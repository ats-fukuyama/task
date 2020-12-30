module fowcomm

  USE commpi
  USE plcomm_parm
  USE fpcomm,ONLY:rkind,nrmax,nthmax,npmax,nsamax

  implicit none

  public

  integer :: model_obload   ! 0          : exec TASK/OB anyway and do not save orbit_x to binary files
                            ! 1[default] : If binary file of 'type(orbit) orbit_x' exists then load member of orbit_x, else exec TASK/OB

  integer :: model_mkcsv    ! 0[default] : do not output to csv files
                            ! 1          : output to csv files
  integer :: max_stp        ! maximum step number for bounce average

  integer:: nthpmax,&                                        ! number of poloidal angle grid points
            nthm1,&                                          ! number of theta_m grid points for 0 <= theta_m <= theta_pnc
            nthm2,&                                          ! number of theta_m grid points for theta_pnc <= theta_m <= theta_co_stg
            nthm3                                            ! number of theta_m grid points for theta_cnt_stg <= theta_m <= pi

  real(rkind),allocatable :: JI(:,:,:,:),& ! dxdydzd(vx)d(vy)d(vz) = JI * dpd(thetam)d(psim)
                             JIR(:,:,:,:)  ! for integrate over velocity space
  real(rkind),allocatable  :: fnorm(:)
! COM --------------------------------------------------------------------------------------------------          
  real(rkind),allocatable,dimension(:) :: psim,&                ! maximum poloidal magnetic flux in an orbit, value at half integer grid points
                                          psimg,&               ! psim at integer grid points
                                          delps

  real(rkind),allocatable :: psim_local(:,:,:,:,:), &
                            thetam_local(:,:,:,:,:), &
                            time_loss(:,:,:,:,:)

  real(rkind),allocatable,dimension(:,:,:,:) :: thetam,&              ! pitch angle along orbit in psim, value at half integer grid points
                                                thetamg,&             ! theta_m at integer grid points
                                                thetam_rg,&           ! theta_m for given pm(np), psimg(nr)
                                                thetam_pg             ! theta_m for given pg(np), psim(nr)

! Diffsion coefficients extend FP's difference equation ------------------------------------------------
  real(rkind),allocatable,dimension(:,:,:,:) :: Dppfow, Dptfow, Dprfow,&
                                                Dtpfow, Dttfow, Dtrfow,&
                                                Drpfow, Drtfow, Drrfow,&
                                                Fppfow,  Fthfow,  Frrfow

! equilibrium variables --------------------------------------------------------------------------------
  real(rkind):: psi0                                             ! poloidal flux at the plasma edge
  real(rkind),allocatable,dimension(:) :: Fpsi,&                 ! poloidai current at psim
                                          Bout,&                 ! abs(magnetic field) at psim when psim is outside the axis
                                          Bin,&                  ! abs(magnetic field) at psim when psim is inside the axis
                                          Fpsig,&                ! Fpsi for grid points
                                          Boutg,&                ! Bout for grid points
                                          Bing                   ! Bin for grid points
  
  real(rkind),allocatable,dimension(:,:) :: Babs                 ! B(psip,thetap)
  real(rkind),allocatable,dimension(:) :: theta_p                 ! poloidal angle

! use for boundary conditions --------------------------------------------------------------------------
  real(rkind),allocatable,dimension(:,:,:) :: theta_pnc,&         ! theta_m of pinch orbit for given pm(np) and psim(nr)
                                              theta_co_stg,&      ! theta_m of co-stagnation orbit for pm(np) and psim(nr)
                                              theta_cnt_stg,&     ! theta_m of counter-stagnation orbit for pm(np) and psim(nr)
                                              psip_pnc_point      ! poloidal flux of pinch orbit given by I = (pm(np), theta_pnc(p,psi_m), psim(nr))
                                              

  real(rkind),allocatable,dimension(:,:,:) :: theta_pnc_pg,&         ! theta_m of pinch orbit for given pg(np) and psim(nr)
                                              theta_co_stg_pg,&      ! theta_m of co-stagnation orbit for given pg(np) and psim(nr)
                                              theta_cnt_stg_pg,&     ! theta_m of counter-stagnation orbit for given pg(np) and psim(nr)
                                              psip_pnc_point_pg      ! poloidal flux of pinch orbit given by I = (pg(np), theta_pnc(p,psim), psim(nr))
                                              

  real(rkind),allocatable,dimension(:,:,:) :: theta_pnc_rg,&         ! theta_m of pinch orbit for given pm(np) and psimg(nr)
                                              theta_co_stg_rg,&      ! theta_m of co-stagnation orbit for given pm(np) and psimg(nr)
                                              theta_cnt_stg_rg,&     ! theta_m of counter-stagnation orbit for given pm(np) and psimg(nr)
                                              psip_pnc_point_rg      ! poloidal flux of pinch orbit given by I = (pm(np), theta_pnc(p,psi_m), psimg(nr))

  real(rkind),allocatable,dimension(:,:,:,:) :: delthm_rg, delthm_pg, delthm

  integer,allocatable,dimension(:) :: nth_stg,&        ! thetamg(nth_stg,(nsa),np,nr,nsa) = theta_stg
                                      nth_pnc           ! thetamg(nth_pnc,(nsa),np,nr,nsa) = theta_pnc

  real(rkind),allocatable,dimension(:,:,:) :: IBCflux_ratio         ! ration between the flux from thetam(nth_pnc-1) to thetam(nth_pnc) and the flux from thetam(nth_pnc-1) to thetam(nth_co_stg)
  integer,allocatable,dimension(:,:,:) ::  nr_pnc_point          ! nr of pinch point of theta_pnc(:,:,:)

  real(rkind),parameter :: NO_PINCH_ORBIT = 19960610d0 ! if theta_pnc(np,nr,nsa) = NO_PINCH_ORBIT then no pinch orbit exists with pm(np) and psi_m = psim(nr)

  type Xstg_as_pnc_point
    integer :: number                   ! number of pinch orbits whose pinch point is at psim of this X-type stagnation orbit(nth,np,nr,nsa)
                                        ! the orbit at (nth,np,nr,nsa) is not X-type stagnation orbit or no pinch orbit has pinch point at psim of X-type stagnation orbit or -> number = 0
    integer,allocatable :: nr(:) ! np and nr of pinch orbits
  end type Xstg_as_pnc_point

  type(Xstg_as_pnc_point),allocatable,dimension(:,:,:) :: Xstg_as_pncp

! use for bounce average -------------------------------------------------------------------------------
  type :: orbit                                                  ! quantities along orbits culcurated by TASK/OB
    integer :: nstp_max
    real(rkind),allocatable,dimension(:) :: time, &              ! time
                                            psip, &              ! poloidal magnetic flux
                                            Babs, &              ! absolute value of magnetic field
                                            theta,&              ! pitch angle
                                            thetap               ! poloidal angle
  end type orbit

  type(orbit),allocatable,dimension(:,:,:,:) :: orbit_p,&    ! (np,nth,nr)=(-0.5,  1 ,  1 ) to (npmax+0.5,nthmax,nrmax)
                                                orbit_th,&   ! (np,nth,nr)=(  1 ,-0.5,  1 ) to (npmax,nthmax+0.5,nrmax)
                                                orbit_r,&    ! (np,nth,nr)=(  1 ,  1 ,-0.5) to (npmax,nthmax,nrmax+0.5)
                                                orbit_m      ! (np,nth,nr)=(1,1,1) to (npmax,nthmax,nrmax)

! ------------------------------------------------------------------------------------------------------

  real(rkind),allocatable,dimension(:) :: xi,&                 ! cos(pitch angle) at psim
                                          xig                  ! xi for grid points

contains

  subroutine fow_allocate
    use fpcomm, only:npmax,nthmax,nrmax,nsamax
    allocate(JI(nthmax,npmax,nrmax,nsamax),JIR(nthmax,npmax,nrmax,nsamax))
    allocate(fnorm(nsamax))
    allocate(psim(nrmax),psimg(nrmax+1))
    allocate(delps(nrmax))
    allocate(psim_local(nthmax,npmax,nrmax,nthpmax,nsamax), thetam_local(nthmax,npmax,nrmax,nthpmax,nsamax))
    allocate(time_loss(nthmax,npmax,nrmax,nthpmax,nsamax))

    allocate(thetam(nthmax,npmax,nrmax,nsamax),thetamg(nthmax+1,npmax,nrmax,nsamax))
    allocate(thetam_pg(nthmax,npmax+1,nrmax,nsamax), thetam_rg(nthmax,npmax,nrmax+1,nsamax))
    ! 
    allocate(Fpsi(nrmax),Bout(nrmax),Bin(nrmax))
    allocate(Fpsig(nrmax+1),Boutg(nrmax+1),Bing(nrmax+1))
    ! 
    allocate(theta_pnc(npmax,nrmax,nsamax),theta_co_stg(npmax,nrmax,nsamax),theta_cnt_stg(npmax,nrmax,nsamax))
    allocate(psip_pnc_point(npmax,nrmax,nsamax))

    allocate(theta_pnc_pg(npmax+1,nrmax,nsamax),theta_co_stg_pg(npmax+1,nrmax,nsamax),theta_cnt_stg_pg(npmax+1,nrmax,nsamax))
    allocate(psip_pnc_point_pg(npmax+1,nrmax,nsamax))

    allocate(theta_pnc_rg(npmax,nrmax+1,nsamax),theta_co_stg_rg(npmax,nrmax+1,nsamax),theta_cnt_stg_rg(npmax,nrmax+1,nsamax))
    allocate(psip_pnc_point_rg(npmax,nrmax+1,nsamax))

    allocate(delthm(nthmax,npmax,nrmax,nsamax),delthm_pg(nthmax,npmax+1,nrmax,nsamax),delthm_rg(nthmax,npmax,nrmax+1,nsamax))
    allocate(nth_stg(nsamax),nth_pnc(nsamax))
    allocate(IBCflux_ratio(npmax,nrmax,nsamax))
    allocate(nr_pnc_point(npmax,nrmax,nsamax))

    allocate(Xstg_as_pncp(npmax,nrmax,nsamax))

    allocate(Babs(nrmax+1,nthpmax))
    allocate(theta_p(nthpmax))
    ! 
    allocate(orbit_p(nthmax,npmax+1,nrmax,nsamax),orbit_th(nthmax+1,npmax,nrmax,nsamax)&
    ,orbit_r(nthmax,npmax,nrmax+1,nsamax),orbit_m(nthmax,npmax,nrmax,nsamax))
    !
    allocate(Dppfow(nthmax,npmax+1,nrmax,nsamax),Dptfow(nthmax,npmax+1,nrmax,nsamax),Dprfow(nthmax,npmax+1,nrmax,nsamax))
    allocate(Dtpfow(nthmax+1,npmax,nrmax,nsamax),Dttfow(nthmax+1,npmax,nrmax,nsamax),Dtrfow(nthmax+1,npmax,nrmax,nsamax))
    allocate(Drpfow(nthmax,npmax,nrmax+1,nsamax),Drtfow(nthmax,npmax,nrmax+1,nsamax),Drrfow(nthmax,npmax,nrmax+1,nsamax))
    allocate(Fppfow(nthmax,npmax+1,nrmax,nsamax),Fthfow(nthmax+1,npmax,nrmax,nsamax),Frrfow(nthmax,npmax,nrmax+1,nsamax))
    !
    allocate(xi(nthmax))
    allocate(xig(nthmax+1))
  end subroutine fow_allocate

  subroutine fow_deallocate

    !
    deallocate(JI)
    deallocate(thetam,thetamg)
    deallocate(thetam_pg, thetam_rg)
    ! 
    deallocate(psim,Fpsi,Bout,Bin)
    deallocate(psimg,Fpsig,Boutg,Bing)
    ! 
    deallocate(theta_pnc,theta_co_stg,theta_cnt_stg)
    deallocate(psip_pnc_point, delthm)

    deallocate(theta_pnc_pg,theta_co_stg_pg,theta_cnt_stg_pg)
    deallocate(psip_pnc_point_pg, delthm_pg)

    deallocate(theta_pnc_rg,theta_co_stg_rg,theta_cnt_stg_rg)
    deallocate(psip_pnc_point_rg, delthm_rg)
    ! 
    deallocate(orbit_p,orbit_th,orbit_r,orbit_m)
    !
    deallocate(xi)
    deallocate(xig)

  end subroutine fow_deallocate

  subroutine fow_read_namelist
    namelist /fow/nthpmax, max_stp, model_obload, model_mkcsv

    open(11,file="fpparm",status='old',action='read')
    read(11,nml=fow)
    close(11)
  end subroutine fow_read_namelist

  subroutine mesh_to_grid1D(f,g)
    use fpcomm,only:rkind
    implicit none
    real(rkind) :: f(:),g(:)
    real(rkind),allocatable :: x(:),fx(:),U(:,:)
    integer :: i,j,imax,jmax,IERR = 0,k

    imax = size(f)
    jmax = size(g)

    if(imax/=jmax-1)then
      write(*,*)"imax/ = jmax-1 at subroutine mesh_to_grid1D"
      STOP
    end if

    allocate(x(imax),fx(imax),U(4,imax))

    do i = 1,imax
      x(i) = i*1.d0
    end do

    call SPL1D(x,f,fx,U,imax,0,IERR)

    do j = 2,jmax-1
      g(j) = 0.d0
      do k = 1,4
        g(j) = g(j)+U(k,j)*0.5d0**(k-1)
      end do
    end do

    g(1) = 0.d0
    g(jmax) = 0.d0
    do k = 1,4
      g(1) = g(1)+U(k,2)*(-0.5d0)**(k-1)
      g(jmax) = g(jmax)+U(k,imax)*(1.5d0)**(k-1)
    end do

    deallocate(x,fx,U)

  end subroutine mesh_to_grid1D

  subroutine solve_quadratic_equation(z,C)
    ! solve C(3)*z**2+C(2)*z+C(1) = 0
    ! z(1) : -sqrt(D)
    ! z(2) : +sqrt(D)
    use fpcomm,only:rkind

    implicit none
    real(rkind),intent(in) :: C(3)
    complex(rkind),intent(out) :: z(2)
    real(rkind) :: D
    complex(rkind) :: ei=(0.d0,1.d0)

    z(1) = (0.d0, 0.d0)
    z(2) = (0.d0, 0.d0)

    D = C(2)**2-4.d0*C(1)*C(3)

    if ( D >= 0.d0 ) then
      z(1) = (-C(2)-sqrt(D))/(2.d0*C(3))+0.d0*ei
      z(2) = (-C(2)+sqrt(D))/(2.d0*C(3))+0.d0*ei
    else
      z(1) = (-C(2)-ei*sqrt(-D))/(2.d0*C(3))
      z(2) = (-C(2)+ei*sqrt(-D))/(2.d0*C(3))
    end if

  end subroutine solve_quadratic_equation

  subroutine first_order_derivative(dfdx, f, x)
    ! calcurate dfdx, the first order derivative of f
    ! error order is O(dx**2)
    use fpcomm,only:rkind

    implicit none
    real(rkind),intent(out) :: dfdx(:)
    real(rkind),intent(in)  :: f(:), x(:)
    real(rkind) :: s, t
    integer :: imax, i

    imax = size(f)
    if ( imax /= size(x) .or. imax /= size(dfdx) ) then
      write(*,*)"ERROR at subroutine first_order_derivative in fowprep.f90"
    end if

    do i = 1, imax
      if ( i /= imax .and. i /= 1) then
        t = x(i+1)-x(i)
        s = x(i)-x(i-1)
        dfdx(i) = (s**2*f(i+1)+(t**2-s**2)*f(i)-t**2*f(i-1))/(s*t*(s+t))
      else if ( i == 1 ) then
        s = x(2)-x(1)
        t = x(3)-x(2)
        dfdx(1) = ((s**2-(s+t)**2)*f(1)+(s+t)**2*f(2)-s**2*f(3))/(s*t*(s+t))
      else if ( i == imax ) then
        t = x(imax)-x(imax-1)
        s = x(imax-1)-x(imax-2)
        dfdx(imax) = (((s+t)**2-t**2)*f(imax)-(s+t)**2*f(imax-1)+t**2*f(imax-2))/(s*t*(s+t))
      end if
    end do

  end subroutine first_order_derivative

  subroutine second_order_derivative(d2fdx2, f, x)
    ! calcurate d2fdx2, the second order derivative of f
    ! error order is O(dx**2)

    use fpcomm,only:rkind

    implicit none
    real(rkind),intent(out) :: d2fdx2(:)
    real(rkind),intent(in)  :: f(:), x(:)
    real(rkind) :: s, t, r, v, w, A(3,3), B(3)
    integer :: imax, i

    imax = size(f)

    do i = 1, imax
      if ( i /= imax .and. i /= 1) then
        t = x(i+1)-x(i)
        s = x(i)-x(i-1)
        d2fdx2(i) = (s*f(i+1)-(s+t)*f(i)+t*f(i-1))/(0.5d0*s*t*(s+t))
      else if ( i == 1 ) then
        r = x(2)-x(1)
        s = x(3)-x(2)
        t = x(4)-x(3)
        v = r+s
        w = r+s+t

        A(1,1) = r
        A(1,2) = r**2/2
        A(1,3) = r**3/6

        A(2,1) = v
        A(2,2) = v**2/2
        A(2,3) = v**3/6

        A(3,1) = w
        A(3,2) = w**2/2
        A(3,3) = w**3/6

        B(1) = f(2)-f(1)
        B(2) = f(3)-f(1)
        B(3) = f(4)-f(1)

        call gauss_jordan(A, B, 3)

        d2fdx2(1) = B(2)

      else if ( i == imax ) then
        r = x(imax-2)-x(imax-3)
        s = x(imax-1)-x(imax-2)
        t = x(imax)-x(imax-1)
        v = s+t
        w = r+s+t

        A(1,1) = -1.d0*t
        A(1,2) = t**2/2
        A(1,3) = -1.d0*t**3/6

        A(2,1) = -1.d0*v
        A(2,2) = v**2/2
        A(2,3) = -1.d0*v**3/6

        A(3,1) = -1.d0*w
        A(3,2) = w**2/2
        A(3,3) = -1.d0*w**3/6

        B(1) = f(imax-1)-f(imax)
        B(2) = f(imax-2)-f(imax)
        B(3) = f(imax-3)-f(imax)

        call gauss_jordan(A, B, 3)

        d2fdx2(imax) = B(2)
        
      end if
    end do

  end subroutine second_order_derivative

  subroutine gauss_jordan(A, B, n)
    use fpcomm,only:rkind

    implicit none
    real(rkind) :: A(n,n), B(n)
    integer n,i,j,k
  
    do k = 1, n
      do j = k + 1, n
        a(k,j) = a(k,j) / a(k,k)
      end do
      b(k) = b(k) / a(k,k)
  
      do i = 1, n
        if ( i .ne. k ) then
          do j = k + 1, n
            a(i,j) = a(i,j) - a(i,k) * a(k,j)
          end do
          b(i) = b(i) - a(i,k) * b(k)
        end if
      end do
  
    end do
  
  end subroutine gauss_jordan
  
  recursive function func_kaijou(n) result(m)
    implicit none
    integer,intent(in) :: n
    integer :: m

    if(n == 1) then
      m = 1
    else
      m = n*func_kaijou(n-1)
    end if

  end function func_kaijou

  subroutine fow_cal_spl(f_out, x_in, f, x)
    ! Calculate spline coefficient y = f(x),
    ! then return f_out = f(x_in)
    use fpcomm,only:rkind

    implicit none
    real(rkind),intent(out) :: f_out
    real(rkind),intent(in) :: x_in, f(:), x(:)
    integer :: i,imax,ierr=0
    real(rkind),allocatable :: U(:,:), fx(:)

    imax = size(f,1)

    allocate(U(4,imax), fx(imax))

    call first_order_derivative(fx, f, x)

    call SPL1D(x,f,fx,U,imax,3,IERR)

    call SPL1DF(x_in,f_out,x,U,imax,IERR)

    return

  end subroutine fow_cal_spl

  subroutine fow_cal_spl2D(f_out, x_in, y_in, f, x, y)
    ! 2D-version of fow_cal_spl
    use fpcomm,only:rkind

    implicit none
    real(rkind),intent(out) :: f_out
    real(rkind),intent(in) :: x_in, y_in, f(:,:), x(:), y(:)
    integer :: i, j, imax, jmax, ierr = 0
    real(rkind),allocatable :: U(:,:,:,:), fx(:,:), fy(:,:) , fxy(:,:)

    imax = size(x)
    jmax = size(y)

    allocate(U(4,4,imax,jmax),fx(imax,jmax),fy(imax,jmax),fxy(imax,jmax))

    if ( size(f,1) /= imax &
        .or. size(f,2) /= jmax ) then
      write(*,*)"error at fow_cal_spl2D"
      STOP
    end if

    call SPL2D(x,y,f,fx,fy,fxy,U,imax,imax,jmax,0,0,IERR)

    call SPL2DF(x_in,y_in,f_out,x,y,U,imax,imax,jmax,ierr)

  end subroutine

end module fowcomm
