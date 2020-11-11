module fowcomm

  USE commpi
  USE plcomm_parm
  USE fpcomm,ONLY:rkind,nrmax,nthmax,npmax,nsamax

  implicit none

  public

  integer:: nthpmax,&                                        ! number of poloidal angle grid points
            nthm1,&                                          ! number of theta_m grid points for 0 <= theta_m <= theta_pnc
            nthm2,&                                          ! number of theta_m grid points for theta_pnc <= theta_m <= theta_co_stg
            nthm3                                            ! number of theta_m grid points for theta_cnt_stg <= theta_m <= pi

  real(rkind),allocatable :: FNSI(:,:,:,:),&     ! distribution function in I=(p,thetam,psim) space
                             Jacobian_I(:,:,:,:) ! dxdydzd(vx)d(vy)d(vz) = Jacobian_I * dpd(thetam)d(psim)
                             
! COM --------------------------------------------------------------------------------------------------          
  real(rkind),allocatable,dimension(:) :: psim,&                ! maximum poloidal magnetic flux in an orbit, value at half integer grid points
                                          psimg                 ! psim at integer grid points
  
  real(rkind),allocatable,dimension(:,:,:,:) :: thetam,&              ! pitch angle along orbit in psim, value at half integer grid points
                                                thetamg,&             ! theta_m at integer grid points
                                                thetam_rg,&           ! theta_m for given pm(np), psimg(nr)
                                                thetam_pg             ! theta_m for given pg(np), psim(nr)

! Diffsion coefficients extend FP's difference equation ------------------------------------------------
  real(rkind),allocatable,dimension(:,:,:,:) :: Dpr,&           ! 
                                                Dtr,&           !
                                                Drp,&           !
                                                Drt             !

! equilibrium variables --------------------------------------------------------------------------------
  real(rkind):: psi0                                             ! poloidal flux at the plasma edge
  real(rkind),allocatable,dimension(:) :: Fpsi,&                 ! poloidai current at psim
                                          Bout,&                 ! abs(magnetic field) at psim when psim is outside the axis
                                          Bin,&                  ! abs(magnetic field) at psim when psim is inside the axis
                                          Fpsig,&                ! Fpsi for grid points
                                          Boutg,&                ! Bout for grid points
                                          Bing                   ! Bin for grid points
  
  real(rkind),allocatable,dimension(:,:) :: Babs                 ! B(psip,thetap)

! use for boundary conditions --------------------------------------------------------------------------
  real(rkind),allocatable,dimension(:,:,:) :: theta_pnc,&         ! theta_m of pinch orbit for given pm(np) and psim(nr)
                                              theta_co_stg,&      ! theta_m of co-stagnation orbit for pm(np) and psim(nr)
                                              theta_cnt_stg,&     ! theta_m of counter-stagnation orbit for pm(np) and psim(nr)
                                              psip_pnc_point,&    ! poloidal flux of pinch orbit given by I = (pm(np), theta_pnc(p,psi_m), psim(nr))
                                              delth1,&            ! grid width of theta_m for 0 <= theta_m <= theta_pnc
                                              delth2,&            ! grid width of theta_m for theta_pnc <= theta_m <= theta_co_stg
                                              delth3              ! grid width of theta_m for theta_cnt_stg <= theta_m <= pi

  real(rkind),allocatable,dimension(:,:,:) :: theta_pnc_pg,&         ! theta_m of pinch orbit for given pg(np) and psim(nr)
                                              theta_co_stg_pg,&      ! theta_m of co-stagnation orbit for given pg(np) and psim(nr)
                                              theta_cnt_stg_pg,&     ! theta_m of counter-stagnation orbit for given pg(np) and psim(nr)
                                              psip_pnc_point_pg,&    ! poloidal flux of pinch orbit given by I = (pg(np), theta_pnc(p,psim), psim(nr))
                                              delth1_pg,&            ! grid width of theta_m for 0 <= theta_m <= theta_pnc
                                              delth2_pg,&            ! grid width of theta_m for theta_pnc <= theta_m <= theta_co_stg
                                              delth3_pg              ! grid width of theta_m for theta_cnt_stg <= theta_m <= pi

  real(rkind),allocatable,dimension(:,:,:) :: theta_pnc_rg,&         ! theta_m of pinch orbit for given pm(np) and psimg(nr)
                                              theta_co_stg_rg,&      ! theta_m of co-stagnation orbit for given pm(np) and psimg(nr)
                                              theta_cnt_stg_rg,&     ! theta_m of counter-stagnation orbit for given pm(np) and psimg(nr)
                                              psip_pnc_point_rg,&    ! poloidal flux of pinch orbit given by I = (pm(np), theta_pnc(p,psi_m), psimg(nr))
                                              delth1_rg,&            ! grid width of theta_m for 0 <= theta_m <= theta_pnc
                                              delth2_rg,&            ! grid width of theta_m for theta_pnc <= theta_m <= theta_co_stg
                                              delth3_rg              ! grid width of theta_m for theta_cnt_stg <= theta_m <= pi

  real(rkind),parameter :: NO_PINCH_ORBIT = 19960610d0 ! if theta_pnc(np,nr,nsa) = NO_PINCH_ORBIT then no pinch orbit exists with pm(np) and psi_m = psim(nr)


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
    allocate(FNSI(nthmax,npmax,nrmax,nsamax),Jacobian_I(nthmax,npmax,nrmax,nsamax))
    allocate(psim(nrmax),psimg(nrmax+1))
    allocate(thetam(nthmax,npmax,nrmax,nsamax),thetamg(nthmax+1,npmax,nrmax,nsamax))
    allocate(thetam_pg(nthmax,npmax+1,nrmax,nsamax), thetam_rg(nthmax,npmax,nrmax+1,nsamax))
    ! 
    allocate(Fpsi(nrmax),Bout(nrmax),Bin(nrmax))
    allocate(Fpsig(nrmax+1),Boutg(nrmax+1),Bing(nrmax+1))
    ! 
    allocate(theta_pnc(npmax,nrmax,nsamax),theta_co_stg(npmax,nrmax,nsamax),theta_cnt_stg(npmax,nrmax,nsamax))
    allocate(psip_pnc_point(npmax,nrmax,nsamax),delth1(npmax,nrmax,nsamax),delth2(npmax,nrmax,nsamax),delth3(npmax,nrmax,nsamax))

    allocate(theta_pnc_pg(npmax+1,nrmax,nsamax),theta_co_stg_pg(npmax+1,nrmax,nsamax),theta_cnt_stg_pg(npmax+1,nrmax,nsamax))
    allocate(psip_pnc_point_pg(npmax+1,nrmax,nsamax),delth1_pg(npmax+1,nrmax,nsamax)&
            ,delth2_pg(npmax+1,nrmax,nsamax),delth3_pg(npmax+1,nrmax,nsamax))

    allocate(theta_pnc_rg(npmax,nrmax+1,nsamax),theta_co_stg_rg(npmax,nrmax+1,nsamax),theta_cnt_stg_rg(npmax,nrmax+1,nsamax))
    allocate(psip_pnc_point_rg(npmax,nrmax+1,nsamax),delth1_rg(npmax,nrmax+1,nsamax)&
            ,delth2_rg(npmax,nrmax+1,nsamax),delth3_rg(npmax,nrmax+1,nsamax))
    allocate(Babs(nrmax+1,nthpmax))
    ! 
    allocate(orbit_p(nthmax,npmax+1,nrmax,nsamax),orbit_th(nthmax+1,npmax,nrmax,nsamax)&
    ,orbit_r(nthmax,npmax,nrmax+1,nsamax),orbit_m(nthmax,npmax,nrmax,nsamax))
    !
    allocate(xi(nthmax))
    allocate(xig(nthmax+1))
  end subroutine fow_allocate

  subroutine fow_deallocate

    !
    deallocate(FNSI,Jacobian_I)
    deallocate(thetam,thetamg)
    deallocate(thetam_pg, thetam_rg)
    ! 
    deallocate(psim,Fpsi,Bout,Bin)
    deallocate(psimg,Fpsig,Boutg,Bing)
    ! 
    deallocate(theta_pnc,theta_co_stg,theta_cnt_stg)
    deallocate(psip_pnc_point,delth1,delth2,delth3)

    deallocate(theta_pnc_pg,theta_co_stg_pg,theta_cnt_stg_pg)
    deallocate(psip_pnc_point_pg,delth1_pg,delth2_pg,delth3_pg)

    deallocate(theta_pnc_rg,theta_co_stg_rg,theta_cnt_stg_rg)
    deallocate(psip_pnc_point_rg,delth1_rg,delth2_rg,delth3_rg)
    ! 
    deallocate(orbit_p,orbit_th,orbit_r,orbit_m)
    !
    deallocate(xi)
    deallocate(xig)

  end subroutine fow_deallocate

  subroutine fow_read_namelist
    namelist /fow/nthpmax

    open(11,file="fpparm",status='old',action='read')
    read(11,nml=fow)
    close(11)
  end subroutine fow_read_namelist

end module fowcomm
