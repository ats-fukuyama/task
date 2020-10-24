module foworbit
  implicit none
  private
  real(8),allocatable :: penergym(:,:),penergyg(:,:)

  public :: fow_orbit_construct, fow_orbit_jacobian, fow_orbit_get_mode

contains

  subroutine fow_orbit_construct(orbit_out)

    use fowcomm,only : orbit
    use fpcomm, only : npmax, nthmax, nrmax, nsamax, rkind

    use plinit
    use obcomm
    use obinit
    use obparm
    use obprep
    use obcalc

    type(orbit),intent(out) :: orbit_out(:,:,:,:)
    integer :: ierr,nobt,nth,np,nr,nsa,mode(3), nobt_in_max
    integer,allocatable :: nobt_in(:,:,:)

    call fow_orbit_get_mode(mode, orbit_out)

    ierr = 0

    ! prepare to execute OB-----------------------------------
    call pl_init
    call EQINIT
    call ob_init
    call ob_parm(1,'../fp.ota/fpparm',ierr)

    call ob_prep(ierr)
    ! --------------------------------------------------------

    ! execute OB ---------------------------------------------
    do nsa = 1, nsamax
      call fow_orbit_prep(nobt_in_max, nobt_in, nsa, mode)
      nobt_max = nobt_in_max
      if ( nobt_in_max>=nobt_m ) then
        write(*,*)"ERROR at fow_orbit_construct : nobt_max must be less than nobt_m"
        write(*,*)"nobt_max = ",nobt_in_max
        write(*,*)"nobt_m   = ",nobt_m
        STOP
      end if
      ns_ob = nsa
      write(*,*)"nsa=",nsa
      call ob_allocate
      call fow_orbit_initial_value(penergy_in, pcangle_in, zeta_in, psipn_in, theta_in, nobt_in, nsa, mode)
      call ob_calc(ierr)
      call fow_construct_orbit(orbit_out, nobt_in, nsa, mode)
      deallocate(nobt_in)
    end do
    ! --------------------------------------------------------

    deallocate(penergym, penergyg)

  end subroutine fow_orbit_construct

  subroutine fow_orbit_jacobian(Jacobian_out, orbit_in)
    use fowcomm
    use fpcomm
    use foworbitclassify
    use fpwrite

    real(rkind),intent(out) :: Jacobian_out(:,:,:,:)
    type(orbit),intent(in) :: orbit_in(:,:,:,:)
    integer :: np, nth, nr, nsa, mode(3)
    real(rkind) :: tau_p, J_c, trans_matrix, Bml, dBmldpsi
    real(rkind),allocatable :: Bml_in(:), Bml_out(:), gml(:), dBm_indpsi(:), dBm_outdpsi(:)&
                            , dBmdpsi(:), dfdpsi(:), lorentz_factor(:,:), xil(:), pl(:,:)

    call fow_orbit_get_mode(mode, orbit_in)
    
    allocate( lorentz_factor(npmax+mode(2),nsamax) )
    allocate( xil(nthmax+mode(1)), pl(npmax+mode(2),nsamax) )
    allocate( Bml_in(nrmax+mode(3)), Bml_out(nrmax+mode(3)), gml(nrmax+mode(3))&
    , dBm_indpsi(nrmax+mode(3)), dBm_outdpsi(nrmax+mode(3)), dfdpsi(nrmax+mode(3)))

    select case(mode(1))
    case(0)
      xil(:) = xi(:)
    case(1)
      xil(:) = xig(:)
    end select

    select case(mode(2))
    case(0)
      do nsa = 1, nsamax
        do np = 1, npmax
          pl(np,nsa) = pm(np,nsa)*ptfp0(nsa)
        end do  
      end do
      
      do nsa = 1, nsamax
        do np = 1, npmax
          lorentz_factor(np,nsa) = sqrt(1.d0+theta0(nsa)*pm(np,nsa)**2)
        end do
      end do

    case(1)
      do nsa = 1, nsamax
        do np = 1, npmax+1
          pl(np,nsa) = pg(np,nsa)*ptfp0(nsa)
        end do  
      end do

      do nsa = 1, nsamax
        do np = 1, npmax+1
          lorentz_factor(np,nsa) = sqrt(1.d0+theta0(nsa)*pg(np,nsa)**2)
        end do
      end do
    end select

    select case(mode(3))
    case(0)
      Bml_in(:) = Bin(:)
      Bml_out(:) = Bout(:)
      gml(:) = Fpsi(:)
      do nr = 1, nrmax+mode(3)
        if ( nr==1 ) then
          dBm_outdpsi(nr) = (-3*Bout(1)+4*Bout(2)-Bout(3))/(-3*psim(1)+4*psim(2)-psim(3))
          dBm_indpsi(nr) = (-3*Bin(1)+4*Bin(2)-Bin(3))/(-3*psim(1)+4*psim(2)-psim(3))
          dFdpsi(nr) = (-3*Fpsi(1)+4*Fpsi(2)-Fpsi(3))/(-3*psim(1)+4*psim(2)-psim(3))
        else if ( nr==nrmax ) then
          dBm_outdpsi(nr) = (3*Bout(nrmax)-4*Bout(nrmax-1)+Bout(nrmax-2))&
                            /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
          dBm_indpsi(nr) = (3*Bin(nrmax)-4*Bin(nrmax-1)+Bin(nrmax-2))&
                            /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
          dFdpsi(nr) = (3*Fpsi(nrmax)-4*Fpsi(nrmax-1)+Fpsi(nrmax-2))&
                        /(3*psim(nrmax)-4*psim(nrmax-1)+psim(nrmax-2))
        else
          dBm_outdpsi(nr) = (Bout(nr+1)-Bout(nr-1))/(psim(nr+1)-psim(nr-1))
          dBm_indpsi(nr) = (Bin(nr+1)-Bin(nr-1))/(psim(nr+1)-psim(nr-1))
          dFdpsi(nr) = (Fpsi(nr+1)-Fpsi(nr-1))/(psim(nr+1)-psim(nr-1))
        end if
      end do
    case(1)
      Bml_in(:) = Bing(:)
      Bml_out(:) = Boutg(:)
      gml(:) = Fpsig(:)
      do nr = 1, nrmax+mode(3)
        if ( nr==1 ) then
          dBm_outdpsi(nr) = (-3*Boutg(1)+4*Boutg(2)-Boutg(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
          dBm_indpsi(nr) = (-3*Bing(1)+4*Bing(2)-Bing(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
          dFdpsi(nr) = (-3*Fpsig(1)+4*Fpsig(2)-Fpsig(3))/(-3*psimg(1)+4*psimg(2)-psimg(3))
        else if ( nr==nrmax+1 ) then
          dBm_outdpsi(nr) = (3*Boutg(nrmax+1)-4*Boutg(nrmax)+Boutg(nrmax-1))&
                            /(3*psimg(nrmax+1)-4*psimg(nrmax)+psimg(nrmax-1))
          dBm_indpsi(nr) = (3*Bing(nrmax+1)-4*Bing(nrmax)+Bing(nrmax-1))&
                            /(3*psimg(nrmax+1)-4*psimg(nrmax)+psimg(nrmax-1))
          dFdpsi(nr) = (3*Fpsig(nrmax+1)-4*Fpsig(nrmax)+Fpsig(nrmax-1))&
                        /(3*psimg(nrmax+1)-4*psimg(nrmax)+psimg(nrmax-1))
        else
          dBm_outdpsi(nr) = (Boutg(nr+1)-Boutg(nr-1))/(psimg(nr+1)-psimg(nr-1))
          dBm_indpsi(nr) = (Bing(nr+1)-Bing(nr-1))/(psimg(nr+1)-psimg(nr-1))
          dFdpsi(nr) = (Fpsig(nr+1)-Fpsig(nr-1))/(psimg(nr+1)-psimg(nr-1))
        end if
      end do
    end select

    call fpcsv1D(dBm_outdpsi,"./csv/dBm_outdpsi.csv")
    call fpcsv1D(dBm_indpsi,"./csv/dBm_indpsi.csv")
    call fpcsv1D(dFdpsi,"./csv/dFdpsi.csv")

    do nsa = 1, nsamax
      do nr = 1, nrmax+mode(3)
        do np = 1, npmax+mode(2)
          do nth = 1, nthmax+mode(1)
            if ( forbitten(nth,np,nr,nsa,mode) ) then
              Jacobian_out(nth,np,nr,nsa) = 0.d0
            else
              if ( xil(nth)>=0.d0 ) then
                Bml = Bml_out(nr)
                dBmldpsi = dBm_outdpsi(nr)
              else
                Bml = Bml_in(nr)
                dBmldpsi = dBm_indpsi(nr)
              end if
              tau_p = orbit_in(nth,np,nr,nsa)%time(orbit_in(nth,np,nr,nsa)%nstp_max)
              write(*,*)
              J_c = 4.d0*pi**2*tau_p/(amfp(nsa)**2*abs(aefp(nsa)))
              trans_matrix = pl(np,nsa)**3/(amfp(nsa)**2*lorentz_factor(np,nsa)*Bml)&
                            *((-Bml*dFdpsi(nr)*xil(nth)**2+0.5d0*(1.d0+xil(nth)**2)*dBmldpsi*gml(nr))*pl(np,nsa)/Bml**2&
                            +aefp(nsa)*xil(nth))
              Jacobian_out(nth,np,nr,nsa) = J_c*trans_matrix  
            end if
          end do
        end do
      end do
    end do

    deallocate( lorentz_factor )
    deallocate( xil, pl )
    deallocate( Bml_in, Bml_out, gml, dBm_indpsi, dBm_outdpsi, dfdpsi)

  end subroutine fow_orbit_jacobian

  subroutine fow_orbit_prep(nobt_in_max, nobt_in, nsa_in, mode)

    use fowcomm
    use fpcomm

    integer,intent(in) :: mode(3), nsa_in
    integer,intent(out) :: nobt_in_max
    integer,allocatable,intent(out) :: nobt_in(:,:,:)
    integer :: nth,np,nr,nsa,i
    real(rkind) :: PVM,PVG

    if ( (.not.allocated(penergym)) .or. (.not.allocated(penergyg)) ) then
      allocate(penergym(npmax,nsamax))
      allocate(penergyg(npmax+1,nsamax))
      do nsa=1,nsamax
        do np=1,npmax+1
          PVG = SQRT(1.D0+THETA0(nsa)*PG(np,nsa)**2)
          penergyg(np,nsa) = (PVG-1.D0)*PTFP0(nsa)**2/(THETA0(nsa)*AMFP(nsa))/(1.d3*aee)
          if( np/=npmax+1 ) then
            PVM = SQRT(1.D0+THETA0(nsa)*PM(np,nsa)**2)
            penergym(np,nsa) = (PVM-1.D0)*PTFP0(nsa)**2/(THETA0(nsa)*AMFP(nsa))/(1.d3*aee)
          end if
        end do
      end do  
    end if

    if ( allocated(nobt_in) ) deallocate(nobt_in)
    allocate(nobt_in(nthmax+mode(1),npmax+mode(2),nrmax+mode(3)))

    i = 0
    do nr = 1, nrmax+mode(3)
      do np = 1, npmax+mode(2)
        do nth = 1, nthmax+mode(1)
          i = i+1
          nobt_in(nth,np,nr) = i
        end do
      end do
    end do

    nobt_in_max = i

  end subroutine fow_orbit_prep

  subroutine fow_orbit_initial_value(penergy_in, pcangle_in, zeta_in, psipn_in, theta_in, nobt_in, nsa_in, mode)

    use fowcomm
    use fpcomm

    integer,intent(in) :: mode(3), nobt_in(:,:,:), nsa_in
    real(rkind),intent(inout) :: penergy_in(:), pcangle_in(:), zeta_in(:), psipn_in(:), theta_in(:)
    integer :: i,nth,np,nr,nsa

    do nr = 1, nrmax+mode(3)
      do np = 1, npmax+mode(2)
        do nth = 1, nthmax+mode(1)

          i = nobt_in(nth,np,nr)

          select case(mode(1))
          case(0)
            pcangle_in(i) = cos(thetam(nth,np,nr,nsa_in))
          case(1)
            pcangle_in(i) = cos(thetamg(nth,np,nr,nsa_in))
          end select

          select case(mode(2))
          case(0)
            penergy_in(i) = penergym(np,nsa_in)
          case(1)
            penergy_in(i) = penergyg(np,nsa_in)
          end select

          select case(mode(3))
          case(0)
            psipn_in(i) = psim(nr)/psi0
          case(1)
            psipn_in(i) = psimg(nr)/psi0
          end select

          if ( pcangle_in(i) >= 0.d0 ) then
            theta_in(i) = 0.d0
          else
            theta_in(i) = pi
          end if

          zeta_in(i) = 0.d0

        end do
      end do
    end do
  
  end subroutine fow_orbit_initial_value

  subroutine fow_construct_orbit(orbit_in, nobt_in, nsa_in, mode)
    use obcomm
    use fowcomm,only : orbit
    use fpcomm, only : npmax, nthmax, nrmax, nsamax, rkind

    integer,intent(in) :: mode(3), nobt_in(:,:,:), nsa_in
    type(orbit),intent(inout) :: orbit_in(:,:,:,:)
    real(rkind), allocatable :: construct_input(:,:)
    integer :: nsa,nr,nth,np
    integer :: i, j

    allocate(construct_input(4,nstp_max+1))
    ! time_ob(0:nstp_max,nobt_max)
    do nr = 1, nrmax+mode(3)
      do np = 1, npmax+mode(2)
        do nth = 1, nthmax+mode(1)

          i = nobt_in(nth,np,nr)

          do j = 1, nstp_max+1
            construct_input(1,j)  = time_ob(j-1,i)
            construct_input(2,j)  = psip_ob(j-1,i)
            construct_input(3,j)  = babs_ob(j-1,i)
            construct_input(4,j)  = acos(vpara_ob(j-1,i)/sqrt(vperp_ob(j-1,i)**2+vpara_ob(j-1,i)**2))
          end do

          orbit_in(nth,np,nr,nsa_in) = construct_orbit(nstp_max_nobt(i), construct_input)  

        end do
      end do
    end do

  end subroutine

  subroutine fow_orbit_get_mode(mode, orbit_in)
    use fowcomm 
    use fpcomm

    integer, intent(out) :: mode(3)
    type(orbit),intent(in) :: orbit_in(:,:,:,:)
    integer :: npm, nthm, nrm

    nthm = size(orbit_in,1)
    npm = size(orbit_in,2)
    nrm = size(orbit_in,3)

    if ( nthm == nthmax ) then
      mode(1) = 0
    else if ( nthm == nthmax+1 )then
      mode(1) = 1
    else
      write(*,*)"ERROR : size(orbit_in,1) is not nthm = nthmax or nthmax+1 in fow_orbit_get_mode"
    end if

    if ( npm == npmax ) then
      mode(2) = 0
    else if ( npm == npmax+1 )then
      mode(2) = 1
    else
      write(*,*)"ERROR : size(orbit_in,2) is not npm = npmax or npmax+1 in fow_orbit_get_mode"
    end if

    if ( nrm == nrmax ) then
      mode(3) = 0
    else if ( nrm == nrmax+1 )then
      mode(3) = 1
    else
      write(*,*)"ERROR : size(orbit_in,3) is not nrm = nrmax or nrmax+1 in fow_orbit_get_mode"
    end if

  end subroutine fow_orbit_get_mode

  function func_orbit_F(orbit_in, nstp_in, nr_in) result(F_ret)
    use fpcomm
    use fowcomm
    implicit none
    real(rkind) :: F_ret
    type(orbit),intent(in) :: orbit_in
    integer,intent(in) :: nstp_in, nr_in

    real(rkind) :: C(3), rr_l
    integer :: nr_max, nr_min, ierr, nr_l

    nr_max = max(nr_in, 4)
    nr_min = max(1, nr_in-3)

    call leastSquareMethodForQuadric(psimg,nrmax+1,nr_min,nr_max,C)
    C(3) = C(3)-orbit_in%psip(nstp_in)
    
    rr_l = C(2)**2-4.d0*C(1)*C(3)
    if ( rr_l < 0.d0 ) then
      F_ret = Fpsi(nr_in)
      return
    else
      rr_l = (-1.d0*C(2)+sqrt(rr_l))/(2*C(1))
    end if

    if ( rr_l <= 0.d0 ) then
      F_ret = Fpsig(1)
    else if ( rr_l >= (nrmax+1)*1.d0 ) then
      F_ret = Fpsig(nrmax+1)
    else 
      nr_l = int(rr_l)
      F_ret = Fpsig(nr_l+1)+(Fpsig(nr_l+1)-Fpsig(nr_l))*(rr_l-dble(nr_l+1))
    end if

    return

  end function func_orbit_F

  function construct_orbit(n,input) result(ret)
    use fowcomm
    type(orbit) :: ret
    integer, intent(in) :: n
    real(rkind), intent(in) :: input(:,:)
    integer :: i

    allocate(ret%time(n+1), ret%psip(n+1), ret%Babs(n+1), ret%theta(n+1))

    ret%nstp_max=n+1
    do i=1,n+1
      ret%time(i) = input(1,i) 
      ret%psip(i) = input(2,i)
      ret%Babs(i) = input(3,i)
      ret%theta(i) = input(4,i)
    end do

  end function construct_orbit

  subroutine destruct_orbit(orbit_in)
    use fowcomm
    type(orbit) :: orbit_in
    deallocate(orbit_in%time, orbit_in%psip, orbit_in%babs, orbit_in%theta)
  end subroutine destruct_orbit

  subroutine leastSquareMethodForQuadric(y,nxmax,nxstart,nxend,reval)
    ! y=reval(1)*nx^2+reval(2)*nx+reval(3)
    use fowprep,only:gauss_jordan
    implicit none
    integer,intent(in) :: nxmax,nxstart,nxend
    double precision,intent(in) :: y(nxmax)
    double precision,intent(out) :: reval(3)
    double precision :: A(3,3)
    integer :: i,j,n
  
    reval(:)=0
    A(:,:)=0
  
    do i=1,3
      do n=nxstart,nxend
        reval(i)=reval(i)+n**(3-i)*y(n)
      end do
      do j=1,3
        do n=nxstart,nxend
          A(i,j)=A(i,j)+n**(6-i-j)
        end do
      end do
    end do
  
    call gauss_jordan(A,reval,3)
  
  end subroutine leastSquareMethodForQuadric
  
end module foworbit
