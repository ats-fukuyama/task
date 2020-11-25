module foworbit
  implicit none
  private
  real(8),allocatable :: penergym(:,:),penergyg(:,:)

  public :: fow_orbit_construct, fow_orbit_save
  public :: func_orbit_F, func_orbit_mean_ra

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
          penergyg(np,nsa) = amfp(nsa)*vc**2*(PVG-1.D0)/(1.d3*aee)
          if( np /= npmax+1 ) then
            PVM = SQRT(1.D0+THETA0(nsa)*PM(np,nsa)**2)
            penergym(np,nsa) = amfp(nsa)*vc**2*(PVM-1.D0)/(1.d3*aee)
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
          ! exclude forbitten region
          if ( mode(1) == 0 .and. nth == nth_forbitten(nsa_in) ) cycle
          if ( mode(2) == 1 .and. np == 1 ) cycle
          if ( mode(3) == 1 .and. nr == 1 ) cycle

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
          ! exclude forbitten region
          if ( mode(1) == 0 .and. nth == nth_forbitten(nsa_in) ) cycle
          if ( mode(2) == 1 .and. np == 1 ) cycle
          if ( mode(3) == 1 .and. nr == 1 ) cycle

          i = nobt_in(nth,np,nr)

          select case(mode(1))
          case(0)
            if ( mode(2) == 0 .and. mode(3) == 0 ) pcangle_in(i) = cos(thetam(nth,np,nr,nsa_in))
            if ( mode(2) == 1 .and. mode(3) == 0 ) pcangle_in(i) = cos(thetam_pg(nth,np,nr,nsa_in))
            if ( mode(2) == 0 .and. mode(3) == 1 ) pcangle_in(i) = cos(thetam_rg(nth,np,nr,nsa_in))            
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

          if ( pcangle_in(i)*aefp(nsa_in) >= 0.d0 ) then
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
    use fowcomm
    use fpcomm, only : npmax, nthmax, nrmax, nsamax, rkind

    integer,intent(in) :: mode(3), nobt_in(:,:,:), nsa_in
    type(orbit),intent(inout) :: orbit_in(:,:,:,:)
    real(rkind), allocatable :: construct_input(:,:)
    integer :: nsa,nr,nth,np
    integer :: i, j

    allocate(construct_input(5,nstp_max+1))
    ! time_ob(0:nstp_max,nobt_max)
    do nr = 1, nrmax+mode(3)
      do np = 1, npmax+mode(2)
        do nth = 1, nthmax+mode(1)
          ! exclude forbitten region
          if ( mode(1) == 0 .and. nth == nth_forbitten(nsa_in) ) then
            orbit_in(nth,np,nr,nsa_in) = construct_orbit(0, construct_input)
            cycle
          end if
          if ( mode(2) == 1 .and. np == 1 ) then
            orbit_in(nth,np,nr,nsa_in) = construct_orbit(0, construct_input)
            cycle
          end if
          if ( mode(3) == 1 .and. nr == 1 ) then
            orbit_in(nth,np,nr,nsa_in) = construct_orbit(0, construct_input)
            cycle
          end if

          i = nobt_in(nth,np,nr)

          do j = 1, nstp_max_nobt(i)+1
            construct_input(1,j)  = time_ob(j-1,i)
            construct_input(2,j)  = psip_ob(j-1,i)
            construct_input(3,j)  = babs_ob(j-1,i)
            construct_input(4,j)  = acos( vpara_ob(j-1,i)/sqrt( vperp_ob(j-1,i)**2+vpara_ob(j-1,i)**2) )
            construct_input(5,j)  = theta_ob(j-1,i)
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
    integer :: nr_max, nr_min, ierr, nr_l, mode(3)

    nr_max = max(nr_in+1, 5)
    nr_max = min(nr_max, nrmax+1)

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

    if ( rr_l <= 1.d0 ) then
      F_ret = Fpsig(1)
    else if ( rr_l >= (nrmax+1)*1.d0 ) then
      F_ret = Fpsig(nrmax+1)
    else 
      nr_l = int(rr_l)
      F_ret = Fpsig(nr_l+1)+(Fpsig(nr_l+1)-Fpsig(nr_l))*(rr_l-dble(nr_l+1))
    end if

    return

  end function func_orbit_F

  function func_orbit_mean_ra(orbit_in, U_rm) result(ra_ret)
    ! calculate mean minor radius along an orbit
    ! U_rm : spline coefficient,  RM(psim)
    use fpcomm
    use fowcomm

    implicit none 
    real(rkind) :: ra_ret
    type(orbit),intent(in) :: orbit_in
    real(rkind),dimension(4,nrmax),intent(in) :: U_rm(:,:)
    integer :: nstp, nstpmax, ierr = 0
    real(rkind) :: sum_psip, mean_psip
 
    nstpmax = orbit_in%nstp_max

    sum_psip = 0.d0

    do nstp = 1, nstpmax
      sum_psip = sum_psip + orbit_in%psip(nstp)
    end do

    mean_psip = sum_psip/dble(nstpmax)

    call SPL1DF(mean_psip, ra_ret, psim, U_rm, nrmax, ierr)

  end function

  function construct_orbit(n,input) result(ret)
    use fowcomm
    implicit none
    type(orbit) :: ret
    integer, intent(in) :: n
    real(rkind), intent(in) :: input(:,:)
    integer :: i

    allocate(ret%time(n+1), ret%psip(n+1), ret%Babs(n+1), ret%theta(n+1), ret%thetap(n+1))

    if ( n == 0 ) then
      ret%nstp_max=n+1
      ret%time(1) = 0.d0
      ret%psip(1) = 0.d0
      ret%Babs(1) = 0.d0
      ret%theta(1) = 0.d0
      ret%thetap(1) = 0.d0
    else
      ret%nstp_max=n+1
      do i=1,n+1
        ret%time(i) = input(1,i) 
        ret%psip(i) = input(2,i)
        ret%Babs(i) = input(3,i)
        ret%theta(i) = input(4,i)
        ret%thetap(i) = input(5,i)
      end do
    end if

  end function construct_orbit

  subroutine destruct_orbit(orbit_in)
    use fowcomm
    type(orbit) :: orbit_in
    deallocate(orbit_in%time, orbit_in%psip, orbit_in%babs, orbit_in%theta)
  end subroutine destruct_orbit

  subroutine fow_orbit_save
    use fpcomm
    use fowcomm
    implicit none
    character(30) :: BIN_DIR, filename
    integer :: nm ,nth, np, nr, nsa, nstp, nstpmax, nmmmax, nmpmax,nmtmax, nmrmax&
            , nstp_sum_m, nstp_sum_p, nstp_sum_t, nstp_sum_r, i, len
    integer,allocatable,dimension(:,:,:,:) :: nma_m, nma_p, nma_r, nma_t
    real(rkind),allocatable,dimension(:) :: time_m, psip_m, theta_m, thetap_m, Babs_m, nstpmax_m
    real(rkind),allocatable,dimension(:) :: time_p, psip_p, theta_p, thetap_p, Babs_p, nstpmax_p
    real(rkind),allocatable,dimension(:) :: time_r, psip_r, theta_r, thetap_r, Babs_r, nstpmax_r
    real(rkind),allocatable,dimension(:) :: time_t, psip_t, theta_t, thetap_t, Babs_t, nstpmax_t

    BIN_DIR = "../fp.ota/bin/"

    filename = TRIM(BIN_DIR)//"fpparm.dat"
    open(10, file=filename, access='direct', recl=rkind*7+4*4, form='unformatted' )
    write(10, rec=1)RR,RA,RKAP,RDLT,RB,BB,RIP,nthmax,npmax,nrmax,nsamax
    close(10)

    allocate(nma_m(nthmax,npmax,nrmax,nsamax))
    allocate(nma_p(nthmax,npmax+1,nrmax,nsamax))
    allocate(nma_r(nthmax,npmax,nrmax+1,nsamax))
    allocate(nma_t(nthmax+1,npmax,nrmax,nsamax))
    nmmmax = nthmax*npmax*nrmax*nsamax
    nmpmax = nthmax*(npmax+1)*nrmax*nsamax
    nmtmax = (nthmax+1)*npmax*nrmax*nsamax
    nmrmax = nthmax*npmax*(nrmax+1)*nsamax
    allocate(nstpmax_m(nmmmax))
    allocate(nstpmax_p(nmpmax))
    allocate(nstpmax_t(nmtmax))
    allocate(nstpmax_r(nmrmax))

    nstp_sum_m = 0
    nstp_sum_p = 0
    nstp_sum_t = 0
    nstp_sum_r = 0    
    do nsa = 1, nsamax
      do nr = 1, nrmax+1
        do np = 1, npmax+1
          do nth = 1, nthmax+1
            if ( nth <= nthmax .and. np <= npmax .and. nr <= nrmax ) then
              nm = (nsa-1)*nthmax*npmax*nrmax&
                  +(nr-1)*nthmax*npmax&
                  +(np-1)*nthmax + nth
              nma_m(nth,np,nr,nsa) = nm
              nstpmax_m(nm) = orbit_m(nth,np,nr,nsa)%nstp_max
              nstp_sum_m = nstp_sum_m + orbit_m(nth,np,nr,nsa)%nstp_max
            end if
            if ( np <= npmax .and. nr <= nrmax ) then
              nm = (nsa-1)*nthmax*npmax*nrmax&
                  +(nr-1)*nthmax*npmax&
                  +(np-1)*nthmax + nth
              nma_t(nth,np,nr,nsa) = nm
              nstpmax_t(nm) = orbit_th(nth,np,nr,nsa)%nstp_max
              nstp_sum_t = nstp_sum_t + orbit_th(nth,np,nr,nsa)%nstp_max
            end if
            if ( nth <= nthmax .and. nr <= nrmax ) then
              nm = (nsa-1)*nthmax*(npmax+1)*nrmax&
                  +(nr-1)*nthmax*(npmax+1)&
                  +(np-1)*nthmax + nth
              nma_p(nth,np,nr,nsa) = nm
              nstpmax_p(nm) = orbit_p(nth,np,nr,nsa)%nstp_max
              nstp_sum_p = nstp_sum_p + orbit_p(nth,np,nr,nsa)%nstp_max
            end if
            if ( np <= npmax .and. nth <= nthmax ) then
              nm = (nsa-1)*nthmax*npmax*(nrmax+1)&
                  +(nr-1)*nthmax*npmax&
                  +(np-1)*nthmax + nth
              nma_r(nth,np,nr,nsa) = nm
              nstpmax_r(nm) = orbit_r(nth,np,nr,nsa)%nstp_max
              nstp_sum_r = nstp_sum_r + orbit_r(nth,np,nr,nsa)%nstp_max
            end if
          end do
        end do
      end do
    end do

    filename = TRIM(BIN_DIR)//"ob_index.dat"
    open(12, file=filename, access='direct', recl=4*8, form='unformatted' )
    write(12,rec=1)nmmmax,nmpmax,nmtmax,nmrmax,nstp_sum_m,nstp_sum_p,nstp_sum_t,nstp_sum_r
    close(12)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_m.dat"
    open(20, file=filename, access='direct', recl=rkind*nmmmax, form='unformatted' )
    write(20,rec=1)nstpmax_m
    close(20)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_p.dat"
    open(30, file=filename, access='direct', recl=rkind*nmpmax, form='unformatted' )
    write(30,rec=1)nstpmax_p
    close(30)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_t.dat"
    open(40, file=filename, access='direct', recl=rkind*nmtmax, form='unformatted' )
    write(40,rec=1)nstpmax_t
    close(40)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_r.dat"
    open(50, file=filename, access='direct', recl=rkind*nmrmax, form='unformatted' )
    write(50,rec=1)nstpmax_r
    close(50)

    allocate(time_m(nstp_sum_m), psip_m(nstp_sum_m))
    allocate(theta_m(nstp_sum_m), thetap_m(nstp_sum_m), Babs_m(nstp_sum_m))

    allocate(time_p(nstp_sum_p), psip_p(nstp_sum_p))
    allocate(theta_p(nstp_sum_p), thetap_p(nstp_sum_p), Babs_p(nstp_sum_p))

    allocate(time_t(nstp_sum_t), psip_t(nstp_sum_t))
    allocate(theta_t(nstp_sum_t), thetap_t(nstp_sum_t), Babs_t(nstp_sum_t))

    allocate(time_r(nstp_sum_r), psip_r(nstp_sum_r))
    allocate(theta_r(nstp_sum_r), thetap_r(nstp_sum_r), Babs_r(nstp_sum_r))


    do nsa = 1, nsamax
      do nr = 1, nrmax+1
        do np = 1, npmax+1
          do nth = 1, nthmax+1
            if ( np <= npmax .and. nth <= nthmax .and. nr <= nrmax ) then
              nstpmax = orbit_m(nth,np,nr,nsa)%nstp_max
              do nstp = 1, nstpmax
                i = (nsa-1)*nstpmax*nthmax*npmax*nrmax&
                    +(nr-1)*nstpmax*nthmax*npmax&
                    +(np-1)*nstpmax*nthmax&
                    +(nth-1)*nstpmax + nstp
                time_m(i) = orbit_m(nth,np,nr,nsa)%time(nstp)
                psip_m(i) = orbit_m(nth,np,nr,nsa)%psip(nstp)
                Babs_m(i) = orbit_m(nth,np,nr,nsa)%Babs(nstp)
                theta_m(i) = orbit_m(nth,np,nr,nsa)%theta(nstp)
                thetap_m(i) = orbit_m(nth,np,nr,nsa)%thetap(nstp)
              end do  
            end if
            if ( nth <= nthmax .and. nr <= nrmax ) then
              nstpmax = orbit_p(nth,np,nr,nsa)%nstp_max
              do nstp = 1, nstpmax
                i = (nsa-1)*nstpmax*nthmax*npmax*nrmax&
                    +(nr-1)*nstpmax*nthmax*npmax&
                    +(np-1)*nstpmax*nthmax&
                    +(nth-1)*nstpmax + nstp
                time_p(i) = orbit_p(nth,np,nr,nsa)%time(nstp)
                psip_p(i) = orbit_p(nth,np,nr,nsa)%psip(nstp)
                Babs_p(i) = orbit_p(nth,np,nr,nsa)%Babs(nstp)
                theta_p(i) = orbit_p(nth,np,nr,nsa)%theta(nstp)
                thetap_p(i) = orbit_p(nth,np,nr,nsa)%thetap(nstp)
              end do
            end if
            if ( np <= npmax .and. nr <= nrmax ) then
              nstpmax = orbit_th(nth,np,nr,nsa)%nstp_max
              do nstp = 1, nstpmax
                i = (nsa-1)*nstpmax*nthmax*npmax*nrmax&
                    +(nr-1)*nstpmax*nthmax*npmax&
                    +(np-1)*nstpmax*nthmax&
                    +(nth-1)*nstpmax + nstp
                time_t(i) = orbit_th(nth,np,nr,nsa)%time(nstp)
                psip_t(i) = orbit_th(nth,np,nr,nsa)%psip(nstp)
                Babs_t(i) = orbit_th(nth,np,nr,nsa)%Babs(nstp)
                theta_t(i) = orbit_th(nth,np,nr,nsa)%theta(nstp)
                thetap_t(i) = orbit_th(nth,np,nr,nsa)%thetap(nstp)
              end do  
            end if
            if ( np <= npmax .and. nth <= nthmax ) then
              nstpmax = orbit_r(nth,np,nr,nsa)%nstp_max
              do nstp = 1, nstpmax
                i = (nsa-1)*nstpmax*nthmax*npmax*nrmax&
                    +(nr-1)*nstpmax*nthmax*npmax&
                    +(np-1)*nstpmax*nthmax&
                    +(nth-1)*nstpmax + nstp
                time_r(i) = orbit_r(nth,np,nr,nsa)%time(nstp)
                psip_r(i) = orbit_r(nth,np,nr,nsa)%psip(nstp)
                Babs_r(i) = orbit_r(nth,np,nr,nsa)%Babs(nstp)
                theta_r(i) = orbit_r(nth,np,nr,nsa)%theta(nstp)
                thetap_r(i) = orbit_r(nth,np,nr,nsa)%thetap(nstp)
              end do  
            end if
    
          end do
        end do
      end do
    end do

    len = rkind*(nstp_sum_m+nstp_sum_p+nstp_sum_t+nstp_sum_r)*5 

    filename = TRIM(BIN_DIR)//"ob_quantities.dat"
    open(60, file=filename, access='direct', recl=rkind*nstp_sum_m*5, form='unformatted' )
    write(60,rec=1)time_m, psip_m, theta_m, thetap_m, Babs_m
    write(60,rec=rkind*nstp_sum_m*5+1)time_p, psip_p, theta_p, thetap_p, Babs_p
    write(60,rec=rkind*(nstp_sum_m+nstp_sum_p)*5+1)time_t, psip_t, theta_t, thetap_t, Babs_t
    write(60,rec=rkind*(nstp_sum_m+nstp_sum_p+nstp_sum_t)*5+1)time_r, psip_r, theta_r, thetap_r, Babs_r
    close(60)

  end subroutine fow_orbit_save

  subroutine fow_orbit_load(flag)
    use fpcomm
    use fowcomm
    implicit none
    integer,intent(out) :: flag
    character(30) :: BIN_DIR, eqin, mesh, filename
    integer :: nm ,nth, np, nr, nsa, nstp, nstpmax, nmmmax, nmpmax,nmtmax, nmrmax&
            , nstp_sum_m, nstp_sum_p, nstp_sum_t, nstp_sum_r, i
    integer,allocatable,dimension(:,:,:,:) :: nma_m, nma_p, nma_r, nma_t
    real(rkind),allocatable,dimension(:) :: time_m, psip_m, theta_m, thetap_m, Babs_m, nstpmax_m
    real(rkind),allocatable,dimension(:) :: time_p, psip_p, theta_p, thetap_p, Babs_p, nstpmax_p
    real(rkind),allocatable,dimension(:) :: time_r, psip_r, theta_r, thetap_r, Babs_r, nstpmax_r
    real(rkind),allocatable,dimension(:) :: time_t, psip_t, theta_t, thetap_t, Babs_t, nstpmax_t

    real(rkind) :: RR_,RA_,RKAP_,RDLT_,RB_,BB_,RIP_
    integer :: nthm_,npm_,nrm_,nsam_

    flag = 0

    BIN_DIR = "../fp.ota/bin"

    filename = TRIM(BIN_DIR)//"fpparm.dat"
    open(10, file=filename, access='direct', recl=rkind*7+4*4, form='unformatted' )
    read(10,rec=1) RR_,RR_,RA_,RKAP_,RDLT_,RB_,BB_,RIP_
    read(10,rec=rkind*7+1)nthm_,npm_,nrm_,nsam_
    close(10)

    if ( RR /= RR_ .or. RA /= RA_ .or. RKAP /= RKAP_ &
        .or. RDLT /= RDLT_ .or. RB /= RB_ .or. BB /= BB_ .or. RIP /= RIP_) then
      flag = 1
      return
    end if
    if ( nthmax /= nthm_ .or. npmax /= npm_ &
        .or. nrmax /= nrm_ .or. nsamax /= nsam_ ) then
      flag = 2
      return    
    end if
    
    filename = TRIM(BIN_DIR)//"ob_index.dat"
    open(11, file=filename, access='direct', recl=4*8, form='unformatted' )
    read(11,rec=1)nmmmax,nmpmax,nmtmax,nmrmax,nstp_sum_m,nstp_sum_p,nstp_sum_t,nstp_sum_r
    close(11)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_m.dat"
    open(20, file=filename, access='direct', recl=rkind*nmmmax, form='unformatted' )
    read(20,rec=1)nstpmax_m
    close(20)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_p.dat"
    open(30, file=filename, access='direct', recl=rkind*nmpmax, form='unformatted' )
    read(30,rec=1)nstpmax_p
    close(30)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_t.dat"
    open(40, file=filename, access='direct', recl=rkind*nmtmax, form='unformatted' )
    read(40,rec=1)nstpmax_t
    close(40)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_r.dat"
    open(50, file=filename, access='direct', recl=rkind*nmrmax, form='unformatted' )
    read(50,rec=1)nstpmax_r
    close(50)


  end subroutine fow_orbit_load

  subroutine leastSquareMethodForQuadric(y,nxmax,nxstart,nxend,reval)
    ! y=reval(1)*nx^2+reval(2)*nx+reval(3)
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
  
end module foworbit
