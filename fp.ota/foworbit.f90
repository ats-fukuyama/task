module foworbit
  
  use obcomm
  use obprep
  private
  public :: fow_orbit, fow_cal_local_COMs, fow_set_obparm
  public :: quantities_at_Bminimum, mean_ra_quantities

  double precision,allocatable :: UR(:,:), Ups(:,:), UF(:,:), UB(:,:,:,:), dBdrt(:,:)

contains

  subroutine fow_set_obparm(ierr)

    use obparm
    use obinit
    use fowcomm
    use fpcomm

    implicit none
    integer,intent(out) :: ierr
    logical,save :: called = .false.
    integer :: nthp
    real(rkind),allocatable :: dradpsi(:)

    ierr = 0

    if ( called ) then
      return
    else
      called = .true.
    end if

    call ob_init
    call ob_parm(1,'../fp.ota/fpparm',ierr)
    nobt_max = 1
    call ob_prep(ierr)
    call ob_allocate

    allocate(Ups(4,nrmax))
    allocate(UF(4,nrmax))
    allocate(UB(4,4,nrmax,nthpmax))
    allocate(dBdrt(nrmax,nthpmax))
    allocate(dradpsi(nrmax))
    allocate(UR(4,nrmax))

    call spl1D(rm,psim,dpsimdr,Ups,nrmax,3,ierr)
    call spl1D(rm,Fpsi,dFdr,UF,nrmax,3,ierr)

    do nthp = 1, nthpmax
      call first_order_derivative(dBdrt(:,nthp),dBdthp(:,nthp),rm)
    end do
    call spl2D(rm,theta_p,Babs,dBdr,dBdthp,dBdrt,UB,nrmax,nrmax,nthpmax,3,3,ierr)

    call first_order_derivative(dradpsi, rm, psim)
    call SPL1D(psim,rm,dradpsi,UR,nrmax,3,ierr)

  end subroutine fow_set_obparm

  subroutine fow_orbit(flag, ierr)
    use fpcomm
    use fowcomm
    use obparm
    use obinit

    implicit none
    integer,intent(out) :: flag, ierr
    integer :: nth,np,nr,nsa,n,nmax
    real(rkind) :: end_time, begin_time, thetap_in

    ierr = 0
    flag = 0

    if ( model_obload >= 1 ) then
      call cpu_time(begin_time)
      call load_orbit(ierr)
      call cpu_time(end_time)
      write(6,'("Load orbit time :",ES10.3,"[sec]")')end_time-begin_time

      if ( ierr == 0 ) then
        return
      else
        write(6,'("IOSTAT = ",I4," in load_orbit")'),ierr
      end if

    end if

    ierr = 0
    call fow_set_obparm(ierr)

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax

            if ( aefp(nsa) * COS( thetam(nth,np,nr,nsa) ) >= 0.d0 ) then
              thetap_in = 0.d0
            else 
              thetap_in = pi
            end if

            call construct_orbit(orbit_m(nth,np,nr,nsa), pm(np,nsa)*ptfp0(nsa)&
                                ,thetam(nth,np,nr,nsa), thetap_in, psim(nr), nsa, ierr)

          end do
        end do
      end do
    end do

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax+1
          do nth = 1, nthmax

            if ( np == 1 ) then
              call construct_orbit_zero(orbit_p(nth,np,nr,nsa))

            else
              if ( aefp(nsa) * COS( thetam_pg(nth,np,nr,nsa) ) >= 0.d0 ) then
                thetap_in = 0.d0
              else 
                thetap_in = pi
              end if
  
              call construct_orbit(orbit_p(nth,np,nr,nsa), pg(np,nsa)*ptfp0(nsa)&
                                  ,thetam_pg(nth,np,nr,nsa), thetap_in, psim(nr), nsa, ierr)

            end if

          end do
        end do
      end do
    end do

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax+1

            if ( aefp(nsa) * COS( thetamg(nth,np,nr,nsa) ) >= 0.d0 ) then
              thetap_in = 0.d0
            else 
              thetap_in = pi
            end if

            call construct_orbit(orbit_th(nth,np,nr,nsa), pm(np,nsa)*ptfp0(nsa)&
                                ,thetamg(nth,np,nr,nsa), thetap_in, psim(nr), nsa, ierr)

          end do
        end do
      end do
    end do

    do nsa = 1, nsamax
      do nr = 1, nrmax+1
        do np = 1, npmax
          do nth = 1, nthmax
            if ( nr == 1 ) then
              call construct_orbit_zero(orbit_r(nth,np,nr,nsa))

            else
              if ( aefp(nsa) * COS( thetam_rg(nth,np,nr,nsa) ) >= 0.d0 ) then
                thetap_in = 0.d0
              else 
                thetap_in = pi
              end if
  
              call construct_orbit(orbit_r(nth,np,nr,nsa), pm(np,nsa)*ptfp0(nsa)&
                                  ,thetam_rg(nth,np,nr,nsa), thetap_in, psimg(nr), nsa, ierr)
  
            end if

          end do
        end do
      end do
    end do

    if ( model_obload /= 0 ) then
      call cpu_time(begin_time)
      call save_orbit(ierr)
      call cpu_time(end_time)
      write(6,'("TASK/OB time:",ES10.3,"[sec]")')end_time-begin_time
      if ( ierr /= 0 ) then
        write(6,'("IOSTAT = ",I4," in save_orbit")'),ierr
      end if
    end if

    flag = 1

  end subroutine fow_orbit

  subroutine fow_cal_local_COMs(thetaml, rhoml, tau_loss, momentum, pitch_angle, theta_pol, psi_pol, ns)
    use fpcomm
    use fowcomm
    implicit none
    real(rkind),intent(out) :: thetaml, rhoml, tau_loss
    real(rkind),intent(in) :: momentum, pitch_angle, theta_pol, psi_pol
    type(orbit) :: ob
    integer,intent(in) :: ns
    integer :: nstp, nstpmax, loc_max,i, im, ierr, ir, il
    real(rkind) :: psip_max, diff_thetap_min,  diff_thetap, thetap_eq(6), A

    call construct_orbit(ob, momentum, pitch_angle, theta_pol, psi_pol, ns, ierr)

    nstpmax = ob%nstp_max 

    if ( ob%psip(nstpmax) > psi0 ) then ! loss to wall
      rhoml = ob%psip(nstpmax)
      thetaml = ACOS( ob%costh(nstpmax) )
      tau_loss = ob%time(nstpmax)
      return
    end if
    
    psip_max = -1.d0
    do nstp = 1, nstpmax
      if ( psip_max <= ob%psip(nstp) ) then
        psip_max = ob%psip(nstp)
        loc_max  = nstp
      end if
    end do

    thetap_eq = [-2.d0*pi, -1.d0*pi, 0.d0, pi, 2.d0*pi, 3.d0*pi] ! opion of thetap at psim

    ! thetap_eq(im) is thetap(psim)
    diff_thetap_min = 1.d8
    do i = 1, 6
      diff_thetap = ABS( thetap_eq(i)-ob%thetap(loc_max) )
      if ( diff_thetap < diff_thetap_min ) then
        diff_thetap_min = diff_thetap
        im = i
      end if
    end do

    ! if ( diff_thetap_min <= pi/1.d2 ) then
      rhoml  = ob%r(loc_max)
      thetaml= ACOS( ob%costh(loc_max) )
    ! else
    ! end if

    tau_loss = 0.d0

  end subroutine fow_cal_local_COMs

  subroutine construct_orbit(ob, momentum, pitch_angle, theta_pol, psi_pol, ns, ierr)
    use fpcomm
    use fowcomm

    use obcalc

    implicit none
    type(orbit),intent(out) :: ob
    real(rkind),intent(in) :: momentum, pitch_angle, theta_pol, psi_pol
    integer,intent(out) :: ierr
    integer,intent(in) :: ns
    integer :: nstp, nstpmax, nthp
    real(rkind) :: pv, mass0, dummy
    logical :: is_allocated
    integer :: A,B,i

    ierr = 0

    mass0 = pa(ns)*amp
    pv = SQRT( 1.d0+momentum**2/(mass0*vc)**2 )
    ns_ob = ns

    penergy_in(1) = mass0*vc**2*(pv-1.d0)/(aee*1.d3)
    pcangle_in(1) = COS( pitch_angle )
    theta_in(1)   = theta_pol!*180.d0/pi
    psipn_in(1)   = psi_pol/psi0
    zeta_in(1)    = 0.d0

    call ob_calc(ierr)

    nstpmax = nstp_max_nobt(1)+1

    if ( nstpmax <= max_stp ) then
      ob%nstp_max = nstpmax
    else
      ob%nstp_max = max_stp
      A = nstpmax / (max_stp-1)
      B = MOD(nstpmax, max_stp-1)  
    end if

    allocate(ob%time(ob%nstp_max))
    allocate(ob%psip(ob%nstp_max))
    allocate(ob%Babs(ob%nstp_max))
    allocate(ob%costh(ob%nstp_max))
    allocate(ob%sinth(ob%nstp_max))
    allocate(ob%thetap(ob%nstp_max))
    allocate(ob%F(ob%nstp_max))
    allocate(ob%r(ob%nstp_max))
    allocate(ob%dpsipdr(ob%nstp_max))
    allocate(ob%dFdr(ob%nstp_max))
    allocate(ob%dBdr(ob%nstp_max))
    allocate(ob%dBdthp(ob%nstp_max))

    do nstp = 1, ob%nstp_max
      if ( nstpmax <= max_stp ) then
        i = nstp-1
      else
        if ( nstp == 1 ) then
          i = 0
        else if ( 2 <= nstp .and. nstp <= max_stp-B ) then
          i = A * (nstp-1)-1
        else
          i = A*(max_stp-B-1)+(A+1)*(nstp-max_stp+B)-1
        end if
      end if

      ob%time(nstp)  = time_ob(i,1)
      ob%psip(nstp)  = psip_ob(i,1)
      ob%Babs(nstp)  = Babs_ob(i,1)
      ob%thetap(nstp)= theta_ob(i,1)

      call spl1DF(psip_ob(i,1),ob%r(nstp),psim,UR,nrmax,ierr)

      ob%costh(nstp) = vpara_ob(i,1)/SQRT( vpara_ob(i,1)**2+vperp_ob(i,1)**2 )
      ob%sinth(nstp) = SQRT( 1.d0-ob%costh(nstp)**2 )

      call spl1DD(ob%r(nstp),dummy,ob%dpsipdr(nstp),rm,Ups,nrmax,ierr)
      call spl1DD(ob%r(nstp),ob%F(nstp),ob%dFdr(nstp),rm,UF,nrmax,ierr)
      call spl2DD(ob%r(nstp),theta_ob(i,1),dummy,ob%dBdr(nstp),ob%dBdthp(nstp),rm,theta_p,UB,nrmax,nrmax,nthpmax,ierr)

    end do

    return

  end subroutine construct_orbit

  subroutine construct_orbit_zero(ob)
    use fpcomm
    use fowcomm
    implicit none
    type(orbit),intent(out) :: ob
    
    ob%nstp_max = 1
    allocate(ob%psip(1))
    allocate(ob%thetap(1))
    allocate(ob%Babs(1))
    allocate(ob%time(1))
    allocate(ob%sinth(1))
    allocate(ob%costh(1))
    allocate(ob%F(1))
    allocate(ob%r(1))
    allocate(ob%dpsipdr(1))
    allocate(ob%dFdr(1))
    allocate(ob%dBdr(1))
    allocate(ob%dBdthp(1))

    ob%psip(1)     = 0.d0
    ob%thetap(1)   = 0.d0
    ob%Babs(1)     = 0.d0
    ob%time(1)     = 0.d0
    ob%sinth(1)    = 0.d0
    ob%costh(1)    = 0.d0
    ob%F(1)        = 0.d0
    ob%r(1)        = 0.d0
    ob%dpsipdr(1)  = 0.d0
    ob%dFdr(1)     = 0.d0
    ob%dBdr(1)     = 0.d0
    ob%dBdthp(1)   = 0.d0

  end subroutine construct_orbit_zero

  subroutine save_orbit(ierr)
    use fpcomm
    use fowcomm
    implicit none
    integer,intent(out) :: ierr
    character(30) :: BIN_DIR, filename
    integer :: nm ,nth, np, nr, nsa, nstp, nstpmax, it, ip, ir, im, i

    ierr = 0

    BIN_DIR = "bin/"

    filename = TRIM(BIN_DIR)//"eqparm.bin"
    open(10,file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)
    write(10,iostat=ierr)RR,RA,RKAP,RDLT,RB,BB,RIP
    close(10)

    filename = TRIM(BIN_DIR)//"fpparm.bin"
    open(11,file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)
    write(11,iostat=ierr)nthmax,npmax,nrmax,nsamax,nthpmax
    close(11)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_m.bin"
    open(50,file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_t.bin"
    open(51,file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_p.bin"
    open(52,file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_r.bin"
    open(53,file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_m.bin"
    open(60,file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_t.bin"
    open(61,file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_p.bin"
    open(62,file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_r.bin"
    open(63,file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)

    do nsa = 1, nsamax
      do nr = 1, nrmax+1
        do np = 1, npmax+1
          do nth = 1, nthmax+1

            if ( nth <= nthmax .and. np <= npmax .and. nr <= nrmax ) then
              write(50,iostat=ierr)orbit_m(nth,np,nr,nsa)%nstp_max
              nstpmax = orbit_m(nth,np,nr,nsa)%nstp_max

              do nstp = 1, nstpmax
                write(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%time(nstp)
                write(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%psip(nstp)
                write(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%thetap(nstp)
                write(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%Babs(nstp)
                write(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%costh(nstp)
                write(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%sinth(nstp)
                write(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%F(nstp)
                write(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%r(nstp)
                write(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%dpsipdr(nstp)
                write(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%dFdr(nstp)
                write(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%dBdr(nstp)
                write(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%dBdthp(nstp)
              end do
              
            end if

            if ( np <= npmax .and. nr <= nrmax ) then
              write(51,iostat=ierr)orbit_th(nth,np,nr,nsa)%nstp_max
              nstpmax = orbit_th(nth,np,nr,nsa)%nstp_max

              do nstp = 1, nstpmax
                write(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%time(nstp)
                write(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%psip(nstp)
                write(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%thetap(nstp)
                write(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%Babs(nstp)
                write(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%costh(nstp)
                write(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%sinth(nstp)
                write(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%F(nstp)
                write(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%r(nstp)
                write(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%dpsipdr(nstp)
                write(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%dFdr(nstp)
                write(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%dBdr(nstp)
                write(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%dBdthp(nstp)
              end do
            end if
    
            if ( nth <= nthmax .and. nr <= nrmax ) then
              write(52,iostat=ierr)orbit_p(nth,np,nr,nsa)%nstp_max
              nstpmax = orbit_p(nth,np,nr,nsa)%nstp_max

              do nstp = 1, nstpmax
                write(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%time(nstp)
                write(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%psip(nstp)
                write(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%thetap(nstp)
                write(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%Babs(nstp)
                write(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%costh(nstp)
                write(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%sinth(nstp)
                write(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%F(nstp)
                write(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%r(nstp)
                write(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%dpsipdr(nstp)
                write(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%dFdr(nstp)
                write(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%dBdr(nstp)
                write(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%dBdthp(nstp)
              end do
            end if

            if ( nth <= nthmax .and. np <= npmax ) then
              write(53,iostat=ierr)orbit_r(nth,np,nr,nsa)%nstp_max
              nstpmax = orbit_r(nth,np,nr,nsa)%nstp_max

              do nstp = 1, nstpmax
                write(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%time(nstp)
                write(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%psip(nstp)
                write(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%thetap(nstp)
                write(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%Babs(nstp)
                write(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%costh(nstp)
                write(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%sinth(nstp)
                write(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%F(nstp)
                write(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%r(nstp)
                write(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%dpsipdr(nstp)
                write(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%dFdr(nstp)
                write(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%dBdr(nstp)
                write(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%dBdthp(nstp)
              end do
            end if

          end do
        end do
      end do
    end do

    do i = 0, 3
      close(50+i)
      close(60+i)
    end do

  end subroutine save_orbit

  subroutine load_orbit(ierr)
    use fpcomm
    use fowcomm
    implicit none
    integer,intent(out) :: ierr
    character(30) :: BIN_DIR, eqin, mesh, filename
    integer :: nm ,nth, np, nr, nsa, nstp, nstpmax, im, ip, it, ir, i

    real(rkind) :: RR_,RA_,RKAP_,RDLT_,RB_,BB_,RIP_
    integer :: nthm_,npm_,nrm_,nsam_,nthpm_
    integer :: access

    ierr = 0

    call system('mkdir -p bin')

    BIN_DIR = "bin/"

    if ( access( TRIM(BIN_DIR)//"fpparm.bin", " ")       /= 0 &
    .or. access( TRIM(BIN_DIR)//"eqparm.bin", " ")       /= 0 &
    .or. access( TRIM(BIN_DIR)//"ob_nstpmax_m.bin", " ") /= 0 &
    .or. access( TRIM(BIN_DIR)//"ob_nstpmax_p.bin", " ") /= 0 &
    .or. access( TRIM(BIN_DIR)//"ob_nstpmax_t.bin", " ") /= 0 &
    .or. access( TRIM(BIN_DIR)//"ob_nstpmax_r.bin", " ") /= 0 &
    .or. access( TRIM(BIN_DIR)//"ob_m.bin", " ")         /= 0 &
    .or. access( TRIM(BIN_DIR)//"ob_p.bin", " ")         /= 0 &
    .or. access( TRIM(BIN_DIR)//"ob_t.bin", " ")         /= 0 &
    .or. access( TRIM(BIN_DIR)//"ob_r.bin", " ")         /= 0 &
    ) then
      ierr = 1000
      return
    end if

    filename = TRIM(BIN_DIR)//"eqparm.bin"
    open(10,file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)
    read(10,iostat=ierr)RR_,RA_,RKAP_,RDLT_,RB_,BB_,RIP_
    close(10)

    filename = TRIM(BIN_DIR)//"fpparm.bin"
    open(11,file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)
    read(11,iostat=ierr)nthm_,npm_,nrm_,nsam_,nthpm_
    close(11)

    if ( RR   /= RR_   &
    .or. RA   /= RA_   &
    .or. RKAP /= RKAP_ &
    .or. RDLT /= RDLT_ &
    .or. RB   /= RB_   &
    .or. BB   /= BB_   &
    .or. RIP  /= RIP_  &
    ) then
      ierr = 1001
    end if
    
    if ( nthmax  /= nthm_  &
    .or. npmax   /= npm_   &
    .or. nrmax   /= nrm_   &
    .or. nsamax  /= nsam_  &
    ) then
      ierr = 1002
    end if

    if ( ierr /= 0 ) return

    filename = TRIM(BIN_DIR)//"ob_nstpmax_m.bin"
    open(50,file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_t.bin"
    open(51,file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_p.bin"
    open(52,file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_r.bin"
    open(53,file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_m.bin"
    open(60,file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_t.bin"
    open(61,file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_p.bin"
    open(62,file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_r.bin"
    open(63,file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)

    do nsa = 1, nsamax
      do nr = 1, nrmax+1
        do np = 1, npmax+1
          do nth = 1, nthmax+1
            if ( nth <= nthmax .and. np <= npmax .and. nr <= nrmax ) then
              read(50,iostat=ierr)nstpmax
              orbit_m(nth,np,nr,nsa)%nstp_max = nstpmax

              allocate(orbit_m(nth,np,nr,nsa)%time(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%psip(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%Babs(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%thetap(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%costh(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%sinth(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%F(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%r(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%dpsipdr(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%dFdr(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%dBdr(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%dBdthp(nstpmax))

              do nstp = 1, nstpmax
                read(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%time(nstp)
                read(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%psip(nstp)
                read(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%thetap(nstp)
                read(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%Babs(nstp)
                read(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%costh(nstp)
                read(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%sinth(nstp)
                read(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%F(nstp)
                read(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%r(nstp)
                read(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%dpsipdr(nstp)
                read(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%dFdr(nstp)
                read(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%dBdr(nstp)
                read(60,iostat=ierr)orbit_m(nth,np,nr,nsa)%dBdthp(nstp)
              end do
              
            end if

            if ( np <= npmax .and. nr <= nrmax ) then
              read(51,iostat=ierr)nstpmax
              orbit_th(nth,np,nr,nsa)%nstp_max = nstpmax

              allocate(orbit_th(nth,np,nr,nsa)%time(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%psip(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%Babs(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%thetap(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%costh(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%sinth(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%F(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%r(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%dpsipdr(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%dFdr(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%dBdr(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%dBdthp(nstpmax))


              do nstp = 1, nstpmax
                read(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%time(nstp)
                read(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%psip(nstp)
                read(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%thetap(nstp)
                read(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%Babs(nstp)
                read(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%costh(nstp)
                read(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%sinth(nstp)
                read(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%F(nstp)
                read(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%r(nstp)
                read(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%dpsipdr(nstp)
                read(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%dFdr(nstp)
                read(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%dBdr(nstp)
                read(61,iostat=ierr)orbit_th(nth,np,nr,nsa)%dBdthp(nstp)
              end do
            end if
    
            if ( nth <= nthmax .and. nr <= nrmax ) then
              read(52,iostat=ierr)nstpmax
              orbit_p(nth,np,nr,nsa)%nstp_max = nstpmax

              allocate(orbit_p(nth,np,nr,nsa)%time(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%psip(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%Babs(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%thetap(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%costh(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%sinth(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%F(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%r(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%dpsipdr(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%dFdr(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%dBdr(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%dBdthp(nstpmax))


              do nstp = 1, nstpmax
                read(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%time(nstp)
                read(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%psip(nstp)
                read(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%thetap(nstp)
                read(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%Babs(nstp)
                read(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%costh(nstp)
                read(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%sinth(nstp)
                read(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%F(nstp)
                read(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%r(nstp)
                read(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%dpsipdr(nstp)
                read(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%dFdr(nstp)
                read(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%dBdr(nstp)
                read(62,iostat=ierr)orbit_p(nth,np,nr,nsa)%dBdthp(nstp)
              end do
            end if

            if ( nth <= nthmax .and. np <= npmax ) then
              read(53,iostat=ierr)nstpmax
              orbit_r(nth,np,nr,nsa)%nstp_max = nstpmax

              allocate(orbit_r(nth,np,nr,nsa)%time(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%psip(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%Babs(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%thetap(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%costh(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%sinth(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%F(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%r(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%dpsipdr(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%dFdr(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%dBdr(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%dBdthp(nstpmax))


              do nstp = 1, nstpmax
                read(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%time(nstp)
                read(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%psip(nstp)
                read(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%thetap(nstp)
                read(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%Babs(nstp)
                read(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%costh(nstp)
                read(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%sinth(nstp)
                read(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%F(nstp)
                read(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%r(nstp)
                read(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%dpsipdr(nstp)
                read(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%dFdr(nstp)
                read(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%dBdr(nstp)
                read(63,iostat=ierr)orbit_r(nth,np,nr,nsa)%dBdthp(nstp)
              end do
            end if

          end do
        end do
      end do
    end do

    do i = 0, 3
      close(50+i)
      close(60+i)
    end do

  end subroutine load_orbit

  subroutine quantities_at_Bminimum(ra_Bmin, theta_Bmin, ob)
    use fpcomm
    use fowcomm
    real(rkind),intent(out) :: ra_Bmin, theta_Bmin
    type(orbit),intent(in) :: ob
    integer :: nstp, nstpmax, nstpmin, mm(1)
    real(rkind) :: B_min, thetap_option(7), thetap_Bmin, psip_Bmin
    real(rkind) :: nstp_near, diff_thetap(7)
    real(rkind),allocatable :: dradpsi(:)

    if ( .not.allocated(UR) ) then
      allocate(dradpsi(nrmax))
      allocate(UR(4,nrmax))
      call first_order_derivative(dradpsi, rm, psim)
      call SPL1D(psimg,rm,dradpsi,UR,nrmax,3,ierr)
    end if
    
    nstpmax = ob%nstp_max
    B_min = 1.0d10
    do nstp = 1, nstpmax
      if ( ob%Babs(nstp) < B_min ) then
        B_min = ob%Babs(nstp)
        nstpmin = nstp
      end if
    end do

    ! assume vertiucally symmmetry system, i.e. MIN(B) is at the equator. 
    thetap_option = [-3.d0*pi, -2.d0*pi, -1.d0*pi, 0.d0, pi, 2.d0*pi, 3.d0*pi]

    diff_thetap(1) = ABS( thetap_option(1) - ob%thetap(nstpmin) )
    diff_thetap(2) = ABS( thetap_option(2) - ob%thetap(nstpmin) )
    diff_thetap(3) = ABS( thetap_option(3) - ob%thetap(nstpmin) )
    diff_thetap(4) = ABS( thetap_option(4) - ob%thetap(nstpmin) )
    diff_thetap(5) = ABS( thetap_option(5) - ob%thetap(nstpmin) )
    diff_thetap(6) = ABS( thetap_option(6) - ob%thetap(nstpmin) )
    diff_thetap(7) = ABS( thetap_option(7) - ob%thetap(nstpmin) )

    mm(1) = MINLOC(diff_thetap, 1)
    thetap_Bmin = thetap_option(mm(1))

    if ( diff_thetap(mm(1)) < twopi/dble(nthpmax)*1.0d-2 ) then
      psip_Bmin  = ob%psip(nstpmin)
      theta_Bmin = ACOS( ob%costh(nstpmin) )

    else if ( nstpmin == 1 .or. nstpmin == nstpmax ) then
      psip_Bmin  = ob%psip(nstpmin)
      theta_Bmin = ACOS( ob%costh(nstpmin) )

    else ! linear interpolate
      if ( ABS( theta_Bmin - ob%thetap(nstpmin-1) ) &
          < ABS( theta_Bmin - ob%thetap(nstpmin+1) ) ) then
        nstp_near = nstpmin-1
      else
        nstp_near = nstpmin+1
      end if

      A = ( ob%psip(nstpmin) - ob%psip(nstp_near) ) &
        /( ob%thetap(nstpmin) - ob%thetap(nstp_near) )

      psip_Bmin = ob%psip(nstpmin) + A * ( thetap_Bmin - ob%thetap(nstpmin) )


      A = ( ACOS( ob%costh(nstpmin) ) - ACOS( ob%costh(nstp_near) ) ) &
        /( ob%thetap(nstpmin) - ob%thetap(nstp_near) )

      theta_Bmin = ACOS( ob%costh(nstpmin) ) + A * ( thetap_Bmin - ob%thetap(nstpmin) )

    end if

    call SPL1DF(psip_Bmin,ra_Bmin,psim,UR,nrmax,IERR)

  end subroutine quantities_at_Bminimum

  subroutine mean_ra_quantities(ob, r0, psip0, costh0, sinth0, B0, F0, dBdr0, dFdr0, dpsipdr0)
    use fpcomm
    use fowcomm
    implicit none
    type(orbit),intent(in) :: ob
    real(rkind),intent(out) :: r0, psip0, costh0, sinth0, B0, F0, dBdr0, dFdr0, dpsipdr0
    integer :: nstp, nstpmax, nstp0
    real(rkind) :: int, dt, sumr

    nstpmax = ob%nstp_max

    sumr = 0.d0
    do nstp = 2, nstpmax
      dt  = ob%time(nstp)-ob%time(nstp-1)
      int = ( ob%r(nstp)+ob%r(nstp-1) ) *0.5d0
      sumr = sumr + int*dt
    end do
    r0 = sumr/ob%time(nstpmax)

    do nstp = 2, nstpmax
      if ( ob%r(nstp-1) <= r0 .and. r0 <= ob%r(nstp) ) then
        nstp0 = nstp
        exit
      end if
    end do

    psip0  = ob%psip(nstp0) + &
          (ob%psip(nstp0)-ob%psip(nstp0-1))/(ob%r(nstp0)-ob%r(nstp0-1))*(ob%r(nstp0)-r0)

    costh0 = ob%costh(nstp0) + &
          (ob%costh(nstp0)-ob%costh(nstp0-1))/(ob%r(nstp0)-ob%r(nstp0-1))*(ob%r(nstp0)-r0)

    sinth0 = ob%sinth(nstp0) + &
          (ob%sinth(nstp0)-ob%sinth(nstp0-1))/(ob%r(nstp0)-ob%r(nstp0-1))*(ob%r(nstp0)-r0)

    B0     = ob%Babs(nstp0) + &
          (ob%Babs(nstp0)-ob%Babs(nstp0-1))/(ob%r(nstp0)-ob%r(nstp0-1))*(ob%r(nstp0)-r0)

    F0     = ob%F(nstp0) + &
          (ob%F(nstp0)-ob%F(nstp0-1))/(ob%r(nstp0)-ob%r(nstp0-1))*(ob%r(nstp0)-r0)

    dBdr0  = ob%dBdr(nstp0) + &
          (ob%dBdr(nstp0)-ob%dBdr(nstp0-1))/(ob%r(nstp0)-ob%r(nstp0-1))*(ob%r(nstp0)-r0)

    dFdr0  = ob%dFdr(nstp0) + &
          (ob%dFdr(nstp0)-ob%dFdr(nstp0-1))/(ob%r(nstp0)-ob%r(nstp0-1))*(ob%r(nstp0)-r0)

    dpsipdr0= ob%dpsipdr(nstp0) + &
          (ob%dpsipdr(nstp0)-ob%dpsipdr(nstp0-1))/(ob%r(nstp0)-ob%r(nstp0-1))*(ob%r(nstp0)-r0)

  end subroutine mean_ra_quantities
  
end module foworbit
