module foworbit
  
  use obcomm
  use obprep
  private
  public :: fow_orbit, fow_cal_local_COMs, fow_set_obparm
  public :: quantities_at_Bminimum, get_F_nstp, get_mean_ra, get_mean_psip

  double precision,allocatable :: UF(:,:), UR(:,:)

contains

  subroutine fow_set_obparm(ierr)

    use obparm
    use obinit

    implicit none
    integer,intent(out) :: ierr

    ierr = 0

    call ob_init
    call ob_parm(1,'../fp.ota/fpparm',ierr)
    nobt_max = 1
    call ob_prep(ierr)
    call ob_allocate

  end subroutine fow_set_obparm

  subroutine fow_orbit(ierr)
    use fpcomm
    use fowcomm
    use obparm
    use obinit

    implicit none
    integer,intent(out) :: ierr
    integer :: nth,np,nr,nsa,n,nmax
    real(rkind) :: end_time, begin_time, thetap_in

    ierr = 0

    if ( model_obload >= 1 ) then
      call cpu_time(begin_time)
      call load_orbit(ierr)
      call cpu_time(end_time)
      write(*,*)"Load orbit_x time :",end_time-begin_time,"[sec]"  

      if ( ierr == 0 ) then
        return
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
      write(*,*)"TASK/OB time:",end_time-begin_time,"[sec]"  
    end if

  end subroutine fow_orbit

  subroutine fow_cal_local_COMs(thetaml, psiml, tau_loss, momentum, pitch_angle, theta_pol, psi_pol, ns)
    use fpcomm
    use fowcomm
    implicit none
    real(rkind),intent(out) :: thetaml, psiml, tau_loss
    real(rkind),intent(in) :: momentum, pitch_angle, theta_pol, psi_pol
    type(orbit) :: orbit_out
    integer,intent(in) :: ns
    integer :: nstp, nstpmax, loc_max,i, im, ierr, ir, il
    real(rkind) :: psip_max, diff_thetap_min,  diff_thetap, thetap_eq(6), A

    call construct_orbit(orbit_out, momentum, pitch_angle, theta_pol, psi_pol, ns, ierr)

    nstpmax = orbit_out%nstp_max 

    if ( orbit_out%psip(nstpmax) > psi0 ) then ! loss to wall
      psiml = orbit_out%psip(nstpmax)
      thetaml = orbit_out%psip(nstpmax)
      tau_loss = orbit_out%time(nstpmax)
      return
    end if
    
    psip_max = -1.d0
    do nstp = 1, nstpmax
      if ( psip_max <= orbit_out%psip(nstp) ) then
        psip_max = orbit_out%psip(nstp)
        loc_max  = nstp
      end if
    end do

    thetap_eq = [-2.d0*pi, -1.d0*pi, 0.d0, pi, 2.d0*pi, 3.d0*pi] ! opion of thetap(psim)

    ! thetap_eq(im) is thetap(psim)
    diff_thetap_min = 1.d8
    do i = 1, 6
      diff_thetap = ABS( thetap_eq(i)-orbit_out%thetap(loc_max) )
      if ( diff_thetap < diff_thetap_min ) then
        diff_thetap_min = diff_thetap
        im = i
      end if
    end do

    ! if ( diff_thetap_min <= pi/1.d2 ) then
      psiml  = orbit_out%psip(loc_max)
      thetaml= orbit_out%theta(loc_max)
    
    ! else
    !   if ( loc_max == 1 ) then
    !     il = 1
    !     ir = 2
    !   else if ( loc_max == nstpmax ) then
    !     il = nstpmax-1
    !     ir = nstpmax
    !   else if ( ABS( orbit_out%thetap(loc_max) - orbit_out%thetap(loc_max-1) ) >= pi ) then
    !     il = loc_max
    !     ir = loc_max+1
    !   else
    !     il = loc_max-1
    !     ir = loc_max
    !   end if
  
    !   if ( orbit_out%thetap(ir) == orbit_out%thetap(il) ) then
    !     psiml = orbit_out%psip(loc_max)
    !     thetaml = orbit_out%theta(loc_max)
    !   else
    !     A = (orbit_out%psip(ir)-orbit_out%psip(il))/(orbit_out%thetap(ir)-orbit_out%thetap(il))
    !     psiml = A*(thetap_eq(im)-orbit_out%thetap(il))+orbit_out%psip(il)
    
    !     A = (orbit_out%theta(ir)-orbit_out%theta(il))/(orbit_out%thetap(ir)-orbit_out%thetap(il))
    !     psiml = A*(thetap_eq(im)-orbit_out%thetap(il))+orbit_out%theta(il)  
    !   end if
    
    ! end if

    tau_loss = 0.d0

  end subroutine fow_cal_local_COMs

  subroutine construct_orbit(orbit_out, momentum, pitch_angle, theta_pol, psi_pol, ns, ierr)
    use fpcomm
    use fowcomm

    use obcalc

    implicit none
    type(orbit),intent(out) :: orbit_out
    real(rkind),intent(in) :: momentum, pitch_angle, theta_pol, psi_pol
    integer,intent(out) :: ierr
    integer,intent(in) :: ns
    integer :: nstp, nstpmax
    real(rkind) :: pv, mass0
    logical :: is_allocated
    integer :: A,B,i

    ierr = 0

    mass0 = pa(ns)*amp
    pv = SQRT( 1.d0+momentum**2/(mass0*vc)**2 )
    ns_ob = ns

    penergy_in(1) = mass0*vc**2*(pv-1.d0)/(aee*1.d3)
    pcangle_in(1) = COS( pitch_angle )
    theta_in(1)   = theta_pol
    psipn_in(1)   = psi_pol
    zeta_in(1)    = 0.d0

    call ob_calc(ierr)

    nstpmax = nstp_max_nobt(1)+1

    if ( max_stp >= nstpmax ) then

      orbit_out%nstp_max = nstpmax

      allocate(orbit_out%time(nstpmax))
      allocate(orbit_out%psip(nstpmax))
      allocate(orbit_out%Babs(nstpmax))
      allocate(orbit_out%theta(nstpmax))
      allocate(orbit_out%thetap(nstpmax))
  
      do nstp = 1, nstpmax
        orbit_out%time(nstp)  = time_ob(nstp-1,1)
        orbit_out%psip(nstp)  = psip_ob(nstp-1,1)/psipa
        orbit_out%Babs(nstp)  = Babs_ob(nstp-1,1)
        orbit_out%thetap(nstp)= theta_ob(nstp-1,1)
  
        if ( vperp_ob(nstp-1,1) == 0.d0 ) then
          if ( vpara_ob(nstp-1,1) >= 0.d0 ) then
            orbit_out%theta(nstp) = 0.d0
          else
            orbit_out%theta(nstp) = pi
          end if
        else
          orbit_out%theta(nstp) = ACOS( vpara_ob(nstp-1,1)&
                                  /SQRT( vpara_ob(nstp-1,1)**2+vperp_ob(nstp-1,1)**2 ) )
        end if
  
      end do

      return
  
    else
      A = nstpmax / (max_stp-1)
      B = MOD(nstpmax, max_stp-1)

      orbit_out%nstp_max = max_stp
      allocate(orbit_out%time(max_stp))
      allocate(orbit_out%psip(max_stp))
      allocate(orbit_out%Babs(max_stp))
      allocate(orbit_out%theta(max_stp))
      allocate(orbit_out%thetap(max_stp))

      orbit_out%time(1)  = time_ob(0,1)
      orbit_out%psip(1)  = psip_ob(0,1)/psi0
      orbit_out%Babs(1)  = Babs_ob(0,1)
      orbit_out%thetap(1)= theta_ob(0,1)

      if ( vperp_ob(0,1) == 0.d0 ) then
        if ( vpara_ob(0,1) >= 0.d0 ) then
          orbit_out%theta(1) = 0.d0
        else
          orbit_out%theta(1) = pi
        end if
      else
        orbit_out%theta(1) = ACOS( vpara_ob(0,1)&
                                /SQRT( vpara_ob(0,1)**2+vperp_ob(0,1)**2 ) )
      end if

      do nstp = 2, max_stp-B
        i = A * (nstp-1)-1

        orbit_out%time(nstp)  = time_ob(i,1)
        orbit_out%psip(nstp)  = psip_ob(i,1)/psi0
        orbit_out%Babs(nstp)  = Babs_ob(i,1)
        orbit_out%thetap(nstp)= theta_ob(i,1)
  
        if ( vperp_ob(i,1) == 0.d0 ) then
          if ( vpara_ob(i,1) >= 0.d0 ) then
            orbit_out%theta(nstp) = 0.d0
          else
            orbit_out%theta(nstp) = pi
          end if
        else
          orbit_out%theta(nstp) = ACOS( vpara_ob(i,1)&
                                  /SQRT( vpara_ob(i,1)**2+vperp_ob(i,1)**2 ) )
        end if

      end do

      do nstp = max_stp-B+1, max_stp
        i = A*(max_stp-B-1)+(A+1)*(nstp-max_stp+B)-1

        orbit_out%time(nstp)  = time_ob(i,1)
        orbit_out%psip(nstp)  = psip_ob(i,1)/psi0
        orbit_out%Babs(nstp)  = Babs_ob(i,1)
        orbit_out%thetap(nstp)= theta_ob(i,1)
  
        if ( vperp_ob(i,1) == 0.d0 ) then
          if ( vpara_ob(i,1) >= 0.d0 ) then
            orbit_out%theta(nstp) = 0.d0
          else
            orbit_out%theta(nstp) = pi
          end if
        else
          orbit_out%theta(nstp) = ACOS( vpara_ob(i,1)&
                                  /SQRT( vpara_ob(i,1)**2+vperp_ob(i,1)**2 ) )
        end if

      end do

      return

    end if
  end subroutine construct_orbit

  subroutine construct_orbit_zero(orbit_out)
    use fpcomm
    use fowcomm
    implicit none
    type(orbit),intent(out) :: orbit_out
    
    orbit_out%nstp_max = 1
    allocate(orbit_out%psip(1))
    allocate(orbit_out%theta(1))
    allocate(orbit_out%thetap(1))
    allocate(orbit_out%Babs(1))
    allocate(orbit_out%time(1))

    orbit_out%psip(1)  = 0.d0
    orbit_out%theta(1) = 0.d0
    orbit_out%thetap(1)= 0.d0
    orbit_out%Babs(1)  = 0.d0
    orbit_out%time(1)  = 0.d0

  end subroutine construct_orbit_zero

  subroutine save_orbit(ierr)
    use fpcomm
    use fowcomm
    implicit none
    integer,intent(out) :: ierr
    character(30) :: BIN_DIR, filename
    integer :: nm ,nth, np, nr, nsa, nstp, nstpmax, it, ip, ir, im, i

    BIN_DIR = "../fp.ota/bin/"

    filename = TRIM(BIN_DIR)//"eqparm.bin"
    open(10,file=filename,access='direct',recl=rkind,form='unformatted')
    write(10,rec=1,iostat=ierr)RR
    write(10,rec=2,iostat=ierr)RA
    write(10,rec=3,iostat=ierr)RKAP
    write(10,rec=4,iostat=ierr)RDLT
    write(10,rec=5,iostat=ierr)RB
    write(10,rec=6,iostat=ierr)BB
    write(10,rec=7,iostat=ierr)RIP
    close(10)

    filename = TRIM(BIN_DIR)//"fpparm.bin"
    open(11,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=ierr)
    write(11,rec=1,iostat=ierr)nthmax
    write(11,rec=2,iostat=ierr)npmax
    write(11,rec=3,iostat=ierr)nrmax
    write(11,rec=4,iostat=ierr)nsamax
    write(11,rec=5,iostat=ierr)nthpmax
    close(11)

    filename = TRIM(BIN_DIR)//"ob_nstpmax_m.bin"
    open(50,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=ierr)
    filename = TRIM(BIN_DIR)//"ob_nstpmax_t.bin"
    open(51,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=ierr)
    filename = TRIM(BIN_DIR)//"ob_nstpmax_p.bin"
    open(52,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=ierr)
    filename = TRIM(BIN_DIR)//"ob_nstpmax_r.bin"
    open(53,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_m.bin"
    open(60,file=filename,access='direct',recl=rkind,form='unformatted')
    filename = TRIM(BIN_DIR)//"ob_t.bin"
    open(61,file=filename,access='direct',recl=rkind,form='unformatted')
    filename = TRIM(BIN_DIR)//"ob_p.bin"
    open(62,file=filename,access='direct',recl=rkind,form='unformatted')
    filename = TRIM(BIN_DIR)//"ob_r.bin"
    open(63,file=filename,access='direct',recl=rkind,form='unformatted')

    it = 1
    ip = 1
    ir = 1
    im = 1
    do nsa = 1, nsamax
      do nr = 1, nrmax+1
        do np = 1, npmax+1
          do nth = 1, nthmax+1

            if ( nth <= nthmax .and. np <= npmax .and. nr <= nrmax ) then
              nm = (nsa-1)*nthmax*npmax*nrmax&
                  +(nr-1)*nthmax*npmax&
                  +(np-1)*nthmax + nth
              write(50,rec=nm,iostat=ierr)orbit_m(nth,np,nr,nsa)%nstp_max
              nstpmax = orbit_m(nth,np,nr,nsa)%nstp_max

              do nstp = 1, nstpmax
                write(60,rec=im,iostat=ierr)orbit_m(nth,np,nr,nsa)%time(nstp)
                write(60,rec=im+1,iostat=ierr)orbit_m(nth,np,nr,nsa)%psip(nstp)
                write(60,rec=im+2,iostat=ierr)orbit_m(nth,np,nr,nsa)%theta(nstp)
                write(60,rec=im+3,iostat=ierr)orbit_m(nth,np,nr,nsa)%thetap(nstp)
                write(60,rec=im+4,iostat=ierr)orbit_m(nth,np,nr,nsa)%Babs(nstp)
                im = im+5
              end do
              
            end if

            if ( np <= npmax .and. nr <= nrmax ) then
              nm = (nsa-1)*(nthmax+1)*npmax*nrmax&
                  +(nr-1)*(nthmax+1)*npmax&
                  +(np-1)*(nthmax+1) + nth

              write(51,rec=nm,iostat=ierr)orbit_th(nth,np,nr,nsa)%nstp_max
              nstpmax = orbit_th(nth,np,nr,nsa)%nstp_max

              do nstp = 1, nstpmax
                write(61,rec=it,iostat=ierr)orbit_th(nth,np,nr,nsa)%time(nstp)
                write(61,rec=it+1,iostat=ierr)orbit_th(nth,np,nr,nsa)%psip(nstp)
                write(61,rec=it+2,iostat=ierr)orbit_th(nth,np,nr,nsa)%theta(nstp)
                write(61,rec=it+3,iostat=ierr)orbit_th(nth,np,nr,nsa)%thetap(nstp)
                write(61,rec=it+4,iostat=ierr)orbit_th(nth,np,nr,nsa)%Babs(nstp)
                it = it+5
              end do
            end if
    
            if ( nth <= nthmax .and. nr <= nrmax ) then
              nm = (nsa-1)*nthmax*(npmax+1)*nrmax&
                  +(nr-1)*nthmax*(npmax+1)&
                  +(np-1)*nthmax + nth

              write(52,rec=nm,iostat=ierr)orbit_p(nth,np,nr,nsa)%nstp_max
              nstpmax = orbit_p(nth,np,nr,nsa)%nstp_max

              do nstp = 1, nstpmax
                write(62,rec=ip,iostat=ierr)orbit_p(nth,np,nr,nsa)%time(nstp)
                write(62,rec=ip+1,iostat=ierr)orbit_p(nth,np,nr,nsa)%psip(nstp)
                write(62,rec=ip+2,iostat=ierr)orbit_p(nth,np,nr,nsa)%theta(nstp)
                write(62,rec=ip+3,iostat=ierr)orbit_p(nth,np,nr,nsa)%thetap(nstp)
                write(62,rec=ip+4,iostat=ierr)orbit_p(nth,np,nr,nsa)%Babs(nstp)
                ip = ip+5
              end do
            end if

            if ( nth <= nthmax .and. np <= npmax ) then
              nm = (nsa-1)*nthmax*npmax*(nrmax+1)&
                  +(nr-1)*nthmax*npmax&
                  +(np-1)*nthmax + nth

              write(53,rec=nm)orbit_r(nth,np,nr,nsa)%nstp_max
              nstpmax = orbit_r(nth,np,nr,nsa)%nstp_max

              do nstp = 1, nstpmax
                write(63,rec=ir,iostat=ierr)orbit_r(nth,np,nr,nsa)%time(nstp)
                write(63,rec=ir+1,iostat=ierr)orbit_r(nth,np,nr,nsa)%psip(nstp)
                write(63,rec=ir+2,iostat=ierr)orbit_r(nth,np,nr,nsa)%theta(nstp)
                write(63,rec=ir+3,iostat=ierr)orbit_r(nth,np,nr,nsa)%thetap(nstp)
                write(63,rec=ir+4,iostat=ierr)orbit_r(nth,np,nr,nsa)%Babs(nstp)
                ir = ir+5
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
    integer :: nthm_,npm_,nrm_,nsam_,access

    ierr = 0

    BIN_DIR = "../fp.ota/bin/"

    if ( access( TRIM(BIN_DIR)//"fpparm.bin", " ") /= 0 .OR. &
         access( TRIM(BIN_DIR)//"eqparm.bin", " ") /= 0) then
      ierr = 3
      return
    end if

    filename = TRIM(BIN_DIR)//"eqparm.bin"
    open(10,file=filename,access='direct',recl=rkind,form='unformatted', &
         status='old',iostat=ierr)
    read(10,rec=1,iostat=ierr)RR_
    read(10,rec=2,iostat=ierr)RA_
    read(10,rec=3,iostat=ierr)RKAP_
    read(10,rec=4,iostat=ierr)RDLT_
    read(10,rec=5,iostat=ierr)RB_
    read(10,rec=6,iostat=ierr)BB_
    read(10,rec=7,iostat=ierr)RIP_
    close(10)

    filename = TRIM(BIN_DIR)//"fpparm.bin"
    open(11,file=filename,access='direct',recl=4,form='unformatted', &
         status='old',iostat=ierr)
    read(11,rec=1,iostat=ierr)nthm_ 
    read(11,rec=2,iostat=ierr)npm_
    read(11,rec=3,iostat=ierr)nrm_ 
    read(11,rec=4,iostat=ierr)nsam_
    close(11)

    if ( RR /= RR_ .or. RA /= RA_ .or. RKAP /= RKAP_ &
        .or. RDLT /= RDLT_ .or. RB /= RB_ .or. BB /= BB_ .or. RIP /= RIP_) then
      ierr = 1
    end if
    if ( nthmax /= nthm_ .or. npmax /= npm_ &
        .or. nrmax /= nrm_ .or. nsamax /= nsam_ ) then
      ierr = 2
    end if

    if ( ierr /= 0 ) return

    filename = TRIM(BIN_DIR)//"ob_nstpmax_m.bin"
    open(50,file=filename,access='direct',recl=4,form='unformatted',status='old',iostat=ierr)
    filename = TRIM(BIN_DIR)//"ob_nstpmax_t.bin"
    open(51,file=filename,access='direct',recl=4,form='unformatted',status='old',iostat=ierr)
    filename = TRIM(BIN_DIR)//"ob_nstpmax_p.bin"
    open(52,file=filename,access='direct',recl=4,form='unformatted',status='old',iostat=ierr)
    filename = TRIM(BIN_DIR)//"ob_nstpmax_r.bin"
    open(53,file=filename,access='direct',recl=4,form='unformatted',status='old',iostat=ierr)

    filename = TRIM(BIN_DIR)//"ob_m.bin"
    open(60,file=filename,access='direct',recl=rkind,form='unformatted',status='old',iostat=ierr)
    filename = TRIM(BIN_DIR)//"ob_t.bin"
    open(61,file=filename,access='direct',recl=rkind,form='unformatted',status='old',iostat=ierr)
    filename = TRIM(BIN_DIR)//"ob_p.bin"
    open(62,file=filename,access='direct',recl=rkind,form='unformatted',status='old',iostat=ierr)
    filename = TRIM(BIN_DIR)//"ob_r.bin"
    open(63,file=filename,access='direct',recl=rkind,form='unformatted',status='old',iostat=ierr)

    im = 1
    ip = 1
    it = 1
    ir = 1
    do nsa = 1, nsamax
      do nr = 1, nrmax+1
        do np = 1, npmax+1
          do nth = 1, nthmax+1
            if ( nth <= nthmax .and. np <= npmax .and. nr <= nrmax ) then
              nm = (nsa-1)*nthmax*npmax*nrmax&
                  +(nr-1)*nthmax*npmax&
                  +(np-1)*nthmax + nth
              read(50,rec=nm)nstpmax
              orbit_m(nth,np,nr,nsa)%nstp_max = nstpmax

              allocate(orbit_m(nth,np,nr,nsa)%time(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%psip(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%Babs(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%theta(nstpmax))
              allocate(orbit_m(nth,np,nr,nsa)%thetap(nstpmax))

              do nstp = 1, nstpmax
                read(60,rec=im,iostat=ierr)orbit_m(nth,np,nr,nsa)%time(nstp)
                read(60,rec=im+1,iostat=ierr)orbit_m(nth,np,nr,nsa)%psip(nstp)
                read(60,rec=im+2,iostat=ierr)orbit_m(nth,np,nr,nsa)%theta(nstp)
                read(60,rec=im+3,iostat=ierr)orbit_m(nth,np,nr,nsa)%thetap(nstp)
                read(60,rec=im+4,iostat=ierr)orbit_m(nth,np,nr,nsa)%Babs(nstp)
                im = im+5
              end do
              
            end if

            if ( np <= npmax .and. nr <= nrmax ) then
              nm = (nsa-1)*(nthmax+1)*npmax*nrmax&
                  +(nr-1)*(nthmax+1)*npmax&
                  +(np-1)*(nthmax+1) + nth
              read(51,rec=nm)nstpmax
              orbit_th(nth,np,nr,nsa)%nstp_max = nstpmax

              allocate(orbit_th(nth,np,nr,nsa)%time(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%psip(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%Babs(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%theta(nstpmax))
              allocate(orbit_th(nth,np,nr,nsa)%thetap(nstpmax))

              do nstp = 1, nstpmax
                read(61,rec=it,iostat=ierr)orbit_th(nth,np,nr,nsa)%time(nstp)
                read(61,rec=it+1,iostat=ierr)orbit_th(nth,np,nr,nsa)%psip(nstp)
                read(61,rec=it+2,iostat=ierr)orbit_th(nth,np,nr,nsa)%theta(nstp)
                read(61,rec=it+3,iostat=ierr)orbit_th(nth,np,nr,nsa)%thetap(nstp)
                read(61,rec=it+4,iostat=ierr)orbit_th(nth,np,nr,nsa)%Babs(nstp)
                it = it+5
              end do
            end if
    
            if ( nth <= nthmax .and. nr <= nrmax ) then
              nm = (nsa-1)*nthmax*(npmax+1)*nrmax&
                  +(nr-1)*nthmax*(npmax+1)&
                  +(np-1)*nthmax + nth

              read(52,rec=nm)nstpmax
              orbit_p(nth,np,nr,nsa)%nstp_max = nstpmax

              allocate(orbit_p(nth,np,nr,nsa)%time(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%psip(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%Babs(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%theta(nstpmax))
              allocate(orbit_p(nth,np,nr,nsa)%thetap(nstpmax))

              do nstp = 1, nstpmax
                read(62,rec=ip,iostat=ierr)orbit_p(nth,np,nr,nsa)%time(nstp)
                read(62,rec=ip+1,iostat=ierr)orbit_p(nth,np,nr,nsa)%psip(nstp)
                read(62,rec=ip+2,iostat=ierr)orbit_p(nth,np,nr,nsa)%theta(nstp)
                read(62,rec=ip+3,iostat=ierr)orbit_p(nth,np,nr,nsa)%thetap(nstp)
                read(62,rec=ip+4,iostat=ierr)orbit_p(nth,np,nr,nsa)%Babs(nstp)
                ip = ip+5
              end do
            end if

            if ( nth <= nthmax .and. np <= npmax ) then
              nm = (nsa-1)*nthmax*npmax*(nrmax+1)&
                  +(nr-1)*nthmax*npmax&
                  +(np-1)*nthmax + nth

              read(53,rec=nm)nstpmax
              orbit_r(nth,np,nr,nsa)%nstp_max = nstpmax

              allocate(orbit_r(nth,np,nr,nsa)%time(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%psip(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%Babs(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%theta(nstpmax))
              allocate(orbit_r(nth,np,nr,nsa)%thetap(nstpmax))

              do nstp = 1, nstpmax
                read(63,rec=ir,iostat=ierr)orbit_r(nth,np,nr,nsa)%time(nstp)
                read(63,rec=ir+1,iostat=ierr)orbit_r(nth,np,nr,nsa)%psip(nstp)
                read(63,rec=ir+2,iostat=ierr)orbit_r(nth,np,nr,nsa)%theta(nstp)
                read(63,rec=ir+3,iostat=ierr)orbit_r(nth,np,nr,nsa)%thetap(nstp)
                read(63,rec=ir+4,iostat=ierr)orbit_r(nth,np,nr,nsa)%Babs(nstp)
                ir = ir+5
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

  subroutine quantities_at_Bminimum(ra_Bmin, theta_Bmin, orbit_in)
    use fpcomm
    use fowcomm
    real(rkind),intent(out) :: ra_Bmin, theta_Bmin
    type(orbit),intent(in) :: orbit_in
    real(rkind),allocatable :: dradpsi(:)
    integer :: nstp, nstpmax, nstpmin, mm(1)
    real(rkind) :: B_min, thetap_option(7), thetap_Bmin, psip_Bmin
    real(rkind) :: nstp_near, diff_thetap(7)

    if ( .not.allocated(UR) ) then
      allocate(dradpsi(nrmax+1))
      allocate(UR(4,nrmax+1))
      call first_order_derivative(dradpsi, rg, psimg)
      call SPL1D(psimg,rg,dradpsi,UR,nrmax+1,3,ierr)
    end if
    
    nstpmax = orbit_in%nstp_max
    B_min = 1.0d10
    do nstp = 1, nstpmax
      if ( orbit_in%Babs(nstp) < B_min ) then
        B_min = orbit_in%Babs(nstp)
        nstpmin = nstp
      end if
    end do

    ! assume vertiucally symmmetry system, i.e. MIN(B) is at the equator. 
    thetap_option = [-3.d0*pi, -2.d0*pi, -1.d0*pi, 0.d0, pi, 2.d0*pi, 3.d0*pi]

    diff_thetap(1) = ABS( thetap_option(1) - orbit_in%thetap(nstpmin) )
    diff_thetap(2) = ABS( thetap_option(2) - orbit_in%thetap(nstpmin) )
    diff_thetap(3) = ABS( thetap_option(3) - orbit_in%thetap(nstpmin) )
    diff_thetap(4) = ABS( thetap_option(4) - orbit_in%thetap(nstpmin) )
    diff_thetap(5) = ABS( thetap_option(5) - orbit_in%thetap(nstpmin) )
    diff_thetap(6) = ABS( thetap_option(6) - orbit_in%thetap(nstpmin) )
    diff_thetap(7) = ABS( thetap_option(7) - orbit_in%thetap(nstpmin) )

    mm(1) = MINLOC(diff_thetap, 1)
    thetap_Bmin = thetap_option(mm(1))

    if ( diff_thetap(mm(1)) < twopi/dble(nthpmax)*1.0d-2 ) then
      psip_Bmin  = orbit_in%psip(nstpmin)
      theta_Bmin = orbit_in%theta(nstpmin)

    else if ( nstpmin == 1 .or. nstpmin == nstpmax ) then
      psip_Bmin  = orbit_in%psip(nstpmin)
      theta_Bmin = orbit_in%theta(nstpmin)

    else ! linear interpolate
      if ( ABS( theta_Bmin - orbit_in%thetap(nstpmin-1) ) &
          < ABS( theta_Bmin - orbit_in%thetap(nstpmin+1) ) ) then
        nstp_near = nstpmin-1
      else
        nstp_near = nstpmin+1
      end if

      A = ( orbit_in%psip(nstpmin) - orbit_in%psip(nstp_near) ) &
        /( orbit_in%thetap(nstpmin) - orbit_in%thetap(nstp_near) )

      psip_Bmin = orbit_in%psip(nstpmin) + A * ( thetap_Bmin - orbit_in%thetap(nstpmin) )


      A = ( orbit_in%theta(nstpmin) - orbit_in%theta(nstp_near) ) &
        /( orbit_in%thetap(nstpmin) - orbit_in%thetap(nstp_near) )

      theta_Bmin = orbit_in%theta(nstpmin) + A * ( thetap_Bmin - orbit_in%thetap(nstpmin) )

    end if

    call SPL1DF(psip_Bmin,ra_Bmin,psimg,UR,nrmax+1,IERR)

  end subroutine

  function get_F_nstp(orbit_in, nstp) result(F_out)
    use fpcomm
    use fowcomm
    implicit none
    real(rkind) :: F_out
    type(orbit),intent(in) :: orbit_in
    integer,intent(in) :: nstp
    real(rkind),allocatable :: dFdpsi(:)
    integer :: ierr = 0

    if ( .not.allocated(UF) ) then
      allocate(dFdpsi(nrmax))
      allocate(UF(4,nrmax))
      call first_order_derivative(dFdpsi, Fpsi, psim)
      call SPL1D(psimg,Fpsig,dFdpsi,UF,nrmax,3,ierr)
    end if

    call SPL1DF(orbit_in%psip(nstp),F_out,psim,UF,nrmax,IERR)
    
  end function get_F_nstp

  function get_mean_ra(orbit_in) result(ra_out)
    use fpcomm
    use fowcomm
    implicit none
    real(rkind) :: ra_out
    type(orbit),intent(in) :: orbit_in
    integer :: nstp, nstpmax, ierr=0
    real(rkind),allocatable :: dradpsi(:)
    real(rkind) :: mean_psip, dt, psipl

    if ( .not.allocated(UR) ) then
      allocate(dradpsi(nrmax))
      allocate(UR(4,nrmax))
      call first_order_derivative(dradpsi, rm, psim)
      call SPL1D(psim,rm,dradpsi,UR,nrmax,3,ierr)
    end if

    nstpmax =  orbit_in%nstp_max
    mean_psip = 0.d0
    do nstp = 2, nstpmax
      dt = orbit_in%time(nstp)-orbit_in%time(nstp-1)
      psipl = ( orbit_in%psip(nstp) + orbit_in%psip(nstp-1) )*0.5d0
      mean_psip = mean_psip + psipl*dt
    end do
    mean_psip = mean_psip/orbit_in%time(nstpmax)

    call SPL1DF(mean_psip,ra_out,psim,UR,nrmax,IERR)

  end function get_mean_ra 

  function get_mean_psip(orbit_in) result(mean_psip)
    use fpcomm
    use fowcomm
    implicit none
    real(rkind) :: mean_psip
    type(orbit),intent(in) :: orbit_in
    integer :: nstp, nstpmax
    real(rkind) :: dt, psipl

    nstpmax =  orbit_in%nstp_max
    mean_psip = 0.d0
    do nstp = 2, nstpmax
      dt = orbit_in%time(nstp)-orbit_in%time(nstp-1)
      psipl = ( orbit_in%psip(nstp) + orbit_in%psip(nstp-1) )*0.5d0
      mean_psip = mean_psip + psipl*dt
    end do
    mean_psip = mean_psip/orbit_in%time(nstpmax)

  end function get_mean_psip 

end module foworbit
