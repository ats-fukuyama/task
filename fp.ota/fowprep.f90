module fowprep
  private
  public :: fow_prep

contains

  subroutine fow_prep

    use fowcomm
    use fpcomm
    use foworbit
    use orbit_classify

    implicit none

    type(orbit) :: orbit_pnc,ob
    integer :: ns, nth, np, nr, nsa, ierr = 0, flag, nstp, nstpmax, ir, nthp
    real(rkind) :: dummy, begin_time, end_time
    real(rkind) :: time_v_co, time_v_cnt, dt, epspsi
    real(rkind) :: thetaml, rhoml, tau_loss, momentum, pitch_angle, theta_pol, psi_pol
    logical :: isCo
    real(rkind) :: summ

    ! load equiliblium variable
    call fow_eqload(ierr)

    ! call output_orbit_classify

    ! calculate thetam
    nthm3 = nthmax/2
    nthm2 = (nthmax-nthm3)/2
    nthm1 = nthmax-nthm2-nthm3

    epspsi = psim(1)*1.d-8

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          ! calculate theta_pnc and psip_pnc_point
          call bisection_method(get_pinch_point, flag, pm(np,nsa)&
                    , theta_pnc(np,nr,nsa), psip_pnc_point(np,nr,nsa), epspsi, psim(nr)-epspsi, psim(nr), nsa)
          if ( flag /= 0 ) theta_pnc(np,nr,nsa) = NO_PINCH_ORBIT

          ! calculate theta_co_stg
          !if orbit with all theta_m are in forbitten region given pm(np) and psim(nr),then theta_co_stg(np,nr,nsa) = 0.d0
          call bisection_method(get_p_stg, flag, pm(np,nsa)&
                                      , dummy, theta_co_stg(np,nr,nsa), 0.d0, pi/2.d0, psim(nr), nsa)
          if ( flag > 0 ) theta_co_stg(np,nr,nsa) = 0.d0

          ! calculate theta_cnt_stg 
          ! if orbit with all theta_m are in forbitten region given pm(np) and psim(nr),then theta_cnt_stg(np,nr,nsa) = pi
          call bisection_method(get_p_stg, flag, pm(np,nsa)&
                                      , dummy, theta_cnt_stg(np,nr,nsa), pi/2.d0, pi, psim(nr), nsa)
          if ( flag > 0 ) theta_cnt_stg(np,nr,nsa) = pi

          ! calculate nr_pnc_point
          if ( theta_pnc(np,nr,nsa) /= NO_PINCH_ORBIT ) then
            do ir = 1, nr
              if ( psimg(ir) <=  psip_pnc_point(np,nr,nsa) &
                  .and. psip_pnc_point(np,nr,nsa) <= psimg(ir+1) ) then
                nr_pnc_point(np,nr,nsa) = ir
                exit
              end if
            end do
          end if
          ! write(6,'(A,3I3,2ES12.4,I2)')"pinch",nsa,nr,np,theta_pnc(np,nr,nsa),psip_pnc_point(np,nr,nsa),nr_pnc_point(np,nr,nsa)

          ! define thetam mesh
          if ( aefp(nsa) >= 0.d0 ) then

            if ( theta_pnc(np,nr,nsa) /= NO_PINCH_ORBIT ) then
              do nth = 1, nthm1
                delthm(nth,np,nr,nsa) = theta_pnc(np,nr,nsa)/dble(nthm1)
              end do

              do nth = nthm1+1, nthm1+nthm2
                delthm(nth,np,nr,nsa) = (theta_co_stg(np,nr,nsa)-theta_pnc(np,nr,nsa))/dble(nthm2)
              end do

              do nth = nthm1+nthm2+1, nthmax
                delthm(nth,np,nr,nsa) = (pi-theta_cnt_stg(np,nr,nsa))/dble(nthm3)
              end do

            else
              do nth = 1, nthm1+nthm2
                delthm(nth,np,nr,nsa) = theta_co_stg(np,nr,nsa)/dble(nthm1+nthm2)
              end do

              do nth = nthm1+nthm2+1, nthmax
                delthm(nth,np,nr,nsa) = (pi-theta_cnt_stg(np,nr,nsa))/dble(nthm3)
              end do

            end if
  
          else 

            if ( theta_pnc(np,nr,nsa) /= NO_PINCH_ORBIT ) then
              do nth = 1, nthm3
                delthm(nth,np,nr,nsa) = theta_co_stg(np,nr,nsa)/dble(nthm3)
              end do

              do nth = nthm3+1, nthm3+nthm2
                delthm(nth,np,nr,nsa) = (theta_pnc(np,nr,nsa)-theta_cnt_stg(np,nr,nsa))/dble(nthm2)
              end do

              do nth = nthm3+nthm2+1, nthmax
                delthm(nth,np,nr,nsa) = (pi-theta_pnc(np,nr,nsa))/dble(nthm1)
              end do

            else
              do nth = 1, nthm3
                delthm(nth,np,nr,nsa) = theta_co_stg(np,nr,nsa)/dble(nthm3)
              end do

              do nth = nthm3+1, nthmax
                delthm(nth,np,nr,nsa) = (pi-theta_cnt_stg(np,nr,nsa))/dble(nthm1+nthm2)
              end do

            end if
  
          end if

          thetam(1,np,nr,nsa) = 0.5d0*delthm(1,np,nr,nsa)
          thetamg(1,np,nr,nsa) = 0.d0

          do nth = 2, nthm1+nthm2
            thetam(nth,np,nr,nsa) = thetam(nth-1,np,nr,nsa) &
                                  +0.5d0*delthm(nth-1,np,nr,nsa)+0.5d0*delthm(nth,np,nr,nsa)
            thetamg(nth,np,nr,nsa) = thetamg(nth-1,np,nr,nsa)+delthm(nth-1,np,nr,nsa)
          end do

          thetam(nthm1+nthm2+1,np,nr,nsa) = theta_cnt_stg(np,nr,nsa) + 0.5d0*delthm(nthm1+nthm2+1,np,nr,nsa)
          thetamg(nthm1+nthm2+1,np,nr,nsa) = theta_cnt_stg(np,nr,nsa)

          do nth = nthm1+nthm2+2, nthmax
            thetam(nth,np,nr,nsa) = thetam(nth-1,np,nr,nsa) &
                                  +0.5d0*delthm(nth-1,np,nr,nsa)+0.5d0*delthm(nth,np,nr,nsa)
            thetamg(nth,np,nr,nsa) = thetamg(nth-1,np,nr,nsa)+delthm(nth-1,np,nr,nsa)
          end do

          thetamg(nthmax+1,np,nr,nsa) = pi 

        end do
      end do

      if ( aefp(nsa) >= 0.d0 ) then
        nth_pnc(nsa) = nthm1+1
        nth_stg(nsa) = nthm1+nthm2+1
      else
        nth_pnc(nsa) = nthm3+nthm2+1
        nth_stg(nsa) = nthm3+1  
      end if

    end do

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax+1
          ! calculate theta_pnc_pg and psip_pnc_point_pg
          if ( np == 1 ) then
            theta_pnc_pg(np,nr,nsa) = NO_PINCH_ORBIT
          else
            call bisection_method(get_pinch_point, flag, pg(np,nsa)&
                      , theta_pnc_pg(np,nr,nsa), psip_pnc_point_pg(np,nr,nsa), epspsi, psim(nr)-epspsi, psim(nr), nsa)

            if ( flag /= 0 ) theta_pnc_pg(np,nr,nsa) = NO_PINCH_ORBIT  
          end if
          
          ! calculate theta_co_stg_pg
          !if orbit with all theta_m are in forbitten region given pm(np) and psim(nr),then theta_co_stg_pg(np,nr,nsa) = 0.d0
          call bisection_method(get_p_stg, flag, pg(np,nsa)&
                                      , dummy, theta_co_stg_pg(np,nr,nsa), 0.d0, pi/2.d0, psim(nr), nsa)
          if ( flag > 0 ) theta_co_stg_pg(np,nr,nsa) = 0.d0

          ! calculate theta_cnt_stg _pg
          ! if orbit with all theta_m are in forbitten region given pm(np) and psim(nr),then theta_cnt_stg_pg(np,nr,nsa) = pi
          call bisection_method(get_p_stg, flag, pg(np,nsa)&
                                      , dummy, theta_cnt_stg_pg(np,nr,nsa), pi/2.d0, pi, psim(nr), nsa)
          if ( flag > 0 ) theta_cnt_stg_pg(np,nr,nsa) = pi

          ! define thetam mesh
          if ( aefp(nsa) >= 0.d0 ) then

            if ( theta_pnc_pg(np,nr,nsa) /= NO_PINCH_ORBIT ) then
              do nth = 1, nthm1
                delthm_pg(nth,np,nr,nsa) = theta_pnc_pg(np,nr,nsa)/dble(nthm1)
              end do

              do nth = nthm1+1, nthm1+nthm2
                delthm_pg(nth,np,nr,nsa) = (theta_co_stg_pg(np,nr,nsa)-theta_pnc_pg(np,nr,nsa))/dble(nthm2)
              end do

              do nth = nthm1+nthm2+1, nthmax
                delthm_pg(nth,np,nr,nsa) = (pi-theta_cnt_stg_pg(np,nr,nsa))/dble(nthm3)
              end do

            else
              do nth = 1, nthm1+nthm2
                delthm_pg(nth,np,nr,nsa) = theta_co_stg_pg(np,nr,nsa)/dble(nthm1+nthm2)
              end do

              do nth = nthm1+nthm2+1, nthmax
                delthm_pg(nth,np,nr,nsa) = (pi-theta_cnt_stg_pg(np,nr,nsa))/dble(nthm3)
              end do

            end if
  
          else 

            if ( theta_pnc_pg(np,nr,nsa) /= NO_PINCH_ORBIT ) then
              do nth = 1, nthm3
                delthm_pg(nth,np,nr,nsa) = theta_co_stg_pg(np,nr,nsa)/dble(nthm3)
              end do

              do nth = nthm3+1, nthm3+nthm2
                delthm_pg(nth,np,nr,nsa) = (theta_pnc_pg(np,nr,nsa)-theta_cnt_stg_pg(np,nr,nsa))/dble(nthm2)
              end do

              do nth = nthm3+nthm2+1, nthmax
                delthm_pg(nth,np,nr,nsa) = (pi-theta_pnc_pg(np,nr,nsa))/dble(nthm1)
              end do

            else
              do nth = 1, nthm3
                delthm_pg(nth,np,nr,nsa) = theta_co_stg_pg(np,nr,nsa)/dble(nthm3)
              end do

              do nth = nthm3+1, nthmax
                delthm_pg(nth,np,nr,nsa) = (pi-theta_cnt_stg_pg(np,nr,nsa))/dble(nthm1+nthm2)
              end do

            end if
  
          end if

          thetam_pg(1,np,nr,nsa) = 0.5d0*delthm_pg(1,np,nr,nsa)

          do nth = 2, nthm1+nthm2
            thetam_pg(nth,np,nr,nsa) = thetam_pg(nth-1,np,nr,nsa) &
                                  +0.5d0*delthm_pg(nth-1,np,nr,nsa)+0.5d0*delthm_pg(nth,np,nr,nsa)
          end do

          thetam_pg(nthm1+nthm2+1,np,nr,nsa) = theta_cnt_stg_pg(np,nr,nsa) + 0.5d0*delthm_pg(nthm1+nthm2+1,np,nr,nsa)

          do nth = nthm1+nthm2+2, nthmax
            thetam_pg(nth,np,nr,nsa) = thetam_pg(nth-1,np,nr,nsa) &
                                  +0.5d0*delthm_pg(nth-1,np,nr,nsa)+0.5d0*delthm_pg(nth,np,nr,nsa)
          end do

        end do
      end do
    end do

    do nsa = 1, nsamax
      do nr = 1, nrmax+1
        do np = 1, npmax
          ! calculate theta_pnc_rg and psip_pnc_point_rg
          if ( nr == 1 ) then
            theta_pnc_rg(np,nr,nsa) = NO_PINCH_ORBIT
          else
            call bisection_method(get_pinch_point, flag, pm(np,nsa)&
                      , theta_pnc_rg(np,nr,nsa), psip_pnc_point_rg(np,nr,nsa), epspsi, psimg(nr)-epspsi, psimg(nr), nsa)

            if ( flag /= 0 ) theta_pnc_rg(np,nr,nsa) = NO_PINCH_ORBIT  

          end if

          ! calculate theta_co_stg_rg
          !if orbit with all theta_m are in forbitten region given pm(np) and psim(nr),then theta_co_stg_rg(np,nr,nsa) = 0.d0
          call bisection_method(get_p_stg, flag, pm(np,nsa)&
                                      , dummy, theta_co_stg_rg(np,nr,nsa), 0.d0, pi/2.d0, psimg(nr), nsa)
          if ( flag > 0 ) theta_co_stg_rg(np,nr,nsa) = 0.d0

          ! calculate theta_cnt_stg_rg 
          ! if orbit with all theta_m are in forbitten region given pm(np) and psim(nr),then theta_cnt_stg_rg(np,nr,nsa) = pi
          call bisection_method(get_p_stg, flag, pm(np,nsa)&
                                      , dummy, theta_cnt_stg_rg(np,nr,nsa), pi/2.d0, pi, psimg(nr), nsa)
          if ( flag > 0 ) theta_cnt_stg_rg(np,nr,nsa) = pi

          ! define theta mesh
          if ( aefp(nsa) >= 0.d0 ) then

            if ( theta_pnc_rg(np,nr,nsa) /= NO_PINCH_ORBIT ) then
              do nth = 1, nthm1
                delthm_rg(nth,np,nr,nsa) = theta_pnc_rg(np,nr,nsa)/dble(nthm1)
              end do

              do nth = nthm1+1, nthm1+nthm2
                delthm_rg(nth,np,nr,nsa) = (theta_co_stg_rg(np,nr,nsa)-theta_pnc_rg(np,nr,nsa))/dble(nthm2)
              end do

              do nth = nthm1+nthm2+1, nthmax
                delthm_rg(nth,np,nr,nsa) = (pi-theta_cnt_stg_rg(np,nr,nsa))/dble(nthm3)
              end do

            else
              do nth = 1, nthm1+nthm2
                delthm_rg(nth,np,nr,nsa) = theta_co_stg_rg(np,nr,nsa)/dble(nthm1+nthm2)
              end do

              do nth = nthm1+nthm2+1, nthmax
                delthm_rg(nth,np,nr,nsa) = (pi-theta_cnt_stg_rg(np,nr,nsa))/dble(nthm3)
              end do

            end if
  
          else 

            if ( theta_pnc_rg(np,nr,nsa) /= NO_PINCH_ORBIT ) then
              do nth = 1, nthm3
                delthm_rg(nth,np,nr,nsa) = theta_co_stg_rg(np,nr,nsa)/dble(nthm3)
              end do

              do nth = nthm3+1, nthm3+nthm2
                delthm_rg(nth,np,nr,nsa) = (theta_pnc_rg(np,nr,nsa)-theta_cnt_stg_rg(np,nr,nsa))/dble(nthm2)
              end do

              do nth = nthm3+nthm2+1, nthmax
                delthm_rg(nth,np,nr,nsa) = (pi-theta_pnc_rg(np,nr,nsa))/dble(nthm1)
              end do

            else
              do nth = 1, nthm3
                delthm_rg(nth,np,nr,nsa) = theta_co_stg_rg(np,nr,nsa)/dble(nthm3)
              end do

              do nth = nthm3+1, nthmax
                delthm_rg(nth,np,nr,nsa) = (pi-theta_cnt_stg_rg(np,nr,nsa))/dble(nthm1+nthm2)
              end do

            end if
  
          end if

          thetam_rg(1,np,nr,nsa) = 0.5d0*delthm_rg(1,np,nr,nsa)

          do nth = 2, nthm1+nthm2
            thetam_rg(nth,np,nr,nsa) = thetam_rg(nth-1,np,nr,nsa) &
                                  +0.5d0*delthm_rg(nth-1,np,nr,nsa)+0.5d0*delthm_rg(nth,np,nr,nsa)
          end do

          thetam_rg(nthm1+nthm2+1,np,nr,nsa) = theta_cnt_stg_rg(np,nr,nsa) + 0.5d0*delthm_rg(nthm1+nthm2+1,np,nr,nsa)

          do nth = nthm1+nthm2+2, nthmax
            thetam_rg(nth,np,nr,nsa) = thetam_rg(nth-1,np,nr,nsa) &
                                  +0.5d0*delthm_rg(nth-1,np,nr,nsa)+0.5d0*delthm_rg(nth,np,nr,nsa)
          end do

        end do
      end do
    end do

    call fow_orbit(flag, ierr)

    ! IBC
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          if ( theta_pnc(np,nr,nsa) == NO_PINCH_ORBIT ) cycle

          orbit_pnc = orbit_th(nth_pnc(nsa),np,nr,nsa)
          nstpmax = orbit_pnc%nstp_max

          if ( 0.d0 <= orbit_pnc%costh(1) .and. orbit_pnc%costh(1) <= 1.d0 ) then
            isCo = .true.
          else 
            isCo = .false.
          end if

          time_v_co = 0.d0
          time_v_cnt = 0.d0
          do nstp = 2, nstpmax
            dt = orbit_pnc%time(nstp)-orbit_pnc%time(nstp-1)
            if ( isCo ) then
              if ( 0.d0 <= orbit_pnc%costh(nstp) .and. orbit_pnc%costh(nstp) <= 1.d0 ) then
                time_v_co = time_v_co + dt
              else
                time_v_co = time_v_co + dt/2.d0
                time_v_cnt = time_v_cnt + dt/2.d0
                isCo = .false.
              end if
            else
              if ( 0.d0 <= orbit_pnc%costh(nstp) .and. orbit_pnc%costh(nstp) <= 1.d0 ) then
                time_v_co = time_v_co + dt/2.d0
                time_v_cnt = time_v_cnt + dt/2.d0
                isCo = .true.
              else
                time_v_cnt = time_v_cnt + dt
              end if
            end if
          end do
          IBCflux_ratio(np,nr,nsa) = time_v_cnt/(time_v_co+time_v_cnt)
          if ( IBCflux_ratio(np,nr,nsa) == 1.d0 .or. IBCflux_ratio(np,nr,nsa) == 0.d0 ) then
            IBCflux_ratio(np,nr,nsa) = 0.5d0
          end if
          ! write(6,'(3I3,ES12.4)')np,nr,nsa,IBCflux_ratio(np,nr,nsa)
        end do
      end do
    end do

    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        do np = 1, npmax
          call search_pinch_orbit(nr_rhom_pinch(np,nr,nsa),np,nr,ns)
        end do
      end do
    end do

    call calculate_jacobian

    ! calculate local COMs
    if ( SUM( modelc ) >= 2*nsamax ) then
      if ( model_obload >= 1 .and. flag == 0) then
        call cpu_time(begin_time)
        call load_local_COM(ierr)
        call cpu_time(end_time)
        write(6,'("Load local COMs time :",ES10.3,"[sec]")')end_time-begin_time
      end if
  
      if ( model_obload == 0 .or. ierr /= 0 .or. flag /= 0) then
        write(6,'(A)')"TASK/OB is runnnig for local COMs ..."
        call cpu_time(begin_time)
        do nsa = 1, nsamax
          do nthp = 1, nthpmax
            do nr = 1, nrmax
              do np = 1, npmax
                do nth = 1, nthmax
                  momentum = pm(np,nsa)*ptfp0(nsa)
                  pitch_angle = thm(nth)
                  theta_pol = theta_p(nthp)
                  psi_pol = psim(nr)
    
                  call fow_cal_local_COMs(thetaml, rhoml, tau_loss, momentum, pitch_angle, theta_pol, psi_pol, nsa)
    
                  thetam_local(nth,np,nr,nthp,nsa) = thetaml
                  rhom_local(nth,np,nr,nthp,nsa) = rhoml
                  time_loss(nth,np,nr,nthp,nsa) = tau_loss
    
                end do
              end do
            end do
          end do
        end do
        call cpu_time(end_time)
        write(6,'("TASK/OB time:",ES10.3,"[sec]")')end_time-begin_time
        if ( model_obload > 0 ) call save_local_COM(ierr)
      end if  
    end if
    
  end subroutine fow_prep

  subroutine fow_eqload(ierr)
    use fowcomm
    use fpcomm,only:rkind,nrmax,nthmax,npmax,rm,rg
    USE libspl2d
    implicit none
    integer,intent(out):: ierr
    character(len = 80) :: line
    real(rkind) :: rr_axis,zz_axis,psit0,qaxis,qsurf
    real(rkind),allocatable,dimension(:) :: ppsi,qpsi,vpsi,rlen,ritpsi,rhotg,rhot
    real(rkind),allocatable,dimension(:,:) :: Br,Bz,Bp,Bt
    real(rkind),allocatable,dimension(:) :: temp,thpa
    integer :: nr,nthp

    real(rkind),allocatable,dimension(:,:) :: Cps, CF, Cq
    real(rkind),allocatable,dimension(:)   :: dpsim, dFpsi, dqdpsi
    real(rkind),allocatable :: dBdrho(:,:), dBdth(:,:), dBdrt(:,:), CB(:,:,:,:)

    allocate(ppsi(nrmax+1),qpsi(nrmax+1),vpsi(nrmax+1),rlen(nrmax+1),ritpsi(nrmax+1)&
            ,rhotg(nrmax+1),rhot(nrmax))
    allocate(Br(nthpmax,nrmax+1),Bz(nthpmax,nrmax+1),Bp(nthpmax,nrmax+1),Bt(nthpmax,nrmax+1))
    allocate(temp(nrmax+1),thpa(nthpmax))

    allocate(Cps(4,nrmax), CF(4,nrmax), Cq(4,nrmax), CB(4,4,nrmax,nthpmax))
    allocate(dpsim(nrmax), dFpsi(nrmax), dqdpsi(nrmax))
    allocate(dBdrho(nrmax,nthpmax), dBdth(nrmax,nthpmax), dBdrt(nrmax,nthpmax))

    ierr = 0

    modelg=3
    call eqload(modelg,knameq,ierr)
    if(ierr.ne.0) return

    write(line,'(a)') 'mdleqc = 1'   ! set boozer poloidal angle
    call eqparm(2,line,ierr)
    if(ierr.ne.0) return

    write(line,'(a,i5)') 'nrmax = ',nrmax+1
    call eqparm(2,line,ierr)
    if(ierr.ne.0) return

    write(line,'(a,i5)') 'nthmax = ',nthpmax
    call eqparm(2,line,ierr)
    if(ierr.ne.0) return

    write(line,'(a,i5)') 'nsumax = ',nrmax+1
    call eqparm(2,line,ierr)
    if(ierr.ne.0) return

    call eqcalq(ierr)
    if(ierr.ne.0) return

    call eqgetp(rhotg,psimg,nrmax+1)                       ! normalized psit radius , use only psimg
    call eqgetqn(ppsi,qpsi,Fpsig,vpsi,rlen,ritpsi,nrmax+1) ! flux functions         , use only Fpsig
    call eqgetbb(Br,Bz,Bp,Bt,nthpmax,nthpmax,nrmax+1)      ! mag field              , use only Bt and Bp
    call eqgeta(rr_axis,zz_axis,psi0,psit0,qaxis,qsurf)    ! axis and mag parameters, use only psi0

    do nthp = 1, nthpmax
      theta_p(nthp) = (dble(nthp)-0.5d0)*2.d0*pi/nthpmax
    end do
    
    psi0 = psi0/(2*pi)
    do nr = 1, nrmax+1
      Fpsig(nr) = Fpsig(nr)/(2*pi)
      psimg(nr) = psimg(nr)/(2*pi)
      if ( nr <= nrmax ) then
        rhot(nr) = rhotg(nr)
        Fpsi(nr) = Fpsig(nr)
        psim(nr) = psimg(nr)
        safety_factor(nr) = qpsi(nr)
        do nthp = 1, nthpmax
          Babs(nr,nthp) = sqrt(Bt(nthp,nr)**2+Bp(nthp,nr)**2)
        end do  
      end if
    end do

    call first_order_derivative(dpsim,psim,rhot)
    call first_order_derivative(dFpsi,Fpsi,rhot)
    call first_order_derivative(dqdpsi,safety_factor,rhot)

    call spl1D(rhot,psim,dpsim,Cps,nrmax,3,ierr)
    call spl1D(rhot,Fpsi,dFpsi,CF,nrmax,3,ierr)
    call spl1D(rhot,safety_factor,dqdpsi,Cq,nrmax,3,ierr)

    call spl2D(rhot,theta_p,Babs,dBdrho,dBdth,dBdrt,CB,nrmax,nrmax,nthpmax,0,0,ierr)
    
    ! re-define
    do nr = 1, nrmax+1
      if ( nr <= nrmax ) then
        call spl1DF(rm(nr),psim(nr),rhot,Cps,nrmax,ierr)
        call spl1DF(rm(nr),Fpsi(nr),rhot,CF,nrmax,ierr)
        call spl1DF(rm(nr),safety_factor(nr),rhot,Cq,nrmax,ierr)
        call spl2DF(rm(nr),0.d0,Bout(nr),rhot,theta_p,CB,nrmax,nrmax,nthpmax,ierr)
        call spl2DF(rm(nr),pi,Bin(nr),rhot,theta_p,CB,nrmax,nrmax,nthpmax,ierr)
      end if
      call spl1DF(rg(nr),psimg(nr),rhot,Cps,nrmax,ierr)
      call spl1DF(rg(nr),Fpsig(nr),rhot,CF,nrmax,ierr)
      call spl2DF(rg(nr),0.d0,Boutg(nr),rhot,theta_p,CB,nrmax,nrmax,nthpmax,ierr)
      call spl2DF(rg(nr),pi,Bing(nr),rhot,theta_p,CB,nrmax,nrmax,nthpmax,ierr)
    end do

    ! calculate derivatives of eq variables
    call first_order_derivative(dpsimdr,psim,rm)
    call first_order_derivative(dFdr,Fpsi,rm)
    call first_order_derivative(dBindr,Bin,rm)
    call first_order_derivative(dBoutdr,Bout,rm)

    call first_order_derivative(dpsimgdr,psimg,rg)
    call first_order_derivative(dFgdr,Fpsig,rg)
    call first_order_derivative(dBingdr,Bing,rg)
    call first_order_derivative(dBoutgdr,Boutg,rg)
  
    do nthp = 1, nthpmax
      call first_order_derivative(dBdr(:,nthp),Babs(:,nthp),rm)
    end do
    do nr = 1, nrmax
      call first_order_derivative(dBdthp(nr,:),Babs(nr,:),theta_p)
    end do

  end subroutine fow_eqload

  subroutine bisection_method(routine, convergence_flag, f, g, x, x1, x2, psim_in, nsa_in)
    ! resolve f = func(x), x1 < x < x2, where psim_in and nsa_in are constant
    !         func(x) is z1 of "routine".
    ! use bisection method if calculation has converged then return convergence_flag = 0,
    !                                                   else return convergence_flag = 1,
    !                      if func(x1)*func(x2) > 0     then return convergence_flag = 2.
    use fpcomm,only:rkind,npmax
    implicit none
    real(rkind),intent(out) :: x, g
    integer,intent(out) :: convergence_flag
    real(rkind),intent(in) :: x1, x2, f, psim_in
    integer,intent(in) :: nsa_in

    interface
      subroutine routine(z1, z2, y1, y2, m)
        double precision , intent(out) :: z1, z2
        double precision, intent(in) :: y1, y2
        integer, intent(in) :: m
      end subroutine
    end interface
    
    real(rkind) :: x_min, x_max, x_mid, f_min, f_max, f_mid, g_min, g_max, g_mid, eps
    integer :: i, imax = 1000

    convergence_flag = 0

    x_min = x1
    x_max = x2
    x_mid = (x_min+x_max)/2.d0
    eps = f*1.0d-4/npmax

    call routine(f_min, g_min, x_min, psim_in, nsa_in)
    f_min = f_min - f
    call routine(f_max, g_max, x_max, psim_in, nsa_in)
    f_max = f_max - f
    call routine(f_mid, g_mid, x_mid, psim_in, nsa_in)
    f_mid = f_mid - f

    if ( f_min*f_max > 0.d0 ) then
      convergence_flag = 2
      x = 0
      return
    end if

    do i = 1, imax

      if ( abs(f_mid) < eps ) then
        x = x_mid
        g = g_mid
        return

      else if ( f_mid*f_max < 0.d0 ) then
        x_min = x_mid
        x_mid = (x_min+x_max)/2.d0
        call routine(f_min, g_min, x_min, psim_in, nsa_in)
        f_min = f_min - f
        call routine(f_mid, g_mid, x_mid, psim_in, nsa_in)
        f_mid = f_mid - f

      else 
        x_max = x_mid
        x_mid = (x_min+x_max)/2.d0
        call routine(f_max, g_max, x_max, psim_in, nsa_in)
        f_max = f_max - f
        call routine(f_mid, g_mid, x_mid, psim_in, nsa_in)
        f_mid = f_mid - f
    
      end if

    end do

    convergence_flag = 1
    x = 0
    ! write(*,*)"do not converge bisection method",nsa_in,psim_in,f_mid
    return

  end subroutine bisection_method

  subroutine get_pinch_point(p_ret, theta_pncp, psip_in, psim_in, nsa_in)
    ! input psip_in that is pinch point of pinch orbit with psi_m = psim_in.
    ! psi_p must be less than psim(nr_in)
    ! return momentum of pinch orbit with psi_m = psim(nr_in), psi_pnc = psip_in
    use fowcomm
    use fpcomm
    
    implicit none
    real(rkind),intent(out) :: p_ret, theta_pncp
    real(rkind),intent(in):: psip_in, psim_in
    integer,intent(in) :: nsa_in

    real(rkind) :: F_pncp, Bin_pncp, BFFB, ps_ratio, dFdpsi_pncp, dBdpsi_pncp, G_m, C(3), w, FB_prime, xi_pncp, xi2
    real(rkind) :: F_m, B_m
    complex(rkind) :: z(2)
    real(rkind),allocatable :: dFdpsi(:), dBdpsi(:)
    integer :: nr

    allocate(dFdpsi(nrmax+1), dBdpsi(nrmax+1))

    call first_order_derivative(dFdpsi, Fpsig, psimg)
    call first_order_derivative(dBdpsi, Bing, psimg)

    call fow_cal_spl(F_pncp, psip_in, Fpsig, psimg)
    call fow_cal_spl(dFdpsi_pncp, psip_in, dFdpsi, psimg)
    call fow_cal_spl(F_m, psim_in, Fpsig, psimg)
    call fow_cal_spl(B_m, psim_in, Boutg, psimg)    
    
    if ( psip_in == 0.d0 ) then
      Bin_pncp = Bing(1)
      dBdpsi_pncp = dBdpsi(1)
    else
      call fow_cal_spl(Bin_pncp, psip_in, Bing, psimg)
      call fow_cal_spl(dBdpsi_pncp, psip_in, dBdpsi, psimg)  
    end if

    ! dFdpsi_pncp = dFdpsi_pncp*psip_in
    ! dBdpsi_pncp = dBdpsi_pncp*psip_in
    dFdpsi_pncp = dFdpsi_pncp*psim_in
    dBdpsi_pncp = dBdpsi_pncp*psim_in

    G_m = aefp(nsa_in)*B_m*psim_in/(amfp(nsa_in)*vc*F_m)

    ps_ratio = 1.d0-psip_in/psim_in

    BFFB = 2.d0*Bin_pncp*dFdpsi_pncp - F_pncp*dBdpsi_pncp

    C(3) = -4.d0*B_m*Bin_pncp**3*F_m**2&
          +B_m**2*(ps_ratio*BFFB+2.d0*Bin_pncp*F_pncp)**2

    C(2) = -4.d0*(Bin_pncp-B_m)*Bin_pncp**3*F_m**2&
           -2.d0*B_m**2*dBdpsi_pncp*F_pncp*ps_ratio&
           *(ps_ratio*BFFB+2.d0*Bin_pncp*F_pncp)

    C(1) = (B_m*dBdpsi_pncp*F_pncp*ps_ratio)**2

    call solve_quadratic_equation(z, C)

    w = real(z(1))

    xi2 = 1.d0-(1.d0-w)*B_m/Bin_pncp
    ! write(*,*)"pnc",xi2,w
    if ( xi2 >= 1.d0 ) xi2 = 1.d0
    if ( xi2 <= 0.d0 ) then
      if( ABS(psip_in-psim_in) >= ABS(psip_in-0.d0) )p_ret = pmax(nsa_in)
      if( ABS(psip_in-psim_in) < ABS(psip_in-0.d0) ) p_ret = 0.d0
      return
    end if

    if ( aefp(nsa_in) >= 0.d0 ) then
      xi_pncp = sqrt(xi2)
    else 
      xi_pncp = -sqrt(xi2)
    end if

    FB_prime = (dFdpsi_pncp*Bin_pncp-F_pncp*dBdpsi_pncp)/Bin_pncp**2*B_m/F_m
    p_ret = G_m*sqrt(w)/(w*FB_prime-0.5d0*(1.d0-xi_pncp**2)*F_pncp*dBdpsi_pncp/F_m/Bin_pncp) ! LHS = gamma*beta
    p_ret = vc*sqrt(p_ret**2/(1.d0+p_ret**2))             ! LHS = velocity of stagnation orbit
    if ( p_ret >= vc ) then
      p_ret = vc-1.d-5
    end if
    p_ret = amfp(nsa_in)*p_ret/(1.d0-p_ret**2/vc**2)      ! momentum of stagnation orbit
    p_ret = p_ret/ptfp0(nsa_in)                           ! normalize
    theta_pncp = acos(xi_pncp)

  end subroutine get_pinch_point

  subroutine get_p_stg(p_ret, dummy, theta_in, psim_in, nsa_in)
    ! return momentum of stagnation orbit for given (theta_m = theta_in, psi_m = psim(nr_in)) and species, nsa_in
    use fowcomm
    use fpcomm
    
    implicit none
    real(rkind),intent(out):: p_ret, dummy
    real(rkind),intent(in):: theta_in, psim_in
    integer,intent(in) :: nsa_in

    integer :: nr
    real(rkind) :: xil, B_p, F_p, G_ml
    real(rkind),allocatable :: G_m(:), B_m(:), dFdpsi(:), dBmdpsi(:)

    allocate(G_m(nrmax), B_m(nrmax),  dFdpsi(nrmax), dBmdpsi(nrmax))


    dummy = 0.d0
    xil = cos(theta_in)
    if ( xil == 0.d0 ) then
      p_ret = 0.d0
      return

    else if ( xil*aefp(nsa_in) > 0.d0 ) then
      do nr = 1, nrmax
        B_m(nr) = Bout(nr)
        G_m(nr) = aefp(nsa_in)*B_m(nr)*psim(nr)/(amfp(nsa_in)*vc*Fpsi(nr))
      end do

    else if ( xil*aefp(nsa_in) < 0.d0 ) then
      do nr = 1, nrmax
        B_m(nr) = Bin(nr)
        G_m(nr) = aefp(nsa_in)*B_m(nr)*psim(nr)/(amfp(nsa_in)*vc*Fpsi(nr))
      end do
    end if

    call first_order_derivative(dFdpsi, Fpsi, psim)
    call first_order_derivative(dBmdpsi, B_m, psim)

    do nr = 1, nrmax
      dFdpsi(nr) = dFdpsi(nr)*psim(nr)/Fpsi(nr)
      dBmdpsi(nr) = dBmdpsi(nr)*psim(nr)/B_m(nr)
    end do

    call fow_cal_spl(F_p, psim_in, dFdpsi, psim)
    call fow_cal_spl(B_p, psim_in, dBmdpsi, psim)
    call fow_cal_spl(G_ml, psim_in, G_m, psim)

    p_ret = G_ml*xil/(xil**2*F_p-0.5d0*(1.d0+xil**2)*B_p) ! LHS = gamma*beta
    p_ret = vc*sqrt(p_ret**2/(1.d0+p_ret**2))             ! LHS = velocity of stagnation orbit
    p_ret = amfp(nsa_in)*p_ret/(1.d0-p_ret**2/vc**2)      ! momentum of stagnation orbit
    p_ret = p_ret/ptfp0(nsa_in)                           ! normalize

  end subroutine get_p_stg

  subroutine save_local_COM(ierr)
    use fpcomm
    use fowcomm
    implicit none
    integer,intent(out) :: ierr
    character(30) :: BIN_DIR, filename
    integer :: nm, nth, np, nr, nsa, nthp

    ierr = 0

    BIN_DIR = "bin/"

    filename = TRIM(BIN_DIR)//"rhom_local.bin"
    open(20, file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)

    filename = TRIM(BIN_DIR)//"thetam_local.bin"
    open(21, file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)

    filename = TRIM(BIN_DIR)//"time_loss.bin"
    open(22, file=filename,access='sequential',form='unformatted',status='replace',iostat=ierr)

    do nsa = 1, nsamax
      do nthp = 1, nthpmax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              write(20,iostat=ierr)rhom_local(nth,np,nr,nthp,nsa)
              write(21,iostat=ierr)thetam_local(nth,np,nr,nthp,nsa)
              write(22,iostat=ierr)time_loss(nth,np,nr,nthp,nsa)
            end do
          end do
        end do
      end do
    end do

    close(20)
    close(21)
    close(22)

  end subroutine save_local_COM

  subroutine load_local_COM(ierr)
    use fpcomm
    use fowcomm
    implicit none
    integer,intent(out) :: ierr
    character(30) :: BIN_DIR, filename
    integer :: nth, np, nr, nsa, nthp, access
    integer :: nthm_, npm_, nrm_, nsam_, nthpm_
    real(rkind) :: RR_, RA_, RKAP_, RDLT_, RB_, BB_, RIP_

    ierr = 0

    call system('mkdir -p bin')

    BIN_DIR = "bin/"
    if ( access( TRIM(BIN_DIR)//"fpparm.bin"      , " ") /= 0 &
    .or. access( TRIM(BIN_DIR)//"eqparm.bin"      , " ") /= 0 &
    .or. access( TRIM(BIN_DIR)//"rhom_local.bin"  , " ") /= 0 &
    .or. access( TRIM(BIN_DIR)//"thetam_local.bin", " ") /= 0 &
    .or. access( TRIM(BIN_DIR)//"time_loss.bin"   , " ") /= 0 &
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
      return
    end if

    if ( nthmax  /= nthm_  &
    .or. npmax   /= npm_   &
    .or. nrmax   /= nrm_   &
    .or. nsamax  /= nsam_  &
    .or. nthpmax /= nthpm_ &
    ) then
      ierr = 1002
      return
    end if

    filename = TRIM(BIN_DIR)//"rhom_local.bin"
    open(20, file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)

    filename = TRIM(BIN_DIR)//"thetam_local.bin"
    open(21, file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)

    filename = TRIM(BIN_DIR)//"time_loss.bin"
    open(22, file=filename,access='sequential',form='unformatted',status='old',iostat=ierr)

    if ( ierr /= 0 ) return

    do nsa = 1, nsamax
      do nthp = 1, nthpmax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              read(20,iostat=ierr)rhom_local(nth,np,nr,nthp,nsa)
              read(21,iostat=ierr)thetam_local(nth,np,nr,nthp,nsa)
              read(22,iostat=ierr)time_loss(nth,np,nr,nthp,nsa)
            end do
          end do
        end do
      end do
    end do

    close(20)
    close(21)
    close(22)

  end subroutine load_local_COM

  subroutine search_pinch_orbit(nr_pinch,np_x,nr_x,ns_x)
    use fpcomm
    use fowcomm
    implicit none
    integer,intent(out) :: nr_pinch
    integer,intent(in) :: np_x,nr_x,ns_x ! index of X point
    real(rkind) :: p_pnc,p_pnc_prev, psim_pinch, epspsi, dummy, psix, p_left, p_right

    epspsi = psim(1)*1.d-8
    psim_pinch = psim(nr_x)+epspsi
    psix = psim(nr_x)
    p_left  = pg(np_x,ns_x)
    p_right = pg(np_x+1,ns_x)
    p_pnc_prev = p_left

    nr_pinch = nr_x
    do while( nr_pinch <= nrmax+1 )

      ! calcurate momentum of X point
      call get_pinch_point(p_pnc, dummy, psix, psim_pinch, ns_x)

      if ( p_left <= p_pnc .and. p_pnc <= p_right ) then
        return

      else if ( p_pnc > p_right .and. ABS( p_pnc-p_right ) < ABS( p_pnc_prev-p_right ) ) then
        return

      else if ( p_pnc > p_right .and. ABS( p_pnc-p_right ) >= ABS( p_pnc_prev-p_right ) ) then
        nr_pinch = nr_pinch - 1
        return

      else
        psim_pinch = psim_pinch + psim(1)
        nr_pinch = nr_pinch + 1
      end if

    end do

  end subroutine search_pinch_orbit

  subroutine calculate_jacobian
    use fowcomm
    use fpcomm
    use foworbit

    implicit none 
    integer :: np, nth, nr, nsa, nstpmax, ierr
    ! for JI
    real(rkind) :: J_x2z, J_z2I, tau_p, dBmdrl, pv, Bml, pl, m0, sthm, cthm
    real(rkind) :: dedp, dmudthm, dmudrm, dPzdthm, dPzdrm
    ! for JIR
    real(rkind) :: r0, psip0, costh0, sinth0, B0, F0, dBdr0, dFdr0, dpsipdr0, drmdr0
    real(rkind) :: A(2,2), b(2), detA
    real(rkind) :: sumU, sumI, normalize

    ierr = 0

    ! calculate JI 
    do nsa = 1, nsamax
      m0 = amfp(nsa)
      do nr = 1, nrmax
        do np = 1, npmax
          pl = pm(np,nsa)*ptfp0(nsa)
          do nth = 1, nthmax
            cthm = COS( thetam(nth,np,nr,nsa) )
            sthm = SIN( thetam(nth,np,nr,nsa) )
            if ( cthm*aefp(nsa) >= 0.d0 ) then
              dBmdrl = dBoutdr(nr)
              Bml = Bout(nr)
            else
              dBmdrl = dBindr(nr)
              Bml = Bin(nr)
            end if
            pv = SQRT(1.d0+theta0(nsa)*pm(np,nsa)**2)
            nstpmax = orbit_m(nth,np,nr,nsa)%nstp_max
            tau_p = orbit_m(nth,np,nr,nsa)%time(nstpmax)

            J_x2z = 4.d0*pi**2*tau_p/(m0**2*ABS(aefp(nsa)))

            dedp    = pl/(m0*pv)
            dmudthm = pl**2*cthm*sthm/(m0*Bml)
            dmudrm  = -1.d0*pl**2*sthm**2/(2.d0*m0)*dBmdrl/Bml**2
            dPzdthm = -1.d0*Fpsi(nr)/Bml*pl*sthm
            dPzdrm  = (Bml*dFdr(nr)-dBmdrl*Fpsi(nr))/Bml**2*pl*cthm-aefp(nsa)*dpsimdr(nr)

            J_z2I = ABS( dedp * (dmudthm*dPzdrm - dmudrm*dPzdthm) )

            JI(nth,np,nr,nsa) = J_x2z * J_z2I
            
          end do
        end do
      end do
    end do


    ! calculate JIR
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          pl = pm(np,nsa)*ptfp0(nsa)
          do nth = 1, nthmax 
            cthm = COS( thetam(nth,np,nr,nsa) )
            sthm = SIN( thetam(nth,np,nr,nsa) )
            if ( cthm*aefp(nsa) >= 0.d0 ) then
              Bml = Bout(nr)
              dBmdrl = dBoutdr(nr)
            else
              Bml = Bin(nr)
              dBmdrl = dBindr(nr)
            end if

            call mean_ra_quantities(orbit_m(nth,np,nr,nsa), r0, psip0, costh0, sinth0, B0, F0, dBdr0, dFdr0, dpsipdr0)

            A(1,1) = 2.d0*sthm*cthm/Bml
            A(1,2) = -1.d0*sthm**2*dBmdrl/Bml**2
            A(2,1) = Fpsi(nr)/Bml*pl*sthm
            A(2,2) = aefp(nsa)*dpsimdr(nr)-( dFdr(nr)*Bml-Fpsi(nr)*dBmdrl )/Bml**2*pl*cthm
            detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)

            if ( detA /= 0.d0 ) then
              b(1) = -1.d0*sinth0**2/B0**2*dBdr0
              b(2) = aefp(nsa)*dpsipdr0 - ( dFdr0*B0-F0*dBdr0 )/B0**2*pl*costh0
              drmdr0 = (A(1,1)*b(2)-A(2,1)*b(1))/detA
            else
              drmdr0 = 1.d0
            end if

            JIR(nth,np,nr,nsa) = JI(nth,np,nr,nsa)*ABS( drmdr0 )/r0
          end do
        end do
      end do

      sumU = 0.d0
      sumI = 0.d0
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            sumI = sumI + delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)
            sumU = sumU + volp(nth,np,nsa)
          end do
        end do
      end do

      normalize = sumU/sumI

      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            JI(nth,np,nr,nsa)  = JI(nth,np,nr,nsa)  * normalize
            JIR(nth,np,nr,nsa) = JIR(nth,np,nr,nsa) * normalize
          end do
        end do
      end do

    end do

  end subroutine calculate_jacobian

end module fowprep
