module fowprep
  private
  public :: fow_prep

  double precision,allocatable :: psim0(:), psimg0(:)

contains

  subroutine fow_prep

    use fowcomm
    use fpcomm
    use foworbit
    use fowdistribution

    use fpwrite

    implicit none

    type(orbit) :: orbit_pnc
    integer :: nth, np, nr, nsa, ierr = 0, flag, nstp, nstpmax, ir, nthp
    real(rkind) :: dummy, begin_time, end_time
    real(rkind) :: time_v_co, time_v_cnt, dt, epspsi
    real(rkind) :: thetaml, psiml, tau_loss, momentum, pitch_angle, theta_pol, psi_pol
    logical :: isCo
    real(rkind) :: summ


    do nth = 1,nthmax
      xi(nth) = (nth-0.5d0)/nthmax*(-2.d0)+1.d0
      xig(nth) = (nth-1.d0)/nthmax*(-2.d0)+1.d0
    end do
    xig(nthmax+1) = -1.d0

    ! load equiliblium variable
    call fow_eqload(ierr)

    do nr = 1, nrmax
      delps(nr) = psimg(nr+1)-psimg(nr)
    end do

    do nthp = 1, nthpmax
      theta_p(nthp) = (nthp-1)*2.d0*pi/nthpmax
    end do

    ! calculate thetam
    nthm3 = nthmax/2
    nthm2 = (nthmax-nthm3)/2
    nthm1 = nthmax-nthm2-nthm3

    epspsi = psim(1)*psi0*1.d-2

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          ! calculate theta_pnc and psip_pnc_point
          call bisection_method_for_IBC(get_pinch_point, flag, pm(np,nsa)&
                    , theta_pnc(np,nr,nsa), psip_pnc_point(np,nr,nsa), epspsi, psim(nr)*psi0-epspsi, psim(nr)*psi0, nsa)

          psip_pnc_point(np,nr,nsa) = psip_pnc_point(np,nr,nsa) / psi0
          if ( flag /= 0 ) theta_pnc(np,nr,nsa) = NO_PINCH_ORBIT

          ! calculate theta_co_stg
          !if orbit with all theta_m are in forbitten region given pm(np) and psim(nr),then theta_co_stg(np,nr,nsa) = 0.d0
          call bisection_method_for_IBC(get_p_stg, flag, pm(np,nsa)&
                                      , dummy, theta_co_stg(np,nr,nsa), 0.d0, pi/2.d0, psim(nr)*psi0, nsa)
          if ( flag > 0 ) theta_co_stg(np,nr,nsa) = 0.d0

          ! calculate theta_cnt_stg 
          ! if orbit with all theta_m are in forbitten region given pm(np) and psim(nr),then theta_cnt_stg(np,nr,nsa) = pi
          call bisection_method_for_IBC(get_p_stg, flag, pm(np,nsa)&
                                      , dummy, theta_cnt_stg(np,nr,nsa), pi/2.d0, pi, psim(nr)*psi0, nsa)
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
            call bisection_method_for_IBC(get_pinch_point, flag, pg(np,nsa)&
                      , theta_pnc_pg(np,nr,nsa), psip_pnc_point_pg(np,nr,nsa), epspsi, psim(nr)*psi0-epspsi, psim(nr)*psi0, nsa)

            psip_pnc_point_pg(np,nr,nsa) = psip_pnc_point_pg(np,nr,nsa) / psi0

            if ( flag /= 0 ) theta_pnc_pg(np,nr,nsa) = NO_PINCH_ORBIT  
          end if
          
          ! calculate theta_co_stg_pg
          !if orbit with all theta_m are in forbitten region given pm(np) and psim(nr),then theta_co_stg_pg(np,nr,nsa) = 0.d0
          call bisection_method_for_IBC(get_p_stg, flag, pg(np,nsa)&
                                      , dummy, theta_co_stg_pg(np,nr,nsa), 0.d0, pi/2.d0, psim(nr)*psi0, nsa)
          if ( flag > 0 ) theta_co_stg_pg(np,nr,nsa) = 0.d0

          ! calculate theta_cnt_stg _pg
          ! if orbit with all theta_m are in forbitten region given pm(np) and psim(nr),then theta_cnt_stg_pg(np,nr,nsa) = pi
          call bisection_method_for_IBC(get_p_stg, flag, pg(np,nsa)&
                                      , dummy, theta_cnt_stg_pg(np,nr,nsa), pi/2.d0, pi, psim(nr)*psi0, nsa)
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
            call bisection_method_for_IBC(get_pinch_point, flag, pm(np,nsa)&
                      , theta_pnc_rg(np,nr,nsa), psip_pnc_point_rg(np,nr,nsa), epspsi, psimg(nr)*psi0-epspsi, psimg(nr)*psi0, nsa)

            psip_pnc_point_rg(np,nr,nsa) = psip_pnc_point_rg(np,nr,nsa) / psi0

            if ( flag /= 0 ) theta_pnc_rg(np,nr,nsa) = NO_PINCH_ORBIT  

          end if

          ! calculate theta_co_stg_rg
          !if orbit with all theta_m are in forbitten region given pm(np) and psim(nr),then theta_co_stg_rg(np,nr,nsa) = 0.d0
          call bisection_method_for_IBC(get_p_stg, flag, pm(np,nsa)&
                                      , dummy, theta_co_stg_rg(np,nr,nsa), 0.d0, pi/2.d0, psimg(nr)*psi0, nsa)
          if ( flag > 0 ) theta_co_stg_rg(np,nr,nsa) = 0.d0

          ! calculate theta_cnt_stg_rg 
          ! if orbit with all theta_m are in forbitten region given pm(np) and psim(nr),then theta_cnt_stg_rg(np,nr,nsa) = pi
          call bisection_method_for_IBC(get_p_stg, flag, pm(np,nsa)&
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

    call fow_orbit(ierr)

    ! calculate the boundary flux partition
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          if ( theta_pnc(np,nr,nsa) == NO_PINCH_ORBIT ) cycle

          orbit_pnc = orbit_th(nth_pnc(nsa),np,nr,nsa)
          nstpmax = orbit_pnc%nstp_max

          if ( 0.d0 <= orbit_pnc%theta(1) .and. orbit_pnc%theta(1) <= pi/2.d0 ) then
            isCo = .true.
          else 
            isCo = .false.
          end if

          time_v_co = 0.d0
          time_v_cnt = 0.d0
          do nstp = 2, nstpmax
            dt = orbit_pnc%time(nstp)-orbit_pnc%time(nstp-1)
            if ( isCo ) then
              if ( 0.d0 <= orbit_pnc%theta(nstp) .and. orbit_pnc%theta(nstp) <= pi/2.d0 ) then
                time_v_co = time_v_co + dt
              else
                time_v_co = time_v_co + dt/2.d0
                time_v_cnt = time_v_cnt + dt/2.d0
                isCo = .false.
              end if
            else
              if ( 0.d0 <= orbit_pnc%theta(nstp) .and. orbit_pnc%theta(nstp) <= pi/2.d0 ) then
                time_v_co = time_v_co + dt/2.d0
                time_v_cnt = time_v_cnt + dt/2.d0
                isCo = .true.
              else
                time_v_cnt = time_v_cnt + dt
              end if
            end if
          end do
          IBCflux_ratio(np,nr,nsa) = time_v_cnt/(time_v_co+time_v_cnt)
        end do
      end do
    end do

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          call construst_Xstg_as_pncp(Xstg_as_pncp(np,nr,nsa))

          do ir = 1, nrmax
            if ( nr_pnc_point(np,ir,nsa) == nr ) then
              call add_pnc_orbit(Xstg_as_pncp(np,nr,nsa), ir)
            end if
          end do

        end do
      end do
    end do

    call calculate_jacobian

    ! calculate local COMs
    if ( model_obload >= 1 ) then
      call cpu_time(begin_time)
      call load_local_COM(ierr)
      call cpu_time(end_time)
      write(6,'("Load local COMs time :",ES10.3,"[sec]")')end_time-begin_time
    end if

    if ( model_obload == 0 .or. ierr /= 0) then
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
  
                call fow_cal_local_COMs(thetaml, psiml, tau_loss, momentum, pitch_angle, theta_pol, psi_pol, nsa)
  
                thetam_local(nth,np,nr,nthp,nsa) = thetaml
                psim_local(nth,np,nr,nthp,nsa) = psiml
                time_loss(nth,np,nr,nthp,nsa) = tau_loss
  
              end do
            end do
          end do
        end do
      end do
      call cpu_time(end_time)
      write(6,'("TASK/OB time:",ES10.3,"[sec]")')end_time-begin_time
      if ( model_obload >= 0 ) call save_local_COM(ierr)
    end if
    
    do nsa = 1, nsamax
      call fI_Maxwellian_sub(fnsp(:,:,:,nsa), nsa)
    end do

    do nsa = 1, nsamax
      fnorm(nsa) = 1.d0!1.d40*RNFP0(NSA)
    end do

  end subroutine fow_prep

  subroutine fow_eqload(ierr)
    use fowcomm
    use fpcomm,only:rkind,nrmax,nthmax,npmax
    USE libgrf
    implicit none
    integer,intent(out):: ierr
    character(len = 80) :: line
    real(rkind) :: rr_axis,zz_axis,psit0,qaxis,qsurf
    real(rkind),allocatable,dimension(:) :: ppsi,qpsi,vpsi,rlen,ritpsi,rhot
    real(rkind),allocatable,dimension(:,:) :: Br,Bz,Bp,Bt
    real(rkind),allocatable,dimension(:) :: temp,thpa
    integer :: nr,nthp

    allocate(ppsi(nrmax+1),qpsi(nrmax+1),vpsi(nrmax+1),rlen(nrmax+1),ritpsi(nrmax+1)&
            ,rhot(nrmax+1))
    allocate(Br(nthpmax,nrmax+1),Bz(nthpmax,nrmax+1),Bp(nthpmax,nrmax+1),Bt(nthpmax,nrmax+1))
    ALLOCATE(temp(nrmax+1),thpa(nthpmax))

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

    call eqgetp(rhot,psimg,nrmax+1)                        ! normalized psit radius , use only psimg
    call eqgetqn(ppsi,qpsi,Fpsig,vpsi,rlen,ritpsi,nrmax+1) ! flux functions         , use only Fpsig
    call eqgetbb(Br,Bz,Bp,Bt,nthpmax,nthpmax,nrmax+1)      ! mag field              , use only Bt and Bp
    call eqgeta(rr_axis,zz_axis,psi0,psit0,qaxis,qsurf)    ! axis and mag parameters, use only psi0

    ! CALL pages
    ! CALL grd1d(1,rhot,psimg,nrmax+1,nrmax+1,1,'@psimg@')
    ! CALL grd1d(2,rhot,Fpsig,nrmax+1,nrmax+1,1,'@Fpsig@')
    ! CALL grd1d(3,rhot,ppsi, nrmax+1,nrmax+1,1,'@ppsi@')
    ! CALL grd1d(4,rhot,qpsi, nrmax+1,nrmax+1,1,'@qpsi@')
    ! CALL pagee

    ! CALL pages
    ! DO nthp=1,nthpmax
    !    thpa(nthp)=2.D0*PI*(nthp-1)/DBLE(nthpmax)
    ! END DO
    ! CALL grd2d(1,thpa,rhot,Br,nthpmax,nthpmax,nrmax+1,'@Br@')
    ! CALL grd2d(2,thpa,rhot,Bz,nthpmax,nthpmax,nrmax+1,'@Bz@')
    ! CALL grd2d(3,thpa,rhot,Bp,nthpmax,nthpmax,nrmax+1,'@Bp@')
    ! CALL grd2d(4,thpa,rhot,Bt,nthpmax,nthpmax,nrmax+1,'@Bt@')
    ! CALL pagee

    ! CALL pages
    ! CALL grd1d(1,rhot,psimg,nrmax+1,nrmax+1,1,'@psimg@')
    ! temp(1)=0.D0
    ! DO nr=1,nrmax
    !    temp(nr+1)=(psimg(nr+1)-psimg(nr))/(rhot(nr+1)-rhot(nr))
    ! END DO
    ! CALL grd1d(2,rhot,temp,nrmax+1,nrmax+1,1,'@dpsimg@')
    ! DO nr=2,nrmax-1
    !    temp(nr)=(psimg(nr+1)-2*psimg(nr)+psimg(nr-1)) &
    !               /(rhot(nr+1)-rhot(nr))**2
    ! END DO
    ! temp(1)=2.D0*temp(2)-temp(3)
    ! temp(nrmax+1)=temp(nrmax)
    ! CALL grd1d(3,rhot,temp,nrmax+1,nrmax+1,1,'@ddpsimg@')
    ! CALL pagee
    
    psi0 = psi0/(2*pi)
    do nr = 1, nrmax+1
      Fpsig(nr) = Fpsig(nr)/(2*pi)
      psimg(nr) = psimg(nr)/(2*pi)
      do nthp = 1, nthpmax
        Babs(nr,nthp) = sqrt(Bt(nthp,nr)**2+Bp(nthp,nr)**2)
      end do
    end do

    call calculate_equator_variable(ierr)

  end subroutine

  subroutine calculate_equator_variable(ierr)
    use fowcomm
    use fpcomm,only:rkind,nrmax,nthmax,npmax

    implicit none
    integer,intent(inout) :: ierr
    real(rkind),allocatable,dimension(:) :: x,f,fx,g,gx,h1,h2,h1x,h2x
    real(rkind),allocatable,dimension(:,:) :: U,V,W1,W2
    real(rkind) :: dps0,dps
    integer :: nr,i,j

    ierr = 0

    allocate(x(nrmax+1),f(nrmax+1),fx(nrmax+1),U(4,nrmax+1))
    allocate(g(nrmax+1),gx(nrmax+1),V(4,nrmax+1))
    allocate(h1(nrmax+1),h1x(nrmax+1),W1(4,nrmax+1))
    allocate(h2(nrmax+1),h2x(nrmax+1),W2(4,nrmax+1))

    do nr = 1,nrmax+1
      x(nr) = (nr-1)*1.d0
      f(nr) = psimg(nr)
      g(nr) = Fpsig(nr)
      h1(nr) = Babs(nr,1)
      h2(nr) = Babs(nr,(nthpmax+1)/2)
    end do

    call SPL1D(x,f,fx,U,nrmax+1,0,IERR)
    call SPL1D(x,g,gx,V,nrmax+1,0,IERR)
    call SPL1D(x,h1,h1x,W1,nrmax+1,0,IERR)
    call SPL1D(x,h2,h2x,W2,nrmax+1,0,IERR)

    Boutg(1) = h1(1)
    Bing(1)  = h2(1)

    do nr = 2,nrmax+1
      if ( psi0 <= cal_spl1D(U(:,nr),1.d0) ) then
        dps0 = 1.d0
        call newton_spl1D(U(:,nr),psi0,dps0)
        dps = (x(nr-1)+dps0)/nrmax

        ! define psimg and Fpsimg
        do i = 2, nrmax+1 ! i is label of psimg and Fpsig
          do j = 2, nrmax+1 
            if ( x(j-1) < (i-1)*dps .and. (i-1)*dps <= x(j) ) then
              psimg(i) = cal_spl1D(U(:,j), (i-1)*dps-x(j-1)) / psi0
              Fpsig(i) = cal_spl1D(V(:,j), (i-1)*dps-x(j-1))
              Boutg(i) = cal_spl1D(W1(:,j), (i-1)*dps-x(j-1))
              Bing(i)  = cal_spl1D(W2(:,j), (i-1)*dps-x(j-1))
              exit
            end if
          end do
        end do

        ! define psim and Fpsim
        do i = 1, nrmax ! i is label of psim and Fpsi
          do j = 2, nrmax+1 
            if ( x(j-1) < (i-0.5d0)*dps .and. (i-0.5d0)*dps <= x(j) ) then
              psim(i) = cal_spl1D( U(:,j), (i-0.5d0)*dps-x(j-1)) / psi0
              Fpsi(i) = cal_spl1D(V(:,j), (i-0.5d0)*dps-x(j-1))
              Bout(i) = cal_spl1D(W1(:,j), (i-0.5d0)*dps-x(j-1))
              Bin(i)  = cal_spl1D(W2(:,j), (i-0.5d0)*dps-x(j-1))
              exit
            end if
          end do
        end do

        exit
      end if
    end do

  end subroutine calculate_equator_variable

  subroutine bisection_method_for_IBC(routine, convergence_flag, f, g, x, x1, x2, psim_in, nsa_in)
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
    write(*,*)"do not converge bisection method",nsa_in,psim_in
    return

  end subroutine bisection_method_for_IBC

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


    if ( .not.allocated(psim0) .or. .not.allocated(psimg0) ) then
      allocate(psim0(nrmax), psimg0(nrmax+1))
      do nr = 1, nrmax
        psim0(nr) = psim(nr)*psi0
        psimg0(nr) = psimg(nr)*psi0
      end do
      psimg0(nrmax+1) = psimg(nrmax+1)*psi0
    end if


    call first_order_derivative(dFdpsi, Fpsig, psimg0)
    call first_order_derivative(dBdpsi, Bing, psimg0)

    call fow_cal_spl(F_pncp, psip_in, Fpsig, psimg0)
    call fow_cal_spl(dFdpsi_pncp, psip_in, dFdpsi, psimg0)
    call fow_cal_spl(F_m, psim_in, Fpsig, psimg0)
    call fow_cal_spl(B_m, psim_in, Boutg, psimg0)    
    
    if ( psip_in == 0.d0 ) then
      Bin_pncp = Bing(1)
      dBdpsi_pncp = dBdpsi(1)
    else
      call fow_cal_spl(Bin_pncp, psip_in, Bing, psimg0)
      call fow_cal_spl(dBdpsi_pncp, psip_in, dBdpsi, psimg0)  
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


    if ( .not.allocated(psim0) .or. .not.allocated(psimg0) ) then
      allocate(psim0(nrmax), psimg0(nrmax+1))
      do nr = 1, nrmax
        psim0(nr) = psim(nr)*psi0
        psimg0(nr) = psimg(nr)*psi0
      end do
      psimg0(nrmax+1) = psimg(nrmax+1)*psi0
    end if


    dummy = 0.d0
    xil = cos(theta_in)
    if ( xil == 0.d0 ) then
      p_ret = 0.d0
      return

    else if ( xil*aefp(nsa_in) > 0.d0 ) then
      do nr = 1, nrmax
        B_m(nr) = Bout(nr)
        G_m(nr) = aefp(nsa_in)*B_m(nr)*psim0(nr)/(amfp(nsa_in)*vc*Fpsi(nr))
      end do

    else if ( xil*aefp(nsa_in) < 0.d0 ) then
      do nr = 1, nrmax
        B_m(nr) = Bin(nr)
        G_m(nr) = aefp(nsa_in)*B_m(nr)*psim0(nr)/(amfp(nsa_in)*vc*Fpsi(nr))
      end do
    end if

    call first_order_derivative(dFdpsi, Fpsi, psim0)
    call first_order_derivative(dBmdpsi, B_m, psim0)

    do nr = 1, nrmax
      dFdpsi(nr) = dFdpsi(nr)*psim0(nr)/Fpsi(nr)
      dBmdpsi(nr) = dBmdpsi(nr)*psim0(nr)/B_m(nr)
    end do

    call fow_cal_spl(F_p, psim_in, dFdpsi, psim0)
    call fow_cal_spl(B_p, psim_in, dBmdpsi, psim0)
    call fow_cal_spl(G_ml, psim_in, G_m, psim0)

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

    filename = TRIM(BIN_DIR)//"psim_local.bin"
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
              write(20,iostat=ierr)psim_local(nth,np,nr,nthp,nsa)
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
    .or. access( TRIM(BIN_DIR)//"psim_local.bin"  , " ") /= 0 &
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

    filename = TRIM(BIN_DIR)//"psim_local.bin"
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
              read(20,iostat=ierr)psim_local(nth,np,nr,nthp,nsa)
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

  subroutine calculate_jacobian
    ! calculate Jacabian_I(nth,np,nr,nsa)

    use fowcomm
    use fpcomm
    use foworbit

    implicit none 
    integer :: np, nth, nr, nsa, nstpmax
    ! for JI
    real(rkind) :: J_x2z, J_z2I, tau_p, dBdpsil, pv, Bml, pl, m0, sthm, cthm
    real(rkind) :: dedp, dmudthm, dmudpsm, dPzdthm, dPzdpsm
    real(rkind) :: sumJ
    real(rkind),allocatable :: dFdpsi(:), dBdpsi(:,:)
    ! for JIR
    integer :: nr0, nstp, ierr, sumi
    real(rkind) :: Vtube, dpsmdr0, psip0, dth0dr0, dFB0dr0, dBinvdr0, dps0dr0
    real(rkind) :: r0pls, r0mns, F0pls, F0mns, dr0, dpsmsr0, r0
    real(rkind) :: F0, B0, cth0, sth0, dPzdr0, dmudr0, fact1, fact2, fact, normalize
    type(orbit) :: ob
    real(rkind),allocatable :: UR(:,:), UPS(:,:)
    real(rkind),allocatable :: dradpsi(:), dpsdra(:)

    ierr = 0

    allocate(dFdpsi(nrmax), dBdpsi(nrmax,2))
    call first_order_derivative(dFdpsi,Fpsi,psim)
    call first_order_derivative(dBdpsi(:,1),Bout,psim)
    call first_order_derivative(dBdpsi(:,2),Bin,psim)

    ! calculate JI 
    do nsa = 1, nsamax
      m0 = amfp(nsa)
      sumJ = 0.d0
      do nr = 1, nrmax
        do np = 1, npmax
          pl = pm(np,nsa)*ptfp0(nsa)
          do nth = 1, nthmax
            cthm = COS( thetam(nth,np,nr,nsa) )
            sthm = SIN( thetam(nth,np,nr,nsa) )
            if ( cthm*aefp(nsa) >= 0.d0 ) then
              dBdpsil = dBdpsi(nr,1)
              Bml = Bout(nr)
            else
              dBdpsil = dBdpsi(nr,2)
              Bml = Bin(nr)
            end if
            pv = SQRT(1.D0+THETA0(NSA)*PM(NP,NSA)**2)
            nstpmax = orbit_m(nth,np,nr,nsa)%nstp_max
            tau_p = orbit_m(nth,np,nr,nsa)%time(nstpmax)

            J_x2z = 4.d0*pi**2*tau_p/(m0**2*abs(aefp(nsa)))

            dedp    = pl/(m0*pv)
            dmudthm = pl**2*cthm*sthm/(m0*Bml)
            dmudpsm = -1.d0*pl**2*sthm**2/(2.d0*m0)*dBdpsil/Bml**2
            dPzdthm = -1.d0*Fpsi(nr)/Bml*pl*sthm
            dPzdpsm = (Bml*dFdpsi(nr)-dBdpsil*Fpsi(nr))/Bml**2*pl*cthm-aefp(nsa)
            J_z2I = ABS( dedp * (dmudthm*dPzdpsm - dmudpsm*dPzdthm) )

            JI(nth,np,nr,nsa) = J_x2z * J_z2I * ptfp0(nsa)*psi0
            sumJ = sumJ + JI(nth,np,nr,nsa)*delps(nr)*delthm(nth,np,nr,nsa)*delp(nsa)
            
          end do
        end do
      end do      
    end do


    ! calculate JIR
    allocate(UR(4,nrmax))
    allocate(dradpsi(nrmax))
    call first_order_derivative(dradpsi, rm, psim)
    call SPL1D(psim,rm,dradpsi,UR,nrmax,3,ierr)

    allocate(dpsdra(nrmax))
    call first_order_derivative(dpsdra, psim, rm)

    do nsa = 1, nsamax
      sumi = 0.d0
      sumJ = 0.d0
      do nr = 1, nrmax
        do np = 1, npmax
          pl = pm(np,nsa)*ptfp0(nsa)
          do nth = 1, nthmax
            ob = orbit_m(nth,np,nr,nsa)
            r0 = get_mean_ra(ob)
            call fow_cal_spl(dpsmsr0, r0, dpsdra, rm)
            ! psip0 = get_mean_psip(ob)
            ! nstpmax = ob%nstp_max

            ! do nstp = 1, nstpmax-1
            !   if ( ob%psip(nstp) <= psip0 .and. psip0 < ob%psip(nstp+1) ) then
            !     call SPL1DF(ob%psip(nstp),r0mns,psim,UR,nrmax,IERR)
            !     call SPL1DF(ob%psip(nstp+1),r0pls,psim,UR,nrmax,IERR)
            !     dr0 = r0pls-r0mns
            !     F0mns = get_F_nstp(ob, nstp)
            !     F0pls = get_F_nstp(ob, nstp+1)
          
            !     dth0dr0  = (ob%theta(nstp+1)-ob%theta(nstp))/dr0
            !     dFB0dr0  = (F0pls/ob%Babs(nstp+1)-F0mns/ob%Babs(nstp))/dr0
            !     dBinvdr0 = (1.d0/ob%Babs(nstp+1)-1.d0/ob%Babs(nstp))/dr0
            !     dps0dr0  = (ob%psip(nstp+1)-ob%psip(nstp))/dr0

            !     F0 = 0.5d0*F0pls + 0.5d0*F0mns
            !     B0 = 0.5d0*ob%Babs(nstp+1) + 0.5d0*ob%Babs(nstp)
            !     sth0 = 0.5d0*SIN( ob%theta(nstp+1) ) + 0.5d0*SIN( ob%theta(nstp) )
            !     cth0 = 0.5d0*COS( ob%theta(nstp+1) ) + 0.5d0*COS( ob%theta(nstp) )
            !     exit
            !   end if
            ! end do

            ! cthm = COS( thetam(nth,np,nr,nsa) )
            ! sthm = SIN( thetam(nth,np,nr,nsa) )
            ! if ( cthm*aefp(nsa) >= 0.d0 ) then
            !   dBdpsil = dBdpsi(nr,1)
            !   Bml = Bout(nr)
            ! else
            !   dBdpsil = dBdpsi(nr,2)
            !   Bml = Bin(nr)
            ! end if
            ! dPzdr0 = dFB0dr0*pl*cth0-F0/B0*pl*sth0*dth0dr0-aefp(nsa)*dps0dr0
            ! dmudr0 = dBinvdr0*sth0**2+dth0dr0*sth0*cth0/B0

            ! fact1 = 2.d0*cthm*( (Bml*dFdpsi(nr)-dBdpsil*Fpsi(nr))/Bml**2*pl*cthm-aefp(nsa) )
            ! fact2 = -1.d0*pl*Fpsi(nr)*dBdpsil/Bml**2*sthm**2

            ! fact = fact1+fact2
            ! dpsmsr0 = ( 2.d0*dPzdr0*cthm + dmudr0*pl*Fpsi(nr) )/fact
            ! write(*,*)dpsmsr0
            Vtube = 2.d0*pi*RR*2.d0*pi*RA*RA*r0
            JIR(nth,np,nr,nsa) = JI(nth,np,nr,nsa)*ABS( dpsmsr0 )*2.d0*pi/Vtube
            sumJ = sumJ+JIR(nth,np,nr,nsa)*delp(nsa)*delthm(nth,np,nr,nsa)
            sumi = sumi + 1
          end do
        end do
      end do

      normalize = dble(sumi)/sumj
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            JI(nth,np,nr,nsa) = JI(nth,np,nr,nsa) * normalize
            JIR(nth,np,nr,nsa) = JIR(nth,np,nr,nsa) * normalize
          end do
        end do
      end do

    end do

  end subroutine calculate_jacobian

  subroutine construst_Xstg_as_pncp(xap)
    use fpcomm
    use fowcomm
    type(Xstg_as_pnc_point),intent(out) :: xap

    xap%number = 0

  end subroutine construst_Xstg_as_pncp

  subroutine add_pnc_orbit(xap, nr_in)
    use fpcomm
    use fowcomm
    type(Xstg_as_pnc_point),intent(inout) :: xap
    integer,intent(in) ::  nr_in
    type(Xstg_as_pnc_point) :: xap_old
    integer :: i

    xap_old = xap
    xap%number = xap%number+1

    if ( xap_old%number >= 1 ) then
      deallocate(xap%nr)
    end if

    allocate(xap%nr(xap%number))

    if ( xap_old%number >= 1 ) then
      do i = 1, xap_old%number
        xap%nr(i) = xap_old%nr(i)
      end do  
    end if

    xap%nr(xap%number) = nr_in

  end subroutine add_pnc_orbit

  function cal_spl1D(U,dx) result(cal)
    use fpcomm,only:rkind
    implicit none
    real(rkind) :: cal
    real(rkind),intent(in) :: U(4),dx
    integer :: k
    cal = 0.d0
    do k = 1, 4
      cal = cal+U(k)*dx**(k-1)
    end do
    
  end function

  subroutine newton_spl1D(U,a,x)
    ! U4*x**3+U3*x**2+U2*x+U1 = a
    use fpcomm,only:rkind
    implicit none
    real(rkind),intent(in) :: U(4),a
    real(rkind),intent(inout):: x
    integer :: i = 0,j
    real(rkind) :: x0,e = 1.d10,eps = 1.d-18,V(4),d

    V = U
    V(1) = V(1)-a

    do
      d = (V(4)*x**3+V(3)*x**2+V(2)*x+V(1))/(3.d0*V(4)*x**2+2.d0*V(3)*x+V(2))
      if(abs(d)<eps)exit

      x0 = x
      x = x0-d

      if(i>e)then
        write(*,*)"Do not converge at subroutine solve_spl1D"
        stop
      end if
    end do

  end subroutine newton_spl1D

end module fowprep
