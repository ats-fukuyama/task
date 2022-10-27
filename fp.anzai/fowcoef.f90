! fowcoef.f90
! [2022/3/30]
! *************************
!  Coefficient calculation
! *************************
! made by ota / modified by anzai
! ver.0.5

module fowcoef
  use fpcomm,only:rkind

  real(rkind),allocatable,dimension(:,:,:,:,:),private :: Dppl, Dptl, Fppl,&
                                                  Dtpl, Dttl, Fthl
  
  real(rkind),allocatable,private :: FNSBL(:,:,:,:,:)
  real(rkind),allocatable,private :: check_zeroD(:,:,:), check_zeroF(:,:)

contains

!=========================================
! Main program of coefficient calculation
!=========================================
  subroutine fow_coef
  !---------------------------------------------
  ! Diffusion Coef and Advection Vec calcul
  !---------------------------------------------

    use fpcomm
    use fowcomm
    use fowdistribution,only:convert_fI_to_fu
    use fpcalc
    use fpcoef
    use fpcalw

    implicit none

    integer :: nthp, nth, np, nr, nsa, i, j
    double precision :: begin_time, end_time

    call cpu_time(begin_time)

    allocate(check_zeroD(3,3,nsamax), check_zeroF(3,nsamax))

    allocate(FNSBL(nthmax,npmax,nrmax,nthpmax,nsamax))
    allocate(Dppl(nthmax,npmax+1,nrmax,nthpmax,nsamax), &
             Dptl(nthmax,npmax+1,nrmax,nthpmax,nsamax))
    allocate(Dtpl(nthmax+1,npmax,nrmax,nthpmax,nsamax), &
             Dttl(nthmax+1,npmax,nrmax,nthpmax,nsamax))
    allocate(Fppl(nthmax,npmax+1,nrmax,nthpmax,nsamax), &
             Fthl(nthmax+1,npmax,nrmax,nthpmax,nsamax))


    !**** Initialization
    Dppl(:,:,:,:,:) = 0.d0
    Dptl(:,:,:,:,:) = 0.d0
    Dtpl(:,:,:,:,:) = 0.d0
    Dttl(:,:,:,:,:) = 0.d0
    Fppl(:,:,:,:,:) = 0.d0
    Fthl(:,:,:,:,:) = 0.d0

    check_zeroF(:,:) = 0.d0
    check_zeroD(:,:,:) = 0.d0

    ! get local distribution function f(p,theta,r/a,thetap)
    call convert_fI_to_fu(FNSBL, FNSP)

    ! calculate coefficients
    MODELA = 0      ! bounce average is not executed in FP_CALX
    ! call FP_CALE    ! FEXX
    ! call FP_CALW    ! DWXX
    call FP_CALC    ! DCXX, FCXX

    !**** calculate local coefficient Dxxl(nth,np,nr,nthp,nsa)
    do nthp = 1, nthpmax

      !**** calculate back graound distribution FNSB(nth,np,nr,nsa) for each nthp
      do nsa = 1, nsamax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              FNSB(nth,np,nr,nsa) = FNSBL(nth,np,nr,nthp,nsa)
            end do
          end do
        end do
      end do

      !**** substitute Dxx at nthp to Dxxl
      do nsa = 1, nsamax
        do nr = 1, nrmax
          do np = 1, npmax+1
            do nth = 1, nthmax+1

              if ( nth /= nthmax+1 ) then 
                Dppl(nth,np,nr,nthp,nsa) = Dppl(nth,np,nr,nthp,nsa) + DCPP(nth,np,nr,nsa) ! + DWPP(nth,np,nr,nsa)
                Dptl(nth,np,nr,nthp,nsa) = Dptl(nth,np,nr,nthp,nsa) + DCPT(nth,np,nr,nsa) ! + DWPT(nth,np,nr,nsa)    
                Fppl(nth,np,nr,nthp,nsa) = Fppl(nth,np,nr,nthp,nsa) + FCPP(nth,np,nr,nsa) ! + FEPP(nth,np,nr,nsa)

                check_zeroD(1,1,nsa) = check_zeroD(1,1,nsa) + ABS( Dppl(nth,np,nr,nthp,nsa) )
                check_zeroD(1,2,nsa) = check_zeroD(1,2,nsa) + ABS( Dptl(nth,np,nr,nthp,nsa) )
                check_zeroF(1,nsa)   = check_zeroF(1,nsa)   + ABS( Fppl(nth,np,nr,nthp,nsa) )
              end if
              if ( np /= npmax+1 ) then
                Dtpl(nth,np,nr,nthp,nsa) = Dtpl(nth,np,nr,nthp,nsa) + DCTP(nth,np,nr,nsa) ! + DWTP(nth,np,nr,nsa) 
                Dttl(nth,np,nr,nthp,nsa) = Dttl(nth,np,nr,nthp,nsa) + DCTT(nth,np,nr,nsa) ! + DWTT(nth,np,nr,nsa)
                Fthl(nth,np,nr,nthp,nsa) = Fthl(nth,np,nr,nthp,nsa) + FCTH(nth,np,nr,nsa) ! + FETH(nth,np,nr,nsa)

                check_zeroD(2,1,nsa) = check_zeroD(2,1,nsa) + ABS( Dtpl(nth,np,nr,nthp,nsa) )
                check_zeroD(2,2,nsa) = check_zeroD(2,2,nsa) + ABS( Dttl(nth,np,nr,nthp,nsa) )
                check_zeroF(2,nsa)   = check_zeroF(2,nsa)   + ABS( Fthl(nth,np,nr,nthp,nsa) )

              end if

            end do
          end do
        end do
      end do

    end do

    call bounce_average

    deallocate(Dppl, Dptl, Fppl,Dtpl, Dttl, Fthl)
    deallocate(FNSBL, check_zeroD, check_zeroF)

    call cpu_time(end_time)
    write(6,'(A,ES10.3,A)')'fowcoef time : ',end_time-begin_time,'[sec]'

  end subroutine fow_coef

!================================
! Bounce average module
!================================
  subroutine bounce_average

    use fpcomm
    use fowcomm

    implicit none
    integer :: nth, np, nr, nsa, nthp, mode(3), nstp, nstpmax, ierr = 0
    real(rkind),dimension(3,3,max_stp) :: dIdu
    real(rkind),allocatable:: U_Dpp(:,:,:,:,:,:,:,:),&
                              U_Dpt(:,:,:,:,:,:,:,:),& 
                              U_Fpp(:,:,:,:,:,:,:,:),& 
                              U_Dtp(:,:,:,:,:,:,:,:),& 
                              U_Dtt(:,:,:,:,:,:,:,:),& 
                              U_Fth(:,:,:,:,:,:,:,:)   
    !** spline coefficient of Dppl in (\theta, \psi_p, \theta_p) plane
    real(rkind) :: Dpp_ob, Dpt_ob, Fpp_ob, Dtp_ob, Dtt_ob, Fth_ob, D_pls, D_mns, dt
    real(rkind) :: cpitch_ob, thetap_ob, psip_ob, JIl
    real(rkind) :: sumt


    allocate(U_Dpp(4,4,4,nthmax,nrmax,nthpmax,npmax+1,nsamax))
    allocate(U_Dpt(4,4,4,nthmax,nrmax,nthpmax,npmax+1,nsamax))
    allocate(U_Fpp(4,4,4,nthmax,nrmax,nthpmax,npmax+1,nsamax))
    allocate(U_Dtp(4,4,4,nthmax+1,nrmax,nthpmax,npmax,nsamax))
    allocate(U_Dtt(4,4,4,nthmax+1,nrmax,nthpmax,npmax,nsamax))
    allocate(U_Fth(4,4,4,nthmax+1,nrmax,nthpmax,npmax,nsamax))

    !**** Initialization
    Dppfow(:,:,:,:) = 0.d0
    Dptfow(:,:,:,:) = 0.d0
    Dprfow(:,:,:,:) = 0.d0
    Dtpfow(:,:,:,:) = 0.d0
    Dttfow(:,:,:,:) = 0.d0
    Dtrfow(:,:,:,:) = 0.d0
    Drpfow(:,:,:,:) = 0.d0
    Drtfow(:,:,:,:) = 0.d0
    Drrfow(:,:,:,:) = 0.d0
    Fppfow(:,:,:,:) = 0.d0
    Fthfow(:,:,:,:) = 0.d0
    Frrfow(:,:,:,:) = 0.d0

    do nsa = 1, nsamax
      if ( check_zeroD(1,1,nsa) >= 1.d-70 ) call make_U_Dxy(U_Dpp, Dppl, 'p', nsa)
      if ( check_zeroD(1,2,nsa) >= 1.d-70 ) call make_U_Dxy(U_Dpt, Dptl, 'p', nsa)
      if ( check_zeroD(2,1,nsa) >= 1.d-70 ) call make_U_Dxy(U_Dtp, Dtpl, 't', nsa)
      if ( check_zeroD(2,2,nsa) >= 1.d-70 ) call make_U_Dxy(U_Dtt, Dttl, 't', nsa)
      if ( check_zeroF(1,nsa) >= 1.d-70 )   call make_U_Dxy(U_Fpp, Fppl, 'p', nsa)
      if ( check_zeroF(2,nsa) >= 1.d-70 )   call make_U_Dxy(U_Fth, Fthl, 't', nsa)  
    end do

    !**** calculate Dpp, Dpt, Dpr, Fp
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax+1
          do nth = 1, nthmax

            if ( np == 1 ) then
              cycle
            end if

            nstpmax = orbit_p(nth,np,nr,nsa)%nstp_max

            mode = [0,1,0]
            call transformation_matrix(dIdu, orbit_p(nth,np,nr,nsa), nth, np, nr, nsa, mode)

            !**** time-integral over an orbit then divide poloidal period
            sumt = 0.d0
            do nstp = 2, nstpmax
              cpitch_ob = orbit_p(nth,np,nr,nsa)%costh(nstp)
              psip_ob   = orbit_p(nth,np,nr,nsa)%psip(nstp)
              thetap_ob = orbit_p(nth,np,nr,nsa)%thetap(nstp)
              
              !****[ calucurate local coefficient along orbit
              ! Dxx(p_orbit, theta_orbit, psip_orbit, thetap_orbit) ]
              call interpolate_D_unlessZero(Dpp_ob, U_Dpp(:,:,:,:,:,:,np,nsa), check_zeroD(1,1,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(Dpt_ob, U_Dpt(:,:,:,:,:,:,np,nsa), check_zeroD(1,2,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(Fpp_ob, U_Fpp(:,:,:,:,:,:,np,nsa), check_zeroF(1,nsa), cpitch_ob, psip_ob, thetap_ob)


              dt = orbit_p(nth,np,nr,nsa)%time(nstp)-orbit_p(nth,np,nr,nsa)%time(nstp-1)

              !**** Dxxfow = int_0^tau_p (integrand) dt
              Dppfow(nth,np,nr,nsa) = Dppfow(nth,np,nr,nsa)&
                                    + Dpp_ob * dt
              Dptfow(nth,np,nr,nsa) = Dptfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,2,nstp) + Dpt_ob * dIdu(2,2,nstp) ) *dt
              Dprfow(nth,np,nr,nsa) = Dprfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,3,nstp) + Dpt_ob * dIdu(2,3,nstp) ) *dt
              Fppfow(nth,np,nr,nsa) = Fppfow(nth,np,nr,nsa)&
                                    + Fpp_ob * dt

            end do

            if ( np == npmax+1 ) then
              JIl = JI(nth,np-1,nr,nsa)
            else if ( np == 1 ) then
              JIl = JI(nth,np,nr,nsa)
            else
              JIl = (JI(nth,np-1,nr,nsa)+JI(nth,np,nr,nsa))*0.50
            end if

            Dppfow(nth,np,nr,nsa) = Dppfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax) * JIl
            Dptfow(nth,np,nr,nsa) = Dptfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax) * JIl
            Dprfow(nth,np,nr,nsa) = Dprfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax) * JIl
            Fppfow(nth,np,nr,nsa) = Fppfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax) * JIl

          end do
        end do
      end do
    end do

    !**** calculate Dtp, Dtt, Dtr, Fth
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax+1

            if ( nth == 1 .or. nth == nthmax+1 .or. nth == nth_stg(nsa) ) then
              cycle
            end if

            nstpmax = orbit_th(nth,np,nr,nsa)%nstp_max

            mode = [1,0,0]
            call transformation_matrix(dIdu, orbit_th(nth,np,nr,nsa), nth, np, nr, nsa, mode)

            !**** time-integral over an orbit then divide poloidal period

            do nstp = 2, nstpmax
              cpitch_ob = orbit_th(nth,np,nr,nsa)%costh(nstp)
              psip_ob   = orbit_th(nth,np,nr,nsa)%psip(nstp)
              thetap_ob = orbit_th(nth,np,nr,nsa)%thetap(nstp)

              !****[ calucurate local coefficient along orbit
              ! Dxx(p_orbit, theta_orbit, psip_orbit, thetap_orbit) ]

              call interpolate_D_unlessZero(D_pls, U_Dpp(:,:,:,:,:,:,np+1,nsa), check_zeroD(1,1,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(D_mns, U_Dpp(:,:,:,:,:,:,np,nsa), check_zeroD(1,1,nsa), cpitch_ob, psip_ob, thetap_ob)
              Dpp_ob = (D_pls+D_mns)/2.d0

              call interpolate_D_unlessZero(D_pls, U_Dpt(:,:,:,:,:,:,np+1,nsa), check_zeroD(1,2,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(D_mns, U_Dpt(:,:,:,:,:,:,np,nsa), check_zeroD(1,2,nsa), cpitch_ob, psip_ob, thetap_ob)
              Dpt_ob = (D_pls+D_mns)/2.d0

              call interpolate_D_unlessZero(D_pls, U_Fpp(:,:,:,:,:,:,np+1,nsa), check_zeroF(1,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(D_mns, U_Fpp(:,:,:,:,:,:,np,nsa), check_zeroF(1,nsa), cpitch_ob, psip_ob, thetap_ob)
              Fpp_ob = (D_pls+D_mns)/2.d0


              call interpolate_D_unlessZero(Dtp_ob, U_Dtp(:,:,:,:,:,:,np,nsa), check_zeroD(2,1,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(Dtt_ob, U_Dtt(:,:,:,:,:,:,np,nsa), check_zeroD(2,2,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(Fth_ob, U_Fth(:,:,:,:,:,:,np,nsa), check_zeroF(2,nsa), cpitch_ob, psip_ob, thetap_ob)

              dt = orbit_th(nth,np,nr,nsa)%time(nstp)-orbit_th(nth,np,nr,nsa)%time(nstp-1)

              !**** Dxxfow = int_0^tau_p (integrand) dt
              Dtpfow(nth,np,nr,nsa) = Dtpfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,2,nstp) + Dpt_ob * dIdu(2,2,nstp) ) * dt

              Dttfow(nth,np,nr,nsa) = Dttfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,2,nstp)**2 + Dpt_ob * dIdu(1,2,nstp)*dIdu(2,2,nstp) &
                                    + Dtp_ob * dIdu(1,2,nstp)*dIdu(2,2,nstp) + Dtt_ob * dIdu(2,2,nstp)**2 ) * dt

              Dtrfow(nth,np,nr,nsa) = Dtrfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,2,nstp)*dIdu(1,3,nstp) + Dpt_ob * dIdu(1,2,nstp)*dIdu(2,3,nstp) &
                                    + Dtp_ob * dIdu(2,2,nstp)*dIdu(1,3,nstp) + Dtt_ob * dIdu(2,2,nstp)*dIdu(2,3,nstp)) * dt

              Fthfow(nth,np,nr,nsa) = Fthfow(nth,np,nr,nsa)&
                                    + ( Fpp_ob*dIdu(1,2,nstp) + Fth_ob*dIdu(2,2,nstp) ) * dt                      
            end do

            if ( nth == nth_pnc(nsa) .and. theta_pnc(np,nr,nsa) /= NO_PINCH_ORBIT ) then
              ! JIl = ( JI(nth-1,np,nr,nsa)*(1.d0 - IBCflux_ratio(np,nr,nsa))&
              !       + JI(nth,np,nr,nsa)*IBCflux_ratio(np,nr,nsa))
              ! JIl = ( JI(nth-1,np,nr,nsa)*IBCflux_ratio(np,nr,nsa)&
              !       * orbit_m(nth-1,np,nr,nsa)%time(orbit_m(nth-1,np,nr,nsa)%nstp_max)&
              !       + JI(nth,np,nr,nsa)*(1.d0 - IBCflux_ratio(np,nr,nsa)) &
              !       * orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max))&
              !       / (orbit_m(nth-1,np,nr,nsa)%time(orbit_m(nth-1,np,nr,nsa)%nstp_max) &
              !       + orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max))
              ! JIl = ( JI(nth-1,np,nr,nsa)*IBCflux_ratio(np,nr,nsa)&
              !       * orbit_m(nth-1,np,nr,nsa)%time(orbit_m(nth-1,np,nr,nsa)%nstp_max)&
              !       + JI(nth,np,nr,nsa)*(1.d0 - IBCflux_ratio(np,nr,nsa)) &
              !       * orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max))&
              !       * orbit_m(nth,np,nr,nsa)%time(nstpmax)
              ! JIl = (JI(nth-1,np,nr,nsa)*orbit_m(nth-1,np,nr,nsa)%time(orbit_m(nth-1,np,nr,nsa)%nstp_max)&
              !       +JI(nth,np,nr,nsa)*orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max))&
              !       *0.5d0/orbit_m(nth,np,nr,nsa)%time(nstpmax)
              ! JIl = (JI(nth-1,np,nr,nsa)&
              !       +JI(nth,np,nr,nsa))&
              !       *0.5d0
              JIl = (JI(nth-1,np,nr,nsa)*orbit_m(nth-1,np,nr,nsa)%time(orbit_m(nth-1,np,nr,nsa)%nstp_max)&
                    +JI(nth,np,nr,nsa)*orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max))&
                    /(orbit_m(nth-1,np,nr,nsa)%time(orbit_m(nth-1,np,nr,nsa)%nstp_max)+ &
                    orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max))
              ! JIl = (JI(nth-1,np,nr,nsa)*orbit_m(nth-1,np,nr,nsa)%time(orbit_m(nth-1,np,nr,nsa)%nstp_max)&
              !       +JI(nth,np,nr,nsa)*orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max))&
              !       /(orbit_m(nth-1,np,nr,nsa)%time(orbit_m(nth-1,np,nr,nsa)%nstp_max)+ &
              !       orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max))
              ! JIl = (JI(nth-1,np,nr,nsa)*orbit_m(nth-1,np,nr,nsa)%time(orbit_m(nth-1,np,nr,nsa)%nstp_max)&
              !       +JI(nth,np,nr,nsa)*orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max))&
              !       *0.5d0*orbit_th(nth,np,nr,nsa)%time(nstpmax)

            else if ( nth == 1 ) then
              JIl = JI(nth,np,nr,nsa)
            else if ( nth == nthmax+1 ) then
              JIl = JI(nth-1,np,nr,nsa)  
            else 
              JIl = (JI(nth,np,nr,nsa)+JI(nth-1,np,nr,nsa))*0.5d0
            end if

            Dtpfow(nth,np,nr,nsa) = Dtpfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax) * JIl
            Dttfow(nth,np,nr,nsa) = Dttfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax) * JIl
            Dtrfow(nth,np,nr,nsa) = Dtrfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax) * JIl      
            Fthfow(nth,np,nr,nsa) = Fthfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax) * JIl

          end do
        end do
      end do
    end do

    !**** calculate Drp, Drt, Drr, Frr
    do nsa = 1, nsamax
      do nr = 1, nrmax+1
        do np = 1, npmax
          do nth = 1, nthmax

            if ( nr == 1 ) then
              cycle    
            end if

            nstpmax = orbit_r(nth,np,nr,nsa)%nstp_max

            mode = [0,0,1]
            call transformation_matrix(dIdu, orbit_r(nth,np,nr,nsa), nth, np, nr, nsa, mode)

            !**** time-integral over an orbit then divide poloidal period

            do nstp = 2, nstpmax
              cpitch_ob = orbit_r(nth,np,nr,nsa)%costh(nstp)
              psip_ob   = orbit_r(nth,np,nr,nsa)%psip(nstp)
              thetap_ob = orbit_r(nth,np,nr,nsa)%thetap(nstp)

              !****[ calucurate local coefficient along orbit
              ! Dxx(p_orbit, theta_orbit, psip_orbit, thetap_orbit) ]
              call interpolate_D_unlessZero(D_pls, U_Dpp(:,:,:,:,:,:,np+1,nsa), check_zeroD(1,1,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(D_mns, U_Dpp(:,:,:,:,:,:,np,nsa), check_zeroD(1,1,nsa), cpitch_ob, psip_ob, thetap_ob)
              Dpp_ob = (D_pls+D_mns)/2.d0

              call interpolate_D_unlessZero(D_pls, U_Dpt(:,:,:,:,:,:,np+1,nsa), check_zeroD(1,2,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(D_mns, U_Dpt(:,:,:,:,:,:,np,nsa), check_zeroD(1,2,nsa), cpitch_ob, psip_ob, thetap_ob)
              Dpt_ob = (D_pls+D_mns)/2.d0

              call interpolate_D_unlessZero(D_pls, U_Fpp(:,:,:,:,:,:,np+1,nsa), check_zeroF(1,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(D_mns, U_Fpp(:,:,:,:,:,:,np,nsa), check_zeroF(1,nsa), cpitch_ob, psip_ob, thetap_ob)
              Fpp_ob = (D_pls+D_mns)/2.d0


              call interpolate_D_unlessZero(Dtp_ob, U_Dtp(:,:,:,:,:,:,np,nsa), check_zeroD(2,1,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(Dtt_ob, U_Dtt(:,:,:,:,:,:,np,nsa), check_zeroD(2,2,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(Fth_ob, U_Fth(:,:,:,:,:,:,np,nsa), check_zeroF(2,nsa), cpitch_ob, psip_ob, thetap_ob)

              dt = orbit_r(nth,np,nr,nsa)%time(nstp)-orbit_r(nth,np,nr,nsa)%time(nstp-1)

              !**** Dxxfow = int_0^tau_p (integrand) dt
              Drpfow(nth,np,nr,nsa) = Drpfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,3,nstp) + Dpt_ob * dIdu(2,3,nstp) ) * dt

              Drtfow(nth,np,nr,nsa) = Drtfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,2,nstp)*dIdu(1,3,nstp) + Dpt_ob * dIdu(1,3,nstp)*dIdu(2,2,nstp) &
                                    + Dtp_ob * dIdu(1,2,nstp)*dIdu(2,3,nstp) + Dtt_ob * dIdu(2,2,nstp)*dIdu(2,3,nstp)) * dt

              Drrfow(nth,np,nr,nsa) = Drrfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,3,nstp)**2 + Dpt_ob * dIdu(1,3,nstp)*dIdu(2,3,nstp) &
                                    + Dtp_ob * dIdu(2,3,nstp)*dIdu(1,3,nstp) + Dtt_ob * dIdu(2,3,nstp)**2) * dt

              Frrfow(nth,np,nr,nsa) = Frrfow(nth,np,nr,nsa)&
                                    + ( Fpp_ob*dIdu(1,3,nstp) + Fth_ob*dIdu(2,3,nstp) ) * dt
                      
            end do

            if ( nr == 1 ) then
              JIl = JI(nth,np,nr,nsa)
            else if ( nr == nrmax+1 ) then
              JIl = JI(nth,np,nr-1,nsa)
            else
              JIl = (JI(nth,np,nr,nsa)+JI(nth,np,nr-1,nsa))*0.5d0
            end if

            Drpfow(nth,np,nr,nsa) = Drpfow(nth,np,nr,nsa) / orbit_r(nth,np,nr,nsa)%time(nstpmax) * JIl
            Drtfow(nth,np,nr,nsa) = Drtfow(nth,np,nr,nsa) / orbit_r(nth,np,nr,nsa)%time(nstpmax) * JIl
            Drrfow(nth,np,nr,nsa) = Drrfow(nth,np,nr,nsa) / orbit_r(nth,np,nr,nsa)%time(nstpmax) * JIl
            Frrfow(nth,np,nr,nsa) = Frrfow(nth,np,nr,nsa) / orbit_r(nth,np,nr,nsa)%time(nstpmax) * JIl

          end do
        end do
      end do
    end do
    

  end subroutine bounce_average

!===============================================
! For matrix calculation modules
!===============================================
  subroutine transformation_matrix(dIdu, ob, nth, np, nr, nsa, mode)
  !----------------------------------------------------------------
  ! Calcuration of transform matrix  from Energy, momentum space
  ! to momentum, pitch angle, the largest magnetic flux
  !----------------------------------------------------------------
    use fpcomm
    use fowcomm
    use foworbit

    implicit none

    real(rkind),dimension(3,3,max_stp),intent(out) :: dIdu
    type(orbit),intent(in) :: ob
    integer,intent(in) :: nth, np, nr, nsa, mode(3)

    ! elements of transformation matrix, dIdu
    real(rkind) :: dpdp,  dthmdp,  drmdp,&
                   dpdth, dthmdth, drmdth,&
                   dpdr,  dthmdr,  drmdr
    real(rkind) :: pl, Bml, dBmdrl, Fml, dFmdrl, cthm, sthm, dpsimdrl
    real(rkind) :: Bob, Fob, cthob, sthob, dBobdr, dpsiobdr, dFobdr
    real(rkind) :: A(2,2), b(2), detA
    integer :: nstp, nstpmax

    select case(mode(1))
    case(0)
      if ( mode(2) == 0 .and. mode(3) == 0 ) cthm = COS( thetam(nth,np,nr,nsa) )
      if ( mode(2) == 1 .and. mode(3) == 0 ) cthm = COS( thetam_pg(nth,np,nr,nsa) )
      if ( mode(2) == 0 .and. mode(3) == 1 ) cthm = COS( thetam_rg(nth,np,nr,nsa) )
    case(1)
      cthm = COS( thetamg(nth,np,nr,nsa) )
    end select
    sthm = SQRT( 1.d0-cthm**2 )

    select case(mode(2))
    case(0)
      pl = pm(np,nsa)*ptfp0(nsa)
    case(1)
      pl = pg(np,nsa)*ptfp0(nsa)
    end select

    select case(mode(3))
    case(0)
      if ( cthm*aefp(nsa) >= 0.d0 ) then
        Bml = Bout(nr)
        dBmdrl = dBoutdr(nr)
      else
        Bml = Bin(nr)
        dBmdrl = dBindr(nr)
      end if
      dpsimdrl = dpsimdr(nr)
      Fml = Fpsi(nr)
      dFmdrl = dFdr(nr)
    case(1)
      if ( cthm*aefp(nsa) >= 0.d0 ) then
        Bml = Boutg(nr)
        dBmdrl = dBoutgdr(nr)
      else
        Bml = Bing(nr)
        dBmdrl = dBingdr(nr)
      end if
      dpsimdrl = dpsimgdr(nr)
      Fml = Fpsig(nr)
      dFmdrl = dFgdr(nr)
    end select

    A(1,1) = 2.d0*sthm*cthm/Bml
    A(1,2) = -1.d0*sthm**2*dBmdrl/Bml**2
    A(2,1) = Fml/Bml*pl*sthm
    A(2,2) = aefp(nsa)*dpsimdrl-(dFmdrl*Bml-Fml*dBmdrl)/Bml**2*pl*cthm
    detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)

    nstpmax = ob%nstp_max

    if ( detA /= 0.d0 ) then
      do nstp = 1, nstpmax
        Bob = ob%Babs(nstp)
        Fob = ob%F(nstp)
        cthob = ob%costh(nstp)
        sthob = ob%sinth(nstp)
        dBobdr = ob%dBdr(nstp)
        dpsiobdr = ob%dpsipdr(nstp)
        dFobdr = ob%dFdr(nstp)

        ! dX/dp
        b(1) = 0.d0
        b(2) = Fml/Bml*cthm-Fob/Bob*cthob
        dpdp   = 1.d0
        dthmdp = (A(2,2)*b(1)-A(1,2)*b(2))/detA
        drmdp  = (A(1,1)*b(2)-A(2,1)*b(1))/detA

        ! dX/dtheta
        b(1) = 2.d0*sthob*cthob/Bob
        b(2) = Fob/Bob*pl*sthob
        dpdth   = 0.d0
        dthmdth = (A(2,2)*b(1)-A(1,2)*b(2))/detA
        drmdth  = (A(1,1)*b(2)-A(2,1)*b(1))/detA

        ! dX/drho
        b(1) = -1.d0*sthob**2/Bob**2*dBobdr
        b(2) = aefp(nsa)*dpsiobdr - ( dFobdr*Bob-Fob*dBobdr )/Bob**2*pl*cthob
        dpdr   = 0.d0
        dthmdr = (A(2,2)*b(1)-A(1,2)*b(2))/detA
        drmdr  = (A(1,1)*b(2)-A(2,1)*b(1))/detA

        dIdu(1,1,nstp) = dpdp
        dIdu(1,2,nstp) = dthmdp*ptfp0(nsa)
        dIdu(1,3,nstp) = drmdp*ptfp0(nsa)
        dIdu(2,1,nstp) = dpdth
        dIdu(2,2,nstp) = dthmdth
        dIdu(2,3,nstp) = drmdth
        dIdu(3,1,nstp) = dpdr
        dIdu(3,2,nstp) = dthmdr
        dIdu(3,3,nstp) = drmdr

      end do

    else
      do nstp = 1, nstpmax
        dIdu(1,1,nstp) = 1.d0
        dIdu(1,2,nstp) = 0.d0
        dIdu(1,3,nstp) = 0.d0
        dIdu(2,1,nstp) = 0.d0
        dIdu(2,2,nstp) = 1.d0
        dIdu(2,3,nstp) = 0.d0
        dIdu(3,1,nstp) = 0.d0
        dIdu(3,2,nstp) = 0.d0
        dIdu(3,3,nstp) = 1.d0
      end do

    end if

    ! do nstp = 1, nstpmax
    !   write(6,'(I4)')nstp
    !   write(6,'(3ES12.4)')dIdu(1,1,nstp),dIdu(1,2,nstp),dIdu(1,3,nstp)
    !   write(6,'(3ES12.4)')dIdu(2,1,nstp),dIdu(2,2,nstp),dIdu(2,3,nstp)
    !   write(6,'(3ES12.4)')dIdu(3,1,nstp),dIdu(3,2,nstp),dIdu(3,3,nstp)
    ! end do

  end subroutine transformation_matrix

  subroutine make_U_Dxy(U_Dxy, Dxyl, x, nsa)
    use fpcomm
    use fowcomm
    USE libspl3d

    implicit none

    real(rkind),intent(out) :: U_Dxy(:,:,:,:,:,:,:,:)
    real(rkind),intent(in) :: Dxyl(:,:,:,:,:)
    character(*),intent(in) :: x
    intent(in) :: nsa
    real(rkind),allocatable :: Dxyl_tmp(:,:,:), U_Dxy_tmp(:,:,:,:,:,:) &
                              , FX(:,:,:), FY(:,:,:), FZ(:,:,:), FXY(:,:,:), FYZ(:,:,:) &
                              , FZX(:,:,:), FXYZ(:,:,:), Xtmp(:)

    integer :: nth, np, nr, nsa, nthp, i, j ,k 
    integer :: nxmax, nymax, nzmax, p, t, ierr = 0

    if ( x == 'p' ) then
      p = 1
      t = 0

      allocate(Xtmp(nthmax))
      do nth = 1, nthmax
        Xtmp(nth) = cosm(nth)
      end do
        
    else if ( x == 't' ) then
      p = 0
      t = 1

      allocate(Xtmp(nthmax+1))
      do nth = 1, nthmax+1
        Xtmp(nth) = cosg(nth)
      end do

    else if ( x == 'm' ) then
      p = 0
      t = 0

      allocate(Xtmp(nthmax))
      do nth = 1, nthmax
        Xtmp(nth) = cosm(nth)
      end do

    end if

    allocate( Dxyl_tmp(nthmax+t,nrmax,nthpmax), U_Dxy_tmp(4,4,4,nthmax+t,nrmax,nthpmax) )
    allocate( FX(nthmax+t,nrmax,nthpmax), FY(nthmax+t,nrmax,nthpmax), FZ(nthmax+t,nrmax,nthpmax), FXY(nthmax+t,nrmax,nthpmax) )
    allocate( FZX(nthmax+t,nrmax,nthpmax), FYZ(nthmax+t,nrmax,nthpmax), FXYZ(nthmax+t,nrmax,nthpmax) )

    do np = 1, npmax+p

      do nthp = 1, nthpmax
        do nr = 1, nrmax
          do nth = 1, nthmax+t
            Dxyl_tmp(nth,nr,nthp) = Dxyl(nth,np,nr,nthp,nsa)
          end do
        end do
      end do

      call SPL3D(Xtmp,psim,theta_p,Dxyl_tmp,FX,FY,FZ,FXY,FYZ,FZX,FXYZ,&
                U_Dxy_tmp,nthmax+t,nrmax,nthmax+t,nrmax,nthpmax,0,0,0,IERR)

      do nthp = 1, nthpmax
        do nr = 1, nrmax
          do nth = 1, nthmax+t

            do k = 1, 4
              do j = 1, 4
                do i = 1, 4
                  U_Dxy(i,j,k,nth,nr,nthp,np,nsa) = U_Dxy_tmp(i,j,k,nth,nr,nthp)
                end do
              end do
            end do

          end do
        end do
      end do

    end do

  end subroutine make_U_Dxy

  subroutine interpolate_D_unlessZero(C_out, U, check0, cpitch_in, psip_in, thetap_in)
    use fpcomm
    use fowcomm
    USE libspl3d
    implicit none

    real(rkind),intent(out) :: C_out
    real(rkind),intent(in) :: U(:,:,:,:,:,:), check0, cpitch_in, psip_in, thetap_in
    integer :: nxmax, nymax, nzmax, ierr
    ierr = 0

    nxmax = size(U,4)
    nymax = size(U,5)
    nzmax = size(U,6)

    if ( check0 < 1.d-70 ) then
      C_out = 0.d0
      return
    end if

    if ( nxmax == nthmax .and. nymax == nrmax ) then
      call SPL3DF(cpitch_in,psip_in,thetap_in,C_out,cosm,psim,theta_p&
                  ,U,nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)

    else if (nxmax == nthmax+1 .and. nymax == nrmax ) then
      call SPL3DF(cpitch_in,psip_in,thetap_in,C_out,cosg,psim,theta_p&
                  ,U,nthmax+1,nrmax,nthmax+1,nrmax,nthpmax,IERR)

    else if (nxmax == nthmax .and. nymax == nrmax+1 ) then
      call SPL3DF(cpitch_in,psip_in,thetap_in,C_out,cosm,psimg,theta_p&
                  ,U,nthmax,nrmax+1,nthmax,nrmax+1,nthpmax,IERR)
                                        
    end if


  end subroutine interpolate_D_unlessZero

end module fowcoef
