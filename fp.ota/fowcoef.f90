module fowcoef
  use fpcomm,only:rkind

  private

  public :: fow_coef

  real(rkind),allocatable,dimension(:,:,:,:,:) :: Dppl, Dptl, Fppl,&
                                                  Dtpl, Dttl, Fthl
  
  real(rkind),allocatable :: FNSBL(:,:,:,:,:)
  real(rkind),allocatable :: dBmdpsi(:,:), dFdpsi(:), dBmgdpsi(:,:), dFgdpsi(:), psim0(:), psimg0(:)
  real(rkind),allocatable :: check_zeroD(:,:,:), check_zeroF(:,:)

contains

  subroutine fow_coef

    use fpcomm
    use fowcomm
    use fowdistribution,only:convert_fI_to_fu
    use fpcalc
    use fpcoef
    use fpcalw

    use fpwrite

    implicit none

    integer :: nthp, nth, np, nr, nsa, i, j


    allocate(check_zeroD(3,3,nsamax), check_zeroF(3,nsamax))

    allocate(FNSBL(nthmax,npmax,nrmax,nthpmax,nsamax))
    allocate(Dppl(nthmax,npmax+1,nrmax,nthpmax,nsamax),Dptl(nthmax,npmax+1,nrmax,nthpmax,nsamax))
    allocate(Dtpl(nthmax+1,npmax,nrmax,nthpmax,nsamax),Dttl(nthmax+1,npmax,nrmax,nthpmax,nsamax))
    allocate(Fppl(nthmax,npmax+1,nrmax,nthpmax,nsamax),Fthl(nthmax+1,npmax,nrmax,nthpmax,nsamax))

    do nsa = 1, nsamax
      do nthp = 1, nthpmax
        do nr = 1, nrmax
          do np = 1, npmax+1
            do nth = 1, nthmax+1
              if ( nth <= nthmax ) then
                Dppl(nth,np,nr,nthp,nsa) = 0.d0
                Dptl(nth,np,nr,nthp,nsa) = 0.d0
                Fppl(nth,np,nr,nthp,nsa) = 0.d0  
              end if
              if ( np <= npmax ) then
                Dtpl(nth,np,nr,nthp,nsa) = 0.d0
                Dttl(nth,np,nr,nthp,nsa) = 0.d0
                Fthl(nth,np,nr,nthp,nsa) = 0.d0  
              end if
            end do
          end do
        end do
      end do  
    end do

    ! get local distribution function f(p,theta,r/a,thetap)
    call convert_fI_to_fu(FNSBL, FNSP)


    ! calculate local coefficient Dxxl(nth,np,nr,nthp,nsa)
    do nthp = 1, nthpmax

      ! calculate back graound distribution FNSB(nth,np,nr,nsa) for each nthp
      do nsa = 1, nsamax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              FNSB(nth,np,nr,nsa) = FNSBL(nth,np,nr,nthp,nsa)
            end do
          end do
        end do
      end do

      ! calculate coefficients
      MODELA = 0      ! bounce average is not executed in FP_CALX
      ! call FP_CALE    ! FEXX
      ! call FP_CALW    ! DWXX
      call FP_CALC    ! DCXX, FCXX

      ! substitute Dxx at nthp to Dxxl
      do nsa = 1, nsamax
        do j = 1, 3
          check_zeroF(j,nsa) = 0.d0
          do i = 1, 3
            check_zeroD(i,j,nsa) = 0.d0
          end do
        end do  
      end do

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
                ! write(*,*)"check zero",check_zeroD(1,1,nsa),check_zeroD(1,2,nsa)
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

    open(100,file="txt/fow_DPPL.txt")
    open(101,file="txt/fow_DTTL.txt")
    open(102,file="txt/fow_DPTL.txt")
    open(103,file="txt/fow_DTPL.txt")
    do nsa = 1, nsamax
      do nthp = 1, nthpmax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              write(100,'(A,5I3,ES12.4)')"DPPL",nth,np,nr,nthp,nsa,DPPL(nth,np,nr,nthp,nsa)
              write(101,'(A,5I3,ES12.4)')"DTTL",nth,np,nr,nthp,nsa,DTTL(nth,np,nr,nthp,nsa)
              write(102,'(A,5I3,ES12.4)')"DPTL",nth,np,nr,nthp,nsa,DPTL(nth,np,nr,nthp,nsa)
              write(103,'(A,5I3,ES12.4)')"DTPL",nth,np,nr,nthp,nsa,DTPL(nth,np,nr,nthp,nsa)
            end do
          end do
        end do
      end do
    end do
    close(100)
    close(101)
    close(102)
    close(103)

    call bounce_average

    deallocate(Dppl, Dptl, Fppl,Dtpl, Dttl, Fthl)
    deallocate(FNSBL, check_zeroD, check_zeroF)

  end subroutine fow_coef

  subroutine bounce_average

    use fpcomm
    use fowcomm

    implicit none
    integer :: nth, np, nr, nsa, nthp, mode(3), nstp, nstpmax, ierr = 0
    real(rkind),allocatable :: dIdu(:,:,:)
    real(rkind),allocatable:: U_Dpp(:,:,:,:,:,:,:,:),& ! spline coefficient of Dppl in (\theta, \psi_p, \theta_p) plane
                              U_Dpt(:,:,:,:,:,:,:,:),& 
                              U_Fpp(:,:,:,:,:,:,:,:),& 
                              U_Dtp(:,:,:,:,:,:,:,:),& 
                              U_Dtt(:,:,:,:,:,:,:,:),& 
                              U_Fth(:,:,:,:,:,:,:,:)   
    real(rkind) :: Dpp_ob, Dpt_ob, Fpp_ob, Dtp_ob, Dtt_ob, Fth_ob, D_pls, D_mns, dt
    real(rkind) :: cpitch_ob, thetap_ob, psip_ob, J_I
    real(rkind) :: sumt

    if ( .not.allocated(psim0) ) then
      allocate(psim0(nrmax), psimg0(nrmax+1))
      do nr = 1, nrmax
        psim0(nr) = psim(nr)*psi0
        psimg0(nr) = psimg(nr)*psi0
      end do
      psimg0(nrmax+1) = psimg(nrmax+1)*psi0
    end if

    if ( .not.allocated(dBmdpsi) ) then
      allocate(dBmdpsi(nrmax,2), dFdpsi(nrmax))
      allocate(dBmgdpsi(nrmax+1,2), dFgdpsi(nrmax+1))
      
      call first_order_derivative(dFdpsi, Fpsi, psim0)
      call first_order_derivative(dBmdpsi(:,1), Bout, psim0)
      call first_order_derivative(dBmdpsi(:,2), Bin, psim0)
      call first_order_derivative(dFgdpsi, Fpsig, psimg0)
      call first_order_derivative(dBmgdpsi(:,1), Boutg, psimg0)
      call first_order_derivative(dBmgdpsi(:,2), Bing, psimg0)
  
    end if

    allocate(U_Dpp(4,4,4,nthmax,nrmax,nthpmax,npmax+1,nsamax))
    allocate(U_Dpt(4,4,4,nthmax,nrmax,nthpmax,npmax+1,nsamax))
    allocate(U_Fpp(4,4,4,nthmax,nrmax,nthpmax,npmax+1,nsamax))
    allocate(U_Dtp(4,4,4,nthmax+1,nrmax,nthpmax,npmax,nsamax))
    allocate(U_Dtt(4,4,4,nthmax+1,nrmax,nthpmax,npmax,nsamax))
    allocate(U_Fth(4,4,4,nthmax+1,nrmax,nthpmax,npmax,nsamax))


    do nsa = 1, nsamax
      if ( check_zeroD(1,1,nsa) >= 1.d-50 ) call make_U_Dxy(U_Dpp, Dppl, 'p', nsa)
      if ( check_zeroD(1,2,nsa) >= 1.d-50 ) call make_U_Dxy(U_Dpt, Dptl, 'p', nsa)
      if ( check_zeroD(2,1,nsa) >= 1.d-50 ) call make_U_Dxy(U_Dtp, Dtpl, 't', nsa)
      if ( check_zeroD(2,2,nsa) >= 1.d-50 ) call make_U_Dxy(U_Dtt, Dttl, 't', nsa)
      if ( check_zeroF(1,nsa) >= 1.d-50 )   call make_U_Dxy(U_Fpp, Fppl, 'p', nsa)
      if ( check_zeroF(2,nsa) >= 1.d-50 )   call make_U_Dxy(U_Fth, Fthl, 't', nsa)  
    end do

    ! calculate Dpp, Dpt, Dpr, Fp
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax+1
          do nth = 1, nthmax

            Dppfow(nth,np,nr,nsa) = 0.d0
            Dptfow(nth,np,nr,nsa) = 0.d0
            Dprfow(nth,np,nr,nsa) = 0.d0
            Fppfow(nth,np,nr,nsa) = 0.d0

            if ( np == 1 ) then
              cycle
            end if

            nstpmax = orbit_p(nth,np,nr,nsa)%nstp_max

            allocate(dIdu(3,3,nstpmax))

            mode = [0,1,0]
            call transformation_matrix(dIdu, orbit_p(nth,np,nr,nsa), nth, np, nr, nsa, mode)

            ! time-integral over an orbit then divide poloidal period 
            sumt = 0.d0
            do nstp = 2, nstpmax
              cpitch_ob = cos(orbit_p(nth,np,nr,nsa)%theta(nstp))
              psip_ob   = orbit_p(nth,np,nr,nsa)%psip(nstp)
              thetap_ob = orbit_p(nth,np,nr,nsa)%thetap(nstp)
              
              ! calucurate local coefficient along orbit Dxx(p_orbit, theta_orbit, psip_orbit, thetap_orbit)
              call interpolate_D_unlessZero(Dpp_ob, U_Dpp(:,:,:,:,:,:,np,nsa), check_zeroD(1,1,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(Dpt_ob, U_Dpt(:,:,:,:,:,:,np,nsa), check_zeroD(1,2,nsa), cpitch_ob, psip_ob, thetap_ob)
              call interpolate_D_unlessZero(Fpp_ob, U_Fpp(:,:,:,:,:,:,np,nsa), check_zeroF(1,nsa), cpitch_ob, psip_ob, thetap_ob)


              dt = orbit_p(nth,np,nr,nsa)%time(nstp)-orbit_p(nth,np,nr,nsa)%time(nstp-1)
              sumt = sumt + dt

              ! Dxxfow = int_0^tau_p (integrand) dt
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
              J_I = Jacobian_I(nth,np-1,nr,nsa)
            else if ( np == 1 ) then
              J_I = Jacobian_I(nth,np,nr,nsa)
            else
              J_I = (Jacobian_I(nth,np-1,nr,nsa)+Jacobian_I(nth,np,nr,nsa))*0.50
            end if

            Dppfow(nth,np,nr,nsa) = Dppfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax) * J_I
            Dptfow(nth,np,nr,nsa) = Dptfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax) * J_I
            Dprfow(nth,np,nr,nsa) = Dprfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax) * J_I
            Fppfow(nth,np,nr,nsa) = Fppfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax) * J_I
            ! Dppfow(nth,np,nr,nsa) = Dppfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax) * J_I
            ! Dptfow(nth,np,nr,nsa) = 0.d0
            ! Dprfow(nth,np,nr,nsa) = 0.d0
            ! Fppfow(nth,np,nr,nsa) = Fppfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax) * J_I

            deallocate(dIdu)

          end do
        end do
      end do
    end do

    ! calculate Dtp, Dtt, Dtr, Fth
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax+1

            Dtpfow(nth,np,nr,nsa) = 0.d0
            Dttfow(nth,np,nr,nsa) = 0.d0
            Dtrfow(nth,np,nr,nsa) = 0.d0
            Fthfow(nth,np,nr,nsa) = 0.d0

            if ( nth == 1 .or. nth == nthmax+1 .or. nth == nth_stg(nsa) ) then
              cycle
            end if

            nstpmax = orbit_th(nth,np,nr,nsa)%nstp_max

            allocate(dIdu(3,3,nstpmax))

            mode = [1,0,0]
            call transformation_matrix(dIdu, orbit_th(nth,np,nr,nsa), nth, np, nr, nsa, mode)

            ! time-integral over an orbit then divide poloidal period 

            do nstp = 2, nstpmax
              cpitch_ob = cos(orbit_th(nth,np,nr,nsa)%theta(nstp))
              psip_ob   = orbit_th(nth,np,nr,nsa)%psip(nstp)
              thetap_ob = orbit_th(nth,np,nr,nsa)%thetap(nstp)

              ! calucurate local coefficient along orbit Dxx(p_orbit, theta_orbit, psip_orbit, thetap_orbit)

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

              ! Dxxfow = int_0^tau_p (integrand) dt
              Dtpfow(nth,np,nr,nsa) = Dtpfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,2,nstp) + Dpt_ob * dIdu(2,2,nstp) ) * dt

              Dttfow(nth,np,nr,nsa) = Dttfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,2,nstp)**2 + Dpt_ob * dIdu(1,2,nstp)*dIdu(2,2,nstp) &
                                    + Dtp_ob * dIdu(1,2,nstp)*dIdu(2,2,nstp) + Dtt_ob * dIdu(2,2,nstp)**2 ) * dt
              ! write(*,'(A,4ES12.4)')"loop",Dpp_ob,Dpt_ob,Dtp_ob,Dtt_ob
              ! write(*,'(A,3ES12.4)')"loop2",dIdu(1,2,nstp),dIdu(2,2,nstp),dIdu(1,3,nstp)
              Dtrfow(nth,np,nr,nsa) = Dtrfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,2,nstp)*dIdu(1,3,nstp) + Dpt_ob * dIdu(1,2,nstp)*dIdu(2,3,nstp) &
                                    + Dtp_ob * dIdu(2,2,nstp)*dIdu(1,3,nstp) + Dtt_ob * dIdu(2,2,nstp)*dIdu(2,3,nstp)) * dt

              Fthfow(nth,np,nr,nsa) = Fthfow(nth,np,nr,nsa)&
                                    + ( Fpp_ob*dIdu(1,2,nstp) + Fth_ob*dIdu(2,2,nstp) ) * dt                      
            end do

            if ( nth == nth_pnc(nsa) .and. theta_pnc(np,nr,nsa) /= NO_PINCH_ORBIT ) then
              J_I = (Jacobian_I(nth-1,np,nr,nsa)*orbit_m(nth-1,np,nr,nsa)%time(orbit_m(nth-1,np,nr,nsa)%nstp_max)&
                    +Jacobian_I(nth,np,nr,nsa)*orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max))&
                    *0.5d0*orbit_th(nth,np,nr,nsa)%time(nstpmax)

            else if ( nth == 1 ) then
              J_I = Jacobian_I(nth,np,nr,nsa)
            else if ( nth == nthmax+1 ) then
              J_I = Jacobian_I(nth-1,np,nr,nsa)  
            else 
              J_I = (Jacobian_I(nth,np,nr,nsa)+Jacobian_I(nth-1,np,nr,nsa))*0.5d0
            end if

            Dtpfow(nth,np,nr,nsa) = Dtpfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax) * J_I
            Dttfow(nth,np,nr,nsa) = Dttfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax) * J_I
            Dtrfow(nth,np,nr,nsa) = Dtrfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax) * J_I      
            Fthfow(nth,np,nr,nsa) = Fthfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax) * J_I
            ! Dtpfow(nth,np,nr,nsa) = 0.d0
            ! Dttfow(nth,np,nr,nsa) = Dttfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax) * J_I
            ! Dtrfow(nth,np,nr,nsa) = 0.d0
            ! Fthfow(nth,np,nr,nsa) = 0.d0
            ! write(*,*)"Dttfow 2",Dttfow(nth,np,nr,nsa),orbit_th(nth,np,nr,nsa)%time(nstpmax),J_I

            deallocate(dIdu)

          end do
        end do
      end do
    end do

    ! calculate Drp, Drt, Drr, Frr
    do nsa = 1, nsamax
      do nr = 1, nrmax+1
        do np = 1, npmax
          do nth = 1, nthmax
            Drpfow(nth,np,nr,nsa) = 0.d0
            Drtfow(nth,np,nr,nsa) = 0.d0
            Drrfow(nth,np,nr,nsa) = 0.d0
            Frrfow(nth,np,nr,nsa) = 0.d0

            if ( nr == 1 ) then
              cycle    
            end if

            nstpmax = orbit_r(nth,np,nr,nsa)%nstp_max

            allocate(dIdu(3,3,nstpmax))

            mode = [0,0,1]
            call transformation_matrix(dIdu, orbit_r(nth,np,nr,nsa), nth, np, nr, nsa, mode)

            ! time-integral over an orbit then divide poloidal period 

            do nstp = 2, nstpmax
              cpitch_ob = cos(orbit_r(nth,np,nr,nsa)%theta(nstp))
              psip_ob   = orbit_r(nth,np,nr,nsa)%psip(nstp)
              thetap_ob = orbit_r(nth,np,nr,nsa)%thetap(nstp)

              ! calucurate local coefficient along orbit Dxx(p_orbit, theta_orbit, psip_orbit, thetap_orbit)
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

              ! Dxxfow = int_0^tau_p (integrand) dt
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
              if ( nth == nthmax/4 .and. np == npmax/2 .and. nr == nrmax/2 .and. nsa == 1 ) then
                write(55,'(5I4,6ES12.4)')nstp,nth,np,nr,nsa,Dpp_ob,Dpt_ob,Fpp_ob,Dtp_ob,Dtt_ob,Fth_ob
              end if
                      
            end do

            if ( nr == 1 ) then
              J_I = Jacobian_I(nth,np,nr,nsa)
            else if ( nr == nrmax+1 ) then
              J_I = Jacobian_I(nth,np,nr-1,nsa)
            else
              J_I = (Jacobian_I(nth,np,nr,nsa)+Jacobian_I(nth,np,nr-1,nsa))*0.5d0
            end if

            Drpfow(nth,np,nr,nsa) = Drpfow(nth,np,nr,nsa) / orbit_r(nth,np,nr,nsa)%time(nstpmax) * J_I
            Drtfow(nth,np,nr,nsa) = Drtfow(nth,np,nr,nsa) / orbit_r(nth,np,nr,nsa)%time(nstpmax) * J_I
            Drrfow(nth,np,nr,nsa) = Drrfow(nth,np,nr,nsa) / orbit_r(nth,np,nr,nsa)%time(nstpmax) * J_I
            Frrfow(nth,np,nr,nsa) = Frrfow(nth,np,nr,nsa) / orbit_r(nth,np,nr,nsa)%time(nstpmax) * J_I
            ! Drpfow(nth,np,nr,nsa) = 0.d0
            ! Drtfow(nth,np,nr,nsa) = 0.d0
            ! Drrfow(nth,np,nr,nsa) = 0.d0
            ! Frrfow(nth,np,nr,nsa) = 0.d0
            if ( nth == nthmax/4 .and. np == npmax/2 .and. nr == nrmax/2 .and. nsa == 1 ) then
              write(55,'(4ES12.4)')Drpfow(nth,np,nr,nsa)/J_I,Drtfow(nth,np,nr,nsa)/J_I,Drrfow(nth,np,nr,nsa)/J_I,Frrfow(nth,np,nr,nsa)/J_I
            end if

            deallocate(dIdu)

          end do
        end do
      end do
    end do
    
    open(90,file="./txt/fow_DPP.txt")
    open(91,file="./txt/fow_DTT.txt")
    open(92,file="./txt/fow_DPT.txt")
    open(93,file="./txt/fow_DTP.txt")
    open(94,file="./txt/fow_DRR.txt")

    open(95,file="./txt/fow_DPR.txt")
    open(96,file="./txt/fow_DTR.txt")
    open(97,file="./txt/fow_DRP.txt")
    open(98,file="./txt/fow_DRT.txt")

    open(100,file="./txt/fow_FPP.txt")
    open(101,file="./txt/fow_FTH.txt")
    open(102,file="./txt/fow_FRR.txt")
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            write(90,'(A,4I4,ES12.4)')"DPP",nth,np,nr,nsa,dppfow(nth,np,nr,nsa)/Jacobian_I(nth,np,nr,nsa)
            write(91,'(A,4I4,ES12.4)')"DTT",nth,np,nr,nsa,dttfow(nth,np,nr,nsa)/Jacobian_I(nth,np,nr,nsa)
            write(92,'(A,4I4,ES12.4)')"DPT",nth,np,nr,nsa,dptfow(nth,np,nr,nsa)/Jacobian_I(nth,np,nr,nsa)
            write(93,'(A,4I4,ES12.4)')"DTP",nth,np,nr,nsa,dtpfow(nth,np,nr,nsa)/Jacobian_I(nth,np,nr,nsa)
            write(94,'(A,4I4,ES12.4)')"DRR",nth,np,nr,nsa,drrfow(nth,np,nr,nsa)/Jacobian_I(nth,np,nr,nsa)

            write(95,'(A,4I4,ES12.4)')"DPR",nth,np,nr,nsa,dprfow(nth,np,nr,nsa)/Jacobian_I(nth,np,nr,nsa)
            write(96,'(A,4I4,ES12.4)')"DTR",nth,np,nr,nsa,dttfow(nth,np,nr,nsa)/Jacobian_I(nth,np,nr,nsa)
            write(97,'(A,4I4,ES12.4)')"DRP",nth,np,nr,nsa,drpfow(nth,np,nr,nsa)/Jacobian_I(nth,np,nr,nsa)
            write(98,'(A,4I4,ES12.4)')"DRT",nth,np,nr,nsa,drrfow(nth,np,nr,nsa)/Jacobian_I(nth,np,nr,nsa)

            write(100,'(A,4I4,ES12.4)')"FPP",nth,np,nr,nsa,fppfow(nth,np,nr,nsa)/Jacobian_I(nth,np,nr,nsa)
            write(101,'(A,4I4,ES12.4)')"FTH",nth,np,nr,nsa,fthfow(nth,np,nr,nsa)/Jacobian_I(nth,np,nr,nsa)
            write(102,'(A,4I4,ES12.4)')"FRR",nth,np,nr,nsa,frrfow(nth,np,nr,nsa)/Jacobian_I(nth,np,nr,nsa)
          end do
        end do
      end do
    end do
    close(90)
    close(91)
    close(92)
    close(93)
    close(94)
    close(95)
    close(96)
    close(97)
    close(98)
    close(100)
    close(101)
    close(102)

  end subroutine bounce_average

  subroutine transformation_matrix(dIdu, orbit_in, nth_in, np_in, nr_in, nsa_in, mode)

    use fpcomm
    use fowcomm
    use foworbit

    implicit none

    real(rkind),intent(out) :: dIdu(:,:,:)
    type(orbit),intent(in) :: orbit_in
    integer,intent(in) :: nth_in, np_in, nr_in, nsa_in, mode(3)

    ! elements of transformation matrix, dIdu
    real(rkind) ::  dpdp,  dxidp,  dpsdp,&
                    dpdth, dxidth, dpsdth,&
                    dpdr,  dxidr,  dpsdr
    real(rkind) :: A,B,C,D,E
    real(rkind) :: Fl, Bl, xil, pl, dFoBdpsi, Fob, Bob, thetaob, eps
    real(rkind) :: dFl, dBl, sinob, cosob
    integer :: nstp, nstpmax

    select case(mode(1))
    case(0)
      if ( mode(2) == 0 .and. mode(3) == 0 ) xil = cos(thetam(nth_in,np_in,nr_in,nsa_in))
      if ( mode(2) == 1 .and. mode(3) == 0 ) xil = cos(thetam_pg(nth_in,np_in,nr_in,nsa_in))
      if ( mode(2) == 0 .and. mode(3) == 1 ) xil = cos(thetam_rg(nth_in,np_in,nr_in,nsa_in))
    case(1)
      xil = cos(thetamg(nth_in,np_in,nr_in,nsa_in))
    end select

    select case(mode(2))
    case(0)
      pl = pm(np_in,nsa_in)*ptfp0(nsa_in)
    case(1)
      pl = pg(np_in,nsa_in)*ptfp0(nsa_in)
    end select

    select case(mode(3))
    case(0)
      Fl = Fpsi(nr_in)
      dFl = dFdpsi(nr_in)
      if ( xil * aefp(nsa_in) >= 0.d0 ) then
        Bl = Bout(nr_in)
        dBl = dBmdpsi(nr_in,1)
      else
        Bl = Bin(nr_in)
        dBl = dBmdpsi(nr_in,2)
      end if
    case(1)
      Fl = Fpsig(nr_in)
      dFl = dFgdpsi(nr_in)
      if ( xil * aefp(nsa_in) >= 0.d0 ) then
        Bl = Boutg(nr_in)
        dBl = dBmgdpsi(nr_in,1)
      else
        Bl = Bing(nr_in)
        dBl = dBmgdpsi(nr_in,2)
      end if
    end select

    nstpmax = orbit_in%nstp_max

    dFoBdpsi = (dFl*Bl-Fl*dBl)/Bl**2

    do nstp = 1, nstpmax
      Fob = get_F_nstp(orbit_in, nstp)
      Bob = orbit_in%Babs(nstp)
      sinob = SIN( orbit_in%theta(nstp) )
      cosob = COS( orbit_in%theta(nstp) )
    
      A = Fob/Bob*cosob-Fl/Bl*xil
      B = -(1.d0-xil**2)/(2.d0*xil)*dBl/Bl
      C = dFoBdpsi*pl*xil-AEFP(nsa_in)
      D = pl*sinob/Bob*(Fl*cosob/xil-Fob)
      E = Bl/Bob*sinob*cosob/xil
  
      dpdp =1.d0
      dpdth=0.d0
      dpdr =0.d0
      dxidr=0.d0
      dpsdr=0.d0
      dpsdp = A/(B*Fl/Bl*pl+C)    *ptfp0(nsa_in)/psi0
      dxidp = A/(B*Fl/Bl*pl+C)*B  *ptfp0(nsa_in)
      dpsdth= D/(B*Fl/Bl*pl+C)    /psi0
      dxidth= D/(B*Fl/Bl*pl+C)*B-E

      if ( abs(xil) >= 1.d0 ) then
        eps = 1.d-3/dble(nthmax) 
      else
        eps = 0.d0
      end if
  
      dIdu(1,1,nstp) = dpdp
      dIdu(1,2,nstp) = dxidp  / sqrt(1.d0-(xil**2-eps))*(-1.d0) ! convert dxi/dp  to dthetam/dp, avoid 1.d0-xil**2 = 0 by eps
      dIdu(1,3,nstp) = dpsdp
      dIdu(2,1,nstp) = dpdth
      dIdu(2,2,nstp) = dxidth / sqrt(1.d0-(xil**2-eps))*(-1.d0) ! convert dxi/dth to dthetam/dth
      dIdu(2,3,nstp) = dpsdth
      dIdu(3,1,nstp) = dpdr
      dIdu(3,2,nstp) = dxidr
      dIdu(3,3,nstp) = dpsdr 

      write(19,*)"dxidp  ",dIdu(1,2,nstp)
      write(19,*)"dpsdp  ",dIdu(1,3,nstp)
      write(19,*)"dthdpth",dIdu(2,2,nstp)
      write(19,*)"dpsdth ",dIdu(2,3,nstp)
      
    end do

  end subroutine transformation_matrix

  subroutine make_U_Dxy(U_Dxy, Dxyl, x, nsa)
    use fpcomm
    use fowcomm

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
    implicit none

    real(rkind),intent(out) :: C_out
    real(rkind),intent(in) :: U(:,:,:,:,:,:), check0, cpitch_in, psip_in, thetap_in
    integer :: nxmax, nymax, nzmax, ierr
    ierr = 0

    nxmax = size(U,4)
    nymax = size(U,5)
    nzmax = size(U,6)

    if ( check0 < 1.d-80 ) then
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