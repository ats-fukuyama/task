module fowcoef
  use fpcomm,only:rkind

  private

  public :: fow_coef

  real(rkind),allocatable,dimension(:,:,:,:,:) :: Dppl, Dptl, Fppl,&
                                                  Dtpl, Dttl, Fthl

  real(rkind),allocatable,dimension(:,:,:,:) :: Dppfow, Dptfow, Dprfow,&
                                                Dtpfow, Dttfow, Dtrfow,&
                                                Drpfow, Drtfow, Drrfow,&
                                                Fppfow,  Fthfow,  Frrfow
  
  real(rkind),allocatable :: FNSBL(:,:,:,:,:)

contains

  subroutine fow_coef

    use fpcomm
    use fowcomm
    use fowdistribution,only:convert_fI_to_fu
    use fpcalc
    use fpcoef
    use fpcalw

    implicit none

    integer :: nthp, nth, np, nr, nsa

    if ( .not.allocated(FNSBL) ) then
      allocate(FNSBL(nthmax,npmax,nrmax,nsamax,nthpmax))
      allocate(Dppl(nthmax,npmax+1,nrmax,nsamax,nthpmax),Dptl(nthmax,npmax+1,nrmax,nsamax,nthpmax))
      allocate(Dtpl(nthmax+1,npmax,nrmax,nsamax,nthpmax),Dttl(nthmax+1,npmax,nrmax,nsamax,nthpmax))
      allocate(Fppl(nthmax,npmax+1,nrmax,nsamax,nthpmax),Fthl(nthmax+1,npmax,nrmax,nsamax,nthpmax))

      allocate(Dppfow(nthmax,npmax+1,nrmax,nsamax),Dptfow(nthmax,npmax+1,nrmax,nsamax),Dprfow(nthmax,npmax+1,nrmax,nsamax))
      allocate(Dtpfow(nthmax+1,npmax,nrmax,nsamax),Dttfow(nthmax+1,npmax,nrmax,nsamax),Dtrfow(nthmax+1,npmax,nrmax,nsamax))
      allocate(Drpfow(nthmax,npmax,nrmax+1,nsamax),Drtfow(nthmax,npmax,nrmax+1,nsamax),Drrfow(nthmax,npmax,nrmax+1,nsamax))
      allocate(Fppfow(nthmax,npmax+1,nrmax,nsamax),Fthfow(nthmax,npmax+1,nrmax,nsamax),Frrfow(nthmax,npmax+1,nrmax,nsamax))
    end if


    do nthp = 1, nthpmax
      do nsa = 1, nsamax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              Dppl(nth,np,nr,nsa,nthp) = 0.d0
              Dptl(nth,np,nr,nsa,nthp) = 0.d0
              Dtpl(nth,np,nr,nsa,nthp) = 0.d0
              Dttl(nth,np,nr,nsa,nthp) = 0.d0
              Fppl(nth,np,nr,nsa,nthp) = 0.d0
              Fthl(nth,np,nr,nsa,nthp) = 0.d0
            end do
          end do
        end do
      end do  
    end do

    ! get local distribution function f(p,theta,r/a,thetap)
    call convert_fI_to_fu(FNSBL, FNSI)

    ! calculate local coefficient Dxxl(nth,np,nr,nsa,nthp)
    do nthp = 1, nthpmax

      ! calculate back graound distribution FNSB(nth,np,nr,nsa) for each nthp
      do nsa = 1, nsamax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              FNSB(nth,np,nr,nsa) = FNSBL(nth,np,nr,nsa,nthp)
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
        do nr = 1, nrmax
          do np = 1, npmax+1
            do nth = 1, nthmax+1

              if ( nth /= nthmax+1 ) then 
                Dppl(nth,np,nr,nsa,nthp) = Dppl(nth,np,nr,nsa,nthp) + DCPP(nth,np,nr,nsa) ! + DWPP(nth,np,nr,nsa)
                Dptl(nth,np,nr,nsa,nthp) = Dptl(nth,np,nr,nsa,nthp) + DCPT(nth,np,nr,nsa) ! + DWPT(nth,np,nr,nsa)    
                Fppl(nth,np,nr,nsa,nthp) = Fppl(nth,np,nr,nsa,nthp) + FCPP(nth,np,nr,nsa) ! + FEPP(nth,np,nr,nsa)
              end if
              if ( np /= npmax+1 ) then
                Dtpl(nth,np,nr,nsa,nthp) = Dtpl(nth,np,nr,nsa,nthp) + DCTP(nth,np,nr,nsa) ! + DWTP(nth,np,nr,nsa) 
                Dttl(nth,np,nr,nsa,nthp) = Dttl(nth,np,nr,nsa,nthp) + DCTT(nth,np,nr,nsa) ! + DWTT(nth,np,nr,nsa)
                Fthl(nth,np,nr,nsa,nthp) = Fthl(nth,np,nr,nsa,nthp) + FCTH(nth,np,nr,nsa) ! + FETH(nth,np,nr,nsa)
              end if

            end do
          end do
        end do
      end do

    end do

    ! execute bounce average
    call fow_bounce_average

  end subroutine fow_coef

  subroutine fow_bounce_average

    use fpcomm
    use fowcomm
    use fowprep,only:first_order_derivative

    implicit none
    integer :: nth, np, nr, nsa, nthp, mode(3), nstp, nstpmax, ierr = 0
    real(rkind),allocatable :: dBmdpsi(:,:), dFdpsi(:), dBmgdpsi(:,:), dFgdpsi(:)
    real(rkind),allocatable :: dIdu(:,:,:)
    real(rkind),allocatable:: U_Dpp(:,:,:,:,:,:,:,:),& ! spline coefficient of Dppl in (\theta, \psi_p, \theta_p) plane
                              U_Dpt(:,:,:,:,:,:,:,:),& ! spline coefficient of Dptl in (\theta, \psi_p, \theta_p) plane
                              U_Fpp(:,:,:,:,:,:,:,:),& ! spline coefficient of Fppl in (\theta, \psi_p, \theta_p) plane
                              U_Dtp(:,:,:,:,:,:,:,:),& ! spline coefficient of Dtpl in (\theta, \psi_p, \theta_p) plane
                              U_Dtt(:,:,:,:,:,:,:,:),& ! spline coefficient of Dttl in (\theta, \psi_p, \theta_p) plane
                              U_Fth(:,:,:,:,:,:,:,:),& ! spline coefficient of Fthl in (\theta, \psi_p, \theta_p) plane
                              theta_p(:), FX(:,:,:), FY(:,:,:), FZ(:,:,:), FXY(:,:,:), FYZ(:,:,:), FZX(:,:,:), FXYZ(:,:,:)
    real(rkind) :: Dpp_ob, Dpt_ob, Fpp_ob, Dtp_ob, Dtt_ob, Fth_ob, D_pls, D_mns, dt
    real(rkind) :: cpitch_ob, thetap_ob, psip_ob


    ! if xi > 0, then second dimension of dBmdpsi = 1 and Bm = Bout
    ! if xi > 0, then second dimension of dBmdpsi = 2 and Bm = Bin
    allocate(dBmdpsi(nrmax,2), dFdpsi(nrmax))
    allocate(dBmgdpsi(nrmax+1,2), dFgdpsi(nrmax+1))
    
    call first_order_derivative(dFdpsi, Fpsi, psim)
    call first_order_derivative(dBmdpsi(:,1), Bout, psim)
    call first_order_derivative(dBmdpsi(:,2), Bin, psim)
    call first_order_derivative(dFgdpsi, Fpsig, psimg)
    call first_order_derivative(dBmdpsi(:,1), Boutg, psimg)
    call first_order_derivative(dBmdpsi(:,2), Bing, psimg)

    ! set spline variable
    allocate(theta_p(nthpmax))
    allocate(FX(nthmax,nrmax,nthpmax), FY(nthmax,nrmax,nthpmax), FZ(nthmax,nrmax,nthpmax), FXY(nthmax,nrmax,nthpmax))
    allocate(FZX(nthmax,nrmax,nthpmax), FYZ(nthmax,nrmax,nthpmax), FXYZ(nthmax,nrmax,nthpmax))

    do nthp = 1, nthpmax
      theta_p(nthp) = (nthp-1)*2.d0*pi/nthpmax
    end do

    allocate(U_Dpp(4,4,4,nthmax,nrmax,nthpmax,npmax+1,nsamax))
    allocate(U_Dpt(4,4,4,nthmax,nrmax,nthpmax,npmax+1,nsamax))
    allocate(U_Fpp(4,4,4,nthmax,nrmax,nthpmax,npmax+1,nsamax))
    allocate(U_Dtp(4,4,4,nthmax+1,nrmax,nthpmax,npmax,nsamax))
    allocate(U_Dtt(4,4,4,nthmax+1,nrmax,nthpmax,npmax,nsamax))
    allocate(U_Fth(4,4,4,nthmax+1,nrmax,nthpmax,npmax,nsamax))

    do nsa = 1, nsamax
      do np = 1, npmax+1
        call SPL3D(cosm,psim,theta_p,Dppl(:,np,:,nsa,:),FX,FY,FZ,FXY,FYZ,FZX,FXYZ,&
                  U_Dpp(:,:,:,:,:,:,np,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,0,0,0,IERR)
        call SPL3D(cosm,psim,theta_p,Dptl(:,np,:,nsa,:),FX,FY,FZ,FXY,FYZ,FZX,FXYZ,&
                  U_Dpt(:,:,:,:,:,:,np,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,0,0,0,IERR)
        call SPL3D(cosm,psim,theta_p,Fppl(:,np,:,nsa,:),FX,FY,FZ,FXY,FYZ,FZX,FXYZ,&
                  U_Fpp(:,:,:,:,:,:,np,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,0,0,0,IERR)

        if ( np /= npmax+1 ) then
          call SPL3D(cosg,psim,theta_p,Dtpl(:,np,:,nsa,:),FX,FY,FZ,FXY,FYZ,FZX,FXYZ,&
                    U_Dtp(:,:,:,:,:,:,np,nsa),nthmax+1,nrmax,nthmax+1,nrmax,nthpmax,0,0,0,IERR)
          call SPL3D(cosg,psim,theta_p,Dttl(:,np,:,nsa,:),FX,FY,FZ,FXY,FYZ,FZX,FXYZ,&
                    U_Dtt(:,:,:,:,:,:,np,nsa),nthmax+1,nrmax,nthmax+1,nrmax,nthpmax,0,0,0,IERR)
          call SPL3D(cosg,psim,theta_p,Fthl(:,np,:,nsa,:),FX,FY,FZ,FXY,FYZ,FZX,FXYZ,&
                    U_Fth(:,:,:,:,:,:,np,nsa),nthmax+1,nrmax,nthmax+1,nrmax,nthpmax,0,0,0,IERR)
        end if

      end do
    end do

    ! calculate Dpp, Dpt, Dpr, Fp
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax+1
          do nth = 1, nthmax

            nstpmax = orbit_p(nth,np,nr,nsa)%nstp_max

            allocate(dIdu(3,3,nstpmax))

            if ( thetam(nth,np,nr,nsa) <= pi/2.d0 ) then
              call transformation_matrix(dIdu, orbit_p(nth,np,nr,nsa), nth, np, nr, nsa, [0,1,0], dBmdpsi(:,1), dFdpsi)
            else
              call transformation_matrix(dIdu, orbit_p(nth,np,nr,nsa), nth, np, nr, nsa, [0,1,0], dBmdpsi(:,2), dFdpsi)
            end if

            ! time-integral over an orbit then divide poloidal period 
            Dppfow(nth,np,nr,nsa) = 0.d0
            Dptfow(nth,np,nr,nsa) = 0.d0
            Dprfow(nth,np,nr,nsa) = 0.d0
            Fppfow(nth,np,nr,nsa) = 0.d0

            do nstp = 2, nstpmax
              cpitch_ob = cos(orbit_p(nth,np,nr,nsa)%theta(nstp))
              psip_ob   = orbit_p(nth,np,nr,nsa)%psip(nstp)
              thetap_ob = orbit_p(nth,np,nr,nsa)%thetap(nstp)

              ! calucurate local coefficient along orbit Dxx(p_orbit, theta_orbit, psip_orbit, thetap_orbit)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,Dpp_ob,cosm,psim,theta_p&
                          ,U_Dpp(:,:,:,:,:,:,np,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,Dpt_ob,cosm,psim,theta_p&
                          ,U_Dpt(:,:,:,:,:,:,np,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,Fpp_ob,cosm,psim,theta_p&
                          ,U_Fpp(:,:,:,:,:,:,np,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)

              dt = orbit_p(nth,np,nr,nsa)%time(nstp)-orbit_p(nth,np,nr,nsa)%time(nstp-1)

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

            Dppfow(nth,np,nr,nsa) = Dppfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax)
            Dptfow(nth,np,nr,nsa) = Dptfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax)
            Dprfow(nth,np,nr,nsa) = Dprfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax)
            Fppfow(nth,np,nr,nsa) = Fppfow(nth,np,nr,nsa) / orbit_p(nth,np,nr,nsa)%time(nstpmax)

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

            nstpmax = orbit_th(nth,np,nr,nsa)%nstp_max

            allocate(dIdu(3,3,nstpmax))

            if ( thetamg(nth,np,nr,nsa) <= pi/2.d0 ) then
              call transformation_matrix(dIdu, orbit_th(nth,np,nr,nsa), nth, np, nr, nsa, [1,0,0], dBmdpsi(:,1), dFdpsi)
            else
              call transformation_matrix(dIdu, orbit_th(nth,np,nr,nsa), nth, np, nr, nsa, [1,0,0], dBmdpsi(:,2), dFdpsi)
            end if

            ! time-integral over an orbit then divide poloidal period 
            Dtpfow(nth,np,nr,nsa) = 0.d0
            Dttfow(nth,np,nr,nsa) = 0.d0
            Dtrfow(nth,np,nr,nsa) = 0.d0
            Fthfow(nth,np,nr,nsa) = 0.d0

            do nstp = 2, nstpmax
              cpitch_ob = cos(orbit_p(nth,np,nr,nsa)%theta(nstp))
              psip_ob   = orbit_p(nth,np,nr,nsa)%psip(nstp)
              thetap_ob = orbit_p(nth,np,nr,nsa)%thetap(nstp)

              ! calucurate local coefficient along orbit Dxx(p_orbit, theta_orbit, psip_orbit, thetap_orbit)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,D_pls,cosm,psim,theta_p&
                          ,U_Dpp(:,:,:,:,:,:,np+1,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,D_mns,cosm,psim,theta_p&
                          ,U_Dpp(:,:,:,:,:,:,np,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              Dpp_ob = (D_pls+D_mns)/2.d0

              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,D_pls,cosm,psim,theta_p&
                          ,U_Dpt(:,:,:,:,:,:,np+1,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,D_mns,cosm,psim,theta_p&
                          ,U_Dpt(:,:,:,:,:,:,np,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              Dpt_ob = (D_pls+D_mns)/2.d0

              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,D_pls,cosm,psim,theta_p&
                          ,U_Fpp(:,:,:,:,:,:,np+1,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,D_mns,cosm,psim,theta_p&
                          ,U_Fpp(:,:,:,:,:,:,np,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              Fpp_ob = (D_pls+D_mns)/2.d0

              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,Dtp_ob,cosg,psim,theta_p&
                          ,U_Dtp(:,:,:,:,:,:,np,nsa),nthmax+1,nrmax,nthmax+1,nrmax,nthpmax,IERR)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,Dtt_ob,cosg,psim,theta_p&
                          ,U_Dtt(:,:,:,:,:,:,np,nsa),nthmax+1,nrmax,nthmax+1,nrmax,nthpmax,IERR)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,Fth_ob,cosg,psim,theta_p&
                          ,U_Fth(:,:,:,:,:,:,np,nsa),nthmax+1,nrmax,nthmax+1,nrmax,nthpmax,IERR)


              dt = orbit_th(nth,np,nr,nsa)%time(nstp)-orbit_th(nth,np,nr,nsa)%time(nstp-1)

              ! Dxxfow = int_0^tau_p (integrand) dt
              Dtpfow(nth,np,nr,nsa) = Dtpfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,2,nstp) + Dpt_ob * dIdu(2,2,nstp) ) * dt

              Dttfow(nth,np,nr,nsa) = Dttfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,2,nstp)**2 + Dpt_ob * dIdu(1,2,nstp)*dIdu(2,2,nstp) &
                                    + Dtp_ob * dIdu(1,2,nstp)*dIdu(2,2,nstp) + Dtt_ob * dIdu(2,2,nstp)**2) * dt

              Dtrfow(nth,np,nr,nsa) = Dtrfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,2,nstp)*dIdu(1,3,nstp) + Dpt_ob * dIdu(1,2,nstp)*dIdu(2,3,nstp) &
                                    + Dtp_ob * dIdu(2,2,nstp)*dIdu(1,3,nstp) + Dtt_ob * dIdu(2,2,nstp)*dIdu(2,3,nstp)) * dt

              Fthfow(nth,np,nr,nsa) = Fthfow(nth,np,nr,nsa)&
                                    + ( Fpp_ob*dIdu(1,2,nstp) + Fth_ob*dIdu(2,2,nstp) ) * dt

            end do

            Dtpfow(nth,np,nr,nsa) = Dtpfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax)
            Dttfow(nth,np,nr,nsa) = Dttfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax)
            Dtrfow(nth,np,nr,nsa) = Dtrfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax)            
            Fthfow(nth,np,nr,nsa) = Fthfow(nth,np,nr,nsa) / orbit_th(nth,np,nr,nsa)%time(nstpmax)

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

            nstpmax = orbit_r(nth,np,nr,nsa)%nstp_max

            allocate(dIdu(3,3,nstpmax))

            if ( thetam(nth,np,nr,nsa) <= pi/2.d0 ) then
              call transformation_matrix(dIdu, orbit_r(nth,np,nr,nsa), nth, np, nr, nsa, [0,0,1], dBmdpsi(:,1), dFdpsi)
            else
              call transformation_matrix(dIdu, orbit_r(nth,np,nr,nsa), nth, np, nr, nsa, [0,0,1], dBmdpsi(:,2), dFdpsi)
            end if

            ! time-integral over an orbit then divide poloidal period 
            Drpfow(nth,np,nr,nsa) = 0.d0
            Drtfow(nth,np,nr,nsa) = 0.d0
            Drrfow(nth,np,nr,nsa) = 0.d0
            Frrfow(nth,np,nr,nsa) = 0.d0

            do nstp = 2, nstpmax
              cpitch_ob = cos(orbit_r(nth,np,nr,nsa)%theta(nstp))
              psip_ob   = orbit_r(nth,np,nr,nsa)%psip(nstp)
              thetap_ob = orbit_r(nth,np,nr,nsa)%thetap(nstp)

              ! calucurate local coefficient along orbit Dxx(p_orbit, theta_orbit, psip_orbit, thetap_orbit)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,D_pls,cosm,psim,theta_p&
                          ,U_Dpp(:,:,:,:,:,:,np+1,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,D_mns,cosm,psim,theta_p&
                          ,U_Dpp(:,:,:,:,:,:,np,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              Dpp_ob = (D_pls+D_mns)/2.d0

              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,D_pls,cosm,psim,theta_p&
                          ,U_Dpt(:,:,:,:,:,:,np+1,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,D_mns,cosm,psim,theta_p&
                          ,U_Dpt(:,:,:,:,:,:,np,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              Dpt_ob = (D_pls+D_mns)/2.d0

              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,D_pls,cosm,psim,theta_p&
                          ,U_Fpp(:,:,:,:,:,:,np+1,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,D_mns,cosm,psim,theta_p&
                          ,U_Fpp(:,:,:,:,:,:,np,nsa),nthmax,nrmax,nthmax,nrmax,nthpmax,IERR)
              Fpp_ob = (D_pls+D_mns)/2.d0

              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,Dtp_ob,cosg,psim,theta_p&
                          ,U_Dtp(:,:,:,:,:,:,np,nsa),nthmax+1,nrmax,nthmax+1,nrmax,nthpmax,IERR)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,Dtt_ob,cosg,psim,theta_p&
                          ,U_Dtt(:,:,:,:,:,:,np,nsa),nthmax+1,nrmax,nthmax+1,nrmax,nthpmax,IERR)
              call SPL3DF(cpitch_ob,psip_ob,thetap_ob,Fth_ob,cosg,psim,theta_p&
                          ,U_Fth(:,:,:,:,:,:,np,nsa),nthmax+1,nrmax,nthmax+1,nrmax,nthpmax,IERR)


              dt = orbit_r(nth,np,nr,nsa)%time(nstp)-orbit_r(nth,np,nr,nsa)%time(nstp-1)

              ! Dxxfow = int_0^tau_p (integrand) dt
              Drpfow(nth,np,nr,nsa) = Drpfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,3,nstp) + Dpt_ob * dIdu(2,3,nstp) ) * dt

              Drtfow(nth,np,nr,nsa) = Drrfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,2,nstp)*dIdu(1,3,nstp) + Dpt_ob * dIdu(1,3,nstp)*dIdu(2,2,nstp) &
                                    + Dtp_ob * dIdu(1,2,nstp)*dIdu(2,3,nstp) + Dtt_ob * dIdu(2,2,nstp)*dIdu(2,3,nstp)) * dt

              Drtfow(nth,np,nr,nsa) = Drtfow(nth,np,nr,nsa)&
                                    + ( Dpp_ob * dIdu(1,3,nstp)**2 + Dpt_ob * dIdu(1,3,nstp)*dIdu(2,3,nstp) &
                                    + Dtp_ob * dIdu(2,3,nstp)*dIdu(1,3,nstp) + Dtt_ob * dIdu(2,3,nstp)**2) * dt

              Frrfow(nth,np,nr,nsa) = Frrfow(nth,np,nr,nsa)&
                                    + ( Fpp_ob*dIdu(1,3,nstp) + Fth_ob*dIdu(2,3,nstp) ) * dt

            end do

            Drpfow(nth,np,nr,nsa) = Drpfow(nth,np,nr,nsa) / orbit_r(nth,np,nr,nsa)%time(nstpmax)
            Drtfow(nth,np,nr,nsa) = Drtfow(nth,np,nr,nsa) / orbit_r(nth,np,nr,nsa)%time(nstpmax)
            Drrfow(nth,np,nr,nsa) = Drrfow(nth,np,nr,nsa) / orbit_r(nth,np,nr,nsa)%time(nstpmax)            
            Frrfow(nth,np,nr,nsa) = Frrfow(nth,np,nr,nsa) / orbit_r(nth,np,nr,nsa)%time(nstpmax)

            deallocate(dIdu)

          end do
        end do
      end do
    end do

  end subroutine fow_bounce_average

  subroutine transformation_matrix(dIdu, orbit_in, nth_in, np_in, nr_in, nsa_in, mode, dBmdpsi, dFdpsi)

    use fpcomm
    use fowcomm
    use foworbit,only:func_orbit_F

    implicit none

    real(rkind),intent(out) :: dIdu(:,:,:)
    type(orbit),intent(in) :: orbit_in
    integer,intent(in) :: nth_in, np_in, nr_in, nsa_in, mode(3)
    real(rkind),intent(in) :: dBmdpsi(:), dFdpsi(:)

    ! elements of transformation matrix, dIdu
    real(rkind) ::  dpdp,  dxidp,  dpsdp,&
                dpdth, dxidth, dpsdth,&
                dpdr,  dxidr,  dpsdr
    real(rkind) :: A,B,C,D,E
    real(rkind) :: Fl, Bl, xil, pl, dFoBdpsi, Fob, Bob, thetaob
    integer :: nstp, nstpmax

    select case(mode(1))
    case(0)
      xil = cos(thetam(nth_in,np_in,nr_in,nsa_in))
    case(1)
      xil = cos(thetamg(nth_in,np_in,nr_in,nsa_in))
    end select

    select case(mode(2))
    case(0)
      pl = pm(np_in,nsa_in)
    case(1)
      pl = pg(np_in,nsa_in)
    end select

    select case(mode(3))
    case(0)
      Fl = Fpsi(nr_in)
      if ( xil >= 0.d0 ) then
        Bl = Bout(nr_in)
      else
        Bl = Bin(nr_in)
      end if
    case(1)
      Fl = Fpsig(nr_in)
      if ( xil >= 0.d0 ) then
        Bl = Boutg(nr_in)
      else
        Bl = Bing(nr_in)
      end if
    end select

    nstpmax = orbit_in%nstp_max

    do nstp = 1, nstpmax
      Fob = func_orbit_F(orbit_in, nstp, nr_in)
      Bob = orbit_in%Babs(nstp)
      thetaob = orbit_in%theta(nstp)
  
      dFoBdpsi = (dFdpsi(nr_in)*Bl-Fl*dBmdpsi(nr_in))/Bl**2
  
      A =Fob/Bob*cos(thetaob)-Fl/Bl*xil
      B = -(1.d0-xil**2)/(2.d0*xil)*dBmdpsi(nr_in)/Bl
      C = dFoBdpsi*pl*xil-AEFP(nsa_in)
      D = pl*sin(thetaob)/Bob*(Fl*cos(thetaob)/xil-Fob)
      E = Bl/Bob*sin(thetaob)*cos(thetaob)/xil
  
      dpdp =1.d0
      dpdth=0.d0
      dpdr =0.d0
      dxidr=0.d0
      dpsdr=0.d0
      dpsdp = A/(B*dFoBdpsi*pl+C)
      dxidp = A/(B*dFoBdpsi*pl+C)*B
      dpsdth= D/(B*dFoBdpsi*pl+C)
      dxidth= D/(B*dFoBdpsi*pl+C)*B-E
  
      dIdu(1,1,nstp) = dpdp
      dIdu(1,2,nstp) = dxidp  / sqrt(1.d0-xil**2)*(-1.d0) ! convert dxi/dp  to dthetam/dp
      dIdu(1,3,nstp) = dpsdp
      dIdu(2,1,nstp) = dpdth
      dIdu(2,2,nstp) = dxidth / sqrt(1.d0-xil**2)*(-1.d0) ! convert dxi/dth to dthetam/dth
      dIdu(2,3,nstp) = dpsdth
      dIdu(3,1,nstp) = dpdr
      dIdu(3,2,nstp) = dxidr
      dIdu(3,3,nstp) = dpsdr  
      
    end do

  end subroutine transformation_matrix

end module fowcoef