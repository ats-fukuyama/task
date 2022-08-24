! fowevalNC.f90

MODULE fowevalNC

  PRIVATE

  PUBLIC fow_fout_NC
  PUBLIC eval_analytic
  PUBLIC eval_rw
  PUBLIC eval_fow

CONTAINS

  SUBROUTINE fow_gout_NC

    USE fpcomm
    USE fowcomm
    IMPLICIT NONE

    REAL(rkind),DIMENSION(nrmax,nsamax):: &
         Dn_banana,Dn_plateau,chi_banana,chi_plateau, &
         Dn_rw,Dn_rw_bav,chi_rw,chi_rw_bav, &
         Dn_fow,chi_fow
    INTEGER:: nth, np, nr, nsa

    CALL eval_analytic(Dn_banana,Dn_plateau,chi_banana,chi_plateau)
  
    CALL eval_rw(Dn_rw,Dn_rw_bav,chi_rw,chi_rw_bav)

    CALL eval_fow(Dn_fow,chi_fow)

    CALL PAGES
    CALL grd1d(1,rm,Dn_banana(:,1),nrmax,nrmax,1,'@Dn_banana_e vs r@')
    CALL grd1d(2,rm,Dn_banana(:,2),nrmax,nrmax,1,'@Dn_banana_i vs r@')
    CALL grd1d(3,rm,Dn_plateau(:,1),nrmax,nrmax,1,'@Dn_plateau_e vs r@')
    CALL grd1d(4,rm,Dn_plateau(:,2),nrmax,nrmax,1,'@Dn_plateau_i vs r@')
    CALL PAGEE
  
    CALL PAGES
    CALL grd1d(1,rm,Dn_rw(:,1),nrmax,nrmax,1,'@Dn_rw_e vs r@')
    CALL grd1d(2,rm,Dn_rw(:,2),nrmax,nrmax,1,'@Dn_rw_i vs r@')
    CALL grd1d(3,rm,Dn_rw_bav(:,1),nrmax,nrmax,1,'@Dn_rw_bav_e vs r@')
    CALL grd1d(4,rm,Dn_rw_bav(:,2),nrmax,nrmax,1,'@Dn_rw_bav_i vs r@')
    CALL PAGEE
  
    CALL PAGES
    CALL grd1d(1,rm,Dn_fow(:,1),nrmax,nrmax,1,'@Dn_fow_e vs r@')
    CALL grd1d(2,rm,Dn_fow(:,2),nrmax,nrmax,1,'@Dn_fow_i vs r@')
    CALL PAGEE

  END SUBROUTINE fow_gout_NC
  
  SUBROUTINE fow_fout_NC
  !----------------------------------------
  ! module for file out put as *.txt
  !----------------------------------------
    USE fpcomm
    USE fowcomm
    USE fpwrite

    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax):: &
         Dn_banana,Dn_plateau,chi_banana,chi_plateau, &
         Dn_rw,Dn_rw_bav,chi_rw,chi_rw_bav, &
         Dn_fow,chi_fow

    !**** make data folder

    CALL system('mkdir -p dat')

    !**** calculation of physical values

    CALL eval_analytic(Dn_banana,Dn_plateau,chi_banana,chi_plateau)

    CALL eval_rw(Dn_rw,Dn_rw_bav,chi_rw,chi_rw_bav)

    CALL eval_fow(Dn_fow,chi_fow)

    !**** write txt file
    CALL fptxt2D(Dn_banana,"dat/Dn_banana.txt")
    CALL fptxt2D(Dn_plateau,"dat/Dn_plateau.txt")
    CALL fptxt2D(chi_banana,"dat/chi_banana.txt")
    CALL fptxt2D(chi_plateau,"dat/chi_plateau.txt")

    CALL fptxt2D(Dn_rw,"dat/Dn_rw.txt")
    CALL fptxt2D(Dn_rw_bav,"dat/Dn_rw_bav.txt")
    CALL fptxt2D(chi_rw,"dat/chi_rw.txt")
    CALL fptxt2D(chi_rw_bav,"dat/chi_rw_bav.txt")

    CALL fptxt2D(Dn_fow,"dat/Dn_fow.txt")
    CALL fptxt2D(chi_fow,"dat/chi_fow.txt")

  END SUBROUTINE fow_fout_NC
  

  SUBROUTINE eval_analytic(Dn_banana,Dn_plateau,chi_banana,chi_plateau)

    !------------------------------------------------------------
    ! SUBROUTINE for neoclassical particle and heat diffusion coefficients
    ! in banana and plateau regime
    ! Formulae from "Tokamaks fourth edition"
    !------------------------------------------------------------

    USE fpcomm
    USE fowcomm
    USE fowlib

    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: &
         Dn_banana,Dn_plateau,chi_banana,chi_plateau
    REAL(rkind),DIMENSION(nrmax,nsamax)  :: &
         pn_nr_nsa,pt_nr_nsa,dpnpr_nr_nsa,dptdr_nr_nsa
    REAL(rkind) fact,fact_ba,fact_pla,tau_ele,rho_e,eps_t,Baxis,B_p
    INTEGER nr, nsa

    fact_taue = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bing(1) ! approximation on B by Baxis

    ! --- radial postion and particle species ---
    
    DO nsa=1,nsamax
       DO nr = 1, nrmax
          pn_nr_nsa(nr,nsa)=rnsl(nr,nsa)*1.D20   ! particle density in m^{-3}
          pt_nr_nsa(nr,nsa)=rwsl(nr,nsa)*1.D6/(1.5D0*pn_nr_nsa(nr_nsa)) !T in J
       END DO
    END DO
       
    ! --- radial derivative of density and temperature ---
    DO nsa = 1, nsamax
      CALL first_order_derivative(dpndr_nr_nsa(:,nsa), pn_nr_nsa(:,nsa), rm)
      CALL first_order_derivative(dptdr_nr_nsa(:,nsa), pt_nr_nsa(:,nsa), rm)
    END DO

    ! --- calculate diffusion coefficient ---

    DO nr = 1, nrmax
       tau_ele = fact*sqrt(AMFP(1))*(SQRT(pt_nr_nsa(nr,1))**3) &
             (rnsl(nr,2)*1.d20*AEFP(2)**2*lnlam(nr,2,1)*AEE**2)
        rho_e = sqrt(2*Ta(nr,1)/AMFP(1))*AMFP(1)/(AEE*Baxis)
        eps_t = rm(nr)*RA/RR
        B_p = rm(nr)*RA*BB/(safety_factor(nr)*RR)
        fact_ba = (safety_factor(nr)**2)*1.d20*rnsl(nr,nsa) &
             *(rho_e**2)/((eps_t**1.5)*tau_ele)*RKEV ! modified by anzai[*RKEV]
        fact_pla = -sqrt(pi)*(eps_t**2)*Ta(nr,1)*RKEV &! modified by anzai [*RKEV]
             *rho_e*1.d20*rnsl(nr,nsa)/(2*AEE*B_p*rm(nr)*RA)

    DO nsa = 1, nsamax
        Dn_banana(nr,nsa) = 0.D0
        Dn_plateau(nr,nsa) = 0.D0
        chi_banana(nr,nsa) = 0.D0
        chi_plateau(nr,nsa) = 0.D0
     END DO
   END DO

 END SUBROUTINE eval_analytic

 SUBROUTINE eval_rw(Dn_rw,Dn_rw_bav,chi_rw,chi_rw_bav)
    USE fpcomm
    USE fowcomm
    USE fowlib

    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: &
         Dn_rw,Dn_rw_bavplateau,chi_banana,chi_plateau
    REAL(rkind),DIMENSION(nrmax,nsamax)  :: Ta,dTadr,dndr
    REAL(rkind) fact,fact_ba,fact_pla,tau_ele,rho_e,eps_t,Baxis,B_p
    INTEGER nr, nsa



    SUBROUTINE Heatneoclass(heatba,heatpla,chi_neo_ba,chi_neo_pla)
  !---------------------------------------------------------------
  !SUBROUTINE for neoclassical heat flux in banana, plateau region
  !calculate heat flux and diffusion coefficient
  ! main calcul using[J]
  ! Formulae from "Tokamaks fourth edition"
  !--------------------------------------------------------------

    USE fpcomm
    USE fowcomm
    USE fowlib

    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: heatba, heatpla
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: chi_neo_ba, chi_neo_pla
    REAL(rkind),DIMENSION(nrmax,nsamax) :: Sr_ba, Sr_pla
    REAL(rkind),DIMENSION(nrmax,nsamax) :: Ta, dTadr, dndr
    REAL(rkind), DIMENSION(nsamax) :: rho_a
    REAL(rkind) fact, fact_s, tau_ele,tau_i, eps_t, Baxis, B_p
    INTEGER nr, nsa

    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bing(1) ! approximation on B by Baxis

    !****temperature make
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)
        !Ta(temperature)[J]
      END DO
    END DO

    !****first order derivation
    DO nsa = 1, nsamax
      CALL first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
      CALL first_order_derivative(dndr(:,nsa), rnsl(:, nsa), rm)
    END DO

    !****calculate flux and diffusion coefficient
    DO nsa = 1, nsamax
      DO nr = 1, nrmax

        tau_ele = fact*sqrt(AMFP(1))*((Ta(nr,1))**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**2*AEE**2*lnlam(nr,2,1))
        tau_i   = fact*sqrt(2.d0)*sqrt(AMFP(2))*((Ta(nr,2))**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2))
        !**** *RKEV converts [keV] to [J]

        rho_a(nsa) = sqrt(Ta(nr,nsa)/AMFP(nsa))*AMFP(nsa)/(AEFP(nsa)*Baxis)
        !**** *RKEV converts [J] to [keV]
        eps_t = rm(nr)*RA/RR
        B_p = rm(nr)*RA*BB/(safety_factor(nr)*RR)
        !****calculate heat flux
        IF (nsa == 1) THEN
         fact_s = (safety_factor(nr)**2) &
             *(rho_a(nsa)**2)/((eps_t**1.5)*tau_ele)
         Sr_ba(nr,nsa) = fact_s*rnsl(nr,nsa)*1.D20*Ta(nr,1) &
            * ((1.53d0*(1+Ta(nr,2)/Ta(nr,1))*dndr(nr, nsa))/rnsl(nr,nsa) &
              - 1.81d0*dTadr(nr,1)/Ta(nr,1) &
              - 0.27d0*dTadr(nr,2)/Ta(nr,1))/RKEV!test
        ! neglect E_para term

         Sr_pla(nr,nsa) = - sqrt(PI)/4*eps_t**2*Ta(nr,nsa) &
               / (AEE*B_p)*rho_a(nsa)*rnsl(nr,nsa)*1.d20*Ta(nr,nsa) &
               / rm(nr)*( (1+Ta(nr,2)/Ta(nr,1))*dndr(nr,nsa) &
               / (rnsl(nr,nsa)*1.d20) &
               + 7.5d0*dTadr(nr,nsa)/Ta(nr,nsa) &
               + 1.5d0*dTadr(nr,2)/Ta(nr,nsa))/RKEV!test
        ! neglect E_para term
        ELSE
         fact_s = (safety_factor(nr)**2) &
              *(rho_a(nsa)**2)/((eps_t**1.5)*tau_i)
         Sr_ba(nr,nsa) = - 0.68d0*fact_s*(1 + 0.48d0*sqrt(eps_t)) &
              *rnsl(nr,nsa)*1.d20*dTadr(nr,2)/RKEV !test
         Sr_pla(nr,nsa) = - 1.5d0*sqrt(pi) &
                        * eps_t**2*Ta(nr,nsa)/(AEE*B_p) &
                        * rho_a(nsa)/rm(nr)*rnsl(nr,nsa) &
                        * 1.d20*dTadr(nr,nsa)/RKEV !test
        END IF

        heatba(nr,nsa)  = Sr_ba(nr, nsa)
        heatpla(nr,nsa) = Sr_pla(nr, nsa)
        chi_neo_ba      = -Sr_ba(nr,nsa)/(dTadr(nr,nsa)/RKEV)
        chi_neo_pla     = -Sr_pla(nr,nsa)/(dTadr(nr,nsa)/RKEV)
      END DO
    END DO

  END SUBROUTINE Heatneoclass


  SUBROUTINE bounce_average_for_Drw(Dout, Din)

    USE fpcomm
    USE fowcomm
    USE fowcoef

    IMPLICIT NONE
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax),INTENT(OUT) :: DOut
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax),INTENT(IN)  :: Din
    REAL(rkind),DIMENSION(3,3,max_stp) :: dIdul
    REAL(rkind),ALLOCATABLE :: U(:,:,:,:,:,:,:,:), Drwl(:,:,:,:,:)
    REAL(rkind) dt, taup, cpitch_ob, psip_ob, thetap_ob, Drw_ob
    type(orbit) ob
    INTEGER nth, np, nr, nsa, nstp, nstpmax, nthp, mode(3)

    allocate(U(4,4,4,nthmax,nrmax,nthpmax,npmax,nsamax))
    allocate(Drwl(nthmax,npmax,nrmax,nthpmax,nsamax))

    mode = [0,0,0]

    DO nsa = 1, nsamax
      DO nthp = 1, nthpmax
        DO nr = 1, nrmax
          DO np = 1, npmax
            DO nth = 1, nthmax
              Drwl(nth,np,nr,nthp,nsa) = Din(nth,np,nr,nsa)
            END DO
          END DO
        END DO
      END DO
    END DO

    DO nsa = 1, nsamax
      CALL make_U_Dxy(U, Drwl, 'm', nsa)
    END DO

    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        DO np = 1, npmax
          DO nth = 1, nthmax
            ob = orbit_m(nth,np,nr,nsa)
            nstpmax = ob%nstp_max
            taup = ob%time(nstpmax)

            CALL transformation_matrix(dIdul, ob, nth, np, nr, nsa, mode)

            DOut(nth,np,nr,nsa) = 0.d0

            DO nstp = 2, nstpmax
              dt = ob%time(nstp)-ob%time(nstp-1)
              cpitch_ob = ob%costh(nstp)
              psip_ob   = ob%psip(nstp)
              thetap_ob = ob%thetap(nstp)
              CALL interpolate_D_unlessZero(Drw_ob, U(:,:,:,:,:,:,np,nsa), &
                    1.d0, cpitch_ob, psip_ob, thetap_ob)

              DOut(nth,np,nr,nsa) = DOut(nth,np,nr,nsa)&
                                  + Drw_ob*dIdul(3,3,nstp)**2*dt
            END DO
            DOut(nth,np,nr,nsa) = DOut(nth,np,nr,nsa)/taup
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE bounce_average_for_Drw

  SUBROUTINE D_random_walk_baverage(Drwav, Drw, Drweff, Drweffav)
  !-----------------------------------------------------
  ! For particle diffusion in banana region
  ! modified and added by anzai [2022/3/2]
  !-----------------------------------------------------
    USE fpcomm
    USE fowcomm
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: Drw, Drwav
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: Drweff, Drweffav
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: dfdrhom
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: fnsp_l
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: Drweff_l, Drweff_ba
    REAL(rkind),DIMENSION(npmax,nrmax,nsamax) :: nu
    REAL(rkind),DIMENSION(nrmax,nsamax) :: Ta
    REAL(rkind),DIMENSION(nrmax,nsamax) :: nu_eff
    REAL(rkind) step_len, step_time, step_time_eff, eps_t
    REAL(rkind) rho_theta, Baxis, pv, B_theta, rho
    REAL(rkind) dVI, sumVI
    INTEGER nth, np, nr, nsa, ns

!     !**** initialization ****
!     Baxis = Bing(1)
!     Drw(:,:) = 0.d0
!     Drwav(:,:) = 0.d0
!     sumVI = 0.d0
!     CALL nu_deflection(nu)
!
! !********** for new dfdrhom module [2022/3/4]
!     fnsp_l(:,:,:,:)=0.d0
!     DO nsa = 1, nsamax
!       ns = ns_nsa(nsa)
!       DO nr= 1, nrmax
!         DO np = 1, npmax
!           IF ( pm(np,ns) > fact_bulk ) exit
!           DO nth = 1, nthmax
!           fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa)
!           END DO
!         END DO
!       END DO
!     END DO
!
!     !**** first order derivative
!     DO nsa = 1, nsamax
!       DO np = 1, npmax
!         DO nth = 1, nthmax
!           CALL first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp_l(nth,np,:,nsa), rm)
!         END DO
!       END DO
!     END DO
! !!**** END of new dfdrhom module
!
!     !************************by anzai[2022/2/15]
!     DO nsa = 1, nsamax
!       DO nr = 1, nrmax
!         Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
!         !** Ta(temperature)[keV]
!       END DO
!     END DO
!
!     CALL nu_effective(nu_eff)
!
!     DO nsa = 1, nsamax
!       DO nr = 1, nrmax
!         eps_t = rm(nr)*RA/RR
!         DO np = 1, npmax
!
!           step_time = 1.d0/nu(np, nr, nsa)
!           step_time_eff = 1.d0/nu_eff(nr, nsa)
!           rho = pm(np,nsa)*ptfp0(nsa)/(aefp(nsa)*Baxis)
!           !** [keV]verocity's length
!
!           DO nth = 1, nthmax
!             step_len = rho * safety_factor(nr)/sqrt(eps_t)
!             Drwlocal(nth,np,nr,nsa) = step_len**2/step_time*sqrt(eps_t)
!             Drweff_l(nth,np,nr,nsa) = step_len**2/step_time_eff*sqrt(eps_t)
!           END DO
!         END DO
!       END DO
!     END DO
!
!     CALL bounce_average_for_Drw(Drwba, Drwlocal)
!     CALL bounce_average_for_Drw(Drweff_ba, Drweff_l)
!
!     DO nsa = 1, nsamax
!       DO nr = 1, nrmax
!         DO np = 1, npmax
!           DO nth = 1, nthmax
!             dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JI(nth,np,nr,nsa)
!
!             Drw(nr,nsa) = Drw(nr,nsa) - Drwlocal(nth,np,nr,nsa) &
!                         * dfdrhom(nth,np,nr,nsa)*dVI
!             Drwav(nr,nsa) = Drwav(nr,nsa) - Drwba(nth,np,nr,nsa) &
!                           * dfdrhom(nth,np,nr,nsa)*dVI
!             Drweff(nr,nsa) = Drweff(nr,nsa) - Drweff_l(nth,np,nr,nsa) &
!                            * dfdrhom(nth,np,nr,nsa)*dVI
!             Drweffav(nr,nsa) = Drweffav(nr,nsa) - Drweff_ba(nth,np,nr,nsa) &
!                              * dfdrhom(nth,np,nr,nsa)*dVI
!             sumVI = sumVI - dfdrhom(nth,np,nr,nsa)*dVI
!           END DO
!         END DO
!         Drw(nr,nsa) = Drw(nr,nsa)/sumVI
!         Drwav(nr,nsa) = Drwav(nr,nsa)/sumVI
!         Drweff(nr,nsa) = Drweff(nr,nsa)/sumVI
!         Drweffav(nr,nsa) = Drweffav(nr,nsa)/sumVI
!       END DO
!     END DO

  END SUBROUTINE D_random_walk_baverage

  SUBROUTINE Heat_rw_baverage(heatrwav, heatrw)
  !---------------------------------------------
  ! heat diffusion coef calcul module
  ! This module corresponds to integral_Heatdiff
  !---------------------------------------------

    USE fpcomm
    USE fowcomm
    USE fowlib

    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT) :: heatrw, heatrwav
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: dfdrhom, fnsp_l
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
    REAL(rkind),DIMENSION(nrmax,nsamax) :: nu, Ta, Drw, vth
    REAL(rkind) step_len, step_time, eps_t, rho_theta, Baxis, pv, rho
    REAL(rkind) dVI, sumVI
    INTEGER nth, np, nr, nsa, ns

    !**** initialization ****
    Baxis = Bing(1)
    heatrw(:,:) = 0.d0
    heatrwav(:,:) = 0.d0
    sumVI = 0.d0
 !**** for new dfdrhom module [2022/3/4]
     fnsp_l(:,:,:,:)=0.d0
     DO nsa = 1, nsamax
       ns = ns_nsa(nsa)
       DO nr= 1, nrmax
         DO np = 1, npmax
           IF ( pm(np,ns) > fact_bulk ) exit
           DO nth = 1, nthmax
           fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa)
           END DO
         END DO
       END DO
     END DO

     !**** first order derivative
     DO nsa = 1, nsamax
       DO np = 1, npmax
         DO nth = 1, nthmax
           CALL first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp_l(nth,np,:    ,nsa), rm)
         END DO
       END DO
     END DO
 !**** END of new dfdrhom module

    !**** Temperature and flecency make
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)
        !**** Ta(temperature)[J]
        ! /RKEV = 1/(AEE*1.D3) unit converter [keV]to [J]
      END DO
    END DO

   CALL nu_effective(nu)

    !**** Calculation of particle diff coef
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        eps_t = rm(nr)*RA/RR
        DO np = 1, npmax
          step_time = 1.d0/nu(nr,nsa)
          !** thermal rho
!          vth(nr,nsa) = SQRT( 2.d0*Ta(nr,nsa)/AMFP(nsa)*RKEV)
!          rho = vth(nr,nsa)*AMFP(nsa)/(aefp(nsa)*Baxis)
          !** original rho
          rho = (pm(np,nsa)*ptfp0(nsa) ) &!/sqrt(AEE*1.d3)) &
              / (aefp(nsa)*Baxis)
          !****Unit [J]

          DO nth = 1, nthmax
            step_len = rho * safety_factor(nr)/sqrt(eps_t)
            Drwlocal(nth,np,nr,nsa) = step_len**2/step_time*sqrt(eps_t)
          END DO
        END DO
      END DO
    END DO

    CALL bounce_average_for_Drw(Drwba, Drwlocal)

    !**** Heat diffusion coef calcul
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        DO np = 1, npmax
          DO nth = 1, nthmax
            !**** Heat flux for FOW momentum not heat velocity [2022/3/22]
            heatrw(nr,nsa) = heatrw(nr,nsa)&
                           - (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                           * Drwlocal(nth,np,nr,nsa) &
                           * JI(nth,np,nr,nsa) &
                           * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa)*1.d20 &
                           / (AEE*1.D3)&  !** Unit convert [J] to [keV]
                           * delp(ns)*delthm(nth,np,nr,nsa)
            heatrwav(nr,nsa) = heatrwav(nr,nsa) &
                             - (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                             * Drwba(nth,np,nr,nsa) &
                             * JI(nth,np,nr,nsa) &
                             * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa)*1.d20 &
                             / (AEE*1.D3)& !** Unit convert [J] to [keV]
                             *delp(ns)* delthm(nth,np,nr,nsa)

!            sumVI = sumVI - dfdrhom(nth,np,nr,nsa)&
!                  * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa))/(AEE*1.D3)
!*** Heat flux of heat velocity [2022/3/22]
!            heatrw(nr,nsa) = heatrw(nr,nsa)&
!                           + (vth(nr,nsa)*AMFP(nsa))**2/(2*AMFP(nsa)) &
!                           * Drwlocal(nth,np,nr,nsa) &
!                           * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
!                           / (AEE*1.D3)  !**** Unit convert [J] to [keV]
!            heatrwav(nr,nsa) = heatrwav(nr,nsa) &
!                             + (vth(nr,nsa)*AMFP(nsa))**2/(2*AMFP(nsa)) &
!                             * Drwba(nth,np,nr,nsa) &
!                             * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
!                             / (AEE*1.D3)  !****Unit convert [J] to [keV]
!
          END DO
        END DO
!        heatrw(nr, nsa)  = heatrw(nr, nsa)/sumVI
!        heatrwav(nr,nsa) = heatrwav(nr,nsa)/sumVI
      END DO
    END DO

  END SUBROUTINE Heat_rw_baverage

!===================================================
! SUBROUTINEs for particle diffusion coef
!===================================================
  SUBROUTINE integral_ParticleDif(Sr_part,Dr)
  !---------------------------------------------
  ! calculation of FOW diff coef
  !---------------------------------------------
    USE fpcomm
    USE fowcomm
    USE fowlib
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT) :: Sr_part, Dr
    REAL(rkind),DIMENSION(nrmax,nsamax) :: Na, dNadr
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: dfdrhom,dfdp
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: dfdthm
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: fnsp_l
    REAL(rkind) Drr_, dVI
    INTEGER nth,np,nr,nsa,ns

    !**** initialization ****
    Sr_part(:,:) = 0.d0
    Dr(:,:) = 0.d0
    !****temperature make
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Na(nr,nsa) = rnsl(nr,nsa)*1.d20
      END DO
    END DO
    !****first order derivation
    DO nsa = 1, nsamax
      CALL first_order_derivative(dNadr(:,nsa), Na(:,nsa), rm)
    END DO

    !********** for new dfdrhom module [2022/3/4]
    fnsp_l(:,:,:,:)=0.d0
    DO nsa = 1, nsamax
      ns = ns_nsa(nsa)
      DO nr= 1, nrmax
        DO np = 1, npmax
          IF ( pm(np,ns) > fact_bulk ) exit
          DO nth = 1, nthmax
          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa)
          END DO
        END DO
      END DO
    END DO

    !**** first order derivative
    DO nsa = 1, nsamax
      DO np = 1, npmax
        DO nth = 1, nthmax
          CALL first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp_l(nth,np,:,nsa), rm)
        END DO
      END DO
    END DO
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        DO np = 1, npmax
          CALL first_order_derivative(dfdthm(:,np,nr,nsa), &
                 fnsp_l(:,np,nr,nsa), thetam(:,np,nr,nsa))
        END DO
      END DO
    END DO
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        DO nth = 1, nthmax
          CALL first_order_derivative(dfdp(nth,:,nr,nsa), &
                 fnsp_l(nth,:,nr,nsa), pm(:,nsa))
        END DO
      END DO
    END DO
    !**** END of new dfdrhom module

    !**** Source calculation
    DO nsa = 1, nsamax
      ns = ns_nsa(nsa)
      DO nr = 1, nrmax
        DO np = 1, npmax
          IF ( pm(np,ns) > fact_bulk ) exit
          DO nth = 1, nthmax
            Sr_part(nr,nsa) = Sr_part(nr,nsa) &
                            - 1.d0&
                            * (Drrfow(nth,np,nr,nsa)&
                            * dfdrhom(nth,np,nr,nsa)*1.d20 &
                            + Drpfow(nth,np,nr,nsa) &
                            * dfdp(nth,np,nr,nsa)*1.d20 &
                            + Drtfow(nth,np,nr,nsa) &
                            * dfdthm(nth,np,nr,nsa)*1.d20&
                            ) &
                            * JI(nth,np,nr,nsa) &
                            !* JI(nth,np,nr,nsa) &
                            * delp(ns) * delthm(nth,np,nr,nsa) &

                          + 1.d0&
                          * Frrfow(nth,np,nr,nsa) &
                          * fnsp_l(nth,np,nr,nsa)*1.d20 &
                          * JI(nth,np,nr,nsa) &
                          !* JI(nth,np,nr,nsa) &
                          * delp(ns) * delthm(nth,np,nr,nsa)
          END DO
        END DO
        Dr(nr,nsa) = -Sr_part(nr,nsa)/dNadr(nr,nsa)
      END DO
    END DO

  END SUBROUTINE integral_ParticleDif

!======================================================================
!SUBROUTINEs for heat fluxes
!======================================================================

  SUBROUTINE integral_Heatdiff(heatfow_out,chi_a)
  !--------------------------------------------
  !SUBROUTINE for FOW heat diffuion coefficient
  !calculate only heat diffusion coefficient
  !--------------------------------------------

    USE fpcomm
    USE fowcomm
    USE fowlib

    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT) :: heatfow_out
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT) :: chi_a
    REAL(rkind),DIMENSION(nrmax,nsamax) :: heatfow_l,sumVI_l,sumVI
    REAL(rkind),DIMENSION(nrmax,nsamax) :: Ta, dTadr
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) ::fnsp_l, dfdrhom
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) ::dfdthm, dfdp
    INTEGER nth,np,nr,nsa,ns

    !**** initialization ****
    heatfow_out(:,:) = 0.d0
    sumVI_l(:,:) = 0.d0
    fnsp_l(:,:,:,:)=0.d0

    !****temperature make
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
        !Ta(temperature)[keV] rwsl[1/m^3*MJ]
      END DO
    END DO
    !****first order derivation
    DO nsa = 1, nsamax
      CALL first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
    END DO
!********** for new dfdrhom module [2022/3/4]
    DO nsa = 1, nsamax
      ns = ns_nsa(nsa)
      DO nr= 1, nrmax
        DO np = 1, npmax
          IF ( pm(np,ns) > fact_bulk ) exit
          DO nth = 1, nthmax
          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa)
          END DO
        END DO
      END DO
    END DO

    !**** first order derivative
    DO nsa = 1, nsamax
      DO np = 1, npmax
        DO nth = 1, nthmax
          CALL first_order_derivative(dfdrhom(nth,np,:,nsa), &
                 fnsp_l(nth,np,:,nsa), rm)
        END DO
      END DO
    END DO
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        DO np = 1, npmax
          CALL first_order_derivative(dfdthm(:,np,nr,nsa), &
                 fnsp_l(:,np,nr,nsa), thetam(:,np,nr,nsa))
        END DO
      END DO
    END DO
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        DO nth = 1, nthmax
          CALL first_order_derivative(dfdp(nth,:,nr,nsa), &
                 fnsp_l(nth,:,nr,nsa), pm(:,nsa))
        END DO
      END DO
    END DO
!!**** END of new dfdrhom module

    !****Integration over moment(np) and pitch angle(nth)
    DO nsa = 1, nsamax
      ns = ns_nsa(nsa)
      DO nr = 1, nrmax
        DO np = 1, npmax
          IF ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
          DO nth = 1, nthmax

           !**** for dfdrhom is made of fnsp*dVI
           heatfow_out(nr,nsa) = heatfow_out(nr,nsa) &
                               - (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa)) &
                               / (AEE*1.D3)*2.d0/3.d0 * 1.d0&
                               *(Drrfow(nth,np,nr,nsa)&
                               * dfdrhom(nth,np,nr,nsa)*1.d20 &
                               + Drpfow(nth,np,nr,nsa) &
                               * dfdp(nth,np,nr,nsa)*1.d20 &
                               + Drtfow(nth,np,nr,nsa) &
                               * dfdthm(nth,np,nr,nsa)*1.d20&
                               ) &
                               * JI(nth,np,nr,nsa) &
                               !* JI(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa) &

                               + (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa) )&
                               / (AEE*1.D3)*2.d0/3.d0 *1.d0&
                               ! unit converter [J] to [keV] &
                               * Frrfow(nth,np,nr,nsa) &
                               * fnsp_l(nth,np,nr,nsa)*1.d20 &
                               * JI(nth,np,nr,nsa) &
                               !* JI(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
                               !** unit [keV] [22/6/6]
          END DO
        END DO
        chi_a(nr,nsa) = -heatfow_out(nr,nsa)/dTadr(nr,nsa)
      END DO
    END DO

  END SUBROUTINE integral_Heatdiff

!===============================================
! SUBROUTINE for factor check
!===============================================
  SUBROUTINE check_factorneo(eps_t,cyclo_rho, dTadr,Ta,tau_ele,tau_i, spVI,spVI2, sTI)
  !-------------------------------------
  ! calculate main factors in neoclassical theory
  !--------------------------------------

    USE fpcomm
    USE fowcomm
    USE fowlib

    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: cyclo_rho, dTadr
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT) :: Ta, spVI,spVI2, sTI
    REAL(rkind),DIMENSION(nrmax),INTENT(OUT) :: eps_t, tau_ele, tau_i
    REAL(rkind) Baxis, fact
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: dfdrhom,dfdrhom2, fnsp_l,fnsp_l2
    INTEGER nth,np,nr,nsa,ns

    Baxis = Bing(1)
    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    !****temperature make
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
        !Ta(temperature)[keV] rwsl[1/m^3*MJ]
        cyclo_rho(nr,nsa) = sqrt(Ta(nr,nsa)*RKEV/AMFP(nsa)) &
                          * AMFP(nsa)/(AEFP(nsa)*Baxis)
        !****Unit [J]
      END DO
    END DO

    !****first order derivation
    DO nsa = 1, nsamax
      CALL first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
    !  CALL first_order_derivative(dndr(:,nsa), rnsl(:, nsa), rm)
    END DO


    !**** Calculation of fators
    DO nr = 1, nrmax
        eps_t(nr) = rm(nr)*RA/RR
        tau_ele(nr) = fact*sqrt(AMFP(1))*((Ta(nr,1)*RKEV)**1.5d0) &
                    / (rnsl(nr,2)*1.d20*AEFP(2)**2*AEE**2*lnlam(nr,2,1))
        tau_i(nr)   = fact*sqrt(2.d0)*sqrt(AMFP(2)) &
                    * ((Ta(nr,2)*RKEV)**1.5d0) &
                    / (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2))
        !**** *RKEV converts [keV] to [J]
    END DO

    !**** initialization ****
    spVI(:,:) = 0.d0

    !********** for new dfdrhom module [2022/3/4]
    fnsp_l(:,:,:,:)=0.d0
    DO nsa = 1, nsamax
      ns = ns_nsa(nsa)
      DO nr= 1, nrmax
        DO np = 1, npmax
          IF ( pm(np,ns) > fact_bulk ) exit
          DO nth = 1, nthmax
          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa) !* JI(nth,np,nr,nsa)
          fnsp_l2(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa)
          END DO
        END DO
      END DO
    END DO

    !**** first order derivative
    DO nsa = 1, nsamax
      DO np = 1, npmax
        DO nth = 1, nthmax
          CALL first_order_derivative(dfdrhom(nth,np,:,nsa), &
                 fnsp_l(nth,np,:,nsa), rm)
          CALL first_order_derivative(dfdrhom2(nth,np,:,nsa), &
                 fnsp_l2(nth,np,:,nsa), rm)
        END DO
      END DO
    END DO
    !!**** END of new dfdrhom module

    !****Integration over moment(np) and pitch angle(nth)
    DO nsa = 1, nsamax
      ns = ns_nsa(nsa)
      DO nr = 1, nrmax
        DO np = 1, npmax
          IF ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
          DO nth = 1, nthmax
            !**** dfdrhom is made of fnsp*dVI
            spVI(nr,nsa) = spVI(nr,nsa) &
                          + 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
!                          * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                          /AEE/1.d3 &![keV]
                          *JI(nth,np,nr,nsa)&
                          !** unit [J][2022/5/29]
                          * delp(ns) * delthm(nth,np,nr,nsa)
            spVI2(nr,nsa) = spVI2(nr,nsa) &
                          + 2.d0/3.d0*dfdrhom2(nth,np,nr,nsa) &
                          * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                          /AEE/1.d3 &![keV]
                          *JI(nth,np,nr,nsa)&
                          !** unit [J][2022/5/29]
                          * delp(ns) * delthm(nth,np,nr,nsa)
            sTI(nr,nsa) = sTI(nr,nsa) &
                          + 2.d0/3.d0*fnsp_l2(nth,np,nr,nsa) &
                          * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                          /AEE/1.d3 &![keV]
                          *JI(nth,np,nr,nsa)&
                          !** unit [J][2022/5/29]
                          * delp(ns) * delthm(nth,np,nr,nsa)
          END DO
        END DO
!        heatfow_out(nr,nsa) = heatfow_out(nr,nsa)/sumVI(nr,nsa)
      END DO
    END DO

  END SUBROUTINE check_factorneo

  SUBROUTINE ch_intHD(hfowout_r,hfowout_p,hfowout_t,hfowout_f,&
                                      chi_Dr,chi_Dp, chi_Dt, chi_Fr)
  !--------------------------------------------
  !SUBROUTINE for factor check of FOW heat currents
  !calculate only heat diffusion coefficient elements
  !--------------------------------------------

    USE fpcomm
    USE fowcomm

    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT) :: hfowout_p, &
                                      hfowout_r,hfowout_t, hfowout_f
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT) :: chi_Dp, &
                                      chi_Dt, chi_Dr, chi_Fr
    REAL(rkind),DIMENSION(nrmax,nsamax) :: heatfow_l,sumVI_l
    REAL(rkind),DIMENSION(nrmax,nsamax) :: Ta, dTadr
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) ::fnsp_l, dfdrhom
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) ::dfdthm, dfdp
    INTEGER nth,np,nr,nsa,ns

    !**** initialization ****
    hfowout_p(:,:) = 0.d0
    hfowout_r(:,:) = 0.d0
    hfowout_t(:,:) = 0.d0
    hfowout_f(:,:) = 0.d0
    chi_Dp(:,:) = 0.d0
    chi_Dr(:,:) = 0.d0
    chi_Dt(:,:) = 0.d0
    chi_Fr(:,:) = 0.d0
    !****temperature make
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
        !Ta(temperature)[keV] rwsl[1/m^3*MJ]
      END DO
    END DO
    !****first order derivation
    DO nsa = 1, nsamax
      CALL first_order_derivative(dTadr(:,nsa), Ta(:,nsa), rm)
    END DO

    !********** for new dfdrhom module [2022/3/4]
    fnsp_l(:,:,:,:)=0.d0
    DO nsa = 1, nsamax
      ns = ns_nsa(nsa)
      DO nr= 1, nrmax
        DO np = 1, npmax
          IF ( pm(np,ns) > fact_bulk ) exit
          DO nth = 1, nthmax
!          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa) * JI(nth,np,nr,nsa)
          fnsp_l(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa) * 1.d20
          END DO
        END DO
      END DO
    END DO

    !**** first order derivative
    DO nsa = 1, nsamax
      DO np = 1, npmax
        DO nth = 1, nthmax
          CALL first_order_derivative(dfdrhom(nth,np,:,nsa), &
                 fnsp_l(nth,np,:,nsa), rm)
        END DO
      END DO
    END DO
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        DO np = 1, npmax
          CALL first_order_derivative(dfdthm(:,np,nr,nsa), &
                 fnsp_l(:,np,nr,nsa), thetam(:,np,nr,nsa))
        END DO
      END DO
    END DO
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        DO nth = 1, nthmax
          CALL first_order_derivative(dfdp(nth,:,nr,nsa), &
                 fnsp_l(nth,:,nr,nsa), pm(:,nsa))
        END DO
      END DO
    END DO
    !!**** END of new dfdrhom module

    !****Integration over moment(np) and pitch angle(nth)
    DO nsa = 1, nsamax
      ns = ns_nsa(nsa)
      DO nr = 1, nrmax
        DO np = 1, npmax
          IF ( pm(np,ns) > fact_bulk ) exit !** fact_balk=5
          DO nth = 1, nthmax

           !**** for dfdrhom is made of fnsp*dVI
           hfowout_r(nr,nsa) = hfowout_r(nr,nsa) &
                               - (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa)) &
                               /(AEE*1.D3) &
                               * Drrfow(nth,np,nr,nsa) &
                               * 2.d0/3.d0*dfdrhom(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
                               !** unit [keV]
           hfowout_p(nr,nsa) = hfowout_p(nr,nsa) &
                               - (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa)) &
                               /(AEE*1.D3) &
                               * Drpfow(nth,np,nr,nsa) &
                               * 2.d0/3.d0*dfdp(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
           hfowout_t(nr,nsa) = hfowout_t(nr,nsa) &
                               - (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa)) &
                               /(AEE*1.D3) &
                               * Drtfow(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa)&
                               * JI(nth,np,nr,nsa) &
                               * 2.d0/3.d0*dfdthm(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
           hfowout_f(nr,nsa) = hfowout_f(nr,nsa) &
                               + (pm(np,nsa)*ptfp0(nsa))**2 &
                               / (2*AMFP(nsa) )&
                               /(AEE*1.D3) &
                               ! unit converter [J] to [keV] &
                               * Frrfow(nth,np,nr,nsa) &
                               * 2.d0/3.d0*fnsp_l(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa) &
                               * JI(nth,np,nr,nsa) &
                               * delp(ns) * delthm(nth,np,nr,nsa)
          END DO
        END DO
        chi_Dr(nr,nsa) = hfowout_r(nr,nsa)!/dTadr(nr,nsa)
        chi_Dp(nr,nsa) = hfowout_p(nr,nsa)!/dTadr(nr,nsa)
        chi_Dt(nr,nsa) = hfowout_t(nr,nsa)!/dTadr(nr,nsa)
        chi_Fr(nr,nsa) = hfowout_f(nr,nsa)!/dTadr(nr,nsa)
      END DO
    END DO

  END SUBROUTINE ch_intHD

  ! *** Modules of collision flecuencies ***

  SUBROUTINE cal_nu_ei(nu)
    USE fpcomm
    USE fowcomm
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT) :: nu
    REAL(rkind) fact, Te, tau
    INTEGER nsa, nr

    fact = 12.d0*pi**1.5d0/sqrt(2.d0)*eps0**2/aee**2*SQRT(AME)
    DO nsa = 1, nsamax
       DO nr = 1, nrmax
          Te = rwsl(nr,1)*1.d6/( 1.5d0*rnsl(nr,1)*1.d20 )
          tau = fact*Te**1.5d0/( lnlam(nr,2,1)*rnsl(nr,2)*1.d20*aefp(2)**2 )
          nu(nr,nsa) = 1.d0/tau !** [J]
       END DO
    END DO

  END SUBROUTINE cal_nu_ei

  SUBROUTINE nu_effective(nu_eff)
  !------------------------------------------------------------------
  ! Calculation of effective collisional flecency[keV] in "Tokamaks"
  ! Formulae from "Tokamaks fourth edition"
  !------------------------------------------------------------------
    USE fpcomm
    USE fowcomm

    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT) :: nu_eff
    REAL(rkind),DIMENSION(nrmax,nsamax) :: Ta
    REAL(rkind) :: fact, tau_e, tau_i
    INTEGER nr, nsa

    !**** initialization
    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)

    !**** Temperature and flecency make
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/( 1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
        !**** Ta(temperature)[keV] /RKEV =*1/( AEE*1.D3)
        IF (nsa == 1) THEN
          nu_eff(nr,nsa) = (1.D0/6.4D14)*rnsl(nr,nsa)*1.D20/(Ta(nr,nsa)**1.5D0)
          nu_eff(nr,nsa) = nu_eff(nr, nsa)/(rm(nr)*RA/RR)
        ELSE
          nu_eff(nr,nsa) = (1.D0/6.6D17)*rnsl(nr,nsa) &
                         * 1.D20/(Ta(nr,nsa)**1.5D0) &
                         * sqrt(AMP/AMFP(nsa)) &
                         * lnlam(nr,nsa,nsa)
                         !** AMP : proton mass
          nu_eff(nr,nsa) = nu_eff(nr, nsa)/(rm(nr)*RA/RR)
        END IF
        !**** Effective collisional flecuency nu
      END DO
    END DO
  END SUBROUTINE nu_effective

  SUBROUTINE nu_deflection(nud)
    USE fpcomm
    USE fowcomm

    IMPLICIT NONE
    REAL(rkind),DIMENSION(npmax,nrmax,nsamax),INTENT(OUT) :: nud
    REAL(rkind) x, fx, gx
    REAL(rkind) nudb, fact, vthb, Tb, va, Ta, pv
    REAL(rkind) nuBra, factBra, numerator, denominator
    INTEGER nth,np,nr,nsa,nsb,nssa,nssb

    fact = 3.d0/2.d0*SQRT( pi/2.d0 )
    factBra = 1.d0/( 12.d0*pi**1.5d0*eps0**2 )

    DO nsa = 1, nsamax
      nssa = ns_nsa(nsa)
      DO nr = 1, nrmax
        DO np = 1, npmax
          nud(np,nr,nsa) = 0.d0
          pv = SQRT(1.d0+theta0(nsa)*pm(np,nsa)**2)
          !** particle verocity [m/s]

          va = vc * SQRT( 1.d0-pv**(-2) )
          !** vc is speed of light[m/s]

          DO nsb = 1, nsbmax
            nssb = ns_nsb(nsa)
            Tb = rwsl(nr,nsb)*1.d6/( 1.5d0*rnsl(nr,nsb)*1.d20 )
            !** Tb[J]
            vthb = SQRT( 2.d0*Tb/amfp(nsb) )

            numerator  = rnsl(nr,nsb)*1.d20 &
                 *AEFP(nsa)**2*AEFP(nsb)**2 * lnlam(nr,nsb,nsa)
            denominator= SQRT( amfp(nsb) )*Tb**1.5d0*SQRT(amfp(nsa)/amfp(nsb))
            nuBra = factBra*numerator/denominator

            x = va/vthb
            CALL gosakn(x,fx,gx)

            nudb = fact*nuBra*(fx-gx)/x**3
            nud(np,nr,nsa) = nud(np,nr,nsa)+nudb

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE nu_deflection

  ! *** Utilities ***

  SUBROUTINE gosakn(x,fx,gx)
    USE fpcomm,only:pi,rkind
    IMPLICIT NONE
    REAL(rkind),INTENT(IN)  :: x
    REAL(rkind),INTENT(OUT) :: fx,gx
    REAL(rkind) :: f,fx1,fx2
    REAL(rkind) :: rh
    DATA RH / 0.70710678118654752440D+00/

    f= DERF(X)*0.5D0
    fx=2.d0*f
    fx1=2.d0/sqrt(pi)*exp(-x**2)
    fx2=-2.d0*x*fx1
    gx=0.d0

    IF(abs(x) .gt. 1.e-10) THEN
    !** .gt. is >,1.e0 is single precision

         gx=(fx-x*fx1)/(2.d0*x**2)
    END IF

    return

  END SUBROUTINE gosakn

end module fowevalNC
