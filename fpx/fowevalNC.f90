MODULE fowevalNC

  PRIVATE

  PUBLIC :: fow_evaluate_nc

CONTAINS

  SUBROUTINE fow_evaluate_nc(Drhmrhm,Drw,Drwav,Dbanana,nud_mono, &
                             Dnewba,Dnewpla,heatfow,heatnewba,heatrw,heatrwav)
    USE fpcomm
    USE fowcomm
    USE fpwrite

    IMPLICIT NONE
    REAL(rkind),INTENT(OUT),dimension(nrmax,nsamax) :: &
         Drhmrhm,Drw,Drwav,Dbanana,nud_mono,Dbanana, &
         Dnewba,Dnewpla,heatfow,heatnewba,heatrw,heatrwav

    REAL(rkind),dimension(nrmax,nsamax) :: nu_D_theory, nu_star_theory,nu_star
    REAL(rkind),dimension(nrmax,nsamax) :: D_mono
    REAL(rkind),dimension(nrmax,nsamax) :: Dnewba,Dnewpla
    REAL(rkind),dimension(nrmax,nsamax) :: jaceffect,heatfow,heatnewba
    INTEGER nth, np, nr, nsa, nsb, mode(3)
   
    !**** make data folder
    CALL system('mkdir -p dat')

    !**** calculation of physical values

    CALL integral_Drr(Drhmrhm)
    CALL integral_Heatdiff(heatfow)

    CALL D_random_walk_baverage(Drwav, Drw)
!    CALL nu_deflection_monoenergy(nud_mono)
!    CALL D_banana(Dbanana)
    CALL newneoclass_ba(Dnewba) ![2022/1/31] editted by anzai
    CALL newneoclass_pla(Dnewpla) ![2022/1/31] editted by anzai
!    CALL check_jeffect(jaceffect) ![2022/2/5] editted by anzai
    CALL newneo_heat_ba(heatnewba) ![2022/2/19] editted by anzai
    CALL Heat_rw_baverage(heatrwav, heatrw) ![2022/2/23] editted by anzai
    
    !**** write txt file
    CALL fptxt2D(Drhmrhm,"dat/Drhmrhm.txt")
    CALL fptxt2D(Drw,"dat/Drw.txt")
    CALL fptxt2D(Drwav,"dat/Drwav.txt")
!    CALL fptxt2D(Dbanana,"dat/Dbanana.txt")
    CALL fptxt2D(Dnewba, "dat/Dnewba.txt") ![2022/1/31] edited by anzai
    CALL fptxt2D(Dnewpla, "dat/Dnewpla.txt") ![2022/1/31] editted by anzai
!    CALL fptxt2D(jaceffect, "dat/jaceffect.txt")![2022/2/5] editted by anzai
    CALL fptxt2D(heatfow, "dat/heatfow.txt")![2022/2/19] editted by anzai
    CALL fptxt2D(heatnewba, "dat/heatnewba.txt")![2022/2/19] editted by anzai
    CALL fptxt2D(heatrw, "dat/heatrw.txt")![2022/2/23] editted by anzai
    CALL fptxt2D(heatrwav, "dat/heatrwav.txt")![2022/2/23] editted by anzai
!    CALL fptxt2D(nud_mono,"dat/nud_mono.txt")

  END SUBROUTINE fow_evaluate_nc

  SUBROUTINE integral_Drr(Drr_out)
    USE fpcomm
    USE fowcomm
    USE fowlib
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT) :: Drr_out
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: dfdrhom
    REAL(rkind) Drr_, dVI, sumVI, JIl_p, JIl_m
    INTEGER nth,np,nr,nsa,ns

    DO nsa = 1, nsamax
      DO np = 1, npmax
        DO nth = 1, nthmax
          CALL first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
        END DO
      END DO
    END DO

    DO nsa = 1, nsamax
      ns = ns_nsa(nsa)
      DO nr = 1, nrmax
        Drr_out(nr,nsa) = 0.d0
        sumVI = 0.d0
        DO np = 1, npmax
          IF ( pm(np,ns) > fact_bulk ) exit
          DO nth = 1, nthmax

            IF ( nr == 1 ) THEN
              JIl_m = JI(nth,np,1,nsa)
              JIl_p = ( JI(nth,np,2,nsa)+JI(nth,np,1,nsa) )*0.5d0
            ELSE IF ( nr == nrmax ) THEN
              JIl_m = ( JI(nth,np,nrmax,nsa)+JI(nth,np,nrmax-1,nsa) )*0.5d0
              JIl_p = JI(nth,np,nrmax,nsa)
            ELSE
              JIl_m = ( JI(nth,np,nr,nsa)+JI(nth,np,nr-1,nsa) )*0.5d0
              JIl_p = ( JI(nth,np,nr,nsa)+JI(nth,np,nr+1,nsa) )*0.5d0
            END IF

            dVI = delp(ns)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)
            Drr_ = ( Drrfow(nth,np,nr+1,nsa)/JIl_p+Drrfow(nth,np,nr,nsa)/JIl_m )*0.5d0

            Drr_out(nr,nsa) = Drr_out(nr,nsa) + Drr_*dfdrhom(nth,np,nr,nsa)*dVI
            sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
          END DO
        END DO
        Drr_out(nr,nsa) = Drr_out(nr,nsa)/sumVI
      END DO
    END DO

  END SUBROUTINE integral_Drr

!  SUBROUTINE D_banana(Dbanana)
!    USE fpcomm
!    USE fowcomm
!    USE fowlib
!    IMPLICIT NONE
!    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: Dbanana
!    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: dfdrhom
!    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
!    REAL(rkind),DIMENSION(nrmax,nsamax) :: nud
!    REAL(rkind) step_len, step_time, eps_t, rho, Baxis, pv, Ta, vtha
!    REAL(rkind) dVI, sumVI
!    INTEGER nth, np, nr, nsa
!
!    Baxis = Bing(1)
!
!    CALL cal_nu_ei(nud)
!
!    DO nsa = 1, nsamax
!      DO np = 1, npmax
!        DO nth = 1, nthmax
!          CALL first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
!        END DO
!      END DO
!    END DO
!
!    DO nsa = 1, nsamax
!      DO nr = 1, nrmax
!        eps_t = rm(nr)*RA/RR
!        Ta = rwsl(nr,nsa)*1.d6/( 1.5d0*rnsl(nr,nsa)*1.d20 )
!        vtha = SQRT( 2.d0*Ta/amfp(nsa) )
!        rho = amfp(nsa)*vtha/(ABS(aefp(nsa))*Baxis)
!        Dbanana(nr,nsa) = safety_factor(nr)**2*rho**2*nud(nr,nsa)/eps_t**1.5d0
!        WRITE(6,'(A,I4,6ES12.4)') 'nr:',nr,eps_t,Ta,vtha,rho,nud(nr,nsa)
!      END DO
!    END DO
!
!  END SUBROUTINE

  SUBROUTINE cal_nu_ei(nu)
  !**************************************
  ![2022/2/8] Modified by anzai
  !**************************************
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
        nu(nr,nsa) = 1.d0/tau
!        WRITE(6,'(A,I6,6ES12.4)') 'nu:',nr,Te,lnlam(nr,2,1),tau,nu(nr,nsa)
      END DO
    END DO
    
  END SUBROUTINE

  SUBROUTINE D_random_walk_baverage(Drwav, Drw)
    USE fpcomm
    USE fowcomm
    USE fowlib
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: Drw, Drwav
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: dfdrhom
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
!    REAL(rkind),DIMENSION(npmax,nrmax,nsamax) :: nud
    REAL(rkind),DIMENSION(nrmax,nsamax) :: nu, Ta
    REAL(rkind) step_len, step_time, eps_t, rho_theta, Baxis, pv, B_theta, rho
    REAL(rkind) dVI, sumVI
    INTEGER nth, np, nr, nsa

    Baxis = Bing(1)

    !CALL nu_deflection(nud)
 !   CALL ca_nu_ei(nu)

    DO nsa = 1, nsamax
      DO np = 1, npmax
        DO nth = 1, nthmax
           CALL first_order_derivative(dfdrhom(nth,np,:,nsa), &
                fnsp(nth,np,:,nsa), rm)
        END DO
      END DO
    END DO
    
    !************************by anzai[2022/2/15]
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/( 1.5d0*rnsl(nr,nsa)*1.d20)  / AEE/1.D3
        !Ta(temperature)[keV]
        nu(nr,nsa) = (1.D0/6.4D14)*rnsl(nr,nsa)*1.D20/(Ta(nr,nsa)**1.5D0)
        nu(nr,nsa) = nu(nr, nsa)/(rm(nr)*RA/RR)*(AMFP(1)/AMFP(nsa)) !Effective nu
      END DO
    END DO
    
    DO nsa = 1, nsamax 
      DO nr = 1, nrmax
        eps_t = rm(nr)*RA/RR
        ! B_theta = safety_factor(nr)/RR*dpsimdr(nr)/RA
        DO np = 1, npmax
          !step_time = 1.d0/nud(np,nr,nsa)
          step_time = 1.d0/nu(nr,nsa)
          rho = pm(np,nsa)*ptfp0(nsa)/(aefp(nsa)*Baxis)
          DO nth = 1, nthmax
            !step_len = rho
            step_len = rho * safety_factor(nr)/sqrt(eps_t)
            !Drwlocal(nth,np,nr,nsa) = step_len**2/step_time*1.46d0*sqrt(eps_t)
            Drwlocal(nth,np,nr,nsa) = step_len**2/step_time*sqrt(eps_t)
          END DO
        END DO
      END DO
    END DO

    CALL bounce_average_for_Drw(Drwba, Drwlocal)

    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Drw(nr,nsa) = 0.d0
        Drwav(nr,nsa) = 0.d0
        sumVI = 0.d0
        DO np = 1, npmax
          DO nth = 1, nthmax
            dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)

            Drw(nr,nsa) = Drw(nr,nsa)&
                          + Drwlocal(nth,np,nr,nsa)*dfdrhom(nth,np,nr,nsa)*dVI
            Drwav(nr,nsa) = Drwav(nr,nsa)&
                          + Drwba(nth,np,nr,nsa)*dfdrhom(nth,np,nr,nsa)*dVI
            sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
          END DO
        END DO
        Drw(nr,nsa) = Drw(nr,nsa)/sumVI
        Drwav(nr,nsa) = Drwav(nr,nsa)/sumVI
      END DO
    END DO

  END SUBROUTINE D_random_walk_baverage

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
          va = vc * SQRT( 1.d0-pv**(-2) )
          DO nsb = 1, nsbmax
            nssb = ns_nsb(nsa)
            Tb = rwsl(nr,nsb)*1.d6/( 1.5d0*rnsl(nr,nsb)*1.d20 )
            vthb = SQRT( 2.d0*Tb/amfp(nsb) )

            numerator  = rnsl(nr,nsb)*1.d20 &
                 *aefp(nsa)**2*aefp(nsb)**2 * lnlam(nr,nsb,nsa)
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

  SUBROUTINE nu_deflection_monoenergy(nud_mono)
    USE fpcomm
    USE fowcomm
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT) :: nud_mono
    REAL(rkind) x, fx, gx
    REAL(rkind) nudb, fact, vthb, Tb, vtha, Ta
    REAL(rkind) nuBra, factBra, numerator, denominator
    INTEGER nr,nsa,nsb,nssa,nssb

    fact = 3.d0/2.d0*SQRT( pi/2.d0 )
    factBra = 1.d0/( 12.d0*pi**1.5d0*eps0**2 )

    DO nsa = 1, nsamax
      nssa = ns_nsa(nsa)
      DO nr = 1, nrmax
        nud_mono(nr,nsa) = 0.d0
        Ta = rwsl(nr,nsa)*1.d6/( 1.5d0*rnsl(nr,nsa)*1.d20 )
        vtha = SQRT( 2.d0*Ta/amfp(nsa) )
        DO nsb = 1, nsbmax
          nssb = ns_nsb(nsa)
          Tb = rwsl(nr,nsb)*1.d6/( 1.5d0*rnsl(nr,nsb)*1.d20 )
          vthb = SQRT( 2.d0*Tb/amfp(nsb) )

          numerator  = rnsl(nr,nsb)*1.d20 * ( aefp(nsa)*aefp(nsb) )**2 * lnlam(nr,nsb,nsa)
          denominator= SQRT( amfp(nsb) ) * Tb**1.5d0 * ( amfp(nsa)/amfp(nsb) )**2
          nuBra = factBra*numerator/denominator

          x = vtha/vthb
          CALL gosakn(x,fx,gx)

          nudb = fact*nuBra*(fx-gx)/x**3

          nud_mono(nr,nsa) = nud_mono(nr,nsa)+nudb

        END DO

      END DO
    END DO

  END SUBROUTINE nu_deflection_monoenergy

  SUBROUTINE gosakn(x,fx,gx)
    USE fpcomm,ONLY:pi,rkind
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
    IF(abs(x).gt.1.e-10) THEN
        gx=(fx-x*fx1)/(2.d0*x**2)
    END IF

    return

  END SUBROUTINE gosakn

  SUBROUTINE bounce_average_for_Drw(Dout, Din)
    USE fpcomm
    USE fowcomm
    USE fowcoef
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax),INTENT(OUT) :: Dout
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax),INTENT(IN)  :: Din
    REAL(rkind),DIMENSION(3,3,max_stp) :: dIdul
    REAL(rkind),allocatable :: U(:,:,:,:,:,:,:,:), Drwl(:,:,:,:,:)
    REAL(rkind) dt, taup, cpitch_ob, psip_ob, thetap_ob, Drw_ob
    TYPE(orbit) ob
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

            Dout(nth,np,nr,nsa) = 0.d0

            DO nstp = 2, nstpmax
              dt = ob%time(nstp)-ob%time(nstp-1)
              cpitch_ob = ob%costh(nstp)
              psip_ob   = ob%psip(nstp)
              thetap_ob = ob%thetap(nstp)
              CALL interpolate_D_unlessZero(Drw_ob, U(:,:,:,:,:,:,np,nsa), &
                    1.d0, cpitch_ob, psip_ob, thetap_ob)
              
              Dout(nth,np,nr,nsa) = Dout(nth,np,nr,nsa)&
                                  + Drw_ob*dIdul(3,3,nstp)**2*dt
            END DO

            Dout(nth,np,nr,nsa) = Dout(nth,np,nr,nsa)/taup
        
          END DO
        END DO
      END DO
    END DO


  END SUBROUTINE
!**********************[2022/1/31] added by anzai***************
  SUBROUTINE newneoclass_ba(Dnewba)
    !----------------------------------
    !SUBROUTINE for neoclassical particle flux in banana region
    !calculate particle flux and diffusion coefficient
    !----------------------------------
    
    USE fpcomm
    USE fowcomm
    USE fowlib
    
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: Dnewba
    REAL(rkind),DIMENSION(nrmax,nsamax) :: Sr_nba, Ta, dTadr, dndr
    !REAL(rkind),DIMENSION(nrmax) :: B_r
    REAL(rkind) fact, fact_s, tau_ele, rho_e, eps_t, Baxis, B_p
    INTEGER nr, nsa

    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bing(1) ! approximation on B by Baxis
 
    !****temperature make
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)
        !**** Ta [J]
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
             (rnsl(nr,2)*1.d20*AEFP(2)**2*lnlam(nr,2,1)*AEE**2)
        rho_e = sqrt(2*Ta(nr,1)/AMFP(1))*AMFP(1)/(AEE*Baxis)
        eps_t = rm(nr)*RA/RR
        !B_p = dpsimdr(nr)*RA/(safety_factor(nr)*RR) ! careful
        fact_s = (safety_factor(nr)**2) &
             *(rho_e**2)/((eps_t**1.5)*tau_ele)
        Sr_nba(nr,nsa) = fact_s*(-1.22d0*(1+Ta(nr,2)/Ta(nr,1))*dndr(nr, nsa)) 
            ! + 4.3d-1*rnsl(nr, nsa)*dTadr(nr,1)/Ta(nr,1) &
            ! + 1.9d-1*rnsl(nr, nsa)*dTadr(nr,2)/Ta(nr,1) ) 

        ! neglect E_para term
        
        Dnewba(nr,nsa) = - Sr_nba(nr, nsa)/dndr(nr, nsa)
!        Dnewba(nr,nsa) = fact_s*1.22d0*(1+Ta(nr,2)/Ta(nr,1))

       ! WRITE(6,'(A,I4,6ES12.4)') 'Dn:',nr, &
            ! TA(nr,1)/AEE,safety_factor(nr),tau_ele,rho_e, &
            ! rho_e**2/tau_ele,Dnewba(nr,nsa)
     END DO
   END DO

  END SUBROUTINE newneoclass_ba

  SUBROUTINE newneoclass_pla(Dnewpla)
    !---------------------------------
    !SUBROUTINE for neoclassical particle flux in plateau region
    !calculate particle flux and diffusion coefficient
    !---------------------------------
    
    USE fpcomm
    USE fowcomm
    USE fowlib
    
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: Dnewpla
    REAL(rkind),DIMENSION(nrmax,nsamax) :: Sr_npla, Ta, dTadr, dndr
    !REAL(rkind),DIMENSION(nrmax) :: B_r
    REAL(rkind) fact, fact_s, tau_ele, rho_e, eps_t, Baxis, B_p
    INTEGER nr, nsa

    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bing(1) ! approximation on B by Baxis

    !****temperature make
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)
        !**** Ta [J]
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
        tau_ele = fact*sqrt(AMFP(1))*(Ta(nr,1)*RKEV)**1.5d0/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**2*AEE**2*lnlam(nr,2,1))
        rho_e = sqrt(Ta(nr,1)*RKEV/AMFP(1))*AMFP(1)/(AEE*Baxis)
        !**** unit [J]
        eps_t = rm(nr)*RA/RR
!        B_p = dpsimdr(nr)*RA/(safety_factor(nr)*RR) ! careful
        B_p = rm(nr)*RA*BB/(safety_factor(nr)*RR) ! careful
        fact_s = -sqrt(pi)*(eps_t**2)*(Ta(nr,1)*RKEV) &
             *rho_e*rnsl(nr,nsa)/(2*AEE*B_p*rm(nr)*RA)        
        
        Sr_npla(nr,nsa) = fact_s*((1+Ta(nr,2)/Ta(nr,1))*dndr(nr,nsa) &
             /(rnsl(nr,nsa)))! &
            ! + 1.5d0*dTadr(nr,1)/Ta(nr,1) + 1.5d0*dTadr(nr,2)/Ta(nr,1) )
        !neglect B term
   
        Dnewpla(nr,nsa) = - Sr_npla(nr, nsa)/dndr(nr, nsa)
      END DO 
    END DO

  END SUBROUTINE newneoclass_pla

  SUBROUTINE check_jeffect(jaceffect)
  !-------------------------------------
  !Integrrate jacobian over theta_m and p
  !--------------------------------------
  
    USE fpcomm
    USE fowcomm
    USE fowlib

    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: jaceffect
     REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: dfdrhom
    REAL(rkind) dVI, sumVI, JIl_p, JIl_m
    INTEGER nth,np,nr,nsa,ns

     !****make dfdrhom
     DO nsa = 1, nsamax
      DO np = 1, npmax
        DO nth = 1, nthmax
          CALL first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
        END DO
      END DO
    END DO
   
    !****integrate
    DO nsa = 1, nsamax
      ns = ns_nsa(nsa)
      DO nr = 1, nrmax
        jaceffect(nr,nsa) = 0.d0
        sumVI = 0.d0
        DO np = 1, npmax
          IF ( pm(np,ns) > fact_bulk ) exit
          DO nth = 1, nthmax
          !-------------------------------------
          !p means plus, m means minus
          !-------------------------------------
            IF ( nr == 1 ) THEN
              JIl_m = JI(nth,np,1,nsa)
              JIl_p = ( JI(nth,np,2,nsa)+JI(nth,np,1,nsa) )*0.5d0
            ELSE IF ( nr == nrmax ) THEN
              JIl_m = ( JI(nth,np,nrmax,nsa)+JI(nth,np,nrmax-1,nsa) )*0.5d0
              JIl_p = JI(nth,np,nrmax,nsa)
            ELSE
              JIl_m = ( JI(nth,np,nr,nsa)+JI(nth,np,nr-1,nsa) )*0.5d0
              JIl_p = ( JI(nth,np,nr,nsa)+JI(nth,np,nr+1,nsa) )*0.5d0
            END IF

            dVI = delp(ns)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)

            jaceffect(nr,nsa) = jaceffect(nr,nsa) + dfdrhom(nth,np,nr,nsa)*dVI
            sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
            
           ! WRITE(*,*) "jacobian:", jaceffect(nr,nsa), sumVI, dfdrhom(nth,np,nr,nsa)
          END DO
        END DO
        jaceffect(nr,nsa) = jaceffect(nr,nsa)/sumVI
      END DO
    END DO

  END SUBROUTINE check_jeffect

!======================================================================
!SUBROUTINEs for heat fluxes
!======================================================================

  SUBROUTINE integral_Heatdiff(heatfow_out)
  !--------------------------------------------
  !SUBROUTINE for FOW heat diffuion coefficient
  !calculate only heat diffusion coefficient
  !--------------------------------------------
    
    USE fpcomm
    USE fowcomm
    USE fowlib
    
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT) :: heatfow_out
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: dfdrhom
    REAL(rkind) Drr_, dVI, sumVI, JIl_p, JIl_m
    INTEGER nth,np,nr,nsa,ns
   
    !****first order derivative
    DO nsa = 1, nsamax
      DO np = 1, npmax
        DO nth = 1, nthmax
          CALL first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
        END DO
      END DO
    END DO

    !****Integration over moment(np) and pitch angle(nth)
    DO nsa = 1, nsamax
      ns = ns_nsa(nsa)
      DO nr = 1, nrmax
        heatfow_out(nr,nsa) = 0.d0
        sumVI = 0.d0
        DO np = 1, npmax
          IF ( pm(np,ns) > fact_bulk ) exit
          DO nth = 1, nthmax

            IF ( nr == 1 ) THEN
              JIl_m = JI(nth,np,1,nsa)
              JIl_p = ( JI(nth,np,2,nsa)+JI(nth,np,1,nsa) )*0.5d0
            ELSE IF ( nr == nrmax ) THEN
              JIl_m = ( JI(nth,np,nrmax,nsa)+JI(nth,np,nrmax-1,nsa) )*0.5d0
              JIl_p = JI(nth,np,nrmax,nsa)
            ELSE
              JIl_m = ( JI(nth,np,nr,nsa)+JI(nth,np,nr-1,nsa) )*0.5d0
              JIl_p = ( JI(nth,np,nr,nsa)+JI(nth,np,nr+1,nsa) )*0.5d0
            END IF

            dVI = delp(ns)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)
            Drr_ = ( Drrfow(nth,np,nr+1,nsa)/JIl_p & 
                 + Drrfow(nth,np,nr,nsa)/JIl_m )*0.5d0

            !**** calculation of heat diffusion coef[keV*m^2/s]
            heatfow_out(nr,nsa) = heatfow_out(nr,nsa) &
                                + (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) & 
                                * Drr_*dfdrhom(nth,np,nr,nsa)*dVI &
                                / (AEE*1.D3)  !****unit convert [J] to [keV]

            sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI &
                 *(pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)*AEE*1.D3) ! AF220220
          END DO
        END DO
        heatfow_out(nr,nsa) = heatfow_out(nr,nsa)/sumVI
      END DO
    END DO

  END SUBROUTINE integral_Heatdiff

  SUBROUTINE newneo_heat_ba(heatnewba)
    !----------------------------------
    !SUBROUTINE for neoclassical heat flux in banana region
    !calculate heat flux and diffusion coefficient
    ! main calcul using[J]
    !----------------------------------

    USE fpcomm
    USE fowcomm
    USE fowlib

    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: heatnewba
    REAL(rkind),DIMENSION(nrmax,nsamax) :: Sr_nba, Ta, dTadr, dndr
    !REAL(rkind),DIMENSION(nrmax) :: B_r
    REAL(rkind), DIMENSION(nsamax) :: rho_a
    REAL(rkind) fact, fact_s, tau_ele,tau_i, eps_t, Baxis, B_p
    INTEGER nr, nsa

    fact = 12.d0*pi**1.5d0*EPS0**2/sqrt(2.d0)
    Baxis = Bing(1) ! approximation on B by Baxis

    !****temperature make
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/(1.5d0*rnsl(nr,nsa)*1.d20)/AEE/1.D3
        !Ta(temperature)[keV]
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

        tau_ele = fact*sqrt(AMFP(1))*((Ta(nr,1)*RKEV)**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**2*AEE**2*lnlam(nr,2,1))
        tau_i   = fact*sqrt(2.d0)*sqrt(AMFP(2))*((Ta(nr,2)*RKEV)**1.5d0)/ &
             (rnsl(nr,2)*1.d20*AEFP(2)**4*lnlam(nr,2,2))
        !**** *RKEV converts [keV] to [J] 
      
        rho_a(nsa) = sqrt(Ta(nr,nsa)*RKEV/AMFP(nsa))*AMFP(nsa)/(AEE*Baxis)
        !**** *RKEV converts [keV] to [J] 
        eps_t = rm(nr)*RA/RR

        !****calculate heat flux
        IF (nsa == 1) THEN
         fact_s = (safety_factor(nr)**2) &
             *(rho_a(nsa)**2)/((eps_t**1.5)*tau_ele)
         Sr_nba(nr,nsa) = fact_s*rnsl(nr,nsa)*1.D20*Ta(nr,1)*RKEV & 
            * ((1.53d0*(1+Ta(nr,2)/Ta(nr,1))*dndr(nr, nsa))/rnsl(nr,nsa) &
              - 1.81d0*dTadr(nr,1)/Ta(nr,1) &
              - 0.27d0*dTadr(nr,2)/Ta(nr,1))
        ! neglect E_para term

        ELSE
         fact_s = (safety_factor(nr)**2) &
              *(rho_a(nsa)**2)/((eps_t**1.5)*tau_i)
         Sr_nba(nr,nsa) = - 0.68d0*fact_s*(1 + 0.48d0*sqrt(eps_t)) &
              *rnsl(nr,nsa)*1.d20*dTadr(nr,2)*RKEV
        END IF  

        heatnewba(nr,nsa) = - Sr_nba(nr, nsa) &
             /(rnsl(nr,nsa)*1.d20*dTadr(nr, nsa)*RKEV)
      END DO
    END DO

  END SUBROUTINE newneo_heat_ba

  SUBROUTINE Heat_rw_baverage(heatrwav, heatrw)
  !---------------------------------------------
  ! heat diffusion coef calcul module
  ! This module corresponds to integral_Heatdiff
  !--------------------------------------------- 
    
    USE fpcomm
    USE fowcomm
    USE fowlib
    
    IMPLICIT NONE
    REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: heatrw, heatrwav
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: dfdrhom
    REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
    REAL(rkind),DIMENSION(nrmax,nsamax) :: nu, Ta, Drw
    REAL(rkind) step_len, step_time, eps_t, rho_theta, Baxis, pv, rho
    REAL(rkind) dVI, sumVI
    INTEGER nth, np, nr, nsa

    Baxis = Bing(1)

    DO nsa = 1, nsamax
      DO np = 1, npmax
        DO nth = 1, nthmax
           CALL first_order_derivative(dfdrhom(nth,np,:,nsa), &
                fnsp(nth,np,:,nsa), rm)
        END DO
      END DO
    END DO

    !**** Temperature and flecency make
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        Ta(nr,nsa) = rwsl(nr,nsa)*1.d6/( 1.5d0*rnsl(nr,nsa)*1.d20) / AEE/1.D3
        !**** Ta(temperature)[keV] RKEV = AEE*1.D3

        nu(nr,nsa) = (1.D0/6.4D14)*rnsl(nr,nsa)*1.D20/(Ta(nr,nsa)**1.5D0)
        nu(nr,nsa) = nu(nr, nsa)/(rm(nr)*RA/RR)*(AMFP(1)/AMFP(nsa)) 
        !**** Effective collisional flecuency nu
      END DO
    END DO
   
    !**** Calculation of particle diff coef
    DO nsa = 1, nsamax
      DO nr = 1, nrmax
        eps_t = rm(nr)*RA/RR
        DO np = 1, npmax
          step_time = 1.d0/nu(nr,nsa)

          rho = pm(np,nsa)*ptfp0(nsa)/(aefp(nsa)*Baxis)
          !****Unit [keV]

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
        heatrw(nr, nsa) = 0.d0
        heatrwav(nr,nsa) = 0.d0
        sumVI = 0.d0
        DO np = 1, npmax
          DO nth = 1, nthmax
            dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)

            heatrw(nr,nsa) = heatrw(nr,nsa)&
                           + (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                           * Drwlocal(nth,np,nr,nsa) & 
                           * dfdrhom(nth,np,nr,nsa)*dVI &
                           / (AEE*1.D3)  !**** Unit convert [J] to [keV]
            heatrwav(nr,nsa) = heatrwav(nr,nsa) &
                             + (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)) &
                             * Drwba(nth,np,nr,nsa) & 
                             * dfdrhom(nth,np,nr,nsa)*dVI &
                             / (AEE*1.D3)  !****Unit convert [J] to [keV]

            sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI &
                  * (pm(np,nsa)*ptfp0(nsa))**2/(2*AMFP(nsa)*AEE*1.D3)
          END DO
        END DO
        heatrw(nr, nsa)  = heatrw(nr, nsa)/sumVI
        heatrwav(nr,nsa) = heatrwav(nr,nsa)/sumVI
      END DO
    END DO

  END SUBROUTINE Heat_rw_baverage

  ! SUBROUTINE banana_width_and_omega_bounce(w_b, omega_b)
  !   USE fpcomm
  !   USE fowcomm
  !   IMPLICIT NONE
  !   REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax),INTENT(OUT) :: w_b, omega_b
  !   REAL(rkind),allocatable :: drdt(:), U(:,:)
  !   TYPE(orbit) ob
  !   INTEGER nth, np, nr, nsa, nstp, nstpmax, ierr
  !   REAL(rkind) r_min

  !   ierr = 0

  !   DO nsa = 1, nsamax
  !     DO nr = 1, nrmax
  !       DO np = 1, npmax
  !         DO nth = 1, nthmax
  !           ob = orbit_m(nth,np,nr,nsa)
  !           nstpmax = ob%nstp_max
  !           omega_b(nth,np,nr,nsa) = 2.d0*pi/ob%time(nstpmax)
  !           allocate(drdt(nstpmax))
  !           allocate(U(4, nstpmax))
  !           CALL spl1D(ob%time, ob%r, drdt, U, nstpmax, 0, ierr)
  !           CALL spl1DF(ob%time(nstpmax)*0.5d0, r_min, ob%time, U ,nstpmax, ierr)
  !           w_b(nth,np,nr,nsa) = ( rm(nr)-r_min )*RA
  !           deallocate(drdt)
  !           deallocate(U)
  !         END DO
  !       END DO
  !     END DO
  !   END DO

  ! END SUBROUTINE banana_width_and_omega_bounce

  ! SUBROUTINE banana_width_and_omega_bounce_model(w_b, omega_b)
  !   USE fpcomm
  !   USE fowcomm
  !   IMPLICIT NONE
  !   REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax),INTENT(OUT) :: w_b, omega_b
  !   REAL(rkind) :: rho, eps,sign_prev,sign,Baxis,vl,pv,B_theta,thetapl
  !   TYPE(orbit) :: ob
  !   logical :: isTrapped
  !   INTEGER :: nth,np,nr,nsa,ns,nstp,nstpmax

  !   Baxis = Boutg(1)

  !   DO nsa = 1, nsamax
  !     ns = ns_nsa(nsa)
  !     DO nr = 1, nrmax
  !       eps = rm(nr)*RA/RR
  !       DO np = 1, npmax
  !         pv = SQRT(1.d0+theta0(nsa)*pm(np,nsa)**2)
  !         vl = vc*SQRT(1.d0-1.d0/pv**2)
  !         DO nth = 1, nthmax
  !           ob = orbit_m(nth,np,nr,nsa)
  !           nstpmax = ob%nstp_max

  !           isTrapped = .false.
  !           IF ( ob%costh(1) < 0.d0 ) THEN
  !             sign_prev = -1.d0
  !           ELSE IF ( ob%costh(1) > 0.d0 ) THEN
  !             sign_prev = 1.d0
  !           ELSE
  !             sign_prev = 0.d0
  !           END IF

  !           DO nstp = 2, nstpmax
  !             IF ( ob%costh(nstp) /= 0.d0 ) THEN
  !               sign = ob%costh(nstp) / ABS( ob%costh(nstp) )
  !             ELSE
  !               sign = 0.d0
  !             END IF

  !             IF ( sign_prev*sign > 0.d0 ) THEN
  !               cycle
  !             ELSE IF ( sign_prev*sign < 0.d0 ) THEN
  !               isTrapped = .True.
  !             ELSE
  !               sign_prev = sign
  !             END IF
  !           END DO

  !           IF ( isTrapped ) THEN
  !             rho = pm(np,ns)*ptfp0(ns)/ABS( pz(ns)*aee*Baxis )
  !             w_b(nth,np,nr,nsa) = safety_factor(nr)*rho/SQRT(eps)
  !             omega_b(nth,np,nr,nsa) = SQRT(eps)*vl/( safety_factor(nr)*RR )
  !           ELSE
  !             IF ( pz(ns)*COS( thetam(nth,np,nr,nsa) ) > 0.d0 ) THEN
  !               thetapl = 0.d0
  !             ELSE
  !               thetapl = pi
  !             END IF
  !             B_theta = safety_factor(nr)/( RR+RA*rm(nr)*COS( thetapl ) )*dpsimdr(nr)
  !             w_b(nth,np,nr,nsa) = 0.d0
  !             omega_b(nth,np,nr,nsa) = B_theta*vl*COS( thetam(nth,np,nr,nsa) )&
  !                                     /( Baxis*rm(nr)*RA )
  !             omega_b(nth,np,nr,nsa) = ABS( omega_b(nth,np,nr,nsa) )
  !             ! omega_b(nth,np,nr,nsa) = 0.d0
  !           END IF

  !         END DO
  !       END DO
  !     END DO
  !   END DO

  ! END SUBROUTINE banana_width_and_omega_bounce_model

  ! SUBROUTINE mean_transmatrix(dIdu,mode)
  !   USE fpcomm
  !   USE fowcomm
  !   USE fowcoef
  !   IMPLICIT NONE
  !   REAL(rkind),DIMENSION(:,:,:,:,:,:),INTENT(OUT) :: dIdu
  !   INTEGER,DIMENSION(3),INTENT(IN) :: mode
  !   REAL(rkind),DIMENSION(3,3,max_stp) :: dIdul
  !   REAL(rkind) :: dt, taup
  !   INTEGER :: nth,np,nr,nsa,nstp,nstpmax,i,j
  !   TYPE(orbit) :: ob

  !   DO nsa = 1, nsamax
  !     DO nr = 1, nrmax+mode(3)
  !       DO np = 1, npmax+mode(2)
  !         DO nth = 1, nthmax+mode(1)
  !           DO j = 1, 3
  !             DO i = 1, 3
  !               dIdu(i,j,nth,np,nr,nsa) = 0.d0
  !             END DO
  !           END DO
  !         END DO
  !       END DO
  !     END DO
  !   END DO

  !   DO nsa = 1, nsamax
  !     DO nr = 1+mode(3), nrmax+mode(3)
  !       DO np = 1+mode(2), npmax+mode(2)
  !         DO nth = 1, nthmax+mode(1)

  !           IF ( mode(1) == 0 .and. mode(2) == 0 .and. mode(3) == 0 ) THEN
  !             ob = orbit_m(nth,np,nr,nsa)
  !           ELSE IF ( mode(1) == 1 .and. mode(2) == 0 .and. mode(3) == 0 ) THEN
  !             ob = orbit_th(nth,np,nr,nsa)
  !             IF ( nth == nth_stg(nsa) ) cycle
  !           ELSE IF ( mode(1) == 0 .and. mode(2) == 1 .and. mode(3) == 0 ) THEN
  !             ob = orbit_p(nth,np,nr,nsa)
  !           ELSE IF ( mode(1) == 0 .and. mode(2) == 0 .and. mode(3) == 1 ) THEN
  !             ob = orbit_r(nth,np,nr,nsa)
  !           END IF

  !           nstpmax = ob%nstp_max
  !           taup = ob%time(nstpmax)
  !           CALL transformation_matrix(dIdul, ob, nth, np, nr, nsa, mode)

  !           DO nstp = 2, nstpmax
  !             dt = ( ob%time(nstp)-ob%time(nstp-1) )/taup
  !             DO j = 1, 3
  !               DO i = 1, 3
  !                 dIdu(i,j,nth,np,nr,nsa) = dIdu(i,j,nth,np,nr,nsa) + dIdul(i,j,nstp)*dt
  !               END DO
  !             END DO

  !           END DO

  !         END DO
  !       END DO
  !     END DO
  !   END DO

  ! END SUBROUTINE  mean_transmatrix

  ! SUBROUTINE D_random_walk_baverage_old(Drw)
  !   USE fpcomm
  !   USE fowcomm
  !   USE fowlib
  !   IMPLICIT NONE
  !   REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: Drw
  !   REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: w_b, omega_b, dfdrhom
  !   REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: Drwlocal, Drwba
  !   REAL(rkind),DIMENSION(npmax,nrmax,nsamax) :: nud
  !   REAL(rkind) step_len, step_time
  !   REAL(rkind) dVI, sumVI
  !   INTEGER nth, np, nr, nsa

  !   CALL banana_width_and_omega_bounce(w_b, omega_b)
  !   CALL nu_deflection(nud)

  !   DO nsa = 1, nsamax
  !     DO np = 1, npmax
  !       DO nth = 1, nthmax
  !         CALL first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
  !       END DO
  !     END DO
  !   END DO

  !   DO nsa = 1, nsamax
  !     DO nr = 1, nrmax
  !       DO np = 1, npmax
  !         step_time = 1.d0/nud(np,nr,nsa)
  !         DO nth = 1, nthmax
  !           step_len = w_b(nth,np,nr,nsa)
  !           Drwlocal(nth,np,nr,nsa) = step_len**2/step_time
  !         END DO
  !       END DO
  !     END DO
  !   END DO

  !   CALL bounce_average_for_Drw(Drwba, Drwlocal)

  !   DO nsa = 1, nsamax
  !     DO nr = 1, nrmax
  !       Drw(nr,nsa) = 0.d0
  !       sumVI = 0.d0
  !       DO np = 1, npmax
  !         DO nth = 1, nthmax
  !           dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)

  !           Drw(nr,nsa) = Drw(nr,nsa)&
  !                         + Drwba(nth,np,nr,nsa)*dfdrhom(nth,np,nr,nsa)*dVI
  !           sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
  !         END DO
  !       END DO
  !       Drw(nr,nsa) = Drw(nr,nsa)/sumVI
  !     END DO
  !   END DO

  ! END SUBROUTINE D_random_walk_baverage_old

  ! SUBROUTINE ionHeatFlux(q_ion)
  !   USE fpcomm
  !   USE fowcomm
  !   USE fowlib
  !   IMPLICIT NONE
  !   REAL(rkind),DIMENSION(nrmax),INTENT(OUT) :: q_ion
  !   REAL(rkind),DIMENSION(nthmax,npmax,nrmax) :: dfdrhom
  !   REAL(rkind) :: Drr_, dVI, K, pv
  !   INTEGER :: nth,np,nr,nsa,ns

  !   DO np = 1, npmax
  !     DO nth = 1, nthmax
  !       CALL first_order_derivative(dfdrhom(nth,np,:), fnsp(nth,np,:,2), rm)
  !       DO nr = 1, nrmax
  !         dfdrhom(nth,np,nr) = dfdrhom(nth,np,nr)/RA
  !       END DO
  !     END DO
  !   END DO

  !   DO nr = 1, nrmax
  !     q_ion(nr) = 0.d0
  !     DO np = 1, npmax
  !       IF ( pm(np,2) > fact_bulk ) exit
  !       DO nth = 1, nthmax
  !         dVI = delp(2)*delthm(nth,np,nr,2)*JIR(nth,np,nr,2)
  !         Drr_ = ( Drrfow(nth,np,nr+1,2)+Drrfow(nth,np,nr,2) )*0.5d0/JI(nth,np,nr,2)
  !         pv = SQRT(1.d0+theta0(2)*pm(np,2)**2)
  !         K = amfp(2)*vc**2*(pv-1.d0)

  !         q_ion(nr) = q_ion(nr) - K*Drr_*dfdrhom(nth,np,nr)*dVI*1.d20
  !       END DO
  !     END DO
  !     write(*,*)"ion",q_ion(nr)
  !   END DO

  ! END SUBROUTINE ionHeatFlux

  ! SUBROUTINE D_neoclass()
  !   !delete D_banana, D_plateau, nu_ei, nu_p, nu_b as variables [2022/1/17]
  !   USE fowcomm
  !   USE fpcomm

  !   IMPLICIT NONE
  !   REAL(rkind),DIMENSION(nrmax),INTENT(OUT) :: D_banana, D_plateau, nu_ei, nu_p, nu_b
  !   INTEGER :: nr, nsb, nssb
  !   REAL(rkind) :: dens ,Te, Vt, omega_e, Baxis, rho_e, epst
  !   REAL(rkind) :: factSpitz

  !   Baxis = Bing(1)
  !   omega_e = aee*Baxis/ame
  !   factSpitz = aee**4/( 93.d0*eps0**2*SQRT(ame) )
  !   DO nr = 1, nrmax
  !     Te = rwsl(nr,1)*1.d6/( 1.5d0*rnsl(nr,1)*1.d20 ) ! [J]
  !     nu_ei(nr) = 0.d0
  !     DO nsb = 1, nsbmax
  !       nssb = ns_nsb(nsb)
  !       IF ( nssb == 1 ) cycle
  !       dens = rnsl(nr,nsb)*1.d20
  !       nu_ei(nr) = nu_ei(nr) + pz(nssb)**2*dens*lnlam(nr,nsb,1)
  !     END DO
  !     nu_ei(nr) = nu_ei(nr)*factSpitz/Te**1.5d0
  !     ! nu_ei = rnud(nr,2,1)

  !     Vt = SQRT( Te/ame )
  !     rho_e = Vt/omega_e
  !     epst = rm(nr)*RA/RR
  !     nu_p(nr) = Vt/( safety_factor(nr)*RR )
  !     nu_b(nr) = epst**1.5d0*nu_p(nr)

  !     !D_banana (nr) = ( safety_factor(nr)*rho_e )**2 * nu_ei(nr) / epst**1.5d0
  !     D_plateau(nr) = ( safety_factor(nr)*rho_e )**2 * nu_p(nr)
  !   END DO

  ! END SUBROUTINE D_neoclass

  ! SUBROUTINE D_random_walk(Drw)
  !   USE fpcomm
  !   USE fowcomm
  !   USE fowlib
  !   IMPLICIT NONE
  !   REAL(rkind),DIMENSION(nrmax,nsamax),INTENT(OUT)  :: Drw
  !   REAL(rkind),DIMENSION(nthmax,npmax,nrmax,nsamax) :: w_b, omega_b, dfdrhom
  !   REAL(rkind),DIMENSION(npmax,nrmax,nsamax)        :: nud
  !   REAL(rkind) step_len, step_time, Drwlocal
  !   REAL(rkind) dVI, sumVI
  !   INTEGER nth, np, nr, nsa

  !   CALL banana_width_and_omega_bounce(w_b, omega_b)
  !   CALL nu_deflection(nud)

  !   DO nsa = 1, nsamax
  !     DO np = 1, npmax
  !       DO nth = 1, nthmax
  !         CALL first_order_derivative(dfdrhom(nth,np,:,nsa), fnsp(nth,np,:,nsa), rm)
  !       END DO
  !     END DO
  !   END DO

  !   DO nsa = 1, nsamax
  !     DO nr = 1, nrmax
  !       Drw(nr,nsa) = 0.d0
  !       sumVI = 0.d0

  !       DO np = 1, npmax
  !         step_time = 1.d0/nud(np,nr,nsa)
  !         DO nth = 1, nthmax
  !           step_len = w_b(nth,np,nr,nsa)/RA
  !           Drwlocal = step_len**2/step_time
            
  !           dVI = delp(nsa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)
  !           sumVI = sumVI + dfdrhom(nth,np,nr,nsa)*dVI
  !           Drw(nr,nsa) = Drw(nr,nsa) + Drwlocal*dfdrhom(nth,np,nr,nsa)*dVI

  !         END DO
  !       END DO

  !       Drw(nr,nsa) = Drw(nr,nsa)/sumVI

  !     END DO
  !   END DO

  ! END SUBROUTINE D_random_walk

END MODULE fowevalNC
