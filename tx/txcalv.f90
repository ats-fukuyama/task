!     $Id$
module variables
  use commons
  implicit none
  public

contains

!***************************************************************
!
!        Calculate variables
!
!***************************************************************

  SUBROUTINE TXCALV(XL)

    use libraries, only : INTG_P, DERIV3SB
    REAL(8), DIMENSION(NQM,0:NRMAX), INTENT(INOUT) :: XL
    INTEGER :: NR
    real(8) :: SUML

    ErV  (0:NRMAX) = XL(LQm1,0:NRMAX)
    EthV (0:NRMAX) = XL(LQm2,0:NRMAX)
    EphV (0:NRMAX) = XL(LQm3,0:NRMAX)
    BthV (0)       = 0.D0
    SUML = 0.D0
    DO NR = 1, NRMAX
       SUML = SUML + INTG_P(XL(LQm4,0:NRMAX),NR)
       BthV(NR) = SUML / R(NR)
    END DO
    BphV (0:NRMAX) = XL(LQm5,0:NRMAX)
    PNeV (0:NRMAX) = XL(LQe1,0:NRMAX)
    UerV (0:NRMAX) = XL(LQe2,0:NRMAX)/PNeV(0:NRMAX)
    UethV(0:NRMAX) = XL(LQe3,0:NRMAX)/PNeV(0:NRMAX)*R(0:NRMAX)
    UephV(0:NRMAX) = XL(LQe4,0:NRMAX)/PNeV(0:NRMAX)
    PTeV (0:NRMAX) = XL(LQe5,0:NRMAX)/PNeV(0:NRMAX)
    PNiV (0:NRMAX) = XL(LQi1,0:NRMAX)
    UirV (0:NRMAX) = XL(LQi2,0:NRMAX)/PNiV(0:NRMAX)
    UithV(0:NRMAX) = XL(LQi3,0:NRMAX)/PNiV(0:NRMAX)*R(0:NRMAX)
    UiphV(0:NRMAX) = XL(LQi4,0:NRMAX)/PNiV(0:NRMAX)
    PTiV (0:NRMAX) = XL(LQi5,0:NRMAX)/PNiV(0:NRMAX)
    PNbV (0:NRMAX) = XL(LQb1,0:NRMAX)
    DO NR = 0, NRMAX
       IF (ABS(PNbV(NR)) < 1.D-7) THEN
          UbthV(NR) = 0.D0
          UbphV(NR) = 0.D0
       ELSE
          UbthV(NR) = XL(LQb3,NR)/PNbV(NR)*R(NR)
          UbphV(NR) = XL(LQb4,NR)/PNbV(NR)
       END IF
    END DO
    PN01V(0:NRMAX) = XL(LQn1,0:NRMAX)
    PN02V(0:NRMAX) = XL(LQn2,0:NRMAX)
    do nr=0,nrmax
!       write(6,'(I3,4F18.10)') nr,XL(LQm1,NR),XL(LQm3,NR),XL(LQe3,NR),XL(LQe3,NR)
!       write(6,'(I3,2F20.10)') nr,XL(LQe4,NR),UephV(NR)
!       write(6,'(I3,3F18.10)') nr,XL(LQi1,NR),XL(LQi4,NR),UiphV(NR)
!       write(6,'(I3,3F18.10)') nr,XL(LQm2,NR),XL(LQm3,NR),XL(LQm4,NR)
!       write(6,'(I3,2F18.10)') nr,XL(LQe1,NR),XL(LQi1,NR)
    end do

    RdBthV(0:NRMAX) = XL(LQm4,0:NRMAX) * R(0:NRMAX) - BthV(0:NRMAX)

    Q(1:NRMAX) = ABS(R(1:NRMAX) * BphV(1:NRMAX) / (RR * BthV(1:NRMAX)))
    Q(0) = (4.D0 * Q(1) - Q(2)) / 3.D0

    RETURN
  END SUBROUTINE TXCALV

!***************************************************************
!
!   Calculate coefficients
!
!***************************************************************

  SUBROUTINE TXCALC

    USE physical_constants, only : AEE, AME, VC, PI, rMU0, EPS0, rKEV
    use libraries, only : EXPV, VALINT_SUB

    INTEGER :: NR, NP, NR1
    REAL(8) :: Sigma0, QL, SL, PNB0, PRFe0, PRFi0, Vte, Vti,  &
         &     rLnLam, EION, XXX, SiV, ScxV, Wte, Wti, EpsL, &
         &     rNuAsE_inv, rNuAsI_inv, BBL, Va, Wpe2, rGC, dQdr, SP, rGBM, &
         &     rGIC, rH, dErdr, dBetadr, &
         &     DCDBM, DeL, ETA, AJPH, AJTH, BN, AJPARA, EPARA, Vcr, Y, &
         &     Ubst, Cs, RhoIT, ExpArg, AiP, DISTAN, &
         &     SiLCL, SiLCthL, SiLCphL
    real(8) :: DERIV3
    real(8), dimension(0:NRMAX) :: Beta, Vexbr, SL1, SL2

    !     *** Constants ***

    !     Neutral cross section

    Sigma0 = 8.8D-21

    !     NBI beam velocity

    Vb =  SQRT(2.D0 * Eb * rKEV / AMB)

    !     Poloidal magnetic field on wall

    IF(FSHL == 0.D0) THEN
       Bthb = rMU0 * rIP * 1.D6 / (2.D0 * PI * RB)
    ELSE
       QL=(Q0-QA)*(1.D0-(RB/RA)**2)+QA
       Bthb = BB*RB/(QL*RR)
    END IF

    !     *** Normalization factor for heating profile ***
    !
    !    SL is a normalization factor for a given heating profile.
    !    It is assumed that the heating profile has a shape of 
    !    exp(-r^2/r_NB^2)*(1-(r/a)^4), in order to renormalize the
    !    heating profile we therefore need to integrate that profile with 
    !    maximum value of unity at axis in advance and calculate the 
    !    normalized factor (i.e. PNB0) which allows an integration
    !    value of renormalized profile to equal designated value (i.e. PNBH).
    !    The factor (1-(r/a)^4), which is an arbitrary function, is
    !    required to damp the profile of exp(-r^2/r_NB^2) rapidly
    !    because if this factor is not imposed unrealistically large current
    !    driven by NBI is generated near the edge where electron and bulk ion
    !    density are dilute although in this case fast ion density from 
    !    NBI is relatively large.

    !  For NBI heating
    SL1(0:NRA) = EXP(- (R(0:NRA) / RNB)**2) * (1.D0 - (R(0:NRA) / RA)** 4)
    CALL VALINT_SUB(SL1,NRA,SL)
    SL = 2.D0 * PI * SL

    PNB0 = PNBH * 1.D6 / (2.D0 * Pi * RR * SL)

    !  For RF heating
    SL2(0:NRA) = EXP(- (R(0:NRA) / RRF)**2) * (1.D0 - (R(0:NRA) / RA)** 4)
    CALL VALINT_SUB(SL2,NRA,SL)
    SL = 2.D0 * PI * SL

    PRFe0 = 0.5D0 * PRFH * 1.D6 / (2.D0 * Pi * RR * SL)
    PRFi0 = 0.5D0 * PRFH * 1.D6 / (2.D0 * Pi * RR * SL)

    Beta(0:NRMAX) = (PNeV(0:NRMAX) * PTeV(0:NRMAX) + PNiV(0:NRMAX) * PTiV(0:NRMAX)) &
         &        * 1.D20 * rKeV /((BphV(0:NRMAX)**2 + BthV(0:NRMAX)**2) / (2.D0 * rMU0))
    Vexbr(1:NRMAX) = ErV(1:NRMAX) &
         &         / (R(1:NRMAX) * SQRT(BphV(1:NRMAX)**2 + BthV(1:NRMAX)**2))

    !     ***** Integer mesh *****

    L_NR: DO NR = 0, NRMAX

       Vte = SQRT(2.D0 * ABS(PTeV(NR)) * rKeV / AME)
       Vti = SQRT(2.D0 * ABS(PTiV(NR)) * rKeV / AMI)
       rLnLam = 15.D0 - LOG(ABS(PNeV(NR))) / 2.D0 + LOG(ABS(PTeV(NR)))

       !     *** Ionization frequency ***

       EION = 13.64D0
       XXX = MAX(PTeV(NR) * 1.D3 / EION, 1.D-2)
       SiV = 1.D-11 * SQRT(XXX) * EXP(- 1.D0 / XXX) &
            &              / (EION**1.5D0 * (6.D0 + XXX))
       rNuION(NR) = FSION * SiV * (PN01V(NR) + PN02V(NR)) * 1.D20

       !     *** Slow neutral diffusion coefficient ***

       D01(NR) = FSD0 * V0**2 &
            &   / (Sigma0 * (PN01V(NR) * V0  + (PNiV(NR) + PN02V(NR)) &
            &      * Vti) * 1.D20)

       !     *** Fast neutral diffusion coefficient ***

       D02(NR) = FSD0 * Vti**2 &
            &   / (Sigma0 * (PNiV(NR) + PN01V(NR) + PN02V(NR)) &
            &      * Vti * 1.D20)

       !     *** Charge exchange frequency ***

       XXX = LOG10(MAX(PTiV(NR) * 1.D3, 50.D0))
       ScxV = 1.57D-16 * SQRT(PTiV(NR) * 1.D3) &
            &          * (XXX * XXX - 14.63D0 * XXX + 53.65D0)
       rNuiCX(NR) = FSCX * ScxV * (PN01V(NR) + PN02V(NR)) * 1.D20

       !     *** Collision frequency (with neutral) ***

       rNu0e(NR) = (PN01V(NR) + PN02V(NR)) * 1.D20 * Sigma0 * Vte
       rNu0i(NR) = (PN01V(NR) + PN02V(NR)) * 1.D20 * Sigma0 * Vti

       !     *** Collision frequency (momentum transfer) ***

       rNuei(NR) = PNiV(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
            &     / (3.D0 * (2.D0 * PI)**1.5D0 &
            &     * EPS0**2 * SQRT(AME) &
            &     * (ABS(PTeV(NR)) * rKeV)**1.5D0)
       rNuii(NR) = PNiV(NR) * 1.D20 * PZ**4 * AEE**4 * rLnLam &
            &     / (3.D0 * SQRT(2.D0) * (2.D0 * PI)**1.5D0 &
            &     * EPS0**2 * SQRT(AMI) &
            &     * (ABS(PTiV(NR)) * rKeV)**1.5D0)
       rNuTei(NR) = PNiV(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
            &     / (3.D0 * SQRT(2.D0 * PI) * PI * EPS0**2 * AME * AMI &
            &     * (  ABS(MAX(PTeV(NR),PTeDIV)) * rKeV / AME &
            &        + ABS(PTiV(NR)) * rKeV / AMI)**1.5D0)

       !     *** Toroidal neoclassical viscosity ***

       Wte = Vte / (Q(NR) * RR) ! Omega_Te
       Wti = Vti / (Q(NR) * RR) ! Omega_Ti
       EpsL = R(NR) / RR        ! Inverse aspect ratio
       rNuAsE_inv = EpsL**1.5D0 * Wte / rNuei(NR)
       rNuAsI_inv = EpsL**1.5D0 * Wti / rNuii(NR)
       BBL = SQRT(BphV(NR)**2 + BthV(NR)**2)
       rNueNC(NR) = FSNC * SQRT(PI) * Q(NR)**2 &
            &     * Wte * 1.78D0 / (rNuAsE_inv + 1.78D0)
       rNuiNC(NR) = FSNC * SQRT(PI) * Q(NR)**2 &
            &     * Wti * 1.78D0 / (rNuAsI_inv + 1.78D0)
       !     &               * (1 + EpsL**1.5D0 * rNuAsI)
       !     &               / (1 + 1.44D0
       !     &                      * ((EpsL**1.5D0 * rNuAsI)**2
       !     &                         + ( ErV(NR)
       !     &                             / ( Vti * BthV(NR)) )**2))
       !!     &                         + ( ErV(NR) * BBL
       !!     &                             / ( Vti * BthV(NR)**2) )**2))

       !     *** Helical neoclassical viscosity ***

       Wte = Vte * NCphi / RR
       Wti = Vti * NCphi / RR
       EpsL = EpsH * R(NR) / RA
       rNuAsE_inv = EpsL**1.5D0 * Wte / rNuei(NR)
       rNuAsI_inv = EpsL**1.5D0 * Wti / rNuii(NR)
       rNueHL(NR) = FSHL * SQRT(PI) &
            &     * Wte * 1.78D0 / (rNuAsE_inv + 1.78D0)
       rNuiHL(NR) = FSHL * SQRT(PI) &
            &     * Wti * 1.78D0 / (rNuAsI_inv + 1.78D0)
       !         WRITE(6,'(I5,1P4E12.4)') 
       !     &        NR,rNueNC(NR),rNuiNC(NR),rNueHL(NR),rNuiHL(NR)
       !         rNueHL(NR) = 0.D0
       !         rNuiHL(NR) = 0.D0
       !     &               * (1 + EpsL**1.5D0 * rNuAsI)
       !     &               / (1 + 1.44D0
       !     &                      * ((EpsL**1.5D0 * rNuAsI)**2
       !     &                         + ( ErV(NR)
       !     &                             / ( Vti * BthV(NR)) )**2))
       !!     &                         + ( ErV(NR) * BBL
       !!     &                             / ( Vti * BthV(NR)**2) )**2))

       !  Derivatives (beta, safety factor, mock ExB velocity)
       dQdr = DERIV3(NR,R,Q,NRMAX,NRM,0)
       S(NR) = R(NR) / Q(NR) * dQdr
       dBetadr = DERIV3(NR,R,Beta,NRMAX,NRM,0)
       Alpha(NR) = - Q(NR)**2 * RR * dBetadr

       !     *** Wave-particle interaction ***

       IF (ABS(FSCDBM) > 0.D0) THEN
          ! Alfven velocity
          Va = SQRT(BBL**2 / (rMU0 * PNiV(NR) * 1.D20 * AMI))
          ! Squared plasma frequency
          Wpe2 = PNeV(NR) * 1.D20 * AEE**2 / (AME * EPS0)
          ! Arbitrary coefficient for CDBM model
          rGC = 8.D0
          ! s-alpha parameter
          IF (Alpha(NR) > 0.D0) THEN
             SP = S(NR) - Alpha(NR)
          ELSE
             SP = Alpha(NR) - S(NR)
          END IF
          ! Fitting function: F for ballooning instability
          IF (SP < 0.D0) THEN
             rGBM = 1.D0 / SQRT(2.D0 * (1.D0 - 2.D0 * SP) &
                  &             * (1.D0 - 2.D0 * SP + 3.D0 * SP**2))
          ELSE
             rGBM = (1.D0 + 9.D0 * SQRT(2.D0) * SP**2.5D0) &
                  &/ (SQRT(2.D0) * (  1.D0 - 2.D0 * SP &
                  &                 + 3.D0 * SP**2 + 2.D0 * SP**3))
          END IF
          ! Magnetic curvature
          rKappa(NR) = - R(NR) / RR * (1.D0 - 1.D0 / Q(NR)**2)
          ! Fitting function: F for interchange instability
          If (Alpha(NR) * rKappa(NR) < 0.D0) THEN
             rGIC = 0.D0
          ELSE
             IF(NR == 0) THEN
                rGIC = 0.D0
             ELSE
                rGIC = ABS(rKappa(NR))**1.5D0 / S(NR)**2
             END IF
          END IF
          ! Select larger effect between ballooning and interchange modes
          FCDBM(NR) = MAX(rGBM, rGIC)
          ! ExB rotational shear
          IF(NR == 0) THEN
             rH=0.D0
          ELSE
             dErdr = DERIV3(NR,R,Vexbr,NRMAX,NRM,0)
             rH = Q(NR) * RR * R(NR) *  dErdr / (Va * S(NR))
          END IF
          
          ! Coefficient of fitting formula
          rG1h2(NR) = 1.D0 / (1.D0 + rG1 * rH**2)
          ! Turbulent transport coefficient calculated by CDBM model
          DCDBM = rGC * FCDBM(NR) * rG1h2(NR) * ABS(Alpha(NR))**1.5D0 &
               &              * VC**2 / Wpe2 * Va / (Q(NR) * RR)
          !write(6,*)DCDBM
          !DCDBM = MAX(DCDBM,1.D-05)
       ELSE
          rG1h2(NR)  = 0.D0
          FCDBM(NR)  = 0.D0
          DCDBM      = 0.D0
       END IF

       IF (R(NR) < RA) THEN
          DeL = FSDFIX * (1.D0 + (PROFD -1) * (R(NR) / RA)**2) &
               &            + FSCDBM * DCDBM
       ELSE
          DeL = FSDFIX * PROFD &
               &       + FSBOHM * PTeV(NR) * rKEV / (16.D0 * AEE * BBL) &
               &       + FSPSCL
       END IF
       ! Particle diffusivity
       De(NR)   = De0   * DeL
       Di(NR)   = Di0   * DeL
       ! Viscosity
       rMue(NR) = rMue0 * DeL
       rMui(NR) = rMui0 * DeL
       ! Thermal diffusivity
       Chie(NR) = Chie0 * DeL
       Chii(NR) = Chii0 * DeL

       ! <omega/m> on half mesh
       WPM(NR) = WPM0 * PTeV(NR) * rKeV / (RA**2 * AEE * BphV(NR))
       ! Force induced by drift wave (eq.(8),(13)) on mesh
       FWthe(NR) = AEE**2         * BphV(NR)**2 * De(NR) &
            &            / (PTeV(NR) * rKeV)
       FWthi(NR) = AEE**2 * PZ**2 * BphV(NR)**2 * Di(NR) &
            &            / (PTiV(NR) * rKeV)

       !     *** Heating profile ***

       IF (R(NR) < RA) THEN
          PNB(NR) = PNB0  * EXP(- R(NR)**2 / RNB**2) * (1 - (R(NR) / RA)** 4)
          SNB(NR) = PNB(NR) / (Eb * rKEV * 1.D20)
          PRFe(NR)= PRFe0 * EXP(- R(NR)**2 / RRF**2) * (1 - (R(NR) / RA)** 4)
          PRFi(NR)= PRFi0 * EXP(- R(NR)**2 / RRF**2) * (1 - (R(NR) / RA)** 4)
       ELSE
          PNB(NR) =0.D0
          SNB(NR) =0.D0
          PRFe(NR)=0.D0
          PRFi(NR)=0.D0
       END IF

       !     *** Current density profile ***

       ! Resistivity
       ETA=AME*rNuei(NR)/(PNeV(NR)*1.D20*AEE**2)
       ! Poloidal current density
       AJPH  = -      AEE * PNeV(NR) * 1.D20 * UephV(NR) &
            &  + PZ * AEE * PNiV(NR) * 1.D20 * UiphV(NR) &
            &  + PZ * AEE * PNbV(NR) * 1.D20 * UbphV(NR)
       ! Toroidal current density
       AJTH  = -      AEE * PNeV(NR) * 1.D20 * UethV(NR) &
            &  + PZ * AEE * PNiV(NR) * 1.D20 * UithV(NR) &
            &  + PZ * AEE * PNbV(NR) * 1.D20 * UbthV(NR)

       BN=SQRT(BthV(NR)**2+BphV(NR)**2)
       ! Parallel current density
       AJPARA=(BthV(NR)*AJTH     + BphV(NR)*AJPH     )/BN
       ! Parallel electric field
       EPARA =(BthV(NR)*EthV(NR) + BphV(NR)*EphV(NR) )/BN
       ! Total current density = parallel current density
       AJ(NR)   = AJPARA
       ! Ohmic current density
       AJOH(NR) = EPARA/ETA
       ! Ohmic heating power
       !POH(NR)  = EPARA*AJPARA
       POH(NR)  = EthV(NR)*AJTH + EphV(NR)*AJPH
       ! NB induced current density
       AJNB(NR) = PZ * AEE * PNbV(NR) * 1.D20 * UbphV(NR)

       !     *** Collision frequency (momentum transfer with beam) ***
       ! reference : 92/04/02, 92/04/21

       Vcr = (3.D0 * SQRT(PI / 2.D0) * PNiV(NR) * PZ**2 / PNeV(NR) * AME / AMI &
            &   * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)**(1.D0/3.D0)
       Y = PNBCD * Vb / Vcr
       IF (Y > 0.D0) THEN
          Ubst = 3.D0 / LOG(1.D0 + Y**3) * Vb
       ELSE
          Ubst = 0.D0
       END IF
       rNube(NR) = PNeV(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
            &     / (3.D0 * (2.D0 * PI)**1.5D0 * EPS0**2 * AMB * AME &
            &             * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)
       rNubi(NR) = PNiV(NR) * 1.D20 * PZ**2 * PZ**2 * AEE**4 * rLnLam &
            &     / (4.D0 * PI * EPS0**2 * AMB) &
            &     * (1.D0 / AMB + 1.D0 / AMI) &
            &     * 1.D0 / ( Ubst**3 + 9.D0 * SQRT(3.D0 * PI) / 4.D0 &
            &     * (ABS(PTiV(NR)) * rKeV / AMI)**1.5D0)
       IF (Y > 0.D0) THEN
          rNuB(NR) = rNube(NR) * 3.D0 / LOG(1.D0 + Y**3)
       ELSE
          rNuB(NR) = 0.D0
       END IF

       !     *** Loss to divertor ***

       IF (R(NR) > RA) THEN
          Cs = SQRT(PTeV(NR) * rKeV / AMI)
!!!!        rNuL(NR) = FSLP * Cs / (2.D0 * PI * Q(NR) * RR &
!!!!     &                  * LOG(0.3D0 / (R(NR) - RA)))
          rNuL(NR) = FSLP * Cs / (2.D0 * PI * Q(NR) * RR &
               &              * (1.D0 + LOG(1.D0 + rLT / (R(NR) - RA)))) &
               &              * (R(NR) - RA) / rLT
       ELSE
          rNuL(NR) = 0.D0
       END IF

    END DO L_NR

    !     ***** Ion Orbit Loss *****

    SiLC  (0:NRMAX) = 0.D0
    SiLCth(0:NRMAX) = 0.D0
    SiLCph(0:NRMAX) = 0.D0
    IF (ABS(FSLC) > 0.D0) THEN
       DO NR = 0, NRA
          EpsL = R(NR) / RR
          Vti = SQRT(2.D0 * PTeV(NR) * rKeV / AMI)
          RhoIT = Vti * AMI / (PZ * AEE * BthV(NR))
          RhoIT = MIN(RhoIT,0.1D0)
          rNuAsI_inv = EpsL**1.5D0 * Vti / rNuii(NR) * Q(NR) * RR
          ExpArg = 2.D0 * EpsL / Vti**2 * (ErV(NR) / BthV(NR))**2
          AiP = rNuii(NR) * SQRT(EpsL) * rNuAsI_inv / (1.D0 + rNuAsI_inv) * EXPV(- ExpArg)
          DO NR1 = NRA, NRMAX
             DISTAN = (R(NR1) - R(NR)) / RhoIT
             SiLCL = AiP * EXPV( - DISTAN**2) * PNiV(NR)
             SiLC(NR) = SiLC(NR) - SiLCL
             SiLC(NR1) = SiLC(NR1) + SiLCL * R(NR) / R(NR1)
             SiLCthL = SiLCL * AMI * UithV(NR)
             SiLCth(NR) = SiLCth(NR) - SiLCthL
             SiLCth(NR1) = SiLCth(NR1) + SiLCthL * R(NR) / R(NR1)
             SiLCphL = SiLCL * AMI * UiphV(NR)
             SiLCph(NR) = SiLCph(NR) - SiLCphL
             SiLCph(NR1) = SiLCph(NR1) + SiLCphL * R(NR) / R(NR1)
          END DO
       END DO
       SiLC  (0:NRMAX) = FSLC * SiLC  (0:NRMAX)
       SiLCth(0:NRMAX) = FSLC * SiLCth(0:NRMAX)
       SiLCph(0:NRMAX) = FSLC * SiLCph(0:NRMAX)
    END IF

    RETURN
  END SUBROUTINE TXCALC
end module variables
