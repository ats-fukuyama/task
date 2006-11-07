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

  SUBROUTINE TXCALV(XL,ID)

    use libraries, only : INTG_P
    REAL(8), DIMENSION(1:NQM,0:NRMAX), INTENT(INOUT) :: XL
    integer, intent(in), optional :: ID
    INTEGER :: NR
    real(8) :: SUML, DERIV3

    Phi  (0:NRMAX) =   XL(LQm1,0:NRMAX)
    DO NR = 0, NRMAX
       ErV (NR) = - 2.D0 * R(NR) * DERIV3(NR,PSI,XL(LQm1,0:NRMAX),NRMAX,NRM,0)
    END DO
    EthV (0)       =   0.D0
    EthV (1:NRMAX) = - XL(LQm2,1:NRMAX) / R(1:NRMAX)
    EphV (0:NRMAX) = - XL(LQm3,0:NRMAX)
    AphV (0:NRMAX) =   XL(LQm4,0:NRMAX)
    RAthV(0:NRMAX) =   XL(LQm5,0:NRMAX)
    DO NR = 0, NRMAX
       BthV(NR) = - 2.D0 * R(NR) * DERIV3(NR,PSI,XL(LQm4,0:NRMAX),NRMAX,NRM,0)
       BphV(NR) =   2.D0         * DERIV3(NR,PSI,XL(LQm5,0:NRMAX),NRMAX,NRM,0)
    END DO
    PNeV (0:NRMAX) =   XL(LQe1,0:NRMAX)
    UerV (0)       =   0.D0
    UerV (1:NRMAX) =   XL(LQe2,1:NRMAX )/ PNeV(1:NRMAX) / R(1:NRMAX)
    UethV(0)       =   0.D0
    UethV(1:NRMAX) =   XL(LQe3,1:NRMAX) / PNeV(1:NRMAX) / R(1:NRMAX)
    UephV(0:NRMAX) =   XL(LQe4,0:NRMAX) / PNeV(0:NRMAX)
    PTeV (0:NRMAX) =   XL(LQe5,0:NRMAX) / PNeV(0:NRMAX)
    PNiV (0:NRMAX) =   XL(LQi1,0:NRMAX)
    UirV (0)       =   0.D0
    UirV (1:NRMAX) =   XL(LQi2,1:NRMAX) / PNiV(1:NRMAX) / R(1:NRMAX)
    UithV(0)       =   0.D0
    UithV(1:NRMAX) =   XL(LQi3,1:NRMAX) / PNiV(1:NRMAX) / R(1:NRMAX)
    UiphV(0:NRMAX) =   XL(LQi4,0:NRMAX) / PNiV(0:NRMAX)
    PTiV (0:NRMAX) =   XL(LQi5,0:NRMAX) / PNiV(0:NRMAX)
    PNbV (0:NRMAX) =   XL(LQb1,0:NRMAX)
    DO NR = 0, NRA-1
       IF(PNbV(NR) == 0.D0) THEN
          UbthV(NR) = 0.D0
          UbphV(NR) = 0.D0
       ELSE
          IF(PNBH == 0.D0 .AND. PNbV(NR) < 1.D-8) THEN
             UbthV(NR) = 0.D0
             UbphV(NR) = 0.D0
          ELSE
             IF(NR == 0) THEN
                UbthV(NR) = 0.D0
                UbphV(NR) = XL(LQb4,NR) / PNbV(NR)
             ELSE
                UbthV(NR) = XL(LQb3,NR) / PNbV(NR) / R(NR)
                UbphV(NR) = XL(LQb4,NR) / PNbV(NR)
             END IF
          END IF
       END IF
    END DO
    IF(PRESENT(ID).AND.ID==1) write(6,*) UbthV(0:NRMAX),UbphV(0:NRMAX)
    UbthV(NRA:NRMAX) = 0.D0
    UbphV(NRA:NRMAX) = 0.D0
!!$    DO NR = 0, NRMAX
!!$       IF (ABS(PNbV(NR)) < 1.D-6) THEN
!!$          UbthV(NR) = 0.D0
!!$          UbphV(NR) = 0.D0
!!$       ELSE
!!$          IF (NR == 0) THEN
!!$             UbthV(NR) = 0.D0
!!$          ELSE
!!$             UbthV(NR) = XL(LQb3,NR) / PNbV(NR) / R(NR)
!!$          END IF
!!$          UbphV(NR) =    XL(LQb4,NR) / PNbV(NR)
!!$       END IF
!!$    END DO
    PN01V(0:NRMAX) =   XL(LQn1,0:NRMAX)
    PN02V(0:NRMAX) =   XL(LQn2,0:NRMAX)

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

    USE physical_constants, only : AEE, AME, VC, PI, rMU0, EPS0, rKeV
    use libraries, only : EXPV, VALINT_SUB, TRCOFS
    use nclass_mod
    
    INTEGER :: NR, NP, NR1, IER, NRPLTE, NRPLTI
    REAL(8) :: Sigma0, QL, SL, PNB0, PRFe0, PRFi0, Vte, Vti, Vtb, &
         &     rLnLam, EION, XXX, SiV, ScxV, Wte, Wti, EpsL, &
         &     rNuAsE_inv, rNuAsI_inv, BBL, Va, Wpe2, rGC, dQdr, SP, rGBM, &
         &     rGIC, rH, dErdr, dpdr, PROFDL, &
         &     DCDBM, DeL, AJPH, AJTH, AJPARA, EPARA, Vcr, &
         &     Cs, RhoIT, ExpArg, AiP, DISTAN, DeLa, &
         &     SiLCL, SiLCthL, SiLCphL, Wbane, Wbani, RL, ALFA, DBW, PTiVA, &
         &     KAPPA, rNuBAR, NGRADB2, K11PSe, K11Be,  K11Pe, K11PSi, K11Bi, K11Pi
    real(8) :: Ce = 0.733D0, Ci = 1.365D0
    real(8) :: FTL, FCL, EFT, ETAS, CR
    real(8) :: DERIV3, AITKEN2P
    real(8), dimension(0:NRMAX) :: p, Vexbr, SL1, SL2

    !     *** Constants ***

    !     Neutral cross section

    Sigma0 = 8.8D-21

    !     NBI beam velocity

    Vb =  SQRT(2.D0 * Eb * rKeV / AMB)

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

    PNB0 = ABS(PNBH) * 1.D6 / (2.D0 * Pi * RR * SL)

    !  For RF heating
    SL2(0:NRA) = EXP(- (R(0:NRA) / RRF)**2) * (1.D0 - (R(0:NRA) / RA)** 4)
    CALL VALINT_SUB(SL2,NRA,SL)
    SL = 2.D0 * PI * SL

    PRFe0 = 0.5D0 * PRFH * 1.D6 / (2.D0 * Pi * RR * SL)
    PRFi0 = 0.5D0 * PRFH * 1.D6 / (2.D0 * Pi * RR * SL)

    p(0:NRMAX) = (PNeV(0:NRMAX) * PTeV(0:NRMAX) + PNiV(0:NRMAX) * PTiV(0:NRMAX)) &
         &        * 1.D20 * rKeV
    Vexbr(1:NRMAX) = ErV(1:NRMAX) &
         &         / (R(1:NRMAX) * SQRT(BphV(1:NRMAX)**2 + BthV(1:NRMAX)**2))

    IF(PROFD == 0.D0.AND.FSDFIX /= 0.D0) THEN
       PROFDL = (PTeV(NRA) * rKeV / (16.D0 * AEE * &
            &    SQRT(BphV(NRA)**2 + BthV(NRA)**2))) / FSDFIX
    ELSE
       PROFDL = PROFD
    END IF

    ! Banana width
    Wbane = (Q0 * SQRT(RR * AME * PTeV(0) * rKeV) / (AEE * BphV(0)))**(2.D0/3.D0)
    Wbani = (Q0 * SQRT(RR * AMI * PTiV(0) * rKeV) / (PZ * AEE * BphV(0)))**(2.D0/3.D0)

    ! Banana width at separatrix
!    PTiVA = PTiV(NRA)
    PTiVA = 0.5D0 * PTiV(0)
    DBW = 3.D0 * SQRT(PTiVA * rKEV * AMI) * Q(NRA) / (PZ * AEE * BphV(NRA)) &
         & / SQRT(R(NRA) / RR)

    !  Coefficients

!    write(6,*) "++++++++++++++++++++++"
    L_NR: DO NR = 0, NRMAX

       Vte = SQRT(2.D0 * ABS(PTeV(NR)) * rKeV / AME)
       Vti = SQRT(2.D0 * ABS(PTiV(NR)) * rKeV / AMI)
       Vtb = SQRT(2.D0 * ABS(PTiV(NR)) * rKeV / AMB)
       rLnLam = 15.2D0 - LOG(ABS(PNeV(NR))) / 2.D0 + LOG(ABS(PTeV(NR)))

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
       rNu0b(NR) = (PN01V(NR) + PN02V(NR)) * 1.D20 * Sigma0 * Vtb

       !     *** Collision frequency (momentum transfer) ***

       ! Braginskii's collision time
       rNuee(NR) = PNeV(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AME) &
            &     * (ABS(PTeV(NR)) * rKeV)**1.5D0)
       rNuei(NR) = PNiV(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AME) &
            &     * (ABS(PTeV(NR)) * rKeV)**1.5D0)
       rNuie(NR) = PNeV(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AMI) &
            &     * (ABS(PTiV(NR)) * rKeV)**1.5D0)
       rNuii(NR) = PNiV(NR) * 1.D20 * PZ**4 * AEE**4 * rLnLam &
            &     / (12.D0 * PI * SQRT(PI) * EPS0**2 * SQRT(AMI) &
            &     * (ABS(PTiV(NR)) * rKeV)**1.5D0)
       ! Energy relaxation time between two kinds of particles
       rNuTei(NR) = PNiV(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
            &     / (3.D0 * SQRT(2.D0 * PI) * PI * EPS0**2 * AME * AMI &
            &     * (  ABS(MAX(PTeV(NR),PTeDIV)) * rKeV / AME &
            &        + ABS(PTiV(NR)) * rKeV / AMI)**1.5D0)

       !     *** Toroidal neoclassical viscosity ***
       !    (Hirshman and Sigmar, Nucl. Fusion 21 (1981) 1079)

       CALL TX_NCLASS(NR,rNueNC(NR),rNuiNC(NR),ETA3(NR),AJBS3(NR),IER)
       IF(IER /= 0) IERR = IER

       Wte = Vte / (Q(NR) * RR) ! Omega_te; transit frequency for electrons
       Wti = Vti / (Q(NR) * RR) ! Omega_ti; transit frequency for ions
       EpsL = R(NR) / RR        ! Inverse aspect ratio
       
       rNuAsE_inv = EpsL**1.5D0 * Wte / (SQRT(2.D0) * rNuei(NR))
       rNuAsI_inv = EpsL**1.5D0 * Wti / (SQRT(2.D0) * rNuii(NR))
       IF(NR /= 0) THEN
          rNuAse(NR) = 1.D0 / rNuAsE_inv
          rNuAsi(NR) = 1.D0 / rNuAsI_inv
       END IF
       BBL = SQRT(BphV(NR)**2 + BthV(NR)**2)
!!$       IF(R(NR) < RA) THEN
!!$          rNueNC(NR) = FSNC * (BphV(0) / BphV(NR))**2 / (2.D0*(1.D0-EpsL**2)**1.5D0) &
!!$               &     * SQRT(PI) * Q(NR)**2 * Wte &
!!$               &     * (1.53D0 / (rNuAsE_inv + SQRT(3.48D0 * rNuAsE_inv) + 1.52D0)) &
!!$               &     / (1.D0 + SQRT(0.37D0 * SQRT(2.D0) * rNuei(NR) / Wte) &
!!$               &                  + 1.25D0 * SQRT(2.D0) * rNuei(NR) / Wte)
!!$          rNuiNC(NR) = FSNC * (BphV(0) / BphV(NR))**2 / (2.D0*(1.D0-EpsL**2)**1.5D0) &
!!$               &     * SQRT(PI) * Q(NR)**2 * Wti &
!!$               &     * (0.53D0 / (rNuAsI_inv + SQRT(0.52D0 * rNuAsI_inv) + 0.56D0)) &
!!$               &     / (1.D0 + SQRT(0.14D0 * SQRT(2.D0) * rNuii(NR) / Wti) &
!!$               &                  + 0.70D0 * SQRT(2.D0) * rNuii(NR) / Wti)
!!$       ELSE
!!$          RL = (R(NR) - RA) / DBW
!!$          rNueNC(NR) = FSNC * (BphV(0) / BphV(NR))**2 / (2.D0*(1.D0-EpsL**2)**1.5D0) &
!!$               &     * SQRT(PI) * Q(NR)**2 * Wte &
!!$               &     * (1.53D0 / (rNuAsE_inv + SQRT(3.48D0 * rNuAsE_inv) + 1.52D0)) &
!!$               &     / (1.D0 + SQRT(0.37D0 * SQRT(2.D0) * rNuei(NR) / Wte) &
!!$               &                  + 1.25D0 * SQRT(2.D0) * rNuei(NR) / Wte) &
!!$               &     * 1.D0 / (1.D0 + RL**2)
!!$          rNuiNC(NR) = FSNC * (BphV(0) / BphV(NR))**2 / (2.D0*(1.D0-EpsL**2)**1.5D0) &
!!$               &     * SQRT(PI) * Q(NR)**2 * Wti &
!!$               &     * (0.53D0 / (rNuAsI_inv + SQRT(0.52D0 * rNuAsI_inv) + 0.56D0)) &
!!$               &     / (1.D0 + SQRT(0.14D0 * SQRT(2.D0) * rNuii(NR) / Wti) &
!!$               &                  + 0.70D0 * SQRT(2.D0) * rNuii(NR) / Wti) &
!!$               &     * 1.D0 / (1.D0 + RL**2)
!!$       END IF

!!$       IF(NR /= 0) THEN
!!$          NGRADB2 = EpsL**2 / (2.D0 * (1.D0 - EpsL**2)**1.5D0) &
!!$               &  * (BphV(0) / (RR * Q(NR)))**2
!!$          FTL  = 1.46D0 * SQRT(EpsL) - 0.46D0 * SQRT(EpsL) * EpsL
!!$          K11PSe = 0.5D0 * Ce * Vte**2 / rNuee(NR)
!!$          FCL = 1.D0 - FTL
!!$          K11Be  = BphV(0)**2 * FTL / FCL / (3.D0 * NGRADB2) &
!!$               & * (NUD(1.D0) * rNuee(NR) + NUD(Vte/Vti) * rNuei(NR))
!!$          K11Pe  = SQRT(PI) * Q(NR) * RR * Vte / 3.D0
!!$          K11PSi = 0.5D0 * Ci * Vti**2 / rNuii(NR)
!!$          K11Bi  = BphV(0)**2 * FTL / FCL / (3.D0 * NGRADB2) &
!!$               & * (NUD(1.D0) * rNuii(NR) + NUD(Vti/Vte) * rNuie(NR))
!!$          K11Pi  = SQRT(PI) * Q(NR) * RR * Vti / 3.D0
!!$          rNueNC(NR) = FSNC*3.D0/(2.D0*(1.D0-EpsL**2)**1.5D0)*(BphV(0)/BphV(NR)/RR)**2 &
!!$            &        *(K11Be * K11PSe) / (K11Be + K11PSe)
!!$          rNuiNC(NR) = FSNC*3.D0/(2.D0*(1.D0-EpsL**2)**1.5D0)*(BphV(0)/BphV(NR)/RR)**2 &
!!$            &        *(K11Bi * K11PSi) / (K11Bi + K11PSi)
!!$       END IF

!!$       rNueNC(NR) = FSNC * SQRT(PI) * Q(NR)**2 &
!!$            &     * Wte * 1.78D0 / (rNuAsE_inv + 1.78D0)
!!$       rNuiNC(NR) = FSNC * SQRT(PI) * Q(NR)**2 &
!!$            &     * Wti * 1.78D0 / (rNuAsI_inv + 1.78D0)
       !     &               * (1 + EpsL**1.5D0 * rNuAsI)
       !     &               / (1 + 1.44D0
       !     &                      * ((EpsL**1.5D0 * rNuAsI)**2
       !     &                         + ( ErV(NR)
       !     &                             / ( Vti * BthV(NR)) )**2))
       !!     &                         + ( ErV(NR) * BBL
       !!     &                             / ( Vti * BthV(NR)**2) )**2))

       !     *** Collision frequency (momentum transfer with beam) ***
       ! reference : 92/04/02, 92/04/21

       Vcr = (3.D0 * SQRT(PI / 2.D0) * PNiV(NR) * PZ**2 / PNeV(NR) * AME / AMI &
            &   * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)**(1.D0/3.D0)
       IF(PNBH == 0.D0 .AND. PNbV(NR) < 1.D-8) THEN
          rNube(NR) = 0.D0
          rNubi(NR) = 0.D0
          rNuB (NR) = 0.D0
       ELSE
          IF(PNbV(NR) /= 0.D0 .AND. NR < NRA) THEN
             RL = (R(NR) - RA) / DBW
             rNube(NR) = PNeV(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
                  &     / (3.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * AMB * AME &
                  &             * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)! &
!                  &     * RL**2 / (1.D0 + RL**2)
             rNubi(NR) = PNiV(NR) * 1.D20 * PZ**2 * PZ**2 * AEE**4 * rLnLam &
                  &     / (4.D0 * PI * EPS0**2 * AMB) &
                  &     * (1.D0 / AMB + 1.D0 / AMI) &
                  &     * 1.D0 / ( ABS(UbphV(NR))**3 + 3.D0 * SQRT(PI) / 4.D0 &
                  &     * Vti**1.5D0)! &
!                  &     * RL**2 / (1.D0 + RL**2)
             rNuB (NR) = rNube(NR) * 3.D0 / LOG(1.D0 + (Vb / Vcr)**3)
          ELSE
             rNube(NR) = 0.D0
             rNubi(NR) = 0.D0
             rNuB (NR) = 0.D0
          END IF
       END IF
!!$       IF(ABS(PNbV(NR)) < 1.D-6) THEN
!!$          rNube(NR) = 0.D0
!!$          rNubi(NR) = 0.D0
!!$          rNuB (NR) = 0.D0
!!$       ELSE
!!$          rNube(NR) = PNeV(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
!!$               &     / (3.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * AMB * AME &
!!$               &             * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)
!!$          rNubi(NR) = PNiV(NR) * 1.D20 * PZ**2 * PZ**2 * AEE**4 * rLnLam &
!!$               &     / (4.D0 * PI * EPS0**2 * AMB) &
!!$               &     * (1.D0 / AMB + 1.D0 / AMI) &
!!$               &     * 1.D0 / ( ABS(UbphV(NR))**3 + 3.D0 * SQRT(PI) / 4.D0 * Vti**1.5D0)
!!$          rNuB (NR) = rNube(NR) * 3.D0 / LOG(1.D0 + (Vb / Vcr)**3)
!!$       END IF

       !     *** Resistivity ***

       ! +++ Original model +++
       rNuBAR = rNuei(NR)+rNube(NR)*AMB*PNbV(NR)/(AME*PNeV(NR))+rNuL(NR)+rNu0e(NR)
       ALFA = (1.D0+rNueNC(NR)/rNuBAR)*(BthV(NR)/BphV(NR))**2
       ETA1(NR) = CORR(Zeff) * AME * (1.D0 + ALFA) * rNuBAR / (PNeV(NR)*1.D20 * AEE**2)

       ! +++ Hirshman, Hawryluk and Birge model +++
       ! Inverse aspect ratio
       EpsL    = R(NR) / RR
       ! Trapped particle fraction
       FTL  = 1.46D0 * SQRT(EpsL) - 0.46D0 * EpsL**1.5D0
       EFT  = FTL * rNuAsE_inv / (rNuAsE_inv + (0.58D0 + 0.20D0 * Zeff))
       ! Spitzer resistivity
       ETAS = CORR(Zeff) * AME * rNuei(NR) / (PNeV(NR) * 1.D20 * AEE**2)
       CR   = 0.56D0 * (3.D0 - Zeff) / ((3.D0 + Zeff) * Zeff)
       ETA2(NR) = ETAS * Zeff * (1.D0 + 0.27D0 * (Zeff - 1.D0)) &
            &   /((1.D0 - EFT) * (1.D0 - CR * EFT) * (1.D0 + 0.47D0 * (Zeff - 1.D0)))
       IF(MDLETA == 0) THEN
          ETA(NR) = ETA1(NR)
       ELSE
          ETA(NR) = ETA2(NR)
       END IF

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
       dQdr = 2.D0 * R(NR) * DERIV3(NR,PSI,Q,NRMAX,NRM,0)
       S(NR) = R(NR) / Q(NR) * dQdr
       dpdr = 2.D0 * R(NR) * DERIV3(NR,PSI,p,NRMAX,NRM,0)
       Alpha(NR) = - Q(NR)**2 * RR * dpdr * 2.D0 * rMU0 / (BphV(NR)**2 + BthV(NR)**2)

       !     *** Wave-particle interaction ***

       IF (ABS(FSCDBM) > 0.D0) THEN
          ! Alfven velocity
          Va = SQRT(BBL**2 / (rMU0 * PNiV(NR) * 1.D20 * AMI))
          ! Squared plasma frequency
          Wpe2 = PNeV(NR) * 1.D20 * AEE**2 / (AME * EPS0)
          ! Arbitrary coefficient for CDBM model
          rGC = 8.D0
          ! Magnetic curvature
          rKappa(NR) = - R(NR) / RR * (1.D0 - 1.D0 / Q(NR)**2)
          ! Calculate CDBM coefficient
          FCDBM(NR) = TRCOFS(S(NR),Alpha(NR),rKappa(NR))
          ! ExB rotational shear
          IF(NR == 0) THEN
             rH=0.D0
          ELSE
             dErdr = 2.D0 * R(NR) * DERIV3(NR,PSI,Vexbr,NRMAX,NRM,0)
             rH = Q(NR) * RR * R(NR) *  dErdr / (Va * S(NR))
          END IF
          
          ! Turbulence suppression by ExB shear
          rG1h2(NR) = 1.D0 / (1.D0 + rG1 * rH**2)
          ! Turbulent transport coefficient calculated by CDBM model
          DCDBM = rGC * FCDBM(NR) * rG1h2(NR) * ABS(Alpha(NR))**1.5D0 &
               &              * VC**2 / Wpe2 * Va / (Q(NR) * RR)
!          write(6,'(I3,3F15.7)') NR,ABS(Alpha(NR))**1.5D0,VC**2 / Wpe2 * Va / (Q(NR) * RR),Chie0*DCDBM
          !DCDBM = MAX(DCDBM,1.D-05)
       ELSE
          rG1h2(NR)  = 0.D0
          FCDBM(NR)  = 0.D0
          DCDBM      = 0.D0
       END IF

       DeL = 1.D0
!!$       DeL = FSDFIX * (1.D0 + (PROFDL -1.D0) * (R(NR) / RA)**2) + FSCDBM * DCDBM
       ! Particle diffusivity
       De(NR)   = De0   * DeL
       Di(NR)   = Di0   * DeL
       IF (R(NR) < RA) THEN
          DeL = FSDFIX * (1.D0 + (PROFDL -1.D0) * (R(NR) / RA)**2) + FSCDBM * DCDBM
       ELSE
          IF(NR == NRA) DeLa = FSDFIX * PROFDL + FSCDBM * DCDBM
          IF(FSPSCL == 0.D0) THEN
             DeL =  (1.D0 - MOD(FSBOHM,2.D0)) * FSDFIX * PROFDL &
                  &+ FSBOHM * PTeV(NR) * rKeV / (16.D0 * AEE * BBL)
          ELSE
             DeL = FSPSCL
          END IF
       END IF
       ! Viscosity
       rMue(NR) = rMue0 * DeL
       rMui(NR) = rMui0 * DeL
       ! Thermal diffusivity
       Chie(NR) = Chie0 * DeL
       Chii(NR) = Chii0 * DeL

       ! <omega/m>
       WPM(NR) = WPM0 * PTeV(NR) * rKeV / (RA**2 * AEE * BphV(NR))
       ! Force induced by drift wave (e.q.(8),(13))
       FWthe(NR)   = AEE**2         * BphV(NR)**2 * De(NR) / (PTeV(NR) * rKeV)
       FWthi(NR)   = AEE**2 * PZ**2 * BphV(NR)**2 * Di(NR) / (PTiV(NR) * rKeV)
       FWthphe(NR) = AEE**2         * BphV(NR)    * De(NR) / (PTeV(NR) * rKeV)
       FWthphi(NR) = AEE**2 * PZ**2 * BphV(NR)    * Di(NR) / (PTiV(NR) * rKeV)
!!$       IF(NR == NRA) THEN
!!$          FWthea = AEE**2         * BphV(NR)**2 * De0 * DeLa / (PTeV(NR) * rKeV)
!!$          FWthia = AEE**2 * PZ**2 * BphV(NR)**2 * Di0 * DeLa / (PTiV(NR) * rKeV)
!!$       END IF
!!$       IF(R(NR) < RA) THEN
!!$          FWthi(NR) = AEE**2 * PZ**2 * BphV(NR)**2 * Di0 * R(NR) &
!!$               &            / (PTiV(NR) * rKeV)
!!$       ELSE
!!$          FWthi(NR) = AEE**2 * PZ**2 * BphV(NR)**2 * Di0 * RA * 2.d0 &
!!$               &            / (PTiV(NR) * rKeV)
!!$       END IF
!!$       IF(NR == NRA) THEN
!!$          FWthia = AEE**2 * PZ**2 * BphV(NR)**2 * Di0 * RA * 2.d0 / (PTiV(NR) * rKeV)
!!$       END IF

       ! Work induced by drift wave
       IF(MDLWTB == 1) THEN
          WNthe (NR) = WPE0 *         AEE    * BphV(NR)    * De(NR)
          WEMthe(NR) = WPE0 *         AEE**2 * BphV(NR)    * De(NR) &
               &                       / (PTeV(NR) * rKeV)
          WWthe (NR) = WPE0 *         AEE**2 * BphV(NR)**2 * De(NR) &
               &            * WPM (NR) / (PTeV(NR) * rKeV)
          WT1the(NR) = WPE0 *         AEE    * BphV(NR)    *(rMue(NR) - 0.5D0 * De(NR)) &
               &                       / (PTeV(NR) * rKeV)
          WT2the(NR) = WPE0 *         AEE    * BphV(NR)    *(rMue(NR) - 0.5D0 * De(NR))
          WNthi (NR) = WPI0 * PZ    * AEE    * BphV(NR)    * Di(NR)
          WEMthi(NR) = WPI0 * PZ**2 * AEE**2 * BphV(NR)    * Di(NR) &
               &                       / (PTiV(NR) * rKeV)
          WWthi (NR) = WPI0 * PZ**2 * AEE**2 * BphV(NR)**2 * Di(NR) &
               &            * WPM (NR) / (PTiV(NR) * rKeV)
          WT1thi(NR) = WPI0 * PZ    * AEE    * BphV(NR)    *(rMui(NR) - 0.5D0 * Di(NR)) &
               &            / (PTiV(NR) * rKeV)
          WT2thi(NR) = WPI0 * PZ    * AEE    * BphV(NR)    *(rMui(NR) - 0.5D0 * Di(NR))
       END IF

       !     *** Heating profile ***

       IF (R(NR) < RA) THEN
          PNB(NR) = PNB0  * EXP(- R(NR)**2 / RNB**2) * (1.D0 - (R(NR) / RA)** 4)
          SNB(NR) = PNB(NR) / (Eb * rKeV * 1.D20)
          PRFe(NR)= PRFe0 * EXP(- R(NR)**2 / RRF**2) * (1.D0 - (R(NR) / RA)** 4)
          PRFi(NR)= PRFi0 * EXP(- R(NR)**2 / RRF**2) * (1.D0 - (R(NR) / RA)** 4)
       ELSE
          PNB(NR) =0.D0
          SNB(NR) =0.D0
          PRFe(NR)=0.D0
          PRFi(NR)=0.D0
       END IF

       !     *** Current density profile ***

       ! Poloidal current density
       AJPH  = -      AEE * PNeV(NR) * 1.D20 * UephV(NR) &
            &  + PZ * AEE * PNiV(NR) * 1.D20 * UiphV(NR) &
            &  + PZ * AEE * PNbV(NR) * 1.D20 * UbphV(NR)
       ! Toroidal current density
       AJTH  = -      AEE * PNeV(NR) * 1.D20 * UethV(NR) &
            &  + PZ * AEE * PNiV(NR) * 1.D20 * UithV(NR) &
            &  + PZ * AEE * PNbV(NR) * 1.D20 * UbthV(NR)

       BBL=SQRT(BthV(NR)**2+BphV(NR)**2)
       ! Parallel current density
       AJPARA=(BthV(NR)*AJTH     + BphV(NR)*AJPH    )/BBL
       ! Parallel electric field
       EPARA =(BthV(NR)*EthV(NR) + BphV(NR)*EphV(NR))/BBL
       ! Total current density = parallel current density(?)
!       AJ(NR)   = AJPARA
       AJ(NR)   = AJPH
       ! Ohmic current density
       AJOH(NR) = EphV(NR) / ETA(NR)
!       AJOH(NR) = EPARA / ETA(NR)
       ! Ohmic heating power
!       POH(NR)  = EPARA*AJPARA
!       POH(NR)  = EthV(NR)*AJTH + EphV(NR)*AJPH
       POH(NR)  = EphV(NR) * AJPH
       ! NB induced current density
       AJNB(NR) = PZ * AEE * PNbV(NR) * 1.D20 * UbphV(NR)

       !     *** Bremsstraulung loss ***

       PBr(NR) = 5.35D-37 * PZ**2 * PNeV(NR) * PNiV(NR) * 1.D40 * SQRT(PTeV(NR))

       !     *** Loss to divertor ***

!       IF (R(NR) + DBW > RA) THEN
       IF (R(NR) > RA) THEN
          Cs = SQRT(2.D0 * PTeV(NR) * rKeV / AMI)
          RL = (R(NR) - RA) / DBW! / 2.D0
          rNuL  (NR) = FSLP  * Cs / (2.D0 * PI * Q(NR) * RR) &
               &             * RL**2 / (1.D0 + RL**2)
          KAPPA = (4.D0*PI*EPS0)**2/(SQRT(AME)*rLnLam*AEE**4*Zeff)*AEE**2.5D0
          rNuLTe(NR) = FSLTE * KAPPA * (PTeV_FIX(NR)*1.D3)**2.5D0 &
               &                  /((2.D0 * PI * Q(NR) * RR)**2 * PNeV_FIX(NR)*1.D20) &
               &             * RL**2 / (1.D0 + RL**2)
!!$          rNuLTe(NR) = FSLTE * Cs / (2.D0 * PI * Q(NR) * RR) &
!!$               &             * RL**2 / (1.D0 + RL**2)
          rNuLTi(NR) = FSLTI * Cs / (2.D0 * PI * Q(NR) * RR) &
               &             * RL**2 / (1.D0 + RL**2)
!          write(6,*) rNuL(NR), rNuLTe(NR), rNuLTi(NR)
!          write(6,*) (4.D0*PI*EPS0)**2/(SQRT(AME)*rLnLam*AEE**4*PZ)*AEE**2.5D0
       ELSE
          rNuL(NR) = 0.D0
          rNuLTe(NR) = 0.D0
          rNuLTi(NR) = 0.D0
       END IF

    END DO L_NR

    rNuAse(0) = AITKEN2P(PSI(0),rNuAse(1),rNuAse(2),rNuAse(3),PSI(1),PSI(2),PSI(3))
    rNuAsi(0) = AITKEN2P(PSI(0),rNuAsi(1),rNuAsi(2),rNuAsi(3),PSI(1),PSI(2),PSI(3))
!    rNueNC(0) = AITKEN2P(PSI(0),rNueNC(1),rNueNC(2),rNueNC(3),PSI(1),PSI(2),PSI(3))
!    rNuiNC(0) = AITKEN2P(PSI(0),rNuiNC(1),rNuiNC(2),rNuiNC(3),PSI(1),PSI(2),PSI(3))
    rNueNC(0) = rNueNC(1)
    rNuiNC(0) = rNuiNC(1)

!###########################
!!$    DO NR = 0, NRA
!!$       RL = (R(NR) - RA) / DBW
!!$       rNueNC(NR) = FSNC * 5.D6 * RL**2 / (1.D0 + RL**2)
!!$       rNuiNC(NR) = FSNC * 1.5D5 * RL**2 / (1.D0 + RL**2)
!!$    END DO

!!$    DO NR=0,NRMAX
!!$       rNueNC(NR) = 20.d6
!!$       rNuiNC(NR) = 1.D6
!!$    END DO
!###########################

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

  REAL(8) FUNCTION CORR(X)
    ! X is the effective charge number
    real(8), intent(in) :: X

    CORR = (1.D0 + 1.198D0 * X + 0.222D0 * X**2) &
    &    / (1.D0 + 2.966D0 * X + 0.753D0 * X**2)

  END FUNCTION CORR

  REAL(8) FUNCTION NUD(X)
    real(8), intent(in) :: X

    NUD = SQRT(1.D0 + X**2) + X**2 * LOG(X / (1.D0 + SQRT(1.D0 + X**2)))

  END FUNCTION NUD

end module variables
