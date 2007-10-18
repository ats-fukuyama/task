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

    use physical_constants, only : rMU0, rKeV
    use libraries, only : INTG_P, DERIVF, VALINT_SUB
    REAL(8), DIMENSION(1:NQM,0:NRMAX), INTENT(INOUT) :: XL
    integer, intent(in), optional :: ID
    INTEGER :: NR
    real(8) :: DERIV3, FCTR

    Phi  (0:NRMAX) =   XL(LQm1,0:NRMAX)
    DO NR = 0, NRMAX
       ErV (NR) = - 2.D0 * R(NR) * DERIVF(NR,PSI,XL,LQm1,NQMAX,NRMAX)
    END DO
    EthV (0)       =   0.D0
    EthV (1:NRMAX) = - XL(LQm2,1:NRMAX) / R(1:NRMAX)
    EphV (0:NRMAX) = - XL(LQm3,0:NRMAX)
    AphV (0:NRMAX) =   XL(LQm4,0:NRMAX) * rMUb2
    RAthV(0:NRMAX) =   XL(LQm5,0:NRMAX) * rMU0
    DO NR = 0, NRMAX
       BthV(NR) = - 2.D0 * R(NR) * DERIVF(NR,PSI,XL,LQm4,NQMAX,NRMAX) * rMUb2
       BphV(NR) =   2.D0         * DERIVF(NR,PSI,XL,LQm5,NQMAX,NRMAX) * rMU0
    END DO
    PNeV (0:NRMAX) =   XL(LQe1,0:NRMAX)
    UerV (0)       =   0.D0
    UerV (1:NRMAX) =   XL(LQe2,1:NRMAX )/ PNeV(1:NRMAX) / R(1:NRMAX)
    RUethV(0:NRMAX)=   XL(LQe3,0:NRMAX) / PNeV(0:NRMAX)
    UethV(0)       =   0.D0
    UethV(1:NRMAX) =   XL(LQe3,1:NRMAX) / PNeV(1:NRMAX) / R(1:NRMAX)
    UephV(0:NRMAX) =   XL(LQe4,0:NRMAX) / PNeV(0:NRMAX) * AMPe4
    IF(MDFIXT == 0) THEN
       PeV  (0:NRMAX) =   XL(LQe5,0:NRMAX)
       PTeV (0:NRMAX) =   XL(LQe5,0:NRMAX) / PNeV(0:NRMAX)
    ELSE
       PeV  (0:NRMAX) =   XL(LQe5,0:NRMAX) * PNeV(0:NRMAX)
       PTeV (0:NRMAX) =   XL(LQe5,0:NRMAX)
    END IF
    PNiV (0:NRMAX) =   XL(LQi1,0:NRMAX)
    UirV (0)       =   0.D0
    UirV (1:NRMAX) =   XL(LQi2,1:NRMAX) / PNiV(1:NRMAX) / R(1:NRMAX)
    RUithV(0:NRMAX)=   XL(LQi3,0:NRMAX) / PNiV(0:NRMAX)
    UithV(0)       =   0.D0
    UithV(1:NRMAX) =   XL(LQi3,1:NRMAX) / PNiV(1:NRMAX) / R(1:NRMAX)
    UiphV(0:NRMAX) =   XL(LQi4,0:NRMAX) / PNiV(0:NRMAX)
    IF(MDFIXT == 0) THEN
       PiV  (0:NRMAX) =   XL(LQi5,0:NRMAX)
       PTiV (0:NRMAX) =   XL(LQi5,0:NRMAX) / PNiV(0:NRMAX)
    ELSE
       PiV  (0:NRMAX) =   XL(LQi5,0:NRMAX) * PNiV(0:NRMAX)
       PTiV (0:NRMAX) =   XL(LQi5,0:NRMAX)
    END IF
    PNbV (0:NRMAX) =   XL(LQb1,0:NRMAX)
    IF(ABS(FSRP) > 0.D0) THEN
       DO NR = 0, NRMAX
          IF(PNbV(NR) == 0.D0) THEN ! The region without beam particles
             UbthV(NR) = 0.D0
             UbphV(NR) = 0.D0
          ELSE ! The region with beam particles
             IF(NR == 0) THEN ! On axis, poloidal beam velocity is assumed to be zero.
                UbthV(NR) = 0.D0
                UbphV(NR) = XL(LQb4,NR) / PNbV(NR)
             ELSE ! The region except the magnetic axis
                UbthV(NR) = XL(LQb3,NR) / PNbV(NR) / R(NR)
                UbphV(NR) = XL(LQb4,NR) / PNbV(NR)
             END IF
          END IF
       END DO
    ELSE
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
       UbthV(NRA:NRMAX) = 0.D0
       UbphV(NRA:NRMAX) = 0.D0
    END IF

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

    PNbRPV(0:NRMAX)=   XL(LQr1,0:NRMAX)

    PT01V(0:NRMAX) =   0.5D0 * AMI * V0**2 / rKeV
    PT02V(0:NRMAX) =   PTiV(0:NRMAX)

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

    use physical_constants, only : AEE, AME, VC, PI, rMU0, EPS0, rKeV, EION
    use libraries, only : EXPV, VALINT_SUB, TRCOFS, INTG_P
    use core_module, only : inv_int
    use nclass_mod
    use sauter_mod

    INTEGER :: NR, NP, NR1, IER, i, imax, nrl, test
    REAL(8) :: Sigma0, QL, SL, SLP1, SLP2, PNBP0, PNBT10, PNBT20, PRFe0, PRFi0, &
         &     Vte, Vti, Vtb, XXX, SiV, ScxV, Wte, Wti, EpsL, rNuPara, &
         &     rNuAsE_inv, rNuAsI_inv, BBL, Va, Wpe2, rGC, SP, rGBM, &
         &     rGIC, rH, dErdr, dpdr, PROFDL, PROFDDL, &
         &     DCDBM, DeL, AJPH, AJTH, AJPARA, EPARA, Vcr, &
         &     Cs, RhoIT, ExpArg, AiP, DISTAN, UbparaL, &
         &     SiLCL, SiLCthL, SiLCphL, Wbane, Wbani, RL, ALFA, DBW, PTiVA, &
         &     KAPPA, rNuBAR, Ecr, factor_bohm, rNuAsIL, VAL1, VAL2, VAL, &
         &     rhob, rNueff, rNubnc, DCB, DRP, Dltcr, DltR2, &
         &     theta1, theta2, dlt, width0, width1, ARC, &
         &     EbL, logEbL, Scx, Vave, Sion, Left, Right, RV0, tmp, RLOSS !&
!!         &     NGRADB2, K11PSe, K11Be,  K11Pe, K11PSi, K11Bi, K11Pi
    real(8) :: rnubarth, rnubarph
!!    real(8) :: Ce = 0.733D0, Ci = 1.365D0
    real(8) :: FCL, EFT, CR, dPTeV, dPTiV, dPPe, dPPi, SUML
    real(8) :: DERIV3, AITKEN2P
    real(8), dimension(0:NRMAX) :: p, Vexbr, dQdr, SP0, SP1, SP2, SN0, SN1, SN2, &
         &                         SP3, th1, th2, Ubpara!,PNbrpL, DERIV
!    real(8), dimension(1:4,0:NRMAX) :: U

    !     *** Constants ***

    !     Neutral cross section
    !     (NRL Plasma Formulary p52 Eq. (1) (2002))

    Sigma0 = 8.8D-21

    !     NBI beam speed

    Vb =  SQRT(2.D0 * Eb * rKeV / AMB)

    !     Poloidal magnetic field on wall

    IF(FSHL == 0.D0) THEN
       Bthb = rMU0 * rIP * 1.D6 / (2.D0 * PI * RB)
    ELSE
       QL=(Q0-QA)*(1.D0-(RB/RA)**2)+QA
       Bthb = BB*RB/(QL*RR)
    END IF

    !   NBI total input power (MW)
    PNBH = PNBHP + PNBHT1 + PNBHT2

    !     *** Normalization factor for heating profile ***
    !
    !    SL is a normalization factor for a given heating profile.
    !    It is assumed that the heating profile has a shape of 
    !    exp(-r^2/r_NB^2)*(1-(r/a)^4), in order to renormalize the
    !    heating profile we therefore need to integrate that profile with 
    !    maximum value of unity at axis in advance and calculate the 
    !    normalized factor (i.e. PNBP0, PNBT0) which allows an integration
    !    value of renormalized profile to equal designated value (i.e. PNBH).
    !    The factor (1-(r/a)^4), which is an arbitrary function, is
    !    required to damp the profile of exp(-r^2/r_NB^2) rapidly
    !    because if this factor is not imposed unrealistically large current
    !    driven by NBI is generated near the edge where electron and bulk ion
    !    density are dilute although in this case fast ion density from 
    !    NBI is relatively large.

    !  For NBI heating
    !  *** Perpendicular
    IF(ABS(FSRP) > 0.D0) THEN
       SP0(0:NRMAX) = EXP(- ((R(0:NRMAX) - RNBP0) / RNBP)**2) * (1.D0 - (R(0:NRMAX) / RB)** 4)
       CALL VALINT_SUB(SP0,NRMAX,SL)
       SN0(0:NRMAX) = SP0(0:NRMAX)
    ELSE
       SP0(0:NRA) = EXP(- ((R(0:NRA) - RNBP0) / RNBP)**2) * (1.D0 - (R(0:NRA) / RA)** 4)
       SP0(NRA+1:NRMAX) = 0.D0
       CALL VALINT_SUB(SP0,NRA,SL)
       SN0(0:NRA) = SP0(0:NRA)
       SN0(NRA+1:NRMAX) = 0.D0
    END IF
    SL = 2.D0 * PI * SL

    PNBP0 = ABS(PNBHP) * 1.D6 / (2.D0 * Pi * RR * SL)

    !  *** Tangential
    IF(ABS(FSRP) > 0.D0) THEN
       SP1(0:NRMAX) = EXP(- ((R(0:NRMAX) - RNBT10) / RNBT1)**2) * (1.D0 - (R(0:NRMAX) / RB)** 4)
       CALL VALINT_SUB(SP1,NRMAX,SL)
       SN1(0:NRMAX) = SP1(0:NRMAX)
    ELSE
       SP1(0:NRA) = EXP(- ((R(0:NRA) - RNBT10) / RNBT1)**2) * (1.D0 - (R(0:NRA) / RA)** 4)
       SP1(NRA+1:NRMAX) = 0.D0
       CALL VALINT_SUB(SP1,NRA,SL)
       SN1(0:NRA) = SP1(0:NRA)
       SN1(NRA+1:NRMAX) = 0.D0
    END IF
    SLP1 = 2.D0 * PI * SL

    IF(ABS(FSRP) > 0.D0) THEN
       SP2(0:NRMAX) = EXP(- ((R(0:NRMAX) - RNBT20) / RNBT2)**2) * (1.D0 - (R(0:NRMAX) / RB)** 2)
       CALL VALINT_SUB(SP2,NRMAX,SL)
       SN2(0:NRMAX) = SP2(0:NRMAX)
    ELSE
       SP2(0:NRA) = EXP(- ((R(0:NRA) - RNBT20) / RNBT2)**2) * (1.D0 - (R(0:NRA) / RA)** 2)
       SP2(NRA+1:NRMAX) = 0.D0
       CALL VALINT_SUB(SP2,NRA,SL)
       SN2(0:NRA) = SP2(0:NRA)
       SN2(NRA+1:NRMAX) = 0.D0
    END IF
    SLP2 = 2.D0 * PI * SL

    PNBT10 = ABS(PNBHT1) * 1.D6 / (2.D0 * Pi * RR * SLP1)
    PNBT20 = ABS(PNBHT2) * 1.D6 / (2.D0 * Pi * RR * SLP2)

    !  For RF heating
!!$    SP3(0:NRMAX) = 0.D0
!!$    SP3(0:NRA) = EXP(- ((R(0:NRA) - RRF0) / RRF)**2) * (1.D0 - (R(0:NRA) / RA)** 4)
!!$    CALL VALINT_SUB(SP3,NRA,SL)
    SP3(0:NRMAX) = EXP(- ((R(0:NRMAX) - RRF0) / RRF)**2) * (1.D0 - (R(0:NRMAX) / RB)** 4)
    CALL VALINT_SUB(SP3,NRMAX,SL)
    SL = 2.D0 * PI * SL

    PRFe0 = 0.5D0 * PRFH * 1.D6 / (2.D0 * Pi * RR * SL)
    PRFi0 = 0.5D0 * PRFH * 1.D6 / (2.D0 * Pi * RR * SL)

    p(0:NRMAX) = (PeV(0:NRMAX) + PiV(0:NRMAX)) * 1.D20 * rKeV
    Vexbr(1:NRMAX) = ErV(1:NRMAX) &
         &         / (R(1:NRMAX) * SQRT(BphV(1:NRMAX)**2 + BthV(1:NRMAX)**2))

    IF(PROFD == 0.D0.AND.FSDFIX /= 0.D0) THEN
       PROFDL = (PTeV(NRA) * rKeV / (16.D0 * AEE * &
            &    SQRT(BphV(NRA)**2 + BthV(NRA)**2))) / FSDFIX
    ELSE
       PROFDL = PROFD
    END IF

    ! Banana width
!    Wbane = (Q(0) * SQRT(RR * AME * PTeV(0) * rKeV) / (AEE * BphV(0)))**(2.D0/3.D0)
!    Wbani = (Q(0) * SQRT(RR * AMI * PTiV(0) * rKeV) / (PZ * AEE * BphV(0)))**(2.D0/3.D0)

    ! Banana width at separatrix
!    PTiVA = PTiV(NRA)
    PTiVA = 0.5D0 * PTiV(0)
    DBW = 3.D0 * SQRT(PTiVA * rKEV * AMI) * Q(NRA) / (PZ * AEE * BphV(NRA)) &
         & / SQRT(R(NRA) / RR)

    ! Ripple amplitude
    DltRP(0:NRMAX) = DltRP0 * (  ((RR + R(0:NRMAX)) / (RR + RA))**(NTCOIL-1) &
         &                     + ((RR - RA) / (RR + R(0:NRMAX)))**(NTCOIL+1) * DIN)

    PNbrpLV(0:NRMAX) = 0.D0

    EbL = Eb * 1.D3 / PA
    logEbL = log10(EbL)

    !  Coefficients

!    write(6,*) "++++++++++++++++++++++"
    L_NR: DO NR = 0, NRMAX

       Vte = SQRT(2.D0 * ABS(PTeV(NR)) * rKeV / AME)
       Vti = SQRT(2.D0 * ABS(PTiV(NR)) * rKeV / AMI)
       Vtb = SQRT(2.D0 * ABS(PTiV(NR)) * rKeV / AMB)

       !     *** Coulomb logarithms ***
       !     (NRL Plasma Formulary p34,35 (2002))
       !     rlnLee = rlnLei = rlnLie 
       !            = 24 + 6*ln10 - ln(sqrt(PNeV [10^20 m^-3] / PTeV [keV]))
       !     rlnLii = 23 + 15/2*ln10 - ln(PZ**2/PTiV*sqrt(2*PNiV*PZ**2/PTiV))

       rlnLe(NR) = 37.8d0 - LOG(SQRT(PNeV(NR)*1.D20)/(PTeV(NR)))
       rlnLi(NR) = 40.3d0 - LOG(PZ**2/PTiV(NR)*SQRT(2.D0*PNiV(NR)*1.D20*PZ**2/PTiV(NR)))

       !     *** Ionization rate ***
       !     (NRL Plasma Formulary p54 Eq. (12) (2002))
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

       !     *** Charge exchange rate ***
!!$       !     (Amano and Okamoto, JAERI-M 8420)
!!$
!!$       XXX = LOG10(MAX(PTiV(NR) * 1.D3, 50.D0))
!!$       ScxV = 1.57D-16 * SQRT(PTiV(NR) * 1.D3) &
!!$            &          * (XXX * XXX - 14.63D0 * XXX + 53.65D0)
!!$       rNuiCX(NR) = FSCX * ScxV * (PN01V(NR) + PN02V(NR)) * 1.D20

       !     (Riviere, NF 11 (1971) 363)
       !  For thermal ions
       Scx = 6.937D-19 * (1.D0 - 0.155D0 * LOG10(PTiV(NR)*1.D3/PA))**2
       Vave = SQRT(8.D0 * PTiV(NR) * rKeV / (PI * AMI))
       rNuiCX(NR) = FSCX * Scx * Vave * (PN01V(NR) + PN02V(NR)) * 1.D20

       !  For beam ions
       Scx = 6.937D-19 * (1.D0 - 0.155D0 * logEbL)**2 &
            & / (1.D0 + 1.112D-15 * EbL**3.3D0)
       Vave = SQRT(8.D0 * Eb * rKeV / (PI * AMI))
       rNubCX(NR) = FSCX * Scx * Vave * (PN01V(NR) + PN02V(NR)) * 1.D20

       !     *** Collision frequency (with neutral) ***

       rNu0e(NR) = (PN01V(NR) + PN02V(NR)) * 1.D20 * Sigma0 * Vte
       rNu0i(NR) = (PN01V(NR) + PN02V(NR)) * 1.D20 * Sigma0 * Vti
       rNu0b(NR) = (PN01V(NR) + PN02V(NR)) * 1.D20 * Sigma0 * Vtb

       !     *** Trapped particle fraction ***
       !     (Y. B. Kim, et al., Phys. Fluids B 3 (1990) 2050)
       EpsL = R(NR) / RR        ! Inverse aspect ratio
       ft(NR) = 1.46D0 * SQRT(EpsL) - 0.46 * EpsL * SQRT(EpsL)

       !     *** Collision frequency (momentum transfer) ***

       ! Braginskii's collision time
       rNuee(NR) = PNeV(NR) * 1.D20 *         AEE**4 * rlnLe(NR) &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AME) &
            &     * (ABS(PTeV(NR)) * rKeV)**1.5D0)
       rNuei(NR) = PNiV(NR) * 1.D20 * PZ**2 * AEE**4 * rlnLe(NR) &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AME) &
            &     * (ABS(PTeV(NR)) * rKeV)**1.5D0)
       rNuie(NR) = (AME / AMI) * rNuei(NR)
       ! Caution!: rNuii is not the same as rNui, rNui = sqrt(2) * rNuii
       rNuii(NR) = PNiV(NR) * 1.D20 * PZ**4 * AEE**4 * rlnLi(NR) &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AMI) &
            &     * (ABS(PTiV(NR)) * rKeV)**1.5D0)
       ! Energy relaxation time between two kinds of particles
!       rNuTei(NR) = PNiV(NR) * 1.D20 * PZ**2 * AEE**4 * rlnLe(NR) &
!            &     / (3.D0 * SQRT(2.D0 * PI) * PI * EPS0**2 * AME * AMI &
!            &     * (  ABS(MAX(PTeV(NR),PTeDIV)) * rKeV / AME &
!            &        + ABS(PTiV(NR)) * rKeV / AMI)**1.5D0)
       rNuTei(NR) = rNuei(NR) * (2.D0 * AME / AMI)

       BBL = SQRT(BphV(NR)**2 + BthV(NR)**2)

       rNuPara = CORR(Zeff) * rNuei(NR)
       rNuei1(NR)  =(BthV(NR)**2 * rNuPara + BphV(NR)**2 * rNuei(NR)) / BBL**2
       rNuei2(NR)  = BphV(NR) / BBL**2 * (rNuPara - rNuei(NR))
       rNuei3(NR)  =(BphV(NR)**2 * rNuPara + BthV(NR)**2 * rNuei(NR)) / BBL**2

       !     *** Toroidal neoclassical viscosity ***
       !     (W. A. Houlberg, et al., Phys. Plasmas 4 (1997) 3230)

       CALL TX_NCLASS(NR,rNueNC(NR),rNuiNC(NR),ETA2(NR),AJBS2(NR),IER)
       IF(IER /= 0) IERR = IER

       Wte = Vte / (Q(NR) * RR) ! Omega_te; transit frequency for electrons
       Wti = Vti / (Q(NR) * RR) ! Omega_ti; transit frequency for ions
       
       rNuAsE_inv = EpsL**1.5D0 * Wte / (SQRT(2.D0) * rNuei(NR))
       rNuAsI_inv = EpsL**1.5D0 * Wti / (SQRT(2.D0) * rNuii(NR))
       IF(NR /= 0) THEN
          rNuAse(NR) = 1.D0 / rNuAsE_inv
          rNuAsi(NR) = 1.D0 / rNuAsI_inv
       END IF
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

       !     *** Beam slowing down time (momentum transfer with beam) ***
       ! reference : memo (92/04/02, 92/04/21)
       !             Tokamaks 3rd pp.246 - 252

!!!       Vcr = (3.D0 * SQRT(PI / 2.D0) * PNiV(NR) * PZ**2 / PNeV(NR) * AME / AMI &
!!!            &   * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)**(1.D0/3.D0)
!       Ecr = (9.D0 * PI / 16.D0 * AMI / AME)**(1.D0/3.D0) * AMB / AMI * PTeV(NR) ! in keV
       Ecr = 14.8D0 * (PA / PA**(2.D0/3.D0)) * PTeV(NR) ! in keV
       PNBcol_i(NR) = NBIi_ratio(Eb/Ecr)
       PNBcol_e(NR) = 1.d0 - PNBcol_i(NR)
       IF(PNBH == 0.D0 .AND. PNbV(NR) < 1.D-8) THEN
          rNube(NR) = 0.D0
          rNubi(NR) = 0.D0
          rNuB (NR) = 0.D0
       ELSE
          rNube(NR) = PNeV(NR) * 1.D20 * PZ**2 * AEE**4 * rlnLe(NR) &
               &     / (3.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * AMB * AME &
               &             * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)
          rNubi(NR) = PNiV(NR) * 1.D20 * PZ**2 * PZ**2 * AEE**4 * rlnLi(NR) &
               &     / (4.D0 * PI * EPS0**2 * AMB) &
               &     * (1.D0 / AMB + 1.D0 / AMI) &
               &     * 1.D0 / (Vb**3 + 0.75D0 * SQRT(PI) * Vti**3)

          rNuPara = CORR(Zeff) * rNube(NR)
          rNube1(NR)  =(BthV(NR)**2 * rNuPara + BphV(NR)**2 * rNube(NR)) / BBL**2
          rNube2(NR)  = BphV(NR) / BBL**2 * (rNuPara - rNube(NR))
          rNube3(NR)  =(BphV(NR)**2 * rNuPara + BthV(NR)**2 * rNube(NR)) / BBL**2

          ! deflection time of beam ions against bulk ions
          ! (Takamura (3.26), Tokamaks 3rd p64)
          ! ** rNuD is similar to rNubi. **
          rNuD(NR) = PNiV(NR) *1.D20 * PZ**2 * PZ**2 * AEE**4 * rlnLi(NR) &
               &   / (2.D0 * PI * EPS0**2 * AMB**2 * Vb**3) &
               &   * (0.5D0 * Vb * (  2.D0 / (Vb + 0.5D0 * SQRT(PI) * Vti) &
               &                    - Vti**2 / (Vb**3 + 0.75D0 * SQRT(PI) * Vti**3)))
          ! effective time of detrapping (Takamura (5.31))
          IF(DltRP(NR) /= 0.D0) THEN
             rNubrp1(NR) = FSRP * rNuD(NR) / DltRP(NR)
             rNubrp2(NR) = FSRP * rNuD(NR) / SQRT(DltRP(NR))! * SQRT(2.D0)
          ELSE
             rNubrp1(NR) = 0.D0
             rNubrp2(NR) = 0.D0
          END IF

          ! The definition given below is a "energy slowing down time" and is not
          ! "particle slowing down time". At present, we assume the former is the
          ! same as the latter.
          rNuB (NR) = rNube(NR) * 1.5D0 / LOG(1.D0 + (Eb / Ecr)**1.5D0)
!!!          rNuB (NR) = rNube(NR) * 3.D0 / LOG(1.D0 + (Vb / Vcr)**3)
!!          rNuB (NR) = rNube(NR) + rNubi(NR)
       END IF
!!$       IF(ABS(PNbV(NR)) < 1.D-6) THEN
!!$          rNube(NR) = 0.D0
!!$          rNubi(NR) = 0.D0
!!$          rNuB (NR) = 0.D0
!!$       ELSE
!!$          rNube(NR) = PNeV(NR) * 1.D20 * PZ**2 * AEE**4 * rlnLe(NR) &
!!$               &     / (3.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * AMB * AME &
!!$               &             * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)
!!$          rNubi(NR) = PNiV(NR) * 1.D20 * PZ**2 * PZ**2 * AEE**4 * rlnLi(NR) &
!!$               &     / (4.D0 * PI * EPS0**2 * AMB) &
!!$               &     * (1.D0 / AMB + 1.D0 / AMI) &
!!$               &     * 1.D0 / ( ABS(UbphV(NR))**3 + 3.D0 * SQRT(PI) / 4.D0 * Vti**1.5D0)
!!$          rNuB (NR) = rNube(NR) * 3.D0 / LOG(1.D0 + (Vb / Vcr)**3)
!!$       END IF

       !     *** Helical neoclassical viscosity ***

       IF(ABS(FSHL) > 0.D0) THEN
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
       END IF

       !  Derivatives (beta, safety factor, mock ExB velocity)
       dQdr(NR) = 2.D0 * R(NR) * DERIV3(NR,PSI,Q,NRMAX,NRM,0)
       S(NR) = R(NR) / Q(NR) * dQdr(NR)
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

!!$       IF (R(NR) < RA) THEN
!!$          DeL = FSDFIX * (1.D0 + 4.D0 * (R(NR) / RA)**2) + FSCDBM * DCDBM
!!$       ELSE
!!$          DeL = 0.2D0 * FSPSCL
!!$       END IF
!!$       DeL = FSDFIX * (1.D0 + (PROFDL - 1.D0) * (R(NR) / RA)**2) + FSCDBM * DCDBM
!       PROFDDL = 8.D0
       PROFDDL = 3.D0
!       PROFDDL = 2.D0
       IF (R(NR) < RA) THEN
          DeL = FSDFIX * (1.D0 + (PROFDDL - 1.D0) * (R(NR) / RA)**6) + FSCDBM * DCDBM
!          DeL = FSDFIX * (1.D0 + (PROFDDL - 1.D0) * (R(NR) / RA)**2) + FSCDBM * DCDBM
       ELSE
          IF(FSPSCL == 0.D0) THEN
             factor_bohm = (FSDFIX * PROFDDL + FSCDBM * DCDBM) &
                  &  / (PTeV(NRA) * rKeV / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
             DeL = factor_bohm * PTeV(NR) * rKeV / (16.D0 * AEE * BBL)
          ELSE
!             DeL = FSPSCL
             DeL = FSPSCL * 0.15d0
          END IF
!!$          DeL = FSDFIX * PROFDDL + FSCDBM * DCDBM
       END IF
       DeL = 1.D0
       ! Particle diffusivity
       De(NR)   = De0   * DeL
       Di(NR)   = Di0   * DeL

       IF (R(NR) < RA) THEN
          DeL = FSDFIX * (1.D0 + (PROFDL -1.D0) * (R(NR) / RA)**2) + FSCDBM * DCDBM
       ELSE
          IF(FSPSCL == 0.D0) THEN
             factor_bohm = (FSDFIX * PROFDL + FSCDBM * DCDBM) &
                  &  / (PTeV(NRA) * rKeV / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
!!$             DeL =  (1.D0 - MOD(FSBOHM,2.D0)) * FSDFIX * PROFDL &
!!$                  &+ FSBOHM * PTeV(NR) * rKeV / (16.D0 * AEE * BBL)
             DeL = factor_bohm * PTeV(NR) * rKeV / (16.D0 * AEE * BBL)
          ELSE
             DeL = FSPSCL
          END IF
       END IF
!!$       ! Particle diffusivity
!!$       De(NR)   = De0   * DeL
!!$       Di(NR)   = Di0   * DeL
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

!!$       IF (R(NR) < RA) THEN
          PNBPD(NR) = PNBP0 * SP0(NR)
          PNBTG(NR) = PNBT10 * SP1(NR) + PNBT20 * SP2(NR)
          PNB(NR)   = PNBPD(NR) + PNBTG(NR)
          SNB(NR)   = (PNBP0 * SN0(NR) + PNBT10 * SN1(NR) + PNBT20 * SN2(NR)) &
               &    / (Eb * rKeV * 1.D20)
          MNB(NR)   = (PNBT10 * SN1(NR) + PNBT20 * SN2(NR)) / (Eb * rKeV * 1.D20)
          PRFe(NR)  = PRFe0 * SP3(NR)
          PRFi(NR)  = PRFi0 * SP3(NR)
!!$       ELSE
!!$          PNBPD(NR) = 0.D0
!!$          PNBTG(NR) = 0.D0
!!$          PNB(NR)   = 0.D0
!!$          SNB(NR)   = 0.D0
!!$          MNB(NR)   = 0.D0
!!$          PRFe(NR)  = 0.D0
!!$          PRFi(NR)  = 0.D0
!!$       END IF

       !     *** Loss to divertor ***

!       IF (R(NR) + DBW > RA) THEN
       IF (R(NR) > RA) THEN
!          Cs = SQRT(2.D0 * PTeV(NR) * rKeV / AMI)
          Cs = SQRT((PZ * PTeV(NR) + 3.D0 * PTiV(NR)) * rKeV / AMI)
          RL = (R(NR) - RA) / DBW! / 2.D0
          rNuL  (NR) = FSLP  * Cs / (2.D0 * PI * Q(NR) * RR) &
               &             * RL**2 / (1.D0 + RL**2)
          KAPPA = (4.D0*PI*EPS0)**2/(SQRT(AME)*rlnLe(NR)*AEE**4*Zeff)*AEE**2.5D0
          rNuLTe(NR) = FSLTE * KAPPA * (PTeV_FIX(NR)*1.D3)**2.5D0 &
               &                  /((2.D0 * PI * Q(NR) * RR)**2 * PNeV_FIX(NR)*1.D20) &
               &             * RL**2 / (1.D0 + RL**2)
!!$          rNuLTe(NR) = FSLTE * Cs / (2.D0 * PI * Q(NR) * RR) &
!!$               &             * RL**2 / (1.D0 + RL**2)
          rNuLTi(NR) = FSLTI * Cs / (2.D0 * PI * Q(NR) * RR) &
               &             * RL**2 / (1.D0 + RL**2)
          Ubpara(NR) = (BphV(NR) * UbphV(NR) + BthV(NR) * UbthV(NR)) / BBL
          IF(NR == NRMAX) Ubpara(NR) = AITKEN2P(R(NRMAX), &
               & Ubpara(NRMAX-1),Ubpara(NRMAX-2),Ubpara(NRMAX-3),&
               & R(NRMAX-1),R(NRMAX-2),R(NRMAX-3))
          UbparaL = max(Ubpara(NR), FSLP*Cs)
          rNuLB(NR) = FSRP * UbparaL / (2.D0 * PI * Q(NR) * RR) &
               &             * RL**2 / (1.D0 + RL**2)
       ELSE
          rNuL(NR) = 0.D0
          rNuLTe(NR) = 0.D0
          rNuLTi(NR) = 0.D0
          rNuLB(NR) = 0.D0
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

       ! Parallel current density
       AJPARA=(BthV(NR)*AJTH     + BphV(NR)*AJPH    )/BBL
       ! Parallel electric field
       EPARA =(BthV(NR)*EthV(NR) + BphV(NR)*EphV(NR))/BBL
       ! Total current density = parallel current density(?)
       AJ(NR)   = AJPARA
!       AJ(NR)   = AJPH
       ! Ohmic heating power
       POH(NR)  = EPARA*AJPARA
!       POH(NR)  = EthV(NR)*AJTH + EphV(NR)*AJPH
!       POH(NR)  = EphV(NR) * AJPH
       ! NB induced current density
       AJNB(NR) = (  (PZ * AEE * PNbV(NR) * 1.D20 * UbphV(NR)) * BphV(NR) &
            &      + (PZ * AEE * PNbV(NR) * 1.D20 * UbthV(NR)) * BthV(NR))/BBL! &
!            &    *(1.D0 - PZ / Zeff)
!       write(6,*) r(nr)/ra,AJ(NR),AJNB(NR)

       !     *** NBI power deposition ***

       PNBe(NR) = Eb * SNB(NR) * PNBcol_e(NR) * (1.D20 * rKeV)
       PNBi(NR) = Eb * SNB(NR) * PNBcol_i(NR) * (1.D20 * rKeV) &
            &   + AMb * Vb * MNB(NR) * (BthV(NR)*UithV(NR)+BphV(NR)*UiphV(NR))/BBL * 1.D20

       !     *** Equipartition power ***

       PEQe(NR)  = - 1.5D0 * rNuTei(NR) * PNeV(NR) * 1.D20 * (PTeV(NR) - PTiV(NR)) * rKeV
       PEQi(NR)  = - 1.5D0 * rNuTei(NR) * PNeV(NR) * 1.D20 * (PTiV(NR) - PTeV(NR)) * rKeV

       !     *** Ohmic power from Equations ***

       POHe(NR) = - AEE * EthV(NR) * PNeV(NR) * 1.D20 * UethV(NR) &
            &     - AEE * EphV(NR) * PNeV(NR) * 1.D20 * UephV(NR)
       POHi(NR) =   AEE * EthV(NR) * PNiV(NR) * 1.D20 * UithV(NR) &
            &     + AEE * EphV(NR) * PNiV(NR) * 1.D20 * UiphV(NR)

       !     *** Bremsstraulung loss ***
       !     (NRL Plasma Formulary p57 Eq. (30) (2002))

       PBr(NR) = 5.35D-37 * PZ**2 * PNeV(NR) * PNiV(NR) * 1.D40 * SQRT(PTeV(NR))
    END DO L_NR

    rNuAse(0) = AITKEN2P(R(0),rNuAse(1),rNuAse(2),rNuAse(3),R(1),R(2),R(3))
    rNuAsi(0) = AITKEN2P(R(0),rNuAsi(1),rNuAsi(2),rNuAsi(3),R(1),R(2),R(3))
    rNueNC(0) = AITKEN2P(R(0),rNueNC(1),rNueNC(2),rNueNC(3),R(1),R(2),R(3))
    rNuiNC(0) = AITKEN2P(R(0),rNuiNC(1),rNuiNC(2),rNuiNC(3),R(1),R(2),R(3))
    if(rNueNC(0) < 0.d0) rNueNC(0) = 0.d0
    if(rNuiNC(0) < 0.d0) rNuiNC(0) = 0.d0

    !     *** Ratio of CX deposition rate to IZ deposition rate ***
    !     (Riviere, NF 11 (1971) 363)

    IF(PNBH == 0.D0) THEN
       RatCX = 0.D0
    ELSE
       IF(Eb > 150.D0) THEN
          Sion = 3.6D-16 / EbL * (- 0.7783D0 + logEbL)
       ELSE
          Sion = 10.D0**(-0.8712D0 * logEbL**2 + 8.156D0 * logEbL - 38.833D0)
       END IF
       RatCX = Scx / (Scx + Sion)
    END IF

    !     *** Resistivity ***

    DO NR = 0, NRMAX
       ! +++ Original model +++
       EpsL = R(NR) / RR
       ALFA = (rNuei1(NR)+rNueNC(NR))/rNuei3(NR)*(BthV(NR)/BphV(NR))**2 &
            & + 2.D0*rNuei2(NR)/rNuei3(NR)*BthV(NR)/BphV(NR)
       ETA1(NR) = AME * (1.D0 + ALFA) * rNuei3(NR) / (PNeV(NR)*1.D20 * AEE**2) &
            &   * BphV(NR)**2 / (BphV(NR)**2 + BthV(NR)**2)

       ! +++ Sauter model +++
       ! Inverse aspect ratio
       EpsL = R(NR) / RR
       dPTeV = DERIV3(NR,R,PTeV,NRMAX,NRM,0) * RA
       dPTiV = DERIV3(NR,R,PTiV,NRMAX,NRM,0) * RA
       dPPe  = DERIV3(NR,R,PeV,NRMAX,NRM,0) * RA
       dPPi  = DERIV3(NR,R,PiV,NRMAX,NRM,0) * RA
       CALL SAUTER(PNeV(NR),PTeV(NR),dPTeV,dPPe,PNiV(NR),PTiV(NR),dPTiV,dPPi, &
            &      Q(NR),BphV(NR),RR*RA*BthV(NR),RR*BphV(NR),EpsL,RR,PZ,Zeff,ft(nr), &
            &      rlnLe_IN=rlnLe(NR),rlnLi_IN=rlnLi(NR),JBS=AJBS3(NR),ETA=ETA3(NR))
       IF(NR == 0) AJBS3(NR) = 0.D0

        ! +++ Hirshman, Hawryluk and Birge model +++
       Vte = SQRT(2.D0 * ABS(PTeV(NR)) * rKeV / AME)
       Wte = Vte / (Q(NR) * RR) ! Omega_te; transit frequency for electrons
       rNuAsE_inv = EpsL**1.5D0 * Wte / (SQRT(2.D0) * rNuei(NR))
       EFT  = ft(NR) * rNuAsE_inv / (rNuAsE_inv + (0.58D0 + 0.2D0 * Zeff))
       ! Spitzer resistivity for hydrogen plasma (parallel direction)
       ETAS(NR) = CORR(1.D0) * AME * rNuei(NR) / (PNeV(NR) * 1.D20 * AEE**2)
       CR   = 0.56D0 * (3.D0 - Zeff) / ((3.D0 + Zeff) * Zeff)
       ETA4(NR) = ETAS(NR) * Zeff * (1.D0 + 0.27D0 * (Zeff - 1.D0)) &
            &   /((1.D0 - EFT) * (1.D0 - CR * EFT) * (1.D0 + 0.47D0 * (Zeff - 1.D0)))
       IF(FSNC /= 0) THEN
          select case(MDLETA)
          case(0)
             ETA(NR) = ETA1(NR)
          case(1)
             ETA(NR) = ETA2(NR)
          case(2)
             ETA(NR) = ETA3(NR)
          case(3)
             ETA(NR) = ETA4(NR)
          case default
             ETA(NR) = ETA2(NR)
          end select
       ELSE
          ! Spitzer resistivity when no neoclassical effects
          ETA(NR) = ETAS(NR)
       END IF

       ! Ohmic current density
!       AJOH(NR) = EphV(NR) / ETA(NR)
       EPARA =(BthV(NR)*EthV(NR) + BphV(NR)*EphV(NR))/BBL
       AJOH(NR) = EPARA / ETA(NR)
!       if(nt==100.or.nt==200) write(6,*) r(nr)/ra,epara,eta(nr)
    END DO

    !     ***** Ion Orbit Loss *****

    SiLC  (0:NRMAX) = 0.D0
    SiLCth(0:NRMAX) = 0.D0
    SiLCph(0:NRMAX) = 0.D0
    rNuOL (0:NRMAX) = 0.D0
    IF (ABS(FSLC) > 0.D0) THEN
       IF(MDLC == 1) THEN
          ! K. C. Shaing, Phys. Fluids B 4 (1992) 3310
          do nr=1,nrmax
             EpsL = R(NR) / RR
             Vti = SQRT(2.D0 * PTiV(NR) * rKeV / AMI)
             Wti = Vti / (Q(NR) * RR)
             rNuAsIL = SQRT(2.D0) * rNuii(NR) / (EpsL**1.5D0 * Wti)
             BBL = sqrt(BphV(NR)**2 + BthV(NR)**2)
             IF(FSLC == 1.D0) THEN
                rNuOL(NR) = 2.25D0 * rNuii(NR) / (sqrt(PI) * sqrt(2.D0 * EpsL)) &
                     &   * EXP(-(rNuAsIL**0.25D0 + AEE * BBL / (BphV(NR) * Vti * AMI) &
                     &         * ABS(- AphV(NR) + AphV(NRA)) / sqrt(2.D0 * EpsL))**2)
             ELSE
                SiLC(NR) = - 2.25D0 * PNiV(NR) * rNuii(NR) / (sqrt(PI) * sqrt(2.D0 * EpsL)) &
                     &   * EXP(-(rNuAsIL**0.25D0 + AEE * BBL / (BphV(NR) * Vti * AMI) &
                     &         * ABS(- AphV(NR) + AphV(NRA)) / sqrt(2.D0 * EpsL))**2)
                SiLCth(NR) = SiLC(NR) * AMI * UithV(NR) * R(NR)
                SiLCph(NR) = SiLC(NR) * AMI * UiphV(NR)
             END IF
          end do
       ELSEIF(MDLC == 2) THEN
          !     S. -I. Itoh and K. Itoh, Nucl. Fusion 29 (1989) 1031
          IF(FSLC == 1.D0) THEN
             ! RLOSS : Numerical coefficient proportional to the relative number of ions
             !         in the loss cone in velocity space
             RLOSS = 0.1D0
             rNuOL(0) = 0.D0
             DO NR = 1, NRMAX
                EpsL = R(NR) / RR
                Vti = SQRT(PTiV(NR) * rKeV / AMI)
                RhoIT = Vti * AMI / (PZ * AEE * BthV(NR))
                RL = (R(NR) - (RA - 1.5D0 * RhoIT)) / DBW ! Alleviation factor
                IF(R(NR) > (RA - RhoIT)) THEN
!                IF(ABS(RA - R(NR)) <= RhoIT .AND. R(NR) < RA) THEN
                   ExpArg = -2.D0 * EpsL * (ErV(NR) / BthV(NR))**2 / Vti**2
                   ExpArg = ExpArg * (R(NR) / RA)**2
                   rNuOL(NR) = RLOSS * rNuii(NR) / SQRT(EpsL) * EXP(ExpArg) &
                        &    * RL**2 / (1.D0 + RL**2)
                ELSE
                   rNuOL(NR) = 0.D0
                END IF
             END DO
          ELSE
             DO NR = 1, NRA
                EpsL = R(NR) / RR
                Vti = SQRT(PTiV(NR) * rKeV / AMI)
                RhoIT = Vti * AMI / (PZ * AEE * BthV(NR))
                RhoIT = MIN(RhoIT,0.1D0)
                Wti = Vti / (Q(NR) * RR)
                rNuAsI_inv = EpsL**1.5D0 * Wti / (SQRT(2.D0) * rNuii(NR))
                ExpArg = 2.D0 * EpsL / Vti**2 * (ErV(NR) / BthV(NR))**2
                AiP = rNuii(NR) * SQRT(EpsL) * rNuAsI_inv / (1.D0 + rNuAsI_inv) &
                     & * EXPV(- ExpArg)
                DO NR1 = NRA, NRMAX
                   DISTAN = (R(NR1) - R(NR)) / RhoIT
                   SiLCL = AiP * EXPV( - DISTAN**2) * PNiV(NR)
                   SiLC(NR) = SiLC(NR) - SiLCL
                   SiLC(NR1) = SiLC(NR1) + SiLCL * R(NR) / R(NR1)
                   SiLCthL = SiLCL * AMI * UithV(NR) * R(NR)
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

!!$          ! *** SiLC correction (int_0^b r * SiLC dr = 0) ***
!!$
!!$          CALL VALINT_SUB(SiLC,NRA-1,VAL1)
!!$          CALL VALINT_SUB(SiLC,NRMAX,VAL2,NRA+2)
!!$          VAL = VAL1 + VAL2
!!$          CALL INV_INT(NRA,VAL,SiLC(NRA-1),SiLC(NRA+1),VAL1)
!!$          SiLC(NRA) = FSLC * VAL1
!!$
!!$          SiLCth(0:NRMAX) = FSLC * SiLC(0:NRMAX) * AMI * UithV(0:NRMAX) * R(0:NRMAX)
!!$          SiLCph(0:NRMAX) = FSLC * SiLC(0:NRMAX) * AMI * UiphV(0:NRMAX)
!!$
!!$          ! *************************************************
          END IF
       END IF
    END IF

    !     ***** Ripple loss transport *****

    IF(ABS(FSRP) > 0.D0) THEN
       ! Ripple well region
       DO NR = 1, NRMAX
          EpsL = R(NR) / RR
          theta1  = 0.d0
          i = 0
          imax = 101
          dlt = 1.d0 / (imax - 1)
          do 
             i = i + 1
             if(i == imax) then
                write(6,*) "ERROR! Ripple Amplitude."
                exit
             end if
             theta1 = theta1 + PI * dlt
             width0 = ripple(NR,theta1)
             width1 = EpsL * sin(theta1) / (NTCOIL * Q(NR))
             if(abs(width0 - width1) < 1.d-6) exit
             if(width0 < width1) then
                theta1  = theta1 - PI * dlt
                dlt = 0.1d0 * dlt
                i = 0
                cycle
             end if
          end do
          ARC = 2.d0 * (theta1 / PI)
          th1(nr) = theta1

          theta2 = PI
          i = 0
          imax = 101
          dlt = 1.d0 / (imax - 1)
          do 
             i = i + 1
             if(i == imax) then
                write(6,*) "ERROR! Ripple Amplitude."
                exit
             end if
             theta2 = theta2 - PI * dlt
             width0 = ripple(NR,theta2)
             width1 = EpsL * sin(theta2) / (NTCOIL * Q(NR))
             if(abs(width0 - width1) < 1.d-6) exit
             if(width0 < width1) then
                theta2  = theta2 + PI * dlt
                dlt = 0.1d0 * dlt
                i = 0
                cycle
             end if
          end do
          ARC = ARC + 2.d0 * ((PI - theta2) / PI)
          RATIO(NR) = ARC / (2.d0 * PI)
          th2(nr) = theta2

          !  Convective loss (the product of the fraction of ripple trapped particle times
          !                   vertical drift velocity)
          Ubrp(NR) = FSRP * 0.5D0 &
               &*(  0.5D0 * AMb * Vb**2 / (PZ * AEE * RR * SQRT(BphV(NR)**2 + BthV(NR)**2)) &
               &  * (theta1*sin(theta1) + (PI - theta2)*sin(theta2)) / (PI + theta1 - theta2))
!!$          rNubL(NR) = Ubrp(NR) / (R(NR) * sin(theta1))
       END DO
!!$       RV0 = AITKEN2P(R(0),r(1)*(pi-th2(1)),r(1)*th1(1),r(2)*th1(2),-R(1),R(1),R(2))
!!$       Ubrp(0) = AITKEN2P(R(0),Ubrp(1),Ubrp(2),Ubrp(3),R(1),R(2),R(3))
!!$       rNubL(0) = Ubrp(0) / RV0
!       Ubrp(0) = 0.5D0 * 0.5D0 * AMb * Vb**2 / (PZ * AEE * RR * SQRT(BphV(0)**2 + BthV(0)**2))
!       Ubrp(0) = 0.D0
       Ubrp(0) = AITKEN2P(R(0),Ubrp(1),Ubrp(2),Ubrp(3),R(1),R(2),R(3))

!!$       CALL SPL1D(R,PNbrpV,DERIV,U,NRMAX+1,0,IER)
!!$       do nr = 0, nrmax
!!$          if(R(NR) <= RB * cos(th1(nr))) then
!!$             tmp = r(nr)*cos(th1(nr))
!!$             call wherenr(r,tmp,nrl,Left,Right)
!!$             CALL SPL1DF(tmp,PNbrpL(NR),R,U,NRMAX+1,IER)
!!$             PNbrpLV(NRL-1) = PNbrpLV(NRL-1) + Left  * PNbrpL(NR)
!!$             PNbrpLV(NRL)   = PNbrpLV(NRL)   + Right * PNbrpL(NR)
!!$!             write(6,*) nrl,sngl(tmp),sngl(r(nrl-1)),sngl(r(nrl)),sngl(PNbrpL(NR)),sngl(Left  * PNbrpL(NR)),sngl(Right * PNbrpL(NR)),sngl(PNbrpLV(NRL-1)),sngl(PNbrpLV(NRL))
!!$          end if
!!$       end do

       !  Diffusive loss
       !  -- Collisional diffusion of trapped fast particles
       do nr = 1, nrmax
          EpsL = R(NR) / RR

          ! rhob : Larmor radius of beam ions
          rhob = AMb * Vb / (PZ * AEE * SQRT(BphV(NR)**2 + BthV(NR)**2))
          ! DltR2 : Square step size of banana particle
          DltR2 = PI * NTCOIL * (Q(NR) / EpsL)**3 * DltRP(NR)**2 * rhob**2

!!$          ! effective collisional frequency
!!$          rNueff = Q(NR)**2 * NTCOIL**2 / EpsL * rNuD(NR)
          ! rNubnc : bounce frequency of beam ions
!!          rNubnc = SQRT(EpsL) * Vb / (2.D0 * PI * Q(NR) * RR)
          rNubnc = SQRT(EpsL) * Vb / (10.5D0 * Q(NR) * RR)
          ! DCB : confined banana diffusion coefficient
          DCB = NTCOIL**2.25D0*Q(NR)**3.25D0*RR*rhob*DltRP(NR)**1.5d0*rNuD(NR)/EpsL**2.5D0
          ! DRP : ripple-plateau diffusion coefficient
!!          DRP = FSRP * rNubnc * DltR2
          DRP = FSRP * 0.25D0 * DltR2 * rNubnc

          ! Collisional ripple well diffusion
!!$          if (rNueff < rNubnc) then
!!$             Dbrp(NR) = DCB
!!$          else
!!$             Dbrp(NR) = DRP
!!$          end if
          Dbrp(NR) = DCB * DRP / (DCB + DRP)

          ! Dltcr : criterion of stochastic diffusion
          Dltcr = (EpsL / (PI * NTCOIL * Q(NR)))**1.5D0 / (rhob * dQdr(NR))
          ! Collisionless stochastic (ergodic) diffusion (whose value is almost
          ! equivalent to that of ripple-plateau diffusion)
          if(DltRP(NR) > Dltcr) Dbrp(NR) = DRP
!          Dbrp(NR) = 0.D0
       end do
       Dbrp(0) = AITKEN2P(R(0),Dbrp(1),Dbrp(2),Dbrp(3),R(1),R(2),R(3))
    ELSE
       Ubrp(0:NRMAX) = 0.D0
       Dbrp(0:NRMAX) = 0.D0
    END IF

    RETURN
  END SUBROUTINE TXCALC

  pure REAL(8) FUNCTION CORR(X)
    ! X is the effective charge number
    real(8), intent(in) :: X

    CORR = (1.D0 + 1.198D0 * X + 0.222D0 * X**2) &
    &    / (1.D0 + 2.966D0 * X + 0.753D0 * X**2) * X

  END FUNCTION CORR

  pure REAL(8) FUNCTION NUD(X)
    real(8), intent(in) :: X

    NUD = SQRT(1.D0 + X**2) + X**2 * LOG(X / (1.D0 + SQRT(1.D0 + X**2)))

  END FUNCTION NUD

  pure real(8) function NBIi_ratio(x) result(f)
    use physical_constants, only : PI
    real(8), intent(in) :: x

    if (x == 0.d0) then
       f = 1.d0
    else
       f = 1.d0 / x * (  1.d0 / 3.d0 * log((1.d0 - sqrt(x) + x) / (1.d0 + sqrt(x))**2) &
            &          + 2.d0 / sqrt(3.d0) * (atan((2.d0 * sqrt(x) - 1.d0) / sqrt(3.d0)) &
            &          + PI / 6.d0))
    end if

  end function NBIi_ratio

  pure real(8) function ripple(NR,theta)
    integer, intent(in) :: NR
    real(8), intent(in) :: theta

    ripple = DltRP0 * (       ((RR + R(NR) * cos(theta)) / (RR + RA))**(NTCOIL-1) &
         &             + DIN *((RR - RA) / (RR + R(NR) * cos(theta)))**(NTCOIL+1))

  end function ripple

!!$  ! Search minimum radial number NR satisfying R(NR) > X.
!!$
!!$  subroutine wherenr(R,X,NR,Left,Right)
!!$    real(8), dimension(0:NRMAX), intent(in) :: R
!!$    real(8), intent(in) :: X
!!$    integer, intent(out) :: NR
!!$    real(8), intent(out) :: Left, Right
!!$    integer :: NRL
!!$
!!$    if(X == 0.d0) then
!!$       NR = 1
!!$       Left  = 0.d0
!!$       Right = 0.d0
!!$       return
!!$    end if
!!$
!!$    do nrl = 1, nrmax
!!$       if(r(nrl) > x) then
!!$          NR = nrl
!!$          Right = (x - r(nr-1)) / (r(nr) - r(nr-1))
!!$          Left  = 1.d0 - Right
!!$          exit
!!$       end if
!!$    end do
!!$
!!$  end subroutine wherenr
end module variables
