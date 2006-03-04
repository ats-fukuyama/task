!     $Id$
module variables
contains

!***************************************************************
!
!        Calculate variables
!
!***************************************************************

  SUBROUTINE TXCALV(XL)

    INCLUDE 'txcomm.inc'

    INTEGER :: NR
    REAL(8), DIMENSION(NQM,0:NRM) :: XL

    ! Half integer variables

    EphHI(0:NRMAX-1)  = XL(LQm3,0:NRMAX-1)
    BphHI(0:NRMAX-1)  = XL(LQm5,0:NRMAX-1)
    PNeHI(0:NRMAX-1)  = XL(LQe1,0:NRMAX-1)
    UephHI(0:NRMAX-1) = XL(LQe4,0:NRMAX-1)/PNeHI(0:NRMAX-1)
    PTeHI(0:NRMAX-1)  = XL(LQe5,0:NRMAX-1)/PNeHI(0:NRMAX-1)
    PNiHI(0:NRMAX-1)  = XL(LQi1,0:NRMAX-1)
    UiphHI(0:NRMAX-1) = XL(LQi4,0:NRMAX-1)/PNiHI(0:NRMAX-1)
    PTiHI(0:NRMAX-1)  = XL(LQi5,0:NRMAX-1)/PNiHI(0:NRMAX-1)
    PNbHI(0:NRMAX-1)  = XL(LQb1,0:NRMAX-1)
    DO NR = 0, NRMAX - 1
       IF(PNbHI(NR) <= 0.D0) THEN
          UbphHI(NR)=0.D0
       ELSE
          UbphHI(NR) = XL(LQb4,NR)/PNbHI(NR)
       END IF
    END DO
    PN01HI(0:NRMAX-1) = XL(LQn1,0:NRMAX-1)
    PN02HI(0:NRMAX-1) = XL(LQn2,0:NRMAX-1)

    ! Interpolation to integer variables

    BphI(1:NRMAX-1) = 0.5D0 * (BphHI(0:NRMAX-2) + BphHI(1:NRMAX-1))
    PNeI(1:NRMAX-1) = 0.5D0 * (PNeHI(0:NRMAX-2) + PNeHI(1:NRMAX-1))
    UephI(1:NRMAX-1)= 0.5D0 * (UephHI(0:NRMAX-2)+ UephHI(1:NRMAX-1))
    PTeI(1:NRMAX-1) = 0.5D0 * (PTeHI(0:NRMAX-2) + PTeHI(1:NRMAX-1))
    PNiI(1:NRMAX-1) = 0.5D0 * (PNiHI(0:NRMAX-2) + PNiHI(1:NRMAX-1))
    UiphI(1:NRMAX-1)= 0.5D0 * (UiphHI(0:NRMAX-2)+ UiphHI(1:NRMAX-1))
    PTiI(1:NRMAX-1) = 0.5D0 * (PTiHI(0:NRMAX-2) + PTiHI(1:NRMAX-1))
    PNbI(1:NRMAX-1) = 0.5D0 * (PNbHI(0:NRMAX-2) + PNbHI(1:NRMAX-1))
    UbphI(1:NRMAX-1)= 0.5D0 * (UbphHI(0:NRMAX-2)+ UbphHI(1:NRMAX-1))

    ! Extraporation for central value

    !     f'(0) = 0, F(0)=f(0.5*DR), F(1)=f(1.5DR)
    !     then f(0) = (9*F(0)-F(1))/8

    BphI(0) = (9.D0 * BphHI(0) - BphHI(1)) / 8.D0
    PNeI(0) = (9.D0 * PNeHI(0) - PNeHI(1)) / 8.D0
    UephI(0)= (9.D0 * UephHI(0)- UephHI(1))/ 8.D0
    PTeI(0) = (9.D0 * PTeHI(0) - PTeHI(1)) / 8.D0
    PNiI(0) = (9.D0 * PNiHI(0) - PNiHI(1)) / 8.D0
    UiphI(0)= (9.D0 * UiphHI(0)- UiphHI(1))/ 8.D0
    PTiI(0) = (9.D0 * PTiHI(0) - PTiHI(1)) / 8.D0
    PNbI(0) = (9.D0 * PNbHI(0) - PNbHI(1)) / 8.D0
    UbphI(0)= (9.D0 * UbphHI(0)- UbphHI(1))/ 8.D0

    ! Extraporation for boundary value

    !     f'(RB) = 0, F(NR-2)=f(RB-1.5*DR), F(NR-1)=f(RB-0.5DR)
    !     then f(RB) = (9*F(NR-1)-F(NR-2))/8

    BphI(NRMAX) = (9.D0 * BphHI(NRMAX-1) - BphHI(NRMAX-2)) / 8.D0
    PNeI(NRMAX) = (9.D0 * PNeHI(NRMAX-1) - PNeHI(NRMAX-2)) / 8.D0
    UephI(NRMAX)= 0.D0
    PTeI(NRMAX) = (9.D0 * PTeHI(NRMAX-1) - PTeHI(NRMAX-2)) / 8.D0
    PNiI(NRMAX) = (9.D0 * PNiHI(NRMAX-1) - PNiHI(NRMAX-2)) / 8.D0
    UiphI(NRMAX)= 0.D0
    PTiI(NRMAX) = (9.D0 * PTiHI(NRMAX-1) - PTiHI(NRMAX-2)) / 8.D0
    PNbI(NRMAX) = (9.D0 * PNbHI(NRMAX-1) - PNbHI(NRMAX-2)) / 8.D0
    UbphI(NRMAX)= 0.D0

    ! Integer variables

    ErI(0:NRMAX)  = XL(LQm1,0:NRMAX)
    EthI(0:NRMAX) = XL(LQm2,0:NRMAX)
!    do nr=0,nrmax
!       write(6,'(I3,4F18.10)') nr,XL(LQm1,NR),XL(LQe1,NR),XL(LQi1,NR),XL(LQb1,NR)
!       write(6,'(I3,3F18.10)') nr,XL(LQe2,NR),XL(LQn1,NR),XL(LQn2,NR)
!       write(6,'(I3,2F18.10)') nr,XL(LQe5,NR),XL(LQi5,NR)
!       write(6,'(I3,2F18.10)') nr,XL(LQm2,NR),XL(LQm5,NR)
!    end do
    BthI(0:NRMAX) = XL(LQm4,0:NRMAX)
    UerI(0:NRMAX) = XL(LQe2,0:NRMAX)/PNeI(0:NRMAX)
    UethI(0:NRMAX)= XL(LQe3,0:NRMAX)/PNeI(0:NRMAX)
    UirI(0:NRMAX) = XL(LQi2,0:NRMAX)/PNiI(0:NRMAX)
    UithI(0:NRMAX)= XL(LQi3,0:NRMAX)/PNiI(0:NRMAX)
    DO NR = 0, NRMAX
       IF (ABS(PNbHI(NR)) < 1.D-30) THEN
          UbthI(NR) = 0.D0
       ELSE
          UbthI(NR) = XL(LQb3,NR)/PNbI(NR)
       END IF
    END DO

    ! Interpolartion to half-integer variables

    ErHI(0:NRMAX-1)  = 0.5D0 * (  ErI(0:NRMAX-1) +   ErI(1:NRMAX))
    EthHI(0:NRMAX-1) = 0.5D0 * ( EthI(0:NRMAX-1) +  EthI(1:NRMAX))
    BthHI(0:NRMAX-1) = 0.5D0 * ( BthI(0:NRMAX-1) +  BthI(1:NRMAX))
    UerHI(0:NRMAX-1) = 0.5D0 * ( UerI(0:NRMAX-1) +  UerI(1:NRMAX))
    UethHI(0:NRMAX-1)= 0.5D0 * (UethI(0:NRMAX-1) + UethI(1:NRMAX))
    UirHI(0:NRMAX-1) = 0.5D0 * ( UirI(0:NRMAX-1) +  UirI(1:NRMAX))
    UithHI(0:NRMAX-1)= 0.5D0 * (UithI(0:NRMAX-1) + UithI(1:NRMAX))
    UbthHI(0:NRMAX-1)= 0.5D0 * (UbthI(0:NRMAX-1) + UbthI(1:NRMAX))
    QHI(0:NRMAX-1)   = ABS(RHI(0:NRMAX-1) * BphHI(0:NRMAX-1) &
         &                / (RR * BthHI(0:NRMAX-1)))

    Q(1:NRMAX) = ABS(R(1:NRMAX) * BphI(1:NRMAX) / (RR * BthI(1:NRMAX)))
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
    use libraries, only : EXPV
    INCLUDE 'txcomm.inc'

    INTEGER :: NR, NP, NR1
    REAL(8) :: Sigma0, QL, SL, PNB0, PRFe0, PRFi0, Vte, Vti,  &
         &     rLnLam, EION, XXX, SiV, ScxV, Wte, Wti, EpsL, &
         &     rNuAsE, rNuAsI, BBL, Va, Wpe2, rGC, dQdr, SP, rGBM, &
         &     rGIC, rH, dErdr, BB1, BB2, Beta1, Beta2, dBetadr, &
         &     DCDBM, DeL, ETA, AJPH, AJTH, BN, AJPARA, EPARA, Vcr, Y, &
         &     Ubst, Cs, RhoIT, ExpArg, AiP, Uith, Uiph,DISTAN, &
         &     SiLCL, SiLCthL, SiLCphL

    !     *** Constants ***

    !     Neutral cross section

    Sigma0 = 8.8D-21

    !     NBI beam velocity

    Vb =  SQRT(2.D0 * Eb * rKEV / AMb)

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
    !
    !  For NBI heating
    SL = 0.D0
    DO NR = 0, NRMAX - 1
       IF (RHI(NR) < RA) THEN
          SL = SL + 2.D0 * PI * RHI(NR) * DR * EXP(- (RHI(NR) / RNB)**2) &
               &  * (1.D0 - (RHI(NR) / RA)** 4)
       END IF
    END DO

    PNB0 = PNBH * 1.D6 / (2.D0 * Pi * RR * SL)

    !  For RF heating
    SL = 0.D0
    DO NR = 0, NRMAX - 1
       IF (RHI(NR) < RA) THEN
          SL = SL + 2.D0 * PI * RHI(NR) * DR * EXP(- (RHI(NR) / RRF)**2) &
               &  * (1.D0 - (RHI(NR) / RA)** 4)
       END IF
    END DO

    PRFe0 = 0.5D0 * PRFH * 1.D6 / (2.D0 * Pi * RR * SL)
    PRFi0 = 0.5D0 * PRFH * 1.D6 / (2.D0 * Pi * RR * SL)

    !     ***** Half integer mesh *****

    L_NR: DO NR = 0, NRMAX-1

       Vte = SQRT(2.D0 * ABS(PTeHI(NR)) * rKeV / AME)
       Vti = SQRT(2.D0 * ABS(PTiHI(NR)) * rKeV / AMI)
       rLnLam = 15.D0 - LOG(ABS(PNeHI(NR))) / 2.D0 + LOG(ABS(PTeHI(NR)))

       !     *** Ionization frequency ***

       EION = 13.64D0
       XXX = MAX(PTeHI(NR) * 1.D3 / EION, 1.D-2)
       SiV = 1.D-11 * SQRT(XXX) * EXP(- 1.D0 / XXX) &
            &              / (EION**1.5D0 * (6.D0 + XXX))
       rNuION(NR) = FSION * SiV * (PN01HI(NR) + PN02HI(NR)) * 1.D20

       !     *** Slow neutral diffusion coefficient ***

       D01(NR) = FSD0 * V0**2 &
            &   / (Sigma0 * (PN01HI(NR) * V0  + (PNiHI(NR) + PN02HI(NR)) &
            &      * Vti) * 1.D20)

       !     *** Fast neutral diffusion coefficient ***

       D02(NR) = FSD0 * Vti**2 &
            &   / (Sigma0 * (PNiHI(NR) + PN01HI(NR) + PN02HI(NR)) &
            &      * Vti * 1.D20)

       !     *** Charge exchange frequency ***

       XXX = LOG10(MAX(PTiHI(NR) * 1.D3, 50.D0))
       ScxV = 1.57D-16 * SQRT(PTiHI(NR) * 1.D3) &
            &          * (XXX * XXX - 14.63D0 * XXX + 53.65D0)
       rNuiCX(NR) = FSCX * ScxV * (PN01HI(NR) + PN02HI(NR)) * 1.D20

       !     *** Collision frequency (with neutral) ***

       rNu0e(NR) = (PN01HI(NR) + PN02HI(NR)) * 1.D20 * Sigma0 * Vte
       rNu0i(NR) = (PN01HI(NR) + PN02HI(NR)) * 1.D20 * Sigma0 * Vti

       !     *** Collision frequency (momentum transfer) ***

       rNuei(NR) = PNiHI(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
            &     / (3.D0 * (2.D0 * PI)**1.5D0 &
            &     * EPS0**2 * SQRT(AME) &
            &     * (ABS(PTeHI(NR)) * rKeV)**1.5D0)
       rNuii(NR) = PNiHI(NR) * 1.D20 * PZ**4 * AEE**4 * rLnLam &
            &     / (3.D0 * SQRT(2.D0) * (2.D0 * PI)**1.5D0 &
            &     * EPS0**2 * SQRT(AMI) &
            &     * (ABS(PTiHI(NR)) * rKeV)**1.5D0)
       rNuTei(NR) = PNiHI(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
            &     / (3.D0 * SQRT(2.D0 * PI) * PI * EPS0**2 * AME * AMI &
            &     * (  ABS(MAX(PTeHI(NR),PTeDIV)) * rKeV / AME &
            &        + ABS(PTiHI(NR)) * rKeV / AMI)**1.5D0)

       !     *** Toroidal neoclassical viscosity ***

       Wte = Vte / (QHI(NR) * RR) ! Omega_Te
       Wti = Vti / (QHI(NR) * RR) ! Omega_Ti
       EpsL = RHI(NR) / RR        ! Inverse aspect ratio
       rNuAsE = rNuei(NR) / (EpsL**1.5D0 * Wte)
       rNuAsI = rNuii(NR) / (EpsL**1.5D0 * Wti)
       BBL = SQRT(BphHI(NR)**2 + BthHI(NR)**2)
       rNueNC(NR) = FSNC * SQRT(PI) * QHI(NR)**2 &
            &     * Wte * 1.78D0 * rNuAsE / (1 + 1.78D0 * rNuAsE)
       rNuiNC(NR) = FSNC * SQRT(PI) * QHI(NR)**2 &
            &     * Wti * 1.78D0 * rNuAsI / (1 + 1.78D0 * rNuAsI)
       !     &               * (1 + EpsL**1.5D0 * rNuAsI)
       !     &               / (1 + 1.44D0
       !     &                      * ((EpsL**1.5D0 * rNuAsI)**2
       !     &                         + ( ErHI(NR)
       !     &                             / ( Vti * BthHI(NR)) )**2))
       !!     &                         + ( ErHI(NR) * BBL
       !!     &                             / ( Vti * BthHI(NR)**2) )**2))

       !     *** Helical neoclassical viscosity ***

       Wte = Vte * NCphi / RR
       Wti = Vti * NCphi / RR
       EpsL = EpsH * RHI(NR) / RA
       rNuAsE = rNuei(NR) / (EpsL**1.5D0 * Wte)
       rNuAsI = rNuii(NR) / (EpsL**1.5D0 * Wti)
       rNueHL(NR) = FSHL * SQRT(PI) &
            &     * Wte * 1.78D0 * rNuAsE / (1.D0 + 1.78D0 * rNuAsE)
       rNuiHL(NR) = FSHL * SQRT(PI) &
            &     * Wti * 1.78D0 * rNuAsI / (1.D0 + 1.78D0 * rNuAsI)
       !         WRITE(6,'(I5,1P4E12.4)') 
       !     &        NR,rNueNC(NR),rNuiNC(NR),rNueHL(NR),rNuiHL(NR)
       !         rNueHL(NR) = 0.D0
       !         rNuiHL(NR) = 0.D0
       !     &               * (1 + EpsL**1.5D0 * rNuAsI)
       !     &               / (1 + 1.44D0
       !     &                      * ((EpsL**1.5D0 * rNuAsI)**2
       !     &                         + ( ErHI(NR)
       !     &                             / ( Vti * BthHI(NR)) )**2))
       !!     &                         + ( ErHI(NR) * BBL
       !!     &                             / ( Vti * BthHI(NR)**2) )**2))

       !     *** Wave-particle interaction ***

       IF (ABS(FSCDBM) > 0.D0) THEN
          ! Alfven velocity
          Va = SQRT(BBL**2 / (rMU0 * PNiHI(NR) * 1.D20 * AMI))
          ! Squared plasma frequency
          Wpe2 = PNeHI(NR) * 1.D20 * AEE**2 / (AME * EPS0)
          ! Arbitrary coefficient for CDBM model
          rGC = 8.D0
          dQdr = (Q(NR+1) - Q(NR)) / DR
          ! Magnetic shear
          S(NR) = RHI(NR) / QHI(NR) * dQdr
          BB1 = SQRT(BphI(NR  )**2 + BthI(NR  )**2)
          BB2 = SQRT(BphI(NR+1)**2 + BthI(NR+1)**2)
          Beta1 = (  PNeI(NR  ) * 1.D20 * PTeI(NR  ) * rKeV &
               &   + PNiI(NR  ) * 1.D20 * PTiI(NR  ) * rKeV) &
               &  / (BB1**2 / (2.D0 * rMU0))
          Beta2 = (  PNeI(NR+1) * 1.D20 * PTeI(NR+1) * rKeV &
               &   + PNiI(NR+1) * 1.D20 * PTiI(NR+1) * rKeV) &
               &  / (BB2**2 / (2.D0 * rMU0))
          dBetadr = (Beta2 - Beta1) / DR
          ! Normalized pressure gradient
          Alpha(NR) = - QHI(NR)**2 * RR * dBetadr
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
          rKappa(NR) = - RHI(NR) / RR * (1.D0 - 1.D0 / QHI(NR)**2)
          ! Fitting function: F for interchange instability
          If (Alpha(NR) * rKappa(NR) < 0.D0) THEN
             rGIC = 0.D0
          ELSE
             rGIC = ABS(rKappa(NR))**1.5D0 / S(NR)**2
          END IF
          ! Select larger effect between ballooning and interchange modes
          FCDBM(NR) = MAX(rGBM, rGIC)
          ! ExB rotational shear
          IF(NR == 0) THEN
             rH=0.D0
          ELSE
             dErdr = (  ErI(NR+1) / (R(NR+1) * BB2) &
                  &   - ErI(NR)   / (R(NR)   * BB1)) / DR
             rH = QHI(NR) * RR * RHI(NR) *  dErdr / (Va * S(NR))
          END IF
          ! Coefficient of fitting formula
          rG1h2(NR) = 1.D0 / (1.D0 + rG1 * rH**2)
          ! Turbulent transport coefficient calculated by CDBM model
          DCDBM = rGC * FCDBM(NR) * rG1h2(NR) * ABS(Alpha(NR))**1.5D0 &
               &              * VC**2 / Wpe2 * Va / (QHI(NR) * RR)
          !write(6,*)DCDBM
          !DCDBM = MAX(DCDBM,1.D-05)
       ELSE
          rG1h2(NR)  = 0.D0
          FCDBM(NR)  = 0.D0
          S(NR)      = 0.D0
          Alpha(NR)  = 0.D0
          rKappa(NR) = 0.D0
          DCDBM      = 0.D0
       END IF

       IF (RHI(NR) < RA) THEN
          DeL = FSDFIX * (1.D0 + (PROFD -1) * (RHI(NR) / RA)**2) &
               &            + FSCDBM * DCDBM
       ELSE
          DeL = FSDFIX * PROFD &
               &       + FSBOHM * PTeHI(NR) * rKEV / (16.D0 * AEE * BBL) &
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
       WPM(NR) = WPM0 * PTeHI(NR) * rKeV / (RA**2 * AEE * BphHI(NR))
       ! Force induced by drift wave (eq.(8),(13)) on half mesh
       FWthe(NR) = AEE**2         * BphHI(NR)**2 * De(NR) &
            &            / (PTeHI(NR) * rKeV)
       FWthi(NR) = AEE**2 * PZ**2 * BphHI(NR)**2 * Di(NR) &
            &            / (PTiHI(NR) * rKeV)

       !     *** Heating profile ***

       IF (RHI(NR) < RA) THEN
          PNB(NR) = PNB0  * EXP(- RHI(NR)**2 / RNB**2) * (1 - (RHI(NR) / RA)** 4)
          SNB(NR) = PNB(NR) / (Eb * rKEV * 1.D20)
          PRFe(NR)= PRFe0 * EXP(- RHI(NR)**2 / RRF**2) * (1 - (RHI(NR) / RA)** 4)
          PRFi(NR)= PRFi0 * EXP(- RHI(NR)**2 / RRF**2) * (1 - (RHI(NR) / RA)** 4)
       ELSE
          PNB(NR) =0.D0
          SNB(NR) =0.D0
          PRFe(NR)=0.D0
          PRFi(NR)=0.D0
       END IF

       !     *** Current density profile ***

       ! Resistivity
       ETA=AME*rNuei(NR)/(PNeHI(NR)*1.D20*AEE**2)
       ! Poloidal current density
       AJPH  = -      AEE * PNeHI(NR) * 1.D20 * UephHI(NR) &
            &  + PZ * AEE * PNiHI(NR) * 1.D20 * UiphHI(NR) &
            &  + PZ * AEE * PNbHI(NR) * 1.D20 * UbphHI(NR)
       ! Toroidal current density
       AJTH  = -      AEE * PNeHI(NR) * 1.D20 * UethHI(NR) &
            &  + PZ * AEE * PNiHI(NR) * 1.D20 * UithHI(NR) &
            &  + PZ * AEE * PNbHI(NR) * 1.D20 * UbthHI(NR)

       BN=SQRT(BthHI(NR)**2+BphHI(NR)**2)
       ! Parallel current density
       AJPARA=(BthHI(NR)*AJTH      + BphHI(NR)*AJPH     )/BN
       ! Parallel electric field
       EPARA =(BthHI(NR)*EthHI(NR) + BphHI(NR)*EphHI(NR))/BN
       ! Total current density = parallel current density
       AJ(NR)   = AJPARA
       ! Ohmic current density
       AJOH(NR) = EPARA/ETA
       ! Ohmic heating power
       !POH(NR)  = EPARA*AJPARA
       POH(NR)  = EthHI(NR)*AJTH + EphHI(NR)*AJPH    
       ! NB induced current density
       AJNB(NR) = PZ * AEE * PNbHI(NR) * 1.D20 * UbphI(NR)

       !     *** Collision frequency (momentum transfer with beam) ***
       ! reference : 92/04/02, 92/04/21

       Vcr = (3.D0 * SQRT(PI / 2.D0) * PNiHI(NR) * PZ**2 / PNeHI(NR) &
            &   * AME / AMI &
            &   * (ABS(PTeHI(NR)) * rKeV / AME)**1.5D0)**(1.D0/3.D0)
       Y = Vb / Vcr
       IF (Y > 0.D0) THEN
          Ubst = 3.D0 / LOG(1.D0 + Y**3) * Vb
       ELSE
          Ubst = 0.D0
       END IF
       rNube(NR) = PNeHI(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
            &     / (3.D0 * (2.D0 * PI)**1.5D0 * EPS0**2 * AMb * AME &
            &             * (ABS(PTeHI(NR)) * rKeV / AME)**1.5D0)
       rNubi(NR) = PNiHI(NR) * 1.D20 * PZ**2 * PZ**2 * AEE**4 * rLnLam &
            &     / (4.D0 * PI * EPS0**2 * AMb) &
            &     * (1.D0 / AMb + 1.D0 / AMI) &
            &     * 1.D0 / ( Ubst**3 + 9.D0 * SQRT(3.D0 * PI) / 4.D0 &
            &     * (ABS(PTiHI(NR)) * rKeV / AMI)**1.5D0)
       IF (Y > 0.D0) THEN
          rNuB(NR) = rNube(NR) * 3.D0 / LOG(1.D0 + Y**3)
       ELSE
          rNuB(NR) = 0.D0
       END IF

       !     *** Loss to divertor ***

       IF (RHI(NR) > RA) THEN
          Cs = SQRT(PTeHI(NR) * rKeV / AMI)
!!!!        rNuL(NR) = FSLP * Cs / (2.D0 * PI * QHI(NR) * RR &
!!!!     &                  * LOG(0.3D0 / (RHI(NR) - RA)))
          rNuL(NR) = FSLP * Cs / (2.D0 * PI * QHI(NR) * RR &
               &              * (1.D0 + LOG(1.D0 + rLT / (RHI(NR) - RA))))
       ELSE
          rNuL(NR) = 0.D0
       END IF

    END DO L_NR

    !     ***** Ion Orbit Loss *****

    SiLC(0:NRMAX-1)   = 0.D0
    SiLCth(0:NRMAX-1) = 0.D0
    SiLCph(0:NRMAX-1) = 0.D0
    IF (ABS(FSLC) > 0.D0) THEN
       NP = NINT(RA / DR)
       DO NR = 1, NP - 1
          EpsL = RHI(NR) / RR
          Vti = SQRT(2 * PTeHI(NR) * rKeV / AMI)
          RhoIT = Vti * AMI / (PZ * AEE * BthHI(NR))
          RhoIT = MIN(RhoIT,0.1D0)
          rNuAsI = rNuii(NR) * QHI(NR) * RR / (EpsL**1.5D0 * Vti)
          ExpArg = 2.D0 * EpsL / Vti**2 * (ErI(NR) / BthI(NR))**2
          AiP = rNuii(NR) * SQRT(EpsL) / (1.D0 + rNuAsI) * EXPV(- ExpArg)
          Uith = 0.5D0*(UithI(NR)+UithI(NR+1))
          Uiph = 0.5D0*(UiphI(NR)+UiphI(NR+1))
          DO NR1 = NP, NRMAX - 1
             DISTAN = (RHI(NR1) - RHI(NR)) / RhoIT
             SiLCL = AiP * EXPV( - DISTAN**2) * PNiHI(NR)
             SiLC(NR) = SiLC(NR) - SiLCL
             SiLC(NR1) = SiLC(NR1) + SiLCL * RHI(NR) / RHI(NR1)
             SiLCthL = SiLCL * AMi * Uith
             SiLCth(NR) = SiLCth(NR) - SiLCthL
             SiLCth(NR1) = SiLCth(NR1) + SiLCthL * RHI(NR) / RHI(NR1)
             SiLCphL = SiLCL * AMi * Uiph
             SiLCph(NR) = SiLCph(NR) - SiLCphL
             SiLCph(NR1) = SiLCph(NR1) + SiLCphL * RHI(NR) / RHI(NR1)
          END DO
       END DO
       SiLC(0:NRMAX-1)   = FSLC * SiLC(0:NRMAX-1)
       SiLCth(0:NRMAX-1) = FSLC * SiLCth(0:NRMAX-1)
       SiLCph(0:NRMAX-1) = FSLC * SiLCph(0:NRMAX-1)
    END IF
    RETURN
  END SUBROUTINE TXCALC
end module variables
