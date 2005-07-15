!     $Id$
!
!     ***************************************************************
!
!        Calculate mesh, etc.
!
!     ***************************************************************
!
      SUBROUTINE TXCALM

      INCLUDE 'txcomm.inc'

      INTEGER :: NR

      AMi = PA * AMP
      AMb = AMi
      DR = RB / NRMAX
      NQMAX = NQM
      UHth=1.D0/SQRT(1.D0+DBLE(NCPHI)*2)
      UHph=DBLE(NCPHI)*UHth

!  Integer mesh

      DO NR = 0, NRMAX
         R(NR) = NR * DR
      ENDDO

!  Half integer mesh

      DO NR = 0, NRMAX
         RHI(NR) = (NR + 0.5D0) * DR
      ENDDO

      RETURN
      END
!
!     ***************************************************************
!
!        Calculate variables
!
!     ***************************************************************
!
      SUBROUTINE TXCALV(XL)

      INCLUDE 'txcomm.inc'

      INTEGER :: NR
      REAL(8), DIMENSION(NQM,0:NRM) :: XL

! Half Integer Variables

      DO NR = 0, NRMAX - 1
         EphHI(NR)  = XL(LQm3,NR)
         BphHI(NR)  = XL(LQm5,NR)
         PNeHI(NR)  = XL(LQe1,NR)
         UephHI(NR) = XL(LQe4,NR)/PNeHI(NR)
         PTeHI(NR)  = XL(LQe5,NR)/PNeHI(NR)
         PNiHI(NR)  = XL(LQi1,NR)
         UiphHI(NR) = XL(LQi4,NR)/PNiHI(NR)
         PTiHI(NR)  = XL(LQi5,NR)/PNiHI(NR)
         PNbHI(NR)  = XL(LQb1,NR)
         IF(PNbHI(NR).LE.0.D0) THEN
            UbphHI(NR)=0.D0
         ELSE
            UbphHI(NR) = XL(LQb4,NR)/PNbHI(NR)
         ENDIF
         PN01HI(NR) = XL(LQn1,NR)
         PN02HI(NR) = XL(LQn2,NR)
      ENDDO

! Interpolation to Integer Variables

      DO NR = 1, NRMAX - 1
         PNeI(NR) = (PNeHI(NR-1) + PNeHI(NR)) / 2
         UephI(NR)= (UephHI(NR-1)+ UephHI(NR))/ 2
         PNiI(NR) = (PNiHI(NR-1) + PNiHI(NR)) / 2
         UiphI(NR)= (UiphHI(NR-1)+ UiphHI(NR))/ 2
         PNbI(NR) = (PNbHI(NR-1) + PNbHI(NR)) / 2
         UbphI(NR)= (UbphHI(NR-1)+ UbphHI(NR))/ 2
         PTeI(NR) = (PTeHI(NR-1) + PTeHI(NR)) / 2
         PTiI(NR) = (PTiHI(NR-1) + PTiHI(NR)) / 2
         BphI(NR) = (BphHI(NR-1) + BphHI(NR)) / 2
      ENDDO

! Extraporation for central value

!     f'(0) = 0, F(0)=f(0.5*DR), F(1)=f(1.5DR)
!     then f(0) = (9*F(0)-F(1))/8

      PNeI(0) = (9.D0 * PNeHI(0) - PNeHI(1)) / 8.D0
      UephI(0)= (9.D0 * UephHI(0)- UephHI(1))/ 8.D0
      PNiI(0) = (9.D0 * PNiHI(0) - PNiHI(1)) / 8.D0
      UiphI(0)= (9.D0 * UiphHI(0)- UiphHI(1))/ 8.D0
      PNbI(0) = (9.D0 * PNbHI(0) - PNbHI(1)) / 8.D0
      UbphI(0)= (9.D0 * UbphHI(0)- UbphHI(1))/ 8.D0
      PTeI(0) = (9.D0 * PTeHI(0) - PTeHI(1)) / 8.D0
      PTiI(0) = (9.D0 * PTiHI(0) - PTiHI(1)) / 8.D0
      BphI(0) = (9.D0 * BphHI(0) - BphHI(1)) / 8.D0

! Extraporation for boundary value

!     f'(RB) = 0, F(NR-2)=f(RB-1.5*DR), F(NR-1)=f(RB-0.5DR)
!     then f(RB) = (9*F(NR-1)-F(NR-2))/8

      PNeI(NRMAX) = (9.D0 * PNeHI(NRMAX-1) - PNeHI(NRMAX-2)) / 8.D0
      UephI(NRMAX)= 0.D0
      PNiI(NRMAX) = (9.D0 * PNiHI(NRMAX-1) - PNiHI(NRMAX-2)) / 8.D0
      UiphI(NRMAX)= 0.D0
      PNbI(NRMAX) = (9.D0 * PNbHI(NRMAX-1) - PNbHI(NRMAX-2)) / 8.D0
      UbphI(NRMAX)= 0.D0
      PTeI(NRMAX) = (9.D0 * PTeHI(NRMAX-1) - PTeHI(NRMAX-2)) / 8.D0
      PTiI(NRMAX) = (9.D0 * PTiHI(NRMAX-1) - PTiHI(NRMAX-2)) / 8.D0
      BphI(NRMAX) = (9.D0 * BphHI(NRMAX-1) - BphHI(NRMAX-2)) / 8.D0

! Integer Variables

      DO NR = 0, NRMAX
         ErI(NR)  = XL(LQm1,NR)
         EthI(NR) = XL(LQm2,NR)
         BthI(NR) = XL(LQm4,NR)
         UerI(NR) = XL(LQe2,NR)/PNeI(NR)
         UethI(NR)= XL(LQe3,NR)/PNeI(NR)
         UirI(NR) = XL(LQi2,NR)/PNiI(NR)
         UithI(NR)= XL(LQi3,NR)/PNiI(NR)
         IF (ABS(PNbI(NR)).LT.1.D-30) THEN
            UbthI(NR) = 0.D0
         ELSE
            UbthI(NR) = XL(LQb3,NR)/PNbI(NR)
         ENDIF
      ENDDO

! Interpolartion to Half-Integer Variables

      DO NR = 0, NRMAX - 1
          ErHI(NR) = (  ErI(NR) +   ErI(NR+1)) / 2
         EthHI(NR) = ( EthI(NR) +  EthI(NR+1)) / 2
         BthHI(NR) = ( BthI(NR) +  BthI(NR+1)) / 2
         UerHI(NR) = ( UerI(NR) +  UerI(NR+1)) / 2
        UethHI(NR) = (UethI(NR) + UethI(NR+1)) / 2
         UirHI(NR) = ( UirI(NR) +  UirI(NR+1)) / 2
        UithHI(NR) = (UithI(NR) + UithI(NR+1)) / 2
        UbthHI(NR) = (UbthI(NR) + UbthI(NR+1)) / 2
           QHI(NR) = ABS(RHI(NR) * BphHI(NR) / (RR * BthHI(NR)))
      ENDDO

      DO NR = 1, NRMAX
         Q(NR) = ABS(R(NR) * BphI(NR) / (RR * BthI(NR)))
      ENDDO
      Q(0) = (4 * Q(1) - Q(2)) / 3.D0
      RETURN
      END
!
!     ***************************************************************
!
!        Calculate coefficients
!
!     ***************************************************************
!
      SUBROUTINE TXCALC

      INCLUDE 'txcomm.inc'

      INTEGER :: NR, NP, NR1
      REAL(8) :: Sigma0, RL, QL, SL, PNB0, PRFe0, PRFi0, Vte, Vti, 
     &           rLnLam, EION, XXX, SiV, ScxV, Wte, Wti, EpsL,
     &           rNuAsE, rNuAsI, BBL, Va, Wpe2, rGC, dQdr, SP, rGBM,
     &           rGIC, rH, DErDr, BB1, BB2, Beta1, Beta2, DBetaDr,
     &           DCDBM, DeL, ETA, AJPH, AJTH, BN, AJPARA, EPARA, Vcr, Y,
     &           Ubst, Cs, RhoIT, ExpArg, AiP, EXPV, Uith, Uiph,DISTAN,
     &           SiLCL, SiLCthL, SiLCphL
    

!     *** Constants ***

!     Neutral corsssection

      Sigma0 = 8.8D-21

!     NBI beam velocity

      Vb =  SQRT(2 * Eb * rKEV / AMb)

!     Poloidal magnetic field on wall

      IF(FSHL.EQ.0.D0) THEN
         Bthb = rMU0 * rIP * 1.D6 / (2.D0 * PI * RB)
      ELSE
         RL=RB
         QL=(Q0-QA)*(1-(RL/RA)**2)+QA
         Bthb = BB*RL/(QL*RR)
      ENDIF

!     *** Normalization factor for heating profile ***

      SL = 0.D0
      DO NR = 0, NRMAX - 1
         IF (RHI(NR) .LT. RA) THEN
            SL = SL + 2 * PI * RHI(NR) * DR
     &                * EXP(- (RHI(NR) / RNB)**2)
     &                * (1 - (RHI(NR) / RA)** 4)
         ENDIF
      ENDDO

      PNB0 = PNBH * 1.D6 / (2 * Pi * RR * SL)

      SL = 0.D0
      DO NR = 0, NRMAX - 1
         IF (RHI(NR) .LT. RA) THEN
            SL = SL + 2 * PI * RHI(NR) * DR
     &                * EXP(- (RHI(NR) / RRF)**2)
     &                * (1 - (RHI(NR) / RA)** 4)
         ENDIF
      ENDDO

      PRFe0 = 0.5D0 * PRFH * 1.D6 / (2 * Pi * RR * SL)
      PRFi0 = 0.5D0 * PRFH * 1.D6 / (2 * Pi * RR * SL)

!     ***** Half Integer Mesh *****

      DO NR = 0, NRMAX-1

         Vte = SQRT(2 * ABS(PTeHI(NR)) * rKeV / AME)
         Vti = SQRT(2 * ABS(PTiHI(NR)) * rKeV / AMI)
         rLnLam = + 15 - LOG(ABS(PNeHI(NR))) / 2 + LOG(ABS(PTeHI(NR)))

!     *** Ionization frequency ***

         EION = 13.64D0
         XXX = MAX(PTeHI(NR) * 1.D3 / EION, 1.D-2)
         SiV = 1.D-11 * SQRT(XXX) * EXP(- 1.D0 / XXX)
     &              / (EION**1.5D0 * (6.D0 + XXX))
         rNuION(NR) = FSION * SiV * (PN01HI(NR) + PN02HI(NR)) * 1.D20

!     *** Slow neutral diffusion coefficient ***

         D01(NR) = FSD0 * V0**2
     &         / (Sigma0 * (PN01HI(NR) * V0
     &         + (PNiHI(NR) + PN02HI(NR)) * Vti) * 1.D20)

!     *** Fast neutral diffusion coefficient ***

         D02(NR) = FSD0 * Vti**2
     &         / (Sigma0 * (PNiHI(NR)
     &         + PN01HI(NR) + PN02HI(NR)) * Vti * 1.D20)

!     *** Charge exchange frequency ***

         XXX = LOG10(MAX(PTiHI(NR) * 1.D3, 50.D0))
         ScxV = 1.57D-16 * SQRT(PTiHI(NR) * 1.D3)
     &          * (XXX * XXX - 14.63D0 * XXX + 53.65D0)
         rNuiCX(NR) = FSCX * ScxV * (PN01HI(NR) + PN02HI(NR)) * 1.D20

!     *** Collision frequency (with neutral) ***

         rNu0e(NR) = (PN01HI(NR) + PN02HI(NR)) * 1.D20 * Sigma0 * Vte
         rNu0i(NR) = (PN01HI(NR) + PN02HI(NR)) * 1.D20 * Sigma0 * Vti

!     *** Collision frequency (momentum-transfer) ***

         rNuei(NR) = PNeHI(NR) * 1.D20 * Zeff  * AEE**4 * rLnLam
     &              / (3 * (2 * PI)**1.5D0
     &                   * EPS0**2 * SQRT(AME)
     &                   * (ABS(PTeHI(NR)) * rKeV)**1.5D0)
         rNuii(NR) = PNiHI(NR) * 1.D20 * PZ**4 * AEE**4 * rLnLam
     &              / (3 * SQRT(2.D0) * (2 * PI)**1.5D0
     &                   * EPS0**2 * SQRT(AMI)
     &                   * (ABS(PTiHI(NR)) * rKeV)**1.5D0)
         rNuTei(NR) = PNiHI(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam
     &               / (3 * SQRT(2 * PI) * PI * EPS0**2 * AME * AMI
     &                  * (  ABS(MAX(PTeHI(NR),PTeDIV)) * rKeV / AME
     &                     + ABS(PTiHI(NR)) * rKeV / AMI)**1.5D0)
!!!     &                  * (  ABS(PTeHI(NR)) * rKeV / AME

!     *** Toroidal Neoclassical viscosity ***

         Wte = Vte / (QHI(NR) * RR)
         Wti = Vti / (QHI(NR) * RR)
         EpsL = RHI(NR) / RR
         rNuAsE = rNuei(NR) / (EpsL**1.5D0 * Wte)
         rNuAsI = rNuii(NR) / (EpsL**1.5D0 * Wti)
         BBL = SQRT(BphHI(NR)**2 + BthHI(NR)**2)
         rNueNC(NR) = FSNC * SQRT(PI) * QHI(NR)**2
     &               * Wte * 1.78D0 * rNuAsE / (1 + 1.78D0 * rNuAsE)
         rNuiNC(NR) = FSNC * SQRT(PI) * QHI(NR)**2
     &               * Wti * 1.78D0 * rNuAsI / (1 + 1.78D0 * rNuAsI)
!     &               * (1 + EpsL**1.5D0 * rNuAsI)
!     &               / (1 + 1.44D0
!     &                      * ((EpsL**1.5D0 * rNuAsI)**2
!     &                         + ( ErHI(NR)
!     &                             / ( Vti * BthHI(NR)) )**2))
!!     &                         + ( ErHI(NR) * BBL
!!     &                             / ( Vti * BthHI(NR)**2) )**2))

!     *** Helical Neoclassical viscosity ***

         Wte = Vte * NCphi / RR
         Wti = Vti * NCphi / RR
         EpsL = EpsH * RHI(NR) / RA
         rNuAsE = rNuei(NR) / (EpsL**1.5D0 * Wte)
         rNuAsI = rNuii(NR) / (EpsL**1.5D0 * Wti)
         rNueHL(NR) = FSHL * SQRT(PI)
     &               * Wte * 1.78D0 * rNuAsE / (1 + 1.78D0 * rNuAsE)
         rNuiHL(NR) = FSHL * SQRT(PI)
     &               * Wti * 1.78D0 * rNuAsI / (1 + 1.78D0 * rNuAsI)
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

         IF (ABS(FSCDBM) .GT. 0.D0) THEN
            Va = SQRT(BBL**2 / (rMU0 * PNiHI(NR) * 1.D20 * AMI))
            Wpe2 = PNeHI(NR) * 1.D20 * AEE**2 / (AME * EPS0)
            rGC = 8.D0
            dQdr = (Q(NR+1) - Q(NR)) / DR
            S(NR) = RHI(NR) / QHI(NR) * dQdr
            BB1 = SQRT(BphI(NR  )**2 + BthI(NR  )**2)
            BB2 = SQRT(BphI(NR+1)**2 + BthI(NR+1)**2)
            Beta1 = (  PNeI(NR  ) * 1.D20 * PTeI(NR  ) * rKeV
     &               + PNiI(NR  ) * 1.D20 * PTiI(NR  ) * rKeV)
     &              / (BB1**2 / (2 * rMU0))
            Beta2 = (  PNeI(NR+1) * 1.D20 * PTeI(NR+1) * rKeV
     &               + PNiI(NR+1) * 1.D20 * PTiI(NR+1) * rKeV)
     &              / (BB2**2 / (2 * rMU0))
            DBetaDr = (Beta2 - Beta1) / DR
            Alpha(NR) = - QHI(NR)**2 * RR * DBetaDr
            IF (Alpha(NR) .GT. 0.D0) THEN
               SP = S(NR) - Alpha(NR)
            ELSE
               SP = Alpha(NR) - S(NR)
            ENDIF
            IF (SP .LT. 0.D0) THEN
               rGBM = 1 / SQRT(2 * (1 - 2 * SP)
     &                         * (1 - 2 * SP + 3 * SP**2))
            ELSE
               rGBM = (1 + 9 * SQRT(2.D0) * SP**2.5D0)
     &                / (SQRT(2.D0)
     &                   * (1 - 2 * SP + 3 * SP**2 + 2 * SP**3))
            ENDIF
            rKappa(NR) = - RHI(NR) / RR * (1 - 1 / QHI(NR)**2)
            If (Alpha(NR) * rKappa(NR) .LT. 0.D0) THEN
               rGIC = 0.D0
            ELSE
               rGIC = ABS(rKappa(NR))**1.5D0 / S(NR)**2
            ENDIF
            FCDBM(NR) = DMAX1(rGBM, rGIC)
            IF(NR.EQ.0) THEN
               rH=0.D0
            ELSE
               DErDr = (ErI(NR+1) / (R(NR+1) * BB2)
     &                 - ErI(NR) / (R(NR) * BB1)) / DR
               rH = QHI(NR) * RR * RHI(NR) *  DErDr / (Va * S(NR))
            ENDIF
            rG1h2(NR) = 1.D0 / (1.D0 + rG1 * rH**2)
            DCDBM = rGC * FCDBM(NR) * rG1h2(NR) * ABS(Alpha(NR))**1.5D0
     &              * VC**2 / Wpe2 * Va / (QHI(NR) * RR)
!            write(6,*)DCDBM
!            DCDBM = MAX(DCDBM,1.D-05)
         ELSE
            rG1h2(NR)  = 0.D0
            FCDBM(NR)  = 0.D0
            S(NR)      = 0.D0
            Alpha(NR)  = 0.D0
            rKappa(NR) = 0.D0
            DCDBM      = 0.D0
         ENDIF
         IF (RHI(NR) .LT. RA) THEN
            DeL = FSDFIX * (1 + (PROFD -1) * (RHI(NR) / RA)**2)
     &            + FSCDBM * DCDBM
         ELSE
            DeL = FSDFIX * PROFD
     &            + FSBOHM * PTeHI(NR) * rKEV / (16 * AEE * BBL)
     &            + FSPSCL
         ENDIF
         De(NR)   = De0   * DeL
         Di(NR)   = Di0   * DeL
         rMue(NR) = rMue0 * DeL
         rMui(NR) = rMui0 * DeL
         Chie(NR) = Chie0 * DeL
         Chii(NR) = Chii0 * DeL

         WPM(NR) = WPM0 * PTeHI(NR) * rKeV / (RA**2 * AEE * BphHI(NR))
         FWthe(NR) = AEE**2         * BphHI(NR)**2 * De(NR)
     &            / (PTeHI(NR) * rKeV)
         FWthi(NR) = AEE**2 * PZ**2 * BphHI(NR)**2 * Di(NR)
     &            / (PTiHI(NR) * rKeV)

!     *** Heating profile ***

         IF (RHI(NR) .LT. RA) THEN
            PNB(NR) = PNB0 * EXP(- RHI(NR)**2 / RNB**2)
     &                    * (1 - (RHI(NR) / RA)** 4)
            SNB(NR) = PNB(NR) / (Eb * rKEV * 1.D20)
            PRFe(NR)= PRFe0 * EXP(- RHI(NR)**2 / RRF**2)
     &                     * (1 - (RHI(NR) / RA)** 4)
            PRFi(NR)= PRFi0 * EXP(- RHI(NR)**2 / RRF**2)
     &                     * (1 - (RHI(NR) / RA)** 4)
         ELSE
            PNB(NR) =0.D0
            SNB(NR) =0.D0
            PRFe(NR)=0.D0
            PRFi(NR)=0.D0
         ENDIF

!     Current

         ETA=AME*rNuei(NR)/(PNeHI(NR)*1.D20*AEE**2)

         AJPH  = -      AEE * PNeHI(NR) * 1.D20 * UephHI(NR)
     &           + PZ * AEE * PNiHI(NR) * 1.D20 * UiphHI(NR)
     &           + PZ * AEE * PNbHI(NR) * 1.D20 * UbphHI(NR)
         AJTH  = -      AEE * PNeHI(NR) * 1.D20 * UethHI(NR)
     &           + PZ * AEE * PNiHI(NR) * 1.D20 * UithHI(NR)
     &           + PZ * AEE * PNbHI(NR) * 1.D20 * UbthHI(NR)

         BN=SQRT(BthHI(NR)**2+BphHI(NR)**2)
         AJPARA=(BthHI(NR)*AJTH     + BphHI(NR)*AJPH    )/BN
         EPARA =(BthHI(NR)*EthHI(NR) + BphHI(NR)*EphHI(NR))/BN
         AJ(NR)   = AJPARA
         AJOH(NR) = EPARA/ETA
!         POH(NR)  = EPARA*AJPARA
         POH(NR)  = EthHI(NR)*AJTH + EphHI(NR)*AJPH    
         AJNB(NR) = PZ * AEE * PNbHI(NR) * 1.D20 * UbphI(NR)

!     *** Collision frequency (momentum-transfer with beam) ***

         Vcr = (3 * SQRT(PI / 2) * PNiHI(NR) * PZ**2 / PNeHI(NR)
     &          * AME / AMI
     &          * (ABS(PTeHI(NR)) * rKeV / AME)**1.5D0
     &         )**(1/3.D0)
         Y = Vb / Vcr
         IF (Y.GT.0.D0) THEN
            Ubst = 3.D0 / LOG(1 + Y**3) * Vb
         ELSE
            Ubst = 0.D0
         ENDIF
         rNube(NR) = PNeHI(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam
     &              / (3 * (2 * PI)**1.5D0
     &                   * EPS0**2 * AMb * AME
     &                   * (ABS(PTeHI(NR)) * rKeV / AME)**1.5D0)
         rNubi(NR) = PNiHI(NR) * 1.D20 * PZ**2 * PZ**2 * AEE**4 * rLnLam
     &              / (4 * PI * EPS0**2 * AMb)
     &              * (1 / AMb + 1 / AMI)
     &              * 1 / (+ Ubst**3
     &                     + 9 * SQRT(3 * PI) / 4
     &                   * (ABS(PTiHI(NR)) * rKeV / AMI)**1.5D0)
         IF (Y.GT.0.D0) THEN
            rNuB(NR) = rNube(NR) * 3.D0 / LOG(1 + Y**3)
         ELSE
            rNuB(NR) = 0.D0
         ENDIF

!     *** Loss to divertor ***

         IF (RHI(NR) .GT. RA) THEN
            Cs = SQRT(PTeHI(NR) * rKeV / AMI)
!!!!            rNuL(NR) = FSLP * Cs / (2 * PI * QHI(NR) * RR
!!!!     &                        * LOG(0.3D0 / (RHI(NR) - RA)))
            rNuL(NR) = FSLP * Cs
     &              / (2 * PI * QHI(NR) * RR
     &                 * (1.D0 + LOG(1.D0 + rLT / (RHI(NR) - RA))))
         ELSE
            rNuL(NR) = 0.D0
         ENDIF

      ENDDO

!     ***** Ion Orbit Loss *****

      DO NR = 0, NRMAX - 1
         SiLC(NR)   = 0.D0
         SiLCth(NR) = 0.D0
         SiLCph(NR) = 0.D0
      ENDDO
      IF (ABS(FSLC) .GT. 0.D0) THEN
         NP = NINT(RA / DR)
         DO NR = 1, NP - 1
            EpsL = RHI(NR) / RR
            Vti = SQRT(2 * PTeHI(NR) * rKeV / AMI)
            RhoIT = Vti * AMI / (PZ * AEE * BthHI(NR))
            RhoIT = MIN(RhoIT,0.1D0)
            rNuAsI = rNuii(NR) * QHI(NR) * RR / (EpsL**1.5D0 * Vti)
            ExpArg = 2 * EpsL / Vti**2 * (ErI(NR) / BthI(NR))**2
            AiP = rNuii(NR) * SQRT(EpsL) / (1 + rNuAsI) * EXPV(- ExpArg)
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
            ENDDO
         ENDDO
         DO NR = 0, NRMAX - 1
            SiLC(NR)   = FSLC * SiLC(NR)
            SiLCth(NR) = FSLC * SiLCth(NR)
            SiLCph(NR) = FSLC * SiLCph(NR)
         ENDDO
      ENDIF
      RETURN
      END
