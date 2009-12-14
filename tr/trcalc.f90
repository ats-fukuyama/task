!     $Id$
!     ***********************************************************

!           CALCULATE TRANSPORT COEFFICIENTS AND SOURCE

!     ***********************************************************

      SUBROUTINE TRCALC(IERR)

      USE TRCOMM, ONLY : AJBS, AJRF, AR1RHOG, ARRHOG, BP, DT, DVRHOG, MDLEQ0, &
           MDLEQB, MDLJBS, MDLUF, MDLPR, MDNCLS, NRAMAX, NRM, NRMAX, &
           NROMAX, NSM, NSMAX, PBCL, PBIN, PCX, PELTIM, PEX, PFCL, PFIN, &
           PI, PIE, PIN, PN, PNB, PNF, POH, PRB, PRC, PRF, PRL, PRSUM, &
           Q0, QP, RDP, RG, &
           RHOA, RR, SCX, SEX, SIE, SNB, SNF, SPE, SSIN, T, TTRHOG, RDPVRHOG
      USE tr_cytran_mod
      IMPLICIT NONE
      INTEGER(4),INTENT(OUT)    :: IERR
      INTEGER(4)                :: NR, NS, NWM, NWMMAX, NWR, NWRMAX
      REAL(8)                   :: FCTR
      REAL(8),DIMENSION(NRMAX)      :: AJWMPL,AJWRPL
      REAL(8),DIMENSION(NRMAX,NSMAX):: PWMPL, PWRPL
      CHARACTER(LEN=80)         :: KID
      real :: g1,g2

      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      IERR=0

      SIE(1:NRMAX)=0.D0
      SNF(1:NRMAX)=0.D0
      SNB(1:NRMAX)=0.D0
      POH(1:NRMAX)=0.D0
      PIE(1:NRMAX)=0.D0
      PCX(1:NRMAX)=0.D0
      PRB(1:NRMAX)=0.D0
      PRC(1:NRMAX)=0.D0
      PRL(1:NRMAX)=0.D0
      PRSUM(1:NRMAX)=0.D0
      PNB(1:NRMAX)=0.D0
      PNF(1:NRMAX)=0.D0
      PBIN(1:NRMAX)=0.D0
      PFIN(1:NRMAX)=0.D0
!      AJNB(1:NRMAX)=0.D0
      AJRF(1:NRMAX)=0.D0
      AJBS(1:NRMAX)=0.D0
      SPE(1:NRMAX,1:NSMAX)=0.D0
      PBCL(1:NRMAX,1:NSMAX)=0.D0
      PFCL(1:NRMAX,1:NSMAX)=0.D0
      IF(MDLUF.NE.0) THEN
         PRF(1:NRMAX,3:NSM)=0.D0
      ELSE
         PRF(1:NRMAX,1:NSM)=0.D0
      ENDIF

      BP(1:NRMAX)=AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR
      QP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)/(4.D0*PI**2*RDPVRHOG(1:NRMAX))
      Q0=FCTR(RG(1),RG(2),QP(1),QP(2))

!     *** RADIAL ELECTRIC FIELD ***

      CALL TRERAD

      IF(T.LT.PELTIM+0.5D0*DT.AND. T.GE.PELTIM-0.5D0*DT) CALL TRPELT
      CALL TRZEFF
      IF(MDLPR.GT.0) CALL TR_CYTRAN

      IF(MDNCLS.NE.0) THEN
         CALL TR_NCLASS(IERR)
         IF(IERR.NE.0) RETURN
      ENDIF

      CALL TRCOEF
      CALL TRLOSS
      IF(MDLUF.NE.1.AND.MDLUF.NE.3) CALL TRPWRF
      CALL TRPWNB

      IF(MDNCLS.NE.0) THEN
         CALL TRAJBS_NCLASS
      ELSE
         select case(MDLJBS)
         case(1)
            CALL TRAJBS
         case(2)
            CALL TRAJBS
         case(3)
            CALL TRAJBS
         case(4)
            CALL TRAJBSNEW
         case(5)
            CALL TRAJBSSAUTER
         case default
            CALL TRAJBS
         end select
      ENDIF

      CALL TRALPH
      CALL TRAJOH

      DO NR=1,NRMAX
         IF(MDLEQ0.EQ.0) THEN
            SSIN(NR,1)=SIE(NR)                            +SNB(NR)+SEX(NR,1)
            SSIN(NR,2)=PN(2)*SIE(NR)/(PN(2)+PN(3))-SNF(NR)+SNB(NR)+SEX(NR,2)
            SSIN(NR,3)=PN(3)*SIE(NR)/(PN(2)+PN(3))-SNF(NR)        +SEX(NR,3)
            SSIN(NR,4)=                            SNF(NR)        +SEX(NR,4)
            SSIN(NR,7)=-SIE(NR)                                   -SCX(NR)
            SSIN(NR,8)=                                    SNB(NR)+SCX(NR)
         ELSEIF(MDLEQ0.EQ.1) THEN
            SSIN(NR,1)=                                    SNB(NR)+SEX(NR,1)
            SSIN(NR,2)=                           -SNF(NR)+SNB(NR)+SEX(NR,2)
            SSIN(NR,3)=                           -SNF(NR)        +SEX(NR,3)
            SSIN(NR,4)=                            SNF(NR)        +SEX(NR,4)
            SSIN(NR,7)=0.D0
            SSIN(NR,8)=                                    SNB(NR)
         ENDIF
         PIN(NR,1)=PBCL(NR,1)+PFCL(NR,1)+PRF(NR,1) &
              &   +POH(NR)-PRSUM(NR)-PIE(NR)+PEX(NR,1)
         PIN(NR,2)=PBCL(NR,2)+PFCL(NR,2)+PRF(NR,2) &
              &   -PN(2)*PCX(NR)/(PN(2)+PN(3))+PEX(NR,2)
         PIN(NR,3)=PBCL(NR,3)+PFCL(NR,3)+PRF(NR,3) &
              &   -PN(3)*PCX(NR)/(PN(2)+PN(3))+PEX(NR,3)
         PIN(NR,4)=PBCL(NR,4)+PFCL(NR,4)+PRF(NR,4)+PEX(NR,4)
      ENDDO

      IF(RHOA.NE.1.D0) NRMAX=NRAMAX

      RETURN
      END SUBROUTINE TRCALC

!     ***********************************************************

!           RADIAL ELECTRIC FIELD

!     **********************************************************

      SUBROUTINE TRERAD

      USE TRCOMM, ONLY : AEE, AMM, BB, BP, DR, EPSRHO, ER, MDLER, NRMAX, &
           & PA, PADD, PBM, PNSS, PTS, PZ, QP, RHOG, RHOM, &
           & RJCB, RKEV, RN, RNF, RT, SUMPBM, VPOL, VTOR
      IMPLICIT NONE
      INTEGER(4):: NR
      REAL(8)   :: ALPHA_NEO, CS, DPD, DRL, EPS, F_UNTRAP, RHO_S, RLNI, &
           & RLTI, TEL, TERM_DP, TIL
      REAL(8)   :: DERIV3P

      IF(SUMPBM.EQ.0.D0) THEN
         PADD(1:NRMAX)=0.D0
      ELSE
         PADD(1:NRMAX)=PBM(1:NRMAX)*1.D-20/RKEV-RNF(1:NRMAX,1)*RT(1:NRMAX,2)
      ENDIF
      DO NR=1,NRMAX
         DRL=RJCB(NR)/DR
         IF(NR.EQ.NRMAX) THEN
            DPD = DERIV3P(PNSS(2)*PTS(2),RN(NR  ,2)*RT(NR  ,2)-PADD(NR  ), &
     &                    RN(NR-1,2)*RT(NR-1,2)-PADD(NR-1),RHOG(NR),RHOM(NR),RHOM(NR-1))
            TERM_DP = DPD*RKEV/(PZ(2)*AEE*PNSS(2))
         ELSE
            DPD =(  RN(NR+1,2)*RT(NR+1,2)-PADD(NR+1)-(RN(NR  ,2)*RT(NR  ,2)-PADD(NR  )))*DRL
            TERM_DP = DPD*RKEV/(PZ(2)*AEE*0.5D0*(RN(NR+1,2)+RN(NR,2)))
         ENDIF
         IF(MDLER.EQ.0) THEN
!     pressure gradient only
            ER(NR) = TERM_DP
         ELSEIF(MDLER.EQ.1) THEN
!     nabla p + toroidal rotation
            ER(NR) = TERM_DP+VTOR(NR)*BP(NR)
         ELSEIF(MDLER.EQ.2) THEN
!     nabla p + V_tor + poloidal rotation
            ER(NR) = TERM_DP+VTOR(NR)*BP(NR)-VPOL(NR)*BB
         ELSEIF(MDLER.EQ.3) THEN
!     Waltz definition
            EPS = EPSRHO(NR)
            F_UNTRAP = 1.D0-1.46D0*SQRT(EPS)+0.46D0*EPS**1.5D0
            ALPHA_NEO = 1.D0-0.8839D0*F_UNTRAP/(0.3477D0+0.4058D0*F_UNTRAP)
            IF(NR.EQ.NRMAX) THEN
               TEL = PTS(1)
               TIL = PTS(2)
               RLNI = -DERIV3P(PNSS(2),RN(NR,2),RN(NR-1,2),RHOG(NR),RHOM(NR),RHOM(NR-1))/PNSS(2)
               RLTI = -DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2),RHOG(NR),RHOM(NR),RHOM(NR-1))/PTS(2)
            ELSE
               TEL = 0.5D0*(RT(NR,1)+RT(NR+1,1))
               TIL = 0.5D0*(RT(NR,2)+RT(NR+1,2))
               RLNI = -(LOG(RN(NR+1,2))-LOG(RN(NR,2)))*DRL
               RLTI = -(LOG(RT(NR+1,2))-LOG(RT(NR,2)))*DRL
            ENDIF
            CS = SQRT(TEL*RKEV/(PA(2)*AMM))
            RHO_S = CS*PA(2)*AMM/(PZ(2)*AEE*BB)
            ER(NR) =-BB*( (TIL/TEL)*RHO_S*CS*(RLNI+ALPHA_NEO*RLTI)-EPS/QP(NR)*VTOR(NR))
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE TRERAD

!     ***********************************************************

!           CALCULATE Z-C

!     **********************************************************

      REAL(8) FUNCTION TRZEC(TE)

      IMPLICIT NONE
      REAL(8):: TE,TEL
      REAL(8):: BC00=2.093736D+03, BC01=5.153766D+03,BC02=5.042105D+03, BC03=2.445345D+03, BC04=5.879207D+02, BC05=5.609128D+01
      REAL(8):: BC10=7.286051D+00, BC11=4.506650D+01,BC12=1.595483D+02, BC13=2.112702D+02, BC14=1.186473D+02, BC15=2.400040D+01
      REAL(8):: BC20=5.998782D+00, BC21=6.808206D-03,BC22=-3.096242D-02,BC23=-4.794194D-02,BC24=1.374732D-01, BC25=5.157571D-01
      REAL(8):: BC30=5.988153D+00, BC31=7.335329D-02,BC32=-1.754858D-01,BC33=2.034126D-01, BC34=-1.145930D-01,BC35=2.517150D-02
      REAL(8):: BC40=-3.412166D+01,BC41=1.250454D+02,BC42=-1.550822D+02,BC43=9.568297D+01, BC44=-2.937297D+01,BC45=3.589667D+00


      TEL = LOG10(TE)
      IF(TE.LE.3.D-3) THEN
         TRZEC=0.D0
      ELSEIF(TE.LE.2.D-2) THEN
         TRZEC= BC00+(BC01*TEL)+(BC02*TEL**2)+(BC03*TEL**3)+(BC04*TEL**4)+(BC05*TEL**5)
      ELSEIF(TE.LE.0.2D0) THEN
         TRZEC= BC10+(BC11*TEL)+(BC12*TEL**2)+(BC13*TEL**3)+(BC14*TEL**4)+(BC15*TEL**5)
      ELSEIF(TE.LE.2.D0) THEN
         TRZEC= BC20+(BC21*TEL)+(BC22*TEL**2)+(BC23*TEL**3)+(BC24*TEL**4)+(BC25*TEL**5)
      ELSEIF(TE.LE.20.D0) THEN
         TRZEC= BC30+(BC31*TEL)+(BC32*TEL**2)+(BC33*TEL**3)+(BC34*TEL**4)+(BC35*TEL**5)
      ELSEIF(TE.LE.100.D0) THEN
         TRZEC= BC40+(BC41*TEL)+(BC42*TEL**2)+(BC43*TEL**3)+(BC44*TEL**4)+(BC45*TEL**5)
      ELSE
         TRZEC= 6.D0
      ENDIF

      RETURN
      END FUNCTION TRZEC

!     ***********************************************************

!           CALCULATE Z-FE

!     ***********************************************************

      REAL(8) FUNCTION TRZEFE(TE)

      IMPLICIT NONE
      REAL(8):: TE, TEL
      REAL(8):: BF10=8.778318D+00, BF11=-5.581412D+01,BF12=-1.225124D+02,BF13=-1.013985D+02,BF14=-3.914244D+01, BF15=-5.851445D+00
      REAL(8):: BF20=1.959496D+01, BF21=1.997852D+01, BF22=9.593373D+00, BF23=-7.609962D+01,BF24=-1.007190D+02, BF25=-1.144046D+01
      REAL(8):: BF30=1.478150D+01, BF31=6.945311D+01, BF32=-1.991202D+02,BF33=2.653804D+02, BF34=-1.631423D+02, BF35=3.770026D+01
      REAL(8):: BF40=2.122081D+01, BF41=6.891607D+00, BF42=-4.076853D+00,BF43=1.577171D+00, BF44=-5.139468D-01, BF45=8.934176D-02


      TEL = LOG10(TE)

      IF(TE.LE.3.D-3) THEN
         TRZEFE=0.D0
      ELSEIF(TE.LE.0.2D0) THEN
         TRZEFE = BF10+(BF11*TEL)+(BF12*TEL**2)+(BF13*TEL**3)+(BF14*TEL**4)+(BF15*TEL**5)
      ELSEIF(TE.LE.2.D0) THEN
         TRZEFE = BF20+(BF21*TEL)+(BF22*TEL**2)+(BF23*TEL**3)+(BF24*TEL**4)+(BF25*TEL**5)
      ELSEIF(TE.LE.20.D0) THEN
         TRZEFE = BF30+(BF31*TEL)+(BF32*TEL**2)+(BF33*TEL**3)+(BF34*TEL**4)+(BF35*TEL**5)
      ELSEIF(TE.LE.100.D0) THEN
         TRZEFE = BF40+(BF41*TEL)+(BF42*TEL**2)+(BF43*TEL**3)+(BF44*TEL**4)+(BF45*TEL**5)
      ELSE
         TRZEFE = 26.D0
      ENDIF

      RETURN
      END FUNCTION TRZEFE

!     ***********************************************************

!           CALCULATE Z-EFFF

!     ***********************************************************

      SUBROUTINE TRZEFF

      USE TRCOMM, ONLY : ANC, ANFE, MDLEQN, MDLUF, NRMAX, PZ, PZC, PZFE, RN, RNF, RT, ZEFF
      IMPLICIT NONE
      INTEGER(4):: NR
      REAL(8)   :: TE, TRZEC, TRZEFE

      IF(MDLUF.EQ.0) THEN
         DO NR=1,NRMAX
            TE =RT(NR,1)
            ZEFF(NR) =(PZ(2)  *PZ(2)  *RN(NR,2) +PZ(3)  *PZ(3)  *RN(NR,3) +PZ(4)  *PZ(4)  *RN(NR,4) &
     &                +TRZEC(TE)**2   *ANC (NR) +TRZEFE(TE)**2  *ANFE(NR))/RN(NR,1)
         ENDDO
      ELSE
         IF(MDLEQN.EQ.0) THEN ! fixed density
            DO NR=1,NRMAX
               ZEFF(NR) =(PZ(2)  *PZ(2)  *RN(NR,2) +PZ(3)  *PZ(3)  *RN(NR,3) +PZ(2)  *PZ(2)  *RNF(NR,1))/RN(NR,1)
            ENDDO
         ENDIF
      ENDIF

      DO NR=1,NRMAX
         TE=RT(NR,1)
         PZC(NR)=TRZEC(TE)
         PZFE(NR)=TRZEFE(TE)
      ENDDO
!
      RETURN
      END SUBROUTINE TRZEFF

!     ***********************************************************

!           RADIATION POWER - C

!     ***********************************************************

      REAL(8) FUNCTION TRRPC(TE)

      IMPLICIT NONE
      REAL(8):: TE,TEL,ARG
      REAL(8):: AC00=1.965300D+03, AC01=4.572039D+03, AC02=4.159590D+03, AC03=1.871560D+03, AC04=4.173889D+02, AC05=3.699382D+01
      REAL(8):: AC10=7.467599D+01, AC11=4.549038D+02, AC12=8.372937D+02, AC13=7.402515D+02, AC14=3.147607D+02, AC15=5.164578D+01
      REAL(8):: AC20=-2.120151D+01,AC21=-3.668933D-01,AC22=7.295099D-01, AC23=-1.944827D-01,AC24=-1.263576D-01,AC25=-1.491027D-01
      REAL(8):: AC30=-2.121979D+01,AC31=-2.346986D-01,AC32=4.093794D-01, AC33=7.874548D-02, AC34=-1.841379D-01,AC35=5.590744D-02
      REAL(8):: AC40=-2.476796D+01,AC41=9.408181D+00, AC42=-9.657446D+00,AC43=4.999161D+00, AC44=-1.237382D+00,AC45=1.160610D-01

      TEL = LOG10(TE)
      IF(TE.LE.3.D-3) THEN
         TEL=LOG10(3.D-3)
         ARG=AC00+(AC01*TEL)+(AC02*TEL**2)+(AC03*TEL**3)+(AC04*TEL**4)+(AC05*TEL**5)
      ELSEIF(TE.LE.2.D-2) THEN
         ARG=AC00+(AC01*TEL)+(AC02*TEL**2)+(AC03*TEL**3)+(AC04*TEL**4)+(AC05*TEL**5)
      ELSEIF(TE.LE.0.2D0) THEN
         ARG=AC10+(AC11*TEL)+(AC12*TEL**2)+(AC13*TEL**3)+(AC14*TEL**4)+(AC15*TEL**5)
      ELSEIF(TE.LE.2.D0) THEN
         ARG=AC20+(AC21*TEL)+(AC22*TEL**2)+(AC23*TEL**3)+(AC24*TEL**4)+(AC25*TEL**5)
      ELSEIF(TE.LE.20.D0) THEN
         ARG=AC30+(AC31*TEL)+(AC32*TEL**2)+(AC33*TEL**3)+(AC34*TEL**4)+(AC35*TEL**5)
      ELSEIF(TE.LE.1.D2) THEN
         ARG=AC40+(AC41*TEL)+(AC42*TEL**2)+(AC43*TEL**3)+(AC44*TEL**4)+(AC45*TEL**5)
      ENDIF
      TRRPC = 10.D0**(ARG-13.D0)

      RETURN
      END FUNCTION TRRPC

!     ***********************************************************

!           RADIATION POWER - FE

!     ***********************************************************

      REAL(8) FUNCTION TRRPFE(TE)

      IMPLICIT NONE
      REAL(8)::TE,TEL,ARG
      REAL(8)::AF10=-2.752599D+01,AF11=-3.908228D+01,AF12=-6.469423D+01,AF13=-5.555048D+01,AF14=-2.405568D+01,AF15=-4.093160D+00
      REAL(8)::AF20=-1.834973D+01,AF21=-1.252028D+00,AF22=-7.533115D+00,AF23=-3.289693D+00,AF24=2.866739D+01, AF25=2.830249D+01
      REAL(8)::AF30=-1.671042D+01,AF31=-1.646143D+01,AF32=3.766238D+01, AF33=-3.944080D+01,AF34=1.918529D+01, AF35=-3.509238D+00
      REAL(8)::AF40=-2.453957D+01,AF41=1.795222D+01, AF42=-2.356360D+01,AF43=1.484503D+01, AF44=-4.542323D+00,AF45=5.477462D-01


      TEL = LOG10(TE)
      IF(TE.LE.3.D-3) THEN
         TEL=LOG10(3.D-3)
         ARG=AF10+(AF11*TEL)+(AF12*TEL**2)+(AF13*TEL**3)+(AF14*TEL**4)+(AF15*TEL**5)
      ELSEIF(TE.LE.0.2D0) THEN
         ARG=AF10+(AF11*TEL)+(AF12*TEL**2)+(AF13*TEL**3)+(AF14*TEL**4)+(AF15*TEL**5)
      ELSEIF(TE.LE.2.D0) THEN
         ARG=AF20+(AF21*TEL)+(AF22*TEL**2)+(AF23*TEL**3)+(AF24*TEL**4)+(AF25*TEL**5)
      ELSEIF(TE.LE.20.D0) THEN
         ARG=AF30+(AF31*TEL)+(AF32*TEL**2)+(AF33*TEL**3)+(AF34*TEL**4)+(AF35*TEL**5)
      ELSEIF(TE.LE.1.D2) THEN
         ARG=AF40+(AF41*TEL)+(AF42*TEL**2)+(AF43*TEL**3)+(AF44*TEL**4)+(AF45*TEL**5)
      ENDIF
      TRRPFE = 10.D0**(ARG-13.D0)

      RETURN
      END FUNCTION TRRPFE

!     ***********************************************************

!           POWER LOSS

!     ***********************************************************

      SUBROUTINE TRLOSS

      USE TRCOMM, ONLY : AEE, ANC, ANFE, ANNU, DT, KUFDEV, MDLUF, MDLPR, &
           NRMAX, NT, NTUM, PCX, PIE, PRB, PRC, PRSUM, PRL, PRLU, RKEV, RN, &
           RT, SCX, SIE, TSCX, TSIE
      USE TRCOM1, ONLY : NTAMAX, NTXMAX, PNBI, TMU
      USE tr_cytran_mod, ONLY: tr_cytran
      IMPLICIT NONE
      INTEGER(4):: IERR, NR
      REAL(8):: ANDX, ANE, ANHE, ANT, EION, PLC, PLD, PLFE, PLHE, PLTT, PRLL, SCH, SION, TD, TE, TN, TNU, TRRPC, TRRPFE, TSL


      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         IF(NT.EQ.0) THEN
            TSL=DT*DBLE(1)
            DO NR=1,NRMAX
               PRSUM(NR)=PRLU(1,NR)
            ENDDO
         ELSE
            TSL=DT*DBLE(NT)
            IF(KUFDEV.EQ.'X') THEN
               DO NR=1,NRMAX
                  CALL TIMESPL(TSL,PRLL,TMU,PRLU(1,NR),NTXMAX,NTUM,IERR)
                  PRSUM(NR)=PRLL
               ENDDO
            ELSE
               DO NR=1,NRMAX
                  IF(PNBI.LT.12.D6) THEN
                     PRSUM(NR)=PRLU(1,NR)
                  ELSE
                     PRSUM(NR)=PRLU(2,NR)
                     IF(NT.EQ.NTAMAX) PRSUM(NR)=PRLU(3,NR)
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ELSEIF(MDLUF.EQ.2) THEN
         IF(NT.EQ.0) THEN
            TSL=DT*DBLE(1)
            DO NR=1,NRMAX
               PRSUM(NR)=PRLU(1,NR)
            ENDDO
         ELSE
            TSL=DT*DBLE(NT)
            DO NR=1,NRMAX
               PRSUM(NR)=PRLU(1,NR)
            ENDDO
         ENDIF
      ELSE
         DO NR=1,NRMAX
            ANE =RN(NR,1)
            ANDX=RN(NR,2)
            ANT =RN(NR,3)
            ANHE=RN(NR,4)
            TE  =RT(NR,1)
!     Radiation loss caused by impurities
            PLFE  = ANE*ANFE(NR)*TRRPFE(TE)*1.D40
            PLC   = ANE*ANC (NR)*TRRPC (TE)*1.D40
            PRL(NR)=PLFE+PLC
!     Bremsstrahlung
            PLD   = ANE*ANDX*5.35D-37*1.D0**2*SQRT(ABS(TE))*1.D40
            PLTT  = ANE*ANT *5.35D-37*1.D0**2*SQRT(ABS(TE))*1.D40
            PLHE  = ANE*ANHE*5.35D-37*2.D0**2*SQRT(ABS(TE))*1.D40
            PRB(NR)= PLD+PLTT+PLHE
!     Sumup
            SELECT CASE(MDLPR)
            CASE(0)
               PRSUM(NR)=PRL(NR)+PRB(NR)
            CASE(1,2)
               PRSUM(NR)=PRL(NR)+PRB(NR)+PRC(NR)
            END SELECT
         ENDDO
      ENDIF

!     ****** IONIZATION LOSS ******

      DO NR=1,NRMAX
         ANE=RN(NR,1)
         TE =RT(NR,1)
         EION  = 13.64D0
!         SIONS  = 1.76D-13/SQRT(TE*1.D3)*FKN(EION/(TE*1.D3))
         TN    = MAX(TE*1.D3/EION,1.D-2)
         SION  = 1.D-11*SQRT(TN)*EXP(-1.D0/TN)/(EION**1.5D0*(6.D0+TN))
         PIE(NR) = ANE*ANNU(NR)*SION*1.D40*EION*AEE
         SIE(NR) = ANE*ANNU(NR)*SION*1.D20
         TSIE(NR)= ANE         *SION*1.D20
      ENDDO

!     ****** CHARGE EXCHANGE LOSS ******

      DO NR=1,NRMAX
         ANE =RN(NR,1)
         ANDX=RN(NR,2)
         TE  =RT(NR,1)
         TD  =RT(NR,2)
         TNU = 0.D0

         TN    = LOG10(MAX(TD*1.D3,50.D0))
         SCH= 1.57D-16*SQRT(ABS(TD)*1.D3)*(TN*TN-14.63D0*TN+53.65D0)
         SCX(NR) =ANDX*ANNU(NR)*SCH*1.D20
         TSCX(NR)=ANDX         *SCH*1.D20

         PCX(NR)=1.5D0*ANDX*ANNU(NR)*SCH*(TD-TNU)*RKEV*1.D40
      ENDDO

      RETURN
      END SUBROUTINE TRLOSS

      REAL(8) FUNCTION FKN(X)

      IMPLICIT NONE
      REAL(8) X,GAMMA

      GAMMA=0.577215664901532D0
      FKN=-(LOG10(X)+GAMMA-X+X**2/4.D0-X**3/1.8D1)

      RETURN
      END FUNCTION FKN

!     ***********************************************

!         BOOTSTRAP CURRENT (NCLASS)

!     ***********************************************

      SUBROUTINE TRAJBS_NCLASS

      USE TRCOMM, ONLY : AJBS, AJBSNC, BB, CJBSP, CJBST, DR, NRM, NRMAX, NSMAX, PBSCD, PNSS, PTS, RG, RM, RN, RT
      IMPLICIT NONE
      INTEGER(4):: NR, NS, NSW
      REAL(8)   :: DRPNW, DRTNW, RPNW, RTNW, SUML
      REAL(8),DIMENSION(NRMAX)::  AJBSL
      REAL(8)   :: DERIV3P


      IF(PBSCD.LE.0.D0) RETURN

      NSW=1
      IF(NSW.EQ.0) THEN
         AJBSL(1:NRMAX)=PBSCD*AJBSNC(1:NRMAX)
      ELSE
         DO NR=1,NRMAX-1
            SUML=0.D0
            DO NS=1,NSMAX
               RTNW=0.5D0*(RT(NR+1,NS)+RT(NR,NS))
               RPNW=0.5D0*(RN(NR+1,NS)*RT(NR+1,NS)+RN(NR  ,NS)*RT(NR  ,NS))
               DRTNW=(RT(NR+1,NS)-RT(NR,NS))/DR
               DRPNW=(RN(NR+1,NS)*RT(NR+1,NS)-RN(NR,NS)*RT(NR,NS))/DR
               SUML=SUML+CJBST(NR,NS)*DRTNW/RTNW +CJBSP(NR,NS)*DRPNW/RPNW
            ENDDO
            AJBSL(NR)=-PBSCD*SUML/BB
         ENDDO

         NR=NRMAX
         SUML=0.D0
         DO NS=1,NSMAX
            RTNW=PTS(NS)
            RPNW=PNSS(NS)*PTS(NS)
            DRTNW=DERIV3P(PTS(NS),RT(NR,NS),RT(NR-1,NS),RG(NR),RM(NR),RM(NR-1))
            DRPNW=DERIV3P(PNSS(NS)*PTS(NS),RN(NR  ,NS)*RT(NR  ,NS),RN(NR-1,NS)*RT(NR-1,NS),RG(NR),RM(NR),RM(NR-1))
            SUML=SUML+CJBST(NR,NS)*DRTNW/RTNW +CJBSP(NR,NS)*DRPNW/RPNW
!CC            if(ns.eq.3) write(6,*) NR,PNSS(NS),PTS(NS)
         ENDDO
         AJBSL(NR)=-PBSCD*SUML/BB
      ENDIF
!
      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO

      RETURN
      END SUBROUTINE TRAJBS_NCLASS

!     ***********************************************

!         BOOTSTRAP CURRENT (O. SAUTER)

!     ***********************************************

      SUBROUTINE TRAJBSSAUTER

      USE TRCOMM, ONLY : AJBS, BB, DR, EPSRHO, MDLTPF, NRM, NRMAX, NSM, PADD, PBSCD, PNSS, PTS, PZ, QP, RDP, RHOG, RHOM, RKEV, &
     &                   RN, RPE, RR, RT, RW, TTRHOG, ZEFF
      IMPLICIT NONE
      INTEGER(4):: NR, NS
      REAL(8)   :: ANE, DPE, DPI, DRL, DTE, DTI, EPS, EPSS, F31TEFF, F32EETEFF, F32EITEFF, F34TEFF, FT, FTPF, PE, PPI, QL, RL31,&
     &             RL32, RL34, RLNLAME, RLNLAMII, RNM, RNP, RNTM, RNTP, RNUE, RNUI, RPIM, RPIP, SALFA, SALFA0, TE, TI, ZEFFL
      REAL(8),DIMENSION(NRMAX):: AJBSL, ANI
      REAL(8)   :: DERIV3P, F31, F32EE, F32EI


      IF(PBSCD.LE.0.D0) RETURN

      DO NR=1,NRMAX-1

         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         QL=ABS(QP(NR))
         ZEFFL=0.5D0*(ZEFF(NR)+ZEFF(NR+1))
         DRL=1.D0/DR

         RNTP=0.D0
         RNP =0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSM
            RNTP=RNTP+RN(NR+1,NS)*RT(NR+1,NS)
            RNP =RNP +RN(NR+1,NS)
            RNTM=RNTM+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR  ,NS)
         ENDDO
         RNTP=RNTP+RW(NR+1,1)+RW(NR+1,2)
         RNTM=RNTM+RW(NR  ,1)+RW(NR  ,2)
         RPIP=RNTP+PADD(NR+1)
         RPIM=RNTM+PADD(NR  )

!     ****** ION PARAMETER ******

!     *** ANI  is the the ion density (ni) ***
!     *** TI   is the ion temperature (Ti) ***
!     *** DTI  is the derivative of ion temperature (dTi/dr) ***
!     *** PPI  is the ion pressure (Pi) ***
!     *** DPI  is the derivative of ion pressure (dPi/dr) ***

         ANI(NR)=0.D0
         DO NS=2,NSM
            ANI(NR)=ANI(NR)+0.5D0*PZ(NS)*(RN(NR+1,NS)+RN(NR,NS))
         ENDDO
         PPI=0.5D0*(RPIP+RPIM)
         DPI=(RPIP-RPIM)*DRL
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         DTI=(RNTP/RNP-RNTM/RNM)*DRL

         rLnLamii=30.D0-LOG(PZ(2)**3*SQRT(ANI(NR)*1.D20)/(ABS(TI*1.D3)**1.5D0))
         RNUI=4.90D-18*QL*RR*ANI(NR)*1.D20*PZ(2)**4*rLnLamii /(ABS(TI*1.D3)**2*EPSS)
!
!     ****** ELECTORON PARAMETER ******

!     *** ANE  is the the electron density (ne) ***
!     *** TE   is the electron temperature (Te) ***
!     *** PE   is the electron pressure (Pe) ***
!     *** DTE  is the derivative of electron temperature (dTe/dr) ***
!     *** DPE  is the derivative of electron pressure (dPe/dr) ***

         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
         TE= 0.5D0*(RT(NR+1,1)+RT(NR,1))
         PE= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
         DTE=(RT(NR+1,1)-RT(NR,1))*DRL
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL

         rLnLame=31.3D0-LOG(SQRT(ANE*1.D20)/ABS(TE*1.D3))
         RNUE=6.921D-18*QL*RR*ANE*1.D20*ZEFFL*rLnLame /(ABS(TE*1.D3)**2*EPSS)

         RPE=PE/(PE+PPI)
         FT=FTPF(MDLTPF,EPS)

!         F33TEFF=FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE)+0.45D0*(1.D0-FT)*RNUE/ZEFFL**1.5)
         F31TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)   +0.5D0*(1.D0-FT)*RNUE/ZEFFL)
         F32EETEFF=FT/(1.D0+0.26D0*(1.D0-FT)*SQRT(RNUE)+0.18D0*(1.D0-0.37D0*FT)*RNUE/SQRT(ZEFFL))
         F32EITEFF=FT/(1.D0+(1.D0+0.6D0*FT)*SQRT(RNUE) +0.85D0*(1.D0-0.37D0*FT)*RNUE*(1.D0+ZEFFL))
         F34TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)   +0.5D0*(1.D0-0.5D0*FT)*RNUE/ZEFFL)

         SALFA0=-1.17D0*(1.D0-FT)/(1.D0-0.22D0*FT-0.19D0*FT**2)
         SALFA=((SALFA0+0.25D0*(1.D0-FT**2)*SQRT(RNUI))/(1.D0+0.5D0*SQRT(RNUI))+0.315D0*RNUI**2*FT**6) &
     &        /(1.D0+0.15D0*RNUI**2*FT**6)

!         RNZ=0.58D0+0.74D0/(0.76D0+ZEFFL)
!         SGMSPTZ=1.9012D4*(TE*1.D3)**1.5/(ZEFFL*RNZ*rLnLame)
!         SGMNEO=SGMSPTZ*F33(F33TEFF,ZEFFL)

         RL31=F31(F31TEFF,ZEFFL)
         RL32=F32EE(F32EETEFF,ZEFFL)+F32EI(F32EITEFF,ZEFFL)
         RL34=F31(F34TEFF,ZEFFL)

         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV *( RL31*(DPE/PE+DPI/PE)+RL32*DTE/TE &
     &             +RL34*SALFA*(1.D0-RPE)/RPE*DTI/TI)/RDP(NR)/BB
      ENDDO

      NR=NRMAX
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         QL=ABS(QP(NR))
         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
         DRL=1.D0/DR

!     In the following, we assume that
!        1. pressures of beam and fusion at rho=1 are negligible,
!        2. calibration of Pbeam (from UFILE) is not necessary at rho=1.

         RNTP=0.D0
         RNP =0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSM
            RNTP=RNTP+PNSS(NS)*PTS(NS)
            RNP =RNP +PNSS(NS)
            RNTM=RNTM+RN(NR-1,NS)*RT(NR-1,NS)+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR-1,NS)+RN(NR  ,NS)
         ENDDO
         RNTM=RNTM+RW(NR-1,1)+RW(NR-1,2)+RW(NR  ,1)+RW(NR  ,2)
         RPIP=RNTP
         RPIM=RNTM+PADD(NR-1)+PADD(NR  )
         RNTM=0.5D0*RNTM
         RNM =0.5D0*RNM
         RPIM=0.5D0*RPIM

!     ****** ION PARAMETER ******

!     *** ANI  is the the ion density (ni) ***
!     *** TI   is the ion temperature (Ti) ***
!     *** DTI  is the derivative of ion temperature (dTi/dr) ***
!     *** PPI  is the ion pressure (Pi) ***
!     *** DPI  is the derivative of ion pressure (dPi/dr) ***

         ANI(NR)=0.D0
         DO NS=2,NSM
            ANI(NR)=ANI(NR)+PZ(NS)*PNSS(NS)
         ENDDO
         PPI=RPIP
         DPI=(RPIP-RPIM)*DRL
         TI =RNTP/RNP
         DTI=(RNTP/RNP-RNTM/RNM)*DRL

         rLnLamii=30.D0-LOG(PZ(2)**3*SQRT(ANI(NR)*1.D20)/(ABS(TI*1.D3)**1.5D0))
         RNUI=4.90D-18*QL*RR*ANI(NR)*1.D20*PZ(2)**4*rLnLamii /(ABS(TI*1.D3)**2*EPSS)
!
!     ****** ELECTORON PARAMETER ******

!     *** ANE  is the the electron density (ne) ***
!     *** TE   is the electron temperature (Te) ***
!     *** PE   is the electron pressure (Pe) ***
!     *** DTE  is the derivative of electron temperature (dTe/dr) ***
!     *** DPE  is the derivative of electron pressure (dPe/dr) ***

         ANE=PNSS(1)
         TE =PTS(1)
         PE =PNSS(1)*PTS(1)
         DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPE=DERIV3P(PNSS(1)*PTS(1),RN(NR,1)*RT(NR,1),RN(NR-1,1)*RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
!
         rLnLame=31.3D0-LOG(SQRT(ANE*1.D20)/ABS(TE*1.D3))
         RNUE=6.921D-18*QL*RR*ANE*1.D20*ZEFFL*rLnLame /(ABS(TE*1.D3)**2*EPSS)
!
         RPE=PE/(PE+PPI)
         FT=FTPF(MDLTPF,EPS)

!     F33TEFF=FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE)+0.45D0*(1.D0-FT)*RNUE/ZEFFL**1.5)
         F31TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)   +0.5D0*(1.D0-FT)*RNUE/ZEFFL)
         F32EETEFF=FT/(1.D0+0.26D0*(1.D0-FT)*SQRT(RNUE)+0.18D0*(1.D0-0.37D0*FT)*RNUE/SQRT(ZEFFL))
         F32EITEFF=FT/(1.D0+(1.D0+0.6D0*FT)*SQRT(RNUE) +0.85D0*(1.D0-0.37D0*FT)*RNUE*(1.D0+ZEFFL))
         F34TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)   +0.5D0*(1.D0-0.5D0*FT)*RNUE/ZEFFL)

         SALFA0=-1.17D0*(1.D0-FT)/(1.D0-0.22D0*FT-0.19D0*FT**2)
         SALFA=((SALFA0+0.25D0*(1.D0-FT**2)*SQRT(RNUI)) /(1.D0+0.5D0*SQRT(RNUI))+0.315D0*RNUI**2*FT**6) &
     &        /(1.D0+0.15D0*RNUI**2*FT**6)

!     RNZ=0.58D0+0.74D0/(0.76D0+ZEFFL)
!     SGMSPTZ=1.9012D4*(TE*1.D3)**1.5/(ZEFFL*RNZ*rLnLame)
!     SGMNEO=SGMSPTZ*F33(F33TEFF,ZEFFL)
!
         RL31=F31(F31TEFF,ZEFFL)
         RL32=F32EE(F32EETEFF,ZEFFL)+F32EI(F32EITEFF,ZEFFL)
         RL34=F31(F34TEFF,ZEFFL)

         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV*( RL31*(DPE/PE+DPI/PE)+RL32*DTE/TE &
     &             +RL34*SALFA*(1.D0-RPE)/RPE*DTI/TI)/RDP(NR)/BB

      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO

      RETURN
      END SUBROUTINE TRAJBSSAUTER

!     *********************
!     *  Fitting Function *
!     *********************

      REAL(8) FUNCTION F33(X,Z)

      IMPLICIT NONE
      REAL(8) X,Z

      F33 = 1.D0-(1.D0+0.36D0/Z)*X+0.59D0/Z*X**2-0.23D0/Z*X**3

      RETURN
      END FUNCTION F33

      REAL(8) FUNCTION F31(X,Z)

      IMPLICIT NONE
      REAL(8) X,Z

      F31 = (1.D0+1.4D0/(Z+1.D0))*X-1.9D0/(Z+1.D0)*X**2 +0.3D0/(Z+1.D0)*X**3+0.2D0/(Z+1.D0)*X**4

      RETURN
      END FUNCTION F31

      REAL(8) FUNCTION F32EE(X,Z)

      IMPLICIT NONE
      REAL(8) X,Z

      F32EE = (0.05D0+0.62D0*Z)/(Z*(1.D0+0.44D0*Z))*(X-X**4)&
     &       +1.D0/(1.D0+0.22D0*Z)*(X**2-X**4-1.2D0*(X**3-X**4))+1.2D0/(1.D0+0.5D0*Z)*X**4

      RETURN
      END FUNCTION F32EE

      REAL(8) FUNCTION F32EI(X,Z)

      IMPLICIT NONE
      REAL(8) X,Z

      F32EI =-(0.56D0+1.93D0*Z)/(Z*(1.D0+0.44D0*Z))*(X-X**4) &
     &       +4.95D0/(1.D0+2.48D0*Z)*(X**2-X**4-0.55D0*(X**3-X**4))-1.2D0/(1.D0+0.5D0*Z)*X**4

      RETURN
      END FUNCTION F32EI

!     ************************************************

!         BOOTSTRAP CURRENT (Hirshman (Wilson))

!     ************************************************

      SUBROUTINE TRAJBSNEW

      USE TRCOMM, ONLY : AJBS, BB, DR, EPSRHO, NRM, NRMAX, NSM, PADD, PBSCD, PNSS, PTS, PZ, RDP, RHOG, RHOM, RKEV, RN, RT, &
     &                   RW, TTRHOG
      IMPLICIT NONE
      INTEGER(4):: NR, NS
      REAL(8) :: DDD, DDX, DPE, DPI, DRL, DTE, DTI, EPS, FT, PE, PPI, RL31, RL32, RNM, RNP, RNTM, RNTP, RPIM, RPIP, TE, TI
      REAL(8),DIMENSION(NRMAX)::  AJBSL, ANI
      REAL(8)   :: DERIV3P


      IF(PBSCD.LE.0.D0) RETURN

      DO NR=1,NRMAX-1

         EPS=EPSRHO(NR)
!         EPSS=SQRT(EPS)**3
!         QL=ABS(QP(NR))
!         ZEFFL=0.5D0*(ZEFF(NR)+ZEFF(NR+1))
         DRL=1.D0/DR

         RNTP=0.D0
         RNP= 0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSM
            RNTP=RNTP+RN(NR+1,NS)*RT(NR+1,NS)
            RNP =RNP +RN(NR+1,NS)
            RNTM=RNTM+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR  ,NS)
         ENDDO
         RNTP=RNTP+RW(NR+1,1)+RW(NR+1,2)
         RNTM=RNTM+RW(NR  ,1)+RW(NR  ,2)
         RPIP=RNTP+PADD(NR+1)
         RPIM=RNTM+PADD(NR  )

!     ****** ION PARAMETER ******

!     ***** ANI  is the the ion density (ni) *****
!     ***** TI   is the ion temperature (Ti) *****
!     ***** DTI  is the derivative of ion temperature (dTi/dr) *****
!     ***** PPI  is the ion pressure (Pi) *****
!     ***** DPI  is the derivative of ion pressure (dPi/dr) *****
!     ***** VTI  is the ion velocity (VTi) *****

         ANI(NR)=0.D0
         DO NS=2,NSM
            ANI(NR)=ANI(NR)+0.5D0*PZ(NS)*(RN(NR+1,NS)+RN(NR,NS))
         ENDDO
         PPI=0.5D0*(RPIP+RPIM)
         DPI=(RPIP-RPIM)*DRL
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         DTI=(RNTP/RNP-RNTM/RNM)*DRL
!         VTI=SQRT(ABS(TI)*RKEV/AMM)

!         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
!         rLnLam=17.3D0-DLOG(ANE)*0.5D0+DLOG(ABS(TI))*1.5D0
!         TAUI=12.D0*PI*SQRT(PI)*EPS0**2*SQRT(AMM)
!     &             *(ABS(TI)*RKEV)**1.5D0/(ANI(NR)*1.D20
!     &             *ZEFFL**4*AEE**4*rLnLam)

!         RNUI=QL*RR/(TAUI*VTI*EPSS)

!     ****** ELECTORON PARAMETER ******

!     ***** ANE  is the the electron density (ne) *****
!     ***** TE   is the electron temperature (Te) *****
!     ***** DTE  is the derivative of electron temperature (dTe/dr) ****
!     ***** PE   is the electron pressure (Pe) *****
!     ***** DPE  is the derivative of electron pressure (dPe/dr) *****
!     ***** VTE  is the electron velocity (VTe) *****

         TE= 0.5D0*(RT(NR+1,1)+RT(NR,1))
         PE= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
         DTE=(RT(NR+1,1)-RT(NR,1))*DRL
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL
!         VTE=SQRT(ABS(TE)*RKEV/AME)

!         rLnLam=15.2D0-DLOG(ANE)*0.5D0+DLOG(ABS(TE))
!         TAUE=6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)*(ABS(TE)*RKEV)**1.5D0/(ANI*1.D20*ZEFFL**2*AEE**4*rLnLam)

!         RNUE=QL*RR/(TAUE*VTE*EPSS)

!         FT=FTPF(MDLTPF,EPS)
         FT=(1.46D0*SQRT(EPS)+2.4D0*EPS)/(1.D0-EPS)**1.5D0

!         DDX=2.4D0+5.4D0*FT+2.6D0*FT**2
!         DDD=-1.17D0/(1.D0+0.46D0*FT)
!         C1=(4.D0+2.6D0*FT)/((1.D0+1.02D0*DSQRT(RNUE)+1.07D0*RNUE)*(1.D0+1.07D0*EPSS*RNUE))
!         C2=C1*TI/TE
!         C3=(7.D0+6.5D0*FT)/((1.D0+0.57D0*DSQRT(RNUE)+0.61D0*RNUE)*(1.D0+0.61D0*EPSS*RNUE))-2.5D0*C1
!         C4=((DDD+0.35D0*DSQRT(RNUI))/(1.D0+0.7D0*DSQRT(RNUI))
!     &      +2.1D0*EPS**3*RNUI**2)*C2/((1.D0-EPS**3*RNUI**2)*(1.D0+EPS**3*RNUE**2))

!         AJBSL(NR)=-PBSCD*FT*PE*1.D20*RKEV/(DDX*BP(NR))*(C1*(DPE/PE)
!     &         +C2*(DPI/PPI)+C3*(DTE/TE)+C4*(DTI/TI))

!     *** S. P. Hirshman, Phys Fluids 31, 1988 3150 ***
!     *** cited (H.R. WILSON, Nucl.Fusion 32,no.2,1992 259-263) ***

         DDX=1.414D0*PZ(2)+PZ(2)**2+FT*(0.754D0+2.657D0*PZ(2) &
     &        +2.D0*PZ(2)**2)+FT**2*(0.348D0+1.243D0*PZ(2)+PZ(2)**2)
         RL31= FT*(0.754D0+2.210D0*PZ(2)+PZ(2)**2 +FT*(0.348D0+1.243D0*PZ(2)+PZ(2)**2))/DDX
         RL32=-FT*(0.884D0+2.074D0*PZ(2))/DDX
         DDD=-1.172D0/(1.D0+0.462D0*FT)

         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV*(RL31*((DPE/PE)+(TI/(PZ(2)*TE)) &
     &        *((DPI/PPI)+DDD*(DTI/TI)))+RL32*(DTE/TE))/RDP(NR)/BB
      ENDDO

      NR=NRMAX
         EPS=EPSRHO(NR)
!         EPSS=SQRT(EPS)**3
!         QL=ABS(QP(NR))
!         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
         DRL=1.D0/DR

!     In the following, we assume that
!        1. pressures of beam and fusion at rho=1 are negligible,
!        2. calibration of Pbeam (from UFILE) is not necessary at rho=1.

         RNTP=0.D0
         RNP= 0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSM
            RNTP=RNTP+PNSS(NS)*PTS(NS)
            RNP =RNP +PNSS(NS)
            RNTM=RNTM+RN(NR-1,NS)*RT(NR-1,NS)+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR-1,NS)+RN(NR  ,NS)
         ENDDO
         RNTM=RNTP+RW(NR-1,1)+RW(NR-1,2)+RW(NR  ,1)+RW(NR  ,2)
         RPIP=RNTP
         RPIM=RNTM+PADD(NR-1)+PADD(NR  )
         RNTM=0.5D0*RNTM
         RNM =0.5D0*RNM
         RPIM=0.5D0*RPIM

!     ****** ION PARAMETER ******

!     ***** ANI  is the the ion density (ni) *****
!     ***** TI   is the ion temperature (Ti) *****
!     ***** DTI  is the derivative of ion temperature (dTi/dr) *****
!     ***** PPI  is the ion pressure (Pi) *****
!     ***** DPI  is the derivative of ion pressure (dPi/dr) *****
!     ***** VTI  is the ion velocity (VTi) *****

         ANI(NR)=0.D0
         DO NS=2,NSM
            ANI(NR)=ANI(NR)+PZ(NS)*PNSS(NS)
         ENDDO
         PPI=0.5D0*(RPIP+RPIM)
         DPI=(RPIP-RPIM)*DRL
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         DTI=(RNTP/RNP-RNTM/RNM)*DRL
!         VTI=SQRT(ABS(TI)*RKEV/AMM)

!         ANE=PNSS(1)
!         rLnLam=17.3D0-DLOG(ANE)*0.5D0+DLOG(ABS(TI))*1.5D0
!         TAUI=12.D0*PI*SQRT(PI)*EPS0**2*SQRT(AMM)
!     &             *(ABS(TI)*RKEV)**1.5D0/(ANI(NR)*1.D20
!     &             *ZEFFL**4*AEE**4*rLnLam)

!         RNUI=QL*RR/(TAUI*VTI*EPSS)

!     ****** ELECTORON PARAMETER ******

!     ***** ANE  is the the electron density (ne) *****
!     ***** TE   is the electron temperature (Te) *****
!     ***** DTE  is the derivative of electron temperature (dTe/dr) ****
!     ***** PE   is the electron pressure (Pe) *****
!     ***** DPE  is the derivative of electron pressure (dPe/dr) *****
!     ***** VTE  is the electron velocity (VTe) *****

         TE =PTS(1)
         PE =PNSS(1)*PTS(1)
         DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPE=DERIV3P(PNSS(1)*PTS(1),RN(NR,1)*RT(NR,1),RN(NR-1,1)*RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
!         VTE=SQRT(ABS(TE)*RKEV/AME)

!         rLnLam=15.2D0-DLOG(ANE)*0.5D0+DLOG(ABS(TE))
!         TAUE=6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)
!     &             *(ABS(TE)*RKEV)**1.5D0/(ANI*1.D20*ZEFFL**2*AEE**4*rLnLam)

!         RNUE=QL*RR/(TAUE*VTE*EPSS)

!         FT=FTPF(MDLTPF,EPS)
         FT=(1.46D0*SQRT(EPS)+2.4D0*EPS)/(1.D0-EPS)**1.5D0

!         DDX=2.4D0+5.4D0*FT+2.6D0*FT**2
!         DDD=-1.17D0/(1.D0+0.46D0*FT)
!         C1=(4.D0+2.6D0*FT)/((1.D0+1.02D0*DSQRT(RNUE)+1.07D0*RNUE)*(1.D0+1.07D0*EPSS*RNUE))
!         C2=C1*TI/TE
!         C3=(7.D0+6.5D0*FT)/((1.D0+0.57D0*DSQRT(RNUE)+0.61D0*RNUE)*(1.D0+0.61D0*EPSS*RNUE))
!     &      -2.5D0*C1
!         C4=((DDD+0.35D0*DSQRT(RNUI))/(1.D0+0.7D0*DSQRT(RNUI))+2.1D0*EPS**3*RNUI**2)*C2
!     &      /((1.D0-EPS**3*RNUI**2)*(1.D0+EPS**3*RNUE**2))

!         AJBSL(NR)=-PBSCD*FT*PE*1.D20*RKEV/(DDX*BP(NR))*(C1*(DPE/PE)
!     &         +C2*(DPI/PPI)+C3*(DTE/TE)+C4*(DTI/TI))

!     *** S. P. Hirshman, Phys Fluids 31, 1988 3150 ***
!     *** cited (H.R. WILSON, Nucl.Fusion 32,no.2,1992 259-263) ***

         DDX=1.414D0*PZ(2)+PZ(2)**2+FT*(0.754D0+2.657D0*PZ(2) &
     &        +2.D0*PZ(2)**2)+FT**2*(0.348D0+1.243D0*PZ(2)+PZ(2)**2)
         RL31=FT*( 0.754D0+2.21D0*PZ(2)+PZ(2)**2+FT*(0.348D0+1.243D0 &
     &            *PZ(2)+PZ(2)**2))/DDX
         RL32=-FT*(0.884D0+2.074D0*PZ(2))/DDX
         DDD=-1.172D0/(1.D0+0.462D0*FT)

         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV &
     &        *(RL31*((DPE/PE)+(TI/(PZ(2)*TE))*((DPI/PPI)+DDD*(DTI/TI)))+RL32*(DTE/TE))/RDP(NR)/BB

      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO

      RETURN
      END SUBROUTINE TRAJBSNEW

!     ***********************************************************

!           BOOTSTRAP CURRENT (Hinton & Hazeltine)

!     ***********************************************************

      SUBROUTINE TRAJBS

      USE TRCOMM, ONLY : AJBS, AME, AMM, BB, BP, DR, EPSRHO, NRM, NRMAX, NSMAX, PA, PBSCD, PNSS, PTS, PZ, QP, RHOG, RHOM, &
     &                   RJCB, RKEV, RN, RR, RT, ZEFF
      IMPLICIT NONE
      INTEGER(4):: NR
      REAL(8)   :: A, AMA, AMD, AMT, ANA, ANDX, ANE, ANT, BPL, DPA, DPD, DPE, DPT, DRL, DTA, DTD, DTE, DTT, EPS, EPSS, FACT, &
     &             FTAUE, FTAUI, H, PAL, PDL, PEL, PTL, RK13E, RK23E, RK3A, RK3D, RK3T, RNUA, RNUD, RNUE, RNUT, TAL, TAUA,   &
     &             TAUD, TAUE, TAUT, TDL, TEL, TTL, VTA, VTD, VTE, VTT, ZEFFL
      REAL(8)   :: DERIV3P
      REAL(8),DIMENSION(NRMAX):: AJBSL

      REAL(8):: RK13=2.30D0, RA13=1.02D0, RB13=1.07D0, RC13=1.07D0, RK23=4.19D0, RA23=0.57D0, RB23=0.61D0, RC23=0.61D0

!     ZEFF=1

!      DATA RK11,RA11,RB11,RC11/1.04D0,2.01D0,1.53D0,0.89D0/
!      DATA RK12,RA12,RB12,RC12/1.20D0,0.76D0,0.67D0,0.56D0/
!      DATA RK22,RA22,RB22,RC22/2.55D0,0.45D0,0.43D0,0.43D0/
!!      DATA RK13,RA13,RB13,RC13/2.30D0,1.02D0,1.07D0,1.07D0/
!!      DATA RK23,RA23,RB23,RC23/4.19D0,0.57D0,0.61D0,0.61D0/
!      DATA RK33,RA33,RB33,RC33/1.83D0,0.68D0,0.32D0,0.66D0/
!      DATA RK2 ,RA2 ,RB2 ,RC2 /0.66D0,1.03D0,0.31D0,0.74D0/

      IF(PBSCD.LE.0.D0) RETURN

      AMD=PA(2)*AMM
      AMT=PA(3)*AMM
      AMA=PA(4)*AMM

      DO NR=1,NRMAX-1

         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
         ANDX=0.5D0*(RN(NR+1,2)+RN(NR,2))
         ANT =0.5D0*(RN(NR+1,3)+RN(NR,3))
         ANA =0.5D0*(RN(NR+1,4)+RN(NR,4))
         TEL=ABS(0.5D0*(RT(NR+1,1)+RT(NR,1)))
         TDL=ABS(0.5D0*(RT(NR+1,2)+RT(NR,2)))
         TTL=ABS(0.5D0*(RT(NR+1,3)+RT(NR,3)))
         TAL=ABS(0.5D0*(RT(NR+1,4)+RT(NR,4)))
         PEL=0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
         PDL=0.5D0*(RN(NR+1,2)*RT(NR+1,2)+RN(NR,2)*RT(NR,2))
         PTL=0.5D0*(RN(NR+1,3)*RT(NR+1,3)+RN(NR,3)*RT(NR,3))
         PAL=0.5D0*(RN(NR+1,4)*RT(NR+1,4)+RN(NR,4)*RT(NR,4))
         ZEFFL=0.5D0*(ZEFF(NR+1)+ZEFF(NR))

         TAUE = FTAUE(ANE,ANDX,TEL,ZEFFL)
         TAUD = FTAUI(ANE,ANDX,TDL,PZ(2),PA(2))
         TAUT = FTAUI(ANE,ANT ,TTL,PZ(3),PA(3))
         TAUA = FTAUI(ANE,ANA ,TAL,PZ(4),PA(4))

         VTE=SQRT(TEL*RKEV/AME)
         VTD=SQRT(TDL*RKEV/AMD)
         VTT=SQRT(TTL*RKEV/AMT)
         VTA=SQRT(TAL*RKEV/AMA)

         RNUE=ABS(QP(NR))*RR/(TAUE*VTE*EPSS)
         RNUD=ABS(QP(NR))*RR/(TAUD*VTD*EPSS)
         RNUT=ABS(QP(NR))*RR/(TAUT*VTT*EPSS)
         RNUA=ABS(QP(NR))*RR/(TAUA*VTA*EPSS)

!         RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)+(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
!         RK12E=RK12*(1.D0/(1.D0+RA12*SQRT(RNUE)+RB12*RNUE)+(EPSS*RC12)**2/RB12*RNUE/(1.D0+RC12*RNUE*EPSS))
!         RK22E=RK22*(1.D0/(1.D0+RA22*SQRT(RNUE)+RB22*RNUE)+(EPSS*RC22)**2/RB22*RNUE/(1.D0+RC22*RNUE*EPSS))
         RK13E=RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)/(1.D0+RC13*RNUE*EPSS)
         RK23E=RK23/(1.D0+RA23*SQRT(RNUE)+RB23*RNUE)/(1.D0+RC23*RNUE*EPSS)
!         RK33E=RK33/(1.D0+RA33*SQRT(RNUE)+RB33*RNUE)/(1.D0+RC33*RNUE*EPSS)

!         RK2D =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUD)+RB2 *RNUD)+(EPSS*RC2 )**2/RB2 *RNUD/(1.D0+RC2 *RNUD*EPSS))
!         RK2T =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUT)+RB2 *RNUT)+(EPSS*RC2 )**2/RB2 *RNUT/(1.D0+RC2 *RNUT*EPSS))
!         RK2A =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUA)+RB2 *RNUA)+(EPSS*RC2 )**2/RB2 *RNUA/(1.D0+RC2 *RNUA*EPSS))
         RK3D=((1.17D0-0.35D0*SQRT(RNUD))/(1.D0+0.7D0*SQRT(RNUD))-2.1D0*(RNUD*EPSS)**2)/(1.D0+(RNUD*EPSS)**2)
         RK3T=((1.17D0-0.35D0*SQRT(RNUT))/(1.D0+0.7D0*SQRT(RNUT))-2.1D0*(RNUT*EPSS)**2)/(1.D0+(RNUT*EPSS)**2)
         RK3A=((1.17D0-0.35D0*SQRT(RNUA))/(1.D0+0.7D0*SQRT(RNUA))-2.1D0*(RNUA*EPSS)**2)/(1.D0+(RNUA*EPSS)**2)

         DRL=RJCB(NR)/DR
         DTE=(RT(NR+1,1)-RT(NR,1))*DRL
         DTD=(RT(NR+1,2)-RT(NR,2))*DRL
         DTT=(RT(NR+1,3)-RT(NR,3))*DRL
         DTA=(RT(NR+1,4)-RT(NR,4))*DRL
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL
         DPD=(RN(NR+1,2)*RT(NR+1,2)-RN(NR,2)*RT(NR,2))*DRL
         DPT=(RN(NR+1,3)*RT(NR+1,3)-RN(NR,3)*RT(NR,3))*DRL
         DPA=(RN(NR+1,4)*RT(NR+1,4)-RN(NR,4)*RT(NR,4))*DRL
         BPL=BP(NR)

         FACT=1.D0/(1.D0+(RNUE*EPS)**2)
         IF(NSMAX.EQ.2) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
         ELSEIF(NSMAX.EQ.3) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL &
     &       +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)+(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)
         ELSEIF(NSMAX.EQ.4) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL) &
     &       +(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)+(PAL/PEL)*(DPA/PAL-RK3A*FACT*DTA/TAL)
         ENDIF
         H=BB/(BB+BP(NR))
         AJBSL(NR)=-PBSCD*H*SQRT(EPS)*ANE*1.D20*TEL*RKEV/BPL *(RK13E*A+RK23E*DTE/TEL)
!
!        WRITE(6,'(I3,1P6E12.4)') NR,RNUE,RK13E,RK23E,A,BPL,AJBSL(NR)
!         IF(NR.EQ.NRMAX) THEN
!            WRITE(6,'(1P6E12.4)') DN,DTE,ANE,PNSS(1),TE,PTS(1)
!            WRITE(6,'(1P6E12.4)') RK13E*(RN(NR,1)*RT(NR,1)
!     &                   +RN(NR,2)*RT(NR,2)
!     &                   +RN(NR,3)*RT(NR,3)
!     &                   +RN(NR,4)*RT(NR,4)
!     &                   +RW(NR,1)
!     &                   +RW(NR,2)         )*DN/ANE,
!     &            +(RK23E-1.5D0*RK13E)*RN(NR,1)*DTE,
!     &            +RK13D*(RK23D-1.5D0)*RN(NR,2)*DTD,
!     &            BPL,AJBS(NR-1),AJBS(NR)
!         ENDIF
      ENDDO

      NR=NRMAX
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         ANE=PNSS(1)
         ANDX=PNSS(2)
!         ANT =PNSS(3)
!         ANA =PNSS(4)
         TEL=ABS(PTS(1))
         TDL=ABS(PTS(2))
         TTL=ABS(PTS(3))
         TAL=ABS(PTS(4))
         PEL=PNSS(1)*PTS(1)
         PDL=PNSS(2)*PTS(2)
         PTL=PNSS(3)*PTS(3)
         PAL=PNSS(4)*PTS(4)
         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)

         TAUE = FTAUE(ANE,ANDX,TEL,ZEFFL)
         TAUD = FTAUI(ANE,ANDX,TDL,PZ(2),PA(2))
         TAUT = FTAUI(ANE,ANT ,TTL,PZ(3),PA(3))
         TAUA = FTAUI(ANE,ANA ,TAL,PZ(4),PA(4))

         VTE=SQRT(TEL*RKEV/AME)
         VTD=SQRT(TDL*RKEV/AMD)
         VTT=SQRT(TTL*RKEV/AMT)
         VTA=SQRT(TAL*RKEV/AMA)

         RNUE=ABS(QP(NR))*RR/(TAUE*VTE*EPSS)
         RNUD=ABS(QP(NR))*RR/(TAUD*VTD*EPSS)
         RNUT=ABS(QP(NR))*RR/(TAUT*VTT*EPSS)
         RNUA=ABS(QP(NR))*RR/(TAUA*VTA*EPSS)

!         RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)+(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
!         RK12E=RK12*(1.D0/(1.D0+RA12*SQRT(RNUE)+RB12*RNUE)+(EPSS*RC12)**2/RB12*RNUE/(1.D0+RC12*RNUE*EPSS))
!         RK22E=RK22*(1.D0/(1.D0+RA22*SQRT(RNUE)+RB22*RNUE)+(EPSS*RC22)**2/RB22*RNUE/(1.D0+RC22*RNUE*EPSS))
         RK13E=RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE) /(1.D0+RC13*RNUE*EPSS)
         RK23E=RK23/(1.D0+RA23*SQRT(RNUE)+RB23*RNUE) /(1.D0+RC23*RNUE*EPSS)
!         RK33E=RK33/(1.D0+RA33*SQRT(RNUE)+RB33*RNUE)/(1.D0+RC33*RNUE*EPSS)

!         RK2D =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUD)+RB2 *RNUD)+(EPSS*RC2 )**2/RB2 *RNUD/(1.D0+RC2 *RNUD*EPSS))
!         RK2T =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUT)+RB2 *RNUT)+(EPSS*RC2 )**2/RB2 *RNUT/(1.D0+RC2 *RNUT*EPSS))
!         RK2A =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUA)+RB2 *RNUA)+(EPSS*RC2 )**2/RB2 *RNUA/(1.D0+RC2 *RNUA*EPSS))
         RK3D=((1.17D0-0.35D0*SQRT(RNUD))/(1.D0+0.7D0*SQRT(RNUD))-2.1D0*(RNUD*EPSS)**2)/(1.D0+(RNUD*EPSS)**2)
         RK3T=((1.17D0-0.35D0*SQRT(RNUT))/(1.D0+0.7D0*SQRT(RNUT))-2.1D0*(RNUT*EPSS)**2)/(1.D0+(RNUT*EPSS)**2)
         RK3A=((1.17D0-0.35D0*SQRT(RNUA))/(1.D0+0.7D0*SQRT(RNUA))-2.1D0*(RNUA*EPSS)**2)/(1.D0+(RNUA*EPSS)**2)

         DRL=RJCB(NR)/DR
         DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DTD=DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DTT=DERIV3P(PTS(3),RT(NR,3),RT(NR-1,3),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DTA=DERIV3P(PTS(4),RT(NR,4),RT(NR-1,4),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPE=DERIV3P(PNSS(1)*PTS(1),RN(NR,1)*RT(NR,1),RN(NR-1,1)*RT(NR-1,1),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPD=DERIV3P(PNSS(2)*PTS(2),RN(NR,2)*RT(NR,2),RN(NR-1,2)*RT(NR-1,2),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPT=DERIV3P(PNSS(3)*PTS(3),RN(NR,3)*RT(NR,3),RN(NR-1,3)*RT(NR-1,3),RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPA=DERIV3P(PNSS(4)*PTS(4),RN(NR,4)*RT(NR,4),RN(NR-1,4)*RT(NR-1,4),RHOG(NR),RHOM(NR),RHOM(NR-1))
         BPL=BP(NR)

         FACT=1.D0/(1.D0+(RNUE*EPS)**2)
         IF(NSMAX.EQ.2) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
         ELSEIF(NSMAX.EQ.3) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL &
     &       +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)+(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)
         ELSEIF(NSMAX.EQ.4) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL) &
     &       +(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)+(PAL/PEL)*(DPA/PAL-RK3A*FACT*DTA/TAL)
         ENDIF
         H=BB/(BB+BP(NR))
         AJBSL(NR)=-PBSCD*H*SQRT(EPS)*ANE*1.D20*TEL*RKEV/BPL*(RK13E*A+RK23E*DTE/TEL)
!
!        WRITE(6,'(I3,1P6E12.4)') NR,RNUE,RK13E,RK23E,A,BPL,AJBSL(NR)
!         IF(NR.EQ.NRMAX) THEN
!            WRITE(6,'(1P6E12.4)') DN,DTE,ANE,PNSS(1),TE,PTS(1)
!            WRITE(6,'(1P6E12.4)') RK13E*(RN(NR,1)*RT(NR,1)
!     &                   +RN(NR,2)*RT(NR,2)
!     &                   +RN(NR,3)*RT(NR,3)
!     &                   +RN(NR,4)*RT(NR,4)
!     &                   +RW(NR,1)
!     &                   +RW(NR,2)         )*DN/ANE,
!     &            +(RK23E-1.5D0*RK13E)*RN(NR,1)*DTE,
!     &            +RK13D*(RK23D-1.5D0)*RN(NR,2)*DTD,
!     &            BPL,AJBS(NR-1),AJBS(NR)
!         ENDIF

      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO

      RETURN
      END SUBROUTINE TRAJBS

!     ***********************************************************

!           CALCULATE AJ, AJOH, POH, EZOH

!     ***********************************************************

      SUBROUTINE TRAJOH

      USE TRCOMM, ONLY : ABRHOG, AJ, AJBS, AJNB, AJOH, AJRF, AJTOR, BB, DR, DVRHO, DVRHOG, ETA, EZOH, KUFDEV, MDLEQB, &
     &                   MDLJQ, MDLUF, NRMAX, POH, RDP, RMU0, RR, TTRHO, TTRHOG, rm, ABVRHOG, RDPVRHOG, PI, rm, rg
      IMPLICIT NONE
      INTEGER(4):: NR
      REAL(8)   :: FACTOR0, FACTORM, FACTORP

      IF(MDLEQB.EQ.1.OR.MDLJQ.EQ.1.OR.(MDLUF.EQ.0.OR.MDLUF.EQ.3)) THEN
      NR=1
         FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
         FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
         AJ(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
      DO NR=2,NRMAX
         FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
         FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
         FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
         AJ(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
      ENDDO
      NR=1
         FACTOR0=RR/(RMU0*DVRHO(NR))
         FACTORP=ABVRHOG(NR  )
         AJTOR(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
      DO NR=2,NRMAX
         FACTOR0=RR/(RMU0*DVRHO(NR))
         FACTORM=ABVRHOG(NR-1)
         FACTORP=ABVRHOG(NR  )
         AJTOR(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
      ENDDO
      ENDIF

      AJOH(1:NRMAX) = AJ(1:NRMAX)-(AJNB(1:NRMAX)+AJRF(1:NRMAX)+AJBS(1:NRMAX))
      EZOH(1:NRMAX) = ETA(1:NRMAX)*AJOH(1:NRMAX)
!!!      IF(KUFDEV.EQ.'lhd') THEN
!!!         POH(1:NRMAX)  = 0.D0
!!!      ELSE
         POH(1:NRMAX)  = EZOH(1:NRMAX)*AJOH(1:NRMAX)
!!!      ENDIF

      RETURN
      END SUBROUTINE TRAJOH

!     ***********************************************************

!           SAWTOOTH OSCILLATION

!     ***********************************************************

      SUBROUTINE TRSAWT

      USE TRCOMM, ONLY : AR1RHOG, ARRHOG, BP, DR, DVRHO, DVRHOG, MDLST, NRM, NRMAX, NSM, PI, QP, RDP, RG, RM, RN, RR, RT, &
     &                   T, TTRHOG, RDPVRHOG
      IMPLICIT NONE
      INTEGER(4):: IONE, IZEROX, LN, LQ, LT, NR, NS
      REAL(8)   :: RNN, RTN, SUML, SUML1, SUML2
      REAL(8),DIMENSION(NRMAX):: QONE


      IF(MDLST.EQ.0) RETURN

      IF(MOD(MDLST,2).EQ.1) THEN
         LT = 1
      ELSE
         LT = 0
      ENDIF
      IF(MOD(MDLST/2,2).EQ.1) THEN
         LN = 1
      ELSE
         LN = 0
      ENDIF
      IF(MOD(MDLST/4,2).EQ.1) THEN
         LQ = 1
      ELSE
         LQ = 0
      ENDIF

      WRITE(6,601) MDLST,T
  601 FORMAT(' ','# SAWTOOTH OSCILLATION -TYPE ',I1, ' AT ',F7.3,' SEC')

      SUML=0.D0
      DO NR=1,NRMAX
         SUML = SUML+(1.D0/ABS(QP(NR))-1.D0)*(DVRHOG(NR)/2.D0*PI*RR)*DR
         IF(SUML.LT.0.D0) GOTO 1000
      ENDDO
      NR=NRMAX

 1000 IZEROX=NR

      DO NR=1,NRMAX
         IF(QP(NR).GE.1.D0) GOTO 2000
      ENDDO
      NR=NRMAX

 2000 IONE=NR

      SUML = SUM((1.D0/ABS(QP(IONE:IZEROX))-1.D0)*DVRHOG(IONE:IZEROX))/(2.D0*PI*RR)*DR

      QONE(1:IZEROX) = 1.D0+SUML*4.D0*RG(1:IZEROX)**2/RG(IZEROX)**4

      IF(LT.EQ.1) THEN
         DO NS=1,NSM
            SUML1 = SUM(RN(1:IZEROX,NS)                *DVRHO(1:IZEROX))
            SUML2 = SUM(RN(1:IZEROX,NS)*RT(1:IZEROX,NS)*DVRHO(1:IZEROX))
            RTN = SUML2/SUML1
            RT(1:IZEROX,NS) = RTN
         ENDDO
      ENDIF

      IF(LN.EQ.1) THEN
         DO NS=1,NSM
            SUML1 = SUM(                DVRHO(1:IZEROX))
            SUML2 = SUM(RN(1:IZEROX,NS)*DVRHO(1:IZEROX))
            RNN = SUML2/SUML1
            RN(1:IZEROX,NS) = RNN
         ENDDO
      ENDIF

      IF(LQ.EQ.1) THEN
         QP(1:IZEROX) = 1.D0/QONE(1:IZEROX)
         RDPVRHOG(1:IZEROX)=TTRHOG(1:IZEROX)*ARRHOG(1:IZEROX) &
     &                /(4.D0*PI**2*QP(1:IZEROX))
         RDP(1:IZEROX)=RDPVRHOG(1:IZEROX)*DVRHOG(1:IZEROX)
         BP(1:IZEROX) =AR1RHOG(1:IZEROX)*RDP(1:IZEROX)/RR
      ENDIF

      WRITE(6,602) RM(IONE),RM(IZEROX),RTN,RNN
602   FORMAT(' ',' R-ONE,R-ZERO,RTN,RNN = ',4F8.3)

      DO NR=1,IZEROX+2
         WRITE(6,'(A,I5,1P5E12.4)') 'NR,R,Q,B,Q,S=', &
                       NR,RM(NR),QP(NR),BP(NR),QONE(NR),SUML
      ENDDO
      RETURN
      END SUBROUTINE TRSAWT

!     ***********************************************************

!           TRAPPED PARTICLE FRACTION

!     ***********************************************************

      REAL(8) FUNCTION FTPF(ID,EPS)

      IMPLICIT NONE
      INTEGER(4)                  :: ID, IERR, N
      INTEGER(4),PARAMETER        :: IMAX=20
      REAL(8)                     :: EPS, EPSC, FTLL, FTUL, OMEGA, PI, S
      REAL(8),DIMENSION(IMAX,IMAX):: TABLE
      EXTERNAL FTL, FTU

      IF(ID.EQ.1) THEN
!  Y. R. Lin-Liu and R. L. Miller, PoP 2 1666 (1995), eqs(7)(13)(18)(19)
         PI=3.14159265358979323846D0
         EPSC=1.D-9
         FTUL=1.D0-(1.D0-1.5D0*SQRT(EPS)+0.5D0*EPS**1.5D0)/SQRT(1-EPS**2)
         CALL RMBRG(0.D0,2.D0*PI,EPSC,S,IMAX,N,IERR,TABLE,EPS,FTL)
         FTLL=1.D0-(1.D0-EPS)**1.5D0/SQRT(1.D0+EPS)*(S/(2.D0*PI))
         OMEGA=(3.D0*SQRT(2.D0)/2.D0*0.69D0-3.D0*SQRT(2.D0)/PI)/(1.5D0-3.D0*SQRT(2.D0)/PI)
         FTPF=OMEGA*FTUL+(1.D0-OMEGA)*FTLL
      ELSEIF(ID.EQ.2) THEN
!  S. P. Hirshman et al., NF 17 611 (1977)
         FTPF=1.D0-(1.D0-EPS)**2.D0/(DSQRT(1.D0-EPS**2)*(1.D0+1.46D0*DSQRT(EPS)))
      ELSEIF(ID.EQ.3) THEN
!  Y. R. Lin-Liu and R. L. Miller, PoP 2 1666 (1995), eqs(16)(17)(18)
         PI=3.14159265358979323846D0
         FTUL=1.5D0*SQRT(EPS)
         FTLL=3.D0*SQRT(2.D0)/PI*SQRT(EPS)
         FTPF=0.75D0*FTUL+0.25D0*FTLL
      ELSEIF(ID.EQ.4) THEN
!  M. N. Rosenbluth et al., PoF 15 116 (1972)
         FTPF=1.46D0*SQRT(EPS)
      ELSE
!  Y. B. Kim et al., PoF B 3 2050 (1991) eq(C18), default
         FTPF=1.46D0*SQRT(EPS)-0.46D0*(EPS)**1.5D0
      ENDIF

      RETURN
      END FUNCTION FTPF

!     *********************************************

!           FUNCTION FOR ROMBERG INTEGRATION

!     *********************************************

      REAL(8) FUNCTION FTU(X,EPS)

      IMPLICIT NONE
      REAL(8):: EPS, X

      FTU = X/SQRT(1.D0-X*(1.D0-EPS))

      RETURN
      END FUNCTION FTU

      REAL(8) FUNCTION FTL(X,EPS)
!
      IMPLICIT NONE
      REAL(8):: X, EPS
      REAL(8):: H

      H = (1.D0 - EPS) / (1.D0 + EPS * COS(X))
      FTL = (1.D0 - SQRT(1.D0 - H) * (1.D0 + 0.5D0 * H)) / H**2

      RETURN
      END FUNCTION FTL
!
!     *********************************************

!           ROMBERG INTEGRATION METHOD

!     *********************************************

      SUBROUTINE RMBRG(A,B,EPS,S,IMAX,N,IERR,T,ARG,F)

!     <input>
!        A     : lower bound
!        B     : upper bound
!        EPS   : stopping criterion
!        IMAX  : maximum division number
!        F     : formula of integrand
!     <output>
!        S     : integration value
!        N     : division number
!        IERR  : error indicator
!        T     : Romberg T table

      IMPLICIT NONE
      REAL(8),INTENT(IN)      ::  A, B, ARG, EPS
      REAL(8),INTENT(OUT)     ::  S
      INTEGER(4),INTENT(INOUT):: IMAX
      INTEGER(4),INTENT(OUT)  :: IERR, N
      REAL(8),DIMENSION(IMAX,IMAX),INTENT(OUT)::  T
      INTEGER(4)::  I, J, K, N2
      REAL(8)   ::  F, X, H, S1, Y1, Y2

      EXTERNAL F

      DO K=1,IMAX
         N=2**(K-1)
         N2=N/2
         H=(B-A)/N
         Y1=0
         IF(N.EQ.1) THEN
            Y2=(F(A,ARG)+F(B,ARG))/2
         ELSE
            DO I=1,N2
               X=A+(2*I-1)*H
               Y1=Y1+F(X,ARG)
            ENDDO
            Y2=Y2+Y1
            Y1=0.D0
         ENDIF
         S=H*Y2
         T(K,1)=S
         IF(K.LE.1) GOTO 10
         DO J=2,K
            T(K,J)=T(K,J-1)+(T(K,J-1)-T(K-1,J-1))/(4**(J-1)-1)
         ENDDO
         S=T(K,K)
         S1=T(K,K-1)
         IF(ABS(S-S1).LT.EPS) THEN
            IERR=0
            IMAX=K
            RETURN
         ENDIF
 10      CONTINUE
      ENDDO
      IERR=1

      RETURN
      END SUBROUTINE RMBRG

!     ***********************************************************

!           COULOMB LOGARITHM

!     ***********************************************************

      REAL(8) FUNCTION COULOG(NS1,NS2,ANEL,TL)

!     ANEL : electron density [10^20 /m^3]
!     TL   : electron or ion temperature [keV]
!            in case of ion-ion collision, TL becomes ion temp.

      IMPLICIT NONE
      INTEGER(4):: NS1,NS2
      REAL(8)   :: ANEL,TL

      IF(NS1.EQ.1.AND.NS2.EQ.1) THEN
         COULOG=14.9D0-0.5D0*LOG(ANEL)+LOG(TL)
      ELSE
         IF(NS1.EQ.1.OR.NS2.EQ.1) THEN
            COULOG=15.2D0-0.5D0*LOG(ANEL)+LOG(TL)
         ELSE
            COULOG=17.3D0-0.5D0*LOG(ANEL)+1.5D0*LOG(TL)
         ENDIF
      ENDIF

      RETURN
      END FUNCTION COULOG

!     ***********************************************************

!           COLLISION TIME

!     ***********************************************************

!     between electrons and ions

      REAL(8) FUNCTION FTAUE(ANEL,ANIL,TEL,ZL)

!     ANEL : electron density [10^20 /m^3]
!     ANIL : ion density [10^20 /m^3]
!     TEL  : electron temperature [kev]
!     ZL   : ion charge number

      USE TRCOMM, ONLY : AEE, AME, EPS0, PI, PZ, RKEV
      IMPLICIT NONE
      REAL(8) :: ANEL, ANIL, TEL, ZL
      REAL(8) :: COEF, COULOG

      COEF = 6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)/(AEE**4*1.D20)
      IF(ZL-PZ(2).LE.1.D-7) THEN
         FTAUE = COEF*(TEL*RKEV)**1.5D0/(ANIL*ZL**2*COULOG(1,2,ANEL,TEL))
      ELSE
!     If the plasma contains impurities, we need to consider the
!     effective charge number instead of ion charge number.
!     From the definition of Zeff=sum(n_iZ_i^2)/n_e,
!     n_iZ_i^2 is replaced by n_eZ_eff at the denominator of tau_e.
         FTAUE = COEF*(TEL*RKEV)**1.5D0/(ANEL*ZL*COULOG(1,2,ANEL,TEL))
      ENDIF

      RETURN
      END FUNCTION FTAUE

!     between ions and ions

      REAL(8) FUNCTION FTAUI(ANEL,ANIL,TIL,ZL,PAL)

!     ANEL : electron density [10^20 /m^3]
!     ANIL : ion density [10^20 /m^3]
!     TIL  : ion temperature [kev]
!     ZL   : ion charge number
!     PAL  : ion atomic number

      USE TRCOMM, ONLY : AEE, AMM, EPS0, PI, RKEV
      IMPLICIT NONE
      REAL(8):: ANEL, ANIL, PAL, TIL, ZL
      REAL(8):: COEF, COULOG

      COEF = 12.D0*PI*SQRT(PI)*EPS0**2*SQRT(PAL*AMM)/(AEE**4*1.D20)
      FTAUI = COEF*(TIL*RKEV)**1.5D0/(ANIL*ZL**4*COULOG(2,2,ANEL,TIL))

      RETURN
      END FUNCTION FTAUI
