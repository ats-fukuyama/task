C     $Id$
C  
C     ***********************************************************
C
C           CALCULATE TRANSPORT COEFFICIENTS AND SOURCE
C
C     ***********************************************************
C
      SUBROUTINE TRCALC(IERR)
C
      INCLUDE 'trcomm.inc'
      DIMENSION PWRPL(NRM,NSM),AJWRPL(NRM)
      DIMENSION PWMPL(NRM,NSM),AJWMPL(NRM)
      CHARACTER KID*80
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      IERR=0
C
      DO NR=1,NRMAX
         SIE(NR)=0.D0
         SNF(NR)=0.D0
         SNB(NR)=0.D0
         POH(NR)=0.D0
         PIE(NR)=0.D0
         PCX(NR)=0.D0
         PRL(NR)=0.D0
         PNB(NR)=0.D0
         PNF(NR)=0.D0
         PBIN(NR)=0.D0
         PFIN(NR)=0.D0
C         AJNB(NR)=0.D0
         AJRF(NR)=0.D0
         AJBS(NR)=0.D0
      DO NS=1,NSM
         SPE(NR,NS)=0.D0
         IF(MDLUF.NE.0) THEN
            IF(NS.GE.3) PRF(NR,NS)=0.D0
         ELSE
            PRF(NR,NS)=0.D0
         ENDIF
         PBCL(NR,NS)=0.D0
         PFCL(NR,NS)=0.D0
      ENDDO
      ENDDO
C
      IF(MODELG.NE.3) THEN
         IF(MDLEQB.NE.0) THEN
            DO NR=1,NRMAX
               QP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)
     &               /(4.D0*PI**2*RDP(NR))
            ENDDO
         ENDIF
      ENDIF
C     calculate q_axis
      Q0=FCTR(RG(1),RG(2),QP(1),QP(2))
C
C     *** RADIAL ELECTRIC FIELD ***
C
      CALL TRERAD
C
      IF(T.LT.PELTIM+0.5D0*DT.AND.
     &   T.GE.PELTIM-0.5D0*DT) THEN
         CALL TRPELT
      ENDIF
      CALL TRZEFF
C
      IF(MDNCLS.NE.0) THEN
         CALL TR_NCLASS(IERR)
         IF(IERR.NE.0) RETURN
      ENDIF
C
      CALL TRCOEF
      CALL TRLOSS
      IF(MDLUF.NE.1.AND.MDLUF.NE.3) CALL TRPWRF
      CALL TRPWNB
C
      CALL PLDATA_GETNW(NWRMAX,NWMMAX)
      DO NWR=1,NWRMAX
         CALL PLDATA_GETWR(NWR,KID,PWRPL,AJWRPL)
         DO NS=1,NSMAX
            DO NR=1,NRMAX
               PRF(NR,NS)=PRF(NR,NS)+PWRPL(NR,NS)
            ENDDO
         ENDDO
         DO NR=1,NRMAX
            AJRF(NR)=AJRF(NR)+AJWRPL(NR)
         ENDDO
      ENDDO
C
      DO NWM=1,NWMMAX
         CALL PLDATA_GETWM(NWM,KID,PWMPL,AJWMPL)
         DO NS=1,NSMAX
            DO NR=1,NRMAX
               PRF(NR,NS)=PRF(NR,NS)+PWMPL(NR,NS)
            ENDDO
         ENDDO
         DO NR=1,NRMAX
            AJRF(NR)=AJRF(NR)+AJWMPL(NR)
         ENDDO
      ENDDO
C
      IF(MDNCLS.NE.0) THEN
         CALL TRAJBS_NCLASS
      ELSE
         IF(MDLJBS.EQ.1) THEN
            CALL TRAJBS
         ELSEIF(MDLJBS.EQ.2) THEN
            CALL TRAJBS
         ELSEIF(MDLJBS.EQ.3) THEN
            CALL TRAJBS
         ELSEIF(MDLJBS.EQ.4) THEN
            CALL TRAJBSNEW
         ELSEIF(MDLJBS.EQ.5) THEN
            CALL TRAJBSSAUTER
         ELSE
            CALL TRAJBS
         ENDIF
      ENDIF
C
      CALL TRALPH
      CALL TRAJOH
C
      DO NR=1,NRMAX
         IF(MDLEQ0.EQ.0) THEN
         SSIN(NR,1)=SIE(NR)                            +SNB(NR)
     &             +SEX(NR,1)
         SSIN(NR,2)=PN(2)*SIE(NR)/(PN(2)+PN(3))-SNF(NR)+SNB(NR)
     &             +SEX(NR,2)
         SSIN(NR,3)=PN(3)*SIE(NR)/(PN(2)+PN(3))-SNF(NR)        
     &             +SEX(NR,3)
         SSIN(NR,4)=                            SNF(NR)        
     &             +SEX(NR,4)
         SSIN(NR,7)=-SIE(NR)                                   -SCX(NR)
         SSIN(NR,8)=                                    SNB(NR)+SCX(NR)
         ELSEIF(MDLEQ0.EQ.1) THEN
         SSIN(NR,1)=                                    SNB(NR)
     &             +SEX(NR,1)
         SSIN(NR,2)=                           -SNF(NR)+SNB(NR)
     &             +SEX(NR,2)
         SSIN(NR,3)=                           -SNF(NR)
     &             +SEX(NR,3)
         SSIN(NR,4)=                            SNF(NR)
     &             +SEX(NR,4)
         SSIN(NR,7)=0.D0
         SSIN(NR,8)=                                    SNB(NR)
         ENDIF
         PIN(NR,1)=PBCL(NR,1)+PFCL(NR,1)+PRF(NR,1)
     &            +POH(NR)-PRL(NR)-PIE(NR)+PEX(NR,1)
         PIN(NR,2)=PBCL(NR,2)+PFCL(NR,2)+PRF(NR,2)
     &            -PN(2)*PCX(NR)/(PN(2)+PN(3))+PEX(NR,2)
         PIN(NR,3)=PBCL(NR,3)+PFCL(NR,3)+PRF(NR,3)
     &            -PN(3)*PCX(NR)/(PN(2)+PN(3))+PEX(NR,3)
         PIN(NR,4)=PBCL(NR,4)+PFCL(NR,4)+PRF(NR,4)+PEX(NR,4)
      ENDDO
C
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
C
      RETURN
      END
C
C     ***********************************************************
C
C           RADIAL ELECTRIC FIELD
C
C     **********************************************************
C
      SUBROUTINE TRERAD
C
      INCLUDE 'trcomm.inc'
C
      DO NR=1,NRMAX
         IF(SUMPBM.EQ.0.D0) THEN
            PADD(NR)=0.D0
         ELSE
            PADD(NR)=PBM(NR)*1.D-20/RKEV-RNF(NR,1)*RT(NR,2)
         ENDIF
      ENDDO
      DO NR=1,NRMAX
         DRL=RJCB(NR)/DR
         IF(NR.EQ.NRMAX) THEN
            DPD = DERIV3P(PNSS(2)*PTS(2),
     &                    RN(NR  ,2)*RT(NR  ,2)-PADD(NR  ),
     &                    RN(NR-1,2)*RT(NR-1,2)-PADD(NR-1),
     &                    RHOG(NR),RHOM(NR),RHOM(NR-1))
            TERM_DP = DPD*RKEV/(PZ(2)*AEE*PNSS(2))
         ELSE
            DPD =(  RN(NR+1,2)*RT(NR+1,2)-PADD(NR+1)
     &            -(RN(NR  ,2)*RT(NR  ,2)-PADD(NR  )))*DRL
            TERM_DP = DPD*RKEV/(PZ(2)*AEE*0.5D0*(RN(NR+1,2)+RN(NR,2)))
         ENDIF
         IF(MDLER.EQ.0) THEN
C     pressure gradient only
            ER(NR) = TERM_DP
         ELSEIF(MDLER.EQ.1) THEN
C     nabla p + toroidal rotation
            ER(NR) = TERM_DP+VTOR(NR)*BP(NR)
         ELSEIF(MDLER.EQ.2) THEN
C     nabla p + V_tor + poloidal rotation
            ER(NR) = TERM_DP+VTOR(NR)*BP(NR)-VPOL(NR)*BB
         ELSEIF(MDLER.EQ.3) THEN
C     Waltz definition
            EPS = EPSRHO(NR)
            F_UNTRAP = 1.D0-1.46D0*SQRT(EPS)+0.46D0*EPS**1.5D0
            ALPHA_NEO = 1.D0-0.8839D0*F_UNTRAP
     &                 /(0.3477D0+0.4058D0*F_UNTRAP)
            IF(NR.EQ.NRMAX) THEN
               TEL = PTS(1)
               TIL = PTS(2)
               RLNI = -DERIV3P(PNSS(2),RN(NR,2),RN(NR-1,2),
     &                         RHOG(NR),RHOM(NR),RHOM(NR-1))/PNSS(2)
               RLTI = -DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2),
     &                         RHOG(NR),RHOM(NR),RHOM(NR-1))/PTS(2)
            ELSE
               TEL = 0.5D0*(RT(NR,1)+RT(NR+1,1))
               TIL = 0.5D0*(RT(NR,2)+RT(NR+1,2))
               RLNI = -(LOG(RN(NR+1,2))-LOG(RN(NR,2)))*DRL
               RLTI = -(LOG(RT(NR+1,2))-LOG(RT(NR,2)))*DRL
            ENDIF
            CS = SQRT(TEL*RKEV/(PA(2)*AMM))
            RHO_S = CS*PA(2)*AMM/(PZ(2)*AEE*BB)
            ER(NR) =-BB*( (TIL/TEL)*RHO_S*CS*(RLNI+ALPHA_NEO*RLTI)
     &                   -EPS/QP(NR)*VTOR(NR))
         ENDIF
      ENDDO
C
      RETURN
      END
C
C     ***********************************************************
C
C           CALCULATE Z-C
C
C     **********************************************************
C
      REAL*8 FUNCTION TRZEC(TE)
C
      IMPLICIT NONE
      REAL*8 TE,TEL
      REAL*8 BC00,BC01,BC02,BC03,BC04,BC05,BC10,BC11,BC12,BC13,BC14,BC15
      REAL*8 BC20,BC21,BC22,BC23,BC24,BC25,BC30,BC31,BC32,BC33,BC34,BC35
      REAL*8 BC40,BC41,BC42,BC43,BC44,BC45
C
      DATA BC00,BC01,BC02,BC03,BC04,BC05/
     &     2.093736D+03, 5.153766D+03, 5.042105D+03,
     &     2.445345D+03, 5.879207D+02, 5.609128D+01/
C
      DATA BC10,BC11,BC12,BC13,BC14,BC15/
     &     7.286051D+00, 4.506650D+01, 1.595483D+02,
     &     2.112702D+02, 1.186473D+02, 2.400040D+01/
C
      DATA BC20,BC21,BC22,BC23,BC24,BC25/
     &     5.998782D+00, 6.808206D-03,-3.096242D-02,
     &    -4.794194D-02, 1.374732D-01, 5.157571D-01/
C
      DATA BC30,BC31,BC32,BC33,BC34,BC35/
     &     5.988153D+00, 7.335329D-02,-1.754858D-01,
     &     2.034126D-01,-1.145930D-01, 2.517150D-02/
C
      DATA BC40,BC41,BC42,BC43,BC44,BC45/
     &    -3.412166D+01, 1.250454D+02,-1.550822D+02,
     &     9.568297D+01,-2.937297D+01, 3.589667D+00/
C
         TEL = LOG10(TE)
         IF(TE.LE.3.D-3) THEN
            TRZEC=0.D0
         ELSEIF(TE.LE.2.D-2) THEN
            TRZEC= BC00+(BC01*TEL)+(BC02*TEL**2)+(BC03*TEL**3)
     &                 +(BC04*TEL**4)+(BC05*TEL**5)
         ELSEIF(TE.LE.0.2D0) THEN
            TRZEC= BC10+(BC11*TEL)+(BC12*TEL**2)+(BC13*TEL**3)
     &                 +(BC14*TEL**4)+(BC15*TEL**5)
         ELSEIF(TE.LE.2.D0) THEN
            TRZEC= BC20+(BC21*TEL)+(BC22*TEL**2)+(BC23*TEL**3)
     &                 +(BC24*TEL**4)+(BC25*TEL**5)
         ELSEIF(TE.LE.20.D0) THEN
            TRZEC= BC30+(BC31*TEL)+(BC32*TEL**2)+(BC33*TEL**3)
     &                 +(BC34*TEL**4)+(BC35*TEL**5)
         ELSEIF(TE.LE.100.D0) THEN
            TRZEC= BC40+(BC41*TEL)+(BC42*TEL**2)+(BC43*TEL**3)
     &                 +(BC44*TEL**4)+(BC45*TEL**5)
         ELSE
            TRZEC= 6.D0
         ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
C           CALCULATE Z-FE
C
C     ***********************************************************
C
      REAL*8 FUNCTION TRZEFE(TE)
C
      IMPLICIT NONE
      REAL*8 TE,TEL
      REAL*8 BF10,BF11,BF12,BF13,BF14,BF15,BF20,BF21,BF22,BF23,BF24,BF25
      REAL*8 BF30,BF31,BF32,BF33,BF34,BF35,BF40,BF41,BF42,BF43,BF44,BF45
C
      DATA BF10,BF11,BF12,BF13,BF14,BF15/
     &     8.778318D+00,-5.581412D+01,-1.225124D+02,
     &    -1.013985D+02,-3.914244D+01,-5.851445D+00/
C
      DATA BF20,BF21,BF22,BF23,BF24,BF25/
     &     1.959496D+01, 1.997852D+01, 9.593373D+00,
     &    -7.609962D+01,-1.007190D+02,-1.144046D+01/
C
      DATA BF30,BF31,BF32,BF33,BF34,BF35/
     &     1.478150D+01, 6.945311D+01,-1.991202D+02,
     &     2.653804D+02,-1.631423D+02, 3.770026D+01/
C
      DATA BF40,BF41,BF42,BF43,BF44,BF45/
     &     2.122081D+01, 6.891607D+00,-4.076853D+00,
     &     1.577171D+00,-5.139468D-01, 8.934176D-02/
C
         TEL = LOG10(TE)
C
         IF(TE.LE.3.D-3) THEN
            TRZEFE=0.D0
         ELSEIF(TE.LE.0.2D0) THEN
            TRZEFE = BF10+(BF11*TEL)+(BF12*TEL**2)+(BF13*TEL**3)
     &                   +(BF14*TEL**4)+(BF15*TEL**5)
         ELSEIF(TE.LE.2.D0) THEN
            TRZEFE = BF20+(BF21*TEL)+(BF22*TEL**2)+(BF23*TEL**3)
     &                   +(BF24*TEL**4)+(BF25*TEL**5)
         ELSEIF(TE.LE.20.D0) THEN
            TRZEFE = BF30+(BF31*TEL)+(BF32*TEL**2)+(BF33*TEL**3)
     &                   +(BF34*TEL**4)+(BF35*TEL**5)
         ELSEIF(TE.LE.100.D0) THEN
            TRZEFE = BF40+(BF41*TEL)+(BF42*TEL**2)+(BF43*TEL**3)
     &                   +(BF44*TEL**4)+(BF45*TEL**5)
         ELSE
            TRZEFE = 26.D0
         ENDIF
      RETURN
      END
C
C     ***********************************************************
C
C           CALCULATE Z-EFFF
C
C     ***********************************************************
C
      SUBROUTINE TRZEFF
C
      INCLUDE 'trcomm.inc'
C
      IF(MDLUF.EQ.0) THEN
         DO NR=1,NRMAX
            TE =RT(NR,1)
            ZEFF(NR) =(PZ(2)  *PZ(2)  *RN(NR,2)
     &                +PZ(3)  *PZ(3)  *RN(NR,3)
     &                +PZ(4)  *PZ(4)  *RN(NR,4)
     &                +TRZEC(TE)**2   *ANC (NR)
     &                +TRZEFE(TE)**2  *ANFE(NR))/RN(NR,1)
         ENDDO
      ELSE
         IF(MDLEQN.EQ.0) THEN ! fixed density
            DO NR=1,NRMAX
               ZEFF(NR) =(PZ(2)  *PZ(2)  *RN(NR,2)
     &                   +PZ(3)  *PZ(3)  *RN(NR,3)
     &                   +PZ(2)  *PZ(2)  *RNF(NR,1))/RN(NR,1)
            ENDDO
         ENDIF
      ENDIF
C
      DO NR=1,NRMAX
         TE=RT(NR,1)
         PZC(NR)=TRZEC(TE)
         PZFE(NR)=TRZEFE(TE)
      ENDDO
C     
      RETURN
      END
C
C     ***********************************************************
C
C           RADIATION POWER - C
C
C     ***********************************************************
C
      REAL*8 FUNCTION TRRPC(TE)
C
      IMPLICIT NONE
      REAL*8 TE,TEL,ARG
      REAL*8 AC00,AC01,AC02,AC03,AC04,AC05,AC10,AC11,AC12,AC13,AC14,AC15
      REAL*8 AC20,AC21,AC22,AC23,AC24,AC25,AC30,AC31,AC32,AC33,AC34,AC35
      REAL*8 AC40,AC41,AC42,AC43,AC44,AC45
C
      DATA AC00,AC01,AC02,AC03,AC04,AC05/
     &     1.965300D+03, 4.572039D+03, 4.159590D+03,
     &     1.871560D+03, 4.173889D+02, 3.699382D+01/
C
      DATA AC10,AC11,AC12,AC13,AC14,AC15/
     &     7.467599D+01, 4.549038D+02, 8.372937D+02,
     &     7.402515D+02, 3.147607D+02, 5.164578D+01/
C
      DATA AC20,AC21,AC22,AC23,AC24,AC25/
     &    -2.120151D+01,-3.668933D-01, 7.295099D-01,
     &    -1.944827D-01,-1.263576D-01,-1.491027D-01/
C
      DATA AC30,AC31,AC32,AC33,AC34,AC35/
     &    -2.121979D+01,-2.346986D-01, 4.093794D-01,
     &     7.874548D-02,-1.841379D-01, 5.590744D-02/
C
      DATA AC40,AC41,AC42,AC43,AC44,AC45/
     &    -2.476796D+01, 9.408181D+00,-9.657446D+00,
     &     4.999161D+00,-1.237382D+00, 1.160610D-01/
C
         TEL = LOG10(TE)
         IF(TE.LE.3.D-3) THEN
            TEL=LOG10(3.D-3)
            ARG=AC00+(AC01*TEL)+(AC02*TEL**2)+(AC03*TEL**3)
     &              +(AC04*TEL**4)+(AC05*TEL**5)
         ELSEIF(TE.LE.2.D-2) THEN
            ARG=AC00+(AC01*TEL)+(AC02*TEL**2)+(AC03*TEL**3)
     &              +(AC04*TEL**4)+(AC05*TEL**5)
         ELSEIF(TE.LE.0.2D0) THEN
            ARG=AC10+(AC11*TEL)+(AC12*TEL**2)+(AC13*TEL**3)
     &              +(AC14*TEL**4)+(AC15*TEL**5)
         ELSEIF(TE.LE.2.D0) THEN
            ARG=AC20+(AC21*TEL)+(AC22*TEL**2)+(AC23*TEL**3)
     &              +(AC24*TEL**4)+(AC25*TEL**5)
         ELSEIF(TE.LE.20.D0) THEN
            ARG=AC30+(AC31*TEL)+(AC32*TEL**2)+(AC33*TEL**3)
     &              +(AC34*TEL**4)+(AC35*TEL**5)
         ELSEIF(TE.LE.1.D2) THEN
            ARG=AC40+(AC41*TEL)+(AC42*TEL**2)+(AC43*TEL**3)
     &              +(AC44*TEL**4)+(AC45*TEL**5)
         ENDIF
         TRRPC = 10.D0**(ARG-13.D0)
      RETURN
      END
C
C     ***********************************************************
C
C           RADIATION POWER - FE
C
C     ***********************************************************
C
      REAL*8 FUNCTION TRRPFE(TE)
C
      IMPLICIT NONE
      REAL*8 TE,TEL,ARG
      REAL*8 AF10,AF11,AF12,AF13,AF14,AF15,AF20,AF21,AF22,AF23,AF24,AF25
      REAL*8 AF30,AF31,AF32,AF33,AF34,AF35,AF40,AF41,AF42,AF43,AF44,AF45
C
      DATA AF10,AF11,AF12,AF13,AF14,AF15/
     &    -2.752599D+01,-3.908228D+01,-6.469423D+01,
     &    -5.555048D+01,-2.405568D+01,-4.093160D+00/
C
      DATA AF20,AF21,AF22,AF23,AF24,AF25/
     &    -1.834973D+01,-1.252028D+00,-7.533115D+00,
     &    -3.289693D+00, 2.866739D+01, 2.830249D+01/
C
      DATA AF30,AF31,AF32,AF33,AF34,AF35/
     &    -1.671042D+01,-1.646143D+01, 3.766238D+01,
     &    -3.944080D+01, 1.918529D+01,-3.509238D+00/
C
      DATA AF40,AF41,AF42,AF43,AF44,AF45/
     &    -2.453957D+01, 1.795222D+01,-2.356360D+01,
     &     1.484503D+01,-4.542323D+00, 5.477462D-01/
C
         TEL = LOG10(TE)
         IF(TE.LE.3.D-3) THEN
            TEL=LOG10(3.D-3)
            ARG=AF10+(AF11*TEL)+(AF12*TEL**2)+(AF13*TEL**3)
     &              +(AF14*TEL**4)+(AF15*TEL**5)
         ELSEIF(TE.LE.0.2D0) THEN
            ARG=AF10+(AF11*TEL)+(AF12*TEL**2)+(AF13*TEL**3)
     &              +(AF14*TEL**4)+(AF15*TEL**5)
         ELSEIF(TE.LE.2.D0) THEN
            ARG=AF20+(AF21*TEL)+(AF22*TEL**2)+(AF23*TEL**3)
     &              +(AF24*TEL**4)+(AF25*TEL**5)
         ELSEIF(TE.LE.20.D0) THEN
            ARG=AF30+(AF31*TEL)+(AF32*TEL**2)+(AF33*TEL**3)
     &              +(AF34*TEL**4)+(AF35*TEL**5)
         ELSEIF(TE.LE.1.D2) THEN
            ARG=AF40+(AF41*TEL)+(AF42*TEL**2)+(AF43*TEL**3)
     &              +(AF44*TEL**4)+(AF45*TEL**5)
         ENDIF
         TRRPFE = 10.D0**(ARG-13.D0)
      RETURN
      END
C
C     ***********************************************************
C
C           POWER LOSS
C
C     ***********************************************************
C
      SUBROUTINE TRLOSS
C
      INCLUDE 'trcomm.inc'
      COMMON /TMSLC1/ TMU(NTUM),TMU1(NTUM)
      COMMON /TMSLC2/ NTAMAX
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
      COMMON /TMSLC4/ PNBI
C
      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         IF(NT.EQ.0) THEN
            TSL=DT*DBLE(1)
            DO NR=1,NRMAX
C               CALL TIMESPL(TSL,PRLL,TMU,PRLU(1,NR),NTXMAX,NTUM,IERR)
C               PRL(NR)=PRLL
               PRL(NR)=PRLU(1,NR)
            ENDDO
         ELSE
            TSL=DT*DBLE(NT)
            IF(KUFDEV.EQ.'X') THEN
               DO NR=1,NRMAX
                  CALL TIMESPL(TSL,PRLL,TMU,PRLU(1,NR),NTXMAX,NTUM,
     &                          IERR)
                  PRL(NR)=PRLL
               ENDDO
            ELSE
               DO NR=1,NRMAX
                  IF(PNBI.LT.12.D6) THEN
                     PRL(NR)=PRLU(1,NR)
                  ELSE
                     PRL(NR)=PRLU(2,NR)
                     IF(NT.EQ.NTAMAX) PRL(NR)=PRLU(3,NR)
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ELSEIF(MDLUF.EQ.2) THEN
         IF(NT.EQ.0) THEN
            TSL=DT*DBLE(1)
            DO NR=1,NRMAX
               PRL(NR)=PRLU(1,NR)
            ENDDO
         ELSE
            TSL=DT*DBLE(NT)
            DO NR=1,NRMAX
               PRL(NR)=PRLU(1,NR)
            ENDDO
         ENDIF
      ELSE
         DO NR=1,NRMAX
            ANE =RN(NR,1)
            ANDX=RN(NR,2)
            ANT =RN(NR,3)
            ANHE=RN(NR,4)
            TE  =RT(NR,1)
C     Radiation loss caused by impurities
            PLFE  = ANE*ANFE(NR)*TRRPFE(TE)*1.D40
            PLC   = ANE*ANC (NR)*TRRPC (TE)*1.D40
C     Bremsstrahlung
            PLD   = ANE*ANDX*5.35D-37*1.D0**2*SQRT(ABS(TE))*1.D40
            PLTT  = ANE*ANT *5.35D-37*1.D0**2*SQRT(ABS(TE))*1.D40
            PLHE  = ANE*ANHE*5.35D-37*2.D0**2*SQRT(ABS(TE))*1.D40
            PRL(NR)= PLFE+PLC+PLD+PLTT+PLHE
         ENDDO
      ENDIF
C
C     ****** IONIZATION LOSS ******
C
      DO NR=1,NRMAX
         ANE=RN(NR,1)
         TE =RT(NR,1)
         EION  = 13.64D0
C         SIONS  = 1.76D-13/SQRT(TE*1.D3)*FKN(EION/(TE*1.D3))
         TN    = MAX(TE*1.D3/EION,1.D-2)
         SION  = 1.D-11*SQRT(TN)*EXP(-1.D0/TN)
     &          /(EION**1.5D0*(6.D0+TN))
         PIE(NR) = ANE*ANNU(NR)*SION*1.D40*EION*AEE
         SIE(NR) = ANE*ANNU(NR)*SION*1.D20
         TSIE(NR)= ANE         *SION*1.D20
      ENDDO
C
C     ****** CHARGE EXCHANGE LOSS ******
C
      DO NR=1,NRMAX
         ANE =RN(NR,1)
         ANDX=RN(NR,2)
         TE  =RT(NR,1)
         TD  =RT(NR,2)
         TNU = 0.D0
C
         TN    = LOG10(MAX(TD*1.D3,50.D0))
         SCH= 1.57D-16*SQRT(ABS(TD)*1.D3)*(TN*TN-14.63D0*TN+53.65D0)
         SCX(NR) =ANDX*ANNU(NR)*SCH*1.D20
         TSCX(NR)=ANDX         *SCH*1.D20
C
         PCX(NR)=1.5D0*ANDX*ANNU(NR)*SCH*(TD-TNU)*RKEV*1.D40
      ENDDO
C
      RETURN
      END
C
      REAL*8 FUNCTION FKN(X)
C
      IMPLICIT NONE
      REAL*8 X,GAMMA
C
      GAMMA=0.577215664901532D0
      FKN=-(LOG10(X)+GAMMA-X+X**2/4.D0-X**3/1.8D1)
C
      RETURN
      END
C
C     ***********************************************
C
C         BOOTSTRAP CURRENT (NCLASS)
C
C     ***********************************************
C
      SUBROUTINE TRAJBS_NCLASS
C
      INCLUDE 'trcomm.inc'
      DIMENSION AJBSL(NRM)
C
      IF(PBSCD.LE.0.D0) RETURN
C
      NSW=1
      IF(NSW.EQ.0) THEN
         DO NR=1,NRMAX
            AJBSL(NR)=PBSCD*AJBSNC(NR)
         ENDDO
      ELSE
         DO NR=1,NRMAX-1
            SUM=0.D0
            DO NS=1,NSMAX
               RTNW=0.5D0*(RT(NR+1,NS)+RT(NR,NS))
               RPNW=0.5D0*(RN(NR+1,NS)*RT(NR+1,NS)
     &                    +RN(NR  ,NS)*RT(NR  ,NS))
               DRTNW=(RT(NR+1,NS)-RT(NR,NS))/DR
               DRPNW=(RN(NR+1,NS)*RT(NR+1,NS)-RN(NR,NS)*RT(NR,NS))/DR
               SUM=SUM+CJBST(NR,NS)*DRTNW/RTNW
     &                +CJBSP(NR,NS)*DRPNW/RPNW
            ENDDO
            AJBSL(NR)=-PBSCD*SUM/BB
         ENDDO
C
         NR=NRMAX
         SUM=0.D0
         DO NS=1,NSMAX
            RTNW=PTS(NS)
            RPNW=PNSS(NS)*PTS(NS)
            DRTNW=DERIV3P(PTS(NS),RT(NR,NS),RT(NR-1,NS),
     &                    RG(NR),RM(NR),RM(NR-1))
            DRPNW=DERIV3P(PNSS(NS)*PTS(NS),
     &                    RN(NR  ,NS)*RT(NR  ,NS),
     &                    RN(NR-1,NS)*RT(NR-1,NS),
     &                    RG(NR),RM(NR),RM(NR-1))
            SUM=SUM+CJBST(NR,NS)*DRTNW/RTNW
     &             +CJBSP(NR,NS)*DRPNW/RPNW
CCC            if(ns.eq.3) write(6,*) NR,PNSS(NS),PTS(NS)
         ENDDO
         AJBSL(NR)=-PBSCD*SUM/BB
      ENDIF
C     
      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO
C
      RETURN
      END
C
C     ***********************************************
C
C         BOOTSTRAP CURRENT (O. SAUTER)
C
C     ***********************************************
C
      SUBROUTINE TRAJBSSAUTER
C
      INCLUDE 'trcomm.inc'
      DIMENSION ANI(NRM),AJBSL(NRM)
C
      IF(PBSCD.LE.0.D0) RETURN
C
      DO NR=1,NRMAX-1
C
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         QL=ABS(QP(NR))
         ZEFFL=0.5D0*(ZEFF(NR)+ZEFF(NR+1))
         DRL=1.D0/DR
C
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
C
C     ****** ION PARAMETER ******
C
C     *** ANI  is the the ion density (ni) ***
C     *** TI   is the ion temperature (Ti) ***
C     *** DTI  is the derivative of ion temperature (dTi/dr) ***
C     *** PPI  is the ion pressure (Pi) ***
C     *** DPI  is the derivative of ion pressure (dPi/dr) ***
C
         ANI(NR)=0.D0
         DO NS=2,NSM
            ANI(NR)=ANI(NR)+0.5D0*PZ(NS)*(RN(NR+1,NS)+RN(NR,NS))
         ENDDO
         PPI=0.5D0*(RPIP+RPIM)
         DPI=(RPIP-RPIM)*DRL
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         DTI=(RNTP/RNP-RNTM/RNM)*DRL
C
         rLnLamii=30.D0-LOG(PZ(2)**3*SQRT(ANI(NR)*1.D20)
     &           /(ABS(TI*1.D3)**1.5D0))
         RNUI=4.90D-18*QL*RR*ANI(NR)*1.D20*PZ(2)**4*rLnLamii
     &       /(ABS(TI*1.D3)**2*EPSS)
C     
C     ****** ELECTORON PARAMETER ******
C
C     *** ANE  is the the electron density (ne) ***
C     *** TE   is the electron temperature (Te) ***
C     *** PE   is the electron pressure (Pe) ***
C     *** DTE  is the derivative of electron temperature (dTe/dr) ***
C     *** DPE  is the derivative of electron pressure (dPe/dr) ***
C
         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
         TE= 0.5D0*(RT(NR+1,1)+RT(NR,1))
         PE= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
         DTE=(RT(NR+1,1)-RT(NR,1))*DRL
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL
C
         rLnLame=31.3D0-LOG(SQRT(ANE*1.D20)/ABS(TE*1.D3))
         RNUE=6.921D-18*QL*RR*ANE*1.D20*ZEFFL*rLnLame
     &       /(ABS(TE*1.D3)**2*EPSS)
C
         RPE=PE/(PE+PPI)
         FT=FTPF(MDLTPF,EPS)
C
C         F33TEFF=FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE)
C     &          +0.45D0*(1.D0-FT)*RNUE/ZEFFL**1.5)
         F31TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)
     &          +0.5D0*(1.D0-FT)*RNUE/ZEFFL)
         F32EETEFF=FT/(1.D0+0.26D0*(1.D0-FT)*SQRT(RNUE)
     &            +0.18D0*(1.D0-0.37D0*FT)*RNUE/SQRT(ZEFFL))
         F32EITEFF=FT/(1.D0+(1.D0+0.6D0*FT)*SQRT(RNUE)
     &            +0.85D0*(1.D0-0.37D0*FT)*RNUE*(1.D0+ZEFFL))
         F34TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)
     &          +0.5D0*(1.D0-0.5D0*FT)*RNUE/ZEFFL)
C
         SALFA0=-1.17D0*(1.D0-FT)/(1.D0-0.22D0*FT-0.19D0*FT**2)
         SALFA=((SALFA0+0.25D0*(1.D0-FT**2)*SQRT(RNUI))
     &        /(1.D0+0.5D0*SQRT(RNUI))+0.315D0*RNUI**2*FT**6)
     &        /(1.D0+0.15D0*RNUI**2*FT**6)
C
C         RNZ=0.58D0+0.74D0/(0.76D0+ZEFFL)
C         SGMSPTZ=1.9012D4*(TE*1.D3)**1.5/(ZEFFL*RNZ*rLnLame)
C         SGMNEO=SGMSPTZ*F33(F33TEFF,ZEFFL)
C
         RL31=F31(F31TEFF,ZEFFL)
         RL32=F32EE(F32EETEFF,ZEFFL)+F32EI(F32EITEFF,ZEFFL)
         RL34=F31(F34TEFF,ZEFFL)
C
         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV
     &           *( RL31*(DPE/PE+DPI/PE)+RL32*DTE/TE
     &             +RL34*SALFA*(1.D0-RPE)/RPE*DTI/TI)/RDP(NR)/BB
      ENDDO
C
      NR=NRMAX
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         QL=ABS(QP(NR))
         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
         DRL=1.D0/DR
C
C     In the following, we assume that
C        1. pressures of beam and fusion at rho=1 are negligible,
C        2. calibration of Pbeam (from UFILE) is not necessary at rho=1.
C
         RNTP=0.D0
         RNP =0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSM
            RNTP=RNTP+PNSS(NS)*PTS(NS)
            RNP =RNP +PNSS(NS)
            RNTM=RNTM+RN(NR-1,NS)*RT(NR-1,NS)
     &               +RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR-1,NS)+RN(NR  ,NS)
         ENDDO
         RNTM=RNTM+RW(NR-1,1)+RW(NR-1,2)
     &            +RW(NR  ,1)+RW(NR  ,2)
         RPIP=RNTP
         RPIM=RNTM+PADD(NR-1)+PADD(NR  )
         RNTM=0.5D0*RNTM
         RNM =0.5D0*RNM
         RPIM=0.5D0*RPIM
C
C     ****** ION PARAMETER ******
C
C     *** ANI  is the the ion density (ni) ***
C     *** TI   is the ion temperature (Ti) ***
C     *** DTI  is the derivative of ion temperature (dTi/dr) ***
C     *** PPI  is the ion pressure (Pi) ***
C     *** DPI  is the derivative of ion pressure (dPi/dr) ***
C
         ANI(NR)=0.D0
         DO NS=2,NSM
            ANI(NR)=ANI(NR)+PZ(NS)*PNSS(NS)
         ENDDO
         PPI=RPIP
         DPI=(RPIP-RPIM)*DRL
         TI =RNTP/RNP
         DTI=(RNTP/RNP-RNTM/RNM)*DRL
C
         rLnLamii=30.D0-LOG(PZ(2)**3*SQRT(ANI(NR)*1.D20)
     &        /(ABS(TI*1.D3)**1.5D0))
         RNUI=4.90D-18*QL*RR*ANI(NR)*1.D20*PZ(2)**4*rLnLamii
     &        /(ABS(TI*1.D3)**2*EPSS)
C     
C     ****** ELECTORON PARAMETER ******
C
C     *** ANE  is the the electron density (ne) ***
C     *** TE   is the electron temperature (Te) ***
C     *** PE   is the electron pressure (Pe) ***
C     *** DTE  is the derivative of electron temperature (dTe/dr) ***
C     *** DPE  is the derivative of electron pressure (dPe/dr) ***
C
         ANE=PNSS(1)
         TE =PTS(1)
         PE =PNSS(1)*PTS(1)
         DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),
     &               RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPE=DERIV3P(PNSS(1)*PTS(1),
     &               RN(NR,1)*RT(NR,1),
     &               RN(NR-1,1)*RT(NR-1,1),
     &               RHOG(NR),RHOM(NR),RHOM(NR-1))
C     
         rLnLame=31.3D0-LOG(SQRT(ANE*1.D20)/ABS(TE*1.D3))
         RNUE=6.921D-18*QL*RR*ANE*1.D20*ZEFFL*rLnLame
     &        /(ABS(TE*1.D3)**2*EPSS)
C     
         RPE=PE/(PE+PPI)
         FT=FTPF(MDLTPF,EPS)
C
C     F33TEFF=FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE)
C     &          +0.45D0*(1.D0-FT)*RNUE/ZEFFL**1.5)
         F31TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)
     &        +0.5D0*(1.D0-FT)*RNUE/ZEFFL)
         F32EETEFF=FT/(1.D0+0.26D0*(1.D0-FT)*SQRT(RNUE)
     &        +0.18D0*(1.D0-0.37D0*FT)*RNUE/SQRT(ZEFFL))
         F32EITEFF=FT/(1.D0+(1.D0+0.6D0*FT)*SQRT(RNUE)
     &        +0.85D0*(1.D0-0.37D0*FT)*RNUE*(1.D0+ZEFFL))
         F34TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)
     &        +0.5D0*(1.D0-0.5D0*FT)*RNUE/ZEFFL)
C
         SALFA0=-1.17D0*(1.D0-FT)/(1.D0-0.22D0*FT-0.19D0*FT**2)
         SALFA=((SALFA0+0.25D0*(1.D0-FT**2)*SQRT(RNUI))
     &        /(1.D0+0.5D0*SQRT(RNUI))+0.315D0*RNUI**2*FT**6)
     &        /(1.D0+0.15D0*RNUI**2*FT**6)
C
C     RNZ=0.58D0+0.74D0/(0.76D0+ZEFFL)
C     SGMSPTZ=1.9012D4*(TE*1.D3)**1.5/(ZEFFL*RNZ*rLnLame)
C     SGMNEO=SGMSPTZ*F33(F33TEFF,ZEFFL)
C     
         RL31=F31(F31TEFF,ZEFFL)
         RL32=F32EE(F32EETEFF,ZEFFL)+F32EI(F32EITEFF,ZEFFL)
         RL34=F31(F34TEFF,ZEFFL)
C
         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV
     &           *( RL31*(DPE/PE+DPI/PE)+RL32*DTE/TE
     &             +RL34*SALFA*(1.D0-RPE)/RPE*DTI/TI)/RDP(NR)/BB
C
      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO
C
      RETURN
      END
C
C     *********************
C     *  Fitting Function *
C     *********************
C
      REAL*8 FUNCTION F33(X,Z)
C
      IMPLICIT NONE
      REAL*8 X,Z
C
      F33 = 1.D0-(1.D0+0.36D0/Z)*X+0.59D0/Z*X**2-0.23D0/Z*X**3
C
      RETURN
      END
C
      REAL*8 FUNCTION F31(X,Z)
C
      IMPLICIT NONE
      REAL*8 X,Z
C
      F31 = (1.D0+1.4D0/(Z+1.D0))*X-1.9D0/(Z+1.D0)*X**2
     &     +0.3D0/(Z+1.D0)*X**3+0.2D0/(Z+1.D0)*X**4
C
      RETURN
      END
C
      REAL*8 FUNCTION F32EE(X,Z)
C
      IMPLICIT NONE
      REAL*8 X,Z
C
      F32EE = (0.05D0+0.62D0*Z)/(Z*(1.D0+0.44D0*Z))*(X-X**4)
     &       +1.D0/(1.D0+0.22D0*Z)*(X**2-X**4-1.2D0*(X**3-X**4))
     &       +1.2D0/(1.D0+0.5D0*Z)*X**4
C
      RETURN
      END
C
      REAL*8 FUNCTION F32EI(X,Z)
C
      IMPLICIT NONE
      REAL*8 X,Z
C
      F32EI =-(0.56D0+1.93D0*Z)/(Z*(1.D0+0.44D0*Z))*(X-X**4)
     &       +4.95D0/(1.D0+2.48D0*Z)*(X**2-X**4-0.55D0*(X**3-X**4))
     &       -1.2D0/(1.D0+0.5D0*Z)*X**4
C
      RETURN
      END
C
C     ************************************************
C
C         BOOTSTRAP CURRENT (Hirshman (Wilson))
C
C     ************************************************
C
      SUBROUTINE TRAJBSNEW
C
      INCLUDE 'trcomm.inc'
      DIMENSION ANI(NRM),AJBSL(NRM)
C
      IF(PBSCD.LE.0.D0) RETURN
C
      DO NR=1,NRMAX-1
C
         EPS=EPSRHO(NR)
C         EPSS=SQRT(EPS)**3
C         QL=ABS(QP(NR))
C         ZEFFL=0.5D0*(ZEFF(NR)+ZEFF(NR+1))
         DRL=1.D0/DR
C
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
C
C     ****** ION PARAMETER ******
C
C     ***** ANI  is the the ion density (ni) *****
C     ***** TI   is the ion temperature (Ti) *****
C     ***** DTI  is the derivative of ion temperature (dTi/dr) *****
C     ***** PPI  is the ion pressure (Pi) *****
C     ***** DPI  is the derivative of ion pressure (dPi/dr) *****
C     ***** VTI  is the ion velocity (VTi) *****
C
         ANI(NR)=0.D0
         DO NS=2,NSM
            ANI(NR)=ANI(NR)+0.5D0*PZ(NS)*(RN(NR+1,NS)+RN(NR,NS))
         ENDDO
         PPI=0.5D0*(RPIP+RPIM)
         DPI=(RPIP-RPIM)*DRL
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         DTI=(RNTP/RNP-RNTM/RNM)*DRL
C         VTI=SQRT(ABS(TI)*RKEV/AMM)
C
C         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
C         rLnLam=17.3D0-DLOG(ANE)*0.5D0+DLOG(ABS(TI))*1.5D0
C         TAUI=12.D0*PI*SQRT(PI)*EPS0**2*SQRT(AMM)
C     &             *(ABS(TI)*RKEV)**1.5D0/(ANI(NR)*1.D20
C     &             *ZEFFL**4*AEE**4*rLnLam)
C
C         RNUI=QL*RR/(TAUI*VTI*EPSS)
C
C     ****** ELECTORON PARAMETER ******
C
C     ***** ANE  is the the electron density (ne) *****
C     ***** TE   is the electron temperature (Te) *****
C     ***** DTE  is the derivative of electron temperature (dTe/dr) ****
C     ***** PE   is the electron pressure (Pe) *****
C     ***** DPE  is the derivative of electron pressure (dPe/dr) *****
C     ***** VTE  is the electron velocity (VTe) *****
C
         TE= 0.5D0*(RT(NR+1,1)+RT(NR,1))
         PE= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
         DTE=(RT(NR+1,1)-RT(NR,1))*DRL
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*DRL
C         VTE=SQRT(ABS(TE)*RKEV/AME)
C
C         rLnLam=15.2D0-DLOG(ANE)*0.5D0+DLOG(ABS(TE))
C         TAUE=6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)
C     &             *(ABS(TE)*RKEV)**1.5D0/(ANI*1.D20
C     &             *ZEFFL**2*AEE**4*rLnLam)
C
C         RNUE=QL*RR/(TAUE*VTE*EPSS)
C
C         FT=FTPF(MDLTPF,EPS)
         FT=(1.46D0*SQRT(EPS)+2.4D0*EPS)/(1.D0-EPS)**1.5D0
C
C         DDX=2.4D0+5.4D0*FT+2.6D0*FT**2
C         DDD=-1.17D0/(1.D0+0.46D0*FT)
C         C1=(4.D0+2.6D0*FT)
C     &      /((1.D0+1.02D0*DSQRT(RNUE)+1.07D0*RNUE)
C     &       *(1.D0+1.07D0*EPSS*RNUE))
C         C2=C1*TI/TE
C         C3=(7.D0+6.5D0*FT)
C     &      /((1.D0+0.57D0*DSQRT(RNUE)+0.61D0*RNUE)
C     &       *(1.D0+0.61D0*EPSS*RNUE))
C     &      -2.5D0*C1
C         C4=((DDD+0.35D0*DSQRT(RNUI))
C     &      /(1.D0+0.7D0*DSQRT(RNUI))
C     &      +2.1D0*EPS**3*RNUI**2)*C2
C     &      /((1.D0-EPS**3*RNUI**2)
C     &       *(1.D0+EPS**3*RNUE**2))
C
C         AJBSL(NR)=-PBSCD*FT*PE*1.D20*RKEV/(DDX*BP(NR))
C     &        *(C1*(DPE/PE)
C     &         +C2*(DPI/PPI)
C     &         +C3*(DTE/TE)
C     &         +C4*(DTI/TI))
C
C     *** S. P. Hirshman, Phys Fluids 31, 1988 3150 ***
C     *** cited (H.R. WILSON, Nucl.Fusion 32,no.2,1992 259-263) ***
C
         DDX=1.414D0*PZ(2)+PZ(2)**2+FT*(0.754D0+2.657D0*PZ(2)
     &        +2.D0*PZ(2)**2)+FT**2*(0.348D0+1.243D0*PZ(2)+PZ(2)**2)
         RL31= FT*(0.754D0+2.210D0*PZ(2)+PZ(2)**2
     &        +FT*(0.348D0+1.243D0*PZ(2)+PZ(2)**2))/DDX
         RL32=-FT*(0.884D0+2.074D0*PZ(2))/DDX
         DDD=-1.172D0/(1.D0+0.462D0*FT)
C
         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV
     &        *(RL31*((DPE/PE)+(TI/(PZ(2)*TE))
     &        *((DPI/PPI)+DDD*(DTI/TI)))+RL32*(DTE/TE))/RDP(NR)/BB
      ENDDO
C
      NR=NRMAX
         EPS=EPSRHO(NR)
C         EPSS=SQRT(EPS)**3
C         QL=ABS(QP(NR))
C         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
         DRL=1.D0/DR
C
C     In the following, we assume that
C        1. pressures of beam and fusion at rho=1 are negligible,
C        2. calibration of Pbeam (from UFILE) is not necessary at rho=1.
C
         RNTP=0.D0
         RNP= 0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSM
            RNTP=RNTP+PNSS(NS)*PTS(NS)
            RNP =RNP +PNSS(NS)
            RNTM=RNTM+RN(NR-1,NS)*RT(NR-1,NS)
     &               +RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR-1,NS)+RN(NR  ,NS)
         ENDDO
         RNTM=RNTP+RW(NR-1,1)+RW(NR-1,2)
     &            +RW(NR  ,1)+RW(NR  ,2)
         RPIP=RNTP
         RPIM=RNTM+PADD(NR-1)+PADD(NR  )
         RNTM=0.5D0*RNTM
         RNM =0.5D0*RNM
         RPIM=0.5D0*RPIM
C
C     ****** ION PARAMETER ******
C
C     ***** ANI  is the the ion density (ni) *****
C     ***** TI   is the ion temperature (Ti) *****
C     ***** DTI  is the derivative of ion temperature (dTi/dr) *****
C     ***** PPI  is the ion pressure (Pi) *****
C     ***** DPI  is the derivative of ion pressure (dPi/dr) *****
C     ***** VTI  is the ion velocity (VTi) *****
C
         ANI(NR)=0.D0
         DO NS=2,NSM
            ANI(NR)=ANI(NR)+PZ(NS)*PNSS(NS)
         ENDDO
         PPI=0.5D0*(RPIP+RPIM)
         DPI=(RPIP-RPIM)*DRL
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         DTI=(RNTP/RNP-RNTM/RNM)*DRL
C         VTI=SQRT(ABS(TI)*RKEV/AMM)
C
C         ANE=PNSS(1)
C         rLnLam=17.3D0-DLOG(ANE)*0.5D0+DLOG(ABS(TI))*1.5D0
C         TAUI=12.D0*PI*SQRT(PI)*EPS0**2*SQRT(AMM)
C     &             *(ABS(TI)*RKEV)**1.5D0/(ANI(NR)*1.D20
C     &             *ZEFFL**4*AEE**4*rLnLam)
C
C         RNUI=QL*RR/(TAUI*VTI*EPSS)
C
C     ****** ELECTORON PARAMETER ******
C
C     ***** ANE  is the the electron density (ne) *****
C     ***** TE   is the electron temperature (Te) *****
C     ***** DTE  is the derivative of electron temperature (dTe/dr) ****
C     ***** PE   is the electron pressure (Pe) *****
C     ***** DPE  is the derivative of electron pressure (dPe/dr) *****
C     ***** VTE  is the electron velocity (VTe) *****
C
         TE =PTS(1)
         PE =PNSS(1)*PTS(1)
         DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),
     &               RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPE=DERIV3P(PNSS(1)*PTS(1),
     &               RN(NR,1)*RT(NR,1),
     &               RN(NR-1,1)*RT(NR-1,1),
     &               RHOG(NR),RHOM(NR),RHOM(NR-1))
C         VTE=SQRT(ABS(TE)*RKEV/AME)
C
C         rLnLam=15.2D0-DLOG(ANE)*0.5D0+DLOG(ABS(TE))
C         TAUE=6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)
C     &             *(ABS(TE)*RKEV)**1.5D0/(ANI*1.D20
C     &             *ZEFFL**2*AEE**4*rLnLam)
C
C         RNUE=QL*RR/(TAUE*VTE*EPSS)
C
C         FT=FTPF(MDLTPF,EPS)
         FT=(1.46D0*SQRT(EPS)+2.4D0*EPS)/(1.D0-EPS)**1.5D0
C
C         DDX=2.4D0+5.4D0*FT+2.6D0*FT**2
C         DDD=-1.17D0/(1.D0+0.46D0*FT)
C         C1=(4.D0+2.6D0*FT)
C     &      /((1.D0+1.02D0*DSQRT(RNUE)+1.07D0*RNUE)
C     &       *(1.D0+1.07D0*EPSS*RNUE))
C         C2=C1*TI/TE
C         C3=(7.D0+6.5D0*FT)
C     &      /((1.D0+0.57D0*DSQRT(RNUE)+0.61D0*RNUE)
C     &       *(1.D0+0.61D0*EPSS*RNUE))
C     &      -2.5D0*C1
C         C4=((DDD+0.35D0*DSQRT(RNUI))
C     &      /(1.D0+0.7D0*DSQRT(RNUI))
C     &      +2.1D0*EPS**3*RNUI**2)*C2
C     &      /((1.D0-EPS**3*RNUI**2)
C     &       *(1.D0+EPS**3*RNUE**2))
C
C         AJBSL(NR)=-PBSCD*FT*PE*1.D20*RKEV/(DDX*BP(NR))
C     &        *(C1*(DPE/PE)
C     &         +C2*(DPI/PPI)
C     &         +C3*(DTE/TE)
C     &         +C4*(DTI/TI))
C
C     *** S. P. Hirshman, Phys Fluids 31, 1988 3150 ***
C     *** cited (H.R. WILSON, Nucl.Fusion 32,no.2,1992 259-263) ***
C
         DDX=1.414D0*PZ(2)+PZ(2)**2+FT*(0.754D0+2.657D0*PZ(2)
     &        +2.D0*PZ(2)**2)+FT**2*(0.348D0+1.243D0*PZ(2)+PZ(2)**2)
         RL31=FT*( 0.754D0+2.21D0*PZ(2)+PZ(2)**2+FT*(0.348D0+1.243D0
     &            *PZ(2)+PZ(2)**2))/DDX
         RL32=-FT*(0.884D0+2.074D0*PZ(2))/DDX
         DDD=-1.172D0/(1.D0+0.462D0*FT)
C
         AJBSL(NR)=-PBSCD*TTRHOG(NR)*PE*1.D20*RKEV
     &        *(RL31*((DPE/PE)+(TI/(PZ(2)*TE))
     &        *((DPI/PPI)+DDD*(DTI/TI)))+RL32*(DTE/TE))/RDP(NR)/BB
C
      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO
C
      RETURN
      END
C
C     ***********************************************************
C
C           BOOTSTRAP CURRENT (Hinton & Hazeltine)
C
C     ***********************************************************
C
      SUBROUTINE TRAJBS
C
      INCLUDE 'trcomm.inc'
      DIMENSION AJBSL(NRM)
C     
C     ZEFF=1
C
C      DATA RK11,RA11,RB11,RC11/1.04D0,2.01D0,1.53D0,0.89D0/
C      DATA RK12,RA12,RB12,RC12/1.20D0,0.76D0,0.67D0,0.56D0/
C      DATA RK22,RA22,RB22,RC22/2.55D0,0.45D0,0.43D0,0.43D0/
      DATA RK13,RA13,RB13,RC13/2.30D0,1.02D0,1.07D0,1.07D0/
      DATA RK23,RA23,RB23,RC23/4.19D0,0.57D0,0.61D0,0.61D0/
C      DATA RK33,RA33,RB33,RC33/1.83D0,0.68D0,0.32D0,0.66D0/
C      DATA RK2 ,RA2 ,RB2 ,RC2 /0.66D0,1.03D0,0.31D0,0.74D0/
C
      IF(PBSCD.LE.0.D0) RETURN
C
      AMD=PA(2)*AMM
      AMT=PA(3)*AMM
      AMA=PA(4)*AMM
C
      DO NR=1,NRMAX-1
C
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
C
         TAUE = FTAUE(ANE,ANDX,TEL,ZEFFL)
         TAUD = FTAUI(ANE,ANDX,TDL,PZ(2),PA(2))
         TAUT = FTAUI(ANE,ANT ,TTL,PZ(3),PA(3))
         TAUA = FTAUI(ANE,ANA ,TAL,PZ(4),PA(4))
C
         VTE=SQRT(TEL*RKEV/AME)
         VTD=SQRT(TDL*RKEV/AMD)
         VTT=SQRT(TTL*RKEV/AMT)
         VTA=SQRT(TAL*RKEV/AMA)
C
         RNUE=ABS(QP(NR))*RR/(TAUE*VTE*EPSS)
         RNUD=ABS(QP(NR))*RR/(TAUD*VTD*EPSS)
         RNUT=ABS(QP(NR))*RR/(TAUT*VTT*EPSS)
         RNUA=ABS(QP(NR))*RR/(TAUA*VTA*EPSS)
C
C         RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)
C     &              +(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
C         RK12E=RK12*(1.D0/(1.D0+RA12*SQRT(RNUE)+RB12*RNUE)
C     &              +(EPSS*RC12)**2/RB12*RNUE/(1.D0+RC12*RNUE*EPSS))
C         RK22E=RK22*(1.D0/(1.D0+RA22*SQRT(RNUE)+RB22*RNUE)
C     &              +(EPSS*RC22)**2/RB22*RNUE/(1.D0+RC22*RNUE*EPSS))
         RK13E=RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)
     &             /(1.D0+RC13*RNUE*EPSS)
         RK23E=RK23/(1.D0+RA23*SQRT(RNUE)+RB23*RNUE)
     &             /(1.D0+RC23*RNUE*EPSS)
C         RK33E=RK33/(1.D0+RA33*SQRT(RNUE)+RB33*RNUE)
C     &             /(1.D0+RC33*RNUE*EPSS)
C
C         RK2D =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUD)+RB2 *RNUD)
C     &              +(EPSS*RC2 )**2/RB2 *RNUD/(1.D0+RC2 *RNUD*EPSS))
C         RK2T =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUT)+RB2 *RNUT)
C     &              +(EPSS*RC2 )**2/RB2 *RNUT/(1.D0+RC2 *RNUT*EPSS))
C         RK2A =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUA)+RB2 *RNUA)
C     &              +(EPSS*RC2 )**2/RB2 *RNUA/(1.D0+RC2 *RNUA*EPSS))
         RK3D=((1.17D0-0.35D0*SQRT(RNUD))/(1.D0+0.7D0*SQRT(RNUD))
     &         -2.1D0*(RNUD*EPSS)**2)/(1.D0+(RNUD*EPSS)**2)
         RK3T=((1.17D0-0.35D0*SQRT(RNUT))/(1.D0+0.7D0*SQRT(RNUT))
     &         -2.1D0*(RNUT*EPSS)**2)/(1.D0+(RNUT*EPSS)**2)
         RK3A=((1.17D0-0.35D0*SQRT(RNUA))/(1.D0+0.7D0*SQRT(RNUA))
     &         -2.1D0*(RNUA*EPSS)**2)/(1.D0+(RNUA*EPSS)**2)
C
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
C
         FACT=1.D0/(1.D0+(RNUE*EPS)**2)
         IF(NSMAX.EQ.2) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL
     &       +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
         ELSEIF(NSMAX.EQ.3) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL
     &       +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
     &       +(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)
         ELSEIF(NSMAX.EQ.4) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL
     &       +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
     &       +(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)
     &       +(PAL/PEL)*(DPA/PAL-RK3A*FACT*DTA/TAL)
         ENDIF
         H=BB/(BB+BP(NR))
         AJBSL(NR)=-PBSCD*H*SQRT(EPS)*ANE*1.D20*TEL*RKEV/BPL
     &            *(RK13E*A+RK23E*DTE/TEL)
C     
C        WRITE(6,'(I3,1P6E12.4)') NR,RNUE,RK13E,RK23E,A,BPL,AJBSL(NR)
C         IF(NR.EQ.NRMAX) THEN
C            WRITE(6,'(1P6E12.4)') DN,DTE,ANE,PNSS(1),TE,PTS(1)
C            WRITE(6,'(1P6E12.4)') RK13E*(RN(NR,1)*RT(NR,1)
C     &                   +RN(NR,2)*RT(NR,2)
C     &                   +RN(NR,3)*RT(NR,3)
C     &                   +RN(NR,4)*RT(NR,4)
C     &                   +RW(NR,1)
C     &                   +RW(NR,2)         )*DN/ANE,
C     &            +(RK23E-1.5D0*RK13E)*RN(NR,1)*DTE,
C     &            +RK13D*(RK23D-1.5D0)*RN(NR,2)*DTD,
C     &            BPL,AJBS(NR-1),AJBS(NR)
C         ENDIF
      ENDDO
C
      NR=NRMAX
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         ANE=PNSS(1)
         ANDX=PNSS(2)
C         ANT =PNSS(3)
C         ANA =PNSS(4)
         TEL=ABS(PTS(1))
         TDL=ABS(PTS(2))
         TTL=ABS(PTS(3))
         TAL=ABS(PTS(4))
         PEL=PNSS(1)*PTS(1)
         PDL=PNSS(2)*PTS(2)
         PTL=PNSS(3)*PTS(3)
         PAL=PNSS(4)*PTS(4)
         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
C
         TAUE = FTAUE(ANE,ANDX,TEL,ZEFFL)
         TAUD = FTAUI(ANE,ANDX,TDL,PZ(2),PA(2))
         TAUT = FTAUI(ANE,ANT ,TTL,PZ(3),PA(3))
         TAUA = FTAUI(ANE,ANA ,TAL,PZ(4),PA(4))
C
         VTE=SQRT(TEL*RKEV/AME)
         VTD=SQRT(TDL*RKEV/AMD)
         VTT=SQRT(TTL*RKEV/AMT)
         VTA=SQRT(TAL*RKEV/AMA)
C
         RNUE=ABS(QP(NR))*RR/(TAUE*VTE*EPSS)
         RNUD=ABS(QP(NR))*RR/(TAUD*VTD*EPSS)
         RNUT=ABS(QP(NR))*RR/(TAUT*VTT*EPSS)
         RNUA=ABS(QP(NR))*RR/(TAUA*VTA*EPSS)
C
C         RK11E=RK11*(1.D0/(1.D0+RA11*SQRT(RNUE)+RB11*RNUE)
C     &              +(EPSS*RC11)**2/RB11*RNUE/(1.D0+RC11*RNUE*EPSS))
C         RK12E=RK12*(1.D0/(1.D0+RA12*SQRT(RNUE)+RB12*RNUE)
C     &              +(EPSS*RC12)**2/RB12*RNUE/(1.D0+RC12*RNUE*EPSS))
C         RK22E=RK22*(1.D0/(1.D0+RA22*SQRT(RNUE)+RB22*RNUE)
C     &              +(EPSS*RC22)**2/RB22*RNUE/(1.D0+RC22*RNUE*EPSS))
         RK13E=RK13/(1.D0+RA13*SQRT(RNUE)+RB13*RNUE)
     &             /(1.D0+RC13*RNUE*EPSS)
         RK23E=RK23/(1.D0+RA23*SQRT(RNUE)+RB23*RNUE)
     &             /(1.D0+RC23*RNUE*EPSS)
C         RK33E=RK33/(1.D0+RA33*SQRT(RNUE)+RB33*RNUE)
C     &             /(1.D0+RC33*RNUE*EPSS)
C
C         RK2D =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUD)+RB2 *RNUD)
C     &              +(EPSS*RC2 )**2/RB2 *RNUD/(1.D0+RC2 *RNUD*EPSS))
C         RK2T =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUT)+RB2 *RNUT)
C     &              +(EPSS*RC2 )**2/RB2 *RNUT/(1.D0+RC2 *RNUT*EPSS))
C         RK2A =RK2 *(1.D0/(1.D0+RA2 *SQRT(RNUA)+RB2 *RNUA)
C     &              +(EPSS*RC2 )**2/RB2 *RNUA/(1.D0+RC2 *RNUA*EPSS))
         RK3D=((1.17D0-0.35D0*SQRT(RNUD))/(1.D0+0.7D0*SQRT(RNUD))
     &         -2.1D0*(RNUD*EPSS)**2)/(1.D0+(RNUD*EPSS)**2)
         RK3T=((1.17D0-0.35D0*SQRT(RNUT))/(1.D0+0.7D0*SQRT(RNUT))
     &         -2.1D0*(RNUT*EPSS)**2)/(1.D0+(RNUT*EPSS)**2)
         RK3A=((1.17D0-0.35D0*SQRT(RNUA))/(1.D0+0.7D0*SQRT(RNUA))
     &         -2.1D0*(RNUA*EPSS)**2)/(1.D0+(RNUA*EPSS)**2)
C
         DRL=RJCB(NR)/DR
         DTE=DERIV3P(PTS(1),RT(NR,1),RT(NR-1,1),
     &               RHOG(NR),RHOM(NR),RHOM(NR-1))
         DTD=DERIV3P(PTS(2),RT(NR,2),RT(NR-1,2),
     &               RHOG(NR),RHOM(NR),RHOM(NR-1))
         DTT=DERIV3P(PTS(3),RT(NR,3),RT(NR-1,3),
     &               RHOG(NR),RHOM(NR),RHOM(NR-1))
         DTA=DERIV3P(PTS(4),RT(NR,4),RT(NR-1,4),
     &               RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPE=DERIV3P(PNSS(1)*PTS(1),
     &               RN(NR,1)*RT(NR,1),
     &               RN(NR-1,1)*RT(NR-1,1),
     &               RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPD=DERIV3P(PNSS(2)*PTS(2),
     &               RN(NR,2)*RT(NR,2),
     &               RN(NR-1,2)*RT(NR-1,2),
     &               RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPT=DERIV3P(PNSS(3)*PTS(3),
     &               RN(NR,3)*RT(NR,3),
     &               RN(NR-1,3)*RT(NR-1,3),
     &               RHOG(NR),RHOM(NR),RHOM(NR-1))
         DPA=DERIV3P(PNSS(4)*PTS(4),
     &               RN(NR,4)*RT(NR,4),
     &               RN(NR-1,4)*RT(NR-1,4),
     &               RHOG(NR),RHOM(NR),RHOM(NR-1))
         BPL=BP(NR)
C
         FACT=1.D0/(1.D0+(RNUE*EPS)**2)
         IF(NSMAX.EQ.2) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL
     &       +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
         ELSEIF(NSMAX.EQ.3) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL
     &       +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
     &       +(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)
         ELSEIF(NSMAX.EQ.4) THEN
            A=           DPE/PEL-    2.5D0*DTE/TEL
     &       +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
     &       +(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)
     &       +(PAL/PEL)*(DPA/PAL-RK3A*FACT*DTA/TAL)
         ENDIF
         H=BB/(BB+BP(NR))
         AJBSL(NR)=-PBSCD*H*SQRT(EPS)*ANE*1.D20*TEL*RKEV/BPL
     &            *(RK13E*A+RK23E*DTE/TEL)
C     
C        WRITE(6,'(I3,1P6E12.4)') NR,RNUE,RK13E,RK23E,A,BPL,AJBSL(NR)
C         IF(NR.EQ.NRMAX) THEN
C            WRITE(6,'(1P6E12.4)') DN,DTE,ANE,PNSS(1),TE,PTS(1)
C            WRITE(6,'(1P6E12.4)') RK13E*(RN(NR,1)*RT(NR,1)
C     &                   +RN(NR,2)*RT(NR,2)
C     &                   +RN(NR,3)*RT(NR,3)
C     &                   +RN(NR,4)*RT(NR,4)
C     &                   +RW(NR,1)
C     &                   +RW(NR,2)         )*DN/ANE,
C     &            +(RK23E-1.5D0*RK13E)*RN(NR,1)*DTE,
C     &            +RK13D*(RK23D-1.5D0)*RN(NR,2)*DTD,
C     &            BPL,AJBS(NR-1),AJBS(NR)
C         ENDIF
C
      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO
C
      RETURN
      END
C
C     ***********************************************************
C
C           CALCULATE AJ, AJOH, POH, EZOH
C
C     ***********************************************************
C
      SUBROUTINE TRAJOH
C
      INCLUDE 'trcomm.inc'
C
      IF(MDLEQB.EQ.1.OR.MDLJQ.EQ.1.OR.(MDLUF.EQ.0.OR.MDLUF.EQ.3)) THEN
      NR=1
         FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
         FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
         AJ(NR) =FACTOR0*FACTORP*RDP(NR)/DR
      DO NR=2,NRMAX
         FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
         FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
         FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
         AJ(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
      ENDDO
      NR=1
         FACTOR0=RR/(RMU0*DVRHO(NR))
         FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
         AJTOR(NR) =FACTOR0*FACTORP*RDP(NR)/DR
      DO NR=2,NRMAX
         FACTOR0=RR/(RMU0*DVRHO(NR))
         FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
         FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
         AJTOR(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
      ENDDO
      ENDIF
C
      DO NR=1,NRMAX
         AJOH(NR) = AJ(NR)-(AJNB(NR)+AJRF(NR)+AJBS(NR))
         EZOH(NR) = ETA(NR)*AJOH(NR)
C         write(6,'(I3,3F17.7)') NR,AJ(NR),AJBS(NR)
         IF(KUFDEV.EQ.'lhd') THEN
            POH(NR)  = 0.D0
         ELSE
            POH(NR)  = EZOH(NR)*AJOH(NR)
         ENDIF
      ENDDO
C
      RETURN
      END
C
C     ***********************************************************
C
C           SAWTOOTH OSCILLATION
C
C     ***********************************************************
C
      SUBROUTINE TRSAWT
C
      INCLUDE 'trcomm.inc'
C
      DIMENSION QONE(NRM)
C
      IF(MDLST.EQ.0) RETURN
C
      IF(MOD(MDLST,2).EQ.1) THEN
         LQ = 1
      ELSE
         LQ = 0
      ENDIF
      IF(MOD(MDLST/2,2).EQ.1) THEN
         LT = 1
      ELSE
         LT = 0
      ENDIF
      IF(MOD(MDLST/4,2).EQ.1) THEN
         LN = 1
      ELSE
         LN = 0
      ENDIF
C
      WRITE(6,601) MDLST,T
  601 FORMAT(' ','# SAWTOOTH OSCILLATION -TYPE ',I1,
     &           ' AT ',F7.3,' SEC')
C
      IF(MDLST.EQ.8) THEN
         LQ = 0
         LT = 1
         LN = 1
         MDLST = 0
      ENDIF
C
      SUM=0.D0
      DO NR=1,NRMAX
         SUM = SUM+(1.D0/ABS(QP(NR))-1.D0)*(DVRHOG(NR)/2.D0*PI*RR)*DR
         IF(SUM.LT.0.D0) GOTO 1000
      ENDDO
      NR=NRMAX
C
 1000 IZEROX=NR
C
      DO NR=1,NRMAX
         IF(QP(NR).GE.1.D0) GOTO 2000
      ENDDO
      NR=NRMAX
C
 2000 IONE=NR
C
      SUM = 0.D0
      DO NR=IONE,IZEROX
         SUM = SUM+(1.D0/ABS(QP(NR))-1.D0)*(DVRHOG(NR)/2.D0*PI*RR)*DR
      ENDDO
C
      DO NR=1,IZEROX
         QONE(NR) = 1.D0+SUM*4.D0*RG(NR)**2/RG(IZEROX)**4
      ENDDO
C
      IF(LT.EQ.1) THEN
         DO NS=1,NSM
            SUM1 = 0.D0
            SUM2 = 0.D0
            DO NR=1,IZEROX
               SUM1 = SUM1+RN(NR,NS)          *DVRHO(NR)
               SUM2 = SUM2+RN(NR,NS)*RT(NR,NS)*DVRHO(NR)
            ENDDO
            RTN = SUM2/SUM1
C
            DO NR=1,IZEROX
               RT(NR,NS) = RTN
            ENDDO
         ENDDO
      ENDIF
C
      IF(LN.EQ.1) THEN
         DO NS=1,NSM
            SUM1 = 0.D0
            SUM2 = 0.D0
            DO NR=1,IZEROX
               SUM1 = SUM1+          DVRHO(NR)
               SUM2 = SUM2+RN(NR,NS)*DVRHO(NR)
            ENDDO
            RNN = SUM2/SUM1
            DO NR=1,IZEROX
               RN(NR,NS) = RNN
            ENDDO
         ENDDO
      ENDIF
C
      IF(LQ.EQ.1) THEN
         DO NR=1,IZEROX
            QP(NR) = 1.D0/QONE(NR)
            RDP(NR)=TTRHOG(NR)*ARRHOG(NR)*DVRHOG(NR)/(4.D0*PI**2*QP(NR))
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
      ENDIF
C
         WRITE(6,602) RM(IONE),RM(IZEROX),RTN,RNN
  602    FORMAT(' ',' R-ONE,R-ZERO,RTN,RNN = ',4F8.3)
      RETURN
      END
C
C     ***********************************************************
C
C           TRAPPED PARTICLE FRACTION
C
C     ***********************************************************
C
      REAL*8 FUNCTION FTPF(ID,EPS)
C
      IMPLICIT NONE
      INTEGER ID,IMAX,N,IERR
      PARAMETER (IMAX=20)
      REAL*8 EPS,PI,FTUL,FTLL,EPSC,S,OMEGA
      REAL*8 TABLE(IMAX,IMAX)
      EXTERNAL FTU,FTL
C
      IF(ID.EQ.1) THEN
C  Y. R. Lin-Liu and R. L. Miller, PoP 2 1666 (1995), eqs(7)(13)(18)(19)
         PI=3.14159265358979323846D0
         EPSC=1.D-9
         FTUL=1.D0-(1.D0-1.5D0*SQRT(EPS)+0.5D0*EPS**1.5D0)/SQRT(1-EPS**2
     &        )
         CALL RMBRG(0.D0,2.D0*PI,EPSC,S,IMAX,N,IERR,TABLE,EPS,FTL)
         FTLL=1.D0-(1.D0-EPS)**1.5D0/SQRT(1.D0+EPS)*(S/(2.D0*PI))
         OMEGA=(3.D0*SQRT(2.D0)/2.D0*0.69D0-3.D0*SQRT(2.D0)/PI)
     &        /(1.5D0-3.D0*SQRT(2.D0)/PI)
         FTPF=OMEGA*FTUL+(1.D0-OMEGA)*FTLL
      ELSEIF(ID.EQ.2) THEN
C  S. P. Hirshman et al., NF 17 611 (1977)
         FTPF=1.D0-(1.D0-EPS)**2.D0
     &        /(DSQRT(1.D0-EPS**2)*(1.D0+1.46D0*DSQRT(EPS)))
      ELSEIF(ID.EQ.3) THEN
C  Y. R. Lin-Liu and R. L. Miller, PoP 2 1666 (1995), eqs(16)(17)(18)
         PI=3.14159265358979323846D0
         FTUL=1.5D0*SQRT(EPS)
         FTLL=3.D0*SQRT(2.D0)/PI*SQRT(EPS)
         FTPF=0.75D0*FTUL+0.25D0*FTLL
      ELSEIF(ID.EQ.4) THEN
C  M. N. Rosenbluth et al., PoF 15 116 (1972)
         FTPF=1.46D0*SQRT(EPS)
      ELSE
C  Y. B. Kim et al., PoF B 3 2050 (1991) eq(C18), default
         FTPF=1.46D0*SQRT(EPS)-0.46D0*(EPS)**1.5D0
      ENDIF
C
      RETURN
      END
C
C     *********************************************
C
C           FUNCTION FOR ROMBERG INTEGRATION
C
C     *********************************************
C
      REAL*8 FUNCTION FTU(X,EPS)
C
      IMPLICIT NONE
      REAL*8 X, EPS
C
      FTU = X/SQRT(1.D0-X*(1.D0-EPS))
C
      RETURN
      END
C
      REAL*8 FUNCTION FTL(X,EPS)
C    
      IMPLICIT NONE
      REAL*8 X, EPS
      REAL*8 H
C
      H = (1.D0 - EPS) / (1.D0 + EPS * COS(X))
      FTL = (1.D0 - SQRT(1.D0 - H) * (1.D0 + 0.5D0 * H)) / H**2
C
      RETURN
      END
C      
C     *********************************************
C
C           ROMBERG INTEGRATION METHOD     
C
C     *********************************************
C
      SUBROUTINE RMBRG(A,B,EPS,S,IMAX,N,IERR,T,ARG,F)
C
C     <input>
C        A     : lower bound
C        B     : upper bound
C        EPS   : stopping criterion
C        IMAX  : maximum division number
C        F     : formula of integrand
C     <output>
C        S     : integration value
C        N     : division number
C        IERR  : error indicator
C        T     : Romberg T table
C
      IMPLICIT NONE
      INTEGER IMAX,N,IERR,K,N2,I,J
      REAL*8 A,B,EPS,S,S1,H,Y1,Y2,F,X,ARG
      REAL*8 T(IMAX,IMAX)
      EXTERNAL F
C
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
C
      RETURN
      END
C
C     ***********************************************************
C
C           COULOMB LOGARITHM
C
C     ***********************************************************
C
      REAL*8 FUNCTION COULOG(NS1,NS2,ANEL,TL)
C
C     ANEL : electron density [10^20 /m^3]
C     TL   : electron or ion temperature [keV]
C            in case of ion-ion collision, TL becomes ion temp.
C
      IMPLICIT NONE
      INTEGER NS1,NS2
      REAL*8 ANEL,TL
C
      IF(NS1.EQ.1.AND.NS2.EQ.1) THEN
         COULOG=14.9D0-0.5D0*LOG(ANEL)+LOG(TL)
      ELSE
         IF(NS1.EQ.1.OR.NS2.EQ.1) THEN
            COULOG=15.2D0-0.5D0*LOG(ANEL)+LOG(TL)
         ELSE
            COULOG=17.3D0-0.5D0*LOG(ANEL)+1.5D0*LOG(TL)
         ENDIF
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
C           COLLISION TIME 
C
C     ***********************************************************
C
C     between electrons and ions
C
      FUNCTION FTAUE(ANEL,ANIL,TEL,ZL)
C
C     ANEL : electron density [10^20 /m^3]
C     ANIL : ion density [10^20 /m^3]
C     TEL  : electron temperature [kev]
C     ZL   : ion charge number
C
      INCLUDE 'trcomm.inc'
C
      COEF = 6.D0*PI*SQRT(2.D0*PI)*EPS0**2*SQRT(AME)/(AEE**4*1.D20)
      IF(ZL-PZ(2).LE.1.D-7) THEN
         FTAUE = COEF*(TEL*RKEV)**1.5D0
     &          /(ANIL*ZL**2*COULOG(1,2,ANEL,TEL))
      ELSE
C     If the plasma contains impurities, we need to consider the
C     effective charge number instead of ion charge number.
C     From the definition of Zeff=sum(n_iZ_i^2)/n_e, 
C     n_iZ_i^2 is replaced by n_eZ_eff at the denominator of tau_e.
         FTAUE = COEF*(TEL*RKEV)**1.5D0
     &          /(ANEL*ZL*COULOG(1,2,ANEL,TEL))
      ENDIF
C
      RETURN
      END
C
C     between ions and ions
C
      FUNCTION FTAUI(ANEL,ANIL,TIL,ZL,PAL)
C
C     ANEL : electron density [10^20 /m^3]
C     ANIL : ion density [10^20 /m^3]
C     TIL  : ion temperature [kev]
C     ZL   : ion charge number
C     PAL  : ion atomic number
C
      INCLUDE 'trcomm.inc'
C
      COEF = 12.D0*PI*SQRT(PI)*EPS0**2*SQRT(PAL*AMM)/(AEE**4*1.D20)
      FTAUI = COEF*(TIL*RKEV)**1.5D0/(ANIL*ZL**4*COULOG(2,2,ANEL,TIL))
C
      RETURN
      END
