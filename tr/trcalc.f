C     $Id$
C  
C     ***********************************************************
C
C           CALCULATE TRANSPORT COEFFICIENTS AND SOURCE
C
C     ***********************************************************
C
      SUBROUTINE TRCALC
C
      INCLUDE 'trcomm.h'
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
         AJNB(NR)=0.D0
         AJRF(NR)=0.D0
         AJBS(NR)=0.D0
      DO NS=1,NSM
         SPE(NR,NS)=0.D0
         PRF(NR,NS)=0.D0
         PBCL(NR,NS)=0.D0
         PFCL(NR,NS)=0.D0
      ENDDO
      ENDDO
C
      IF(MODELG.EQ.3) THEN
         DO NR=1,NRMAX
            QP(NR)=QRHO(NR)*BPRHO(NR)/BP(NR)
         ENDDO
      ELSE
         DO NR=1,NRMAX
            QP(NR)=FKAP*RG(NR)*RA*BB/(RR*BP(NR))
CCC            QP(NR)=RG(NR)*RA*BB/(RR*BP(NR))
C            if(nr.le.7) write(6,*) NR,BP(NR)
         ENDDO
      ENDIF
      Q0  = (4.D0*QP(1) -QP(2) )/3.D0
C
      IF(T.LT.PELTIM+0.5D0*DT.AND.
     &   T.GE.PELTIM-0.5D0*DT) THEN
         CALL TRPELT
      ENDIF
      CALL TRZEFF
C
      IF(MDNCLS.NE.0) CALL TR_NCLASS
C
      CALL TRCOEF
      CALL TRLOSS
      CALL TRPWRF
      CALL TRPWNB
C
      IF(MDNCLS.NE.0) THEN
         CALL TRAJBS_NCLASS
      ELSE
         IF(MDLJBS.EQ.1) THEN
            CALL OLDTRAJBS
         ELSEIF(MDLJBS.EQ.2) THEN
            CALL TRAJBSO
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
C      IF(MDLUF.EQ.3) THEN
C         DO NR=1,NRMAX
C            AJBS(NR)=0.D0
C         ENDDO
C      ENDIF
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
      RETURN
      END
C
C     ***********************************************************
C
C           CALCULATE Z-C
C
C     ***********************************************************
C
      FUNCTION TRZEC(TE)
C
      IMPLICIT REAL*8 (A-H,O-Z)
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
      FUNCTION TRZEFE(TE)
C
      IMPLICIT REAL*8 (A-H,O-Z)
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
      INCLUDE 'trcomm.h'
C
      IF(MDLUF.EQ.3) THEN
         DO NR=1,NRMAX
            ZEFF(NR)=2.D0
         ENDDO
      ELSE
         DO NR=1,NRMAX
            TE =RT(NR,1)
            ZEFF(NR) =(PZ(2)  *PZ(2)  *RN(NR,2)
     &                +PZ(3)  *PZ(3)  *RN(NR,3)
     &                +PZ(4)  *PZ(4)  *RN(NR,4)
     &                +TRZEC(TE)**2   *ANC (NR)
     &                +TRZEFE(TE)**2  *ANFE(NR))/RN(NR,1)
         ENDDO
      ENDIF
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
      FUNCTION TRRPC(TE)
C
      IMPLICIT REAL*8 (A-H,O-Z)
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
      FUNCTION TRRPFE(TE)
C
      IMPLICIT REAL*8 (A-H,O-Z)
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
      INCLUDE 'trcomm.h'
C
      DO NR=1,NRMAX
         ANE =RN(NR,1)
         ANDX=RN(NR,2)
         ANT =RN(NR,3)
         ANHE=RN(NR,4)
         TE  =RT(NR,1)
         PLFE  = ANE*ANFE(NR)*TRRPFE(TE)*1.D40
         PLC   = ANE*ANC (NR)*TRRPC (TE)*1.D40
         PLD   = ANE*ANDX*5.35D-37*1.D0**2*SQRT(ABS(TE))*1.D40
         PLTT  = ANE*ANT *5.35D-37*1.D0**2*SQRT(ABS(TE))*1.D40
         PLHE  = ANE*ANHE*5.35D-37*2.D0**2*SQRT(ABS(TE))*1.D40
         PRL(NR)= PLFE+PLC+PLD+PLTT+PLHE
      ENDDO
C
C     ****** IONIZATION LOSS ******
C
      DO NR=1,NRMAX
         ANE=RN(NR,1)
         TE =RT(NR,1)
         EION  = 13.64D0
         TN    = MAX(TE*1.D3/EION,1.D-2)
         SION  = 1.D-11*SQRT(TN)*EXP(-1.D0/TN)
     &          /(EION**1.5D0*(6.D0+TN))
         PIE(NR) = ANE*ANNU(NR)*SION*1.D40*EION*AEE
         SIE(NR) = ANE*ANNU(NR)*SION*1.D20
         TSIE(NR)= ANE*SION*1.D20
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
         SCX(NR)=ANE*ANNU(NR)*SCH*1.D20
         TSCX(NR)=ANE*SCH*1.D20
C
         PCX(NR)=(-1.5D0*ANE*ANNU(NR)*SION*TNU
     &         +  1.5D0*ANDX*ANNU(NR)*SCH*(TD-TNU))*RKEV*1.D40
      ENDDO
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
      INCLUDE 'trcomm.h'
      DIMENSION AJBSL(NRM)
C
      IF(PBSCD.LE.0.D0) RETURN
C
      NSW=1
      IF(NSW.EQ.0) THEN
         DO NR=1,NRMAX
            AJBS(NR)=PBSCD*AJBSNC(NR)
         ENDDO
      ELSE
         DO NR=1,NRMAX-1
            SUM=0.D0
            DO NS=1,NSMAX
               RTNW=0.5D0*(RT(NR+1,NS)+RT(NR,NS))
               RPNW=0.5D0*(RN(NR+1,NS)*RT(NR+1,NS)
     &                    +RN(NR  ,NS)*RT(NR  ,NS))
               DRTNW=(RT(NR+1,NS)-RT(NR,NS))*AR1RHO(NR)/DR
               DRPNW=(RN(NR+1,NS)*RT(NR+1,NS)-RN(NR,NS)*RT(NR,NS))
     &              *AR1RHO(NR)/DR
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
            DRTNW=2.D0*(PTS(NS)-RT(NR,NS))*AR1RHO(NR)/DR
            DRPNW=2.D0*(PNSS(NS)*PTS(NR)-RN(NR,NS)*RT(NR,NS))
     &                *AR1RHO(NR)/DR
            SUM=SUM+CJBST(NR,NS)*DRTNW/RTNW
     &             +CJBSP(NR,NS)*DRPNW/RPNW
         ENDDO
         AJBSL(NR)=-PBSCD*SUM/BB
C     
         AJBS(1)=0.5D0*AJBSL(1)
         DO NR=2,NRMAX
            AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
         ENDDO
      ENDIF
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
      INCLUDE 'trcomm.h'
      DIMENSION ANI(NRM),AJBSL(NRM)
C
      IF(PBSCD.LE.0.D0) RETURN
C
      DO NR=1,NRMAX-1
C
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         QL=ABS(QP(NR))
         ZEFFL=0.5D0*(ZEFF(NR-1)+ZEFF(NR))
C
         RNTP=0.D0
         RNP =0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+RN(NR+1,NS)*RT(NR+1,NS)
            RNP =RNP +RN(NR+1,NS)
            RNTM=RNTM+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR  ,NS)
         ENDDO
         RNTP=RNTP+RW(NR+1,1)+RW(NR+1,2)
         RNTM=RNTM+RW(NR  ,1)+RW(NR  ,2)
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
         DO NS=2,NSMAX
            ANI(NR)=ANI(NR)+0.5D0*PZ(NS)*(RN(NR+1,NS)+RN(NR,NS))
         ENDDO
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         PPI=0.5D0*(RNTP+RNTM)
         DTI=(RNTP/RNP-RNTM/RNM)*AR1RHO(NR)/DR
         DPI=(RNTP-RNTM)*AR1RHO(NR)/DR
C
         rLnLamii=30.D0-LOG(PZ(2)**3*SQRT(ANI(NR)*1.D20)
     &           /(ABS(TI*1.D3)**1.5))
         RNUI=4.90D-18*QL*RR*ANI(NR)*1.D20*PZ(2)**4*rLnLamii
     &       /(ABS(TI*1.D3)**2*EPSS)
C     
C     ****** ELECTORON PARAMETER ******
C
C     *** ANE  is the the electron density (ne) ***
C     *** TE   is the electron temperature (Te) ***
C     *** PE   is the electron pressure (Pe) ***
C     *** DNE  is the derivative of electron density (dne/dr) ***
C     *** DTE  is the derivative of electron temperature (dTe/dr) ***
C     *** DPE  is the derivative of electron pressure (dPe/dr) ***
C
         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
         TE= 0.5D0*(RT(NR+1,1)+RT(NR,1))
         PE= 0.5D0*(RN(NR+1,1)*RT(NR+1,1)+RN(NR,1)*RT(NR,1))
         DNE=(RN(NR+1,1)-RN(NR,1))*AR1RHO(NR)/DR
         DTE=(RT(NR+1,1)-RT(NR,1))*AR1RHO(NR)/DR
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*AR1RHO(NR)/DR
C
         rLnLame=31.3D0-LOG(SQRT(ANE*1.D20)/ABS(TE*1.D3))
         RNUE=6.921D-18*QL*RR*ANE*1.D20*ZEFFL*rLnLame
     &       /(ABS(TE*1.D3)**2*EPSS)
C
         RPE=PE/(PE+PPI)
C     <1>
         FT1=1.D0-(1.D0-EPS)**2.D0
     &         /(DSQRT(1.D0-EPS**2)*(1.D0+1.46D0*DSQRT(EPS)))
C     <2>
         FT2=1.D0-(1.D0-1.5D0*SQRT(EPS)+0.5D0*EPSS)/SQRT(1-EPS**2)
C     <3>
         FT3=1.46D0*SQRT(EPS)
C     <4>
         FTU=1.5D0*SQRT(EPS)
         FTL=3.D0*SQRT(2.D0)/PI*SQRT(EPS)
         FT4=0.75D0*FTU+0.25D0*FTL
         FT=FT1
C         write(6,'I2,4E20.12') NR,FT1,FT2,FT3,FT4
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
C         AJBSL(NR-1)=-PBSCD*PE*1.D20*RKEV
C     &             *( RL31*(DPE/PE+DPI/PE)+RL32*DTE/TE
C     &               +RL34*SALFA*(1.D0-RPE)/RPE*DTI/TI)/BP(NR-1)
         AJBSL(NR)=-PBSCD*(PE+PPI)*1.D20*RKEV
     &            *( RL31*DNE/ANE
     &              +RPE*(RL31+RL32)*DTE/TE
     &             +(1.D0-RPE)*(1.D0+RL34/RL31*SALFA)*RL31*DTI/TI)
     &            /BP(NR)
      ENDDO
C
      NR=NRMAX
      EPS=EPSRHO(NR)
      EPSS=SQRT(EPS)**3
      QL=ABS(QP(NR))
      ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
C
      RNTP=0.D0
      RNP =0.D0
      RNTM=0.D0
      RNM =0.D0
      DO NS=2,NSMAX
         RNTP=RNTP+PNSS(NS)*PTS(NS)
         RNP =RNP +PNSS(NS)
         RNTM=RNTM+RN(NR,NS)*RT(NR,NS)
         RNM =RNM +RN(NR,NS)
      ENDDO
      RNTP=RNTP+RW(NR  ,1)+RW(NR  ,2)
      RNTM=RNTM+RW(NR  ,1)+RW(NR  ,2)
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
      DO NS=2,NSMAX
         ANI(NR)=ANI(NR)+PZ(NS)*PNSS(NS)
      ENDDO
      TI =RNTP/RNP
      PPI=RNTP
      DTI=2.D0*(RNTP/RNP-RNTM/RNM)*AR1RHO(NR)/DR
      DPI=2.D0*(RNTP-RNTM)*AR1RHO(NR)/DR
C
      rLnLamii=30.D0-LOG(PZ(2)**3*SQRT(ANI(NR)*1.D20)
     &     /(ABS(TI*1.D3)**1.5))
      RNUI=4.90D-18*QL*RR*ANI(NR)*1.D20*PZ(2)**4*rLnLamii
     &     /(ABS(TI*1.D3)**2*EPSS)
C     
C     ****** ELECTORON PARAMETER ******
C
C     *** ANE  is the the electron density (ne) ***
C     *** TE   is the electron temperature (Te) ***
C     *** PE   is the electron pressure (Pe) ***
C     *** DNE  is the derivative of electron density (dne/dr) ***
C     *** DTE  is the derivative of electron temperature (dTe/dr) ***
C     *** DPE  is the derivative of electron pressure (dPe/dr) ***
C
      ANE=PNSS(1)
      TE =PTS(1)
      PE =PNSS(1)*PTS(1)
      DNE=2.D0*(PNSS(1)-RN(NR,1))*AR1RHO(NR)/DR
      DTE=2.D0*(PTS (1)-RT(NR,1))*AR1RHO(NR)/DR
      DPE=2.D0*(PNSS(1)*PTS(1)-RN(NR,1)*RT(NR,1))*AR1RHO(NR)/DR
C     
      rLnLame=31.3D0-LOG(SQRT(ANE*1.D20)/ABS(TE*1.D3))
      RNUE=6.921D-18*QL*RR*ANE*1.D20*ZEFFL*rLnLame
     &     /(ABS(TE*1.D3)**2*EPSS)
C     
      RPE=PE/(PE+PPI)
C     <1>
      FT1=1.D0-(1.D0-EPS)**2.D0
     &     /(DSQRT(1.D0-EPS**2)*(1.D0+1.46D0*DSQRT(EPS)))
C     <2>
      FT2=1.D0-(1.D0-1.5D0*SQRT(EPS)+0.5D0*EPSS)/SQRT(1-EPS**2)
C     <3>
      FT3=1.46D0*SQRT(EPS)
C     <4>
      FTU=1.5D0*SQRT(EPS)
      FTL=3.D0*SQRT(2.D0)/PI*SQRT(EPS)
      FT4=0.75D0*FTU+0.25D0*FTL
      FT=FT1
C         write(6,'I2,4E20.12') NR,FT1,FT2,FT3,FT4
C
C     F33TEFF=FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE)
C     &          +0.45D0*(1.D0-FT)*RNUE/ZEFFL**1.5)
      F31TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)
     &     +0.5D0*(1.D0-FT)*RNUE/ZEFFL)
      F32EETEFF=FT/(1.D0+0.26D0*(1.D0-FT)*SQRT(RNUE)
     &     +0.18D0*(1.D0-0.37D0*FT)*RNUE/SQRT(ZEFFL))
      F32EITEFF=FT/(1.D0+(1.D0+0.6D0*FT)*SQRT(RNUE)
     &     +0.85D0*(1.D0-0.37D0*FT)*RNUE*(1.D0+ZEFFL))
      F34TEFF=FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(RNUE)
     &     +0.5D0*(1.D0-0.5D0*FT)*RNUE/ZEFFL)
C
      SALFA0=-1.17D0*(1.D0-FT)/(1.D0-0.22D0*FT-0.19D0*FT**2)
      SALFA=((SALFA0+0.25D0*(1.D0-FT**2)*SQRT(RNUI))
     &     /(1.D0+0.5D0*SQRT(RNUI))+0.315D0*RNUI**2*FT**6)
     &     /(1.D0+0.15D0*RNUI**2*FT**6)
C
C     RNZ=0.58D0+0.74D0/(0.76D0+ZEFFL)
C     SGMSPTZ=1.9012D4*(TE*1.D3)**1.5/(ZEFFL*RNZ*rLnLame)
C     SGMNEO=SGMSPTZ*F33(F33TEFF,ZEFFL)
C     
      RL31=F31(F31TEFF,ZEFFL)
      RL32=F32EE(F32EETEFF,ZEFFL)+F32EI(F32EITEFF,ZEFFL)
      RL34=F31(F34TEFF,ZEFFL)
C
C     AJBSL(NR)=-PBSCD*PE*1.D20*RKEV
C     &             *( RL31*(DPE/PE+DPI/PE)+RL32*DTE/TE
C     &               +RL34*SALFA*(1.D0-RPE)/RPE*DTI/TI)/BP(NR)
      AJBSL(NR)=-PBSCD*(PE+PPI)*1.D20*RKEV
     &                *( RL31*DNE/ANE
     &                  +RPE*(RL31+RL32)*DTE/TE
     &                  +(1.D0-RPE)*(1.D0+RL34/RL31*SALFA)*RL31*DTI/TI
     &                  )/BP(NR)
C
      AJBS(1)=0.5D0*AJBSL(1)
      DO NR=2,NRMAX
         AJBS(NR)=0.5D0*(AJBSL(NR)+AJBSL(NR-1))
      ENDDO
C
      RETURN
      END
C
C     *****
C
      FUNCTION F33(X,Z)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      F33 = 1.D0-(1.D0+0.36D0/Z)*X+0.59D0/Z*X**2-0.23D0/Z*X**3
C
      RETURN
      END
C
      FUNCTION F31(X,Z)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      F31 = (1.D0+1.4D0/(Z+1.D0))*X-1.9D0/(Z+1.D0)*X**2
     &     +0.3D0/(Z+1.D0)*X**3+0.2D0/(Z+1.D0)*X**4
C
      RETURN
      END
C
      FUNCTION F32EE(X,Z)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      F32EE = (0.05D0+0.62D0*Z)/(Z*(1.D0+0.44D0*Z))*(X-X**4)
     &       +1.D0/(1.D0+0.22D0*Z)*(X**2-X**4-1.2D0*(X**3-X**4))
     &       +1.2D0/(1.D0+0.5D0*Z)*X**4
C
      RETURN
      END
C
      FUNCTION F32EI(X,Z)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
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
C         BOOTSTRAP CURRENT (TOKAMAKS, H.R. Wilson)
C
C     ************************************************
C
      SUBROUTINE TRAJBSNEW
C
      INCLUDE 'trcomm.h'
      DIMENSION ANI(NRM),AJBSL(NRM)
C
      IF(PBSCD.LE.0.D0) RETURN
C
      DO NR=1,NRMAX-1
C
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         QL=ABS(QP(NR))
         ZEFFL=0.5D0*(ZEFF(NR-1)+ZEFF(NR))
C
         RNTP=0.D0
         RNP= 0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+RN(NR+1,NS)*RT(NR+1,NS)
            RNP =RNP +RN(NR+1,NS)
            RNTM=RNTM+RN(NR  ,NS)*RT(NR  ,NS)
            RNM =RNM +RN(NR  ,NS)
         ENDDO
         RNTP=RNTP+RW(NR+1,1)+RW(NR+1,2)
         RNTM=RNTM+RW(NR  ,1)+RW(NR  ,2)
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
         DO NS=2,NSMAX
            ANI(NR)=ANI(NR)+0.5D0*PZ(NS)*(RN(NR+1,NS)+RN(NR,NS))
         ENDDO
         TI =0.5D0*(RNTP/RNP+RNTM/RNM)
         PPI=0.5D0*(RNTP+RNTM)
         DTI=(RNTP/RNP-RNTM/RNM)*AR1RHO(NR)/DR
         DPI=(RNTP-RNTM)*AR1RHO(NR)/DR
         VTI=SQRT(ABS(TI)*RKEV/AMM)
C
         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
         rLnLam=17.3D0-DLOG(ANE)*0.5D0+DLOG(ABS(TI))*1.5D0
         TAUI=12.D0*PI*SQRT(PI)*AEPS0**2*SQRT(AMM)
     &             *(ABS(TI)*RKEV)**1.5D0/(ANI(NR)*1.D20
     &             *ZEFFL**4*AEE**4*rLnLam)
C
         RNUI=QL*RR/(TAUI*VTI*EPSS)
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
         DTE=(RT(NR+1,1)-RT(NR,1))*AR1RHO(NR)/DR
         DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR,1)*RT(NR,1))*AR1RHO(NR)/DR
         VTE=SQRT(ABS(TE)*RKEV/AME)
C
         rLnLam=15.2D0-DLOG(ANE)*0.5D0+DLOG(ABS(TE))
         TAUE=6.D0*PI*SQRT(2.D0*PI)*AEPS0**2*SQRT(AME)
     &             *(ABS(TE)*RKEV)**1.5D0/(ANE*1.D20
     &             *ZEFFL**2*AEE**4*rLnLam)
C
         RNUE=QL*RR/(TAUE*VTE*EPSS)
C
         FT=1.D0-(1.D0-EPS)**2.D0
     &         /(DSQRT(1.D0-EPS**2)*(1.D0+1.46D0*DSQRT(EPS)))
C         FT=1.46D0*SQRT(EPS)
         DDX=2.4D0+5.4D0*FT+2.6D0*FT**2
C
         DDD=-1.17D0/(1.D0+0.46D0*FT)
C
         C1=(4.D0+2.6D0*FT)
     &      /((1.D0+1.02D0*DSQRT(RNUE)+1.07D0*RNUE)
     &       *(1.D0+1.07D0*EPSS*RNUE))
         C2=C1*TI/TE
         C3=(7.D0+6.5D0*FT)
     &      /((1.D0+0.57D0*DSQRT(RNUE)+0.61D0*RNUE)
     &       *(1.D0+0.61D0*EPSS*RNUE))
     &      -2.5D0*C1
         C4=((DDD+0.35D0*DSQRT(RNUI))
     &      /(1.D0+0.7D0*DSQRT(RNUI))
     &      +2.1D0*EPS**3*RNUI**2)*C2
     &      /((1.D0+EPS**3*RNUI**2)
     &       *(1.D0+EPS**3*RNUE**2))
C
         AJBSL(NR)=-PBSCD*FT*PE*1.D20*RKEV/(DDX*BP(NR))
     &        *(C1*(DPE/PE)
     &         +C2*(DPI/PPI)
     &         +C3*(DTE/TE)
     &         +C4*(DTI/TI))
C
C     *** H.R. WILSON, Nucl.Fusion 32,no.2,1992 259-263 ***
C
C         DDX=1.414D0*ZEFFL+ZEFFL**2+FT*(0.754D0+2.657D0*ZEFFL
C     &        +2.D0*ZEFFL**2)+FT**2*(0.348D0+1.243D0*ZEFFL+ZEFFL**2)
C         RL31=FT*( 0.754D0+2.21D0*ZEFFL+ZEFFL**2+FT*(0.348D0+1.243D0
C     &            *ZEFFL+ZEFFL**2))/DDX
C         RL32=-FT*(0.884D0+2.074D0*ZEFFL)/DDX
C         DDD=-1.172D0/(1.D0+0.462D0*FT)
C
C         AJBSL(NR)=-PBSCD*PE*1.D20*RKEV*(RL31*(DPE/PE)+(TI/(ZEFFL*TE))
C     &        *((DPI/PPI)+DDD*(DTI/TI))+RL32*(DTE/TE))/BP(NR)
      ENDDO
C
      NR=NRMAX
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         QL=ABS(QP(NR))
         ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
C
         RNTP=0.D0
         RNP= 0.D0
         RNTM=0.D0
         RNM =0.D0
         DO NS=2,NSMAX
            RNTP=RNTP+PNSS(NS)*PTS(NS)
            RNP =RNP +PNSS(NS)
            RNTM=RNTM+RN(NR,NS)*RT(NR,NS)
            RNM =RNM +RN(NR,NS)
         ENDDO
         RNTP=RNTP+RW(NR  ,1)+RW(NR  ,2)
         RNTM=RNTM+RW(NR  ,1)+RW(NR  ,2)
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
         DO NS=2,NSMAX
            ANI(NR)=ANI(NR)+PZ(NS)*PNSS(NS)
         ENDDO
         TI =RNTP/RNP
         PPI=RNTP
         DTI=2.D0*(RNTP/RNP-RNTM/RNM)*AR1RHO(NR)/DR
         DPI=2.D0*(RNTP-RNTM)*AR1RHO(NR)/DR
         VTI=SQRT(ABS(TI)*RKEV/AMM)
C
         ANE=PNSS(1)
         rLnLam=17.3D0-DLOG(ANE)*0.5D0+DLOG(ABS(TI))*1.5D0
         TAUI=12.D0*PI*SQRT(PI)*AEPS0**2*SQRT(AMM)
     &             *(ABS(TI)*RKEV)**1.5D0/(ANI(NR)*1.D20
     &             *ZEFFL**4*AEE**4*rLnLam)
C
         RNUI=QL*RR/(TAUI*VTI*EPSS)
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
         DTE=2.D0*(PTS(1)-RT(NR,1))*AR1RHO(NR)/DR
         DPE=2.D0*(PNSS(1)*PTS(1)-RN(NR,1)*RT(NR,1))*AR1RHO(NR)/DR
         VTE=SQRT(ABS(TE)*RKEV/AME)
C
         rLnLam=15.2D0-DLOG(ANE)*0.5D0+DLOG(ABS(TE))
         TAUE=6.D0*PI*SQRT(2.D0*PI)*AEPS0**2*SQRT(AME)
     &             *(ABS(TE)*RKEV)**1.5D0/(ANE*1.D20
     &             *ZEFFL**2*AEE**4*rLnLam)
C
         RNUE=QL*RR/(TAUE*VTE*EPSS)
C
         FT=1.D0-(1.D0-EPS)**2.D0
     &         /(DSQRT(1.D0-EPS**2)*(1.D0+1.46D0*DSQRT(EPS)))
C         FT=1.46D0*SQRT(EPS)
         DDX=2.4D0+5.4D0*FT+2.6D0*FT**2
C
         DDD=-1.17D0/(1.D0+0.46D0*FT)
C
         C1=(4.D0+2.6D0*FT)
     &      /((1.D0+1.02D0*DSQRT(RNUE)+1.07D0*RNUE)
     &       *(1.D0+1.07D0*EPSS*RNUE))
         C2=C1*TI/TE
         C3=(7.D0+6.5D0*FT)
     &      /((1.D0+0.57D0*DSQRT(RNUE)+0.61D0*RNUE)
     &       *(1.D0+0.61D0*EPSS*RNUE))
     &      -2.5D0*C1
         C4=((DDD+0.35D0*DSQRT(RNUI))
     &      /(1.D0+0.7D0*DSQRT(RNUI))
     &      +2.1D0*EPS**3*RNUI**2)*C2
     &      /((1.D0+EPS**3*RNUI**2)
     &       *(1.D0+EPS**3*RNUE**2))
C
         AJBSL(NR)=-PBSCD*FT*PE*1.D20*RKEV/(DDX*BP(NR))
     &        *(C1*(DPE/PE)
     &         +C2*(DPI/PPI)
     &         +C3*(DTE/TE)
     &         +C4*(DTI/TI))
C
C     *** H.R. WILSON, Nucl.Fusion 32,no.2,1992 259-263 ***
C
C         DDX=1.414D0*ZEFFL+ZEFFL**2+FT*(0.754D0+2.657D0*ZEFFL
C     &        +2.D0*ZEFFL**2)+FT**2*(0.348D0+1.243D0*ZEFFL+ZEFFL**2)
C         RL31=FT*( 0.754D0+2.21D0*ZEFFL+ZEFFL**2+FT*(0.348D0+1.243D0
C     &            *ZEFFL+ZEFFL**2))/DDX
C         RL32=-FT*(0.884D0+2.074D0*ZEFFL)/DDX
C         DDD=-1.172D0/(1.D0+0.462D0*FT)
C
C         AJBSL(NR)=-PBSCD*PE*1.D20*RKEV*(RL31*(DPE/PE)+(TI/(ZEFFL*TE))
C     &        *((DPI/PPI)+DDD*(DTI/TI))+RL32*(DTE/TE))/BP(NR)
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
C           BOOTSTRAP CURRENT (HINTON & HAZELTINE)
C
C     ***********************************************************
C
      SUBROUTINE TRAJBS
C
      INCLUDE 'trcomm.h'
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
      DO 100 NR=1,NRMAX-1
C
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         ANE=0.5D0*(RN(NR+1,1)+RN(NR,1))
C         ANDX=0.5D0*(RN(NR+1,2)+RN(NR,2))
C         ANT =0.5D0*(RN(NR+1,3)+RN(NR,3))
C         ANA =0.5D0*(RN(NR+1,4)+RN(NR,4))
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
         COEF = 12.D0*PI*SQRT(PI)*AEPS0**2
     &         /(ANE*1.D20*ZEFFL*AEE**4*15.D0)
         TAUE = COEF*SQRT(AME)*(TEL*RKEV)**1.5D0/SQRT(2.D0)
         TAUD = COEF*SQRT(AMD)*(TDL*RKEV)**1.5D0/PZ(2)**2
         TAUT = COEF*SQRT(AMT)*(TTL*RKEV)**1.5D0/PZ(3)**2
         TAUA = COEF*SQRT(AMA)*(TAL*RKEV)**1.5D0/PZ(4)**2
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
         DRL=AR1RHO(NR)/DR
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
         AJBSL(NR)=-PBSCD*SQRT(EPS)*ANE*1.D20*TEL*RKEV/BPL
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
  100 CONTINUE
C
      NR=NRMAX
         EPS=EPSRHO(NR)
         EPSS=SQRT(EPS)**3
         ANE=PNS(1)
C         ANDX=0.5D0*(RN(NR+1,2)+RN(NR,2))
C         ANT =0.5D0*(RN(NR+1,3)+RN(NR,3))
C         ANA =0.5D0*(RN(NR+1,4)+RN(NR,4))
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
         COEF = 12.D0*PI*SQRT(PI)*AEPS0**2
     &         /(ANE*1.D20*ZEFFL*AEE**4*15.D0)
         TAUE = COEF*SQRT(AME)*(TEL*RKEV)**1.5D0/SQRT(2.D0)
         TAUD = COEF*SQRT(AMD)*(TDL*RKEV)**1.5D0/PZ(2)**2
         TAUT = COEF*SQRT(AMT)*(TTL*RKEV)**1.5D0/PZ(3)**2
         TAUA = COEF*SQRT(AMA)*(TAL*RKEV)**1.5D0/PZ(4)**2
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
         DRL=AR1RHO(NR)/DR
         DTE=2.D0*(PTS(1)-RT(NR,1))*DRL
         DTD=2.D0*(PTS(2)-RT(NR,2))*DRL
         DTT=2.D0*(PTS(3)-RT(NR,3))*DRL
         DTA=2.D0*(PTS(4)-RT(NR,4))*DRL
         DPE=2.D0*(PNSS(1)*PTS(1)-RN(NR,1)*RT(NR,1))*DRL
         DPD=2.D0*(PNSS(2)*PTS(2)-RN(NR,2)*RT(NR,2))*DRL
         DPT=2.D0*(PNSS(3)*PTS(3)-RN(NR,3)*RT(NR,3))*DRL
         DPA=2.D0*(PNSS(4)*PTS(4)-RN(NR,4)*RT(NR,4))*DRL
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
         AJBSL(NR)=-PBSCD*SQRT(EPS)*ANE*1.D20*TEL*RKEV/BPL
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
C           BOOTSTRAP CURRENT (KEPT ON 96/08/13)
C
C     ***********************************************************
C
      SUBROUTINE TRAJBSO
C
      INCLUDE 'trcomm.h'
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
      DO 100 NR=1,NRMAX
C
         EPS=RA*RM(NR)/RR
         EPSS=SQRT(EPS)**3
         ANE=RN(NR,1)
         TEL=ABS(RT(NR,1))
         TDL=ABS(RT(NR,2))
         TTL=ABS(RT(NR,3))
         TAL=ABS(RT(NR,4))
         PEL=RN(NR,1)*RT(NR,1)
         PDL=RN(NR,2)*RT(NR,2)
         PTL=RN(NR,3)*RT(NR,3)
         PAL=RN(NR,4)*RT(NR,4)
         ZEFFL=ZEFF(NR)
C
         COEF = 12.D0*PI*SQRT(PI)*AEPS0**2
     &         /(ANE*1.D20*ZEFFL*AEE**4*15.D0)
         TAUE = COEF*SQRT(AME)*(TEL*RKEV)**1.5D0/SQRT(2.D0)
         TAUD = COEF*SQRT(AMD)*(TDL*RKEV)**1.5D0/PZ(2)**2
         TAUT = COEF*SQRT(AMT)*(TTL*RKEV)**1.5D0/PZ(3)**2
         TAUA = COEF*SQRT(AMA)*(TAL*RKEV)**1.5D0/PZ(4)**2
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
         RK3D=((1.17D0-0.35D0*SQRT(RNUD))
     &        /(1.D0+0.7D0*SQRT(RNUD))
     &         -2.1D0*(RNUD*EPSS)**2)/(1.D0+(RNUD*EPSS)**2)
         RK3T=((1.17D0-0.35D0*SQRT(RNUT))
     &        /(1.D0+0.7D0*SQRT(RNUT))
     &         -2.1D0*(RNUT*EPSS)**2)/(1.D0+(RNUT*EPSS)**2)
         RK3A=((1.17D0-0.35D0*SQRT(RNUA))
     &        /(1.D0+0.7D0*SQRT(RNUA))
     &         -2.1D0*(RNUA*EPSS)**2)/(1.D0+(RNUA*EPSS)**2)
C
         IF(NR.EQ.1) THEN
            DRL=AR1RHO(NR)/(2.D0*DR)
            DTE=(RT(2,1)-RT(1,1))*DRL
            DTD=(RT(2,2)-RT(1,2))*DRL
            DTT=(RT(2,3)-RT(1,3))*DRL
            DTA=(RT(2,4)-RT(1,4))*DRL
            DPE=(RN(2,1)*RT(2,1)-RN(1,1)*RT(1,1))*DRL
            DPD=(RN(2,2)*RT(2,2)-RN(1,2)*RT(1,2))*DRL
            DPT=(RN(2,3)*RT(2,3)-RN(1,3)*RT(1,3))*DRL
            DPA=(RN(2,4)*RT(2,4)-RN(1,4)*RT(1,4))*DRL
            BPL=0.5D0*BP(NR)
         ELSEIF(NR.EQ.NRMAX) THEN
            DRL=AR1RHO(NR)/DR
            DTE=(PTS(1) -0.5D0*(RT(NR,1)+RT(NR-1,1)))*DRL
            DTD=(PTS(2) -0.5D0*(RT(NR,2)+RT(NR-1,2)))*DRL
            DTT=(PTS(3) -0.5D0*(RT(NR,3)+RT(NR-1,3)))*DRL
            DTA=(PTS(4) -0.5D0*(RT(NR,4)+RT(NR-1,4)))*DRL
            DPE=(PNSS(1)*PTS(1) 
     &           -0.5D0*(RN(NR,1)*RT(NR,1)+RN(NR-1,1)*RT(NR-1,1)))*DRL
            DPD=(PNSS(2)*PTS(2) 
     &           -0.5D0*(RN(NR,2)*RT(NR,2)+RN(NR-1,2)*RT(NR-1,2)))*DRL
            DPT=(PNSS(3)*PTS(3) 
     &           -0.5D0*(RN(NR,3)*RT(NR,3)+RN(NR-1,3)*RT(NR-1,3)))*DRL
            DPA=(PNSS(4)*PTS(4) 
     &           -0.5D0*(RN(NR,4)*RT(NR,4)+RN(NR-1,4)*RT(NR-1,4)))*DRL
            BPL=0.5D0*(BP(NR-1)+BP(NR))
         ELSE
            DRL=AR1RHO(NR)/(2.D0*DR)
            DTE=(RT(NR+1,1)-RT(NR-1,1))*DRL
            DTD=(RT(NR+1,2)-RT(NR-1,2))*DRL
            DTT=(RT(NR+1,3)-RT(NR-1,3))*DRL
            DTA=(RT(NR+1,4)-RT(NR-1,4))*DRL
            DPE=(RN(NR+1,1)*RT(NR+1,1)-RN(NR-1,1)*RT(NR-1,1))*DRL
            DPD=(RN(NR+1,2)*RT(NR+1,2)-RN(NR-1,2)*RT(NR-1,2))*DRL
            DPT=(RN(NR+1,3)*RT(NR+1,3)-RN(NR-1,3)*RT(NR-1,3))*DRL
            DPA=(RN(NR+1,4)*RT(NR+1,4)-RN(NR-1,4)*RT(NR-1,4))*DRL
            BPL=0.5D0*(BP(NR-1)+BP(NR))
         ENDIF
C
         FACT=1.D0/(1.D0+(RNUE*EPS)**2)
         A=                   DPE/PEL-    2.5D0*DTE/TEL
     &    +(PDL/PEL)*(DPD/PDL-RK3D*FACT*DTD/TDL)
     &    +(PTL/PEL)*(DPT/PTL-RK3T*FACT*DTT/TTL)
     &    +(PAL/PEL)*(DPA/PAL-RK3A*FACT*DTA/TAL)
         AJBS(NR)=-PBSCD*SQRT(EPS)*RN(NR,1)*1.D20*RT(NR,1)*RKEV/BPL
     &            *(RK13E*A+RK23E*DTE/TEL)
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
  100 CONTINUE
C
      RETURN
      END
C
C     ***********************************************************
C
C           BOOTSTRAP CURRENT (OLD VERSION)
C
C     ***********************************************************
C
      SUBROUTINE OLDTRAJBS
C
      INCLUDE 'trcomm.h'
C
      IF(PBSCD.LE.0.D0) RETURN
C
      AMD=PA(2)*AMM
      AMT=PA(3)*AMM
      AMA=PA(4)*AMM
C
      DO 100 NR=1,NRMAX
C
         ANE=RN(NR,1)
         TEL=ABS(RT(NR,1))
         TDL=ABS(RT(NR,2))
         TTL=ABS(RT(NR,3))
         TAL=ABS(RT(NR,4))
         ZEFFL=ZEFF(NR)
C
         COEF = 6.D0*PI*SQRT(2.D0*PI)*AEPS0**2/(1.D20*AEE**4*15.D0)
         TAUE = COEF*SQRT(AME)*(TEL*RKEV)**1.5D0/ANE
         TAUD = COEF*SQRT(AMD)*(TDL*RKEV)**1.5D0/ANE
         TAUT = COEF*SQRT(AMT)*(TTL*RKEV)**1.5D0/ANE
         TAUA = COEF*SQRT(AMA)*(TAL*RKEV)**1.5D0/ANE
         TAUBE = TAUE*2.0D0/(1.D0+ZEFFL)
         TAUBD = TAUD/ZEFFL
         TAUBT = TAUT/ZEFFL
         TAUBA = TAUA/ZEFFL
C
C         COEF=6.D0*AEPS0**2*(PI*RKEV)**1.5D0
C     &       /(ANE*1.D20*ZEFF(NR)*AEE**4*15.D0)
C         TAUBE=     COEF*TEL**1.5D0*SQRT(AME)
C         TAUBD=2.D0*COEF*TDL**1.5D0*SQRT(AMD)/PZ(2)**2
C         TAUBT=2.D0*COEF*TTL**1.5D0*SQRT(AMT)/PZ(3)**2
C         TAUBA=2.D0*COEF*TAl**1.5D0*SQRT(AMA)/PZ(4)**2
C
         VTE=SQRT(TEL*RKEV/AME)
         VTD=SQRT(TDL*RKEV/AMD)
         VTT=SQRT(TTL*RKEV/AMT)
         VTA=SQRT(TAL*RKEV/AMA)
C
         EPS=EPSRHO(NR)
         RNUES=QP(NR)*RR/(TAUBE*VTE*EPS**1.5D0)
         RNUDS=QP(NR)*RR/(TAUBD*VTD*EPS**1.5D0)
         RNUTS=QP(NR)*RR/(TAUBT*VTT*EPS**1.5D0)
         RNUAS=QP(NR)*RR/(TAUBA*VTA*EPS**1.5D0)
C
         RK13E=2.44D0/(1.D0+0.85D0*RNUES)
         RK13D=2.44D0/(1.D0+0.85D0*RNUDS)
         RK13T=2.44D0/(1.D0+0.85D0*RNUTS)
         RK13A=2.44D0/(1.D0+0.85D0*RNUAS)
         RK23E=4.35D0/(1.D0+0.40D0*RNUES)
         RK23D=(1.33D0+3.D0*RNUDS)/(1.D0+RNUDS)
         RK23T=(1.33D0+3.D0*RNUTS)/(1.D0+RNUTS)
         RK23A=(1.33D0+3.D0*RNUAS)/(1.D0+RNUAS)
C
         IF(NR.EQ.1) THEN
            DRL=AR1RHO(NR)/(2.D0*DR)
            DN =(RN(2,1)-RN(1,1))*DRL
            DTE=(RT(2,1)-RT(1,1))*DRL
            DTD=(RT(2,2)-RT(1,2))*DRL
            DTT=(RT(2,3)-RT(1,3))*DRL
            DTA=(RT(2,4)-RT(1,4))*DRL
            BPL=0.5D0*BP(NR)
         ELSEIF(NR.EQ.NRMAX) THEN
            DRL=AR1RHO(NR)/DR
            DN =(PNSS(1)-0.5D0*(RN(NR,1)+RN(NR-1,1)))*DRL
            DTE=(PTS(1) -0.5D0*(RT(NR,1)+RT(NR-1,1)))*DRL
            DTD=(PTS(2) -0.5D0*(RT(NR,2)+RT(NR-1,2)))*DRL
            DTT=(PTS(3) -0.5D0*(RT(NR,3)+RT(NR-1,3)))*DRL
            DTA=(PTS(4) -0.5D0*(RT(NR,4)+RT(NR-1,4)))*DRL
            BPL=0.5D0*(BP(NR-1)+BP(NR))
         ELSE
            DRL=AR1RHO(NR)/(2.D0*DR)
            DN =(RN(NR+1,1)-RN(NR-1,1))*DRL
            DTE=(RT(NR+1,1)-RT(NR-1,1))*DRL
            DTD=(RT(NR+1,2)-RT(NR-1,2))*DRL
            DTT=(RT(NR+1,3)-RT(NR-1,3))*DRL
            DTA=(RT(NR+1,4)-RT(NR-1,4))*DRL
            BPL=0.5D0*(BP(NR-1)+BP(NR))
         ENDIF
         AJBS(NR)=-SQRT(RM(NR))
     &           *(RK13E*(RN(NR,1)*RT(NR,1)
     &                   +RN(NR,2)*RT(NR,2)
     &                   +RN(NR,3)*RT(NR,3)
     &                   +RN(NR,4)*RT(NR,4)
     &                   +RW(NR,1)
     &                   +RW(NR,2)         )*DN/ANE
     &            +(RK23E-1.5D0*RK13E)*RN(NR,1)*DTE
     &            +RK13D*(RK23D-1.5D0)*RN(NR,2)*DTD
     &            +RK13T*(RK23T-1.5D0)*RN(NR,3)*DTT
     &            +RK13A*(RK23A-1.5D0)*RN(NR,4)*DTA)
     &           *PBSCD*1.D20*RKEV/(SQRT(RR)*BPL)
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
  100 CONTINUE
C
      RETURN
      END
C     ***********************************************************
C
C           CALCULATE AJ, AJOH, POH, EZOH
C
C     ***********************************************************
C
      SUBROUTINE TRAJOH
C
      INCLUDE 'trcomm.h'
C
      IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
c$$$         NR=1
c$$$            AJ(NR) = (RG(NR)*RA*BP(NR)                  )
c$$$     &              /(RM(NR)*RA**2*DR)*FKAP/(RKAP*AMYU0)
c$$$CCC            AJ(NR) = (RG(NR)*RA*BP(NR)                  )
c$$$CCC     &              /(RM(NR)*RA**2*DR)/(RKAP*AMYU0)
c$$$         DO NR=2,NRMAX
c$$$            AJ(NR) = (RG(NR)*RA*BP(NR)-RG(NR-1)*RA*BP(NR-1))
c$$$     &              /(RM(NR)*RA**2*DR)*FKAP/(RKAP*AMYU0)
c$$$CCC            AJ(NR) = (RG(NR)*RA*BP(NR)-RG(NR-1)*RA*BP(NR-1))
c$$$CCC     &              /(RM(NR)*RA**2*DR)/(RKAP*AMYU0)
c$$$         ENDDO
         NR=1
            FACTOR0=TTRHO(NR)/(ARRHO(NR)*AMYU0*DVRHO(NR))
            FACTOR2=DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  )
            FACTOR3=DVRHO(NR+1)*ABRHO(NR+1)/TTRHO(NR+1)
            FACTORP=0.5D0*(FACTOR2+FACTOR3)
            AJ(NR)= FACTOR0*FACTORP*BP(NR)/DR*RA
     &             *(FKAP/RKAP)
         DO NR=2,NRMAX-1
            FACTOR0=TTRHO(NR)/(ARRHO(NR)*AMYU0*DVRHO(NR))
            FACTOR1=DVRHO(NR-1)*ABRHO(NR-1)/TTRHO(NR-1)
            FACTOR2=DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  )
            FACTOR3=DVRHO(NR+1)*ABRHO(NR+1)/TTRHO(NR+1)
            FACTORM=0.5D0*(FACTOR1+FACTOR2)
            FACTORP=0.5D0*(FACTOR2+FACTOR3)
            AJ(NR)= FACTOR0*(FACTORP*BP(NR)-FACTORM*BP(NR-1))/DR*RA
     &             *(FKAP/RKAP)
C            if(nr.ge.nrmax-4)
C     &     write(6,'(I3,1P4E19.12)') NR,FACTORP,BP(NR),FACTORM,BP(NR-1)
         ENDDO
         NR=NRMAX
            FACTOR0=TTRHO(NR)/(ARRHO(NR)*AMYU0*DVRHO(NR))
            FACTOR1=DVRHO(NR-1)*ABRHO(NR-1)/TTRHO(NR-1)
            FACTOR2=DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  )
            FACTORM=0.5D0*(FACTOR1+FACTOR2)
            FACTORP=0.5D0*(3.D0*FACTOR2-FACTOR1)
            AJ(NR)= FACTOR0*(FACTORP*BP(NR)-FACTORM*BP(NR-1))/DR*RA
     &             *(FKAP/RKAP)
C            write(6,'(I3,1P4E19.12)') NR,FACTORP,BP(NR),FACTORM,BP(NR-1)
C            write(6,*) "*****"
      ELSE
         NR=1
            FACTOR0=TTRHO(NR)/(ARRHO(NR)*AMYU0*DVRHO(NR))
            FACTOR2=DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  )
            FACTOR3=DVRHO(NR+1)*ABRHO(NR+1)/TTRHO(NR+1)
            FACTORP=0.5D0*(FACTOR2+FACTOR3)
            AJ(NR)= FACTOR0*FACTORP*BP(NR)/DR*RA
         DO NR=2,NRMAX-1
            FACTOR0=TTRHO(NR)/(ARRHO(NR)*AMYU0*DVRHO(NR))
            FACTOR1=DVRHO(NR-1)*ABRHO(NR-1)/TTRHO(NR-1)
            FACTOR2=DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  )
            FACTOR3=DVRHO(NR+1)*ABRHO(NR+1)/TTRHO(NR+1)
            FACTORM=0.5D0*(FACTOR1+FACTOR2)
            FACTORP=0.5D0*(FACTOR2+FACTOR3)
            AJ(NR)= FACTOR0*(FACTORP*BP(NR)-FACTORM*BP(NR-1))/DR*RA
         ENDDO
         NR=NRMAX
            FACTOR0=TTRHO(NR)/(ARRHO(NR)*AMYU0*DVRHO(NR))
            FACTOR1=DVRHO(NR-1)*ABRHO(NR-1)/TTRHO(NR-1)
            FACTOR2=DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  )
            FACTORM=0.5D0*(FACTOR1+FACTOR2)
            FACTORP=0.5D0*(3.D0*FACTOR2-FACTOR1)
            AJ(NR)= FACTOR0*(FACTORP*BP(NR)-FACTORM*BP(NR-1))/DR*RA
      ENDIF
C
      DO NR=1,NRMAX
         AJOH(NR) = AJ(NR)-(AJNB(NR  )+AJRF(NR  )+AJBS(NR ))
         EZOH(NR) = ETA(NR)*AJOH(NR)
         POH(NR)  = EZOH(NR)*AJOH(NR)
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
      INCLUDE 'trcomm.h'
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
      DO 10 NR=1,NRMAX
         SUM = SUM+(1.D0/ABS(QP(NR))-1.D0)*DSRHO(NR)*DR
         IF(SUM.LT.0.D0) GOTO 1000
   10 CONTINUE
      NR=NRMAX
C
 1000 IZEROX=NR
C
      DO 20 NR=1,NRMAX
         IF(QP(NR).GE.1.D0) GOTO 2000
   20 CONTINUE
      NR=NRMAX
C
 2000 IONE=NR
C
      SUM = 0.D0
      DO 30 NR=IONE,IZEROX
         SUM = SUM+(1.D0/ABS(QP(NR))-1.D0)*DSRHO(NR)*DR
   30 CONTINUE
C
      DO 40 NR=1,IZEROX
         QONE(NR) = 1.D0+SUM*4.D0*RG(NR)**2/RG(IZEROX)**4
   40 CONTINUE
C
      IF(LT.EQ.1) THEN
         DO 70 NS=1,NSM
            SUM1 = 0.D0
            SUM2 = 0.D0
            DO 50 NR=1,IZEROX
               SUM1 = SUM1+RN(NR,NS)          *DVRHO(NR)
               SUM2 = SUM2+RN(NR,NS)*RT(NR,NS)*DVRHO(NR)
   50       CONTINUE
            RTN = SUM2/SUM1
C
            DO 60 NR=1,IZEROX
               RT(NR,NS) = RTN
   60       CONTINUE
   70    CONTINUE
      ENDIF
C
      IF(LN.EQ.1) THEN
         DO 100 NS=1,NSM
            SUM1 = 0.D0
            SUM2 = 0.D0
            DO 80 NR=1,IZEROX
               SUM1 = SUM1+          DVRHO(NR)
               SUM2 = SUM2+RN(NR,NS)*DVRHO(NR)
   80       CONTINUE
            RNN = SUM2/SUM1
            DO 90 NR=1,IZEROX
               RN(NR,NS) = RNN
   90       CONTINUE
  100    CONTINUE
      ENDIF
C
      IF(LQ.EQ.1) THEN
         IF(MODELG.EQ.3) THEN
            DO NR=1,IZEROX
               QP(NR) = 1.D0/QONE(NR)
               BP(NR) = BPRHO(NR)*QRHO(NR)/QP(NR)
            ENDDO
         ELSE
            DO NR=1,IZEROX
               QP(NR) = 1.D0/QONE(NR)
               BP(NR)  = FKAP*RA*RG(NR)*BB/(RR*QP(NR))
CCC               BP(NR)  = RA*RG(NR)*BB/(RR*QP(NR))
C               write(6,*) NR,BP(NR)
            ENDDO
         ENDIF
      ENDIF
C
         WRITE(6,602) RM(IONE),RM(IZEROX),RTN,RNN
  602    FORMAT(' ',' R-ONE,R-ZERO,RTN,RNN = ',4F8.3)
      RETURN
      END
