C     $Id$
C
C     ***********************************************************
C
C           PRODUCTION OF ALPHA PARTICLES
C
C     ***********************************************************
C
      SUBROUTINE TRALPH
C
      INCLUDE 'trcomm.h'
C
      IF(MDLNF.EQ.0) THEN
         DO 5 NR=1,NRMAX
            TAUF(NR)=1.D0
    5    CONTINUE
         RETURN
      ENDIF
C
      AMD=PA(2)*AMM
      AMT=PA(3)*AMM
      AMA=PA(4)*AMM
      VF =SQRT(2.D0*3.5D3 *RKEV/AMA)
C
      DO 10 NR=1,NRMAX
         ANE= RN(NR,1)
         TE = RT(NR,1)
         TD = RT(NR,2)
         TT = RT(NR,3)
         SS = SIGMAM(TD,TT)
         ZEFFM = (PZ(2)  *PZ(2)  *RN(NR,2)/PA(2)
     &           +PZ(3)  *PZ(3)  *RN(NR,3)/PA(3)
     &           +PZ(4)  *PZ(4)  *RN(NR,4)/PA(4)
     &           +PZC(NR) *PZC(NR) *ANC(NR) /12.D0
     &           +PZFE(NR)*PZFE(NR)*ANFE(NR)/52.D0)/ANE
         EC  = 14.8D0*TE*PA(2)*ZEFFM**(2.D0/3.D0)
         TAUS= 0.2D0*PA(2)*ABS(TE)**1.5D0/(PZ(2)**2*ANE*15.D0)
         PTNT= PBIN(NR)*TAUS/(RN(NR,2)*1.D20*PNBENG*RKEV)
         SSB = SIGMAB(PNBENG,EC,TT,PTNT)
         SNF(NR) = (SS+SSB)*RN(NR,2)*RN(NR,3)*1.D20
         PNF(NR) = SNF(NR)*3.5D3*RKEV*1.D20
         IF(MDLNF.LE.1) SNF(NR) = 0.D0
   10 CONTINUE
C
      DO 20 NR=1,NRMAX
         ANE= RN(NR,1)
         TE = RT(NR,1)
         WF = 1.5D0*RW(NR,2)
         P1   = 3.D0*SQRT(0.5D0*PI)*AME/ANE
     &         *(ABS(TE)*RKEV/AME)**1.5D0
         VCD3 = P1*RN(NR,2)*PZ(2)**2/AMD
         VCT3 = P1*RN(NR,3)*PZ(3)**2/AMT
         VCA3 = P1*RN(NR,4)*PZ(4)**2/AMA
         VC3  = VCD3+VCT3+VCA3
         VCR  = VC3**(1.D0/3.D0)
         HYF=HY(VF/VCR)
         TAUS = 0.2D0*PA(4)*ABS(TE)**1.5D0/(PZ(4)**2*ANE*15.D0)
         TAUF(NR)= 0.5D0*TAUS*(1.D0-HYF)
         RNF(NR,2)= 2.D0*LOG(1.D0+(VF/VCR)**3)*WF
     &             /(3.D0*(1.D0-HYF)*3.5D3)
         IF(RNF(NR,2).GT.0.D0) THEN
            RTF(NR,2)= WF/(1.5D0*RNF(NR,2))
         ELSE
            RTF(NR,2)= 0.D0
         ENDIF
         PFIN(NR) = WF*RKEV*1.D20/TAUF(NR)
         PFCL(NR,1)=    (1.D0-HYF)*PFIN(NR)
         PFCL(NR,2)=(VCD3/VC3)*HYF*PFIN(NR)
         PFCL(NR,3)=(VCT3/VC3)*HYF*PFIN(NR)
         PFCL(NR,4)=(VCA3/VC3)*HYF*PFIN(NR)
   20 CONTINUE
C
      RETURN
      END
C
C     ***********************************************************
C
C           REACTION CROSS SECTION (MAXELLIAN)
C
C     ***********************************************************
C
      FUNCTION SIGMAM(TD,TT)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      TI = (3.D0*ABS(TD)+2.D0*ABS(TT))/5.D0
      H  = TI/37.D0
     &     +5.45D0/(3.D0+TI*(1.D0+(TI/37.5D0)**2.8D0))
      ARG= -20.D0/TI**(1.D0/3.D0)
      IF(ARG.GE.-100.D0)  THEN
         SIGMAM = 3.7D-18*TI**(-2.D0/3.D0)*EXP(ARG)/H
      ELSE
         SIGMAM = 0.D0
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
C      REACTION RATE : TAIL
C
C     ***********************************************************
C
      FUNCTION SIGMAB(EB,EC,TI,PTNT)
C
C      APPROXIMATE FORMULA OF FUSION REACTION RATE
C         FOR SLOWING DOWN ION DISTRIBUTION
C      REF. TAKIZUKA AND YAMAGIWA, JAERI-M 87-066
C
C      EB   : BEAM ENERGY (KEV)
C      EC   : CRITICAL ENERGY (KEV)
C      TI   : TRITIUM TEMPERATURE (KEV)
C      PTNT : PB * TAUS / (ND * EB)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      XB=SQRT(EB/127.D0)
      XC=SQRT(EC/127.D0)
C
      AG1= 1.06D0-0.058D0*SQRT(ABS(TI))
      AG2= 1.06D0-0.058D0*SQRT(ABS(TI))
      AG3= 0.33D0
      AL1= 1.D0/(0.40D0+0.032D0*TI)
      AL2=-1.D0/(0.91D0+0.016D0*SQRT(ABS(TI)))
      AL3=-0.11D0
      X1=0.97D0-AG1/AL1
      X2=0.97D0
      X3=0.97D0+(AG2-AG3)/(AL3-AL2)
      X4=0.97D0+3.D0
C
      IF(XB.LT.X1) THEN
         SA=0.D0
      ELSE
         SA=-SIGMBS(X1,AG1,AL1,XC)
         IF(XB.LT.X2) THEN
            SA=SA+SIGMBS(XB,AG1,AL1,XC)
         ELSE
            SA=SA+SIGMBS(X2,AG1,AL1,XC)-SIGMBS(X2,AG2,AL2,XC)
            IF(XB.LT.X3) THEN
               SA=SA+SIGMBS(XB,AG2,AL2,XC)
            ELSE
               SA=SA+SIGMBS(X3,AG2,AL2,XC)-SIGMBS(X3,AG3,AL3,XC)
               IF(XB.LT.X4) THEN
                  SA=SA+SIGMBS(XB,AG3,AL3,XC)
               ELSE
                  SA=SA+SIGMBS(X4,AG3,AL3,XC)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      SIGMAB=PTNT*1.67D-21*SA
      RETURN
      END
C
      FUNCTION SIGMBS(XX,RGG,RGL,XC)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      X=XX/XC
      SIGMBS=((RGG-0.97D0*RGL)/3.D0+XC*RGL/6.D0)*LOG(X*X*X+1.D0)
     &      +XC*RGL*(X-LOG(X+1.D0)/2.D0
     &               -ATAN((2.D0*X-1.D0)/SQRT(3.D0))/SQRT(3.D0))
      RETURN
      END
