C     $Id$
C
C     ***********************************************************
C
C           CALCULATE GLOBAL QUANTITIES
C
C     ***********************************************************
C
      SUBROUTINE TXGLOB
C
      INCLUDE 'txcomm.inc'
C
      DIMENSION BP(NRM),BETA(NRM),BETAP(NRM)
      DIMENSION BETAL(NRM),BETAPL(NRM),BETAQ(NRM)
C
C     Line Averaged Density and Temperature
C     Core Density and Temperature
C     Electoron and ion Work Quantities
C
      RKAP=1.D0
      FKAP=1.D0
C
      CALL TXSUMD(PNeHI,RHI,NRMAX,RNSUM)
      CALL TXSUMT(PTeHI,PNeHI,RHI,NRMAX,RTSUM)
      ANSAV(1) = RNSUM*2.D0*PI*DR/(PI*RA*RA)
      ANS0(1)  = PNeI(0)
      IF(RNSUM.GT.0.D0) THEN
         TSAV(1) = RTSUM/RNSUM
      ELSE
         TSAV(1)=0.D0
      ENDIF
      TS0(1) = PTeI(0)
      WST(1) = RTSUM*1.5D0*2.D0*PI*RR*2.D0*PI*DR*RKEV*1.D14
C
      CALL TXSUMD(PNiHI,RHI,NRMAX,RNSUM)
      CALL TXSUMT(PTiHI,PNeHI,RHI,NRMAX,RTSUM)
      ANSAV(2) = RNSUM*2.D0*PI*DR/(PI*RA*RA)
      ANS0(2)  = PNiI(0)
      IF(RNSUM.GT.0.D0) THEN
         TSAV(2) = RTSUM/RNSUM
      ELSE
         TSAV(2)=0.D0
      ENDIF
      TS0(2) = PTiI(0)
      WST(2) = RTSUM*1.5D0*2.D0*PI*RR*2.D0*PI*DR*RKEV*1.D14
C
         CALL TXSUMD(PNbHI,RHI,NRMAX,ANFSUM)
         CALL TXSUMD(SNB,RHI,NRMAX,RWSUM)
         WFT(1) = 0.5D0*AMi*RWSUM**2.D0*1.5D0
     &           *2.D0*PI*RR*2.D0*PI*DR*RKAP*RKEV*1.D14
         ANFAV(1) = ANFSUM*2.D0*PI*DR/(PI*RA*RA)
         ANF0(1)  = PNbI(0)
         IF(ANFSUM.GT.0.D0) THEN
            TFAV(1)  = RWSUM/ANFSUM
         ELSE
            TFAV(1)  = 0.D0
         ENDIF
         IF(ANF0(1).GT.0.D0) THEN
            TF0(1)  = (9.D0*SNB(1)-SNB(2))/8.D0
     &                 /ANF0(1)
         ELSE
            TF0(1)  = 0.D0
         ENDIF
C   20 CONTINUE
C
C     Input power
C
      CALL TXSUMD(POH,RHI,NRMAX,POHSUM)
      CALL TXSUMD(PNB,RHI,NRMAX,PNBSUM)
C      CALL TXSUMD(PNF,RHI,NRMAX,PNFSUM)
      PNFSUM=0.D0
      POHT = POHSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
      PNBT = PNBSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
      PNFT = PNFSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
C
      CALL TXSUMD(PRFe,RHI,NRMAX,PRFeSUM)
      CALL TXSUMD(PRFi,RHI,NRMAX,PRFiSUM)
      PRFeTOT=PRFeSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
      PRFiTOT=PRFiSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
      PRFT=PRFeTOT+PRFiTOT
C
C      DO 30 NS=1,NSM
C        CALL TXSUMD(PRF,RHI,NRMAX,PRFSUM)
C        PRFT = PRFSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
C   30 CONTINUE
C
C      CALL TXSUMD(PBIN,RHI,NRMAX,PBSUM)
C      PBINT = PBSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
CC      DO 40 NS=1,NSM
C        CALL TXSUMD(PBCL,RHI,NRMAX,PBSUM)
C        PBCLT = PBSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
C   40 CONTINUE
C
C      CALL TXSUMD(PFIN,RHI,NRMAX,PFSUM)
C      PFINT = PFSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
C      DO 50 NS=1,NSM
C        CALL TXSUMD(PFCL(1,NS),RHI,NRMAX,PFSUM)
C        PFCLT(NS) = PFSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
C   50 CONTINUE
C
C      CALL TXSUMD(PRL,RHI,NRMAX,PRLSUM)
C
C     Output power
C
      PIESUM=0.D0
      SIESUM=0.D0
      PCXSUM=0.D0
C
      DO I=1,NRMAX
         EION  = 13.64D0
         PIE=PNeHI(I)*rNuION(I)*1.D20*EION*AEE
         SIE=PNeHI(I)*rNuION(I)*1.D20
         PIESUM=PIESUM+PIE*RHI(I)
         SIESUM=SIESUM+SIE*RHI(I)
C
         TNU=0.D0
         PCX=(-1.5D0*PNeHI(I)*rNuION(I)*TNU/1.D20
     &             +1.5D0*PNiHI(I)*rNuiCX(I)
     &             *(PTiHI(I)-TNU))*RKEV*1.D20
         PCXSUM=PCXSUM+PCX*RHI(I)
      ENDDO
C
C      CALL TXSUMD(PCX,RHI,NRMAX,PCXSUM)
C      CALL TXSUMD(PIE,RHI,NRMAX,PIESUM)
C      PRLT = PRLSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
      PCXT = PCXSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
      PIET = PIESUM*2.D0*PI*RR*2.D0*PI*DR*RKAP/1.D6
C
C     Current
C
      CALL TXSUMD(AJ  ,RHI,NRMAX,AJTSUM)
      CALL TXSUMD(AJOH,RHI,NRMAX,AOHSUM)
      CALL TXSUMD(AJNB,RHI,NRMAX,ANBSUM)
C
      AJT   = AJTSUM*            2.D0*PI*DR*RKAP/1.D6
      AJOHT = AOHSUM*            2.D0*PI*DR*RKAP/1.D6
      AJNBT = ANBSUM*            2.D0*PI*DR*RKAP/1.D6
C
C      DRH=0.5D0*DR
C      DO NS=1,NSM
C         VNP=AV(NRMAX,NS)
C         DNP=AD(NRMAX,NS)
C
C         VTP=(AVK(NRMAX,NS)/1.5D0+AV(NRMAX,NS))
C         DTP= AK(NRMAX,NS)/(1.5D0)
C
C         VXP=0.D0
C         DXP=(1.5D0*AD(NRMAX,NS)-AK(NRMAX,NS))*PTS(NS)
C
C         SLT(NS) =((     DNP/DRH)*RN(NRMAX,NS)
C     &            +( VNP-DNP/DRH)*PNSS(NS))
C     &            *2.D0*PI*RR*2.D0*PI*RA*FKAP
C
C         PLT(NS) =((     DXP/DRH)*RN(NRMAX,NS)
C     &            +(     DTP/DRH)*RN(NRMAX,NS)*RT(NRMAX,NS)*1.5D0
C     &            +( VXP-DXP/DRH)*PNSS(NS)
C     &            +( VTP-DTP/DRH)*PNSS(NS)*PTS(NS)*1.5D0)
C     &            *2.D0*PI*RR*2.D0*PI*RA*FKAP*RKEV*1.D14
C      ENDDO
C
C      CALL TXSUMD(SIE,RHI,NRMAX,SIESUM)
C      CALL TXSUMD(SNF,RHI,NRMAX,SNFSUM)
      CALL TXSUMD(SNB,RHI,NRMAX,SNBSUM)
      SIET = SIESUM*2.D0*PI*RR*2.D0*PI*DR*RKAP
C      SNFT = SNFSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP
      SNBT = SNBSUM*2.D0*PI*RR*2.D0*PI*DR*RKAP
C
C      DO NS=1,NSM
C         CALL TXSUMD(SPE(1,NS),RHI,NRMAX,SPESUM)
C         SPET(NS) = SPESUM*2.D0*PI*RR*2.D0*PI*DR
C      ENDDO
C
      WBULKT=0.D0
      DO NS=1,NSM
         WBULKT=WBULKT+WST(NS)
      ENDDO
C
      WTAILT=0.D0
C      DO NF=1,NFM
         NF=1
         WTAILT=WTAILT+WFT(NF)
C      ENDDO
C
      WPT =WBULKT+WTAILT
      PINT=PNBT+PRFT+POHT
      POUT=PCXT+PIET
      SINT=SIET+SNBT
C      SOUT=SLST
C
      IF(ABS(TIME-TPRE).LE.1.D-70) THEN
         WPDOT=0.D0
      ELSE
         WPDOT=(WPT-WPPRE)/(TIME-TPRE)
      ENDIF
      WPPRE=WPT
      TPRE=TIME
C
      IF(PINT.LE.0.D0) THEN
         TAUE1=0.D0
         TAUE2=0.D0
      ELSE
         TAUE1=WPT/PINT
         TAUE2=WPT/(PINT-WPDOT)
      ENDIF
C
      SUM=0.D0
      SUP=0.D0
      SUL=0.D0
      FACT=1.D0
      FACT=RKAP/FKAP**2
C
         DO I=1,NRMAX
            BP(I)=BthI(I)
         ENDDO
C
      DO I=1,NRMAX-1
         BBL  = SQRT(BthI(I)**2+BphI(I)**2)
         SUML = PNeHI(I)*PTeHI(I)*RKEV*1.D20
     &        + PNiHI(I)*PTiHI(I)*RKEV*1.D20
         SUPL = (PNeHI(I+1)*PTeHI(I+1)-PNeHI(I  )*PTeHI(I  )
     &          +PNiHI(I+1)*PTiHI(I+1)-PNiHI(I  )*PTiHI(I  ))
     &          *RKEV*1.D20/DR
C
         SUM = SUM + SUML*2.D0*PI*RHI(I)*DR
         SUP = SUP + 0.5D0*SUPL*PI*RHI(I)*RHI(I)*DR
         SUL = SUL + BP(I)**2*R(I)
         BETA(I)   = 2.D0*SUM *rMU0/(PI*(R(I)*BBL)**2)
         BETAL(I)  = 2.D0*SUML*rMU0/(           BBL **2)
         BETAP(I)  = 2.D0*SUM *rMU0/(PI*(R(I)*BP(NRMAX))**2)*FACT
         BETAPL(I) = 2.D0*SUML*rMU0/(           BP(NRMAX)**2)*FACT
         BETAQ(I)  =-2.D0*SUP *rMU0/(PI*(R(I)*BP(I))**2)*FACT
         SUP = SUP + 0.5D0*SUPL*PI*RHI(I+1)*RHI(I+1)*DR
      ENDDO
C
      I=NRMAX
         BBL  = SQRT(BthI(I)**2+BphI(I)**2)
         SUML =(PNeHI(I)*PTeHI(I)
     &         +PNiHI(I)*PTiHI(I))
     &         *RHI(I)*RKEV*1.D20
         PNES = PNa * EXP(-(RB-RA) / rLn)
         SUPL = (PNES*PTea-PNeHI(I)*PTeHI(I)
     &          +PNES/PZ*PTia-PNiHI(I)*PTiHI(I))
     &          *RKEV*1.D20/DR*2.D0
C
C         DO NF=1,NFM
C            SUML = SUML +SNB(I)*RHI(I)*RKEV*1.D20
C            SUPL = SUPL +(0.D0-SNB(I))*RKEV*1.D20/DR*2.D0
C         ENDDO
C
         SUM = SUM + SUML*2.D0*PI*RHI(I)*DR
         SUP = SUP + 0.5D0*SUPL*PI*RHI(I)*RHI(I)*DR
         SUL = SUL + 0.5D0*BP(I)**2*R(I)
         BETA(I)   = 2.D0*SUM *rMU0/(PI*(R(I)*BBL)**2)
         BETAL(I)  = 2.D0*SUML*rMU0/(           BBL **2)
         BETAP(I)  = 2.D0*SUM *rMU0/(PI*(R(I)*BP(NRMAX))**2)*FACT
         BETAPL(I) = 2.D0*SUML*rMU0/(           BP(NRMAX)**2)*FACT
         BETAQ(I)  =-2.D0*SUP *rMU0/(PI*(R(I)*BP(I))**2)*FACT
C
      BETA0 =(4.D0*BETA(1)  -BETA(2)  )/3.D0
      BETAP0=(4.D0*BETAPL(1)-BETAPL(2))/3.D0
      BETAQ0=0.D0
C
      BETAPA=BETAP(NRMAX)
      BETAA =BETA(NRMAX)
C
      ALI=8.D0*PI**2*DR*SUL*FKAP**2/((rMU0*rIp*1.D6)**2)
      VLOOP = EphHI(NRMAX-1)*2*PI*RR
C
C      PAI=(PA(2)*PN(2)+PA(3)*PN(3)+PA(4)*PN(4))/(PN(2)+PN(3)+PN(4))
      PAI=PA
C
      IF(PINT.LE.0.D0) THEN
         TAUEP=0.D0
         QF=0.D0
      ELSE
      TAUEP=4.8D-2*(rIP**0.85D0)
     &            *(RR**1.2D0)
     &            *(RA**0.3D0)
     &            *(RKAP**0.5D0)
     &            *(ANSAV(1)**0.1D0)
     &            *(BB**0.2D0)
     &            *(PAI**0.5D0)
     &            *(PINT**(-0.5D0))
         QF=5.D0*PNFT/PINT
      ENDIF
C
      IF(Q(0).GE.1.D0) THEN
         RQ1=0.D0
      ELSE
         IF(Q(1).GT.1.D0) THEN
            RQ1=SQRT( (1.D0-Q(0) )/(Q(1)-Q(0)) )*DR
            GOTO 310
         ENDIF
         DO I=2,NRMAX
            IF(Q(I).GT.1.D0) THEN
               RQ1=(R(I)-R(I-1))*(1.D0-Q(I-1))/(Q(I)-Q(I-1))
     &             +R(I-1)
               GOTO 310
            ENDIF
         ENDDO
         RQ1=RA
  310    CONTINUE
      ENDIF
C
C      ZEFF0=(4.D0*ZEFF(1)-ZEFF(2))/3.D0
      ZEFF0=Zeff
C
      RETURN
      END
C
C     ********************************************************
C
C           RADIAL INTEGRATION  (DOUBLE VARIABLES)
C
C     ********************************************************
C
      SUBROUTINE TXSUMD(A,B,NMAX,SUM)
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NMAX
      REAL(8), DIMENSION(NMAX), INTENT(IN) :: A, B
      REAL(8), INTENT(OUT) :: SUM
      INTEGER :: N
C
      SUM=0.D0
      DO N=1,NMAX
         SUM=SUM+A(N)*B(N)
      ENDDO
C
      RETURN
      END
C
C     ********************************************************
C
C           RADIAL INTEGRATION  (TRIPLE VARIABLES)
C
C     ********************************************************
C
      SUBROUTINE TXSUMT(A,B,C,NMAX,SUM)
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NMAX
      REAL(8), DIMENSION(NMAX), INTENT(IN) :: A, B, C
      REAL(8), INTENT(OUT) :: SUM
      INTEGER :: N
C
      SUM=0.D0
      DO N=1,NMAX
         SUM=SUM+A(N)*B(N)*C(N)
      ENDDO
      RETURN
      END
