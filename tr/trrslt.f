C     $Id$
C
C     ***********************************************************
C
C           CALCULATE GLOBAL QUANTITIES
C
C     ***********************************************************
C
      SUBROUTINE TRGLOB
C
      INCLUDE 'trcomm.h'
C
      IF (MODELG.EQ.3) THEN
         FKAP=1.D0
         RKAP=1.D0
      ELSE
         FKAP=1.D0
      ENDIF
C
      VOL=0.D0
      DO NR=1,NRMAX
         VOL=VOL+DVRHO(NR)*DR
      ENDDO
C
      DO NS=1,NSM
         CALL TRSUMD(RN(1,NS),DVRHO,NRMAX,RNSUM)
         CALL TRSUMT(RT(1,NS),RN(1,NS),DVRHO,NRMAX,RTSUM)
         ANSAV(NS) = RNSUM*DR/VOL
         ANS0(NS) = (9.D0*RN(1,NS)-RN(2,NS))/8.D0
         IF(RNSUM.GT.0.D0) THEN
            TSAV(NS) = RTSUM/RNSUM
         ELSE
            TSAV(NS) = 0.D0
         ENDIF
         TS0(NS) = (9.D0*RT(1,NS)-RT(2,NS))/8.D0
         WST(NS) = 1.5D0*RTSUM*DR*RKEV*1.D14
      ENDDO
C
      DO NF=1,NFM
         CALL TRSUMD(RNF(1,NF),DVRHO,NRMAX,ANFSUM)
         CALL TRSUMD(RW(1,NF),DVRHO,NRMAX,RWSUM)
         WFT(NF) = 1.5D0*RWSUM*DR*RKEV*1.D14
         ANFAV(NF) = ANFSUM*DR/VOL
         ANF0(NF)  = (9.D0*RNF(1,NF)-RNF(2,NF))/8.D0
         IF(ANFSUM.GT.0.D0) THEN
            TFAV(NF)  = RWSUM/ANFSUM
         ELSE
            TFAV(NF)  = 0.D0
         ENDIF
         IF(RNF(1,NF).GT.0.D0) THEN
            TF0(NF)  = (9.D0*RW(1,NF)-RW(2,NF))/8.D0
     &                 /ANF0(NF)
         ELSE
            TF0(NF)  = 0.D0
         ENDIF
      ENDDO
C
C      DO NR=1,NRMAX
C         write(6,*) NR,POH(NR),DVRHO(NR)
C      ENDDO
      CALL TRSUMD(POH,DVRHO,NRMAX,POHSUM)
      CALL TRSUMD(PNB,DVRHO,NRMAX,PNBSUM)
      CALL TRSUMD(PNF,DVRHO,NRMAX,PNFSUM)
      POHT = POHSUM*DR/1.D6
      PNBT = PNBSUM*DR/1.D6
      PNFT = PNFSUM*DR/1.D6
      DO NS=1,NSM
         CALL TRSUMD(PEX(1,NS),DVRHO,NRMAX,PEXSUM)
         PEXT(NS) = PEXSUM*DR/1.D6
      ENDDO
C
      DO 30 NS=1,NSM
        CALL TRSUMD(PRF(1,NS),DVRHO,NRMAX,PRFSUM)
        PRFT(NS) = PRFSUM*DR/1.D6
   30 CONTINUE
C
      CALL TRSUMD(PBIN,DVRHO,NRMAX,PBSUM)
      PBINT = PBSUM*DR/1.D6
      DO 40 NS=1,NSM
        CALL TRSUMD(PBCL(1,NS),DVRHO,NRMAX,PBSUM)
        PBCLT(NS) = PBSUM*DR/1.D6
   40 CONTINUE
C
      CALL TRSUMD(PFIN,DVRHO,NRMAX,PFSUM)
      PFINT = PFSUM*DR/1.D6
      DO 50 NS=1,NSM
        CALL TRSUMD(PFCL(1,NS),DVRHO,NRMAX,PFSUM)
        PFCLT(NS) = PFSUM*DR/1.D6
   50 CONTINUE
C
      CALL TRSUMD(PRL,DVRHO,NRMAX,PRLSUM)
      CALL TRSUMD(PCX,DVRHO,NRMAX,PCXSUM)
      CALL TRSUMD(PIE,DVRHO,NRMAX,PIESUM)
      PRLT = PRLSUM*DR/1.D6
      PCXT = PCXSUM*DR/1.D6
      PIET = PIESUM*DR/1.D6
C
      CALL TRSUMD(AJNB,DSRHO,NRMAX,ANBSUM)
      CALL TRSUMD(AJRF,DSRHO,NRMAX,ARFSUM)
      CALL TRSUMD(AJBS,DSRHO,NRMAX,ABSSUM)
C      AJNBT = ANBSUM*DR/1.D6
C      AJRFT = ARFSUM*DR/1.D6
C      AJBST = ABSSUM*DR/1.D6
      AJNBT = ANBSUM*DR/1.D6*(RKAP/FKAP)
      AJRFT = ARFSUM*DR/1.D6*(RKAP/FKAP)
      AJBST = ABSSUM*DR/1.D6*(RKAP/FKAP)
C
      CALL TRSUMD(AJ  ,DSRHO,NRMAX,AJTSUM)
      CALL TRSUMD(AJOH,DSRHO,NRMAX,AOHSUM)
C
C      AJT   = AJTSUM*DR/1.D6
C      AJOHT = AOHSUM*DR/1.D6
      AJT   = AJTSUM*DR/1.D6*(RKAP/FKAP)
      AJOHT = AOHSUM*DR/1.D6*(RKAP/FKAP)
C
C      DO NR=1,NRMAX
C         write(6,*) NR,AJ(NR)
C      ENDDO
C      write(6,*) AJT,DR
C
      DRH=0.5D0*DR*RA
      DO 70 NS=1,NSM
         VNP=AV(NRMAX,NS)
         DNP=AD(NRMAX,NS)
C
         VTP=(AVK(NRMAX,NS)/1.5D0+AV(NRMAX,NS))
         DTP= AK(NRMAX,NS)/(1.5D0)
C
         VXP=0.D0
         DXP=(1.5D0*AD(NRMAX,NS)-AK(NRMAX,NS))*PTS(NS)
C
         SLT(NS) =((     DNP/DRH)*RN(NRMAX,NS)
     &            +( VNP-DNP/DRH)*PNSS(NS))
     &            *DVRHO(NRMAX)
     &            *(FKAP/RKAP)/(RM(NRMAX)*RA)
C
         PLT(NS) =((     DXP/DRH)*RN(NRMAX,NS)
     &            +(     DTP/DRH)*RN(NRMAX,NS)*RT(NRMAX,NS)*1.5D0
     &            +( VXP-DXP/DRH)*PNSS(NS)
     &            +( VTP-DTP/DRH)*PNSS(NS)*PTS(NS)*1.5D0)
     &            *DVRHO(NRMAX)*RKEV*1.D14
     &            *(FKAP/RKAP)/(RM(NRMAX)*RA)
   70 CONTINUE
C
      CALL TRSUMD(SIE,DVRHO,NRMAX,SIESUM)
      CALL TRSUMD(SNF,DVRHO,NRMAX,SNFSUM)
      CALL TRSUMD(SNB,DVRHO,NRMAX,SNBSUM)
      SIET = SIESUM*DR
      SNFT = SNFSUM*DR
      SNBT = SNBSUM*DR
C
      DO 80 NS=1,NSM
         CALL TRSUMD(SPE(1,NS),DVRHO,NRMAX,SPESUM)
         SPET(NS) = SPESUM*DR/RKAP
   80 CONTINUE
C
      WBULKT=0.D0
      PEXST =0.D0
      PRFST =0.D0
      PLST  =0.D0
      SLST  =0.D0
      DO 100 NS=1,NSM
         WBULKT=WBULKT+WST(NS)
         PEXST =PEXST +PEXT(NS)
         PRFST =PRFST +PRFT(NS)
         PLST  =PLST  +PLT(NS)
         SLST  =SLST  +SLT(NS)
  100 CONTINUE
      WTAILT=0.D0
      DO 110 NF=1,NFM
         WTAILT=WTAILT+WFT(NF)
  110 CONTINUE
C
      WPT =WBULKT+WTAILT
      PINT=POHT+PNBT+PRFST+PNFT+PEXST
      POUT=PLST+PCXT+PIET+PRLT
      SINT=SIET+SNBT
      SOUT=SLST
C
      IF(ABS(T-TPRE).LE.1.D-70) THEN
         WPDOT=0.D0
      ELSE
         WPDOT=(WPT-WPPRE)/(T-TPRE)
      ENDIF
      WPPRE=WPT
      TPRE=T
C
      TAUE1=WPT/PINT
      TAUE2=WPT/(PINT-WPDOT)
C
      SUM=0.D0
      SUP=0.D0
      SUL=0.D0
C
C      SUMS=0.D0
      DO NR=1,NRMAX-1
C         IF(NR.EQ.1) THEN
C            SUMS=SUMS+DVRHO(NR)*DR
C         ELSE
C            SUMS=SUMS+0.5D0*(DVRHO(NR-1)+DVRHO(NR))*DR
C         ENDIF
         SUMS=DVRHO(NR)
         SUML=0.D0
         SUPL=0.D0
         DO NS=1,NSM
            SUML = SUML +RN(NR,NS)*RT(NR,NS)*RKEV*1.D20
            SUPL = SUPL +(RN(NR+1,NS)*RT(NR+1,NS)
     &                   -RN(NR  ,NS)*RT(NR  ,NS))*RKEV*1.D20/(DR*RA)
         ENDDO
         DO NF=1,NFM
            SUML = SUML +RW(NR,NF)*RKEV*1.D20
            SUPL = SUPL +(RW(NR+1,NF)
     &                   -RW(NR  ,NF))*RKEV*1.D20/(DR*RA)
         ENDDO
C
C         SUM = SUM + SUML*DVRHO(NR)*DR
C         SUP = SUP + 0.5D0*SUPL*SUMS*DR
C         SUL = SUL + BP(NR)**2*DSRHO(NR)*DR
C         BETA(NR)   = 2.D0*SUM *AMYU0/(SUMS*BB**2)
C         BETAP(NR)  = 2.D0*SUM *AMYU0/(SUMS*BP(NRMAX)**2)
C         BETAPL(NR) = 2.D0*SUML*AMYU0/(     BP(NRMAX)**2)
C         BETAQ(NR)  =-2.D0*SUP *AMYU0/(SUMS*BP(NR   )**2)
C         SUP = SUP + 0.5D0*SUPL*SUMS*DR
         SUM = SUM + SUML*DVRHO(NR)*DR/(2.D0*PI*RR*RKAP)
         SUP = SUP + 0.5D0*SUPL*SUMS*DR*RM(NR)*RA/(2.D0*RKAP*2.D0*PI*RR)
         SUL = SUL + BP(NR)**2*DSRHO(NR)*DR
     &       *RG(NR)/(2.D0*PI*FKAP*RA*RM(NR)*DR)
         BETA(NR)   = 2.D0*SUM *AMYU0/(SUMS*BB**2
     &                *RG(NR)**2/(2.D0*RKAP*2.D0*PI*RR*RM(NR)))
         BETAL(NR)  = 2.D0*SUML*AMYU0/(     BB**2)
         BETAP(NR)  = 2.D0*SUM *AMYU0/(SUMS*BP(NRMAX)**2
     &                *RG(NR)**2/(2.D0*RKAP*2.D0*PI*RR*RM(NR)))
     &                *RKAP/FKAP**2
         BETAPL(NR) = 2.D0*SUML*AMYU0/(     BP(NRMAX)**2)*RKAP/FKAP**2
         BETAQ(NR)  =-2.D0*SUP *AMYU0/(SUMS*BP(NR   )**2
     &                *RG(NR)**2/(2.D0*RKAP*2.D0*PI*RR*RM(NR)))
     &                *RKAP/FKAP**2
         SUP = SUP + 0.5D0*SUPL*PI*RM(NR+1)*RM(NR+1)*DR*RA**3
      ENDDO
C

      NR=NRMAX
C         SUMS=SUMS+0.5D0*(DVRHO(NR-1)+DVRHO(NR))*DR
         SUMS=DVRHO(NR)
         SUML=0.D0
         SUPL=0.D0
         DO NS=1,NSM
C            SUML = SUML +RN(NR,NS)*RT(NR,NS)*RKEV*1.D20
            SUML = SUML +RN(NR,NS)*RT(NR,NS)*RM(NR)*RA*RKEV*1.D20
            SUPL = SUPL +(PNSS(NS)   *PTS(NS)
     &                   -RN(NR  ,NS)*RT(NR  ,NS))*RKEV*1.D20
     &                   /(DR*RA)*2.D0
         ENDDO
         DO NF=1,NFM
C            SUML = SUML +RW(NR,NF)*RKEV*1.D20
            SUML = SUML +RW(NR,NF)*RM(NR)*RA*RKEV*1.D20
            SUPL = SUPL +(0.D0-RW(NR  ,NF))*RKEV*1.D20/(DR*RA)*2.D0
         ENDDO
C
C         SUM = SUM + SUML*DVRHO(NR)*DR
C         SUP = SUP + 0.5D0*SUPL*SUMS*DR
C         SUL = SUL + 0.5D0*BP(NR)**2*DSRHO(NR)*DR
C         BETA(NR)   = 2.D0*SUM *AMYU0/(SUMS*BB**2)
C         BETAL(NR)  = 2.D0*SUML*AMYU0/(     BB **2)
C         BETAP(NR)  = 2.D0*SUM *AMYU0/(SUMS*BP(NRMAX)**2)
C         BETAPL(NR) = 2.D0*SUML*AMYU0/(     BP(NRMAX)**2)
C         BETAQ(NR)  =-2.D0*SUP *AMYU0/(SUMS*BP(NR   )**2)
         SUM = SUM + SUML*DVRHO(NR)*DR/(2.D0*PI*RR*RKAP)
         SUP = SUP + 0.5D0*SUPL*SUMS*DR*RM(NR)*RA/(2.D0*RKAP*2.D0*PI*RR)
         SUL = SUL + 0.5D0*BP(NR)**2*DSRHO(NR)*DR
     &       *RG(NR)/(2.D0*PI*FKAP*RA*RM(NR)*DR)
         BETA(NR)   = 2.D0*SUM *AMYU0/(SUMS*BB**2
     &                *RG(NR)**2/(2.D0*RKAP*2.D0*PI*RR*RM(NR)))
         BETAL(NR)  = 2.D0*SUML*AMYU0/(     BB**2)
         BETAP(NR)  = 2.D0*SUM *AMYU0/(SUMS*BP(NRMAX)**2
     &                *RG(NR)**2/(2.D0*RKAP*2.D0*PI*RR*RM(NR)))
     &                *RKAP/FKAP**2
         BETAPL(NR) = 2.D0*SUML*AMYU0/(     BP(NRMAX)**2)*RKAP/FKAP**2
         BETAQ(NR)  =-2.D0*SUP *AMYU0/(SUMS*BP(NR   )**2
     &                *RG(NR)**2/(2.D0*RKAP*2.D0*PI*RR*RM(NR)))
     &                *RKAP/FKAP**2
C
      BETA0 =(4.D0*BETA(1)  -BETA(2)  )/3.D0
      BETAP0=(4.D0*BETAPL(1)-BETAPL(2))/3.D0
      BETAQ0=0.D0
C
      BETAPA=BETAP(NRMAX)
      BETAA =BETA(NRMAX)
      BETAN =BETAA*1.D2/(RIP/(RA*BB))
C
C      ALI=4.D0*PI*SUL/((AMYU0*AJT*1.D6)**2)
      ALI=8.D0*PI**2*DR*RA*SUL*FKAP**2/((AMYU0*AJT*1.D6)**2)
      VLOOP = EZOH(NRMAX)*2.D0*PI*RR
C
      PAI=(PA(2)*PN(2)+PA(3)*PN(3)+PA(4)*PN(4))/(PN(2)+PN(3)+PN(4))
C
      TAUEP=4.8D-2*(RIP**0.85D0)
     &            *(RR**1.2D0)
     &            *(RA**0.3D0)
     &            *(RKAP**0.5D0)
     &            *(ANSAV(1)**0.1D0)
     &            *(BB**0.2D0)
     &            *(PAI**0.5D0)
     &            *(PINT**(-0.5D0))
C
      QF=5.D0*PNFT/(POHT+PNBT+PRFST+PEXST)
C
      IF(Q0.GE.1.D0) THEN
         RQ1=0.D0
      ELSE
         IF(QP(1).GT.1.D0) THEN
            RQ1=SQRT((1.D0-Q0)/(QP(1)-Q0))*DR
            GOTO 310
         ENDIF
         DO 300 NR=2,NRMAX
            IF(QP(NR).GT.1.D0) THEN
              RQ1=(RG(NR)-RG(NR-1))*RA*(1.D0-QP(NR-1))/(QP(NR)-QP(NR-1))
     &            +RG(NR-1)*RA
               GOTO 310
            ENDIF
  300    CONTINUE
         RQ1=RA
  310    CONTINUE
      ENDIF
C
      ZEFF0=(4.D0*ZEFF(1)-ZEFF(2))/3.D0
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
      SUBROUTINE TRSUMD(A,B,NMAX,SUM)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION A(NMAX),B(NMAX)
C
      SUM=0.D0
      DO 100 N=1,NMAX
         SUM=SUM+A(N)*B(N)
  100 CONTINUE
      RETURN
      END
C
C     ********************************************************
C
C           RADIAL INTEGRATION  (TRIPLE VARIABLES)
C
C     ********************************************************
C
      SUBROUTINE TRSUMT(A,B,C,NMAX,SUM)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION A(NMAX),B(NMAX),C(NMAX)
C
      SUM=0.D0
      DO 100 N=1,NMAX
         SUM=SUM+A(N)*B(N)*C(N)
  100 CONTINUE
      RETURN
      END
C
C     ***********************************************************
C
C           SAVE DATA TO DATAT FOR GRAPHICS
C
C     ***********************************************************
C
      SUBROUTINE TRATOT
C
      INCLUDE 'trcomm.h'
C
      IF(NGT.GE.NTM) RETURN
      NGT=NGT+1
C
      GT    (NGT) = GCLIP(T)
C     
      GVT(NGT, 1) = GCLIP(ANS0(1))
      GVT(NGT, 2) = GCLIP(ANS0(2))
      GVT(NGT, 3) = GCLIP(ANS0(3))
      GVT(NGT, 4) = GCLIP(ANS0(4))
      GVT(NGT, 5) = GCLIP(ANSAV(1))
      GVT(NGT, 6) = GCLIP(ANSAV(2))
      GVT(NGT, 7) = GCLIP(ANSAV(3))
      GVT(NGT, 8) = GCLIP(ANSAV(4))
C
      GVT(NGT, 9) = GCLIP(TS0(1))
      GVT(NGT,10) = GCLIP(TS0(2))
      GVT(NGT,11) = GCLIP(TS0(3))
      GVT(NGT,12) = GCLIP(TS0(4))
      GVT(NGT,13) = GCLIP(TSAV(1))
      GVT(NGT,14) = GCLIP(TSAV(2))
      GVT(NGT,15) = GCLIP(TSAV(3))
      GVT(NGT,16) = GCLIP(TSAV(4))
C
      GVT(NGT,17) = GCLIP(WST(1))
      GVT(NGT,18) = GCLIP(WST(2))
      GVT(NGT,19) = GCLIP(WST(3))
      GVT(NGT,20) = GCLIP(WST(4))
C
      GVT(NGT,21) = GCLIP(ANF0(1))
      GVT(NGT,22) = GCLIP(ANF0(2))
      GVT(NGT,23) = GCLIP(ANFAV(1))
      GVT(NGT,24) = GCLIP(ANFAV(2))
      GVT(NGT,25) = GCLIP(TF0(1))
      GVT(NGT,26) = GCLIP(TF0(2))
      GVT(NGT,27) = GCLIP(TFAV(1))
      GVT(NGT,28) = GCLIP(TFAV(2))
C
      GVT(NGT,29) = GCLIP(WFT(1))
      GVT(NGT,30) = GCLIP(WFT(2))
      GVT(NGT,31) = GCLIP(WBULKT)
      GVT(NGT,32) = GCLIP(WTAILT)
      GVT(NGT,33) = GCLIP(WPT)
C
      GVT(NGT,34) = GCLIP(AJT)
      GVT(NGT,35) = GCLIP(AJOHT)
      GVT(NGT,36) = GCLIP(AJNBT)
      GVT(NGT,37) = GCLIP(AJRFT)
      GVT(NGT,38) = GCLIP(AJBST)
C
      GVT(NGT,39) = GCLIP(PINT)
      GVT(NGT,40) = GCLIP(POHT)
      GVT(NGT,41) = GCLIP(PNBT)
      GVT(NGT,42) = GCLIP(PRFT(1))
      GVT(NGT,43) = GCLIP(PRFT(2))
      GVT(NGT,44) = GCLIP(PRFT(3))
      GVT(NGT,45) = GCLIP(PRFT(4))
      GVT(NGT,46) = GCLIP(PNFT)
C
      GVT(NGT,47) = GCLIP(PBINT)
      GVT(NGT,48) = GCLIP(PBCLT(1))
      GVT(NGT,49) = GCLIP(PBCLT(2))
      GVT(NGT,50) = GCLIP(PBCLT(3))
      GVT(NGT,51) = GCLIP(PBCLT(4))
      GVT(NGT,52) = GCLIP(PFINT)
      GVT(NGT,53) = GCLIP(PFCLT(1))
      GVT(NGT,54) = GCLIP(PFCLT(2))
      GVT(NGT,55) = GCLIP(PFCLT(3))
      GVT(NGT,56) = GCLIP(PFCLT(4))
C
      GVT(NGT,57) = GCLIP(POUT)
      GVT(NGT,58) = GCLIP(PCXT)
      GVT(NGT,59) = GCLIP(PIET)
      GVT(NGT,60) = GCLIP(PRLT)
      GVT(NGT,61) = GCLIP(PLT(1))
      GVT(NGT,62) = GCLIP(PLT(2))
      GVT(NGT,63) = GCLIP(PLT(3))
      GVT(NGT,64) = GCLIP(PLT(4))
C
      GVT(NGT,65) = GCLIP(SINT)
      GVT(NGT,66) = GCLIP(SIET)
      GVT(NGT,67) = GCLIP(SNBT)
      GVT(NGT,68) = GCLIP(SNFT)
      GVT(NGT,69) = GCLIP(SOUT)
      GVT(NGT,70) = GCLIP(SLT(1))
      GVT(NGT,71) = GCLIP(SLT(2))
      GVT(NGT,72) = GCLIP(SLT(3))
      GVT(NGT,73) = GCLIP(SLT(4))
C
      GVT(NGT,74) = GCLIP(VLOOP)
      GVT(NGT,75) = GCLIP(ALI)
      GVT(NGT,76) = GCLIP(RQ1)
      GVT(NGT,77) = GCLIP(Q0)
C
      GVT(NGT,78) = GCLIP(WPDOT)
      GVT(NGT,79) = GCLIP(TAUE1)
      GVT(NGT,80) = GCLIP(TAUE2)
      GVT(NGT,81) = GCLIP(TAUEP)
C
      GVT(NGT,82) = GCLIP(BETAP0)
      GVT(NGT,83) = GCLIP(BETAPA)
      GVT(NGT,84) = GCLIP(BETA0)
      GVT(NGT,85) = GCLIP(BETAA)
C
      GVT(NGT,86) = GCLIP(ZEFF0)
      GVT(NGT,87) = GCLIP(QF)
      GVT(NGT,88) = GCLIP(RIP)
      GVT(NGT,89) = GCLIP(PEXST)
C
C     *** FOR 3D ***
C
      DO NR=1,NRMAX
         G3D(NR,NGT, 1) = GCLIP(RT(NR,1))
         G3D(NR,NGT, 2) = GCLIP(RT(NR,2))
         G3D(NR,NGT, 3) = GCLIP(RT(NR,3))
         G3D(NR,NGT, 4) = GCLIP(RT(NR,4))
C
         G3D(NR,NGT, 5) = GCLIP(RN(NR,1))
         G3D(NR,NGT, 6) = GCLIP(RN(NR,2))
         G3D(NR,NGT, 7) = GCLIP(RN(NR,3))
         G3D(NR,NGT, 8) = GCLIP(RN(NR,4))
C
         G3D(NR,NGT, 9) = GCLIP(AJ  (NR))
         G3D(NR,NGT,10) = GCLIP(AJOH(NR))
         G3D(NR,NGT,11) = GCLIP(AJNB(NR))
         G3D(NR,NGT,12) = GCLIP(AJRF(NR))
         G3D(NR,NGT,13) = GCLIP(AJBS(NR))
C
         G3D(NR,NGT,14) = GCLIP(POH(NR)+PNB(NR)+PNF(NR)
     &                         +PEX(NR,1)+PEX(NR,2)+PEX(NR,3)+PEX(NR,4)
     &                         +PRF(NR,1)+PRF(NR,2)+PRF(NR,3)+PRF(NR,4))
         G3D(NR,NGT,15) = GCLIP(POH(NR))
         G3D(NR,NGT,16) = GCLIP(PNB(NR))
         G3D(NR,NGT,17) = GCLIP(PNF(NR))
         G3D(NR,NGT,18) = GCLIP(PRF(NR,1))
         G3D(NR,NGT,19) = GCLIP(PRF(NR,2))
         G3D(NR,NGT,20) = GCLIP(PRF(NR,3))
         G3D(NR,NGT,21) = GCLIP(PRF(NR,4))
         G3D(NR,NGT,22) = GCLIP(PRL(NR))
         G3D(NR,NGT,23) = GCLIP(PCX(NR))
         G3D(NR,NGT,24) = GCLIP(PIE(NR))
         G3D(NR,NGT,25) = GCLIP(PEX(NR,1))
         G3D(NR,NGT,26) = GCLIP(PEX(NR,2))
         G3D(NR,NGT,27) = GCLIP(PEX(NR,3))
         G3D(NR,NGT,28) = GCLIP(PEX(NR,4))
C
         IF (NR.EQ.1) THEN
            G3D(NR,NGT,27) = GCLIP(Q0)
         ELSE
            G3D(NR,NGT,27) = GCLIP(QP(NR))
         ENDIF
         G3D(NR,NGT,28) = GCLIP(EZOH(NR))
         G3D(NR,NGT,29) = GCLIP(BETA(NR))
         G3D(NR,NGT,30) = GCLIP(BETAP(NR))
         G3D(NR,NGT,31) = GCLIP(EZOH(NR)*2.D0*PI*RR)
         G3D(NR,NGT,32) = GCLIP(ETA(NR))
      ENDDO
C      
      RETURN
      END
C
C     ***********************************************************
C
C           SAVE DATA TO DATA FOR GRAPHICS
C
C     ***********************************************************
C
      SUBROUTINE TRATOG
C
      INCLUDE 'trcomm.h'
C
      IF(NGR.GE.NGM) RETURN
      NGR=NGR+1
      GTR(NGR)=GCLIP(T)
C
      DO 10 NR=1,NRMAX
         GVR(NR,NGR, 1)  = GCLIP(RN(NR,1))
         GVR(NR,NGR, 2)  = GCLIP(RN(NR,2))
         GVR(NR,NGR, 3)  = GCLIP(RN(NR,3))
         GVR(NR,NGR, 4)  = GCLIP(RN(NR,4))
         GVR(NR,NGR, 5)  = GCLIP(RT(NR,1))
         GVR(NR,NGR, 6)  = GCLIP(RT(NR,2))
         GVR(NR,NGR, 7)  = GCLIP(RT(NR,3))
         GVR(NR,NGR, 8)  = GCLIP(RT(NR,4))
         GVR(NR+1,NGR, 9)  = GCLIP(QP(NR))
         GVR(NR,NGR,10)  = GCLIP(AJ(NR)*1.D-6)
         GVR(NR,NGR,11)  = GCLIP(EZOH(NR))
         GVR(NR,NGR,12)  = GCLIP(AJOH(NR)*1.D-6)
         GVR(NR,NGR,13)  = GCLIP((AJNB(NR)+AJRF(NR))*1.D-6)
         GVR(NR,NGR,14)  = GCLIP(AJBS(NR)*1.D-6)
         GVR(NR,NGR,15)  = GCLIP((PIN(NR,1)+PIN(NR,2)
     &                           +PIN(NR,3)+PIN(NR,4))*1.D-6)
         GVR(NR,NGR,16)  = GCLIP(POH(NR)*1.D-6)
         GVR(NR,NGR,17)  = GCLIP(VGR1(NR,2))
         GVR(NR,NGR,18)  = GCLIP(VGR1(NR,1))
         GVR(NR,NGR,19)  = GCLIP(VGR1(NR,3))
C         GVR(NR,NGR,19)  = GCLIP(VGR3(NR,1))
         GVR(NR,NGR,20)  = GCLIP(AK(NR,2))
C         GVR(NR,NGR,17)  = GCLIP(RW(NR,1)*1.D-6*1.5D0)
C         GVR(NR,NGR,18)  = GCLIP(RW(NR,2)*1.D-6*1.5D0)
C         GVR(NR,NGR,19)  = GCLIP(PNB(NR)*1.D-6)
C         GVR(NR,NGR,20)  = GCLIP(PNF(NR)*1.D-6)
   10 CONTINUE
         GVR(1,NGR, 9)  = GCLIP(Q0)
C
      RETURN
      END
C
C     ***********************************************************
C
C           PRINT GLOBAL QUANTITIES
C
C     ***********************************************************
C
      SUBROUTINE TRPRNT(KID)
C
      INCLUDE 'trcomm.h'
C
      CHARACTER KID*1
      CHARACTER K1*3,K2*3,K3*3,K4*3,K5*3,K6*3
      CHARACTER KCOM*40
C
      IF(KID.EQ.'1') THEN
         WRITE(6,601) T,
     &                WPT,WBULKT,WTAILT,WPDOT,
     &                TAUE1,TAUE2,TAUEP,QF,
     &                BETAP0,BETAPA,BETA0,BETAA,
     &                Q0,RQ1,ZEFF0,BETAN
  601    FORMAT(' ','# TIME : ',F7.3,' SEC'/
     &          ' ',3X,'WPT   =',1PD10.3,'  WBULKT=',1PD10.3,
     &               '  WTAILT=',1PD10.3,'  WPDOT =',1PD10.3/
     &          ' ',3X,'TAUE1 =',1PD10.3,'  TAUE2 =',1PD10.3,
     &               '  TAUEP =',1PD10.3,'  QF    =',1PD10.3/
     &          ' ',3X,'BETAP0=',1PD10.3,'  BETAPA=',1PD10.3,
     &               '  BETA0 =',1PD10.3,'  BETAA =',1PD10.3/
     &          ' ',3X,'Q0    =',1PD10.3,'  RQ1   =',1PD10.3,
     &               '  ZEFF0 =',1PD10.3,'  BETAN =',1PD10.3)
C
         WRITE(6,602) WST(1),TS0(1),TSAV(1),ANSAV(1),
     &                WST(2),TS0(2),TSAV(2),ANSAV(2),
     &                WST(3),TS0(3),TSAV(3),ANSAV(3),
     &                WST(4),TS0(4),TSAV(4),ANSAV(4),
     &                WFT(1),TF0(1),TFAV(1),ANFAV(1),
     &                WFT(2),TF0(2),TFAV(2),ANFAV(2)
  602    FORMAT(' ',3X,'WE    =',1PD10.3,'  TE0   =',1PD10.3,
     &               '  TEAVE =',1PD10.3,'  NEAVE =',1PD10.3/
     &          ' ',3X,'WD    =',1PD10.3,'  TD0   =',1PD10.3,
     &               '  TDAVE =',1PD10.3,'  NDAVE =',1PD10.3/
     &          ' ',3X,'WT    =',1PD10.3,'  TT0   =',1PD10.3,
     &               '  TTAVE =',1PD10.3,'  NTAVE =',1PD10.3/
     &          ' ',3X,'WA    =',1PD10.3,'  TA0   =',1PD10.3,
     &               '  TAAVE =',1PD10.3,'  NAAVE =',1PD10.3/
     &          ' ',3X,'WB    =',1PD10.3,'  TB0   =',1PD10.3,
     &               '  TBAVE =',1PD10.3,'  NBAVE =',1PD10.3/
     &          ' ',3X,'WF    =',1PD10.3,'  TF0   =',1PD10.3,
     &               '  TFAVE =',1PD10.3,'  NFAVE =',1PD10.3)
C
         WRITE(6,603) AJT,VLOOP,ALI,VSEC,
     &                AJOHT,AJNBT,AJRFT,AJBST
  603    FORMAT(' ',3X,'AJT   =',1PD10.3,'  VLOOP =',1PD10.3,
     &               '  ALI   =',1PD10.3,'  VSEC  =',1PD10.3/
     &          ' ',3X,'AJOHT =',1PD10.3,'  AJNBT =',1PD10.3,
     &               '  AJRFT =',1PD10.3,'  AJBST =',1PD10.3)
C
         WRITE(6,604) PINT,POHT,PNBT,PNFT,
     &                PRFT(1),PRFT(2),PRFT(3),PRFT(4),
     &                PBINT,PFINT,
     &                PBCLT(1),PBCLT(2),PBCLT(3),PBCLT(4),
     &                PFCLT(1),PFCLT(2),PFCLT(3),PFCLT(4),
     &                POUT,PRLT,PCXT,PIET,
     &                PLT(1),PLT(2),PLT(3),PLT(4)
  604    FORMAT(' ',3X,'PINT  =',1PD10.3,'  POHT  =',1PD10.3,
     &               '  PNBT  =',1PD10.3,'  PNFTE =',1PD10.3/
     &          ' ',3X,'PRFTE =',1PD10.3,'  PRFTD =',1PD10.3,
     &               '  PRFTT =',1PD10.3,'  PRFTA =',1PD10.3/
     &          ' ',3X,'PBIN  =',1PD10.3,'  PFIN  =',1PD10.3/
     &          ' ',3X,'PBCLE =',1PD10.3,'  PBCLD =',1PD10.3,
     &               '  PBCLT =',1PD10.3,'  PBCLA =',1PD10.3/
     &          ' ',3X,'PFCLE =',1PD10.3,'  PFCLD =',1PD10.3,
     &               '  PFCLT =',1PD10.3,'  PFCLA =',1PD10.3/
     &          ' ',3X,'POUT  =',1PD10.3,'  PRLT  =',1PD10.3,
     &               '  PCXT  =',1PD10.3,'  PIETE =',1PD10.3/
     &          ' ',3X,'PLTE  =',1PD10.3,'  PLTD  =',1PD10.3,
     &               '  PLTTE =',1PD10.3,'  PLTA  =',1PD10.3)
C
         WRITE(6,605) SINT,SIET,SNBT,SNFT,
     &                SOUT,ZEFF(1),ANC(1),ANFE(1),
     &                SLT(1),SLT(2),SLT(3),SLT(4)
  605    FORMAT(' ',3X,'SINT  =',1PD10.3,'  SIET  =',1PD10.3,
     &               '  SNBT  =',1PD10.3,'  SNFT  =',1PD10.3/
     &          ' ',3X,'SOUT  =',1PD10.3,'  ZEFF0 =',1PD10.3,
     &               '  ANC0  =',1PD10.3,'  ANFE0 =',1PD10.3/
     &          ' ',3X,'SLTET =',1PD10.3,'  SLTD  =',1PD10.3,
     &               '  SLTTT =',1PD10.3,'  SLTA  =',1PD10.3)
      ENDIF
C
      IF(KID.EQ.'2') THEN
         WRITE(6,611) Q0,(QP(I),I=1,NRMAX)
  611    FORMAT(' ','* Q PROFILE *'/
     &         (' ',5F7.3,2X,5F7.3))
      ENDIF
C
      IF(KID.EQ.'3') THEN
         CALL GUTIME(GTCPU2)
         WRITE(6,621) GTCPU2-GTCPU1
  621    FORMAT(' ','# CPU TIME = ',F8.3,' S')
         RETURN
      ENDIF
C
      IF(KID.EQ.'4') THEN
         WRITE(6,631)
  631    FORMAT(' ','#',12X,'FIRST',7X,'MAX',9X,'MIN',9X,'LAST')
         CALL TRMXMN( 1,'  NE0  ')
C         CALL TRMXMN( 2,'  ND0  ')
C         CALL TRMXMN( 3,'  NT0  ')
C         CALL TRMXMN( 4,'  NA0  ')
         CALL TRMXMN( 5,'  NEAV ')
C         CALL TRMXMN( 6,'  NDAV ')
C         CALL TRMXMN( 7,'  NTAV ')
C         CALL TRMXMN( 8,'  NAAV ')
         CALL TRMXMN( 9,'  TE0  ')
         CALL TRMXMN(10,'  TD0  ')
         CALL TRMXMN(11,'  TT0  ')
C         CALL TRMXMN(12,'  TA0  ')
         CALL TRMXMN(13,'  TEAV ')
         CALL TRMXMN(14,'  TDAV ')
C         CALL TRMXMN(15,'  TTAV ')
C         CALL TRMXMN(16,'  TAAV ')
C         CALL TRMXMN(17,'  WE   ')
C         CALL TRMXMN(18,'  WD   ')
C         CALL TRMXMN(19,'  WT   ')
C         CALL TRMXMN(20,'  WA   ')
C
C         CALL TRMXMN(21,'  NB0  ')
C         CALL TRMXMN(22,'  NF0  ')
C         CALL TRMXMN(23,'  NBAV ')
C         CALL TRMXMN(24,'  NFAV ')
C         CALL TRMXMN(25,'  TB0  ')
C         CALL TRMXMN(26,'  TF0  ')
C         CALL TRMXMN(27,'  TBAV ')
C         CALL TRMXMN(28,'  TFAV ')
         CALL TRMXMN(29,'  WB   ')
         CALL TRMXMN(30,'  WF   ')
         CALL TRMXMN(31,' WBULK ')
C         CALL TRMXMN(32,' WTAIL ')
         CALL TRMXMN(33,'  WP   ')
C
C         CALL TRMXMN(34,'  IP   ')
         CALL TRMXMN(35,'  IOH  ')
         CALL TRMXMN(36,'  INB  ')
C         CALL TRMXMN(37,'  IRF  ')
         CALL TRMXMN(38,'  IBS  ')
C
C         CALL TRMXMN(39,'  PIN  ')
         CALL TRMXMN(40,'  POH  ')
         CALL TRMXMN(41,'  PNB  ')
C         CALL TRMXMN(42,'  PRFE ')
C         CALL TRMXMN(43,'  PRFD ')
C         CALL TRMXMN(44,'  PRFT ')
C         CALL TRMXMN(45,'  PRFA ')
         CALL TRMXMN(46,'  PNF  ')
C         CALL TRMXMN(47,'  PBINT')
C         CALL TRMXMN(48,'  PBCLE')
C         CALL TRMXMN(49,'  PBCLD')
C         CALL TRMXMN(50,'  PBCLT')
C         CALL TRMXMN(51,'  PBCLA')
C         CALL TRMXMN(52,'  PFINT')
C         CALL TRMXMN(53,'  PFCLE')
C         CALL TRMXMN(54,'  PFCLD')
C         CALL TRMXMN(55,'  PFCLT')
C         CALL TRMXMN(56,'  PFCLA')
C         CALL TRMXMN(57,'  POUT ')
         CALL TRMXMN(58,'  PCX  ')
         CALL TRMXMN(59,'  PIE  ')
         CALL TRMXMN(60,'  PRL  ')
         CALL TRMXMN(61,'  PLE  ')
         CALL TRMXMN(62,'  PLD  ')
         CALL TRMXMN(63,'  PLT  ')
         CALL TRMXMN(64,'  PLA  ')
C
C         CALL TRMXMN(65,'  SIN  ')
C         CALL TRMXMN(66,'  SIE  ')
C         CALL TRMXMN(67,'  SNB  ')
C         CALL TRMXMN(68,'  SNF  ')
C         CALL TRMXMN(69,'  SOUT ')
C         CALL TRMXMN(70,'  SLE  ')
C         CALL TRMXMN(71,'  SLD  ')
C         CALL TRMXMN(72,'  SLT  ')
C         CALL TRMXMN(73,'  SLA  ')
C
         CALL TRMXMN(74,' VLOOP ')
         CALL TRMXMN(75,'  LI   ')
C         CALL TRMXMN(76,'  RQ1  ')
C         CALL TRMXMN(77,'   Q0  ')
C         CALL TRMXMN(78,' WPDOT ')
C         CALL TRMXMN(79,' TAUE1 ')
C         CALL TRMXMN(80,' TAUE2 ')
C         CALL TRMXMN(81,' TAUEP ')
C         CALL TRMXMN(82,' BETAP0')
         CALL TRMXMN(83,' BETAPA')
C         CALL TRMXMN(84,' BETA0 ')
C         CALL TRMXMN(85,' BETAA ')
C         CALL TRMXMN(86,' ZEFF0 ')
         CALL TRMXMN(87,'   QF  ')
C         CALL TRMXMN(88,'   IP  ')
         CALL TRMXMN(89,'  PEX  ')
      ENDIF
C
      IF(KID.EQ.'5') THEN
         WRITE(6,641)NGR,NGT,DT,NTMAX
  641    FORMAT(' ','# PARAMETER INFORMATION',/
     &          ' ','  NGR   =',I3,'    NGT   =',I3,/
     &          ' ','  DT    =',1F5.3,'  NTMAX =',I3)
      ENDIF
C
      IF(KID.EQ.'6') THEN
         WRITE(6,651)T,TAUE1,TAUE2,TAUEP,PINT
 651     FORMAT(' ','# TIME : ',F7.3,' SEC'/
     &          ' ',3X,'TAUE1 =',1PD10.3,'  TAUE2 =',1PD10.3,
     &              '  TAUEP =',1PD10.3,'  PINT  =',1PD10.3)
      ENDIF
C
      IF(KID.EQ.'7'.OR.KID.EQ.'8') THEN
         WRITE(6,671) T,
     &                WPT,TAUE1,TAUE2,TAUEP,
     &                BETAN,BETAPA,BETA0,BETAA
  671    FORMAT(' ','# TIME : ',F7.3,' SEC'/
     &          ' ',3X,'WPT   =',1PD10.3,'  TAUE1 =',1PD10.3,
     &               '  TAUE2 =',1PD10.3,'  TAUEP =',1PD10.3/
     &          ' ',3X,'BETAN =',1PD10.3,'  BETAPA=',1PD10.3,
     &               '  BETA0 =',1PD10.3,'  BETAA =',1PD10.3)
C
         WRITE(6,672) WST(1),TS0(1),TSAV(1),ANSAV(1),
     &                WST(2),TS0(2),TSAV(2),ANSAV(2)
  672    FORMAT(' ',3X,'WE    =',1PD10.3,'  TE0   =',1PD10.3,
     &               '  TEAVE =',1PD10.3,'  NEAVE =',1PD10.3/
     &          ' ',3X,'WD    =',1PD10.3,'  TD0   =',1PD10.3,
     &               '  TDAVE =',1PD10.3,'  NDAVE =',1PD10.3)
C
         WRITE(6,673) AJT,VLOOP,ALI,Q0,
     &                AJOHT,AJNBT,AJRFT,AJBST
  673    FORMAT(' ',3X,'AJT   =',1PD10.3,'  VLOOP =',1PD10.3,
     &               '  ALI   =',1PD10.3,'  Q0    =',1PD10.3/
     &          ' ',3X,'AJOHT =',1PD10.3,'  AJNBT =',1PD10.3,
     &               '  AJRFT =',1PD10.3,'  AJBST =',1PD10.3)
C
         WRITE(6,674) PINT,POHT,PNBT,
     &                PRFT(1)+PRFT(2)+PRFT(3)+PRFT(4),
     &                POUT,PRLT,PCXT,PIET
  674    FORMAT(' ',3X,'PINT  =',1PD10.3,'  POHT  =',1PD10.3,
     &               '  PNBT  =',1PD10.3,'  PRFT  =',1PD10.3/
     &          ' ',3X,'POUT  =',1PD10.3,'  PRLT  =',1PD10.3,
     &               '  PCXT  =',1PD10.3,'  PIETE =',1PD10.3)
C
      IF(KID.EQ.'8') THEN
 1600    WRITE(6,*) '## INPUT COMMENT FOR trn.data (A40)'
         READ(5,'(A40)',END=9000,ERR=1600) KCOM
C
C         OPEN(16,POSITION='APPEND',FILE=KFNLOG)
         OPEN(16,ACCESS='APPEND',FILE=KFNLOG)
C
         CALL GUDATE(NDY,NDM,NDD,NTH1,NTM1,NTS1)
         WRITE(K1,'(I3)') 100+NDY
         WRITE(K2,'(I3)') 100+NDM
         WRITE(K3,'(I3)') 100+NDD
         WRITE(K4,'(I3)') 100+NTH1
         WRITE(K5,'(I3)') 100+NTM1
         WRITE(K6,'(I3)') 100+NTS1
         WRITE(16,1670) K1(2:3),K2(2:3),K3(2:3),K4(2:3),K5(2:3),K6(2:3),
     &                  KCOM,
     &                  RIPS,RIPE,PN(1),PN(2),BB,PICTOT,PLHTOT,PLHNPR
 1670    FORMAT(' '/
     &          ' ','## DATE : ',
     &              A2,'-',A2,'-',A2,'  ',A2,':',A2,':',A2,' : ',A40/
     &          ' ',3X,'RIPS  =',1PD10.3,'  RIPE  =',1PD10.3,
     &               '  PNE   =',1PD10.3,'  PNI   =',1PD10.3/
     &          ' ',3X,'BB    =',1PD10.3,'  PICTOT=',1PD10.3,
     &               '  PLHTOT=',1PD10.3,'  PLHNPR=',1PD10.3)
         WRITE(16,1671) T,
     &                WPT,TAUE1,TAUE2,TAUEP,
     &                BETAN,BETAPA,BETA0,BETAA
 1671    FORMAT(' ','# TIME : ',F7.3,' SEC'/
     &          ' ',3X,'WPT   =',1PD10.3,'  TAUE  =',1PD10.3,
     &               '  TAUED =',1PD10.3,'  TAUEP =',1PD10.3/
     &          ' ',3X,'BETAN =',1PD10.3,'  BETAPA=',1PD10.3,
     &               '  BETA0 =',1PD10.3,'  BETAA =',1PD10.3)
C
         WRITE(16,1672) WST(1),TS0(1),TSAV(1),ANSAV(1),
     &                WST(2),TS0(2),TSAV(2),ANSAV(2)
 1672    FORMAT(' ',3X,'WE    =',1PD10.3,'  TE0   =',1PD10.3,
     &               '  TEAVE =',1PD10.3,'  NEAVE =',1PD10.3/
     &          ' ',3X,'WD    =',1PD10.3,'  TD0   =',1PD10.3,
     &               '  TDAVE =',1PD10.3,'  NDAVE =',1PD10.3)
C
         WRITE(16,1673) AJT,VLOOP,ALI,Q0,
     &                AJOHT,AJNBT,AJRFT,AJBST
 1673    FORMAT(' ',3X,'AJT   =',1PD10.3,'  VLOOP =',1PD10.3,
     &               '  ALI   =',1PD10.3,'  Q0    =',1PD10.3/
     &          ' ',3X,'AJOHT =',1PD10.3,'  AJNBT =',1PD10.3,
     &               '  AJRFT =',1PD10.3,'  AJBST =',1PD10.3)
C
         WRITE(16,1674) PINT,POHT,PNBT,
     &                PRFT(1)+PRFT(2)+PRFT(3)+PRFT(4),
     &                POUT,PRLT,PCXT,PIET
 1674    FORMAT(' ',3X,'PINT  =',1PD10.3,'  POHT  =',1PD10.3,
     &               '  PNBT  =',1PD10.3,'  PRFT  =',1PD10.3/
     &          ' ',3X,'POUT  =',1PD10.3,'  PRLT  =',1PD10.3,
     &               '  PCXT  =',1PD10.3,'  PIETE =',1PD10.3)
         CLOSE(16)
      ENDIF
      ENDIF
C
      IF(KID.EQ.'9') THEN
         CALL TRDATA
      ENDIF
C
 9000 RETURN
      END
C
C     ***********************************************************
C
C           PRINT LOCAL DATA
C
C     ***********************************************************
C
      SUBROUTINE TRDATA
C
      INCLUDE 'trcomm.h'
C
    1 WRITE(6,*) '## INPUT MODE : 1:GVT(NT)  2:GVR(NR)  3:GVR(NG)'
      READ(5,*,END=9000,ERR=1) NID
      IF(NID.EQ.0) GOTO 9000
C
      IF(NID.EQ.1) THEN
   10    WRITE(6,*) '## INPUT NID,NTMIN,NTMAX,NTSTEP'
         READ(5,*,END=1,ERR=10) MID,MTMIN,MTMAX,MTSTEP
         IF(MID.EQ.0) GOTO 1
         DO 1000 NT=MTMIN,MTMAX,MTSTEP
            WRITE(6,601) NT,GT(NT),GVT(NT,MID)
 1000    CONTINUE
         GOTO 10
      ELSEIF(NID.EQ.2) THEN
   20    WRITE(6,*) '## INPUT NID,NG,NRMIN,NRMAX,NRSTEP'
         READ(5,*,END=1,ERR=20) MID,NG,MRMIN,MRMAX,MRSTEP
         IF(MID.EQ.0) GOTO 1
         DO 2000 NR=MRMIN,MRMAX,MRSTEP
            WRITE(6,602) NG,NR,GRM(NR),GVR(NR,NG,MID)
 2000    CONTINUE
         GOTO 20
      ELSEIF(NID.EQ.3) THEN
   30    WRITE(6,*) '## INPUT NID,NR,NGMIN,NGMAX,NGSTEP'
         READ(5,*,END=1,ERR=30) MID,NR,MGMIN,MGMAX,MGSTEP
         IF(MID.EQ.0) GOTO 1
         DO 3000 NG=MGMIN,MGMAX,MGSTEP
            WRITE(6,602) NG,NR,GRM(NR),GVR(NR,NG,MID)
 3000    CONTINUE
         GOTO 30
      ENDIF
      GOTO 1
C
 9000 RETURN
  601 FORMAT(' ','  NT=',I3,'  T=',1PE12.4,'  DATA=',1PE12.4)
  602 FORMAT(' ','  NG=',I3,'  NR=',I3,
     &                      '  R=',1PE12.4,    '  DATA=',1PE12.4)
      END
C
C     ***********************************************************
C
C          PEAK-VALUE WO SAGASE
C
C     ***********************************************************
C
      SUBROUTINE TRMXMN(N,STR)
C
      INCLUDE 'trcomm.h'
C
      CHARACTER STR*7
C
      GVMAX=GVT(1,N)
      GVMIN=GVT(1,N)
C
      DO 100 NT=2,NGT
         GVMAX=MAX(GVMAX,GVT(NT,N))
         GVMIN=MIN(GVMIN,GVT(NT,N))
  100 CONTINUE
C
      WRITE(6,600) STR,GVT(1,N),GVMAX,GVMIN,GVT(NGT,N)
  600 FORMAT(' ',A8,5X,1PD10.3,2X,1PD10.3,2X,1PD10.3,2X,1PD10.3)
C
      RETURN
      END
C
C     ***********************************************************
C
C           SIMPLE STATUS REPORT
C
C     ***********************************************************
C
      SUBROUTINE TRSNAP
C
      INCLUDE 'trcomm.h'
C
      WRITE(6,601) T,WPT,TAUE1,Q0,RT(1,1),RT(1,2),RT(1,3),RT(1,4)
  601 FORMAT(' ','# T: ',F7.3,'(S)     WP:',F7.2,'(MJ)  ',
     &           '  TAUE:',F7.3,'(S)   Q0:',F7.3,/
     &       ' ','  TE:',F7.3,'(KEV)   TD:',F7.3,'(KEV) ',
     &           '  TT:',F7.3,'(KEV)   TA:',F7.3,'(KEV)')
      RETURN
      END
C
C
C     ***********************************************************
C
C           SAVE PROFILE DATA
C
C     ***********************************************************
C
      SUBROUTINE TRXOUT
C
      INCLUDE 'trcomm.h'
C      INCLUDE 'trxcom.f'
C
      COMMON /TRXDT1/ KXNDEV,KXNDCG,KXNID
      COMMON /TRKID2/ KDIRW1,KDIRW2
      CHARACTER KXNDEV*80,KXNDCG*80,KXNID*80
      CHARACTER KDIRW*80,KDIRW1*80,KDIRW2*80,KFID*80
      DIMENSION GF1(NTM),GF2(NRMP,NTM),GRX(NRM),GRIN(NRMP)
      DIMENSION DIN(NRM),DRM(NRM),DERIV(NRM)
      DIMENSION UTEOUT(4,NRM),UTIOUT(4,NRM),UQOUT(4,NRM)
C
      NRXMAX=NRMAX+1
C
      KXNDEV='X'
C      KXNDCG='11'
      KXNDCG='test'
      KXNID ='sim'
C
      CALL KTRIM(KXNDEV,IKNDEV)
      CALL KTRIM(KXNDCG,IKNDCG)
      CALL KTRIM(KXNID ,IKNID )
      KDIRW='../../tr.new/data/'//KXNDEV(1:IKNDEV)//'/'
     &                          //KXNDCG(1:IKNDCG)//'/'
     &                          //KXNID (1:IKNID )//'/'
      CALL KTRIM(KDIRW,IKDIRW)
      KDIRW1=KDIRW(1:IKDIRW)//KXNDEV(1:IKNDEV)
     &       //'1d'//KXNDCG(1:IKNDCG)//'.'
      KDIRW2=KDIRW(1:IKDIRW)//KXNDEV(1:IKNDEV)
     &       //'2d'//KXNDCG(1:IKNDCG)//'.'
C
c$$$      KFID='LI'
c$$$      DO 1000 NT=1,NGT
c$$$         GF1(NT)=GVT(NT,75)
c$$$ 1000 CONTINUE
c$$$      CALL TRXW1D(KFID,GT,GF1,NTM,NGT)
c$$$C
c$$$      KFID='POHM'
c$$$      DO 1100 NT=1,NGT
c$$$         GF1(NT)=GVT(NT,40)*1.E6
c$$$ 1100 CONTINUE
c$$$      CALL TRXW1D(KFID,GT,GF1,NTM,NGT)
c$$$C
c$$$      KFID='TE0'
c$$$      DO 1200 NT=1,NGT
c$$$         GF1(NT)=GVT(NT,9)*1.E3
c$$$ 1200 CONTINUE
c$$$      CALL TRXW1D(KFID,GT,GF1,NTM,NGT)
c$$$C
c$$$      KFID='TI0'
c$$$      DO 1300 NT=1,NGT
c$$$         GF1(NT)=GVT(NT,10)*1.E3
c$$$ 1300 CONTINUE
c$$$      CALL TRXW1D(KFID,GT,GF1,NTM,NGT)
c$$$C
c$$$      KFID='VSURF'
c$$$      DO 1400 NT=1,NGT
c$$$         GF1(NT)=GVT(NT,74)
c$$$ 1400 CONTINUE
c$$$      CALL TRXW1D(KFID,GT,GF1,NTM,NGT)
c$$$C
c$$$      KFID='WTH'
c$$$      DO 1500 NT=1,NGT
c$$$         GF1(NT)=GVT(NT,31)*1.E6
c$$$ 1500 CONTINUE
c$$$      CALL TRXW1D(KFID,GT,GF1,NTM,NGT)
c$$$C
c$$$      KFID='WTOT'
c$$$      DO 1600 NT=1,NGT
c$$$         GF1(NT)=GVT(NT,33)*1.E6
c$$$ 1600 CONTINUE
c$$$      CALL TRXW1D(KFID,GT,GF1,NTM,NGT)
c$$$C
c$$$      DGR=1.D0/DBLE(NRMAX)
c$$$      DO 2000 NR=1,NRMAX
c$$$         GRX(NR)=DGR*(NR-0.5D0)
c$$$ 2000 CONTINUE
C
      KFID='TE'
      CALL TR_UFILE2D_CREATE(KFID, 1,1.D3 ,NRXMAX,IERR)
C
      KFID='TI'
      CALL TR_UFILE2D_CREATE(KFID, 2,1.D3 ,NRXMAX,IERR)
C
      KFID='NE'
      CALL TR_UFILE2D_CREATE(KFID, 5,1.D20,NRXMAX,IERR)
C
      KFID='CUR'
      CALL TR_UFILE2D_CREATE(KFID, 9,1.D0 ,NRXMAX,IERR)
C
      KFID='CURBS'
      CALL TR_UFILE2D_CREATE(KFID,13,1.D0 ,NRXMAX,IERR)
C
      KFID='POH'
      CALL TR_UFILE2D_CREATE(KFID,15,1.D0 ,NRXMAX,IERR)
C
      KFID='PCX'
      CALL TR_UFILE2D_CREATE(KFID,23,1.D0 ,NRXMAX,IERR)
C
      KFID='PIE'
      CALL TR_UFILE2D_CREATE(KFID,24,1.D0 ,NRXMAX,IERR)
C
      KFID='RFHE'
      CALL TR_UFILE2D_CREATE(KFID,25,1.D0 ,NRXMAX,IERR)
C
      KFID='RFHI'
      CALL TR_UFILE2D_CREATE(KFID,26,1.D0 ,NRXMAX,IERR)
C
      KFID='Q'
      CALL TR_UFILE2D_CREATE(KFID,27,1.D0 ,NRXMAX,IERR)
C
      KFID='V'
      CALL TR_UFILE2D_CREATE(KFID,31,1.D0 ,NRXMAX,IERR)
C
      KFID='ETA_NC'
      CALL TR_UFILE2D_CREATE(KFID,32,1.D0 ,NRXMAX,IERR)
C
      RETURN
      END
C
C     *****
C
      SUBROUTINE TR_UFILE2D_CREATE(KFID,NUM,AMP,NRXMAX,IERR)
C
      INCLUDE 'trcomm.h'
      DIMENSION GF2(NRMP,NTM),GRIN(NRMP)
      DIMENSION DIN(NRM),DRM(NRM),DERIV(NRM),UOUT(4,NRM)
      CHARACTER KFID*80,KERR*80,KERRF*80
C
      DO NT=1,NGT
         DO NR=1,NRMAX
            DRM(NR)=DBLE(GRM(NR))
            DIN(NR)=DBLE(G3D(NR,NT,NUM))
         ENDDO
         DO NR=1,NRM
            DERIV(NR)=0.D0
         ENDDO
         CALL SPL1D(DRM,DIN,DERIV,UOUT,NRMAX,0,IERR)
         WRITE(KERR,'(A,I2,A,I2)')
     &        'XX TRXOUT: SPL1D G3D(',NUM,'): IERR=',IERR
         IF(IERR.NE.0) WRITE(6,*) KERR
         DO NR=1,NRXMAX
            RIN=DBLE(NR-1)/DBLE(NRMAX)
            CALL SPL1DF(RIN,DATOUT,DRM,UOUT,NRMAX,IERR)
            WRITE(KERRF,'(A,I2,A,I2)')
     &           'XX TRXOUT: SPL1DF G3D(',NUM,'): IERR=',IERR
            IF(IERR.NE.0) WRITE(6,*) KERRF
            GF2(NR,NT)=GUCLIP(DATOUT*AMP)
            GRIN(NR)=GUCLIP(RIN)
         ENDDO
      ENDDO
      CALL TRXW2D(KFID,GT,GRIN,GF2,NRMP,NTM,NRXMAX,NGT)
C
      RETURN
      END
C
C     *****
C
      SUBROUTINE TRXW1D(KFID,GT,GF,NTM,NTXMAX)
C
      IMPLICIT REAL*8(A-F,H,O-Z)
      DIMENSION GT(NTM),GF(NTM)
      COMMON /TRXDT1/ KXNDEV,KXNDCG,KXNID
      COMMON /TRKID2/ KDIRW1,KDIRW2
      CHARACTER KXNDEV*80,KXNDCG*80,KXNID*80
      CHARACTER KDIRW1*80,KDIRW2*80,KFID*80,KFILE*80
C
      CALL KTRIM(KDIRW1,KL1)
      KFILE=KDIRW1(1:KL1)//KFID
      WRITE(6,*) '- OPEN FILE:',KFILE(1:55)
C
      OPEN(16,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',ERR=10)
C
      WRITE(16,'(1X,A8,A8,A14,A18)') KXNDCG(1:8),KXNDEV(1:8),
     &     '               ',
     &     ';-SHOT #- DEVICE -'
      WRITE(16,'(1X,A30,A29)')
     &     'TIME          SECONDS         ',
     &     ';-INDEPENDENT VARIABLE LABEL-'
      WRITE(16,'(1X,A30,A27)') KFID,
     &     ';-DEPENDENT VARIABLE LABEL-'
      WRITE(16,'(1X,I30,A33)') NTXMAX,
     &     ';-# OF PTS-  X, F(X) DATA FOLLOW:'
C
      WRITE(16,'(1X,1P6E13.6)') (GT(NTX),NTX=1,NTXMAX)
      WRITE(16,'(1X,1P6E13.6)') (GF(NTX),NTX=1,NTXMAX)
C
      CLOSE(16)
      RETURN
C
   10 WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
      RETURN
      END
C
C     *****
C
      SUBROUTINE TRXW2D(KFID,GT,GR,GF,NRM,NTM,NRXMAX,NTXMAX)
C
      IMPLICIT REAL*8(A-F,H,O-Z)
      DIMENSION GT(NTM),GR(NRM),GF(NRM,NTM)
      COMMON /TRXDT1/ KXNDEV,KXNDCG,KXNID
      COMMON /TRKID2/ KDIRW1,KDIRW2
      CHARACTER KXNDEV*80,KXNDCG*80,KXNID*80
      CHARACTER KDIRW1*80,KDIRW2*80,KFID*80,KFILE*80
C
      CALL KTRIM(KDIRW2,KL2)
      KFILE=KDIRW2(1:KL2)//KFID
      WRITE(6,*) '- OPEN FILE:',KFILE(1:55)
C
      OPEN(16,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',ERR=10)
C
C
      WRITE(16,'(1X,A8,A8,A14,A18)') KXNDCG(1:8),KXNDEV(1:8),
     &     '               ',
     &     ';-SHOT #- DEVICE -'
      WRITE(16,'(1X,A30,A32)')
     &     'RHO                           ',
     &     ';-INDEPENDENT VARIABLE LABEL: X-'
      WRITE(16,'(1X,A30,A32)')
     &     'TIME          SECONDS         ',
     &     ';-INDEPENDENT VARIABLE LABEL: Y-'
      WRITE(16,'(1X,A30,A27)') KFID,
     &     ';-DEPENDENT VARIABLE LABEL-'
      WRITE(16,'(1X,I30,A12)') NRXMAX,
     &     ';-# OF X PTS-:'
      WRITE(16,'(1X,I30,A35)') NTXMAX,
     &     ';-# OF Y PTS-  X, F(X) DATA FOLLOW:'
C
      WRITE(16,'(1X,1P6E13.6)') (GR(NRX),NRX=1,NRXMAX)
      WRITE(16,'(1X,1P6E13.6)') (GT(NTX),NTX=1,NTXMAX)
      WRITE(16,'(1X,1P6E13.6)') ((GF(NRX,NTX),NRX=1,NRXMAX),
     &                                        NTX=1,NTXMAX)
C
      CLOSE(16)
      RETURN
C
   10 WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
      RETURN
      END
