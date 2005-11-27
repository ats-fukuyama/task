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
      INCLUDE 'trcomm.inc'
      DIMENSION DSRHO(NRM)!,FACT(NRM)
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
C
      VOL=0.D0
      DO NR=1,NRMAX
         VOL=VOL+DVRHO(NR)*DR
         DSRHO(NR)=DVRHO(NR)/(2.D0*PI*RMJRHO(NR))
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
      IF(MDLUF.NE.0) THEN
         CALL TRSUMD(RNF(1,1),DVRHO,NRMAX,ANFSUM)
         CALL TRSUMT(RNF(1,1),RT(1,2),DVRHO,NRMAX,RNTSUM)
         CALL TRSUMD(PBM,DVRHO,NRMAX,RWSUM)
         NF=1
         WFT(NF) = RWSUM*DR*1.D-6-1.5D0*RNTSUM*DR*RKEV*1.D14
         ANFAV(NF) = ANFSUM*DR/VOL
         ANF0(NF)  = (9.D0*RNF(1,1)-RNF(2,1))/8.D0
         IF(ANFSUM.GT.0.D0) THEN
            TFAV(NF)  = RWSUM/(RKEV*1.D20)/ANFSUM
         ELSE
            TFAV(NF)  = 0.D0
         ENDIF
         IF(RNF(1,1).GT.0.D0) THEN
            TF0(NF)  = (9.D0*PBM(1)-PBM(2))/8.D0/(RKEV*1.D20)
     &                 /ANF0(NF)
         ELSE
            TF0(NF)  = 0.D0
         ENDIF
         NF=2
         WFT(NF) = 0.D0
         ANFAV(NF) = 0.D0
         ANF0(NF)  = 0.D0
         TFAV(NF)  = 0.D0
         TF0(NF)   = 0.D0
      ELSE
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
      ENDIF
C
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
      DO NS=1,NSM
         CALL TRSUMD(PRFV(1,NS,1),DVRHO,NRMAX,PRFV1SUM)
         CALL TRSUMD(PRFV(1,NS,2),DVRHO,NRMAX,PRFV2SUM)
         CALL TRSUMD(PRFV(1,NS,3),DVRHO,NRMAX,PRFV3SUM)
         CALL TRSUMD(PRF (1,NS  ),DVRHO,NRMAX,PRFSUM  )
         PRFVT(NS,1) = PRFV1SUM*DR/1.D6
         PRFVT(NS,2) = PRFV2SUM*DR/1.D6
         PRFVT(NS,3) = PRFV3SUM*DR/1.D6
         PRFT (NS  ) = PRFSUM  *DR/1.D6
      ENDDO
C
      CALL TRSUMD(PBIN,DVRHO,NRMAX,PBSUM)
      PBINT = PBSUM*DR/1.D6
      DO NS=1,NSM
         CALL TRSUMD(PBCL(1,NS),DVRHO,NRMAX,PBSUM)
         PBCLT(NS) = PBSUM*DR/1.D6
      ENDDO
C
      CALL TRSUMD(PFIN,DVRHO,NRMAX,PFSUM)
      PFINT = PFSUM*DR/1.D6
      DO NS=1,NSM
         CALL TRSUMD(PFCL(1,NS),DVRHO,NRMAX,PFSUM)
         PFCLT(NS) = PFSUM*DR/1.D6
      ENDDO
C     
      CALL TRSUMD(PRL,DVRHO,NRMAX,PRLSUM)
      CALL TRSUMD(PCX,DVRHO,NRMAX,PCXSUM)
      CALL TRSUMD(PIE,DVRHO,NRMAX,PIESUM)
      PRLT = PRLSUM*DR/1.D6
      PCXT = PCXSUM*DR/1.D6
      PIET = PIESUM*DR/1.D6
C
      CALL TRSUMD(AJNB,DSRHO,NRMAX,ANBSUM)
      CALL TRSUMD(AJRFV(1,1),DSRHO,NRMAX,ARF1SUM)
      CALL TRSUMD(AJRFV(1,2),DSRHO,NRMAX,ARF2SUM)
      CALL TRSUMD(AJRFV(1,3),DSRHO,NRMAX,ARF3SUM)
      CALL TRSUMD(AJRF,DSRHO,NRMAX,ARFSUM)
      CALL TRSUMD(AJBS,DSRHO,NRMAX,ABSSUM)
c$$$      DO NR=1,NRMAX
c$$$         FACT(NR)=BB/(BB+BP(NR))
c$$$      ENDDO
c$$$      CALL TRSUMT(AJNB,DSRHO,FACT,NRMAX,ANBSUM)
c$$$      CALL TRSUMT(AJRFV(1,1),DSRHO,FACT,NRMAX,ARF1SUM)
c$$$      CALL TRSUMT(AJRFV(1,2),DSRHO,FACT,NRMAX,ARF2SUM)
c$$$      CALL TRSUMT(AJRFV(1,3),DSRHO,FACT,NRMAX,ARF3SUM)
c$$$      CALL TRSUMT(AJRF,DSRHO,FACT,NRMAX,ARFSUM)
c$$$      CALL TRSUMT(AJBS,DSRHO,FACT,NRMAX,ABSSUM)
      AJNBT = ANBSUM*DR/1.D6
      AJRFVT(1) = ARF1SUM*DR/1.D6
      AJRFVT(2) = ARF2SUM*DR/1.D6
      AJRFVT(3) = ARF3SUM*DR/1.D6
      AJRFT = ARFSUM*DR/1.D6
      AJBST = ABSSUM*DR/1.D6
C
      CALL TRSUMD(AJ   ,DSRHO,NRMAX,AJTSUM )
      CALL TRSUMD(AJTOR,DSRHO,NRMAX,AJTTSUM)
      CALL TRSUMD(AJOH ,DSRHO,NRMAX,AOHSUM )
c$$$      CALL TRSUMT(AJOH,DSRHO,FACT,NRMAX,AOHSUM)
C
      AJT    = AJTSUM *DR/1.D6
      AJTTOR = AJTTSUM*DR/1.D6
      AJOHT  = AOHSUM *DR/1.D6
C
      IF(RHOA.EQ.1.D0) THEN
         NRL=NRMAX
      ELSE
         NRL=NRAMAX
      ENDIF
      NSW=3
      NRMAX=NRAMAX
      CALL TR_COEF_DECIDE(NRL,NSW,DV53)
      NRMAX=NROMAX
      NMK=2
      DRH=DR/DVRHO(NRL)**(2.D0/3.D0)
      DO NEQ=1,NEQMAX
         NSSN=NSS(NEQ)
         NSVN=NSV(NEQ)
         SUM=0.D0
         DO NW=1,NEQMAX
            NSSN1=NSS(NW)
            NSVN1=NSV(NW)
            IF(NSVN1.EQ.1) THEN
               SUM=SUM+  DD(NEQ,NW,NMK,NSW)*2.D0 *DRH *RN(NRL,NSSN1)
     &                +( VV(NEQ,NW,NMK,NSW)      *DRH
     &                  -DD(NEQ,NW,NMK,NSW)*2.D0)*DRH *PNSS(NSSN1)
            ELSEIF(NSVN1.EQ.2) THEN
               SUM=SUM+  DD(NEQ,NW,NMK,NSW)*2.D0 *DRH *RN(NRL,NSSN1)
     &                                                *RT(NRL,NSSN1)
     &                +( VV(NEQ,NW,NMK,NSW)      *DRH
     &                  -DD(NEQ,NW,NMK,NSW)*2.D0 *DRH)*PNSS(NSSN1)
     &                                                *PTS (NSSN1)
            ENDIF
         ENDDO
         IF(NSVN.EQ.1) THEN
            SLT(NSSN)=SUM*RKEV*1.D14
         ELSEIF(NSVN.EQ.2) THEN
            PLT(NSSN)=SUM*RKEV*1.D14
         ENDIF
      ENDDO
C
      CALL TRSUMD(SIE,DVRHO,NRMAX,SIESUM)
      CALL TRSUMD(SNF,DVRHO,NRMAX,SNFSUM)
      CALL TRSUMD(SNB,DVRHO,NRMAX,SNBSUM)
      SIET = SIESUM*DR
      SNFT = SNFSUM*DR
      SNBT = SNBSUM*DR
C
      DO NS=1,NSM
         CALL TRSUMD(SPE(1,NS),DVRHO,NRMAX,SPESUM)
         SPET(NS) = SPESUM*DR/RKAP
      ENDDO
C
      WBULKT=0.D0
      PEXST =0.D0
      PRFST =0.D0
      PLST  =0.D0
      SLST  =0.D0
      DO NS=1,NSM
         WBULKT=WBULKT+WST(NS)
         PEXST =PEXST +PEXT(NS)
         PRFST =PRFST +PRFT(NS)
         PLST  =PLST  +PLT(NS)
         SLST  =SLST  +SLT(NS)
      ENDDO
      WTAILT=0.D0
      DO NF=1,NFM
         WTAILT=WTAILT+WFT(NF)
      ENDDO
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
      VOL=0.D0
      DO NR=1,NRMAX-1
         SUML=0.D0
         DO NS=1,NSM
            SUML = SUML +RN(NR,NS)*RT(NR,NS)*RKEV*1.D20
         ENDDO
         DO NF=1,NFM
            SUML = SUML +RW(NR,NF)*RKEV*1.D20
         ENDDO
C
         VOL = VOL +         DVRHO(NR)*DR
         SUM = SUM + SUML   *DVRHO(NR)*DR
         SUP = SUP + SUML**2*DVRHO(NR)*DR
         BETA(NR)   = 2.D0*RMU0*SUM      /(     VOL *BB**2)
         BETAL(NR)  = 2.D0*RMU0*SUML     /(          BB**2)
         BETAP(NR)  = 2.D0*RMU0*SUM      /(     VOL *BP(NRMAX)**2)
         BETAPL(NR) = 2.D0*RMU0*SUML     /(          BP(NRMAX)**2)
         BETAQ(NR)  = 2.D0*RMU0*SQRT(SUP)/(SQRT(VOL)*BB**2)
      ENDDO
C

      NR=NRMAX
         SUML=0.D0
         DO NS=1,NSM
            SUML = SUML +RN(NR,NS)*RT(NR,NS)*RKEV*1.D20
         ENDDO
         DO NF=1,NFM
            SUML = SUML +RW(NR,NF)*RKEV*1.D20
         ENDDO
C
         VOL = VOL +         DVRHO(NR)*DR
         SUM = SUM + SUML   *DVRHO(NR)*DR
         SUP = SUP + SUML**2*DVRHO(NR)*DR
         BETA(NR)   = 2.D0*RMU0*SUM      /(     VOL *BB**2)
         BETAL(NR)  = 2.D0*RMU0*SUML     /(          BB**2)
         BETAP(NR)  = 2.D0*RMU0*SUM      /(     VOL *BP(NRMAX)**2)
         BETAPL(NR) = 2.D0*RMU0*SUML     /(          BP(NRMAX)**2)
         BETAQ(NR)  = 2.D0*RMU0*SQRT(SUP)/(SQRT(VOL)*BB**2)
C
      BETA0 =(4.D0*BETA(1)  -BETA(2)  )/3.D0
      BETAP0=(4.D0*BETAPL(1)-BETAPL(2))/3.D0
      BETAQ0=(4.D0*BETAQ(1) -BETAQ(2) )/3.D0
C
      BETAPA=BETAP(NRMAX)
      BETAA =BETA(NRMAX)
      BETAN =BETAA*1.D2/(RIP/(RA*BB))
C
      WPOL=0.D0
      DO NR=1,NRMAX
         WPOL=WPOL+DVRHOG(NR)*ABRHOG(NR)*RDP(NR)**2*DR/(2.D0*RMU0)
      ENDDO
      ALI=4.D0*WPOL/(RMU0*RR*(AJT*1.D6)**2)
C
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
         DO NR=2,NRMAX
            IF(QP(NR).GT.1.D0) THEN
              RQ1=(RG(NR)-RG(NR-1))*RA*(1.D0-QP(NR-1))/(QP(NR)-QP(NR-1))
     &            +RG(NR-1)*RA
               GOTO 310
            ENDIF
         ENDDO
         RQ1=RA
  310    CONTINUE
      ENDIF
C
      ZEFF0=(4.D0*ZEFF(1)-ZEFF(2))/3.D0
C
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
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
      DO N=1,NMAX
         SUM=SUM+A(N)*B(N)
      ENDDO
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
      DO N=1,NMAX
         SUM=SUM+A(N)*B(N)*C(N)
      ENDDO
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
      INCLUDE 'trcomm.inc'
      DIMENSION DERIV(NRM),U(4,NRM),U0(NRM)
C
      IF(NGT.GE.NTM) RETURN
      NGT=NGT+1
C
      GT    (NGT) = GUCLIP(T)
C     
      GVT(NGT, 1) = GUCLIP(ANS0(1))
      GVT(NGT, 2) = GUCLIP(ANS0(2))
      GVT(NGT, 3) = GUCLIP(ANS0(3))
      GVT(NGT, 4) = GUCLIP(ANS0(4))
      GVT(NGT, 5) = GUCLIP(ANSAV(1))
      GVT(NGT, 6) = GUCLIP(ANSAV(2))
      GVT(NGT, 7) = GUCLIP(ANSAV(3))
      GVT(NGT, 8) = GUCLIP(ANSAV(4))
C
      GVT(NGT, 9) = GUCLIP(TS0(1))
      GVT(NGT,10) = GUCLIP(TS0(2))
      GVT(NGT,11) = GUCLIP(TS0(3))
      GVT(NGT,12) = GUCLIP(TS0(4))
      GVT(NGT,13) = GUCLIP(TSAV(1))
      GVT(NGT,14) = GUCLIP(TSAV(2))
      GVT(NGT,15) = GUCLIP(TSAV(3))
      GVT(NGT,16) = GUCLIP(TSAV(4))
C
      GVT(NGT,17) = GUCLIP(WST(1))
      GVT(NGT,18) = GUCLIP(WST(2))
      GVT(NGT,19) = GUCLIP(WST(3))
      GVT(NGT,20) = GUCLIP(WST(4))
C
      GVT(NGT,21) = GUCLIP(ANF0(1))
      GVT(NGT,22) = GUCLIP(ANF0(2))
      GVT(NGT,23) = GUCLIP(ANFAV(1))
      GVT(NGT,24) = GUCLIP(ANFAV(2))
      GVT(NGT,25) = GUCLIP(TF0(1))
      GVT(NGT,26) = GUCLIP(TF0(2))
      GVT(NGT,27) = GUCLIP(TFAV(1))
      GVT(NGT,28) = GUCLIP(TFAV(2))
C
      GVT(NGT,29) = GUCLIP(WFT(1))
      GVT(NGT,30) = GUCLIP(WFT(2))
      GVT(NGT,31) = GUCLIP(WBULKT)
      GVT(NGT,32) = GUCLIP(WTAILT)
      GVT(NGT,33) = GUCLIP(WPT)
C
      GVT(NGT,34) = GUCLIP(AJT)
      GVT(NGT,35) = GUCLIP(AJOHT)
      GVT(NGT,36) = GUCLIP(AJNBT)
      GVT(NGT,37) = GUCLIP(AJRFT)
      GVT(NGT,38) = GUCLIP(AJBST)
C
      GVT(NGT,39) = GUCLIP(PINT)
      GVT(NGT,40) = GUCLIP(POHT)
      GVT(NGT,41) = GUCLIP(PNBT)
      GVT(NGT,42) = GUCLIP(PRFT(1))
      GVT(NGT,43) = GUCLIP(PRFT(2))
      GVT(NGT,44) = GUCLIP(PRFT(3))
      GVT(NGT,45) = GUCLIP(PRFT(4))
      GVT(NGT,46) = GUCLIP(PNFT)
C
      GVT(NGT,47) = GUCLIP(PBINT)
      GVT(NGT,48) = GUCLIP(PBCLT(1))
      GVT(NGT,49) = GUCLIP(PBCLT(2))
      GVT(NGT,50) = GUCLIP(PBCLT(3))
      GVT(NGT,51) = GUCLIP(PBCLT(4))
      GVT(NGT,52) = GUCLIP(PFINT)
      GVT(NGT,53) = GUCLIP(PFCLT(1))
      GVT(NGT,54) = GUCLIP(PFCLT(2))
      GVT(NGT,55) = GUCLIP(PFCLT(3))
      GVT(NGT,56) = GUCLIP(PFCLT(4))
C
      GVT(NGT,57) = GUCLIP(POUT)
      GVT(NGT,58) = GUCLIP(PCXT)
      GVT(NGT,59) = GUCLIP(PIET)
      GVT(NGT,60) = GUCLIP(PRLT)
      GVT(NGT,61) = GUCLIP(PLT(1))
      GVT(NGT,62) = GUCLIP(PLT(2))
      GVT(NGT,63) = GUCLIP(PLT(3))
      GVT(NGT,64) = GUCLIP(PLT(4))
C
      GVT(NGT,65) = GUCLIP(SINT)
      GVT(NGT,66) = GUCLIP(SIET)
      GVT(NGT,67) = GUCLIP(SNBT)
      GVT(NGT,68) = GUCLIP(SNFT)
      GVT(NGT,69) = GUCLIP(SOUT)
      GVT(NGT,70) = GUCLIP(SLT(1))
      GVT(NGT,71) = GUCLIP(SLT(2))
      GVT(NGT,72) = GUCLIP(SLT(3))
      GVT(NGT,73) = GUCLIP(SLT(4))
C
      GVT(NGT,74) = GUCLIP(VLOOP)
      GVT(NGT,75) = GUCLIP(ALI)
      GVT(NGT,76) = GUCLIP(RQ1)
      GVT(NGT,77) = GUCLIP(Q0)
C
      GVT(NGT,78) = GUCLIP(WPDOT)
      GVT(NGT,79) = GUCLIP(TAUE1)
      GVT(NGT,80) = GUCLIP(TAUE2)
      GVT(NGT,81) = GUCLIP(TAUEP)
C
      GVT(NGT,82) = GUCLIP(BETAP0)
      GVT(NGT,83) = GUCLIP(BETAPA)
      GVT(NGT,84) = GUCLIP(BETA0)
      GVT(NGT,85) = GUCLIP(BETAA)
C
      GVT(NGT,86) = GUCLIP(ZEFF0)
      GVT(NGT,87) = GUCLIP(QF)
      GVT(NGT,88) = GUCLIP(RIP)
C 
      GVT(NGT,89) = GUCLIP(PEXT(1))
      GVT(NGT,90) = GUCLIP(PEXT(2))
      GVT(NGT,91) = GUCLIP(PRFVT(1,1)) ! ECH  to electron
      GVT(NGT,92) = GUCLIP(PRFVT(2,1)) ! ECH  to ions
      GVT(NGT,93) = GUCLIP(PRFVT(1,2)) ! LH   to electron
      GVT(NGT,94) = GUCLIP(PRFVT(2,2)) ! LH   to ions
      GVT(NGT,95) = GUCLIP(PRFVT(1,3)) ! ICRH to electron
      GVT(NGT,96) = GUCLIP(PRFVT(2,3)) ! ICRH to ions
C
      GVT(NGT,97) = GUCLIP(RR)
      GVT(NGT,98) = GUCLIP(RA)
      GVT(NGT,99) = GUCLIP(BB)
      GVT(NGT,100)= GUCLIP(RKAP)
      GVT(NGT,101)= GUCLIP(AJTTOR)
C
C     *** FOR 3D ***
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      CALL SPL1D  (RM,DVRHO,DERIV,U,NRMAX,0,IERR)
      CALL SPL1DI0(RM,U,U0,NRMAX,IERR)
      DO NR=1,NRMAX
         G3D(NR,NGT, 1) = GUCLIP(RT(NR,1))
         G3D(NR,NGT, 2) = GUCLIP(RT(NR,2))
         G3D(NR,NGT, 3) = GUCLIP(RT(NR,3))
         G3D(NR,NGT, 4) = GUCLIP(RT(NR,4))
C
         G3D(NR,NGT, 5) = GUCLIP(RN(NR,1))
         G3D(NR,NGT, 6) = GUCLIP(RN(NR,2))
         G3D(NR,NGT, 7) = GUCLIP(RN(NR,3))
         G3D(NR,NGT, 8) = GUCLIP(RN(NR,4))
C
         G3D(NR,NGT, 9) = GUCLIP(AJ  (NR))
         G3D(NR,NGT,10) = GUCLIP(AJOH(NR))
         G3D(NR,NGT,11) = GUCLIP(AJNB(NR))
         G3D(NR,NGT,12) = GUCLIP(AJRF(NR))
         G3D(NR,NGT,13) = GUCLIP(AJBS(NR))
C
         G3D(NR,NGT,14) = GUCLIP(POH(NR)+PNB(NR)+PNF(NR)
     &                         +PEX(NR,1)+PEX(NR,2)+PEX(NR,3)+PEX(NR,4)
     &                         +PRF(NR,1)+PRF(NR,2)+PRF(NR,3)+PRF(NR,4))
         G3D(NR,NGT,15) = GUCLIP(POH(NR))
         G3D(NR,NGT,16) = GUCLIP(PNB(NR))
         G3D(NR,NGT,17) = GUCLIP(PNF(NR))
         G3D(NR,NGT,18) = GUCLIP(PRF(NR,1))
         G3D(NR,NGT,19) = GUCLIP(PRF(NR,2))
         G3D(NR,NGT,20) = GUCLIP(PRF(NR,3))
         G3D(NR,NGT,21) = GUCLIP(PRF(NR,4))
         G3D(NR,NGT,22) = GUCLIP(PRL(NR))
         G3D(NR,NGT,23) = GUCLIP(PCX(NR))
         G3D(NR,NGT,24) = GUCLIP(PIE(NR))
         G3D(NR,NGT,25) = GUCLIP(PEX(NR,1))
         G3D(NR,NGT,26) = GUCLIP(PEX(NR,2))
C
C         IF (NR.EQ.1) THEN
C            G3D(NR,NGT,27) = GUCLIP(Q0)
C         ELSE
            G3D(NR,NGT,27) = GUCLIP(QP(NR))
C         ENDIF
         G3D(NR,NGT,28) = GUCLIP(EZOH(NR))
         G3D(NR,NGT,29) = GUCLIP(BETA(NR))
         G3D(NR,NGT,30) = GUCLIP(BETAP(NR))
         G3D(NR,NGT,31) = GUCLIP(EZOH(NR)*2.D0*PI*RR)
         G3D(NR,NGT,32) = GUCLIP(ETA(NR))
         G3D(NR,NGT,33) = GUCLIP(ZEFF(NR))
         G3D(NR,NGT,34) = GUCLIP(AK(NR,1))
         G3D(NR,NGT,35) = GUCLIP(AK(NR,2))
C
         G3D(NR,NGT,36) = GUCLIP(PRFV(NR,1,1))
         G3D(NR,NGT,37) = GUCLIP(PRFV(NR,1,2))
         G3D(NR,NGT,38) = GUCLIP(PRFV(NR,1,3))
         G3D(NR,NGT,39) = GUCLIP(PRFV(NR,2,1))
         G3D(NR,NGT,40) = GUCLIP(PRFV(NR,2,2))
         G3D(NR,NGT,41) = GUCLIP(PRFV(NR,2,3))
C
         G3D(NR,NGT,42) = GUCLIP(AJRFV(NR,1))
         G3D(NR,NGT,43) = GUCLIP(AJRFV(NR,2))
         G3D(NR,NGT,44) = GUCLIP(AJRFV(NR,3))
C
         G3D(NR,NGT,45) = GUCLIP(RW(NR,1)+RW(NR,2))
         G3D(NR,NGT,46) = GUCLIP(ANC(NR)+ANFE(NR))
         G3D(NR,NGT,47) = GUCLIP(BP(NR))
         G3D(NR,NGT,48) = GUCLIP(RPSI(NR))
C
         G3D(NR,NGT,49) = GUCLIP(RMJRHO(NR))
         G3D(NR,NGT,50) = GUCLIP(RMNRHO(NR))
         RMN=(DBLE(NR)-0.5D0)*DR
         CALL SPL1DI(RMN,F0D,RM,U,U0,NRMAX,IERR)
         G3D(NR,NGT,51) = GUCLIP(F0D)
         G3D(NR,NGT,52) = GUCLIP(RKPRHO(NR))
         G3D(NR,NGT,53) = GUCLIP(1.D0) ! DELTAR
         G3D(NR,NGT,54) = GUCLIP(AR1RHO(NR))
         G3D(NR,NGT,55) = GUCLIP(AR2RHO(NR))
         G3D(NR,NGT,56) = GUCLIP(AKDW(NR,1))
         G3D(NR,NGT,57) = GUCLIP(AKDW(NR,2))
         G3D(NR,NGT,58) = GUCLIP(RN(NR,1)*RT(NR,1))
         G3D(NR,NGT,59) = GUCLIP(RN(NR,2)*RT(NR,2))
C
         G3D(NR,NGT,60) = GUCLIP(VTOR(NR))
         G3D(NR,NGT,61) = GUCLIP(VPOL(NR))
C
         G3D(NR,NGT,62) = GUCLIP(SALPHA(NR))
C
      ENDDO
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
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
      INCLUDE 'trcomm.inc'
C
      IF(NGR.GE.NGM) RETURN
      NGR=NGR+1
      GTR(NGR)=GUCLIP(T)
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      DO NR=1,NRMAX
         GVR(NR,NGR, 1)  = GUCLIP(RN(NR,1))
         GVR(NR,NGR, 2)  = GUCLIP(RN(NR,2))
         GVR(NR,NGR, 3)  = GUCLIP(RN(NR,3))
         GVR(NR,NGR, 4)  = GUCLIP(RN(NR,4))
         GVR(NR,NGR, 5)  = GUCLIP(RT(NR,1))
         GVR(NR,NGR, 6)  = GUCLIP(RT(NR,2))
         GVR(NR,NGR, 7)  = GUCLIP(RT(NR,3))
         GVR(NR,NGR, 8)  = GUCLIP(RT(NR,4))
         GVR(NR+1,NGR, 9)  = GUCLIP(QP(NR))
         GVR(NR,NGR,10)  = GUCLIP(AJ(NR)*1.D-6)
         GVR(NR,NGR,11)  = GUCLIP(EZOH(NR))
         GVR(NR,NGR,12)  = GUCLIP(AJOH(NR)*1.D-6)
         GVR(NR,NGR,13)  = GUCLIP((AJNB(NR)+AJRF(NR))*1.D-6)
         GVR(NR,NGR,14)  = GUCLIP(AJBS(NR)*1.D-6)
         GVR(NR,NGR,15)  = GUCLIP((PIN(NR,1)+PIN(NR,2)
     &                           +PIN(NR,3)+PIN(NR,4))*1.D-6)
         GVR(NR,NGR,16)  = GUCLIP(POH(NR)*1.D-6)
         GVR(NR,NGR,17)  = GUCLIP(VGR1(NR,2))
         GVR(NR,NGR,18)  = GUCLIP(VGR1(NR,1))
         GVR(NR,NGR,19)  = GUCLIP(VGR1(NR,3))
C         GVR(NR,NGR,19)  = GUCLIP(VGR3(NR,1))
         GVR(NR,NGR,20)  = GUCLIP(AK(NR,1))
         GVR(NR,NGR,21)  = GUCLIP(AK(NR,2))
C         GVR(NR,NGR,17)  = GUCLIP(RW(NR,1)*1.D-6*1.5D0)
C         GVR(NR,NGR,18)  = GUCLIP(RW(NR,2)*1.D-6*1.5D0)
C         GVR(NR,NGR,19)  = GUCLIP(PNB(NR)*1.D-6)
C         GVR(NR,NGR,20)  = GUCLIP(PNF(NR)*1.D-6)
         GVR(NR,NGR,22)  = GUCLIP(BP(NR))
         GVR(NR,NGR,23)  = GUCLIP(RPSI(NR))
      ENDDO
         GVR(1,NGR, 9)  = GUCLIP(Q0)
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
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
      INCLUDE 'trcomm.inc'
C
      CHARACTER KID*1
      CHARACTER K1*3,K2*3,K3*3,K4*3,K5*3,K6*3
      CHARACTER KCOM*40
C
      IF(KID.EQ.'N') THEN
         CALL TRNLIN(-29,IST,IERR)
      ELSEIF(KID.EQ.'1') THEN
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
     &                PBINT,PFINT,AJ(1)*1.D-6,
     &                PBCLT(1),PBCLT(2),PBCLT(3),PBCLT(4),
     &                PFCLT(1),PFCLT(2),PFCLT(3),PFCLT(4),
     &                POUT,PRLT,PCXT,PIET,
     &                PLT(1),PLT(2),PLT(3),PLT(4)
  604    FORMAT(' ',3X,'PINT  =',1PD10.3,'  POHT  =',1PD10.3,
     &               '  PNBT  =',1PD10.3,'  PNFTE =',1PD10.3/
     &          ' ',3X,'PRFTE =',1PD10.3,'  PRFTD =',1PD10.3,
     &               '  PRFTT =',1PD10.3,'  PRFTA =',1PD10.3/
     &          ' ',3X,'PBIN  =',1PD10.3,'  PFIN  =',1PD10.3,
     &               '  AJ0   =',1PD10.3/
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
C         OPEN(16,ACCESS='APPEND',FILE=KFNLOG)
         OPEN(16,ACCESS='SEQUENTIAL',FILE=KFNLOG)
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
      INCLUDE 'trcomm.inc'
C
    1 WRITE(6,*) '## INPUT MODE : 1:GVT(NT)  2:GVR(NR)  3:GVR(NG)'
      WRITE(6,*) '                NOW NGR=',NGR
      READ(5,*,END=9000,ERR=1) NID
      IF(NID.EQ.0) GOTO 9000
C
      IF(NID.EQ.1) THEN
   10    WRITE(6,*) '## INPUT NID,NTMIN,NTMAX,NTSTEP'
         READ(5,*,END=1,ERR=10) MID,MTMIN,MTMAX,MTSTEP
         IF(MID.EQ.0) GOTO 1
         DO NT=MTMIN,MTMAX,MTSTEP
            WRITE(6,601) NT,GT(NT),GVT(NT,MID)
         ENDDO
         GOTO 10
      ELSEIF(NID.EQ.2) THEN
   20    WRITE(6,*) '## INPUT NID,NG,NRMIN,NRMAX,NRSTEP,G or H(1 or 2)'
         READ(5,*,END=1,ERR=20) MID,NG,MRMIN,MRMAX,MRSTEP,MGH
         IF(MID.EQ.0) GOTO 1
         DO NR=MRMIN,MRMAX,MRSTEP
            IF(MGH.EQ.1) THEN
               WRITE(6,602) NG,NR,GRG(NR+1),GVR(NR,NG,MID)
            ELSEIF(MGH.EQ.2) THEN
               WRITE(6,602) NG,NR,GRM(NR),GVR(NR,NG,MID)
            ELSE
               GOTO 1
            ENDIF
         ENDDO
         GOTO 20
      ELSEIF(NID.EQ.3) THEN
   30    WRITE(6,*) '## INPUT NID,NR,NGMIN,NGMAX,NGSTEP,G or H(1 or 2)'
         READ(5,*,END=1,ERR=30) MID,NR,MGMIN,MGMAX,MGSTEP,MGH
         IF(MID.EQ.0) GOTO 1
         DO NG=MGMIN,MGMAX,MGSTEP
            IF(MGH.EQ.1) THEN
               WRITE(6,602) NG,NR,GRG(NR+1),GVR(NR,NG,MID)
            ELSEIF(MGH.EQ.2) THEN
               WRITE(6,602) NG,NR,GRM(NR),GVR(NR,NG,MID)
            ELSE
               GOTO 1
            ENDIF
         ENDDO
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
      INCLUDE 'trcomm.inc'
C
      CHARACTER STR*7
C
      GVMAX=GVT(1,N)
      GVMIN=GVT(1,N)
C
      DO NT=2,NGT
         GVMAX=MAX(GVMAX,GVT(NT,N))
         GVMIN=MIN(GVMIN,GVT(NT,N))
      ENDDO
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
      INCLUDE 'trcomm.inc'
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
      INCLUDE 'trcomm.inc'
C      INCLUDE 'trxcom.f'
C
      COMMON /TRXDT1/ KXNDEV,KXNDCG,KXNID
      COMMON /TRKID2/ KDIRW1,KDIRW2
      CHARACTER KXNDEV*80,KXNDCG*80,KXNID*80
      CHARACTER KDIRW*80,KDIRW1*80,KDIRW2*80,KFID*80
C
      KXNDEV='X'
      KXNDCG='test'
      KXNID ='in'
C
      CALL KTRIM(KXNDEV,IKNDEV)
      CALL KTRIM(KXNDCG,IKNDCG)
      CALL KTRIM(KXNID ,IKNID )
      KDIRW='../../tr.new/data/'//KXNDEV(1:IKNDEV)//'/'
     &                          //KXNDCG(1:IKNDCG)//'/'
     &                          //KXNID (1:IKNID )//'/'
C      KDIRW='../../../profile/profile_data/'//KXNDEV(1:IKNDEV)//'/'
C     &                          //KXNDCG(1:IKNDCG)//'/'
C     &                          //KXNID (1:IKNID )//'/'
      CALL KTRIM(KDIRW,IKDIRW)
      KDIRW1=KDIRW(1:IKDIRW)//KXNDEV(1:IKNDEV)
     &       //'1d'//KXNDCG(1:IKNDCG)//'.'
      KDIRW2=KDIRW(1:IKDIRW)//KXNDEV(1:IKNDEV)
     &       //'2d'//KXNDCG(1:IKNDCG)//'.'
C
C     *** 1D DATA ***
C
      KFID='IP'
      CALL TR_UFILE1D_CREATE(KFID,34,1.D6 ,IERR)
C
      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 2,1.D0 ,IERR)
C
      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 3,1.D0 ,IERR)
C
      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 4,1.D0 ,IERR)
C
      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 5,1.D0 ,IERR)
C
      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 6,1.D0 ,IERR)
C
      IF(MDLUF.NE.0) THEN
      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 8,1.D0 ,IERR)
      ELSE
      KFID='PNBI'
      CALL TR_UFILE1D_CREATE(KFID,41,1.D0 ,IERR)
      ENDIF
C
      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 9,1.D0 ,IERR)
C
      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID,10,1.D0 ,IERR)
C
      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID,11,1.D0 ,IERR)
C
      KFID='PRAD'
      CALL TR_UFILE1D_CREATE(KFID,60,1.D0 ,IERR)
C
      KFID='ZEFF'
      CALL TR_UFILE1D_CREATE(KFID,86,1.D0 ,IERR)
C
      KFID='VSURF'
      CALL TR_UFILE1D_CREATE(KFID,74,1.D0 ,IERR)
C
      KFID='LI'
      CALL TR_UFILE1D_CREATE(KFID,75,1.D0 ,IERR)
C
      KFID='WTH'
      CALL TR_UFILE1D_CREATE(KFID,31,1.D6 ,IERR)
C
      KFID='WTOT'
      CALL TR_UFILE1D_CREATE(KFID,33,1.D6 ,IERR)
C
      KFID='TE0'
      CALL TR_UFILE1D_CREATE(KFID, 9,1.D3 ,IERR)
C
      KFID='TI0'
      CALL TR_UFILE1D_CREATE(KFID,10,1.D3 ,IERR)
C
      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID,26,1.D0 ,IERR)
C
      KFID='POHM'
      CALL TR_UFILE1D_CREATE(KFID,40,1.D0 ,IERR)
C
      KFID='IBOOT'
      CALL TR_UFILE1D_CREATE(KFID,38,1.D6 ,IERR)
C
      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID,29,1.D0 ,IERR)
C
      KFID='PFUSION'
      CALL TR_UFILE1D_CREATE(KFID,46,1.D0 ,IERR)
C
C     *** 2D DATA ***
C
      KFID='TE'
      CALL TR_UFILE2D_CREATE(KFID, 1,1.D3 ,0,IERR)
C
      KFID='TI'
      CALL TR_UFILE2D_CREATE(KFID, 2,1.D3 ,0,IERR)
C
      KFID='NE'
      CALL TR_UFILE2D_CREATE(KFID, 5,1.D20,0,IERR)
C
      IF(MDLUF.NE.0) THEN
      KFID='QNBIE'
      CALL TR_UFILE2D_CREATE(KFID,89,1.D0 ,0,IERR)
      ENDIF
C
      KFID='QICRHE'
      CALL TR_UFILE2D_CREATE(KFID,38,1.D0 ,0,IERR)
C
      KFID='QECHE'
      CALL TR_UFILE2D_CREATE(KFID,36,1.D0 ,0,IERR)
C
      KFID='QLHE'
      CALL TR_UFILE2D_CREATE(KFID,37,1.D0 ,0,IERR)
C
      IF(MDLUF.NE.0) THEN
      KFID='QNBII'
      CALL TR_UFILE2D_CREATE(KFID,90,1.D0 ,0,IERR)
      ENDIF
C
      KFID='QICRHI'
      CALL TR_UFILE2D_CREATE(KFID,41,1.D0 ,0,IERR)
C
      KFID='QECHI'
      CALL TR_UFILE2D_CREATE(KFID,39,1.D0 ,0,IERR)
C
      KFID='QLHI'
      CALL TR_UFILE2D_CREATE(KFID,40,1.D0 ,0,IERR)
C
      KFID='CURNBI'
      CALL TR_UFILE2D_CREATE(KFID,11,1.D0 ,0,IERR)
C
      KFID='CURICRH'
      CALL TR_UFILE2D_CREATE(KFID,44,1.D0 ,0,IERR)
C
      KFID='CURECH'
      CALL TR_UFILE2D_CREATE(KFID,42,1.D0 ,0,IERR)
C
      KFID='CURLH'
      CALL TR_UFILE2D_CREATE(KFID,43,1.D0 ,0,IERR)
C
      KFID='NFAST'
      CALL TR_UFILE2D_CREATE(KFID,45,1.D20,0,IERR)
C
      KFID='QRAD'
      CALL TR_UFILE2D_CREATE(KFID,22,1.D0 ,0,IERR)
C
      KFID='ZEFFR'
      CALL TR_UFILE2D_CREATE(KFID,33,1.D0 ,0,IERR)
C
      KFID='Q'
      CALL TR_UFILE2D_CREATE(KFID,27,1.D0 ,1,IERR)
C
      KFID='CHIE'
      CALL TR_UFILE2D_CREATE(KFID,34,1.D0 ,1,IERR)
C
      KFID='CHII'
      CALL TR_UFILE2D_CREATE(KFID,35,1.D0 ,1,IERR)
C
      KFID='NM1'
      CALL TR_UFILE2D_CREATE(KFID, 6,1.D20,0,IERR)
C
      KFID='CURTOT'
      CALL TR_UFILE2D_CREATE(KFID, 9,1.D0 ,0,IERR)
C
      KFID='NIMP'
      CALL TR_UFILE2D_CREATE(KFID,46,1.D20,0,IERR)
C
      KFID='QOHM'
      CALL TR_UFILE2D_CREATE(KFID,15,1.D0 ,0,IERR)
C
      KFID='BPOL'
      CALL TR_UFILE2D_CREATE(KFID,47,1.D0 ,1,IERR)
C
      KFID='RMAJOR'
      CALL TR_UFILE2D_CREATE(KFID,49,1.D0 ,0,IERR)
C
      KFID='RMINOR'
      CALL TR_UFILE2D_CREATE(KFID,50,1.D0 ,0,IERR)
C
      KFID='VOLUME'
      CALL TR_UFILE2D_CREATE(KFID,51,1.D0 ,0,IERR)
C
      KFID='KAPPAR'
      CALL TR_UFILE2D_CREATE(KFID,52,1.D0 ,0,IERR)
C
      KFID='DELTAR'
      CALL TR_UFILE2D_CREATE(KFID,53,1.D0 ,0,IERR)
C
      KFID='GRHO1'
      CALL TR_UFILE2D_CREATE(KFID,54,1.D0 ,0,IERR)
C
      KFID='GRHO2'
      CALL TR_UFILE2D_CREATE(KFID,55,1.D0 ,0,IERR)
C
      KFID='CURBS'
      CALL TR_UFILE2D_CREATE(KFID,13,1.D0 ,0,IERR)
C
      KFID='CHITBE'
      CALL TR_UFILE2D_CREATE(KFID,56,1.D0 ,1,IERR)
C
      KFID='CHITBI'
      CALL TR_UFILE2D_CREATE(KFID,57,1.D0 ,1,IERR)
C
      KFID='ETAR'
      CALL TR_UFILE2D_CREATE(KFID,32,1.D0 ,0,IERR)
C
      RETURN
      END
C
C     *****
C
      SUBROUTINE TR_UFILE1D_CREATE(KFID,NUM,AMP,IERR)
C
      INCLUDE 'trcomm.inc'
      DIMENSION GTL(NTM),GF1(NTM),TF(NTM),F1(NRM)
      DIMENSION DGT(NTM),DIN(NTM),DERIV(NTM),UOUT(4,NTM)
      DIMENSION DERIVQ(NRM),UQ95(4,NRM)
      CHARACTER KFID*80,KERRF*80
C
      IF(KFID.EQ.'DIRECT') THEN
         IF(NUM.EQ.2) THEN
            KFID='BT'
            DO NTL=1,NGT
               TF(NTL)=BB
            ENDDO
         ELSEIF(NUM.EQ.3) THEN
            KFID='AMIN'
            DO NTL=1,NGT
               TF(NTL)=RA
            ENDDO
         ELSEIF(NUM.EQ.4) THEN
            KFID='RGEO'
            DO NTL=1,NGT
               TF(NTL)=RR
            ENDDO
         ELSEIF(NUM.EQ.5) THEN
            KFID='KAPPA'
            DO NTL=1,NGT
               TF(NTL)=RKAP
            ENDDO
         ELSEIF(NUM.EQ.6) THEN
            KFID='DELTA'
            DO NTL=1,NGT
               TF(NTL)=0.D0
            ENDDO
         ELSEIF(NUM.EQ.8) THEN
            KFID='PNBI'
            DO NTL=1,NGT
               TF(NTL)=DBLE(GVT(NTL,89)+GVT(NTL,90))
            ENDDO
         ELSEIF(NUM.EQ.9) THEN
            KFID='PECH'
            DO NTL=1,NGT
               TF(NTL)=DBLE(GVT(NTL,91)+GVT(NTL,92))
            ENDDO
         ELSEIF(NUM.EQ.10) THEN
            KFID='PICRH'
            DO NTL=1,NGT
               TF(NTL)=DBLE(GVT(NTL,95)+GVT(NTL,96))
            ENDDO
         ELSEIF(NUM.EQ.11) THEN
            KFID='PLH'
            DO NTL=1,NGT
               TF(NTL)=DBLE(GVT(NTL,93)+GVT(NTL,94))
            ENDDO
         ELSEIF(NUM.EQ.26) THEN
            KFID='Q95'
            ID=0
            IF(MDLUF.NE.0.AND.NRMAX.NE.NROMAX) THEN
               ID=1
               NRMAX=NROMAX
            ENDIF
            DO NTL=1,NGT
               DO NRL=1,NRMAX
                  F1(NRL)=DBLE(G3D(NRL,NTL,27))
               ENDDO
               CALL SPL1D (RG,F1,DERIVQ,UQ95,NRMAX,0,IERR)
               CALL SPL1DF(0.95D0,FQ95,RG,UQ95,NRMAX,IERR)
               TF(NTL)=FQ95
            ENDDO
            IF(ID.NE.0) NRMAX=NRAMAX
         ELSEIF(NUM.EQ.29) THEN
            KFID='RHOA'
            DO NTL=1,NGT
               TF(NTL)=RHOA
            ENDDO
         ENDIF
C
         DO NTL=1,NGT
            DGT(NTL)=DBLE(GT(NTL))
            DIN(NTL)=DBLE(TF(NTL))
            DERIV(NTL)=0.D0
         ENDDO
      ELSE
         DO NTL=1,NGT
            DGT(NTL)=DBLE(GT(NTL))
            DIN(NTL)=DBLE(GVT(NTL,NUM))
            DERIV(NTL)=0.D0
         ENDDO
      ENDIF
C
      CALL SPL1D(DGT,DIN,DERIV,UOUT,NGT,0,IERR)
      IF(IERR.NE.0) THEN
         IF(KFID.EQ.'DIRECT') THEN
            WRITE(6,'(A,I2,A,I2)')
     &           'XX TRXOUT: SPL1D DIRECT(',NUM,'): IERR=',IERR
         ELSE
            WRITE(6,'(A,I2,A,I2)')
     &           'XX TRXOUT: SPL1D GVT(',NUM,'): IERR=',IERR
         ENDIF
      ENDIF
C
      DTL=0.05D0
      NTLMAX=INT((GT(NGT)-GT(1))/SNGL(DTL))+1
C
      DO NTL=1,NTLMAX
         TIN=DBLE(GT(1))+DTL*DBLE(NTL-1)
         CALL SPL1DF(TIN,DATOUT,DGT,UOUT,NGT,IERR)
         WRITE(KERRF,'(A,I2,A,I2)')
     &           'XX TRXOUT: SPL1DF GVT(',NUM,'): IERR=',IERR
         IF(IERR.NE.0) WRITE(6,*) KERRF
         GTL(NTL)=GUCLIP(TIN)
         GF1(NTL)=GUCLIP(DATOUT*AMP)
      ENDDO
C
      CALL TRXW1D(KFID,GTL,GF1,NTM,NTLMAX)
C
      RETURN
      END
C
C     *****
C
      SUBROUTINE TR_UFILE2D_CREATE(KFID,NUM,AMP,ID,IERR)
C
      INCLUDE 'trcomm.inc'
      DIMENSION GF2(NRMP,NTM),GRL(NRMP),GTL(NTM)
      DIMENSION DGT(NTM),DIN(NTM)
      DIMENSION DERIV(NTM),U(4,NTM)
      CHARACTER KFID*80
C
      DO NTL=1,NGT
         DGT(NTL)=DBLE(GT(NTL))
      ENDDO
C
      NRLMAX=NRMAX
      DTL=0.05D0
      NTLMAX=INT((GT(NGT)-GT(1))/SNGL(DTL))+1
      IF(ID.EQ.0) THEN
         DO NRL=1,NRLMAX
            GRL(NRL)=GRM(NRL)
            DO NTL=1,NGT
               DIN(NTL)=DBLE(G3D(NRL,NTL,NUM))
            ENDDO
            CALL SPL1D(DGT,DIN,DERIV,U,NGT,0,IERR)
            IF(IERR.NE.0) WRITE(6,'(A,I2,A,I2)')
     &           'XX TRXOUT: SPL1D G3D(',NUM,'): IERR=',IERR
C     
            DO NTL=1,NTLMAX
               TIN=DBLE(GT(1))+DTL*DBLE(NTL-1)
               CALL SPL1DF(TIN,F0,DGT,U,NGT,IERR)
               IF(IERR.NE.0) WRITE(*,'(A,I2,A,I2)')
     &              'XX TRXOUT: SPL1DF G3D(',NUM,'): IERR=',IERR
               GTL(NTL)    =GUCLIP(TIN)
               GF2(NRL,NTL)=GUCLIP(F0*AMP)
            ENDDO
         ENDDO
      ELSEIF(ID.EQ.1) THEN
         NRLMAX=NRMAX+1
         NRL=1
            GRL(NRL)=GRG(NRL)
            IF(KFID.EQ.'Q') THEN
               DO NTL=1,NGT
                  DIN(NTL)=(4.D0*DBLE(G3D(NRL  ,NTL,NUM))
     &                          -DBLE(G3D(NRL+1,NTL,NUM)))/3.D0
               ENDDO
            ELSEIF(KFID.EQ.'BPOL') THEN
               DO NTL=1,NGT
                  DIN(NTL)=0.D0
               ENDDO
            ELSE
               DO NTL=1,NGT
                  R1=DBLE(GRL(NRL))
                  R2=DBLE(GRL(NRL+1))
                  F1=DBLE(G3D(NRL  ,NTL,NUM))
                  F2=DBLE(G3D(NRL+1,NTL,NUM))
                  DIN(NTL)=FCTR(R1,R2,F1,F2)
               ENDDO
            ENDIF
            CALL SPL1D(DGT,DIN,DERIV,U,NGT,0,IERR)
            IF(IERR.NE.0) WRITE(6,'(A,I2,A,I2)')
     &           'XX TRXOUT: SPL1D G3D(',NUM,'): IERR=',IERR
            DO NTL=1,NTLMAX
               TIN=DBLE(GT(1))+DTL*DBLE(NTL-1)
               CALL SPL1DF(TIN,F0,DGT,U,NGT,IERR)
               IF(IERR.NE.0) WRITE(*,'(A,I2,A,I2)')
     &              'XX TRXOUT: SPL1DF G3D(',NUM,'): IERR=',IERR
               GTL(NTL)    =GUCLIP(TIN)
               GF2(NRL,NTL)=GUCLIP(F0*AMP)
            ENDDO
C
         DO NRL=2,NRLMAX
            GRL(NRL)=GRG(NRL)
            DO NTL=1,NGT
               DIN(NTL)=DBLE(G3D(NRL-1,NTL,NUM))
            ENDDO
            CALL SPL1D(DGT,DIN,DERIV,U,NGT,0,IERR)
            IF(IERR.NE.0) WRITE(6,'(A,I2,A,I2)')
     &           'XX TRXOUT: SPL1D G3D(',NUM,'): IERR=',IERR
            DO NTL=1,NTLMAX
               TIN=DBLE(GT(1))+DTL*DBLE(NTL-1)
               CALL SPL1DF(TIN,F0,DGT,U,NGT,IERR)
               IF(IERR.NE.0) WRITE(*,'(A,I2,A,I2)')
     &              'XX TRXOUT: SPL1DF G3D(',NUM,'): IERR=',IERR
               GTL(NTL)    =GUCLIP(TIN)
               GF2(NRL,NTL)=GUCLIP(F0*AMP)
            ENDDO
         ENDDO
      ENDIF
      CALL TRXW2D(KFID,GTL,GRL,GF2,NRMP,NTM,NRLMAX,NTLMAX)
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
      CHARACTER KXNDEV*80,KXNDCG*80,KXNID*80,CDATE*9
      CHARACTER KDIRW1*80,KDIRW2*80,KFID*80,KFILE*80
C
      CALL GET_DATE(CDATE)
C
      CALL KTRIM(KDIRW1,KL1)
      KFILE=KDIRW1(1:KL1)//KFID
      WRITE(6,*) '- OPEN FILE:',KFILE(1:55)
C
      OPEN(16,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',ERR=10)
C
      WRITE(16,'(1X,A8,A8,A14,A29,A9)') KXNDCG(1:8),KXNDEV(1:8),
     &     '               ',
     &     ';-SHOT #- F(X) DATA -UF1DWR- ',CDATE
      WRITE(16,'(1X,A10,A20,A38)') 'TR:/tasktr','                    ',
     &     ';-SHOT DATE-  UFILES ASCII FILE SYSTEM'
      WRITE(16,'(1X,A30,A29)')
     &     'TIME          SECONDS         ',
     &     ';-INDEPENDENT VARIABLE LABEL-'
      WRITE(16,'(1X,A30,A27)') KFID,
     &     ';-DEPENDENT VARIABLE LABEL-'
      WRITE(16,'(1X,I1,A29,A39)') 2,'                             ',
     &     ';-PROC CODE- 0:RAW 1:AVG 2:SM. 3:AVG+SM'
      WRITE(16,'(1X,I11,A19,A33)') NTXMAX,'                   ',
     &     ';-# OF PTS-  X, F(X) DATA FOLLOW:'
C
      WRITE(16,'(1X,1P6E13.6)') (GT(NTX),NTX=1,NTXMAX)
      WRITE(16,'(1X,1P6E13.6)') (GF(NTX),NTX=1,NTXMAX)
C
      WRITE(16,'(A52)')
     &     ';----END-OF-DATA-----------------COMMENTS:-----------'
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
      CHARACTER KXNDEV*80,KXNDCG*80,KXNID*80,CDATE*9
      CHARACTER KDIRW1*80,KDIRW2*80,KFID*80,KFILE*80
C
      CALL GET_DATE(CDATE)
C
      CALL KTRIM(KDIRW2,KL2)
      KFILE=KDIRW2(1:KL2)//KFID
      WRITE(6,*) '- OPEN FILE:',KFILE(1:55)
C
      OPEN(16,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',ERR=10)
C
      WRITE(16,'(1X,A8,A8,A14,A29,A9)') KXNDCG(1:8),KXNDEV(1:8),
     &     '               ',
     &     ';-SHOT #- F(X) DATA -UF1DWR- ',CDATE
      WRITE(16,'(1X,A10,A20,A38)') 'TR:/tasktr','                    ',
     &     ';-SHOT DATE-  UFILES ASCII FILE SYSTEM'
      WRITE(16,'(1X,A30,A32)')
     &     'RHO                           ',
     &     ';-INDEPENDENT VARIABLE LABEL: X-'
      WRITE(16,'(1X,A30,A32)')
     &     'TIME          SECONDS         ',
     &     ';-INDEPENDENT VARIABLE LABEL: Y-'
      WRITE(16,'(1X,A30,A27)') KFID,
     &     ';-DEPENDENT VARIABLE LABEL-'
      WRITE(16,'(1X,I1,A29,A39)') 2,'                             ',
     &     ';-PROC CODE- 0:RAW 1:AVG 2:SM. 3:AVG+SM'
      WRITE(16,'(1X,I11,A19,A12)') NRXMAX,'                   ',
     &     ';-# OF X PTS-:'
      WRITE(16,'(1X,I11,A19,A35)') NTXMAX,'                   ',
     &     ';-# OF Y PTS-  X, F(X) DATA FOLLOW:'
C
      WRITE(16,'(1X,1P6E13.6)') (GR(NRX),NRX=1,NRXMAX)
      WRITE(16,'(1X,1P6E13.6)') (GT(NTX),NTX=1,NTXMAX)
      WRITE(16,'(1X,1P6E13.6)') ((GF(NRX,NTX),NRX=1,NRXMAX),
     &                                        NTX=1,NTXMAX)
C
      WRITE(16,'(A52)')
     &     ';----END-OF-DATA-----------------COMMENTS:-----------'
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
      SUBROUTINE GET_DATE(CDATE)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
      CHARACTER CDD*2,CDM*3,CDY*2,CDATE*9,CDATA(12)*3
      DATA (CDATA(I),I=1,12)
     &     /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
     &      'Oct','Nov','Dec'/
C
      CALL GUDATE(NDY,NDM,NDD,NTIH,NTIM,NTIS)
      CDM=CDATA(NDM)
      IF(NDD.LT.10) THEN
         WRITE(CDD,'(A1,I1)') ' ',NDD
      ELSE
         WRITE(CDD,'(I2)') NDD
      ENDIF
      NDY=NDY-100
      IF(NDY.LT.10) THEN
         WRITE(CDY,'(I1,I1)') 0,NDY
      ELSE
         WRITE(CDY,'(I2)') NDY
      ENDIF
      WRITE(CDATE,'(A2,A1,A3,A1,A2)') CDD,'-',CDM,'-',CDY
C
      RETURN
      END
