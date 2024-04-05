!     ***********************************************************

!           CALCULATE GLOBAL QUANTITIES

!     ***********************************************************

      SUBROUTINE TRGLOB

      USE trcomm
      USE trexec
      USE libitp
      IMPLICIT NONE
      INTEGER:: NEQ, NF, NMK, NR, NRL, NS, NSSN, NSSN1, NSVN, NSVN1, NSW, NW
      INTEGER:: NNB
      REAL(rkind)   :: ANFSUM, C83, DRH, DV53, PAI, PLST, RNSUM, RNTSUM, &
           & RTSUM, RWSUM, SLST, SUMM, SUML, SUMP, VOL, WPOL
      REAL(rkind),DIMENSION(NRMAX):: DSRHO

      IF(RHOA.NE.1.D0) NRMAX=NROMAX

!     *** Local beta ***
!        BETAL : toroidal beta
!        BETA  : volume-averaged toroidal beta
!        BETAPL: poloidal beta
!        BETAP : volume-averaged poloidal beta
!        BETAQ : toroidal beta for reaction rate
!               (ref. TOKAMAKS 3rd, p115)

      SUMM=0.D0
      SUMP=0.D0
      VOL =0.D0
      DO NR=1,NRMAX-1
         SUML = (SUM(RN(NR,1:NSM)*RT(NR,1:NSM))  &
               + SUM(RW(NR,1:NFMAX)))*RKEV*1.D20

!!         DSRHO(NR)=DVRHO(NR)/(2.D0*PI*RMJRHO(NR))
         DSRHO(NR)=DVRHO(NR)/(2.D0*PI*RR)
         VOL  = VOL  +         DVRHO(NR)*DR
         SUMM = SUMM + SUML   *DVRHO(NR)*DR
         SUMP = SUMP + SUML**2*DVRHO(NR)*DR
         BETA(NR)   = 2.D0*RMU0*SUMM     /(     VOL *BB**2)
         BETAL(NR)  = 2.D0*RMU0*SUML     /(          BB**2)
         BETAP(NR)  = 2.D0*RMU0*SUMM     /(     VOL *BP(NRMAX)**2)
         BETAPL(NR) = 2.D0*RMU0*SUML     /(          BP(NRMAX)**2)
         BETAQ(NR)  = 2.D0*RMU0*SQRT(SUMP)/(SQRT(VOL)*BB**2)
      ENDDO

      NR=NRMAX
      SUML = (SUM(RN(NR,1:NSM)*RT(NR,1:NSM)) &
            + SUM(RW(NR,1:NFMAX)))*RKEV*1.D20

!!         DSRHO(NR)=DVRHO(NR)/(2.D0*PI*RMJRHO(NR))
         DSRHO(NR)=DVRHO(NR)/(2.D0*PI*RR)
         VOL  = VOL  +         DVRHO(NR)*DR
         SUMM = SUMM + SUML   *DVRHO(NR)*DR
         SUMP = SUMP + SUML**2*DVRHO(NR)*DR
         BETA(NR)   = 2.D0*RMU0*SUMM     /(     VOL *BB**2)
         BETAL(NR)  = 2.D0*RMU0*SUML     /(          BB**2)
         BETAP(NR)  = 2.D0*RMU0*SUMM     /(     VOL *BP(NRMAX)**2)
         BETAPL(NR) = 2.D0*RMU0*SUML     /(          BP(NRMAX)**2)
         BETAQ(NR)  = 2.D0*RMU0*SQRT(SUMP)/(SQRT(VOL)*BB**2)

      BETA0 =(4.D0*BETA(1)  -BETA(2)  )/3.D0
      BETAP0=(4.D0*BETAPL(1)-BETAPL(2))/3.D0
      BETAQ0=(4.D0*BETAQ(1) -BETAQ(2) )/3.D0

!     *** Global beta ***
!        BETAPA: poloidal beta at separatrix
!        BETAA : toroidal beta at separatrix
!        BETAN : normalized toroidal beta (Troyon beta)

      BETAPA=BETAP(NRMAX)
      BETAA =BETA(NRMAX)
      BETAN =BETAA*1.D2/(RIP/(RA*BB))

!     *** Volume-averaged density and temperature ***
!     *** Central density and temperature         ***
!     *** Stored energy                           ***

!     +++ for electron and bulk ions +++
      DO NS=1,NSM
         RNSUM = SUM(RN(1:NRMAX,NS)*DVRHO(1:NRMAX))
         RTSUM = SUM(RN(1:NRMAX,NS)*RT(1:NRMAX,NS)*DVRHO(1:NRMAX))
         ANSAV(NS) = RNSUM*DR/VOL
         ANS0(NS) = FCTR(RM(1),RM(2),RN(1,NS),RN(2,NS))
         IF(RNSUM.GT.0.D0) THEN
            TSAV(NS) = RTSUM/RNSUM
         ELSE
            TSAV(NS) = 0.D0
         ENDIF
         TS0(NS) = FCTR(RM(1),RM(2),RT(1,NS),RT(2,NS))
         WST(NS) = 1.5D0*RTSUM*DR*RKEV*1.D14
      ENDDO

!     +++ for fast particles +++
      IF(MDLUF.NE.0) THEN
         ANFSUM = SUM(RNF(1:NRMAX,1)*DVRHO(1:NRMAX))
         RNTSUM = SUM(RNF(1:NRMAX,1)*RT(1:NRMAX,2)*DVRHO(1:NRMAX))
         RWSUM  = SUM(PBM(1:NRMAX)  *DVRHO(1:NRMAX))
         NF=1
            WFT(NF) = RWSUM*DR*1.D-6-1.5D0*RNTSUM*DR*RKEV*1.D14
            ANFAV(NF) = ANFSUM*DR/VOL
            ANF0(NF)  = FCTR(RM(1),RM(2),RNF(1,1),RNF(2,1))
            IF(ANFSUM.GT.0.D0) THEN
               TFAV(NF)  = RWSUM/(RKEV*1.D20)/ANFSUM
            ELSE
               TFAV(NF)  = 0.D0
            ENDIF
            IF(RNF(1,1).GT.0.D0) THEN
               TF0(NF)  = FCTR(RM(1),RM(2),PBM(1),PBM(2))/(RKEV*1.D20)/ANF0(NF)
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
         DO NF=1,NFMAX
            ANFSUM=0.D0
            RWSUM=0.D0
            DO NR=1,NRMAX
               ANFSUM = ANFSUM+RNF(NR,NF)*DVRHO(NR)
               RWSUM  = RWSUM +RW (NR,NF)*DVRHO(NR)
            ENDDO
            WFT(NF) = 1.5D0*RWSUM*DR*RKEV*1.D14
            ANFAV(NF) = ANFSUM*DR/VOL
            ANF0(NF)  = FCTR(RM(1),RM(2),RNF(1,NF),RNF(2,NF))
            IF(ANFSUM.GT.0.D0) THEN
               TFAV(NF)  = RWSUM/ANFSUM
            ELSE
               TFAV(NF)  = 0.D0
            ENDIF
            IF(RNF(1,NF).GT.0.D0) THEN
               TF0(NF)  = FCTR(RM(1),RM(2),RW(1,NF),RW(2,NF))/ANF0(NF)
            ELSE
               TF0(NF)  = 0.D0
            ENDIF
         ENDDO
      ENDIF

!     *** Line-averaged densities ***

      DO NS=1,NSM
         ANLAV(NS)=SUM(RN(1:NRMAX,NS))*DR
      ENDDO

!     *** Particle source ***

      DO NS=1,NSMAX
         SPSCT(NS) = SUM(SPSC(1:NRMAX,NS)*DVRHO(1:NRMAX))*DR
      END DO

!     *** Ohmic, NBI and fusion powers ***

      DO NR=1,NRMAX
         PNBNR(NR)=0.D0
         PNFNR(NR)=0.D0
         DO NS=1,NSMAX
            PNBNR(NR)=PNBNR(NR)+PNB(NS,NR)
            PNFNR(NR)=PNFNR(NR)+PNF(NS,NR)
         END DO
      END DO
            
      POHT = SUM(POH(1:NRMAX)*DVRHO(1:NRMAX))*DR/1.D6
      PNBT = SUM(PNBNR(1:NRMAX)*DVRHO(1:NRMAX))*DR/1.D6
      PNFT = SUM(PNFNR(1:NRMAX)*DVRHO(1:NRMAX))*DR/1.D6

!     *** External power typically for NBI from exp. data ***

      DO NS=1,NSM
         PEXT(NS) = SUM(PEX(1:NRMAX,NS)*DVRHO(1:NRMAX))*DR/1.D6
      ENDDO

!     *** RF power ***

      DO NS=1,NSM
         PRFVT(NS,1) = SUM(PRFV(1:NRMAX,NS,1)*DVRHO(1:NRMAX))*DR/1.D6
         PRFVT(NS,2) = SUM(PRFV(1:NRMAX,NS,2)*DVRHO(1:NRMAX))*DR/1.D6
         PRFVT(NS,3) = SUM(PRFV(1:NRMAX,NS,3)*DVRHO(1:NRMAX))*DR/1.D6
         PRFT (NS  ) = SUM(PRF (1:NRMAX,NS)  *DVRHO(1:NRMAX))*DR/1.D6
      ENDDO

!     *** Total NBI power distributed on electrons and bulk ions ***

      PBINT = SUM(PBIN(1:NRMAX)*DVRHO(1:NRMAX))*DR/1.D6
      DO NS=1,NSM
         PBCLT(NS) = SUM(PBCL(1:NRMAX,NS)*DVRHO(1:NRMAX))*DR/1.D6
      ENDDO

!     *** Total RF power distributed on electrons and bulk ions ***

      PFINT = SUM(PFIN(1:NRMAX)*DVRHO(1:NRMAX))*DR/1.D6
      DO NS=1,NSM
         PFCLT(NS) = SUM(PFCL(1:NRMAX,NS)*DVRHO(1:NRMAX))*DR/1.D6
      ENDDO

!     *** Radiation, charge exchange and ionization losses ***

      PRBT  = SUM(PRB(1:NRMAX)*DVRHO(1:NRMAX))*DR/1.D6
      PRCT  = SUM(PRC(1:NRMAX)*DVRHO(1:NRMAX))*DR/1.D6
      PRLT  = SUM(PRL(1:NRMAX)*DVRHO(1:NRMAX))*DR/1.D6
      PRSUMT= SUM(PRSUM(1:NRMAX)*DVRHO(1:NRMAX))*DR/1.D6
      PCXT  = SUM(PCX(1:NRMAX)*DVRHO(1:NRMAX))*DR/1.D6
      PIET  = SUM(PIE(1:NRMAX)*DVRHO(1:NRMAX))*DR/1.D6

!     *** Current densities ***

      AJNBT     = SUM(AJNB (1:NRMAX)  *DSRHO(1:NRMAX))*DR/1.D6
      AJRFVT(1) = SUM(AJRFV(1:NRMAX,1)*DSRHO(1:NRMAX))*DR/1.D6
      AJRFVT(2) = SUM(AJRFV(1:NRMAX,2)*DSRHO(1:NRMAX))*DR/1.D6
      AJRFVT(3) = SUM(AJRFV(1:NRMAX,3)*DSRHO(1:NRMAX))*DR/1.D6
      AJRFT     = SUM(AJRF (1:NRMAX)  *DSRHO(1:NRMAX))*DR/1.D6
      AJBST     = SUM(AJBS (1:NRMAX)  *DSRHO(1:NRMAX))*DR/1.D6

      AJT    = SUM(AJ   (1:NRMAX)*DSRHO(1:NRMAX))*DR/1.D6
      ! Plasma current, AJTTOR
      AJTTOR = SUM(AJTOR(1:NRMAX)*DSRHO(1:NRMAX))*DR/1.D6
      AJOHT  = SUM(AJOH (1:NRMAX)*DSRHO(1:NRMAX))*DR/1.D6

!     *** Output source and power ***
!        evaluate output power using flux at NR=NRMAX

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
      C83=8.D0/3.D0
      DO NEQ=1,NEQMAX
         NSSN=NSS(NEQ)
         NSVN=NSV(NEQ)
         SUML=0.D0
         DO NW=1,NEQMAX
            NSSN1=NSS(NW)
            NSVN1=NSV(NW)
            IF(NSVN1.EQ.1) THEN
               SUML=SUML+(-DD(NEQ,NW,NMK,NSW)/3.D0)*DRH *RN(NRL-1,NSSN1) &
     &                  +  DD(NEQ,NW,NMK,NSW)*3.D0 *DRH *RN(NRL  ,NSSN1) &
     &                  +( VV(NEQ,NW,NMK,NSW)      *DRH &
     &                    -DD(NEQ,NW,NMK,NSW)*C83 )*DRH *PNSS(NSSN1)
            ELSEIF(NSVN1.EQ.2) THEN
               SUML=SUML+(-DD(NEQ,NW,NMK,NSW)/3.D0)*DRH *RN(NRL-1,NSSN1)*RT(NRL-1,NSSN1) &
     &                  +  DD(NEQ,NW,NMK,NSW)*3.D0 *DRH *RN(NRL  ,NSSN1)*RT(NRL  ,NSSN1) &
     &                  +( VV(NEQ,NW,NMK,NSW)      *DRH &
     &                    -DD(NEQ,NW,NMK,NSW)*C83  *DRH)*PNSS(NSSN1)*PTS (NSSN1)
            ENDIF
         ENDDO
         IF(NSVN.EQ.1) THEN
            SLT(NSSN)=SUML*RKEV*1.D14
         ELSEIF(NSVN.EQ.2) THEN
            PLT(NSSN)=SUML*RKEV*1.D14
         ENDIF
      ENDDO

!     *** Ionization, fusion and NBI fuelling ***

      SIET = SUM(SIE(1:NRMAX)*DVRHO(1:NRMAX))*DR
      IF(NSMAX.GE.4) THEN
         SNFT=SUM(SNF(4,1:NRMAX)*DVRHO(1:NRMAX))*DR
      ELSE
         SNFT=0.D0
      END IF
      SNBT=0.D0
      DO NNB=1,NNBMAX
         SNBT = SNBT+SUM(SNB(NSPNB(NNB),1:NRMAX)*DVRHO(1:NRMAX))*DR
      END DO

!     *** Pellet injection fuelling ***

      DO NS=1,NSM
         SPET(NS) = SUM(SPE(1:NRMAX,NS)*DVRHO(1:NRMAX))*DR/RKAP
      ENDDO

!     *** Input and output sources and powers ***

      WBULKT=SUM(WST(1:NSM))
      PEXST =SUM(PEXT(1:NSM))
      PRFST =SUM(PRFT(1:NSM))
      PLST  =SUM(PLT(1:NSM))
      SLST  =SUM(SLT(1:NSM))
      WTAILT=SUM(WFT(1:NFMAX))

      WPT =WBULKT+WTAILT
      PINT=POHT+PNBT+PRFST+PNFT+PEXST
      POUT=PLST+PCXT+PIET+PRBT+PRCT+PRLT
      SINT=SIET+SNBT
      SOUT=SLST

!     *** Energy confinement times ***
!        TAUE1: steady state
!        TAUE2: transient

      IF(ABS(T-TPRE).LE.1.D-70) THEN
         WPDOT=0.D0
      ELSE
         WPDOT=(WPT-WPPRE)/(T-TPRE)
      ENDIF
      WPPRE=WPT
      TPRE=T

      TAUE1=WPT/PINT
      TAUE2=WPT/(PINT-WPDOT)

!     *** Inductance and one-turn voltage ***

      WPOL=SUM(ABVRHOG(1:NRMAX)*RDPVRHOG(1:NRMAX)**2*DVRHOG(1:NRMAX))*DR/(2.D0*RMU0)
      ALI=4.D0*WPOL/(RMU0*RR*(AJTTOR*1.D6)**2)

      VLOOP = EZOH(NRMAX)*2.D0*PI*RR

!     *** Confinement scalings (TOKAMAKS 3rd p183,184) ***
!        TAUE89: ITER89-P L-mode scaling
!        TAUE98: IPB98(y,2) H-mode scaling with ELMs
!        H98Y2: H-mode factor

!     volume-averaged isotopic mass number
      PAI = (PA(2)*ANSAV(2)+PA(3)*ANSAV(3)+PA(4)*ANSAV(4))  /(ANSAV(2)+ANSAV(3)+ANSAV(4))

      TAUE89=4.8D-2*(ABS(RIP)**0.85D0)    *(RR**1.2D0) *(RA**0.3D0)  *(RKAP**0.5D0) &
     &             *(ANLAV(1)**0.1D0)*(ABS(BB)**0.2D0) *(PAI**0.5D0) *(PINT**(-0.5D0))
      TAUE98=0.145D0*(ABS(RIP)**0.93D0)   *(RR**1.39D0)*(RA**0.58D0) *(RKAP**0.78D0) &
     &            *(ANLAV(1)**0.41D0)*(ABS(BB)**0.15D0)*(PAI**0.19D0)*(PINT**(-0.69D0))
      H98Y2=TAUE2/TAUE98

!     *** Fusion production rate ***

      QF=5.D0*PNFT/(POHT+PNBT+PRFST+PEXST)

!     *** Distance of q=1 surface from magnetic axis ***

      IF(Q0.GE.1.D0) THEN
         RQ1=0.D0
      ELSE
         IF(QP(1).GT.1.D0) THEN
            RQ1=SQRT((1.D0-Q0)/(QP(1)-Q0))*DR
            GOTO 310
         ENDIF
         DO NR=2,NRMAX
            IF(QP(NR).GT.1.D0) THEN
              RQ1=(RG(NR)-RG(NR-1))*RA*(1.D0-QP(NR-1))/(QP(NR)-QP(NR-1))+RG(NR-1)*RA
               GOTO 310
            ENDIF
         ENDDO
         RQ1=RA
  310    CONTINUE
      ENDIF

!     *** Effective charge number at axis ***

      ZEFF0=(4.D0*ZEFF(1)-ZEFF(2))/3.D0

      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
      RETURN
      END SUBROUTINE TRGLOB

!     ***********************************************************

!           SAVE DATA TO DATAT FOR GRAPHICS

!     ***********************************************************

      SUBROUTINE TRATOT

        USE TRCOMM
        USE libspl1d
      IMPLICIT NONE
      INTEGER:: IERR, NR
      REAL(rkind)   :: RMN, F0D
      REAL(rkind),DIMENSION(NRMAX):: DERIV, U0
      REAL(rkind),DIMENSION(4,NRMAX):: U
      REAL(rkind)   :: TRCOFS
      REAL   :: GUCLIP


      IF(NGT.GE.NTM) RETURN
      NGT=NGT+1

      GT    (NGT) = GUCLIP(T)
!
      GVT(NGT, 1) = GUCLIP(ANS0(1))
      GVT(NGT, 2) = GUCLIP(ANS0(2))
      GVT(NGT, 3) = GUCLIP(ANS0(3))
      GVT(NGT, 4) = GUCLIP(ANS0(4))
      GVT(NGT, 5) = GUCLIP(ANSAV(1))
      GVT(NGT, 6) = GUCLIP(ANSAV(2))
      GVT(NGT, 7) = GUCLIP(ANSAV(3))
      GVT(NGT, 8) = GUCLIP(ANSAV(4))

      GVT(NGT, 9) = GUCLIP(TS0(1))
      GVT(NGT,10) = GUCLIP(TS0(2))
      GVT(NGT,11) = GUCLIP(TS0(3))
      GVT(NGT,12) = GUCLIP(TS0(4))
      GVT(NGT,13) = GUCLIP(TSAV(1))
      GVT(NGT,14) = GUCLIP(TSAV(2))
      GVT(NGT,15) = GUCLIP(TSAV(3))
      GVT(NGT,16) = GUCLIP(TSAV(4))

      GVT(NGT,17) = GUCLIP(WST(1))
      GVT(NGT,18) = GUCLIP(WST(2))
      GVT(NGT,19) = GUCLIP(WST(3))
      GVT(NGT,20) = GUCLIP(WST(4))

      GVT(NGT,21) = GUCLIP(ANF0(1))
      GVT(NGT,22) = GUCLIP(ANF0(2))
      GVT(NGT,23) = GUCLIP(ANFAV(1))
      GVT(NGT,24) = GUCLIP(ANFAV(2))
      GVT(NGT,25) = GUCLIP(TF0(1))
      GVT(NGT,26) = GUCLIP(TF0(2))
      GVT(NGT,27) = GUCLIP(TFAV(1))
      GVT(NGT,28) = GUCLIP(TFAV(2))

      GVT(NGT,29) = GUCLIP(WFT(1))
      GVT(NGT,30) = GUCLIP(WFT(2))
      GVT(NGT,31) = GUCLIP(WBULKT)
      GVT(NGT,32) = GUCLIP(WTAILT)
      GVT(NGT,33) = GUCLIP(WPT)

      GVT(NGT,34) = GUCLIP(AJT)
      GVT(NGT,35) = GUCLIP(AJOHT)
      GVT(NGT,36) = GUCLIP(AJNBT)
      GVT(NGT,37) = GUCLIP(AJRFT)
      GVT(NGT,38) = GUCLIP(AJBST)

      GVT(NGT,39) = GUCLIP(PINT)
      GVT(NGT,40) = GUCLIP(POHT)
      GVT(NGT,41) = GUCLIP(PNBT)
      GVT(NGT,42) = GUCLIP(PRFT(1))
      GVT(NGT,43) = GUCLIP(PRFT(2))
      GVT(NGT,44) = GUCLIP(PRFT(3))
      GVT(NGT,45) = GUCLIP(PRFT(4))
      GVT(NGT,46) = GUCLIP(PNFT)

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

      GVT(NGT,57) = GUCLIP(POUT)
      GVT(NGT,58) = GUCLIP(PCXT)
      GVT(NGT,59) = GUCLIP(PIET)
      GVT(NGT,60) = GUCLIP(PRSUMT)
      GVT(NGT,61) = GUCLIP(PLT(1))
      GVT(NGT,62) = GUCLIP(PLT(2))
      GVT(NGT,63) = GUCLIP(PLT(3))
      GVT(NGT,64) = GUCLIP(PLT(4))

      GVT(NGT,65) = GUCLIP(SINT)
      GVT(NGT,66) = GUCLIP(SIET)
      GVT(NGT,67) = GUCLIP(SNBT)
      GVT(NGT,68) = GUCLIP(SNFT)
      GVT(NGT,69) = GUCLIP(SOUT)
      GVT(NGT,70) = GUCLIP(SLT(1))
      GVT(NGT,71) = GUCLIP(SLT(2))
      GVT(NGT,72) = GUCLIP(SLT(3))
      GVT(NGT,73) = GUCLIP(SLT(4))

      GVT(NGT,74) = GUCLIP(VLOOP)
      GVT(NGT,75) = GUCLIP(ALI)
      GVT(NGT,76) = GUCLIP(RQ1)
      GVT(NGT,77) = GUCLIP(Q0)

      GVT(NGT,78) = GUCLIP(WPDOT)
      GVT(NGT,79) = GUCLIP(TAUE1)
      GVT(NGT,80) = GUCLIP(TAUE2)
      GVT(NGT,81) = GUCLIP(TAUE89)

      GVT(NGT,82) = GUCLIP(BETAP0)
      GVT(NGT,83) = GUCLIP(BETAPA)
      GVT(NGT,84) = GUCLIP(BETA0)
      GVT(NGT,85) = GUCLIP(BETAA)

      GVT(NGT,86) = GUCLIP(ZEFF0)
      GVT(NGT,87) = GUCLIP(QF)
      GVT(NGT,88) = GUCLIP(RIP)
!
      GVT(NGT,89) = GUCLIP(PEXT(1))
      GVT(NGT,90) = GUCLIP(PEXT(2))
      GVT(NGT,91) = GUCLIP(PRFVT(1,1)) ! ECH  to electron
      GVT(NGT,92) = GUCLIP(PRFVT(2,1)) ! ECH  to ions
      GVT(NGT,93) = GUCLIP(PRFVT(1,2)) ! LH   to electron
      GVT(NGT,94) = GUCLIP(PRFVT(2,2)) ! LH   to ions
      GVT(NGT,95) = GUCLIP(PRFVT(1,3)) ! ICRH to electron
      GVT(NGT,96) = GUCLIP(PRFVT(2,3)) ! ICRH to ions

      GVT(NGT,97) = GUCLIP(RR)
      GVT(NGT,98) = GUCLIP(RA)
      GVT(NGT,99) = GUCLIP(BB)
      GVT(NGT,100)= GUCLIP(RKAP)
      GVT(NGT,101)= GUCLIP(AJTTOR)

      GVT(NGT,102)= GUCLIP(TAUE98)
      GVT(NGT,103)= GUCLIP(H98Y2)
      GVT(NGT,104)= GUCLIP(ANLAV(1))
      GVT(NGT,105)= GUCLIP(ANLAV(2))
      GVT(NGT,106)= GUCLIP(ANLAV(3))
      GVT(NGT,107)= GUCLIP(ANLAV(4))

      GVT(NGT,108)= GUCLIP(PRBT)
      GVT(NGT,109)= GUCLIP(PRCT)
      GVT(NGT,110)= GUCLIP(PRLT)

!     *** FOR 3D ***

      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      CALL SPL1D  (RM,DVRHO,DERIV,U,NRMAX,0,IERR)
      CALL SPL1DI0(RM,U,U0,NRMAX,IERR)
      DO NR=1,NRMAX
         GVRT(NR,NGT, 1) = GUCLIP(RT(NR,1))
         GVRT(NR,NGT, 2) = GUCLIP(RT(NR,2))
         GVRT(NR,NGT, 3) = GUCLIP(RT(NR,3))
         GVRT(NR,NGT, 4) = GUCLIP(RT(NR,4))

         GVRT(NR,NGT, 5) = GUCLIP(RN(NR,1))
         GVRT(NR,NGT, 6) = GUCLIP(RN(NR,2))
         GVRT(NR,NGT, 7) = GUCLIP(RN(NR,3))
         GVRT(NR,NGT, 8) = GUCLIP(RN(NR,4))

         GVRT(NR,NGT, 9) = GUCLIP(AJ  (NR))
         GVRT(NR,NGT,10) = GUCLIP(AJOH(NR))
         GVRT(NR,NGT,11) = GUCLIP(AJNB(NR))
         GVRT(NR,NGT,12) = GUCLIP(AJRF(NR))
         GVRT(NR,NGT,13) = GUCLIP(AJBS(NR))

         GVRT(NR,NGT,14) = GUCLIP(POH(NR)+PNBNR(NR)+PNFNR(NR) &
     &                         +PEX(NR,1)+PEX(NR,2)+PEX(NR,3)+PEX(NR,4) &
     &                         +PRF(NR,1)+PRF(NR,2)+PRF(NR,3)+PRF(NR,4))
         GVRT(NR,NGT,15) = GUCLIP(POH(NR))
         GVRT(NR,NGT,16) = GUCLIP(PNBNR(NR))
         GVRT(NR,NGT,17) = GUCLIP(PNFNR(NR))
         GVRT(NR,NGT,18) = GUCLIP(PRF(NR,1))
         GVRT(NR,NGT,19) = GUCLIP(PRF(NR,2))
         GVRT(NR,NGT,20) = GUCLIP(PRF(NR,3))
         GVRT(NR,NGT,21) = GUCLIP(PRF(NR,4))
         GVRT(NR,NGT,22) = GUCLIP(PRL(NR))
         GVRT(NR,NGT,23) = GUCLIP(PCX(NR))
         GVRT(NR,NGT,24) = GUCLIP(PIE(NR))
         GVRT(NR,NGT,25) = GUCLIP(PEX(NR,1))
         GVRT(NR,NGT,26) = GUCLIP(PEX(NR,2))

         GVRT(NR,NGT,27) = GUCLIP(QP(NR))
         GVRT(NR,NGT,28) = GUCLIP(EZOH(NR))
         GVRT(NR,NGT,29) = GUCLIP(BETA(NR))
         GVRT(NR,NGT,30) = GUCLIP(BETAP(NR))
         GVRT(NR,NGT,31) = GUCLIP(EZOH(NR)*2.D0*PI*RR)
         GVRT(NR,NGT,32) = GUCLIP(ETA(NR))
         GVRT(NR,NGT,33) = GUCLIP(ZEFF(NR))
         GVRT(NR,NGT,34) = GUCLIP(AK(NR,1))
         GVRT(NR,NGT,35) = GUCLIP(AK(NR,2))

         GVRT(NR,NGT,36) = GUCLIP(PRFV(NR,1,1))
         GVRT(NR,NGT,37) = GUCLIP(PRFV(NR,1,2))
         GVRT(NR,NGT,38) = GUCLIP(PRFV(NR,1,3))
         GVRT(NR,NGT,39) = GUCLIP(PRFV(NR,2,1))
         GVRT(NR,NGT,40) = GUCLIP(PRFV(NR,2,2))
         GVRT(NR,NGT,41) = GUCLIP(PRFV(NR,2,3))

         GVRT(NR,NGT,42) = GUCLIP(AJRFV(NR,1))
         GVRT(NR,NGT,43) = GUCLIP(AJRFV(NR,2))
         GVRT(NR,NGT,44) = GUCLIP(AJRFV(NR,3))

         GVRT(NR,NGT,45) = GUCLIP(RW(NR,1)+RW(NR,2))
         GVRT(NR,NGT,46) = GUCLIP(ANC(NR)+ANFE(NR))
         GVRT(NR,NGT,47) = GUCLIP(BP(NR))
         GVRT(NR,NGT,48) = GUCLIP(RPSI(NR))

         GVRT(NR,NGT,49) = GUCLIP(RMJRHO(NR))
         GVRT(NR,NGT,50) = GUCLIP(RMNRHO(NR))
         RMN=(DBLE(NR)-0.5D0)*DR
         CALL SPL1DI(RMN,F0D,RM,U,U0,NRMAX,IERR)
         GVRT(NR,NGT,51) = GUCLIP(F0D)
         GVRT(NR,NGT,52) = GUCLIP(RKPRHO(NR))
         GVRT(NR,NGT,53) = GUCLIP(1.D0) ! DELTAR
         GVRT(NR,NGT,54) = GUCLIP(AR1RHO(NR))
         GVRT(NR,NGT,55) = GUCLIP(AR2RHO(NR))
         GVRT(NR,NGT,56) = GUCLIP(AKDW(NR,1))
         GVRT(NR,NGT,57) = GUCLIP(AKDW(NR,2))
         GVRT(NR,NGT,58) = GUCLIP(RN(NR,1)*RT(NR,1))
         GVRT(NR,NGT,59) = GUCLIP(RN(NR,2)*RT(NR,2))

         GVRT(NR,NGT,60) = GUCLIP(VTOR(NR))
         GVRT(NR,NGT,61) = GUCLIP(VPOL(NR))

         GVRT(NR,NGT,62) = GUCLIP(S(NR)-ALPHA(NR))
         GVRT(NR,NGT,63) = GUCLIP(ER(NR))
         GVRT(NR,NGT,64) = GUCLIP(S(NR))
         GVRT(NR,NGT,65) = GUCLIP(ALPHA(NR))
         GVRT(NR,NGT,66) = GUCLIP(TRCOFS(S(NR),ALPHA(NR),RKCV(NR)))
         GVRT(NR,NGT,67) = GUCLIP(2.D0*PI/QP(NR))
         
         IF(NSMAX.GE.1) GVRT(NR,NGT,68) = GUCLIP(PNB(1,NR))
         IF(NSMAX.GE.2) GVRT(NR,NGT,69) = GUCLIP(PNB(2,NR))
         IF(NSMAX.GE.3) GVRT(NR,NGT,70) = GUCLIP(PNB(3,NR))
         IF(NSMAX.GE.4) GVRT(NR,NGT,71) = GUCLIP(PNB(4,NR))
         IF(NSMAX.GE.1) GVRT(NR,NGT,72) = GUCLIP(PNF(1,NR))
         IF(NSMAX.GE.2) GVRT(NR,NGT,73) = GUCLIP(PNF(2,NR))
         IF(NSMAX.GE.3) GVRT(NR,NGT,74) = GUCLIP(PNF(3,NR))
         IF(NSMAX.GE.4) GVRT(NR,NGT,75) = GUCLIP(PNF(4,NR))

      ENDDO
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
!
      RETURN
      END SUBROUTINE TRATOT

!     ***********************************************************

!           SAVE DATA TO DATA FOR GRAPHICS

!     ***********************************************************

      SUBROUTINE TRATOG

      USE TRCOMM, ONLY : AJ, AJBS, AJNB, AJOH, AJRF, AK, BP, EZOH, GTR, GVR, NGM, NGR, NRAMAX, NRMAX , NROMAX, &
     &                   PIN, POH, Q0, QP, RHOA, RN, RPSI, RT, T, VGR1
      IMPLICIT NONE
      INTEGER:: NR
      REAL   :: GUCLIP


      IF(NGR.GE.NGM) RETURN
      NGR=NGR+1
      GTR(NGR)=GUCLIP(T)

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
         GVR(NR,NGR,15)  = GUCLIP((PIN(NR,1)+PIN(NR,2) &
     &                           +PIN(NR,3)+PIN(NR,4))*1.D-6)
         GVR(NR,NGR,16)  = GUCLIP(POH(NR)*1.D-6)
         GVR(NR,NGR,17)  = GUCLIP(VGR1(NR,2))
         GVR(NR,NGR,18)  = GUCLIP(VGR1(NR,1))
         GVR(NR,NGR,19)  = GUCLIP(VGR1(NR,3))
!         GVR(NR,NGR,19)  = GUCLIP(VGR3(NR,1))
         GVR(NR,NGR,20)  = GUCLIP(AK(NR,1))
         GVR(NR,NGR,21)  = GUCLIP(AK(NR,2))
!         GVR(NR,NGR,17)  = GUCLIP(RW(NR,1)*1.D-6*1.5D0)
!         GVR(NR,NGR,18)  = GUCLIP(RW(NR,2)*1.D-6*1.5D0)
!         GVR(NR,NGR,19)  = GUCLIP(PNB(NR)*1.D-6)
!         GVR(NR,NGR,20)  = GUCLIP(PNF(NR)*1.D-6)
         GVR(NR,NGR,22)  = GUCLIP(BP(NR))
         GVR(NR,NGR,23)  = GUCLIP(RPSI(NR))
      ENDDO
         GVR(1,NGR, 9)  = GUCLIP(Q0)
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX

      RETURN
      END SUBROUTINE TRATOG

!     ***********************************************************

!           PRINT GLOBAL QUANTITIES

!     ***********************************************************

      SUBROUTINE TRPRNT(KID)

      USE trcomm
      USE trparm
      IMPLICIT NONE
      CHARACTER(LEN=1),INTENT(IN):: KID
      INTEGER:: I, IERR, IST, NDD, NDM, NDY, NTH1, NTM1, NTS1
      REAL   :: GTCPU2
      CHARACTER(LEN=3) :: K1, K2, K3, K4, K5, K6
      CHARACTER(LEN=40):: KCOM


      IF(KID.EQ.'N') THEN
         CALL tr_nlin(-29,IST,IERR)
      ELSEIF(KID.EQ.'1') THEN
         WRITE(6,601) T,WPT,WBULKT,WTAILT,WPDOT,TAUE1,TAUE2,TAUE89,TAUE98, &
     &                QF,BETAP0,BETAPA,BETA0,BETAA,Q0,RQ1,ZEFF0,BETAN
  601    FORMAT(' ','# TIME : ',F7.3,' SEC'/ &
     &          ' ',3X,'WPT   =',1PD10.3,'  WBULKT=',1PD10.3, &
     &               '  WTAILT=',1PD10.3,'  WPDOT =',1PD10.3/ &
     &          ' ',3X,'TAUE1 =',1PD10.3,'  TAUE2 =',1PD10.3, &
     &               '  TAUE89=',1PD10.3,'  TAUE98=',1PD10.3/ &
     &          ' ',3X,'QF    =',1PD10.3/ &
     &          ' ',3X,'BETAP0=',1PD10.3,'  BETAPA=',1PD10.3, &
     &               '  BETA0 =',1PD10.3,'  BETAA =',1PD10.3/ &
     &          ' ',3X,'Q0    =',1PD10.3,'  RQ1   =',1PD10.3, &
     &               '  ZEFF0 =',1PD10.3,'  BETAN =',1PD10.3)

         WRITE(6,602) WST(1),TS0(1),TSAV(1),ANSAV(1), &
     &                WST(2),TS0(2),TSAV(2),ANSAV(2), &
     &                WST(3),TS0(3),TSAV(3),ANSAV(3), &
     &                WST(4),TS0(4),TSAV(4),ANSAV(4), &
     &                WFT(1),TF0(1),TFAV(1),ANFAV(1), &
     &                WFT(2),TF0(2),TFAV(2),ANFAV(2)
  602    FORMAT(' ',3X,'WE    =',1PD10.3,'  TE0   =',1PD10.3, &
     &               '  TEAVE =',1PD10.3,'  NEAVE =',1PD10.3/ &
     &          ' ',3X,'WD    =',1PD10.3,'  TD0   =',1PD10.3, &
     &               '  TDAVE =',1PD10.3,'  NDAVE =',1PD10.3/ &
     &          ' ',3X,'WT    =',1PD10.3,'  TT0   =',1PD10.3, &
     &               '  TTAVE =',1PD10.3,'  NTAVE =',1PD10.3/ &
     &          ' ',3X,'WA    =',1PD10.3,'  TA0   =',1PD10.3, &
     &               '  TAAVE =',1PD10.3,'  NAAVE =',1PD10.3/ &
     &          ' ',3X,'WB    =',1PD10.3,'  TB0   =',1PD10.3, &
     &               '  TBAVE =',1PD10.3,'  NBAVE =',1PD10.3/ &
     &          ' ',3X,'WF    =',1PD10.3,'  TF0   =',1PD10.3, &
     &               '  TFAVE =',1PD10.3,'  NFAVE =',1PD10.3)

         WRITE(6,603) AJT,VLOOP,ALI,VSEC, &
     &                AJOHT,AJNBT,AJRFT,AJBST
  603    FORMAT(' ',3X,'AJT   =',1PD10.3,'  VLOOP =',1PD10.3, &
     &               '  ALI   =',1PD10.3,'  VSEC  =',1PD10.3/ &
     &          ' ',3X,'AJOHT =',1PD10.3,'  AJNBT =',1PD10.3, &
     &               '  AJRFT =',1PD10.3,'  AJBST =',1PD10.3)

!         WRITE(16,603) AJTTOR,VLOOP,ALI,VSEC, &
!     &                AJT,AJOHT,AJNBT,AJBST
!  603    FORMAT(' ',3X,'AJTTOR=',1PD10.3,'  VLOOP =',1PD10.3, &
!     &               '  ALI   =',1PD10.3,'  VSEC  =',1PD10.3/ &
!     &          ' ',3X,'AJT   =',1PD10.3,'  AJOHT =',1PD10.3, &
!     &               '  AJNBT =',1PD10.3,'  AJBST =',1PD10.3)

         WRITE(6,604) PINT,POHT,PNBT,PNFT, &
     &                PRFT(1),PRFT(2),PRFT(3),PRFT(4), &
     &                PBINT,PFINT,AJ(1)*1.D-6, &
     &                PBCLT(1),PBCLT(2),PBCLT(3),PBCLT(4), &
     &                PFCLT(1),PFCLT(2),PFCLT(3),PFCLT(4), &
     &                POUT,PRSUMT,PCXT,PIET, &
     &                PLT(1),PLT(2),PLT(3),PLT(4), &
                      PRBT,PRCT,PRLT
  604    FORMAT(' ',3X,'PINT  =',1PD10.3,'  POHT  =',1PD10.3, &
     &               '  PNBT  =',1PD10.3,'  PNFTE =',1PD10.3/ &
     &          ' ',3X,'PRFTE =',1PD10.3,'  PRFTD =',1PD10.3, &
     &               '  PRFTT =',1PD10.3,'  PRFTA =',1PD10.3/ &
     &          ' ',3X,'PBIN  =',1PD10.3,'  PFIN  =',1PD10.3, &
     &               '  AJ0   =',1PD10.3/ &
     &          ' ',3X,'PBCLE =',1PD10.3,'  PBCLD =',1PD10.3, &
     &               '  PBCLT =',1PD10.3,'  PBCLA =',1PD10.3/ &
     &          ' ',3X,'PFCLE =',1PD10.3,'  PFCLD =',1PD10.3, &
     &               '  PFCLT =',1PD10.3,'  PFCLA =',1PD10.3/ &
     &          ' ',3X,'POUT  =',1PD10.3,'  PRSUMT=',1PD10.3, &
     &               '  PCXT  =',1PD10.3,'  PIETE =',1PD10.3/ &
     &          ' ',3X,'PLTE  =',1PD10.3,'  PLTD  =',1PD10.3, &
     &               '  PLTTE =',1PD10.3,'  PLTA  =',1PD10.3/ &
     &          ' ',3X,'PRBT  =',1PD10.3,'  PRCT  =',1PD10.3, &
     &               '  PRLT  =',1PD10.3)

         WRITE(6,605) SINT,SIET,SNBT,SNFT, &
     &                SOUT,ZEFF(1),ANC(1),ANFE(1), &
     &                SLT(1),SLT(2),SLT(3),SLT(4)
  605    FORMAT(' ',3X,'SINT  =',1PD10.3,'  SIET  =',1PD10.3, &
     &               '  SNBT  =',1PD10.3,'  SNFT  =',1PD10.3/ &
     &          ' ',3X,'SOUT  =',1PD10.3,'  ZEFF0 =',1PD10.3, &
     &               '  ANC0  =',1PD10.3,'  ANFE0 =',1PD10.3/ &
     &          ' ',3X,'SLTET =',1PD10.3,'  SLTD  =',1PD10.3, &
     &               '  SLTTT =',1PD10.3,'  SLTA  =',1PD10.3)
      ENDIF

      IF(KID.EQ.'2') THEN
         WRITE(6,611) Q0,(QP(I),I=1,NRMAX)
  611    FORMAT(' ','* Q PROFILE *'/ &
     &         (' ',5F7.3,2X,5F7.3))
      ENDIF

      IF(KID.EQ.'3') THEN
         CALL GUTIME(GTCPU2)
         WRITE(6,621) GTCPU2-GTCPU1
  621    FORMAT(' ','# CPU TIME = ',F8.3,' S')
         RETURN
      ENDIF

      IF(KID.EQ.'4') THEN
         WRITE(6,631)
  631    FORMAT(' ','#',12X,'FIRST',7X,'MAX',9X,'MIN',9X,'LAST')
         CALL TRMXMN( 1,'  NE0  ')
!         CALL TRMXMN( 2,'  ND0  ')
!         CALL TRMXMN( 3,'  NT0  ')
!         CALL TRMXMN( 4,'  NA0  ')
         CALL TRMXMN( 5,'  NEAV ')
!         CALL TRMXMN( 6,'  NDAV ')
!         CALL TRMXMN( 7,'  NTAV ')
!         CALL TRMXMN( 8,'  NAAV ')
         CALL TRMXMN( 9,'  TE0  ')
         CALL TRMXMN(10,'  TD0  ')
         CALL TRMXMN(11,'  TT0  ')
!         CALL TRMXMN(12,'  TA0  ')
         CALL TRMXMN(13,'  TEAV ')
         CALL TRMXMN(14,'  TDAV ')
!         CALL TRMXMN(15,'  TTAV ')
!         CALL TRMXMN(16,'  TAAV ')
!         CALL TRMXMN(17,'  WE   ')
!         CALL TRMXMN(18,'  WD   ')
!         CALL TRMXMN(19,'  WT   ')
!         CALL TRMXMN(20,'  WA   ')

!         CALL TRMXMN(21,'  NB0  ')
!         CALL TRMXMN(22,'  NF0  ')
!         CALL TRMXMN(23,'  NBAV ')
!         CALL TRMXMN(24,'  NFAV ')
!         CALL TRMXMN(25,'  TB0  ')
!         CALL TRMXMN(26,'  TF0  ')
!         CALL TRMXMN(27,'  TBAV ')
!         CALL TRMXMN(28,'  TFAV ')
         CALL TRMXMN(29,'  WB   ')
         CALL TRMXMN(30,'  WF   ')
         CALL TRMXMN(31,' WBULK ')
!         CALL TRMXMN(32,' WTAIL ')
         CALL TRMXMN(33,'  WP   ')

!         CALL TRMXMN(34,'  IP   ')
         CALL TRMXMN(35,'  IOH  ')
         CALL TRMXMN(36,'  INB  ')
!         CALL TRMXMN(37,'  IRF  ')
         CALL TRMXMN(38,'  IBS  ')

!         CALL TRMXMN(39,'  PIN  ')
         CALL TRMXMN(40,'  POH  ')
         CALL TRMXMN(41,'  PNB  ')
!         CALL TRMXMN(42,'  PRFE ')
!         CALL TRMXMN(43,'  PRFD ')
!         CALL TRMXMN(44,'  PRFT ')
!         CALL TRMXMN(45,'  PRFA ')
         CALL TRMXMN(46,'  PNF  ')
!         CALL TRMXMN(47,'  PBINT')
!         CALL TRMXMN(48,'  PBCLE')
!         CALL TRMXMN(49,'  PBCLD')
!         CALL TRMXMN(50,'  PBCLT')
!         CALL TRMXMN(51,'  PBCLA')
!         CALL TRMXMN(52,'  PFINT')
!         CALL TRMXMN(53,'  PFCLE')
!         CALL TRMXMN(54,'  PFCLD')
!         CALL TRMXMN(55,'  PFCLT')
!         CALL TRMXMN(56,'  PFCLA')
!         CALL TRMXMN(57,'  POUT ')
         CALL TRMXMN(58,'  PCX  ')
         CALL TRMXMN(59,'  PIE  ')
         CALL TRMXMN(60,'  PRL  ')
         CALL TRMXMN(61,'  PLE  ')
         CALL TRMXMN(62,'  PLD  ')
         CALL TRMXMN(63,'  PLT  ')
         CALL TRMXMN(64,'  PLA  ')

!         CALL TRMXMN(65,'  SIN  ')
!         CALL TRMXMN(66,'  SIE  ')
!         CALL TRMXMN(67,'  SNB  ')
!         CALL TRMXMN(68,'  SNF  ')
!         CALL TRMXMN(69,'  SOUT ')
!         CALL TRMXMN(70,'  SLE  ')
!         CALL TRMXMN(71,'  SLD  ')
!         CALL TRMXMN(72,'  SLT  ')
!         CALL TRMXMN(73,'  SLA  ')

         CALL TRMXMN(74,' VLOOP ')
         CALL TRMXMN(75,'  LI   ')
!         CALL TRMXMN(76,'  RQ1  ')
!         CALL TRMXMN(77,'   Q0  ')
!         CALL TRMXMN(78,' WPDOT ')
!         CALL TRMXMN(79,' TAUE1 ')
!         CALL TRMXMN(80,' TAUE2 ')
!         CALL TRMXMN(81,' TAUE89')
!         CALL TRMXMN(82,' BETAP0')
         CALL TRMXMN(83,' BETAPA')
!         CALL TRMXMN(84,' BETA0 ')
!         CALL TRMXMN(85,' BETAA ')
!         CALL TRMXMN(86,' ZEFF0 ')
         CALL TRMXMN(87,'   QF  ')
!         CALL TRMXMN(88,'   IP  ')
         CALL TRMXMN(89,'  PEX  ')
      ENDIF

      IF(KID.EQ.'5') THEN
         WRITE(6,641)NGR,NGT,DT,NTMAX
  641    FORMAT(' ','# PARAMETER INFORMATION',/ &
     &          ' ','  NGR   =',I3,'    NGT   =',I3,/ &
     &          ' ','  DT    =',1F5.3,'  NTMAX =',I3)
      ENDIF

      IF(KID.EQ.'6') THEN
         WRITE(6,651)T,TAUE1,TAUE2,TAUE89,PINT
 651     FORMAT(' ','# TIME : ',F7.3,' SEC'/ &
     &          ' ',3X,'TAUE1 =',1PD10.3,'  TAUE2 =',1PD10.3, &
     &               '  TAUE89=',1PD10.3,'  PINT  =',1PD10.3)
      ENDIF

      IF(KID.EQ.'7'.OR.KID.EQ.'8') THEN
         WRITE(6,671) T,WPT,TAUE1,TAUE2,TAUE89,BETAN,BETAPA,BETA0,BETAA
  671    FORMAT(' ','# TIME : ',F7.3,' SEC'/ &
     &          ' ',3X,'WPT   =',1PD10.3,'  TAUE1 =',1PD10.3, &
     &               '  TAUE2 =',1PD10.3,'  TAUE89=',1PD10.3/ &
     &          ' ',3X,'BETAN =',1PD10.3,'  BETAPA=',1PD10.3, &
     &               '  BETA0 =',1PD10.3,'  BETAA =',1PD10.3)

         WRITE(6,672) WST(1),TS0(1),TSAV(1),ANSAV(1), &
     &                WST(2),TS0(2),TSAV(2),ANSAV(2)
  672    FORMAT(' ',3X,'WE    =',1PD10.3,'  TE0   =',1PD10.3, &
     &               '  TEAVE =',1PD10.3,'  NEAVE =',1PD10.3/ &
     &          ' ',3X,'WD    =',1PD10.3,'  TD0   =',1PD10.3, &
     &               '  TDAVE =',1PD10.3,'  NDAVE =',1PD10.3)

         WRITE(6,673) AJTTOR,VLOOP,ALI,Q0,AJOHT,AJNBT,AJRFT,AJBST
  673    FORMAT(' ',3X,'AJTTOR=',1PD10.3,'  VLOOP =',1PD10.3, &
     &               '  ALI   =',1PD10.3,'  Q0    =',1PD10.3/ &
     &          ' ',3X,'AJOHT =',1PD10.3,'  AJNBT =',1PD10.3, &
     &               '  AJRFT =',1PD10.3,'  AJBST =',1PD10.3)

         WRITE(6,674) PINT,POHT,PNBT, &
     &                PRFT(1)+PRFT(2)+PRFT(3)+PRFT(4),POUT,PRLT,PCXT,PIET
  674    FORMAT(' ',3X,'PINT  =',1PD10.3,'  POHT  =',1PD10.3, &
     &               '  PNBT  =',1PD10.3,'  PRFT  =',1PD10.3/ &
     &          ' ',3X,'POUT  =',1PD10.3,'  PRLT  =',1PD10.3, &
     &               '  PCXT  =',1PD10.3,'  PIETE =',1PD10.3)

      IF(KID.EQ.'8') THEN
 1600    WRITE(6,*) '## INPUT COMMENT FOR trn.data (A40)'
         READ(5,'(A40)',END=9000,ERR=1600) KCOM

!         OPEN(16,POSITION='APPEND',FILE=KFNLOG)
!         OPEN(16,ACCESS='APPEND',FILE=KFNLOG)
         OPEN(16,ACCESS='SEQUENTIAL',FILE=KFNLOG)

         CALL GUDATE(NDY,NDM,NDD,NTH1,NTM1,NTS1)
         WRITE(K1,'(I3)') 100+NDY
         WRITE(K2,'(I3)') 100+NDM
         WRITE(K3,'(I3)') 100+NDD
         WRITE(K4,'(I3)') 100+NTH1
         WRITE(K5,'(I3)') 100+NTM1
         WRITE(K6,'(I3)') 100+NTS1
         WRITE(16,1670) K1(2:3),K2(2:3),K3(2:3),K4(2:3),K5(2:3),K6(2:3), &
     &                  KCOM, &
     &                  RIPS,RIPE,PN(1),PN(2),BB,PICTOT,PLHTOT,PLHNPR
 1670    FORMAT(' '/ &
     &          ' ','## DATE : ', &
     &              A2,'-',A2,'-',A2,'  ',A2,':',A2,':',A2,' : ',A40/ &
     &          ' ',3X,'RIPS  =',1PD10.3,'  RIPE  =',1PD10.3, &
     &               '  PNE   =',1PD10.3,'  PNI   =',1PD10.3/ &
     &          ' ',3X,'BB    =',1PD10.3,'  PICTOT=',1PD10.3, &
     &               '  PLHTOT=',1PD10.3,'  PLHNPR=',1PD10.3)
         WRITE(16,1671) T, &
     &                WPT,TAUE1,TAUE2,TAUE89, &
     &                BETAN,BETAPA,BETA0,BETAA
 1671    FORMAT(' ','# TIME : ',F7.3,' SEC'/ &
     &          ' ',3X,'WPT   =',1PD10.3,'  TAUE  =',1PD10.3, &
     &               '  TAUED =',1PD10.3,'  TAUE89=',1PD10.3/ &
     &          ' ',3X,'BETAN =',1PD10.3,'  BETAPA=',1PD10.3, &
     &               '  BETA0 =',1PD10.3,'  BETAA =',1PD10.3)

         WRITE(16,1672) WST(1),TS0(1),TSAV(1),ANSAV(1), &
     &                WST(2),TS0(2),TSAV(2),ANSAV(2)
 1672    FORMAT(' ',3X,'WE    =',1PD10.3,'  TE0   =',1PD10.3, &
     &               '  TEAVE =',1PD10.3,'  NEAVE =',1PD10.3/ &
     &          ' ',3X,'WD    =',1PD10.3,'  TD0   =',1PD10.3, &
     &               '  TDAVE =',1PD10.3,'  NDAVE =',1PD10.3)

         WRITE(16,1673) AJTTOR,VLOOP,ALI,Q0, &
     &                AJOHT,AJNBT,AJRFT,AJBST
 1673    FORMAT(' ',3X,'AJTTOR=',1PD10.3,'  VLOOP =',1PD10.3, &
     &               '  ALI   =',1PD10.3,'  Q0    =',1PD10.3/ &
     &          ' ',3X,'AJOHT =',1PD10.3,'  AJNBT =',1PD10.3, &
     &               '  AJRFT =',1PD10.3,'  AJBST =',1PD10.3)

         WRITE(16,1674) PINT,POHT,PNBT, &
     &                PRFT(1)+PRFT(2)+PRFT(3)+PRFT(4), &
     &                POUT,PRLT,PCXT,PIET
 1674    FORMAT(' ',3X,'PINT  =',1PD10.3,'  POHT  =',1PD10.3, &
     &               '  PNBT  =',1PD10.3,'  PRFT  =',1PD10.3/ &
     &          ' ',3X,'POUT  =',1PD10.3,'  PRLT  =',1PD10.3, &
     &               '  PCXT  =',1PD10.3,'  PIETE =',1PD10.3)
         CLOSE(16)
      ENDIF
      ENDIF

      IF(KID.EQ.'9') THEN
         CALL TRDATA
      ENDIF

 9000 RETURN
      END SUBROUTINE TRPRNT

!     ***********************************************************

!           PRINT LOCAL DATA

!     ***********************************************************

      SUBROUTINE TRDATA

      USE TRCOMM, ONLY : GRG, GRM, GT, GVR, GVT, NGR, NT
      IMPLICIT NONE
      INTEGER:: MGH, MGMAX, MGMIN, MGSTEP, MID, MRMAX, MRMIN, MRSTEP, MTMAX, MTMIN, MTSTEP, NG, NID, NR


    1 WRITE(6,*) '## INPUT MODE : 1:GVT(NT)  2:GVR(NR)  3:GVR(NG)'
      WRITE(6,*) '                NOW NGR=',NGR
      READ(5,*,END=9000,ERR=1) NID
      IF(NID.EQ.0) GOTO 9000

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

 9000 RETURN
  601 FORMAT(' ','  NT=',I3,'  T=',1PE12.4,'  DATA=',1PE12.4)
  602 FORMAT(' ','  NG=',I3,'  NR=',I3,'  R=',1PE12.4,    '  DATA=',1PE12.4)
      END SUBROUTINE TRDATA

!     ***********************************************************

!          PEAK-VALUE WO SAGASE

!     ***********************************************************

      SUBROUTINE TRMXMN(N,STR)

      USE TRCOMM, ONLY : GVT, NGT, NT
      IMPLICIT NONE
      INTEGER:: N
      REAL   :: GVMAX, GVMIN

      CHARACTER STR*7

      GVMAX=GVT(1,N)
      GVMIN=GVT(1,N)

      DO NT=2,NGT
         GVMAX=MAX(GVMAX,GVT(NT,N))
         GVMIN=MIN(GVMIN,GVT(NT,N))
      ENDDO

      WRITE(6,600) STR,GVT(1,N),GVMAX,GVMIN,GVT(NGT,N)
  600 FORMAT(' ',A8,5X,1PD10.3,2X,1PD10.3,2X,1PD10.3,2X,1PD10.3)

      RETURN
      END SUBROUTINE TRMXMN

!     ***********************************************************

!           SIMPLE STATUS REPORT

!     ***********************************************************

      SUBROUTINE TRSNAP

      USE TRCOMM, ONLY : Q0, RT, T, TAUE1, WPT
      IMPLICIT NONE


      WRITE(6,601) T,WPT,TAUE1,Q0,RT(1,1),RT(1,2),RT(1,3),RT(1,4)
  601 FORMAT(' ','# T: ',F8.3,'(S)    WP:',F7.2,'(MJ)  ', &
     &           '  TAUE:',F7.3,'(S)   Q0:',F7.3,/ &
     &       ' ','  TE:',F7.3,'(KEV)   TD:',F7.3,'(KEV) ', &
     &           '  TT:',F7.3,'(KEV)   TA:',F7.3,'(KEV)')
      RETURN
      END SUBROUTINE TRSNAP


!     ***********************************************************

!           SAVE PROFILE DATA

!     ***********************************************************

      SUBROUTINE TRXOUT

        USE TRCOMM, ONLY : MDLUF,KXNDEV,KXNDCG,KXNID,KDIRW1,KDIRW2
        IMPLICIT NONE
        INTEGER:: IERR, IKDIRW, IKNDCG, IKNDEV, IKNID
        CHARACTER(LEN=80):: KDIRW, KFID


      KXNDEV='X'
      KXNDCG='test'
      KXNID ='in'

      IKNDEV=len_trim(KXNDEV)
      IKNDCG=len_trim(KXNDCG)
      IKNID =len_trim(KXNID )
      KDIRW='./profile_data/'//KXNDEV(1:IKNDEV)//'/' &
     &                          //KXNDCG(1:IKNDCG)//'/' &
     &                          //KXNID (1:IKNID )//'/'
!      KDIRW='../../tr.new/data/'//KXNDEV(1:IKNDEV)//'/' &
!     &                          //KXNDCG(1:IKNDCG)//'/' &
!     &                          //KXNID (1:IKNID )//'/'
!!      KDIRW='../../../profile/profile_data/'//KXNDEV(1:IKNDEV)//'/'
!!     &                          //KXNDCG(1:IKNDCG)//'/'
!!     &                          //KXNID (1:IKNID )//'/'
      IKDIRW=len_trim(KDIRW)
      KDIRW1=KDIRW(1:IKDIRW)//KXNDEV(1:IKNDEV) &
     &       //'1d'//KXNDCG(1:IKNDCG)//'.'
      KDIRW2=KDIRW(1:IKDIRW)//KXNDEV(1:IKNDEV) &
     &       //'2d'//KXNDCG(1:IKNDCG)//'.'

!     *** 1D DATA ***

      KFID='IP'
      CALL TR_UFILE1D_CREATE(KFID,34,1.D6 ,IERR)

      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 2,1.D0 ,IERR)

      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 3,1.D0 ,IERR)

      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 4,1.D0 ,IERR)

      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 5,1.D0 ,IERR)

      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 6,1.D0 ,IERR)

      IF(MDLUF.NE.0) THEN
      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 8,1.D0 ,IERR)
      ELSE
      KFID='PNBI'
      CALL TR_UFILE1D_CREATE(KFID,41,1.D0 ,IERR)
      ENDIF

      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID, 9,1.D0 ,IERR)

      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID,10,1.D0 ,IERR)

      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID,11,1.D0 ,IERR)

      KFID='PRAD'
      CALL TR_UFILE1D_CREATE(KFID,60,1.D0 ,IERR)

      KFID='ZEFF'
      CALL TR_UFILE1D_CREATE(KFID,86,1.D0 ,IERR)

      KFID='VSURF'
      CALL TR_UFILE1D_CREATE(KFID,74,1.D0 ,IERR)

      KFID='LI'
      CALL TR_UFILE1D_CREATE(KFID,75,1.D0 ,IERR)

      KFID='WTH'
      CALL TR_UFILE1D_CREATE(KFID,31,1.D6 ,IERR)

      KFID='WTOT'
      CALL TR_UFILE1D_CREATE(KFID,33,1.D6 ,IERR)

      KFID='TE0'
      CALL TR_UFILE1D_CREATE(KFID, 9,1.D3 ,IERR)

      KFID='TI0'
      CALL TR_UFILE1D_CREATE(KFID,10,1.D3 ,IERR)

      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID,26,1.D0 ,IERR)

      KFID='POHM'
      CALL TR_UFILE1D_CREATE(KFID,40,1.D0 ,IERR)

      KFID='IBOOT'
      CALL TR_UFILE1D_CREATE(KFID,38,1.D6 ,IERR)

      KFID='DIRECT'
      CALL TR_UFILE1D_CREATE(KFID,29,1.D0 ,IERR)

      KFID='PFUSION'
      CALL TR_UFILE1D_CREATE(KFID,46,1.D0 ,IERR)

!     *** 2D DATA ***

      KFID='TE'
      CALL TR_UFILE2D_CREATE(KFID, 1,1.D3 ,0,IERR)

      KFID='TI'
      CALL TR_UFILE2D_CREATE(KFID, 2,1.D3 ,0,IERR)

      KFID='NE'
      CALL TR_UFILE2D_CREATE(KFID, 5,1.D20,0,IERR)

      IF(MDLUF.NE.0) THEN
      KFID='QNBIE'
      CALL TR_UFILE2D_CREATE(KFID,89,1.D0 ,0,IERR)
      ENDIF

      KFID='QICRHE'
      CALL TR_UFILE2D_CREATE(KFID,38,1.D0 ,0,IERR)

      KFID='QECHE'
      CALL TR_UFILE2D_CREATE(KFID,36,1.D0 ,0,IERR)

      KFID='QLHE'
      CALL TR_UFILE2D_CREATE(KFID,37,1.D0 ,0,IERR)

      IF(MDLUF.NE.0) THEN
      KFID='QNBII'
      CALL TR_UFILE2D_CREATE(KFID,90,1.D0 ,0,IERR)
      ENDIF

      KFID='QICRHI'
      CALL TR_UFILE2D_CREATE(KFID,41,1.D0 ,0,IERR)

      KFID='QECHI'
      CALL TR_UFILE2D_CREATE(KFID,39,1.D0 ,0,IERR)

      KFID='QLHI'
      CALL TR_UFILE2D_CREATE(KFID,40,1.D0 ,0,IERR)

      KFID='CURNBI'
      CALL TR_UFILE2D_CREATE(KFID,11,1.D0 ,0,IERR)

      KFID='CURICRH'
      CALL TR_UFILE2D_CREATE(KFID,44,1.D0 ,0,IERR)

      KFID='CURECH'
      CALL TR_UFILE2D_CREATE(KFID,42,1.D0 ,0,IERR)

      KFID='CURLH'
      CALL TR_UFILE2D_CREATE(KFID,43,1.D0 ,0,IERR)

      KFID='NFAST'
      CALL TR_UFILE2D_CREATE(KFID,45,1.D20,0,IERR)

      KFID='QRAD'
      CALL TR_UFILE2D_CREATE(KFID,22,1.D0 ,0,IERR)

      KFID='ZEFFR'
      CALL TR_UFILE2D_CREATE(KFID,33,1.D0 ,0,IERR)

      KFID='Q'
      CALL TR_UFILE2D_CREATE(KFID,27,1.D0 ,1,IERR)

      KFID='CHIE'
      CALL TR_UFILE2D_CREATE(KFID,34,1.D0 ,1,IERR)

      KFID='CHII'
      CALL TR_UFILE2D_CREATE(KFID,35,1.D0 ,1,IERR)

      KFID='NM1'
      CALL TR_UFILE2D_CREATE(KFID, 6,1.D20,0,IERR)

      KFID='CURTOT'
      CALL TR_UFILE2D_CREATE(KFID, 9,1.D0 ,0,IERR)

      KFID='NIMP'
      CALL TR_UFILE2D_CREATE(KFID,46,1.D20,0,IERR)

      KFID='QOHM'
      CALL TR_UFILE2D_CREATE(KFID,15,1.D0 ,0,IERR)

      KFID='BPOL'
      CALL TR_UFILE2D_CREATE(KFID,47,1.D0 ,1,IERR)

      KFID='RMAJOR'
      CALL TR_UFILE2D_CREATE(KFID,49,1.D0 ,0,IERR)

      KFID='RMINOR'
      CALL TR_UFILE2D_CREATE(KFID,50,1.D0 ,0,IERR)

      KFID='VOLUME'
      CALL TR_UFILE2D_CREATE(KFID,51,1.D0 ,0,IERR)

      KFID='KAPPAR'
      CALL TR_UFILE2D_CREATE(KFID,52,1.D0 ,0,IERR)

      KFID='DELTAR'
      CALL TR_UFILE2D_CREATE(KFID,53,1.D0 ,0,IERR)

      KFID='GRHO1'
      CALL TR_UFILE2D_CREATE(KFID,54,1.D0 ,0,IERR)

      KFID='GRHO2'
      CALL TR_UFILE2D_CREATE(KFID,55,1.D0 ,0,IERR)

      KFID='CURBS'
      CALL TR_UFILE2D_CREATE(KFID,13,1.D0 ,0,IERR)

      KFID='CHITBE'
      CALL TR_UFILE2D_CREATE(KFID,56,1.D0 ,1,IERR)

      KFID='CHITBI'
      CALL TR_UFILE2D_CREATE(KFID,57,1.D0 ,1,IERR)

      KFID='ETAR'
      CALL TR_UFILE2D_CREATE(KFID,32,1.D0 ,0,IERR)

      RETURN
      END SUBROUTINE TRXOUT

!     *****

      SUBROUTINE TR_UFILE1D_CREATE(KFID,NUM,AMP,IERR)

      USE TRCOMM, ONLY : BB, GVRT, GT, GVT, MDLUF, NGT, NRAMAX, NRMAX, NROMAX, NTM, RA, RG, RHOA, RKAP, RR, rkind
      USE libspl1d
      IMPLICIT NONE
      CHARACTER(LEN=80),INTENT(INOUT):: KFID
      INTEGER,INTENT(IN) :: NUM
      REAL(rkind),   INTENT(IN) :: AMP
      INTEGER,INTENT(OUT):: IERR
      INTEGER:: ID, NTL, NTLMAX
      REAL(rkind)   :: DATOUT, DTL, FQ95, TIN
      REAL,DIMENSION(NTM)  :: GTL, GF1
      REAL(rkind),DIMENSION(NTM)  :: TF, DGT, DIN, DERIV
      REAL(rkind),DIMENSION(NRMAX)  :: F1, DERIVQ
      REAL(rkind),DIMENSION(4,NTM):: UOUT
      REAL(rkind),DIMENSION(4,NRMAX):: UQ95
      CHARACTER(LEN=80)::KERRF
      REAL :: GUCLIP

      IF(KFID.EQ.'DIRECT') THEN
         IF(NUM.EQ.2) THEN
            KFID='BT'
            TF(1:NGT)=BB
         ELSEIF(NUM.EQ.3) THEN
            KFID='AMIN'
            TF(1:NGT)=RA
         ELSEIF(NUM.EQ.4) THEN
            KFID='RGEO'
            TF(1:NGT)=RR
         ELSEIF(NUM.EQ.5) THEN
            KFID='KAPPA'
            TF(1:NGT)=RKAP
         ELSEIF(NUM.EQ.6) THEN
            KFID='DELTA'
            TF(1:NGT)=0.D0
         ELSEIF(NUM.EQ.8) THEN
            KFID='PNBI'
            TF(1:NGT)=DBLE(GVT(1:NGT,89)+GVT(1:NGT,90))
         ELSEIF(NUM.EQ.9) THEN
            KFID='PECH'
            TF(1:NGT)=DBLE(GVT(1:NGT,91)+GVT(1:NGT,92))
         ELSEIF(NUM.EQ.10) THEN
            KFID='PICRH'
            TF(1:NGT)=DBLE(GVT(1:NGT,95)+GVT(1:NGT,96))
         ELSEIF(NUM.EQ.11) THEN
            KFID='PLH'
            TF(1:NGT)=DBLE(GVT(1:NGT,93)+GVT(1:NGT,94))
         ELSEIF(NUM.EQ.26) THEN
            KFID='Q95'
            ID=0
            IF(MDLUF.NE.0.AND.NRMAX.NE.NROMAX) THEN
               ID=1
               NRMAX=NROMAX
            ENDIF
            DO NTL=1,NGT
               F1(1:NRMAX)=DBLE(GVRT(1:NRMAX,NTL,27))
               CALL SPL1D (RG,F1,DERIVQ,UQ95,NRMAX,0,IERR)
               CALL SPL1DF(0.95D0,FQ95,RG,UQ95,NRMAX,IERR)
               TF(NTL)=FQ95
            ENDDO
            IF(ID.NE.0) NRMAX=NRAMAX
         ELSEIF(NUM.EQ.29) THEN
            KFID='RHOA'
            TF(1:NGT)=RHOA
         ENDIF

         DGT(1:NGT)=DBLE(GT(1:NGT))
         DIN(1:NGT)=DBLE(TF(1:NGT))
         DERIV(1:NGT)=0.D0
      ELSE
         DGT(1:NGT)=DBLE(GT(1:NGT))
         DIN(1:NGT)=DBLE(GVT(1:NGT,NUM))
         DERIV(1:NGT)=0.D0
      ENDIF

      CALL SPL1D(DGT,DIN,DERIV,UOUT,NGT,0,IERR)
      IF(IERR.NE.0) THEN
         IF(KFID.EQ.'DIRECT') THEN
            WRITE(6,'(A,I2,A,I2)') 'XX TRXOUT: SPL1D DIRECT(',NUM,'): IERR=',IERR
         ELSE
            WRITE(6,'(A,I2,A,I2)') 'XX TRXOUT: SPL1D GVT(',NUM,'): IERR=',IERR
         ENDIF
      ENDIF

      DTL=0.05D0
      NTLMAX=INT((GT(NGT)-GT(1))/SNGL(DTL))+1

      DO NTL=1,NTLMAX
         TIN=DBLE(GT(1))+DTL*DBLE(NTL-1)
         CALL SPL1DF(TIN,DATOUT,DGT,UOUT,NGT,IERR)
         WRITE(KERRF,'(A,I2,A,I2)') 'XX TRXOUT: SPL1DF GVT(',NUM,'): IERR=',IERR
         IF(IERR.NE.0) WRITE(6,*) KERRF
         GTL(NTL)=GUCLIP(TIN)
         GF1(NTL)=GUCLIP(DATOUT*AMP)
      ENDDO

      CALL TRXW1D(KFID,GTL,GF1,NTM,NTLMAX)

      RETURN
      END SUBROUTINE TR_UFILE1D_CREATE

!     *****

      SUBROUTINE TR_UFILE2D_CREATE(KFID,NUM,AMP,ID,IERR)

      USE TRCOMM, ONLY : GVRT, GRG, GRM, GT, NGT, NRMAX, NRMP, NTM, rkind
      USE libitp  
      USE libspl1d
      IMPLICIT NONE
      CHARACTER(LEN=80),INTENT(IN) :: KFID
      INTEGER       ,INTENT(IN) :: NUM, ID
      REAL(rkind)          ,INTENT(IN) :: AMP
      INTEGER       ,INTENT(OUT):: IERR
      INTEGER::NTL, NRLMAX, NTLMAX, NRL
      REAL(rkind)   ::DTL, TIN, F0, R1, R2, F1, F2
      REAL,DIMENSION(NRMP)    :: GRL
      REAL,DIMENSION(NTM)     :: GTL
      REAL,DIMENSION(NRMP,NTM):: GF2
      REAL(rkind),DIMENSION(NTM)     :: DGT,DIN,DERIV
      REAL(rkind),DIMENSION(4,NTM)   :: U
      REAL   :: GUCLIP


      DGT(1:NGT)=DBLE(GT(1:NGT))

      NRLMAX=NRMAX
      DTL=0.05D0
      NTLMAX=INT((GT(NGT)-GT(1))/SNGL(DTL))+1
      IF(ID.EQ.0) THEN
         DO NRL=1,NRLMAX
            GRL(NRL)=GRM(NRL)
            DIN(1:NGT)=DBLE(GVRT(NRL,1:NGT,NUM))
            CALL SPL1D(DGT,DIN,DERIV,U,NGT,0,IERR)
            IF(IERR.NE.0) WRITE(6,'(A,I2,A,I2)') 'XX TRXOUT: SPL1D GVRT(',NUM,'): IERR=',IERR
!
            DO NTL=1,NTLMAX
               TIN=DBLE(GT(1))+DTL*DBLE(NTL-1)
               CALL SPL1DF(TIN,F0,DGT,U,NGT,IERR)
               IF(IERR.NE.0) WRITE(*,'(A,I2,A,I2)') 'XX TRXOUT: SPL1DF GVRT(',NUM,'): IERR=',IERR
               GTL(NTL)    =GUCLIP(TIN)
               GF2(NRL,NTL)=GUCLIP(F0*AMP)
            ENDDO
         ENDDO
      ELSEIF(ID.EQ.1) THEN
         NRLMAX=NRMAX+1
         NRL=1
            GRL(NRL)=GRG(NRL)
            IF(KFID.EQ.'Q') THEN
               DIN(1:NGT)=(4.D0*DBLE(GVRT(NRL,1:NGT,NUM))-DBLE(GVRT(NRL+1,1:NGT,NUM)))/3.D0
            ELSEIF(KFID.EQ.'BPOL') THEN
               DIN(1:NGT)=0.D0
            ELSE
               DO NTL=1,NGT
                  R1=DBLE(GRL(NRL))
                  R2=DBLE(GRL(NRL+1))
                  F1=DBLE(GVRT(NRL  ,NTL,NUM))
                  F2=DBLE(GVRT(NRL+1,NTL,NUM))
                  DIN(NTL)=FCTR(R1,R2,F1,F2)
               ENDDO
            ENDIF
            CALL SPL1D(DGT,DIN,DERIV,U,NGT,0,IERR)
            IF(IERR.NE.0) WRITE(6,'(A,I2,A,I2)') 'XX TRXOUT: SPL1D GVRT(',NUM,'): IERR=',IERR
            DO NTL=1,NTLMAX
               TIN=DBLE(GT(1))+DTL*DBLE(NTL-1)
               CALL SPL1DF(TIN,F0,DGT,U,NGT,IERR)
               IF(IERR.NE.0) WRITE(*,'(A,I2,A,I2)') 'XX TRXOUT: SPL1DF GVRT(',NUM,'): IERR=',IERR
               GTL(NTL)    =GUCLIP(TIN)
               GF2(NRL,NTL)=GUCLIP(F0*AMP)
            ENDDO

         DO NRL=2,NRLMAX
            GRL(NRL)=GRG(NRL)
            DIN(1:NGT)=DBLE(GVRT(NRL-1,1:NGT,NUM))
            CALL SPL1D(DGT,DIN,DERIV,U,NGT,0,IERR)
            IF(IERR.NE.0) WRITE(6,'(A,I2,A,I2)') 'XX TRXOUT: SPL1D GVRT(',NUM,'): IERR=',IERR
            DO NTL=1,NTLMAX
               TIN=DBLE(GT(1))+DTL*DBLE(NTL-1)
               CALL SPL1DF(TIN,F0,DGT,U,NGT,IERR)
               IF(IERR.NE.0) WRITE(*,'(A,I2,A,I2)') 'XX TRXOUT: SPL1DF GVRT(',NUM,'): IERR=',IERR
               GTL(NTL)    =GUCLIP(TIN)
               GF2(NRL,NTL)=GUCLIP(F0*AMP)
            ENDDO
         ENDDO
      ENDIF
      CALL TRXW2D(KFID,GTL,GRL,GF2,NRMP,NTM,NRLMAX,NTLMAX)

      RETURN
      END SUBROUTINE TR_UFILE2D_CREATE

!     *****

      SUBROUTINE TRXW1D(KFID,GT,GF,NTM,NTXMAX)

      USE TRCOMM, ONLY : KDIRW1, KXNDCG, KXNDEV
      IMPLICIT NONE
      CHARACTER(LEN=80)     ,INTENT(IN):: KFID
      INTEGER            ,INTENT(IN):: NTM, NTXMAX
      REAL,DIMENSION(NTM),INTENT(IN):: GT, GF
      INTEGER:: KL1, IST, NTX
      CHARACTER(LEN=9) :: CDATE
      CHARACTER(LEN=80):: KFILE


      CALL GET_DATE(CDATE)

      KL1=len_trim(KDIRW1)
      KFILE=KDIRW1(1:KL1)//KFID
      WRITE(6,*) '- OPEN FILE:',KFILE(1:55)

      OPEN(16,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',ERR=10)

      WRITE(16,'(1X,A8,A8,A14,A29,A9)') KXNDCG(1:8),KXNDEV(1:8), &
     &     '               ', &
     &     ';-SHOT #- F(X) DATA -UF1DWR- ',CDATE
      WRITE(16,'(1X,A10,A20,A38)') 'TR:/tasktr','                    ', &
     &     ';-SHOT DATE-  UFILES ASCII FILE SYSTEM'
      WRITE(16,'(1X,A30,A29)') &
     &     'TIME          SECONDS         ', &
     &     ';-INDEPENDENT VARIABLE LABEL-'
      WRITE(16,'(1X,A30,A27)') KFID, &
     &     ';-DEPENDENT VARIABLE LABEL-'
      WRITE(16,'(1X,I1,A29,A39)') 2,'                             ', &
     &     ';-PROC CODE- 0:RAW 1:AVG 2:SM. 3:AVG+SM'
      WRITE(16,'(1X,I11,A19,A33)') NTXMAX,'                   ', &
     &     ';-# OF PTS-  X, F(X) DATA FOLLOW:'

      WRITE(16,'(1X,1P6E13.6)') (GT(NTX),NTX=1,NTXMAX)
      WRITE(16,'(1X,1P6E13.6)') (GF(NTX),NTX=1,NTXMAX)

      WRITE(16,'(A52)') &
     &     ';----END-OF-DATA-----------------COMMENTS:-----------'

      CLOSE(16)
      RETURN

   10 WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
      RETURN
      END SUBROUTINE TRXW1D

!     *****

      SUBROUTINE TRXW2D(KFID,GT,GR,GF,NRM,NTM,NRXMAX,NTXMAX)

      USE TRCOMM, ONLY :KDIRW2, KXNDCG, KXNDEV
      IMPLICIT NONE
      CHARACTER(LEN=80)         ,INTENT(IN):: KFID
      INTEGER                ,INTENT(IN):: NRM, NTM, NRXMAX, NTXMAX
      REAL,DIMENSION(NTM)    ,INTENT(IN):: GT
      REAL,DIMENSION(NRM)    ,INTENT(IN):: GR
      REAL,DIMENSION(NRM,NTM),INTENT(IN):: GF
      INTEGER:: KL2,IST,NRX,NTX
      CHARACTER(LEN=9) :: CDATE
      CHARACTER(LEN=80):: KFILE


      CALL GET_DATE(CDATE)

      KL2=len_trim(KDIRW2)
      KFILE=KDIRW2(1:KL2)//KFID
      WRITE(6,*) '- OPEN FILE:',KFILE(1:55)

      OPEN(16,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',ERR=10)

      WRITE(16,'(1X,A8,A8,A14,A29,A9)') KXNDCG(1:8),KXNDEV(1:8), &
     &     '               ', &
     &     ';-SHOT #- F(X) DATA -UF1DWR- ',CDATE
      WRITE(16,'(1X,A10,A20,A38)') 'TR:/tasktr','                    ', &
     &     ';-SHOT DATE-  UFILES ASCII FILE SYSTEM'
      WRITE(16,'(1X,A30,A32)') &
     &     'RHO                           ', &
     &     ';-INDEPENDENT VARIABLE LABEL: X-'
      WRITE(16,'(1X,A30,A32)') &
     &     'TIME          SECONDS         ', &
     &     ';-INDEPENDENT VARIABLE LABEL: Y-'
      WRITE(16,'(1X,A30,A27)') KFID, &
     &     ';-DEPENDENT VARIABLE LABEL-'
      WRITE(16,'(1X,I1,A29,A39)') 2,'                             ', &
     &     ';-PROC CODE- 0:RAW 1:AVG 2:SM. 3:AVG+SM'
      WRITE(16,'(1X,I11,A19,A12)') NRXMAX,'                   ', &
     &     ';-# OF X PTS-:'
      WRITE(16,'(1X,I11,A19,A35)') NTXMAX,'                   ', &
     &     ';-# OF Y PTS-  X, F(X) DATA FOLLOW:'

      WRITE(16,'(1X,1P6E13.6)') (GR(NRX),NRX=1,NRXMAX)
      WRITE(16,'(1X,1P6E13.6)') (GT(NTX),NTX=1,NTXMAX)
      WRITE(16,'(1X,1P6E13.6)') ((GF(NRX,NTX),NRX=1,NRXMAX),NTX=1,NTXMAX)

      WRITE(16,'(A52)')  ';----END-OF-DATA-----------------COMMENTS:-----------'

      CLOSE(16)
      RETURN

   10 WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
      RETURN
      END SUBROUTINE TRXW2D

!     *****

      SUBROUTINE GET_DATE(CDATE)

      IMPLICIT NONE
      CHARACTER(LEN=9),INTENT(OUT):: CDATE
      INTEGER      :: NDD, NDM, NDY, NTIH, NTIM, NTIS
      CHARACTER(LEN=2):: CDD, CDY
      CHARACTER(LEN=3):: CDM
      CHARACTER(LEN=3),DIMENSION(12):: CDATA = &
     &     (/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/)


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

      RETURN
      END SUBROUTINE GET_DATE

!   *** setup variable name strings ***

      SUBROUTINE tr_setup_kv

      USE TRCOMM
      IMPLICIT NONE

      KVT( 1) = 'ANS0(1)   '
      KVT( 2) = 'ANS0(2)   '
      KVT( 3) = 'ANS0(3)   '
      KVT( 4) = 'ANS0(4)   '
      KVT( 5) = 'ANSAV(1)  '
      KVT( 6) = 'ANSAV(2)  '
      KVT( 7) = 'ANSAV(3)  '
      KVT( 8) = 'ANSAV(4)  '

      KVT( 9) = 'TS0(1)    '
      KVT(10) = 'TS0(2)    '
      KVT(11) = 'TS0(3)    '
      KVT(12) = 'TS0(4)    '
      KVT(13) = 'TSAV(1)   '
      KVT(14) = 'TSAV(2)   '
      KVT(15) = 'TSAV(3)   '
      KVT(16) = 'TSAV(4)   '

      KVT(17) = 'WST(1)    '
      KVT(18) = 'WST(2)    '
      KVT(19) = 'WST(3)    '
      KVT(20) = 'WST(4)    '

      KVT(21) = 'ANF0(1)   '
      KVT(22) = 'ANF0(2)   '
      KVT(23) = 'ANFAV(1)  '
      KVT(24) = 'ANFAV(2)  '
      KVT(25) = 'TF0(1)    '
      KVT(26) = 'TF0(2)    '
      KVT(27) = 'TFAV(1)   '
      KVT(28) = 'TFAV(2)   '

      KVT(29) = 'WFT(1)    '
      KVT(30) = 'WFT(2)    '
      KVT(31) = 'WBULKT    '
      KVT(32) = 'WTAILT    '
      KVT(33) = 'WPT       '

      KVT(34) = 'AJT       '
      KVT(35) = 'AJOHT     '
      KVT(36) = 'AJNBT     '
      KVT(37) = 'AJRFT     '
      KVT(38) = 'AJBST     '

      KVT(39) = 'PINT      '
      KVT(40) = 'POHT      '
      KVT(41) = 'PNBT      '
      KVT(42) = 'PRFT(1)   '
      KVT(43) = 'PRFT(2)   '
      KVT(44) = 'PRFT(3)   '
      KVT(45) = 'PRFT(4)   '
      KVT(46) = 'PNFT      '

      KVT(47) = 'PBINT     '
      KVT(48) = 'PBCLT(1)  '
      KVT(49) = 'PBCLT(2)  '
      KVT(50) = 'PBCLT(3)  '
      KVT(51) = 'PBCLT(4)  '
      KVT(52) = 'PFINT     '
      KVT(53) = 'PFCLT(1)  '
      KVT(54) = 'PFCLT(2)  '
      KVT(55) = 'PFCLT(3)  '
      KVT(56) = 'PFCLT(4)  '

      KVT(57) = 'POUT      '
      KVT(58) = 'PCXT      '
      KVT(59) = 'PIET      '
      KVT(60) = 'PRSUMT    '
      KVT(61) = 'PLT(1)    '
      KVT(62) = 'PLT(2)    '
      KVT(63) = 'PLT(3)    '
      KVT(64) = 'PLT(4)    '

      KVT(65) = 'SINT      '
      KVT(66) = 'SIET      '
      KVT(67) = 'SNBT      '
      KVT(68) = 'SNFT      '
      KVT(69) = 'SOUT      '
      KVT(70) = 'SLT(1)    '
      KVT(71) = 'SLT(2)    '
      KVT(72) = 'SLT(3)    '
      KVT(73) = 'SLT(4)    '

      KVT(74) = 'VLOOP     '
      KVT(75) = 'ALI       '
      KVT(76) = 'RQ1       '
      KVT(77) = 'Q0        '

      KVT(78) = 'WPDOT     '
      KVT(79) = 'TAUE1     '
      KVT(80) = 'TAUE2     '
      KVT(81) = 'TAUE89    '

      KVT(82) = 'BETAP0    '
      KVT(83) = 'BETAPA    '
      KVT(84) = 'BETA0     '
      KVT(85) = 'BETAA     '

      KVT(86) = 'ZEFF0     '
      KVT(87) = 'QF        '
      KVT(88) = 'RIP       '
!
      KVT(89) = 'PEXT(1)   '
      KVT(90) = 'PEXT(2)   '
      KVT(91) = 'PRFVT(1,1)' ! ECH  to electron
      KVT(92) = 'PRFVT(2,1)' ! ECH  to ions
      KVT(93) = 'PRFVT(1,2)' ! LH   to electron
      KVT(94) = 'PRFVT(2,2)' ! LH   to ions
      KVT(95) = 'PRFVT(1,3)' ! ICRH to electron
      KVT(96) = 'PRFVT(2,3)' ! ICRH to ions

      KVT(97) = 'RR        '
      KVT(98) = 'RA        '
      KVT(99) = 'BB        '
      KVT(100)= 'RKAP      '
      KVT(101)= 'AJTTOR    '

      KVT(102)= 'TAUE98    '
      KVT(103)= 'H98Y2     '
      KVT(104)= 'ANLAV(1)  '
      KVT(105)= 'ANLAV(2)  '
      KVT(106)= 'ANLAV(3)  '
      KVT(107)= 'ANLAV(4)  '

      KVT(108)= 'PRBT      '
      KVT(109)= 'PRCT      '
      KVT(110)= 'PRLT      '

!     *** FOR 3D ***

      KVRT( 1) = 'RT(1)     '
      KVRT( 2) = 'RT(2)     '
      KVRT( 3) = 'RT(3)     '
      KVRT( 4) = 'RT(4)     '

      KVRT( 5) = 'RN(1)     '
      KVRT( 6) = 'RN(2)     '
      KVRT( 7) = 'RN(3)     '
      KVRT( 8) = 'RN(4)     '

      KVRT( 9) = 'AJ        '
      KVRT(10) = 'AJOH      '
      KVRT(11) = 'AJNB      '
      KVRT(12) = 'AJRF      '
      KVRT(13) = 'AJBS      '

      KVRT(14) = 'PTOT      '
      KVRT(15) = 'POH       '
      KVRT(16) = 'PNB       '
      KVRT(17) = 'PNF       '
      KVRT(18) = 'PRF(1)    '
      KVRT(19) = 'PRF(2)    '
      KVRT(20) = 'PRF(3)    '
      KVRT(21) = 'PRF(4)    '
      KVRT(22) = 'PRL       '
      KVRT(23) = 'PCX       '
      KVRT(24) = 'PIE       '
      KVRT(25) = 'PEX(1)    '
      KVRT(26) = 'PEX(2)    '
      KVRT(27) = 'QP        '
      KVRT(28) = 'EZOH      '
      KVRT(29) = 'BETA      '
      KVRT(30) = 'BETAP     '
      KVRT(31) = 'EZOH*2PIRR'
      KVRT(32) = 'ETA       '
      KVRT(33) = 'ZEFF      '
      KVRT(34) = 'AK(1)     '
      KVRT(35) = 'AK(2)     '

      KVRT(36) = 'PRFV(1,1) '
      KVRT(37) = 'PRFV(1,2) '
      KVRT(38) = 'PRFV(1,3) '
      KVRT(39) = 'PRFV(2,1) '
      KVRT(40) = 'PRFV(2,2) '
      KVRT(41) = 'PRFV(2,3) '

      KVRT(42) = 'AJRFV(1)  '
      KVRT(43) = 'AJRFV(2)  '
      KVRT(44) = 'AJRFV(3)  '

      KVRT(45) = 'RW(1+2)   '
      KVRT(46) = 'ANC+ANFE  '
      KVRT(47) = 'BP        '
      KVRT(48) = 'RPSI      '

      KVRT(49) = 'RMJRHO    '
      KVRT(50) = 'RMNRHO    '
      KVRT(51) = 'F0D       '
      KVRT(52) = 'RKPRHO    '
      KVRT(53) = 'DELTAR    '
      KVRT(54) = 'AR1RHO    '
      KVRT(55) = 'AR2RHO    '
      KVRT(56) = 'AKDW(1)   '
      KVRT(57) = 'AKDW(2)   '
      KVRT(58) = 'RN*RT(1)  '
      KVRT(59) = 'RN*RT(2)  '

      KVRT(60) = 'VTOR      '
      KVRT(61) = 'VPOL      '

      KVRT(62) = 'S-ALPHA   '
      KVRT(63) = 'ER        '
      KVRT(64) = 'S         '
      KVRT(65) = 'ALPHA     '
      KVRT(66) = 'TRCOFS    '
      KVRT(67) = '2PI/QP    '

      RETURN
    END SUBROUTINE tr_setup_kv
