C     $Id$
C
C     ***********************************************************
C
C           MAIN ROUTINE FOR TRANSPORT CALCULATION
C
C     ***********************************************************
C
      SUBROUTINE TRLOOP
C
      INCLUDE 'trcomm.h'
      PARAMETER(NURM=51,NUTM=61)
      COMMON /TRBCH1/ RAD(NURM),FRFHE(NURM,NUTM),FRFHI(NURM,NUTM)
      COMMON /TRBCH3/ FUT(NUTM),NUFMAX,NTXMAX
      COMMON /PRETREAT2/ NTAMAX
      DIMENSION DERIVR(NURM),URFHEPH(4,NURM),URFHIPH(4,NURM)
      DIMENSION DERIVT(NUTM),UFUTE(4,NUTM),UFUTI(4,NUTM)
      DIMENSION FRFHET(NUTM),FRFHIT(NUTM),FRFHER(NURM),FRFHIR(NURM)
C
      DIMENSION XX(MLM),YY(2,NRM)
C
      IF(MDLUF.NE.1) THEN
         NT=0
         RIP=RIPS
         IF(NTMAX.NE.0) DIP=(RIPE-RIPS)/DBLE(NTMAX)
      ELSE
         IF(NTMAX.GT.NTAMAX) NTMAX=NTAMAX
      ENDIF
      ICHCK=0
C
C     *****
C
      IF(MDLUF.EQ.3) THEN
         DO NUR=1,NURM
            DERIVR(NUR)=0.D0
         ENDDO
         DO NUT=1,NUTM
            DERIVT(NUT)=0.D0
         ENDDO
         DO NUF=1,NUFMAX
            DO NTX=1,NTXMAX
               FRFHET(NTX)=FRFHE(NUF,NTX)
               FRFHIT(NTX)=FRFHI(NUF,NTX)
            ENDDO
            CALL SPL1D(FUT,FRFHET,DERIVT,UFUTE,NTXMAX,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TRLOOP: SPL1D FUTE: IERR=',IERR
            CALL SPL1D(FUT,FRFHIT,DERIVT,UFUTI,NTXMAX,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TRLOOP: SPL1D FUTI: IERR=',IERR
C     
            CALL SPL1DF(T,PRFHEA,FUT,UFUTE,NTXMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRLOOP: SPL1DF PRFHEA: IERR=',IERR
            CALL SPL1DF(T,PRFHIA,FUT,UFUTI,NTXMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRLOOP: SPL1DF PRFHIA: IERR=',IERR
            FRFHER(NUF)=PRFHEA
            FRFHIR(NUF)=PRFHIA
         ENDDO
C
         CALL SPL1D(RAD,FRFHER,DERIVR,URFHEPH,NUFMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRLOOP: SPL1D FRFHER: IERR=',IERR
         CALL SPL1D(RAD,FRFHIR,DERIVR,URFHIPH,NUFMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRLOOP: SPL1D FRFHIR: IERR=',IERR
         DO NR=1,NRMAX
            RMNOW=RM(NR)
            CALL SPL1DF(RMNOW,PRFHE,RAD,URFHEPH,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRLOOP: SPL1DF PRFHE: IERR=',IERR
            CALL SPL1DF(RMNOW,PRFHI,RAD,URFHIPH,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRLOOP: SPL1DF PRFHI: IERR=',IERR
            PEX(NR,1)=PRFHE
            PEX(NR,2)=PRFHI
            PEX(NR,3)=0.D0
            PEX(NR,4)=0.D0
         ENDDO
      ENDIF
C
C     *****
C
      CALL TRCALC(IERR)
      IF(IERR.NE.0) RETURN
      CALL TRGLOB
      CALL TRSNAP
      IF(NGT.EQ.0) THEN
         CALL TRATOT
         CALL TRATOTN
      ENDIF
      IF(NGR.EQ.0) CALL TRATOG
      IF(NT.GE.NTMAX) GOTO 9000
C
 1000 L=0
C     /* Making New Variables */
      CALL TRATOX
C
C     /* Stored Variables for Convergence Check */
      DO J=1,NEQMAX
      DO NR=1,NRMAX
         XX(NEQMAX*(NR-1)+J)=XV(J,NR)
      ENDDO
      ENDDO
      DO J=1,NFM
      DO NR=1,NRMAX
         YY(J,NR)=YV(J,NR)
      ENDDO
      ENDDO
C
 2000 CONTINUE
C
C     /* Matrix Producer */
      CALL TRMTRX(NEQRMAX)
C
C     /* Matrix Solver */
      MWRMAX=4*NEQRMAX-1
      CALL BANDRD(AX,X,NEQRMAX*NRMAX,MWRMAX,MWM,IERR)
      IF(IERR.EQ.30000) THEN
         WRITE(6,*) 'XX ERROR IN TRLOOP : MATRIX AA IS SINGULAR ',
     &              ' AT ',NT,' STEP.'
         GOTO 9000
      ENDIF
C     
      DO J=1,NFM
      DO NR=1,NRMAX
         Y(J,NR) = Y(J,NR)/AY(J,NR)
      ENDDO
      ENDDO
C
C     /* Convergence check */
      DO I=1,NEQRMAX*NRMAX
         IF (ABS(X(I)-XX(I)).GT.EPSLTR*ABS(X(I))) GOTO 3000
      ENDDO
      DO J=1,NFM
      DO NR=1,NRMAX
         IF (ABS(Y(J,NR)-YY(J,NR)).GT.EPSLTR*ABS(Y(J,NR))) GOTO 3000
      ENDDO
      ENDDO
C
C     *****
C
      IF(MDLUF.EQ.3) THEN
C         write(6,*) "T=",T
         DO NUF=1,NUFMAX
            DO NTX=1,NTXMAX
               FRFHET(NTX)=FRFHE(NUF,NTX)
               FRFHIT(NTX)=FRFHI(NUF,NTX)
            ENDDO
            CALL SPL1D(FUT,FRFHET,DERIVT,UFUTE,NTXMAX,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TRLOOP: SPL1D FUTE: IERR=',IERR
            CALL SPL1D(FUT,FRFHIT,DERIVT,UFUTI,NTXMAX,0,IERR)
            IF(IERR.NE.0) WRITE(6,*) 'XX TRLOOP: SPL1D FUTI: IERR=',IERR
C     
            CALL SPL1DF(T,PRFHEA,FUT,UFUTE,NTXMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRLOOP: SPL1DF PRFHEA: IERR=',IERR
            CALL SPL1DF(T,PRFHIA,FUT,UFUTI,NTXMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRLOOP: SPL1DF PRFHIA: IERR=',IERR
            FRFHER(NUF)=PRFHEA
            FRFHIR(NUF)=PRFHIA
         ENDDO
C
         CALL SPL1D(RAD,FRFHER,DERIVR,URFHEPH,NUFMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRLOOP: SPL1D FRFHER: IERR=',IERR
         CALL SPL1D(RAD,FRFHIR,DERIVR,URFHIPH,NUFMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TRLOOP: SPL1D FRFHIR: IERR=',IERR
         DO NR=1,NRMAX
            RMNOW=RM(NR)
            CALL SPL1DF(RMNOW,PRFHE,RAD,URFHEPH,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRLOOP: SPL1DF PRFHE: IERR=',IERR
            CALL SPL1DF(RMNOW,PRFHI,RAD,URFHIPH,NUFMAX,IERR)
            IF(IERR.NE.0)
     &           WRITE(6,*) 'XX TRLOOP: SPL1DF PRFHI: IERR=',IERR
            PEX(NR,1)=PRFHE
            PEX(NR,2)=PRFHI
            PEX(NR,3)=0.D0
            PEX(NR,4)=0.D0
         ENDDO
      ENDIF
C
C     *****
C
      GOTO 4000
C
 3000 L=L+1
      IF(L.GE.LMAXTR) GOTO 4000
C
C     /* Stored Variables for Convergence Check */
      DO I=1,NEQRMAX*NRMAX
         XX(I) = X(I)
      ENDDO
      DO J=1,NFM
      DO NR=1,NRMAX
         YY(J,NR) = Y(J,NR)
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            NSTN=NST(NEQ)
            IF(NSVN.EQ.0) THEN
               IF(NSTN.EQ.0) THEN
                  BP(NR) = XV(NEQ,NR)
               ELSE
                  BP(NR) = 0.5D0*(XV(NEQ,NR)+X(NEQRMAX*(NR-1)+NSTN))
               ENDIF
            ELSEIF(NSVN.EQ.1) THEN
               IF(NSSN.EQ.1.AND.MDLEQE.EQ.0) THEN
                  RN(NR,NSSN) = 0.D0
                  DO NEQ1=1,NEQMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     NSTN1=NST(NEQ1)
                     IF(NSVN1.EQ.1.AND.NSSN1.NE.1) THEN
                        IF(NSTN1.EQ.0) THEN
                           RN(NR,NSSN) = RN(NR,NSSN)+PZ(NSSN1)
     &                                  *XV(NEQ1,NR)
                        ELSE
                           RN(NR,NSSN) = RN(NR,NSSN)+PZ(NSSN1)
     &                      *0.5D0*(XV(NEQ1,NR)+X(NEQRMAX*(NR-1)+NSTN1))
                        ENDIF
                     ENDIF
                  ENDDO
                  RN(NR,NSSN) = RN(NR,NSSN)
     &                         +PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
                  IF(MDLEQE.EQ.0) XV(NEQ,NR) = RN(NR,NSSN)
               ELSEIF(NSSN.EQ.1.AND.MDLEQE.EQ.1) THEN
                  RN(NR,NSSN) = 0.5D0*(XV(NEQ,NR)+X(NEQRMAX*(NR-1)+NEQ))
               ELSE
                  IF(NSTN.EQ.0) THEN
                     RN(NR,NSSN) = XV(NEQ,NR)
                  ELSE
                     RN(NR,NSSN) = 0.5D0*(XV(NEQ,NR)
     &                                  +  X(NEQRMAX*(NR-1)+NSTN))
                  ENDIF
               ENDIF
            ELSEIF(NSVN.EQ.2) THEN
               IF(NSTN.EQ.0) THEN
                  IF(NSSN.NE.NSM) THEN
                     RT(NR,NSSN) = XV(NEQ,NR)/RN(NR,NSSN)
                  ELSE
                     IF(RN(NR,NSM).LT.1.D-70) THEN
                        RT(NR,NSM) = 0.D0
                     ELSE
                        RT(NR,NSM) = XV(NEQ,NR)/RN(NR,NSM)
                     ENDIF
                  ENDIF
               ELSE
                  IF(NSSN.NE.NSM) THEN
                     RT(NR,NSSN) = 0.5D0*(XV(NEQ,NR)
     &                    +X(NEQRMAX*(NR-1)+NSTN))/RN(NR,NSSN)
                  ELSE
                     IF(RN(NR,NSM).LT.1.D-70) THEN
                        RT(NR,NSM) = 0.D0
                     ELSE
                        RT(NR,NSM) = 0.5D0*(XV(NEQ,NR)
     &                       +X(NEQRMAX*(NR-1)+NSTN))/RN(NR,NSM)
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF(NSVN.EQ.3) THEN
               IF(NSTN.EQ.0) THEN
                  RU(NR,NSSN) = XV(NEQ,NR)/(PA(NSSN)*AMM*RN(NR,NSSN))
               ELSE
                  RU(NR,NSSN) = 0.5D0*(XV(NEQ,NR)
     &                 +  X(NEQRMAX*(NR-1)+NSTN))
     &                 /(PA(NSSN)*AMM*RN(NR,NSSN))
               ENDIF
            ENDIF
         ENDDO
         ANNU(NR)=RN(NR,7)+RN(NR,8)
      ENDDO
C
      CALL TRCHCK(ICHCK)
      IF(ICHCK.EQ.1) GOTO 4000
C
      CALL TRCALC(IERR)
      IF(IERR.NE.0) RETURN
      GOTO 2000
C
 4000 NT=NT+1
      T=T+DT
      VSEC=VSEC+VLOOP*DT
      IF(Q0.LT.1.D0) TST=TST+DT
      IF(MDLUF.EQ.1) THEN
         RIP=RIPU(NT)
      ELSE
         RIP=RIP+DIP
      ENDIF
C      write(6,'(A,1P3E12.5)') "RIP,RIPE,DIP=",RIP,RIPE,DIP
C
C     /* Making new XV(NEQ,NR) and YV(NF,NR) */
      DO NEQ=1,NEQMAX
         NSTN=NST(NEQ)
         IF(NSTN.NE.0) THEN
            DO NR=1,NRMAX
               XV(NEQ,NR) = X(NEQRMAX*(NR-1)+NSTN)
            ENDDO
         ENDIF
      ENDDO
C
      DO NF=1,NFM
      DO NR=1,NRMAX
         YV(NF,NR) = Y(NF,NR)
      ENDDO
      ENDDO
C
C     /* Making New Physical Variables */
      CALL TRXTOA
      IF(MDLUF.EQ.1) CALL TR_UFREAD
C
      IF(ICHCK.EQ.0) CALL TRCHCK(ICHCK)
      IF(ICHCK.EQ.1) THEN
         CALL TRGLOB
         CALL TRATOT
         CALL TRATOTN
         CALL TRATOG
         GOTO 9000
      ENDIF
C
C     /* Sawtooth Oscillation */
      IF(TST+0.5D0*DT.GT.TPRST) THEN
         CALL TRSAWT
         TST=0.D0
      ENDIF
C
      CALL TRCALC(IERR)
      IF(IERR.NE.0) RETURN
C
      IDGLOB=0
      IF(MOD(NT,NTSTEP).EQ.0) THEN
         IF(IDGLOB.EQ.0) CALL TRGLOB
         IDGLOB=1
         CALL TRSNAP
      ENDIF
      IF(MOD(NT,NGTSTP).EQ.0) THEN
         IF(IDGLOB.EQ.0) CALL TRGLOB
         IDGLOB=1
         CALL TRATOT
         CALL TRATOTN
      ENDIF
      IF(MOD(NT,NGRSTP).EQ.0) THEN
         IF(IDGLOB.EQ.0) CALL TRGLOB
         IDGLOB=1
         CALL TRATOG
      ENDIF
      IF(MODELG.EQ.3.AND.MOD(NT,NTEQIT).EQ.0) THEN
C         CALL TRCONV(L,IERR)
C         WRITE(6,*) "L=",L
         CALL TRSETG
         IF(IERR.NE.0) RETURN
      ENDIF
      IF(NT.LT.NTMAX) GOTO 1000
C
 9000 IF(MDLUF.NE.1) RIPS=RIPE
      RETURN
      END
C
C     ***********************************************************
C
C           COMPUTE MATRIX ELEMENTS
C
C     ***********************************************************
C
      SUBROUTINE TRMTRX(NEQRMAX)
C
      INCLUDE 'trcomm.h'
C
      COMMON /TRLCL1/ A(NVM,NVM,NRM),B(NVM,NVM,NRM),C(NVM,NVM,NRM)
      COMMON /TRLCL2/ D(NVM,NRM)
      COMMON /TRLCL3/ RD(NEQM,NRM)
      COMMON /TRLCL4/ PPA(NEQM,NRM),PPB(NEQM,NRM),PPC(NEQM,NRM)
C     
      IF(MDLCD.EQ.0) THEN
         BPS= AMYU0*RIP*1.D6/(2.D0*PI*RA*RKAPS)
      ELSE
         NEQ=1
         NSVN=NSS(NEQ)
         IF(NSVN.EQ.0) THEN
            BPA=XV(NEQ,NRMAX)
         ELSE
            BPA=0.D0
         ENDIF
         RLP=RA*(LOG(8.D0*RR/RA)-2.D0)
         BPS=BPA-AMYU0*DT*EZOH(NRMAX)/RLP
      ENDIF
      IF(MDLUF.EQ.1) THEN
         IF(NT.EQ.0) THEN
            RKAP=RKAPU(1)
            RKAPS=SQRT(RKAP)
            BPS=AMYU0*RIPU(1)*1.D6/(2.D0*PI*RA*RKAPS)
         ELSE
            RKAP=RKAPU(NT)
            RKAPS=SQRT(RKAP)
            BPS=AMYU0*RIPU(NT)*1.D6/(2.D0*PI*RA*RKAPS)
         ENDIF
      ENDIF
C      IF(MODELG.EQ.3) THEN
C         BPS=BPSEQ*RIP/RIPEQ
C      ENDIF
C
      COEF = AEE**4*5.D0*1.D20/(SQRT(2.D0*PI)*PI*AEPS0**2)
C
      DO NR=1,NRMAX
      DO NW=1,NEQMAX
      DO NV=1,NEQMAX
         A(NV,NW,NR)=0.D0
         B(NV,NW,NR)=0.D0
         C(NV,NW,NR)=0.D0
      ENDDO
      ENDDO
      ENDDO
      DO NR=1,NRMAX
      DO NEQ=1,NEQMAX
         D(NEQ,NR)=0.D0
         PPA(NEQ,NR)=0.D0
         PPB(NEQ,NR)=0.D0
         PPC(NEQ,NR)=0.D0
      ENDDO
      ENDDO
      DO NW=1,NEQMAX
      DO NV=1,NEQMAX
      DO NX=1,4
      DO NSW=1,3
         IF (NX.LT.3) THEN
            VV(NV,NW,NX,NSW)=0.D0
            DD(NV,NW,NX,NSW)=0.D0
            VI(NV,NW,NX,NSW)=0.D0
            DI(NV,NW,NX,NSW)=0.D0
         ELSE 
            VV(NV,NW,NX,NSW)=0.D0
            DD(NV,NW,NX,NSW)=0.D0
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      MWMAX=4*NEQMAX-1
      DO MV=1,NEQMAX*NRMAX
      DO MW=1,MWMAX
         AX(MW,MV) = 0.D0
      ENDDO
      ENDDO
C
C      FADV=0.5D0
      FADV=1.0D0
C
      PRV=(1.D0-FADV)*DT
      ADV=FADV*DT
C
C          /----------\
C    ***   |   NR=1   |   ***
C          \----------/
C
      NR=1
      NSW=1
      CALL TR_COEF_DECIDE(NR,NSW,DV53)
C
      DO NV=1,NEQMAX
      DO NW=1,NEQMAX
         A(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)+DI(NV,NW,1,NSW)
         B(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)-DI(NV,NW,1,NSW)
     &                -0.5D0*VI(NV,NW,2,NSW)-DI(NV,NW,2,NSW)
         C(NV,NW,NR) =-0.5D0*VI(NV,NW,2,NSW)+DI(NV,NW,2,NSW)
      ENDDO
      ENDDO
C
      DO NEQ =1,NEQMAX
      DO NEQ1=1,NEQMAX
         NSVN1=NSV(NEQ1)
         NSVN =NSV(NEQ )
         IF(NSVN1.EQ.2.AND.NSVN.EQ.2) THEN
            NSSN1=NSS(NEQ1)
            NSSN =NSS(NEQ )
            IF(NSSN1.NE.NSSN) THEN
               C1=COEF/((RTM(NSSN)+RTM(NSSN1))**1.5D0*AMZ(NSSN)
     &           *AMZ(NSSN1))*DV53*1.5D0
               B(NEQ,NEQ, NR)=B(NEQ,NEQ, NR)-RN(NR,NSSN1)*C1
               B(NEQ,NEQ1,NR)=B(NEQ,NEQ1,NR)+RN(NR,NSSN )*C1
            ENDIF
         ENDIF
      ENDDO
      ENDDO
C
      CALL TR_IONIZATION(NR)
      CALL TR_CHARGE_EXCHANGE(NR)
C
C     ***** RHS Vector *****
C
      DO NEQ=1,NEQMAX
         X(NEQMAX*(NR-1)+NEQ) = RD(NEQ,NR)*XV(NEQ,NR)+DT*D(NEQ,NR)
      ENDDO
C
      DO NW=1,NEQMAX
      DO NV=1,NEQMAX
         X(NEQMAX*(NR-1)+NV) = X(NEQMAX*(NR-1)+NV)
     &                        +PRV*(
     &                              +B(NV,NW,NR)*XV(NW,NR  )
     &                              +C(NV,NW,NR)*XV(NW,NR+1))
      ENDDO
      ENDDO
C
C     ***** Evolution of fast ion components *****
C
      Y(1,NR)=(1.5D0-PRV/TAUB(NR))*YV(1,NR)+PNB(NR)*DT/(RKEV*1.D20)
      Y(2,NR)=(1.5D0-PRV/TAUF(NR))*YV(2,NR)+PNF(NR)*DT/(RKEV*1.D20)
      AY(1,NR)=1.5D0+ADV/TAUB(NR)
      AY(2,NR)=1.5D0+ADV/TAUF(NR)
C
C          /---------------------\
C    ***   |   NR=2 to NRMAX-1   |   ***
C          \---------------------/
C
      NSW=2
      DO NR=2,NRMAX-1
         CALL TR_COEF_DECIDE(NR,NSW,DV53)
C     
         DO NV=1,NEQMAX
         DO NW=1,NEQMAX
            A(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)+DI(NV,NW,1,NSW)
            B(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)-DI(NV,NW,1,NSW)
     &                   -0.5D0*VI(NV,NW,2,NSW)-DI(NV,NW,2,NSW)
            C(NV,NW,NR) =-0.5D0*VI(NV,NW,2,NSW)+DI(NV,NW,2,NSW)
         ENDDO
         ENDDO
C
         DO NEQ =1,NEQMAX
         DO NEQ1=1,NEQMAX
            NSVN1=NSV(NEQ1)
            NSVN =NSV(NEQ )
            IF(NSVN1.EQ.2.AND.NSVN.EQ.2) THEN
               NSSN1=NSS(NEQ1)
               NSSN =NSS(NEQ )
               IF(NSSN1.NE.NSSN) THEN
                  C1=COEF/((RTM(NSSN)+RTM(NSSN1))**1.5D0*AMZ(NSSN)
     &                 *AMZ(NSSN1))*DV53*1.5D0
                  B(NEQ,NEQ, NR)=B(NEQ,NEQ, NR)-RN(NR,NSSN1)*C1
                  B(NEQ,NEQ1,NR)=B(NEQ,NEQ1,NR)+RN(NR,NSSN )*C1
               ENDIF
            ENDIF
         ENDDO
         ENDDO
C
         CALL TR_IONIZATION(NR)
         CALL TR_CHARGE_EXCHANGE(NR)
C
C     ***** RHS Vector *****
C
         DO NEQ=1,NEQMAX
            X(NEQMAX*(NR-1)+NEQ) = RD(NEQ,NR)*XV(NEQ,NR)+DT*D(NEQ,NR)
         ENDDO
C
         DO NW=1,NEQMAX
         DO NV=1,NEQMAX
            X(NEQMAX*(NR-1)+NV) = X(NEQMAX*(NR-1)+NV)
     &                           +PRV*(A(NV,NW,NR)*XV(NW,NR-1)
     &                                +B(NV,NW,NR)*XV(NW,NR  )
     &                                +C(NV,NW,NR)*XV(NW,NR+1))
         ENDDO
         ENDDO
C
C     ***** Evolution of fast ion components *****
C
         Y(1,NR)=(1.5D0-PRV/TAUB(NR))*YV(1,NR)+PNB(NR)*DT/(RKEV*1.D20)
         Y(2,NR)=(1.5D0-PRV/TAUF(NR))*YV(2,NR)+PNF(NR)*DT/(RKEV*1.D20)
         AY(1,NR)=1.5D0+ADV/TAUB(NR)
         AY(2,NR)=1.5D0+ADV/TAUF(NR)
C 
      ENDDO
C
C          /--------------\
C    ***   |   NR=NRMAX   |   ***
C          \--------------/
C
      NR=NRMAX
      NSW=3
      CALL TR_COEF_DECIDE(NR,NSW,DV53)
C     
      DO NV=1,NEQMAX
      DO NW=1,NEQMAX
         A(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)+DI(NV,NW,1,NSW)
         B(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)-DI(NV,NW,1,NSW)
     &                -0.0D0*VI(NV,NW,2,NSW)-DI(NV,NW,2,NSW)
     &                                      -DI(NV,NW,2,NSW)
         C(NV,NW,NR) =-0.0D0
      ENDDO
      ENDDO
C
      DO NEQ =1,NEQMAX
      DO NEQ1=1,NEQMAX
         NSVN1=NSV(NEQ1)
         NSVN =NSV(NEQ )
         IF(NSVN1.EQ.2.AND.NSVN.EQ.2) THEN
            NSSN1=NSS(NEQ1)
            NSSN =NSS(NEQ )
            IF(NSSN1.NE.NSSN) THEN
               C1=COEF/((RTM(NSSN)+RTM(NSSN1))**1.5D0*AMZ(NSSN)
     &              *AMZ(NSSN1))*DV53*1.5D0
               B(NEQ,NEQ, NR)=B(NEQ,NEQ, NR)-RN(NR,NSSN1)*C1
               B(NEQ,NEQ1,NR)=B(NEQ,NEQ1,NR)+RN(NR,NSSN )*C1
            ENDIF
         ENDIF
      ENDDO
      ENDDO
C
      CALL TR_IONIZATION(NR)
      CALL TR_CHARGE_EXCHANGE(NR)
C     
C     ***** RHS Vector *****
C
      DO NEQ=1,NEQMAX
         X(NEQMAX*(NR-1)+NEQ) = RD(NEQ,NR)*XV(NEQ,NR)+DT*D(NEQ,NR)
      ENDDO
C
      DO NW=1,NEQMAX
      DO NV=1,NEQMAX
         X(NEQMAX*(NR-1)+NV) = X(NEQMAX*(NR-1)+NV)
     &                        +PRV*(A(NV,NW,NR)*XV(NW,NR-1)
     &                             +B(NV,NW,NR)*XV(NW,NR  )
     &                                                     )
      ENDDO
      ENDDO
C
C     ***** Evolution of fast ion components *****
C
      Y(1,NR)=(1.5D0-PRV/TAUB(NR))*YV(1,NR)+PNB(NR)*DT/(RKEV*1.D20)
      Y(2,NR)=(1.5D0-PRV/TAUF(NR))*YV(2,NR)+PNF(NR)*DT/(RKEV*1.D20)
      AY(1,NR)=1.5D0+ADV/TAUB(NR)
      AY(2,NR)=1.5D0+ADV/TAUF(NR)
C
C     ***** Band Matrix *****
C
      DO NR=1,NRMAX
         CALL TR_BAND_GEN(NR,NEQRMAX,ADV)
      ENDDO
C
C     ***** RHS Vector Reduce *****
C
      DO NR=1,NRMAX
      DO NEQ=1,NEQMAX
         NSTN=NST(NEQ)
         IF(NSTN.NE.0) THEN
            X(NEQRMAX*(NR-1)+NSTN) = X(NEQMAX*(NR-1)+NEQ)
     &                         +PPA(NEQ,NR)+PPB(NEQ,NR)+PPC(NEQ,NR)
         ENDIF
      ENDDO
      ENDDO
C
C     ***** Surface Boundary Condition for Bp *****
C
      NEQ=1
      NSVN=NSV(NEQ)
      IF(NSVN.EQ.0) THEN
         MV=NEQRMAX*(NRMAX-1)+NEQ
         DO MW=1,MWMAX
            AX(MW,MV)=0.D0
         ENDDO
         AX(2*NEQRMAX,MV)=RD(NEQ,NRMAX)
         X(MV)=BPS
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
C           BAND MATRIX GENERATOR
C
C     ***********************************************************
C
      SUBROUTINE TR_BAND_GEN(NR,NEQRMAX,ADV)
C
      INCLUDE 'trcomm.h'
C
      COMMON /TRLCL1/ A(NVM,NVM,NRM),B(NVM,NVM,NRM),C(NVM,NVM,NRM)
      COMMON /TRLCL3/ RD(NEQM,NRM)
      COMMON /TRLCL4/ PPA(NEQM,NRM),PPB(NEQM,NRM),PPC(NEQM,NRM)
C
      DO NV=1,NEQMAX
      DO NW=1,NEQMAX
         A(NV,NW,NR) = -ADV*A(NV,NW,NR)
         B(NV,NW,NR) = -ADV*B(NV,NW,NR) 
         C(NV,NW,NR) = -ADV*C(NV,NW,NR)
      ENDDO
      ENDDO
C
      DO NEQ=1,NEQMAX
         B(NEQ,NEQ,NR)=B(NEQ,NEQ,NR)+RD(NEQ,NR)
      ENDDO
C
      IF(NR.EQ.1) THEN
         CALL TR_BANDREDUCE(B(1,1,NR),XV(1,NR  ),PPB(1,NR),NEQRMAX)
         CALL TR_BANDREDUCE(C(1,1,NR),XV(1,NR+1),PPC(1,NR),NEQRMAX)
      ELSEIF(NR.EQ.NRMAX) THEN
         CALL TR_BANDREDUCE(A(1,1,NR),XV(1,NR-1),PPA(1,NR),NEQRMAX)
         CALL TR_BANDREDUCE(B(1,1,NR),XV(1,NR  ),PPB(1,NR),NEQRMAX)
      ELSE
         CALL TR_BANDREDUCE(A(1,1,NR),XV(1,NR-1),PPA(1,NR),NEQRMAX)
         CALL TR_BANDREDUCE(B(1,1,NR),XV(1,NR  ),PPB(1,NR),NEQRMAX)
         CALL TR_BANDREDUCE(C(1,1,NR),XV(1,NR+1),PPC(1,NR),NEQRMAX)
      ENDIF
C
      DO NV=1,NEQRMAX
      DO NW=1,NEQRMAX
         AX(  NEQRMAX+NW-NV,NEQRMAX*(NR-1)+NV) = A(NV,NW,NR)
         AX(2*NEQRMAX+NW-NV,NEQRMAX*(NR-1)+NV) = B(NV,NW,NR)
         AX(3*NEQRMAX+NW-NV,NEQRMAX*(NR-1)+NV) = C(NV,NW,NR)
      ENDDO
      ENDDO
C
      RETURN
      END
C
C     ***********************************************************
C
C           BAND REDUCER
C
C     ***********************************************************
C
      SUBROUTINE TR_BANDREDUCE(A,XR,P,NEQRMAX)
C
      INCLUDE 'trcomm.h'
C     
      DIMENSION A(NEQM,NEQM),AA(NEQM,NEQM)
      DIMENSION XR(NEQM),P(NEQM),AL(NEQM,NEQM),AM(NEQM)
C
      DO NV=1,NEQMAX
         DO NW=1,NEQMAX
            AA(NV,NW)=A(NV,NW)
         ENDDO
         AM(NV)=0.D0
      ENDDO
C
      NBSIZE=NEQMAX
      LOOP=0
      DO NEQ=1,NEQM
         NNSN=NNS(NEQ)
         NNSN1=NNS(NEQ)
         IF(NNSN.EQ.0) THEN
            NEQRMAX=NEQMAX-LOOP
            GOTO 1000
         ENDIF
         IF(LOOP.EQ.0) THEN
            NNSOLD=NNSN
         ELSE
            IF(NNSN.GT.NNSOLD) THEN
               NNSN=NNSN-LOOP
               NNSOLD=NNSN
            ELSE
               STOP 'ERROR'
            ENDIF
         ENDIF
         LOOP=LOOP+1
C
C     /* Obtaining correction terms */
C
         DO NW=1,NEQMAX
            DO NV=1,NEQMAX
               IF(NW.EQ.NNSN1) AL(LOOP,NV)=AA(NV,NW)*XR(NNSN1)
            ENDDO
         ENDDO
C
C     /* Band reducer*/
C
         IF(NNSN.NE.NBSIZE) THEN
            NBSIZE=NBSIZE-1
            DO NV=1,NBSIZE+1
               DO NW=1,NBSIZE
                  IF(NW.GE.NNSN) A(NV,NW)=A(NV,NW+1)
               ENDDO
               A(NV,NBSIZE+1)=0.D0
            ENDDO
C
            DO NW=1,NBSIZE+1
               DO NV=1,NBSIZE
                  IF(NV.GE.NNSN) A(NV,NW)=A(NV+1,NW)
               ENDDO
               A(NBSIZE+1,NW)=0.D0
            ENDDO
         ELSE
            NBSIZE=NBSIZE-1
            DO NV=1,NBSIZE+1
               A(NV,NBSIZE+1)=0.D0
            ENDDO
            DO NW=1,NBSIZE+1
               A(NBSIZE+1,NW)=0.D0
            ENDDO
         ENDIF
      ENDDO
C
 1000 CONTINUE
C
C     /* Consummation */
C
      DO L=1,LOOP
         DO NEQ=1,NEQMAX
            AM(NEQ)=AM(NEQ)+AL(L,NEQ)
         ENDDO
      ENDDO
      DO NEQ=1,NEQMAX
         IF(NST(NEQ).NE.0) P(NEQ)=P(NEQ)-AM(NEQ)
      ENDDO
C
      RETURN
      END
C
C     ***********************************************************
C
C           CONVERSION FROM ANT TO XV
C
C     ***********************************************************
C
      SUBROUTINE TRATOX
C
      INCLUDE 'trcomm.h'
C
      DO NR=1,NRMAX
         DO NEQ=1,NEQMAX
            NSVN=NSV(NEQ)
            NSSN=NSS(NEQ)
            IF(NSVN.EQ.0) THEN
               XV(NEQ,NR) = BP(NR)
            ELSEIF(NSVN.EQ.1) THEN
               XV(NEQ,NR) = RN(NR,NSSN)
            ELSEIF(NSVN.EQ.2) THEN
               XV(NEQ,NR) = RN(NR,NSSN)*RT(NR,NSSN)
            ELSEIF(NSVN.EQ.3) THEN
               XV(NEQ,NR) = PA(NSSN)*AMM*RN(NR,NSSN)*RU(NR,NSSN)
            ENDIF
         ENDDO
         YV(  1,NR) = RW(NR,1)
         YV(  2,NR) = RW(NR,2)
      ENDDO
C
      RETURN
      END
C
C     ***********************************************************
C
C           CONVERT FROM XV TO ANT
C
C     ***********************************************************
C
      SUBROUTINE TRXTOA
C
      INCLUDE 'trcomm.h'
C
      DO NR=1,NRMAX
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            IF(NSVN.EQ.0) THEN
               BP(NR) = XV(NEQ,NR)
            ELSEIF(NSVN.EQ.1) THEN
               IF(NSSN.EQ.1.AND.MDLEQE.EQ.0) THEN
                  RN(NR,NSSN) = 0.D0
                  DO NEQ1=1,NEQMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1.AND.NSSN1.NE.1) THEN
                        RN(NR,NSSN) = RN(NR,NSSN)+PZ(NSSN1)*XV(NEQ1,NR)
                     ENDIF
                  ENDDO
                  RN(NR,NSSN) = RN(NR,NSSN)
     &                         +PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
C     I think the above expression is obsolete.
               ELSE
                  RN(NR,NSSN) = XV(NEQ,NR)
               ENDIF
            ELSEIF(NSVN.EQ.2) THEN
               IF(NSSN.NE.NSM) THEN
                  RT(NR,NSSN) = XV(NEQ,NR)/RN(NR,NSSN)
               ELSE
                  IF(RN(NR,NSM).LT.1.D-70) THEN
                     RT(NR,NSM) = 0.D0
                  ELSE
                     RT(NR,NSM) = XV(NEQ,NR)/RN(NR,NSM)
                  ENDIF
               ENDIF
            ELSEIF(NSVN.EQ.3) THEN
               RU(NR,NSSN) = XV(NEQ,NR)/(PA(NSSN)*AMM*RN(NR,NSSN))
            ENDIF
         ENDDO
         RW(NR,1) = YV(1,NR)
         RW(NR,2) = YV(2,NR)
         ANNU(NR) = RN(NR,7)+RN(NR,8)
C         if(NR.EQ.40.OR.NR.EQ.45) write(6,*) NT,NR,RN(NR,7),RN(NR,8)
C         write(6,*) ANNU(NR),RN(NR,7),RN(NR,8)
      ENDDO
C
      RETURN
      END
C
C     ***********************************************************
C
C           SAVE TRANSIENT DATA FOR GRAPHICS
C
C     ***********************************************************
C
      SUBROUTINE TRATOTN
C
      INCLUDE 'trcomm.h'
C     
      IF(NGT.GE.NTM) RETURN
C
      IF(MDLST.EQ.0) RETURN
C      
      IF(TSST.GT.T) RETURN
C      
      NGST=NGST+1
      K=0
      L=IZERO
      IX=INT((NRMAX+1-IZERO)/(NGPST-1))
C
      GTS(NGST) = GCLIP(T)
C
  100 GVT(NGST,K+89) = GCLIP(RT(L,1))
      GVT(NGST,K+90) = GCLIP(RT(L,2))
      GVT(NGST,K+91) = GCLIP(RT(L,3))
      GVT(NGST,K+92) = GCLIP(RT(L,4))
C
      K=K+4
      L=L+IX
      IF(L.LE.IZERO+(NGPST-1)*IX) GOTO 100
C
      RETURN
      END
C
C     ***********************************************************
C
C           CHECK NEGATIVE DENSITY OR TEMPERATURE
C
C     ***********************************************************
C
      SUBROUTINE TRCHCK(ICHCK)
C
      INCLUDE 'trcomm.h'
C
      ICHCK = 0
      DO NEQ=1,NEQMAX
         NSSN=NSS(NEQ)
         IF(NSSN.NE.0) THEN
            DO NR=1,NRMAX
               IF(RT(NR,NSSN).LT.0.D0) GOTO 100
            ENDDO
         ENDIF
      ENDDO
      GOTO 9000
C
  100 CONTINUE
      IND=0
      WRITE(6,*)
     &   'XX ERROR : NEGATIVE TEMPERATURE AT STEP ',NT
      DO NEQ=1,NEQMAX
         NSSN=NSS(NEQ)
         IF(NSSN.EQ.1.AND.IND.NE.1) THEN
            WRITE(6,*) '     TE (',NR,')=',RT(NR,NSSN)
            IND=1
         ELSEIF(NSSN.EQ.2.AND.IND.NE.2) THEN
            WRITE(6,*) '     TD (',NR,')=',RT(NR,NSSN)
            IND=2
         ELSEIF(NSSN.EQ.3.AND.IND.NE.3) THEN
            WRITE(6,*) '     TT (',NR,')=',RT(NR,NSSN)
            IND=3
         ELSEIF(NSSN.EQ.4.AND.IND.NE.4) THEN
            WRITE(6,*) '     TA (',NR,')=',RT(NR,NSSN)
            IND=4
         ELSEIF(NSSN.EQ.5.AND.IND.NE.5) THEN
            WRITE(6,*) '     TI1(',NR,')=',RT(NR,NSSN)
            IND=5
         ELSEIF(NSSN.EQ.6.AND.IND.NE.6) THEN
            WRITE(6,*) '     TI2(',NR,')=',RT(NR,NSSN)
            IND=6
C         ELSEIF(NSSN.EQ.7.AND.IND.NE.7) THEN
C            WRITE(6,*) '     TNC(',NR,')=',RT(NR,NSSN)
C            IND=7
C         ELSEIF(NSSN.EQ.8.AND.IND.NE.8) THEN
C            WRITE(6,*) '     TNH(',NR,')=',RT(NR,NSSN)
C            IND=8
         ENDIF
      ENDDO
      ICHCK=1
C
 9000 RETURN
      END
C
C     ***********************************************************
C
C           DECIDE COEEFICIENTS FOR EQUATIONS
C
C     ***********************************************************
C
      SUBROUTINE TR_COEF_DECIDE(NR,NSW,DV53)
C
      INCLUDE 'trcomm.h'
      COMMON /TRLCL1/ A(NVM,NVM,NRM),B(NVM,NVM,NRM),C(NVM,NVM,NRM)
      COMMON /TRLCL2/ D(NVM,NRM)
      COMMON /TRLCL3/ RD(NEQM,NRM)
      DIMENSION F2C(2),SIG(4),FCB(4)
C
      DV11=DVRHO(NR)
      DV23=DVRHO(NR)**(2.D0/3.D0)
      DV53=DVRHO(NR)**(5.D0/3.D0)
C
      IF (MODELG.EQ.0) THEN
         FVL = RKAPS/RKAP
C         FVL = 1.D0
      ELSEIF (MODELG.EQ.3) THEN
         FVL = 1.D0
      ENDIF
C
      DO NEQ=1,NEQMAX
C
         NSSN=NSS(NEQ)
         NSVN=NSV(NEQ)
         IF(NSSN.NE.0) RTM(NSSN)=RT(NR,NSSN)*RKEV/(PA(NSSN)*AMM)
C
C     /* Coefficients of left hand side variables */
C
         IF(NSVN.EQ.0) THEN
            RD(NEQ,NR)=1.D0
         ELSEIF(NSVN.EQ.1) THEN
            RD(NEQ,NR)=DV11
         ELSEIF(NSVN.EQ.2) THEN
            RD(NEQ,NR)=1.5D0*DV53
         ELSEIF(NSVN.EQ.3) THEN
            RD(NEQ,NR)=DV11
         ENDIF
C     
C     /* Coefficients of main variables */
C     
         DO NF=1,4
            NRF=NR+(NF-2)
            IF(NRF.EQ.0.OR.NRF.GT.NRMAX) THEN
               FCB(NF)=0.D0
            ELSE
               FCB(NF)=DVRHO(NRF)*ABRHO(NRF)/TTRHO(NRF)
            ENDIF
         ENDDO
         IF(NR.EQ.NRMAX-1) FCB(4)=2.D0*FCB(3)-FCB(2)
C     
         DO NMK=1,3
            IF (NSW.EQ.2.OR.NSW.NE.NMK) THEN
               FA(NMK,NSW)=DVRHO(NR+(NMK-2))*AR1RHO(NR+(NMK-2))*FVL/DR
               FB(NMK,NSW)=DVRHO(NR+(NMK-2))*AR2RHO(NR+(NMK-2))*FVL
     &                    /(DR*DR)
            ELSE
               FA(NMK,NSW)=0.D0
               FB(NMK,NSW)=0.D0
            ENDIF
            IF(NSW.NE.3) THEN
               FC(NMK,NSW)=0.5D0*(FCB(NMK)+FCB(NMK+1))*FVL/(DR*DR)
            ELSE
               FC(NMK,NSW)=0.D0
            ENDIF
         ENDDO
C     
         SIG(1)= 2.D0
         SIG(2)= 2.D0
         SIG(3)=-2.D0
         SIG(4)=-2.D0
         DO NMK=1,4
            IF(NMK.EQ.1) THEN
               NI=1
               NJ=1
            ELSEIF(NMK.EQ.2) THEN
               NI=2
               NJ=2
            ELSEIF(NMK.EQ.3) THEN
               NI=2
               NJ=1
            ELSE
               NI=3
               NJ=2
            ENDIF
            NRJ=NR+(NJ-2)
C
            IF(MDNCLS.EQ.0) THEN
C
            IF(NSVN.EQ.0) THEN
               IF(NSW.NE.3) THEN
                  F2C(NJ)=ETA(NRJ+1)*TTRHO(NRJ+1)
     &                   /(AMYU0*ARRHO(NRJ+1)*DVRHO(NRJ+1))
                  VV(NEQ,NEQ  ,NMK,NSW)= SIG(NMK)*F2C(NJ)*FC(NI,NSW)
                  DD(NEQ,NEQ  ,NMK,NSW)=          F2C(NJ)*FC(NI,NSW)
               ENDIF
            ELSEIF(NSVN.EQ.1) THEN
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*AV(NRJ,NSSN)
                  DD(NEQ,NEQ  ,NMK,NSW)= FB(NI,NSW)*AD(NRJ,NSSN)
               ENDIF
            ELSEIF(NSVN.EQ.2) THEN
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*(AVK(NRJ,NSSN)
     &                                  +AV(NRJ,NSSN)*1.5D0)*DV23
                  DD(NEQ,NEQ  ,NMK,NSW)= FB(NI,NSW)*AK(NRJ,NSSN)*DV23
                  VV(NEQ,NEQ-1,NMK,NSW)= 0.D0*DV23
                  IF(NRJ.EQ.NRMAX) THEN
                  DD(NEQ,NEQ-1,NMK,NSW)= FB(NI,NSW)*(AD(NRJ,NSSN)
     &                                  *1.5D0-AK(NRJ,NSSN))*DV23
     &                                  *PTS(NSSN)
                  ELSE
                  DD(NEQ,NEQ-1,NMK,NSW)= FB(NI,NSW)*(AD(NRJ,NSSN)
     &                                  *1.5D0-AK(NRJ,NSSN))*DV23*0.5D0
     &                                  *(RT(NRJ,NSSN)+RT(NRJ+1,NSSN))
                  ENDIF
               ENDIF
            ELSEIF(NSVN.EQ.3) THEN
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*VOID
                  DD(NEQ,NEQ  ,NMK,NSW)= FB(NI,NSW)*VOID
               ENDIF
            ENDIF
C
            ELSE
C
C     *** NCLASS ***
C     *
C     This expression below is the interim way to solve the problem.
C     That is because this simple expression cannot be used if we
C       consider impurity transports.
            IF(MDLEQ0.EQ.1) THEN
               NEQLMAX=NEQMAX-2 ! This means "Local NEQMAX".
            ELSE
               NEQLMAX=NEQMAX
            ENDIF
C
            IF(NSVN.EQ.0) THEN
               IF(NSW.NE.3) THEN
                  F2C(NJ)=ETA(NRJ+1)*TTRHO(NRJ+1)
     &                   /(AMYU0*ARRHO(NRJ+1)*DVRHO(NRJ+1))
                  VV(NEQ,NEQ  ,NMK,NSW)= SIG(NMK)*F2C(NJ)*FC(NI,NSW)
                  DD(NEQ,NEQ  ,NMK,NSW)=          F2C(NJ)*FC(NI,NSW)
               ENDIF
            ELSEIF(NSVN.EQ.1) THEN
               IF(NSSN.EQ.7.OR.NSSN.EQ.8) THEN
                  IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                     VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*AV(NRJ,NSSN)
                     DD(NEQ,NEQ  ,NMK,NSW)= FB(NI,NSW)*AD(NRJ,NSSN)
                  ENDIF
               ELSE
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  IF(NRJ.EQ.NRMAX) THEN
                     VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)
     &                      *(AV(NRJ,NSSN)+RGFLS(NRJ,5,NSSN)/PNSS(NSSN))
                  DO NEQ1=1,NEQLMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)
     &                       * PNSS(NSSN)
     &                       * ADNCLA(NRJ,NSSN1,NSSN,NSVN)
     &                       / PNSS(NSSN1)
                     ELSEIF(NSVN1.EQ.2) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)
     &                       * PNSS(NSSN)
     &                       * ADNCLB(NRJ,NSSN1,NSSN,NSVN)
     &                       /(PNSS(NSSN1)*PTS(NSSN1))
                     ENDIF
                  ENDDO
                  ELSE
                     VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)
     &                       *(AV(NRJ,NSSN)+RGFLS(NRJ,5,NSSN)
     &                        /(0.5D0*(RN(NRJ,NSSN)+RN(NRJ+1,NSSN))))
                  DO NEQ1=1,NEQLMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)
     &                       * 0.5D0*(RN(NRJ,NSSN )+RN(NRJ+1,NSSN ))
     &                       * ADNCLA(NRJ,NSSN1,NSSN,NSVN)
     &                       /(0.5D0*(RN(NRJ,NSSN1)+RN(NRJ+1,NSSN1)))
                     ELSEIF(NSVN1.EQ.2) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)
     &                       * 0.5D0*(RN(NRJ  ,NSSN )+RN(NRJ+1,NSSN ))
     &                       * ADNCLB(NRJ,NSSN1,NSSN,NSVN)
     &                       /(0.5D0*(RN(NRJ  ,NSSN1)*RT(NRJ  ,NSSN1)
     &                               +RN(NRJ+1,NSSN1)*RT(NRJ+1,NSSN1)))
                     ENDIF
                  ENDDO
                  ENDIF
               ENDIF
               ENDIF
            ELSEIF(NSVN.EQ.2) THEN
               CC=1.5D0
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  IF(NRJ.EQ.NRMAX) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*DV23
     &                    *( AVK(NRJ,NSSN)+AV(NRJ,NSSN)*CC
     &                      +RQFLS(NRJ,5,NSSN)
     &                      /(PNSS(NSSN)*PTS(NSSN))
     &                      +RGFLS(NRJ,5,NSSN)*CC
     &                      / PNSS(NSSN))
                  DO NEQ1=1,NEQLMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)*DV23
     &                       * PNSS(NSSN)*PTS(NSSN)
     &                       *( AKNCLB(NRJ,NSSN1,NSSN)
     &                         +ADNCLA(NRJ,NSSN1,NSSN,NSVN)*CC)
     &                       /PNSS(NSSN1)
                     ELSEIF(NSVN1.EQ.2) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)*DV23
     &                       * PNSS(NSSN)*PTS(NSSN)
     &                       *( AKNCLA(NRJ,NSSN1,NSSN)
     &                         +ADNCLB(NRJ,NSSN1,NSSN,NSVN)*CC)
     &                       /(PNSS(NSSN1)*PTS(NSSN1))
                     ENDIF
                  ENDDO
                  ELSE
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*DV23
     &                    *( AVK(NRJ,NSSN)+AV(NRJ,NSSN)*CC
     &                      +RQFLS(NRJ,5,NSSN)
     &                      /(0.5D0*( RN(NRJ  ,NSSN)*RT(NRJ  ,NSSN)
     &                               +RN(NRJ+1,NSSN)*RT(NRJ+1,NSSN)))
     &                      +RGFLS(NRJ,5,NSSN)*CC
     &                      /(0.5D0*( RN(NRJ,NSSN)  +RN(NRJ+1,NSSN))))
                  DO NEQ1=1,NEQLMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)*DV23
     &                       * 0.5D0*( RN(NRJ  ,NSSN )*RT(NRJ  ,NSSN )
     &                                +RN(NRJ+1,NSSN )*RT(NRJ+1,NSSN ))
     &                       *( AKNCLB(NRJ,NSSN1,NSSN)
     &                         +ADNCLA(NRJ,NSSN1,NSSN,NSVN)*CC)
     &                       /(0.5D0*( RN(NRJ  ,NSSN1)+RN(NRJ+1,NSSN1)))
                     ELSEIF(NSVN1.EQ.2) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)*DV23
     &                       * 0.5D0*( RN(NRJ  ,NSSN )*RT(NRJ  ,NSSN )
     &                                +RN(NRJ+1,NSSN )*RT(NRJ+1,NSSN ))
     &                       *( AKNCLA(NRJ,NSSN1,NSSN)
     &                         +ADNCLB(NRJ,NSSN1,NSSN,NSVN)*CC)
     &                       /(0.5D0*( RN(NRJ  ,NSSN1)*RT(NRJ  ,NSSN1)
     &                                +RN(NRJ+1,NSSN1)*RT(NRJ+1,NSSN1)))
                     ENDIF
                  ENDDO
                  ENDIF
               ENDIF
            ELSEIF(NSVN.EQ.3) THEN
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*VOID
                  DD(NEQ,NEQ  ,NMK,NSW)= FB(NI,NSW)*VOID
               ENDIF
            ENDIF
C     *
C     **************
C
            ENDIF
         ENDDO
C
      ENDDO
C
      IF(MDLEQE.EQ.1) THEN
         IND=0
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            IF(NSSN.EQ.1.AND.NSVN.EQ.1) THEN
               DO NMK=1,4
                  VV(NEQ,NEQ,NMK,NSW)=0.D0
                  DD(NEQ,NEQ,NMK,NSW)=0.D0
               ENDDO
               DO NEQ1=1,NEQMAX
                  NSSN1=NSS(NEQ1)
                  NSVN1=NSV(NEQ1)
                  IF(NSSN1.NE.1.AND.NSVN1.EQ.1) THEN
                     DO NMK=1,4
                        VV(NEQ,NEQ1,NMK,NSW) = VV(NEQ1,NEQ1,NMK,NSW)
     &                                        *PZ(NSSN1)
                        DD(NEQ,NEQ1,NMK,NSW) = DD(NEQ1,NEQ1,NMK,NSW)
     &                                        *PZ(NSSN1)
                     ENDDO
                  ENDIF
               ENDDO
               IND=IND+1
            ELSEIF(NSSN.EQ.1.AND.NSVN.EQ.2) THEN
               DO NMK=1,4
                  VV(NEQ,NEQ-1,NMK,NSW)=0.D0
                  DD(NEQ,NEQ-1,NMK,NSW)=0.D0
               ENDDO
               DO NEQ1=1,NEQMAX
                  NSSN1=NSS(NEQ1)
                  NSVN1=NSV(NEQ1)
                  IF(NSSN1.NE.1.AND.NSVN1.EQ.2) THEN
                     DO NMK=1,4
                        VV(NEQ,NEQ1-1,NMK,NSW) = VV(NEQ1,NEQ1-1,NMK,NSW)
     &                                          *PZ(NSSN1)
                        DD(NEQ,NEQ1-1,NMK,NSW) = DD(NEQ1,NEQ1-1,NMK,NSW)
     &                                          *PZ(NSSN1)
                     ENDDO
                  ENDIF
               ENDDO
               IND=IND+1
            ENDIF
            IF(IND.EQ.2) GOTO 1000
         ENDDO
      ENDIF
 1000 CONTINUE
C
C     /* Interim Parameter */
C
      IF(NR.EQ.NRMAX) THEN
      DO NV=1,NEQMAX
      DO NW=1,NEQMAX
         NO=1
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO,NSW)+VV(NV,NW,NO+2,NSW))
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO,NSW)+DD(NV,NW,NO+2,NSW))
         NO=2
         VI(NV,NW,NO,NSW)=VV(NV,NW,NO,NSW)
         DI(NV,NW,NO,NSW)=DD(NV,NW,NO,NSW)
      ENDDO
      ENDDO
      ELSE
      DO NV=1,NEQMAX
      DO NW=1,NEQMAX
      DO NO=1,2
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO,NSW)+VV(NV,NW,NO+2,NSW))
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO,NSW)+DD(NV,NW,NO+2,NSW))
      ENDDO
      ENDDO
      ENDDO
      ENDIF
C     
C     /* Coefficients of source term */
C     
      DO NEQ=1,NEQMAX
         NSSN=NSS(NEQ)
         NSVN=NSV(NEQ)
         IF(NR.EQ.NRMAX.AND.MDNCLS.EQ.0) THEN
            IF(NSSN.EQ.0) THEN
               D(NEQ,NR) = 0.D0
            ELSEIF(NSSN.EQ.5) THEN
               D(NEQ,NR) = VOID
            ELSEIF(NSSN.EQ.6) THEN
               D(NEQ,NR) = VOID
            ELSEIF(NSSN.EQ.7.OR.NSSN.EQ.8) THEN
               D(NEQ,NR) = SSIN(NR,NSSN)*DV11
     &                 +(-VI(NEQ,NEQ,2,NSW)+DI(NEQ,NEQ,2,NSW)*2.D0)
     &                 *PNSS(NSSN)
            ELSE
               IF(NSVN.EQ.1) THEN
                  D(NEQ,NR) = (SSIN(NR,NSSN)+SPE(NR,NSSN)/DT)*DV11
     &                 +(-VI(NEQ,NEQ,2,NSW)+DI(NEQ,NEQ,2,NSW)*2.D0)
     &                 *PNSS(NSSN)
               ELSEIF(NSVN.EQ.2) THEN
                  D(NEQ,NR) = (PIN(NR,NSSN)/(RKEV*1.D20)  )*DV53
     &                 +(-VI(NEQ,NEQ-1,2,NSW)+DI(NEQ,NEQ-1,2,NSW)*2.D0)
     &                 *PNSS(NSSN)
     &                 +(-VI(NEQ,NEQ  ,2,NSW)+DI(NEQ,NEQ  ,2,NSW)*2.D0)
     &                 *PNSS(NSSN)*PTS(NSSN)
               ELSEIF(NSVN.EQ.3) THEN
                  D(NEQ,NR) = VOID
               ELSE
                  STOP 'ERROR: must be NSSV=0 if NSSN=0'
               ENDIF
            ENDIF
C     *** NCLASS ***
C     *
         ELSEIF(NR.EQ.NRMAX.AND.MDNCLS.NE.0) THEN
            IF(NSSN.EQ.0) THEN
               D(NEQ,NR) = 0.D0
            ELSEIF(NSSN.EQ.5) THEN
               D(NEQ,NR) = VOID
            ELSEIF(NSSN.EQ.6) THEN
               D(NEQ,NR) = VOID
            ELSEIF(NSSN.EQ.7.OR.NSSN.EQ.8) THEN
               D(NEQ,NR) = SSIN(NR,NSSN)*DV11
               DO NEQ1=1,NEQMAX
                  NSSN1=NSS(NEQ1)
                  NSVN1=NSV(NEQ1)
                  IF(NSVN1.EQ.1) THEN
                     D(NEQ,NR) = D(NEQ,NR)
     &                    +(-VI(NEQ,NEQ1,2,NSW)+DI(NEQ,NEQ1,2,NSW)*2.D0)
     &                    *PNSS(NSSN1)
                  ELSEIF(NSVN1.EQ.2) THEN
                     D(NEQ,NR) = D(NEQ,NR)
     &                    +(-VI(NEQ,NEQ1,2,NSW)+DI(NEQ,NEQ1,2,NSW)*2.D0)
     &                    *PNSS(NSSN1)*PTS(NSSN1)
                  ENDIF
               ENDDO
            ELSE
               IF(NSVN.EQ.1) THEN
                  D(NEQ,NR) = (SSIN(NR,NSSN)+SPE(NR,NSSN)/DT)*DV11
                  DO NEQ1=1,NEQMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1) THEN
                        D(NEQ,NR) = D(NEQ,NR)
     &                 +(-VI(NEQ,NEQ1,2,NSW)+DI(NEQ,NEQ1,2,NSW)*2.D0)
     &                 *PNSS(NSSN1)
                     ELSEIF(NSVN1.EQ.2) THEN
                        D(NEQ,NR) = D(NEQ,NR)
     &                 +(-VI(NEQ,NEQ1,2,NSW)+DI(NEQ,NEQ1,2,NSW)*2.D0)
     &                 *PNSS(NSSN1)*PTS(NSSN1)
                     ENDIF
                  ENDDO
               ELSEIF(NSVN.EQ.2) THEN
                  D(NEQ,NR) = (PIN(NR,NSSN)/(RKEV*1.D20)  )*DV53
                  DO NEQ1=1,NEQMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1) THEN
                        D(NEQ,NR) = D(NEQ,NR)
     &                 +(-VI(NEQ,NEQ1,2,NSW)+DI(NEQ,NEQ1,2,NSW)*2.D0)
     &                 *PNSS(NSSN1)
                     ELSEIF(NSVN1.EQ.2) THEN
                        D(NEQ,NR) = D(NEQ,NR)
     &                 +(-VI(NEQ,NEQ1,2,NSW)+DI(NEQ,NEQ1,2,NSW)*2.D0)
     &                 *PNSS(NSSN1)*PTS(NSSN1)
                     ENDIF
                  ENDDO
               ELSEIF(NSVN.EQ.3) THEN
                  D(NEQ,NR) = VOID
               ELSE
                  STOP 'ERROR: must be NSSV=0 if NSSN=0'
               ENDIF
            ENDIF
C     *
C     *************
C
         ELSE
            IF(NSSN.EQ.0) THEN
               D(NEQ,NR) = ETA(NR  )*AR1RHO(NR  )*BB/(TTRHO(NR  )*RR
     &                    *ARRHO(NR  )*DR)*(AJ(NR  )-AJOH(NR  ))
     &                    -ETA(NR+1)*AR1RHO(NR+1)*BB/(TTRHO(NR+1)*RR
     &                    *ARRHO(NR+1)*DR)*(AJ(NR+1)-AJOH(NR+1))
            ELSEIF(NSSN.EQ.5) THEN
               D(NEQ,NR) = VOID
            ELSEIF(NSSN.EQ.6) THEN
               D(NEQ,NR) = VOID
            ELSEIF(NSSN.EQ.7.OR.NSSN.EQ.8) THEN
               D(NEQ,NR) = SSIN(NR,NSSN)*DV11
            ELSE
               IF(NSVN.EQ.1) THEN
                  D(NEQ,NR) = (SSIN(NR,NSSN)+SPE(NR,NSSN)/DT)*DV11
               ELSEIF(NSVN.EQ.2) THEN
                  D(NEQ,NR) = (PIN(NR,NSSN)/(RKEV*1.D20)    )*DV53
               ELSEIF(NSVN.EQ.3) THEN
                  D(NEQ,NR) = VOID
               ELSE
                  STOP 'ERROR: must be NSSV=0 if NSSN=0'
               ENDIF
            ENDIF
         ENDIF
      ENDDO
C
      IF(MDLEQE.EQ.1) THEN
C
      IF(NR.EQ.NRMAX) THEN
         IND=0
         VISUMN=0.D0
         DISUMN=0.D0
         VISUMT1=0.D0
         DISUMT1=0.D0
         VISUMT2=0.D0
         DISUMT2=0.D0
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            IF(NSSN.EQ.1.AND.NSVN.EQ.1) THEN
               D(NEQ,NR) = 0.D0
               DO NEQ1=1,NEQMAX
                  NSSN1=NSS(NEQ1)
                  NSVN1=NSV(NEQ1)
                  IF(NSSN1.NE.1.AND.NSVN1.EQ.1) THEN
                     D(NEQ,NR) = D(NEQ,NR)+PZ(NSSN1)
     &                          *(SSIN(NR,NSSN1)+SPE(NR,NSSN1)/DT)*DV11
                     VISUMN = VISUMN+PZ(NSSN1)*VI(NEQ1,NEQ1,2,NSW)
     &                              *PNSS(NSSN1)
                     DISUMN = DISUMN+PZ(NSSN1)*DI(NEQ1,NEQ1,2,NSW)
     &                              *PNSS(NSSN1)
                  ENDIF
               ENDDO
               D(NEQ,NR) = D(NEQ,NR)+(-VISUMN+DISUMN*2.D0)
               IND=IND+1
            ELSEIF(NSSN.EQ.1.AND.NSVN.EQ.2) THEN
               D(NEQ,NR) = 0.D0
               DO NEQ1=1,NEQMAX
                  NSSN1=NSS(NEQ1)
                  NSVN1=NSV(NEQ1)
                  IF(NSSN1.NE.1.AND.NSVN1.EQ.2) THEN
                     D(NEQ,NR) = D(NEQ,NR)+PZ(NSSN1)
     &                          *(PIN(NR,NSSN1)/(RKEV*1.D20)  )*DV53
                     VISUMT1 = VISUMT1+PZ(NSSN1)*VI(NEQ1,NEQ1-1,2,NSW)
     &                                *PNSS(NSSN1)
                     DISUMT1 = DISUMT1+PZ(NSSN1)*DI(NEQ1,NEQ1-1,2,NSW)
     &                                *PNSS(NSSN1)
                     VISUMT2 = VISUMT2+PZ(NSSN1)*VI(NEQ1,NEQ1  ,2,NSW)
     &                                *PNSS(NSSN1)*PTS(NSSN1) 
                     DISUMT2 = DISUMT2+PZ(NSSN1)*DI(NEQ1,NEQ1  ,2,NSW)
     &                                *PNSS(NSSN1)*PTS(NSSN1) 
                  ENDIF
               ENDDO
               D(NEQ,NR) = D(NEQ,NR)+(-VISUMT1+DISUMT1*2.D0)
     &                              +(-VISUMT2+DISUMT2*2.D0)
               IND=IND+1
            ENDIF
           IF(IND.EQ.2) GOTO 2000
         ENDDO
      ELSE
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            IF(NSSN.EQ.1.AND.NSVN.EQ.1) THEN
               D(NEQ,NR) = 0.D0
               DO NEQ1=1,NEQMAX
                  NSSN1=NSS(NEQ1)
                  NSVN1=NSV(NEQ1)
                  IF(NSSN1.NE.1.AND.NSVN1.EQ.1) THEN
                     D(NEQ,NR)=D(NEQ,NR)+PZ(NSSN1)*D(NEQ1,NR)
                  ENDIF
               ENDDO
               GOTO 2000
            ENDIF
         ENDDO
      ENDIF
C
      ENDIF
 2000 CONTINUE
C
      RETURN
      END
C
C     ***********************************************************
C
C           IONIZATION 
C
C     ***********************************************************
C
      SUBROUTINE TR_IONIZATION(NR)
C
      INCLUDE 'trcomm.h'
      COMMON /TRLCL1/ A(NVM,NVM,NRM),B(NVM,NVM,NRM),C(NVM,NVM,NRM)
C
      IF(MDLEQ0.EQ.1) THEN
         DO NV=1,NEQMAX
            NSSV=NSS(NV)
            NSVV=NSV(NV)
            IF(NSSV.EQ.1.AND.NSVV.EQ.1) THEN
               DO NW=1,NEQMAX
                  NSSW=NSS(NW)
                  IF(NSSW.EQ.7) THEN
                     B(NV,NW,NR)=B(NV,NW,NR)+TSIE(NR)*DVRHO(NR)
                  ENDIF
               ENDDO
            ELSEIF(NSSV.EQ.2.AND.NSVV.EQ.1) THEN
               DO NW=1,NEQMAX
                  NSSW=NSS(NW)
                  IF(NSSW.EQ.7) THEN
                     B(NV,NW,NR)=B(NV,NW,NR)+(PN(2)/(PN(2)+PN(3)))
     &                                      *TSIE(NR)*DVRHO(NR)
                  ENDIF
               ENDDO
            ELSEIF(NSSV.EQ.3.AND.NSVV.EQ.1) THEN
               DO NW=1,NEQMAX
                  NSSW=NSS(NW)
                  IF(NSSW.EQ.7) THEN
                     B(NV,NW,NR)=B(NV,NW,NR)+(PN(3)/(PN(2)+PN(3)))
     &                                      *TSIE(NR)*DVRHO(NR)
                  ENDIF
               ENDDO
            ELSE
               DO NW=1,NEQMAX
                  NSSW=NSS(NW)
                  IF(NV.EQ.NW.AND.NSSW.EQ.7) THEN
                     B(NV,NW,NR)=B(NV,NW,NR)-TSIE(NR)*DVRHO(NR)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
C           CHARGE EXCHANGE 
C
C     ***********************************************************
C
      SUBROUTINE TR_CHARGE_EXCHANGE(NR)
C
      INCLUDE 'trcomm.h'
      COMMON /TRLCL1/ A(NVM,NVM,NRM),B(NVM,NVM,NRM),C(NVM,NVM,NRM)
C
      IF(MDLEQ0.EQ.1) THEN
         DO NV=1,NEQMAX
            NSSV=NSS(NV)
            DO NW=1,NEQMAX
               NSSW=NSS(NW)
               IF(NV.EQ.NW.AND.NSSW.EQ.7) THEN
                  B(NV,NW,NR)=B(NV,NW,NR)-TSCX(NR)*DVRHO(NR)
               ELSEIF(NSSV.EQ.8.AND.NSSW.EQ.7) THEN
                  B(NV,NW,NR)=B(NV,NW,NR)+TSCX(NR)*DVRHO(NR)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
