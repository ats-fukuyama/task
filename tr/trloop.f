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
      PARAMETER(NURM=NRM+1,NUTM=61)
      COMMON /TRBCH1/ RAD(NURM),FRFHE(NURM,NUTM),FRFHI(NURM,NUTM)
      COMMON /TRBCH3/ FUT(NUTM),NUFMAX,NTXMAX
      DIMENSION DERIVR(NURM),URFHEPH(4,NURM),URFHIPH(4,NURM)
      DIMENSION DERIVT(NUTM),UFUTE(4,NUTM),UFUTI(4,NUTM)
      DIMENSION FRFHET(NUTM),FRFHIT(NUTM),FRFHER(NURM),FRFHIR(NURM)
C
      DIMENSION XX(MLM),YY(2,NRM)
C
      NT=0
      RIP=RIPS
      IF(NTMAX.NE.0) DIP=(RIPE-RIPS)/NTMAX
      ICHCK=0
C
      NQM=NVM
      MWMAX=MWM
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
      CALL TRCALC
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
      CALL TRATOX(NQM)
C
      DO J=1,NQM
      DO NR=1,NRMAX
         XX(NQM*(NR-1)+J)=XV(J,NR)
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
      CALL TRMTRX(NQM)
C
      CALL BANDRD (AX,X,NQM*NRMAX,MWMAX,MWM,IERR)
      IF(IERR.EQ.30000) THEN
         WRITE(6,*) 'XX ERROR IN TRLOOP : MATRIX AA IS SINGULAR ',
     &              ' AT ',NT,' STEP.'
         GOTO 9000
      ENDIF
      DO J=1,NFM
      DO NR=1,NRMAX
         Y(J,NR) = Y(J,NR)/AY(J,NR)
      ENDDO
      ENDDO
C
      DO I=1,NQM*NRMAX
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
      DO I=1,NQM*NRMAX
         XX(I) = X(I)
      ENDDO
      DO J=1,NFM
      DO NR=1,NRMAX
         YY(J,NR) = Y(J,NR)
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         RN(NR,1) = 0.D0
         DO NS=2,NSM
            RN(NR,NS) = 0.5D0*(XV(NS,NR)+X(NQM*(NR-1)+NS))
            RN(NR,1)  = RN(NR,1)+PZ(NS)*RN(NR,NS)
         ENDDO
         RN(NR,1) = RN(NR,1)+PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
         DO NS=1,NSM
            IF(NS.NE.NSM) THEN
               RT(NR,NS) = 0.5D0*(XV(NSM+NS,NR)+X(NQM*(NR-1)+(NSM+NS)))
     &                          /RN(NR,NS)
            ELSE
            IF(RN(NR,NSM).LT.1.D-70) THEN
               RT(NR,NS) = 0.D0
            ELSE
               RT(NR,NS) = 0.5D0*(XV(NSM+NS,NR)+X(NQM*(NR-1)+(NSM+NS)))
     &                          /RN(NR,NS)
            ENDIF
            ENDIF
         ENDDO
         BP(NR)    = 0.5D0*(XV(NQM,NR)+X(NQM*NR))
         RW(NR,1)  = 0.5D0*(YV(1,NR)+Y(1,NR))
         RW(NR,2)  = 0.5D0*(YV(2,NR)+Y(2,NR))
      ENDDO
C
      CALL TRCHCK(ICHCK)
      IF(ICHCK.EQ.1) GOTO 4000
C
      CALL TRCALC
      GOTO 2000
C
 4000 NT=NT+1
      T=T+DT
      VSEC=VSEC+VLOOP*DT
      IF(Q0.LT.1.D0) TST=TST+DT
      RIP=RIP+DIP
C      write(6,'(A,1P3E12.5)') "RIP,RIPE,DIP=",RIP,RIPE,DIP
C
      DO J=1,NQM
      DO NR=1,NRMAX
         XV(J,NR) = X(NQM*(NR-1)+J)
      ENDDO
      ENDDO
      DO NF=1,NFM
      DO NR=1,NRMAX
         YV(NF,NR) = Y(NF,NR)
      ENDDO
      ENDDO
C
      CALL TRXTOA(NQM)
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
      IF(TST+0.5D0*DT.GT.TPRST) THEN
         CALL TRSAWT
         TST=0.D0
      ENDIF
C
      CALL TRCALC
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
         CALL TRCONV(L,IERR)
         write(6,*) "L=",L
         IF(IERR.NE.0) RETURN
      ENDIF
      IF(NT.LT.NTMAX) GOTO 1000
C
 9000 RIPS=RIPE
      RETURN
      END
C
C     ***********************************************************
C
C           COMPUTE MATRIX ELEMENTS
C
C     ***********************************************************
C
      SUBROUTINE TRMTRX(NQM)
C
      INCLUDE 'trcomm.h'
C
      COMMON /TRLCL1/ A(NVM,NVM,NRM),B(NVM,NVM,NRM),C(NVM,NVM,NRM)
      COMMON /TRLCL2/ D(NVM,NRM)
C
      RKAPX=(RKAP-1.D0)/(RKAP+1.D0)
      FKAP=0.5D0*(RKAP+1.D0)
     &          *(1.D0+RKAPX/4.D0+RKAPX*RKAPX/64.D0)
      DO NS=1,NSM
         AMZ(NS)=PA(NS)*AMM/PZ(NS)**2
      ENDDO
C
      IF(MDLCD.EQ.0) THEN
C         BPS= AMYU0*RIP*1.D6/(2.D0*PI*RA*FKAP)
         BPS= AMYU0*RIP*1.D6/(2.D0*PI*RA*RKAP)
      ELSE
         BPA=XV(NQM,NRMAX)
         RLP=RA*(LOG(8.D0*RR/RA)-2.D0)
         BPS=BPA-AMYU0*DT*EZOH(NRMAX)/RLP
      ENDIF
C
      COEF = AEE**4*5.D0*1.D20/(SQRT(2.D0*PI)*PI*AEPS0**2)
C
      DO NR=1,NRMAX
      DO NW=1,NVM
      DO NV=1,NVM
         A(NV,NW,NR)=0.D0
         B(NV,NW,NR)=0.D0
         C(NV,NW,NR)=0.D0
      ENDDO
      ENDDO
      ENDDO
      DO NR=1,NRMAX
      DO NV=1,NVM
         D(NV,NR)=0.D0
      ENDDO
      ENDDO
      DO NW=1,NVM
      DO NV=1,NVM
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
C
C      NN=2*NSM+1
      NN=NQM
C
      DO NR=2,NRMAX-1
C
      NSW=2
      CALL TRIPTC(DV11,DV53,NSW,NR)
C
      DO NV=1,NQM
      DO NW=1,NQM
      DO NO=1,2
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
      ENDDO
         A(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)+DI(NV,NW,1,NSW)
         B(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)-DI(NV,NW,1,NSW)
     &                -0.5D0*VI(NV,NW,2,NSW)-DI(NV,NW,2,NSW)
         C(NV,NW,NR) =-0.5D0*VI(NV,NW,2,NSW)+DI(NV,NW,2,NSW)
      ENDDO
      ENDDO
C
      DO NS=1,NSM
         D(    NS,NR) = (SSIN(NR,NS)+SPE(NR,NS)/DT)*DV11
         D(NSM+NS,NR) = (PIN(NR,NS)/(RKEV*1.D20)  )*DV53
         D(   NQM,NR) = ETA(NR  )*AR1RHO(NR  )*BB/(TTRHO(NR  )*RR
     &                 *ARRHO(NR  )*DR)*(AJ(NR  )-AJOH(NR  ))
     &                 -ETA(NR+1)*AR1RHO(NR+1)*BB/(TTRHO(NR+1)*RR
     &                 *ARRHO(NR+1)*DR)*(AJ(NR+1)-AJOH(NR+1))
C
      DO NS1=1,NSM
         IF(NS1.NE.NS) THEN
            C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1))*DV53
     &             *1.5D0
            B(NSM+NS,NSM+NS, NR)=B(NSM+NS,NSM+NS, NR)-RN(NR,NS1)*C1
            B(NSM+NS,NSM+NS1,NR)=B(NSM+NS,NSM+NS1,NR)+RN(NR,NS )*C1
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      NR=1
C     
      NSW=1
      CALL TRIPTC(DV11,DV53,NSW,NR)
C
C      NO=2
      DO NV=1,NQM
      DO NW=1,NQM
      DO NO=1,2
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
      ENDDO
         A(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)+DI(NV,NW,1,NSW)
         B(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)-DI(NV,NW,1,NSW)
     &                -0.5D0*VI(NV,NW,2,NSW)-DI(NV,NW,2,NSW)
         C(NV,NW,NR) =-0.5D0*VI(NV,NW,2,NSW)+DI(NV,NW,2,NSW)
      ENDDO
      ENDDO
C
      DO NS=1,NSM
         D(    NS,NR) = (SSIN(NR,NS)+SPE(NR,NS)/DT)*DV11
         D(NSM+NS,NR) = (PIN(NR,NS)/(RKEV*1.D20)  )*DV53
         D(   NQM,NR) = ETA(NR  )*AR1RHO(NR  )*BB/(TTRHO(NR  )*RR
     &                 *ARRHO(NR  )*DR)*(AJ(NR  )-AJOH(NR  ))
     &                 -ETA(NR+1)*AR1RHO(NR+1)*BB/(TTRHO(NR+1)*RR
     &                 *ARRHO(NR+1)*DR)*(AJ(NR+1)-AJOH(NR+1))
C
      DO NS1=1,NSM
         IF(NS1.NE.NS) THEN
            C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1))*DV53
     &             *1.5D0
            B(NSM+NS,NSM+NS, NR)=B(NSM+NS,NSM+NS, NR)-RN(NR,NS1)*C1
            B(NSM+NS,NSM+NS1,NR)=B(NSM+NS,NSM+NS1,NR)+RN(NR,NS )*C1
         ENDIF
      ENDDO
      ENDDO
C
      NR=NRMAX
C
      NSW=3
      CALL TRIPTC(DV11,DV53,NSW,NR)
C
      DO NV=1,NQM
      DO NW=1,NQM
         NO=1
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
C
         NO=2
         VI(NV,NW,NO,NSW)=VV(NV,NW,NO  ,NSW)
         DI(NV,NW,NO,NSW)=DD(NV,NW,NO  ,NSW)
C
C     *** ATTENTION PLEASE!! **************************************
C     *                                                           *
C     * NO=2 section is slightly diffrent from that of tr.020319. *
C     *                                                           *
C     *************************************************************
C
         A(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)+DI(NV,NW,1,NSW)
         B(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)-DI(NV,NW,1,NSW)
     &                -0.0D0*VI(NV,NW,2,NSW)-DI(NV,NW,2,NSW)
     &                                      -DI(NV,NW,2,NSW)
         C(NV,NW,NR) =-0.0D0
      ENDDO
      ENDDO
C
      DO NS=1,NSM
         D(    NS,NR) = (SSIN(NR,NS)+SPE(NR,NS)/DT)*DV11
     &                 +(-VI(NS,NS,2,NSW)+DI(NS,NS,2,NSW)*2.D0)*PNSS(NS)
         D(NSM+NS,NR) = (PIN(NR,NS)/(RKEV*1.D20)  )*DV53
     &                 +(-VI(NSM+NS,NS,2,NSW)+DI(NSM+NS,NS,2,NSW)*2.D0)
     &                 *PNSS(NS)
     &                 +(-VI(NSM+NS,NSM+NS,2,NSW)
     &                 +  DI(NSM+NS,NSM+NS,2,NSW)*2.D0)
     &                 *PNSS(NS)*PTS(NS)
         D(   NQM,NR) = 0.D0
C
      DO NS1=1,NSM
         IF(NS1.NE.NS) THEN
            C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1))*DV53
     &             *1.5D0
            B(NSM+NS,NSM+NS, NR)=B(NSM+NS,NSM+NS, NR)-RN(NR,NS1)*C1
            B(NSM+NS,NSM+NS1,NR)=B(NSM+NS,NSM+NS1,NR)+RN(NR,NS )*C1
         ENDIF
      ENDDO
      ENDDO
C
C      FADV=0.5D0
      FADV=1.0D0
C
      PRV=(1.D0-FADV)*DT
      ADV=FADV*DT
C
      DO NR=1,NRMAX
         DV11=DVRHO(NR)
         DV53=DVRHO(NR)**(5.D0/3.D0)
         DO NS=1,NSM
            NV=NS
            X(NQM*(NR-1)+NV) =       DV11*XV(NV,NR)+DT*D(NV,NR)
            NV=NSM+NS
            X(NQM*(NR-1)+NV) = 1.5D0*DV53*XV(NV,NR)+DT*D(NV,NR)
         ENDDO
         NV=NQM
         X(NQM*(NR-1)+NV) = XV(NV,NR)+DT*D(NV,NR)
      ENDDO
C
      NR=1
      DO NW=1,NQM
      DO NV=1,NQM
         X(NQM*(NR-1)+NV)=X(NQM*(NR-1)+NV)+PRV*(
     &                                         +B(NV,NW,NR)*XV(NW,NR  )
     &                                         +C(NV,NW,NR)*XV(NW,NR+1))
      ENDDO
      ENDDO
C
      DO NR=2,NRMAX-1
      DO NW=1,NQM
      DO NV=1,NQM
         X(NQM*(NR-1)+NV)=X(NQM*(NR-1)+NV)+PRV*(A(NV,NW,NR)*XV(NW,NR-1)
     &                                         +B(NV,NW,NR)*XV(NW,NR  )
     &                                         +C(NV,NW,NR)*XV(NW,NR+1))
      ENDDO
      ENDDO
      ENDDO
C
      NR=NRMAX
      DO NW=1,NQM
      DO NV=1,NQM
         X(NQM*(NR-1)+NV)=X(NQM*(NR-1)+NV)+PRV*(A(NV,NW,NR)*XV(NW,NR-1)
     &                                         +B(NV,NW,NR)*XV(NW,NR  )
     &                                                                 )
      ENDDO
      ENDDO
C
      DO MV=1,NVM*NRMAX
      DO MW=1,MWM
         AX(MW,MV) = 0.D0
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX
      DO NV=1,NQM
      DO NW=1,NQM
         AX(  NQM+NW-NV,NQM*(NR-1)+NV) = -ADV*A(NV,NW,NR)
         AX(2*NQM+NW-NV,NQM*(NR-1)+NV) = -ADV*B(NV,NW,NR)
         AX(3*NQM+NW-NV,NQM*(NR-1)+NV) = -ADV*C(NV,NW,NR)
      ENDDO
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         DV11=DVRHO(NR)
         DV53=DVRHO(NR)**(5.D0/3.D0)
         DO NS=1,NSM
            AX(2*NQM,NQM*(NR-1)+NS) 
     &    = AX(2*NQM,NQM*(NR-1)+NS)     + DV11
            AX(2*NQM,NQM*(NR-1)+NSM+NS) 
     &    = AX(2*NQM,NQM*(NR-1)+NSM+NS) + 1.5D0*DV53
         ENDDO
            AX(2*NQM,NQM*(NR-1)+NN) 
     &    = AX(2*NQM,NQM*(NR-1)+NN)     + 1.D0
      ENDDO
C
C     ***** Surface Boundary Condition for Bp *****
C
      MV=NQM*NRMAX
      DO MW=1,MWM
         AX(MW,MV)=0.D0
      ENDDO
      AX(2*NQM,MV)=1.D0
      X(MV)=BPS
C
C     ***** Evolution of fast ion components *****
C
      DO NR=1,NRMAX
         Y(1,NR)=(1.5D0-PRV/TAUB(NR))*YV(1,NR)+PNB(NR)*DT/(RKEV*1.D20)
         Y(2,NR)=(1.5D0-PRV/TAUF(NR))*YV(2,NR)+PNF(NR)*DT/(RKEV*1.D20)
         AY(1,NR)=1.5D0+ADV/TAUB(NR)
         AY(2,NR)=1.5D0+ADV/TAUF(NR)
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
      SUBROUTINE TRATOX(NQM)
C
      INCLUDE 'trcomm.h'
C
      DO NR=1,NRMAX
      DO NS=1,NSM
         XV(    NS,NR) = RN(NR,NS)
         XV(NSM+NS,NR) = RN(NR,NS)*RT(NR,NS)
      ENDDO
         XV(NQM,NR) = BP(NR)
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
      SUBROUTINE TRXTOA(NQM)
C
      INCLUDE 'trcomm.h'
C
      DO NR=1,NRMAX
         RN(NR,1) = 0.D0
         DO NS=2,NSM
            RN(NR,NS) = XV(NS,NR)
            RN(NR,1)  = RN(NR,1)+PZ(NS)*RN(NR,NS)
         ENDDO
         RN(NR,1) = RN(NR,1)+PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
         DO NS=1,NSM
            IF(NS.NE.NSM) THEN
               RT(NR,NS) = XV(NSM+NS,NR)/RN(NR,NS)
            ELSE
               IF(RN(NR,NSM).LT.1.D-70) THEN
                  RT(NR,NS) = 0.D0
               ELSE
                  RT(NR,NS) = XV(NSM+NS,NR)/RN(NR,NS)
               ENDIF
            ENDIF
         ENDDO
         BP(NR)   = XV(NQM,NR)
         RW(NR,1) = YV(1,NR)
         RW(NR,2) = YV(2,NR)
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
      DO NS=1,NSM
      DO NR=1,NRMAX
         IF(RT(NR,NS).LT.0.D0) THEN
            GOTO 100
         ENDIF
      ENDDO
      ENDDO
      GOTO 9000
C
  100 CONTINUE
      WRITE(6,*)
     &   'XX ERROR : NEGATIVE TEMPERATURE AT STEP ',NT
      WRITE(6,*) '     TE(',NR,')=',RT(NR,1)
      WRITE(6,*) '     TD(',NR,')=',RT(NR,2)
      WRITE(6,*) '     TT(',NR,')=',RT(NR,3)
      WRITE(6,*) '     TA(',NR,')=',RT(NR,4)
      ICHCK=1
C
 9000 RETURN
      END
C
C     ***********************************************************
C
C           INPUT COEFFECIENTS OF TRANSPORT EQUATIONS
C
C     ***********************************************************
C
      SUBROUTINE TRIPTC(DV11,DV53,NSW,NR)
C
      INCLUDE 'trcomm.h'
C
      DIMENSION F2C(2),SIG(4),FCB(4)
C
      DO NS=1,NSM
         RTM(NS)=RT(NR,NS)*RKEV/(PA(NS)*AMM)
      ENDDO
C
      DV11=DVRHO(NR)
      DV23=DVRHO(NR)**(2.D0/3.D0)
      DV53=DVRHO(NR)**(5.D0/3.D0)
C
      IF (MODELG.EQ.0) THEN
CCC         FVL = FKAP/RKAP
         FVL = 1.D0
      ELSEIF (MODELG.EQ.3) THEN
         FVL = 1.D0
      ENDIF
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
            FB(NMK,NSW)=DVRHO(NR+(NMK-2))*AR2RHO(NR+(NMK-2))*FVL/(DR*DR)
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
      NSX=2*NSM+1
      DO NMK=1,4
         IF (NMK.EQ.1) THEN
            NI=1
            NJ=1
         ELSEIF (NMK.EQ.2) THEN
            NI=2
            NJ=2
         ELSEIF (NMK.EQ.3) THEN
            NI=2
            NJ=1
         ELSE
            NI=3
            NJ=2
         ENDIF
         NRJ=NR+(NJ-2)
         IF (NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
            DO NS=1,NSM
            VV(    NS,    NS,NMK,NSW) = FA(NI,NSW)*AV(NRJ,NS)
            DD(    NS,    NS,NMK,NSW) = FB(NI,NSW)*AD(NRJ,NS)
            VV(NSM+NS,NSM+NS,NMK,NSW) = FA(NI,NSW)*(AVK(NRJ,NS)
     &                                 +AV(NRJ,NS)*1.5D0)*DV23
            DD(NSM+NS,NSM+NS,NMK,NSW) = FB(NI,NSW)*AK(NRJ,NS)*DV23
            VV(NSM+NS,    NS,NMK,NSW) = 0.D0*DV23
            IF(NRJ.EQ.NRMAX) THEN
            DD(NSM+NS,    NS,NMK,NSW) = FB(NI,NSW)*(AD(NRJ,NS)
     &                                       *1.5D0-AK(NRJ,NS))*DV23
     &                                       *PTS(NS)
            ELSE
            DD(NSM+NS,    NS,NMK,NSW) = FB(NI,NSW)*(AD(NRJ,NS)
     &                                       *1.5D0-AK(NRJ,NS))*DV23
     &                                 *0.5D0*(RT(NRJ,NS)+RT(NRJ+1,NS))
            ENDIF
            ENDDO
         ENDIF
         IF(NSW.NE.3) THEN
            IF(NR.EQ.NRMAX) THEN
               F2C(NJ)=0.D0
            ELSE
               F2C(NJ)=ETA(NRJ+1)*TTRHO(NRJ+1)
     &               /(AMYU0*ARRHO(NRJ+1)*DVRHO(NRJ+1))
            ENDIF
            VV(   NSX,   NSX,NMK,NSW) = SIG(NMK)*F2C(NJ)*FC(NI,NSW)
            DD(   NSX,   NSX,NMK,NSW) =          F2C(NJ)*FC(NI,NSW)
         ENDIF
      ENDDO
C
      RETURN
      END
