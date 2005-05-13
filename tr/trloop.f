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
      INCLUDE 'trcomm.inc'
      COMMON /TMSLC1/ TMU(NTUM),TMU1(NTUM)
      COMMON /TMSLC2/ NTAMAX
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
C
      DIMENSION XX(MLM),YY(2,NRM),ZZ(NSM,NRM)
C
      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         IF(NTMAX.GT.NTAMAX) NTMAX=NTAMAX
      ELSE
         NT=0
         RIP=RIPS
         IF(NTMAX.NE.0) DIP=(RIPE-RIPS)/DBLE(NTMAX)
      ENDIF
      ICHCK=0
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
      IF(MDTC.NE.0) THEN
         DO J=1,NSMAX
         DO NR=1,NRMAX
            ZZ(J,NR)=ZV(J,NR)
         ENDDO
         ENDDO
      ENDIF
C
 2000 CONTINUE
      CALL TR_EDGE_SELECTOR(0)
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
      IF(MDTC.NE.0) THEN
         DO J=1,NSMAX
         DO NR=1,NRMAX
            Z(J,NR) = Z(J,NR)/AZ(J,NR)
         ENDDO
         ENDDO
      ENDIF
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
      IF(MDTC.NE.0) THEN
         DO J=1,NSMAX
         DO NR=1,NRMAX
            IF (ABS(Z(J,NR)-ZZ(J,NR)).GT.EPSLTR*ABS(Z(J,NR))) GOTO 3000
         ENDDO
         ENDDO
      ENDIF
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
      IF(MDTC.NE.0) THEN
         DO J=1,NSMAX
         DO NR=1,NRMAX
            ZZ(J,NR) = Z(J,NR)
         ENDDO
         ENDDO
      ENDIF
C
      DO NR=1,NRMAX
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            NSTN=NST(NEQ)
            IF(NSVN.EQ.0) THEN
               IF(NSTN.EQ.0) THEN
                  RDP(NR) = XV(NEQ,NR)
               ELSE
                  RDP(NR) = 0.5D0*(XV(NEQ,NR)+X(NEQRMAX*(NR-1)+NSTN))
               ENDIF
            ELSEIF(NSVN.EQ.1) THEN
               IF(MDLEQN.NE.0) THEN !!!
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
               ELSEIF(NSSN.EQ.1.AND.MDLEQE.EQ.1) THEN
                  RN(NR,NSSN) = 0.5D0*(XV(NEQ,NR)+X(NEQRMAX*(NR-1)+NEQ))
               ELSEIF(NSSN.EQ.2.AND.MDLEQE.EQ.2) THEN
                  RN(NR,NSSN) = 0.D0
                  DO NEQ1=1,NEQMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     NSTN1=NST(NEQ1)
                     IF(NSVN1.EQ.1.AND.NSSN1.NE.2) THEN
                        IF(NSTN1.EQ.0) THEN
                           RN(NR,NSSN) = RN(NR,NSSN)+ABS(PZ(NSSN1))
     &                                  *XV(NEQ1,NR)
                        ELSE
                           RN(NR,NSSN) = RN(NR,NSSN)+ABS(PZ(NSSN1))
     &                      *0.5D0*(XV(NEQ1,NR)+X(NEQRMAX*(NR-1)+NSTN1))
                        ENDIF
                     ENDIF
                  ENDDO
                  RN(NR,NSSN) = RN(NR,NSSN)
     &                         +PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
               ELSE
                  IF(NSTN.EQ.0) THEN
                     RN(NR,NSSN) = XV(NEQ,NR)
                  ELSE
                     RN(NR,NSSN) = 0.5D0*(XV(NEQ,NR)
     &                                  +  X(NEQRMAX*(NR-1)+NSTN))
                  ENDIF
               ENDIF
            ELSE !!!
               IF(NSTN.EQ.0) THEN
                  RN(NR,NSSN) = XV(NEQ,NR)
               ELSE
                  RN(NR,NSSN) = 0.5D0*(XV(NEQ,NR)
     &                              +  X(NEQRMAX*(NR-1)+NSTN))
               ENDIF
            ENDIF !!!
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
         BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         RW(NR,1)  = 0.5D0*(YV(1,NR)+Y(1,NR))
         RW(NR,2)  = 0.5D0*(YV(2,NR)+Y(2,NR))
      ENDDO
C      write(6,*) L,XV(2,1)/RN(1,1),X(1)/RN(1,1)
C
      CALL TRCHCK(ICHCK)
      IF(ICHCK.EQ.1) GOTO 4000
C
      CALL TR_EDGE_SELECTOR(1)
      CALL TRCALC(IERR)
      IF(IERR.NE.0) RETURN
      GOTO 2000
C
 4000 NT=NT+1
      T=T+DT
      VSEC=VSEC+VLOOP*DT
      IF(Q0.LT.1.D0) TST=TST+DT
      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         TSL=DT*DBLE(NT)
         CALL LAGLANGE(TSL,RIP,TMU1,RIPU,NTXMAX1,NTUM,IERR)
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
      IF(MDTC.NE.0) THEN
         DO NS=1,NSMAX
         DO NR=1,NRMAX
            ZV(NS,NR) = Z(NS,NR)
         ENDDO
         ENDDO
      ENDIF
C
C     /* Making New Physical Variables */
      CALL TRXTOA
      CALL TR_EDGE_SELECTOR(1)
      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) CALL TR_UFREAD
C      IF(MDLUF.EQ.2.AND.MODEP.EQ.3) CALL TR_UFREAD_S
      IF(MDLUF.EQ.2) CALL TR_UFREAD_S
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
      IF(MDTC.NE.0) THEN
         CALL TRXTOA_AKDW
         CALL TRCFDW_AKDW
      ENDIF
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
      IF(MODELQ.EQ.3.AND.MOD(NT,NTEQIT).EQ.0) THEN
C         CALL TRCONV(L,IERR)
C         WRITE(6,*) "L=",L
         CALL TRSETG
         IF(IERR.NE.0) RETURN
      ENDIF
      IF(IDGLOB.EQ.0) CALL TRGLOB
      IF(NT.LT.NTMAX) GOTO 1000
C
 9000 IF(MDLUF.NE.1.AND.MDLUF.NE.3) RIPS=RIPE
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
      INCLUDE 'trcomm.inc'
C
      COMMON /TRLCL1/ A(NVM,NVM,NRM),B(NVM,NVM,NRM),C(NVM,NVM,NRM)
      COMMON /TRLCL2/ D(NVM,NRM)
      COMMON /TRLCL3/ RD(NEQM,NRM)
      COMMON /TRLCL4/ PPA(NEQM,NRM),PPB(NEQM,NRM),PPC(NEQM,NRM)
      COMMON /TMSLC1/ TMU(NTUM),TMU1(NTUM)
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
C
      IF(MDLCD.EQ.0) THEN
         RDPS=2.D0*PI*RMU0*RIP*1.D6/(DVRHOG(NRMAX)*ABRHOG(NRMAX))
      ELSE
         NEQ=1
         NSVN=NSS(NEQ)
         IF(NSVN.EQ.0) THEN
            RDPA=XV(NEQ,NRMAX)
         ELSE
            RDPA=0.D0
         ENDIF
         RLP=RA*(LOG(8.D0*RR/RA)-2.D0)
         RDPS=RDPA-4.D0*PI*PI*RMU0*RA/(RLP*DVRHOG(NRMAX)*ABRHOG(NRMAX))
      ENDIF
      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         IF(NT.EQ.0) THEN
            TSL=DT*DBLE(1)
            CALL LAGLANGE(TSL,RIPL,TMU1,RIPU,NTXMAX1,NTUM,IERR)
            RDPS=2.D0*PI*RMU0*RIPL*1.D6
     &           /(DVRHOG(NRMAX)*ABRHOG(NRMAX))
         ELSE
            TSL=DT*DBLE(NT)
            CALL LAGLANGE(TSL,RIPL,TMU1,RIPU,NTXMAX1,NTUM,IERR)
            RDPS=2.D0*PI*RMU0*RIPL*1.D6
     &           /(DVRHOG(NRMAX)*ABRHOG(NRMAX))
         ENDIF
      ENDIF
C      IF(MODELQ.EQ.3) THEN
C         BPS=BPSEQ*RIP/RIPEQ
C      ENDIF
C
      COEF = AEE**4*1.D20/(3.D0*SQRT(2.D0*PI)*PI*EPS0**2)
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
     &           *AMZ(NSSN1))*DV53
     &           *COULOG(NSSN,NSSN1,RN(NR,1),RT(NR,1))
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
      Y(1,NR)=(1.D0-PRV/TAUB(NR))*YV(1,NR)+PNB(NR)*DT/(RKEV*1.D20)
      Y(2,NR)=(1.D0-PRV/TAUF(NR))*YV(2,NR)+PNF(NR)*DT/(RKEV*1.D20)
      AY(1,NR)=1.D0+ADV/TAUB(NR)
      AY(2,NR)=1.D0+ADV/TAUF(NR)
      IF(MDTC.NE.0) THEN
         DO NS=1,NSMAX
            Z(NS,NR)=(1.D0-PRV/TAUK(NR))*ZV(NS,NR)
     &              +AKDW(NR,NS)*DT/TAUK(NR)
            AZ(NS,NR)=1.D0+ADV/TAUK(NR)
         ENDDO
      ENDIF
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
     &                 *AMZ(NSSN1))*DV53
     &                 *COULOG(NSSN,NSSN1,RN(NR,1),RT(NR,1))
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
         Y(1,NR)=(1.D0-PRV/TAUB(NR))*YV(1,NR)+PNB(NR)*DT/(RKEV*1.D20)
         Y(2,NR)=(1.D0-PRV/TAUF(NR))*YV(2,NR)+PNF(NR)*DT/(RKEV*1.D20)
         AY(1,NR)=1.D0+ADV/TAUB(NR)
         AY(2,NR)=1.D0+ADV/TAUF(NR)
         IF(MDTC.NE.0) THEN
            DO NS=1,NSMAX
               Z(NS,NR)=(1.D0-PRV/TAUK(NR))*ZV(NS,NR)
     &                 +AKDW(NR,NS)*DT/TAUK(NR)
               AZ(NS,NR)=1.D0+ADV/TAUK(NR)
            ENDDO
         ENDIF
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
     &              *AMZ(NSSN1))*DV53
     &              *COULOG(NSSN,NSSN1,RN(NR,1),RT(NR,1))
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
      Y(1,NR)=(1.D0-PRV/TAUB(NR))*YV(1,NR)+PNB(NR)*DT/(RKEV*1.D20)
      Y(2,NR)=(1.D0-PRV/TAUF(NR))*YV(2,NR)+PNF(NR)*DT/(RKEV*1.D20)
      AY(1,NR)=1.D0+ADV/TAUB(NR)
      AY(2,NR)=1.D0+ADV/TAUF(NR)
      IF(MDTC.NE.0) THEN
         DO NS=1,NSMAX
            Z(NS,NR)=(1.D0-PRV/TAUK(NR))*ZV(NS,NR)
     &              +AKDW(NR,NS)*DT/TAUK(NR)
            AZ(NS,NR)=1.D0+ADV/TAUK(NR)
         ENDDO
      ENDIF
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
         X(MV)=RDPS
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
      INCLUDE 'trcomm.inc'
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
      INCLUDE 'trcomm.inc'
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
      INCLUDE 'trcomm.inc'
C
      DO NR=1,NRMAX
         DO NEQ=1,NEQMAX
            NSVN=NSV(NEQ)
            NSSN=NSS(NEQ)
            IF(NSVN.EQ.0) THEN
               XV(NEQ,NR) = RDP(NR)
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
      IF(MDTC.NE.0) THEN
         DO NR=1,NRMAX
            DO NS=1,NSMAX
               ZV(NS,NR) = AKDW(NR,NS)
            ENDDO
         ENDDO
      ENDIF
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
      INCLUDE 'trcomm.inc'
C
      DO NR=1,NRMAX
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            IF(NSVN.EQ.0) THEN
               RDP(NR) = XV(NEQ,NR)
            ELSEIF(NSVN.EQ.1) THEN
               IF(MDLEQN.NE.0) THEN
                  IF(NSSN.EQ.1.AND.MDLEQE.EQ.0) THEN
                     RN(NR,NSSN) = 0.D0
                     DO NEQ1=1,NEQMAX
                        NSSN1=NSS(NEQ1)
                        NSVN1=NSV(NEQ1)
                        IF(NSVN1.EQ.1.AND.NSSN1.NE.1) THEN
                           RN(NR,NSSN) = RN(NR,NSSN)
     &                                  +ABS(PZ(NSSN1))*XV(NEQ1,NR)
                        ENDIF
                     ENDDO
                     RN(NR,NSSN) = RN(NR,NSSN)
     &                    +PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
C     I think the above expression is obsolete.
                  ELSEIF(NSSN.EQ.2.AND.MDLEQE.EQ.2) THEN
                     RN(NR,NSSN) = 0.D0
                     DO NEQ1=1,NEQMAX
                        NSSN1=NSS(NEQ1)
                        NSVN1=NSV(NEQ1)
                        IF(NSVN1.EQ.1.AND.NSSN1.NE.2) THEN
                           RN(NR,NSSN) = RN(NR,NSSN)
     &                                  +ABS(PZ(NSSN1))*XV(NEQ1,NR)
                        ENDIF
                     ENDDO
                     RN(NR,NSSN) = RN(NR,NSSN)
     &                    +PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
                  ELSE
                     RN(NR,NSSN) = XV(NEQ,NR)
                  ENDIF
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
         BP(NR)   = AR1RHOG(NR)*RDP(NR)/RR
      ENDDO
C
      SUM=0.D0
      DO NR=1,NROMAX
         SUM=SUM+RDP(NR)*DR
         RPSI(NR)=SUM
      ENDDO
C
      N=NRM
      CALL PLDATA_SETP(N,RN,RT,RU)
C
      ENTRY TRXTOA_AKDW
C
      IF(MDTC.NE.0) THEN
         DO NR=1,NRMAX
            DO NS=1,NSMAX
               AKDW(NR,NS) = ZV(NS,NR)
            ENDDO
         ENDDO
      ENDIF
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
      INCLUDE 'trcomm.inc'
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
      GTS(NGST) = GUCLIP(T)
C
  100 GVT(NGST,K+89) = GUCLIP(RT(L,1))
      GVT(NGST,K+90) = GUCLIP(RT(L,2))
      GVT(NGST,K+91) = GUCLIP(RT(L,3))
      GVT(NGST,K+92) = GUCLIP(RT(L,4))
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
      INCLUDE 'trcomm.inc'
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
      INCLUDE 'trcomm.inc'
      COMMON /TRLCL1/ A(NVM,NVM,NRM),B(NVM,NVM,NRM),C(NVM,NVM,NRM)
      COMMON /TRLCL2/ D(NVM,NRM)
      COMMON /TRLCL3/ RD(NEQM,NRM)
      DIMENSION F2C(2),SIG(4),FCB(4)
C
      DV11=DVRHO(NR)
      DV23=DVRHO(NR)**(2.D0/3.D0)
      DV53=DVRHO(NR)**(5.D0/3.D0)
      CC=1.5D0
C
      DO NEQ=1,NEQMAX  ! *** NEQ ***
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
               FA(NMK,NSW)=DVRHO(NR+(NMK-2))*AR1RHO(NR+(NMK-2))/DR
               FB(NMK,NSW)=DVRHO(NR+(NMK-2))*AR2RHO(NR+(NMK-2))/(DR*DR)
            ELSE
               FA(NMK,NSW)=0.D0
               FB(NMK,NSW)=0.D0
            ENDIF
            IF(NSW.NE.3) THEN
               FC(NMK,NSW)=0.5D0*(FCB(NMK)+FCB(NMK+1))/(DR*DR)
            ELSE
               FC(NMK,NSW)=0.D0
            ENDIF
         ENDDO
C     
         SIG(1)= 2.D0
         SIG(2)= 2.D0
         SIG(3)=-2.D0
         SIG(4)=-2.D0
         DO NMK=1,4  ! *** NMK ***
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
            IF(MDDIAG.EQ.0) THEN  ! *** MDDIAG ***
C
            IF(NSVN.EQ.0) THEN
               IF(NSW.NE.3) THEN
                  F2C(NJ)=ETA(NRJ+1)*TTRHO(NRJ+1)
     &                   /(RMU0*ARRHO(NRJ+1)*DVRHO(NRJ+1))
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
                  IF(MDLFLX.EQ.0) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*(AVK(NRJ,NSSN)
     &                                  +AV(NRJ,NSSN)*CC)*DV23
                  ELSE
                  VV(NEQ,NEQ  ,NMK,NSW)=(FA(NI,NSW)*AVK(NRJ,NSSN)
     &                                  +CC*RGFLX(NRJ,NSSN)
     &                                  /(RNV(NRJ,NSSN)*DR))*DV23
                  ENDIF
                  DD(NEQ,NEQ  ,NMK,NSW)= FB(NI,NSW)*AK(NRJ,NSSN)*DV23
                  VV(NEQ,NEQ-1,NMK,NSW)= 0.D0*DV23
                  DD(NEQ,NEQ-1,NMK,NSW)= FB(NI,NSW)*(AD(NRJ,NSSN)
     &                                  *CC-AK(NRJ,NSSN))
     &                                  *RTV(NRJ,NSSN)*DV23
               ENDIF
            ELSEIF(NSVN.EQ.3) THEN
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*VOID
                  DD(NEQ,NEQ  ,NMK,NSW)= FB(NI,NSW)*VOID
               ENDIF
            ENDIF
C
            ELSE  ! *** MDDIAG ***
C
C     *** OFF-DIAGONAL TRANSPORT ***
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
     &                   /(RMU0*ARRHO(NRJ+1)*DVRHO(NRJ+1))
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
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)
     &                 *(AV(NRJ,NSSN)+RGFLS(NRJ,5,NSSN)/RNV(NRJ,NSSN))
                  DO NEQ1=1,NEQLMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)
     &                       *RNV(NRJ,NSSN)*ADLD(NRJ,NSSN1,NSSN)
     &                       /RNV(NRJ,NSSN1)
                     ELSEIF(NSVN1.EQ.2) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)
     &                       *RNV(NRJ,NSSN)*ADLP(NRJ,NSSN1,NSSN)
     &                       /RPV(NRJ,NSSN1)
                     ENDIF
                  ENDDO
               ENDIF
               ENDIF
            ELSEIF(NSVN.EQ.2) THEN
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  IF(MDLFLX.EQ.0) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*DV23
     &                    *( AVK(NRJ,NSSN)+AV(NRJ,NSSN)*CC
     &                      +RQFLS(NRJ,5,NSSN)   /RPV(NRJ,NSSN)
     &                      +RGFLS(NRJ,5,NSSN)*CC/RNV(NRJ,NSSN))
                  ELSE
                  VV(NEQ,NEQ  ,NMK,NSW)= DV23*(FA(NI,NSW)
     &                    *( AVK(NRJ,NSSN)
     &                      +RQFLS(NRJ,5,NSSN)   / RPV(NRJ,NSSN)
     &                      +RGFLS(NRJ,5,NSSN)*CC/ RNV(NRJ,NSSN))
     &                     +(RGFLX(NRJ,NSSN)  *CC/(RNV(NRJ,NSSN)*DR)))
                  ENDIF
                  DO NEQ1=1,NEQLMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)*DV23
     &                       * RPV(NRJ,NSSN)
     &                       *( AKLD(NRJ,NSSN1,NSSN)
     &                         +ADLD(NRJ,NSSN1,NSSN)*CC)
     &                       / RNV(NRJ,NSSN1)
                     ELSEIF(NSVN1.EQ.2) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)*DV23
     &                       * RPV(NRJ,NSSN)
     &                       *( AKLP(NRJ,NSSN1,NSSN)
     &                         +ADLP(NRJ,NSSN1,NSSN)*CC)
     &                       / RPV(NRJ,NSSN1)
                     ENDIF
                  ENDDO
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
            ENDIF  ! *** MDDIAG ***
         ENDDO  ! *** NMK ***
C
      ENDDO  ! *** NEQ ***
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
         IF(NR.EQ.NRMAX.AND.MDDIAG.EQ.0) THEN
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
C
C     *** OFF-DIAGONAL TRANSPORT ***
C     *
         ELSEIF(NR.EQ.NRMAX.AND.MDDIAG.NE.0) THEN
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
               D(NEQ,NR) = ETA(NR  )*BB/(TTRHO(NR  )*ARRHO(NR  )*DR)
     &                    *(AJ(NR  )-AJOH(NR  ))
     &                    -ETA(NR+1)*BB/(TTRHO(NR+1)*ARRHO(NR+1)*DR)
     &                    *(AJ(NR+1)-AJOH(NR+1))
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
      INCLUDE 'trcomm.inc'
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
      INCLUDE 'trcomm.inc'
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
C
C     ***********************************************************
C
C           LOCAL VARIABLES FOR DENSITY, TEMPERATURE, PRESSURE
C
C     ***********************************************************
C
      FUNCTION RNV(NR,NS)
C
      INCLUDE 'trcomm.inc'
C
      IF(NR.EQ.NRMAX) THEN
         RNV=PNSS(NS)
      ELSEIF(NR.NE.0.AND.NR.NE.NRMAX) THEN
         RNV=0.5D0*(RN(NR,NS)+RN(NR+1,NS))
      ENDIF
C
      RETURN
      END
C
      FUNCTION RTV(NR,NS)
C
      INCLUDE 'trcomm.inc'
C
      IF(NR.EQ.NRMAX) THEN
         RTV=PTS(NS)
      ELSEIF(NR.NE.0.AND.NR.NE.NRMAX) THEN
         RTV=0.5D0*(RT(NR,NS)+RT(NR+1,NS))
      ENDIF
C
      RETURN
      END
C
      FUNCTION RPV(NR,NS)
C
      INCLUDE 'trcomm.inc'
C
      IF(NR.EQ.NRMAX) THEN
         RPV=PNSS(NS)*PTS(NS)
      ELSEIF(NR.NE.0.AND.NR.NE.NRMAX) THEN
         RPV=0.5D0*(RN(NR,NS)*RT(NR,NS)+RN(NR+1,NS)*RT(NR+1,NS))
      ENDIF
C
      RETURN
      END

