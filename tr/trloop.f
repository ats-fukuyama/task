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
C
      DIMENSION XX(MLM),YY(2,NRM)
C
      NT=0
      RIP=RIPS
      IF(NTMAX.NE.0) DIP=(RIPE-RIPS)/NTMAX
      ICHCK=0
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (MDLEQ.EQ.0) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         NQM=NVM
         MWMAX=MWM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF (MDLEQ.EQ.10) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         NQM=NSM+1
         MWMAX=4*NQM-1
      ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (MDLEQ.EQ.0) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF (MDLEQ.EQ.10) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO NR=1,NRMAX
         RN(NR,1) = 0.D0
         DO NS=2,NSM
            RN(NR,1)  = RN(NR,1)+PZ(NS)*RN(NR,NS)
         ENDDO
         RN(NR,1) = RN(NR,1)+PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
         DO NS=1,NSM
            IF(NS.NE.NSM) THEN
               RT(NR,NS) = 0.5D0*(XV(NS,NR)+X(NQM*(NR-1)+NS))/RN(NR,NS)
            ELSE
            IF(RN(NR,NSM).LT.1.D-70) THEN
               RT(NR,NS) = 0.D0
            ELSE
               RT(NR,NS) = 0.5D0*(XV(NS,NR)+X(NQM*(NR-1)+NS))/RN(NR,NS)
            ENDIF
            ENDIF
         ENDDO
         BP(NR)    = 0.5D0*(XV(NQM,NR)+X(NQM*NR))
         RW(NR,1)  = 0.5D0*(YV(1,NR)+Y(1,NR))
         RW(NR,2)  = 0.5D0*(YV(2,NR)+Y(2,NR))
      ENDDO
      ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C 
      CALL TRCHCK(ICHCK)
      IF(ICHCK.EQ.1) GOTO 4000
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
      IF(MOD(NT,NTEQIT).EQ.0) THEN
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
C      DIMENSION DNDR(NRM,NSM,3)
C
      RKAPX=(RKAP-1.D0)/(RKAP+1.D0)
      FKAP=0.5D0*(RKAP+1.D0)
     &          *(1.D0+RKAPX/4.D0+RKAPX*RKAPX/64.D0)
      DO NS=1,NSM
         AMZ(NS)=PA(NS)*AMM/PZ(NS)**2
      ENDDO
C
      IF(MDLCD.EQ.0) THEN
         BPS= AMYU0*RIP*1.D6/(2.D0*PI*RA*FKAP)
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
      DO NW=1,NVM-1
      DO NV=1,NVM-1
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
         FAM=ETA(NR  )/(RR*TTRHO(NR  )*ARRHO(NR  ))
         FAP=ETA(NR+1)/(RR*TTRHO(NR+1)*ARRHO(NR+1))
         FBM=RR*TTRHO(NR  )**2/(AMYU0*DVRHO(NR  ))
         FBP=RR*TTRHO(NR+1)**2/(AMYU0*DVRHO(NR+1))
         FC1=DVRHO(NR-1)*ABRHO(NR-1)/TTRHO(NR-1)
         FC2=DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  )
         FC3=DVRHO(NR+1)*ABRHO(NR+1)/TTRHO(NR+1)
         IF(NR.EQ.NRMAX-1) THEN
            FC4=2.D0*FC3-FC2
         ELSE
            FC4=DVRHO(NR+2)*ABRHO(NR+2)/TTRHO(NR+2)
         ENDIF
         F1M=FAM*AR1RHO(NR  )*BB/DR
         F1P=FAP*AR1RHO(NR+1)*BB/DR
         F2M=FAM*FBM/(DR*DR)
         F2P=FAP*FBP/(DR*DR)
         FCM=0.5D0*(FC1+FC2)
         FC0=0.5D0*(FC2+FC3)
         FCP=0.5D0*(FC3+FC4)
         IF (MODELG.EQ.0) THEN
            FVL = FKAP/RKAP
         ELSEIF (MODELG.EQ.3) THEN
            FVL = 1.D0
         ENDIF
C
         A(NN,NN,NR) = F2M*FCM*FVL
         B(NN,NN,NR) =-F2M*FC0*FVL
     &                -F2P*FC0*FVL
         C(NN,NN,NR) = F2P*FCP*FVL
         D(NN   ,NR) = F1M*(AJ(NR  )-AJOH(NR  ))
     &                -F1P*(AJ(NR+1)-AJOH(NR+1))
      ENDDO
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (MDLEQ.EQ.0) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO NR=2,NRMAX-1
C
         DV11=DVRHO(NR)
         DV23=DVRHO(NR)**(2.D0/3.D0)
         DV53=DVRHO(NR)**(5.D0/3.D0)
C
      NSW=2
      CALL TRIPTC(NSW,NR)
C
      DO NV=1,NQM-1
      DO NW=1,NQM-1
      DO NO=1,2
         IF (NV.LE.NSM) THEN
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
         ELSEIF (NW.GT.NSM) THEN
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
     &                         *DV23
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
     &                         *DV23
         ELSE
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
     &                         *DV23
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
     &                   *0.5D0*(RT(NR+(NO-2),NW)+RT(NR+(NO-1),NW))
     &                         *DV23
         ENDIF
      ENDDO
C      IF (NV.GT.NSM.AND.NV.EQ.NW) write(6,*) NR,DI(NV,NW,2,NSW)*2.D0
C     &     /(3.D0*DV53)
C      IF (NV.GT.NSM.AND.NV.EQ.NW+NSM) write(6,*) NR,DI(NV,NW,2,NSW)
C     &     /DV53
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF (MDLEQ.EQ.10) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO NR=2,NRMAX-1
C
         DV11=DVRHO(NR)
         DV23=DVRHO(NR)**(2.D0/3.D0)
         DV53=DVRHO(NR)**(5.D0/3.D0)
C
      NSW=2
      CALL TRIPTC(NSW,NR)
C      CALL RKDNDR(DNDR,NR)
C
      DO NV=1,NQM-1
      DO NW=1,NQM-1
      DO NO=1,2
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
     &                         *DV23
         IF (NR.EQ.2) THEN
         VI(NV,NW,NO,NSW)=VI(NV,NW,NO,NSW)
     &                   +0.5D0*( VV(NSM+NV,NSM+NW,NO  ,NSW)
     &                           +VV(NSM+NV,NSM+NW,NO+2,NSW))
     &                  /(0.5D0*(RN(NR+(NO-2),NW)+RN(NR+(NO-1),NW)))
C     &                   *0.5D0*(DNDR(NR+(NO-2),NW,NO)
C     &                          +DNDR(NR+(NO-1),NW,NO+1))
     &                         *(RN(NR+ NO   ,NW)-RN(NR+(NO-2),NW))
     &                         *DV23
         ELSEIF (NR.EQ.NRMAX-1) THEN
         VI(NV,NW,NO,NSW)=VI(NV,NW,NO,NSW)
     &                   +0.5D0*( VV(NSM+NV,NSM+NW,NO  ,NSW)
     &                           +VV(NSM+NV,NSM+NW,NO+2,NSW))
     &                  /(0.5D0*(RN(NR+(NO-2),NW)+RN(NR+(NO-1),NW)))
C     &                   *0.5D0*(DNDR(NR+(NO-2),NW,NO)
C     &                          +DNDR(NR+(NO-1),NW,NO+1))
     &                         *(RN(NR+(NO-1),NW)-RN(NR+(NO-3),NW))
     &                         *DV23
         ELSE
          VI(NV,NW,NO,NSW)=VI(NV,NW,NO,NSW)
     &                    +0.5D0*( VV(NSM+NV,NSM+NW,NO  ,NSW)
     &                            +VV(NSM+NV,NSM+NW,NO+2,NSW))
     &                   /(0.5D0*(RN(NR+(NO-2),NW)+RN(NR+(NO-1),NW)))
C     &                   *0.5D0*(DNDR(NR+(NO-2),NW,NO)
C     &                          +DNDR(NR+(NO-1),NW,NO+1))
     &                    *0.5D0*(RN(NR+ NO   ,NW)+RN(NR+(NO-1),NW)
     &                           -RN(NR+(NO-2),NW)-RN(NR+(NO-3),NW))
     &                          *DV23
         ENDIF
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
     &                         *DV23
      ENDDO
         A(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)+DI(NV,NW,1,NSW)
         B(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)-DI(NV,NW,1,NSW)
     &                -0.5D0*VI(NV,NW,2,NSW)-DI(NV,NW,2,NSW)
         C(NV,NW,NR) =-0.5D0*VI(NV,NW,2,NSW)+DI(NV,NW,2,NSW)
      ENDDO
      ENDDO
C
      DO NS=1,NSM
         D(NS,NR) = (PIN(NR,NS)/(RKEV*1.D20))*DV53
C
      DO NS1=1,NSM
         IF(NS1.NE.NS) THEN
            C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1))*DV53
     &             *1.5D0
            B(NS,NS, NR)=B(NS,NS, NR)-RN(NR,NS1)*C1
            B(NS,NS1,NR)=B(NS,NS1,NR)+RN(NR,NS )*C1
         ENDIF
      ENDDO
      ENDDO
C
      ENDDO
      ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      NR=1
         FAM=ETA(NR  )/(RR*TTRHO(NR  )*ARRHO(NR  ))
         FAP=ETA(NR+1)/(RR*TTRHO(NR+1)*ARRHO(NR+1))
         FBM=RR*TTRHO(NR  )**2/(AMYU0*DVRHO(NR  ))
         FBP=RR*TTRHO(NR+1)**2/(AMYU0*DVRHO(NR+1))
         FC1=0.D0
         FC2=DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  )
         FC3=DVRHO(NR+1)*ABRHO(NR+1)/TTRHO(NR+1)
         FC4=DVRHO(NR+2)*ABRHO(NR+2)/TTRHO(NR+2)
         F1M=FAM*AR1RHO(NR  )*BB/DR
         F1P=FAP*AR1RHO(NR+1)*BB/DR
         F2M=FAM*FBM/(DR*DR)
         F2P=FAP*FBP/(DR*DR)
         FCM=0.5D0*(FC1+FC2)
         FC0=0.5D0*(FC2+FC3)
         FCP=0.5D0*(FC3+FC4)
         IF (MODELG.EQ.0) THEN
            FVL = FKAP/RKAP
         ELSEIF (MODELG.EQ.3) THEN
            FVL = 1.D0
         ENDIF
C
C         A(NN,NN,NR) = F2M*FCM*FVL
         A(NN,NN,NR) = 0.D0*FVL
         B(NN,NN,NR) =-F2M*FC0*FVL
     &                -F2P*FC0*FVL
         C(NN,NN,NR) = F2P*FCP*FVL
         D(NN   ,NR) = F1M*(AJ(NR  )-AJOH(NR  ))
     &                -F1P*(AJ(NR+1)-AJOH(NR+1))
C     
         DV11=DVRHO(NR)
         DV23=DVRHO(NR)**(2.D0/3.D0)
         DV53=DVRHO(NR)**(5.D0/3.D0)
C
      NSW=1
      CALL TRIPTC(NSW,NR)
C
      NO=2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (MDLEQ.EQ.0) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO NV=1,NQM-1
      DO NW=1,NQM-1
         IF (NV.LE.NSM) THEN
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
         ELSEIF (NW.GT.NSM) THEN
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
     &                         *DV23
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
     &                         *DV23
         ELSE
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
     &                         *DV23
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
     &                   *0.5D0*(RT(NR+(NO-2),NW)+RT(NR+(NO-1),NW))
     &                         *DV23
         ENDIF
C      IF (NV.GT.NSM.AND.NV.EQ.NW) write(6,*) DI(NV,NW,2,NSW)*2.D0
C     &     /(3.D0*DV53)
C      IF (NV.GT.NSM.AND.NV.EQ.NW+NSM) write(6,*) DI(NV,NW,2,NSW)
C     &     /DV53
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF (MDLEQ.EQ.10) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C      CALL RKDNDR(DNDR,NR)
      DO NV=1,NQM-1
      DO NW=1,NQM-1
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
     &                         *DV23
         VI(NV,NW,NO,NSW)=VI(NV,NW,NO,NSW)
     &                   +0.5D0*( VV(NSM+NV,NSM+NW,NO  ,NSW)
     &                           +VV(NSM+NV,NSM+NW,NO+2,NSW))
     &                  /(0.5D0*(RN(NR+(NO-2),NW)+RN(NR+(NO-1),NW)))
C     &                   *0.5D0*(DNDR(NR+(NO-2),NW,NO)
C     &                          +DNDR(NR+(NO-1),NW,NO+1))
     &                         *(RN(NR+ NO   ,NW)-RN(NR+(NO-2),NW))
     &                         *DV23
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
     &                         *DV23
C
         A(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)+DI(NV,NW,1,NSW)
         B(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)-DI(NV,NW,1,NSW)
     &                -0.5D0*VI(NV,NW,2,NSW)-DI(NV,NW,2,NSW)
         C(NV,NW,NR) =-0.5D0*VI(NV,NW,2,NSW)+DI(NV,NW,2,NSW)
      ENDDO
      ENDDO
C
      DO NS=1,NSM
         D(NS,NR) = (PIN(NR,NS)/(RKEV*1.D20))*DV53
C
      DO NS1=1,NSM
         IF(NS1.NE.NS) THEN
            C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1))*DV53
     &             *1.5D0
            B(NS,NS, NR)=B(NS,NS, NR)-RN(NR,NS1)*C1
            B(NS,NS1,NR)=B(NS,NS1,NR)+RN(NR,NS )*C1
         ENDIF
      ENDDO
      ENDDO
      ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      NR=NRMAX
         A(NN,NN,NR) = 0.D0
         B(NN,NN,NR) = 0.D0
         C(NN,NN,NR) = 0.D0
         D(NN   ,NR) = 0.D0
C
         DV11=DVRHO(NR)
         DV23=DVRHO(NR)**(2.D0/3.D0)
         DV53=DVRHO(NR)**(5.D0/3.D0)
C
      NSW=3
      CALL TRIPTC(NSW,NR)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (MDLEQ.EQ.0) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO NV=1,NQM-1
      DO NW=1,NQM-1
         NO=1
         IF (NV.LE.NSM) THEN
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
         ELSEIF (NW.GT.NSM) THEN
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
     &                         *DV23
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
     &                         *DV23
         ELSE
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
     &                         *DV23
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
     &                   *0.5D0*(RT(NR+(NO-2),NW)+RT(NR+(NO-1),NW))
     &                         *DV23
         ENDIF
C
         NO=2
         IF (NV.LE.NSM) THEN
         VI(NV,NW,NO,NSW)=VV(NV,NW,NO  ,NSW)
         DI(NV,NW,NO,NSW)=DD(NV,NW,NO  ,NSW)
         ELSEIF (NW.GT.NSM) THEN
         VI(NV,NW,NO,NSW)=VV(NV,NW,NO  ,NSW)        *DV23
         DI(NV,NW,NO,NSW)=DD(NV,NW,NO  ,NSW)        *DV23
         ELSE
         VI(NV,NW,NO,NSW)=VV(NV,NW,NO  ,NSW)        *DV23
         DI(NV,NW,NO,NSW)=DD(NV,NW,NO  ,NSW)*PTS(NW)*DV23
         ENDIF
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF (MDLEQ.EQ.10) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C      CALL RKDNDR(DNDR,NR)
      DO NV=1,NQM-1
      DO NW=1,NQM-1
         NO=1
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO  ,NSW)+VV(NV,NW,NO+2,NSW))
     &                         *DV23
         VI(NV,NW,NO,NSW)=VI(NV,NW,NO,NSW)
     &                   +0.5D0*( VV(NSM+NV,NSM+NW,NO  ,NSW)
     &                           +VV(NSM+NV,NSM+NW,NO+2,NSW))
     &                  /(0.5D0*(RN(NR+(NO-2),NW)+RN(NR+(NO-1),NW)))
C     &                   *0.5D0*(DNDR(NR+(NO-2),NW,NO)
C     &                          +DNDR(NR+(NO-1),NW,NO+1))
     &                         *(RN(NR+(NO-1),NW)-RN(NR+(NO-3),NW))
     &                         *DV23
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO  ,NSW)+DD(NV,NW,NO+2,NSW))
     &                         *DV23
C
         NO=2
         VI(NV,NW,NO,NSW)=VV(NV,NW,NO  ,NSW)*DV23
         VI(NV,NW,NO,NSW)=VI(NV,NW,NO,NSW)
     &                   +VV(NSM+NV,NSM+NW,NO  ,NSW)/PNSS(NW)
C     &                   *DNDR(NR+(NO-2),NW,NO)
     &                   *(2.D0*PNSS(NW)-RN(NR+(NO-2),NW)
     &                                  -RN(NR+(NO-3),NW))     
     &                   *DV23
         DI(NV,NW,NO,NSW)=DD(NV,NW,NO  ,NSW)*DV23
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
         D(NS,NR) = (PIN(NR,NS)/(RKEV*1.D20)  )*DV53
     &             +(-VI(NS,NS,2,NSW)+DI(NS,NS,2,NSW)*2.D0)
     &             *PNSS(NS)*PTS(NS)
C
      DO NS1=1,NSM
         IF(NS1.NE.NS) THEN
            C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1))*DV53
     &             *1.5D0
            B(NS,NS, NR)=B(NS,NS, NR)-RN(NR,NS1)*C1
            B(NS,NS1,NR)=B(NS,NS1,NR)+RN(NR,NS )*C1
         ENDIF
      ENDDO
      ENDDO
      ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C      FADV=0.5D0
      FADV=1.0D0
C
      PRV=(1.D0-FADV)*DT
      ADV=FADV*DT
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (MDLEQ.EQ.0) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF (MDLEQ.EQ.10) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO NR=1,NRMAX
         DV53=DVRHO(NR)**(5.D0/3.D0)
         DO NS=1,NSM
            X(NQM*(NR-1)+NS) = 1.5D0*DV53*XV(NS,NR)+DT*D(NS,NR)
         ENDDO
         NS=NQM
         X(NQM*(NR-1)+NS) = XV(NS,NR)+DT*D(NS,NR)
      ENDDO
      ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (MDLEQ.EQ.0) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF (MDLEQ.EQ.10) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO NR=1,NRMAX
         DV53=DVRHO(NR)**(5.D0/3.D0)
         DO NS=1,NSM
            AX(2*NQM,NQM*(NR-1)+NS) 
     &    = AX(2*NQM,NQM*(NR-1)+NS) + 1.5D0*DV53
         ENDDO
            AX(2*NQM,NQM*(NR-1)+NN) 
     &    = AX(2*NQM,NQM*(NR-1)+NN) + 1.D0
      ENDDO
      ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (MDLEQ.EQ.0) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO NR=1,NRMAX
      DO NS=1,NSM
         XV(    NS,NR) = RN(NR,NS)
         XV(NSM+NS,NR) = RN(NR,NS)*RT(NR,NS)
      ENDDO
         XV(NQM,NR) = BP(NR)
         YV(  1,NR) = RW(NR,1)
         YV(  2,NR) = RW(NR,2)
      ENDDO
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF (MDLEQ.EQ.10) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO NR=1,NRMAX
      DO NS=1,NSM
         XV(NS,NR) = RN(NR,NS)*RT(NR,NS)
      ENDDO
         XV(NQM,NR) = BP(NR)
         YV(  1,NR) = RW(NR,1)
         YV(  2,NR) = RW(NR,2)
      ENDDO
      ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (MDLEQ.EQ.0) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF (MDLEQ.EQ.10) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO NR=1,NRMAX
         RN(NR,1) = 0.D0
         DO NS=2,NSM
            RN(NR,1) = RN(NR,1)+PZ(NS)*RN(NR,NS)
         ENDDO
         RN(NR,1) = RN(NR,1)+PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
         DO NS=1,NSM
         IF(NS.NE.NSM) THEN
            RT(NR,NS) = XV(NS,NR)/RN(NR,NS)
         ELSE
            IF(RN(NR,NSM).LT.1.D-70) THEN
               RT(NR,NS) = 0.D0
            ELSE
               RT(NR,NS) = XV(NS,NR)/RN(NR,NS)
            ENDIF
         ENDIF
      ENDDO
      BP(NR)   = XV(NQM,NR)
      RW(NR,1) = YV(1,NR)
      RW(NR,2) = YV(2,NR)
      ENDDO
      ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
      SUBROUTINE TRIPTC(NSW,NR)
C
      INCLUDE 'trcomm.h'
C
      DIMENSION SUM(3),DDM(NSM),DNDR(NRM,NSM,3)
C
      DO NS=1,NSM
         RTM(NS)=RT(NR,NS)*RKEV/(PA(NS)*AMM)
      ENDDO
C
      IF (MODELG.EQ.0) THEN
         FVL = FKAP/RKAP
      ELSEIF (MODELG.EQ.3) THEN
         FVL = 1.D0
      ENDIF
      DO NMK=1,3
         IF (NSW.EQ.2.OR.NSW.NE.NMK) THEN
            FA(NMK,NSW)=DVRHO(NR+(NMK-2))*AR1RHO(NR+(NMK-2))*FVL/DR
            FB(NMK,NSW)=DVRHO(NR+(NMK-2))*AR2RHO(NR+(NMK-2))*FVL/(DR*DR)
         ELSE
            FA(NMK,NSW)=0.D0
            FB(NMK,NSW)=0.D0
         ENDIF
      ENDDO
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (MDLEQ.EQ.0) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (MDCHG.EQ.0) THEN
      DO NS=1,NSM
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
         IF (NSW.EQ.2.OR.NSW.NE.NI.AND.NR+(NJ-2).NE.0) THEN
            VV(    NS,    NS,NMK,NSW) = FA(NI,NSW)*AV(NR+(NJ-2),NS)
            DD(    NS,    NS,NMK,NSW) = FB(NI,NSW)*AD(NR+(NJ-2),NS)
            VV(NSM+NS,NSM+NS,NMK,NSW) = FA(NI,NSW)*(AVK(NR+(NJ-2),NS)
     &                                 +AV(NR+(NJ-2),NS)*1.5D0)
            DD(NSM+NS,NSM+NS,NMK,NSW) = FB(NI,NSW)*AK(NR+(NJ-2),NS)
            VV(NSM+NS,    NS,NMK,NSW) = 0.D0
            DD(NSM+NS,    NS,NMK,NSW) = FB(NI,NSW)*(AD(NR+(NJ-2),NS)
     &                                       *1.5D0-AK(NR+(NJ-2),NS))
         ENDIF
      ENDDO
      ENDDO

C   /////////////////////////////////////////////////////////////////////
C
      ELSEIF (MDCHG.EQ.1) THEN
C
      DO NMK=1,3
C     ** BEFORE **
      IF (NSW.EQ.2.OR.NSW.NE.NMK) THEN
      DO NS=1,NSM
         VV(    NS,    NS,NMK,NSW) = FA(NMK,NSW)*AV(NR+(NMK-2),1)
         VV(NSM+NS,NSM+NS,NMK,NSW) = FA(NMK,NSW)*(AVK(NR+(NMK-2),NS)
     &                       +AV(NR+(NMK-2),NS)*1.5D0)*RN(NR+(NMK-2),NS)
         DD(NSM+NS,NSM+NS,NMK,NSW) = FB(NMK,NSW)*RN(NR+(NMK-2),NS)
     &                         *AK(NR+(NMK-2),NS)
         VV(NSM+NS,    NS,NMK,NSW) = 0.D0
         DD(NSM+NS,    NS,NMK,NSW) = FB(NMK,NSW)*1.5D0*RT(NR+(NMK-2),NS)
     &                         *AD(NR+(NMK-2),NS)
      ENDDO
C
      DO NW=1,NSM
      DO NS=1,NSM
         DD(    NS,    NW,NMK,NSW) = FB(NMK,NSW)*AD(NR+(NMK-2),1)
     &        *RN(NR+(NMK-2),NS)*RT(NR+(NMK-2),NW)
     &        /(RN(NR+(NMK-2),1)*RT(NR+(NMK-2),1))
         DD(    NS,NSM+NW,NMK,NSW) = FB(NMK,NSW)*AD(NR+(NMK-2),1)
     &        *3.D0*RN(NR+(NMK-2),NS)*RN(NR+(NMK-2),NW)
     &        /(2.D0*RN(NR+(NMK-2),1)*RT(NR+(NMK-2),1))
      ENDDO
      ENDDO
      ENDIF
C     ** AFTER **
      DO NS=1,NSM
         VV(NSM+NS,NSM+NS,NMK,NSW) = VV(NSM+NS,NSM+NS,NMK,NSW)
     &                              /RN(NR+(NMK-2),NS)
         DD(NSM+NS,NSM+NS,NMK,NSW) = DD(NSM+NS,NSM+NS,NMK,NSW)
     &                              /RN(NR+(NMK-2),NS)
         DD(NSM+NS,    NS,NMK,NSW) = DD(NSM+NS,    NS,NMK,NSW)
     &                    -(RT(NR+(NMK-2),NS)/RN(NR+(NMK-2),NS))
     &                    *(DD(NSM+NS,NSM+NS,NMK,NSW)*RN(NR+(NMK-2),NS))
      ENDDO
C
      DO NW=1,NSM
      DO NS=1,NSM
         DD(    NS,    NW,NMK,NSW) = DD(    NS,    NW,NMK,NSW)
     &                              -DD(    NS,NSM+NW,NMK,NSW)
     &                            *(RT(NR+(NMK-2),NW)/RN(NR+(NMK-2),NW))
         DD(    NS,NSM+NW,NMK,NSW) = DD(    NS,NSM+NW,NMK,NSW)
     &                              /RN(NR+(NMK-2),NW)
      ENDDO
      ENDDO
C
      ENDDO
C
C   /////////////////////////////////////////////////////////////////////
C
      ELSEIF (MDCHG.EQ.2) THEN
C
      DO NMK=1,3
         SUM(NMK)=0.D0
      ENDDO
      DO NMK=1,3
      DO NS=2,NSM
         SUM(NMK)=SUM(NMK)+PZ(NS)*RN(NR+(NMK-2),NS)*RT(NR+(NMK-2),NS)
      ENDDO
      ENDDO
C
      DO NMK=1,3
C     ** BEFORE **
      IF (NSW.EQ.2.OR.NSW.NE.NMK) THEN
      DO NS=1,NSM
         VV(    NS,    NS,NMK,NSW) = FA(NMK,NSW)*AV(NR+(NMK-2),NS)
         DD(    NS,    NS,NMK,NSW) = FB(NMK,NSW)*AD(NR+(NMK-2),NS)
     &                              *(1.D0+SUM(NMK)/(RN(NR+(NMK-2),1)
     &                              *RT(NR+(NMK-2),1)))
         VV(    NS,NSM+NS,NMK,NSW) = FA(NMK,NSW)*AV(NR+(NMK-2),NS)
         VV(NSM+NS,NSM+NS,NMK,NSW) = FA(NMK,NSW)*(AVK(NR+(NMK-2),NS)
     &                              +AV(NR+(NMK-2),NS)*1.5D0)
     &                              *RN(NR+(NMK-2),NS)
         DD(NSM+NS,NSM+NS,NMK,NSW) = FB(NMK,NSW)*RN(NR+(NMK-2),NS)
     &                              *AK(NR+(NMK-2),NS)
         VV(NSM+NS,    NS,NMK,NSW) = 0.D0
         DD(NSM+NS,    NS,NMK,NSW) = FB(NMK,NSW)*1.5D0*RT(NR+(NMK-2),NS)
     &                              *AD(NR+(NMK-2),NS)
         DDM(NS) = DD(    NS,    NS,NMK,NSW)
      ENDDO
C
      DO NW=1,NSM
      DO NS=1,NSM
      DD(    NS,NSM+NW,NMK,NSW) = FB(NMK,NSW)*AD(NR+(NMK-2),NS)
     &        *3.D0*RN(NR+(NMK-2),NS)*ABS(PZ(NW))*RN(NR+(NMK-2),NW)
     &        /(2.D0*RN(NR+(NMK-2),1)*RT(NR+(NMK-2),1))
      ENDDO
      ENDDO
      ENDIF
C     ** AFTER **
      DO NS=1,NSM
         VV(NSM+NS,NSM+NS,NMK,NSW) = VV(NSM+NS,NSM+NS,NMK,NSW)
     &                              /RN(NR+(NMK-2),NS)
         DD(NSM+NS,NSM+NS,NMK,NSW) = DD(NSM+NS,NSM+NS,NMK,NSW)
     &                              /RN(NR+(NMK-2),NS)
         DD(NSM+NS,    NS,NMK,NSW) = DD(NSM+NS,    NS,NMK,NSW)
     &                            -(RT(NR+(NMK-2),NS)/RN(NR+(NMK-2),NS))
     &                    *(DD(NSM+NS,NSM+NS,NMK,NSW)*RN(NR+(NMK-2),NS))
         DD(    NS,     1,NMK,NSW) =-DD(    NS,NSM+ 1,NMK,NSW)
     &                             *(RT(NR+(NMK-2),1)/RN(NR+(NMK-2),1))
         DD(    NS,NSM +1,NMK,NSW) = DD(    NS,NSM +1,NMK,NSW)
     &                              /RN(NR+(NMK-2),1)
      ENDDO
C
      DO NW=2,NSM
      DO NS=1,NSM
         DD(    NS,    NW,NMK,NSW) =-DD(    NS,NSM+NW,NMK,NSW)
     &                            *(RT(NR+(NMK-2),NW)/RN(NR+(NMK-2),NW))
         DD(    NS,NSM+NW,NMK,NSW) = DD(    NS,NSM+NW,NMK,NSW)
     &                              /RN(NR+(NMK-2),NW)
      ENDDO
      ENDDO
      DO NS=1,NSM
         DD(    NS,    NS,NMK,NSW) = DD(    NS,    NS,NMK,NSW)+DDM(NS)
      ENDDO
C
      ENDDO
      ENDIF
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF (MDLEQ.EQ.10) THEN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      CALL RKDNDR(DNDR,NR)
C
      DO NS=1,NSM
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
      IF (DNDR(NR,NS,NI).NE.0.D0.AND.NR+(NJ-2).NE.0) THEN
         VV(    NS,    NS,NMK,NSW)= FA(NI,NSW)*(AVK(NR+(NJ-2),NS)
     &                                        +  AV(NR+(NJ-2),NS)*1.5D0)
         VV(NSM+NS,NSM+NS,NMK,NSW)=-FA(NI,NSW)*( AD(NR+(NJ-2),NS)*1.5D0
     &                                        -  AK(NR+(NJ-2),NS))
         DD(    NS,    NS,NMK,NSW)= FB(NI,NSW)*  AK(NR+(NJ-2),NS)
      ENDIF
      ENDDO
      ENDDO
C
      ENDIF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      RETURN
      END
C
C     ***********************************************************
C
C           RECKON NUMERATOR OF DN/DR
C
C     ***********************************************************
C
      SUBROUTINE RKDNDR(DNDR,NR)
C
      INCLUDE 'trcomm.h'
C
      DIMENSION DNDR(NRM,NSM,3)
C
      DO NS=1,NSM
      IF (NR.EQ.1) THEN
         NMK=1
         DNDR(NR,NS,NMK) = 0.D0
         DO NMK=2,3
         DNDR(NR,NS,NMK) = RN(NR+(NMK-2)+1,NS)-RN(NR+(NMK-2),NS)
         ENDDO
      ELSEIF (NR.EQ.NRMAX-1) THEN
         DO NMK=1,2
         DNDR(NR,NS,NMK) = RN(NR+(NMK-2)+1,NS)-RN(NR+(NMK-2),NS)
         ENDDO
         NMK=3
         DNDR(NR,NS,NMK) = 2.D0*(PNSS(NS)-RN(NR+(NMK-2),NS))
      ELSEIF (NR.EQ.NRMAX) THEN
         NMK=1
         DNDR(NR,NS,NMK) = RN(NR+(NMK-2)+1,NS)-RN(NR+(NMK-2),NS)
         NMK=2
         DNDR(NR,NS,NMK) = 2.D0*(PNSS(NS)-RN(NR+(NMK-2),NS))
         NMK=3
         DNDR(NR,NS,NMK) = 0.D0
      ELSE
         DO NMK=1,3
         DNDR(NR,NS,NMK) = RN(NR+(NMK-2)+1,NS)-RN(NR+(NMK-2),NS)
         ENDDO
      ENDIF
      ENDDO
C
      RETURN
      END
