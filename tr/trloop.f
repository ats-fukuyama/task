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
      CALL TRATOX
C
      DO 100 J=1,NVM
      DO 100 NR=1,NRMAX
         XX(NVM*(NR-1)+J)=XV(J,NR)
  100 CONTINUE
      DO 110 J=1,NFM
      DO 110 NR=1,NRMAX
         YY(J,NR)=YV(J,NR)
  110 CONTINUE
C
 2000 CONTINUE
C
      CALL TRMTRX
C
      CALL BANDRD (AX,X,NVM*NRMAX,MWM,MWM,IERR)
      IF(IERR.EQ.30000) THEN
         WRITE(6,*) 'XX ERROR IN TRLOOP : MATRIX AA IS SINGULAR ',
     &              ' AT ',NT,' STEP.'
         GOTO 9000
      ENDIF
      DO 210 J=1,NFM
      DO 210 NR=1,NRMAX
         Y(J,NR) = Y(J,NR)/AY(J,NR)
  210 CONTINUE
C
      DO 220 I=1,NVM*NRMAX
         IF (ABS(X(I)-XX(I)).GT.EPSLTR*ABS(X(I))) GOTO 3000
  220 CONTINUE
      DO 230 J=1,NFM
      DO 230 NR=1,NRMAX
         IF (ABS(Y(J,NR)-YY(J,NR)).GT.EPSLTR*ABS(Y(J,NR))) GOTO 3000
  230 CONTINUE
      GOTO 4000
C
 3000 L=L+1
      IF(L.GE.LMAXTR) GOTO 4000
C
      DO 310 I=1,NVM*NRMAX
         XX(I) = X(I)
  310 CONTINUE
      DO 320 J=1,NFM
      DO 320 NR=1,NRMAX
         YY(J,NR) = Y(J,NR)
  320 CONTINUE
C
      DO 330 NR=1,NRMAX
         RN(NR,2) = 0.5D0*(XV(2,NR)+X(NVM*NR-7))
         RN(NR,3) = 0.5D0*(XV(3,NR)+X(NVM*NR-6))
         RN(NR,4) = 0.5D0*(XV(4,NR)+X(NVM*NR-5))
         RN(NR,1) = PZ(2)*RN(NR,2) +PZ(3)*RN(NR,3)+PZ(4)*RN(NR,4)
     &             +PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
         RT(NR,1)  = 0.5D0*(XV(5,NR)+X(NVM*NR-4))/(1.5D0*RN(NR,1))
         RT(NR,2)  = 0.5D0*(XV(6,NR)+X(NVM*NR-3))/(1.5D0*RN(NR,2))
         RT(NR,3)  = 0.5D0*(XV(7,NR)+X(NVM*NR-2))/(1.5D0*RN(NR,3))
         IF(RN(NR,4).LT.1.D-70) THEN
            RT(NR,4)  = 0.D0
         ELSE
            RT(NR,4)  = 0.5D0*(XV(8,NR)+X(NVM*NR-1))/(1.5D0*RN(NR,4))
         ENDIF
         BP(NR)    = 0.5D0*(XV(9,NR)+X(NVM*NR))
         QP(NR)    =  FKAP*RG(NR)*BB/(RR*BP(NR))
         RW(NR,1)  = 0.5D0*(YV(1,NR)+Y(1,NR))
         RW(NR,2)  = 0.5D0*(YV(2,NR)+Y(2,NR))
  330 CONTINUE
C
      Q0  = (4.D0*QP(1) -QP(2) )/3.D0
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
C
      DO 400 J=1,NVM
      DO 400 NR=1,NRMAX
         XV(J,NR) = X(NVM*(NR-1)+J)
  400 CONTINUE
      DO 410 NF=1,NFM
      DO 410 NR=1,NRMAX
         YV(NF,NR) = Y(NF,NR)
  410 CONTINUE
C
      CALL TRXTOA
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
      SUBROUTINE TRMTRX
C
      INCLUDE 'trcomm.h'
C
      COMMON /TRLCL1/ A(NVM,NVM,NRM),B(NVM,NVM,NRM),C(NVM,NVM,NRM)
      COMMON /TRLCL2/ D(NVM,NRM)
      DIMENSION RTM(NSM),AMZ(NSM)
C
      RKAPX=(RKAP-1.D0)/(RKAP+1.D0)
      FKAP=0.5D0*(RKAP+1.D0)
     &     *(1.D0+RKAPX/4.D0+RKAPX*RKAPX/64.D0)
      DO 5 NS=1,NSM
        AMZ(NS)=PA(NS)*AMM/PZ(NS)**2
    5 CONTINUE
C
      IF(MDLCD.EQ.0) THEN
         BPS= AMYU0*RIP*1.D6/(2.D0*PI*RA*FKAP)
      ELSE
         BPA=XV(9,NRMAX)
         RLP=RA*LOG(8.D0*RR/RA-2.D0)
         BPS=BPA-DT*EZOH(NRMAX)/RLP
      ENDIF
C
      COEF = AEE**4*5.D0*1.D20/(SQRT(2.D0*PI)*PI*AEPS0**2)
      FVL  = FKAP/RKAP
      DR1  = FVL/DR
      DR2  = FVL/(DR*DR)
C
      DO 10 NR=1,NRMAX
      DO 10 NW=1,NVM
      DO 10 NV=1,NVM
         A(NV,NW,NR)=0.D0
         B(NV,NW,NR)=0.D0
         C(NV,NW,NR)=0.D0
   10 CONTINUE
      DO 15 NR=1,NRMAX
      DO 15 NV=1,NVM
         D(NV,NR)=0.D0
   15 CONTINUE
         
C
      DO 20 NR=2,NRMAX-1
         NN=2*NSM+1
         A(NN,NN,NR) = ETA(NR  )*RG(NR-1)/(RM(NR  )*AMYU0)  *DR2
         B(NN,NN,NR) =-ETA(NR  )*RG(NR  )/(RM(NR  )*AMYU0)  *DR2
     &                -ETA(NR+1)*RG(NR  )/(RM(NR+1)*AMYU0)  *DR2
         C(NN,NN,NR) = ETA(NR+1)*RG(NR+1)/(RM(NR+1)*AMYU0)  *DR2
         D(NN   ,NR) = ETA(NR  )*(AJ(NR  )-AJOH(NR  ))      /DR
     &                -ETA(NR+1)*(AJ(NR+1)-AJOH(NR+1))      /DR
C
         DO 18 NS=1,NSM
            RTM(NS)=RT(NR,NS)*RKEV/(PA(NS)*AMM)
   18    CONTINUE
      DO 20 NS=1,NSM
         VNM=RG(NR-1)*AV(NR-1,NS)/RM(NR)
         VNP=RG(NR  )*AV(NR  ,NS)/RM(NR)
         DNM=RG(NR-1)*AD(NR-1,NS)/RM(NR)
         DNP=RG(NR  )*AD(NR  ,NS)/RM(NR)
C
         VTM=RG(NR-1)*(AVK(NR-1,NS)/1.5D0+AV(NR-1,NS))/RM(NR)
         VTP=RG(NR  )*(AVK(NR  ,NS)/1.5D0+AV(NR  ,NS))/RM(NR)
         DTM=RG(NR-1)* AK(NR-1,NS)/(1.5D0*RM(NR))
         DTP=RG(NR  )* AK(NR  ,NS)/(1.5D0*RM(NR))
C
         VXM=0.D0
         VXP=0.D0
         RTL=0.5D0*(RT(NR-1,NS)+RT(NR  ,NS))
         DXM=RG(NR-1)*(1.5D0*AD(NR-1,NS)-AK(NR-1,NS))*RTL/RM(NR)
         RTL=0.5D0*(RT(NR  ,NS)+RT(NR+1,NS))
         DXP=RG(NR  )*(1.5D0*AD(NR  ,NS)-AK(NR  ,NS))*RTL/RM(NR)
C
         A(    NS,    NS,NR) = 0.5D0*VNM*DR1+DNM*DR2
         B(    NS,    NS,NR) = 0.5D0*VNM*DR1-DNM*DR2
     &                        -0.5D0*VNP*DR1-DNP*DR2
         C(    NS,    NS,NR) =-0.5D0*VNP*DR1+DNP*DR2
C
         A(NSM+NS,NSM+NS,NR) = 0.5D0*VTM*DR1+DTM*DR2
         B(NSM+NS,NSM+NS,NR) = 0.5D0*VTM*DR1-DTM*DR2
     &                        -0.5D0*VTP*DR1-DTP*DR2
         C(NSM+NS,NSM+NS,NR) =-0.5D0*VTP*DR1+DTP*DR2
C
         A(NSM+NS,    NS,NR) = 0.5D0*VXM*DR1+DXM*DR2
         B(NSM+NS,    NS,NR) = 0.5D0*VXM*DR1-DXM*DR2
     &                        -0.5D0*VXP*DR1-DXP*DR2
         C(NSM+NS,    NS,NR) =-0.5D0*VXP*DR1+DXP*DR2
C
         D(    NS,NR) = SSIN(NR,NS)+SPE(NR,NS)/DT
         D(NSM+NS,NR) = PIN(NR,NS)/(RKEV*1.D20)
C
      DO 20 NS1=1,NSM
         IF(NS1.NE.NS) THEN
            C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1))
            B(NSM+NS,NSM+NS, NR)=B(NSM+NS,NSM+NS, NR)-RN(NR,NS1)*C1
            B(NSM+NS,NSM+NS1,NR)=B(NSM+NS,NSM+NS1,NR)+RN(NR,NS )*C1
         ENDIF
   20 CONTINUE
C
      NR=1
         NN=2*NSM+1
         A(NN,NN,NR) = 0.D0
         B(NN,NN,NR) =-ETA(NR  )*RG(NR  )/(RM(NR  )*AMYU0)  *DR2
     &                -ETA(NR+1)*RG(NR  )/(RM(NR+1)*AMYU0)  *DR2
         C(NN,NN,NR) = ETA(NR+1)*RG(NR+1)/(RM(NR+1)*AMYU0)  *DR2
         D(NN   ,NR) = ETA(NR  )*(AJ(NR  )-AJOH(NR  ))      /DR
     &                -ETA(NR+1)*(AJ(NR+1)-AJOH(NR+1))      /DR
C
         DO 28 NS=1,NSM
            RTM(NS)=RT(NR,NS)*RKEV/(PA(NS)*AMM)
   28    CONTINUE
      DO 30 NS=1,NSM
         VNM=0.D0
         VNP=RG(NR  )*AV(NR  ,NS)/RM(NR)
         DNM=0.D0
         DNP=RG(NR  )*AD(NR  ,NS)/RM(NR)
         VTM=0.D0
         VTP=RG(NR  )*(AVK(NR  ,NS)/1.5D0+AV(NR  ,NS))/RM(NR)
         DTM=0.D0
         DTP=RG(NR  )* AK(NR  ,NS)/(1.5D0*RM(NR))
         VXM=0.D0
         VXP=0.D0
         DXM=0.D0
         RTL=0.5D0*(RT(NR  ,NS)+RT(NR+1,NS))
         DXP=RG(NR  )*(1.5D0*AD(NR  ,NS)-AK(NR  ,NS))*RTL/RM(NR)
C
         A(    NS,    NS,NR) = 0.5D0*VNM*DR1+DNM*DR2
         B(    NS,    NS,NR) = 0.5D0*VNM*DR1-DNM*DR2
     &                        -0.5D0*VNP*DR1-DNP*DR2
         C(    NS,    NS,NR) =-0.5D0*VNP*DR1+DNP*DR2
C
         A(NSM+NS,NSM+NS,NR) = 0.5D0*VTM*DR1+DTM*DR2
         B(NSM+NS,NSM+NS,NR) = 0.5D0*VTM*DR1-DTM*DR2
     &                        -0.5D0*VTP*DR1-DTP*DR2
         C(NSM+NS,NSM+NS,NR) =-0.5D0*VTP*DR1+DTP*DR2
C
         A(NSM+NS,    NS,NR) = 0.5D0*VXM*DR1+DXM*DR2
         B(NSM+NS,    NS,NR) = 0.5D0*VXM*DR1-DXM*DR2
     &                        -0.5D0*VXP*DR1-DXP*DR2
         C(NSM+NS,    NS,NR) =-0.5D0*VXP*DR1+DXP*DR2
C
         D(    NS,NR) = SSIN(NR,NS)+SPE(NR,NS)/DT
         D(NSM+NS,NR) = PIN(NR,NS)/(RKEV*1.D20)
C
      DO 30 NS1=1,NSM
         IF(NS1.NE.NS) THEN
            C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1))
            B(NSM+NS,NSM+NS, NR)=B(NSM+NS,NSM+NS, NR)-RN(NR,NS1)*C1
            B(NSM+NS,NSM+NS1,NR)=B(NSM+NS,NSM+NS1,NR)+RN(NR,NS )*C1
         ENDIF
   30 CONTINUE
C
      NR=NRMAX
         NN=2*NSM+1
         A(NN,NN,NR) = 0.D0
         B(NN,NN,NR) = 0.D0
         C(NN,NN,NR) = 0.D0
         D(NN   ,NR) = 0.D0
C
         DO 38 NS=1,NSM
            RTM(NS)=RT(NR,NS)*RKEV/(PA(NS)*AMM)
   38    CONTINUE
      DO 40 NS=1,NSM
         VNM=RG(NR-1)*AV(NR-1,NS)/RM(NR)
         VNP=RG(NR  )*AV(NR  ,NS)/RM(NR)
         DNM=RG(NR-1)*AD(NR-1,NS)/RM(NR)
         DNP=RG(NR  )*AD(NR  ,NS)/RM(NR)
C
         VTM=RG(NR-1)*(AVK(NR-1,NS)/1.5D0+AV(NR-1,NS))/RM(NR)
         VTP=RG(NR  )*(AVK(NR  ,NS)/1.5D0+AV(NR  ,NS))/RM(NR)
         DTM=RG(NR-1)* AK(NR-1,NS)/(1.5D0*RM(NR))
         DTP=RG(NR  )* AK(NR  ,NS)/(1.5D0*RM(NR))
C
         VXM=0.D0
         VXP=0.D0
         RTL=0.5D0*(RT(NR-1,NS)+RT(NR,NS))
         DXM=RG(NR-1)*(1.5D0*AD(NR-1,NS)-AK(NR-1,NS))*RTL/RM(NR)
         RTL=PTS(NS)
         DXP=RG(NR  )*(1.5D0*AD(NR  ,NS)-AK(NR  ,NS))*RTL/RM(NR)
C
         A(    NS,    NS,NR) = 0.5D0*VNM*DR1+DNM*DR2
         B(    NS,    NS,NR) = 0.5D0*VNM*DR1-DNM*DR2
     &                        -0.0D0*VNP*DR1-DNP*DR2*2.D0
         C(    NS,    NS,NR) = 0.D0
C
         A(NSM+NS,NSM+NS,NR) = 0.5D0*VTM*DR1+DTM*DR2
         B(NSM+NS,NSM+NS,NR) = 0.5D0*VTM*DR1-DTM*DR2
     &                        -0.0D0*VTP*DR1-DTP*DR2*2.D0
         C(NSM+NS,NSM+NS,NR) = 0.0D0
C
         A(NSM+NS,    NS,NR) = 0.5D0*VXM*DR1+DXM*DR2
         B(NSM+NS,    NS,NR) = 0.5D0*VXM*DR1-DXM*DR2
     &                        -0.0D0*VXP*DR1-DXP*DR2*2.D0
         C(NSM+NS,    NS,NR) = 0.0D0
C
         D(    NS,NR) = SSIN(NR,NS)+SPE(NR,NS)/DT
     &                 +(-VNP*DR1+DNP*DR2*2.D0)*PNSS(NS)
         D(NSM+NS,NR) = PIN(NR,NS)/(RKEV*1.D20)
     &                 +(-VXP*DR1+DXP*DR2*2.D0)*PNSS(NS)
     &                 +(-VTP*DR1+DTP*DR2*2.D0)*PNSS(NS)*PTS(NS)*1.5D0

      DO 40 NS1=1,NSM
         IF(NS1.NE.NS) THEN
            C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1))
            B(NSM+NS,NSM+NS, NR)=B(NSM+NS,NSM+NS, NR)-RN(NR,NS1)*C1
            B(NSM+NS,NSM+NS1,NR)=B(NSM+NS,NSM+NS1,NR)+RN(NR,NS )*C1
         ENDIF
   40 CONTINUE
C
C      FADV=0.5D0
      FADV=1.0D0
C
      PRV=(1.D0-FADV)*DT
      ADV=FADV*DT
C
      DO 80 NR=1,NRMAX
      DO 80 NV=1,NVM
         X(NVM*(NR-1)+NV) = XV(NV,NR)+DT*D(NV,NR)
   80 CONTINUE
C
      NR=1
      DO 90 NW=1,NVM
      DO 90 NV=1,NVM
         X(NVM*(NR-1)+NV)=X(NVM*(NR-1)+NV)+PRV*(
     &                                         +B(NV,NW,NR)*XV(NW,NR  )
     &                                         +C(NV,NW,NR)*XV(NW,NR+1))
   90 CONTINUE
C
      DO 92 NR=2,NRMAX-1
      DO 92 NW=1,NVM
      DO 92 NV=1,NVM
         X(NVM*(NR-1)+NV)=X(NVM*(NR-1)+NV)+PRV*(A(NV,NW,NR)*XV(NW,NR-1)
     &                                         +B(NV,NW,NR)*XV(NW,NR  )
     &                                         +C(NV,NW,NR)*XV(NW,NR+1))
   92 CONTINUE
C
      NR=NRMAX
      DO 94 NW=1,NVM
      DO 94 NV=1,NVM
         X(NVM*(NR-1)+NV)=X(NVM*(NR-1)+NV)+PRV*(A(NV,NW,NR)*XV(NW,NR-1)
     &                                         +B(NV,NW,NR)*XV(NW,NR  )
     &                                         )
   94 CONTINUE
C
      DO 100 MV=1,NVM*NRMAX
      DO 100 MW=1,MWM
         AX(MW,MV) = 0.D0
  100 CONTINUE
C
      DO 110 NV=1,NVM
      DO 110 NW=1,NVM
      DO 110 NR=1,NRMAX
         AX(  NVM+NW-NV,NVM*(NR-1)+NV) = -ADV*A(NV,NW,NR)
         AX(2*NVM+NW-NV,NVM*(NR-1)+NV) = -ADV*B(NV,NW,NR)
         AX(3*NVM+NW-NV,NVM*(NR-1)+NV) = -ADV*C(NV,NW,NR)
  110 CONTINUE
C
      DO 120 NV=1,NVM
      DO 120 NR=1,NRMAX
         AX(2*NVM,NVM*(NR-1)+NV) = AX(2*NVM,NVM*(NR-1)+NV) + 1.D0
  120 CONTINUE
C
C     ***** Surface Boundary Condition for Bp *****
C
      MV=NVM*NRMAX
      DO 130 MW=1,MWM
         AX(MW,MV)=0.D0
  130 CONTINUE
      AX(2*NVM,MV)=1.D0
      X(MV)=BPS
C
C     ***** Evolution of fast ion components *****
C
      DO 150 NR=1,NRMAX
         Y(1,NR)=(1.D0-PRV/TAUB(NR))*YV(1,NR)+PNB(NR)*DT/(RKEV*1.D20)
         Y(2,NR)=(1.D0-PRV/TAUF(NR))*YV(2,NR)+PNF(NR)*DT/(RKEV*1.D20)
         AY(1,NR)=1.D0+ADV/TAUB(NR)
         AY(2,NR)=1.D0+ADV/TAUF(NR)
  150 CONTINUE
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
      DO 100 NR=1,NRMAX
         XV(1,NR) = RN(NR,1)
         XV(2,NR) = RN(NR,2)
         XV(3,NR) = RN(NR,3)
         XV(4,NR) = RN(NR,4)
         XV(5,NR) = 1.5D0*RN(NR,1)*RT(NR,1)
         XV(6,NR) = 1.5D0*RN(NR,2)*RT(NR,2)
         XV(7,NR) = 1.5D0*RN(NR,3)*RT(NR,3)
         XV(8,NR) = 1.5D0*RN(NR,4)*RT(NR,4)
         XV(9,NR) = BP(NR)
         YV(1,NR) = 1.5D0*RW(NR,1)
         YV(2,NR) = 1.5D0*RW(NR,2)
  100 CONTINUE
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
      DO 100 NR=1,NRMAX
         RN(NR,2) = XV(2,NR)
         RN(NR,3) = XV(3,NR)
         RN(NR,4) = XV(4,NR)
         RN(NR,1) = PZ(2)*RN(NR,2)+PZ(3)*RN(NR,3)+PZ(4)*RN(NR,4)
     &             +PZC(NR)*ANC(NR)+PZFE(NR)*ANFE(NR)
         RT(NR,1) = XV(5,NR)/(1.5D0*RN(NR,1))
         RT(NR,2) = XV(6,NR)/(1.5D0*RN(NR,2))
         RT(NR,3) = XV(7,NR)/(1.5D0*RN(NR,3))
         IF(RN(NR,4).LT.1.D-70)THEN
            RT(NR,4) = 0.D0
         ELSE
            RT(NR,4) = XV(8,NR)/(1.5D0*RN(NR,4))
         ENDIF
         BP(NR)  = XV(9,NR)
         QP(NR)  = FKAP*RG(NR)*BB/(RR*BP(NR))
         RW(NR,1)  = YV(1,NR)/1.5D0
         RW(NR,2)  = YV(2,NR)/1.5D0
  100 CONTINUE
C
      Q0  = (4.D0*QP(1) -QP(2) )/3.D0
C
      RETURN
      END
C
C     ***********************************************************
C
C           SAVE TRANSIENT DATA FOR GRAPHICS
C
C     ******************************************************
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
      DO 10 NS=1,NSM
      DO 10 NR=1,NRMAX
         IF(RT(NR,NS).LT.0.D0) THEN
            GOTO 100
         ENDIF
   10 CONTINUE
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
