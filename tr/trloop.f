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
         RW(NR,1)  = 0.5D0*(YV(1,NR)+Y(1,NR))
         RW(NR,2)  = 0.5D0*(YV(2,NR)+Y(2,NR))
  330 CONTINUE
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
      DO NS=1,NSM
         AMZ(NS)=PA(NS)*AMM/PZ(NS)**2
      ENDDO
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
C
      NN=2*NSM+1
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
            FC4=FC3
         ELSE
            FC4=DVRHO(NR+2)*ABRHO(NR+2)/TTRHO(NR+2)
         ENDIF
         F1M=FAM*AR1RHO(NR  )/DR
         F1P=FAP*AR1RHO(NR+1)/DR
         F2M=FAM*FBM/(DR*DR)
         F2P=FAP*FBP/(DR*DR)
         FCM=0.5D0*(FC1+FC2)
         FC0=0.5D0*(FC2+FC3)
         FCP=0.5D0*(FC3+FC4)
C
         A(NN,NN,NR) = F2M*FCM
         B(NN,NN,NR) =-F2M*FC0
     &                -F2P*FC0
         C(NN,NN,NR) = F2P*FCP
         D(NN   ,NR) = F1M*(AJ(NR  )-AJOH(NR  ))
     &                -F1P*(AJ(NR+1)-AJOH(NR+1))
C
         FA1=DVRHO(NR-1)*AR1RHO(NR-1)/DR
         FA2=DVRHO(NR  )*AR1RHO(NR  )/DR
         FA3=DVRHO(NR+1)*AR1RHO(NR+1)/DR
         FB1=DVRHO(NR-1)*AR2RHO(NR-1)/(DR*DR)
         FB2=DVRHO(NR  )*AR2RHO(NR  )/(DR*DR)
         FB3=DVRHO(NR+1)*AR2RHO(NR+1)/(DR*DR)
         DO NS=1,NSM
            RTM(NS)=RT(NR,NS)*RKEV/(PA(NS)*AMM)
         ENDDO
C
      DO NS=1,NSM
         VNM=0.5D0*(FA1*AV(NR-1,NS)+FA2*AV(NR  ,NS))
         VNP=0.5D0*(FA2*AV(NR  ,NS)+FA3*AV(NR+1,NS))
         DNM=0.5D0*(FB1*AD(NR-1,NS)+FB2*AD(NR  ,NS))
         DNP=0.5D0*(FB2*AD(NR  ,NS)+FB3*AD(NR+1,NS))
C
         DV23=DVRHO(NR)**(2.D0/3.D0)
         DV53=DVRHO(NR)**(5.D0/3.D0)
         AVK1=AVK(NR-1,NS)+AV(NR-1,NS)*2.5D0
         AVK2=AVK(NR  ,NS)+AV(NR  ,NS)*2.5D0
         AVK3=AVK(NR+1,NS)+AV(NR+1,NS)*2.5D0
         ADK1=AK(NR-1,NS)
         ADK2=AK(NR  ,NS)
         ADK3=AK(NR+1,NS)
         VTM=0.5D0*(FA1*AVK1+FA2*AVK2)*DV23
         VTP=0.5D0*(FA2*AVK2+FA3*AVK3)*DV23
         DTM=0.5D0*(FB1*ADK1+FB2*ADK2)*DV23
         DTP=0.5D0*(FB2*ADK2+FB3*ADK3)*DV23
C
         ADX1=(2.5D0*AD(NR-1,NS)-AK(NR-1,NS))*RT(NR-1,NS)
         ADX2=(2.5D0*AD(NR  ,NS)-AK(NR  ,NS))*RT(NR  ,NS)
         ADX3=(2.5D0*AD(NR+1,NS)-AK(NR+1,NS))*RT(NR+1,NS)
         VXM=0.D0*DV23
         VXP=0.D0*DV23
         DXM=0.5D0*(FB1*ADX1+FB2*ADX2)*DV23
         DXP=0.5D0*(FB2*ADX2+FB3*ADX3)*DV23
C
         A(    NS,    NS,NR) = 0.5D0*VNM+DNM
         B(    NS,    NS,NR) = 0.5D0*VNM-DNM
     &                        -0.5D0*VNP-DNP
         C(    NS,    NS,NR) =-0.5D0*VNP+DNP
C
         A(NSM+NS,NSM+NS,NR) = 0.5D0*VTM+DTM
         B(NSM+NS,NSM+NS,NR) = 0.5D0*VTM-DTM
     &                        -0.5D0*VTP-DTP
         C(NSM+NS,NSM+NS,NR) =-0.5D0*VTP+DTP
C
         A(NSM+NS,    NS,NR) = 0.5D0*VXM+DXM
         B(NSM+NS,    NS,NR) = 0.5D0*VXM-DXM
     &                        -0.5D0*VXP-DXP
         C(NSM+NS,    NS,NR) =-0.5D0*VXP+DXP
C
         D(    NS,NR) = (SSIN(NR,NS)+SPE(NR,NS)/DT)*DV53
         D(NSM+NS,NR) = (PIN(NR,NS)/(RKEV*1.D20)  )*DV53
C
      DO NS1=1,NSM
         IF(NS1.NE.NS) THEN
            C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1))*DV53
            B(NSM+NS,NSM+NS, NR)=B(NSM+NS,NSM+NS, NR)-RN(NR,NS1)*C1
            B(NSM+NS,NSM+NS1,NR)=B(NSM+NS,NSM+NS1,NR)+RN(NR,NS )*C1
         ENDIF
      ENDDO
      ENDDO
      ENDDO
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
         F1M=FAM*AR1RHO(NR  )/DR
         F1P=FAP*AR1RHO(NR+1)/DR
         F2M=FAM*FBM/(DR*DR)
         F2P=FAP*FBP/(DR*DR)
         FCM=0.5D0*(FC1+FC2)
         FC0=0.5D0*(FC2+FC3)
         FCP=0.5D0*(FC3+FC4)
C
         A(NN,NN,NR) = F2M*FCM
         B(NN,NN,NR) =-F2M*FC0
     &                -F2P*FC0
         C(NN,NN,NR) = F2P*FCP
         D(NN   ,NR) = F1M*(AJ(NR  )-AJOH(NR  ))
     &                -F1P*(AJ(NR+1)-AJOH(NR+1))
C
         FA1=0.D0
         FA2=DVRHO(NR  )*AR1RHO(NR  )/DR
         FA3=DVRHO(NR+1)*AR1RHO(NR+1)/DR
         FB1=0.D0
         FB2=DVRHO(NR  )*AR2RHO(NR  )/(DR*DR)
         FB3=DVRHO(NR+1)*AR2RHO(NR+1)/(DR*DR)
         DO NS=1,NSM
            RTM(NS)=RT(NR,NS)*RKEV/(PA(NS)*AMM)
         ENDDO
C
      DO NS=1,NSM
         VNM=0.D0
         VNP=0.5D0*(FA2*AV(NR  ,NS)+FA3*AV(NR+1,NS))
         DNM=0.D0
         DNP=0.5D0*(FB2*AD(NR  ,NS)+FB3*AD(NR+1,NS))
C
         DV23=DVRHO(NR)**(2.D0/3.D0)
         DV53=DVRHO(NR)**(5.D0/3.D0)
         AVK1=0.D0
         AVK2=AVK(NR  ,NS)+AV(NR  ,NS)*2.5D0
         AVK3=AVK(NR+1,NS)+AV(NR+1,NS)*2.5D0
         ADK1=0.D0
         ADK2=AK(NR  ,NS)
         ADK3=AK(NR+1,NS)
         VTM=0.D0
         VTP=0.5D0*(FA2*AVK2+FA3*AVK3)*DV23
         DTM=0.D0
         DTP=0.5D0*(FB2*ADK2+FB3*ADK3)*DV23
C
         ADX2=(2.5D0*AD(NR  ,NS)-AK(NR  ,NS))*RT(NR  ,NS)
         ADX3=(2.5D0*AD(NR+1,NS)-AK(NR+1,NS))*RT(NR+1,NS)
         VXM=0.D0
         VXP=0.D0*DV23
         DXM=0.0D0
         DXP=0.5D0*(FB2*ADX2+FB3*ADX3)*DV23
C
         A(    NS,    NS,NR) = 0.5D0*VNM+DNM
         B(    NS,    NS,NR) = 0.5D0*VNM-DNM
     &                        -0.5D0*VNP-DNP
         C(    NS,    NS,NR) =-0.5D0*VNP+DNP
C
         A(NSM+NS,NSM+NS,NR) = 0.5D0*VTM+DTM
         B(NSM+NS,NSM+NS,NR) = 0.5D0*VTM-DTM
     &                        -0.5D0*VTP-DTP
         C(NSM+NS,NSM+NS,NR) =-0.5D0*VTP+DTP
C
         A(NSM+NS,    NS,NR) = 0.5D0*VXM+DXM
         B(NSM+NS,    NS,NR) = 0.5D0*VXM-DXM
     &                        -0.5D0*VXP-DXP
         C(NSM+NS,    NS,NR) =-0.5D0*VXP+DXP
C
         D(    NS,NR) = (SSIN(NR,NS)+SPE(NR,NS)/DT)*DV53
         D(NSM+NS,NR) = (PIN(NR,NS)/(RKEV*1.D20)  )*DV53
C
      DO NS1=1,NSM
         IF(NS1.NE.NS) THEN
            C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1))*DV53
            B(NSM+NS,NSM+NS, NR)=B(NSM+NS,NSM+NS, NR)-RN(NR,NS1)*C1
            B(NSM+NS,NSM+NS1,NR)=B(NSM+NS,NSM+NS1,NR)+RN(NR,NS )*C1
         ENDIF
      ENDDO
      ENDDO
C
      NR=NRMAX
         A(NN,NN,NR) = 0.D0
         B(NN,NN,NR) = 0.D0
         C(NN,NN,NR) = 0.D0
         D(NN   ,NR) = 0.D0
C
         FA1=DVRHO(NR-1)*AR1RHO(NR-1)/DR
         FA2=DVRHO(NR  )*AR1RHO(NR  )/DR
         FA3=DVRHO(NR+1)*AR1RHO(NR+1)/DR
         FB1=DVRHO(NR-1)*AR2RHO(NR-1)/(DR*DR)
         FB2=DVRHO(NR  )*AR2RHO(NR  )/(DR*DR)
         FB3=DVRHO(NR+1)*AR2RHO(NR+1)/(DR*DR)
         DO NS=1,NSM
            RTM(NS)=RT(NR,NS)*RKEV/(PA(NS)*AMM)
         ENDDO
C
      DO NS=1,NSM
         VNM=0.5D0*(FA1*AV(NR-1,NS)+FA2*AV(NR  ,NS))
         VNP=       FA2*AV(NR  ,NS)
         DNM=0.5D0*(FB1*AD(NR-1,NS)+FB2*AD(NR  ,NS))
         DNP=       FB2*AD(NR  ,NS)
C
         DV23=DVRHO(NR)**(2.D0/3.D0)
         DV53=DVRHO(NR)**(5.D0/3.D0)
         AVK1=AVK(NR-1,NS)+AV(NR-1,NS)*2.5D0
         AVK2=AVK(NR  ,NS)+AV(NR  ,NS)*2.5D0
         ADK1=AK(NR-1,NS)
         ADK2=AK(NR  ,NS)
         VTM=0.5D0*(FA1*AVK1+FA2*AVK2)*DV23
         VTP=       FA2*AVK2          *DV23
         DTM=0.5D0*(FB1*ADK1+FB2*ADK2)*DV23
         DTP=       FB2*ADK2          *DV23
C
         ADX1=(2.5D0*AD(NR-1,NS)-AK(NR-1,NS))*RT(NR-1,NS)
         ADX2=(2.5D0*AD(NR  ,NS)-AK(NR  ,NS))*RT(NR  ,NS)
         VXM=0.D0*DV23
         VXP=0.D0*DV23
         DXM=0.5D0*(FB1*ADX1+FB2*ADX2)*DV23
         DXP=       FB2*ADX2          *DV23
C
         A(    NS,    NS,NR) = 0.5D0*VNM+DNM
         B(    NS,    NS,NR) = 0.5D0*VNM-DNM
     &                        -0.0D0*VNP-DNP-DNP
         C(    NS,    NS,NR) =-0.0D0
C
         A(NSM+NS,NSM+NS,NR) = 0.5D0*VTM+DTM
         B(NSM+NS,NSM+NS,NR) = 0.5D0*VTM-DTM
     &                        -0.0D0*VTP-DTP-DTP
         C(NSM+NS,NSM+NS,NR) =-0.0D0
C
         A(NSM+NS,    NS,NR) = 0.5D0*VXM+DXM
         B(NSM+NS,    NS,NR) = 0.5D0*VXM-DXM
     &                        -0.0D0*VXP-DXP-DXP
         C(NSM+NS,    NS,NR) =-0.0D0
C
         D(    NS,NR) = (SSIN(NR,NS)+SPE(NR,NS)/DT)*DV53
     &                 +(-VNP+DNP*2.D0)*PNSS(NS)
         D(NSM+NS,NR) = (PIN(NR,NS)/(RKEV*1.D20)  )*DV53
     &                 +(-VXP+DXP*2.D0)*PNSS(NS)
     &                 +(-VTP+DTP*2.D0)*PNSS(NS)*PTS(NS)*1.5D0
C
      DO NS1=1,NSM
         IF(NS1.NE.NS) THEN
            C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1))*DV53
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
      DO NV=1,NVM
         X(NVM*(NR-1)+NV) = XV(NV,NR)+DT*D(NV,NR)
      ENDDO
      ENDDO
C
      NR=1
      DO NW=1,NVM
      DO NV=1,NVM
         X(NVM*(NR-1)+NV)=X(NVM*(NR-1)+NV)+PRV*(
     &                                         +B(NV,NW,NR)*XV(NW,NR  )
     &                                         +C(NV,NW,NR)*XV(NW,NR+1))
      ENDDO
      ENDDO
C
      DO NR=2,NRMAX-1
      DO NW=1,NVM
      DO NV=1,NVM
         X(NVM*(NR-1)+NV)=X(NVM*(NR-1)+NV)+PRV*(A(NV,NW,NR)*XV(NW,NR-1)
     &                                         +B(NV,NW,NR)*XV(NW,NR  )
     &                                         +C(NV,NW,NR)*XV(NW,NR+1))
      ENDDO
      ENDDO
      ENDDO
C
      NR=NRMAX
      DO NW=1,NVM
      DO NV=1,NVM
         X(NVM*(NR-1)+NV)=X(NVM*(NR-1)+NV)+PRV*(A(NV,NW,NR)*XV(NW,NR-1)
     &                                         +B(NV,NW,NR)*XV(NW,NR  )
     &                                         )
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
      DO NV=1,NVM
      DO NW=1,NVM
         AX(  NVM+NW-NV,NVM*(NR-1)+NV) = -ADV*A(NV,NW,NR)
         AX(2*NVM+NW-NV,NVM*(NR-1)+NV) = -ADV*B(NV,NW,NR)
         AX(3*NVM+NW-NV,NVM*(NR-1)+NV) = -ADV*C(NV,NW,NR)
      ENDDO
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         DV53=DVRHO(NR)**(5.D0/3.D0)
         DO NS=1,NSM
            AX(2*NVM,NVM*(NR-1)+NS) 
     &           = AX(2*NVM,NVM*(NR-1)+NS) + 1.D0
            AX(2*NVM,NVM*(NR-1)+NSM+NS) 
     &           = AX(2*NVM,NVM*(NR-1)+NSM+NS) + 1.5D0*DV53
         ENDDO
            AX(2*NVM,NVM*(NR-1)+NN) 
     &           = AX(2*NVM,NVM*(NR-1)+NN) + 1.D0
      ENDDO
C
C     ***** Surface Boundary Condition for Bp *****
C
      MV=NVM*NRMAX
      DO MW=1,MWM
         AX(MW,MV)=0.D0
      ENDDO
      AX(2*NVM,MV)=1.D0
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
      SUBROUTINE TRATOX
C
      INCLUDE 'trcomm.h'
C
      DO NR=1,NRMAX
         XV(1,NR) = RN(NR,1)
         XV(2,NR) = RN(NR,2)
         XV(3,NR) = RN(NR,3)
         XV(4,NR) = RN(NR,4)
         XV(5,NR) = RN(NR,1)*RT(NR,1)
         XV(6,NR) = RN(NR,2)*RT(NR,2)
         XV(7,NR) = RN(NR,3)*RT(NR,3)
         XV(8,NR) = RN(NR,4)*RT(NR,4)
         XV(9,NR) = BP(NR)
         YV(1,NR) = RW(NR,1)
         YV(2,NR) = RW(NR,2)
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
         RN(NR,2) = XV(2,NR)
         RN(NR,3) = XV(3,NR)
         RN(NR,4) = XV(4,NR)
         RN(NR,1) = PZ(2)*RN(NR,2)+PZ(3)*RN(NR,3)+PZ(4)*RN(NR,4)
     &             +PZC(NR)*ANC(NR)+PZFE(NR)*ANFE(NR)
         RT(NR,1) = XV(5,NR)/RN(NR,1)
         RT(NR,2) = XV(6,NR)/RN(NR,2)
         RT(NR,3) = XV(7,NR)/RN(NR,3)
         IF(RN(NR,4).LT.1.D-70)THEN
            RT(NR,4) = 0.D0
         ELSE
            RT(NR,4) = XV(8,NR)/RN(NR,4)
         ENDIF
         BP(NR)  = XV(9,NR)
         RW(NR,1)  = YV(1,NR)
         RW(NR,2)  = YV(2,NR)
      ENDDO
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
