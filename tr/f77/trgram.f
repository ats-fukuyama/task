C
C     ***********************************************************
C
C           GRAPHIC : CONTROL ROUTINE FOR COMPARISON
C
C     ***********************************************************
C
      SUBROUTINE TRCOMP(K2,INQ)
C
      INCLUDE 'trcomm.inc'
      CHARACTER K2*1
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
C
      IF(K2.EQ.'1') CALL TRCMP1(INQ)
      IF(K2.EQ.'2') CALL TRCMP2(INQ)
      IF(K2.EQ.'3') CALL TRCMP3(INQ)
      IF(K2.EQ.'4') CALL TRCMP4(INQ)
      IF(K2.EQ.'5') CALL TRCMP5(INQ)
      IF(K2.EQ.'6') CALL TRCMP6(INQ)
      IF(K2.EQ.'7') CALL TRCMP7(INQ)
C
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
C
      RETURN
      END
C
C     ***********************************
C
C        COMPARE WITH DIFFERENT MODELS
C
C     ***********************************
C
      SUBROUTINE TRCMP1(INQ)
C
      INCLUDE 'trcomm.inc'
      DIMENSION RGFLSUM(NRMP,NSM),RQFLSUM(NRMP,NSM)
      DIMENSION RNN(NRM,NSM),DNN(NRM,NSM),DTN(NRM,NSM)
      DIMENSION AKNCL(NRM,NSM),ADNCL(NRM,NSM)
      DIMENSION AJBSSTCK(NRM),ETASTCK(NRM)
C
      DO NR=1,NRMAX
         ETASTCK(NR)=ETA(NR)
         AJBSSTCK(NR)=AJBS(NR)
      ENDDO
C
C     *** Bootstrap Current and Neoclassical Resistivity ***
C
      CALL TRAJBS
      DO NR=1,NRMAX
         GJB(NR,1)=GUCLIP(AJBS(NR)*1.D-6)
      ENDDO
      CALL TRAJBSNEW
      DO NR=1,NRMAX
         GJB(NR,2)=GUCLIP(AJBS(NR)*1.D-6)
      ENDDO
      CALL TRAJBSSAUTER
      DO NR=1,NRMAX
         GJB(NR,3)=GUCLIP(AJBS(NR)*1.D-6)
      ENDDO
C
      MDLETASTCK=MDLETA
      MDNCLSSTCK=MDNCLS
      IF(MDNCLS.EQ.1) MDNCLS=0
      DO MDLETA=1,4
         CALL TRCFET
         DO NR=1,NRMAX
            GET(NR,MDLETA)=GLOG(ETA(NR),1.D-10,1.D0)
         ENDDO
      ENDDO
      MDLETA=MDLETASTCK
      MDNCLS=MDNCLSSTCK
C
      IF(MDNCLS.EQ.0) MDNCLS=1
      CALL TR_NCLASS(IERR)
      CALL TRAJBS_NCLASS
      DO NR=1,NRMAX
         GJB(NR,4)=GUCLIP(AJBS(NR)*1.D-6)
         GET(NR,5)=GLOG(ETANC(NR),1.D-10,1.D0)
      ENDDO
      MDNCLS=MDNCLSSTCK
      DO NR=1,NRMAX
         ETA(NR)=ETASTCK(NR)
         AJBS(NR)=AJBSSTCK(NR)
      ENDDO
C
C     *** Neoclassical Particle and Heat Flux Diffusivity ***
C
      MODE=0
      IF(MODE.EQ.0) THEN
         DO NR=1,NRMAX
            DO NS=1,2
               AKNCL(NR,NS)=AKNCP(NR,NS,NS)+AKNCT(NR,NS,NS)
               ADNCL(NR,NS)=ADNC(NR,NS)
            ENDDO
         ENDDO
      ELSE
         DO NR=1,NRMAX
            DO NS=1,NSMAX
               RGFLSUM(NR,NS)=0.D0
               RQFLSUM(NR,NS)=0.D0
               DO NA=1,5
                  RGFLSUM(NR,NS)=RGFLSUM(NR,NS)+RGFLS(NR,NA,NS)
                  RQFLSUM(NR,NS)=RQFLSUM(NR,NS)+RQFLS(NR,NA,NS)
               ENDDO
            ENDDO
         ENDDO
         DO NR=1,NRMAX-1
            DO NS=1,NSMAX
               RNN(NR,NS)=(RN(NR+1,NS)+RN(NR,NS))*0.5D0
               DNN(NR,NS)=(RN(NR+1,NS)-RN(NR,NS))     *RJCB(NR)/DR
               DTN(NR,NS)=(RT(NR+1,NS)-RT(NR,NS))*RKEV*RJCB(NR)/DR
            ENDDO
         ENDDO
         NR=NRMAX
         DO NS=1,NSMAX
            RNN(NR,NS)=PNSS(NS)
            DNN(NR,NS)=DERIV3P(PNSS(NS),RN(NR,NS),RN(NR-1,NS),
     &                        RHOG(NR),RHOM(NR),RHOM(NR-1))
            DTN(NR,NS)=DERIV3P(PTS(NS),RT(NR,NS),RT(NR-1,NS),
     &                        RHOG(NR),RHOM(NR),RHOM(NR-1))*RKEV
         ENDDO
         DO NR=1,NRMAX
            DO NS=1,NSMAX
               ADNCL(NR,NS)=-RGFLSUM(NR,NS)/DNN(NR,NS)
               AKNCL(NR,NS)=-RQFLSUM(NR,NS)/(RNN(NR,NS)*DTN(NR,NS))
            ENDDO
         ENDDO
      ENDIF
C
      MDLKNCSTCK=MDLKNC
      DO MDLKNC=1,3,2
         CALL TRCFNC
         DO NR=1,NRMAX
            GAK(NR+1,MDLKNC+2) = GUCLIP(AKNC(NR,1))
            GAK(NR+1,MDLKNC+3) = GUCLIP(AKNC(NR,2))
         ENDDO
         GAK(1,MDLKNC+2) = GUCLIP(AKNC(1,1))
         GAK(1,MDLKNC+3) = GUCLIP(AKNC(1,2))
      ENDDO
      MDLKNC=MDLKNCSTCK
      DO NR=1,NRMAX
         GAD(NR+1,1) = GUCLIP(ADNCL(NR,1))
         GAD(NR+1,2) = GUCLIP(ADNCL(NR,2))
         GAK(NR+1,1) = GUCLIP(AKNCL(NR,1))
         GAK(NR+1,2) = GUCLIP(AKNCL(NR,2))
         GAD(NR+1,3) = GUCLIP(ADNC(NR,1))
         GAD(NR+1,4) = GUCLIP(ADNC(NR,2))
      ENDDO
         GAD(1,1) = GUCLIP(ADNCL(1,1))
         GAD(1,2) = GUCLIP(ADNCL(1,2))
         GAK(1,1) = GUCLIP(AKNCL(1,1))
         GAK(1,2) = GUCLIP(AKNCL(1,2))
         GAD(1,3) = GUCLIP(ADNC(1,1))
         GAD(1,4) = GUCLIP(ADNC(1,2))
C
C     *** Graphic Routine ***
C
      CALL PAGES
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GJB,NRMP,NRMAX,4,
     &            '@JBS [MA/m$+2$=]  vs r@',2+INQ)
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GET,NRMP,NRMAX,5,
     &            '@LOG:ETA  vs r @',11+INQ)
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GAD,NRMP,NRMAX+1,4,
     &            '@ADNCE, ADNCD [m$+2$=/s]  vs r@',2+INQ)
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GAK,NRMP,NRMAX+1,6,
     &            '@AKNCE, AKNCD [m$+2$=/s]  vs r @',2+INQ)
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
C
C     ****************************************
C
C        COMPARE WITH UFILE DATA (TEMPORAL)
C
C     ****************************************
C
      SUBROUTINE TRCMP2(INQ)
C
      INCLUDE 'trcomm.inc'
      DIMENSION TMU(NTUM)
      DIMENSION TE0(NTUM),TI0(NTUM),WTOT(NTUM),RIBS(NTUM),RIPL(NTUM)
      DIMENSION PICRH(NTUM),PNBI(NTUM)
      CHARACTER KFID*10
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
C
      IF(MDLUF.EQ.0) RETURN
C
      IF(MDLUF.EQ.2) THEN
         AMP=1.D-3
         KFID='TE0'
         CALL UF1DT(KFID,KUFDEV,KUFDCG,NTS,ATE0 ,AMP,MDLXP,IERR)
         KFID='TI0'
         CALL UF1DT(KFID,KUFDEV,KUFDCG,NTS,ATI0 ,AMP,MDLXP,IERR)
         AMP=1.D-6
         KFID='WTOT'
         CALL UF1DT(KFID,KUFDEV,KUFDCG,NTS,AWTOT,AMP,MDLXP,IERR)
         KFID='IBOOT'
         CALL UF1DT(KFID,KUFDEV,KUFDCG,NTS,ARIBS,AMP,MDLXP,IERR)
         KFID='IP'
         CALL UF1DT(KFID,KUFDEV,KUFDCG,NTS,ARIPL,AMP,MDLXP,IERR)
         KFID='PICRH'
         CALL UF1DT(KFID,KUFDEV,KUFDCG,NTS,AICRH,AMP,MDLXP,IERR)
         KFID='PNBI'
         CALL UF1DT(KFID,KUFDEV,KUFDCG,NTS,ANBI ,AMP,MDLXP,IERR)
C
         CALL PAGES
         CALL SETLIN(0,0,7)
         CALL SETCHS(0.3,0.0)
C         CALL SETFNT(32)
         CALL GTEXT( 2.0, 17.0, '< CALCULATED >', 14, 0)
         CALL GTEXT(12.0, 17.0, '< EXPERIMENT >', 14, 0)
C
         CALL MOVE(12.0, 16.0)
         CALL TEXT('Te0 =', 5)
         CALL NUMBD(ATE0,'(1F7.3)',7)
         CALL TEXT('[keV]', 5)
C
         CALL MOVE(12.0, 15.0)
         CALL TEXT('Ti0 =', 5)
         CALL NUMBD(ATI0,'(1F7.3)',7)
         CALL TEXT('[keV]', 5)
C
         CALL MOVE(12.0, 14.0)
         CALL TEXT('Wtot=', 5)
         CALL NUMBD(AWTOT,'(1F7.3)',7)
         CALL TEXT('[MJ]', 4)
C
         CALL MOVE(12.0, 13.0)
         CALL TEXT('IBS =', 5)
         CALL NUMBD(ARIBS,'(1F7.3)',7)
         CALL TEXT('[MA]', 4)
C
         CALL MOVE(12.0, 12.0)
         CALL TEXT('IP  =', 5)
         CALL NUMBD(ARIPL,'(1F7.3)',7)
         CALL TEXT('[MA]', 4)
C
         CALL MOVE(12.0, 11.0)
         CALL TEXT('PICH=', 5)
         CALL NUMBD(AICRH,'(1F7.3)',7)
         CALL TEXT('[MW]', 4)
C
         CALL MOVE(12.0, 10.0)
         CALL TEXT('PNBI=', 5)
         CALL NUMBD(ANBI,'(1F7.3)',7)
         CALL TEXT('[MW]', 4)
C
C     ++++++++++++++++++++++++++++++++++++++
C
         CALL MOVE( 2.0, 16.0)
         CALL TEXT('Te0 =', 5)
         CALL NUMBD(TS0(1),'(1F7.3)',7)
         CALL TEXT('[keV]', 5)
C
         CALL MOVE( 2.0, 15.0)
         CALL TEXT('Ti0 =', 5)
         CALL NUMBD(TS0(2),'(1F7.3)',7)
         CALL TEXT('[keV]', 5)
C
         CALL MOVE( 2.0, 14.0)
         CALL TEXT('Wtot=', 5)
         CALL NUMBD(WPT,'(1F7.3)',7)
         CALL TEXT('[MJ]', 4)
C
         CALL MOVE( 2.0, 13.0)
         CALL TEXT('IBS =', 5)
         CALL NUMBD(AJBST,'(1F7.3)',7)
         CALL TEXT('[MA]', 4)
C
         CALL MOVE( 2.0, 12.0)
         CALL TEXT('IP  =', 5)
         CALL NUMBD(AJT,'(1F7.3)',7)
         CALL TEXT('[MA]', 4)
C
         CALL MOVE( 2.0, 11.0)
         CALL TEXT('PICH=', 5)
         CALL NUMBD(PRFST,'(1F7.3)',7)
         CALL TEXT('[MW]', 4)
C
         CALL MOVE( 2.0, 10.0)
         CALL TEXT('PNBI=', 5)
         CALL NUMBD(PEXST,'(1F7.3)',7)
         CALL TEXT('[MW]', 4)
C
         CALL TRGRTM
         CALL PAGEE
      ELSE
      ICK=2
      TMUMAX=0.D0
C
      AMP=1.D-3
      KFID='TE0'
      CALL UF1DG(KFID,KUFDEV,KUFDCG,GT,TMU,TE0  ,AMP,NGT,MDLXP,IERR)
      KFID='TI0'
      CALL UF1DG(KFID,KUFDEV,KUFDCG,GT,TMU,TI0  ,AMP,NGT,MDLXP,IERR)
      AMP=1.D-6
      KFID='WTOT'
      CALL UF1DG(KFID,KUFDEV,KUFDCG,GT,TMU,WTOT ,AMP,NGT,MDLXP,IERR)
      KFID='IBOOT'
      CALL UF1DG(KFID,KUFDEV,KUFDCG,GT,TMU,RIBS ,AMP,NGT,MDLXP,IERR)
      KFID='IP'
      CALL UF1DG(KFID,KUFDEV,KUFDCG,GT,TMU,RIPL ,AMP,NGT,MDLXP,IERR)
      IF(MDLUF.EQ.3) THEN
         KFID='PLH'
         CALL UF1DG(KFID,KUFDEV,KUFDCG,GT,TMU,PICRH,AMP,NGT,MDLXP,IERR)
      ELSE
         KFID='PICRH'
         CALL UF1DG(KFID,KUFDEV,KUFDCG,GT,TMU,PICRH,AMP,NGT,MDLXP,IERR)
      ENDIF
      KFID='PNBI'
      CALL UF1DG(KFID,KUFDEV,KUFDCG,GT,TMU,PNBI ,AMP,NGT,MDLXP,IERR)
C
      CALL PAGES
C
      DO NG=1,NGT
         GYT(NG,1)=GVT(NG, 9)
         GYT(NG,2)=GUCLIP(TE0(NG))
         TSL=DBLE(GT(NG))
         NR=NRMAX/2
         CALL TIMESPL(TSL,PTE,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
         GYT(NG,3)=G3D(NR,NG,1)
         GYT(NG,4)=GUCLIP(PTE)
      ENDDO
      CALL TRGR1D( 3.0,12.0,14.0,17.0,GT,GYT,NTM,NGT,4,
     &            '@TE0(TR),TE0(XP) [keV]  vs t@',2+INQ)
C
      DO NG=1,NGT
         GYT(NG,1)=GVT(NG,10)
         GYT(NG,2)=GUCLIP(TI0(NG))
         TSL=DBLE(GT(NG))
         NR=NRMAX/2
         CALL TIMESPL(TSL,PTI,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
         GYT(NG,3)=G3D(NR,NG,2)
         GYT(NG,4)=GUCLIP(PTI)
      ENDDO
      CALL TRGR1D(15.0,24.0,14.0,17.0,GT,GYT,NTM,NGT,4,
     &            '@TI0(TR),TI0(XP) [keV]  vs t@',2+INQ)
C
      DO NG=1,NGT
         GYT(NG,1)=GVT(NG,33)
         GYT(NG,2)=GUCLIP(WTOT(NG))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 9.7,12.7,GT,GYT,NTM,NGT,2,
     &            '@WTOT(TR),WTOT(XP) [MJ]  vs t@',2+INQ)
C
      DO NG=1,NGT
         GYT(NG,1)=GVT(NG,42)+GVT(NG,43)
         GYT(NG,2)=GUCLIP(PICRH(NG))
      ENDDO
      CALL TRGR1D(15.0,24.0, 9.7,12.7,GT,GYT,NTM,NGT,2,
     &            '@PICRH(TR),PICRH(XP) [MW]  vs t@',2+INQ)
C
      DO NG=1,NGT
         GYT(NG,1)=GVT(NG,38)
         GYT(NG,2)=GUCLIP(RIBS(NG))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 5.4, 8.4,GT,GYT,NTM,NGT,2,
     &            '@IBS(TR),IBS(XP) [MA]  vs t@',2+INQ)
C
      DO NG=1,NGT
         GYT(NG,1)=GVT(NG,89)
         GYT(NG,2)=GVT(NG,90)
         GYT(NG,3)=GVT(NG,89)+GVT(NG,90)
         GYT(NG,4)=GUCLIP(PNBI(NG))
      ENDDO
      CALL TRGR1D(15.0,24.0, 5.4, 8.4,GT,GYT,NTM,NGT,4,
     &            '@PNBIE;I;TOT(TR),PNBI(XP) [MW]  vs t@',2+INQ)
C
      DO NG=1,NGT
         GYT(NG,1)=GVT(NG,34)
         GYT(NG,2)=GUCLIP(RIPL(NG))
      ENDDO
      CALL TRGR1D( 3.0,12.0, 1.1, 4.1,GT,GYT,NTM,NGT,2,
     &            '@IP(TR),IP(XP) [MA]  vs t@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
      ENDIF
C
      RETURN
      END
C
C     **********************************************
C
C        COMPARE WITH UFILE DATA (RADIAL PROFILE)
C
C     **********************************************
C
      SUBROUTINE TRCMP3(INQ)
C
      INCLUDE 'trcomm.inc'
      COMMON /TMSLC1/ TMU(NTUM),TMU1(NTUM)
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
C
      TSL=DT*DBLE(NT)
      IF(MDLUF.EQ.1) THEN
         CALL PAGES
C
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(RT(NR,1))
            CALL TIMESPL(TSL,RTL,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
            GYR(NR,2) = GUCLIP(RTL)
         ENDDO
         CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &               '@TE(TR),TE(XP) [keV]  vs r@',2+INQ)
C
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(RT(NR,2))
            CALL TIMESPL(TSL,RTL,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
            GYR(NR,2) = GUCLIP(RTL)
          ENDDO
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &               '@TI(TR),TI(XP) [keV]  vs r@',2+INQ)
C
         IF(MDLJQ.NE.1) THEN
            DO NR=1,NRMAX
               GYR(NR+1,1) = GUCLIP(QP(NR))
               CALL TIMESPL(TSL,QPL,TMU,QPU(1,NR),NTXMAX,NTUM,IERR)
               GYR(NR+1,2) = GUCLIP(QPL)
            ENDDO
            GYR(1,1) = GUCLIP((4.D0*QP(1)-QP(2))/3.D0)
            GYR(1,2) = GUCLIP((4.0*GYR(2,2)-GYR(3,2))/3.0)
            CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,
     &                  '@QP(TR),QP(XP) vs r@',2+INQ)
         ELSE
            DO NR=1,NRMAX
               GYR(NR,1) = GUCLIP(AJ(NR)    *1.D-6)
               CALL TIMESPL(TSL,AJL,TMU,AJU(1,NR),NTXMAX,NTUM,IERR)
               GYR(NR,2) = GUCLIP(AJL*1.D-6)
            ENDDO
            CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,
     &                  '@AJ(TR),AJ(XP) [MA/m$+2$=]  vs r@',2+INQ)
         ENDIF
C
         DO NR=1,NRMAX
            GYR(NR+1,1) = GUCLIP(BP(NR))
            CALL TIMESPL(TSL,BPL,TMU,BPU(1,NR),NTXMAX,NTUM,IERR)
            GYR(NR+1,2) = GUCLIP(BPL)
         ENDDO
         GYR(1,1) = 0.0
         GYR(1,2) = 0.0
         CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,
     &               '@BP(TR),BP(XP) [T] vs r@',2+INQ)
C
         CALL TRGRTM
         CALL PAGEE
      ELSEIF(MDLUF.EQ.2) THEN
         CALL PAGES
C
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(RT(NR,1))
            GYR(NR,2) = GUCLIP(RTU(1,NR,1))
         ENDDO
         CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &               '@TE(TR),TE(XP) [keV]  vs r@',2+INQ)
C
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(RT(NR,2))
            GYR(NR,2) = GUCLIP(RTU(1,NR,2))
         ENDDO
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &               '@TI(TR),TI(XP) [keV]  vs r@',2+INQ)
C
         IF(MDLJQ.NE.1) THEN
            DO NR=1,NRMAX
               GYR(NR+1,1) = GUCLIP(QP(NR))
               GYR(NR+1,2) = GUCLIP(QPU(1,NR))
            ENDDO
            GYR(1,1) = GUCLIP((4.D0*QP(1)-QP(2))/3.D0)
            GYR(1,2) = GUCLIP((4.0*GYR(2,2)-GYR(3,2))/3.0)
            CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,
     &                  '@QP(TR),QP(XP)  vs r@',2+INQ)
         ELSE
            DO NR=1,NRMAX
               GYR(NR,1) = GUCLIP(AJ(NR)   *1.D-6)
               GYR(NR,2) = GUCLIP(AJU(1,NR)*1.D-6)
            ENDDO
            CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,
     &                  '@AJ(TR),AJ(XP) [MA/m$+2$=]  vs r@',2+INQ)
         ENDIF
C
         DO NR=1,NRMAX
            GYR(NR+1,1) = GUCLIP(BP(NR))
            GYR(NR+1,2) = GUCLIP(BPU(1,NR))
         ENDDO
         GYR(1,1) = 0.0
         GYR(1,2) = 0.0
         CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,
     &               '@BP(TR),BP(XP) [T] vs r@',2+INQ)
C
         CALL TRGRTM
         CALL PAGEE
      ELSEIF(MDLUF.EQ.3) THEN
         CALL PAGES
C
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(RT(NR,1))
            CALL TIMESPL(TSL,RTL,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
            GYR(NR,2) = GUCLIP(RTL)
         ENDDO
         CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &               '@TE(TR),TE(XP) [keV]  vs r@',2+INQ)
C
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(RT(NR,2))
            CALL TIMESPL(TSL,RTL,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
            GYR(NR,2) = GUCLIP(RTL)
         ENDDO
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &               '@TI(TR),TI(XP) [keV]  vs r@',2+INQ)
C
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(AJ(NR)    *1.D-6)
            CALL TIMESPL(TSL,AJL,TMU,AJU(1,NR),NTXMAX,NTUM,IERR)
            GYR(NR,2) = GUCLIP(AJL*1.D-6)
         ENDDO
         CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,
     &               '@AJ(TR),AJ(XP) [MA/m$+2$=]  vs r@',2+INQ)
C
         IF(KUFDEV.EQ.'X'.AND.KUFDCG.EQ.'14') THEN
            DO NR=1,NRMAX
               GYR(NR,1) = GUCLIP(AJBS(NR)    *1.D-6)
               CALL TIMESPL(TSL,AJL,TMU,AJBSU(1,NR),NTXMAX,NTUM,IERR)
               GYR(NR,2) = GUCLIP(AJL*1.D-6)
            ENDDO
            CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,
     &                  '@AJBS(TR),AJBS(XP) [MA/m$+2$=]  vs r@',2+INQ)
         ELSE
            DO NR=1,NRMAX
               GYR(NR+1,1) = GUCLIP(QP(NR))
               CALL TIMESPL(TSL,QPL,TMU,QPU(1,NR),NTXMAX,NTUM,IERR)
               GYR(NR+1,2) = GUCLIP(QPL)
            ENDDO
            GYR(1,1) = GUCLIP((4.D0*QP(1)-QP(2))/3.D0)
            GYR(1,2) = GUCLIP((4.0*GYR(2,2)-GYR(3,2))/3.0)
            CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,
     &                  '@QP(TR),QP(XP)  vs r@',2+INQ)
         ENDIF
         CALL TRGRTM
         CALL PAGEE
      ELSE
         RETURN
      ENDIF
C
      RETURN
      END
C
C     **********************************************
C
C        COMPARE WITH UFILE DATA (ERROR BAR)
C
C     **********************************************
C
      SUBROUTINE TRCMP4(INQ)
C
      INCLUDE 'trcomm.inc'
      COMMON /TRUFC4/ NREMAX(2),GRE(NRM,2)
      COMMON /TRERU1/ RTEXU(NTUM,NRMU),   RTIXU(NTUM,NRMU),
     &                RNEXU(NTUM,NRMU)
      COMMON /TRERU2/ RTEXEU(NTUM,NRMU),  RTIXEU(NTUM,NRMU),
     &                RNEXEU(NTUM,NRMU)
C
      IF(MDLUF.EQ.2) THEN
         CALL PAGES
         DO NG=1,2
            DO NR=1,NRMP
               GYR(NR,NG)=0.0
            ENDDO
         ENDDO
         DO NR=1,NREMAX(1)
            GYR(NR,1) = GUCLIP(RTEXU(1,NR))
            GER(NR,1) = GUCLIP(RTEXEU(1,NR))
         ENDDO
         DO NR=1,NRMAX
            GYR(NR,2) = GUCLIP(RT(NR,1))
         ENDDO
         CALL TRGR1DE( 3.0,12.0,11.0,17.0,GRM,GYR,GRE(1,1),GER,NRMP,
     &                NRMAX,NREMAX(1),2,
     &                '@TE(XP),TE(TR) [keV]  vs r@',2+INQ)
C     
         DO NR=1,NREMAX(1)
            GYR(NR,1) = GUCLIP(RTIXU(1,NR))
            GER(NR,1) = GUCLIP(RTIXEU(1,NR))
         ENDDO
         DO NR=1,NRMAX
            GYR(NR,2) = GUCLIP(RT(NR,2))
         ENDDO
         CALL TRGR1DE(15.5,24.5,11.0,17.0,GRM,GYR,GRE(1,1),GER,NRMP,
     &                NRMAX,NREMAX(1),2,
     &                '@TI(XP),TI(TR) [keV]  vs r@',2+INQ)
C     
         DO NR=1,NREMAX(2)
            GYR(NR,1) = GUCLIP(RNEXU(1,NR))
            GER(NR,1) = GUCLIP(RNEXEU(1,NR))
         ENDDO
         DO NR=1,NRMAX
            GYR(NR,2) = GUCLIP(RN(NR,1))
         ENDDO            
         CALL TRGR1DE( 3.0,12.0, 2.0, 8.0,GRM,GYR,GRE(1,2),GER,NRMP,
     &                NRMAX,NREMAX(2),2,
     &                '@NE(XP),NE(TR) [10$+20$=/m$+3$=]  vs r@',2+INQ)
         CALL TRGRTM
         CALL PAGEE
      ENDIF
C
      RETURN
      END
C
C     **********************************************
C
C        COMPARE WITH UFILE DATA (OTHER)
C
C     **********************************************
C
      SUBROUTINE TRCMP5(INQ)
C
      INCLUDE 'trcomm.inc'
      COMMON /TMSLC1/ TMU(NTUM),TMU1(NTUM)
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
C
      IF(MDLUF.EQ.0) RETURN
C
      CALL PAGES
      TSL=DT*DBLE(NT)
      IF(MDLUF.EQ.2) THEN
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(POH(NR)*1.D-6)
            GYR(NR,2) = GUCLIP(POHU(1,NR)*1.D-6)
         ENDDO
      ELSE
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(POH(NR)*1.D-6)
            CALL TIMESPL(TSL,RTL,TMU,POHU(1,NR),NTXMAX,NTUM,IERR)
            GYR(NR,2) = GUCLIP(RTL*1.D-6)
         ENDDO
      ENDIF
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &            '@POH(TR),POH(XP) [MW/m$+3$=]  vs r@',2+INQ)
C
      IF(MDLUF.EQ.2) THEN
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(ZEFF(NR))
            GYR(NR,2) = GUCLIP(ZEFFU_ORG(1,NR))
         ENDDO
      ELSE
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(ZEFF(NR))
            CALL TIMESPL(TSL,RTL,TMU,ZEFFU_ORG(1,NR),NTXMAX,NTUM,IERR)
            GYR(NR,2) = GUCLIP(RTL)
         ENDDO
      ENDIF
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,
     &            '@ZEFF(TR),ZEFF(XP) [keV]  vs r@',2+INQ)
C
      IF(MDLUF.EQ.2) THEN
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(RN(NR,2))
            GYR(NR,2) = GUCLIP(RNU_ORG(1,NR,2)+RNU_ORG(1,NR,4))
         ENDDO
      ELSE
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(RN(NR,2))
            CALL TIMESPL(TSL,RN2L,TMU,RNU_ORG(1,NR,2),NTXMAX,NTUM,IERR)
            CALL TIMESPL(TSL,RN4L,TMU,RNU_ORG(1,NR,4),NTXMAX,NTUM,IERR)
            GYR(NR,2) = GUCLIP(RN2L+RN4L)
         ENDDO
      ENDIF
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,
     &     '@NI(TR), NI(XP) [10$+20$=/m$+3$=]  vs r@',2+INQ)
C
      IF(MDLUF.EQ.2) THEN
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(RN(NR,3))
            GYR(NR,2) = GUCLIP(RNU_ORG(1,NR,3))
         ENDDO
      ELSE
         DO NR=1,NRMAX
            GYR(NR,1) = GUCLIP(RN(NR,3))
            CALL TIMESPL(TSL,RNL,TMU,RNU_ORG(1,NR,3),NTXMAX,NTUM,IERR)
            GYR(NR,2) = GUCLIP(RNL)
         ENDDO
      ENDIF
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,
     &     '@Nimp(TR), Nimp(XP) [10$+20$=/m$+3$=]  vs r@',2+INQ)
C
      CALL PAGEE
      CALL TRGRTM
C
      RETURN
      END
C
C     ****************************************************
C
C        COMPARE WITH UFILE DATA (TEMPORAL; TEMPERATURE)
C
C     ****************************************************
C
      SUBROUTINE TRCMP6(INQ)
C
      INCLUDE 'trcomm.inc'
      CHARACTER KFID*10
      COMMON /TMSLC1/ TMU(NTUM),TMU1(NTUM)
      COMMON /TMSLC3/ NTXMAX,NTXMAX1
C
      IF(MDLUF.EQ.0.OR.MDLUF.EQ.2) RETURN
C
      CALL PAGES
C
      DO NG=1,NGT
         GYT(NG,1)=GVT(NG, 9)
         TSL=DBLE(GT(NG))
         NR=1
         CALL TIMESPL(TSL,PTE,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
         GYT(NG,2)=GUCLIP(PTE)
         NR=INT(0.3D0*NRMAX)
         CALL TIMESPL(TSL,PTE,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
         GYT(NG,3)=G3D(NR,NG,1)
         GYT(NG,4)=GUCLIP(PTE)
         NR=INT(0.6D0*NRMAX)
         CALL TIMESPL(TSL,PTE,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
         GYT(NG,5)=G3D(NR,NG,1)
         GYT(NG,6)=GUCLIP(PTE)
      ENDDO
      CALL TRGR1D( 3.0,24.0, 9.7,17.0,GT,GYT,NTM,NGT,6,
     &            '@TE0(TR),TE0(XP) [keV]  vs t@',2+INQ)
C
      DO NG=1,NGT
         GYT(NG,1)=GVT(NG,10)
         TSL=DBLE(GT(NG))
         NR=1
         CALL TIMESPL(TSL,PTI,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
         GYT(NG,2)=GUCLIP(PTI)
         NR=INT(0.3D0*NRMAX)
         CALL TIMESPL(TSL,PTI,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
         GYT(NG,3)=G3D(NR,NG,2)
         GYT(NG,4)=GUCLIP(PTI)
         NR=INT(0.6D0*NRMAX)
         CALL TIMESPL(TSL,PTI,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
         GYT(NG,5)=G3D(NR,NG,2)
         GYT(NG,6)=GUCLIP(PTI)
      ENDDO
      CALL TRGR1D( 3.0,24.0, 1.1, 8.4,GT,GYT,NTM,NGT,6,
     &            '@TI0(TR),TI0(XP) [keV]  vs t@',2+INQ)
C
      CALL TRGRTM
      CALL PAGEE
C
      RETURN
      END
C
C     ****************************************************
C
C        COMPARE WITH RADIAL ELECTRIC FIELD MODELS
C
C     ****************************************************
C
      SUBROUTINE TRCMP7(INQ)
C
      INCLUDE 'trcomm.inc'
      DIMENSION ER_ORG(NRM)
C
      CALL PAGES
C
      MDLER_ORG=MDLER
      DO NR=1,NRMAX
         ER_ORG(NR)=ER(NR)
      ENDDO
      DO MDLER=1,4
         CALL TRERAD
         DO NR=1,NRMAX
            GYR(NR,MDLER)=GUCLIP(ER(NR))
         ENDDO
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX,4,
     &            '@ER  vs r@',2+INQ)
C
      MDLER=MDLER_ORG
      DO NR=1,NRMAX
         ER(NR)=ER_ORG(NR)
      ENDDO
C
      CALL PAGEE
      CALL TRGRTM
C
      RETURN
      END
