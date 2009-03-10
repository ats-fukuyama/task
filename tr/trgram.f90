!     ***********************************************************

!           GRAPHIC : CONTROL ROUTINE FOR COMPARISON

!     ***********************************************************

      SUBROUTINE TRCOMP(K2,INQ)

      USE TRCOMM, ONLY : NRAMAX, NRMAX, NROMAX, RHOA
      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN):: K2
      INTEGER(4),       INTENT(IN):: INQ

      IF(RHOA.NE.1.D0) NRMAX=NROMAX

      IF(K2.EQ.'1') CALL TRCMP1(INQ)
      IF(K2.EQ.'2') CALL TRCMP2(INQ)
      IF(K2.EQ.'3') CALL TRCMP3(INQ)
      IF(K2.EQ.'4') CALL TRCMP4(INQ)
      IF(K2.EQ.'5') CALL TRCMP5(INQ)
      IF(K2.EQ.'6') CALL TRCMP6(INQ)
      IF(K2.EQ.'7') CALL TRCMP7(INQ)

      IF(RHOA.NE.1.D0) NRMAX=NRAMAX

      RETURN
      END SUBROUTINE TRCOMP

!     ***********************************

!        COMPARE WITH DIFFERENT MODELS

!     ***********************************

      SUBROUTINE TRCMP1(INQ)

      USE TRCOMM, ONLY : ADNC, AJBS, AKNC, AKNCP, AKNCT, DR, ETA, ETANC, GAD, GAK, GET, GJB, GRG, GRM, MDLETA, MDLKNC, &
     &                   MDNCLS, NRM, NRMAX, NRMP, NSM, NSMAX, PNSS, PTS, RGFLS, RHOG, RHOM, RJCB, RKEV, RN, RQFLS, RT
      IMPLICIT NONE
      INTEGER(4),INTENT(IN):: INQ
      INTEGER(4)                  :: IERR, MDLETASTCK, MDLKNCSTCK, MDNCLSSTCK, MODE, NR, NS
      REAL(8), DIMENSION(NRMAX)     :: AJBSSTCK, ETASTCK
      REAL(8), DIMENSION(NRMAX,NSMAX) :: ADNCL, AKNCL, DNN, DTN, RNN
      REAL(8), DIMENSION(NRMP,NSMAX):: RGFLSUM, RQFLSUM
      REAL(8) :: DERIV3P
      REAL(4) :: GLOG


      ETASTCK(1:NRMAX)=ETA(1:NRMAX)
      AJBSSTCK(1:NRMAX)=AJBS(1:NRMAX)

!     *** Bootstrap Current and Neoclassical Resistivity ***

      CALL TRAJBS
      GJB(1:NRMAX,1)=SNGL(AJBS(1:NRMAX))*1.E-6
      CALL TRAJBSNEW
      GJB(1:NRMAX,2)=SNGL(AJBS(1:NRMAX))*1.E-6
      CALL TRAJBSSAUTER
      GJB(1:NRMAX,3)=SNGL(AJBS(1:NRMAX))*1.E-6

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

      IF(MDNCLS.EQ.0) MDNCLS=1
      CALL TR_NCLASS(IERR)
      CALL TRAJBS_NCLASS
      GJB(1:NRMAX,4)=SNGL(AJBS(1:NRMAX))*1.E-6
      DO NR=1,NRMAX
         GET(NR,5)=GLOG(ETANC(NR),1.D-10,1.D0)
      ENDDO
      MDNCLS=MDNCLSSTCK
      ETA(1:NRMAX)=ETASTCK(1:NRMAX)
      AJBS(1:NRMAX)=AJBSSTCK(1:NRMAX)

!     *** Neoclassical Particle and Heat Flux Diffusivity ***

      MODE=0
      IF(MODE.EQ.0) THEN
         DO NS=1,2
            AKNCL(1:NRMAX,NS)=AKNCP(1:NRMAX,NS,NS)+AKNCT(1:NRMAX,NS,NS)
         ENDDO
         ADNCL(1:NRMAX,1:2)=ADNC(1:NRMAX,1:2)
      ELSE
         DO NR=1,NRMAX
            DO NS=1,NSMAX
               RGFLSUM(NR,NS)=SUM(RGFLS(NR,1:5,NS))
               RQFLSUM(NR,NS)=SUM(RQFLS(NR,1:5,NS))
            ENDDO
         ENDDO
         DO NS=1,NSMAX
            RNN(1:NRMAX-1,NS)=(RN(1:NRMAX-1+1,NS)+RN(1:NRMAX-1,NS))*0.5D0
            DNN(1:NRMAX-1,NS)=(RN(1:NRMAX-1+1,NS)-RN(1:NRMAX-1,NS))     *RJCB(1:NRMAX-1)/DR
            DTN(1:NRMAX-1,NS)=(RT(1:NRMAX-1+1,NS)-RT(1:NRMAX-1,NS))*RKEV*RJCB(1:NRMAX-1)/DR
         ENDDO
         NR=NRMAX
         DO NS=1,NSMAX
            RNN(NR,NS)=PNSS(NS)
            DNN(NR,NS)=DERIV3P(PNSS(NS),RN(NR,NS),RN(NR-1,NS),RHOG(NR),RHOM(NR),RHOM(NR-1))
            DTN(NR,NS)=DERIV3P(PTS(NS),RT(NR,NS),RT(NR-1,NS), RHOG(NR),RHOM(NR),RHOM(NR-1))*RKEV
         ENDDO
         ADNCL(1:NRMAX,1:NSMAX)=-RGFLSUM(1:NRMAX,1:NSMAX)/DNN(1:NRMAX,1:NSMAX)
         AKNCL(1:NRMAX,1:NSMAX)=-RQFLSUM(1:NRMAX,1:NSMAX)/(RNN(1:NRMAX,1:NSMAX)*DTN(1:NRMAX,1:NSMAX))
      ENDIF

      MDLKNCSTCK=MDLKNC
      DO MDLKNC=1,3,2
         CALL TRCFNC
         GAK(2:NRMAX+1,MDLKNC+2) = SNGL(AKNC(1:NRMAX,1))
         GAK(2:NRMAX+1,MDLKNC+3) = SNGL(AKNC(1:NRMAX,2))
         GAK(1,MDLKNC+2) = SNGL(AKNC(1,1))
         GAK(1,MDLKNC+3) = SNGL(AKNC(1,2))
      ENDDO
      MDLKNC=MDLKNCSTCK
      GAD(2:NRMAX+1,1) = SNGL(ADNCL(1:NRMAX,1))
      GAD(2:NRMAX+1,2) = SNGL(ADNCL(1:NRMAX,2))
      GAK(2:NRMAX+1,1) = SNGL(AKNCL(1:NRMAX,1))
      GAK(2:NRMAX+1,2) = SNGL(AKNCL(1:NRMAX,2))
      GAD(2:NRMAX+1,3) = SNGL(ADNC(1:NRMAX,1))
      GAD(2:NRMAX+1,4) = SNGL(ADNC(1:NRMAX,2))
      GAD(1,1) = SNGL(ADNCL(1,1))
      GAD(1,2) = SNGL(ADNCL(1,2))
      GAK(1,1) = SNGL(AKNCL(1,1))
      GAK(1,2) = SNGL(AKNCL(1,2))
      GAD(1,3) = SNGL(ADNC(1,1))
      GAD(1,4) = SNGL(ADNC(1,2))

!     *** Graphic Routine ***

      CALL PAGES
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GJB,NRMP,NRMAX,4,'@JBS [MA/m$+2$=]  vs r@',2+INQ)
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GET,NRMP,NRMAX,5,'@LOG:ETA  vs r @',11+INQ)
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GAD,NRMP,NRMAX+1,4,'@ADNCE, ADNCD [m$+2$=/s]  vs r@',2+INQ)
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GAK,NRMP,NRMAX+1,6,'@AKNCE, AKNCD [m$+2$=/s]  vs r @',2+INQ)
      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRCMP1

!     ****************************************

!        COMPARE WITH UFILE DATA (TEMPORAL)

!     ****************************************

      SUBROUTINE TRCMP2(INQ)

      USE TRCOMM, ONLY : AJBST, AJTTOR, G3D, GT, GVT, GYT, KUFDCG, KUFDEV, MDLUF, MDLXP, NGT, NRMAX, NTM, NTS, &
     &                   NTUM, PEXST, PRFST, RTU, TS0, WPT
      USE TRCOM1, ONLY : NTXMAX
      IMPLICIT NONE
      INTEGER(4), INTENT(IN):: INQ
      REAL(8)    :: AICRH, AMP, ANBI, ARIBS, ARIPL, ATE0, ATI0, AWTOT, PTE, PTI, TMUMAX, TSL
      REAL(8),DIMENSION(NTUM) :: PICRH, PNBI, RIBS, RIPL, TE0, TI0, TMU, WTOT
      INTEGER(4) :: ICK, IERR, NG, NR
      CHARACTER(LEN=10) ::  KFID


      IF(MDLUF.EQ.0) RETURN

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

         CALL PAGES
         CALL SETLIN(0,0,7)
         CALL SETCHS(0.3,0.0)
!         CALL SETFNT(32)
         CALL GTEXT( 2.0, 17.0, '< CALCULATED >', 14, 0)
         CALL GTEXT(12.0, 17.0, '< EXPERIMENT >', 14, 0)

         CALL MOVE(12.0, 16.0)
         CALL TEXT('Te0 =', 5)
         CALL NUMBD(ATE0,'(1F7.3)',7)
         CALL TEXT('[keV]', 5)

         CALL MOVE(12.0, 15.0)
         CALL TEXT('Ti0 =', 5)
         CALL NUMBD(ATI0,'(1F7.3)',7)
         CALL TEXT('[keV]', 5)

         CALL MOVE(12.0, 14.0)
         CALL TEXT('Wtot=', 5)
         CALL NUMBD(AWTOT,'(1F7.3)',7)
         CALL TEXT('[MJ]', 4)

         CALL MOVE(12.0, 13.0)
         CALL TEXT('IBS =', 5)
         CALL NUMBD(ARIBS,'(1F7.3)',7)
         CALL TEXT('[MA]', 4)

         CALL MOVE(12.0, 12.0)
         CALL TEXT('IP  =', 5)
         CALL NUMBD(ARIPL,'(1F7.3)',7)
         CALL TEXT('[MA]', 4)

         CALL MOVE(12.0, 11.0)
         CALL TEXT('PICH=', 5)
         CALL NUMBD(AICRH,'(1F7.3)',7)
         CALL TEXT('[MW]', 4)

         CALL MOVE(12.0, 10.0)
         CALL TEXT('PNBI=', 5)
         CALL NUMBD(ANBI,'(1F7.3)',7)
         CALL TEXT('[MW]', 4)

!     ++++++++++++++++++++++++++++++++++++++

         CALL MOVE( 2.0, 16.0)
         CALL TEXT('Te0 =', 5)
         CALL NUMBD(TS0(1),'(1F7.3)',7)
         CALL TEXT('[keV]', 5)

         CALL MOVE( 2.0, 15.0)
         CALL TEXT('Ti0 =', 5)
         CALL NUMBD(TS0(2),'(1F7.3)',7)
         CALL TEXT('[keV]', 5)

         CALL MOVE( 2.0, 14.0)
         CALL TEXT('Wtot=', 5)
         CALL NUMBD(WPT,'(1F7.3)',7)
         CALL TEXT('[MJ]', 4)

         CALL MOVE( 2.0, 13.0)
         CALL TEXT('IBS =', 5)
         CALL NUMBD(AJBST,'(1F7.3)',7)
         CALL TEXT('[MA]', 4)

         CALL MOVE( 2.0, 12.0)
         CALL TEXT('IP  =', 5)
!!         CALL NUMBD(AJT,'(1F7.3)',7)
         CALL NUMBD(AJTTOR,'(1F7.3)',7)
         CALL TEXT('[MA]', 4)

         CALL MOVE( 2.0, 11.0)
         CALL TEXT('PICH=', 5)
         CALL NUMBD(PRFST,'(1F7.3)',7)
         CALL TEXT('[MW]', 4)

         CALL MOVE( 2.0, 10.0)
         CALL TEXT('PNBI=', 5)
         CALL NUMBD(PEXST,'(1F7.3)',7)
         CALL TEXT('[MW]', 4)

         CALL TRGRTM
         CALL PAGEE
      ELSE
      ICK=2
      TMUMAX=0.D0

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

      CALL PAGES

      GYT(1:NGT,1)=GVT(1:NGT,9)
      GYT(1:NGT,2)=SNGL(TE0(1:NGT))
      DO NG=1,NGT
         TSL=DBLE(GT(NG))
         NR=NRMAX/2
         CALL TIMESPL(TSL,PTE,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
         GYT(NG,3)=G3D(NR,NG,1)
         GYT(NG,4)=SNGL(PTE)
      ENDDO
      CALL TRGR1D( 3.0,12.0,14.0,17.0,GT,GYT,NTM,NGT,4,'@TE0(TR),TE0(XP) [keV]  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,10)
      GYT(1:NGT,2)=SNGL(TI0(1:NGT))
      DO NG=1,NGT
         TSL=DBLE(GT(NG))
         NR=NRMAX/2
         CALL TIMESPL(TSL,PTI,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
         GYT(NG,3)=G3D(NR,NG,2)
         GYT(NG,4)=SNGL(PTI)
      ENDDO
      CALL TRGR1D(15.0,24.0,14.0,17.0,GT,GYT,NTM,NGT,4,'@TI0(TR),TI0(XP) [keV]  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,33)
      GYT(1:NGT,2)=SNGL(WTOT(1:NGT))
      CALL TRGR1D( 3.0,12.0, 9.7,12.7,GT,GYT,NTM,NGT,2,'@WTOT(TR),WTOT(XP) [MJ]  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,42)+GVT(1:NGT,43)
      GYT(1:NGT,2)=SNGL(PICRH(1:NGT))
      CALL TRGR1D(15.0,24.0, 9.7,12.7,GT,GYT,NTM,NGT,2,'@PICRH(TR),PICRH(XP) [MW]  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,38)
      GYT(1:NGT,2)=SNGL(RIBS(1:NGT))
      CALL TRGR1D( 3.0,12.0, 5.4, 8.4,GT,GYT,NTM,NGT,2,'@IBS(TR),IBS(XP) [MA]  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,89)
      GYT(1:NGT,2)=GVT(1:NGT,90)
      GYT(1:NGT,3)=GVT(1:NGT,89)+GVT(1:NGT,90)
      GYT(1:NGT,4)=SNGL(PNBI(1:NGT))
      CALL TRGR1D(15.0,24.0, 5.4, 8.4,GT,GYT,NTM,NGT,4,'@PNBIE;I;TOT(TR),PNBI(XP) [MW]  vs t@',2+INQ)

      GYT(1:NGT,1)=GVT(1:NGT,34)
      GYT(1:NGT,2)=SNGL(RIPL(1:NGT))
      CALL TRGR1D( 3.0,12.0, 1.1, 4.1,GT,GYT,NTM,NGT,2,'@IP(TR),IP(XP) [MA]  vs t@',2+INQ)

      CALL TRGRTM
      CALL PAGEE
      ENDIF

      RETURN
      END SUBROUTINE TRCMP2

!     **********************************************

!        COMPARE WITH UFILE DATA (RADIAL PROFILE)

!     **********************************************

      SUBROUTINE TRCMP3(INQ)

      USE TRCOMM, ONLY : AJ, AJBS, AJBSU, AJU, BP, BPU, DT, GRG, GRM, GYR, KUFDCG, KUFDEV, MDLJQ, MDLUF, NRMAX, NRMP, &
     &                   NT, NTUM, QP, QPU, RT, RTU
      USE TRCOM1, ONLY : NTXMAX, TMU
      IMPLICIT NONE
      INTEGER(4), INTENT(IN):: INQ
      INTEGER(4):: IERR, NR
      REAL(8)   :: AJL, BPL, QPL, RTL, TSL


      TSL=DT*DBLE(NT)
      select case(MDLUF)
      case(1)
         CALL PAGES

         GYR(1:NRMAX,1) = SNGL(RT(1:NRMAX,1))
         DO NR=1,NRMAX
            CALL TIMESPL(TSL,RTL,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
            GYR(NR,2) = SNGL(RTL)
         ENDDO
         CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,'@TE(TR),TE(XP) [keV]  vs r@',2+INQ)

         GYR(1:NRMAX,1) = SNGL(RT(1:NRMAX,2))
         DO NR=1,NRMAX
            CALL TIMESPL(TSL,RTL,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
            GYR(NR,2) = SNGL(RTL)
          ENDDO
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,'@TI(TR),TI(XP) [keV]  vs r@',2+INQ)

         IF(MDLJQ.NE.1) THEN
            GYR(2:NRMAX+1,1) = SNGL(QP(1:NRMAX))
            DO NR=1,NRMAX
               CALL TIMESPL(TSL,QPL,TMU,QPU(1,NR),NTXMAX,NTUM,IERR)
               GYR(NR+1,2) = SNGL(QPL)
            ENDDO
            GYR(1,1) = SNGL((4.D0*QP(1)-QP(2))/3.D0)
            GYR(1,2) = (4.0*GYR(2,2)-GYR(3,2))/3.0
            CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,'@QP(TR),QP(XP) vs r@',2+INQ)
         ELSE
            GYR(1:NRMAX,1) = SNGL(AJ(1:NRMAX))*1.E-6
            DO NR=1,NRMAX
               CALL TIMESPL(TSL,AJL,TMU,AJU(1,NR),NTXMAX,NTUM,IERR)
               GYR(NR,2) = SNGL(AJL*1.D-6)
            ENDDO
            CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,'@AJ(TR),AJ(XP) [MA/m$+2$=]  vs r@',2+INQ)
         ENDIF

         GYR(2:NRMAX+1,1) = SNGL(BP(1:NRMAX))
         DO NR=1,NRMAX
            CALL TIMESPL(TSL,BPL,TMU,BPU(1,NR),NTXMAX,NTUM,IERR)
            GYR(NR+1,2) = SNGL(BPL)
         ENDDO
         GYR(1,1) = 0.0
         GYR(1,2) = 0.0
         CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,'@BP(TR),BP(XP) [T] vs r@',2+INQ)

         CALL TRGRTM
         CALL PAGEE
      case(2)
         CALL PAGES

         GYR(1:NRMAX,1) = SNGL(RT(1:NRMAX,1))
         GYR(1:NRMAX,2) = SNGL(RTU(1,1:NRMAX,1))
         CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,'@TE(TR),TE(XP) [keV]  vs r@',2+INQ)

         GYR(1:NRMAX,1) = SNGL(RT(1:NRMAX,2))
         GYR(1:NRMAX,2) = SNGL(RTU(1,1:NRMAX,2))
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,'@TI(TR),TI(XP) [keV]  vs r@',2+INQ)

         IF(MDLJQ.NE.1) THEN
            GYR(2:NRMAX+1,1) = SNGL(QP(1:NRMAX))
            GYR(2:NRMAX+1,2) = SNGL(QPU(1,1:NRMAX))
            GYR(1,1) = SNGL((4.D0*QP(1)-QP(2))/3.D0)
            GYR(1,2) = (4.0*GYR(2,2)-GYR(3,2))/3.0
            CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,'@QP(TR),QP(XP)  vs r@',2+INQ)
         ELSE
            GYR(1:NRMAX,1) = SNGL(AJ(1:NRMAX)   *1.D-6)
            GYR(1:NRMAX,2) = SNGL(AJU(1,1:NRMAX)*1.D-6)
            CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,'@AJ(TR),AJ(XP) [MA/m$+2$=]  vs r@',2+INQ)
         ENDIF

         GYR(2:NRMAX+1,1) = SNGL(BP(1:NRMAX))
         GYR(2:NRMAX+1,2) = SNGL(BPU(1,1:NRMAX))
         GYR(1,1) = 0.0
         GYR(1,2) = 0.0
         CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,'@BP(TR),BP(XP) [T] vs r@',2+INQ)

         CALL TRGRTM
         CALL PAGEE
      case(3)
         CALL PAGES

         GYR(1:NRMAX,1) = SNGL(RT(1:NRMAX,1))
         DO NR=1,NRMAX
            CALL TIMESPL(TSL,RTL,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
            GYR(NR,2) = SNGL(RTL)
         ENDDO
         CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,'@TE(TR),TE(XP) [keV]  vs r@',2+INQ)

         GYR(1:NRMAX,1) = SNGL(RT(1:NRMAX,2))
         DO NR=1,NRMAX
            CALL TIMESPL(TSL,RTL,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
            GYR(NR,2) = SNGL(RTL)
         ENDDO
         CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,'@TI(TR),TI(XP) [keV]  vs r@',2+INQ)

         GYR(1:NRMAX,1) = SNGL(AJ(1:NRMAX))*1.E-6
         DO NR=1,NRMAX
            CALL TIMESPL(TSL,AJL,TMU,AJU(1,NR),NTXMAX,NTUM,IERR)
            GYR(NR,2) = SNGL(AJL*1.D-6)
         ENDDO
         CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,'@AJ(TR),AJ(XP) [MA/m$+2$=]  vs r@',2+INQ)

         IF(KUFDEV.EQ.'X'.AND.KUFDCG.EQ.'14') THEN
            GYR(1:NRMAX,1) = SNGL(AJBS(1:NRMAX))*1.E-6
            DO NR=1,NRMAX
               CALL TIMESPL(TSL,AJL,TMU,AJBSU(1,NR),NTXMAX,NTUM,IERR)
               GYR(NR,2) = SNGL(AJL*1.D-6)
            ENDDO
            CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2,'@AJBS(TR),AJBS(XP) [MA/m$+2$=]  vs r@',2+INQ)
         ELSE
            GYR(2:NRMAX+1,1) = SNGL(QP(1:NRMAX))
            DO NR=1,NRMAX
               CALL TIMESPL(TSL,QPL,TMU,QPU(1,NR),NTXMAX,NTUM,IERR)
               GYR(NR+1,2) = SNGL(QPL)
            ENDDO
            GYR(1,1) = SNGL((4.D0*QP(1)-QP(2))/3.D0)
            GYR(1,2) = (4.0*GYR(2,2)-GYR(3,2))/3.0
            CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRG,GYR,NRMP,NRMAX+1,2,'@QP(TR),QP(XP)  vs r@',2+INQ)
         ENDIF
         CALL TRGRTM
         CALL PAGEE
      end select

      RETURN
      END SUBROUTINE TRCMP3

!     **********************************************

!        COMPARE WITH UFILE DATA (ERROR BAR)

!     **********************************************

      SUBROUTINE TRCMP4(INQ)

      USE TRCOMM, ONLY : GER, GRM, GRM, GYR, MDLUF, NRM, NRMAX, NRMP, NRUM, NTUM, RN, RT
      USE TRCOM1, ONLY : GRE, NREMAX, RNEXEU, RNEXU, RTEXEU, RTEXU, RTIXEU, RTIXU
      IMPLICIT NONE
      INTEGER(4), INTENT(IN):: INQ
      INTEGER(4)  :: NR

      IF(MDLUF.EQ.2) THEN
         CALL PAGES
         GYR(1:NRMP,1:2)=0.0
         GYR(1:NREMAX(1),1) = SNGL(RTEXU(1,1:NREMAX(1)))
         GER(1:NREMAX(1),1) = SNGL(RTEXEU(1,1:NREMAX(1)))
         GYR(1:NRMAX,2) = SNGL(RT(1:NRMAX,1))
         CALL TRGR1DE( 3.0,12.0,11.0,17.0,GRM,GYR,GRE(1,1),GER,NRMP, &
     &                NRMAX,NREMAX(1),2,'@TE(XP),TE(TR) [keV]  vs r@',2+INQ)
!
         GYR(1:NREMAX(1),1) = SNGL(RTIXU(1,1:NREMAX(1)))
         GER(1:NREMAX(1),1) = SNGL(RTIXEU(1,1:NREMAX(1)))
         GYR(1:NRMAX,2) = SNGL(RT(1:NRMAX,2))
         CALL TRGR1DE(15.5,24.5,11.0,17.0,GRM,GYR,GRE(1,1),GER,NRMP, &
     &                NRMAX,NREMAX(1),2,'@TI(XP),TI(TR) [keV]  vs r@',2+INQ)
!
         GYR(1:NREMAX(2),1) = SNGL(RNEXU(1,1:NREMAX(2)))
         GER(1:NREMAX(2),1) = SNGL(RNEXEU(1,1:NREMAX(2)))
         GYR(1:NRMAX,2) = SNGL(RN(1:NRMAX,1))
         CALL TRGR1DE( 3.0,12.0, 2.0, 8.0,GRM,GYR,GRE(1,2),GER,NRMP, &
     &                NRMAX,NREMAX(2),2,'@NE(XP),NE(TR) [10$+20$=/m$+3$=]  vs r@',2+INQ)
         CALL TRGRTM
         CALL PAGEE
      ENDIF

      RETURN
      END SUBROUTINE TRCMP4

!     **********************************************

!        COMPARE WITH UFILE DATA (OTHER)

!     **********************************************

      SUBROUTINE TRCMP5(INQ)

      USE TRCOMM, ONLY : DT, GRM, GYR, MDLUF, NRMAX, NRMP, NT, NTUM, POH, POHU, RN, RNU_ORG, ZEFF, ZEFFU_ORG
      USE TRCOM1, ONLY : NTXMAX, TMU
      IMPLICIT NONE
      INTEGER(4), INTENT(IN):: INQ
      INTEGER(4):: NR, IERR
      REAL(8)   :: RN2L, RN4L, RNL, RTL, TSL

      IF(MDLUF.EQ.0) RETURN

      CALL PAGES
      TSL=DT*DBLE(NT)
      IF(MDLUF.EQ.2) THEN
         GYR(1:NRMAX,1) = SNGL(POH(1:NRMAX))*1.E-6
         GYR(1:NRMAX,2) = SNGL(POHU(1,1:NRMAX))*1.E-6
      ELSE
         GYR(1:NRMAX,1) = SNGL(POH(1:NRMAX))*1.E-6
         DO NR=1,NRMAX
            CALL TIMESPL(TSL,RTL,TMU,POHU(1,NR),NTXMAX,NTUM,IERR)
            GYR(NR,2) = SNGL(RTL*1.D-6)
         ENDDO
      ENDIF
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,'@POH(TR),POH(XP) [MW/m$+3$=]  vs r@',2+INQ)

      IF(MDLUF.EQ.2) THEN
         GYR(1:NRMAX,1) = SNGL(ZEFF(1:NRMAX))
         GYR(1:NRMAX,2) = SNGL(ZEFFU_ORG(1,1:NRMAX))
      ELSE
         GYR(1:NRMAX,1) = SNGL(ZEFF(1:NRMAX))
         DO NR=1,NRMAX
            CALL TIMESPL(TSL,RTL,TMU,ZEFFU_ORG(1,NR),NTXMAX,NTUM,IERR)
            GYR(NR,2) = SNGL(RTL)
         ENDDO
      ENDIF
      CALL TRGR1D(15.5,24.5,11.0,17.0,GRM,GYR,NRMP,NRMAX,2,'@ZEFF(TR),ZEFF(XP) [keV]  vs r@',2+INQ)

      IF(MDLUF.EQ.2) THEN
         GYR(1:NRMAX,1) = SNGL(RN(1:NRMAX,2))
         GYR(1:NRMAX,2) = SNGL(RNU_ORG(1,1:NRMAX,2)+RNU_ORG(1,1:NRMAX,4))
      ELSE
         GYR(1:NRMAX,1) = SNGL(RN(1:NRMAX,2))
         DO NR=1,NRMAX
            CALL TIMESPL(TSL,RN2L,TMU,RNU_ORG(1,NR,2),NTXMAX,NTUM,IERR)
            CALL TIMESPL(TSL,RN4L,TMU,RNU_ORG(1,NR,4),NTXMAX,NTUM,IERR)
            GYR(NR,2) = SNGL(RN2L+RN4L)
         ENDDO
      ENDIF
      CALL TRGR1D( 3.0,12.0, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2, &
     &     '@NI(TR), NI(XP) [10$+20$=/m$+3$=]  vs r@',2+INQ)

      IF(MDLUF.EQ.2) THEN
         GYR(1:NRMAX,1) = SNGL(RN(1:NRMAX,3))
         GYR(1:NRMAX,2) = SNGL(RNU_ORG(1,1:NRMAX,3))
      ELSE
         GYR(1:NRMAX,1) = SNGL(RN(1:NRMAX,3))
         DO NR=1,NRMAX
            CALL TIMESPL(TSL,RNL,TMU,RNU_ORG(1,NR,3),NTXMAX,NTUM,IERR)
            GYR(NR,2) = SNGL(RNL)
         ENDDO
      ENDIF
      CALL TRGR1D(15.5,24.5, 2.0, 8.0,GRM,GYR,NRMP,NRMAX,2, &
     &     '@Nimp(TR), Nimp(XP) [10$+20$=/m$+3$=]  vs r@',2+INQ)

      CALL PAGEE
      CALL TRGRTM

      RETURN
      END SUBROUTINE TRCMP5

!     ****************************************************

!        COMPARE WITH UFILE DATA (TEMPORAL; TEMPERATURE)

!     ****************************************************

      SUBROUTINE TRCMP6(INQ)

      USE TRCOMM, ONLY : G3D, GT, GVT, GYT, MDLUF, NGT, NRMAX, NTM, NTUM, RTU
      USE TRCOM1, ONLY : NTXMAX, TMU
      IMPLICIT NONE
      INTEGER(4), INTENT(IN):: INQ
      INTEGER(4):: IERR, NG, NR
      REAL(8)   :: PTE, PTI, TSL

      IF(MDLUF.EQ.0.OR.MDLUF.EQ.2) RETURN

      CALL PAGES

      DO NG=1,NGT
         GYT(NG,1)=GVT(NG, 9)
         TSL=DBLE(GT(NG))
         NR=1
         CALL TIMESPL(TSL,PTE,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
         GYT(NG,2)=SNGL(PTE)
         NR=INT(0.3D0*NRMAX)
         CALL TIMESPL(TSL,PTE,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
         GYT(NG,3)=G3D(NR,NG,1)
         GYT(NG,4)=SNGL(PTE)
         NR=INT(0.6D0*NRMAX)
         CALL TIMESPL(TSL,PTE,TMU,RTU(1,NR,1),NTXMAX,NTUM,IERR)
         GYT(NG,5)=G3D(NR,NG,1)
         GYT(NG,6)=SNGL(PTE)
      ENDDO
      CALL TRGR1D( 3.0,24.0, 9.7,17.0,GT,GYT,NTM,NGT,6,'@TE0(TR),TE0(XP) [keV]  vs t@',2+INQ)

      DO NG=1,NGT
         GYT(NG,1)=GVT(NG,10)
         TSL=DBLE(GT(NG))
         NR=1
         CALL TIMESPL(TSL,PTI,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
         GYT(NG,2)=SNGL(PTI)
         NR=INT(0.3D0*NRMAX)
         CALL TIMESPL(TSL,PTI,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
         GYT(NG,3)=G3D(NR,NG,2)
         GYT(NG,4)=SNGL(PTI)
         NR=INT(0.6D0*NRMAX)
         CALL TIMESPL(TSL,PTI,TMU,RTU(1,NR,2),NTXMAX,NTUM,IERR)
         GYT(NG,5)=G3D(NR,NG,2)
         GYT(NG,6)=SNGL(PTI)
      ENDDO
      CALL TRGR1D( 3.0,24.0, 1.1, 8.4,GT,GYT,NTM,NGT,6,'@TI0(TR),TI0(XP) [keV]  vs t@',2+INQ)

      CALL TRGRTM
      CALL PAGEE

      RETURN
      END SUBROUTINE TRCMP6

!     ****************************************************

!        COMPARE WITH RADIAL ELECTRIC FIELD MODELS

!     ****************************************************

      SUBROUTINE TRCMP7(INQ)

      USE TRCOMM, ONLY : ER, GRG, GYR, MDLER, NRMP, NRMAX
      IMPLICIT NONE
      INTEGER(4), INTENT(IN):: INQ
      INTEGER(4):: MDLER_ORG, NR
      REAL(8),DIMENSION(NRMAX):: ER_ORG

      CALL PAGES

      MDLER_ORG=MDLER
      ER_ORG(1:NRMAX)=ER(1:NRMAX)
      DO MDLER=1,4
         CALL TRERAD
         GYR(1:NRMAX,MDLER)=SNGL(ER(1:NRMAX))
      ENDDO
      CALL TRGR1D( 3.0,12.0,11.0,17.0,GRG,GYR,NRMP,NRMAX,4,'@ER  vs r@',2+INQ)

      MDLER=MDLER_ORG
      ER(1:NRMAX)=ER_ORG(1:NRMAX)

      CALL PAGEE
      CALL TRGRTM

      RETURN
      END SUBROUTINE TRCMP7
