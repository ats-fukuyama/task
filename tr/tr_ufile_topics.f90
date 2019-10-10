!     *****************************************************************

!           TIME EVOLUTIONAL UFILE DATA INPUT ROUTINE FOR TOPICS

!     *****************************************************************

      SUBROUTINE TR_TIME_UFILE_TOPICS

      USE TRCOMM, ONLY : ABRHOU, AJBSU, AJU, AR1RHOU, AR2RHOU, ARRHOU, BB, BBU, DR, DT, DVRHOU, EPSRHO, KUFDCG, KUFDEV, &
     &                   MDLXP, MDPHIA, NRAMAX, NRMAX, NRMP, NRUM, NTUM, PECU, PHIA, PI, PICU, PNBU, PNSU, PNSUA, POHU, &
     &                   PRLU, PZ, QPU, RA, RAU, RG, RHOA, RIPE, RIPS, RIPU, RKAP, RKAPU, RKPRHO, RMJRHOU, RMNRHOU, RNU, &
     &                   RR, RRU, RTU, TTRHOU, VOLAU, ZEFFU
      USE TRCOM1, ONLY : NTAMAX, NTXMAX, NTXMAX1, TMU, TMU1
      IMPLICIT NONE
      INTEGER(4):: ICK, IERR, MDAMIN, MDBT, MDIP, MDKAPPA, MDRGEO, MDVOL, NTAMAX1, NTX
      INTEGER(4):: MDALL, MDCHK, NR, NRLMAX
      REAL(8)   :: AMP, TMUMAX
      REAL(8),DIMENSION(NTUM)     :: PV, PVA, TL, VOL
      REAL(8),DIMENSION(NRUM)     :: RL
      REAL(8),DIMENSION(NTUM,NRMP):: FAT
      REAL(8),DIMENSION(NTUM,NRUM):: F2
      CHARACTER(LEN=10):: KFID
      CHARACTER(LEN=80):: KDIRX


      ICK=0
      TMUMAX=0.D0

      MDRGEO=0
      MDAMIN=0
      MDIP=0
      MDBT=0
      MDKAPPA=0
      MDVOL=0
      MDPHIA=0

!     *** 1D VALUE ***

      KFID='RGEO'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,RRU  ,NTAMAX1,NTXMAX1,TMUMAX,ICK,MDLXP,IERR)
      MDRGEO=IERR

      KFID='AMIN'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,RAU  ,NTAMAX1,NTXMAX1,TMUMAX,ICK,MDLXP,IERR)
      MDAMIN=IERR

      KFID='IP'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,RIPU ,NTAMAX1,NTXMAX1,TMUMAX,ICK,MDLXP,IERR)
      RIPU(1:NTXMAX1)=ABS(RIPU(1:NTXMAX1)*1.D-6)
      MDIP=IERR

      KFID='BT'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,BBU  ,NTAMAX1,NTXMAX1,TMUMAX,ICK,MDLXP,IERR)
      MDBT=IERR

      KFID='KAPPA'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,RKAPU,NTAMAX1,NTXMAX1,TMUMAX,ICK,MDLXP,IERR)
      MDKAPPA=IERR

      KFID='VOL'
      CALL UF1D(KFID,KUFDEV,KUFDCG,DT,TMU1,VOL, NTAMAX1,NTXMAX1,TMUMAX,ICK,MDLXP,IERR)
      MDVOL=IERR

      RR    = RRU(1)
      RA    = RAU(1)
      RIPS  = RIPU(1)
      RIPE  = RIPU(1)
      BB    = BBU(1)
      RKAP  = RKAPU(1)
      PHIA  = 0.D0

!     *** 2D VALUE ***

      AMP=1.D-3
      KFID='TE'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)
      RTU(1:NTXMAX,1:NRMAX,1)=FAT(1:NTXMAX,1:NRMAX)

      KFID='TI'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)
      RTU(1:NTXMAX,1:NRMAX,2)=FAT(1:NTXMAX,1:NRMAX)
      RTU(1:NTXMAX,1:NRMAX,3)=FAT(1:NTXMAX,1:NRMAX)

      AMP=1.D-20
      KFID='NE'
      CALL UF2DTP(KFID,KUFDEV,KUFDCG,DR,DT,PV,PVA,TMU,FAT,AMP, &
     &          NTAMAX,NTXMAX,RHOA,NRAMAX,NRMAX,TMUMAX,0,ICK,MDLXP,IERR)
      RNU(1:NTXMAX,1:NRMAX,1)=FAT(1:NTXMAX,1:NRMAX)
!      RNU(1:NTXMAX,1:NRMAX,2)=FAT(1:NTXMAX,1:NRMAX)
      PNSU (1:NTXMAX,1)=PV (1:NTXMAX)
!      PNSU (1:NTXMAX,2)=PV (1:NTXMAX)
      IF(RHOA.NE.1.D0) THEN
         PNSUA(1:NTXMAX,1)=PVA(1:NTXMAX)
!         PNSUA(1:NTXMAX,2)=PVA(1:NTXMAX)
      ENDIF

      KFID='NIMP'
      CALL UF2DTP(KFID,KUFDEV,KUFDCG,DR,DT,PV,PVA,TMU,FAT,AMP,NTAMAX,NTXMAX,RHOA,NRAMAX,NRMAX,TMUMAX,0,ICK,MDLXP,IERR)
      RNU(1:NTXMAX,1:NRMAX,2)=(RNU(1:NTXMAX,1:NRMAX,1)-PZ(3)*FAT(1:NTXMAX,1:NRMAX))/PZ(2)
      RNU(1:NTXMAX,1:NRMAX,3)=FAT(1:NTXMAX,1:NRMAX)
      PNSU(1:NTXMAX,2)=(PNSU(1:NTXMAX,1)-PZ(3)*PV(1:NTXMAX))/PZ(2)
      PNSU(1:NTXMAX,3)=PV(1:NTXMAX)
      IF(RHOA.NE.1.D0) THEN
         PNSUA(1:NTXMAX,2)=(PNSUA(1:NTXMAX,1)-PZ(3)*PVA(1:NTXMAX))/PZ(2)
         PNSUA(1:NTXMAX,3)=PVA(1:NTXMAX)
      ENDIF

      KFID='NM1'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)
      RNU(1:NTXMAX,1:NRMAX,4)=FAT(1:NTXMAX,1:NRMAX)

      AMP=1.D0
      KFID='ZEFFR'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,1,1,ICK,MDLXP,IERR)
      ZEFFU(1:NTXMAX,1:NRMAX)=FAT(1:NTXMAX,1:NRMAX)

      KFID='Q'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,QPU,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,1,1,ICK,MDLXP,IERR)
      KFID='CURTOT'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,AJU,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)
      IF(KUFDEV.EQ.'X'.AND.KUFDCG.EQ.'14') THEN
         KFID='CURBS'
         CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,AJBSU,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)
      ENDIF

      KFID='QNBIE'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)
      forall(NTX=1:NTXMAX,NR=1:NRMAX,FAT(NTX,NR) < 0.D0) FAT(NTX,NR) = 0.D0
      PNBU(1:NTXMAX,1:NRMAX,1) = FAT(1:NTXMAX,1:NRMAX)

      KFID='QNBII'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)
      forall(NTX=1:NTXMAX,NR=1:NRMAX,FAT(NTX,NR) < 0.D0) FAT(NTX,NR) = 0.D0
      PNBU(1:NTXMAX,1:NRMAX,2) = FAT(1:NTXMAX,1:NRMAX)

      PICU(1:NTXMAX,1:NRMAX,1)=0.D0
      PICU(1:NTXMAX,1:NRMAX,2)=0.D0
      PECU(1:NTXMAX,1:NRMAX  )=0.D0
      POHU(1:NTXMAX,1:NRMAX  )=0.D0
!      KFID='QICRHE'
      KFID='QLHE'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)
      forall(NTX=1:NTXMAX,NR=1:NRMAX,FAT(NTX,NR) < 0.D0) FAT(NTX,NR) = 0.D0
      PICU(1:NTXMAX,1:NRMAX,1) = FAT(1:NTXMAX,1:NRMAX)

!      KFID='QICRHI'
      KFID='QLHI'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)
      forall(NTX=1:NTXMAX,NR=1:NRMAX,FAT(NTX,NR) < 0.D0) FAT(NTX,NR) = 0.D0
      PICU(1:NTXMAX,1:NRMAX,2) = FAT(1:NTXMAX,1:NRMAX)

      KFID='QECH'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)
      forall(NTX=1:NTXMAX,NR=1:NRMAX,FAT(NTX,NR) < 0.D0) FAT(NTX,NR) = 0.D0
      PECU(1:NTXMAX,1:NRMAX) = FAT(1:NTXMAX,1:NRMAX)

      KFID='QRAD'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,PRLU,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)

      KFID='QOHM'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,POHU,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)

!     *** GEOMETRY FACTORS ***

      KFID='RMAJOR'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,RMJRHOU,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,1,1,ICK,MDLXP,IERR)

      KFID='RMINOR'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,RMNRHOU,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,1,0,ICK,MDLXP,IERR)

      KFID='GRHO1'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,AR1RHOU,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)

      KFID='GRHO2'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,AR2RHOU,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)

      KFID='VOLUME'
      CALL UFREAD2_TIME(KDIRX,KUFDEV,KUFDCG,KFID,RL,TL,F2,NRLMAX,NTXMAX,NRUM,NTUM,MDCHK,IERR)
      VOLAU(1:NTXMAX)=F2(1:NTXMAX,NRLMAX)

      KFID='VRO'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,DVRHOU,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,0,ICK,MDLXP,IERR)

      KFID='AAT'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,ARRHOU,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)

      KFID='HDT'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)
      TTRHOU(1:NTXMAX,1:NRMAX)=FAT(1:NTXMAX,1:NRMAX)/ARRHOU(1:NTXMAX,1:NRMAX)

      KFID='CKT'
      CALL UF2DT(KFID,KUFDEV,KUFDCG,DR,DT,TMU,FAT,AMP,NTAMAX,NTXMAX,NRMAX,TMUMAX,0,1,ICK,MDLXP,IERR)
      ABRHOU(1:NTXMAX,1:NRMAX)=FAT(1:NTXMAX,1:NRMAX)/DVRHOU(1:NTXMAX,1:NRMAX)**2

      DO NR=1,NRMAX
         RG(NR)    =DBLE(NR)*DR
      ENDDO
      EPSRHO(1:NRMAX)=RMNRHOU(1,1:NRMAX)/RMJRHOU(1,1:NRMAX)
      RKPRHO(1:NRMAX)=RKAP

!     *** 1D VALUE ***

      MDALL=MDRGEO+MDAMIN+MDIP+MDBT+MDKAPPA
      IF(MDALL.NE.0) THEN
         NTXMAX1=NTXMAX
         IF(MDRGEO.NE.0)  RRU(1:NTXMAX1)=RR
         IF(MDAMIN.NE.0)  RAU(1:NTXMAX1)=RA
         IF(MDIP.NE.0)    RIPU(1:NTXMAX1)=RIPS
         IF(MDBT.NE.0)    BBU(1:NTXMAX1)=BB
         IF(MDKAPPA.NE.0) RKAPU(1:NTXMAX1)=RKAP
         IF(MDVOL.NE.0)   VOL(1:NTXMAX1)=PI*RKAP*RA**2*2.D0*PI*RR
         TMU1(1:NTXMAX1)=TMU(1:NTXMAX1)
      ENDIF

!     ***

      RETURN
      END SUBROUTINE TR_TIME_UFILE_TOPICS
