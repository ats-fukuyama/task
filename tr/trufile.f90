! trufile.f90
!     ***********************************************************

!           CONTROL ROUTINE FOR UFILE READING

!     ***********************************************************

      SUBROUTINE TR_UFILE_CONTROL(NSW)

!     MDNI is a parameter that can control which data is used
!     to determine bulk density, impurity density or effective
!     charge number among those data.
!     If all the data above do not exist, MDNI is set to zero
!     automatically regardless of original MDNI.
!           0 : NSMAX=2, ne=ni
!           1 : calculate nimp and ni profiles from NE, ZIMP and ZEFFR
!           2 : calculate nimp and zeff profiles from NE, ZIMP and NM1
!           3 : calculate zeff and ni profiles from NE, ZIMP and NIMP

!     For the time being, we deal with transport calculation with
!     three particles so that densities of 4th argument requires only
!     finite small values and their values can be ignored in the context
!     of charge neutrality.

      USE TRCOMM, ONLY : &
           AD0, CDP, CNP, KUFDCG, KUFDEV, MDLFLX, MDNI, NEQMAX, NRAMAX, &
           NRMAX, NROMAX, NSMAX, NST
      
      IMPLICIT NONE
      INTEGER, INTENT(IN):: NSW
      INTEGER:: MDSLCT,  NEQ, NEQL

      IF(NSW.NE.0) THEN
         CALL CHECK_IMPURITY(MDSLCT)

         select case(MDSLCT)
         case(0)
            MDNI=0
            NEQL=0
            DO NEQ=1,NEQMAX
               IF(NST(NEQ).NE.0) NEQL=NEQL+1
            ENDDO
            IF(NEQL.EQ.1) THEN
               NSMAX=1
            ELSE
               NSMAX=2
            ENDIF
         case(1)
            IF(MDNI.EQ.2.OR.MDNI.EQ.3) MDNI=1
         case(2)
            IF(MDNI.EQ.1.OR.MDNI.EQ.3) MDNI=2
         case(3)
            IF(MDNI.EQ.3) MDNI=1
         case(4)
            IF(MDNI.EQ.1.OR.MDNI.EQ.2) MDNI=3
         case(5)
            IF(MDNI.EQ.2) MDNI=1
         case(6)
            IF(MDNI.EQ.1) MDNI=2
         end select

         IF(MDLFLX.NE.0) THEN
            CNP=0.D0
            CDP=0.D0
            AD0=0.D0
         ENDIF

         WRITE(6,'(A7,A8,A6,A15)') 'DEVICE=',KUFDEV,'SHOT#=',KUFDCG

         select case(NSW)
         case(1)
            CALL TR_TIME_UFILE
         case(2)
            CALL TR_STEADY_UFILE
         case(3)
            CALL TR_TIME_UFILE_TOPICS
         end select
      ENDIF

      NROMAX = NRMAX
      NRMAX  = NRAMAX

      RETURN
      END SUBROUTINE TR_UFILE_CONTROL

!     ***********************************************************

!           CHECKING WHETHER IMPURITY EXISTS

!     ***********************************************************

      SUBROUTINE CHECK_IMPURITY(MDSLCT)

      USE TRCOMM, ONLY : KUFDCG, KUFDEV, MDNI
      USE TRCOM1, ONLY : KDIRX, NMCHK
      IMPLICIT NONE
      INTEGER, INTENT(OUT):: MDSLCT
      INTEGER             :: IKDIRX, IKNDCG, IKNDEV, KL2
      CHARACTER(LEN=80)      :: KDIRR2, KFILE
      CHARACTER(LEN=20)      :: KFID
      LOGICAL LEX


      IKNDEV = len_trim(KUFDEV)
      IKNDCG = len_trim(KUFDCG)

      IKDIRX = len_trim(KDIRX)
      KDIRR2=KDIRX(1:IKDIRX)//KUFDEV(1:IKNDEV)//'2d'//KUFDCG(1:IKNDCG)//'.'
      KL2 = len_trim(KDIRR2)

      IF(MDNI.LT.0.OR.MDNI.GT.3) MDNI=0
      MDSLCT=0

      KFID='ZEFFR'
      KFILE=KDIRR2(1:KL2)//KFID
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
      IF(LEX) MDSLCT=MDSLCT+1
      KFID='NM1'
      KFILE=KDIRR2(1:KL2)//KFID
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
      IF(LEX) THEN
         NMCHK=0
         MDSLCT=MDSLCT+2
         KFID='NM2'
         KFILE=KDIRR2(1:KL2)//KFID
         INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
         IF(LEX) NMCHK=1
         KFID='NM3'
         KFILE=KDIRR2(1:KL2)//KFID
         INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
         IF(LEX) NMCHK=2
      ENDIF
      KFID='NIMP'
      KFILE=KDIRR2(1:KL2)//KFID
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
      IF(LEX) MDSLCT=MDSLCT+4

 9000 RETURN
      END SUBROUTINE CHECK_IMPURITY

!     ***********************************************************

!           STEADY STATE UFILE DATA INPUT ROUTINE

!     ***********************************************************

      SUBROUTINE TR_STEADY_UFILE

      USE TRCOMM, ONLY : ABRHOU, AJNBU, AJU, AR1RHOU, AR2RHOU, ARRHOU, BB, BBU, BPU, DR, DVRHOU, KUFDCG, KUFDEV, MDLJQ, MDLXP, &
     &                   MDNI, MDPHIA, NRAMAX, NRMAX, NRMP, NRUM, NTS, NTUM, PBMU, PECU, PHIA, PHIAU, PICU, PN, PNBU, PNS,&
     &                   PNSA, POHU, PRLU, PT, PTS, PTSA, PZ, QPU, RA, RAU, RHOA, RIPE, RIPS, RIPU, RKAP, RKAPU, RKPRHOU,      &
     &                   RMJRHOU, RMNRHOU, RNFU, RNU, RNU_ORG, RR, RRU, RT, RTU, SNBU, SWLU, TIME_INT, TTRHOU, VOLAU, WROTU,   &
     &                   ZEFFU, ZEFFU_ORG, rkind
      USE TRCOM1, ONLY : GRE, KDIRX, NMCHK, NREMAX, NTXMAX, NTXMAX1, RNEXEU, RNEXU, RTEXEU, RTEXU, RTIXEU, RTIXU, TMU, TMU1
      IMPLICIT NONE
      INTEGER:: IERR, MDAMIN, MDBT, MDCHK, MDIP, MDKAPPA, MDRGEO, MDSUM, NR, NRFMAX, NRLMAX, NTX, NTX1
      REAL(rkind)   :: AMP, DT_MIN, PNF1, PNF2, PNIMP1, PNIMP2, PNM1, PNM2, PNM3, PNM4, PNM5, PNM6, PTMP1, PTMP2, PZEF1, PZEF2
      REAL(rkind),DIMENSION(NRUM) :: RUF, RL
      REAL(rkind),DIMENSION(NTUM) :: F1, TL
      REAL(rkind),DIMENSION(NRMP) :: FAS
      REAL(rkind),DIMENSION(NTUM,NRUM)::  F2
      CHARACTER(LEN=10)       :: KFID


      MDRGEO=0
      MDAMIN=0
      MDIP=0
      MDBT=0
      MDKAPPA=0
      MDPHIA=0

      NTS=0

!     *** 1D VALUE ***

      KFID='RGEO'
      IF(MDLXP.EQ.0) THEN
         CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TMU1,F1,NTXMAX1,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_MDS1(KUFDEV,KUFDCG,KFID,NTUM,TMU1,F1,NTXMAX1,IERR)
      ENDIF
      IF(IERR.EQ.1) THEN
         MDRGEO=IERR
         IERR=0
      ELSE
         RRU(1:NTXMAX1)=F1(1:NTXMAX1)
      ENDIF

      KFID='AMIN'
      IF(MDLXP.EQ.0) THEN
         CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TMU1,F1,NTXMAX1,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_MDS1(KUFDEV,KUFDCG,KFID,NTUM,TMU1,F1,NTXMAX1,IERR)
      ENDIF
      IF(IERR.EQ.1) THEN
         MDAMIN=IERR
         IERR=0
      ELSE
         RAU(1:NTXMAX1)=F1(1:NTXMAX1)
      ENDIF

      KFID='IP'
      IF(MDLXP.EQ.0) THEN
         CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TMU1,F1,NTXMAX1,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_MDS1(KUFDEV,KUFDCG,KFID,NTUM,TMU1,F1,NTXMAX1,IERR)
      ENDIF
      IF(IERR.EQ.1) THEN
         MDIP=IERR
         IERR=0
      ELSE
         RIPU(1:NTXMAX1)=ABS(F1(1:NTXMAX1)*1.D-6)
      ENDIF

      KFID='BT'
      IF(MDLXP.EQ.0) THEN
         CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TMU1,F1,NTXMAX1,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_MDS1(KUFDEV,KUFDCG,KFID,NTUM,TMU1,F1,NTXMAX1,IERR)
      ENDIF
      IF(IERR.EQ.1) THEN
         MDBT=IERR
         IERR=0
      ELSE
         BBU(1:NTXMAX1)=ABS(F1(1:NTXMAX1))
      ENDIF

      IF(KUFDEV.EQ.'tftr'.AND.(KUFDCG.EQ.'50862' .OR. KUFDCG.EQ.'50921'.OR.KUFDCG.EQ.'52527')) THEN
         RKAPU(1:NTXMAX1)=1.D0
      ELSE
         KFID='KAPPA'
         IF(MDLXP.EQ.0) THEN
            CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TMU1,F1,NTXMAX1,NTUM,MDCHK,IERR)
         ELSE
            CALL IPDB_MDS1(KUFDEV,KUFDCG,KFID,NTUM,TMU1,F1,NTXMAX1,IERR)
         ENDIF
         IF(IERR.EQ.1) THEN
            MDKAPPA=IERR
            IERR=0
         ELSE
            RKAPU(1:NTXMAX1)=F1(1:NTXMAX1)
         ENDIF
      ENDIF
!      IF(IERR.NE.0) STOP 'SOME 1D UFILES DO NOT EXIST.'

      KFID='PHIA'
      IF(MDLXP.EQ.0) THEN
         CALL UFREAD_TIME(KDIRX,KUFDEV,KUFDCG,KFID,TMU1,F1,NTXMAX1,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_MDS1(KUFDEV,KUFDCG,KFID,NTUM,TMU1,F1,NTXMAX1,IERR)
      ENDIF
      IF(IERR.EQ.1) THEN
         MDPHIA=IERR
         IERR=0
      ELSE
         PHIAU(1:NTXMAX1)=F1(1:NTXMAX1)
      ENDIF

!     *** 2D VALUE ***

      AMP=1.D-3
      KFID='TE'
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PTMP1,PTMP2,FAS,AMP,RHOA,NRAMAX,NRMAX,MDLXP,IERR)
      RTU(1,1:NRMAX,1)=FAS(1:NRMAX)
      RT(1:NRMAX,1)=FAS(1:NRMAX)
      PT(1)=RT(1,1)
      PTS(1)=PTMP1
      IF(RHOA.NE.1.D0) PTSA(1)=PTMP2

      KFID='TEXP'
      IF(MDLXP.EQ.0) THEN
         CALL UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,RTEXU,NRFMAX,NTXMAX,NRUM,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_RAW2(KUFDEV,KUFDCG,KFID,NRUM,NTUM,RUF,TMU,RTEXU,NRFMAX,NTXMAX,IERR)
      ENDIF
      IF(IERR.EQ.0) THEN
         NREMAX(1)=NRFMAX
         GRE(1:NRFMAX,1)=SNGL(RUF(1:NRFMAX))
      ENDIF
      KFID='TEXPEB'
      IF(MDLXP.EQ.0) THEN
         CALL UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,RTEXEU,NRFMAX,NTXMAX,NRUM,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_RAW2(KUFDEV,KUFDCG,KFID,NRUM,NTUM,RUF,TMU,RTEXEU,NRFMAX,NTXMAX,IERR)
      ENDIF
!
      KFID='TI'
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PTMP1,PTMP2,FAS,AMP,RHOA,NRAMAX,NRMAX,MDLXP,IERR)
      RTU(1,1:NRMAX,2)=FAS(1:NRMAX)
      RTU(1,1:NRMAX,3)=FAS(1:NRMAX)
      RT(1:NRMAX,2)=FAS(1:NRMAX)
      RT(1:NRMAX,3)=FAS(1:NRMAX)
      PT  (2)=RT(1,2)
      PT  (3)=RT(1,3)
      PT  (4)=RT(1,2)
      PTS (2)=PTMP1
      PTS (3)=PTMP1
      PTS (4)=PTMP1
      IF(RHOA.NE.1.D0) THEN
         PTSA(2)=PTMP2
         PTSA(3)=PTMP2
         PTSA(4)=PTMP2
      ENDIF

      KFID='TIXP'
      IF(MDLXP.EQ.0) THEN
         CALL UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,RTIXU,NRFMAX,NTXMAX,NRUM,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_RAW2(KUFDEV,KUFDCG,KFID,NRUM,NTUM,RUF,TMU,RTIXU,NRFMAX,NTXMAX,IERR)
      ENDIF
      KFID='TIXPEB'
      IF(MDLXP.EQ.0) THEN
         CALL UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,RTIXEU,NRFMAX,NTXMAX,NRUM,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_RAW2(KUFDEV,KUFDCG,KFID,NRUM,NTUM,RUF,TMU,RTIXEU,NRFMAX,NTXMAX,IERR)
      ENDIF

      AMP=1.D-20
      KFID='NE'
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PTMP1,PTMP2,FAS,AMP,RHOA,NRAMAX,NRMAX,MDLXP,IERR)
      RNU(1,1:NRMAX,1)     = FAS(1:NRMAX)
      RNU(1,1:NRMAX,2)     = FAS(1:NRMAX)
      RNU_ORG(1,1:NRMAX,1) = FAS(1:NRMAX)
      RNU_ORG(1,1:NRMAX,2) = FAS(1:NRMAX)
      ZEFFU(1,1:NRMAX)     = 1.D0 ! in case of MDNI=0
      ZEFFU_ORG(1,1:NRMAX) = 1.D0 ! in case of MDNI=0
      PN(1)  =RNU(1,1,1)
      PN(2)  =RNU(1,1,2)
      PN(3)  =1.D-7
      PN(4)  =1.D-7
      PNS(1) =PTMP1
      PNS(2) =PNS(1)
      PNS(3) =1.D-8
      PNS(4) =1.D-8
      IF(RHOA.NE.1.D0) THEN
         PNSA(1)=PTMP2
         PNSA(2)=PNSA(1)
         PNSA(3)=1.D-8
         PNSA(4)=1.D-8
      ENDIF

      KFID='NEXP'
      IF(MDLXP.EQ.0) THEN
         CALL UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,RNEXU,NRFMAX,NTXMAX,NRUM,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_RAW2(KUFDEV,KUFDCG,KFID,NRUM,NTUM,RUF,TMU,RNEXU,NRFMAX,NTXMAX,IERR)
      ENDIF
      IF(IERR.EQ.0) THEN
         NREMAX(2)=NRFMAX
         GRE(1:NRFMAX,2)=SNGL(RUF(1:NRFMAX))
      ENDIF
      KFID='NEXPEB'
      IF(MDLXP.EQ.0) THEN
         CALL UFREAD2_ERROR(KDIRX,KUFDEV,KUFDCG,KFID,RUF,TMU,RNEXEU,NRFMAX,NTXMAX,NRUM,NTUM,MDCHK,IERR)
      ELSE
         CALL IPDB_RAW2(KUFDEV,KUFDCG,KFID,NRUM,NTUM,RUF,TMU,RNEXEU,NRFMAX,NTXMAX,IERR)
      ENDIF

      KFID='NFAST'
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PNF1,PNF2,FAS,AMP,RHOA,NRMAX,NRMAX,MDLXP,IERR)
      DO NR=1,NRMAX
         IF(IERR.EQ.0.AND.FAS(NR).LE.0.D0) FAS(NR)=1.D-10
         RNFU(1,NR)=FAS(NR)
      ENDDO

 100  CONTINUE
      KFID='NM1'
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PNM1,PNM2,FAS,AMP,RHOA,NRAMAX,NRMAX,MDLXP,IERR)
      DO NR=1,NRMAX
         IF(NMCHK.EQ.0) THEN
            IF(MDNI.NE.0) RNU(1,NR,2)=FAS(NR)
            IF(IERR.EQ.0) RNU_ORG(1,NR,2)=FAS(NR)
            RNU_ORG(1,NR,4)=0.D0
         ELSE
            RNU_ORG(1,NR,4)=FAS(NR)
         ENDIF
      ENDDO
      IF(NMCHK.GE.1) THEN
         KFID='NM2'
         CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PNM3,PNM4,FAS,AMP,RHOA,NRAMAX,NRMAX,MDLXP,IERR)
         IF(MDNI.NE.0) RNU(1,1:NRMAX,2)=FAS(1:NRMAX)+RNU_ORG(1,1:NRMAX,4)
         RNU_ORG(1,1:NRMAX,2)=FAS(1:NRMAX)
         PNM1=PNM1+PNM3
         PNM2=PNM2+PNM4
      ENDIF
      IF(NMCHK.GE.2) THEN
         KFID='NM3'
         CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PNM5,PNM6,FAS,AMP,RHOA,NRAMAX,NRMAX,MDLXP,IERR)
         IF(MDNI.NE.0) RNU(1,1:NRMAX,2)=RNU(1,1:NRMAX,2)+FAS(1:NRMAX)
         RNU_ORG(1,1:NRMAX,4)=RNU_ORG(1,1:NRMAX,4)+FAS(1:NRMAX)
         PNM1=PNM1+PNM5
         PNM2=PNM2+PNM6
      ENDIF

      AMP=1.D0
      KFID='ZEFFR'
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PZEF1,PZEF2,FAS,AMP,RHOA,NRAMAX,NRMAX,MDLXP,IERR)
      IF(MDNI.NE.0) ZEFFU(1,1:NRMAX)=FAS(1:NRMAX)
      ZEFFU_ORG(1,1:NRMAX)=FAS(1:NRMAX)

      AMP=1.D-20
      KFID='NIMP'
      CALL UF2DSP(KFID,KUFDEV,KUFDCG,DR,PNIMP1,PNIMP2,FAS,AMP,RHOA,NRAMAX,NRMAX,MDLXP,IERR)
      IF(MDNI.NE.0) RNU(1,1:NRMAX,3)=FAS(1:NRMAX)
      RNU_ORG(1,1:NRMAX,3)=FAS(1:NRMAX)

      select case(MDNI)
      case(1)
         
      DO NR=1,NRMAX
         RNU(1,NR,2)=( (PZ(3)-ZEFFU(1,NR))*RNU(1,NR,1)-(PZ(3)-PZ(2))*PZ(2)*RNFU(1,NR))/(PZ(2)*(PZ(3)-PZ(2)))
         RNU(1,NR,3)=(PZ(2)-ZEFFU(1,NR))/(PZ(3)*(PZ(2)-PZ(3)))*RNU(1,NR,1)
      ENDDO
      PN(2)=RNU(1,1,2)
      PN(3)=RNU(1,1,3)
      PN(4)=1.D-7
      PNS(2)=( (PZ(3)-PZEF1)*PNS(1)-(PZ(3)-PZ(2))*PZ(2)*PNF1)/(PZ(2)*(PZ(3)-PZ(2)))
      PNS(3)=(PZ(2)-PZEF1)/(PZ(3)*(PZ(2)-PZ(3)))*PNS(1)
      PNS(4)=1.D-8
      IF(RHOA.NE.1.D0) THEN
         PNSA(2)=( (PZ(3)-PZEF2)*PNSA(1)-(PZ(3)-PZ(2))*PZ(2)*PNF2)/(PZ(2)*(PZ(3)-PZ(2)))
         PNSA(3)=(PZ(2)-PZEF2)/(PZ(3)*(PZ(2)-PZ(3)))*PNSA(1)
         PNSA(4)=1.D-8
      ENDIF

      case(2)

      DO NR=1,NRMAX
         RNU(1,NR,3)=(RNU(1,NR,1)-( RNU(1,NR,2)+RNFU(1,NR))*PZ(2))/ PZ(3)
         ZEFFU(1,NR)=PZ(2)*(PZ(2)-PZ(3))*(RNU(1,NR,2)+RNFU(1,NR)) / RNU(1,NR,1)+PZ(3)
      ENDDO
      PN(2) =RNU(1,1,2)
      PN(3) =RNU(1,1,3)
      PN(4) =1.D-7
      IF(PNM1.LE.1.D-3) THEN
         PNS(1)=1.D-2
         PNM1=1.D-2-1.D-3
      ENDIF
      PNS(2)=PNM1
      PNS(3)=(PNS(1)-(PNM1+PNF1)*PZ(2))/PZ(3)
      PNS(4)=1.D-8
      IF(RHOA.NE.1.D0) THEN
         PNSA(2)=PNM2
         PNSA(3)=(PNSA(1)-(PNM2+PNF2)*PZ(2))/PZ(3)
         PNSA(4)=1.D-8
      ENDIF

      case(3)

      DO NR=1,NRMAX
         RNU(1,NR,2)=(RNU(1,NR,1)-PZ(2)*RNFU(1,NR)-PZ(3)*RNU(1,NR,3)) /PZ(2)
         ZEFFU(1,NR)=PZ(2)+PZ(3)*(PZ(3)-PZ(2))*(RNU(1,NR,3)/RNU(1,NR,1))
      ENDDO
      PNS(2)=(PNS(1)-PZ(2)*PNF1-PZ(3)*PNIMP1)/PZ(2)
      PNS(3)=PNIMP1
      PNS(4)=1.D-8
      IF(RHOA.NE.1.D0) THEN
         PNSA(2)=(PNSA(1)-PZ(2)*PNF2-PZ(3)*PNIMP2)/PZ(2)
         PNSA(3)=PNIMP2
         PNSA(4)=1.D-8
      ENDIF

      end select
      IF(PNS(3).LE.1.D-8) PNS(3)=1.D-8

      IF(MDNI.NE.0) THEN
!     check whether density developed is appropriate (positive) or not
      DO NR=1,NRMAX
         IF(RNU(1,NR,2).LE.0.D0.OR.RNU(1,NR,3).LE.0.D0) THEN
            WRITE(6,*)'XX TR_STEADY_UFILE: DENSITY NEGATIVE: WRONG MDNI=',MDNI
            MDNI=MDNI+1
            IF(MDNI.LE.3) THEN
               GOTO 100
            ELSE
               STOP
            ENDIF
         ENDIF
      ENDDO
      ENDIF

      AMP=1.D0
      KFID='PBEAM'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      PBMU(1,1:NRMAX)=FAS(1:NRMAX)

      KFID='Q'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,1,1,MDLXP,IERR)
      IF(IERR.NE.0.AND.MDLJQ.EQ.1) MDLJQ=0
      QPU(1,1:NRMAX)=FAS(1:NRMAX)

      KFID='CURTOT'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      IF(IERR.NE.0.AND.MDLJQ.EQ.0) MDLJQ=1
      AJU(1,1:NRMAX)=FAS(1:NRMAX)

      KFID='CURNBI'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      AJNBU(1,1:NRMAX)=FAS(1:NRMAX)

      KFID='BPOL'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,1,0,MDLXP,IERR)
      BPU(1,1:NRMAX)=FAS(1:NRMAX)

      KFID='QNBIE'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      forall(NR=1:NRMAX,FAS(NR) < 0.d0) FAS(NR) = 0.D0
      PNBU(1,1:NRMAX,1) = FAS(1:NRMAX)

      KFID='QNBII'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      forall(NR=1:NRMAX,FAS(NR) < 0.d0) FAS(NR) = 0.D0
      PNBU(1,1:NRMAX,2) = FAS(1:NRMAX)

      PICU(1,1:NRMAX,1)=0.D0
      PICU(1,1:NRMAX,2)=0.D0
      PECU(1,1:NRMAX  )=0.D0
      POHU(1,1:NRMAX  )=0.D0
      WROTU(1,1:NRMAX )=0.D0
      KFID='QICRHE'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      forall(NR=1:NRMAX,FAS(NR) < 0.d0) FAS(NR) = 0.D0
      PICU(1,1:NRMAX,1) = FAS(1:NRMAX)

      KFID='QICRHI'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      forall(NR=1:NRMAX,FAS(NR) < 0.d0) FAS(NR) = 0.D0
      PICU(1,1:NRMAX,2) = FAS(1:NRMAX)

      KFID='QECH'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      forall(NR=1:NRMAX,FAS(NR) < 0.d0) FAS(NR) = 0.D0
      PECU(1,1:NRMAX) = FAS(1:NRMAX)

      KFID='QECHE'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      forall(NR=1:NRMAX,FAS(NR) < 0.d0) FAS(NR) = 0.D0
      PECU(1,1:NRMAX) = FAS(1:NRMAX)

      KFID='QECHI'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      forall(NR=1:NRMAX,FAS(NR) < 0.d0) FAS(NR) = 0.D0
      PICU(1,1:NRMAX,2) = PICU(1,1:NRMAX,2) + FAS(1:NRMAX)

      KFID='QRAD'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      PRLU(1,1:NRMAX)=FAS(1:NRMAX)

      KFID='QOHM'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      POHU(1,1:NRMAX)=FAS(1:NRMAX)

      AMP=1.D-20
      KFID='SNBIE'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      SNBU(1,1:NRMAX,1)=FAS(1:NRMAX)

      KFID='SNBII'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      SNBU(1,1:NRMAX,2)=FAS(1:NRMAX)

      KFID='SWALL'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      SWLU(1,1:NRMAX)=FAS(1:NRMAX)

      AMP=1.D0
      KFID='VROT'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,1,1,MDLXP,IERR)
      WROTU(1,1:NRMAX)=FAS(1:NRMAX)

!     *****

      MDSUM=MDRGEO+MDAMIN+MDIP+MDBT+MDKAPPA+MDPHIA
      IF(MDSUM.NE.0) THEN
         IF(MDRGEO .NE.0) RRU  (1:NTXMAX)=RR
         IF(MDAMIN .NE.0) RAU  (1:NTXMAX)=RA
         IF(MDIP   .NE.0) RIPU (1:NTXMAX)=RIPS
         IF(MDBT   .NE.0) BBU  (1:NTXMAX)=BB
         IF(MDKAPPA.NE.0) RKAPU(1:NTXMAX)=RKAP
         IF(MDPHIA .NE.0) PHIAU(1:NTXMAX)=0.D0
      ENDIF

!     *** in case of discharge with only one time slice ****
!     *     obtain time point designated by 2d ufiles      *
!     *     from time data array attained by 1d ufiles     *
!     *** in case of discharge with mulitple time slices ***
!     *     obtain time point designated by user           *
!     *     from time data array attained by 2d ufiles     *
!     ******************************************************

      IF(TMU(2).EQ.0.D0) THEN
         DT_MIN=ABS(TMU1(1)-TMU(1))
         DO NTX1=2,NTXMAX1
            IF(DT_MIN.GT.ABS(TMU1(NTX1)-TMU(1))) THEN
               DT_MIN=ABS(TMU1(NTX1)-TMU(1))
               NTS=NTX1
            ENDIF
         ENDDO
         IF(NTS.EQ.0) NTS=1
      ELSE
         DT_MIN=ABS(TMU(1)-TIME_INT)
         DO NTX=2,NTUM
            IF(DT_MIN.GT.ABS(TMU(NTX)-TIME_INT)) THEN
               DT_MIN=ABS(TMU(NTX)-TIME_INT)
               NTS=NTX
            ENDIF
         ENDDO
         IF(NTS.EQ.0) NTS=1
      ENDIF

!     *** GEOMETRY FACTORS ***

      KFID='RMAJOR'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,1,1,MDLXP,IERR)
      RMJRHOU(1,1:NRMAX)=FAS(1:NRMAX)
!      ARRHOU(1,1:NRMAX)=1.D0/RMJRHOU(1,1:NRMAX)**2
      ARRHOU(1,1:NRMAX)=1.D0/RRU(NTS)**2
      TTRHOU(1,1:NRMAX)=BBU(NTS)*RRU(NTS)

      KFID='RMINOR'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,1,0,MDLXP,IERR)
      RMNRHOU(1,1:NRMAX)=FAS(1:NRMAX)

      KFID='KAPPAR'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      IF(IERR.EQ.0) THEN
         RKPRHOU(1,1:NRMAX)=FAS(1:NRMAX)
      ELSE
         RKPRHOU(1,1:NRMAX)=RKAPU(NTS)
      ENDIF

      KFID='GRHO1'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      AR1RHOU(1,1:NRMAX)=FAS(1:NRMAX)

      KFID='GRHO2'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,1,MDLXP,IERR)
      AR2RHOU(1,1:NRMAX)=FAS(1:NRMAX)
      ABRHOU(1,1:NRMAX)=AR2RHOU(1,1:NRMAX)*ARRHOU(1,1:NRMAX)

      KFID='SURF'
      CALL UF2DS(KFID,KUFDEV,KUFDCG,DR,TMU,FAS,AMP,NRMAX,0,0,MDLXP,IERR)
      DVRHOU(1,1:NRMAX)=FAS(1:NRMAX)/AR1RHOU(1,1:NRMAX)

      KFID='VOLUME'
      CALL UFREAD2_TIME(KDIRX,KUFDEV,KUFDCG,KFID,RL,TL,F2,NRLMAX,NTXMAX,NRUM,NTUM,MDCHK,IERR)
      VOLAU(1:NTXMAX)=F2(1:NTXMAX,NRLMAX)

!     *****

      RR    = RRU(NTS)
      RA    = RAU(NTS)
      RIPS  = RIPU(NTS)
      RIPE  = RIPU(NTS)
      BB    = BBU(NTS)
      RKAP  = RKAPU(NTS)
      PHIA  = PHIAU(NTS)

!$$$      IF(NTXMAX1.EQ.0) THEN
!$$$         NTXMAX1=NTXMAX
!$$$         TMU1(1:NTXMAX)=TMU(1:NTXMAX)
!$$$      ENDIF
!$$$C
!$$$      DO NTX=1,NTXMAX1
!$$$         IF(ABS(TMU1(NTX)-TSLC).LE.1.D-5) THEN
!$$$            RR   = RRU(NTX)
!$$$            RA   = RAU(NTX)
!$$$            RIPS = RIPU(NTX)
!$$$            BB   = BBU(NTX)
!$$$            RKAP = RKAPU(NTX)
!$$$            PHIA = PHIAU(NTX)
!$$$            GOTO 2000
!$$$         ENDIF
!$$$      ENDDO
!$$$      WRITE(6,*) 'XX 1D UFILES NO SLICE TIME THE SAME AS ',
!$$$     &           'THAT OF 2D UFILES!'
!$$$      WRITE(6,*) '## THE PROCESS CONTINUES THROUGH USING THE TIMESPL ',
!$$$     &           'INTERPOLATION...'
!$$$      CALL TIMESPL(TSLC,RR  ,TMU1,RRU  ,NTXMAX1,NTUM,IERR)
!$$$      CALL TIMESPL(TSLC,RA  ,TMU1,RAU  ,NTXMAX1,NTUM,IERR)
!$$$      CALL TIMESPL(TSLC,RIPS,TMU1,RIPU ,NTXMAX1,NTUM,IERR)
!$$$      CALL TIMESPL(TSLC,BB  ,TMU1,BBU  ,NTXMAX1,NTUM,IERR)
!$$$      CALL TIMESPL(TSLC,RKAP,TMU1,RKAPU,NTXMAX1,NTUM,IERR)
!$$$      CALL TIMESPL(TSLC,PHIA,TMU1,PHIAU,NTXMAX1,NTUM,IERR)
!$$$ 2000 CONTINUE

      RETURN
      END SUBROUTINE TR_STEADY_UFILE

!     ***********************************************************

!           READING DATA FROM UFILES WITH NT INCREMENT

!     ***********************************************************

      SUBROUTINE TR_UFREAD

      USE TRCOMM, ONLY : &
           NTUM, DT, NT, RKAP, RKAPU, RR, RRU, RA, RAU, BB, BBU, PHIA, &
           PHIAU, PNBIU, VOLAU, RMJRHOU, NRMAX, RHOA, NROMAX, MDLUF, MDLEQN, &
           RNU, RN, MDNI, MDLEOI, RTU, RT, QPU, QP, KUFDEV, PEX, PNBU, PICU, &
           PECU, PRF, WROTU, TTRHOU, DVRHOU, ABRHOU, ARRHOU, AR1RHOU, &
           AR2RHOU, RMNRHOU, RKPRHOU, WROT, VTOR, TTRHO, DVRHO, ABRHO, &
           ARRHO, AR1RHO, AR2RHO, RMJRHO, RMNRHO, RKPRHO, MDPHIA, PI, RJCB, &
           RHOM, RM, RHOG, RG, PBMU, RNFU, PBM, RNF, MDLEQB, MDLJQ, AJU, &
           AJNBU, AJ, AJNB, RMU0, DVRHOG, TTRHOG, RDP, DR, BP, &
           AR1RHOG, ARRHOG, AJTOR, NRAMAX, RIPA, PNSU, PTSU, PNSUA, PTSUA, &
           PNS, PTS, PNSA, PTSA, ALP, PROFT1, PROFT2, ANC, PNC, ANFE, PNFE, &
           NSM, PZ, PZFE, PZC, PNSS, PNSSA, ABVRHOG, RDPVRHOG, rkind
      USE TRCOM1, ONLY : INS, NTAMAX, NTXMAX, NTXMAX1, PNBI, TMU, TMU1
      USE trmetric
      USE libitp
      IMPLICIT NONE
      INTEGER:: IERR, NR, NS, NTL, NTLL
      REAL(rkind)   :: &
           ABRL, AJL, AJNBL, ANEAVE, ANI, ANZ, AR1RL, AR2RL, ARRL, DILUTE, &
           DVRL, FACTOR0, FACTORM, FACTORP, PBML, PECL, PICDL, PICEL, PNBDL, &
           PNBEL, PNSAL, PNSL, PROF, PTSAL, PTSL, QPL, RHO_A, RKAPL, RMJRL, &
           RMJRXL, RMNRL, RNDL, RNEL, RNFL, RNIL, RTDL, RTDM, RTEL, RTIL, &
           TSL, TTRL, VOL, VOLM, WROTL

      REAL(rkind),DIMENSION(NRMAX):: AJTMP, DSRHO


      TSL=DT*DBLE(NT)
      CALL TIMESPL(TSL,RKAP,TMU1,RKAPU,NTXMAX1,NTUM,IERR)
      CALL TIMESPL(TSL,RR  ,TMU1,RRU  ,NTXMAX1,NTUM,IERR)
      CALL TIMESPL(TSL,RA  ,TMU1,RAU  ,NTXMAX1,NTUM,IERR)
      CALL TIMESPL(TSL,BB  ,TMU1,BBU  ,NTXMAX1,NTUM,IERR)
      CALL TIMESPL(TSL,PHIA,TMU1,PHIAU,NTXMAX1,NTUM,IERR)
      CALL TIMESPL(TSL,PNBI,TMU1,PNBIU,NTXMAX1,NTUM,IERR)

      CALL TIMESPL(TSL,VOLM,TMU,VOLAU,NTXMAX,NTUM,IERR)
      CALL TIMESPL(TSL,RMJRXL,TMU,RMJRHOU(1:NTUM,NRMAX),NTXMAX,NTUM,IERR)

      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      select case(MDLUF)
      case(1)
      DO NR=1,NRMAX
         IF(MDLEQN.EQ.0) THEN
            CALL TIMESPL(TSL,RNEL ,TMU,RNU(1:NTUM,NR,1),NTXMAX,NTUM,IERR)
            CALL TIMESPL(TSL,RNDL ,TMU,RNU(1:NTUM,NR,2),NTXMAX,NTUM,IERR)
            RN(NR,1)=RNEL
            RN(NR,2)=RNDL
            IF(MDNI.NE.0) THEN
               CALL TIMESPL(TSL,RNIL ,TMU,RNU(1:NTUM,NR,3),NTXMAX,NTUM,IERR)
               RN(NR,3)=RNIL
            ENDIF
         ENDIF
         IF(INS.NE.0) THEN
            IF(MDLEOI.EQ.1) THEN
               CALL TIMESPL(TSL,RTDL ,TMU,RTU(1:NTUM,NR,2),NTXMAX,NTUM,IERR)
               RT(NR,2)=RTDL
            ELSEIF(MDLEOI.EQ.2) THEN
               CALL TIMESPL(TSL,RTEL ,TMU,RTU(1:NTUM,NR,1),NTXMAX,NTUM,IERR)
               RT(NR,1)=RTEL
            ENDIF
            CALL TIMESPL(TSL,RTIL ,TMU,RTU(1:NTUM,NR,3),NTXMAX,NTUM,IERR)
            IF(INS.EQ.2) RT(NR,3)=RTIL
         ENDIF
         CALL TIMESPL(TSL,QPL ,TMU,QPU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         QP(NR)=QPL

         IF(KUFDEV.EQ.'X') THEN
            IF(PNBI.LT.12.D6) THEN
               PEX(NR,1)=PNBU(1,NR,1)
               PEX(NR,2)=PNBU(1,NR,2)
            ELSE
               PEX(NR,1)=PNBU(2,NR,1)
               PEX(NR,2)=PNBU(2,NR,2)
               IF(NT.EQ.NTAMAX) THEN
                  PEX(NR,1)=PNBU(3,NR,1)
                  PEX(NR,2)=PNBU(3,NR,2)
               ENDIF
            ENDIF
         ELSE
            CALL TIMESPL(TSL,PNBEL,TMU,PNBU(1:NTUM,NR,1),NTXMAX,NTUM,IERR)
            CALL TIMESPL(TSL,PNBDL,TMU,PNBU(1:NTUM,NR,2),NTXMAX,NTUM,IERR)
            IF(PNBEL.LT.0.D0) PNBEL=0.D0
            IF(PNBDL.LT.0.D0) PNBDL=0.D0
            PEX(NR,1)=PNBEL
            PEX(NR,2)=PNBDL
         ENDIF

         CALL TIMESPL(TSL,PICEL,TMU,PICU(1:NTUM,NR,1),NTXMAX,NTUM,IERR)
         IF(PICEL.LT.0.D0) PICEL=0.D0
         CALL TIMESPL(TSL,PICDL,TMU,PICU(1:NTUM,NR,2),NTXMAX,NTUM,IERR)
         IF(PICDL.LT.0.D0) PICDL=0.D0
         CALL TIMESPL(TSL,PECL ,TMU,PECU(1:NTUM,NR  ),NTXMAX,NTUM,IERR)
         IF(PECL.LT.0.D0) PECL=0.D0
         PRF(NR,1)=PICEL+PECL
         PRF(NR,2)=PICDL
         CALL TIMESPL(TSL,WROTL,TMU,WROTU  (1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,TTRL ,TMU,TTRHOU (1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,DVRL ,TMU,DVRHOU (1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,ABRL ,TMU,ABRHOU (1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,ARRL ,TMU,ARRHOU (1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,AR1RL,TMU,AR1RHOU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,AR2RL,TMU,AR2RHOU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,RMJRL,TMU,RMJRHOU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,RMNRL,TMU,RMNRHOU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,RKAPL,TMU,RKPRHOU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         WROT(NR)=WROTL
         VTOR(NR)=WROTL*RMJRL
         TTRHO(NR)=TTRL
         DVRHO(NR)=DVRL
         ABRHO(NR)=ABRL
         ARRHO(NR)=ARRL
         AR1RHO(NR)=AR1RL
         AR2RHO(NR)=AR2RL
         RMJRHO(NR)=RMJRL
         RMNRHO(NR)=RMNRL
         RKPRHO(NR)=RKAPL
         IF(MDPHIA.EQ.0) THEN
            RHO_A=SQRT(PHIA/(PI*BB))
            RJCB(NR)=1.D0/RHO_A
            RHOM(NR)=RM(NR)*RHO_A
            RHOG(NR)=RG(NR)*RHO_A
         ELSE
            RHO_A=SQRT(VOLM/(2.D0*PI**2*RMJRXL))
            RJCB(NR)=1.D0/RHO_A
            RHOM(NR)=RM(NR)/RJCB(NR)
            RHOG(NR)=RG(NR)/RJCB(NR)
         ENDIF
         CALL TIMESPL(TSL,PBML,TMU,PBMU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,RNFL,TMU,RNFU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         PBM(NR)=PBML
         RNF(NR,1)=RNFL
      ENDDO
      CALL TRGFRG
      select case(MDLEQB)
      case(0)
         IF(MDLJQ.EQ.0) THEN ! *** MDLJQ ***
            NR=1
               CALL TIMESPL(TSL,AJL  ,TMU,AJU  (1:NTUM,NR),NTXMAX,NTUM,IERR)
               CALL TIMESPL(TSL,AJNBL,TMU,AJNBU(1:NTUM,NR),NTXMAX,NTUM,IERR)
               AJ(NR)=AJL
               AJNB(NR)=AJNBL
               FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
               FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
               RDPVRHOG(NR)=FACTOR0*DR/FACTORP
               RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
               BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
               QP(NR)=TTRHOG(NR)*ARRHOG(NR)/(4.D0*PI**2*RDPVRHOG(NR))
            DO NR=2,NRMAX
               CALL TIMESPL(TSL,AJL  ,TMU,AJU  (1:NTUM,NR),NTXMAX,NTUM,IERR)
               CALL TIMESPL(TSL,AJNBL,TMU,AJNBU(1:NTUM,NR),NTXMAX,NTUM,IERR)
               AJ(NR)=AJL
               AJNB(NR)=AJNBL
               FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
               FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
               FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
               RDPVRHOG(NR)=(FACTOR0*DR+FACTORM*RDPVRHOG(NR-1))/FACTORP
               RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
               BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
               QP(NR)=TTRHOG(NR)*ARRHOG(NR)/(4.D0*PI**2*RDPVRHOG(NR))
            ENDDO
            NR=1
               FACTOR0=RR/(RMU0*DVRHO(NR))
               FACTORP=ABVRHOG(NR)
               AJTOR(NR)=FACTOR0*FACTORP*RDPVRHOG(NR)/DR
            DO NR=2,NRMAX
               FACTOR0=RR/(RMU0*DVRHO(NR))
               FACTORM=ABVRHOG(NR-1)
               FACTORP=ABVRHOG(NR  )
               AJTOR(NR)=FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
            ENDDO
         ELSE ! *** MDLJQ ***
            DO NR=1,NRMAX
               RDPVRHOG(NR)=TTRHOG(NR)*ARRHOG(NR)/(4.D0*PI**2*QP(NR))
               RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
               BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
            ENDDO
         ENDIF ! *** MDLJQ ***
      case default
!        boundary condition for polidal flux at rhoa defined by exp. data
         VOL=0.D0
         DO NR=1,NRAMAX
            CALL TIMESPL(TSL,AJL ,TMU,AJU  (1:NTUM,NR),NTXMAX,NTUM,IERR)
            VOL=VOL+DVRHO(NR)*DR
            DSRHO(NR)=DVRHO(NR)/(2.D0*PI*RR)
            AJTMP(NR)=AJL
         ENDDO
         RIPA=SUM(AJTMP(1:NRMAX)*DSRHO(1:NRMAX))*DR/1.D6

         DO NR=1,NRMAX
            CALL TIMESPL(TSL,AJNBL,TMU,AJNBU(1:NTUM,NR),NTXMAX,NTUM,IERR)
            AJNB(NR)=AJNBL
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
            QP(NR)=TTRHOG(NR)*ARRHOG(NR)/(4.D0*PI**2*RDPVRHOG(NR))
         ENDDO
      end select
      case(3)
      DO NR=1,NRMAX
         CALL TIMESPL(TSL,RNEL ,TMU,RNU(1:NTUM,NR,1),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,RNDL ,TMU,RNU(1:NTUM,NR,2),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,RNIL ,TMU,RNU(1:NTUM,NR,3),NTXMAX,NTUM,IERR)
         RN(NR,1)=RNEL
         RN(NR,2)=RNDL
         RN(NR,3)=RNIL
!         CALL TIMESPL(TSL,QPL ,TMU,QPU(1,NR),NTXMAX,NTUM,IERR)
!         QP(NR)=QPL
         CALL TIMESPL(TSL,PNBEL,TMU,PNBU(1:NTUM,NR,1),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,PNBDL,TMU,PNBU(1:NTUM,NR,2),NTXMAX,NTUM,IERR)
!         IF(PNBEL.LT.0.D0) PNBEL=0.D0
!         IF(PNBDL.LT.0.D0) PNBDL=0.D0
         PEX(NR,1)=PNBEL
         PEX(NR,2)=PNBDL
         CALL TIMESPL(TSL,PICEL,TMU,PICU(1:NTUM,NR,1),NTXMAX,NTUM,IERR)
!         IF(PICEL.LT.0.D0) PICEL=0.D0
         CALL TIMESPL(TSL,PICDL,TMU,PICU(1:NTUM,NR,2),NTXMAX,NTUM,IERR)
!         IF(PICDL.LT.0.D0) PICDL=0.D0
         CALL TIMESPL(TSL,PECL ,TMU,PECU(1:NTUM,NR  ),NTXMAX,NTUM,IERR)
!         IF(PECL.LT.0.D0) PECL=0.D0
         IF(TSL.LE.1.D0) THEN
            PRF(NR,1)=0.D0
            PRF(NR,2)=0.D0
         ELSE
            DO NTL=1,NTXMAX
               IF(TMU(NTL).GT.TSL) THEN
                  NTLL=NTL
                  GOTO 1000
               ENDIF
            ENDDO
            NTLL=NTXMAX
 1000       PRF(NR,1)=PICU(NTLL,NR,1)+PECU(NTLL,NR)
            PRF(NR,2)=PICU(NTLL,NR,2)
         ENDIF
         CALL TIMESPL(TSL,TTRL ,TMU,TTRHOU (1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,DVRL ,TMU,DVRHOU (1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,ABRL ,TMU,ABRHOU (1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,ARRL ,TMU,ARRHOU (1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,AR1RL,TMU,AR1RHOU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,AR2RL,TMU,AR2RHOU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,RMJRL,TMU,RMJRHOU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,RMNRL,TMU,RMNRHOU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         CALL TIMESPL(TSL,RKAPL,TMU,RKPRHOU(1:NTUM,NR),NTXMAX,NTUM,IERR)
         TTRHO(NR)=TTRL
         DVRHO(NR)=DVRL
         ABRHO(NR)=ABRL
         ARRHO(NR)=ARRL
         AR1RHO(NR)=AR1RL
         AR2RHO(NR)=AR2RL
         RMJRHO(NR)=RMJRL
         RMNRHO(NR)=RMNRL
         RKPRHO(NR)=RKAPL
         IF(MDPHIA.EQ.0) THEN
            RHO_A=SQRT(PHIA/(PI*BB))
            RJCB(NR)=1.D0/RHO_A
            RHOM(NR)=RM(NR)*RHO_A
            RHOG(NR)=RG(NR)*RHO_A
         ELSE
            RHO_A=SQRT(VOLM/(2.D0*PI**2*RMJRXL))
            RJCB(NR)=1.D0/RHO_A
            RHOM(NR)=RM(NR)/RJCB(NR)
            RHOG(NR)=RG(NR)/RJCB(NR)
         ENDIF
      ENDDO
      end select
      CALL TRGFRG

      IF(MDLUF.EQ.1) THEN
         DO NS=1,2
            CALL TIMESPL(TSL,PNSL ,TMU,PNSU (1:NTUM,NS),NTXMAX,NTUM,IERR)
            CALL TIMESPL(TSL,PTSL ,TMU,PTSU (1:NTUM,NS),NTXMAX,NTUM,IERR)
            CALL TIMESPL(TSL,PNSAL,TMU,PNSUA(1:NTUM,NS),NTXMAX,NTUM,IERR)
            CALL TIMESPL(TSL,PTSAL,TMU,PTSUA(1:NTUM,NS),NTXMAX,NTUM,IERR)
            PNS (NS)=PNSL
            PTS (NS)=PTSL
            PNSA(NS)=PNSAL
            PTSA(NS)=PTSAL
         ENDDO
         IF(MDNI.NE.0) THEN
            NS=3
            CALL TIMESPL(TSL,PNSL ,TMU,PNSU (1:NTUM,NS),NTXMAX,NTUM,IERR)
            CALL TIMESPL(TSL,PTSL ,TMU,PTSU (1:NTUM,NS),NTXMAX,NTUM,IERR)
            CALL TIMESPL(TSL,PNSAL,TMU,PNSUA(1:NTUM,NS),NTXMAX,NTUM,IERR)
            CALL TIMESPL(TSL,PTSAL,TMU,PTSUA(1:NTUM,NS),NTXMAX,NTUM,IERR)
            PNS (NS)=PNSL
            PTS (NS)=PTSL
            PNSA(NS)=PNSAL
            PTSA(NS)=PTSAL
         ENDIF
         IF(RHOA.NE.1.D0) THEN
            CALL TIMESPL(TSL,RTDM,TMU,RTU(1:NTUM,NRMAX,2),NTXMAX,NTUM,IERR)
            DO NR=NRAMAX+1,NRMAX
               PROF    =(1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               CALL TIMESPL(TSL,RTEL,TMU,RTU(1:NTUM,NR,1),NTXMAX,NTUM,IERR)
               CALL TIMESPL(TSL,RTDL,TMU,RTU(1:NTUM,NR,2),NTXMAX,NTUM,IERR)
               RT(NR,1)= RTEL
               RT(NR,2)= RTDL
               RT(NR,3)= RTDL
               RT(NR,4)=(RTDL-RTDM)*PROF+RTDM
            ENDDO
            IF(MDNI.EQ.0) THEN
               DO NR=NRAMAX+1,NRMAX
                  CALL TIMESPL(TSL,RTDL,TMU,RTU(1:NTUM,NR,2),NTXMAX,NTUM,IERR)
                  PROF    =(1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
                  RT(NR,3)=(RTDL-RTDM)*PROF+RTDM
               ENDDO
            ENDIF
         ENDIF
      ENDIF

!     *** CALCULATE PZC,PZFE ***

      CALL TRZEFF

!     *** CALCULATE ANEAVE ***

      ANEAVE=SUM(RN(1:NRMAX,1)*RM(1:NRMAX))*2.D0*DR

!     *** CALCULATE IMPURITY DENSITY
!                ACCORDING TO ITER PHYSICS DESIGN GUIDELINE ***

      IF(MDLUF.NE.3) THEN
         DO NR=1,NRMAX
            ANC (NR)= (0.9D0+0.60D0*(0.7D0/ANEAVE)**2.6D0)*PNC*1.D-2*RN(NR,1)
            ANFE(NR)= (0.0D0+0.05D0*(0.7D0/ANEAVE)**2.3D0)*PNFE*1.D-2*RN(NR,1)
            ANI = SUM(PZ(2:NSM)*RN(NR,2:NSM))
            ANZ = PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
            DILUTE = 1.D0-ANZ/ANI
            RN(NR,2:NSM) = RN(NR,2:NSM)*DILUTE
         ENDDO
         PNSS(1)=PNS(1)
         PNSS(2:NSM)=PNS(2:NSM)*DILUTE
         PNSS(7)=PNS(7)
         PNSS(8)=PNS(8)
         IF(RHOA.NE.1.D0) THEN
            PNSSA(1)=PNSA(1)
            DO NS=2,NSM
               PNSSA(NS)=PNSA(NS)*DILUTE
            ENDDO
            PNSSA(7)=PNSA(7)
            PNSSA(8)=PNSA(8)
         ENDIF
      ELSE
!$$$         ANC(1:NRMAX)=1.D0/3.D1*RN(1:NRMAX,1)
         PNSS(1:NSM)=PNS(1:NSM)
         PNSS(7)=PNS(7)
         PNSS(8)=PNS(8)
         IF(RHOA.NE.1.D0) THEN
!            PNSSA(1:NSM)=PNSA(1:NSM)*DILUTE      ! DILUTE is not defined
            PNSSA(1:NSM)=PNSA(1:NSM)
            PNSSA(7)=PNSA(7)
            PNSSA(8)=PNSA(8)
         ENDIF
      ENDIF

      IF(RHOA.NE.1.D0) NRMAX=NRAMAX

      RETURN
      END SUBROUTINE TR_UFREAD

!     ******

      SUBROUTINE TR_UFREAD_S

      USE TRCOMM, ONLY : ABRHO, ABRHOU, ANC, ANFE, AR1RHO, AR1RHOG, AR1RHOU, AR2RHO, AR2RHOU, ARRHO, ARRHOG, ARRHOU, BB, BBU, &
     &                   BP, DR, DVRHO, DVRHOG, DVRHOU, MDLEOI, MDLEQN, MDNI, MDPHIA, NRAMAX, NRMAX, NROMAX, NSM, NTS, PECU,  &
     &                   PEX, PHIA, PHIAU, PI, PICU, PNBU, PNC, PNFE, PNS, PNSA, PNSS, PNSSA, PRF, PZ, PZC, PZFE, QP, RA, RAU,&
     &                   RDP, RG, RHOA, RHOG, RHOM, RJCB, RKAP, RKAPU, RKPRHO, RKPRHOU, RM, RMJRHO, RMJRHOU, RMNRHO, RMNRHOU, &
     &                   RN, RNF, RNFU, RNU, RR, RRU, RT, RTU, TTRHO, TTRHOG, TTRHOU, VOLAU, VTOR, WROT, WROTU, RDPVRHOG, rkind
      USE TRCOM1, ONLY : INS
      USE trmetric
      IMPLICIT NONE
      INTEGER:: NR
      REAL(rkind)   :: ANEAVE, ANI, ANZ, DILUTE, RHO_A, RMJRXL, VOLM


      RKAP=RKAPU(NTS)
      RR=RRU(NTS)
      RA=RAU(NTS)
      BB=BBU(NTS)
      PHIA=PHIAU(NTS)

      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      VOLM=VOLAU(1)
      RMJRXL=RMJRHOU(1,NRMAX)
      DO NR=1,NRMAX
         IF(MDLEQN.EQ.0) THEN
            RN(NR,1)=RNU(1,NR,1)
            RN(NR,2)=RNU(1,NR,2)
            IF(MDNI.NE.0) RN(NR,3)=RNU(1,NR,3)
         ENDIF
         IF(INS.NE.0) THEN
            IF(MDLEOI.EQ.1) THEN
               RT(NR,2)=RTU(1,NR,2)
            ELSEIF(MDLEOI.EQ.2) THEN
               RT(NR,1)=RTU(1,NR,1)
            ENDIF
            IF(INS.EQ.2) RT(NR,3)=RTU(1,NR,3)
         ENDIF
         PEX(NR,1)=PNBU(1,NR,1)
         PEX(NR,2)=PNBU(1,NR,2)
         PRF(NR,1)=PICU(1,NR,1)+PECU(1,NR)
         PRF(NR,2)=PICU(1,NR,2)
         WROT(NR) =WROTU(1,NR)
         VTOR(NR) =WROTU(1,NR)*RMJRHOU(1,NR)
         TTRHO(NR)=TTRHOU(1,NR)
         DVRHO(NR)=DVRHOU(1,NR)
         ABRHO(NR)=ABRHOU(1,NR)
         ARRHO(NR)=ARRHOU(1,NR)
         AR1RHO(NR)=AR1RHOU(1,NR)
         AR2RHO(NR)=AR2RHOU(1,NR)
         RMJRHO(NR)=RMJRHOU(1,NR)
         RMNRHO(NR)=RMNRHOU(1,NR)
         RKPRHO(NR)=RKPRHOU(1,NR)
         IF(MDPHIA.EQ.0) THEN
            RHO_A=SQRT(PHIA/(PI*BB))
            RJCB(NR)=1.D0/RHO_A
            RHOM(NR)=RM(NR)*RHO_A
            RHOG(NR)=RG(NR)*RHO_A
         ELSE
            RHO_A=SQRT(VOLM/(2.D0*PI**2*RMJRXL))
            RJCB(NR)=1.D0/RHO_A
            RHOM(NR)=RM(NR)/RJCB(NR)
            RHOG(NR)=RG(NR)/RJCB(NR)
         ENDIF
         RNF(NR,1)=RNFU(1,NR)
      ENDDO
      CALL TRGFRG
      DO NR=1,NRMAX
         RDPVRHOG(NR)=TTRHOG(NR)*ARRHOG(NR)/(4.D0*PI**2*QP(NR))
         RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
         BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
      ENDDO

!     *** CALCULATE PZC,PZFE ***

      CALL TRZEFF

!     *** CALCULATE ANEAVE ***

      ANEAVE=SUM(RN(1:NRMAX,1)*RM(1:NRMAX))*2.D0*DR

!     *** CALCULATE IMPURITY DENSITY
!                ACCORDING TO ITER PHYSICS DESIGN GUIDELINE ***

      DO NR=1,NRMAX
         ANC (NR)= (0.9D0+0.60D0*(0.7D0/ANEAVE)**2.6D0)*PNC*1.D-2*RN(NR,1)
         ANFE(NR)= (0.0D0+0.05D0*(0.7D0/ANEAVE)**2.3D0)*PNFE*1.D-2*RN(NR,1)
         ANI = SUM(PZ(2:NSM)*RN(NR,2:NSM))
         ANZ = PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
         DILUTE = 1.D0-ANZ/ANI
         RN(NR,2:NSM) = RN(NR,2:NSM)*DILUTE
      ENDDO
      PNSS(1)=PNS(1)
      PNSS(2:NSM)=PNS(2:NSM)*DILUTE
      PNSS(7)=PNS(7)
      PNSS(8)=PNS(8)
      IF(RHOA.NE.1.D0) THEN
         PNSSA(1)=PNSA(1)
         PNSSA(2:NSM)=PNSA(2:NSM)*DILUTE
         PNSSA(7)=PNSA(7)
         PNSSA(8)=PNSA(8)
      ENDIF

      IF(RHOA.NE.1.D0) NRMAX=NRAMAX

      RETURN
      END SUBROUTINE TR_UFREAD_S
