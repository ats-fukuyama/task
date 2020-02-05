MODULE dpmenu

PRIVATE

PUBLIC dp_menu

CONTAINS

!     ***** TASK/DP MENU *****

  SUBROUTINE dp_menu

    USE dpcomm_local
    USE plprof
    USE plprofw
    USE plparm,ONLY: pl_view
    USE dpparm,ONLY: dp_parm,dp_view
    USE dproot,ONLY: dp_root,dpgrp1,dpgrp0
    USE dpcont,ONLY: dp_cont2,dp_cont3
    USE dpcont4,ONLY: dp_cont4
    USE dptnsr0,ONLY: dp_tnsr0
    IMPLICIT NONE
      
    CHARACTER(LEN=1):: KID
    CHARACTER(LEN=80):: LINE
    INtEGER:: MODE,IERR,NID
    COMPLEX(rkind):: CD1(6),CD2(6),CD3(6),CD4(6)
    INTEGER:: I
    INTEGER:: NSMAX_SAVE,NRMAX_SAVE,MODELV_SAVE,MODELP_SAVE
    REAL(rkind):: RF0_SAVE,RFI0_SAVE,RKZ0_SAVE,RKX0_SAVE,RL,RHON
    COMPLEX(rkind):: CW,CKPR,CKPP
    TYPE(pl_mag_type):: mag
    TYPE(pl_plfw_type),DIMENSION(nsmax):: plfw
    TYPE(pl_grd_type),DIMENSION(nsmax):: grd

1   CONTINUE
    WRITE(6,*) '## DP MENU: P,V/PARM  ', &
               'D0,D1,D2,D3,D4,D5/DISP  F/ROOT  T,S,K/TEST  Q/QUIT'

    CALL TASK_KLIN(LINE,KID,MODE,DP_PARM)
    IF(MODE.NE.1) GOTO 1
    CALL GUCPTL(KID)

    IF(KID.EQ.'P') THEN
       CALL DP_PARM(0,'DP',IERR)
    ELSEIF(KID.EQ.'V') THEN
       CALL PL_VIEW
       CALL DP_VIEW
    ELSEIF(KID.EQ.'D') THEN
       READ(LINE(2:),*,ERR=1,END=1) NID
       SELECT CASE(NID)
       CASE(0)
          CALL DPGRP0
       CASE(1)
          CALL DPGRP1
       CASE(2)
          CALL DP_CONT2
       CASE(3)
          CALL DP_CONT3
       CASE(4,5)
          CALL DP_CONT4(NID)
       CASE DEFAULT
          WRITE(6,*) 'XX DPMENU: unknown NID'
       END SELECT
    ELSEIF(KID.EQ.'F') THEN
       CALL DP_ROOT
    ELSEIF(KID.EQ.'T') THEN
       NSMAX_SAVE=NSMAX
       NRMAX_SAVE=NRMAX_DP
       MODELV_SAVE=MODELV(1)
       MODELP_SAVE=MODELP(1)
       RF0_SAVE=RF0
       RFI0_SAVE=RFI0
       RKZ0_SAVE=RKZ0
       RKX0_SAVE=RKX0
       RL=RR
1001   WRITE(6,*) '# INPUT: RL,RF0,RFI0,RKZ0,RKX0 ='
       WRITE(6,'(10X,1P5E12.4)') RL,RF0,RFI0,RKZ0,RKX0
       READ(5,*,ERR=1001,END=1002) RL,RF0,RFI0,RKZ0,RKX0
       IF(RF0.LE.0.D0) GOTO 1002
       CW=2.D0*PI*DCMPLX(RF0,RFI0)*1.D6
       CKPR=MAX(RKZ0,1.D-8)
       CKPP=RKX0
       CALL PL_MAG(RL,0.D0,0.D0,mag)
       RHON=mag%rhon
       CALL PL_PROFW(RHON,plfw)
       CALL PL_GRAD(RHON,grd)

       NSMAX=1
       NRMAX_DP=1

       MODELV(1)=0
       MODELP(1)=6
       CALL DP_TNSR0(CW,CKPR,CKPP,1,mag,plfw,grd,CD1)
       MODELV(1)=1
       MODELP(1)=6
       CALL DP_TNSR0(CW,CKPR,CKPP,1,mag,plfw,grd,CD2)
       MODELV(1)=0
       MODELP(1)=9
       CALL DP_TNSR0(CW,CKPR,CKPP,1,mag,plfw,grd,CD3)
       MODELV(1)=3
       MODELP(1)=6
       CALL DP_TNSR0(CW,CKPR,CKPP,1,mag,plfw,grd,CD4)

       WRITE(6,602) 
602    FORMAT(6X,'MODELV=0',13X,'MODELV=1', &
             13X,'MODELV=0',13X,'MODELV=3')
       WRITE(6,603) 
603    FORMAT(6X,'MODELP=6',13X,'MODELP=6', &
             13X,'MODELP=9',13X,'MODELP=6')
       WRITE(6,604) (CD1(I),CD2(I),CD3(I),CD4(I),I=1,6)
604    FORMAT((1PE9.2,1P7E10.2))
       GOTO 1001

1002   CONTINUE
       NSMAX=NSMAX_SAVE
       NRMAX_DP=NRMAX_SAVE
       MODELV(1)=MODELV_SAVE
       MODELP(1)=MODELP_SAVE
       RF0=RF0_SAVE
       RFI0=RFI0_SAVE
       RKZ0=RKZ0_SAVE
       RKX0=RKX0_SAVE
       GOTO 1
    ELSEIF(KID.EQ.'K') THEN
       CONTINUE
    ELSEIF(KID.EQ.'S') THEN
       CONTINUE
    ELSEIF(KID.EQ.'Q') THEN
       GOTO 9000
    ENDIF
    GOTO 1

9000 RETURN
  END SUBROUTINE dp_menu

END MODULE dpmenu
