MODULE dpmenu

PRIVATE

PUBLIC dp_menu

CONTAINS

!     ***** TASK/DP MENU *****

  SUBROUTINE dp_menu

    USE dpcomm
    USE plprof
    USE plparm,ONLY: pl_view
    USE dpparm,ONLY: dp_parm,dp_view
    USE dproot,ONLY: dp_root,dpgrp1,dpgrp0
    USE dpcont,ONLY: dp_cont,dp_contx
    USE dpcont2,ONLY: dp_cont2
    USE dptens,ONLY: dp_tens
    IMPLICIT NONE
      
    CHARACTER(LEN=1):: KID
    CHARACTER(LEN=80):: LINE
    INtEGER:: MODE,IERR,NID
    COMPLEX(rkind):: CD4(6),CD5(6),CD6(6)
    INTEGER:: MODELP_SAVE,I
    REAL(rkind):: RF0_SAVE,RKZ0_SAVE,RKX0_SAVE,RL,RHON
    COMPLEX(rkind):: CW,CKPR,CKPP
    TYPE(pl_mag_type):: mag
    TYPE(pl_plf_type),DIMENSION(nsmax):: plf
    TYPE(pl_grd_type),DIMENSION(nsmax):: grd

1   CONTINUE
    WRITE(6,*) '## DP MENU: P,V/PARM  ', &
               'D0,D1,D2,D3,D4/DISP  F/ROOT  T,S,K/TEST  Q/QUIT'

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
       IF(NID.EQ.0) THEN
          CALL DPGRP0
       ELSEIF(NID.EQ.1) THEN
          CALL DPGRP1
       ELSEIF(NID.EQ.2) THEN
          CALL DP_CONT
       ELSEIF(NID.EQ.3) THEN
          CALL DP_CONTX
       ELSEIF(NID.EQ.4) THEN
          CALL DP_CONT2
       ELSE
          WRITE(6,*) 'XX DPMENU: unknown NID'
       ENDIF
    ELSEIF(KID.EQ.'F') THEN
       CALL DP_ROOT
    ELSEIF(KID.EQ.'T') THEN
       MODELP_SAVE=MODELP(1)
       RF0_SAVE=RF0
       RKZ0_SAVE=RKZ0
       RKX0_SAVE=RKX0
       RL=RR
1001   WRITE(6,*) '# INPUT: RL,RF0,RKZ0,RKX0 ='
       READ(5,*,ERR=1001,END=1002) RL,RF0,RKZ0,RKX0
       IF(RF0.LE.0.D0) GOTO 1002
       CW=2.D0*PI*DCMPLX(RF0,RFI0)*1.D6
       CKPR=MAX(RKZ0,1.D-4)
       CKPP=RKX0
       CALL PL_MAG(RL,0.D0,0.D0,mag)
       RHON=mag%rhon
       CALL PL_PROF(RHON,plf)
       CALL PL_GRAD(RHON,grd)
       MODELP(1)=1
       CALL DP_TENS(CW,CKPR,CKPP,1,mag,plf,grd,CD4)
       MODELP(1)=4
       CALL DP_TENS(CW,CKPR,CKPP,1,mag,plf,grd,CD5)
       MODELP(1)=5
       CALL DP_TENS(CW,CKPR,CKPP,1,mag,plf,grd,CD6)
       WRITE(6,602) 
602    FORMAT(8X,'MODELP=1',16X,'MODELP=4',16X,'MODELP=5')
       WRITE(6,603) (I,CD4(I),CD5(I),CD6(I),I=1,6)
603    FORMAT(I5,1P6E12.4)
       GOTO 1001

1002   MODELP(1)=MODELP_SAVE
       RF0=RF0_SAVE
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
