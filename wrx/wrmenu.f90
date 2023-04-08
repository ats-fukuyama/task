! wrmenu.f90

MODULE wrmenu

  PRIVATE
  PUBLIC wr_menu

CONTAINS

!     ***** TASK/WR MENU *****

  SUBROUTINE WR_MENU

    USE wrcomm
    USE plparm,ONLY: pl_view
    USE dpparm,ONLY: dp_view
    USE dproot,ONLY: dp_root,dpgrp1
    USE dpcont,ONLY: dp_cont2,dp_cont3
    USE wrparm,ONLY: wr_parm
    USE wrview,ONLY: wr_view
    USE wrsetup,ONLY: wr_setup
    USE wrexec,ONLY: wr_exec
    USE wrgout,ONLY: wr_gout
    USE wrfile,ONLY: wr_save,wr_load,wr_write
    USE libkio
    IMPLICIT NONE
    CHARACTER(LEN=1):: KID
    CHARACTER(LEN=80):: LINE
    INTEGER:: NSTAT,IERR,MODE,NID

      NSTAT=0

    1 CONTINUE
         IERR=0
         WRITE(6,601)
  601    FORMAT('## WR MENU: P,V/PARM  R,B/RAY  G/GRAPH  S,L,W/FILE', &
                '  Dn/DISP  F/ROOT  Q/QUIT')
         CALL TASK_KLIN(LINE,KID,MODE,WR_PARM)
         IF(MODE.EQ.3) STOP
      IF(MODE.NE.1) GOTO 1

      IF(KID.EQ.'P') THEN
         CALL WR_PARM(0,'WR',IERR)
      ELSEIF(KID.EQ.'V') THEN
         CALL PL_VIEW
         CALL DP_VIEW
         CALL WR_VIEW
      ELSEIF(KID.EQ.'D') THEN
         READ(LINE(2:),*,ERR=1,END=1) NID
         IF(NID.EQ.1) THEN
            CALL DPGRP1
         ELSEIF(NID.EQ.2) THEN
            CALL DP_CONT2
         ELSEIF(NID.EQ.3) THEN
            CALL DP_CONT3
         ELSE
            WRITE(6,*) 'XX WRMENU: unknown NID'
         ENDIF
      ELSEIF(KID.EQ.'F') THEN
         CALL DP_ROOT
      ELSEIF(KID.EQ.'R') THEN
         CALL wr_allocate
         CALL wr_setup(ierr)
         IF(ierr.NE.0) GO TO 1
         CALL wr_exec(nstat,ierr)
      ELSEIF(KID.EQ.'G') THEN
         CALL WR_GOUT(NSTAT)
      ELSEIF(KID.EQ.'S') THEN
         CALL wr_save
      ELSEIF(KID.EQ.'L') THEN
         CALL wr_load(NSTAT)
      ELSEIF(KID.EQ.'W') THEN
         CALL wr_write
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ELSE
         WRITE(6,*) 'XX WRMENU: UNKNOWN KID'
      ENDIF
      GOTO 1

 9000 RETURN
  END SUBROUTINE WR_MENU
END MODULE wrmenu
