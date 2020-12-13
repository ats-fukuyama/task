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
    USE wrparm,ONLY: wr_parm,wr_view
    USE wrcalc,ONLY: wr_calc
    USE wrbeam,ONLY: wr_beam
    USE wrgout,ONLY: wr_gout
    USE wrfile,ONLY: wrsave,wrload
    IMPLICIT NONE
    CHARACTER(LEN=1):: KID
    CHARACTER(LEN=80):: LINE
    INTEGER:: NSTAT,IERR,MODE,NID
    EXTERNAL TASK_KLIN

      NSTAT=0

    1 CONTINUE
         IERR=0
         WRITE(6,601)
  601    FORMAT('## WR MENU: P,V/PARM  R,B/RAY  G/GRAPH  S,L/FILE', &
                '  Dn/DISP  F/ROOT  Q/QUIT')
         CALL TASK_KLIN(LINE,KID,MODE,WR_PARM)
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
         CALL WR_CALC(IERR)
         NSTAT=1
      ELSEIF(KID.EQ.'B') THEN
         CALL wr_allocate
         CALL WR_BEAM
         NSTAT=2
      ELSEIF(KID.EQ.'G') THEN
         CALL WR_GOUT(NSTAT)
      ELSEIF(KID.EQ.'S') THEN
         CALL WRSAVE
      ELSEIF(KID.EQ.'L') THEN
         CALL WRLOAD(NSTAT)
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ELSE
         WRITE(6,*) 'XX WRMENU: UNKNOWN KID'
      ENDIF
      GOTO 1

 9000 RETURN
  END SUBROUTINE WR_MENU
END MODULE wrmenu
