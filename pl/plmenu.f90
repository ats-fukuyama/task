!  $Id$

!     ***** TASK/PL MENU *****

  MODULE plmenu

  CONTAINS

    SUBROUTINE pl_menu

      USE plcomm
      USE plparm,ONLY: pl_parm,pl_view
      USE plgout,ONLY: pl_gout
      USE plload,ONLY: pl_load

      IMPLICIT NONE
      INTEGER(ikind)    :: IERR, MODE
      CHARACTER         :: KID
      CHARACTER(LEN=80) :: LINE

    1 CONTINUE
         IERR=0
         WRITE(6,601)
  601    FORMAT('## PL MENU: P,V/PARM  G/graph  L/LOAD  Q/QUIT')

         CALL TASK_KLIN(LINE,KID,MODE,pl_parm)
      IF(MODE.NE.1) GOTO 1

      IF(KID.EQ.'P') THEN
         CALL pl_parm(0,'PL',IERR)
      ELSEIF(KID.EQ.'V') THEN
         CALL pl_view
      ELSEIF(KID.EQ.'G') THEN
         CALL pl_gout
      ELSEIF(KID.EQ.'L') THEN
         CALL pl_load(ierr)
         IF(ierr.ne.0) GO TO 1
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ELSE
         WRITE(6,*) 'XX PLMENU: UNKNOWN KID'
      ENDIF
      GOTO 1

 9000 RETURN
    END SUBROUTINE pl_menu
  END MODULE plmenu
