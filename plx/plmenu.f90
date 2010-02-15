!  $Id$

!     ***** TASK/PL MENU *****

  MODULE plmenu

  CONTAINS

    SUBROUTINE pl_menu

      USE plcomm
      USE plinit

      IMPLICIT NONE
      INTEGER(ikind)    :: IERR, MODE
      CHARACTER         :: KID
      CHARACTER(LEN=80) :: LINE

    1 CONTINUE
         IERR=0
         WRITE(6,601)
  601    FORMAT('## PL MENU: P,V/PARM  Q/QUIT')

         CALL TASK_KLIN(LINE,KID,MODE,pl_parm)
      IF(MODE.NE.1) GOTO 1

      IF(KID.EQ.'P') THEN
         CALL pl_parm(0,'PL',IERR)
      ELSEIF(KID.EQ.'V') THEN
         CALL pl_view
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ELSE
         WRITE(6,*) 'XX PLMENU: UNKNOWN KID'
      ENDIF
      GOTO 1

 9000 RETURN
    END SUBROUTINE pl_menu
  END MODULE plmenu
