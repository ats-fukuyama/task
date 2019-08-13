!     ***** TASK/PL MENU *****

      SUBROUTINE PLMENU

!      INCLUDE 'plcomm.inc'
      IMPLICIT NONE
      EXTERNAL PLPARM
      INTEGER(4)        :: IERR, MODE
      CHARACTER         :: KID
      CHARACTER(LEN=80) :: LINE

    1 CONTINUE
         IERR=0
         WRITE(6,601)
  601    FORMAT('## PL MENU: P,V/PARM  Q/QUIT')

         CALL TASK_KLIN(LINE,KID,MODE,PLPARM)
      IF(MODE.NE.1) GOTO 1

      IF(KID.EQ.'P') THEN
         CALL PLPARM(0,'PL',IERR)
      ELSEIF(KID.EQ.'V') THEN
         CALL PLVIEW
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ELSE
         WRITE(6,*) 'XX EQMENU: UNKNOWN KID'
      ENDIF
      GOTO 1

 9000 RETURN
      END
