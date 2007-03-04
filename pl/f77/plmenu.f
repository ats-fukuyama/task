C     $Id$
C
C     ***** TASK/PL MENU *****
C
      SUBROUTINE PLMENU
C
      INCLUDE 'plcomm.inc'
C
      EXTERNAL PLPARM
      CHARACTER KID*1,LINE*80
C
    1 CONTINUE
         IERR=0
         WRITE(6,601)
  601    FORMAT('## PL MENU: P,V/PARM  Q/QUIT')
C
         CALL TASK_KLIN(LINE,KID,MODE,PLPARM)
      IF(MODE.NE.1) GOTO 1
C
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
C
 9000 RETURN
      END
