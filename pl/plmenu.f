C     $Id$
C
C     ***** TASK/PL MENU *****
C
      SUBROUTINE PLMENU
C
      INCLUDE 'plcomm.inc'
C
      CHARACTER KID*1
C
    1 WRITE(6,601)
  601 FORMAT('## PL MENU: P,V/PARM  Q/QUIT')
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL GUCPTL(KID)
C
      IF(KID.EQ.'P') THEN
         CALL PLPARM(IERR)
      ELSEIF(KID.EQ.'V') THEN
         CALL PLVIEW
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ENDIF
      GOTO 1
C
 9000 RETURN
      END
