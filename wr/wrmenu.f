C     $Id$
C
C     ***** TASK/WR MENU *****
C
      SUBROUTINE WRMENU
C
      INCLUDE 'wrcomm.inc'
C
      CHARACTER KID*1,LINE*80
C
    1 WRITE(6,601)
  601 FORMAT('## WR MENU: P,V/PARM  R,B/RAY  G/GRAPH  S/SAVE',
     &       '  1,2,3/DISP  F/ROOT  Q/QUIT')
      CALL WRKLIN(LINE,KID,MODE)
      IF(MODE.NE.1) GOTO 1
C
      IF(KID.EQ.'P') THEN
         CALL WRPARM(KID)
      ELSEIF(KID.EQ.'V') THEN
         CALL PLVIEW
         CALL DPVIEW
         CALL WRVIEW
      ELSEIF(KID.EQ.'1') THEN
         CALL DPGRP1
      ELSEIF(KID.EQ.'2') THEN
         CALL DPCONT
      ELSEIF(KID.EQ.'3') THEN
         CALL DPCONTX
      ELSEIF(KID.EQ.'F') THEN
         CALL DPROOT
      ELSEIF(KID.EQ.'R') THEN
         CALL WRCALC
      ELSEIF(KID.EQ.'B') THEN
         CALL WRBEAM
      ELSEIF(KID.EQ.'G') THEN
         CALL WRGOUT
      ELSEIF(KID.EQ.'S') THEN
         CALL WRSAVE
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ENDIF
      GOTO 1
C
 9000 RETURN
      END
C
C     ***** INPUT KID or LINE *****
C                   MODE=0: LINE INPUT 
C                        1: KID INPUT
C                        2: PARM INPUT
C                        3: NEW PROMPT
C
      SUBROUTINE WRKLIN(LINE,KID,MODE)
C
      CHARACTER LINE*80,KID*1
C
      READ(5,'(A80)',ERR=2,END=3) LINE
C
      ID=0
      DO I=1,80
         IF(LINE(I:I).EQ.'=') ID=1
      ENDDO
      IF(ID.EQ.1) THEN
         CALL WRPARL(LINE)
         MODE=2
         RETURN
      ENDIF
C
      KID=LINE(1:1)
      CALL GUCPTL(KID)
      IF((KID.GE.'0'.AND.KID.LE.'9').OR.
     &   (KID.GE.'A'.AND.KID.LE.'Z')) THEN
         MODE=1
         RETURN
      ENDIF
C
      KID=' '
      MODE=0
      RETURN
C
    2 WRITE(6,*) 'XX INPUT ERROR !'
      MODE=3
      RETURN
C
    3 KID='Q'
      MODE=1
      RETURN
      END
