C     $Id$
C
C     ***** TASK/WR MENU *****
C
      SUBROUTINE WRMENU
C
      INCLUDE 'wrcomm.inc'
C
      EXTERNAL WRPARM
      CHARACTER KID*1,LINE*80
C
    1 CONTINUE
         IERR=0
         WRITE(6,601)
  601    FORMAT('## WR MENU: P,V/PARM  R,B/RAY  G/GRAPH  S/SAVE',
     &          '  Dn/DISP  F/ROOT  Q/QUIT')
         CALL TASK_KLIN(LINE,KID,MODE,WRPARM)
      IF(MODE.NE.1) GOTO 1
C
      IF(KID.EQ.'P') THEN
         CALL WRPARM(0,'WR',IERR)
      ELSEIF(KID.EQ.'V') THEN
         CALL PLVIEW
         CALL DPVIEW
         CALL WRVIEW
      ELSEIF(KID.EQ.'D') THEN
         READ(LINE(2:),*,ERR=1,END=1) NID
         IF(NID.EQ.1) THEN
            CALL DPGRP1
         ELSEIF(NID.EQ.2) THEN
            CALL DPCONT
         ELSEIF(NID.EQ.3) THEN
            CALL DPCONTX
         ELSE
            WRITE(6,*) 'XX WRMENU: unknown NID'
         ENDIF
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
      ELSE
         WRITE(6,*) 'XX WRMENU: UNKNOWN KID'
      ENDIF
      GOTO 1
C
 9000 RETURN
      END
