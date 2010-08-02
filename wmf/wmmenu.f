C     $Id$
C
C     ***** TASK/WM MENU *****
C
      SUBROUTINE WMMENU
C
      INCLUDE 'wmcomm.inc'
C
      EXTERNAL WMPARM
      CHARACTER KID*1,LINE*80
      SAVE INIT,NFILEINI
      DATA INIT/0/
C
      IF(INIT.EQ.0) THEN
         NFILEINI=0
         INIT=1
      ENDIF
C
    1 CONTINUE
         IF(MYRANK.EQ.0) THEN
            WRITE(6,601)
  601       FORMAT('## WM MENU: P,V/PARM R/RUN D0-3/AMP F/ROOT ',
     &      'G/GRAPH T/TAE O/OUT S,W/SAVE Q/QUIT')
            CALL TASK_KLIN(LINE,KID,MODE,WMPARM)
         ENDIF
         CALL MPBCIA(MODE)
         IF(MODE.EQ.2) CALL WMPRBC
      IF(MODE.NE.1) GOTO 1
C
    2 CONTINUE
         CALL MPBCKA(KID)
C
         IF (KID.EQ.'P') THEN
            IF(MYRANK.EQ.0) CALL WMPARM(0,'WM',IERR)
            CALL MPSYNC
            CALL WMPRBC
            KID=' '
         ELSE IF(KID.EQ.'V') THEN
            IF(MYRANK.EQ.0) CALL WMVIEW
            CALL MPSYNC
            KID=' '
C
C        *** WAVE CALCULATION ***
C
         ELSEIF (KID.EQ.'R') THEN
            CALL WMLOOP(IERR)
            CALL MPSYNC
            IF(IERR.NE.0) GOTO 1
            KID=' '
C
C        *** AMPLITUDE SURVEY ***
C
      ELSEIF(KID.EQ.'D') THEN
         READ(LINE(2:),*,ERR=1,END=1) NID
         IF(NID.EQ.0) THEN
            CALL WMAM0D(KID,LINE)
         ELSEIF(NID.EQ.1) THEN
            CALL WMAM1D(KID,LINE)
         ELSEIF(NID.EQ.2) THEN
            CALL WMAM2D(KID,LINE)
         ELSEIF(NID.EQ.3) THEN
            CALL WMSCAN(KID,LINE)
         ELSE
            WRITE(6,*) 'XX WMMENU: unknown NID'
         ENDIF
C
C        *** FIND ROOT ***
C
         ELSE IF (KID.EQ.'F') THEN
            CALL WMEIGN(KID,LINE)
C
C        *** GRAPHICS ***
C
         ELSE IF (KID.EQ.'G') THEN
            IF(MYRANK.EQ.0) CALL WMGOUT
            CALL MPSYNC
            KID=' '
C
C        *** FILE OUTPUT ***
C
         ELSE IF (KID.EQ.'S') THEN
            IF(MYRANK.EQ.0) THEN
               CALL WMSAVE
            ENDIF
            CALL MPSYNC
            KID=' '
         ELSE IF (KID.EQ.'W') THEN
            IF(MYRANK.EQ.0) THEN
               IF(NFILEINI.EQ.0) THEN
                  REWIND(26)
                  NFILEINI=1
               ENDIF
               CALL WMWOUT
            ENDIF
            CALL MPSYNC
            KID=' '
C
C        *** TAE FREQUENCY ***
C
         ELSE IF (KID.EQ.'T') THEN
            MODEEG=0
            CALL wmfem_setg(IERR)
               IF(IERR.NE.0) GOTO 1
            CALL wmfem_setj(IERR)
               IF(IERR.NE.0) GOTO 1
            IF(MYRANK.EQ.0) CALL WMTAE
            CALL MPSYNC
            KID=' '
C
C        *** Pabs(r,s) output for TOPICS ***
C
         ELSE IF (KID.EQ.'O') THEN
            CALL WMLOOP(IERR,0)
            CALL MPSYNC
            IF(IERR.NE.0) GOTO 1
            KID=' '
c$$$            CALL WM_TOPICS(IERR)
C
         ELSE IF(KID.EQ.'H') THEN
            WRITE(6,*) '# KID: P: PARAMETER INPUT (VARNAME = VALUE)'
            WRITE(6,*) '       V:  VIEW PARAMETERS'
            WRITE(6,*) '       R:  WAVE EXCITED BY EXTERNAL ANTENNA'
            WRITE(6,*) '       D0: AMPLITUDE OF INTERNALLY EXCITED WAVE'
            WRITE(6,*) '       D1: FREQUENCY SCAN OF AMPLITUDE'
            WRITE(6,*) '       D2: COMPLEX FREQUENCY SCAN OF AMPLITUDE'
            WRITE(6,*) '       D3: PARAMETER SCAN OF EIGEN VALUE'
            WRITE(6,*) '       F:  FIND EIGEN VALUE'
            WRITE(6,*) '       G: GRAPHICS'
            WRITE(6,*) '       W: FILE 26 OUTPUT: FIELD DATA'
            WRITE(6,*) '       H: HELP'
            WRITE(6,*) '       Q: QUIT'
            KID=' '
C
         ELSE IF(KID.EQ.'X') THEN
            KID=' '
C
         ELSE IF(KID.EQ.'Q') THEN
            GOTO 9000
         ELSE
            IF(MYRANK.EQ.0) WRITE(6,*) 'XX WMMENU: UNKNOWN KID'
            KID=' '
         END IF
C
      IF(KID.NE.' '.AND.KID.NE.'Q') GOTO 2
      GO TO 1
C
 9000 RETURN
      END
