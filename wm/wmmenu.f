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
  601       FORMAT('## WM MENU: P,V/PARM  R/RUN  A,F,C/AMP  E,S/SCAN  ',
     &      'G/GRAPH T/TAE W/WRITE Q/QUIT')
            CALL TASK_KLIN(LINE,KID,MODE,WMPARM)
         ENDIF
         CALL MPBCIA(MODE)
         IF(MODE.EQ.2) CALL WMPRBC
      IF(MODE.NE.1) GOTO 1
C
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
            CALL WMEXEC(IERR)
            CALL MPSYNC
            IF(IERR.NE.0) GOTO 1
            KID=' '
C
C        *** AMPLITUDE SURVEY ***
C
         ELSEIF (KID.EQ.'A') THEN
            CALL WMAM0D(KID)
         ELSE IF (KID.EQ.'F') THEN
            CALL WMAM1D(KID)
         ELSE IF (KID.EQ.'C') THEN
            CALL WMAM2D(KID)
C
C        *** EIGENMODE ***
C
         ELSE IF (KID.EQ.'E') THEN
            CALL WMEIGN(KID)
         ELSE IF (KID.EQ.'S') THEN
            CALL WMSCAN(KID)
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
            CALL WMSETG(IERR)
               IF(IERR.NE.0) GOTO 1
            CALL WMSETJ(IERR)
               IF(IERR.NE.0) GOTO 1
            IF(MYRANK.EQ.0) CALL WMTAE
            CALL MPSYNC
            KID=' '
C
         ELSE IF(KID.EQ.'H') THEN
            WRITE(6,*) '# KID:  P: PARAMETER INPUT (VARNAME = VALUE)'
            WRITE(6,*) '        V: VIEW PARAMETERS'
            WRITE(6,*) '        R: WAVE EXCITED BY EXTERNAL ANTENNA'
            WRITE(6,*) '        A: AMPLITUDE OF INTERNALLY EXCITED WAVE'
            WRITE(6,*) '        F: REAL FREQUENCY SCAN OF AMPLITUDE'
            WRITE(6,*) '        C: COMPLEX FREQUENCY SCAN OF AMPLITUDE'
            WRITE(6,*) '        E: EIGEN VALUE SEARCH'
            WRITE(6,*) '        S: PARAMETER SCAN OF EIGEN VALUE'
            WRITE(6,*) '        G: GRAPHICS'
            WRITE(6,*) '        W: FILE 26 OUTPUT: FIELD DATA'
            WRITE(6,*) '        H: HELP'
            WRITE(6,*) '        Q: QUIT'
            KID=' '
C
         ELSEIF (KID.EQ.'D') THEN
            CALL WMDEBUG(IERR)
            CALL MPSYNC
            IF(IERR.NE.0) GOTO 1
            KID=' '
C
         ELSE IF(KID.EQ.'Q') THEN
            GOTO 9000
C
         ELSE IF(KID.EQ.'X'.OR.KID.EQ.'#') THEN
            KID=' '
         ELSE
            IF(MYRANK.EQ.0) WRITE(6,*) 'XX WMMAIN: UNKNOWN KID'
            KID=' '
         END IF
C
      IF(KID.NE.' '.AND.KID.NE.'Q') GOTO 2
      GO TO 1
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
      SUBROUTINE WMKLIN(LINE,KID,MODE)
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
         CALL WMPARL(LINE)
         MODE=2
         RETURN
      ENDIF
C
      KID=LINE(1:1)
      CALL GUCPTL(KID)
      IF(KID.GE.'A'.AND.
     &   KID.LE.'Z') THEN
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
