C     $Id$
C
C     ***** TASK/FP MENU *****
C
      SUBROUTINE FPMENU
C
      INCLUDE 'fpcomm.inc'
C
      CHARACTER KID*1,LINE*80
C
    1 IERR=0
      WRITE(6,601)
  601 FORMAT('## FP MENU: R:RUN C:CONT P,V:PARAM G,F:GRAPH',
     &                  ' I:RESET W:WRITE Y:COEF Q:QUIT')
C
      CALL FPKLIN(LINE,KID,MODE)
      IF(MODE.NE.1) GOTO 1
C
      IF (KID.EQ.'R') THEN
         TIMEFP=0.D0
         NTG1=0
         NTG2=0
         CALL FPPREP(IERR)
         IF(IERR.NE.0) GOTO 1
         CALL FPFINI
         CALL FPSAVI
         CALL FPLOOP
      ELSEIF (KID.EQ.'C') THEN
         IF(KID.EQ.'F') THEN
            NTG1=0
            NTG2=0
         ENDIF
         CALL FPLOOP
      ELSEIF (KID.EQ.'P') THEN
         CALL FPPARM(KID)
      ELSEIF (KID.EQ.'V') THEN
         CALL PLVIEW
         CALL FPVIEW
      ELSEIF (KID.EQ.'G') THEN
         CALL FPGRAF
      ELSEIF (KID.EQ.'F') THEN
         CALL FPFOUT
      ELSEIF (KID.EQ.'W') THEN
         CALL FPSGLB
         CALL FPWRIT
      ELSEIF (KID.EQ.'Y') THEN
         TIMEFP=0.D0
         CALL FPPREP(IERR)
         IF(IERR.NE.0) GOTO 1
         CALL FPFINI
         CALL FPCOEF
         CALL FPSGLB
         CALL FPWRIT
      ELSEIF (KID.EQ.'I') THEN
         NTG1=0
         NTG2=0
      ELSEIF (KID.EQ.'S') THEN
         CALL FPSAVE
      ELSEIF (KID.EQ.'L') THEN
         CALL FPLOAD
      ELSEIF (KID.EQ.'Q') THEN
         GO TO 9000
C
      ELSE IF(KID.EQ.'X'.OR.KID.EQ.'#') THEN
         KID=' '
      ELSE
         WRITE(6,*) 'XX FPMAIN: UNKNOWN KID'
         KID=' '
      END IF
C
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
      SUBROUTINE FPKLIN(LINE,KID,MODE)
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
         CALL FPPARL(LINE)
         MODE=2
         RETURN
      ENDIF
C
      KID=LINE(1:1)
      CALL GUCPTL(KID)
      IF(KID.EQ.'P'.OR.
     &   KID.EQ.'V'.OR.
     &   KID.EQ.'R'.OR.
     &   KID.EQ.'C'.OR.
     &   KID.EQ.'G'.OR.
     &   KID.EQ.'F'.OR.
     &   KID.EQ.'W'.OR.
     &   KID.EQ.'Y'.OR.
     &   KID.EQ.'I'.OR.
     &   KID.EQ.'I'.OR.
     &   KID.EQ.'S'.OR.
     &   KID.EQ.'L'.OR.
     &   KID.EQ.'Q') THEN
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
