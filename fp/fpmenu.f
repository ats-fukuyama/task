C     $Id$
C
C     ***** TASK/FP MENU *****
C
      SUBROUTINE FPMENU
C
      INCLUDE 'fpcomm.inc'
C
      CHARACTER KID*1,LINE*80
      EXTERNAL FPPARM
      INTEGER NTG1,NTG2
      DATA NTG1/0/,NTG2/0/
C
    1 CONTINUE
         IERR=0
         WRITE(6,601)
  601    FORMAT('## FP MENU: R:RUN C:CONT P,V:PARAM G,F:GRAPH',
     &                     ' I:RESET L,S:FILE Y:COEF Q:QUIT')
C
         CALL TASK_KLIN(LINE,KID,MODE,FPPARM)
      IF(MODE.NE.1) GOTO 1
C
      IF (KID.EQ.'R') THEN
         TIMEFP=0.D0
         NTG1=0
         NTG2=0
         CALL FPPREP(IERR)
         IF(IERR.NE.0) GOTO 1
         CALL FPLOOP
      ELSEIF (KID.EQ.'C') THEN
         CALL FPLOOP
      ELSEIF (KID.EQ.'P') THEN
         CALL FPPARM(0,'FP',IERR)
      ELSEIF (KID.EQ.'V') THEN
         CALL PLVIEW
         CALL FPVIEW
      ELSEIF (KID.EQ.'G') THEN
         CALL FPGRAF
      ELSEIF (KID.EQ.'F') THEN
         NGRAPH_SAVE=NGRAPH
         NGRAPH=0
         CALL FPGRAF
         NGRAPH=NGRAPH_SAVE
      ELSEIF (KID.EQ.'W') THEN
         IF(NTG1.NE.0) THEN
            CALL FPSGLB
            CALL FPWRT2
            CALL FPSPRF
            CALL FPWRT1
         ELSE
            WRITE(6,*) 'XX no data to write'
         ENDIF
      ELSEIF (KID.EQ.'Y') THEN
         TIMEFP=0.D0
         CALL FPPREP(IERR)
         IF(IERR.NE.0) GOTO 1
         DO NSA=1,NSAMAX
            CALL FPCOEF(NSA)
         END DO
         CALL FPSGLB
         CALL FPWRT2
         CALL FPSPRF
         CALL FPWRT1
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
         CONTINUE
      ELSE
         WRITE(6,*) 'XX FPMENU: UNKNOWN KID'
      END IF
C
      GO TO 1
C
 9000 RETURN
      END
