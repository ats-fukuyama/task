C     $Id$
C
C ************************************************************
C
C                    3D FOKKER-PLANCK CODE
C
C                  TASK/FP  V2.0  1992/12/12
C                           V2.1  1993/09/16
C                           V2.2  1997/03/18
C                           V3.0  1997/08/05
C
C                        PROGRAMMED BY
C                         A. FUKUYAMA
C                      OKAYAMA UNIVERSITY
C
C ************************************************************
C
      INCLUDE 'fpcomm.inc'
C
      CHARACTER KID*1,LINE*80
C
C     ------ INITIALIZATION ------
C
      WRITE(6,*) '*** TASK/FP V3.2 [04/03/21] ***'
      OPEN(33,STATUS='SCRATCH',FORM='FORMATTED')
C
      CALL GSOPEN
      CALL GUTIME(GTCPU1)
      CALL PLINIT
      CALL FPINIT
      CALL PLPARF(IERR)
      CALL FPPARF
C
      IERR=0
      NTSMAX=NTMAX
C
    1 WRITE(6,601)
  601 FORMAT('#FP> SELECT : R:RUN C:CONT P,V:PARAM G,F:GRAPH',
     &                    ' I:RESET W:WRITE Y:COEF Q:QUIT')
C
C
      CALL FPKLIN(LINE,KID,MODE,NTMAX,NTSMAX,IERR)
      IF(MODE.EQ.1.AND.IERR.EQ.2) GOTO 100
      IF(MODE.EQ.1.AND.IERR.EQ.3) GOTO 9000
      IF(MODE.NE.1) GOTO 1
C
  100 CONTINUE
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
 9000 CALL GSCLOS
      CLOSE(33)
      CALL GUTIME(GTCPU2)
      WRITE(6,666) GTCPU2-GTCPU1
  666 FORMAT(' ','#      CPU TIME :   ',F8.3,' SEC ')
      STOP
      END
C
C     ***** INPUT KID or LINE *****
C                   MODE=0: LINE INPUT 
C                        1: KID INPUT
C                        2: PARM INPUT
C                        3: NEW PROMPT
C
      SUBROUTINE FPKLIN(LINE,KID,MODE,NTMAX,NTSMAX,IERR)
C
      CHARACTER LINE*80,KID*1
C
      READ(5,'(A80)',ERR=2,END=3) LINE
C
      IF(IERR.NE.0) GOTO 4
      ID=0
      DO I=1,80
         IF(LINE(I:I).EQ.'=') ID=1
      ENDDO
      IF(ID.EQ.1) THEN
         CALL FPPARL(LINE)
         IF(NTSMAX.NE.NTMAX) NTSMAX=NTMAX
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
C
    4 KID=LINE(1:1)
      CALL GUCPTL(KID)
      IF(KID.EQ.'G') THEN
         MODE=1
         IERR=2
      ELSEIF (KID.EQ.'Q') THEN
         MODE=1
         IERR=3
      ENDIF
      RETURN
      END
