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
      CHARACTER KID*1
C
C      ON REAL*8 UNDERFLOW CALL ERUFL8
C
      WRITE(6,*) '*** TASK/FP V3.1 [02/09/20] ***'
      CALL GSOPEN
      CALL PLINIT
      CALL FPINIT
      CALL PLPARF(IERR)
      CALL FPPARF
C
    1 WRITE(6,601)
  601 FORMAT('#FP> SELECT : R:RUN  C:CONT  P,V:PARAM  G,F:GRAPH'/
     &       '              I:RESET  W:WRITE  Y:COEF  Q:QUIT')
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL GUCPTL(KID)
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
         CALL FPPARM
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
      END IF
      GO TO 1
C
 9000 CALL GSCLOS
      STOP
      END
