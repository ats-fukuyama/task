C     $Id$
C
C               ############# TASK/WR #############
C
C                        RAY TRACING
C                     BASED ON TASK/DP
C
C                           A. Fukuyama
C
C                Department of Nuclear Engineering
C                       Kyoto Univerisity
C                     Kyoto 606-8501, Japan
C
C                      V1.00  : 1988 FEB 10
C                      V2.00  : 1996 FEB 04
C                      V3.00  : 1997 AUG 05
C                      V3.10  : 2000 NOV 25
C     --------------------------------------------------------
C
      INCLUDE 'wrcomm.inc'
C
      CHARACTER KID*1
C
      CALL GSOPEN
      CALL PLINIT
      CALL PLPARF(IERR)
      CALL DPINIT
      CALL DPPARF
      CALL WRINIT
      CALL WRPARF
C
    1 WRITE(6,601)
  601 FORMAT('#WR> SELECT : P,V/PARM  R,B/RAY  G/GRAPH  S/SAVE',
     &       '  1,2,3/DISP  F/ROOT  Q/QUIT')
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL GUCPTL(KID)
C
      IF(KID.EQ.'P') THEN
         CALL WRPARM
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
 9000 CALL GSCLOS
      STOP
      END
