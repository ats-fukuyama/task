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
C
      WRITE(6,*) '## TASK/WR 2004/11/08'
C
      CALL GSOPEN
      OPEN(7,STATUS='SCRATCH')
      CALL PLINIT
      CALL EQINIT
      CALL DPINIT
      CALL WRINIT
      CALL PLPARM(1,'plparm',IERR)
      CALL EQPARM(1,'eqparm',IERR)
      CALL DPPARM(1,'dpparm',IERR)
      CALL WRPARM(1,'wrparm',IERR)
C
      CALL WRMENU
C
      CLOSE(7)
      CALL GSCLOS
      STOP
      END
