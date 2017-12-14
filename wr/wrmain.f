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
      USE plcomm
      USE plinit,ONLY: PL_INIT,PL_PARM
      INCLUDE 'wrcomm.inc'

      CALL mtx_initialize

      IF(NRANK.EQ.0) THEN
         WRITE(6,*) '## TASK/WR 2017/12/06'
         CALL GSOPEN
      END IF
      CALL MPSYNC

      CALL PL_INIT
      CALL EQINIT
      CALL DPINIT
      CALL WRINIT
      IF(NRANK.EQ.0) THEN
         OPEN(7,STATUS='SCRATCH')
         CALL PL_PARM(1,'plparm',IERR)
         CALL EQPARM(1,'eqparm',IERR)
         CALL DPPARM(1,'dpparm',IERR)
         CALL WRPARM(1,'wrparm',IERR)
      END IF
      CALL MPSYNC
      CALL WRPRBC

      CALL WRMENU

      IF(NRANK.EQ.0) THEN
         CLOSE(7)
         CALL GSCLOS
      END IF

      CALL mtx_finalize

      STOP
      END
