C     $Id$
C
C               ############# TASK/DP #############
C
C                       DISPERSION RELATION
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
C
C-----------------------------------------------------------------------
C
      INCLUDE 'dpcomm.inc'
C
      WRITE(6,*) '## TASK/DP 2004/11/08'
      CALL GSOPEN
      OPEN(7,STATUS='SCRATCH')
      CALL PLINIT
      CALL DPINIT
      CALL PLPARM(1,'plparm',IERR)
      CALL DPPARM(1,'dpparm',IERR)
C
      CALL DPMENU
C
      CLOSE(7)
      CALL GSCLOS
      STOP
      END
