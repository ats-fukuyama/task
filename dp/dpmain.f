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
      CHARACTER KPNAME*80
C
      WRITE(6,*) '## TASK/DP 2004/11/08'
      CALL GSOPEN
      CALL PLINIT
      CALL DPINIT
      KPNAME='plparm'
      CALL PLPARF(KPNAME)
      KPNAME='dpparm'
      CALL DPPARF(KPNAME)
C
      CALL DPMENU
C
      CALL GSCLOS
      STOP
      END
