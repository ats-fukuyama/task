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
      CHARACTER KPNAME*80
C
      WRITE(6,601) 
  601 FORMAT('## TASK/WR 2004/11/08')
      CALL GSOPEN
      CALL PLINIT
      CALL DPINIT
      CALL WRINIT
      KPNAME='plparm'
      CALL PLPARF(KPNAME)
      KPNAME='dpparm'
      CALL DPPARF(KPNAME)
      KPNAME='wrparm'
      CALL WRPARF(KPNAME)
      CALL PLDATA_SETN(50,1)
C
      CALL WRMENU
C
      CALL GSCLOS
      STOP
      END
