C     $Id$
C
C               ############# TASK/PL #############
C
C             BASIC PARAMETERS OF DEVICE AND PLASMAS
C
C                           A. Fukuyama
C
C                Department of Nuclear Engineering
C                       Kyoto Univerisity
C                     Kyoto 606-8501, Japan
C
C                     V1.00  : 1997 AUG 05
C                     V1.10  : 2000 NOV 25
C
C-----------------------------------------------------------------------
C
      INCLUDE 'plcomm.inc'
C
      WRITE(6,*) '## TASK/PL 2004/11/08'
      OPEN(7,STATUS='SCRATCH')
      CALL PLINIT
      CALL PLPARM(1,'plparm',IER)
C
      CALL PLMENU
C
      CLOSE(7)
      STOP
      END
