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
      CHARACTER KID*1
C
      CALL PLINIT
      CALL PLPARF(IERR)
C
    1 WRITE(6,601)
  601 FORMAT('#PL> SELECT : P,V/PARM  Q/QUIT')
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL GUCPTL(KID)
C
      IF(KID.EQ.'P') THEN
         CALL PLPARM(IERR)
      ELSEIF(KID.EQ.'V') THEN
         CALL PLVIEW
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ENDIF
      GOTO 1
C
 9000 STOP
      END
