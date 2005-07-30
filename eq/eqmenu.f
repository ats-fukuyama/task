C     $Id$
C
C   ***** TASK/EQ MENU *****
C
      SUBROUTINE EQMENU
C
      INCLUDE 'eqcomc.inc'
C
      EXTERNAL EQPARM
      CHARACTER KNAMEQ1*80,KNAM*80,KPNAME*80
      CHARACTER KID*1,LINE*80
      SAVE INIT,MSTAT,KPNAME
      DATA INIT/0/
C
      IF(INIT.EQ.0) THEN
         MSTAT=0
         KPNAME='eqparm'
         INIT=1
      ENDIF
C
    1 CONTINUE
         IERR=0
         WRITE(6,601) 
  601    FORMAT('## EQ MENU: R/RUN  C/CONT  P,V,I/PARM  G/GRAPH',
     &                  '  M/MULT  S,L,K/FILE  Q/QUIT')
C
         CALL TASK_KLIN(LINE,KID,MODE,EQPARM)
      IF(MODE.NE.1) GOTO 1
C
      IF(KID.EQ.'R') THEN
         MSTAT=1
         CALL EQCALC(IERR)
            IF(IERR.NE.0) GOTO 1
         NRMAX1=NRMAX
         NTHMAX1=NTHMAX
         NSUMAX1=NSUMAX
         CALL EQCALQ(NRMAX1,NTHMAX1,NSUMAX1,IERR)
            IF(IERR.NE.0) GOTO 1
         MSTAT=1
C
      ELSEIF(KID.EQ.'C') THEN
         IF(MSTAT.EQ.1) THEN 
  101       WRITE(6,*) '#EQ> INPUT PP0,PP1,PP2,PJ0,PJ1,PJ2,RIP,HM:'
            READ(5,*,ERR=101,END=1) PP0,PP1,PP2,PJ0,PJ1,PJ2,RIP,HM
            CALL EQLOOP(IERR)
               IF(IERR.NE.0) GOTO 1
            CALL EQTORZ
            CALL EQCALP
            NRMAX1=NRMAX
            NTHMAX1=NTHMAX
            NSUMAX1=NSUMAX
            CALL EQCALQ(NRMAX1,NTHMAX1,NSUMAX1,IERR)
            IF(IERR.NE.0) GOTO 1
         ELSE
            WRITE(6,*) 'XX No data for continuing calculation!'
         ENDIF
C
      ELSEIF(KID.EQ.'P') THEN
         CALL EQPARM(0,'EQ',IERR)
C
      ELSEIF(KID.EQ.'V') THEN
         CALL EQVIEW
C
      ELSEIF(KID.EQ.'I') THEN
   20    WRITE(6,'(A,A)') '#EQ> INPUT : EQPARM FILE NAME : ',KPNAME
         READ(5,'(A80)',ERR=20,END=9000) KPNAME
         CALL EQPARM(1,KPNAME,IERR)
C
      ELSEIF(KID.EQ.'G') THEN
         CALL EQGOUT(MSTAT)
C
      ELSEIF(KID.EQ.'M') THEN
  102    WRITE(6,*) '#EQ> INPUT WHAT IS YOUR NUMBER OF TIMES? (1-5):'
         READ(5,*,ERR=102,END=1) NTIMES
         IF (NTIMES.LT.5) THEN
  103    WRITE(6,*) '#EQ> INPUT NEXT PV0:'
         READ(5,*,ERR=103,END=1) PV0
         ENDIF
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            PSIO(NTG,NSG,NTIMES) = PSI(NTG,NSG)
         ENDDO
         ENDDO
C
      ELSEIF(KID.EQ.'S') THEN
         CALL EQSAVE
C
      ELSEIF(KID.EQ.'L') THEN
         KNAMEQ1=KNAMEQ
         CALL KTRIM(KNAMEQ1,KL)
   10    WRITE(6,*) '#EQ> INPUT : EQDATA FILE NAME : ',KNAMEQ1(1:KL)
         READ(5,'(A80)',ERR=10,END=9000) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMEQ1=KNAM
C
         CALL EQLOAD(3,KNAMEQ1,IERR)
         IF(IERR.EQ.1) GOTO 10
         NRMAX1=NRMAX
         NTHMAX1=NTHMAX
         NSUMAX1=NSUMAX
         CALL EQCALQ(NRMAX1,NTHMAX1,NSUMAX1,IERR)
         MSTAT=2
C
      ELSEIF(KID.EQ.'K') THEN
         KNAMEQ1=KNAMEQ
         CALL KTRIM(KNAMEQ1,KL)
   11    WRITE(6,*) '#EQ> INPUT : EQDATA FILE NAME : ',KNAMEQ1(1:KL)
         READ(5,'(A80)',ERR=11,END=9000) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMEQ1=KNAM
C
         CALL EQLOAD(5,KNAMEQ1,IERR)
         IF(IERR.NE.0) GOTO 11
         NRMAX1=NRMAX
         NTHMAX1=NTHMAX
         NSUMAX1=NSUMAX
         CALL EQCALQ(NRMAX1,NTHMAX1,NSUMAX1,IERR)
         MSTAT=2
C
      ELSE IF(KID.EQ.'X'.OR.KID.EQ.'#') THEN
         CONTINUE
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ELSE
         WRITE(6,*) 'XX EQMENU: UNKNOWN KID'
      ENDIF
      GOTO 1
C
 9000 RETURN
      END
