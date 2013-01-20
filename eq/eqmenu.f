C     $Id$
C
C   ***** TASK/EQ MENU *****
C
      SUBROUTINE EQMENU
C
      INCLUDE 'eqcomc.inc'
C
      EXTERNAL EQPARM
      CHARACTER KNAM*80,KPNAME*80
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
         WRITE(6,602) 
  601    FORMAT('## EQ MENU: R/RUN  C/CONT  P,V,I/PARM  G/GRAPH',
     &                    '  M/MULT  S,L,K,F/FILE')
 602     FORMAT('            H/RPPL  Q/QUIT')
C
         CALL TASK_KLIN(LINE,KID,MODE,EQPARM)
      IF(MODE.NE.1) GOTO 1
C
      IF(KID.EQ.'R') THEN
         MSTAT=0
         CALL EQCALC(IERR)
            IF(IERR.NE.0) GOTO 1
         MSTAT=1
         CALL EQCALQ(IERR)
            IF(IERR.NE.0) GOTO 1
C
      ELSEIF(KID.EQ.'C') THEN
         IF(MSTAT.GE.1) THEN 
  101       WRITE(6,*) '#EQ> INPUT PP0,PP1,PP2,PJ0,PJ1,PJ2,RIP,HM:'
            READ(5,*,ERR=101,END=1) PP0,PP1,PP2,PJ0,PJ1,PJ2,RIP,HM

            CALL EQLOOP(IERR)
               IF(IERR.NE.0) GO TO 101
            CALL EQTORZ
            CALL EQCALP
            MSTAT=1
            CALL EQCALQ(IERR)
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
      write(6,'(A)') KNAMEQ
      write(6,'(A)') KNAMEQ2
         CALL KTRIM(KNAMEQ,KL)
   10    WRITE(6,*) '#EQ> INPUT : EQDATA FILE NAME : ',KNAMEQ(1:KL)
         READ(5,'(A80)',ERR=10,END=9000) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMEQ=KNAM
C
         CALL EQREAD(IERR)
         IF(IERR.NE.0) GOTO 10
         CALL EQCALQ(IERR)
         MSTAT=2
C
      ELSEIF(KID.EQ.'K') THEN
         CALL KTRIM(KNAMEQ,KL)
   11    WRITE(6,*) '#EQ> INPUT : EQDSK FILE NAME : ',KNAMEQ(1:KL)
         READ(5,'(A80)',ERR=11,END=9000) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMEQ=KNAM
C
         MODELG=5
         CALL EQREAD(IERR)
         IF(IERR.NE.0) GOTO 11
         CALL EQCALQ(IERR)
         MSTAT=2
C
      ELSEIF(KID.EQ.'F') THEN
         IF(MSTAT.NE.0) CALL EQMETRIC(IERR)
C
      ELSEIF(KID.EQ.'H') THEN
         if(MSTAT .eq. 0) return
         call read_rppl(ierr)
         if(ierr .eq. 0) call eqrppl(IERR)
         MSTAT=3
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
