C     $Id$
C
C   ************************************************
C   **  Grad-Shafranov Equation (Fixed Boundary)  **
C   ************************************************
C
C   PSIN: 0 on axis, 1 on boundary
C
      INCLUDE 'eqcomc.inc'
C
      CHARACTER KNAMEQ1*72,KNAM*72,KPNAME*32
      CHARACTER KID*1,LINE*80
C
      WRITE(6,*) '## TASK/EQ 2003/09/28'
      MSTAT=0
      CALL GSOPEN
      OPEN(7,STATUS='SCRATCH')
      CALL EQINIT
      KPNAME='eqparm'
      CALL EQPARF(KPNAME)
C
    1 WRITE(6,*) '#EQ> SELECT : R/RUN C/CONT P,V,I/PARM G/GRAPH',
     &           ' M/MULTI S,L/FILE U/UFILE Q/QUIT'
C
      CALL EQKLIN(LINE,KID,MODE,IERR)
      IF(MODE.EQ.1.AND.IERR.EQ.2) GOTO 100
      IF(MODE.EQ.1.AND.IERR.EQ.3) GOTO 9000
      IF(MODE.NE.1) GOTO 1
C
  100 CONTINUE
      IF(KID.EQ.'R') THEN
         CALL EQMESH
         CALL EQPSIN
         CALL EQLOOP(IERR)
            IF(IERR.NE.0) GOTO 1
         CALL EQTORZ
         CALL EQCALP
         CALL EQSETP
         NRMAX1=NRMAX
         NTHMAX1=NTHMAX
         NSUMAX1=NSUMAX
         CALL EQPSIC(NRMAX1,NTHMAX1,NSUMAX1,IERR)
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
            CALL EQSETP
            NRMAX1=NRMAX
            NTHMAX1=NTHMAX
            NSUMAX1=NSUMAX
            CALL EQPSIC(NRMAX1,NTHMAX1,NSUMAX1,IERR)
         ELSE
            WRITE(6,*) 'XX No data for continuing calculation!'
         ENDIF
C
      ELSEIF(KID.EQ.'P') THEN
         CALL EQPARM(KID)
         IF(KID.EQ.'Q') GOTO 9000
C
      ELSEIF(KID.EQ.'V') THEN
         CALL EQVIEW
C
      ELSEIF(KID.EQ.'I') THEN
   20    WRITE(6,'(A,A)') '#EQ> INPUT : EQPARM FILE NAME : ',KPNAME
         READ(5,'(A72)',ERR=20,END=9000) KPNAME
         CALL EQPARF(KPNAME)
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
   10    WRITE(6,*) '#EQ> INPUT : EQDATA FILE NAME : ',KNAMEQ1
         READ(5,'(A72)',ERR=10,END=9000) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMEQ1=KNAM
C
         CALL EQLOAD(3,KNAMEQ1,IERR)
         IF(IERR.EQ.1) GOTO 10
         CALL EQSETP
         NRMAX1=NRMAX
         NTHMAX1=NTHMAX
         NSUMAX1=NSUMAX
         CALL EQPSIC(NRMAX1,NTHMAX1,NSUMAX1,IERR)
         MSTAT=2
C
      ELSEIF(KID.EQ.'U') THEN
         KNAMEQ1=KNAMEQ
   30    WRITE(6,*) '#EQ> INPUT : EQDATA FILE NAME : ',KNAMEQ1
         READ(5,'(A72)',ERR=30,END=9000) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMEQ1=KNAM
C
         CALL EQUFIN(KNAMEQ1,IERR)
         IF(IERR.EQ.1) GOTO 30
         CALL EQSETP
         NRMAX1=NRMAX
         NTHMAX1=NTHMAX
         NSUMAX1=NSUMAX
         CALL EQPSIC(NRMAX1,NTHMAX1,NSUMAX1,IERR)
         MSTAT=2
C
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ELSE
         WRITE(6,*) 'X UNKNOWN COMMAND CHAR'
      ENDIF
      GOTO 1
C
 9000 CALL GSCLOS
      STOP
      END
C
C     ***** INPUT KID or LINE *****
C                   MODE=0: LINE INPUT 
C                        1: KID INPUT
C                        2: PARM INPUT
C                        3: NEW PROMPT
C
      SUBROUTINE EQKLIN(LINE,KID,MODE,IERR)
C
      CHARACTER LINE*80,KID*1
C
      READ(5,'(A80)',ERR=2,END=3) LINE
C
      IERR=0
      ID=0
      DO I=1,80
         IF(LINE(I:I).EQ.'=') ID=1
      ENDDO
      IF(ID.EQ.1) THEN
         CALL EQPARL(LINE)
         MODE=2
         RETURN
      ENDIF
C
      KID=LINE(1:1)
      CALL GUCPTL(KID)
      IF(KID.EQ.'R'.OR.
     &   KID.EQ.'C'.OR.
     &   KID.EQ.'P'.OR.
     &   KID.EQ.'V'.OR.
     &   KID.EQ.'I'.OR.
     &   KID.EQ.'G'.OR.
     &   KID.EQ.'M'.OR.
     &   KID.EQ.'S'.OR.
     &   KID.EQ.'L'.OR.
     &   KID.EQ.'U'.OR.
     &   KID.EQ.'Q') THEN
         MODE=1
         RETURN
      ENDIF
C
      KID=' '
      MODE=0
      RETURN
C
    2 WRITE(6,*) 'XX INPUT ERROR !'
      MODE=3
      RETURN
C
    3 KID='Q'
      MODE=1
      RETURN
C
    4 KID=LINE(1:1)
      CALL GUCPTL(KID)
      IF(KID.EQ.'G') THEN
         MODE=1
         IERR=2
      ELSEIF (KID.EQ.'Q') THEN
         MODE=1
         IERR=3
      ENDIF
      RETURN
      END
