C     $Id$
C   ************************************************
C   **  Grad-Shafranov Equation (Fixed Boundary)  **
C   ************************************************
C
C      INCLUDE 'eqcomm.h'
      INCLUDE 'eqcomc.h'
C
      CHARACTER KNAMEQ1*72,KNAM*72,KPNAM*72
      CHARACTER KID*1
C
      WRITE(6,*) '## TASK/EQ 2001/08/25'
      MODE=0
      CALL GSOPEN
      CALL EQINIT
      CALL EQPARF
C
    1 WRITE(6,*) '#EQ> SELECT : R/RUN C/CONT P,V,I/PARM G/GRAPH',
     &           ' M/MULTI S,L/FILE U/UFILE Q/QUIT'
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL GUCPTL(KID)
C
      IF(KID.EQ.'R') THEN
         CALL EQMESH
         CALL EQPSIN
         CALL EQLOOP(IERR)
            IF(IERR.NE.0) GOTO 1
         CALL EQTORZ
         CALL EQSETP
         NRMAX1=NRMAX
         NTHMAX1=NTHMAX
         NSUMAX1=NSUMAX
         CALL EQPSIC(NRMAX1,NTHMAX1,NSUMAX1,IERR)
         MODE=1
C
      ELSEIF(KID.EQ.'C') THEN
         IF(MODE.EQ.1) THEN 
  101       WRITE(6,*) '#EQ> INPUT PP0,PP1,PP2,PJ0,PJ1,PJ2,RIP,HM:'
            READ(5,*,ERR=101,END=1) PP0,PP1,PP2,PJ0,PJ1,PJ2,RIP,HM
            CALL EQLOOP(IERR)
               IF(IERR.NE.0) GOTO 1
            CALL EQTORZ
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
         CALL EQPARM
      ELSEIF(KID.EQ.'V') THEN
         CALL EQVIEW
C
      ELSEIF(KID.EQ.'G') THEN
         CALL EQGOUT(MODE)
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
      ELSEIF(KID.EQ.'S') THEN
         CALL EQSAVE
      ELSEIF(KID.EQ.'L') THEN
         KNAMEQ1=KNAMEQ
   10    WRITE(6,*) '#EQ> INPUT : EQDATA FILE NAME : ',KNAMEQ1
         READ(5,'(A72)',ERR=10,END=9000) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMEQ1=KNAM
C
         CALL EQLOAD(1,KNAMEQ1,IERR)
         IF(IERR.EQ.1) GOTO 10
         CALL EQSETP
         NRMAX1=NRMAX
         NTHMAX1=NTHMAX
         NSUMAX1=NSUMAX
         CALL EQPSIC(NRMAX1,NTHMAX1,NSUMAX1,IERR)
         MODE=2
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ELSEIF(KID.EQ.'I') THEN
   20    WRITE(6,*) '#EQ> INPUT : EQPARM FILE NAME : ','eqparm'
         READ(5,'(A72)',ERR=20,END=9000) KPNAM
         CALL EQPARG(KPNAM)
c$$$      ELSEIF(KID.EQ.'U') THEN
c$$$         CALL EQUFRD
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
         MODE=2
      ELSE
         WRITE(6,*) 'X UNKNOWN COMMAND CHAR'
      ENDIF
      IF (KMODE.EQ.1) CALL EQPARF
      GOTO 1
C
 9000 CALL GSCLOS
      STOP
      END
