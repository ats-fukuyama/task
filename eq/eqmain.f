C     $Id$
C   ************************************************
C   **  Grad-Shafranov Equation (Fixed Boundary)  **
C   ************************************************
C
      INCLUDE 'eqcomm.h'
C
      CHARACTER KNAM*32
      CHARACTER KID*1
C
      WRITE(6,*) '## TASK/EQ V2.0 1999/06/10'
      MODE=0
      CALL GSOPEN
      CALL EQINIT
      CALL EQPARF
C
    1 WRITE(6,*) '#EQ> SELECT : R/RUN C/CONT P,V/PARM G/GRAPH S,L/FILE',
     &           ' Q/QUIT'
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL GUCPTL(KID)
C
      IF(KID.EQ.'R') THEN
         CALL EQMESH
         CALL EQPSIN
         CALL EQLOOP
         CALL EQTORZ
         CALL EQSETP
         CALL EQPSIC(NRMAX,NTHMAX,NSUMAX)
         MODE=1
C
      ELSEIF(KID.EQ.'C') THEN
         IF(MODE.EQ.1) THEN 
  101       WRITE(6,*) '# INPUT PP0,PP1,PP2,PJ0,PJ1,PJ2 :'
            READ(5,*,ERR=101,END=1) PP0,PP1,PP2,PJ0,PJ1,PJ2
            CALL EQLOOP
            CALL EQTORZ
            CALL EQSETP
            CALL EQPSIC(NRMAX,NTHMAX,NSUMAX)
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
         IF(MODE.EQ.0) THEN
            WRITE(6,*) 'XX DATA UNDEFINED!'
         ELSE
            IF(MODE.EQ.1) THEN
               CALL EQGRAP
            ENDIF
            CALL SHOWPS
            CALL EQDRAW(1)
            CALL EQDRAW(2)
            CALL EQDRAW(3)
            CALL EQDRAW(4)
         ENDIF
C
      ELSEIF(KID.EQ.'S') THEN
         CALL EQSAVE
      ELSEIF(KID.EQ.'L') THEN
   10    WRITE(6,*) '# INPUT : EQDATA FILE NAME : ',KNAMEQ
         READ(5,'(A32)',ERR=10,END=9000) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMEQ=KNAM
C
         CALL EQLOAD(1,KNAMEQ,IERR)
         IF(IERR.EQ.1) GOTO 10
         CALL EQSETP
         CALL EQPSIC(NRMAX,NTHMAX,NSUMAX)
         MODE=2
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
