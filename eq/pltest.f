C     $Id$
C
      IMPLICIT COMPLEX*16(C),REAL*8(A,B,D-F,H,O-Z)
      CHARACTER*80 KNAMEQ
C
      CALL GSOPEN
    1 WRITE(6,*) '#EQ> INPUT 3:TASK/EQ 5:EQDSK'
      READ(5,*,ERR=1,END=9000) MODELG
      IF(MODELG.EQ.3) THEN
         KNAMEQ='eqdata'
      ELSEIF(MODELG.EQ.5) THEN
         KNAMEQ='eqdskdata'
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELG FOR EQLOAD'
         GOTO 1
      ENDIF
      CALL EQLOAD(MODELG,KNAMEQ,IERR)
      IF(IERR.EQ.1) GOTO 1
      CALL EQSETP
      NRMAX1=50
      NTHMAX1=16
      NSUMAX1=41
      CALL EQPSIC(NRMAX1,NTHMAX1,NSUMAX1,IERR)
C
    2 WRITE(6,*) '#EQ> INPUT PSIN'
      READ(5,*,ERR=2,END=9000) PSIN
      IF(PSIN.LT.0.D0) GOTO 9000
      CALL GETPP(PSIN,PP)
      CALL GETQP(PSIN,QP)
      CALL GETRMN(PSIN,RRMINL)
      CALL GETRMX(PSIN,RRMAXL)
      WRITE(6,'(A,1P4E12.4)') 'PP,QP,RRMIN,RRMAX =',PP,QP,RRMINL,RRMAXL
      GOTO 2
C
 9000 CALL GSCLOS
      STOP
      END
