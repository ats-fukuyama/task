C     $Id$
C
      IMPLICIT COMPLEX*16(C),REAL*8(A,B,D-F,H,O-Z)
      CHARACTER*32 KNAMEQ
C
      NRMAX=50
      NTHMAX=16
      NSUMAX=41
C
      CALL GSOPEN
      KNAMEQ='eqdata'
C
      CALL AKEQIN(KNAMEQ,NRMAX,NTHMAX,NSUMAX,
     &            ALPMAX,RAXIS,ZAXIS,IERR)
      IF(IERR.NE.0) GOTO 9000
      WRITE(6,'(A,1PE12.4)') 'ALPMAX=',ALPMAX
      WRITE(6,'(A,1PE12.4)') 'RAXIS =',RAXIS
      WRITE(6,'(A,1PE12.4)') 'ZAXIS =',ZAXIS
C
      CALL EQGETB(BB,RR,RIP,RA,RKAP,RDEL,RB)
      WRITE(6,'(A,1PE12.4)') 'BB    =',BB
      WRITE(6,'(A,1PE12.4)') 'RR    =',RR
      WRITE(6,'(A,1PE12.4)') 'RIP   =',RIP
      WRITE(6,'(A,1PE12.4)') 'RA    =',RA
      WRITE(6,'(A,1PE12.4)') 'RKAP  =',RKAP
      WRITE(6,'(A,1PE12.4)') 'RDEL  =',RDEL
C
    1 WRITE(6,*) '# INPUT R,Z'
      READ(5,*,ERR=1,END=9000) R,Z
C
      CALL RZTOAB(R,Z,ALPHA,BETA)
      WRITE(6,'(A,1P2E12.4)') 'ALPHA,BETA=',ALPHA,BETA
      CALL ABTORZ(ALPHA,BETA,R1,Z1,IERR)
      WRITE(6,'(A,1P2E12.4,I5)') 'R1,Z1,IERR=',R1,Z1,IERR
      ER=R1-R
      EZ=Z1-Z
      ET=SQRT(ER**2+EZ**2)/RAXIS
      WRITE(6,'(A,1P3E12.4,I5)') 'ERRORS=',ER,EZ,ET
C
C      CALL SUBQPA(ALPHA,Q,DQDA)
C      WRITE(6,'(A,1P2E12.4)') 'Q,DQ/DA=',Q,DQDA
C      CALL SUBBMX(ALPHA,BMX,DBMXDA)
C      WRITE(6,'(A,1P2E12.4)') 'BMAX,DBMAX/DA=',BMX,DBMXDA
C      CALL SUBBMN(ALPHA,BMN,DBMNDA)
C      WRITE(6,'(A,1P2E12.4)') 'BMIN,DBMIN/DA=',BMN,DBMNDA
      CALL SUBMAG(ALPHA,BETA,BR,BZ,BT,BTOT)
      WRITE(6,'(A,1P4E12.4)') 'BR,BZ,BT,BTOT=',BR,BZ,BT,BTOT
      GOTO 1
C
 9000 CALL GSCLOS
      STOP
      END
