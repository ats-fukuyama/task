C
      IMPLICIT REAL*8 (A-G,H,O-Z)
      PARAMETER (NM=101)
      DIMENSION A(3,NM),RHS(NM),X(NM)
C
      NMAX=21
      NS=11
      D=0.1D0
C
    1 CONTINUE
      WRITE(6,'(A,1PE12.4,2I5)') 'D,NS,NMAX=',D,NS,NMAX
      READ(5,*,ERR=1,END=9000) D,NS,NMAX
C
      DO N=1,NMAX
         A(1,N)=   -D
         A(2,N)=1+2*D
         A(3,N)=   -D
         RHS(N)=0.D0
      ENDDO
      A(1,1)=0.D0
      A(2,1)=1+D
      A(2,NMAX)=1+D
      A(3,NMAX)=0.D0
      RHS(NS)=1.D0
      CALL TDMSRD(A,RHS,NMAX,X,IERR)
      IF(IERR.NE.0) WRITE(6,'(A,I5)') 'XX IERR=',IERR
C
      WRITE(6,'(1P5E12.4)') (X(N),N=1,NMAX)
      SUM=0.D0
      DO N=1,NMAX
         SUM=SUM+X(N)
      ENDDO
      WRITE(6,'(A,1PE12.4)') 'SUM=',SUM
C
      DO N=1,NMAX
         A(1,N)=   -D
         A(2,N)=1+2*D
         A(3,N)=   -D
         RHS(N)=0.D0
      ENDDO
      RHS(NS)=1.D0
      CALL TDMPRD(A,RHS,NMAX,X,IERR)
      IF(IERR.NE.0) WRITE(6,'(A,I5)') 'XX IERR=',IERR
C
      WRITE(6,'(1P5E12.4)') (X(N),N=1,NMAX)
      SUM=0.D0
      DO N=1,NMAX
         SUM=SUM+X(N)
      ENDDO
      WRITE(6,'(A,1PE12.4)') 'SUM=',SUM
      GOTO 1
C
 9000 CONTINUE
      STOP
      END
      
