C
      IMPLICIT REAL*8 (A-G,H,O-Z)
      PARAMETER (NM=101)
      DIMENSION A(3,NM),RHS(NM),X(NM),A0(3,NM)
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
      DO N =1,NMAX
         RHS(N)=1.D0/N
      ENDDO
C      RHS(NS)=1.D0
C
      DO N=1,NMAX
         A0(1,N)=A(1,N)
         A0(2,N)=A(2,N)
         A0(3,N)=A(3,N)
      ENDDO
C
      CALL TDMSRD(A,RHS,NMAX,X,IERR)
      IF(IERR.NE.0) WRITE(6,'(A,I5)') 'XX IERR=',IERR
C
      N=1
      XX=A0(2,N)*X(N)+A0(3,N)*X(N+1)
      WRITE(6,*) N,RHS(N),XX
      DO N=2,NMAX-1
         XX=A0(1,N)*X(N-1)+A0(2,N)*X(N)+A0(3,N)*X(N+1)
         WRITE(6,*) N,RHS(N),XX
      ENDDO
      N=NMAX
      XX=A0(1,N)*X(N-1)+A0(2,N)*X(N)
      WRITE(6,*) N,RHS(N),XX
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
C
      DO N=1,NMAX
         A0(1,N)=A(1,N)
         A0(2,N)=A(2,N)
         A0(3,N)=A(3,N)
      ENDDO
C
      CALL TDMPRD(A,RHS,NMAX,X,IERR)
      IF(IERR.NE.0) WRITE(6,'(A,I5)') 'XX IERR=',IERR
C
      N=1
      XX=A0(1,N)*X(NMAX)+A0(2,N)*X(N)+A0(3,N)*X(N+1)
      WRITE(6,*) N,RHS(N),XX
      DO N=2,NMAX-1
         XX=A0(1,N)*X(N-1)+A0(2,N)*X(N)+A0(3,N)*X(N+1)
         WRITE(6,*) N,RHS(N),XX
      ENDDO
      N=NMAX
      XX=A0(1,N)*X(N-1)+A0(2,N)*X(N)+A0(3,N)*X(1)
      WRITE(6,*) N,RHS(N),XX
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
      
