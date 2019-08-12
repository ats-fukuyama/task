C
C     ***** TRI-DIAGONAL MATRIX SOLVER *****
C
      SUBROUTINE TDMSRD(A,RHS,NMAX,X,IERR)
C
C     +++++ INPUT +++++
C        A(3,NMAX) : D : MATRIX COEEFICIENS (CONTENTS TO BE DESTROYED)
C        RHS(NMAX) : D : RIGHT HAND SIDE VECTOR
C             A(1,N)*X(N-1)+A(2,N)*X(N)+A(3,N)*X(N+1)=RHS(N)
C        NMAX      : I : NUMBER OF EQUATIONS
C     +++++ OUTPUT ++++
C        X(NMAX)   : D : SOLUTION VECTOR
C        IERR      : I : ERROR INDICATOR (0 for naormal)
C     +++++++++++++++++
C
      IMPLICIT NONE
      REAL*8 A,RHS,X
      REAL*8 AAN,BBN,AN,BN,CN,DN,FACT,XN
      INTEGER NMAX,IERR,N
      DIMENSION A(3,NMAX),RHS(NMAX),X(NMAX)
C
      AAN=0.D0
      BBN=0.D0
      DO N=1,NMAX
         AN=A(1,N)
         BN=A(2,N)
         CN=A(3,N)
         DN=RHS(N)
         FACT=AN*AAN+BN
         IF(FACT.EQ.0.D0) GOTO 9001
         FACT=1.D0/FACT
         AAN=-CN*FACT
         BBN= (DN-AN*BBN)*FACT
         A(1,N)=AAN
         A(2,N)=BBN
      ENDDO
      XN=BBN
C
      X(NMAX)=XN
      DO N=NMAX-1,1,-1
         AAN=A(1,N)
         BBN=A(2,N)
         X(N)=AAN*X(N+1)+BBN
      ENDDO
      IERR=0
      RETURN
C
 9001 IERR=9001
      RETURN
      END
C
C     ***** PERIODIC TRI-DIAGONAL MATRIX SOLVER *****
C
      SUBROUTINE TDMPRD(A,RHS,NMAX,X,IERR)
C
C     +++++ INPUT +++++
C        A(3,NMAX) : D : MATRIX COEEFICIENS (CONTENTS TO BE DESTROYED)
C        RHS(NMAX) : D : RIGHT HAND SIDE VECTOR
C             A(1,N)*X(N-1)+A(2,N)*X(N)+A(3,N)*X(N+1)=RHS(N)
C        NMAX      : I : NUMBER OF EQUATIONS
C     +++++ OUTPUT ++++
C        X(NMAX)   : D : SOLUTION VECTOR
C        IERR      : I : ERROR INDICATOR (0 for naormal)
C     +++++ COMMENT +++
C        X(0)  =X(N)
C        X(N+1)=X(1)
C     +++++++++++++++++
C
      IMPLICIT NONE
      REAL*8 A,RHS,X
      REAL*8 AAN,BBN,CCN,DDN,EEN,DDNN,EENN
      REAL*8 AN,BN,CN,DN,FACT,X1,XN
      INTEGER NMAX,IERR,N
      DIMENSION A(3,NMAX),RHS(NMAX),X(NMAX)
C
      AAN=0.D0
      BBN=0.D0
      CCN=1.D0
      DO N=1,NMAX
         AN=A(1,N)
         BN=A(2,N)
         CN=A(3,N)
         DN=RHS(N)
         FACT=AN*AAN+BN
         IF(FACT.EQ.0.D0) GOTO 9001
         FACT=1.D0/FACT
         AAN=-CN*FACT
         BBN= (DN-AN*BBN)*FACT
         CCN=-AN*CCN*FACT
         A(1,N)=AAN
         A(2,N)=BBN
         A(3,N)=CCN
      ENDDO
      IF(CCN.EQ.1.D0) GOTO 9002
      DDN=BBN/(1.D0-CCN)
      EEN=AAN/(1.D0-CCN)
      DDNN=DDN
      EENN=EEN
C
      DO N=NMAX-1,1,-1
         AAN=A(1,N)
         BBN=A(2,N)
         CCN=A(3,N)
         DDN=AAN*DDN+BBN+CCN*DDNN
         EEN=AAN*EEN    +CCN*EENN
      ENDDO
      IF(EEN.EQ.1.D0) GOTO 9003
      X1=DDN/(1.D0-EEN)
      XN=DDNN+EENN*X1
C
      X(NMAX)=XN
      DO N=NMAX-1,1,-1
         AAN=A(1,N)
         BBN=A(2,N)
         CCN=A(3,N)
         X(N)=AAN*X(N+1)+BBN+CCN*XN
      ENDDO
      IERR=0
      RETURN
C
 9001 IERR=9001
      RETURN
 9002 IERR=9002
      RETURN
 9003 IERR=9003
      RETURN
      END
C
C     ***** TRI-DIAGONAL MATRIX SOLVER *****
C
      SUBROUTINE TDMSCD(A,RHS,NMAX,X,IERR)
C
C     +++++ INPUT +++++
C        A(3,NMAX) : Z : MATRIX COEEFICIENS (CONTENTS TO BE DESTROYED)
C        RHS(NMAX) : Z : RIGHT HAND SIDE VECTOR
C             A(1,N)*X(N-1)+A(2,N)*X(N)+A(3,N)*X(N+1)=RHS(N)
C        NMAX      : I : NUMBER OF EQUATIONS
C     +++++ OUTPUT ++++
C        X(NMAX)   : D : SOLUTION VECTOR
C        IERR      : I : ERROR INDICATOR (0 for naormal)
C     +++++++++++++++++
C
      IMPLICIT NONE
      COMPLEX*16 A,RHS,X
      COMPLEX*16 AAN,BBN,AN,BN,CN,DN,FACT,XN
      INTEGER NMAX,IERR,N
      DIMENSION A(3,NMAX),RHS(NMAX),X(NMAX)
C
      AAN=0.D0
      BBN=0.D0
      DO N=1,NMAX
         AN=A(1,N)
         BN=A(2,N)
         CN=A(3,N)
         DN=RHS(N)
         FACT=AN*AAN+BN
         IF(FACT.EQ.0.D0) GOTO 9001
         FACT=1.D0/FACT
         AAN=-CN*FACT
         BBN= (DN-AN*BBN)*FACT
         A(1,N)=AAN
         A(2,N)=BBN
      ENDDO
      XN=BBN
C
      X(NMAX)=XN
      DO N=NMAX-1,1,-1
         AAN=A(1,N)
         BBN=A(2,N)
         X(N)=AAN*X(N+1)+BBN
      ENDDO
      IERR=0
      RETURN
C
 9001 IERR=9001
      RETURN
      END
C
C     ***** PERIODIC TRI-DIAGONAL MATRIX SOLVER *****
C
      SUBROUTINE TDMPCD(A,RHS,NMAX,X,IERR)
C
C     +++++ INPUT +++++
C        A(3,NMAX) : Z : MATRIX COEEFICIENS (CONTENTS TO BE DESTROYED)
C        RHS(NMAX) : Z : RIGHT HAND SIDE VECTOR
C             A(1,N)*X(N-1)+A(2,N)*X(N)+A(3,N)*X(N+1)=RHS(N)
C        NMAX      : I : NUMBER OF EQUATIONS
C     +++++ OUTPUT ++++
C        X(NMAX)   : D : SOLUTION VECTOR
C        IERR      : I : ERROR INDICATOR (0 for naormal)
C     +++++ COMMENT +++
C        X(0)  =X(N)
C        X(N+1)=X(1)
C     +++++++++++++++++
C
      IMPLICIT NONE
      COMPLEX*16 A,RHS,X
      COMPLEX*16 AAN,BBN,CCN,DDN,EEN,DDNN,EENN
      COMPLEX*16 AN,BN,CN,DN,FACT,X1,XN
      INTEGER NMAX,IERR,N
      DIMENSION A(3,NMAX),RHS(NMAX),X(NMAX)
C
      AAN=0.D0
      BBN=0.D0
      CCN=1.D0
      DO N=1,NMAX
         AN=A(1,N)
         BN=A(2,N)
         CN=A(3,N)
         DN=RHS(N)
         FACT=AN*AAN+BN
         IF(FACT.EQ.0.D0) GOTO 9001
         FACT=1.D0/FACT
         AAN=-CN*FACT
         BBN= (DN-AN*BBN)*FACT
         CCN=-AN*CCN*FACT
         A(1,N)=AAN
         A(2,N)=BBN
         A(3,N)=CCN
      ENDDO
      IF(CCN.EQ.1.D0) GOTO 9002
      DDN=BBN/(1.D0-CCN)
      EEN=AAN/(1.D0-CCN)
      DDNN=DDN
      EENN=EEN
C
      DO N=NMAX-1,1,-1
         AAN=A(1,N)
         BBN=A(2,N)
         CCN=A(3,N)
         DDN=AAN*DDN+BBN+CCN*DDNN
         EEN=AAN*EEN    +CCN*EENN
      ENDDO
      IF(EEN.EQ.1.D0) GOTO 9003
      X1=DDN/(1.D0-EEN)
      XN=DDNN+EENN*X1
C
      X(NMAX)=XN
      DO N=NMAX-1,1,-1
         AAN=A(1,N)
         BBN=A(2,N)
         CCN=A(3,N)
         X(N)=AAN*X(N+1)+BBN+CCN*XN
      ENDDO
      IERR=0
      RETURN
C
 9001 IERR=9001
      RETURN
 9002 IERR=9002
      RETURN
 9003 IERR=9003
      RETURN
      END
