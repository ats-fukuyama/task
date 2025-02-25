C     $Id$
C
      CALL TEST1
      CALL TEST2
      STOP
      END
C
      SUBROUTINE TEST1
C
      REAL*8 X1,X2,X,FUNC,FBRENT
      EXTERNAL FUNC
C
      X1=0.D0
      X2=2.D0
C      TOL=1.D-8
      TOL=1.D-15
C      TOL=0.D0
C
      X=FBRENT(FUNC,X1,X2,TOL)
C
      WRITE(6,'(A,1P3E18.10)') 'X,F(X)=',X,FUNC(X)
C
      RETURN
      END
C
      SUBROUTINE TEST2
C
C     THIS IS A SAMPLE PROGRAM FOR THE EASY-TO-USE VERSION BRENT1
C     OF SUBROUTINE BRENTM. THIS PROGRAM SOLVES THE DISCRETE
C     BOUNDARY VALUE PROBLEM DEFINED BY THE SYSTEM OF NONLINEAR
C     EQUATIONS
C
C     2*X(I) - X(I-1) - X(I+1)
C
C            + 0.5*(H**2)*(X(I) + I*H + 1)**3 = 0 , I = 1,...,N
C
C     WHERE H = 1/(N+1), AND X(0) = X(N+1) = 0.
C
C     **********
      INTEGER I,INT,INFO,LWA,N,NFCALL,NFEV,NWRITE
      DOUBLE PRECISION FNORM1,FNORM2,H,TEMP,TOL
      DOUBLE PRECISION X(10),FVEC(10),WA(130)
      DOUBLE PRECISION DFLOAT
      REAL*8 Q(10,10),F(10),Y(10),W1(10),W2(10),EPS
      EXTERNAL BVP,FBVP
      COMMON /REFNUM/ NFCALL
C
C     LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.
C
      DATA NWRITE /6/
C
      DFLOAT(INT) = INT
      LWA = 130
      TOL = 1.D-10
      N = 10
C
C     STARTING VALUES.
C
      H = 1.D0/DFLOAT(N+1)
      DO 10 I = 1, N
         TEMP = DFLOAT(I)*H
         X(I) = TEMP*(TEMP - 1.D0)
   10     CONTINUE
C
C     INITIAL MAX-NORM OF THE RESIDUALS.
C
      FNORM1 = 0.D0
      DO I = 1, N
         CALL BVP(N,X,FVEC,I)
         FNORM1 = DMAX1(FNORM1,DABS(FVEC(I)))
      ENDDO
C
      NFCALL = 0
      CALL FBRENTN(BVP,N,X,FVEC,TOL,INFO,WA,LWA)
      NFEV = NFCALL/N
C
C     FINAL MAX-NORM OF THE RESIDUALS.
C
      FNORM2 = 0.D0
      DO I = 1, N
         CALL BVP(N,X,FVEC,I)
         FNORM2 = DMAX1(FNORM2,DABS(FVEC(I)))
      ENDDO
C
      WRITE (NWRITE,1000) N,FNORM1,FNORM2,NFEV,INFO,(X(I),I=1,N)
C
      RETURN
 1000  FORMAT (5X,10H DIMENSION,I5,5X //
     1     5X,34H INITIAL MAX-NORM OF THE RESIDUALS,D15.7 //
     2     5X,34H FINAL MAX-NORM OF THE RESIDUALS  ,D15.7 //
     3     5X,33H NUMBER OF FUNCTION EVALUATIONS  ,I10 //
     4     5X,15H EXIT PARAMETER ,18X,I10 //
     5     5X,27H FINAL APPROXIMATE SOLUTION  // (5X,5D15.7))
C
C     LAST CARD OF SAMPLE PROGRAM.
C
      END
C
      SUBROUTINE BVP(N,X,FVEC,IFLAG)
      INTEGER N,IFLAG
      DOUBLE PRECISION X(N),FVEC(N)
C     **********
C
C     SUBROUTINE BVP DEFINES THE BOUNDARY VALUE PROBLEM.
C
C     **********
      INTEGER INT,NFCALL
      DOUBLE PRECISION H,TEMP,TEMP1,TEMP2
      DOUBLE PRECISION DFLOAT
      COMMON /REFNUM/ NFCALL
      DFLOAT(INT) = INT
      H = 1.D0/DFLOAT(N+1)
      TEMP = 0.5D0*((X(IFLAG) + DFLOAT(IFLAG)*H + 1.D0)**3)
      TEMP1 = 0.D0
      IF (IFLAG .NE. 1) TEMP1 = X(IFLAG-1)
      TEMP2 = 0.D0
      IF (IFLAG .NE. N) TEMP2 = X(IFLAG+1)
      FVEC(IFLAG) = 2.D0*X(IFLAG) - TEMP1 - TEMP2 + TEMP*H**2
      NFCALL = NFCALL + 1
      RETURN
C
C     LAST CARD OF SUBROUTINE BVP.
C
      END
C
      DOUBLE PRECISION FUNCTION FUNC(X)
C
      REAL*8 X
      FUNC=-0.5D0+SIN(X)
      RETURN
      END
