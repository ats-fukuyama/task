C     $Id$
C############################################################
C
C          ブレント法のサブルーチンプログラム
C
C############################################################
C
      SUBROUTINE BRENT(N,X,FUNC,EPS,M,Q,F,Y,W1,W2,IERR)
C
      IMPLICIT REAL*8 (A-B,D-H,O-Z) 
      IMPLICIT COMPLEX*16 (C)
C
      EXTERNAL FUNC
C
      DIMENSION X(N),F(N),Y(N),W1(N),W2(N),Q(N,N)
C
      MMAX = M
      M    = 0
      IC2  = 0
      EPS2 = EPS**2
      IERR = 0
C
      DO 10 I=1,N
         F(I)=FUNC(N,X,I)
   10 CONTINUE
      FNORM=0.D0
      DO 11 I=1,N
         FNORM=FNORM+F(I)**2
   11 CONTINUE
      IF (FNORM.LT.EPS2) RETURN
C
  100 M=M+1
      DO 15 I=1,N
         Y(I)=X(I)
   15 CONTINUE
      DO 50 I=1,N
      DO 50 J=1,N
         Q(I,J)=0.D0
         IF (I.EQ.J) Q(I,J)=Q(I,J)+1
   50 CONTINUE
      DO 60 K=1,N
         FF=FUNC(N,Y,K)
         DELH=0.D0
         DO 61 I=1,N
            DELH=DELH+Y(I)**2
   61    CONTINUE
         DELH=SQRT(DELH)
         IF (DELH.LE.(1.D-8)) DELH=1.D0
         DELH=1.D-8*DELH
         DO 62 I=K,N
            DO 63 J=1,N
               W1(J)=Y(J)+DELH*Q(J,I)
   63       CONTINUE
            W2(I)=(FUNC(N,W1,K)-FF)/DELH
   62    CONTINUE
         ANORM=0.D0
         DO 64 I=K,N
            ANORM=ANORM+W2(I)**2
   64    CONTINUE
         ANORM=SQRT(ANORM)
         IF (ANORM.LT.(1.D-50)) GO TO 9003
         ALPHA=(ANORM*(ANORM+W2(K)))/2.D0
         ALPHA=SQRT(ABS(ALPHA))
         W1(K)=ALPHA/ANORM
         DO 65 I=K+1,N
            W1(I)=W2(I)/(2.D0*ALPHA)
   65    CONTINUE
         SIGMA=W2(K)-2.D0*W1(K)*ALPHA
         DO 66 I=1,K-1
            QSUM=0.D0
            DO 68 L=K,N
               QSUM=QSUM+Q(I,L)*W1(L)
   68       CONTINUE
            DO 67 J=K,N
               Q(I,J)=Q(I,J)-2.D0*W1(J)*QSUM
   67       CONTINUE
   66    CONTINUE
         DO 71 I=K,N
            QSUM=0.D0
            DO 73 L=K,N
               QSUM=QSUM+Q(I,L)*W1(L)
   73       CONTINUE
            DO 72 J=K,N
               Q(I,J)=Q(I,J)-2.D0*W1(J)*QSUM
   72       CONTINUE
   71    CONTINUE
         DO 75 I=1,N
            Y(I)=Y(I)-(FF/SIGMA)*Q(I,K)
   75    CONTINUE
   60 CONTINUE
C
      DO 80 I=1,N
         Y(I)=Y(I)-X(I)
   80 CONTINUE
      T=1
      IC=1
 1000 CONTINUE
         DO 81 I=1,N
            W1(I)=X(I)+T*Y(I)
   81    CONTINUE
         DO 82 I=1,N
            W2(I)=FUNC(N,W1,I)
   82    CONTINUE
         FNORMN=0.D0
         DO 83 I=1,N
            FNORMN=FNORMN+W2(I)**2
   83    CONTINUE
         IF (FNORMN.LT.(((1.D0-T/2.D0)**2)*FNORM)) GO TO 1100
         IF (IC.EQ.30) GO TO 1200
         T=T/2.D0
         IC=IC+1
      GO TO 1000
C
 1100 IC2=0
      GO TO 1300
 1200 IC2=IC2+1
      IF (IC2.GE.3) GO TO 9002
 1300 DO 84 I=1,N
         X(I)=W1(I)
   84 CONTINUE
      FNORM=FNORMN
      IF (FNORM.LT.EPS2) RETURN
      IF (M.GE.MMAX) GO TO 9001
      GO TO 100
C    
 9001 WRITE(6,101)
      IERR=1
      RETURN
 9002 WRITE(6,102)
      IERR=2
      RETURN
 9003 WRITE(6,103)
      IERR=3
      RETURN
  101 FORMAT(1H ,25HFAILED TO CONVERGW(CASE1))
  102 FORMAT(1H ,25HFAILED TO CONVERGW(CASE2))
  103 FORMAT(1H ,25HFAILED TO CONVERGW(CASE3))
      END
C   ***********************************************
C   **               ブレント法                  **
C   ***********************************************
C
      FUNCTION ZBRENT(FUNC,X1,X2,TOL)
C
      INTEGER ITMAX
      REAL*8 ZBRENT,TOL,X1,X2,FUNC,EPS
      EXTERNAL FUNC
      PARAMETER (ITMAX=100,EPS=1.D-15)
      INTEGER ITER
      REAL*8 A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF((FA.GT.0..AND.FB.GT.0.).OR.(FA.LT.0..AND.FB.LT.0.)) THEN
         WRITE(6,'(A,1P3E12.4)') 'XX ZBRENT: ROOT MUST BE BETWEEN',X1,X2
         ZBRENT=A
         RETURN
      ENDIF
      C=B
      FC=FB
      DO 11 ITER=1,ITMAX
        IF((FB.GT.0..AND.FC.GT.0.).OR.(FB.LT.0..AND.FC.LT.0.))THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.D0*EPS*ABS(B)+0.5D0*TOL
        XM=0.5D0*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.D0)THEN
          ZBRENT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.D0*XM*S
            Q=1.D0-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.D0*XM*Q*(Q-R)-(B-A)*(R-1.D0))
            Q=(Q-1.D0)*(R-1.D0)*(S-1.D0)
          ENDIF
          IF(P.GT.0.D0) Q=-Q
          P=ABS(P)
          IF(2.D0*P .LT. MIN(3.D0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      WRITE(6,*) 'ZBRENT EXCEEDING MAXIMUM ITERATIONS'
      ZBRENT=B
      RETURN
      END
