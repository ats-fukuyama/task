C     $Id$
C############################################################
C
C          BRENT METHOD
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
      DO I=1,N
         F(I)=FUNC(N,X,I)
      ENDDO
      FNORM=0.D0
      DO I=1,N
         FNORM=FNORM+F(I)**2
      ENDDO
      IF (FNORM.LT.EPS2) RETURN
C
  100 M=M+1
      DO I=1,N
         Y(I)=X(I)
      ENDDO
      DO I=1,N
      DO J=1,N
         Q(I,J)=0.D0
         IF (I.EQ.J) Q(I,J)=Q(I,J)+1
      ENDDO
      ENDDO
      DO K=1,N
         FF=FUNC(N,Y,K)
         DELH=0.D0
         DO I=1,N
            DELH=DELH+Y(I)**2
         ENDDO
         DELH=SQRT(DELH)
         IF (DELH.LE.(1.D-8)) DELH=1.D0
         DELH=1.D-8*DELH
         DO I=K,N
            DO J=1,N
               W1(J)=Y(J)+DELH*Q(J,I)
            ENDDO
            W2(I)=(FUNC(N,W1,K)-FF)/DELH
         ENDDO
         ANORM=0.D0
         DO I=K,N
            ANORM=ANORM+W2(I)**2
         ENDDO
         ANORM=SQRT(ANORM)
         IF (ANORM.LT.(1.D-50)) GO TO 9003
         ALPHA=(ANORM*(ANORM+W2(K)))/2.D0
         ALPHA=SQRT(ABS(ALPHA))
         W1(K)=ALPHA/ANORM
         DO I=K+1,N
            W1(I)=W2(I)/(2.D0*ALPHA)
         ENDDO
         SIGMA=W2(K)-2.D0*W1(K)*ALPHA
         DO I=1,K-1
            QSUM=0.D0
            DO L=K,N
               QSUM=QSUM+Q(I,L)*W1(L)
            ENDDO
            DO J=K,N
               Q(I,J)=Q(I,J)-2.D0*W1(J)*QSUM
            ENDDO
         ENDDO
         DO I=K,N
            QSUM=0.D0
            DO L=K,N
               QSUM=QSUM+Q(I,L)*W1(L)
            ENDDO
            DO J=K,N
               Q(I,J)=Q(I,J)-2.D0*W1(J)*QSUM
            ENDDO
         ENDDo
         DO I=1,N
            Y(I)=Y(I)-(FF/SIGMA)*Q(I,K)
         ENDDO
      ENDDO
C
      DO I=1,N
         Y(I)=Y(I)-X(I)
      ENDDO
      T=1
      IC=1
 1000 CONTINUE
         DO I=1,N
            W1(I)=X(I)+T*Y(I)
         ENDDO
         DO I=1,N
            W2(I)=FUNC(N,W1,I)
         ENDDO
         FNORMN=0.D0
         DO I=1,N
            FNORMN=FNORMN+W2(I)**2
         ENDDO
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
 1300 DO I=1,N
         X(I)=W1(I)
      ENDDO
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
  101 FORMAT(' FAILED TO CONVERGE (CASE1)')
  102 FORMAT(' FAILED TO CONVERGE (CASE2)')
  103 FORMAT(' FAILED TO CONVERGE (CASE3)')
      END
C   ***********************************************
C   **               BRENT METHOD                **
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
         WRITE(6,'(A,1P2E12.4)') 'XX ZBRENT: ROOT MUST BE BETWEEN',X1,X2
         WRITE(6,'(A,1P2E12.4)') '         : F(X1),F(X2)=',FA,FB
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
