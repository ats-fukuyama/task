C     $Id$
C************************************************************************
C  INCOMPLETE LU DECOMPOSITION CONJUGATED GRADIENT SQUARED METHOD FOR  
C  FINITE DIFFERENCE METHOD.(for COMPLEX MATRIX VERSION)
C
C  INPUT  : A(LA,N) : COEFFICIENT MATRIX (complex * 16)
C           B(N)    : RIGHT-HAND-SIDE VECTOR (complex * 16)
C           N       : MATRIX LENGTH
C           M       : DISTANCE FROM CENTER TO REIGHT AND LEFT EDGE
C           LA      : SIZE OF ARRAY A
C  OUTPUT : X(N)    : SOLUTION VECTOR (complex * 16)
C           IERR    ' ERROR CODE
C  WORK   : W(N),P(N),Q(N),R(N),D(N),R0(N),E(N),H(N) (complex * 16)
C
C  COPYRIGHT	TSUTOMU OGUNI      SEP. 1 1991      VER. 1
C  ARRANGED BY TOMOYA TOHNAI   JUN. 16 1998 VER. 3
C***********************************************************************
C
      SUBROUTINE CGSCD(A,N,M,LA,B,EPS,ITR,X,D,P,Q,R,R0,E,H,W,IERR)
C
      COMPLEX * 16 A,B,X,W
      COMPLEX * 16 P,Q,R,D,R0,E,H
      COMPLEX * 16 SS
      COMPLEX * 16 C1,C2,C3
      COMPLEX * 16 Y
      COMPLEX * 16 ALPHA
      COMPLEX * 16 BETA
      REAL    *  8 EPS,RES
      REAL    *  8 X1,X2
C
      DIMENSION A(LA,N),B(N),X(N)
      DIMENSION W (N)
      DIMENSION P (N),Q(N)
      DIMENSION R (N),D(N)
      DIMENSION R0(N),E(N),H(N)
C
      EPS=1.D-8
      ITR=2000
C
      IERR = 0
C
      NRP=0
      DO 12 MS=1,N
         B(MS)=(0.D0,0.D0)
         DO 10 MB=1,2*M+1
            A(MB,MS)=(0.D0,0.D0)
   10    CONTINUE
         CALL WMSETM(A(1,MS),B(MS),MS,LA,NRP)
   12 CONTINUE
C
      DO I=1,N
         X(I)=(0.D0,0.D0)
      ENDDO
C
      DO 15 I=1,N
         D(I) = (0.D0,0.D0)
   15 CONTINUE
C
C  INCOMPLETE LDU DECOMPOSITION
C
      DO 20 I=1,N
         SS=A(M+1,I)
         DO K=I-1,MAX(I-M,1),-1
            J=I-K
            SS=SS-A(M+1-J,I)*A(M+1+J,K)*D(K)
         ENDDO
         D(I) = 1.D0 / SS
   20 CONTINUE
C
      DO 30 I=1,N
         Q(I)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q(I)=Q(I)+A(M+1+J,I)*X(K)
         ENDDO
   30 CONTINUE
C
      DO 40 I=1,N
         R(I) = B(I) - Q(I)
   40 CONTINUE
C
      DO 50 I=1,N
         DO K=I-1,MAX(I-M,1),-1
            J=I-K
            R(I)=R(I)-A(M+1-J,I)*R(K)
         ENDDO
         R(I)=D(I)*R(I)
   50 CONTINUE
C
      DO 60 I=N,1,-1
         DO K=I+1,MIN(I+M,N)
            J=K-I
            R(I)=R(I)-A(M+1+J,I)*R(K)*D(I)
         ENDDO
   60 CONTINUE
C
      C1 = (0.D0,0.D0)
      DO 70 I=1,N
         R0(I) = R(I)
         P(I) = R(I)
         E(I) = R(I)
         C1 = C1 + R(I) * R(I)
   70 CONTINUE
C
      WRITE(6,'(A8,1PE13.4)') 'SC1=',CDABS(C1)
C
C     ITERATION PHASE
C
      KK=0
  200 KK=KK+1
C
      DO 80 I=1,N
         Q(I)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q(I)=Q(I)+A(M+1+J,I)*P(K)
         ENDDO
   80 CONTINUE
C
      DO 90 I=1,N
         DO K=I-1,MAX(I-M,1),-1
            J=I-K
            Q(I)=Q(I)-A(M+1-J,I)*Q(K)
         ENDDO
         Q(I)=D(I)*Q(I)
   90 CONTINUE
C
      DO 100 I=N,1,-1
         DO K=I+1,MIN(I+M,N)
            J=K-I
            Q(I)=Q(I)-D(I)*A(M+1+J,I)*Q(K)
         ENDDO
  100 CONTINUE
C
      C2 = (0.D0,0.D0)
      DO 110 I=1,N
         C2 = C2 + Q(I) * R0(I)
  110 CONTINUE
C
      IF (CDABS(C2) .EQ. 0.D0) THEN
         IERR = 5
         EPS  = RES
         ITR  = KK
         RETURN
      ENDIF
C
      ALPHA = C1 / C2
      C3 = (0.D0,0.D0)
      X1 = 0.D0
      X2 = 0.D0
C
      DO 120 I=1,N
         H(I) = E(I) - ALPHA * Q(I)
  120 CONTINUE
C
      DO 130 I=1,N
         W(I) = E(I) + H(I)
  130 CONTINUE
C
      DO 140 I=1,N
         Q(I)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q(I)=Q(I)+A(M+1+J,I)*W(K)
         ENDDO
  140 CONTINUE
C
      DO 150 I=1,N
         DO K=I-1,MAX(I-M,1),-1
            J=I-K
            Q(I)=Q(I)-A(M+1-J,I)*Q(K)
         ENDDO
         Q(I)=D(I)*Q(I)
  150 CONTINUE
C
      DO 160 I=N,1,-1
         DO K=I+1,MIN(I+M,N)
            J=K-I
            Q(I)=Q(I)-D(I)*A(M+1+J,I)*Q(K)
         ENDDO
  160 CONTINUE
C
      DO 170 I=1,N
         Y = X(I)
         R(I) = R(I) - ALPHA * Q(I)
         X(I) = X(I) + ALPHA * W(I)
         C3 = C3 + R(I) * R0(I)
         X1 = X1 + DCONJG(Y) * Y
         X2 = X2 + DCONJG(X(I) - Y)*(X(I) - Y)
  170 CONTINUE
C
      IF (X1 .NE. 0.D0) THEN
         RES = DSQRT(X2/X1)
C         write(6,*) 'EPS=',kk,res
         IF(RES.LE.EPS) THEN
            ITR  = KK
            IERR = 0
            EPS  = RES
            RETURN
         ENDIF
      ENDIF
C
      BETA = C3 / C1
      C1  = C3
C
      DO 180 I=1,N
         E(I) = R(I) + BETA * H(I)
         P(I) = E(I) + BETA * (H(I) + BETA * P(I))
  180 CONTINUE
C
      IF(KK.LT.ITR) GOTO 200
C
      IERR = 1
      WRITE(6,*) '(SUBR. ILUCGS) NO CONVERGENCE. '
      EPS = RES
C
      RETURN
C  END OF ILUCGS
      END




