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
      SUBROUTINE CGSCDB(A,N,M,LA,LD,B,EPS,ITR,X,D,P,Q,R,R0,E,H,W,
     &                 TEMP1,TEMP2,NM,NBSIZ,IERR)
C
      COMPLEX * 16 A,B,X,W
      COMPLEX * 16 P,Q,R,D,R0,E,H
      COMPLEX * 16 C1,C2,C3
      COMPLEX * 16 Y
      COMPLEX * 16 ALPHA
      COMPLEX * 16 BETA
      COMPLEX * 16 TEMP1,TEMP2
      REAL    *  8 EPS,RES
      REAL    *  8 X1,X2
C
      DIMENSION A(LA,N),B(N),X(N)
      DIMENSION W(N)
      DIMENSION P(N),Q(N)
      DIMENSION R(N),D(LD,N)
      DIMENSION R0(N),E(N),H(N)
      DIMENSION TEMP1(LD,NBSIZ),TEMP2(NBSIZ)
      DIMENSION NM(NBSIZ)
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
      DO J=1,N
         DO I=1,NBSIZ
            D(I,J)=(0.D0,0.D0)
         ENDDO
      ENDDO
C
C  INCOMPLETE LDU DECOMPOSITION
C
      ICF =2*NBSIZ-NBSIZ
      IDF =2*NBSIZ
      IEF =2*NBSIZ+NBSIZ
C
      DO IBN=1,N/NBSIZ
         IBLC=NBSIZ*(IBN-1)+1
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               TEMP1(I,J)=(0.D0,0.D0)
               DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
                  K=IC+NBSIZ-IBLC+1
                  TEMP1(I,J)=TEMP1(I,J)
     &                      +A(ICF +(K-1)-(J-1),IBLC+(J-1))
     &                      *D(I               ,IC        )
               ENDDO
            ENDDO
         ENDDO
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               D(I,IBLC+(J-1))=(0.D0,0.D0)
               DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
                  K=IC+NBSIZ-IBLC+1
                  D(I,IBLC+(J-1))=D(I,IBLC+(J-1))
     &                           +TEMP1(K,J)
     &                           *A(IEF+(I-1)-(K-1),IC)
               ENDDO
               D(I,IBLC+(J-1))=A(IDF+(I-1)-(J-1),IBLC+(J-1))
     &                        -D(I,IBLC+(J-1))
            ENDDO
         ENDDO
C
         CALL INVMCD(D(1,IBLC),NBSIZ,LD,NM,IERR)
C
      ENDDO
C
      DO 30 I=1,N
         Q(I)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q(I)=Q(I)+A(M+1+J,I)*X(K)
         ENDDO
   30 CONTINUE
      DO 40 I=1,N
         R(I) = B(I) - Q(I)
   40 CONTINUE
C
      DO IBN=1,N/NBSIZ
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ 
            TEMP2(I)=R(IBLC+(I-1))
            DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)-A(ICF+(K-1)-(I-1),IBLC+(I-1))
     &                          *R(IC)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            R(IBLC+(I-1))=(0.D0,0.D0)
            DO K=1,NBSIZ
               R(IBLC+(I-1))=R(IBLC+(I-1))
     &                      +D(K,IBLC+(I-1))*TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO IBN=N/NBSIZ,1,-1
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=(0.D0,0.D0)
            DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,N)
               K=IC-IBLC-NBSIZ+1
               TEMP2(I)=TEMP2(I)+A(IEF+(K-1)-(I-1),IBLC+(I-1))
     &                          *R(IC)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO K=1,NBSIZ
               R(IBLC+(I-1))=R(IBLC+(I-1))
     &                      -D(K,IBLC+(I-1))
     &                      *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      C1 = (0.D0,0.D0)
      DO 70 I=1,N
         R0(I) = R(I)
         P(I) = R(I)
         E(I) = R(I)
         C1 = C1 + R(I) * R(I)
   70 CONTINUE
C
C      WRITE(6,'(A8,1PE13.4)') 'SC1=',CDABS(C1)
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
      DO IBN=1,N/NBSIZ
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=Q(IBLC+(I-1))
            DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)-A(ICF+(K-1)-(I-1),IBLC+(I-1))
     &                          *Q(IBLC-NBSIZ+(K-1))
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            Q(IBLC+(I-1))=(0.D0,0.D0)
            DO K=1,NBSIZ
               Q(IBLC+(I-1))=Q(IBLC+(I-1))
     &                      +D(K,IBLC+(I-1))*TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO IBN=N/NBSIZ,1,-1
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=(0.D0,0.D0)
            DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,N)
               K=IC-IBLC-NBSIZ+1
               TEMP2(I)=TEMP2(I)+A(IEF+(K-1)-(I-1),IBLC+(I-1))
     &                          *Q(IBLC+NBSIZ+(K-1))
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO K=1,NBSIZ
               Q(IBLC+(I-1))=Q(IBLC+(I-1))
     &                      -D(K,IBLC+(I-1))
     &                      *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
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
      DO IBN=1,N/NBSIZ
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=Q(IBLC+(I-1))
            DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)-A(ICF+(K-1)-(I-1),IBLC+(I-1))
     &                          *Q(IBLC-NBSIZ+(K-1))
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            Q(IBLC+(I-1))=(0.D0,0.D0)
            DO K=1,NBSIZ
               Q(IBLC+(I-1))=Q(IBLC+(I-1))
     &                      +D(K,IBLC+(I-1))*TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO IBN=N/NBSIZ,1,-1
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=(0.D0,0.D0)
            DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,N)
               K=IC-IBLC-NBSIZ+1
               TEMP2(I)=TEMP2(I)+A(IEF+(K-1)-(I-1),IBLC+(I-1))
     &                          *Q(IBLC+NBSIZ+(K-1))
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO K=1,NBSIZ
               Q(IBLC+(I-1))=Q(IBLC+(I-1))
     &                      -D(K,IBLC+(I-1))
     &                      *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO 170 I=1,N
         Y = X(I)
         R(I) = R(I) - ALPHA * Q(I)
         X(I) = X(I) + ALPHA * W(I)
         C3 = C3 + R(I) * R0(I)
         X1 = X1 + Y * DCONJG(Y)
         X2 = X2 + (X(I) - Y)*DCONJG(X(I) - Y)
  170 CONTINUE
C
      IF (X1 .NE. 0.D0) THEN
         RES = DSQRT(X2/X1)
C         write(6,*) 'EPS=',kk,res
         IF(RES.LE.EPS) THEN
            ITR  = KK
            IERR = 0
            EPS  = RES
            DO I=1,N
               B(I)=X(I)
            ENDDO
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
      RETURN
C  END OF ILUCGS
      END




