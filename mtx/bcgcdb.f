C     $Id$
C**********************************************************************
C  BI-CONJUGATE GRADIENT METHOD WITH INCOMPLETE LU DECOMPOSITION.     
C  MATRIX FOR COMPLEX
C  COPYRIGHT:     TSUTOMU OGUNI       SEP. 1 1991               
C  ARRANGED BY TOMOYA TOHNAI    JUN. 19 1998 VER.1
C                                                                   
C  INPUT  : A(LA,N) : COEFFICIENT MATRIX (COMPLEX)
C           B(N)    : RIGHT-HAND-SIDE VECTOR (COMPLEX) 
C           N       : MATRIX SIZE1
C           M       : DISTANCE FROM CENTER TO RIGHT AND LEFT SIDE
C  OUTPUT : X(N)    : SOLUTION VECTOR (COMPLEX)
C           IERR    : ERROR CODE
C    
C  WORK   : R1(N),R2(N),P1(N),P2(N),Q1(N),Q2(N)
C
C**********************************************************************
C
      SUBROUTINE BCGCDB(A,N,M,LA,LD,B,EPS,ITR,X,D,R1,R2,P1,P2,Q1,Q2,
     &                  TEMP1,TEMP2,NM,NBSIZ,IERR)
C
      COMPLEX * 16 A,B,X,D
      COMPLEX * 16 P1,P2,Q1,Q2,R1,R2
      COMPLEX * 16 C1,C2,C3
      COMPLEX * 16 ALPHA
      COMPLEX * 16 BETA
      COMPLEX * 16 Y
      COMPLEX * 16 TEMP1,TEMP2
      REAL    *  8 X1,X2
      REAL    *  8 EPS
      REAL    *  8 RES
C
      DIMENSION A(LA,N),B(N),X(N)
      DIMENSION D(LD,N)
      DIMENSION P1(N),P2(N)
      DIMENSION Q1(N),Q2(N)
      DIMENSION R1(N),R2(N)
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
      ICF=2*NBSIZ-NBSIZ
      IDF=2*NBSIZ
      IEF=2*NBSIZ+NBSIZ
C      do ibn=1,n/nbsiz
C         iblc=nbsiz*(ibn-1)+1
C         do j=1,nbsiz
C            do i=1,icf-j
C               write(21,*) i,iblc+j-1,a(i,iblc+j-1)
C            enddo
C         enddo
C      enddo
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
     &                      *D(I,IC)
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
      DO I=1,N
         Q1(I)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q1(I)=Q1(I)+A(M+1+J,I)*X(K)
         ENDDO
C            WRITE(20,*) 'I=,Q(I)=',I,Q1(I)
      ENDDO
C
      DO I=1,N
         R1(I)=B(I)-Q1(I)
C         WRITE(20,'(A,I3,1PE13.4)') 'I,R=',I,R1(I)
      ENDDO
C
      DO IBN=1,N/NBSIZ
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=R1(IBLC+(I-1))
            DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)-A(ICF+(K-1)-(I-1),IBLC+(I-1))
     &                          *R1(IC)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            R1(IBLC+(I-1))=(0.D0,0.D0)
            DO K=1,NBSIZ
               R1(IBLC+(I-1))=R1(IBLC+(I-1))
     &                       +D(K,IBLC+(I-1))*TEMP2(K)
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
               TEMP2(I)=TEMP2(I)+A (IEF+(K-1)-(I-1),IBLC+(I-1))
     &                          *R1(IC)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO K=1,NBSIZ
               R1(IBLC+(I-1))=R1(IBLC+(I-1))
     &                       -D(K,IBLC+(I-1))
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
            
C
      DO I=1,N
         Q2(I)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q2(I)=Q2(I)+A(M+1-J,K)*DCONJG(X(K))
         ENDDO
      ENDDO
C
      DO I=1,N
         R2(I)=DCONJG(B(I))-Q2(I)
      ENDDO
C
      DO IBN=1,N/NBSIZ
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=(0.D0,0.D0)
            DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)+D (I,IC)
     &                          *R2(IC)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
               K=IC+NBSIZ-IBLC+1
               R2(IBLC+(I-1))=R2(IBLC+(I-1))
     &                       -A(IEF+(I-1)-(K-1),IC)
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO IBN=N/NBSIZ,1,-1
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=R2(IBLC+(I-1))
            DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,N)
               K=IC-IBLC-NBSIZ+1
               TEMP2(I)=TEMP2(I)
     &                 -A (ICF+(I-1)-(K-1),IC)
     &                 *R2(IC)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            R2(IBLC+(I-1))=(0.D0,0.D0)
            DO K=1,NBSIZ
               R2(IBLC+(I-1))=R2(IBLC+(I-1))
     &                       +D(I,IBLC+(K-1))
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      C1 = (0.D0,0.D0)
      DO I=1,N
         P1(I)=R1(I)
         P2(I)=R2(I)
         C1=C1+R1(I)*R2(I)
      ENDDO
C      WRITE(6,'(A8,1PE13.4)') 'SC1=',CDABS(C1)
C
C  ITERATION PHASE
C
      KK=0
  200 KK=KK+1
C
      DO IBN=1,N/NBSIZ
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            Q1(IBLC+(I-1))=P2(IBLC+(I-1))
            TEMP2(I)=(0.D0,0.D0)
            DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)+D(I,IC)
     &                          *Q1(IC)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
               K=IC+NBSIZ-IBLC+1
               Q1(IBLC+(I-1))=Q1(IBLC+(I-1))
     &                       -A(IEF+(I-1)-(K-1),IC)
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO IBN=N/NBSIZ,1,-1
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=Q1(IBLC+(I-1))
            DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,N)
               K=IC-IBLC-NBSIZ+1
               TEMP2(I)=TEMP2(I)
     &                 -A(ICF+(I-1)-(K-1),IC)
     &                 *Q1(IC)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            Q1(IBLC+(I-1))=(0.D0,0.D0)
            DO K=1,NBSIZ
               Q1(IBLC+(I-1))=Q1(IBLC+(I-1))
     &                       +D(I,IBLC+(K-1))
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO I=1,N
         Q2(I)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q2(I)=Q2(I)+A(M+1-J,K)*Q1(K)
         ENDDO
C     WRITE(22,*) 'I=,Q2=',I,Q2(I)
      ENDDO
C
      DO I=1,N
         Q1(I)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q1(I)=Q1(I)+A(M+1+J,I)*P1(I+J)
         ENDDO
C     WRITE(22,*) 'I=,Q1(I)',I,Q1(I)
      ENDDO
C
      DO IBN=1,N/NBSIZ
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=Q1(IBLC+(I-1))
            DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)
     &                 -A(ICF+(K-1)-(I-1),IBLC+(I-1))
     &                 *Q1(IC)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            Q1(IBLC+(I-1))=(0.D0,0.D0)
            DO K=1,NBSIZ
               Q1(IBLC+(I-1))=Q1(IBLC+(I-1))
     &                       +D(K,IBLC+(I-1))
     &                       *TEMP2(K)
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
     &                          *Q1(IC)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO K=1,NBSIZ
               Q1(IBLC+(I-1))=Q1(IBLC+(I-1))
     &                       -D(K,IBLC+(I-1))
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      C2 = (0.D0,0.D0)
      DO 150 I=1,N
         C2 = C2+Q1(I)*P2(I)
  150 CONTINUE
C
      IF(CDABS(C2).EQ.0.D0) THEN
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
      DO 160 I=1,N
         Y = X(I)
         X(I) = X(I) + ALPHA *P1(I)
         R1(I) = R1(I) - ALPHA * Q1(I)
         R2(I) = R2(I) - ALPHA * Q2(I)
         C3 = C3 + R1(I) * R2(I)
         X1 = X1 + Y * DCONJG(Y)
         X2 = X2 + (X(I) - Y)*DCONJG(X(I) - Y)
  160 CONTINUE
C
      IF(X1.NE.0.D0) THEN
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
      C1   = C3
C
      DO 170 I=1,N
         P1(I) = R1(I) + BETA * P1(I)
         P2(I) = R2(I) + BETA * P2(I)
  170 CONTINUE
C
      IF(KK.LE.ITR) GOTO 200
C
      IERR = 1
      WRITE(6,*) ' (SUBR. ILUBCG) NO CONVERGENCE.'
      EPS = RES
C
      RETURN
C  END OF ILUBCG
      END




