C     $Id$
C************************************************************************
C  INCOMPLETE LU DECOMPOSITION CONJUGATED GRADIENT SQUARED METHOD FOR  
C  FINITE DIFFERENCE METHOD.(for COMPLEX MATRIX AND MPI VERSION)
C                           (BLOCK CALCULATION)
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
      SUBROUTINE CGSCDBM(A,N,M,LA,LD,B,EPS,ITR,X,D,P,Q,R,R0,E,H,W,F,
     &                   TEMP1,TEMP2,NM,NBSIZ,NFST,NEND,IERR)
C
      INCLUDE '../mpi/mpif.inc'
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
      COMMON /MPIVR3/ NBST,NBED
C
      COMPLEX * 16 A,B,X,W
      COMPLEX * 16 P,Q,R,D,R0,E,H
      COMPLEX * 16 C1,C2,C3
      COMPLEX * 16 Y
      COMPLEX * 16 ALPHA
      COMPLEX * 16 BETA
      COMPLEX * 16 TEMP1,TEMP2
      COMPLEX * 16 CTMP
      COMPLEX * 16 F
      REAL    *  8 EPS,RES
      REAL    *  8 X1,X2
      REAL    *  8 TMP
      INTEGER      ISTATUS(MPI_STATUS_SIZE)
C
      DIMENSION A(LA,NFST:NEND),B(NFST:NEND),X(NFST:NEND)
      DIMENSION W(NFST:NEND)
      DIMENSION P(NFST:NEND),Q(NFST:NEND)
      DIMENSION R(NFST:NEND),D(LD,NFST:NEND)
      DIMENSION R0(NFST:NEND),E(NFST:NEND),H(NFST:NEND)
      DIMENSION TEMP1(LD,NBSIZ),TEMP2(NBSIZ)
      DIMENSION F(M*(2*M+1))
      DIMENSION NM(NBSIZ)
C
      EPS=1.D-8
      ITR=2000
C
      IERR = 0
C
      NRP=0
      DO 12 MS=ISTA,IEND
         B(MS-ISTA+1)=(0.D0,0.D0)
         DO 10 MB=1,2*M+1
            A(MB,MS-ISTA+1)=(0.D0,0.D0)
   10    CONTINUE
         CALL WMSETM(A(1,MS-ISTA+1),B(MS-ISTA+1),MS,LA,NRP)
   12 CONTINUE
C
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      DO I=ISTA-M,IEND+M
         X(I-ISTA+1)=(0.D0,0.D0)
      ENDDO
C
      DO J=ISTA-NBSIZ,IEND+NBSIZ
         DO I=1,NBSIZ
            D(I,J-ISTA+1)=(0.D0,0.D0)
         ENDDO
      ENDDO
C
C  INCOMPLETE LDU DECOMPOSITION
C
      CALL SENDRCVC2(A,F,N,M,LA,NFST,NEND)
C
      ICF =2*NBSIZ-NBSIZ
      IDF =2*NBSIZ
      IEF =2*NBSIZ+NBSIZ
C
      DO NRANK=0,NPROCS-1
C
         IF(MYRANK.EQ.NRANK-1) THEN
            MN=0
            DO J=1,NBSIZ
               DO I=1,NBSIZ
                  MN=MN+1
                  F(MN)=D(I,IEND-NBSIZ+1+(J-1)-ISTA+1)
               ENDDO
            ENDDO
            CALL MPI_SEND(F(1),MN,MPI_DOUBLE_COMPLEX,
     &           MYRANK+1,1,MPI_COMM_WORLD,IERR)
         ELSE IF(MYRANK.EQ.NRANK) THEN
            IF(MYRANK.NE.0) THEN
               MLTI=NBSIZ*NBSIZ
               CALL MPI_RECV(F(1),MLTI,MPI_DOUBLE_COMPLEX,
     &              MYRANK-1,1,MPI_COMM_WORLD,ISTATUS,IERR)
               MN=0
               DO J=1,NBSIZ
                  DO I=1,NBSIZ
                     MN=MN+1
                     D(I,1-NBSIZ+(J-1))=F(MN)
                  ENDDO
               ENDDO
            ENDIF
C     
            DO IBN=NBST,NBED
               IBLC=NBSIZ*(IBN-1)+1
C
               DO J=1,NBSIZ
                  DO I=1,NBSIZ
                     TEMP1(I,J)=(0.D0,0.D0)
                     DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
                        K=IC+NBSIZ-IBLC+1
                        TEMP1(I,J)=TEMP1(I,J)
     &                       +A(ICF +(K-1)-(J-1),IBLC+(J-1)-ISTA+1)
     &                       *D(I               ,IC        -ISTA+1)
                     ENDDO
                  ENDDO
               ENDDO
C
               DO J=1,NBSIZ
                  DO I=1,NBSIZ
                     D(I,IBLC+(J-1)-ISTA+1)=(0.D0,0.D0)
                     DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
                        K=IC+NBSIZ-IBLC+1
                        D(I,IBLC+(J-1)-ISTA+1)=D(I,IBLC+(J-1)-ISTA+1)
     &                       +TEMP1(K,J)
     &                       *A(IEF+(I-1)-(K-1),IC-ISTA+1)
                     ENDDO
                     D(I,IBLC+(J-1)-ISTA+1)
     &                    =A(IDF+(I-1)-(J-1),IBLC+(J-1)-ISTA+1)
     &                    -D(I,IBLC+(J-1)-ISTA+1)
                  ENDDO
               ENDDO
C
               CALL INVMCD(D(1,IBLC-ISTA+1),NBSIZ,LD,NM,IERR)
C
            ENDDO
C
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      ENDDO
C
      DO 30 I=ISTA,IEND
         Q(I-ISTA+1)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q(I-ISTA+1)=Q(I-ISTA+1)+A(M+1+J,I-ISTA+1)*X(K-ISTA+1)
         ENDDO
   30 CONTINUE
C
      DO 40 I=ISTA,IEND
         R(I-ISTA+1) = B(I-ISTA+1) - Q(I-ISTA+1)
   40 CONTINUE
C
      DO NRANK=0,NPROCS-1
C
         IF(MYRANK.EQ.NRANK-1) THEN
            CALL MPI_SEND(R(IEND-NBSIZ+1-ISTA+1),NBSIZ,
     &           MPI_DOUBLE_COMPLEX,MYRANK+1,1,MPI_COMM_WORLD,IERR)
         ELSE IF(MYRANK.EQ.NRANK) THEN
            IF(MYRANK.NE.0) THEN
               CALL MPI_RECV(R(1-NBSIZ),NBSIZ,MPI_DOUBLE_COMPLEX,
     &              MYRANK-1,1,MPI_COMM_WORLD,ISTATUS,IERR)
            ENDIF
C
            DO IBN=NBST,NBED
               IBLC=NBSIZ*(IBN-1)+1
               DO I=1,NBSIZ 
                  TEMP2(I)=R(IBLC+(I-1)-ISTA+1)
                  DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
                     K=IC+NBSIZ-IBLC+1
                     TEMP2(I)=TEMP2(I)
     &                    -A(ICF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                    *R(IC-ISTA+1)
                  ENDDO
               ENDDO
               DO I=1,NBSIZ
                  R(IBLC+(I-1)-ISTA+1)=(0.D0,0.D0)
                  DO K=1,NBSIZ
                     R(IBLC+(I-1)-ISTA+1)=R(IBLC+(I-1)-ISTA+1)
     &                    +D(K,IBLC+(I-1)-ISTA+1)*TEMP2(K)
                  ENDDO
               ENDDO
            ENDDO
C
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      ENDDO
C
      DO NRANK=NPROCS-1,0,-1
C
         IF(MYRANK.EQ.NRANK+1) THEN
            CALL MPI_SEND(R(1),NBSIZ,MPI_DOUBLE_COMPLEX,
     &           MYRANK-1,1,MPI_COMM_WORLD,IERR)
         ELSE IF(MYRANK.EQ.NRANK) THEN
            IF(MYRANK.NE.NPROCS-1) THEN
               CALL MPI_RECV(R(IEND+1-ISTA+1),NBSIZ,
     &              MPI_DOUBLE_COMPLEX,MYRANK+1,1,
     &              MPI_COMM_WORLD,ISTATUS,IERR)
            ENDIF
C
            DO IBN=NBED,NBST,-1
               IBLC=NBSIZ*(IBN-1)+1
               DO I=1,NBSIZ
                  TEMP2(I)=(0.D0,0.D0)
                  DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,N)
                     K=IC-IBLC-NBSIZ+1
                     TEMP2(I)=TEMP2(I)
     &                    +A(IEF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                    *R(IC-ISTA+1)
                  ENDDO
               ENDDO
               DO I=1,NBSIZ
                  DO K=1,NBSIZ
                     R(IBLC+(I-1)-ISTA+1)=R(IBLC+(I-1)-ISTA+1)
     &                    -D(K,IBLC+(I-1)-ISTA+1)
     &                    *TEMP2(K)
                  ENDDO
               ENDDO
            ENDDO
C
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      ENDDO
C
      C1 = (0.D0,0.D0)
      DO 70 I=ISTA,IEND
         R0(I-ISTA+1) = R(I-ISTA+1)
         P (I-ISTA+1) = R(I-ISTA+1)
         E (I-ISTA+1) = R(I-ISTA+1)
         C1 = C1 + R(I-ISTA+1) * R(I-ISTA+1)
   70 CONTINUE
      CALL MPI_REDUCE(C1,CTMP,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      C1=CTMP
      CALL MPI_BCAST(C1,1,MPI_DOUBLE_COMPLEX,0,
     &     MPI_COMM_WORLD,IERR)
C
C      IF(MYRANK.EQ.0) WRITE(6,'(A8,1PE13.4)') 'SC1=',CDABS(C1)
C
C     ITERATION PHASE
C
      KK=0
  200 KK=KK+1
C
      CALL SENDRCVC(P,N,M,NFST,NEND)
C
      DO 80 I=ISTA,IEND
         Q(I-ISTA+1)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q(I-ISTA+1)=Q(I-ISTA+1)+A(M+1+J,I-ISTA+1)*P(K-ISTA+1)
         ENDDO
   80 CONTINUE
C
      DO NRANK=0,NPROCS-1
C
         IF(MYRANK.EQ.NRANK-1) THEN
            CALL MPI_SEND(Q(IEND-NBSIZ+1-ISTA+1),NBSIZ,
     &           MPI_DOUBLE_COMPLEX,MYRANK+1,1,MPI_COMM_WORLD,IERR)
         ELSE IF(MYRANK.EQ.NRANK) THEN
            IF(MYRANK.NE.0) THEN
               CALL MPI_RECV(Q(1-NBSIZ),NBSIZ,MPI_DOUBLE_COMPLEX,
     &              MYRANK-1,1,MPI_COMM_WORLD,ISTATUS,IERR)
            ENDIF
C
            DO IBN=NBST,NBED
               IBLC=NBSIZ*(IBN-1)+1
               DO I=1,NBSIZ
                  TEMP2(I)=Q(IBLC+(I-1)-ISTA+1)
                  DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
                     K=IC+NBSIZ-IBLC+1
                     TEMP2(I)=TEMP2(I)
     &                    -A(ICF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                    *Q(IBLC-NBSIZ+(K-1)-ISTA+1)
                  ENDDO
               ENDDO
               DO I=1,NBSIZ
                  Q(IBLC+(I-1)-ISTA+1)=(0.D0,0.D0)
                  DO K=1,NBSIZ
                     Q(IBLC+(I-1)-ISTA+1)=Q(IBLC+(I-1)-ISTA+1)
     &                    +D(K,IBLC+(I-1)-ISTA+1)*TEMP2(K)
                  ENDDO
               ENDDO
            ENDDO
C
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      ENDDO
C
      DO NRANK=NPROCS-1,0,-1
C
         IF(MYRANK.EQ.NRANK+1) THEN
            CALL MPI_SEND(Q(1),NBSIZ,MPI_DOUBLE_COMPLEX,
     &           MYRANK-1,1,MPI_COMM_WORLD,IERR)
         ELSE IF(MYRANK.EQ.NRANK) THEN
            IF(MYRANK.NE.NPROCS-1) THEN
               CALL MPI_RECV(Q(IEND+1-ISTA+1),NBSIZ,
     &              MPI_DOUBLE_COMPLEX,MYRANK+1,1,
     &              MPI_COMM_WORLD,ISTATUS,IERR)
            ENDIF
C
            DO IBN=NBED,NBST,-1
               IBLC=NBSIZ*(IBN-1)+1
               DO I=1,NBSIZ
                  TEMP2(I)=(0.D0,0.D0)
                  DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,N)
                     K=IC-IBLC-NBSIZ+1
                     TEMP2(I)=TEMP2(I)
     &                    +A(IEF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                    *Q(IBLC+NBSIZ+(K-1)-ISTA+1)
                  ENDDO
               ENDDO
               DO I=1,NBSIZ
                  DO K=1,NBSIZ
                     Q(IBLC+(I-1)-ISTA+1)=Q(IBLC+(I-1)-ISTA+1)
     &                    -D(K,IBLC+(I-1)-ISTA+1)
     &                    *TEMP2(K)
                  ENDDO
               ENDDO
            ENDDO
C
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      ENDDO
C
      C2 = (0.D0,0.D0)
      DO 110 I=ISTA,IEND
         C2 = C2 + Q(I-ISTA+1) * R0(I-ISTA+1)
  110 CONTINUE
      CALL MPI_REDUCE(C2,CTMP,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      C2=CTMP
      CALL MPI_BCAST(C2,1,MPI_DOUBLE_COMPLEX,0,
     &     MPI_COMM_WORLD,IERR)
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
      DO 120 I=ISTA,IEND
         H(I-ISTA+1) = E(I-ISTA+1) - ALPHA * Q(I-ISTA+1)
  120 CONTINUE
C
      DO 130 I=ISTA,IEND
         W(I-ISTA+1) = E(I-ISTA+1) + H(I-ISTA+1)
  130 CONTINUE
C
      CALL SENDRCVC(W,N,M,NFST,NEND)
C
      DO 140 I=ISTA,IEND
         Q(I-ISTA+1)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q(I-ISTA+1)=Q(I-ISTA+1)+A(M+1+J,I-ISTA+1)*W(K-ISTA+1)
         ENDDO
  140 CONTINUE
C
      DO NRANK=0,NPROCS-1
C
         IF(MYRANK.EQ.NRANK-1) THEN
            CALL MPI_SEND(Q(IEND-NBSIZ+1-ISTA+1),NBSIZ,
     &           MPI_DOUBLE_COMPLEX,MYRANK+1,1,MPI_COMM_WORLD,IERR)
         ELSE IF(MYRANK.EQ.NRANK) THEN
            IF(MYRANK.NE.0) THEN
               CALL MPI_RECV(Q(1-NBSIZ),NBSIZ,MPI_DOUBLE_COMPLEX,
     &              MYRANK-1,1,MPI_COMM_WORLD,ISTATUS,IERR)
            ENDIF
C
            DO IBN=NBST,NBED
               IBLC=NBSIZ*(IBN-1)+1
               DO I=1,NBSIZ
                  TEMP2(I)=Q(IBLC+(I-1)-ISTA+1)
                  DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
                     K=IC+NBSIZ-IBLC+1
                     TEMP2(I)=TEMP2(I)
     &                    -A(ICF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                    *Q(IBLC-NBSIZ+(K-1)-ISTA+1)
                  ENDDO
               ENDDO
               DO I=1,NBSIZ
                  Q(IBLC+(I-1)-ISTA+1)=(0.D0,0.D0)
                  DO K=1,NBSIZ
                     Q(IBLC+(I-1)-ISTA+1)=Q(IBLC+(I-1)-ISTA+1)
     &                    +D(K,IBLC+(I-1)-ISTA+1)*TEMP2(K)
                  ENDDO
               ENDDO
            ENDDO
C
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      ENDDO
C
      DO NRANK=NPROCS-1,0,-1
C
         IF(MYRANK.EQ.NRANK+1) THEN
            CALL MPI_SEND(Q(1),NBSIZ,MPI_DOUBLE_COMPLEX,
     &           MYRANK-1,1,MPI_COMM_WORLD,IERR)
         ELSE IF(MYRANK.EQ.NRANK) THEN
            IF(MYRANK.NE.NPROCS-1) THEN
               CALL MPI_RECV(Q(IEND+1-ISTA+1),NBSIZ,
     &              MPI_DOUBLE_COMPLEX,MYRANK+1,1,
     &              MPI_COMM_WORLD,ISTATUS,IERR)
            ENDIF
C
            DO IBN=NBED,NBST,-1
               IBLC=NBSIZ*(IBN-1)+1
               DO I=1,NBSIZ
                  TEMP2(I)=(0.D0,0.D0)
                  DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,N)
                     K=IC-IBLC-NBSIZ+1
                     TEMP2(I)=TEMP2(I)
     &                    +A(IEF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                    *Q(IBLC+NBSIZ+(K-1)-ISTA+1)
                  ENDDO
               ENDDO
               DO I=1,NBSIZ
                  DO K=1,NBSIZ
                     Q(IBLC+(I-1)-ISTA+1)=Q(IBLC+(I-1)-ISTA+1)
     &                    -D(K,IBLC+(I-1)-ISTA+1)
     &                    *TEMP2(K)
                  ENDDO
               ENDDO
            ENDDO
C
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      ENDDO
C
      DO 170 I=ISTA,IEND
         Y = X(I-ISTA+1)
         R(I-ISTA+1) = R(I-ISTA+1)-ALPHA*Q(I-ISTA+1)
         X(I-ISTA+1) = X(I-ISTA+1)+ALPHA*W(I-ISTA+1)
         C3 = C3 + R(I-ISTA+1) * R0(I-ISTA+1)
         X1 = X1 + Y * DCONJG(Y)
         X2 = X2 + (X(I-ISTA+1) - Y)*DCONJG(X(I-ISTA+1) - Y)
  170 CONTINUE
      CALL MPI_REDUCE(C3,CTMP,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      C3=CTMP
      CALL MPI_BCAST(C3,1,MPI_DOUBLE_COMPLEX,0,
     &     MPI_COMM_WORLD,IERR)      
C
      CALL MPI_REDUCE(X1,TMP,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      X1=TMP
      CALL MPI_BCAST(X1,1,MPI_DOUBLE_PRECISION,0,
     &     MPI_COMM_WORLD,IERR)
C
      CALL MPI_REDUCE(X2,TMP,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      X2=TMP
      CALL MPI_BCAST(X2,1,MPI_DOUBLE_PRECISION,0,
     &     MPI_COMM_WORLD,IERR)
C
      IF (X1 .NE. 0.D0) THEN
         RES = DSQRT(X2/X1)
C         IF(MYRANK.EQ.0) write(6,*) 'EPS=',kk,res
         IF(RES.LE.EPS) THEN
            ITR  = KK
            IERR = 0
            EPS  = RES
            DO I=ISTA,IEND
C               B(I-M-ISTA+1)=X(I-ISTA+1)
                B(I-2*NBSIZ-ISTA+1)=X(I-ISTA+1)
            ENDDO
            RETURN
         ENDIF
      ENDIF
C
      BETA = C3 / C1
      C1  = C3
C
      DO 180 I=ISTA,IEND
         E(I-ISTA+1)=R(I-ISTA+1)+BETA*H(I-ISTA+1)
         P(I-ISTA+1)=E(I-ISTA+1)
     &              +BETA*(H(I-ISTA+1)+BETA*P(I-ISTA+1))
  180 CONTINUE
C
      IF(KK.LT.ITR) GOTO 200
C
      IERR = 1
      IF(MYRANK.EQ.0) WRITE(6,*) '(SUBR. ILUCGS) NO CONVERGENCE. '
      EPS = RES
      RETURN
C  END OF ILUCGS
      END




