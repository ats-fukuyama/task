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
      SUBROUTINE BCGCDBMA(A,N,M,LA,LD,B,EPS,ITR,X,D,R1,R2,P1,P2,Q1,Q2,
     &                  F,TEMP1,TEMP2,NM,NBSIZ,NFST,NEND,IERR)
C
      INCLUDE '../mpi/mpif.inc'
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
      COMMON /MPIVR3/ NBST,NBED
C
      COMPLEX * 16 A,B,X,D
      COMPLEX * 16 P1,P2,Q1,Q2,R1,R2
      COMPLEX * 16 C1,C2,C3
      COMPLEX * 16 ALPHA
      COMPLEX * 16 BETA
      COMPLEX * 16 Y
      COMPLEX * 16 TEMP1,TEMP2
      COMPLEX * 16 CTMP
      COMPLEX * 16 F
      REAL    *  8 EPS,RES
      REAL    *  8 X1,X2
      REAL    *  8 TMP
C
      DIMENSION A(LA,NFST:NEND),B(NFST:NEND),X(NFST:NEND)
      DIMENSION D(LD,NFST:NEND)
      DIMENSION P1(NFST:NEND),P2(NFST:NEND)
      DIMENSION Q1(NFST:NEND),Q2(NFST:NEND)
      DIMENSION R1(NFST:NEND),R2(NFST:NEND)
      DIMENSION TEMP1(LD,NBSIZ),TEMP2(NBSIZ)
      DIMENSION F(M*(2*M+1))
      DIMENSION NM(NBSIZ)
C
      EPS=1.D-8
      ITR=200
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
      ISCNF=NBSIZ*(NBST-1)+1
      IECNF=NBSIZ*NBED
C
C  INCOMPLETE LDU DECOMPOSITION
C
      CALL SENDRCVC2(A,F,N,M,LA,NFST,NEND)
C
      ICF=2*NBSIZ-NBSIZ
      IDF=2*NBSIZ
      IEF=2*NBSIZ+NBSIZ
C
      DO IBN=NBST,NBED
         IBLC=NBSIZ*(IBN-1)+1
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               TEMP1(I,J)=(0.D0,0.D0)
               DO IC=MAX(IBLC-NBSIZ,ISCNF),IBLC-1
                  K=IC+NBSIZ-IBLC+1
                  TEMP1(I,J)=TEMP1(I,J)
     &                      +A(ICF +(K-1)-(J-1),IBLC+(J-1)-ISTA+1)
     &                      *D(I,IC-ISTA+1)
               ENDDO
            ENDDO
         ENDDO
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               D(I,IBLC+(J-1)-ISTA+1)=(0.D0,0.D0)
               DO IC=MAX(IBLC-NBSIZ,ISCNF),IBLC-1
                  K=IC+NBSIZ-IBLC+1
                  D(I,IBLC+(J-1)-ISTA+1)=D(I,IBLC+(J-1)-ISTA+1)
     &                           +TEMP1(K,J)
     &                           *A(IEF+(I-1)-(K-1),IC-ISTA+1)
               ENDDO
               D(I,IBLC+(J-1)-ISTA+1)
     &                 =A(IDF+(I-1)-(J-1),IBLC+(J-1)-ISTA+1)
     &                 -D(I,IBLC+(J-1)-ISTA+1)
            ENDDO
         ENDDO
C
         CALL INVMCD(D(1,IBLC-ISTA+1),NBSIZ,LD,NM,IERR)
C
      ENDDO
C
 1000 CALL SENDRCVC(X,N,M,NFST,NEND) 
C
      DO I=ISTA,IEND
         Q1(I-ISTA+1)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q1(I-ISTA+1)=Q1(I-ISTA+1)
     &                  +A(M+1+J,I-ISTA+1)*X(K-ISTA+1)
         ENDDO
C            WRITE(20,*) 'I=,Q(I)=',I,Q1(I)
      ENDDO
C
      DO I=ISTA,IEND
         R1(I-ISTA+1)=B(I-ISTA+1)-Q1(I-ISTA+1)
C         WRITE(20,'(A,I3,1PE13.4)') 'I,R=',I,R1(I)
      ENDDO
C
      DO IBN=NBST,NBED
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=R1(IBLC+(I-1)-ISTA+1)
            DO IC=MAX(IBLC-NBSIZ,ISCNF),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)
     &                 -A(ICF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                 *R1(IC-ISTA+1)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            R1(IBLC+(I-1)-ISTA+1)=(0.D0,0.D0)
            DO K=1,NBSIZ
               R1(IBLC+(I-1)-ISTA+1)=R1(IBLC+(I-1)-ISTA+1)
     &                +D(K,IBLC+(I-1)-ISTA+1)*TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO IBN=NBED,NBST,-1
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=(0.D0,0.D0)
            DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,IECNF)
               K=IC-IBLC-NBSIZ+1
               TEMP2(I)=TEMP2(I)
     &                 +A (IEF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                 *R1(IC-ISTA+1)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO K=1,NBSIZ
               R1(IBLC+(I-1)-ISTA+1)=R1(IBLC+(I-1)-ISTA+1)
     &                       -D(K,IBLC+(I-1)-ISTA+1)
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO I=ISTA,IEND
         Q2(I-ISTA+1)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q2(I-ISTA+1)=Q2(I-ISTA+1)
     &                  +A(M+1-J,K-ISTA+1)*DCONJG(X(K-ISTA+1))
         ENDDO
      ENDDO
C
      DO I=ISTA,IEND
         R2(I-ISTA+1)=DCONJG(B(I-ISTA+1))-Q2(I-ISTA+1)
      ENDDO

      DO IBN=NBST,NBED
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=(0.D0,0.D0)
            DO IC=MAX(IBLC-NBSIZ,ISCNF),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)+D (I,IC-ISTA+1)
     &                          *R2(IC-ISTA+1)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO IC=MAX(IBLC-NBSIZ,ISCNF),IBLC-1
               K=IC+NBSIZ-IBLC+1
               R2(IBLC+(I-1)-ISTA+1)=R2(IBLC+(I-1)-ISTA+1)
     &                       -A(IEF+(I-1)-(K-1),IC-ISTA+1)
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO IBN=NBED,NBST,-1
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=R2(IBLC+(I-1)-ISTA+1)
            DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,IECNF)
               K=IC-IBLC-NBSIZ+1
               TEMP2(I)=TEMP2(I)
     &                 -A (ICF+(I-1)-(K-1),IC-ISTA+1)
     &                 *R2(IC-ISTA+1)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            R2(IBLC+(I-1)-ISTA+1)=(0.D0,0.D0)
            DO K=1,NBSIZ
               R2(IBLC+(I-1)-ISTA+1)=R2(IBLC+(I-1)-ISTA+1)
     &                       +D(I,IBLC+(K-1)-ISTA+1)
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      C1 = (0.D0,0.D0)
      DO I=ISTA,IEND
         P1(I-ISTA+1)=R1(I-ISTA+1)
C         P2(I-ISTA+1)=R1(I-ISTA+1)
C         R2(I-ISTA+1)=R1(I-ISTA+1)
         P2(I-ISTA+1)=R2(I-ISTA+1)
         C1=C1+R1(I-ISTA+1)*R2(I-ISTA+1)
      ENDDO
      CALL MPI_REDUCE(C1,CTMP,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      C1=CTMP
      CALL MPI_BCAST(C1,1,MPI_DOUBLE_COMPLEX,0,
     &     MPI_COMM_WORLD,IERR)
C
      IF(MYRANK.EQ.0) WRITE(6,'(A8,1PE13.4)') 'SC1=',CDABS(C1)
C
C  ITERATION PHASE
C
      KK=0
  200 KK=KK+1
C
      DO IBN=NBST,NBED
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            Q1(IBLC+(I-1)-ISTA+1)=P2(IBLC+(I-1)-ISTA+1)
            TEMP2(I)=(0.D0,0.D0)
            DO IC=MAX(IBLC-NBSIZ,ISCNF),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)+D(I,IC-ISTA+1)
     &                          *Q1(IC-ISTA+1)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO IC=MAX(IBLC-NBSIZ,ISCNF),IBLC-1
               K=IC+NBSIZ-IBLC+1
               Q1(IBLC+(I-1)-ISTA+1)=Q1(IBLC+(I-1)-ISTA+1)
     &                       -A(IEF+(I-1)-(K-1),IC-ISTA+1)
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO IBN=NBED,NBST,-1
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=Q1(IBLC+(I-1)-ISTA+1)
            DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,IECNF)
               K=IC-IBLC-NBSIZ+1
               TEMP2(I)=TEMP2(I)
     &                 -A(ICF+(I-1)-(K-1),IC-ISTA+1)
     &                 *Q1(IC-ISTA+1)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            Q1(IBLC+(I-1)-ISTA+1)=(0.D0,0.D0)
            DO K=1,NBSIZ
               Q1(IBLC+(I-1)-ISTA+1)=Q1(IBLC+(I-1)-ISTA+1)
     &                       +D(I,IBLC+(K-1)-ISTA+1)
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      CALL SENDRCVC(Q1,N,M,NFST,NEND)
C
      DO I=ISTA,IEND
         Q2(I-ISTA+1)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q2(I-ISTA+1)=Q2(I-ISTA+1)
     &                  +A(M+1-J,K-ISTA+1)*Q1(K-ISTA+1)
         ENDDO
C     WRITE(22,*) 'I=,Q2=',I,Q2(I)
      ENDDO
C
      CALL SENDRCVC(P1,N,M,NFST,NEND)
C
      DO I=ISTA,IEND
         Q1(I-ISTA+1)=(0.D0,0.D0)
         DO K=MAX(I-M,1),MIN(I+M,N)
            J=K-I
            Q1(I-ISTA+1)=Q1(I-ISTA+1)
     &                  +A(M+1+J,I-ISTA+1)*P1(I+J-ISTA+1)
         ENDDO
C     WRITE(22,*) 'I=,Q1(I)',I,Q1(I)
      ENDDO
C
      DO IBN=NBST,NBED
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=Q1(IBLC+(I-1)-ISTA+1)
            DO IC=MAX(IBLC-NBSIZ,ISCNF),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)
     &                 -A(ICF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                 *Q1(IC-ISTA+1)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            Q1(IBLC+(I-1)-ISTA+1)=(0.D0,0.D0)
            DO K=1,NBSIZ
               Q1(IBLC+(I-1)-ISTA+1)=Q1(IBLC+(I-1)-ISTA+1)
     &                       +D(K,IBLC+(I-1)-ISTA+1)
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      DO IBN=NBED,NBST,-1
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=(0.D0,0.D0)
            DO IC=IBLC+NBSIZ,MIN(IBLC+2*NBSIZ-1,IECNF)
               K=IC-IBLC-NBSIZ+1
               TEMP2(I)=TEMP2(I)
     &                 +A(IEF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                 *Q1(IC-ISTA+1)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO K=1,NBSIZ
               Q1(IBLC+(I-1)-ISTA+1)=Q1(IBLC+(I-1)-ISTA+1)
     &                       -D(K,IBLC+(I-1)-ISTA+1)
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      C2 = (0.D0,0.D0)
      DO 150 I=ISTA,IEND
         C2 = C2+Q1(I-ISTA+1)*P2(I-ISTA+1)
  150 CONTINUE
      CALL MPI_REDUCE(C2,CTMP,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     &     MPI_COMM_WORLD,IERR)
      C2=CTMP
      CALL MPI_BCAST(C2,1,MPI_DOUBLE_COMPLEX,0,
     &     MPI_COMM_WORLD,IERR)
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
      DO 160 I=ISTA,IEND
         Y = X(I-ISTA+1)
         X(I-ISTA+1) = X(I-ISTA+1) + ALPHA *P1(I-ISTA+1)
         R1(I-ISTA+1) = R1(I-ISTA+1) - ALPHA * Q1(I-ISTA+1)
         R2(I-ISTA+1) = R2(I-ISTA+1) - ALPHA * Q2(I-ISTA+1)
         C3 = C3 + R1(I-ISTA+1) * R2(I-ISTA+1)
         X1 = X1 + Y * DCONJG(Y)
         X2 = X2 + (X(I-ISTA+1) - Y)*DCONJG(X(I-ISTA+1) - Y)
  160 CONTINUE
C
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
      IF(X1.NE.0.D0) THEN
         RES = DSQRT(X2/X1)
         if(myrank.eq.0) write(6,*) kk,res
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
      C1   = C3
C
      DO 170 I=ISTA,IEND
         P1(I-ISTA+1) = R1(I-ISTA+1) + BETA * P1(I-ISTA+1)
         P2(I-ISTA+1) = R2(I-ISTA+1) + BETA * P2(I-ISTA+1)
  170 CONTINUE
C
      IF(KK.LE.ITR) GOTO 200
C
C      go to 1000
C
      IERR = 1
      IF(MYRANK.EQ.0) WRITE(6,*) ' (SUBR. ILUBCG) NO CONVERGENCE.'
      EPS = RES
C
      RETURN
C  END OF ILUBCG
      END




