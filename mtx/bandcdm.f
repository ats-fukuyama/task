C     $Id$
C
C     ****** SOLUTION OF BAND MATRIX (GAUSSIAN ELIMINATION) ******
C
      SUBROUTINE BANDCDM(A,X,N,L,LA,F,NFST,NEND,IERR)
C
      INCLUDE '../mpi/mpif.inc'
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
C
      INTEGER         ISTATUS(MPI_STATUS_SIZE)
      COMPLEX * 16    A(LA,NFST:NEND),X(NFST:NEND),TEMP
      COMPLEX * 16    F(L*(L-1)/2)
      REAL    *  8    EPS , ABS1 , ABS2
      DATA EPS/ 1.D-70 /
C
C      CALL GUTIME(GT1)
      NRP=0
      DO J=ISTA,IEND
         X(J-ISTA+1)=(0.D0,0.D0)
         DO I=1,L
            A(I,J-ISTA+1)=(0.D0,0.D0)
         ENDDO
         CALL WMSETM(A(1,J-ISTA+1),X(J-ISTA+1),J,LA,NRP)
      ENDDO
C      CALL GUTIME(GT2)
C
      IF( MOD(L,2) .EQ. 0 ) GO TO 9000
      LH  = (L+1)/2
      LHM = LH-1
C
      IF(MYRANK.EQ.0) THEN
C
         DO 30 K = 1 , LHM
            LHMK = LH-K
            NPMK = N+1-K
            DO 30 I = 1 , LHMK
               LPMI = L+1-I
               DO 40 J = 2 , L
                  A( J-1 , K-ISTA+1 ) = A( J , K-ISTA+1 )
   40          CONTINUE
               A( L    , K-ISTA+1   ) = ( 0.D0 , 0.D0 )
   30    CONTINUE
C
      ENDIF
C
      IF(MYRANK.EQ.NPROCS-1) THEN
C
         DO K=1,LHM
            LHMK = LH-K
            NPMK = N+1-K
            DO I = 1 , LHMK
               LPMI = L+1-I
               A( LPMI , NPMK-ISTA+1 ) = ( 0.D0 , 0.D0 )
            ENDDO
         ENDDO
C
      ENDIF
C
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      IF(MYRANK.EQ.NPROCS-1) THEN
         IENDD=IEND-1
      ELSE
         IENDD=IEND
      ENDIF
C
      CALL SENDRCVC2(A,F,N,(L-1)/2,LA,NFST,NEND)
      CALL SENDRCVC (X,N,(L-1)/2,NFST,NEND)
C
      DO NRANK=0,NPROCS-1
C
         IF(MYRANK.EQ.NRANK-1) THEN
C
            MN=0
            DO J=IEND+1,IEND+(L-1)/2
               DO I=1,L
                  MN=MN+1
                  F(MN)=A(I,J-ISTA+1)
               ENDDO
            ENDDO
C
            CALL MPI_SEND(F(1),MN,MPI_DOUBLE_COMPLEX,
     &           MYRANK+1,1,MPI_COMM_WORLD,IERR)
            CALL MPI_SEND(X(IEND+1-ISTA+1),(L-1)/2,
     &           MPI_DOUBLE_COMPLEX,MYRANK+1,1,MPI_COMM_WORLD,IERR)
            CALL MPI_SEND(LH,1,MPI_INTEGER,MYRANK+1,1,
     &           MPI_COMM_WORLD,IERR)
C
         ELSE IF(MYRANK.EQ.NRANK) THEN
C
            IF(MYRANK.NE.0) THEN
               CALL MPI_RECV(F(1),L*(L-1)/2,MPI_DOUBLE_COMPLEX,
     &              MYRANK-1,1,MPI_COMM_WORLD,ISTATUS,IERR)
               CALL MPI_RECV(X(1),(L-1)/2,MPI_DOUBLE_COMPLEX,
     &              MYRANK-1,1,MPI_COMM_WORLD,ISTATUS,IERR)
               CALL MPI_RECV(LH,1,MPI_INTEGER,MYRANK-1,1,
     &              MPI_COMM_WORLD,ISTATUS,IERR)
C
               MN=0
               DO J=1,(L-1)/2
                  DO I=1,L
                     MN=MN+1
                     A(I,J)=F(MN)
                  ENDDO
               ENDDO
C
            ENDIF
C
            DO 50 I = ISTA,IENDD
               IPIVOT = I
               IP     = I+1
               ABS2   = CDABS( A(1,IPIVOT-ISTA+1) )
               DO 60 K = IP , LH
                  ABS1 = CDABS( A(1,K-ISTA+1) )
                  IF( ABS1 .GT. ABS2 ) THEN
                     IPIVOT = K
                     ABS2 = ABS1
                  ENDIF
   60          CONTINUE
C
               IF( ABS2 .LT. EPS ) GO TO 9002
               IF( IPIVOT .NE. I ) THEN
                  TEMP             = X( I-ISTA+1 )
                  X(I-ISTA+1     ) = X( IPIVOT-ISTA+1 )
                  X(IPIVOT-ISTA+1) = TEMP
                  DO 90 J = 1 , L
                     TEMP                 = A( J, I-ISTA+1)
                     A( J, I-ISTA+1     ) = A( J, IPIVOT-ISTA+1)
                     A( J, IPIVOT-ISTA+1) = TEMP
   90             CONTINUE
               END IF
C
               TEMP   = 1.D0   / A( 1, I-ISTA+1)
               X( I-ISTA+1 ) = X( I-ISTA+1 ) * TEMP
C
               DO 120 J = 2 , L
                  A( J , I-ISTA+1 ) = A( J , I-ISTA+1 ) * TEMP
  120          CONTINUE
C
               DO 130 K = IP , LH
                  TEMP   = A( 1 , K-ISTA+1 )
                  X(K-ISTA+1)=X(K-ISTA+1)-X(I-ISTA+1)*TEMP
                  DO 140 J = 2 , L
                     A(J-1,K-ISTA+1)=A(J,K-ISTA+1)
     &                              -A(J,I-ISTA+1)*TEMP
  140             CONTINUE
C
                  A(L,K-ISTA+1)=(0.D0,0.D0)
  130          CONTINUE
               IF( LH .LT. N ) LH = LH + 1
   50       CONTINUE
C
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      ENDDO
C
      IF(MYRANK.EQ.NPROCS-1) THEN
         IF( CDABS(A(1,N-ISTA+1)) .LT. EPS ) GO TO 9002
         X(N-ISTA+1) = X(N-ISTA+1)/A(1,N-ISTA+1)
C         WRITE(6,*) 'A(1',N,')',A(1,N-ISTA+1)
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      JJ = 2
C
      DO NRANK=NPROCS-1,0,-1
C
         IF(MYRANK.EQ.NRANK+1) THEN
            CALL MPI_SEND(X(1),L,MPI_DOUBLE_COMPLEX,
     &           MYRANK-1,1,MPI_COMM_WORLD,IERR)
            CALL MPI_SEND(JJ,1,MPI_INTEGER,MYRANK-1,1,
     &           MPI_COMM_WORLD,IERR)
         ELSE IF(MYRANK.EQ.NRANK) THEN
            IF(MYRANK.NE.NPROCS-1) THEN
               CALL MPI_RECV(X(IEND+1-ISTA+1),L,MPI_DOUBLE_COMPLEX,
     &              MYRANK+1,1,MPI_COMM_WORLD,ISTATUS,IERR)
               CALL MPI_RECV(JJ,1,MPI_INTEGER,MYRANK+1,1,
     &              MPI_COMM_WORLD,ISTATUS,IERR)
            ENDIF
C
            DO 160 K = IENDD,ISTA,-1
               TEMP = ( 0.D0 , 0.D0 )
               DO 170 J = 2 , JJ
                  TEMP = A(J,K-ISTA+1)*X(K-1+J-ISTA+1)+TEMP
  170          CONTINUE
               X(K-ISTA+1)=X(K-ISTA+1)-TEMP
               IF( JJ .LT. L ) JJ = JJ + 1
  160       CONTINUE
C
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      ENDDO
C
      DO I=ISTA,IEND
C         X(I-(L-1)/2-ISTA+1)=X(I-ISTA+1)
         X(I-(L+1)/2-ISTA+1)=X(I-ISTA+1)
      ENDDO
C      CALL GUTIME(GT3)
C      WRITE(6,'(A,1P2E12.4)') 'TIME1,TIME2=',GT2-GT1,GT3-GT2
C
      IERR = 0
      RETURN
C
 9000 IERR = 10000
      RETURN
 9002 IERR = 30000 + I
      RETURN
C
      END


