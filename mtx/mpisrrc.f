C     $Id$
C
C     ### SEND AND RECIEVE ROUTINE FOR COMPLEX MATRIX AND VECTOR
C         IN BLOCK DIVISION ###
C
C     ######## SEND AND RECIEVE ROUTINE ########
C     
      SUBROUTINE SENDRCVC(B,NMAX,M,NFST,NEND)
C     
      INCLUDE 'mpif.h'
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
C
      INTEGER ISTATUS(MPI_STATUS_SIZE)
C
      COMPLEX * 16 B
      DIMENSION B(NFST:NEND)
C
      ISUB=M
      IF(MOD(MYRANK,2).EQ.0) THEN
C
         IF(MYRANK.NE.NPROCS-1) THEN
            ISTAS=IEND-M+1-ISTA+1
            CALL MPI_SEND(B(ISTAS),ISUB,MPI_DOUBLE_COMPLEX,
     &           MYRANK+1,1,MPI_COMM_WORLD,IERR)
         ENDIF
C
         IF(MYRANK.NE.0) THEN
C
            ISTAR=1-M
            CALL MPI_RECV(B(ISTAR),ISUB,MPI_DOUBLE_COMPLEX,
     &           MYRANK-1,1,MPI_COMM_WORLD,ISTATUS,IERR)
C
            IF(MYRANK.EQ.NPROCS-1.AND.IEND-ISTA+1.LT.M) THEN
               ISUB=IEND-ISTA+1
            ENDIF
            CALL MPI_SEND(B(1),ISUB,MPI_DOUBLE_COMPLEX,
     &           MYRANK-1,1,MPI_COMM_WORLD,IERR)
         ENDIF
C
         IF(MYRANK.NE.NPROCS-1) THEN
C            ISTAR=(MYRANK+1)*IEND-ISTA+1+1
            ISTAR=IEND+1-ISTA+1
            IF(MYRANK.EQ.NPROCS-2.AND.IEND-ISTA+1.LT.M) THEN
               ISUB=NMAX-ISTAR+1
            ENDIF
            CALL MPI_RECV(B(ISTAR),ISUB,MPI_DOUBLE_COMPLEX,
     &           MYRANK+1,1,MPI_COMM_WORLD,ISTATUS,IERR)
         ENDIF
C
      ENDIF
C
      IF(MOD(MYRANK,2).EQ.1) THEN
C
         ISTAR=1-M
         CALL MPI_RECV(B(ISTAR),ISUB,MPI_DOUBLE_COMPLEX,
     &        MYRANK-1,1,MPI_COMM_WORLD,ISTATUS,IERR)
C
         ISTAS=IEND-M+1-ISTA+1
         IF(MYRANK.NE.NPROCS-1) THEN
C
            CALL MPI_SEND(B(ISTAS),ISUB,MPI_DOUBLE_COMPLEX,
     &           MYRANK+1,1,MPI_COMM_WORLD,IERR)
C
            ISTAR=IEND+1-ISTA+1
            IF(MYRANK.EQ.NPROCS-2.AND.IEND-ISTA+1.LT.M) THEN
               ISUB=NMAX-ISTAR+1
            ENDIF
            CALL MPI_RECV(B(ISTAR),ISUB,MPI_DOUBLE_COMPLEX,
     &           MYRANK+1,1,MPI_COMM_WORLD,ISTATUS,IERR)
C
         ENDIF
C
         IF(MYRANK.EQ.NPROCS-1.AND.IEND-ISTA+1.LT.M) THEN
            ISUB=IEND-ISTA+1
         ENDIF
         CALL MPI_SEND(B(1),ISUB,MPI_DOUBLE_COMPLEX,
     &           MYRANK-1,1,MPI_COMM_WORLD,IERR)
C
      ENDIF
C
      RETURN
      END
C
C     ######## SEND AND RECIEVE ROUTINE ########
C
      SUBROUTINE SENDRCVC2(SR,F,NMAX,M,LA,NFST,NEND)
C
      INCLUDE 'mpif.h'
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
C
      INTEGER ISTATUS(MPI_STATUS_SIZE)
C
      COMPLEX * 16 SR,F
      DIMENSION SR(LA,NFST:NEND),F(M*(2*M+1))
C
      ISUB=M
      IF(MOD(MYRANK,2).EQ.0) THEN
C     
         IF(MYRANK.NE.NPROCS-1) THEN
            ISTAS=IEND-M+1
            MN=0
            DO I=ISTAS,IEND
               DO J=1,2*M+1
                  MN=MN+1
                  F(MN)=SR(J,I-ISTA+1)
               ENDDO
            ENDDO
            CALL MPI_SEND(F(1),MN,MPI_DOUBLE_COMPLEX,
     &           MYRANK+1,1,MPI_COMM_WORLD,IERR)
         ENDIF
C     
         IF(MYRANK.NE.0) THEN
C     
            ISTAR=ISTA-M
            MN=M*(2*M+1)
            CALL MPI_RECV(F(1),MN,MPI_DOUBLE_COMPLEX,
     &           MYRANK-1,1,MPI_COMM_WORLD,ISTATUS,IERR)
            MN=0
            DO I=ISTAR,ISTA-1
               DO J=1,2*M+1
                  MN=MN+1
                  SR(J,I-ISTA+1)=F(MN)
               ENDDO
            ENDDO
C
            IF(MYRANK.EQ.NPROCS-1.AND.IEND-ISTA+1.LT.M) THEN
               ISUB=IEND-ISTA+1
            ENDIF
            MN=0
            DO I=ISTA,ISTA+ISUB-1
               DO J=1,2*M+1
                  MN=MN+1
                  F(MN)=SR(J,I-ISTA+1)
               ENDDO
            ENDDO
            CALL MPI_SEND(F(1),MN,MPI_DOUBLE_COMPLEX,
     &           MYRANK-1,1,MPI_COMM_WORLD,IERR)
         ENDIF
C
         IF(MYRANK.NE.NPROCS-1) THEN
C     
            ISTAR=IEND+1
            IF(MYRANK.EQ.NPROCS-2.AND.IEND-ISTA+1.LT.M) THEN
               ISUB=NMAX-ISTAR+1
            ENDIF
            MN=ISUB*(2*M+1)
            CALL MPI_RECV(F(1),MN,MPI_DOUBLE_COMPLEX,
     &           MYRANK+1,1,MPI_COMM_WORLD,ISTATUS,IERR)
            MN=0
            DO I=ISTAR,ISTAR+ISUB-1
               DO J=1,2*M+1
                  MN=MN+1
                  SR(J,I-ISTA+1)=F(MN)
               ENDDO
            ENDDO
C     
         ENDIF
C
      ENDIF
C
      IF(MOD(MYRANK,2).EQ.1) THEN
C
         ISTAR=ISTA-M
         MN=M*(2*M+1)
         CALL MPI_RECV(F(1),MN,MPI_DOUBLE_COMPLEX,
     &        MYRANK-1,1,MPI_COMM_WORLD,ISTATUS,IERR)
         MN=0
         DO I=ISTAR,ISTA-1
            DO J=1,2*M+1
               MN=MN+1
               SR(J,I-ISTA+1)=F(MN)
            ENDDO
         ENDDO
C     
         ISTAS=IEND-M+1
         IF(MYRANK.NE.NPROCS-1) THEN
C     
            MN=0
            DO I=ISTAS,IEND
               DO J=1,2*M+1
                  MN=MN+1
                  F(MN)=SR(J,I-ISTA+1)
               ENDDO
            ENDDO
            CALL MPI_SEND(F(1),MN,MPI_DOUBLE_COMPLEX,
     &           MYRANK+1,1,MPI_COMM_WORLD,IERR)
C     
            ISTAR=IEND+1
            IF(MYRANK.EQ.NPROCS-2.AND.IEND-ISTA+1.LT.M) THEN
               ISUB=NMAX-ISTAR+1
            ENDIF
            MN=ISUB*(2*M+1)
            CALL MPI_RECV(F(1),MN,MPI_DOUBLE_COMPLEX,
     &           MYRANK+1,1,MPI_COMM_WORLD,ISTATUS,IERR)
            MN=0
            DO I=ISTAR,ISTAR+ISUB-1
               DO J=1,2*M+1
                  MN=MN+1
                  SR(J,I-ISTA+1)=F(MN)
               ENDDO
            ENDDO
C     
         ENDIF
C
         IF(MYRANK.EQ.NPROCS-1.AND.IEND-ISTA+1.LT.M) THEN
            ISUB=IEND-ISTA+1
         ENDIF
         MN=0
         DO I=ISTA,ISTA+ISUB-1
            DO J=1,2*M+1
               MN=MN+1
               F(MN)=SR(J,I-ISTA+1)
            ENDDO
         ENDDO
         CALL MPI_SEND(F(1),MN,MPI_DOUBLE_COMPLEX,
     &        MYRANK-1,1,MPI_COMM_WORLD,IERR)
C
      ENDIF
C     
  100 CONTINUE
C
      RETURN
      END
C
