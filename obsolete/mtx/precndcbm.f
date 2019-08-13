C     $Id$
C
C     ######## CALCULATE D INVERCE ########
C
      SUBROUTINE CALDINV(A,D,NM,LA,LD,N,M,NBSIZ,NFST,NEND,IERR)
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
      COMMON /MPIVR3/ NBST,NBED
C
      COMMON /WMBLPR/ ICF,IDF,IEF,ISCNF,IECNF
C
      COMPLEX * 16 A(LA,NFST:NEND),D(LD,NFST:NEND)
      COMPLEX * 16 NM(NBSIZ)
C
      DO IBN=NBST,NBED
         IBLC=NBSIZ*(IBN-1)+1
C
         DO I=1,NBSIZ
            DO J=1,NBSIZ
               D(J,IBLC+(I-1)-ISTA+1)
     &              =A(IDF+(J-1)-(I-1),IBLC+(I-1)-ISTA+1)
            ENDDO
         ENDDO
C
         CALL INVMCD(D(1,IBLC-ISTA+1),NBSIZ,LD,IERR)
C
      ENDDO
C
      RETURN
      END
C
C
C     ######## PRECONDITIONING ########
C                  (JACOBI)
C
      SUBROUTINE CALJACBM(A,D,TEMP2,X,LA,LD,N,M,NBSIZ,
     &                    NFST,NEND,IERR)
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
      COMMON /MPIVR3/ NBST,NBED
C
      COMMON /WMBLPR/ ICF,IDF,IEF,ISCNF,IECNF
C
      COMPLEX * 16 A(LA,NFST:NEND),D(LD,NFST:NEND)
      COMPLEX * 16 X(NFST:NEND),TEMP2(NBSIZ)
C
      IF(N.LE.0.OR.NBSIZ.LE.0) THEN
         IERR=1
      ENDIF
C
      DO IBN=NBST,NBED
         IBLC=NBSIZ*(IBN-1)+1
C
         DO I=1,NBSIZ
            TEMP2(I)=(0.D0,0.D0)
            DO K=1,NBSIZ
               TEMP2(I)=TEMP2(I)
     &              +D(K,IBLC+(I-1)-ISTA+1)*X(IBLC+(K-1)-ISTA+1)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
           X(IBLC+(K-1)-ISTA+1)=TEMP2(I)
        ENDDO
C
      ENDDO
C
      RETURN
      END
C
C     ######## PRECONDITIONING ROUTINE ########
C       (INCOMPLETE LDU DICOMPOSITION VERSION)
C     
      SUBROUTINE CALCUDBM(A,D,TEMP1,NM,LA,LD,N,M,NBSIZ,
     &                    NFST,NEND,IERR)
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
      COMMON /MPIVR3/ NBST,NBED
C
      COMMON /WMBLPR/ ICF,IDF,IEF,ISCNF,IECNF
C
      COMPLEX * 16 A(LA,NFST:NEND),D(LD,NFST:NEND)
      COMPLEX * 16 TEMP1(LD,NBSIZ),NM(NBSIZ)
C
      IF(N.LE.0.OR.NBSIZ.LE.0) THEN
         IERR=1
      ENDIF
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
         CALL INVMCD(D(1,IBLC-ISTA+1),NBSIZ,LD,IERR)
C
      ENDDO
C
      RETURN
      END
C
C
C     ######## CALCULATE x=(LDU)^{-1}*y ########
C
      SUBROUTINE CALLDUBM(A,D,TEMP2,X,LA,LD,N,M,NBSIZ,
     &                    NFST,NEND,IERR)
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
      COMMON /MPIVR3/ NBST,NBED
C
      COMMON /WMBLPR/ ICF,IDF,IEF,ISCNF,IECNF
C
      COMPLEX * 16 A(LA,NFST:NEND),D(LD,NFST:NEND)
      COMPLEX * 16 X(NFST:NEND),TEMP2(NBSIZ)
C
      IF(N.LE.0.OR.NBSIZ.LE.0) THEN
         IERR=1
      ENDIF
C
      DO IBN=NBST,NBED
         IBLC=NBSIZ*(IBN-1)+1
         DO I=1,NBSIZ
            TEMP2(I)=X(IBLC+(I-1)-ISTA+1)
            DO IC=MAX(IBLC-NBSIZ,ISCNF),IBLC-1
               K=IC+NBSIZ-IBLC+1
               TEMP2(I)=TEMP2(I)
     &                 -A(ICF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                 *X(IC-ISTA+1)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            X(IBLC+(I-1)-ISTA+1)=(0.D0,0.D0)
            DO K=1,NBSIZ
               X(IBLC+(I-1)-ISTA+1)=X(IBLC+(I-1)-ISTA+1)
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
     &                 *X(IC-ISTA+1)
            ENDDO
         ENDDO
         DO I=1,NBSIZ
            DO K=1,NBSIZ
               X(IBLC+(I-1)-ISTA+1)=X(IBLC+(I-1)-ISTA+1)
     &                       -D(K,IBLC+(I-1)-ISTA+1)
     &                       *TEMP2(K)
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C
C     ######## LDU DECOMPOSITION ########
C            (KEEP INDEPENDANCE)
C
      SUBROUTINE CALCUDBMS(A,D,TEMP1,NM,LA,LD,N,M,NBSIZ,F,
     &                    NFST,NEND,IERR)
C
      INCLUDE '../mpi/mpif.inc'
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
      COMMON /MPIVR3/ NBST,NBED
C
      COMMON /WMBLPR/ ICF,IDF,IEF,ISCNF,IECNF
C
      COMPLEX * 16 A(LA,NFST:NEND),D(LD,NFST:NEND)
      COMPLEX * 16 TEMP1(LD,NBSIZ),NM(NBSIZ)
      COMPLEX * 16 F(M*(2*M+1))
      INTEGER      ISTATUS(MPI_STATUS_SIZE)
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
               CALL INVMCD(D(1,IBLC-ISTA+1),NBSIZ,LD,IERR)
C
            ENDDO
C
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      ENDDO
C
      RETURN
      END
C
C     ######## y=(LDU)^{-1}x ##########
C          (KEEP INDEPENDENCE)
C
      SUBROUTINE CALLDUBMS(A,D,TEMP2,X,LA,LD,N,M,NBSIZ,F,
     &                    NFST,NEND,IERR)
C
      INCLUDE '../mpi/mpif.inc'
C
      COMMON /MPIVR1/ NPROCS,MYRANK
      COMMON /MPIVR2/ ISTA,IEND
      COMMON /MPIVR3/ NBST,NBED
C
      COMMON /WMBLPR/ ICF,IDF,IEF,ISCNF,IECNF
C
      COMPLEX * 16 A(LA,NFST:NEND),D(LD,NFST:NEND)
      COMPLEX * 16 X(NFST:NEND),TEMP2(NBSIZ)
      COMPLEX * 16 F(M*(2*M+1))
      INTEGER      ISTATUS(MPI_STATUS_SIZE)
C
      DO NRANK=0,NPROCS-1
C
         IF(MYRANK.EQ.NRANK-1) THEN
            CALL MPI_SEND(X(IEND-NBSIZ+1-ISTA+1),NBSIZ,
     &           MPI_DOUBLE_COMPLEX,MYRANK+1,1,MPI_COMM_WORLD,IERR)
         ELSE IF(MYRANK.EQ.NRANK) THEN
            IF(MYRANK.NE.0) THEN
               CALL MPI_RECV(X(1-NBSIZ),NBSIZ,MPI_DOUBLE_COMPLEX,
     &              MYRANK-1,1,MPI_COMM_WORLD,ISTATUS,IERR)
            ENDIF
C
            DO IBN=NBST,NBED
               IBLC=NBSIZ*(IBN-1)+1
               DO I=1,NBSIZ 
                  TEMP2(I)=X(IBLC+(I-1)-ISTA+1)
                  DO IC=MAX(IBLC-NBSIZ,1),IBLC-1
                     K=IC+NBSIZ-IBLC+1
                     TEMP2(I)=TEMP2(I)
     &                    -A(ICF+(K-1)-(I-1),IBLC+(I-1)-ISTA+1)
     &                    *X(IC-ISTA+1)
                  ENDDO
               ENDDO
               DO I=1,NBSIZ
                  X(IBLC+(I-1)-ISTA+1)=(0.D0,0.D0)
                  DO K=1,NBSIZ
                     X(IBLC+(I-1)-ISTA+1)=X(IBLC+(I-1)-ISTA+1)
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
            CALL MPI_SEND(X(1),NBSIZ,MPI_DOUBLE_COMPLEX,
     &           MYRANK-1,1,MPI_COMM_WORLD,IERR)
         ELSE IF(MYRANK.EQ.NRANK) THEN
            IF(MYRANK.NE.NPROCS-1) THEN
               CALL MPI_RECV(X(IEND+1-ISTA+1),NBSIZ,
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
     &                    *X(IC-ISTA+1)
                  ENDDO
               ENDDO
               DO I=1,NBSIZ
                  DO K=1,NBSIZ
                     X(IBLC+(I-1)-ISTA+1)=X(IBLC+(I-1)-ISTA+1)
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
      RETURN
      END
