C     $Id$
C
C     ****** SOLUTION OF BAND MATRIX (GAUSSIAN ELIMINATION IN MPI) ******
C
      SUBROUTINE BANDCDB(A,D,X,N,L,LA,LD,NM,TMPA,TMPG,TMPB,
     &                   TINV1,TEMP1,TEMP2,NBSIZ,IERR)
C
      COMPLEX * 16 A (LA,N),D(N),X(N)
      COMPLEX * 16 TMPA (LD,NBSIZ)
      COMPLEX * 16 TMPG (LD,NBSIZ)
      COMPLEX * 16 TMPB (NBSIZ)
      COMPLEX * 16 TINV1(LD,NBSIZ)
      COMPLEX * 16 TEMP1(LD,NBSIZ)
      COMPLEX * 16 TEMP2(NBSIZ)
      INTEGER      NM(NBSIZ)
C
      NRP=0
      DO J=1,N
         D(J)=(0.D0,0.D0)
         DO I=1,L
            A(I,J)=(0.D0,0.D0)
         ENDDO
         CALL WMSETM(A(1,J),D(J),J,LA,NRP)
      ENDDO
C
      IAF=2*NBSIZ-NBSIZ
      IBF=2*NBSIZ
      ICF=2*NBSIZ+NBSIZ
C
C     ###### UPWARD SCANNING ######
C
C     ###### SET INITIAL VALUE ######
C
      IBLC=NBSIZ*(N/NBSIZ-1)+1
C
      DO J=1,NBSIZ
         DO I=1,NBSIZ
            TINV1(I,J)=A(IBF+(I-1)-(J-1),IBLC+(J-1))
         ENDDO
      ENDDO
C
      CALL INVMCD(TINV1,NBSIZ,LD,NM,IERR)
      if(ierr.ne.0) write(6,*) '1 ierr=',ierr
C
      DO J=1,NBSIZ
C
         TMPB(J)=(0.D0,0.D0)
         DO K=1,NBSIZ
            TMPB(J)=TMPB(J)
     &             +TINV1(K,J)*D(IBLC+(K-1))
         ENDDO
C
         DO I=1,NBSIZ
C
            TMPA(I,J)=(0.D0,0.D0)
            DO K=1,NBSIZ
               TMPA(I,J)=TMPA(I,J)
     &                  -TINV1(K,J)
     &                  *A(IAF+(I-1)-(K-1),IBLC+(K-1))
            ENDDO
C
            TMPG(I,J)=(0.D0,0.D0)
C
         ENDDO
C
      ENDDO
C
      DO J=1,NBSIZ
         DO I=1,NBSIZ
            A(IAF+(I-1)-(J-1),IBLC+(J-1))=TMPA(I,J)
            A(IBF+(I-1)-(J-1),IBLC+(J-1))=TMPG(I,J)
         ENDDO
         D(IBLC+(J-1))=TMPB(J)
      ENDDO
C
C     ###### UPWARD SCANNING MAIN ######
C
      DO IBN=N/NBSIZ-1,2,-1
         IBLC=NBSIZ*(IBN-1)+1
C
         CALL CALNCFN(A,LA,N,IBLC,ICF,IBF,TMPA,TINV1,
     &                LD,NBSIZ,IERR)
C
C     ###### CALCULATE BETA ######
C
         DO J=1,NBSIZ
            TEMP2(J)=D(IBLC+(J-1))
            DO K=1,NBSIZ
               TEMP2(J)=TEMP2(J)
     &                 -A(ICF+(K-1)-(J-1),IBLC+(J-1))
     &                 *TMPB(K)
            ENDDO
         ENDDO
C
         DO J=1,NBSIZ
            TMPB(J)=(0.D0,0.D0)
            DO K=1,NBSIZ
               TMPB(J)=TMPB(J)+TINV1(K,J)*TEMP2(K)
            ENDDO
         ENDDO
C
C     ###### CALCULATE ALPHA ######
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               TMPA(I,J)=(0.D0,0.D0)
               DO K=1,NBSIZ
                  TMPA(I,J)=TMPA(I,J)
     &                     -TINV1(K,J)
     &                     *A(IAF+(I-1)-(K-1),IBLC+(K-1))
               ENDDO
            ENDDO
         ENDDO
C
C     ###### CALCULATE GAMMA ######
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               TEMP1(I,J)=(0.D0,0.D0)
               DO K=1,NBSIZ
                  TEMP1(I,J)=TEMP1(I,J)
     &                      -A(ICF+(K-1)-(J-1),IBLC+(J-1))
     &                      *TMPG(I,K)
               ENDDO
            ENDDO
         ENDDO
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               TMPG(I,J)=(0.D0,0.D0)
               DO K=1,NBSIZ
                  TMPG(I,J)=TMPG(I,J)+TINV1(K,J)*TEMP1(I,K)
               ENDDO
            ENDDO
         ENDDO
C
         DO J=1,NBSIZ
            DO I=1,NBSIZ
               A(IAF+(I-1)-(J-1),IBLC+(J-1))=TMPA(I,J)
               A(IBF+(I-1)-(J-1),IBLC+(J-1))=TMPG(I,J)
            ENDDO
            D(IBLC+(J-1))=TMPB(J)
         ENDDO
C
      ENDDO
C
C     ####### SOLUTION VECTOR ######
C
      DO J=1,NBSIZ
         DO I=1,NBSIZ
            TINV1(I,J)=A(IBF+(I-1)-(J-1),J)
            DO K=1,NBSIZ
               TINV1(I,J)=TINV1(I,J)+A(ICF+(K-1)-(J-1),J)
     &                              *TMPA(I,K)
            ENDDO
         ENDDO
      ENDDO
C
      CALL INVMCD(TINV1,NBSIZ,LD,NM,IERR)
C
      DO J=1,NBSIZ
         TEMP2(J)=D(J)
         DO K=1,NBSIZ
            TEMP2(J)=TEMP2(J)-A(ICF+(K-1)-(J-1),J)*TMPB(K)
         ENDDO
      ENDDO
C
      DO J=1,NBSIZ
         X(J)=(0.D0,0.D0)
         DO K=1,NBSIZ
            X(J)=X(J)+TINV1(K,J)*TEMP2(K)
         ENDDO
      ENDDO
C
      DO IBN=2,N/NBSIZ
         IBLC=NBSIZ*(IBN-1)+1
C
         DO I=1,NBSIZ
            X(IBLC+(I-1))=D(IBLC+(I-1))
            DO K=1,NBSIZ
               X(IBLC+(I-1))=X(IBLC+(I-1))
     &                      +A(IAF+(K-1)-(I-1),IBLC+(I-1))
     &                      *X(IBLC-NBSIZ+(K-1))
            ENDDO
         ENDDO
C
      ENDDO
C
      DO I=1,N
         D(I)=X(I)
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE NEXT BLOCK COEFFICIENT ******
C
      SUBROUTINE CALNCFN(A,LA,N,IBLC,IAF,IBF,TEMP1,TINV1,
     &                   LD,NBSIZ,IERR)
C
      COMPLEX * 16 A(LA,N)
      COMPLEX * 16 TEMP1(LD,NBSIZ)
      COMPLEX * 16 TINV1(LD,NBSIZ)
C
      DO J=1,NBSIZ
         DO I=1,NBSIZ
            TINV1(I,J)=A(IBF+(I-1)-(J-1),IBLC+(J-1))
            DO K=1,NBSIZ
               TINV1(I,J)=TINV1(I,J)
     &                   +A(IAF+(K-1)-(J-1),IBLC+(J-1))
     &                   *TEMP1(I,K)
            ENDDO
         ENDDO
      ENDDO
C
      CALL INVMCD(TINV1,NBSIZ,LD,NM,IERR)
      if(ierr.ne.0) write(6,*) '3 ierr=',ierr
C
      RETURN
      END














