C
      SUBROUTINE INVEMATRIX(NDIM,NDIMM,MATRX,INMATRX,WORK)
C
      IMPLICIT NONE
      INTEGER I,J,K,NDIM,NDIMM
      REAL*8 MATRX(NDIMM,*),INMATRX(NDIMM,*)
      REAL*8 WORK(NDIMM,*)
C
C      WRITE(*,*) 'INPUT DIMENSION OF MATRIX?'
C      READ(*,*) NDIM
C
      DO I=1,NDIM
         DO J=1,NDIM
C
C            WRITE(*,*) 'INPUT NUMBER OF MATRIX(',I,',',J,')?'
C            READ(*,*) WORK(I,J)
C
            WORK(I,J)=MATRX(I,J)
            IF (I .EQ. J) THEN
	       WORK(I,J+NDIM)=1
            ELSE
	       WORK(I,J+NDIM)=0
            END IF
         END DO
      END DO
C 
      DO K=1,NDIM
         DO I=1,NDIM
            DO J=1,NDIM*2
               IF (K .NE. I) THEN
                  IF (K .NE. J) THEN
                     WORK(I,J)=WORK(I,J)
     &                        -WORK(K,J)*(WORK(I,K)/WORK(K,K))
                  END IF
               END IF
            END DO
         END DO
C
         DO J=1,NDIM*2
            IF (K .NE. J) THEN
               WORK(K,J)=WORK(K,J)/WORK(K,K)
            END IF
         END DO
         DO I=1,NDIM
            WORK(I,K)=0
         END DO
         WORK(K,K)=1
      END DO
C
      DO I=1,NDIM
      DO J=1,NDIM
         INMATRX(I,J)=WORK(I,NDIM+J)
      END DO
      END DO
C
 2000 FORMAT(1X,100F7.2)
C      DO I=1,NDIM
C            WRITE(*,2000) (WORK(I,J),J=NDIM+1,DIM*2)
C      END DO
      RETURN
      END
