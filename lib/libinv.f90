!     $Id$

!     ****** REAL MATRIX INVERSION ******

      SUBROUTINE INVMRD(A,N,NA,ILL)

      implicit none
      integer(4),                intent(in)    :: N, NA
      real(8), dimension(NA,NA), intent(inout) :: A
      integer(4),                intent(out)   :: ILL
      integer(4) :: NN, NW, IP, I, J
      integer(4), dimension(N) :: NM
      real(8) :: W, P, EPS
      logical :: LEX
      DATA EPS/1.D-14/

!!$      IF(N > 128) THEN
!!$         WRITE(6,*) 'XX INVMRD: N > 128: N=',N
!!$         ILL=999
!!$         RETURN
!!$      END IF

      DO NN=1,N
         NM(NN) = NN
      END DO

      DO NN=1,N
         P=0.D0
         IP=NN
         DO I=NN,N
            IF(ABS(A(I,1)) > P) THEN
               P=ABS(A(I,1))
               IP=I
            ENDIF
         END DO
         IF(P <= EPS) THEN
            ILL=900
            RETURN
         END IF
         IF(IP /= NN) THEN
            NW=NM(IP)
            NM(IP)=NM(NN)
            NM(NN)=NW
            DO J=1,N
               W=A(IP,J)
               A(IP,J)=A(NN,J)
               A(NN,J)=W
            END DO
         ENDIF
         W=1.D0/A(NN,1)
         A(NN,1:N-1)=A(NN,2:N)*W
         A(NN,N)=W
         DO I=1,N
            IF(I /= NN) THEN
               W=A(I,1)
               A(I,1:N-1)=A(I,2:N)-W*A(NN,1:N-1)
               A(I,N)=-W*A(NN,N)
            ENDIF
         END DO
      END DO

      OUT:DO NN=1,N
         LEX = .FALSE.
         IN:DO J=NN,N
            IF(NM(J) == NN) THEN
               LEX = .TRUE.
               EXIT IN
            END IF
         END DO IN
         IF(LEX.EQV..FALSE.) THEN
            ILL=900
            RETURN
         END IF

         IF(J /= NN) THEN
            NM(J)=NM(NN)
            DO I=1,N
               W=A(I,J)
               A(I,J)=A(I,NN)
               A(I,NN)=W
            END DO
         ENDIF
      END DO OUT
      ILL=0
      RETURN

      END SUBROUTINE INVMRD

!     ***** COMPLEX MATRIX INVERSION *****

      SUBROUTINE INVMCD(A,N,NA,ILL)

      implicit none
      integer(4),                   intent(in)    :: N, NA
      complex(8), dimension(NA,NA), intent(inout) :: A
      integer(4),                   intent(out)   :: ILL
      integer(4) :: NN, NW, IP, I, J
      integer(4), dimension(N) :: NM
      real(8)    :: P, EPS
      complex(8) :: W
      logical    :: LEX
      DATA EPS/1.D-14/

      DO NN=1,N
         NM(NN)=NN
      END DO

      DO NN=1,N
         P=0.D0
         IP=NN
         DO I=NN,N
            IF(ABS(A(I,1)) > P) THEN
               P=ABS(A(I,1))
               IP=I
            ENDIF
         END DO
         IF(P <= EPS) THEN
            ILL=NN
            RETURN
         END IF
         IF(IP /= NN) THEN
            NW=NM(IP)
            NM(IP)=NM(NN)
            NM(NN)=NW
            DO J=1,N
               W=A(IP,J)
               A(IP,J)=A(NN,J)
               A(NN,J)=W
            END DO
         ENDIF
         W=(1.D0,0.D0)/A(NN,1)
         A(NN,1:N-1)=A(NN,2:N)*W
         A(NN,N)=W
         DO I=1,N
            IF(I /= NN) THEN
               W=A(I,1)
               A(I,1:N-1)=A(I,2:N)-W*A(NN,1:N-1)
               A(I,N)=-W*A(NN,N)
            ENDIF
         END DO
      END DO

      OUT:DO NN=1,N
         IN:DO J=NN,N
            IF(NM(J) == NN) THEN
               LEX = .TRUE.
               EXIT IN
            END IF
         END DO IN
         IF(LEX.EQV..FALSE.) THEN
            ILL=900
            RETURN
         END IF

         IF(J /= NN) THEN
            NM(J)=NM(NN)
            DO I=1,N
               W=A(I,J)
               A(I,J)=A(I,NN)
               A(I,NN)=W
            END DO
         ENDIF
      END DO OUT
      ILL=0
      RETURN

      END SUBROUTINE INVMCD
