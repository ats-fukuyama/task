! librk.f90

MODULE librk

  PRIVATE
  PUBLIC oderk
  PUBLIC oderkn

CONTAINS

! *** Ordinary differential equation ***

! --- old version ---

      SUBROUTINE ODERK(NEQ,SUB,X0,XE,N,Y0,YN,WORK)

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER, INTENT(IN)   :: NEQ, N
      REAL(dp)   , INTENT(INOUT):: X0
      REAL(dp)   , INTENT(IN)   :: XE
      REAL(dp), INTENT(INOUT), DIMENSION(NEQ) :: Y0, YN
      REAL(dp), INTENT(INOUT), DIMENSION(NEQ,2) :: WORK
      REAL(dp)     :: H
      INTEGER  :: I,J

      EXTERNAL SUB

      H = (XE - X0) / N
      DO I = 1,N
         CALL RKSTEP(NEQ,SUB,X0,H,Y0,YN,WORK(1,1),WORK(1,2))
         X0 = X0 + H
         DO J = 1,NEQ
            Y0(J) = YN(J)
         ENDDO
      ENDDO
      X0 = XE
      RETURN
      END SUBROUTINE ODERK

! --- old version ---

      SUBROUTINE ODERKN(neq_max,sub,xs,xe,nt_max,ys,ye)

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER, INTENT(IN)   :: neq_max,nt_max
      REAL(dp)   , INTENT(IN)   :: xs
      REAL(dp)   , INTENT(IN)   :: xe
      REAL(dp), INTENT(IN), DIMENSION(neq_max) :: ys
      REAL(dp), INTENT(OUT), DIMENSION(neq_max) :: ye
      REAL(dp), ALLOCATABLE:: y(:),work(:,:)
      REAL(dp)     :: xstep,x
      INTEGER  :: neq,nt

      EXTERNAL sub

      ALLOCATE(y(neq_max),work(neq_max,2))

      x = xs
      xstep = (xe - xs) / nt_max
      DO neq = 1,neq_max
         y(neq)=ys(neq)
      END DO
      DO nt= 1,nt_max
         CALL RKSTEP(neq_max,sub,x,xstep,y,ye, &
              work(1:neq_max,1),work(1:neq_max,2))
         x = x + xstep
         DO neq = 1,neq_max
            y(neq) = ye(neq)
         ENDDO
      ENDDO

      DEALLOCATE(y,work)
      RETURN
      END SUBROUTINE ODERKN

! --- 4th-order Runge-Kutta one stap ---

      SUBROUTINE RKSTEP(NEQ,SUB,X,H,Y0,YN,AK,W)

      USE task_kinds,ONLY: dp
      IMPLICIT NONE
      INTEGER, INTENT(IN)   :: NEQ
      REAL(dp)   , INTENT(INOUT):: X, H
      REAL(dp), INTENT(INOUT), DIMENSION(NEQ) :: Y0, YN, AK, W
      INTEGER  :: I
      REAL(dp), PARAMETER :: A2 = 0.5D0, A3 = A2, B2 = 0.5D0, B3 = B2
      REAL(dp), PARAMETER :: C1 = 1.0D0/6.D0, C2 = 1.0D0/3.D0, C3 = C2, C4 = C1
      EXTERNAL SUB

      CALL SUB(X,Y0,AK)
      DO I = 1,NEQ
         YN(I) = Y0(I) + H * C1 * AK(I)
      ENDDO
      DO I = 1,NEQ
         W(I) = Y0(I) + H * B2 * AK(I)
      ENDDO

      CALL SUB(X + A2 * H,W,AK)
      DO I = 1,NEQ
         YN(I) = YN(I) + H * C2 * AK(I)
      ENDDO
      DO I = 1,NEQ
         W(I) = Y0(I) + H * B3 * AK(I)
      ENDDO

      CALL SUB(X + A3 * H,W,AK)
      DO I = 1,NEQ
         YN(I) = YN(I) + H * C3 * AK(I)
      ENDDO
      DO I = 1,NEQ
         W(I) = Y0(I) + H * AK(I)
      ENDDO

      CALL SUB(X + H,W,AK)
      DO I = 1,NEQ
         YN(I) = YN(I) + H * C4 * AK(I)
      ENDDO
      RETURN
      END SUBROUTINE RKSTEP
END MODULE librk
