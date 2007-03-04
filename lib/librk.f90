      SUBROUTINE ODERK(NEQ,SUB,X0,XE,N,Y0,YN,WORK)

      IMPLICIT NONE
      INTEGER(4), INTENT(IN)   :: NEQ, N
      REAL(8)   , INTENT(INOUT):: X0
      REAL(8)   , INTENT(IN)   :: XE
      REAL(8), INTENT(INOUT), DIMENSION(NEQ) :: Y0, YN
      REAL(8), INTENT(INOUT), DIMENSION(NEQ,2) :: WORK
      REAL(8)     :: H
      INTEGER(4)  :: I,J

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

      SUBROUTINE RKSTEP(NEQ,SUB,X,H,Y0,YN,AK,W)

      IMPLICIT NONE
      INTEGER(4), INTENT(IN)   :: NEQ
      REAL(8)   , INTENT(INOUT):: X, H
      REAL(8), INTENT(INOUT), DIMENSION(NEQ) :: Y0, YN, AK, W
      INTEGER(4)  :: I
      REAL(8), PARAMETER :: A2 = 0.5D0, A3 = A2, B2 = 0.5D0, B3 = B2
      REAL(8), PARAMETER :: C1 = 1.0D0/6.D0, C2 = 1.0D0/3.D0, C3 = C2, C4 = C1
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
