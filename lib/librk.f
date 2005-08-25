C     $Id$
C
      SUBROUTINE RK(NEQ,SUB,X0,XE,N,Y0,YN,WORK)
C
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL SUB
      DIMENSION Y0(NEQ),YN(NEQ),WORK(NEQ,2)
C
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
      END
C
      SUBROUTINE RKSTEP(NEQ,SUB,X,H,Y0,YN,AK,W)
C
      IMPLICIT REAL * 8(A-H,O-Z)
      PARAMETER(A2 = 0.5D0, A3 = A2)
      PARAMETER(B2 = 0.5D0, B3 = B2)
      PARAMETER(C1 = 1.0D0/6.D0, C4 = C1)
      PARAMETER(C2 = 1.0D0/3.D0, C3 = C2)
      DIMENSION Y0(NEQ),YN(NEQ),AK(NEQ),W(NEQ)
      EXTERNAL SUB
C
      CALL SUB(X,Y0,AK)
      DO I = 1,NEQ
         YN(I) = Y0(I) + H * C1 * AK(I)
      ENDDO
      DO I = 1,NEQ
         W(I) = Y0(I) + H * B2 * AK(I)
      ENDDO
C
      CALL SUB(X + A2 * H,W,AK)
      DO I = 1,NEQ
         YN(I) = YN(I) + H * C2 * AK(I)
      ENDDO
      DO I = 1,NEQ
         W(I) = Y0(I) + H * B3 * AK(I)
      ENDDO
C
      CALL SUB(X + A3 * H,W,AK)
      DO I = 1,NEQ
         YN(I) = YN(I) + H * C3 * AK(I)
      ENDDO
      DO I = 1,NEQ
         W(I) = Y0(I) + H * AK(I)
      ENDDO
C
      CALL SUB(X + H,W,AK)
      DO I = 1,NEQ
         YN(I) = YN(I) + H * C4 * AK(I)
      ENDDO
      RETURN
      END
