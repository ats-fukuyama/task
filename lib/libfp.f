C     $Id$
C
C
C ***********************************************************
C
C     ERROR FUNCTION
C
C ***********************************************************
C
      FUNCTION ERF0(X)
      REAL*8 ERF0,X,DERF
C
      ERF0=DERF(X)
      RETURN
      END
C
C ***********************************************************
C
C     DERIVATIVE OF ERROR FUNCTION
C
C ***********************************************************
C
      FUNCTION ERF1(U)
C
      IMPLICIT REAL*8(A-H,O-R,T-Z)
C
      PI=3.14159265358979323846D0
      U2=U**2
      IF (U2.GT.100.D0)THEN
         ERF1=0.D0
      ELSE
         ERF1=2/SQRT(PI)*EXP(-U2)
      END IF
C
      RETURN
      END
C
C ***************************************************************
C
C     LEGENDRE POLYNOMIALS
C
C ***************************************************************
C
      SUBROUTINE DPLEG(X,N,Y,IERR)
C
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION Y(*)
C
      IERR = 0
      IF(N .LT. 0) THEN
         WRITE(6,*) 'XX DPLEG: N.LT.0: N=',N
         IERR = 1
         RETURN
      ENDIF
C
      Y(1) = 1.D0
      IF(N .EQ. 0) RETURN
      Y(2) = X
      IF(N .EQ. 1) RETURN
      DO I=2,N
         W  = X * Y(I)
         WY = W - Y(I-1)
         Y(I+1) = WY + W - WY / DBLE(I)
      ENDDO
      RETURN
      END
