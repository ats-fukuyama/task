C     $Id$
C
C **************************************
C     REAL*8 UNDERFLOW ERROR ROUTINE
C **************************************
C
      SUBROUTINE ERUFL8(A)
      REAL*8 A
      A=0.D0
      RETURN
      END
C
C ***********************************************************
C
C                    ERROR FUNCTION(1)
C
C ***********************************************************
C
      FUNCTION ERF1(U)
C
      IMPLICIT REAL*8(A-H,O-R,T-Z)
C$ALIAS DERF='erf'(%VAL)
C
      ERF1=ERF(U)
C
      RETURN
      END
C
C ***********************************************************
C
C                    ERROR FUNCTION(2)
C
C ***********************************************************
C
      FUNCTION ERF2(U)
C
      IMPLICIT REAL*8(A-H,O-R,T-Z)
C
      PI=3.14159265358979323846D0
      U2=U**2
      IF (U2.GT.100.D0)THEN
         ERF2=0.D0
      ELSE
         ERF2=2/SQRT(PI)*EXP(-U2)
      END IF
C
      RETURN
      END
C
C ***************************************************************
C
      SUBROUTINE DPLEG( X, N, Y, IER)
C
C ***************************************************************
C
C ***************************************************************
C *   SUBROUTINE DPLEG(X,N,Y,IER)                               *
C *   DPLEG CALCULATES LEGENDRE POLYNOMIALS P(N,X)  (DOUBLE)    *
C *      INPUT PARAMETERS:                                      *
C *        X = ARGUMENT OF P (REAL*8)                           *
C *        N = MAXIMUM DEGREE OF P (INTEGER)                    *
C *      OUTPUT DATA:                                           *
C *        Y(N) = VALUE OF P(N,X) (REAL*8)                      *
C *        IER = 0  (N.GE.0) :ERROR CONDITION CODE              *
C *            = -1 (N.LT.0)                                    *
C ***************************************************************
C
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION Y(*)
C
C     TEST OF (N)
      IF(N .LT. 0) THEN
         WRITE(6,1000) N
         IER = -1
         RETURN
      ENDIF
      IER = 0
      Y(1) = 1.D0
      IF(N .EQ. 0) RETURN
      Y(2) = X
      IF(N .EQ. 1) RETURN
      DO 20 I=2,N
         W  = X * Y(I)
         WY = W - Y(I-1)
   20 Y(I+1) = WY + W - WY / DBLE(I)
      RETURN
 1000 FORMAT(1H0,
     &       '(SUBR. DPLEG) N=',I5,', N MUST BE NON-NEGATIVE.')
      END
