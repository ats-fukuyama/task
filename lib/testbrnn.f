C
      REAL*8 X1,X2,FUNC
      EXTERNAL FUNC
C
      X1=0.D0
      X2=2.D0
      TOL=1.D-8
C
      F1=FBRENT(X1,X2,FUNC,TOL)
C
      F2=zeroin(X1,X2,FUNC,TOL)
C
      WRITE(6,*) F1,F2
C
      STOP
      END
C
      DOUBLE PRECISION FUNCTION FUNC(X)
C
      REAL*8 X
      FUNC=-0.5D0+SIN(X)
      RETURN
      END
