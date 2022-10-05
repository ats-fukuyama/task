!     libplog.f90

MODULE libplog

CONTAINS

!     ***********************************************************

!           CEILING FUNCTION FOR LOG10 PLOT

!     ***********************************************************

  FUNCTION PLOG(X,XMIN,XMAX)

    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    REAL(dp), INTENT(IN) ::  X, XMIN, XMAX
    REAL(dp):: PLOG

    IF(X.LT.XMIN) THEN
       PLOG=LOG10(XMIN)
    ELSEIF(X.GT.XMAX) THEN
       PLOG=LOG10(XMAX)
    ELSE
       PLOG=LOG10(X)
    ENDIF

    RETURN
  END FUNCTION PLOG
END MODULE libplog
