! wmlib.f90

MODULE wmlib

  PRIVATE
  PUBLIC adjust_power2

CONTAINS

  ! *** adjust an integer n to the nearest power of two ***
  
  SUBROUTINE adjust_power2(n,n_adjusted)
    USE wmcomm,ONLY: rkind
    IMPLICIT NONE
    INTEGER,INTENT(IN):: n
    INTEGER,INTENT(OUT):: n_adjusted
    REAL(rkind):: power2
    INTEGER:: ipower2

    power2=LOG(DBLE(n))/LOG(2.D0)
    ipower2=NINT(power2)
    IF(ABS(power2-DBLE(ipower2)).GE.1.D-6) THEN
       n_adjusted=2**ipower2
       IF(n_adjusted.LT.n) n_adjusted=2**(ipower2+1)
    ELSE
       n_adjusted=n
    END IF
    RETURN
  END SUBROUTINE adjust_power2
END MODULE wmlib
