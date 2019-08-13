MODULE w1sub

CONTAINS

!     ****** INTERFACE FOR FFT ******

  SUBROUTINE W1FFTL(CAF,N,KEY)
    USE libfft
    USE w1comm,ONLY: rkind
    IMPLICIT NONE

    COMPLEX(rkind),INTENT(INOUT):: CAF(N)
    COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CBF
    INTEGER,INTENT(IN):: N,KEY

    ALLOCATE(CBF(N))
    CALL CFFT1D(CAF,CBF,N,KEY)
    CAF(1:N)=CBF(1:N)
    DEALLOCATE(CBF)

    RETURN
  END SUBROUTINE W1FFTL

!     H(X) FOR SLOWING-DOWN DISTRIBUTION

  FUNCTION W1_HBEAM(X)
    USE w1comm,ONLY: rkind,PI
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X
    REAL(rkind):: w1_HBEAM,SQR3

    SQR3=SQRT(3.D0)
    W1_HBEAM=(LOG((1.D0+X**3)/(1.D0+X)**3)/3.D0 &
           +(ATAN((2.D0*X-1.D0)/SQR3)+PI/6.D0)/SQR3)/X**2
    RETURN
  END FUNCTION W1_HBEAM

END MODULE w1sub
