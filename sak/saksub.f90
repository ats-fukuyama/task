! saksub.f90

MODULE saksub
  USE sakcomm
  REAL(dp):: sg2_local,rk2_local

  PRIVATE
  PUBLIC cfeps,subeps,sg2_local,rk2_local

CONTAINS

  FUNCTION cfeps(cw,rk2,sg2)
    IMPLICIT NONE
    COMPLEX(dp),INTENT(IN):: cw
    REAL(dp),INTENT(IN):: rk2,sg2
    COMPLEX(dp):: cfeps

    cfeps= 1.D0-(1.D0+sg2+3.D0*(1+sg2)**2*rk2/cw**2)/cw**2 &
         +CI*SQRT(0.5D0*Pi)*cw/(RK2**1.5D0*SQRT(1.D0+sg2)) &
         *EXP(-0.5D0*cw**2/(rk2*(1.D0+sg2)))
  END FUNCTION cfeps

  SUBROUTINE subeps(x,y,f)
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: x,y
    REAL(dp),INTENT(OUT):: f
    COMPLEX(dp):: cw,cf

    cw=CMPLX(x,y)
    cf=cfeps(cw,rk2_local,sg2_local)
    f=REAL(cf)**2+AIMAG(cf)**2
    RETURN
  END SUBROUTINE subeps
    
END MODULE saksub
