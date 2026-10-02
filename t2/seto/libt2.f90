!    $Id$

MODULE LIBT2

  USE bpsd,ONLY: rkind

  PRIVATE
  PUBLIC integ_f,func_1,x_ab,func_phi,func_G

  REAL(rkind):: x_ab

CONTAINS

  FUNCTION func_phi(x)
    USE bpsd,ONLY: rkind
    USE libspf,ONLY: ERF0
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: x
    REAL(rkind):: func_phi
    func_phi=ERF0(x)
  END FUNCTION func_phi
    
  FUNCTION func_G(x)
    USE bpsd,ONLY: rkind
    USE libspf,ONLY: ERF0,ERF1
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: x
    REAL(rkind):: func_G
    IF(x.LE.0.D0) THEN
       func_G=0.D0
    ELSE
       func_G=(ERF0(x)-x*ERF1(x))/(2*x*x)
    END IF
  END FUNCTION func_G
    
  FUNCTION func_1(x2)
    USE bpsd,ONLY: rkind,pi
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: x2
    REAL(rkind):: func_1,x_a,x_b,f

    x_a=SQRT(x2) ! x_a=v/v_a   x_ab=v_b/v_a
    x_b=x_a/x_ab
    f=(func_phi(x_b)-func_G(x_b))/(x_a)**3
    func_1=(4.D0/(3.D0*SQRT(pi)))*f*EXP(-x_a**2)*x_a**3
    print*,'x_a=',x_a,'x_ab=',x_ab
    RETURN
  END FUNCTION func_1

  SUBROUTINE integ_f(FUNC,FINT,FERR,EPS,H0,KID,ILST)

    USE bpsd,ONLY: rkind
    USE LIBDE,ONLY: DEHIFE
    IMPLICIT NONE
    REAL(rkind),INTENT(OUT):: FINT   ! Integral
    REAL(rkind),INTENT(OUT):: FERR   ! Estimated error 
    REAL(rkind),INTENT(IN),OPTIONAL::  EPS  ! Convergence thrshold
    REAL(rkind),INTENT(IN),OPTIONAL::  H0   ! Initial step size
    INTEGER,INTENT(IN),OPTIONAL::  ILST ! print out control: 0 for no print out
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: KID   ! function identifier string
    INTERFACE
       FUNCTION FUNC(X2) ! X2 = x^2, use X=SQRT(X2)
         USE bpsd,ONLY: rkind
         REAL(rkind):: FUNC
         REAL(rkind),INTENT(IN):: X2
       END FUNCTION FUNC
    END INTERFACE
    REAL(rkind):: EPS_,H0_
    INTEGER:: ILST_

    IF(PRESENT(EPS)) THEN
       EPS_=EPS
    ELSE
       EPS_=1.D-14
    END IF
    IF(PRESENT(H0)) THEN
       H0_=H0
    ELSE
       H0_=0.5D0
    END IF
    IF(PRESENT(ILST)) THEN
       ILST_=ILST
    ELSE
       ILST_=0
    END IF
    IF(PRESENT(KID)) THEN
       CALL DEHIFE(FINT,FERR,H0_,EPS_,ILST_,FUNC,KID=KID)
    ELSE
       CALL DEHIFE(FINT,FERR,H0_,EPS_,ILST_,FUNC,KID='DEHIFT-FUNC')
    ENDIF
    RETURN
  END SUBROUTINE integ_f
END MODULE LIBT2
