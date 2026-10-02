! testsigma.f90

MODULE sigmav_nrl_dt_comm
  USE bpsd_kinds
  REAL(rkind):: T_bulk_dt
END MODULE sigmav_nrl_dt_comm

PROGRAM testsigma
  USE bpsd_kinds
  USE bpsd_constants
  USE libsigma
  USE libgrf
  IMPLICIT NONE

  REAL(rkind),DIMENSION(10):: x
  REAL(rkind),DIMENSION(10,4):: y
  INTEGER:: i

  CALL GSOPEN
  
  CALL set_usigmavmal_dt
  
  DO i=1,10
     x(i)=LOG10(tempa(i))
     y(i,1)=LOG10(sva(i))
     y(i,2)=LOG10(sigmavm_int_dt(tempa(i)))
     y(i,3)=LOG10(sigmavm_low_dt(tempa(i)))
     y(i,4)=LOG10(sigmavm_nrl_dt(tempa(i)))
     WRITE(6,'(I4,5ES12.4)') i,x(i),y(i,1),y(i,2),y(i,3),y(i,4)
  END DO

  CALL PAGES
  CALL grd1d(0,x,y,10,10,4,'@temp vs sigmav@',3)
  CALL PAGEE

  CALL GSCLOS
  STOP

CONTAINS

  ! analytic sigmavm

  FUNCTION sigmavm_low_dt(temp)

    IMPLICIT NONE
    REAL(rkind), INTENT(IN):: temp
    REAL(rkind):: sigmavm_low_dt
    REAL(rkind):: temp3

    IF(temp.LE.0.D0) THEN
       sigmavm_low_dt=0.D0
       RETURN
    END IF
    temp3=1.D0/temp**(1.D0/3.D0)
!    sigmavm_low_dd=2.33D-14*temp3**2*EXP(-18.76D0*temp3)
    sigmavm_low_dt=3.68D-12*temp3**2*EXP(-19.94D0*temp3)
    RETURN
  END FUNCTION sigmavm_low_dt

END PROGRAM testsigma
      
