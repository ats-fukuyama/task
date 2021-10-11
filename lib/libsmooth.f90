! libsmooth.f90

MODULE libsmooth

  ! *** library of smmoth funcutions ***

  USE task_kinds
  
  PRIVATE
  PUBLIC smooth_1

CONTAINS

  FUNCTION smooth_1(x)
    REAL(rkind),INTENT(IN):: x
    REAL(rkind):: smooth_1

    smooth_1=x**2/(x**2+(1.D0-x)**2)
    RETURN
  END FUNCTION smooth_1

END MODULE libsmooth
       

  
