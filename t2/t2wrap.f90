!C
!C
!C 
!C  TASK/PL <=> TASK/T2
!C
!C                                                                      1

MODULE T2WRAP

  IMPLICIT NONE
CONTAINS

  SUBROUTINE T2WRAP_PL_TO_T2
    
    USE PLINIT
    
    CALL PL_INIT
       pl_prof(
       RETURN
  END SUBROUTINE T2WRAP_PL_TO_T2

  
END MODULE T2WRAP
