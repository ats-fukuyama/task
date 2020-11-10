! ptview.f90

MODULE ptview

CONTAINS

  !     ****** SHOW PARAMETERS ******

  SUBROUTINE pt_view

    USE ptcomm_parm
    IMPLICIT NONE

    WRITE(6,'(A,ES12.4)') 'bb                     =',bb
    WRITE(6,'(A,ES12.4)') 'rr                     =',rr
    WRITE(6,'(A,ES12.4)') 'ra                     =',ra
    WRITE(6,'(A,ES12.4)') 'rkap                   =',rkap
    WRITE(6,'(A,ES12.4)') 'rdlt                   =',rdlt
    WRITE(6,'(A,I8)')     'ngxmax                 =',ngxmax
    WRITE(6,'(A,I8)')     'ngymax                 =',ngymax
    WRITE(6,'(A,I8)')     'nthmax                 =',nthmax
    RETURN
  END SUBROUTINE pt_view
END MODULE ptview
