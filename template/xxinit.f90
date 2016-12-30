MODULE xxinit

PUBLIC xx_init

CONTAINS

! ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE xx_init

    USE xxcomm
    IMPLICIT NONE

    nxmax=100    ! number of points
    x0=0.D0      ! initial points
    dx=0.1D0     ! distance between points
    am=1.D0      ! amplitude

    RETURN
  END SUBROUTINE xx_init
END MODULE xxinit
