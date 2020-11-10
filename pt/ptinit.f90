! ptinit.f90

MODULE ptinit

  PRIVATE
  PUBLIC pt_init

CONTAINS

  SUBROUTINE pt_init
    USE ptcomm_parm
    IMPLICIT NONE

    bb=5.3D0
    rr=6.2D0
    ra=2.0D0
    rkap=1.7D0
    rdlt=0.33D0

    ngxmax=201
    ngymax=201
    nthmax=201

    RETURN
  END SUBROUTINE pt_init
END MODULE ptinit
