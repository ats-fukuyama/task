!     $Id$

MODULE wiinit

CONTAINS

! ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE wi_init

    USE wicomm
    IMPLICIT NONE

    MODEWI=1           ! calculation mode: 1:unmag
    NXMAX=400          ! number of grid points in the x direction
    NWMAX=200          ! number of grid points for integration of kernel
    XMAX=200.D0        ! maximum value of x (minimum value is 0.0)
    PN0=2.0D0          ! normalized plasma density (omega_pe^2/omega^2) at x=0
    ALFA=0.005D0       ! normalized density gradient [ n=n_0 exp(-ALFA*x) ]
    CFYN=(1.D0,0.D0)   ! Electric field of the incident wave at nx=nxmax
    AKY=0.2D0          ! refractive index in the y direction (k_y c/omega)
    BETA=0.1D0         ! ratio of electron thermal velocity to light velocity

    RETURN
  END SUBROUTINE wi_init
END MODULE wiinit
