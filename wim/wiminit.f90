MODULE wiminit

CONTAINS

! ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE wim_init

    USE wimcomm
    IMPLICIT NONE

    NZMAX=200           ! Number of mesh in parallel direction
    NWMAX=200           ! Number of mesh for integration
    NTMAX=100           ! Number of mesh for kernel function tabulation
    ZMIN=-300.D0        ! Minimum z position  (omega=omega_ce at z=0)
    ZMAX= 500.D0        ! Maximum z position
    TMAX=20.D0          ! Maximum argument of kernel functiion tabulation
    PN0=0.5D0           ! Plasma density (omega_pe^2/omega_ce^2)
    DBDZ=0.0005D0       ! Gradient of Magnetic field DB/B(0) DZ
    CER1=(0.D0,0.D0)    ! RHS wave electric field at z=0
    CEL1=(0.D0,0.D0)    ! LHS wave electric field at z=0
    CER2=(1.D0,0.D0)    ! RHS wave electric field at z=zmax
    CEL2=(1.D0,0.D0)    ! LHS wave electric field at z=zmax
    ANX=0.1D0           ! Perpendicular refractive index (k c / oemga)
    BETA=0.01D0         ! Plasma temperature (v_te/c)
    DZMAX=5.D0          ! Mesh accumulation factor
    DZWID=50.D0         ! Mesh accumulation width
    MODELW=0            ! 0: absorption boundary
                        ! 1: reflection boundary

    NTSET=0             ! tabulation setup (do no change)
    TMAXST=0.D0         ! tabulation setup (do no change)

    RETURN
  END SUBROUTINE wim_init
END MODULE wiminit
