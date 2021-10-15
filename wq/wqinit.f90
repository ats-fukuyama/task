! wqinit.f90

MODULE wqinit

  PRIVATE
  PUBLIC wq_init

CONTAINS

  SUBROUTINE wq_init

    USE wqcomm_parm
    IMPLICIT NONE

    freq     = 5.0d9     ! 5GHz [Hz]
    dtfactor = 1.0d-4    ! ratio between time step and wave period
    dxfactor = 1.0d-1    ! ratio between mesh size and wave length in x
    dyfactor = 1.0d-1    ! ratio between mesh size and wave length in y
    nufactor = 1.0d-3    ! ratio between collision freq. and wave freq.
    B0       = 0.072d0   ! magnetic field [T]
    RR       = 0.22d0    ! major radius [m]
    RA       = 0.16d0    ! minor radius [m]
    q0       = 1.d0      ! safety factor on magnetic axis
    qa       = 3.d0      ! safety factor on plasma surface
    n0       = 1.0d17    ! plasma density [m^{-3}]
    ntmax    = 10        ! maximum time step
    nxmax    = 61        ! number of meshes in x direction
    nymax    = 61        ! number of meshes in y direction
    INMODE   = 1         !
    TMN      = 1.0d0     !

    source_width     = 2.D0   ! length of pulse [wavelength=VC/freq]
    model_pulse      = 0      ! 0:constant, 1:one pulse
    pulse_length     = 1.D0   ! length of pulse [period=1/freq]
    model_ramp       = 0      ! 0:step function, 1:linear, 2:smooth
    ramp_length      = 1.D0   ! length of pulse  [period=1/freq]
    model_dielectric = 1      ! 0:vacuum, 1:plasma, 2:de const, 3:de resonance
    model_solver     = 0      ! 0:original, 1:explicit, 2:implicit
    model_plot       = 0      ! 0:no plot save, 1:plot save
    dielectric_2     = 1.D0   ! dielectric constant in a layer for model_de=2
    dielectric_3     = 1.D0   ! dielectric constant in a layer for model_de=3
    freq_resonance   = 1.D0   ! resonance frequency for model_de=3 [freq]
    freq_collision   = 0.003  ! collision frequency for model_de=3 [freq]
    ntplot_interval  = 1000   ! nt interval for plot save
    ntplot_max       = 5      ! total number of plot save 
    
    RETURN
  END SUBROUTINE wq_init
END MODULE wqinit

