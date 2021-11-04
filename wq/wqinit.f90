! wqinit.f90

MODULE wqinit

  PRIVATE
  PUBLIC wq_init

CONTAINS

  SUBROUTINE wq_init

    USE wqcomm_parm
    IMPLICIT NONE
    INTEGER:: i

    model_geometry = 0   ! id of geometry
                         !   0: rectangular, 1: cyrindrical
    xnmin    =  0.D0     ! minimum of normalized x: xn [in wave length]
    xnmax    = 10.D0     ! maximum of normalized x: xn [in wave length]
    ynmin    =  0.D0     ! minimum of normalized y: yn [in wave length]
    ynmax    = 10.D0     ! maximum of normalised y: xn [in wave length]

    B0       = 0.072d0   ! magnetic field [T]
    RR       = 0.22d0    ! major radius [m]
    RA       = 0.16d0    ! minor radius [m]
    q0       = 1.d0      ! safety factor on magnetic axis
    qa       = 3.d0      ! safety factor on plasma surface

    freq     = 5.0d9     ! 5GHz [Hz]    vacuum wave length = VC/freq [m]
    rkz = 0.D0           ! wave number in z direction [unit=vacuum wave length]
    nph = 0              ! mode number in phi direction (model_geometry=1)
    
    model_source = 1          ! id of source model
                              !   0:   no source
                              !   1-4: source on wall (-x,-y,+x,+y) S (Ez=0)
                              !   5-8: source on wall (-x,-y,+x,+y) P (Ex=Ey=0)
                              !   10:  initial source profile S (Ez=0)
                              !   10:  initial source profile P (Ex=Ey=0)
    source_position_xn = 0.5D0*(xnmin+xnmax) ! [unit=vacuum wave length]
    source_position_yn = 0.5D0*(ynmin+ynmax) ! [unit=vacuum wave length]
    source_width       = 4.D0   ! length of pulse [unit=vacuum wave length]
    source_thickness   = 2.D0   ! length of pulse [unit=vacuum wave length]
    source_angle       = 0.D0                ! [unit=degree, 0=direction of x]
    
    model_pulse      = 0      ! 0:flat, 1:pulse
    pulse_length     = 1.D0   ! length of pulse [unit=period=1/freq]
    model_ramp       = 0      ! 0:step function, 1:linear, 2:smooth
    ramp_length      = 1.D0   ! length of pulse  [period=1/freq]

    medium_max   = 1          ! number of medium (max=medium_m)
    id_medium(1) = 1          ! 0:vacuum, 1:plasma, 2:de const, 3:de resonance
    xnmin_medium(1) = xnmin   ! minimum xn of medium
    xnmax_medium(1) = xnmax   ! maximum xn of medium
    ynmin_medium(1) = ynmin   ! minimum yn of medium
    ynmax_medium(1) = ynmax   ! maximum yn of medium
    dielectric_medium(1)=1.D0 ! dielectric constant of medium
    res_freq_medium(1)=1.D0   ! resonance frequency of medium [in wave freq.]
    res_coll_medium(1)=0.1D0  ! collision frequency at resonenc [in wave freq]
    density_medium(1)=1.0D17  ! plasma density [m^{-3}]
    collision_medium(1)=0.003D0! collsion frequency normalize by wave freq.

    model_solver     = 1      ! 1:original, 2:matrix
    model_plot       = 0      ! 0:no plot save, 1:plot save

    fimplicit        = 1.D0   ! implicit factor (0: expl., 0.5:C-N, 1:impl.)
    ntype_mat        = 1      ! type of matrix solver
    eps_mat          = 1.D-10 ! convergence limit for matrix solver

!    dtfactor = 1.0d-4    ! ratio between time step and wave period
    dtfactor = 1.0d-2    ! ratio between time step and wave period
    dxfactor = 1.0d-1    ! ratio between mesh size and vacuum wave length in x
    dyfactor = 1.0d-1    ! ratio between mesh size and vacuum wave length in y

    ntmax    =    10     ! maximum of time step
                         !   maximum time = ntmax*dt_factor
                         !   normalized by wave period
    ntstep           = 1      ! nt interval for print status
    ngtstep          = 1      ! nt interval for storing global var to plot
    ngrstep          = 10     ! nt interval for storing profile to plot

    DO i=1,idebug_m
       idebuga(i)=0           ! debug options (41: nrank data, 61: mtx coef)
    END DO
    
    RETURN
  END SUBROUTINE wq_init
END MODULE wqinit

