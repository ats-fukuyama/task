! obinit.f90

MODULE obinit

CONTAINS

!     ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE ob_init

    USE obcomm_parm
    IMPLICIT NONE
    INTEGER:: nobt

!    modelg=3
    
    nobt_max=1                ! number of orbits
    nstp_max=10000            ! maximum number of orbit step
    ns_ob=2                   ! id of particle species
    lmax_nw=20                ! maximum number of iteration (initial condition)

    mdlobp=0                  ! model id of equation of motion
                              !     0: Eq of Motion  with Boozer coordinates
                              !     1: Eq of Motion  with Cylindrical coord.
    mdlobi=0                  ! model id of input scheme of initial parameters
                              !     0: penergy,pcangle,zeta,psipn,theta
                              !     1: penergy,pcangle,zeta,rr,zz
                              !   100: line input with psipn,theta
                              !   101: line input with rr,zz
    mdlobq=0                  ! model id of initial input parameters
                              !     0: 4th-order Runge-Kutta-Gill 
                              !     1: universal ODE solver (TBI)
                              !     2: symplectic solver (TBI)
    mdlobt=1                  ! model id of time normalization
                              !     0: real time
                              !     1: normalized by approximate bounce time
    mdlobc=0                  ! model id of one cycle calculation
                              !     0: independent of cycle, until tmax_ob
                              !     1: one cycle for trapped and untrapped
    mdlobw=0                  ! model id of output interval
                              !     0: no output
                              !     1: every step
                              !     2: every 10 step
                              !     3: every 100 step
                              !     4: every 1000 step
                              !     5: every 10000 step
    mdlobg=0                  ! model id of graphics
                              !     0: default
    mdlobx=1                  ! model id of wall
                              !     0: calculate only inside the wall
                              !            i.e. psip_ob<=psipa
                              !     1: continue Runge-Kutta
                              !            while satifying psip_ob>psipa
    tmax_ob=10.D0             ! maximum of orbit following time in omega_bounce
    delt_ob=0.1d0             ! time step size in omega_bounce
                              ! t_bounce = 2 Pi/ omega_bounce
                              ! omega bounce = (v_perp/qR) SQRT(r/2R)
                              ! omega_bounce^2 = (mu B /m)*(r/q^2 R^3)
    eps_ob=1.D-6              ! convergence criterion of orbit solution
    del_ob=1.D-4              ! step size of iteration (initial condition)
    eps_nw=1.D-6              ! convergence criterion of iteration (initial c.)

    penergy_ob_in(1)=1.D0        ! initial particle energy (mdlobi=0,1) [keV]
    pcangle_ob_in(1)=0.5D0       ! initial cosine of pitch angle (mdlobi=0,1)
    zeta_ob_in(1)=0.D0           ! initial toroidal angle (mdlobi=0,1) [degree]
    psipn_ob_in(1)=0.5D0         ! initial normalized poloidal flux (mdlobi=0)
    theta_ob_in(1)=0.D0          ! initial poloidal angle (mdlobi=0) [deg]
    rr_ob_in(1)=4.D0             ! initial major radius (mdlobi=1) [m]
    zz_ob_in(1)=0.D0             ! initial vertical position (mdlobi=1) [m]

! --- equilibirum parameters ---

    nrmax_ob=100   ! number of radial meshes
    nthmax_ob=64   ! number of radial meshes
    nsumax_ob=100  ! number of radial meshes

    RETURN
  END SUBROUTINE ob_init
END MODULE obinit
