! obinit.f90

MODULE obinit

CONTAINS

!     ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE ob_init

    USE obcomm_parm
    IMPLICIT NONE
    INTEGER:: nobt

    nobt_max=1                ! number of orbits
    nstp_max=100              ! maximum number of orbits
    ns_ob=2                   ! id of particle species
    lmax_nw=20                ! maximum number of iteration (initial condition)

    mdlobp=0                  ! model id of equation of motion
                              !     0: Eq of Motion  with Boozer coordinates
                              !     1: Eq of Motion  with Cylindrical coord.
    mdlobi=0                  ! model id of input scheme of initial parameters
                              !     0: penergy,pangle,zeta,pzeta,theta
                              !     1: penergy,pangle,zeta,rr,zz
                              !   100: line input with pzeta,theta
                              !   101: line input with rr,zz
    mdlobq=0                  ! model id of initial input parameters
                              !     0: 4th -order Runge-Kutta-Gill 
    mdlobg=0                  ! model id of graphics
                              !     0: default

    smax=10.D0                ! maximum of orbit length
    dels=0.01d0               ! step size of orbit length
    eps_obt=1.D-6             ! convergence criterion of orbit solution
    del_obt=1.D-4             ! step size of iteration (initial condition)
    eps_nw=1.D-6              ! convergence criterion of iteration (initial c.)

    penergy_in(1)=1.D0        ! initial particle energy (mdlobi=0,1) [keV]
    pangle_in(1)=1.D0         ! initial sine of pitch angle (mdlobi=0,1) [mu/E]
    zeta_in(1)=0.D0           ! initial toroidal angle (mdlobi=0,1) [degree]
    pzeta_in(1)=-0.5D0        ! initial toroidal momentum (mdlobi=0) [P/E]
    theta_in(1)=0.D0          ! initial poloidal angle (mdlobi=0) [deg]
    rr_in(1)=4.D0             ! initial major radius (mdlobi=1) [m]
    zz_in(1)=0.D0             ! initial vertical position (mdlobi=1) [m]

    RETURN
  END SUBROUTINE ob_init
END MODULE obinit
