! wrinit.f90

MODULE wrinit

CONTAINS

!     ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE WR_INIT

    USE wrcomm_parm
    IMPLICIT NONE
    INTEGER:: NRAY,i

!     *********** Ray/Beam inital parameters ***********

!     NRAYMAX : Number of rays

!     RFIN    : WAVE FREQUENCY FOR RAY TRACING [MHZ]
!     RPIN    : INITIAL MAJOR RADIUS R [M]
!     ZPIN    : INITIAL VERTICAL POSITION Z [M]
!     PHIIN   : INITIAL TOROIDAL ANGLE [RADIAN]
!     RNKPIN  : INITIAL poloidal REFRACTIVE INDEX
!     RNKTIN  : INITIAL toroidal REFRACTIVE INDEX
!     ANGPIN  : INITIAL poloidal injection angle [degree] from horizonal plane
!     ANGTIN  : INITIAL toroidal injection angle [degree] from equi-phi plane
!     MODEWIN : Selection of initial wave mode
!              0: slow wave: larger k  (negative, inward for forward wave)
!              1: fast wave: smaller k (negative, inward for forward wave)
!              2: slow wave: larger k  (positive, inward for backward wave)
!              3: fast wave: smaller k (positive, inward for backward wave) 
!     RNKIN   : INITIAL estimate of inirial refractive index (1.D0 for vaccum)
!     UUIN    : INITIAL WAVE ENERGY

!     RCURVA : INITIAL WAVE-FRONT CURVATURE RADIUS (0 for Plane wave)
!     RCURVB : INITIAL WAVE-FRONT CURVATURE RADIUS (0 for Plane wave)
!                 RCURVA perp to k and B
!                 RCURVB perp to k and in kxB plane
!     RBRADA  : INITIAL BEAM RADIUS
!     RBRADB  : INITIAL BEAM RADIUS
!                 RBRADA perp to k and B
!                 RBRADB perp to k and in kxB plane

      NRAYMAX  = 1

      DO NRAY=1,NRAYM
         RFIN(NRAY)     = 170.D3
         RPIN(NRAY)     = 3.95D0
         ZPIN(NRAY)     = 0.D0
         PHIIN(NRAY)    = 0.D0
         ANGPIN(NRAY)   = 0.D0
         ANGTIN(NRAY)   = 0.D0
         RNKPIN(NRAY)   = 0.0D0
         RNKTIN(NRAY)   = 0.5D0
         MODEWIN(NRAY)  = 1
         rnkin(nray)    = 1.d0
         UUIN(NRAY)     = 1.D0
         RCURVAIN(NRAY) = 0.D0
         RCURVBIN(NRAY) = 0.D0
         RBRADAIN(NRAY) = 0.03D0
         RBRADBIN(NRAY) = 0.03D0
      ENDDO

!     *********** Ray/Beam control parameters ***********

!     NSTPMAX  : Maxmum number of steps 
!     NRSMAX   : Number of minor radius division for absorbed power
!     NRLMAX   : Number of major radius division for absorbed power

!     SMAX   : MAXIMUM RAY LENGTH
!     DELS   : INCREMENTAL LENGTH OF RAY
!     UUMIN  : MINIMUM POWER TO ADVANCE RAY

!     EPSRAY : CONVERGENCE CRITEIRION IN RAY TRACING
!     DELRAY : MINIMUM STEP SIZE IN RAY TRACING
!     DELDER : STEP SIZE TO CALCULATE DERIVATIVES IN RAY TRACING

!     DELKR  : STEP SIZE TO ESTIMATE D/DKR IN NEWTON METHOD
!     EPSNW  : CONVERGENCE CRITEIRION IN NEWTON METHOD
!     LMAXNW : MAXIMUM ITERATION COUNT IN NEWTON METHOD

!     mode_beam : 0 for ray tracing, nonzero for beam tracing

      mode_beam=0

!     MDLWRI : INPUT TYPE OF WAVE PARAMETERS
!        0,1 : RFIN,RPIN,ZPIN,PHIIN,ANGPIN,ANGTIN,MODEWIN,UUIN: interactive
!        2,3 : RFIN,RPIN,ZPIN,PHIIN,RNPIN,RNTIN,MODEWIN,UUIN:   interactive
!    100,101 : RFIN,RPIN,ZPIN,PHIIN,ANGPIN,ANGTIN,MODEWIN,UUIN: namelist
!    102,103 : RFIN,RPIN,ZPIN,PHIIN,RNPIN,RNTIN,MODEWIN,UUIN:   namelist
!    0,2,100,102: poloidal first definition: k_p = k sin angp
!    1,3,101,103: toridal first definition:  k_t = k sin angt
      
!     MDLWRG : TYPE OF GRAPHICS
!              0 : FULL TORUS, FULL RADIUS FOR DEPOSITION
!              1 : PARTIAL TORUS, FULL RADIUS FOR DEPOSITION
!              2 : FULL TORUS, PARTIAL RADIUS FOR DEPOSITION
!              3 : PARTIAL TORUS, PARTIAL RADIUS FOR DEPOSITION
!             11 : 2D plane

!     MDLWRP : TYPE OF Power deoposition profile GRAPHICS
!              0 : Sum up of power deposition of each species
!              1 : Power flux with respect to minor radius
!              2 : Power flud with respect to major radius

!     MDLWRQ : TYPE OF DIFFERENTIAL EQUATION IN RAY TRACING
!              0 : RUNGE-KUTTA, FIXED STEPSIZE
!              1 : RUNGE-KUTTA, FIXED STEPSIZE, with k_X adjust to satisfy D=0
!              2 : RUNGE-KUTTA, FIXED STEPSIZE, with mode conversion
!              3 : RUNGE-KUTTA, VARIABLE STEPSIZE
!              4 : SYMPLECTIC METHOD, FIXED STEPSIZE (not completed)

!     MDLWRW : Level of PRINT OUTPUT
!              0 : NO output
!             -1 : Write initial kr calculation
!              1 : Write data every step
!              2 : Write data every 10 steps
!              3 : Write data every 100 step
!              4 : Write data every 1000 step
!              5 : Write data every 10000 step

!     nres_type : plot type of resonance curves
!              0 : power abs density (max, 50% for nres_max=3)
!              1 : power flux        (25%, 50%, 75% for nres_max=3)
!              2 : ray length        (25%, 50%, 75% for nres_max=3)
!     nres_max : number of resonance curves

      SMAX   = 1.0D0
      DELS   = 0.05D0
      UUMIN  = 1.D-4

      EPSRAY = 1.D-4
      DELRAY = 1.D-3
      DELDER = 1.D-4

      DELKR  = 1.D0
      EPSNW  = 1.D-6
      LMAXNW = 100

      NSTPMAX  = 10000
      NRSMAX   = 100
      NRLMAX   = 200

      MDLWRI = 1
      MDLWRG = 0
      MDLWRP = 1
      MDLWRQ = 1
      MDLWRW = 0

      nres_type = 0
      nres_max = 3

      ! pne_threshold: threshold value to identify to be in plasma [10^20 m^-3]
      ! bdr_threshold: threshold value to boundary factor: F in below

      pne_threshold=1.D-6
      bdr_threshold=1.2

      ! Rmax_wr: maximum of R for ray calculation, if 0, set RR+F*RA
      ! Rmin_wr: minimum of R for ray calculation, if 0, set RR-F*RA >1.D-6
      ! Zmax_wr: maximum of Z for ray calculation, if 0, set RR+F*rkap*RA
      ! Zmin_wr: minimum of Z for ray calculation, if 0, set RR-F*rkap*RA

      Rmax_wr=0.D0
      Rmin_wr=0.D0
      Zmax_wr=0.D0
      Zmin_wr=0.D0

      DO i=1,idebug_max
         idebug_wr(i)=0
      END DO

! --- defined in dp ---

      NPMAX_DP=50
      NTHMAX_DP=50
      NRMAX_DP=50

      RETURN
  END SUBROUTINE WR_INIT
END MODULE wrinit
