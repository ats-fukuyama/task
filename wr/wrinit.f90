! wrinit.f90

MODULE wrinit

CONTAINS

!     ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE WR_INIT

    USE wrcomm_parm
    IMPLICIT NONE
    INTEGER:: NRAY

!     *********** Ray/Beam inital parameters ***********

!     RF     : WAVE FREQUENCY FOR RAY TRACING [MHZ]
!     RPI    : INITIAL MAJOR RADIUS R [M]
!     ZPI    : INITIAL VERTICAL POSITION Z [M]
!     PHII   : INITIAL TOROIDAL ANGLE [RADIAN]
!     RNZI   : INITIAL VERTICAL REFRACTIVE INDEX
!     RNPHII : INITIAL TOROIDAL REFRACTIVE INDEX
!     RKR0   : SPECULATED INITIAL RADIAL WAVE NUMBER [1/M]
!     MODEW  : Selection of initial k_r
!              0: found RKR near RKR0
!              1: larger  RKR (negative, inward for forward wave)    
!              2: smaller RKR (negative, inward for forward wave)    
!             -1: larger RKR  (positive, inward for backward wave)    
!             -2: smaller RKR (positive, inward for backward wave)    
!     UUI    : INITIAL WAVE ENERGY

!     RCURVA : INITIAL WAVE-FRONT CURVATURE RADIUS (0 for Plane wave)
!     RCURVB : INITIAL WAVE-FRONT CURVATURE RADIUS (0 for Plane wave)
!                 RCURVA perp to k and B
!                 RCURVB perp to k and in kxB plane
!     RBRADA  : INITIAL BEAM RADIUS
!     RBRADB  : INITIAL BEAM RADIUS
!                 RBRADA perp to k and B
!                 RBRADB perp to k and in kxB plane

      RF     = 170.D3
      RPI    = 3.95D0
      ZPI    = 0.D0
      PHII   = 0.D0
      RNZI   = 0.D0
      RNPHII = 0.5D0
      RKR0   = -1000.D0
      MODEW  = 0
      UUI    = 1.D0
      RCURVA = 0.D0
      RCURVB = 0.D0
      RBRADA = 0.03D0
      RBRADB = 0.03D0

      DO NRAY=1,NRAYM
         RFIN(NRAY)    = RF
         RPIN(NRAY)    = RPI
         ZPIN(NRAY)    = ZPI
         PHIIN(NRAY)   = PHII
         RNZIN(NRAY)   = RNZI
         RNPHIIN(NRAY) = RNPHII
         RKRIN(NRAY)   = RKR0
         MODEWIN(NRAY) = MODEW
         UUIN(NRAY)    = UUI
         RCURVAIN(NRAY)= RCURVA
         RCURVBIN(NRAY)= RCURVB
         RBRADAIN(NRAY)= RBRADA
         RBRADBIN(NRAY)= RBRADB
      ENDDO

!     *********** Ray/Beam control parameters ***********

!     NRAYMAX  : Number of rays
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
!              0 : RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU
!              1 : RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU
!             11 : RF,RP,ZP,RKR0,RNZ,RNPHI,UU
!            100 : RFIN,RPIN,ZPIN,PHIIN,RKRIN,RNZIN,RNPHIIN,UUIN: namelist
!            101 : RFIN,RPIN,ZPIN,PHIIN,RKRIN,ANGZIN,ANGPHIN,UUIN: namelist

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

      NRAYMAX  = 1
      NSTPMAX  = 10000
      NRSMAX   = 100
      NRLMAX   = 200

      MDLWRI = 0
      MDLWRG = 0
      MDLWRP = 1
      MDLWRQ = 1
      MDLWRW = 0

      nres_type = 0
      nres_max = 3

      ! pne_threshold: threshold value to identify to be in plasma [10^20 m^-3]

      pne_threshold=1.D-6

      ! Rmax_wr: maximum of R for ray calculation, if 0, set RR+1.2*RA
      ! Rmin_wr: minimum of R for ray calculation, if 0, set RR-1.2*RA >1.D-6
      ! Zmax_wr: maximum of Z for ray calculation, if 0, set RR+1.2*rkap*RA
      ! Zmin_wr: minimum of Z for ray calculation, if 0, set RR-1.2*rkap*RA

      Rmax_wr=0.D0
      Rmin_wr=0.D0
      Zmax_wr=0.D0
      Zmin_wr=0.D0

! --- defined in dp ---

      NPMAX_DP=50
      NTHMAX_DP=50
      NRMAX_DP=50

      RETURN
  END SUBROUTINE WR_INIT
END MODULE wrinit
