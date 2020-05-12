! obinit.f90

MODULE obinit

CONTAINS

!     ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE OB_INIT

    USE obcomm_parm
    IMPLICIT NONE
    INTEGER:: NOBT

!     *********** Obt/Beam inital parameters ***********

!     RF     : WAVE FREQUENCY FOR OBT TRACING [MHZ]
!     RPI    : INITIAL MAJOR RADIUS R [M]
!     ZPI    : INITIAL VERTICAL POSITION Z [M]
!     PHII   : INITIAL TOROIDAL ANGLE [RADIAN]
!     RNZI   : INITIAL VERTICAL REFRACTIVE INDEX
!     RNPHII : INITIAL TOROIDAL REFRACTIVE INDEX
!     RKR0   : SPECULATED INITIAL RADIAL WAVE NUMBER [1/M]
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
      UUI    = 1.D0
      RCURVA = 0.D0
      RCURVB = 0.D0
      RBRADA = 0.03D0
      RBRADB = 0.03D0

      DO NOBT=1,NOBTM
         RFIN(NOBT)    = RF
         RPIN(NOBT)    = RPI
         ZPIN(NOBT)    = ZPI
         PHIIN(NOBT)   = PHII
         RNZIN(NOBT)   = RNZI
         RNPHIIN(NOBT) = RNPHII
         RKRIN(NOBT)   = RKR0
         MODEWIN(NOBT) = MODEW
         UUIN(NOBT)    = UUI
         RCURVAIN(NOBT)= RCURVA
         RCURVBIN(NOBT)= RCURVB
         RBRADAIN(NOBT)= RBRADA
         RBRADBIN(NOBT)= RBRADB
      ENDDO

!     *********** Obt/Beam control parameters ***********

!     NOBTMAX  : Number of obts
!     NSTPMAX  : Maxmum number of steps 
!     NRSMAX   : Number of minor radius division for absorbed power
!     NRRMAX   : Number of major radius division for absorbed power

!     SMAX   : MAXIMUM OBT LENGTH
!     DELS   : INCREMENTAL LENGTH OF OBT
!     UUMIN  : MINIMUM POWER TO ADVANCE OBT

!     EPSOBT : CONVERGENCE CRITEIRION IN OBT TRACING
!     DELOBT : MINIMUM STEP SIZE IN OBT TRACING
!     DELDER : STEP SIZE TO CALCULATE DERIVATIVES IN OBT TRACING

!     DELKR  : STEP SIZE TO ESTIMATE D/DKR IN NEWTON METHOD
!     EPSNW  : CONVERGENCE CRITEIRION IN NEWTON METHOD
!     LMAXNW : MAXIMUM ITERATION COUNT IN NEWTON METHOD

!     MDLOBI : INPUT TYPE OF WAVE PARAMETERS
!              0 : RF,RP,ZP,PHI,RKR0,RNZ,RNPHI,UU
!              1 : RF,RP,ZP,PHI,RKR0,ANGZ,ANGPH,UU
!             11 : RF,RP,ZP,RKR0,RNZ,RNPHI,UU
!            100 : RFIN,RPIN,ZPIN,PHIIN,RKRIN,RNZIN,RNPHIIN,UUIN: namelist
!            101 : RFIN,RPIN,ZPIN,PHIIN,RKRIN,ANGZIN,ANGPHIN,UUIN: namelist

!     MDLOBG : TYPE OF GRAPHICS
!              0 : FULL TORUS, FULL RADIUS FOR DEPOSITION
!              1 : PARTIAL TORUS, FULL RADIUS FOR DEPOSITION
!              2 : FULL TORUS, PARTIAL RADIUS FOR DEPOSITION
!              3 : PARTIAL TORUS, PARTIAL RADIUS FOR DEPOSITION
!             11 : 2D plane

!     MDLOBP : TYPE OF Power deoposition profile GRAPHICS
!              0 : Sum up of power deposition of each species
!              1 : Power flux with respect to minor radius
!              2 : Power flud with respect to major radius

!     MDLOBQ : TYPE OF DIFFERENTIAL EQUATION IN OBT TRACING
!              0 : RUNGE-KUTTA, FIXED STEPSIZE
!              1 : RUNGE-KUTTA, FIXED STEPSIZE, with k_X adjust to satisfy D=0
!              2 : RUNGE-KUTTA, FIXED STEPSIZE, with mode conversion
!              3 : RUNGE-KUTTA, VARIABLE STEPSIZE
!              4 : SYMPLECTIC METHOD, FIXED STEPSIZE (not completed)

!     MDLOBW : Level of PRINT OUTPUT
!              0 : NO output
!             -1 : Write initial kr calculation
!              1 : Write data every step
!              2 : Write data every 10 steps
!              3 : Write data every 100 step
!              4 : Write data every 1000 step
!              5 : Write data every 10000 step

      SMAX   = 1.0D0
      DELS   = 0.05D0
      UUMIN  = 1.D-4

      EPSOBT = 1.D-4
      DELOBT = 1.D-3
      DELDER = 1.D-4

      DELKR  = 1.D0
      EPSNW  = 1.D-6
      LMAXNW = 100

      NOBTMAX  = 1
      NSTPMAX  = 10000
      NRSMAX   = 100
      NRRMAX   = 200

      MDLOBI = 0
      MDLOBG = 0
      MDLOBP = 1
      MDLOBQ = 1
      MDLOBW = 0

      RETURN
  END SUBROUTINE OB_INIT
END MODULE obinit
