! wrinit.f90

MODULE wrinit

CONTAINS

!     ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE WR_INIT

    USE wrcomm_parm
    IMPLICIT NONE
    INTEGER:: NRAY,i,nsa

!     *********** Ray/Beam inital parameters ***********

!     NRAYMAX : Number of rays

!     RFIN    : WAVE FREQUENCY FOR RAY TRACING [MHZ]
!     RPIN    : INITIAL MAJOR RADIUS R [M]
!     ZPIN    : INITIAL VERTICAL POSITION Z [M]
!     PHIIN   : INITIAL TOROIDAL ANGLE [RADIAN]
!     ANGTIN  : INITIAL toroidal injection angle [degree] from equi-phi plane
!     ANGPIN  : INITIAL poloidal injection angle [degree] from horizonal plane
!     RNPHIN  : INITIAL poloidal REFRACTIVE INDEX
!     MODEWIN : Selection of initial wave mode
!               1: slow wave: larger k  (negative, inward for forward wave)
!               2: fast wave: smaller k (negative, inward for forward wave)
!     RNKIN   : INITIAL estimate of initial refractive index (1.D0 for vaccum)
!               positive for forward wave propagating inward 
!               negative for backward wave propagating inward 
!     UUIN    : INITIAL WAVE POWER (default 1.D0)

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
         ANGTIN(NRAY)   = 0.D0
         ANGPIN(NRAY)   = 0.D0
         RNPHIN(NRAY)   = 0.0D0
         RNZIN(NRAY)    = 0.0D0
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

!     EPSRAY : CONVERGENCE CRITERION IN RAY TRACING
!     DELRAY : MINIMUM STEP SIZE IN RAY TRACING
!     DELDER : STEP SIZE TO CALCULATE DERIVATIVES IN RAY TRACING
!     EPSD0  : CONVERGENCE CRITERIION for D=0 in WRMODNWTN

!     DELKR  : STEP SIZE TO ESTIMATE D/DKR IN NEWTON METHOD
!     EPSNW  : CONVERGENCE CRITEIRION IN NEWTON METHOD
!     LMAXNW : MAXIMUM ITERATION COUNT IN NEWTON METHOD

!     mode_beam : 0 for ray tracing, nonzero for beam tracing

      mode_beam=0

!     MDLWRI : INPUT TYPE OF WAVE PARAMETERS
!                   RFIN,RPIN,ZPIN,PHIIN,...,MODEWIN,RNKIN,UUIN
!         1, 2, 3 : ANGTIN,ANGPIN (toroidal, poloidal incident angles) 
!        11,12,13 : ANGTIN,ANGPIN (toroidal, poloidal absolute angles : TRAVIS)
!        21,22,23 : RNPHIN,ANGPIN (angt=ARCSIN(rnph/(rnk*COS(angp))))
!        31,32,33 : RNPHIN,ANGPIN (angt=ARCSIN(rnph/(rnk*COS(angp))) : TRAVIS)
!        41,42,43 : RNPHIN,RNZIN  (angp=   ARCSIN(rnz/rnk) LFS)
!        51,52,53 : RNPHIN,RNZIN  (angp=pi+ARCSIN(rnz/rnk) HFS)
!               1 : poloidal first angle  k_p = k sin angp
!               2 : toroidal first angle  k_t = k sin angt
!               3 : Intuitive angle       k_p = k sin angp, k_t = k sin angt
!            +100 : interactive parameter input for each ray
!      
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

      !     model_fdrv: type of driver to calculate del D/del r, del D/del k
      !        1 : original driver
      !        2 : new driver on 2023/09/03
      !        3 : new driver on 2023/11/21
      
      !     model_fdrv_ds: definition of ds in driver wrfdrv
      !        0 : original    DS
      !        1 : corrected   DOMG
      
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
      DELDER = 1.D-6

      EPSD0  = 1.D-4

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

      model_fdrv = 3
      model_fdrv_ds = 0

      nres_type = 0
      nres_max = 3

      ! pne_threshold: threshold value to identify to be in plasma
      !                          [omegape^2/omega^2]
      ! bdr_threshold: threshold value to boundary factor: F in below

      pne_threshold=1.D-3
      bdr_threshold=1.5

      ! For F=bdr_threshold.GT.0.D0
      !    Rmax_wr: maximum of R for ray calculation, rmax_wr=RR+F*RA
      !    Rmin_wr: minimum of R for ray calculation, rmin_wr,RR-F*RA
      !    Zmax_wr: maximum of Z for ray calculation, zmax_wr, F*rkap*RA
      !    Zmin_wr: minimum of Z for ray calculation, zmin_wr,-F*rkap*RA
      ! For F=bdr_threshold.EQ.0.D0
      !    given Rmax_wr,Rmin_wr,Zmax_wr,Zmin_wr are used

      Rmax_wr=RR+bdr_threshold*RA
      Rmin_wr=RR-bdr_threshold+RA
      Zmax_wr= bdr_threshold*RKAP*RA
      Zmin_wr=-bdr_threshold*RKAP*RA

      ! ra_wr:   maximum minor radius for plot
      ra_wr=1.1D0

      ! wr ascii data file

      KNAMWRW='wr-text-data'

      ! --- debug output contral ---
      !        idebug_wr =  1: initial position, wave number, denisty
      !        idebug_wr =  2: vacuum step position, density
      !        idebug_wr =  3: initial plasma position, wave number
      !        idebug_wr =  4: initial plasma position, k_para, kperp 
      !        idebug_wr =  5: initial plasma position, k convergence check
      !        idebug_wr =  6: wr_newton: iteration rk
      !        idebug_wr =  7: wr_newton: polarization check
      !        idebug_wr =  8: plasma step: position, wave number
      !        idebug_wr = 10: wr_write_line: output Y(1:6)
      !        idebug_wr = 11: wrrkft: initial and each step
      !        idebug_wr = 12: wrfdrv: functions and derivatives
      !        idebug_wr = 13: wrmodenwtn: iteration data
      !        idebug_wr = 90: wrgout: wrgrf5: R,cexyz,ceoxp, error
      !        idebug_wr = 91: wrgout: wrgrf6: nres_max
      !        idebug_wr = 92: wrgout: wrgrf2: gxorg,gxstep,gystep
      
      
      DO i=1,idebug_max
         idebug_wr(i)=0
      END DO

! --- defined in dp ---

      NPMAX_DP=50
      NTHMAX_DP=50
      NRMAX_DP=50
      NSAMAX_WR=2
      DO NSA=1,NSM
         NS_NSA_WR(NSA)=NSA
      END DO
      nsa_grf=1

      RETURN
  END SUBROUTINE WR_INIT
END MODULE wrinit
