! dpinit.f90

MODULE dpinit

  PRIVATE
  PUBLIC dp_init

CONTAINS

  SUBROUTINE dp_init

    USE dpcomm_parm_local
    IMPLICIT NONE
    INTEGER:: NS

!
!     *********** INPUT PARAMETERS fof DP library ***********
!
!     --- PLASMA MODEL PARAMETERS ---
!
!     MODELP(NS): TYPE OF ANALYTIC DIELECTRIC TENSOR
!                 0 : COLLISIONLESS COLD MODEL
!                 1 : COLLISIONAL COLD MODEL
!                 2 : IDEAL MHD MODEL
!                 3 : RESISTIVE MHD MODEL
!                 4 : KINETIC MODEL WITHOUT FLR
!                 5 : KINETIC MODEL WITH FLR (symmetric)
!                 6 : KINETIC MODEL WITH FLR
!                 7 : Kinetic with parallel drift and FLR by Swanson
!                 8 : Kinetic with parallel drift and FLR by T. Watanabe
!                 9 : KINETIC MODEL WITH RELATIVISTIC EFFECTS (test)
!                11 : (WM) MHD plasma
!                12 : (WM) Cold plasma
!                13 : (WM) Hot plasma (No FLR)
!                14 : (WM) Hot plasma (Cold FLR)
!                15 : (WM) Hot plasma (FLR)
!                16 : (WM) Drift kinetic plasma
!                21 : DRIFTKINETIC MODEL (coming)
!                22 : GYROKINETIC MODEL (coming)
!                31 : Cold beam model (FLR)
!
!             0- 99 : PROPAGATION  = GIVEN MODEL
!                     POLARIZATION = GIVEN MODEL
!                     ABSORPTION   = GIVEN MODEL
!
!           100-199 : PROPAGATION  = COLD
!                     POLARIZATION = GIVEN MODEL
!                     ABSORPTION   = GIVEN MODEL
!
!           200-299 : PROPAGATION  = COLD
!                     POLARIZATION = COLD
!                     ABSORPTION   = GIVEN MODEL
!
!           300-399 : PROPAGATION  = KINETIC 4
!                     POLARIZATION = GIVEN MODEL
!                     ABSORPTION   = GIVEN MODEL
!
!           400-499 : PROPAGATION  = KINETIC 4
!                     POLARIZATION = KINETIC 4
!                     ABSORPTION   = GIVEN MODEL
!
!           500-599 : PROPAGATION  = KINETIC 6
!                     POLARIZATION = GIVEN MODEL
!                     ABSORPTION   = GIVEN MODEL
!
!           600-699 : PROPAGATION  = KINETIC 6
!                     POLARIZATION = KINETIC 6
!                     ABSORPTION   = GIVEN MODEL
!
!     MODELV(NS) : NUMERICAL MODEL (*: not yet implemented)
!              0 : ANALYTIC MODEL
!              1 : KINETIC: ANALYTIC MAXWELLIAN DISTRIBUTION
!              2 : KINETIC: READ FPDATA DISTRIBUTION
!              3 : KINETIC: ANALYTIC MAXWELLIAN DISTRIBUTUION (RELATIVISTIC)
!              4 : KINETIC: READ FPDATA DISTRIBUTION (RELATIVISTIC)
!              5 : *DRIFTKINETIC: ANALYTIC MAXWELLIAN DISTRIBUTION 
!              6 : *DRIFTKINETIC: READ FPDATA DISTRIBUTION
!              9 : *LOCAL MODEL (MODELV locally specified by MODELVR)
!
!     NCMIN(NS): MINIMUM HARMONIC NUMBER
!     NCMAX(NS): MAXMUM  HARMONIC NUMBER

    DO NS=1,NSM
       MODELP(NS)= 0
       MODELV(NS)= 0
       NCMIN(NS)=-2
       NCMAX(NS)= 2
    ENDDO

!     MODEL_ES : 0: Electromagnetic wave, 1:Electrostatic wave
    MODEL_ES=0
    
!     EPSRT  : CONVERGENCE CRITERION OF ROOT FINDING
!     LMAXRT : MAXIMUM ITERATION COUNT OF ROOT FINDING

    EPSRT  = 1.D-12
    LMAXRT = 20

!     --- Velocity distribution function parameters ---
!             --- usually read from fpfile ---
!
!     NS_NSA_DP(NSA): particle species of NSA
!     PMAX_dp(NSA)  : maximum momentum normalized by p_thermal
!     EMAX_dp(NSA)  : maximum energy in keV, if EMAX is not zero
!     rhon_min(NSA) : minimum radius of velocity distribution function (r/a)
!     rhon_max(NSA) : maximum radius of velocity distribution function (r/a)
!
!     NPMAX_DP : number of momentum magnitude mesh
!     NTHMAX_DP: number of momentum angle mesh
!     NRMAX_DP : number of radial mesh
!     NSAMAX_DP: number of test particle species

    DO NS=1,NSM
       NS_NSA_DP(NS)=NS
       PMAX_dp(NS)= 7.D0
       EMAX_dp(NS)= 7.D0
       rhon_min(NS)=0.D0
       rhon_max(NS)=1.D0
    ENDDO

    NPMAX_DP=100
    NTHMAX_DP=100
    NRMAX_DP=3
    NSAMAX_DP=2

!     --- Root finding and dispersion plot parameters ---
!
!     RF0,RFI0,RKX0,RKY0,RKZ0 : STANDARD PARAMETER FOR ROOT FINDING
!     RX0,RY0,RZ0             : STANDARD POSITION FOR ROOT FINDING
!     RK0,RKANG0              : STANDARD WAVE NUMBER AND ANGLE TO B (B:y,k:x,y)

    RF0    = 160.D3
    RFI0   =   0.D0
    RKX0   = 800.D0
    RKY0   = 160.D0
    RKZ0   =   0.D0
    RX0    = RR
    RY0    = 0.D0
    RZ0    = 0.D0
    RK0    = 100.D0
    RKANG0 = 89.70

!
!     *********** INPUT PARAMETERS fof DP program ***********
!
!     --- dispersin range parameters ---
!
!     RF1,RF2                 : SCAN RANGE OF REAL FREQUENCY (MHZ)
!     RFI1,RFI2               : SCAN RANGE OF IMAGINARY FREQUENCY (MHZ)
!     RKX1,RKX2               : SCAN RANGE OF WAVE NUMBER KX (1/M)
!     RKY1,RKY2               : SCAN RANGE OF WAVE NUMBER KY (1/M)
!     RKZ1,RKZ2               : SCAN RANGE OF WAVE NUMBER KZ (1/M)
!     RX1,RX2                 : SCAN RANGE OF POSITION X (M)
!     RY1,RY2                 : SCAN RANGE OF POSITION Y (M)
!     RZ1,RZ2                 : SCAN RANGE OF POSITION Z (M)
!     RK1,RK2                 : SCAN RANGE OF WAVE NUMBER K (1/M) 

    RF1    = 80000.D0
    RF2    = 16000.D0
    RFI1   =-5.D0
    RFI2   = 5.D0
    RKX1   = 0.D0
    RKX2   = 1600.D0
    RKY1   = 0.D0
    RKY2   = 1600.D0
    RKZ1   = 0.D0
    RKZ2   = 1600.D0
    RX1    = RR
    RX2    = RR+RA
    RY1    = -0.5D0
    RY2    =  0.5D0
    RZ1    = -RKAP*RA
    RZ2    =  RKAP*RA
    RK1    = 0.D0
    RK2    = 1600.D0

!     NGXMAX : NUMBER OF 1D SCAN POINTS
!     NGYMAX : NUMBER OF 2D SCAN POINTS
!     NGPMAX : NUMBER OF PARAMETER SCAN POINTS

    NGXMAX  = 21
    NGYMAX  = 21
    NGPMAX  = 21

!     EPSDP : Convergence torelance for root finding
!     EPSRF : Torelance for dumped mode (f_i < EPSRF*ABS(f_r) removed)

    EPSDP  = 1.D0
    EPSRF  = 3.D-3

!     --- Graphic parameters ---
!                 WC: absolute value of angular cyclotron frequency
!                 VT: thermal velocity
!                 VA: Alfven velocity
!     NORMF  : frequency normalization
!                 0: no normalization (Hz)
!                 positive: normalized by WC of NS=NORMF
!     NORMK  : wave number normalization
!                 0: no normalization (Hz)
!                 positive: normalized by VT/WC of NS=NORMK
!                 negative: normalized by VA/WC of NS=NORMK

    NORMF=0
    NORMK=0

!     NFLOUT: file output control
!             0: no file output
!            21: D4 and D5 point data

    NFLOUT=0

    RETURN
  END SUBROUTINE dp_init
END MODULE dpinit


