! dpinit.f90

MODULE dpinit

  PRIVATE
  PUBLIC dp_init

  CONTAINS

  SUBROUTINE dp_init

    USE dpcomm_parm
    IMPLICIT NONE
    INTEGER:: NS

!
!     *********** INPUT PARAMETERS ***********
!
!     --- PLASMA MODEL PARAMETERS ---
!
!     MODELP: TYPE OF ANALYTIC DIELECTRIC TENSOR
!                 0 : COLLISIONLESS COLD MODEL
!                 1 : COLLISIONAL COLD MODEL
!                 2 : IDEAL MHD MODEL
!                 3 : RESISTIVE MHD MODEL
!                 4 : KINETIC MODEL WITHOUT FLR
!                 5 : KINETIC MODEL WITH FLR (symmetric)
!                 6 : KINETIC MODEL WITH FLR
!                 7 : KINETIC MODEL WITH RELATIVISTIC EFFECTS (test)
!                 8 : DRIFTKINETIC MODEL (coming)
!                 9 : GYROKINETIC MODEL (coming)
!                11 : (WM) MHD plasma
!                12 : (WM) Cold plasma
!                13 : (WM) Hot plasma (No FLR)
!                14 : (WM) Hot plasma (Cold FLR)
!                15 : (WM) Hot plasma (FLR)
!                16 : (WM) Drift kinetic plasma
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
!     MODELV : NUMERICAL MODEL (*: not yet implemented)
!              0 : ANALYTIC MODEL
!              1 : KINETIC: ANALYTIC MAXWELLIAN DISTRIBUTION
!              2 : KINETIC: READ FPDATA DISTRIBUTION
!              3 : KINETIC: ANALYTIC MAXWELLIAN DISTRIBUTUION (RELATIVISTIC)
!              4 : KINETIC: READ FPDATA DISTRIBUTION (RELATIVISTIC)
!              5 : DRIFTKINETIC: ANALYTIC MAXWELLIAN DISTRIBUTION
!              6 : DRIFTKINETIC: READ FPDATA DISTRIBUTION
!              9 : LOCAL MODEL (MODELV locally specified by MODELVR)
!
!     NDISP1: MINIMUM HARMONIC NUMBER
!     NDISP2: MAXMUM  HARMONIC NUMBER

    DO NS=1,NSM
       MODELP(NS)= 0
       MODELV(NS)= 0
       NDISP1(NS)=-2
       NDISP2(NS)= 2
       PMAX(NS)= 7.D0
       NS_NSA(NS)=NS
    ENDDO

!     --- Root finding and dispersion plot parameters ---
!
!     RF0,RFI0,RKX0,RKY0,RKZ0 : STANDARD PARAMETER FOR ROOT FINDING
!     RX0,RY0,RZ0             : STANDARD POSITION FOR ROOT FINDING
!     RF1,RF2                 : SCAN RANGE OF REAL FREQUENCY (MHZ)
!     RFI1,RFI2               : SCAN RANGE OF IMAGINARY FREQUENCY (MHZ)
!     RKX1,RKX2               : SCAN RANGE OF WAVE NUMBER KX (1/M)
!     RKY1,RKY2               : SCAN RANGE OF WAVE NUMBER KY (1/M)
!     RKZ1,RKZ2               : SCAN RANGE OF WAVE NUMBER KZ (1/M)
!     RX1,RX2                 : SCAN RANGE OF POSITION X (M)
!     RY1,RY2                 : SCAN RANGE OF POSITION Y (M)
!     RZ1,RZ2                 : SCAN RANGE OF POSITION Z (M)
!
!     NGXMAX : NUMBER OF 1D SCAN POINTS
!     NGYMAX : NUMBER OF 2D SCAN POINTS
!     NGPMAX : NUMBER OF PARAMETER SCAN POINTS
!
!     EPSRT  : CONVERGENCE CRITERION OF ROOT FINDING
!     LMAXRT : MAXIMUM ITERATION COUNT OF ROOT FINDING

    RF0    = 160.D3
    RFI0   =   0.D0
    RKX0   = 800.D0
    RKY0   = 160.D0
    RKZ0   =   0.D0
    RX0    = RR
    RY0    = 0.D0
    RZ0    = 0.D0
!
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

    NGXMAX  = 21
    NGYMAX  = 21
    NGPMAX  = 21

    EPSRT  = 1.D-8
    LMAXRT = 20

!     --- Velocity distribution function parameters ---
!
!     NPMAX  : number of momentum magnitude mesh
!     NTHMAX : number of momentum angle mesh
!     NRMAX  : number of radial mesh
!     NSAMAX : number of test particle species
!     MODEFA : Type of fast particle contribution
!     RMIN   : minimum of radial mesh (r/a)
!     RMAX   : maximum of radial mesh (r/a)

    NPMAX=50
    NTHMAX=50
    NRMAX=3
    NSMAX=2
    NSAMAX=2
    RMIN=0.1D0
    RMAX=0.3D0
    MODEFA   = 0
    RMIN=0.D0
    RMAX=1.D0
!
    RETURN
  END SUBROUTINE dp_init
END MODULE dpinit
