!     $Id$
!***************************************************************
!
!   Set constants and initial parameters
!
!***************************************************************

SUBROUTINE TXINIT
  use tx_commons

  implicit none

  !   ***** Configuration parameters *****

  !   Plasma minor radius (m)
  RA = 0.35D0

  !   Wall radius (m)
  RB = 0.38D0

  !   Partition radius (m) available
  RC = RA

  !   Plasma major radius (m)
  RR = 1.3D0

  !   Toroidal magnetic field (T)
  BB = 1.3D0

  !   Plasma current start (MA)
  rIPs= 0.150D0

  !   Plasma current end (MA)
  rIPe= 0.150D0

  !   ***** Plasma components *****

  !   Atomic number of ion
  PA = 1.D0

  !   Charge number of ion
  PZ = 1.D0

  !   Effective charge
  Zeff = 2.D0

  !   ***** Initial plasma parameters *****

  !   Initial electron density at r = 0 (10^20 m^-3)
  PN0 = 0.4D0

  !   Initial electron density at r = a (10^20 m^-3)
  PNa = 0.05D0

  !   Electron density in diverter region (Minimum density in SOL)
  PNeDIV = 10.D-3

  !   Initial electron temperature at r = 0 (keV)
  PTe0 = 700.D-3

  !   Initial electron temperature at r = a (keV)
  PTea =  50.D-3

  !   Electron temperature in diverter region (Minimum Te in SOL)
  PTeDIV = 10.D-3

  !   Initial ion temperature  at r = 0 (keV)
  PTi0 = 700.D-3

  !   Initial ion temperature  at r = a (keV)
  PTia =  50.D-3

  !   Ion temperature in diverter region (Minimum Ti in SOL)
  PTiDIV = 10.D-3

  !   Initial current profile parameter
  PROFJ = 2.D0

  !   Initial density profile parameters
  PROFN1 = 3.D0
  PROFN2 = 1.D0

  !   Initial temperature profile parameters
  PROFT1 = 2.D0
  PROFT2 = 2.D0

  !   ***** Particle diffusivity and viscosity parameters *****

  !   Turbulent pinch velocity parameter
  VWpch0 = 0.D0

  !   Electron-driven diffusion parameter
  De0 = 0.1D0

  !   Ion-driven diffusion parameter
  Di0 = 0.D0

  !   Electron viscosity parameter
  rMue0 = 3.D0

  !   Ion viscosity parameter
  rMui0 = 3.D0

  !   Wave-particle interaction parameter
  WPE0 = 0.D0
  WPI0 = 1.D0

  !   Drift frequency parameter (omega/omega*e)
  WPM0 = 0.D0

  !  ***** Thermal diffusivity parameters *****

  !   Electron thermal diffusivity parameter (Chie/D)
  !     0 for fixed temperature profile
  Chie0 = 3.D0

  !   Ion thermal diffusivity parameter (Chie/D)
  !     0 for fixed temperature profile
  Chii0 = 3.D0

  !   ***** Turbulent transport control parameters *****

  !   Fixed transport coefficient parameter
  !   Suitable between 0.01 and 0.05
  FSDFIX = 0.05D0

  !   CDBM transport coefficient parameter
  !     (Current diffusive ballooning mode)
  !      The finer time step size, typically less than or equal to DT=5.D-4,
  !         is anticipated when using CDBM.
  FSCDBM = 0.D0

  !   Bohm transport coefficient parameter in SOL
  FSBOHM = 0.D0

  !   Pseud-classical particle transport coefficient parameter in SOL
  FSPCLD = 0.D0

  !   Pseud-classical heat & mom transport coefficient parameter in SOL
  FSPCLC = 0.D0

  !   Particle diffusion coefficient profile parameter (D(r=a)/D(r=0))
  PROFD =  3.D0

  !   Heat & mom diffusion coefficient profile parameter (D(r=a)/D(r=0))
  PROFC = 10.D0

  !   ***** Other transport parameters *****

  !   Charge exchange parameter
  FSCX = 1.D0

  !   Orbit loss parameter
  !     FSLC = 0 : No orbit loss included.
  !            1 : Loss term is expressed as damping term.
  !            2 : Loss term is expressed as source and sink term.
  FSLC = 0.D0

  !   Orbit loss model
  !     MDLC = 1 : Shaing model
  !            2 : Itoh model
  MDLC = 1

  !   Ripple loss parameter
  !     FSRP = 0.D0 : No ripple effect
  !            1.D0 : Full ripple effect
  !            2.D0 : Ripple effect but no ripple trapped ions
  FSRP = 0.D0

  !   Alpha heating paramter
  FSNF = 0.D0

  !   Neoclassical thermal diffusivity parameter
  ChiNC = 0.D0

  !   Poloidal neoclassical viscosity parameter
  FSNC = 1.D0

  !   Helical neoclassical viscosity parameter
  FSHL = 0.D0

  !   Particle loss to divertor parameter (default = 0.3)
  FSLP = 0.3D0

  !   Heat loss to divertor parameter (default = 1.0)
  FSLTE = 1.D0
  FSLTI = 1.D0

  !   Ionization parameter
  FSION = 1.D0

  !   Slow neutral diffusion factor
  FSD01 = 1.D0

  !   Fast neutral diffusion factor
  FSD02 = 1.D0

  !   Factor of E x B rotation shear
  rG1 = 24.D0

  !   ***** initial parameters *****

  !   Initial Density scale length in SOL (m), valid if MDITSN /= 0
  rLn = 0.03D0

  !   Initail Temperature scale length in SOL (m), valid if MDITST /= 0
  rLT = 0.030D0

  !   ***** Heating parameters *****

  !   NBI beam energy (keV)
  Eb = 32.D0

  !   Heating radius of perp NBI heating (m)
  RNBP  = 0.175D0

  !   Heating center of perp NBI heating (m)
  RNBP0 = 0.D0

  !   Heating radius of first tangential NBI heating (m)
  RNBT1  = 0.175D0

  !   Heating radius of second tangential NBI heating (m)
  RNBT2  = 0.175D0

  !   Heating center of first tangential NBI heating (m)
  RNBT10 = 0.D0

  !   Heating center of second tangential NBI heating (m)
  RNBT20 = 0.D0

  !   Perpendicular NBI input power (MW)
  PNBHP = 0.D0

  !   First tangential NBI input power (MW)
  PNBHT1 = 0.D0

  !   Second tangential NBI input power (MW)
  PNBHT2 = 0.D0

  !   NBI current drive parameter
  PNBCD = 1.D0

  !   Rate of the collisional slowing down part of the perpendicular NBI
  PNBMPD = 0.D0

  !   Different NBI deposition profiles for electrons and ions due to banana orbit effect
  !     MDLNBD = 0 : No charge separation
  !              1 : Orbit effect for banana particles only
  !              2 : Orbit effect for all beam ions including passing particles
  MDLNBD = 0

  !   Momentum input model for beam ions
  !     MDLMOM = 0 : Parallel velocity = initial beam velocity, Toroidal torque only
  !              1 : Toroidal and poloidal torques are given distinctly from "Vbpara".
  MDLMOM = 0 

  !   Refractive index of RF waves
  rNRFe = 0.D0
  rNRFi = 0.D0

  !   Heating width of RF heating (m)
  RRFew = 0.175D0
  RRFiw = 0.175D0

  !   Heating center radius of RF heating (m)
  RRFe0 = 0.D0
  RRFi0 = 0.D0

  !   RF input power (MW)
  PRFHe = 0.D0
  PRFHi = 0.D0

  !   Virtual torque input (N m)
  Tqi0  = 0.D0

  !   ***** Neutral parameters *****

  !   Initial Neutral density (10^20 m^-3)
  PN0s = 1.D-8

  !   Neutral thermal velocity (m/s)
  !     V0 = SQRT(2.D0*X[eV]*AEE/(PA*AMP))
  V0 = 1.5D3

  !   Recycling rate in SOL
  rGamm0 = 0.8D0

  !   Gas-puff particle flux (10^20 m^-2 1/s)
  !      If you input Gamma_0 [particles/sec], you translate it into
  !      rGASPF [1/(m^2 s)]:
  !        rGASPF = Gamma_0 / (2.D0*PI*RR*2.D0*PI*RB)
  rGASPF = 0.1D0

  !   ***** Ripple loss parameters *****

  ! Number of toroidal field coils in JT-60U
  NTCOIL = 18

  ! Minimum ripple amplitude in the system (R = Rmag0 = 2.4m on JT-60U)
  DltRPn = 0.0002D0 * 1.D-2

  ! Effective ellipticity for ripple amplitude estimation
  kappa = 1.2D0

  !   ***** Parameters for toroidal neoclassical viscosity *****

  ! Maximum poloidal mode number of error field (should be power of two)
  m_pol  = 32

  ! Toroidal mode number of error field
  n_tor  = NTCOIL

  ! Magnetic braiding parameters ***AF (2008-06-08)
  DMAG0   = 0.D0  ! Magnetic field line diffusivity : (Delta r)**2/Delta z [m]
  RMAGMN  = RA    ! minimum radius of the braiding region [m]
  RMAGMX  = RA    ! maximum radius of the braiding region [m]

  !   Helical ripple amplitude at r=a, linear in (r/a)
  EpsH = 0.1D0    ! amplitude 
  NCph = 5        ! toroidal pitch number
  NCth = 2        ! poloidal pitch number
  Q0 = 3.D0       ! q(0) by external coils
  QA = 2.D0       ! q(a) by external coils


  !   ***** Numerical parameters *****

  !   Time step size(s)
!!!!  DT = 1.D-4
  DT = 1.D-3

  !   Convergence parameter
  EPS = 1.D-2

  !   Iteration
  ICMAX = 10

  !   Time-advancing method
  !     ADV = 0     : Explicit scheme       (Not usable)
  !           0.5   : Crank-Nicolson scheme (Not usable)
  !           2/3   : Galerkin scheme       (Not usable)
  !           0.878 : Liniger scheme
  !           1     : Implicit scheme       (Recommended)
  ADV = 1.D0

  !   Mode of Backward Differential (Gear's) Formula (Second-order accuracy)
  !   (http://www.scholarpedia.org/article/Backward_differentiation_formulas)
  !   (Bibun Houteishiki no Suuchi Kaihou I, Taketomo Mitsui, Iwanami, 1993, p29,30)
  !   (P.M.Gresho and R.L.Sani, "Incompressible Flow and the Finite Element Method,
  !    Vol.1 Advection-Diffusion", John Wiley and Sons, 2000, p.263)
  !     IGBDF = 0   : not used
  !             1   : Use BDF2
  IGBDF = 0

  !   Lower bound of dependent variables
  tiny_cap = 1.d-14

  !   Amplitude for numerical convergence
  AMPe4 = 1.D3
  AMPm5 = 1.D3

  !   Permittivity switch for numerical convergence
  rMUb1 = rMU0
  rMUb2 = 1.d0

  !   ***** Mesh number parameters *****

  !   Magnitude of mesh peakness
  CMESH0 =  2.D0
  CMESH  = 30.D0

  !   Width of mesh peakness
  WMESH0 = 0.2D0
  WMESH  = 5.D-2

  !   Number of nodes
  NRMAX  = 50

  !   Number of elements
  NEMAX = NRMAX

  !   ***** Time parameter *****

  !   Number of time step
  NTMAX = 10

  !   ***** Graphic parameters *****

  !   Time step interval between print output
  NTSTEP = 10

  !   Time step interval between lines in f(r) graph
  NGRSTP = 1

  !   Time step interval between points in f(t) graph
  NGTSTP = 1

  !   Time step interval between points in f(t) graph
  NGVSTP = 1

  !   Mode of Graph
  !   0 : for Display (with grid, w/o power)
  !   1 : for Display (with grid and power)
  !   2 : for Print Out (w/o grid, with power)
  MODEG = 1

  !   MODE of Graph Line
  !   0 : Change Line Color (Last Color Fixed)
  !   1 : Change Line Color and Style
  !   2 : Change Line Color and Style (With Legend)
  !   3 : Change Line Color, Style and Mark
  !   4 : Change Line Color, Style and Mark (With Legend)
  MODEGL = 1

  !   Mode of AV
  !   0 : OFF
  !   n : Number of Display
  MODEAV = 0

  !   Diagnostic parameter
  !   0 : OFF
  !   1 : debug message output (ntstep == 1)
  !   2 : debug message output (few)
  !   3 : debug message output (many)
  !  -1 : message for steady state check at NGTSTP intervals
  IDIAG = 0

  !   ***** Model parameters *****

  !   Mode of LAPACK
  !   0    : Use BANDRD
  !   else : Use LAPACK_DGBSV
  MDLPCK = 0

  !   Mode of Wave-particle interaction model
  !   0    : Non-ambipolar model
  !   1    : Shaing model
  MDLWTB = 0

  !   Mode of neoclassical resistivity model
  !   0    : original
  !   1    : NCLASS
  !   2    : Sauter
  !   3    : Hirshman, Hawryluk and Birge
  MDLETA = 0

  !   Mode of fixed temperature profile
  !   0    : not fixed
  !   1    : fixed
  MDFIXT = 0

  !   Mode of nonlinear iteration method
  !   0    : fixed point method (or Picard method), linear convergence
  !   1    : Secant method, golden ratio (1.618...) convergence (NOT working now)
  !   2    : Steffensen's method, quadratic convergence (Validity unconfirmed)
  MDSOLV = 0

  !   Mode of initial density and temperature profiles in the SOL
  !   0    : polynominal model
  !   1    : exponential decay model
  MDITSN = 0  ! for density profiles
  MDITST = 0  ! for temperature profiles

  !   Mode of initial temperature profiles
  !   0    : original
  !   1    : pedestal model
  MDINTT = 0

  !   Mode of initial plasma profils
  !   0    : many profiles are analytically calculated 
  !   1    : minimal profiles are calculated
  MDINIT = 0

  !   Multiplication factor for graphic in the radial direction
  !   default : 1.0
  !
  gDIV(1:NGYRM) = 1.0
  gDIV(1)  = 1.E20
  gDIV(2)  = 1.E14
  gDIV(4)  = 1.E3
  gDIV(5)  = 1.E3
  gDIV(7)  = 1.E3
  gDIV(8)  = 1.E3
  gDIV(9)  = 1.E3
  gDIV(12) = 1.E18
  gDIV(13) = 1.E3
  gDIV(16) = 1.E14
  gDIV(17) = 1.E-3
  gDIV(19) = 1.E3
  gDIV(21) = 1.E6
  gDIV(22) = 1.E6
  gDIV(23) = 1.E3
  gDIV(24) = 1.E3
  gDIV(25) = 1.E3
  gDIV(26) = 1.E3
  gDIV(27) = 1.E3
  gDIV(28) = 1.E3
  gDIV(35) = 1.E14
  gDIV(36) = 1.E12
  gDIV(37) = 1.E20
  gDIV(38) = 1.E3
  gDIV(41) = 1.E6
  gDIV(42) = 1.E3
  gDIV(45) = 1.E3
  gDIV(46) = 1.E3
  gDIV(50) = 1.E3
  gDIV(51) = 1.E3
  gDIV(53) = 1.E-2
  gDIV(55) = 1.E3
  gDIV(56) = 1.E-4
  gDIV(57) = 1.E-8
  gDIV(60) = 1.E6
  gDIV(64) = 1.E-20
  gDIV(65) = 1.E-20
  gDIV(72) = 1.E3
  gDIV(73) = 1.E6
  gDIV(74) = 1.E6
  gDIV(75) = 1.E6
  gDIV(89) = 1.E3
  gDIV(90) = 1.E20
  gDIV(96) = 1.E6
  gDIV(97) = 1.E6
  gDIV(98) = 1.E6
  gDIV(100) = 1.E6
  gDIV(101) = 1.E6
  gDIV(102) = 1.E6
  gDIV(103) = 1.E6
  gDIV(104) = 1.E6
  gDIV(109) = 1.E15
  gDIV(110) = 1.E-4
!  gDIV(115) = 1.E15
  gDIV(121) = 1.E6
  gDIV(123) = 1.E6
  gDIV(124) = 1.E6

  !   *** Obsolete or not used parameter ***

  !   Radius where density increase by command DEL
  DelR = 0.175D0

  !   Amount of increase of density by command DEL
  DelN = 5.D-1

  !   Index for graphic save interval

  NGR=-1

  RETURN
END SUBROUTINE TXINIT

!***************************************************************
!
!        Calculate mesh, etc.
!
!***************************************************************

SUBROUTINE TXCALM

  use tx_commons
  use tx_interface, only : LORENTZ, LORENTZ_PART, BISECTION

  implicit none
  INTEGER(4) :: NR, NRL, NR_RC_NEAR
  real(8)    :: DR, MAXAMP, CL, WL, C1L, C2L, W1L, W2L, RL, RCL, CLNEW

  !   Ion mass number
  AMI   = PA * AMP
  !   Beam ion mass number
  AMB   = AMI
  !   Radial step width
  !    DR    = RB / NRMAX
  !   Number of equations
  NQMAX = NQM

  !   Helical system
  UHth  = DBLE(NCth) / DBLE(NCph)
  UHph  = 1.D0

  !   Square root permittivity for LQm1
  sqeps0 = sqrt(EPS0)

  !  Mesh

  !  As a trial, generate mesh using given CL and WL and seek the position in the
  !  original coordinate, which becomes the nearest mesh of separatrix after mapping.
  C1L = CMESH0
  C2L = CMESH
  W1L = WMESH0
  W2L = WMESH
  MAXAMP = LORENTZ(RB,C1L,C2L,W1L,W2L,0.D0,RC) / RB
  NR_RC_NEAR = 0
  R(0) = 0.D0
  DO NR = 1, NRMAX - 1
     RL = DBLE(NR) / DBLE(NRMAX) * RB
     CALL BISECTION(LORENTZ,C1L,C2L,W1L,W2L,0.D0,RC,MAXAMP,RL,RB,R(NR))
     IF(ABS(R(NR)-RC) <= ABS(R(NR)-R(NR_RC_NEAR))) NR_RC_NEAR = NR
  END DO
  R(NRMAX) = RB

  !  Construct new CL value that separatrix is just on mesh.
  !  New CL is chosen in order not to be settle so far from given CL.
  !  The mesh finally obtained is well-defined.
  RCL = DBLE(NR_RC_NEAR) / DBLE(NRMAX) * RB
  CLNEW = ( (RC - RCL) * RB - RCL * C1L * LORENTZ_PART(RB,W1L,W2L,0.D0,RC,0) &
       &   + RB * C1L * LORENTZ_PART(RC,W1L,W2L,0.D0,RC,0)) &
       &  / (  RCL * LORENTZ_PART(RB,W1L,W2L,0.D0,RC,1) &
       &     - RB  * LORENTZ_PART(RC,W1L,W2L,0.D0,RC,1))
  MAXAMP = LORENTZ(RB,C1L,CLNEW,W1L,W2L,0.D0,RC) / RB
  R(0) = 0.D0
  DO NR = 1, NRMAX - 1
     RL = DBLE(NR) / DBLE(NRMAX) * RB
     CALL BISECTION(LORENTZ,C1L,CLNEW,W1L,W2L,0.D0,RC,MAXAMP,RL,RB,R(NR))
  END DO
  R(NRMAX) = RB

  !  Maximum NR till RA

  DO NR = 0, NRMAX-1
     IF(R(NR) <= RA.AND.R(NR+1) >= RA) THEN
        NRL = NR
        EXIT
     END IF
  END DO

  !  Adjust RC on mesh

  IF(ABS(R(NRL)-RC) < ABS(R(NRL+1)-RC)) THEN
     NRA = NRL
     RL = 0.5d0 * ABS(R(NRA) - RC)
     R(NRA  ) = RC
     R(NRA-1) = R(NRA-1) + RL
  ELSE
     NRA = NRL + 1
     RL = 0.5d0 * ABS(R(NRA) - RC)
     R(NRA  ) = RC
     R(NRA+1) = R(NRA+1) - RL
  END IF

  !  Mesh number at the center of the core plasma

  DO NR = 0, NRMAX-1
     IF(R(NR) <= 0.5d0*RA.AND.R(NR+1) >= 0.5d0*RA) THEN
        IF(ABS(R(NR)-0.5d0*RA) < ABS(R(NR+1)-0.5d0*RA)) THEN
           NRC = NR
        ELSE
           NRC = NR + 1
        END IF
        EXIT
     END IF
  END DO

  !  Mesh coordinate

  RHO(0:NRMAX) = R(0:NRMAX) / RA
  PSI(0:NRMAX) = R(0:NRMAX)**2

  !  Mesh interval

  H   (1:NEMAX) = R  (1:NRMAX) - R  (0:NRMAX-1)
  HPSI(1:NEMAX) = PSI(1:NRMAX) - PSI(0:NRMAX-1)

  RETURN
END SUBROUTINE TXCALM

!***************************************************************
!
!   Initialize profiles
!
!***************************************************************

SUBROUTINE TXPROF

  use tx_commons
  use tx_variables
  use tx_interface, only : INTG_P, INTDERIV3, detect_datatype, INTG_F, dfdx

  implicit none
  INTEGER(4) :: NR, NQ, I, IER, ifile
  REAL(8) :: RL, PROF, PROFN, PROFT, PTePROF, PTiPROF, QL, RIP1, RIP2, dRIP, SSN, SSPe, SSPi
  REAL(8) :: AJFCT, SUM_INT
  REAL(8) :: ALP, dPe, dPi, DR1, DR2
  REAL(8) :: EpsL, Vte, Wte, rNuAsE_inv, FTL, EFT, CR
  real(8) :: FACT, PBA, dPN, CfN1, CfN2, pea, pia, pediv, pidiv, dpea, dpia, &
       &     Cfpe1, Cfpe2, Cfpi1, Cfpi2
  REAL(8) :: DERIV4, FCTR ! External functions
  real(8), dimension(:), allocatable :: AJPHL, TMP, RHSV, dPedr, dPidr
  real(8), dimension(:,:), allocatable :: CMTX

  !  Read spline table for neoclassical toroidal viscosity
  IF(FSRP /= 0.D0) CALL Wnm_spline(fmnq, wnm, umnq, nmnqm)

  NEMAX = NRMAX

  !  Define basic quantities like mass of particles, mesh, etc.

  CALL TXCALM

  !  Contribution of perturbed magnetic field
!  IF(DltRPn /= 0.D0) CALL perturb_mag

  !  Initialize variable vector

  X(1:NQMAX,0:NRMAX) = 0.D0

  !  Variables

  allocate(AJPHL(0:NRMAX))
  PBA   = RB - RA
  dPN   = - 3.D0 * (PN0 - PNa) / RA
  CfN1  = - (3.D0 * PBA * dPN + 4.D0 * (PNa - PNeDIV)) / PBA**3
  CfN2  =   (2.D0 * PBA * dPN + 3.D0 * (PNa - PNeDIV)) / PBA**4
  IF(MDFIXT == 0) THEN
     pea   = PNa    * PTea   ;  pia  = PNa    / PZ * PTia
     dpea  = dPN    * PTea   ; dpia  = dPN    / PZ * PTia
     pediv = PNeDIV * PTeDIV ; pidiv = PNeDIV / PZ * PTiDIV
  ELSE
     pea   = PTea   ;  pia   = PTia
     dpea  = 0.d0   ; dpia   = 0.d0
     pediv = PTeDIV ; pidiv = PTiDIV
  END IF
  Cfpe1 = - (3.D0 * PBA * dpea + 4.D0 * (pea - pediv)) / PBA**3
  Cfpe2 =   (2.D0 * PBA * dpea + 3.D0 * (pea - pediv)) / PBA**4
  Cfpi1 = - (3.D0 * PBA * dpia + 4.D0 * (pia - pidiv)) / PBA**3
  Cfpi2 =   (2.D0 * PBA * dpia + 3.D0 * (pia - pidiv)) / PBA**4
  DO NR = 0, NRMAX
     RL=R(NR)
     IF (RL < RA) THEN
        PROFN = (1.D0 - RHO(NR)**PROFN1)**PROFN2
        PROFT = (1.D0 - RHO(NR)**PROFT1)**PROFT2
        X(LQe1,NR) = (PN0 - PNa) * PROFN + PNa ! Ne
        X(LQi1,NR) = X(LQe1,NR) / PZ           ! Ni
        IF(MDINTT == 0) THEN
           PTePROF = (PTe0 - PTea) * PROFT + PTea
           PTiPROF = (PTi0 - PTia) * PROFT + PTia
           IF(MDFIXT == 0) THEN
              X(LQe5,NR) = PTePROF * X(LQe1,NR) ! Ne*Te
              X(LQi5,NR) = PTiPROF * X(LQi1,NR) ! Ni*Ti
           ELSE 
              X(LQe5,NR) = PTePROF ! Te
              X(LQi5,NR) = PTiPROF ! Ti
           END IF
        ELSE
           PTePROF = (PTe0 - PTea) * (0.8263d0 * (1.d0 - RHO(NR)**2 )**1.5d0 &
                &                   + 0.167d0  * (1.d0 - RHO(NR)**30)**1.25d0) + PTea
           PTiPROF = (PTi0 - PTia) * (0.8263d0 * (1.d0 - RHO(NR)**2 )**1.5d0 &
                &                   + 0.167d0  * (1.d0 - RHO(NR)**30)**1.25d0) + PTia
           IF(MDFIXT == 0) THEN
              X(LQe5,NR) = PTePROF * X(LQe1,NR) ! Ne*Te
              X(LQi5,NR) = PTiPROF * X(LQi1,NR) ! Ne*Te
           ELSE
              X(LQe5,NR) = PTePROF ! Te
              X(LQi5,NR) = PTiPROF ! Ti
           END IF
        END IF
     ELSE
        IF(MDITSN == 0) THEN
           X(LQe1,NR) = PNa + dPN * (RL - RA) + CfN1 * (RL - RA)**3 &
                &                             + CfN2 * (RL - RA)**4
           X(LQi1,NR) = X(LQe1,NR) / PZ
        ELSE
           X(LQe1,NR) = PNa * EXP(- (RL - RA) / rLN)
           X(LQi1,NR) = X(LQe1,NR) / PZ
        END IF
        IF(MDITST == 0) THEN
           X(LQe5,NR) = pea + dpea * (RL - RA) + Cfpe1 * (RL - RA)**3 &
                &                              + Cfpe2 * (RL - RA)**4
           X(LQi5,NR) = pia + dpia * (RL - RA) + Cfpi1 * (RL - RA)**3 &
                &                              + Cfpi2 * (RL - RA)**4
        ELSE
           PTePROF = PTea * EXP(- (RL - RA) / rLT)
           PTiPROF = PTia * EXP(- (RL - RA) / rLT)
           IF(MDFIXT == 0) THEN
              X(LQe5,NR) = PTePROF * X(LQe1,NR)
              X(LQi5,NR) = PTiPROF * X(LQi1,NR)
           ELSE
              X(LQe5,NR) = PTePROF
              X(LQi5,NR) = PTiPROF
           END IF
        END IF
     END IF
     ! N0_1 (slow neutrals)
     X(LQn1,NR) = PN0s
     ! N0_2 (fast neutrals)
     X(LQn2,NR) = 0.D0
     ! Bphi
     X(LQm5,NR) = 0.5D0 * PSI(NR) * BB / rMU0 / AMPm5
     BphV(NR)   = BB
     ! Fixed densities to keep them constant during iterations
     PNeV_FIX(NR) = X(LQe1,NR)
     PNiV_FIX(NR) = X(LQi1,NR)
     IF(MDFIXT == 0) THEN
        PTeV_FIX(NR) = X(LQe5,NR) / X(LQe1,NR)
        PTiV_FIX(NR) = X(LQi5,NR) / X(LQi1,NR)
     ELSE
        PTeV_FIX(NR) = X(LQe5,NR)
        PTiV_FIX(NR) = X(LQi5,NR)
     END IF
  END DO

  ! Poloidal magnetic field

  IF(FSHL == 0.D0) THEN
     BthV(0) = 0.D0
     DO NR = 1, NRMAX
        RL = R(NR)
        IF(RL < RA) THEN
           PROF = 1.D0 - (RL / RA)**2
           ! Btheta
           BthV(NR) = rMUb1 * rIPs * 1.D6 / (2.D0 * PI * RL) * (1.D0 - PROF**(PROFJ+1))
        ELSE
           BthV(NR) = rMUb1 * rIPs * 1.D6 / (2.D0 * PI * RL)
        END IF
     END DO
  ELSE
     BthV(0) = 0.D0
     DO NR = 1, NRMAX
        RL = R(NR)
        QL = (Q0 - QA) * (1.D0 - (RL / RA)**2) + QA
        BthV(NR) = BB * RL / (QL * RR)
     END DO
  END IF
  Bthb = BthV(NRMAX)

  ! Toroidal electron current

  ifile = detect_datatype('LQe4')
  if(ifile == 0) then
     DO NR = 0, NRMAX
!!$        AJPHL(NR) = 2.D0 / rMUv1 * DERIV4(NR,PSI,R(0:NRMAX)*BthV(0:NRMAX),NRMAX,0)
!!$        X(LQe4,NR) = - AJPHL(NR) / (AEE * 1.D20) / AMPe4

        IF((1.D0-(R(NR)/RA)**2) <= 0.D0) THEN
           PROF= 0.D0    
        ELSE             
           PROF= 1.D0-(R(NR)/RA)**2
        END IF
        IF(FSHL == 0.D0) THEN
           ! Ne*UePhi
           AJPHL(NR) = rIPs * 1.D6 / (PI * RA**2) * (PROFJ + 1.D0) * PROF**PROFJ
           X(LQe4,NR) = - AJPHL(NR) / (AEE * 1.D20) / AMPe4
           AJOH(NR)= AJPHL(NR)
!           AJOH(NR)= PROF
        ELSE
           AJPHL(NR) = 0.D0
           X(LQe4,NR) = 0.D0
           AJOH(NR)= 0.D0
        END IF
     END DO
  else
     call inexpolate(infiles(ifile)%nol,infiles(ifile)%r,infiles(ifile)%data,NRMAX,RHO,5,AJPHL)
     AJFCT = rIPs * 1.D6 / (2.D0 * PI * INTG_F(AJPHL))
     DO NR = 0, NRMAX
        AJPHL(NR) = AJFCT * AJPHL(NR)
        X(LQe4,NR) = - AJPHL(NR) / (AEE * 1.D20) / AMPe4
        AJOH(NR) = AJPHL(NR)
     END DO
     BthV(0) = 0.d0
     SUM_INT = 0.d0
     DO NR = 1, NRMAX
        SUM_INT = SUM_INT + INTG_P(AJPHL,NR,0)
        BthV(NR) = rMUb1 * SUM_INT / R(NR)
     END DO
  end if

  ! Inverse matrix of derivative formula for integration

  allocate(CMTX(1:NRMAX,1:NRMAX),RHSV(1:NRMAX))
  CMTX(:,:) = 0.D0
  DO NR = 1, NRMAX
     IF(NR == 1) THEN
        DR1 = PSI(NR-1) - PSI(NR)
        DR2 = PSI(NR+1) - PSI(NR)
        CMTX(NR,NR  ) = - (DR1 + DR2) / (DR1 * DR2)
        CMTX(NR,NR+1) = - DR1**2 / (DR1 * DR2 * (DR2 - DR1))
     ELSEIF(NR == NRMAX) THEN
        DR1 = PSI(NR-1) - PSI(NR)
        DR2 = PSI(NR-2) - PSI(NR)
        CMTX(NR,NR-2) = - DR1**2 / (DR1 * DR2 * (DR2 - DR1))
        CMTX(NR,NR-1) =   DR2**2 / (DR1 * DR2 * (DR2 - DR1))
        CMTX(NR,NR  ) = - (DR1 + DR2) / (DR1 * DR2)
     ELSE
        DR1 = PSI(NR-1) - PSI(NR)
        DR2 = PSI(NR+1) - PSI(NR)
        CMTX(NR,NR-1) =   DR2**2 / (DR1 * DR2 * (DR2 - DR1))
        CMTX(NR,NR  ) = - (DR1 + DR2) / (DR1 * DR2)
        CMTX(NR,NR+1) = - DR1**2 / (DR1 * DR2 * (DR2 - DR1))
     END IF
  END DO
  CALL INVMRD(CMTX,NRMAX,NRMAX,IER)

  ! Numerical solution for AphV

  IF((FSHL == 0.0) .AND. (ifile == 0) .AND. &
     (PROFJ == 1 .OR. PROFJ == 2 .OR. PROFJ == 3 .OR. &
      PROFJ == 4 .OR. PROFJ == 5)) THEN

  ! Analytic solution for AphV (Valid for PROFJ=1,2,3,4,5, otherwise use below.)

     DO NR = 0, NRMAX
        IF(R(NR) < RA) THEN
           RL = R(NR) / RA
           IF(PROFJ == 1) THEN
              FACT = RL**2*(2.d0-0.5d0*RL**2)
           ELSE IF(PROFJ == 2) THEN
              FACT = RL**2*(3.d0-1.5D0*RL**2+RL**4/3.D0)
           ELSE IF(PROFJ == 3) THEN
              FACT = RL**2*(4.d0-3.d0*RL**2+(4.d0/3.d0)*RL**4-0.25d0*RL**6)
           ELSE IF(PROFJ == 4) THEN
              FACT = RL**2*(5.d0-5.d0*RL**2+(10.d0/3.d0)*RL**4-1.25d0*RL**6+0.2d0*RL**8)
           ELSE IF(PROFJ == 5) THEN
              FACT = RL**2*(6.d0-7.5d0*RL**2+(20.d0/3.d0)*RL**4-3.75d0*RL**6+1.2d0*RL**8 &
                   &-RL**10/6.d0)
           END IF
           X(LQm4,NR) = - rMUb1 * rIPs * 1.D6 / (4.D0 * PI) * FACT
        ELSE
           RL = 1.D0
           IF(PROFJ == 1) THEN
              FACT = RL**2*(2.d0-0.5d0*RL**2)
           ELSE IF(PROFJ == 2) THEN
              FACT = RL**2*(3.d0-1.5D0*RL**2+RL**4/3.D0)
           ELSE IF(PROFJ == 3) THEN
              FACT = RL**2*(4.d0-3.d0*RL**2+(4.d0/3.d0)*RL**4-0.25d0*RL**6)
           ELSE IF(PROFJ == 4) THEN
              FACT = RL**2*(5.d0-5.d0*RL**2+(10.d0/3.d0)*RL**4-1.25d0*RL**6+0.2d0*RL**8)
           ELSE IF(PROFJ == 5) THEN
              FACT = RL**2*(6.d0-7.5d0*RL**2+(20.d0/3.d0)*RL**4-3.75d0*RL**6+1.2d0*RL**8 &
                   &-RL**10/6.d0)
           END IF
           X(LQm4,NR) = - rMUb1 * rIPs * 1.D6 / (4.D0 * PI) * FACT &
                &       - rMUb1 * rIPs * 1.D6 / (4.D0 * PI) * LOG(PSI(NR)/RA**2)
        END IF
     END DO
  ELSE
     RHSV(1:NRMAX) = - 0.5D0 * BthV(1:NRMAX) / R(1:NRMAX)
     X(LQm4,0) = 0.D0
     X(LQm4,1:NRMAX) = matmul(CMTX,RHSV)
  END IF

  ! Poloidal current density (Virtual current for helical system)

  allocate(TMP(0:NRMAX))
  IF(FSHL == 0.D0) THEN
     AJV(0:NRMAX)=0.D0
  ELSE
     TMP(0:NRMAX) = R(0:NRMAX) * BthV(0:NRMAX)
     DO NR = 1, NRMAX
        dRIP = DERIV4(NR,R,TMP,NRMAX,0) * 2.D0 * PI / rMUb1
        AJV(NR)=dRIP / (2.D0 * PI * R(NR))
     END DO
     AJV(0)=FCTR(R(1),R(2),AJV(1),AJV(2))
!     write(6,'(I5,1P2E12.4)') (NR,R(NR),AJV(NR),NR=0,NRMAX)
  END IF
  deallocate(TMP)

  ! Toroidal electric field for initial NCLASS calculation

  Q(1:NRMAX) = ABS(R(1:NRMAX) * BB / (RR * BthV(1:NRMAX)))
  Q(0) = FCTR(R(1),R(2),Q(1),Q(2))

  DO NR = 0, NRMAX
     ! +++ Hirshman, Hawryluk and Birge model +++
     PNeV(NR) = X(LQe1,NR)
     PNiV(NR) = X(LQi1,NR)
     IF(MDFIXT == 0) THEN
        PTeV(NR) = X(LQe5,NR) / X(LQe1,NR)
     ELSE
        PTeV(NR) = X(LQe5,NR)
     END IF
     ! Inverse aspect ratio
     EpsL = R(NR) / RR
     Vte = SQRT(2.D0 * ABS(PTeV(NR)) * rKeV / AME) ! Thermal velocity for ions
     Wte = Vte / (Q(NR) * RR) ! Omega_te; transit frequency for electrons
     rlnLei(NR) = 37.8d0 - LOG(SQRT(PNeV(NR)*1.D20)/(PTeV(NR)))
     rNuei(NR) = PNiV(NR) * 1.D20 * PZ**2 * AEE**4 * rlnLei(NR) &
          &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AME) &
          &     * (ABS(PTeV(NR)) * rKeV)**1.5D0)
     rNuAsE_inv = EpsL**1.5D0 * Wte / (SQRT(2.D0) * rNuei(NR))
     ! Trapped particle fraction
     FTL  = 1.46D0 * SQRT(EpsL) - 0.46D0 * EpsL**1.5D0
     EFT  = FTL * rNuAsE_inv / (rNuAsE_inv + (0.58D0 + 0.20D0 * Zeff))
     ! Spitzer resistivity for hydrogen plasma
     ETAS(NR) = CORR(1.D0) * AME * rNuei(NR) / (PNeV(NR) * 1.D20 * AEE**2)
     CR   = 0.56D0 * (3.D0 - Zeff) / ((3.D0 + Zeff) * Zeff)
     IF(ABS(FSNC) > 0.D0) THEN
        ETA(NR) = ETAS(NR) * Zeff * (1.D0 + 0.27D0 * (Zeff - 1.D0)) &
             &  /((1.D0 - EFT) * (1.D0 - CR * EFT) * (1.D0 + 0.47D0 * (Zeff - 1.D0)))
     ELSE
        ETA(NR) = ETAS(NR)
     END IF
     IF(FSHL == 0.D0) THEN
        ! Ephi
        X(LQm3,NR) = - ETA(NR) *  AJPHL(NR)
     ELSE
        X(LQm3,NR) = 0.D0
     END IF
     IF(X(LQm3,NR) == 0.D0) X(LQm3,NR) = -1.D-4
  END DO
  deallocate(AJPHL)

  !   NBI total input power (MW)
  PNBH = PNBHP + PNBHT1 + PNBHT2

  T_TX=0.D0
  NGT=-1
  NGR=-1
  NGVV=-1
  rIP=rIPs

  !  Define physical variables from X

  CALL TXCALV(X)

  !  Calculate various physical quantities

  CALL TXCALC

  !  Initial condition Part II

  IF(MDINIT == 0) THEN

     allocate(dPedr(0:NRMAX),dPidr(0:NRMAX))
     dPedr(0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PeV,NRMAX,0) * rKeV
     dPidr(0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PiV,NRMAX,0) * rKeV

     DO NR = 0, NRMAX
        dPe = dPedr(NR)
        dPi = dPidr(NR)
        IF(NR == NRMAX) THEN
           dPe = 0.D0
           dPi = 0.D0
        END IF
        IF(rNueNC(NR) == 0.D0) THEN
           ALP = 0.D0
        ELSE
           ALP = (AMI / AME) * (rNuiNC(NR) / rNueNC(NR))
        END IF
        X(LQi3,NR) = (- BthV(NR) / BphV(NR) * X(LQe4,NR) * AMPe4 + (dPe + dPi) &
             &     / (AEE * BphV(NR))) / (PZ + ALP) * R(NR)
        X(LQe3,NR) =- ALP * X(LQi3,NR)
        ErV(NR)    =- BphV(NR) / PNiV(NR) * (- BthV(NR) / BphV(NR) * X(LQe4,NR) * AMPe4 &
             &      + (dPe + dPi) / (AEE * BphV(NR))) / (PZ + ALP) &
             &      + dPi / (PZ * AEE * PNiV(NR))
        X(LQe2,NR) =- (AMI * rNuiNC(NR) /(AEE * BphV(NR))) * X(LQi3,NR)
        X(LQi2,NR) = X(LQe2,NR) / PZ
        X(LQm3,NR) = BthV(NR) / PNeV(NR) * (-(AMI * rNuiNC(NR) /(AEE * BphV(NR))) &
             &     / (PZ + ALP) * (- BthV(NR) / BphV(NR) * X(LQe4,NR) * AMPe4 &
             &     +(dPe + dPi) / (AEE * BphV(NR)))) &
             &     + AME * rNuei3(NR) / (AEE * PNeV(NR)) * X(LQe4,NR) * AMPe4
     END DO

     deallocate(dPedr,dPidr)

     ! Scalar potential

     allocate(TMP(0:NRMAX))
     ! TMP(0) is an arbitrary value (INTDERIV3 does not require the value at axis node
     !                               in case of (LAST ARGUMENT)=1.)
     TMP(1:NRMAX) = - 0.5D0 * ErV(1:NRMAX) / R(1:NRMAX)
     TMP(0) = FCTR(PSI(1),PSI(2),TMP(1),TMP(2))
     CALL INTDERIV3(TMP,PSI,X(LQm1,0:NRMAX),0.D0,NRMAX,1)
     
     ! AthV

     TMP(1:NRMAX) = 0.5D0 * rMU0 * AEE * (X(LQe3,1:NRMAX) - PZ * X(LQi3,1:NRMAX)) * 1.D20 &
          &       / PSI(1:NRMAX)
     TMP(0) = FCTR(PSI(1),PSI(2),TMP(1),TMP(2))
     CALL INTDERIV3(TMP,PSI,BphV,BB,NRMAX,1)
     RHSV(1:NRMAX) = 0.5D0 * BphV(1:NRMAX)
     X(LQm5,0) = 0.D0
     X(LQm5,1:NRMAX) = matmul(CMTX,RHSV) / rMU0 / AMPm5
     deallocate(CMTX,RHSV,TMP)

     CALL TXCALV(X)
     CALL TXCALC

  END IF

  !  Calculate global quantities for storing and showing initial status

  CALL TXGLOB

  RETURN
END SUBROUTINE TXPROF

!*****************************************************************************************

module tx_parameter_control
  use tx_commons
  implicit none
  public
  NAMELIST /TX/ &
       & RA,RB,RC,RR,BB, &
       & PA,PZ,Zeff, &
       & PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ,PROFN1,PROFN2,PROFT1,PROFT2, &
       & De0,Di0,VWpch0,rMue0,rMui0,WPM0,WPE0,WPI0, &
       & Chie0,Chii0,ChiNC, &
       & FSDFIX,FSCDBM,FSBOHM,FSPCLD,FSPCLC,PROFD,PROFC, &
       & FSCX,FSLC,FSRP,FSNF,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD01,FSD02,MDLC, &
       & rLn,rLT, &
       & Eb,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBHP,PNBHT1,PNBHT2,PNBCD,PNBMPD, &
       & rNRFe,RRFew,RRFe0,PRFHe,Tqi0,rNRFi,RRFiw,RRFi0,PRFHi, &
       & PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV, &
       & NTCOIL,DltRPn,kappa,m_pol,n_tor, &
       & DT,EPS,ICMAX,ADV,tiny_cap,CMESH0,CMESH,WMESH0,WMESH, &
       & NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP, &
       & DelR,DelN, &
       & DMAG0,RMAGMN,RMAGMX,EpsH,NCph,NCth,&
       & rG1,FSHL,Q0,QA, &
       & rIPs,rIPe, &
       & MODEG,gDIV,MODEAV,MODEGL,MDLPCK,MDLWTB, &
       & MDLETA,MDFIXT,MDITSN,MDITST,MDINTT,MDINIT,IDIAG,IGBDF,MDSOLV,MDLNBD,MDLMOM
  private :: TXPLST

contains
!***************************************************************
!
!   Change input parameters
!
!***************************************************************

  SUBROUTINE TXPARM(KID)

    INTEGER(4) :: IST
    LOGICAL :: LEX
    character(len=*)  :: KID

    DO 
       WRITE(6,*) '# INPUT &TX :'
       READ(5,TX,IOSTAT=IST)
       IF(IST > 0) THEN
          CALL TXPLST
          CYCLE
       ELSE IF(IST < 0) THEN
          KID='Q'
          EXIT
       ELSE
          KID=' '
          EXIT
       END IF
    END DO
    IERR=0

  END SUBROUTINE TXPARM

  SUBROUTINE TXPARL(KLINE)

    integer(4) :: IST
    logical :: LEX
    character(len=*) :: KLINE
    character(len=90) :: KNAME

    KNAME=' &TX '//KLINE//' &END'
    WRITE(7,'(A90)') KNAME
    REWIND(7)
    READ(7,TX,IOSTAT=IST)
    IF(IST /= 0) THEN
       CALL TXPLST
    ELSE
       WRITE(6,'(A)') ' ## PARM INPUT ACCEPTED.'
    END IF
    REWIND(7)
    IERR=0

  END SUBROUTINE TXPARL

  SUBROUTINE TXPARF(KPNAME)

    integer(4) :: IST, KL
    logical :: LEX
    character(len=*) :: KPNAME

    INQUIRE(FILE=KPNAME,EXIST=LEX)
    IF(.NOT.LEX) RETURN

    OPEN(25,FILE=KPNAME,IOSTAT=IST,STATUS='OLD')
    IF(IST > 0) THEN
       WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
       RETURN
    END IF
    READ(25,TX,IOSTAT=IST)
    IF(IST > 0) THEN
       WRITE(6,*) 'XX PARM FILE READ ERROR'
       RETURN
    ELSE IF(IST < 0) THEN
       WRITE(6,*) 'XX PARM FILE EOF ERROR'
       RETURN
    END IF
    CALL KTRIM(KPNAME,KL)
    WRITE(6,*) &
         &     '## FILE (',KPNAME(1:KL),') IS ASSIGNED FOR PARM INPUT'
    IERR=0

  END SUBROUTINE TXPARF

  SUBROUTINE TXPARM_CHECK

    DO 
       ! System integers
       IF(NRMAX > NRM .OR. NRMAX < 0) EXIT
       IF(NQMAX > NQM .OR. NQMAX < 0) EXIT
       IF(NTMAX < 0) EXIT
       IF(NTSTEP < 0 .OR. NGRSTP < 0 .OR. NGTSTP < 0 .OR. NGVSTP < 0) EXIT
       ! Physical variables
       IF(RA >= RB) EXIT
       IF(RC < 0.D0 .OR. RC > RB) EXIT
       IF(RR <= RB) EXIT
       IF(rIPs < 0.D0 .OR. rIPe < 0.D0) EXIT
       IF(PA < 0.D0 .OR. PZ < 0.D0 .OR. Zeff < 1.D0) EXIT
       IF(PN0 < 0.D0 .OR. PNa < 0.D0) EXIT
       IF(PTe0 < 0.D0 .OR. PTea < 0.D0) EXIT
       IF(PTi0 < 0.D0 .OR. PTia < 0.D0) EXIT
       IF(De0 < 0.D0 .OR. Di0 < 0.D0) EXIT
       IF(rMue0 < 0.D0 .OR. rMui0 < 0.D0) EXIT
       IF(WPE0 < 0.D0 .OR. WPI0 < 0.D0 .OR. WPM0 < 0.D0) EXIT
       IF(Chie0 < 0.D0 .OR. Chii0 < 0.D0 .OR. ChiNC < 0.D0) EXIT
       IF(FSDFIX < 0.D0 .OR. FSCDBM < 0.D0 .OR. FSBOHM < 0.D0) EXIT
       IF(FSPCLD < 0.D0 .OR. FSPCLC < 0.D0) EXIT
       IF(FSLC < 0.D0 .OR. FSRP < 0.D0 .OR. FSNC < 0.D0) EXIT
       IF(FSHL < 0.D0 .OR. FSNF < 0.D0) EXIT
       IF(FSLP < 0.D0 .OR. FSLTE < 0.D0 .OR. FSLTI < 0.D0) EXIT
       IF(FSION < 0.D0 .OR. FSD01 < 0.D0 .OR. FSD02 < 0.D0) EXIT
       IF(MDLC /= 1 .AND. MDLC /= 2) EXIT
       IF(rG1 < 0.D0) EXIT
       IF(Eb < 0.D0 .OR. PNBHP < 0.D0 .OR. PNBHT1 < 0.D0 .OR. PNBHT2 < 0.D0) EXIT
       IF(ABS(PNBCD) > 1.D0 .OR. ABS(PNBMPD) > 1.D0) EXIT
       IF(RNBP0 > RB .OR. RNBP0 < 0.D0) EXIT
       IF(RNBT10 > RB .OR. RNBT10 < 0.D0) EXIT
       IF(RNBT20 > RB .OR. RNBT20 < 0.D0) EXIT
       IF(rNRFe < 0.D0 .OR. PRFHe < 0.D0) EXIT
       IF(rNRFi < 0.D0 .OR. PRFHi < 0.D0) EXIT
       IF(RRFe0 > RB .OR. RRFe0 < 0.D0) EXIT
       IF(RRFi0 > RB .OR. RRFi0 < 0.D0) EXIT
       IF(PN0s < 0.D0 .OR. V0 < 0.D0) EXIT
       IF(rGamm0 < 0.D0 .OR. rGamm0 > 1.D0) EXIT
       IF(rGASPF < 0.D0) EXIT
       IF(PNeDIV < 0.D0 .OR. PTeDIV < 0.D0 .OR. PTiDIV < 0.D0) EXIT
       IF(NTCOIL <= 0 .OR. DltRPn < 0.D0 .OR. DltRPn > 1.D0 .OR. kappa < 0.D0) EXIT
       IF(DT < 0.D0 .OR. EPS < 0.D0) EXIT
       IF(ICMAX < 0) EXIT
       IF(ADV < 0.D0 .OR. ADV > 1.D0) EXIT
       IF(tiny_cap < 0.D0) EXIT
       IF(CMESH0 < 0.D0 .OR. CMESH < 0.D0) EXIT
       IF(WMESH0 < 0.D0 .OR. WMESH < 0.D0) EXIT
       RETURN
    END DO

    WRITE(6,*) 'XX CONSISTENCY ERROR: PLEASE CHECK CONSISTENCY OF INPUT PARAMETERS.'
    STOP
    
  END SUBROUTINE TXPARM_CHECK

  !***** INPUT PARAMETER LIST *****

  SUBROUTINE TXPLST

    WRITE(6,601)
    RETURN

601 FORMAT(' ','# &TX : RA,RB,RC,RR,BB,PA,PZ,Zeff,'/ &
         &       ' ',8X,'PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ,,PROFN1,PROFN2,PROFT1,PROFT2,'/ &
         &       ' ',8X,'De0,Di0,VWpch0,rMue0,rMui0,WPM0,WPE0,WPI0,'/ &
         &       ' ',8X,'Chie0,Chii0,ChiNC,'/ &
         &       ' ',8X,'FSDFIX,FSCDBM,FSBOHM,FSPCLD,FSPCLC,PROFD,PROFC,'/ &
         &       ' ',8X,'FSCX,FSLC,FSRP,FSNF,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD01,FSD02,'/&
         &       ' ',8X,'MDLC,rLn,rLT,'/ &
         &       ' ',8X,'Eb,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBHP,PNBHT1,PNBHT2,'/ &
         &       ' ',8X,'PNBCD,PNBMPD,rNRFe,RRFew,RRFe0,PRFHe,Tqi0,rNRFe,RRFew,RRFe0,PRFHe,'/&
         &       ' ',8X,'PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV,'/ &
         &       ' ',8X,'NTCOIL,DltRPn,kappa,m_pol,n_tor,'/ &
         &       ' ',8X,'DT,EPS,ICMAX,ADV,tiny_cap,CMESH0,CMESH,WMESH0,WMESH,'/ &
         &       ' ',8X,'NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,'/ &
         &       ' ',8X,'DelR,DelN,'/ &
         &       ' ',8X,'Dmag0,RMAGMN,RMAGMX,EpsH,NCph,NCth,'/ &
         &       ' ',8X,'rG1,FSHL,Q0,QA,'/ &
         &       ' ',8X,'rIPs,rIPe,'/ &
         &       ' ',8X,'MODEG,gDIV,MODEAV,MODEGL,MDLPCK,MDLWTB'/ &
         &       ' ',8X,'MDLETA,MDFIXT,MDITSN,MDITST,MDINTT,MDINIT,IDIAG,IGBDF,MDSOLV' / & 
         &       ' ',8X,'MDLNBD,MDLMOM')
  END SUBROUTINE TXPLST

!***************************************************************
!
!   View input parameters
!
!***************************************************************

  SUBROUTINE TXVIEW

    WRITE(6,'((1X,A6," =",1PD9.2,3(2X,A6," =",1PD9.2)))') &
         &   'RA    ', RA    ,  'RB    ', RB    ,  &
         &   'RR    ', RR    ,  'BB    ', BB    ,  &
         &   'PA    ', PA    ,  'PZ    ', PZ    ,  &
         &   'PN0   ', PN0   ,  'PNa   ', PNa   ,  &
         &   'PTe0  ', PTe0  ,  'PTea  ', PTea  ,  &
         &   'PTi0  ', PTi0  ,  'PTia  ', PTia  ,  &
         &   'rIP   ', rIP   ,  'Zeff  ', Zeff  ,  &
         &   'PROFJ ', PROFJ ,  'PROFN1', PROFN1,  &
         &   'PROFN2', PROFN2,  'PROFT1', PROFT1,  &
         &   'PROFT2', PROFT2,  'CMESH0', CMESH0,  &
         &   'WMESH0', WMESH0,  'CMESH ', CMESH ,  &
         &   'WMESH ', WMESH ,  'ADV   ', ADV   ,  &
         &   'De0   ', De0   ,  'Di0   ', Di0   ,  &
         &   'rMue0 ', rMue0 ,  'rMui0 ', rMui0 ,  &
         &   'VWpch0', VWpch0,  'WPM0  ', WPM0  ,  &
         &   'WPE0  ', WPE0  ,  'WPI0  ', WPI0  ,  &
         &   'PROFD ', PROFD ,  'PROFC ', PROFC , &
         &   'Chie0 ', Chie0 ,  'Chii0 ', Chii0 ,  &
         &   'ChiNC ', ChiNC , &
         &   'FSDFIX', FSDFIX,  'FSCDBM', FSCDBM,  &
         &   'FSBOHM', FSBOHM,  'FSPCLD', FSPCLD,  &
         &   'FSPCLC', FSPCLC,  'FSCX  ', FSCX  ,  &
         &   'FSLC  ', FSLC  ,  'FSRP  ', FSRP  ,  &
         &   'FSNF  ', FSNF  ,  'FSNC  ', FSNC  ,  &
         &   'FSLP  ', FSLP  ,  'FSLTE ', FSLTE ,  &
         &   'FSLTI ', FSLTI ,  'FSION ', FSION ,  &
         &   'FSD01 ', FSD01 ,  'FSD02 ', FSD02 ,  &
         &   'rLn   ', rLn   ,  'rLT   ', rLT   ,  &
         &   'Eb    ', Eb    ,  &
         &   'RNBP  ', RNBP  ,  'RNBP0 ', RNBP0 ,  &
         &   'RNBT1 ', RNBT1 ,  'RNBT10', RNBT10,  &
         &   'RNBT2 ', RNBT2 ,  'RNBT20', RNBT20,  &
         &   'PNBHP ', PNBHP ,  'PNBMPD', PNBMPD,  &
         &   'PNBHT1', PNBHT1,  'PNBHT2', PNBHT2,  &
         &   'rNRFe ', rNRFe ,  'RRFew ', RRFew ,  &
         &   'RRFe0 ', RRFe0 ,  'PRFHe ', PRFHe ,  &
         &   'rNRFi ', rNRFe ,  'RRFiw ', RRFiw ,  &
         &   'RRFi0 ', RRFi0 ,  'PRFHi ', PRFHi ,  &
         &   'Tqi0  ', Tqi0  , &
         &   'rGamm0', rGamm0,  'V0    ', V0    ,  &
         &   'rGASPF', rGASPF,  'PNeDIV', PNeDIV,  &
         &   'PTeDIV', PTeDIV,  'PTiDIV', PTiDIV,  &
         &   'DltRPn', DltRPn,  'kappa ', kappa ,  &
         &   'PN0s  ', PN0s  ,  'EPS   ', EPS   ,  &
         &   'tiny  ', tiny_cap,'DT    ', DT    ,  &
         &   'rG1   ', rG1   ,  'Zeff  ', Zeff  ,  &
         &   'rIPs  ', rIPs  ,  'rIPe  ', rIPe  ,  &
         &   'FSHL  ', FSHL  ,  'EpsH  ', EpsH  ,  &
         &   'DMAG0 ', DMAG0 ,  'RMAGMN', RMAGMN,  &
         &   'RMAGMX', RMAGMX,  &
         &   'FSHL  ', FSHL  ,  'EpsH  ', EpsH  ,  &
         &   'Q0    ', Q0    ,  'QA    ', QA
    WRITE(6,'((" ",A6," =",I5,3(6X,A6," =",I5)))') &
         &   'NRMAX ', NRMAX ,  &
         &   'NTMAX ', NTMAX ,  'NTSTEP', NTSTEP,  &
         &   'NGRSTP', NGRSTP,  'NGTSTP', NGTSTP,  &
         &   'NGVSTP', NGVSTP,  'ICMAX ', ICMAX ,  &
         &   'MODEG ', MODEG ,  'MODEAV', MODEAV,  &
         &   'MODEGL', MODEGL,  'MDLPCK', MDLPCK,  &
         &   'MDLWTB', MDLWTB,  'MDLETA', MDLETA,  &
         &   'MDFIXT', MDFIXT,  'MDITSN', MDITSN,  &
         &   'MDITST', MDITST,  'MDINTT', MDINTT,  &
         &   'MDINIT', MDINIT,  'IDIAG ', IDIAG ,  &
         &   'IGBDF ', IGBDF,   'MDSOLV', MDSOLV,  &
         &   'NTCOIL', NTCOIL,  'MDLC  ', MDLC,    &
         &   'm_pol ', m_pol ,  'n_tor ', n_tor,   &
         &   'MDLNBD', MDLNBD,  'MDLMOM', MDLMOM,  &
         &   'NCph  ', NCph  ,  'NCth  ', NCth

    RETURN
  END SUBROUTINE TXVIEW
end module tx_parameter_control
