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

  !   Mesh accumulation radius (m)
  !     RACCUM coincides with RA when RACCUM is minus. (-1: default)
  !     RACCUM is valid when RACCUM is plus.
  RACCUM = - 1.d0

  !   Plasma major radius (m)
  RR = 1.3D0

  !   Toroidal magnetic field (T)
  BB = 1.3D0

  !   Plasma current start (MA)
  rIPs= 0.15D0

  !   Plasma current end (MA)
  rIPe= 0.15D0

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
  De0 = 0.05D0

  !   Ion-driven diffusion parameter
  Di0 = 0.D0

  !   Electron viscosity parameter
  rMue0 = 0.3D0

  !   Ion viscosity parameter
  rMui0 = 0.3D0

  !   Drift frequency parameter (omega/omega*e)
  WPM0 = 0.D0

  !  ***** Thermal diffusivity parameters *****

  !   Electron thermal diffusivity parameter (Chie/D)
  !     0 for fixed temperature profile
  Chie0 = 0.3D0

  !   Ion thermal diffusivity parameter (Chie/D)
  !     0 for fixed temperature profile
  Chii0 = 0.3D0

  !   ***** Turbulent transport control parameters *****

  !   Fixed transport coefficient parameter
  !     1 : particle transport
  !     2 : momentum transport
  !     3 : thermal transport
  FSDFIX(1:3) = 1.D0

  !   Anomalous turbulent transport models
  !     1 : CDBM model (Current diffusive ballooning mode)
  !     2 : CDIM model (Current diffusive interchange mode)
  !     3 : MMM95 (Multi-mode model)
  !   If MDANOM is negative, smoothing of the pressure gradient in the direction of time is off.
  MDANOM = 1

  !   Switch and amplification of turbulence coefficients (typically 0 or 1)
  !     1 : particle transport
  !     2 : momentum transport
  !     3 : thermal transport
  FSANOM(1:3) = 0.D0

  !   Effect of ExB shear stabilization
  FSCBSH = 0.D0

  !   Position of ETB shoulder (valid only when MDLETB /= 0)
  !     1 : particle transport
  !     2 : momentum transport
  !     3 : thermal transport
  RhoETB(1:3) = 0.9D0

  !   ==== CDBM transport coefficient parameters ============================
  !      The finer time step size, typically less than or equal to DT=5.D-4,
  !         is anticipated when using CDBM.

  !   Effect of magnetic curvature
  !     (It could destabilize a numerical robustness when FS2 > FS1 (see TRCOFS).)
  FSCBKP = 0.D0

  !   Effect of elongation (CDBM05 model)
  !      [M.Honda and A.Fukuyama NF 46 (2006) 580]
  FSCBEL = 1.D0

  !   Factor of E x B rotation shear
  rG1 = 24.D0

  !   =======================================================================

  !   Bohm transport coefficient parameter in SOL
  FSBOHM = 0.D0

  !   Pseud-classical particle transport coefficient parameter in SOL
  FSPCLD = 0.D0

  !   Pseud-classical mom. transport coefficient parameter in SOL
  FSPCLM = 0.D0

  !   Pseud-classical heat transport coefficient parameter in SOL
  FSPCLC = 0.D0

  !   Controller for the degree of the effect of Te'
  !     on the turbulent particle transport
  !   Only valid if MDVAHL = 1 or 2.
  FSVAHL = 0.D0

  !   Particle diffusion coefficient profile parameter (D(r=a)/D(r=0))
  PROFD =  3.D0

  !   Exponent of particle diffusion coefficient profile
  PROFD1 = 3.D0

  !   Gaussian modification of particle diffusion coefficient profile
  PROFD2 = 0.D0

  !   Mom. diffusion coefficient profile parameter (D(r=a)/D(r=0))
  PROFM = 10.D0

  !   Exponent of mom. diffusion coefficient profile
  PROFM1 = 2.D0

  !   Heat diffusion coefficient profile parameter (D(r=a)/D(r=0))
  PROFC = 10.D0

  !   Exponent of heat diffusion coefficient profile
  PROFC1 = 2.D0

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

  !   Thermal neutral diffusion factor
  FSD02 = 1.D0

  !   Halo neutral diffusion factor
  FSD03 = 1.D0

  !   Poloidal torque paramter due to neoclassical heat flux
  FSNCPL = 0.D0

  !   ***** initial parameters *****

  !   Initial Density scale length in SOL (m), valid if MDITSN /= 0
  rLn = 0.03D0

  !   Initail Temperature scale length in SOL (m), valid if MDITST /= 0
  rLT = 0.03D0

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

  !   NBI input power from input file (MW)
  PNBHex = 0.D0

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
  Tqt0  = 0.D0  ! Toroidal
  Tqp0  = 0.D0  ! Poloidal

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
!  EpsH = 0.1D0     ! amplitude 
  EpsHM(1:NHFMmx,0:3) = 0.d0
!  EpsHM(1,0:3) = 0.D0, 0.D0, 0.1D0, 0.D0     ! amplitude 
  EpsHM(1,2) = 0.1D0
!  NCph = 5        ! toroidal pitch number
!  NCth = 2        ! poloidal pitch number
  HPN(1:NHFMmx, 1:2) = 0
  HPN(1,1) = 2    ! poloidal pitch number
  HPN(1,2) = 5    ! toroidal pitch number
  Q0 = 3.D0       ! q(0) by external coils
  QA = 2.D0       ! q(a) by external coils

  !   ***** Numerical parameters *****

  !   Time step size(s)
!!!!  DT = 1.D-4
  DT = 1.D-3

  !   Convergence parameter
  !      Convergence calculation is going on until IC = ICMAX when EPS is negative.
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

  !   Permittivity switch for numerical convergence
  rMUb1 = rMU0
  rMUb2 = 1.d0

  !   Convergence model
  !     (Definition: ERROR means |X^j - X^{j-1}|, where j denotes the number of iteration.)
  !     0: Relative error 
  !          ERROR divided by root-mean-square,
  !          evaluated for each equation and each grid point.
  !     1: Relative error
  !          L2 norm of ERROR divided by L2 norm of X^{j-1},
  !          evaluated for each equation.
  MODECV = 0

  !   ***** Mesh number parameters *****

  !   Magnitude of mesh peakness
  CMESH0 =  2.D0
  CMESH  = 10.D0

  !   Width of mesh peakness
  WMESH0 = 0.2D0
  WMESH  = 5.D-2

  !   Number of nodes
  NRMAX  = 60

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

  !   Mode of AV (Diagnostic message in terms of convergence)
  !   0 : OFF (recommended)
  !   n : Interval of displaying diagnostic message
  MODEAV = 0

  !   Diagnostic parameter
  !   0 : OFF
  !   1 : debug message output (ntstep == 1)
  !   2 : debug message output (few)
  !   3 : debug message output (many)
  !   4 : debug message output (many,each prof)
  ! +10 : convergence output for each equation
  !  -1 : message for steady state check at NGTSTP intervals
  IDIAG = 0

  !   ***** Model parameters *****

  !   Mode of LAPACK
  !   0    : Use BANDRD
  !   else : Use LAPACK_DGBSV or LA_GBSV
  MDLPCK = 0

  !   Mode of orbit squeezing effect in NCLASS
  !   0    : No orbit squeezing effect
  !   1    : Orbit squeezing effect
  !   2    : Orbit squeezing effect, fixed during iteration
  MDOSQZ = 2

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

  !   Mode of initial density profiles
  !   -2   : read from file and smooth
  !   -1   : read from file
  !   0    : original
  MDINTN = 0

  !   Mode of initial temperature profiles
  !   -2   : read from file and smooth
  !   -1   : read from file
  !   0    : original
  !   1    : pedestal model
  !   2    : empirical steady state temperature profile
  MDINTT = 0

  !   Mode of initial current density profiles
  !   -2   : read from file and smooth
  !   -1   : read from file
  !   0    : original
  MDINTC = 0

  !   Mode of initial density profile in the SOL
  !   0    : polynominal model
  !   1    : exponential decay model
  MDITSN = 1

  !   Mode of initial temperature profile in the SOL
  !   0    : polynominal model
  !   1    : exponential decay model
  !   2    : exponential decay model 2
  !          This should be chosen if MDINTT=2.
  MDITST = 1

  !   Mode of initial plasma profils
  !   0    : minimal profiles are calculated
  !   1    : many profiles are analytically calculated 
  MDINIT = 1

  !   Mode of inherent convection annihilator
  !   0    : Nothing to do
  !   1    : Use parameter FSVAHL controlling the contribution of d(lnT)/dln(n)
  !          (Te gradient term is only changed according to FSVAHL.)
  !   2    : Annihilate the contribution of Er to the flux
  !   3    : Annihilate the inherent convection (= obtain pure diffusion)
  !          when FSNC = 0
  MDVAHL = 0

  !   Mode of Edge Transport barrier
  !   0    : Nothing to do
  !   1    : Turn off anomalous effect outside RhoETB
  MDLETB = 0

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
  gDIV(16) = 1.E16
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
  gDIV(35) = 1.E15
  gDIV(36) = 1.E13
  gDIV(37) = 1.E20
  gDIV(38) = 1.E3
  gDIV(39) = 1.E-3
  gDIV(41) = 1.E3
  gDIV(42) = 1.E3
  gDIV(45) = 1.E3
  gDIV(46) = 1.E3
  gDIV(50) = 1.E3
  gDIV(51) = 1.E3
  gDIV(53) = 1.E-2
  gDIV(55) = 1.E3
  gDIV(56) = 1.E-4
  gDIV(57) = 1.E-6
  gDIV(60) = 1.E6
  gDIV(64) = 1.E-20
  gDIV(65) = 1.E-20
  gDIV(72) = 1.E3
  gDIV(73) = 1.E3
  gDIV(74) = 1.E6
  gDIV(75) = 1.E6
  gDIV(85) = 1.E6
  gDIV(86) = 1.E6
  gDIV(89) = 1.E3
  gDIV(90) = 1.E20
  gDIV(92) = 1.E3
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
  gDIV(133) = 1.E20
  gDIV(136) = 1.E3
  gDIV(137) = 1.E-6
  gDIV(138) = 1.E13
  gDIV(139) = 1.E3
  gDIV(142) = 1.E27
  gDIV(143) = 1.E27

  !   *** Density perturbation technique ***

  !   Normalized radius where density increase
  DelRho = 0.2D0

  !   Amount of increase in density (10^20 m^-3)
  DelN = 1.D-1

  !   ***

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
  real(8)    :: MAXAMP, C1L, C2L, W1L, W2L, RL, RC, RCL, CLNEW

  !   Ion mass number
  AMI   = PA * AMP
  !   Beam ion mass number
  AMB   = AMI
  !   Number of equations
  NQMAX = NQM

  !   Square root permittivity for LQm1
  !     for the sake of acceleration of convergence
  sqeps0 = sqrt(EPS0)

  !  Mesh

  if(RACCUM < 0.d0) then
     RC = RA
  else
     RC = RACCUM
  end if

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

!!$  !  Maximum NR till RA
!!$
!!$  DO NR = 0, NRMAX-1
!!$     IF(R(NR) <= RA .AND. R(NR+1) >= RA) THEN
!!$        NRL = NR
!!$        EXIT
!!$     END IF
!!$  END DO

  !  Maximum NR till RC

  DO NR = 0, NRMAX-1
     IF(R(NR) <= RC .AND. R(NR+1) >= RC) THEN
        NRL = NR
        EXIT
     END IF
  END DO

  !  Adjust RC on mesh (for the time being RC is assumed to be equivalent to RA)

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

  !  RACCUM doesn't coincide with the plasma surface RA

  IF(RACCUM >= 0.D0) THEN

     !  Maximum NR till RA

     DO NR = 0, NRMAX-1
        IF(R(NR) <= RA .AND. R(NR+1) >= RA) THEN
           NRL = NR
           EXIT
        END IF
     END DO

     !  Adjust RA on mesh

     IF(ABS(R(NRL)-RA) < ABS(R(NRL+1)-RA)) THEN
        NRA = NRL
        RL = 0.5d0 * ABS(R(NRA) - RA)
        R(NRA  ) = RA
        R(NRA-1) = R(NRA-1) + RL
     ELSE
        NRA = NRL + 1
        RL = 0.5d0 * ABS(R(NRA) - RA)
        R(NRA  ) = RA
        R(NRA+1) = R(NRA+1) - RL
     END IF
     
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
  use tx_interface, only : INTG_P, INTDERIV3, detect_datatype, INTG_F, dfdx, &
       &                   initprof_input, moving_average

  implicit none
  INTEGER(4) :: NR, IER, ifile, NHFM, NR_smt, NR_smt_start = 10
  REAL(8) :: RL, PROF, PROFN, PROFT, PTePROF, PTiPROF!, QL, dRIP
  REAL(8) :: AJFCT, SUM_INT
  REAL(8) :: ALP, dPe, dPi, DR1, DR2
  REAL(8) :: EpsL, Vte, Wte, rNuAsE_inv, FTL, EFT, CR
  real(8) :: FACT, PBA, dPN, CfN1, CfN2, pea, pia, pediv, pidiv, dpea, dpia, &
       &     Cfpe1, Cfpe2, Cfpi1, Cfpi2, sigma, fexp, PN0L, PNaL, PNeDIVL, &
       &     PTe0L, PTi0L, PTeaL, PTiaL, PTeDIVL, PTiDIVL
  REAL(8) :: DERIV4, FCTR ! External functions
  real(8), dimension(:), allocatable :: AJPHL, TMP, RHSV, dPedr, dPidr, Prof1, Prof2
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

  if(MDINTN < 0 .or. MDINTT < 0 .or. ABS(MDINTC) /= 0) call initprof_input
  if(MDINTN < 0) then ! density at the boundaries
     call initprof_input(  0,1,PN0L)
     call initprof_input(NRA,1,PNaL)
     PNeDIVL = 0.25d0 * PNaL
  else
     PN0L = PN0
     PNaL = PNa
     PNeDIVL = PNeDIV
  end if
  if(MDINTT < 0) then ! temperature at the boundaries
     call initprof_input(  0,2,PTe0L)
     call initprof_input(  0,3,PTi0L)
     call initprof_input(NRA,2,PTeaL)
     call initprof_input(NRA,3,PTiaL)
     PTeDIVL = 0.25d0 * PTeaL
     PTiDIVL = 0.25d0 * PTiaL
  else
     PTe0L   = PTe0
     PTi0L   = PTi0
     PTeaL   = PTea
     PTiaL   = PTia
     PTeDIVL = PTeDIV
     PTiDIVL = PTiDIV
  end if
  PBA   = RB - RA
  dPN   = - 3.D0 * (PN0L - PNaL) / RA
  CfN1  = - (3.D0 * PBA * dPN + 4.D0 * (PNaL - PNeDIVL)) / PBA**3
  CfN2  =   (2.D0 * PBA * dPN + 3.D0 * (PNaL - PNeDIVL)) / PBA**4
  IF(MDFIXT == 0) THEN
     pea   = PNaL    * PTeaL   ; pia   = PNaL    / PZ * PTiaL
     dpea  = dPN     * PTeaL   ; dpia  = dPN     / PZ * PTiaL
     pediv = PNeDIVL * PTeDIVL ; pidiv = PNeDIVL / PZ * PTiDIVL
  ELSE
     pea   = PTeaL   ; pia   = PTiaL
     dpea  = 0.d0    ; dpia  = 0.d0
     pediv = PTeDIVL ; pidiv = PTiDIVL
  END IF
  Cfpe1 = - (3.D0 * PBA * dpea + 4.D0 * (pea - pediv)) / PBA**3
  Cfpe2 =   (2.D0 * PBA * dpea + 3.D0 * (pea - pediv)) / PBA**4
  Cfpi1 = - (3.D0 * PBA * dpia + 4.D0 * (pia - pidiv)) / PBA**3
  Cfpi2 =   (2.D0 * PBA * dpia + 3.D0 * (pia - pidiv)) / PBA**4

  DO NR = 0, NRMAX
     RL = R(NR)
     IF (RL <= RA) THEN ! +++ Core +++
        IF(MDINTN < 0) THEN
           call initprof_input(NR,1,X(LQe1,NR)) ! Ne
           X(LQi1,NR) = X(LQe1,NR) / PZ         ! Ni
        ELSE
           PROFN = (1.D0 - RHO(NR)**PROFN1)**PROFN2
           X(LQe1,NR) = (PN0L - PNaL) * PROFN + PNaL ! Ne
           X(LQi1,NR) = X(LQe1,NR) / PZ           ! Ni
        END IF
        IF(MDINTT < 0) THEN
           call initprof_input(NR,2,PTePROF) ! Te
           call initprof_input(NR,3,PTiPROF) ! Ti
        ELSE IF(MDINTT == 0) THEN
           PROFT = (1.D0 - RHO(NR)**PROFT1)**PROFT2
           PTePROF = (PTe0L - PTeaL) * PROFT + PTeaL
           PTiPROF = (PTi0L - PTiaL) * PROFT + PTiaL
        ELSE IF(MDINTT == 1) THEN
           PTePROF = (PTe0L - PTeaL) * (0.8263d0 * (1.d0 - RHO(NR)**2 )**1.5d0 &
                &                     + 0.167d0  * (1.d0 - RHO(NR)**30)**1.25d0) + PTeaL
           PTiPROF = (PTi0L - PTiaL) * (0.8263d0 * (1.d0 - RHO(NR)**2 )**1.5d0 &
                &                     + 0.167d0  * (1.d0 - RHO(NR)**30)**1.25d0) + PTiaL
        ELSE
           PTePROF = PTe0L * exp(- log(PTe0L / PTeaL) * RHO(NR)**2)
           PTiPROF = PTi0L * exp(- log(PTi0L / PTiaL) * RHO(NR)**2)
        END IF
        IF(MDFIXT == 0) THEN
           X(LQe5,NR) = PTePROF * X(LQe1,NR) ! Ne*Te
           X(LQi5,NR) = PTiPROF * X(LQi1,NR) ! Ni*Ti
        ELSE 
           X(LQe5,NR) = PTePROF ! Te
           X(LQi5,NR) = PTiPROF ! Ti
        END IF
     ELSE ! +++ SOL +++
        ! density
        IF(MDITSN == 0) THEN
           X(LQe1,NR) = PNaL + dPN * (RL - RA) + CfN1 * (RL - RA)**3 &
                &                              + CfN2 * (RL - RA)**4
           X(LQi1,NR) = X(LQe1,NR) / PZ
        ELSE
           X(LQe1,NR) = PNaL * EXP(- (RL - RA) / rLn)
           X(LQi1,NR) = X(LQe1,NR) / PZ
        END IF
        IF(MDITST == 0) THEN
           X(LQe5,NR) = pea + dpea * (RL - RA) + Cfpe1 * (RL - RA)**3 &
                &                              + Cfpe2 * (RL - RA)**4
           X(LQi5,NR) = pia + dpia * (RL - RA) + Cfpi1 * (RL - RA)**3 &
                &                              + Cfpi2 * (RL - RA)**4
        ELSE
           IF(MDITST == 1) THEN
              PTePROF = PTeaL * EXP(- (RL - RA) / rLT)
              PTiPROF = PTiaL * EXP(- (RL - RA) / rLT)
           ELSE
              sigma = 0.3d0 * (RHO(NRMAX) - 1.d0)
              fexp  = exp(- (RHO(NR) - 1.d0)**2 / (2.d0 * sigma**2))
              PTePROF = PTe0L   * exp(- log(PTe0L / PTeaL) * RHO(NR)**2) *         fexp &
                   &  + PTeDIVL                                          * (1.d0 - fexp)
              PTiPROF = PTi0L   * exp(- log(PTi0L / PTiaL) * RHO(NR)**2) *         fexp &
                   &  + PTiDIVL                                          * (1.d0 - fexp)
           END IF
           ! pressure (MDFIXT=0) or temperature (MDFIXT=1)
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
     ! N0_2 (thermal neutrals)
!     X(LQn2,NR) = 0.D0
     X(LQn2,NR) = 1.D-20 ! when ThntSW = 0.D0
     ! N0_3 (halo neutrals)
     X(LQn3,NR) = 0.D0
     ! Bphi
     X(LQm5,NR) = 0.5D0 * PSI(NR) * BB / rMU0
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

  IF(MDINTT == -2) THEN ! Smoothing temperatures
     allocate(Prof1(0:NRMAX),Prof2(0:NRMAX))
     Prof1(0:NRMAX) = X(LQe5,0:NRMAX) / X(LQe1,0:NRMAX)
     Prof2(0:NRMAX) = X(LQi5,0:NRMAX) / X(LQi1,0:NRMAX)
     NR_smt = NRA - NR_smt_start ! smoothing data only in the edge region
     DO NR = NR_smt, NRMAX
        X(LQe5,NR) = moving_average(NR,Prof1,NRMAX) * X(LQe1,NR)
        X(LQi5,NR) = moving_average(NR,Prof2,NRMAX) * X(LQi1,NR)
     END DO
     deallocate(Prof1,Prof2)
  END IF
  IF(MDINTN == -2) THEN ! Smoothing densities
     allocate(Prof1(0:NRMAX))
     X(LQe5,0:NRMAX) = X(LQe5,0:NRMAX) / X(LQe1,0:NRMAX)
     X(LQi5,0:NRMAX) = X(LQi5,0:NRMAX) / X(LQi1,0:NRMAX)
     Prof1(0:NRMAX) = X(LQe1,0:NRMAX)
     NR_smt = NRA - NR_smt_start ! smoothing data only in the edge region
     DO NR = NR_smt, NRMAX
        X(LQe1,NR) = moving_average(NR,Prof1,NRMAX)
        X(LQi1,NR) = X(LQe1,NR) / PZ         ! Ni
     END DO
     X(LQe5,0:NRMAX) = X(LQe5,0:NRMAX) * X(LQe1,0:NRMAX)
     X(LQi5,0:NRMAX) = X(LQi5,0:NRMAX) * X(LQi1,0:NRMAX)
     deallocate(Prof1)
  END IF

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
     Q(0) = Q0
     DO NR = 1, NRMAX
        RL = R(NR)
!        QL = (Q0 - QA) * (1.D0 - (RL / RA)**2) + QA
        Q(NR) = (Q0 - QA) * (1.D0 - (RL / RA)**2) + QA
        BthV(NR) = BB * RL / (Q(NR) * RR)
     END DO
  END IF
  Bthb = BthV(NRMAX)

  ! Toroidal electron current

  allocate(AJPHL(0:NRMAX))
  ifile = detect_datatype('LQe4')
  if(ifile == 0) then
     IF(MDINTC <= -1) THEN
        DO NR = 0, NRA
           call initprof_input(NR,4,AJPHL(NR)) ! NeUeph
        END DO
        AJPHL(NRA+1:NRMAX) = 0.D0
        AJFCT = rIPs * 1.D6 / (2.D0 * PI * INTG_F(AJPHL))
        ! Artificially extrapolate a current density in the SOL for numerical stability
        DO NR = NRA+1, NRMAX
           AJPHL(NR) = AJPHL(NRA) * EXP(- (R(NR) - RA) / (0.5d0 * rLn))
        END DO
        
        IF(MDINTC == -2) THEN ! Smoothing current density
           allocate(Prof1(0:NRMAX))
           Prof1(0:NRMAX) = AJPHL(0:NRMAX)
           NR_smt = NRA - NR_smt_start ! smoothing data only in the edge region
           DO NR = NR_smt, NRMAX
              AJPHL(NR) = moving_average(NR,Prof1,NRMAX)
           END DO
           deallocate(Prof1)
        END IF

        DO NR = 0, NRMAX
           AJPHL(NR) = AJFCT * AJPHL(NR)
           X(LQe4,NR) = - AJPHL(NR) / (AEE * 1.D20)
           AJOH(NR) = AJPHL(NR)
        END DO
        BthV(0) = 0.d0
        SUM_INT = 0.d0
        DO NR = 1, NRMAX
           SUM_INT = SUM_INT + INTG_P(AJPHL,NR,0)
           BthV(NR) = rMUb1 * SUM_INT / R(NR)
        END DO
     ELSE ! (MDINTC == 0)
        DO NR = 0, NRMAX
!!$           AJPHL(NR) = 2.D0 / rMUv1 * DERIV4(NR,PSI,R(0:NRMAX)*BthV(0:NRMAX),NRMAX,0)
!!$           X(LQe4,NR) = - AJPHL(NR) / (AEE * 1.D20)

           IF((1.D0-(R(NR)/RA)**2) <= 0.D0) THEN
              PROF= 0.D0    
           ELSE             
              PROF= 1.D0-(R(NR)/RA)**2
           END IF
           IF(FSHL == 0.D0) THEN
              ! Ne*UePhi
              AJPHL(NR)  = rIPs * 1.D6 / (PI * RA**2) * (PROFJ + 1.D0) * PROF**PROFJ
              X(LQe4,NR) = - AJPHL(NR) / (AEE * 1.D20)
              AJOH(NR)   = AJPHL(NR)
!              AJOH(NR)   = PROF
           ELSE
              AJPHL(NR)  = 0.D0
              X(LQe4,NR) = 0.D0
              AJOH(NR)   = 0.D0
           END IF
        END DO
     END IF
  else
     call inexpolate(infiles(ifile)%nol,infiles(ifile)%r,infiles(ifile)%data,NRMAX,RHO,5,AJPHL)
     AJFCT = rIPs * 1.D6 / (2.D0 * PI * INTG_F(AJPHL))
     DO NR = 0, NRMAX
        AJPHL(NR) = AJFCT * AJPHL(NR)
        X(LQe4,NR) = - AJPHL(NR) / (AEE * 1.D20)
        AJOH(NR) = AJPHL(NR)
     END DO
     BthV(0) = 0.d0
     SUM_INT = 0.d0
     DO NR = 1, NRMAX
        SUM_INT = SUM_INT + INTG_P(AJPHL,NR,0)
        BthV(NR) = rMUb1 * SUM_INT / R(NR)
     END DO
  end if

  if(MDINTN < 0 .or. MDINTT < 0 .or. ABS(MDINTC) /= 0) call initprof_input(idx = 0)

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

  ! LQm4

  IF((FSHL == 0.0) .AND. (ifile == 0) .AND. (MDINTC == 0) .AND. &
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

  ! Numerical solution for AphV

     RHSV(1:NRMAX) = - 0.5D0 * BthV(1:NRMAX) / R(1:NRMAX)
     X(LQm4,0) = 0.D0
     X(LQm4,1:NRMAX) = matmul(CMTX,RHSV)
  END IF

  ! Poloidal current density (Virtual current for helical system)

  IF(FSHL == 0.D0) THEN
     AJV(0:NRMAX)=0.D0

     Q(1:NRMAX) = ABS(R(1:NRMAX) * BB / (RR * BthV(1:NRMAX)))
     Q(0) = FCTR(R(1),R(2),Q(1),Q(2))
  ELSE
     ! Integrate 1 / (r * rMU0) * d/dr (r * BthV) to obtain AJV
!!$     allocate(TMP(0:NRMAX))
!!$     TMP(0:NRMAX) = R(0:NRMAX) * BthV(0:NRMAX)
!!$     DO NR = 1, NRMAX
!!$        dRIP = DERIV4(NR,R,TMP,NRMAX,0) * 2.D0 * PI / rMUb1
!!$        AJV(NR)=dRIP / (2.D0 * PI * R(NR))
!!$     END DO
!!$     AJV(0)=FCTR(R(1),R(2),AJV(1),AJV(2))
!!$     deallocate(TMP)
     AJV(0:NRMAX) = BB / (RR * rMU0) * 2.d0 * Q0 / Q(0:NRMAX)**2
  END IF

  ! Toroidal electric field for initial NCLASS calculation

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

  !  Check whether (m=0, n>0) Fourier component exists or not.   miki_m 10-09-07
  UHphSwitch = 0
  do NHFM = 1, NHFMmx 
     if (HPN(NHFM,1) == 0 .and. HPN(NHFM,2) > 0) UHphSwitch = 1
  enddo

  !  Define physical variables from X

  CALL TXCALV(X,0) ! Set variables as well as pres0 and ErV0
  !  --- Fixed Er to keep it constant during iterations ---
  ErV_FIX(0:NRMAX) = ErV(0:NRMAX)

  !  Calculate various physical quantities

  CALL TXCALC(0)

  !  Initial condition Part II

  IF(MDINIT == 1) THEN

     allocate(dPedr(0:NRMAX),dPidr(0:NRMAX))
     dPedr(0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PeV,NRMAX,0) * rKeV
     dPidr(0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PiV,NRMAX,0) * rKeV

     DO NR = 1, NRMAX
        dPe = dPedr(NR)
        dPi = dPidr(NR)
        IF(NR == NRMAX) THEN
           dPe = 0.D0
           dPi = 0.D0
        END IF
        IF(rNueNC(NR) == 0.D0) THEN
           X(LQe3,NR) = 0.D0
           X(LQi3,NR) = 0.D0
        ELSE
           ALP = (AMI / AME) * (rNuiNC(NR) / rNueNC(NR))
           X(LQi3,NR) =( (- BthV(NR) * X(LQe4,NR) + (dPe + dPi) / AEE) / BphV(NR) &
                &       +(FQeth(NR) + AMI / AME * FQith(NR)) /(rNueNC(NR) * 1.D20)) &
                &     / (PZ + ALP) * R(NR)
           X(LQe3,NR) =- ALP * X(LQi3,NR) + (FQeth(NR) + AMI / AME * FQith(NR)) &
                &                         /(rNueNC(NR) * 1.D20) * R(NR)
        END IF
        ErV(NR)    =- BphV(NR) / PNiV(NR) * X(LQi3,NR) / R(NR) + dPi / (PZ * AEE * PNiV(NR))
        X(LQe2,NR) = (rNueNC(NR) * X(LQe3,NR) - FQeth(NR) * 1.D-20 * R(NR)) &
             &     * AME / (AEE * BphV(NR))
        X(LQi2,NR) = X(LQe2,NR) / PZ
        X(LQm3,NR) = BthV(NR) / PNeV(NR) * X(LQe2,NR) / R(NR) &
             &     + AME * rNuei3(NR) / (AEE * PNeV(NR)) * X(LQe4,NR)
     END DO
     X(LQe3,0) = 0.D0 ; X(LQi3,0) = 0.D0
     ErV(0)    = 0.D0
     X(LQe2,0) = 0.D0 ; X(LQi2,0) = 0.D0
     X(LQm3,0) = AME * rNuei3(0) / (AEE * PNeV(0)) * X(LQe4,0)

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
     X(LQm5,1:NRMAX) = matmul(CMTX,RHSV) / rMU0
     deallocate(CMTX,RHSV,TMP)

     CALL TXCALV(X,0) ! Set variables as well as ErV0
     ErV_FIX(0:NRMAX) = ErV(0:NRMAX)

     CALL TXCALC(0)

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
       & RA,RB,RACCUM,RR,BB, &
       & PA,PZ,Zeff, &
       & PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ,PROFN1,PROFN2,PROFT1,PROFT2, &
       & De0,Di0,VWpch0,rMue0,rMui0,WPM0, &
       & Chie0,Chii0,ChiNC, &
       & FSDFIX,FSANOM,FSCBKP,FSCBEL,FSCBSH,FSBOHM,FSPCLD,FSPCLM,FSPCLC,FSVAHL,MDANOM,RhoETB, &
       & PROFD,PROFD1,PROFD2,PROFM,PROFM1,PROFC,PROFC1, &
       & FSCX,FSLC,FSRP,FSNF,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD01,FSD02,FSD03,FSNCPL,MDLC, &
       & rLn,rLT, &
       & Eb,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBHP,PNBHT1,PNBHT2,PNBCD,PNBMPD, &
       & rNRFe,RRFew,RRFe0,PRFHe,Tqt0,Tqp0,rNRFi,RRFiw,RRFi0,PRFHi, &
       & PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV, &
       & NTCOIL,DltRPn,kappa,m_pol,n_tor, &
       & DT,EPS,ICMAX,ADV,tiny_cap,CMESH0,CMESH,WMESH0,WMESH, &
       & NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP, &
       & DelRho,DelN, &
!       & DMAG0,RMAGMN,RMAGMX,EpsH,NCph,NCth, &
       & DMAG0,RMAGMN,RMAGMX, &
       & rG1,FSHL,Q0,QA, &
       & rIPs,rIPe, &
       & MODEG,gDIV,MODEAV,MODEGL,MDLPCK,MODECV, &
       & MDOSQZ,MDLETA,MDFIXT,MDITSN,MDITST,MDINTN,MDINTT,MDINTC,MDINIT,MDVAHL,MDLETB, &
       & IDIAG,IGBDF,MDSOLV,MDLNBD,MDLMOM, & ! 09/06/17~ miki_m
       & EpsHM, HPN  ! 10/08/06 miki_m
  private :: TXPLST

contains
!***************************************************************
!
!   Change input parameters
!
!***************************************************************

  SUBROUTINE TXPARM(KID)

    INTEGER(4) :: IST
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

    integer(4) :: idx

    idx = 1
    DO 
       ! /// idx = 1 - 10 ///
       ! System integers
       IF(NRMAX < 0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(NTMAX < 0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(NTSTEP < 0 .OR. NGRSTP < 0 .OR. NGTSTP < 0 .OR. NGVSTP < 0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       ! Physical variables
       IF(RA >= RB) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(RACCUM > RB) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(RR <= RB) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rIPs < 0.D0 .OR. rIPe < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(PA < 0.D0 .OR. PZ < 0.D0 .OR. Zeff < 1.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(PN0 < 0.D0 .OR. PNa < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(PTe0 < 0.D0 .OR. PTea < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       ! /// idx = 11 - 20 ///
       IF(PTi0 < 0.D0 .OR. PTia < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(De0 < 0.D0 .OR. Di0 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rMue0 < 0.D0 .OR. rMui0 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(Chie0 < 0.D0 .OR. Chii0 < 0.D0 .OR. ChiNC < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(minval(FSDFIX) < 0.D0 .OR. minval(FSANOM) < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(minval(RhoETB) < 0.D0 .OR. maxval(RhoETB) > 1.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(FSCBKP < 0.D0 .OR. FSCBEL < 0.D0 .OR. FSCBSH < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(FSBOHM < 0.D0 .OR. FSPCLD < 0.D0 .OR. FSPCLM < 0.D0 .OR. FSPCLC < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(FSLC < 0.D0 .OR. FSRP < 0.D0 .OR. FSNC < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(FSHL < 0.D0 .OR. FSNF < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       ! /// idx = 21 - 30 ///
       IF(FSLP < 0.D0 .OR. FSLTE < 0.D0 .OR. FSLTI < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(FSION < 0.D0)  THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(FSD01 < 0.D0 .OR. FSD02 < 0.D0 .OR. FSD03 < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(FSNCPL < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(MDLC /= 1 .AND. MDLC /= 2) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rG1 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(Eb < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(PNBHP < 0.D0 .OR. PNBHT1 < 0.D0 .OR. PNBHT2 < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(ABS(PNBCD) > 1.D0 .OR. ABS(PNBMPD) > 1.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(RNBP0 > RB .OR. RNBP0 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       ! /// idx = 31 - 40 ///
       IF(RNBT10 > RB .OR. RNBT10 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(RNBT20 > RB .OR. RNBT20 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rNRFe < 0.D0 .OR. PRFHe < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rNRFi < 0.D0 .OR. PRFHi < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(RRFe0 > RB .OR. RRFe0 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(RRFi0 > RB .OR. RRFi0 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(PN0s < 0.D0 .OR. V0 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rGamm0 < 0.D0 .OR. rGamm0 > 1.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rGASPF < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(PNeDIV < 0.D0 .OR. PTeDIV < 0.D0 .OR. PTiDIV < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       ! /// idx = 41 - 46 ///
       IF(NTCOIL <= 0 .OR. DltRPn < 0.D0 .OR. DltRPn > 1.D0 .OR. kappa < 0.D0) THEN
          EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(DT < 0.D0 .OR. EPS == 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(ICMAX < 0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(ADV < 0.D0 .OR. ADV > 1.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(tiny_cap < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(CMESH0 < 0.D0 .OR. CMESH < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(WMESH0 < 0.D0 .OR. WMESH < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF

       RETURN
    END DO

    WRITE(6,'(A,I3)') 'XX INPUT ERROR: Please check consistency of INPUT PARAMETERS. idx =',idx
    STOP
    
  END SUBROUTINE TXPARM_CHECK

  !***** INPUT PARAMETER LIST *****

  SUBROUTINE TXPLST

    WRITE(6,601)
    RETURN

601 FORMAT(' ','# &TX : RA,RB,RACCUM,RR,BB,PA,PZ,Zeff,'/ &
         &       ' ',8X,'PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ,,PROFN1,PROFN2,PROFT1,PROFT2,'/ &
         &       ' ',8X,'De0,Di0,VWpch0,rMue0,rMui0,WPM0,'/ &
         &       ' ',8X,'Chie0,Chii0,ChiNC,'/ &
         &       ' ',8X,'FSDFIX,FSANOM,FSCBKP,FSCBEL,FSCBSH,FSBOHM,FSPCLD,FSPCLM,FSPCLC,FSVAHL,'/ &
         &       ' ',8X,'MDANOM,RhoETB,'/ &
         &       ' ',8X,'PROFD,PROFD1,PROFD2,PROFM,PROFM1,PROFC,PROFC1,'/ &
         &       ' ',8X,'FSCX,FSLC,FSRP,FSNF,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD01,FSD02,FSD03,FSNCPL,'/&
         &       ' ',8X,'MDLC,rLn,rLT,'/ &
         &       ' ',8X,'Eb,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBHP,PNBHT1,PNBHT2,'/ &
         &       ' ',8X,'PNBCD,PNBMPD,rNRFe,RRFew,RRFe0,PRFHe,Tqt0,Tqp0,rNRFe,RRFew,RRFe0,PRFHe,'/&
         &       ' ',8X,'PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV,'/ &
         &       ' ',8X,'NTCOIL,DltRPn,kappa,m_pol,n_tor,'/ &
         &       ' ',8X,'DT,EPS,ICMAX,ADV,tiny_cap,CMESH0,CMESH,WMESH0,WMESH,'/ &
         &       ' ',8X,'NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,'/ &
         &       ' ',8X,'DelRho,DelN,'/ &
!         &       ' ',8X,'Dmag0,RMAGMN,RMAGMX,EpsH,NCph,NCth,'/ & 
         &       ' ',8X,'Dmag0,RMAGMN,RMAGMX,EpsHM,HPN,'/ &   ! 10/08/06 miki_m
         &       ' ',8X,'rG1,FSHL,Q0,QA,'/ &
         &       ' ',8X,'rIPs,rIPe,'/ &
         &       ' ',8X,'MODEG,gDIV,MODEAV,MODEGL,MDLPCK,MODECV,'/ &
         &       ' ',8X,'MDOSQZ,MDLETA,MDFIXT,MDITSN,MDITST,MDINTN,MDINTT,MDINIT,MDVAHL,MDLETB,' / & 
         &       ' ',8X,'IDIAG,IGBDF,MDSOLV,MDLNBD,MDLMOM')
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
         &   'PROFD ', PROFD ,  'PROFD1', PROFD1,  &
         &   'PROFD2', PROFD2,  'PROFM ', PROFM ,  &
         &   'PROFM1', PROFM1,  'PROFC ', PROFC ,  &
         &   'PROFC1', PROFC1,  'ChiNC ', ChiNC ,  &
         &   'Chie0 ', Chie0 ,  'Chii0 ', Chii0 ,  &
         &   'FSDFX1', FSDFIX(1),  'FSDFX2', FSDFIX(2),  &
         &   'FSDFX3', FSDFIX(3),  'FANOM1', FSANOM(1),  &
         &   'FANOM2', FSANOM(2),  'FANOM3', FSANOM(3),  &
         &   'RoETB1', RhoETB(1),  'RoETB2', RhoETB(2),  &
         &   'RoETB3', RhoETB(3),  'FSCBKP', FSCBKP,  &
         &   'FSCBEL', FSCBEL,  'FSCBSH', FSCBSH,  &
         &   'FSBOHM', FSBOHM,  'FSPCLD', FSPCLD,  &
         &   'FSPCLM', FSPCLM,  'FSPCLC', FSPCLC,  &
         &   'FSVAHL', FSVAHL,  'FSCX  ', FSCX  ,  &
         &   'FSLC  ', FSLC  ,  'FSRP  ', FSRP  ,  &
         &   'FSNF  ', FSNF  ,  'FSNC  ', FSNC  ,  &
         &   'FSLP  ', FSLP  ,  'FSLTE ', FSLTE ,  &
         &   'FSLTI ', FSLTI ,  'FSION ', FSION ,  &
         &   'FSD01 ', FSD01 ,  'FSD02 ', FSD02 ,  &
         &   'FSD03 ', FSD03 ,  'FSNCPL', FSNCPL,  &
         &   'rLn   ', rLn   ,  'rLT   ', rLT   ,  &
         &   'Eb    ', Eb    ,  &
         &   'RNBP  ', RNBP  ,  'RNBP0 ', RNBP0 ,  &
         &   'RNBT1 ', RNBT1 ,  'RNBT10', RNBT10,  &
         &   'RNBT2 ', RNBT2 ,  'RNBT20', RNBT20,  &
         &   'PNBHP ', PNBHP ,  'PNBHT1', PNBHT1,  &
         &   'PNBHT2', PNBHT2,  'PNBHex', PNBHex,  &
         &   'PNBMPD', PNBMPD,  'rNRFe ', rNRFe ,  &
         &   'RRFew ', RRFew ,  'RRFe0 ', RRFe0 ,  &
         &   'PRFHe ', PRFHe ,  'rNRFi ', rNRFe ,  &
         &   'RRFiw ', RRFiw ,  'RRFi0 ', RRFi0 ,  &
         &   'PRFHi ', PRFHi ,  &
         &   'Tqt0  ', Tqt0  ,  'Tqp0  ', Tqp0  ,  &
         &   'rGamm0', rGamm0,  'V0    ', V0    ,  &
         &   'rGASPF', rGASPF,  'PNeDIV', PNeDIV,  &
         &   'PTeDIV', PTeDIV,  'PTiDIV', PTiDIV,  &
         &   'DltRPn', DltRPn,  'kappa ', kappa ,  &
         &   'PN0s  ', PN0s  ,  'EPS   ', EPS   ,  &
         &   'tiny  ', tiny_cap,'DT    ', DT    ,  &
         &   'rG1   ', rG1   ,  'Zeff  ', Zeff  ,  &
         &   'rIPs  ', rIPs  ,  'rIPe  ', rIPe  ,  &
         &   'DMAG0 ', DMAG0 ,  'RMAGMN', RMAGMN,  &
         &   'RMAGMX', RMAGMX,  &
         &   'FSHL  ', FSHL  ,    &  ! Too many elements of EpsHM to show miki_m 10-08-11
         &   'Q0    ', Q0    ,  'QA    ', QA    ,  &
!!$         &   'EpsHM(1,:) ', EpsHM(1,0:3) ,         &
!!$         &   'EpsHM(2,:) ', EpsHM(2,0:3) ,         &
!!$         &   'EpsHM(3,:) ', EpsHM(3,0:3)
         &   'EpsH10 ', EpsHM(1,0), 'EpsH11 ', EpsHM(1,1), &
         &   'EpsH12 ', EpsHM(1,2), 'EpsH13 ', EpsHM(1,3), &
         &   'EpsH20 ', EpsHM(2,0), 'EpsH21 ', EpsHM(2,1), &
         &   'EpsH22 ', EpsHM(2,2), 'EpsH23 ', EpsHM(2,3), &
         &   'EpsH30 ', EpsHM(3,0), 'EpsH31 ', EpsHM(3,1), &
         &   'EpsH32 ', EpsHM(3,2), 'EpsH33 ', EpsHM(3,3), &
         &   'EpsH40 ', EpsHM(4,0), 'EpsH41 ', EpsHM(4,1), &
         &   'EpsH42 ', EpsHM(4,2), 'EpsH43 ', EpsHM(4,3)
    WRITE(6,'((" ",A6," =",I5,3(6X,A6," =",I5)))') &
         &   'NRMAX ', NRMAX ,  &
         &   'NTMAX ', NTMAX ,  'NTSTEP', NTSTEP,  &
         &   'NGRSTP', NGRSTP,  'NGTSTP', NGTSTP,  &
         &   'NGVSTP', NGVSTP,  'ICMAX ', ICMAX ,  &
         &   'MODEG ', MODEG ,  'MODEAV', MODEAV,  &
         &   'MODEGL', MODEGL,  'MDLPCK', MDLPCK,  &
         &   'MODECV', MODECV,  &
         &   'MDOSQZ', MDOSQZ,  'MDLETA', MDLETA,  &
         &   'MDANOM', MDANOM,  'MDFIXT', MDFIXT,  &
         &   'MDITSN', MDITSN,  'MDITST', MDITST,  &
         &   'MDINTN', MDINTN,  'MDINTT', MDINTT,  &
         &   'MDINTC', MDINTC,  'MDINIT', MDINIT,  &
         &   'MDVAHL', MDVAHL,  'MDLETB', MDLETB,  &
         &   'IDIAG ', IDIAG ,  'IGBDF ', IGBDF,   &
         &   'MDSOLV', MDSOLV,  &
         &   'NTCOIL', NTCOIL,  'MDLC  ', MDLC,    &
         &   'm_pol ', m_pol ,  'n_tor ', n_tor,   &
         &   'MDLNBD', MDLNBD,  'MDLMOM', MDLMOM ,  &
!         &   'NCph  ', NCph  ,  'NCth  ', NCth,    &
         &   'HPNth1   ', HPN(1,1), 'HPNph1   ', HPN(1,2), &
         &   'HPNth2   ', HPN(2,1), 'HPNph2   ', HPN(2,2), &
         &   'HPNth3   ', HPN(3,1), 'HPNph3   ', HPN(3,2), &
         &   'HPNth4   ', HPN(4,1), 'HPNph4   ', HPN(4,2)
    RETURN
  END SUBROUTINE TXVIEW
end module tx_parameter_control
