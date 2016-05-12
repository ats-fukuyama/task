!     $Id: txinit.f90,v 1.94 2011/06/13 07:53:20 honda Exp $
!***************************************************************
!
!   Set constants and initial parameters
!
!***************************************************************

SUBROUTINE TXINIT
  use tx_commons
  use tx_graphic, only : MODEG, MODEGL, NGYRM, NGR, NGRSTP, NGTSTP, NGVSTP, gDIV

  implicit none

  !   Number of equations
  NQMAX = NQM

  !   ***** Configuration parameters *****

  !   Plasma minor radius (m), geometrically defined by (Rmax-Rmin)/2
  ra = 0.8D0

  !   Plasma minor radius (m), defined by sqrt(V/(2 Pi Pi R))
  ravl = 0.8D0

  !   Virtual wall radius in rho coordinate (-) for defining rb
  !     The position of the virtual wall follows the change in that of the separatrix
  !     when an equilibrium evolves.
  !     That is, the virtual wall always locates at rhob distance from the separatrix.
  rhob = 1.1d0

  !   Mesh accumulation radius (-)
  !     rhoaccum coincides with rho=1.d0 when rhoaccum is minus. (-1: default)
  !     rhoaccum is valid when rhoaccum is plus.
  rhoaccum = - 1.d0

  !   Plasma major radius (m), geometrically defined by (Rmax+Rmin)/2
  RR = 3.2D0

  !   Toroidal magnetic field (T) at R=RR
  BB = 2.68D0

  !   Poloidal current function at the virtual wall (Tm)
  rbvt = rr * bb

  !   Plasma current start (MA)
  rIPs= 1.D0

  !   Plasma current end (MA)
  rIPe= 1.D0

  !   ***** Plasma components *****

  !   Atomic number
  amas(1) = aep     ! electron
  amas(2) = 2.d0    ! bulk ion
  amb     = amas(2) ! beam ion

  !   Charge number
  achg(1) = -1.d0 
  achg(2) =  1.d0 
  achgb   =  1.d0 

  !   Effective charge
  Zeff = 1.D0

  !   ***** Initial plasma parameters *****

  !   Initial electron density at rho = 0 (10^20 m^-3)
  PN0 = 0.2D0

  !   Initial electron density at rho = a (10^20 m^-3)
  PNa = 0.05D0

  !   Electron density in divertor region (Minimum density in SOL)
  PNeDIV = 0.01D0

  !   Ion density in divertor region (Minimum density in SOL)
  PNiDIV = PNeDIV

  !   Initial electron temperature at rho = 0 (keV)
  PTe0 = 2.D0

  !   Initial electron temperature at rho = a (keV)
  PTea = 0.2D0

  !   Electron temperature in divertor region (Minimum Te in SOL)
  PTeDIV = 0.05D0

  !   Initial ion temperature  at rho = 0 (keV)
  PTi0 = 2.D0

  !   Initial ion temperature  at rho = a (keV)
  PTia = 0.2D0

  !   Ion temperature in divertor region (Minimum Ti in SOL)
  PTiDIV = 0.05D0

  !   Initial current profile parameter
  PROFJ = 2.D0

  !   Initial density profile parameters
  PROFN1 = 2.D0
  PROFN2 = 1.D0

  !   Initial temperature profile parameters
  PROFT1 = 2.D0
  PROFT2 = 2.D0

  !   Initial ion toroidal rotation velocity at rho = 0 (m/s)
  Uiph0 = 0.d0

  !   Equilibrium parameters
  !     ieqread = 0 : Large aspect ratio limit
  !             = 1 : Large aspect ratio approximation
  !             = 2 : Read equilibrium parameters (not yet)
  ieqread = 1

  !   ***** Particle diffusivity and viscosity parameters *****

  !   Turbulent pinch velocity parameter
  VWpch0 = 0.D0

  !   Electron-driven diffusion parameter
  De0 = 0.1D0

  !   Ion-driven diffusion parameter (usually inactive)
  Di0 = 0.D0

  !   Electron viscosity parameter
  rMue0 = 0.5D0

  !   Ion viscosity parameter
  rMui0 = 0.5D0

  !   Drift frequency parameter (omega/omega*e)
  WPM0 = 0.D0

  !  ***** Thermal diffusivity parameters *****

  !   Electron thermal diffusivity parameter (Chie/D)
  !     0 for fixed temperature profile
  Chie0 = 0.5D0

  !   Ion thermal diffusivity parameter (Chie/D)
  !     0 for fixed temperature profile
  Chii0 = 0.5D0

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
  !     0 : No ExB shear stabilization
  !     1 : Lorentz-type ExB shear stabilization
  !    -1 : Exponent-type ExB shear stabilization [M. Honda (2007) dissertation, Chapter 4]
  FSCBSH = 0.D0

  !   Position of ETB shoulder (valid only when MDLETB /= 0)
  !     1 : particle transport
  !     2 : momentum transport
  !     3 : thermal transport
  RhoETB(1:3) = 0.9D0

  !   Momentum pinch parameter
  !     0 : No momentum pinch
  !     1 : Momentum pinch effective
  FSMPCH(1:NSM) = 1.d0

  !   ==== CDBM transport coefficient parameters ============================
  !      The finer time step size, typically less than or equal to DT=5.D-4,
  !         is anticipated when using CDBM.

  !   Effect of magnetic curvatures-alpha
  FSCBAL = 1.D0

  !   Effect of magnetic curvature
  !     (It could destabilize a numerical robustness when FS2 > FS1. )
  FSCBKP = 0.D0

  !   Effect of elongation (CDBM05 model)
  !      [M.Honda and A.Fukuyama NF 46 (2006) 580]
  FSCBEL = 1.D0

  !   Factor of E x B rotation shear
  rG1 = 10.D0

  !   ==== Neoclassical transport coefficient parameters ====================

  !   Projection of toroidal momentum flux onto parallel flow equations
  !     These terms preserve consistency that has to be satisfied by LQ*2 with
  !     LQ*3 and LQ*4 at the magnetic axis, where terms of grad P and grad Phi
  !     become nil.
  FSPARV(1:NSM) = 1.d0

  !   Neoclassical thermal diffusivity parameter
  ChiNC = 1.D0

  !   Neoclassical viscosity parameter
  FSNC = 1.D0

  !   Beam neoclassical viscosity parameter
  FSNCB = 1.D0

  !   Helical neoclassical viscosity parameter
  FSHL = 0.D0

  !   =======================================================================

  !   Bohm transport coefficient parameter in SOL
  FSBOHM = 1.D0

  !   Pseud-classical particle transport coefficient parameter in SOL
  FSPCLD = 1.D0

  !   Pseud-classical mom. transport coefficient parameter in SOL
  FSPCLM = 0.D0

  !   Pseud-classical heat transport coefficient parameter in SOL
  FSPCLC = 0.D0

  !   Controller for thermodiffusive pinch term of turbulent particle flux
  FSVAHL =-0.5D0

  !   Particle diffusion coefficient profile parameter (D(r=a)/D(r=0))
  PROFD =  5.D0

  !   Exponent of particle diffusion coefficient profile
  PROFD1 = 2.D0

  !   Gaussian modification of particle diffusion coefficient profile
  PROFD2 = 0.D0

  !   Pedestal of particle diffusion coefficient profile
  PROFDB = 0.D0

  !   Mom. diffusion coefficient profile parameter (D(r=a)/D(r=0))
  PROFM = 10.D0

  !   Exponent of mom. diffusion coefficient profile
  PROFM1 = 2.D0

  !   Pedestal of mom. diffusion coefficient profile
  PROFMB = 0.D0

  !   Heat diffusion coefficient profile parameter (chi(r=a)/chi(r=0))
  PROFC = 10.D0

  !   Exponent of heat diffusion coefficient profile
  PROFC1 = 2.D0

  !   Pedestal of heat diffusion coefficient profile
  PROFCB = 0.D0

  !   ***** Other transport parameters *****

  !   Advection parameter
  FSADV = 1.d0

  !   Beam advection parameter (finer time step required when on)
  FSADVB = 0.d0

  !   Velocity of flux surfaces
  FSUG = 1.d0

  !   Charge exchange parameter
  FSCX = 1.D0

  !   Orbit loss parameter
  !     FSLC = 0 : No orbit loss included.
  !            1 : Loss term is expressed as damping term.
  !          +10 : Loss term is expressed as source and sink term.
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

  !   ***** initial parameters *****

  !   Initial Density scale length in SOL normalized by minor radius (-), valid if MDITSN /= 0
  rLn = 0.0857D0

  !   Initail Temperature scale length in SOL normalized by minor radius (-), valid if MDITST /= 0
  rLT = 0.0857D0

  !   ***** Heating parameters *****

  !   Maximum NBI beam energy (keV)
  Ebmax = 80.D0

  !   Fraction of particles with Ebmax, Ebmax/2 and Ebmax/3 energies
  !      Positive-ion-source NBI : 0.75, 0.15, 0.10 (typically)
  !      Negative-ion-source NBI : 1.00, 0.00, 0.00
  esps(1) = 0.75d0
  esps(2) = 0.15d0
  esps(3) = 0.10d0

  !   Heating radius of perp NBI heating (-)
  RNBP  = 0.5d0

  !   Heating center of perp NBI heating (-)
  RNBP0 = 0.d0

  !   Heating radius of first tangential NBI heating (-)
  RNBT1  = 0.5d0

  !   Heating radius of second tangential NBI heating (-)
  RNBT2  = 0.5d0

  !   Heating center of first tangential NBI heating (-)
  RNBT10 = 0.2d0

  !   Heating center of second tangential NBI heating (-)
  RNBT20 = 0.d0

  !   Perpendicular NBI input power (MW)
  PNBHP = 0.D0

  !   First tangential NBI input power (MW)
  PNBHT1 = 0.D0

  !   Second tangential NBI input power (MW)
  PNBHT2 = 0.D0

  !   NBI input power from input file (MW)
  PNBHex = 0.D0

  !   NBI current drive parameter ; Direction of injection
  !     1.0 : co, -1.0 : counter
  PNBCD = 1.D0

  !   Fraction of the collisional slowing down part of the perpendicular NBI
  !     -1.0 (counter) <= PNBMPD <= 1.0 (co)
  !     0.0 : Exact perpendicular injection without momentum input
  PNBMPD = 0.D0

  !   Pitch angle between the tanjential NBI chord and the magnetic field line
  PNBPTC = 0.8d0

  !   Different NBI deposition profiles for electrons and ions due to banana orbit effect
  !     MDLNBD = 0 : No charge separation
  !              1 : Orbit effect for banana particles only
  !              2 : Orbit effect for all beam ions including passing particles
  MDLNBD = 0

  !   Refractive index of RF waves
  rNRFe = 0.D0
  rNRFi = 0.D0

  !   Heating width of RF heating (-)
  RRFew = 0.5d0
  RRFiw = 0.5d0

  !   Heating center radius of RF heating (-)
  RRFe0 = 0.d0
  RRFi0 = 0.d0

  !   RF input power (MW)
  PRFHe = 0.D0
  PRFHi = 0.D0

  !   Additional torque input (N m)
  Tqt0  = 0.D0  ! Toroidal
  Tqp0  = 0.D0  ! Poloidal

  !   ***** Neutral parameters *****

  !   Initial Neutral density (10^20 m^-3)
  PN0s = 1.D-8

  !   Neutral thermal velocity (m/s)
  !     V0 = SQRT(2.D0*X[eV]*AEE/(amas(2)*AMP))
  V0 = 1.6954D4

  !   Recycling rate in SOL
  rGamm0 = 0.8D0

  !   Gas-puff particle flux (10^20 m^-2 1/s)
  !      If you input Gamma_0 [particles/sec], you translate it into
  !      rGASPF [1/(m^2 s)]:
  !        rGASPF = Gamma_0 / <|nabla V|>
  !               = Gamma_0 / surt(NRMAX)
  !               ~ Gamma_0 / (2.D0*PI*RR*2.D0*PI*RB)
  rGASPF = 0.2D0

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

  !   Cumulative number of time steps
  NTCUM = 0

  !   Time step size(s)
!!!!  DT = 1.D-4
  DT = 1.D-3

  !   Convergence parameter
!!  !      Convergence calculation is going on until IC = ICMAX when EPS is negative.
  !     Computation continues even if not converged at IC = ICMAX when EPS is positive.
  EPS = 1.D-3

  !   Iteration
  ICMAX = 50

  !   Time-advancing method
  !     ADV = 0     : Explicit scheme       (Not usable)
  !           0.5   : Crank-Nicolson scheme (Not usable)
  !           2/3   : Galerkin scheme       (Not usable)
  !           0.729 : Minimum value at which a simulation can be run.
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

  !   Mode of the scheme to evaluate dPhi/dV and dPs/dpsi for suppressing wiggles near the axis
  !     ISMTHD = 0  : dPhiV/dV and dPs/dpsi
  !              1  : dPhiV/drho with dPhiV/drho = 0 at rho = 0 and dPs/drho with dPs/drho at rho = 0
  !              11 : dPhiV/drho(1) and dPs/drho(1) are interpolated by "replace_interpolate_value".
  !                   Smooth profiles regarding diamag. flow would be obtained.
  ISMTHD = 11

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

  !   Ratio of previous dependent values to next dependent values for convergence
  !   Increasing this value generally stabilizes oscilllation during convergence, 
  !   but requires more iteration numbers even when it is well behaved.
  !     0 <= oldmix <= 1
  oldmix = 0.d0

  !   SUPG on/off switch
  !     iSUPG2 for LQe2CC and LQi2CC
  !     iSUPG3 for LQe3CC and LQi3CC; LQb3CC only when iSUPG3 = -1
  !     iSUPG6 for LQe6CC and LQi6CC
  iSUPG2 = 0 ! Should be on when the turbulent particle flux is off (De0=0.0).
  iSUPG3 = 0
  iSUPG6 = 0

  !   Enhancement factor for numerical stability in advection problems
  !     p = SUPGstab / sqrt(15.d0)
  !     where SUPGstab = 1.d0 : [Raymond and Garder, Monthly Weather Review (1976)]
  !                    = 2.d0 : [Ganjoo and Tezduyar, NASA report CR-180124 (1986)]
  !                    = 4.d0 : Seemingly best, but no ground in the light of numerical analysis
  SUPGstab = 2.d0

  !   ***** Mesh number parameters *****

  !   Magnitude of mesh peakness
  CMESH0 =  2.D0
  CMESH  = 10.D0

  !   Width of mesh peakness (-)
!  WMESH0 = 0.2D0
!  WMESH  = 5.D-2
  WMESH0 = 0.25d0
  WMESH  = 0.0625d0

  !   Number of nodes
  NRMAX  = 60

  !   Number of elements
  NEMAX = NRMAX

  !   ***** Time parameter *****

  !   Number of time step
  NTMAX = 100

  !   ***** Diagnostics parameters *****

  !   Time step interval between print output
  NTSTEP = 50

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

  !   ***** Graphic parameters (module tx_graphic) *****

  !   Time step interval between lines in f(r) graph
  NGRSTP = 20

  !   Time step interval between points in f(t) graph
  NGTSTP = 5

  !   Time step interval between points in f(t) graph
  NGVSTP = 5

  !   Mode of Graph
  !   0 : for Display (with grid, w/o power)
  !   1 : for Display (with grid and power)
  !   2 : for Print Out (w/o grid, with power)
  MODEG = 2

  !   MODE of Graph Line
  !   0 : Change Line Color (Last Color Fixed)
  !   1 : Change Line Color and Style
  !   2 : Change Line Color and Style (With Legend)
  !   3 : Change Line Color, Style and Mark
  !   4 : Change Line Color, Style and Mark (With Legend)
  MODEGL = 1

  !   ***** Model parameters *****

  !   Mode of LAPACK
  !   0    : Use BANDRD
  !   else : Use LAPACK_DGBSV or LA_GBSV
  MDLPCK = 1

  !   Mode of fixed temperature profile
  !   0    : not fixed
  !   1    : fixed
  MDFIXT = 0

  !   Mode of equations solved for beam ions
  !   0    : minimum ; LQb1,       LQb3, LQb4, LQb7
  !   1    : more    ; LQb1, LQb2, LQb3, LQb4, LQb7
  MDBEAM = 1

  !   Mode of how to calculate the derivative of Phi w.r.t. psi for orbit squeezing effect
  !   0    : Calculate the derivative in an usual manner
  !   1    : The derivative is fixed during iteration
  !   +10  : Smoothing d/dpsi(dPhi/dpsi) profile using moving_average
  MDOSQZ = 0

  !   Model of neoclassical resistivity model (mainly for graphics)
  !   1    : depending upon MDLNEO
  !   2    : Sauter model
  !   3    : Hirshman, Hawryluk and Birge
  MDLETA = 1

  !   Model of neoclassical transport model, espcially for calculating
  !     friction coefficients and viscosities
  !   1    : Matrix Inversion
  !   2    : NCLASS
  !   +10  : Both (Users cannot specify)
  MDLNEO  = 1

  !   Mode of orbit squeezing effect for neoclassical transport solver
  !   0    : not considered
  !   1    : considered
  MDOSQZN = 0

  !   Choice of taking whether the bootstrap (BS) current or the resistivity
  !     from the neoclassical transport solver determined by MDLNEO
  !   This parameter is used for estimating the fraction of the BS current
  !     and the ohmic current, because TASK/TX cannot decompose the components
  !     of the current.
  !   0    : BS current is calculated by the external module.
  !          Ohmic current is subserviently determined. (strongly recommended)
  !   else : Ohmic current is calculated using the resistivity calculated by
  !          the external module. BS current is subserviently determined.
  !          NOTE: This option may cause the estimate of the negative BS current
  !                when NBs are injected.
  MDBSETA  = 0

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
  MDINTT = 2

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
  MDITST = 2

  !   Mode of Edge Transport barrier
  !   0    : Nothing to do
  !   1    : Turn off anomalous effect outside RhoETB
  MDLETB = 0

  !   Multiplication factor for graphic in the radial direction (module tx_graphic)
  !   default : 1.0
  !
  gDIV(1:NGYRM) = 1.0
  gDIV(1)  = 1.E20
  gDIV(2)  = 1.E14
  gDIV(4)  = 1.E3
  gDIV(5)  = 1.E3
  gDIV(7:9)  = 1.E3
  gDIV(12) = 1.E18
  gDIV(13) = 1.E3
  gDIV(16) = 1.E16
  gDIV(17) = 1.E-3
  gDIV(19) = 1.E3
  gDIV(21) = 1.E6
  gDIV(22) = 1.E6
  gDIV(23:28) = 1.E3
  gDIV(35) = 1.E15
  gDIV(36) = 1.E13
  gDIV(37) = 1.E20
  gDIV(38) = 1.E3
  gDIV(39) = 1.E-3
  gDIV(42) = 1.E3
  gDIV(45) = 1.E3
  gDIV(46) = 1.E3
  gDIV(48) = 1.E3
  gDIV(49) = 1.E3
  gDIV(52) = 1.E3
  gDIV(53) = 1.E3
  gDIV(55) = 1.E3
  gDIV(56) = 1.E3
  gDIV(57) = 1.E-2
  gDIV(58) = 1.E3
  gDIV(59) = 1.E-4
  gDIV(60) = 1.E-6
  gDIV(61) = 1.E-6
  gDIV(62) = 1.E-4
  gDIV(64) = 1.E-20
  gDIV(65) = 1.E-20
  gDIV(72) = 1.E3
  gDIV(73) = 1.E3
  gDIV(85) = 1.E6
  gDIV(86) = 1.E6
  gDIV(89) = 1.E3
  gDIV(90) = 1.E20
  gDIV(92) = 1.E3
  gDIV(96:98) = 1.E6
  gDIV(100:104) = 1.E6
  gDIV(109) = 1.E3
  gDIV(110) = 1.E15
!  gDIV(115) = 1.E15
  gDIV(121) = 1.E6
  gDIV(123) = 1.E6
  gDIV(124) = 1.E6
  gDIV(133) = 1.E20
  gDIV(136) = 1.E20
  gDIV(137) = 1.E20
  gDIV(138) = 1.E13
  gDIV(139) = 1.E3
  gDIV(141) = 1.E6
  gDIV(142:145) = 1.E3
  gDIV(146:153) = 1.E-6
  gDIV(154:155) = 1.E3
  gDIV(158:159) = 1.E3
  gDIV(162:165) = 1.E3
  gDIV(166:173) = 1.E3
!  gDIV(184) = 1.E-9
!  gDIV(186) = 1.E-9
!  gDIV(187) = 1.E-9
  gDIV(188) = 1.E-6
  gDIV(196:197) = 1.E3

  !   *** Density perturbation technique ***

  !   Normalized radius where density increase
  DelRho = 0.2D0

  !   Amount of increase in density (10^20 m^-3)
  DelN = 1.D-1

  !   ***

  !   Index for graphic save interval (module tx_graphic)

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
  INTEGER(4) :: NR, NRL, nr_rhoc_near
  real(8)    :: MAXAMP, C1L, C2L, W1L, W2L, CLNEW
  real(8)    :: rhoc, rhol, rhocl

  !  Mesh

  if(rhoaccum < 0.d0) then
     rhoc = 1.d0
  else
     rhoc = rhoaccum
  end if

  !  As a trial, generate mesh using given CL and WL and seek the position in the
  !  original coordinate, which becomes the nearest mesh of separatrix after mapping.
  C1L = CMESH0
  C2L = CMESH
  W1L = WMESH0
  W2L = WMESH
  MAXAMP = LORENTZ(rhob,C1L,C2L,W1L,W2L,0.D0,rhoc) / rhob ! normalized parameter C_0
  nr_rhoc_near = 0
  R(0) = 0.D0
  DO NR = 1, NRMAX - 1
     rhol = NR * rhob / NRMAX
     CALL BISECTION(LORENTZ,C1L,C2L,W1L,W2L,0.D0,rhoc,MAXAMP,rhol,rhob,rho(NR))
     IF(ABS(rho(NR)-rhoc) <= ABS(rho(NR)-rho(nr_rhoc_near))) nr_rhoc_near = NR
  END DO
  rho(NRMAX) = rhob

  !  Construct new CL value that separatrix is just on mesh.
  !  New CL is chosen in order not to be settle so far from given CL.
  !  The mesh finally obtained is well-defined.
  rhocl = nr_rhoc_near * rhob / NRMAX
  CLNEW = ( (rhoc - rhocl) * rhob - rhocl * C1L * LORENTZ_PART(rhob,W1L,W2L,0.D0,rhoc,0) &
       &   + rhob * C1L * LORENTZ_PART(rhoc,W1L,W2L,0.D0,rhoc,0)) &
       &  / (  rhocl * LORENTZ_PART(rhob,W1L,W2L,0.D0,rhoc,1) &
       &     - rhob  * LORENTZ_PART(rhoc,W1L,W2L,0.D0,rhoc,1))
  MAXAMP = LORENTZ(rhob,C1L,CLNEW,W1L,W2L,0.D0,rhoc) / rhob
  rho(0) = 0.D0
  DO NR = 1, NRMAX - 1
     rhol = NR * rhob / NRMAX
     CALL BISECTION(LORENTZ,C1L,CLNEW,W1L,W2L,0.D0,rhoc,MAXAMP,rhol,rhob,rho(NR))
  END DO
  rho(NRMAX) = rhob

  !  Maximum NR till rhoc

  DO NR = 0, NRMAX-1
     IF(rho(NR) <= rhoc .AND. rho(NR+1) >= rhoc) THEN
        NRL = NR
        EXIT
     END IF
  END DO

  !  Adjust rhoc on mesh (for the time being rhoc is assumed to be equivalent to 1.d0)

  IF(ABS(rho(NRL)-rhoc) < ABS(rho(NRL+1)-rhoc)) THEN
     NRA = NRL
     rhol = 0.5d0 * ABS(rho(NRA) - rhoc)
     rho(NRA  ) = rhoc
     rho(NRA-1) = rho(NRA-1) + rhol
  ELSE
     NRA = NRL + 1
     rhol = 0.5d0 * ABS(rho(NRA) - rhoc)
     rho(NRA  ) = rhoc
     rho(NRA+1) = rho(NRA+1) - rhol
  END IF

  !  rhoaccum doesn't coincide with the plasma surface 1.d0

  IF(rhoaccum >= 0.D0) THEN

     !  Maximum NR till rho=1.d0

     DO NR = 0, NRMAX-1
        IF(rho(NR) <= 1.d0 .AND. rho(NR+1) >= 1.d0) THEN
           NRL = NR
           EXIT
        END IF
     END DO

     !  Adjust rho=1.d0 on mesh

     IF(ABS(rho(NRL)-1.d0) < ABS(rho(NRL+1)-1.d0)) THEN
        NRA = NRL
        rhol = 0.5d0 * ABS(rho(NRA) - 1.d0)
        rho(NRA  ) = 1.d0
        rho(NRA-1) = rho(NRA-1) + rhol
     ELSE
        NRA = NRL + 1
        rhol = 0.5d0 * ABS(rho(NRA) - 1.d0)
        rho(NRA  ) = 1.d0
        rho(NRA+1) = rho(NRA+1) - rhol
     END IF
     
  END IF

  !  Mesh number at the center of the core plasma

  DO NR = 0, NRMAX-1
     IF(rho(NR) <= 0.5d0.AND.rho(NR+1) >= 0.5d0) THEN
        IF(ABS(rho(NR)-0.5d0) < ABS(rho(NR+1)-0.5d0)) THEN
           NRC = NR
        ELSE
           NRC = NR + 1
        END IF
        EXIT
     END IF
  END DO

  !  Equilibrium
  call txequ

  RETURN
END SUBROUTINE TXCALM

!***************************************************************
!
!   Initialize profiles
!
!***************************************************************

SUBROUTINE TXPROF

  use tx_commons
  use tx_graphic, only : NGR, NGT, NGVV
  use tx_variables
  use tx_interface, only : INTDERIV3, detect_datatype, dfdx, &
       &                   initprof_input, moving_average, CORR, coulog, inexpolate
  use tx_core_module, only : intg_area, intg_area_p, intg_vol_p 
  use sauter_mod
  use tx_ntv, only : perturb_mag, Wnm_spline
  use eqread_mod, only : AJphVRL
#ifdef laself
  ! for self-compiled lapack
  use f95_lapack, only : GESV => LA_GESV
#else
  ! for intel mkl LAPACK95, 
  !  Note: This module file includes "ptsv" subroutine, 
  !        whose name conflicts with PTsV defined in TASK/TX.  
  use lapack95, only : GESV
#endif

  implicit none
  INTEGER(4) :: NR, IER, ifile, NHFM, NR_smt, NR_smt_start = 10
  REAL(8) :: rhol, PROFN, PROFT, PTePROF, PTiPROF!, RL, QL, dRIP
  REAL(8) :: AJFCT, SUM_INT, DR1, DR2
  REAL(8) :: EpsL, FTL, PBA, dPN, CfN1, CfN2, pea, pia, pediv, pidiv, dpea, dpia, &
       &     Cfpe1, Cfpe2, Cfpi1, Cfpi2, sigma, fexp, PN0L, PNaL, PNeDIVL, &
       &     PTe0L, PTi0L, PTeaL, PTiaL, PTeDIVL, PTiDIVL
  real(8) :: BCLQm3, etanc, etaspz, dum=0.d0, tmp
  REAL(8) :: aitken2p!, DERIV4
  real(8), dimension(:), allocatable :: AJPHL, tmpa, RHSV, Prof1, Prof2, &
       & Profsdt, dProfsdt
  real(8), dimension(:,:), allocatable :: CMTX!, dPsV

  !  Read spline table for neoclassical toroidal viscosity
  IF(FSRP /= 0.D0) CALL Wnm_spline

  NEMAX = NRMAX

  !  Define basic quantities like mass of particles, mesh, etc.

  CALL TXCALM

  !  Contribution of perturbed magnetic field
!  IF(DltRPn /= 0.D0) CALL perturb_mag

  !  Initialize variable array

  X = 0.D0

  !  Variables

  ! ********************************************************************
  !      B.C. and profiles of density, temperature, pressure
  ! ********************************************************************

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
  PBA   = rhob - 1.d0
  dPN   = - 3.D0 * (PN0L - PNaL) / ravl
  CfN1  = - (3.D0 * PBA * dPN + 4.D0 * (PNaL - PNeDIVL)) / PBA**3
  CfN2  =   (2.D0 * PBA * dPN + 3.D0 * (PNaL - PNeDIVL)) / PBA**4
  IF(MDFIXT == 0) THEN
     pea   = PNaL    * PTeaL   ; pia   = PNaL    / achg(2) * PTiaL
     dpea  = dPN     * PTeaL   ; dpia  = dPN     / achg(2) * PTiaL
     pediv = PNeDIVL * PTeDIVL ; pidiv = PNeDIVL / achg(2) * PTiDIVL
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
     rhol = rho(NR)
     ! === Density, temperature, pressure ===
     ! +++ Core +++
     IF (rhol <= 1.d0) THEN
        IF(MDINTN < 0) THEN
           call initprof_input(NR,1,X(NR,LQe1)) ! Ne
           X(NR,LQi1) = X(NR,LQe1) / achg(2)    ! Ni
        ELSE
           PROFN = (1.D0 - rho(NR)**PROFN1)**PROFN2
           X(NR,LQe1) = (PN0L - PNaL) * PROFN + PNaL ! Ne
           X(NR,LQi1) = X(NR,LQe1) / achg(2)         ! Ni
        END IF
        IF(MDINTT < 0) THEN
           call initprof_input(NR,2,PTePROF) ! Te
           call initprof_input(NR,3,PTiPROF) ! Ti
        ELSE IF(MDINTT == 0) THEN
           PROFT = (1.D0 - rho(NR)**PROFT1)**PROFT2
           PTePROF = (PTe0L - PTeaL) * PROFT + PTeaL
           PTiPROF = (PTi0L - PTiaL) * PROFT + PTiaL
        ELSE IF(MDINTT == 1) THEN
           PTePROF = (PTe0L - PTeaL) * (0.8263d0 * (1.d0 - rho(NR)**2 )**1.5d0 &
                &                     + 0.167d0  * (1.d0 - rho(NR)**30)**1.25d0) + PTeaL
           PTiPROF = (PTi0L - PTiaL) * (0.8263d0 * (1.d0 - rho(NR)**2 )**1.5d0 &
                &                     + 0.167d0  * (1.d0 - rho(NR)**30)**1.25d0) + PTiaL
        ELSE
           PTePROF = PTe0L * exp(- log(PTe0L / PTeaL) * rho(NR)**2)
           PTiPROF = PTi0L * exp(- log(PTi0L / PTiaL) * rho(NR)**2)
        END IF
        IF(MDFIXT == 0) THEN
           X(NR,LQe5) = PTePROF * X(NR,LQe1) ! Ne*Te
           X(NR,LQi5) = PTiPROF * X(NR,LQi1) ! Ni*Ti
        ELSE 
           X(NR,LQe5) = PTePROF ! Te
           X(NR,LQi5) = PTiPROF ! Ti
        END IF
     ELSE
        ! +++ SOL +++
        ! Density
        IF(MDITSN == 0) THEN
           X(NR,LQe1) = PNaL + dPN * (rhol - 1.d0) + CfN1 * (rhol - 1.d0)**3 &
                &                                  + CfN2 * (rhol - 1.d0)**4
           X(NR,LQi1) = X(NR,LQe1) / achg(2)
        ELSE
           X(NR,LQe1) = PNaL * EXP(- (rhol - 1.d0) / rLn)
           X(NR,LQi1) = X(NR,LQe1) / achg(2)
        END IF
        IF(MDITST == 0) THEN
           X(NR,LQe5) = pea + dpea * (rhol - 1.d0) + Cfpe1 * (rhol - 1.d0)**3 &
                &                                  + Cfpe2 * (rhol - 1.d0)**4
           X(NR,LQi5) = pia + dpia * (rhol - 1.d0) + Cfpi1 * (rhol - 1.d0)**3 &
                &                                  + Cfpi2 * (rhol - 1.d0)**4
        ELSE
           IF(MDITST == 1) THEN
              PTePROF = PTeaL * EXP(- (rhol - 1.d0) / rLT)
              PTiPROF = PTiaL * EXP(- (rhol - 1.d0) / rLT)
           ELSE
              sigma = 0.3d0 * (rho(NRMAX) - 1.d0)
              fexp  = exp(- (rho(NR) - 1.d0)**2 / (2.d0 * sigma**2))
              PTePROF = PTe0L   * exp(- log(PTe0L / PTeaL) * rho(NR)**2) *         fexp &
                   &  + PTeDIVL                                          * (1.d0 - fexp)
              PTiPROF = PTi0L   * exp(- log(PTi0L / PTiaL) * rho(NR)**2) *         fexp &
                   &  + PTiDIVL                                          * (1.d0 - fexp)
           END IF
           ! Pressure (MDFIXT=0) or temperature (MDFIXT=1)
           IF(MDFIXT == 0) THEN
              X(NR,LQe5) = PTePROF * X(NR,LQe1)
              X(NR,LQi5) = PTiPROF * X(NR,LQi1)
           ELSE
              X(NR,LQe5) = PTePROF
              X(NR,LQi5) = PTiPROF
           END IF
        END IF
     END IF

     ! === Neutrals ===
     ! N0_1 (slow neutrals)
     X(NR,LQn1) = PN0s
     ! N0_2 (thermal neutrals)
!     X(NR,LQn2) = 0.D0
     X(NR,LQn2) = 1.D-20 ! when ThntSW = 0.D0
     ! N0_3 (halo neutrals)
     X(NR,LQn3) = 0.D0

     ! === Fixed densities to keep them constant during iterations ===
     PNsV_FIX(NR,1) = X(NR,LQe1)
     PNsV_FIX(NR,2) = X(NR,LQi1)
     IF(MDFIXT == 0) THEN
        PTsV_FIX(NR,1) = X(NR,LQe5) / X(NR,LQe1)
        PTsV_FIX(NR,2) = X(NR,LQi5) / X(NR,LQi1)
     ELSE
        PTsV_FIX(NR,1) = X(NR,LQe5)
        PTsV_FIX(NR,2) = X(NR,LQi5)
     END IF
  END DO

  ! === Smoothing profiles ===
  IF(MDINTT == -2) THEN ! Smoothing temperatures
     allocate(Prof1(0:NRMAX),Prof2(0:NRMAX))
     Prof1(:) = X(:,LQe5) / X(:,LQe1)
     Prof2(:) = X(:,LQi5) / X(:,LQi1)
     NR_smt = NRA - NR_smt_start ! smoothing data only in the edge region
     DO NR = NR_smt, NRMAX
        X(NR,LQe5) = moving_average(NR,Prof1,NRMAX) * X(NR,LQe1)
        X(NR,LQi5) = moving_average(NR,Prof2,NRMAX) * X(NR,LQi1)
     END DO
     deallocate(Prof1,Prof2)
  END IF
  IF(MDINTN == -2) THEN ! Smoothing densities
     allocate(Prof1(0:NRMAX))
     X(:,LQe5) = X(:,LQe5) / X(:,LQe1)
     X(:,LQi5) = X(:,LQi5) / X(:,LQi1)
     Prof1(:) = X(:,LQe1)
     NR_smt = NRA - NR_smt_start ! smoothing data only in the edge region
     DO NR = NR_smt, NRMAX
        X(NR,LQe1) = moving_average(NR,Prof1,NRMAX)
        X(NR,LQi1) = X(NR,LQe1) / achg(2)       ! Ni
     END DO
     X(:,LQe5) = X(:,LQe5) * X(:,LQe1)
     X(:,LQi5) = X(:,LQi5) * X(:,LQi1)
     deallocate(Prof1)
  END IF

  ! ********************************************************************
  !      Poloidal current function, fipol = R B_t
  !      Poloidal magnetic field, BthV
  !      dpsi/dV, sdt
  ! ********************************************************************

  if( ieqread >= 2 ) then
     ! PsitV, fipol and sdt have already been determined in intequ.

     X(:,LQm5) = PsitV(:) / rMU0 ! PsitV / rMU0

     BthV(0)       = 0.d0
     BthV(1:NRMAX) = sqrt(ckt(1:NRMAX)) * sdt(1:NRMAX)
  else
     sum_int = 0.d0
     fipol(:) = rbvt ! flat fipol assumed at initial
     do NR = 0, NRMAX
        sum_int = sum_int + intg_vol_p(aat,nr) * fipol(NRMAX) / (4.d0 * Pisq)
        X(NR,LQm5) = sum_int / rMU0 ! PsitV / rMU0
     end do

     ! === Poloidal magnetic field ===

     allocate(Profsdt(0:NRMAX),dProfsdt(0:NRMAX))
     BCLQm3 = 2.d0 * Pi * rMUb1 * rIps * 1.d6
     !  dPsi/dV
     do NR = 1, NRMAX
        if(rho(nr) < 1.d0) then ! Core
           ! Profsdt is used for AJPHL. dProfsdt = d Profsdt/d (vv/vv(NRA))
           Profsdt(NR) = 1.d0 - (1.d0 - vv(NR)/vv(NRA))**(PROFJ+1)
           dProfsdt(NR) = (PROFJ+1) * (1.d0 - vv(NR)/vv(NRA))**PROFJ
           sdt(NR)  = BCLQm3 / ckt(NR) * Profsdt(NR)
        else ! SOL, no current assumption
           Profsdt(NR) = 1.d0
           dProfsdt(NR) = 0.d0
           sdt(NR)  = BCLQm3 / ckt(NR)
        end if
        BthV(NR) = sqrt(ckt(NR)) * sdt(NR)
     end do
     BthV(0) = 0.d0
     sdt(0)  = aitken2p(0.d0,sdt(1),sdt(2),sdt(3),vv(1),vv(2),vv(3))
     Profsdt(0) = 1.d0
     dProfsdt(0) = PROFJ+1
  end if

!  IF(FSHL /= 0.D0) THEN
!     BthV(0) = 0.D0
!     Q(0) = Q0
!     DO NR = 1, NRMAX
!        RL = R(NR)
!        Q(NR) = (Q0 - QA) * (1.D0 - (RL / RA)**2) + QA
!        BthV(NR) = BB * RL / (Q(NR) * RR)
!     END DO
!  END IF

  ! ********************************************************************
  !      Toroidal electron current 
  !
  !   electron current: AJPHL   = - e ne <ueph/R> / <1/R>
  !                   : AJphVRL = - e ne <ueph/R>
  !
  !      dpsi/dV, sdt (in some case)
  ! ********************************************************************

  if( ieqread >= 2 ) then
     ! AJphVRL have already been allocated and determined in intequ.

     DO NR = 0, NRMAX
        X(NR,LQe7) =-AJphVRL(NR) / (AEE * 1.D20)
        X(NR,LQe4) = X(NR,LQe7) / aat(NR) ! approx
        AJOH(NR)   = X(NR,LQe7) / ait(NR)
        X(NR,LQm4) = PsiV(NR)
     END DO

  else
     
     ! === File input ===

     allocate(AJPHL(0:NRMAX),AJphVRL(0:NRMAX))
     ifile = detect_datatype('LQe4')
     if(ifile == 0) then ! No toroidal current (LQe4) data in a structured type
        IF(MDINTC <= -1) THEN ! Current density read from file
           DO NR = 0, NRA
              call initprof_input(NR,4,AJPHL(NR)) ! electron current read from an external file
           END DO
           AJPHL(NRA+1:NRMAX) = 0.D0
           AJFCT = rIPs * 1.D6 / intg_area(AJPHL)
           ! Artificially extrapolate a current density in the SOL for numerical stability
           DO NR = NRA+1, NRMAX
              AJPHL(NR) = AJPHL(NRA) * EXP(- (rho(NR) - 1.d0) / (0.5d0 * rLn))
           END DO

           IF(MDINTC == -2) THEN ! Smoothing current density
              allocate(Prof1(0:NRMAX))
              Prof1(:) = AJPHL(:)
              NR_smt = NRA - NR_smt_start ! smoothing data only in the edge region
              DO NR = NR_smt, NRMAX
                 AJPHL(NR) = moving_average(NR,Prof1,NRMAX)
              END DO
              deallocate(Prof1)
           END IF

           DO NR = 0, NRMAX
              AJPHL(NR)   = AJFCT * AJPHL(NR)
              AJphVRL(NR) = AJPHL(NR) * ait(NR)
              X(NR,LQe4)  =-AJphVRL(NR) / (AEE * 1.D20) * rrt(NR) ! approx
              X(NR,LQe7)  =-AJphVRL(NR) / (AEE * 1.D20)
              AJOH(NR)    = AJPHL(NR)
           END DO

           BthV(0) = 0.d0
           sum_int = 0.d0
           DO NR = 1, NRMAX
              sum_int  = sum_int + intg_vol_p(AJphVRL,NR)
              BthV(NR) = rMUb1 * sum_int / sqrt(ckt(NR))
              sdt(NR)  = rMUb1 * sum_int / ckt(NR)
           END DO
           sdt(0)  = aitken2p(0.d0,sdt(1),sdt(2),sdt(3),vv(1),vv(2),vv(3))
        ELSE ! (MDINTC == 0); Current density constructed
           DO NR = 0, NRMAX
              AJphVRL(NR) = BCLQm3 / (rMUb1 * vlt(nra)) * dProfsdt(NR) ! <j_zeta/R>
              AJPHL(NR)   = AJphVRL(NR) / ait(NR) ! <j_zeta/R>/<1/R>
              X(NR,LQi7)  = (Uiph0 * ait(NR)) * X(NR,LQi1) * (dProfsdt(NR) / dProfsdt(0))
              X(NR,LQi4)  = X(NR,LQi7) * rrt(NR) ! n <R u_zeta> = n <u_zeta/R>/<R^2>
              X(NR,LQe7)  =-AJphVRL(NR) / (AEE * 1.D20) + achg(2) * X(NR,LQi7)
              X(NR,LQe4)  = X(NR,LQe7) * rrt(NR) ! n <R u_zeta> = n <u_zeta/R>/<R^2>
              AJOH(NR)    = AJPHL(NR)
!              IF(FSHL == 0.D0) THEN
!                 AJphVRL(NR) = 0.D0
!                 X(NR,LQe4)  = 0.D0
!                 X(NR,LQe7)  = 0.D0
!                 AJOH(NR)    = 0.D0
!              END IF
           END DO
        END IF

     else ! Detect toroidal current data in a structured type
        call inexpolate(infiles(ifile)%nol,infiles(ifile)%r,infiles(ifile)%data,NRMAX,RHO,5,AJPHL)
        AJFCT = rIPs * 1.D6 / intg_area(AJPHL)
        DO NR = 0, NRMAX
           AJPHL(NR)   = AJFCT * AJPHL(NR)
           AJphVRL(NR) = AJPHL(NR) * ait(NR)
           X(NR,LQe4)  =-AJPHL(NR) / (AEE * 1.D20) * rrt(NR) ! approx
           X(NR,LQe7)  =-AJPHL(NR) / (AEE * 1.D20)
           AJOH(NR)    = AJPHL(NR)
        END DO

        BthV(0) = 0.d0
        sum_int = 0.d0
        DO NR = 1, NRMAX
           sum_int  = sum_int + intg_vol_p(AJphVRL,NR)
           BthV(NR) = rMUb1 * sum_int / sqrt(ckt(NR))
           sdt(NR)  = rMUb1 * sum_int / ckt(NR)
        END DO
        sdt(0)  = aitken2p(0.d0,sdt(1),sdt(2),sdt(3),vv(1),vv(2),vv(3))
     end if

     deallocate(Profsdt,dProfsdt,AJPHL)

     sum_int = 0.d0
     X(0,LQm4) = 0.d0
     do NR = 1, NRMAX
        sum_int = sum_int + intg_vol_p(sdt,nr)
        X(NR,LQm4) = sum_int
     end do
  end if

  ! deallocate arrays in initprof_input
  if(MDINTN < 0 .or. MDINTT < 0 .or. ABS(MDINTC) /= 0) call initprof_input(idx = 0)

  ! ********************************************************************
  !      Electrostatic potential, Phi
  ! ********************************************************************

  ! Inverse matrix of derivative formula for integration

  allocate(CMTX(1:NRMAX,1:NRMAX),RHSV(1:NRMAX))
  CMTX(:,:) = 0.D0
  NR = 1
     DR1 = vv(NR-1) - vv(NR)
     DR2 = vv(NR+1) - vv(NR)
     CMTX(NR,NR  ) = - (DR1 + DR2) / (DR1 * DR2)
     CMTX(NR,NR+1) = - DR1**2 / (DR1 * DR2 * (DR2 - DR1))
  DO NR = 2, NRMAX-1
     DR1 = vv(NR-1) - vv(NR)
     DR2 = vv(NR+1) - vv(NR)
     CMTX(NR,NR-1) =   DR2**2 / (DR1 * DR2 * (DR2 - DR1))
     CMTX(NR,NR  ) = - (DR1 + DR2) / (DR1 * DR2)
     CMTX(NR,NR+1) = - DR1**2 / (DR1 * DR2 * (DR2 - DR1))
  END DO
  NR = NRMAX
     DR1 = vv(NR-1) - vv(NR)
     DR2 = vv(NR-2) - vv(NR)
     CMTX(NR,NR-2) = - DR1**2 / (DR1 * DR2 * (DR2 - DR1))
     CMTX(NR,NR-1) =   DR2**2 / (DR1 * DR2 * (DR2 - DR1))
     CMTX(NR,NR  ) = - (DR1 + DR2) / (DR1 * DR2)

  ! Numerical solution for LQm1: integrate Phi' = - p_i'/(Z_i e n_i) over V

  allocate(tmpa(0:NRMAX))
  if(MDFIXT == 0) then
     tmpa(:) = dfdx(vv,X(:,LQi5),NRMAX,0) * rKilo ! p'
  else
     tmpa(:) = dfdx(vv,X(:,LQi1)*X(:,LQi5),NRMAX,0) * rKilo ! p'
  end if

  RHSV(1:NRMAX) = - tmpa(0:NRMAX-1) / (X(0:NRMAX-1,LQi1) * achg(2))
  if( MDLPCK /= 0 ) then
     call gesv(CMTX,RHSV)
  else
     CALL INVMRD(CMTX,NRMAX,NRMAX,IER)
     RHSV(1:NRMAX) = matmul(CMTX,RHSV)
  end if

  do NR = NRMAX, 1, -1
     X(NR,LQm1) = RHSV(NR) - RHSV(NRMAX)
  end do
  X(0,LQm1) = - RHSV(NRMAX)
  deallocate(tmpa,CMTX,RHSV)

  ! ********************************************************************
  !      Parallel flows
  !      Safety factor, Q
  ! ********************************************************************

  ! === Parallel flow ===
  if( ieqread >= 2 ) then
     X(:,LQe3) = X(:,LQe4) / X(:,LQe1) * (bbt(:) / fipol(:))
     X(:,LQi3) = X(:,LQi4) / X(:,LQi1) * (bbt(:) / fipol(:))
  else
     do NR = 0, NRMAX
        bbt(NR) = bb * bb
        tmp = bbt(NR) / fipol(NR)
        X(NR,LQe3) = X(NR,LQe4) / X(NR,LQe1) * tmp
        X(NR,LQi3) = X(NR,LQi4) / X(NR,LQi1) * tmp
     end do
  end if

  ! === Safety factor ===

  Q(:) = fipol(:) * aat(:) / (4.d0 * Pisq * sdt(:))

  ! === Poloidal current density (Virtual current for helical system) ===

  AJV(:)=0.D0

!  IF(FSHL == 0.D0) THEN
!     ! Integrate 1 / (r * rMU0) * d/dr (r * BthV) to obtain AJV
!     AJV(:) = BB / (RR * rMU0) * 2.d0 * Q0 / Q(:)**2
!  END IF

  ! ********************************************************************
  !      Toroidal electric field for initial NCLASS calculation
  ! ********************************************************************

  DO NR = 0, NRMAX
     Var(NR,1)%n = X(NR,LQe1)
     Var(NR,2)%n = X(NR,LQi1)
     IF(MDFIXT == 0) THEN
        Var(NR,1)%T = X(NR,LQe5) / X(NR,LQe1)
        Var(NR,2)%T = X(NR,LQi5) / X(NR,LQi1)
     ELSE
        Var(NR,1)%T = X(NR,LQe5)
        Var(NR,2)%T = X(NR,LQi5)
     END IF
     ! Inverse aspect ratio
     EpsL = epst(NR)
     ! Trapped particle fraction
     FTL  = 1.46D0 * SQRT(EpsL) - 0.46D0 * EpsL**1.5D0
     ! Estimating parallel resistivity
     call sauter(Var(NR,1)%n,Var(NR,1)%T,dum,dum,Var(NR,2)%n,Var(NR,2)%T,dum,dum, &
       &         Q(NR),sdt(NR),fipol(NR),EpsL,RR,achg(2),Zeff,ftl, &
       &         rlnLei_IN=coulog(1.d0,Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amas(1),achg(1),amas(2),achg(2)), &
       &         rlnLii_IN=coulog(1.d0,Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amas(2),achg(2),amas(2),achg(2)), &
       &         BJBS=dum,ETA=etanc,ETAS=etaspz)
     if(abs(FSNC) > 0.d0) then
        eta(NR) = etanc
     else
        eta(NR) = etaspz
     end if
     ! R Et = PsitdotV
     !    <B E//> = eta// <B j//> and neglecting poloidal components yields <Bt Et> = eta// <Bt j//>
     !    <Bt Et> = I<1/R^2> R E_t and <Bt jt> = I<jt/R> give R Et = eta// <jt/R> / <1/R^2>
     X(NR,LQm3) = eta(NR) * AJphVRL(NR) / aat(NR)
     IF(X(NR,LQm3) == 0.D0) X(NR,LQm3) = 1.D-4
  END DO
!  IF(FSHL == 0.D0) X(:,LQm3) = 0.D0
  if( allocated(AJphVRL) ) deallocate(AJphVRL)

  ! ********************************************************************
  !      Miscellaneous
  ! ********************************************************************

  !   NBI total input power (MW)
  PNBH = PNBHP + PNBHT1 + PNBHT2

  !     Averaged beam injection energy
  Eb =  (esps(1) + 0.5d0 * esps(2) + esps(3) / 3.d0) * Ebmax

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
  ErV_FIX(:) = ErV(:)

  !  Calculate various physical quantities

  CALL TXCALC(0)

  !  Calculate global quantities for storing and showing initial status

  CALL TXGLOB

  RETURN
END SUBROUTINE TXPROF

!*****************************************************************************************

module tx_parameter_control
  use tx_commons
  use tx_graphic, only : MODEG, MODEGL, NGRSTP, NGTSTP, NGVSTP, gDIV
  implicit none
  public
  NAMELIST /TX/ &
       & RA,rhob,rhoaccum,RR,BB,rbvt, &
       & amas,achg,Zeff, &
       & PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ,PROFN1,PROFN2,PROFT1,PROFT2,Uiph0, &
       & De0,Di0,VWpch0,rMue0,rMui0,WPM0, &
       & Chie0,Chii0,ChiNC, &
       & FSDFIX,FSANOM,FSCBAL,FSCBKP,FSCBEL,FSCBSH,FSBOHM,FSPCLD,FSPCLM,FSPCLC,FSVAHL,MDANOM,RhoETB, &
       & FSMPCH,FSPARV,PROFD,PROFD1,PROFD2,PROFDB,PROFM,PROFM1,PROFMB,PROFC,PROFC1,PROFCB, &
       & FSCX,FSLC,FSRP,FSNF,FSADV,FSADVB,FSUG,FSNC,FSNCB,FSLP,FSLTE,FSLTI,FSION,FSD01,FSD02,FSD03,MDLC, &
       & rLn,rLT, &
       & Ebmax,esps,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBHP,PNBHT1,PNBHT2,PNBCD,PNBMPD,PNBPTC, &
       & rNRFe,RRFew,RRFe0,PRFHe,Tqt0,Tqp0,rNRFi,RRFiw,RRFi0,PRFHi, &
       & PN0s,V0,rGamm0,rGASPF,PNeDIV,PNiDIV,PTeDIV,PTiDIV, &
       & NTCOIL,DltRPn,kappa,m_pol,n_tor, &
       & DT,EPS,ICMAX,ADV,tiny_cap,CMESH0,CMESH,WMESH0,WMESH, &
       & NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,ieqread, &
       & DelRho,DelN, &
!       & DMAG0,RMAGMN,RMAGMX,EpsH,NCph,NCth, &
       & DMAG0,RMAGMN,RMAGMX, &
       & rG1,FSHL,Q0,QA, &
       & rIPs,rIPe, &
       & MODEG,gDIV,MODEAV,MODEGL,MDLPCK,MODECV,oldmix,iSUPG2,iSUPG3,iSUPG6,SUPGstab, &
       & MDFIXT,MDBEAM,MDOSQZ,MDOSQZN,MDLETA,MDLNEO,MDBSETA,MDITSN,MDITST,MDINTN,MDINTT,MDINTC,MDLETB, &
       & IDIAG,IGBDF,ISMTHD,MDLNBD, & ! 09/06/17~ miki_m
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
       WRITE(6,*) '  INPUT PARAMETERS MUST BE BRACKETED BETWEEN &TX AND &END; CTRL-D FOR EXIT.'
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
       IF(rhob <= 1.d0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rhoaccum > rhob) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(RR <= rhob) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rIPs < 0.D0 .OR. rIPe < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(amas(1) < 0.D0 .OR. amas(2) < 0.D0 .OR. achg(1) > 0.D0 .OR. achg(2) < 0.D0 .OR. Zeff < 1.D0) THEN
          EXIT ; ELSE ; idx = idx + 1 ; ENDIF
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
       IF(minval(FSMPCH) < 0.d0 .or. maxval(FSMPCH) > 1.d0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(minval(FSPARV) < 0.d0 .or. maxval(FSPARV) > 1.d0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(FSCBAL < 0.D0 .OR. FSCBKP < 0.D0 .OR. FSCBEL < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(FSBOHM < 0.D0 .OR. FSPCLD < 0.D0 .OR. FSPCLM < 0.D0 .OR. FSPCLC < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       ! /// idx = 21 - 30 ///
       IF(FSLC < 0.D0 .OR. FSRP < 0.D0 .OR. FSNC < 0.D0 .OR. FSNCB < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(FSHL < 0.D0 .OR. FSNF < 0.D0 .OR. FSADV < 0.d0 .OR. FSADVB < 0.d0 .OR. FSUG < 0.d0) THEN
          EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(FSLP < 0.D0 .OR. FSLTE < 0.D0 .OR. FSLTI < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(FSION < 0.D0)  THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(FSD01 < 0.D0 .OR. FSD02 < 0.D0 .OR. FSD03 < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(MDLC /= 1 .AND. MDLC /= 2) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rG1 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(Ebmax < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(sum(esps) /= 1.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(PNBHP < 0.D0 .OR. PNBHT1 < 0.D0 .OR. PNBHT2 < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       ! /// idx = 31 - 40 ///
       IF(ABS(PNBCD) > 1.D0 .OR. ABS(PNBMPD) > 1.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(PNBPTC < 0.D0 .OR. PNBPTC > 1.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(RNBP0 > rhob*ra .OR. RNBP0 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(RNBT10 > rhob*ra .OR. RNBT10 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(RNBT20 > rhob*ra .OR. RNBT20 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rNRFe < 0.D0 .OR. PRFHe < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rNRFi < 0.D0 .OR. PRFHi < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(RRFe0 > rhob*ra .OR. RRFe0 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(RRFi0 > rhob*ra .OR. RRFi0 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(PN0s < 0.D0 .OR. V0 < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(rGamm0 < 0.D0 .OR. rGamm0 > 1.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       ! /// idx = 41 - 50 ///
       IF(rGASPF < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(PNeDIV < 0.D0 .OR. PNiDIV < 0.D0 .OR. PTeDIV < 0.D0 .OR. PTiDIV < 0.D0) THEN ; EXIT ; ELSE
          idx = idx + 1 ; ENDIF
       IF(NTCOIL <= 0 .OR. DltRPn < 0.D0 .OR. DltRPn > 1.D0 .OR. kappa < 0.D0) THEN
          EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(DT < 0.D0 .OR. EPS == 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(ICMAX < 0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(ADV < 0.D0 .OR. ADV > 1.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(tiny_cap < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(oldmix < 0.D0 .OR. oldmix > 1.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(CMESH0 < 0.D0 .OR. CMESH < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(WMESH0 < 0.D0 .OR. WMESH < 0.D0) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       ! /// idx = 51 - 52 ///
       IF(MDLNEO /= 1 .AND. MDLNEO /= 2) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       IF(MDBEAM /= 0 .AND. MDBEAM /= 1) THEN ; EXIT ; ELSE ; idx = idx + 1 ; ENDIF
       RETURN
    END DO

    WRITE(6,'(A,I3)') 'XX INPUT ERROR: Please check consistency and/or sanity of INPUT PARAMETERS. idx =',idx
    STOP
    
  END SUBROUTINE TXPARM_CHECK

  !***** INPUT PARAMETER LIST *****

  SUBROUTINE TXPLST

    WRITE(6,601)
    RETURN

601 FORMAT(' ','# &TX : RA,rhob,rhoaccum,RR,BB,rbvt,amas,achg,Zeff,'/ &
         &       ' ',8X,'PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ,,PROFN1,PROFN2,PROFT1,PROFT2,Uiph0,'/ &
         &       ' ',8X,'De0,Di0,VWpch0,rMue0,rMui0,WPM0,'/ &
         &       ' ',8X,'Chie0,Chii0,ChiNC,'/ &
         &       ' ',8X,'FSDFIX,FSANOM,FSCBAL,FSCBKP,FSCBEL,FSCBSH,FSBOHM,FSPCLD,FSPCLM,FSPCLC,FSVAHL,'/ &
         &       ' ',8X,'MDANOM,RhoETB,'/ &
         &       ' ',8X,'PROFD,PROFD1,PROFD2,PROFDB,PROFM,PROFM1,PROFMB,PROFC,PROFC1,PROFCB,'/ &
         &       ' ',8X,'FSCX,FSLC,FSRP,FSNF,FSADV,FSADVB,FSUG,FSNC,FSNCB,FSLP,FSLTE,FSLTI,FSION,FSD01,FSD02,FSD03,'/&
         &       ' ',8X,'MDLC,rLn,rLT,'/ &
         &       ' ',8X,'Ebmax,esps,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBHP,PNBHT1,PNBHT2,'/ &
         &       ' ',8X,'PNBCD,PNBMPD,PNBPTC,rNRFe,RRFew,RRFe0,PRFHe,Tqt0,Tqp0,rNRFe,RRFew,RRFe0,PRFHe,'/&
         &       ' ',8X,'PN0s,V0,rGamm0,rGASPF,PNeDIV,PNiDIV,PTeDIV,PTiDIV,'/ &
         &       ' ',8X,'NTCOIL,DltRPn,kappa,m_pol,n_tor,'/ &
         &       ' ',8X,'DT,EPS,ICMAX,ADV,tiny_cap,CMESH0,CMESH,WMESH0,WMESH,'/ &
         &       ' ',8X,'NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,ieqread,'/ &
         &       ' ',8X,'DelRho,DelN,'/ &
!         &       ' ',8X,'Dmag0,RMAGMN,RMAGMX,EpsH,NCph,NCth,'/ & 
         &       ' ',8X,'Dmag0,RMAGMN,RMAGMX,EpsHM,HPN,'/ &   ! 10/08/06 miki_m
         &       ' ',8X,'rG1,FSHL,Q0,QA,'/ &
         &       ' ',8X,'rIPs,rIPe,'/ &
         &       ' ',8X,'MODEG,gDIV,MODEAV,MODEGL,MDLPCK,MODECV,oldmix,iSUPG2,iSUPG3,iSUPG6,SUPGstab,'/ &
         &       ' ',8X,'MDFIXT,MDBEAM,MDOSQZ,MDOSQZN,MDLETA,MDLNEO,MDBSETA,MDITSN,MDITST,MDINTN,MDINTT,MDLETB,' / & 
         &       ' ',8X,'IDIAG,IGBDF,ISMTHD,MDLNBD')
  END SUBROUTINE TXPLST

!***************************************************************
!
!   View input parameters
!
!***************************************************************

  SUBROUTINE TXVIEW

    WRITE(6,'((1X,A10," =",1PD11.4,2(2X,A10," =",1PD11.4)))') &
         &   'RA        ', RA       , 'RHOB      ', RHOB     ,  &
         &   'RR        ', RR       , 'BB        ', BB       ,  &
         &   'RAVL      ', RAVL     , 'RBVL      ', RBVL     ,  &
         &   'rbvt      ', rbvt     , &
         &   'amas(2)   ', amas(2)  , 'achg(2)   ', achg(2)  ,  &
         &   'PN0       ', PN0      , 'PNa       ', PNa      ,  &
         &   'PTe0      ', PTe0     , 'PTea      ', PTea     ,  &
         &   'PTi0      ', PTi0     , 'PTia      ', PTia     ,  &
         &   'Zeff      ', Zeff     , 'rIP       ', rIP      ,  &
         &   'rIPs      ', rIPs     , 'rIPe      ', rIPe     ,  &
         &   'PROFJ     ', PROFJ    , 'PROFN1    ', PROFN1   ,  &
         &   'PROFN2    ', PROFN2   , 'PROFT1    ', PROFT1   ,  &
         &   'PROFT2    ', PROFT2   , 'Uiph0     ', Uiph0    ,  &
         &   'CMESH0    ', CMESH0   , &
         &   'WMESH0    ', WMESH0   , 'CMESH     ', CMESH    ,  &
         &   'WMESH     ', WMESH    , 'ADV       ', ADV      ,  &
         &   'De0       ', De0      , 'Di0       ', Di0      ,  &
         &   'rMue0     ', rMue0    , 'rMui0     ', rMui0    ,  &
         &   'FSMPCH(1) ', FSMPCH(1), 'FSMPCH(2) ', FSMPCH(2),  &
         &   'FSPARV(1) ', FSPARV(1), 'FSPARV(2) ', FSPARV(2),  &
         &   'VWpch0    ', VWpch0   , 'WPM0      ', WPM0     ,  &
         &   'PROFD     ', PROFD    , 'PROFD1    ', PROFD1   ,  &
         &   'PROFD2    ', PROFD2   , 'PROFDB    ', PROFDB   ,  &
         &   'PROFM     ', PROFM    , 'PROFM1    ', PROFM1   ,  &
         &   'PROFMB    ', PROFMB   , 'PROFC     ', PROFC    ,  &
         &   'PROFC1    ', PROFC1   , 'PROFCB    ', PROFCB   ,  &
         &   'ChiNC     ', ChiNC    , &
         &   'Chie0     ', Chie0    , 'Chii0     ', Chii0    ,  &
         &   'FSDFIX(1) ', FSDFIX(1), 'FSDFIX(2) ', FSDFIX(2),  &
         &   'FSDFIX(3) ', FSDFIX(3), 'FSANOM(1) ', FSANOM(1),  &
         &   'FSANOM(2) ', FSANOM(2), 'FSANOM(3) ', FSANOM(3),  &
         &   'RhoETB(1) ', RhoETB(1), 'RhoETB(2) ', RhoETB(2),  &
         &   'RhoETB(3) ', RhoETB(3),  &
         &   'FSCBAL    ', FSCBAL   , 'FSCBKP    ', FSCBKP   ,  &
         &   'FSCBEL    ', FSCBEL   , 'FSCBSH    ', FSCBSH   ,  &
         &   'FSBOHM    ', FSBOHM   , 'FSPCLD    ', FSPCLD   ,  &
         &   'FSPCLM    ', FSPCLM   , 'FSPCLC    ', FSPCLC   ,  &
         &   'FSVAHL    ', FSVAHL   , 'FSCX      ', FSCX     ,  &
         &   'FSLC      ', FSLC     , 'FSRP      ', FSRP     ,  &
         &   'FSNF      ', FSNF     , 'FSNC      ', FSNC     ,  &
         &   'FSADV     ', FSADV    , 'FSADVB    ', FSADVB   ,  &
         &   'FSUG      ', FSUG     ,  &
         &   'FSLP      ', FSLP     , 'FSLTE     ', FSLTE    ,  &
         &   'FSLTI     ', FSLTI    , 'FSION     ', FSION    ,  &
         &   'FSD01     ', FSD01    , 'FSD02     ', FSD02    ,  &
         &   'FSD03     ', FSD03    ,  &
         &   'rLn       ', rLn      , 'rLT       ', rLT      ,  &
         &   'Ebmax     ', Ebmax    , 'esps(1)   ', esps(1)  ,  &
         &   'esps(2)   ', esps(2)  , 'esps(3)   ', esps(3)  ,  &
         &   'FSNCB     ', FSNCB    ,  &
         &   'RNBP      ', RNBP     , 'RNBP0     ', RNBP0    ,  &
         &   'RNBT1     ', RNBT1    , 'RNBT10    ', RNBT10   ,  &
         &   'RNBT2     ', RNBT2    , 'RNBT20    ', RNBT20   ,  &
         &   'PNBHP     ', PNBHP    , 'PNBHT1    ', PNBHT1   ,  &
         &   'PNBHT2    ', PNBHT2   , 'PNBHex    ', PNBHex   ,  &
         &   'PNBMPD    ', PNBMPD   , 'PNBPTC    ', PNBPTC   ,  &
         &   'rNRFe     ', rNRFe    , 'RRFew     ', RRFew    ,  &
         &   'RRFe0     ', RRFe0    , 'PRFHe     ', PRFHe    ,  &
         &   'rNRFi     ', rNRFe    , 'RRFiw     ', RRFiw    ,  &
         &   'RRFi0     ', RRFi0    , 'PRFHi     ', PRFHi    ,  &
         &   'Tqt0      ', Tqt0     , 'Tqp0      ', Tqp0     ,  &
         &   'rGamm0    ', rGamm0   , 'V0        ', V0       ,  &
         &   'rGASPF    ', rGASPF   ,  &
         &   'PNeDIV    ', PNeDIV   , 'PNiDIV    ', PNiDIV   ,  &
         &   'PTeDIV    ', PTeDIV   , 'PTiDIV    ', PTiDIV   ,  &
         &   'DltRPn    ', DltRPn   , 'kappa     ', kappa    ,  &
         &   'PN0s      ', PN0s     , 'EPS       ', EPS      ,  &
         &   'tiny_cap  ', tiny_cap , 'DT        ', DT       ,  &
         &   'rG1       ', rG1      ,  &
         &   'DMAG0     ', DMAG0    , 'RMAGMN    ', RMAGMN   ,  &
         &   'RMAGMX    ', RMAGMX   ,  &
         &   'FSHL      ', FSHL     ,  &  ! Too many elements of EpsHM to show miki_m 10-08-11
         &   'Q0        ', Q0       , 'QA        ', QA       ,  &
!!$         &   'EpsHM(1,:)', EpsHM(1,0:3) , &
!!$         &   'EpsHM(2,:)', EpsHM(2,0:3) , &
!!$         &   'EpsHM(3,:)', EpsHM(3,0:3) , &
         &   'SUPGstb   ', SUPGstab,  'oldmix    ', oldmix   , &
         &   'EpsHM(1,0)', EpsHM(1,0),'EpsHM(1,1)', EpsHM(1,1), &
         &   'EpsHM(1,2)', EpsHM(1,2),'EpsHM(1,3)', EpsHM(1,3), &
         &   'EpsHM(2,0)', EpsHM(2,0),'EpsHM(2,1)', EpsHM(2,1), &
         &   'EpsHM(2,2)', EpsHM(2,2),'EpsHM(2,3)', EpsHM(2,3), &
         &   'EpsHM(3,0)', EpsHM(3,0),'EpsHM(3,1)', EpsHM(3,1), &
         &   'EpsHM(3,2)', EpsHM(3,2),'EpsHM(3,3)', EpsHM(3,3), &
         &   'EpsHM(4,0)', EpsHM(4,0),'EpsHM(4,1)', EpsHM(4,1), &
         &   'EpsHM(4,2)', EpsHM(4,2),'EpsHM(4,3)', EpsHM(4,3)
    WRITE(6,'((" ",A10," =",I5,2(8X,A10," =",I5)))') &
         &   'NRMAX     ', NRMAX    , &
         &   'NTMAX     ', NTMAX    , 'NTSTEP    ', NTSTEP   ,  &
         &   'NGRSTP    ', NGRSTP   , 'NGTSTP    ', NGTSTP   ,  &
         &   'NGVSTP    ', NGVSTP   , 'ieqread   ', ieqread  ,  &
         &   'ICMAX     ', ICMAX    , 'MODEG     ', MODEG    ,  &
         &   'MODEAV    ', MODEAV   , 'MODEGL    ', MODEGL   ,  &
         &   'MDLPCK    ', MDLPCK   , 'MODECV    ', MODECV   ,  &
         &   'iSUPG2    ', iSUPG2   , 'iSUPG3    ', iSUPG3   ,  &
         &   'iSUPG6    ', iSUPG6   , 'MDFIXT    ', MDFIXT   ,  &
         &   'MDBEAM    ', MDBEAM   , &
         &   'MDOSQZ    ', MDOSQZ   , 'MDOSQZN   ', MDOSQZN  ,  &
         &   'MDLETA    ', MDLETA   , 'MDLNEO    ', MDLNEO   ,  &
         &   'MDBSETA   ', MDBSETA  , 'MDANOM    ', MDANOM   ,  &
         &   'MDITSN    ', MDITSN   , 'MDITST    ', MDITST   ,  &
         &   'MDINTN    ', MDINTN   , 'MDINTT    ', MDINTT   ,  &
         &   'MDINTC    ', MDINTC   , &
         &   'MDLETB    ', MDLETB   , 'IDIAG     ', IDIAG    ,  &
         &   'IGBDF     ', IGBDF    , 'ISMTHD    ', ISMTHD   ,  &
         &   'NTCOIL    ', NTCOIL   , 'MDLC      ', MDLC     ,  &
         &   'm_pol     ', m_pol    , 'n_tor     ', n_tor    ,  &
         &   'MDLNBD    ', MDLNBD   ,  &
!         &   'NCph      ', NCph     , 'NCth      ', NCth     ,  &
         &   'HPN(1,1)  ', HPN(1,1) , 'HPN(1,2)  ', HPN(1,2) , &
         &   'HPN(2,1)  ', HPN(2,1) , 'HPN(2,2)  ', HPN(2,2) , &
         &   'HPN(3,1)  ', HPN(3,1) , 'HPN(3,2)  ', HPN(3,2) , &
         &   'HPN(4,1)  ', HPN(4,1) , 'HPN(4,2)  ', HPN(4,2)
    RETURN
  END SUBROUTINE TXVIEW
end module tx_parameter_control
