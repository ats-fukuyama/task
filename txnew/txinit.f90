!!***************************************************************
!
!   Set constants and initial parameters
!
!***************************************************************

subroutine TXINIT
  use tx_commons
  use tx_graphic, only : MODEG, MODEGL, NGYRM, NGR, NGRSTP, NGTSTP, NGVSTP, gDIV

  implicit none
  integer(4) :: i, j

  !   Number of equations
  NQMAX = NQM

  !   ***** Configuration parameters *****

  !   Plasma minor radius (m), geometrically defined by (Rmax-Rmin)/2
  ra = 0.8d0

  !   Plasma minor radius (m), defined by sqrt(V/(2 Pi Pi R))
  ravl = 0.8d0

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
  RR = 3.2d0

  !   Toroidal magnetic field (T) at R=RR
  BB = 2.68d0

  !   Poloidal current function at the virtual wall (Tm)
  rbvt = rr * bb

  !   Plasma current start (MA)
  rIPs= 1.d0

  !   Plasma current end (MA)
  rIPe= 1.d0

  !   ***** Plasma components *****

  !   Atomic number
  amas(1) = aep     ! electron
  amas(2) = 2.d0    ! bulk ion (deuterium)
  amas(3) = 12.d0   ! impurity (carbon)
  amb     = amas(2) ! beam ion

  !   Charge number
  achg(1) = -1.d0
  achg(2) =  1.d0
  achg(3) =  6.d0
  achgb   =  1.d0 

  !   Effective charge
  Zeffin = 1.5d0

  !   ***** Initial plasma parameters *****

  !   Initial electron density at rho = 0 (10^20 m^-3)
  PN0 = 0.2d0

  !   Initial electron density at rho = 1 (10^20 m^-3)
  PNa = 0.05d0

  !   Electron density in divertor region (Minimum density in SOL)
  PNsDIV(1) = 0.01d0

  !   Ion density in divertor region (Minimum density in SOL)
  PNsDIV(2) = (achg(3) - Zeffin) / (achg(2) * (achg(3) - achg(2))) * PNsDIV(1)

  !   Impurity density in divertor region (Minimum density in SOL)
  PNsDIV(3) = (Zeffin - achg(2)) / (achg(3) * (achg(3) - achg(2))) * PNsDIV(1)

  !   Initial electron temperature at rho = 0 (keV)
  PTe0 = 2.d0

  !   Initial electron temperature at rho = 1 (keV)
  PTea = 0.2d0

  !   Electron temperature in divertor region (Minimum Te in SOL)
  PTsDIV(1) = 0.05d0

  !   Initial ion temperature  at rho = 0 (keV)
  PTi0 = 2.d0

  !   Initial ion temperature  at rho = 1 (keV)
  PTia = 0.2d0

  !   Ion temperature in divertor region (Minimum Ti in SOL)
  PTsDIV(2) = 0.05d0

  !   Initial impurity temperature  at rho = 0 (keV)
  PTz0 = 2.d0

  !   Initial impurity temperature  at rho = 1 (keV)
  PTza = 0.2d0

  !   Impurity temperature in divertor region (Minimum Ti in SOL)
  PTsDIV(3) = 0.05d0

  !   Initial current profile parameter
  PROFJ = 2.d0

  !   Initial density profile parameters
  PROFN1 = 2.d0
  PROFN2 = 1.d0

  !   Initial temperature profile parameters
  PROFT1 = 2.d0
  PROFT2 = 2.d0

  !   Initial ion toroidal rotation velocity at rho = 0 (m/s)
  Uiph0 = 0.d0

  !   Equilibrium parameters
  !     ieqread = 0 : Large aspect ratio limit
  !             = 1 : Large aspect ratio approximation
  !             = 2 : Read equilibrium parameters
  ieqread = 1

  !    irktrc  = 0           : bilinear interpolation
  !            = 1           : Runge-Kutta interpolation (rkck)
  !            = 2 (default) : Runge-Kutta interpolation (eqrk4)
  !            = negative    : Does not allow to enlarge the trace step size if tracing fails.
  irktrc = 2

  !   ***** Particle diffusivity and viscosity parameters *****

  !   Turbulent pinch velocity parameter
  VWpch0 = 0.d0

  !   Electron-driven diffusion parameter
  Dfs0(1) = 0.1d0

  !   Ion and impurity driven diffusion parameter (usually inactive)
  Dfs0(2:NSM) = 0.d0

  !   Electron, ion and impurity viscosity parameters
  rMus0(:) = 0.5d0

  !   Drift frequency parameter (omega/omega*e)
  WPM0 = 0.d0

  !  ***** Thermal diffusivity parameters *****

  !   Electron, ion and impurity thermal diffusivity parameter (Chi/D)
  !     0 for fixed temperature profile
  Chis0(1:NSM) = 0.5d0

  !   ***** Turbulent transport control parameters *****

  !   Fixed transport coefficient parameter
  !     1 : particle transport
  !     2 : momentum transport
  !     3 : thermal transport
  FSDFIX(1:3) = 1.d0

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
  FSANOM(1:3) = 0.d0

  !   Model of ExB shear stabilization
  !     0 : No ExB shear stabilization
  !     1 : Lorentz-type ExB shear stabilization (very weak)
  !     2 : Exponent-type ExB shear stabilization [M. Honda (2007) dissertation, Chapter 4]
  !     3 : Lorentz-type ExB shear stabilization [M. Yagi (2012) CPP p.372]
  FSCBSH = 0.d0

  !   Position of ETB shoulder (valid only when MDLETB /= 0)
  !     1 : particle transport
  !     2 : momentum transport
  !     3 : thermal transport
  RhoETB(1:3) = 0.9d0

  !   Momentum pinch parameter
  !     0 : No momentum pinch
  !     1 : Momentum pinch effective
  FSMPCH(1:NSM) = 1.d0

  !   Implementation of turbulent particle transport model
  !     0 : Electron density gradient term remains as it is.
  !     1 : Electron density gradient term is replaced using the radial force balance equation
  !         for the sake of numerical stability.
  !     2 : Electron density gradient term is replaced using the dimagnetic flow
  !         for the sake of numerical stability.
  FSTPTM(1:NSM) = 2

  !   ==== CDBM transport coefficient parameters ============================
  !      The finer time step size, typically less than or equal to DT=5.D-4,
  !         is anticipated when using CDBM.
  !      See also iprestab

  !   Effect of magnetic curvatures-alpha
  FSCBAL = 1.d0

  !   Effect of magnetic curvature
  !     (It could destabilize a numerical robustness when FS2 > FS1. )
  FSCBKP = 0.d0

  !   Effect of elongation (CDBM05 model)
  !      [M.Honda and A.Fukuyama NF 46 (2006) 580]
  FSCBEL = 1.d0

  !   Factor of E x B rotation shear
  rG1 = 10.d0

  !   ==== Neoclassical transport coefficient parameters ====================

  !   Projection of toroidal momentum flux onto parallel flow equations
  !   * These terms preserve consistency that has to be satisfied by LQ*2 with
  !     LQ*3 and LQ*4 at the magnetic axis, where terms of grad P and grad Phi
  !     become nil.
  !   * They are necessary not only for numerical reasons but also physical reasons.
  !     Hence, FSPARV should always be unity.
  FSPARV(1:NSM) = 1.d0

  !   Neoclassical thermal diffusivity parameter
  ChiNC = 1.d0

  !   Neoclassical transport parameters
  !     FSNC(1) : neoclassical viscosity
  !     FSNC(2) : neoclassical friction force
  FSNC(1:2) = 1.d0

  !   Beam neoclassical transport parameter
  !     FSNCB(1) : neoclassical viscosity
  !     FSNCB(2) : neoclassical friction force
  FSNCB(1:2) = 1.d0

  !   Helical neoclassical viscosity parameter
  FSHL = 0.d0

  !   =======================================================================

  !   Bohm transport coefficient parameter in SOL
  FSBOHM = 1.d0

  !   Pseudo-classical transport coefficient parameter in SOL
  !     1 : particle transport
  !     2 : momentum transport
  !     3 : thermal transport
  FSPCL(1:3) = 1.d0

  !   Controller for thermodiffusive pinch term of turbulent particle flux
  FSVAHL =-0.5d0

  !   Particle diffusion coefficient profile parameter (D(r=a)/D(r=0))
  PROFD =  5.d0

  !   Exponent of particle diffusion coefficient profile
  PROFD1 = 2.d0

  !   Gaussian modification of particle diffusion coefficient profile
  PROFD2 = 0.d0

  !   Pedestal of particle diffusion coefficient profile
  PROFDB = 0.d0

  !   Mom. diffusion coefficient profile parameter (D(r=a)/D(r=0))
  PROFM = 10.d0

  !   Exponent of mom. diffusion coefficient profile
  PROFM1 = 2.d0

  !   Pedestal of mom. diffusion coefficient profile
  PROFMB = 0.d0

  !   Heat diffusion coefficient profile parameter (chi(r=a)/chi(r=0))
  PROFC = 10.d0

  !   Exponent of heat diffusion coefficient profile
  PROFC1 = 2.d0

  !   Pedestal of heat diffusion coefficient profile
  PROFCB = 0.d0

  !   ***** Other transport parameters *****

  !   Advection parameter
  FSADV = 1.d0

  !   Beam advection parameter (finer time step required when on)
  FSADVB = 0.d0

  !   Velocity of flux surfaces
  FSUG = 1.d0

  !   Charge exchange parameter
  FSCX = 1.d0

  !   Orbit loss parameter
  !     FSLC = 0 : No orbit loss included.
  !            1 : Loss term is expressed as damping term.
  !          +10 : Loss term is expressed as source and sink term.
  FSLC = 0.d0

  !   Orbit loss model
  !     MDLC = 1 : Shaing model
  !            2 : Itoh model
  MDLC = 1

  !   Ripple loss parameter
  !     FSRP = 0.d0 : No ripple effect
  !            1.d0 : Full ripple effect
  !            2.d0 : Ripple effect but no ripple trapped ions
  FSRP = 0.d0

  !   Alpha heating paramter
  FSNF = 0.d0

  !   Particle loss to divertor parameter (default = 0.3)
  FSLP = 0.3d0

  !   Heat loss to divertor parameter (default = 1.0)
  FSLTs(1:NSM) = 1.d0

  !   Beam ion loss to divertor parameter (default = 0.3)
  FSLPB = 1.d0

  !   Ionization parameter
  FSION = 1.d0

  !   Neutral collision parameter
  FSNCOL = 1.d0

  !   Slow neutral diffusion factor
  FSD01 = 1.d0

  !   Thermal neutral diffusion factor
  FSD02 = 1.d0

  !   Halo neutral diffusion factor
  FSD03 = 1.d0

  !   Impurity neutral diffusion factor
  FSD0z = 1.d0

  !   ***** initial parameters *****

  !   Initial Density scale length in SOL normalized by minor radius (-), valid if MDITSN /= 0
  rLn = 0.0857d0

  !   Initail Temperature scale length in SOL normalized by minor radius (-), valid if MDITST /= 0
  rLT = 0.0857d0

  !   ***** Heating parameters *****

  !   Maximum NBI beam energy (keV)
  Ebmax = 80.d0

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
  PNBHP = 0.d0

  !   First tangential NBI input power (MW)
  PNBHT1 = 0.d0

  !   Second tangential NBI input power (MW)
  PNBHT2 = 0.d0

  !   NBI input power from input file (MW)
  PNBHex = 0.d0

  !   NBI current drive parameter ; Direction of injection
  !     1.0 : co, -1.0 : counter
  PNBCD = 1.d0

  !   Fraction of the collisional slowing down part of the perpendicular NBI
  !     -1.0 (counter) <= PNBMPD <= 1.0 (co)
  !     0.0 : Exact perpendicular injection without momentum input
  PNBMPD = 0.d0

  !   Pitch angle between the tanjential NBI chord and the magnetic field line
  PNBPTC = 0.8d0

  !   Different NBI deposition profiles for electrons and ions due to banana orbit effect
  !     MDLNBD = 0 : No charge separation
  !              1 : Orbit effect for banana particles only
  !              2 : Orbit effect for all beam ions including passing particles
  MDLNBD = 0

  !   Refractive index of RF waves
  rNRFe = 0.d0
  rNRFi = 0.d0

  !   Heating width of RF heating (-)
  RRFew = 0.5d0
  RRFiw = 0.5d0

  !   Heating center radius of RF heating (-)
  RRFe0 = 0.d0
  RRFi0 = 0.d0

  !   RF input power (MW)
  PRFHe = 0.d0
  PRFHi = 0.d0

  !   Additional torque input (N m)
  Tqt0  = 0.d0  ! Toroidal
  Tqp0  = 0.d0  ! Sheared torque

  !   ***** Neutral parameters *****

  !   Initial Neutral density (10^20 m^-3)
  PN0s  = 1.D-8 ! Hydrogen species
  PN0zs = 1.D-8 ! Impurity species

  !   Neutral thermal velocity (m/s)
  !     V0 = SQRT(2.d0*X[eV]*AEE/(amas(2)*AMP)),
  !     where X=3eV is assumed (see p. 529 of [4])
  V0  = 1.6954D4 ! Hydrogen species
  V0z = 6.9206D3 ! Impurity species

  !   Recycling rate in SOL
  rGamm0  = 0.8d0 ! Hydrogen species
  rGamm0z = 0.8d0 ! Impurity species

  !   Gas-puff particle flux (10^20 m^-2 1/s)
  !      If you input Gamma_0 [particles/sec], you translate it into
  !      rGASPF [1/(m^2 s)]:
  !        rGASPF = Gamma_0 / <|nabla V|>
  !               = Gamma_0 / surt(NRMAX)
  !               ~ Gamma_0 / (2.d0*PI*RR*2.d0*PI*RB)
  rGASPF  = 0.2d0 ! Hydrogen species
  rGASPFz = 0.0d0 ! Impurity species

  !   ***** Ripple loss parameters *****

  ! Number of toroidal field coils in JT-60U
  NTCOIL = 18

  ! Minimum ripple amplitude in the system (R = Rmag0 = 2.4m on JT-60U)
  DltRPn = 0.0002d0 * 1.D-2

  ! Effective ellipticity for ripple amplitude estimation
  kappa = 1.2d0

  !   ***** Parameters for toroidal neoclassical viscosity *****

  ! Maximum poloidal mode number of error field (should be power of two)
  m_pol  = 32

  ! Toroidal mode number of error field
  n_tor  = NTCOIL

  ! Magnetic braiding parameters ***AF (2008-06-08)
  DMAG0   = 0.d0  ! Magnetic field line diffusivity : (Delta r)**2/Delta z [m]
  RMAGMN  = RA    ! minimum radius of the braiding region [m]
  RMAGMX  = RA    ! maximum radius of the braiding region [m]

  !   Helical ripple amplitude at r=a, linear in (r/a)
  !  EpsH = 0.1d0     ! amplitude 
  EpsHM(1:NHFMmx,0:3) = 0.d0
  !  EpsHM(1,0:3) = 0.d0, 0.d0, 0.1d0, 0.d0     ! amplitude 
  EpsHM(1,2) = 0.1d0
  !  NCph = 5        ! toroidal pitch number
  !  NCth = 2        ! poloidal pitch number
  HPN(1:NHFMmx, 1:2) = 0
  HPN(1,1) = 2    ! poloidal pitch number
  HPN(1,2) = 5    ! toroidal pitch number
  Q0 = 3.d0       ! q(0) by external coils
  QA = 2.d0       ! q(a) by external coils

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
  ICMAX = 100

  !   Time-advancing method
  !     ADV = 0     : Explicit scheme       (Not usable)
  !           0.5   : Crank-Nicolson scheme (Not usable)
  !           2/3   : Galerkin scheme       (Not usable)
  !           0.729 : Minimum value at which a simulation can be run.
  !           0.878 : Liniger scheme
  !           1     : Implicit scheme       (Recommended)
  ADV = 1.d0

  !   Mode of Backward Differential (Gear's) Formula (Second-order accuracy)
  !   (http://www.scholarpedia.org/article/Backward_differentiation_formulas)
  !   (Bibun Houteishiki no Suuchi Kaihou I, Taketomo Mitsui, Iwanami, 1993, p29,30)
  !   (P.M.Gresho and R.L.Sani, "Incompressible Flow and the Finite Element Method,
  !    Vol.1 Advection-Diffusion", John Wiley and Sons, 2000, p.263)
  !     IGBDF = 0   : not used
  !             1   : Use BDF2
  IGBDF = 1

  !   Mode of the scheme to evaluate Er and dPs/dpsi for suppressing wiggles near the axis
  !     ISMTHD = 0  : use dPhi/dV to calculate Er and differentiate p w.r.t. psi to obtain dPs/dpsi
  !              1  : use dPhi/drho with dPhi/drho = 0 at rho = 0 to calculate Er
  !                   and differentiate p w.r.t. rho to obtain dPs/drho with dPs/drho at rho = 0
  !              11 : dPhi/drho(1) and dPs/drho(1) are interpolated by "replace_interpolate_value".
  !                   Smooth profiles regarding diamag. flow would be obtained.
  ISMTHD = 1

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
  !     iSUPG3 for LQ*3CC; LQb3CC only when iSUPG3 = -1
  !     iSUPG6 for LQ*6CC
  !     iSUPG8 for LQ*8CC
  iSUPG3 = 0 ! Able to be activated only when FSPARV(:) = 0.0
  iSUPG6 = 0
  iSUPG8 = 1

  !   Enhancement factor for numerical stability in advection problems
  !     p = SUPGstab / sqrt(15.d0)
  !     where SUPGstab = 1.d0 : [Raymond and Garder, Monthly Weather Review (1976)] suitable for the case with beam
  !                    = 2.d0 : [Ganjoo and Tezduyar, NASA report CR-180124 (1986)]
  !                    = 4.d0 : Seemingly best, but no ground in the light of numerical analysis
  SUPGstab = 1.d0

  !   Model for stabilizing oscillation due to the nonlinearity of the pressure gradient
  !   when typically using CDBM model.
  !     0 : Doing nothing (no stabilization; often not work)
  !     1 : Mixing pressure in that instant and that at IC = 1 (good)
  !     2 : Mixing pressure in that instant and that converged at previous time
  !         (slightly unstable when ExB shearing is on with larger time step size.)
  !     3 : Mixing pressure at IC = 1 and that converged at previous time (excessive assumption)
  !
  !   N.B. Dt=1.e-3 is ok, but dt=1.e-4 realizes more stable calculation.
  !        Choosing iprestab=3 is that diffusivities are estimated by using the pressure gradient
  !        at almost previous time, which is fixed during iteration.
  iprestab = 1

  !   ***** Mesh number parameters *****

  !   Magnitude of mesh peakness
  CMESH0 =  2.d0
  CMESH  = 10.d0

  !   Width of mesh peakness (-)
  !  WMESH0 = 0.2d0
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

  !   Debug output control parameter for Matrix Inversion
  !   0 : OFF, 1 : ON
  !     midbg(1) : check the symmetry of friction coefficients
  !                simulation will stop right after output if midbg(1)=0.
  !     midbg(2) : check the ambipolarity of the flux
  midbg(1) = 0
  midbg(2) = 0

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
  !   1    : more    ; LQb1, LQb2, LQb3, LQb4, LQb7, LQb8
  !   2    : no flow ; LQb1
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
  !   1    : Matrix Inversion (booth9)
  !   2    : Matrix Inversion (nccoe)
  !   3    : NCLASS
  !   +10  : Both (Users cannot specify)
  MDLNEO  = 2

  !   Options of neoclassical transport model
  !   Valid only when MDLNEO = 1 or 2
  imodel_neo(1) = 0 ! Fast ion viscosity
  imodel_neo(2) = 1 ! Required for NBCD
  imodel_neo(3) = 2 ! Fast ion contribution to friction forces
  imodel_neo(4) = 0 ! Unused
  imodel_neo(5) = 0 ! PS contribution nil when 1
  imodel_neo(6) = 3 ! Higher-order flow contribution, valid only for nccoe

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
  !   +-10 : smoothing profiles to ensure that the first-order derivative of them
  !          is continuous and smooth
  MDINTN = 10

  !   Mode of initial temperature profiles
  !   -2   : read from file and smooth
  !   -1   : read from file
  !   0    : original
  !   1    : pedestal model
  !   2    : empirical steady state temperature profile
  !   +-10 : smoothing profiles to ensure that the first-order derivative of them
  !          is continuous and smooth
  MDINTT = 12

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

  !   Poynting flux calculation
  !  << Details are shown in [M. Honda, Nucl. Fusion 58 (2018) 026006] >>
  !    iPoyntpol and iPoynttor are flags to choose Poynting flux equations:
  !    iPoyntpol = 1 and iPoynttor = 0  ==>  eq. (35) in the reference
  !    iPoyntpol = 0 and iPoynttor = 1  ==>  eq. (36) in the reference
  !    iPoyntpol = 1 and iPoynttor = 1  ==>  eq. (25) in the reference
  iPoyntpol = 1
  iPoynttor = 1

!!$  !   Multiplication factor for graphic in the radial direction (module tx_graphic)
!!$  !   default : 1.0
!!$  !
!!$  gDIV(:)   = 1.0
!!$  gDIV(  1) = gkilo
!!$  gDIV(  2) = 1.e-3
!!$  gDIV(  5) = gkilo
!!$  gDIV(  8) = gkilo
!!$  gDIV(  9) = gkilo
!!$  gDIV( 11) = gkilo
!!$  gDIV( 12) = gkilo
!!$  gDIV( 13) = gkilo
!!$  gDIV( 16) = gkilo
!!$  gDIV( 17) = gkilo
!!$  gDIV( 19) = gkilo
!!$  gDIV( 20) = gkilo
!!$  gDIV( 21) = gkilo
!!$  gDIV( 24) = gkilo
!!$  gDIV( 25) = gkilo
!!$  gDIV( 27) = gkilo
!!$  gDIV( 28) = gkilo
!!$  gDIV( 29) = gkilo
!!$  gDIV( 30) = 1.e-2
!!$  gDIV( 32) = gkilo
!!$  gDIV( 33) = gkilo
!!$  gDIV( 34) = gkilo
!!$  gDIV( 35) = gkilo
!!$  gDIV( 36) = 1.e-4
!!$  gDIV( 37) = 1.e-6
!!$  gDIV( 38) = 1.e-6
!!$  gDIV( 39) = 1.e-4
!!$  gDIV( 40) = 1.e-4
!!$
!!$  do i = 1, NSM
!!$     j = (i - 1) * 13
!!$     gDIV( 52+j) = 1.e20
!!$     gDIV( 54+j) = gmega
!!$     gDIV( 56+j) = gkilo
!!$     gDIV( 57+j) = gkilo
!!$     gDIV( 58+j) = gkilo
!!$     gDIV( 59+j) = gkilo
!!$     gDIV( 60+j) = gkilo
!!$     gDIV( 62+j) = gkilo
!!$     gDIV( 63+j) = gkilo
!!$     gDIV( 64+j) = gkilo
!!$  end do
!!$  gDIV(80) = gkilo
!!$
!!$  gDIV( 91) = gkilo
!!$  gDIV( 92) = 1.e-3
!!$  gDIV(99:109) = gmega
!!$  gDIV(110:112) = gkilo
!!$
!!$  gDIV(114) = 1.e18
!!$  gDIV(116:121) = gkilo
!!$  gDIV(122) = 1.e15
!!$  gDIV(125) = 1.e15
!!$  gDIV(126:128) = 1.e13
!!$
!!$  gDIV(129:137) = gkilo
!!$  gDIV(138) = 1.e14
!!$  gDIV(157:159) = gkilo
!!$  gDIV(179) = gkilo
!!$  gDIV(192) = gmega
!!$
!!$  gDIV(208:219) = 1.e-6
!!$  gDIV(220:231) = gkilo
!!$  gDIV(256) = 1.e-10
!!$
!!$  gDIV(271:272) = 1.e20
!!$  gDIV(274:283) = gmega
!!$  gDIV(286:288) = gmega
!!$  gDIV(290) = gkilo
!!$  gDIV(292) = gmega
!!$  gDIV(295) = gkilo
!!$
!!$  gDIV(338:341) = 1.e20
!!$  gDIV(342:347) = gkilo

  !   *** Density perturbation technique ***

  !   Normalized radius where density increase
  DelRho = 0.2d0

  !   Amount of increase in density (10^20 m^-3)
  DelN = 1.D-1

  !   ***

  !   Index for graphic save interval (module tx_graphic)

  NGR=-1

end subroutine TXINIT

!***************************************************************
!
!        Calculate mesh, etc.
!
!***************************************************************

subroutine TXCALM

  use tx_commons
  use tx_interface, only : LORENTZ, LORENTZ_PART, BISECTION

  implicit none
  integer(4) :: NR, NRL, nr_rhoc_near
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
  MAXAMP = LORENTZ(rhob,C1L,C2L,W1L,W2L,0.d0,rhoc) / rhob ! normalized parameter C_0
  nr_rhoc_near = 0
!  R(0) = 0.d0
  do NR = 1, NRMAX - 1
     rhol = NR * rhob / NRMAX
     call BISECTION(LORENTZ,C1L,C2L,W1L,W2L,0.d0,rhoc,MAXAMP,rhol,rhob,rho(NR))
     if(abs(rho(NR)-rhoc) <= abs(rho(NR)-rho(nr_rhoc_near))) nr_rhoc_near = NR
  end do
  rho(NRMAX) = rhob

  !  Construct new CL value that separatrix is just on mesh.
  !  New CL is chosen in order not to be settle so far from given CL.
  !  The mesh finally obtained is well-defined.
  rhocl = nr_rhoc_near * rhob / NRMAX
  CLNEW = ( (rhoc - rhocl) * rhob - rhocl * C1L * LORENTZ_PART(rhob,W1L,W2L,0.d0,rhoc,0) &
       &   + rhob * C1L * LORENTZ_PART(rhoc,W1L,W2L,0.d0,rhoc,0)) &
       &  / (  rhocl * LORENTZ_PART(rhob,W1L,W2L,0.d0,rhoc,1) &
       &     - rhob  * LORENTZ_PART(rhoc,W1L,W2L,0.d0,rhoc,1))
  MAXAMP = LORENTZ(rhob,C1L,CLNEW,W1L,W2L,0.d0,rhoc) / rhob
  rho(0) = 0.d0
  do NR = 1, NRMAX - 1
     rhol = NR * rhob / NRMAX
     call BISECTION(LORENTZ,C1L,CLNEW,W1L,W2L,0.d0,rhoc,MAXAMP,rhol,rhob,rho(NR))
  end do
  rho(NRMAX) = rhob

  !  Maximum NR till rhoc

  do NR = 0, NRMAX-1
     if(rho(NR) <= rhoc .and. rho(NR+1) >= rhoc) then
        NRL = NR
        exit
     end if
  end do

  !  Adjust rhoc on mesh (for the time being rhoc is assumed to be equivalent to 1.d0)

  if(abs(rho(NRL)-rhoc) < abs(rho(NRL+1)-rhoc)) then
     NRA = NRL
     rhol = 0.5d0 * abs(rho(NRA) - rhoc)
     rho(NRA  ) = rhoc
     rho(NRA-1) = rho(NRA-1) + rhol
  else
     NRA = NRL + 1
     rhol = 0.5d0 * abs(rho(NRA) - rhoc)
     rho(NRA  ) = rhoc
     rho(NRA+1) = rho(NRA+1) - rhol
  end if

  !  rhoaccum doesn't coincide with the plasma surface 1.d0

  if(rhoaccum >= 0.d0) then

     !  Maximum NR till rho=1.d0

     do NR = 0, NRMAX-1
        if(rho(NR) <= 1.d0 .and. rho(NR+1) >= 1.d0) then
           NRL = NR
           exit
        end if
     end do

     !  Adjust rho=1.d0 on mesh

     if(abs(rho(NRL)-1.d0) < abs(rho(NRL+1)-1.d0)) then
        NRA = NRL
        rhol = 0.5d0 * abs(rho(NRA) - 1.d0)
        rho(NRA  ) = 1.d0
        rho(NRA-1) = rho(NRA-1) + rhol
     else
        NRA = NRL + 1
        rhol = 0.5d0 * abs(rho(NRA) - 1.d0)
        rho(NRA  ) = 1.d0
        rho(NRA+1) = rho(NRA+1) - rhol
     end if
     
  end if

  !  Mesh number at the center of the core plasma

  do NR = 0, NRMAX-1
     if(rho(NR) <= 0.5d0.and.rho(NR+1) >= 0.5d0) then
        if(abs(rho(NR)-0.5d0) < abs(rho(NR+1)-0.5d0)) then
           NRC = NR
        else
           NRC = NR + 1
        end if
        exit
     end if
  end do

  !  Equilibrium
  call txequ

end subroutine TXCALM

!***************************************************************
!
!   Initialize profiles
!
!***************************************************************

subroutine TXPROF

  use tx_commons
  use tx_graphic, only : NGR, NGT, NGVV
  use tx_variables
  use tx_interface, only : INTDERIV3, detect_datatype, dfdx, &
       &                   initprof_input, moving_average, CORR, coulog, inexpolate
  use tx_core_module, only : intg_area, intg_area_p, intg_vol_p 
  use tx_glob, only : TXGLOB
  use sauter_mod
  use tx_ntv, only : perturb_mag, Wnm_spline
  use eqread_mod, only : AJphVRL
  use mod_eqneo, only : wrap_eqneo
  use mod_cross_section, only : spline_table_carbon_rate_coef_adas, spline_table_beam_rate_coef
  use libitp, only : FCTR4pt
  use libspl1d, only : spl1d,spl1df

  implicit none
  integer(4) :: NR, IER, ifile, NHFM, NR_smt, NR_smt_start = 10
  integer(4) :: MDINTN1, MDINTN2, MDINTT1, MDINTT2
  real(8) :: rhol, PROFN, PROFT, PTePROF, PTiPROF, PTzPROF!, RL, QL, dRIP
  real(8) :: AJFCT, SUM_INT, DR1, DR2
  real(8) :: EpsL, FTL, PBA, dPN, sigma, fexp, CfN1, CfN2, PN0L, PNaL, PNeDIVL 
  real(8) :: pea, pia, pza, pediv, pidiv, pzdiv, dpea, dpia, dpza, faci, facz, &
       &     Cfpe1, Cfpe2, Cfpi1, Cfpi2, Cfpz1, Cfpz2, &
       &     PTe0L, PTi0L, PTz0L, PTeaL, PTiaL, PTzaL, PTeDIVL, PTiDIVL, PTzDIVL
  real(8) :: BCLQm3, etanc, etaspz, dum=0.d0, tmp
  real(8), dimension(:), allocatable :: AJPHL, ProfR, Profsdt, dProfsdt
  real(8), dimension(:,:), allocatable :: ProfRS

  MDINTN2 = abs(MDINTN / 10)      ! 2nd digit of MDINTN
  MDINTN1 = MDINTN - MDINTN2 * 10 ! 1st digit of MDINTN
  MDINTT2 = abs(MDINTT / 10)      ! 2nd digit of MDINTT
  MDINTT1 = MDINTT - MDINTT2 * 10 ! 1st digit of MDINTT

  !  Read spline table for rate coefficients
  call spline_table_carbon_rate_coef_adas
  call spline_table_beam_rate_coef

  !  Read spline table for neoclassical toroidal viscosity
  if(FSRP /= 0.d0) call Wnm_spline

  NEMAX = NRMAX

  !  Define basic quantities like mass of particles, mesh, etc.

  call TXCALM

  !  Contribution of perturbed magnetic field
!  if(DltRPn /= 0.d0) call perturb_mag

  !  Variables

  ! ********************************************************************
  !      B.C. and profiles of density, temperature, pressure
  ! ********************************************************************

  if(MDINTN1 < 0 .or. MDINTT1 < 0 .or. abs(MDINTC) /= 0) call initprof_input
  ! === Density at the boundaries ===
  if(MDINTN1 < 0) then
     ! --> Read from the file
     call initprof_input(  0,1,PN0L)
     call initprof_input(NRA,1,PNaL)
     PNeDIVL = 0.25d0 * PNaL
  else
     ! --> Read from the namelist
     PN0L = PN0
     PNaL = PNa
     PNeDIVL = PNsDIV(1)
  end if
  ! === Temperature at the boundaries ===
  if(MDINTT1 < 0) then
     ! --> Read from the file
     call initprof_input(  0,2,PTe0L)
     call initprof_input(  0,3,PTi0L)
     call initprof_input(NRA,2,PTeaL)
     call initprof_input(NRA,3,PTiaL)
     PTz0L   = PTi0L
     PTzaL   = PTiaL
     PTsDIV(1) = 0.25d0 * PTeaL
     PTsDIV(2) = 0.25d0 * PTiaL
     PTsDIV(3) = PTsDIV(2)
  else
     ! --> Read from the namelist
     PTe0L   = PTe0
     PTi0L   = PTi0
     PTz0L   = PTz0
     PTeaL   = PTea
     PTiaL   = PTia
     PTzaL   = PTza
  end if
  PTeDIVL = PTsDIV(1)
  PTiDIVL = PTsDIV(2)
  PTzDIVL = PTsDIV(3)

  PBA   = rhob - 1.d0
  dPN   = - 3.d0 * (PN0L - PNaL) / ravl
  CfN1  = - (3.d0 * PBA * dPN + 4.d0 * (PNaL - PNeDIVL)) / PBA**3
  CfN2  =   (2.d0 * PBA * dPN + 3.d0 * (PNaL - PNeDIVL)) / PBA**4
  faci  = (achg(3)-Zeffin)/(achg(2)*(achg(3)-achg(2)))
  facz  = (Zeffin-achg(2))/(achg(3)*(achg(3)-achg(2)))
  if(MDFIXT == 0) then
     pea   = PNaL    * PTeaL
     dpea  = dPN     * PTeaL
     pediv = PNeDIVL * PTeDIVL
     pia   = faci*PNaL    * PTiaL
     dpia  = faci*dPN     * PTiaL
     pidiv = faci*PNeDIVL * PTiDIVL
     pza   = facz*PNaL    * PTzaL
     dpza  = facz*dPN     * PTzaL
     pzdiv = facz*PNeDIVL * PTzDIVL
  else
     pea   = PTeaL   ; pia   = PTiaL   ; pza   = PTzaL
     dpea  = 0.d0    ; dpia  = 0.d0    ; dpza  = 0.d0
     pediv = PTeDIVL ; pidiv = PTiDIVL ; pzdiv = PTzDIVL
  end if
  Cfpe1 = - (3.d0 * PBA * dpea + 4.d0 * (pea - pediv)) / PBA**3
  Cfpe2 =   (2.d0 * PBA * dpea + 3.d0 * (pea - pediv)) / PBA**4
  Cfpi1 = - (3.d0 * PBA * dpia + 4.d0 * (pia - pidiv)) / PBA**3
  Cfpi2 =   (2.d0 * PBA * dpia + 3.d0 * (pia - pidiv)) / PBA**4
  Cfpz1 = - (3.d0 * PBA * dpza + 4.d0 * (pza - pzdiv)) / PBA**3
  Cfpz2 =   (2.d0 * PBA * dpza + 3.d0 * (pza - pzdiv)) / PBA**4

  ! === Density profile ===
  do NR = 0, NRMAX
     rhol = rho(NR)
     if( rhol <= 1.d0 ) then
        ! +++ Core +++
        if( MDINTN1 < 0 ) then
           ! --> Read from the file
           call initprof_input(NR,1,X(NR,LQe1)) ! Ne [10^{20}m^{-3}]
           X(NR,LQi1) = faci * X(NR,LQe1)       ! Ni [10^{20}m^{-3}]
           X(NR,LQz1) = facz * X(NR,LQe1)       ! Nz [10^{20}m^{-3}]
        else
           ! --> Read from the namelist
           PROFN = (1.d0 - rho(NR)**PROFN1)**PROFN2
           X(NR,LQe1) = (PN0L - PNaL) * PROFN + PNaL ! Ne [10^{20}m^{-3}]
           X(NR,LQi1) = faci * X(NR,LQe1)            ! Ni [10^{20}m^{-3}]
           X(NR,LQz1) = facz * X(NR,LQe1)            ! Nz [10^{20}m^{-3}]
        end if

     else
        ! +++ SOL +++
        if( MDITSN == 0 ) then
           X(NR,LQe1) = PNaL + dPN * (rhol - 1.d0) + CfN1 * (rhol - 1.d0)**3 &
                &                                  + CfN2 * (rhol - 1.d0)**4
           X(NR,LQi1) = faci * X(NR,LQe1)
           X(NR,LQz1) = facz * X(NR,LQe1)
        else
           X(NR,LQe1) = PNaL * exp(- (rhol - 1.d0) / rLn)
           X(NR,LQi1) = faci * X(NR,LQe1)
           X(NR,LQz1) = facz * X(NR,LQe1)
        end if
     end if

     ! === Fixed densities to keep them constant during iterations ===
     PNsV_FIX(NR,1) = X(NR,LQe1)
     PNsV_FIX(NR,2) = X(NR,LQi1)
     PNsV_FIX(NR,3) = X(NR,LQz1)
  end do

  ! === Temperature and pressure profiles ===
  do NR = 0, NRMAX
     rhol = rho(NR)
     if( rhol <= 1.d0 ) then
        ! +++ Core +++
        if( MDINTT1 < 0 ) then
           ! --> Read from the file
           call initprof_input(NR,2,PTePROF) ! Te
           call initprof_input(NR,3,PTiPROF) ! Ti
           PTzPROF = PTiPROF
        else if( MDINTT1 == 0 ) then
           ! --> Read from the namelist
           PROFT = (1.d0 - rho(NR)**PROFT1)**PROFT2
           PTePROF = (PTe0L - PTeaL) * PROFT + PTeaL
           PTiPROF = (PTi0L - PTiaL) * PROFT + PTiaL
           PTzPROF = (PTz0L - PTzaL) * PROFT + PTzaL
        else if( MDINTT1 == 1 ) then
           ! --> Read from the namelist
           PTePROF = (PTe0L - PTeaL) * (0.8263d0 * (1.d0 - rho(NR)**2 )**1.5d0 &
                &                     + 0.167d0  * (1.d0 - rho(NR)**30)**1.25d0) + PTeaL
           PTiPROF = (PTi0L - PTiaL) * (0.8263d0 * (1.d0 - rho(NR)**2 )**1.5d0 &
                &                     + 0.167d0  * (1.d0 - rho(NR)**30)**1.25d0) + PTiaL
           PTzPROF = (PTz0L - PTzaL) * (0.8263d0 * (1.d0 - rho(NR)**2 )**1.5d0 &
                &                     + 0.167d0  * (1.d0 - rho(NR)**30)**1.25d0) + PTzaL
        else
           ! --> Read from the namelist
           PTePROF = PTe0L * exp(- log(PTe0L / PTeaL) * rho(NR)**2)
           PTiPROF = PTi0L * exp(- log(PTi0L / PTiaL) * rho(NR)**2)
           PTzPROF = PTz0L * exp(- log(PTz0L / PTzaL) * rho(NR)**2)
        end if

        if( MDFIXT == 0 ) then
           ! In case that LQ*5 represents pressure variables
           X(NR,LQe5) = PTePROF * X(NR,LQe1) ! Ne*Te [10^{20}m^{-3}*keV]
           X(NR,LQi5) = PTiPROF * X(NR,LQi1) ! Ni*Ti [10^{20}m^{-3}*keV]
           X(NR,LQz5) = PTzPROF * X(NR,LQz1) ! Nz*Tz [10^{20}m^{-3}*keV]
        else
           ! In case that LQ*5 represents temperature variables
           X(NR,LQe5) = PTePROF ! Te
           X(NR,LQi5) = PTiPROF ! Ti
           X(NR,LQz5) = PTzPROF ! Tz
        end if
        
     else
        ! +++ SOL +++
        if( MDITST == 0 ) then
           X(NR,LQe5) = pea + dpea * (rhol - 1.d0) + Cfpe1 * (rhol - 1.d0)**3 &
                &                                  + Cfpe2 * (rhol - 1.d0)**4
           X(NR,LQi5) = pia + dpia * (rhol - 1.d0) + Cfpi1 * (rhol - 1.d0)**3 &
                &                                  + Cfpi2 * (rhol - 1.d0)**4
           X(NR,LQz5) = pza + dpza * (rhol - 1.d0) + Cfpz1 * (rhol - 1.d0)**3 &
                &                                  + Cfpz2 * (rhol - 1.d0)**4
        else
           if( MDITST == 1 ) then
              PTePROF = PTeaL * exp(- (rhol - 1.d0) / rLT)
              PTiPROF = PTiaL * exp(- (rhol - 1.d0) / rLT)
              PTzPROF = PTzaL * exp(- (rhol - 1.d0) / rLT)
           else
              sigma = 0.3d0 * (rho(NRMAX) - 1.d0)
              fexp  = exp(- (rho(NR) - 1.d0)**2 / (2.d0 * sigma**2))
              PTePROF = PTe0L   * exp(- log(PTe0L / PTeaL) * rho(NR)**2) *         fexp &
                   &  + PTeDIVL                                          * (1.d0 - fexp)
              PTiPROF = PTi0L   * exp(- log(PTi0L / PTiaL) * rho(NR)**2) *         fexp &
                   &  + PTiDIVL                                          * (1.d0 - fexp)
              PTzPROF = PTz0L   * exp(- log(PTz0L / PTzaL) * rho(NR)**2) *         fexp &
                   &  + PTzDIVL                                          * (1.d0 - fexp)
           end if
           ! Pressure (MDFIXT=0) or temperature (MDFIXT=1)
           if( MDFIXT == 0 ) then
              X(NR,LQe5) = PTePROF * X(NR,LQe1)
              X(NR,LQi5) = PTiPROF * X(NR,LQi1)
              X(NR,LQz5) = PTzPROF * X(NR,LQz1)
           else
              X(NR,LQe5) = PTePROF
              X(NR,LQi5) = PTiPROF
              X(NR,LQz5) = PTzPROF
           end if
        end if
     end if

     ! === Fixed temperatures to keep them constant during iterations ===
     PTsV_FIX(NR,1) = PTePROF
     PTsV_FIX(NR,2) = PTiPROF
     PTsV_FIX(NR,3) = PTzPROF
  end do

  ! === Neutrals ===
  do NR = 0, NRMAX
     ! N0_1 (slow neutrals)
     X(NR,LQn1) = PN0s
     ! N0_2 (thermal neutrals)
     X(NR,LQn2) = 1.D-20 ! when ThntSW = 0.d0
     ! N0_3 (halo neutrals)
     X(NR,LQn3) = 0.d0
     ! N0z  (Impurity neutrals)
     X(NR,LQnz) = PN0zs
  end do

  ! ********************************************************************
  !      Smoothing profiles
  ! ********************************************************************

  ! === Temperature and pressure profiles ===
  if( MDINTT1 == -2 ) then
     allocate(ProfRS, mold=array_init_NRNS)
     if( MDFIXT == 0 ) then
        ProfRS(:,1) = X(:,LQe5) / X(:,LQe1)
        ProfRS(:,2) = X(:,LQi5) / X(:,LQi1)
        ProfRS(:,3) = X(:,LQz5) / X(:,LQz1)
     else
        ProfRS(:,1) = X(:,LQe5)
        ProfRS(:,2) = X(:,LQi5)
        ProfRS(:,3) = X(:,LQz5)
     end if
     NR_smt = NRA - NR_smt_start ! smoothing data only in the edge region
     if( MDFIXT == 0 ) then
        do NR = NR_smt, NRMAX
           X(NR,LQe5) = moving_average(NR,ProfRS(:,1),NRMAX) * X(NR,LQe1)
           X(NR,LQi5) = moving_average(NR,ProfRS(:,2),NRMAX) * X(NR,LQi1)
           X(NR,LQz5) = moving_average(NR,ProfRS(:,3),NRMAX) * X(NR,LQz1)
        end do
     else
        do NR = NR_smt, NRMAX
           X(NR,LQe5) = moving_average(NR,ProfRS(:,1),NRMAX)
           X(NR,LQi5) = moving_average(NR,ProfRS(:,2),NRMAX)
           X(NR,LQz5) = moving_average(NR,ProfRS(:,3),NRMAX)
        end do
     end if
     deallocate(ProfRS)
  end if

  if( MDINTT2 == 1 ) then
     ! Render the first derivative of the profile continuous and smooth using spline
     block
       integer(4) :: i, ierr
       real(8), dimension(:), allocatable :: rhoeq, deriv
       real(8), dimension(:,:), allocatable :: u

       allocate(rhoeq, deriv, ProfR, source=array_init_NR)
       do NR = 0, NRMAX
          rhoeq(NR) = NR / real(NRMAX,8) * rhob
       end do

       allocate(ProfRS, mold=array_init_NRNS)
       if( MDFIXT == 0 ) then
          ProfRS(:,1) = X(:,LQe5) / X(:,LQe1)
          ProfRS(:,2) = X(:,LQi5) / X(:,LQi1)
          ProfRS(:,3) = X(:,LQz5) / X(:,LQz1)
       else
          ProfRS(:,1) = X(:,LQe5)
          ProfRS(:,2) = X(:,LQi5)
          ProfRS(:,3) = X(:,LQz5)
       end if

       allocate(u(0:NRMAX,4))
       do i = 1, NSM
          call spl1d(rho,ProfRS(:,i),deriv,u,NRMAX+1,0,ierr)
          ! Map ProfRS on unequally-spaced rho onto equally-spaced rho
          do NR = 0, NRMAX
             call spl1df(rhoeq(NR),ProfR(NR),rho,u,NRMAX+1,ierr)
          end do
          call spl1d(rhoeq,ProfR,deriv,u,NRMAX+1,0,ierr)
          ! Map ProfR on equally-spaced rho back onto unequally-spaced rho
          ! This process can render the first-derivative of ProfR continuous and smooth.
          do NR = 0, NRMAX
             call spl1df(rho(NR),ProfRS(NR,i),rhoeq,u,NRMAX+1,ierr)
          end do
       end do

       if( MDFIXT == 0 ) then
          X(:,LQe5) = ProfRS(:,1) * X(:,LQe1)
          X(:,LQi5) = ProfRS(:,2) * X(:,LQi1)
          X(:,LQz5) = ProfRS(:,3) * X(:,LQz1)
       else
          X(:,LQe5) = ProfRS(:,1)
          X(:,LQi5) = ProfRS(:,2)
          X(:,LQz5) = ProfRS(:,3)
       end if

       deallocate(rhoeq,deriv,ProfR,ProfRS,u)
     end block
  end if
  
  ! === Density profile ===
  if(MDINTN1 == -2) then
     allocate(ProfR, mold=array_init_NR)
     if( MDFIXT == 0 ) then
        allocate(ProfRS, mold=array_init_NRNS)
        ProfRS(:,1) = X(:,LQe5) / X(:,LQe1)
        ProfRS(:,2) = X(:,LQi5) / X(:,LQi1)
        ProfRS(:,3) = X(:,LQz5) / X(:,LQz1)
     end if

     ProfR(:) = X(:,LQe1)
     NR_smt = NRA - NR_smt_start ! smoothing data only in the edge region
     do NR = NR_smt, NRMAX
        X(NR,LQe1) = moving_average(NR,ProfR,NRMAX)
        X(NR,LQi1) = faci * X(NR,LQe1) ! Ni
        X(NR,LQz1) = facz * X(NR,LQe1) ! Nz
     end do

     if( MDFIXT == 0 ) then
        X(:,LQe5) = ProfRS(:,1) * X(:,LQe1)
        X(:,LQi5) = ProfRS(:,2) * X(:,LQi1)
        X(:,LQz5) = ProfRS(:,3) * X(:,LQz1)
        deallocate(ProfRS)
     end if
     deallocate(ProfR)
  end if

  if( MDINTN2 == 1 ) then
     ! Render the first derivative of the profile continuous and smooth using spline
     block
       integer(4) :: i, ierr
       real(8), dimension(:), allocatable :: rhoeq, deriv
       real(8), dimension(:,:), allocatable :: u

       allocate(rhoeq, deriv, ProfR, source=array_init_NR)
       do NR = 0, NRMAX
          rhoeq(NR) = NR / real(NRMAX,8) * rhob
       end do

       if( MDFIXT == 0 ) then
          ! Keep temperature profiles in ProfRS temporarily
          allocate(ProfRS, mold=array_init_NRNS)
          ProfRS(:,1) = X(:,LQe5) / X(:,LQe1)
          ProfRS(:,2) = X(:,LQi5) / X(:,LQi1)
          ProfRS(:,3) = X(:,LQz5) / X(:,LQz1)
       end if

       allocate(u(0:NRMAX,4))
       ProfR(:) = X(:,LQe1)
       call spl1d(rho,ProfR,deriv,u,NRMAX+1,0,ierr)
       do NR = 0, NRMAX
          call spl1df(rhoeq(NR),ProfR(NR),rho,u,NRMAX+1,ierr)
       end do
       call spl1d(rhoeq,ProfR,deriv,u,NRMAX+1,0,ierr)
       do NR = 0, NRMAX
          call spl1df(rho(NR),X(NR,LQe1),rhoeq,u,NRMAX+1,ierr)
       end do
       X(:,LQi1) = faci * X(:,LQe1)
       X(:,LQz1) = facz * X(:,LQe1)

       if( MDFIXT == 0 ) then
          X(:,LQe5) = ProfRS(:,1) * X(:,LQe1)
          X(:,LQi5) = ProfRS(:,2) * X(:,LQi1)
          X(:,LQz5) = ProfRS(:,3) * X(:,LQz1)
          deallocate(ProfRS)
       end if

       deallocate(rhoeq,deriv,ProfR,u)
     end block
  end if

  ! ********************************************************************
  !      Effective charge
  ! ********************************************************************

  Zeff(:) = (  achg(2)**2 * X(:,LQi1) &
       &     + achg(3)**2 * X(:,LQz1) ) / X(:,LQe1)

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
        sum_int = sum_int + intg_vol_p(aat,NR) * fipol(NRMAX) / (4.d0 * Pisq)
        X(NR,LQm5) = sum_int / rMU0 ! PsitV / rMU0
     end do

     ! === Poloidal magnetic field ===

     allocate(Profsdt, dProfsdt, mold=vv)
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
     sdt(0) = FCTR4pt(rho(1),rho(2),rho(3),sdt(1),sdt(2),sdt(3))
     Profsdt(0) = 1.d0
     dProfsdt(0) = PROFJ+1
  end if

!  if(FSHL /= 0.d0) then
!     BthV(0) = 0.d0
!     Q(0) = Q0
!     do NR = 1, NRMAX
!        RL = rpt(NR)
!        Q(NR) = (Q0 - QA) * (1.d0 - (RL / RA)**2) + QA
!        BthV(NR) = BB * RL / (Q(NR) * RR)
!     end do
!  end if

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

     do NR = 0, NRMAX
        X(NR,LQe7) =-AJphVRL(NR) / (AEE * 1.D20)
        X(NR,LQe4) = X(NR,LQe7) / aat(NR) ! approx
        AJOH(NR)   = X(NR,LQe7) / ait(NR)
        X(NR,LQm4) = PsiV(NR)
     end do

  else
     
     ! === File input ===

     allocate(AJPHL, AJphVRL, source=array_init_NR)
     ifile = detect_datatype('LQe4')
     if(ifile == 0) then ! No toroidal current (LQe4) data in a structured type
        if(MDINTC <= -1) then ! Current density read from file
           do NR = 0, NRA
              call initprof_input(NR,4,AJPHL(NR)) ! electron current read from an external file
           end do
           AJFCT = rIPs * 1.D6 / intg_area(AJPHL)
           ! Artificially extrapolate a current density in the SOL for numerical stability
           do NR = NRA+1, NRMAX
              AJPHL(NR) = AJPHL(NRA) * EXP(- (rho(NR) - 1.d0) / (0.5d0 * rLn))
           end do

           if(MDINTC == -2) then ! Smoothing current density
              allocate(ProfR, source=AJPHL)
              NR_smt = NRA - NR_smt_start ! smoothing data only in the edge region
              do NR = NR_smt, NRMAX
                 AJPHL(NR) = moving_average(NR,ProfR,NRMAX)
              end do
              deallocate(ProfR)
           end if

           do NR = 0, NRMAX
              AJPHL(NR)   = AJFCT * AJPHL(NR)
              AJphVRL(NR) = AJPHL(NR) * ait(NR)
              X(NR,LQe4)  =-AJphVRL(NR) / (AEE * 1.D20) * rrt(NR) ! approx
              X(NR,LQe7)  =-AJphVRL(NR) / (AEE * 1.D20)
              AJOH(NR)    = AJPHL(NR)
           end do

           BthV(0) = 0.d0
           sum_int = 0.d0
           do NR = 1, NRMAX
              sum_int  = sum_int + intg_vol_p(AJphVRL,NR)
              BthV(NR) = rMUb1 * sum_int / sqrt(ckt(NR))
              sdt(NR)  = rMUb1 * sum_int / ckt(NR)
           end do
           sdt(0) = FCTR4pt(rho(1),rho(2),rho(3),sdt(1),sdt(2),sdt(3))
        else ! (MDINTC == 0); Current density constructed
           do NR = 0, NRMAX
              AJphVRL(NR) = BCLQm3 / (rMUb1 * vlt(nra)) * dProfsdt(NR) ! <j_zeta/R>
              AJPHL(NR)   = AJphVRL(NR) / ait(NR) ! <j_zeta/R>/<1/R>
              X(NR,LQi7)  = (Uiph0 * ait(NR)) * X(NR,LQi1) * (dProfsdt(NR) / dProfsdt(0))
              X(NR,LQi4)  = X(NR,LQi7) * rrt(NR) ! n <R u_zeta> = n <u_zeta/R>/<R^2>
              X(NR,LQe7)  =-AJphVRL(NR) / (AEE * 1.D20) + achg(2) * X(NR,LQi7)
              X(NR,LQe4)  = X(NR,LQe7) * rrt(NR) ! n <R u_zeta> = n <u_zeta/R>/<R^2>
              AJOH(NR)    = AJPHL(NR)
!              if(FSHL == 0.d0) then
!                 AJphVRL(NR) = 0.d0
!                 X(NR,LQe4)  = 0.d0
!                 X(NR,LQe7)  = 0.d0
!                 AJOH(NR)    = 0.d0
!              end if
           end do
        end if

     else ! Detect toroidal current data in a structured type
        call inexpolate(infiles(ifile)%nol,infiles(ifile)%r,infiles(ifile)%data,NRMAX,RHO,5,AJPHL)
        AJFCT = rIPs * 1.D6 / intg_area(AJPHL)
        do NR = 0, NRMAX
           AJPHL(NR)   = AJFCT * AJPHL(NR)
           AJphVRL(NR) = AJPHL(NR) * ait(NR)
           X(NR,LQe4)  =-AJPHL(NR) / (AEE * 1.D20) * rrt(NR) ! approx
           X(NR,LQe7)  =-AJPHL(NR) / (AEE * 1.D20)
           AJOH(NR)    = AJPHL(NR)
        end do

        BthV(0) = 0.d0
        sum_int = 0.d0
        do NR = 1, NRMAX
           sum_int  = sum_int + intg_vol_p(AJphVRL,NR)
           BthV(NR) = rMUb1 * sum_int / sqrt(ckt(NR))
           sdt(NR)  = rMUb1 * sum_int / ckt(NR)
        end do
        sdt(0) = FCTR4pt(rho(1),rho(2),rho(3),sdt(1),sdt(2),sdt(3))
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
  if(MDINTN1 < 0 .or. MDINTT1 < 0 .or. abs(MDINTC) /= 0) call initprof_input(idx = 0)

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

!  if(FSHL == 0.d0) then
!     ! Integrate 1 / (r * rMU0) * d/dr (r * BthV) to obtain AJV
!     AJV(:) = BB / (RR * rMU0) * 2.d0 * Q0 / Q(:)**2
!  end if

  ! ********************************************************************
  !      Toroidal electric field for initial NCLASS calculation
  ! ********************************************************************

  do NR = 0, NRMAX
     Var(NR,1)%n = X(NR,LQe1)
     Var(NR,2)%n = X(NR,LQi1)
     Var(NR,3)%n = X(NR,LQz1)
     if(MDFIXT == 0) then
        Var(NR,1)%T = X(NR,LQe5) / X(NR,LQe1)
        Var(NR,2)%T = X(NR,LQi5) / X(NR,LQi1)
        Var(NR,3)%T = X(NR,LQz5) / X(NR,LQz1)
     else
        Var(NR,1)%T = X(NR,LQe5)
        Var(NR,2)%T = X(NR,LQi5)
        Var(NR,3)%T = X(NR,LQz5)
     end if
     ! Inverse aspect ratio
     EpsL = epst(NR)
     ! Trapped particle fraction
     FTL  = 1.46d0 * SQRT(EpsL) - 0.46d0 * EpsL**1.5d0
     ! Estimating parallel resistivity
     call sauter(Var(NR,1)%n,Var(NR,1)%T,dum,dum,Var(NR,2)%n,Var(NR,2)%T,dum,dum, &
       &         Q(NR),sdt(NR),fipol(NR),EpsL,RR,achg(2),Zeff(NR),ftl, &
       &         rlnLei_IN=coulog(Zeff(NR),Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amas(1),achg(1),amas(2),achg(2)), &
       &         rlnLii_IN=coulog(Zeff(NR),Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amas(2),achg(2),amas(2),achg(2)), &
       &         BJBS=dum,ETA=etanc,ETAS=etaspz)
     if(abs(FSNC(1)) > 0.d0) then
        eta(NR) = etanc
     else
        eta(NR) = etaspz
     end if
     ! R Et = PsitdotV
     !    <B E//> = eta// <B j//> and neglecting poloidal components yields <Bt Et> = eta// <Bt j//>
     !    <Bt Et> = I<1/R^2> R E_t and <Bt jt> = I<jt/R> give R Et = eta// <jt/R> / <1/R^2>
     X(NR,LQm3) = eta(NR) * AJphVRL(NR) / aat(NR)
     if(X(NR,LQm3) == 0.d0) X(NR,LQm3) = 1.D-4
  end do
!  if(FSHL == 0.d0) X(:,LQm3) = 0.d0
  if( allocated(AJphVRL) ) deallocate(AJphVRL)

  ! ********************************************************************
  !      Electrostatic potential, Phi
  ! ********************************************************************

  !  dPhi/dV is computed by summing LQ*2 over species

  block
    integer(4) :: i, j, NE
    real(8) :: phipre_NRMAX
    real(8), dimension(:), allocatable :: sumden, sumdpre, sumupara, sumuph, kernel, phipre
    real(8), dimension(:,:), allocatable :: dpre

    allocate(sumden, sumdpre, sumupara, sumuph, kernel, source=array_init_NR)
    allocate(dpre, source=array_init_NRNS)

    do i = 1, NSM
       dpre(:,i) = dfdx(vv,Var(:,i)%n*Var(:,i)%T,NRMAX,0) * rKilo ! dp/dV
    end do

    do NR = 0, NRMAX
       do i = 1, NSM
          j = (i - 1) * 8 ! 8 indicates the number of equations for each species.
          sumden(NR)   = sumden(NR)   + achg(i) * Var(NR,i)%n / amas(i)
          sumdpre(NR)  = sumdpre(NR)  + dpre(NR,i) / amas(i)
          sumupara(NR) = sumupara(NR) + achg(i) * Var(NR,i)%n / amas(i) * X(NR,LQe3+j)
          sumuph(NR)   = sumuph(NR)   + achg(i) / amas(i) * X(NR,LQe4+j)
       end do
    end do

    !  Evaluate bri
    !    Note: qhatsq(0) remains indefinite temporarily because it is not used here.

    !  Square of the poloidal magnetic field: <B_p^2> = ckt * (dpsi/dV)^2
    Bpsq(:) = ckt(:) * sdt(:)*sdt(:)
    bri(0) = 0.d0
    do NR = 1, NRMAX
       !  qhat square: q^^2 = I^2/(2<B_p^2>)(<1/R^2>-1/<R^2>)
       qhatsq(NR) = 0.5d0 * fipol(NR)*fipol(NR) / Bpsq(NR) * (aat(NR) - 1.d0 / rrt(NR))
       !  Metric coefficient: <B^2><R^2>-I^2
       bri(NR) = rrt(NR) * Bpsq(NR) * (1.d0 + 2.d0 * qhatsq(NR))
    end do

    !  kernel = dPhi/dV

    kernel(1:) = (- bri(1:) * sumdpre(1:) + fipol(1:) * sdt(1:) * sumupara(1:) &
         &        - bbt(1:) * sdt(1:) * sumuph(1:)) / (bri(1:) * sumden(1:))
    kernel(0) = FCTR4pt(rho(1),rho(2),rho(3),kernel(1),kernel(2),kernel(3))

    !  Integrate dPhi/dV w.r.t. V from the boundary where Phi=0 to the axis

    sum_int = 0.d0
    X(NRMAX,LQm1) = 0.d0
    do NE = NEMAX, 1, -1
       NR = NE - 1
       sum_int = sum_int + intg_vol_p(kernel,NE)
       X(NR,LQm1) = sum_int ! Phi
    end do

    ! Diamagnetic particle flow
    !   Note: The sign of kernel(=dPhi/dV) have to flip so as to be consistent with
    !         the V-derivative of X(:,LQm1).

    do NR = 1, NRMAX
       X(NR,LQe8) = - fipol(NR) &
            & * ( dpre(NR,1) / ( achg(1) * Var(NR,1)%n ) - kernel(NR) ) / sdt(NR)
       X(NR,LQi8) = - fipol(NR) &
            & * ( dpre(NR,2) / ( achg(2) * Var(NR,2)%n ) - kernel(NR) ) / sdt(NR)
       X(NR,LQz8) = - fipol(NR) &
            & * ( dpre(NR,3) / ( achg(3) * Var(NR,3)%n ) - kernel(NR) ) / sdt(NR)
    end do
    X(0,LQe8) = FCTR4pt(rho(1),rho(2),rho(3),X(1,LQe8),X(2,LQe8),X(3,LQe8))
    X(0,LQi8) = FCTR4pt(rho(1),rho(2),rho(3),X(1,LQi8),X(2,LQi8),X(3,LQi8))
    X(0,LQz8) = FCTR4pt(rho(1),rho(2),rho(3),X(1,LQz8),X(2,LQz8),X(3,LQz8))

    deallocate(dpre, sumden, sumdpre, sumupara, sumuph, kernel)
  end block

  ! ********************************************************************
  !      Coefficients for Pfirsch-Schluter viscosity (nccoe)
  ! ********************************************************************

  call wrap_eqneo

!!$  if( ieqread >= 2 ) then
!!$!     if( MDLNEOL == 2 ) call eqneo
!!$     call eqneo
!!$  else
!!$     do NR = 1, NRMAX
!!$        gamneo(NR) = 4.d0 * Pisq * sdt(NR) / bbrt(NR) ! = bthco(NR) / bbrt(NR)
!!$        mxneo(NR) = 3
!!$        epsl = epst(NR)
!!$        coefmneo = 1.d0 - epsl**2
!!$        do i = 1, mxneo(NR)
!!$           fmneo(i,NR) = real(i,8)* ( (1.d0-sqrt(coefmneo))/epsl)**(2*i) &
!!$                &                 * (1.d0+real(i,8)*sqrt(coefmneo)) &
!!$                &                 / (coefmneo*sqrt(coefmneo)*(q(NR)*RR)**2)
!!$        end do
!!$        fmneo(mxneo(NR)+1:10,NR) = 0.d0
!!$     end do
!!$  end if
!!$  ! Even at axis, fmneo=0 should be avoided because it would cause K_PS = 0.
!!$  NR = 0
!!$     gamneo(NR) = 4.d0 * Pisq * sdt(NR) / bbrt(NR) ! = bthco(NR) / bbrt(NR)
!!$     mxneo(NR) = 3
!!$!     fmneo(1:10,NR) = 0.d0
!!$     epsl = smallvalue
!!$     coefmneo = 1.d0 - epsl**2
!!$     do i = 1, mxneo(NR)
!!$        fmneo(i,NR) = real(i,8)* ( (1.d0-sqrt(coefmneo))/epsl)**(2*i) &
!!$             &                 * (1.d0+real(i,8)*sqrt(coefmneo)) &
!!$             &                 / (coefmneo*sqrt(coefmneo)*(q(NR)*RR)**2)
!!$     end do
!!$     fmneo(mxneo(NR)+1:10,NR) = 0.d0

  ! ********************************************************************
  !      Miscellaneous
  ! ********************************************************************

  !   NBI total input power (MW)
  PNBH = PNBHP + PNBHT1 + PNBHT2

  !     Averaged beam injection energy
  Eb =  (esps(1) + 0.5d0 * esps(2) + esps(3) / 3.d0) * Ebmax

  T_TX=0.d0
  NGT=-1
  NGR=-1
  NGVV=-1
  rIP=rIPs

  !  Check whether (m=0, n>0) Fourier component exists or not.   miki_m 10-09-07
  UHphSwitch = 0
  do NHFM = 1, NHFMmx 
     if (HPN(NHFM,1) == 0 .and. HPN(NHFM,2) > 0) UHphSwitch = 1
  end do

  !  Define physical variables from X

  call TXCALV(X,0) ! Set variables as well as pres0 and ErV0
  !  --- Fixed Er to keep it constant during iterations ---
  ErV_FIX(:) = ErV(:)

  !  Calculate various physical quantities

  call TXCALC(0)

  !  Calculate global quantities for storing and showing initial status

  call TXGLOB

end subroutine TXPROF

!*****************************************************************************************

module tx_parameter_control
  use tx_commons
  use tx_graphic, only : MODEG, MODEGL, NGRSTP, NGTSTP, NGVSTP, gDIV
  implicit none
  public
  namelist /TX/ &
       & RA,rhob,rhoaccum,RR,BB,rbvt, &
       & amas,achg,Zeffin, &
       & PN0,PNa,PTe0,PTea,PTi0,PTia,PTz0,PTza,PROFJ,PROFN1,PROFN2,PROFT1,PROFT2,Uiph0, &
       & Dfs0,VWpch0,rMus0,WPM0,Chis0,ChiNC, &
       & FSDFIX,FSANOM,FSCBAL,FSCBKP,FSCBEL,FSCBSH,FSBOHM,FSPCL,FSVAHL,MDANOM,RhoETB, &
       & FSMPCH,FSTPTM,FSPARV,PROFD,PROFD1,PROFD2,PROFDB,PROFM,PROFM1,PROFMB,PROFC,PROFC1,PROFCB, &
       & FSCX,FSLC,FSRP,FSNF,FSADV,FSADVB,FSUG,FSNC,FSNCB,FSLP,FSLPB,FSLTs,FSION, &
       & FSNCOL,FSD01,FSD02,FSD03,FSD0z,MDLC,rLn,rLT, &
       & Ebmax,esps,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBHP,PNBHT1,PNBHT2,PNBCD,PNBMPD,PNBPTC, &
       & rNRFe,RRFew,RRFe0,PRFHe,Tqt0,Tqp0,rNRFi,RRFiw,RRFi0,PRFHi, &
       & PN0s,V0,rGamm0,rGASPF,PN0zs,V0z,rGamm0z,rGASPFz,PNsDIV,PTsDIV, &
       & NTCOIL,DltRPn,kappa,m_pol,n_tor, &
       & DT,EPS,ICMAX,ADV,tiny_cap,CMESH0,CMESH,WMESH0,WMESH, &
       & NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,ieqread,irktrc, &
       & DelRho,DelN, &
!       & DMAG0,RMAGMN,RMAGMX,EpsH,NCph,NCth, &
       & DMAG0,RMAGMN,RMAGMX, &
       & rG1,FSHL,Q0,QA, &
       & rIPs,rIPe, &
       & MODEG,gDIV,MODEAV,MODEGL,MDLPCK,MODECV,oldmix,iSUPG3,iSUPG6,iSUPG8,SUPGstab,iprestab, &
       & MDFIXT,MDBEAM,MDOSQZ,MDOSQZN,MDLETA,MDLNEO,MDBSETA,MDITSN,MDITST,MDINTN,MDINTT,MDINTC,MDLETB, &
       & iPoyntpol,iPoynttor,IDIAG,midbg,IGBDF,ISMTHD,MDLNBD,imodel_neo, & ! 09/06/17~ miki_m
       & EpsHM, HPN  ! 10/08/06 miki_m
  private :: TXPLST

contains
!***************************************************************
!
!   Change input parameters
!
!***************************************************************

  subroutine TXPARM(KID)

    integer(4) :: IST
    character(len=*)  :: KID

    do 
       write(6,*) '# INPUT &TX :'
       write(6,*) '  INPUT PARAMETERS MUST BE BRACKETED BETWEEN &TX AND &END; CTRL-D FOR EXIT.'
       read(5,TX,IOSTAT=IST)
       if(IST > 0) then
          call TXPLST
          cycle
       else if(IST < 0) then
          KID='Q'
          exit
       else
          KID=' '
          exit
       end if
    end do
    IERR=0

  end subroutine TXPARM

  subroutine TXPARL(KLINE)

    integer(4) :: IST
    character(len=*) :: KLINE
    character(len=90) :: KNAME

    KNAME=' &TX '//KLINE//' &END'
    write(7,'(A90)') KNAME
    rewind(7)
    read(7,TX,IOSTAT=IST)
    if(IST /= 0) then
       call TXPLST
    else
       write(6,'(A)') ' ## PARM INPUT ACCEPTED.'
    end if
    rewind(7)
    IERR=0

  end subroutine TXPARL

  subroutine TXPARF(KPNAME)

    use libchar, only : ktrim
    integer(4) :: IST, KL
    logical :: LEX
    character(len=*) :: KPNAME

    inquire(file=KPNAME,exist=LEX)
    if(.not. LEX) return

    open(25,file=KPNAME,iostat=IST,status='old')
    if(IST > 0) then
       write(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
       return
    end if
    read(25,TX,IOSTAT=IST)
    if(IST > 0) then
       write(6,*) 'XX PARM FILE READ ERROR'
       return
    else if(IST < 0) then
       write(6,*) 'XX PARM FILE EOF ERROR'
       return
    end if
    call KTRIM(KPNAME,KL)
    write(6,*) &
         &     '## FILE (',KPNAME(1:KL),') IS ASSIGNED FOR PARM INPUT'
    IERR=0
    close(25)

  end subroutine TXPARF

  subroutine TXPARM_CHECK

    integer(4) :: idx

    idx = 1
    do 
       ! /// idx = 1 - 10 ///
       ! System integers
       if(NRMAX < 0) then ; exit ; else ; idx = idx + 1 ; end if
       if(NTMAX < 0) then ; exit ; else ; idx = idx + 1 ; end if
       if(NTSTEP < 0 .or. NGRSTP < 0 .or. NGTSTP < 0 .or. NGVSTP < 0) then ; exit ; else
          idx = idx + 1 ; end if
       ! Physical variables
       if(rhob <= 1.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(rhoaccum > rhob) then ; exit ; else ; idx = idx + 1 ; end if
       if(RR <= rhob * ra) then ; exit ; else ; idx = idx + 1 ; end if
       if(rIPs < 0.d0 .or. rIPe < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(minval(amas) < 0.d0) then ; exit ; else ; idx = idx + 1 ; endif
       if(achg(1) > 0.d0 .or. minval(achg(2:3)) < 0.d0 .or. Zeffin < 1.d0 .or. Zeffin >= maxval(achg(2:3))) then
          exit ; else ; idx = idx + 1 ; end if
       if(PN0 < 0.d0 .or. PNa < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       ! /// idx = 11 - 20 ///
       if(PTe0 < 0.d0 .or. PTea < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(PTi0 < 0.d0 .or. PTia < 0.d0 .or. PTz0 < 0.d0 .or. PTza < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(minval(Dfs0) < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(minval(rMus0) < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(minval(Chis0) < 0.d0 .or. ChiNC < 0.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(minval(FSDFIX) < 0.d0 .or. minval(FSANOM) < 0.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(minval(RhoETB) < 0.d0 .or. maxval(RhoETB) > 1.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(minval(FSMPCH) < 0.d0 .or. maxval(FSMPCH) > 1.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(minval(FSTPTM) < 0.d0 .or. maxval(FSTPTM) > 2.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(minval(FSPARV) < 0.d0 .or. maxval(FSPARV) > 1.d0) then ; exit ; else
          idx = idx + 1 ; end if
       ! /// idx = 21 - 30 ///
       if(FSCBAL < 0.d0 .or. FSCBKP < 0.d0 .or. FSCBEL < 0.d0 .or. FSCBSH < 0.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(FSBOHM < 0.d0 .or. minval(FSPCL) < 0.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(FSLC < 0.d0 .or. FSRP < 0.d0 .or. minval(FSNC) < 0.d0 .or. minval(FSNCB) < 0.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(FSHL < 0.d0 .or. FSNF < 0.d0 .or. FSADV < 0.d0 .or. FSADVB < 0.d0 .or. FSUG < 0.d0) then
          exit ; else ; idx = idx + 1 ; end if
       if(FSLP < 0.d0 .or. minval(FSLTs) < 0.d0 .or. FSLPB < 0.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(FSION < 0.d0 .or. FSNCOL < 0.d0)  then ; exit ; else ; idx = idx + 1 ; end if
       if(FSD01 < 0.d0 .or. FSD02 < 0.d0 .or. FSD03 < 0.d0 .or. FSD0z < 0.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(MDLC /= 1 .and. MDLC /= 2) then ; exit ; else ; idx = idx + 1 ; end if
       if(rG1 < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(Ebmax < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       ! /// idx = 31 - 40 ///
       if(sum(esps) /= 1.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(PNBHP < 0.d0 .or. PNBHT1 < 0.d0 .or. PNBHT2 < 0.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(abs(PNBCD) > 1.d0 .or. abs(PNBMPD) > 1.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(PNBPTC < 0.d0 .or. PNBPTC > 1.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(RNBP0 > rhob .or. RNBP0 < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(RNBT10 > rhob .or. RNBT10 < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(RNBT20 > rhob .or. RNBT20 < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(rNRFe < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(rNRFi < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(RRFe0 > rhob .or. RRFe0 < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       ! /// idx = 41 - 50 ///
       if(RRFi0 > rhob .or. RRFi0 < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(PN0s < 0.d0 .or. V0 < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(rGamm0 < 0.d0 .or. rGamm0 > 1.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(rGASPF < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(minval(PNsDIV) < 0.d0 .or. minval(PTsDIV) < 0.d0) then ; exit ; else
          idx = idx + 1 ; end if
       if(NTCOIL <= 0 .or. DltRPn < 0.d0 .or. DltRPn > 1.d0 .or. kappa < 0.d0) then
          exit ; else ; idx = idx + 1 ; end if
       if(DT < 0.d0 .or. EPS == 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(ICMAX < 0) then ; exit ; else ; idx = idx + 1 ; end if
       if(ADV < 0.d0 .or. ADV > 1.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(tiny_cap < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       ! /// idx = 51 - 58 ///
       if(oldmix < 0.d0 .or. oldmix > 1.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(CMESH0 < 0.d0 .or. CMESH < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(WMESH0 < 0.d0 .or. WMESH < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(MDLNEO /= 1 .and. MDLNEO /= 2 .and. MDLNEO /= 3) then ; exit ; else ; idx = idx + 1 ; end if
       if(MDBEAM <  0 .or. MDBEAM > 2) then ; exit ; else ; idx = idx + 1 ; end if
       if(PN0zs < 0.d0 .or. V0z < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(rGamm0z < 0.d0 .or. rGamm0z > 1.d0) then ; exit ; else ; idx = idx + 1 ; end if
       if(rGASPFz < 0.d0) then ; exit ; else ; idx = idx + 1 ; end if
       return
    end do

    write(6,'(A,I3)') 'XX INPUT ERROR: Please check consistency and/or sanity of INPUT PARAMETERS. idx =',idx
    stop
    
  end subroutine TXPARM_CHECK

  !***** INPUT PARAMETER LIST *****

  subroutine TXPLST

    write(6,601)
    return

601 format(' ','# &TX : RA,rhob,rhoaccum,RR,BB,rbvt,amas,achg,Zeffin,'/ &
         &       ' ',8X,'PN0,PNa,PTe0,PTea,PTi0,PTia,PTz0,PTza,PROFJ,PROFN1,PROFN2,PROFT1,PROFT2,Uiph0,'/ &
         &       ' ',8X,'Dfs0,VWpch0,rMus0,WPM0,Chis0,ChiNC,'/ &
         &       ' ',8X,'FSDFIX,FSANOM,FSCBAL,FSCBKP,FSCBEL,FSCBSH,FSBOHM,FSPCL,FSVAHL,'/ &
         &       ' ',8X,'MDANOM,RhoETB,FSMPCH,FSTPTM,FSPARV'/ &
         &       ' ',8X,'PROFD,PROFD1,PROFD2,PROFDB,PROFM,PROFM1,PROFMB,PROFC,PROFC1,PROFCB,'/ &
         &       ' ',8X,'FSCX,FSLC,FSRP,FSNF,FSADV,FSADVB,FSUG,FSNC,FSNCB,FSLP,FSLTs,FSLPB,FSION,'/ &
         &       ' ',8X,'FSNCOL,FSD01,FSD02,FSD03,FSD0z,MDLC,rLn,rLT,'/ &
         &       ' ',8X,'Ebmax,esps,RNBP,RNBP0,RNBT1,RNBT2,RNBT10,RNBT20,PNBHP,PNBHT1,PNBHT2,'/ &
         &       ' ',8X,'PNBCD,PNBMPD,PNBPTC,rNRFe,RRFew,RRFe0,PRFHe,Tqt0,Tqp0,rNRFe,RRFew,RRFe0,PRFHe,'/&
         &       ' ',8X,'PN0s,V0,rGamm0,rGASPF,PN0zs,V0z,rGamm0z,rGASPFz,PNsDIV,PTsDIV,'/ &
         &       ' ',8X,'NTCOIL,DltRPn,kappa,m_pol,n_tor,'/ &
         &       ' ',8X,'DT,EPS,ICMAX,ADV,tiny_cap,CMESH0,CMESH,WMESH0,WMESH,'/ &
         &       ' ',8X,'NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,ieqread,irktrc,'/ &
         &       ' ',8X,'DelRho,DelN,'/ &
!         &       ' ',8X,'Dmag0,RMAGMN,RMAGMX,EpsH,NCph,NCth,'/ & 
         &       ' ',8X,'Dmag0,RMAGMN,RMAGMX,EpsHM,HPN,'/ &   ! 10/08/06 miki_m
         &       ' ',8X,'rG1,FSHL,Q0,QA,'/ &
         &       ' ',8X,'rIPs,rIPe,'/ &
         &       ' ',8X,'MODEG,gDIV,MODEAV,MODEGL,MDLPCK,MODECV,oldmix,iSUPG3,iSUPG6,iSUPG8,SUPGstab,iprestab,'/ &
         &       ' ',8X,'MDFIXT,MDBEAM,MDOSQZ,MDOSQZN,MDLETA,MDLNEO,MDBSETA,MDITSN,MDITST,MDINTN,MDINTT,MDLETB,' / & 
         &       ' ',8X,'IDIAG,midbg,IGBDF,ISMTHD,MDLNBD,imodel_neo')
  end subroutine TXPLST

!***************************************************************
!
!   View input parameters
!
!***************************************************************

  subroutine TXVIEW

    write(6,'((1X,A10," =",ES11.4,2(2X,A10," =",ES11.4)))') &
         &   'RA        ', RA       , 'RHOB      ', RHOB     ,  &
         &   'RR        ', RR       , 'BB        ', BB       ,  &
         &   'RAVL      ', RAVL     , 'RBVL      ', RBVL     ,  &
         &   'rbvt      ', rbvt     , &
         &   'amas(2)   ', amas(2)  , 'achg(2)   ', achg(2)  ,  &
         &   'amas(3)   ', amas(3)  , 'achg(3)   ', achg(3)  ,  &
         &   'PN0       ', PN0      , 'PNa       ', PNa      ,  &
         &   'PTe0      ', PTe0     , 'PTea      ', PTea     ,  &
         &   'PTi0      ', PTi0     , 'PTia      ', PTia     ,  &
         &   'Zeffin    ', Zeffin   , 'rIP       ', rIP      ,  &
         &   'rIPs      ', rIPs     , 'rIPe      ', rIPe     ,  &
         &   'PROFJ     ', PROFJ    , 'PROFN1    ', PROFN1   ,  &
         &   'PROFN2    ', PROFN2   , 'PROFT1    ', PROFT1   ,  &
         &   'PROFT2    ', PROFT2   , 'Uiph0     ', Uiph0    ,  &
         &   'CMESH0    ', CMESH0   , &
         &   'WMESH0    ', WMESH0   , 'CMESH     ', CMESH    ,  &
         &   'WMESH     ', WMESH    , 'ADV       ', ADV      ,  &
         &   'Dfs0(1)   ', Dfs0(1)  , 'Dfs0(2)   ', Dfs0(2)  ,  &
         &   'Dfs0(3)   ', Dfs0(3)  , 'rMus0(1)  ', rMus0(1) ,  &
         &   'rMus0(2)  ', rMus0(2) , 'rMus0(3)  ', rMus0(3) ,  &
         &   'ChiNC     ', ChiNC    , 'Chis0(1)  ', Chis0(1) ,  &
         &   'Chis0(2)  ', Chis0(2) , 'Chis0(3)  ', Chis0(3) ,  &
         &   'FSMPCH(1) ', FSMPCH(1), 'FSMPCH(2) ', FSMPCH(2),  &
         &   'FSMPCH(3) ', FSMPCH(3), 'FSTPTM(1) ', FSTPTM(1),  &
         &   'FSTPTM(2) ', FSTPTM(2), 'FSTPTM(3) ', FSTPTM(3),  &
         &   'FSPARV(1) ', FSPARV(1), 'FSPARV(2) ', FSPARV(2),  &
         &   'FSPARV(3) ', FSPARV(3),  &
         &   'VWpch0    ', VWpch0   , 'WPM0      ', WPM0     ,  &
         &   'PROFD     ', PROFD    , 'PROFD1    ', PROFD1   ,  &
         &   'PROFD2    ', PROFD2   , 'PROFDB    ', PROFDB   ,  &
         &   'PROFM     ', PROFM    , 'PROFM1    ', PROFM1   ,  &
         &   'PROFMB    ', PROFMB   , 'PROFC     ', PROFC    ,  &
         &   'PROFC1    ', PROFC1   , 'PROFCB    ', PROFCB   ,  &
         &   'FSDFIX(1) ', FSDFIX(1), 'FSDFIX(2) ', FSDFIX(2),  &
         &   'FSDFIX(3) ', FSDFIX(3), 'FSANOM(1) ', FSANOM(1),  &
         &   'FSANOM(2) ', FSANOM(2), 'FSANOM(3) ', FSANOM(3),  &
         &   'RhoETB(1) ', RhoETB(1), 'RhoETB(2) ', RhoETB(2),  &
         &   'RhoETB(3) ', RhoETB(3),  &
         &   'FSCBAL    ', FSCBAL   , 'FSCBKP    ', FSCBKP   ,  &
         &   'FSCBEL    ', FSCBEL   , 'FSCBSH    ', FSCBSH   ,  &
         &   'FSBOHM    ', FSBOHM   , 'FSPCL(1)  ', FSPCL(1) ,  &
         &   'FSPCL(2)  ', FSPCL(2) , 'FSPCL(3)  ', FSPCL(3) ,  &
         &   'FSVAHL    ', FSVAHL   , 'FSCX      ', FSCX     ,  &
         &   'FSLC      ', FSLC     , 'FSRP      ', FSRP     ,  &
         &   'FSNF      ', FSNF     , 'FSNC(1)   ', FSNC(1)  ,  &
         &   'FSNC(2)   ', FSNC(2)  , 'FSNCB(1)  ', FSNCB(1) ,  &
         &   'FSNCB(2)  ', FSNCB(2) , 'FSADV     ', FSADV    ,  &
         &   'FSADVB    ', FSADVB   , 'FSUG      ', FSUG     ,  &
         &   'FSLP      ', FSLP     , 'FSLTs(1)  ', FSLTs(1) ,  &
         &   'FSLTs(2)  ', FSLTs(2) , 'FSLTs(3)  ', FSLTs(3) ,  &
         &   'FSLPB     ', FSLPB    , 'FSION     ', FSION    ,  &
         &   'FSNCOL    ', FSNCOL   , 'FSD01     ', FSD01    ,  &
         &   'FSD02     ', FSD02    , 'FSD03     ', FSD03    ,  &
         &   'FSD0z     ', FSD0z    ,  &
         &   'rLn       ', rLn      , 'rLT       ', rLT      ,  &
         &   'Ebmax     ', Ebmax    , 'esps(1)   ', esps(1)  ,  &
         &   'esps(2)   ', esps(2)  , 'esps(3)   ', esps(3)  ,  &
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
         &   'rGASPF    ', rGASPF   , 'rGamm0z   ', rGamm0z  ,  &
         &   'V0z       ', V0z      , 'rGASPFz   ', rGASPFz  ,  &
         &   'PNeDIV    ', PNsDIV(1), 'PNiDIV    ', PNsDIV(2),  &
         &   'PNzDIV    ', PNsDIV(3), 'PTeDIV    ', PTsDIV(1),  &
         &   'PTiDIV    ', PTsDIV(2), 'PTzDIV    ', PTsDIV(3),  &
         &   'DltRPn    ', DltRPn   , 'kappa     ', kappa    ,  &
         &   'PN0s      ', PN0s     , 'PN0zs     ', PN0zs    ,  &
         &   'EPS       ', EPS      , 'tiny_cap  ', tiny_cap ,  &
         &   'DT        ', DT       , 'rG1       ', rG1      ,  &
         &   'DMAG0     ', DMAG0    , 'RMAGMN    ', RMAGMN   ,  &
         &   'RMAGMX    ', RMAGMX   ,  &
         &   'FSHL      ', FSHL     ,  &
         &   'Q0        ', Q0       , 'QA        ', QA       ,  &
         &   'SUPGstb   ', SUPGstab,  'oldmix    ', oldmix   !, &
!!$         &   'EpsHM(1,0)', EpsHM(1,0),'EpsHM(1,1)', EpsHM(1,1), &
!!$         &   'EpsHM(1,2)', EpsHM(1,2),'EpsHM(1,3)', EpsHM(1,3), &
!!$         &   'EpsHM(2,0)', EpsHM(2,0),'EpsHM(2,1)', EpsHM(2,1), &
!!$         &   'EpsHM(2,2)', EpsHM(2,2),'EpsHM(2,3)', EpsHM(2,3), &
!!$         &   'EpsHM(3,0)', EpsHM(3,0),'EpsHM(3,1)', EpsHM(3,1), &
!!$         &   'EpsHM(3,2)', EpsHM(3,2),'EpsHM(3,3)', EpsHM(3,3), &
!!$         &   'EpsHM(4,0)', EpsHM(4,0),'EpsHM(4,1)', EpsHM(4,1), &
!!$         &   'EpsHM(4,2)', EpsHM(4,2),'EpsHM(4,3)', EpsHM(4,3)
    write(6,'((" ",A10," =",I5,2(8X,A10," =",I5)))') &
         &   'NRMAX     ', NRMAX    , &
         &   'NTMAX     ', NTMAX    , 'NTSTEP    ', NTSTEP   ,  &
         &   'NGRSTP    ', NGRSTP   , 'NGTSTP    ', NGTSTP   ,  &
         &   'NGVSTP    ', NGVSTP   , 'ieqread   ', ieqread  ,  &
         &   'ICMAX     ', ICMAX    , 'MODEG     ', MODEG    ,  &
         &   'MODEAV    ', MODEAV   , 'MODEGL    ', MODEGL   ,  &
         &   'MDLPCK    ', MDLPCK   , 'MODECV    ', MODECV   ,  &
         &   'IDIAG     ', IDIAG    , 'midbg(1)  ', midbg(1) ,  &
         &   'midbg(2)  ', midbg(2) ,  &
         &   'iSUPG3    ', iSUPG3   , 'iSUPG6    ', iSUPG6   ,  &
         &   'iSUPG8    ', iSUPG8   , 'IGBDF     ', IGBDF    ,  &
         &   'ISMTHD    ', ISMTHD   , 'iprestab  ', iprestab ,  &
         &   'MDFIXT    ', MDFIXT   , 'MDBEAM    ', MDBEAM   ,  &
         &   'MDOSQZ    ', MDOSQZ   , 'MDOSQZN   ', MDOSQZN  ,  &
         &   'MDLETA    ', MDLETA   , 'MDLNEO    ', MDLNEO   ,  &
         &   'MDBSETA   ', MDBSETA  , 'MDANOM    ', MDANOM   ,  &
         &   'MDITSN    ', MDITSN   , 'MDITST    ', MDITST   ,  &
         &   'MDINTN    ', MDINTN   , 'MDINTT    ', MDINTT   ,  &
         &   'MDINTC    ', MDINTC   , 'iPoyntpol ',iPoyntpol ,  &
         &   'iPoynttor ',iPoynttor ,  &
         &   'MDLETB    ', MDLETB   , 'MDLNBD    ', MDLNBD   ,  &
         &   'NTCOIL    ', NTCOIL   , 'MDLC      ', MDLC     ,  &
         &   'imodelneo1', imodel_neo(1), 'imodelneo2', imodel_neo(2),  &
         &   'imodelneo3', imodel_neo(3), 'imodelneo4', imodel_neo(4),  &
         &   'imodelneo5', imodel_neo(5), 'imodelneo6', imodel_neo(6)!,  &
!!$         &   'm_pol     ', m_pol    , 'n_tor     ', n_tor    ,  &
!!$         &   'NCph      ', NCph     , 'NCth      ', NCth     ,  &
!!$         &   'HPN(1,1)  ', HPN(1,1) , 'HPN(1,2)  ', HPN(1,2) , &
!!$         &   'HPN(2,1)  ', HPN(2,1) , 'HPN(2,2)  ', HPN(2,2) , &
!!$         &   'HPN(3,1)  ', HPN(3,1) , 'HPN(3,2)  ', HPN(3,2) , &
!!$         &   'HPN(4,1)  ', HPN(4,1) , 'HPN(4,2)  ', HPN(4,2)

  end subroutine TXVIEW

end module tx_parameter_control
