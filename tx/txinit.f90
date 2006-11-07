!     $Id$
module init_prof
  use commons
  implicit none
  public

contains

!***************************************************************
!
!   Set constants and initial parameters
!
!***************************************************************

  SUBROUTINE TXINIT

    !   ***** Configuration parameters *****

    !   Plasma minor radius (m)
    RA = 0.35D0

    !   Wall radius (m)
    RB = 0.4D0

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

    !   Initial electron temperature at r = 0 (keV)
    PTe0 = 700.D-3

    !   Initial electron temperature at r = a (keV)
    PTea =  50.D-3

    !   Initial ion temperature  at r = 0 (keV)
    PTi0 = 700.D-3

    !   Initial ion temperature  at r = a (keV)
    PTia =  50.D-3

    !   Initail current profile parameter
    PROFJ = 2.D0

    !   ***** Particle diffusivity and viscosity parameters *****

    !   Electron-driven diffusion parameter
    De0 = 0.D0

    !   Ion-driven diffusion parameter
    Di0 = 0.2D0

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
    Chie0 = 0.D0

    !   Ion thermal diffusivity parameter (Chie/D)
    !     0 for fixed temperature profile
    Chii0 = 0.D0

    !   ***** Turbulent transport control parameters *****

    !   Fixed transport coefficient parameter
    !   Suitable between 0.01 and 0.05
    FSDFIX = 0.05D0

    !   CDBM transport coefficient parameter
    !     (Current diffusive ballooning mode)
    FSCDBM = 0.D0

    !   Bohm transport coefficient parameter in SOL
    FSBOHM = 0.D0

    !   Pseud-classical transport coefficient parameter in SOL
    FSPSCL = 0.D0

    !   Diffusion coefficient profile parameter (D(r=a)/D(r=0))
    PROFD = 10.D0

    !   ***** Other transport parameters *****

    !   Charge exchange parameter
    FSCX = 1.D0

    !   Orbit loss parameter
!!!!    FSLC = 1.D0
    FSLC = 0.D0

    !   Toroidal neoclassical viscosity parameter
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

    !   Neutral diffusion factor
    FSD0 = 1.D0

    !   Factor of E x B rotation shear
    rG1 = 24.D0

    !   ***** initial parameters *****

    !   Initial Density scale length in SOL (m)
    rLn = 0.03D0

    !   Initail Temperature scale length in SOL (m)
    rLT = 0.030D0

    !   ***** Heating parameters *****

    !   NBI beam energy (keV)
    Eb = 32.D0

    !   Heating radius of NBI heating (m)
    RNB = 0.175D0

    !   NBI input power (MW)
    PNBH = 0.D0

    !   NBI current drive parameter
    PNBCD= 0.D0

    !   Refractive index of RF waves
    rNRF = 0.D0

    !   Heating radius of RF heating (m)
    RRF = 0.175D0

    !   RF input power (MW)
    PRFH = 0.D0

    !   ***** Neutral parameters *****

    !   Initial Neutral density (10^20 m^-3)
    PN0s = 1.D-8

    !   Neutral thermal velocity (m/s)
    V0 = 1.5D3

    !   Recycling rate in SOL
    rGamm0 = 0.8D0

    !   Gas-puff particle flux (10^20 m^-2 1/s)
    rGASPF = 0.1D0

    !   Electron density in diverter region (Minimum density in SOL)
    PNeDIV = 1.D-2

    !   Electron temperature in diverter region (Minimum Te in SOL)
    PTeDIV = 1.D-2

    !   Ion temperature in diverter region (Minimum Ti in SOL)
    PTiDIV = 1.D-2

    !   ***** Numerical parameters *****

    !   Implicitness parameter (Not used now)
    DLT = 1.0D0

    !   Time step size(s)
!!!!    DT = 1.D-4
    DT = 1.D-3

    !   Convergence parameter
    EPS = 1.D-2

    !   Iteration
    ICMAX = 10

    !   Time-advancing method
    !     ADV = 0     : Explicit scheme
    !           0.5   : Crank-Nicolson scheme
    !           2/3   : Galerkin scheme
    !           0.878 : Liniger scheme
    !           1     : Implicit scheme
    ADV = 1.D0

    !   ***** Mesh number parameters *****

    !   Magnitude of mesh peakness
    CMESH  = 30.D0

    !   Width of mesh peakness
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
    MODEGL=1

    !   Mode of AV
    !   0 : OFF
    !   n : Number of Display
    MODEAV = 0

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
    !   0    : original model
    !   1    : Hirshman, Hawryluk and Birge model
    MDLETA = 0

    !   Mode of fixed temperature profile
    !   0    : not fixed
    !   1    : fixed
    MDFIXT = 0

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
    gDIV(18) = 1.E6
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
    gDIV(36) = 1.E10
    gDIV(37) = 1.E20
    gDIV(38) = 1.E20
    gDIV(39) = 1.E23
    gDIV(40) = 1.E23
    gDIV(41) = 1.E20
    gDIV(42) = 1.E20
    gDIV(43) = 1.E20
    gDIV(44) = 1.E20
    gDIV(45) = 1.E20
    gDIV(46) = 1.E20
    gDIV(47) = 1.E23
    gDIV(60) = 1.E6
    gDIV(64) = 1.E-20
    gDIV(65) = 1.E-20
    gDIV(72) = 1.E3
    gDIV(73) = 1.E6
    gDIV(74) = 1.E6
    gDIV(75) = 1.E6
    gDIV(88) = 1.E3
    gDIV(89) = 1.E6

    !   Radius where density increase by command DEL
    DelR = 0.175D0

    !   Amount of increase of density by command DEL
    DelN = 5.D-1

    !   Helical ripple amplitude at r=a, Bhelical/Btoroidal, linear to r/a
    EpsH = 0.1D0

    !   Helical pitch number
    NCphi = 10

    !   Safety factor for helical
    Q0 = 3.D0
    QA = 2.D0

    NGR=-1

    RETURN
  END SUBROUTINE TXINIT

!***************************************************************
!
!        Calculate mesh, etc.
!
!***************************************************************

  SUBROUTINE TXCALM

    use physical_constants, only : AMP
    use libraries, only : LORENTZ, LORENTZ_PART, BISECTION

    INTEGER :: NR, NRL, NR_RC_NEAR
    real(8) :: DR, MAXAMP, CL, WL, RL, RCL, CLNEW

    !   Ion mass number
    AMI   = PA * AMP
    !   Beam ion mass number
    AMB   = AMI
    !   Radial step width
    DR    = RB / NRMAX
    !   Number of equations
    NQMAX = NQM
    !   Helical system
    UHth  = 1.D0 / SQRT(1.D0 + DBLE(NCphi)*2)
    UHph  = DBLE(NCphi) * UHth

    !  Mesh

!!$    CL  = 4.5D0
!!$    WL  = 5.D-2
!!$    MAXAMP = LORENTZ(RB,CL,WL,RC)
!!$    R(0) = 0.D0
!!$    DO NR = 1, NRMAX - 1
!!$       RL = DBLE(NR) / DBLE(NRMAX)
!!$       CALL BISECTION(LORENTZ,CL,WL,RC,MAXAMP,RL,RB,R(NR))
!!$    END DO
!!$    R(NRMAX) = RB

    !  As a trial, generate mesh using given CL and WL and seek the position in the
    !  original coordinate, which becomes the nearest mesh of separatrix after mapping.
    CL  = CMESH
    WL  = WMESH
    MAXAMP = LORENTZ(RB,CL,WL,RC) / RB
    NR_RC_NEAR = 0
    R(0) = 0.D0
    DO NR = 1, NRMAX - 1
       RL = DBLE(NR) / DBLE(NRMAX) * RB
       CALL BISECTION(LORENTZ,CL,WL,RC,MAXAMP,RL,RB,R(NR))
       IF(ABS(R(NR)-RC) <= ABS(R(NR)-R(NR_RC_NEAR))) NR_RC_NEAR = NR
    END DO
    R(NRMAX) = RB

    !  Construct new CL value that separatrix is just on mesh.
    !  New CL is chosen in order not to be settle so far from given CL.
    !  The mesh finally obtained is well-defined.
    RCL = DBLE(NR_RC_NEAR) / DBLE(NRMAX) * RB
    CLNEW = (RC - RCL) / (RCL * LORENTZ_PART(RB,WL,RC) / RB - LORENTZ_PART(RC,WL,RC))
    MAXAMP = LORENTZ(RB,CLNEW,WL,RC) / RB
    R(0) = 0.D0
    DO NR = 1, NRMAX - 1
       RL = DBLE(NR) / DBLE(NRMAX) * RB
       CALL BISECTION(LORENTZ,CLNEW,WL,RC,MAXAMP,RL,RB,R(NR))
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

    use physical_constants, only : AEE, AME, PI, rMU0, EPS0, rKeV
    use results
    use variables
    use libraries, only : INTG_P, DERIVS, INTDERIV3

    INTEGER :: NR, NQ, I, IER
    REAL(8) :: RL, PROF, PROFT, QL, RIP1, RIP2, dRIP, SSN, SSPe, SSPi
    REAL(8) :: ALP, dPe, dPi, DR1, DR2
    REAL(8) :: EpsL, Vte, Wte, rLnLam, rNuAsE_inv, FTL, EFT, ETAS, CR
    REAL(8) :: DERIV3, FCTR ! External functions
    real(8), dimension(:), allocatable :: AJPHL, TMP, RHSV
    real(8), dimension(:,:), allocatable :: CMTX

    NEMAX = NRMAX

    !  Define basic quantities like mass of particles, mesh, etc.

    CALL TXCALM

    !  Initialize variable vector

    X(1:NQMAX,0:NRMAX) = 0.D0

    !  Variables

    allocate(AJPHL(0:NRMAX))
    DO NR = 0, NRMAX
       RL=R(NR)
       IF (RL < RA) THEN
!          PROF  =(1.D0 - (RL / RA)**2)**2
          PROF  = 1.D0 - (RL / RA)**2
          PROFT = PROF**2
          ! Ne
          X(LQe1,NR) = (PN0 - PNa) * PROF + PNa
          ! Ni
          X(LQi1,NR) = X(LQe1,NR) / PZ
          ! Ne*Te
          X(LQe5,NR) = ((PTe0 - PTea) * PROFT + PTea) * X(LQe1,NR)
          ! Ni*Ti
          X(LQi5,NR) = ((PTi0 - PTia) * PROFT + PTia) * X(LQi1,NR)
       ELSE
          SSN = 2.D0 * (RB - RA) * (PN0 - PNa) / (RA * (PNa - PNeDIV))
          X(LQe1,NR) = (PNa - PNeDIV) * ((RB - RL) / (RB - RA))**SSN + PNeDIV!PNa * EXP(-(RL-RA) / rLn)!
          X(LQi1,NR) = X(LQe1,NR) / PZ
          SSPe = 2.D0 * (RB - RA) * (PN0 - PNa)       * PTea &
               & / (RA * (PNa*PTea - PNeDIV*PTeDIV))
          SSPi = 2.D0 * (RB - RA) *((PN0 - PNa) / PZ) * PTia &
               & / (RA * (PNa*PTia - PNeDIV*PTiDIV) / PZ)
          X(LQe5,NR) = ((PNa*PTea - PNeDIV*PTeDIV)      * ((RB - RL) / (RB - RA))**SSPe) &
               &     + PNeDIV*PTeDIV!PTea*X(LQe1,NR)
          X(LQi5,NR) = ((PNa*PTia - PNeDIV*PTiDIV) / PZ * ((RB - RL) / (RB - RA))**SSPi) &
               &     + PNeDIV*PTiDIV / PZ!PTia*X(LQi1,NR)
       END IF
       ! N0_1 (slow neutrals)
       X(LQn1,NR) = PN0s
       ! N0_2 (fast neutrals)
       X(LQn2,NR) = 0.D0
       ! Bphi
       X(LQm5,NR) = 0.5D0 * PSI(NR) * BB
       BphV(NR)   = BB
       ! Fixed densities to keep them constant during iterations
       PNeV_FIX(NR) = X(LQe1,NR)
       PTeV_FIX(NR) = X(LQe5,NR) / X(LQe1,NR)
    END DO
    CALL DERIVS(PSI,X(LQe1,0:NRMAX),dPNeV_FIX,NRMAX)

    ! Poloidal magnetic field

    IF(FSHL == 0.D0) THEN
       BthV(0) = 0.D0
       DO NR = 1, NRMAX
          RL = R(NR)
          IF(RL < RA) THEN
             PROF = 1.D0 - (RL / RA)**2
             ! Btheta
             BthV(NR) = rMU0 * rIPs * 1.D6 / (2.D0 * PI * RL) * (1.D0 - PROF**(PROFJ+1))
          ELSE
             BthV(NR) = rMU0 * rIPs * 1.D6 / (2.D0 * PI * RL)
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

    DO NR = 0, NRMAX
!!$       AJPHL(NR) = 2.D0 / rMU0 * DERIV3(NR,PSI,R(0:NRMAX)*BthV(0:NRMAX),NRMAX,NRM,0)
!!$       X(LQe4,NR) = - AJPHL(NR) / (AEE * 1.D20)

       IF((1.D0-(R(NR)/RA)**2) <= 0.D0) THEN
          PROF= 0.D0    
       ELSE             
          PROF= 1.D0-(R(NR)/RA)**2
       END IF
       IF(FSHL == 0.D0) THEN
          ! Ne*UePhi
          AJPHL(NR) = rIPs * 1.D6 / (PI * RA**2) * (PROFJ + 1.D0) * PROF**PROFJ
          X(LQe4,NR) = - AJPHL(NR) / (AEE * 1.D20)
          AJOH(NR)= PROF
       ELSE
          AJPHL(NR) = 0.D0
          X(LQe4,NR) = 0.D0
          AJOH(NR)= 0.D0
       END IF
    END DO

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

!!$    ! Numerical solution for AphV
!!$
!!$    RHSV(1:NRMAX) = - 0.5D0 * BthV(1:NRMAX) / R(1:NRMAX)
!!$    X(LQm4,0) = 0.D0
!!$    X(LQm4,1:NRMAX) = matmul(CMTX,RHSV)

    ! Analytic solution for AphV

    DO NR = 0, NRMAX
       IF(R(NR) < RA) THEN
          X(LQm4,NR) = - rMU0 * rIPs * 1.D6 / (4.D0 * PI * RA**2) &
               & *(3.d0*PSI(NR)-1.5D0*PSI(NR)**2/RA**2+PSI(NR)**3/(3.D0*RA**4))
       ELSE
          X(LQm4,NR) = - rMU0 * rIPs * 1.D6 / (4.D0 * PI * RA**2) &
               & *(3.d0*RA**2-1.5D0*RA**2+RA**2/3.D0) &
               & - rMU0 * rIPs * 1.D6 / (4.D0 * PI) * LOG(PSI(NR)/RA**2)
       END IF
    END DO

    ! Poloidal current density (Virtual current for helical system)

    allocate(TMP(0:NRMAX))
    IF(FSHL == 0.D0) THEN
       AJV(0:NRMAX)=0.D0
    ELSE
       TMP(0:NRMAX) = R(0:NRMAX) * BthV(0:NRMAX)
       DO NR = 0, NRMAX
          dRIP = DERIV3(NR,R,TMP,NRMAX,NRM,0) * 2.D0 * PI / rMU0
          AJV(NR)=dRIP / (2.D0 * PI * R(NR))
       END DO
    END IF
    deallocate(TMP)

    ! Toroidal electric field for initial NCLASS calculation

    Q(1:NRMAX) = ABS(R(1:NRMAX) * BB / (RR * BthV(1:NRMAX)))
    Q(0) = (4.D0 * Q(1) - Q(2)) / 3.D0

    DO NR = 0, NRMAX
       ! +++ Hirshman, Hawryluk and Birge model +++
       PNeV(NR) = X(LQe1,NR)
       PNiV(NR) = X(LQi1,NR)
       PTeV(NR) = X(LQe5,NR) / X(LQe1,NR)
       Vte = SQRT(2.D0 * ABS(PTeV(NR)) * rKeV / AME)
       Wte = Vte / (Q(NR) * RR)
       rLnLam = 15.2D0 - LOG(ABS(PNeV(NR))) / 2.D0 + LOG(ABS(PTeV(NR)))
       rNuei(NR) = PNiV(NR) * 1.D20 * PZ**2 * AEE**4 * rLnLam &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AME) &
            &     * (ABS(PTeV(NR)) * rKeV)**1.5D0)
       rNuAsE_inv = EpsL**1.5D0 * Wte / (SQRT(2.D0) * rNuei(NR))
       ! Inverse aspect ratio
       EpsL    = R(NR) / RR
       ! Trapped particle fraction
       FTL  = 1.46D0 * SQRT(EpsL) - 0.46D0 * EpsL**1.5D0
       EFT  = FTL * rNuAsE_inv / (rNuAsE_inv + (0.58D0 + 0.20D0 * Zeff))
       ! Spitzer resistivity
       ETAS = CORR(Zeff) * AME * rNuei(NR) / (PNeV(NR) * 1.D20 * AEE**2)
       CR   = 0.56D0 * (3.D0 - Zeff) / ((3.D0 + Zeff) * Zeff)
       ETA(NR) = ETAS * Zeff * (1.D0 + 0.27D0 * (Zeff - 1.D0)) &
            &  /((1.D0 - EFT) * (1.D0 - CR * EFT) * (1.D0 + 0.47D0 * (Zeff - 1.D0)))
       IF(FSHL == 0.D0) THEN
          ! Ephi
          X(LQm3,NR) = - ETA(NR) *  AJPHL(NR)
       ELSE
          X(LQm3,NR) = 0.D0
       END IF
       IF(X(LQm3,NR) == 0.D0) X(LQm3,NR) = -1.D-4
    END DO
    deallocate(AJPHL)

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

    DO NR = 0, NRMAX
       dPe = 2.D0 * R(NR) * DERIV3(NR,PSI,X(LQe5,0:NRMAX),NRMAX,NRM,0) * rKeV
       dPi = 2.D0 * R(NR) * DERIV3(NR,PSI,X(LQi5,0:NRMAX),NRMAX,NRM,0) * rKeV
       IF(rNueNC(NR) == 0.D0) THEN
          ALP = 0.D0
       ELSE
          ALP = (AMI / AME) * (rNuiNC(NR) / rNueNC(NR))
       END IF
       X(LQi3,NR) = (- BthV(NR) / BphV(NR) * X(LQe4,NR) + (dPe + dPi) / (AEE * BphV(NR)))&
            &     / (PZ + ALP) * R(NR)
       X(LQe3,NR) =- ALP * X(LQi3,NR)
       ErV(NR)    =- BphV(NR) / PNiV(NR) * (- BthV(NR) / BphV(NR) * X(LQe4,NR) &
            &      + (dPe + dPi) / (AEE * BphV(NR))) / (PZ + ALP) &
            &      + dPi / (PZ * AEE * PNiV(NR))
       X(LQe2,NR) =- (AMI * rNuiNC(NR) /(AEE * BphV(NR))) * X(LQi3,NR)
       X(LQi2,NR) = X(LQe2,NR) / PZ
       X(LQm3,NR) = BthV(NR) / PNeV(NR) * (-(AMI * rNuiNC(NR) /(AEE * BphV(NR))) &
            &     / (PZ + ALP) * (- BthV(NR) / BphV(NR) * X(LQe4,NR) &
            &     +(dPe + dPi) / (AEE * BphV(NR)))) &
            &     + AME * rNuei(NR) / (AEE * PNeV(NR)) * X(LQe4,NR)
    END DO

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
    X(LQm5,1:NRMAX) = matmul(CMTX,RHSV)
    deallocate(CMTX,RHSV,TMP)

    CALL TXCALV(X)
    CALL TXCALC

    !  Calculate global quantities for storing and showing initial status

    CALL TXGLOB

    RETURN
  END SUBROUTINE TXPROF

end module init_prof

module parameter_control
  use commons
  implicit none
  public
  NAMELIST /TX/ &
       & RA,RB,RC,RR,BB, &
       & PA,PZ,Zeff, &
       & PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ, &
       & De0,Di0,rMue0,rMui0,WPM0,WPE0,WPI0, &
       & Chie0,Chii0, &
       & FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD, &
       & FSCX,FSLC,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD0, &
       & rLn,rLT, &
       & Eb,RNB,PNBH,PNBCD,rNRF,RRF,PRFH, &
       & PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV, &
       & DLT,DT,EPS,ICMAX,ADV,CMESH,WMESH &
       & NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP, &
       & DelR,DelN, &
       & rG1,EpsH,FSHL,NCphi,Q0,QA, &
       & rIPs,rIPe, &
       & MODEG, gDIV, MODEAV, MODEGL, MDLPCK, MDLWTB, &
       & MDLETA, MDFIXT
  private :: TXPLST

contains
!***************************************************************
!
!   Change input parameters
!
!***************************************************************

  SUBROUTINE TXPARM(KID)

    INTEGER :: IST
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

    integer :: IST
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

    integer :: IST, KL
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

  !***** INPUT PARAMETER LIST *****

  SUBROUTINE TXPLST

    WRITE(6,601)
    RETURN

601 FORMAT(' ','# &TX : RA,RB,RC,RR,BB,PA,PZ,Zeff,'/ &
         &       ' ',8X,'PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ,'/ &
         &       ' ',8X,'De0,Di0,rMue0,rMui0,WPM0,WPE0,WPI0,'/ &
         &       ' ',8X,'Chie0,Chii0,'/ &
         &       ' ',8X,'FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD,'/ &
         &       ' ',8X,'FSCX,FSLC,FSNC,FSLP,FSLTE,FSLTI,FSION,FSD0,'/ &
         &       ' ',8X,'rLn,rLT,'/ &
         &       ' ',8X,'Eb,RNB,PNBH,PNBCD,rNRF,RRF,PRFH,'/ &
         &       ' ',8X,'PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV,'/ &
         &       ' ',8X,'DLT,DT,EPS,ICMAX,ADV,CMESH,WMESH,'/ &
         &       ' ',8X,'NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,'/ &
         &       ' ',8X,'DelR,DelN,'/ &
         &       ' ',8X,'rG1,EpsH,FSHL,NCphi,Q0,QA,'/ &
         &       ' ',8X,'rIPs,rIPe,'/ &
         &       ' ',8X,'MODEG, gDIV, MODEAV, MODEGL, MDLPCK,'/ &
         &       ' ',8X,'MDLWTB')
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
         &   'PROFJ ', PROFJ ,  'CMESH ', CMESH ,  &
         &   'WMESH ', WMESH ,  'ADV   ', ADV   ,  &
         &   'De0   ', De0   ,  'Di0   ', Di0   ,  &
         &   'rMue0 ', rMue0 ,  'rMui0 ', rMui0 ,  &
         &   'WPM0  ', WPM0  ,  'WPE0  ', WPE0  ,  &
         &   'WPI0  ', WPI0  ,  'PROFD ', PROFD ,  &
         &   'Chie0 ', Chie0 ,  'Chii0 ', Chii0 ,  &
         &   'FSDFIX', FSDFIX,  'FSCDBM', FSCDBM,  &
         &   'FSBOHM', FSBOHM,  'FSPSCL', FSPSCL,  &
         &   'FSCX  ', FSCX  ,  'FSLC  ', FSLC  ,  &
         &   'FSNC  ', FSNC  ,  'FSLP  ', FSLP  ,  &
         &   'FSLTE ', FSLTE ,  'FSLTI ', FSLTI ,  &
         &   'FSION ', FSION ,  'FSD0  ', FSD0  ,  &
         &   'rLn   ', rLn   ,  'rLT   ', rLT   ,  &
         &   'Eb    ', Eb    ,  'RNB   ', RNB   ,  &
         &   'PNBH  ', PNBH  ,  'rNRF  ', rNRF  ,  &
         &   'RRF   ', RRF   ,  'PRFH  ', PRFH  ,  &
         &   'rGamm0', rGamm0,  'V0    ', V0    ,  &
         &   'rGASPF', rGASPF,  'PNeDIV', PNeDIV,  &
         &   'PTeDIV', PTeDIV,  'PTiDIV', PTiDIV,  &
         &   'PN0s  ', PN0s  ,  'DLT   ', DLT   ,  &
         &   'EPS   ', EPS   ,  'DT    ', DT    ,  &
         &   'rG1   ', rG1   ,  'Zeff  ', Zeff  ,  &
         &   'rIPs  ', rIPs  ,  'rIPe  ', rIPe  ,  &
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
         &   'MDFIXT', MDFIXT,  'NCphi ', NCphi

    RETURN
  END SUBROUTINE TXVIEW
end module parameter_control
