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

    !   Partition radius (m) available if NRCMAX /= 0
    RC = 0.D0

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

    !   Particle loss to divertor parameter
    FSLP = 1.D0

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

    !   Gas-puff particle flux (10^20 1/s)
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

    !   ***** Mesh number parameters *****

    !   Number of nodes
    !     NRCMAX: first section (mainly core or whole region)
    !             No partitioning is carried out if NRCMAX=0 regardless of any RC values.
    NRCMAX = 0
    NRMAX  = 50

    !   Whether partitioning or not
    IDVD  = 0

    !   Number of elements
    NEMAX = NRMAX

    !   Number of time step
    NTMAX = 10

    !   Time step interval between print output
    NTSTEP = 10

    !   Time step interval between lines in f(r) graph
    NGRSTP = 1

    !   Time step interval between points in f(t) graph
    NGTSTP = 1

    !   Time step interval between points in f(t) graph
    NGVSTP = 1

    !   Mode of Graph
    !   1 : for Display
    !   2 : for Print Out
    MODEG = 2

    !   MODE of Graph Line
    !   0 : Change Line Color (Last Color Fixed)
    !   1 : Change Line Color and Style
    !   2 : Change Line Color, Style and Mark
    !   3 : Change Line Color, Style and Mark (With Legend)
    MODEL=1

    !   Mode of AV
    !   0 : OFF
    !   n : Number of Display
    MODEAV = 0

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
    gDIV(19) = 1.E6
    gDIV(21) = 1.E6
    gDIV(22) = 1.E6
    gDIV(23) = 1.E3
    gDIV(24) = 1.E3
    gDIV(25) = 1.E3
    gDIV(26) = 1.E3
    gDIV(27) = 1.E3
    gDIV(28) = 1.E3
    gDIV(35) = 1.E14
    gDIV(36) = 1.E14
    gDIV(37) = 1.E20
    gDIV(38) = 1.E20
    gDIV(39) = 1.E23
    gDIV(40) = 1.E23
    gDIV(41) = 1.E20
    gDIV(42) = 1.E20
    gDIV(43) = 1.E20
    gDIV(44) = 1.E20
    gDIV(45) = 1.E20
    gDIV(46) = 1.E23
    gDIV(47) = 1.E23
    gDIV(79) = 1.E6

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

    USE physical_constants, only : AMP

    INTEGER :: NR
    real(8) :: DR, DR1, DR2

    !   Ion mass number
    AMI   = PA * AMP
    !   Beam ion mass number
    AMB   = AMI
    !   Radial step width
    DR    = RB / NRMAX
    IF(NRCMAX /= 0.AND.RC /= 0) THEN
       IF(NRCMAX > NRMAX) STOP 'XX TXCALM: NRCMAX must be smaller than NRMAX!'
       DR1   = RC / NRCMAX
       DR2   =(RB - RC) / (NRMAX - NRCMAX)
       IDVD  = 1
    END IF
    !   Number of equations
    NQMAX = NQM
    !   Helical system
    UHth  = 1.D0 / SQRT(1.D0 + DBLE(NCphi)*2)
    UHph  = DBLE(NCphi) * UHth

    !  Mesh

    IF(IDVD == 0) THEN
       DO NR = 0, NRMAX
          R(NR) = DBLE(NR) * DR
       END DO
    ELSE
       DO NR = 0, NRCMAX
          R(NR) = DBLE(NR) * DR1
       END DO
       DO NR = 1, NRMAX - NRCMAX
          R(NR+NRCMAX) = RC + DBLE(NR) * DR2
       END DO
    END IF

    !  Maximum NR till RA

    DO NR = 0, NRMAX-1
       IF(R(NR) <= RA.AND.R(NR+1) >= RA) THEN
          NRA = NR
          EXIT
       END IF
    END DO

    !  Mesh interval

    H(1:NEMAX) = R(1:NRMAX) - R(0:NRMAX-1)

    RETURN
  END SUBROUTINE TXCALM

!***************************************************************
!
!   Initialize profiles
!
!***************************************************************

  SUBROUTINE TXPROF

    use physical_constants, only : AEE, AME, PI, rMU0, EPS0, rKEV
    use results
    use variables
    use libraries, only : F33, VALINT_SUB

    INTEGER :: NR, NQ
    REAL(8) :: RL, PROF, PROFT, QL, RIP1, RIP2, rLnLam, ETA, rJP, dRIP
    REAL(8) :: EpsL, rLnLame, RNZ, SGMSPTZ, FT, RNUE, F33TEFF
    REAL(8) :: ALP, dPe, dPi, VALE, VALI, BINCA
    REAL(8) :: DERIV3, FCTR ! External functions
    real(8), dimension(:), allocatable :: AJPHL, TMP, BINC

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
          PROF  = 1.D0 - (RL / RA)**2
          PROFT = PROF**2
          ! Ne
          X(LQe1,NR)  = (PN0 - PNa) * PROF + PNa
          ! Ni
          X(LQi1,NR)  = X(LQe1,NR) / PZ
          ! Ne*Te
          X(LQe5,NR) = ((PTe0 - PTea) * PROFT + PTea) * X(LQe1,NR)
          ! Ni*Ti
          X(LQi5,NR) = ((PTi0 - PTia) * PROFT + PTia) * X(LQi1,NR)
       ELSE
          X(LQe1,NR)  = PNa * EXP(-(RL-RA) / rLn)!(RL - RB)**2 + PNa - (RA - RB)**2
          X(LQi1,NR)  = X(LQe1,NR) / PZ
!!!!        X(LQe5,NR) = PTea * EXP(-(RL-RA) / rLT)
!!!!        X(LQi5,NR) = PTia * EXP(-(RL-RA) / rLT)
          X(LQe5,NR) = PTea*X(LQe1,NR)
          X(LQi5,NR) = PTia*X(LQi1,NR)
       END IF
       ! N0_1 (slow neutrals)
       X(LQn1,NR) = PN0s
       ! N0_2 (fast neutrals)
       X(LQn2,NR) = 0.D0
       ! Bphi
       X(LQm5,NR) = BB

       IF((1.D0-(R(NR)/RA)**PROFJ) <= 0.D0) THEN
          PROF= 0.D0    
       ELSE             
          PROF= 1.D0-(R(NR)/RA)**PROFJ
       END IF
       IF(FSHL == 0.D0) THEN
          ! Ne*UePhi
          AJPHL(NR) = rIPs * 1.D6 / (PI * RA**2) * (PROFJ + 1.D0) * PROF**PROFJ
          X(LQe4,NR) = - AJPHL(NR) / (AEE * 1.D20)
          AJOH(NR)= PROF
       ELSE
          X(LQe4,NR) = 0.D0
          AJOH(NR)= 0.D0
       END IF
    END DO

    ! Poloidal magnetic field

    IF(FSHL == 0.D0) THEN
       BthV(0) = 0.D0
       DO NR = 1, NRMAX
          RL=R(NR)
          IF (RL < RA) THEN
             PROF = 1.D0 - (RL / RA)**2
             ! Btheta
             BthV(NR) = rMU0 * rIPs * 1.D6 / (2.D0 * PI * RL) &
                  &              * (1.D0 - PROF**(PROFJ+1))
          ELSE
             BthV(NR) = rMU0 * rIPs * 1.D6 / (2.D0 * PI * RL)
          END IF
       END DO
    ELSE
       BthV(0) = 0.D0
       DO NR = 1, NRMAX
          RL=R(NR)
          QL=(Q0-QA)*(1.D0-(RL/RA)**2)+QA
          BthV(NR) = BB*RL/(QL*RR)
       END DO
    END IF
    DO NR = 0, NRMAX
       RL=R(NR)
       IF (RL < RA) THEN
          X(LQm4,NR) = - rMu0 * AEE * X(LQe4,NR) * 1.D20
       ELSE
          X(LQm4,NR) =0.D0
       END IF
    END DO

    ! Poloidal current density (Virtual current for helical system)

    allocate(TMP(0:NRMAX))
    IF(FSHL == 0.D0) THEN
       AJV(0:NRMAX)=0.D0
    ELSE
       TMP(0:NRMAX) = R(0:NRMAX) * X(LQm4,0:NRMAX)
       DO NR = 0, NRMAX
          dRIP = DERIV3(NR,R,TMP,NRMAX,NRM,0) * 2.D0 * PI / rMU0
          AJV(NR)=dRIP/(2.D0*PI*R(NR))
       END DO
    END IF
    deallocate(TMP)

    ! Toroidal electric field

    DO NR = 0, NRMAX
       RL=R(NR)
       IF (RL < RA) THEN
          PROF = (1.D0 - (RL / RA)**2)**PROFJ
       ELSE
          PROF = 0.D0
       END IF
       ! Ne on Integer Mesh
       PNeV(NR)=X(LQe1,NR)
       ! Te on Integer Mesh
       PTeV(NR)=X(LQe5,NR)/X(LQe1,NR)
       ! Neoclassical resistivity by Sauter et al.
       EpsL    = R(NR) / RR        ! Inverse aspect ratio
       rLnLame = 31.3D0-LOG(SQRT(PNeV(NR)*1.D20)/ABS(PTeV(NR)*1.D3)) ! Coulomb logarithm
       RNZ     = 0.58D0+0.74D0/(0.76D0+Zeff)
       SGMSPTZ = 1.9012D4*(PTeV(NR)*1.D3)**1.5D0/(Zeff*RNZ*rLnLame) ! Spitzer resistivity
       FT      = 1.46D0*SQRT(EpsL)-0.46D0*(EpsL)**1.5D0 ! Trapped particle fraction
       IF(NR == 0) THEN
          F33TEFF = 0.D0
       ELSE
          RNUE    = 6.921D-18*ABS(Q(NR))*RR*PNeV(NR)*1.D20 &
               &   *Zeff*rLnLame/((PTeV(NR)*1.D3)**2*SQRT(EpsL)**3)
          F33TEFF = FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(RNUE) &
               &   +0.45D0*(1.D0-FT)*RNUE/Zeff**1.5D0)
       END IF
       ETA     = 1.D0/(SGMSPTZ*F33(F33TEFF,Zeff))
       IF(FSHL == 0.D0) THEN
          ! Ephi
          X(LQm3,NR) = ETA *  AJPHL(NR)
       ELSE
          X(LQm3,NR) = 0.D0
       END IF
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
       dPe = DERIV3(NR,R,X(LQe5,0:NRMAX),NRMAX,NRM,0) * rKEV
       dPi = DERIV3(NR,R,X(LQi5,0:NRMAX),NRMAX,NRM,0) * rKEV
       ALP = PZ * (AME / AMI) * (rNueNC(NR) / rNuiNC(NR))
       IF(NR /= 0) THEN
          X(LQe3,NR) =(- PNiV(NR) * dPe - PNeV(NR) * dPi &
               &       + AEE * PNiV(NR) * BthV(NR) * X(LQe4,NR) &
               &       - AEE * PNeV(NR) * BthV(NR) * X(LQi4,NR)) &
               &     / ( AEE * PNiV(NR) * BphV(NR) * (1.D0 + ALP) * R(NR))
          X(LQi3,NR) =(  PNiV(NR) * dPe + PNeV(NR) * dPi &
               &       - AEE * PNiV(NR) * BthV(NR) * X(LQe4,NR) &
               &       + AEE * PNeV(NR) * BthV(NR) * X(LQi4,NR)) &
               &     * ALP / (AEE * PNiV(NR) * BphV(NR) * (1.D0 + ALP) * R(NR))
       END IF
       X(LQm2,NR) =(  PNiV(NR) * dPe + PNeV(NR) * dPi &
            &       - AEE * PNiV(NR) * BthV(NR) * X(LQe4,NR) &
            &       + AEE * PNeV(NR) * BthV(NR) * X(LQi4,NR)) &
            &     * rNueNC(NR) * AME &
            &     / (AEE**2 * PNiV(NR) * PNeV(NR) * BphV(NR) * (1.D0 + ALP))
       X(LQm1,NR) =(- PNiV(NR) * dPe * ALP + PNeV(NR) * dPi &
            &       + AEE * PNeV(NR) * BthV(NR) * X(LQi4,NR) &
            &       + AEE * PNiV(NR) * BthV(NR) * X(LQe4,NR) * ALP) &
            &     / ( AEE * PNeV(NR) * PNiV(NR) * (1.D0 + ALP))
!       write(6,'(I3,4F18.10)') NR,X(LQe3,NR),X(LQi3,NR),X(LQm3,NR),X(LQm1,NR)
    END DO
    X(LQe3,0) = FCTR(R(1),R(2),X(LQe3,1),X(LQe3,2))
    X(LQi3,0) = FCTR(R(1),R(2),X(LQi3,1),X(LQi3,2))

    allocate(BINC(0:NRMAX))
    DO NR = 0, NRMAX
       CALL VALINT_SUB(X(LQe3,0:NRMAX),NR,VALE)
       CALL VALINT_SUB(X(LQi3,0:NRMAX),NR,VALI)
       BINC(NR) = rMU0 * AEE * ( VALE * 1.D20 - PZ * VALI * 1.D20)
    END DO
    BINCA = BINC(NRMAX)
    BINC(0:NRMAX) = BINC(0:NRMAX) - BINCA
    X(LQm5,0:NRMAX) = BB + BINC(0:NRMAX)
    deallocate(BINC)

    X(LQe2,0:NRMAX) =((        X(LQm5,0:NRMAX) * De(0:NRMAX) / PTeV(0:NRMAX) / rKEV) &
         &          * (X(LQe3,0:NRMAX) - WPM(0:NRMAX) * R(0:NRMAX) * PNeV(0:NRMAX)) &
         &          - (PZ**2 * X(LQm5,0:NRMAX) * Di(0:NRMAX) / PTiV(0:NRMAX) / rKEV) &
         &          * (X(LQi3,0:NRMAX) - WPM(0:NRMAX) * R(0:NRMAX) * PNiV(0:NRMAX))) * AEE
    X(LQi2,0:NRMAX) = X(LQe2,0:NRMAX) / PZ
    X(LQm3,0:NRMAX) = X(LQm3,0:NRMAX) - BthV(0:NRMAX) * (X(LQe2,0:NRMAX) / PNeV(0:NRMAX) &
         &                                              +X(LQi2,0:NRMAX) / PNiV(0:NRMAX))

    CALL TXCALV(X)

    !  Calculate coefficient matrices for differential equations

    !  CALL TXCALA

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
       & De0,Di0,rMue0,rMui0,WPM0, &
       & Chie0,Chii0, &
       & FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD, &
       & FSCX,FSLC,FSNC,FSLP,FSION,FSD0, &
       & rLn,rLT, &
       & Eb,RNB,PNBH,PNBCD,rNRF,RRF,PRFH, &
       & PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV, &
       & DLT,DT,EPS,ICMAX, &
       & NRMAX,NRCMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP, &
       & DelR,DelN, &
       & rG1,EpsH,FSHL,NCphi,Q0,QA, &
       & rIPs,rIPe, &
       & MODEG, gDIV, MODEAV, MODEl
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
         &       ' ',8X,'De0,Di0,rMue0,rMui0,WPM0,'/ &
         &       ' ',8X,'Chie0,Chii0,'/ &
         &       ' ',8X,'FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD,'/ &
         &       ' ',8X,'FSCX,FSLC,FSNC,FSLP,FSION,FSD0,'/ &
         &       ' ',8X,'rLn,rLT,'/ &
         &       ' ',8X,'Eb,RNB,PNBH,PNBCD,rNRF,RRF,PRFH,'/ &
         &       ' ',8X,'PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV,'/ &
         &       ' ',8X,'DLT,DT,EPS,ICMAX,'/ &
         &       ' ',8X,'NRMAX,NRCMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,'/ &
         &       ' ',8X,'DelR,DelN,'/ &
         &       ' ',8X,'rG1,EpsH,FSHL,NCphi,Q0,QA,'/ &
         &       ' ',8X,'rIPs,rIPe,'/ &
         &       ' ',8X,'MODEG, gDIV, MODEAV, MODEl,')
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
         &   'PROFJ ', PROFJ , &
         &   'De0   ', De0   ,  'Di0   ', Di0   ,  &
         &   'rMue0 ', rMue0 ,  'rMui0 ', rMui0 ,  &
         &   'WPM0  ', WPM0  ,  'PROFD ', PROFD ,  &
         &   'Chie0 ', Chie0 ,  'Chii0 ', Chii0 ,  &
         &   'FSDFIX', FSDFIX,  'FSCDBM', FSCDBM,  &
         &   'FSBOHM', FSBOHM,  'FSPSCL', FSPSCL,  &
         &   'FSCX  ', FSCX  ,  'FSLC  ', FSLC  ,  &
         &   'FSNC  ', FSNC  ,  'FSLP  ', FSLP  ,  &
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
         &   'MODEl ', MODEl ,  'NCphi ', NCphi

    RETURN
  END SUBROUTINE TXVIEW
end module parameter_control
