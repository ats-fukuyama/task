!     $Id$
module init_prof
contains

!***************************************************************
!
!   Set constants and initial parameters
!
!***************************************************************

  SUBROUTINE TXINIT

    INCLUDE 'txcomm.inc'

    !   ***** Configuration parameters *****

    !   Plasma minor radius (m)
    RA = 0.35D0

    !   Wall radius (m)
    RB = 0.4D0

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
    ICMAX=10

    !   ***** Mesh number parameters *****

    !   Radial step number
    NRMAX = 50

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
    INCLUDE 'txcomm.inc'

    INTEGER :: NR

    !   Ion mass number
    AMi   = PA * AMP
    !   Beam ion mass number
    AMb   = AMi
    !   Radial step width
    DR    = RB / NRMAX
    !   Number of equations
    NQMAX = NQM
!    NQMAX = 1
    !   ???
    UHth  = 1.D0 / SQRT(1.D0 + DBLE(NCphi)*2)
    UHph  = DBLE(NCphi) * UHth

    !  Integer mesh

    DO NR = 0, NRMAX
       R(NR) = DBLE(NR) * DR
    END DO

    !  Half integer mesh

    DO NR = 0, NRMAX
       RHI(NR) = (DBLE(NR) + 0.5D0) * DR
    END DO

    RETURN
  END SUBROUTINE TXCALM

!***************************************************************
!
!   Initialize profiles
!
!***************************************************************

  SUBROUTINE TXPROF

    USE physical_constants, only : AEE, AME, PI, rMU0, EPS0, rKEV
    use results
    use variables
    INCLUDE 'txcomm.inc'

    INTEGER :: NR, NQ
    REAL(8) :: RL, PROF, PROFT, QL, RIP1, RIP2, rLnLam, ETA, rJP

    !  Define basic quantities like mass of particles, mesh, etc.

    CALL TXCALM

    !  Define radial coordinate for graph

    CALL TXPRFG

    !  Initialize variable vector

    X(1:NQMAX,0:NRMAX) = 0.D0

    !  Half integer mesh variables

    DO NR = 0, NRMAX - 1
       RL=RHI(NR)
       IF (RL < RA) THEN
          PROF  = 1.D0 - (RL / RA)**2
          PROFT = PROF**2
          ! Ne
          X(LQe1,NR)  = (PN0 - PNa) * PROF + PNa
          ! Ni
          X(LQi1,NR)  = X(LQe1,NR) / PZ
          ! Ne*Te
          X(LQe5,NR) = ((PTe0 - PTea) * PROFT + PTea)*X(LQe1,NR)
          ! Ni*Ti
          X(LQi5,NR) = ((PTi0 - PTia) * PROFT + PTia)*X(LQi1,NR)
       ELSE
          X(LQe1,NR)  = PNa * EXP(-(RL-RA) / rLn)
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

       IF((1.D0-(RHI(NR)/RA)**PROFJ) <= 0.D0) THEN
          PROF=0.D0    
       ELSE             
          PROF= (1.D0-(RHI(NR)/RA)**PROFJ)
       END IF
       IF(FSHL == 0.D0) THEN
          ! Ne*UePhi
          X(LQe4,NR) = - rIPs * 1.D6 / (AEE * PI * RA**2 * 1.D20) &
               &                     * (PROFJ + 1) * PROF**PROFJ
          AJOH(NR)= PROF
       ELSE
          X(LQe4,NR) = 0.D0
          AJOH(NR)= 0.D0
       END IF
    END DO

    ! Integer mesh variables

    IF(FSHL == 0.D0) THEN
       X(LQm4,0) = 0.D0
       DO NR = 1, NRMAX
          RL=R(NR)
          IF (RL < RA) THEN
             PROF = 1.D0 - (RL / RA)**2
             ! Btheta
             X(LQm4,NR) = rMU0 * rIPs * 1.D6 / (2 * PI * RL) &
                  &              * (1 - PROF**(PROFJ+1))
          ELSE
             X(LQm4,NR) = rMU0 * rIPs * 1.D6 / (2 * PI * RL)
          END IF
       END DO
    ELSE
       X(LQm4,0) = 0.D0
       DO NR = 1, NRMAX
          RL=R(NR)
          QL=(Q0-QA)*(1-(RL/RA)**2)+QA
          X(LQm4,NR) = BB*RL/(QL*RR)
       END DO
    END IF

    IF(FSHL == 0.D0) THEN
       AJV(0:NRMAX)=0.D0
    ELSE
       DO NR=0,NRMAX
          RL=R(NR)
          RIP1=2.D0*PI*RL*X(LQm4,NR)/rMU0
          RL=R(NR+1)
          RIP2=2.D0*PI*RL*X(LQm4,NR+1)/rMU0
          RL=RHI(NR)
          ! Poloidal current density
          AJV(NR)=(RIP2-RIP1)/(2.D0*PI*RL*DR)
       END DO
    END IF

    DO NR = 0, NRMAX - 1
       RL=RHI(NR)
       IF (RL < RA) THEN
          PROF = (1.D0 - (RL / RA)**2)**PROFJ
       ELSE
          PROF = 0.D0
       END IF
       ! Ne on Half Integer Mesh
       PNeHI(NR)=X(LQe1,NR)
       ! Te on Half Integer Mesh
       PTeHI(NR)=X(LQe5,NR)
       ! Coulomb logarithm
       rLnLam = 15.D0 - LOG(ABS(PNeHI(NR))) / 2.D0 + LOG(ABS(PTeHI(NR)))
       ! Resistivity
       ETA =  SQRT(AME) * Zeff  * AEE**2 * rLnLam &
            &           / (3.D0 * (2.D0 * PI)**1.5D0 * EPS0**2 &
            &                   * (ABS(PTeHI(NR)) * rKeV)**1.5D0)
       ! Jphi
       rJP = rIPs * 1.D6 / (PI * RA**2) * (PROFJ + 1.D0) * PROF
       IF(FSHL == 0.D0) THEN
          ! Ephi
          X(LQm3,NR) = ETA * rJP
       ELSE
          X(LQm3,NR) = 0.D0
       END IF
    END DO

    T_TX=0.D0
    NGT=-1
    NGR=-1
    NGVV=-1
    rIP=rIPs

    !  Define physical variables from X on half and integer mesh

    CALL TXCALV(X)

    !  Calculate various physical quantities

    CALL TXCALC

    !  Calculate coefficient matrices for differential equations

    !  CALL TXCALA

    !  Calculate global quantities for storing and showing initial status

    CALL TXGLOB

    !  Store center or edge values of variables for showing time-evolution graph

    CALL TXSTGT(SNGL(T_TX))

    !  Store global quantities for showing time-evolution graph

    CALL TXSTGV(SNGL(T_TX))

    !  Store profile data for showing graph

    CALL TXSTGR

    RETURN
  END SUBROUTINE TXPROF

end module init_prof

module parameters
contains
!***************************************************************
!
!   Change input parameters
!
!***************************************************************

  SUBROUTINE TXPARM(KID)

    INCLUDE 'txcomm.inc'

    NAMELIST /TX/ &
         & RA,RB,RR,BB, &
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
         & NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP, &
         & DelR,DelN, &
         & rG1,EpsH,FSHL,NCphi,Q0,QA, &
         & rIPs,rIPe, &
         & MODEG, gDIV, MODEAV, MODEl

    INTEGER :: MODE, IST, KL
    LOGICAL :: LEX
    CHARACTER(80) :: KPNAME, KLINE, KNAME*90, KID*1

    MODE=0
1   CONTINUE
    WRITE(6,*) '# INPUT &TX :'
    READ(5,TX,IOSTAT=IST)
    IF(IST > 0) THEN
       CALL TXPLST
       GO TO 1
    ELSE IF(IST < 0) THEN
       KID='Q'
    ELSE
       KID=' '
    END IF
    GOTO 3000

    ENTRY TXPARL(KLINE)

    MODE=1
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
    GOTO 3000

    ENTRY TXPARF(KPNAME)

    MODE=2
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

3000 IERR=0

    !    ERROR CHECK

    IF(IERR /= 0.AND.MODE == 0) GO TO 1

    RETURN

  END SUBROUTINE TXPARM

  !***** INPUT PARAMETER LIST *****

  SUBROUTINE TXPLST

    WRITE(6,601)
    RETURN

601 FORMAT(' ','# &TX : RA,RB,RR,BB,PA,PZ,Zeff,'/ &
         &       ' ',8X,'PN0,PNa,PTe0,PTea,PTi0,PTia,PROFJ,'/ &
         &       ' ',8X,'De0,Di0,rMue0,rMui0,WPM0,'/ &
         &       ' ',8X,'Chie0,Chii0,'/ &
         &       ' ',8X,'FSDFIX,FSCDBM,FSBOHM,FSPSCL,PROFD,'/ &
         &       ' ',8X,'FSCX,FSLC,FSNC,FSLP,FSION,FSD0,'/ &
         &       ' ',8X,'rLn,rLT,'/ &
         &       ' ',8X,'Eb,RNB,PNBH,PNBCD,rNRF,RRF,PRFH,'/ &
         &       ' ',8X,'PN0s,V0,rGamm0,rGASPF,PNeDIV,PTeDIV,PTiDIV,'/ &
         &       ' ',8X,'DLT,DT,EPS,ICMAX,'/ &
         &       ' ',8X,'NRMAX,NTMAX,NTSTEP,NGRSTP,NGTSTP,NGVSTP,'/ &
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

    INCLUDE 'txcomm.inc'

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
end module parameters
