C     $Id$
C
C     ***********************************************************
C
C           INITIALIZE CONSTANTS AND DEFAULT VALUES
C
C     ***********************************************************
C
      SUBROUTINE TRINIT
C
      INCLUDE 'trcomm.inc'
C
C     *** CONSTANTS ****
C
C        PI    : Pi
C        AEE   : Elementaty charge
C        AME   : Electron mass
C        AMM   : Proton mass
C        VC    : Speed of light in vacuum
C        RMU0  : Permeability of free space
C        EPS0  : Permittivity of free space
C        VOID  : 0.D0
C
      PI      = ASIN(1.D0)*2.D0
      AEE     = 1.60217733D-19
      AME     = 9.1093897D-31
      AMM     = 1.6726231D-27
      VC      = 2.99792458D8
      RMU0    = 4.D0*PI*1.D-7
      EPS0    = 1.D0/(VC*VC*RMU0)
      RKEV    = AEE*1.D3
      VOID    = 0.D0
C
C     ==== DEVICE PARAMETERS ====
C
C        RR     : PLASMA MAJOR RADIUS (M)
C        RA     : PLASMA MINOR RADIUS (M)
C        RKAP   : ELIPTICITY OF POLOIDAL CROSS SECTION
C        RDLT   : TRIANGULARITY OF POLOIDAL CROSS SECTION
C        BB     : TOROIDAL MAGNETIC FIELD ON PLASMA AXIS (T)
C        RIPS   : INITIAL VALUE OF PLASMA CURRENT (MA)
C        RIPE   : FINAL VALUE OF PLASMA CURRENT (MA)
C        RHOA   : EDGE OF CALCULATE REGION (NORMALIZED SMALL RADIUS)
C
      RR      = 3.0D0
      RA      = 1.2D0
      RKAP    = 1.5D0
      RDLT    = 0.0D0
      BB      = 3.D0
      RIPS    = 3.D0
      RIPE    = 3.D0
      RHOA    = 1.D0
C
C     ==== PLASMA PARAMETERS ====
C
C        NSMAX  : NUMBER OF MAIN PARTICLE SPECIES (NS=1:ELECTRON)
C        NSZMAX : NUMBER OF IMPURITIES SPECIES
C        NSNMAX : NUMBER OF NEUTRAL SPECIES
C
C        PA(NS) : ATOMIC NUMBER
C        PZ(NS) : CHARGE NUMBER
C        PN(NS) : INITIAL NUMBER DENSITY ON AXIS (1.E20 M**-3)
C        PNS(NS): INITIAL NUMBER DENSITY ON SURFACE (1.E20 M**-3)
C        PT(NS) : INITIAL TEMPERATURE ON AXIS (KEV)
C        PTS(IS): INITIAL TEMPERATURE ON SURFACE (KEV)
C
      NSMAX=2
      NSZMAX=0  ! the number of impurities
      NSNMAX=0  ! the number of neutrals, 0 or 2 fixed
C
      PA(1)   = AME/AMM
      PZ(1)   =-1.D0
      PN(1)   = 0.5D0
      PT(1)   = 1.5D0
      PTS(1)  = 0.05D0
      PNS(1)  = 0.05D0
C
      PA(2)   = 2.D0
      PZ(2)   = 1.D0
      PN(2)   = 0.5D0-2.D-7
      PT(2)   = 1.5D0
      PTS(2)  = 0.05D0
      PNS(2)  = 0.05D0-2.D-8
C
      PA(3)   = 3.D0
      PZ(3)   = 1.D0
      PN(3)   = 1.D-7
      PT(3)   = 1.5D0
      PTS(3)  = 0.05D0
      PNS(3)  = 1.D-8
C
      PA(4)   = 4.D0
      PZ(4)   = 2.D0
      PN(4)   = 1.D-7
      PT(4)   = 1.5D0
      PTS(4)  = 0.05D0
      PNS(4)  = 1.D-8
C
      PA(5)   = 12.D0
      PZ(5)   = 2.D0
      PN(5)   = VOID
      PT(5)   = 0.D0
      PTS(5)  = 0.D0
      PNS(5)  = VOID
C
      PA(6)   = 12.D0
      PZ(6)   = 4.D0
      PN(6)   = VOID
      PT(6)   = 0.D0
      PTS(6)  = 0.D0
      PNS(6)  = VOID
C
      PA(7)   = 2.D0
      PZ(7)   = 0.D0
      PN(7)   = 1.D-15
      PT(7)   = 0.D0
      PTS(7)  = 0.D0
      PNS(7)  = 2.D-4
C
      PA(8)   = 2.D0
      PZ(8)   = 0.D0
      PN(8)   = 1.D-15
      PT(8)   = 0.D0
      PTS(8)  = 0.D0
      PNS(8)  = 1.D-15
C
C     ==== IMPURITY PARAMETERS ====
C
C        PNC    : CARBON DENSITY FACTOR
C        PNFE   : IRON DENSITY FACTOR
C                      COMPARED WITH ITER PHYSICS DESIGN GUIDELINE
C        PNNU   : NEUTRAL NUMBER DENSITY ON AXIS (1.E20 M**-3)
C        PNNUS  :                        ON SURFACE (1.E20 M**-3)
C
      PNC     = 0.D0
      PNFE    = 0.D0
      PNNU    = 0.D0
      PNNUS   = 0.D0
C
C     ==== PROFILE PARAMETERS ====
C
C        PROFN*: PROFILE PARAMETER OF INITIAL DENSITY
C        PROFT*: PROFILE PARAMETER OF INITIAL TEMPERATURE
C        PROFU*: PROFILE PARAMETER OF INITIAL NEUTRAL DENSITY
C        PROFJ*: PROFILE PARAMETER OF INITIAL CURRENT DENSITY
C                    (X0-XS)(1-RHO**PROFX1)**PROFX2+XS
C
C        ALP   : ADDITIONAL PARAMETERS
C           ALP(1): RADIUS REDUCTION FACTOR
C           ALP(2): MASS WEIGHTING FACTOR FOR NC
C           ALP(3): CHARGE WEIGHTING FACTOR FOR NC
C
      PROFN1 = 2.D0
      PROFN2 = 0.5D0
      PROFT1 = 2.D0
      PROFT2 = 1.D0
      PROFU1 =12.D0
      PROFU2 = 1.D0
      PROFJ1 =-2.D0
      PROFJ2 = 1.D0
C
      ALP(1) = 1.0D0
      ALP(2) = 0.D0
      ALP(3) = 0.D0
C
C     ==== TRANSPORT PARAMETERS ====
C
C        AV0    : INWARD PARTICLE PINCH FACTOR
C        AD0    : PARTICLE DIFFUSION FACTOR
C        CNP    : COEFFICIENT FOR NEOCLASICAL PARTICLE DIFFUSION
C        CNH    : COEFFICIENT FOR NEOCLASICAL HEAT DIFFUSION
C        CDP    : COEFFICIENT FOR TURBULENT PARTICLE DIFFUSION
C        CDH    : COEFFICIENT FOR TURBLUENT HEAT DIFFUSION
C        CDW(8) : COEFFICIENTS FOR DW MODEL
C
      AV0    = 0.5D0
      AD0    = 0.5D0
C
      CNP    = 1.D0
      CNH    = 1.D0
      CDP    = 1.D0
      CDH    = 1.D0
      CDW(1) = 0.04D0
      CDW(2) = 0.04D0
      CDW(3) = 0.04D0
      CDW(4) = 0.04D0
      CDW(5) = 0.04D0
      CDW(6) = 0.04D0
      CDW(7) = 0.04D0
      CDW(8) = 0.04D0
C
C     ==== TRANSPORT MODEL ====
C
C        MDLKAI: TURBULENT TRANSPORT MODEL
C
C   *************************************************************
C   ***  0.GE.MDLKAI.LT.10 : CONSTANT COEFFICIENT MODEL       ***
C   *** 10.GE.MDLKAI.LT.20 : DRIFT WAVE (+ITG +ETG) MODEL     ***
C   *** 20.GE.MDLKAI.LT.30 : REBU-LALLA MODEL                 ***
C   *** 30.GE.MDLKAI.LT.40 : CURRENT-DIFFUSIVITY DRIVEN MODEL ***
C   *** 40.GE.MDLKAI.LT.60 : DRIFT WAVE BALLOONING MODEL      ***
C   ***       MDLKAI.GE.60 : ITG(/TEM, ETG) MODEL ETC         ***
C   *************************************************************
C
C      ***  MDLKAI.EQ. 0   : CONSTANT                              ***
C      ***  MDLKAI.EQ. 1   : CONSTANT/(1-A*rho^2)                  ***
C      ***  MDLKAI.EQ. 2   : CONSTANT*(dTi/drho)^B/(1-A*rho^2)     ***
C      ***  MDLKAI.EQ. 3   : CONSTANT*(dTi/drho)^B*Ti^C            ***
C
C      ***  MDLKAI.EQ. 10  : etac=1                                ***
C      ***  MDLKAI.EQ. 11  : etac=1 1/(1+exp)                      ***
C      ***  MDLKAI.EQ. 12  : etac=1 1/(1+exp) *q                   ***
C      ***  MDLKAI.EQ. 13  : etac=1 1/(1+exp) *(1+q^2)             ***
C      ***  MDLKAI.EQ. 14  : etac=1+2.5*(Ln/RR-0.2) 1/(1+exp)      ***
C      ***  MDLKAI.EQ. 15  : etac=1 1/(1+exp) func(q,eps,Ln)       ***
C
C      ***  MDLKAI.EQ. 20  : Rebu-Lalla model                      ***
C
C      ***  MDLKAI.EQ. 30  : CDBM 1/(1+s)                          ***
C      ***  MDLKAI.EQ. 31  : CDBM F(s,alpha,kappaq)                ***
C      ***  MDLKAI.EQ. 32  : CDBM F(s,alpha,kappaq)/(1+WE1^2)      ***
C      ***  MDLKAI.EQ. 33  : CDBM F(s,0,kappaq)                    ***
C      ***  MDLKAI.EQ. 34  : CDBM F(s,0,kappaq)/(1+WE1^2)          ***
C      ***  MDLKAI.EQ. 35  : CDBM (s-alpha)^2/(1+s^2.5)            ***
C      ***  MDLKAI.EQ. 36  : CDBM (s-alpha)^2/(1+s^2.5)/(1+WE1^2)  ***
C      ***  MDLKAI.EQ. 37  : CDBM s^2/(1+s^2.5)                    ***
C      ***  MDLKAI.EQ. 38  : CDBM s^2/(1+s^2.5)/(1+WE1^2)          ***
C      ***  MDLKAI.EQ. 39  : CDBM F2(s,alpha,kappaq,a/R)           ***
C      ***  MDLKAI.EQ. 40  : CDBM F3(s,alpha,kappaq,a/R)/(1+WS1^2) ***
C
C      ***  MDLKAI.EQ. 60  : GLF23 model                           ***
C      ***  MDLKAI.EQ. 61  : GLF23 (stability enhanced version)    ***
C      ***  MDLKAI.EQ. 62  : IFS/PPPL model                        ***
C      ***  MDLKAI.EQ. 63  : Weiland model                         ***
C      ***  MDLKAI.EQ. 64  : Modified Weiland model                ***
C      ***  MDLKAI.EQ. 65  : Bohm/Gyro-Bohm model                  ***
C
C     +++++ WARNING +++++++++++++++++++++++++++++++++++++++++++
C     +  Parameters below are valid only if MDNCLS /= 0,      +
C     +  that is, one do not use NCLASS,                      +
C     +  otherwise NCLASS automatically calculates all        +
C     +  variables in the following:                          +
C     +                                                       +
C     +    MDLETA: RESISTIVITY MODEL                          +
C     +               1: Hinton and Hazeltine                 +
C     +               2: Hirshman, Hawryluk                   +
C     +               3: Sauter                               +
C     +               4: Hirshman, Sigmar                     +
C     +               else: CLASSICAL                         +
C     +    MDLAD : PARTICLE DIFFUSION MODEL                   +
C     +               1: CONSTANT D                           +
C     +               2: TURBULENT EFFECT                     +
C     +               3: Hinton and Hazeltine                 +
C     +               4: Hinton and Hazeltine w/ TURBULENT    +
C     +               else: NO PARTICLE TRANSPORT             +
C     +    MDLAVK: HEAT PINCH MODEL                           +
C     +               1: Arbitrary amplitude                  +
C     +               2: Arbitrary amplitude w/ pressure dep. +
C     +               3: Hinton and Hazeltine                 +
C     +               else: NO HEAT PINCH                     +
C     +    MDLJBS: BOOTSTRAP CURRENT MODEL                    +
C     +               1-3: Hinton and Hazeltine               +
C     +               4: Hirshman, Sigmar                     +
C     +               5: Sauter                               +
C     +               else: Hinton and Hazeltine              +
C     +    MDLKNS: NEOCLASSICAL TRANSPORT MODEL               +
C     +            0    : Hinton and Hazeltine                +
C     +            else : Chang and Hinton                    +
C     +                                                       +
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C        MDLTPF: TRAPPED PARTICLE FRACTION MODEL
C
      MDLKAI = 31
      MDLETA = 3
      MDLAD  = 3
      MDLAVK = 3
      MDLJBS = 5
      MDLKNC = 1
      MDLTPF = 0
C
C        MDLWLD : Weiland model mode selector
C            0    : using effective transport coefficients
C            else : using transport coefficients' matrices
C
      MDLWLD=0
C
C        MDDW : mode selector for anom. particle transport coefficient
C            you must NOT modify this parameter.
C            0    : if MDDW=0 from start to finish when you choose
C                   a certain transport model (MDLKAI),
C                   you could control a ratio of anomalous particle
C                   transport to total particle transport to manipulate
C                   the factor of AD0.
C            else : this is because you chose MDLKAI=60, 61, or 63
C                   which assign the transport models that can calculate
C                   an anomalous particle transport coefficient
C                   on their own.
      MDDW=0
C
C     ==== Semi-Empirical Parameter for Anomalous Transport ====
C
      CHP    = 0.D0
      CK0    = 12.D0 ! for electron
      CK1    = 12.D0 ! for ions
      CWEB   = 0.D0  ! for omega ExB
      CALF   = 1.D0  ! for s-alpha
      CKALFA = 0.D0
      CKBETA = 0.D0
      CKGUMA = 0.D0
C
C     ==== CONTROL PARAMETERS ====
C
C        DT     : SIZE OF TIME STEP
C        NRMAX  : NUMBER OF RADIAL MESH POINTS
C        NTMAX  : NUMBER OF TIME STEP
C        NTSTEP : INTERVAL OF SNAP DATA PRINT
C        NGRSTP : INTERVAL OF RADIAL PROFILE SAVE
C        NGTSTP : INTERVAL OF TIME EVOLUTION SAVE
C        NGPST  : ???
C        TSST   : ???
C
      DT     = 0.01D0 
      NRMAX  = 50
      NTMAX  = 100
      NTSTEP = 10
      NGRSTP = 100
      NGTSTP = 2
      NGPST  = 4
      TSST   = 1.D9
C
C     ==== Convergence Parameter ====
C
C        EPSLTR : CONVERGENCE CRITERION OF ITERATION
C        LMAXTR : MAXIMUM COUNT OF ITERATION
C
      EPSLTR = 0.001D0
C      EPSLTR = 1.D99
      LMAXTR = 10
C
C     ==== SAWTOOTH PARAMETERS ====
C
C        TPRST  : SAWTOOTH PERIOD (S)
C        MDLST  : SAWTOOTH MODEL TYPE
C                    0:OFF
C                    1:ON
C        IZERO  : SAWTOOTH CRASH TYPE
C
      TPRST  = 0.1D0
      MDLST  = 0
      IZERO  = 3
C
C     ==== FUSION REACTION PARAMETERS ====
C
C        MDLNF  : FUSION REACTION MODEL TYPE
C                    0:OFF
C                    1:ON
C
      MDLNF  = 0
C
C     ==== NBI HEATING PARAMETERS ====
C
C        PNBTOT : NBI TOTAL INPUT POWER (MW)
C        PNBR0  : RADIAL POSITION OF NBI POWER DEPOSITION (M)
C        PNBRW  : RADIAL WIDTH OF NBI POWER DEPOSITION (M)
C        PNBENG : NBI BEAM ENERGY (keV)
C        PNBRTG : TANGENTIAL RADIUS OF NBI BEAM (M)
C        PNBCD  : CURRENT DRIVE FACTOR
C        MDLNB  : NBI MODEL TYPE
C                    0:OFF
C                    1:GAUSSIAN
C                    2:PENCIL BEAM
C
      PNBTOT = 0.D0
      PNBR0  = 0.D0
      PNBRW  = 0.5D0
      PNBENG = 80.D0
      PNBRTG = 3.D0
      PNBCD  = 1.D0
      MDLNB  = 1
C
C     ==== ECRF PARAMETERS ====
C
C        PECTOT : ECRF INPUT POWER (MW)
C        PECR0  : RADIAL POSITION OF POWER DEPOSITION (M)
C        PECRW  : RADIAL WIDTH OF POWER DEPOSITION (M)
C        PECTOE : POWER PARTITION TO ELECTRON
C        PECNPR : PARALLEL REFRACTIVE INDEX
C        PECCD  : CURRENT DRIVE FACTOR
C        MDLEC  : ECRF MODEL
C
      PECTOT = 0.D0
      PECR0  = 0.D0
      PECRW  = 0.2D0
      PECTOE = 1.D0
      PECNPR = 0.D0
      PECCD  = 0.D0
      MDLEC  = 0
C
C     ==== LHRF PARAMETERS ====
C
C        PLHTOT : LHRF INPUT POWER (MW)
C        PLHR0  : RADIAL POSITION OF POWER DEPOSITION (M)
C        PLHRW  : RADIAL WIDTH OF POWER DEPOSITION (M)
C        PLHTOE : POWER PARTITION TO ELECTRON
C        PLHNPR : PARALLEL REFRACTIVE INDEX
C        PLHCD  : CURRENT DRIVE FACTOR
C        MDLLH  : LHRF MODEL
C
      PLHTOT = 0.D0
      PLHR0  = 0.D0
      PLHRW  = 0.2D0
      PLHTOE = 1.D0
      PLHNPR = 2.D0
      MDLLH  = 0
C
C     ==== ICRF PARAMETERS ====
C
C        PICTOT : ICRF INPUT POWER (MW)
C        PICR0  : RADIAL POSITION OF POWER DEPOSITION (M)
C        PICRW  : RADIAL WIDTH OF POWER DEPOSITION (M)
C        PICTOE : POWER PARTITION TO ELECTRON
C        PICNPR : PARALLEL REFRACTIVE INDEX
C        PICCD  : CURRENT DRIVE FACTOR
C        MDLIC  : ICRF MODEL
C
      PICTOT = 0.D0
      PICR0  = 0.D0
      PICRW  = 0.5D0
      PICTOE = 0.5D0
      PICNPR = 2.D0
      PICCD  = 0.D0
      MDLIC  = 0
C
C     ==== CURRENT DRIVE PARAMETERS ====
C
C        PBSCD : BOOTSTRAP CURRENT DRIVE FACTOR
C        MDLCD : CURRENT DRIVE OPERATION MODEL
C                  0: TOTAL PLASMA CURRENT FIXED
C                  1: TOTAL PLASMA CURRENT VARIABLE
C
      PBSCD  = 1.D0
      MDLCD  = 0
C
C     ==== PELLET INJECTION PARAMETERS ====
C
C        MDLPEL : PELLET INJECTION MODEL TYPE
C                    0:OFF  1:GAUSSIAN  2:NAKAMURA  3:HO
C        PELTOT : TOTAL NUMBER OF PARTICLES IN PELLET
C        PELR0  : RADIAL POSITION OF PELLET DEPOSITION (M)
C        PELRW  : RADIAL WIDTH OF PELLET DEPOSITION (M)
C        PELRAD : RADIUS OF PELLET (M)
C        PELVEL : PELLET INJECTION VELOCITY (M/S)
C        PELTIM : TIME FOR PELLET TO BE INJECTED
C        PELPAT : PARTICLE RATIO IN PELLET'
C
      MDLPEL = 1
      PELTOT = 0.D0
      PELR0  = 0.D0
      PELRW  = 0.5D0
      PELRAD = 0.D0
      PELVEL = 0.D0
      PELTIM = -10.D0
C
      DO NS=1,NSMAX
         PELPAT(NS) = 1.0D0
      ENDDO
C
C     ==== DEVICE NAME AND SHOT NUMBER IN UFILE DATA ====
C        KUFDEV : DEVICE NAME
C        KUFDCG : DISCHARGE NUMBER
C
      KUFDEV=''
      KUFDCG=''
C
C     ==== FILE NAME ====
C
C        KNAMEQ: Filename of equilibrium data
C        KNAMTR: Filename of transport data
C        KFNLOG : LOG FILE NAME
C
      KNAMEQ='eqdata'
      KNAMTR='trdata'
      KFNLOG='trf.log'
C
C     ==== INTERACTION WITH EQ ====
C
C        MODELG: 2 : TOROIDAL GEOMETRY
C                3 : READ TASK/EQ FILE
C                5 : READ EQDSK FILE
C                9 : CALCULATE TASK/EQ
C        NTEQIT: STEP INTERVAL OF EQ CALCULATION
C                0 : INITIAL EQUILIBRIUM ONLY
C
      MODELG=2
      NTEQIT=0
C
C     ==== INPUT FROM EXPERIMENTAL DATA ====
C
C        MDLXP :
C           0 : from ufiles
C        else : MDSplus
C
C        MDLUF :
C           0 : not used
C           1 : time evolution
C           2 : steady state
C           3 : compared with TOPICS
C
      MDLXP=0
      MDLUF=0
C
C     ==== IMPURITY TREATMENT ====
C
C        MDNI  :
C           0 : NSMAX=2, ne=ni
C           1 : calculate nimp and zeff profiles from NE, ZIMP and NM1
C           2 : calculate nimp and ni profiles from NE, ZIMP and ZEFFR
C           3 : calculate zeff and ni profiles from NE, ZIMP and NIMP
C
      MDNI=0
C
C     ==== INITIAL PROFILE SWITCH ====
C
C        MODEP : initial profile selector for steady-state simulation
C
      MODEP=3
C
C     ==== INITIAL CURRENT PROFILE SWITCH ====
C
C        MDLJQ : 
C
C           0 : create AJ(NR) profile from experimental Q profile
C           1 : create QP(NR) profile from experimental CURTOT profile
C
      MDLJQ=0
C
C     ==== FLUX SWITCH ====
C
C        MDLFLX :
C
C           0 : use diffusion and convection coefficients
C           1 : use FLUX term made from source files
C
      MDLFLX=0
C
C     ==== Eqs. Selection Parameter ====
C
      MDLEQB=1  ! 0/1 for B_theta
      MDLEQN=0  ! 0/1 for density
      MDLEQT=1  ! 0/1 for heat
      MDLEQU=0  ! 0/1 for rotation
      MDLEQZ=0  ! 0/1 for impurity
      MDLEQ0=0  ! 0/1 for neutral
      MDLEQE=0  ! 0/1/2 for electron density
C               ! 0: electron only, 1: both, 2: ion only
      MDLEOI=0  ! 0/1/2 for electron only or bulk ion only if NSMAX=1
C
C     ==== RADIAL ELECTRIC FIELD SWITCH ====
C
C        0: depend on only pressure gradient
C        1: depend on "0" + toroidal velocity
C        2: depend on "1" + poloidal velocity
C        3: simple circular formula for real geometory
C           (R. Waltz et al, PoP 4 2482 (1997)
C
      MDLER=3
C
C     ==== NCLASS SWITCH ====
C
C        0    : off
C        else : on
C        NSLMAX : the number of species for NCLASS
C                 this parameter is of advantage when you'd like to 
C                 solve only one particle but the other particle
C                 effects are included in the calculation.
C                 i.e. NSLMAX never be less than 2.
      MDNCLS=0
      NSLMAX=NSMAX
C
C     ==== MODERATE TIME EVOLUTION FOR ANOM. TRANSPORT COEFFICIENTS ====
C        0    : off
C        else : multiplyer for TAUK (which is the required time of
C               averaging magnetic surface)
      MDTC=0
C
C     ==== LAPACK ====
C
C        0    : off (using BANDRD for band matrix solver)
C        1    : on  (using DGBTRF and DGBTRS for band matrix solver)
C        else : on  (using DGBSV for band matrix solver)
C
      MDLPCK=0
C
C     ====== INTIALIZE ======
C
      NT=0
      TIME_INT=0.D0
      CNB=1.D0
      SUMPBM=0.D0
      IREAD=0
C
      RETURN
      END
C
C     ***********************************************************
C
C           PARAMETER INPUT
C
C     ***********************************************************
C
      SUBROUTINE TRPARM(MODE,KIN,IERR)
C
C     MODE=0 : standard namelinst input
C     MODE=1 : namelist file input
C     MODE=2 : namelist line input
C
C     IERR=0 : normal end
C     IERR=1 : namelist standard input error
C     IERR=2 : namelist file does not exist
C     IERR=3 : namelist file open error
C     IERR=4 : namelist file read error
C     IERR=5 : namelist file abormal end of file
C     IERR=6 : namelist line input error
C     IERR=7 : unknown MODE
C     IERR=10X : input parameter out of range
C
      INCLUDE 'trcomm.inc'
C
      EXTERNAL TRNLIN,TRPLST
      CHARACTER KIN*(*)
C
    1 CALL TASK_PARM(MODE,'TR',KIN,TRNLIN,TRPLST,IERR)
      IF(IERR.NE.0) RETURN
C
      CALL TRCHEK(IERR)
      NTMAX_SAVE=NTMAX
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100
C
      RETURN
      END
C
C     ****** INPUT NAMELIST ******
C
      SUBROUTINE TRNLIN(NID,IST,IERR)
C
      INCLUDE 'trcomm.inc'
C
      NAMELIST /TR/ RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RHOA,
     &              PA,PZ,PN,PNS,PT,PTS,PNC,PNFE,PNNU,PNNUS,
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,
     &              PROFJ1,PROFJ2,ALP,AD0,AV0,CNP,CNH,CDP,CDH,CDW,
     &              CWEB,CALF,CNB,
     &              MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,MDLTPF,
     &              DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST,
     &              EPSLTR,LMAXTR,CHP,CK0,CK1,CKALFA,CKBETA,CKGUMA,
     &              TPRST,CDW,
     &              MDLST,MDLNF,IZERO,MODELG,NTEQIT,
     &              MDLXP,MDLUF,MDNCLS,MDLWLD,MDLFLX,MDLER,
     &              PNBTOT,PNBR0,PNBRW,PNBENG,PNBRTG,MDLNB,
     &              PECTOT,PECR0,PECRW,PECTOE,PECNPR,MDLEC,
     &              PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,MDLLH,
     &              PICTOT,PICR0,PICRW,PICTOE,PICNPR,MDLIC,
     &              PNBCD,PECCD,PLHCD,PICCD,PBSCD,MDLCD,
     &              PELTOT,PELR0,PELRW,PELRAD,PELVEL,MDLPEL,
     &              PELTIM,PELPAT,KNAMEQ,KNAMTR,KFNLOG,
     &              MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE,
     &              MDLEOI,NSMAX,NSZMAX,NSNMAX,
     &              KUFDEV,KUFDCG,TIME_INT,MODEP,MDNI,MDLJQ,MDTC,MDLPCK
C
      IF(NID.LT.0) THEN
         WRITE(-NID,TR,IOSTAT=IST,ERR=9800)
      ELSE
         READ(NID,TR,IOSTAT=IST,ERR=9800,END=9900)
         NTMAX_SAVE=NTMAX
      ENDIF
      IERR=0
      RETURN
C
 9800 IERR=8
      RETURN
 9900 IERR=9
      RETURN
      END
C
C     ***** INPUT PARAMETER LIST *****
C
      SUBROUTINE TRPLST
C
      WRITE(6,601)
      RETURN
C
  601 FORMAT(' ','# &TR : RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RHOA'/
     &       ' ',8X,'(PA,PZ,PN,PNS,PT,PTS:NSM)'/
     &       ' ',8X,'PNC,PNFE,PNNU,PNNUS'/
     &       ' ',8X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2'/
     &       ' ',8X,'PROFJ1,PROFJ2,ALP'/
     &       ' ',8X,'CK0,CK1,CNP,CNH,CDP,CDH,CDW,CNB'/
     &       ' ',8X,'CWEB,CALF,CKALFA,CKBETA,MDLKNC,MDLTPF'/
     &       ' ',8X,'AD0,CHP,MDLAD,MDLAVK,CKGUMA,MDLKAI,MDLETA,MDLJBS'/
     &       ' ',8X,'DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST'/
     &       ' ',8X,'EPSLTR,LMAXTR,PRST,MDLST,MDLNF,IZERO,PBSCD,MDLCD'/
     &       ' ',8X,'PNBTOT,PNBR0,PNBRW,PNBENG,PNBRTG,PNBCD,MDLNB'/
     &       ' ',8X,'PECTOT,PECR0,PECRW,PECTOE,PECNPR,PECCD,MDLEC'/
     &       ' ',8X,'PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,PLHCD,MDLLH'/
     &       ' ',8X,'PICTOT,PICR0,PICRW,PICTOE,PICNPR,PICCD,MDLIC'/
     &       ' ',8X,'PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM,MDLPEL'/
     &       ' ',8X,'PELTIM,PELPAT,MODELG,NTEQIT'/
     &       ' ',8X,'MDLXP,MDLUF,MDNCLS,MDLWLD,MDLFLX,MDLER'/
     &       ' ',8X,'MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0'/
     &       ' ',8X,'MDLEQE,MDLEOI,NSMAX,NSZMAX,NSNMAX,KUFDEV,KUFDCG'/
     &       ' ',8X,'TIME_INT,MODEP,MDNI,MDLJQ,MDTC,MDLPCK'/
     &       ' ',8X,'KNAMEQ,KNAMTR,KFNLOG')
      END
C
C     ***** CHECK INPUT PARAMETERS *****
C
      SUBROUTINE TRCHEK(IERR)
C
      INCLUDE 'trcomm.inc'
C
      IERR=0
C
      IF(NRMAX.LT.1.OR.NRMAX.GT.NRM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NRMAX'
         WRITE(6,*) '                  NRMAX,NRM =',NRMAX,NRM
         NRMAX=NRM
         IERR=1
      ENDIF
C
      IF(NTMAX.LT.1.OR.NTMAX/NGTSTP.GT.NTM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NTMAX'
         WRITE(6,*) '                  NTMAX,NTM =',NTMAX,NTM
         NTMAX=NTM*NGTSTP
         IERR=1
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
C           VIEW INPUT PARAMETER
C
C     ***********************************************************
C
      SUBROUTINE TRVIEW(ID)
C
      INCLUDE 'trcomm.inc'
C
      WRITE(6,*) '** TRANSPORT **'
      WRITE(6,602) 'MDLEQB',MDLEQB,
     &             'MDLEQN',MDLEQN,
     &             'MDLEQT',MDLEQT,
     &             'MDLEQU',MDLEQU
      WRITE(6,602) 'MDLEQZ',MDLEQZ,
     &             'MDLEQ0',MDLEQ0,
     &             'MDLEQE',MDLEQE,
     &             'MDLEOI',MDLEOI
      WRITE(6,602) 'NSMAX ',NSMAX,
     &             'NSZMAX',NSZMAX,
     &             'NSNMAX',NSNMAX
      WRITE(6,601) 'RR    ',RR,
     &             'RA    ',RA,
     &             'RKAP  ',RKAP,
     &             'RDLT  ',RDLT
      WRITE(6,601) 'RIPS  ',RIPS,
     &             'RIPE  ',RIPE,
     &             'BB    ',BB
C
      WRITE(6,611)
  611 FORMAT(' ','NS',2X,'PA           PZ    PN(E20)  PNS(E20) ',
     &                   'PT(KEV)  PTS(KEV)  PELPAT')
      DO NS=1,NSMAX
         WRITE(6,612) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PT(NS),PTS(NS),
     &                PELPAT(NS)
  612    FORMAT(' ',I2,1PD12.4,0P,F6.1,5F9.4)
      ENDDO
C
      WRITE(6,601) 'PNC   ',PNC,
     &             'PNFE  ',PNFE,
     &             'PNNU  ',PNNU,
     &             'PNNUS ',PNNUS
C
      WRITE(6,601) 'PROFN1',PROFN1,
     &             'PROFT1',PROFT1,
     &             'PROFU1',PROFU1,
     &             'PROFJ1',PROFJ1
C
      WRITE(6,601) 'PROFN2',PROFN2,
     &             'PROFT2',PROFT2,
     &             'PROFU2',PROFU2,
     &             'PROFJ2',PROFJ2
C
      WRITE(6,601) 'ALP(1)',ALP(1),
     &             'ALP(2)',ALP(2),
     &             'ALP(3)',ALP(3),
     &             'PBSCD ',PBSCD
C
      WRITE(6,602) 'MDLKAI',MDLKAI,
     &             'MDLETA',MDLETA,
     &             'MDLAD ',MDLAD,
     &             'MDLAVK',MDLAVK
C
      WRITE(6,602) 'MDLJBS',MDLJBS,
     &             'MDLKNC',MDLKNC,
     &             'MDLTPF',MDLTPF,
     &             'MDNCLS',MDNCLS
C
      WRITE(6,604) 'MDLUF ',MDLUF,
     &             'KUFDEV',KUFDEV,
     &             'KUFDCG',KUFDCG,
     &             'MDNI  ',MDNI
C
      WRITE(6,605) 'MDLJQ ',MDLJQ,
     &             'MDLFLX',MDLFLX,
     &             'MDTC  ',MDTC,
     &             'RHOA  ',RHOA
C
      WRITE(6,602) 'MDLWLD',MDLWLD,
     &             'MDLER ',MDLER,
     &             'MODELG',MODELG,
     &             'NTEQIT',NTEQIT
C
      WRITE(6,601) 'CK0   ',CK0,
     &             'CK1   ',CK1
C
      WRITE(6,601) 'CNP   ',CNP,
     &             'CNH   ',CNH,
     &             'CDP   ',CDP,
     &             'CDH   ',CDH
C
      WRITE(6,601) 'AD0   ',AD0,
     &             'CHP   ',CHP,
     &             'CWEB  ',CWEB,
     &             'CALF  ',CALF
C
      IF((MDLKAI.GE.1.AND.MDLKAI.LT.10).OR.ID.EQ.1)
     &   WRITE(6,601) 'CKALFA',CKALFA,
     &                'CKBETA',CKBETA,
     &                'CKGUMA',CKGUMA
C
      IF((MDLKAI.GE.10.AND.MDLKAI.LT.20).OR.ID.EQ.1) 
     &   WRITE(6,613) CDW(1),CDW(2),CDW(3),CDW(4),
     &                CDW(5),CDW(6),CDW(7),CDW(8)
  613 FORMAT(' ','    AKDW(E) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/
     &       ' ','    AKDW(D) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/
     &       ' ','    AKDW(T) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/
     &       ' ','    AKDW(A) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW')
C
      WRITE(6,601) 'DT    ',DT,
     &             'EPSLTR',EPSLTR,
     &             'TSST  ',TSST,
     &             'TPRST ',TPRST
      WRITE(6,602) 'LMAXTR',LMAXTR,
     &             'NRMAX ',NRMAX,
     &             'NTMAX ',NTMAX,
     &             'NTSTEP',NTSTEP
      WRITE(6,602) 'NGRSTP',NGRSTP,
     &             'NGTSTP',NGTSTP,
     &             'NGPST ',NGPST,
     &             'IZERO ',IZERO
C
      WRITE(6,602) 'MDLST ',MDLST,
     &             'MDLCD ',MDLCD,
     &             'MDLNF ',MDLNF
C
      IF((PNBTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PNBTOT',PNBTOT,
     &                'PNBR0 ',PNBR0,
     &                'PNBRW ',PNBRW,
     &                'PNBENG',PNBENG
         WRITE(6,603) 'MDLNB ',MDLNB,
     &                'PNBRTG',PNBRTG,
     &                'PNBCD ',PNBCD
      ENDIF
C
      IF((PECTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PECTOT',PECTOT,
     &                'PECR0 ',PECR0,
     &                'PECRW ',PECRW,
     &                'PECTOE',PECTOE
         WRITE(6,603) 'MDLEC ',MDLEC,
     &                'PECNPR',PECNPR,
     &                'PECCD ',PECCD
      ENDIF
C
      IF((PLHTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PLHTOT',PLHTOT,
     &                'PLHR0 ',PLHR0,
     &                'PLHRW ',PLHRW,
     &                'PLHTOE',PLHTOE
         WRITE(6,603) 'MDLLH ',MDLLH,
     &                'PLHNPR',PLHNPR,
     &                'PLHCD ',PLHCD
      ENDIF
C
      IF((PICTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PICTOT',PICTOT,
     &                'PICR0 ',PICR0,
     &                'PICRW ',PICRW,
     &                'PICTOE',PICTOE
         WRITE(6,603) 'MDLIC ',MDLIC,
     &                'PICNPR',PICNPR,
     &                'PICCD ',PICCD
      ENDIF
C
      IF((PELTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PELTOT',PELTOT,
     &                'PELR0 ',PELR0,
     &                'PELRW ',PELRW
         WRITE(6,603) 'MDLPEL',MDLPEL,
     &                'PELRAD',PELRAD,
     &                'PELVEL',PELVEL,
     &                'PELTIM',PELTIM
      ENDIF
C
      WRITE(6,601) 'CNB   ',CNB 
      RETURN
C
  601 FORMAT(' ',A6,'=',1PE11.3 :2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
  602 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X   :2X,A6,'=',I7)
  603 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
  604 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',1X,A6,4X:
     &        2X,A6,'=',1X,A6,4X:2X,A6,'=',I7)
  605 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X   :2X,A6,'=',1PE11.3)
      END
C
C     ***********************************************************
C
C           MODEL SELECTOR
C
C     ***********************************************************
C
      SUBROUTINE TR_EQS_SELECT(INIT)
C
      INCLUDE 'trcomm.inc'
      COMMON /TRINS1/ INS
      SAVE NSSMAX
C
C     If INS is zero, all particles designated by NSMAX are fully 
C     calculated. If INS=1, we handle only electrons and bulk ions
C     IF INS=2, three species are employed in the simulation but
C     just one of them is calculated.
C
      NSCMAX=NSMAX+NSZMAX ! the number of charged particles
      NSTMAX=NSMAX+NSZMAX+NSNMAX ! the number of all particles
C
      IF(NSMAX.EQ.1.AND.MDLEOI.EQ.0) MDLEOI=1
      IF(INIT.EQ.0) THEN
         INS=0
      ELSE
         IF(NSSMAX.NE.NSMAX) THEN
            INS=0
         ENDIF
      ENDIF
      NSSMAX=NSMAX
      CALL CHECK_IMPURITY(MDSLCT)
      IF(MDLUF.NE.0) THEN
         IF(MDSLCT.EQ.0) THEN
            IF(NSMAX.EQ.1) THEN
               INS=1
               NSMAX=2
            ENDIF
         ELSE
            IF(NSMAX.EQ.1) INS=2
            NSMAX=3
            PA(3)=12.D0
C            PZ(3)=6.D0
         ENDIF
      ENDIF
C     *** for NCLASS ***
      IF(NSMAX.EQ.1) THEN
         NSLMAX=2
      ELSE
         NSLMAX=NSMAX
      ENDIF
C     ***
C
      IF(MDLEQT.EQ.0) THEN
         NEQMAX=MDLEQB+(MDLEQN+MDLEQT+MDLEQU)*NSMAX
     &         +MDLEQ0*NSNMAX+MDLEQZ*NSZMAX
      ELSEIF(MDLEQT.EQ.1) THEN
         NEQMAX=MDLEQB+NSMAX+(MDLEQT+MDLEQU)*NSMAX
     &         +MDLEQ0*NSNMAX+MDLEQZ*NSZMAX
      ENDIF
C
      IF(MDLEQN.EQ.0.AND.(MDLEQE.EQ.1.OR.MDLEQE.EQ.2)) THEN
         WRITE(6,*) "ERROR! : MDLEQE can be 1 or 2 when MDLEQN is 1."
         STOP
      ENDIF
C
      DO NEQ=1,NEQM
         NSS(NEQ)=-1
         NSV(NEQ)=-1
         NNS(NEQ)=0
         NST(NEQ)=0
      ENDDO
      NEQ=0
      IF(MDLEQB.EQ.1) THEN
         NEQ=NEQ+1
         NSS(NEQ)=0
         NSV(NEQ)=0
      ENDIF
      IND=0
      INDH=0
      INDHD=0
      DO NS=1,NSM
         IF(MDLEQN.EQ.1.OR.MDLEQT.EQ.1) THEN
            IF(MDLEQN.EQ.1) IND=1
            CALL TR_TABLE(NS,NEQ,1,IND,INDH,INDHD)
            IF(IND.EQ.-1) CONTINUE
         ENDIF
         IF(MDLEQT.EQ.1) THEN
            CALL TR_TABLE(NS,NEQ,2,IND,INDH,INDHD)
            IF(IND.EQ.-1) CONTINUE
         ENDIF
         IF(MDLEQU.EQ.1) THEN
            CALL TR_TABLE(NS,NEQ,3,IND,INDH,INDHD)
            IF(IND.EQ.-1) CONTINUE
         ENDIF
      ENDDO
      IF(MDLEQZ.EQ.1) THEN
         IF(NSZMAX.EQ.1) THEN
            NEQ=NEQ+1
            NSS(NEQ  )=5
            NSV(NEQ  )=1
         ELSEIF(NSZMAX.EQ.2) THEN
            NEQ=NEQ+2
            NSS(NEQ-1)=5
            NSV(NEQ-1)=1
            NSS(NEQ  )=6
            NSV(NEQ  )=1
         ENDIF
      ENDIF
      IF(MDLEQ0.EQ.1) THEN
         NEQ=NEQ+2
         NSS(NEQ-1)=7
         NSV(NEQ-1)=1
         NSS(NEQ  )=8
         NSV(NEQ  )=1
      ENDIF
C
      IF(INS.NE.0) THEN
         NEQI=0
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            IF(NSSN.NE.MDLEOI.OR.NSVN.NE.2) THEN
               NEQI=NEQI+1
               NNS(NEQI)=NEQ
            ENDIF
         ENDDO
      ENDIF
C
      NNSC=0
      DO NEQ=1,NEQMAX
         NNSN=NNS(NEQ)
         IF(NNSN.NE.0) THEN
            NNSC=NNSC+1
         ELSE
            GOTO 1000
         ENDIF
      ENDDO
 1000 CONTINUE
      NNSMAX=NNSC
      NEQRMAX=NEQMAX-NNSMAX
C
      NEQS=1
      NEQT=1
      DO NEQ=1,NEQMAX
         NNSN=NNS(NEQ)
         DO NEQ1=1,NEQMAX
            IF(NEQ1.GE.NEQS) THEN
               IF(NEQ1.LT.NNSN.OR.NNSN.EQ.0) THEN
                  NST(NEQ1)=NEQT
                  NEQT=NEQT+1
                  IF(NEQT.GT.NEQRMAX) GOTO 1200
               ELSEIF(NNSN.EQ.NEQ1) THEN
                  NST(NEQ1)=0
                  NEQS=NNSN+1
                  GOTO 1100
               ENDIF
            ENDIF
         ENDDO
 1100    CONTINUE
      ENDDO
 1200 CONTINUE
C
      DO NS=1,8 ! 8 means all particles we should consider.
         IF(PZ(NS).NE.0.D0) THEN
            AMZ(NS)=PA(NS)*AMM/PZ(NS)**2
         ELSE
            AMZ(NS)=0.D0
C     I don't know whether this representation is true or not.
         ENDIF
      ENDDO
C
      WRITE(6,600) 'NEQ','NSS','NSV','NNS','NST'
      DO NEQ=1,NEQMAX
         WRITE(6,610) NEQ,NSS(NEQ),NSV(NEQ),NNS(NEQ),NST(NEQ)
      ENDDO
 600  FORMAT(' ',5(' ',A))
 610  FORMAT(' ',5I4)
C
C     CHECK WHETHER TURBULENT TRANSPORT MODEL HAS OFF-DIAGONAL PARTS
C
      IF(MDLKAI.EQ.61.OR.(MDLKAI.EQ.63.AND.MDLWLD.EQ.1)) THEN
         MDANOM=1
      ELSE
         MDANOM=0
      ENDIF
C
C     |-----------------------|
C     |MDDIAG |NCLASS |MDANOM |
C     |-------|-------|-------|
C     |   0   |   *   |   *   |
C     |   1   |   o   |   *   |
C     |   2   |   *   |   o   |
C     |   3   |   o   |   o   |
C     |-----------------------|
C
      IF(MDNCLS.EQ.0) THEN
         IF(MDANOM.EQ.1) THEN
            MDDIAG=2
         ELSE
            MDDIAG=0
         ENDIF
         DO NR=1,NRMAX
            DO NA=1,5
               DO NS=1,NSM
                  RGFLS(NR,NA,NS)=0.D0
                  RQFLS(NR,NA,NS)=0.D0
               ENDDO
            ENDDO
         ENDDO
      ELSE
         IF(MDANOM.EQ.1) THEN
            MDDIAG=3
         ELSE
            MDDIAG=1
         ENDIF
      ENDIF
C
      RETURN
      END
C
C     ***********************************************************
C
C           SORTER AS MAIN PART OF MODEL SELECTOR
C
C     ***********************************************************
C
      SUBROUTINE TR_TABLE(NS,NEQ,NSW,IND,INDH,INDHD)
C
      INCLUDE 'trcomm.inc'
C
      REM=AME/AMM
      IF(NS.LE.NSMAX) THEN
         IF(ABS(PA(NS)-REM).LE.1.D-10) THEN
C     electron
            NEQ=NEQ+1
            NSS(NEQ)=1
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1.AND.MDLEQE.EQ.0) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 100
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 100           CONTINUE
            ELSEIF(NSW.EQ.2.AND.IND.EQ.1) THEN
               IF(MDLEQE.EQ.0) THEN
                  IF(NSS(1).EQ.0) THEN
                     NNS(1)=2
                  ELSEIF(NSS(1).EQ.1) THEN
                     NNS(1)=1
                  ENDIF
               ELSEIF(MDLEQE.EQ.2) THEN
                  NNS(1)=3
               ENDIF
            ENDIF
         ELSEIF(PA(NS).EQ.1.D0.OR.ABS(PA(NS)-2.D0).LT.0.5D0) THEN
C     regard the particle whose mass is 1.0 as HYDROGEN
C     regard the particle whose mass is between 1.5 and 2.5 as DEUTERIUM
C     If bulk particle is hydrogen, INDH=1
            IF(PA(NS).EQ.1.D0) INDH=1
            NEQ=NEQ+1
C     If plasma is composed of hydrogen and deuterium, INDHD=1
            IF(INDH.EQ.1.AND.ABS(PA(NS)-2.D0).LT.0.5D0) THEN
               NSS(NEQ)=3
               INDHD=1
            ELSE
               NSS(NEQ)=2
            ENDIF
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 200
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 200           CONTINUE
            ENDIF
         ELSEIF(ABS(PA(NS)-3.D0).LT.0.5D0) THEN
C     regard the particle whose mass is between 2.5 and 3.5 as TRITIUM
            NEQ=NEQ+1
            IF(INDHD.EQ.1) THEN
               NSS(NEQ)=4
            ELSE
               NSS(NEQ)=3
            ENDIF
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 300
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 300           CONTINUE
            ENDIF
         ELSEIF(ABS(PA(NS)-4.D0).LT.0.5D0) THEN
C     regard the particle whose mass is between 3.5 and 4.5 as HELIUM
            NEQ=NEQ+1
            NSS(NEQ)=4
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 400
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 400           CONTINUE
            ENDIF
         ELSEIF(ABS(PA(NS)-12.D0).LT.2.D0.AND.NSMAX.EQ.3) THEN
C     regard the particle whose mass is between 10.0 and 14.0 as CARBON
            NEQ=NEQ+1
            NSS(NEQ)=3
            NSV(NEQ)=NSW
            IF(NSW.EQ.2.AND.IND.EQ.0) THEN
               DO NEQI=NEQ-1,NEQMAX
                  NSVN=NSV(NEQI)
                  IF(NSVN.EQ.1) THEN
                     DO NEQII=1,NEQMAX
                        NNSN=NNS(NEQII)
                        IF(NNSN.EQ.0) THEN
                           NNS(NEQII)=NEQI
                           GOTO 500
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
 500           CONTINUE
            ENDIF
         ELSEIF(PA(NS).EQ.0.D0) THEN
            IND=-1
         ENDIF
      ENDIF
C     
      RETURN
      END
