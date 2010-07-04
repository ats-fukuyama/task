!     ***********************************************************

!           INITIALIZE CONSTANTS AND DEFAULT VALUES

!     ***********************************************************

      SUBROUTINE TRINIT

      USE TRCOMM,ONLY : AD0, AEE, ALP, AME, AMM, AV0, BB, CALF, CDH, CDP, CDW, CHP, CK0, CK1, CKALFA, CKBETA, CKGUMA,          &
     &                  CNB, CNH, CNN, CNP, CSPRS, CWEB, DT, EPS0, EPSLTR, IREAD, IZERO, KFNLOG, KNAMEQ, KNAMTR, KUFDCG,       &
     &                  KUFDEV, LMAXTR, MDCD05, MDDW, MDEDGE, MDLAD, MDLAVK, MDLCD, MDLEC, MDLEOI, MDLEQ0, MDLEQB, MDLEQE,     &
     &                  MDLEQN, MDLEQT, MDLEQU, MDLEQZ, MDLER, MDLETA, MDLFLX, MDLIC, MDLJBS, MDLJQ, MDLKAI, MDLKNC, MDLLH,    &
     &                  MDLNB, MDLNF, MDLPCK, MDLPEL, MDLST, MDLTPF, MDLUF, MDLWLD, MDLXP, MDNCLS, MDNI, MDTC, MODELG, MODEP,  &
     &                  NGPST, NGRSTP, NGTSTP, NRMAX, NSLMAX, NSMAX, NSNMAX, NSTM, NSZMAX, NT, NTEQIT, NTMAX, NTSTEP, PA, PBSCD, &
     &                  PECCD, PECNPR, PECR0, PECRW, PECTOE, PECTOT, PELPAT, PELR0, PELRAD, PELRW, PELTIM, PELTOT, PELVEL, PI, &
     &                  PICCD, PICNPR, PICR0, PICRW, PICTOE, PICTOT, PLHNPR, PLHR0, PLHRW, PLHTOE, PLHTOT, PN, PNBCD, PNBENG,  &
     &                  PNBR0, PNBRTG, PNBRW, PNBTOT, PNBVW, PNBVY, PNC, PNFE, PNNU, PNNUS, PNS, PROFJ1, PROFJ2, PROFN1,       &
     &                  PROFN2, PROFT1, PROFT2, PROFU1, PROFU2, PT, PTS, PZ, RA, RDLT, RHOA, RIPE, RIPS, RKAP, RKEV, RMU0, RR, &
     &                  SUMPBM, TIME_INT, TPRST, TSST, VC, VOID, KUFDIR, &
     MDLPR,SYNCABS,SYNCSELF, PU, PUS, PROFNU1, PROFNU2, &
     ELMWID, ELMDUR, ELMNRD, ELMTRD, ELMENH, NSMM, MDLELM
      IMPLICIT NONE
      INTEGER(4) NS, IERR

!     ==== DEVICE PARAMETERS ====

!        RR     : PLASMA MAJOR RADIUS (M)
!        RA     : PLASMA MINOR RADIUS (M)
!        RKAP   : ELIPTICITY OF POLOIDAL CROSS SECTION
!        RDLT   : TRIANGULARITY OF POLOIDAL CROSS SECTION
!        BB     : TOROIDAL MAGNETIC FIELD ON PLASMA AXIS (T)
!        RIPS   : INITIAL VALUE OF PLASMA CURRENT (MA)
!        RIPE   : FINAL VALUE OF PLASMA CURRENT (MA)
!        RHOA   : EDGE OF CALCULATE REGION (NORMALIZED SMALL RADIUS)

      RR      = 3.0D0
      RA      = 1.2D0
      RKAP    = 1.5D0
      RDLT    = 0.0D0
      BB      = 3.D0
      RIPS    = 3.D0
      RIPE    = 3.D0
      RHOA    = 1.D0

!     ==== PLASMA PARAMETERS ====

!        NSMAX  : NUMBER OF MAIN PARTICLE SPECIES (NS=1:ELECTRON)
!        NSZMAX : NUMBER OF IMPURITIES SPECIES
!        NSNMAX : NUMBER OF NEUTRAL SPECIES

!        PA(NS) : ATOMIC NUMBER
!        PZ(NS) : CHARGE NUMBER
!        PN(NS) : INITIAL NUMBER DENSITY ON AXIS (1.E20 M**-3)
!        PNS(NS): INITIAL NUMBER DENSITY ON SURFACE (1.E20 M**-3)
!        PT(NS) : INITIAL TEMPERATURE ON AXIS (KEV)
!        PTS(NS): INITIAL TEMPERATURE ON SURFACE (KEV)
!        PU(NS) : INITIAL TOROIDAL VELOCITY ON AXIS (KEV)
!        PUS(US): INITIAL TOROIDAL VELOCITY ON SURFACE (KEV)

      NSMAX=2
      NSZMAX=0  ! the number of impurities
      NSNMAX=2  ! the number of neutrals, 0 or 2 fixed

      PA(1)   = AME/AMM
      PZ(1)   =-1.D0
      PN(1)   = 0.5D0
      PT(1)   = 1.5D0
      PTS(1)  = 0.05D0
      PNS(1)  = 0.05D0
      PU(1)   = 0.D0
      PUS(1)  = 0.D0

      PA(2)   = 2.D0
      PZ(2)   = 1.D0
      PN(2)   = 0.5D0-2.D-7
      PT(2)   = 1.5D0
      PTS(2)  = 0.05D0
      PNS(2)  = 0.05D0-2.D-8
      PU(2)   = 0.D0
      PUS(2)  = 0.D0

      PA(3)   = 3.D0
      PZ(3)   = 1.D0
      PN(3)   = 1.D-7
      PT(3)   = 1.5D0
      PTS(3)  = 0.05D0
      PNS(3)  = 1.D-8
      PU(3)   = 0.D0
      PUS(3)  = 0.D0

      PA(4)   = 4.D0
      PZ(4)   = 2.D0
      PN(4)   = 1.D-7
      PT(4)   = 1.5D0
      PTS(4)  = 0.05D0
      PNS(4)  = 1.D-8
      PU(4)   = 0.D0
      PUS(4)  = 0.D0

      PA(5)   = 12.D0
      PZ(5)   = 2.D0
      PN(5)   = VOID
      PT(5)   = 0.D0
      PTS(5)  = 0.D0
      PNS(5)  = VOID
      PU(5)   = 0.D0
      PUS(5)  = 0.D0

      PA(6)   = 12.D0
      PZ(6)   = 4.D0
      PN(6)   = VOID
      PT(6)   = 0.D0
      PTS(6)  = 0.D0
      PNS(6)  = VOID
      PU(6)   = 0.D0
      PUS(6)  = 0.D0

      DO NS=7,NSMM
         PA(NS)  = 2.D0
         PZ(NS)  = 0.D0
         PN(NS)  = 0.D0
         PT(NS)  = 0.D0
         PTS(NS) = 0.D0
         PNS(NS) = 0.D0
         PU(NS)  = 0.D0
         PUS(NS) = 0.D0
      END DO

!     ==== IMPURITY PARAMETERS ====

!        PNC    : CARBON DENSITY FACTOR
!        PNFE   : IRON DENSITY FACTOR
!                      COMPARED WITH ITER PHYSICS DESIGN GUIDELINE
!        PNNU   : NEUTRAL NUMBER DENSITY ON AXIS (1.E20 M**-3)
!        PNNUS  :                        ON SURFACE (1.E20 M**-3)

      PNC     = 0.D0
      PNFE    = 0.D0
      PNNU    = 0.D0
      PNNUS   = 0.D0

!     ==== PROFILE PARAMETERS ====

!        PROFN* : PROFILE PARAMETER OF INITIAL DENSITY
!        PROFT* : PROFILE PARAMETER OF INITIAL TEMPERATURE
!        PROFU* : PROFILE PARAMETER OF INITIAL TOROIDAL VELOCITY
!        PROFNU*: PROFILE PARAMETER OF INITIAL NEUTRAL DENSITY
!        PROFJ* : PROFILE PARAMETER OF INITIAL CURRENT DENSITY
!                    (X0-XS)(1-RHO**PROFX1)**PROFX2+XS

!        ALP   : ADDITIONAL PARAMETERS
!           ALP(1): RADIUS REDUCTION FACTOR
!           ALP(2): MASS WEIGHTING FACTOR FOR NC
!           ALP(3): CHARGE WEIGHTING FACTOR FOR NC

      PROFN1 = 2.D0
      PROFN2 = 0.5D0
      PROFT1 = 2.D0
      PROFT2 = 1.D0
      PROFU1 = 2.D0
      PROFU2 = 1.D0
      PROFNU1=12.D0
      PROFNU2= 1.D0
      PROFJ1 =-2.D0
      PROFJ2 = 1.D0

      ALP(1) = 1.0D0
      ALP(2) = 0.D0
      ALP(3) = 0.D0

!     ==== TRANSPORT PARAMETERS ====

!        AV0    : INWARD PARTICLE PINCH FACTOR
!        AD0    : PARTICLE DIFFUSION FACTOR
!        CNP    : COEFFICIENT FOR NEOCLASICAL PARTICLE DIFFUSION
!        CNH    : COEFFICIENT FOR NEOCLASICAL HEAT DIFFUSION
!        CDP    : COEFFICIENT FOR TURBULENT PARTICLE DIFFUSION
!        CDH    : COEFFICIENT FOR TURBLUENT HEAT DIFFUSION
!        CNN    : COEFFICIENT FOR NEUTRAL DIFFUSION
!        CDW(8) : COEFFICIENTS FOR DW MODEL

      AV0    = 0.5D0
      AD0    = 0.5D0

      CNP    = 1.D0
      CNH    = 1.D0
      CDP    = 1.D0
      CDH    = 1.D0
      CNN    = 1.D0
      CDW(1) = 0.04D0
      CDW(2) = 0.04D0
      CDW(3) = 0.04D0
      CDW(4) = 0.04D0
      CDW(5) = 0.04D0
      CDW(6) = 0.04D0
      CDW(7) = 0.04D0
      CDW(8) = 0.04D0

!     ==== TRANSPORT MODEL ====

!        MDLKAI: TURBULENT TRANSPORT MODEL

!   *************************************************************
!   ***  0.GE.MDLKAI.LT.10 : CONSTANT COEFFICIENT MODEL       ***
!   *** 10.GE.MDLKAI.LT.20 : DRIFT WAVE (+ITG +ETG) MODEL     ***
!   *** 20.GE.MDLKAI.LT.30 : REBU-LALLA MODEL                 ***
!   *** 30.GE.MDLKAI.LT.40 : CURRENT-DIFFUSIVITY DRIVEN MODEL ***
!   *** 40.GE.MDLKAI.LT.60 : DRIFT WAVE BALLOONING MODEL      ***
!   ***       MDLKAI.GE.60 : ITG(/TEM, ETG) MODEL ETC         ***
!   *************************************************************

!      ***  MDLKAI.EQ. 0   : CONSTANT*(1+A*rho^2)              ***
!      ***  MDLKAI.EQ. 1   : CONSTANT/(1-A*rho^2)                  ***
!      ***  MDLKAI.EQ. 2   : CONSTANT*(dTi/drho)^B/(1-A*rho^2)     ***
!      ***  MDLKAI.EQ. 3   : CONSTANT*(dTi/drho)^B*Ti^C            ***

!      ***  MDLKAI.EQ. 10  : etac=1                                ***
!      ***  MDLKAI.EQ. 11  : etac=1 1/(1+exp)                      ***
!      ***  MDLKAI.EQ. 12  : etac=1 1/(1+exp) *q                   ***
!      ***  MDLKAI.EQ. 13  : etac=1 1/(1+exp) *(1+q^2)             ***
!      ***  MDLKAI.EQ. 14  : etac=1+2.5*(Ln/RR-0.2) 1/(1+exp)      ***
!      ***  MDLKAI.EQ. 15  : etac=1 1/(1+exp) func(q,eps,Ln)       ***

!      ***  MDLKAI.EQ. 20  : Rebu-Lalla model                      ***

!      ***  MDLKAI.EQ. 30  : CDBM 1/(1+s)                          ***
!      ***  MDLKAI.EQ. 31  : CDBM F(s,alpha,kappaq)                ***
!      ***  MDLKAI.EQ. 32  : CDBM F(s,alpha,kappaq)/(1+WE1^2)      ***
!      ***  MDLKAI.EQ. 33  : CDBM F(s,0,kappaq)                    ***
!      ***  MDLKAI.EQ. 34  : CDBM F(s,0,kappaq)/(1+WE1^2)          ***
!      ***  MDLKAI.EQ. 35  : CDBM (s-alpha)^2/(1+s^2.5)            ***
!      ***  MDLKAI.EQ. 36  : CDBM (s-alpha)^2/(1+s^2.5)/(1+WE1^2)  ***
!      ***  MDLKAI.EQ. 37  : CDBM s^2/(1+s^2.5)                    ***
!      ***  MDLKAI.EQ. 38  : CDBM s^2/(1+s^2.5)/(1+WE1^2)          ***
!      ***  MDLKAI.EQ. 39  : CDBM F2(s,alpha,kappaq,a/R)           ***
!      ***  MDLKAI.EQ. 40  : CDBM F3(s,alpha,kappaq,a/R)/(1+WS1^2) ***

!      ***  MDLKAI.EQ. 60  : GLF23 model                           ***
!      ***  MDLKAI.EQ. 61  : GLF23 (stability enhanced version)    ***
!      ***  MDLKAI.EQ. 62  : IFS/PPPL model                        ***
!      ***  MDLKAI.EQ. 63  : Weiland model                         ***
!      ***  MDLKAI.EQ. 64  : Modified Weiland model                ***


!      ***  MDLKAI.EQ. 130 : CDBM model                            ***
!      ***  MDLKAI.EQ. 131 : CDBM05 model                          ***
!      ***  MDLKAI.EQ. 132 : CDBM model with ExB shear             ***
!      ***  MDLKAI.EQ. 134 : CDBM05 model with ExB shear           ***

!      ***  MDLKAI.EQ. 140 : mBgB (mixed Bonm and gyro-Boahm) model ***
!      ***  MDLKAI.EQ. 141 : mBgB model with suppresion by Tara     ***
!      ***  MDLKAI.EQ. 142 : mBgB model with suppresion by Pacher (EXB)  ***
!      ***  MDLKAI.EQ. 143 : mBgB model with suppresion by Pacher (EXB+SHAR)***

!     +++++ WARNING +++++++++++++++++++++++++++++++++++++++++++
!     +  Parameters below are valid only if MDNCLS /= 0,      +
!     +  that is, one do not use NCLASS,                      +
!     +  otherwise NCLASS automatically calculates all        +
!     +  variables in the following:                          +
!     +                                                       +
!     +    MDLETA: RESISTIVITY MODEL                          +
!     +               1: Hinton and Hazeltine                 +
!     +               2: Hirshman, Hawryluk                   +
!     +               3: Sauter                               +
!     +               4: Hirshman, Sigmar                     +
!     +               else: CLASSICAL                         +
!     +    MDLAD : PARTICLE DIFFUSION MODEL                   +
!     +               1: CONSTANT D                           +
!     +               2: TURBULENT EFFECT                     +
!     +               3: Hinton and Hazeltine                 +
!     +               4: Hinton and Hazeltine w/ TURBULENT    +
!     +               else: NO PARTICLE TRANSPORT             +
!     +    MDLAVK: HEAT PINCH MODEL                           +
!     +               1: Arbitrary amplitude                  +
!     +               2: Arbitrary amplitude w/ pressure dep. +
!     +               3: Hinton and Hazeltine                 +
!     +               else: NO HEAT PINCH                     +
!     +    MDLJBS: BOOTSTRAP CURRENT MODEL                    +
!     +               1-3: Hinton and Hazeltine               +
!     +               4: Hirshman, Sigmar                     +
!     +               5: Sauter                               +
!     +               else: Hinton and Hazeltine              +
!     +    MDLKNC: NEOCLASSICAL TRANSPORT MODEL               +
!     +               0    : Hinton and Hazeltine             +
!     +               else : Chang and Hinton                 +
!     +                                                       +
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!        MDLTPF: TRAPPED PARTICLE FRACTION MODEL
!                   1: Y. R. Lin-Liu and R. L. Miller (numerical)
!                   2: S. P. Hirshman et al.
!                   3: Y. R. Lin-Liu and R. L. Miller (analytic)
!                   4: C  M. N. Rosenbluth et al.
!                   else: Y. B. Kim et al. (default)

      MDLKAI = 31
      MDLETA = 3
      MDLAD  = 3
      MDLAVK = 3
      MDLJBS = 5
      MDLKNC = 1
      MDLTPF = 0

!        MDLWLD : Weiland model mode selector
!            0    : using effective transport coefficients
!            else : using transport coefficients' matrices

      MDLWLD=0

!        MDCD05 : choose either original CDBM or CDBM05 model
!            0    : original CDBM model
!            else : CDBM05 model with the elongation effect

      MDCD05=0

!        MDDW : mode selector for anom. particle transport coefficient
!            you must NOT modify this parameter.
!            0    : if MDDW=0 from start to finish when you choose
!                   a certain transport model (MDLKAI),
!                   you could control a ratio of anomalous particle
!                   transport to total particle transport to manipulate
!                   the factor of AD0.
!            else : this is because you chose MDLKAI=60, 61, or 63
!                   which assign the transport models that can calculate
!                   an anomalous particle transport coefficient
!                   on their own.
      MDDW=0

!     ==== Semi-Empirical Parameter for Anomalous Transport ====

      CHP    = 0.D0
      CK0    = 12.D0 ! for electron
      CK1    = 12.D0 ! for ions
      CWEB   = 1.D0  ! for omega ExB
      CALF   = 1.D0  ! for s-alpha
      CKALFA = 0.D0
      CKBETA = 0.D0
      CKGUMA = 0.D0

!     ==== CONTROL PARAMETERS ====

!        DT     : SIZE OF TIME STEP
!        NRMAX  : NUMBER OF RADIAL MESH POINTS
!        NTMAX  : NUMBER OF TIME STEP
!        NTSTEP : INTERVAL OF SNAP DATA PRINT
!        NGRSTP : INTERVAL OF RADIAL PROFILE SAVE
!        NGTSTP : INTERVAL OF TIME EVOLUTION SAVE
!        NGPST  : ???
!        TSST   : ???

      DT     = 0.01D0
      NRMAX  = 50
      NTMAX  = 100
      NTSTEP = 10
      NGRSTP = 100
      NGTSTP = 2
      NGPST  = 4
      TSST   = 1.D9

!     ==== Convergence Parameter ====

!        EPSLTR : CONVERGENCE CRITERION OF ITERATION
!        LMAXTR : MAXIMUM COUNT OF ITERATION

      EPSLTR = 0.001D0
!      EPSLTR = 1.D99
      LMAXTR = 10

!     ==== SAWTOOTH PARAMETERS ====

!        TPRST  : SAWTOOTH PERIOD (S)
!        MDLST  : SAWTOOTH MODEL TYPE
!                    0:OFF
!                    1:ON
!        IZERO  : SAWTOOTH CRASH TYPE

      TPRST  = 0.1D0
      MDLST  = 0
      IZERO  = 3

!     ==== FUSION REACTION PARAMETERS ====

!        MDLNF  : FUSION REACTION MODEL TYPE
!                    0:OFF
!                    1:ON (DT) without particle source
!                    2:ON (DT) with particle source
!                    3:ON (DT) with NB beam component without particle source
!                    4:ON (DT) with NB beam component particle source
!                    5:ON (DHe3) without particle source
!                    6:ON (DHe3) with particle source

      MDLNF  = 0

!     ==== NBI HEATING PARAMETERS ====

!        PNBTOT : NBI TOTAL INPUT POWER (MW)
!        PNBR0  : RADIAL POSITION OF NBI POWER DEPOSITION (M)
!        PNBRW  : RADIAL WIDTH OF NBI POWER DEPOSITION (M)
!        PNBVY  : VERTICAL POSITION OF NBI (M)
!                 VALID FOR MDLNB=3 or 4
!        PNBVW  : VERTICAL WIDTH OF NBI (M)
!                 VALID FOR MDLNB=3 or 4
!        PNBENG : NBI BEAM ENERGY (keV)
!        PNBRTG : TANGENTIAL RADIUS OF NBI BEAM (M)
!                 VALID FOR MDLNB=3 or 4
!        PNBCD  : CURRENT DRIVE FACTOR
!                 0: off
!                 0 to 1: ratio of v_parallel to v
!        MDLNB  : NBI MODEL TYPE
!                    0:OFF
!                    1:GAUSSIAN (NO PARTICLE SOURCE)
!                    2:GAUSSIAN
!                    3:PENCIL BEAM (NO PARTICLE SOURCE)
!                    4:PENCIL BEAM

      PNBTOT = 0.D0
      PNBR0  = 0.D0
      PNBRW  = 0.5D0
      PNBVY  = 0.D0
      PNBVW  = 0.5D0
      PNBENG = 80.D0
      PNBRTG = 3.D0
      PNBCD  = 1.D0
      MDLNB  = 1

!     ==== ECRF PARAMETERS ====

!        PECTOT : ECRF INPUT POWER (MW)
!        PECR0  : RADIAL POSITION OF POWER DEPOSITION (M)
!        PECRW  : RADIAL WIDTH OF POWER DEPOSITION (M)
!        PECTOE : POWER PARTITION TO ELECTRON
!        PECNPR : PARALLEL REFRACTIVE INDEX
!        PECCD  : CURRENT DRIVE FACTOR
!        MDLEC  : ECRF MODEL

      PECTOT = 0.D0
      PECR0  = 0.D0
      PECRW  = 0.2D0
      PECTOE = 1.D0
      PECNPR = 0.D0
      PECCD  = 0.D0
      MDLEC  = 0

!     ==== LHRF PARAMETERS ====

!        PLHTOT : LHRF INPUT POWER (MW)
!        PLHR0  : RADIAL POSITION OF POWER DEPOSITION (M)
!        PLHRW  : RADIAL WIDTH OF POWER DEPOSITION (M)
!        PLHTOE : POWER PARTITION TO ELECTRON
!        PLHNPR : PARALLEL REFRACTIVE INDEX
!        PLHCD  : CURRENT DRIVE FACTOR
!        MDLLH  : LHRF MODEL

      PLHTOT = 0.D0
      PLHR0  = 0.D0
      PLHRW  = 0.2D0
      PLHTOE = 1.D0
      PLHNPR = 2.D0
      MDLLH  = 0

!     ==== ICRF PARAMETERS ====

!        PICTOT : ICRF INPUT POWER (MW)
!        PICR0  : RADIAL POSITION OF POWER DEPOSITION (M)
!        PICRW  : RADIAL WIDTH OF POWER DEPOSITION (M)
!        PICTOE : POWER PARTITION TO ELECTRON
!        PICNPR : PARALLEL REFRACTIVE INDEX
!        PICCD  : CURRENT DRIVE FACTOR
!        MDLIC  : ICRF MODEL

      PICTOT = 0.D0
      PICR0  = 0.D0
      PICRW  = 0.5D0
      PICTOE = 0.5D0
      PICNPR = 2.D0
      PICCD  = 0.D0
      MDLIC  = 0

!     ==== CURRENT DRIVE PARAMETERS ====

!        PBSCD : BOOTSTRAP CURRENT DRIVE FACTOR
!        MDLCD : CURRENT DRIVE OPERATION MODEL
!                  0: TOTAL PLASMA CURRENT FIXED
!                  1: TOTAL PLASMA CURRENT VARIABLE

      PBSCD  = 1.D0
      MDLCD  = 0

!     ==== PELLET INJECTION PARAMETERS ====

!        MDLPEL : PELLET INJECTION MODEL TYPE
!                    0:OFF  1:GAUSSIAN  2:NAKAMURA  3:HO
!        PELTOT : TOTAL NUMBER OF PARTICLES IN PELLET
!        PELR0  : RADIAL POSITION OF PELLET DEPOSITION (M)
!        PELRW  : RADIAL WIDTH OF PELLET DEPOSITION (M)
!        PELRAD : RADIUS OF PELLET (M)
!        PELVEL : PELLET INJECTION VELOCITY (M/S)
!        PELTIM : TIME FOR PELLET TO BE INJECTED
!        PELPAT : PARTICLE RATIO IN PELLET'

      MDLPEL = 1
      PELTOT = 0.D0
      PELR0  = 0.D0
      PELRW  = 0.5D0
      PELRAD = 0.D0
      PELVEL = 0.D0
      PELTIM = -10.D0

      PELPAT(1) = 1.0D0
      PELPAT(2) = 1.0D0
      DO NS=3,NSMM
         PELPAT(NS) = 0.0D0
      ENDDO

!     ==== RADIATION ====

!     MDLPR  : MODEL OF RADIATION
!               0: PRSUM=PRB+PRL
!               1: PRSUM=PRB+PRL+PRC

!     SYNCABS : fraction of cyclotron radiation absorption by walls
!     SYNCSELF: fraction of x (o) mode reflected as x (o) mode

      MDLPR=0
      SYNCABS=0.2D0
      SYNCSELF=0.95D0

!     ==== EDGE MODEL ====

!      MDEDGE: model parameter for plasma edge model0
!        0 : off
!        1 : simple edge model (set to zero outside NREDGE)

      MDEDGE=0
      CSPRS=0.D0

!     ==== DEVICE NAME AND SHOT NUMBER IN UFILE DATA ====
!        KUFDIR : UFILE DIRECTORY
!        KUFDEV : DEVICE NAME
!        KUFDCG : DISCHARGE NUMBER

      KUFDIR='../../../profiledb/profile_data/'
      KUFDEV='jt60u'
      KUFDCG='21695'

!     ==== FILE NAME ====

!        KNAMEQ: Filename of equilibrium data
!        KNAMTR: Filename of transport data
!        KFNLOG : LOG FILE NAME

      KNAMEQ='eqdata'
      KNAMTR='trdata'
      KFNLOG='trf.log'

!     ==== INTERACTION WITH EQ ====

!        MODELG: 2 : TOROIDAL GEOMETRY
!                3 : READ TASK/EQ FILE
!                5 : READ EQDSK FILE
!                9 : CALCULATE TASK/EQ
!        NTEQIT: STEP INTERVAL OF EQ CALCULATION
!                0 : INITIAL EQUILIBRIUM ONLY

      MODELG=2
      NTEQIT=0

!     ==== INPUT FROM EXPERIMENTAL DATA ====

!        MDLXP :
!           0 : from ufiles
!        else : MDSplus

!        MDLUF :
!           0 : not used
!           1 : time evolution
!           2 : steady state
!           3 : compared with TOPICS

      MDLXP=0
      MDLUF=0

!     ==== IMPURITY TREATMENT ====

!        MDNI  :
!           0 : NSMAX=2, ne=ni
!           1 : Use exp. ZEFFR profile if available
!           2 : Use exp. NM1 (or NM2) profile if available
!           3 : Use exp. NIMP profile if available

      MDNI=0

!     ==== INITIAL PROFILE SWITCH ====

!        MODEP : initial profile selector for steady-state simulation

      MODEP=3

!     ==== INITIAL CURRENT PROFILE SWITCH ====

!        MDLJQ :

!           0 : create QP(NR) profile from experimental CURTOT profile
!           1 : create AJ(NR) profile from experimental Q profile
!           2 : no current with experimental Q profile for helical

      MDLJQ=0

!     ==== FLUX SWITCH ====

!        MDLFLX :

!           0 : use diffusion and convection coefficients
!           1 : use FLUX term made from source files

      MDLFLX=0

!     ==== Eqs. Selection Parameter ====

      MDLEQB=1  ! 0/1 for B_theta
      MDLEQN=0  ! 0/1 for density
      MDLEQT=1  ! 0/1 for heat
      MDLEQU=0  ! 0/1 for rotation
      MDLEQZ=0  ! 0/1 for impurity
      MDLEQ0=0  ! 0/1 for neutral
      MDLEQE=0  ! 0/1/2 for electron density
!               ! 0: electron only, 1: both, 2: ion only
      MDLEOI=0  ! 0/1/2 for electron only or bulk ion only if NSMAX=1
!               ! 0: both, 1: electron, 2: ion

!     ==== RADIAL ELECTRIC FIELD SWITCH ====

!        0: depend on only pressure gradient
!        1: depend on "0" + toroidal velocity
!        2: depend on "1" + poloidal velocity
!        3: simple circular formula for real geometory
!           (R. Waltz et al, PoP 4 2482 (1997)

      MDLER=3

!     ==== NCLASS SWITCH ====

!        0    : off
!        else : on
!        NSLMAX : the number of species for NCLASS
!                 this parameter is of advantage when you'd like to
!                 solve only one particle but the other particle
!                 effects are included in the calculation.
!                 i.e. NSLMAX never be less than 2.
      MDNCLS=0
      NSLMAX=NSMAX

!     ==== ELM MODEL ====
!        MDLELM: 0: no ELM
!                1: simple reduction model: ELMWID, ELMTRD, ELMTRD
!        ELMWID: Factor of Radial width from plasma surface at ELM
!        ELMDUR: Factor of Time duration of ELM
!        ELMNRD: Density reduction factor at ELM for species NS
!        ELMTRD: Temperature reduction factor at ELM for species NS
!        ELMENH: Transport enhancement factor at ELM for species NS

      MDLELM=0
      ELMWID=1.D0
      ELMDUR=1.D0
      DO NS=1,NSMM
         ELMNRD(NS)=1.D0
         ELMTRD(NS)=1.D0
         ELMENH(NS)=1.D0
      END DO

!     ==== MODERATE TIME EVOLUTION FOR ANOM. TRANSPORT COEFFICIENTS ====
!        0    : off
!        else : multiplyer for TAUK (which is the required time of
!               averaging magnetic surface)
      MDTC=0

!     ==== LAPACK ====

!        0    : off (using BANDRD for band matrix solver)
!        1    : on  (using DGBTRF and DGBTRS for band matrix solver)
!        else : on  (using DGBSV for band matrix solver)

      MDLPCK=0

!     ====== INTIALIZE ======

      NT=0
      TIME_INT=0.D0
      CNB=1.D0
      SUMPBM=0.D0
      IREAD=0

      RETURN
      END SUBROUTINE TRINIT

!     ***********************************************************

!           PARAMETER INPUT

!     ***********************************************************

      SUBROUTINE TRPARM(MODE,KIN,IERR)

!     MODE=0 : standard namelinst input
!     MODE=1 : namelist file input
!     MODE=2 : namelist line input

!     IERR=0 : normal end
!     IERR=1 : namelist standard input error
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE
!     IERR=10X : input parameter out of range

      USE TRCOMM, ONLY : NTMAX, NTMAX_SAVE
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: MODE
      CHARACTER(LEN=*),INTENT(IN)::  KIN
      INTEGER(4),INTENT(OUT):: IERR
      EXTERNAL TRNLIN,TRPLST

    1 CALL TASK_PARM(MODE,'TR',KIN,TRNLIN,TRPLST,IERR)
      IF(IERR.NE.0) RETURN

      CALL TRCHEK(IERR)
      NTMAX_SAVE=NTMAX
      IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
      IF(IERR.NE.0) IERR=IERR+100

      RETURN
      END SUBROUTINE TRPARM

!     ****** INPUT NAMELIST ******

      SUBROUTINE TRNLIN(NID,IST,IERR)

      USE TRCOMM, ONLY : RR, RA, RKAP, RDLT, BB, RIPS, RIPE, RHOA, PNC, PNFE, PNNU, PNNUS,            &
     &                   PROFN1, PROFN2, PROFT1, PROFT2, PROFU1, PROFU2, PROFJ1, PROFJ2, ALP, AD0, AV0, CNP, CNH, CDP,          &
     &                   CDH, CNN, CDW, CWEB, CALF, CNB, CSPRS, MDLKAI, MDLETA, MDLAD,                                          &
     &                   MDLAVK, MDLJBS, MDLKNC, MDLTPF, DT, NRMAX, NTMAX, NTSTEP, NGTSTP, NGRSTP, NGPST, TSST,                 &
     &                   EPSLTR, LMAXTR, CHP, CK0, CK1, CKALFA, CKBETA, CKGUMA, TPRST,                                          &
     &                   MDLST, MDLNF, IZERO, MODELG, NTEQIT, MDEDGE, MDLXP, MDLUF, MDNCLS, MDLWLD, MDLFLX, MDLER, MDCD05,      &
     &                   PNBTOT, PNBR0, PNBRW, PNBVY, PNBVW, PNBENG, PNBRTG, MDLNB, PECTOT, PECR0, PECRW, PECTOE, PECNPR, MDLEC,&
     &                   PLHTOT, PLHR0, PLHRW, PLHTOE, PLHNPR, MDLLH, PICTOT, PICR0, PICRW, PICTOE, PICNPR, MDLIC,              &
     &                   PNBCD, PECCD, PLHCD, PICCD, PBSCD, MDLCD, PELTOT, PELR0, PELRW, PELRAD, PELVEL, MDLPEL,                &
     MDLPR, SYNCABS, SYNCSELF,  &
     &                   PELTIM, KNAMEQ, KNAMTR, KFNLOG, MDLEQB, MDLEQN, MDLEQT, MDLEQU, MDLEQZ, MDLEQ0, MDLEQE,        &
     &                   MDLEOI, NSMAX, NSZMAX, NSNMAX, KUFDIR, KUFDEV, KUFDCG, TIME_INT, MODEP, MDNI, MDLJQ,  &
     &                   MDTC, MDLPCK, NTMAX_SAVE
      USE TRCOM3
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NID
      INTEGER(4),INTENT(OUT):: IST, IERR

      NAMELIST /TR/ RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RHOA, &
     &              PA,PZ,PN,PNS,PT,PTS,PNC,PNFE,PNNU,PNNUS, &
     &              PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
     &              PROFJ1,PROFJ2,ALP,AD0,AV0,CNP,CNH,CDP,CDH,CNN,CDW, &
     &              CWEB,CALF,CNB,CSPRS, &
     &              MDLKAI,MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,MDLTPF, &
     &              DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST, &
     &              EPSLTR,LMAXTR,CHP,CK0,CK1,CKALFA,CKBETA,CKGUMA, &
     &              TPRST,CDW, &
     &              MDLST,MDLNF,IZERO,MODELG,NTEQIT,MDEDGE, &
     &              MDLXP,MDLUF,MDNCLS,MDLWLD,MDLFLX,MDLER,MDCD05, &
     &              PNBTOT,PNBR0,PNBRW,PNBVY,PNBVW,PNBENG,PNBRTG,MDLNB, &
     &              PECTOT,PECR0,PECRW,PECTOE,PECNPR,MDLEC, &
     &              PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,MDLLH, &
     &              PICTOT,PICR0,PICRW,PICTOE,PICNPR,MDLIC, &
     &              PNBCD,PECCD,PLHCD,PICCD,PBSCD,MDLCD, &
     &              PELTOT,PELR0,PELRW,PELRAD,PELVEL,MDLPEL, &
                    MDLPR,SYNCABS,SYNCSELF, &
     &              PELTIM,PELPAT,KNAMEQ,KNAMTR,KFNLOG, &
     &              MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE, &
     &              MDLEOI,NSMAX,NSZMAX,NSNMAX, &
     &              KUFDIR,KUFDEV,KUFDCG,TIME_INT,MODEP,MDNI,MDLJQ,MDTC,MDLPCK

      IF(NID.LT.0) THEN
         IF(NSTMAX2 .GE. NSMAX+NSZMAX+NSNMAX) THEN
            CALL COM3COPY1(NSTMAX2,PA,PZ,PN,PNS,PT,PTS,PELPAT)
            WRITE(-NID,TR,IOSTAT=IST,ERR=9800)
         ELSE
            WRITE(6,*) "TRNLIN error NSMAX+NSZMAX+NSNMAX too large"
            WRITE(6,*) "TRNLIN       NSTMAX2,NSMAX,NSZMAX,NSNMAX:", &
                                     NSTMAX2,NSMAX,NSZMAX,NSNMAX
            GOTO 9800
         ENDIF
      ELSE
            CALL COM3COPY1(NSTMAX2,PA,PZ,PN,PNS,PT,PTS,PELPAT)
         READ(NID,TR,IOSTAT=IST,ERR=9800,END=9900)
         NTMAX_SAVE=NTMAX
         IF(NSTMAX2 .GE. NSMAX+NSZMAX+NSNMAX) THEN
            CALL COM3COPY2(NSTMAX2,PA,PZ,PN,PNS,PT,PTS,PELPAT)
         ELSE
            WRITE(6,*) "TRNLIN error NSMAX+NSZMAX+NSNMAX too large"
            WRITE(6,*) "TRNLIN       NSTMAX2,NSMAX,NSZMAX,NSNMAX:",NSTMAX2,NSMAX,NSZMAX,NSNMAX
            GOTO 9800
         ENDIF
      ENDIF
      IERR=0
      RETURN

 9800 IERR=8
      RETURN
 9900 IERR=9
      RETURN
      END SUBROUTINE TRNLIN

      SUBROUTINE COM3COPY1(NSTMAX2,PA_L,PZ_L,PN_L,PNS_L,PT_L,PTS_L,PELPAT_L)
      USE TRCOMM,ONLY : NSTM,PA,PZ,PN,PNS,PT,PTS,PELPAT
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSTMAX2
      REAL(8),DIMENSION(NSTMAX2),INTENT(OUT) :: PA_L,PZ_L,PN_L,PNS_L,PT_L,PTS_L,PELPAT_L
      INTEGER(4) :: NS

      DO NS=1,NSTM
         PA_L(NS)     = PA(NS)
         PZ_L(NS)     = PZ(NS)
         PN_L(NS)     = PN(NS)
         PNS_L(NS)    = PNS(NS)
         PT_L(NS)     = PT(NS)
         PTS_L(NS)    = PTS(NS)
         PELPAT_L(NS) = PELPAT(NS)
      ENDDO

      RETURN
      END SUBROUTINE COM3COPY1


      SUBROUTINE COM3COPY2(NSTMAX2,PA_L,PZ_L,PN_L,PNS_L,PT_L,PTS_L,PELPAT_L)
      USE TRCOMM,ONLY : NSTM,PA,PZ,PN,PNS,PT,PTS,PELPAT
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: NSTMAX2
      REAL(8),DIMENSION(NSTMAX2),INTENT(IN) :: PA_L,PZ_L,PN_L,PNS_L,PT_L,PTS_L,PELPAT_L
      INTEGER(4) :: NS

!      IF(NSTMAX .NE. NSMAX+NSZMAX+NSNMAX) THEN
!        NSTMAX=NSMAX+NSZMAX+NSNMAX
!        DEALLOCATE(PA,PZ,PN,PNS,PT,PTS,PELPAT)
!        ALLOCATE(PA(NSTMAX),PZ(NSTMAX),PN(NSTMAX),PNS(NSTMAX),PT(NSTMAX),PTS(NSTMAX),PELPAT(NSTMAX))
!      ENDIF

      DO NS=1,NSTM
         PA(NS)     = PA_L(NS)
         PZ(NS)     = PZ_L(NS)
         PN(NS)     = PN_L(NS)
         PNS(NS)    = PNS_L(NS)
         PT(NS)     = PT_L(NS)
         PTS(NS)    = PTS_L(NS)
         PELPAT(NS) = PELPAT_L(NS)
      ENDDO

      RETURN
      END SUBROUTINE COM3COPY2


!     ***** INPUT PARAMETER LIST *****

      SUBROUTINE TRPLST

      WRITE(6,601)
      RETURN

  601 FORMAT(' ','# &TR : RR,RA,RKAP,RDLT,BB,RIPS,RIPE,RHOA'/ &
     &       ' ',8X,'(PA,PZ,PN,PNS,PT,PTS:NSM)'/ &
     &       ' ',8X,'PNC,PNFE,PNNU,PNNUS'/ &
     &       ' ',8X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2'/ &
     &       ' ',8X,'PROFJ1,PROFJ2,ALP'/ &
     &       ' ',8X,'CK0,CK1,CNP,CNH,CDP,CDH,CNN,CDW,CNB,CSPRS'/ &
     &       ' ',8X,'CWEB,CALF,CKALFA,CKBETA,MDLKNC,MDLTPF'/ &
     &       ' ',8X,'AD0,CHP,MDLAD,MDLAVK,CKGUMA,MDLKAI,MDLETA,MDLJBS'/ &
     &       ' ',8X,'DT,NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP,NGPST,TSST'/ &
     &       ' ',8X,'EPSLTR,LMAXTR,PRST,MDLST,MDLNF,IZERO,PBSCD,MDLCD'/ &
     &       ' ',8X,'PNBTOT,PNBR0,PNBRW,PNBVY,PNBVW,PNBENG,PNBRTG'/ &
     &       ' ',8X,'PNBCD,MDLNB'/ &
     &       ' ',8X,'PECTOT,PECR0,PECRW,PECTOE,PECNPR,PECCD,MDLEC'/ &
     &       ' ',8X,'PLHTOT,PLHR0,PLHRW,PLHTOE,PLHNPR,PLHCD,MDLLH'/ &
     &       ' ',8X,'PICTOT,PICR0,PICRW,PICTOE,PICNPR,PICCD,MDLIC'/ &
     &       ' ',8X,'PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM,MDLPEL'/ &
     &       ' ',8X,'PELTIM,PELPAT,MDLPR,SYNCABS,SYNCSELF,MODELG,NTEQIT'/&
     &       ' ',8X,'MDEDGE'/ &
     &       ' ',8X,'MDLXP,MDLUF,MDNCLS,MDLWLD,MDLFLX,MDLER,MDCD05'/ &
     &       ' ',8X,'MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0'/ &
     &       ' ',8X,'MDLEQE,MDLEOI,NSMAX,NSZMAX,NSNMAX,KUFDIR,KUFDEV,KUFDCG'/ &
     &       ' ',8X,'TIME_INT,MODEP,MDNI,MDLJQ,MDTC,MDLPCK'/ &
     &       ' ',8X,'KNAMEQ,KNAMTR,KFNLOG')
      END SUBROUTINE TRPLST

!     ***** CHECK INPUT PARAMETERS *****

      SUBROUTINE TRCHEK(IERR)

      USE TRCOMM, ONLY : NGTSTP, NRMAX, NTM, NTMAX
      IMPLICIT NONE
      INTEGER(4), INTENT(OUT):: IERR


      IERR=0

      IF(NRMAX.LT.1) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NRMAX'
         WRITE(6,*) '                  NRMAX =',NRMAX
         IERR=1
      ENDIF

      IF(NTMAX.LT.0.OR.NTMAX/NGTSTP.GT.NTM) THEN
         WRITE(6,*) 'XXX INPUT ERROR : ILLEGAL NTMAX'
         WRITE(6,*) '                  NTMAX,NTM =',NTMAX,NTM
         NTMAX=NTM*NGTSTP
         IERR=1
      ENDIF

      RETURN
      END SUBROUTINE TRCHEK

!     ***********************************************************

!           VIEW INPUT PARAMETER

!     ***********************************************************

      SUBROUTINE TRVIEW(ID)

      USE TRCOMM, ONLY : &
           AD0, ALP, BB, CALF, CDH, CDP, CDW, CHP, CK0, CK1, CKALFA, CKBETA, &
           CKGUMA, CNB, CNH, CNN, CNP, CSPRS, CWEB, DT, EPSLTR, IZERO, &
           KUFDCG, KUFDIR,KUFDEV, LMAXTR, MDCD05, MDEDGE, MDLAD, MDLAVK, &
           MDLCD, MDLEC, MDLEOI, MDLEQ0, MDLEQB, MDLEQE, MDLEQN, MDLEQT, &
           MDLEQU, MDLEQZ, MDLER, MDLETA, MDLFLX, MDLIC, MDLJBS, MDLJQ, &
           MDLKAI, MDLKNC, MDLLH, MDLNB, MDLNF, MDLPEL, MDLPR, SYNCABS, &
           SYNCSELF, MDLST, MDLTPF, MDLUF, MDLWLD, MDNCLS, MDNI, MDTC, &
           MODELG, NGPST, NGRSTP, NGTSTP, NRMAX, NSMAX, NSNMAX, NSZMAX, &
           NTEQIT, NTMAX, NTSTEP, PA, PBSCD, PECCD, PECNPR, PECR0, PECRW, &
           PECTOE, PECTOT, PELPAT, PELR0, PELRAD, PELRW, PELTIM, PELTOT, &
           PELVEL, PICCD, PICNPR, PICR0, PICRW, PICTOE, PICTOT, PLHCD, &
           PLHNPR, PLHR0, PLHRW,  PLHTOE, PLHTOT, PN, PNBCD, PNBENG, PNBR0, &
           PNBRTG, PNBRW, PNBTOT, PNBVW, PNBVY, PNC, PNFE, PNNU, PNNUS, &
           PNS, PROFJ1, PROFJ2, PROFN1, PROFN2, PROFT1, PROFT2, PROFU1, &
           PROFU2, PT, PTS, PZ, RA, RDLT, RHOA, RIPE, RIPS, RKAP, RR, TPRST, &
           TSST
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: ID
      INTEGER(4) :: NS


      WRITE(6,*) '** TRANSPORT **'
      WRITE(6,602) 'MDLEQB',MDLEQB,'MDLEQN',MDLEQN,'MDLEQT',MDLEQT,'MDLEQU',MDLEQU
      WRITE(6,602) 'MDLEQZ',MDLEQZ,'MDLEQ0',MDLEQ0,'MDLEQE',MDLEQE,'MDLEOI',MDLEOI
      WRITE(6,602) 'NSMAX ',NSMAX, 'NSZMAX',NSZMAX,'NSNMAX',NSNMAX
      WRITE(6,601) 'RR    ',RR,    'RA    ',RA,    'RKAP  ',RKAP,  'RDLT  ',RDLT
      WRITE(6,601) 'RIPS  ',RIPS,  'RIPE  ',RIPE,  'BB    ',BB

      WRITE(6,611)
  611 FORMAT(' ','NS',2X,'PA           PZ      PN(E20)  PNS(E20) ','PT(KEV)  PTS(KEV) PELPAT')
      DO NS=1,NSMAX
         WRITE(6,612) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PT(NS),PTS(NS),PELPAT(NS)
  612    FORMAT(' ',I2,1PD12.4,0P,F8.3,5F9.4)
      ENDDO

      WRITE(6,601) 'PNC   ',PNC,   'PNFE  ',PNFE,  'PNNU  ',PNNU,  'PNNUS ',PNNUS
      WRITE(6,601) 'PROFN1',PROFN1,'PROFT1',PROFT1,'PROFU1',PROFU1,'PROFJ1',PROFJ1
      WRITE(6,601) 'PROFN2',PROFN2,'PROFT2',PROFT2,'PROFU2',PROFU2,'PROFJ2',PROFJ2
      WRITE(6,601) 'ALP(1)',ALP(1),'ALP(2)',ALP(2),'ALP(3)',ALP(3),'PBSCD ',PBSCD
      WRITE(6,602) 'MDLKAI',MDLKAI,'MDLETA',MDLETA,'MDLAD ',MDLAD, 'MDLAVK',MDLAVK
      WRITE(6,602) 'MDLJBS',MDLJBS,'MDLKNC',MDLKNC,'MDLTPF',MDLTPF,'MDNCLS',MDNCLS
      WRITE(6,604) 'MDLUF ',MDLUF, 'KUFDEV',KUFDEV,'KUFDCG',KUFDCG,'MDNI  ',MDNI
      WRITE(6,605) 'MDLJQ ',MDLJQ, 'MDLFLX',MDLFLX,'MDTC  ',MDTC,  'RHOA  ',RHOA
      WRITE(6,602) 'MDLWLD',MDLWLD,'MDLER ',MDLER, 'MODELG',MODELG,'NTEQIT',NTEQIT
      WRITE(6,603) 'MDCD05',MDCD05,'CK0   ',CK0,   'CK1   ',CK1
      WRITE(6,603) 'MDEDGE',MDEDGE,'CSPRS ',CSPRS, 'CNN   ',CNN
      WRITE(6,601) 'CNP   ',CNP,   'CNH   ',CNH,   'CDP   ',CDP,   'CDH   ',CDH
      WRITE(6,601) 'AD0   ',AD0,   'CHP   ',CHP,   'CWEB  ',CWEB,  'CALF  ',CALF
      IF((MDLKAI.GE.1.AND.MDLKAI.LT.10).OR.ID.EQ.1) &
     &   WRITE(6,601) 'CKALFA',CKALFA,'CKBETA',CKBETA,'CKGUMA',CKGUMA

      IF((MDLKAI.GE.10.AND.MDLKAI.LT.20).OR.ID.EQ.1)  &
     &   WRITE(6,613) CDW(1),CDW(2),CDW(3),CDW(4),CDW(5),CDW(6),CDW(7),CDW(8)
  613 FORMAT(' ','    AKDW(E) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/ &
     &       ' ','    AKDW(D) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/ &
     &       ' ','    AKDW(T) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/ &
     &       ' ','    AKDW(A) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW')

      WRITE(6,601) 'DT    ',DT,    'EPSLTR',EPSLTR,'TSST  ',TSST,  'TPRST ',TPRST
      WRITE(6,602) 'LMAXTR',LMAXTR,'NRMAX ',NRMAX, 'NTMAX ',NTMAX, 'NTSTEP',NTSTEP
      WRITE(6,602) 'NGRSTP',NGRSTP,'NGTSTP',NGTSTP,'NGPST ',NGPST, 'IZERO ',IZERO
      WRITE(6,602) 'MDLST ',MDLST, 'MDLCD ',MDLCD, 'MDLNF ',MDLNF

      IF((PNBTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PNBTOT',PNBTOT,'PNBR0 ',PNBR0,'PNBRW ',PNBRW,'PNBENG',PNBENG
         WRITE(6,603) 'MDLNB ',MDLNB, 'PNBRTG',PNBRTG,'PNBCD ',PNBCD,'PNBVY ',PNBVY
         WRITE(6,601) 'PNBVW ',PNBVW
      ENDIF

      IF((PECTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PECTOT',PECTOT, 'PECR0 ',PECR0,'PECRW ',PECRW,'PECTOE',PECTOE
         WRITE(6,603) 'MDLEC ',MDLEC,  'PECNPR',PECNPR,'PECCD ',PECCD
      ENDIF

      IF((PLHTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PLHTOT',PLHTOT, 'PLHR0 ',PLHR0,'PLHRW ',PLHRW, 'PLHTOE',PLHTOE
         WRITE(6,603) 'MDLLH ',MDLLH,  'PLHNPR',PLHNPR,'PLHCD ',PLHCD
      ENDIF

      IF((PICTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PICTOT',PICTOT, 'PICR0 ',PICR0, 'PICRW ',PICRW,'PICTOE',PICTOE
         WRITE(6,603) 'MDLIC ',MDLIC,  'PICNPR',PICNPR,'PICCD ',PICCD
      ENDIF

      IF((PELTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PELTOT',PELTOT,'PELR0 ',PELR0,'PELRW ',PELRW
         WRITE(6,603) 'MDLPEL',MDLPEL,'PELRAD',PELRAD,'PELVEL',PELVEL,'PELTIM',PELTIM
      ENDIF

      IF((MDLPR.GE.1).OR.(ID.EQ.1)) THEN
         WRITE(6,623) 'MDLPR   ',MDLPR,   'SYNCABS ',SYNCABS, &
                      'SYNCSELF',SYNCSELF
      ENDIF

      WRITE(6,601) 'CNB   ',CNB
      RETURN

  601 FORMAT(' ',A6,'=',1PE11.3 :2X,A6,'=',1PE11.3: &
     &        2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
  602 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',I7,4X  : &
     &        2X,A6,'=',I7,4X   :2X,A6,'=',I7)
  603 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',1PE11.3: &
     &        2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
  604 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',1X,A6,4X: &
     &        2X,A6,'=',1X,A6,4X:2X,A6,'=',I7)
  605 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',I7,4X  : &
     &        2X,A6,'=',I7,4X   :2X,A6,'=',1PE11.3)
  623 FORMAT(' ',A8,'=',I7,4X   :2X,A8,'=',1PE11.3: &
     &        2X,A8,'=',1PE11.3)
      END SUBROUTINE TRVIEW

!     ***********************************************************

!           MODEL SELECTOR

!     ***********************************************************

      SUBROUTINE TR_EQS_SELECT(INIT)

      USE TRCOMM, ONLY : AMM, AMZ, MDDIAG, MDLEOI, MDLEQ0, MDLEQB, MDLEQE, MDLEQN, MDLEQT, MDLEQU, MDLEQZ, MDLKAI, MDLUF, &
     &                   MDLWLD, MDNCLS, NEA, NEQM, NEQMAX, NEQMAXM, NNS, NREDGE, NRMAX, NSCMAX, NSLMAX, NSM, NSMAX,      &
     &                   NSNMAX, NSS, NST, NSTM, NSTMAX, NSV, NSZMAX, PA, PZ, RGFLS, RQFLS
      USE TRCOM1, ONLY : INS
      IMPLICIT NONE
      INTEGER(4),INTENT(IN) :: INIT
      INTEGER(4) :: IND, INDH, INDHD, MDANOM, MDSLCT, NEQ, NEQ1, NEQI, NEQRMAX, NEQS, NEQT, NNSC, NNSMAX, NNSN, &
     &              NS, NSSN, NSVN
      INTEGER(4), SAVE :: NSSMAX


!     If INS is zero, all particles designated by NSMAX are fully
!     calculated. If INS=1, we handle only electrons and bulk ions
!     IF INS=2, three species are employed in the simulation but
!     just one of them is calculated.

      NSCMAX=NSMAX+NSZMAX          ! the number of charged particles
      NSTMAX=NSMAX+NSZMAX+NSNMAX   ! the number of all particles

      IF(NSMAX.EQ.1.AND.MDLEOI.EQ.0) MDLEOI=1
      IF(INIT.EQ.0) THEN
         INS=0
      ELSE
         IF(NSSMAX.NE.NSMAX) THEN
            INS=0
         ENDIF
      ENDIF
      NSSMAX=NSMAX
      IF(MDLUF.NE.0) THEN
         CALL CHECK_IMPURITY(MDSLCT)
         IF(MDSLCT.EQ.0) THEN
            IF(NSMAX.EQ.1) THEN
               INS=1
               NSMAX=2
            ENDIF
         ELSE
            IF(NSMAX.EQ.1) INS=2
            NSMAX=3
!            PA(3)=12.D0
!            PZ(3)=6.D0
         ENDIF
      ENDIF
!     *** for NCLASS ***
      IF(NSMAX.EQ.1) THEN
         NSLMAX=2
      ELSE
         NSLMAX=NSMAX
      ENDIF
!     ***

!      IF(MDLEQT.EQ.0) THEN
!         NEQMAX=MDLEQB+(MDLEQN+MDLEQT+MDLEQU)*2+MDLEQ0*NSNMAX+MDLEQZ*NSZMAX
!      ELSEIF(MDLEQT.EQ.1) THEN
!         NEQMAX=MDLEQB+2+(MDLEQT+MDLEQU)*2+MDLEQ0*NSNMAX+MDLEQZ*NSZMAX
!      ENDIF
      IF(MDLEQT.EQ.0) THEN
         NEQMAX=MDLEQB+(MDLEQN+MDLEQT+MDLEQU)*NSMAX+MDLEQ0*NSNMAX+MDLEQZ*NSZMAX
      ELSEIF(MDLEQT.EQ.1) THEN
         NEQMAX=MDLEQB+NSMAX+(MDLEQT+MDLEQU)*NSMAX+MDLEQ0*NSNMAX+MDLEQZ*NSZMAX
      ENDIF

      IF(MDLEQN.EQ.0.AND.(MDLEQE.EQ.1.OR.MDLEQE.EQ.2)) THEN
         WRITE(6,*)  'XX TR_EQS_SELECT : MDLEQE can be 1 or 2 when MDLEQN is 1.'
         STOP
      ENDIF

      NSS(1:NEQMAXM)=-1
      NSV(1:NEQMAXM)=-1
      NNS(1:NEQMAXM)=0
      NST(1:NEQMAXM)=0
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

      DO NS=1,8 ! 8 means all particles we should consider.
         IF(PZ(NS).NE.0.D0) THEN
            AMZ(NS)=PA(NS)*AMM/PZ(NS)**2
         ELSE
            AMZ(NS)=0.D0
!     I don't know whether this representation is true or not.
         ENDIF
      ENDDO

      WRITE(6,600) 'NEQ','NSS','NSV','NNS','NST'
      DO NEQ=1,NEQMAX
         WRITE(6,610) NEQ,NSS(NEQ),NSV(NEQ),NNS(NEQ),NST(NEQ)
      ENDDO
 600  FORMAT(' ',5(' ',A))
 610  FORMAT(' ',5I4)

!     *** EQUATION SELECTOR ***

!     Format : NEA(species,equation) for all equations

      NEA(0:NSTM,0:3)=0
      DO NEQ=1,NEQMAX
         NEA(NSS(NEQ),NSV(NEQ))=NEQ
      ENDDO

!     CHECK WHETHER TURBULENT TRANSPORT MODEL HAS OFF-DIAGONAL PARTS

      IF(MDLKAI.EQ.61.OR.(MDLKAI.EQ.63.AND.MDLWLD.EQ.1)) THEN
         MDANOM=1
      ELSE
         MDANOM=0
      ENDIF

!     |-----------------------|
!     |MDDIAG |NCLASS |MDANOM |
!     |-------|-------|-------|
!     |   0   |   *   |   *   |
!     |   1   |   o   |   *   |
!     |   2   |   *   |   o   |
!     |   3   |   o   |   o   |
!     |-----------------------|

      IF(MDNCLS.EQ.0) THEN
         IF(MDANOM.EQ.1) THEN
            MDDIAG=2
         ELSE
            MDDIAG=0
         ENDIF
         RGFLS(1:NRMAX,1:5,1:NSMAX)=0.D0
         RQFLS(1:NRMAX,1:5,1:NSMAX)=0.D0
      ELSE
         IF(MDANOM.EQ.1) THEN
            MDDIAG=3
         ELSE
            MDDIAG=1
         ENDIF
      ENDIF

!     *** GRID POINT OF EDGE REGION ***

      NREDGE=NINT(0.93*NRMAX)

      RETURN
      END SUBROUTINE TR_EQS_SELECT

!     ***********************************************************

!           SORTER AS MAIN PART OF MODEL SELECTOR

!     ***********************************************************

      SUBROUTINE TR_TABLE(NS,NEQ,NSW,IND,INDH,INDHD)

      USE TRCOMM, ONLY : AME, AMM, MDLEQE, NEQMAX, NNS, NSMAX, NSS, NSV, PA
      IMPLICIT NONE
      INTEGER(4),INTENT(IN)   :: NS, NSW
      INTEGER(4),INTENT(INOUT):: NEQ, IND, INDH, INDHD
      INTEGER(4) :: NEQI, NEQII, NNSN, NSVN
      REAL(8)    :: REM

      REM=AME/AMM
      IF(NS.LE.NSMAX) THEN
         IF(ABS(PA(NS)-REM).LE.1.D-10) THEN
!     electron
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
!     regard the particle whose mass is 1.0 as HYDROGEN
!     regard the particle whose mass is between 1.5 and 2.5 as DEUTERIUM
!     If bulk particle is hydrogen, INDH=1
            IF(PA(NS).EQ.1.D0) INDH=1
            NEQ=NEQ+1
!     If plasma is composed of hydrogen and deuterium, INDHD=1
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
!     regard the particle whose mass is between 2.5 and 3.5 as TRITIUM
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
!     regard the particle whose mass is between 3.5 and 4.5 as HELIUM
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
         ELSEIF(ABS(PA(NS)-12.D0).LT.3.D0.AND.NSMAX.EQ.3) THEN
!     regard the particle whose mass is between 9.0 and 15.0 as CARBON
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

      RETURN
      END SUBROUTINE TR_TABLE
