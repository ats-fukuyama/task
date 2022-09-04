! trinit.f90

MODULE trinit

  PRIVATE
  PUBLIC tr_init

CONTAINS

!     ***********************************************************

!           INITIALIZE CONSTANTS AND DEFAULT VALUES

!     ***********************************************************

      SUBROUTINE tr_init

      USE trcomm
      IMPLICIT NONE
      INTEGER NS, NPSC

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
      PN(2)   = 0.5D0
      PT(2)   = 1.5D0
      PTS(2)  = 0.05D0
      PNS(2)  = 0.05D0
      PU(2)   = 0.D0
      PUS(2)  = 0.D0

      PA(3)   = 3.D0
      PZ(3)   = 1.D0
      PN(3)   = 0.D0
      PT(3)   = 1.5D0
      PTS(3)  = 0.05D0
      PNS(3)  = 0.D0
      PU(3)   = 0.D0
      PUS(3)  = 0.D0

      PA(4)   = 4.D0
      PZ(4)   = 2.D0
      PN(4)   = 0.D0
      PT(4)   = 1.5D0
      PTS(4)  = 0.05D0
      PNS(4)  = 0.D0
      PU(4)   = 0.D0
      PUS(4)  = 0.D0

      PA(5)   = 12.D0
      PZ(5)   = 2.D0
      PN(5)   = 0.D0
      PT(5)   = 0.D0
      PTS(5)  = 0.D0
      PNS(5)  = 0.D0
      PU(5)   = 0.D0
      PUS(5)  = 0.D0

      PA(6)   = 12.D0
      PZ(6)   = 4.D0
      PN(6)   = 0.D0
      PT(6)   = 0.D0
      PTS(6)  = 0.D0
      PNS(6)  = 0.D0
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

!     MDLIMP : MODEL IMPURITY TREATMENT with PNC and PNFE
!              0 : PNC and PNFE are not used
!              1 : Initial Impurity density according to ITER PHYS GD for 1.D0
!              2 : Initial Impurity density factor: ANC=PNC*ANE
!              3 : Electron density changes with Te through PZC/PZFE for case 1
!              4 : Electron density changes with Te through PZC/PZFE for case 2

      MDLIMP=0

!        PNC    : CARBON DENSITY FACTOR (1.D0 for Guideline, reduce for low PN)
!        PNFE   : IRON DENSITY FACTOR   (1.D0 for Guideline, reduce for low PN)
!                      COMPARED WITH ITER PHYSICS DESIGN GUIDELINE
!        PNNU   : NEUTRAL NUMBER DENSITY ON AXIS (1.E20 M**-3)    : not used
!        PNNUS  :                        ON SURFACE (1.E20 M**-3) : not used

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

!   *** TRANSPORT MODEL ***

!   model_chi_tb: TURBULENT TRANSPORT MODEL

!        [  0..  9] : CONSTANT COEFFICIENT MODEL
!        [ 10.. 19] : DRIFT WAVE (+ITG +ETG) MODEL
!        [ 20.. 29] : REBU-LALLA MODEL
!        [ 30.. 39] : CURRENT-DIFFUSIVITY DRIVEN MODEL
!        [ 40.. 50] : DRIFT WAVE BALLOONING MODEL
!        [ 60.. 64] : GLF23,IFS/PPPL,Weiland models
!        [130..134] : CDBM model
!        [140..143] : Mixed Boam and gyroBohm model
!        [150..151] : mmm95 model
!        [160..162] : mmm7_1 model (ETG not included)

!      ***  model_chi_tb.EQ. 0   : CONSTANT*(1+A*rho^2)                  ***
!      ***  model_chi_tb.EQ. 1   : CONSTANT/(1-A*rho^2)                  ***
!      ***  model_chi_tb.EQ. 2   : CONSTANT*(dTi/drho)^B/(1-A*rho^2)     ***
!      ***  model_chi_tb.EQ. 3   : CONSTANT*(dTi/drho)^B*Ti^C            ***
!      ***  model_chi_tb.EQ. 4   : PROP. TO CUBIC FUNC. with Bohm  ***
!      ***  model_chi_tb.EQ. 5   : PROP. TO CUBIC FUNC.            ***

!      ***  model_chi_tb.EQ. 10  : etac=1                                ***
!      ***  model_chi_tb.EQ. 11  : etac=1 1/(1+exp)                      ***
!      ***  model_chi_tb.EQ. 12  : etac=1 1/(1+exp) *q                   ***
!      ***  model_chi_tb.EQ. 13  : etac=1 1/(1+exp) *(1+q^2)             ***
!      ***  model_chi_tb.EQ. 14  : etac=1+2.5*(Ln/RR-0.2) 1/(1+exp)      ***
!      ***  model_chi_tb.EQ. 15  : etac=1 1/(1+exp) func(q,eps,Ln)       ***
!      ***  model_chi_tb.EQ. 16  : (model_chi_tb=15) + ZONAL FLOW        ***

!      ***  model_chi_tb.EQ. 20  : Rebu-Lalla model                      ***

!      ***  model_chi_tb.EQ. 30  : CDBM 1/(1+s)                          ***
!      ***  model_chi_tb.EQ. 31  : CDBM F(s,alpha,kappaq)                ***
!      ***  model_chi_tb.EQ. 32  : CDBM F(s,alpha,kappaq)/(1+WE1^2)      ***
!      ***  model_chi_tb.EQ. 33  : CDBM F(s,0,kappaq)                    ***
!      ***  model_chi_tb.EQ. 34  : CDBM F(s,0,kappaq)/(1+WE1^2)          ***
!      ***  model_chi_tb.EQ. 35  : CDBM (s-alpha)^2/(1+s^2.5)            ***
!      ***  model_chi_tb.EQ. 36  : CDBM (s-alpha)^2/(1+s^2.5)/(1+WE1^2)  ***
!      ***  model_chi_tb.EQ. 37  : CDBM s^2/(1+s^2.5)                    ***
!      ***  model_chi_tb.EQ. 38  : CDBM s^2/(1+s^2.5)/(1+WE1^2)          ***
!      ***  model_chi_tb.EQ. 39  : CDBM F2(s,alpha,kappaq,a/R)           ***

!      ***  model_chi_tb.EQ. 60  : GLF23 model                           ***
!      ***  model_chi_tb.EQ. 61  : GLF23 (stability enhanced version)    ***
!      ***  model_chi_tb.EQ. 62  : IFS/PPPL model                        ***
!      ***  model_chi_tb.EQ. 63  : Weiland model                         ***
!      ***  model_chi_tb.EQ. 64  : Modified Weiland model                ***

!      ***  model_chi_tb.EQ. 130 : CDBM model                            ***
!      ***  model_chi_tb.EQ. 131 : CDBM05 model                          ***
!      ***  model_chi_tb.EQ. 132 : CDBM model with ExB shear             ***
!      ***  model_chi_tb.EQ. 134 : CDBM05 model with ExB shear           ***

!      ***  model_chi_tb.EQ. 140 : mBgB (mixed Bonm and gyro-Boahm) model ***
!      ***  model_chi_tb.EQ. 141 : mBgB model suppressed by Tara          ***
!      ***  model_chi_tb.EQ. 142 : mBgB model suppressed by Pacher (EXB)  ***
!      ***  model_chi_tb.EQ. 143 : mBgB model supp. by Pacher (EXB+SHEAR) ***

!      ***  model_chi_tb.EQ. 150 : mmm95 (Multi-Mode tr Model) (no ExB) ***
!      ***  model_chi_tb.EQ. 151 : mmm95 (Multi-Mode tr Model) (with ExB) ***

!      ***  model_chi_tb.EQ. 160 : mmm7_1 (Multi-Mode tr Model) (no ExB) ***
!      ***  model_chi_tb.EQ. 161 : mmm7_1 (Multi-Mode tr Model) (with ExB) ***
!     ------------------------------------------------------------------------
      
!     +    model_dp_tb: turbulent particle diffusion model
!     +               0: None
!     +               1: same as thermal diffusion model (multiplied by factor)
!     +               2: constatn (factor)
!     +    model_vk_tb: turbulent heat pinch model
!     +    model_vp_tb: turbulent particle pinch model
!     +               0: None
!     +               1: Arbitrary amplitude (multiplied by -factor rD/a^2)
!     +    model_chi_nc: neocllasical heat diffusion model
!     +    model_dp_nc : neocllasical particle diffusion model
!     +    model_vk_nc : neocllasical heat pinch model
!     +    model_vp_nc : neocllasical model
!     +               1: nclass
!     +               2: Hinton and Hazeltine                 +
!     +               3: Chang and Hinton
!     +    model_eta: RESISTIVITY MODEL
!     +               1: nclass
!     +               2: Hinton and Hazeltine
!     +               3: Hirshman, Hawryluk
!     +               4: Hirshman, Sigmar
!     +               5: Sauter
!     +    model_bs: BOOTSTRAP CURRENT MODEL
!     +               1: nclass
!     +               2-3: Hinton and Hazeltine
!     +               4: Hirshman, Sigmar
!     +               5: Sauter
!          model_tpfrac
!                   0: Y. B. Kim et al. (default)
!                   1: Y. R. Lin-Liu and R. L. Miller (numerical)
!                   2: S. P. Hirshman et al.
!                   3: Y. R. Lin-Liu and R. L. Miller (analytic)
!                   4: C  M. N. Rosenbluth et al.

      model_chi_tb=31  ! model for turbulent thermal diffusion 
      model_dp_tb=1    ! model for turbulent particle diffusion 
      model_vk_tb=0    ! model for turbulent thermal pinch
      model_vp_tb=0    ! model for turbulent particle pinch
      model_chi_nc=2   ! model for neoclassical therma diffusion
      model_dp_nc=2    ! model for neoclassical particle diffusion
      model_vk_nc=2    ! model for neocllasical thermal pinch
      model_vp_nc=2    ! model for neocllasical particle pinch
      model_eta=2      ! model for resistivity
      model_bs=5       ! model for bootstrap current
      model_tpfrac=0   ! model for trapped particle fraction

      factor_vk_tb=1.D0   ! factor for turbulent thermal pinch
      factor_vk_nc=1.D0   ! factor for neocllasical thermal pinch
      factor_vp_tb=1.D0   ! factor for turbulent particle pinch
      factor_vp_nc=1.D0   ! factor for neocllasical particle pinch
      factor_eta=1.D0     ! factor for resistivity
      factor_bs=1.D0      ! factor for bootstrap current
       
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
!                   a certain transport model (model_chi_tb),
!                   you could control a ratio of anomalous particle
!                   transport to total particle transport to manipulate
!                   the factor of AD0.
!            else : this is because you chose model_chi_tb=60, 61, 63, or 150:159 
!                   which assign the transport models that can calculate
!                   an anomalous particle transport coefficient
!                   on their own.
      MDDW=0

      CNN    = 1.D0
      CDW(1) = 0.04D0
      CDW(2) = 0.04D0
      CDW(3) = 0.04D0
      CDW(4) = 0.04D0
      CDW(5) = 0.04D0
      CDW(6) = 0.04D0
      CDW(7) = 0.04D0
      CDW(8) = 0.04D0

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

!     ==== CONTINUOUS PARTICLE SOURCE PARAMETERS ====

!        MDLPSC : Partcicle SOurce MODEL TYPE
!                    0:OFF  1:GAUSSIAN
!        NPSCMAX: Number of particle sources (Max: NPSCM)
!        PSCTOT : Time rate of particle source [10^{20} /s]
!        PSCR0  : Radial position of particle source [m]
!        PSCRW  : Radial width of particle source [m]
!        NSPSC  : Ion species of particle source

      MDLPSC = 0
      NPSCMAX = 1
      DO NPSC=1,NPSCM
         PSCTOT(NPSC) = 0.D0
         PSCR0(NPSC)  = 0.D0
         PSCRW(NPSC)  = 0.5D0
         NSPSC(NPSC)  = 2
      END DO

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

!        KNAMEQ:  file name of equilibrium data
!        KNAMEQ2: file name of equilibrium data
!        KNAMTR:  file name of transport data
!        KFNLOG : file name of log data
!        KFNTXT : file name of text output
!        KFNCVS : file name of cvs output

      KNAMEQ='eqdata'
      KNAMEQ2=''
      KNAMTR='trdata'
      KFNLOG='tr.log'
      KFNTXT='tr.txt'
      KFNCVS='tr.cvs'

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

      MDLXP=0

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
      END SUBROUTINE tr_init
END MODULE trinit
