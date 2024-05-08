! trinit.f90

MODULE trinit

  PRIVATE
  PUBLIC tr_init

CONTAINS

!     ***********************************************************

!           INITIALIZE CONSTANTS AND DEFAULT VALUES

!     ***********************************************************

      SUBROUTINE tr_init

      USE trcomm_parm
      IMPLICIT NONE
      INTEGER NS,NNB,NEC,NLH,NIC,NPEL,NPSC,nnf

      !  ==== Configuration parameters ====

      !  MODELG: Configuration switch
      !          2 : TOROIDAL GEOMETRY
      !          3 : READ TASK/EQ FILE
      !          5 : READ EQDSK FILE
      !          9 : CALCULATE with TASK/EQ
      !  NTEQIT: STEP INTERVAL OF EQ CALCULATION
      !          0 : INITIAL EQUILIBRIUM ONLY

      MODELG=2
      NTEQIT=0

      !  ==== DEVICE PARAMETERS for MODELG=2 ====

      !  RR     : PLASMA MAJOR RADIUS (m)
      !  RA     : PLASMA MINOR RADIUS (m)
      !  RB     : WALL MINOR RADIUS (m)
      !  RKAP   : ELIPTICITY OF POLOIDAL CROSS SECTION
      !  RDLT   : TRIANGULARITY OF POLOIDAL CROSS SECTION
      !  BB     : TOROIDAL MAGNETIC FIELD ON PLASMA AXIS (T)
      !  RIPS   : INITIAL VALUE OF PLASMA CURRENT (MA)
      !  RIPE   : FINAL VALUE OF PLASMA CURRENT (MA)

      RR      = 3.0D0
      RA      = 1.2D0
      RB      = 1.3D0
      RKAP    = 1.5D0
      RDLT    = 0.0D0
      BB      = 3.D0
      RIPS    = 3.D0
      RIPE    = 3.D0

      !  ==== PLASMA PARAMETERS ====

      !  NSMAX  : NUMBER OF MAIN PARTICLE SPECIES (NS=1:ELECTRON)
      !  NSZMAX : NUMBER OF IMPURITIES SPECIES
      !  NSNMAX : NUMBER OF NEUTRAL SPECIES

      !  KID_NS(NS): Particle name in two characters
      !  ID_NS(NS): -1: electron, 0: neutral, 1: ion, 2: fast ion
      !  NPA(NS): Atomic number (0 for electron)
      !  PM(NS) : Mass NUMBER
      !  PZ(NS) : Charge NUMBER
      !  PN(NS) : INITIAL NUMBER DENSITY ON AXIS (1.E20 M**-3)
      !  PNS(NS): INITIAL NUMBER DENSITY ON SURFACE (1.E20 M**-3)
      !  PT(NS) : INITIAL TEMPERATURE ON AXIS (KEV)
      !  PTS(NS): INITIAL TEMPERATURE ON SURFACE (KEV)
      !  PU(NS) : INITIAL TOROIDAL VELOCITY ON AXIS (KEV)
      !  PUS(NS): INITIAL TOROIDAL VELOCITY ON SURFACE (KEV)

      !  Single ion : NS=1:e 2:D
      !  DT         : NS=1:e 2:D 3:T 4:He4
      !  DD+DD+DHe3 : NS=1:e 2:D 3:T 4:He4 5:H 6:He3

      NSMAX=2
      NSZMAX=0  ! the number of impurities
      NSNMAX=2  ! the number of neutrals, 0 or 2 fixed

      NS_e=1
      KID_NS(1)= ' e'
      ID_NS(1) = -1
      NPA(1)   = 0
      PM(1)    = AME/AMP
      PZ(1)    =-1.D0

      NS_D=2
      KID_NS(2)= ' D'
      ID_NS(2) = 1
      NPA(2)   = 1
      PM(2)    = 2.D0
      PZ(2)    = 1.D0

      NS_T=3
      KID_NS(3)= ' T'
      ID_NS(3) = 1
      NPA(3)   = 1
      PM(3)    = 3.D0
      PZ(3)    = 1.D0

      NS_A=4
      KID_NS(4)= 'He'
      ID_NS(4) = 1
      NPA(4)   = 2
      PM(4)    = 4.D0
      PZ(4)    = 2.D0

      NS_H=5
      KID_NS(5)= ' H'
      ID_NS(5) = 1
      NPA(5)   = 1
      PM(5)    = 1.D0
      PZ(5)    = 1.D0

      NS_He3=6
      KID_NS(6)= 'He'
      ID_NS(6) = 1
      NPA(6)   = 2
      PM(6)    = 3.D0
      PZ(6)    = 2.D0

      NS_C=7
      KID_NS(7)= ' C'
      ID_NS(7) = 1
      NPA(7)   = 6
      PM(7)    = 12.D0
      PZ(7)    = 6.D0

      NS_Fe=8
      KID_NS(8)= 'Fe'
      ID_NS(8) = 1
      NPA(8)   = 26
      PM(8)    = 56.D0
      PZ(8)    = 26.D0

      DO NS=9,NSM
         KID_NS(NS)= ' H'
         ID_NS(NS) = 1
         NPA(NS)   = 1
         PM(NS)    = 1.D0
         PZ(NS)    = 1.D0
      END DO

      DO NS=1,2
         PN(NS)   = 0.5D0
         PNS(NS)  = 0.05D0
         PT(NS)   = 1.5D0
         PTS(NS)  = 0.05D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
      END DO

      DO NS=3,NSM
         PN(NS)   = 0.D0
         PNS(NS)  = 0.0D0
         PT(NS)   = 1.5D0
         PTS(NS)  = 0.05D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
      END DO

      !    ==== PROFILE PARAMETERS ====

      !    model_prof: profile model
      !    knam_prof: profile data file name
      
      !    PROFN* : PROFILE PARAMETER OF INITIAL DENSITY
      !    PROFT* : PROFILE PARAMETER OF INITIAL TEMPERATURE
      !    PROFU* : PROFILE PARAMETER OF INITIAL TOROIDAL VELOCITY
      !    PROFNU*: PROFILE PARAMETER OF INITIAL NEUTRAL DENSITY
      !    PROFJ* : PROFILE PARAMETER OF INITIAL CURRENT DENSITY
      !             (X0-XS)(1-RHO**PROFX1)**PROFX2+XS

      !    ALP   : ADDITIONAL PARAMETERS
      !       ALP(1): RADIUS REDUCTION FACTOR
      !       ALP(2): MASS WEIGHTING FACTOR FOR NC
      !       ALP(3): CHARGE WEIGHTING FACTOR FOR NC
      !       ALP(4): ADDW factor for D
      !       ALP(5): ADDW factor for T
      !       ALP(6): ADDW factor for He

      model_prof=0
      knam_prof='prof.data'

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
      ALP(4) = 1.0D0
      ALP(5) = 1.0D0
      ALP(6) = 1.0D0
      
      !  ====== Fixed profile =====      

      !  model_nfixed: for density profile
      !  model_tfixed: for temperature profile
      !     0: no fixed profile      
      !     1: fixed profile (x=n for density, x=t for temperatrue) 
      !     2: fixed profile for rho_min_xfixed <= rho <= rho_max_xfixed
      !     read parameters from file 'xprof_coef_data'
      !        ntime_xfixed_max:         number of time points
      !        ndata_xfixed_max:         number of coefficients
      !        rho_min_xfixed:           rho minimum of fixed profile
      !        rho_max_xfixed:           rho maximum of fixed profile
      !        time_xfixed(ntime):       start time of fixed profile   
      !        coef_xfixed(ndata,ntime): coefficients of fixed profile
      !        f(nr)=coef(0) &
      !             +0.5D0*coef(1) &
      !             *(tanh((1.D0-coef(2)*coef(3)-rho(nr))/coef(3))+1.D0) &
      !             +coef(4)*(1.D0-rho(nr)*rho(nr))**coef(5) &
      !             +0.5D0*coef(8)*(1.D0-erf((rho(nr)-coef(9))
      !              /SQRT(2.D0*coef(10))))
      !
      !  knam_nfixed: density fixed profile file name
      !  knam_tfixed: temperature fixed profile file name
      
      model_nfixed=0
      model_tfixed=0
      knam_nfixed='nprof_coef_data'
      knam_tfixed='tprof_coef_data'


      !  ==== IMPURITY ans neutral PARAMETERS ====
      
      !  MDLIMP : MODEL IMPURITY TREATMENT with PNC and PNFE
      !       0 : PNC and PNFE are not used
      !       1 : Initial Impurity density according to ITER PHYS GD
      !       2 : Initial Impurity density factor: ANC=PNC*ANE
      !       3 : n_e changes with Te through PZC/PZFE for case 1
      !       4 : n_e changes with Te through PZC/PZFE for case 2

      !  MDLNI  :
      !       0 : NSMAX=2, ne=ni
      !       1 : Use exp. ZEFFR profile if available
      !       2 : Use exp. NM1 (or NM2) profile if available
      !       3 : Use exp. NIMP profile if available


      !  PNC    : CARBON DENSITY FACTOR (1.D0 for Guideline)
      !  PNFE   : IRON DENSITY FACTOR   (1.D0 for Guideline,)
      !                COMPARED WITH ITER PHYSICS DESIGN GUIDELINE
      !  PNNU   : NEUTRAL NUMBER DENSITY ON AXIS (1.E20 M**-3)    : not used
      !  PNNUS  :                        ON SURFACE (1.E20 M**-3) : not used

      MDLIMP=0
      MDLNI=0

      PNC     = 0.D0
      PNFE    = 0.D0

      PNNU    = 0.D0
      PNNUS   = 0.D0

      !  ==== TRANSPORT MODEL PARAMETERS ====

      !  MDLKAI: TURBULENT TRANSPORT MODEL

      !  ***   0. . 9 : CONSTANT COEFFICIENT MODEL       ***
      !  ***  10.. 19 : DRIFT WAVE (+ITG +ETG) MODEL     ***
      !  ***  20.. 29 : REBU-LALLA MODEL                 ***
      !  ***  30.. 39 : CURRENT-DIFFUSIVITY DRIVEN MODEL ***
      !  ***  40.. 49 : DRIFT WAVE BALLOONING MODEL      ***
      !  ***  60.. 64 : CLF23,IFS/PPPL,Weiland models    ***
      !  *** 130..134 : CDBM model                       ***
      !  *** 140..143 : Mixed Boam and gyroBohm model    ***
      !  *** 150..151 : mmm95 model                      ***
      !  *** 160..162 : mmm7_1 model (ETG not included)  ***

      !  ***  MDLKAI.EQ. 0   : CONSTANT*(1+A*rho^2)
      !                   1   : CONSTANT/(1-A*rho^2)
      !                   2   : CONSTANT*(dTi/drho)^B/(1-A*rho^2)
      !                   3   : CONSTANT*(dTi/drho)^B*Ti^C

      !  ***  MDLKAI.EQ. 10  : etac=1
      !                  11  : etac=1 1/(1+exp)
      !                  12  : etac=1 1/(1+exp) *q
      !                  13  : etac=1 1/(1+exp) *(1+q^2)
      !                  14  : etac=1+2.5*(Ln/RR-0.2) 1/(1+exp)
      !                  15  : etac=1 1/(1+exp) func(q,eps,Ln)
      !                  16  : (MDLKAI=15) + ZONAL FLOW

      !  ***  MDLKAI.EQ. 20  : Rebu-Lalla model

      !  ***  MDLKAI.EQ. 30  : CDBM 1/(1+s)
      !                  31  : CDBM F(s,alpha,kappaq)
      !                  32  : CDBM F(s,alpha,kappaq)/(1+WE1^2)
      !                  33  : CDBM F(s,0,kappaq)
      !                  34  : CDBM F(s,0,kappaq)/(1+WE1^2)
      !                  35  : CDBM (s-alpha)^2/(1+s^2.5)
      !                  36  : CDBM (s-alpha)^2/(1+s^2.5)/(1+WE1^2)
      !                  37  : CDBM s^2/(1+s^2.5)
      !                  38  : CDBM s^2/(1+s^2.5)/(1+WE1^2)
      !                  39  : CDBM F2(s,alpha,kappaq,a/R)
      !                  40  : CDBM F3(s,alpha,kappaq,a/R)/(1+WS1^2)

      !  ***  MDLKAI.EQ. 60  : GLF23 model
      !                  61  : GLF23 (stability enhanced version)
      !                  62  : IFS/PPPL model
      !                  63  : Weiland model
      !                  64  : Modified Weiland model


      !  ***  MDLKAI.EQ. 130 : CDBM model (subroutine tr_cdbm)
      !                  131 : CDBM05 model
      !                  132 : CDBM model with ExB shear
      !                  134 : CDBM05 model with ExB shear

      !  ***  MDLKAI.EQ. 140 : mBgB (mixed Bonm and gyro-Boahm) model
      !                  141 : mBgB model with suppresion by Tara
      !                  142 : mBgB model with suppresion by Pacher (EXB)
      !                  143 : mBgB model with suppresion by Pacher (EXB+SHAR)

      !  ***  MDLKAI.EQ. 150 : mmm95 (Multi-Mode transport Model) (no ExB)
      !                  151 : mmm95 (Multi-Mode transport Model) (with ExB)

      !  ***  MDLKAI.EQ. 160 : mmm7_1 (Multi-Mode transport Model) (no ExB)
      !                  161 : mmm7_1 (Multi-Mode transport Model) (with ExB)

      !     +++++ WARNING +++++++++++++++++++++++++++++++++++++++++++
      !     +  Parameters below are valid only if MDLNCL /= 0,      +
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
      !     +               1: CONSTANT D with pinch AV0            +
      !     +               2: TURBULENT EFFECT with pinch AV0      +
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
      !     +               1    : Hinton and Hazeltine             +
      !     +               else : Chang and Hinton                 +
      !     +                                                       +
      !     +    MDLTPF: TRAPPED PARTICLE FRACTION MODEL
      !     +               1: Y. R. Lin-Liu and R. L. Miller (numerical)
      !     +               2: S. P. Hirshman et al.
      !     +               3: Y. R. Lin-Liu and R. L. Miller (analytic)
      !     +               4: C  M. N. Rosenbluth et al.
      !     +               else: Y. B. Kim et al. (default)
      !     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      MDLKAI = 31
      MDLETA = 3
      MDLAD  = 3
      MDLAVK = 3
      MDLJBS = 5
      MDLKNC = 1
      MDLTPF = 0

      !     ==== NCLASS SWITCH ====

      !     MDLNCL   0    : off
      !              else : on
      !     NSLMAX : the number of species for NCLASS
      !              this parameter is of advantage when you'd like to
      !              solve only one particle but the other particle
      !              effects are included in the calculation.
      !              i.e. NSLMAX never be less than 2.

      MDLNCL=0
      NSLMAX=NSMAX

      !     ==== Turbulence model SWITCH ====

      !   MDLWLD  : Weiland model mode selector
      !     0     : using effective transport coefficients
      !     else  : using transport coefficients' matrices
      !   MDLCD05 : choose either original CDBM or CDBM05 model
      !      0    : original CDBM model
      !      else : CDBM05 model with the elongation effect
      !   MDLDW   : mode selector for anom. particle transport coefficient
      !             you must NOT modify this parameter.
      !      0    : if MDLDW=0 from start to finish when you choose
      !             a certain transport model (MDLKAI),
      !             you could control a ratio of anomalous particle
      !             transport to total particle transport to manipulate
      !             the factor of AD0.
      !      else : because you chose MDLKAI=60, 61, 63, or 150:159 
      !             which assign the transport models that can calculate
      !             an anomalous particle transport coefficient
      !             on their own.

      MDLWLD=0
      MDLCD05=0
      MDLDW=0

      !  ==== MODERATE TIME EVOLUTION FOR ANOM. TRANSPORT COEFFICIENTS ====
      !           MDLTC= 0    : off
      !               else : multiplyer for TAUK (which is the required time of
      !                      averaging magnetic surface)

      MDLTC=0

      !  ==== TRANSPORT PARAMETERS ====

      !  AD0    : PARTICLE DIFFUSION FACTOR
      !  AV0    : INWARD PARTICLE PINCH FACTOR
      
      !  CNP    : COEFFICIENT FOR NEOCLASICAL PARTICLE DIFFUSION
      !  CNH    : COEFFICIENT FOR NEOCLASICAL HEAT DIFFUSION
      !  CDP    : COEFFICIENT FOR TURBULENT PARTICLE DIFFUSION
      !  CDH    : COEFFICIENT FOR TURBLUENT HEAT DIFFUSION
      !  CNN    : COEFFICIENT FOR NEUTRAL DIFFUSION
      !  CDW(8) : COEFFICIENTS FOR DW MODEL

      AD0    = 0.5D0
      AV0    = 0.0D0

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

      !  ==== Semi-Empirical Parameter for Anomalous Transport ====

      CHP    = 0.D0
      CK0    = 12.D0 ! for electron
      CK1    = 12.D0 ! for ions
      CWEB   = 1.D0  ! for omega ExB
      CALF   = 1.D0  ! for s-alpha
      CKALFA = 0.D0
      CKBETA = 0.D0
      CKGUMA = 0.D0

      !  ==== RADIAL ELECTRIC FIELD SWITCH ====

      !  MDLER: 0: depend on only pressure gradient
      !         1: depend on "0" + toroidal velocity
      !         2: depend on "1" + poloidal velocity
      !         3: simple circular formula for real geometory
      !              (R. Waltz et al, PoP 4 2482 (1997)

      MDLER=3

      !  ==== EDGE MODEL ====

      !  MDLEDGE: model parameter for plasma edge model0
      !      0 : off
      !      1 : simple edge model (set factor to CSPRS outside NREDGE)
      !  RHOA  : EDGE OF CALCULATE REGION (NORMALIZED SMALL RADIUS)

      MDLEDGE=0
      NREDGE=0
      CSPRS=0.D0
      RHOA=1.D0

      !   ==== ELM MODEL ====
      
      !     MDLELM: 0: no ELM
      !             1: simple reduction model: ELMWID, ELMTRD, ELMTRD
      !     ELMWID: Factor of Radial width from plasma surface at ELM
      !     ELMDUR: Factor of Time duration of ELM
      !     ELMNRD: Density reduction factor at ELM for species NS
      !     ELMTRD: Temperature reduction factor at ELM for species NS
      !     ELMENH: Transport enhancement factor at ELM for species NS

      MDLELM=0
      ELMWID=1.D0
      ELMDUR=1.D0
      DO NS=1,NSM
         ELMNRD(NS)=1.D0
         ELMTRD(NS)=1.D0
         ELMENH(NS)=1.D0
      END DO

      !  ==== FUSION REACTION PARAMETERS ====

      !  nnfmax : number of fusion product id
      !  model_nnf  : FUSION REACTION MODEL TYPE
      !        0:OFF
      !        1:ON He4 (DT) without particle source
      !        2:ON He4 (DT) with particle source
      !        3:ON He4 (DT) with NB beam component without particle source
      !        4:ON He4 (DT) with NB beam component particle source
      !       11:ON H   (DD) without particle source
      !       12:ON H   (DD) with particle source
      !       13:ON H   (DD) with NB beam component without particle source
      !       14:ON H   (DD) with NB beam component particle source
      !       21:ON D   (DD) without particle source
      !       22:ON D   (DD) with particle source
      !       23:ON D   (DD) with NB beam component without particle source
      !       24:ON D   (DD) with NB beam component particle source
      !       31:ON He3 (DD) without particle source
      !       32:ON He3 (DD) with particle source
      !       33:ON He3 (DD) with NB beam component without particle source
      !       34:ON He3 (DD) with NB beam component particle source
      !       41:ON H   (DHe3) without particle source
      !       42:ON H   (DHe3) with particle source
      !       43:ON H   (DHe3) with NB beam component without particle source
      !       44:ON H   (DHe3) with NB beam component particle source
      !       51:ON He4 (DHe3) without particle source
      !       52:ON He4 (DHe3) with particle source
      !       53:ON He4 (DHe3) with NB beam component without particle source
      !       54:ON He4 (DHe3) with NB beam component particle source

      nnfmax=0
      DO nnf=1,nnfm
         model_nnf(nnf)  = 1
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

      !  *** Source parameters ***

      !  ==== NBI HEATING PARAMETERS ====

      !  NNBMAX       : number of NBI beam (MAX=nnbm)
      !  model_nnb(nnbm)  : NBI MODEL TYPE
      !       0 : OFF
      !       1 : GAUSSIAN (NO PARTICLE SOURCE)
      !       2 : GAUSSIAN
      !       3 : PENCIL BEAM (NO PARTICLE SOURCE)
      !       4 : PENCIL BEAM
      !  PNBIN(nnbm) : NBI TOTAL INPUT POWER (MW)
      !  PNBR0(nnbm)  : RADIAL POSITION OF NBI POWER DEPOSITION (M)
      !  PNBRW(nnbm)  : RADIAL WIDTH OF NBI POWER DEPOSITION (M)
      !  PNBVY(nnbm)  : VERTICAL POSITION OF NBI (M)
      !                 VALID FOR MDLNB=3 or 4
      !  PNBVW(nnbm)  : VERTICAL WIDTH OF NBI (M)
      !                 VALID FOR MDLNB=3 or 4
      !  PNBENG(nnbm) : NBI BEAM ENERGY (keV)
      !  PNBRTG(nnbm) : TANGENTIAL RADIUS OF NBI BEAM (M)
      !                 VALID FOR MDLNB=3 or 4
      !  PNBCD(nnbm)  : CURRENT DRIVE FACTOR
      !             0 : off
      !             0 to 1: ratio of v_parallel to v
      !  ns_nnb(nnbm)  : Particle species number
      !  nrmax_nnb(nnbm): number of division of NB profile      

      NNBMAX=0
      DO NNB=1,NNBM
         model_nnb(nnb) = 1    ! 
         ns_nnb(nnb)    = 2
         nrmax_nnb(nnb) = 10
         PNBIN(NNB)  = 0.D0
         PNBR0(NNB)  = 0.D0
         PNBRW(NNB)  = 0.5D0
         PNBVY(NNB)  = 0.D0
         PNBVW(NNB)  = 0.5D0
         PNBENG(NNB) = 80.D0
         PNBRTG(NNB) = 3.D0
         PNBCD(NNB)  = 1.D0
      END DO

      !  ==== ECRF PARAMETERS ====

      !  NECMAX       : number of EC source   (MAX=necm)
      !  MDLEC(necm)  : ECRF MODEL
      !  PECIN(necm)  : ECRF INPUT POWER (MW)
      !  PECR0(necm)  : RADIAL POSITION OF POWER DEPOSITION (M)
      !  PECRW(necm)  : RADIAL WIDTH OF POWER DEPOSITION (M)
      !  PECTOE(necm) : POWER PARTITION TO ELECTRON
      !  PECNPR(necm) : PARALLEL REFRACTIVE INDEX
      !  PECCD(necm)  : CURRENT DRIVE FACTOR

      NECMAX=1
      DO NEC=1,NECM
         MDLEC(NEC)  = 0
         PECIN(NEC)  = 0.D0
         PECR0(NEC)  = 0.D0
         PECRW(NEC)  = 0.2D0
         PECTOE(NEC) = 1.D0
         PECNPR(NEC) = 0.D0
         PECCD(NEC)  = 0.D0
      END DO

      !  ==== LHRF PARAMETERS ====

      !  NLHMAX       : number of LH source   (MAX=nlhm)
      !  MDLLH(nlhm)  : LHRF MODEL
      !  PLHIN(nlhm)  : LHRF INPUT POWER (MW)
      !  PLHR0(nlhm)  : RADIAL POSITION OF POWER DEPOSITION (M)
      !  PLHRW(nlhm)  : RADIAL WIDTH OF POWER DEPOSITION (M)
      !  PLHTOE(nlhm) : POWER PARTITION TO ELECTRON
      !  PLHNPR(nlhm) : PARALLEL REFRACTIVE INDEX
      !  PLHCD(nlhm)  : CURRENT DRIVE FACTOR

      NLHMAX=1
      DO NLH=1,NLHM
         MDLLH(NLH)  = 0
         PLHIN(NLH)  = 0.D0
         PLHR0(NLH)  = 0.D0
         PLHRW(NLH)  = 0.2D0
         PLHTOE(NLH) = 1.D0
         PLHNPR(NLH) = 2.D0
      END DO

      !  ==== ICRF PARAMETERS ====

      !  NICMAX       : number of IC source   (MAX=nicm)
      !  MDLIC(nicm)  : ICRF MODEL
      !  PICIN(nicm)  : ICRF INPUT POWER (MW)
      !  PICR0(nicm)  : RADIAL POSITION OF POWER DEPOSITION (M)
      !  PICRW(nicm)  : RADIAL WIDTH OF POWER DEPOSITION (M)
      !  PICTOE(nicm) : POWER PARTITION TO ELECTRON
      !  PICNPR(nicm) : PARALLEL REFRACTIVE INDEX
      !  PICCD(nicm)  : CURRENT DRIVE FACTOR

      NICMAX=1
      DO NIC=1,NICM
         MDLIC(NIC)  = 0
         PICIN(NIC)  = 0.D0
         PICR0(NIC)  = 0.D0
         PICRW(NIC)  = 0.5D0
         PICTOE(NIC) = 0.5D0
         PICNPR(NIC) = 2.D0
         PICCD(NIC)  = 0.D0
      END DO

      !  ==== PELLET INJECTION PARAMETERS ====

      !  NPELMAX       : number of pellet source   (MAX=npelm)
      !  MDLPEL(npelm) : PELLET INJECTION MODEL TYPE
      !              0:OFF  1:GAUSSIAN  2:NAKAMURA  3:HO
      !  PELTOT(npelm) : TOTAL NUMBER OF PARTICLES IN PELLET
      !  PELR0(npelm)  : RADIAL POSITION OF PELLET DEPOSITION (M)
      !  PELRW(npelm)  : RADIAL WIDTH OF PELLET DEPOSITION (M)
      !  PELRAD(npelm) : RADIUS OF PELLET (M)
      !  PELVEL(npelm) : PELLET INJECTION VELOCITY (M/S)
      !  PELTIM(npelm) : TIME FOR PELLET TO BE INJECTED
      !  PELPAT(npelm) : PARTICLE RATIO IN PELLET'
      !  pellet_time_start(npelm) : start time of pellet injection repeat (s)
      !  pellet_time_interval(npelm) : interval of pellet injection repeat (s)
      !  number_of_pellet_repeat(npelm) : max. repeat count of pellet injection

      NPELMAX=1
      DO npel=1,npelm
         PELIN(npel) = 0.D0
         PELR0(npel) = 0.D0
         PELRW(npel) = 0.5D0
         PELRAD(npel)= 0.D0
         PELVEL(npel)= 0.D0
         PELTIM(npel)= -10.D0
         pellet_time_start(npel)=0.D0
         pellet_time_interval(npel)=1.D0
         DO NS=1,NSM
            PELPAT(NS,npel) = 0.0D0
         END DO
         PELPAT(2,npel) = 1.0D0
         MDLPEL(npel)= 1
         number_of_pellet_repeat(npel)=0
      END DO
      
!      NPELMAX=1
!      DO NPEL=1,NPELM
!         MDLPEL(NPEL) = 1
!         PELIN(NPEL) = 0.D0
!         PELR0(NPEL)  = 0.D0
!         PELRW(NPEL)  = 0.5D0
!         PELRAD(NPEL) = 0.D0
!         PELVEL(NPEL) = 0.D0
!         PELTIM(NPEL) = -10.D0

!         PELPAT(1,NPEL) = 1.0D0
!         PELPAT(2,NPEL) = 1.0D0
!         DO NS=3,NSMM
!            PELPAT(NS,NPEL) = 0.0D0
!         ENDDO
!
!         pellet_time_start(NPEL)=0.D0
!         pellet_time_interval(NPEL)=1.D0
!         number_of_pellet_repeat(NPEL)=0
!      END DO

      !  ==== CONTINUOUS PARTICLE SOURCE PARAMETERS ====

      !  NPSCMAX       : Number of particle sources (Max=npscm)
      !  MDLPSC(npscm) : Partcicle SOurce MODEL TYPE
      !              0:OFF  1:GAUSSIAN
      !  NSPSC(npscm)  : Ion species of particle source
      !  PSCIN(npscm) : Time rate of particle source [10^{20} /s]
      !  PSCR0(npscm)  : Radial position of particle source [m]
      !  PSCRW(npscm)  : Radial width of particle source [m]

      NPSCMAX = 1
      DO NPSC=1,NPSCM
         MDLPSC(NPSC) = 0
         PSCIN(NPSC)  = 0.D0
         PSCR0(NPSC)  = 0.D0
         PSCRW(NPSC)  = 0.5D0
         NSPSC(NPSC)  = 2
      END DO

      !  ==== CURRENT DRIVE PARAMETERS ====

      !  MDLCD : CURRENT DRIVE OPERATION MODEL
      !            0: TOTAL PLASMA CURRENT FIXED
      !            1: TOTAL PLASMA CURRENT VARIABLE
      !  PBSCD : BOOTSTRAP CURRENT DRIVE FACTOR

      MDLCD  = 0
      PBSCD  = 1.D0

      !  ==== SAWTOOTH PARAMETERS ====

      !  MDLST  : SAWTOOTH MODEL TYPE
      !              0:OFF
      !              1:ON
      !  IZERO  : SAWTOOTH CRASH TYPE
      !  NGPST  : graphic parameter
      !  TPRST  : SAWTOOTH PERIOD (S)
      !  TSST   : artifical sawtooth time

      MDLST  = 0
      IZERO  = 3
      NGPST  = 4
      TPRST  = 0.1D0
      TSST   = 1.D9

      !  ==== INPUT FROM EXPERIMENTAL DATA ====

      !  MODEP : initial profile selector for steady-state simulation

      MODEP=3

      !  ==== INITIAL CURRENT PROFILE SWITCH ====

      !  MDLJQ :
      !      0 : create QP(NR) profile from experimental CURTOT profile
      !      1 : create AJ(NR) profile from experimental Q profile
      !      2 : no current with experimental Q profile for helical

      MDLJQ=0

      !  ==== Simulateion CONTROL PARAMETERS ====

      !  DT     : SIZE OF TIME STEP
      !  NRMAX  : NUMBER OF RADIAL MESH POINTS
      !  NTMAX  : NUMBER OF TIME STEP
      !  NTSTEP : INTERVAL OF SNAP DATA PRINT
      !  NGTSTP : INTERVAL OF TIME EVOLUTION SAVE
      !  NGRSTP : INTERVAL OF RADIAL PROFILE SAVE

      DT     = 0.01D0
      NRMAX  = 50
      NTMAX  = 100
      NTSTEP = 10
      NGTSTP = 2
      NGRSTP = 100

      !  ==== Convergence Parameter ====

      !  EPSLTR : CONVERGENCE CRITERION OF ITERATION
      !  LMAXTR : MAXIMUM COUNT OF ITERATION

      EPSLTR = 0.001D0
      LMAXTR = 10

      !  ==== Eqs. Selection Parameter ====

      MDLEQB=1  ! 0/1 for B_theta
      MDLEQN=0  ! 0/1 for density
      MDLEQT=1  ! 0/1 for heat
      MDLEQU=0  ! 0/1 for rotation
      MDLEQZ=0  ! 0/1 for impurity
      MDLEQ0=0  ! 0/1 for neutral
      MDLEQE=0  ! 0/1/2 for electron density
      !           0: electron only, 1: both, 2: ion only
      MDLEOI=0  ! 0/1/2 for electron only or bulk ion only if NSMAX=1
      !           0: both, 1: electron, 2: ion

      !  ==== LAPACK parameter ====

      !  MDLPCK : 0    : off (using BANDRD for band matrix solver)
      !           1    : on  (using DGBTRF and DGBTRS for band matrix solver)
      !           else : on  (using DGBSV for band matrix solver)

      MDLPCK=0

      !  ==== FILE NAME ====

      !  KNAMEQ  : file name of equilibrium data
      !  KNAMEQ2 : file name of equilibrium data
      !  KNAMTR  : file name of transport data
      !  KFNLOG  : file name of log data
      !  KFNTXT  : file name of text output
      !  KFNCVS  : file name of cvs output

      KNAMEQ='eqdata'
      KNAMEQ2=''
      KNAMTR='trdata'
      KFNLOG='tr.log'
      KFNTXT='tr.txt'
      KFNCVS='tr.cvs'

      RETURN
      END SUBROUTINE tr_init
END MODULE trinit
