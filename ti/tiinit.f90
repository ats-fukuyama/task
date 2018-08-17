!
MODULE tiinit

  PRIVATE
  PUBLIC ti_init

CONTAINS

!     ***********************************************************

!           INITIALIZE CONSTANTS AND DEFAULT VALUES

!     ***********************************************************

  SUBROUTINE ti_init

    USE ticomm_parm
    IMPLICIT NONE
    INTEGER(ikind) NS

!     ==== DEVICE PARAMETERS ==== * defined in plinit.f90

!        RR*    : PLASMA MAJOR RADIUS (M)
!        RA*    : PLASMA MINOR RADIUS (M)
!        RKAP*  : ELIPTICITY OF POLOIDAL CROSS SECTION
!        RDLT*  : TRIANGULARITY OF POLOIDAL CROSS SECTION
!        BB*    : TOROIDAL MAGNETIC FIELD ON PLASMA AXIS (T)

!        RIP    : PLASMA CURRENT (MA)

      RIP     = 3.D0

!     ==== PLASMA PARAMETERS ====

!        NSMAX* : NUMBER OF MAIN PARTICLE SPECIES (NS=1:ELECTRON)

!        PM(NS)*: Atomic mass number
!        PZ(NS)*: Atomic number (Charge number in full ionized state)
!        PN(NS)*: INITIAL NUMBER DENSITY ON AXIS (1.E20 M**-3)
!        PNS(NS)*:INITIAL NUMBER DENSITY ON SURFACE (1.E20 M**-3)
!        PTPR(NS)*: INITIAL TEMPERATURE ON AXIS (KEV)
!        PTPP(NS)*: INITIAL TEMPERATURE ON AXIS (KEV)
!        PTS(NS)*: INITIAL TEMPERATURE ON SURFACE (KEV)
!        PU(NS)*: INITIAL TOROIDAL VELOCITY ON AXIS (KEV)
!        PUS(NS)* INITIAL TOROIDAL VELOCITY ON SURFACE (KEV)
!        PROFN(NS)* : PROFILE PARAMETER OF INITIAL DENSITY
!        PROFT(NS)* : PROFILE PARAMETER OF INITIAL TEMPERATURE
!        PROFU(NS)* : PROFILE PARAMETER OF INITIAL TOROIDAL VELOCITY

!        KID_NS(NS)*  : Name of particle species
!        ID_NS(NS)*   : Type of particle species
!        NZMIN_NS(NS)*: minimum of particle charge number: 0 for neutral
!        NZMAX_NS(NS)*: maximum of particle charge number: PS for fully-ionized
!        NZINI_NS(NS)*: charge number of initial profile: NZINI=-1 for all 0
      
!        PM(NS),PZ(NS) : definition of element for NS
!        IF(ID_NS==-1) (electron)
!           NSA.range=1
!   　　    PMA(NSA)=PM(NS)
!           PZA(NSA)=PZ(NS)
!        ELSEIF(ID_NS==0) Not solved
!           NSA.range=0
!        ELSE IF(ID_NS==1) (thermal ion)
!           NSA.range=1
!   　　    PMA(NSA)=PM(NS)
!           PZA(NSA)=PZ(NS)
!        ELSE IF(ID_NS==2) (beam ion)
!           NSA.range=1
!   　　    PMA(NSA)=PM(NS)
!           PZA(NSA)=PZ(NS)
!        ELSE IF(ID_NS==5) average-ionization [ADPOST]
!           NSA.range=1
!           PMA(NSA)=PM(NS)
!           PZA(NSA)=PZ(NS)
!        ELSE IF(ID_NS==6) average-ionization [OPEN-ADAS]
!           NSA.range=1
!           PMA(NSA)=PM(NS)
!           PZA(NSA)=PZ(NS)
!        ELSE IF(ID_NS==10) multi-level ionization [OPEN-ADAS]
!           NSA.range=NZMAX_NS(NS)-NZMIN_NS(NS)+1 (max PZ+1)
!           PMA(NSA)=PM(NS)
!           PZA(NSA)=NS-NSAMIN+NZMIN_NS(NS)
!           NZINI_NS: NZ for initial profile, 0 for all zero
!           RNA(NSAMIN:NSAMIN+NSA.range)
!           RTA(NSAMIN)
!           RUA(NSAMIN)
!        ELSE IF(ID_NS==11) multi-level ionization [OPEN-ADAS]
!           RNA(NSAMIN:NSAMIN+NSA.range)
!           RTA(NSAMIN:NSAMIN+NSA.range)
!           RUA(NSAMIN)
!        ELSE IF(ID_NS==12) multi-level ionization [OPEN-ADAS]
!           RNA(NSAMIN:NSAMIN+NSA.range)
!           RTA(NSAMIN:NSAMIN+NSA.range)
!           RUA(NSAMIN:NSAMIN+NSA.range)
!        END IF

!        PT(NS)=(PTPR(NS)+2*PTPP(NS))/3
    NS=1
       KID_NS(NS)= ' e'
       ID_NS(NS)=-1
       NZMIN_NS(NS) = 0
       NZMAX_NS(NS) = 0
       NZINI_NS(NS) = 0
    DO NS =2, NSM
       KID_NS(NS)= ' H'
       ID_NS(NS)= 1
       NZMIN_NS(NS) = 0
       NZMAX_NS(NS) = 1
       NZINI_NS(NS) = 0
    END DO

!     ==== PROFILE PARAMETERS ====

!        PROFJ  : PROFILE PARAMETER OF INITIAL CURRENT DENSITY
!                    (X0-XS)(1-RHO**PROFX1)**PROFX2+XS

      PROFJ1 =-2.D0
      PROFJ2 = 1.D0

!     ==== CONTROL PARAMETERS ====

!        DT     : SIZE OF TIME STEP [s]
!        NRMAX  : NUMBER OF RADIAL MESH POINTS
!        NTMAX  : NUMBER OF TIME STEP
!        NTSTEP : INTERVAL OF SNAP DATA PRINT
!        NGTSTEP: INTERVAL OF TIME EVOLUTION SAVE
!        NGRSTEP: INTERVAL OF RADIAL PROFILE SAVE
!        NGTMAX_INIT: INITIAL SIZE OF NGTSTEP DATA
!        NGRMAX_INIT: INITIAL SIZE OF NGRSTEP DATA


      DT     = 0.01D0
      NRMAX  = 50
      NTMAX  = 100
      NTSTEP = 10
      NGTSTEP = 2
      NGRSTEP = 100
      NGTMAX_INIT = 10000
      NGRMAX_INIT = 1000

!     ==== Convergence Parameter ====

!        EPSLTI : CONVERGENCE CRITERION OF ITERATION
!        LMAXTI : MAXIMUM COUNT OF ITERATION

      EPSLTI = 0.001D0
      LMAXTI = 20

!     ==== Solver Parameter ====
!        MODEL_EQB: Solve poloidal magnetic field (0:Off, 1:On)
!        MODEL_EQN: Solve density (0:Off, 1:On)
!        MODEL_EQT: Solve temperature (0:Off, 1:On)
!        MODEL_EQU: Solve toroidal velocity (0:Off, 1:On)

      MODEL_EQB=0
      MODEL_EQN=1
      MODEL_EQT=0
      MODEL_EQU=0

!     ==== TRANSPORT MODEL ====

!        MODEL_KAI: TURBULENT THERMAL TRANSPORT MODEL
!                     0: fixed
!                    10: CDBM model
!        MODEL_DRR: TURBULENT PARTICLE TRANSPORT MODEL
!                     0: fixed
!                    10: CDBM model
!        MODEL_VR:  TURBULENT PARTICLE PINCH MODEL
!                     0: fixed
!                    10: CDBM model
!        MODEL_NC:  NEOCLASSICAL TRANSPORT MODEL
!                     0: None
!                    10: NCLASS model

      MODEL_KAI= 0
      MODEL_DRR= 0
      MODEL_VR = 0
      MODEL_NC =10

!        MODEL_NF:  NUCULER FUSION REACTION MODEL
!                     0: None
!                    10: On
!        MODEL_NB:  NEUTRAL BEAM INJECTION H/CD MODEL
!                     0: None
!                    10: On
!        MODEL_EC:  ELECTRON CYCLOTRON WAVE H/CD MODEL
!                     0: None
!                    10: On
!        MODEL_LH:  LOWER HYBRID WAVE H/CD MODEL
!                     0: None
!                    10: On
!        MODEL_IC:  ION CYCLOTRON WAVE H/CD MODEL
!                     0: None
!                    10: On
!        MODEL_CD:  NON-INDUCTIVE CURRENT DRIVE MODEL
!                     0: None
!                    10: On
!        MODEL_SYNC: SYNCHROTRON RADIATION MODEL
!                     0: None
!                    10: On
!        MODEL_PEL: PELLET INJECTION MODEL
!                     0: None
!                    10: On
!        MODEL_PSC: PARTICLE SOURCE MODEL
!                     0: None
!                    10: On
!        MODEL_PROF*: PROFILE MODEL
!                     0: PROFX1/2(NS)=PROFX1/2(1)
!                     1: PROFX1/2(NS): separately defined

      MODEL_NF  =0
      MODEL_NB  =0
      MODEL_EC  =0
      MODEL_LH  =0
      MODEL_IC  =0
      MODEL_CD  =0
      MODEL_SYNC=0
      MODEL_PEL =0
      MODEL_PSC =0

!     ==== TRANSPORT PARAMETERS ====

!        AK0    : THERMAL DIFFUSION FACTOR
!        AD0    : PARTICLE DIFFUSION FACTOR
!        AV0    : INWARD PARTICLE PINCH FACTOR
!        DK0    : CENTERAL VALUE of KAI-fixed
!        DKS    : SURFACE VALUE of KAI-fixed

      AK0    = 1.0D0
      AD0    = 0.1D0
      AV0    = 0.5D0
      DK0    = 0.3D0
      DKS    = 3.0D0

!     ==== Source Parameters ====

!     PNBIN: NBI heating power [MW}
!     PNBR0:  NBI deposition center [m]
!     PNBRW:  NBI deposition width [m]
!     PNBENG: NBI beam energy [MW]
!     PNBRTG: NBI tangentical radius [m]

!     PECIN: EC heating power [MW}
!     PECR0:  EC deposition center [m]
!     PECRW:  EC deposition width [m]
!     PECTOE: EC power partition to electron
!     PECNPR: EC parallel refractive index

!     PLHIN: LH heating power [MW}
!     PLHR0:  LH deposition center [m]
!     PLHRW:  LH deposition width [m]
!     PLHTOE: LH power partition to electron
!     PLHNPR: LH parallel refractive index

!     PICIN: IC heating power [MW}
!     PICR0:  IC deposition center [m]
!     PICRW:  IC deposition width [m]
!     PICTOE: IC power partition to electron
!     PICNPR: IC parallel refractive index

!     PNBCD: NB current drive (0: off, 1: on)
!     PECCD: EC current drive (0: off, 1: on)
!     PLHCD: LH current drive (0: off, 1: on)
!     PICCD: IC current drive (0: off, 1: on)
!     PBSCD: Bootstrap current drive (0: off, 1: on)

!     PELIN: Pellet total number of partilces in a pellet [10^20]
!     PELR0:  Pellet deposition center [m]
!     PELRW:  Pellet deposition width [m]
!     PELPAT(NSM): Pellet partition of particle species
!     PRLTS:  Pellet injection start time [s]
!     PRLDE:  Pellet injection time interval [s]
!     PRLTE:  Pellet injection end time [s]
!     PELRAD: Pellet radius [m]
!     PELVEL: Pellet velocity [m/s]

!     PSCIN: Particle source injection rate [10^20/s]
!     PSCR0:  Particle source deposition center [m]
!     PSCRW:  Particle source deposition width [m]

!     SYNC_WALL: fraction of synchrotron radiation absorption by walls
!     SYNC_CONV: fraction of x (o) mode reflected as x (o) mode

      PNBIN = 0.D0
      PNBR0  = 0.D0
      PNBRW  = 0.5D0
      PNBENG = 80.D0
      PNBRTG = 3.D0
      PNBCD  = 1.D0

      PECIN = 0.D0
      PECR0  = 0.D0
      PECRW  = 0.2D0
      PECTOE = 1.D0
      PECNPR = 0.D0
      PECCD  = 0.D0

      PLHIN = 0.D0
      PLHR0  = 0.D0
      PLHRW  = 0.2D0
      PLHTOE = 1.D0
      PLHNPR = 2.D0
      PLHCD  = 1.D0

      PICIN = 0.D0
      PICR0  = 0.D0
      PICRW  = 0.5D0
      PICTOE = 0.5D0
      PICNPR = 2.D0
      PICCD  = 0.D0

      PBSCD  = 1.D0

      PELIN = 0.D0
      PELR0  = 0.D0
      PELRW  = 0.5D0
      PELPAT(1)=1.D0
      PELPAT(2)=1.D0
      DO NS=3,NSM
         PELPAT(NS) = 0.0D0
      ENDDO
      PELRAD = 0.D0
      PELVEL = 0.D0
      PELTS  = 0.D0
      PELDT  = 0.D0
      PELTE  = 0.D0
      PELRAD = 0.D0
      PELVEL = 0.D0

      PSCIN = 0.D0
      PSCR0  = 0.D0
      PSCRW  = 0.5D0
      PSCPAT(1)=1.D0
      PSCPAT(2)=1.D0
      DO NS=3,NSM
         PSCPAT(NS) = 0.0D0
      ENDDO

      SYNC_WALL=0.2D0
      SYNC_CONV=0.95D0

!     ==== FILE NAME ====

!        KNAMEQ: Filename of equilibrium data*
!        KNAMTR: Filename of transport data*
!        KNAMLOG: LOG FILE NAME

      KNAMLOG='ti.log'

      DO NS=1,NSM
         PT(NS)=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      END DO

    RETURN
  END SUBROUTINE ti_init
END MODULE tiinit
