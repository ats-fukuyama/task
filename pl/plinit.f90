  MODULE plinit

  CONTAINS

!     ****** INITIALIZE INPUT PARAMETERS ******

    SUBROUTINE pl_init

      USE plcomm
      IMPLICIT NONE
      INTEGER:: NS

      CALL pl_allocate_ns

!     ======( DEVICE PARAMETERS )======

!        RR    : Plasma major radius                             (m)
!        RA    : Plasma minor radius                             (m)
!        RB    : Wall minor radius                               (m)
!        RKAP  : Plasma shape elongation
!        RDEL  : Plasma shape triangularity *
!        BB    : Magnetic field at center                        (T)
!        Q0    : Safety factor at center
!        QA    : Safety factor on plasma surface
!        RIP   : Plasma current                                 (MA)
!        PROFJ : Curren density profile parameter (power of (1 - rho^2))

      RR    = 3.D0
      RA    = 1.D0
      RB    = 1.2D0
      RKAP  = 1.D0
      RDLT  = 0.D0

      BB    = 3.D0
      Q0    = 1.D0
      QA    = 3.D0
      RIP   = 3.D0
      PROFJ = 2.D0

!     *** CONFIGURATION PARAMETERS (MIRROR: MODELG=0) ***

!        RMIR  : Mirror ratio
!        ZBB   : Periodic length along magnetic axis           (m)

      RMIR   = 2.0D0
      ZBB    = 0.15D0

!     *** CIRCULAR COIL PARAMETERS (MODELG=0 MODELB=1) ***

!        NCMAX : Number of coil
!        RC    : Radial position of coil current               (m)
!        ZC    : Axial position of coil current                (m)
!        BC    : Magnetic field on axis, center of coil        (T)

      NCOILMAX  = 3
      RCOIL(1)  = 0.35D0
      ZCOIL(1)  = 0.D0
      BCOIL(1)  = 0.001D0
      RCOIL(2)  = 0.35D0
      ZCOIL(2)  = 0.05D0
      BCOIL(2)  =-0.001D0
      RCOIL(3)  = 0.35D0
      ZCOIL(3)  =-0.05D0
      BCOIL(3)  =-0.001D0
      

!     *** CONFIGURATION PARAMETERS (HELICAL: MODELG=11) ***

!        Hpitch1: Helical pitch (2*pi/L) for B profile
!        Hpitch2: Helical pitch (2*pi/L) for metric
!        RRCH  : Helical Coil radius                                   (m)

      Hpitch1= 1.25D0
      Hpitch2= 1.25D0
      RRCH   = 0.95D0

!     ======( PLASMA PARAMETERS )======

!        NSMAX : Number of particle species
!        NPA   : Atomic number (0 for electron)
!        PA    : Mass number (to be replaced by PM)
!        PZ    : Charge number
!        PN    : Density at center                     (1.0E20/m**3)
!        PNS   : Density on plasma surface             (1.0E20/m**3)
!        PTPR  : Parallel temperature at center                (keV)
!        PTPP  : Perpendicular temperature at center           (keV)
!        PTS   : Temperature on surface                        (keV)
!        PU    : Toroidal rotation velocity at center          (m/s)
!        PUS   : Toroidal rotation velocity on surface         (m/s)
!        PUPR  : typical parallel velocity                     (m/s)
!        PUPP  : typical perpendicular velocity                (m/s)
!        RHOITB: rho at ITB (0 for no ITB)
!        PNITB : Density increment at ITB              (1.0E20/Mm*3)
!        PTITB : Temperature increment at ITB                  (keV)
!        PUITB : Toroidal rotation velocity increment at ITB   (m/s)
!        PNUC  : Factor for collision frequency
!                   model_coll=1:  RNUC (no nueutral)
!                   model_coll!=1: RNUC (with nueutral)
!        PZCL  : Factor for collisional damping
!                   model_coll=0: nu/omega=PZCL
!                   model_coll=1: nu/omega=RNUC (no nueutral)
!                   model_coll=2: nu/omega=RNUC (with nueutral)
!        KID_NS: index of particle species
!        ID_NS : -1 : electron
!                 0 : neutral
!                 1 : ion
!                 2 : fast ion

      NSMAX = 2                  ! Default number of particle species

!     *** electron ***
         NS = 1

         KID_NS(NS)= ' e'
         ID_NS(NS) = -1
         NPA(NS)  = 0
         PA(NS)   = AME/AMP
         PZ(NS)   =-1.0D0
         PN(NS)   = 1.0D0
         PNS(NS)  = 0.0D0
         PTPR(NS) = 5.0D0
         PTPP(NS) = 5.0D0
         PTS(NS)  = 0.05D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         PUPR(NS) = 0.D0
         PUPP(NS) = 0.D0
         RHOITB(NS)=0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0
         PNUC(NS) = 0.D0
         PZCL(NS) = 0.D0

!     *** deuteron ***
         NS = 2

         KID_NS(NS)= ' D'
         ID_NS(NS) = 1
         NPA(NS)  = 1
         PA(NS)   = 2.0D0
         PZ(NS)   = 1.0D0
         PN(NS)   = 1.0D0
         PNS(NS)  = 0.0D0
         PTPR(NS) = 5.0D0
         PTPP(NS) = 5.0D0
         PTS(NS)  = 0.05D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         PUPR(NS) = 0.D0
         PUPP(NS) = 0.D0
         RHOITB(NS)=0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0
         PNUC(NS) = 0.D0
         PZCL(NS) = 0.D0

!     *** triton ***

         NS = 3

         KID_NS(NS)= ' T'
         ID_NS(NS) = 1
         NPA(NS)  = 1
         PA(NS)   = 3.0D0
         PZ(NS)   = 1.0D0
         PN(NS)   = 1.0D0
         PNS(NS)  = 0.0D0
         PTPR(NS) = 5.0D0
         PTPP(NS) = 5.0D0
         PTS(NS)  = 0.05D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         PUPR(NS) = 0.D0
         PUPP(NS) = 0.D0
         RHOITB(NS)=0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0
         PNUC(NS) = 0.D0
         PZCL(NS) = 0.D0

!     *** Helium ion ***
         NS = 4

         KID_NS(NS)= 'He'
         ID_NS(NS) = 1
         NPA(NS)  = 2
         PA(NS)   = 4.0D0
         PZ(NS)   = 2.0D0
         PN(NS)   = 1.0D0
         PNS(NS)  = 0.0D0
         PTPR(NS) = 5.0D0
         PTPP(NS) = 5.0D0
         PTS(NS)  = 0.05D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         PUPR(NS) = 0.D0
         PUPP(NS) = 0.D0
         RHOITB(NS)=0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0
         PNUC(NS) = 0.D0
         PZCL(NS) = 0.D0

         ! *** dummy ***
      DO NS = 5, NSM
         KID_NS(NS)= ' H'
         ID_NS(NS)= 1
         NPA(NS)  = 1
         PA(NS)   = 1.0D0
         PZ(NS)   = 1.0D0
         PN(NS)   = 0.0D0
         PNS(NS)  = 0.0D0
         PTPR(NS) = 5.0D0
         PTPP(NS) = 5.0D0
         PTS(NS)  = 0.0D0
         PU(NS)   = 0.D0
         PUS(NS)  = 0.D0
         PUPR(NS) = 0.D0
         PUPP(NS) = 0.D0
         RHOITB(NS)=0.D0
         PNITB(NS)= 0.D0
         PTITB(NS)= 0.D0
         PUITB(NS)= 0.D0
         PNUC(NS) = 0.D0
         PZCL(NS) = 0.D0
      ENDDO

!     *** PLANE  PARAMETERS ***

  r_corner(1)=0.D0
  r_corner(2)=RA
  r_corner(3)=0.D0
  z_corner(1)=0.D0
  z_corner(2)=0.D0
  z_corner(3)=2.D0*PI*RR

  br_corner(1:3)=0.D0
  bz_corner(1:3)=0.D0
  bt_corner(1:3)=BB

  DO ns=1,nsm
     pn_corner(1,ns)=pn(ns)
     pn_corner(2,ns)=pns(ns)
     pn_corner(3,ns)=pn(ns)
     ptpr_corner(1,ns)=ptpr(ns)
     ptpr_corner(2,ns)=pts(ns)
     ptpr_corner(3,ns)=ptpr(ns)
     ptpp_corner(1,ns)=ptpp(ns)
     ptpp_corner(2,ns)=pts(ns)
     ptpp_corner(3,ns)=ptpp(ns)
  END DO

!     ======( PROFILE PARAMETERS )======

!        PROFN1: Density profile parameter (power of rho)
!        PROFN2: Density profile parameter (power of (1 - rho^PROFN1))
!        PROFT1: Temperature profile parameter (power of rho)
!        PROFT2: Temperature profile parameter (power of (1 - rho^PROFN1))
!        PROFU1: Rotation profile parameter (power of rho)
!        PROFU2: Rotation profile parameter (power of (1 - rho^PROFN1))

  DO NS=1,NSM
     PROFN1(NS)= 2.D0
     PROFN2(NS)= 0.5D0
     PROFT1(NS)= 2.D0
     PROFT2(NS)= 1.D0
     PROFU1(NS)= 2.D0
     PROFU2(NS)= 1.D0
  END DO

!     ======( TRAVIS PROFILE PARAMETERS )======

!        PROF=g-h+(1-g+h)(1-x^p)^q + h(1-EXP(-x^2/w^2))

  profn_travis_g=0.D0
  profn_travis_h=0.D0
  profn_travis_p=2.D0
  profn_travis_q=1.D0
  profn_travis_w=1.D0
  proft_travis_g=0.D0
  proft_travis_h=0.D0
  proft_travis_p=2.D0
  proft_travis_q=1.D0
  proft_travis_w=1.D0

  !     ======( MODEL PARAMETERS )======

!        MODELG: Control plasma geometry model
!              0: XYZ Slab geometry
!                 MODELB: 
!                      0: Translational mirror
!                      1: Translational mirror with curent rod (NCMAX.GT.0)
!                      2: XY 2D plane profile (defined at corners)
!                               (MODELGX: 0:linear, 1,2:parabolic)
!                      3: XY 2D plane Read 2D mag file 
!              1: RZphi Cylindrical geometry
!                 MODEL B:
!                      0: Axisymmetric mirror
!                      1: Axisymmetric mirror with circular coils (NCMAX.GT.0)
!                      2: RZ 2D toroidal profile (defined at corners)
!                          (MODELGX: 0:linear,1,2:parabolic)
!                      3: RZ 2D toroidal Read 2D mag file 
!              2: RZphi Toroidal geometry
!              3: RZphi Read TASK/EQ output geometry
!              4: RZphi Read VMEC output geometry
!              5: RZphi Read EQDSK output geometry
!              6: RZphi Read Boozer output geometry
!              7: RZphi Read new VMEC output geometry
!              8: RZphi call TASK/EQU
!              9: RZphi call TASK/EQ
!             10: reserved for GAMMA-10
!             11: Straight helical geometry
!             12: 2D plane profile (B read from file)
!             13: 2D plane profile (simple parabolic cylinder)

!        MODELN: Control plasma profile
!                   0: Calculated from PN,PNS,PTPR,PTPP,PTS,PU,PUS; PN=0 in SOL
!                   1: Calculated from PN,PNS,PTPR,PTPP,PTS,PU,PUS; PNS in SOL
!                   2: n,T from pressure profile; u from PU,PUS; PNS in SOL
!                   3: with RHOEDG from PN,PNS,PTPR,PTPP,PTS,PU,PUS; PNS in SOL
!                   8: Read from file through WMXPRF (JT-60)
!                   9: Read from bpsd_plasmaf
!                  21: Read from trdata
!                  31: Calculated from profn_travis and proft_travis
!        MODELQ: Control safety factor profile (for MODELG=1,2)
!                   0: Parabolic q profile (Q0,QA,RHOMIN)
!                   1: Given current profile (RIP,PROFJ)
!        model_coll: collision parameter
!                   0: nu_coll=PZCL*omega
!                   1: nu_coll=Coulomb collision freq.
!                   2: nu_coll=Coulomb collision + neutral collision freq.
!        MODEL_PROF: profile parameter
!                   0: PROFX1(NS)=PROFX1(1),PROFX2(NS)=PROFX2(1): compatibility
!                   1: PROFX1(NS),PROFX2(NS): defined separately
!        MODEL_NPROF: neutral profile parameter
!                   0: Flat profile
!                   1: Flat only in plasma, 0 outside
!                   2: (1-psi) dependence, 0 outside

      MODELG= 2
      MODELB= 0
      MODELN= 0
      MODELQ= 0
      model_coll=0
      MODEL_PROF=0
      MODEL_NPROF=0

!        RHOMIN: rho at minimum q (0 for positive shear)
!        QMIN  : q minimum for reversed shear
!        RHOITB: rho at ITB (0 for no ITB)
!        RHOEDG: rho at EDGE for smoothing (1 for no smooth)

      RHOMIN = 0.D0
      QMIN   = 1.5D0
      RHOEDG = 1.D0

!        BAXIS_SCALED: if not 0.D0, magnetic field is scaled
!                      by Baxis_scaled/Baxis_real

      BAXIS_SCALED=0.D0
      
!        PPN0: Neutral pressure [Pa] 1 Torr = 1 mmHg = 133.322 Pa
!        PTN0: Neutral temperature [eV]
!        RF_PL: wave frequency [MHz], usually set in wave code

      PPN0 = 3.0D0
      PTN0 = 0.03D0
      RF_PL = 1.D6    ! 1THz to reduce pzcl

!     ======( GRAPHIC PARAMETERS )======

!        RHOGMN: minimum rho in radial profile
!        RHOGMX: maximum rho in radial profile

      RHOGMN = 0.D0
      RHOGMX = 1.D0

!     ======( profile for wave parameter )======

!        mdlplw: 0  rupr=pupr, rupp=pupp
!                1  rupr=ru*bnt+rupl*bnp  rupp=-ru*bnp+rupl*bnt
!                        ru:   toroidal velocity, rupl: poloidal velocity 
!                        rupr: parallel velocity, rupp: perpendicular velocity 
      mdlplw=0


!     ======( IO FILE NAMES )======

!        KNAMEQ: Filename of equilibrium data
!        KNAMWR: Filename of ray tracing data
!        KNAMWM: Filename of full wave data
!        KNAMFP: Filename of Fokker-Planck data
!        KNAMFO: Filename of File output
!        KNAMPF: Filename of profile data
!        KNAMTR: Filename of transport data

      KNAMEQ = 'eqdata'
      KNAMWR = 'wrdata'
      KNAMWM = 'wmdata'
      KNAMFP = 'fpdata'
      KNAMFO = 'fodata'
      KNAMPF = 'pfdata'
      KNAMTR = 'trdata'

!     ======( FILE IO MODES )======

!        MODEFR: File name interaction at reading file
!                 0 : WITHOUT PROMPT
!                 1 : WITH FILE NAME INPUT
!        MODEFW: File name interaction at writing file
!                 0 : WITHOUT PROMPT, ALWAYS OVERWRITE
!                 1 : WITHOUT PROMPT, CONFIRM, IF FILE EXISTS
!                 2 : WITHOUT PROMPT, ASK NEW NAME, IF FILE EXISTS
!                 3 : WITHOUT PROMPT, ERROR, IF FILE EXISTS
!                 4 : WITH FILE NAME INPUT, ALWAYS OVERWRITE
!                 5 : WITH FILE NAME INPUT, CONFIRM, IF FILE EXISTS
!                 6 : WITH FILE NAME INPUT, ASK NEW NAME, IF FILE EXISTS
!                 7 : WITH FILE NAME INPUT, ERROR, IF FILE EXISTS

!     ======( Default values )======
!             eq_init sets MODEFR=0 and MODEFW=0 due to histrical reason
      
      MODEFR = 1
      MODEFW = 5

      IDEBUG = 0

      RETURN
    END SUBROUTINE pl_init

END MODULE plinit
