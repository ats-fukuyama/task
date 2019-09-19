MODULE w1init

CONTAINS

! ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE w1_init

    USE w1comm_parm
    IMPLICIT NONE
    INTEGER:: ns, na

!     ******* INPUT PARAMETERS THROUGH THE NAMELIST /W1/ *******
!
!     ======( MACHINE PARAMETERS )======

!     BB    : MAGNETIC FIELD AT X=0   (T)
!     RR    : PLASMA MAJOR RADIUS     (M)
!     RA    : PLASMA MINOR RADIUS     (M)
!     RD    : ANTENNA RADIUS          (M)
!     RB    : WALL RADIUS             (M)
!     RZ    : PERIODIC LENGTH         (M)  (IN Z DIRECTION)
!                           *** IF RZ.EQ.0. THEN RZ=2*PI*(RR+RB) ***
!     WALLR : WALL RESISTIVITY        (OHM-M)
!     EPSH  : HELICAL RIPPLE IN HELICAL SYSTEM

    BB    = 3.5D0
    RR    = 3.000D0
    RA    = 1.2D0
    RD    = 1.25D0
    RB    = 1.3D0
    RZ    = 2.D0*PI*(RR+RB)
    WALLR = 0.0D-5
    EPSH  = -0.3D0

!     ======( ANTENNA PARAMETERS )======

!     RF    : WAVE FREQUENCY          (MHZ)
!     RKZ   : WAVE NUMBER IN Z-DIRECTION (/M) (VALID FOR NZP=1)
!     NAMAX : NUMBER OF ANTENNAS
!     AJYH  : ANTENNA CURRENT Y AT -RD  (KA)
!     AJYL  :                 Y AT  RD  (KA)
!     AJZH  : ANTENNA CURRENT DENSITY Z AT -RD  (KA/m)
!     AJZL  :                         Z AT  RD  (KA/m)
!     APYH  : TOROIDAL ANGLE CURRENT Y AT -RD ( DEGREE )  -180..180
!     APYL  :                        Y AT  RD ( DEGREE )  -180..180
!     APZH  : TOROIDAL ANGLE CURRENT Z AT -RD ( DEGREE )  -180..180
!     APZL  :                        Z AT  RD ( DEGREE )  -180..180
!     ALZH  : TOROIDAL LENGTH CURRENT Z AT -RD ( DEGREE ) APZH..APZH+ALZH
!     ALZL  :                         Z AT  RD ( DEGREE ) APZL..APZL+ALZL
!     APHH  : PHASE OF CURRENT AT -RD ( DEGREE )
!     APHL  :                  AT  RD ( DEGREE )

    RF    = 52.5D0
    RKZ   = 8.00D0
    NAMAX =   1
    DO NA=1,NAM
       AJYH( NA )  =   0.D0
       AJZH( NA )  =   0.D0
       APYH( NA )  =   0.D0
       APZH( NA )  =   0.D0
       ALZH( NA )  =   0.D0
       APHH( NA )  =   0.D0
       AJYL( NA )  =   0.D0
       AJZL( NA )  =   0.D0
       APYL( NA )  =   0.D0
       APZL( NA )  =   0.D0
       ALZL( NA )  =   0.D0
       APHL( NA )  =   0.D0
    END DO
    AJYL( 1 )  =   1.D0

!     ======( PLASMA PARAMETERS )======

!     NSMAX : NUMBER OF PARTICLE SPECIES
!     PA    : ATOMIC NUMBER
!     PZ    : CHARGE NUMBER
!     PN    : DENSITY AT X=0             (1.E20 M**-3)
!     PTPP  : PERPENDICULAR TEMPERATURE AT X=0   (KEV)
!     PTPR  : PARALLEL TEMPERATURE      AT X=0   (KEV)
!     PU    : PARALLEL DRIFT VELOCITY   AT X=0   (KEV)
!     PNS   : DENSITY    ON SURFACE X=RA (1.E20 M**-3)
!     PTS   : TEPERATURE ON SURFACE X=RA         (KEV)
!     PU    : PARALLEL DRIFT VELOCITY   AT X=0   (KEV)
!     PZCL  : NUMERICAL FACTOR OF COLLISION  (NOT USED IN THIS VERSION)
!     IHARM : MAXIMUM HARMONIC NUMBER (-IHARM TO IHARM)
!               (APPLICABLE FOR NMODEL GE 2.
!                WHEN NMODEL EQ 4 OR NMODEL EQ 5,
!                   IF IHARM GE 3 THEN FAST WAVE FLR FOR 3..IHARM
!                   IF IHARM LT 0 THEN FAST WAVE FLR FOR 0..ABS(IHARM))
!     IELEC : 1 for electron, 0 for else other
!     IELEC : 1 for electron, 0 for else other

    NSMAX = 3

    PA  ( 1 ) =  5.4466D-4
    PZ  ( 1 ) = -1.D0
    PN  ( 1 ) =  1.000D0
    PTPP( 1 ) =  2.000D0
    PTPR( 1 ) =  2.000D0
    PNS ( 1 ) =  0.1D0
    PTS ( 1 ) =  0.1D0
    PU  ( 1 ) =  0.D0
    PZCL( 1 ) =  0.D0
    IHARM(1 ) =  2
    IELEC(1 ) =  1

    PA  ( 2 ) =  2.D0
    PZ  ( 2 ) =  1.D0
    PN  ( 2 ) =  0.950D0
    PTPP( 2 ) =  2.000D0
    PTPR( 2 ) =  2.000D0
    PNS ( 2 ) =  0.095D0
    PTS ( 2 ) =  0.1D0
    PU  ( 2 ) =  0.D0
    PZCL( 2 ) =  0.D0
    IHARM(2 ) =  2
    IELEC(2 ) =  0

    PA  ( 3 ) =  1.D0
    PZ  ( 3 ) =  1.D0
    PN  ( 3 ) =  0.050D0
    PTPP( 3 ) =  2.000D0
    PTPR( 3 ) =  2.000D0
    PNS ( 3 ) =  0.005D0
    PTS ( 3 ) =  0.1D0
    PU  ( 3 ) =  0.D0
    PZCL( 3 ) =  0.D0
    IHARM(3 ) =  2
    IELEC(3 ) =  0

    DO NS=4,NSM
       PA  ( NS ) =  1.D0
       PZ  ( NS ) =  1.D0
       PN  ( NS ) =  0.D0
       PTPP( NS ) =  1.D0
       PTPR( NS ) =  1.D0
       PNS ( NS ) =  0.D0
       PTS ( NS ) =  0.1D0
       PU  ( NS ) =  0.D0
       PZCL( NS ) =  0.D0
       IHARM(NS ) =  2
       IELEC(NS ) =  0
    END DO

!   ======( PROFILE CONTROL PARAMETERS )======

!     APRFPN: PROFILE FACTOR FOR PLASMA DENSITY
!     APRFTR: PROFILE FACTOR FOR PARALLEL TEMPERATURE
!     APRFTP: PROFILE FACTOR FOR PERPENDICULAR TEMPERATURE

    APRFPN= 0.5D0
    APRFTR= 1.D0
    APRFTP= 1.D0

!   ======( PROGRAM CONTROL )======
!
!     NXMAX: NUMBER OF X-MESHES
!            NXPMAX = NXMAX*RA/RB : NUMBER OF X-MESHES IN PLASMA
!            NXVMAX = NXMAX-NXPMAX: NUMBER OF X-MESHES IN VACUUM (EVEN NUMBER)
!     NZMAX: NUMBER OF Z-MESHES  (POWER OF 2, 2**N), 1 for single mode
!
!     NMODEL: 0 : NO FLR CONDUCTIVITY MODEL        (FINITE ELEMENT)
!             1X: NO FLR CONDUCTIVITY MODEL        (MULTI LAYER)
!             2 : FAST WAVE FLR CONDUCTIVITY MODEL (FINITE ELEMENT)
!             3X: FAST WAVE FLR CONDUCTIVITY MODEL (MULTIPLE LAYER)
!             4 : DIFFERENTIAL CONDUCTIVITY MODEL  (FINITE ELEMENT)
!             5X: DIFFERENTIAL CONDUCTIVITY MODEL  (MULTIPLE LAYER)
!             6 : COLD CONDUCTIVITY MODEL          (FULL FINITE ELEMENT) REFL
!             7 : WARM CONDUCTIVITY MODEL          (FULL FINITE ELEMENT) REFL
!             8 : INTEGRAL CONDUCTIVITY MODEL      (FULL FINITE ELEMENT) REFL
!             9 : COLD CONDUCTIVITY MODEL          (FULL FINITE ELEMENT) ABS
!            10 : WARM CONDUCTIVITY MODEL          (FULL FINITE ELEMENT) ABS
!            11 : INTEGRAL CONDUCTIVITY MODEL      (FULL FINITE ELEMENT) ABS
!            12 : INTEGRAL CONDUCTIVITY MODEL      (FINITE ELEMENT)
!
!     NPRINT: 0 :    LP OUTPUT (GLOBAL DATA)
!             1 :    LP OUTPUT (1-D DATA K-DEPENDENCE)
!             2 :    LP OUTPUT (1-D DATA X-DEPENDENCE)
!             3 :    LP OUTPUT (2-D DATA ELECTRIC FIELD)
!
!     NFILE : 0 : NO FILE OUTPUT
!             1 :    FILE OUTPUT (1-D DATA)
!             2 :    FILE OUTPUT (2-D DATA)
!
!     NGRAPH: 0 : NO GRAPHICS
!             1 :    GRAPHIC (E,J,P: 5 FIGS)
!             2 :    GRAPHIC (E,P: 4 FIGS)
!             3 :    GRAPHIC (E,P/TWICE HEIGHT: 4 FIGS)
!            +4 :    GRAPHIC (B FIELD, VECTOR POTENTIAL: 6 FIGS)
!            +8 :    GRAPHIC (HELICITY AND FORCE: 4 FIGS)
!           +16 :    GRAPHIC (HELICITY AND CURRENT: 4 FIGS)
!           +32 :    GRAPHIC (E,P,J: 4+1 FIGS)
!           +64 :    GRAPHIC (HELICITY, FORCE AND CURRENT: 8FIGS)
!
!     NLOOP : NUMBER OF PARAMETER SCAN LOOP
!     DRF  = INCREMENT OF RF
!     DRKZ = INCREMENT OF RKZ
!
!     NSYM  : 0 : WITHOUT SYMMETRY IN Z-DIRECTION
!             1 : WITH    SYMMETRY
!            -1 : WITH   ASYMMETRY
!
!     MODELN:  0: parabolic with pn=0 in SOL
!              1: parabolic with pn=pns in SOL
!             10: linear with pn=0 in SOL
!             11: linear with pn=pns in SOL
!
!     NALPHA: 0 : NO ALPHA PARTICLE EFFECT
!             1 : SLOWING DOWN DISTRIBUTION FOR ALPHA (IS=4, NMODEL=1)
!             2 : ALPHA DENSITY BY FUSION REACTION (IS 2:D 3:T 4:ALPHA)
!             3 : ALPHA DENSITY AND SLOWING DOWN DISTRIBUTION
!
!     NSYS  : 0 : TOKAMAK CONFIGURATION
!             1 : HELICAL CONFIGURATION
!
!     NDISP : 0 : FIELD CALCULATION
!             1 : WAVE NUMBER DISPLAY
!
!     XDMAX : MAXIMUM VALUE OF GYRORADIUS IN INTEGRATION (NMODEL=6)
!     NDMAX : MAXIMUM NUMBER OF INTEGRAL TABLE (NMODEL=6)
!     DXFACT: MESH ACCUMULATION FACTOR (DEFAULT 0.D0)
!     DXWDTH:                   WIDTH  (DEFAULT 3.D0 : 3*LARMOR RADIUS)
!     NXABS : 0 : NO ABSORBING WALL
!            .GT.0 :    ABSORBING WALL AT NX=NXABS
!            .LT.0 :    ABSORBING WALL AT NX=ABS(NXABS) FOR BACKWARD WAVE

    NXMAX= 130
    NZMAX= 1

    NMODEL= 5
    MODELN= 0

    NPRINT= 0
    NFILE = 0
    NGRAPH= 1
    NLOOP = 1
    DRF   = 0.D0
    DRKZ  = 0.D0
    NSYM  = 0
    NALPHA= 0
    NSYS  = 0
    NDISP = 0
    XDMAX = 10.D0
    NDMAX = 100
    DXFACT= 0.D0
    DXWDTH= 4.D0
    NXABS = 0

!     ======( CURRENT DRIVE PARAMETER ) =====

!     NCDTYP: CURRENT DRIVE TYPE (0:LANDAU (LH,EC), 1:TTMP (IC))
!     ZEFF  : Z EFFECTIVE FOR CURRENT DRIVE EFFICIENCY
!     WVYSIZ: Vertical size of wave guide

    NCDTYP=1
    ZEFF=2.D0
    WVYSIZ=0.D0

!     ======( WG Parameters )======
!     MDLWG: Wave guide excitation parameter
!            0  no WG
!            1  LFS O-mode
!            2  LFS X-mode
!            3  HFS O-mode
!            4  HFS X-mode
!     MDLWGS: Wave guide shape parameter
!            0: Rectangular shape
!            1: Parabolic shape
!            2: Gaussian shape
!     WGZ1   : Lower end position of WG [m]
!     WGZ2   : Upper end position of WG [m]
!     WGAMP  : Waveguide electric field amplitude [V/m]
!     WGNZ   : Parallel refractive index

      MDLWG  = 0
      MDLWGS = 0
      WGZ1   = 0.05D0
      WGZ2   = 0.15D0
      WGAMP  = 1.D0
      WGNZ   = 0.2D0

    RETURN
  END SUBROUTINE w1_init
END MODULE w1init
