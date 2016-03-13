MODULE w1init

CONTAINS

! ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE w1_init

!     ******* INPUT PARAMETERS THROUGH THE NAMELIST /W1/ *******
!
!     BB    : MAGNETIC FIELD AT X=0   (T)
!     RR    : PLASMA MAJOR RADIUS     (M)
!     RZ    : PERIODIC LENGTH         (M)  (IN Z DIRECTION)
!                           *** IF RZ.EQ.0. THEN RZ=2*PI*RR ***
!     RA    : PLASMA MINOR RADIUS     (M)
!     RD    : ANTENNA RADIUS          (M)
!     RB    : WALL RADIUS             (M)
!     EPSH  : HELICAL RIPPLE IN HELICAL SYSTEM
!
!     WALLR : WALL RESISTIVITY        (OHM-M)
!     NXABS : 0 : NO ABSORBING WALL
!         .GT.0 :    ABSORBING WALL AT NX=NXABS
!         .LT.0 :    ABSORBING WALL AT NX=ABS(NXABS) FOR BACKWARD WAVE
!
!     NCDTYP: CURRENT DRIVE TYPE (0:LANDAU, 1:TTMP)
!     ZEFF  : Z EFFECTIVE FOR CURRENT DRIVE EFFICIENCY
!
!     RF    : WAVE FREQUENCY          (MHZ)
!     RKZ   : WAVE NUMBER IN Z-DIRECTION (/M) (VALID FOR NZP=1)
!     IAMAX : NUMBER OF ANTENNAS
!     AJYH  : ANTENNA     CURRENT Y AT -RD  (KA)
!     AJYL  :                     Y AT  RD  (KA)
!     ALYH  : POSITION OF CURRENT Y AT -RD ( DEGREE )  -180..180
!     ALYL  :                     Y AT  RD ( DEGREE )  -180..180
!     APYH  : PHASE    OF CURRENT Y AT -RD ( DEGREE )
!     APYL  :                     Y AT  RD ( DEGREE )
!     AJZH  : ANTENNA     CURRENT Z AT -RD  (KA)
!     AJZL  :                     Z AT  RD  (KA)
!     DXFACT: MESH ACCUMULATION FACTOR (DEFAULT 0.D0)
!     DXWDTH:                   WIDTH  (DEFAULT 3.D0 : 3*LARMOR RADIUS)
!
!     ISMAX : NUMBER OF PARTICLE SPECIES
!     PA    : ATOMIC NUMBER
!     PZ    : CHARGE NUMBER
!     PN    : DENSITY AT X=0             (1.E20 M**-3)
!     PTPP  : PERPENDICULAR TEMPERATURE AT X=0   (KEV)
!     PTPR  : PARALLEL TEMPERATURE      AT X=0   (KEV)
!     PU    : PARALLEL DRIFT VELOCITY   AT X=0   (KEV)
!     PNS   : DENSITY    ON SURFACE X=RA (1.E20 M**-3)
!     PTS   : TEPERATURE ON SURFACE X=RA         (KEV)
!     PZCL  : NUMERICAL FACTOR OF COLLISION  (NOT USED IN THIS VERSION)
!     IHARM : MAXIMUM HARMONIC NUMBER (-IHARM TO IHARM)
!               (APPLICABLE FOR NMODEL GE 2.
!                WHEN NMODEL EQ 4 OR NMODEL EQ 5,
!                   IF IHARM GE 3 THEN FAST WAVE FLR FOR 3..IHARM
!                   IF IHARM LT 0 THEN FAST WAVE FLR FOR 0..ABS(IHARM))
!     APRFPN: PROFILE FACTOR FOR PLASMA DENSITY
!     APRFTR: PROFILE FACTOR FOR PARALLEL TEMPERATURE
!     APRFTP: PROFILE FACTOR FOR PERPENDICULAR TEMPERATURE
!
!     NXP   : NUMBER OF X-MESHES IN PLASMA
!     NXV   : NUMBER OF X-MESHES IN VACUUM (EVEN NUMBER)
!     NZP   : NUMBER OF Z-MESHES  (POWER OF 2, 2**N)
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
!     NLOOP : NUMBER OF LOOP : DRF  = INCREMENT OF RF
!                              DRKZ = INCREMENT OF RKZ
!
!     NSYM  : 0 : WITHOUT SYMMETRY IN Z-DIRECTION
!             1 : WITH    SYMMETRY
!            -1 : WITH   ASYMMETRY
!
!     NMODEL: 0 : NO FLR CONDUCTIVITY MODEL        (FINITE ELEMENT)
!             1 : NO FLR CONDUCTIVITY MODEL        (MULTI LAYER)
!             2 : FAST WAVE FLR CONDUCTIVITY MODEL (FINITE ELEMENT)
!             3 : FAST WAVE FLR CONDUCTIVITY MODEL (MULTIPLE LAYER)
!             4 : DIFFERENTIAL CONDUCTIVITY MODEL  (FINITE ELEMENT)
!             5 : DIFFERENTIAL CONDUCTIVITY MODEL  (MULTIPLE LAYER)
!             6 : INTEGRAL CONDUCTIVITY MODEL      (FINITE ELEMENT)
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
!
    USE w1comm
    IMPLICIT NONE
    INTEGER:: ns, na

!   ======( PROGRAM CONTROL )======
    NPRINT= 0
    NFILE = 0
    NGRAPH= 1
    NLOOP = 1
    DRF   = 0.D0
    DRKZ  = 0.D0
    NSYM  = 0
    DXFACT= 0.D0
    DXWDTH= 4.D0
    NALPHA= 0
    APRFPN= 0.5D0
    APRFTR= 1.D0
    APRFTP= 1.D0
    NDMAX = 100
    XDMAX = 10.D0
    NMODEL= 5
    NSYS  = 0
    NDISP = 0
    NXABS = 0
!     ======( CURRENT DRIVE PARAMETER ) =====
    ZEFF=2.D0
    NCDTYP=1
    WVYSIZ=0.D0
!     ======( MESH SIZE )======
    NXPMAX= 100
    NXVMAX= 10
    NZPMAX= 1
!     ======( MACHINE PARAMETER )======
    RR    = 3.000D0
    RB    = 1.3D0
    BB    = 3.5D0
    RA    = 1.2D0
    WALLR = 0.0D-5
    EPSH  = -0.3D0
!     ======( ANTENNA PARAMETER )======
    RF    = 52.5D0
    RD    = 1.25D0
    NAMAX =   1
    RKZ   = 8.00D0

    DO NA=1,NAM
       AJYH( NA )  =   0.D0
       ALYH( NA )  =   0.D0
       APYH( NA )  =   0.D0
       AJYL( NA )  =   0.D0
       ALYL( NA )  =   0.D0
       APYL( NA )  =   0.D0
       AJZH( NA )  =   0.D0
       AJZL( NA )  =   0.D0
    END DO
    AJYL( 1 )  =   1.D0
!     ======( PLASMA PARAMETER )======
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
    RETURN
  END SUBROUTINE w1_init
END MODULE w1init
