C
C     ****** TASK/W1 ******
C
C     ****** DEFAULT INPUT PARAMETERS ******
C
      SUBROUTINE W1INIT
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM3/ PA(ISM),PZ(ISM),PN(ISM),PTPP(ISM),PTPR(ISM),
     &                PNS(ISM),PTS(ISM),PU(ISM),PZCL(ISM)
      COMMON /W1PRM4/ RZ,DZ,NZP
      COMMON /W1PRM6/ IELEC(ISM)
      COMMON /W1CTRL/ DRF,DRKZ,DXFACT,DXWDTH,APRFPN,APRFTR,APRFTP,
     &                NPRINT,NFILE,NGRAPH,NLOOP,NSYM,NFLR,NALPHA,
     &                NSYS,NDISP,NXABS
      COMMON /W1QCTL/ XDMAX,DXD,NDMAX,NHARM,IHARM(ISM),MATL,NMODEL
      COMMON /W1ANT1/ AJYL(IAM),AJYH(IAM),AJZL(IAM),AJZH(IAM),
     &                ALYL(IAM),ALYH(IAM),APYL(IAM),APYH(IAM)
      COMMON /W1PRM5/ EPSH
      COMMON /W1CDDT/ AJCDX(NXPM),AJCDK(NZPM),AJCDT,ZEFF,WVYSIZ,ICDTYP
C
C     ======( PHYSICAL CONSTANTS )======
      AEE   = 1.6021892  D-19
      AME   = 9.109534   D-31
      AMM   = 1.6726485  D-27
      VC    = 2.99792458 D  8
      EPS0  = 8.854187818D-12
      AMYU0 = 1.256637061D- 6
      PI    = 3.141592654D  0
C     ======( PROGRAM CONTROL )======
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
C     ======( CURRENT DRIVE PARAMETER ) =====
      ZEFF=2.D0
      ICDTYP=1
      WVYSIZ=0.D0
C     ======( MESH SIZE )======
      NXP   = 100
      NXV   = 10
      NZP   = 1
C     ======( MACHINE PARAMETER )======
      RR    = 3.000 D 0
      RB    = 1.3   D 0
      BB    = 3.5   D 0
      RA    = 1.2   D 0
      WALLR = 0.    D-5
      EPSH  = -0.3D0
C     ======( ANTENNA PARAMETER )======
      RF    = 52.5  D 0
      RD    = 1.25  D 0
      IAMAX =   1
      RKZ   = 8.00  D0
      DO 10 IA=1,IAM
         AJYH( IA )  =   0.    D 0
         ALYH( IA )  =   0.    D 0
         APYH( IA )  =   0.    D 0
         AJYL( IA )  =   0.    D 0
         ALYL( IA )  =   0.    D 0
         APYL( IA )  =   0.    D 0
         AJZH( IA )  =   0.    D 0
         AJZL( IA )  =   0.    D 0
   10 CONTINUE
      AJYL( 1 )  =   1.    D 0
C     ======( PLASMA PARAMETER )======
      ISMAX = 3
C
      PA  ( 1 ) =  5.4466 D-4
      PZ  ( 1 ) = -1.     D 0
      PN  ( 1 ) =  1.000  D 0
      PTPP( 1 ) =  2.000  D 0
      PTPR( 1 ) =  2.000  D 0
      PNS ( 1 ) =  0.1    D 0
      PTS ( 1 ) =  0.1    D 0
      PU  ( 1 ) =  0.     D 0
      PZCL( 1 ) =  0.     D 0
      IHARM(1 ) =  2
      IELEC(1 ) =  1
C
      PA  ( 2 ) =  2.     D 0
      PZ  ( 2 ) =  1.     D 0
      PN  ( 2 ) =  0.950  D 0
      PTPP( 2 ) =  2.000  D 0
      PTPR( 2 ) =  2.000  D 0
      PNS ( 2 ) =  0.095  D 0
      PTS ( 2 ) =  0.1    D 0
      PU  ( 2 ) =  0.     D 0
      PZCL( 2 ) =  0.     D0
      IHARM(2 ) =  2
      IELEC(2 ) =  0
C
      PA  ( 3 ) =  1.     D 0
      PZ  ( 3 ) =  1.     D 0
      PN  ( 3 ) =  0.050  D 0
      PTPP( 3 ) =  2.000  D 0
      PTPR( 3 ) =  2.000  D 0
      PNS ( 3 ) =  0.005  D 0
      PTS ( 3 ) =  0.1    D 0
      PU  ( 3 ) =  0.     D 0
      PZCL( 3 ) =  0.     D 0
      IHARM(3 ) =  2
      IELEC(3 ) =  0
C
      DO 20 IS=4,ISM
         PA  ( IS ) =  1.     D 0
         PZ  ( IS ) =  1.     D 0
         PN  ( IS ) =  0.     D 0
         PTPP( IS ) =  1.     D 0
         PTPR( IS ) =  1.     D 0
         PNS ( IS ) =  0.     D 0
         PTS ( IS ) =  0.1    D 0
         PU  ( IS ) =  0.     D 0
         PZCL( IS ) =  0.     D 0
         IHARM(IS ) =  2
         IELEC(IS ) =  0
   20 CONTINUE
      RETURN
      END
