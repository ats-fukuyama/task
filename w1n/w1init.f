!     $Id$

C
C     ****** TASK/W1 ******
C
C     ****** DEFAULT INPUT PARAMETERS ******
C
      SUBROUTINE W1INIT
C
      USE w1comm
      IMPLICIT NONE
      INTEGER NS,NA
C
C     ======( PROGRAM CONTROL )======
      NPRINT= 0
      NFILE = 0
      NGRAPH= 1
      NLOOP = 1
      NSYM  = 0
      NALPHA= 0
      DRF   = 0.D0
      DRKZ  = 0.D0
      DXFACT= 0.D0
      DXWDTH= 4.D0
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
      NCDTYP=1
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
      NAMAX =   1
      RKZ   = 8.00  D0
      DO 10 NA=1,NAM
         AJYH( NA )  =   0.    D 0
         ALYH( NA )  =   0.    D 0
         APYH( NA )  =   0.    D 0
         AJYL( NA )  =   0.    D 0
         ALYL( NA )  =   0.    D 0
         APYL( NA )  =   0.    D 0
         AJZH( NA )  =   0.    D 0
         AJZL( NA )  =   0.    D 0
   10 CONTINUE
      AJYL( 1 )  =   1.    D 0
C     ======( PLASMA PARAMETER )======
      NSMAX = 3
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
      DO 20 NS=4,NSM
         PA  ( NS ) =  1.     D 0
         PZ  ( NS ) =  1.     D 0
         PN  ( NS ) =  0.     D 0
         PTPP( NS ) =  1.     D 0
         PTPR( NS ) =  1.     D 0
         PNS ( NS ) =  0.     D 0
         PTS ( NS ) =  0.1    D 0
         PU  ( NS ) =  0.     D 0
         PZCL( NS ) =  0.     D 0
         IHARM(NS ) =  2
         IELEC(NS ) =  0
   20 CONTINUE
      RETURN
      END
