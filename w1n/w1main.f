C
C     ************************* TASK.W1 *************************
C
C     ICRF WAVE PROPAGATION AND ABSORPTION
C
C          DIFFERENTIAL OR INTEGRO-DIFFERENTIAL ANALYSIS
C          SLAB PLASMA MODEL
C          MULTIPLE LAYER OR FINITE ELEMENT METHOD
C
C     PROGRAMED BY
C
C          A. FUKUYAMA, S. NISHIYAMA, T. MATSUISHI AND k. HAMAMATSU(*)
C
C          DEPARTMENT OF ELECTRICAL AND ELECTRONIC ENGINEERING
C          FACULTY OF ENGINEERING
C          OKAYAMA UNIVERSITY
C          OKAYA 700, JAPAN
C
C      (*) DEPARTMENT OF LARGE TOKAMAK RESEARCH
C          NAKA FUSION RESEARCH ESTABLISHMENT
C          JAPAN ATOMIC ENERGY RESEARCH INSTITUTE
C          NAKA, IBARAKI 311-01, JAPAN
C
C     PROGRAM VERSION
C
C          V1.68 : 90/02/17
C          V1.70 : 90/03/14
C          V1.71 : 90/06/08 : NCDTYP ADDED
C          V1.72 : 90/06/16 : WVYSIZ ADDED
C          V1.73 : 90/06/18 : W1CLCD CORRECTED, WVYSIZ CORRECTED
C          V2.00 : 90/09/08 : B, A, HELICITY, ABSORBED BOUNDARY ADDED
C          V2.01 : 90/09/12 : HELICITY CURRENT ADDED
C          V2.02 : 90/09/14 : CORRECTION
C          V2.03 : 90/09/20 : TALPHA,TAU CORRECTED
C          V2.04 : 90/09/22 : NEW GRAPH, QH RECORRECTED
C          V2.05 : 90/10/16 : TRAPPED PARTICLE EFFECTS
C          V2.06 : 91/03/29 : NONRESONANT FORCE MODIFIED
C          V2.10 : 91/06/18 : UNDERFLOW ADJUSTMENT FOR SX
C          V2.11 : 92/03/18 : WORKING ON HP WS
C          V2.12 : 94/01/12 : SPLIT SOURCE FILES, INCLUDE FILE, ERF
C
C     ******* INPUT PARAMETERS THROUGH THE NAMELIST /W1/ *******
C
C     BB    : MAGNETIC FIELD AT X=0   (T)
C     RR    : PLASMA MAJOR RADIUS     (M)
C     RZ    : PERIODIC LENGTH         (M)  (IN Z DIRECTION)
C                           *** IF RZ.EQ.0. THEN RZ=2*PI*RR ***
C     RA    : PLASMA MINOR RADIUS     (M)
C     RD    : ANTENNA RADIUS          (M)
C     RB    : WALL RADIUS             (M)
C     EPSH  : HELICAL RIPPLE IN HELICAL SYSTEM
C
C     WALLR : WALL RESISTIVITY        (OHM-M)
C     NXABS : 0 : NO ABSORBING WALL
C         .GT.0 :    ABSORBING WALL AT NX=NXABS
C         .LT.0 :    ABSORBING WALL AT NX=ABS(NXABS) FOR BACKWARD WAVE
C
C     NCDTYP: CURRENT DRIVE TYPE (0:LANDAU, 1:TTMP)
C     ZEFF  : Z EFFECTIVE FOR CURRENT DRIVE EFFICIENCY
C
C     RF    : WAVE FREQUENCY          (MHZ)
C     RKZ   : WAVE NUMBER IN Z-DIRECTION (/M) (VALID FOR NZP=1)
C     NAMAX : NUMBER OF ANTENNAS
C     AJYH  : ANTENNA     CURRENT Y AT -RD  (KA)
C     AJYL  :                     Y AT  RD  (KA)
C     ALYH  : POSITION OF CURRENT Y AT -RD ( DEGREE )  -180..180
C     ALYL  :                     Y AT  RD ( DEGREE )  -180..180
C     APYH  : PHASE    OF CURRENT Y AT -RD ( DEGREE )
C     APYL  :                     Y AT  RD ( DEGREE )
C     AJZH  : ANTENNA     CURRENT Z AT -RD  (KA)
C     AJZL  :                     Z AT  RD  (KA)
C     DXFACT: MESH ACCUMULATION FACTOR (DEFAULT 0.D0)
C     DXWDTH:                   WIDTH  (DEFAULT 3.D0 : 3*LARMOR RADIUS)
C
C     NSMAX : NUMBER OF PARTICLE SPECIES
C     PA    : ATOMIC NUMBER
C     PZ    : CHARGE NUMBER
C     PN    : DENSITY AT X=0             (1.E20 M**-3)
C     PTPP  : PERPENDICULAR TEMPERATURE AT X=0   (KEV)
C     PTPR  : PARALLEL TEMPERATURE      AT X=0   (KEV)
C     PU    : PARALLEL DRIFT VELOCITY   AT X=0   (KEV)
C     PNS   : DENSITY    ON SURFACE X=RA (1.E20 M**-3)
C     PTS   : TEPERATURE ON SURFACE X=RA         (KEV)
C     PZCL  : NUMERICAL FACTOR OF COLLISION  (NOT USED IN THIS VERSION)
C     IHARM : MAXIMUM HARMONIC NUMBER (-IHARM TO IHARM)
C               (APPLICABLE FOR NMODEL GE 2.
C                WHEN NMODEL EQ 4 OR NMODEL EQ 5,
C                   IF IHARM GE 3 THEN FAST WAVE FLR FOR 3..IHARM
C                   IF IHARM LT 0 THEN FAST WAVE FLR FOR 0..ABS(IHARM))
C     APRFPN: PROFILE FACTOR FOR PLASMA DENSITY
C     APRFTR: PROFILE FACTOR FOR PARALLEL TEMPERATURE
C     APRFTP: PROFILE FACTOR FOR PERPENDICULAR TEMPERATURE
C
C     NXP   : NUMBER OF X-MESHES IN PLASMA
C     NXV   : NUMBER OF X-MESHES IN VACUUM (EVEN NUMBER)
C     NZP   : NUMBER OF Z-MESHES  (POWER OF 2, 2**N)
C
C     NPRINT: 0 :    LP OUTPUT (GLOBAL DATA)
C             1 :    LP OUTPUT (1-D DATA K-DEPENDENCE)
C             2 :    LP OUTPUT (1-D DATA X-DEPENDENCE)
C             3 :    LP OUTPUT (2-D DATA ELECTRIC FIELD)
C
C     NFILE : 0 : NO FILE OUTPUT
C             1 :    FILE OUTPUT (1-D DATA)
C             2 :    FILE OUTPUT (2-D DATA)
C
C     NGRAPH: 0 : NO GRAPHICS
C             1 :    GRAPHIC (E,J,P: 5 FIGS)
C             2 :    GRAPHIC (E,P: 4 FIGS)
C             3 :    GRAPHIC (E,P/TWICE HEIGHT: 4 FIGS)
C            +4 :    GRAPHIC (B FIELD, VECTOR POTENTIAL: 6 FIGS)
C            +8 :    GRAPHIC (HELICITY AND FORCE: 4 FIGS)
C           +16 :    GRAPHIC (HELICITY AND CURRENT: 4 FIGS)
C           +32 :    GRAPHIC (E,P,J: 4+1 FIGS)
C           +64 :    GRAPHIC (HELICITY, FORCE AND CURRENT: 8FIGS)
C
C     NLOOP : NUMBER OF LOOP : DRF  = INCREMENT OF RF
C                              DRKZ = INCREMENT OF RKZ
C
C     NSYM  : 0 : WITHOUT SYMMETRY IN Z-DIRECTION
C             1 : WITH    SYMMETRY
C            -1 : WITH   ASYMMETRY
C
C     NMODEL: 0 : NO FLR CONDUCTIVITY MODEL        (FINITE ELEMENT)
C             1 : NO FLR CONDUCTIVITY MODEL        (MULTI LAYER)
C             2 : FAST WAVE FLR CONDUCTIVITY MODEL (FINITE ELEMENT)
C             3 : FAST WAVE FLR CONDUCTIVITY MODEL (MULTIPLE LAYER)
C             4 : DIFFERENTIAL CONDUCTIVITY MODEL  (FINITE ELEMENT)
C             5 : DIFFERENTIAL CONDUCTIVITY MODEL  (MULTIPLE LAYER)
C             6 : INTEGRAL CONDUCTIVITY MODEL      (FINITE ELEMENT)
C
C     NALPHA: 0 : NO ALPHA PARTICLE EFFECT
C             1 : SLOWING DOWN DISTRIBUTION FOR ALPHA (NS=4, NMODEL=1)
C             2 : ALPHA DENSITY BY FUSION REACTION (NS 2:D 3:T 4:ALPHA)
C             3 : ALPHA DENSITY AND SLOWING DOWN DISTRIBUTION
C
C     NSYS  : 0 : TOKAMAK CONFIGURATION
C             1 : HELICAL CONFIGURATION
C
C     NDISP : 0 : FIELD CALCULATION
C             1 : WAVE NUMBER DISPLAY
C
C     XDMAX : MAXIMUM VALUE OF GYRORADIUS IN INTEGRATION (NMODEL=6)
C     NDMAX : MAXIMUM NUMBER OF INTEGRAL TABLE (NMODEL=6)
C
C     *************************************************************
C
      USE w1comm
      USE w1comm_alf1
      USE w1comm_dec
      USE w1comm_gdat

      IMPLICIT NONE
      INTEGER,SAVE:: IGRAPH=0
      INTEGER:: NS,NA,NX,NZ,NL,IC,ICL,IERR
      REAL(rkind):: PKXMIN,PKXMAX,TT1,TT2,XMIN,XMAX
      REAL(rkind):: RFSAVE,RKSAVE
C
C     ************************************************************
C
      NAMELIST /W1/ BB,RR,RZ,RA,RD,RB,RF,WALLR,APRFPN,APRFTR,APRFTP,
     &              RKZ,DRF,DRKZ,DXFACT,DXWDTH,
     &              AJYH,ALYH,APYH,AJYL,ALYL,APYL,AJZH,AJZL,NAMAX,
     &              PA,PZ,PN,PTPP,PTPR,PU,PNS,PTS,PZCL,NSMAX,
     &              NXP,NXV,NZP,NPRINT,NFILE,NGRAPH,NLOOP,NSYM,
     &              NMODEL,NALPHA,NDMAX,XDMAX,IHARM,NSYS,NDISP,
     &              EPSH,ZEFF,WVYSIZ,NCDTYP,NXABS,IELEC
C
C     ******* INITNALIZATION *******
      WRITE(6,600)
      WRITE(6,601) 'TASK/W1 --- V3.00 : 2014/09/18',
     &             NXPM,NXVM,NZPM,NSM,NAM,MATLM,NHARMM,NDM
      CALL w1comm_init
      CALL w1comm_alf1_init
      CALL w1comm_dec_init
      CALL w1comm_gdat_init
C
      CALL W1INIT
      PKXMIN=-1000.D0
      PKXMAX= 1000.D0
C
C     ******* PARAMETER INPUT *******
C
      OPEN(15,FILE='w1parm',STATUS='old',FORM='formatted',ERR=1)
      READ(15,W1,ERR=1,END=1)
      WRITE(6,*) '## &W1 UPDATED BY w1data'
C
    1 WRITE( 6,602)
    2 WRITE( 6,*) '## INPUT : &W1'
      RZ=0.D0
      READ(5,W1,ERR=1,END=9000)
         IF(NXP.LE.0) GO TO 9000
C
      IF(MOD(NXV,2).NE.0) NXV=NXV+1
      IF(ABS(RZ).LE.1.D-8) RZ=2.D0*PI*RR
      IF(NDMAX.GT.NDM-1) NDMAX=NDM-1
      NHARM=0
      IF(NMODEL.EQ.6) THEN
         DO 11 NS=1,NSMAX
            IF(ABS(IHARM(NS)).GT.NHARMM) IHARM(NS)=NHARMM
            IF(ABS(IHARM(NS)).GT.NHARM)  NHARM=IHARM(NS)
   11    CONTINUE
         DXD=XDMAX/NDMAX
         CALL W1QTBL(NDMAX,XDMAX,NHARM)
      ENDIF
C
      CALL FCLOCK(TT1)
      WRITE(6,609) BB,RR,RZ,RA,RD,RB
      WRITE(6,610) RF,RKZ,WALLR,ZEFF,WVYSIZ,NCDTYP,DXFACT,DXWDTH,EPSH
      WRITE(6,623) APRFPN,APRFTR,APRFTP
      WRITE(6,611) NXP,NXV,NZP,NSYM,NPRINT,NFILE,NGRAPH,NLOOP,
     &             NMODEL,NALPHA,NSYS,NXABS
      WRITE(6,620)
      DO 20 NA = 1 , NAMAX
         WRITE(6,621) NA,AJYH(NA),AJZH(NA),ALYH(NA),APYH(NA)
   20 CONTINUE
      WRITE(6,622)
      DO 30 NA = 1 , NAMAX
         WRITE(6,621) NA,AJYL(NA),AJZL(NA),ALYL(NA),APYL(NA)
   30 CONTINUE
      WRITE(6,612)
      DO 10 NS=1,NSMAX
         WRITE(6,613) NS,PA(NS),PZ(NS),PN(NS),PTPP(NS),PTPR(NS),PU(NS)
   10 CONTINUE
      WRITE(6,614)
      DO 12 NS=1,NSMAX
         WRITE(6,615) PNS(NS),PTS(NS),PZCL(NS),IHARM(NS),IELEC(NS)
   12 CONTINUE
C
C     ******* DISPERSION RELATION *******
C
      IF(NDISP.EQ.1) THEN
         XMIN=-RA
         XMAX= RA
         IF(NGRAPH.GT.0) THEN
            IF(IGRAPH.EQ.0) CALL GSOPEN
            IGRAPH=1
         ENDIF
  100    WRITE(6,*) '## INPUT : XMIN,XMAX,KXMIN,KXMAX,NXP ?'
         READ(5,*,ERR=100,END=2) XMIN,XMAX,PKXMIN,PKXMAX,NXP
         IF(ABS(XMAX-XMIN).LE.1.D-32) GOTO 2
         DO 110 NX=1,NXP
            XAM(NX)=XMIN+(XMAX-XMIN)*(NX-0.5)/NXP
  110    CONTINUE
         CALL W1PROF
         NZ=1
         CALL W1DSPA(NALPHA)
         IF(NMODEL.LE.3) THEN
            CALL W1WKXB
         ELSE
            CALL W1WKXD
         ENDIF
         CALL W1GDSP(PKXMIN,PKXMAX)
         GOTO 100
      ENDIF
C
C     ******* 2-DIMENSIONAL ANALYSIS *******
C
      RFSAVE=RF
      RKSAVE=RKZ
      DO 2000 NL=1 , NLOOP
         IF(NLOOP.NE.1) WRITE(6,630) RF,RKZ
         CALL W1SETZ(IERR)
            IF(IERR.NE.0) GOTO 2000
         CALL W1ANTS
         CALL W1SETX(IERR)
            IF(IERR.NE.0) GOTO 2000
         CALL W1PROF
         CALL W1PWRI
C
C     ******* FOURIER TRANSFORM OF ANTENNA CURRENT *******
C
         CALL W1FFTL(CJ1,NZP,0)
         CALL W1FFTL(CJ2,NZP,0)
         CALL W1FFTL(CJ3,NZP,0)
         CALL W1FFTL(CJ4,NZP,0)
C
C     ******* CALCULATION FOR EACH KZ *******
C
         DO 1000 NZ = 1 , NZP
            IF(NZ.LE.(NZP/2+1).OR.NSYM.EQ.0) THEN
               RKZ   = AKZ(NZ)
               CFJY1 = CJ1(NZ)
               CFJY2 = CJ2(NZ)
               CFJZ1 = CJ3(NZ)
               CFJZ2 = CJ4(NZ)
               CALL W1BCND
               IF(NMODEL.LE.5) THEN
                  CALL W1DSPA(NALPHA)
               ELSE
                  CALL W1DSPQ(ICL)
               ENDIF
               IF(NMODEL.EQ.0.OR.NMODEL.EQ.2) THEN
                  CALL W1BNDA(IERR)
                     IF(IERR.NE.0) GOTO 2000
                  CALL W1EPWA(NZ)
               ELSEIF(NMODEL.EQ.1.OR.NMODEL.EQ.3) THEN
                  CALL W1WKXB
                  CALL W1BNDB(IERR,NXABS)
                     IF(IERR.NE.0) GOTO 2000
                  CALL W1EPWB(NZ)
                  CALL W1HELD(4)
               ELSEIF(NMODEL.EQ.4) THEN
                  CALL W1BNDC(IERR)
                     IF(IERR.NE.0) GOTO 2000
                  CALL W1EPWC(NZ)
               ELSEIF(NMODEL.EQ.5) THEN
                  CALL W1WKXD
                  CALL W1BNDD(IERR,NXABS)
                     IF(IERR.NE.0) GOTO 2000
                  CALL W1EPWD(NZ)
                  CALL W1HELD(6)
               ELSEIF(NMODEL.EQ.6) THEN
                  CALL W1BNDQ(IERR)
                     IF(IERR.NE.0) GOTO 2000
                  CALL W1EPWQ(NZ)
               ENDIF
               CALL W1EVAC(NZ,NSYM)
               CALL W1CLCD(NZ)
               CALL W1CLPW(NZ,NSYM)
            ELSE
               CALL W1SYMS(NZ,NSYM)
            ENDIF
 1000    CONTINUE
         IF(NMODEL.EQ.6) WRITE(6,616) MATL,MATLM,ICL,NCLM,NDMAX,XDMAX
C
C     ******* INVERSE FOURIER TRANSFORM *******
C
         CALL W1FFTL(CJ1,NZP,1)
         CALL W1FFTL(CJ2,NZP,1)
         CALL W1FFTL(CJ3,NZP,1)
         CALL W1FFTL(CJ4,NZP,1)
C
         DO 1100 NX=1,NXT
         DO 1100 IC=1,3
            CALL W1FFTL(CE2DA(1,NX,IC),NZP,1)
 1100    CONTINUE
C
C     ******* POWER ABSORPTION AND OUTPUT *******
C
         CALL W1PWRS
         CALL W1PRNT(NPRINT)
         CALL W1FILE(NFILE)
         IF(NGRAPH.GT.0) THEN
            IF(IGRAPH.EQ.0) CALL GSOPEN
            IF(NZP.EQ.1) THEN
               CALL W1GR1D(MOD(NGRAPH,   4))
               CALL W1GR1B(MOD(NGRAPH/4, 2))
               CALL W1GR1F(MOD(NGRAPH/8, 2))
               CALL W1GR1H(MOD(NGRAPH/16,2))
               CALL W1GRUD(MOD(NGRAPH/32,2))
               CALL W1GRUF(MOD(NGRAPH/64,2))
            ELSE
               CALL W1GR2D(MOD(NGRAPH,   4))
               CALL W1GR1D(MOD(NGRAPH,   4))
            ENDIF
            IGRAPH=1
         ENDIF
         IF(NLOOP.NE.1) THEN
            RF   = RF   + DRF
            RKZ  = RKZ  + DRKZ
         ENDIF
 2000 CONTINUE
C
      RF =RFSAVE
      RKZ=RKSAVE
      CALL FCLOCK(TT2)
      WRITE(6,699) (TT2-TT1)
  699 FORMAT(1H ,'## CPU TIME = ',F15.3,' SEC')
      GOTO 2
C
 9000 CONTINUE
      IF(IGRAPH.NE.0) CALL GSCLOS
C      WRITE(6,*) 'UNDERFLOW COUNT (4,8) = ',IUNFL4,IUNFL8
      STOP
C
  600 FORMAT(1H1)
  601 FORMAT(1H ,'## ',A28,' ',
     &           '## NXPM,NXVM,NZPM,NSM,NAM,MATLM,NHARMM,NDM ##'/
     &       1H ,35X,I4,',',I4,',',I4,',',I3,',',I3,',',
     &               I5,',',I6,',',I3)
  602 FORMAT(1H ,'## NAMELIST PARM : BB,RR,RZ,RA,RD,RB,WALLR,EPSH,',
     &           'APRFPN,APRFTR,APRFTP'/
     &       1H ,19X,'AJYH,ALYH,APYH,AJYL,ALYL,APYL,AJZH,AJZL:(NAMAX)'/
     &       1H ,19X,'PA,PZ,PN,PTPP,PTPR,PU,PNS,PTS,PZCL,IHARM:(NSMAX)'/
     &       1H ,19X,'RF,RKZ,DRF,DRKZ,DXFACT,DXWDTH,NXABS '/
     &       1H ,19X,'NXP,NXV,NZP,NPRINT,NFILE,NGRAPH,NLOOP,NSYM,',
     &           'NSYS,NDISP,IELEC'/
     &       1H ,19X,'NMODEL,NALPHA,NDMAX,XDMAX,NCDTYP,WVYSIZ,ZEFF')
  609 FORMAT(1H ,'BB    =',1PD12.4,' (T)   ',
     &           'RR    =',  D12.4,' (M)   ',
     &           'RZ    =',  D12.4,' (M)'/
     &       1H ,'RA    =',1PD12.4,' (M)   ',
     &           'RD    =',  D12.4,' (M)   ',
     &           'RB    =',  D12.4,' (M)')
  610 FORMAT(1H ,'RF    =',1PD12.4,' (MHZ) ',
     &           'RKZ   =',1PD12.4,' (/M)  ',
     &           'WALLR =',1PD12.4,'(OHM M)'/
     &       1H ,'ZEFF  =',1PD12.4,'       ',
     &           'WVYSIZ=',1PD12.4,' (M)   ',
     &           'NCDTYP=',I6,6X,  '       '/
     &       1H ,'DXFACT=',1PD12.4,'       ',
     &           'DXWDTH=',1PD12.4,' (M)   ',
     &           'EPSH  =',1PD12.4,'       ')
  620 FORMAT(1H ,'HFS  ANTY(KA)  ANTZ(KA)  POS(DEG)  PHA(DEG)     ')
  621 FORMAT(1H ,I3,4(F10.4))
  622 FORMAT(1H ,'LHS  ANTY(KA)  ANTZ(KA)  POS(DEG)  PHA(DEG)     ')
  623 FORMAT(1H ,'APRFPN=',1PD12.4,'       ',
     &           'APRFTR=',1PD12.4,'       ',
     &           'APRFTP=',1PD12.4,'       ')
  611 FORMAT(1H ,'NXP   =',I6,6X,'NXV   =',I6,
     &        6X,'NZP   =',I6,6X,'NSYM  =',I6/
     &       1H ,'NPRINT=',I6,6X,'NFILE =',I6,
     &        6X,'NGRAPH=',I6,6X,'NLOOP =',I6/
     &       1H ,'NMODEL=',I6,6X,'NALPHA=',I6,
     &        6X,'NSYS  =',I6,6X,'NXABS =',I6)
  612 FORMAT(1H ,'NS',2X,'PA(MASS)    PZ(CHARGE)  PN(E20/M3)  ',
     &                   'PTPP(KEV)   PTPR(KEV)   PU(KEV)     ')
  613 FORMAT(1H ,I2,1P6D12.4)
  614 FORMAT(1H ,    28X,'PNS(E20/M3) PTS(KEV)    PZCL    IHARM',
     &                   '   IELEC')
  615 FORMAT(1H ,26X,1P3D12.4,I3,5X,I3)
  616 FORMAT(1H ,'## INTEGRO-DIFF EQ. : ',
     &           'MATL  MATLM  ICL   NCLM  NDMAX XDMAX'/
     &       1H ,22X,I3,3X,I3,3X,I5,1X,I6,I3,1P1D12.4)
  630 FORMAT(1H /
     &       1H ,'** RF = ',1PD12.4,' (MHZ) , RKZ = ',1PD12.4,
     &           ' (/M) **')
      END
