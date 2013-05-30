C     $Id$
C
C     ****** INITIALIZE PARAMETERS ******
C
      SUBROUTINE WFINIT
C
      INCLUDE 'wfcomm.inc'
C
C     *** CONSTANTS ***
C
C        CI    : Imaginary unit
C        PI    : Pi
C        AEE   : Elementaty charge
C        AME   : Electron mass
C        AMP   : Proton mass
C        VC    : Speed of light in vacuum
C        RMU0  : Permeability of free space
C        EPS0  : Permittivity of free space
C
      CI     = (0.D0,1.D0)
      PI     = ACOS(0.D0)*2.D0
      AEE    = 1.60217733D-19
      AME    = 9.1093897D-31
      AMP    = 1.6726231D-27
      VC     = 2.997924580D8
      RMU0   = 4.D0*PI*1.D-7
      EPS0   = 1.D0/(VC*VC*RMU0)
C
C     *** CONTROL PARAMETERS ***
C
C        KFNAME: File name of element data
C        KFNAMA: File name of antenna data
C        KFNAMF: File name of field data
C        KFNAMZ: File name of zone data
C
      KFNAME = 'elm-data'
      KFNAMA = 'ant-data'
      KFNAMF = 'fld-data'
      KFNAMZ = 'zone-data'
C
C     *** CONFIGURATION PARAMETERS ***
C
C        BB    : Magnetic field at center                        (T)
C        RA    : Plasma minor radius                             (m)
C
      BB     = 0.08D0
      RA     = 0.08D0
C
C     *** CONFIGURATION PARAMETERS (MIRROR: MODELB=3,4) ***
C
C        RMIR  : Mirror ratio
C        ZBB   : Periodic length along magnetic axis           (m)
C
      RMIR   = 2.0D0
      ZBB    = 0.15D0
C
C     *** CONFIGURATION PARAMETERS (TOKAMAK: MODELB=5) ***
C
C        RR    : Plasma major radius                             (m)
C        Q0    : Safety factor at center
C        QA    : Safety factor on plasma surface
C        RKAP  : Plasma shape elongation
C        RDEL  : Plasma shape triangularity *
C
      RR     = 3.D0
      Q0     = 1.D0
      QA     = 3.D0
      RKAP   = 1.D0
      RDEL   = 0.D0
C
C     *** CONFIGURATION PARAMETERS (HELICAL: MODELB=6) ***
C
C        H1    : Helical pitch (2*pi/L) for B profile
C        H2    : Helical pitch (2*pi/L) for metric
C        RRC   : Coil radius                                   (m)
C
      H1     = 1.25D0
      H2     = 1.25D0
      RRC    = 0.95D0
C
C     *** CIRCULAR COIL PARAMETERS (MODELB=7) ***
C
C        NCMAX : Number of coil
C        RC    : Radial position of coil current               (m)
C        ZC    : Axial position of coil current                (m)
C        B2    : Magnetic filed on axis, center of coil        (T)
C
      NCMAX  = 3
      RC(1)  = 0.35D0
      ZC(1)  = 0.D0
      BC(1)  = 0.001D0
      RC(2)  = 0.35D0
      ZC(2)  = 0.05D0
      BC(2)  =-0.001D0
      RC(3)  = 0.35D0
      ZC(3)  =-0.05D0
      BC(3)  =-0.001D0
C
C     *** RF PARAMETERS ***
C
C        RF    : Wave frequency                               (MHz)
C        RKZ   : Wave number in (z or Z) direction          (m**-1)
C        NPHI  : Mode number in (phi) direction
C        NZMAX : Number of Fourier modes in (z, phi or Z) direction
C        RZ    : Periodic length in (z or Z) direction          (m)
C
      RF     = 2450.D0
      RKZ    = 2.5D0
      NPHI   = -1
C
      NZMAX  = 1
      RZ     = 1.D0
C
C     *** ANTENNA PARAMETERS ***
C
C        NAMAX : Number of antennae
C        AJ    : Antenna current density                       (A/m)
C        APH   : Antenna phase                              (degree)
C        AWD   : Antenna width in (z, phi, Z) direction     (degree)
C        APOS  : Antenna position in (z, phi, Z) direction  (degree)
C        APOS  : Antenna position in (z, phi, Z) direction  (degree)
C
      NAMAX=0
      DO 10 NA=1,NAM
         AJ(NA)   = 1.D0
         IF(MOD(NA,2).EQ.1) THEN
            APH(NA)  = 0.D0
         ELSE
            APH(NA)  = 180.D0
         ENDIF
         AWD(NA)  = 0.D0
         APOS(NA) = 0.D0
   10 CONTINUE
C
C     *** PLASMA PARAMETERS ***
C
C        NSMAX : Number of particle species
C        PA    : Mass number
C        PZ    : Charge number
C        PN    : Density at center                     (1.0E20/Mm*3)
C        PNS   : Density on plasma surface             (1.0E20/m**3)
C        PTPR  : Parallel temperature at center                (keV)
C        PTPP  : Perpendicular temperature at center           (keV)
C        PTS   : Temperature on surface                        (keV)
C        PZCL  : Ratio of collision frequency to wave frequency
C
      NSMAX= 2
      IF(NSMAX.GT.NSM) NSMAX=NSM
C
         PA(1)= AME/AMP
         PZ(1)=-1.0D0
         PN(1)= 0.0002D0
         PNS(1)= 0.D0
         PTPR(1)=0.01D0
         PTPP(1)=0.01D0
         PTS(1)=0.D0
         PZCL(1)= 0.02D0
C
      IF(NSM.GE.2) THEN
         PA(2)= 39.9480D0
         PZ(2)= 1.0D0
         PN(2)= 0.0002D0
         PNS(2)= 0.0D0
         PTPR(2)=0.01D0
         PTPP(2)=0.01D0
         PTS(2)=0.D0
         PZCL(2)= 0.02D0
      ENDIF
C
      IF(NSM.GE.3) THEN
         PA(3)= 1.0D0
         PZ(3)= 1.0D0
         PN(3)= 0.0D0
         PNS(3)= 0.0D0
         PTPR(3)=0.1D0
         PTPP(3)=0.1D0
         PTS(3)=0.D0
         PZCL(3)= 0.001D0
      ENDIF
C
C     *** CONTROL PARAMETERS ***
C
C        MODELS: 0 : No symmetry
C                1 : Axial symmetry (Y axis)
C                2 : Axial symmetry (Y axis -RR)
C                3 : Helical symmetry (Z axis)
C        MODELB: 0 : X axis
C                1 : Y axis
C                2 : Z axis
C                3 : Axisymmetric mirror
C                4 : Translational mirror
C                5 : Tokamak
C                6 : Helical
C                7 : Circular coils
C        MODELD: 0 : Cold plasma model
C                1 : Warm plasma model
C                2 : Hot plasma model
C        MODELP: 0 : Flat profile
C                1 : Radially parabolic profile
C                2 : Axially exponential profile
C                3 : Radially and axially parabolic profile
C                4 : Temporal use
C                5 : Radially parabolic and axially quartic profile
C
C        MODELW: 0 : Fixed density and fixed temperature on boundary
C                1 : Free density and fixed temperature on boundary
C                2 : Free Density and free temperature on boundary
C        MODELT: 0 : Fixed temperature model
C                1 : Density gradient model
C                2 : Pressure gradient model
C
C        MODELN: 0 : Fixed crosssection for electron collision with neutrals
C                1 : Mometum transder collision data
C
C        MODELV:   : Type of divide model
C
         MODELS = 1
         MODELB = 3
         MODELD = 0
         MODELP = 1
         MODELW = 0
         MODELT = 2
         MODELN = 0
         MODELV = 0
C
C     *** OUTPUT PARAMETERS ***
C
C        NPRINT: Print output parameter
C                0 - No output
C              * 1 - Parameter and global field data
C                2 - Local field data
C                3 - Element data
C        NDRAWD: Drawing parameter for elemendt divider
C                0 : Boundary shape
C              * 1 : Element shape
C                2 : Element shape + Element number
C                3 : Element shape + Element number + Node number
C        NDRAWA: Drawing parameter for antenna generater
C                0 : Antenna primary data
C                1 : Antenna secondary data
C              * 2 : Antenna secondary data + Element shape
C        NRMAX : Number of radial mesh points
C
      NPRINT = 1
      NDRAWD = 1
      NDRAWA = 2
      NRMAX  = 101
C
C     *** DIVIDER PARAMETERS ***
C
C        BXMIN : Minimum x (m)
C        BXMAX : Maximum x (m)
C        BYMIN : Minimum y (m)
C        BYMAX : Maximum y (m)
C        RB    : Boundary radius (m)
C        BKAP  : Boundary elongation
C        BDEL  : Boundary triangularity *
C        DELX  : Typical element size in x direction (m)
C        DELY  : Typical element size in y direction (m)
C
      BXMIN  = 0.0D0
      BXMAX  = 0.1D0
      BYMIN  =-0.15D0
      BYMAX  = 0.15D0
      RB     = 0.25D0
      BKAP   = 1.D0
      BDEL   = 0.D0
      DELX   = 0.01D0
      DELY   = 0.01D0
C
C     *** ANTENNA SHAPE PARAMETERS ***
C
C        PIN   : Input Power (W); Set 0.0 to calculate from antenna current
C        PHIW  : Potential of wave electrode                    (V)
C        RD    : Antenna radius (m)
C        THETJ1: Start angle of arc antenna (degree)
C        THETJ2: End angle of arc antenna (degree)
C        ZJH1  : Axial start position of helical antenna (m)
C        ZJH2  : Axial end position of helical antenna (m)
C        PHJH  : Rotation angle of helical antenna (degree)
C        NTYPJH: Type of helical antenna
C                0 : Loop antenas on both ends
C                1 : Loop antena on the second end
C                2 : Loop antena on the first end
C                3 : No loop antena
C        NJMAX : Number of primary grid points of antenna
C
      PIN    = 0.D0
      PHIW   = 0.D0
      RD     = 0.22D0
      THETJ1 = 40.D0
      THETJ2 =-40.D0
      ZJH1   = 0.01D0
      ZJH2   = 0.19D0
      PHJH   = 360.D0
      NTYPJH = 0
      NJMAX  = 41
C
C     *** TRANSPORT PARAMETERS ***
C
C     DT     : Time step size
C     NTMAX  : Iteration number
C     NSTEP  : Number of transport calculations after one wave calculation
C     NFMAX  : Number of particle species in transport calculation
C
      DT     = 1.D-6
      NTMAX  = 1
      NTSTEP = 10
      NFMAX  = 2
C
C     *** TRANSPORT PLASMA PARAMETERS ***
C
C     PPN0   : Neutral pressure (Pa)
C                 1 torr = 133.322 Pa
C     PTN0   : Initial neutral temperarure (eV)
C     PNE0   : Initial electron density (1.D20/m^3)
C     PTE0   : Initial electron temperature (eV)
C     PTI0   : Initial ion temperature (eV)
C     PNES   : Edge electron density (1.D20/m^3)
C     PTES   : Edge electron temperature (eV)
C     PTIS   : Edge ion temperature (eV)
C
      PPN0   = 1.D0
      PTN0   = 0.03D0
      PNE0   = 1.D-6
      PTE0   = 0.03D0
      PTI0   = 0.03D0
      PNES   = 1.D-6
      PTES   = 0.03D0
      PTIS   = 0.03D0
C
C     *** BOHM DIFFUSION ***
C
C     DC     : FACTOR OF BOHM DIFFUSION COEFFICIENT
C
      DC     = 1.D0
C
C     *** COMPUTATION PARAMETERS ***
C
C     EPSIMP : CONVERGENCE CRITERION FOR IMPLICIT TIME EVOLUTION
C     EPSSUM : BOUNDARY BETWEEN RELATIVE ERROR AND ABSOLUTE ERROR
C     MAXIMP : MAXIMUM LOOP COUNT FOR IMPLICIT TIME EVOLUTION
C     FACIMP : IMPLICIT FACTOR
C
      EPSIMP = 1.D-4
      EPSSUM = 1.D0
      MAXIMP = 1
      FACIMP = 1.D0
C
C     *** ARTIFICAL SOURCE ***
C
C     PGIVEN : MAXIMUM POWER DENSITY
C     SGIVEN : MAXIMUM PLASMA SOURCE DENSITY
C     XGIVEN : X COORDINATE OF THE CENTER OF THE SOURCE
C     YGIVEN : Y COORDINATE OF THE CENTER OF THE SOURCE
C     RGIVEN : DECAY LENGTH OF THE SOURCE
C
      PGIVEN = 0.D0
      SGIVEN = 0.D0
      XGIVEN = 0.05D0
      YGIVEN = 0.D0
      RGIVEN = 0.05D0
C
C     *** ELECTRODE PARAMETERS ***
C
C     RFES   : FREQUENCY OF ELECTRODE POTENTIAL
C     PHIES  : AMPLITUDE OF ELECTRODE POTENTIAL
C
      RFES   = 13.56D0
      PHIES  = 0.D0
C
C     *** GRAPHICS CONTROL PARAMETERS ***
C
C     KGINX : GRAPHIC CONTROL STRINGS
C     KGINV : GRAPHIC CONTROL STRINGS
C
      KGINX(0)='EXI EYI EZR PP1C'
      KGINX(1)='EXR EXI EYR EYI'
      KGINX(2)='EZR EZI PP1C PP2C'
      KGINX(3)='EXR EXI EYR EYI EZR EZI PP1C PP2C'
      KGINX(4)='AXR AXI AYR AYI AZR AZI AFR AFI'
      KGINX(5)='PF0C PNEC PTEC PTIC PIOC PC0C'
      KGINX(6)='PF0Y0 PFOY0.03 PF0X0 PF0X0.05'
      KGINX(7)='PNEY0 PNEY0.03 PNEX0 PNEX0.05'
      KGINX(8)='PTEY0 PTEY0.03 PTEX0 PTEX0.05'
      KGINX(9)='L2 TF0 TNE TTE TTI'
C
      KGINV(0)='EXR,EXI,EYR,EYI,EZR,EZI,PP1C,PP2C,PNEC,PTEC,TNE,TTE'
      KGINV(1)='PC1C PC2C PR1C PR2C'
      KGINV(2)='PU1C PU2C PV1C PV2C'
      KGINV(3)='PD1C PD2C PE1C PE2C'
      KGINV(4)='PH1C PH2C PK1C PK2C'
      KGINV(5)='PI0C'
      KGINV(6)='PI0C'
      KGINV(7)='PI0C'
      KGINV(8)='PI0C'
      KGINV(9)='PI0C'
C
C     *** GRAPHICS CONTROL PARAMETERS ***
C
C     NGRAPH: Type of 1D graphic output
C             1 : Autoscale plot
C             2 : Symmetric scale plot
C             3 : Amplitude and phase plot
C     NGRAPH: Type of 2D graphic output
C             1 : Contour plot (in element mesh)      
C             2 : Color-painted contour plot (in rectangular mesh)
C             3 : Bird's eye view
C             4 : Contour plot (in rectangular mesh)
C     FRATIO: Horizontal expansion factor for 2D graphics
C
      NGRAPH=1
      FRATIO=1.D0
C
C     *** 3D GRAPHICS CONTROL PARAMETERS ***
C
      GA=-25.0
      GB=  0.0
      GC=-30.0
      GD=  0.0
      GE= 1000.0
      GXN= 6.0
      GYN= 9.0
      GZN= 3.0
      GXN1=-5.0
      GXN2= 5.0
      GYN1= 0.0
      GYN2=10.0
      IXY= 3
      IDN=-3
C
      MODIFY=0
C
C     *** INITIALIZATION PARAMETERS (DO NOT MODIFY) ***
C
      NNOD   = 0
      NFOPEN = 0
      XDMIN  = 0.D0
      XDMAX  = 0.D0
C
      NDMAX=0
      NBWMAX=0
      NBPMAX=0
      NMMAX=1
      EPSDM(1)=1.D0
      RMUDM(1)=1.D0
      SIGDM(1)=0.D0
      NVWMAX=0
      NVPMAX=0
C
      RETURN
      END
C
C     ***** INPUT PARAMETER LIST *****
C
      SUBROUTINE WFPLST
C
C     --- show valid namelist variables ---
C
C
      WRITE(6,*) '&WF: BB,RA,RZ,RR,Q0,QA,RKAP,RDEL,RMIR,ZBB,H1,H2,RRC,'
      WRITE(6,*) '     RF,RKZ,NPHI,PHIW,AJ,APH,AWD,APOS,'
      WRITE(6,*) '     PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,'
      WRITE(6,*) '     NSMAX,NAMAX,NRMAX,NZMAX,NPRINT,NDRAWD,NDRAWA,'
      WRITE(6,*) '     KFNAME,KFNAMA,KFNAMF,'
      WRITE(6,*) '     BXMIN,BXMAX,BYMIN,BYMAX,'
      WRITE(6,*) '     RB,BKAP,BDEL,DELX,DELY,'
      WRITE(6,*) '     PIN,EPSIMP,FACIMP,EPSSUM,MAXIMP,'
      WRITE(6,*) '     RD,THETJ1,THETJ2,ZJH1,ZJH2,PHJH,NTYPJH,NJMAX,'
      WRITE(6,*) '     DT,NTMAX,NTSTEP,NFMAX,RFES,PHIES,MODELT,MODELN,'
      WRITE(6,*) '     PNE0,PTE0,PTI0,PNES,PTES,PTIS,PPN0,PTN0'
      WRITE(6,*) '     DC,PGIVEN,SGIVEN,XGIVEN,YGIVEN,RGIVEN'
      WRITE(6,*) '     NGRAPH,FRATIO,GA,GB,GC,GD,GE,GXN,GYN,GZN,'
      WRITE(6,*) '     GXN1,GXN2,GYN1,GYN2,IXY,IDN'
      RETURN
      END
C
C     ****** INPUT PARAMETERS ******
C
      SUBROUTINE WFPARM(KID)
C
      INCLUDE 'wfcomm.inc'
C
      NAMELIST /WF/ BB,RA,RZ,RR,Q0,QA,RKAP,RDEL,
     &              RMIR,ZBB,H1,H2,RRC,RC,ZC,BC,NCMAX,
     &              RF,RKZ,PHIW,NPHI,AJ,APH,AWD,APOS,
     &              PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,
     &              NSMAX,NAMAX,NRMAX,NZMAX,
     &              NPRINT,NDRAWD,NDRAWA,
     &              KFNAME,KFNAMA,KFNAMF,KFNAMZ,
     &              BXMIN,BXMAX,BYMIN,BYMAX,
     &              RB,BKAP,BDEL,DELX,DELY,
     &              PIN,EPSIMP,EPSSUM,FACIMP,MAXIMP,
     &              RD,THETJ1,THETJ2,ZJH1,ZJH2,PHJH,NTYPJH,NJMAX,
     &              DT,NTMAX,NTSTEP,NFMAX,RFES,PHIES,MODELT,
     &              PNE0,PTE0,PTI0,PNES,PTES,PTIS,PPN0,PTN0,
     &              DC,PGIVEN,SGIVEN,XGIVEN,YGIVEN,RGIVEN,
     &              MODELS,MODELB,MODELD,MODELP,MODELW,MODELT,MODELN,
     &              KGINX,KGINV,NGRAPH,FRATIO,
     &              GA,GB,GC,GD,GE,GXN,GYN,GZN,GXN1,GXN2,GYN1,GYN2,
     &              IXY,IDN,MODIFY
      CHARACTER KSNAME*32,KSNAMZ*32,KSNAMA*32,KSNAMF*32
      CHARACTER KPNAME*32,KLINE*70,KNAME*80,KID*1
      LOGICAL LEX
C
      MODE=0
 1000 CONTINUE
C
    1 CONTINUE
         WRITE(6,*) '## INPUT : &WF '
C
         KSNAME=KFNAME
         KSNAMZ=KFNAMZ
         KSNAMA=KFNAMA
         KSNAMF=KFNAMF
         READ(5,WF,IOSTAT=IST,ERR=2,END=3)
         IF(KFNAME.NE.KSNAME) THEN
            CALL WFRELM
            CALL WFRZON
         ENDIF
         IF(KFNAMZ.NE.KSNAMZ) THEN
            CALL WFRZON
         ENDIF
         IF(KFNAMA.NE.KSNAMA) NAMAX=0
         IF(KFNAMF.NE.KSNAMF) THEN
            IF(NFOPEN.NE.0) THEN
               CLOSE(26)
               NFOPEN=0
            ENDIF
         ENDIF
         KID=' '
         GOTO 4
C
    2    WRITE(6,*) 'XX IO ERROR: IOSTAT=',IST
         CALL WFPLST
         GOTO 1
    3    KID='Q'
    4    CONTINUE
C
      GOTO 3000
C
C     --- process input line including '=' ---
C
      ENTRY WFPARL(KLINE)
C
      MODE=1
      KNAME=' &WF '//KLINE//' &END'
C      KNAME=' &WF '//KLINE//' /'
C
C     --- when internal file does not accept namelist ---
C
      WRITE(7,'(A80)') KNAME
      REWIND(7)
      READ(7,WF,IOSTAT=IST,ERR=8,END=8)
C
C     --- when internal file accepts namelist ---
C
C      READ(KNAME,IN,ERR=8,END=8)
C
      WRITE(6,*) '#### PARM INPUT ACCEPTED.'
      GOTO 9
    8 WRITE(6,*) 'XX IO ERROR: IOSTAT=',IST
      CALL WFPLST
    9 REWIND(7)
      GOTO 3000
C
C     --- process default input file ---
C
      ENTRY WFPARF
C
      MODE=2
      KPNAME='wfparm'
      INQUIRE(FILE=KPNAME,EXIST=LEX)
      IF(.NOT.LEX) RETURN
C
      OPEN(8,FILE=KPNAME,IOSTAT=IST,STATUS='OLD',ERR=9100)
      READ(8,WF,IOSTAT=IST,ERR=9800,END=9900)
      WRITE(6,*) '## FILE (',KPNAME,') IS ASSIGNED FOR PARM INPUT'
C
C     --- check input data are valid---
C
 3000 IERR=0
      IF(NSMAX.GT.NSM) THEN
         WRITE(6,*) '## NSMAX .GT. NSM : NSMAX =',NSMAX,'  NSM = ',NSM
         NSMAX=NSM
         IERR=1
      ENDIF
      IF(NFMAX.GT.NFM) THEN
         WRITE(6,*) '## NFMAX .GT. NFM : NFMAX =',NFMAX,'  NFM = ',NFM
         NFMAX=NFM
         IERR=1
      ENDIF
      IF(NAMAX.GT.NAM) THEN
         WRITE(6,*) '## NAMAX .GT. NAM : NAMAX =',NAMAX,'  NAM = ',NAM
         NAMAX=NAM
         IERR=1
      ENDIF
      IF(NRMAX.GT.NRM) THEN
         WRITE(6,*) '## NRMAX .GT. NRM : NRMAX =',NRMAX,'  NRM = ',NRM
         NRMAX=NRM
         IERR=1
      ENDIF
      IF(NZMAX.GT.NZM) THEN
         WRITE(6,*) '## NZMAX .GT. NZM : NZMAX =',NZMAX,'  NZM = ',NZM
         NZMAX=NZM
         IERR=1
      ENDIF
C
      IF(IERR.NE.0.AND.MODE.EQ.0) GOTO 1000
      RETURN
C
 9100 WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
      RETURN
 9800 WRITE(6,*) 'XX PARM FILE READ ERROR : IOSTAT = ',IST
      RETURN
 9900 WRITE(6,*) 'XX PARM FILE EOF ERROR'
      RETURN
      END
C
C     ****** SHOW PARAMETERS ******
C
      SUBROUTINE WFVIEW
C
      INCLUDE 'wfcomm.inc'
C
      WRITE(6,602) 'MODELS',MODELS,'MODELB',MODELB,
     &             'MODELD',MODELD,'MODELP',MODELP
      WRITE(6,602) 'MODELW',MODELW,'MODELT',MODELT,
     &             'MODELN',MODELN,'MODIFY',MODIFY
      IF(MODELB.LE.2.OR.MODELB.EQ.4) THEN
         WRITE(6,601) 'BB    ',BB    ,'RA    ',RA
         IF(NZMAX.EQ.1) THEN
            WRITE(6,601) 'RF    ',RF    ,'RKZ   ',RKZ   ,
     &                   'PHIW  ',PHIW
         ELSE
            WRITE(6,605) 'RF    ',RF    ,'NZMAX ',NZMAX ,
     &                   'RZ    ',RZ    ,'PHIW  ',PHIW
         ENDIF
      ELSEIF(MODELB.EQ.3) THEN
         WRITE(6,601) 'BB    ',BB    ,'RA    ',RA    ,
     &                'RMIR  ',RMIR  ,'ZBB   ',ZBB
         IF(NZMAX.EQ.1) THEN
            WRITE(6,603) 'RF    ',RF    ,'NPHI  ',NPHI,
     &                   'PHIW  ',PHIW
         ELSE
            WRITE(6,603) 'RF    ',RF    ,'NZMAX ',NZMAX,
     &                   'PHIW  ',PHIW
         ENDIF
      ELSEIF(MODELB.EQ.5) THEN
         WRITE(6,601) 'BB    ',BB    ,'RA    ',RA    ,
     &                'RR    ',RR
         WRITE(6,601) 'Q0    ',Q0    ,'QA    ',QA    ,
     &                'RKAP  ',RKAP  ,'RDEL  ',RDEL
         IF(NZMAX.EQ.1) THEN
            WRITE(6,603) 'RF    ',RF    ,'NPHI  ',NPHI,
     &                   'PHIW  ',PHIW
         ELSE
            WRITE(6,603) 'RF    ',RF    ,'NZMAX ',NZMAX,
     &                   'PHIW  ',PHIW
         ENDIF
      ELSEIF(MODELB.EQ.6) THEN
         WRITE(6,601) 'BB    ',BB    ,'RA    ',RA
         WRITE(6,601) 'H1    ',H1    ,'H2    ',H2    ,
     &                'RRC   ',RRC
         IF(NZMAX.EQ.1) THEN
            WRITE(6,601) 'RF    ',RF    ,'RKZ   ',RKZ   ,
     &                   'PHIW  ',PHIW
         ELSE
            WRITE(6,605) 'RF    ',RF    ,'NZMAX ',NZMAX ,
     &                   'RZ    ',RZ    ,'PHIW  ',PHIW
         ENDIF
      ENDIF
C
      IF(NNOD.GT.0) THEN
         WRITE(6,601) 'DELX  ',DELX  ,'DELY  ',DELY
         WRITE(6,604) 'NNOD  ',NNOD  ,'NELM  ',NELM 
         WRITE(6,604) 'NNODM ',NNODM ,'NELMM ',NELMM
      ENDIF
C
      WRITE(6,*) '      AJ         APH        AWD        APOS       PIN'
      DO NA=1,NAMAX
         WRITE(6,610) NA,AJ(NA),APH(NA),AWD(NA),APOS(NA),PIN
      ENDDO
C
      IF(MODELD.EQ.0) THEN
         WRITE(6,698)
  698    FORMAT(1H ,'NS    PA',9X,'PZ',9X,'PN',9X,'PNS',
     &                         8X,'PZCL')
         DO NS=1,NSMAX
            WRITE(6,610) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PZCL(NS)
         ENDDO
      ELSE
         WRITE(6,699)
  699    FORMAT(1H ,'NS    PA',9X,'PZ',9X,'PN',9X,'PNS',
     &                         8X,'PTPR',7X,'PTPP',7X,'PTS')
         DO NS=1,NSMAX
            WRITE(6,610) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),
     &                   PTPR(NS),PTPP(NS),PTS(NS)
         ENDDO
      ENDIF
      WRITE(6,604) 'NPRINT',NPRINT,'NDRAWD',NDRAWD,
     &             'NDRAWA',NDRAWA,'NRMAX ',NRMAX
      WRITE(6,603) 'DT    ',DT    ,'NTMAX ',NTMAX ,
     &             'DC    ',DC    ,'NTSTEP',NTSTEP
      WRITE(6,601) 'PPN0  ',PPN0  ,'PTN0  ',PTN0  ,
     &             'PNE0  ',PNE0  ,'PTE0  ',PTE0
      WRITE(6,601) 'PTI0  ',PTI0  ,'PNES  ',PNES  ,
     &             'PTES  ',PTES  ,'PTIS  ',PTIS
      WRITE(6,603) 'EPSIMP',EPSIMP,'MAXIMP',MAXIMP,
     &             'EPSSUM',EPSSUM,'MODELT',MODELT
      WRITE(6,601) 'SGIVEN',SGIVEN,'PGIVEN',PGIVEN,
     &             'XGIVEN',XGIVEN,'YGIVEN',YGIVEN
      WRITE(6,601) 'RGIVEN',RGIVEN,'RFES  ',RFES  ,
     &             'PHIES ',PHIES ,'FACIMP',FACIMP
      WRITE(6,*)
      RETURN
C
  601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  602 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  :
     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
  603 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',I7,4X  :
     &        2X,A6,'=',1PE11.3:2X,A6,'=',I7)
  604 FORMAT(1H ,A6,'=',I6     :2X,A6,'=',I6     :
     &        2X,A6,'=',I6     :2X,A6,'=',I6     :
     &        2X,A6,'=',I6)
  605 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',I7,4X  :
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
  610 FORMAT(1H ,I1,7(1PE11.3))
      END
