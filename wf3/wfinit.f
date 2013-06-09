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
C        CII   : Imaginary unit
C        PI    : Pi
C        AEE   : Elementaty charge
C        AME   : Electron mass
C        AMP   : Atomic unit mass
C        VC    : Speed of light in vacuum
C        AMU0  : Permeability of free space
C        EPS0  : Permittivity of free space
C
      CII    = (0.D0,1.D0)
      PI     = ACOS(0.D0)*2.D0
      AEE    = 1.60217733D-19
      AME    = 9.1093897D-31
      AMP    = 1.6605402D-27
      VC     = 2.997924580D8
      AMU0   = 4.D0*PI*1.D-7
      EPS0   = 1.D0/(VC*VC*AMU0)
C
C     *** CONTROL PARAMETERS ***
C
C        MODELI = Definition of imaginary unit
C               '0' : exp(-i omega t)
C               '1' : exp( j omega t)
C        MODELG = Model geometry
C               '2' : 3D cartesian coordinates
C
C        MODELB = Magnetic Field Configuration
C               '0' : xyz : Bx only
C               '1' : xyz : By only
C               '2' : xyz : Bz only
C
C        MODELD = Dielectric Tensor Model
C               '0' : Cold
C
C        MODELP = Density Profile
C               '0' : Flat
C               '1' : Step function with radius RA
C               '2' : Parabolic with radius RA
C               '3' : Decay in -Z
C
C        MODELN = nas-data Process | Waveguide Boundary
C               '0' : Neglect WG, Absorbing boundary and material
C               '1' : Neglect WG and Absorbing boundary
C               '2' : Include WG but Neglect Absorbing boundary and material
C               '3' : Include WG but Neglect Absorbing boundary
C               '4' : Include Absorbing boundary, but neglect WG and material
C               '5' : Include Absorbing boundary and material, but neglect WG
C               '6' : Include WG and Absorbing boundary, but neglect material
C               '7' : Include WG, Absorbing boundary and material
C             'XX0' : WG mode XX, input from Z=ZNDMAX
C             'XX1' : WG mode XX, input from Z=0
C             'XX2' : WG mode XX, input from Z=ZNDMAX, output Z=0
C             'XX3' : WG mode XX, input from Z=0, output Z=ZNDMAX
C
C        MODELS = Sorting Order
C               '0' : X+Y+Z
C               '1' : YZ plane, then X axis
C               '2' : ZX plane, then Y axis
C               '3' : XY plane, then Z axis
C               '4' : X axis, then YZ plane
C               '5' : Y axis, then ZX plane
C               '6' : Z axis, then XY plane
C
C        MODELX = Equation Type
C               '0' : original
C               '1' : Nabla (Nabla dot A) added
C               '2' : nabla dot A = 0 used
C               '3' : '1'+'2'
C
C        MODELA = Absorpation Type
C               '0' : No absorption
C               '1' : Absorption layer in x direction
C               '2' : Absorption layer in y direction
C               '3' : Absorption layer in z direction
C               '4' : Absorption layer in radial (x-y) direction
C        POSRES: position of resonance surface
C        POSABS: position of absorption layer surface
C        EPSABS: typical value of real part of epsilon
C        DLTABS: typical value of imag part of epsilon
C
C            eps=1.D0+EPSABS*(X-POSABS)**2
C                           /((POSRES-POSABS)*(POSRES-X+CI*DLTEPS))
C
      MODELG=2
      MODELB=0
      MODELD=0
      MODELP=0
      MODELN=0
      MODELS=6
      MODELX=0
      MODELA=0
      POSRES=0.2D0
      POSABS=0.15D0
      EPSABS=1.0D0
      DLTABS=0.1D0
C
C        KFNAME: File name of element data
C        KFNAMA: File name of antenna data
C        KFNAMF: File name of field data
C        KFNAMN: File name of nas data
C        KFNAMB: File name of buffer
C
      KFNAME = 'elm-data'
      KFNAMA = 'ant-data'
      KFNAMF = 'fld-data'
      KFNAMN = 'nas-data'
      KFNAMB = '/tmp/wfx-buffer'
C
C     *** CONFIGURATION PARAMETERS ***
C
C        BB    : Magnetic field at center                        (T)
C        RA    : Plasma minor radius                             (m)
C        ZPMAX : Plasma maximum z                                (m)
C        ZPMIN : Plasma minimum z                                (m)
C        ZPLEN : Plasma decay lenth                              (m)
C
      BB     = 1.D-8
      RA     = 0.25D0
      ZPMAX  = 0.20D0
      ZPMIN  = 0.00D0
      ZPLEN  = 0.05D0
C
C     *** RF PARAMETERS ***
C
C        RF    : Wave frequency                               (MHz)
C
      RF     = 13.56D0
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
      DO NA=1,NAM
         AJ(NA)   = 1.D0
         APH(NA)  = 0.D0
         AWD(NA)  = 0.D0
         APOS(NA) = 0.D0
      ENDDO
C
C     *** PLASMA PARAMETERS ***
C
C        NSMAX : Number of particle species
C        PA    : Mass number
C        PZ    : Charge number
C        PN    : Density at center                     (1.0E20/m**3)
C        PNS   : Density on plasma surface             (1.0E20/m**3)
C        PTPR  : Parallel temperature at center                (keV)
C        PTPP  : Perpendicular temperature at center           (keV)
C        PTS   : Temperature on surface                        (keV)
C        PZCL  : Ratio of collision frequency to wave frequency
C
      NSMAX= 0
      IF(NSMAX.GT.NSM) NSMAX=NSM
C
         PA(1)= AME/AMP
         PZ(1)=-1.0D0
         PN(1)= 0.0002D0
         PNS(1)= 0.D0
         PTPR(1)=0.005D0
         PTPP(1)=0.005D0
         PTS(1)=0.D0
         PZCL(1)= 0.02D0
C
      IF(NSM.GE.2) THEN
         PA(2)= 39.948D0
         PZ(2)= 1.0D0
         PN(2)= 0.0002D0
         PNS(2)= 0.0D0
         PTPR(2)=0.005D0
         PTPP(2)=0.005D0
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
C     *** MEDIUM PARAMETERS ***
C
      NMMAX=0
      DO NM=1,NMM
         EPSDM(NM)=1.D0
         AMUDM(NM)=1.D0
         SIGDM(NM)=0.D0
      ENDDO
C
      NBMAX=0
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
C        NGRAPH: Drawing parameter for field data
C                0 : text file output
C              * 1 : Contour plot
C                2 : Paint plot
C                3-6 : Bird'seye view from four directions
C
      NPRINT = 1
      NDRAWD = 1
      NDRAWA = 2
      NGRAPH = 1
C
C     *** DIVIDER PARAMETERS ***
C
C        BXMIN : Minimum x (m)
C        BXMAX : Maximum x (m)
C        BYMIN : Minimum y (m)
C        BYMAX : Maximum y (m)
C        BZMIN : Minimum z (m)
C        BZMAX : Maximum z (m)
C        RB    : Boundary radius (m)
C        RBAX  : Inner Boundary radius (m)
C        DELX  : Typical element size in x direction (m)
C        DELY  : Typical element size in y direction (m)
C        DELZ  : Typical element size in z direction (m)
C
      BXMIN  =-0.35D0
      BXMAX  = 0.35D0
      BYMIN  =-0.35D0
      BYMAX  = 0.35D0
      RB     = 0.35D0
      RBAX   = 0.00D0
      BZMIN  = 0D0
      BZMAX  = 0.3D0
      DELX   = 0.1D0
      DELY   = 0.1D0
      DELZ   = 0.1D0
C
C     *** ANTENNA SHAPE PARAMETERS ***
C
C        PIN   : Input Power (W)
C        RD    : Antenna radius (m)
C        THETJ1: Start angle of arc antenna (degree)
C        THETJ2: End angle of arc antenna (degree)
C        NJMAX : Number of primary grid points of antenna
C        ZANT  : Antenna position in z axis (m)
C
      PIN    = 1.D0
      RD     = 0.22D0
      THETJ1 = 40.D0
      THETJ2 =-40.D0
      NJMAX  = 41
      THETS1 = 0.D0
      THETS2 = 1440.D0
      RD1    = 0.08D0
      RD2    = 0.22D0
      ZANT   = 0.23D0
      ZWALL  = 0.299999D0
C
C     PPN0   : Neutral pressure (Pa)  1 Torr = 1 mmHg = 133.322 Pa
C     PTN0   : Initial neutral temperarure (eV)
C
      PPN0   = 3.D0
      PTN0   = 0.03D0
C
C     *** GRAPHICS PARAMETER ***
C
      NGXMAX = 31
      NGYMAX = 31
      NGVMAX = 31
C
      KGINX(0)='EXRZ0.0 EXIZ0.0 EYRZ0.0 EYIZ0.0 EZRZ0.0 EZIZ0.0 '//
     &         'P1CZ0.0 P2CZ0.0'
      KGINX(1)='EXRZ0.1 EXIZ0.1 EYRZ0.1 EYIZ0.1 EZRZ0.1 EZIZ0.1 '//
     &         'P1CZ0.1 P2CZ0.1'
      KGINX(2)='EXRZ0.2 EXIZ0.2 EYRZ0.2 EYIZ0.2 EZRZ0.2 EZIZ0.2 '//
     &         'P1CZ0.2 P2CZ0.2'
      KGINX(3)='EXRZ0.3 EXIZ0.3 EYRZ0.3 EYIZ0.3 EZRZ0.3 EZIZ0.3 '//
     &         'P1CZ0.3 P2CZ0.3'
      KGINX(4)='AXRZ0.2 AXIZ0.2 AYRZ0.2 AYIZ0.2 AZRZ0.2 AZIZ0.2 '//
     &         'AFRZ0.2 AFIZ0.2'
      KGINX(5)='EXRX0.0 EXIX0.0 EYRX0.0 EYIX0.0 EZRX0.0 EZIX0.0 '//
     &         'P1CX0.0 P2CX0.0'
      KGINX(6)='EXRY0.0 EXIY0.0 EYRY0.0 EYIY0.0 EZRY0.0 EZIY0.0 '//
     &         'P1CY0.0 P2CY0.0'
      KGINX(7)='EXRZ0.25 EXIZ0.25 EYRZ0.25 EYIZ0.25 '//
     &         'EZRZ0.25 EZIZ0.25 P1CZ0.25 P2CZ0.25'
      KGINX(8)=' '
      KGINX(9)=' '
C
      KGINV(0)=' '
      KGINV(1)=' '
      KGINV(2)=' '
      KGINV(3)=' '
      KGINV(4)=' '
      KGINV(5)=' '
      KGINV(6)=' '
      KGINV(7)=' '
      KGINV(8)=' '
      KGINV(9)=' '
C
C     *** DEBUG CONTROL PARAMETER ***
C
      IDEBUG = 0
C
C     *** INITIALIZATION PARAMETERS (DO NOT MODIFY) ***
C
      NNMAX  = 0
      NFOPEN = 0
      XNDMIN = 0.D0
      XNDMAX = 0.D0
C
      NDFILE=25
C
      RETURN
      END
C
C     ***** INPUT PARAMETER LIST *****
C
      SUBROUTINE WFPLST
C
      WRITE(6,*) '&WF: BB,RA,RF,ZPMIN.ZPMAX,ZPLEN,AJ,APH,AWD,APOS,'
      WRITE(6,*) '     PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PPN0,PTPN0,'
      WRITE(6,*) '     NSMAX,NAMAX,MODELI,'
      WRITE(6,*) '     MODELG,MODELB,MODELD,MODELP,MODELN,MODELS,'
      WRITE(6,*) '     POSRES,POSABS,EPSABS,DLTABS,MODELA,MODELX,'
      WRITE(6,*) '     NPRINT,NDRAWD,NDRAWA,NGRAPH,'
      WRITE(6,*) '     KFNAME,KFNAMA,KFNAMF,KFNAMN,KFNAMB,'
      WRITE(6,*) '     BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX,'
      WRITE(6,*) '     DELX,DELY,DELZ,'
      WRITE(6,*) '     PIN,RD,THETJ1,THETJ2,NJMAX,ZANT,ZWALL,'
      WRITE(6,*) '     THETS1,THETS2,RD1,RD2,'
      WRITE(6,*) '     NGXMAX,NGYMAX,NGZMAX,KGINX,KGINV,IDEBUG'
      RETURN
      END
C
C     ****** INPUT PARAMETERS ******
C
      SUBROUTINE WFPARM(KID)
C
      INCLUDE 'wfcomm.inc'
C
      LOGICAL LEX
      CHARACTER KPNAME*32,KLINE*70,KNAME*80,KID*1
C
      NAMELIST /WF/ BB,RA,RF,ZPMIN,ZPMAX,ZPLEN,
     &              AJ,APH,AWD,APOS,
     &              PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,
     &              NSMAX,NAMAX,MODELI,
     &              MODELG,MODELB,MODELD,MODELP,MODELN,MODELS,
     &              POSRES,POSABS,EPSABS,DLTABS,MODELA,MODELX,
     &              NPRINT,NDRAWD,NDRAWA,NGRAPH,
     &              KFNAME,KFNAMA,KFNAMF,KFNAMN,KFNAMB,
     &              BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX,
     &              DELX,DELY,DELZ,
     &              PIN,RD,THETJ1,THETJ2,NJMAX,ZANT,ZWALL,
     &              PPN0,PTN0,
     &              THETS1,THETS2,RD1,RD2,
     &              NGXMAX,NGYMAX,NGVMAX,KGINX,KGINV,IDEBUG
C
      MODE=0
 1000 CONTINUE
C
    1    CONTINUE
         WRITE(6,*) '## INPUT : &WF '
         READ(5,WF,ERR=2,END=3)
         KID=' '
         GOTO 4
C
    2    CALL WFPLST
         GOTO 1
    3    KID='Q'
    4    CONTINUE
C
      GOTO 3000
C
      ENTRY WFPARL(KLINE)
C
      MODE=1
      KNAME=' &WF '//KLINE//' &END'
      WRITE(33,'(A80)') KNAME
      REWIND(33)
      READ(33,WF,ERR=8,END=8)
      WRITE(6,'(A)') ' ## PARM INPUT ACCEPTED.'
      GOTO 9
    8 CALL WFPLST
    9 REWIND(33)
      GOTO 3000
C
      ENTRY WFPARF
C
      MODE=2
      KPNAME='wfparm'
      INQUIRE(FILE=KPNAME,EXIST=LEX)
      IF(.NOT.LEX) RETURN
C
      OPEN(25,FILE=KPNAME,IOSTAT=IST,STATUS='OLD',ERR=9100)
      READ(25,WF,ERR=9800,END=9900)
      WRITE(6,*) '## FILE (',KPNAME,') IS ASSIGNED FOR PARM INPUT'
C
 3000 IERR=0
      IF(MODELI.EQ.0) THEN
         CI=CII
      ELSE
         CI=-CII
      ENDIF
      IF(NSMAX.GT.NSM) THEN
         WRITE(6,*) '## NSMAX .GT. NSM : NSMAX =',NSMAX,'  NSM = ',NSM
         NSMAX=NSM
         IERR=1
      ENDIF
      IF(NAMAX.GT.NAM) THEN
         WRITE(6,*) '## NAMAX .GT. NAM : NAMAX =',NAMAX,'  NAM = ',NAM
         NAMAX=NAM
         IERR=1
      ENDIF
C
      IF(IERR.NE.0.AND.MODE.EQ.0) GOTO 1000
      RETURN
C
 9100 WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
      RETURN
 9800 WRITE(6,*) 'XX PARM FILE READ ERROR'
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
      IF(MODELN.EQ.0) THEN
         WRITE(6,601) 'BB    ',BB    ,'RA    ',RA    ,
     &                'RF    ',RF
      ELSE
         WRITE(6,601) 'RF    ',RF
      ENDIF
C
      WRITE(6,'(6X,6A10)')  
     &     '     NNMAX','     NEMAX','    NSDMAX',
     &     '    NSFMAX','      MLEN','      MBND'
      WRITE(6,'(A4,2X,6I10)') 
     &     'USED',NNMAX,NEMAX,NSDMAX,NSFMAX,MLEN ,MBND
      WRITE(6,'(A4,2X,6I10)') 
     &     'SIZE',NNM  ,NEM  ,NSDM  ,NSFM  ,MLENM,MBNDM
C
      IF(NAMAX.GT.0) THEN
         WRITE(6,*) '*** ANT ***'
         WRITE(6,*) '      AJ         APH        AWD        APOS'
         DO NA=1,NAMAX
            WRITE(6,610) NA,AJ(NA),APH(NA),AWD(NA),APOS(NA)
         ENDDO
      ENDIF
C
      IF(NSMAX.GT.0) THEN
         WRITE(6,*) '*** COLD ***'
         WRITE(6,698)
  698    FORMAT(' ','NS    PA',9X,'PZ',9X,'PN',9X,'PNS',
     &                      8X,'PZCL')
         DO NS=1,NSMAX
            WRITE(6,610) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PZCL(NS)
         ENDDO
      ENDIF
C
      WRITE(6,*) '*** CONTROL ***'
      WRITE(6,604) 'MODELG',MODELG,'MODELB',MODELB,
     &             'MODELD',MODELD,'MODELP',MODELP
      WRITE(6,604) 'MODELN',MODELN,'MODELS',MODELS,
     &             'MODELX',MODELX,'MODELA',MODELA
      WRITE(6,604) 'NPRINT',NPRINT,'NDRAWD',NDRAWD,
     &             'NDRAWA',NDRAWA,'NGRAPH',NGRAPH
      WRITE(6,604) 'NGXMAX',NGXMAX,'NGYMAX',NGYMAX,
     &             'NGVMAX',NGVMAX,'MODELI',MODELI
      IF(MODELA.NE.0) THEN
      WRITE(6,601) 'POSRES',POSRES,'POSABS',POSABS,
     &             'EPSABS',EPSABS,'DLTABS',DLTABS
      ENDIF
      WRITE(6,601) 'PPN0  ',PPN0  ,'PTN0  ',PTN0  
      WRITE(6,*)
      RETURN
C
  601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3:
     &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
C  602 FORMAT(' ',A6,'=',I7,4X  :2X,A6,'=',I7,4X  :
C     &        2X,A6,'=',I7,4X  :2X,A6,'=',I7)
C  603 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',I7,4X  :
C     &        2X,A6,'=',1PE11.3:2X,A6,'=',I7)
  604 FORMAT(' ',A6,'=',I6     :2X,A6,'=',I6     :
     &        2X,A6,'=',I6     :2X,A6,'=',I6     :
     &        2X,A6,'=',I6)
  610 FORMAT(' ',I1,7(1PE11.3))
C  614 FORMAT(' ',A6,'=',I11    :2X,A6,'=',I11    :
C     &        2X,A6,'=',I11)
      END
