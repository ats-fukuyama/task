!     $Id$
!
!     ****** INITIALIZE PARAMETERS ******

SUBROUTINE WFINIT
  
  use wfcomm
  implicit none
  integer :: NA,NM

!     *** CONTROL PARAMETERS ***
!
!        MODELI = Definition of imaginary unit
!               '0' : exp(-i omega t)
!               '1' : exp( j omega t)
!        MODELG = Model geometry
!               '2' : 3D cartesian coordinates
!
!        MODELB = Magnetic Field Configuration
!               '0' : xyz : Bx only
!               '1' : xyz : By only
!               '2' : xyz : Bz only
!
!        MODELD = Dielectric Tensor Model
!               '0' : Cold
!
!        MODELP = Density Profile
!               '0' : Flat
!               '1' : Step function with radius RA
!               '2' : Parabolic with radius RA
!               '3' : Decay in -Z
!
!        MODELN = nas-data Process | Waveguide Boundary
!               '0' : Neglect WG, Absorbing boundary and material
!               '1' : Neglect WG and Absorbing boundary
!               '2' : Include WG but Neglect Absorbing boundary and material
!               '3' : Include WG but Neglect Absorbing boundary
!               '4' : Include Absorbing boundary, but neglect WG and material
!               '5' : Include Absorbing boundary and material, but neglect WG
!               '6' : Include WG and Absorbing boundary, but neglect material
!               '7' : Include WG, Absorbing boundary and material
!             'XX0' : WG mode XX, input from Z=ZNDMAX
!             'XX1' : WG mode XX, input from Z=0
!             'XX2' : WG mode XX, input from Z=ZNDMAX, output Z=0
!             'XX3' : WG mode XX, input from Z=0, output Z=ZNDMAX
!
!        MODELS = Sorting Order
!               '0' : X+Y+Z
!               '1' : YZ plane, then X axis
!               '2' : ZX plane, then Y axis
!               '3' : XY plane, then Z axis
!               '4' : X axis, then YZ plane
!               '5' : Y axis, then ZX plane
!               '6' : Z axis, then XY plane
!
!        MODELX = Equation Type
!               '0' : original
!               '1' : Nabla (Nabla dot A) added
!               '2' : nabla dot A = 0 used
!               '3' : '1'+'2'
!
!        MODELA = Absorpation Type
!               '0' : No absorption
!               '1' : Absorption layer in x direction
!               '2' : Absorption layer in y direction
!               '3' : Absorption layer in z direction
!               '4' : Absorption layer in radial (x-y) direction
!        POSRES: position of resonance surface
!        POSABS: position of absorption layer surface
!        EPSABS: typical value of real part of epsilon
!        DLTABS: typical value of imag part of epsilon
!
!            eps=1.D0+EPSABS*(X-POSABS)**2
!                           /((POSRES-POSABS)*(POSRES-X+CI*DLTEPS))

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
  
!        KFNAME: File name of element data
!        KFNAMA: File name of antenna data
!        KFNAMF: File name of field data
!        KFNAMN: File name of nas data
!        KFNAMB: File name of buffer

  KFNAME = 'elm-data'
  KFNAMA = 'ant-data'
  KFNAMF = 'fld-data'
  KFNAMN = 'nas-data'
  KFNAMB = '/tmp/wfx-buffer'

!     *** CONFIGURATION PARAMETERS ***
!
!        BB    : Magnetic field at center                        (T)
!        RA    : Plasma minor radius                             (m)
!        ZPMAX : Plasma maximum z                                (m)
!        ZPMIN : Plasma minimum z                                (m)
!        ZPLEN : Plasma decay lenth                              (m)

  BB     = 1.D-8
  RA     = 0.25D0
  ZPMAX  = 0.20D0
  ZPMIN  = 0.00D0
  ZPLEN  = 0.05D0

!     *** RF PARAMETERS ***
!
!        RF    : Wave frequency                               (MHz)

  RF     = 13.56D0

!     *** ANTENNA PARAMETERS ***
!
!        NAMAX : Number of antennae
!        AJ    : Antenna current density                       (A/m)
!        APH   : Antenna phase                              (degree)
!        AWD   : Antenna width in (z, phi, Z) direction     (degree)
!        APOS  : Antenna position in (z, phi, Z) direction  (degree)
!        APOS  : Antenna position in (z, phi, Z) direction  (degree)

  NAMAX=0
  DO NA=1,NAM
     AJ(NA)   = 1.D0
     APH(NA)  = 0.D0
     AWD(NA)  = 0.D0
     APOS(NA) = 0.D0
  ENDDO

!     *** PLASMA PARAMETERS ***
!
!        NSMAX : Number of particle species
!        PA    : Mass number
!        PZ    : Charge number
!        PN    : Density at center                     (1.0E20/m**3)
!        PNS   : Density on plasma surface             (1.0E20/m**3)
!        PTPR  : Parallel temperature at center                (keV)
!        PTPP  : Perpendicular temperature at center           (keV)
!        PTS   : Temperature on surface                        (keV)
!        PZCL  : Ratio of collision frequency to wave frequency

  NSMAX= 0
  IF(NSMAX.GT.NSM) NSMAX=NSM
  
  PA(1)  = AME/AMP
  PZ(1)  =-1.0D0
  PN(1)  = 0.0002D0
  PNS(1) = 0.0D0
  PTPR(1)= 0.005D0
  PTPP(1)= 0.005D0
  PTS(1) = 0.0D0
  PZCL(1)= 0.02D0
  
  IF(NSM.GE.2) THEN
     PA(2)  = 39.948D0
     PZ(2)  = 1.0D0
     PN(2)  = 0.0002D0
     PNS(2) = 0.0D0
     PTPR(2)= 0.005D0
     PTPP(2)= 0.005D0
     PTS(2) = 0.0D0
     PZCL(2)= 0.02D0
  ENDIF
  
  IF(NSM.GE.3) THEN
     PA(3)  = 1.0D0
     PZ(3)  = 1.0D0
     PN(3)  = 0.0D0
     PNS(3) = 0.0D0
     PTPR(3)= 0.1D0
     PTPP(3)= 0.1D0
     PTS(3) = 0.0D0
     PZCL(3)= 0.001D0
  ENDIF

!     *** MEDIUM PARAMETERS ***

  call wfmed_allocate
  call wfmed_allocate

  NMMAX=0
  DO NM=1,NMM
     EPSDM(NM)=1.D0
     AMUDM(NM)=1.D0
     SIGDM(NM)=0.D0
  ENDDO
  
  NBMAX=0

!     *** OUTPUT PARAMETERS ***
!
!        NPRINT: Print output parameter
!                0 - No output
!              * 1 - Parameter and global field data
!                2 - Local field data
!                3 - Element data
!        NDRAWD: Drawing parameter for elemendt divider
!                0 : Boundary shape
!              * 1 : Element shape
!                2 : Element shape + Element number
!                3 : Element shape + Element number + Node number
!        NDRAWA: Drawing parameter for antenna generater
!               0 : Antenna primary data
!                1 : Antenna secondary data
!              * 2 : Antenna secondary data + Element shape
!        NGRAPH: Drawing parameter for field data
!                0 : text file output
!              * 1 : Contour plot
!                2 : Paint plot
!                3-6 : Bird'seye view from four directions

  NPRINT = 1
  NDRAWD = 1
  NDRAWA = 2
  NGRAPH = 1

!     *** DIVIDER PARAMETERS ***
!
!        BXMIN : Minimum x (m)
!        BXMAX : Maximum x (m)
!        BYMIN : Minimum y (m)
!        BYMAX : Maximum y (m)
!        BZMIN : Minimum z (m)
!        BZMAX : Maximum z (m)
!        RB    : Boundary radius (m)
!        RBAX  : Inner Boundary radius (m)
!        DELX  : Typical element size in x direction (m)
!        DELY  : Typical element size in y direction (m)
!        DELZ  : Typical element size in z direction (m)

  BXMIN  =-0.35D0
  BXMAX  = 0.35D0
  BYMIN  =-0.35D0
  BYMAX  = 0.35D0
  RB     = 0.35D0
  RBAX   = 0.00D0
  BZMIN  = 0.0D0
  BZMAX  = 0.3D0
  DELX   = 0.1D0
  DELY   = 0.1D0
  DELZ   = 0.1D0
  
!     *** ANTENNA SHAPE PARAMETERS ***
!
!        PIN   : Input Power (W)
!        RD    : Antenna radius (m)
!        THETJ1: Start angle of arc antenna (degree)
!        THETJ2: End angle of arc antenna (degree)
!        NJMAX : Number of primary grid points of antenna
!        ZANT  : Antenna position in z axis (m)

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
  
!     PPN0   : Neutral pressure (Pa)  1 Torr = 1 mmHg = 133.322 Pa
!     PTN0   : Initial neutral temperarure (eV)

  PPN0   = 3.0D0
  PTN0   = 0.03D0
      
!     *** GRAPHICS PARAMETER ***

  NGXMAX = 31
  NGYMAX = 31
  NGVMAX = 31
  
  KGINX(0)='EXRZ0.0 EXIZ0.0 EYRZ0.0 EYIZ0.0 EZRZ0.0 EZIZ0.0 '//&
       &   'P1CZ0.0 P2CZ0.0'
  KGINX(1)='EXRZ0.1 EXIZ0.1 EYRZ0.1 EYIZ0.1 EZRZ0.1 EZIZ0.1 '//&
       &   'P1CZ0.1 P2CZ0.1'
  KGINX(2)='EXRZ0.2 EXIZ0.2 EYRZ0.2 EYIZ0.2 EZRZ0.2 EZIZ0.2 '//&
       &   'P1CZ0.2 P2CZ0.2'
  KGINX(3)='EXRZ0.3 EXIZ0.3 EYRZ0.3 EYIZ0.3 EZRZ0.3 EZIZ0.3 '//&
       &   'P1CZ0.3 P2CZ0.3'
  KGINX(4)='AXRZ0.2 AXIZ0.2 AYRZ0.2 AYIZ0.2 AZRZ0.2 AZIZ0.2 '//&
       &   'AFRZ0.2 AFIZ0.2'
  KGINX(5)='EXRX0.0 EXIX0.0 EYRX0.0 EYIX0.0 EZRX0.0 EZIX0.0 '//&
       &   'P1CX0.0 P2CX0.0'
  KGINX(6)='EXRY0.0 EXIY0.0 EYRY0.0 EYIY0.0 EZRY0.0 EZIY0.0 '//&
       &   'P1CY0.0 P2CY0.0'
  KGINX(7)='EXRZ0.25 EXIZ0.25 EYRZ0.25 EYIZ0.25 '//&
       &   'EZRZ0.25 EZIZ0.25 P1CZ0.25 P2CZ0.25'
  KGINX(8)=' '
  KGINX(9)=' '
  
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
  
!     *** DEBUG CONTROL PARAMETER ***

  IDEBUG = 0

!     *** INITIALIZATION PARAMETERS (DO NOT MODIFY) ***
  
  NNMAX  = 0
  NFOPEN = 0
  XNDMIN = 0.D0
  XNDMAX = 0.D0
  
  NDFILE=25
  
  RETURN
END SUBROUTINE WFINIT

!     ***** INPUT PARAMETER LIST *****

SUBROUTINE WFPLST

  use wfcomm

  if (nrank.eq.0) then
     WRITE(6,*) '&WF: BB,RA,RF,ZPMIN,ZPMAX,ZPLEN,AJ,APH,AWD,APOS,'
     WRITE(6,*) '     PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PPN0,PTN0,'
     WRITE(6,*) '     NSMAX,NAMAX,MODELI,'
     WRITE(6,*) '     MODELG,MODELB,MODELD,MODELP,MODELN,MODELS,'
     WRITE(6,*) '     POSRES,POSABS,EPSABS,DLTABS,MODELA,MODELX,'
     WRITE(6,*) '     NPRINT,NDRAWD,NDRAWA,NGRAPH,'
     WRITE(6,*) '     KFNAME,KFNAMA,KFNAMF,KFNAMN,KFNAMB,'
     WRITE(6,*) '     BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX,'
     WRITE(6,*) '     DELX,DELY,DELZ,'
     WRITE(6,*) '     PIN,RD,THETJ1,THETJ2,NJMAX,ZANT,ZWALL,'
     WRITE(6,*) '     THETS1,THETS2,RD1,RD2,'
     WRITE(6,*) '     NGXMAX,NGYMAX,NGVMAX,KGINX,KGINV,IDEBUG'
  end if
  RETURN
END SUBROUTINE WFPLST

!     ****** INPUT PARAMETERS ******

SUBROUTINE WFPARM(KID)
  
  use wfcomm
  implicit none
  integer :: MODE,IST,IERR
  LOGICAL LEX
  CHARACTER KPNAME*32,KLINE*70,KNAME*80,KID*1

  NAMELIST /WF/ BB,RA,RF,ZPMIN,ZPMAX,ZPLEN,AJ,APH,AWD,APOS,&
                PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PPN0,PTN0,&
                NSMAX,NAMAX,MODELI,&
                MODELG,MODELB,MODELD,MODELP,MODELN,MODELS,&
                POSRES,POSABS,EPSABS,DLTABS,MODELA,MODELX,&
                NPRINT,NDRAWD,NDRAWA,NGRAPH,&
                KFNAME,KFNAMA,KFNAMF,KFNAMN,KFNAMB,&
                BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX,&
                DELX,DELY,DELZ,&
                PIN,RD,THETJ1,THETJ2,NJMAX,ZANT,ZWALL,&
                THETS1,THETS2,RD1,RD2,&
                NGXMAX,NGYMAX,NGVMAX,KGINX,KGINV,IDEBUG
  
  MODE=0
1000 continue

1 CONTINUE
  WRITE(6,*) '## INPUT : &WF '
  READ(5,WF,ERR=2,END=3)
  KID=' '
  GOTO 4
  
2 CALL WFPLST
  GOTO 1
3 KID='Q'
4 CONTINUE
  GOTO 3000

  ENTRY WFPARL(KLINE)
  
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
    
  ENTRY WFPARF

  MODE=2
  KPNAME='wfparm'
  INQUIRE(FILE=KPNAME,EXIST=LEX)
  IF(.NOT.LEX) RETURN
  
  OPEN(25,FILE=KPNAME,IOSTAT=IST,STATUS='OLD',ERR=9100)
  READ(25,WF,ERR=9800,END=9900)
  WRITE(6,*) '## FILE (',KPNAME,') IS ASSIGNED FOR PARM INPUT'

3000 IERR=0
  
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
  IF(IERR.NE.0.AND.MODE.EQ.0) GOTO 1000

  RETURN
  
9100 WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
  RETURN
9800 WRITE(6,*) 'XX PARM FILE READ ERROR'
  RETURN
9900 WRITE(6,*) 'XX PARM FILE EOF ERROR'
  RETURN
END SUBROUTINE WFPARM

!     ****** SHOW PARAMETERS ******

SUBROUTINE WFVIEW
  
  use wfcomm
  implicit none
  integer :: NA,NS

  IF(MODELN.EQ.0) THEN
     WRITE(6,601) 'BB    ',BB    ,'RA    ',RA    ,&
                  'RF    ',RF
  ELSE
     WRITE(6,601) 'RF    ',RF
  ENDIF
  
  WRITE(6,'(6X,5A10)') &  
          '     NNMAX','     NEMAX','    NSDMAX',&
          '    NSFMAX','      MLEN'
  WRITE(6,'(A4,2X,5I10)') &
          'USED',NNMAX,NEMAX,NSDMAX,NSFMAX,MLEN
  
  IF(NAMAX.GT.0) THEN
     WRITE(6,*) '*** ANT ***'
     WRITE(6,*) '      AJ         APH        AWD        APOS'
     DO NA=1,NAMAX
        WRITE(6,610) NA,AJ(NA),APH(NA),AWD(NA),APOS(NA)
     ENDDO
  ENDIF
  
  IF(NSMAX.GT.0) THEN
     WRITE(6,*) '*** COLD ***'
     WRITE(6,698)
698     FORMAT(' ','NS    PA',9X,'PZ',9X,'PN',9X,'PNS',&
             &                      8X,'PZCL')
     DO NS=1,NSMAX
        WRITE(6,610) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PZCL(NS)
     ENDDO
  ENDIF
  
  WRITE(6,*) '*** CONTROL ***'
  WRITE(6,604) 'MODELG',MODELG,'MODELB',MODELB,&
               'MODELD',MODELD,'MODELP',MODELP
  WRITE(6,604) 'MODELN',MODELN,'MODELS',MODELS,&
               'MODELX',MODELX,'MODELA',MODELA
  WRITE(6,604) 'NPRINT',NPRINT,'NDRAWD',NDRAWD,&
               'NDRAWA',NDRAWA,'NGRAPH',NGRAPH
  WRITE(6,604) 'NGXMAX',NGXMAX,'NGYMAX',NGYMAX,&
               'NGVMAX',NGVMAX,'MODELI',MODELI
  IF(MODELA.NE.0) THEN
     WRITE(6,601) 'POSRES',POSRES,'POSABS',POSABS,&
                  'EPSABS',EPSABS,'DLTABS',DLTABS
  ENDIF
  WRITE(6,601) 'PPN0  ',PPN0  ,'PTN0  ',PTN0  
  WRITE(6,*)
  RETURN
  
601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3:&
         &        2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
604 FORMAT(' ',A6,'=',I6     :2X,A6,'=',I6     :&
         &        2X,A6,'=',I6     :2X,A6,'=',I6     :&
          &        2X,A6,'=',I6)
610 FORMAT(' ',I1,7(1PE11.3))

END SUBROUTINE WFVIEW
! --------------------------------------------------------
subroutine wfparm_broadcast

  use libmpi
  use libmtx
  use wfcomm
  implicit none

  integer,dimension(20) :: idata
  real(8),dimension(31) :: ddata
  
! ---  broadcast integer data -----

  if(nrank.eq.0) then
     idata(1) =NSMAX
     idata(2) =NAMAX
     idata(3) =MODELI
     idata(4) =MODELG
     idata(5) =MODELB
     idata(6) =MODELD
     idata(7) =MODELP
     idata(8) =MODELN
     idata(9) =MODELS
     idata(10)=MODELA
     idata(11)=MODELX
     idata(12)=NPRINT
     idata(13)=NDRAWD
     idata(14)=NDRAWA
     idata(15)=NGRAPH
     idata(16)=NJMAX
     idata(17)=NGXMAX
     idata(18)=NGYMAX
     idata(19)=NGVMAX
     idata(20)=IDEBUG
  end if
  
  call mtx_broadcast_integer(idata,20)
  
  NSMAX =idata(1)
  NAMAX =idata(2)
  MODELI=idata(3)
  MODELG=idata(4)
  MODELB=idata(5)
  MODELD=idata(6)
  MODELP=idata(7)
  MODELN=idata(8)
  MODELS=idata(9)
  MODELA=idata(10)
  MODELX=idata(11)
  NPRINT=idata(12)
  NDRAWD=idata(13)
  NDRAWA=idata(14)
  NGRAPH=idata(15)
  NJMAX =idata(16)
  NGXMAX=idata(17)
  NGYMAX=idata(18)
  NGVMAX=idata(19)
  IDEBUG=idata(20)

! ----- broadcast real(8) data ------
  if(nrank.eq.0) then
     ddata(1) =BB
     ddata(2) =RA
     ddata(3) =RF
     ddata(4) =ZPMIN
     ddata(5) =ZPMAX
     ddata(6) =ZPLEN
     ddata(7) =POSRES
     ddata(8) =POSABS
     ddata(9) =EPSABS
     ddata(10)=DLTABS
     ddata(11)=BXMIN
     ddata(12)=BXMAX
     ddata(13)=BYMIN
     ddata(14)=BYMAX
     ddata(15)=BZMIN
     ddata(16)=BZMAX
     ddata(17)=DELX
     ddata(18)=DELY
     ddata(19)=DELZ
     ddata(20)=PIN
     ddata(21)=RD
     ddata(22)=THETJ1
     ddata(23)=THETJ2
     ddata(24)=ZANT
     ddata(25)=ZWALL
     ddata(26)=PPN0
     ddata(27)=PTN0
     ddata(28)=THETS1
     ddata(29)=THETS2
     ddata(30)=RD1
     ddata(31)=RD2
  end if

  call mtx_broadcast_real8(ddata,34)
  
  BB    =ddata(1)
  RA    =ddata(2)
  RF    =ddata(3)
  ZPMIN =ddata(4)
  ZPMAX =ddata(5)
  ZPLEN =ddata(6)
  POSRES=ddata(7) 
  POSABS=ddata(8)
  EPSABS=ddata(9)
  DLTABS=ddata(10)
  BXMIN =ddata(11)
  BXMAX =ddata(12)
  BYMIN =ddata(13)
  BYMAX =ddata(14)
  BZMIN =ddata(15)
  BZMAX =ddata(16)
  DELX  =ddata(17)
  DELY  =ddata(18)
  DELZ  =ddata(19)
  PIN   =ddata(20)
  RD    =ddata(21)
  THETJ1=ddata(22)
  THETJ2=ddata(23)
  ZANT  =ddata(24)
  ZWALL =ddata(25)
  PPN0  =ddata(26)
  PTN0  =ddata(27)
  THETS1=ddata(28)
  THETS2=ddata(29)
  RD1   =ddata(30)
  RD2   =ddata(31)

  call mtx_broadcast_real8(AJ  ,8)
  call mtx_broadcast_real8(APH ,8)
  call mtx_broadcast_real8(AWD ,8)
  call mtx_broadcast_real8(APOS,8)
  call mtx_broadcast_real8(PA  ,3)
  call mtx_broadcast_real8(PZ  ,3)
  call mtx_broadcast_real8(PN  ,3)
  call mtx_broadcast_real8(PNS ,3)
  call mtx_broadcast_real8(PZCL,3)
  call mtx_broadcast_real8(PTPR,3)
  call mtx_broadcast_real8(PTPP,3)
  call mtx_broadcast_real8(PTS ,3)

! ------ broadcast character ------

  call mtx_broadcast_character(KFNAME,32)
  call mtx_broadcast_character(KFNAMA,32)
  call mtx_broadcast_character(KFNAMF,32)
  call mtx_broadcast_character(KFNAMN,32)
  call mtx_broadcast_character(KFNAMB,32)

  call mtx_broadcast_character(KGINX(1),80)
  call mtx_broadcast_character(KGINX(2),80)
  call mtx_broadcast_character(KGINX(3),80)
  call mtx_broadcast_character(KGINX(4),80)
  call mtx_broadcast_character(KGINX(5),80)
  call mtx_broadcast_character(KGINX(6),80)
  call mtx_broadcast_character(KGINX(7),80)
  call mtx_broadcast_character(KGINX(8),80)
  call mtx_broadcast_character(KGINX(9),80)
  
  call mtx_broadcast_character(KGINV(1),80)
  call mtx_broadcast_character(KGINV(2),80)
  call mtx_broadcast_character(KGINV(3),80)
  call mtx_broadcast_character(KGINV(4),80)
  call mtx_broadcast_character(KGINV(5),80)
  call mtx_broadcast_character(KGINV(6),80)
  call mtx_broadcast_character(KGINV(7),80)
  call mtx_broadcast_character(KGINV(8),80)
  call mtx_broadcast_character(KGINV(9),80)

  IF(MODELI.EQ.0) THEN
     CII=CI
  ELSE
     CII=-CI
  ENDIF

  return
end subroutine wfparm_broadcast

