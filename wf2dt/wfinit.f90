!     $Id: wfinit.f90,v 1.23 2012/03/05 06:29:02 maruyama Exp $
!
!     ****** INITIALIZE PARAMETERS ******

SUBROUTINE WFINIT
  
  use wfcomm
  implicit none
  integer :: NA,NM

!     *** CONTROL PARAMETERS ***
!
!        MODELI = Definition of imaginary unit
!              * 0 : exp(-i omega t)
!                1 : exp( j omega t)
!
!        MODELG = Model geometry
!                0 : slab geometry
!              * 2 : toroidal geometry
!
!        MODELB = Magnetic Field Configuration
!              * 0 : rzp : tokamak
!
!        MODELD = Dielectric Tensor Model
!              * 0 : Cold
!
!        MODELP = Density Profile
!                0 : Flat
!                1 : Step function with radius RA 
!              * 2 : Parabolic with radius RA
!
!

  CALL pl_allocate_ns

  MODELI=0
  MODELG=2
  MODELB=0
  MODELD=0
  MODELP=2
  
!        KFNAME: File name of element data
!        KFNAMA: File name of antenna data
!        KFNAMF: File name of field data
!        KFNAMB: File name of buffer

  KFNAME = 'elm-data'
  KFNAMA = 'ant-data'
  KFNAMF = 'fld-data'
  KFNAMB = '/tmp/wfx-buffer'

  !     *** INITIAL PARAMETERS FOR DIVIDE ***
  !
  !     NRM,NZM: PARAMETER FOR WFDIV
                                                 
  NRM = 1001
  NZM = 1001

!     *** CONFIGURATION PARAMETERS ***
!
!        BB    : Magnetic field at center                        (T)
!        RR    : Plasma major radius                             (m)
!        RA    : Plasma minor radius                             (m)

  BB     = 0.072D0
  RR     = 0.22D0
  RA     = 0.16D0

!     *** RF PARAMETERS ***
!
!        RF    : Wave frequency                               (MHz)
!        NPH   : Toroidal Mode Number

  RF     = 5000.D0
  NPH    = 0

!     *** ANTENNA PARAMETERS ***
!
!        NAMAX : Number of antennae
!        AJ    : Antenna current density                       (A/m)
!        APH   : Antenna phase                              (degree)
!        AWD   : Antenna width in (z, phi, Z) direction     (degree)
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

  NSMAX= 2
  IF(NSMAX.GT.NSM) NSMAX=NSM
  
  PA(1)  = AME/AMP
  PZ(1)  =-1.0D0
  PN(1)  = 1.0D-3
  PNS(1) = 0.0D0
  PTPR(1)= 1.0D0
  PTPP(1)= 1.0D0
  PTS(1) = 0.0D0
  PZCL(1)= 3.0D-3
  
  IF(NSM.GE.2) THEN
     PA(2)  = 1.D0
     PZ(2)  = 1.0D0
     PN(2)  = 1.D-3
     PNS(2) = 0.0D0
     PTPR(2)= 1.D0
     PTPP(2)= 1.D0
     PTS(2) = 0.0D0
     PZCL(2)= 3.0D-3
  ENDIF
  
  IF(NSM.GE.3) THEN
     PA(3)  = 1.0D0
     PZ(3)  = 1.0D0
     PN(3)  = 0.0D0
     PNS(3) = 0.0D0
     PTPR(3)= 1.0D0
     PTPP(3)= 1.0D0
     PTS(3) = 0.0D0
     PZCL(3)= 3.0D-3
  ENDIF

!     *** MEDIUM PARAMETERS ***

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
!                0 : No output
!              * 1 : Parameter and global field data
!                2 : Local field data
!                3 : Element data
!        NDRAWD: Drawing parameter for elemendt divider
!                0 : Boundary shape
!              * 1 : Element shape
!                2 : Element shape + Element number
!                3 : Element shape + Element number + Node number
!        NDRAWA: Drawing parameter for antenna generater
!                0 : Antenna primary data
!                1 : Antenna secondary data
!              * 2 : Antenna secondary data + Element shape
!        NDRAWE: Drawing parameter for electric field profile
!              * 0 : cylinder (r,phi,z)
!                1 : toroidal (r,theta,phi)
!        NDRAWV: Vector field output parameter
!              * 0 : No output
!                1 : 2D Vector field in poloidal cross section
!        NGRAPH: Drawing parameter for field data
!                0 : text file output
!              * 1 : Contour plot
!                2 : Paint plot
!                3-6 : Bird's eye view from four directions

  NPRINT = 1
  NDRAWD = 1
  NDRAWA = 2
  NDRAWE = 0
  NDRAWV = 0
  NGRAPH = 1

!     *** DIVIDER PARAMETERS ***
!
!        RB    : Boundary radius (m)
!        DELR  : Typical element size in r direction (m)
!        DELZ  : Typical element size in z direction (m)

  RB     = 0.2D0
  DELR   = 0.05D0
  DELZ   = 0.05D0
  BRMIN  = 0.01d0
  BRMAX  = 0.5d0
  BZMIN  =-0.5d0
  BZMAX  = 0.5d0

!     *** ANTENNA SHAPE PARAMETERS ***
!
!        PIN   : Input Power (W)
!        RD    : Antenna radius (m)
!        THETJ1: Start angle of arc antenna (degree)
!        THETJ2: End angle of arc antenna (degree)
!        NJMAX : Number of primary grid points of antenna

  PIN    = 1.D0
  RD     = 0.1855D0
  THETJ1 = 40.D0
  THETJ2 =-40.D0
  NJMAX  = 41
  
!     PPN0   : Neutral pressure (Pa)  1 Torr = 1 mmHg = 133.322 Pa
!     PTN0   : Initial neutral temperarure (eV)

  PPN0   = 3.0D0
  PTN0   = 0.03D0
      
!     *** GRAPHICS PARAMETER ***

  NGXMAX = 31
  NGYMAX = 31
  NGVMAX = 31  
  
!     *** Numerical computation parameter ***

  tolerance = 1.D-8

!     *** DEBUG CONTROL PARAMETER ***

  IDEBUG = 0

!     *** INITIALIZATION PARAMETERS (DO NOT MODIFY) ***
  
  NNMAX  = 0
  NFOPEN = 0
  RNDMIN = 0.D0
  RNDMAX = 0.D0
  
  NDFILE=25
  
  RETURN
END SUBROUTINE WFINIT

!     ***** INPUT PARAMETER LIST *****

SUBROUTINE WFPLST

  use wfcomm

  if (nrank.eq.0) then
     WRITE(6,*) '&WF: BB,RA,RR,RF,AJ,APH,AWD,APOS,NPH,'
     WRITE(6,*) '     PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PPN0,PTN0,'
     WRITE(6,*) '     NSMAX,NAMAX,MODELI,'
     WRITE(6,*) '     MODELG,MODELB,MODELD,MODELP,MODELN,'
     WRITE(6,*) '     NPRINT,NDRAWD,NDRAWA,NDRAWE,NGRAPH,NDRAWV,'
     WRITE(6,*) '     KFNAME,KFNAMA,KFNAMF,KFNAMN,KFNAMB,'
     WRITE(6,*) '     BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX,'
     WRITE(6,*) '     DELR,DELZ,'
     WRITE(6,*) '     PIN,RD,THETJ1,THETJ2,NJMAX,ZANT,ZWALL,'
     WRITE(6,*) '     NGXMAX,NGYMAX,NGVMAX,IDEBUG,'
     WRITE(6,*) '     tolerance'
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

  NAMELIST /WF/ BB,RA,RR,RF,AJ,APH,AWD,APOS,NPH,&
                PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PPN0,PTN0,&
                NSMAX,NAMAX,MODELI,&
                MODELG,MODELB,MODELD,MODELP,MODELN,&
                NPRINT,NDRAWD,NDRAWA,NDRAWE,NGRAPH,NDRAWV,&
                KFNAME,KFNAMA,KFNAMF,KFNAMN,KFNAMB,&
                BRMIN,BRMAX,BZMIN,BZMAX,&
                DELR,DELZ,&
                PIN,RD,THETJ1,THETJ2,NJMAX,&
                NGXMAX,NGYMAX,NGVMAX,IDEBUG, &
                tolerance
  
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
  CLOSE(25)

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

  write(6,*) '*** USED PARAMETERS ***'
  write(6,'(A8,1PE11.3:2X,A7,1PE11.3:2X,A7,1PE11.3)')&
       & ' BB    =',BB    ,'RA    =',RA    ,'RR    ',RR    
  write(6,'(A8,1PE11.3:2X,A7,I11)')&
       & ' RF    =',RF    ,'NPH   =',NPH
  
  write(6,*) '***** DIV *****'
  write(6,'(5A10)') &  
          '     NNMAX','     NEMAX',&
          '    NSDMAX','      MLEN'
  write(6,'(5I10)') &
          NNMAX,NEMAX,NSDMAX,MLEN
  
  IF(NAMAX.GT.0) THEN
     WRITE(6,*) '***** ANT *****'
     WRITE(6,*) '      AJ         APH        AWD        APOS'
     DO NA=1,NAMAX
        WRITE(6,610) NA,AJ(NA),APH(NA),AWD(NA),APOS(NA)
     ENDDO
  ENDIF
  
  IF(NSMAX.GT.0) THEN
     WRITE(6,*) '***** COLD *****'
     WRITE(6,698)
698     FORMAT(' ','NS    PA',9X,'PZ',9X,'PN',9X,'PNS',&
             &                8X,'PZCL')
     DO NS=1,NSMAX
        WRITE(6,610) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PZCL(NS)
     ENDDO
  ENDIF
  
  WRITE(6,*) '***** CONTROL *****'
  WRITE(6,604) 'MODELG',MODELG,'MODELB',MODELB,&
               'MODELD',MODELD,'MODELP',MODELP
  WRITE(6,604) 'MODELN',MODELN,'MODELI',MODELI
  WRITE(6,604) 'NPRINT',NPRINT,'NDRAWD',NDRAWD,&
               'NDRAWA',NDRAWA,'NDRAWE',NDRAWE
  WRITE(6,604) 'NGXMAX',NGXMAX,'NGYMAX',NGYMAX,&
               'NGVMAX',NGVMAX,'NGRAPH',NGRAPH
  WRITE(6,601) 'PPN0  ',PPN0  ,'PTN0  ',PTN0  
  WRITE(6,602) 'tolerance',tolerance
  RETURN
  
601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3:&
         &  2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
602 FORMAT(' ',A10,'=',1PE12.4:2X,A10,'=',1PE12.4:&
         &  2X,A10,'=',1PE12.4)
604 FORMAT(' ',A6,'=',I6     :2X,A6,'=',I6     :&
          & 2X,A6,'=',I6     :2X,A6,'=',I6     :&
          & 2X,A6,'=',I6)
610 FORMAT(' ',I1,7(1PE11.3))

END SUBROUTINE WFVIEW
! --------------------------------------------------------
subroutine wfparm_broadcast

  use libmpi
  use wfcomm
  implicit none

  integer,dimension(19) :: idata
  real(8),dimension(17) :: ddata
  
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
     idata(9) =NPRINT
     idata(10)=NDRAWD
     idata(11)=NDRAWA
     idata(12)=NDRAWE
     idata(13)=NGRAPH
     idata(14)=NJMAX
     idata(15)=NGXMAX
     idata(16)=NGYMAX
     idata(17)=NGVMAX
     idata(18)=IDEBUG
     idata(19)=NPH
  end if
  
  call mtx_broadcast_integer(idata,19)
  
  NSMAX =idata(1)
  NAMAX =idata(2)
  MODELI=idata(3)
  MODELG=idata(4)
  MODELB=idata(5)
  MODELD=idata(6)
  MODELP=idata(7)
  MODELN=idata(8)
  NPRINT=idata(9)
  NDRAWD=idata(10)
  NDRAWA=idata(11)
  NDRAWE=idata(12)
  NGRAPH=idata(13)
  NJMAX =idata(14)
  NGXMAX=idata(15)
  NGYMAX=idata(16)
  NGVMAX=idata(17)
  IDEBUG=idata(18)
  NPH   =idata(19)

! ----- broadcast real(8) data ------

  if(nrank.eq.0) then
     ddata(1) =BB
     ddata(2) =RA
     ddata(3) =RF
     ddata(4) =BRMIN
     ddata(5) =BRMAX
     ddata(6) =BZMIN
     ddata(7) =BZMAX
     ddata(8) =DELR
     ddata(9) =DELZ
     ddata(10)=PIN
     ddata(11)=RD
     ddata(12)=THETJ1
     ddata(13)=THETJ2
     ddata(14)=PPN0
     ddata(15)=PTN0
     ddata(16)=RR
     ddata(17)=tolerance
  end if

  call mtx_broadcast_real8(ddata,17)
  
  BB    =ddata(1)
  RA    =ddata(2)
  RF    =ddata(3)
  BRMIN =ddata(4)
  BRMAX =ddata(5)
  BZMIN =ddata(6)
  BZMAX =ddata(7)
  DELR  =ddata(8)
  DELZ  =ddata(9)
  PIN   =ddata(10)
  RD    =ddata(11)
  THETJ1=ddata(12)
  THETJ2=ddata(13)
  PPN0  =ddata(14)
  PTN0  =ddata(15)
  RR    =ddata(16)
  tolerance =ddata(17)

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

  IF(MODELI.EQ.0) THEN
     CII=CI
  ELSE
     CII=-CI
  ENDIF

  return
end subroutine wfparm_broadcast

