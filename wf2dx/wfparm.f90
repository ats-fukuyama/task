! wfparm.f90

MODULE WFPARM

  PRIVATE
  PUBLIC wf_parm
  PUBLIC wf_view
  PUBLIC wfparm_broadcast

CONTAINS

!     ****** INPUT PARAMETERS ******

  SUBROUTINE wf_parm(MODE,KIN,IERR)

!     MODE=0 : standard namelinst input
!     MODE=1 : namelist file input
!     MODE=2 : namelist line input

!     IERR=0 : normal end
!     IERR=1 : namelist standard input error
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE
!     IERR=10X : input parameter out of range

    USE libkio
    IMPLICIT NONE
    INTEGER,INTENT(IN):: MODE
    CHARACTER(LEN=*),INTENT(IN):: KIN
    INTEGER,INTENT(OUT):: IERR

    IERR=0

1   CALL TASK_PARM(MODE,'WF',KIN,wf_nlin,wf_namelist,IERR)
    IF(IERR.NE.0) RETURN

    CALl EQCHEK(IERR)
    CALl wf_check(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1

    RETURN
  END SUBROUTINE wf_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE wf_nlin(NID,IST,IERR)
  
    USE wfcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NID
    INTEGER,INTENT(OUT):: IST,IERR

    NAMELIST /WF/ &
           RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
           RMIR,ZBB,Hpitch1,Hpitch2,RRCH,RCOIL,ZCOIL,BCOIL,NCOILMAX, &
           NSMAX,NPA,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PUPR,PUPP,PNUC,PZCL, &
           ID_NS,KID_NS, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG, &
           PPN0,PTN0,RF_PL,BAXIS_SCALED, &
           r_corner,z_corner, &
           br_corner,bz_corner,bt_corner, &
           pn_corner,ptpr_corner,ptpp_corner, &
           profn_travis_g,profn_travis_h,profn_travis_p,profn_travis_q, &
           profn_travis_w,proft_travis_g,proft_travis_h,proft_travis_p, &
           proft_travis_q,proft_travis_w, &
           MODELG,MODELB,MODELN,MODELQ,model_coll,MODEL_PROF,MODEL_NPROF, &
           RHOGMN,RHOGMX, &
           KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF, &
           MODEFR,MODEFW,IDEBUG,mdlplw, &
           !           
           MODELP,MODELV,NCMIN,NCMAX,MODEL_ES,EPSRT,LMAXRT, &
           NS_NSA_DP,PMAX_dp,EMAX_dp,RHON_MIN,RHON_MAX, &
           NPMAX_DP,NTHMAX_DP,NRMAX_DP,NSAMAX_DP, &
           !
           MODELI,KFNAME,KFNAMA,KFNAMF,KFNAMN,KFNAMB,KNAMWG, &
           RF,RKZ,NPH,nant_max,AJ,APH,AWD,APOS,PIN,RD,THETJ1,THETJ2, &
           modelwg,modelwf,&
           x1wg,y1wg,x2wg,y2wg,ph1WG,ph2wg,ampwg,angwg,elpwg,dphwg, &
           MDAMP,WDAMP,FDAMP,rdamp_min,rdamp_max,zdamp_min,zdamp_max, &
           model_coll_enhance,factor_coll_enhance, &
           xpos_coll_enhance,xwidth_coll_enhance, &
           ypos_coll_enhance,ywidth_coll_enhance, &
           model_interpolation,nmmax,nkmax,epsdm,amudm,sigdm,nmka, &
!           nbmax,nbpmax,KABDY,PHIBDY,RESBDY,PWRBDY,PHABDY, &
!           XGBDY,YGBDY,ZGBDY,XNBDY,YNBDY,ZNBDY,XPBDY,YPBDY,ZPBDY,SZBDY, &
!           NDBDY,NMBDY,NENBP,NDNBP, &
           nprint,ngraph,ndrawd,ndraws,ndrawa,ndrawe,ndrawv, &
           sort_weight_x,sort_weight_y,delr,delz,bdrmin,bdrmax,bdzmin,bdzmax, &
           ngxmax,ngymax,ngvmax,gfactor,nxzone_max,nyzone_max,tolerance,idebuga

    IERR=0
    
    READ(NID,WF,IOSTAT=IST,ERR=9800,END=9900)
    IF(IST.NE.0) GO TO 9100
     
    CALL wf_check(IERR)
    RETURN
  
9100 WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
    RETURN
9800 WRITE(6,*) 'XX PARM FILE READ ERROR'
    RETURN
9900 WRITE(6,*) 'XX PARM FILE EOF ERROR'
    RETURN
  END SUBROUTINE wf_nlin

!     ****** INPUT PARAMETER CHECK ******

  SUBROUTINE wf_check(IERR)
  
    USE wfcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR

    IERR=0
    
    IF(NSMAX.GT.NSM) THEN
       WRITE(6,*) '## NSMAX .GT. NSM : NSMAX =',NSMAX,'  NSM = ',NSM
       NSMAX=NSM
       IERR=1
    ENDIF

    IF(nant_max.GT.NAM) THEN
       WRITE(6,*) '## nant_max .GT. NAM : nant_max =',nant_max,'  NAM = ',NAM
       nant_max=NAM
       IERR=1
    ENDIF

    IF(NCOILMAX.GT.NCOILM) THEN
       WRITE(6,*) '## NCOILMAX .GT. NCOILM : NCOILMAX =',NCOILMAX, &
                                          '  NCOILM = ',NCOILM
       NCOILMAX=NCOILM
       IERR=1
    ENDIF

    RETURN
  END SUBROUTINE wf_check

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE wf_namelist

    use wfcomm

    WRITE(6,*) '&WF : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'
    WRITE(6,*) 'RMIR,ZBB,Hpitch1,Hpitch2,RRCH,RCOI,ZCOIL,BCOIL,NCOILMAX,'
    WRITE(6,*) 'NSMAX,PA,PZ,PN,PNS,PTPR,PTPP,PTS,'
    WRITE(6,*) 'PU,PUS,PUPR,PUPP,PNUC,PZCL,ID_NS,KID_NS,'
    WRITE(6,*) 'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'
    WRITE(6,*) 'r_corner,z_corner,br_corner,bz_corner,bt_corner,'
    WRITE(6,*) 'pn_corner,ptpr_corner,ptpp_corner,'
    WRITE(6,*) 'profn_travis_g,profn_travis_h,profn_travis_p,'
    WRITE(6,*) 'profn_travis_q,profn_travis_w,proft_travis_g,'
    WRITE(6,*) 'proft_travis_h,proft_travis_p,proft_travis_q,'
    WRITE(6,*) 'proft_travis_w,'
    WRITE(6,*) 'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,'
    WRITE(6,*) 'PPN0,PTN0,RFCL,BAXIS_SCALED,'
    WRITE(6,*) 'MODELG,MODELB,MODELN,MODELQ,'
    WRITE(6,*) 'model_coll,MODEL_PROF,MODEL_NPROF,RHOGMN,RHOGMX,'
    WRITE(6,*) 'KNAMEQ,KNAMWR,KNAMFP,KNAMFO,KNAMEQ2'
    WRITE(6,*) 'MODEFW,MODEFR,IDEBUG,mdlplw,'
    
    WRITE(6,*) 'MODELP,MODELV,NCMIN,NCMAX,'
    WRITE(6,*) 'MODEL_ES,EPSRT,LMAXRT,'
    WRITE(6,*) 'NS_NSA_DP,PMAX_dp,EMAX_dp,ROHN_MIN,ROHN_MAX,'
    WRITE(6,*) 'NPMAX_DP,NTHMAX_DP,NRMAX_DP,NSAMAX_DP,'

    WRITE(6,*) 'MODELI,KFNAME,KFNAMA,KFNAMF,KFNAMB,KNAMWG,'
    WRITE(6,*) 'RF,RKZ,NPH,nant_max,AJ,APH,AWD,APOS,modelwg,modelwf,'
    WRITE(6,*) 'PIN,RD,THETJ1,THETJ2,'
    WRITE(6,*) 'x1wg,y1wg,x2wg,y2wg,ph1WG,ph2wg,ampwg,angwg,elpwg,dphwg,'
    WRITE(6,*) 'MDAMP,WDAMP,FDAMP,rdamp_min,rdamp_max,zdamp_min,zdamp_max,'
    WRITE(6,*) 'model_coll_enhance,factor_coll_enhance,'
    WRITE(6,*) 'xpos_coll_enhance,xwidth_coll_enhance,'
    WRITE(6,*) 'ypos_coll_enhance,ywidth_coll_enhance,'
    WRITE(6,*) 'model_interpolation,nmmax,nkmax,epsdm,amudm,sigdm,nmka,'
!    WRITE(6,*) 'nbmax,nbpmax,KABDY,PHIBDY,RESBDY,PWRBDY,PHABDY,'
!    WRITE(6,*) 'XGBDY,YGBDY,ZGBDY,XNBDY,YNBDY,ZNBDY,XPBDY,YPBDY,ZPBDY,'
!    WRITE(6,*) 'SZBDY,NDBDY,NMBDY,NENBP,NDNBP'
    WRITE(6,*) 'nprint,ngraph,ndrawd,ndraws,ndrawa,ndrawe,ndrawv,'
    WRITE(6,*) 'sort_weight_x,sort_weight_y,'
    WRITE(6,*) 'delr,delz,bdrmin,bdrmax,bdzmin,bdzmax,'
    WRITE(6,*) 'ngxmax,ngymax,ngvmax,gfactor,nxzone_max,nyzone_max,'
    WRITE(6,*) 'tolerance,idebuga'
    RETURN
  END SUBROUTINE wf_namelist

!     ****** SHOW PARAMETERS ******

  SUBROUTINE wf_view
  
    use wfcomm
    implicit none
    integer :: NA,NM

    write(6,*) '*** WF2D PARAMETERS ***'

    WRITE(6,'(A10,I6)') 'MODELI=   ',MODELI
    WRITE(6,'(A10,A)')  'KFNAME=   ',TRIM(KFNAME)
    WRITE(6,'(A10,A)')  'KFNAMA=   ',TRIM(KFNAMA)
    WRITE(6,'(A10,A)')  'KFNAMF=   ',TRIM(KFNAMF)
    WRITE(6,'(A10,A)')  'KFNAMN=   ',TRIM(KFNAMN)
    WRITE(6,'(A10,A)')  'KFNAMB=   ',TRIM(KFNAMB)
    WRITE(6,'(A10,A)')  'KNAMWG=   ',TRIM(KNAMWG)

    WRITE(6,'(A10,ES12.4)')  'RF=       ',RF
    WRITE(6,'(A10,ES12.4)')  'RKZ=      ',RKZ
    WRITE(6,'(A10,I12)')     'NPH=      ',NPH

    IF(nant_max.GT.0) THEN
       WRITE(6,*) '***** ANT *****'
       WRITE(6,'(A10,I12)')     'nant_max=    ',nant_max
       WRITE(6,*) '      AJ         APH        AWD        APOS'
       WRITE(6,*) '     PIN          RD     THETJ1      THETJ2'
       DO NA=1,nant_max
          WRITE(6,610) NA,AJ(NA),APH(NA),AWD(NA),APOS(NA)
       ENDDO
       DO NA=1,nant_max
          WRITE(6,610) NA,PIN(NA),RD(NA),THETJ1(NA),THETJ2(NA)
       ENDDO
    ENDIF

    IF(AMPWG.GT.0.D0) THEN
       WRITE(6,601) 'X1WG  ',X1WG  ,'Y1WG  ',Y1WG  , &
                    'X2WG  ',X2WG  ,'X2WG  ',Y2WG
       WRITE(6,601) 'PH1WG ',PH1WG ,'PH2WG ',PH2WG , &
                    'AMPWG ',AMPWG ,'ANGWG ',ANGWG
       WRITE(6,601) 'ELPWG ',ELPWG ,'DPHWG ',DPHWG
       WRITE(6,605) 'MODELWG',MODELWG,'MODELWF',MODELWF
    END IF
  
    WRITE(6,'(A10,I12)')     'MDAMP     ',MDAMP
    WRITE(6,'(A10,ES12.4)')  'WDAMP=    ',WDAMP
    WRITE(6,'(A10,ES12.4)')  'FDAMP=    ',FDAMP
    WRITE(6,'(A10,ES12.4)')  'rdamp_min=',rdamp_min
    WRITE(6,'(A10,ES12.4)')  'rdamp_max=',rdamp_max
    WRITE(6,'(A10,ES12.4)')  'zdamp_min=',zdamp_min
    WRITE(6,'(A10,ES12.4)')  'zdamp_max=',zdamp_max

    WRITE(6,'(A24,I4,8X,A24,ES12.4)') &
            'model_coll_enhance     =',model_coll_enhance, &
            'factor_coll_enhance    =',factor_coll_enhance
    WRITE(6,'(A24,ES12.4,A24,ES12.4)') &
            'xpos_coll_enhance      =',xpos_coll_enhance, &
            'xwidth_coll_enhance    =',xwidth_coll_enhance
    WRITE(6,'(A24,ES12.4,A24,ES12.4)') &
            'ypos_coll_enhance      =',ypos_coll_enhance, &
            'ywidth_coll_enhance    =',ywidth_coll_enhance

    WRITE(6,'(A24,I8)') 'model_interpolation    =',model_interpolation

    DO NM=1,NMMAX
       IF(NM.EQ.1) &
            WRITE(6,'(A)') '  NM       epsdm       amudm       sigdm    nmka'
       WRITE(6,'(I4,3ES12.4,I8)') nm,epsdm(nm),amudm(nm),sigdm(nm),nmka(nm)
    END DO

    
    WRITE(6,*) '***** output control *****'
    
    WRITE(6,604) 'NPRINT',NPRINT,'NGRAPH',NGRAPH,&
                 'NDRAWD',NDRAWD,'NDRAWS',NDRAWS
    WRITE(6,604) 'NDRAWA',NDRAWA,'NDRAWE',NDRAWE,&
                 'NDRAWV',NDRAWV
    WRITE(6,604) 'NGXMAX',NGXMAX,'NGYMAX',NGYMAX,&
                 'NGVMAX',NGVMAX
    WRITE(6,602) 'gfactor   ',gfactor

    WRITE(6,*) '***** control *****'
    
    WRITE(6,602) 'tolerance ',tolerance
    WRITE(6,603) 'sort_weight_x   ',sort_weight_x, &
                 'sort_weight_y   ',sort_weight_y
    WRITE(6,602) 'delr      ',delr,      'delz      ',delz
    WRITE(6,602) 'bdrmin    ',bdrmin,    'bdrmax    ',bdrmax
    WRITE(6,602) 'bdzmin    ',bdzmin,    'bdzmax    ',bdzmax
    WRITE(6,'(A24,I12,4X,A24,I12)') &
         'nxzone_max             =',nxzone_max, &
         'nyzone_max             =',nyzone_max
  RETURN
  
601 FORMAT(' ',A6,'=',ES12.4:1X,A6,'=',ES12.4:&
         &  1X,A6,'=',ES12.4:1X,A6,'=',ES12.4)
602 FORMAT(' ',A10,'=',ES12.4:2X,A10,'=',ES12.4:&
         &  2X,A10,'=',ES12.4)
603 FORMAT(' ',A16,'=',1PE12.4:2X,A16,'=',1PE12.4)
604 FORMAT(' ',A6,'=',I6     :2X,A6,'=',I6     :&
          & 2X,A6,'=',I6     :2X,A6,'=',I6     :&
          & 2X,A6,'=',I6)
605 FORMAT(' ',A10,'=',I6    :2X,A10,'=',I6    :&
          & 2X,A10,'=',I6    :2X,A10,'=',I6)
610 FORMAT(' ',I2,7(1PE11.3))

END SUBROUTINE wf_view

! --------------------------------------------------------

SUBROUTINE wfparm_broadcast

  USE libmpi
  USE wfcomm
  IMPLICIT NONE

  INTEGER,DIMENSION(23) :: idata
  REAL(rkind),DIMENSION(36) :: ddata
  
! ---  broadcast integer data -----

  idata(1) =MODELI
  idata(2) =NPH
  idata(3) =nant_max
  idata(5) =MODELWG
  idata(6) =MODELWF
  idata(7) =MDAMP
  idata(8) =model_coll_enhance
  idata(9) =model_interpolation
  idata(10)=nmmax
  idata(11)=nkmax
  idata(12)=nprint
  idata(13)=ngraph
  idata(14)=ndrawd
  idata(15)=ndraws
  idata(16)=ndrawa
  idata(17)=ndrawe
  idata(18)=ndrawv
  idata(19)=ngxmax
  idata(20)=ngymax
  idata(21)=ngvmax
  idata(22)=nxzone_max
  idata(23)=nyzone_max
  
  call mtx_broadcast_integer(idata,23)
  
  MODELI=idata(1)
  NPH=idata(2)
  nant_max=idata(3)
  MODELWG=idata(5)
  MODELWF=idata(6)
  MDAMP=idata(7)
  model_coll_enhance=idata(8)
  model_interpolation=idata(9)
  nmmax=idata(10)
  nkmax=idata(11)
  nprint=idata(12)
  ngraph=idata(13)
  ndrawd=idata(14)
  ndraws=idata(15)
  ndrawa=idata(16)
  ndrawe=idata(17)
  ndrawv=idata(18)
  ngxmax=idata(19)
  ngymax=idata(20)
  ngvmax=idata(21)
  nxzone_max=idata(22)
  nyzone_max=idata(23)

! ----- broadcast real(8) data ------

  ddata(1) =RF
  ddata(2) =RKZ
  ddata(3) =x1wg
  ddata(4) =y1wg
  ddata(5) =x2wg
  ddata(6) =y2wg
  ddata(7) =ph1wg
  ddata(8) =ph2wg
  ddata(9) =ampwg
  ddata(10)=angwg
  ddata(11)=elpwg
  ddata(12)=dphwg
  ddata(13)=wdamp
  ddata(14)=fdamp
  ddata(15)=rdamp_min
  ddata(16)=rdamp_max
  ddata(17)=zdamp_min
  ddata(18)=zdamp_max
  ddata(19)=factor_coll_enhance
  ddata(20)=xpos_coll_enhance
  ddata(21)=xwidth_coll_enhance
  ddata(22)=ypos_coll_enhance
  ddata(23)=sort_weight_x
  ddata(24)=sort_weight_y
  ddata(25)=delr
  ddata(26)=delz
  ddata(27)=bdrmin
  ddata(28)=bdrmax
  ddata(29)=bdzmin
  ddata(30)=bdzmax
  ddata(31)=gfactor
  ddata(32)=tolerance

  call mtx_broadcast_real8(ddata,32)
  
  RF=ddata(1)
  RKZ=ddata(2)
  x1wg=ddata(3)
  y1wg=ddata(4)
  x2wg=ddata(5)
  y2wg=ddata(6)
  ph1wg=ddata(7)
  ph2wg=ddata(8)
  ampwg=ddata(9)
  angwg=ddata(10)
  elpwg=ddata(11)
  dphwg=ddata(12)
  wdamp=ddata(13)
  fdamp=ddata(14)
  rdamp_min=ddata(15)
  rdamp_max=ddata(16)
  zdamp_min=ddata(17)
  zdamp_max=ddata(18)
  factor_coll_enhance=ddata(19)
  xpos_coll_enhance=ddata(20)
  xwidth_coll_enhance=ddata(21)
  ypos_coll_enhance=ddata(22)
  sort_weight_x=ddata(23)
  sort_weight_y=ddata(24)
  delr=ddata(25)
  delz=ddata(26)
  bdrmin=ddata(27)
  bdrmax=ddata(28)
  bdzmin=ddata(29)
  bdzmax=ddata(30)
  gfactor=ddata(31)
  tolerance=ddata(32)

  call mtx_broadcast_real8(AJ,nant_max)
  call mtx_broadcast_real8(APH,nant_max)
  call mtx_broadcast_real8(AWD,nant_max)
  call mtx_broadcast_real8(APOS,nant_max)
  call mtx_broadcast_real8(PIN,nant_max)
  call mtx_broadcast_real8(RD,nant_max)
  call mtx_broadcast_real8(THETJ1,nant_max)
  call mtx_broadcast_real8(THETJ2,nant_max)
  
  call mtx_broadcast_real8(epsdm,nmmax)
  call mtx_broadcast_real8(amudm,nmmax)
  call mtx_broadcast_real8(sigdm,nmmax)
  
!  call mtx_broadcast_integer(nbpmax,nbmax)
!  call mtx_broadcast_integer(kabdy,nbmax)
!  call mtx_broadcast_real8(phibdy,nbmax)
!  call mtx_broadcast_real8(resbdy,nbmax)
!  call mtx_broadcast_real8(pwrbdy,nbmax)
!  call mtx_broadcast_real8(phabdy,nbmax)
!  call mtx_broadcast_real8(xgbdy,nbmax)
!  call mtx_broadcast_real8(ygbdy,nbmax)
!  call mtx_broadcast_real8(zgbdy,nbmax)
!  call mtx_broadcast_real8(xnbdy,3*nbmax)
!  call mtx_broadcast_real8(ynbdy,3*nbmax)
!  call mtx_broadcast_real8(znbdy,3*nbmax)
!  call mtx_broadcast_real8(xpbdy,nbmax)
!  call mtx_broadcast_real8(ypbdy,nbmax)
!  call mtx_broadcast_real8(zpbdy,nbmax)
!  call mtx_broadcast_real8(szbdy,2*nbmax)
!  call mtx_broadcast_integer(ndbdy,nbmax)
!  call mtx_broadcast_integer(nmbdy,nbmax)
!  DO nb=1,nbmax
!     call mtx_broadcast_integer(nenbp(1,nb),nbpmax)
!     call mtx_broadcast_integer(ndnbp(1,nb),nbpmax)
!  END DO
  
  call mtx_broadcast_integer(idebuga,idebuga_max)

! ------ broadcast character ------

  call mtx_broadcast_character(KFNAME,80)
  call mtx_broadcast_character(KFNAMA,80)
  call mtx_broadcast_character(KFNAMF,80)
  call mtx_broadcast_character(KFNAMN,80)
  call mtx_broadcast_character(KFNAMB,80)

  call mtx_broadcast_character(KNAMWG,80)

  IF(MODELI.EQ.0) THEN
     CII=CI
  ELSE
     CII=-CI
  ENDIF

  RETURN
END SUBROUTINE wfparm_broadcast
END MODULE WFPARM
