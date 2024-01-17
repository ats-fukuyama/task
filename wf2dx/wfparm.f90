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
           model_config,model_shape, &
           xdiv_min,xdiv_max,ydiv_min,ydiv_max,delx,dely, &
           rdiv_min,rdiv_max,thdiv_min,thdiv_max, &
           RF,RKZ,nph,nant_max,AJ,APH,AWD,APOS,PIN,RD,THETJ1,THETJ2, &
           model_wg,model_wf,&
           xwg_min,xwg_max,ywg_min,ywg_max, &
           phase_wg_min,phase_wg_cen,phase_wg_max, &
           amp_wg,angle_wg,ellip_wg, &
           model_damp,xdamp_min,xdamp_max,ydamp_min,ydamp_max, &
           thdamp_min,thdamp_max,width_damp,factor_damp, &
           model_coll_enhance,factor_coll_enhance, &
           xpos_coll_enhance,xwidth_coll_enhance, &
           ypos_coll_enhance,ywidth_coll_enhance, &
           model_interpolation, &
           nmed_max,model_nmed,epsilon_nmed,amu_nmed,sigma_nmed, &
           xmin_nmed,xmax_nmed,ymin_nmed,ymax_nmed, &
           rmin_nmed,rmax_nmed,thmin_nmed,thmax_nmed, &
           nprint,ngraph,ndrawd,ndraws,ndrawa,ndrawe,ndrawv, &
           sort_weight_x,sort_weight_y, &
           ngxmax,ngymax,ngvmax,gaspect,nxzone_max,nyzone_max, &
           tolerance,modeli, &
           KFNAME,KFNAMA,KFNAMF,KFNAMB,KNAMWG, &
           idebuga

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

    IF(nant_max.GT.nantm) THEN
       WRITE(6,*) '## nant_max .GT. nantm : nant_max =',nant_max, &
            '  nantm = ',nantm
       nant_max=nantm
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
    WRITE(6,*) 'RF,RKZ,NPH,nant_max,AJ,APH,AWD,APOS,model_wg,model_wf,'
    WRITE(6,*) 'PIN,RD,THETJ1,THETJ2,'
    WRITE(6,*) 'xwg_min,xwg_max,ywg_min,ywg_max,phwg_min,phwg_cen,phwg_max,'
    WRITE(6,*) 'amp_wg,angle_wg,ellip_wg,'
    WRITE(6,*) 'model_damp,xdamp_min,xdamp_max,ydamp_min,ydamp_max,'
    WRITE(6,*) 'thdamp_min,thdamp_max,width_damp,factor_damp,'
    WRITE(6,*) 'model_coll_enhance,factor_coll_enhance,'
    WRITE(6,*) 'xpos_coll_enhance,xwidth_coll_enhance,'
    WRITE(6,*) 'ypos_coll_enhance,ywidth_coll_enhance,'
    WRITE(6,*) 'model_interpolation,'
    WRITE(6,*) 'nmed_max,model_nmed,epsilon_nmed,amu_nmed,sigma_nmed,'
    WRITE(6,*) 'xmin_nmed,xmax_nmed,ymin_nmed,ymax_nmed,'
    WRITE(6,*) 'rmin_nmed,rmax_nmed,thmin_nmed,thmax_nmed,'
    WRITE(6,*) 'nprint,ngraph,ndrawd,ndraws,ndrawa,ndrawe,ndrawv,'
    WRITE(6,*) 'sort_weight_x,sort_weight_y,'
    WRITE(6,*) 'delr,delz,bdrmin,bdrmax,bdzmin,bdzmax,'
    WRITE(6,*) 'ngxmax,ngymax,ngvmax,gaspect,nxzone_max,nyzone_max,'
    WRITE(6,*) 'tolerance,idebuga'
    RETURN
  END SUBROUTINE wf_namelist

!     ****** SHOW PARAMETERS ******

  SUBROUTINE wf_view
  
    use wfcomm
    implicit none
    integer :: nant,nmed

    write(6,*) '*** WF2D PARAMETERS ***'

    WRITE(6,'(A10,I6)') 'MODELI=   ',MODELI
    WRITE(6,'(A10,A)')  'KFNAME=   ',TRIM(KFNAME)
    WRITE(6,'(A10,A)')  'KFNAMA=   ',TRIM(KFNAMA)
    WRITE(6,'(A10,A)')  'KFNAMF=   ',TRIM(KFNAMF)
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
       DO nant=1,nant_max
          WRITE(6,610) nant,AJ(nant),APH(nant),AWD(nant),APOS(nant)
       ENDDO
       DO nant=1,nant_max
          WRITE(6,610) nant,PIN(nant),RD(nant),THETJ1(nant),THETJ2(nant)
       ENDDO
    ENDIF

    IF(amp_wg.GT.0.D0) THEN
       WRITE(6,611) 'xwg_min   ',xwg_min, 'xwg_max   ',xwg_max
       WRITE(6,611) 'ywg_min   ',ywg_min, 'ywg_max   ',ywg_max
       WRITE(6,603) 'phase_wg_min',phase_wg_min, &
                    'phase_wg_cen',phase_wg_cen
       WRITE(6,603) 'phase_wg_max',phase_wg_max
       WRITE(6,611) 'amp_wg    ',amp_wg,   'angle_wg ',angle_wg
       WRITE(6,611) 'ellip_wg  ',ellip_wg
       WRITE(6,605) 'model_wg  ',model_wg,'model_wf  ',model_wf
    END IF
  
    WRITE(6,'(A12,I12)')     'model_damp =',model_damp
    WRITE(6,'(A12,ES12.4)')  'xdamp_min  =',xdamp_min
    WRITE(6,'(A12,ES12.4)')  'xdamp_max  =',xdamp_max
    WRITE(6,'(A12,ES12.4)')  'ydamp_min  =',ydamp_min
    WRITE(6,'(A12,ES12.4)')  'ydamp_max  =',ydamp_max
    WRITE(6,'(A12,ES12.4)')  'thdamp_min =',thdamp_min
    WRITE(6,'(A12,ES12.4)')  'thdamp_max =',thdamp_max
    WRITE(6,'(A10,ES12.4)')  'width_damp =',width_damp
    WRITE(6,'(A10,ES12.4)')  'factor_damp=',factor_damp

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

    DO nmed=1,nmed_max
       IF(nmed.EQ.1) THEN
          WRITE(6,'(A)') 'nmed  model_nmed     epsilon         amu       sigma'
          WRITE(6,'(A)') 'nmed        xmin        xmax        ymin        ymax'
          WRITE(6,'(A)') 'nmed        rmin        rmax       thmin       thmax'
       END IF
       WRITE(6,'(I4,I8,4X,4ES12.4)') &
            nmed,model_nmed(nmed),epsilon_nmed(nmed),amu_nmed(nmed), &
            sigma_nmed(nmed)
       WRITE(6,'(I4,4ES12.4)') &
           nmed,xmin_nmed(nmed),xmax_nmed(nmed),ymin_nmed(nmed),ymax_nmed(nmed)
       WRITE(6,'(I4,4ES12.4)') &
         nmed,rmin_nmed(nmed),rmax_nmed(nmed),thmin_nmed(nmed),thmax_nmed(nmed)
    END DO

    
    WRITE(6,*) '***** output control *****'
    
    WRITE(6,604) 'NPRINT',NPRINT,'NGRAPH',NGRAPH,&
                 'NDRAWD',NDRAWD,'NDRAWS',NDRAWS
    WRITE(6,604) 'NDRAWA',NDRAWA,'NDRAWE',NDRAWE,&
                 'NDRAWV',NDRAWV
    WRITE(6,604) 'NGXMAX',NGXMAX,'NGYMAX',NGYMAX,&
                 'NGVMAX',NGVMAX
    WRITE(6,602) 'gaspect   ',gaspect

    WRITE(6,*) '***** control *****'
    
    WRITE(6,602) 'tolerance ',tolerance
    WRITE(6,603) 'sort_weight_x   ',sort_weight_x, &
                 'sort_weight_y   ',sort_weight_y
    WRITE(6,602) 'gaspect   ',gaspect
    WRITE(6,'(A24,I12,4X,A24,I12)') &
         'nxzone_max             =',nxzone_max, &
         'nyzone_max             =',nyzone_max
  RETURN
  
!601 FORMAT(' ',A6,'=',ES12.4:1X,A6,'=',ES12.4:&
!         &  1X,A6,'=',ES12.4:1X,A6,'=',ES12.4)
602 FORMAT(' ',A10,'=',ES12.4:2X,A10,'=',ES12.4:&
         &  2X,A10,'=',ES12.4)
603 FORMAT(' ',A16,'=',1PE12.4:2X,A16,'=',1PE12.4)
604 FORMAT(' ',A6,'=',I6     :2X,A6,'=',I6     :&
          & 2X,A6,'=',I6     :2X,A6,'=',I6     :&
          & 2X,A6,'=',I6)
605 FORMAT(' ',A10,'=',I6    :2X,A10,'=',I6    :&
          & 2X,A10,'=',I6    :2X,A10,'=',I6)
610 FORMAT(' ',I2,7(1PE11.3))
611 FORMAT(' ',A10,'=',ES12.4:1X,A10,'=',ES12.4:&
            1X,A10,'=',ES12.4)

END SUBROUTINE wf_view

! --------------------------------------------------------

SUBROUTINE wfparm_broadcast

  USE libmpi
  USE wfcomm
  IMPLICIT NONE

  INTEGER,DIMENSION(50) :: idata
  REAL(rkind),DIMENSION(100) :: ddata
  
! ---  broadcast integer data -----

  idata(1) =model_config
  idata(2) =model_shape
  idata(3) =nph
  idata(4) =nant_max
  idata(5) =model_wg
  idata(6) =model_wf
  idata(7) =model_damp
  idata(8) =model_coll_enhance
  idata(9) =model_interpolation
  idata(10)=nmed_max
  idata(11)=nprint
  idata(12)=ngraph
  idata(13)=ndrawd
  idata(14)=ndraws
  idata(15)=ndrawa
  idata(16)=ndrawe
  idata(17)=ndrawv
  idata(18)=ngxmax
  idata(19)=ngymax
  idata(20)=ngvmax
  idata(21)=nxzone_max
  idata(22)=nyzone_max
  idata(23)=modeli
  
  call mtx_broadcast_integer(idata,23)
  
  model_config=idata(1)
  model_shape=idata(2)
  nph=idata(3)
  nant_max=idata(4)
  model_wg=idata(5)
  model_wf=idata(6)
  model_damp=idata(7)
  model_coll_enhance=idata(8)
  model_interpolation=idata(9)
  nmed_max=idata(10)
  nprint=idata(11)
  ngraph=idata(12)
  ndrawd=idata(13)
  ndraws=idata(14)
  ndrawa=idata(15)
  ndrawe=idata(16)
  ndrawv=idata(17)
  ngxmax=idata(18)
  ngymax=idata(19)
  ngvmax=idata(20)
  nxzone_max=idata(21)
  nyzone_max=idata(22)
  modeli=idata(23)

! ----- broadcast real(8) data ------

  ddata(1) =xdiv_min
  ddata(2) =xdiv_max
  ddata(3) =ydiv_min
  ddata(4) =ydiv_max
  ddata(5) =delx
  ddata(6) =dely
  ddata(7) =rdiv_min
  ddata(8) =rdiv_max
  ddata(9) =thdiv_min
  ddata(10) =thdiv_max
  ddata(11) =RF
  ddata(12) =RKZ
  ddata(13) =xwg_min
  ddata(14) =xwg_max
  ddata(15) =ywg_min
  ddata(16) =ywg_max
  ddata(17) =phase_wg_min
  ddata(18) =phase_wg_cen
  ddata(19) =phase_wg_max
  ddata(20)=amp_wg
  ddata(21)=angle_wg
  ddata(22)=ellip_wg
  ddata(23)=xdamp_min
  ddata(24)=xdamp_max
  ddata(25)=ydamp_min
  ddata(26)=ydamp_max
  ddata(27)=thdamp_min
  ddata(28)=thdamp_max
  ddata(29)=width_damp
  ddata(30)=factor_damp
  ddata(31)=factor_coll_enhance
  ddata(32)=xpos_coll_enhance
  ddata(33)=xwidth_coll_enhance
  ddata(34)=ypos_coll_enhance
  ddata(35)=ywidth_coll_enhance
  ddata(36)=sort_weight_x
  ddata(37)=sort_weight_y
  ddata(38)=gaspect
  ddata(39)=tolerance

  call mtx_broadcast_real8(ddata,39)
  
  xdiv_min=ddata(1)
  xdiv_max=ddata(2)
  ydiv_min=ddata(3)
  ydiv_max=ddata(4)
  delx=ddata(5)
  dely=ddata(6)
  rdiv_min=ddata(7)
  rdiv_max=ddata(8)
  thdiv_min=ddata(9)
  thdiv_max=ddata(10)
  RF=ddata(11)
  RKZ=ddata(12)
  xwg_min=ddata(13)
  xwg_max=ddata(14)
  ywg_min=ddata(15)
  ywg_max=ddata(16)
  phase_wg_min=ddata(17)
  phase_wg_cen=ddata(18)
  phase_wg_max=ddata(19)
  amp_wg=ddata(20)
  angle_wg=ddata(21)
  ellip_wg=ddata(21)
  xdamp_min=ddata(23)
  xdamp_max=ddata(24)
  ydamp_min=ddata(25)
  ydamp_max=ddata(26)
  thdamp_min=ddata(27)
  thdamp_max=ddata(28)
  width_damp=ddata(29)
  factor_damp=ddata(30)
  factor_coll_enhance=ddata(31)
  xpos_coll_enhance=ddata(32)
  xwidth_coll_enhance=ddata(33)
  ypos_coll_enhance=ddata(34)
  ywidth_coll_enhance=ddata(35)
  sort_weight_x=ddata(36)
  sort_weight_y=ddata(37)
  gaspect=ddata(38)
  tolerance=ddata(39)

  call mtx_broadcast_real8(AJ,nant_max)
  call mtx_broadcast_real8(APH,nant_max)
  call mtx_broadcast_real8(AWD,nant_max)
  call mtx_broadcast_real8(APOS,nant_max)
  call mtx_broadcast_real8(PIN,nant_max)
  call mtx_broadcast_real8(RD,nant_max)
  call mtx_broadcast_real8(THETJ1,nant_max)
  call mtx_broadcast_real8(THETJ2,nant_max)
  
  call mtx_broadcast_integer(model_nmed,nmed_max)
  call mtx_broadcast_real8(epsilon_nmed,nmed_max)
  call mtx_broadcast_real8(amu_nmed,nmed_max)
  call mtx_broadcast_real8(sigma_nmed,nmed_max)
  call mtx_broadcast_real8(xmin_nmed,nmed_max)
  call mtx_broadcast_real8(xmax_nmed,nmed_max)
  call mtx_broadcast_real8(ymin_nmed,nmed_max)
  call mtx_broadcast_real8(ymax_nmed,nmed_max)
  call mtx_broadcast_real8(rmin_nmed,nmed_max)
  call mtx_broadcast_real8(rmax_nmed,nmed_max)
  call mtx_broadcast_real8(thmin_nmed,nmed_max)
  call mtx_broadcast_real8(thmax_nmed,nmed_max)
  
  call mtx_broadcast_integer(idebuga,idebuga_max)

! ------ broadcast character ------

  call mtx_broadcast_character(KFNAME,80)
  call mtx_broadcast_character(KFNAMA,80)
  call mtx_broadcast_character(KFNAMF,80)
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
