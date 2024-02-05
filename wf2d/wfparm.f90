! wfparm.f90

MODULE WFPARM

  PRIVATE
  PUBLIC wf_parm,wf_view,wfparm_broadcast

CONTAINS

!     ****** INPUT PARAMETERS ******

  SUBROUTINE WF_PARM(MODE,KIN,IERR)

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

1   CALL TASK_PARM(MODE,'WF',KIN,WFNLIN,WFPLST,IERR)
    IF(IERR.NE.0) RETURN

    CALl EQCHEK(IERR)
    CALl WFCHEK(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1

    RETURN
  END SUBROUTINE WF_PARM

!     ****** INPUT NAMELIST ******

  SUBROUTINE WFNLIN(NID,IST,IERR)
  
    USE wfcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NID
    INTEGER,INTENT(OUT):: IST,IERR

    NAMELIST /WF/ BB,RA,RB,RR,RF,AJ,APH,AWD,APOS,NPH,RKZ,&
                  PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PPN0,PTN0,&
                  NSMAX,NAMAX,MODELI,&
                  PROFN1,PROFN2,PROFT1,PROFT2,NCMIN,NCMAX, &
                  MODELG,MODELB,MODELD,MODELP,MODELN,&
                  model_coll_enhance,factor_coll_enhance, &
                  xpos_coll_enhance,xwidth_coll_enhance, &
                  ypos_coll_enhance,ywidth_coll_enhance, &
                  NPRINT,NDRAWD,NDRAWA,NDRAWE,NGRAPH,NDRAWV,&
                  KFNAME,KFNAMA,KFNAMF,KFNAMN,KFNAMB,&
                  KNAMPF,KNAMWG, &
                  BDRMIN,BDRMAX,BDZMIN,BDZMAX,&
                  DELR,DELZ,&
                  rmin_div,rmax_div,thmin_div,thmax_div,delr_div,delth_div, &
                  nlayer_max,ch_layer_mode,posl_nlayer,thickness_nlayer, &
                  pos_min,pos_max,step_size, &
                  PIN,RD,THETJ1,THETJ2,NJMAX,&
                  R1WG,Z1WG,R2WG,Z2WG,PH1WG,PH2WG, &
                  AMPWG,ANGWG,ELPWG,DPHWG,MODELWG, &
                  th_wg_min,th_wg_max,gauss_wg, &
                  phase_wg_min,phase_wg_cen,phase_wg_max, &
                  MODELWF, &
                  NGXMAX,NGYMAX,NGVMAX,IDEBUG, &
                  sort_weight_x,sort_weight_y, &
                  br_corner,bz_corner,bt_corner, &
                  pn_corner,ptpr_corner,ptpp_corner, &
                  tolerance,wdamp,fdamp,gfactor, &
                  mdamp,rdamp_min,rdamp_max,zdamp_min,zdamp_max, &
                  thdamp_min,thdamp_max, &
                  NCOILMAX,RCOIL,ZCOIL,BCOIL, &
                  nxzone_max,nyzone_max,idebug_wf
    
    IERR=0
    
    READ(NID,WF,IOSTAT=IST,ERR=9800,END=9900)
    IF(IST.NE.0) GO TO 9100
     
    CALL WFCHEK(IERR)
    RETURN
  
9100 WRITE(6,*) 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
    IERR=1
    RETURN
9800 WRITE(6,*) 'XX PARM FILE READ ERROR'
    IERR=2
    RETURN
9900 WRITE(6,*) 'XX PARM FILE EOF ERROR'
    IERR=3
    RETURN
  END SUBROUTINE WFNLIN

!     ****** INPUT PARAMETER CHECK ******

  SUBROUTINE WFCHEK(IERR)
  
    USE wfcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR

    IERR=0
    
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

    IF(NCOILMAX.GT.NCOILM) THEN
       WRITE(6,*) '## NCOILMAX .GT. NCOILM : NCOILMAX =',NCOILMAX, &
                                          '  NCOILM = ',NCOILM
       NCOILMAX=NCOILM
       IERR=1
    ENDIF

    RETURN
  END SUBROUTINE WFCHEK

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE WFPLST

    use wfcomm

    if (nrank.eq.0) then
       WRITE(6,*) '&WF: BB,RA,RB,RR,RF,AJ,APH,AWD,APOS,NPH,RKZ,'
       WRITE(6,*) '     PA,PZ,PN,PNS,PZCL,PTPR,PTPP,PTS,PPN0,PTN0,'
       WRITE(6,*) '     NSMAX,NAMAX,MODELI,'
       WRITE(6,*) '     PROFN1,PROFN2,PROFT1,PROFT2,NCMIN,NCMAX,'
       WRITE(6,*) '     MODELG,MODELB,MODELD,MODELP,MODELN,'
       WRITE(6,*) '     NPRINT,NDRAWD,NDRAWA,NDRAWE,NGRAPH,NDRAWV,'
       WRITE(6,*) '     model_coll_enhance,factor_coll_enhance,'
       WRITE(6,*) '     xpos_coll_enhance,xwidth_coll_enhance,'
       WRITE(6,*) '     ypos_coll_enhance,ywidth_coll_enhance,'
       WRITE(6,*) '     KFNAME,KFNAMA,KFNAMF,KFNAMN,KFNAMB,'
       WRITE(6,*) '     KNAMPF,KNAMWG,'
       WRITE(6,*) '     BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX,'
       WRITE(6,*) '     DELR,DELZ,delr_div,delth_div,'
       WRITE(6,*) '     rmin_div,rmax_div,thmin_div,thmax_div,'
       WRITE(6,*) '     nlayer_max,ch_layer_mode,posl_nlayer,'
       WRITE(6,*) '     thickness_nlayer,pos_min,pos_max,step_size,'
       WRITE(6,*) '     PIN,RD,THETJ1,THETJ2,NJMAX,ZANT,ZWALL,'
       WRITE(6,*) '     R1WG,Z1WG,R2WG,Z2WG,PH1WG,PH2WG,'
       WRITE(6,*) '     AMPWG,ANGWG,ELPWG,DPHWG,MODELWG,'
       WRITE(6,*) '     th_wg_min,th_wg_max,gauss_wg,'
       WRITE(6,*) '     phase_wg_min,phase_wg_cen,phase_wg_max,'
       WRITE(6,*) '     MODELWF,'
       WRITE(6,*) '     NGXMAX,NGYMAX,NGVMAX,IDEBUG,'
       WRITE(6,*) '     sort_weight_x,sort_weight_y'
       WRITE(6,*) '     br_corner,bz_corner,bt_corner,'
       WRITE(6,*) '     pn_corner,ptpr_corner,ptpp_corner,'
       WRITE(6,*) '     tolerance,wdamp,fdamp,gfactor,'
       WRITE(6,*) '     mdamp,rdamp_min,rdamp_max,zdamp_min,zdamp_max,'
       WRITE(6,*) '     thrdamp_min,thdamp_max,'
       WRITE(6,*) '     NCOILMAX,RCOIL,ZCOIL,BCOIL,'
       WRITE(6,*) '     nxzone_max,nyzone_max,idebug_wf '
    end if
    RETURN
  END SUBROUTINE WFPLST

!     ****** SHOW PARAMETERS ******

  SUBROUTINE WF_VIEW
  
    use wfcomm
    implicit none
    integer :: NA,NS,NC,ID,i
    REAL(rkind):: RZCL(NSMAX)

    write(6,*) '*** USED PARAMETERS ***'
    write(6,'(A8,1PE11.3:2X,A7,1PE11.3:2X,A7,1PE11.3)') &
         ' BB    =',BB    ,'RA    =',RA    ,'RR    ',RR    
    write(6,'(A8,1PE11.3:2X,A7,I11    :2X,A7,1PE11.3)') &
         ' RF    =',RF    ,'NPH   =',NPH   ,'RKZ   ',RKZ
    write(6,'(A8,1PE11.3:2X,A7,1PE11.3:2X,A7,1PE11.3)') &
         ' RB    =',RB
  
    write(6,*) '***** DIV *****'
    write(6,'(5A10)') &  
         '  node_max','  nelm_max', &
         '    NSDMAX','      MLEN'
    write(6,'(5I10)') &
         node_max,nelm_max,NSDMAX,MLEN

    IF(NCOILMAX.GT.0) THEN
       WRITE(6,'(A,I5)') 'NCOILMAX=',NCOILMAX
       WRITE(6,'((A,I5,3(A,1PE12.4)/))') &
            ('  NC=',NC,'  RC=',RCOIL(NC), &
             '  ZC=',ZCOIL(NC),'  BC=',BCOIL(NC), &
             NC=1,NCOILMAX)
    END IF
  
    IF(NAMAX.GT.0) THEN
       WRITE(6,*) '***** ANT *****'
       WRITE(6,*) '       AJ          APH         AWD         APOS'
       DO NA=1,NAMAX
          WRITE(6,610) NA,AJ(NA),APH(NA),AWD(NA),APOS(NA)
       ENDDO
    ENDIF

    IF(AMPWG.GT.0.D0) THEN
       WRITE(6,601) 'R1WG  ',R1WG  ,'Z1WG  ',Z1WG  , &
                    'R2WG  ',R2WG  ,'Z2WG  ',Z2WG
       WRITE(6,601) 'PH1WG ',PH1WG ,'PH2WG ',PH2WG , &
                    'AMPWG ',AMPWG ,'ANGWG ',ANGWG
       WRITE(6,601) 'ELPWG ',ELPWG ,'DPHWG ',DPHWG
       WRITE(6,605) 'MODELWG',MODELWG,'MODELWF',MODELWF
    END IF
  
    IF(NSMAX.GT.0) THEN
       WRITE(6,*) '***** COLD *****'
       WRITE(6,698)
698    FORMAT(' ','NS     PA',10X,'PZ',10X,'PN',10X,'PNS',9X, &
             &           'PZCL',8X,'MODELP')
       SELECT CASE(MODELG)
       CASE(0,12)
          DO NS=1,NSMAX
             WRITE(6,614) NS,PA(NS),PZ(NS),pn_corner(1,NS),pn_corner(2,NS), &
                          PZCL(NS),MODELP(NS)
          ENDDO
       CASE(2)
          DO NS=1,NSMAX
             WRITE(6,614) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PZCL(NS),MODELP(NS)
          ENDDO
       END SELECT

       ID=0
       DO NS=1,NSMAX
          IF(PZCL(NS).EQ.0.D0) ID=1
       END DO
       IF(ID.NE.0) CALL WFCOLL(PN,PTPR,PTPP,RZCL,1)

       WRITE(6,'(A24,I4,8X,A24,ES12.4)') &
            'model_coll_enhance     =',model_coll_enhance, &
            'factor_coll_enhance    =',factor_coll_enhance
       WRITE(6,'(A24,ES12.4,A24,ES12.4)') &
            'xpos_coll_enhance      =',xpos_coll_enhance, &
            'xwidth_coll_enhance    =',xwidth_coll_enhance
       WRITE(6,'(A24,ES12.4,A24,ES12.4)') &
            'ypos_coll_enhance      =',ypos_coll_enhance, &
            'ywidth_coll_enhance    =',ywidth_coll_enhance
    ENDIF

  IF(MODELG.EQ.0) THEN
     WRITE(6,601) 'br_c1 ',br_corner(1),'br_c2 ',br_corner(2), &
                  'br_c3 ',br_corner(3)
     WRITE(6,601) 'bz_c1 ',bz_corner(1),'bz_c2 ',bz_corner(2), &
                  'bz_c3 ',bz_corner(3)
     WRITE(6,601) 'bt_c1 ',bt_corner(1),'bt_c2 ',bt_corner(2), &
                  'bt_c3 ',bt_corner(3)
     DO NS=1,NSMAX
        WRITE(6,611) 'PN  ',ns,pn_corner(1,ns),pn_corner(2,ns), &
                               pn_corner(3,ns)
        WRITE(6,611) 'PTPR',ns,ptpr_corner(1,ns),ptpr_corner(2,ns), &
                               ptpr_corner(3,ns)
        WRITE(6,611) 'PTPP',ns,ptpp_corner(1,ns),ptpp_corner(2,ns), &
                               ptpp_corner(3,ns)
     ENDDO
  END IF
  
  WRITE(6,*) '***** CONTROL *****'
  WRITE(6,604) 'MODELG',MODELG,'MODELB',MODELB,&
               'MODELD',MODELD
  WRITE(6,604) 'MODELN',MODELN,'MODELI',MODELI
  WRITE(6,604) 'NPRINT',NPRINT,'NDRAWD',NDRAWD,&
               'NDRAWA',NDRAWA,'NDRAWE',NDRAWE
  WRITE(6,604) 'NGXMAX',NGXMAX,'NGYMAX',NGYMAX,&
               'NGVMAX',NGVMAX,'NGRAPH',NGRAPH
  WRITE(6,601) 'PPN0  ',PPN0  ,'PTN0  ',PTN0  
  WRITE(6,602) 'tolerance ',tolerance,'wdamp     ',wdamp, &
               'fdamp     ',fdamp
  WRITE(6,603) 'sort_weight_x   ',sort_weight_x, &
               'sort_weight_y   ',sort_weight_y
  WRITE(6,604) 'mdamp ',mdamp
  WRITE(6,602) 'rdamp_min ',rdamp_min, 'rdamp_max ',rdamp_max
  WRITE(6,602) 'zdamp_min ',zdamp_min, 'zdamp_max ',zdamp_max
  WRITE(6,602) 'thdamp_min',thdamp_min,'thdamp_max',thdamp_max
  WRITE(6,605) 'nxzone_max',nxzone_max,'nyzone_max',nyzone_max

  DO i=1,100
     IF(idebug_wf(I).NE.0) &
          WRITE(6,'(A,I3,A,I4)') 'idebug_wf(',i,')=',idebug_wf(i)
  END DO
  RETURN
  
601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3:&
         &  2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
602 FORMAT(' ',A10,'=',1PE12.4:2X,A10,'=',1PE12.4:&
         &  2X,A10,'=',1PE12.4)
603 FORMAT(' ',A16,'=',1PE12.4:2X,A16,'=',1PE12.4)
604 FORMAT(' ',A6,'=',I6     :2X,A6,'=',I6     :&
          & 2X,A6,'=',I6     :2X,A6,'=',I6     :&
          & 2X,A6,'=',I6)
605 FORMAT(' ',A10,'=',I6    :2X,A10,'=',I6    :&
          & 2X,A10,'=',I6    :2X,A10,'=',I6)
610 FORMAT(' ',I2,6(ES12.4))
611 FORMAT(' ',A4,I2,3(ES12.4))
614 FORMAT(' ',I2,5(ES12.4),I6)

END SUBROUTINE WF_VIEW

! --------------------------------------------------------

SUBROUTINE wfparm_broadcast

  USE libmpi
  USE wfcomm
  IMPLICIT NONE

  INTEGER,DIMENSION(24) :: idata
  REAL(rkind),DIMENSION(58) :: ddata
  
! ---  broadcast integer data -----

  if(nrank.eq.0) then
     idata(1) =NSMAX
     idata(2) =NAMAX
     idata(3) =MODELI
     idata(4) =MODELG
     idata(5) =MODELB
     idata(6) =MODELD
     idata(7) =nlayer_max
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
     idata(20)=MODELWG
     idata(21)=MDAMP
     idata(22)=NCOILMAX
     idata(23)=MODELWF
     idata(24)=model_coll_enhance
  end if
  
  call mtx_broadcast_integer(idata,24)
  
  NSMAX =idata(1)
  NAMAX =idata(2)
  MODELI=idata(3)
  MODELG=idata(4)
  MODELB=idata(5)
  MODELD=idata(6)
  nlayer_max=idata(7)
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
  MODELWG =idata(20)
  MDAMP =idata(21)
  NCOILMAX =idata(22)
  MODELWF =idata(23)
  model_coll_enhance=idata(24)

! ----- broadcast real(8) data ------

  if(nrank.eq.0) then
     ddata(1) =BB
     ddata(2) =RA
     ddata(3) =RF
     ddata(4) =BDRMIN
     ddata(5) =BDRMAX
     ddata(6) =BDZMIN
     ddata(7) =BDZMAX
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
     ddata(18)=wdamp
     ddata(19)=fdamp
     ddata(20)=r1wg
     ddata(21)=z1wg
     ddata(22)=r2wg
     ddata(23)=z2wg
     ddata(24)=ph1wg
     ddata(25)=ph2wg
     ddata(26)=ampwg
     ddata(27)=angwg
     ddata(28)=elpwg
     ddata(29)=dphwg
     ddata(30)=gfactor
     ddata(31)=rdamp_min
     ddata(32)=rdamp_max
     ddata(33)=zdamp_min
     ddata(34)=zdamp_max
     ddata(35)=rkz
     ddata(36)=factor_coll_enhance
     ddata(37)=xpos_coll_enhance
     ddata(38)=xwidth_coll_enhance
     ddata(39)=ypos_coll_enhance
     ddata(40)=ywidth_coll_enhance
     ddata(41)=thdamp_min
     ddata(42)=thdamp_max
     ddata(43)=RB
     ddata(44)=th_wg_min
     ddata(45)=th_wg_max
     ddata(46)=phase_wg_min
     ddata(47)=phase_wg_cen
     ddata(48)=phase_wg_max
     ddata(49)=gauss_wg
     ddata(50)=rmin_div
     ddata(51)=rmax_div
     ddata(52)=thmin_div
     ddata(53)=thmax_div
     ddata(54)=delr_div
     ddata(55)=delth_div
     ddata(56)=pos_min
     ddata(57)=pos_max
     ddata(58)=step_size
  end if

  call mtx_broadcast_real8(ddata,58)
  
  BB    =ddata(1)
  RA    =ddata(2)
  RF    =ddata(3)
  BDRMIN=ddata(4)
  BDRMAX=ddata(5)
  BDZMIN=ddata(6)
  BDZMAX=ddata(7)
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
  wdamp =ddata(18)
  fdamp =ddata(19)
  r1wg  =ddata(20)
  z1wg  =ddata(21)
  r2wg  =ddata(22)
  z2wg  =ddata(23)
  ph1wg =ddata(24)
  ph2wg =ddata(25)
  ampwg =ddata(26)
  angwg =ddata(27)
  elpwg =ddata(28)
  dphwg =ddata(29)
  gfactor=ddata(30)
  rdamp_min=ddata(31)
  rdamp_max=ddata(32)
  zdamp_min=ddata(33)
  zdamp_max=ddata(34)
  rkz   =ddata(35)
  factor_coll_enhance=ddata(36)
  xpos_coll_enhance=ddata(37)
  xwidth_coll_enhance=ddata(38)
  ypos_coll_enhance=ddata(39)
  ywidth_coll_enhance=ddata(40)
  thdamp_min=ddata(41)
  thdamp_max=ddata(42)
  RB=ddata(43)
  th_wg_min=ddata(44)
  th_wg_max=ddata(45)
  phase_wg_min=ddata(46)
  phase_wg_cen=ddata(47)
  phase_wg_max=ddata(48)
  gauss_wg=ddata(49)
  rmin_div=ddata(50)
  rmax_div=ddata(51)
  thmin_div=ddata(52)
  thmax_div=ddata(53)
  delr_div=ddata(54)
  delth_div=ddata(55)
  pos_min=ddata(56)
  pos_max=ddata(57)
  step_size=ddata(58)

  call mtx_broadcast_integer(MODELP,nsmax)
  call mtx_broadcast_real8(AJ  ,8)
  call mtx_broadcast_real8(APH ,8)
  call mtx_broadcast_real8(AWD ,8)
  call mtx_broadcast_real8(APOS,8)
  call mtx_broadcast_real8(PA  ,nsmax)
  call mtx_broadcast_real8(PZ  ,nsmax)
  call mtx_broadcast_real8(PN  ,nsmax)
  call mtx_broadcast_real8(PNS ,nsmax)
  call mtx_broadcast_real8(PZCL,nsmax)
  call mtx_broadcast_real8(PTPR,nsmax)
  call mtx_broadcast_real8(PTPP,nsmax)
  call mtx_broadcast_real8(PTS ,nsmax)
  call mtx_broadcast_real8(br_corner,3)
  call mtx_broadcast_real8(bz_corner,3)
  call mtx_broadcast_real8(bt_corner,3)
  call mtx_broadcast_real8(pn_corner,3*nsmax)
  call mtx_broadcast_real8(ptpr_corner,3*nsmax)
  call mtx_broadcast_real8(ptpp_corner,3*nsmax)
  call mtx_broadcast_real8(RCOIL,NCOILMAX)
  call mtx_broadcast_real8(ZCOIL,NCOILMAX)
  call mtx_broadcast_real8(BCOIL,NCOILMAX)
  call mtx_broadcast_real8(posl_nlayer,NLM+1)
  call mtx_broadcast_real8(thickness_nlayer,NLM)

! ------ broadcast character ------

  call mtx_broadcast_character(KFNAME,32)
  call mtx_broadcast_character(KFNAMA,32)
  call mtx_broadcast_character(KFNAMF,32)
  call mtx_broadcast_character(KFNAMN,32)
  call mtx_broadcast_character(KFNAMB,32)

  call mtx_broadcast_character(KNAMPF,80)
  call mtx_broadcast_character(KNAMWG,80)
  call mtx_broadcast_character(ch_layer_mode,1)

  IF(MODELI.EQ.0) THEN
     CII=CI
  ELSE
     CII=-CI
  ENDIF

  RETURN
END SUBROUTINE wfparm_broadcast
END MODULE WFPARM
