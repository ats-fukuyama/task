! obparm.f90

MODULE obparm

  PRIVATE
  PUBLIC ob_parm
  PUBLIC ob_view
  PUBLIC ob_broadcast

CONTAINS

!     ****** INPUT PARAMETERS ******

  SUBROUTINE ob_parm(mode,kin,ierr)

!     mode=0 : standard namelinst input
!     mode=1 : namelist file input
!     mode=2 : namelist line input

!     ierr=0 : normal end
!     ierr=1 : namelist standard input error
!     ierr=2 : namelist file does not exist
!     ierr=3 : namelist file open error
!     ierr=4 : namelist file read error
!     ierr=5 : namelist file abormal end of file
!     ierr=6 : namelist line input error
!     ierr=7 : unknown MODE
!     ierr=10X : input parameter out of range

    USE libkio
    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode
    CHARACTER(LEN=*),INTENT(IN):: kin
    INTEGER,INTENT(OUT):: ierr

    ierr=0

1   CALL task_parm(mode,'OB',kin,ob_nlin,ob_plst,ierr)
    IF(ierr.NE.0) RETURN

    CALl eqchek(ierr)
    IF(ierr.NE.0) GOTO 1
    CALl ob_chek(ierr)
    IF(mode.EQ.0.AND.ierr.NE.0) GOTO 1

    RETURN
  END SUBROUTINE ob_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE ob_nlin(nid,ist,ierr)

    USE obcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nid
    INTEGER,INTENT(OUT):: ist,ierr
    INTEGER:: NS

    NAMELIST /ob/ &
         RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
         RMIR,ZBB,Hpitch1,Hpitch2,RRCH,RCOIL,ZCOIL,BCOIL,NCOILMAX, &
         NSMAX,NPA,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PUPR,PUPP,PZCL, &
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
         MODELG,MODELB,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF, &
         RHOGMN,RHOGMX, &
         KNAMEQ,KNAMWR,KNAMWM,KNAMFP,KNAMFO,KNAMPF, &
         MODEFR,MODEFW,IDEBUG,mdlplw, &
    !
         nobt_max,nstp_max,ns_ob,lmax_nw, &
         mdlobp,mdlobi,mdlobq,mdlobt,mdlobc,mdlobw,mdlobg,mdlobx, &
         tmax_ob,delt_ob,eps_ob,del_ob,eps_nw, &
         penergy_ob_in,pcangle_ob_in,zeta_ob_in,psipn_ob_in, &
         theta_ob_in,rr_ob_in,zz_ob_in, &
         nrmax_ob,nthmax_ob,nsumax_ob

    READ(nid,ob,IOSTAT=ist,ERR=9800,END=9900)
    
    IF(model_prof.EQ.0) THEN
       DO ns=1,nsmax
          profn1(ns)=profn1(1)
          profn2(ns)=profn2(1)
          proft1(ns)=proft1(1)
          proft2(ns)=proft2(1)
          profu1(ns)=profu1(1)
          profu2(ns)=profu2(1)
       END DO
    END IF
    ierr=0
    RETURN

9800 CONTINUE
    ierr=8
    RETURN
9900 CONTINUE
    ierr=9
    RETURN
  END SUBROUTINE ob_nlin

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE ob_plst

    WRITE(6,601)
    RETURN

  601 FORMAT(' ','# &OB : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'/ &
             9X,'RMIR,ZBB,Hpitch1,Hpitch2,RRCH,RCOI,ZCOIL,BCOIL,NCOILMAX,'/ &
             9X,'NSMAX,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PUPR,PUPP,PZCL,'/ &
             9X,'ID_NS,KID_NS,'/ &
             9X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'/ &
             9X,'r_corner,z_corner,br_corner,bz_corner,bt_corner,'/ &
             9X,'pn_corner,ptpr_corner,ptpp_corner,'/ &
             9X,'profn_travis_g,profn_travis_h,profn_travis_p,'/ &
             9X,'profn_travis_q,profn_travis_w,proft_travis_g,'/ &
             9X,'proft_travis_h,proft_travis_p,proft_travis_q,'/ &
             9X,'proft_travis_w,'/ &
             9X,'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,'/ &
             9X,'PPN0,PTN0,RFCL,BAXIS_SCALED,'/ &
             9X,'MODELG,MODELB,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF,'/ &
             9X,'RHOGMN,RHOGMX,'/ &
             9X,'KNAMEQ,KNAMWR,KNAMFP,KNAMFO,KNAMEQ2'/ &
             9X,'MODEFW,MODEFR,IDEBUG,mdlplw,'/ &
             9X,'nobt_max,nstp_max,ns_ob,lmax_nw,'/ &
             9X,'eps_ob,del_ob,eps_nw,tmax_ob,delt_ob,'/ &
             9X,'mdlobp,mdlobi,mdlobq,mdlobt,mdlobc,'/ &
             9X,'mdlobw,mdlobg,mdlobx,max,'/ &
             9X,'penergy_ob_in,spcangle_ob_in,zeta_ob_in,psipn_ob_in,'/&
             9X,'theta_ob_in,rr_ob_in,zz_ob_in,nrmax_ob,nthmax_ob,nsumax_ob')
  END SUBROUTINE ob_plst

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE ob_chek(ierr)

    USE obcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    ierr=0

    RETURN
  END SUBROUTINE ob_chek

  !     ***** BROADCAST INPUT PARAMETERS *****

  SUBROUTINE ob_broadcast

    Use obcomm_parm
    USE libmpi
    USE libmtx
    IMPLICIT NONE
    INTEGER,DIMENSION(99):: idata
    REAL(rkind),DIMENSION(99):: ddata
    
    idata( 1)=nobt_max
    idata( 2)=nstp_max
    idata( 3)=ns_ob
    idata( 4)=lmax_nw
    idata( 5)=mdlobp
    idata( 6)=mdlobi
    idata( 7)=mdlobq
    idata( 8)=mdlobt
    idata( 9)=mdlobc
    idata(10)=mdlobw
    idata(11)=mdlobg
    idata(12)=mdlobx
    idata(13)=nrmax_ob
    idata(14)=nthmax_ob
    idata(15)=nsumax_ob
    CALL mtx_broadcast_integer(idata,15)
    nobt_max=idata( 1)
    nstp_max=idata( 2)
    ns_ob=idata( 3)
    lmax_nw=idata( 4)
    mdlobp=idata( 5)
    mdlobi=idata( 6)
    mdlobq=idata( 7)
    mdlobt=idata( 8)
    mdlobc=idata( 9)
    mdlobw=idata(10)
    mdlobg=idata(11)
    mdlobx=idata(12)
    nrmax_ob=idata(13)
    nthmax_ob=idata(14)
    nsumax_ob=idata(15)

    ddata(1)=tmax_ob
    ddata(2)=delt_ob
    ddata(3)=eps_ob
    ddata(4)=del_ob
    ddata(5)=eps_nw
    CALL mtx_broadcast_real8(ddata, 5)
    tmax_ob=ddata( 1)
    delt_ob=ddata( 2)
    eps_ob=ddata( 3)
    del_ob=ddata( 4)
    eps_nw=ddata( 5)

    CALL mtx_broadcast_real8(penergy_ob_in,nobt_max)
    CALL mtx_broadcast_real8(pcangle_ob_in,nobt_max)
    CALL mtx_broadcast_real8(zeta_ob_in,nobt_max)
    CALL mtx_broadcast_real8(psipn_ob_in,nobt_max)
    CALL mtx_broadcast_real8(theta_ob_in,nobt_max)
    CALL mtx_broadcast_real8(rr_ob_in,nobt_max)
    CALL mtx_broadcast_real8(zz_ob_in,nobt_max)

    RETURN
  END SUBROUTINE ob_broadcast

END MODULE obparm

