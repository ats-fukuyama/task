! obparm.f90

MODULE obparm

  PRIVATE
  PUBLIC ob_parm,ob_view

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
         tmax,delt,eps_obt,del_obt,eps_nw, &
         penergy_in,pcangle_in,zeta_in,psipn_in,theta_in,rr_in,zz_in, &
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
             9X,'eps_obt,del_obt,eps_nw,tmax,delt,'/ &
             9X,'mdlobp,mdlobi,mdlobq,mdlobt,mdlobc,mdlobw,mdlobg,mdlobx,max,'/ &
             9X,'penergy_in,pcangle_in,zeta_in,psipn_in,theta_in,'/ &
             9X,'rr_in,zz_in,nrmax_ob,nthmax_ob,nsumax_ob')
  END SUBROUTINE ob_plst

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE ob_chek(ierr)

    USE obcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    ierr=0

    RETURN
  END SUBROUTINE ob_chek
END MODULE obparm

