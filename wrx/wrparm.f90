! wrparm.f90

MODULE wrparm

  PRIVATE
  PUBLIC wr_parm

CONTAINS

!     ****** INPUT PARAMETERS ******

  SUBROUTINE WR_PARM(MODE,KIN,IERR)

!     MODE=0 : standard namelist input
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

    USE dpparm,ONLY: dp_chek
    USE dpcomm,ONLY: nsamax_dp
!    USE wrcomm,ONLY: nsamax_wr,nsmax
    USE wrcomm
    USE libkio
    IMPLICIT NONE
    INTEGER,INTENT(IN):: MODE
    CHARACTER(LEN=*),INTENT(IN):: KIN
    INTEGER,INTENT(OUT):: IERR
    EXTERNAL EQCHEK

    IERR=0
1   CALL TASK_PARM(MODE,'WR',KIN,WRNLIN,WRPLST,IERR)
    IF(IERR.NE.0) RETURN

    CALl EQCHEK(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
    IF(NSAMAX_WR.GT.NSMAX) NSAMAX_WR=NSMAX
    NSAMAX_DP=NSAMAX_WR
    CALl DP_CHEK(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1

    RETURN
  END SUBROUTINE WR_PARM

!     ****** INPUT NAMELIST ******

  SUBROUTINE WRNLIN(NID,IST,IERR)

    USE wrcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NID
    INTEGER,INTENT(OUT):: IST,IERR
    INTEGER:: NS

    NAMELIST /WR/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
                  NSMAX,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PZCL, &
                  r_corner,z_corner, &
                  br_corner,bz_corner,bt_corner, &
                  pn_corner,ptpr_corner,ptpp_corner, &
                  PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
                  RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG, &
                  PPN0,PTN0,RF_PL, &
                  MODELG,MODELN,MODELQ,model_coll,MODEL_PROF,MODEL_NPROF, &
                  RHOGMN,RHOGMX, &
                  KNAMEQ,KNAMWR,KNAMFP,KNAMFO,KNAMEQ2, &
                  MODEFW,MODEFR,IDEBUG, &
                  MODELP,MODELV,NCMIN,NCMAX,NS_NSA_DP,PMAX_DP,EMAX_DP, &
                  RFIN,RPIN,ZPIN,PHIIN,RNPHIN,MODEWIN,RNKIN,UUIN, &
                  ANGTIN,ANGPIN,RCURVAIN,RCURVBIN,RBRADAIN,RBRADBIN, &
                  NRAYMAX,NSTPMAX,NRSMAX,NRLMAX,LMAXNW, &
                  NPMAX_DP,NTHMAX_DP,NRMAX_DP,NSAMAX_WR,NS_NSA_WR,nsa_grf, &
                  MDLWRI,MDLWRG,MDLWRP,MDLWRQ,MDLWRW,nres_max,nres_type, &
                  model_fdrv,model_fdrv_ds, &
                  SMAX,DELS,UUMIN,EPSRAY,DELRAY,DELDER,DELKR,EPSNW,EPSD0, &
                  mode_beam,pne_threshold,bdr_threshold, &
                  Rmax_wr,Rmin_wr,Zmax_wr,Zmin_wr,ra_wr,idebug_wr,KNAMWRW

    READ(NID,WR,IOSTAT=IST,ERR=9800,END=9900)
    
    IF(MODEL_PROF.EQ.0) THEN
       DO NS=1,NSMAX
          PROFN1(NS)=PROFN1(1)
          PROFN2(NS)=PROFN2(1)
          PROFT1(NS)=PROFT1(1)
          PROFT2(NS)=PROFT2(1)
          PROFU1(NS)=PROFU1(1)
          PROFU2(NS)=PROFU2(1)
       END DO
    END IF
    IERR=0
      RETURN

 9800 IERR=8
      RETURN
 9900 IERR=9
      RETURN
  END SUBROUTINE WRNLIN

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE WRPLST

    WRITE(6,601)
    RETURN

  601 FORMAT(' ','# &WR : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'/ &
             9X,'NSMAX,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PZCL,'/ &
             9X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'/ &
             9X,'r_corner,z_corner,br_corner,bz_corner,bt_corner,'/ &
             9X,'pn_corner,ptpr_corner,ptpp_corner,'/ &
             9X,'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,'/ &
             9X,'PPN0,PTN0,RFCL,'/ &
             9X,'MODELG,MODELN,MODELQ,model_coll,MODEL_PROF,MODEL_NPROF,'/ &
             9X,'RHOGMN,RHOGMX,'/ &
             9X,'KNAMEQ,KNAMWR,KNAMFP,KNAMFO,KNAMEQ2'/ &
             9X,'MODEFW,MODEFR,IDEBUG'/ &
             9X,'MODELP,MODELV,NCMIN,NCMAX,NS_NSA_DP,PMAX_DP,EMAX_DP,'/ &
             9X,'RFIN,RPIN,ZPIN,PHIIN,RNPHIN,MODEWIN,RNKIN,UUIN,'/ &
             9X,'ANGTIN,ANGPIN,RCURVAIN,RCURVBIN,RBRADAIN,RBRADBIN,'/ &
             9X,'NRAYMAX,NSTPMAX,NRSMAX,NRLMAX,LMAXNW,'/ &
             9X,'NPMAX_DP,NTHMAX_DP,NRMAX_DP,NSAMAX_WR,NS_NSA_WR,nsa_grf'/ &
             9X,'MDLWRI,MDLWRG,MDLWRP,MDLWRQ,MDLWRW,nres_max,nres_type,'/ &
             9X,'model_fdrv,model_fdrv_ds,'/ &
             9X,'SMAX,DELS,UUMIN,EPSRAY,DELRAY,DELDER,DELKR,EPSNW,EPSD0,'/ &
             9X,'mode_beam,pne_threshold,bdr_thershold'/ &
             9X,'Rmax_wr,Rmin_wr,Zmax_wr,Zmin_wr,ra_wr,idebug_wr'/ &
             9X,'KNAMWRW')
  END SUBROUTINE WRPLST

END MODULE wrparm
