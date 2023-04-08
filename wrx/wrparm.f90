! wrparm.f90

MODULE wrparm

  PRIVATE
  PUBLIC wr_parm
  PUBLIC wr_chek

CONTAINS

!     ****** INPUT PARAMETERS ******

  SUBROUTINE WR_PARM(MODE,KIN,IERR)

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

    USE dpparm,ONLY: dp_chek
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
    CALl DP_CHEK(IERR)
    CALl WR_CHEK(IERR)
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
                  NPMAX_DP,NTHMAX_DP,NRMAX_DP, &
                  MDLWRI,MDLWRG,MDLWRP,MDLWRQ,MDLWRW,nres_max,nres_type, &
                  MDLWRF, &
                  SMAX,DELS,UUMIN,EPSRAY,DELRAY,DELDER,DELKR,EPSNW, &
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
             9X,'NPMAX_DP,NTHMAX_DP,NRMAX_DP,'/ &
             9X,'MDLWRI,MDLWRG,MDLWRP,MDLWRQ,MDLWRW,nres_max,nres_type,'/ &
             9X,'MDLWRF,'/ &
             9X,'SMAX,DELS,UUMIN,EPSRAY,DELRAY,DELDER,DELKR,EPSNW'/ &
             9X,'mode_beam,pne_threshold,bdr_thershold'/ &
             9X,'Rmax_wr,Rmin_wr,Zmax_wr,Zmin_wr,ra_wr,idebug_wr'/ &
             9X,'KNAMWRW')
  END SUBROUTINE WRPLST

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE WR_CHEK(IERR)

    USE wrcomm_parm
    USE dpparm,ONLY: dpprep_local
    USE equnit
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    CHARACTER(LEN=80):: LINE
    INTEGER,SAVE:: INITEQ=0
    INTEGER:: nsu
    REAL(rkind):: dth
    EXTERNAL EQCALQ,EQGETB  !,eqget_rzsu
    INTERFACE
       SUBROUTINE eqget_rzsu(rsu,zsu,nsumax)
         USE task_kinds,ONLY: dp
         REAL(dp),ALLOCATABLE,INTENT(OUT):: rsu(:),zsu(:)
         INTEGER,INTENT(OUT):: nsumax
       END SUBROUTINE eqget_rzsu
    END INTERFACE

    IERR=0

    IF(MODELG.EQ.3.OR.MODELG.EQ.5) THEN ! task/eq and eqdsk
       IF(INITEQ.EQ.0) THEN
          CALL eq_load(MODELG,KNAMEQ,IERR)
          IF(IERR.EQ.0) THEN
             WRITE(LINE,'(A,I5)') 'NRMAX =',51
             CALL eq_parm(2,LINE,IERR)
             WRITE(LINE,'(A,I5)') 'NTHMAX=',64
             CALL eq_parm(2,LINE,IERR)
             WRITE(LINE,'(A,I5)') 'NSUMAX=',64
             CALL eq_parm(2,LINE,IERR)
             CALL EQCALQ(IERR)
             CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
             WRITE(6,*) ALLOCATED(rsu_wr)
             CALL eqget_rzsu(rsu_wr,zsu_wr,nsumax)
             rmax_eq=rsu_wr(1)
             rmin_eq=rsu_wr(1)
             zmax_eq=zsu_wr(1)
             zmin_eq=zsu_wr(1)
             DO nsu=2,nsumax
                IF(rsu_wr(nsu).GT.rmax_eq) rmax_eq=rsu_wr(nsu)
                IF(rsu_wr(nsu).LT.rmin_eq) rmin_eq=rsu_wr(nsu)
                IF(zsu_wr(nsu).GT.zmax_eq) zmax_eq=zsu_wr(nsu)
                IF(zsu_wr(nsu).LT.zmin_eq) zmin_eq=zsu_wr(nsu)
             END DO
             raxis_eq=0.5D0*(rmin_eq+rmax_eq)
             zaxis_eq=0.5D0*(zmin_eq+zmax_eq)
             rmax_wr=raxis_eq+bdr_threshold*(rmax_eq-raxis_eq)
             rmin_wr=raxis_eq+bdr_threshold*(rmin_eq-raxis_eq)
             zmax_wr=zaxis_eq+bdr_threshold*(zmax_eq-zaxis_eq)
             zmin_wr=zaxis_eq+bdr_threshold*(zmin_eq-zaxis_eq)
             IF(rmin_wr.LT.0.D0) rmin_wr=0.D0

             WRITE(6,'(A,4ES12.4)') '_eq:',rmin_eq,rmax_eq,zmin_eq,zmax_eq
             WRITE(6,'(A,2ES12.4)') '_ax:',raxis_eq,zaxis_eq
             WRITE(6,'(A,4ES12.4)') '_wr:',rmin_wr,rmax_wr,zmin_wr,zmax_wr
             INITEQ=1
          ELSE
             WRITE(6,*) 'XX EQLOAD: IERR=',IERR
             INITEQ=0
          ENDIF
       ENDIF
    ELSE IF(MODELG.EQ.8) THEN ! task/equ
       IF(INITEQ.EQ.0) THEN
          CALL eq_read(IERR)
          IF(IERR.EQ.0) THEN
             CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
             CALL eqget_rzsu(rsu_wr,zsu_wr,nsumax)
             rmax_eq=rsu_wr(1)
             rmin_eq=rsu_wr(1)
             zmax_eq=zsu_wr(1)
             zmin_eq=zsu_wr(1)
             DO nsu=2,nsumax
                IF(rsu_wr(nsu).GT.rmax_eq) rmax_eq=rsu_wr(nsu)
                IF(rsu_wr(nsu).LT.rmin_eq) rmin_eq=rsu_wr(nsu)
                IF(zsu_wr(nsu).GT.zmax_eq) zmax_eq=rsu_wr(nsu)
                IF(zsu_wr(nsu).LT.zmin_eq) zmin_eq=rsu_wr(nsu)
             END DO
             raxis_eq=0.5D0*(rmin_eq+rmax_eq)
             zaxis_eq=0.5D0*(zmin_eq+zmax_eq)
             rmax_wr=raxis_eq+bdr_threshold*(rmax_eq-raxis_eq)
             rmin_wr=raxis_eq+bdr_threshold*(rmin_eq-raxis_eq)
             zmax_wr=zaxis_eq+bdr_threshold*(zmax_eq-zaxis_eq)
             zmin_wr=zaxis_eq+bdr_threshold*(zmin_eq-zaxis_eq)
             IF(rmin_wr.LT.0.D0) rmin_wr=0.D0
             INITEQ=1
          ELSE
             WRITE(6,*) 'XX EQREAD: IERR=',IERR
             INITEQ=0
          ENDIF
       ENDIF
    ELSE
       INITEQ=0
       raxis_eq=RR
       zaxis_eq=0.D0
       rmax_eq=RR+RA
       rmin_wr=RR-RA
       zmax_wr= RKAP*RA
       zmin_wr=-RKAP*RA
       rmax_wr=raxis_eq+bdr_threshold*(rmax_eq-raxis_eq)
       rmin_wr=raxis_eq+bdr_threshold*(rmin_eq-raxis_eq)
       zmax_wr=zaxis_eq+bdr_threshold*(zmax_eq-zaxis_eq)
       zmin_wr=zaxis_eq+bdr_threshold*(zmin_eq-zaxis_eq)
       IF(rmin_wr.LT.0.D0) rmin_wr=0.D0
       nsumax=128
       IF(ALLOCATED(rsu_wr)) DEALLOCATE(rsu_wr)
       IF(ALLOCATED(zsu_wr)) DEALLOCATE(zsu_wr)
       ALLOCATE(rsu_wr(nsumax+1),zsu_wr(nsumax+1))
       dth=2.D0*PI/nsumax
       DO nsu=1,nsumax+1
          rsu_wr(nsu)=RR+RA*COS((nsu-1)*dth)
          zsu_wr(nsu)=RKAP*RA*SIN((nsu-1)*dth)
       END DO
    ENDIF

    CALL DPPREP_LOCAL(IERR)

    RETURN
  END SUBROUTINE WR_CHEK
END MODULE wrparm
