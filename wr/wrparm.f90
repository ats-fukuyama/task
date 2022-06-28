! wrparm.f90

MODULE wrparm

  PRIVATE
  PUBLIC wr_parm

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
                  MODELG,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF, &
                  RHOGMN,RHOGMX, &
                  KNAMEQ,KNAMWR,KNAMFP,KNAMFO,KNAMEQ2, &
                  MODEFW,MODEFR,IDEBUG, &
                  MODELP,MODELV,NCMIN,NCMAX,NS_NSA_DP,PMAX,EMAX, &
                  RF,RPI,ZPI,PHII,RNZI,RNPHII,RKR0,MODEW,UUI, &
                  RCURVA,RCURVB,RBRADA,RBRADB, &
                  RFIN,RPIN,ZPIN,PHIIN,RNZIN,RNPHIIN,RKRIN,MODEWIN,UUIN, &
                  ANGZIN,ANGPHIN,RCURVAIN,RCURVBIN,RBRADAIN,RBRADBIN, &
                  NRAYMAX,NSTPMAX,NRSMAX,NRLMAX,LMAXNW, &
                  NPMAX_DP,NTHMAX_DP,NRMAX_DP, &
                  MDLWRI,MDLWRG,MDLWRP,MDLWRQ,MDLWRW,nres_max,nres_type, &
                  SMAX,DELS,UUMIN,EPSRAY,DELRAY,DELDER,DELKR,EPSNW, &
                  mode_beam

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
             9X,'MODELG,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF,'/ &
             9X,'RHOGMN,RHOGMX,'/ &
             9X,'KNAMEQ,KNAMWR,KNAMFP,KNAMFO,KNAMEQ2'/ &
             9X,'MODEFW,MODEFR,IDEBUG'/ &
             9X,'MODELP,MODELV,NCMIN,NCMAX,NS_NSA_DP,PMAX,EMAX,'/ &
             9X,'RF,RPI,ZPI,PHII,RNZI,RNPHII,RKR0,MODEW,UUI,'/ &
             9X,'RCURVA,RCURVB,RBRADA,RBRADB,'/ &
             9X,'RFIN,RPIN,ZPIN,PHIIN,RNZIN,RNPHIIN,RKR0IN,MODEWIN,UUIN,'/ &
             9X,'ANGZIN,ANGPHIN,RCURVAIN,RCURVBIN,RBRADAIN,RBRADBIN,'/ &
             9X,'NRAYMAX,NSTPMAX,NRSMAX,NRLMAX,LMAXNW,'/ &
             9X,'NPMAX_DP,NTHMAX_DP,NRMAX_DP,'/ &
             9X,'MDLWRI,MDLWRG,MDLWRP,MDLWRQ,MDLWRW,nres_max,nres_type,'/ &
             9X,'SMAX,DELS,UUMIN,EPSRAY,DELRAY,DELDER,DELKR,EPSNW'/ &
             9X,'mode_beam')
  END SUBROUTINE WRPLST

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE WR_CHEK(IERR)

    USE wrcomm_parm
    USE dpparm,ONLY: dpprep_local
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    CHARACTER(LEN=80):: LINE
    INTEGER,SAVE:: INITEQ=0
    EXTERNAL EQLOAD,EQPARM,EQCALC,EQCALQ,EQGETB,EQREAD

    IERR=0

    IF(MODELG.EQ.3.OR.MODELG.EQ.5) THEN
       IF(INITEQ.EQ.0) THEN
          CALL EQLOAD(MODELG,KNAMEQ,IERR)
          IF(IERR.EQ.0) THEN
             WRITE(LINE,'(A,I5)') 'NRMAX =',51
             CALL EQPARM(2,LINE,IERR)
             WRITE(LINE,'(A,I5)') 'NTHMAX=',64
             CALL EQPARM(2,LINE,IERR)
             WRITE(LINE,'(A,I5)') 'NSUMAX=',64
             CALL EQPARM(2,LINE,IERR)
             CALL EQCALQ(IERR)
             CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
             INITEQ=1
          ELSE
             WRITE(6,*) 'XX EQLOAD: IERR=',IERR
             INITEQ=0
          ENDIF
       ENDIF
    ELSE IF(MODELG.EQ.8) THEN
       IF(INITEQ.EQ.0) THEN
          CALL EQREAD(IERR)
          IF(IERR.EQ.0) THEN
             CALL EQGETB(BB,RR,RIP,RA,RKAP,RDLT,RB)
             INITEQ=1
          ELSE
             WRITE(6,*) 'XX EQREAD: IERR=',IERR
             INITEQ=0
          ENDIF
       ENDIF
    ELSE
       INITEQ=0
    ENDIF

    CALL DPPREP_LOCAL(IERR)

    RETURN
  END SUBROUTINE WR_CHEK
END MODULE wrparm
