! obparm.f90

MODULE obparm

  PRIVATE
  PUBLIC ob_parm,ob_view

CONTAINS

!     ****** INPUT PARAMETERS ******

  SUBROUTINE OB_PARM(MODE,KIN,IERR)

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

    IMPLICIT NONE
    INTEGER,INTENT(IN):: MODE
    CHARACTER(LEN=*),INTENT(IN):: KIN
    INTEGER,INTENT(OUT):: IERR

    IERR=0

1   CALL TASK_PARM(MODE,'OB',KIN,OBNLIN,OBPLST,IERR)
    IF(IERR.NE.0) RETURN

    CALl EQCHEK(IERR)
    CALl OB_CHEK(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1

    RETURN
  END SUBROUTINE OB_PARM

!     ****** INPUT NAMELIST ******

  SUBROUTINE OBNLIN(NID,IST,IERR)

    USE obcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NID
    INTEGER,INTENT(OUT):: IST,IERR
    INTEGER:: NS

    NAMELIST /OB/ RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
                  NSMAX,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PZCL, &
                  r_corner,z_corner, &
                  br_corner,bz_corner,bt_corner, &
                  pn_corner,ptpr_corner,ptpp_corner, &
                  PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
                  RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG, &
                  PPN0,PTN0,RF_PL, &
                  MODELG,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF, &
                  RHOGMN,RHOGMX, &
                  KNAMEQ,KNAMFP,KNAMFO,KNAMEQ2, &
                  MODEFW,MODEFR,IDEBUG, &
                  RF,RPI,ZPI,PHII,RNZI,RNPHII,RKR0,MODEW,UUI, &
                  RCURVA,RCURVB,RBRADA,RBRADB, &
                  RFIN,RPIN,ZPIN,PHIIN,RNZIN,RNPHIIN,RKRIN,MODEWIN,UUIN, &
                  ANGZIN,ANGPHIN,RCURVAIN,RCURVBIN,RBRADAIN,RBRADBIN, &
                  NOBTMAX,NSTPMAX,NRSMAX,NRRMAX,LMAXNW, &
                  MDLOBI,MDLOBG,MDLOBP,MDLOBQ,MDLOBW, &
                  SMAX,DELS,UUMIN,EPSOBT,DELOBT,DELDER,DELKR,EPSNW

    READ(NID,OB,IOSTAT=IST,ERR=9800,END=9900)
    
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
  END SUBROUTINE OBNLIN

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE OBPLST

    WRITE(6,601)
    RETURN

  601 FORMAT(' ','# &OB : RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ,'/ &
             9X,'NSMAX,PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PZCL,'/ &
             9X,'PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2,'/ &
             9X,'r_corner,z_corner,br_corner,bz_corner,bt_corner,'/ &
             9X,'pn_corner,ptpr_corner,ptpp_corner,'/ &
             9X,'RHOMIN,QMIN,RHOITB,PNITB,PTITB,PUITB,RHOEDG,'/ &
             9X,'PPN0,PTN0,RFCL,'/ &
             9X,'MODELG,MODELN,MODELQ,MODEL_PROF,MODEL_NPROF,'/ &
             9X,'RHOGMN,RHOGMX,'/ &
             9X,'KNAMEQ,KNAMOB,KNAMFP,KNAMFO,KNAMEQ2'/ &
             9X,'MODEFW,MODEFR,IDEBUG'/ &
             9X,'RF,RPI,ZPI,PHII,RNZI,RNPHII,RKR0,MODEW,UUI,'/ &
             9X,'RCURVA,RCURVB,RBRADA,RBRADB,'/ &
             9X,'RFIN,RPIN,ZPIN,PHIIN,RNZIN,RNPHIIN,RKR0IN,MODEWIN,UUIN,'/ &
             9X,'ANGZIN,ANGPHIN,RCURVAIN,RCURVBIN,RBRADAIN,RBRADBIN,'/ &
             9X,'NOBTMAX,NSTPMAX,NRSMAX,NRRMAX,LMAXNW,'/ &
             9X,'MDLOBI,MDLOBG,MDLOBP,MDLOBQ,MDLOBW,'/ &
             9X,'SMAX,DELS,UUMIN,EPSOBT,DELOBT,DELDER,DELKR,EPSNW')
  END SUBROUTINE OBPLST

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE OB_CHEK(IERR)

    USE obcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    CHARACTER(LEN=80):: LINE
    INTEGER,SAVE:: INITEQ=0

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

    RETURN
  END SUBROUTINE OB_CHEK

!     ****** SHOW PARAMETERS ******

  SUBROUTINE OB_VIEW

    USE obcomm_parm
    IMPLICIT NONE

    WRITE(6,603) 'NOBTMAX     ',NOBTMAX, &
                 'NSTPMAX     ',NSTPMAX
    WRITE(6,603) 'NRSMAX      ',NRSMAX, &
                 'NRRMAX      ',NRRMAX
    WRITE(6,601) 'SMAX  ',SMAX  ,'DELS  ',DELS  , &
                 'UUMIN ',UUMIN
    WRITE(6,601) 'EPSOBT',EPSOBT,'DELOBT',DELOBT, &
                 'DELDER',DELDER
    WRITE(6,601) 'DELKR ',DELKR, 'EPSNW ',EPSNW
    WRITE(6,602) 'LMAXNW',LMAXNW
    WRITE(6,602) 'MDLOBI',MDLOBI,'MDLOBG',MDLOBG, &
                 'MDLOBP',MDLOBP
    WRITE(6,602) 'MDLOBQ',MDLOBQ,'MDLOBW',MDLOBW
    RETURN

601 FORMAT(1H ,A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
            2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
602 FORMAT(1H ,A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
            2X,A6,'=',I7,4X  :2X,A6,'=',I7)
603 FORMAT(1H ,A12,'=',I7,4X  :2X,A12,'=',I7,4X  : &
            2X,A12,'=',I7)
  END SUBROUTINE OB_VIEW
END MODULE obparm
