! cvmenu.f90

MODULE cvmenu

  PRIVATE
  PUBLIC cv_menu

CONTAINS

!     ***** TASK/CV MENU *****

  SUBROUTINE cv_menu

    USE cvcomm
    USE cvparm,ONLY: cv_parm
    USE cvview,ONLY: cv_view
    USE cvfile,ONLY: cv_load,cv_global
    USE cvgout,ONLY: cv_gout
    USE cvregion,ONLY: cv_region
    USE cvselect,ONLY: cv_select
    USE cvpopulation,ONLY: cv_population_load
    USE cvrank,ONLY: cv_rank
    USE cvlib
    USE libcharx
    IMPLICIT NONE
    CHARACTER(LEN=1):: kid
    CHARACTER(LEN=80):: line
    INTEGER:: ierr,mode,nid,ndate,ndate1
    CHARACTER(LEN=10):: date_id
    CHARACTER(LEN=10):: kidn
    INTEGER,PARAMETER:: lword=10
    INTEGER:: nword,nword_max
    CHARACTER(LEN=lword),ALLOCATABLE:: kworda(:)

    CALL cv_load(ierr)
    CALL cv_population_load(ierr)

1   CONTINUE
    ierr=0
    WRITE(6,'(A,A)') &
         '## CV MENU: P,V/parm G/graph L/load R/region A/all ', &
         'S/select N/rank Q/quit'
    CALL task_klin(line,kid,mode,cv_parm)
    IF(mode.NE.1) GOTO 1

    IF(kid.EQ.'P') THEN
       CALL cv_parm(0,'CV',ierr)
    ELSEIF(kid.EQ.'V') THEN
       CALL cv_view
    ELSEIF(kid.EQ.'G') THEN
       CALL cv_gout
    ELSEIF(kid.EQ.'R') THEN
       CALL cv_region
    ELSEIF(kid.EQ.'A') THEN
       CALL cv_global(ierr)
    ELSEIF(kid.EQ.'S') THEN
       CALL cv_select
    ELSEIF(kid.EQ.'L') THEN
       CALL cv_load(ierr)
       CALL cv_population_load(ierr)
    ELSEIF(kid.EQ.'N') THEN
       CALL cv_rank
    ELSEIF(kid.EQ.'T') THEN
7001   CONTINUE
       WRITE(6,'(A)') '## cvtest: input ndate:'
       READ(5,*,END=1,ERR=7001) ndate
       CALL convert_ndate_to_date_id(ndate,date_id)
       WRITE(6,'(A,A10)') '  date_id=',date_id
       CALL convert_date_id_to_ndate(date_id,ndate1)
       WRITE(6,'(A,2I10)') '  ndate=',ndate,ndate1
       GO TO 7001
    ELSEIF(kid.EQ.'U') THEN
7002   CONTINUE
       WRITE(6,'(A)') '## ksplitn test: input line?'
       READ(5,'(A)',END=1,ERR=7002) line
       kidn=' ,'
       CALL ksplitn(line,kidn,lword,nword_max,kworda)
       WRITE(6,'(A,I8)') '  nword_max=',nword_max
       DO nword=1,nword_max
          WRITE(6,'(2I8,2X,A)') nword,LEN_TRIM(kworda(nword)),kworda(nword)
       END DO
       GO TO 7002
    ELSEIF(kid.EQ.'Q') THEN
       GO TO 9000
    ELSE
       WRITE(6,*) 'XX cv_menu: unknown kid'
    ENDIF
    GO TO 1

9000  CONTINUE
    RETURN
  END SUBROUTINE cv_menu
END MODULE cvmenu
