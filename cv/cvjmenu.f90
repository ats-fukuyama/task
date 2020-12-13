! cvjmenu.f90

MODULE cvjmenu

  PRIVATE
  PUBLIC cvj_menu

CONTAINS

!     ***** TASK/CVJ MENU *****

  SUBROUTINE cvj_menu

    USE cvcomm
    USE cvparm,ONLY: cv_parm
    USE cvview,ONLY: cv_view
    USE cvjfile,ONLY: cvj_load
    USE cvgout,ONLY: cv_gout
    USE cvjpopulation,ONLY: cvj_population_load
    USE cvrank,ONLY: cv_rank
    USE cvlib
    USE libcharx
    IMPLICIT NONE
    CHARACTER(LEN=1):: kid
    CHARACTER(LEN=80):: line
    INTEGER:: ierr,mode,ndate,ndate1
    CHARACTER(LEN=10):: date_id
    CHARACTER(LEN=10):: kidn
    INTEGER,PARAMETER:: lword=10
    INTEGER:: nword,nword_max
    CHARACTER(LEN=lword),ALLOCATABLE:: kworda(:)
    EXTERNAL TASK_KLIN

    CALL cvj_load(ierr)
    CALL cvj_population_load(ierr)
    population_min_rank=0.D0

1   CONTINUE
    ierr=0
    WRITE(6,'(A,A)') &
         '## CVJ MENU: P,V/parm G/graph L/load ', &
         'N/rank Q/quit'
    CALL task_klin(line,kid,mode,cv_parm)
    IF(mode.NE.1) GOTO 1

    IF(kid.EQ.'P') THEN
       CALL cv_parm(0,'CV',ierr)
    ELSEIF(kid.EQ.'V') THEN
       CALL cv_view
    ELSEIF(kid.EQ.'G') THEN
       CALL cv_gout
    ELSEIF(kid.EQ.'L') THEN
       CALL cvj_load(ierr)
       CALL cvj_population_load(ierr)
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
    ELSEIF(kid.EQ.'Q') THEN
       GO TO 9000
    ELSE
       WRITE(6,*) 'XX cvj_menu: unknown kid'
    ENDIF
    GO TO 1

9000  CONTINUE
    RETURN
  END SUBROUTINE cvj_menu
END MODULE cvjmenu
