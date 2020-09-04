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
    USE cvfile,ONLY: cv_read,cv_write
    USE cvgout,ONLY: cv_gout
    USE cvregion,ONLY: cv_region
    USE cvselect,ONLY: cv_select
    USE cvpopulation,ONLY: cv_population_read
    USE cvlib
    IMPLICIT NONE
    CHARACTER(LEN=1):: kid
    CHARACTER(LEN=80):: line
    INTEGER:: nstat,ierr,mode,nid,ndate,ndate1
    CHARACTER(LEN=10):: date_id

    nstat=0

1   CONTINUE
    ierr=0
    WRITE(6,'(A,A)') &
         '## CV MENU: P,V/parm  L/load  C/convert  G/graph  R/region  ', &
         'S/select  Q/quit'
    CALL task_klin(line,kid,mode,cv_parm)
    IF(mode.NE.1) GOTO 1

    IF(kid.EQ.'P') THEN
       CALL cv_parm(0,'CV',ierr)
    ELSEIF(kid.EQ.'V') THEN
       CALL cv_view
    ELSEIF(kid.EQ.'L') THEN
       CALL cv_read(ierr)
       CALL cv_population_read(ierr)
       nstat=1
    ELSEIF(kid.EQ.'C') THEN
       IF(nstat.EQ.0) THEN
          WRITE(6,'(A)') 'XX cvmenu: No date ha been loaded'
       ELSE
          CALL cv_write(ierr)
       END IF
    ELSEIF(kid.EQ.'G') THEN
       IF(nstat.EQ.0) THEN
          WRITE(6,'(A)') 'XX cvmenu: No date ha been loaded'
       ELSE
          CALL cv_gout
       END IF
    ELSEIF(kid.EQ.'R') THEN
       IF(nstat.EQ.0) THEN
          WRITE(6,'(A)') 'XX cvmenu: No date ha been loaded'
       ELSE
          CALL cv_region
       END IF
    ELSEIF(kid.EQ.'S') THEN
       IF(nstat.EQ.0) THEN
          WRITE(6,'(A)') 'XX cvmenu: No date ha been loaded'
       ELSE
          CALL cv_select
       END IF
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
       WRITE(6,*) 'XX cv_menu: unknown kid'
    ENDIF
    GO TO 1

9000  CONTINUE
    RETURN
  END SUBROUTINE cv_menu
END MODULE cvmenu
