!  $Id$

!  ***** TASK/WI MENU *****

MODULE wimenu

CONTAINS

  SUBROUTINE wi_menu

    USE wicomm,ONLY: ikind,rkind,wi_allocate,wi_deallocate,nxmax
    USE wiparm,ONLY: wi_parm,wi_view
    USE wiexec,ONLY: wi_exec
    USE wiscan,ONLY: wi_scan
    USE wigout,ONLY: wi_gout

    IMPLICIT NONE
    INTEGER(ikind)    :: ierr,mode,ind
    CHARACTER         :: kid
    CHARACTER(LEN=80) :: line
    INTEGER(ikind)    :: init=0
    INTEGER(ikind)    :: nxmax_save=0
    REAL(rkind)       :: ratea

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') '## WI MENU: P,V/PARM  R/RUN  S/SCAN  G/GRAF  Q/QUIT'

    CALL TASK_KLIN(line,kid,mode,wi_parm)
    IF(mode /= 1) GOTO 1
    IF(nxmax.NE.nxmax_save) THEN  ! data structure was modified
       INIT=0
    END IF

    IF(kid.EQ.'P') THEN
       CALL wi_parm(0,'WI',ierr)
    ELSEIF(kid.EQ.'V') THEN
       CALL wi_view
    ELSEIF(kid.EQ.'R') THEN
       CALL wi_allocate
       nxmax_save=nxmax
       CALL wi_exec(ierr)
       INIT=1
    ELSEIF(kid.EQ.'S') THEN
       CALL wi_allocate
       nxmax_save=nxmax
       CALL wi_scan(ierr)
       INIT=1
    ELSEIF(kid.EQ.'G') THEN
       IF(INIT.EQ.0) THEN
          WRITE(6,*) 'XX data is not ready or destroyed'
       ELSE
          CALL wi_gout
       END IF
    ELSEIF(kid.EQ.'Q') THEN
       GOTO 9000
    ELSE
       WRITE(6,*) 'XX WIMENU: UNKNOWN kid: kid = ',kid
    ENDIF
    GOTO 1

9000 CALL wi_deallocate
    RETURN
  END SUBROUTINE wi_menu
END MODULE wimenu
