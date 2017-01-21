!  ***** TASK/WIM MENU *****

MODULE wimmenu

CONTAINS

  SUBROUTINE wim_menu

    USE wimcomm,ONLY: ikind,rkind,wim_allocate,wim_deallocate, &
         nzmax,nwmax,modelp
    USE wimparm,ONLY: wim_parm,wim_view
    USE wimexec,ONLY: wim_exec,subfw
    USE wimgout,ONLY: wim_gout

    IMPLICIT NONE
    INTEGER(ikind)    :: ierr,mode,ind
    CHARACTER         :: kid
    CHARACTER(LEN=80) :: line
    INTEGER(ikind)    :: init=0
    INTEGER(ikind)    :: nzmax_save=0
    INTEGER(ikind)    :: nwmax_save=0

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') '## WIM MENU: P,V/PARM  R/RUN  G/GRAF  K/Kernel Q/QUIT'

    CALL TASK_KLIN(line,kid,mode,wim_parm)
    IF(mode /= 1) GOTO 1
    IF(nzmax.NE.nzmax_save) THEN  ! data structure was modified
       INIT=0
    END IF

    IF(kid.EQ.'P') THEN
       CALL wim_parm(0,'WIM',ierr)
    ELSEIF(kid.EQ.'V') THEN
       CALL wim_view
    ELSEIF(kid.EQ.'R') THEN
       NWMAX_SAVE=NWMAX
       IF(MODELP.EQ.0) NWMAX=1
       CALL wim_allocate
       nzmax_save=nzmax
       CALL wim_exec(ierr)
       NWMAX=NWMAX_SAVE
       NZMAX_SAVE=NZMAX
       INIT=1
    ELSEIF(kid.EQ.'G') THEN
       IF(INIT.EQ.0) THEN
          WRITE(6,*) 'WIM data is not ready or destroyed'
       ELSE
          CALL wim_gout
       END IF
    ELSEIF(kid.EQ.'K') THEN
       CALL wim_allocate
       CALL SUBFW(1)
    ELSEIF(kid.EQ.'Q') THEN
       GOTO 9000
    ELSE
       WRITE(6,*) 'WIM WIMMENU: UNKNOWN kid: kid = ',kid
    ENDIF
    GOTO 1

9000 CALL wim_deallocate
    RETURN
  END SUBROUTINE wim_menu
END MODULE wimmenu
