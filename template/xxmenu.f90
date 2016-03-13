!  $Id$

!  ***** TASK/XX MENU *****

MODULE xxmenu

CONTAINS

  SUBROUTINE xx_menu

    USE xxcomm,ONLY: ikind,rkind,xx_allocate,xx_deallocate,nxmax
    USE xxparm,ONLY: xx_parm,xx_view
    USE xxexec,ONLY: xx_exec
    USE xxgout,ONLY: xx_gout

    IMPLICIT NONE
    INTEGER(ikind)    :: ierr,mode,ind
    CHARACTER         :: kid
    CHARACTER(LEN=80) :: line
    INTEGER(ikind)    :: init=0
    INTEGER(ikind)    :: nxmax_save=0

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') '## XX MENU: P,V/PARM  R/RUN  G/GRAF  Q/QUIT'

    CALL TASK_KLIN(line,kid,mode,xx_parm)
    IF(mode /= 1) GOTO 1
    IF(nxmax.NE.nxmax_save) THEN  ! data structure was modified
       INIT=0
    END IF

    IF(kid.EQ.'P') THEN
       CALL xx_parm(0,'XX',ierr)
    ELSEIF(kid.EQ.'V') THEN
       CALL xx_view
    ELSEIF(kid.EQ.'R') THEN
       CALL xx_allocate
       nxmax_save=nxmax
       CALL xx_exec(ierr)
       INIT=1
    ELSEIF(kid.EQ.'G') THEN
       IF(INIT.EQ.0) THEN
          WRITE(6,*) 'XX data is not ready or destroyed'
       ELSE
          CALL xx_gout
       END IF
    ELSEIF(kid.EQ.'Q') THEN
       GOTO 9000
    ELSE
       WRITE(6,*) 'XX XXMENU: UNKNOWN kid: kid = ',kid
    ENDIF
    GOTO 1

9000 CALL xx_deallocate
    RETURN
  END SUBROUTINE xx_menu
END MODULE xxmenu
