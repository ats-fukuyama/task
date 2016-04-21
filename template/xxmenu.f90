MODULE xxmenu

PRIVATE
PUBLIC xx_menu

CONTAINS

  SUBROUTINE xx_menu

    USE xxcomm,ONLY: ikind,rkind,xx_allocate,xx_deallocate
    USE xxparm,ONLY: xx_parm,xx_view
    USE xxexec,ONLY: xx_exec
    USE xxgout,ONLY: xx_gout

    IMPLICIT NONE
    INTEGER(ikind)    :: ierr,mode,ind
    CHARACTER         :: kid
    CHARACTER(LEN=80) :: line

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') '## XX MENU: P,V/PARM  R/RUN  G/GRAF  Q/QUIT'

    CALL TASK_KLIN(line,kid,mode,xx_parm)
    IF(mode /= 1) GOTO 1

    IF(kid.EQ.'P') THEN
       CALL xx_parm(0,'XX',ierr)
    ELSEIF(kid.EQ.'V') THEN
       CALL xx_view
    ELSEIF(kid.EQ.'R') THEN
       CALL xx_allocate
       CALL xx_exec(ierr)
    ELSEIF(kid.EQ.'G') THEN
       CALL xx_gout
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
