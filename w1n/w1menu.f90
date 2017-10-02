MODULE w1menu

CONTAINS

  SUBROUTINE w1_menu

    USE w1comm,ONLY: ikind,rkind,w1_allocate,w1_deallocate
    USE w1parm,ONLY: w1_parm,w1_view
    USE w1exec,ONLY: w1_exec
    USE w1gout,ONLY: w1_gout

    IMPLICIT NONE
    INTEGER(ikind)    :: ierr,mode,ind
    CHARACTER         :: kid
    CHARACTER(LEN=80) :: line
    INTEGER(ikind)    :: init=0

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') '## W1 MENU: P,V/PARM  R/RUN  G/GRAF  D/DISP  Q/QUIT'

    CALL TASK_KLIN(line,kid,mode,w1_parm)
    IF(mode /= 1) GOTO 1
!    IF(nxmax.NE.nxmax_save) THEN  ! data structure was modified
!       INIT=0
!    END IF

    IF(kid.EQ.'P') THEN
       CALL w1_parm(0,'W1',ierr)
    ELSEIF(kid.EQ.'V') THEN
       CALL w1_view
    ELSEIF(kid.EQ.'R') THEN
       CALL w1_allocate
       CALL w1_exec(ierr)
       INIT=1
    ELSEIF(kid.EQ.'G') THEN
       IF(INIT.EQ.0) THEN
          WRITE(6,*) 'W1 data is not ready or destroyed'
       ELSE
          CALL w1_gout
       END IF
    ELSEIF(kid.EQ.'D') THEN
       CALL w1_allocate
       CALL w1_disp(ierr)
       INIT=1
    ELSEIF(kid.EQ.'Q') THEN
       GOTO 9000
    ELSE
       WRITE(6,*) 'W1 W1MENU: UNKNOWN kid: kid = ',kid
    ENDIF
    GOTO 1

9000 CALL w1_deallocate
    RETURN

  END SUBROUTINE w1_menu

END MODULE w1menu
