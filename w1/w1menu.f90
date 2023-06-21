MODULE w1menu

CONTAINS

  SUBROUTINE w1_menu

    USE w1comm,ONLY: ikind,rkind,w1_allocate,w1_deallocate
    USE w1parm,ONLY: w1_parm,w1_view
    USE w1exec,ONLY: w1_exec
    USE w1gdsp,ONLY: w1_gdsp
    USE w1gout,ONLY: w1_gout
    USE w1file,ONLY: w1_save,w1_load
    USE libkio

    IMPLICIT NONE
    INTEGER(ikind)    :: ierr,mode,ind
    CHARACTER         :: kid
    CHARACTER(LEN=80) :: line
    INTEGER(ikind)    :: init=0

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') &
         '## W1 MENU: P,V/PARM  R/RUN  G/GRAF  D/DISP  S,L/FILE  Q/QUIT'

    CALL TASK_KLIN(line,kid,mode,w1_parm)
    IF(mode /= 1) GOTO 1
!    IF(nxmax.NE.nxmax_save) THEN  ! data structure was modified
!       INIT=0
!    END IF

    SELECT CASE(kid)
    CASE('P')
       CALL w1_parm(0,'W1',ierr)
    CASE('V')
       CALL w1_view
    CASE('R')
       CALL w1_allocate
       CALL w1_exec(ierr)
       IF(ierr.EQ.0) THEN
          INIT=1
       ELSE
          INIT=0
       END IF
    CASE('G')
       IF(INIT.EQ.0) THEN
          WRITE(6,*) 'W1 data is not ready or destroyed'
       ELSE
          CALL w1_gout
       END IF
    CASE('D')
       CALL w1_allocate
       CALL w1_gdsp
       INIT=1
    CASE('S')
       CALL w1_save(ierr)
    CASE('L')
       CALL w1_load(ierr)
       IF(ierr.EQ.0) THEN
          INIT=1
       ELSE
          INIT=0
       END IF
    CASE('Q')
       GOTO 9000
    CASE DEFAULT
       WRITE(6,*) 'W1 W1MENU: UNKNOWN kid: kid = ',kid
    END SELECT
    GOTO 1

9000 CALL w1_deallocate
    RETURN

  END SUBROUTINE w1_menu

END MODULE w1menu
