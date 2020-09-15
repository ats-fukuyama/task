! wqmenu.f90

MODULE wqmenu

  PRIVATE
  PUBLIC wq_menu

CONTAINS

!     ***** TASK/WQ MENU *****

  SUBROUTINE wq_menu

    USE wqcomm
    USE wqparm,ONLY: wq_parm
    USE wqview,ONLY: wq_view
    USE wqprep,ONLY: wq_prep
    USE wqexec,ONLY: wq_exec
    USE wqgout,ONLY: wq_gout
    IMPLICIT NONE
    CHARACTER(LEN=1):: kid
    CHARACTER(LEN=80):: line
    INTEGER:: ierr,mode
    INTEGER:: istatus=0

1   CONTINUE
    ierr=0
    WRITE(6,'(A,A)') &
         '## WQ MENU: P,V/parm  R/run  C/cont  G/graph  Q/quit'
    CALL task_klin(line,kid,mode,wq_parm)
    IF(mode.NE.1) GOTO 1

    IF(kid.EQ.'P') THEN
       CALL wq_parm(0,'WQ',ierr)
    ELSEIF(kid.EQ.'V') THEN
       CALL wq_view
    ELSEIF(kid.EQ.'R') THEN
       CALL wq_allocate
       CALL wq_prep
       istatus=1
       CALL wq_exec
    ELSEIF(kid.EQ.'C') THEN
       IF(istatus.EQ.1) THEN
          CALL wq_exec
       ELSE
          WRITE(6,'(A)') 'XX wq_menu: No result to be continued'
       END IF
    ELSEIF(kid.EQ.'G') THEN
       IF(istatus.EQ.1) THEN
          CALL wq_gout
       ELSE
          WRITE(6,'(A)') 'XX wq_menu: No result to be drawn'
       END IF
    ELSEIF(kid.EQ.'Q') THEN
       GO TO 9000
    ELSE
       WRITE(6,*) 'XX wq_menu: unknown kid'
    ENDIF
    GO TO 1

9000  CONTINUE
    RETURN
  END SUBROUTINE wq_menu
END MODULE wqmenu
