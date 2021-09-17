! fowmenu.f90

MODULE fowmenu

  PRIVATE
  PUBLIC fow_menu

CONTAINS

  SUBROUTINE fow_menu
    use fowcomm
    use fowprep
    use foworbit
    use fowdistribution
    use fowloop

    USE fpparm
    USE fpprep,ONLY: fp_prep
    use fpwrite
    USE fpsub,ONLY: FPMXWL
    use equnit_mod
    USE libkio
    use libmpi

    implicit none

    CHARACTER(LEN=1)::  KID
    CHARACTER(LEN=80):: LINE
    INTEGER:: MODE,IERR

1   CONTINUE
    IF(nrank.EQ.0) THEN
       ierr=0
       WRITE(6,601)
601    FORMAT('## FOW MENU: R:RUN Q:QUIT')
       CALL TASK_KLIN(LINE,KID,MODE,fp_parm)
    ENDIF
    CALL mtx_barrier
    CALL mtx_broadcast_character(KID,1)
    CALL mtx_broadcast1_integer(MODE)
    IF(MODE.NE.1) GOTO 1


    SELECT CASE(kid)
    CASE('R')
       call fp_broadcast
       call fp_prep(IERR)

       call fow_allocate
       call fow_prep(ierr)

       call fow_loop

       write(*,*)"end"
       call fow_deallocate
    CASE('Q')
       GO TO 9000
    END SELECT
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE fow_menu

END MODULE fowmenu
