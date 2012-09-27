MODULE trgsub

  USE trcomm, ONLY: ikind, rkind

  PRIVATE
  PUBLIC tr_gr_time

CONTAINS

  SUBROUTINE tr_gr_time
! ***********************************************************************
!            Write time on figure
! ***********************************************************************

    USE trcomm, ONLY : t
    IMPLICIT NONE

    CALL SETLIN(0,0,7)
    CALL SETCHS(0.3,0.0)
    CALL SETFNT(32)
    CALL MOVE(11.8,18.0)
    CALL TEXT('t =',2)
    CALL NUMBD(t,'(1F7.3)',7)
    CALL TEXT(' sec.',4)

    RETURN
  END SUBROUTINE tr_gr_time

END MODULE trgsub
