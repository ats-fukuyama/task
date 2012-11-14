MODULE trgsub

  USE trcomm, ONLY: ikind, rkind

  PRIVATE
  PUBLIC tr_gr_time

CONTAINS

  SUBROUTINE tr_gr_time(idexp)
! ----------------------------------------------------------------------
!            Write time on GSAF pages
!
!  idexp = 0 : simulation time
!        = 1 : time of snap shot of experimental data (mdluf = 1)
!        = 2 : time of snap shot of experimental data (mdluf = 2, 3)
!               (arbitrary time slice every graphic view)
! ----------------------------------------------------------------------       
    USE trcomm, ONLY : t, time_slc, time_snap
    IMPLICIT NONE

    INTEGER(ikind),INTENT(IN) :: idexp

    CALL SETLIN(0,0,7)
    CALL SETCHS(0.3,0.0)
    CALL SETFNT(32)
    CALL MOVE(11.8,18.0)
    CALL TEXT('t =',2)

    SELECT CASE(idexp)
    CASE(0)
       CALL NUMBD(t,'(1F7.3)',7)
    CASE(1)
       CALL NUMBD(time_slc,'(1F7.3)',7)
    CASE(2)
       CALL NUMBD(time_snap,'(1F7.3)',7)
    END SELECT

    CALL TEXT(' sec.',4)

    RETURN
  END SUBROUTINE tr_gr_time

END MODULE trgsub
