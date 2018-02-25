!     *************************** TASK.TI ***************************
!     *                                                             *
!     *  IMPURITY TRANSPORT IN A TOKAMAK                            *
!     *                                                             *
!     *  PROGRAMMED BY A. FUKUYAMA                                  *
!     *                                                             *
!     ***************************************************************

PROGRAM ti

  USE ticomm
  USE bpsd
  USE plinit,ONLY: pl_init,pl_parm
  USE equnit_mod
  USE tiinit,ONLY: ti_init
  USE tiparm,ONLY: ti_parm
  USE timenu,ONLY: ti_menu
  IMPLICIT NONE
  INTEGER(ikind):: IERR

!     ------ INITIALIZATION ------

  WRITE(6,*) '***** TASK/TI  18/02/24 *****'

  CALL GSOPEN

  CALL pl_init
  call eq_init
  CALL ti_init

  CALL pl_parm(1,'plparm',IERR)
  CALL eq_parm(1,'eqparm',IERR)
  CALL ti_parm(1,'tiparm',IERR)

  CALL ti_menu

  CALL GSCLOS
  STOP
END PROGRAM ti
