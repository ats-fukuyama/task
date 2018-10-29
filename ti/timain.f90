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
  USE plinit,ONLY: pl_init
  USE plparm,ONLY: pl_parm
  USE equnit_mod
  USE tiinit,ONLY: ti_init
  USE tiparm,ONLY: ti_parm
  USE timenu,ONLY: ti_menu
  USE ADPOST
  IMPLICIT NONE
  INTEGER(ikind):: IERR

!     ------ INITIALIZATION ------

  WRITE(6,*) '***** TASK/TI  18/02/24 *****'

  CALL GSOPEN
  OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
  CALL read_adpost(IERR)
  IF(IERR.NE.0) WRITE(6,*) 'XX read_adpost: IERR=',IERR

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
