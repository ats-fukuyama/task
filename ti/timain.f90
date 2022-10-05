!     *************************** TASK.TI ***************************
!     *                                                             *
!     *  IMPURITY TRANSPORT IN A TOKAMAK                            *
!     *                                                             *
!     *  PROGRAMMED BY A. FUKUYAMA                                  *
!     *                                                             *
!     ***************************************************************

PROGRAM ti

  USE ticomm
  USE plinit,ONLY: pl_init
  USE equnit,ONLY: eq_init
  USE tiinit,ONLY: ti_init
  USE tiparm,ONLY: ti_parm,ti_broadcast
  USE timenu,ONLY: ti_menu
  USE libmtx
  IMPLICIT NONE
  INTEGER(ikind):: IERR

!     ------ INITIALIZATION ------

  CALL mtx_initialize

  IF(nrank.EQ.0) THEN
     WRITE(6,*) '***** TASK/TI  19/01/11 *****'
     CALL GSOPEN
     OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
  END IF

  CALL pl_init
  call eq_init
  CALL ti_init

  IF(nrank.EQ.0) THEN
     CALL ti_parm(1,'tiparm',IERR)
  END IF
  CALL ti_broadcast

  CALL ti_menu

  IF(nrank.EQ.0) THEN
     CALL GSCLOS
  END IF
  CALL mtx_finalize
  STOP
END PROGRAM ti
