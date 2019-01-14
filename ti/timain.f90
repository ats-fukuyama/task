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
  USE plparm,ONLY: pl_parm
  USE equnit_mod
  USE tiinit,ONLY: ti_init
  USE tiparm,ONLY: ti_parm
  USE timenu,ONLY: ti_menu
  USE ADPOST
  USE libmtx
  IMPLICIT NONE
  INTEGER(ikind):: IERR

!     ------ INITIALIZATION ------

  CALL mtx_initialize

  IF(nrank.EQ.0) THEN
     WRITE(6,*) '***** TASK/TI  19/01/11 *****'
     CALL GSOPEN
     OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
     CALL read_adpost(IERR)
     IF(IERR.NE.0) WRITE(6,*) 'XX read_adpost: IERR=',IERR
  END IF

  CALL pl_init
  call eq_init
  CALL ti_init

  CALL ti_parm(1,'tiparm',IERR)

  CALL ti_menu

  IF(nrank.EQ.0) THEN
     CALL GSCLOS
  END IF
  CALL mtx_finalize
  STOP
END PROGRAM ti
