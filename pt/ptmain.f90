!     *************************** TASK.pt ***************************
!     *                                                             *
!     *  PLOT diagram                                               *
!     *                                                             *
!     *  PROGRAMMED BY A. FUKUYAMA                                  *
!     *                                                             *
!     ***************************************************************

PROGRAM pt

  USE ptcomm
  USE ptinit,ONLY: pt_init
  USE ptparm,ONLY: pt_parm
  USE ptmenu,ONLY: pt_menu
  USE libmtx
  IMPLICIT NONE
  INTEGER:: IERR

!     ------ INITIALIZATION ------

  CALL mtx_initialize

!  IF(nrank.EQ.0) THEN
     WRITE(6,*) '***** TASK/PT  20/11/02 *****'
     CALL GSOPEN
     OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
!  END IF

  CALL pt_init

  CALL pt_parm(1,'ptparm',IERR)

  CALL pt_menu

!  IF(nrank.EQ.0) THEN
     CALL GSCLOS
!  END IF
  CALL mtx_finalize
  STOP
END PROGRAM pt
