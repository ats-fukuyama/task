! wqmain.f90

!               ############# TASK/WQ #############

!       Wave propagation analysis with quasi-optical formulation

!                           A. Fukuyama
!               originally coded by  Y. Maruyama in 2009

!                      V1.00  : 2020 AUGUST 23
!     --------------------------------------------------------

PROGRAM wq

  USE wqcomm
  USE wqinit,ONLY: wq_init
  USE wqparm,ONLY: wq_parm
  USE wqmenu,ONLY: wq_menu
  USE libmtx
  IMPLICIT NONE
  INTEGER:: ierr

  WRITE(6,*) '## TASK/WQ 2020/09/15'
  CALL GSOPEN
  CALL mtx_initialize

  CALL wq_init
  OPEN(7,STATUS='SCRATCH')
  CALL wq_parm(1,'wqparm',IERR)

  CALL wq_menu

  CLOSE(7)
  CALL mtx_finalize
  CALL GSCLOS

  STOP
END PROGRAM wq
