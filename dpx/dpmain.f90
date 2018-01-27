PROGRAM DP

!               ############# TASK/DP #############

!                       DISPERSION RELATION

!                           A. Fukuyama

!                Department of Nuclear Engineering
!                       Kyoto Univerisity
!                     Kyoto 606-8501, Japan

!                      V1.00  : 1988 FEB 10
!                      V2.00  : 1996 FEB 04
!                      V3.00  : 1997 AUG 05
!                      V3.10  : 2000 NOV 25

!-----------------------------------------------------------------------

  USE plinit, ONLY: pl_init,pl_parm
  USE dpinit, ONLY: dp_init
  USE dpparm, ONLY: dp_parm
  USE dpmenu, ONLY: dp_menu
  IMPLICIt NONE
  INTEGER:: IERR

  WRITE(6,*) '## TASK/DP 2017/08/18'
  CALL GSOPEN
  OPEN(7,STATUS='SCRATCH')
  CALL PL_INIT
  CALL DP_INIT
  CALL PL_PARM(1,'plparm',IERR)
  CALL DP_PARM(1,'dpparm',IERR)
  
  CALL DP_MENU

  CLOSE(7)
  CALL GSCLOS
  STOP
END program DP
