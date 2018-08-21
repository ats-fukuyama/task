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
  USE dpcomm, ONLY: BB,nrank
  use libmtx
  IMPLICIt NONE
  INTEGER:: IERR

  CALL mtx_initialize
  IF(nrank.EQ.0) THEN
     WRITE(6,*) '## TASK/DP 2018/08/21'
     CALL GSOPEN
     OPEN(7,STATUS='SCRATCH')
  END IF
  CALL PL_INIT
  CALL DP_INIT
  CALL PL_PARM(1,'plparm',IERR)
  CALL DP_PARM(1,'dpparm',IERR)
  
  CALL DP_MENU

  IF(nrank.EQ.0) THEN
     CLOSE(7)
     CALL GSCLOS
  END IF
  STOP
END program DP
