! fpmain.f90

!     *******************************************************

!                    3D FOKKER-PLANCK CODE

!                  TASK/FP  V2.0  1992/12/12
!                           V2.1  1993/09/16
!                           V2.2  1997/03/18
!                           V3.0  1997/08/05
!                           V4.0  2021/08/08

!                         PROGRAMMED BY
!                      FUKUYAMA Atsushi
!                      NUGA Hideo (Relativistic effects, runaway electrons)
!                      OTA Keigo (Finite orbit width effects)
!                         Kyoto University

!     ********************************************************

PROGRAM fp

  USE fpcomm
  USE plinit
  USE plparm
  USE equnit_mod
  USE fpinit
  USE fpparm
  USE fpmenu
  USE fpwrin
  USE libmtx

  IMPLICIT NONE
  INTEGER:: ierr

  CALL mtx_initialize
  IF(nrank.EQ.0) THEN
     WRITE(6,*) '***** TASK/FP V4.0 2021/08/08 *****'
     CALL GSOPEN
     OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
  ENDIF

  CALL pl_init
  CALL eq_init
  CALL fp_init
  IF(nrank.EQ.0) THEN
     CALL fp_parm(1,'fpparm',IERR)
  ENDIF
  CALL fp_menu

  IF(nrank.EQ.0) CALL GSCLOS
  CALL mtx_finalize
  CALL fp_wr_deallocate
  STOP
END PROGRAM fp
