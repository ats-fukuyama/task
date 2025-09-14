!     $Id: fpmain.f90,v 1.9 2013/01/22 16:21:46 fukuyama Exp $

!     *******************************************************

!                    3D FOKKER-PLANCK CODE

!                  TASK/FP  V2.0  1992/12/12
!                           V2.1  1993/09/16
!                           V2.2  1997/03/18
!                           V3.0  1997/08/05
!                           V4.0  2011/07/14 with Radial transport 
!                           V5.0  2023/01/20 with Finite orbit width effects
!                           V5.1  2024/09/04 combining update by H. Nuga

!                           PROGRAMMED BY
!                  A. Fukuyama, H. Nuga, K. Ota
!                      Okayama/Kyoto University

!     ********************************************************

PROGRAM fp

  USE fpcomm
  USE plinit
  USE plparm
  USE equnit
  USE obinit
  USE obparm
  USE fpinit
  USE fpparm
  USE fpmenu
  USE libmtx

  IMPLICIT NONE
  INTEGER:: IERR

  CALL mtx_initialize
  IF(nrank.EQ.0) THEN
     WRITE(6,*) '***** TASK/FP 2024/09/04 *****'
     CALL GSOPEN
     OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
  ENDIF

  CALL pl_init
  CALL eq_init
  CALL ob_init
  CALL fp_init
  IF(nrank.EQ.0) THEN
     CALL pl_parm(1,'plparm',IERR)
     CALL eqparm(1,'eqparm',IERR)
     CALL ob_parm(1,'obparm',IERR)
     CALL fp_parm(1,'fpparm',IERR)
  ENDIF
  
  CALL fp_menu

  IF(nrank.EQ.0) CALL GSCLOS
  CALL mtx_finalize
  STOP
END PROGRAM fp
