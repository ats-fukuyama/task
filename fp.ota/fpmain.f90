!     $Id: fpmain.f90,v 1.9 2013/01/22 16:21:46 fukuyama Exp $

!     *******************************************************

!                    3D FOKKER-PLANCK CODE

!                  TASK/FP  V2.0  1992/12/12
!                           V2.1  1993/09/16
!                           V2.2  1997/03/18
!                           V3.0  1997/08/05

!                        PROGRAMMED BY
!                         A. FUKUYAMA
!                      OKAYAMA UNIVERSITY

!     ********************************************************
      program fp

      use fpcomm
      use plinit
      use plparm
      use equnit_mod
      use fpinit
      use fpmenu
      use fpwrin
      use libmtx

      implicit none
      integer:: IERR

      CALL mtx_initialize
      IF(nrank.EQ.0) THEN
         WRITE(6,*) '***** TASK/FP 2009/09/18 *****'
         CALL GSOPEN
         OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
      ENDIF
      ierr_g=0
      N_f1=0

      CALL pl_init
      CALL eq_init
      CALL fp_init
      IF(nrank.EQ.0) THEN
         CALL pl_parm(1,'plparm',IERR)
         CALL eqparm(1,'eqparm',IERR)
         CALL fp_parm(1,'fpparm',IERR)
      ENDIF
      CALL fp_menu

      IF(nrank.EQ.0) CALL GSCLOS
      CALL mtx_finalize
      CALL fp_wr_deallocate
      STOP
      END PROGRAM fp
