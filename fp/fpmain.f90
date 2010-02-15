!     $Id$

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
      use equnit_mod
      use fpinit
      use fpmenu
      use libmtx

      implicit none
      integer:: IERR

      CALL mtx_initialize(nrank,nprocs)
      IF(nrank.EQ.0) THEN
         WRITE(6,*) '***** TASK/FP 2009/09/18 *****'
         CALL GSOPEN
         OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
      ENDIF

      CALL pl_init
      CALL eq_init
      CALL fp_init
      IF(nrank.EQ.0) THEN
         CALL pl_parm(1,'plparm',IERR)
         CALL eqparm(1,'eqparm',IERR)
         CALL fp_parm(1,'fpparm',IERR)
      ENDIF
      CALL fp_menu

      CLOSE(7)
      IF(nrank.EQ.0) CALL GSCLOS
      CALL mtx_finalize
      STOP
      END PROGRAM fp
