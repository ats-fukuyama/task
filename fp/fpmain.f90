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

      N_f1=0
         OPEN(9,file="f1_1.dat")
         open(10,file='time_evol.dat') 
         open(11,file='efield_e1.dat') 
         open(12,file='dndt.dat') 
         open(13,file='radial.dat') 
         open(14,file='nth-re.dat')
        open(15,file='re_pitch.dat')
        open(18,file='efield_ref.dat')
      CALL pl_init
      CALL eq_init
      CALL fp_init
      IF(nrank.EQ.0) THEN
         CALL pl_parm(1,'plparm',IERR)
         CALL eqparm(1,'eqparm',IERR)
         CALL fp_parm(1,'fpparm',IERR)
      ENDIF
      CALL fp_menu

         close(9)
         close(10)
         close(11)
         close(12)
         close(13)
         close(14)
         close(15)
         close(18)
      IF(nrank.EQ.0) CALL GSCLOS
      CALL mtx_finalize
      CALL fp_wr_deallocate
      STOP
      END PROGRAM fp
