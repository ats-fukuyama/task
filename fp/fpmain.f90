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
      use fpwrin
      use libmtx

      implicit none
      integer:: IERR

      CALL mtx_initialize
      write(6,*) 'nsize,nrank,ncomm=',nsize,nrank,ncomm
      IF(nrank.EQ.0) THEN
         WRITE(6,*) '***** TASK/FP 2015/02/08 *****'
         CALL GSOPEN
         OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
      ENDIF

      N_f1=0
!         OPEN(9,file="f1_1.dat")
!      IF(model_disrupt.ne.0)THEN
!         open(10,file='time_evol.dat') 
!         open(11,file='efield_e1.dat') 
!         open(12,file='dndt.dat') 
!         open(13,file='radial.dat') 
!         open(14,file='nth-re.dat')
!         open(15,file='re_pitch.dat')
!         OPEN(16,file="FNS01.dat")
!         OPEN(17,file="FNS03.dat")
!         OPEN(18,file="FNS05.dat")
!      END IF

      CALL pl_init
      CALL eq_init
      CALL fp_init
      IF(nrank.EQ.0) THEN
         CALL pl_parm(1,'plparm',IERR)
         CALL eqparm(1,'eqparm',IERR)
         CALL fp_parm(1,'fpparm',IERR)
      ENDIF
      CALL fp_menu

!         close(9)
 !     IF(model_disrupt.ne.0)THEN
         CLOSE(7)
!         close(10)
!         close(11)
!         close(12)
!         close(13)
!         close(14)
!         close(15)
!         close(16)
!         close(17)
!         close(18)
 !     END IF
      IF(nrank.EQ.0) CALL GSCLOS
      CALL mtx_finalize
      CALL fp_wr_deallocate
      STOP
      END PROGRAM fp
