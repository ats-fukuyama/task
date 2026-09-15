! equmain.f90

!     $Id$
!
!   ****************************************
!   **  FREE BOUNDARY EQUILIBRIUM SOLVER  **
!   **     developed by Masafumi Azumi    **
!   ****************************************
!
      program equ
      use eqmenu_mod
      use equunit_mod
      use trunit_mod
!
      write(6,*) '## topics/equ 2006/08/24'
      call gsopen
      open(7,STATUS='SCRATCH',FORM='FORMATTED')
!
      call equ_init
      call tr_init
!
      call equ_parm(1,'equparm',ierr)
!
      call eqmenu
!
      close(7)
      call gsclos
      stop
      end program equ
