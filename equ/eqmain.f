C     $Id$
c
c   ****************************************
c   **  FREE BOUNDARY EQUILIBRIUM SOLVER  **
c   **     developed by Masafumi Azumi    **
c   ****************************************
c
      program equ
      use eqmenu_mod
      use equunit_mod
      use trunit_mod
c
      write(6,*) '## topics/equ 2006/08/24'
      call gsopen
      open(7,STATUS='SCRATCH',FORM='FORMATTED')
      call equ_init
      call tr_init
c
      call equ_parm(1,'equparm',ierr)
c
      call eqmenu
c
      close(7)
      call gsclos
      stop
      end program equ
