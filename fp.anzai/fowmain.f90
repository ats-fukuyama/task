! fowmain.f90
! [2022/3/9]
! *********************************
!  Main program of FOW
! *********************************
!  made by ota / modified by anzai
!

program fow

  use fowcomm
  use fowprep
  use fowloop

  use plinit
  use plparm
  use equnit_mod
  use fpparm
  use fpmenu
  use fpcomm
  use fpprep
  use fpsub

  use fpwrite

  use libmtx

  implicit none

  integer :: IERR=0

  call mtx_initialize

  CALL pl_init
  CALL eq_init
  CALL fp_init
  CALL pl_parm(1,'plparm',IERR)
  CALL eqparm(1,'eqparm',IERR)
  CALL fp_parm(1,'fpparm',IERR)

  call fp_broadcast
  call fp_prep(IERR)
  call fow_read_namelist
  call fow_allocate
  call fow_prep

  call fow_loop

  write(*,*)"END of Program"

  call fow_deallocate

end program fow

