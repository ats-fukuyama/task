program fow

  use eqfp
  use fow_global
  use fpwrite

  use fpinit
  use fpmenu
  use fpcomm
  use fpprep
  use libmtx

  implicit none

  integer :: IERR

  call fow_read_namelist
  call fow_allocate
  call interfaceEQandFP(psi,psirz,Fpsi,Frz,Bout,Bin,Brz,psi0)
  call fow_prep

  call mtx_initialize
  call fp_init
  call fp_parm(1,'fpparm',IERR)
  call fp_broadcast
  call fp_prep(IERR)


  ! call fpcsv1D(psi,"./csv/psi.csv")
  ! call fpcsv1D(F,"./csv/F.csv")
  ! call fpcsv2D(Frz,"./csv/Frz.csv")
  ! call fpcsv2D(psirz,"./csv/psirz.csv")
  ! call fpcsv2D(Brz,"./csv/Brz.csv")
  ! call fpcsv1D(Bout,"./csv/Bout.csv")
  ! call fpcsv1D(Bin,"./csv/Bin.csv")

  call fow_deallocate

  write(*,*)"end"

end program fow
