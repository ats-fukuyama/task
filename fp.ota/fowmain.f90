program fow

  use foweq
  use fowcomm
  use fowprep
  use fowob
  use foworbitclassify

  use fpwrite

  use plinit
  use plparm
  use equnit_mod
  use fpinit
  use fpmenu
  use fpcomm
  use fpprep
  use libmtx

  implicit none

  integer :: IERR=0,nps,nze,np
  
  real(8),allocatable :: check_orbit(:,:),lorentz(:),lorentzg(:),beta(:),betag(:)

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
  call interfaceEQandFP(psim,psirz,Fpsi,Frz,Bout,Bin,Brz,psi0)
  call fow_prep

  call fow_ob_run

  allocate(check_orbit(nzemax,npsmax))


  do nps=1,npsmax
    do nze=1,nzemax
      do np=npmax,1,-1
        if(trapped(np,nze,nps,2,[0,0,0]))then
          check_orbit(nze,nps)=pm(np,2)
          exit
        elseif(np==1.and.(.not.trapped(np,nze,nps,2,[0,0,0])))then
          check_orbit(nze,nps)=0.d0
        end if
      end do
    end do
  end do
  call fpcsv2D(check_orbit,"./csv/check_orbit_trapped.csv")
  do nps=1,npsmax
    do nze=1,nzemax
      do np=npmax,1,-1
        if(.not.forbitten(np,nze,nps,2,[1,1,1]))then
          check_orbit(nze,nps)=pm(np,2)
          exit
        elseif(np==1.and.(forbitten(np,nze,nps,2,[1,1,1])))then
          check_orbit(nze,nps)=0.d0
        end if
      end do
    end do
  end do
  call fpcsv2D(check_orbit,"./csv/check_orbit_forbitten.csv")

  call fpcsv1D(psimg,"./csv/psimg.csv")
  call fpcsv1D(Fpsig,"./csv/Fpsig.csv")
  call fpcsv1D(Boutg,"./csv/Boutg.csv")
  call fpcsv1D(Bing,"./csv/Bing.csv")
  call fpcsv1D(psim,"./csv/psim.csv")
  call fpcsv1D(Fpsi,"./csv/Fpsi.csv")
  call fpcsv1D(Bout,"./csv/Bout.csv")
  call fpcsv1D(Bin,"./csv/Bin.csv")

  deallocate(check_orbit)
  write(*,*)npmax,nzemax,npsmax

  write(*,*)"end"
  call fow_deallocate

end program fow
