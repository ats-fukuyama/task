program fow

  use fowcomm
  use fowprep
  use foworbit
  use fowcoef
  use fowdistribution
  use fowloop

  use fpwrite

  use plinit
  use plparm
  use equnit_mod
  use fpparm
  use fpmenu
  use fpcomm
  use fpprep
  use fpsub

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


  ! if ( model_mkcsv /= 0 ) call fow_csv(nrmax/2)

  write(*,*)"end"
  ! call fow_deallocate

end program fow

subroutine fow_csv(nr_out)
  use fowcomm
  use fowdistribution
  use fpcomm
  use fpwrite

  integer,intent(in) :: nr_out
  real(rkind),allocatable :: gamma_r(:), Deff(:), rmg(:)
  allocate(gamma_r(nrmax))
  allocate(Deff(nrmax))
  allocate(rmg(nrmax))

  do nr = 1, nrmax
    rmg(nr) = rg(nr+1)
  end do
  

  call radial_particle_flux(gamma_r,2)
  call effective_diffusion_cosfficient(Deff, 2)

  call fpcsv1D(gamma_r,"./csv/gamma.csv")
  call fpcsv1D(Deff,"./csv/Deff.csv")

  call fpcsv1D(psimg,"./csv/psimg.csv")
  call fpcsv1D(Fpsig,"./csv/Fpsig.csv")
  call fpcsv1D(Boutg,"./csv/Boutg.csv")
  call fpcsv1D(Bing,"./csv/Bing.csv")
  call fpcsv1D(psim,"./csv/psim.csv")
  call fpcsv1D(Fpsi,"./csv/Fpsi.csv")
  call fpcsv1D(Bout,"./csv/Bout.csv")
  call fpcsv1D(Bin,"./csv/Bin.csv")
  call fpcsv1D(pm(:,1),"./csv/pm_ele.csv")
  call fpcsv1D(pg(:,1),"./csv/pg_ele.csv")
  call fpcsv1D(pm(:,2),"./csv/pm_ion.csv")
  call fpcsv1D(pg(:,2),"./csv/pg_ion.csv")
  call fpcsv1D(rm,"./csv/rm.csv")
  call fpcsv1D(rg,"./csv/rg.csv")
  call fpcsv1D(rmg,"./csv/rmg.csv")


  call fpcsv2D(thetam(:,:,nr_out,2),"./csv/thetam.csv")
  call fpcsv2D(JI(:,:,nr_out,2),"./csv/Jacobian.csv")

  call fpcsv2D(thetam(:,:,nr_out,1),"./csv/thetam_ele.csv")
  call fpcsv2D(JI(:,:,nr_out,1),"./csv/Jacobian_ele.csv")

  call fpcsv2D(theta_pnc(:,:,1),"./csv/theta_pnc_ele.csv")
  call fpcsv2D(theta_pnc(:,:,2),"./csv/theta_pnc_ion.csv")
  call fpcsv2D(theta_co_stg(:,:,1),"./csv/theta_co_ele.csv")
  call fpcsv2D(theta_co_stg(:,:,2),"./csv/theta_co_ion.csv")
  call fpcsv2D(theta_cnt_stg(:,:,1),"./csv/theta_cnt_ele.csv")
  call fpcsv2D(theta_cnt_stg(:,:,2),"./csv/theta_cnt_ion.csv")


  call fpcsv2D(fnsp(:,:,nr_out,2),"./csv/fnsp.csv")
  call fpcsv1D(fnsp(npmax/2,nthmax/4,:,2),"./csv/fnsp_radial_th4.csv")
  call fpcsv1D(fnsp(npmax/2,nthmax/3,:,2),"./csv/fnsp_radial_th3.csv")
  call fpcsv1D(fnsp(npmax/2,nthmax/2,:,2),"./csv/fnsp_radial_th2.csv")

  call fpcsv2D(fnsm(:,:,nr_out,2),"./csv/fnsm.csv")
  call fpcsv1D(fnsm(npmax/2,nthmax/4,:,2),"./csv/fnsm_radial_th4.csv")
  call fpcsv1D(fnsm(npmax/2,nthmax/3,:,2),"./csv/fnsm_radial_th3.csv")
  call fpcsv1D(fnsm(npmax/2,nthmax/2,:,2),"./csv/fnsm_radial_th2.csv")

end subroutine

