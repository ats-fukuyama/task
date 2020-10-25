program fow

  use fowcomm
  use fowprep
  use foworbit
  use foworbitclassify
  use fowdistribution

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

  call fow_orbit_construct(orbit_m)
  call fow_orbit_construct(orbit_th)

  write(*,*)"end"
  call fow_deallocate

end program fow

subroutine fow_debug

  use fowcomm
  use foworbitclassify
  use fowdistribution

  use fpcomm
  use fpwrite

  implicit none

  integer :: IERR=0,nr,nth,np,nsa,nstp
  integer :: mode(3), nr_out
  real(rkind),allocatable :: check_orbit(:,:),lorentz(:),lorentzg(:),beta(:),betag(:),mean_r(:,:,:,:)&
                            ,trapped_boundary(:,:,:), forbitten_boundary(:,:,:), X_boundary(:,:),dummy(:),dummy2(:)
  real(rkind) :: mean_psip, dummy3, dummy4

  if ( nrmax==50 ) then
    nr_out=33
  else
    nr_out=nrmax/2+3
  end if

  write(*,*)"mean"
  allocate(mean_r(nthmax,npmax,nrmax,nsamax))
  do nsa = 1, nsamax
    do nr = 1, nrmax
      do np = 1, npmax
        do nth = 1, nthmax
          if ( forbitten(nth,np,nr,nsa,[0,0,0]) ) then
            mean_r(nth,np,nr,nsa) = 0.d0
          else

            if ( size(orbit_m(nth,np,nr,nsa)%psip)==0 ) then
              write(*,*)"STOP at fowdebug"
              write(*,*)nr,"size(orbit_m(nth,np,nr,nsa)%psip)==0"
              STOP
            end if
 
            mean_psip = sum(orbit_m(nth,np,nr,nsa)%psip)/size(orbit_m(nth,np,nr,nsa)%psip)
            call fow_get_ra_from_psip(mean_r(nth,np,nr,nsa), mean_psip)
            mean_r(nth,np,nr,nsa) = rm(nr)-mean_r(nth,np,nr,nsa)
          end if
        end do
      end do
    end do
  end do

  call fpcsv2D(mean_r(:,:,2,2),"./csv/mean_r_center.csv")
  call fpcsv2D(mean_r(:,:,nr_out,2),"./csv/mean_r_quarter.csv")
  call fpcsv2D(mean_r(:,:,nrmax,2),"./csv/mean_r_edge.csv")

  write(*,*)"equiv_variable"
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
  call fpcsv1D(xi,"./csv/xi.csv")
  call fpcsv1D(xig,"./csv/xig.csv")

  allocate(dummy(nrmax+1))
  call first_order_derivative(dummy, Boutg, psimg)
  call fpcsv1D(dummy,"./csv/dBoutgdpsi.csv")
  call second_order_derivative(dummy, Boutg, psimg)
  call fpcsv1D(dummy,"./csv/d2Boutgdpsi.csv")
  deallocate(dummy)

  allocate(trapped_boundary(nthmax,nrmax,nsamax),forbitten_boundary(nthmax,nrmax,nsamax))
  call fow_trapped_boundary(trapped_boundary)
  call fow_forbitten_boundary(forbitten_boundary)

  ! do nsa = 1, nsamax
  !   do nr = 1, nrmax
  !     do nth = 1, nthmax
  !       if ( trapped_boundary(nth,nr,nsa) >= pm(npmax,nsa) ) then
  !         trapped_boundary(nth,nr,nsa) = pm(npmax,nsa)
  !       end if
  !       if ( forbitten_boundary(nth,nr,nsa) >= pm(npmax,nsa) ) then
  !         forbitten_boundary(nth,nr,nsa) = pm(npmax,nsa)
  !       end if
  !     end do
  !   end do  
  ! end do

  call fpcsv2D(trapped_boundary(:,:,2),"./csv/trapped_boundary.csv")
  call fpcsv1D(trapped_boundary(:,2,2),"./csv/trapped_boundary_center.csv")
  call fpcsv1D(trapped_boundary(:,nr_out,2),"./csv/trapped_boundary_quarter.csv")
  call fpcsv1D(trapped_boundary(:,nrmax,2),"./csv/trapped_boundary_edge.csv")

  call fpcsv2D(forbitten_boundary(:,:,2),"./csv/forbitten_boundary.csv")
  call fpcsv1D(forbitten_boundary(:,2,2),"./csv/forbitten_boundary_center.csv")
  call fpcsv1D(forbitten_boundary(:,nr_out,2),"./csv/forbitten_boundary_quarter.csv")
  call fpcsv1D(forbitten_boundary(:,nrmax,2),"./csv/forbitten_boundary_edge.csv")

  allocate(X_boundary(nrmax,2))
  call fow_stagnation_type(X_boundary)
  ! nrmax = 50 -> nr = 33, r/a = 0.65 gaii

  allocate(dummy(nthmax),dummy2(nthmax))
  do nth = 1, nthmax
    dummy(nth) = 1.d0 - ( 1.d0-X_boundary(nr_out,1) )*(nth-1)/(nthmax-1)
  end do
  call fpcsv1D(dummy,"./csv/xi_x_stagnation_co.csv")
  do nth = 1, nthmax
    call fow_stagnation_orbit_velocity(dummy2(nth),dummy(nth),rm(nr_out),2)
  end do
  call fpcsv1D(dummy2,"./csv/beta_x_stagnation_co.csv")

  do nth = 1, nthmax
    dummy(nth) = X_boundary(nr_out,1) - ( X_boundary(nr_out,1)-0.d0)*(nth-1)/(nthmax-1)
  end do
  call fpcsv1D(dummy,"./csv/xi_o_stagnation_co.csv")
  do nth = 1, nthmax
    call fow_stagnation_orbit_velocity(dummy2(nth),dummy(nth),rm(nr_out),2)
  end do
  call fpcsv1D(dummy2,"./csv/beta_o_stagnation_co.csv")

  do nth = 1, nthmax
    dummy(nth) = 0.d0 - ( 0.d0-X_boundary(nr_out,2) )*(nth-1)/(nthmax-1)
  end do
  call fpcsv1D(dummy,"./csv/xi_x_stagnation_cnt.csv")
  do nth = 1, nthmax
    call fow_stagnation_orbit_velocity(dummy2(nth),dummy(nth),rm(nr_out),2)
  end do
  call fpcsv1D(dummy2,"./csv/beta_x_stagnation_cnt.csv")

  do nth = 1, nthmax
    dummy(nth) = X_boundary(nr_out,2) - ( X_boundary(nr_out,2)+1.d0 )*(nth-1)/(nthmax-1)
  end do
  call fpcsv1D(dummy,"./csv/xi_o_stagnation_cnt.csv")
  do nth = 1, nthmax
    call fow_stagnation_orbit_velocity(dummy2(nth),dummy(nth),rm(nr_out),2)
  end do
  call fpcsv1D(dummy2,"./csv/beta_o_stagnation_cnt.csv")

  deallocate(dummy,dummy2)
  allocate(dummy(nr_out), dummy2(nr_out))
  call  fow_pinch_orbit(dummy, dummy2, nr_out, 2)
  call fpcsv1D(dummy,"./csv/beta_pinch.csv")
  call fpcsv1D(dummy2,"./csv/xi_pinch.csv")

  write(*,*)"E", amfp(2)*vc**2*((1.d0-forbitten_boundary(nthmax*7/8,nr_out,2)**2)**(-0.5d0)-1.d0)/aee*1.0d-3
  write(*,*)"psin", psim(nr_out)/psi0
  write(*,*)"cospicth", xi(nthmax*7/8)

end subroutine
