module fowdebug
  
contains
  subroutine check_Jacobian
    use fowcomm
    use fpcomm

    implicit none 
    integer :: np, nth, nr, nsa
    real(rkind) :: V_v, V_x, J_sum, deltaps, deltath, deltap, pl, PVmax, Vmax

    do nsa = 1, nsamax
      PVmax = SQRT(1.D0+THETA0(NSA)*PM(NPMAX,NSA)**2)
      Vmax = vc * sqrt(1.d0-1.d0/PVmax**2)

      V_v = 4.d0/3.d0*pi*Vmax**3
      
      V_x = (pi*ra**2)*(2.d0*pi*rr)

      J_sum = 0.d0
      do nr = 1, nrmax
        deltaps = psimg(nr+1)-psimg(nr)
        do np = 1, npmax
          deltap = (pg(np+1,nsa)-pg(np,nsa))*ptfp0(nsa)
          do nth = 1, nthmax
            deltath = thetamg(nth+1,np,nr,nsa)-thetamg(nth,np,nr,nsa)
            J_sum = J_sum + deltap*deltath*deltaps*Jacobian_I(nth,np,nr,nsa)
          end do
        end do
      end do
      write(*,*)"nsa = ",nsa
      write(*,*)"Jacobian,V_xv,Vmax",J_sum,V_x*V_v, Vmax,PVmax

    end do

  end subroutine

  subroutine fow_debug

    use fowcomm
    use foworbitclassify
    use fowdistribution
    use fowcoef
    use fowsource
  
    use fpcomm
    use fpwrite
  
    implicit none
  
    integer :: IERR=0,nr,nth,np,nsa,nstp
    integer :: mode(3), nr_out
    real(rkind),allocatable :: check_orbit(:,:),lorentz(:),lorentzg(:),beta(:),betag(:),fi(:,:,:,:),fu(:,:,:,:),ful(:,:,:,:,:)&
                              ,trapped_boundary(:,:,:), forbitten_boundary(:,:,:), X_boundary(:,:),dummy(:),dummy2(:), D2(:,:)
    real(rkind) :: mean_psip, dummy3, dummy4
  
    if ( nrmax==50 ) then
      nr_out=33
    else
      nr_out=nrmax/3
    end if
  
    allocate(fI(nthmax,npmax,nrmax,nsamax))
    allocate(fu(nthmax,npmax,nrmax,nsamax))
    allocate(ful(nthmax,npmax,nrmax,nthpmax,nsamax))
    call fow_Maxwellian_COM(fnsp)
    allocate(D2(nrmax,nsamax))
    call moment_0th_order_COM(D2,fnsp)
    do nr = 1, nrmax
      write(*,*)"M0",D2(nr,2)
    end do
    call convert_fI_to_fu(ful, fnsp)
    call fow_coef
    call fow_calculate_source
  
  
    ! ! call fpcsv2D(fu(:,:,nr_out,2),"./csv/fu.csv")
    ! ! call fpcsv2D(fI(:,:,nr_out,2),"./csv/fI.csv")
    call fpcsv2D(thetam(:,:,nr_out,2),"./csv/thetam.csv")
    call fpcsv2D(Jacobian_I(:,:,nr_out,2),"./csv/Jacobian.csv")
  
  
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
  
    ! allocate(dummy(nrmax+1))
    ! call first_order_derivative(dummy, Boutg, psimg)
    ! call fpcsv1D(dummy,"./csv/dBoutgdpsi.csv")
    ! call second_order_derivative(dummy, Boutg, psimg)
    ! call fpcsv1D(dummy,"./csv/d2Boutgdpsi.csv")
    ! deallocate(dummy)
  
    ! ! allocate(trapped_boundary(nthmax,nrmax,nsamax),forbitten_boundary(nthmax,nrmax,nsamax))
    ! ! call fow_trapped_boundary(trapped_boundary)
    ! ! call fow_forbitten_boundary(forbitten_boundary)
  
    ! ! ! do nsa = 1, nsamax
    ! ! !   do nr = 1, nrmax
    ! ! !     do nth = 1, nthmax
    ! ! !       if ( trapped_boundary(nth,nr,nsa) >= pm(npmax,nsa) ) then
    ! ! !         trapped_boundary(nth,nr,nsa) = pm(npmax,nsa)
    ! ! !       end if
    ! ! !       if ( forbitten_boundary(nth,nr,nsa) >= pm(npmax,nsa) ) then
    ! ! !         forbitten_boundary(nth,nr,nsa) = pm(npmax,nsa)
    ! ! !       end if
    ! ! !     end do
    ! ! !   end do  
    ! ! ! end do
  
    ! ! call fpcsv2D(trapped_boundary(:,:,2),"./csv/trapped_boundary.csv")
    ! ! call fpcsv1D(trapped_boundary(:,2,2),"./csv/trapped_boundary_center.csv")
    ! ! call fpcsv1D(trapped_boundary(:,nr_out,2),"./csv/trapped_boundary_quarter.csv")
    ! ! call fpcsv1D(trapped_boundary(:,nrmax,2),"./csv/trapped_boundary_edge.csv")
  
    ! ! call fpcsv2D(forbitten_boundary(:,:,2),"./csv/forbitten_boundary.csv")
    ! ! call fpcsv1D(forbitten_boundary(:,2,2),"./csv/forbitten_boundary_center.csv")
    ! ! call fpcsv1D(forbitten_boundary(:,nr_out,2),"./csv/forbitten_boundary_quarter.csv")
    ! ! call fpcsv1D(forbitten_boundary(:,nrmax,2),"./csv/forbitten_boundary_edge.csv")
  
    ! ! allocate(X_boundary(nrmax,2))
    ! ! call fow_stagnation_type(X_boundary)
    ! ! ! nrmax = 50 -> nr = 33, r/a = 0.65 gaii
  
    ! ! allocate(dummy(nthmax),dummy2(nthmax))
    ! ! do nth = 1, nthmax
    ! !   dummy(nth) = 1.d0 - ( 1.d0-X_boundary(nr_out,1) )*(nth-1)/(nthmax-1)
    ! ! end do
    ! ! call fpcsv1D(dummy,"./csv/xi_x_stagnation_co.csv")
    ! ! do nth = 1, nthmax
    ! !   call fow_stagnation_orbit_velocity(dummy2(nth),dummy(nth),rm(nr_out),2)
    ! ! end do
    ! ! call fpcsv1D(dummy2,"./csv/beta_x_stagnation_co.csv")
  
    ! ! do nth = 1, nthmax
    ! !   dummy(nth) = X_boundary(nr_out,1) - ( X_boundary(nr_out,1)-0.d0)*(nth-1)/(nthmax-1)
    ! ! end do
    ! ! call fpcsv1D(dummy,"./csv/xi_o_stagnation_co.csv")
    ! ! do nth = 1, nthmax
    ! !   call fow_stagnation_orbit_velocity(dummy2(nth),dummy(nth),rm(nr_out),2)
    ! ! end do
    ! ! call fpcsv1D(dummy2,"./csv/beta_o_stagnation_co.csv")
  
    ! ! do nth = 1, nthmax
    ! !   dummy(nth) = 0.d0 - ( 0.d0-X_boundary(nr_out,2) )*(nth-1)/(nthmax-1)
    ! ! end do
    ! ! call fpcsv1D(dummy,"./csv/xi_x_stagnation_cnt.csv")
    ! ! do nth = 1, nthmax
    ! !   call fow_stagnation_orbit_velocity(dummy2(nth),dummy(nth),rm(nr_out),2)
    ! ! end do
    ! ! call fpcsv1D(dummy2,"./csv/beta_x_stagnation_cnt.csv")
  
    ! ! do nth = 1, nthmax
    ! !   dummy(nth) = X_boundary(nr_out,2) - ( X_boundary(nr_out,2)+1.d0 )*(nth-1)/(nthmax-1)
    ! ! end do
    ! ! call fpcsv1D(dummy,"./csv/xi_o_stagnation_cnt.csv")
    ! ! do nth = 1, nthmax
    ! !   call fow_stagnation_orbit_velocity(dummy2(nth),dummy(nth),rm(nr_out),2)
    ! ! end do
    ! ! call fpcsv1D(dummy2,"./csv/beta_o_stagnation_cnt.csv")
  
    ! ! deallocate(dummy,dummy2)
    ! ! allocate(dummy(nr_out), dummy2(nr_out))
    ! ! call  fow_pinch_orbit(dummy, dummy2, nr_out, 2)
    ! ! call fpcsv1D(dummy,"./csv/beta_pinch.csv")
    ! ! call fpcsv1D(dummy2,"./csv/xi_pinch.csv")
  
    ! ! write(*,*)"E", amfp(2)*vc**2*((1.d0-forbitten_boundary(nthmax*7/8,nr_out,2)**2)**(-0.5d0)-1.d0)/aee*1.0d-3
    ! ! write(*,*)"psin", psim(nr_out)/psi0
    ! ! write(*,*)"cospicth", xi(nthmax*7/8)
  
  end subroutine  
end module fowdebug