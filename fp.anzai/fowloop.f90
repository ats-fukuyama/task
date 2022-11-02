! fowloop.f90
! [2022/3/9]
! ****************************
!  Main loop of FOW
! ****************************
!  made by ota / modified by anzai
!

module fowloop

  private
  public :: fow_loop 

  double precision,allocatable,dimension(:,:,:) :: density, p_flux
  double precision,allocatable,dimension(:,:,:) :: Deff, temperature
  double precision,allocatable,dimension(:,:,:) :: Sr, Sr_Dp
  double precision,allocatable,dimension(:,:,:) :: Sr_Dth, Sr_Dr, Sr_F
  double precision,allocatable,dimension(:,:,:,:,:) :: fI, source

contains

  subroutine fow_loop
  !----------------------------
  ! Main loop of fow program
  !----------------------------

    use fowcomm
    use fowexec
    use fowsource
    use fowcoef
    use fpcomm
    use fowdistribution
    use fpsub
    use fpprep
    use fpwrite
    use check_neoclass
    
    implicit none
    logical :: iteration_flag, coef_is_updated
    real(rkind) :: begin_time, end_time, begin_time_loop, end_time_loop
    real(rkind) :: deps, sumF0, sumFd
    integer :: nt, nth, np, nr, nsa, n_iterate, ierr = 0, its
    
    coef_is_updated = .false.

    !**** initial distribution
    call fI_Maxwellian(fnsp)

    !**** initial coefficient

    call fow_coef    
    call fow_calculate_source(0)

    do nt = 1, ntmax
      write(*,*)"nt=",nt
      call cpu_time(begin_time_loop)

      do nsa = 1, nsamax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              fnsm(nth,np,nr,nsa) = fnsp(nth,np,nr,nsa)
            end do
          end do
        end do
      end do

      n_iterate = 0
      iteration_flag = .true.

      do while( iteration_flag )
        n_iterate = n_iterate + 1

        do nsa = 1, nsamax          
          call fow_exec(nsa, ierr, its)
        end do

        deps = 0.d0
        do nsa = 1, nsamax
          sumF0 = 0.d0
          sumFd = 0.d0
          do nr = 1, nrmax
            do np = 1, npmax
              do nth = 1, nthmax
                sumFd = sumFd + ( fnsp(nth,np,nr,nsa) &
                      - fns0(nth,np,nr,nsa) )**2
                sumF0 = sumF0 + fns0(nth,np,nr,nsa)**2
              end do
            end do
          end do
          deps = MAX( deps, sumFd/sumF0 )
        end do
        write(*,*)"deps,it",deps,n_iterate

        if ( deps <= epsfp .or. n_iterate >= lmaxfp )&
              iteration_flag = .false.

        do nsa = 1, nsamax
          do nr = 1, nrmax
            do np = 1, npmax
              do nth = 1, nthmax
                fnsp(nth,np,nr,nsa) = fns0(nth,np,nr,nsa)
              end do
            end do
          end do
        end do

        do nsa = 1, nsamax
          if ( modelc(nsa) >= 2 ) then
            call update_bulk_temperature
            call coulomb_log
            call fow_coef
            call fow_calculate_source(nt)
            coef_is_updated = .true.
            exit
          end if
        end do

      end do !** end of do while

      if ( .not.coef_is_updated ) then
        call update_bulk_temperature
        call coulomb_log
        call fow_coef
        call fow_calculate_source(nt)
        coef_is_updated = .false.
      end if

      !**** update density and temperature
      call moment_0th_order_COM(rnsl, fnsp)
      call moment_2nd_order_COM(rwsl, fnsp)

      call output_data(nt)
      !call output_orbit_classify

      call cpu_time(end_time_loop)
      write(6,'(A,I0,A,ES10.3,A)'),'time to loop(nt=',nt,'):' &
                       ,end_time_loop-begin_time_loop,'[sec]'
      
    end do !** time loop
      call fptxt1D(rm,"dat/rm.txt")

  end subroutine fow_loop

  subroutine update_bulk_temperature
  !-----------------------------------------
  ! Subroutine of time evolution of bulk Ta
  !-----------------------------------------

    use fpcomm
    use fowcomm
    implicit none
    real(rkind) :: dVI, energy, pv, rns_bulk, rws_bulk, J2keV
    integer :: nth, np, nr, nsa, ns

    J2keV = 1.d-3/aee
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax

        rns_bulk = 0.d0
        rws_bulk = 0.d0

        do np = 1, npmax
          if ( pm(np,ns) > fact_bulk ) exit
          pv = sqrt(1.d0+theta0(ns)*pm(np,ns)**2)
          do nth = 1, nthmax
            energy = amfp(ns)*vc**2*(pv-1.d0)
            dVI = JI(nth,np,nr,nsa)*delp(ns)*delthm(nth,np,nr,nsa)
            rns_bulk = rns_bulk + fnsp(nth,np,nr,nsa)*dVI
            rws_bulk = rws_bulk + energy*fnsp(nth,np,nr,nsa)*dVI
          end do
        end do

        rt_bulk(nr,ns) = rws_bulk*J2keV / (1.5d0*rns_bulk)

      end do
    end do

    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        rn_temp(nr,ns) = rnsl(nr,ns)
        rt_temp(nr,ns) = rt_bulk(nr,ns)
      end do
    end do

  end subroutine update_bulk_temperature

  subroutine output_data(nt)
    use fowcomm
    use fpcomm
    use fpwrite
    use foworbit
    use fowdistribution
    use check_neoclass
    use orbit_classify
    implicit none
    integer,intent(in) :: nt
    double precision,dimension(nrmax) :: rmg
    double precision,dimension(npmax,nsamax) :: momm
    double precision,dimension(npmax+1,nsamax) :: momg
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: taup, r0
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: gammaI
    double precision :: r_, psip0, costh0, sinth0, B0
    double precision :: F0, dBdr0, dFdr0, dpsipdr0, MJ2keV
    integer :: nth,np,nr,nsa

    ! if ( .not.allocated(density) ) allocate(density(nrmax,nsamax,ntmax))
    ! if ( .not.allocated(p_flux ) ) allocate(p_flux (nrmax,nsamax,ntmax))
    ! if ( .not.allocated(Sr     ) ) allocate(Sr     (nrmax,nsamax,ntmax))
    ! if ( .not.allocated(Sr_Dp  ) ) allocate(Sr_Dp  (nrmax,nsamax,ntmax))
    ! if ( .not.allocated(Sr_Dth ) ) allocate(Sr_Dth (nrmax,nsamax,ntmax))
    ! if ( .not.allocated(Sr_Dr  ) ) allocate(Sr_Dr  (nrmax,nsamax,ntmax))
    ! if ( .not.allocated(Sr_F   ) ) allocate(Sr_F   (nrmax,nsamax,ntmax))
    ! if ( .not.allocated(Deff   ) ) allocate(Deff   (nrmax,nsamax,ntmax))
    ! if ( .not.allocated(fI     ) )  &
    !   allocate(fI(nthmax,npmax,nrmax,nsamax,ntmax))
    ! if ( .not.allocated(source ) )  &
    !   allocate(source (nthmax,npmax,nrmax,nsamax,ntmax))

    ! if ( .not.allocated(temperature)) &
    !   allocate(temperature(nrmax,nsamax,ntmax))
    

    ! MJ2keV = 1.d-3/aee*1.d6

    !  if ( nt == 1 ) then
    !  if ( nt == 0 ) then !** modified by anzai
    !   do nr = 1, nrmax
    !     rmg(nr) = rg(nr+1)
    !   end do

    !   do nsa = 1, nsamax
    !     do np = 1, npmax
    !       momm(np,nsa) = pm(np,nsa)*ptfp0(nsa)
    !       momg(np,nsa) = pg(np,nsa)*ptfp0(nsa)
    !     end do
    !     momg(npmax+1,nsa) = pg(npmax+1,nsa)*ptfp0(nsa)
    !   end do
    
      ! do nsa = 1, nsamax
      !   do nr = 1, nrmax
      !     do np = 1, npmax
      !       do nth = 1, nthmax
      !         call mean_ra_quantities(orbit_m(nth,np,nr,nsa), &
      !                r_, psip0, costh0, sinth0, B0, F0, dBdr0,&
      !                dFdr0, dpsipdr0)
      !         r0(nth,np,nr,nsa) = r_
      !         taup(nth,np,nr,nsa) = &
      !         orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max)
      !       end do
      !     end do
      !   end do
      ! end do
      
      !***** make directory
      call system('mkdir -p dat')

      !***** txt file output
      ! call fptxt1D(psimg,"dat/psimg.txt")
      ! call fptxt1D(psim,"dat/psim.txt")
      ! call fptxt1D(Fpsig,"dat/Fpsig.txt")
      ! call fptxt1D(Boutg,"dat/Boutg.txt")
      ! call fptxt1D(Bing,"dat/Bing.txt")
      ! call fptxt1D(Bout,"dat/Bout.txt")
      ! call fptxt1D(Bin,"dat/Bin.txt")
      ! call fptxt1D(Fpsi,"dat/Fpsi.txt")
      call fptxt1D(rm,"dat/rm.txt")
      ! call fptxt1D(rg,"dat/rg.txt")
      ! call fptxt1D(rmg,"dat/rmg.txt")

    !   call fptxt1D(safety_factor,"dat/safety_factor.txt")
    
    !   call fptxt2D(pm,"dat/pm.txt")
    !   call fptxt2D(pg,"dat/pg.txt")
    !   call fptxt2D(momm,"dat/momentum.txt")
    !   call fptxt2D(momg,"dat/momentumg.txt")
    
    !   call fptxt4D(thetam,"dat/thetam.txt")
    !   call fptxt4D(thetamg,"dat/thetamg.txt")
    !   call fptxt4D(thetam_pg,"dat/thetam_pg.txt")
    !   call fptxt4D(thetam_rg,"dat/thetam_rg.txt")
    !   call fptxt4D(JI,"dat/Jacobian.txt")
    !   call fptxt4D(JI,"dat/JI.txt")
    !   call fptxt4D(taup,"dat/taup.txt")
    !   call fptxt4D(r0,"dat/r_mean.txt")
    
    !   call fptxt4D(Dppfow,"dat/Dpp.txt")
    !   call fptxt4D(Dptfow,"dat/Dpt.txt")
    !   call fptxt4D(Dprfow,"dat/Dpr.txt")
    
    !   call fptxt4D(Dtpfow,"dat/Dtp.txt")
    !   call fptxt4D(Dttfow,"dat/Dtt.txt")
    !   call fptxt4D(Dtrfow,"dat/Dtr.txt")
    
    !   call fptxt4D(Drpfow,"dat/Drp.txt")
    !   call fptxt4D(Drtfow,"dat/Drt.txt")
    !   call fptxt4D(Drrfow,"dat/Drr.txt") 
      
    !   call fptxt4D(Fppfow,"dat/Fpp.txt")
    !   call fptxt4D(Fthfow,"dat/Fth.txt")
    !   call fptxt4D(Frrfow,"dat/Frr.txt")

    !   call fptxt4D(Dcpp,"dat/Dpp_fp.txt")
    !   call fptxt4D(Dcpt,"dat/Dpt_fp.txt")
    
    !   call fptxt4D(Dctp,"dat/Dtp_fp.txt")
    !   call fptxt4D(Dctt,"dat/Dtt_fp.txt")

    !   call fptxt4D(Fcpp,"dat/Fpp_fp.txt")
    !   call fptxt4D(Fcth,"dat/Fth_fp.txt")

    !  end if

    ! do nsa = 1, nsamax
    !   do nr = 1, nrmax
    !     do np = 1, npmax
    !       do nth = 1, nthmax
    !         fI(nth,np,nr,nsa,nt)     = fnsp(nth,np,nr,nsa)
    !         source(nth,np,nr,nsa,nt) = SPP(nth,np,nr,nsa)
    !         if ( fnsp(nth,np,nr,nsa) /= fnsp(nth,np,nr,nsa) ) then
    !           write(6,*)"ERROR",nth,np,nr,nsa
    !           STOP
    !         end if
    !       end do
    !     end do
    !   end do
    ! end do

    ! call moment_0th_order_COM(density(:,:,nt), fnsp)
    ! call moment_2nd_order_COM(temperature(:,:,nt), fnsp)

    ! do nsa = 1, nsamax
    !   do nr = 1, nrmax
    !     temperature(nr,nsa,nt) = temperature(nr,nsa,nt) &
    !                            * MJ2keV/(1.5d0*density(nr,nsa,nt)*1.d20)
    !   end do
    ! end do

    ! call particle_flux_element(Sr(:,:,nt), Sr_Dp(:,:,nt),&
    !        Sr_Dth(:,:,nt), Sr_Dr(:,:,nt), Sr_F(:,:,nt))
    ! call effective_diffusion_cosfficient(Deff(:,:,nt))

    ! if ( nt == ntmax ) then
    !   call fptxt5D(fI,"dat/f.txt")
    !   call fptxt5D(source,"dat/source.txt")
    !   call fptxt3D(density,"dat/density.txt")
    !   call fptxt3D(Deff,"dat/Deff.txt")
    !   call fptxt3D(p_flux,"dat/p_flux.txt")
    !   call fptxt3D(Sr,"dat/Sr.txt")
    !   call fptxt3D(Sr_Dp,"dat/Sr_Dp.txt")
    !   call fptxt3D(Sr_Dth,"dat/Sr_Dth.txt")
    !   call fptxt3D(Sr_Dr,"dat/Sr_Dr.txt")
    !   call fptxt3D(Sr_F,"dat/Sr_F.txt")
    !   call fptxt3D(temperature,"dat/temperature.txt")
    ! end if

  end subroutine output_data

end module fowloop