module fowloop
  private
  public :: fow_loop 

  double precision,allocatable,dimension(:,:,:)     :: density, p_flux, Deff, temperature
  double precision,allocatable,dimension(:,:,:)     :: Sr, Sr_Dp, Sr_Dth, Sr_Dr, Sr_F
  double precision,allocatable,dimension(:,:,:,:,:) :: fI, source
  double precision,allocatable,dimension(:,:,:,:) :: nu_collision

  double precision,allocatable,dimension(:,:)       :: BMT, DLT, FMT
  double precision,allocatable,dimension(:,:,:)     :: ALT

contains

  subroutine fow_loop
    use fowcomm
    use fowexec
    use fowsource
    use fowcoef
    use fpcomm
    use fowdistribution
    use fpsub
    use fpprep

    use fpwrite
    
    implicit none
    integer :: nt, nth, np, nr, nsa, n_iterate, ierr = 0, its
    real(rkind) :: deps, sumF0, sumFd
    logical :: iteration_flag, coef_is_updated
    real(rkind) :: begin_time, end_time, begin_time_loop, end_time_loop

    coef_is_updated = .false.

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
                ! write(*,*)fnsp(nth,np,nr,nsa),fns0(nth,np,nr,nsa)
                sumFd = sumFd + ( fnsp(nth,np,nr,nsa)-fns0(nth,np,nr,nsa) )**2
                sumF0 = sumF0 + fns0(nth,np,nr,nsa)**2
              end do
            end do
          end do
          deps = MAX( deps, sumFd/sumF0 )
        end do
        write(*,*)"deps,it",deps,n_iterate

        if ( deps <= epsfp .or. n_iterate >= lmaxfp ) iteration_flag = .false.

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

      end do ! end of do while

      if ( .not.coef_is_updated ) then
        call update_bulk_temperature
        call coulomb_log
        call fow_coef
        call fow_calculate_source(nt)
        coef_is_updated = .false.
      end if

      ! ! update density and temperature
      call moment_0th_order_COM(rnsl, fnsp)
      call moment_2nd_order_COM(rwsl, fnsp)

      call output_data(nt)

      call cpu_time(end_time_loop)
      write(6,'(A,I0,A,ES10.3,A)'),'time to loop(nt=',nt,'):',end_time_loop-begin_time_loop,'[sec]'
      
    end do

  end subroutine fow_loop

  subroutine update_bulk_temperature
    use fpcomm
    use fowcomm
    implicit none
    integer :: nth, np, nr, nsa, ns
    real(rkind) :: dVI, energy, pv, rns_bulk, rws_bulk, J2keV

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
            dVI = JIR(nth,np,nr,nsa)*delp(ns)*delthm(nth,np,nr,nsa)
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
    implicit none
    integer,intent(in) :: nt
    double precision :: r_, psip0, costh0, sinth0, B0, F0, dBdr0, dFdr0, dpsipdr0, MJ2keV
    double precision,dimension(nrmax) :: rmg
    double precision,dimension(nthmax,npmax,nrmax,nsamax) :: taup, r0, gammaI
    integer :: nth,np,nr,nsa,nm,nm2

    if ( .not.allocated(density) ) allocate(density(             nrmax,nsamax,ntmax))
    if ( .not.allocated(p_flux ) ) allocate(p_flux (             nrmax,nsamax,ntmax))
    if ( .not.allocated(Sr     ) ) allocate(Sr     (             nrmax,nsamax,ntmax))
    if ( .not.allocated(Sr_Dp  ) ) allocate(Sr_Dp  (             nrmax,nsamax,ntmax))
    if ( .not.allocated(Sr_Dth ) ) allocate(Sr_Dth (             nrmax,nsamax,ntmax))
    if ( .not.allocated(Sr_Dr  ) ) allocate(Sr_Dr  (             nrmax,nsamax,ntmax))
    if ( .not.allocated(Sr_F   ) ) allocate(Sr_F   (             nrmax,nsamax,ntmax))
    if ( .not.allocated(Deff   ) ) allocate(Deff   (             nrmax,nsamax,ntmax))
    if ( .not.allocated(fI     ) ) allocate(fI     (nthmax,npmax,nrmax,nsamax,ntmax))
    if ( .not.allocated(source ) ) allocate(source (nthmax,npmax,nrmax,nsamax,ntmax))
    if ( .not.allocated(BMT    ) ) allocate(BMT    (nmmax,ntmax                    ))
    if ( .not.allocated(DLT    ) ) allocate(DLT    (nmmax,ntmax                    ))
    if ( .not.allocated(FMT    ) ) allocate(FMT    (nmmax,ntmax                    ))
    if ( .not.allocated(ALT    ) ) allocate(ALT    (nmmax,nlmaxm,ntmax             ))

    if ( .not.allocated(temperature) ) allocate(temperature(nrmax,nsamax,ntmax))
    if ( .not.allocated(nu_collision) ) allocate(nu_collision(nrmax,nsbmax,nsamax,ntmax))
    

    MJ2keV = 1.d-3/aee*1.d6

    if ( nt == 1 ) then
      do nr = 1, nrmax
        rmg(nr) = rg(nr+1)
      end do
    
      do nsa = 1, nsamax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              call mean_ra_quantities(orbit_m(nth,np,nr,nsa), r_, psip0, costh0, sinth0, B0, F0, dBdr0, dFdr0, dpsipdr0)
              r0(nth,np,nr,nsa) = r_
              taup(nth,np,nr,nsa) = orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max)
            end do
          end do
        end do
      end do

      call system('mkdir -p dat')
    
      call fptxt1D(psimg,"dat/psimg.txt")
      call fptxt1D(Fpsig,"dat/Fpsig.txt")
      call fptxt1D(Boutg,"dat/Boutg.txt")
      call fptxt1D(Bing,"dat/Bing.txt")
      call fptxt1D(psim,"dat/psim.txt")
      call fptxt1D(Fpsi,"dat/Fpsi.txt")
      call fptxt1D(Bout,"dat/Bout.txt")
      call fptxt1D(Bin,"dat/Bin.txt")
      call fptxt1D(rm,"dat/rm.txt")
      call fptxt1D(rg,"dat/rg.txt")
      call fptxt1D(rmg,"dat/rmg.txt")
    
      call fptxt2D(pm,"dat/pm.txt")
      call fptxt2D(pg,"dat/pg.txt")
    
      call fptxt4D(thetam,"dat/thetam.txt")
      call fptxt4D(thetamg,"dat/thetamg.txt")
      call fptxt4D(thetam_pg,"dat/thetam_pg.txt")
      call fptxt4D(thetam_rg,"dat/thetam_rg.txt")
      call fptxt4D(JI,"dat/Jacobian.txt")
      call fptxt4D(JIR,"dat/JIR.txt")
      call fptxt4D(taup,"dat/taup.txt")
      call fptxt4D(r0,"dat/r_mean.txt")
    
      call fptxt4D(Dppfow,"dat/Dpp.txt")
      call fptxt4D(Dptfow,"dat/Dpt.txt")
      call fptxt4D(Dprfow,"dat/Dpr.txt")
    
      call fptxt4D(Dtpfow,"dat/Dtp.txt")
      call fptxt4D(Dttfow,"dat/Dtt.txt")
      call fptxt4D(Dtrfow,"dat/Dtr.txt")
    
      call fptxt4D(Drpfow,"dat/Drp.txt")
      call fptxt4D(Drtfow,"dat/Drt.txt")
      call fptxt4D(Drrfow,"dat/Drr.txt") 
      
      call fptxt4D(Fppfow,"dat/Fpp.txt")
      call fptxt4D(Fthfow,"dat/Fth.txt")
      call fptxt4D(Frrfow,"dat/Frr.txt")

      call fptxt4D(Dcpp,"dat/Dpp_fp.txt")
      call fptxt4D(Dcpt,"dat/Dpt_fp.txt")
    
      call fptxt4D(Dctp,"dat/Dtp_fp.txt")
      call fptxt4D(Dctt,"dat/Dtt_fp.txt")

      call fptxt4D(Fcpp,"dat/Fpp_fp.txt")
      call fptxt4D(Fcth,"dat/Fth_fp.txt")

    end if

    do nm = 1, nmmax
      BMT(nm,nt) = BM(nm)
      DLT(nm,nt) = DL(nm)
      FMT(nm,nt) = FM(nm)
      do nm2 = 1, nlmaxm
        ALT(nm,nm2,nt) = AL(nm,nm2)
      end do
    end do

    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            fI(nth,np,nr,nsa,nt)     = fnsp(nth,np,nr,nsa)
            source(nth,np,nr,nsa,nt) = SPP(nth,np,nr,nsa)
            if ( fnsp(nth,np,nr,nsa) /= fnsp(nth,np,nr,nsa) ) then
              write(6,*)"ERROR",nth,np,nr,nsa
              STOP
            end if
          end do
        end do
      end do
    end do

    call moment_0th_order_COM(density(:,:,nt), fnsp)
    call moment_2nd_order_COM(temperature(:,:,nt), fnsp)

    do nsa = 1, nsamax
      do nr = 1, nrmax
        temperature(nr,nsa,nt) = temperature(nr,nsa,nt)*MJ2keV/(1.5d0*density(nr,nsa,nt)*1.d20)
      end do
    end do

    call particle_flux_element(Sr(:,:,nt), Sr_Dp(:,:,nt), Sr_Dth(:,:,nt), Sr_Dr(:,:,nt), Sr_F(:,:,nt))
    call effective_diffusion_cosfficient(Deff(:,:,nt))
    call collision_frequency(nu_collision(:,:,:,nt),.false.)

    if ( nt == ntmax ) then
      call fptxt5D(fI,"dat/f.txt")
      call fptxt5D(source,"dat/source.txt")
      call fptxt3D(density,"dat/density.txt")
      call fptxt3D(Deff,"dat/Deff.txt")
      call fptxt3D(p_flux,"dat/p_flux.txt")
      call fptxt3D(Sr,"dat/Sr.txt")
      call fptxt3D(Sr_Dp,"dat/Sr_Dp.txt")
      call fptxt3D(Sr_Dth,"dat/Sr_Dth.txt")
      call fptxt3D(Sr_Dr,"dat/Sr_Dr.txt")
      call fptxt3D(Sr_F,"dat/Sr_F.txt")
      call fptxt3D(temperature,"dat/temperature.txt")
      call fptxt4D(nu_collision,"dat/nu_coll.txt")
      call fptxt2D(BMT, "dat/bm.txt")
      call fptxt2D(DLT, "dat/dl.txt")
      call fptxt3D(ALT, "dat/al.txt")
      call fptxt2D(FMT, "dat/fm.txt")
    end if

  end subroutine output_data

  subroutine collision_frequency(nu_coll, normalize_)
    use fpcomm
    use fowcomm

    implicit none
    real(rkind),dimension(nrmax,nsbmax,nsamax),intent(out) :: nu_coll
    logical,intent(in),optional :: normalize_

    logical :: normalize
    integer :: nth, np, nr, nsa, nsb, nssa, nssb
    real(rkind) :: mu, numerator, denominator, v0, temp, dens
    real(rkind) :: taup, taup_mean, dVI

    if ( present(normalize_) ) then
      normalize = normalize_
    else
      normalize = .true.
    end if

    do nsa = 1, nsamax
      nssa = ns_nsb(nsa)
      do nsb = 1, nsbmax
        nssb = ns_nsa(nsb)
        mu = pa(nssa)*pa(nssb)/(pa(nssa)+pa(nssb))*amp
        do nr = 1, nrmax
          dens = rnsl(nr,nsb)*1.d20
          temp = rwsl(nr,nsa)/( 1.5d0*rnsl(nr,nsa)*1.d20 ) ! [MJ]
          temp = temp * 1.d6                               ! [J]
          v0 = SQRT( 3.d0*temp/(pa(nssa)*amp) )            ! [m/s]

          numerator   = dens*(pz(nssa)*pz(nssb))**2*aee**4
          denominator = 2.d0*pi*eps0**2*mu**2*v0**3

          nu_coll(nr,nsb,nsa) = numerator/denominator * lnlam(nr,nsb,nsa)
        end do
      end do
    end do

    ! normalize by tau_orbit
    if ( normalize ) then

      do nsa = 1, nsamax
        nssa = ns_nsa(nsa)
        do nr = 1, nrmax
          taup_mean = 0.d0
          do np = 1, npmax
            do nth = 1, nthmax
              dVI = delp(nssa)*delthm(nth,np,nr,nsa)*JIR(nth,np,nr,nsa)
              taup = orbit_m(nth,np,nr,nsa)%time(orbit_m(nth,np,nr,nsa)%nstp_max)
              taup_mean = taup_mean + taup*fnsp(nth,np,nr,nsa)*dVI
            end do
          end do
          taup_mean = taup_mean*rnfp0(nsa)/rnsl(nr,nsa)

          do nsb = 1, nsbmax
            nu_coll(nr,nsb,nsa) = nu_coll(nr,nsb,nsa)*taup_mean
          end do

        end do
      end do
      
    end if
    

  end subroutine collision_frequency
  
end module fowloop