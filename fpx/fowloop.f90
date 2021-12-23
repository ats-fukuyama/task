module fowloop
  USE fpcomm,ONLY: rkind
  
  private
  public :: fow_loop 

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
    use fpfout

    use fpwrite
    
    implicit none
    integer :: nt, nth, np, nr, nsa, n_iterate, ierr = 0, its
    real(rkind) :: deps, sumF0, sumFd
    logical :: iteration_flag, coef_is_updated
    real(rkind) :: begin_time_loop, end_time_loop

    coef_is_updated = .false.

    ! initial distribution
    call fI_Maxwellian(fnsp)
    ! initial coefficient
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

      IF(model_fow_fout.NE.0) call fp_fout_tot(nt)

      call cpu_time(end_time_loop)
      write(6,'(A,I0,A,ES10.3,A)'),'time to loop(nt=',nt,'):',end_time_loop-begin_time_loop,'[sec]'
      
    end do ! time loop

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

end module fowloop
