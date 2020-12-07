module fowloop
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

    use fpwrite
    
    implicit none
    integer :: nt, nth, np, nr, nsa, n_iterate, ierr = 0, its
    real(rkind) :: deps, sumF0, sumFd
    logical :: iteration_flag

    call update_quantities
    call fow_coef
    call fow_calculate_source  
    open(109,file="./txt/fow_FNSP.txt")
    do nsa = 1, nsamax
      do nr = 1, nrmax
        do np = 1, npmax
          do nth = 1, nthmax
            write(109,*)fnsp(nth,np,nr,nsa)
          end do
        end do
      end do
    end do
    close(109)

    call fpcsv1D(rnsl(:,2),"./csv/n0.csv")

    do nt = 1, ntmax
      write(*,*)"nt=",nt

      do nsa = 1, nsamax
        do nr = 1, nrmax
          do np = 1, npmax
            do nth = 1, nthmax
              fnsm(nth,np,nr,nsa)=fnsp(nth,np,nr,nsa)
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
        ! lmaxfp = 1
        if ( deps <= epsfp .or. n_iterate >= lmaxfp ) iteration_flag = .false.

        do nsa = 1, nsamax
          do nr = 1, nrmax
            do np = 1, npmax
              do nth = 1, nthmax
                ! write(*,'(A,4I3,2ES12.4)')"F",nth,np,nr,nsa,fnsp(nth,np,nr,nsa),fns0(nth,np,nr,nsa)
                fnsp(nth,np,nr,nsa) = fns0(nth,np,nr,nsa)
              end do
            end do
          end do
        end do

        call update_quantities

        call fow_coef

      end do ! end of do while
      
    end do

  end subroutine fow_loop

  subroutine update_quantities
    use plprof
    use fpcomm
    use fowcomm
    use fowdistribution
    use foworbit

    implicit none
    type(pl_plf_type),dimension(nsmax):: plf
    integer :: nth, np, nr, nsa, ns, nstp,nstpmax
    real(rkind) :: rhon, Ntot, mean_ra

    ! do nsa = 1, nsamax
    !   fnorm(nsa) = 1.d-40*RNFP0(NSA)
    ! end do

    ! update rn_temp, rt_temp
    do nsa = 1, nsamax
      ns = ns_nsa(nsa)
      do nr = 1, nrmax
        rhon = rm(nr)
        call pl_prof(rhon,plf)
        rn_temp(nr,ns) = plf(ns)%rn
        rt_temp(nr,ns) = (plf(ns)%rtpr+2.d0*plf(ns)%rtpp)/3.d0
      end do
    end do

    ! ! update density
    call moment_0th_order_COM(RNSL, fnsp)

    ! ! update temperature
    ! call moment_2nd_order_COM(RWSL, fnsp)

    ! update bulk temperature

    ! update Coulomb logarithm

    ! -----------------s

    do nsa = 1, nsamax
      do nr = 1, nrmax
            write(*,*)"RNSL",RNSL(nr,nsa)
      end do
    end do

    ! do nsa = 1, nsamax
    !   do nr = 1, nrmax
    !         write(*,*)"RWSL",RWSL(nr,nsa)/aee*1.D3
    !   end do
    ! end do

    call total_N(Ntot, fnsp, 1)
    write(*,*)"electron total", Ntot
    call total_N(Ntot, fnsp, 2)
    write(*,*)"ion total", Ntot

  end subroutine update_quantities

end module fowloop