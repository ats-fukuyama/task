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

    do nt = 1, ntmax
      write(*,*)"nt=",nt

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
        ! lmaxfp = 1
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

        call update_quantities

        call fow_coef

      end do ! end of do while

      if ( model_mkcsv /= 0 ) call save_csv(nt)
      
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

  subroutine save_csv(nt)
    use fowcomm
    use fpcomm
    use fpwrite
    use fowdistribution

    implicit none

    integer,intent(in) :: nt
    character(len=30) :: filef, filen, filegamma, filedeff, filefrad4, filefrad2, filefrad34
    real(rkind),dimension(nrmax) :: gamma_r, Deff, rmg
    real(rkind),dimension(nrmax,nsamax) :: N_prev
    integer :: nr, nr_out

    nr_out=nrmax/2

    if ( nt == 1 ) then
  
      do nr = 1, nrmax
        rmg(nr) = rg(nr+1)
      end do
  
      call moment_0th_order_COM(N_prev, fnsm)
      
      call fpcsv2D(fnsm(:,:,nr_out,2),"csv/fnsp0.csv")
      call fpcsv1D(fnsm(npmax/2,nthmax/2,:,2),"csv/frad20.csv")
      call fpcsv1D(fnsm(npmax/2,nthmax/4,:,2),"csv/frad40.csv")
      call fpcsv1D(fnsm(npmax/2,nthmax/4*3,:,2),"csv/frad340.csv")    
      call fpcsv1D(N_prev(:,2),"csv/n0.csv")

      call fpcsv1D(rm,"./csv/rm.csv")
      call fpcsv1D(rg,"./csv/rg.csv")
      call fpcsv1D(rmg,"./csv/rmg.csv")
      call fpcsv2D(thetam(:,:,nr_out,2),"./csv/thetam.csv")
      call fpcsv2D(JI(:,:,nr_out,2),"./csv/Jacobian.csv")
    
      call fpcsv2D(thetam(:,:,nr_out,1),"./csv/thetam_ele.csv")
      call fpcsv2D(JI(:,:,nr_out,1),"./csv/Jacobian_ele.csv")
      call fpcsv1D(pm(:,1),"./csv/pm_ele.csv")
      call fpcsv1D(pg(:,1),"./csv/pg_ele.csv")
      call fpcsv1D(pm(:,2),"./csv/pm_ion.csv")
      call fpcsv1D(pg(:,2),"./csv/pg_ion.csv")

    end if

    call radial_particle_flux(gamma_r,2)
    call effective_diffusion_cosfficient(Deff,2)  

    write(filef,'(A,I1,A)')"csv/fnsp",nt,".csv"
    write(filen,'(A,I1,A)')"csv/n",nt,".csv"
    write(filefrad2,'(A,I1,A)')"csv/frad2",nt,".csv"
    write(filefrad4,'(A,I1,A)')"csv/frad4",nt,".csv"
    write(filefrad34,'(A,I1,A)')"csv/frad34",nt,".csv"
    write(filegamma,'(A,I1,A)')"csv/gamma",nt,".csv"
    write(filedeff,'(A,I1,A)')"csv/Deff",nt,".csv"

    call fpcsv2D(fnsp(:,:,nr_out,2),filef)
    call fpcsv1D(fnsp(npmax/2,nthmax/2,:,2),filefrad2)
    call fpcsv1D(fnsp(npmax/2,nthmax/4,:,2),filefrad4)
    call fpcsv1D(fnsp(npmax/2,nthmax/4*3,:,2),filefrad34)    
    call fpcsv1D(rnsl(:,2),filen)

    call fpcsv1D(gamma_r,filegamma)
    call fpcsv1D(Deff,filedeff)

  end subroutine 

end module fowloop