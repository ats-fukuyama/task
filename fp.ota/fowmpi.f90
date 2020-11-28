subroutine fow_save_orbit_mpi(ierr)
  use libmtx
  use libmpi
  use fpcomm
  use fowcomm
  implicit none
  integer,intent(out) :: ierr
  integer :: err, err_recv
  character(30) :: BIN_DIR, filename
  integer :: nm ,nth, np, nr, nsa, nstp, nstpmax, dummy
  integer :: nmax_m, nmax_t, nmax_p, nmax_r, nmax, nmax_ml, nmax_tl, nmax_pl, nmax_rl
  integer,allocatable :: nstpm_m(:), nstpm_t(:), nstpm_p(:), nstpm_r(:), sum_m(:), sum_t(:), sum_p(:), sum_r(:)
  integer,allocatable :: nstpm_ml(:), nstpm_tl(:), nstpm_pl(:), nstpm_rl(:)
  INTEGER,DIMENSION(nsize) :: ilena,iposa

  call mtx_reset_communicator

  nmax_m = npmax*nthmax*nrmax*nsamax
  allocate(nstpm_m(nmax),sum_m(nmax+1))
  nmax_t = npmax*(nthmax+1)*nrmax*nsamax
  allocate(nstpm_t(nmax),sum_t(nmax+1))
  nmax_p = (npmax+1)*nthmax*nrmax*nsamax
  allocate(nstpm_p(nmax),sum_p(nmax+1))
  nmax_r = npmax*nthmax*(nrmax+1)*nsamax
  allocate(nstpm_r(nmax),sum_r(nmax+1))

  nmax_ml = nthmax*(npend-npstart+1)*(nrend-nrstart+1)*(nsaend-nsastart+1)
  allocate(nstpm_ml(nmax))
  nmax_tl = (nthmax+1)*(npend-npstart+1)*(nrend-nrstart+1)*(nsaend-nsastart+1)
  allocate(nstpm_tl(nmax))
  nmax_pl = nthmax*(npendwg-npstart+1)*(nrend-nrstart+1)*(nsaend-nsastart+1)
  allocate(nstpm_pl(nmax))
  nmax_rl = nthmax*(npend-npstart+1)*(nrendwg-nrstart+1)*(nsaend-nsastart+1)
  allocate(nstpm_rl(nmax))

  BIN_DIR = "../fp.ota/bin/"

  filename = TRIM(BIN_DIR)//"eqparm.bin"
  open(10,file=filename,access='direct',recl=rkind,form='unformatted',status='replace',iostat=err)
  write(10,rec=1)RR
  write(10,rec=2)RA
  write(10,rec=3)RKAP
  write(10,rec=4)RDLT
  write(10,rec=5)RB
  write(10,rec=6)BB
  write(10,rec=7)RIP
  close(10)
  call mtx_allreduce1_integer(err**2,3,err_recv,dummy)
  if ( err_recv /= 0 ) then
    ierr = err
    return
  end if


  filename = TRIM(BIN_DIR)//"fpparm.bin"
  open(11,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=err)
  write(11,rec=1)nthmax
  write(11,rec=2)npmax
  write(11,rec=3)nrmax
  write(11,rec=4)nsamax
  close(11)
  call mtx_allreduce1_integer(err**2,3,err_recv,dummy)
  if ( err_recv /= 0 ) then
    ierr = err
    return
  end if


  filename = TRIM(BIN_DIR)//"ob_rec_m.bin"
  open(40,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_rec_t.bin"
  open(41,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_rec_p.bin"
  open(42,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_rec_r.bin"
  open(43,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=err)


  filename = TRIM(BIN_DIR)//"ob_nstpmax_m.bin"
  open(50,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_nstpmax_t.bin"
  open(51,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_nstpmax_p.bin"
  open(52,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_nstpmax_r.bin"
  open(53,file=filename,access='direct',recl=4,form='unformatted',status='replace',iostat=err)

  call mtx_allreduce1_integer(err**2,3,err_recv,dummy)
  if ( err_recv /= 0 ) then
    ierr = err
    return
  end if

  
  do nsa = nsastart, nsaend
    do nr = nrstart, nrendwg
      do np = npstart, npendwg
        do nth = 1, nthmax+1

          if ( nth <= nthmax .and. np <= npmax .and. nr <= nrmax ) then
            nm = (nsa-1)*nthmax*(npend-npstart+1)*(nrend-nrstart+1)&
                +(nr-1)*nthmax*(npend-npstart+1)&
                +(np-1)*nthmax + nth
            write(50,rec=nm)orbit_m(nth,np,nr,nsa)%nstp_max
            nstpm_ml(nm) = orbit_m(nth,np,nr,nsa)%nstp_max

          end if

          if ( np <= npmax .and. nr <= nrmax ) then
            nm = (nsa-1)*(nthmax+1)*(npend-npstart+1)*(nrend-nrstart+1)&
                +(nr-1)*(nthmax+1)*(npend-npstart+1)&
                +(np-1)*(nthmax+1) + nth

            write(51,rec=nm)orbit_th(nth,np,nr,nsa)%nstp_max
            nstpm_tl(nm) = orbit_th(nth,np,nr,nsa)%nstp_max

          end if
  
          if ( nth <= nthmax .and. nr <= nrmax ) then
            nm = (nsa-1)*nthmax*(npendwg-npstart+1)*(nrend-nrstart+1)&
                +(nr-1)*nthmax*(npendwg-npstart+1)&
                +(np-1)*nthmax + nth

            write(52,rec=nm)orbit_p(nth,np,nr,nsa)%nstp_max
            nstpm_pl(nm) = orbit_p(nth,np,nr,nsa)%nstp_max

          end if

          if ( nth <= nthmax .and. np <= npmax ) then
            nm = (nsa-1)*nthmax*(npend-npstart+1)*(nrendwg-nrstart+1)&
                +(nr-1)*nthmax*(npend-npstart+1)&
                +(np-1)*nthmax + nth

            write(53,rec=nm)orbit_r(nth,np,nr,nsa)%nstp_max
            nstpm_rl(nm) = orbit_r(nth,np,nr,nsa)%nstp_max

          end if

        end do
      end do
    end do
  end do

  ilena(nrank-1) = nmax_ml
  iposa(nrank-1) = (nsastart-1)*nthmax*npmax*nrmax+&
                  +(nrstart-1)*nthmax*npmax+(npstart-1)*nthmax+nth
  call mtx_gatherv_integer(nstpm_ml,nmax_ml,nstpm_m,nmax_m,ilena,iposa)

  ilena(nrank-1) = nmax_pl
  iposa(nrank-1) = (nsastart-1)*nthmax*(npmax+1)*nrmax+&
                  +(nrstart-1)*nthmax*(npmax+1)+(npstart-1)*nthmax+nth
  call mtx_gatherv_integer(nstpm_pl,nmax_pl,nstpm_p,nmax_p,ilena,iposa)

  ilena(nrank-1) = nmax_tl
  iposa(nrank-1) = (nsastart-1)*(nthmax+1)*npmax*nrmax+&
                  +(nrstart-1)*(nthmax+1)*npmax+(npstart-1)*(nthmax+1)+nth
  call mtx_gatherv_integer(nstpm_tl,nmax_tl,nstpm_t,nmax_t,ilena,iposa)

  ilena(nrank-1) = nmax_rl
  iposa(nrank-1) = (nsastart-1)*nthmax*npmax*(nrmax+1)&
                  +(nrstart-1)*nthmax*npmax+(npstart-1)*nthmax+nth
  call mtx_gatherv_integer(nstpm_rl,nmax_rl,nstpm_r,nmax_r,ilena,iposa)

  if ( nrank == 0 ) then
    sum_m(1)=1
    sum_t(1)=1
    sum_p(1)=1
    sum_r(1)=1

    nmax = MAX( nmax_m,nmax_t,nmax_p,nmax_r )
    do nm = 1, nmax
      if ( nm <= nmax_m ) then
        write(40,rec=nm,iostat=err)sum_m(nm)
        sum_m(nm+1) = sum_m(nm)+nstpm_m(nm)*5
      end if
      if ( nm <= nmax_t ) then
        write(41,rec=nm,iostat=err)sum_t(nm)
        sum_t(nm+1) = sum_t(nm)+nstpm_t(nm)*5
      end if
      if ( nm <= nmax_p ) then
        write(42,rec=nm,iostat=err)sum_p(nm)
        sum_p(nm+1) = sum_p(nm)+nstpm_p(nm)*5
      end if
      if ( nm <= nmax_r ) then
        write(43,rec=nm,iostat=err)sum_r(nm)
        sum_r(nm+1) = sum_r(nm)+nstpm_r(nm)*5
      end if
    end do

  end if
  
  call mtx_broadcast_integer(sum_m,nmax_m)
  call mtx_broadcast_integer(sum_t,nmax_t)
  call mtx_broadcast_integer(sum_p,nmax_p)
  call mtx_broadcast_integer(sum_t,nmax_r)

  call mtx_allreduce1_integer(err**2,3,err_recv,dummy)
  if ( err_recv /= 0 ) then
    ierr = err
    return
  end if

  filename = TRIM(BIN_DIR)//"ob_m.bin"
  open(60,file=filename,access='direct',recl=rkind,form='unformatted',status='replace',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_t.bin"
  open(61,file=filename,access='direct',recl=rkind,form='unformatted',status='replace',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_p.bin"
  open(62,file=filename,access='direct',recl=rkind,form='unformatted',status='replace',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_r.bin"
  open(63,file=filename,access='direct',recl=rkind,form='unformatted',status='replace',iostat=err)

  call mtx_allreduce1_integer(err**2,3,err_recv,dummy)
  if ( err_recv /= 0 ) then
    ierr = err
    return
  end if

  do nsa = nsastart, nsaend
    do nr = nrstart, nrendwg
      do np = npstart, npendwg
        do nth = 1, nthmax+1

          if ( nth <= nthmax .and. np <= npmax .and. nr <= nrmax ) then
            nm = (nsa-1)*(nthmax+1)*npmax*nrmax&
            +(nr-1)*(nthmax+1)*npmax&
            +(np-1)*(nthmax+1) + nth

            nstpmax = orbit_th(nth,np,nr,nsa)%nstp_max
            do nstp = 1, nstpmax
              write(60,rec=sum_m(nm))orbit_m(nth,np,nr,nsa)%time(nstp)
              write(60,rec=sum_m(nm)+1)orbit_m(nth,np,nr,nsa)%psip(nstp)
              write(60,rec=sum_m(nm)+2)orbit_m(nth,np,nr,nsa)%theta(nstp)
              write(60,rec=sum_m(nm)+3)orbit_m(nth,np,nr,nsa)%thetap(nstp)
              write(60,rec=sum_m(nm)+4)orbit_m(nth,np,nr,nsa)%Babs(nstp)
            end do             
          end if

          if ( np <= npmax .and. nr <= nrmax ) then
            nm = (nsa-1)*(nthmax+1)*npmax*nrmax&
                +(nr-1)*(nthmax+1)*npmax&
                +(np-1)*(nthmax+1) + nth

            nstpmax = orbit_th(nth,np,nr,nsa)%nstp_max
            do nstp = 1, nstpmax
              write(61,rec=sum_t(nm))orbit_th(nth,np,nr,nsa)%time(nstp)
              write(61,rec=sum_t(nm)+1)orbit_th(nth,np,nr,nsa)%psip(nstp)
              write(61,rec=sum_t(nm)+2)orbit_th(nth,np,nr,nsa)%theta(nstp)
              write(61,rec=sum_t(nm)+3)orbit_th(nth,np,nr,nsa)%thetap(nstp)
              write(61,rec=sum_t(nm)+4)orbit_th(nth,np,nr,nsa)%Babs(nstp)
            end do
          end if
  
          if ( nth <= nthmax .and. nr <= nrmax ) then
            nm = (nsa-1)*nthmax*(npmax+1)*nrmax&
                +(nr-1)*nthmax*(npmax+1)&
                +(np-1)*nthmax + nth

            nstpmax = orbit_p(nth,np,nr,nsa)%nstp_max
            do nstp = 1, nstpmax
              write(62,rec=sum_p(nm))orbit_p(nth,np,nr,nsa)%time(nstp)
              write(62,rec=sum_p(nm)+1)orbit_p(nth,np,nr,nsa)%psip(nstp)
              write(62,rec=sum_p(nm)+2)orbit_p(nth,np,nr,nsa)%theta(nstp)
              write(62,rec=sum_p(nm)+3)orbit_p(nth,np,nr,nsa)%thetap(nstp)
              write(62,rec=sum_p(nm)+4)orbit_p(nth,np,nr,nsa)%Babs(nstp)
            end do
          end if

          if ( nth <= nthmax .and. np <= npmax ) then
            nm = (nsa-1)*nthmax*npmax*(nrmax+1)&
                +(nr-1)*nthmax*npmax&
                +(np-1)*nthmax + nth

            nstpmax = orbit_r(nth,np,nr,nsa)%nstp_max
            do nstp = 1, nstpmax
              write(63,rec=sum_r(nm))orbit_r(nth,np,nr,nsa)%time(nstp)
              write(63,rec=sum_r(nm)+1)orbit_r(nth,np,nr,nsa)%psip(nstp)
              write(63,rec=sum_r(nm)+2)orbit_r(nth,np,nr,nsa)%theta(nstp)
              write(63,rec=sum_r(nm)+3)orbit_r(nth,np,nr,nsa)%thetap(nstp)
              write(63,rec=sum_r(nm)+4)orbit_r(nth,np,nr,nsa)%Babs(nstp)
            end do
          end if

          call mtx_allreduce1_integer(err**2,3,err_recv,dummy)
          if ( err_recv /= 0 ) then
            ierr = err
            return
          end if  

        end do      
      end do
    end do
  end do

  do np = 0, 3
    close(40+np)
    close(50+np)
    close(60+np)
  end do

end subroutine fow_save_orbit_mpi

subroutine fow_load_orbit_mpi(ierr)
  use libmtx
  use libmpi
  use fpcomm
  use fowcomm
  implicit none
  integer,intent(out) :: ierr
  integer :: err, err_recv, dummy, record
  character(30) :: BIN_DIR, eqin, mesh, filename
  integer :: nm ,nth, np, nr, nsa, nstp, nstpmax

  real(rkind) :: RR_,RA_,RKAP_,RDLT_,RB_,BB_,RIP_
  integer :: nthm_,npm_,nrm_,nsam_,access

  ierr = 0
  err = 0
  err_recv = 0

  BIN_DIR = "./bin/"

  if ( access( TRIM(BIN_DIR)//"fpparm.bin", " ") /= 0 &
      .and. access( TRIM(BIN_DIR)//"eqparm.bin", " ") /= 0 ) then
    ierr = 1
    return
  end if

  filename = TRIM(BIN_DIR)//"eqparm.bin"
  open(10, file=filename,access='direct',recl=rkind,form='unformatted',status='old',iostat=err)
  read(10,rec=1,iostat=err)RR_
  read(10,rec=2,iostat=err)RA_
  read(10,rec=3,iostat=err)RKAP_
  read(10,rec=4,iostat=err)RDLT_
  read(10,rec=5,iostat=err)RB_
  read(10,rec=6,iostat=err)BB_
  read(10,rec=7,iostat=err)RIP_
  close(10)

  filename = TRIM(BIN_DIR)//"fpparm.bin"
  open(11, file=filename,access='direct',recl=4,form='unformatted',status='old',iostat=err)
  read(11,rec=1,iostat=err)nthm_
  read(11,rec=2,iostat=err)npm_
  read(11,rec=3,iostat=err)nrm_
  read(11,rec=4,iostat=err)nsam_
  close(11)
  call mtx_allreduce1_integer(err**2,3,err_recv,dummy)
  if ( err_recv /= 0 ) then
    ierr = err
    return
  end if


  if ( RR /= RR_ .or. RA /= RA_ .or. RKAP /= RKAP_ &
      .or. RDLT /= RDLT_ .or. RB /= RB_ .or. BB /= BB_ .or. RIP /= RIP_) then
    ierr = 2
    write(*,*)"Equilibrium parameters do not match with binary data."
    return
  end if
  if ( nthmax /= nthm_ .or. npmax /= npm_ &
      .or. nrmax /= nrm_ .or. nsamax /= nsam_ ) then
    ierr = 3
    write(*,*)"Mesh width does not match with binary data."
    return    
  end if
  

  filename = TRIM(BIN_DIR)//"ob_nstpmax_m.bin"
  open(50,file=filename,access='direct',recl=4,form='unformatted',status='old',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_nstpmax_t.bin"
  open(51,file=filename,access='direct',recl=4,form='unformatted',status='old',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_nstpmax_p.bin"
  open(52,file=filename,access='direct',recl=4,form='unformatted',status='old',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_nstpmax_r.bin"
  open(53,file=filename,access='direct',recl=4,form='unformatted',status='old',iostat=err)

  filename = TRIM(BIN_DIR)//"ob_m.bin"
  open(60,file=filename,access='direct',recl=rkind,form='unformatted',status='old',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_t.bin"
  open(61,file=filename,access='direct',recl=rkind,form='unformatted',status='old',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_p.bin"
  open(62,file=filename,access='direct',recl=rkind,form='unformatted',status='old',iostat=err)
  filename = TRIM(BIN_DIR)//"ob_r.bin"
  open(63,file=filename,access='direct',recl=rkind,form='unformatted',status='old',iostat=err)

  call mtx_allreduce1_integer(err**2,3,err_recv,dummy)
  if ( err_recv /= 0 ) then
    ierr = err
    return
  end if

  do nsa = nsastart, nsaend
    do nr = nrstart, nrendwg
      do np = npstart, npendwg
        do nth = 1, nthmax
          if ( nth <= nthmax .and. np <= npmax .and. nr <= nrmax ) then
            nm = (nsa-1)*nthmax*npmax*nrmax&
                +(nr-1)*nthmax*npmax&
                +(np-1)*nthmax + nth
            read(50,rec=nm,iostat=err)nstpmax
            orbit_m(nth,np,nr,nsa)%nstp_max = nstpmax

            allocate(orbit_m(nth,np,nr,nsa)%time(nstpmax))
            allocate(orbit_m(nth,np,nr,nsa)%psip(nstpmax))
            allocate(orbit_m(nth,np,nr,nsa)%Babs(nstpmax))
            allocate(orbit_m(nth,np,nr,nsa)%theta(nstpmax))
            allocate(orbit_m(nth,np,nr,nsa)%thetap(nstpmax))
            
            read(40,rec=nm,iostat=err)record
            do nstp = 1, nstpmax
              read(60,rec=record,iostat=err)orbit_m(nth,np,nr,nsa)%time(nstp)
              read(60,rec=record+1,iostat=err)orbit_m(nth,np,nr,nsa)%psip(nstp)
              read(60,rec=record+2,iostat=err)orbit_m(nth,np,nr,nsa)%theta(nstp)
              read(60,rec=record+3,iostat=err)orbit_m(nth,np,nr,nsa)%thetap(nstp)
              read(60,rec=record+4,iostat=err)orbit_m(nth,np,nr,nsa)%Babs(nstp)
            end do
            
          end if

          if ( np <= npmax .and. nr <= nrmax ) then
            nm = (nsa-1)*(nthmax+1)*npmax*nrmax&
                +(nr-1)*(nthmax+1)*npmax&
                +(np-1)*(nthmax+1) + nth
            read(51,rec=nm,iostat=err)nstpmax
            orbit_th(nth,np,nr,nsa)%nstp_max = nstpmax

            allocate(orbit_th(nth,np,nr,nsa)%time(nstpmax))
            allocate(orbit_th(nth,np,nr,nsa)%psip(nstpmax))
            allocate(orbit_th(nth,np,nr,nsa)%Babs(nstpmax))
            allocate(orbit_th(nth,np,nr,nsa)%theta(nstpmax))
            allocate(orbit_th(nth,np,nr,nsa)%thetap(nstpmax))

            read(41,rec=nm,iostat=err)record
            do nstp = 1, nstpmax
              read(61,rec=record,iostat=err)orbit_th(nth,np,nr,nsa)%time(nstp)
              read(61,rec=record+1,iostat=err)orbit_th(nth,np,nr,nsa)%psip(nstp)
              read(61,rec=record+2,iostat=err)orbit_th(nth,np,nr,nsa)%theta(nstp)
              read(61,rec=record+3,iostat=err)orbit_th(nth,np,nr,nsa)%thetap(nstp)
              read(61,rec=record+4,iostat=err)orbit_th(nth,np,nr,nsa)%Babs(nstp)
            end do
          end if
  
          if ( nth <= nthmax .and. nr <= nrmax ) then
            nm = (nsa-1)*nthmax*(npmax+1)*nrmax&
                +(nr-1)*nthmax*(npmax+1)&
                +(np-1)*nthmax + nth

            read(52,rec=nm,iostat=err)nstpmax
            orbit_p(nth,np,nr,nsa)%nstp_max = nstpmax

            allocate(orbit_p(nth,np,nr,nsa)%time(nstpmax))
            allocate(orbit_p(nth,np,nr,nsa)%psip(nstpmax))
            allocate(orbit_p(nth,np,nr,nsa)%Babs(nstpmax))
            allocate(orbit_p(nth,np,nr,nsa)%theta(nstpmax))
            allocate(orbit_p(nth,np,nr,nsa)%thetap(nstpmax))

            read(42,rec=nm,iostat=err)record
            do nstp = 1, nstpmax
              read(62,rec=record,iostat=err)orbit_p(nth,np,nr,nsa)%time(nstp)
              read(62,rec=record+1,iostat=err)orbit_p(nth,np,nr,nsa)%psip(nstp)
              read(62,rec=record+2,iostat=err)orbit_p(nth,np,nr,nsa)%theta(nstp)
              read(62,rec=record+3,iostat=err)orbit_p(nth,np,nr,nsa)%thetap(nstp)
              read(62,rec=record+4,iostat=err)orbit_p(nth,np,nr,nsa)%Babs(nstp)
            end do
          end if

          if ( nth <= nthmax .and. np <= npmax ) then
            nm = (nsa-1)*nthmax*npmax*(nrmax+1)&
                +(nr-1)*nthmax*npmax&
                +(np-1)*nthmax + nth

            read(53,rec=nm)nstpmax
            orbit_r(nth,np,nr,nsa)%nstp_max = nstpmax

            allocate(orbit_r(nth,np,nr,nsa)%time(nstpmax))
            allocate(orbit_r(nth,np,nr,nsa)%psip(nstpmax))
            allocate(orbit_r(nth,np,nr,nsa)%Babs(nstpmax))
            allocate(orbit_r(nth,np,nr,nsa)%theta(nstpmax))
            allocate(orbit_r(nth,np,nr,nsa)%thetap(nstpmax))

            read(43,rec=nm,iostat=err)record
            do nstp = 1, nstpmax
              read(63,rec=record,iostat=err)orbit_r(nth,np,nr,nsa)%time(nstp)
              read(63,rec=record+1,iostat=err)orbit_r(nth,np,nr,nsa)%psip(nstp)
              read(63,rec=record+2,iostat=err)orbit_r(nth,np,nr,nsa)%theta(nstp)
              read(63,rec=record+3,iostat=err)orbit_r(nth,np,nr,nsa)%thetap(nstp)
              read(63,rec=record+4,iostat=err)orbit_r(nth,np,nr,nsa)%Babs(nstp)
            end do
          end if

        end do
      end do
    end do
  end do

  call mtx_allreduce1_integer(err**2,3,err_recv,dummy)
  if ( err_recv /= 0 ) then
    ierr = err
    return
  end if

  do np = 0, 3
    close(40+np)
    close(50+np)
    close(60+np)
  end do

end subroutine fow_load_orbit_mpi

