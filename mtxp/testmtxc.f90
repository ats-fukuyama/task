!     Test program of libmtx for 1D/2D/3D Poisson equation

!     Input parameters
!        idimen : number of dimension (i or 2 or 3),  0 for quit
!        isiz : number of mesh point in one dimension
!        isource : source position (isource, isource, isource)
!        itype: type of linear solver (0 for default)
!        tolerance : tolerance in iterative method

program testmtxc

  use libmtxc
  implicit none
  integer :: nrank,nprocs
  integer :: isize,itype,idata(2)
  integer :: istart,iend,its
  integer :: imax,jwidth,jsource
  integer :: ndiv,i,j
  integer :: MODE
  real(8) :: tolerance,ddata(1)
  real(8) :: xmin,xmax,dt,dx
  real(4) :: cputime1,cputime2
  complex(8),dimension(:),pointer:: x
  complex(8) :: k,cdata(1)
  character :: character*1
  real(8),dimension(:)  ,pointer :: FX
  real(8),dimension(:,:),pointer :: FY
  character ::STR*80

  call mtx_initialize(nrank,nprocs)

  itype = 0
  tolerance=1.d-7
 
  dt   = 0.001
  xmin = 0.d0
  xmax = 1.d0
  ndiv = 100

1 continue

  if(nrank.eq.0) then
     write(6,*) "########## test program : Schroedinger equation solver ##############"
     write(6,*) "                    coded Y.Maruyama                                 "   
2    write(6,"(A,/E12.4,1X,2I4,D12.4)")&
             &  "#INPUT: dt,ndiv,itype=",&
               &         dt,ndiv,itype
     read(5,*,ERR=2) dt,ndiv,itype

     dx = (xmax-xmin)/DBLE(ndiv)

     k = dt/((0.d0, 1.d0)*(dx**2))

     isize = ndiv

     idata(1)= isize
     idata(2)= itype
     ddata(1)= tolerance
     cdata(1)= k

  end if

  call mtx_broadcast_integer(idata,2)
  call mtx_broadcast_real8(ddata,1)
  call mtx_broadcast_complex8(cdata,1)

  isize     = idata(1)
  itype     = idata(2)
  tolerance = ddata(1)
  k         = cdata(1)

  if(nrank.eq.0) call cpu_time(cputime1)

  !first order solver
  imax = isize
  jwidth = 3
  allocate(x(imax))

  do i=1,imax
     x(i)=exp(-0.1*(REAL(i)*dx-0.5)**2)
  end do

  call mtx_setup(imax,istart,iend,jwidth)

  do i=istart,iend
     if(i.gt.1) call mtx_set_matrix(i,i-1,k)
     call mtx_set_matrix(i,i,1-2*k)
     if(i.lt.isize) call mtx_set_matrix(i,i+1,k)
  end do

100 do i=istart,iend
     call mtx_set_source(i,x(i))
  end do

  call mtx_solve(itype,tolerance,its)
  if(nrank.eq.0) then
     write(6,*) 'Iteration Number=',its
  end if

  call mtx_gather_vector(x)

  if(nrank.eq.0) then
!     open(7,file='out')
!     do i=1,isize
!        write(7,'(E12.4,1X,E12.4)'),dx*i,real(x(i))
!     end do
!     close(7 )
     call gsopen

     allocate(FX(isize),FY(isize,1))
     do i=1,isize
        FX(i) = dx*i
     end do
        do i=1,isize
           FY(i,1) = real(x(i))
        end do
     STR=""
     MODE = 0
     call pages
     call GRD1D(0,FX,FY,isize,isize,1,STR,MODE)
     call pagee

3    write(6,*) "#INPUT: (C)CONTINUE,(Q)QUIT"
     read (5,*,ERR=3) character
     if (character.eq."c")then
        goto 100
     else if(character.eq."q")then
        deallocate(FX,FY)
        goto 4
     else
        go to 3
     end if
  end if

4 call mtx_cleanup
  
  deallocate(x)
  if(nrank.eq.0) then
     call cpu_time(cputime2)
     write(6,'(A,F12.3)') &
          '--cpu time =',cputime2-cputime1
     call gsclos
  end if
  call mtx_finalize

  stop
end program testmtxc
