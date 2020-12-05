!     Test program of libmtxc for 1D Shroedinger equation
!     coded by Y.maruyama

PROGRAM testmtxc

  USE task_kinds,ONLY: dp
  USE libmpi
  USE commpi
  USE libmtx
  USE libgrf,only:grd1d

  IMPLICIT NONE
  integer :: isize,itype,idata(2)
  integer :: istart,iend,its
  integer :: imax,jwidth
  integer :: ndiv,i
  integer :: MODE
  real(dp) :: tolerance,ddata(1)
  real(dp) :: xmin,xmax,dt,dx
  real :: cputime1,cputime2
  complex(dp),dimension(:),pointer:: x
  complex(dp) :: k
  character :: character*1
  real(dp),dimension(:)  ,pointer :: FX
  real(dp),dimension(:,:),pointer :: FY
  character ::STR*80

  call mtx_initialize
  if(nrank.eq.0) then
     write(6,*) 'nrank, nsize = ',nrank,nsize
     call gsopen
  endif

  itype = 0
  tolerance=1.d-7
 
  dt   = 0.05
  xmin = 0.d0
  xmax = 1.d0
  ndiv = 101

1 continue
  if(nrank.eq.0) then
     write(6,'(A)') "#TEST PROGRAM : Schroedinger equation solver"
     write(6,"(A,/E12.4,1X,2I4,D12.4)")&
          &  "#INPUT: dt,ndiv,itype=",&
          &           dt,ndiv,itype
     read(5,*,ERR=1)  dt,ndiv,itype

     isize = ndiv

     idata(1)= isize
     idata(2)= itype
     ddata(1)= tolerance

  end if

  call mtx_broadcast_integer(idata,2)
  call mtx_broadcast_real8(ddata,1)

  isize     = idata(1)
  itype     = idata(2)
  tolerance = ddata(1)

  dx = (xmax-xmin)/DBLE(ndiv-1)
  k = dt/((0.d0, 1.d0)*(dx**2))

  if(nrank.eq.0) call cpu_time(cputime1)

  imax = isize
  jwidth = 3
  allocate(x(imax))

  do i=1,imax
     x(i)=exp(-20.0D0*((i-1)*dx-0.5D0)**2)
  end do

  call mtxc_setup(imax,istart,iend)

  do i=istart,iend
     if(i.gt.1) call mtxc_set_matrix(i,i-1,k)
     call mtxc_set_matrix(i,i,1-2*k)
     if(i.lt.isize) call mtxc_set_matrix(i,i+1,k)
  end do

  do i=istart,iend
     call mtxc_set_source(i,x(i))
  end do

  call mtxc_solve(itype,tolerance,its)
  if(nrank.eq.0) then
     write(6,*) 'Iteration Number=',its
  end if

  call mtxc_gather_vector(x)

2 continue

  if(nrank.eq.0) then

     allocate(FX(isize),FY(isize,2))
     do i=1,isize
        FX(i) = dx*(i-1)
     end do
        do i=1,isize
           FY(i,1) = real(x(i))
           FY(i,2) = aimag(x(i))
           write(6,'(i5,1p3e12.4)') i,FX(i),FY(i,1),FY(i,2)
        end do
     STR="@test1@"
     MODE = 0
     call pages
     call GRD1D(0,FX,FY,isize,isize,2,STR,MODE)
     call pagee
     deallocate(FX,FY)

3    continue
     write(6,*) "#INPUT: (C)CONTINUE,(Q)QUIT"
     read (5,*,ERR=3) character
  endif
  CALL mtx_broadcast_character(character,1)

  call mtxc_cleanup
  if (character.eq."c")then
     go to 1
  else if(character.ne."q")then
     go to 2
  end if
  
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
