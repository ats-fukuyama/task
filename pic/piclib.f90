MODULE piclib

  PRIVATE
  PUBLIC kine,pote,sumdim,sumdim1

  CONTAINS

!***********************************************************************
    subroutine kine(npmax,vx,vy,vz,akin,mass)
!***********************************************************************
      implicit none
      real(8), dimension(npmax) :: vx, vy, vz
      real(8) :: akin, mass
      integer(4) :: npmax, np
      akin = 0.d0
      do np = 1, npmax
         akin = akin + vx(np)**2 + vy(np)**2 + vz(np)**2
      end do

      akin = 0.5 * akin * mass /dble(npmax)
    end subroutine kine

!***********************************************************************
    subroutine pote(nxmax,nymax,ex,ey,ez,bx,by,bz,vcfact,apot)
!***********************************************************************
      implicit none
      real(8), dimension(0:nxmax,0:nymax) :: ex, ey, ez, bx, by, bz
      real(8) :: apot, vcfact
      integer(4) :: nxmax, nymax, nx, ny
      apot = 0.d0

      do ny = 0, nymax-1
      do nx = 0, nxmax-1

         apot = apot + ex(nx,ny)**2 + ey(nx,ny)**2 + ez(nx,ny)**2 &
                     + bx(nx,ny)**2 + by(nx,ny)**2 + bz(nx,ny)**2
            
      end do
      end do
      apot = 0.5 * apot / (dble(nxmax)*dble(nymax))

    end subroutine pote

!***********************************************************************
    subroutine sumdim(nodes,myid,a,b,ndim)
!***********************************************************************
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      integer stat1(mpi_status_size)
      integer stat2(mpi_status_size)
      dimension a(ndim), b(ndim)
!***********************************************************************
      kmod = 1 
!.....................................................
      if(nodes.gt.1) then
!.....................................................
 
    1 continue
 
      kmod1 = kmod
      kmod = kmod*2
 
      if (mod(myid,kmod) .lt. kmod1) then
        idiff = kmod1
      else
        idiff = -kmod1
      endif
 
      call mpi_isend(a,ndim,mpi_real8,myid+idiff,300,  &
                            mpi_comm_world,ireq1,ierr)
      call mpi_irecv(b,ndim,mpi_real8,myid+idiff,300,  &
                            mpi_comm_world,ireq2,ierr)

      call mpi_wait(ireq1,stat1,ierr)
      call mpi_wait(ireq2,stat2,ierr)

       do j = 1, ndim
       a(j) = a(j) + b(j)
       end do

      if (kmod .lt. nodes) goto 1
 
!.....................................................
      endif
!.....................................................
 
    end subroutine sumdim

!***********************************************************************
    subroutine sumdim1(nodes,myid,a1,b1)
!***********************************************************************
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      integer stat1(mpi_status_size)
      integer stat2(mpi_status_size)
      dimension a(1), b(1)
      a(1)=a1
      b(1)=b1
      kmod = 1 
!.....................................................
      if(nodes.gt.1) then
!.....................................................
 
    1 continue
 
      kmod1 = kmod
      kmod = kmod*2
 
      if (mod(myid,kmod) .lt. kmod1) then
        idiff = kmod1
      else
        idiff = -kmod1
      endif
 
      call mpi_isend(a,ndim,mpi_real8,myid+idiff,300,  &
                            mpi_comm_world,ireq1,ierr)
      call mpi_irecv(b,ndim,mpi_real8,myid+idiff,300,  &
                            mpi_comm_world,ireq2,ierr)

      call mpi_wait(ireq1,stat1,ierr)
      call mpi_wait(ireq2,stat2,ierr)

       do j = 1, ndim
       a(j) = a(j) + b(j)
       end do

      if (kmod .lt. nodes) goto 1
 
      endif
!.....................................................
 
    end subroutine sumdim1

END Module piclib
