!  ***** TASK/PIC EXEC *****

Module picexec

  USE piccomm
 
  PRIVATE
  PUBLIC pic_exec

CONTAINS

  SUBROUTINE pic_exec(iout)

    USE picfield,ONLY: poissn,fftpic
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER,INTENT(OUT):: iout

!-----------------------------------------------------------------------
!----- start of main do-loop -------------------------------------------
!-----------------------------------------------------------------------
      do iloop = 1, iend

         time = iloop * dt

!----- output particle positions
!        write(41,1000) xe(100),ye(100),xe(5000),ye(5000), &
!                       xi(100),yi(100)
!1000    format(6e20.4)

         if( myid .eq. 0 ) write(*,*) iloop

         !----- charge assignment
         rho(:,:) = 0.d0
         call source(np,nx,ny,xe,ye,rho,chrge,cfact)
         call source(np,nx,ny,xi,yi,rho,chrgi,cfact)
                     !..... sum charge densities over cores
         call sumdim(nodes,myid,rho,phi,nxy)

         !----- calculate electric field
         !.......... fourier transform rho
         ifset = -1
         call fftpic(nx,ny,nxh1,nx1,ny1,rho,rhof,awk,afwk,ifset)

         !.......... calculate phi from rho in fourier space
         ipssn = 1
         call poissn(nx,ny,nxh1,rhof,phif,cform,ipssn)

         !.......... inverse fourier transform phi
         ifset = 1
         call fftpic(nx,ny,nxh1,nx1,ny1,phi,phif,awk,afwk,ifset)

         !.......... calculate ex and ey
         call efield(nx,ny,phi,ex,ey)

         !..... diagnostics to check energy conservation 
         !.....            before pushing 
         if( mod(iloop,nhmod) .eq. 0 ) then
            iene = iene + 1
            call kine(np,vxe,vye,akine1,me)
            call kine(np,vxi,vyi,akini1,mi)
            call pote(nx,ny,ex,ey,apot,cfacti)
            call sumdim1(nodes,myid,akine1,wkword)
            call sumdim1(nodes,myid,akini1,wkword)
         endif

         !----- push particles
         call push(np,nx,ny,xe,ye,vxe,vye,ex,ey,dt,ctome)
         call push(np,nx,ny,xi,yi,vxi,vyi,ex,ey,dt,ctomi)

         !----- treat particles being out of the boundary
         call bound(np,xe,ye,x1,x2,y1,y2,alx,aly)
         call bound(np,xi,yi,x1,x2,y1,y2,alx,aly)

         !..... diagnostics to check energy conservation
         !.....            after pushing 
         if( mod(iloop,nhmod) .eq. 0 ) then
            call kine(np,vxe,vye,akine2,me)
            call kine(np,vxi,vyi,akini2,mi)
            call sumdim1(nodes,myid,akine2,wkword)
            call sumdim1(nodes,myid,akini2,wkword)
            akine = 0.5d0 * ( akine1 + akine2 )
            akini = 0.5d0 * ( akini1 + akini2 )
            aktot = akine + akini
            atot  = aktot + apot

            if( iene .eq. 1 ) then
               akine0 = akine
               akini0 = akini
               aktot0 = aktot
               apot0  = apot
               atot0  = atot
            endif

            akine = akine - akine0
            akini = akini - akini0
            aktot = aktot - aktot0
            apot  = apot  - apot0
            atot  = atot  - atot0

            if( myid .eq. 0 ) &
               write(21,*) time, akine, akini, aktot, apot, atot
         endif

      end do
!-----------------------------------------------------------------------
!----- end of main do-loop ---------------------------------------------
!-----------------------------------------------------------------------

      !..... output wall clock time
      call mpi_barrier(mpi_comm_world,ierr)
      wtime2 = mpi_wtime()
      wtime  = wtime2 - wtime1

      if( myid .eq. 0 ) write(*,*) '*** wall clock time = ***', wtime
      iout=1

    END SUBROUTINE pic_exec


!***********************************************************************
      subroutine push(np,nx,ny,x,y,vx,vy,ex,ey,dt,ctom)
!***********************************************************************
      implicit none
      real(8), dimension(np) :: x, y, vx, vy
      real(8), dimension(0:nx,0:ny) :: ex, ey
      real(8) :: ctom, dx, dy, dx1, dy1, dt, exx, eyy
      integer :: np, nx, ny, i, ip, jp

      do i = 1, np

! calculate the electric field at the particle position

      ip = x(i)
      jp = y(i)

      dx = x(i) - dble(ip)
      dy = y(i) - dble(jp)
      dx1 = 1.0d0 - dx
      dy1 = 1.0d0 - dy

! electric field
 
      exx =    ex(ip ,jp  )*dx1*dy1 + ex(ip+1,jp  )*dx*dy1   &
             + ex(ip ,jp+1)*dx1*dy  + ex(ip+1,jp+1)*dx*dy
      eyy =    ey(ip ,jp  )*dx1*dy1 + ey(ip+1,jp  )*dx*dy1   &
             + ey(ip ,jp+1)*dx1*dy  + ey(ip+1,jp+1)*dx*dy

! push particles by dt

      vx(i) = vx(i) + ctom * exx * dt
      vy(i) = vy(i) + ctom * eyy * dt
      x(i)  = x(i)  + vx(i) * dt
      y(i)  = y(i)  + vy(i) * dt

      end do

    end subroutine push

!***********************************************************************
    subroutine bound(np,x,y,x1,x2,y1,y2,alx,aly)
!***********************************************************************
      implicit none
      real(8), dimension(np) :: x, y
      real(8) :: alx, aly, x1, x2, y1, y2
      integer :: np, i

      do i = 1, np
         if( x(i) .lt. x1 ) then
            x(i) = x(i) + alx
         elseif( x(i) .gt. x2 ) then
            x(i) = x(i) - alx
         endif
 
         if( y(i) .lt. y1 ) then
             y(i) = y(i) + aly
         elseif( y(i) .gt. y2 ) then
             y(i) = y(i) - aly
         endif
      end do

    end subroutine bound

!***********************************************************************
    subroutine source(np,nx,ny,x,y,rho,chrg,cfact)
!***********************************************************************
      implicit none
      real(8), dimension(np)        :: x,y
      real(8), dimension(0:nx,0:ny) :: rho
      real(8) :: chrg, dx, dy, dx1, dy1, cfact
      integer :: np, nx, ny, i, ip, jp, ix, iy

!*poption parallel, psum(rho)

      do i = 1, np

         ip = x(i)
         jp = y(i)

         dx  = x(i) - dble(ip)
         dy  = y(i) - dble(jp)
         dx1 = 1.0d0 - dx
         dy1 = 1.0d0 - dy

         rho(ip  ,jp  ) = rho(ip  ,jp  ) + dx1 * dy1 * chrg
         rho(ip+1,jp  ) = rho(ip+1,jp  ) + dx  * dy1 * chrg
         rho(ip  ,jp+1) = rho(ip  ,jp+1) + dx1 * dy  * chrg
         rho(ip+1,jp+1) = rho(ip+1,jp+1) + dx  * dy  * chrg

      end do

      !..... set charge densities at the boundary
      if( chrg .gt. 0.d0 ) then
         do iy = 0, ny
            rho(0,iy) = rho(0,iy) + rho(nx,iy)
         end do

         do ix = 0, nx
            rho(ix,0) = rho(ix,0) + rho(ix,ny)
         end do

         do iy = 0, ny
         do ix = 0, nx
            rho(ix,iy) = cfact * rho(ix,iy)
         end do
         end do
      endif

    end subroutine source

!***********************************************************************
    subroutine efield(nx,ny,phi,ex,ey)
!***********************************************************************
      implicit none
      real(8), dimension(0:nx,0:ny) :: phi, ex, ey
      integer :: nx, ny, i, j, im, ip, jm, jp

      do j = 0, ny
      do i = 0, nx

         im = i - 1
         ip = i + 1
         jm = j - 1
         jp = j + 1

         if( i .eq. 0  ) im = nx - 1
         if( i .eq. nx ) ip = 1
         if( j .eq. 0  ) jm = ny - 1
         if( j .eq. ny ) jp = 1

         ex(i,j) = 0.5d0 * ( phi(im,j ) - phi(ip,j ) )
         ey(i,j) = 0.5d0 * ( phi(i ,jm) - phi(i ,jp) )

      end do
      end do

    end subroutine efield

!***********************************************************************
    subroutine d2phi(nx,ny,phi,rho)
!***********************************************************************
!..... this subroutine is used only for check
!..... not used now
      implicit none
      real(8), dimension(0:nx,0:ny) :: phi, rho
      integer :: nx, ny, i, j, im, ip, jm, jp

      do j = 0, ny
      do i = 0, nx

         im = i - 1
         ip = i + 1
         jm = j - 1
         jp = j + 1

         if( i .eq. 0  ) im = nx - 1
         if( i .eq. nx ) ip = 1
         if( j .eq. 0  ) jm = ny - 1
         if( j .eq. ny ) jp = 1

         rho(i,j) = phi(im,j) + phi(ip,j) + phi(i,jm) + phi(i,jp)  &
                      - 4.0d0 * phi(i,j)

      end do
      end do

    end subroutine d2phi

!***********************************************************************
    subroutine kine(np,vx,vy,akin,mass)
!***********************************************************************
      implicit none
      real(8), dimension(np) :: vx, vy
      real(8) :: akin, mass 
      integer(4) :: np, i

      akin = 0.d0
      do i = 1, np
         akin = akin + vx(i)*vx(i) + vy(i)*vy(i)
      end do

      akin = 0.5 * akin * mass

    end subroutine kine

!***********************************************************************
    subroutine pote(nx,ny,ex,ey,apot,cfacti)
!***********************************************************************
      implicit none
      real(8), dimension(0:nx,0:ny) :: ex, ey 
      real(8) :: apot, cfacti 
      integer(4) :: nx, ny, ix, iy 

      apot = 0.d0
      do iy = 0, ny-1
      do ix = 0, nx-1
         apot = apot + ex(ix,iy)*ex(ix,iy) + ey(ix,iy)*ey(ix,iy)
      end do
      end do

      apot = 0.5 * cfacti * apot

    end subroutine pote

!***********************************************************************
    subroutine sumdim(nodes,myid,a,b,ndim)
!***********************************************************************
      implicit real*8(a-h,o-z)
      include'mpif.h'
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
      include'mpif.h'
      integer stat1(mpi_status_size)
      integer stat2(mpi_status_size)
      dimension a(1), b(1)
!***********************************************************************
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
 
!.....................................................
      endif
!.....................................................
 
    end subroutine sumdim1
END Module picexec
