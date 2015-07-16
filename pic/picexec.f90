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
    REAL(8),DIMENSION(:,:),ALLOCATABLE:: work
    INTEGER:: ienemax_old

!-----------------------------------------------------------------------
!----- allocate time hisotry data -------------------------------------------
!-----------------------------------------------------------------------

   IF(iene.eq.0) THEN
      IF(ALLOCATED(timet)) &
         DEALLOCATE(timet,akinet,akinit,aktott,apott,atott)
      ienemax=iend/nhmod
      ALLOCATE(timet(ienemax))
      ALLOCATE(akinet(ienemax))
      ALLOCATE(akinit(ienemax))
      ALLOCATE(aktott(ienemax))
      ALLOCATE(apott(ienemax))
      ALLOCATE(atott(ienemax))
   ELSE
      ALLOCATE(work(ienemax,6))
      work(1:ienemax,1)=timet (1:ienemax)
      work(1:ienemax,2)=akinet(1:ienemax)
      work(1:ienemax,3)=akinit(1:ienemax)
      work(1:ienemax,4)=aktott(1:ienemax)
      work(1:ienemax,5)=apott (1:ienemax)
      work(1:ienemax,6)=atott (1:ienemax)
      DEALLOCATE(timet,akinet,akinit,aktott,apott,atott)
      ienemax_old=ienemax
      ienemax=ienemax+iend/nhmod
      ALLOCATE(timet(ienemax))
      ALLOCATE(akinet(ienemax))
      ALLOCATE(akinit(ienemax))
      ALLOCATE(aktott(ienemax))
      ALLOCATE(apott(ienemax))
      ALLOCATE(atott(ienemax))
      timet (1:ienemax_old)=work(1:ienemax_old,1)
      akinet(1:ienemax_old)=work(1:ienemax_old,2)
      akinit(1:ienemax_old)=work(1:ienemax_old,3)
      aktott(1:ienemax_old)=work(1:ienemax_old,4)
      apott (1:ienemax_old)=work(1:ienemax_old,5)
      atott (1:ienemax_old)=work(1:ienemax_old,6)
      DEALLOCATE(work)
   END IF
      
!-----------------------------------------------------------------------
!----- start of main do-loop -------------------------------------------
!-----------------------------------------------------------------------
      do iloop = 1, iend

!         time = iloop * dt
         time = time + dt

!----- output particle positions
!        write(41,1000) xe(100),ye(100),xe(5000),ye(5000), &
!                       xi(100),yi(100)
!1000    format(6e20.4)

         if( myid .eq. 0 ) write(6,'(2I8,1PE12.4)') iloop, iene+1, time

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

         !.......... calculate ex and ey and ez
         call efield(nx,ny,nz,phi,ex,ey,ez,ezg)

         !.......... calculate bxg and byg and bzg
         call bfield(nx,ny,nz,bx,by,bz,bxg,byg,bzg)
         

         !..... diagnostics to check energy conservation 
         !.....            before pushing 
         if( mod(iloop,nhmod) .eq. 0 ) then
            iene = iene + 1
            call kine(np,vxe,vye,vze,akine1,me)
            call kine(np,vxi,vyi,vzi,akini1,mi)
            call pote(nx,ny,nz,ex,ey,ezg,apot,cfacti)
            call sumdim1(nodes,myid,akine1,wkword)
            call sumdim1(nodes,myid,akini1,wkword)
         endif

         !----- push particles
         call push(np,nx,ny,nz,xe,ye,ze,vxe,vye,vze,ex,ey,ezg,bxg,byg,bzg,dt,&
ctome)
         call push(np,nx,ny,nz,xi,yi,zi,vxi,vyi,vzi,ex,ey,ezg,bxg,byg,bzg,dt,&
ctomi)

         !----- treat particles being out of the boundary
         call bound(np,nx,ny,nz,xe,ye,ze,x1,x2,y1,y2,z1,z2,alx,aly,alz)
         call bound(np,nx,ny,nz,xi,yi,zi,x1,x2,y1,y2,z1,z2,alx,aly,alz)

         !..... diagnostics to check energy conservation
         !.....            after pushing 
         if( mod(iloop,nhmod) .eq. 0 ) then
            call kine(np,vxe,vye,vze,akine2,me)
            call kine(np,vxi,vyi,vzi,akini2,mi)
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
            
            timet(iene)=time
            akinet(iene)=akine
            akinit(iene)=akini
            aktott(iene)=aktot
            apott(iene)=apot
            atott(iene)=atot  
          
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

      ienemax=iene
      iout=1



       END SUBROUTINE pic_exec

!***********************************************************************
      subroutine push(np,nx,ny,nz,x,y,z,vx,vy,vz,ex,ey,ezg,bxg,byg,bzg,dt,ctom)
!***********************************************************************
      implicit none
      real(8), dimension(np) :: x, y, z, vx, vy, vz
      real(8), dimension(0:nx,0:ny,0:nz) :: ex, ey, ezg, bxg, byg, bzg
      real(8) :: ctom, dx, dy, dz, dx1, dy1, dz1, dt, exx, eyy, ezz, ez, bxx,&
 byy, bzz, a, b, c
      integer :: np, nx, ny, nz, i, ip, jp, kp

      do i = 1, np

! calculate the electric field at the particle position

      ip = x(i)
      jp = y(i)
      kp = z(i)

      dx = x(i) - dble(ip)
      dy = y(i) - dble(jp)
      dz = z(i) - dble(kp)
      dx1 = 1.0d0 - dx
      dy1 = 1.0d0 - dy
      dz1 = 1.0d0 - dz

! electric field
 
      exx =  ex(ip  ,jp  ,kp  )*dx1*dy1*dz1 + ex(ip+1,jp  ,kp  )*dx*dy1*dz1 &
           + ex(ip  ,jp+1,kp  )*dx1*dy*dz1  + ex(ip  ,jp  ,kp+1)*dx1*dy1*dz &
           + ex(ip+1,jp+1,kp  )*dx*dy*dz1   + ex(ip  ,jp+1,kp+1)*dx1*dy*dz &
           + ex(ip+1,jp  ,kp+1)*dx*dy1*dz   + ex(ip+1,jp+1,kp+1)*dx*dy*dz

      eyy =  ey(ip  ,jp  ,kp  )*dx1*dy1*dz1 + ey(ip+1,jp  ,kp  )*dx*dy1*dz1 &
           + ey(ip  ,jp+1,kp  )*dx1*dy*dz1  + ey(ip  ,jp  ,kp+1)*dx1*dy1*dz &
           + ey(ip+1,jp+1,kp  )*dx*dy*dz1   + ey(ip  ,jp+1,kp+1)*dx1*dy*dz &
           + ey(ip+1,jp  ,kp+1)*dx*dy1*dz   + ey(ip+1,jp+1,kp+1)*dx*dy*dz

      ezz =  ezg(ip  ,jp  ,kp  )*dx1*dy1*dz1 + ezg(ip+1,jp  ,kp  )*dx*dy1*dz1 &
           + ezg(ip  ,jp+1,kp  )*dx1*dy*dz1  + ezg(ip  ,jp  ,kp+1)*dx1*dy1*dz &
           + ezg(ip+1,jp+1,kp  )*dx*dy*dz1   + ezg(ip  ,jp+1,kp+1)*dx1*dy*dz &
           + ezg(ip+1,jp  ,kp+1)*dx*dy1*dz   + ezg(ip+1,jp+1,kp+1)*dx*dy*dz


! magnetic field
      
    bxx =    bxg(ip  ,jp  ,kp  )*dx1*dy1*dz1 + bxg(ip+1,jp  ,kp  )*dx*dy1*dz1 &
           + bxg(ip  ,jp+1,kp  )*dx1*dy*dz1  + bxg(ip  ,jp  ,kp+1)*dx1*dy1*dz &
           + bxg(ip+1,jp+1,kp  )*dx*dy*dz1   + bxg(ip  ,jp+1,kp+1)*dx1*dy*dz &
           + bxg(ip+1,jp  ,kp+1)*dx*dy1*dz   + bxg(ip+1,jp+1,kp+1)*dx*dy*dz

    byy =    byg(ip  ,jp  ,kp  )*dx1*dy1*dz1 + byg(ip+1,jp  ,kp  )*dx*dy1*dz1 &
           + byg(ip  ,jp+1,kp  )*dx1*dy*dz1  + byg(ip  ,jp  ,kp+1)*dx1*dy1*dz &
           + byg(ip+1,jp+1,kp  )*dx*dy*dz1   + byg(ip  ,jp+1,kp+1)*dx1*dy*dz &
           + byg(ip+1,jp  ,kp+1)*dx*dy1*dz   + byg(ip+1,jp+1,kp+1)*dx*dy*dz

    bzz =    bzg(ip  ,jp  ,kp  )*dx1*dy1*dz1 + bzg(ip+1,jp  ,kp  )*dx*dy1*dz1 &
           + bzg(ip  ,jp+1,kp  )*dx1*dy*dz1  + bzg(ip  ,jp  ,kp+1)*dx1*dy1*dz &
           + bzg(ip+1,jp+1,kp  )*dx*dy*dz1   + bzg(ip  ,jp+1,kp+1)*dx1*dy*dz &
           + bzg(ip+1,jp  ,kp+1)*dx*dy1*dz   + bzg(ip+1,jp+1,kp+1)*dx*dy*dz

    ! push particles by dt
    a = vx(i)
    b = vy(i)
    c = vz(i)
     

      vx(i) = vx(i) + ctom * (exx + b * bzz - c * byy) * dt
      vy(i) = vy(i) + ctom * (eyy + c * bxx - a * bzz) * dt
      vz(i) = vz(i) + ctom * (ezz + a * byy - b * bxx) * dt
      x(i)  = x(i)  + vx(i) * dt
      y(i)  = y(i)  + vy(i) * dt
      z(i)  = z(i)  + vz(i) * dt

      end do

    end subroutine push

!***********************************************************************
    subroutine bound(np,nx,ny,nz,x,y,z,x1,x2,y1,y2,z1,z2,alx,aly,alz)
!***********************************************************************
      implicit none
      real(8), dimension(np) :: x, y, z
      real(8) :: alx, aly, alz, x1, x2, y1, y2, z1, z2
      integer :: np, i, nx, ny, nz

      do i = 1, np
         if( x(i) .lt. x1 ) then
            do while(x(i) .lt. x1 + nx)
               x(i) = x(i) + dble(nx)
            end do
         x(i) = x(i)+ alx
         elseif( x(i) .gt. x2 ) then
            do while(x(i) .gt. x2 - nx)
               x(i) = x(i) - dble(nx)
            end do
         x(i) = x(i) - alx
         endif
 
         if( y(i) .lt. y1 ) then
            do while(y(i) .lt. y1 + ny)
               y(i) = y(i) + dble(ny)
            end do
         y(i) = y(i) - aly
         elseif( y(i) .gt. y2 ) then
            do while(y(i) .gt. y2 - ny)
               y(i) = y(i) - dble(ny)
            end do
         y(i) = y(i) - aly
         endif

         if( z(i) .lt. z1 ) then
            do while(z(i) .lt. z1 + nz)
               z(i) = z(i) + dble(nz)
            end do
         z(i) = z(i) + alz
         elseif( z(i) .gt. z2 ) then
            do while(z(i) .gt. z2 - nz)
               z(i) = z(i) - dble(nz)
            end do
         z(i) = z(i) - alz
         endif
          
      end do

    end subroutine bound

!***********************************************************************
    subroutine source(np,nx,ny,x,y,rho,chrg,cfact)
!***********************************************************************
      implicit none
      real(8), dimension(np)        :: x, y
      real(8), dimension(0:nx,0:ny) :: rho
      real(8) :: chrg, dx, dy, dx1, dy1, cfact
      integer :: np, nx, ny, i, ip, jp, kp, ix, iy

!*poption parallel, psum(rho)

      do i = 1, np

         ip = x(i)
         jp = y(i)
         !kp = z(i)

         dx  = x(i) - dble(ip)
         dy  = y(i) - dble(jp)
         !dz  = z(i) - dble(kp)
         dx1 = 1.0d0 - dx
         dy1 = 1.0d0 - dy
         !dz1 = 1.0d0 - dz

         rho(ip  ,jp  ) = rho(ip  ,jp  ) + dx1 * dy1 * chrg
         rho(ip+1,jp  ) = rho(ip+1,jp  ) + dx  * dy1 * chrg
         rho(ip  ,jp+1) = rho(ip  ,jp+1) + dx1 * dy  * chrg
         rho(ip  ,jp  ) = rho(ip  ,jp  ) + dx1 * dy1 * chrg
         !rho(ip+1,jp+1,kp+1) = rho(ip+1,jp+1,kp+1) + dx  * dy  * dz  * chrg

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
   subroutine efield(nx,ny,nz,phi,ex,ey,ez,ezg)
!***********************************************************************
      implicit none
      real(8), dimension(0:nx,0:ny) :: phi
      real(8), dimension(0:nx,0:ny,0:nz) :: ex, ey, ezg
      integer :: nx, ny, nz, i, j, im, ip, jm, jp, k, km, kp
      real(rkind)::ez

      do j = 0, ny
      do i = 0, nx
      do k = 0, nz

         im = i - 1
         ip = i + 1
         jm = j - 1
         jp = j + 1
         km = k - 1
         kp = k + 1

         if( i .eq. 0  ) im = nx - 1
         if( i .eq. nx ) ip = 1
         if( j .eq. 0  ) jm = ny - 1
         if( j .eq. ny ) jp = 1
         if( k .eq. 0  ) km = nz - 1
         if( k .eq. nz ) kp = 1

         ex(i,j,k) = 0.5d0 * ( phi(im,j ) - phi(ip,j ) )
         ey(i,j,k) = 0.5d0 * ( phi(i ,jm) - phi(i ,jp) )
         ezg(i,j,k) = ez

      end do
      end do
      end do
    end subroutine efield
!***********************************************************************
    subroutine bfield(nx,ny,nz,bx,by,bz,bxg,byg,bzg)
!***********************************************************************
      implicit none
      real(rkind), dimension(0:nx,0:ny,0:nz) :: bxg, byg, bzg
      integer :: nx, ny, nz, i, j, k, im, ip, jm, jp, km, kp
      real(rkind)::bx, by, bz

      do j = 0, ny
      do i = 0, nx
      do k = 0, nz

         im = i - 1
         ip = i + 1
         jm = j - 1
         jp = j + 1
         km = k - 1
         kp = k + 1

         if( i .eq. 0  ) im = nx - 1
         if( i .eq. nx ) ip = 1
         if( j .eq. 0  ) jm = ny - 1
         if( j .eq. ny ) jp = 1
         if( k .eq. 0  ) km = nz - 1
         if( k .eq. ny ) kp = 1
         

         bxg(i,j,k) = bx
         byg(i,j,k) = by
         bzg(i,j,k) = bz

      end do
      end do
      end do

    end subroutine bfield

    

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
    subroutine kine(np,vx,vy,vz,akin,mass)
!***********************************************************************
      implicit none
      real(8), dimension(np) :: vx, vy, vz
      real(8) :: akin, mass 
      integer(4) :: np, i
      akin = 0.d0
      do i = 1, np
         akin = akin + vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i)
      end do

      akin = 0.5 * akin * mass
    end subroutine kine

!***********************************************************************
    subroutine pote(nx,ny,nz,ex,ey,ezg,apot,cfacti)
!***********************************************************************
      implicit none
      real(8), dimension(0:nx,0:ny,0:nz) :: ex, ey, ezg
      real(8) ::  apot, cfacti 
      integer(4) :: nx, ny, nz, ix, iy, iz 

      apot = 0.d0
      do iy = 0, ny-1
      do ix = 0, nx-1
      do iz = 0, nz-1      
         apot = apot + ex(ix,iy,iz)*ex(ix,iy,iz) + ey(ix,iy,iz)*ey(ix,iy,iz)&
+ezg(ix,iy,iz)*ezg(ix,iy,iz)
      end do
      end do
      end do

      apot = 0.5 * cfacti * apot

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
