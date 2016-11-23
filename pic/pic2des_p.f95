!***********************************************************************
      program pic2des_p
!***********************************************************************
!.......................................................................
!.............****** 2 dimensional electrostatic PIC code *******.......
!.............           programmed by Hiroshi Naitou,           .......
!.............                Yamaguchi University    Dec. 2009  .......
!.......................................................................
      implicit none
      include 'mpif.h'
!.......................................................................
!............. npx   : number of particles in x                  .......
!............. npy   : number of particles in y                  .......
!............. np = npx * npy : total number of each particles   .......
!............. nx    : number of grids in x                      .......
!............. ny    : number of grids in y                      .......
!............. iend  : total time steps                          .......
!............. nhmod : output energies in each nhmod steps       .......
!.......................................................................
      integer, parameter :: npx = 1000, npy = 1000,            &
               nx = 128, ny = 128, iend = 1000,  nhmod = 1,    &
               np = npx * npy,                                 &
               nxh1 = nx / 2 + 1, nx1 = nx + 1, ny1 = ny + 1,  &
               nxy = nx1 * ny1  
!.......................................................................

      real(8), dimension(0:nx,0:ny)   :: ex, ey, rho, phi
      real(8), dimension(nx,ny)       :: awk
      real(8), dimension(np)          :: xe, ye, vxe, vye,     &
                                         xi, yi, vxi, vyi
      real(8),     dimension(nxh1,ny) :: cform 
      complex(8),  dimension(nxh1,ny) :: rhof, phif, afwk

      real(8) :: dt, me, mi, chrge, chrgi, ctome, ctomi,       &
                 te, ti, vte, vti, cfact, cfacti,              &
                 akine , akini , aktot , apot , atot ,         &
                 akine0, akini0, aktot0, apot0, atot0,         &
                 akine1, akine2, akini1, akini2, time,         &
                 eps, x1, x2, y1, y2, alx, aly,                &
                 wkword, wtime, wtime1, wtime2
      integer :: iloop, ifset, ipssn, iran, iene
      integer :: ierr, myid, nodes

!----- start parallel calculation
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,nodes,ierr)
      if( myid .eq. 0 ) write(*,*) '*** number of nodes = ***', nodes

!----- set some parameters -------------------------------
      dt     =    0.2d0     !: time step size
      me     =    1.0d0     !: electron mass
      mi     = 1836.0d0     !: ion mass
      chrge  =   -1.0d0     !: electron charge
      chrgi  =    1.0d0     !: ion charge
      te     =    1.0d0     !: electron temperature
      ti     =    1.0d0     !: ion temperature
!---------------------------------------------------------

      ctome  = chrge / me       !: charge to mass ratio of electrons
      ctomi  = chrgi / mi       !: charge to mass ratio of ions
      vte    = sqrt( te / me )  !: electron thermal velocity
      vti    = sqrt( ti / mi )  !: ion thermal velocity

      iene   = 0                !: conuter for energy outputs
      iran   = 14142 * myid     !: initial parameter for random number

      !..... factors for chrage density and field energy 
      cfact  = dble(nx) * dble(ny) / dble(np) / dble(nodes) 
      cfacti = 1.d0 / cfact

      !..... constants to define boundary condition
      eps = 0.00000001d0
      alx = dble(nx)
      aly = dble(ny)
      x1  = eps 
      x2  = alx - eps
      y1  = eps 
      y2  = aly - eps
      alx = alx - 2.d0 * eps
      aly = aly - 2.d0 * eps

      !..... set initial positions and velocities of electrons 
      call iniset(np,npx,npy,nx,ny,xe,ye,vxe,vye,vte,iran)
      call iniset(np,npx,npy,nx,ny,xi,yi,vxi,vyi,vti,iran)

      !..... initialize poisson solver
      ipssn = 0
      call poissn(nx,ny,nxh1,rhof,phif,cform,ipssn)

      !..... initialize FFT 
      ifset = 0
      call fftpic(nx,ny,nxh1,nx1,ny1,rho,rhof,awk,afwk,ifset)

      !..... initialize wall clock time
      call mpi_barrier(mpi_comm_world,ierr)
      wtime1 = mpi_wtime()

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
            call sumdim(nodes,myid,akine1,wkword,1)
            call sumdim(nodes,myid,akini1,wkword,1)
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
            call sumdim(nodes,myid,akine2,wkword,1)
            call sumdim(nodes,myid,akini2,wkword,1)
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

      !----- end parallel calcultion
      call mpi_finalize(ierr)

      end program

!***********************************************************************
      subroutine iniset(np,npx,npy,nx,ny,x,y,vx,vy,vt,iran)
!***********************************************************************
      implicit none
      real(8), dimension(np) :: x, y, vx, vy
      real(8) :: vt, alx, aly, factx, facty, rvx,rvy
      integer :: np, npx, npy, nx, ny, ix, iy, i, iran

      alx = dble(nx)
      aly = dble(ny)

      factx = alx / dble(npx)
      facty = aly / dble(npy)

      i = 0
      do iy = 1, npy
      do ix = 1, npx
         i  = i + 1
         x(i) = ( dble(ix) - 0.5d0 ) * factx
         y(i) = ( dble(iy) - 0.5d0 ) * facty
         call gauss(rvx,rvy,iran)
         vx(i) = rvx * vt 
         vy(i) = rvy * vt 
      end do
      end do

      end subroutine

!***********************************************************************
      subroutine gauss(rvx,rvy,iran)
!***********************************************************************
      implicit none
      real(8) :: rvx, rvy, r1, r2, rv
      real(8) :: pi, twopi, eps, aln 
      real(8) :: rmod = 2147483648.d0, ramda = 65539.d0, wran
      integer :: iran

      pi    = 3.14159265358979d0
      twopi = 2.d0 * pi
      eps   = 0.00247875d0
      aln   = 1.d0 - eps

      !----- generate first random number
      if( iran .lt. 0 ) iran = -iran
      if( iran .eq. 0 ) iran = 3907
      wran     = iran
      wran     = mod( ramda * wran, rmod )
      iran  = wran
      r1    = wran / rmod

      !----- generate second random number
      if( iran .lt. 0 ) iran = -iran
      if( iran .eq. 0 ) iran = 3907
      wran     = iran
      wran     = mod( ramda * wran, rmod )
      iran  = wran
      r2    = wran / rmod

      !----- generate two gaussian random number
      r1  = eps + aln * r1
      rv  = sqrt( -2.d0 * log(r1) )

      rvx = rv * cos( twopi * r2 )
      rvy = rv * sin( twopi * r2 )

      end subroutine

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

      end subroutine

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

      end subroutine

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

      end subroutine

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

      end subroutine

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

      end subroutine

!***********************************************************************
      subroutine poissn(nx,ny,nxh1,rhof,phif,cform,ipssn)
!***********************************************************************
      implicit none
      complex(8), dimension(nxh1,ny) :: rhof, phif
      real(8), dimension(nxh1,ny) :: cform
      integer(4) :: nx, ny, nxh1, nyh, i, j, ipssn
      real(8) :: alxi, alyi, pi, twopi, am, an, amn2, afsp, afsp2

         afsp  = 1.5 !+++ size of finite-size-particle
         afsp2 = afsp * afsp 

      if( ipssn .eq. 0 ) then
         pi    = 3.14159265358979d0
         twopi = 2.d0 * pi 
         alxi   = 1.d0 / dble(nx)
         alyi   = 1.d0 / dble(ny)
         nyh    = ny / 2

         do j = 1, nyh+1
         do i = 1, nxh1
            am = twopi * dble(i-1) * alxi
            an = twopi * dble(j-1) * alyi
            amn2 = am*am + an*an
           
            if( i .eq. 1 .and. j .eq. 1 ) then
               cform(i,j) = 0.0
            else 
               cform(i,j) = 1.d0 / amn2 * exp( - amn2 * afsp2 ) 
            endif

            if( j .ge. 2 .and. j .le. nyh ) then
               cform(i,ny-j+2) = cform(i,j)
            endif
         end do
         end do

      else

         !----- solve poisson equation
         do j = 1, ny
         do i = 1, nxh1
            phif(i,j) = cform(i,j) * rhof(i,j)
         end do
         end do 

      endif
      
      end subroutine

!***********************************************************************
      subroutine fftpic(nx,ny,nxh1,nx1,ny1,a,af,awk,afwk,ifset)
!***********************************************************************
      implicit none
      include 'fftw3.f'
      real(8), dimension(nx1,ny1) :: a
      real(8), dimension(nx,ny) :: awk
      real(8) :: alx, aly 
      complex(8), dimension(nxh1,ny) :: af, afwk
      integer(4) :: nx, ny, nxh1, nx1, ny1, ifset, i, j
      !....integer(8), save :: FFTW_ESTIMATE
      integer(8), save :: plan1, plan2

      alx = dble(nx)
      aly = dble(ny)

      if( ifset .eq. 0 ) then
         !----- initialization of fourier transform
         call dfftw_plan_dft_r2c_2d(plan1,nx,ny,awk,afwk,FFTW_ESTIMATE)
         call dfftw_plan_dft_c2r_2d(plan2,nx,ny,afwk,awk,FFTW_ESTIMATE)

      elseif( ifset .eq. -1 ) then

         !----- fourier transform
         do j = 1, ny
         do i = 1, nx
            awk(i,j) = a(i,j)
         end do
         end do

         call dfftw_execute(plan1)

         do j = 1, ny
         do i = 1, nxh1
         af(i,j) = afwk(i,j) / ( alx * aly )
         end do
         end do

      else

         !----- inverse fourier transform

         do j = 1, ny
         do i = 1, nxh1
         afwk(i,j) = af(i,j)
         end do
         end do

         call dfftw_execute(plan2)

         do j = 1, ny
         do i = 1, nx
            a(i,j) = awk(i,j)
         end do
         end do

         do j = 1, ny
            a(nx1,j) = a(1,j)
         end do

         do i = 1, nx1
            a(i,ny1) = a(i,1)
         end do

      endif

      end subroutine 

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

      end subroutine

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

      end subroutine

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
 
      end subroutine
