!  ***** TASK/PIC PREPARATION *****

Module picprep
  PRIVATE
  PUBLIC pic_prep
 
CONTAINS

  SUBROUTINE pic_prep(iout)

    USE piccomm
    USE picfield,ONLY: poissn,fftpic
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER,INTENT(OUT):: iout

      np = npx * npy
      nxh1 = nx / 2 + 1
      nx1 = nx + 1
      ny1 = ny + 1
      nxy = nx1 * ny1  

      ctome  = chrge / me       !: charge to mass ratio of electrons
      ctomi  = chrgi / mi       !: charge to mass ratio of ions
      vte    = sqrt( te / me )  !: electron thermal velocity
      vti    = sqrt( ti / mi )  !: ion thermal velocity

      iene   = 0                !: counter for energy outputs
      time   = 0.D0             !: time
      iran   = 14142 * myid     !: initial parameter for random number

      !..... factors for chrage density and field energy 
      cfact  = dble(nx) * dble(ny) / dble(np) / dble(nodes) 
      cfacti = 1.d0 / cfact

      !..... constants to define boundary condition
      alx = dble(nx)
      aly = dble(ny)
      x1  = eps 
      x2  = alx - eps
      y1  = eps 
      y2  = aly - eps
      alx = alx - 2.d0 * eps
      aly = aly - 2.d0 * eps

      CALL pic_allocate

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

      iout=ierr
  END SUBROUTINE pic_prep

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

    end subroutine iniset

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

      end subroutine gauss
END Module picprep
