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

      npmax = npxmax * npymax
      nxmaxh1 = nxmax / 2 + 1
      nxmax1 = nxmax + 1
      nymax1 = nymax + 1
      nxymax = nxmax1 * nymax1

      ctome  = chrge / me       !: charge to mass ratio of electrons
      ctomi  = chrgi / mi       !: charge to mass ratio of ions
      vte    = sqrt( te / me )  !: electron thermal velocity
      vti    = sqrt( ti / mi )  !: ion thermal velocity

      iene   = 0                !: counter for energy outputs
      time   = 0.D0             !: time
      iran   = 14142 * myid     !: initial parameter for random number

      !..... constants to define boundary condition
      alx = dble(nxmax)
      aly = dble(nymax)
      alz = 1.D0
      x1  = eps 
      x2  = alx - eps
      y1  = eps 
      y2  = aly - eps
      z1  = eps 
      z2  = alz - eps
      alx = alx - 2.d0 * eps
      aly = aly - 2.d0 * eps
      alz = alz - 2.d0 * eps

      CALL pic_allocate

      !..... set initial positions and velocities of electrons 
      call iniset(npmax,npxmax,npymax,nxmax,nymax,xe,ye,ze,xeb,yeb,zeb,&
           vxe,vye,vze,vte,iran)
      call iniset(npmax,npxmax,npymax,nxmax,nymax,xi,yi,zi,xib,yib,zib,&
           vxi,vyi,vzi,vti,iran)

      !..... initialize poisson solver
      ipssn = 0
      call poissn(nxmax,nymax,nxmaxh1,rhof,phif,cform,ipssn)

      !..... initialize FFT 
      ifset = 0
      call fftpic(nxmax,nymax,nxmaxh1,nxmax1,nymax1,rho,rhof,awk,afwk,ifset)

      !..... initialize wall clock time
      call mpi_barrier(mpi_comm_world,ierr)
      wtime1 = mpi_wtime()

    Axb(:,:) = 0.0d0
    Ayb(:,:) = 0.0d0
    Azb(:,:) = 0.0d0
    Axbb(:,:) = 0.0d0
    Aybb(:,:) = 0.0d0
    Azbb(:,:) = 0.0d0
    phib(:,:) = 0.0d0

      iout=ierr
  END SUBROUTINE pic_prep

!***********************************************************************
      subroutine iniset(npmax,npxmax,npymax,nxmax,nymax,x,y,z,xb,yb,zb,vx,vy,vz,vt,iran)
!***********************************************************************
      implicit none
      real(8), dimension(npmax) :: x, y, z, xb, yb, zb, vx, vy, vz
      real(8) :: vt, alx, aly, factx, facty, rvx, rvy, rvz
      integer :: npmax, npxmax, npymax, nxmax, nymax, ix, iy, iz, i, iran

      alx = dble(nxmax)
      aly = dble(nymax)
      
      factx = alx / dble(npxmax)
      facty = aly / dble(npymax)

      i = 0
      do iy = 1, npymax
      do ix = 1, npxmax      
         i  = i + 1
         x(i) = ( dble(ix) - 0.5d0 ) * factx
         y(i) = ( dble(iy) - 0.5d0 ) * facty
         z(i) = 0.D0

         xb(i) = x(i)
         yb(i) = y(i)
         zb(i) = z(i)
         call gauss(rvx,rvy,rvz,iran)
         vx(i) = rvx * vt 
         vy(i) = rvy * vt
         vz(i) = rvz * vt
      end do
      end do
    end subroutine iniset

!***********************************************************************
      subroutine gauss(rvx,rvy,rvz,iran)
!***********************************************************************
      implicit none
      real(8) :: rvx, rvy, rvz, r1, r2, r3, rv
      real(8) :: pi, twopi, eps, aln ,ab
      real(8) :: rmod = 2147483648.d0, ramda = 65539.d0, wran
      integer :: iran

      pi    = 3.14159265358979d0
      twopi = 2.d0 * pi
      eps   = 0.00247875d0
      aln   = 1.d0 - eps
      call random_number(ab)

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

       !----- generate third random number
      if( iran .lt. 0 ) iran = -iran
      if( iran .eq. 0 ) iran = 3907
      wran     = iran
      wran     = mod( ramda * wran, rmod )
      iran  = wran
      r3    = wran / rmod

      !----- generate two gaussian random number
      r1  = eps + aln * r1
      rv  = sqrt( -2.d0 * log(r1) )
      rvx = rv * cos( twopi * r2 ) * sin( twopi * r3 )
      rvy = rv * sin( twopi * r2 ) * sin( twopi * r3 )
      rvz = rv * cos( twopi * r3 )

    end subroutine gauss
    
END Module picprep
