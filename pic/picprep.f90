!  ***** TASK/PIC PREPARATION *****

Module picprep
  PRIVATE
  PUBLIC pic_prep
 
CONTAINS

  SUBROUTINE pic_prep(iout)

    USE piccomm
    USE picfield,ONLY: poissn,fftpic
    USE piclib
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER,INTENT(OUT):: iout
    INTEGER:: nx,ny
    REAL(8):: factor

      npmax = npxmax * npymax
      nxmaxh1 = nxmax / 2 + 1
      nxmax1 = nxmax + 1
      nymax1 = nymax + 1
      nxymax = nxmax1 * nymax1

      ctome  = chrge / me       !: charge to mass ratio of electrons
      ctomi  = chrgi / mi       !: charge to mass ratio of ions
      vte    = sqrt( te / me )  !: electron thermal velocity
      vti    = sqrt( ti / mi )  !: ion thermal velocity

      time   = 0.D0             !: time
      ntgcount= 0               !: counter for global outputs
      ntpcount= 0               !: counter for profile outputs
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
      call iniset(npmax,npxmax,npymax,nxmax,nymax, &
                  xe,ye,ze,xeb,yeb,zeb,vxe,vye,vze,vte,dt,iran)

      !..... set initial positions and velocities of ions 
      call iniset(npmax,npxmax,npymax,nxmax,nymax, &
                  xi,yi,zi,xib,yib,zib,vxi,vyi,vzi,vti,dt,iran)

      !..... initialize scalar potential by poisson solver
      ipssn = 0
      call poissn(nxmax,nymax,nxmaxh1,rhof,phif,cform,ipssn)

      !..... initialize FFT 
      ifset = 0
      call fftpic(nxmax,nymax,nxmaxh1,nxmax1,nymax1,rho,rhof,awk,afwk,ifset)

      !..... initialize wall clock time
      call mpi_barrier(mpi_comm_world,ierr)
      wtime1 = mpi_wtime()

      Ax(:,:) = 0.0d0
      Ay(:,:) = 0.0d0
      Az(:,:) = 0.0d0
      Axb(:,:) = 0.0d0
      Ayb(:,:) = 0.0d0
      Azb(:,:) = 0.0d0
      phi(:,:) = 0.0d0

      DO nx=0,nxmax
         DO ny=0,nymax
            factor=DBLE(nx)/DBLE(nxmax)
            bxbg(nx,ny)=bxmin+(bxmax-bxmin)*factor
            bybg(nx,ny)=bxmin+(bxmax-bxmin)*factor
            bybg(nx,ny)=bxmin+(bxmax-bxmin)*factor
         END DO
      END DO

      call kine(npmax,vxe,vye,vze,akine0,me)
      call kine(npmax,vxi,vyi,vzi,akini0,mi)
      do ny=1,nymax
      do nx=1,nxmax  
         ex(nx,ny)=0.D0
         ey(nx,ny)=0.D0
         ez(nx,ny)=0.D0
         bx(nx,ny)=bxbg(nx,ny)
         by(nx,ny)=bybg(nx,ny)
         bz(nx,ny)=bzbg(nx,ny)
      end do
      end do
      call pote(nxmax,nymax,ex,ey,ez,bx,by,bz,vcfact,apot0)
      call sumdim1(nodes,myid,akine0,wkword)
      call sumdim1(nodes,myid,akini0,wkword)
      call sumdim1(nodes,myid,apot0,wkword)
      aktot0 = akine0 + akini0
      atot0  = aktot0 + apot0

      IF( myid .eq. 0 ) THEN
         WRITE(6,'(I8,1PE12.4,I8,1P3E12.4)') &
              nt,time,ntgcount,aktot0,apot0,atot0
      END IF

      iout=ierr
  END SUBROUTINE pic_prep

!***********************************************************************
      subroutine iniset(npmax,npxmax,npymax,nxmax,nymax, &
                        x,y,z,xb,yb,zb,vx,vy,vz,vt,dt,iran)
!***********************************************************************
      implicit none
      real(8), dimension(npmax) :: x, y, z, xb, yb, zb, vx, vy, vz
      integer :: npmax, npxmax, npymax, nxmax, nymax, iran
      real(8) :: vt, dt, factx, facty, rvx, rvy, rvz
      integer :: npx, npy, np

      factx = dble(nxmax) / dble(npxmax)
      facty = dble(nymax) / dble(npymax)

      np = 0
      do npy = 1, npymax
      do npx = 1, npxmax      
         np  = np + 1
         x(np) = ( dble(npx) - 0.5d0 ) * factx
         y(np) = ( dble(npy) - 0.5d0 ) * facty
         z(np) = 0.D0

         call gauss(rvx,rvy,rvz,iran)
         vx(np) = rvx * vt 
         vy(np) = rvy * vt
         vz(np) = rvz * vt

         xb(np) = x(np) - vx(np) * dt
         yb(np) = y(np) - vy(np) * dt
         zb(np) = z(np) - vz(np) * dt

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
