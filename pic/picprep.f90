!  ***** TASK/PIC PREPARATION *****

Module picprep
  PRIVATE
  PUBLIC pic_prep

CONTAINS

  SUBROUTINE pic_prep(iout)

    USE piccomm
    USE picsub,ONLY: poisson_f,poisson_m,efield,bfield,kine,pote
    USE piclib
    USE libmpi
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER,INTENT(OUT):: iout
    INTEGER:: nx,ny,np,locv
    REAL(8):: factor,sum

      npmax = npxmax * npymax
      nxmaxh1 = nxmax / 2 + 1
      nxmax1 = nxmax + 1
      nymax1 = nymax + 1
      nxymax = nxmax1 * nymax1
      nzmax  = MIN(nxmax,nymax)

      ctome  = chrge / me       !: charge to mass ratio of electrons
      ctomi  = chrgi / mi       !: charge to mass ratio of ions
      vte    = sqrt( te / me )  !: electron thermal velocity
      vti    = sqrt( ti / mi )  !: ion thermal velocity

      time   = 0.D0             !: time
      ntcount = 0               !: time counter
      ntgcount= 0               !: counter for global outputs
      ntpcount= 0               !: counter for profile outputs
      ntocount= 0               !: counter for orbit outputs
      iran   = 14142 * nrank    !: initial parameter for random number
      !..... constants to define boundary condition
      alx = dble(nxmax)
      aly = dble(nymax)
      alz = dble(nzmax)
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
      call iniset(npmax,npxmax,npymax,nxmax,nymax,densx, &
                  xe,ye,ze,xeb,yeb,zeb,vxe,vye,vze,vte,dt,iran)

      !..... set initial positions and velocities of ions
      call iniset(npmax,npxmax,npymax,nxmax,nymax,densx, &
                  xi,yi,zi,xib,yib,zib,vxi,vyi,vzi,vti,dt,iran)

      !..... initialize scalar potential by poisson solver
      ipssn = 0
       IF(model_boundary.EQ.0) THEN
          call poisson_f(nxmax,nymax,nxmaxh1,nxmax1,nymax1, &
                         rho,phi,rhof,phif,awk,afwk,cform,ipssn)
       END IF

      !..... initialize wall clock time
      call mpi_barrier(mpi_comm_world,ierr)
      wtime1 = mpi_wtime()

      DO nx=0,nxmax
         DO ny=0,nymax
            factor=DBLE(nx)/DBLE(nxmax)
            bxbg(nx,ny)=bxmin+(bxmax-bxmin)*factor
            bybg(nx,ny)=bymin+(bymax-bymin)*factor
            bzbg(nx,ny)=bzmin+(bzmax-bzmin)*factor
         END DO
      END DO

       !.......... calculate ex and ey and ez
       call efield(nxmax,nymax,dt,phi,Ax,Ay,Az,Axb,Ayb,Azb, &
                               ex,ey,ez,esx,esy,esz,emx,emy,emz, &
                               model_push,model_boundary)
       !.......... calculate bx and by and bz
       call bfield(nxmax,nymax,Ax,Ay,Az,Axb,Ayb,Azb, &
                               bx,by,bz,bxbg,bybg,bzbg,bb, &
                               model_push,model_boundary)
      do np=1,npmax
         vparae(np)=vye(np)
         vperpe(np)=SQRT(vxe(np)**2+vze(np)**2)
         vparai(np)=vyi(np)
         vperpi(np)=SQRT(vxi(np)**2+vzi(np)**2)
      end do

      call kine(npmax,vxe,vye,vze,akine0,me)
      call kine(npmax,vxi,vyi,vzi,akini0,mi)
      call pote(nxmax,nymax,ex,ey,ez,bx,by,bz,bxbg,bybg,bzbg,vcfact, &
                apote0,apotm0)
      call mtx_allreduce1_real8(akine0,3,sum,locv) ! sum
      akine0=sum/dble(nsize)
      call mtx_allreduce1_real8(akini0,3,sum,locv) ! sum
      akini0=sum/dble(nsize)
      call mtx_allreduce1_real8(apote0,3,sum,locv) ! sum
      apote0=sum/dble(nsize)
      call mtx_allreduce1_real8(apotm0,3,sum,locv) ! sum
      apotm0=sum/dble(nsize)

      aktot0 = akine0 + akini0
      aptot0 = apote0 + apotm0
      atot0  = aktot0 + aptot0

      IF( nrank .eq. 0 ) THEN
         WRITE(6,'(A)') &
              '      nt        time     ntg    ktot        ptot        Etot'
         WRITE(6,'(I8,1PE12.4,I8,1P3E12.4)') &
              ntcount,time,ntgcount,aktot0,aptot0,atot0
      END IF

      iout=ierr
  END SUBROUTINE pic_prep

!***********************************************************************
      subroutine iniset(npmax,npxmax,npymax,nxmax,nymax,densx,&
                        x,y,z,xb,yb,zb,vx,vy,vz,vt,dt,iran)
!***********************************************************************
      implicit none
      real(8), dimension(npmax) :: x, y, z, xb, yb, zb, vx, vy, vz
      integer :: npmax, npxmax, npymax, nxmax, nymax, iran
      real(8) :: vt, dt, factx, facty, rvx, rvy, rvz, densx, inter, position
      integer :: npx, npy, np

      factx = dble(nxmax) / dble(npxmax)
      facty = dble(nymax) / dble(npymax)
      np = 0
      if(densx .lt. 0.d0) then ! subroutine for uniform density
      do npy = 1, npymax
      do npx = 1, npxmax
        np = np + 1
         x(np) = (dble(npx) - 0.5d0 ) * factx
         y(np) = (dble(npy) - 0.5d0 ) * facty
         call gauss(rvx,rvy,rvz,iran)
         vx(np) = rvx * vt
         vy(np) = rvy * vt
         vz(np) = rvz * vt
         xb(np) = x(np) - vx(np) * dt
         yb(np) = y(np) - vy(np) * dt
         zb(np) = z(np) - vz(np) * dt
      end do
      end do
   else ! subroutine for density gradient
      inter = dble(nxmax) / (dble(npxmax) + 1.0d0 &
                          - densx * (dble(npxmax)+1.0d0)/2.0d0)
      do npy = 1, npymax
        position = 0.d0
      do npx = 1, npxmax
         np = np + 1
         position = position &
                  + inter * (1.0d0 - densx * (dble(npx) - 1.0d0)/dble(npxmax))
         x(np) = position
         y(np) = (dble(npy) - 0.5d0 ) * facty
         
         call gauss(rvx,rvy,rvz,iran)
         vx(np) = rvx * vt
         vy(np) = rvy * vt
         vz(np) = rvz * vt
         xb(np) = x(np) - vx(np) * dt
         yb(np) = y(np) - vy(np) * dt
         zb(np) = z(np) - vz(np) * dt

      end do
      end do
      end if
      end subroutine iniset

!***********************************************************************
      subroutine gauss(rvx,rvy,rvz,iran)
!***********************************************************************
      implicit none
      real(8) :: rvx, rvy, rvz, r1, r2, r3, rv
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
