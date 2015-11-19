!./pic  ***** TASK/PIC EXEC *****

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
    INTEGER:: ienemax_old,i,j
    REAL(8)::abc
    INTEGER::ii,jj
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

!       time = iloop * dt
         time = time + dt

!----- output particle positions
!        write(41,1000) xe(100),ye(100),xe(5000),ye(5000), &
!                       xi(100),yi(100)
!1000    format(6e20.4)

         if( myid .eq. 0 ) write(6,'(2I8,1PE12.4)') iloop, iene+1, time
         
         !----- charge assignment
         rho(:,:)=0.0d0
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
         jx(:,:)=0.d0
         jy(:,:)=0.d0
         jz(:,:)=0.d0
         !.......... calculate current by electrons
         call current(np,nx,ny,xe,ye,vxe,vye,vze,jx,jy,jz,dt,chrge,cfact,&
              jxbg,jybg,jzbg,omega,thetax,thetay,thetaz,time)
         !.......... calculate current by ions
         call current(np,nx,ny,xi,yi,vxi,vyi,vzi,jx,jy,jz,dt,chrgi,cfact,&
              jxbg,jybg,jzbg,omega,thetax,thetay,thetaz,time)
         !.......... calculate vector potential
         call phia(nx,ny,c,dt,phi,phib,jx,jy,jz,Ax,Ay,Az,Axb,Ayb,Azb,&
              Axbb,Aybb,Azbb,cfact)

           !.......... calculate ex and ey and ez
         call efield(nx,ny,phi,ex,ey,ez,Axb,Ayb,Azb,Axbb,Aybb,Azbb)
           !.......... calculate bxg and byg and bzg
         call bfield(nx,ny,bxmin,bxmax,bymin,bymax,bzmin,bzmax,&
              bx,by,bz,Axb,Ayb,Azb)
         
          if( mod(iloop,nhmod) .eq. 0 ) then
            iene = iene + 1
            call kine(np,vxe,vye,vze,akine1,me,cfacti)
            call kine(np,vxi,vyi,vzi,akini1,mi,cfacti)
            call pote(nx,ny,ex,ey,ez,bx,by,bz,c,apot,cfacti)
            call sumdim1(nodes,myid,akine1,wkword)
            call sumdim1(nodes,myid,akini1,wkword)
         endif
         !..... push electrons
         call push(np,nx,ny,xe,ye,ze,vxe,vye,vze,ex,ey,ez,bx,by,bz,dt,&
              ctome,xeb,yeb,zeb,phi,phib,Axb,Ayb,Azb,Axbb,Aybb,Azbb, &
              vparae,vperpe)
         
         !..... push ions
         call push(np,nx,ny,xi,yi,zi,vxi,vyi,vzi,ex,ey,ez,bx,by,bz,dt,&
              ctomi,xib,yib,zib,phi,phib,Axb,Ayb,Azb,Axbb,Aybb,Azbb, &
              vparai,vperpi)

         !----- treat particles being out of the boundary
         call bound(np,xe,ye,ze,x1,x2,y1,y2,z1,z2,alx,aly,alz)
         call bound(np,xi,yi,zi,x1,x2,y1,y2,z1,z2,alx,aly,alz)

         !..... diagnostics to check energy conservation
         !.....            after pushing 
         if( mod(iloop,nhmod) .eq. 0 ) then
            call kine(np,vxe,vye,vze,akine2,me,cfacti)
            call kine(np,vxi,vyi,vzi,akini2,mi,cfacti)
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
       subroutine push(np,nx,ny,x,y,z,vx,vy,vz,ex,ey,ez,bx,by,bz,dt,&
            ctom,xb,yb,zb,phi,phib,Axb,Ayb,Azb,Axbb,Aybb,Azbb,vpara,vperp)
!***********************************************************************
      implicit none
      real(8), dimension(np) :: x, y, z, xb, yb, zb, vx, vy, vz, vpara, vperp
      real(8), dimension(0:nx,0:ny) :: phi,phib,Axb, Ayb, Azb, Axbb, Aybb, Azbb
      real(8), dimension(0:nx,0:ny) :: ex, ey, ez, bx, by, bz
      real(8) :: ctom, dx, dy, dx1, dy1, dt, exx, eyy, ezz, bxx,&
                 byy, bzz, vxn, vyn, vzn, vxzero, vyzero, vzzero, vxp, vyp,&
                 vzp, sx1p, sx1m, sy1p, sy1m, sx2, sy2, sx2p, sx2m, sy2m, sy2p
      real(8) :: btot, vtot
      integer :: np, nx, ny, nz, i, j, ip, jp, ipp, ipm, jpp, jpm, ippp, jppp
     
      do i = 1, np

! calculate the electric field at the particle position

         ip = x(i)
         jp = y(i)
         dx = x(i) - dble(ip)
         dy = y(i) - dble(jp)
         dx1 = 1.0d0 - dx
         dy1 = 1.0d0 - dy

          if(dx .le. 0.5d0) then
            sx2  = 3.0d0/4 - dx ** 2
            sx2p = 1.0d0/2 * (1.0d0/2 + dx) ** 2
            sx2m = 1.0d0/2 * (1.0d0/2 - dx) ** 2
         else
            sx2  = 3.0d0/4 - (dx - 1.0d0) ** 2
            sx2p = 1.0d0/2 * (-1.0d0/2 + dx) ** 2
            sx2m = 1.0d0/2 * (3.0d0/2 - dx) ** 2
         endif
         if(dy .le. 0.5d0) then
            sy2  = 3.0d0/4 - dy ** 2
            sy2p = 1.0d0/2 * (1.0d0/2 + dy) ** 2
            sy2m = 1.0d0/2 * (1.0d0/2 - dy) ** 2
         else
            sy2  = 3.0d0/4 - (dy - 1.0d0) ** 2
            sy2p = 1.0d0/2 * (-1.0d0/2 + dy) ** 2
            sy2m = 1.0d0/2 * (3.0d0/2 - dy) ** 2
         endif

         ipm = ip - 1
         ipp = ip + 1
         jpm = jp - 1
         jpp = jp + 1
         ippp = ip + 2
         jppp = jp + 2
         if( ip .eq. 0  ) ipm = nx - 1
         if( jp .eq. 0  ) jpm = ny - 1
         if( ip .eq. nx-1) ippp = 1
         if( jp .eq. ny-1) jppp = 1
         ! electric field and magnetic field
         if(dx .le. 0.5d0 .and. dy .le. 0.5d0) then
            exx = ex(ipp,jpp)*dx*sy2p + ex(ip ,jpp)*dx1*sy2p &
                + ex(ipp,jp )*dx*sy2  + ex(ip ,jp )*dx1*sy2  &
                + ex(ipp,jpm)*dx*sy2m + ex(ip ,jpm)*dx1*sy2m
        
            eyy = ey(ipp,jpp)*sx2p*dy + ey(ipp,jp )*sx2p*dy1 &
                + ey(ip ,jpp)*sx2 *dy + ey(ip ,jp )*sx2 *dy1 &
                + ey(ipm,jpp)*sx2m*dy + ey(ipm,jp )*sx2m*dy1 
            
            ezz = ez(ipp,jpp)*sx2p*sy2p + ez(ipp,jp )*sx2p*sy2 &
                + ez(ipp,jpm)*sx2p*sy2m + ez(ip ,jpp)*sx2 *sy2p&
                + ez(ip ,jp )*sx2 *sy2  + ez(ip ,jpm)*sx2 *sy2m&
                + ez(ipm,jpp)*sx2m*sy2p + ez(ipm,jp )*sx2m*sy2 &
                + ez(ipm,jpm)*sx2m*sy2m
            
            bxx = bx(ipp,jpp)*dx*sy2p + bx(ip ,jpp)*dx1*sy2p &
                + bx(ipp,jp )*dx*sy2  + bx(ip ,jp )*dx1*sy2  &
                + bx(ipp,jpm)*dx*sy2m + bx(ip ,jpm)*dx1*sy2m
        
            byy = by(ipp,jpp)*sx2p*dy + by(ipp,jp )*sx2p*dy1 &
                + by(ip ,jpp)*sx2 *dy + by(ip ,jp )*sx2 *dy1 &
                + by(ipm,jpp)*sx2m*dy + by(ipm,jp )*sx2m*dy1
        
            bzz = bz(ipp,jpp)*sx2p*sy2p + bz(ip ,jpp)*sx2 *sy2p &
                + bz(ipm,jpp)*sx2m*sy2p + bz(ipp,jp )*sx2p*sy2  &
                + bz(ip ,jp )*sx2 *sy2  + bz(ipm,jp )*sx2m*sy2  &
                + bz(ipp,jpm)*sx2p*sy2m + bz(ip ,jpm)*sx2 *sy2m &
                + bz(ipm,jpm)*sx2m*sy2m

         else if(dx .le. 0.5d0 .and. dy .ge. 0.5d0) then
            exx = ex(ipp,jppp)*dx*sy2p + ex(ip ,jppp)*dx1*sy2p &
                + ex(ipp,jpp )*dx*sy2  + ex(ip ,jpp )*dx1*sy2  &
                + ex(ipp,jp  )*dx*sy2m + ex(ip ,jp  )*dx1*sy2m
        
            eyy = ey(ipp,jpp)*sx2p*dy + ey(ipp,jp )*sx2p*dy1 &
                + ey(ip ,jpp)*sx2 *dy + ey(ip ,jp )*sx2 *dy1 &
                + ey(ipm,jpp)*sx2m*dy + ey(ipm,jp )*sx2m*dy1 
            
            ezz = ez(ipp,jppp)*sx2p*sy2p + ez(ipp,jpp )*sx2p*sy2 &
                + ez(ipp,jp  )*sx2p*sy2m + ez(ip ,jppp)*sx2 *sy2p&
                + ez(ip ,jpp )*sx2 *sy2  + ez(ip ,jp  )*sx2 *sy2m&
                + ez(ipm,jppp)*sx2m*sy2p + ez(ipm,jpp )*sx2m*sy2 &
                + ez(ipm,jp  )*sx2m*sy2m
      
            bxx = bx(ipp,jppp)*dx*sy2p + bx(ip ,jppp)*dx1*sy2p &
                + bx(ipp,jpp )*dx*sy2  + bx(ip ,jpp )*dx1*sy2  &
                + bx(ipp,jp  )*dx*sy2m + bx(ip ,jp  )*dx1*sy2m
        
            byy = by(ipp,jpp)*sx2p*dy + by(ipp,jp )*sx2p*dy1 &
                + by(ip ,jpp)*sx2 *dy + by(ip ,jp )*sx2 *dy1 &
                + by(ipm,jpp)*sx2m*dy + by(ipm,jp )*sx2m*dy1
        
            bzz = bz(ipp,jppp)*sx2p*sy2p + bz(ip ,jppp)*sx2 *sy2p &
                + bz(ipm,jppp)*sx2m*sy2p + bz(ipp,jpp )*sx2p*sy2  &
                + bz(ip ,jpp )*sx2 *sy2  + bz(ipm,jpp )*sx2m*sy2  &
                + bz(ipp,jp  )*sx2p*sy2m + bz(ip ,jp  )*sx2 *sy2m &
                + bz(ipm,jp  )*sx2m*sy2m

            else if(dx .ge. 0.5d0 .and. dy .le. 0.5d0) then
            exx = ex(ipp,jpp)*dx*sy2p + ex(ip ,jpp)*dx1*sy2p &
                + ex(ipp,jp )*dx*sy2  + ex(ip ,jp )*dx1*sy2  &
                + ex(ipp,jpm)*dx*sy2m + ex(ip ,jpm)*dx1*sy2m
        
            eyy = ey(ippp,jpp)*sx2p*dy + ey(ippp,jp )*sx2p*dy1 &
                + ey(ipp ,jpp)*sx2 *dy + ey(ipp ,jp )*sx2 *dy1 &
                + ey(ip  ,jpp)*sx2m*dy + ey(ip  ,jp )*sx2m*dy1 
            
            ezz = ez(ippp,jpp)*sx2p*sy2p + ez(ippp,jp )*sx2p*sy2 &
                + ez(ippp,jpm)*sx2p*sy2m + ez(ipp ,jpp)*sx2 *sy2p&
                + ez(ipp ,jp )*sx2 *sy2  + ez(ipp ,jpm)*sx2 *sy2m&
                + ez(ip  ,jpp)*sx2m*sy2p + ez(ip  ,jp )*sx2m*sy2 &
                + ez(ip  ,jpm)*sx2m*sy2m
            
            bxx = bx(ipp,jpp)*dx*sy2p + bx(ip ,jpp)*dx1*sy2p &
                + bx(ipp,jp )*dx*sy2  + bx(ip ,jp )*dx1*sy2  &
                + bx(ipp,jpm)*dx*sy2m + bx(ip ,jpm)*dx1*sy2m
        
            byy = by(ippp,jpp)*sx2p*dy + by(ippp,jp )*sx2p*dy1 &
                + by(ipp ,jpp)*sx2 *dy + by(ipp ,jp )*sx2 *dy1 &
                + by(ip  ,jpp)*sx2m*dy + by(ip  ,jp )*sx2m*dy1
        
            bzz = bz(ippp,jpp)*sx2p*sy2p + bz(ipp ,jpp)*sx2 *sy2p &
                + bz(ip  ,jpp)*sx2m*sy2p + bz(ippp,jp )*sx2p*sy2  &
                + bz(ipp ,jp )*sx2 *sy2  + bz(ip  ,jp )*sx2m*sy2  &
                + bz(ippp,jpm)*sx2p*sy2m + bz(ipp ,jpm)*sx2 *sy2m &
                + bz(ip  ,jpm)*sx2m*sy2m

            else
            exx = ex(ipp,jppp)*dx*sy2p + ex(ip ,jppp)*dx1*sy2p &
                + ex(ipp,jpp )*dx*sy2  + ex(ip ,jpp )*dx1*sy2  &
                + ex(ipp,jp  )*dx*sy2m + ex(ip ,jp  )*dx1*sy2m
        
            eyy = ey(ippp,jpp)*sx2p*dy + ey(ippp,jp )*sx2p*dy1 &
                + ey(ipp ,jpp)*sx2 *dy + ey(ipp ,jp )*sx2 *dy1 &
                + ey(ip  ,jpp)*sx2m*dy + ey(ip  ,jp )*sx2m*dy1 
            
            ezz = ez(ippp,jppp)*sx2p*sy2p + ez(ippp,jpp )*sx2p*sy2 &
                + ez(ippp,jp  )*sx2p*sy2m + ez(ipp ,jppp)*sx2 *sy2p&
                + ez(ipp ,jpp )*sx2 *sy2  + ez(ipp ,jp  )*sx2 *sy2m&
                + ez(ip  ,jppp)*sx2m*sy2p + ez(ip  ,jpp )*sx2m*sy2 &
                + ez(ip  ,jp  )*sx2m*sy2m
            
            bxx = bx(ipp,jppp)*dx*sy2p + bx(ip ,jppp)*dx1*sy2p &
                + bx(ipp,jpp )*dx*sy2  + bx(ip ,jpp )*dx1*sy2  &
                + bx(ipp,jp  )*dx*sy2m + bx(ip ,jp  )*dx1*sy2m
        
            byy = by(ippp,jpp)*sx2p*dy + by(ippp,jp )*sx2p*dy1 &
                + by(ipp ,jpp)*sx2 *dy + by(ipp ,jp )*sx2 *dy1 &
                + by(ip  ,jpp)*sx2m*dy + by(ip  ,jp )*sx2m*dy1
        
            bzz = bz(ippp,jppp)*sx2p*sy2p + bz(ipp ,jppp)*sx2 *sy2p &
                + bz(ip  ,jppp)*sx2m*sy2p + bz(ippp,jpp )*sx2p*sy2  &
                + bz(ipp ,jpp )*sx2 *sy2  + bz(ip  ,jpp )*sx2m*sy2  &
                + bz(ippp,jp  )*sx2p*sy2m + bz(ipp ,jp  )*sx2 *sy2m &
                + bz(ip  ,jp  )*sx2m*sy2m
         endif
        
            ! push particles by using Buneman-Boris method
         vxn = vx(i) + 1.0d0/2 * ctom * exx * dt 
         vyn = vy(i) + 1.0d0/2 * ctom * eyy * dt
         vzn = vz(i) + 1.0d0/2 * ctom * ezz * dt

         vxzero = vxn + 1.0d0/2 * ctom * (vyn * bzz - vzn * byy) * dt
         vyzero = vyn + 1.0d0/2 * ctom * (vzn * bxx - vxn * bzz) * dt 
         vzzero = vzn + 1.0d0/2 * ctom * (vxn * byy - vyn * bxx) * dt

         vxp = vxn + 1.0d0/(1.0d0 + 0.25d0 * (ctom * dt) ** 2 & 
             * (bxx ** 2 + byy ** 2 + bzz ** 2)) & 
             * ctom * (vyzero * bzz - vzzero * byy) * dt
         vyp = vyn + 1.0d0/(1.0d0 + 0.25d0 * (ctom * dt) ** 2 &
             * (bxx ** 2 + byy ** 2 + bzz ** 2)) &
             * ctom * (vzzero * bxx - vxzero * bzz) * dt 
         vzp = vzn + 1.0d0/(1.0d0 + 0.25d0 * (ctom * dt) ** 2 & 
             * (bxx ** 2 + byy ** 2 + bzz ** 2)) & 
             * ctom * (vxzero * byy - vyzero * bxx) * dt
 
         vx(i) = vxp + 1.0d0/2 * ctom * exx * dt
         vy(i) = vyp + 1.0d0/2 * ctom * eyy * dt
         vz(i) = vzp + 1.0d0/2 * ctom * ezz * dt
         
         xb(i) = x(i)
         yb(i) = y(i)
         zb(i) = z(i)
         x(i) = x(i) + vx(i) * dt
         y(i) = y(i) + vy(i) * dt
         z(i) = z(i) + vz(i) * dt

         btot=SQRT(bxx**2+byy**2+bzz**2)
         IF(btot.EQ.0.D0) THEN
            vpara(i)=vx(i)
            vperp(i)=SQRT(vy(i)**2+vz(i)**2)
         ELSE
            vtot=SQRT(vx(i)**2+vy(i)**2+vz(i)**2)
            vpara(i)=(bxx*vx(i)+byy*vy(i)+bzz*vz(i))/btot
            vperp(i)=(SQRT(vtot**2-vpara(i)**2))
         END IF

      end do

    end subroutine push

!***********************************************************************
    subroutine bound(np,x,y,z,x1,x2,y1,y2,z1,z2,alx,aly,alz)
!***********************************************************************
      implicit none
      real(8), dimension(np) :: x, y, z
      real(8) :: alx, aly, alz, x1, x2, y1, y2, z1, z2
      integer :: np, i
      do i = 1, np
         if( x(i) .lt. x1 ) then
            do while(x(i) .lt. x1)
               x(i) = x(i) + alx
            end do
         elseif( x(i) .gt. x2 ) then
            do while(x(i) .gt. x2)
               x(i) = x(i) - alx
            end do
         endif
 
         if( y(i) .lt. y1 ) then
            do while(y(i) .lt. y1)
               y(i) = y(i) + aly
            end do
         elseif( y(i) .gt. y2 ) then
            do while(y(i) .gt. y2)
               y(i) = y(i) - aly
            end do
         endif

         if( z(i) .lt. z1 ) then
            do while(z(i) .lt. z1)
               z(i) = z(i) + alz
            end do
         elseif( z(i) .gt. z2 ) then
            do while(z(i) .gt. z2)
               z(i) = z(i) - alz
            end do
         endif
      end do

    end subroutine bound

!***********************************************************************
    subroutine source(np,nx,ny,x,y,rho,chrg,cfact)
!***********************************************************************
      implicit none
      real(8), dimension(np)        :: x, y
      real(8), dimension(0:nx,0:ny) :: rho
      real(8) :: chrg, cfact, dx, dy, dz, sx2, sy2, sx2p, sy2p, sx2m, sy2m,&
           dx1,dy1
      integer :: np, nx, ny, i, ip, jp, ix, iy, ipp, ipm, jpp, jpm, ippp, jppp

!*poption parallel, psum(rho)

      do i = 1, np

         ip = x(i)
         jp = y(i)

         dx  = x(i) - dble(ip)
         dy  = y(i) - dble(jp)
         dx1 = 1.0d0 - dx
         dy1 = 1.0d0 - dy
         if(dx .le. 0.5d0) then
            sx2  = 3.0d0/4 - dx ** 2
            sx2p = 1.0d0/2 * (1.0d0/2 + dx) ** 2
            sx2m = 1.0d0/2 * (1.0d0/2 - dx) ** 2
         else
            sx2  = 3.0d0/4 - (1.0d0 - dx) ** 2
            sx2p = 1.0d0/2 * (-1.0d0/2 + dx) ** 2
            sx2m = 1.0d0/2 * (3.0d0/2 - dx) ** 2
         endif
         if(dy .le. 0.5d0) then
            sy2  = 3.0d0/4 - dy ** 2
            sy2p = 1.0d0/2 * (1.0d0/2 + dy) ** 2
            sy2m = 1.0d0/2 * (1.0d0/2 - dy) ** 2
         else
            sy2  = 3.0d0/4 - (1.0d0 - dy) ** 2
            sy2p = 1.0d0/2 * (-1.0d0/2 + dy) ** 2
            sy2m = 1.0d0/2 * (3.0d0/2 - dy) ** 2
         endif
         ipm = ip - 1
         ipp = ip + 1
         jpm = jp - 1
         jpp = jp + 1
         ippp= ip + 2
         jppp= jp + 2

         if( ip .eq. 0  ) ipm = nx - 1
         if( jp .eq. 0  ) jpm = ny - 1
         if( ip .eq. nx-1) ippp=1
         if( jp .eq. ny-1) jppp=1
         if(dx .le. 0.5d0 .and. dy .le. 0.5d0) then
            rho(ipm,jpm) = rho(ipm,jpm) + sx2m * sy2m * chrg
            rho(ipm,jp ) = rho(ipm,jp ) + sx2m * sy2  * chrg
            rho(ipm,jpp) = rho(ipm,jpp) + sx2m * sy2p * chrg
            rho(ip ,jpm) = rho(ip ,jpm) + sx2  * sy2m * chrg
            rho(ip ,jp ) = rho(ip ,jp ) + sx2  * sy2  * chrg
            rho(ip ,jpp) = rho(ip ,jpp) + sx2  * sy2p * chrg
            rho(ipp,jpm) = rho(ipp,jpm) + sx2p * sy2m * chrg
            rho(ipp,jp ) = rho(ipp,jp ) + sx2p * sy2  * chrg
            rho(ipp,jpp) = rho(ipp,jpp) + sx2p * sy2p * chrg
         else if(dx .le. 0.5d0 .and. dy .ge. 0.5d0) then
            rho(ipm ,jp  ) = rho(ipm ,jp  ) + sx2m * sy2m * chrg
            rho(ipm ,jpp ) = rho(ipm ,jpp ) + sx2m * sy2  * chrg
            rho(ipm ,jppp) = rho(ipm ,jppp) + sx2m * sy2p * chrg
            rho(ip  ,jp  ) = rho(ip  ,jp  ) + sx2  * sy2m * chrg
            rho(ip  ,jpp ) = rho(ip  ,jpp ) + sx2  * sy2  * chrg
            rho(ip  ,jppp) = rho(ip  ,jppp) + sx2  * sy2p * chrg
            rho(ipp ,jp  ) = rho(ipp ,jp  ) + sx2p * sy2m * chrg
            rho(ipp ,jpp ) = rho(ipp ,jpp ) + sx2p * sy2  * chrg
            rho(ipp ,jppp) = rho(ipp ,jppp) + sx2p * sy2p * chrg
         else if(dx .ge. 0.5d0 .and. dy .le. 0.5d0) then
            rho(ip  ,jpm) = rho(ip  ,jpm ) + sx2m * sy2m * chrg
            rho(ip  ,jp ) = rho(ip  ,jp  ) + sx2m * sy2  * chrg
            rho(ip  ,jpp) = rho(ip  ,jpp ) + sx2m * sy2p * chrg
            rho(ipp ,jpm) = rho(ipp ,jpm ) + sx2  * sy2m * chrg
            rho(ipp ,jp ) = rho(ipp ,jp  ) + sx2  * sy2  * chrg
            rho(ipp ,jpp) = rho(ipp ,jpp ) + sx2  * sy2p * chrg
            rho(ippp,jpm) = rho(ippp,jpm ) + sx2p * sy2m * chrg
            rho(ippp,jp ) = rho(ippp,jp  ) + sx2p * sy2  * chrg
            rho(ippp,jpp) = rho(ippp,jpp ) + sx2p * sy2p * chrg
         else 
            rho(ip  ,jp  ) = rho(ip  ,jp  ) + sx2m * sy2m * chrg
            rho(ip  ,jpp ) = rho(ip  ,jpp ) + sx2m * sy2  * chrg
            rho(ip  ,jppp) = rho(ip  ,jppp) + sx2m * sy2p * chrg
            rho(ipp ,jp  ) = rho(ipp ,jp  ) + sx2  * sy2m * chrg
            rho(ipp ,jpp ) = rho(ipp ,jpp ) + sx2  * sy2  * chrg
            rho(ipp ,jppp) = rho(ipp ,jppp) + sx2  * sy2p * chrg
            rho(ippp,jp  ) = rho(ippp,jp  ) + sx2p * sy2m * chrg
            rho(ippp,jpp ) = rho(ippp,jpp ) + sx2p * sy2  * chrg
            rho(ippp,jppp) = rho(ippp,jppp) + sx2p * sy2p * chrg
         endif
         end do
      !!..... set charge densities at the boundary
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
    subroutine phia(nx,ny,c,dt,phi,phib,jx,jy,jz,Ax,Ay,Az,Axb,Ayb,Azb,&
                    Axbb,Aybb,Azbb,cfact)
!***********************************************************************
   !original subroutine
      implicit none
      real(8), dimension(0:nx,0:ny) :: phi, phib, jx, jy, jz, Ax, Ay, Az,&
      Axb, Ayb, Azb, Axbb, Aybb, Azbb
      integer :: nx, ny, i, j, im, ip, jp, jm, jpm, ix ,iy
      real(rkind) :: c, dt, cfact

 ! Solution of maxwell equation in the A-phi formulation by difference method
 ! c is the ratio of the light speed to lattice parameter times plasma frequency
      do i = 0, nx
      do j = 0, ny

         im = i - 1
         ip = i + 1
         jm = j - 1
         jp = j + 1

         if( i .eq. 0  ) im = nx - 1
         if( i .eq. nx ) ip = 1
         if( j .eq. 0  ) jm = ny - 1
         if( j .eq. ny ) jp = 1
      
        Ax(i,j) = dt ** 2 * c ** 2 * (Axb(ip,j) + Axb(im,j) &
                + Axb(i,jp) + Axb(i,jm) - 4.0d0 * Axb(i,j)) &
                + dt ** 2 * jx(i,j) - 0.5d0 * dt * (phi(ip,j) - phib(ip,j)&
                - phi(im,j) + phib(im,j)) + 2.0d0 * Axb(i,j) - Axbb(i,j) 

        Ay(i,j) = dt ** 2 * c ** 2 * (Ayb(ip,j) + Ayb(im,j) &
                + Ayb(i,jp) + Ayb(i,jm) - 4.0d0 * Ayb(i,j)) &
                + dt ** 2 * jy(i,j) - 0.5d0 * dt * (phi(i,jp) - phib(i,jp)&
                - phi(i,jm) + phib(i,jm)) + 2.0d0 * Ayb(i,j) - Aybb(i,j) 

        Az(i,j) = dt ** 2 * c ** 2 * (Azb(ip,j) + Azb(im,j) &
                + Azb(i,jp) + Azb(i,jm) - 4.0d0 * Azb(i,j)) &
                + dt ** 2 * jz(i,j) + 2.0d0 * Azb(i,j) - Azbb(i,j)
      
      end do
      end do
      do i = 0, nx
      do j = 0, ny
        Axbb(i,j) = Axb(i,j)
        Aybb(i,j) = Ayb(i,j)
        Azbb(i,j) = Azb(i,j)
        Axb(i,j)  = Ax(i,j)
        Ayb(i,j)  = Ay(i,j)
        Azb(i,j)  = Az(i,j)
        phib(i,j) = phi(i,j)

      end do
   end do
   
      end subroutine phia                            

!***********************************************************************
   subroutine efield(nx,ny,phi,ex,ey,ez,Axb,Ayb,Azb,Axbb,Aybb,Azbb)
!***********************************************************************
     implicit none
      real(8), dimension(0:nx,0:ny) :: phi, ex, ey, ez, Axb, Ayb,&
           Azb, Axbb, Aybb, Azbb
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
         
         esx(i,j) = 0.5d0 * ( phi(im,j) - phi(ip,j))
         emx(i,j) = - ( Axb(i,j) - Axbb(i,j) ) / dt
         esy(i,j) = 0.5d0 * ( phi(i,jm) - phi(i,jp))
         emy(i,j) = - ( Ayb(i,j) - Aybb(i,j) ) / dt
         esz(i,j) = 0.d0
         emz(i,j) = - ( Azb(i,j) - Azbb(i,j) ) / dt
         ex(i,j) = esx(i,j) + emx(i,j)
         ey(i,j) = esy(i,j) + emy(i,j)
         ez(i,j) = esz(i,j) + emz(i,j)

      end do
   end do
    end subroutine efield
!***********************************************************************
    subroutine bfield(nx,ny,bxmin,bxmax,bymin,bymax,bzmin,bzmax,bx,by,bz,&
         Ax,Ay,Az)
!***********************************************************************
      implicit none
      real(rkind), dimension(0:ny) :: bxnab,bznab
      real(rkind), dimension(0:nx) :: bynab
      real(rkind), dimension(0:nx,0:ny) :: bx, by, bz, Ax, Ay, Az
      integer :: nx, ny, i, j, ip, jp, im, jm
      real(rkind)::bxmin, bxmax, bymin, bymax, bzmin, bzmax
      do i=0,nx-1
         bynab(i) = bymin + dble(i) / dble(nx) * (bymax - bymin)
      end do
      bynab(nx) = bymin
      !bxnab(0) = bxnab(nx) = bxmin
      !bxnab(nx-1) = bxmax
      !do j=0,ny-1
      !   bxnab(j) = j / ny * 
      !enddo
      !do i=0,nx-1
      !   bznab(i) = dble(i) / dble(nx) * bzbg
      !enddo

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
         
         bx(i,j) = 0.25d0 * (Azb(i,jp) + Azbb(i,jp) - Azb(i,jm) - Azbb(i,jm))
         by(i,j) = 0.25d0 * (Azb(ip,j) + Azbb(ip,j) - Azb(im,j) - Azbb(im,j)) &
                 + bynab(i)
         bz(i,j) = 0.25d0 * (Ayb(ip,j) + Aybb(ip,j) - Ayb(im,j)-Aybb(im,j) &
                 - (Axb(i,jp) + Axbb(i,jp)-Axb(i,jm)-Axbb(i,jm)))
      end do
      end do
    !  j=ny
    !  do i = 0, nx-1
    !     bx(i,j) = 0.25d0 * (Azb(i,1)+Azbb(i,1)-Azb(i,j)-Azbb(i,j))
    !     by(i,j) = bynab(i) + 0.25d0 * (Azb(i+1,j)+Azbb(i+1,j) &
    !             - Azb(i,j)-Azbb(i,j))
    !     bz(i,j) = 0.25d0 * (Ayb(i+1,j) + Aybb(i+1,j)-Ayb(i,j)-Aybb(i,j) &
    !             - (Axb(i,1)+Aybb(i,1)-Axb(i,j)-Axbb(i,j)))
    !  end do
     
    !  i=nx
    !  do j = 0, ny-1
    !     bx(i,j) = 0.25d0 * (Azb(i,j+1)+Azbb(i,j+1) &
    !             - Azb(i,j)-Azbb(i,j))
    !     by(i,j) = bynab(i) + 0.25d0 * (Azb(1,j)+Azbb(1,j)-Azb(i,j)-Azbb(i,j))
    !     bz(i,j) = 0.25d0 * (Ayb(1,j)+Aybb(1,j)- Ayb(i,j)-Aybb(i,j) &
    !             - (Axb(i,j+1)+Axbb(i,j+1)-Axb(i,j)-Axbb(i,j)))
    !  end do

    !  i=nx
    !  j=ny
    !     bx(i,j) = 0.25d0 * (Azb(i,1)+Azbb(i,1)-Azb(i,j)-Azbb(1,j))
    !     by(i,j) = bynab(i) - 0.25d0 * (Azb(1,j)+Azbb(1,j)-Azb(i,j)-Azbb(i,j))
    !     bz(i,j) = 0.25d0 * (Ayb(1,j)+Aybb(1,j)-Ayb(i,j)-Aybb(i,j) &
    !             - (Axb(i,1)+Axb(i,1)-Axb(i,j)-Axbb(i,j)))
    !  do j = 0, ny
    !  do i = 0, nx
    !     bb(i,j) = SQRT(bx(i,j)**2+by(i,j)**2+bz(i,j)**2)
    !     AA(i,j) = SQRT(Ax(i,j)**2+Ay(i,j)**2+Az(i,j)**2)
    !  end do
    !  end do

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

         rho(i,j) = phi(im,j) + phi(ip,j) + phi(i,jm) + phi(i,jp) &
                  - 4.0d0 * phi(i,j)

      end do
      end do

    end subroutine d2phi

!***********************************************************************
    subroutine kine(np,vx,vy,vz,akin,mass,cfacti)
!***********************************************************************
      implicit none
      real(8), dimension(np) :: vx, vy, vz
      real(8) :: akin, mass ,cfacti
      integer(4) :: np, i
      akin = 0.d0
      do i = 1, np
         akin = akin + vx(i)**2 + vy(i)**2 + vz(i)**2
      end do

      akin = 0.5 * akin * mass
    end subroutine kine

!***********************************************************************
    subroutine pote(nx,ny,ex,ey,ez,bx,by,bz,c,apot,cfacti)
!***********************************************************************
      implicit none
      real(8), dimension(0:nx,0:ny) :: ex, ey, ez, bx, by, bz
      real(8) :: apot, cfacti,c
      integer(4) :: nx, ny, i, j
      apot = 0.d0

      do j = 0, ny-1
      do i = 0, nx-1

         apot = apot + ex(i,j)**2 + ey(i,j)**2 + ez(i,j)**2 &
              + bx(i,j)**2 + by(i,j)**2 + bz(i,j)** 2
            
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

!***********************************************************************
    subroutine current(np,nx,ny,x,y,vx,vy,vz,jx,jy,jz,dt,chrg,cfact,&
         jxbg,jybg,jzbg,omega,thetax,thetay,thetaz,time)
!***********************************************************************
    implicit none

    real(8), dimension(np) :: x, y, vx, vy, vz 
    real(8), dimension(0:nx,0:ny) :: jx, jy, jz
    real(8) :: chrg, dt, dx, dy, dx1, dy1, &
         sx1p, sy1p, sx1m, sy1m, sx2, sy2, sx2p, sy2p, sx2m, sy2m, cfact,time,&
         jxbg, jybg, jzbg,omega,thetax,thetay,thetaz,pi,twopi,jxbgg,jybgg,jzbgg
    integer :: np, nx, ny, i, ip, jp, ipm, jpm, ipp, jpp, ippp, jppp, ix, iy

     do i = 1, np
        
        ip = x(i)
        jp = y(i)
        dx = x(i) - dble(ip)
        dy = y(i) - dble(jp)
        dx1 = 1.0d0 - dx
        dy1 = 1.0d0 - dy 
        if(dx .le. 0.5d0) then
            sx2  = 3.0d0/4 - dx ** 2
            sx2p = 1.0d0/2 * (1.0d0/2 + dx) ** 2
            sx2m = 1.0d0/2 * (1.0d0/2 - dx) ** 2
         else
            sx2  = 3.0d0/4 - (dx - 1.0d0) ** 2
            sx2p = 1.0d0/2 * (-1.0d0/2 + dx) ** 2
            sx2m = 1.0d0/2 * (3.0d0/2 - dx) ** 2
         endif
         if(dy .le. 0.5d0) then
            sy2  = 3.0d0/4 - dy ** 2
            sy2p = 1.0d0/2 * (1.0d0/2 + dy) ** 2
            sy2m = 1.0d0/2 * (1.0d0/2 - dy) ** 2
         else
            sy2  = 3.0d0/4 - (dy - 1.0d0) ** 2
            sy2p = 1.0d0/2 * (-1.0d0/2 + dy) ** 2
            sy2m = 1.0d0/2 * (3.0d0/2 - dy) ** 2
         endif
         ipm = ip - 1
         ipp = ip + 1
         jpm = jp - 1
         jpp = jp + 1
         ippp = ip + 2
         jppp = jp + 2
         if( ip .eq. 0  ) ipm = nx-1
         if( jp .eq. 0  ) jpm = ny-1
         if( ip .eq. nx-1  ) ippp = 1
         if( jp .eq. ny-1  ) jppp = 1
        
         if (dx .le. 0.5d0 .and. dy .le. 0.5d0) then
           jx(ipp,jp ) = jx(ipp,jp ) + chrg * vx(i) * sy2  * dx
           jx(ipp,jpp) = jx(ipp,jpp) + chrg * vx(i) * sy2p * dx
           jx(ipp,jpm) = jx(ipp,jpm) + chrg * vx(i) * sy2m * dx
           jx(ip ,jp ) = jx(ip ,jp ) + chrg * vx(i) * sy2  * dx1
           jx(ip ,jpp) = jx(ip ,jpp) + chrg * vx(i) * sy2p * dx1
           jx(ip ,jpm) = jx(ip ,jpm) + chrg * vx(i) * sy2m * dx1
           
           jy(ip ,jpp) = jy(ip ,jpp) + chrg * vy(i) * sx2  * dy
           jy(ipp,jpp) = jy(ipp,jpp) + chrg * vy(i) * sx2p * dy
           jy(ipm,jpp) = jy(ipm,jpp) + chrg * vy(i) * sx2m * dy
           jy(ip ,jp ) = jy(ip ,jp ) + chrg * vy(i) * sx2  * dy1
           jy(ipp,jp ) = jy(ipp,jp ) + chrg * vy(i) * sx2p * dy1
           jy(ipm,jp ) = jy(ipm,jp ) + chrg * vy(i) * sx2m * dy1
      
           jz(ipm,jpp) = jz(ipm,jpp) + chrg * vz(i) * sx2m * sy2p
           jz(ipm,jp ) = jz(ipm,jp ) + chrg * vz(i) * sx2m * sy2
           jz(ipm,jpm) = jz(ipm,jpm) + chrg * vz(i) * sx2m * sy2m
           jz(ip ,jpp) = jz(ip ,jpp) + chrg * vz(i) * sx2  * sy2p
           jz(ip ,jp ) = jz(ip ,jp ) + chrg * vz(i) * sx2  * sy2
           jz(ip ,jpm) = jz(ip ,jpm) + chrg * vz(i) * sx2  * sy2m
           jz(ipp,jpp) = jz(ipp,jpp) + chrg * vz(i) * sx2p * sy2p
           jz(ipp,jp ) = jz(ipp,jp ) + chrg * vz(i) * sx2p * sy2
           jz(ipp,jpm) = jz(ipp,jpm) + chrg * vz(i) * sx2p * sy2m

        else if(dx .le. 0.5d0 .and. dy .ge. 0.5d0) then
           jx(ipp,jpp ) = jx(ipp,jpp ) + chrg * vx(i) * sy2  * dx
           jx(ipp,jppp) = jx(ipp,jppp) + chrg * vx(i) * sy2p * dx
           jx(ipp,jp  ) = jx(ipp,jp  ) + chrg * vx(i) * sy2m * dx
           jx(ip ,jpp ) = jx(ip ,jpp ) + chrg * vx(i) * sy2  * dx1
           jx(ip ,jppp) = jx(ip ,jppp) + chrg * vx(i) * sy2p * dx1
           jx(ip ,jp  ) = jx(ip ,jp  ) + chrg * vx(i) * sy2m * dx1
           
           jy(ip ,jpp) = jy(ip ,jpp) + chrg * vy(i) * sx2  * dy
           jy(ipp,jpp) = jy(ipp,jpp) + chrg * vy(i) * sx2p * dy
           jy(ipm,jpp) = jy(ipm,jpp) + chrg * vy(i) * sx2m * dy
           jy(ip ,jp ) = jy(ip ,jp ) + chrg * vy(i) * sx2  * dy1
           jy(ipp,jp ) = jy(ipp,jp ) + chrg * vy(i) * sx2p * dy1
           jy(ipm,jp ) = jy(ipm,jp ) + chrg * vy(i) * sx2m * dy1
      
           jz(ipm,jppp) = jz(ipm,jppp) + chrg * vz(i) * sx2m * sy2p
           jz(ipm,jpp ) = jz(ipm,jpp ) + chrg * vz(i) * sx2m * sy2
           jz(ipm,jp  ) = jz(ipm,jp  ) + chrg * vz(i) * sx2m * sy2m
           jz(ip ,jppp) = jz(ip ,jppp) + chrg * vz(i) * sx2  * sy2p
           jz(ip ,jpp ) = jz(ip ,jpp ) + chrg * vz(i) * sx2  * sy2
           jz(ip ,jp  ) = jz(ip ,jp  ) + chrg * vz(i) * sx2  * sy2m
           jz(ipp,jppp) = jz(ipp,jppp) + chrg * vz(i) * sx2p * sy2p
           jz(ipp,jpp ) = jz(ipp,jpp ) + chrg * vz(i) * sx2p * sy2
           jz(ipp,jp  ) = jz(ipp,jp  ) + chrg * vz(i) * sx2p * sy2m

        else if(dx .ge. 0.5d0 .and. dy .le. 0.5d0) then
           jx(ipp,jp ) = jx(ipp,jp ) + chrg * vx(i) * sy2  * dx
           jx(ipp,jpp) = jx(ipp,jpp) + chrg * vx(i) * sy2p * dx
           jx(ipp,jpm) = jx(ipp,jpm) + chrg * vx(i) * sy2m * dx
           jx(ip ,jp ) = jx(ip ,jp ) + chrg * vx(i) * sy2  * dx1
           jx(ip ,jpp) = jx(ip ,jpp) + chrg * vx(i) * sy2p * dx1
           jx(ip ,jpm) = jx(ip ,jpm) + chrg * vx(i) * sy2m * dx1
           
           jy(ipp ,jpp) = jy(ipp ,jpp) + chrg * vy(i) * sx2  * dy
           jy(ippp,jpp) = jy(ippp,jpp) + chrg * vy(i) * sx2p * dy
           jy(ip  ,jpp) = jy(ip  ,jpp) + chrg * vy(i) * sx2m * dy
           jy(ipp ,jp ) = jy(ipp ,jp ) + chrg * vy(i) * sx2  * dy1
           jy(ippp,jp ) = jy(ippp,jp ) + chrg * vy(i) * sx2p * dy1
           jy(ip  ,jp ) = jy(ip  ,jp ) + chrg * vy(i) * sx2m * dy1
      
           jz(ip  ,jpp) = jz(ip  ,jpp) + chrg * vz(i) * sx2m * sy2p
           jz(ip  ,jp ) = jz(ip  ,jp ) + chrg * vz(i) * sx2m * sy2
           jz(ip  ,jpm) = jz(ip  ,jpm) + chrg * vz(i) * sx2m * sy2m
           jz(ipp ,jpp) = jz(ipp ,jpp) + chrg * vz(i) * sx2  * sy2p
           jz(ipp ,jp ) = jz(ipp ,jp ) + chrg * vz(i) * sx2  * sy2
           jz(ipp ,jpm) = jz(ipp ,jpm) + chrg * vz(i) * sx2  * sy2m
           jz(ippp,jpp) = jz(ippp,jpp) + chrg * vz(i) * sx2p * sy2p
           jz(ippp,jp ) = jz(ippp,jp ) + chrg * vz(i) * sx2p * sy2
           jz(ippp,jpm) = jz(ippp,jpm) + chrg * vz(i) * sx2p * sy2m

        else
           jx(ipp,jpp ) = jx(ipp,jpp ) + chrg * vx(i) * sy2  * dx
           jx(ipp,jppp) = jx(ipp,jppp) + chrg * vx(i) * sy2p * dx
           jx(ipp,jp  ) = jx(ipp,jp  ) + chrg * vx(i) * sy2m * dx
           jx(ip ,jpp ) = jx(ip ,jpp ) + chrg * vx(i) * sy2  * dx1
           jx(ip ,jppp) = jx(ip ,jppp) + chrg * vx(i) * sy2p * dx1
           jx(ip ,jp  ) = jx(ip ,jp  ) + chrg * vx(i) * sy2m * dx1
           
           jy(ipp ,jpp) = jy(ipp ,jpp) + chrg * vy(i) * sx2  * dy
           jy(ippp,jpp) = jy(ippp,jpp) + chrg * vy(i) * sx2p * dy
           jy(ip  ,jpp) = jy(ip  ,jpp) + chrg * vy(i) * sx2m * dy
           jy(ipp ,jp ) = jy(ipp ,jp ) + chrg * vy(i) * sx2  * dy1
           jy(ippp,jp ) = jy(ippp,jp ) + chrg * vy(i) * sx2p * dy1
           jy(ip  ,jp ) = jy(ip  ,jp ) + chrg * vy(i) * sx2m * dy1
      
           jz(ip  ,jppp) = jz(ip  ,jppp) + chrg * vz(i) * sx2m * sy2p
           jz(ip  ,jpp ) = jz(ip  ,jpp ) + chrg * vz(i) * sx2m * sy2
           jz(ip  ,jp  ) = jz(ip  ,jp  ) + chrg * vz(i) * sx2m * sy2m
           jz(ipp ,jppp) = jz(ipp ,jppp) + chrg * vz(i) * sx2  * sy2p
           jz(ipp ,jpp ) = jz(ipp ,jpp ) + chrg * vz(i) * sx2  * sy2
           jz(ipp ,jp  ) = jz(ipp ,jp  ) + chrg * vz(i) * sx2  * sy2m
           jz(ippp,jppp) = jz(ippp,jppp) + chrg * vz(i) * sx2p * sy2p
           jz(ippp,jpp ) = jz(ippp,jpp ) + chrg * vz(i) * sx2p * sy2
           jz(ippp,jp  ) = jz(ippp,jp  ) + chrg * vz(i) * sx2p * sy2m

        endif
       
        end do
      if( chrg .gt. 0.d0 ) then
         jxbgg = jxbg * cos (omega * time + thetax)
         jybgg = jybg * cos (omega * time + thetay)
         jzbgg = jzbg * cos (omega * time + thetaz)
         ! inserting background current density
         do iy=1,ny
            jy(0,iy) = jx(ix,0) + 0.5d0 * jybgg 
            jy(1,iy) = jx(ix,1) + 0.5d0 * jybgg
         end do
         !--------------
         do iy = 0, ny
            jx(0,iy) = jx(0,iy) + jx(nx,iy)
            jy(0,iy) = jy(0,iy) + jy(nx,iy)
            jz(0,iy) = jz(0,iy) + jz(nx,iy)
            jx(nx,iy) = jx(0,iy)
            jy(nx,iy) = jy(0,iy)
            jz(nx,iy) = jz(0,iy)
         end do
         do ix = 0, nx
            jx(ix,0) = jx(ix,0) + jx(ix,ny)
            jy(ix,0) = jy(ix,0) + jy(ix,ny) 
            jz(ix,0) = jz(ix,0) + jz(ix,ny)
            jx(ix,ny) = jx(ix,0)
            jy(ix,ny) = jy(ix,0)
            jz(ix,ny) = jz(ix,0)
         enddo
         do iy = 0, ny
         do ix = 0, nx
            jx(ix,iy) = cfact * jx(ix,iy)
            jy(ix,iy) = cfact * jy(ix,iy)
            jz(ix,iy) = cfact * jz(ix,iy)
         end do
         end do
      endif
   end subroutine current
END Module picexec
