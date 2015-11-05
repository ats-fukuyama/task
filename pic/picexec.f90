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
    Axb(:,:) = 0.0d0
    Ayb(:,:) = 0.0d0
    Azb(:,:) = 0.0d0
    Axbb(:,:) = 0.0d0
    Aybb(:,:) = 0.0d0
    Azbb(:,:) = 0.0d0
    phib(:,:) = 0.0d0
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
         call current(np,nx,ny,nz,xe,ye,ze,vxe,vye,vze,jx,jy,jz,dt,chrge,cfact)
                  !.......... calculate current by ions
         call current(np,nx,ny,nz,xi,yi,zi,vxi,vyi,vzi,jx,jy,jz,dt,chrgi,cfact)
         
         !.......... calculate vector potential
         call phia(nx,ny,c,omega,dt,phi,phib,jx,jy,jz,Ax,Ay,Az,Axb,Ayb,Azb,&
              Axbb,Aybb,Azbb,cfact)

           !.......... calculate ex and ey and ez
         call efield(nx,ny,nz,phi,ex,ey,ez,Axb,Ayb,Azb,Axbb,Aybb,Azbb)

           !.......... calculate bxg and byg and bzg
         call bfield(nx,ny,nz,bx,by,bz,bxg,byg,bzg,Axb,Ayb,Azb)

          if( mod(iloop,nhmod) .eq. 0 ) then
            iene = iene + 1
            call kine(np,vxe,vye,vze,akine1,me)
            call kine(np,vxi,vyi,vzi,akini1,mi)
            call pote(nx,ny,nz,phi,ex,ey,ez,Ax,Ay,Az,apot,cfacti)
            call sumdim1(nodes,myid,akine1,wkword)
            call sumdim1(nodes,myid,akini1,wkword)
         endif

         !..... push electrons
         call push(np,nx,ny,nz,xe,ye,ze,vxe,vye,vze,ex,ey,ez,bxg,byg,bzg,dt,&
              ctome,xeb,yeb,zeb,phi,phib,Axb,Ayb,Azb,Axbb,Aybb,Azbb, &
              vparae,vperpe)
         
         !..... push ions
         call push(np,nx,ny,nz,xi,yi,zi,vxi,vyi,vzi,ex,ey,ez,bxg,byg,bzg,dt,&
              ctomi,xib,yib,zib,phi,phib,Axb,Ayb,Azb,Axbb,Aybb,Azbb, &
              vparai,vperpi)

         !----- treat particles being out of the boundary
         call bound(np,xe,ye,ze,x1,x2,y1,y2,z1,z2,alx,aly,alz)
         call bound(np,xi,yi,zi,x1,x2,y1,y2,z1,z2,alx,aly,alz)

         !..... diagnostics to check energy conservation
         !.....            after pushing 
         if( mod(iloop,nhmod) .eq. 0 ) then
            call kine(np,vxe,vye,vze,akine2,me)
            call kine(np,vxi,vyi,vzi,akini2,mi)
            call sumdim1(nodes,myid,akine2,wkword)
            call sumdim1(nodes,myid,akini2,wkword)
            !apot=0.d0
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
       subroutine push(np,nx,ny,nz,x,y,z,vx,vy,vz,ex,ey,ezg,bxg,byg,bzg,dt,&
            ctom,xb,yb,zb,phi,phib,Axb,Ayb,Azb,Axbb,Aybb,Azbb,vpara,vperp)
!***********************************************************************
      implicit none
      real(8), dimension(np) :: x, y, z, xb, yb, zb, vx, vy, vz, vpara, vperp
      real(8), dimension(0:nx,0:ny) :: phi,phib,Axb, Ayb, Azb, Axbb, Aybb, Azbb
      real(8), dimension(0:nx,0:ny,0:nz) :: ex, ey, ezg, bxg, byg, bzg
      real(8) :: ctom, dx, dy, dz, dx1, dy1, dz1, dt, exx, eyy, ezz, bxx,&
                 byy, bzz, vxn, vyn, vzn, vxzero, vyzero, vzzero, vxp, vyp,&
                 vzp, sx1p, sx1m, sy1p, sy1m, sx2, sy2, sx2p, sx2m, sy2m, sy2p
      real(8) :: btot, vtot
      integer :: np, nx, ny, nz, i, j, ip, jp, kp, ipp, ipm, jpp, jpm
     
      do i = 1, np

! calculate the electric field at the particle position

         ip = x(i)
         jp = y(i)
         kp = z(i)
         dx = x(i) - dble(ip)
         dy = y(i) - dble(jp)
         !dz = z(i) - dble(kp)
         dx1 = 1.0d0 - dx
         dy1 = 1.0d0 - dy
         !dz1 = 1.0d0 - dz

         if (dx .le. 0.5d0) then 
            sx1p = 0.5d0 + dx
            sx1m = 0.5d0 - dx
         else
            sx1p = dx - 0.5d0
            sx1m = 1.5d0 - dx
         end if
         if (dy .le. 0.5d0) then 
            sy1p = 0.5d0 + dy
            sy1m = 0.5d0 - dy
         else
            sy1p = dy - 0.5d0
            sy1m = 1.5d0 - dy
         end if

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

         if( ip .eq. 0  ) ipm = nx - 1
         if( ip .eq. nx ) ipp = 1
         if( jp .eq. 0  ) jpm = ny - 1
         if( jp .eq. ny ) jpp = 1       
         
         ! electric field
         exx = ex(ip ,jp  ,0)*dx1*dy1 + ex(ip+1,jp  ,0)*dx*dy1 &
             + ex(ip ,jp+1,0)*dx1*dy  + ex(ip+1,jp+1,0)*dx*dy  
              
             !  ex(ipp,jpp,0) * sx2p * sy2p &
             !+ ex(ipp,jp ,0) * sx2p * sy2  &
             !+ ex(ipp,jpm,0) * sx2p * sy2m &
             !+ ex(ip ,jpp,0) * sx2  * sy2p &
             !+ ex(ip ,jp ,0) * sx2  * sy2  &
             !+ ex(ip ,jpm,0) * sx2  * sy2m &
             !+ ex(ipm,jpp,0) * sx2m * sy2p &
             !+ ex(ipm,jp ,0) * sx2m * sy2  &
             !+ ex(ipm,jpm,0) * sx2m * sy2m
         
         !- phi(ipp,jpp) * sy2p * sx1p &
                !- phi(ipp,jp ) * sy2  * sx1p &
                !- phi(ipp,jpm) * sy2m * sx1p &
                !+ phi(ip ,jpp) * sy2p * (sx1p - sx1m) &
                !+ phi(ip ,jp ) * sy2  * (sx1p - sx1m) &
                !+ phi(ip ,jpm) * sy2m * (sx1p - sx1m) &
                !+ phi(ipm,jpp) * sy2p * sx1m &
                !+ phi(ipm,jp ) * sy2  * sx1m &
                !+ phi(ipm,jpm) * sy2m * sx1m 
                !- (Axb(ip ,jpp) - Axbb(ip ,jpp)) / dt * sx1p * sy2p &
                !- (Axb(ip ,jp ) - Axbb(ip ,jp )) / dt * sx1p * sy2  &
                !- (Axb(ip ,jpm) - Axbb(ip ,jpm)) / dt * sx1p * sy2m &
                !- (Axb(ipm,jpp) - Axbb(ipm,jpp)) / dt * sx1m * sy2p &
                !- (Axb(ipm,jp ) - Axbb(ipm,jp )) / dt * sx1m * sy2  &
                !- (Axb(ipm,jpm) - Axbb(ipm,jpm)) / dt * sx1m * sy2m 
         
         eyy = ey(ip ,jp  ,0)*dx1*dy1 + ey(ip+1,jp  ,0)*dx*dy1 &
             + ey(ip ,jp+1,0)*dx1*dy  + ey(ip+1,jp+1,0)*dx*dy  

             !  ey(ipp,jpp,0) * sx2p * sy2p &
             !+ ey(ip ,jpp,0) * sx2  * sy2p &
             !+ ey(ipm,jpp,0) * sx2m * sy2p &
             !+ ey(ipp,jp ,0) * sx2p * sy2  &
             !+ ey(ip ,jp ,0) * sx2  * sy2  &
             !+ ey(ipm,jp ,0) * sx2m * sy2  &
             !+ ey(ipp,jpm,0) * sx2p * sy2m &
             !+ ey(ip ,jpm,0) * sx2  * sy2m &
             !+ ey(ipm,jpm,0) * sx2m * sy2m  
             
         !- phi(ipp,jpp) * sx2p * sy1p &
                !- phi(ip ,jpp) * sx2  * sy1p &
                !- phi(ipm,jpp) * sx2m * sy1p &
                !+ phi(ipp,jp ) * sx2p * (sy1p - sy1m) &
                !+ phi(ip ,jp ) * sx2  * (sy1p - sy1m) &
                !+ phi(ipm,jp ) * sx2m * (sy1p - sy1m) &
                !+ phi(ipp,jpm) * sx2p * sy1m &
                !+ phi(ip ,jpm) * sx2  * sy1m &
                !+ phi(ipm,jpm) * sx2m * sy1m 
                !- (Ayb(ipp,jp ) - Aybb(ipp,jp )) / dt * sx2p * sy1p &
                !- (Ayb(ip ,jp ) - Aybb(ip ,jp )) / dt * sx2  * sy1p &
                !- (Ayb(ipm,jp ) - Aybb(ipm,jp )) / dt * sx2m * sy1p &
                !- (Ayb(ipp,jpm) - Aybb(ipp,jpm)) / dt * sx2p * sy1m &
                !- (Ayb(ip ,jpm) - Aybb(ip ,jpm)) / dt * sx2  * sy1m &
                !- (Ayb(ipm,jpm) - Aybb(ipm,jpm)) / dt * sx2m * sy1m      
         
             ezz = ezg(ip ,jp  ,0)*dx1*dy1 + ezg(ip+1,jp  ,0)*dx*dy1 &
                 + ezg(ip ,jp+1,0)*dx1*dy  + ezg(ip+1,jp+1,0)*dx*dy  
         !- (Azb(ipp,jpp) - Azbb(ipp,jpp)) * sx2p * sy2p / dt &
                !- (Azb(ipp,jp ) - Azbb(ipp,jp )) * sx2p * sy2  / dt &
                !- (Azb(ipp,jpm) - Azbb(ipp,jpm)) * sx2p * sy2m / dt &
                !- (Azb(ip ,jpp) - Azbb(ip ,jpp)) * sx2  * sy2p / dt &
                !- (Azb(ip ,jp ) - Azbb(ip ,jp )) * sx2  * sy2  / dt &
                !- (Azb(ip ,jpm) - Azbb(ip ,jpm)) * sx2  * sy2m / dt &
                !- (Azb(ipm,jpp) - Azbb(ipm,jpp)) * sx2m * sy2p / dt &
                !- (Azb(ipm,jp ) - Azbb(ipm,jp )) * sx2m * sy2  / dt &
                !- (Azb(ipm,jpm) - Azbb(ipm,jpm)) * sx2m * sy2m / dt         
            
            !exx=0.d0
            !eyy=0.d0
            !ezz=0.d0

           ! magnetic field
             bxx = bxg(ip  ,jp,  0)*dx1*dy1 + bxg(ip+1,jp  ,0)*dx*dy1 &
                 + bxg(ip  ,jp+1,0)*dx1*dy  + bxg(ip+1,jp+1,0)*dx*dy
             
                !(Azb(ipp,jpp) + Azbb(ipp,jpp))/2 * sx2p * sy1p &
                !+ (Azb(ip ,jpp) + Azbb(ip ,jpp))/2 * sx2  * sy1p &
                !+ (Azb(ipm,jpp) + Azbb(ipm,jpp))/2 * sx2m * sy1p &
                !- (Azb(ipp,jp ) + Azbb(ipp,jp ))/2 * sx2p * (sy1p-sy1m)&
                !- (Azb(ip ,jp ) + Azbb(ip ,jp ))/2 * sx2  * (sy1p-sy1m)&
                !- (Azb(ipm,jp ) + Azbb(ipm,jp ))/2 * sx2m * (sy1p-sy1m) &
                !- (Azb(ipp,jpm) + Azbb(ipp,jpm))/2 * sx2p * sy1m &
                !- (Azb(ip ,jpm) + Azbb(ip ,jpm))/2 * sx2  * sy1m &
                !- (Azb(ipm,jpm) + Azbb(ipm,jpm))/2 * sx2m * sy1m 

             byy = byg(ip  ,jp,  0)*dx1*dy1 + byg(ip+1,jp  ,0)*dx*dy1 &
                 + byg(ip  ,jp+1,0)*dx1*dy  + byg(ip+1,jp+1,0)*dx*dy
                !- (Azb(ipp,jpp) + Azbb(ipp,jpp))/2 * sy2p * sx1p &
                !- (Azb(ipp,jp ) + Azbb(ipp,jp ))/2 * sy2  * sx1p &
                !- (Azb(ipp,jpm) + Azbb(ipp,jpm))/2 * sy2m * sx1p &
                !+ (Azb(ip ,jpp) + Azbb(ip ,jpp))/2 * sy2p * (sx1p-sx1m) &
                !+ (Azb(ip ,jp ) + Azbb(ip ,jp ))/2 * sy2  * (sx1p-sx1m) &
                !+ (Azb(ip ,jpm) + Azbb(ip ,jpm))/2 * sy2m * (sx1p-sx1m) &
                !+ (Azb(ipm,jpp) + Azbb(ipm,jpp))/2 * sy2p * sx1m &
                !+ (Azb(ipm,jp ) + Azbb(ipm,jp ))/2 * sy2  * sx1m &
                !+ (Azb(ipm,jpm) + Azbb(ipm,jpm))/2 * sy2m * sx1m
        
             bzz = bzg(ip  ,jp,  0)*dx1*dy1 + bzg(ip+1,jp  ,0)*dx*dy1 &
                 + bzg(ip  ,jp+1,0)*dx1*dy  + bzg(ip+1,jp+1,0)*dx*dy
                !(Ayb(ipp,jp ) + Aybb(ipp,jp )) / 2 * sy1p * sx1p &
                !- (Ayb(ip ,jp ) + Aybb(ip ,jp )) / 2 * sy1p * (sx1p-sx1m) &
                !- (Ayb(ipm,jp ) + Aybb(ipm,jp )) / 2 * sy1p * sx1m &
                !+ (Ayb(ipp,jpm) + Aybb(ipp,jpm)) / 2 * sy1m * sx1p &
                !- (Ayb(ip ,jpm) + Aybb(ip ,jpm)) / 2 * sy1m * (sx1p-sx1m) &
                !- (Ayb(ipm,jpm) + Aybb(ipm,jpm)) / 2 * sy1m * sx1m &
                !- (Axb(ip ,jpp) + Axbb(ip ,jpp)) / 2 * sx1p * sy1p &
                !+ (Axb(ip ,jp ) + Axbb(ip ,jp )) / 2 * sx1p * (sy1p-sy1m) &
                !+ (Axb(ip ,jpm) + Axbb(ip ,jpm)) / 2 * sx1p * sy1m &
                !- (Axb(ipm,jpp) + Axbb(ipm,jpp)) / 2 * sx1m * sy1p &
                !+ (Axb(ipm,jp ) + Axbb(ipm,jp )) / 2 * sx1m * (sy1p-sy1m) &
                !+ (Axb(ipm,jpm) + Axbb(ipm,jpm)) / 2 * sx1m * sy1m 
                
            ! push particles by using Buneman-Boris method
          
         vxn = vx(i) + 1.0d0/2 * ctom * exx * dt 
         vyn = vy(i) + 1.0d0/2 * ctom * eyy * dt
         vzn = vz(i) + 1.0d0/2 * ctom * ezz * dt

         vxzero = vxn + 1.0d0/2 * ctom * (vyn * bzz - vzn * byy) * dt
         vyzero = vyn + 1.0d0/2 * ctom * (vzn * bxx - vxn * bzz) * dt 
         vzzero = vzn + 1.0d0/2 * ctom * (vxn * byy - vyn * bxx) * dt

         vxp = vxn + 2.0d0/(1.0d0 + 0.25d0 * (ctom * dt) ** 2 & 
             * (bxx ** 2 + byy ** 2 + bzz ** 2)) * 1.0d0/2 & 
             * ctom * (vyzero * bzz - vzzero * byy) * dt
         vyp = vyn + 2.0d0/(1.0d0 + 0.25d0 * (ctom * dt) ** 2 &
             * (bxx ** 2 + byy ** 2 + bzz ** 2)) * 1.0d0/2 &
             * ctom * (vzzero * bxx - vxzero * bzz) * dt 
         vzp = vzn + 2.0d0/(1.0d0 + 0.25d0 * (ctom * dt) ** 2 & 
             * (bxx ** 2 + byy ** 2 + bzz ** 2)) * 1.0d0/2 & 
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
      integer :: np, nx, ny, i, ip, jp, kp, ix, iy, ipp, ipm, jpp, jpm

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
         if(dx .le. 0.5d0) then
            sx2  = 3.0d0/4 - dx ** 2
            sx2p = 1.0d0/2 * (1.0d0/2 + dx) ** 2
            sx2m = 1.0d0/2 * (1.0d0/2 - dx) ** 2
         else
            sx2  = 3.0d0/4 - (dx - 1.0d0) ** 2
            sx2p = 1.0d0/2 * (3.0d0/2 - dx) ** 2
            sx2m = 1.0d0/2 * (-1.0d0/2 + dx) ** 2
         endif
         if(dy .le. 0.5d0) then
            sy2  = 3.0d0/4 - dy ** 2
            sy2p = 1.0d0/2 * (1.0d0/2 + dy) ** 2
            sy2m = 1.0d0/2 * (1.0d0/2 - dy) ** 2
         else
            sy2  = 3.0d0/4 - (dy - 1.0d0) ** 2
            sy2p = 1.0d0/2 * (3.0d0/2 - dy) ** 2
            sy2m = 1.0d0/2 * (-1.0d0/2 + dy) ** 2
         endif
         ipm = ip - 1
         ipp = ip + 1
         jpm = jp - 1
         jpp = jp + 1

         if( ip .eq. 0  ) ipm = nx - 1
         if( ip .eq. nx ) ipp = 1
         if( jp .eq. 0  ) jpm = ny - 1
         if( jp .eq. ny ) jpp = 1

         rho(ip  ,jp  ) = rho(ip  ,jp  ) + dx1 * dy1 * chrg
         rho(ip+1,jp  ) = rho(ip+1,jp  ) + dx  * dy1 * chrg
         rho(ip  ,jp+1) = rho(ip  ,jp+1) + dx1 * dy  * chrg
         rho(ip+1,jp+1) = rho(ip+1,jp+1) + dx  * dy  * chrg

           ! rho(ipm,jpm) = rho(ipm,jpm) + sx2m * sy2m * chrg
           ! rho(ipm,jp ) = rho(ipm,jp ) + sx2m * sy2  * chrg
           ! rho(ipm,jpp) = rho(ipm,jpp) + sx2m * sy2p * chrg
           ! rho(ip ,jpm) = rho(ip ,jpm) + sx2  * sy2m * chrg
           ! rho(ip ,jp ) = rho(ip ,jp ) + sx2  * sy2  * chrg
           ! rho(ip ,jpp) = rho(ip ,jpp) + sx2  * sy2p * chrg
           ! rho(ipp,jpm) = rho(ipp,jpm) + sx2p * sy2m * chrg
           ! rho(ipp,jp ) = rho(ipp,jp ) + sx2p * sy2  * chrg
           ! rho(ipp,jpp) = rho(ipp,jpp) + sx2p * sy2p * chrg
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
    subroutine phia(nx,ny,c,omega,dt,phi,phib,jx,jy,jz,Ax,Ay,Az,Axb,Ayb,Azb,&
                    Axbb,Aybb,Azbb,cfact)
!***********************************************************************
   !original subroutine
      implicit none
      real(8), dimension(0:nx,0:ny) :: phi, phib, jx, jy, jz, Ax, Ay, Az,&
      Axb, Ayb, Azb, Axbb, Aybb, Azbb
      integer :: nx, ny, i, j, k, im, ip, jm, jp
      real(rkind) :: c, omega, dt, cfact

 ! Solution of maxwell equation in the A-phi formulation by difference method
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
      
        Ax(i,j) = dt ** 2 * c ** 2 / omega ** 2 * (Axb(ip,j) + Axb(im,j) &
                + Axb(i,jp) + Axb(i,jm) - 4.0d0 * Axb(i,j)) &
                + dt ** 2 * jx(i,j) - dt * (phi(ip,j) - phib(ip,j) &
                - phi(i,j) + phib(i,j)) + 2.0d0 * Axb(i,j) - Axbb(i,j) 

        Ay(i,j) = dt ** 2 * c ** 2 / omega ** 2 * (Ayb(ip,j) + Ayb(im,j) &
                + Ayb(i,jp) + Ayb(i,jm) - 4.0d0 * Ayb(i,j)) &
                + dt ** 2 * jy(i,j) - dt * (phi(i,jp) - phib(i,jp) &
                - phi(i,j) + phib(i,j)) + 2.0d0 * Ayb(i,j) - Aybb(i,j) 

        Az(i,j) = dt ** 2 * c ** 2 / omega ** 2 * (Azb(ip,j) + Azb(im,j) &
                + Azb(i,jp) + Azb(i,jm) - 4.0d0 * Azb(i,j)) &
                + dt ** 2 * jz(i,j) + 2.0d0 * Azb(i,j) - Azbb(i,j)
      
      end do
   end do

   !do i = 0, ny
    !     do j = 0, nx
     !       Ax(i,j) = cfact * Ax(i,j)
     !       Ay(i,j) = cfact * Ay(i,j)
     !       Az(i,j) = cfact * Az(i,j)
     !    end do
     !    end do
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
   subroutine efield(nx,ny,nz,phi,ex,ey,ezg,Axb,Ayb,Azb,Axbb,Aybb,Azbb)
!***********************************************************************
     implicit none
      real(8), dimension(0:nx,0:ny) :: phi,Axb,Ayb,Azb,Axbb,Aybb,Azbb
      real(8), dimension(0:nx,0:ny,0:nz) :: ex, ey, ezg
      integer :: nx, ny, nz, i, j, im, ip, jm, jp, k, km, kp
      real(rkind) :: ez
      ez=0.d0
      do j = 0, ny
      do i = 0, nx

         im = i - 1
         ip = i + 1
         jm = j - 1
         jp = j + 1
         !km = k - 1
         !kp = k + 1

         if( i .eq. 0  ) im = nx - 1
         if( i .eq. nx ) ip = 1
         if( j .eq. 0  ) jm = ny - 1
         if( j .eq. ny ) jp = 1
         !if( k .eq. 0  ) km = nz - 1
         !if( k .eq. nz ) kp = 1

         ex(i,j,0) = 0.5d0 * ( phi(im,j ) - phi(ip,j ) ) &
                   - (Ax(i,j) - Axb(i,j)) / dt
         ey(i,j,0) = 0.5d0 * ( phi(i ,jm) - phi(i ,jp) ) &
                   - (Ay(i,j) - Ayb(i,j)) / dt 
         ezg(i,j,0)= - (Az(i,j) - Azb(i,j)) / dt 
         !ex(i,j,k) = 0
         !ey(i,j,k) = 0

      end do
      end do
    end subroutine efield
!***********************************************************************
    subroutine bfield(nx,ny,nz,bx,by,bz,bxg,byg,bzg,Axb,Ayb,Azb)
!***********************************************************************
      implicit none
      real(rkind), dimension(0:nx,0:ny,0:nz) :: bxg, byg, bzg
      real(rkind), dimension(0:nx,0:ny) :: Axb, Ayb, Azb
      integer :: nx, ny, nz, i, j, k, im, ip, jm, jp, km, kp
      real(rkind)::bx, by, bz

      do j = 0, ny
      do i = 0, nx
      !do k = 0, nz
         im = i - 1
         ip = i + 1
         jm = j - 1
         jp = j + 1
         !km = k - 1
         !kp = k + 1

         if( i .eq. 0  ) im = nx - 1
         if( i .eq. nx ) ip = 1
         if( j .eq. 0  ) jm = ny - 1
         if( j .eq. ny ) jp = 1
         !if( k .eq. 0  ) km = nz - 1
         !if( k .eq. nz ) kp = 1
         

         bxg(i,j,0) = 0.5d0 * (Azb(ip,j) - Azb(im,j))
         byg(i,j,0) = - 0.5d0 * (Azb(i,jp) - Azb(i,jm))
         bzg(i,j,0) = 0.5d0 * (Ayb(ip,j) - Ayb(im,j) - (Axb(i,jp) - Axb(i,jm)))

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

         rho(i,j) = phi(im,j) + phi(ip,j) + phi(i,jm) + phi(i,jp) &
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
    subroutine pote(nx,ny,nz,phi,ex,ey,ez,Ax,Ay,Az,apot,cfacti)
!***********************************************************************
      implicit none
      real(8), dimension(0:nx,0:ny,0:nz) :: ex, ey, ez
      real(8), dimension(0:nx,0:ny) :: phi, Ax, Ay, Az
      real(8) :: apot, cfacti
      integer(4) :: nx, ny, nz, i, j, im, ip, jm, jp

      apot = 0.d0

      do j = 0, ny-1
      do i = 0, nx-1

         im = i - 1
         ip = i + 1
         jm = j - 1
         jp = j + 1

         if( i .eq. 0  ) im = nx - 1
         if( i .eq. nx ) ip = 1
         if( j .eq. 0  ) jm = ny - 1
         if( j .eq. ny ) jp = 1     
         apot = apot + ex(i,j,0)**2 + ey(i,j,0)**2 + ez(i,j,0)**2
              !+ 0.25d0 * (phi(ip,j) - phi(im,j)) ** 2 &
              !+ 0.25d0 * (phi(i,jp) - phi(i,jm)) ** 2 &
              !+ 0.25d0 / c**2 * ((Ayb(ip,j) - Ayb(im,j) - Axb(i,jp) &
              !+ Axb(i,jm))** 2 + (Azb(ip,j) - Azb(im,j)) ** 2 &
              !+ (Azb(i,jp) - Azb(i,jm)) ** 2)
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
 
!.....................................................
      endif
!.....................................................
 
    end subroutine sumdim1

!***********************************************************************
    subroutine current(np,nx,ny,nz,x,y,z,vx,vy,vz,jx,jy,jz,dt,chrg,cfact)
!***********************************************************************
    implicit none

    real(8), dimension(np) :: x, y, z, vx, vy, vz 
    real(8), dimension(0:nx,0:ny) :: jx, jy, jz
    real(8) :: chrg, dt, dx, dy, dz, dx1, dy1, dz1, deltax,deltay,deltaz,&
               cfact, sx1p, sy1p, sx1m, sy1m, sx2, sy2, sx2p, sy2p, sx2m, sy2m
    integer :: np, nx, ny, nz, i, ip, jp, kp, ix, iy, ipp, ipm, jpp, jpm 

     do i = 1, np
        
        ip = x(i)
        jp = y(i)
        kp = z(i)
        dx = x(i) - dble(ip)
        dy = y(i) - dble(jp)
        dz = z(i) - dble(kp)
        dx1 = 1.0d0 - dx
        dy1 = 1.0d0 - dy 
        dz1 = 1.0d0 - dz
        !deltax = x(i) - xb(i)
        !deltay = y(i) - yb(i)
        !deltaz = z(i) - zb(i)
         if (dx .le. 0.5d0) then 
            sx1p = 0.5d0 + dx
            sx1m = 0.5d0 - dx
         else
            sx1p = dx - 0.5d0
            sx1m = 1.5d0 - dx
         end if
         if (dy .le. 0.5d0) then 
            sy1p = 0.5d0 + dy
            sy1m = 0.5d0 - dy
         else
            sy1p = dy - 0.5d0
            sy1m = 1.5d0 - dy
         end if
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

         if( ip .eq. 0  ) ipm = nx - 1
         if( ip .eq. nx ) ipp = 1
         if( jp .eq. 0  ) jpm = ny - 1
         if( jp .eq. ny ) jpp = 1

         jx(ip  ,jp  ) = jx(ip  ,jp  ) + dx1 * dy1 * chrg
         jx(ip+1,jp  ) = jx(ip+1,jp  ) + dx  * dy1 * chrg
         jx(ip  ,jp+1) = jx(ip  ,jp+1) + dx1 * dy  * chrg
         jx(ip+1,jp+1) = jx(ip+1,jp+1) + dx  * dy  * chrg

         jy(ip  ,jp  ) = jy(ip  ,jp  ) + dx1 * dy1 * chrg
         jy(ip+1,jp  ) = jy(ip+1,jp  ) + dx  * dy1 * chrg
         jy(ip  ,jp+1) = jy(ip  ,jp+1) + dx1 * dy  * chrg
         jy(ip+1,jp+1) = jy(ip+1,jp+1) + dx  * dy  * chrg

         jz(ip  ,jp  ) = jz(ip  ,jp  ) + dx1 * dy1 * chrg
         jz(ip+1,jp  ) = jz(ip+1,jp  ) + dx  * dy1 * chrg
         jz(ip  ,jp+1) = jz(ip  ,jp+1) + dx1 * dy  * chrg
         jz(ip+1,jp+1) = jz(ip+1,jp+1) + dx  * dy  * chrg
         
          ! jx(ip ,jp ) = jx(ip ,jp ) + chrg * vx(i) * sy2  * sx1p
          ! jx(ip ,jpp) = jx(ip ,jpp) + chrg * vx(i) * sy2p * sx1p
          ! jx(ip ,jpm) = jx(ip ,jpm) + chrg * vx(i) * sy2m * sx1p
          ! jx(ipm,jp ) = jx(ipm,jp ) + chrg * vx(i) * sy2  * sx1m
          ! jx(ipm,jpp) = jx(ipm,jpp) + chrg * vx(i) * sy2p * sx1m
          ! jx(ipm,jpm) = jx(ipm,jpm) + chrg * vx(i) * sy2m * sx1m
           
          ! jy(ip ,jp ) = jy(ip ,jp ) + chrg * vy(i) * sx2  * sy1p
          ! jy(ipp,jp ) = jy(ipm,jp ) + chrg * vy(i) * sx2p * sy1p
          ! jy(ipm,jp ) = jy(ipp,jp ) + chrg * vy(i) * sx2m * sy1p
          ! jy(ip ,jpm) = jy(ip ,jpm) + chrg * vy(i) * sx2  * sy1m
          ! jy(ipp,jpm) = jy(ipm,jpm) + chrg * vy(i) * sx2p * sy1m
          ! jy(ipm,jpm) = jy(ipp,jpm) + chrg * vy(i) * sx2m * sy1m
      
          ! jz(ipm,jpp) = jz(ipm,jpp) + chrg * vz(i) * sx2m * sy2p
          ! jz(ipm,jp ) = jz(ipm,jp ) + chrg * vz(i) * sx2m * sy2
          ! jz(ipm,jpm) = jz(ipm,jpm) + chrg * vz(i) * sx2m * sy2m
          ! jz(ip ,jpp) = jz(ip ,jpp) + chrg * vz(i) * sx2  * sy2p
          ! jz(ip ,jp ) = jz(ip ,jp ) + chrg * vz(i) * sx2  * sy2
          ! jz(ip ,jpm) = jz(ip ,jpm) + chrg * vz(i) * sx2  * sy2m
          ! jz(ipp,jpp) = jz(ipp,jpp) + chrg * vz(i) * sx2p * sy2p
          ! jz(ipp,jp ) = jz(ipp,jp ) + chrg * vz(i) * sx2p * sy2
          ! jz(ipp,jpm) = jz(ipp,jpm) + chrg * vz(i) * sx2p * sy2m
       
     end do

   end subroutine current
END Module picexec
