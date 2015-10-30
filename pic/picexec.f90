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

         !.......... calculate ex and ey and ez
         !call efield(nx,ny,nz,phi,ex,ey,ez,ezg)

         !.......... calculate bxg and byg and bzg
         ! call bfield(nx,ny,nz,bx,by,bz,bxg,byg,bzg)
        
         !..... diagnostics to check energy conservation 
         !.....            before pushing 
         if( mod(iloop,nhmod) .eq. 0 ) then
            iene = iene + 1
            call kine(np,vxe,vye,vze,akine1,me)
            call kine(np,vxi,vyi,vzi,akini1,mi)
            call pote(nx,ny,nz,ex,ey,ez,apot,cfacti)
            call sumdim1(nodes,myid,akine1,wkword)
            call sumdim1(nodes,myid,akini1,wkword)
         endif
         
         jx(:,:)=0.d0
         jy(:,:)=0.d0
         jz(:,:)=0.d0
         
         !.......... calculate current by electrons
         call current(np,nx,ny,nz,xe,ye,ze,xeb,yeb,zeb,jx,jy,jz,dt,chrge) 
          !.......... calculate current by ions
         call current(np,nx,ny,nz,xi,yi,zi,xib,yib,zib,jx,jy,jz,dt,chrgi)
         !jx(:,:)=0.d0
         !jy(:,:)=0.d0
         !jz(:,:)=0.d0
         !.......... calculate vector potential
         call phia(nx,ny,c,omega,dt,phi,phib,jx,jy,jz,Ax,Ay,Az,Axb,Ayb,Azb,&
              Axbb,Aybb,Azbb)
         ex(:,:,:)=0.d0
         ey(:,:,:)=0.d0
         ez(:,:,:)=0.d0
         !----- push electrons
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
            apot=0.d0
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
       subroutine push(np,nx,ny,nz,x,y,z,vx,vy,vz,ex,ey,ez,bxg,byg,bzg,dt,&
            ctom,xb,yb,zb,phi,phib,Axb,Ayb,Azb,Axbb,Aybb,Azbb,vpara,vperp)
!***********************************************************************
      implicit none
      real(8), dimension(np) :: x, y, z, xb, yb, zb, vx, vy, vz, vpara, vperp
      real(8), dimension(0:nx,0:ny) :: phi,phib,Axb, Ayb, Azb, Axbb, Aybb, Azbb
      real(8), dimension(0:nx,0:ny,0:nz) :: ex, ey, ez, bxg, byg, bzg
      real(8) :: ctom, dx, dy, dz, dx1, dy1, dz1, dt, exx, eyy, ezz, bxx,&
                 byy, bzz, vxn, vyn, vzn, vxzero, vyzero, vzzero, vxp, vyp,&
                 vzp, sx1p, sx1m, sy1p, sy1m, sx2, sy2, sx2p, sx2m, sy2m, sy2p
      real(8) :: btot, vtot
      integer :: np, nx, ny, nz, i, j, ip, jp, kp
     
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
         ! electric field
         if(jp .ge. 1 .and. ip .ge. 1) then
            exx = - phi(ip+1,jp+1) * sy2p * sx1p &
                - phi(ip+1,jp  ) * sy2  * sx1p &
                - phi(ip+1,jp-1) * sy2m * sx1p &
                + phi(ip  ,jp+1) * sy2p * (sx1p - sx1m) &
                + phi(ip  ,jp  ) * sy2  * (sx1p - sx1m) &
                + phi(ip  ,jp-1) * sy2m * (sx1p - sx1m) &
                + phi(ip-1,jp+1) * sy2p * sx1m &
                + phi(ip-1,jp  ) * sy2  * sx1m &
                + phi(ip-1,jp-1) * sy2m * sx1m &
                - (Axb(ip  ,jp+1) - Axbb(ip  ,jp+1)) / dt * sx1p * sy2p &
                - (Axb(ip  ,jp  ) - Axbb(ip  ,jp  )) / dt * sx1p * sy2  &
                - (Axb(ip  ,jp-1) - Axbb(ip  ,jp-1)) / dt * sx1p * sy2m &
                - (Axb(ip-1,jp+1) - Axbb(ip-1,jp+1)) / dt * sx1m * sy2p &
                - (Axb(ip-1,jp  ) - Axbb(ip-1,jp  )) / dt * sx1m * sy2  &
                - (Axb(ip-1,jp-1) - Axbb(ip-1,jp-1)) / dt * sx1m * sy2m 
         else if(ip .ge. 1 .and. jp .eq. 0) then
            exx = - phi(ip+1,jp+1) * sy2p * sx1p &
                - phi(ip+1,jp  ) * sy2  * sx1p &
                - phi(ip+1,ny-1) * sy2m * sx1p &
                + phi(ip  ,jp+1) * sy2p * (sx1p - sx1m) &
                + phi(ip  ,jp  ) * sy2  * (sx1p - sx1m) &
                + phi(ip  ,ny-1) * sy2m * (sx1p - sx1m) &
                + phi(ip-1,jp+1) * sy2p * sx1m &
                + phi(ip-1,jp  ) * sy2  * sx1m &
                + phi(ip-1,ny-1) * sy2m * sx1m &
                - (Axb(ip  ,jp+1) - Axbb(ip  ,jp+1)) / dt * sx1p * sy2p &
                - (Axb(ip  ,jp  ) - Axbb(ip  ,jp  )) / dt * sx1p * sy2  &
                - (Axb(ip  ,ny-1) - Axbb(ip  ,ny-1)) / dt * sx1p * sy2m &
                - (Axb(ip-1,jp+1) - Axbb(ip-1,jp+1)) / dt * sx1m * sy2p &
                - (Axb(ip-1,jp  ) - Axbb(ip-1,jp  )) / dt * sx1m * sy2  &
                - (Axb(ip-1,ny-1) - Axbb(ip-1,ny-1)) / dt * sx1m * sy2m 
         else if(ip .eq. 0 .and. jp .ge. 1) then
            exx = - phi(ip+1,jp+1) * sy2p * sx1p &
                - phi(ip+1,jp  ) * sy2  * sx1p &
                - phi(ip+1,jp-1) * sy2m * sx1p &
                + phi(ip  ,jp+1) * sy2p * (sx1p - sx1m) &
                + phi(ip  ,jp  ) * sy2  * (sx1p - sx1m) &
                + phi(ip  ,jp-1) * sy2m * (sx1p - sx1m) &
                + phi(nx-1,jp+1) * sy2p * sx1m &
                + phi(nx-1,jp  ) * sy2  * sx1m &
                + phi(nx-1,jp-1) * sy2m * sx1m &
                - (Axb(ip  ,jp+1) - Axbb(ip  ,jp+1)) / dt * sx1p * sy2p &
                - (Axb(ip  ,jp  ) - Axbb(ip  ,jp  )) / dt * sx1p * sy2  &
                - (Axb(ip  ,jp-1) - Axbb(ip  ,jp-1)) / dt * sx1p * sy2m &
                - (Axb(nx-1,jp+1) - Axbb(nx-1,jp+1)) / dt * sx1m * sy2p &
                - (Axb(nx-1,jp  ) - Axbb(nx-1,jp  )) / dt * sx1m * sy2  &
                - (Axb(nx-1,jp-1) - Axbb(nx-1,jp-1)) / dt * sx1m * sy2m 
         else
            exx = - phi(ip+1,jp+1) * sy2p * sx1p &
                - phi(ip+1,jp  ) * sy2  * sx1p &
                - phi(ip+1,ny-1) * sy2m * sx1p &
                + phi(ip  ,jp+1) * sy2p * (sx1p - sx1m) &
                + phi(ip  ,jp  ) * sy2  * (sx1p - sx1m) &
                + phi(ip  ,ny-1) * sy2m * (sx1p - sx1m) &
                + phi(nx-1,jp+1) * sy2p * sx1m &
                + phi(nx-1,jp  ) * sy2  * sx1m &
                + phi(nx-1,ny-1) * sy2m * sx1m &
                - (Axb(ip  ,jp+1) - Axbb(ip  ,jp+1)) / dt * sx1p * sy2p &
                - (Axb(ip  ,jp  ) - Axbb(ip  ,jp  )) / dt * sx1p * sy2  &
                - (Axb(ip  ,ny-1) - Axbb(ip  ,ny-1)) / dt * sx1p * sy2m &
                - (Axb(nx-1,jp+1) - Axbb(nx-1,jp+1)) / dt * sx1m * sy2p &
                - (Axb(nx-1,jp  ) - Axbb(nx-1,jp  )) / dt * sx1m * sy2  &
                - (Axb(nx-1,ny-1) - Axbb(nx-1,ny-1)) / dt * sx1m * sy2m         
            
         endif
      
         if(ip .ge. 1 .and. jp .ge. 1) then
            eyy = - phi(ip+1,jp+1) * sx2p * sy1p &
                - phi(ip  ,jp+1) * sx2  * sy1p &
                - phi(ip-1,jp+1) * sx2m * sy1p &
                + phi(ip+1,jp  ) * sx2p * (sy1p - sy1m) &
                + phi(ip  ,jp  ) * sx2  * (sy1p - sy1m) &
                + phi(ip-1,jp  ) * sx2m * (sy1p - sy1m) &
                + phi(ip+1,jp-1) * sx2p * sy1m &
                + phi(ip  ,jp-1) * sx2  * sy1m &
                + phi(ip-1,jp-1) * sx2m * sy1m &
                - (Ayb(ip+1,jp  ) - Aybb(ip+1,jp  )) / dt * sx2p * sy1p &
                - (Ayb(ip  ,jp  ) - Aybb(ip  ,jp  )) / dt * sx2  * sy1p &
                - (Ayb(ip-1,jp  ) - Aybb(ip-1,jp  )) / dt * sx2m * sy1p &
                - (Ayb(ip+1,jp-1) - Aybb(ip+1,jp-1)) / dt * sx2p * sy1m &
                - (Ayb(ip  ,jp-1) - Aybb(ip  ,jp-1)) / dt * sx2  * sy1m &
                - (Ayb(ip-1,jp-1) - Aybb(ip-1,jp-1)) / dt * sx2m * sy1m
                
         else if (ip .eq. 0 .and. jp .ge. 1 ) then
            eyy = - phi(ip+1,jp+1) * sx2p * sy1p &
                - phi(ip  ,jp+1) * sx2  * sy1p &
                - phi(nx-1,jp+1) * sx2m * sy1p &
                + phi(ip+1,jp  ) * sx2p * (sy1p - sy1m) &
                + phi(ip  ,jp  ) * sx2  * (sy1p - sy1m) &
                + phi(nx-1,jp  ) * sx2m * (sy1p - sy1m) &
                + phi(ip+1,jp-1) * sx2p * sy1m &
                + phi(ip  ,jp-1) * sx2  * sy1m &
                + phi(nx-1,jp-1) * sx2m * sy1m &
                - (Ayb(ip+1,jp  ) - Aybb(ip+1,jp  )) / dt * sx2p * sy1p &
                - (Ayb(ip  ,jp  ) - Aybb(ip  ,jp  )) / dt * sx2  * sy1p &
                - (Ayb(nx-1,jp  ) - Aybb(nx-1,jp  )) / dt * sx2m * sy1p &
                - (Ayb(ip+1,jp-1) - Aybb(ip+1,jp-1)) / dt * sx2p * sy1m &
                - (Ayb(ip  ,jp-1) - Aybb(ip  ,jp-1)) / dt * sx2  * sy1m &
                - (Ayb(nx-1,jp-1) - Aybb(nx-1,jp-1)) / dt * sx2m * sy1m
         else if (ip .ge. 1 .and. jp .eq. 0 ) then
            eyy = - phi(ip+1,jp+1) * sx2p * sy1p &
                - phi(ip  ,jp+1) * sx2  * sy1p &
                - phi(ip-1,jp+1) * sx2m * sy1p &
                + phi(ip+1,jp  ) * sx2p * (sy1p - sy1m) &
                + phi(ip  ,jp  ) * sx2  * (sy1p - sy1m) &
                + phi(ip-1,jp  ) * sx2m * (sy1p - sy1m) &
                + phi(ip+1,ny-1) * sx2p * sy1m &
                + phi(ip  ,ny-1) * sx2  * sy1m &
                + phi(ip-1,ny-1) * sx2m * sy1m &
                - (Ayb(ip+1,jp  ) - Aybb(ip+1,jp  )) / dt * sx2p * sy1p &
                - (Ayb(ip  ,jp  ) - Aybb(ip  ,jp  )) / dt * sx2  * sy1p &
                - (Ayb(ip-1,jp  ) - Aybb(ip-1,jp  )) / dt * sx2m * sy1p &
                - (Ayb(ip+1,ny-1) - Aybb(ip+1,ny-1)) / dt * sx2p * sy1m &
                - (Ayb(ip  ,ny-1) - Aybb(ip  ,ny-1)) / dt * sx2  * sy1m &
                - (Ayb(ip-1,ny-1) - Aybb(ip-1,ny-1)) / dt * sx2m * sy1m
         else
            eyy = - phi(ip+1,jp+1) * sx2p * sy1p &
                - phi(ip  ,jp+1) * sx2  * sy1p &
                - phi(nx-1,jp+1) * sx2m * sy1p &
                + phi(ip+1,jp  ) * sx2p * (sy1p - sy1m) &
                + phi(ip  ,jp  ) * sx2  * (sy1p - sy1m) &
                + phi(nx-1,jp  ) * sx2m * (sy1p - sy1m) &
                + phi(ip+1,ny-1) * sx2p * sx1m &
                + phi(ip  ,ny-1) * sx2  * sy1m &
                + phi(nx-1,ny-1) * sx2m * sy1m &
                - (Ayb(ip+1,jp  ) - Aybb(ip+1,jp  )) / dt * sx2p * sy1p &
                - (Ayb(ip  ,jp  ) - Aybb(ip  ,jp  )) / dt * sx2  * sy1p &
                - (Ayb(nx-1,jp  ) - Aybb(nx-1,jp  )) / dt * sx2m * sy1p &
                - (Ayb(ip+1,ny-1) - Aybb(ip+1,ny-1)) / dt * sx2p * sy1m &
                - (Ayb(ip  ,ny-1) - Aybb(ip  ,ny-1)) / dt * sx2  * sy1m &
                - (Ayb(nx-1,ny-1) - Aybb(nx-1,ny-1)) / dt * sx2m * sy1m
            
         endif
         
         if (ip .ge. 1 .and. jp .ge. 1) then
            ezz = - (Azb(ip+1,jp+1) - Azbb(ip+1,jp+1)) * sx2p * sy2p / dt &
                - (Azb(ip+1,jp  ) - Azbb(ip+1,jp  )) * sx2p * sy2  / dt &
                - (Azb(ip+1,jp-1) - Azbb(ip+1,jp-1)) * sx2p * sy2m / dt &
                - (Azb(ip  ,jp+1) - Azbb(ip  ,jp+1)) * sx2  * sy2p / dt &
                - (Azb(ip  ,jp  ) - Azbb(ip  ,jp  )) * sx2  * sy2  / dt &
                - (Azb(ip  ,jp-1) - Azbb(ip  ,jp-1)) * sx2  * sy2m / dt &
                - (Azb(ip-1,jp+1) - Azbb(ip-1,jp+1)) * sx2m * sy2p / dt &
                - (Azb(ip-1,jp  ) - Azbb(ip-1,jp  )) * sx2m * sy2  / dt &
                - (Azb(ip-1,jp-1) - Azbb(ip-1,jp-1)) * sx2m * sy2m / dt         
         elseif(ip .ge. 1 .and. jp .eq. 0) then
            ezz = - (Azb(ip+1,jp+1) - Azbb(ip+1,jp+1)) * sx2p * sy2p / dt &
                - (Azb(ip+1,jp  ) - Azbb(ip+1,jp  )) * sx2p * sy2  / dt &
                - (Azb(ip+1,ny-1) - Azbb(ip+1,ny-1)) * sx2p * sy2m / dt &
                - (Azb(ip  ,jp+1) - Azbb(ip  ,jp+1)) * sx2  * sy2p / dt &
                - (Azb(ip  ,jp  ) - Azbb(ip  ,jp  )) * sx2  * sy2  / dt &
                - (Azb(ip  ,ny-1) - Azbb(ip  ,ny-1)) * sx2  * sy2m / dt &
                - (Azb(ip-1,jp+1) - Azbb(ip-1,jp+1)) * sx2m * sy2p / dt &
                - (Azb(ip-1,jp  ) - Azbb(ip-1,jp  )) * sx2m * sy2  / dt &
                - (Azb(ip-1,ny-1) - Azbb(ip-1,ny-1)) * sx2m * sy2m / dt
         elseif(ip .eq. 0 .and. jp .ge. 1) then
            ezz = - (Azb(ip+1,jp+1) - Azbb(ip+1,jp+1)) * sx2p * sy2p / dt &
                - (Azb(ip+1,jp  ) - Azbb(ip+1,jp  )) * sx2p * sy2  / dt &
                - (Azb(ip+1,jp-1) - Azbb(ip+1,jp-1)) * sx2p * sy2m / dt &
                - (Azb(ip  ,jp+1) - Azbb(ip  ,jp+1)) * sx2  * sy2p / dt &
                - (Azb(ip  ,jp  ) - Azbb(ip  ,jp  )) * sx2  * sy2  / dt &
                - (Azb(ip  ,jp-1) - Azbb(ip  ,jp-1)) * sx2  * sy2m / dt &
                - (Azb(nx-1,jp+1) - Azbb(nx-1,jp+1)) * sx2m * sy2p / dt &
                - (Azb(nx-1,jp  ) - Azbb(nx-1,jp  )) * sx2m * sy2  / dt &
                - (Azb(nx-1,jp-1) - Azbb(nx-1,jp-1)) * sx2m * sy2m / dt
         else
            ezz = - (Azb(ip+1,jp+1) - Azbb(ip+1,jp+1)) * sx2p * sy2p / dt &
                - (Azb(ip+1,jp  ) - Azbb(ip+1,jp  )) * sx2p * sy2  / dt &
                - (Azb(ip+1,ny-1) - Azbb(ip+1,ny-1)) * sx2p * sy2m / dt &
                - (Azb(ip  ,jp+1) - Azbb(ip  ,jp+1)) * sx2  * sy2p / dt &
                - (Azb(ip  ,jp  ) - Azbb(ip  ,jp  )) * sx2  * sy2  / dt &
                - (Azb(ip  ,ny-1) - Azbb(ip  ,ny-1)) * sx2  * sy2m / dt &
                - (Azb(nx-1,jp+1) - Azbb(nx-1,jp+1)) * sx2m * sy2p / dt &
                - (Azb(nx-1,jp  ) - Azbb(nx-1,jp  )) * sx2m * sy2  / dt &
                - (Azb(nx-1,ny-1) - Azbb(nx-1,ny-1)) * sx2m * sy2m / dt
         endif
         ex(ip,jp,kp) = ex(ip,jp,kp) + exx
         ey(ip,jp,kp) = ey(ip,jp,kp) + eyy
         ez(ip,jp,kp) = ez(ip,jp,kp) + ezz
      
     ! exx = ex(ip  ,jp  ,kp  )*dx1*dy1*dz1 + ex(ip+1,jp  ,kp  )*dx*dy1*dz1 &
     !     + ex(ip  ,jp+1,kp  )*dx1*dy*dz1  + ex(ip  ,jp  ,kp+1)*dx1*dy1*dz &
     !     + ex(ip+1,jp+1,kp  )*dx*dy*dz1   + ex(ip  ,jp+1,kp+1)*dx1*dy*dz &
     !     + ex(ip+1,jp  ,kp+1)*dx*dy1*dz   + ex(ip+1,jp+1,kp+1)*dx*dy*dz

     ! eyy = ey(ip  ,jp  ,kp  )*dx1*dy1*dz1 + ey(ip+1,jp  ,kp  )*dx*dy1*dz1 &
     !     + ey(ip  ,jp+1,kp  )*dx1*dy*dz1  + ey(ip  ,jp  ,kp+1)*dx1*dy1*dz &
     !     + ey(ip+1,jp+1,kp  )*dx*dy*dz1   + ey(ip  ,jp+1,kp+1)*dx1*dy*dz &
     !     + ey(ip+1,jp  ,kp+1)*dx*dy1*dz   + ey(ip+1,jp+1,kp+1)*dx*dy*dz

     ! ezz = ezg(ip  ,jp  ,kp  )*dx1*dy1*dz1 + ezg(ip+1,jp  ,kp  )*dx*dy1*dz1 &
     !     + ezg(ip  ,jp+1,kp  )*dx1*dy*dz1  + ezg(ip  ,jp  ,kp+1)*dx1*dy1*dz &
     !     + ezg(ip+1,jp+1,kp  )*dx*dy*dz1   + ezg(ip  ,jp+1,kp+1)*dx1*dy*dz &
     !     + ezg(ip+1,jp  ,kp+1)*dx*dy1*dz   + ezg(ip+1,jp+1,kp+1)*dx*dy*dz

! magnetic field
         if(ip .ge. 1 .and. jp .ge. 1) then
            bxx = (Azb(ip+1,jp+1) + Azbb(ip+1,jp+1))/2 * sx2p * sy1p &
                + (Azb(ip  ,jp+1) + Azbb(ip  ,jp+1))/2 * sx2  * sy1p &
                + (Azb(ip-1,jp+1) + Azbb(ip-1,jp+1))/2 * sx2m * sy1p &
                - (Azb(ip+1,jp  ) + Azbb(ip+1,jp  ))/2 * sx2p * (sy1p-sy1m)&
                - (Azb(ip  ,jp  ) + Azbb(ip  ,jp  ))/2 * sx2  * (sy1p-sy1m)&
                - (Azb(ip-1,jp  ) + Azbb(ip-1,jp  ))/2 * sx2m * (sy1p-sy1m)&
                - (Azb(ip+1,jp-1) + Azbb(ip+1,jp-1))/2 * sx2p * sy1m &
                - (Azb(ip  ,jp-1) + Azbb(ip  ,jp-1))/2 * sx2  * sy1m &
                - (Azb(ip-1,jp-1) + Azbb(ip-1,jp-1))/2 * sx2m * sy1m 

         else if(ip .eq. 0 .and. jp .ge. 1) then
            bxx = (Azb(ip+1,jp+1) + Azbb(ip+1,jp+1))/2 * sx2p * sy1p &
                + (Azb(ip  ,jp+1) + Azbb(ip  ,jp+1))/2 * sx2  * sy1p &
                + (Azb(nx-1,jp+1) + Azbb(nx-1,jp+1))/2 * sx2m * sy1p &
                - (Azb(ip+1,jp  ) + Azbb(ip+1,jp  ))/2 * sx2p * (sy1p-sy1m)&
                - (Azb(ip  ,jp  ) + Azbb(ip  ,jp  ))/2 * sx2  * (sy1p-sy1m)&
                - (Azb(nx-1,jp  ) + Azbb(nx-1,jp  ))/2 * sx2m * (sy1p-sy1m)&
                - (Azb(ip+1,jp-1) + Azbb(ip+1,jp-1))/2 * sx2p * sy1m &
                - (Azb(ip  ,jp-1) + Azbb(ip  ,jp-1))/2 * sx2  * sy1m &
                - (Azb(nx-1,jp-1) + Azbb(nx-1,jp-1))/2 * sx2m * sy1m 

         else if(ip .ge. 1 .and. jp .eq. 0) then
            bxx = (Azb(ip+1,jp+1) + Azbb(ip+1,jp+1))/2 * sx2p * sy1p &
                + (Azb(ip  ,jp+1) + Azbb(ip  ,jp+1))/2 * sx2  * sy1p &
                + (Azb(ip-1,jp+1) + Azbb(ip-1,jp+1))/2 * sx2m * sy1p &
                - (Azb(ip+1,jp  ) + Azbb(ip+1,jp  ))/2 * sx2p * (sy1p-sy1m)&
                - (Azb(ip  ,jp  ) + Azbb(ip  ,jp  ))/2 * sx2  * (sy1p-sy1m)&
                - (Azb(ip-1,jp  ) + Azbb(ip-1,jp  ))/2 * sx2m * (sy1p-sy1m)&
                - (Azb(ip+1,ny-1) + Azbb(ip+1,ny-1))/2 * sx2p * sy1m &
                - (Azb(ip  ,ny-1) + Azbb(ip  ,ny-1))/2 * sx2  * sy1m &
                - (Azb(ip-1,ny-1) + Azbb(ip-1,ny-1))/2 * sx2m * sy1m 

         else
            bxx = (Azb(ip+1,jp+1) + Azbb(ip+1,jp+1))/2 * sx2p * sy1p &
                + (Azb(ip  ,jp+1) + Azbb(ip  ,jp+1))/2 * sx2  * sy1p &
                + (Azb(nx-1,jp+1) + Azbb(nx-1,jp+1))/2 * sx2m * sy1p &
                - (Azb(ip+1,jp  ) + Azbb(ip+1,jp  ))/2 * sx2p * (sy1p-sy1m)&
                - (Azb(ip  ,jp  ) + Azbb(ip  ,jp  ))/2 * sx2  * (sy1p-sy1m)&
                - (Azb(nx-1,jp  ) + Azbb(nx-1,jp  ))/2 * sx2m * (sy1p-sy1m)&
                - (Azb(ip+1,ny-1) + Azbb(ip+1,ny-1))/2 * sx2p * sy1m &
                - (Azb(ip  ,ny-1) + Azbb(ip  ,ny-1))/2 * sx2  * sy1m &
                - (Azb(nx-1,ny-1) + Azbb(nx-1,ny-1))/2 * sx2m * sy1m
         endif
         
         if(ip .ge. 1 .and. jp .ge. 1) then
            byy = - (Azb(ip+1,jp+1) + Azbb(ip+1,jp+1))/2 * sy2p * sx1p &
                - (Azb(ip+1,jp  ) + Azbb(ip+1,jp  ))/2 * sy2  * sx1p &
                - (Azb(ip+1,jp-1) + Azbb(ip+1,jp-1))/2 * sy2m * sx1p &
                + (Azb(ip  ,jp+1) + Azbb(ip  ,jp+1))/2 * sy2p * (sx1p-sx1m) &
                + (Azb(ip  ,jp  ) + Azbb(ip  ,jp  ))/2 * sy2  * (sx1p-sx1m) &
                + (Azb(ip  ,jp-1) + Azbb(ip  ,jp-1))/2 * sy2m * (sx1p-sx1m) &
                + (Azb(ip-1,jp+1) + Azbb(ip-1,jp+1))/2 * sy2p * sx1m &
                + (Azb(ip-1,jp  ) + Azbb(ip-1,jp  ))/2 * sy2  * sx1m &
                + (Azb(ip-1,jp-1) + Azbb(ip-1,jp-1))/2 * sy2m * sx1m
         else if(ip .eq. 0 .and. jp .ge. 1) then
            byy = - (Azb(ip+1,jp+1) + Azbb(ip+1,jp+1))/2 * sy2p * sx1p &
                - (Azb(ip+1,jp  ) + Azbb(ip+1,jp  ))/2 * sy2  * sx1p &
                - (Azb(ip+1,jp-1) + Azbb(ip+1,jp-1))/2 * sy2m * sx1p &
                + (Azb(ip  ,jp+1) + Azbb(ip  ,jp+1))/2 * sy2p * (sx1p-sx1m) &
                + (Azb(ip  ,jp  ) + Azbb(ip  ,jp  ))/2 * sy2  * (sx1p-sx1m) &
                + (Azb(ip  ,jp-1) + Azbb(ip  ,jp-1))/2 * sy2m * (sx1p-sx1m) &
                + (Azb(nx-1,jp+1) + Azbb(nx-1,jp+1))/2 * sy2p * sx1m &
                + (Azb(nx-1,jp  ) + Azbb(nx-1,jp  ))/2 * sy2  * sx1m &
                + (Azb(nx-1,jp-1) + Azbb(nx-1,jp-1))/2 * sy2m * sx1m
         else if(ip .ge. 1 .and. jp .eq. 0) then
            byy = - (Azb(ip+1,jp+1) + Azbb(ip+1,jp+1))/2 * sy2p * sx1p &
                - (Azb(ip+1,jp  ) + Azbb(ip+1,jp  ))/2 * sy2  * sx1p &
                - (Azb(ip+1,ny-1) + Azbb(ip+1,ny-1))/2 * sy2m * sx1p &
                + (Azb(ip  ,jp+1) + Azbb(ip  ,jp+1))/2 * sy2p * (sx1p-sx1m) &
                + (Azb(ip  ,jp  ) + Azbb(ip  ,jp  ))/2 * sy2  * (sx1p-sx1m) &
                + (Azb(ip  ,ny-1) + Azbb(ip  ,ny-1))/2 * sy2m * (sx1p-sx1m) &
                + (Azb(ip-1,jp+1) + Azbb(ip-1,jp+1))/2 * sy2p * sx1m &
                + (Azb(ip-1,jp  ) + Azbb(ip-1,jp  ))/2 * sy2  * sx1m &
                + (Azb(ip-1,ny-1) + Azbb(ip-1,ny-1))/2 * sy2m * sx1m
         else
            byy = - (Azb(ip+1,jp+1) + Azbb(ip+1,jp+1))/2 * sy2p * sx1p &
                - (Azb(ip+1,jp  ) + Azbb(ip+1,jp  ))/2 * sy2  * sx1p &
                - (Azb(ip+1,ny-1) + Azbb(ip+1,ny-1))/2 * sy2m * sx1p &
                + (Azb(ip  ,jp+1) + Azbb(ip  ,jp+1))/2 * sy2p * (sx1p-sx1m) &
                + (Azb(ip  ,jp  ) + Azbb(ip  ,jp  ))/2 * sy2  * (sx1p-sx1m) &
                + (Azb(ip  ,ny-1) + Azbb(ip  ,ny-1))/2 * sy2m * (sx1p-sx1m) &
                + (Azb(nx-1,jp+1) + Azbb(nx-1,jp+1))/2 * sy2p * sx1m &
                + (Azb(nx-1,jp  ) + Azbb(nx-1,jp  ))/2 * sy2  * sx1m &
                + (Azb(nx-1,ny-1) + Azbb(nx-1,ny-1))/2 * sy2m * sx1m 
         endif

         if(ip .ge. 1 .and. jp .ge. 1) then
            bzz = (Ayb(ip+1,jp  ) + Aybb(ip+1,jp  )) / 2 * sy1p * sx1p &
                - (Ayb(ip  ,jp  ) + Aybb(ip  ,jp  )) / 2 * sy1p * (sx1p-sx1m)&
                - (Ayb(ip-1,jp  ) + Aybb(ip-1,jp  )) / 2 * sy1p * sy1m &
                + (Ayb(ip+1,jp-1) + Aybb(ip+1,jp-1)) / 2 * sy1m * sx1p &
                + (Ayb(ip  ,jp-1) + Aybb(ip  ,jp-1)) / 2 * sy1m * (sx1p-sx1m) &
                + (Ayb(ip-1,jp-1) + Aybb(ip-1,jp-1)) / 2 * sy1m * sx1m &
                - (Axb(ip  ,jp+1) + Axbb(ip  ,jp+1)) / 2 * sx1p * sy1p &
                + (Axb(ip  ,jp  ) + Axbb(ip  ,jp  )) / 2 * sx1p * (sy1p-sy1m) &
                + (Axb(ip  ,jp-1) + Axbb(ip  ,jp-1)) / 2 * sx1p * sy1m &
                - (Axb(ip-1,jp+1) + Axbb(ip-1,jp+1)) / 2 * sx1m * sy1p &
                + (Axb(ip-1,jp  ) + Axbb(ip-1,jp  )) / 2 * sx1m * (sy1p-sy1m) &
                + (Axb(ip-1,jp-1) + Axbb(ip-1,jp-1)) / 2 * sx1m * sy1m 
                
         else if(ip .eq. 0 .and. jp .ge. 1) then
            bzz = (Ayb(ip+1,jp  ) + Aybb(ip+1,jp  )) / 2 * sy1p * sx1p &
                - (Ayb(ip  ,jp  ) + Aybb(ip  ,jp  )) / 2 * sy1p * (sx1p-sx1m)&
                - (Ayb(nx-1,jp  ) + Aybb(nx-1,jp  )) / 2 * sy1p * sy1m &
                + (Ayb(ip+1,jp-1) + Aybb(ip+1,jp-1)) / 2 * sy1m * sx1p &
                + (Ayb(ip  ,jp-1) + Aybb(ip  ,jp-1)) / 2 * sy1m * (sx1p-sx1m) &
                + (Ayb(nx-1,jp-1) + Aybb(nx-1,jp-1)) / 2 * sy1m * sx1m &
                - (Axb(ip  ,jp+1) + Axbb(ip  ,jp+1)) / 2 * sx1p * sy1p &
                + (Axb(ip  ,jp  ) + Axbb(ip  ,jp  )) / 2 * sx1p * (sy1p-sy1m) &
                + (Axb(ip  ,jp-1) + Axbb(ip  ,jp-1)) / 2 * sx1p * sy1m &
                - (Axb(nx-1,jp+1) + Axbb(nx-1,jp+1)) / 2 * sx1m * sy1p &
                + (Axb(nx-1,jp  ) + Axbb(nx-1,jp  )) / 2 * sx1m * (sy1p-sy1m) &
                + (Axb(nx-1,jp-1) + Axbb(nx-1,jp-1)) / 2 * sx1m * sy1m 
                
          else if(ip .ge. 1 .and. jp .eq. 0) then
            bzz = (Ayb(ip+1,jp  ) + Aybb(ip+1,jp  )) / 2 * sy1p * sx1p &
                - (Ayb(ip  ,jp  ) + Aybb(ip  ,jp  )) / 2 * sy1p * (sx1p-sx1m)&
                - (Ayb(ip-1,jp  ) + Aybb(ip-1,jp  )) / 2 * sy1p * sy1m &
                + (Ayb(ip+1,ny-1) + Aybb(ip+1,ny-1)) / 2 * sy1m * sx1p &
                + (Ayb(ip  ,ny-1) + Aybb(ip  ,ny-1)) / 2 * sy1m * (sx1p-sx1m) &
                + (Ayb(ip-1,ny-1) + Aybb(ip-1,ny-1)) / 2 * sy1m * sx1m &
                - (Axb(ip  ,jp+1) + Axbb(ip  ,jp+1)) / 2 * sx1p * sy1p &
                + (Axb(ip  ,jp  ) + Axbb(ip  ,jp  )) / 2 * sx1p * (sy1p-sy1m) &
                + (Axb(ip  ,ny-1) + Axbb(ip  ,ny-1)) / 2 * sx1p * sy1m &
                - (Axb(ip-1,jp+1) + Axbb(ip-1,jp+1)) / 2 * sx1m * sy1p &
                + (Axb(ip-1,jp  ) + Axbb(ip-1,jp  )) / 2 * sx1m * (sy1p-sy1m) &
                + (Axb(ip-1,ny-1) + Axbb(ip-1,ny-1)) / 2 * sx1m * sy1m 
          else
            bzz = (Ayb(ip+1,jp  ) + Aybb(ip+1,jp  )) / 2 * sy1p * sx1p &
                - (Ayb(ip  ,jp  ) + Aybb(ip  ,jp  )) / 2 * sy1p * (sx1p-sx1m)&
                - (Ayb(nx-1,jp  ) + Aybb(nx-1,jp  )) / 2 * sy1p * sy1m &
                + (Ayb(ip+1,ny-1) + Aybb(ip+1,ny-1)) / 2 * sy1m * sx1p &
                + (Ayb(ip  ,ny-1) + Aybb(ip  ,ny-1)) / 2 * sy1m * (sx1p-sx1m) &
                + (Ayb(nx-1,ny-1) + Aybb(nx-1,ny-1)) / 2 * sy1m * sx1p &
                - (Axb(ip  ,jp+1) + Axbb(ip  ,jp+1)) / 2 * sx1p * sy1p &
                + (Axb(ip  ,jp  ) + Axbb(ip  ,jp  )) / 2 * sx1p * (sy1p-sy1m) &
                + (Axb(ip  ,ny-1) + Axbb(ip  ,ny-1)) / 2 * sx1p * sy1m &
                - (Axb(nx-1,jp+1) + Axbb(nx-1,jp+1)) / 2 * sx1m * sy1p &
                + (Axb(nx-1,jp  ) + Axbb(nx-1,jp  )) / 2 * sx1m * (sy1p-sy1m) &
                + (Axb(nx-1,ny-1) + Axbb(nx-1,ny-1)) / 2 * sx1m * sy1m 
         end if

    ! bxx = bxg(ip  ,jp  ,kp  )*dx1*dy1*dz1 + bxg(ip+1,jp  ,kp  )*dx*dy1*dz1 &
    !     + bxg(ip  ,jp+1,kp  )*dx1*dy*dz1  + bxg(ip  ,jp  ,kp+1)*dx1*dy1*dz &
    !     + bxg(ip+1,jp+1,kp  )*dx*dy*dz1   + bxg(ip  ,jp+1,kp+1)*dx1*dy*dz &
    !     + bxg(ip+1,jp  ,kp+1)*dx*dy1*dz   + bxg(ip+1,jp+1,kp+1)*dx*dy*dz

    ! byy = byg(ip  ,jp  ,kp  )*dx1*dy1*dz1 + byg(ip+1,jp  ,kp  )*dx*dy1*dz1 &
    !     + byg(ip  ,jp+1,kp  )*dx1*dy*dz1  + byg(ip  ,jp  ,kp+1)*dx1*dy1*dz &
    !     + byg(ip+1,jp+1,kp  )*dx*dy*dz1   + byg(ip  ,jp+1,kp+1)*dx1*dy*dz &
    !     + byg(ip+1,jp  ,kp+1)*dx*dy1*dz   + byg(ip+1,jp+1,kp+1)*dx*dy*dz

    ! bzz = bzg(ip  ,jp  ,kp  )*dx1*dy1*dz1 + bzg(ip+1,jp  ,kp  )*dx*dy1*dz1 &
    !     + bzg(ip  ,jp+1,kp  )*dx1*dy*dz1  + bzg(ip  ,jp  ,kp+1)*dx1*dy1*dz &
    !     + bzg(ip+1,jp+1,kp  )*dx*dy*dz1   + bzg(ip  ,jp+1,kp+1)*dx1*dy*dz &
    !     + bzg(ip+1,jp  ,kp+1)*dx*dy1*dz   + bzg(ip+1,jp+1,kp+1)*dx*dy*dz

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
      real(8) :: chrg, cfact, dx, dy, dz, sx2, sy2, sx2p, sy2p, sx2m, sy2m
      integer :: np, nx, ny, i, ip, jp, kp, ix, iy

!*poption parallel, psum(rho)

      do i = 1, np

         ip = x(i)
         jp = y(i)
         !kp = z(i)

         dx  = x(i) - dble(ip)
         dy  = y(i) - dble(jp)
         !dz  = z(i) - dble(kp)
         !dx1 = 1.0d0 - dx
         !dy1 = 1.0d0 - dy
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

         if( ip .ne. 0  .and. jp .ne. 0) then
            rho(ip-1,jp-1) = rho(ip-1,jp-1) + sx2m * sy2m * chrg
            rho(ip-1,jp  ) = rho(ip-1,jp  ) + sx2m * sy2  * chrg
            rho(ip-1,jp+1) = rho(ip-1,jp+1) + sx2m * sy2p * chrg
            rho(ip  ,jp-1) = rho(ip  ,jp-1) + sx2  * sy2m * chrg
            rho(ip+1,jp-1) = rho(ip+1,jp-1) + sx2p * sy2m * chrg
         else if ( ip .eq. 0 .and. jp .ne. 0) then
            rho(nx-1,jp-1) = rho(nx-1,jp-1) + sx2m * sy2m * chrg
            rho(nx-1,jp  ) = rho(nx-1,jp  ) + sx2m * sy2  * chrg
            rho(nx-1,jp+1) = rho(nx-1,jp+1) + sx2m * sy2p * chrg
            rho(0   ,jp-1) = rho(0   ,jp-1) + sx2  * sy2m * chrg
            rho(1   ,jp-1) = rho(1   ,jp-1) + sx2p * sy2m * chrg
         else if ( ip .ne. 0 .and. jp .eq. 0) then
            rho(ip-1,ny-1) = rho(ip-1,ny-1) + sx2m * sy2m * chrg
            rho(ip-1,0   ) = rho(ip-1,0   ) + sx2m * sy2  * chrg
            rho(ip-1,1   ) = rho(ip-1,1   ) + sx2m * sy2p * chrg
            rho(ip  ,ny-1) = rho(ip  ,ny-1) + sx2  * sy2m * chrg
            rho(ip+1,ny-1) = rho(ip+1,ny-1) + sx2p * sy2m * chrg
         else
            rho(nx-1,ny-1) = rho(nx-1,ny-1) + sx2m * sy2m * chrg
            rho(nx-1,0   ) = rho(nx-1,0   ) + sx2m * sy2  * chrg
            rho(nx-1,1   ) = rho(nx-1,1   ) + sx2m * sy2p * chrg
            rho(0   ,ny-1) = rho(0   ,ny-1) + sx2  * sy2m * chrg
            rho(1   ,ny-1) = rho(1   ,ny-1) + sx2p * sy2m * chrg
         endif
         rho(ip  ,jp  ) = rho(ip  ,jp  ) + sx2  * sy2  * chrg
         rho(ip+1,jp  ) = rho(ip+1,jp  ) + sx2p * sy2  * chrg
         rho(ip  ,jp+1) = rho(ip  ,jp+1) + sx2  * sy2p * chrg
         rho(ip+1,jp+1) = rho(ip+1,jp+1) + sx2p * sy2p * chrg
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
    subroutine phia(nx,ny,c,omega,dt,phi,phib,jx,jy,jz,Ax,Ay,Az,Axb,Ayb,Azb,&
                    Axbb,Aybb,Azbb)
!***********************************************************************
   !original subroutine
      implicit none
      real(8), dimension(0:nx,0:ny) :: phi, phib, jx, jy, jz, Ax, Ay, Az,&
      Axb, Ayb, Azb, Axbb, Aybb, Azbb
      integer :: nx, ny, i, j, k
      real(rkind) :: c, omega, dt
      
 ! Solution of maxwell equation in the A-phi formulation by difference method
      do i = 1, nx-1
      do j = 1, ny-1

      
        Ax(i,j) = dt ** 2 * c ** 2 / omega ** 2 * (Axb(i+1,j) + Axb(i-1,j) &
                + Axb(i,j+1) + Axb(i,j-1) - 4.0d0 * Axb(i,j)) &
                + dt ** 2 * jx(i,j) - dt * (phi(i+1,j) - phib(i+1,j)&
                - phi(i,j) + phib(i,j)) + 2.0d0 * Axb(i,j) - Axbb(i,j) 

        Ay(i,j) = dt ** 2 * c ** 2 / omega ** 2 * (Ayb(i+1,j) + Ayb(i-1,j) &
                + Ayb(i,j+1) + Ayb(i,j-1) - 4.0d0 * Ayb(i,j)) &
                + dt ** 2 * jy(i,j) - dt * (phi(i,j+1) - phib(i,j+1)&
                - phi(i,j) + phib(i,j)) + 2.0d0 * Ayb(i,j) - Aybb(i,j) 

        Az(i,j) = dt ** 2 * c ** 2 / omega ** 2 * (Azb(i+1,j) + Azb(i-1,j) &
                + Azb(i,j+1) + Azb(i,j-1) - 4.0d0 * Azb(i,j)) &
                + dt ** 2 * jz(i,j) + 2.0d0 * Azb(i,j) - Azbb(i,j)
      
      end do
      end do

      do j = 1, ny-1
      
        Ax(0,j) = dt ** 2 * c ** 2 / omega ** 2 * (Axb(1,j) + Axb(nx-1,j) &
                + Axb(0,j+1) + Axb(0,j-1) - 4.0d0 * Axb(0,j)) &
                + dt ** 2 * jx(0,j) - dt * (phi(1,j) - phib(1,j) & 
                - phi(0,j) + phib(0,j)) + 2.0d0 * Axb(0,j) - Axbb(0,j)

        Ax(nx,j) = Ax(0,j)
    
        Ay(0,j) = dt ** 2 * c ** 2 / omega ** 2 * (Ayb(1,j) + Ayb(nx-1,j) &
                + Ayb(0,j+1) + Ayb(0,j-1) - 4.0d0 * Ayb(0,j)) & 
                + dt ** 2 * jy(0,j) - dt * (phi(0,j+1) - phib(0,j+1) & 
                - phi(0,j) + phib(0,j)) + 2.0d0 * Ayb(0,j) - Aybb(0,j)

        Ay(nx,j) = Ay(0,j)
      
        Az(0,j) = dt ** 2 * c ** 2 / omega ** 2 * (Azb(1,j) + Azb(nx-1,j) &
                + Azb(0,j+1) + Azb(0,j-1) - 4.0d0 * Azb(0,j)) &
                + dt ** 2 * jz(0,j) + 2.0d0 * Azb(0,j) - Azbb(0,j)

        Az(nx,j) = Az(0,j)
  
      end do

      do i = 1, nx-1

        Ax(i,0) = dt ** 2 * c ** 2 / omega ** 2 * (Axb(i+1,0) + Axb(i-1,0) &
                + Axb(i,1) + Axb(i,ny-1) - 4.0d0 * Axb(i,0)) &
                + dt ** 2 * jx(i,0) - dt * (phi(i+1,0) - phib(i+1,0)&
                - phi(i,0) + phib(i,0))  + 2.0d0 * Axb(i,0) - Axbb(i,0) 
     
        Ax(i,ny) = Ax(i,0) 

        Ay(i,0) = dt ** 2 * c ** 2 / omega ** 2 * (Ayb(i+1,0) + Ayb(i-1,0) &
                + Ayb(i,1) + Ayb(i,ny-1) - 4.0d0 * Ayb(i,0)) &
                + dt ** 2 * jy(i,0) - dt * (phi(i,1) - phib(i,1) &
                - phi(i,0) + phib(i,0)) + 2.0d0 * Ayb(i,0) - Aybb(i,0)

        Ay(i,ny) = Ay(i,0)
        
        Az(i,0) = dt ** 2 * c ** 2 / omega ** 2 * (Azb(i+1,0) + Azb(i-1,0) &
                + Azb(i,1) + Azb(i,ny-1) - 4.0d0 * Azb(i,0)) &
                + dt ** 2 * jz(i,0) + 2.0d0 * Azb(i,0) - Azbb(i,0)
 
        Az(i,ny) = Az(i,0)
      end do

      Ax(0,0) = dt ** 2 * c ** 2 / omega ** 2 * (Axb(1,0) + Axb(nx-1,0) &
              + Axb(0,1)  + Axb(0,ny-1) - 4.0d0 * Axb(0,0)) &
              + dt ** 2 * jx(0,0) - dt * (phi(1,0) - phib(1,0) & 
              - phi(0,0) + phib(0,0)) + 2.0d0 * Axb(0,0) - Axbb(0,0)
       
      Ax(nx,ny) = Ax(0,0) 
      Ax(0,ny) = Ax(0,0) 
      Ax(nx,0) = Ax(0,0)
      
      Ay(0,0) = dt ** 2 * c ** 2 / omega ** 2 * (Ayb(1,0) + Ayb(nx-1,0) &
              + Ayb(0,1) + Ayb(0,ny-1) - 4.0d0 * Ayb(0,0)) &
              + dt ** 2 * jy(0,0) - dt * (phi(0,1) - phib(0,1) &
              - phi(0,0) + phib(0,0)) + 2.0d0 * Ayb(0,0) - Aybb(0,0)

      Ay(nx,ny) = Ay(0,0)
      Ay(0,ny) = Ay(0,0)
      Ay(nx,0) = Ay(0,0)

      Az(0,0) = dt ** 2 * c ** 2 / omega ** 2 * (Azb(1,0) + Azb(nx-1,0) &
              + Azb(0,1) + Azb(0,ny-1) - 4.0d0 * Azb(0,0)) &
              + dt ** 2 * jz(0,0) + 2.0d0 * Azb(0,0) - Azbb(0,0)

      Az(nx,ny) = Az(0,0) 
      Az(0,ny) = Az(0,0) 
      Az(nx,0) = Az(0,0)
      
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
   subroutine efield(nx,ny,nz,phi,ex,ey,ez,ezg)
!***********************************************************************
     implicit none
      real(8), dimension(0:nx,0:ny) :: phi
      real(8), dimension(0:nx,0:ny,0:nz) :: ex, ey, ezg
      integer :: nx, ny, nz, i, j, im, ip, jm, jp, k, km, kp
      real(rkind) :: ez

      do j = 0, ny
      do i = 0, nx
      do k = 0, nz

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

         ex(i,j,k) = 0.5d0 * ( phi(im,j ) - phi(ip,j ) )
         ey(i,j,k) = 0.5d0 * ( phi(i ,jm) - phi(i ,jp) )
         ezg(i,j,k) = ez
         !ex(i,j,k) = 0
         !ey(i,j,k) = 0

      end do
      end do
      end do
    end subroutine efield
!***********************************************************************
    subroutine bfield(nx,ny,nz,bx,by,bz,bxg,byg,bzg)
!***********************************************************************
      implicit none
      real(rkind), dimension(0:nx,0:ny,0:nz) :: bxg, byg, bzg
      integer :: nx, ny, nz, i, j, k!, im, ip, jm, jp, km, kp
      real(rkind)::bx, by, bz

      do j = 0, ny
      do i = 0, nx
      do k = 0, nz
         !im = i - 1
         !ip = i + 1
         !jm = j - 1
         !jp = j + 1
         !km = k - 1
         !kp = k + 1

         !if( i .eq. 0  ) im = nx - 1
         !if( i .eq. nx ) ip = 1
         !if( j .eq. 0  ) jm = ny - 1
         !if( j .eq. ny ) jp = 1
         !if( k .eq. 0  ) km = nz - 1
         !if( k .eq. nz ) kp = 1
         

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
    subroutine pote(nx,ny,nz,ex,ey,ez,apot,cfacti)
!***********************************************************************
      implicit none
      real(8), dimension(0:nx,0:ny,0:nz) :: ex, ey, ez
      real(8) :: apot, cfacti 
      integer(4) :: nx, ny, nz, ix, iy, iz 

      apot = 0.d0
      do iy = 0, ny-1
      do ix = 0, nx-1
      do iz = 0, nz-1      
         apot = apot + ex(ix,iy,iz)*ex(ix,iy,iz) + ey(ix,iy,iz)*ey(ix,iy,iz)&
              + ez(ix,iy,iz)*ez(ix,iy,iz)
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
    subroutine current(np,nx,ny,nz,x,y,z,xb,yb,zb,jx,jy,jz,dt,chrg)
!***********************************************************************
    implicit none

    real(8), dimension(np) :: x, y, z, xb, yb, zb 
    real(8), dimension(0:nx,0:ny) :: jx, jy, jz
    real(8) :: chrg, dt, dx, dy, dz, dx1, dy1, dz1, deltax, deltay, deltaz,&
               cfact, sx1p, sy1p, sx1m, sy1m, sx2, sy2, sx2p, sy2p, sx2m, sy2m
    integer :: np, nx, ny, nz, i, ip, jp, kp, ix, iy

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
        deltax = x(i) - xb(i)
        deltay = y(i) - yb(i)
        deltaz = z(i) - zb(i)
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
            sx2  = 3.0d0/4 - (dx-1.0d0) ** 2
            sx2p = 1.0d0/2 * (-1.0d0/2 + dx) ** 2
            sx2m = 1.0d0/2 * (3.0d0/2 - dx) ** 2
         endif
         if(dy .le. 0.5d0) then
            sy2  = 3.0d0/4 - dy ** 2
            sy2p = 1.0d0/2 * (1.0d0/2 + dy) ** 2
            sy2m = 1.0d0/2 * (1.0d0/2 - dy) ** 2
         else
            sy2  = 3.0d0/4 - (dy-1.0d0) ** 2
            sy2p = 1.0d0/2 * (-1.0d0/2 + dy) ** 2
            sy2m = 1.0d0/2 * (3.0d0/2 - dy) ** 2
         endif
        if(ip .ne. 0 .and. jp .ne. 0) then
           jx(ip  ,jp  ) = jx(ip  ,jp  ) + chrg / dt * deltax * sy2  * sx1p
           jx(ip  ,jp+1) = jx(ip  ,jp+1) + chrg / dt * deltax * sy2p * sx1p
           jx(ip  ,jp-1) = jx(ip  ,jp-1) + chrg / dt * deltax * sy2m * sx1p
           jx(ip-1,jp  ) = jx(ip-1,jp  ) + chrg / dt * deltax * sy2  * sx1m
           jx(ip-1,jp+1) = jx(ip-1,jp+1) + chrg / dt * deltax * sy2p * sx1m
           jx(ip-1,jp-1) = jx(ip-1,jp-1) + chrg / dt * deltax * sy2m * sx1m
        else if(ip .eq. 0 .and. jp .ne. 0) then
           jx(ip  ,jp  ) = jx(ip  ,jp  ) + chrg / dt * deltax * sy2  * sx1p
           jx(ip  ,jp+1) = jx(ip  ,jp+1) + chrg / dt * deltax * sy2p * sx1p
           jx(ip  ,jp-1) = jx(ip  ,jp-1) + chrg / dt * deltax * sy2m * sx1p
           jx(nx-1,jp  ) = jx(nx-1,jp  ) + chrg / dt * deltax * sy2  * sx1m
           jx(nx-1,jp+1) = jx(nx-1,jp+1) + chrg / dt * deltax * sy2p * sx1m
           jx(nx-1,jp-1) = jx(nx-1,jp-1) + chrg / dt * deltax * sy2m * sx1m
        else if(ip .ne. 0 .and. jp .eq. 0) then 
           jx(ip  ,jp  ) = jx(ip  ,jp  ) + chrg / dt * deltax * sy2  * sx1p
           jx(ip  ,jp+1) = jx(ip  ,jp+1) + chrg / dt * deltax * sy2p * sx1p
           jx(ip  ,ny-1) = jx(ip  ,ny-1) + chrg / dt * deltax * sy2m * sx1p
           jx(ip-1,jp  ) = jx(ip-1,jp  ) + chrg / dt * deltax * sy2  * sx1m
           jx(ip-1,jp+1) = jx(ip-1,jp+1) + chrg / dt * deltax * sy2p * sx1m
           jx(ip-1,ny-1) = jx(ip-1,ny-1) + chrg / dt * deltax * sy2m * sx1m
        else 
           jx(ip  ,jp  ) = jx(ip  ,jp  ) + chrg / dt * deltax * sy2  * sx1p
           jx(ip  ,jp+1) = jx(ip  ,jp+1) + chrg / dt * deltax * sy2p * sx1p
           jx(ip  ,ny-1) = jx(ip  ,ny-1) + chrg / dt * deltax * sy2m * sx1p
           jx(nx-1,jp  ) = jx(nx-1,jp  ) + chrg / dt * deltax * sy2  * sx1m
           jx(nx-1,jp+1) = jx(nx-1,jp+1) + chrg / dt * deltax * sy2p * sx1m
           jx(nx-1,ny-1) = jx(nx-1,ny-1) + chrg / dt * deltax * sy2m * sx1m
        end if 
      
        if(ip .ne. 0 .and. jp .ne. 0) then
           jy(ip  ,jp  ) = jy(ip  ,jp  ) + chrg / dt * deltay * sx2  * sy1p
           jy(ip-1,jp  ) = jy(ip-1,jp  ) + chrg / dt * deltay * sx2m * sy1p
           jy(ip+1,jp  ) = jy(ip+1,jp  ) + chrg / dt * deltay * sx2p * sy1p
           jy(ip  ,jp-1) = jy(ip,  jp-1) + chrg / dt * deltay * sx2  * sy1m
           jy(ip-1,jp-1) = jy(ip-1,jp-1) + chrg / dt * deltay * sx2m * sy1m
           jy(ip+1,jp-1) = jy(ip+1,jp-1) + chrg / dt * deltay * sx2p * sy1m
        else if (ip .eq. 0 .and. jp .ne. 0) then
           jy(ip  ,jp  ) = jy(ip  ,jp  ) + chrg / dt * deltay * sx2  * sy1p
           jy(nx-1,jp  ) = jy(nx-1,jp  ) + chrg / dt * deltay * sx2m * sy1p
           jy(ip+1,jp  ) = jy(ip+1,jp  ) + chrg / dt * deltay * sx2p * sy1p
           jy(ip  ,jp-1) = jy(ip  ,jp-1) + chrg / dt * deltay * sx2  * sy1m
           jy(nx-1,jp-1) = jy(nx-1,jp-1) + chrg / dt * deltay * sx2m * sy1m
           jy(ip+1,jp-1) = jy(ip+1,jp-1) + chrg / dt * deltay * sx2p * sy1m
        else if (ip .ne. 0 .and. jp .eq. 0) then
           jy(ip  ,jp  ) = jy(ip  ,jp  ) + chrg / dt * deltay * sx2  * sy1p
           jy(ip-1,jp  ) = jy(ip-1,jp  ) + chrg / dt * deltay * sx2m * sy1p
           jy(ip+1,jp  ) = jy(ip+1,jp  ) + chrg / dt * deltay * sx2p * sy1p
           jy(ip  ,ny-1) = jy(ip  ,ny-1) + chrg / dt * deltay * sx2  * sy1m
           jy(ip-1,ny-1) = jy(ip-1,ny-1) + chrg / dt * deltay * sx2m * sy1m
           jy(ip+1,ny-1) = jy(ip+1,ny-1) + chrg / dt * deltay * sx2p * sy1m
        else 
           jy(ip  ,jp  ) = jy(ip  ,jp  ) + chrg / dt * deltay * sx2  * sy1p
           jy(nx-1,jp  ) = jy(nx-1,jp  ) + chrg / dt * deltay * sx2m * sy1p
           jy(ip+1,jp  ) = jy(ip+1,jp  ) + chrg / dt * deltay * sx2p * sy1p
           jy(ip  ,ny-1) = jy(ip  ,ny-1) + chrg / dt * deltay * sx2  * sy1m
           jy(nx-1,ny-1) = jy(nx-1,ny-1) + chrg / dt * deltay * sx2m * sy1m
           jy(ip+1,ny-1) = jy(ip+1,ny-1) + chrg / dt * deltay * sx2p * sy1m
        end if
   
        jz(ip  ,jp  ) = jz(ip  ,jp  ) + chrg / dt * deltaz * sx2  * sy2
        jz(ip+1,jp  ) = jz(ip+1,jp  ) + chrg / dt * deltaz * sx2p * sy2
        jz(ip  ,jp+1) = jz(ip  ,jp+1) + chrg / dt * deltaz * sx2  * sy2p
        jz(ip+1,jp+1) = jz(ip+1,jp+1) + chrg / dt * deltaz * sx2p * sy2p

        if(ip .ne. 0 .and. jp .ne. 0) then
           jz(ip-1,jp  ) = jz(ip-1,jp  ) + chrg / dt * deltaz * sx2m * sy2
           jz(ip-1,jp-1) = jz(ip-1,jp-1) + chrg / dt * deltaz * sx2m * sy2m
           jz(ip-1,jp+1) = jz(ip-1,jp+1) + chrg / dt * deltaz * sx2m * sy2p
           jz(ip  ,jp-1) = jz(ip  ,jp-1) + chrg / dt * deltaz * sx2  * sy2m
           jz(ip+1,jp-1) = jz(ip+1,jp-1) + chrg / dt * deltaz * sx2p * sy2m
        else if (ip .eq. 0 .and. jp .ne. 0) then
           jz(nx-1,jp  ) = jz(ip  ,jp  ) + chrg / dt * deltaz * sx2m * sy2
           jz(nx-1,jp-1) = jz(nx-1,jp-1) + chrg / dt * deltaz * sx2m * sy2m
           jz(nx-1,jp+1) = jz(ip+1,jp+1) + chrg / dt * deltaz * sx2m * sy2p
           jz(ip  ,jp-1) = jz(ip  ,jp-1) + chrg / dt * deltaz * sx2  * sy2m
           jz(ip+1,jp-1) = jz(ip+1,jp-1) + chrg / dt * deltaz * sx2p * sy2m
        else if (ip .ne. 0 .and. jp .eq. 0) then
           jz(ip-1,jp  ) = jz(ip-1,jp  ) + chrg / dt * deltaz * sx2m * sy2
           jz(ip-1,ny-1) = jz(ip-1,ny-1) + chrg / dt * deltaz * sx2m * sy2m
           jz(ip-1,jp+1) = jz(ip-1,jp+1) + chrg / dt * deltaz * sx2m * sy2p
           jz(ip  ,ny-1) = jz(ip  ,ny-1) + chrg / dt * deltaz * sx2  * sy2m
           jz(ip+1,ny-1) = jz(ip+1,ny-1) + chrg / dt * deltaz * sx2p * sy2m
        else
           jz(nx-1,jp  ) = jz(nx-1,jp  ) + chrg / dt * deltaz * sx2m * sy2
           jz(nx-1,ny-1) = jz(nx-1,ny-1) + chrg / dt * deltaz * sx2m * sy2m
           jz(nx-1,jp+1) = jz(nx-1,jp+1) + chrg / dt * deltaz * sx2m * sy2p
           jz(ip  ,ny-1) = jz(ip  ,ny-1) + chrg / dt * deltaz * sx2  * sy2m
           jz(ip+1,ny-1) = jz(ip+1,ny-1) + chrg / dt * deltaz * sx2p * sy2m

        end if
     end do
     

   end subroutine current
END Module picexec
