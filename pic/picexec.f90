!./pic  ***** TASK/PIC EXEC *****

Module picexec

  PRIVATE
  PUBLIC pic_exec

CONTAINS

  SUBROUTINE pic_exec(iout)

    USE piccomm
    USE picsub,ONLY: poissn,fftpic,efield,bfield,kine,pote
    USE piclib
    USE libmpi
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER,INTENT(OUT):: iout
    REAL(8)::abc
    INTEGER:: nx,ny,nt

!-----------------------------------------------------------------------
!----- allocate time hisotry data -------------------------------------------
!-----------------------------------------------------------------------

    CALL pic_expand_storage

!-----------------------------------------------------------------------
!----- start of main do-loop -------------------------------------------
!-----------------------------------------------------------------------

    do nt = 1, ntmax

       time = time + dt
       ntcount = ntcount + 1

       do nx = 0, nxmax
       do ny = 0, nymax
          Axbb(nx,ny) = Axb(nx,ny)
          Aybb(nx,ny) = Ayb(nx,ny)
          Azbb(nx,ny) = Azb(nx,ny)
          Axb(nx,ny)  = Ax(nx,ny)
          Ayb(nx,ny)  = Ay(nx,ny)
          Azb(nx,ny)  = Az(nx,ny)
          phib(nx,ny) = phi(nx,ny)
       end do
       end do
   
       !----- charge assignment
       rho(:,:)=0.0d0
       call source(npmax,nxmax,nymax,xe,ye,rho,chrge)
       call source(npmax,nxmax,nymax,xi,yi,rho,chrgi)
       call boundary_rho(nxmax,nymax,rho)
          
       !..... sum charge densities over cores
       call sumdim(nodes,myid,rho,phi,nxymax)

       !----- calculate electric field
       ipssn = 1
       call poissn(nxmax,nymax,nxmaxh1,nxmax1,nymax1, &
                   rho,phi,rhof,phif,awk,afwk,cform,ipssn)

       !----- current assignment
       jx(:,:)=0.d0
       jy(:,:)=0.d0
       jz(:,:)=0.d0

       call current(npmax,nxmax,nymax,xe,ye,vxe,vye,vze,chrge,jx,jy,jz)
       call current(npmax,nxmax,nymax,xi,yi,vxi,vyi,vzi,chrgi,jx,jy,jz)

       call antenna(nxmax,nymax,jxant,jyant,jzant,phxant,phyant,phzant, &
                    omega,time,jx,jy,jz)
       call boundary_j(nxmax,nymax,jx,jy,jz)

       !..... sum chrrent densities over cores
       call sumdim(nodes,myid,jx,Ax,nxymax)
       call sumdim(nodes,myid,jy,Ay,nxymax)
       call sumdim(nodes,myid,jx,Ax,nxymax)

       !.......... calculate vector potential
       call phia(nxmax,nymax,vcfact,dt,phi,phib,jx,jy,jz,Ax,Ay,Az, &
                 Axb,Ayb,Azb,Axbb,Aybb,Azbb)

       !.......... calculate ex and ey and ez
       call efield(nxmax,nymax,dt,phi,Ax,Ay,Az,Axb,Ayb,Azb, &
                               ex,ey,ez,esx,esy,esz,emx,emy,emz)

       !.......... calculate bxg and byg and bzg
       call bfield(nxmax,nymax,Ax,Ay,Az,Axb,Ayb,Azb, &
                               bx,by,bz,bxbg,bybg,bzbg)
         
       if( mod(nt,ntgstep) .eq. 0 ) then
          call kine(npmax,vxe,vye,vze,akine1,me)
          call kine(npmax,vxi,vyi,vzi,akini1,mi)
          call pote(nxmax,nymax,ex,ey,ez,bx,by,bz,vcfact,apot)
          call sumdim1(nodes,myid,akine1,wkword)
          call sumdim1(nodes,myid,akini1,wkword)
          call sumdim1(nodes,myid,apot,wkword)
       endif

       !..... push electrons
       call push(npmax,nxmax,nymax,xe,ye,ze,vxe,vye,vze, &
                 ex,ey,ez,bx,by,bz,dt,ctome,xeb,yeb,zeb, &
                 phi,phib,Axb,Ayb,Azb,Axbb,Aybb,Azbb, &
                 vparae,vperpe)
         
       !..... push ions
       call push(npmax,nxmax,nymax,xi,yi,zi,vxi,vyi,vzi, &
                 ex,ey,ez,bx,by,bz,dt,ctomi,xib,yib,zib, &
                 phi,phib,Axb,Ayb,Azb,Axbb,Aybb,Azbb, &
                 vparai,vperpi)

       !----- treat particles being out of the boundary
       call bound(npmax,xe,ye,ze,x1,x2,y1,y2,z1,z2,alx,aly,alz)
       call bound(npmax,xi,yi,zi,x1,x2,y1,y2,z1,z2,alx,aly,alz)

       !..... diagnostics to check energy conservation
       !.....            after pushing 
       if( mod(nt,ntgstep) .eq. 0 ) then
          call kine(npmax,vxe,vye,vze,akine2,me)
          call kine(npmax,vxi,vyi,vzi,akini2,mi)
          call sumdim1(nodes,myid,akine2,wkword)
          call sumdim1(nodes,myid,akini2,wkword)
          akine = 0.5d0 * ( akine1 + akine2 )
          akini = 0.5d0 * ( akini1 + akini2 )
          aktot = akine + akini
          atot  = aktot + apot

!          if( nt .eq. 1 ) then
!             akine0 = akine
!             akini0 = akini
!             aktot0 = aktot
!             apot0  = apot
!             atot0  = atot
!          endif
            
          akine = akine - akine0
          akini = akini - akini0
          aktot = aktot - aktot0
          apot  = apot  - apot0
          atot  = atot  - atot0

          ntgcount=ntgcount+1
          timet(ntgcount)=time
          akinet(ntgcount)=akine
          akinit(ntgcount)=akini
          aktott(ntgcount)=aktot
          apott(ntgcount)=apot
          atott(ntgcount)=atot
            
       endif
 
       IF( myid .eq. 0 ) THEN
          IF(MOD(ntcount,ntstep).EQ.0) THEN
             WRITE(6,'(I8,1PE12.4,I8,1P3E12.4)') &
                  ntcount,time,ntgcount,aktot,apot,atot
          END IF
       END IF
    end do
!-----------------------------------------------------------------------
!----- end of main do-loop ---------------------------------------------
!-----------------------------------------------------------------------

    !..... output wall clock time
    call mtx_barrier
    wtime2 = mpi_wtime()
    wtime  = wtime2 - wtime1

    if( myid .eq. 0 ) write(*,*) '*** wall clock time = ***', wtime

    ntgmax=ntgcount
    iout=1
  END SUBROUTINE pic_exec

!***********************************************************************
!      expand strage for history data
!
  SUBROUTINE pic_expand_storage
!***********************************************************************
    USE piccomm,ONLY: ntmax,ntgstep,ntgcount,ntgmax, &
                      timet,akinet,akinit,aktott,apott,atott
    IMPLICIT NONE
    REAL(8),DIMENSION(:,:),ALLOCATABLE:: work
    INTEGER:: ntgmax_old

    IF(ntgcount.eq.0) THEN
       IF(ALLOCATED(timet)) &
            DEALLOCATE(timet,akinet,akinit,aktott,apott,atott)
       ntgmax=ntmax/ntgstep
       ALLOCATE(timet(ntgmax))
       ALLOCATE(akinet(ntgmax))
       ALLOCATE(akinit(ntgmax))
       ALLOCATE(aktott(ntgmax))
       ALLOCATE(apott(ntgmax))
       ALLOCATE(atott(ntgmax))
    ELSE
       ALLOCATE(work(ntgmax,6))
       work(1:ntgmax,1)=timet (1:ntgmax)
       work(1:ntgmax,2)=akinet(1:ntgmax)
       work(1:ntgmax,3)=akinit(1:ntgmax)
       work(1:ntgmax,4)=aktott(1:ntgmax)
       work(1:ntgmax,5)=apott (1:ntgmax)
       work(1:ntgmax,6)=atott (1:ntgmax)
       DEALLOCATE(timet,akinet,akinit,aktott,apott,atott)
       ntgmax_old=ntgmax
       ntgmax=ntgmax+ntmax/ntgstep
       ALLOCATE(timet(ntgmax))
       ALLOCATE(akinet(ntgmax))
       ALLOCATE(akinit(ntgmax))
       ALLOCATE(aktott(ntgmax))
       ALLOCATE(apott(ntgmax))
       ALLOCATE(atott(ntgmax))
       timet (1:ntgmax_old)=work(1:ntgmax_old,1)
       akinet(1:ntgmax_old)=work(1:ntgmax_old,2)
       akinit(1:ntgmax_old)=work(1:ntgmax_old,3)
       aktott(1:ntgmax_old)=work(1:ntgmax_old,4)
       apott (1:ntgmax_old)=work(1:ntgmax_old,5)
       atott (1:ntgmax_old)=work(1:ntgmax_old,6)
       DEALLOCATE(work)
    END IF
  END SUBROUTINE pic_expand_storage
      
!***********************************************************************
  subroutine push(npmax,nxmax,nymax,x,y,z,vx,vy,vz,ex,ey,ez,bx,by,bz,dt,&
                  ctom,xb,yb,zb,phi,phib,Axb,Ayb,Azb,Axbb,Aybb,Azbb, &
                  vpara,vperp)
!***********************************************************************
    implicit none
    real(8), dimension(npmax) :: x, y, z, xb, yb, zb
    real(8), dimension(npmax) :: vx, vy, vz, vpara, vperp
    real(8), dimension(0:nxmax,0:nymax) :: phi,phib,Axb,Ayb,Azb,Axbb,Aybb,Azbb
    real(8), dimension(0:nxmax,0:nymax) :: ex, ey, ez, bx, by, bz
    real(8) :: ctom, dx, dy, dx1, dy1, dt, exx, eyy, ezz, bxx,&
               byy, bzz, vxn, vyn, vzn, vxzero, vyzero, vzzero, vxp, vyp,&
               vzp, sx1p, sx1m, sy1p, sy1m, sx2, sy2, sx2p, sx2m, sy2m, sy2p
    real(8) :: btot, vtot
    integer :: npmax, nxmax, nymax
    integer :: np, nx, ny, nxp, nyp, nxpp, nxpm, nypp, nypm, nxppp, nyppp
     
    do np = 1, npmax

! calculate the electric field at the particle position

       nxp = x(np)
       nyp = y(np)
       dx = x(np) - dble(nxp)
       dy = y(np) - dble(nyp)
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

       nxpm  = nxp - 1
       nxpp  = nxp + 1
       nxppp = nxp + 2
       nypm  = nyp - 1
       nypp  = nyp + 1
       nyppp = nyp + 2

       if( nxp .eq. 0  ) nxpm = nxmax - 1
       if( nyp .eq. 0  ) nypm = nymax - 1
       if( nxp .eq. nxmax-2) nxppp=0
       if( nyp .eq. nymax-2) nyppp=0
       if( nxp .eq. nxmax-1) nxpp =0
       if( nyp .eq. nymax-1) nypp =0
       if( nxp .eq. nxmax-1) nxppp = 1
       if( nyp .eq. nymax-1) nyppp = 1

       ! electric field and magnetic field
       if(dx .le. 0.5d0 .and. dy .le. 0.5d0) then
          exx = ex(nxpp,nypp)*dx*sy2p + ex(nxp ,nypp)*dx1*sy2p &
              + ex(nxpp,nyp )*dx*sy2  + ex(nxp ,nyp )*dx1*sy2  &
              + ex(nxpp,nypm)*dx*sy2m + ex(nxp ,nypm)*dx1*sy2m
 
          eyy = ey(nxpp,nypp)*sx2p*dy + ey(nxpp,nyp )*sx2p*dy1 &
              + ey(nxp ,nypp)*sx2 *dy + ey(nxp ,nyp )*sx2 *dy1 &
              + ey(nxpm,nypp)*sx2m*dy + ey(nxpm,nyp )*sx2m*dy1 
            
          ezz = ez(nxpp,nypp)*sx2p*sy2p + ez(nxpp,nyp )*sx2p*sy2 &
              + ez(nxpp,nypm)*sx2p*sy2m + ez(nxp ,nypp)*sx2 *sy2p&
              + ez(nxp ,nyp )*sx2 *sy2  + ez(nxp ,nypm)*sx2 *sy2m&
              + ez(nxpm,nypp)*sx2m*sy2p + ez(nxpm,nyp )*sx2m*sy2 &
              + ez(nxpm,nypm)*sx2m*sy2m
            
          bxx = bx(nxpp,nypp)*dx*sy2p + bx(nxp ,nypp)*dx1*sy2p &
              + bx(nxpp,nyp )*dx*sy2  + bx(nxp ,nyp )*dx1*sy2  &
              + bx(nxpp,nypm)*dx*sy2m + bx(nxp ,nypm)*dx1*sy2m
      
          byy = by(nxpp,nypp)*sx2p*dy + by(nxpp,nyp )*sx2p*dy1 &
              + by(nxp ,nypp)*sx2 *dy + by(nxp ,nyp )*sx2 *dy1 &
              + by(nxpm,nypp)*sx2m*dy + by(nxpm,nyp )*sx2m*dy1
        
          bzz = bz(nxpp,nypp)*sx2p*sy2p + bz(nxp ,nypp)*sx2 *sy2p &
              + bz(nxpm,nypp)*sx2m*sy2p &
              + bz(nxpp,nyp )*sx2p*sy2  + bz(nxp ,nyp )*sx2 *sy2  &
              + bz(nxpm,nyp )*sx2m*sy2  &
              + bz(nxpp,nypm)*sx2p*sy2m + bz(nxp ,nypm)*sx2 *sy2m &
              + bz(nxpm,nypm)*sx2m*sy2m

       else if(dx .le. 0.5d0 .and. dy .ge. 0.5d0) then
          exx = ex(nxpp,nyppp)*dx*sy2p + ex(nxp ,nyppp)*dx1*sy2p &
              + ex(nxpp,nypp )*dx*sy2  + ex(nxp ,nypp )*dx1*sy2  &
              + ex(nxpp,nyp  )*dx*sy2m + ex(nxp ,nyp  )*dx1*sy2m
      
          eyy = ey(nxpp,nypp)*sx2p*dy + ey(nxpp,nyp )*sx2p*dy1 &
              + ey(nxp ,nypp)*sx2 *dy + ey(nxp ,nyp )*sx2 *dy1 &
              + ey(nxpm,nypp)*sx2m*dy + ey(nxpm,nyp )*sx2m*dy1 
            
          ezz = ez(nxpp,nyppp)*sx2p*sy2p + ez(nxpp,nypp )*sx2p*sy2 &
              + ez(nxpp,nyp  )*sx2p*sy2m + ez(nxp ,nyppp)*sx2 *sy2p&
              + ez(nxp ,nypp )*sx2 *sy2  + ez(nxp ,nyp  )*sx2 *sy2m&
              + ez(nxpm,nyppp)*sx2m*sy2p + ez(nxpm,nypp )*sx2m*sy2 &
              + ez(nxpm,nyp  )*sx2m*sy2m
      
          bxx = bx(nxpp,nyppp)*dx*sy2p + bx(nxp ,nyppp)*dx1*sy2p &
              + bx(nxpp,nypp )*dx*sy2  + bx(nxp ,nypp )*dx1*sy2  &
              + bx(nxpp,nyp  )*dx*sy2m + bx(nxp ,nyp  )*dx1*sy2m
        
          byy = by(nxpp,nypp)*sx2p*dy + by(nxpp,nyp )*sx2p*dy1 &
              + by(nxp ,nypp)*sx2 *dy + by(nxp ,nyp )*sx2 *dy1 &
              + by(nxpm,nypp)*sx2m*dy + by(nxpm,nyp )*sx2m*dy1
        
          bzz = bz(nxpp,nyppp)*sx2p*sy2p + bz(nxp ,nyppp)*sx2 *sy2p &
              + bz(nxpm,nyppp)*sx2m*sy2p + bz(nxpp,nypp )*sx2p*sy2  &
              + bz(nxp ,nypp )*sx2 *sy2  + bz(nxpm,nypp )*sx2m*sy2  &
              + bz(nxpp,nyp  )*sx2p*sy2m + bz(nxp ,nyp  )*sx2 *sy2m &
              + bz(nxpm,nyp  )*sx2m*sy2m

       else if(dx .ge. 0.5d0 .and. dy .le. 0.5d0) then
          exx = ex(nxpp,nypp)*dx*sy2p + ex(nxp ,nypp)*dx1*sy2p &
              + ex(nxpp,nyp )*dx*sy2  + ex(nxp ,nyp )*dx1*sy2  &
              + ex(nxpp,nypm)*dx*sy2m + ex(nxp ,nypm)*dx1*sy2m
        
          eyy = ey(nxppp,nypp)*sx2p*dy + ey(nxppp,nyp )*sx2p*dy1 &
              + ey(nxpp ,nypp)*sx2 *dy + ey(nxpp ,nyp )*sx2 *dy1 &
              + ey(nxp  ,nypp)*sx2m*dy + ey(nxp  ,nyp )*sx2m*dy1 
            
          ezz = ez(nxppp,nypp)*sx2p*sy2p + ez(nxppp,nyp )*sx2p*sy2 &
              + ez(nxppp,nypm)*sx2p*sy2m + ez(nxpp ,nypp)*sx2 *sy2p&
              + ez(nxpp ,nyp )*sx2 *sy2  + ez(nxpp ,nypm)*sx2 *sy2m&
              + ez(nxp  ,nypp)*sx2m*sy2p + ez(nxp  ,nyp )*sx2m*sy2 &
              + ez(nxp  ,nypm)*sx2m*sy2m
            
          bxx = bx(nxpp,nypp)*dx*sy2p + bx(nxp ,nypp)*dx1*sy2p &
              + bx(nxpp,nyp )*dx*sy2  + bx(nxp ,nyp )*dx1*sy2  &
              + bx(nxpp,nypm)*dx*sy2m + bx(nxp ,nypm)*dx1*sy2m

          byy = by(nxppp,nypp)*sx2p*dy + by(nxppp,nyp )*sx2p*dy1 &
              + by(nxpp ,nypp)*sx2 *dy + by(nxpp ,nyp )*sx2 *dy1 &
              + by(nxp  ,nypp)*sx2m*dy + by(nxp  ,nyp )*sx2m*dy1
      
          bzz = bz(nxppp,nypp)*sx2p*sy2p + bz(nxpp ,nypp)*sx2 *sy2p &
              + bz(nxp  ,nypp)*sx2m*sy2p + bz(nxppp,nyp )*sx2p*sy2  &
              + bz(nxpp ,nyp )*sx2 *sy2  + bz(nxp  ,nyp )*sx2m*sy2  &
              + bz(nxppp,nypm)*sx2p*sy2m + bz(nxpp ,nypm)*sx2 *sy2m &
              + bz(nxp  ,nypm)*sx2m*sy2m

       else
          exx = ex(nxpp,nyppp)*dx*sy2p + ex(nxp ,nyppp)*dx1*sy2p &
              + ex(nxpp,nypp )*dx*sy2  + ex(nxp ,nypp )*dx1*sy2  &
              + ex(nxpp,nyp  )*dx*sy2m + ex(nxp ,nyp  )*dx1*sy2m
        
          eyy = ey(nxppp,nypp)*sx2p*dy + ey(nxppp,nyp )*sx2p*dy1 &
              + ey(nxpp ,nypp)*sx2 *dy + ey(nxpp ,nyp )*sx2 *dy1 &
              + ey(nxp  ,nypp)*sx2m*dy + ey(nxp  ,nyp )*sx2m*dy1 
            
          ezz = ez(nxppp,nyppp)*sx2p*sy2p + ez(nxppp,nypp )*sx2p*sy2 &
              + ez(nxppp,nyp  )*sx2p*sy2m + ez(nxpp ,nyppp)*sx2 *sy2p&
              + ez(nxpp ,nypp )*sx2 *sy2  + ez(nxpp ,nyp  )*sx2 *sy2m&
              + ez(nxp  ,nyppp)*sx2m*sy2p + ez(nxp  ,nypp )*sx2m*sy2 &
              + ez(nxp  ,nyp  )*sx2m*sy2m
          
          bxx = bx(nxpp,nyppp)*dx*sy2p + bx(nxp ,nyppp)*dx1*sy2p &
              + bx(nxpp,nypp )*dx*sy2  + bx(nxp ,nypp )*dx1*sy2  &
              + bx(nxpp,nyp  )*dx*sy2m + bx(nxp ,nyp  )*dx1*sy2m
        
          byy = by(nxppp,nypp)*sx2p*dy + by(nxppp,nyp )*sx2p*dy1 &
              + by(nxpp ,nypp)*sx2 *dy + by(nxpp ,nyp )*sx2 *dy1 &
              + by(nxp  ,nypp)*sx2m*dy + by(nxp  ,nyp )*sx2m*dy1
        
          bzz = bz(nxppp,nyppp)*sx2p*sy2p + bz(nxpp ,nyppp)*sx2 *sy2p &
              + bz(nxp  ,nyppp)*sx2m*sy2p + bz(nxppp,nypp )*sx2p*sy2  &
              + bz(nxpp ,nypp )*sx2 *sy2  + bz(nxp  ,nypp )*sx2m*sy2  &
              + bz(nxppp,nyp  )*sx2p*sy2m + bz(nxpp ,nyp  )*sx2 *sy2m &
              + bz(nxp  ,nyp  )*sx2m*sy2m
       endif
        
       ! push particles by using Buneman-Boris method

       vxn = vx(np) + 1.0d0/2 * ctom * exx * dt 
       vyn = vy(np) + 1.0d0/2 * ctom * eyy * dt
       vzn = vz(np) + 1.0d0/2 * ctom * ezz * dt

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
 
       vx(np) = vxp + 1.0d0/2 * ctom * exx * dt
       vy(np) = vyp + 1.0d0/2 * ctom * eyy * dt
       vz(np) = vzp + 1.0d0/2 * ctom * ezz * dt
         
       xb(np) = x(np)
       yb(np) = y(np)
       zb(np) = z(np)
       x(np) = x(np) + vx(np) * dt
       y(np) = y(np) + vy(np) * dt
       z(np) = z(np) + vz(np) * dt

       btot=SQRT(bxx**2+byy**2+bzz**2)
       IF(btot.EQ.0.D0) THEN
          vpara(np)=vx(np)
          vperp(np)=SQRT(vy(np)**2+vz(np)**2)
       ELSE
          vtot=SQRT(vx(np)**2+vy(np)**2+vz(np)**2)
          vpara(np)=(bxx*vx(np)+byy*vy(np)+bzz*vz(np))/btot
          vperp(np)=(SQRT(vtot**2-vpara(np)**2))
       END IF

    end do

  end subroutine push

!***********************************************************************
    subroutine bound(npmax,x,y,z,x1,x2,y1,y2,z1,z2,alx,aly,alz)
!***********************************************************************
      implicit none
      real(8), dimension(npmax) :: x, y, z
      real(8) :: alx, aly, alz, x1, x2, y1, y2, z1, z2
      integer :: npmax, np

      do np = 1, npmax
         if( x(np) .lt. x1 ) then
            do while(x(np) .lt. x1)
               x(np) = x(np) + alx
            end do
         elseif( x(np) .gt. x2 ) then
            do while(x(np) .gt. x2)
               x(np) = x(np) - alx
            end do
         endif
 
         if( y(np) .lt. y1 ) then
            do while(y(np) .lt. y1)
               y(np) = y(np) + aly
            end do
         elseif( y(np) .gt. y2 ) then
            do while(y(np) .gt. y2)
               y(np) = y(np) - aly
            end do
         endif

         if( z(np) .lt. z1 ) then
            do while(z(np) .lt. z1)
               z(np) = z(np) + alz
            end do
         elseif( z(np) .gt. z2 ) then
            do while(z(np) .gt. z2)
               z(np) = z(np) - alz
            end do
         endif
      end do

    end subroutine bound

!***********************************************************************
    subroutine source(npmax,nxmax,nymax,x,y,rho,chrg)
!***********************************************************************

      implicit none
      real(8), dimension(npmax)        :: x, y
      real(8), dimension(0:nxmax,0:nymax) :: rho
      real(8) :: chrg, factor, dx, dy, dz, sx2, sy2, sx2p, sy2p, sx2m, sy2m,&
           dx1,dy1
      integer :: npmax, nxmax, nymax
      integer:: np, nxp, nyp, nx, ny, nxpp, nxpm, nypp, nypm, nxppp, nyppp

      factor=chrg*dble(nxmax)*dble(nymax)/dble(npmax)

!*poption parallel, psum(rho)

      do np = 1, npmax

         nxp = x(np)
         nyp = y(np)

         dx  = x(np) - dble(nxp)
         dy  = y(np) - dble(nyp)
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

         nxpm = nxp - 1
         nxpp = nxp + 1
         nxppp= nxp + 2
         nypm = nyp - 1
         nypp = nyp + 1
         nyppp= nyp + 2

         if( nxp .eq. 0  ) nxpm = nxmax - 1
         if( nyp .eq. 0  ) nypm = nymax - 1
         if( nxp .eq. nxmax-2) nxppp=0
         if( nyp .eq. nymax-2) nyppp=0
         if( nxp .eq. nxmax-1) nxpp =0
         if( nyp .eq. nymax-1) nypp =0
         if( nxp .eq. nxmax-1) nxppp=1
         if( nyp .eq. nymax-1) nyppp=1

         if(dx .le. 0.5d0 .and. dy .le. 0.5d0) then
            rho(nxpm,nypm) = rho(nxpm,nypm) + sx2m * sy2m * factor
            rho(nxpm,nyp ) = rho(nxpm,nyp ) + sx2m * sy2  * factor
            rho(nxpm,nypp) = rho(nxpm,nypp) + sx2m * sy2p * factor
            rho(nxp ,nypm) = rho(nxp ,nypm) + sx2  * sy2m * factor
            rho(nxp ,nyp ) = rho(nxp ,nyp ) + sx2  * sy2  * factor
            rho(nxp ,nypp) = rho(nxp ,nypp) + sx2  * sy2p * factor
            rho(nxpp,nypm) = rho(nxpp,nypm) + sx2p * sy2m * factor
            rho(nxpp,nyp ) = rho(nxpp,nyp ) + sx2p * sy2  * factor
            rho(nxpp,nypp) = rho(nxpp,nypp) + sx2p * sy2p * factor
         else if(dx .le. 0.5d0 .and. dy .ge. 0.5d0) then
            rho(nxpm ,nyp  ) = rho(nxpm ,nyp  ) + sx2m * sy2m * factor
            rho(nxpm ,nypp ) = rho(nxpm ,nypp ) + sx2m * sy2  * factor
            rho(nxpm ,nyppp) = rho(nxpm ,nyppp) + sx2m * sy2p * factor
            rho(nxp  ,nyp  ) = rho(nxp  ,nyp  ) + sx2  * sy2m * factor
            rho(nxp  ,nypp ) = rho(nxp  ,nypp ) + sx2  * sy2  * factor
            rho(nxp  ,nyppp) = rho(nxp  ,nyppp) + sx2  * sy2p * factor
            rho(nxpp ,nyp  ) = rho(nxpp ,nyp  ) + sx2p * sy2m * factor
            rho(nxpp ,nypp ) = rho(nxpp ,nypp ) + sx2p * sy2  * factor
            rho(nxpp ,nyppp) = rho(nxpp ,nyppp) + sx2p * sy2p * factor
         else if(dx .ge. 0.5d0 .and. dy .le. 0.5d0) then
            rho(nxp  ,nypm) = rho(nxp  ,nypm ) + sx2m * sy2m * factor
            rho(nxp  ,nyp ) = rho(nxp  ,nyp  ) + sx2m * sy2  * factor
            rho(nxp  ,nypp) = rho(nxp  ,nypp ) + sx2m * sy2p * factor
            rho(nxpp ,nypm) = rho(nxpp ,nypm ) + sx2  * sy2m * factor
            rho(nxpp ,nyp ) = rho(nxpp ,nyp  ) + sx2  * sy2  * factor
            rho(nxpp ,nypp) = rho(nxpp ,nypp ) + sx2  * sy2p * factor
            rho(nxppp,nypm) = rho(nxppp,nypm ) + sx2p * sy2m * factor
            rho(nxppp,nyp ) = rho(nxppp,nyp  ) + sx2p * sy2  * factor
            rho(nxppp,nypp) = rho(nxppp,nypp ) + sx2p * sy2p * factor
         else 
            rho(nxp  ,nyp  ) = rho(nxp  ,nyp  ) + sx2m * sy2m * factor
            rho(nxp  ,nypp ) = rho(nxp  ,nypp ) + sx2m * sy2  * factor
            rho(nxp  ,nyppp) = rho(nxp  ,nyppp) + sx2m * sy2p * factor
            rho(nxpp ,nyp  ) = rho(nxpp ,nyp  ) + sx2  * sy2m * factor
            rho(nxpp ,nypp ) = rho(nxpp ,nypp ) + sx2  * sy2  * factor
            rho(nxpp ,nyppp) = rho(nxpp ,nyppp) + sx2  * sy2p * factor
            rho(nxppp,nyp  ) = rho(nxppp,nyp  ) + sx2p * sy2m * factor
            rho(nxppp,nypp ) = rho(nxppp,nypp ) + sx2p * sy2  * factor
            rho(nxppp,nyppp) = rho(nxppp,nyppp) + sx2p * sy2p * factor
         endif
         end do
     end subroutine source

!***********************************************************************
    subroutine boundary_rho(nxmax,nymax,rho)
!***********************************************************************
      implicit none
      real(8), dimension(0:nxmax,0:nymax) :: rho
      integer:: nxmax,nymax
      integer:: nx,ny

      !!..... set charge densities at the boundary

      do ny = 0, nymax
         rho(0,ny) = rho(0,ny) + rho(nxmax,ny)
         rho(nxmax,ny) = rho(0,ny)
      end do

      do nx = 0, nxmax
         rho(nx,0) = rho(nx,0) + rho(nx,nymax)
         rho(nx,nymax) = rho(nx,0)
      end do
    end subroutine boundary_rho

!***********************************************************************
    subroutine current(npmax,nxmax,nymax,x,y,vx,vy,vz,chrg,jx,jy,jz)
!***********************************************************************
    implicit none

    real(8), dimension(npmax) :: x, y, vx, vy, vz 
    real(8), dimension(0:nxmax,0:nymax) :: jx, jy, jz
    real(8) :: chrg, dt, dx, dy, dx1, dy1, &
         sx1p, sy1p, sx1m, sy1m, sx2, sy2, sx2p, sy2p, sx2m, sy2m, factor
    integer :: npmax, nxmax, nymax
    integer :: np, nxp, nyp, nxpm, nypm, nxpp, nypp, nxppp, nyppp, nx, ny

    factor=chrg*dble(nxmax)*dble(nymax)/dble(npmax)

    do np = 1, npmax

        nxp = x(np)
        nyp = y(np)
        dx = x(np) - dble(nxp)
        dy = y(np) - dble(nyp)
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
         nxpm = nxp - 1
         nxpp = nxp + 1
         nypm = nyp - 1
         nypp = nyp + 1
         nxppp = nxp + 2
         nyppp = nyp + 2
         if( nxp .eq. 0  )       nxpm = nxmax-1
         if( nyp .eq. 0  )       nypm = nymax-1
         if( nxp .eq. nxmax-2  ) nxppp = 0
         if( nyp .eq. nymax-2  ) nyppp = 0
         if( nxp .eq. nxmax-1  ) nxppp = 1
         if( nyp .eq. nymax-1  ) nyppp = 1
        
         if (dx .le. 0.5d0 .and. dy .le. 0.5d0) then
           jx(nxpp,nyp ) = jx(nxpp,nyp ) + factor * vx(np) * sy2  * dx
           jx(nxpp,nypp) = jx(nxpp,nypp) + factor * vx(np) * sy2p * dx
           jx(nxpp,nypm) = jx(nxpp,nypm) + factor * vx(np) * sy2m * dx
           jx(nxp ,nyp ) = jx(nxp ,nyp ) + factor * vx(np) * sy2  * dx1
           jx(nxp ,nypp) = jx(nxp ,nypp) + factor * vx(np) * sy2p * dx1
           jx(nxp ,nypm) = jx(nxp ,nypm) + factor * vx(np) * sy2m * dx1
           
           jy(nxp ,nypp) = jy(nxp ,nypp) + factor * vy(np) * sx2  * dy
           jy(nxpp,nypp) = jy(nxpp,nypp) + factor * vy(np) * sx2p * dy
           jy(nxpm,nypp) = jy(nxpm,nypp) + factor * vy(np) * sx2m * dy
           jy(nxp ,nyp ) = jy(nxp ,nyp ) + factor * vy(np) * sx2  * dy1
           jy(nxpp,nyp ) = jy(nxpp,nyp ) + factor * vy(np) * sx2p * dy1
           jy(nxpm,nyp ) = jy(nxpm,nyp ) + factor * vy(np) * sx2m * dy1
      
           jz(nxpm,nypp) = jz(nxpm,nypp) + factor * vz(np) * sx2m * sy2p
           jz(nxpm,nyp ) = jz(nxpm,nyp ) + factor * vz(np) * sx2m * sy2
           jz(nxpm,nypm) = jz(nxpm,nypm) + factor * vz(np) * sx2m * sy2m
           jz(nxp ,nypp) = jz(nxp ,nypp) + factor * vz(np) * sx2  * sy2p
           jz(nxp ,nyp ) = jz(nxp ,nyp ) + factor * vz(np) * sx2  * sy2
           jz(nxp ,nypm) = jz(nxp ,nypm) + factor * vz(np) * sx2  * sy2m
           jz(nxpp,nypp) = jz(nxpp,nypp) + factor * vz(np) * sx2p * sy2p
           jz(nxpp,nyp ) = jz(nxpp,nyp ) + factor * vz(np) * sx2p * sy2
           jz(nxpp,nypm) = jz(nxpp,nypm) + factor * vz(np) * sx2p * sy2m

        else if(dx .le. 0.5d0 .and. dy .ge. 0.5d0) then
           jx(nxpp,nypp ) = jx(nxpp,nypp ) + factor * vx(np) * sy2  * dx
           jx(nxpp,nyppp) = jx(nxpp,nyppp) + factor * vx(np) * sy2p * dx
           jx(nxpp,nyp  ) = jx(nxpp,nyp  ) + factor * vx(np) * sy2m * dx
           jx(nxp ,nypp ) = jx(nxp ,nypp ) + factor * vx(np) * sy2  * dx1
           jx(nxp ,nyppp) = jx(nxp ,nyppp) + factor * vx(np) * sy2p * dx1
           jx(nxp ,nyp  ) = jx(nxp ,nyp  ) + factor * vx(np) * sy2m * dx1
           
           jy(nxp ,nypp) = jy(nxp ,nypp) + factor * vy(np) * sx2  * dy
           jy(nxpp,nypp) = jy(nxpp,nypp) + factor * vy(np) * sx2p * dy
           jy(nxpm,nypp) = jy(nxpm,nypp) + factor * vy(np) * sx2m * dy
           jy(nxp ,nyp ) = jy(nxp ,nyp ) + factor * vy(np) * sx2  * dy1
           jy(nxpp,nyp ) = jy(nxpp,nyp ) + factor * vy(np) * sx2p * dy1
           jy(nxpm,nyp ) = jy(nxpm,nyp ) + factor * vy(np) * sx2m * dy1
      
           jz(nxpm,nyppp) = jz(nxpm,nyppp) + factor * vz(np) * sx2m * sy2p
           jz(nxpm,nypp ) = jz(nxpm,nypp ) + factor * vz(np) * sx2m * sy2
           jz(nxpm,nyp  ) = jz(nxpm,nyp  ) + factor * vz(np) * sx2m * sy2m
           jz(nxp ,nyppp) = jz(nxp ,nyppp) + factor * vz(np) * sx2  * sy2p
           jz(nxp ,nypp ) = jz(nxp ,nypp ) + factor * vz(np) * sx2  * sy2
           jz(nxp ,nyp  ) = jz(nxp ,nyp  ) + factor * vz(np) * sx2  * sy2m
           jz(nxpp,nyppp) = jz(nxpp,nyppp) + factor * vz(np) * sx2p * sy2p
           jz(nxpp,nypp ) = jz(nxpp,nypp ) + factor * vz(np) * sx2p * sy2
           jz(nxpp,nyp  ) = jz(nxpp,nyp  ) + factor * vz(np) * sx2p * sy2m

        else if(dx .ge. 0.5d0 .and. dy .le. 0.5d0) then
           jx(nxpp,nyp ) = jx(nxpp,nyp ) + factor * vx(np) * sy2  * dx
           jx(nxpp,nypp) = jx(nxpp,nypp) + factor * vx(np) * sy2p * dx
           jx(nxpp,nypm) = jx(nxpp,nypm) + factor * vx(np) * sy2m * dx
           jx(nxp ,nyp ) = jx(nxp ,nyp ) + factor * vx(np) * sy2  * dx1
           jx(nxp ,nypp) = jx(nxp ,nypp) + factor * vx(np) * sy2p * dx1
           jx(nxp ,nypm) = jx(nxp ,nypm) + factor * vx(np) * sy2m * dx1
           
           jy(nxpp ,nypp) = jy(nxpp ,nypp) + factor * vy(np) * sx2  * dy
           jy(nxppp,nypp) = jy(nxppp,nypp) + factor * vy(np) * sx2p * dy
           jy(nxp  ,nypp) = jy(nxp  ,nypp) + factor * vy(np) * sx2m * dy
           jy(nxpp ,nyp ) = jy(nxpp ,nyp ) + factor * vy(np) * sx2  * dy1
           jy(nxppp,nyp ) = jy(nxppp,nyp ) + factor * vy(np) * sx2p * dy1
           jy(nxp  ,nyp ) = jy(nxp  ,nyp ) + factor * vy(np) * sx2m * dy1
      
           jz(nxp  ,nypp) = jz(nxp  ,nypp) + factor * vz(np) * sx2m * sy2p
           jz(nxp  ,nyp ) = jz(nxp  ,nyp ) + factor * vz(np) * sx2m * sy2
           jz(nxp  ,nypm) = jz(nxp  ,nypm) + factor * vz(np) * sx2m * sy2m
           jz(nxpp ,nypp) = jz(nxpp ,nypp) + factor * vz(np) * sx2  * sy2p
           jz(nxpp ,nyp ) = jz(nxpp ,nyp ) + factor * vz(np) * sx2  * sy2
           jz(nxpp ,nypm) = jz(nxpp ,nypm) + factor * vz(np) * sx2  * sy2m
           jz(nxppp,nypp) = jz(nxppp,nypp) + factor * vz(np) * sx2p * sy2p
           jz(nxppp,nyp ) = jz(nxppp,nyp ) + factor * vz(np) * sx2p * sy2
           jz(nxppp,nypm) = jz(nxppp,nypm) + factor * vz(np) * sx2p * sy2m

        else
           jx(nxpp,nypp ) = jx(nxpp,nypp ) + factor * vx(np) * sy2  * dx
           jx(nxpp,nyppp) = jx(nxpp,nyppp) + factor * vx(np) * sy2p * dx
           jx(nxpp,nyp  ) = jx(nxpp,nyp  ) + factor * vx(np) * sy2m * dx
           jx(nxp ,nypp ) = jx(nxp ,nypp ) + factor * vx(np) * sy2  * dx1
           jx(nxp ,nyppp) = jx(nxp ,nyppp) + factor * vx(np) * sy2p * dx1
           jx(nxp ,nyp  ) = jx(nxp ,nyp  ) + factor * vx(np) * sy2m * dx1
           
           jy(nxpp ,nypp) = jy(nxpp ,nypp) + factor * vy(np) * sx2  * dy
           jy(nxppp,nypp) = jy(nxppp,nypp) + factor * vy(np) * sx2p * dy
           jy(nxp  ,nypp) = jy(nxp  ,nypp) + factor * vy(np) * sx2m * dy
           jy(nxpp ,nyp ) = jy(nxpp ,nyp ) + factor * vy(np) * sx2  * dy1
           jy(nxppp,nyp ) = jy(nxppp,nyp ) + factor * vy(np) * sx2p * dy1
           jy(nxp  ,nyp ) = jy(nxp  ,nyp ) + factor * vy(np) * sx2m * dy1
      
           jz(nxp  ,nyppp) = jz(nxp  ,nyppp) + factor * vz(np) * sx2m * sy2p
           jz(nxp  ,nypp ) = jz(nxp  ,nypp ) + factor * vz(np) * sx2m * sy2
           jz(nxp  ,nyp  ) = jz(nxp  ,nyp  ) + factor * vz(np) * sx2m * sy2m
           jz(nxpp ,nyppp) = jz(nxpp ,nyppp) + factor * vz(np) * sx2  * sy2p
           jz(nxpp ,nypp ) = jz(nxpp ,nypp ) + factor * vz(np) * sx2  * sy2
           jz(nxpp ,nyp  ) = jz(nxpp ,nyp  ) + factor * vz(np) * sx2  * sy2m
           jz(nxppp,nyppp) = jz(nxppp,nyppp) + factor * vz(np) * sx2p * sy2p
           jz(nxppp,nypp ) = jz(nxppp,nypp ) + factor * vz(np) * sx2p * sy2
           jz(nxppp,nyp  ) = jz(nxppp,nyp  ) + factor * vz(np) * sx2p * sy2m

        endif
       
     end do

   end subroutine current

!***********************************************************************
    subroutine antenna(nxmax,nymax,jxant,jyant,jzant,phxant,phyant,phzant, &
                       omega,time,jx,jy,jz)
!***********************************************************************
    implicit none

    integer :: nxmax,nymax
    real(8) :: jxant,jyant,jzant,phxant,phyant,phzant,omega,time
    real(8), dimension(0:nxmax,0:nymax) :: jx, jy, jz
    real(8) :: jxt,jyt,jzt
    integer :: nx,ny
     
         jxt = jxant * cos (omega * time + phxant)
         jyt = jyant * cos (omega * time + phyant)
         jzt = jzant * cos (omega * time + phzant)

         ! add antenna current density

         do ny=1,nymax
            jy(0,ny) = jy(0,ny) + 0.5d0 * jyt
            jy(1,ny) = jy(1,ny) + 0.5d0 * jyt
            jz(0,ny) = jz(0,ny) + 0.5d0 * jzt
            jz(1,ny) = jz(1,ny) + 0.5d0 * jzt
         end do

   end subroutine antenna

!***********************************************************************
    subroutine boundary_j(nxmax,nymax,jx,jy,jz)
!***********************************************************************
    implicit none

    integer :: nxmax,nymax
    real(8), dimension(0:nxmax,0:nymax) :: jx, jy, jz
    integer :: nx,ny
     
         !--------------
         do ny = 0, nymax
            jx(0,ny) = jx(0,ny) + jx(nxmax,ny)
            jy(0,ny) = jy(0,ny) + jy(nxmax,ny)
            jz(0,ny) = jz(0,ny) + jz(nxmax,ny)
            jx(nxmax,ny) = jx(0,ny)
            jy(nxmax,ny) = jy(0,ny)
            jz(nxmax,ny) = jz(0,ny)
         end do
         do nx = 0, nxmax
            jx(nx,0) = jx(nx,0) + jx(nx,nymax)
            jy(nx,0) = jy(nx,0) + jy(nx,nymax) 
            jz(nx,0) = jz(nx,0) + jz(nx,nymax)
            jx(nx,nymax) = jx(nx,0)
            jy(nx,nymax) = jy(nx,0)
            jz(nx,nymax) = jz(nx,0)
         enddo

   end subroutine boundary_j

!***********************************************************************
    subroutine phia(nxmax,nymax,vcfact,dt,phi,phib,jx,jy,jz,Ax,Ay,Az, &
                    Axb,Ayb,Azb,Axbb,Aybb,Azbb)
!***********************************************************************
   !original subroutine
      implicit none
      real(8), dimension(0:nxmax,0:nymax) :: phi,phib,jx,jy,jz,Ax,Ay,Az, &
                                             Axb,Ayb,Azb,Axbb,Aybb,Azbb
      integer :: nxmax, nymax, nx, ny, nxm, nxp, nyp, nym, nypm
      real(8) :: vcfact, dt

 ! Solution of maxwell equation in the A-phi formulation by difference method
 ! c is the ratio of the light speed to lattice parameter times plasma frequency
      do nx = 0, nxmax
      do ny = 0, nymax

         nxm = nx - 1
         nxp = nx + 1
         nym = ny - 1
         nyp = ny + 1

         if( nx .eq. 0  )    nxm = nxmax - 1
         if( nx .eq. nxmax ) nxp = 1
         if( ny .eq. 0  )    nym = nymax - 1
         if( ny .eq. nymax ) nyp = 1
      
        Ax(nx,ny) = dt ** 2 * vcfact ** 2 * (Axb(nxp,ny) + Axb(nxm,ny) &
                                           + Axb(nx,nyp) + Axb(nx,nym) &
                                           - 4.0d0 * Axb(nx,ny)) &
                  + dt ** 2 * jx(nx,ny) &
                  - 0.5d0 * dt * (phi(nxp,ny) - phib(nxp,ny) &
                                - phi(nxm,ny) + phib(nxm,ny)) &
                  + 2.0d0 * Axb(nx,ny) - Axbb(nx,ny) 

        Ay(nx,ny) = dt ** 2 * vcfact ** 2 * (Ayb(nxp,ny) + Ayb(nxm,ny) &
                                           + Ayb(nx,nyp) + Ayb(nx,nym) &
                                           - 4.0d0 * Ayb(nx,ny)) &
                  + dt ** 2 * jy(nx,ny) &
                  - 0.5d0 * dt * (phi(nx,nyp) - phib(nx,nyp) &
                                - phi(nx,nym) + phib(nx,nym)) &
                  + 2.0d0 * Ayb(nx,ny) - Aybb(nx,ny) 

        Az(nx,ny) = dt ** 2 * vcfact ** 2 * (Azb(nxp,ny) + Azb(nxm,ny) &
                                           + Azb(nx,nyp) + Azb(nx,nym) &
                                           - 4.0d0 * Azb(nx,ny)) &
                  + dt ** 2 * jz(nx,ny) &
                  + 2.0d0 * Azb(nx,ny) - Azbb(nx,ny)
      
      end do
      end do

    end subroutine phia

END Module picexec
