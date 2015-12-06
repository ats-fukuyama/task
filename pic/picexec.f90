!./pic  ***** TASK/PIC EXEC *****

Module picexec

  PRIVATE
  PUBLIC pic_exec

CONTAINS

  SUBROUTINE pic_exec(iout)

    USE piccomm
    USE picsub,ONLY: poisson_f,poisson_m,efield,bfield,kine,pote
    USE piclib
    USE libmpi
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER,INTENT(OUT):: iout
    REAL(8)::abc,sum
    INTEGER:: nx,ny,nt,locv,npo,np
    INTEGER,DIMENSION(:),ALLOCATABLE:: locva
    REAL(8),DIMENSION(:,:),ALLOCATABLE:: suma

    ALLOCATE(locva(nxymax))
    ALLOCATE(suma(0:nxmax,0:nymax))

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
       call source(npmax,nxmax,nymax,xe,ye,rho,chrge,model_boundary)
       call source(npmax,nxmax,nymax,xi,yi,rho,chrgi,model_boundary)
       call boundary_rho(nxmax,nymax,rho,model_boundary)
       !..... sum charge densities over cores
       call mtx_allreduce_real8(rho,nxymax,3,suma,locva)
       DO ny=0,nymax
       DO nx=0,nxmax
          rho(nx,ny)=suma(nx,ny)
       END DO
       END DO

       !----- calculate electric field
       ipssn = 1
       IF(model_boundary.EQ.0) THEN
          call poisson_f(nxmax,nymax,nxmaxh1,nxmax1,nymax1, &
                         rho,phi,rhof,phif,awk,afwk,cform,ipssn)
       ELSE
          call poisson_m(nxmax1,nymax1,rho,phi,ipssn, &
                         model_matrix0,model_matrix1,model_matrix2, &
                         tolerance_matrix)
       END IF

       !----- current assignment
       jx(:,:)=0.d0
       jy(:,:)=0.d0
       jz(:,:)=0.d0


       call current(npmax,nxmax,nymax,xe,ye,vxe,vye,vze,chrge,jx,jy,jz, &
                    model_boundary)
       call current(npmax,nxmax,nymax,xi,yi,vxi,vyi,vzi,chrgi,jx,jy,jz, &
                    model_boundary)
       if (model_antenna .eq. 1) then
       call antenna(nxmax,nymax,jxant,jyant,jzant,phxant,phyant,phzant, &
            omega,time,jx,jy,jz)
       end if
       
       call boundary_j(nxmax,nymax,jx,jy,jz,model_boundary)
       
       !..... sum current densities over cores
       call mtx_allreduce_real8(jx,nxymax,3,suma,locva)
       DO ny=0,nymax
       DO nx=0,nxmax
          jx(nx,ny)=suma(nx,ny)
       END DO
       END DO
       call mtx_allreduce_real8(jy,nxymax,3,suma,locva)
       DO ny=0,nymax
       DO nx=0,nxmax
          jy(nx,ny)=suma(nx,ny)
       END DO
       END DO
       call mtx_allreduce_real8(jz,nxymax,3,suma,locva)
       DO ny=0,nymax
       DO nx=0,nxmax
          jz(nx,ny)=suma(nx,ny)
       END DO
       END DO

       !.......... calculate vector potential
       if(model_boundary .eq. 0) then
          call phia_periodic(nxmax,nymax,vcfact,dt,phi,phib,jx,jy,jz, &
                             Ax,Ay,Az,Axb,Ayb,Azb,Axbb,Aybb,Azbb)
       else
          call phia_reflective(nxmax,nymax,vcfact,dt,phi,phib,jx,jy,jz, &
                               Ax,Ay,Az,Axb,Ayb,Azb,Axbb,Aybb,Azbb, &
                               model_wg,xmin_wg,xmax_wg,ymin_wg,ymax_wg, &
                               amp_wg,ph_wg,rot_wg,eli_wg,omega,time,pi)
       endif
       !.......... calculate ex and ey and ez
       call efield(nxmax,nymax,dt,phi,Ax,Ay,Az,Axb,Ayb,Azb, &
            ex,ey,ez,esx,esy,esz,emx,emy,emz,model_push,model_boundary)

       !.......... calculate bxg and byg and bzg
       call bfield(nxmax,nymax,Ax,Ay,Az,Axb,Ayb,Azb, &
                               bx,by,bz,bxbg,bybg,bzbg,bb, &
                               model_push,model_boundary)
       if( mod(nt,ntgstep) .eq. 0 ) then
          call kine(npmax,vxe,vye,vze,akine1,me)
          call kine(npmax,vxi,vyi,vzi,akini1,mi)
          call pote(nxmax,nymax,ex,ey,ez,bx,by,bz,bxbg,bybg,bzbg,vcfact, &
                    apote,apotm)
          call mtx_allreduce1_real8(akine1,3,sum,locv)
          akine1=sum
          call mtx_allreduce1_real8(akini1,3,sum,locv)
          akini1=sum
          call mtx_allreduce1_real8(apote,3,sum,locv)
          apote=sum
          call mtx_allreduce1_real8(apotm,3,sum,locv)
          apotm=sum
       endif

       !..... push electrons
       call push(npmax,nxmax,nymax,xe,ye,ze,vxe,vye,vze, &
                 ex,ey,ez,bx,by,bz,dt,ctome,xeb,yeb,zeb, &
                 vparae,vperpe,model_boundary)
         
       !..... push ions
       call push(npmax,nxmax,nymax,xi,yi,zi,vxi,vyi,vzi, &
                 ex,ey,ez,bx,by,bz,dt,ctomi,xib,yib,zib, &
                 vparai,vperpi,model_boundary)

       !----- treat particles being out of the boundary

       if(model_boundary .eq. 0) then
          call bound_periodic(npmax,xe,ye,ze,x1,x2,y1,y2,z1,z2,alx,aly,alz)
          call bound_periodic(npmax,xi,yi,zi,x1,x2,y1,y2,z1,z2,alx,aly,alz)
       else
          call bound_reflective(npmax,xe,ye,ze,vxe,vye, &
                                x1,x2,y1,y2,z1,z2,alx,aly,alz)
          call bound_reflective(npmax,xi,yi,zi,vxi,vyi, &
                                x1,x2,y1,y2,z1,z2,alx,aly,alz)
       endif

       IF(npomax.GT.0) THEN
          ntocount=ntocount+1
          DO npo=1,npomax
             np=1+npostep*(npo-1)
             xpo(npo,ntocount)=xe(np)
             ypo(npo,ntocount)=ye(np)
             zpo(npo,ntocount)=ze(np)
             vxpo(npo,ntocount)=vxe(np)
             vypo(npo,ntocount)=vye(np)
             vzpo(npo,ntocount)=vze(np)
          END DO
       END IF



       !..... diagnostics to check energy conservation
       !.....            after pushing 
       if( mod(nt,ntgstep) .eq. 0 ) then
          call kine(npmax,vxe,vye,vze,akine2,me)
          call kine(npmax,vxi,vyi,vzi,akini2,mi)
          call mtx_allreduce1_real8(akine2,3,sum,locv) ! sum
          akine2=sum
          call mtx_allreduce1_real8(akini2,3,sum,locv) ! sum
          akini2=sum
          akine = 0.5d0 * ( akine1 + akine2 )
          akini = 0.5d0 * ( akini1 + akini2 )
          aktot = akine + akini
          aptot = apote + apotm
          atot  = aktot + aptot

          akine = akine - akine0
          akini = akini - akini0
          aktot = aktot - aktot0
          apote = apote - apote0
          apotm = apotm - apotm0
          aptot = aptot - aptot0
          atot  = atot  - atot0

          ntgcount=ntgcount+1
          timet(ntgcount)=time
          akinet(ntgcount)=akine
          akinit(ntgcount)=akini
          aktott(ntgcount)=aktot
          apotet(ntgcount)=apote
          apotmt(ntgcount)=apotm
          aptott(ntgcount)=aptot
          atott(ntgcount)=atot
       endif

       !..... evaluate profiles

       IF( MOD(nt,ntpstep) .EQ. 0 ) THEN
          ntpcount=ntpcount+1
          CALL profile(npmax,nxmax,nymax, &
                       xe,ye,vxe,vye,vze,vparae,vperpe,me, &
                       profilee(0:nxmax,0:nymax,1:9,ntpcount),model_boundary)
          CALL profile(npmax,nxmax,nymax, &
                       xi,yi,vxi,vyi,vzi,vparai,vperpi,mi, &
                       profilei(0:nxmax,0:nymax,1:9,ntpcount),model_boundary)
       END IF
 
       IF( nrank .eq. 0 ) THEN
          IF(MOD(ntcount,ntstep).EQ.0) THEN
             WRITE(6,'(I8,1PE12.4,I8,1P3E12.4)') &
                  ntcount,time,ntgcount,aktot,aptot,atot
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

    if( nrank .eq. 0 ) write(*,*) '*** wall clock time = ***', wtime

    DEALLOCATE(locva,suma)

    ntgmax=ntgcount
    ntomax=ntocount
    write(6,*) 'ntgmax,ntomax=',ntgmax,ntomax
    iout=1
  END SUBROUTINE pic_exec

!***********************************************************************
!      expand strage for history data
!
  SUBROUTINE pic_expand_storage
!***********************************************************************
    USE piccomm,ONLY: ntmax,ntgstep,ntgcount,ntgmax,ntpstep,ntpcount,ntpmax, &
                      npomax,ntomax,ntostep,ntocount,  &
                      timet,akinet,akinit,aktott,apotet,apotmt,aptott,atott, &
                      xpo,ypo,zpo,vxpo,vypo,vzpo, &
                      nxmax,nymax,timep,profilee,profilei
    IMPLICIT NONE
    REAL(8),DIMENSION(:,:),ALLOCATABLE:: work
    REAL(8),DIMENSION(:),ALLOCATABLE:: workt
    REAL(8),DIMENSION(:,:,:),ALLOCATABLE:: workp
    REAL(8),DIMENSION(:,:,:,:,:),ALLOCATABLE:: work_prof
    INTEGER:: ntgmax_old,ntomax_old,ntpmax_old,npomax_save
    REAL(8):: factor

    IF(ntgcount.eq.0) THEN
       IF(ALLOCATED(timet)) &
            DEALLOCATE(timet,akinet,akinit,aktott,apotet,apotmt,aptott,atott)
       ntgmax=MAX(ntmax/ntgstep,1)
       ALLOCATE(timet(ntgmax))
       ALLOCATE(akinet(ntgmax))
       ALLOCATE(akinit(ntgmax))
       ALLOCATE(aktott(ntgmax))
       ALLOCATE(apotet(ntgmax))
       ALLOCATE(apotmt(ntgmax))
       ALLOCATE(aptott(ntgmax))
       ALLOCATE(atott(ntgmax))
    ELSE
       ALLOCATE(work(ntgmax,8))
       work(1:ntgmax,1)=timet (1:ntgmax)
       work(1:ntgmax,2)=akinet(1:ntgmax)
       work(1:ntgmax,3)=akinit(1:ntgmax)
       work(1:ntgmax,4)=aktott(1:ntgmax)
       work(1:ntgmax,5)=apotet (1:ntgmax)
       work(1:ntgmax,6)=apotmt (1:ntgmax)
       work(1:ntgmax,7)=aptott (1:ntgmax)
       work(1:ntgmax,8)=atott (1:ntgmax)
       DEALLOCATE(timet,akinet,akinit,aktott,apotet,apotmt,aptott,atott)
       ntgmax_old=ntgmax
       ntgmax=ntgmax+ntmax/ntgstep
       ALLOCATE(timet(ntgmax))
       ALLOCATE(akinet(ntgmax))
       ALLOCATE(akinit(ntgmax))
       ALLOCATE(aktott(ntgmax))
       ALLOCATE(apotet(ntgmax))
       ALLOCATE(apotmt(ntgmax))
       ALLOCATE(aptott(ntgmax))
       ALLOCATE(atott(ntgmax))
       timet (1:ntgmax_old)=work(1:ntgmax_old,1)
       akinet(1:ntgmax_old)=work(1:ntgmax_old,2)
       akinit(1:ntgmax_old)=work(1:ntgmax_old,3)
       aktott(1:ntgmax_old)=work(1:ntgmax_old,4)
       apotet(1:ntgmax_old)=work(1:ntgmax_old,5)
       apotmt(1:ntgmax_old)=work(1:ntgmax_old,6)
       aptott(1:ntgmax_old)=work(1:ntgmax_old,7)
       atott (1:ntgmax_old)=work(1:ntgmax_old,8)
       DEALLOCATE(work)
    END IF

    IF(ntpcount.eq.0) THEN
       IF(ALLOCATED(timep)) &
            DEALLOCATE(timep,profilee,profilei)
       ntpmax=MAX(ntmax/ntpstep,1)
       ALLOCATE(timep(ntpmax))
       ALLOCATE(profilee(0:nxmax,0:nymax,9,ntpmax))
       ALLOCATE(profilei(0:nxmax,0:nymax,9,ntpmax))
    ELSE
       ALLOCATE(workt(ntpmax))
       ALLOCATE(work_prof(0:nxmax,0:nymax,9,ntpmax,2))
       workt(1:ntpmax)=timep(1:ntpmax)
       work_prof(0:nxmax,0:nymax,1:9,1:ntpmax,1) &
            =profilee(0:nxmax,0:nymax,1:9,1:ntpmax)
       work_prof(0:nxmax,0:nymax,1:9,1:ntpmax,2) &
            =profilei(0:nxmax,0:nymax,1:9,1:ntpmax)
       DEALLOCATE(timep,profilee,profilei)
       ntpmax_old=ntpmax
       ntpmax=ntpmax+ntmax/ntpstep
       ALLOCATE(timep(ntpmax))
       ALLOCATE(profilee(0:nxmax,0:nymax,9,ntpmax))
       ALLOCATE(profilei(0:nxmax,0:nymax,9,ntpmax))
       timep(1:ntpmax_old) = workt(1:ntpmax_old)
       profilee(0:nxmax,0:nymax,1:9,1:ntpmax_old) &
            =work_prof(0:nxmax,0:nymax,1:9,1:ntpmax_old,1)
       profilei(0:nxmax,0:nymax,1:9,1:ntpmax_old) &
            =work_prof(0:nxmax,0:nymax,1:9,1:ntpmax_old,2)
       DEALLOCATE(workt,work_prof)
    END IF

    IF(npomax.GT.0) THEN
       IF(ntocount.eq.0) THEN
          IF(ALLOCATED(xpo)) &
               DEALLOCATE(xpo,ypo,zpo,vxpo,vypo,vzpo)
          ntomax=MAX(ntmax/ntostep,1)
          ALLOCATE(xpo(npomax,ntomax))
          ALLOCATE(ypo(npomax,ntomax))
          ALLOCATE(zpo(npomax,ntomax))
          ALLOCATE(vxpo(npomax,ntomax))
          ALLOCATE(vypo(npomax,ntomax))
          ALLOCATE(vzpo(npomax,ntomax))
          npomax_save=npomax
       ELSE
          npomax=npomax_save
          ALLOCATE(workp(npomax,ntomax,6))
          workp(1:npomax,1:ntomax,1)=xpo(1:npomax,1:ntomax)
          workp(1:npomax,1:ntomax,2)=ypo(1:npomax,1:ntomax)
          workp(1:npomax,1:ntomax,3)=zpo(1:npomax,1:ntomax)
          workp(1:npomax,1:ntomax,4)=vxpo(1:npomax,1:ntomax)
          workp(1:npomax,1:ntomax,5)=vypo(1:npomax,1:ntomax)
          workp(1:npomax,1:ntomax,6)=vzpo(1:npomax,1:ntomax)
          DEALLOCATE(xpo,ypo,zpo,vxpo,vypo,vzpo)
          ntomax_old=ntomax
          ntomax=ntomax+ntmax/ntostep
          ALLOCATE(xpo(npomax,ntomax))
          ALLOCATE(ypo(npomax,ntomax))
          ALLOCATE(zpo(npomax,ntomax))
          ALLOCATE(vxpo(npomax,ntomax))
          ALLOCATE(vypo(npomax,ntomax))
          ALLOCATE(vzpo(npomax,ntomax))
          xpo(1:npomax,1:ntgmax_old)=workp(1:npomax,1:ntgmax_old,1)
          ypo(1:npomax,1:ntgmax_old)=workp(1:npomax,1:ntgmax_old,2)
          zpo(1:npomax,1:ntgmax_old)=workp(1:npomax,1:ntgmax_old,3)
          vxpo(1:npomax,1:ntgmax_old)=workp(1:npomax,1:ntgmax_old,4)
          vypo(1:npomax,1:ntgmax_old)=workp(1:npomax,1:ntgmax_old,5)
          vzpo(1:npomax,1:ntgmax_old)=workp(1:npomax,1:ntgmax_old,6)
          DEALLOCATE(workp)
       END IF
    END IF
  END SUBROUTINE pic_expand_storage
      
!***********************************************************************
  subroutine push(npmax,nxmax,nymax,x,y,z,vx,vy,vz,ex,ey,ez,bx,by,bz,dt,&
                  ctom,xb,yb,zb,vpara,vperp,model_boundary)
!***********************************************************************
    implicit none
    real(8), dimension(npmax) :: x, y, z, xb, yb, zb
    real(8), dimension(npmax) :: vx, vy, vz, vpara, vperp
    real(8), dimension(0:nxmax,0:nymax) :: ex, ey, ez, bx, by, bz
    real(8) :: ctom, dx, dy, dx1, dy1, dt, exx, eyy, ezz, bxx,&
               byy, bzz, vxn, vyn, vzn, vxzero, vyzero, vzzero, vxp, vyp,&
               vzp, sx1p, sx1m, sy1p, sy1m, sx2, sy2, sx2p, sx2m, sy2m, sy2p
    real(8) :: btot, vtot, bb2
    integer :: npmax, nxmax, nymax, model_boundary
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

       IF(model_boundary.EQ.0) THEN ! periodic
          if( nxp .eq. 0  ) nxpm = nxmax - 1
          if( nyp .eq. 0  ) nypm = nymax - 1
          !if( nxp .eq. nxmax-2) nxppp=0
          !if( nyp .eq. nymax-2) nyppp=0
          !if( nxp .eq. nxmax-1) nxpp =0
          !if( nyp .eq. nymax-1) nypp =0
          if( nxp .eq. nxmax-1) nxppp = 1
          if( nyp .eq. nymax-1) nyppp = 1
       ELSE   ! reflective: 
          if( nxp .eq. 0  ) nxpm = 0
          if( nyp .eq. 0  ) nypm = 0
          !if( nxp .eq. nxmax-2) nxppp=nxmax
          !if( nyp .eq. nymax-2) nyppp=nymax
          !if( nxp .eq. nxmax-1) nxpp =nxmax
          !if( nyp .eq. nymax-1) nypp =nymax
          if( nxp .eq. nxmax-1) nxppp = nxmax
          if( nyp .eq. nymax-1) nyppp = nymax
       END IF

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

       bb2 = bxx ** 2 + byy ** 2 + bzz ** 2
       vxn = vx(np) + 0.5D0 * ctom * exx * dt 
       vyn = vy(np) + 0.5D0 * ctom * eyy * dt
       vzn = vz(np) + 0.5D0 * ctom * ezz * dt

       vxzero = vxn + 0.5D0 * ctom * (vyn * bzz - vzn * byy) * dt
       vyzero = vyn + 0.5D0 * ctom * (vzn * bxx - vxn * bzz) * dt 
       vzzero = vzn + 0.5D0 * ctom * (vxn * byy - vyn * bxx) * dt

       vxp = vxn + 1.0d0/(1.0d0 + 0.25d0 * (ctom * dt) ** 2 * bb2) & 
                   * ctom * (vyzero * bzz - vzzero * byy) * dt
       vyp = vyn + 1.0d0/(1.0d0 + 0.25d0 * (ctom * dt) ** 2 * bb2) &
                   * ctom * (vzzero * bxx - vxzero * bzz) * dt 
       vzp = vzn + 1.0d0/(1.0d0 + 0.25d0 * (ctom * dt) ** 2 * bb2) & 
                   * ctom * (vxzero * byy - vyzero * bxx) * dt
 
       vx(np) = vxp + 0.5D0 * ctom * exx * dt
       vy(np) = vyp + 0.5D0 * ctom * eyy * dt
       vz(np) = vzp + 0.5D0 * ctom * ezz * dt
         
       xb(np) = x(np)
       yb(np) = y(np)
       zb(np) = z(np)
       x(np) = x(np) + vx(np) * dt
       y(np) = y(np) + vy(np) * dt
       z(np) = z(np) + vz(np) * dt

       btot=SQRT(bb2)
       IF(btot.EQ.0.D0) THEN
          vpara(np)=vx(np)
          vperp(np)=SQRT(vy(np)**2+vz(np)**2)
       ELSE
          vtot=SQRT(vx(np)**2+vy(np)**2+vz(np)**2)
          vpara(np)=(bxx*vx(np)+byy*vy(np)+bzz*vz(np))/btot
          vperp(np)=SQRT(vtot**2-vpara(np)**2)
       END IF

    end do

  end subroutine push

!***********************************************************************
    subroutine bound_periodic(npmax,x,y,z,x1,x2,y1,y2,z1,z2,alx,aly,alz)
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

    end subroutine bound_periodic

!***********************************************************************
    subroutine bound_reflective(npmax,x,y,z,vx,vy, &
                                x1,x2,y1,y2,z1,z2,alx,aly,alz)
!***********************************************************************
      implicit none
      real(8), dimension(npmax) :: x, y, z, vx, vy
      real(8) :: x1, x2, y1, y2, z1, z2, alx, aly, alz
      integer :: npmax, np

      do np = 1, npmax
         if( x(np) .lt. x1 ) then
            x(np) = -x(np)
            vx(np) = -vx(np)
         elseif( x(np) .gt. x2 ) then
            x(np) = alx - (x(np) - alx)
            vx(np) = -vx(np)
         endif
 
         if( y(np) .lt. y1 ) then
            y(np) = -y(np)
            vy(np) = -vy(np)
         elseif( y(np) .gt. y2 ) then
            y(np) = aly - (y(np) - aly)
            vy(np) = -vy(np)
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

    end subroutine bound_reflective

!***********************************************************************
    subroutine source(npmax,nxmax,nymax,x,y,rho,chrg,model_boundary)
!***********************************************************************

      implicit none
      real(8), dimension(npmax)        :: x, y
      real(8), dimension(0:nxmax,0:nymax) :: rho
      real(8) :: chrg, factor, dx, dy, dz, sx2, sy2, sx2p, sy2p, sx2m, sy2m,&
           dx1,dy1
      integer :: npmax, nxmax, nymax, model_boundary
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

         IF(model_boundary.EQ.0) THEN ! periodic
            if( nxp .eq. 0  ) nxpm = nxmax - 1
            if( nyp .eq. 0  ) nypm = nymax - 1
            !if( nxp .eq. nxmax-2) nxppp=0
            !if( nyp .eq. nymax-2) nyppp=0
            !if( nxp .eq. nxmax-1) nxpp =0
            !if( nyp .eq. nymax-1) nypp =0
            if( nxp .eq. nxmax-1) nxppp=1
            if( nyp .eq. nymax-1) nyppp=1
         ELSE   ! reflective: 
            if( nxp .eq. 0  ) nxpm = 0
            if( nyp .eq. 0  ) nypm = 0
            !if( nxp .eq. nxmax-2) nxppp=nxmax
            !if( nyp .eq. nymax-2) nyppp=nymax
            !if( nxp .eq. nxmax-1) nxpp =nxmax
            !if( nyp .eq. nymax-1) nypp =nymax
            if( nxp .eq. nxmax-1) nxppp = nxmax
            if( nyp .eq. nymax-1) nyppp = nymax
         END IF

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
    subroutine boundary_rho(nxmax,nymax,rho,model_boundary)
!***********************************************************************
      implicit none
      real(8), dimension(0:nxmax,0:nymax) :: rho
      integer:: nxmax,nymax,model_boundary
      integer:: nx,ny

      !!..... set charge densities at the boundary

      IF(model_boundary.EQ.0) THEN  ! periodic
         do ny = 0, nymax
            rho(0,ny) = rho(0,ny) + rho(nxmax,ny)
            rho(nxmax,ny) = rho(0,ny)
         end do
         do nx = 0, nxmax
            rho(nx,0) = rho(nx,0) + rho(nx,nymax)
            rho(nx,nymax) = rho(nx,0)
         end do
      ELSE                         ! reflecting
         do ny = 0, nymax
            rho(0,ny)     = 2.D0 * rho(0,ny)
            rho(nxmax,ny) = 2.D0 * rho(nxmax,ny)
         end do
         do nx = 1, nxmax-1
            rho(nx,0)     = 2.D0 * rho(nx,0)
            rho(nx,nymax) = 2.D0 * rho(nx,nymax)
         end do
      END IF

    end subroutine boundary_rho

!***********************************************************************
    subroutine current(npmax,nxmax,nymax,x,y,vx,vy,vz,chrg,jx,jy,jz, &
                       model_boundary)
!***********************************************************************
    implicit none

    real(8), dimension(npmax) :: x, y, vx, vy, vz 
    real(8), dimension(0:nxmax,0:nymax) :: jx, jy, jz
    real(8) :: chrg, dt, dx, dy, dx1, dy1, &
         sx1p, sy1p, sx1m, sy1m, sx2, sy2, sx2p, sy2p, sx2m, sy2m, factor
    integer :: npmax, nxmax, nymax, model_boundary
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

         IF(model_boundary.EQ.0) THEN ! periodic
            if( nxp .eq. 0  ) nxpm = nxmax - 1
            if( nyp .eq. 0  ) nypm = nymax - 1
            !if( nxp .eq. nxmax-2) nxppp=0
            !if( nyp .eq. nymax-2) nyppp=0
            !if( nxp .eq. nxmax-1) nxpp =0
            !if( nyp .eq. nymax-1) nypp =0
            if( nxp .eq. nxmax-1) nxppp=1
            if( nyp .eq. nymax-1) nyppp=1
         ELSE   ! reflective: 
            if( nxp .eq. 0  ) nxpm = 0
            if( nyp .eq. 0  ) nypm = 0
            !if( nxp .eq. nxmax-2) nxppp=nxmax
            !if( nyp .eq. nymax-2) nyppp=nymax
            !if( nxp .eq. nxmax-1) nxpp =nxmax
            !if( nyp .eq. nymax-1) nypp =nymax
            if( nxp .eq. nxmax-1) nxppp = nxmax
            if( nyp .eq. nymax-1) nyppp = nymax
         END IF

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

         do ny=5,10
            jy(3,ny) = jy(3,ny) + 0.5d0 * jyt
            jy(4,ny) = jy(4,ny) + 0.5d0 * jyt
            jz(3,ny) = jz(3,ny) + 0.5d0 * jzt
            jz(4,ny) = jz(4,ny) + 0.5d0 * jzt
         end do

   end subroutine antenna

!***********************************************************************
    subroutine boundary_j(nxmax,nymax,jx,jy,jz,model_boundary)
!***********************************************************************
    implicit none

    integer :: nxmax,nymax,model_boundary
    real(8), dimension(0:nxmax,0:nymax) :: jx, jy, jz
    integer :: nx,ny

    IF(model_boundary.EQ.0) THEN     ! periodic
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
      ELSE                           ! reflective
        do ny=0,nymax
           jx(0,ny) = 2.0d0 * jx(0,ny)
           jy(0,ny) = 2.0d0 * jy(0,ny)
           jz(0,ny) = 2.0d0 * jz(0,ny)
           jx(nxmax,ny) = 2.0d0 * jx(nxmax,ny)
           jy(nxmax,ny) = 2.0d0 * jy(nxmax,ny)
           jz(nxmax,ny) = 2.0d0 * jz(nxmax,ny)
        end do
        do nx=1,nxmax-1
           jx(nx,0) = 2.0d0 * jx(nx,0)
           jy(nx,0) = 2.0d0 * jy(nx,0)
           jz(nx,0) = 2.0d0 * jz(nx,0)
           jx(nx,nymax) = 2.0d0 * jx(nx,nymax)
           jy(nx,nymax) = 2.0d0 * jy(nx,nymax)
           jz(nx,nymax) = 2.0d0 * jz(nx,nymax)
        end do
     END IF

   end subroutine boundary_j

!***********************************************************************
    subroutine phia_periodic(nxmax,nymax,vcfact,dt,phi,phib,jx,jy,jz, &
                             Ax,Ay,Az,Axb,Ayb,Azb,Axbb,Aybb,Azbb)
!***********************************************************************
   !original subroutine
      implicit none
      real(8), dimension(0:nxmax,0:nymax) :: phi,phib,jx,jy,jz,Ax,Ay,Az, &
                                             Axb,Ayb,Azb,Axbb,Aybb,Azbb
      integer :: nxmax, nymax, nx, ny, nxm, nxp, nyp, nym, nypm
      real(8) :: vcfact, dt

 ! Solution of maxwell equation in the A-phi formulation by difference method
 ! vcfact is the ratio of the light speed to lattice parameter times 
 ! plasma frequency

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

    end subroutine phia_periodic

!***********************************************************************
    subroutine phia_reflective(nxmax,nymax,vcfact,dt,phi,phib,jx,jy,jz, &
                               Ax,Ay,Az,Axb,Ayb,Azb,Axbb,Aybb,Azbb, &
                               model_wg,xmin_wg,xmax_wg,ymin_wg,ymax_wg, &
                               amp_wg,ph_wg,rot_wg,eli_wg,omega,time,pi)
!***********************************************************************
   !original subroutine
      implicit none
      real(8), dimension(0:nxmax,0:nymax) :: phi,phib,jx,jy,jz,Ax,Ay,Az, &
                                             Axb,Ayb,Azb,Axbb,Aybb,Azbb
      integer :: nxmax, nymax, nx, ny, nxm, nxp, nyp, nym, nypm
      integer :: model_wg
      real(8) :: xmin_wg,xmax_wg,ymin_wg,ymax_wg,amp_wg,ph_wg,rot_wg,eli_wg
      real(8) :: omega,time,pi
      real(8) :: vcfact,dt,dph,x,y

 ! Solution of maxwell equation in the A-phi formulation by difference method
 ! vcfact is the ratio of the light speed to lattice parameter times plasma 
 ! frequency

      do nx = 1, nxmax-1
      do ny = 1, nymax-1

         nxm = nx - 1
         nxp = nx + 1
         nym = ny - 1
         nyp = ny + 1

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
      !boundary condition for reflection

      Ax(0,:)=0.d0
      Ay(0,:)=0.d0
      Az(0,:)=0.d0
      Ax(nxmax,:)=0.d0
      Ay(nxmax,:)=0.d0
      Az(nxmax,:)=0.d0
      Ax(:,0)=0.d0
      Ay(:,0)=0.d0
      Az(:,0)=0.d0
      Ax(:,nymax)=0.d0
      Ay(:,nymax)=0.d0
      Az(:,nymax)=0.d0

      SELECT CASE(model_wg)
      CASE(0)
         dph=ph_wg/(ymax_wg-ymin_wg)
         DO ny=1,nymax
            y=dble(ny)
            IF(y.GE.ymin_wg.AND.y.LE.ymax_wg) &
               Ay(0,ny)=amp_wg*sin(omega*time-2.D0*pi*dph*(y-ymin_wg))
         END DO
      END SELECT

    end subroutine phia_reflective

!***********************************************************************
    subroutine profile(npmax,nxmax,nymax,x,y,vx,vy,vz,vpara,vperp,mass, &
                       profiles,model_boundary)
!***********************************************************************
      IMPLICIT NONE
      INTEGER:: npmax,nxmax,nymax,model_boundary
      REAL(8):: mass
      REAL(8):: x(npmax),y(npmax),vx(npmax),vy(npmax),vz(npmax)
      REAL(8):: vpara(npmax),vperp(npmax)
      REAL(8):: profiles(0:nxmax,0:nymax,9)
      INTEGER:: np,nx,ny
      REAL(8):: xp,yp,fp(9),factor

      profiles(0:nxmax,0:nymax,1:9)=0.D0
      factor=DBLE(nxmax)*DBLE(nymax)/DBLE(npmax)

      DO np=1,npmax
         xp=x(np)
         yp=y(np)
         fp(1)=1.D0
         fp(2)=vx(np)
         fp(3)=vy(np)
         fp(4)=vz(np)
         fp(5)=vpara(np)
         fp(6)=vperp(np)
         fp(7)=0.5D0*mass*vpara(np)**2
         fp(8)=0.5D0*mass*vperp(np)**2
         fp(9)=0.5D0*mass*(vx(np)**2+vy(np)**2+vz(np)**2)
         CALL profile_sub(nxmax,nymax,xp,yp,fp,profiles,9,model_boundary)
      END DO

      CALL profile_boundary(nxmax,nymax,profiles,9,model_boundary)

      profiles(0:nxmax,0:nymax,1:9)=factor*profiles(0:nxmax,0:nymax,1:9)

      DO ny=0,nymax
         DO nx=0,nxmax
            profiles(nx,ny,2)=profiles(nx,ny,2)/profiles(nx,ny,1)
            profiles(nx,ny,3)=profiles(nx,ny,3)/profiles(nx,ny,1)
            profiles(nx,ny,4)=profiles(nx,ny,4)/profiles(nx,ny,1)
            profiles(nx,ny,5)=profiles(nx,ny,5)/profiles(nx,ny,1)
            profiles(nx,ny,6)=profiles(nx,ny,6)/profiles(nx,ny,1)
            profiles(nx,ny,7)=profiles(nx,ny,7)/profiles(nx,ny,1)
            profiles(nx,ny,8)=profiles(nx,ny,8)/profiles(nx,ny,1)
            profiles(nx,ny,9)=profiles(nx,ny,9)/profiles(nx,ny,1)
         END DO
      END DO
    END subroutine profile

!***********************************************************************
    subroutine profile_sub(nxmax,nymax,xp,yp,fp,fxy,imax,model_boundary)
!***********************************************************************

      implicit none
      real(8) :: xp, yp, fp(imax)
      real(8), dimension(0:nxmax,0:nymax,imax) :: fxy
      real(8) :: dx, dy, dz, sx2, sy2, sx2p, sy2p, sx2m, sy2m, dx1,dy1
      integer :: nxmax, nymax, imax, model_boundary
      integer:: nx, ny, nxp, nyp, nxpp, nxpm, nypp, nypm, nxppp, nyppp, i

      nxp = xp
      nyp = yp

      dx  = xp - dble(nxp)
      dy  = yp - dble(nyp)

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

      IF(model_boundary.EQ.0) THEN ! periodic
         if( nxp .eq. 0  ) nxpm = nxmax - 1
         if( nyp .eq. 0  ) nypm = nymax - 1
         if( nxp .eq. nxmax-1) nxppp=1
         if( nyp .eq. nymax-1) nyppp=1
      ELSE   ! reflective: 
         if( nxp .eq. 0  ) nxpm = 0
         if( nyp .eq. 0  ) nypm = 0
         if( nxp .eq. nxmax-1) nxppp = nxmax
         if( nyp .eq. nymax-1) nyppp = nymax
      END IF

      if(dx .le. 0.5d0 .and. dy .le. 0.5d0) then
         DO i=1,imax
            fxy(nxpm,nypm,i) = fxy(nxpm,nypm,i) + sx2m * sy2m * fp(i)
            fxy(nxpm,nyp ,i) = fxy(nxpm,nyp ,i) + sx2m * sy2  * fp(i)
            fxy(nxpm,nypp,i) = fxy(nxpm,nypp,i) + sx2m * sy2p * fp(i)
            fxy(nxp ,nypm,i) = fxy(nxp ,nypm,i) + sx2  * sy2m * fp(i)
            fxy(nxp ,nyp ,i) = fxy(nxp ,nyp ,i) + sx2  * sy2  * fp(i)
            fxy(nxp ,nypp,i) = fxy(nxp ,nypp,i) + sx2  * sy2p * fp(i)
            fxy(nxpp,nypm,i) = fxy(nxpp,nypm,i) + sx2p * sy2m * fp(i)
            fxy(nxpp,nyp ,i) = fxy(nxpp,nyp ,i) + sx2p * sy2  * fp(i)
            fxy(nxpp,nypp,i) = fxy(nxpp,nypp,i) + sx2p * sy2p * fp(i)
         END DO
      else if(dx .le. 0.5d0 .and. dy .ge. 0.5d0) then
         DO i=1,imax
            fxy(nxpm ,nyp  ,i) = fxy(nxpm ,nyp  ,i) + sx2m * sy2m * fp(i)
            fxy(nxpm ,nypp ,i) = fxy(nxpm ,nypp ,i) + sx2m * sy2  * fp(i)
            fxy(nxpm ,nyppp,i) = fxy(nxpm ,nyppp,i) + sx2m * sy2p * fp(i)
            fxy(nxp  ,nyp  ,i) = fxy(nxp  ,nyp  ,i) + sx2  * sy2m * fp(i)
            fxy(nxp  ,nypp ,i) = fxy(nxp  ,nypp ,i) + sx2  * sy2  * fp(i)
            fxy(nxp  ,nyppp,i) = fxy(nxp  ,nyppp,i) + sx2  * sy2p * fp(i)
            fxy(nxpp ,nyp  ,i) = fxy(nxpp ,nyp  ,i) + sx2p * sy2m * fp(i)
            fxy(nxpp ,nypp ,i) = fxy(nxpp ,nypp ,i) + sx2p * sy2  * fp(i)
            fxy(nxpp ,nyppp,i) = fxy(nxpp ,nyppp,i) + sx2p * sy2p * fp(i)
         END DO
      else if(dx .ge. 0.5d0 .and. dy .le. 0.5d0) then
         DO i=1,imax
            fxy(nxp  ,nypm,i) = fxy(nxp  ,nypm ,i) + sx2m * sy2m * fp(i)
            fxy(nxp  ,nyp ,i) = fxy(nxp  ,nyp  ,i) + sx2m * sy2  * fp(i)
            fxy(nxp  ,nypp,i) = fxy(nxp  ,nypp ,i) + sx2m * sy2p * fp(i)
            fxy(nxpp ,nypm,i) = fxy(nxpp ,nypm ,i) + sx2  * sy2m * fp(i)
            fxy(nxpp ,nyp ,i) = fxy(nxpp ,nyp  ,i) + sx2  * sy2  * fp(i)
            fxy(nxpp ,nypp,i) = fxy(nxpp ,nypp ,i) + sx2  * sy2p * fp(i)
            fxy(nxppp,nypm,i) = fxy(nxppp,nypm ,i) + sx2p * sy2m * fp(i)
            fxy(nxppp,nyp ,i) = fxy(nxppp,nyp  ,i) + sx2p * sy2  * fp(i)
            fxy(nxppp,nypp,i) = fxy(nxppp,nypp ,i) + sx2p * sy2p * fp(i)
         END DO
      else 
         DO i=1,imax
            fxy(nxp  ,nyp  ,i) = fxy(nxp  ,nyp  ,i) + sx2m * sy2m * fp(i)
            fxy(nxp  ,nypp ,i) = fxy(nxp  ,nypp ,i) + sx2m * sy2  * fp(i)
            fxy(nxp  ,nyppp,i) = fxy(nxp  ,nyppp,i) + sx2m * sy2p * fp(i)
            fxy(nxpp ,nyp  ,i) = fxy(nxpp ,nyp  ,i) + sx2  * sy2m * fp(i)
            fxy(nxpp ,nypp ,i) = fxy(nxpp ,nypp ,i) + sx2  * sy2  * fp(i)
            fxy(nxpp ,nyppp,i) = fxy(nxpp ,nyppp,i) + sx2  * sy2p * fp(i)
            fxy(nxppp,nyp  ,i) = fxy(nxppp,nyp  ,i) + sx2p * sy2m * fp(i)
            fxy(nxppp,nypp ,i) = fxy(nxppp,nypp ,i) + sx2p * sy2  * fp(i)
            fxy(nxppp,nyppp,i) = fxy(nxppp,nyppp,i) + sx2p * sy2p * fp(i)
         END DO
      endif
    END subroutine profile_sub

!***********************************************************************
    subroutine profile_boundary(nxmax,nymax,fxy,imax,model_boundary)
!***********************************************************************

      implicit none
      INTEGER:: nxmax,nymax,imax,model_boundary
      real(8), dimension(0:nxmax,0:nymax,imax) :: fxy
      integer:: nx, ny, i

      IF(model_boundary.EQ.0) THEN  ! periodic
         DO i=1,imax
            DO ny = 0, nymax
               fxy(0,ny,i) = fxy(0,ny,i) + fxy(nxmax,ny,i)
               fxy(nxmax,ny,i) = fxy(0,ny,i)
            END DO
         END DO
         DO i=1,imax
            DO nx = 0, nxmax
               fxy(nx,0,i) = fxy(nx,0,i) + fxy(nx,nymax,i)
               fxy(nx,nymax,i) = fxy(nx,0,i)
            END DO
         END DO
      ELSE                         ! reflecting
         DO i=1,imax
            DO ny = 0, nymax
               fxy(0,ny,i)     = 2.D0 * fxy(0,ny,i)
               fxy(nxmax,ny,i) = 2.D0 * fxy(nxmax,ny,i)
            END DO
         END DO
         DO nx = 1, nxmax-1
            DO i=1,imax
               fxy(nx,0,i)     = 2.D0 * fxy(nx,0,i)
               fxy(nx,nymax,i) = 2.D0 * fxy(nx,nymax,i)
            END DO
         END DO
      END IF

    END subroutine profile_boundary
    
END Module picexec
