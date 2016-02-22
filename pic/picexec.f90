!./pic  ***** TASK/PIC EXEC *****

MODULE picexec

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
    REAL(8)::sum
    INTEGER:: nx,ny,nt,locv,npo,np
    INTEGER,DIMENSION(:),ALLOCATABLE:: locva
    REAL(8),DIMENSION(:,:),ALLOCATABLE:: suma

    ALLOCATE(locva(nxymax))
    ALLOCATE(suma(0:nxmax,0:nymax))

    !----- start parallel calculation

    !-----------------------------------------------------------------------
    !----- allocate time hisotry data -------------------------------------------
    !-----------------------------------------------------------------------

    CALL pic_expand_storage

    !-----------------------------------------------------------------------
    !----- start of main do-loop -------------------------------------------
    !-----------------------------------------------------------------------
    DO nt = 1, ntmax
       time = time + dt
       ntcount = ntcount + 1
       DO nx = 0, nxmax
          DO ny = 0, nymax
             Axbb(nx,ny) = Axb(nx,ny)
             Aybb(nx,ny) = Ayb(nx,ny)
             Azbb(nx,ny) = Azb(nx,ny)
             Axb(nx,ny)  = Ax(nx,ny)
             Ayb(nx,ny)  = Ay(nx,ny)
             Azb(nx,ny)  = Az(nx,ny)
             phib(nx,ny) = phi(nx,ny)
          END DO
       END DO
       !----- charge assignment
       rho(:,:)=0.0d0
       CALL source(npmax,nxmax,nymax,xe,ye,rho,chrge,model_boundary)
       CALL source(npmax,nxmax,nymax,xi,yi,rho,chrgi,model_boundary)
       CALL boundary_rho(nxmax,nymax,rho,model_boundary)
       !..... sum charge densities over cores
       CALL mtx_allreduce_real8(rho,nxymax,3,suma,locva)
       DO ny=0,nymax
          DO nx=0,nxmax
             rho(nx,ny)=suma(nx,ny)/dble(nsize)
          END DO
       END DO
       !----- calculate electric field
       ipssn = 1
       IF(model_boundary.EQ.0) THEN
          CALL poisson_f(nxmax,nymax,nxmaxh1,nxmax1,nymax1, &
               rho,phi,rhof,phif,awk,afwk,cform,ipssn)
       ELSE
          CALL poisson_m(nxmax1,nymax1,rho,phi,ipssn, &
               model_matrix0,model_matrix1,model_matrix2, &
               tolerance_matrix,model_boundary,dlen)
       END IF
       !----- current assignment
       jx(:,:)=0.d0
       jy(:,:)=0.d0
       jz(:,:)=0.d0
       CALL current(npmax,nxmax,nymax,xe,ye,vxe,vye,vze,chrge,jx,jy,jz, &
            model_boundary)
       CALL current(npmax,nxmax,nymax,xi,yi,vxi,vyi,vzi,chrgi,jx,jy,jz, &
            model_boundary)
       IF (model_antenna .EQ. 1) THEN
          CALL antenna(nxmax,nymax,jxant,jyant,jzant,phxant,phyant,phzant, &
               omega,time,jx,jy,jz)
       END IF

       CALL boundary_j(nxmax,nymax,jx,jy,jz,model_boundary)
       !..... sum current densities over cores
       CALL mtx_allreduce_real8(jx,nxymax,3,suma,locva)
       DO ny=0,nymax
          DO nx=0,nxmax
             jx(nx,ny)=suma(nx,ny)/dble(nsize)
          END DO
       END DO
       CALL mtx_allreduce_real8(jy,nxymax,3,suma,locva)
       DO ny=0,nymax
          DO nx=0,nxmax
             jy(nx,ny)=suma(nx,ny)/dble(nsize)
          END DO
       END DO
       CALL mtx_allreduce_real8(jz,nxymax,3,suma,locva)
       DO ny=0,nymax
          DO nx=0,nxmax
             jz(nx,ny)=suma(nx,ny)/dble(nsize)
          END DO
       END DO

       !.......... calculate vector potential
       IF(model_boundary .EQ. 0) THEN
          CALL phia_periodic(nxmax,nymax,vcfact,dt,phi,phib,jx,jy,jz, &
               Ax,Ay,Az,Axb,Ayb,Azb,Axbb,Aybb,Azbb)
       ELSE
          CALL phia_reflective(nxmax,nymax,vcfact,dt,phi,phib,jx,jy,jz, &
               Ax,Ay,Az,Axb,Ayb,Azb,Axbb,Aybb,Azbb, &
               model_wg,xmin_wg,xmax_wg,ymin_wg,ymax_wg, &
               amp_wg,ph_wg,rot_wg,eli_wg,omega,time,pi,&
               model_boundary,dlen)
       ENDIF
       !.......... calculate ex and ey and ez
       CALL efield(nxmax,nymax,dt,phi,Ax,Ay,Az,Axb,Ayb,Azb, &
            ex,ey,ez,esx,esy,esz,emx,emy,emz,model_push,model_boundary)

       !.......... calculate bxg and byg and bzg
       CALL bfield(nxmax,nymax,Ax,Ay,Az,Axb,Ayb,Azb, &
            bx,by,bz,bxbg,bybg,bzbg,bb, model_push,model_boundary)
       IF( MOD(nt,ntgstep) .EQ. 0 ) THEN
          CALL kine(npmax,vxe,vye,vze,akine1,me)
          CALL kine(npmax,vxi,vyi,vzi,akini1,mi)
          CALL pote(nxmax,nymax,ex,ey,ez,bx,by,bz,bxbg,bybg,bzbg,vcfact, &
               apote,apotm)
          CALL mtx_allreduce1_real8(akine1,3,sum,locv)
          akine1=sum/dble(nsize)
          CALL mtx_allreduce1_real8(akini1,3,sum,locv)
          akini1=sum/dble(nsize)
          CALL mtx_allreduce1_real8(apote,3,sum,locv)
          apote=sum/dble(nsize)
          CALL mtx_allreduce1_real8(apotm,3,sum,locv)
          apotm=sum/dble(nsize)
       ENDIF
       !..... push electrons
       CALL push(npmax,nxmax,nymax,xe,ye,ze,vxe,vye,vze, &
            ex,ey,ez,bx,by,bz,dt,ctome,xeb,yeb,zeb, &
            vparae,vperpe,model_boundary,vcfact)

       !..... push ions
       CALL push(npmax,nxmax,nymax,xi,yi,zi,vxi,vyi,vzi, &
            ex,ey,ez,bx,by,bz,dt,ctomi,xib,yib,zib, &
            vparai,vperpi,model_boundary,vcfact)

       !----- treat particles being out of the boundary

       IF(model_boundary .EQ. 0) THEN
          CALL bound_periodic(npmax,xe,ye,ze,x1,x2,y1,y2,z1,z2,alx,aly,alz)
          CALL bound_periodic(npmax,xi,yi,zi,x1,x2,y1,y2,z1,z2,alx,aly,alz)
       ELSE
          CALL bound_reflective(npmax,xe,ye,ze,vxe,vye, &
               x1,x2,y1,y2,z1,z2,alx,aly,alz)
          CALL bound_reflective(npmax,xi,yi,zi,vxi,vyi, &
               x1,x2,y1,y2,z1,z2,alx,aly,alz)
       ENDIF

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
       IF( MOD(nt,ntgstep) .EQ. 0 ) THEN
          CALL kine(npmax,vxe,vye,vze,akine2,me)
          CALL kine(npmax,vxi,vyi,vzi,akini2,mi)
          CALL mtx_allreduce1_real8(akine2,3,sum,locv) ! sum
          akine2=sum/dble(nsize)
          CALL mtx_allreduce1_real8(akini2,3,sum,locv) ! sum
          akini2=sum/dble(nsize)
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
       ENDIF

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

       IF( nrank .EQ. 0 ) THEN
          IF(MOD(ntcount,ntstep).EQ.0) THEN
             WRITE(6,'(I8,1PE12.4,I8,1P3E12.4)') &
                  ntcount,time,ntgcount,aktot,aptot,atot
          END IF
       END IF
    END DO
    !-----------------------------------------------------------------------
    !----- end of main do-loop ---------------------------------------------
    !-----------------------------------------------------------------------

    !..... output wall clock time
    CALL mtx_barrier
    wtime2 = mpi_wtime()
    wtime  = wtime2 - wtime1

    IF( nrank .EQ. 0 ) WRITE(*,*) '*** wall clock time = ***', wtime

    DEALLOCATE(locva,suma)

    ntgmax=ntgcount
    ntomax=ntocount
    WRITE(6,*) 'ntgmax,ntomax=',ntgmax,ntomax
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

    IF(ntgcount.EQ.0) THEN
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

    IF(ntpcount.EQ.0) THEN
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
       IF(ntocount.EQ.0) THEN
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
  SUBROUTINE push(npmax,nxmax,nymax,x,y,z,vx,vy,vz,ex,ey,ez,bx,by,bz,dt,&
       ctom,xb,yb,zb,vpara,vperp,model_boundary,vcfact)
    !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(npmax) :: x, y, z, xb, yb, zb
    REAL(8), DIMENSION(npmax) :: vx, vy, vz, vpara, vperp
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: ex, ey, ez, bx, by, bz
    REAL(8) :: ctom, dx, dy, dx1, dy1, dt, exx, eyy, ezz, bxx,&
         byy, bzz, vxm, vym, vzm, vxzero, vyzero, vzzero, vxp, vyp,&
         vzp, sx1p, sx1m, sy1p, sy1m, sx2, sy2, sx2p, sx2m, sy2m, sy2p
    REAL(8) :: btot, vtot, bb2 ,vcfact, gamma
    INTEGER :: npmax, nxmax, nymax, model_boundary
    INTEGER :: np, nx, ny, nxp, nyp, nxpp, nxpm, nypp, nypm, nxppp, nyppp

    DO np = 1, npmax

       ! calculate the electric field at the particle position

       nxp = x(np)
       nyp = y(np)
       dx = x(np) - DBLE(nxp)
       dy = y(np) - DBLE(nyp)
       dx1 = 1.0d0 - dx
       dy1 = 1.0d0 - dy

       IF(dx .LE. 0.5d0) THEN
          sx2  = 3.0d0/4 - dx ** 2
          sx2p = 1.0d0/2 * (1.0d0/2 + dx) ** 2
          sx2m = 1.0d0/2 * (1.0d0/2 - dx) ** 2
       ELSE
          sx2  = 3.0d0/4 - (dx - 1.0d0) ** 2
          sx2p = 1.0d0/2 * (-1.0d0/2 + dx) ** 2
          sx2m = 1.0d0/2 * (3.0d0/2 - dx) ** 2
       ENDIF
       IF(dy .LE. 0.5d0) THEN
          sy2  = 3.0d0/4 - dy ** 2
          sy2p = 1.0d0/2 * (1.0d0/2 + dy) ** 2
          sy2m = 1.0d0/2 * (1.0d0/2 - dy) ** 2
       ELSE
          sy2  = 3.0d0/4 - (dy - 1.0d0) ** 2
          sy2p = 1.0d0/2 * (-1.0d0/2 + dy) ** 2
          sy2m = 1.0d0/2 * (3.0d0/2 - dy) ** 2
       ENDIF

       nxpm  = nxp - 1
       nxpp  = nxp + 1
       nxppp = nxp + 2
       nypm  = nyp - 1
       nypp  = nyp + 1
       nyppp = nyp + 2

       IF(model_boundary.EQ.0) THEN ! periodic
          IF( nxp .EQ. 0  ) nxpm = nxmax-1
          IF( nyp .EQ. 0  ) nypm = nymax-1
          IF( nxp .EQ. nxmax-1) nxppp = 1
          IF( nyp .EQ. nymax-1) nyppp = 1
       ELSE   ! reflective:
          IF( nxp .EQ. 0  ) nxpm = 0
          IF( nyp .EQ. 0  ) nypm = 0
          IF( nxp .EQ. nxmax-1) nxppp = nxmax
          IF( nyp .EQ. nymax-1) nyppp = nymax
       END IF

       ! electric field and magnetic field
       IF(dx .LE. 0.5d0 .AND. dy .LE. 0.5d0) THEN
        !  IF(model_boundary .NE. 0 .AND. nxp .EQ. 0)  THEN
        !    ey(nxpm,nypp) = - ey(nxpp,nypp)
        !    ey(nxpm,nyp ) = - ey(nxpp,nyp )
        !    ez(nxpm,nypp) = - ez(nxpp,nypp)
        !    ez(nxpm,nyp ) = - ez(nxpp,nyp )
        !    ez(nxpm,nypm) = - ez(nxpp,nypm)
        !    by(nxpm,nypp) = by(nxpp,nypp)
        !    by(nxpm,nyp ) = by(nxpp,nyp )
        !    bz(nxpm,nypp) = bz(nxpp,nypp)
        !    bz(nxpm,nyp ) = bz(nxpp,nyp )
        !    bz(nxpm,nypm) = bz(nxpp,nypm)
        !  ENDIF
        !  IF(model_boundary .NE. 0 .AND. nyp .EQ. 0) THEN
        !    ex(nxpp,nypm) = - ex(nxpp,nypp)
        !    ex(nxp ,nypm) = - ex(nxp ,nypp)
        !    ez(nxpp,nypm) = - ez(nxpp,nypp)
        !    ez(nxp ,nypm) = - ez(nxp ,nypp)
        !    ez(nxpm,nypm) = - ez(nxpm,nypp)
        !    bx(nxpp,nypm) = bx(nxpp,nypp)
        !    bx(nxp ,nypm) = bx(nxp ,nypp)
        !    bz(nxpp,nypm) = bz(nxpp,nypp)
        !    bz(nxp ,nypm) = bz(nxp ,nypp)
        !    bz(nxpm,nypm) = bz(nxpm,nypp)
        !  ENDIF
         !IF(model_boundary .EQ. 0 .OR. (nxp .NE. 0 .OR. nyp .NE. 0)) then
         IF( nxp .EQ. 0  ) sx2m =0.d0
         IF( nyp .EQ. 0  ) sy2m =0.d0
          exx = ex(nxpp,nypp)*dx*sy2p + ex(nxp ,nypp)*dx1*sy2p &
               + ex(nxpp,nyp )*dx*sy2  + ex(nxp ,nyp )*dx1*sy2  &
               + ex(nxpp,nypm)*dx*sy2m + ex(nxp ,nypm)*dx1*sy2m

          eyy = ey(nxpp,nypp)*sx2p*dy + ey(nxpp,nyp )*sx2p*dy1 &
               + ey(nxp ,nypp)*sx2 *dy + ey(nxp ,nyp )*sx2 *dy1 &
               + ey(nxpm,nypp)*sx2m*dy + ey(nxpm,nyp )*sx2m*dy1

          ezz = ez(nxpp,nypp)*sx2p*sy2p + ez(nxpp,nyp )*sx2p*sy2 &
               + ez(nxpp,nypm)*sx2p*sy2m + ez(nxp ,nypp)*sx2 *sy2p &
               + ez(nxp ,nyp )*sx2 *sy2  + ez(nxp ,nypm)*sx2 *sy2m &
               + ez(nxpm,nypp)*sx2m*sy2p + ez(nxpm,nyp )*sx2m*sy2 &
               + ez(nxpm,nypm)*sx2m*sy2m

          bxx = bx(nxpp,nypp)*dx*sy2p + bx(nxp ,nypp)*dx1*sy2p &
               + bx(nxpp,nyp )*dx*sy2  + bx(nxp ,nyp )*dx1*sy2 &
               + bx(nxpp,nypm)*dx*sy2m + bx(nxp ,nypm)*dx1*sy2m

          byy = by(nxpp,nypp)*sx2p*dy + by(nxpp,nyp )*sx2p*dy1 &
               + by(nxp ,nypp)*sx2 *dy + by(nxp ,nyp )*sx2 *dy1 &
               + by(nxpm,nypp)*sx2m*dy + by(nxpm,nyp )*sx2m*dy1

          bzz = bz(nxpp,nypp)*sx2p*sy2p + bz(nxp ,nypp)*sx2 *sy2p &
               + bz(nxpm,nypp)*sx2m*sy2p &
               + bz(nxpp,nyp )*sx2p*sy2  + bz(nxp ,nyp )*sx2 *sy2 &
               + bz(nxpm,nyp )*sx2m*sy2  &
               + bz(nxpp,nypm)*sx2p*sy2m + bz(nxp ,nypm)*sx2 *sy2m &
               + bz(nxpm,nypm)*sx2m*sy2m

        !  ELSEIF(model_boundary .NE. 0 .AND. (nxp .EQ. 0 .AND. nyp .NE. 0)) then
        !    exx = ex(nxpp,nypp)*dx*sy2p + ex(nxp ,nypp)*dx1*sy2p &
        !         + ex(nxpp,nyp )*dx*sy2  + ex(nxp ,nyp )*dx1*sy2  &
        !         + ex(nxpp,nypm)*dx*sy2m + ex(nxp ,nypm)*dx1*sy2m
         !
        !    eyy = ey(nxpp,nypp)*sx2p*dy + ey(nxpp,nyp )*sx2p*dy1 &
        !         + ey(nxp ,nypp)*sx2 *dy + ey(nxp ,nyp )*sx2 *dy1 &
        !         + ey(nxpm,nypp)*sx2m*dy + ey(nxpm,nyp )*sx2m*dy1
         !
        !    ezz = ez(nxpp,nypp)*sx2p*sy2p + ez(nxpp,nyp )*sx2p*sy2 &
        !         + ez(nxpp,nypm)*sx2p*sy2m + ez(nxp ,nypp)*sx2 *sy2p &
        !         + ez(nxp ,nyp )*sx2 *sy2  + ez(nxp ,nypm)*sx2 *sy2m &
        !         + ez(nxpm,nypp)*sx2m*sy2p + ez(nxpm,nyp )*sx2m*sy2 &
        !         + ez(nxpm,nypm)*sx2m*sy2m
         !
        !    bxx = bx(nxpp,nypp)*dx*sy2p + bx(nxp ,nypp)*dx1*sy2p &
        !         + bx(nxpp,nyp )*dx*sy2  + bx(nxp ,nyp )*dx1*sy2 &
        !         + bx(nxpp,nypm)*dx*sy2m + bx(nxp ,nypm)*dx1*sy2m
         !
        !    byy = by(nxpp,nypp)*sx2p*dy + by(nxpp,nyp )*sx2p*dy1 &
        !         + by(nxp ,nypp)*sx2 *dy + by(nxp ,nyp )*sx2 *dy1 &
        !         + by(nxpm,nypp)*sx2m*dy + by(nxpm,nyp )*sx2m*dy1
         !
        !    bzz = bz(nxpp,nypp)*sx2p*sy2p + bz(nxp ,nypp)*sx2 *sy2p &
        !         + bz(nxpm,nypp)*sx2m*sy2p &
        !         + bz(nxpp,nyp )*sx2p*sy2  + bz(nxp ,nyp )*sx2 *sy2 &
        !         + bz(nxpm,nyp )*sx2m*sy2  &
        !         + bz(nxpp,nypm)*sx2p*sy2m + bz(nxp ,nypm)*sx2 *sy2m &
        !         + bz(nxpm,nypm)*sx2m*sy2m
        !      ENDIF

       ELSE IF(dx .LE. 0.5d0 .AND. dy .GE. 0.5d0) THEN
        !  IF(model_boundary .NE. 0 .AND. nxp .EQ. 0)  THEN
        !    ey(nxpm,nypp ) = - ey(nxpp,nypp )
        !    ey(nxpm,nyp  ) = - ey(nxpp,nyp  )
        !    ez(nxpm,nyppp) = - ez(nxpp,nyppp)
        !    ez(nxpm,nypp ) = - ez(nxpp,nypp )
        !    ez(nxpm,nyp  ) = - ez(nxpp,nyp  )
        !    by(nxpm,nypp ) = by(nxpp,nypp )
        !    by(nxpm,nyp  ) = by(nxpp,nyp  )
        !    bz(nxpm,nyppp) = bz(nxpp,nyppp)
        !    bz(nxpm,nypp ) = bz(nxpp,nypp )
        !    bz(nxpm,nyp  ) = bz(nxpp,nyp  )
        !  ENDIF
        !  IF(model_boundary .NE. 0 .AND. nyp .EQ. nymax-1) THEN
        !    ex(nxpp,nyppp) = - ex(nxpp,nyp)
        !    ex(nxp ,nyppp) = - ex(nxp ,nyp)
        !    ez(nxpp,nyppp) = - ez(nxpp,nyp)
        !    ez(nxp ,nyppp) = - ez(nxp ,nyp)
        !    ez(nxpm,nyppp) = - ez(nxpm,nyp)
        !    bx(nxpp,nyppp) = bx(nxpp,nyp)
        !    bx(nxp ,nyppp) = bx(nxp ,nyp)
        !    bz(nxpp,nyppp) = bz(nxpp,nyp)
        !    bz(nxp ,nyppp) = bz(nxp ,nyp)
        !    bz(nxpm,nyppp) = bz(nxpm,nyp)
        !  ENDIF
        IF( nxp .EQ. 0  ) sx2m =0.d0
        IF( nyp .EQ. nymax-1  ) sy2p =0.d0
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

       ELSE IF(dx .GE. 0.5d0 .AND. dy .LE. 0.5d0) THEN
        !  IF(model_boundary .NE. 0 .AND. nxp .EQ. nxmax-1)  THEN
        !    ey(nxppp,nypp) = - ey(nxp ,nypp)
        !    ey(nxppp,nyp ) = - ey(nxp ,nyp )
        !    ez(nxppp,nypp) = - ez(nxp ,nypp)
        !    ez(nxppp,nyp ) = - ez(nxp ,nyp )
        !    ez(nxppp,nypm) = - ez(nxp ,nypm)
        !    by(nxppp,nypp) = by(nxp ,nypp)
        !    by(nxppp,nyp ) = by(nxp ,nyp )
        !    bz(nxppp,nypp) = bz(nxp ,nypp)
        !    bz(nxppp,nyp ) = bz(nxp ,nyp )
        !    bz(nxppp,nypm) = bz(nxp ,nypm)
        !  ENDIF
        !  IF(model_boundary .NE. 0 .AND. nyp .EQ. 0) THEN
        !    ex(nxpp ,nypm) = - ex(nxpp ,nypp)
        !    ex(nxp  ,nypm) = - ex(nxp  ,nypp)
        !    ez(nxppp,nypm) = - ez(nxppp,nypp)
        !    ez(nxpp ,nypm) = - ez(nxpp ,nypp)
        !    ez(nxp  ,nypm) = - ez(nxp  ,nypp)
        !    bx(nxpp ,nypm) = bx(nxpp ,nypp)
        !    bx(nxp  ,nypm) = bx(nxp  ,nypp)
        !    bz(nxppp,nypm) = bz(nxppp,nypp)
        !    bz(nxpp ,nypm) = bz(nxpp ,nypp)
        !    bz(nxp  ,nypm) = bz(nxp  ,nypp)
        !  ENDIF
        IF( nxp .EQ. nxmax-1  ) sx2p =0.d0
        IF( nyp .EQ. 0  ) sy2m =0.d0
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
       ELSE
          !  IF(model_boundary .NE. 0 .AND. nxp .EQ. nxmax-1)  THEN
          !    ey(nxppp,nypp ) = - ey(nxp ,nypp )
          !    ey(nxppp,nyp  ) = - ey(nxp ,nyp  )
          !    ez(nxppp,nyppp) = - ez(nxp ,nyppp)
          !    ez(nxppp,nypp ) = - ez(nxp ,nypp )
          !    ez(nxppp,nyp  ) = - ez(nxp ,nyp  )
          !    by(nxppp,nypp ) = by(nxp ,nypp )
          !    by(nxppp,nyp  ) = by(nxp ,nyp  )
          !    bz(nxppp,nyppp) = bz(nxp ,nyppp)
          !    bz(nxppp,nypp ) = bz(nxp ,nypp )
          !    bz(nxppp,nyp  ) = bz(nxp ,nyp  )
          !  ENDIF
          !  IF(model_boundary .NE. 0 .AND. nyp .EQ. nymax-1)  THEN
          !    ex(nxpp ,nyppp) = - ex(nxpp ,nyp)
          !    ex(nxp  ,nyppp) = - ex(nxp  ,nyp)
          !    ez(nxppp,nyppp) = - ez(nxppp,nyp)
          !    ez(nxpp ,nyppp) = - ez(nxpp ,nyp)
          !    ez(nxp  ,nyppp) = - ez(nxp  ,nyp)
          !    bx(nxpp ,nyppp) = bx(nxpp ,nyp)
          !    bx(nxp  ,nyppp) = bx(nxp  ,nyp)
          !    bz(nxppp,nyppp) = bz(nxppp,nyp)
          !    bz(nxpp ,nyppp) = bz(nxpp ,nyp)
          !    bz(nxp  ,nyppp) = bz(nxp  ,nyp)
          !  ENDIF
          IF( nxp .EQ. nxmax-1  ) sx2p =0.d0
          IF( nyp .EQ. nymax-1  ) sy2p =0.d0
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
       ENDIF

       ! push particles by using Buneman-Boris method
       ! gamma is lorentz factor
       bb2 = bxx ** 2 + byy ** 2 + bzz ** 2
       gamma = 1.d0 / sqrt(1.d0 - (vx(np)**2+vy(np)**2+vz(np)**2)/vcfact**2)

       vxm = vx(np) * gamma + 0.5D0 * ctom * exx * dt
       vym = vy(np) * gamma + 0.5D0 * ctom * eyy * dt
       vzm = vz(np) * gamma + 0.5D0 * ctom * ezz * dt

       gamma = 1.d0 / sqrt(1.d0 + (vxm**2+vym**2+vzm**2) / vcfact**2)

       vxzero = vxm + 0.5D0 * ctom * (vym * bzz - vzm * byy) * dt * gamma
       vyzero = vym + 0.5D0 * ctom * (vzm * bxx - vxm * bzz) * dt * gamma
       vzzero = vzm + 0.5D0 * ctom * (vxm * byy - vym * bxx) * dt * gamma

       vxp = vxm + 1.0d0/(1.0d0 + 0.25d0 * (ctom * dt * gamma) ** 2 * bb2) &
            * ctom * (vyzero * bzz - vzzero * byy) * dt * gamma
       vyp = vym + 1.0d0/(1.0d0 + 0.25d0 * (ctom * dt * gamma) ** 2 * bb2) &
            * ctom * (vzzero * bxx - vxzero * bzz) * dt * gamma
       vzp = vzm + 1.0d0/(1.0d0 + 0.25d0 * (ctom * dt * gamma) ** 2 * bb2) &
            * ctom * (vxzero * byy - vyzero * bxx) * dt * gamma

       gamma = 1.d0 / sqrt(1.d0 + (vxp**2+vyp**2+vzp**2)/vcfact**2)
       vx(np) = vxp * gamma + 0.5D0 * ctom * exx * dt
       vy(np) = vyp * gamma + 0.5D0 * ctom * eyy * dt
       vz(np) = vzp * gamma + 0.5D0 * ctom * ezz * dt

       !push particles
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

    END DO

  END SUBROUTINE push

  !***********************************************************************
  SUBROUTINE bound_periodic(npmax,x,y,z,x1,x2,y1,y2,z1,z2,alx,aly,alz)
    !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(npmax) :: x, y, z
    REAL(8) :: alx, aly, alz, x1, x2, y1, y2, z1, z2
    INTEGER :: npmax, np

    DO np = 1, npmax
       IF( x(np) .LT. x1 ) THEN
          DO WHILE(x(np) .LT. x1)
             x(np) = x(np) + alx
          END DO
       ELSEIF( x(np) .GT. x2 ) THEN
          DO WHILE(x(np) .GT. x2)
             x(np) = x(np) - alx
          END DO
       ENDIF

       IF( y(np) .LT. y1 ) THEN
          DO WHILE(y(np) .LT. y1)
             y(np) = y(np) + aly
          END DO
       ELSEIF( y(np) .GT. y2 ) THEN
          DO WHILE(y(np) .GT. y2)
             y(np) = y(np) - aly
          END DO
       ENDIF

       IF( z(np) .LT. z1 ) THEN
          DO WHILE(z(np) .LT. z1)
            z(np) = z(np) + alz
          END DO
       ELSEIF( z(np) .GT. z2 ) THEN
          DO WHILE(z(np) .GT. z2)
             z(np) = z(np) - alz
          END DO
       ENDIF
    END DO

  END SUBROUTINE bound_periodic

  !***********************************************************************
  SUBROUTINE bound_reflective(npmax,x,y,z,vx,vy, &
       x1,x2,y1,y2,z1,z2,alx,aly,alz)
    !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(npmax) :: x, y, z, vx, vy
    REAL(8) :: x1, x2, y1, y2, z1, z2, alx, aly, alz
    INTEGER :: npmax, np

    DO np = 1, npmax
       IF( x(np) .LT. x1 ) THEN
          x(np) = -x(np)
          vx(np) = -vx(np)
       ELSEIF( x(np) .GT. x2 ) THEN
          x(np) = alx - (x(np) - alx)
          vx(np) = -vx(np)
       ENDIF
       IF( y(np) .LT. y1 ) THEN
          y(np) = -y(np)
          vy(np) = -vy(np)
       ELSEIF( y(np) .GT. y2 ) THEN
          y(np) = aly - (y(np) - aly)
          vy(np) = -vy(np)
       ENDIF

         IF( z(np) .LT. z1 ) THEN
            DO WHILE(z(np) .LT. z1)
               z(np) = z(np) + alz
          END DO
        ELSEIF( z(np) .GT. z2 ) THEN
            DO WHILE(z(np) .GT. z2)
               z(np) = z(np) - alz
            END DO
         ENDIF
    END DO

  END SUBROUTINE bound_reflective

  !***********************************************************************
  SUBROUTINE source(npmax,nxmax,nymax,x,y,rho,chrg,model_boundary)
    !***********************************************************************

    IMPLICIT NONE
    REAL(8), DIMENSION(npmax)        :: x, y
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: rho
    REAL(8) :: chrg, factor, dx, dy, dz, sx2, sy2, sx2p, sy2p, sx2m, sy2m,&
         dx1,dy1
    INTEGER :: npmax, nxmax, nymax, model_boundary
    INTEGER:: np, nxp, nyp, nx, ny, nxpp, nxpm, nypp, nypm, nxppp, nyppp

    IF(npmax.EQ.0) THEN
       factor=chrg*DBLE(nxmax)*DBLE(nymax)
    ELSE
       factor=chrg*DBLE(nxmax)*DBLE(nymax)/DBLE(npmax)
    END IF

    !*poption parallel, psum(rho)

    DO np = 1, npmax

       nxp = x(np)
       nyp = y(np)

       dx  = x(np) - DBLE(nxp)
       dy  = y(np) - DBLE(nyp)
       dx1 = 1.0d0 - dx
       dy1 = 1.0d0 - dy
       IF(dx .LE. 0.5d0) THEN
          sx2  = 3.0d0/4 - dx ** 2
          sx2p = 1.0d0/2 * (1.0d0/2 + dx) ** 2
          sx2m = 1.0d0/2 * (1.0d0/2 - dx) ** 2
       ELSE
          sx2  = 3.0d0/4 - (1.0d0 - dx) ** 2
          sx2p = 1.0d0/2 * (-1.0d0/2 + dx) ** 2
          sx2m = 1.0d0/2 * (3.0d0/2 - dx) ** 2
       ENDIF
       IF(dy .LE. 0.5d0) THEN
          sy2  = 3.0d0/4 - dy ** 2
          sy2p = 1.0d0/2 * (1.0d0/2 + dy) ** 2
          sy2m = 1.0d0/2 * (1.0d0/2 - dy) ** 2
       ELSE
          sy2  = 3.0d0/4 - (1.0d0 - dy) ** 2
          sy2p = 1.0d0/2 * (-1.0d0/2 + dy) ** 2
          sy2m = 1.0d0/2 * (3.0d0/2 - dy) ** 2
       ENDIF

       nxpm = nxp - 1
       nxpp = nxp + 1
       nxppp= nxp + 2
       nypm = nyp - 1
       nypp = nyp + 1
       nyppp= nyp + 2

       IF(model_boundary.EQ.0) THEN ! periodic
          IF( nxp .EQ. 0  ) nxpm = nxmax - 1
          IF( nyp .EQ. 0  ) nypm = nymax - 1
          IF( nxp .EQ. nxmax - 1) nxppp = 1
          IF( nyp .EQ. nymax - 1) nyppp = 1
       ELSE   ! reflective:
          IF( nxp .EQ. 0  ) nxpm = 0
          IF( nyp .EQ. 0  ) nypm = 0
          IF( nxp .EQ. nxmax-1) nxppp = nxmax
          IF( nyp .EQ. nymax-1) nyppp = nymax
       END IF

       IF(dx .LT. 0.5d0 .AND. dy .LT. 0.5d0) THEN
         ! use Mirror image method
         IF(model_boundary .NE. 0 .AND. nxp .EQ. 0) THEN
           sx2p = sx2p - sx2m
           sx2 =  0.d0
           sx2m = 0.d0!sx2m - sx2p
         ELSEIF(model_boundary .NE. 0 .AND. nxp .EQ. nxmax-1) THEN
           sx2p = 0.d0
           sx2 =  sx2
           sx2m = sx2m!sx2m - sx2p
         ENDIF
         IF(model_boundary .NE. 0 .AND. nyp .EQ. 0) THEN
           sy2p = sy2p - sy2m
           sy2 =  0.d0
           sy2m = 0.d0!sy2m - sy2p
         ELSEIF(model_boundary .NE. 0 .AND. nyp .EQ. nymax-1) THEN
           sy2p = 0.d0
           sy2 =  sy2
           sy2m = sy2m!sx2m - sx2p
         ENDIF

          rho(nxpm,nypm) = rho(nxpm,nypm) + sx2m * sy2m * factor
          rho(nxpm,nyp ) = rho(nxpm,nyp ) + sx2m * sy2  * factor
          rho(nxpm,nypp) = rho(nxpm,nypp) + sx2m * sy2p * factor
          rho(nxp ,nypm) = rho(nxp ,nypm) + sx2  * sy2m * factor
          rho(nxp ,nyp ) = rho(nxp ,nyp ) + sx2  * sy2  * factor
          rho(nxp ,nypp) = rho(nxp ,nypp) + sx2  * sy2p * factor
          rho(nxpp,nypm) = rho(nxpp,nypm) + sx2p * sy2m * factor
          rho(nxpp,nyp ) = rho(nxpp,nyp ) + sx2p * sy2  * factor
          rho(nxpp,nypp) = rho(nxpp,nypp) + sx2p * sy2p * factor
       ELSE IF(dx .LT. 0.5d0 .AND. dy .GT. 0.5d0) THEN
         IF(model_boundary .NE. 0 .AND. nxp .EQ. 0)  THEN
           sx2p = sx2p - sx2m
           sx2  = 0.d0
           sx2m = 0.d0!sx2m - sx2p
       ELSEIF(model_boundary .NE. 0 .AND. nxp .EQ. nxmax-1) THEN
         sx2p = 0.d0
         sx2  = sx2
         sx2m = sx2m!sx2m - sx2p
       ENDIF
         IF(model_boundary .NE. 0 .AND. nyp .EQ. nymax-1) THEN
           sy2m = sy2m - sy2p
           sy2  = 0.d0
           sy2p = 0.d0!sy2p - sy2m
       ELSEIF(model_boundary .NE. 0 .AND. nyp .EQ. 0) THEN
         sy2p = sy2p
         sy2  = sy2
         sy2m = 0.d0!sx2m - sx2p
       ENDIF
          rho(nxpm ,nyp  ) = rho(nxpm ,nyp  ) + sx2m * sy2m * factor
          rho(nxpm ,nypp ) = rho(nxpm ,nypp ) + sx2m * sy2  * factor
          rho(nxpm ,nyppp) = rho(nxpm ,nyppp) + sx2m * sy2p * factor
          rho(nxp  ,nyp  ) = rho(nxp  ,nyp  ) + sx2  * sy2m * factor
          rho(nxp  ,nypp ) = rho(nxp  ,nypp ) + sx2  * sy2  * factor
          rho(nxp  ,nyppp) = rho(nxp  ,nyppp) + sx2  * sy2p * factor
          rho(nxpp ,nyp  ) = rho(nxpp ,nyp  ) + sx2p * sy2m * factor
          rho(nxpp ,nypp ) = rho(nxpp ,nypp ) + sx2p * sy2  * factor
          rho(nxpp ,nyppp) = rho(nxpp ,nyppp) + sx2p * sy2p * factor
       ELSE IF(dx .GT. 0.5d0 .AND. dy .LT. 0.5d0) THEN
         IF(model_boundary .NE. 0 .AND. nxp .EQ. nxmax-1)  THEN
           sx2m = sx2m - sx2p
           sx2  = 0.d0
           sx2p = 0.d0!sx2p - sx2m
         ELSEIF(model_boundary .NE. 0 .AND. nyp .EQ. 0) THEN
           sx2p = sx2p
           sx2  = sx2
           sx2m = 0.d0!sx2m - sx2p
         ENDIF
         IF(model_boundary .NE. 0 .AND. nyp .EQ. 0) THEN
           sy2p = sy2p - sy2m
           sy2 =  0.d0
           sy2m = 0.d0!sy2m - sy2p
       ELSEIF(model_boundary .NE. 0 .AND. nyp .EQ. nymax-1) THEN
         sy2p = 0.d0
         sy2 =  sy2
         sy2m = sy2m!sx2m - sx2p
       ENDIF

          rho(nxp  ,nypm) = rho(nxp  ,nypm ) + sx2m * sy2m * factor
          rho(nxp  ,nyp ) = rho(nxp  ,nyp  ) + sx2m * sy2  * factor
          rho(nxp  ,nypp) = rho(nxp  ,nypp ) + sx2m * sy2p * factor
          rho(nxpp ,nypm) = rho(nxpp ,nypm ) + sx2  * sy2m * factor
          rho(nxpp ,nyp ) = rho(nxpp ,nyp  ) + sx2  * sy2  * factor
          rho(nxpp ,nypp) = rho(nxpp ,nypp ) + sx2  * sy2p * factor
          rho(nxppp,nypm) = rho(nxppp,nypm ) + sx2p * sy2m * factor
          rho(nxppp,nyp ) = rho(nxppp,nyp  ) + sx2p * sy2  * factor
          rho(nxppp,nypp) = rho(nxppp,nypp ) + sx2p * sy2p * factor
       ELSE
         IF(model_boundary .NE. 0 .AND. nxp .EQ. nxmax-1)  THEN
           sx2m = sx2m - sx2p
           sx2  = 0.d0
           sx2p = 0.d0!sx2p - sx2m
         ELSEIF(model_boundary .NE. 0 .AND. nxp .EQ. 0) THEN
           sx2p = sy2p
           sx2  = sx2
           sx2m = 0.d0!sx2m - sx2p
         ENDIF
         IF(model_boundary .NE. 0 .AND. nyp .EQ. nymax-1) THEN
           sy2m = sy2m - sy2p
           sy2  = 0.d0
           sy2p = 0.d0!sy2p - sy2m
         ELSEIF(model_boundary .NE. 0 .AND. nyp .EQ. 0) THEN
           sy2p = sy2p
           sy2 =  sy2
           sy2m = 0.d0!sy2m - sy2p
         ENDIF
          rho(nxp  ,nyp  ) = rho(nxp  ,nyp  ) + sx2m * sy2m * factor
          rho(nxp  ,nypp ) = rho(nxp  ,nypp ) + sx2m * sy2  * factor
          rho(nxp  ,nyppp) = rho(nxp  ,nyppp) + sx2m * sy2p * factor
          rho(nxpp ,nyp  ) = rho(nxpp ,nyp  ) + sx2  * sy2m * factor
          rho(nxpp ,nypp ) = rho(nxpp ,nypp ) + sx2  * sy2  * factor
          rho(nxpp ,nyppp) = rho(nxpp ,nyppp) + sx2  * sy2p * factor
          rho(nxppp,nyp  ) = rho(nxppp,nyp  ) + sx2p * sy2m * factor
          rho(nxppp,nypp ) = rho(nxppp,nypp ) + sx2p * sy2  * factor
          rho(nxppp,nyppp) = rho(nxppp,nyppp) + sx2p * sy2p * factor
       ENDIF
    END DO
    END SUBROUTINE source

  !***********************************************************************
  SUBROUTINE boundary_rho(nxmax,nymax,rho,model_boundary)
    !***********************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: rho
    INTEGER:: nxmax,nymax,model_boundary
    INTEGER:: nx,ny

    !!..... set charge densities at the boundary

    IF(model_boundary.EQ.0) THEN  ! periodic
       DO ny = 0, nymax
          rho(0,ny) = rho(0,ny) + rho(nxmax,ny)
          rho(nxmax,ny) = rho(0,ny)
       END DO
       DO nx = 0, nxmax
          rho(nx,0) = rho(nx,0) + rho(nx,nymax)
          rho(nx,nymax) = rho(nx,0)
       END DO
    !ELSE                         ! reflecting
      !rho(0,:) = 0.d0
      !rho(nxmax,:) = 0.d0
      !rho(:,0) = 0.d0
      !rho(:,nymax) = 0.d0
       !DO ny = 1, nymax-1
        !  rho(0,ny)     = 2.D0 * rho(0,ny)
        !  rho(nxmax,ny) = 2.D0 * rho(nxmax,ny)
       !END DO
       !DO nx = 1, nxmax-1
        !  rho(nx,0)     = 2.D0 * rho(nx,0)
        !  rho(nx,nymax) = 2.D0 * rho(nx,nymax)
       !END DO
       !rho(0,0)         = 4.D0 * rho(0,0)
       !rho(0,nymax)     = 4.D0 * rho(0,nymax)
       !rho(nxmax,0)     = 4.D0 * rho(nxmax,0)
       !rho(nxmax,nymax) = 4.D0 * rho(nxmax,nymax)
    END IF

  END SUBROUTINE boundary_rho

  !***********************************************************************
  SUBROUTINE current(npmax,nxmax,nymax,x,y,vx,vy,vz,chrg,jx,jy,jz, &
       model_boundary)
    !***********************************************************************
    IMPLICIT NONE

    REAL(8), DIMENSION(npmax) :: x, y, vx, vy, vz
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: jx, jy, jz
    REAL(8) :: chrg, dt, dx, dy, dx1, dy1, &
         sx1p, sy1p, sx1m, sy1m, sx2, sy2, sx2p, sy2p, sx2m, sy2m, factor
    INTEGER :: npmax, nxmax, nymax, model_boundary
    INTEGER :: np, nxp, nyp, nxpm, nypm, nxpp, nypp, nxppp, nyppp, nx, ny

    IF(npmax.EQ.0) THEN
       factor=chrg*DBLE(nxmax)*DBLE(nymax)
    ELSE
       factor=chrg*DBLE(nxmax)*DBLE(nymax)/DBLE(npmax)
    END IF

    DO np = 1, npmax

       nxp = x(np)
       nyp = y(np)
       dx = x(np) - DBLE(nxp)
       dy = y(np) - DBLE(nyp)
       dx1 = 1.0d0 - dx
       dy1 = 1.0d0 - dy
       IF(dx .LE. 0.5d0) THEN
          sx2  = 3.0d0/4 - dx ** 2
          sx2p = 1.0d0/2 * (1.0d0/2 + dx) ** 2
          sx2m = 1.0d0/2 * (1.0d0/2 - dx) ** 2
       ELSE
          sx2  = 3.0d0/4 - (dx - 1.0d0) ** 2
          sx2p = 1.0d0/2 * (-1.0d0/2 + dx) ** 2
          sx2m = 1.0d0/2 * (3.0d0/2 - dx) ** 2
       ENDIF
       IF(dy .LE. 0.5d0) THEN
          sy2  = 3.0d0/4 - dy ** 2
          sy2p = 1.0d0/2 * (1.0d0/2 + dy) ** 2
          sy2m = 1.0d0/2 * (1.0d0/2 - dy) ** 2
       ELSE
          sy2  = 3.0d0/4 - (dy - 1.0d0) ** 2
          sy2p = 1.0d0/2 * (-1.0d0/2 + dy) ** 2
          sy2m = 1.0d0/2 * (3.0d0/2 - dy) ** 2
       ENDIF
       nxpm = nxp - 1
       nxpp = nxp + 1
       nypm = nyp - 1
       nypp = nyp + 1
       nxppp = nxp + 2
       nyppp = nyp + 2

       IF(model_boundary.EQ.0) THEN ! periodic
          IF( nxp .EQ. 0  ) nxpm = nxmax - 1
          IF( nyp .EQ. 0  ) nypm = nymax - 1
          IF( nxp .EQ. nxmax-1) nxppp=1
          IF( nyp .EQ. nymax-1) nyppp=1
       ELSE   ! reflective:
          IF( nxp .EQ. 0  ) nxpm=0
          IF( nyp .EQ. 0  ) nypm=0
          IF( nxp .EQ. nxmax-1) nxppp=nxmax
          IF( nyp .EQ. nymax-1) nyppp=nymax
       END IF

       IF (dx .LE. 0.5d0 .AND. dy .LE. 0.5d0) THEN
         ! use Mirror image method
         IF(model_boundary .NE. 0 .AND. nxp .EQ. 0)  THEN
           sx2p = sx2p - sx2m
           sx2m = 0.d0
           sx2 = 0.d0
           dx1 = 2.0d0 * dx1
        ELSE IF(model_boundary .NE. 0 .AND. nxp .EQ. nxmax-1)  THEN
             sx2p = sx2p
             sx2m = 0.d0
             sx2 = sx2
             dx = 2.0d0 * dx
           ENDIF
         IF(model_boundary .NE. 0 .AND. nyp .EQ. 0) THEN
           sy2p = sy2p - sy2m
           sy2m = 0.d0
           sy2 = 0.d0
           dy1 = 2.0d0 * dy1
       ELSE IF(model_boundary .NE. 0 .AND. nyp .EQ. nymax-1) THEN
         sy2p = sy2p
         sy2m = 0.d0
         sy2 = sx2
         dy = 2.0d0 * dy
       ENDIF
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

       ELSE IF(dx .LE. 0.5d0 .AND. dy .GE. 0.5d0) THEN
         IF(model_boundary .NE. 0 .AND. nxp .EQ. 0)  THEN
           sx2p = sx2p - sx2m
           sx2 = 0.d0
           sx2m = 0.d0
           dx1 = 2.0d0 * dx1
         ELSE IF(model_boundary .NE. 0 .AND. nxp .EQ. nxmax-1)  THEN
              sx2p = sx2p
              sx2m = 0.d0
              sx2 = sx2
              dx = 2.0d0 * dx
         ENDIF
         IF(model_boundary .NE. 0 .AND. nyp .EQ. nymax-1) THEN
           sy2m = sy2m - sy2p
           sy2 = 0.d0
           sy2p = 0.d0
           dy = 2.0d0 * dy
         ELSE IF(model_boundary .NE. 0 .AND. nyp .EQ. 0) THEN
           sy2p = sy2p
           sy2m = 0.d0
           sy2 = sx2
           dy1 = 2.0d0 * dy1
         ENDIF
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

       ELSE IF(dx .GE. 0.5d0 .AND. dy .LE. 0.5d0) THEN
         IF(model_boundary .NE. 0 .AND. nxp .EQ. nxmax-1)  THEN
           sx2m = sx2m - sx2p
           sx2 = 0.d0
           sx2p = 0.d0
           dx = 2.0d0 * dx
         ELSE IF(model_boundary .NE. 0 .AND. nxp .EQ. 0) THEN
           sx2p = sx2p
           sx2m = 0.d0
           sx2 = sx2
           dx1 = 2.0d0 * dx1
         ENDIF
         IF(model_boundary .NE. 0 .AND. nyp .EQ. 0) THEN
           sy2p = sy2p - sy2m
           sy2 = 0.d0
           sy2m = 0.d0
           dy1 = 2.0d0 * dy1
         ELSE IF(model_boundary .NE. 0 .AND. nyp .EQ. nymax-1) THEN
           sy2p = sy2p
           sy2m = 0.d0
           sy2 = sx2
           dy = 2.0d0 * dy
         ENDIF
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

       ELSE
         IF(model_boundary .NE. 0 .AND. nxp .EQ. nxmax-1) THEN
           sx2m = sx2m - sx2p
           sx2 = 0.d0
           sx2p = 0.d0
           dx = 2.0d0 * dx
         ELSE IF(model_boundary .NE. 0 .AND. nxp .EQ. 0) THEN
           sx2p = sx2p
           sx2m = 0.d0
           sx2 = sx2
           dx1 = 2.0d0 * dx1
         ENDIF
         IF(model_boundary .NE. 0 .AND. nyp .EQ. nymax-1) THEN
           sy2m = sy2m - sy2p
           sy2 = 0.d0
           sy2p = 0.d0
           dy = 2.0d0 * dy
         ELSE IF(model_boundary .NE. 0 .AND. nyp .EQ. 0) THEN
           sy2p = sy2p
           sy2m = 0.d0
           sy2 = sx2
           dy1 = 2.0d0 * dy1
         ENDIF
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

       ENDIF

    END DO

  END SUBROUTINE current

  !***********************************************************************
  SUBROUTINE antenna(nxmax,nymax,jxant,jyant,jzant,phxant,phyant,phzant, &
       omega,time,jx,jy,jz)
    !***********************************************************************
    IMPLICIT NONE

    INTEGER :: nxmax,nymax
    REAL(8) :: jxant,jyant,jzant,phxant,phyant,phzant,omega,time
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: jx, jy, jz
    REAL(8) :: jxt,jyt,jzt
    INTEGER :: nx,ny

    jxt = jxant * COS (omega * time + phxant)
    jyt = jyant * COS (omega * time + phyant)
    jzt = jzant * COS (omega * time + phzant)

    ! add antenna current density

    DO ny=5,10
       jy(3,ny) = jy(3,ny) + 0.5d0 * jyt
       jy(4,ny) = jy(4,ny) + 0.5d0 * jyt
       jz(3,ny) = jz(3,ny) + 0.5d0 * jzt
       jz(4,ny) = jz(4,ny) + 0.5d0 * jzt
    END DO

  END SUBROUTINE antenna

  !***********************************************************************
  SUBROUTINE boundary_j(nxmax,nymax,jx,jy,jz,model_boundary)
    !***********************************************************************
    IMPLICIT NONE

    INTEGER :: nxmax,nymax,model_boundary
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: jx, jy, jz
    INTEGER :: nx,ny

    IF(model_boundary.EQ.0) THEN     ! periodic
        DO ny = 0, nymax
          jx(0,ny) = jx(0,ny) + jx(nxmax,ny)
          jy(0,ny) = jy(0,ny) + jy(nxmax,ny)
          jz(0,ny) = jz(0,ny) + jz(nxmax,ny)
          jx(nxmax,ny) = jx(0,ny)
          jy(nxmax,ny) = jy(0,ny)
          jz(nxmax,ny) = jz(0,ny)
       END DO
       DO nx = 0, nxmax
          jx(nx,0) = jx(nx,0) + jx(nx,nymax)
          jy(nx,0) = jy(nx,0) + jy(nx,nymax)
          jz(nx,0) = jz(nx,0) + jz(nx,nymax)
          jx(nx,nymax) = jx(nx,0)
          jy(nx,nymax) = jy(nx,0)
          jz(nx,nymax) = jz(nx,0)
       ENDDO
    ELSE                           ! reflective
       DO ny=1,nymax-1
          jx(0,ny) = 2.0d0 * jx(0,ny)
          jy(0,ny) = 2.0d0 * jy(0,ny)
          jz(0,ny) = 2.0d0 * jz(0,ny)
          jx(nxmax,ny) = 2.0d0 * jx(nxmax,ny)
          jy(nxmax,ny) = 2.0d0 * jy(nxmax,ny)
          jz(nxmax,ny) = 2.0d0 * jz(nxmax,ny)
       END DO
       DO nx=1,nxmax-1
          jx(nx,0) = 2.0d0 * jx(nx,0)
          jy(nx,0) = 2.0d0 * jy(nx,0)
          jz(nx,0) = 2.0d0 * jz(nx,0)
          jx(nx,nymax) = 2.0d0 * jx(nx,nymax)
          jy(nx,nymax) = 2.0d0 * jy(nx,nymax)
          jz(nx,nymax) = 2.0d0 * jz(nx,nymax)
       END DO
       jx(0,0) = 4.0d0 * jx(0,0)
       jy(0,0) = 4.0d0 * jy(0,0)
       jz(0,0) = 4.0d0 * jz(0,0)
       jx(0,nymax) = 4.0d0 * jx(0,nymax)
       jy(0,nymax) = 4.0d0 * jy(0,nymax)
       jz(0,nymax) = 4.0d0 * jz(0,nymax)
       jx(nxmax,0) = 4.0d0 * jx(nxmax,0)
       jy(nxmax,0) = 4.0d0 * jy(nxmax,0)
       jz(nxmax,0) = 4.0d0 * jz(nxmax,0)
       jx(nxmax,nymax) = 4.0d0 * jx(nxmax,nymax)
       jy(nxmax,nymax) = 4.0d0 * jy(nxmax,nymax)
       jz(nxmax,nymax) = 4.0d0 * jz(nxmax,nymax)
    END IF

  END SUBROUTINE boundary_j

  !***********************************************************************
  SUBROUTINE phia_periodic(nxmax,nymax,vcfact,dt,phi,phib,jx,jy,jz, &
       Ax,Ay,Az,Axb,Ayb,Azb,Axbb,Aybb,Azbb)
    !***********************************************************************
    !original subroutine
    IMPLICIT NONE
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: phi,phib,jx,jy,jz,Ax,Ay,Az, &
         Axb,Ayb,Azb,Axbb,Aybb,Azbb
    INTEGER :: nxmax, nymax, nx, ny, nxm, nxp, nyp, nym, nypm
    REAL(8) :: vcfact, dt, dph, x, y
    ! Solution of maxwell equation in the A-phi formulation by difference method
    ! vcfact is the ratio of the light speed to lattice parameter times
    ! plasma frequency

    DO nx = 0, nxmax
       DO ny = 0, nymax

          nxm = nx - 1
          nxp = nx + 1
          nym = ny - 1
          nyp = ny + 1

          IF( nx .EQ. 0  )    nxm = nxmax - 1
          IF( nx .EQ. nxmax ) nxp = 1
          IF( ny .EQ. 0  )    nym = nymax - 1
          IF( ny .EQ. nymax ) nyp = 1

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

       END DO
    END DO

  END SUBROUTINE phia_periodic

  !***********************************************************************
  SUBROUTINE phia_reflective(nxmax,nymax,vcfact,dt,phi,phib,jx,jy,jz, &
       Ax,Ay,Az,Axb,Ayb,Azb,Axbb,Aybb,Azbb, &
       model_wg,xmin_wg,xmax_wg,ymin_wg,ymax_wg, &
       amp_wg,ph_wg,rot_wg,eli_wg,omega,time,pi,&
       model_boundary,dlen)
    !***********************************************************************
    !original subroutine
    IMPLICIT NONE
    REAL(8), DIMENSION(0:nxmax,0:nymax) :: phi,phib,jx,jy,jz,Ax,Ay,Az, &
         Axb,Ayb,Azb,Axbb,Aybb,Azbb
    INTEGER :: nxmax, nymax, nx, ny, nxm, nxp, nyp, nym,ilen
    INTEGER :: model_wg,model_boundary
    REAL(8) :: xmin_wg,xmax_wg,ymin_wg,ymax_wg,amp_wg,ph_wg,rot_wg,eli_wg
    REAL(8) :: omega,time,pi,dlen,inv
    REAL(8) :: vcfact,dt,dph,x,y,xmax,ymax,yc,ylen,factor,amp_start

    ! Solution of maxwell equation in the A-phi formulation by difference method
    ! vcfact is the ratio of the light speed to lattice parameter times plasma
    ! frequency

    DO nx = 0, nxmax
       DO ny = 0, nymax

          nxm = nx - 1
          nxp = nx + 1
          nym = ny - 1
          nyp = ny + 1
          IF( nx .EQ. 0  )    nxm = nxmax
          IF( nx .EQ. nxmax ) nxp = 0
          IF( ny .EQ. 0  )    nym = nymax
          IF( ny .EQ. nymax ) nyp = 0

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
                    + dt ** 2 * jz(nx,ny) + 2.0d0 * Azb(nx,ny) - Azbb(nx,ny)

       END DO
    END DO
    IF(model_boundary .EQ. 2) THEN !damping A in absorbing boundary
       ilen = int(dlen)
       inv = 1.0d0 / dlen
       DO ny = 1,nymax
          DO nx = nxmax-ilen,nxmax
             x=DBLE(nx)
             xmax=DBLE(nxmax)
             Ax(nx,ny) = Ax(nx,ny)*(-1.0d0*inv**2*x**2 &
                                    +2.0d0*inv**2*(xmax-dlen)*x&
                                    +1.0d0-1.0d0*inv**2*(xmax-dlen)**2)
             Ay(nx,ny) = Ay(nx,ny)*(-1.0d0*inv**2*x**2 &
                                    +2.0d0*inv**2*(xmax-dlen)*x&
                                    +1.0d0-1.0d0*inv**2*(xmax-dlen)**2)
             Az(nx,ny) = Az(nx,ny)*(-1.0d0*inv**2*x**2 &
                                    +2.0d0*inv**2*(xmax-dlen)*x&
                                    +1.0d0-1.0d0*inv**2*(xmax-dlen)**2)
          ENDDO
       ENDDO
       DO nx = 1,nxmax-ilen
          DO ny = nymax-ilen,nymax
             y=DBLE(ny)
             ymax=DBLE(nymax)
             Ax(nx,ny) = Ax(nx,ny)*(-1.0d0*inv**2*y**2 &
                                    +2.0d0*inv**2*(ymax-dlen)*y &
                                    +1.0d0-1.0d0*inv**2*(ymax-dlen)**2)
             Ay(nx,ny) = Ay(nx,ny)*(-1.0d0*inv**2*y**2 &
                                    +2.0d0*inv**2*(ymax-dlen)*y &
                                    +1.0d0-1.0d0*inv**2*(ymax-dlen)**2)
             Az(nx,ny) = Az(nx,ny)*(-1.0d0*inv**2*y**2 &
                                    +2.0d0*inv**2*(ymax-dlen)*y &
                                    +1.0d0-1.0d0*inv**2*(ymax-dlen)**2)
          ENDDO
       ENDDO
       DO nx = 1, nxmax-ilen
          DO ny = 1, ilen
             y=DBLE(ny)
             Ax(nx,ny) = Ax(nx,ny)*(-1.0d0*inv**2*y**2+2.0d0*inv*y)
             Ay(nx,ny) = Ay(nx,ny)*(-1.0d0*inv**2*y**2+2.0d0*inv*y)
             Az(nx,ny) = Az(nx,ny)*(-1.0d0*inv**2*y**2+2.0d0*inv*y)
          ENDDO
       ENDDO
    ENDIF

    !boundary condition for reflection
     !Ax(0,:)=0.d0
     !Ay(0,:)=0.d0
     !Az(0,:)=0.d0
     !Ax(nxmax,:)=0.d0
     !Ay(nxmax,:)=0.d0
     !Az(nxmax,:)=0.d0
     !Ax(:,0)=0.d0
     !Ay(:,0)=0.d0
     !Az(:,0)=0.d0
     !Ax(:,nymax)=0.d0
     !Ay(:,nymax)=0.d0
     !Az(:,nymax)=0.d0
    SELECT CASE(model_wg)
    CASE(0)
       yc=0.5d0*(ymin_wg+ymax_wg)
       ylen=(ymax_wg-ymin_wg)
       IF(ylen .NE. 0) dph=ph_wg/ylen
       amp_start = time * vcfact/100.d0
       IF(amp_start .GE. 1.d0) amp_start=1.0d0
       DO ny=1,nymax
          y=DBLE(ny)
          IF(y.GE.ymin_wg.AND.y.LE.ymax_wg) THEN
             factor=EXP(-12.D0*(y-yc)**2/(ylen)**2)
             Ay(0,ny)=amp_wg*amp_start &
                  *factor*COS(rot_wg*pi/180.D0) &
                  *SIN(omega*time-pi*dph*(y-ymin_wg)/180.D0)
             Az(0,ny)=amp_wg*amp_start &
                  *factor*SIN(rot_wg*pi/180.D0) &
                  *SIN(omega*time-pi*dph*(y-ymin_wg)/180.D0)
          END IF
       END DO
    END SELECT

  END SUBROUTINE phia_reflective

  !***********************************************************************
  SUBROUTINE profile(npmax,nxmax,nymax,x,y,vx,vy,vz,vpara,vperp,mass, &
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
    IF(npmax.EQ.0) THEN
       factor=DBLE(nxmax-1)*DBLE(nymax-1)
    ELSE
       factor=DBLE(nxmax-1)*DBLE(nymax-1)/DBLE(npmax)
    END IF

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
          IF(profiles(nx,ny,1).NE.0.D0) THEN
             profiles(nx,ny,2)=profiles(nx,ny,2)/profiles(nx,ny,1)
             profiles(nx,ny,3)=profiles(nx,ny,3)/profiles(nx,ny,1)
             profiles(nx,ny,4)=profiles(nx,ny,4)/profiles(nx,ny,1)
             profiles(nx,ny,5)=profiles(nx,ny,5)/profiles(nx,ny,1)
             profiles(nx,ny,6)=profiles(nx,ny,6)/profiles(nx,ny,1)
             profiles(nx,ny,7)=profiles(nx,ny,7)/profiles(nx,ny,1)
             profiles(nx,ny,8)=profiles(nx,ny,8)/profiles(nx,ny,1)
             profiles(nx,ny,9)=profiles(nx,ny,9)/profiles(nx,ny,1)
          END IF
       END DO
    END DO
  END SUBROUTINE profile

  !***********************************************************************
  SUBROUTINE profile_sub(nxmax,nymax,xp,yp,fp,fxy,imax,model_boundary)
    !***********************************************************************

    IMPLICIT NONE
    REAL(8) :: xp, yp, fp(imax)
    REAL(8), DIMENSION(0:nxmax,0:nymax,imax) :: fxy
    REAL(8) :: dx, dy, dz, sx2, sy2, sx2p, sy2p, sx2m, sy2m, dx1, dy1
    INTEGER :: nxmax, nymax, imax, model_boundary
    INTEGER :: nx, ny, nxp, nyp, nxpp, nxpm, nypp, nypm, nxppp, nyppp, i

    nxp = xp
    nyp = yp

    dx  = xp - DBLE(nxp)
    dy  = yp - DBLE(nyp)

    dx1 = 1.0d0 - dx
    dy1 = 1.0d0 - dy

    IF(dx .LE. 0.5d0) THEN
       sx2  = 3.0d0/4 - dx ** 2
       sx2p = 1.0d0/2 * (1.0d0/2 + dx) ** 2
       sx2m = 1.0d0/2 * (1.0d0/2 - dx) ** 2
    ELSE
       sx2  = 3.0d0/4 - (1.0d0 - dx) ** 2
       sx2p = 1.0d0/2 * (-1.0d0/2 + dx) ** 2
       sx2m = 1.0d0/2 * (3.0d0/2 - dx) ** 2
    ENDIF
    IF(dy .LE. 0.5d0) THEN
       sy2  = 3.0d0/4 - dy ** 2
       sy2p = 1.0d0/2 * (1.0d0/2 + dy) ** 2
       sy2m = 1.0d0/2 * (1.0d0/2 - dy) ** 2
    ELSE
       sy2  = 3.0d0/4 - (1.0d0 - dy) ** 2
       sy2p = 1.0d0/2 * (-1.0d0/2 + dy) ** 2
       sy2m = 1.0d0/2 * (3.0d0/2 - dy) ** 2
    ENDIF

    nxpm = nxp - 1
    nxpp = nxp + 1
    nxppp= nxp + 2
    nypm = nyp - 1
    nypp = nyp + 1
    nyppp= nyp + 2

    IF(model_boundary.EQ.0) THEN ! periodic
       IF( nxp .EQ. 0  ) nxpm = nxmax - 1
       IF( nyp .EQ. 0  ) nypm = nymax - 1
       IF( nxp .EQ. nxmax-1) nxppp=1
       IF( nyp .EQ. nymax-1) nyppp=1
    ELSE   ! reflective:
       IF( nxp .EQ. 0  ) nxpm = 0
       IF( nyp .EQ. 0  ) nypm = 0
       IF( nxp .EQ. nxmax-1) nxppp = nxmax
       IF( nyp .EQ. nymax-1) nyppp = nymax
    END IF

    IF(dx .LE. 0.5d0 .AND. dy .LE. 0.5d0) THEN
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
    ELSE IF(dx .LE. 0.5d0 .AND. dy .GE. 0.5d0) THEN
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
    ELSE IF(dx .GE. 0.5d0 .AND. dy .LE. 0.5d0) THEN
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
    ELSE
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
    ENDIF
  END SUBROUTINE profile_sub

  !***********************************************************************
  SUBROUTINE profile_boundary(nxmax,nymax,fxy,imax,model_boundary)
    !***********************************************************************

    IMPLICIT NONE
    INTEGER:: nxmax,nymax,imax,model_boundary
    REAL(8), DIMENSION(0:nxmax,0:nymax,imax) :: fxy
    INTEGER:: nx, ny, i

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
          DO ny = 1, nymax-1
             fxy(0,ny,i)     = 2.D0 * fxy(0,ny,i)
             fxy(nxmax,ny,i) = 2.D0 * fxy(nxmax,ny,i)
          END DO
       END DO
       DO i=1,imax
          DO nx = 1, nxmax-1
             fxy(nx,0,i)     = 2.D0 * fxy(nx,0,i)
             fxy(nx,nymax,i) = 2.D0 * fxy(nx,nymax,i)
          END DO
       END DO
       DO i=1,imax
          fxy(0,0,i)         = 4.D0 * fxy(0,0,i)
          fxy(0,nymax,i)     = 4.D0 * fxy(0,nymax,i)
          fxy(nxmax,0,i)     = 4.D0 * fxy(nxmax,0,i)
          fxy(nxmax,nymax,i) = 4.D0 * fxy(nxmax,nymax,i)
       END DO
    END IF

  END SUBROUTINE profile_boundary
END MODULE picexec
