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
      nxmaxh1 = nxmax/ 2 + 1
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
      ex(:,:) = 0.d0
      ey(:,:) = 0.d0
      ez(:,:) = 0.d0
      exb(:,:) = 0.d0
      eyb(:,:) = 0.d0
      ezb(:,:) = 0.d0
      exbb(:,:) = 0.d0
      eybb(:,:) = 0.d0
      ezbb(:,:) = 0.d0
      esx(:,:) = 0.d0
      esy(:,:) = 0.d0
      esz(:,:) = 0.d0
      emx(:,:) = 0.d0
      emy(:,:) = 0.d0
      emz(:,:) = 0.d0
      bx(:,:) = 0.d0
      by(:,:) = 0.d0
      bz(:,:) = 0.d0
      bxb(:,:) = 0.d0
      byb(:,:) = 0.d0
      bzb(:,:) = 0.d0
      bxbb(:,:) = 0.d0
      bybb(:,:) = 0.d0
      bzbb(:,:) = 0.d0
      Axbb(:,:) = 0.d0
      Axb(:,:) = 0.d0
      Aybb(:,:) = 0.d0
      Ayb(:,:) = 0.d0
      Azbb(:,:) = 0.d0
      Azb(:,:) = 0.d0
      Axb(:,:) = 0.d0
      Ax(:,:) = 0.d0
      Ayb(:,:) = 0.d0
      Ay(:,:) = 0.d0
      Azb(:,:) = 0.d0
      Az(:,:) = 0.d0
      phib(:,:) = 0.d0
      phi(:,:) = 0.d0

      !..... set initial positions and velocities of electrons
      call iniset(npmax,npxmax,npymax,nxmax,nymax,densx, &
                  xe,ye,ze,xeb,yeb,zeb,xemid,yemid,vxe,vye,vze,vte,dt,iran,&
                  x1,x2,y1,y2,alx,aly,model_boundary,vzone)

      !..... set initial positions and velocities of ions
      call iniset(npmax,npxmax,npymax,nxmax,nymax,densx, &
                  xi,yi,zi,xib,yib,zib,ximid,yimid,vxi,vyi,vzi,vti,dt,iran,&
                  x1,x2,y1,y2,alx,aly,model_boundary,vzone)
      
      !..... initialize scalar potential by poisson solver
      ipssn = 0
       IF(model_boundary.EQ.0) THEN
          call poisson_f(nxmax,nymax,nxmaxh1,nxmax1,nymax1, &
                         rho,phi,rhof,phif,awk,afwk,cform,ipssn)
       END IF

      !..... initialize wall clock time
      call mtx_barrier
      wtime1 = mpi_wtime()

      DO nx=vzone,nxmax-vzone
         DO ny=vzone,nymax-vzone
            factor=DBLE(nx)/DBLE(nxmax)
            bxbg(nx,ny)=bxmin+(bxmax-bxmin)*factor
            bybg(nx,ny)=bymin+(bymax-bymin)*factor
            bzbg(nx,ny)=bzmin+(bzmax-bzmin)*factor
         END DO
      END DO
       !.......... calculate ex and ey and ez
         CALL efield(nxmax,nymax,dt,phi,Ax,Ay,Az,Axb,Ayb,Azb, &
                     ex,ey,ez,exb,eyb,ezb,exbb,eybb,ezbb,bxb,byb,bzb,&
                     esx,esy,esz,emx,emy,emz,jx,jy,jz,vcfact,model_push,model_boundary)
         !.......... calculate bx and by and bz
         CALL bfield(nxmax,nymax,dt,Ax,Ay,Az,Axb,Ayb,Azb,ex,ey,ez,&
                     bx,by,bz,bxb,byb,bzb,bxbb,bybb,bzbb,bxbg,bybg,bzbg,bb, &
                     vcfact,model_push,model_boundary)
      do np=1,npmax
         vparae(np)=vye(np)
         vperpe(np)=SQRT(vxe(np)**2+vze(np)**2)
         vparai(np)=vyi(np)
         vperpi(np)=SQRT(vxi(np)**2+vzi(np)**2)
      end do

      call kine(npmax,vxe,vye,vze,akine0,me,vcfact)
      call kine(npmax,vxi,vyi,vzi,akini0,mi,vcfact)
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
                        x,y,z,xb,yb,zb,xmid,ymid,vx,vy,vz,vt,dt,iran,&
                        x1,x2,y1,y2,alx,aly,model_boundary,vzone)
!***********************************************************************
      implicit none
      INTEGER,INTENT(IN) :: &
           npmax, npxmax, npymax, nxmax, nymax, iran, vzone, model_boundary
      REAL(8),DIMENSION(npmax),INTENT(INOUT) :: &
           x, y, z, xb, yb, zb, vx, vy, vz, xmid, ymid
      REAL(8),INTENT(IN) :: densx,dt, vt,x1,x2,y1,y2,alx,aly
      REAL(8):: factx, facty, rvx, rvy, rvz, inter, position,&
           x3,x4,y3,y4,alx1,aly1,vdzone,&
           xl_before,xl_after,yl_before,yl_after
      integer :: npx, npy, np

      factx = dble(nxmax) / dble(npxmax)
      facty = dble(nymax) / dble(npymax)
      vdzone = dble(vzone)
      IF(model_boundary .ne. 0) THEN
        factx = dble(nxmax-2*vzone) / dble(npxmax)
        facty = dble(nymax-2*vzone) / dble(npymax)
        x3 = x1 + vdzone
        y3 = y1 + vdzone
        x4 = x2 - vdzone
        y4 = y2 - vdzone
        alx1 = alx - vdzone
        aly1 = aly - vdzone
      ENDIF
      np = 0
      IF(densx .lt. 0.d0) then ! subroutine for uniform density
      DO npy = 1, npymax
      DO npx = 1, npxmax
        np = np + 1
         x(np) = (dble(npx) - 0.5d0 ) * factx + vdzone
         y(np) = (dble(npy) - 0.5d0 ) * facty + vdzone
         call gauss(rvx,rvy,rvz,iran)
         vx(np) = rvx * vt
         vy(np) = rvy * vt
         vz(np) = rvz * vt
         xb(np) = x(np) - vx(np) * dt
         yb(np) = y(np) - vy(np) * dt
         zb(np) = z(np) - vz(np) * dt
        IF(model_boundary .eq. 0) THEN
           IF( xb(np) .LT. x1 ) THEN
              DO WHILE(xb(np) .LT. x1)
                xb(np) = xb(np) + alx
              END DO
           ELSEIF( xb(np) .GT. x2 ) THEN
              DO WHILE(xb(np) .GT. x2)
                xb(np) = xb(np) - alx
              END DO
           END IF
           xmid(np) = 0.5d0*(xb(np)+x(np))
           IF( yb(np) .LT. y1 ) THEN
              DO WHILE(yb(np) .LT. y1)
                yb(np) = yb(np) + aly
              END DO
           ELSEIF( yb(np) .GT. y2 ) THEN
              DO WHILE(yb(np) .GT. 2)
                yb(np) = yb(np) - aly
              END DO
           ENDIF
           ymid(np) = 0.5d0*(yb(np)+y(np))
         ELSE
           xmid(np)=0.5D0*(xb(np)+x(np))
           IF( x(np) .LT. x3  ) THEN
              xl_before=xb(np) - x3
              x(np) = x3 + (x3  - x(np))
              xl_after=x(np) - x3
              xmid(np)=x3+0.5D0*ABS(xl_after-xl_before)
              vx(np) = -vx(np)
           ELSEIF( x(np) .GT. x4 ) THEN
              xl_before=x4 - xb(np)
              x(np) = x4 - (x(np) - x4)
              xl_after=x4 - x(np)
              xmid(np)=x4-0.5D0*ABS(xl_after-xl_before)
              vx(np) = -vx(np)
           ENDIF
           ymid(np)=0.5D0*(yb(np)+y(np))
           IF( y(np) .LT. y3 ) THEN
             yl_before=yb(np) - y3
             y(np) = y3 + (y3  - y(np))
             yl_after=y(np) - y3
             ymid(np)=y3+0.5D0*ABS(yl_after-yl_before)
             vy(np) = -vy(np)
           ELSEIF( y(np) .GT. y4 ) THEN
             yl_before=y4 - yb(np)
             y(np) = y4 - (y(np) - y4)
             yl_after=y4 - y(np)
             ymid(np)=y4-0.5D0*ABS(yl_after-yl_before)
             vy(np) = -vy(np)
           ENDIF
        ENDIF                                             
      END DO
      END DO
   ELSE ! subroutine for density gradient
      inter = 0.d0
      DO npx = 1 , npxmax+1
      inter = inter + 1.d0/(dble(npx)*(1.d0/(1.d0-densx)-1.d0)/(dble(npxmax-1))&
                           +(dble(npxmax)-1.d0/(1.d0-densx))/dble(npxmax-1))
      END DO
       ! inter = dble(nxmax-2*vdzone)/((dble(npxmax)+1.0d0)*(1.0d0-0.5d0*densx))
      !inter = dble(nxmax-2*vdzone) / inter
      inter = 0.5d0*(1.d0+densx)*(1.d0-densx)
      DO npy = 1, npymax
         position = 0.d0
      DO npx = 1, npxmax
         np = np + 1
         !position = position &
         !         + inter/(dble(npx)*(1.d0/(1.d0-densx)-1.d0)/(dble(npxmax-1))&
         !                  +(dble(npxmax)-1.d0/(1.d0-densx))/dble(npxmax-1))
         !x(np) = position + vdzone
         x(np) =dble(nxmax-2*vzone)*(sqrt((1.d0-0.5*densx)**2+2.d0*densx*dble(npx)/dble(npxmax+1))-1.d0+0.5*densx)/densx+dble(vzone)
         y(np) = (dble(npy) - 0.5d0 ) * facty + vdzone

         call gauss(rvx,rvy,rvz,iran)
         vx(np) = rvx * vt
         vy(np) = rvy * vt
         vz(np) = rvz * vt
         xb(np) = x(np) - vx(np) * dt
         yb(np) = y(np) - vy(np) * dt
         zb(np) = z(np) - vz(np) * dt
         ! particle reflect condition on boundary
        IF(model_boundary .eq. 0) THEN
           IF( xb(np) .LT. x1 ) THEN
              DO WHILE(xb(np) .LT. x1)
                xb(np) = xb(np) + alx
              END DO
           ELSEIF( xb(np) .GT. x2 ) THEN
              DO WHILE(xb(np) .GT. x2)
                xb(np) = xb(np) - alx
              END DO
           END IF
           xmid(np) = 0.5d0*(xb(np)+x(np))
           IF( yb(np) .LT. y1 ) THEN
              DO WHILE(yb(np) .LT. y1)
                yb(np) = yb(np) + aly
              END DO
           ELSEIF( yb(np) .GT. y2 ) THEN
              DO WHILE(yb(np) .GT. 2)
                yb(np) = yb(np) - aly
              END DO
           ENDIF
           ymid(np) = 0.5d0*(yb(np)+y(np))
         ELSE
           xmid(np)=0.5D0*(xb(np)+x(np))
           IF( x(np) .LT. x3  ) THEN
             xl_before=xb(np) - x3
             x(np) = x3 + (x3  - x(np))
             xl_after=x(np) - x3
             xmid(np)=x3+0.5D0*ABS(xl_after-xl_before)
             vx(np) = -vx(np)
           ELSEIF( x(np) .GT. x4 ) THEN
             xl_before=x4 - xb(np)
             x(np) = x4 - (x(np) - x4)
             xl_after=x4 - x(np)
             xmid(np)=x4-0.5D0*ABS(xl_after-xl_before)
             vx(np) = -vx(np)
           ENDIF
           ymid(np)=0.5D0*(yb(np)+y(np))
           IF( y(np) .LT. y3 ) THEN
             yl_before=yb(np) - y3
             y(np) = y3 + (y3  - y(np))
             yl_after=y(np) - y3
             ymid(np)=y3+0.5D0*ABS(yl_after-yl_before)
             vy(np) = -vy(np)
           ELSEIF( y(np) .GT. y4 ) THEN
             yl_before=y4 - yb(np)
             y(np) = y4 - (y(np) - y4)
             yl_after=y4 - y(np)
             ymid(np)=y4-0.5D0*ABS(yl_after-yl_before)
             vy(np) = -vy(np)
           ENDIF
        ENDIF
      END DO
      END DO
      END IF
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
