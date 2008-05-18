!#####################################################################
      SUBROUTINE parab(psin)
!#####################################################################

  USE hfrmod, ONLY : a0, agg, aggps, ai, aii, aiips, aiota, aiott, b0, gg, loop,  &
                     psi, r0
  USE hfreya_all, ONLY : dspln
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: psin
  INTEGER(4) :: i, ii, iii
  INTEGER(4) :: ifst=0
  REAL(8) :: z
  REAL(8) :: y1(100),x1(100)
  REAL(8) :: wai(4,100),wgg(4,100),wio(4,100)

      IF(ifst.eq.0) THEN
        ai(0)=0.0
        CALL dspln(y1,x1,0,ai,psi,loop+1,wai,0)
        CALL dspln(y1,x1,0,gg,psi,loop+1,wgg,0)
        CALL dspln(y1,x1,0,aiota,psi,loop+1,wio,0)
        ifst=1
      ENDIF

      IF(psin.eq.psi(0))THEN
        ii=0
        GO TO 100
      ENDIF

      DO i=1,loop
        IF(psin.gt.psi(i-1).and.psin.le.psi(i)) THEN
          ii=i-1
          GO TO 100
        ENDIF
      END DO

      ii=loop
 100  CONTINUE

      z=psin-psi(ii)
      iii=ii+1
      aii=wai(1,iii)+z*(wai(2,iii)+z*(wai(3,iii)+z*wai(4,iii)))
      agg=wgg(1,iii)+z*(wgg(2,iii)+z*(wgg(3,iii)+z*wgg(4,iii)))
      aiott=wio(1,iii)+z*(wio(2,iii)+z*(wio(3,iii)+z*wio(4,iii)))
      aiips=wai(2,iii)+z*(2.*wai(3,iii)+3.d0*z*wai(4,iii))
      aggps=wgg(2,iii)+z*(2.*wgg(3,iii)+3.d0*z*wgg(4,iii))
      RETURN
      END


!#####################################################################
      SUBROUTINE intarz(psin)
!#####################################################################

  USE hfrmod, ONLY : abmnum, apnm, arnm, aznm, cp, cr, cz, dp, dr, dz, ep, er, ez, &
                     loop, modmax, pnmn, psi, rnmn, znmn
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: psin
  INTEGER(4):: i, ii, k
  REAL(8) ::  x1, xabm, xx, z

      IF(psin.eq.psi(0)) THEN
        ii=0
        GO TO 100
      ENDIF
      DO i=1,loop
        IF(psin.gt.psi(i-1).and.psin.le.psi(i)) THEN
          ii=i-1
          GO TO 100
        ENDIF
      END DO
      ii=loop
 100  CONTINUE
      z=psin-psi(ii)
      xx=z/psi(loop)
      x1=psin/psi(loop)
      DO k=1,modmax
        xabm=x1**abmnum(k)
        arnm(k)=rnmn(k,ii)+(cr(k,ii)+(dr(k,ii)+er(k,ii)*xx)*xx)*xx*xabm
        aznm(k)=znmn(k,ii)+(cz(k,ii)+(dz(k,ii)+ez(k,ii)*xx)*xx)*xx*xabm
        apnm(k)=pnmn(k,ii)+(cp(k,ii)+(dp(k,ii)+ep(k,ii)*xx)*xx)*xx*xabm
      END DO

      RETURN
      END
!
!#####################################################################
      SUBROUTINE derivb(psin,thin,fin)
!#####################################################################

  USE hfrmod, ONLY : abnm, bdfi, bdps, bdth, dabnm, mnumbr, modmax, nnumbr, sb
  IMPLICIT NONE
  REAL(8),INTENT(IN) :: psin, thin, fin
  INTEGER(4) :: iarg, k
  REAL(8) :: a, arg, cpi2, dc, ds


! the purpose of this subroutine is to calculate the quantities of
! sb, bdth, bdfi,  bdps, by interpolation from thier values at
! grid points.

!     data
!ppp
!cc      print * ,' derivb is started.'
!
      cpi2=dacos(-1.d0)*2.d0
      CALL intada(psin)
!ppp
!cc      print * ,' derivb is started.'

      sb=0.d0
      bdth=0.d0
      bdfi=0.d0
      bdps=0.d0
      DO k=2,modmax
        a=abnm(k)
        arg=(mnumbr(k)*thin-nnumbr(k)*fin)/cpi2
        iarg=arg
        arg=(arg-iarg)*cpi2
        dc=dcos(arg)
        ds=dsin(arg)
        sb=sb+a*dc
        bdth=bdth-mnumbr(k)*a*ds
        bdfi=bdfi+nnumbr(k)*a*ds
        bdps=bdps+dabnm(k)*dc
      END DO
      sb=sb*dsqrt(2.*psin)+abnm(1)
      bdps=bdps+dabnm(1)
!
      RETURN
      END

!#####################################################################
      SUBROUTINE intada(psin)
!#####################################################################

  USE hfrmod, ONLY : abmnum, abnm, bnmn, cb, dabnm, db, eb, loop, lxy, modmax, psi
  IMPLICIT NONE

  REAL(8),INTENT(IN) :: psin
  INTEGER(4) :: i, ii, k
  REAL(8) :: bharm, dbharm, dxabm, x1, xabm, xx, z
!c
!ppp
!cc      print * ,' intada is started.'

      IF(psin.eq.psi(0)) THEN
        ii=0
        z=0.d0
!ppp
!cc         print * ,' psin=psi .'
        GO TO 100
      ENDIF

!ppp
!cc         print * ,' psin neq. psi .'
      DO i=1,loop
        IF(psin.gt.psi(i-1).and.psin.le.psi(i)) THEN
          ii=i-1
          GO TO 100
        ENDIF
      END DO

      ii=loop
!ppp

 100  CONTINUE
      z=psin-psi(ii)
      x1=psin/psi(loop)
      xx=z/psi(loop)
!$$$         print * ,' psi(lp) =',psi(loop)
!$$$         print * ,' psin    =',psin
!$$$         print * ,' x1      =',x1
!$$$         print * ,' xx      =',xx
!$$$         print * ,' ii      =',ii

      IF(.not.lxy) THEN
!ppp
!ccc         print * ,' (not lxy) is true.'

        DO k=1,modmax
           xabm=x1**abmnum(k)
           dxabm=x1**(abmnum(k)-1)*abmnum(k)
           bharm=bnmn(k,ii)+(cb(k,ii)+(db(k,ii)+eb(k,ii)*xx)*xx)*xx
           dbharm=cb(k,ii)+(db(k,ii)*2.d0+eb(k,ii)*3.d0*xx)*xx
           abnm(k)=xabm*bharm
           dabnm(k)=(dxabm*bharm+xabm*dbharm)/psi(loop)
        END DO
!ppp
!cc         print * ,' ((not lxy) is true) is end.'
      ELSE
!
!ppp
!cc         print * ,' (not lxy) is not true.'

!     b=b/rho
!     bdps=bdps*rho,bdth=bdth/rho,bdfi=bdfi/rho
!
        DO k=2,modmax
           xabm=(x1+1.d-60)**(abmnum(k)-0.5d0)/dsqrt(2.d0*psi(loop))
           dxabm=(x1+1.d-60)**(abmnum(k)-0.5d0)*abmnum(k)*dsqrt(2.*psi(loop))
           bharm=bnmn(k,ii)+(cb(k,ii)+(db(k,ii)+eb(k,ii)*xx)*xx)*xx
           dbharm=cb(k,ii)+(db(k,ii)*2.d0+eb(k,ii)*3.d0*xx)*xx
           abnm(k)=xabm*bharm
           dabnm(k)=(dxabm*bharm+2.d0*psin*xabm*dbharm)/psi(loop)
         END DO

         bharm=bnmn(1,ii)+(cb(1,ii)+(db(1,ii)+eb(1,ii)*xx)*xx)*xx
         dbharm=cb(1,ii)+(db(1,ii)*2.d0+eb(1,ii)*3.d0*xx)*xx
         abnm(1)=bharm
         dabnm(1)=(2.d0*psin)**0.5d0*dbharm/psi(loop)
!ppp
!cc         print * ,' ((not lxy) is not true) is end.'
      ENDIF
      RETURN
      END
!
!#####################################################################
      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hdid,hmin,nok,nbad,iflg,nst,irunge)
!#####################################################################

  USE hfrmod, ONLY : cmax, dtime, dtx, dxsav, icmax, kmax, kount, xp, yp
  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: nvar, irunge
  INTEGER(4),INTENT(OUT):: nok, nbad, nst,iflg
  REAL(8),INTENT(IN) :: x1, x2, h1, hmin, eps
  REAL(8),INTENT(OUT):: hdid
  REAL(8),INTENT(INOUT):: ystart(nvar)

  INTEGER(4),PARAMETER :: maxstp=10000,nmax=10
  REAL(8)   ,PARAMETER :: two=2.d0,zero=0.d0,tiny=1.d-40
  INTEGER(4) :: i, icount, iorb, nstp
  REAL(8)    :: h, hnext, x, x0, xsav, zps, &
                dydx(nmax), y(nmax), y0(nmax), yscal(nmax)

!      hmax=h1
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      DO i=1,nvar
        y(i)=ystart(i)
      END DO
      xsav=x-dxsav*two
      DO nstp=1,maxstp
        nst=nstp
        CALL derivs(x,y,dydx)
        DO i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+tiny
!c    yscal(i)=1.0
        END DO
        IF(kmax.gt.0)THEN
          IF(abs(x-xsav).gt.abs(dxsav)) THEN
            IF(kount.lt.kmax-1)THEN
               kount=kount+1
               xp(kount)=x
               DO i=1,nvar
                 yp(i,kount)=y(i)
               END DO
               xsav=x
             ENDIF
           ENDIF
         ENDIF
         IF((x+h-x2)*(x+h-x1).gt.zero) h=x2-x
         zps=0.5*(y(2)**2+y(4)**2)
         icount=0
         x0=x
         DO i=1,nvar
           y0(i)=y(i)
         END DO
!
 21      h=dmin1(h,dtx/(100.*((zps/cmax)**4+0.02)))
         IF(irunge.eq.0) THEN
           CALL rkqc(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,iflg)
         ELSE
           CALL bsstep(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,iflg)
         ENDIF
         IF(iflg.lt.0) RETURN
         zps=0.5*(y(2)**2+y(4)**2)
         IF(zps.gt.cmax) THEN
           icount=icount+1
           IF(icount.le.icmax) THEN
             x=x0
             h=hdid/2.
             DO i=1,nvar
               y(i)=y0(i)
             END DO
             CALL derivs(x,y,dydx)
             DO i=1,nvar
               yscal(i)=abs(y(i))+abs(h*dydx(i))+tiny
             END DO
             GO TO 21
           ELSE
             iflg=0
             DO i=1,nvar
               ystart(i)=y(i)
             END DO
             RETURN
           ENDIF
         ENDIF

         iorb=zps/cmax*20+1.1
         dtime(iorb)=dtime(iorb)+hdid
         IF(hdid.eq.h)THEN
           nok=nok+1
         ELSE
           nbad=nbad+1
         ENDIF
         IF((x-x2)*(x2-x1).ge.zero)THEN
           DO i=1,nvar
             ystart(i)=y(i)
           END DO
           IF(kmax.ne.0)THEN
             kount=kount+1
             xp(kount)=x
             DO i=1,nvar
               yp(i,kount)=y(i)
             END DO
           ENDIF
           RETURN
         ENDIF
         h=hnext
!cc   if(h.lt.hmin) then
!cc   iflag=0
!cc   hdid=h
!cc   return
!cc   endif
      END DO
 17   iflg=0

      RETURN
      END


!#####################################################################
      SUBROUTINE rkqc(y,dydx,n,x,htry,eps,yscal,hdid,hnext,iflg)
!#####################################################################

  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: n
  INTEGER(4),INTENT(OUT):: iflg
  REAL(8),INTENT(IN) :: dydx(n),htry,eps, yscal(n)
  REAL(8),INTENT(OUT):: hdid, hnext
  REAL(8),INTENT(INOUT):: y(n), x
  INTEGER(4),PARAMETER :: nmax=10
  REAL(8),PARAMETER :: fcor=.0666666667d0, one=1.d0, safety=9.d-1, errcon=6.d-4
  INTEGER(4):: i
  REAL(8):: errmax, h, hh, pgrow, pshrnk, xsav, ytemp(nmax), ysav(nmax), dysav(nmax)

      pgrow=-0.20
      pshrnk=-0.25
      xsav=x
      DO i=1,n
        ysav(i)=y(i)
        dysav(i)=dydx(i)
      END DO
      h=htry
 1    hh=0.5*h
      CALL rk4(ysav,dysav,n,xsav,hh,ytemp)
      x=xsav+hh
      CALL derivs(x,ytemp,dydx)
      CALL rk4(ytemp,dydx,n,x,hh,y)
      x=xsav+h
      IF(x.eq.xsav) THEN
        hdid=h
        iflg=-1
        RETURN
      ENDIF
      CALL rk4(ysav,dysav,n,xsav,h,ytemp)
!c    call derivs(time,y,f)
!c    errmax=dabs(engy11-engy00)/engy00
      errmax=0.
      DO i=1,n
        ytemp(i)=y(i)-ytemp(i)
        errmax=max(errmax,abs(ytemp(i)/yscal(i)))
      END DO
      errmax=errmax/eps
      IF(errmax.gt.one) THEN
        h=safety*h*(errmax**pshrnk)
        GOTO 1
      ELSE
        hdid=h
        IF(errmax.gt.errcon)THEN
          hnext=safety*h*(errmax**pgrow)
        ELSE
          hnext=4.*h
        ENDIF
      ENDIF
      DO i=1,n
        y(i)=y(i)+ytemp(i)*fcor
      END DO
      RETURN
      END

!#####################################################################
      SUBROUTINE rk4(y,dydx,n,x,h,yout)
!#####################################################################

  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: n
  REAL,INTENT(IN) :: y(n), dydx(n), x, h
  REAL,INTENT(OUT):: yout(n)
  INTEGER(4), PARAMETER :: nmax=10
  INTEGER(4):: i
  REAL:: hh, h6, xh, yt(nmax), dyt(nmax), dym(nmax)

      hh=h*0.5
      h6=h/6.
      xh=x+hh
      DO i=1,n
        yt(i)=y(i)+hh*dydx(i)
      END DO
      CALL derivs(xh,yt,dyt)
      DO i=1,n
        yt(i)=y(i)+hh*dyt(i)
      END DO
      CALL derivs(xh,yt,dym)
      DO i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
      END DO
      CALL derivs(x+h,yt,dyt)
      DO i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
      END DO
      RETURN
      END


!#####################################################################
      SUBROUTINE derivs(time,y,f)
!#####################################################################

  USE hfrmod, ONLY : agg, aggps, aii, aiips, aiott, ambip, amu, bdfi, bdps, bdth,  &
                      cpsie, k1, lxy, sb, vpa
  IMPLICIT NONE
  REAL(8),INTENT(IN) :: y(4)
  REAL(8),INTENT(OUT):: f(4)
  REAL(8) :: time   !dummy
  REAL(8):: d, dphips, dpsdpf, dpsdpz, drodpf, drodpz, rmrb, roai1, rogio, rrb, zf6,&
            zps, zrho,zf4, zt


!     y(1)=fi,y(2)=theta,y(3)=rou(parallel),
!     y(4)=psi.all these quantities are demensionless variablesfor
!     every particles. y(1),y(2),y(3),y(4) are canonical variables.
!     f(1)=dfi/dt,f(2)=dthet/dt,
!     f(3)=drouparall/dt,f(4)=dpsi/dt.

      IF(.not.lxy) THEN
        CALL derivb(y(4),y(2),y(1))
        CALL parab(y(4))

        rogio=y(3)*aggps-aiott
        roai1=y(3)*aiips+1

        rmrb=amu(k1)+y(3)**2*sb
        rrb=y(3)*sb**2
        d=agg*roai1-aii*rogio

        drodpz=-rogio/d
        drodpf=roai1/d
        dpsdpz=agg/d
        dpsdpf=-aii/d
        dphips=-ambip/cpsie

        f(1) = (rmrb*bdps+dphips)*dpsdpf+rrb*drodpf
        f(2) = (rmrb*bdps+dphips)*dpsdpz+rrb*drodpz
        f(3) = -rmrb*(bdth*drodpz+bdfi*drodpf)
        f(4) = -rmrb*(bdth*dpsdpz+bdfi*dpsdpf)

        vpa(k1)=dsqrt((y(3)*sb)**2+2.d0*amu(k1)*sb)
!!         engy11=0.5d0*vpa(k1)**2
!
      ELSE
!
        zt=datan2(y(4),y(2))
        zps=0.5d0*(y(2)**2+y(4)**2)
        zrho=dsqrt(2.*zps)
!
        CALL derivb(zps,zt,y(1))
        CALL parab(zps)

        rogio=y(3)*aggps-aiott
        roai1=y(3)*aiips+1

        rmrb=amu(k1)+y(3)**2*sb
        rrb=y(3)*sb**2

        d=agg*roai1-aii*rogio
        drodpz=-rogio/d
        drodpf=roai1/d
        dpsdpz=agg/d
        dpsdpf=-aii/d
        dphips=-ambip/cpsie
!
        f(1)=(rmrb*bdps/(zrho+1.d-40)+dphips)*dpsdpf+rrb*drodpf


        f(3)=-rmrb*(bdth*drodpz+bdfi*drodpf)*zrho

        zf4=(rmrb*bdps+dphips*zrho)*dpsdpz+rrb*drodpz*zrho
        zf6=-rmrb*(bdth*dpsdpz+bdfi*dpsdpf)

        f(2)=dcos(zt)*zf6-dsin(zt)*zf4
        f(4)=dsin(zt)*zf6+dcos(zt)*zf4
!
        vpa(k1)=dsqrt((y(3)*sb)**2+2.d0*amu(k1)*sb)
!!          engy11=0.5d0*vpa(k1)**2
!
      ENDIF
      RETURN
      END

!
!#####################################################################
      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,iflg)
!#####################################################################

  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: nv
  INTEGER(4),INTENT(OUT):: iflg
  REAL(8),INTENT(IN) :: y(nv), dydx(nv), htry,eps,yscal(nv)
  REAL(8),INTENT(OUT):: hdid, hnext
  REAL(8),INTENT(INOUT):: x

  INTEGER(4),PARAMETER:: nmax=10, imax=11, nuse=7
  REAL(8),PARAMETER:: one=1.d0,shrink=.95d0,grow=1.2d0
  INTEGER(4):: i, j
  REAL(8):: errmax, h, xest, xsav, yerr(nmax), ysav(nmax),dysav(nmax),yseq(nmax)
  INTEGER(4):: nseq(imax)=(/2,4,6,8,12,16,24,32,48,64,96/)

      h=htry
      xsav=x
      DO i=1,nv
        ysav(i)=y(i)
        dysav(i)=dydx(i)
      END DO
 1    DO i=1,imax
        CALL mmid(ysav,dysav,nv,xsav,h,nseq(i),yseq)
        xest=(h/nseq(i))**2
        CALL rzextr(i,xest,yseq,y,yerr,nv,nuse)
        errmax=0.
        DO j=1,nv
          errmax=max(errmax,abs(yerr(j)/yscal(j)))
        END DO
        errmax=errmax/eps
        IF(errmax.lt.one) THEN
          x=x+h
          hdid=h
          IF(i.eq.nuse)THEN
            hnext=h*shrink
          ELSE IF(i.eq.nuse-1)THEN
            hnext=h*grow
          ELSE
            hnext=(h*nseq(nuse-1))/nseq(i)
          ENDIF
          RETURN
        ENDIF
      END DO
      h=0.25*h/2**((imax-nuse)/2)
      IF(x+h.eq.x) THEN
        iflg=-1
        RETURN
      ENDIF
      GOTO 1
      END
!
!
!#####################################################################
      SUBROUTINE rzextr(iest,xest,yest,yz,dy,nv,nuse)
!#####################################################################
!
  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: iest, nv,nuse
  REAL(8),INTENT(IN) :: xest, yest(nv)
  REAL(8),INTENT(OUT):: yz(nv), dy(nv)
  INTEGER(4),PARAMETER:: imax=11, nmax=10, ncol=7
  INTEGER(4):: j, k, m1
  REAL(8):: b, b1, c, ddy, v, x, yy, d(nmax,ncol), fx(ncol)

      x(iest)=xest
      IF(iest.eq.1) THEN
        DO j=1,nv
          yz(j)=yest(j)
          d(j,1)=yest(j)
          dy(j)=yest(j)
        END DO
      ELSE
        m1=min(iest,nuse)
        DO k=1,m1-1
          fx(k+1)=x(iest-k)/xest
        END DO
        DO j=1,nv
          yy=yest(j)
          v=d(j,1)
          c=yy
          d(j,1)=yy
          DO k=2,m1
            b1=fx(k)*v
            b=b1-c
            IF(b.ne.0.) THEN
              b=(c-v)/b
              ddy=c*b
              c=b1*b
            ELSE
              ddy=v
            ENDIF
            v=d(j,k)
            d(j,k)=ddy
            yy=yy+ddy
          END DO
          dy(j)=ddy
          yz(j)=yy
        END DO
      ENDIF
      RETURN
      END
!
!
!#####################################################################
      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv,nuse)
!#####################################################################

  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: iest, nv, nuse
  REAL(8),INTENT(IN) :: xest, yest(nv)
  REAL(8),INTENT(OUT):: yz(nv), dy(nv)
  INTEGER(4),PARAMETER:: imax=11, nmax=10, ncol=7
  INTEGER(4):: j, k1, m1
  REAL(8):: delta, f1, f2, q, x(imax), d(nmax), qcol(nmax,ncol)

      x(iest)=xest
      DO j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
      END DO
      IF(iest.eq.1) THEN
        DO j=1,nv
          qcol(j,1)=yest(j)
        END DO
      ELSE
        m1=min(iest,nuse)
        DO j=1,nv
          d(j)=yest(j)
        END DO
        DO k1=1,m1-1
          delta=1./(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          DO j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
          END DO
        END DO
!         do 16 j=1,nv
!            qcol(j,m1)=dy(j)
! 16      continue
      ENDIF
      RETURN
      END
!
!#####################################################################
      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout)
!#####################################################################

  IMPLICIT NONE
  INTEGER(4),INTENT(IN) ::nvar, nstep
  REAL(8),INTENT(IN) :: y(nvar), dydx(nvar), xs, htot
  REAL(8),INTENT(OUT):: yout(nvar)
  INTEGER(4),PARAMETER:: nmax=10
  INTEGER(4):: i, n
  REAL(8):: h, h2, swap, x, ym(nmax), yn(nmax)

      h=htot/nstep
      DO i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
      END DO
      x=xs+h
      CALL derivs(x,yn,yout)
      h2=2.*h
      DO n=2,nstep
        DO i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
        END DO
        x=x+h
        CALL derivs(x,yn,yout)
      END DO
      DO i=1,nvar
        yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
      END DO
      RETURN
      END


!#####################################################################
      SUBROUTINE coeff
!#####################################################################

  USE hfrmod, ONLY: abmnum, ai, aiota, bnm, bnmn, cb, cp, cpsie, cr, cz, db, &
    dp, dr, dz, eb, ep, er, ez, gg, icut, icut0, loop, lpi2, mm, mnumbr,     &
    modmax, nn, pnm, pnmn, psi, psino, rnm, rnmn, vol, vprime, znm, znmn
  USE hfreya_all,  ONLY: cubspl
  IMPLICIT NONE
  INTEGER(4):: i, ii, k, loop0, loop1, ibcbeg, ibcend
  REAL(8):: cpi2 , cbc(4,nn+1), crc(4,nn+1), czc(4,nn+1), cpc(4,nn+1)
  REAL(8):: fcent

      IF(lpi2.eq.1) THEN
        cpi2=2.d0*dacos(-1.d0)
        DO i=0,loop
          psi(i)=psi(i)/cpi2
        END DO
      ENDIF
!
      WRITE(7,601) (psi(i),i=0,loop)
 601  FORMAT(1h //,'  psi'/(1p10e12.3))
      cpsie=psi(loop)
      loop0=loop
      loop=(loop0-icut0)/icut
      WRITE(7,610) icut,loop0,loop
 610  FORMAT(1h //,'  icut,loop0,loop'/3i5)
      DO k=1,modmax
        DO i=1,loop
          ii=i*icut+icut0
          bnm(k,i)=bnm(k,ii)
          rnm(k,i)=rnm(k,ii)
          znm(k,i)=znm(k,ii)
          pnm(k,i)=pnm(k,ii)
        END DO
      END DO
!
      DO k=1,modmax
        abmnum(k)=0.0d0
!c        abmnum(k)=0.5d0*abs(mnumbr(k))
      END DO
      psino(0)=0.0
      DO i=1,loop
        ii=i*icut+icut0
        psi(i)=psi(ii)
        vol(i)=vol(ii)
        vprime(i)=vprime(ii)
        gg(i)=gg(ii)
        ai(i)=ai(ii)
        aiota(i)=aiota(ii)
      END DO
!      cpie=psi(loop)
      DO i=0,loop
        psino(i)=psi(i)/cpsie
      END DO
      WRITE(7,602) (psi(i),i=0,loop)
 602  FORMAT(1h //,'  psi'/(1p10e12.3))
      DO k=1,modmax
        IF(mnumbr(k).eq.0.0) THEN
          bnm(k,0)=fcent(bnm,psino,k,mm)
          rnm(k,0)=fcent(rnm,psino,k,mm)
        ENDIF
      END DO
      DO k=1,modmax
        DO i=1,loop
          bnmn(k,i)=bnm(k,i)/(psino(i)**abmnum(k))
          rnmn(k,i)=rnm(k,i)/(psino(i)**abmnum(k))
          znmn(k,i)=znm(k,i)/(psino(i)**abmnum(k))
          pnmn(k,i)=pnm(k,i)/(psino(i)**abmnum(k))
        END DO
        bnmn(k,0)=fcent(bnmn,psino,k,mm)
        rnmn(k,0)=fcent(rnmn,psino,k,mm)
        znmn(k,0)=fcent(znmn,psino,k,mm)
        pnmn(k,0)=fcent(pnmn,psino,k,mm)
!

        DO i=0,loop
          cbc(1,i+1)=bnmn(k,i)
          crc(1,i+1)=rnmn(k,i)
          czc(1,i+1)=znmn(k,i)
          cpc(1,i+1)=pnmn(k,i)
        END DO
        loop1=loop+1
        ibcbeg=1
        ibcend=0
        cbc(2,1)=0.0
        crc(2,1)=0.0
        czc(2,1)=0.0
        cpc(2,1)=0.0
        CALL cubspl(psino,cbc,loop1,ibcbeg,ibcend)
        CALL cubspl(psino,crc,loop1,ibcbeg,ibcend)
        CALL cubspl(psino,czc,loop1,ibcbeg,ibcend)
        CALL cubspl(psino,cpc,loop1,ibcbeg,ibcend)
        DO i=0,loop
          cb(k,i)=cbc(2,i+1)
          cr(k,i)=crc(2,i+1)
          cz(k,i)=czc(2,i+1)
          cp(k,i)=cpc(2,i+1)
          db(k,i)=cbc(3,i+1)/2.d0
          dr(k,i)=crc(3,i+1)/2.d0
          dz(k,i)=czc(3,i+1)/2.d0
          dp(k,i)=cpc(3,i+1)/2.d0
          eb(k,i)=cbc(4,i+1)/6.d0
          er(k,i)=crc(4,i+1)/6.d0
          ez(k,i)=czc(4,i+1)/6.d0
          ep(k,i)=cpc(4,i+1)/6.d0
        END DO
      END DO
      RETURN
      END

!#####################################################################
      REAL(8) FUNCTION fcent(z,psii,k,m)
!#####################################################################

  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: k,m
  REAL(8),INTENT(IN):: psii(0:*),z(m,0:*)

      fcent= &
         (psii(2)*psii(3))/(psii(1)-psii(2))/(psii(1)-psii(3))*z(k,1) &
        +(psii(3)*psii(1))/(psii(2)-psii(3))/(psii(2)-psii(1))*z(k,2) &
        +(psii(1)*psii(2))/(psii(3)-psii(1))/(psii(3)-psii(2))*z(k,3)
!
!c    &  (psii(2)**2*z(k,1)-psii(1)**2*z(k,2))/(psii(2)**2-psii(1)**2)

      RETURN
      END


!!!#####################################################################
!!      function fcent0(z,psii)
!!!#####################################################################
!!
!!      implicit real*8(a-h,o-z)
!!      dimension psii(0:*),z(0:*)
!!
!!      fcent0= &
!!     &     (psii(2)*psii(3))/(psii(1)-psii(2))/(psii(1)-psii(3))*z(1) &
!!     &     +(psii(3)*psii(1))/(psii(2)-psii(3))/(psii(2)-psii(1))*z(2) &
!!     &     +(psii(1)*psii(2))/(psii(3)-psii(1))/(psii(3)-psii(2))*z(3)
!!!
!!!c   &  (psii(2)**2*z(k,1)-psii(1)**2*z(k,2))/(psii(2)**2-psii(1)**2)
!!      return
!!      end



!!!#####################################################################
!!      subroutine absmax(a,km,ka,loop,bmax,md)
!!!#####################################################################
!!
!!     implicit real*8 (a-h,o-z)
!!    dimension a(md,*)
!!
!!      bmax=-1.e70
!!      do 100 k=km,km+ka-1
!!         do 100 i=1,loop
!!            bmax=dmax1(dabs(a(k,i)),bmax)
!! 100  continue
!!      return
!!      end


!#####################################################################
      SUBROUTINE coeff1
!#####################################################################

  USE hfrmod, ONLY: abmnum, ai, aiota, bnm, bnmn, cb, cp, cpsie, cr, cz, db, dp,  &
    dr, dz, eb, ep, er, ez, gg, icut, loop, modmax, nn, pnm, pnmn, psi, psino,    &
    rnm, rnmn, vol, vprime, znm, znmn
  USE hfreya_all, ONLY : dspln1
  IMPLICIT NONE
  INTEGER(4)::i, ii, k,  loop0, loop1
  REAL(8):: x0(101), y0(101), x(101), y(101),   &
    cbc(4,nn+1), cpc(4,nn+1), crc(4,nn+1), czc(4,nn+1)


      loop0=loop
      cpsie=psi(loop0)
      loop=loop0/icut
      WRITE(6,666) icut,loop0,loop
 666  FORMAT(1h //,'  icut,loop0,loop'/(3i10))
      DO k=1,modmax
        DO i=0,loop
          ii=i*icut
          bnm(k,i)=bnm(k,ii)
          rnm(k,i)=rnm(k,ii)
          znm(k,i)=znm(k,ii)
          pnm(k,i)=pnm(k,ii)
        END DO
      END DO

      DO k=1,modmax
        abmnum(k)=0.0d0
!c         abmnum(k)=0.5*abs(mnumbr(k))
      END DO
      psino(0)=0.0d0
      DO i=1,loop
        ii=icut*i
        psi(i)=psi(ii)
        vol(i)=vol(ii)
        vprime(i)=vprime(ii)
        gg(i)=gg(ii)
        ai(i)=ai(ii)
        aiota(i)=aiota(ii)
        psino(i)=psi(i)/cpsie
      END DO

      DO k=1,modmax
        DO i=1,loop
          bnmn(k,i)=bnm(k,i)/(psino(i)**abmnum(k))
          rnmn(k,i)=rnm(k,i)/(psino(i)**abmnum(k))
          znmn(k,i)=znm(k,i)/(psino(i)**abmnum(k))
          pnmn(k,i)=pnm(k,i)/(psino(i)**abmnum(k))
        END DO

!         bnmn(k,0)=fcent(bnmn,psino,k,mm)
!         rnmn(k,0)=fcent(rnmn,psino,k,mm)
!         znmn(k,0)=fcent(znmn,psino,k,mm)
!         pnmn(k,0)=fcent(pnmn,psino,k,mm)

        bnmn(k,0)=bnm(k,0)
        rnmn(k,0)=rnm(k,0)
        znmn(k,0)=znm(k,0)
        pnmn(k,0)=pnm(k,0)

!
        DO i=0,loop
          y0(i+1)=bnmn(k,i)
          x0(i+1)=psino(i)
        END DO
        loop1=loop+1
        CALL dspln1(y,x,0,y0,x0,loop1,cbc,0)
        DO i=0,loop
          y0(i+1)=rnmn(k,i)
        END DO
        CALL dspln1(y,x,0,y0,x0,loop1,crc,0)
        DO i=0,loop
          y0(i+1)=znmn(k,i)
        END DO
        CALL dspln1(y,x,0,y0,x0,loop1,czc,0)
        DO i=0,loop
          y0(i+1)=pnmn(k,i)
        END DO
        CALL dspln1(y,x,0,y0,x0,loop1,cpc,0)
        DO i=0,loop
          cb(k,i)=cbc(2,i+1)
          db(k,i)=cbc(3,i+1)
          eb(k,i)=cbc(4,i+1)
          cr(k,i)=crc(2,i+1)
          dr(k,i)=crc(3,i+1)
          er(k,i)=crc(4,i+1)
          cz(k,i)=czc(2,i+1)
          dz(k,i)=czc(3,i+1)
          ez(k,i)=czc(4,i+1)
          cp(k,i)=cpc(2,i+1)
          dp(k,i)=cpc(3,i+1)
          ep(k,i)=cpc(4,i+1)
        END DO
      END DO
!c    do 900 k=1,modmax
!c    write(6,901)(cb(k,j),j=0,loop)
 901  format(1h /,' cb=',1p8e13.3)
!     write(6,902)(cr(k,j),j=0,loop)
 902  format(1h /,' cr=',1p8e13.3)
!     write(6,903)(cz(k,j),j=0,loop)
 903  format(1h /,' cz=',1p8e13.3)
!     write(6,904)(cp(k,j),j=0,loop)
 904  format(1h, /' cp=',1p8e13.3)
 900  CONTINUE
      RETURN
      END

!
!!!***********************************************************************
!!      subroutine sortd1 ( a, l, n, ifg )
!!!***********************************************************************
!!!      sort a vector of real number
!!
!!!     parameter : a    (i) sort suru data
!!!                      (o) sort gono data
!!!                 l    (i) 1 kara n madeno suu
!!!                      (o) a o sort sita jyunban
!!!                 n    (i) data-suu
!!!                 ifg  (i) sort-flg    ifg.eq.1 : up    ifg.ne.1 : down
!!!***********************************************************************
!!      implicit    real*8 (a-h,o-z)
!!      dimension   a(n), l(n)
!!      logical     lsort
!!
!!      ng = n
!!      do 40 while ( ng.gt.1 )
!!         ng = ng/2
!!         imax = n - ng
!!         do 30 while ( k.ne.0 )
!!            k = 0
!!            do 20 i=1,imax
!!               ip = i + ng
!!               if ( ifg.eq.1 ) then
!!                  lsort = a(i).gt.a(ip)
!!               else
!!                  lsort = a(i).lt.a(ip)
!!               endif
!!               if ( lsort ) then
!!                  aw    = a(i)
!!                  a(i)  = a(ip)
!!                  a(ip) = aw
!!                  lw    = l(i)
!!                  l(i)  = l(ip)
!!                  l(ip) = lw
!!                  k     = k + 1
!!               endif
!!   20       continue
!!   30    continue
!!   40 continue
!!
!!      return
!!      end
