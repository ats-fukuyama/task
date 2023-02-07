c     
      module tpxssl_mod
      public
      contains
c=======================================================================
      subroutine deriv(dy,y,x,n)
c=======================================================================
      implicit none
! argument
      integer   n
      real*8    dy(n), y(n), x(n)
! local variables
      integer   m, nm
      real*8    dx0, dx1, dx2, dy0, dy1, dy2
c=======================================================================
      nm=n-1
      dx1=x(2)-x(1)
      dx2=x(3)-x(2)
      dy1=(y(2)-y(1))/dx1
      dy2=(y(3)-y(2))/dx2
      dy(1)=(-dx1*dy2+(dx2+2.*dx1)*dy1)/(dx2+dx1)
      do m=2,nm
      dx0=dx1
      dy0=dy1
      dx1=x(m+1)-x(m)
      dy1=(y(m+1)-y(m))/dx1
      dy(m)=(dx0*dy1+dx1*dy0)/(dx1+dx0)
      enddo
      dy(n)=((2.*dx1+dx0)*dy1-dx1*dy0)/(dx1+dx0)
      return
      end subroutine deriv
c
c=======================================================================
      subroutine spln(y,x,n,y0,x0,n0,w,i)
c=======================================================================
      implicit none
! argument
      integer   n, n0, i
      real*8    y(n), x(n), y0(n0), x0(n0), w(n0,4)
! local variables
      integer   j, k, m, nm
      real*8    dd0, dd1, dx0, dx1, dx2, dy0, dy1, dy2, z
c=======================================================================
      if(i.eq.0)then
      nm=n0-2
      dx1=x0(2)-x0(1)
      dx2=x0(3)-x0(2)
      dy1=(y0(2)-y0(1))/dx1
      dy2=(y0(3)-y0(2))/dx2
      dd1=(-dx1*dy2+(dx2+2.*dx1)*dy1)/(dx2+dx1)
      do m=1,nm
      dx0=dx1
      dy0=dy1
      dd0=dd1
      dx1=x0(m+2)-x0(m+1)
      dy1=(y0(m+2)-y0(m+1))/dx1
      dd1=(dx0*dy1+dx1*dy0)/(dx1+dx0)
      w(m,1)=y0(m)
      w(m,2)=dd0
      w(m,3)=(3.*dy0-dd1-2.*dd0)/dx0
      w(m,4)=(-2.*dy0+dd1+dd0)/(dx0*dx0)
      enddo
      dd0=dd1
      dd1=((2.*dx1+dx0)*dy1-dx1*dy0)/(dx1+dx0)
      m=n0-1
      w(m,1)=y0(m)
      w(m,2)=dd0
      w(m,3)=(3.*dy1-dd1-2.*dd0)/dx1
      w(m,4)=(-2.*dy1+dd1+dd0)/(dx1*dx1)
      w(n0,1)=y0(n0)
      if(n.eq.0)return
      endif
c-----------------------------------------------------------------------
      if(x(1).lt.x0( 1))x(1)=x0( 1)
      if(x(n).gt.x0(n0))x(n)=x0(n0)
c-----
      j=2
      do m=1,n
      do while(x(m).gt.x0(j))
      j=j+1
      if(j.ge.n0)then
      j=n0
      endif
      enddo
      k=j-1
      z=x(m)-x0(k)
      y(m)=w(k,1)+z*(w(k,2)+z*(w(k,3)+z*w(k,4)))
      enddo
      return
      end subroutine spln
c
c=======================================================================
      subroutine prnts(ft06,titl,x,n,x0)
c=======================================================================
      implicit none
! argument
      integer   ft06, n
      character titl*4
      real*8    x(n), x0
! local variables
      integer   i, i0, i1, iy, k, l
      real*8    xmax, y(100), y0, ymax, yymax
c=======================================================================
c-----NORMALIZATION FACTOR
      xmax=0.
      do i=1,n
      if(dabs(x(i)).gt.xmax)xmax=dabs(x(i))
      enddo
c-----
      if(xmax.eq.0.)then
      write(ft06,'(a3,a4,1pd11.4)')'***',titl,xmax
      return
      endif
c-----
      ymax=xmax
      y0=x0
      if(dabs(y0).le.1.d-20)y0=1.
      if(ymax.ne.0.)then
      if(dabs(y0).gt.1.d-20)ymax=ymax/dabs(y0)
      if(ymax.ge.10..or.ymax.lt.0.999999)then
      yymax=0.
      if(xmax.gt.1.d-20)yymax=dlog10(xmax)
      iy=0
      if(yymax.gt.0.)iy=yymax
      if(yymax.lt.0.)iy=yymax-1.
      y0=1.
      if(iy.gt.0)y0=10.**iy
      if(iy.lt.0)y0=1./10.**(-iy)
      endif
      endif
c-----------------------------------------------------------------------
c-----PRINT OUT
      i0=1
      i1=11
      if(n.le.i1)i1=n
      if(n.gt.i1)i1=i1-1
      k=0
      do i=i0,i1
      k=k+1
      y(k)=x(i)/y0
      enddo 
      do while(i1.lt.n)
      i0=i1+1
      i1=i1+11
      if(n.le.i1)i1=n
      if(n.gt.i1)i1=i1-1
      k=0
      do i=i0,i1
      k=k+1
      y(k)=x(i)/y0
      enddo
c      write(ft06,'(19x,11f10.4)')(y(l),l=1,k)
      write(ft06,'(10f8.4)')(y(l),l=1,k)
      enddo
      return
      end subroutine prnts
c
c===================================================================
      subroutine matslv(mmnn,m,n,b,c,err,ill)
c===================================================================
c     SOLVE B*X=D
c===================================================================
      implicit none
! argument
      integer   mmnn, m, n, ill
      real*8    b(mmnn,mmnn), c(mmnn,mmnn), err
! local variables
      integer   i, j, k, l, nn, nnn
      real*8    x
c===================================================================
c     LU DECOMPOSITION
c===================================================================
      do j=1,m
c-------------------------------------------------------------------
c::B(J,K) : K <= J
c-------------------------------------------------------------------
      do k=1,j
      x=0.
      if(k.gt.1)then
      do l=1,k-1
      x=x+b(j,l)*b(l,k)
      enddo
      endif
c-----
      b(j,k)=b(j,k)-x
      enddo
c-------------------------------------------------------------------
c::B(J,K) : K > J
c-------------------------------------------------------------------
      if(j.lt.m)then
      do k=j+1,m
      x=0.
c-----
      do l=1,j-1
      x=x+b(j,l)*b(l,k)
      enddo
c-----
      b(j,k)=(b(j,k)-x)/b(j,j)
      enddo
      endif
c-----------------------------------------------------------------
      enddo
c=================================================================
      if(n.lt.m)return
c=================================================================
      if(n.gt.m)then
      nn=n-m
      do nnn=1,nn
      do j=1,m
      x=0.
c-----
      if(j.gt.1)then
      do k=1,j-1
      x=x+b(j,k)*b(k,m+nnn)
      enddo
      endif
c-----
      b(j,m+nnn)=(b(j,m+nnn)-x)/b(j,j)
      enddo
c-------------------------------------------------------------------
      do j=m,1,-1
      x=0.
c-----
      if(j.lt.m)then
      do k=j+1,m
      x=x+b(j,k)*b(k,m+nnn)
      enddo
      endif
c-----
      b(j,m+nnn)=b(j,m+nnn)-x
      enddo
      enddo
c-----
      ill=0
      err=0.
      return
      endif
c=====================================================================
c  INVERSE MATRIX : C(I,J)
c=====================================================================
      do i=1,m
      do j=1,m
      c(i,j)=0.
      enddo
      c(i,i)=1.
      enddo
c-----
      do nnn=1,m
      do j=1,m
      x=0.
c-----
      if(j.gt.1)then
      do k=1,j-1
      x=x+b(j,k)*c(k,nnn)
      enddo
      endif
c-----
      c(j,nnn)=(c(j,nnn)-x)/b(j,j)
      enddo
c-------------------------------------------------------------------
      do j=m,1,-1
      x=0.
c-----
      if(j.lt.m)then
      do k=j+1,m
      x=x+b(j,k)*c(k,nnn)
      enddo
      endif
c-----
      c(j,nnn)=c(j,nnn)-x
      enddo
      enddo
c=====================================================================
      ill=0
      err=0.
      return
      end subroutine matslv
c
c=====================================================================
      subroutine cpusec(cputim)
c=====================================================================
      implicit none
      real*8    cputim
      return
      end subroutine cpusec
c
c=======================================================================
      double precision function   random(i)
c=======================================================================
      implicit none
! argument
      integer   i
! local variables
      integer   icover, ix
      real*8    rlamda, rmod, w
c
      data ix /1234567/
c=======================================================================
      rmod=2147483648.d0
      rlamda=65539.d0
      if(ix) 2,4,6
   2  ix=-ix
      go to 6
   4  ix=3907
   6  w=ix
      w=rlamda*w
      if(w-rmod) 20,10,10
  10  icover=w/rmod
      w=w-float(icover)*rmod
  20  ix=w
      random=w/rmod
      return
      end function random
c
c=======================================================================
      subroutine fitfun(nn,xdata,ydata,ysm)
c=======================================================================
      implicit none
! argument variables
      integer   nn
      real*8    xdata(nn), ydata(nn), ysm(nn)
! local variables
      integer   i, ii, ill, ismin, j, jj, k, l
      real*8    amat(5,5), bmat(5,5), coef(4), err, rok, rol, x, y
c-----------------------------------------------------------------------
      ismin=1
      do k=1,4
       do l=1,4
        x=0.d0
        do j=0,10
         jj=ismin+j
         rol=1.
         if(l.gt.1)rol=xdata(jj)**(l-1)
         rok=1.
         if(k.gt.1)rok=xdata(jj)**(k-1)
         x=x+rok*rol
        enddO
        amat(k,l)=x
       enddo
       y=0.d0
       do j=0,10
        jj=ismin+j
        rok=1.
        if(k.gt.1)rok=xdata(jj)**(k-1)
        y=y+ydata(jj)*rok
       enddo
       amat(k,5)=y
      enddo
c
      call matslv(5,4,5,amat,bmat,err,ill)
c
      do i=1,4
       coef(i)=amat(i,5)
      enddo
c
      do i=1,11
       ii=ismin+i-1
       ysm(ii)=coef(1)
       do j=1,3
        ysm(ii)=ysm(ii)+coef(j+1)*xdata(ii)**j
       enddo
      enddo
      do i=12,nn
       ii=ismin+i-1
       ysm(ii)=ydata(ii)
      enddo
c
      return
      end subroutine fitfun
c
c=======================================================================
      subroutine recur(xx,aa,bb,cc,dd,nn,mode)
c=======================================================================
      implicit none
! argument
      integer   nn, mode
      real*8    xx(nn), aa(nn), bb(nn), cc(nn), dd(nn)
! local variables
      integer   mm, n, nm
      real*8    ee(1001), ff(1001)
c=======================================================================
      nm=nn-1
      mm=nm-1
      ee(nm)=0.
      ff(nm)=xx(nn)
      do n=mm,1,-1
      ee(n)=-cc(n+1)/(bb(n+1)+aa(n+1)*ee(n+1))
      ff(n)=(dd(n+1)-aa(n+1)*ff(n+1))/(bb(n+1)+aa(n+1)*ee(n+1))
      enddo
c----------------------------------------------------------------------
      xx(1)=0.
c<X'(0)=0>
      if(mode.eq.1)
     >xx(1)=((4.-ee(2))*ff(1)-ff(2))/(3.-(4.-ee(2))*ee(1))
c<X''(0)=0>
      if(mode.eq.2)
     >xx(1)=((2.-ee(2))*ff(1)-ff(2))/(1.+(ee(2)-2.)*ee(1))
c-----
      if(mode.eq.3)
     >xx(1)=(dd(1)-aa(1)*ff(1))/(bb(1)+aa(1)*ee(1))
c----------------------------------------------------------------------
      do n=2,nm
      xx(n)=ee(n-1)*xx(n-1)+ff(n-1)
      enddo
      return
      end subroutine recur
c
c======================================================================
      subroutine csplsf(xfun,yfun,spcoe,sig,sigmax
     >                 ,ilmax,iimax,jjmax,itab)
c======================================================================
c     CUBIC SPLINE WITH LEAST SQUARE FITTING
c======================================================================
      implicit none
! argument
      integer   ilmax, iimax, jjmax, itab(50)
      real*8    xfun(501), yfun(501), spcoe(50,4), sig, sigmax
! local variables
      integer   i, i4, ii, iimax0, iimax1, iimax2, iimax4, iimax5, il
     >        , ill, j, jlmax, k, kkmax, ll, llmax, n,  istop
      real*8    dsig, err, x, xa, xb, xc, xd, xe, xxmat(200,200)
     >        , ymax, yymat(200,200), zfun(501)
c======================================================================
      ymax=0.
      do n=1,ilmax
      if(ymax.lt.dabs(yfun(n))) ymax=dabs(yfun(n))
      enddo
      sig=-1.
      if(ymax.le.0.)return
c-----
      iimax0=iimax-1
      iimax1=jjmax-iimax+1
      if(iimax1.lt.1)iimax1=1
      iimax=iimax-1
      istop=0
      kkmax=0
      do while(istop.eq.0)
      kkmax=kkmax+1
      iimax2=iimax
      iimax=iimax0+kkmax
      if(iimax.le.iimax2)iimax=iimax2+1
      llmax=(ilmax-1)/iimax
      jlmax=llmax*iimax+1
      if(jlmax.lt.ilmax)then
      iimax=iimax+1
      endif
      itab(1)=1
      do i=1,iimax-1
      itab(i+1)=i*llmax+1
      enddo
      itab(iimax+1)=ilmax
c-----
      iimax4=4*iimax
      iimax5=iimax4+1
      do i=1,iimax4
      do j=1,iimax5
      xxmat(i,j)=0.
      enddo
      enddo
c-----
      if(yfun(1).ne.0)then
      xxmat(1,2)=1.
      xxmat(1,1)=1.d-15
      else
      xxmat(1,1)=1.
      endif
      i=1
      do ii=1,iimax
      xa=0.
      xb=0.
      xc=0.
      xd=0.
      xe=0.
      do ll=itab(ii),itab(ii+1)
      x=xfun(ll)
      xa=xa+1.
      xb=xb+x
      xc=xc+x*x
      xd=xd+x*x*x
      xe=xe+yfun(ll)
      enddo
      i=i+1
      xxmat(i,i-1)=xa
      xxmat(i,i )=xb
      xxmat(i,i+1)=xc
      xxmat(i,i+2)=xd
      xxmat(i,iimax5)= xe
c-----
      i=i+1
      xxmat(i,i-2)=1.
      xxmat(i,i-1)=x
      xxmat(i,i )=x*x
      xxmat(i,i+1)=x*x*x
      if(ii.lt.iimax)then
      xxmat(i,i-2+4)=-1
      xxmat(i,i-1+4)=-x
      xxmat(i,i  +4)=-x*x
      xxmat(i,i+1+4)=-x*x*x
      else
      xxmat(i,iimax5)=yfun(ilmax)
      endif
c-----
      if(ii.lt.iimax)then
      i=i+1
      xxmat(i,i-2)=1.
      xxmat(i,i-1)=2.*x
      xxmat(i,i )=3.*x*x
      xxmat(i,i-2+4)=-1
      xxmat(i,i-1+4)=-2.*x
      xxmat(i,i  +4)=-3.*x*x
      endif
c--
      i=i+1
      if(ii.lt.iimax)then
      xxmat(i,i-2)=2.
      xxmat(i,i-1)=6.*x
      xxmat(i,i )=1.d-15
      xxmat(i,i-2+4)=-2.
      xxmat(i,i-1+4)=-6.*x
      else
      xxmat(i,i-1)=2.
      xxmat(i,i )=6.*x
      endif
c--
      enddo
      xxmat(iimax4,2)=1.
      xxmat(iimax4,iimax4)=1.d-15
      call matslv(200,iimax4,iimax5,xxmat,yymat,err,ill)
      do i=1,iimax
      i4=4*i-3
      spcoe(i,  1)=xxmat(i4,iimax5)
      spcoe(i,  2)=xxmat(i4+1,iimax5)
      spcoe(i,  3)=xxmat(i4+2,iimax5)
      spcoe(i,  4)=xxmat(i4+3,iimax5)
      enddo
c-----
      sig=0.
      dsig=0.
      do i=1,iimax
      do k=itab(i),itab(i+1)
      x=xfun(k)
      zfun(k)=spcoe(i,1)+x*(spcoe(i,2)+x*(spcoe(i,3)+x*spcoe(i,4)))
      sig=sig+yfun(k)**2
      dsig=dsig+(zfun(k)-yfun(k))**2
      enddo
      enddo
      i=iimax
      k=ilmax
      x=xfun(k)
      zfun(k)=spcoe(i,1)+x*(spcoe(i,2)+x*(spcoe(i,3)+x*spcoe(i,4)))
      sig=sig+yfun(k)**2
      dsig=dsig+(zfun(k)-yfun(k))**2
      sig=dsig/sig
      if(sig.lt.sigmax)istop=1
      if(kkmax.ge.iimax1)istop=1
      enddo
      do il=1,ilmax
      yfun(il)=zfun(il)
      enddo
      return
      end subroutine csplsf
c
      end module tpxssl_mod
