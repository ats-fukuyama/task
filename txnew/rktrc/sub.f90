module subs
  implicit none

contains

  ! bisection solves f(x) - offset = 0 with respect to x
  ! f(x) is computed by the "sub" subroutine, which is usually set to "spl1df".
  ! x,u and nmax are arguments for sub.
  subroutine bisection(offset,val,sub,x,u,nmax,xmin,xmax,tol)
    integer(4), intent(in)  :: nmax
    real(8),    intent(in)  :: offset, xmin, xmax, tol
    real(8),    intent(in)  :: x(1:nmax), u(1:4,1:nmax)
    real(8),    intent(out) :: val
    integer(4) :: i, ierr, n
    real(8) :: a, b, c, fa, fc

    interface
       subroutine sub(X0,F0,X,U,NXMAX,IERR)
         integer(4),                intent(in) :: NXMAX
         real(8),                   intent(in) :: X0
         real(8),dimension(NXMAX),  intent(in) :: X
         real(8),dimension(4,NXMAX),intent(in) :: U
         real(8),                   intent(out):: F0
         integer(4),                intent(out):: IERR
       end subroutine sub
    end interface

    a = xmin
    b = xmax
    n = log10((b - a) / tol) / log10(2.d0) + 0.5d0
    call sub(a,fa,x,u,nmax,ierr)
    fa = fa - offset
    do i = 1, n
       c = 0.5d0 * (a + b)
       call sub(c,fc,x,u,nmax,ierr)
       fc = fc - offset
       if(fa * fc < 0.d0) then
          b = c
       else
          a = c
       end if
    end do
    val = c

  end subroutine bisection

  ! bisectionf solves f(x) - func = 0 with respect to x
  ! f(x) is computed by the "sub1" subroutine, which is usually set to "spl1df".
  ! func is computed using output from the "sub2" subroutine, which is usually set to "spl1df".
  subroutine bisectionf(func,val,sub1,x,u1,sub2,u2,nmax,xmin,xmax,tol,r1,z1,r2,z2)
    integer(4), intent(in)  :: nmax
    real(8),    intent(in)  :: xmin, xmax, tol, r1, z1, r2, z2
    real(8),    intent(in)  :: x(1:nmax), u1(1:4,1:nmax), u2(1:4,1:nmax)
    real(8),    intent(out) :: val
    integer(4) :: i, ierr, n
    real(8) :: a, b, c, fa, fc, r
    real(8), external :: func

    interface
       subroutine sub1(X0,F0,X,U,NXMAX,IERR)
         integer(4),                intent(in) :: NXMAX
         real(8),                   intent(in) :: X0
         real(8),dimension(NXMAX),  intent(in) :: X
         real(8),dimension(4,NXMAX),intent(in) :: U
         real(8),                   intent(out):: F0
         integer(4),                intent(out):: IERR
       end subroutine sub1
    end interface

    interface
       subroutine sub2(X0,F0,X,U,NXMAX,IERR)
         integer(4),                intent(in) :: NXMAX
         real(8),                   intent(in) :: X0
         real(8),dimension(NXMAX),  intent(in) :: X
         real(8),dimension(4,NXMAX),intent(in) :: U
         real(8),                   intent(out):: F0
         integer(4),                intent(out):: IERR
       end subroutine sub2
    end interface

    a = xmin
    b = xmax
    n = log10((b - a) / tol) / log10(2.d0) + 0.5d0
    call sub2(a,r ,x,u2,nmax,ierr)
    call sub1(a,fa,x,u1,nmax,ierr)
    fa = fa - func(r,r1,z1,r2,z2)
    do i = 1, n
       c = 0.5d0 * (a + b)
       call sub2(c,r ,x,u2,nmax,ierr)
       call sub1(c,fc,x,u1,nmax,ierr)
       fc = fc - func(r,r1,z1,r2,z2)
       if(fa * fc < 0.d0) then
          b = c
       else
          a = c
       end if
    end do
    val = c

  end subroutine bisectionf

  ! bisection solves f'(x) - offset = 0 with respect to x
  ! f'(x) is computed by the "sub" subroutine, which is usually set to "spl1dd".
  ! x,u and nmax are arguments for sub.
  subroutine bisectiond(offset,val,fval,sub,x,u,nmax,xmin,xmax,tol)
    integer(4), intent(in)  :: nmax
    real(8),    intent(in)  :: offset, xmin, xmax, tol
    real(8),    intent(in)  :: x(1:nmax), u(1:4,1:nmax)
    real(8),    intent(out) :: val, fval
    integer(4) :: i, ierr, n
    real(8) :: a, b, c, fa, fda, fc, fdc

    interface
       subroutine sub(X0,F0,FX0,X,U,NXMAX,IERR)
         integer(4),                intent(in) :: NXMAX
         real(8),                   intent(in) :: X0
         real(8),dimension(NXMAX),  intent(in) :: X
         real(8),dimension(4,NXMAX),intent(in) :: U
         real(8),                   intent(out):: F0, FX0
         integer(4),                intent(out):: IERR
       end subroutine sub
    end interface

    a = xmin
    b = xmax
    n = log10((b - a) / tol) / log10(2.d0) + 0.5d0
    call sub(a,fa,fda,x,u,nmax,ierr)
    fda = fda - offset
    do i = 1, n
       c = 0.5d0 * (a + b)
       call sub(c,fc,fdc,x,u,nmax,ierr)
       fdc = fdc - offset
       if(fda * fdc < 0.d0) then
          b = c
       else
          a = c
       end if
    end do
    val  = c
    fval = fc

  end subroutine bisectiond

  !     ***** 1D Newton's method *****
  ! newton1d solves f(x) = 0 or f(x) = vtarget with respect to x
  ! f(x) is computed by the "sub" subroutine, which has to be "spl1dd" in usual.
  ! x,u and nmax are arguments for sub.

  subroutine newton1d(x_in,x,sub,adim,uspl,ndim,tol,vtarget)
    real(8), intent(in)  :: x_in, tol, adim(1:ndim), uspl(1:4,1:ndim)
    real(8), intent(in), optional :: vtarget
    real(8), intent(out) :: x
    integer(4), intent(in)  :: ndim

    real(8) :: fx, fdx, dx
    integer(4) :: n, ierr

    interface
       subroutine sub(X0,F0,FX0,X,U,NXMAX,IERR)
         integer(4),                intent(in) :: NXMAX
         real(8),                   intent(in) :: X0
         real(8),dimension(NXMAX),  intent(in) :: X
         real(8),dimension(4,NXMAX),intent(in) :: U
         real(8),                   intent(out):: F0, FX0
         integer(4),                intent(out):: IERR
       end subroutine sub
    end interface

    n = 0
    x = x_in

    do
       call sub(x,fx,fdx,adim,uspl,ndim,ierr)
       if( ierr /= 0 ) then
          write(6,*) 'newton1d: spline error.'
          stop
       endif
       dx = ( vtarget - fx ) / fdx
       !     write(6,'(I3,4ES15.7)') n,x,fx,fdx,dx
       n = n + 1
       x = x + dx
       if( abs( dx / x ) < tol ) exit
       if( n >= 100 ) then
          write(6,*) 'newton1d: too much iteration. Aborting!'
          stop
       endif
    enddo

    return
  end subroutine newton1d

  !     ***** 1D Newton's method *****
  ! newton1d solves f(x) = 0 or f(x) = vtarget with respect to x
  ! f(x) is computed by the "sub" subroutine, which has to be "spl1dd" in usual.
  ! x,u and nmax are arguments for sub.

  subroutine newton1df(x_in,x,sub,adim,uspl,sub2,uspl2,ndim,tol,func,r1,z1,r2,z2)
    real(8), intent(in)  :: x_in, tol, adim(1:ndim), uspl(1:4,1:ndim), uspl2(1:4,1:ndim), r1, z1, r2, z2
    real(8), intent(out) :: x
    integer(4), intent(in)  :: ndim
    real(8), external :: func

    real(8) :: fx, fdx, dx, ffunc, fdfunc, rl, drl
    integer(4) :: n, ierr

    interface
       subroutine sub(X0,F0,FX0,X,U,NXMAX,IERR)
         integer(4),                intent(in) :: NXMAX
         real(8),                   intent(in) :: X0
         real(8),dimension(NXMAX),  intent(in) :: X
         real(8),dimension(4,NXMAX),intent(in) :: U
         real(8),                   intent(out):: F0, FX0
         integer(4),                intent(out):: IERR
       end subroutine sub
    end interface
    interface
       subroutine sub2(X0,F0,FX0,X,U,NXMAX,IERR)
         integer(4),                intent(in) :: NXMAX
         real(8),                   intent(in) :: X0
         real(8),dimension(NXMAX),  intent(in) :: X
         real(8),dimension(4,NXMAX),intent(in) :: U
         real(8),                   intent(out):: F0, FX0
         integer(4),                intent(out):: IERR
       end subroutine sub2
    end interface

    n = 0
    x = x_in

    fdfunc = (z2-z1)/(r2-r1)
    do
       call sub2(x,rl,drl,adim,uspl2,ndim,ierr) ! get R(l)
       ffunc = func(rl,r1,z1,r2,z2)
       call sub(x,fx,fdx,adim,uspl,ndim,ierr)
       if( ierr /= 0 ) then
          write(6,*) 'newton1df: spline error.'
          stop
       endif
!       write(6,'(I3,5ES15.7)') n,x,rl,fx,ffunc,fx-ffunc
       fx  = fx - ffunc
       fdx = fdx - fdfunc * drl

       dx = - fx / fdx
       !     write(6,'(I3,4ES15.7)') n,x,fx,fdx,dx
       n = n + 1
       x = x + dx
       if( abs( dx / x ) < tol ) exit
       if( n >= 100 ) then
          write(6,*) 'newton1d: too much iteration. Aborting!'
          stop
       endif
    enddo

    return
  end subroutine newton1df

  !     ***** 1D Newton's method *****
  ! newton1dd solves df(x)/dx = 0 with respect to x to obtain a maximum or minimum value
  ! f(x) is computed by the "sub" subroutine, which has to be "spl1ddd" in usual.
  ! x,u and nmax are arguments for sub.

  subroutine newton1dd(x_in,x,fx,sub,adim,uspl,ndim,tol)
    real(8), intent(in)  :: x_in, tol, adim(1:ndim), uspl(1:4,1:ndim)
    real(8), intent(out) :: x, fx
    integer(4), intent(in)  :: ndim

    real(8) :: fdx, fddx, dx
    integer(4) :: n, ierr

    interface
       subroutine sub(X0,F0,FX0,FXX0,X,U,NXMAX,IERR)
         integer(4),                intent(in) :: NXMAX
         real(8),                   intent(in) :: X0
         real(8),dimension(NXMAX),  intent(in) :: X
         real(8),dimension(4,NXMAX),intent(in) :: U
         real(8),                   intent(out):: F0, FX0, FXX0
         integer(4),                intent(out):: IERR
       end subroutine sub
    end interface

    n = 0
    x = x_in

    do
       call sub(x,fx,fdx,fddx,adim,uspl,ndim,ierr)
       if( ierr /= 0 ) then
          write(6,*) 'newton1dd: spline error.'
          stop
       endif
       dx = - fdx / fddx
       n = n + 1
       x = x + dx
!       write(6,'(I3,5ES15.7)') n,x,fx,fdx,fddx,abs( dx / x )
       if( abs( dx / x ) < tol ) then
          call sub(x,fx,fdx,fddx,adim,uspl,ndim,ierr)
          exit
       endif
       if( n >= 100 ) then
          write(6,*) 'newton1dd: too much iteration. Aborting!'
          stop
       endif
    enddo

    return
  end subroutine newton1dd

!     ****** TWO-DIMENSIONAL NEWTON METHOD ******

  subroutine NEWTN(SUB,X,Y,XX,YY,DELT,EPS,ILMAX,LIST,IER)
    real(8)   , intent(IN)   :: DELT, EPS
    real(8)   , intent(INOUT):: X, Y
    real(8)   , intent(OUT)  ::XX,YY
    integer(4), intent(IN)   ::ILMAX,LIST
    integer(4), intent(OUT)  ::IER
    real(8)     :: DET, DF, DFN, DX, DY, FXX, FYY, FXY, FXY1, FXY2, G0, GN, GX, GY, &
         &         H11, H12, H21, H22, HX, HY, S0, SN, SX, SY, TT
    integer(4)  :: ITER
!    external SUB

    interface
       subroutine sub(r,z,dpsidr,dpsidz,zpsi)
         real(8), intent(in)  :: r, z
         real(8), intent(out) :: dpsidr, dpsidz
         real(8), intent(out), optional :: zpsi
       end subroutine sub
    end interface

    IER=0
    ITER=0
    if(abs(X) > 1.D0) then
       HX=DELT*X
    else
       HX=DELT
    endif
    if(abs(Y) > 1.D0) then
       HY=DELT*Y
    else
       HY=DELT
    endif
    call SUB(X,   Y,   S0,G0)
    if(LIST > 0) write(6,600) X,Y,S0,G0
    if(LIST > 0) write(6,*) 'HX=',HX,'HY=',HY
1   call SUB(X+HX,Y,   SX,GX)
    call SUB(X,   Y+HY,SY,GY)
    FXX =(SX-S0)/HX
    FXY1=(GX-G0)/HX
    FYY =(GY-G0)/HY
    FXY2=(SY-S0)/HY
    FXY=0.5D0*(FXY1+FXY2)
    if(LIST > 1) write(6,601) SX,SY,GX,GY
    if(LIST > 1) write(6,601) FXX,FYY,FXY1,FXY2
    DF=sqrt(S0*S0+G0*G0)
    DET=FXX*FYY-FXY*FXY
    H11= FYY/DET
    H12=-FXY/DET
    H21=-FXY/DET
    H22= FXX/DET

    DX=-(H11*S0+H12*G0)
    DY=-(H21*S0+H22*G0)
    TT=1.0D0
2   X=X+TT*DX
    Y=Y+TT*DY
    call SUB(X,   Y,   SN,GN)
    if(LIST > 0) write(6,602) DX,DY
    if(LIST > 0) write(6,600) X,Y,SN,GN
    DFN=sqrt(SN*SN+GN*GN)
    if(DFN > DF) then
       X=X-TT*DX
       Y=Y-TT*DY
       if(DF < EPS) goto 900
       TT=0.5D0*TT
       ITER=ITER+1
       if(TT < 1.D-3) goto 800
       if(ITER < ILMAX) goto 2
    else
       S0=SN
       G0=GN
       DF=DFN
       ITER=ITER+1
       if(DF < EPS) goto 900
       if(ITER < ILMAX) goto 1
    endif

    IER=2
    if(LIST > 0) &
         &write(6,*) 'XX NEWTN: LOOP COUNT EXCEEDS UPPER BOUND.'
    goto 900

800 IER=1
    if(LIST > 0) &
         &write(6,*) 'XX NEWTN: DOES NOT CONVERGE.'
    goto 900

900 XX=X
    YY=Y
    return
600 format(" ",6X,'X,Y,FX,FY = ',4ES15.7)
601 format(" ",6X,'FXX,YY,XY = ',4ES15.7)
602 format(" ",6X,'DX,DY     = ',2ES15.7)
  end subroutine NEWTN

!     ***** One-variable Brent Method *****

  real(8) function fbrent(f,ax,bx,tol)
    real(8), intent(in):: ax,bx,tol
    real(8), external  :: f

!      a zero of the function  f(x)  is computed in the interval ax,bx .

!  input..

!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result (.ge.0.)

!  output..

!  zeroin abscissa approximating a zero of  f  in the interval ax,bx

!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  this is checked, and an error message is printed if this is not
!  satisfied.   zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
!  the  relative machine precision defined as the smallest representable
!  number such that  1.+macheps .gt. 1.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice-hall, inc. (1973).

    real(8) :: a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s
    integer(4) :: it
    real(8),    parameter :: eps=1.D-15
    integer(4), parameter :: itmax=100

    tol1 = eps+1.0d0

    a=ax
    b=bx
    fa=f(a)
    fb=f(b)
!     check that f(ax) and f(bx) have different signs
    if (fa .ne. 0.0d0 .and. fb .ne. 0.0d0) then
       if (fa * (fb/abs(fb)) .gt. 0.0d0) then
          write(6,*) 'XX FBRENT: f(ax) and f(bx) have same sign.'
          write(6,'(A,4ES12.4)') '  ax,bx,f(ax),f(bx)=',a,b,fa,fb
          return
       endif
    endif
    c=b
    fc=fb

    DO IT=1,ITMAX
       if ((fb*(fc/abs(fc))).gt.0.0d0) THEN
          c=a
          fc=fa
          d=b-a
          e=d
       endif

       if (abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       endif
       tol1=2.0d0*eps*abs(b)+0.5d0*tol
       xm = 0.5d0*(c-b)

       if ((abs(xm).le.tol1).or.(fb.eq.0.0d0)) then
          fbrent=b
          return
       endif

       ! see if a bisection is forced

       if ((abs(e).ge.tol1).and.(abs(fa).gt.abs(fb))) then
          s=fb/fa
          if (a.eq.c) then

             ! linear interpolation

             p=2.0d0*xm*s
             q=1.0d0-s
          else

             ! inverse quadratic interpolation

             q=fa/fc
             r=fb/fc
             p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
             q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
          endif
          if (p.gt.0.0d0) then
             q=-q
          else
             p=-p
          endif
          s=e
          e=d
          if (((2.0d0*p).lt.(3.0d0*xm*q-abs(tol1*q))).and. &
               &           (p.lt.abs(0.5d0*s*q))) then
             d=p/q
          else
             d=xm
             e=d
          endif
       else
          d=xm
          e=d
       endif
       a=b
       fa=fb
       if (abs(d).gt.tol1) then
          b=b+d
       else
          if (xm.gt.0.0d0) then
             b=b+tol1
          else
             b=b-tol1
          endif
       endif
       fb=f(b)
    ENDDO

    WRITE(6,*) 'XX FBRENT EXCEEDING MAXIMUM ITERATIONS'
    fbrent=b
    return
  end function fbrent

end module subs
