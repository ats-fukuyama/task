module mod_num_recipe
  implicit none

contains
! Using a combination of Newton-Raphson and bisection, 
! return the root of a function bracketed between x1 and x2.

! Numerical Recipes 3rd ed., p.460
  FUNCTION rtsafe(funcd,x1,x2,xacc,idebug)
    USE nrtype; USE nrutil, ONLY : nrerror
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x1,x2,xacc
    REAL(DP) :: rtsafe
    integer, optional :: idebug
    INTERFACE
       SUBROUTINE funcd(x,fval,fderiv)
         USE nrtype
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x
         REAL(DP), INTENT(OUT) :: fval,fderiv
       END SUBROUTINE funcd
    END INTERFACE
    INTEGER(I4B), PARAMETER :: MAXIT=100
    INTEGER(I4B) :: j
    REAL(DP) :: df,dx,dxold,f,fh,fl,temp,xh,xl
    call funcd(x1,fl,df)
    call funcd(x2,fh,df)
    if(present(idebug)) then
       write(6,*) "input=",x1,x2
       write(6,*) "value=",fl,fh
    end if
    if ((fl > 0.0_dp .and. fh > 0.0_dp) .or. &
         (fl < 0.0_dp .and. fh < 0.0_dp)) &
         call nrerror('root must be bracketed in rtsafe')
    if (fl == 0.0_dp) then
       rtsafe=x1
       RETURN
    else if (fh == 0.0_dp) then
       rtsafe=x2
       RETURN
    else if (fl < 0.0_dp) then
       xl=x1
       xh=x2
    else
       xh=x1
       xl=x2
    end if
    rtsafe=0.5_dp*(x1+x2)
    dxold=abs(x2-x1)
    dx=dxold
    call funcd(rtsafe,f,df)
    do j=1,MAXIT
       if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) > 0.0 .or. &
            abs(2.0_dp*f) > abs(dxold*df) ) then
          dxold=dx
          dx=0.5_dp*(xh-xl)
          rtsafe=xl+dx
          if (xl == rtsafe) RETURN
       else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if (temp == rtsafe) RETURN
       end if
       if (abs(dx) < xacc) RETURN
       call funcd(rtsafe,f,df)
       if (f < 0.0_dp) then
          xl=rtsafe
       else
          xh=rtsafe
       end if
    end do
    call nrerror('rtsafe: exceeded maximum iterations')
  END FUNCTION rtsafe

  ! ---

  ! Cash-Karp type Runge-Kutta method

  SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
    USE nrtype; USE nrutil, ONLY : assert_eq
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: y,dydx
    REAL(DP), INTENT(IN) :: x,h
    REAL(DP), DIMENSION(:), INTENT(OUT) :: yout,yerr
    INTERFACE
       SUBROUTINE derivs(x,y,dydx)
         USE nrtype
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x
         REAL(DP), DIMENSION(:), INTENT(IN) :: y
         REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
       END SUBROUTINE derivs
    END INTERFACE
    INTEGER(I4B) :: ndum
    REAL(DP), DIMENSION(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
    REAL(DP), PARAMETER :: A2=0.2_dp,A3=0.3_dp,A4=0.6_dp,A5=1.0_dp,&
         A6=0.875_dp,B21=0.2_dp,B31=3.0_dp/40.0_dp,B32=9.0_dp/40.0_dp,&
         B41=0.3_dp,B42=-0.9_dp,B43=1.2_dp,B51=-11.0_dp/54.0_dp,&
         B52=2.5_dp,B53=-70.0_dp/27.0_dp,B54=35.0_dp/27.0_dp,&
         B61=1631.0_dp/55296.0_dp,B62=175.0_dp/512.0_dp,&
         B63=575.0_dp/13824.0_dp,B64=44275.0_dp/110592.0_dp,&
         B65=253.0_dp/4096.0_dp,C1=37.0_dp/378.0_dp,&
         C3=250.0_dp/621.0_dp,C4=125.0_dp/594.0_dp,&
         C6=512.0_dp/1771.0_dp,DC1=C1-2825.0_dp/27648.0_dp,&
         DC3=C3-18575.0_dp/48384.0_dp,DC4=C4-13525.0_dp/55296.0_dp,&
         DC5=-277.0_dp/14336.0_dp,DC6=C6-0.25_dp
    ndum=assert_eq(size(y),size(dydx),size(yout),size(yerr),'rkck')
    ytemp=y+B21*h*dydx
    call derivs(x+A2*h,ytemp,ak2)
    ytemp=y+h*(B31*dydx+B32*ak2)
    call derivs(x+A3*h,ytemp,ak3)
    ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
    call derivs(x+A4*h,ytemp,ak4)
    ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
    call derivs(x+A5*h,ytemp,ak5)
    ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
    call derivs(x+A6*h,ytemp,ak6)
    yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
    yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
  END SUBROUTINE rkck

end module mod_num_recipe
