! sak3.f90

MODULE sak3

  PRIVATE
  PUBLIC sak_3

CONTAINS

  SUBROUTINE sak_3

    USE sakcomm
    USE saksub
    USE libgrf
    IMPLICIT NONE
    REAL(dp):: wr,wi    ! real and imag parts of cw
    REAL(dp):: rk2      ! (k vte/omegape)**2 = (k lambda)**2
    REAL(dp):: sg2      ! sigma**2 = (l q_theta/k)**2
    COMPLEX(dp):: cf    ! 
    INTEGER:: nxmax,nymax,nx,ny
    REAL(dp):: delta_nw,eps_nw,wwr,wwi
    INTEGER:: lmax_nw,list_nw,mode_nw,ierr

    wr=1.D0
    wi=0.D0
    rk2=0.1D0
    sg2=0.D0

!        delta_nw: Step size in evaluating derivatives in Newton method
!        eps_nw  : Convergence criterion in Newton method
!        lmax_nw : Maximum iteration count in Newton method
!        list_nw : Listing in Newton method
!        mode_nw : Type of Newton method
    delta_nw = 1.D-6
    eps_nw = 1.D-6
    lmax_nw= 40
    list_nw= 1
    mode_nw= 0

1   CONTINUE
    WRITE(6,'(A/6ES12.4,2I4)') &
         '## INPUT: wr,wi,rk2,sg2,delta_nw,eps_nw,lmax_nw,list_nw?', &
         wr,wi,rk2,sg2,delta_nw,eps_nw,lmax_nw,list_nw
    READ(5,*,ERR=1,END=9000) wr,wi,rk2,sg2,delta_nw,eps_nw,lmax_nw,list_nw

    rk2_local=rk2
    sg2_local=sg2
    CALL newtn0(subeps,wr,wi,wwr,wwi,delta_nw,eps_nw,lmax_nw,list_nw,ierr)

    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE sak_3

!     ****** TWO-DIMENSIONAL NEWTON METHOD (minimum grad f) ****** 

  SUBROUTINE newtn0(sub,xi,yi,xx,yy,delt,eps,ilmax,list,ierr)

    USE sakcomm
    IMPLICIT NONE
    EXTERNAL sub
    REAL(dp),INTENT(IN):: xi,yi    ! initial location
    REAL(dp),INTENT(OUT):: xx,yy  ! final location
    REAL(dp),INTENT(IN):: delt,eps
    INTEGER,INTENT(IN):: ilmax,list
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: il
    REAL(dp):: hx,hy,f00,fm0,fp0,f0m,f0p,fmm,fmp,fpm,fpp
    REAL(dp):: fx,fy,fxx,fyy,fxy,df,det,h11,h12,h21,h22
    REAL(dp):: x,y,dx,dy,v,dv,tt,fxn,fyn,dfn

    ierr=0
    il=0
    x=xi
    y=yi
    IF(ABS(x).GT.1.D0) THEN
       hx=delt*x
    ELSE
       hx=delt
    ENDIF
    IF(ABS(Y).GT.1.D0) THEN
       hy=delt*y
    ELSE
       hy=delt
    ENDIF

    CALL sub(x,   y,   f00)
    IF(LIST.GE.1) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x,y,f00
    IF(f00.LE.2.D-15) GOTO 9001

    CALL sub(x-hx,y,   fm0)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x-hx,y,fm0
    CALL sub(x+hx,y,   fp0)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x+hx,y,fp0
    CALL sub(x,   y-hy,f0m)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x,y-hy,f0m
    CALL sub(x,   y+hy,f0p)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x,y+hy,f0p

    fx=(fp0-fm0)/(2*hx)
    fy=(f0p-f0m)/(2*hy)
    
1   CONTINUE
    CALL sub(x-hx,y-hy,fmm)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x-hx,y-hy,fmm
    CALL sub(x-hx,y+hy,fmp)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x-hx,y+hy,fmp
    CALL sub(x+hx,y-hy,fpm)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x+hx,y-hy,fpm
    CALL sub(x+hx,y+hy,fpp)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x+hx,y+hy,fpp

    fxx=(fp0-2*f00+fm0)/(hx*hx)
    fyy=(f0p-2*f00+f0m)/(hy*hy)
    fxy=(fpp-fpm-fmp+fmm)/(4*hx*hy)
    IF(LIST.GE.2) WRITE(6,'(A,2ES12.4)') 'fx,fy       =',fx,fy
    IF(LIST.GE.2) WRITE(6,'(A,3ES12.4)') 'fxx,fxy,fyy =',fxx,fxy,fyy

    df=SQRT(fx*fx+fy*fy)
    det=fxx*fyy-fxy*fxy
    IF(LIST.GE.2) WRITE(6,'(A,2ES12.4)') 'df,det:0    =',df,det
    IF(df.LE.2.D-15) GO TO 9000

    h11= fyy/det
    h12=-fxy/det
    h21=-fxy/det
    h22= fxx/det
    dx=-(h11*fx+h12*fy)
    dy=-(h21*fx+h22*fy)
    v=SQRT(x*x+y*y+eps)
    dv=SQRT(dx*dx+dy*dy)
    IF(LIST.GE.2) WRITE(6,'(A,3ES12.4)') 'v,dv,dv/v   =',v,dv,dv/v
    IF(dv/v.LE.eps) GO TO 9000

    tt=1.D0
2   CONTINUE
    x=x+tt*dx
    y=y+tt*dy

    CALL sub(x,   y,   f00)
    IF(LIST.GE.1) WRITE(6,'(A,5ES12.4)') 'x,y,f,tt,df =',x,y,f00,tt,df
    IF(f00.LE.2.D-15) GOTO 9001

    CALL sub(x-hx,y,   fm0)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x-hx,y,fm0  =',x-hx,y,fm0
    CALL sub(x+hx,y,   fp0)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x+hx,y,fp0  =',x+hx,y,fp0
    CALL sub(x,   y-hy,f0m)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y-hy,f0m  =',x,y-hy,f0m
    CALL sub(x,   y+hy,f0p)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y+hy,f0p  =',x,y+hy,f0p

    fxn=(fp0-fm0)/(2*hx)
    fyn=(f0p-f0m)/(2*hy)

    dfn=SQRT(fxn*fxn+fyn*fyn)
    IF(dfn.GT.df) THEN
       x=x-tt*dx
       y=y-tt*dy
       tt=0.5d0*tt
       il=il+1
       IF(tt.LE.1.D-3) GOTO 8000
       IF(il.LE.ilmax) GOTO 2
    ELSE
       fx=fxn
       fy=fyn
       df=dfn
       il=il+1
       IF(il.LE.ilmax) GO TO 1
    ENDIF

    ierr=2
    IF(list.GE.1) &
         WRITE(6,*) 'XX NEWTN0: LOOP COUNT EXCEEDS UPPER BOUND.'
    GOTO 9001

8000 ierr=1
    IF(list.GE.1) &
         WRITE(6,*) 'XX NEWTN0: DOES NOT CONVERGE.'
    GOTO 9001

9000 CONTINUE
    CALL sub(x,y,f00)
    IF(LIST.GE.1) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x,y,f00
9001 CONTINUE
    xx=x
    yy=y
    RETURN
  END SUBROUTINE newtn0

!     ****** TWO-DIMENSIONAL NEWTON METHOD (minimum f) ****** 

  SUBROUTINE newtn1(sub,xi,yi,xx,yy,delt,eps,ilmax,list,ierr)

    USE sakcomm
    IMPLICIT NONE
    EXTERNAL sub
    REAL(dp),INTENT(IN):: xi,yi    ! initial location
    REAL(dp),INTENT(OUT):: xx,yy  ! final location
    REAL(dp),INTENT(IN):: delt,eps
    INTEGER,INTENT(IN):: ilmax,list
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: il
    REAL(dp):: hx,hy,f00,fm0,fp0,f0m,f0p,fmm,fmp,fpm,fpp
    REAL(dp):: fx,fy,fxx,fyy,fxy,df,det,h11,h12,h21,h22
    REAL(dp):: x,y,dx,dy,v,dv,tt,fxn,fyn,dfn

    ierr=0
    il=0
    x=xi
    y=yi
    IF(ABS(x).GT.1.D0) THEN
       hx=delt*x
    ELSE
       hx=delt
    ENDIF
    IF(ABS(Y).GT.1.D0) THEN
       hy=delt*y
    ELSE
       hy=delt
    ENDIF
    IF(LIST.GE.1) WRITE(6,'(A,2ES12.4)') 'x,y         =',x,y

    CALL sub(x,   y,   f00)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x,y,f00
    IF(f00.LE.2.D-15) GOTO 9001

    CALL sub(x-hx,y,   fm0)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x-hx,y,fm0
    CALL sub(x+hx,y,   fp0)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x+hx,y,fp0
    CALL sub(x,   y-hy,f0m)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x,y-hy,f0m
    CALL sub(x,   y+hy,f0p)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x,y+hy,f0p

    fx=(fp0-fm0)/(2*hx)
    fy=(f0p-f0m)/(2*hy)
    
1   CONTINUE
    CALL sub(x-hx,y-hy,fmm)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x-hx,y-hy,fmm
    CALL sub(x-hx,y+hy,fmp)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x-hx,y+hy,fmp
    CALL sub(x+hx,y-hy,fpm)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x+hx,y-hy,fpm
    CALL sub(x+hx,y+hy,fpp)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x+hx,y+hy,fpp

    fxx=(fp0-2*f00+fm0)/(hx*hx)
    fyy=(f0p-2*f00+f0m)/(hy*hy)
    fxy=(fpp-fpm-fmp+fmm)/(4*hx*hy)
    IF(LIST.GE.2) WRITE(6,'(A,2ES12.4)') 'fx,fy       =',fx,fy
    IF(LIST.GE.2) WRITE(6,'(A,3ES12.4)') 'fxx,fxy,fyy =',fxx,fxy,fyy

    df=SQRT(fx*fx+fy*fy)
    det=fxx*fyy-fxy*fxy
    IF(LIST.GE.2) WRITE(6,'(A,2ES12.4)') 'df,det:0    =',df,det
    IF(df.LE.2.D-15) GO TO 9000

    h11= fyy/det
    h12=-fxy/det
    h21=-fxy/det
    h22= fxx/det
    dx=-(h11*fx+h12*fy)
    dy=-(h21*fx+h22*fy)
    v=SQRT(x*x+y*y+eps)
    dv=SQRT(dx*dx+dy*dy)
    IF(LIST.GE.2) WRITE(6,'(A,3ES12.4)') 'v,dv,dv/v   =',v,dv,dv/v
    IF(dv/v.LE.eps) GO TO 9000

    tt=1.D0
2   CONTINUE
    x=x+tt*dx
    y=y+tt*dy

    CALL sub(x,   y,   f00)
    IF(LIST.GE.1) WRITE(6,'(A,3ES12.4)') 'x,y,f       =',x,y,f00
    IF(f00.LE.2.D-15) GOTO 9001

    CALL sub(x-hx,y,   fm0)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x-hx,y,fm0  =',x-hx,y,fm0
    CALL sub(x+hx,y,   fp0)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x+hx,y,fp0  =',x+hx,y,fp0
    CALL sub(x,   y-hy,f0m)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y-hy,f0m  =',x,y-hy,f0m
    CALL sub(x,   y+hy,f0p)
    IF(LIST.GE.3) WRITE(6,'(A,3ES12.4)') 'x,y+hy,f0p  =',x,y+hy,f0p

    fxn=(fp0-fm0)/(2*hx)
    fyn=(f0p-f0m)/(2*hy)

    dfn=SQRT(fxn*fxn+fyn*fyn)
    IF(dfn.GT.df) THEN
       x=x-tt*dx
       y=y-tt*dy
       tt=0.5d0*tt
       il=il+1
       IF(tt.LE.1.D-3) GOTO 8000
       IF(il.LE.ilmax) GOTO 2
    ELSE
       fx=fxn
       fy=fyn
       df=dfn
       il=il+1
       IF(il.LE.ilmax) GO TO 1
    ENDIF

    ierr=2
    IF(list.GE.1) &
         WRITE(6,*) 'XX NEWTN0: LOOP COUNT EXCEEDS UPPER BOUND.'
    GOTO 9001

8000 ierr=1
    IF(list.GE.1) &
         WRITE(6,*) 'XX NEWTN0: DOES NOT CONVERGE.'
    GOTO 9001

9000 CONTINUE
    CALL sub(x,y,f00)
9001 CONTINUE
    xx=x
    yy=y
    RETURN
  END SUBROUTINE newtn1

END MODULE sak3
