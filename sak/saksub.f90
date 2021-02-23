! saksub.f90

MODULE saksub
  USE sakcomm
  REAL(dp):: rk_local,sg_local

  PRIVATE
  PUBLIC set_rksg,cfeps,subeps,newtn0,cwaprx

CONTAINS

  SUBROUTINE set_rksg(rk,sg)
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: rk,sg

    rk_local=rk
    sg_local=sg
    RETURN
  END SUBROUTINE set_rksg

  SUBROUTINE subeps(x,y,f)
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: x,y
    REAL(dp),INTENT(OUT):: f
    COMPLEX(dp):: cw,cf

    cw=DCMPLX(x,y)
    cf=cfeps(cw,rk_local,sg_local)
    f=REAL(cf)**2+AIMAG(cf)**2
    RETURN
  END SUBROUTINE subeps
    
  FUNCTION cfeps(cw,rk,sg)
    IMPLICIT NONE
    COMPLEX(dp),INTENT(IN):: cw
    REAL(dp),INTENT(IN):: rk,sg
    COMPLEX(dp):: cfeps,cw2
    REAL(dp):: rk2,sg2

    cw2=cw**2
    rk2=rk**2
    sg2=sg**2
    IF(ABS(rk).LE.1.D-8) THEN
       cfeps=1.D0 &
            -(1.D0+sg2+3.D0*(1.D0+sg2)**2*rk2/cw2)/cw2
    ELSE
       cfeps= 1.D0 &
            -(1.D0+sg2+3.D0*(1.D0+sg2)**2*rk2/cw2)/cw2 &
            +CI*SQRT(0.5D0*Pi)*cw/(rk**3*SQRT(1.D0+sg2)) &
            *EXP(-0.5D0*cw2/(rk2*(1.D0+sg2)))
    END IF
    RETURN
  END FUNCTION cfeps

!     ****** TWO-DIMENSIONAL NEWTON METHOD (minimum grad f) ****** 

  SUBROUTINE newtn0(sub,xi,yi,xx,yy,rd,delt,eps,ilmax,list,ierr)

    USE sakcomm
    IMPLICIT NONE
    EXTERNAL sub
    REAL(dp),INTENT(IN):: xi,yi    ! initial location
    REAL(dp),INTENT(OUT):: xx,yy  ! final location
    REAL(dp),INTENT(OUT):: rd     ! residue
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
    rd=f00
    RETURN
  END SUBROUTINE newtn0

  SUBROUTINE cwaprx(rk,sg,wr1,wi1,wi2)
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: rk,sg
    REAL(dp),INTENT(OUT):: wr1,wi1,wi2
    
    wr1=SQRT((1.D0+sg**2)*(1.D0+3.D0*rk**2))
    IF(ABS(rk).LE.1.D-6) THEN
       wi1=0.D0
       wi2=0.D0
    ELSE
       wi1=-SQRT(0.125D0*Pi)*SQRT(1.D0+sg**2)/rk**3 &
            *EXP(-0.5D0/rk**2)
       wi2=-SQRT(0.125D0*Pi)*SQRT(1.D0+sg**2)*(1.D0+3.D0*rk**2)**2/rk**3 &
            *EXP(-0.5D0/rk**2*(1.D0+3.D0*rk**2))
    END IF
    RETURN
  END SUBROUTINE cwaprx
END MODULE saksub
