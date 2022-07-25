MODULE grdutils

  USE task_kinds,ONLY: dp
  USE grfutils,ONLY: setrgba
  IMPLICIT NONE
  
  PRIVATE
  public gdclip
  PUBLIC ngdlen
  PUBLIC setrgbd
  PUBLIC setrgbda

CONTAINS

  !     ****** AVOID REAL*4 UNDERFLOW ******
  
  FUNCTION gdclip(d)

    REAL(dp),INTENT(IN):: d
    REAL:: gdclip

    IF(ABS(d).LT.1.D-30) then
       gdclip=0.0
    ELSEIF(d.GT. 1.D30) THEN
       gdclip= 1.E30
    ELSEIF(d.LT.-1.D30) THEN
       gdclip=-1.E30
    ELSE
       gdclip=SNGL(d)
    ENDIF
    RETURN
  END FUNCTION gdclip

  !   ***** OPTIMUM NUM LENGTH FOR GVALUE *****

  FUNCTION ngdlen(step)

    REAL(dp),INTENT(IN):: step
    INTEGER:: ngdlen
    REAL(dp):: gxl,gx
    INTEGER:: ngx
    
    gxl=LOG10(step*0.11D0)
    IF(gxl.LT.0.0) THEN
       ngx = -INT(-gxl)
    ELSE
       ngx =  INT( gxl)
    ENDIF
    IF(ngx.LT.-4.OR.ngx.GT.4)THEN
       gx = gxl-ngx
       IF(gx.LT.0.D0) gx=gx+1.D0
       IF(gx.LE.0.15D0) THEN
          ngx=100
       ELSE
          ngx=101
       ENDIF
    ELSEIF(ngx.GE.0) THEN
       ngx=0
    ENDIF
    ngdlen=-ngx
    RETURN
  END FUNCTION ngdlen

  !   ***** setup rgb in dp *****

  SUBROUTINE setrgbd(r,g,b)
    REAL(dp),INTENT(IN):: r,g,b

    CALL setrgb(gdclip(r),gdclip(g),gdclip(b))
    RETURN
  END SUBROUTINE setrgbd

  !   ***** setup rgb in dp *****

  SUBROUTINE setrgbda(rgb)
    REAL(dp),INTENT(IN):: rgb(3)
    REAL:: grgb(3)

    grgb(1)=gdclip(rgb(1))
    grgb(2)=gdclip(rgb(2))
    grgb(3)=gdclip(rgb(3))
    CALL setrgba(grgb)
    RETURN
  END SUBROUTINE setrgbda
END MODULE grdutils
