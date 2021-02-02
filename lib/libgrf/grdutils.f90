MODULE grdutils

  USE task_kinds,ONLY: dp
  USE grfutils,ONLY: setrgba
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC gsclip,ngslen,setrgbd,setrgbda

CONTAINS

  !     ****** AVOID REAL*4 UNDERFLOW ******
  
  FUNCTION gsclip(d)

    REAL(dp),INTENT(IN):: d
    REAL:: gsclip

    IF(ABS(d).LT.1.D-30) then
       gsclip=0.0
    ELSEIF(d.GT. 1.D30) THEN
       gsclip= 1.E30
    ELSEIF(d.LT.-1.D30) THEN
       gsclip=-1.E30
    ELSE
       gsclip=SNGL(d)
    ENDIF
    RETURN
  END FUNCTION gsclip

  !   ***** OPTIMUM NUM LENGTH FOR GVALUE *****

  FUNCTION ngslen(step)

    REAL(dp),INTENT(IN):: step
    INTEGER:: ngslen
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
    ngslen=-ngx
    RETURN
  END FUNCTION ngslen

  !   ***** setup rgb in dp *****

  SUBROUTINE setrgbd(r,g,b)
    REAL(dp),INTENT(IN):: r,g,b

    CALL setrgb(gsclip(r),gsclip(g),gsclip(b))
    RETURN
  END SUBROUTINE setrgbd

  !   ***** setup rgb in dp *****

  SUBROUTINE setrgbda(rgb)
    REAL(dp),INTENT(IN):: rgb(3)
    REAL:: grgb(3)

    grgb(1)=gsclip(rgb(1))
    grgb(2)=gsclip(rgb(2))
    grgb(3)=gsclip(rgb(3))
    CALL setrgba(grgb)
    RETURN
  END SUBROUTINE setrgbda

END MODULE grdutils
