MODULE grfutils

  USE task_kinds,ONLY: dp
  PRIVATE
  PUBLIC grf_title,grf_frame1d,grf_frame2d,grf_frame3d, &
         grf_info,grfut1,grfut2,grfut3,grfut4,setrgba, &
         gsclip,ngslen

CONTAINS

! ***** DRAW GRAPH TITLES *****

  SUBROUTINE GRF_TITLE(A)

    USE grftype, ONLY: grf_attr_type
    IMPLICIT NONE

    TYPE(grf_attr_type),INTENT(IN):: A
    REAL:: SIZE
    INTEGER:: NG,NL

    IF(A%TITLE_LEN.GT.0) THEN
       CALL SETCHS(A%TITLE_SIZE,0.0)
       CALL SETFNT(A%TITLE_FONT)
       CALL SETRGBA(A%TITLE_RGB(1:3))
       CALL INQTSZ(A%TITLE,A%TITLE_LEN,SIZE)
       SELECT CASE(A%TITLE_POS)
       CASE(0)
          CALL MOVE(A%GPXMIN,A%GPYMAX+A%TITLE_SEP)
       CASE(1)
          CALL MOVE(0.5*(A%GPXMIN+A%GPXMAX)-0.5*SIZE,A%GPYMAX+A%TITLE_SEP)
       CASE(2)
          CALL MOVE(A%GPXMAX-SIZE,A%GPYMAX+A%TITLE_SEP)
       CASE DEFAULT
          CALL MOVE(A%GPXMIN,A%GPYMAX+A%TITLE_SEP)
       END SELECT
       CALL TEXT(A%TITLE,A%TITLE_LEN)
    ENDIF

    IF(A%XTITLE_LEN.GT.0) THEN
       CALL SETCHS(A%XTITLE_SIZE,0.0)
       CALL SETFNT(A%XTITLE_FONT)
       CALL SETRGBA(A%XTITLE_RGB(1:3))
       CALL INQTSZ(A%XTITLE,A%XTITLE_LEN,SIZE)
       CALL MOVE(0.5*(A%GPXMIN+A%GPXMAX)-0.5*SIZE,A%GPYMIN-A%XTITLE_SEP)
       CALL TEXT(A%XTITLE,A%XTITLE_LEN)
    ENDIF

    IF(A%YTITLE_LEN.GT.0) THEN
       CALL SETCHS(A%YTITLE_SIZE,90.0)
       CALL SETFNT(A%YTITLE_FONT)
       CALL SETRGBA(A%YTITLE_RGB(1:3))
       CALL INQTSZ(A%YTITLE,A%YTITLE_LEN,SIZE)
       CALL MOVE(A%GPXMIN-A%YTITLE_SEP,0.5*(A%GPYMIN+A%GPYMAX)-0.5*SIZE)
       CALL TEXT(A%YTITLE,A%YTITLE_LEN)
    ENDIF

    RETURN
  END SUBROUTINE GRF_TITLE

! ***** DRAW GRAPH FRAME AND SCALES *****

  SUBROUTINE GRF_FRAME1D(A)

    USE grftype, ONLY: grf_attr_type
    IMPLICIT NONE

    TYPE(grf_attr_type),INTENT(IN):: A
    REAL:: SIZE
    INTEGER:: NG,NL

    CALL GDEFIN(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX, &
                A%XMIN,A%XMAX,A%FMIN,A%FMAX)

    CALL SETRGBA(A%FRAME_RGB(1:3))
    CALL SETLNW(A%FRAME_THICKNESS)
    CALL GFRAME

    CALL SETLNW(A%SCALE_THICKNESS)
    CALL SETCHS(A%VALUE_SIZE,0.0)
    CALL SETFNT(A%VALUE_FONT)
    CALL SETRGBA(A%SCALE_RGB(1:3))

!    WRITE(6,'(A,4ES12.4)') 'GRF_FRAME1D:',A%XORG,A%XMIN,A%XMAX,A%XSCALE_STEP
    
    IF(A%MODE_LS.EQ.0.OR.A%MODE_LS.EQ.2) THEN
       IF(A%XSCALE_ZERO.NE.0) &
            CALL GSCALE(0.0,2*(A%XMAX-A%XMIN),0.0,0.0, &
                        FLOAT(A%XSCALE_ZERO_PAT),0)
       CALL GSCALE(A%XORG,A%XSCALE_STEP,0.0,0.0,A%XSCALE_SIZE,A%XSCALE_TYPE)
       IF(A%XGRID_STEP.NE.0.0) &
            CALL GSCALE(A%XORG,A%XGRID_STEP,0.0,0.0,A%XSCALE_SIZE,0)
       CALL SETRGBA(A%VALUE_RGB(1:3))
       CALL GVALUE(A%XORG,A%XVALUE_STEP*A%XSCALE_STEP,0.0,0.0,A%XVALUE_TYPE)
    ELSE
       CALL GSCALL(A%XORG,A%XSCALE_LTYPE,0.0,0,A%XSCALE_SIZE,A%XSCALE_TYPE)
       IF(A%XGRID_LTYPE.NE.0) &
            CALL GSCALL(A%XORG,A%XGRID_LTYPE,0.0,0,A%XSCALE_SIZE,0)
       CALL SETRGBA(A%VALUE_RGB(1:3))
       CALL GVALUL(A%XORG,A%XVALUE_LTYPE,0.0,0,A%XVALUE_TYPE)
    ENDIF

    CALL SETRGBA(A%SCALE_RGB(1:3))
    IF(A%MODE_LS.EQ.0.OR.A%MODE_LS.EQ.1) THEN
       IF(A%FSCALE_ZERO.NE.0) &
            CALL GSCALE(0.0,0.0,0.0,2*(A%FMAX-A%FMIN), &
                        FLOAT(A%FSCALE_ZERO_PAT),0)
       CALL GSCALE(0.0,0.0,A%FORG,A%FSCALE_STEP,A%FSCALE_SIZE,A%FSCALE_TYPE)
       IF(A%FGRID_STEP.NE.0.0) &
            CALL GSCALE(0.0,0.0,A%FORG,A%FGRID_STEP,A%FSCALE_SIZE,0)
       CALL SETRGBA(A%VALUE_RGB(1:3))
       CALL GVALUE(0.0,0.0,A%FORG,A%FVALUE_STEP*A%FSCALE_STEP,A%FVALUE_TYPE)
    ELSE
       CALL GSCALL(0.0,0,A%FORG,A%FSCALE_LTYPE,A%FSCALE_SIZE,A%FSCALE_TYPE)
       IF(A%FGRID_LTYPE.NE.0) &
            CALL GSCALL(0.0,0,A%FORG,A%FGRID_LTYPE,A%FSCALE_SIZE,0)
       CALL SETRGBA(A%VALUE_RGB(1:3))
       CALL GVALUL(0.0,0,A%FORG,A%FVALUE_LTYPE,A%FVALUE_TYPE)
    ENDIF

  END SUBROUTINE GRF_FRAME1D

  SUBROUTINE GRF_FRAME2D(A)

    USE grftype, ONLY: grf_attr_type
    IMPLICIT NONE

    TYPE(grf_attr_type),INTENT(IN):: A
    REAL:: SIZE
    INTEGER:: NG,NL

    CALL GDEFIN(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX, &
                A%XMIN,A%XMAX,A%YMIN,A%YMAX)

    CALL SETRGBA(A%FRAME_RGB(1:3))
    CALL SETLNW(A%FRAME_THICKNESS)
    CALL GFRAME

    CALL SETLNW(A%SCALE_THICKNESS)
    CALL SETCHS(A%VALUE_SIZE,0.0)
    CALL SETFNT(A%VALUE_FONT)
    CALL SETRGBA(A%SCALE_RGB(1:3))

    IF(A%MODE_LS.EQ.0.OR.A%MODE_LS.EQ.2) THEN
       IF(A%XSCALE_ZERO.NE.0) &
            CALL GSCALE(0.0,2*(A%XMAX-A%XMIN),0.0,0.0, &
                        FLOAT(A%XSCALE_ZERO_PAT),0)
       CALL GSCALE(A%XORG,A%XSCALE_STEP,0.0,0.0,A%XSCALE_SIZE,A%XSCALE_TYPE)
       IF(A%XGRID_STEP.NE.0.0) &
            CALL GSCALE(A%XORG,A%XGRID_STEP,0.0,0.0,A%XSCALE_SIZE,0)
       CALL SETRGBA(A%VALUE_RGB(1:3))
       CALL GVALUE(A%XORG,A%XVALUE_STEP*A%XSCALE_STEP,0.0,0.0,A%XVALUE_TYPE)
    ELSE
       CALL GSCALL(A%XORG,A%XSCALE_LTYPE,0.0,0,A%XSCALE_SIZE,A%XSCALE_TYPE)
       IF(A%XGRID_LTYPE.NE.0) &
            CALL GSCALL(A%XORG,A%XGRID_LTYPE,0.0,0,A%XSCALE_SIZE,0)
       CALL SETRGBA(A%VALUE_RGB(1:3))
       CALL GVALUL(A%XORG,A%XVALUE_LTYPE,0.0,0,A%XVALUE_TYPE)
    ENDIF

    CALL SETRGBA(A%SCALE_RGB(1:3))
    IF(A%MODE_LS.EQ.0.OR.A%MODE_LS.EQ.1) THEN
       IF(A%YSCALE_ZERO.NE.0) &
            CALL GSCALE(0.0,0.0,0.0,2*(A%YMAX-A%YMIN), &
                        FLOAT(A%YSCALE_ZERO_PAT),0)
       CALL GSCALE(0.0,0.0,A%YORG,A%YSCALE_STEP,A%YSCALE_SIZE,A%YSCALE_TYPE)
       IF(A%XGRID_STEP.NE.0.0) &
            CALL GSCALE(0.0,0.0,A%XORG,A%XGRID_STEP,A%XSCALE_SIZE,0)
       CALL SETRGBA(A%VALUE_RGB(1:3))
       CALL GVALUE(0.0,0.0,A%YORG,A%YVALUE_STEP*A%YSCALE_STEP,A%YVALUE_TYPE)
    ELSE
       CALL GSCALL(0.0,0,A%YORG,A%YSCALE_LTYPE,A%YSCALE_SIZE,A%YSCALE_TYPE)
       IF(A%YGRID_LTYPE.NE.0) &
            CALL GSCALL(0.0,0,A%YORG,A%YGRID_LTYPE,A%YSCALE_SIZE,0)
       CALL SETRGBA(A%VALUE_RGB(1:3))
       CALL GVALUL(0.0,0,A%YORG,A%YVALUE_LTYPE,A%YVALUE_TYPE)
    ENDIF

  END SUBROUTINE GRF_FRAME2D

  SUBROUTINE GRF_FRAME3D(A)

    USE grftype, ONLY: grf_attr_type
    IMPLICIT NONE

    TYPE(grf_attr_type),INTENT(IN):: A
    REAL:: SIZE
    INTEGER:: NG,NL

    CALL SETRGBA(A%SCALE_RGB(1:3))
    CALL SETLNW(A%SCALE_THICKNESS)

    CALL GAXIS3D(0)

    CALL GSCALE3DX(A%XORG,A%XSCALE_STEP,A%XSCALE_SIZE,A%XSCALE_TYPE)
    CALL GSCALE3DY(A%YORG,A%YSCALE_STEP,A%YSCALE_SIZE,A%YSCALE_TYPE)
    CALL GSCALE3DZ(A%FORG,A%FSCALE_STEP,A%FSCALE_SIZE,A%FSCALE_TYPE)

    CALL SETFNT(A%VALUE_FONT)
    CALL SETCHS(A%VALUE_SIZE,0.0)
    CALL SETRGBA(A%VALUE_RGB(1:3))

    CALL GVALUE3DX(A%XORG,A%XVALUE_STEP*A%XSCALE_STEP,A%XVALUE_POS, &
                   A%XVALUE_TYPE)
    CALL GVALUE3DY(A%YORG,A%YVALUE_STEP*A%YSCALE_STEP,A%YVALUE_POS, &
                   A%YVALUE_TYPE)
    CALL GVALUE3DZ(A%FORG,A%FVALUE_STEP*A%FSCALE_STEP,A%FVALUE_POS, &
                   A%FVALUE_TYPE)
    RETURN
  END SUBROUTINE GRF_FRAME3D

! ***** DRAW GRAPH INFO *****

  SUBROUTINE GRF_INFO(A)

    USE grftype, ONLY: grf_attr_type
    IMPLICIT NONE

    TYPE(grf_attr_type),INTENT(IN):: A
    REAL:: STEP
    CHARACTER(LEN=64) LINE

    SELECT CASE(A%MODE_2D)
    CASE(1)
       STEP=A%LINE_VALUE(2)-A%LINE_VALUE(1)
    CASE(2)
       STEP=A%PAINT_VALUE(2)-A%PAINT_VALUE(1)
    CASE DEFAULT
       STEP=0.0
    END SELECT
       
    CALL SETCHS(0.5*A%VALUE_SIZE,0.0)
    CALL SETRGB(0.0,0.0,0.0)
!    WRITE(LINE,'(A,1PE12.4,A,1PE12.4,A,1PE12.4)') &
!         'MIN  =',A%FMIN,'     MAX  =',A%FMAX,'     STEP =', STEP
!    CALL MOVE(A%GPXMAX-32*A%VALUE_SIZE/1.5, &
!              A%GPYMAX+A%TITLE_SEP)
!    CALL TEXT(LINE,64)
    WRITE(LINE,'(A,1PE12.4)') 'MAX  =',A%FMAX
    CALL MOVE(A%GPXMAX-6*A%VALUE_SIZE,A%GPYMAX+1.5*A%TITLE_SEP+A%VALUE_SIZE)
    CALL TEXT(LINE,18)
    WRITE(LINE,'(A,1PE12.4)') 'MIN  =',A%FMIN
    CALL MOVE(A%GPXMAX-6*A%VALUE_SIZE,A%GPYMAX+A%TITLE_SEP+0.5*A%VALUE_SIZE)
    CALL TEXT(LINE,18)
    WRITE(LINE,'(A,1PE12.4)') 'STEP =',STEP
    CALL MOVE(A%GPXMAX-6*A%VALUE_SIZE,A%GPYMAX+0.5*A%TITLE_SEP)
    CALL TEXT(LINE,18)
    CALL SETRGB(0.0,0.0,0.0)
    RETURN
  END SUBROUTINE GRF_INFO


!     ***** EVALUATE GXMIN,GXMAX,GYMIN,GYMAX *****

  SUBROUTINE GRFUT1(GX,NXMAX,GXMIN,GXMAX,NXMIN,NXSTEP)

    IMPLICIT NONE

    INTEGER,INTENT(IN):: NXMAX
    REAL,INTENT(IN):: GX(NXMAX)
    REAL,INTENT(OUT):: GXMIN,GXMAX
    INTEGER,INTENT(IN),OPTIONAL:: NXMIN,NXSTEP
    INTEGER:: NXMINL,NXSTEPL

    IF(PRESENT(NXMIN)) THEN
       NXMINL=NXMIN
    ELSE
       NXMINL=1
    ENDIF
    IF(PRESENT(NXSTEP)) THEN
       NXSTEPL=NXSTEP
    ELSE
       NXSTEPL=1
    ENDIF
    CALL GMNMX1(GX,NXMINL,NXMAX,NXSTEPL,GXMIN,GXMAX)
    IF(ABS(GXMAX-GXMIN).LT.1.E-6) THEN
       GXMIN=GXMIN-0.999E-6
       GXMAX=GXMAX+1.000E-6
    ENDIF

    RETURN
  END SUBROUTINE GRFUT1

!     ***** EVALUATE GXMIN,GXMAX,GYMIN,GYMAX *****

  SUBROUTINE GRFUT2(GY,NXM,NXMAX,NYMAX,GYMIN,GYMAX, &
                       NXMIN,NXSTEP,NYMIN,NYSTEP)

    IMPLICIT NONE

    INTEGER,INTENT(IN):: NXM,NXMAX,NYMAX
    REAL,INTENT(IN):: GY(NXM,NYMAX)
    REAL,INTENT(OUT):: GYMIN,GYMAX
    INTEGER,INTENT(IN),OPTIONAL:: NXMIN,NXSTEP,NYMIN,NYSTEP
    INTEGER:: NXMINL,NXSTEPL,NYMINL,NYSTEPL
    REAL:: GYNORM

    IF(PRESENT(NXMIN)) THEN
       NXMINL=NXMIN
    ELSE
       NXMINL=1
    ENDIF
    IF(PRESENT(NXSTEP)) THEN
       NXSTEPL=NXSTEP
    ELSE
       NXSTEPL=1
    ENDIF
    IF(PRESENT(NYMIN)) THEN
       NYMINL=NYMIN
    ELSE
       NYMINL=1
    ENDIF
    IF(PRESENT(NYSTEP)) THEN
       NYSTEPL=NYSTEP
    ELSE
       NYSTEPL=1
    ENDIF

    CALL GMNMX2(GY,NXM,NXMINL,NXMAX,NXSTEPL,NYMINL,NYMAX,NYSTEPL,GYMIN,GYMAX)
    GYNORM=ABS(GYMIN)+ABS(GYMAX)
    IF(ABS(GYMAX-GYMIN).LT.1.E-6*GYNORM) THEN
       GYMIN=GYMIN-0.999E-6*GYNORM
       GYMAX=GYMAX+1.000E-6*GYNORM
    ENDIF

    RETURN
  END SUBROUTINE GRFUT2

! ***** EVALUATE GSMIN,GSMAX,GSTEP,GORG *****

  SUBROUTINE GRFUT3(GMIN,GMAX,GSMIN,GSMAX,GSTEP,GORG)

    IMPLICIT NONE
    REAL,INTENT(IN):: GMIN,GMAX
    REAL,INTENT(OUT):: GSMIN,GSMAX,GSTEP,GORG
    REAL(dp):: TEMP

    CALL GQSCAL(GMIN,GMAX,GSMIN,GSMAX,GSTEP)

    IF(GMIN.LE.0.0.AND.GMAX.GE.0.0) THEN
       GORG=0.0
    ELSE
       IF(ABS(GSTEP).GT.1.E-16) THEN
          TEMP=GSMIN/(2*GSTEP)
          IF(TEMP.LT.-1.D8) THEN
             GORG=GSMIN
          ELSE IF(TEMP.GT.1.D8) THEN
             GORG=GSMIN
          ELSE IF(TEMP.GT.0.D0) THEN
             GORG=(INT(TEMP)+1)*2*GSTEP
          ELSE
             GORG=(INT(TEMP))*2*GSTEP
          ENDIF
          IF(GORG.LE.GSMIN) GORG=GORG+2*GSTEP   ! 2020-09-14 AF: 
          IF(GORG.GE.GSMAX) GORG=GORG-2*GSTEP   !   cope with no scale case
       ELSE
          GORG=GSMIN
       END IF
    ENDIF
!    WRITE(6,'(A,6ES12.4)') 'GRFUT3:',GMIN,GMAX,GSMIN,GSMAX,GSTEP,GORG
    RETURN
  END SUBROUTINE GRFUT3


  SUBROUTINE GRFUT4(NGP,GP)

    IMPLICIT NONE

    INTEGER:: NGP
    REAL:: GP(4)
    INTEGER:: IX,IY
    REAL:: X0,XL,Y0,YL,XS,YS

    IF(NGP.EQ.0) THEN
       GP(1)= 4.0
       GP(2)=24.0
       GP(3)= 2.0
       GP(4)=17.0
    ELSE
       IF(NGP.LE.4) THEN ! 4 figs
          X0= 2.0     ! 3.0
          XL=12.0     ! 12.0
          XS=12.8
          Y0= 1.0     ! 1.75
          YL= 8.5     ! 8.5
          YS= 9.0     ! 8.5
          IX=MOD(NGP-1,2)
          IY=1-(NGP-1)/2
       ELSEIF(NGP.LE.13) THEN ! 9 figs
          X0= 2.0*0.666
          XL=12.0*0.666
          XS=12.8*0.666
          Y0= 1.0*0.666
          YL= 8.5*0.666
          YS= 9.0*0.666
          IX=MOD(NGP-5,3)
          IY=2-(NGP-5)/3
       ELSEIF(NGP.LE.29) THEN ! 16 figs
          X0= 2.0*0.5
          XL=12.0*0.5
          XS=12.8*0.5
          Y0= 1.0*0.5
          YL= 8.5*0.5
          YS= 9.0*0.5
          IX=MOD(NGP-14,4)
          IY=3-(NGP-14)/4
       ELSEIF(NGP.LE.54) THEN ! 25 figs
          X0= 2.0*0.4
          XL=12.0*0.4
          XS=12.8*0.4
          Y0= 1.0*0.4
          YL= 8.5*0.4
          YS= 9.0*0.4
          IX=MOD(NGP-30,5)
          IY=4-(NGP-30)/5
       ELSEIF(NGP.LE.90) THEN ! 36 figs
          X0= 2.0*0.333
          XL=12.0*0.333
          XS=12.8*0.333
          Y0= 1.0*0.333
          YL= 8.5*0.333
          YS= 9.0*0.4
          IX=MOD(NGP-30,5)
          IY=4-(NGP-30)/5
       ENDIF
       GP(1)=IX*XS+X0
       GP(2)=IX*XS+XL
       GP(3)=IY*YS+Y0
       GP(4)=IY*YS+YL
    ENDIF
    RETURN
  END SUBROUTINE GRFUT4

! ***** SET RGB BY ARRAY *****

  SUBROUTINE SETRGBA(RGB)

    IMPLICIT NONE
    REAL,DIMENSION(3),INTENT(IN):: RGB

    CALL SETRGB(RGB(1),RGB(2),RGB(3))
    RETURN
  END SUBROUTINE SETRGBA

!     ****** AVOID REAL*4 UNDERFLOW ******
  
  FUNCTION gsclip(s)

    REAL,INTENT(IN):: s
    REAL:: gsclip

    IF(ABS(s).LT.1.E-30) then
       gsclip=0.0
    ELSEIF(s.GT. 1.E30) THEN
       gsclip= 1.E30
    ELSEIF(s.LT.-1.D30) THEN
       gsclip=-1.E30
    ELSE
       gsclip=s
    ENDIF
    RETURN
  END FUNCTION gsclip

  !   ***** OPTIMUM NUM LENGTH FOR GVALUE *****

  FUNCTION ngslen(step)

    REAL,INTENT(IN):: step
    INTEGER:: ngslen
    REAL:: gxl,gx
    INTEGER:: ngx
    
    gxl=LOG10(step*0.11)
    IF(gxl.LT.0.0) THEN
       ngx = -INT(-gxl)
    ELSE
       ngx =  INT( gxl)
    ENDIF
    IF(ngx.LT.-4.OR.ngx.GT.4)THEN
       gx = gxl-ngx
       IF(gx.LT.0.D0) gx=gx+1.0
       IF(gx.LE.0.15) THEN
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

END MODULE grfutils
