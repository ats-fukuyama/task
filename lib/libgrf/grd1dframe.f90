MODULE grd1dframe

  USE task_kinds,ONLY: dp
  USE grftype, ONLY: grf_attr_type
  PRIVATE
  PUBLIC grd1d_frame_start,grd1d_frame_end

  TYPE(grf_attr_type):: A

CONTAINS

! ***** DRAW 1D GRAPH *****

  SUBROUTINE GRD1D_FRAME_START(NGP,XMIN,XMAX,FMIN,FMAX, &
                   TITLE, & ! optional positional args
                   GPXMIN,GPXMAX,GPYMIN,GPYMAX, & ! top of optional args
                   XORG,YORG,YMIN,YMAX,FORG, &
                   XSPACE_FACTOR,YSPACE_FACTOR, &
                   NXMIN,NYMIN,NXSTEP,NYSTEP,NLMAX,NPMAX, &
                   LINE_VALUE,LINE_RGB,LINE_THICKNESS,LINE_PAT, &
                   LINE_MARK,LINE_MARK_STEP,LINE_MARK_SIZE, &
                   PAINT_VALUE,PAINT_RGB, &
                   ASPECT,BEV_XLEN,BEV_YLEN,BEV_ZLEN, &
                   BEV_PHI,BEV_CHI,BEV_DISTANCE,BEV_TYPE, &
                   XTITLE,YTITLE, &
                   TITLE_LEN,XTITLE_LEN,YTITLE_LEN, &
                   TITLE_SIZE,XTITLE_SIZE,YTITLE_SIZE, &
                   TITLE_FONT,XTITLE_FONT,YTITLE_FONT, &
                   TITLE_RGB,XTITLE_RGB,YTITLE_RGB, &
                   TITLE_SEP,XTITLE_SEP,YTITLE_SEP,TITLE_POS, &
                   FRAME_RGB,FRAME_THICKNESS, &
                   SCALE_RGB,SCALE_ZERO_RGB, &
                   SCALE_THICKNESS,SCALE_ZERO_THICKNESS, &
                   XSCALE_STEP,YSCALE_STEP,FSCALE_STEP, &
                   XSCALE_SIZE,YSCALE_SIZE,FSCALE_SIZE, &
                   XSCALE_TYPE,YSCALE_TYPE,FSCALE_TYPE, &
                   XSCALE_LTYPE,YSCALE_LTYPE,FSCALE_LTYPE, &
                   XSCALE_ZERO,YSCALE_ZERO,FSCALE_ZERO, &
                   XSCALE_ZERO_PAT,YSCALE_ZERO_PAT,FSCALE_ZERO_PAT, &
                   XGRID_STEP,YGRID_STEP,FGRID_STEP, &
                   XGRID_LTYPE,YGRID_LTYPE,FGRID_LTYPE, &
                   VALUE_FONT,VALUE_SIZE,VALUE_RGB, &
                   XVALUE_STEP,YVALUE_STEP,FVALUE_STEP, &
                   XVALUE_POS,YVALUE_POS,FVALUE_POS, &
                   XVALUE_TYPE,YVALUE_TYPE,FVALUE_TYPE, &
                   XVALUE_LTYPE,YVALUE_LTYPE,FVALUE_LTYPE, &
                   MODE_XY,MODE_SPL,FRAME_TYPE, &
                   NOTITLE,NOFRAME,NOINFO, &
                   NOXSCALE,NOYSCALE,NOFSCALE, &
                   NOXVALUE,NOYVALUE,NOFVALUE, &
                   LINE_RGB_SUB,PAINT_RGB_SUB)

    USE grdconvert, ONLY: grd_convert

    IMPLICIT NONE

    INTEGER,INTENT(IN):: NGP                       ! Graph size and postion ID
    REAL(dp),INTENT(IN):: XMIN,XMAX
    REAL(dp),INTENT(IN):: FMIN,FMAX

    REAL(dp),INTENT(IN),OPTIONAL:: GPXMIN,GPXMAX,GPYMIN,GPYMAX
    REAL(dp),INTENT(IN),OPTIONAL:: XORG
    REAL(dp),INTENT(IN),OPTIONAL:: FORG
    REAL(dp),INTENT(IN),OPTIONAL:: YMIN,YMAX,YORG
    REAL(dp),INTENT(IN),OPTIONAL:: XSPACE_FACTOR,YSPACE_FACTOR
    INTEGER,INTENT(IN),OPTIONAL:: NXMIN,NYMIN,NXSTEP,NYSTEP

    INTEGER,INTENT(IN),OPTIONAL:: NLMAX,NPMAX
    REAL(dp),INTENT(IN),OPTIONAL:: LINE_VALUE(:),LINE_RGB(:,:)
    REAL(dp),INTENT(IN),OPTIONAL:: LINE_THICKNESS(:)
    INTEGER,INTENT(IN),OPTIONAL:: LINE_PAT(:)
    INTEGER,INTENT(IN),OPTIONAL:: LINE_MARK(:),LINE_MARK_STEP(:)
    REAL(dp),INTENT(IN),OPTIONAL:: LINE_MARK_SIZE(:)
    REAL(dp),INTENT(IN),OPTIONAL:: PAINT_VALUE(:),PAINT_RGB(:,:)
    REAL(dp),INTENT(IN),OPTIONAL:: ASPECT,BEV_XLEN,BEV_YLEN,BEV_ZLEN
    REAL(dp),INTENT(IN),OPTIONAL:: BEV_PHI,BEV_CHI,BEV_DISTANCE
    INTEGER,INTENT(IN),OPTIONAL:: BEV_TYPE

    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: TITLE,XTITLE,YTITLE
    INTEGER,INTENT(IN),OPTIONAL:: TITLE_LEN,XTITLE_LEN,YTITLE_LEN
    REAL(dp),INTENT(IN),OPTIONAL:: TITLE_SIZE,XTITLE_SIZE,YTITLE_SIZE
    INTEGER,INTENT(IN),OPTIONAL:: TITLE_FONT,XTITLE_FONT,YTITLE_FONT
    REAL(dp),INTENT(IN),OPTIONAL:: TITLE_RGB(3),XTITLE_RGB(3),YTITLE_RGB(3)
    REAL(dp),INTENT(IN),OPTIONAL:: TITLE_SEP,XTITLE_SEP,YTITLE_SEP
    INTEGER,INTENT(IN),OPTIONAL:: TITLE_POS

    REAL(dp),INTENT(IN),OPTIONAL:: FRAME_RGB(3),FRAME_THICKNESS

    REAL(dp),INTENT(IN),OPTIONAL:: SCALE_RGB(3),SCALE_Zero_RGB(3)
    REAL(dp),INTENT(IN),OPTIONAL:: SCALE_THICKNESS,SCALE_ZERO_THICKNESS
    REAL(dp),INTENT(IN),OPTIONAL:: XSCALE_STEP,YSCALE_STEP,FSCALE_STEP
    REAL(dp),INTENT(IN),OPTIONAL:: XSCALE_SIZE,YSCALE_SIZE,FSCALE_SIZE
    INTEGER,INTENT(IN),OPTIONAL:: XSCALE_TYPE,YSCALE_TYPE,FSCALE_TYPE
    INTEGER,INTENT(IN),OPTIONAL:: XSCALE_LTYPE,YSCALE_LTYPE,FSCALE_LTYPE
    INTEGER,INTENT(IN),OPTIONAL:: XSCALE_ZERO,YSCALE_ZERO,FSCALE_ZERO
    INTEGER,INTENT(IN),OPTIONAL:: XSCALE_ZERO_PAT,YSCALE_ZERO_PAT, &
                                  FSCALE_ZERO_PAT
    REAL(dp),INTENT(IN),OPTIONAL:: XGRID_STEP,YGRID_STEP,FGRID_STEP
    INTEGER,INTENT(IN),OPTIONAL:: XGRID_LTYPE,YGRID_LTYPE,FGRID_LTYPE

    INTEGER,INTENT(IN),OPTIONAL:: VALUE_FONT
    REAL(dp),INTENT(IN),OPTIONAL:: VALUE_SIZE,VALUE_RGB(3)
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_STEP,YVALUE_STEP,FVALUE_STEP
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_POS,YVALUE_POS,FVALUE_POS
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_TYPE,YVALUE_TYPE,FVALUE_TYPE
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_LTYPE,YVALUE_LTYPE,FVALUE_LTYPE

    INTEGER,INTENT(IN),OPTIONAL:: MODE_XY
    INTEGER,INTENT(IN),OPTIONAL:: MODE_SPL,FRAME_TYPE
    INTEGER,INTENT(IN),OPTIONAL:: NOTITLE,NOFRAME,NOINFO
    INTEGER,INTENT(IN),OPTIONAL:: NOXSCALE,NOYSCALE,NOFSCALE
    INTEGER,INTENT(IN),OPTIONAL:: NOXVALUE,NOYVALUE,NOFVALUE
    INTEGER:: MODE_2D,MODE_PRD,MODE_LS

    INTERFACE
       SUBROUTINE LINE_RGB_SUB(VALUE,RGB)
         REAL,INTENT(IN):: VALUE
         REAL,DIMENSION(3),INTENT(OUT)::RGB
       END SUBROUTINE LINE_RGB_SUB
       SUBROUTINE PAINT_RGB_SUB(VALUE,RGB)
         REAL,INTENT(IN):: VALUE
         REAL,DIMENSION(3),INTENT(OUT):: RGB
       END SUBROUTINE PAINT_RGB_SUB
    END INTERFACE
    OPTIONAL LINE_RGB_SUB,PAINT_RGB_SUB

    INTEGER:: NX,NY,NXM,NXMAX,NGMAX
    REAL,DIMENSION(:),ALLOCATABLE:: GX,GY
    REAL,DIMENSION(:,:),ALLOCATABLE:: GF
    REAL:: GUCLIP

    NXM=2
    NXMAX=2
    NGMAX=1
    ALLOCATE(GX(NXMAX),GY(NGMAX),GF(NXM,NGMAX))
    A%MODE_2D=0

    GX(1)=GUCLIP(XMIN)
    GX(2)=GUCLIP(XMAX)
    GY(1)=0.0
    GF(1,1)=GUCLIP(FMIN)
    GF(2,1)=GUCLIP(FMAX)
    MODE_LS=0

    CALL GRD_CONVERT(NGP,GX,GY,GF,NXM,NXMAX,NGMAX, &
                     GPXMIN,GPXMAX,GPYMIN,GPYMAX, &
                     XMIN,XMAX,XORG,YMIN,YMAX,YORG, &
                     FMIN,FMAX,FORG, &
                     XSPACE_FACTOR,YSPACE_FACTOR, &
                     NXMIN,NYMIN,NXSTEP,NYSTEP,NLMAX,NPMAX, &
                     LINE_VALUE,LINE_RGB,LINE_THICKNESS,LINE_PAT, &
                     LINE_MARK,LINE_MARK_STEP,LINE_MARK_SIZE, &
                     PAINT_VALUE,PAINT_RGB, &
                     ASPECT,BEV_XLEN,BEV_YLEN,BEV_ZLEN, &
                     BEV_PHI,BEV_CHI,BEV_DISTANCE,BEV_TYPE, &
                     TITLE,XTITLE,YTITLE, &
                     TITLE_LEN,XTITLE_LEN,YTITLE_LEN, &
                     TITLE_SIZE,XTITLE_SIZE,YTITLE_SIZE, &
                     TITLE_FONT,XTITLE_FONT,YTITLE_FONT, &
                     TITLE_RGB,XTITLE_RGB,YTITLE_RGB, &
                     TITLE_SEP,XTITLE_SEP,YTITLE_SEP,TITLE_POS, &
                     FRAME_RGB,FRAME_THICKNESS, &
                     SCALE_RGB,SCALE_ZERO_RGB, &
                     SCALE_THICKNESS,SCALE_ZERO_THICKNESS, &
                     XSCALE_STEP,YSCALE_STEP,FSCALE_STEP, &
                     XSCALE_SIZE,YSCALE_SIZE,FSCALE_SIZE, &
                     XSCALE_TYPE,YSCALE_TYPE,FSCALE_TYPE, &
                     XSCALE_LTYPE,YSCALE_LTYPE,FSCALE_LTYPE, &
                     XSCALE_ZERO,YSCALE_ZERO,FSCALE_ZERO, &
                     XSCALE_ZERO_PAT,YSCALE_ZERO_PAT,FSCALE_ZERO_PAT, &
                     XGRID_STEP,YGRID_STEP,FGRID_STEP, &
                     XGRID_LTYPE,YGRID_LTYPE,FGRID_LTYPE, &
                     VALUE_FONT,VALUE_SIZE,VALUE_RGB, &
                     XVALUE_STEP,YVALUE_STEP,FVALUE_STEP, &
                     XVALUE_POS,YVALUE_POS,FVALUE_POS, &
                     XVALUE_TYPE,YVALUE_TYPE,FVALUE_TYPE, &
                     XVALUE_LTYPE,YVALUE_LTYPE,FVALUE_LTYPE, &
                     MODE_2D,MODE_XY,MODE_PRD,MODE_SPL,MODE_LS,FRAME_TYPE, &
                     NOTITLE,NOFRAME,NOINFO, &
                     NOXSCALE,NOYSCALE,NOFSCALE, &
                     NOXVALUE,NOYVALUE,NOFVALUE, &
                     LINE_RGB_SUB,PAINT_RGB_SUB,A)

    CALL GRF1D_FRAME_SETUP
    DEALLOCATE(GX,GY,GF)
    RETURN
  END SUBROUTINE GRD1D_FRAME_START

  SUBROUTINE GRD1D_FRAME_END

    USE grfutils,ONLY: grf_title,grf_frame1d,grf_info
    IMPLICIT NONE

    CALL OFFCLP

    IF(A%NOFRAME.EQ.0) CALL GRF_FRAME1D(A)
    IF(A%NOTITLE.EQ.0) CALL GRF_TITLE(A)
    IF(A%NOINFO.EQ.0)  CALL GRF_INFO(A)

    CALL SETRGB(0.0,0.0,0.0)
    IF(ALLOCATED(A%LINE_VALUE))     DEALLOCATE(A%LINE_VALUE)
    IF(ALLOCATED(A%LINE_RGB))       DEALLOCATE(A%LINE_RGB)
    IF(ALLOCATED(A%LINE_THICKNESS)) DEALLOCATE(A%LINE_THICKNESS)
    IF(ALLOCATED(A%LINE_PAT))       DEALLOCATE(A%LINE_PAT)
    IF(ALLOCATED(A%LINE_MARK))      DEALLOCATE(A%LINE_MARK)
    IF(ALLOCATED(A%LINE_MARK_STEP)) DEALLOCATE(A%LINE_MARK_STEP)
    IF(ALLOCATED(A%LINE_MARK_SIZE)) DEALLOCATE(A%LINE_MARK_SIZE)
    IF(ALLOCATED(A%PAINT_VALUE))    DEALLOCATE(A%PAINT_VALUE)
    IF(ALLOCATED(A%PAINT_RGB))      DEALLOCATE(A%PAINT_RGB)
    RETURN
  END SUBROUTINE GRD1D_FRAME_END

! ***** DRAW 1D GRAPH FRAME *****

  SUBROUTINE GRF1D_FRAME_SETUP

!   *** after calling this and drawing, callin OFFCLP is required.

    IMPLICIT NONE

    CALL GDEFIN(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX, &
                A%XMIN,A%XMAX,A%FMIN,A%FMAX)
    CALL SETCLP(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX)
    RETURN
  END SUBROUTINE GRF1D_FRAME_SETUP

END MODULE grd1dframe

