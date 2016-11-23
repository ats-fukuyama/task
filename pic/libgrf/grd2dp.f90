MODULE grd2dp_mod

  PRIVATE
  PUBLIC grd2dp

CONTAINS

! ***** DRAW 2D GRAPH *****

  SUBROUTINE GRD2DP(NGP,FX,FY,FF,NFMAX, &   ! indispensable arguments
                   TITLE,MODE_LS,MODE_PRD,MODE_2D, & ! optional positional args
                   GPXMIN,GPXMAX,GPYMIN,GPYMAX, & ! top of optional args
                   XMIN,XMAX,XORG,YMIN,YMAX,YORG,FMIN,FMAX,FORG, &
                   XSPACE_FACTOR,YSPACE_FACTOR, &
                   NXMIN,NYMIN,NXSTEP,NYSTEP,NLMAX,NMMAX, &
                   LINE_VALUE,LINE_RGB,LINE_THICKNESS,LINE_PAT, &
                   LINE_MARK,LINE_MARK_STEP,LINE_MARK_SIZE, &
                   LINE_MARK_THICKNESS, &
                   MARK_VALUE,MARK_RGB, &
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
                   VALUE_FONT,VALUE_SIZE,VALUE_RGB, &
                   XVALUE_STEP,YVALUE_STEP,FVALUE_STEP, &
                   XVALUE_POS,YVALUE_POS,FVALUE_POS, &
                   XVALUE_TYPE,YVALUE_TYPE,FVALUE_TYPE, &
                   XVALUE_LTYPE,YVALUE_LTYPE,FVALUE_LTYPE, &
                   MODE_XY,MODE_SPL,FRAME_TYPE, &
                   NOTITLE,NOFRAME,NOINFO, &
                   NOXSCALE,NOYSCALE,NOFSCALE, &
                   NOXVALUE,NOYVALUE,NOFVALUE, &
                   LINE_RGB_SUB,MARK_RGB_SUB)

    USE grftype, ONLY: grf_attr_type
    USE grdpconvert, ONLY: grdp_convert
    USE grf2dpexec, ONLY: grf2dp_exec

    IMPLICIT NONE

    INTEGER,INTENT(IN):: NGP                       ! Graph size and postion ID
    REAL(8),INTENT(IN):: FX(NFMAX),FY(NFMAX)       ! 1D data
    REAL(8),INTENT(IN):: FF(NFMAX)                 ! 1D data
    INTEGER,INTENT(IN):: NFMAX                     ! Number of data

    REAL(8),INTENT(IN),OPTIONAL:: GPXMIN,GPXMAX,GPYMIN,GPYMAX
    REAL(8),INTENT(IN),OPTIONAL:: XMIN,XMAX,XORG
    REAL(8),INTENT(IN),OPTIONAL:: YMIN,YMAX,YORG
    REAL(8),INTENT(IN),OPTIONAL:: FMIN,FMAX,FORG
    REAL(8),INTENT(IN),OPTIONAL:: XSPACE_FACTOR,YSPACE_FACTOR
    INTEGER,INTENT(IN),OPTIONAL:: NXMIN,NYMIN,NXSTEP,NYSTEP

    INTEGER,INTENT(IN),OPTIONAL:: NLMAX,NMMAX
    REAL(8),INTENT(IN),OPTIONAL:: LINE_VALUE(:),LINE_RGB(:,:)
    REAL(8),INTENT(IN),OPTIONAL:: LINE_THICKNESS(:)
    INTEGER,INTENT(IN),OPTIONAL:: LINE_PAT(:)
    INTEGER,INTENT(IN),OPTIONAL:: LINE_MARK(:),LINE_MARK_STEP(:)
    REAL(8),INTENT(IN),OPTIONAL:: LINE_MARK_SIZE(:)
    REAL(8),INTENT(IN),OPTIONAL:: LINE_MARK_THICKNESS(:)
    REAL(8),INTENT(IN),OPTIONAL:: MARK_VALUE(:),MARK_RGB(:,:)
    REAL(8),INTENT(IN),OPTIONAL:: ASPECT,BEV_XLEN,BEV_YLEN,BEV_ZLEN
    REAL(8),INTENT(IN),OPTIONAL:: BEV_PHI,BEV_CHI,BEV_DISTANCE
    INTEGER,INTENT(IN),OPTIONAL:: BEV_TYPE

    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: TITLE,XTITLE,YTITLE
    INTEGER,INTENT(IN),OPTIONAL:: TITLE_LEN,XTITLE_LEN,YTITLE_LEN
    REAL(8),INTENT(IN),OPTIONAL:: TITLE_SIZE,XTITLE_SIZE,YTITLE_SIZE
    INTEGER,INTENT(IN),OPTIONAL:: TITLE_FONT,XTITLE_FONT,YTITLE_FONT
    REAL(8),INTENT(IN),OPTIONAL:: TITLE_RGB(3),XTITLE_RGB(3),YTITLE_RGB(3)
    REAL(8),INTENT(IN),OPTIONAL:: TITLE_SEP,XTITLE_SEP,YTITLE_SEP
    INTEGER,INTENT(IN),OPTIONAL:: TITLE_POS

    REAL(8),INTENT(IN),OPTIONAL:: FRAME_RGB(3),FRAME_THICKNESS

    REAL(8),INTENT(IN),OPTIONAL:: SCALE_RGB(3),SCALE_Zero_RGB(3)
    REAL(8),INTENT(IN),OPTIONAL:: SCALE_THICKNESS,SCALE_ZERO_THICKNESS
    REAL(8),INTENT(IN),OPTIONAL:: XSCALE_STEP,YSCALE_STEP,FSCALE_STEP
    REAL(8),INTENT(IN),OPTIONAL:: XSCALE_SIZE,YSCALE_SIZE,FSCALE_SIZE
    INTEGER,INTENT(IN),OPTIONAL:: XSCALE_TYPE,YSCALE_TYPE,FSCALE_TYPE
    INTEGER,INTENT(IN),OPTIONAL:: XSCALE_LTYPE,YSCALE_LTYPE,FSCALE_LTYPE
    INTEGER,INTENT(IN),OPTIONAL:: XSCALE_ZERO,YSCALE_ZERO,FSCALE_ZERO
    INTEGER,INTENT(IN),OPTIONAL:: XSCALE_ZERO_PAT,YSCALE_ZERO_PAT, &
                                  FSCALE_ZERO_PAT

    INTEGER,INTENT(IN),OPTIONAL:: VALUE_FONT
    REAL(8),INTENT(IN),OPTIONAL:: VALUE_SIZE,VALUE_RGB(3)
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_STEP,YVALUE_STEP,FVALUE_STEP
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_POS,YVALUE_POS,FVALUE_POS
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_TYPE,YVALUE_TYPE,FVALUE_TYPE
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_LTYPE,YVALUE_LTYPE,FVALUE_LTYPE

    INTEGER,INTENT(IN),OPTIONAL:: MODE_2D,MODE_XY,MODE_PRD,MODE_LS
    INTEGER,INTENT(IN),OPTIONAL:: MODE_SPL,FRAME_TYPE
    INTEGER,INTENT(IN),OPTIONAL:: NOTITLE,NOFRAME,NOINFO
    INTEGER,INTENT(IN),OPTIONAL:: NOXSCALE,NOYSCALE,NOFSCALE
    INTEGER,INTENT(IN),OPTIONAL:: NOXVALUE,NOYVALUE,NOFVALUE

    INTERFACE
       SUBROUTINE LINE_RGB_SUB(VALUE,RGB)
         REAL(4),INTENT(IN):: VALUE
         REAL(4),DIMENSION(3),INTENT(OUT)::RGB
       END SUBROUTINE LINE_RGB_SUB
       SUBROUTINE MARK_RGB_SUB(VALUE,RGB)
         REAL(4),INTENT(IN):: VALUE
         REAL(4),DIMENSION(3),INTENT(OUT)::RGB
       END SUBROUTINE MARK_RGB_SUB
    END INTERFACE
    OPTIONAL LINE_RGB_SUB,MARK_RGB_SUB

    TYPE(grf_attr_type):: A

    INTEGER:: NF,NM
    REAL(4):: GX(NFMAX),GY(NFMAX),GF(NFMAX)
    REAL(4):: GUCLIP
    REAL(4),DIMENSION(:),ALlOCATABLE:: WLM

    IF(PRESENT(MODE_2D)) THEN
       A%MODE_2D=MODE_2D
    ELSE
       A%MODE_2D=1
    ENDIF

    IF(PRESENT(LINE_MARK_THICKNESS)) THEN
       IF(NMMAX.GE.1) THEN
          ALLOCATE(WLM(NMMAX))
          DO NM=1,NMMAX
             WLM(NM)=LINE_MARK_THICKNESS(NM)
          END DO
       ELSE
          ALLOCATE(WLM(1))
          WLM(1)=0.035
       END IF
    ELSE
       ALLOCATE(WLM(1))
       WLM(1)=0.035
    ENDIF

    DO NF=1,NFMAX
       GX(NF)=GUCLIP(FX(NF))
       GY(NF)=GUCLIP(FY(NF))
       GF(NF)=GUCLIP(FF(NF))
    ENDDO

    CALL GRDP_CONVERT(NGP,GX,GY,GF,NFMAX, &
                     GPXMIN,GPXMAX,GPYMIN,GPYMAX, &
                     XMIN,XMAX,XORG,YMIN,YMAX,YORG,FMIN,FMAX,FORG, &
                     XSPACE_FACTOR,YSPACE_FACTOR, &
                     NXMIN,NYMIN,NXSTEP,NYSTEP,NLMAX,NMMAX, &
                     LINE_VALUE,LINE_RGB,LINE_THICKNESS,LINE_PAT, &
                     LINE_MARK,LINE_MARK_STEP,LINE_MARK_SIZE, &
                     MARK_VALUE,MARK_RGB, &
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
                     VALUE_FONT,VALUE_SIZE,VALUE_RGB, &
                     XVALUE_STEP,YVALUE_STEP,FVALUE_STEP, &
                     XVALUE_POS,YVALUE_POS,FVALUE_POS, &
                     XVALUE_TYPE,YVALUE_TYPE,FVALUE_TYPE, &
                     XVALUE_LTYPE,YVALUE_LTYPE,FVALUE_LTYPE, &
                     MODE_2D,MODE_XY,MODE_PRD,MODE_SPL,MODE_LS,FRAME_TYPE, &
                     NOTITLE,NOFRAME,NOINFO, &
                     NOXSCALE,NOYSCALE,NOFSCALE, &
                     NOXVALUE,NOYVALUE,NOFVALUE, &
                     LINE_RGB_SUB,MARK_RGB_SUB,A)

       CALL GRF2DP_EXEC(GX,GY,GF,NFMAX,A,WLM, &
            LINE_RGB_SUB=LINE_RGB_SUB,MARK_RGB_SUB=MARK_RGB_SUB)

    IF(ALLOCATED(A%LINE_VALUE))     DEALLOCATE(A%LINE_VALUE)
    IF(ALLOCATED(A%LINE_RGB))       DEALLOCATE(A%LINE_RGB)
    IF(ALLOCATED(A%LINE_THICKNESS)) DEALLOCATE(A%LINE_THICKNESS)
    IF(ALLOCATED(A%LINE_PAT))       DEALLOCATE(A%LINE_PAT)
    IF(ALLOCATED(A%LINE_MARK))      DEALLOCATE(A%LINE_MARK)
    IF(ALLOCATED(A%LINE_MARK_STEP)) DEALLOCATE(A%LINE_MARK_STEP)
    IF(ALLOCATED(A%LINE_MARK_SIZE)) DEALLOCATE(A%LINE_MARK_SIZE)
    IF(ALLOCATED(A%PAINT_VALUE))    DEALLOCATE(A%PAINT_VALUE)
    IF(ALLOCATED(A%PAINT_RGB))      DEALLOCATE(A%PAINT_RGB)
    DEALLOCATE(WLM)
    RETURN
  END SUBROUTINE GRD2DP

END MODULE grd2dp_mod
