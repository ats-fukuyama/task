MODULE grdconvert

  PRIVATE
  PUBLIC grd_convert

CONTAINS

! ***** Convert input parameters to A: REAL(4) *****

  SUBROUTINE GRD_CONVERT(NGP,GX,GY,GF,NXM,NXMAX,NYMAX, &
                         GPXMIN,GPXMAX,GPYMIN,GPYMAX, &
                         XMIN,XMAX,XORG,YMIN,YMAX,YORG,FMIN,FMAX,FORG, &
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
                         MODE_2D,MODE_XY,MODE_PRD,MODE_SPL, &
                         MODE_LS,FRAME_TYPE, &
                         NOTITLE,NOFRAME,NOINFO, &
                         NOXSCALE,NOYSCALE,NOFSCALE, &
                         NOXVALUE,NOYVALUE,NOFVALUE, &
                         LINE_RGB_SUB,PAINT_RGB_SUB,A)

    USE grftype, ONLY: grf_attr_type
    USE grfutils, ONLY: grfut1,grfut2,grfut3,grfut4
    IMPLICIT NONE

    INTEGER,INTENT(IN):: NGP
    REAL(4),INTENT(IN):: GX(NXMAX),GY(NYMAX),GF(NXM,NYMAX)
    INTEGER,INTENT(IN):: NXM,NXMAX,NYMAX

    REAL(8),INTENT(IN),OPTIONAL:: GPXMIN,GPXMAX,GPYMIN,GPYMAX
    REAL(8),INTENT(IN),OPTIONAL:: XMIN,XMAX,XORG
    REAL(8),INTENT(IN),OPTIONAL:: YMIN,YMAX,YORG
    REAL(8),INTENT(IN),OPTIONAL:: FMIN,FMAX,FORG
    REAL(8),INTENT(IN),OPTIONAL:: XSPACE_FACTOR,YSPACE_FACTOR
    INTEGER,INTENT(IN),OPTIONAL:: NXMIN,NYMIN,NXSTEP,NYSTEP

    INTEGER,INTENT(IN),OPTIONAL:: NLMAX,NPMAX
    REAL(8),INTENT(IN),OPTIONAL:: LINE_VALUE(:),LINE_RGB(:,:)
    REAL(8),INTENT(IN),OPTIONAL:: LINE_THICKNESS(:)
    INTEGER,INTENT(IN),OPTIONAL:: LINE_PAT(:)
    INTEGER,INTENT(IN),OPTIONAL:: LINE_MARK(:),LINE_MARK_STEP(:)
    REAL(8),INTENT(IN),OPTIONAL:: LINE_MARK_SIZE(:)
    REAL(8),INTENT(IN),OPTIONAL:: PAINT_VALUE(:),PAINT_RGB(:,:)
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
    REAL(8),INTENT(IN),OPTIONAL:: XGRID_STEP,YGRID_STEP,FGRID_STEP
    INTEGER,INTENT(IN),OPTIONAL:: XGRID_LTYPE,YGRID_LTYPE,FGRID_LTYPE

    INTEGER,INTENT(IN),OPTIONAL:: VALUE_FONT
    REAL(8),INTENT(IN),OPTIONAL:: VALUE_SIZE,VALUE_RGB(3)
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_STEP,YVALUE_STEP,FVALUE_STEP
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_POS,YVALUE_POS,FVALUE_POS
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_TYPE,YVALUE_TYPE,FVALUE_TYPE
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_LTYPE,YVALUE_LTYPE,FVALUE_LTYPE

    INTEGER,INTENT(IN),OPTIONAL:: MODE_2D,MODE_XY,MODE_PRD,MODE_SPL
    INTEGER,INTENT(IN),OPTIONAL:: MODE_LS,FRAME_TYPE
    INTEGER,INTENT(IN),OPTIONAL:: NOTITLE,NOFRAME,NOINFO
    INTEGER,INTENT(IN),OPTIONAL:: NOXSCALE,NOYSCALE,NOFSCALE
    INTEGER,INTENT(IN),OPTIONAL:: NOXVALUE,NOYVALUE,NOFVALUE

    INTERFACE
       SUBROUTINE LINE_RGB_SUB(VALUE,RGB)
         REAL(4),INTENT(IN):: VALUE
         REAL(4),DIMENSION(3),INTENT(OUT)::RGB
       END SUBROUTINE LINE_RGB_SUB
       SUBROUTINE PAINT_RGB_SUB(VALUE,RGB)
         REAL(4),INTENT(IN):: VALUE
         REAL(4),DIMENSION(3),INTENT(OUT):: RGB
       END SUBROUTINE PAINT_RGB_SUB
       FUNCTION GUCLIP(VALUE)
         REAL(8),INTENT(IN):: VALUE
         REAL(4):: GUCLIP
       END FUNCTION GUCLIP
    END INTERFACE
    OPTIONAL LINE_RGB_SUB,PAINT_RGB_SUB

    TYPE(grf_attr_type):: A

    INTEGER,PARAMETER:: NLM=5
    REAL(4):: GP(4),GL,GFACTOR
    REAL(4):: GXMINL,GXMAXL,GYMINL,GYMAXL,GFMINL,GFMAXL
    REAL(4):: GXMIN,GXMAX,GXSTEP,GXORG
    REAL(4):: GYMIN,GYMAX,GYSTEP,GYORG
    REAL(4):: GFMIN,GFMAX,GFSTEP,GFORG
    REAL(4):: GRGB(3,NLM),FACTOR
    REAL(4):: GPC,GPD
    INTEGER:: NCH,NGULEN,NL,NLL,NP,NPL,I
    INTEGER:: IPAT(NLM)
    CHARACTER(LEN=80):: KTITLE
    CHARACTER(LEN=1):: KDL

    DATA IPAT/0,2,3,4,6/
    DATA GRGB/0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,0.0,1.0, &
              0.0,1.0,0.0, 1.0,0.0,1.0/

    IF(NGP.LT.0.OR.NGP.GT.54) THEN
       WRITE(6,*) 'XX GRD_CONVERT: INVALID NGP'
       RETURN
    ENDIF

    IF(PRESENT(MODE_XY)) THEN
       A%MODE_XY=MODE_XY
    ELSE
       A%MODE_XY=0
    ENDIF
    IF(PRESENT(MODE_PRD)) THEN
       A%MODE_PRD=MODE_PRD
    ELSE
       A%MODE_PRD=0
    ENDIF
    IF(PRESENT(MODE_SPL)) THEN
       A%MODE_SPL=MODE_SPL
    ELSE
       A%MODE_SPL=0
    ENDIF
    IF(PRESENT(MODE_LS)) THEN
       A%MODE_LS=MODE_LS
    ELSE
       A%MODE_LS=0
    ENDIF

!     ----- define graph position -----

    CALL GRFUT4(NGP,GP) ! assign predefined graph position to GP
    IF(PRESENT(GPXMIN)) THEN
       A%GPXMIN=GUCLIP(GPXMIN)
    ELSE
       A%GPXMIN=GP(1)
    ENDIF
    IF(PRESENT(GPXMAX)) THEN
       A%GPXMAX=GUCLIP(GPXMAX)
    ELSE
       A%GPXMAX=GP(2)
    ENDIF
    IF(PRESENT(GPYMIN)) THEN
       A%GPYMIN=GUCLIP(GPYMIN)
    ELSE
       A%GPYMIN=GP(3)
    ENDIF
    IF(PRESENT(GPYMAX)) THEN
       A%GPYMAX=GUCLIP(GPYMAX)
    ELSE
       A%GPYMAX=GP(4)
    ENDIF

    IF(PRESENT(ASPECT)) THEN
       A%ASPECT=GUCLIP(ASPECT)
       IF(A%ASPECT /= 0.0) THEN
          IF(A%ASPECT >= 0.75) THEN  ! 0.75=15/20=(17-2)/(24-4)
             A%GPXMAX=A%GPXMIN+(A%GPYMAX-A%GPYMIN)/A%ASPECT
          ELSE
             A%GPYMAX=A%GPYMIN+(A%GPXMAX-A%GPXMIN)*A%ASPECT
          END IF
       END IF
    ELSE
       A%ASPECT=(A%GPYMAX-A%GPYMIN)/(A%GPXMAX-A%GPXMIN)
    END IF

    SELECT CASE(A%MODE_XY)
    CASE(1)
       IF(A%GPXMAX-A%GPXMIN .GT. A%GPYMAX-A%GPYMIN) THEN
          GPC=0.5*(A%GPXMAX+A%GPXMIN)
          GPD=0.5*(A%GPYMAX-A%GPYMIN)
          A%GPXMIN=GPC-GPD
          A%GPXMAX=GPC+GPD
       ELSE
          GPC=0.5*(A%GPYMAX+A%GPYMIN)
          GPD=0.5*(A%GPXMAX-A%GPXMIN)
          A%GPYMIN=GPC-GPD
          A%GPYMAX=GPC+GPD
       END IF
    CASE(2)
       IF(A%GPXMAX-A%GPXMIN .GT. 2.0*(A%GPYMAX-A%GPYMIN)) THEN
          GPC=0.5*(A%GPXMAX+A%GPXMIN)
          GPD=     A%GPYMAX-A%GPYMIN
          A%GPXMIN=GPC-GPD
          A%GPXMAX=GPC+GPD
       ELSE
          GPC=0.5*(A%GPYMAX+A%GPYMIN)
          GPD=0.25*(A%GPXMAX-A%GPXMIN)
          A%GPYMIN=GPC-GPD
          A%GPYMAX=GPC+GPD
       END IF
    END SELECT

!     ----- define graph size scaling factor -----

    GL=MIN(ABS(A%GPXMAX-A%GPXMIN),ABS(A%GPYMAX-A%GPYMIN))
    GFACTOR=GL/15.0

!     ----- SCALE ADJUSTMENT SECTION -----

    IF(PRESENT(XSPACE_FACTOR)) THEN
       A%XSPACE_FACTOR=GUCLIP(XSPACE_FACTOR)
    ELSE
       A%XSPACE_FACTOR=0.0
    END IF
    IF(PRESENT(YSPACE_FACTOR)) THEN
       A%YSPACE_FACTOR=GUCLIP(YSPACE_FACTOR)
    ELSE
       A%YSPACE_FACTOR=0.0
    END IF

!     ----- DATA SECTION -----

    IF(PRESENT(NXMIN)) THEN
       A%NXMIN=NXMIN
    ELSE
       A%NXMIN=1
    END IF
    IF(PRESENT(NXSTEP)) THEN
       A%NXSTEP=NXSTEP
    ELSE
       A%NXSTEP=1
    END IF

    IF(PRESENT(NYMIN)) THEN
       A%NYMIN=NYMIN
    ELSE
       A%NYMIN=1
    END IF
    IF(PRESENT(NYSTEP)) THEN
       A%NYSTEP=NYSTEP
    ELSE
       A%NYSTEP=1
    END IF

    IF(A%MODE_2D.EQ.0) THEN
       CALL GRFUT1(GX,NXMAX,GXMINL,GXMAXL,NXMIN=A%NXMIN,NXSTEP=A%NXSTEP)
       CALL GRFUT2(GF,NXM,NXMAX,NYMAX,GFMINL,GFMAXL, &
                   NXMIN=A%NXMIN,NXSTEP=A%NXSTEP,NYMIN=A%NYMIN,NYSTEP=A%NYSTEP)
       GYMINL=0.0
       GYMAXL=0.0
    ELSE
       CALL GRFUT1(GX,NXMAX,GXMINL,GXMAXL,NXMIN=A%NXMIN,NXSTEP=A%NXSTEP)
       CALL GRFUT1(GY,NYMAX,GYMINL,GYMAXL,NXMIN=A%NYMIN,NXSTEP=A%NYSTEP)
       CALL GRFUT2(GF,NXM,NXMAX,NYMAX,GFMINL,GFMAXL, &
                   NXMIN=A%NXMIN,NXSTEP=A%NXSTEP,NYMIN=A%NYMIN,NYSTEP=A%NYSTEP)
    ENDIF

    IF(PRESENT(XMIN)) GXMINL=GUCLIP(XMIN)
    IF(PRESENT(XMAX)) GXMAXL=GUCLIP(XMAX)
    IF(PRESENT(YMIN)) GYMINL=GUCLIP(YMIN)
    IF(PRESENT(YMAX)) GYMAXL=GUCLIP(YMAX)
    IF(PRESENT(FMIN)) GFMINL=GUCLIP(FMIN)
    IF(PRESENT(FMAX)) GFMAXL=GUCLIP(FMAX)

    CALL GRFUT3(GXMINL,GXMAXL,GXMIN,GXMAX,GXSTEP,GXORG)
    CALL GRFUT3(GYMINL,GYMAXL,GYMIN,GYMAX,GYSTEP,GYORG)
    CALL GRFUT3(GFMINL,GFMAXL,GFMIN,GFMAX,GFSTEP,GFORG)
    
    SELECT CASE(A%MODE_XY)
    CASE(1)
       GXMIN =-GXMAX
       GYMIN =-GXMAX
       GYMAX = GXMAX
       GYORG = GXORG
       GYSTEP= GXSTEP
       GYMINL=-GXMAXL
       GYMAXL= GXMAXL
    CASE(2)
       GXMIN=-GXMAX
       GYMIN= 0.0
       GYMAX= GXMAX
       GYORG = GXORG
       GYSTEP= GXSTEP
       GYMINL= 0.0
       GYMAXL= GXMAXL
    END SELECT

    IF(PRESENT(XMIN)) THEN
       A%XMIN=GUCLIP(XMIN)
    ELSE
       A%XMIN=GXMIN
    ENDIF
    IF(PRESENT(XMAX)) THEN
       A%XMAX=GUCLIP(XMAX)
    ELSE
       A%XMAX=GXMAX
    ENDIF
    IF(PRESENT(XSCALE_STEP)) THEN
       A%XSCALE_STEP=GUCLIP(XSCALE_STEP)
    ELSE
       A%XSCALE_STEP=GXSTEP
    ENDIF
    IF(PRESENT(XGRID_STEP)) THEN
       A%XGRID_STEP=GUCLIP(XGRID_STEP)
    ELSE
       A%XGRID_STEP=0.0
    ENDIF
    IF(PRESENT(XORG)) THEN
       A%XORG=GUCLIP(XORG)
    ELSE
       A%XORG=GXORG
    ENDIF

    IF(PRESENT(YMIN)) THEN
       A%YMIN=GUCLIP(YMIN)
    ELSE
       A%YMIN=GYMIN
    ENDIF
    IF(PRESENT(YMAX)) THEN
       A%YMAX=GUCLIP(YMAX)
    ELSE
       A%YMAX=GYMAX
    ENDIF
    IF(PRESENT(YSCALE_STEP)) THEN
       A%YSCALE_STEP=GUCLIP(YSCALE_STEP)
    ELSE
       A%YSCALE_STEP=GYSTEP
    ENDIF
    IF(PRESENT(YGRID_STEP)) THEN
       A%YGRID_STEP=GUCLIP(YGRID_STEP)
    ELSE
       A%YGRID_STEP=0.0
    ENDIF
    IF(PRESENT(YORG)) THEN
       A%YORG=GUCLIP(YORG)
    ELSE
       A%YORG=GYORG
    ENDIF

    IF(PRESENT(FMIN)) THEN
       A%FMIN=GUCLIP(FMIN)
    ELSE
       A%FMIN=GFMIN
    ENDIF
    IF(PRESENT(FMAX)) THEN
       A%FMAX=GUCLIP(FMAX)
    ELSE
       A%FMAX=GFMAX
    ENDIF
    IF(PRESENT(FSCALE_STEP)) THEN
       A%FSCALE_STEP=GUCLIP(FSCALE_STEP)
    ELSE
       A%FSCALE_STEP=GFSTEP
    ENDIF
    IF(PRESENT(FGRID_STEP)) THEN
       A%FGRID_STEP=GUCLIP(FGRID_STEP)
    ELSE
       A%FGRID_STEP=0.0
    ENDIF
    IF(PRESENT(FORG)) THEN
       A%FORG=GUCLIP(FORG)
    ELSE
       A%FORG=GFORG
    ENDIF

!     ----- Adjust min max according to space factor at both ends-----

    GL=A%XMAX-A%XMIN
    A%XMAX=A%XMAX+0.5*GL*A%XSPACE_FACTOR
    A%XMIN=A%XMIN-0.5*GL*A%XSPACE_FACTOR
    GL=A%YMAX-A%YMIN
    A%YMAX=A%YMAX+0.5*GL*A%YSPACE_FACTOR
    A%YMIN=A%YMIN-0.5*GL*A%YSPACE_FACTOR

!     ---- Adjust origin to be between min and max 
    
    IF(A%MODE_LS.EQ.0.OR.A%MODE_LS.EQ.2) THEN
       IF(A%XORG.LT.A%XMIN) A%XORG=A%XMIN
       IF(A%XORG.LT.A%XMAX) A%XORG=A%XMAX
    ELSE
       IF(A%XORG.LT.A%XMIN) A%XORG=A%XORG+INT(A%XMIN-A%XORG)+1.0
       IF(A%XORG.LT.A%XMAX) A%XORG=A%XORG-INT(A%XORG-A%XMAX)-1.0
    END IF

    IF(A%MODE_LS.EQ.0.OR.A%MODE_LS.EQ.1) THEN
       IF(A%YORG.LT.A%YMIN) A%YORG=A%YMIN
       IF(A%YORG.LT.A%YMAX) A%YORG=A%YMAX
    ELSE
       IF(A%YORG.LT.A%YMIN) A%YORG=A%YORG+INT(A%YMIN-A%YORG)+1.0
       IF(A%YORG.LT.A%YMAX) A%YORG=A%YORG-INT(A%YORG-A%YMAX)-1.0
    END IF

!     ----- Adjust of graph shape by actual 2D/1D shape -----

    IF(A%MODE_2D >= 1) THEN
       IF(A%ASPECT == 0.0) THEN
          IF(A%XMAX-A%XMIN /= 0.0) THEN
             A%ASPECT=(A%YMAX-A%YMIN)/(A%XMAX-A%XMIN)
             IF(A%ASPECT /= 0.0) THEN
                IF(A%ASPECT >= 0.75) THEN   ! 0.75=15/20=(17-2)/(24-4)
                   A%GPXMAX=A%GPXMIN+(A%GPYMAX-A%GPYMIN)/A%ASPECT
                ELSE
                   A%GPYMAX=A%GPYMIN+(A%GPXMAX-A%GPXMIN)*A%ASPECT
                END IF
             END IF
          END IF
       END IF
    ELSE
       IF(A%ASPECT == 0.0) THEN
          IF(A%XMAX-A%XMIN /= 0.0) THEN
             A%ASPECT=(A%FMAX-A%FMIN)/(A%XMAX-A%XMIN)
             IF(A%ASPECT /= 0.0) THEN
                IF(A%ASPECT >= 0.75) THEN  ! 0.75=15/20=(17-2)/(24-4)
                   A%GPXMAX=A%GPXMIN+(A%GPYMAX-A%GPYMIN)/A%ASPECT
                ELSE
                   A%GPYMAX=A%GPYMIN+(A%GPXMAX-A%GPXMIN)*A%ASPECT
                END IF
             END IF
          END IF
       END IF
    END IF


!     ----- LINE SECTION -----

    SELECT CASE(A%MODE_2D)
    CASE(0,21,22) ! 1D data and X-Y plot
       IF(PRESENT(NLMAX)) THEN
          A%NLMAX=NLMAX
       ELSE
          A%NLMAX=NLM
       ENDIF

       ALLOCATE(A%LINE_RGB(3,A%NLMAX))
       ALLOCATE(A%LINE_THICKNESS(A%NLMAX))
       ALLOCATE(A%LINE_PAT(A%NLMAX))
       ALLOCATE(A%LINE_MARK(A%NLMAX))
       ALLOCATE(A%LINE_MARK_STEP(A%NLMAX))
       ALLOCATE(A%LINE_MARK_SIZE(A%NLMAX))

       IF(PRESENT(LINE_RGB)) THEN
          NLL=SIZE(LINE_RGB,DIM=2)
          DO NL=1,A%NLMAX
             DO I=1,3
                A%LINE_RGB(I,NL)=GUCLIP(LINE_RGB(I,MOD(NL-1,NLL)+1))
             END DO
          ENDDO
       ELSE
          DO NL=1,A%NLMAX
             A%LINE_RGB(1:3,NL)=GRGB(1:3,MOD(NL-1,NLM)+1)
          END DO
       END IF
       IF(PRESENT(LINE_THICKNESS)) THEN
          NLL=SIZE(LINE_THICKNESS,DIM=1)
          DO NL=1,A%NLMAX
             A%LINE_THICKNESS(NL)=GUCLIP(LINE_THICKNESS(MOD(NL-1,NLL)+1)) &
                                 *GFACTOR
          END DO
       ELSE
          A%LINE_THICKNESS(1:A%NLMAX)=0.07*GFACTOR
       ENDIF
       IF(PRESENT(LINE_PAT)) THEN
          NLL=SIZE(LINE_PAT,DIM=1)
          DO NL=1,A%NLMAX
             A%LINE_PAT(NL)=LINE_PAT(MOD(NL-1,NLL)+1)
          END DO
       ELSE
          DO NL=1,A%NLMAX
             A%LINE_PAT(NL)=IPAT(MOD(NL-1,NLM)+1)
          END DO
       ENDIF
       IF(PRESENT(LINE_MARK)) THEN
          NLL=SIZE(LINE_MARK,DIM=1)
          DO NL=1,A%NLMAX
             A%LINE_MARK(NL)=LINE_MARK(MOD(NL-1,NLL)+1)
          END DO
       ELSE
          A%LINE_MARK(1:A%NLMAX)=0
       ENDIF
       IF(PRESENT(LINE_MARK_STEP)) THEN
          NLL=SIZE(LINE_MARK_STEP,DIM=1)
          DO NL=1,A%NLMAX
             A%LINE_MARK_STEP(NL)=LINE_MARK_STEP(MOD(NL-1,NLL)+1)
          END DO
       ELSE
          A%LINE_MARK_STEP(1:A%NLMAX)=0
       ENDIF
       IF(PRESENT(LINE_MARK_SIZE)) THEN
          NLL=SIZE(LINE_MARK_SIZE,DIM=1)
          DO NL=1,A%NLMAX
             A%LINE_MARK_SIZE(NL)=LINE_MARK_SIZE(MOD(NL-1,NLL)+1)*GFACTOR
          END DO
       ELSE
          A%LINE_MARK_SIZE(1:A%NLMAX)=0.3*GFACTOR
       ENDIF

    CASE(1:19) ! 2D contour line and paint and 3D Bird's ey view

!       IF(A%FMIN < 0.0 .AND. A%FMAX > 0.0) THEN
!           A%FMAX=MAX(ABS(A%FMIN),ABS(A%FMAX))
!           A%FMIN=-A%FMAX
!       ENDIF

       IF(PRESENT(NLMAX)) THEN
          A%NLMAX=NLMAX
       ELSE
          A%NLMAX=MIN(NINT(ABS((A%FMAX-A%FMIN)/A%FSCALE_STEP))+2,1000)
          IF(A%FMIN < 0.0 .AND. A%FMAX > 0.0) THEN
!             IF(MOD(A%NLMAX,2) == 0) A%NLMAX=A%NLMAX+1
             IF(MOD(A%NLMAX,2) == 1) A%NLMAX=A%NLMAX+1
          ENDIF
       ENDIF

       IF(PRESENT(NPMAX)) THEN
          A%NPMAX=NPMAX
       ELSE
          A%NPMAX=A%NLMAX
       ENDIF

       ALLOCATE(A%LINE_VALUE(A%NLMAX))
       ALLOCATE(A%LINE_RGB(3,A%NLMAX))
       ALLOCATE(A%LINE_THICKNESS(A%NLMAX))
       ALLOCATE(A%LINE_PAT(A%NLMAX))
       ALLOCATE(A%PAINT_VALUE(A%NPMAX))
       ALLOCATE(A%PAINT_RGB(3,A%NPMAX))

       IF(PRESENT(LINE_VALUE)) THEN
          NLL=SIZE(LINE_VALUE,DIM=1)
          DO NL=1,A%NLMAX
             A%LINE_VALUE(NL)=GUCLIP(LINE_VALUE(MOD(NL-1,NLL)+1))
          END DO
       ELSE
          IF(A%NLMAX == 1) THEN
             FACTOR=0.0
          ELSE
             FACTOR=(A%FMAX-A%FMIN)/(A%NLMAX-1)
          END IF
          DO NL=1,A%NLMAX
             A%LINE_VALUE(NL)=A%FMIN+FACTOR*(NL-1)
          END DO
       ENDIF

       IF(PRESENT(LINE_RGB)) THEN
          NLL=SIZE(LINE_RGB,DIM=2)
          DO NL=1,A%NLMAX
             DO I=1,3
                A%LINE_RGB(I,NL)=GUCLIP(LINE_RGB(I,MOD(NL-1,NLL)+1))
             END DO
          END DO
       ELSE IF(PRESENT(LINE_RGB_SUB)) THEN
          DO NL=1,A%NLMAX
             FACTOR=FLOAT(NL-1)/FLOAT(A%NLMAX-1)
             CALL LINE_RGB_SUB(FACTOR,A%LINE_RGB(1:3,NL))
          END DO
       ELSE
          IF(A%FMIN < 0.0 .AND. A%FMAX > 0.0) THEN
             DO NL=1,A%NLMAX
                IF(A%LINE_VALUE(NL) > 0.0) THEN
                   A%LINE_RGB(1:3,NL)=(/1.0,0.0,0.0/)
                ELSE IF(A%LINE_VALUE(NL) < 0.0) THEN
                   A%LINE_RGB(1:3,NL)=(/0.0,0.0,1.0/)
                ELSE
                   A%LINE_RGB(1:3,NL)=(/0.5,0.5,0.5/)
                END IF
             END DO
          ELSE IF(A%FMIN >= 0.0) THEN
             DO NL=1,A%NLMAX
                A%LINE_RGB(1:3,NL)=(/1.0,0.0,0.0/)
             END DO
          ELSE
             DO NL=1,A%NLMAX
                A%LINE_RGB(1:3,NL)=(/0.0,0.0,1.0/)
             END DO
          ENDIF
       END IF

       IF(PRESENT(PAINT_VALUE)) THEN
          NPL=SIZE(PAINT_VALUE,DIM=1)
          DO NP=1,A%NPMAX
             A%PAINT_VALUE(NP)=GUCLIP(PAINT_VALUE(MOD(NP-1,NPL)+1))
          END DO
       ELSE
          IF(A%NPMAX == 1) THEN
             FACTOR=0.0
          ELSE
             FACTOR=(A%FMAX-A%FMIN)/(A%NPMAX-1)
          END IF
          DO NP=1,A%NPMAX
             A%PAINT_VALUE(NP)=A%FMIN+FACTOR*(NP-0.5)
          END DO
       ENDIF

       IF(PRESENT(PAINT_RGB)) THEN
          NPL=SIZE(PAINT_RGB,DIM=2)
          DO NP=1,A%NPMAX
             DO I=1,3
                A%PAINT_RGB(I,NP)=GUCLIP(PAINT_RGB(I,MOD(NP-1,NPL)+1))
             END DO
          END DO
       ELSE IF(PRESENT(PAINT_RGB_SUB)) THEN
          DO NP=1,A%NPMAX
             FACTOR=FLOAT(NP-1)/FLOAT(A%NPMAX-1)
             CALL PAINT_RGB_SUB(FACTOR,A%PAINT_RGB(1:3,NP))
          END DO
       ELSE
          IF(A%FMIN < 0.0 .AND. A%FMAX > 0.0) THEN
             DO NP=1,A%NPMAX
                FACTOR=FLOAT(NP-1)/FLOAT(A%NPMAX-1)
                CALL R2W2B(FACTOR,A%PAINT_RGB(1:3,NP))
             END DO
          ELSE
             IF(A%FMAX > 0.0) THEN
                DO NP=1,A%NPMAX
                   FACTOR=FLOAT(NP-1)/FLOAT(A%NPMAX-1)
                   CALL R2Y2W(FACTOR,A%PAINT_RGB(1:3,NP))
                END DO
             ELSE
                DO NP=1,A%NPMAX
                   FACTOR=FLOAT(NP-1)/FLOAT(A%NPMAX-1)
                   CALL W2G2B(FACTOR,A%PAINT_RGB(1:3,NP))
                END DO
             END IF
          END IF
       END IF

       IF(PRESENT(LINE_THICKNESS)) THEN
          NLL=SIZE(LINE_THICKNESS,DIM=1)
          DO NL=1,A%NLMAX
             A%LINE_THICKNESS(NL)=GUCLIP(LINE_THICKNESS(MOD(NL-1,NLL)+1)) &
                                 *GFACTOR
          END DO
       ELSE
          A%LINE_THICKNESS(1:A%NLMAX)=0.07*GFACTOR
       ENDIF

       IF(PRESENT(LINE_PAT)) THEN
          NLL=SIZE(LINE_PAT,DIM=1)
          DO NL=1,A%NLMAX
             A%LINE_PAT(NL)=LINE_PAT(MOD(NL-1,NLL)+1)
          END DO
       ELSE
          IF(A%FMIN < 0.0 .AND. A%FMAX > 0.0) THEN
             DO NL=1,A%NLMAX
                IF(A%LINE_VALUE(NL) > 0.0) THEN
                   A%LINE_PAT(NL)=0
                ELSE IF(A%LINE_VALUE(NL) < 0.0) THEN
                   A%LINE_PAT(NL)=3
                ELSE
                   A%LINE_PAT(NL)=5
                END IF
             END DO
          ELSE
             DO NL=1,A%NLMAX
                A%LINE_PAT(NL)=0
             END DO
          END IF
       ENDIF

       IF(PRESENT(XVALUE_POS)) THEN
          A%XVALUE_POS=XVALUE_POS
       ELSE
          A%XVALUE_POS=1
       ENDIF
       IF(PRESENT(YVALUE_POS)) THEN
          A%YVALUE_POS=YVALUE_POS
       ELSE
          A%YVALUE_POS=1
       ENDIF
       IF(PRESENT(FVALUE_POS)) THEN
          A%FVALUE_POS=FVALUE_POS
       ELSE
          A%FVALUE_POS=2
       ENDIF
       
       SELECT CASE(A%MODE_2D)
       CASE(11:12)
          IF(PRESENT(BEV_XLEN)) THEN
             A%BEV_XLEN=GUCLIP(BEV_XLEN)*GFACTOR
          ELSE
             A%BEV_XLEN=15.0*GFACTOR
          END IF
          IF(PRESENT(BEV_YLEN)) THEN
             A%BEV_YLEN=GUCLIP(BEV_YLEN)*GFACTOR
          ELSE
             A%BEV_YLEN=30.0*GFACTOR
          END IF
          IF(PRESENT(BEV_ZLEN)) THEN
             A%BEV_ZLEN=GUCLIP(BEV_ZLEN)*GFACTOR
          ELSE
             A%BEV_ZLEN=15.0*GFACTOR
          END IF
          IF(PRESENT(BEV_PHI)) THEN
             A%BEV_PHI=GUCLIP(BEV_PHI)
          ELSE
             A%BEV_PHI=-60.0
          END IF
          IF(PRESENT(BEV_CHI)) THEN
             A%BEV_CHI=GUCLIP(BEV_CHI)
          ELSE
             A%BEV_CHI=65.0
          END IF
          IF(PRESENT(BEV_DISTANCE)) THEN
             A%BEV_DISTANCE=GUCLIP(BEV_DISTANCE)
          ELSE
             A%BEV_DISTANCE=100.0
          END IF
          
          IF(PRESENT(BEV_TYPE)) THEN
             A%BEV_TYPE=BEV_TYPE
          ELSE
             SELECT CASE(A%MODE_2D)
             CASE(11)
                IF(NXMAX.GT.100.OR.NYMAX.GT.100) THEN
                   A%BEV_TYPE=4
                ELSE
                   A%BEV_TYPE=7
                ENDIF
             CASE(12)
                IF(NXMAX.GT.100.OR.NYMAX.GT.100) THEN
                   A%BEV_TYPE=24
                ELSE
                   A%BEV_TYPE=27
                ENDIF
             END SELECT
          END IF
       END SELECT
    END SELECT

!     ----- TITLE SECTION -----

    IF(PRESENT(TITLE)) THEN
       KDL=TITLE(1:1)
       NCH=2
1      IF(TITLE(NCH:NCH).EQ.KDL.OR.NCH.EQ.LEN(TITLE)) GOTO 2
       KTITLE(NCH-1:NCH-1)=TITLE(NCH:NCH)
       NCH=NCH+1
       GOTO 1
2      CONTINUE
       A%TITLE=KTITLE
       A%TITLE_LEN=NCH-2
    ELSE
       A%TITLE_LEN=0
    ENDIF
       
    IF(PRESENT(TITLE_SIZE)) THEN
       A%TITLE_SIZE=GUCLIP(TITLE_SIZE)*GFACTOR
    ELSE
       A%TITLE_SIZE=0.6*GFACTOR
    END IF

    IF(PRESENT(TITLE_FONT)) THEN
       A%TITLE_FONT=TITLE_FONT
    ELSE
       A%TITLE_FONT=32
    END IF

    IF(PRESENT(TITLE_RGB)) THEN
       DO I=1,3
          A%TITLE_RGB(I)=GUCLIP(TITLE_RGB(I))
       END DO
    ELSE
       A%TITLE_RGB(1:3)=0.0
    END IF

    IF(PRESENT(TITLE_SEP)) THEN
       A%TITLE_SEP=GUCLIP(TITLE_SEP)*GFACTOR
    ELSE
       A%TITLE_SEP=0.3*GFACTOR
    END IF

    IF(PRESENT(TITLE_POS)) THEN
       A%TITLE_POS=TITLE_POS
    ELSE
       A%TITLE_POS=0
    END IF
      
    IF(PRESENT(XTITLE)) THEN
       KDL=XTITLE(1:1)
       NCH=2
3      IF(XTITLE(NCH:NCH).EQ.KDL.OR.NCH.EQ.LEN(XTITLE)) GOTO 4
       KTITLE(NCH-1:NCH-1)=XTITLE(NCH:NCH)
       NCH=NCH+1
       GOTO 3
4      CONTINUE
       A%XTITLE=KTITLE
       A%XTITLE_LEN=NCH-2
    ELSE
       A%XTITLE_LEN=0
    ENDIF
       
    IF(PRESENT(XTITLE_SIZE)) THEN
       A%XTITLE_SIZE=GUCLIP(XTITLE_SIZE)*GFACTOR
    ELSE
       A%XTITLE_SIZE=A%TITLE_SIZE
    END IF

    IF(PRESENT(XTITLE_FONT)) THEN
       A%XTITLE_FONT=XTITLE_FONT
    ELSE
       A%XTITLE_FONT=A%TITLE_FONT
    END IF

    IF(PRESENT(XTITLE_RGB)) THEN
       DO I=1,3
          A%XTITLE_RGB(I)=GUCLIP(XTITlE_RGB(I))
       END DO
    ELSE
       A%XTITLE_RGB(1:3)=A%TITlE_RGB(1:3)
    END IF

    IF(PRESENT(XTITLE_SEP)) THEN
       A%XTITLE_SEP=GUCLIP(XTITLE_SEP)*GFACTOR
    ELSE
       A%XTITLE_SEP=A%TITLE_SIZE*3.0+A%TITLE_SEP
    END IF

    IF(PRESENT(YTITLE)) THEN
       KDL=YTITLE(1:1)
       NCH=2
5      IF(YTITLE(NCH:NCH).EQ.KDL.OR.NCH.EQ.LEN(YTITLE)) GOTO 6
       KTITLE(NCH-1:NCH-1)=YTITLE(NCH:NCH)
       NCH=NCH+1
       GOTO 5
6      CONTINUE
       A%YTITLE=KTITLE
       A%YTITLE_LEN=NCH-2
    ELSE
       A%YTITLE_LEN=0
    ENDIF
       
    IF(PRESENT(YTITLE_SIZE)) THEN
       A%YTITLE_SIZE=GUCLIP(YTITLE_SIZE)*GFACTOR
    ELSE
       A%YTITLE_SIZE=A%TITLE_SIZE
    END IF

    IF(PRESENT(YTITLE_FONT)) THEN
       A%YTITLE_FONT=YTITLE_FONT
    ELSE
       A%YTITLE_FONT=A%TITLE_FONT
    END IF

    IF(PRESENT(YTITLE_RGB)) THEN
       DO I=1,3
          A%YTITLE_RGB(I)=GUCLIP(YTITlE_RGB(I))
       END DO
    ELSE
       A%YTITLE_RGB(1:3)=A%TITlE_RGB(1:3)
    END IF

    IF(PRESENT(YTITLE_SEP)) THEN
       A%YTITLE_SEP=GUCLIP(YTITLE_SEP)*GFACTOR
    ELSE
       A%YTITLE_SEP=A%YTITLE_SIZE*3.0+A%TITLE_SEP
    END IF

!     ----- FRAME SECTION -----

    IF(PRESENT(FRAME_RGB)) THEN
       DO I=1,3
          A%FRAME_RGB(I)=GUCLIP(FRAME_RGB(I))
       END DO
    ELSE
       A%FRAME_RGB(1:3)=0.0
    END IF

    IF(PRESENT(FRAME_THICKNESS)) THEN
       A%FRAME_THICKNESS=GUCLIP(FRAME_THICKNESS)*GFACTOR
    ELSE
       A%FRAME_THICKNESS=0.07*GFACTOR
    END IF

!     ----- SCALE SECTION -----

    IF(PRESENT(SCALE_RGB)) THEN
       DO I=1,3
          A%SCALE_RGB(I)=GUCLIP(SCALE_RGB(I))
       END DO
    ELSE
       A%SCALE_RGB(1:3)=A%FRAME_RGB(1:3)
    END IF

    IF(PRESENT(SCALE_ZERO_RGB)) THEN
       DO I=1,3
          A%SCALE_ZERO_RGB(I)=GUCLIP(SCALE_ZERO_RGB(I))
       END DO
    ELSE
       A%SCALE_ZERO_RGB(1:3)=A%FRAME_RGB(1:3)
    END IF

    IF(PRESENT(SCALE_THICKNESS)) THEN
       A%SCALE_THICKNESS=GUCLIP(SCALE_THICKNESS)*GFACTOR
    ELSE
       A%SCALE_THICKNESS=A%FRAME_THICKNESS
    END IF

    IF(PRESENT(SCALE_ZERO_THICKNESS)) THEN
       A%SCALE_ZERO_THICKNESS=GUCLIP(SCALE_ZERO_THICKNESS)*GFACTOR
    ELSE
       A%SCALE_ZERO_THICKNESS=A%FRAME_THICKNESS
    END IF

    IF(PRESENT(VALUE_FONT)) THEN
       A%VALUE_FONT=VALUE_FONT
    ELSE
       A%VALUE_FONT=A%TITLE_FONT
    END IF

    IF(PRESENT(VALUE_SIZE)) THEN
       A%VALUE_SIZE=GUCLIP(VALUE_SIZE)*GFACTOR
    ELSE
       A%VALUE_SIZE=A%TITLE_SIZE
    END IF

    IF(PRESENT(VALUE_RGB)) THEN
       DO I=1,3
          A%VALUE_RGB(I)=GUCLIP(VALUE_RGB(I))
       END DO
    ELSE
       A%VALUE_RGB(1:3)=A%FRAME_RGB(1:3)
    END IF

!     ----- XSCALE SECTION -----

    IF(PRESENT(XSCALE_SIZE)) THEN
       A%XSCALE_SIZE=GUCLIP(XSCALE_SIZE)*GFACTOR
    ELSE
       A%XSCALE_SIZE=0.2*GFACTOR
    END IF

    IF(PRESENT(XSCALE_TYPE)) THEN
       A%XSCALE_TYPE=XSCALE_TYPE
    ELSE
       SELECT CASE(A%MODE_2D)
       CASE(11:12)
          A%XSCALE_TYPE=1
       CASE DEFAULT
          A%XSCALE_TYPE=9
       END SELECT
    END IF

    IF(PRESENT(XSCALE_LTYPE)) THEN
       A%XSCALE_LTYPE=XSCALE_LTYPE
    ELSE
       A%XSCALE_LTYPE=9
    END IF

    IF(PRESENT(XGRID_LTYPE)) THEN
       A%XGRID_LTYPE=XGRID_LTYPE
    ELSE
       A%XGRID_LTYPE=1
    END IF

    IF(PRESENT(XSCALE_ZERO)) THEN
       A%XSCALE_ZERO=XSCALE_ZERO
    ELSE
       A%XSCALE_ZERO=1
    END IF

    IF(PRESENT(XSCALE_ZERO_PAT)) THEN
       A%XSCALE_ZERO_PAT=XSCALE_ZERO_PAT
    ELSE
       A%XSCALE_ZERO_PAT=0
    END IF

    IF(PRESENT(XVALUE_STEP)) THEN
       A%XVALUE_STEP=XVALUE_STEP
    ELSE
       A%XVALUE_STEP=2
    END IF

    IF(PRESENT(XVALUE_TYPE)) THEN
       A%XVALUE_TYPE=XVALUE_TYPE
    ELSE
       A%XVALUE_TYPE=NGULEN(2*GXSTEP)
    END IF

    IF(PRESENT(XVALUE_LTYPE)) THEN
       A%XVALUE_LTYPE=XVALUE_LTYPE
    ELSE
       A%XVALUE_LTYPE=1
    END IF

!     ----- YSCALE SECTION -----

    IF(PRESENT(YSCALE_SIZE)) THEN
       A%YSCALE_SIZE=GUCLIP(YSCALE_SIZE)*GFACTOR
    ELSE
       A%YSCALE_SIZE=0.2*GFACTOR
    END IF

    IF(PRESENT(YSCALE_TYPE)) THEN
       A%YSCALE_TYPE=YSCALE_TYPE
    ELSE
       SELECT CASE(A%MODE_2D)
       CASE(11:12)
          A%YSCALE_TYPE=1
       CASE DEFAULT
          A%YSCALE_TYPE=9
       END SELECT
    END IF

    IF(PRESENT(YSCALE_LTYPE)) THEN
       A%YSCALE_LTYPE=YSCALE_LTYPE
    ELSE
       A%YSCALE_LTYPE=9
    END IF

    IF(PRESENT(YGRID_LTYPE)) THEN
       A%YGRID_LTYPE=YGRID_LTYPE
    ELSE
       A%YGRID_LTYPE=1
    END IF

    IF(PRESENT(YSCALE_ZERO)) THEN
       A%YSCALE_ZERO=YSCALE_ZERO
    ELSE
       A%YSCALE_ZERO=1
    END IF

    IF(PRESENT(YSCALE_ZERO_PAT)) THEN
       A%YSCALE_ZERO_PAT=YSCALE_ZERO_PAT
    ELSE
       A%YSCALE_ZERO_PAT=0
    END IF

    IF(PRESENT(YVALUE_STEP)) THEN
       A%YVALUE_STEP=YVALUE_STEP
    ELSE
       A%YVALUE_STEP=2
    END IF

    IF(PRESENT(YVALUE_TYPE)) THEN
       A%YVALUE_TYPE=YVALUE_TYPE
    ELSE
       A%YVALUE_TYPE=NGULEN(2*GYSTEP)
    END IF

    IF(PRESENT(YVALUE_LTYPE)) THEN
       A%YVALUE_LTYPE=YVALUE_LTYPE
    ELSE
       A%YVALUE_LTYPE=1
    END IF

!     ----- FSCALE SECTION -----

    IF(PRESENT(FSCALE_SIZE)) THEN
       A%FSCALE_SIZE=GUCLIP(FSCALE_SIZE)*GFACTOR
    ELSE
       A%FSCALE_SIZE=0.2*GFACTOR
    END IF

    IF(PRESENT(FSCALE_TYPE)) THEN
       A%FSCALE_TYPE=FSCALE_TYPE
    ELSE
       SELECT CASE(A%MODE_2D)
       CASE(11:12)
          A%FSCALE_TYPE=1
       CASE DEFAULT
          A%FSCALE_TYPE=9
       END SELECT
    END IF

    IF(PRESENT(FSCALE_LTYPE)) THEN
       A%FSCALE_LTYPE=FSCALE_LTYPE
    ELSE
       A%FSCALE_LTYPE=9
    END IF

    IF(PRESENT(FGRID_LTYPE)) THEN
       A%FGRID_LTYPE=FGRID_LTYPE
    ELSE
       A%FGRID_LTYPE=1
    END IF

    IF(PRESENT(FSCALE_ZERO)) THEN
       A%FSCALE_ZERO=FSCALE_ZERO
    ELSE
       A%FSCALE_ZERO=1
    END IF

    IF(PRESENT(FSCALE_ZERO_PAT)) THEN
       A%FSCALE_ZERO_PAT=FSCALE_ZERO_PAT
    ELSE
       A%FSCALE_ZERO_PAT=0
    END IF

    IF(PRESENT(FVALUE_STEP)) THEN
       A%FVALUE_STEP=FVALUE_STEP
    ELSE
       A%FVALUE_STEP=2
    END IF

    IF(PRESENT(FVALUE_TYPE)) THEN
       A%FVALUE_TYPE=FVALUE_TYPE
    ELSE
       A%FVALUE_TYPE=NGULEN(2*GFSTEP)
    END IF

    IF(PRESENT(FVALUE_LTYPE)) THEN
       A%FVALUE_LTYPE=FVALUE_LTYPE
    ELSE
       A%FVALUE_LTYPE=1
    END IF

!   --- draw selection section ---

    IF(PRESENT(NOTITLE)) THEN
       A%NOTITLE=NOTITLE
    ELSE
       A%NOTITLE=0
    ENDIF

    IF(PRESENT(NOFRAME)) THEN
       A%NOFRAME=NOFRAME
    ELSE
       A%NOFRAME=0
    ENDIF

    IF(PRESENT(NOINFO)) THEN
       A%NOINFO=NOINFO
    ELSE
       A%NOINFO=0
    ENDIF

    RETURN
  END SUBROUTINE GRD_CONVERT
END MODULE grdconvert
