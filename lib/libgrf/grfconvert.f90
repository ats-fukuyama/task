MODULE grfconvert

  PRIVATE
  PUBLIC grf_convert

CONTAINS

! ***** Convert input parameters to A: REAL(4) *****

  SUBROUTINE GRF_CONVERT(NGP,GX,GY,GF,NXM,NXMAX,NYMAX, &
                         GPXMIN,GPXMAX,GPYMIN,GPYMAX, &
                         XMIN,XMAX,XORG,YMIN,YMAX,YORG,FMIN,FMAX,FORG, &
                         XSPACE_FACTOR,YSPACE_FACTOR,FSPACE_FACTOR, &
                         NXMIN,NYMIN,NXSTEP,NYSTEP,NLMAX,NPMAX, &
                         LINE_VALUE,LINE_RGB,LINE_THICKNESS,LINE_PAT, &
                         LINE_MARK,LINE_MARK_STEP,LINE_MARK_SIZE, &
                         PAINT_VALUE,PAINT_RGB, &
                         BEV_XORG,BEV_YORG,BEV_ZORG, &
                         BEV_PHI,BEV_CHI,BEV_DISTANCE, &
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
                         XVALUE_TYPE,YVALUE_TYPE,FVALUE_TYPE, &
                         XVALUE_LTYPE,YVALUE_LTYPE,FVALUE_LTYPE, &
                         MODE_2D,MODE_XY,MODE_PRD,MODE_LS,FRAME_TYPE, &
                         NOTITLE,NOFRAME, &
                         NOXSCALE,NOYSCALE,NOFSCALE, &
                         NOXVALUE,NOYVALUE,NOFVALUE, &
                         LINE_RGB_SUB,PAINT_RGB_SUB,A)

    USE grftype, ONLY: grf_attr_type
    USE grfutils, ONLY: grfut1,grfut2,grfut3,grfut4
    IMPLICIT NONE

    INTEGER,INTENT(IN):: NGP
    REAL(4),INTENT(IN):: GX(NXMAX),GY(NYMAX),GF(NXM,NYMAX)
    INTEGER,INTENT(IN):: NXM,NXMAX,NYMAX

    REAL(4),INTENT(IN),OPTIONAL:: GPXMIN,GPXMAX,GPYMIN,GPYMAX
    REAL(4),INTENT(IN),OPTIONAL:: XMIN,XMAX,XORG
    REAL(4),INTENT(IN),OPTIONAL:: YMIN,YMAX,YORG
    REAL(4),INTENT(IN),OPTIONAL:: FMIN,FMAX,FORG
    REAL(4),INTENT(IN),OPTIONAL:: XSPACE_FACTOR,YSPACE_FACTOR,FSPACE_FACTOR
    INTEGER,INTENT(IN),OPTIONAL:: NXMIN,NYMIN,NXSTEP,NYSTEP

    INTEGER,INTENT(IN),OPTIONAL:: NLMAX,NPMAX
    REAL(4),INTENT(IN),OPTIONAL:: LINE_VALUE(*),LINE_RGB(3,*)
    REAL(4),INTENT(IN),OPTIONAL:: LINE_THICKNESS(*)
    INTEGER,INTENT(IN),OPTIONAL:: LINE_PAT(*)
    INTEGER,INTENT(IN),OPTIONAL:: LINE_MARK(*),LINE_MARK_STEP(*)
    REAL(4),INTENT(IN),OPTIONAL:: LINE_MARK_SIZE(*)
    REAL(4),INTENT(IN),OPTIONAL:: PAINT_VALUE(*),PAINT_RGB(3,*)
    REAL(4),INTENT(IN),OPTIONAL:: BEV_XORG,BEV_YORG,BEV_ZORG
    REAL(4),INTENT(IN),OPTIONAL:: BEV_PHI,BEV_CHI,BEV_DISTANCE

    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: TITLE,XTITLE,YTITLE
    INTEGER,INTENT(IN),OPTIONAL:: TITLE_LEN,XTITLE_LEN,YTITLE_LEN
    REAL(4),INTENT(IN),OPTIONAL:: TITLE_SIZE,XTITLE_SIZE,YTITLE_SIZE
    INTEGER,INTENT(IN),OPTIONAL:: TITLE_FONT,XTITLE_FONT,YTITLE_FONT
    REAL(4),INTENT(IN),OPTIONAL:: TITLE_RGB(3),XTITLE_RGB(3),YTITLE_RGB(3)
    REAL(4),INTENT(IN),OPTIONAL:: TITLE_SEP,XTITLE_SEP,YTITLE_SEP
    INTEGER,INTENT(IN),OPTIONAL:: TITLE_POS

    REAL(4),INTENT(IN),OPTIONAL:: FRAME_RGB(3),FRAME_THICKNESS

    REAL(4),INTENT(IN),OPTIONAL:: SCALE_RGB(3),SCALE_Zero_RGB(3)
    REAL(4),INTENT(IN),OPTIONAL:: SCALE_THICKNESS,SCALE_ZERO_THICKNESS
    REAL(4),INTENT(IN),OPTIONAL:: XSCALE_STEP,YSCALE_STEP,FSCALE_STEP
    REAL(4),INTENT(IN),OPTIONAL:: XSCALE_SIZE,YSCALE_SIZE,FSCALE_SIZE
    INTEGER,INTENT(IN),OPTIONAL:: XSCALE_TYPE,YSCALE_TYPE,FSCALE_TYPE
    INTEGER,INTENT(IN),OPTIONAL:: XSCALE_LTYPE,YSCALE_LTYPE,FSCALE_LTYPE
    INTEGER,INTENT(IN),OPTIONAL:: XSCALE_ZERO,YSCALE_ZERO,FSCALE_ZERO
    INTEGER,INTENT(IN),OPTIONAL:: XSCALE_ZERO_PAT,YSCALE_ZERO_PAT, &
                                  FSCALE_ZERO_PAT

    INTEGER,INTENT(IN),OPTIONAL:: VALUE_FONT
    REAL(4),INTENT(IN),OPTIONAL:: VALUE_SIZE,VALUE_RGB(3)
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_STEP,YVALUE_STEP,FVALUE_STEP
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_TYPE,YVALUE_TYPE,FVALUE_TYPE
    INTEGER,INTENT(IN),OPTIONAL:: XVALUE_LTYPE,YVALUE_LTYPE,FVALUE_LTYPE

    INTEGER,INTENT(IN),OPTIONAL:: MODE_2D,MODE_XY,MODE_PRD,MODE_LS,FRAME_TYPE
    INTEGER,INTENT(IN),OPTIONAL:: NOTITLE,NOFRAME
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
    END INTERFACE
    OPTIONAL LINE_RGB_SUB,PAINT_RGB_SUB

    TYPE(grf_attr_type):: A

    INTEGER,PARAMETER:: NLM=5
    REAL(4):: GP(4),GL,GFACTOR
    REAL(4):: GXMINL,GXMAXL,GYMINL,GYMAXL,GFMINL,GFMAXL
    REAL(4):: GXMIN,GXMAX,GXSTEP,GXORG
    REAL(4):: GYMIN,GYMAX,GYSTEP,GYORG
    REAL(4):: GFMIN,GFMAX,GFSTEP,GFORG
    REAL(4):: GRGB(3,NLM)
    INTEGER:: NL,NCH,NGULEN
    INTEGER:: IPAT(NLM)
    CHARACTER(LEN=80):: KTITLE
    CHARACTER(LEN=1):: KDL

    DATA IPAT/0,2,3,4,6/
    DATA GRGB/0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,0.0,1.0, &
              0.0,1.0,0.0, 1.0,0.0,1.0/

    IF(NGP.LT.0.OR.NGP.GT.29) THEN
       WRITE(6,*) 'XX GRF_CONVERT: INVALID NGP'
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
    IF(PRESENT(MODE_LS)) THEN
       A%MODE_LS=MODE_LS
    ELSE
       A%MODE_LS=0
    ENDIF

    IF(PRESENT(MODE_2D)) THEN
       A%MODE_2D=MODE_2D
    ELSE
       A%MODE_2D=0
    ENDIF

!     ----- define graph position -----

    CALL GRFUT4(NGP,GP) ! assign predefined graph position to GP
    IF(PRESENT(GPXMIN)) THEN
       A%GPXMIN=GPXMIN
    ELSE
       A%GPXMIN=GP(1)
    ENDIF
    IF(PRESENT(GPXMAX)) THEN
       A%GPXMAX=GPXMAX
    ELSE
       A%GPXMAX=GP(2)
    ENDIF
    IF(PRESENT(GPYMIN)) THEN
       A%GPYMIN=GPYMIN
    ELSE
       A%GPYMIN=GP(3)
    ENDIF
    IF(PRESENT(GPYMAX)) THEN
       A%GPYMAX=GPYMAX
    ELSE
       A%GPYMAX=GP(4)
    ENDIF

!     ----- define graph size scaling factor -----

    GL=MIN(ABS(A%GPXMAX-A%GPXMIN),ABS(A%GPYMAX-A%GPYMIN))
    GFACTOR=GL/15.0

!     ----- SCALE ADJUSTMENT SECTION -----

    IF(PRESENT(XSPACE_FACTOR)) THEN
       A%XSPACE_FACTOR=XSPACE_FACTOR
    ELSE
       A%XSPACE_FACTOR=0.0
    END IF
    IF(PRESENT(YSPACE_FACTOR)) THEN
       A%YSPACE_FACTOR=YSPACE_FACTOR
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

    IF(PRESENT(XMIN)) GXMINL=XMIN
    IF(PRESENT(XMAX)) GXMAXL=XMAX
    IF(PRESENT(YMIN)) GYMINL=YMIN
    IF(PRESENT(YMAX)) GYMAXL=YMAX
    IF(PRESENT(FMIN)) GFMINL=FMIN
    IF(PRESENT(FMAX)) GFMAXL=FMAX

    CALL GRFUT3(GXMINL,GXMAXL,GXMIN,GXMAX,GXSTEP,GXORG)
    CALL GRFUT3(GYMINL,GYMAXL,GYMIN,GYMAX,GYSTEP,GYORG)
    CALL GRFUT3(GFMINL,GFMAXL,GFMIN,GFMAX,GFSTEP,GFORG)

    IF(PRESENT(XMIN)) THEN
       A%XMIN=XMIN
    ELSE
       A%XMIN=GXMIN
    ENDIF
    IF(PRESENT(XMAX)) THEN
       A%XMAX=XMAX
    ELSE
       A%XMAX=GXMAX
    ENDIF
    IF(PRESENT(XSCALE_STEP)) THEN
       A%XSCALE_STEP=XSCALE_STEP
    ELSE
       A%XSCALE_STEP=GXSTEP
    ENDIF
    IF(PRESENT(XORG)) THEN
       A%XORG=XORG
    ELSE
       A%XORG=GXORG
    ENDIF

    IF(PRESENT(YMIN)) THEN
       A%YMIN=YMIN
    ELSE
       A%YMIN=GYMIN
    ENDIF
    IF(PRESENT(YMAX)) THEN
       A%YMAX=YMAX
    ELSE
       A%YMAX=GYMAXL
    ENDIF
    IF(PRESENT(YSCALE_STEP)) THEN
       A%YSCALE_STEP=YSCALE_STEP
    ELSE
       A%YSCALE_STEP=GYSTEP
    ENDIF
    IF(PRESENT(YORG)) THEN
       A%YORG=YORG
    ELSE
       A%YORG=GYORG
    ENDIF

    IF(PRESENT(FMIN)) THEN
       A%FMIN=FMIN
    ELSE
       A%FMIN=GFMIN
    ENDIF
    IF(PRESENT(FMAX)) THEN
       A%FMAX=FMAX
    ELSE
       A%FMAX=GFMAXL
    ENDIF
    IF(PRESENT(FSCALE_STEP)) THEN
       A%FSCALE_STEP=FSCALE_STEP
    ELSE
       A%FSCALE_STEP=GFSTEP
    ENDIF
    IF(PRESENT(FORG)) THEN
       A%FORG=FORG
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
    GL=A%FMAX-A%FMIN
    A%FMAX=A%FMAX+0.5*GL*A%FSPACE_FACTOR
    A%FMIN=A%FMIN-0.5*GL*A%FSPACE_FACTOR

!     ----- LINE SECTION -----

    IF(A%MODE_2D.EQ.0) THEN  ! 1D data
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
          A%LINE_RGB(1:3,1:A%NLMAX)=LINE_RGB(1:3,1:A%NLMAX)
       ELSE
          DO NL=1,A%NLMAX
             A%LINE_RGB(1:3,NL)=GRGB(1:3,MOD(NL-1,NLM)+1)
          ENDDO
       END IF
       IF(PRESENT(LINE_THICKNESS)) THEN
          A%LINE_THICKNESS(1:A%NLMAX)=LINE_THICKNESS(1:A%NLMAX)*GFACTOR
       ELSE
          A%LINE_THICKNESS(1:A%NLMAX)=0.07*GFACTOR
       ENDIF
       IF(PRESENT(LINE_PAT)) THEN
          A%LINE_PAT(1:A%NLMAX)=LINE_PAT(1:A%NLMAX)
       ELSE
          DO NL=1,A%NLMAX
             A%LINE_PAT(NL)=IPAT(MOD(NL-1,NLM)+1)
          END DO
       ENDIF
       IF(PRESENT(LINE_MARK)) THEN
          A%LINE_MARK(1:A%NLMAX)=LINE_MARK(1:A%NLMAX)
       ELSE
          A%LINE_MARK(1:A%NLMAX)=0
       ENDIF
       IF(PRESENT(LINE_MARK_STEP)) THEN
          A%LINE_MARK_STEP(1:A%NLMAX)=LINE_MARK_STEP(1:A%NLMAX)
       ELSE
          A%LINE_MARK_STEP(1:A%NLMAX)=0
       ENDIF
       IF(PRESENT(LINE_MARK_SIZE)) THEN
          A%LINE_MARK_SIZE(1:A%NLMAX)=LINE_MARK_SIZE(1:A%NLMAX)*GFACTOR
       ELSE
          A%LINE_MARK_SIZE(1:A%NLMAX)=0.3*GFACTOR
       ENDIF

    END IF

!     ----- PAINT SECTION -----

    IF(PRESENT(NPMAX)) THEN
       A%NPMAX=NPMAX
    ELSE
       A%NPMAX=0
    END IF
    IF(A%NPMAX.GT.0) THEN
       ALLOCATE(A%PAINT_VALUE(A%NPMAX))
       ALLOCATE(A%PAINT_RGB(3,A%NPMAX))
       IF(PRESENT(PAINT_VALUE)) THEN
          A%PAINT_VALUE(1:A%NLMAX)=PAINT_VALUE(1:A%NLMAX)
       END IF
       IF(PRESENT(PAINT_RGB)) THEN
          A%PAINT_RGB(1:3,1:A%NLMAX)=PAINT_RGB(1:3,1:A%NLMAX)
       ELSE
          A%PAINT_RGB(1:3,1:A%NLMAX)=GRGB(1:3,1:A%NLMAX)
       END IF
    END IF

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
       A%TITLE_SIZE=TITLE_SIZE*GFACTOR
    ELSE
       A%TITLE_SIZE=0.6*GFACTOR
    END IF

    IF(PRESENT(TITLE_FONT)) THEN
       A%TITLE_FONT=TITLE_FONT
    ELSE
       A%TITLE_FONT=32
    END IF

    IF(PRESENT(TITLE_RGB)) THEN
       A%TITLE_RGB(1:3)=TITLE_RGB(1:3)
    ELSE
       A%TITLE_RGB(1:3)=0.0
    END IF

    IF(PRESENT(TITLE_SEP)) THEN
       A%TITLE_SEP=TITLE_SEP*GFACTOR
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
       A%XTITLE_SIZE=XTITLE_SIZE*GFACTOR
    ELSE
       A%XTITLE_SIZE=A%TITLE_SIZE
    END IF

    IF(PRESENT(XTITLE_FONT)) THEN
       A%XTITLE_FONT=XTITLE_FONT
    ELSE
       A%XTITLE_FONT=A%TITLE_FONT
    END IF

    IF(PRESENT(XTITLE_RGB)) THEN
       A%XTITLE_RGB(1:3)=XTITlE_RGB(1:3)
    ELSE
       A%XTITLE_RGB(1:3)=A%TITlE_RGB(1:3)
    END IF

    IF(PRESENT(XTITLE_SEP)) THEN
       A%XTITLE_SEP=XTITLE_SEP*GFACTOR
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
       A%YTITLE_SIZE=YTITLE_SIZE*GFACTOR
    ELSE
       A%YTITLE_SIZE=A%TITLE_SIZE
    END IF

    IF(PRESENT(YTITLE_FONT)) THEN
       A%YTITLE_FONT=YTITLE_FONT
    ELSE
       A%YTITLE_FONT=A%TITLE_FONT
    END IF

    IF(PRESENT(YTITLE_RGB)) THEN
       A%YTITLE_RGB(1:3)=YTITlE_RGB(1:3)
    ELSE
       A%YTITLE_RGB(1:3)=A%TITlE_RGB(1:3)
    END IF

    IF(PRESENT(YTITLE_SEP)) THEN
       A%YTITLE_SEP=YTITLE_SEP*GFACTOR
    ELSE
       A%YTITLE_SEP=A%YTITLE_SIZE*3.0+A%TITLE_SEP
    END IF

!     ----- FRAME SECTION -----

    IF(PRESENT(FRAME_RGB)) THEN
       A%FRAME_RGB(1:3)=FRAME_RGB(1:3)
    ELSE
       A%FRAME_RGB(1:3)=0.0
    END IF

    IF(PRESENT(FRAME_THICKNESS)) THEN
       A%FRAME_THICKNESS=FRAME_THICKNESS*GFACTOR
    ELSE
       A%FRAME_THICKNESS=0.07*GFACTOR
    END IF

!     ----- SCALE SECTION -----

    IF(PRESENT(SCALE_RGB)) THEN
       A%SCALE_RGB(1:3)=SCALE_RGB(1:3)
    ELSE
       A%SCALE_RGB(1:3)=A%FRAME_RGB(1:3)
    END IF

    IF(PRESENT(SCALE_ZERO_RGB)) THEN
       A%SCALE_ZERO_RGB(1:3)=SCALE_ZERO_RGB(1:3)
    ELSE
       A%SCALE_ZERO_RGB(1:3)=A%FRAME_RGB(1:3)
    END IF

    IF(PRESENT(SCALE_THICKNESS)) THEN
       A%SCALE_THICKNESS=SCALE_THICKNESS*GFACTOR
    ELSE
       A%SCALE_THICKNESS=A%FRAME_THICKNESS
    END IF

    IF(PRESENT(SCALE_ZERO_THICKNESS)) THEN
       A%SCALE_ZERO_THICKNESS=SCALE_ZERO_THICKNESS*GFACTOR
    ELSE
       A%SCALE_ZERO_THICKNESS=A%FRAME_THICKNESS
    END IF

    IF(PRESENT(VALUE_FONT)) THEN
       A%VALUE_FONT=VALUE_FONT
    ELSE
       A%VALUE_FONT=A%TITLE_FONT
    END IF

    IF(PRESENT(VALUE_SIZE)) THEN
       A%VALUE_SIZE=VALUE_SIZE*GFACTOR
    ELSE
       A%VALUE_SIZE=A%TITLE_SIZE
    END IF

    IF(PRESENT(VALUE_RGB)) THEN
       A%VALUE_RGB(1:3)=VALUE_RGB(1:3)
    ELSE
       A%VALUE_RGB(1:3)=A%FRAME_RGB(1:3)
    END IF

!     ----- XSCALE SECTION -----

    IF(PRESENT(XSCALE_SIZE)) THEN
       A%XSCALE_SIZE=XSCALE_SIZE*GFACTOR
    ELSE
       A%XSCALE_SIZE=0.2*GFACTOR
    END IF

    IF(PRESENT(XSCALE_TYPE)) THEN
       A%XSCALE_TYPE=XSCALE_TYPE
    ELSE
       A%XSCALE_TYPE=9
    END IF

    IF(PRESENT(XSCALE_LTYPE)) THEN
       A%XSCALE_LTYPE=XSCALE_LTYPE
    ELSE
       A%XSCALE_LTYPE=9
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
       A%YSCALE_SIZE=YSCALE_SIZE*GFACTOR
    ELSE
       A%YSCALE_SIZE=0.2*GFACTOR
    END IF

    IF(PRESENT(YSCALE_TYPE)) THEN
       A%YSCALE_TYPE=YSCALE_TYPE
    ELSE
       A%YSCALE_TYPE=9
    END IF

    IF(PRESENT(YSCALE_LTYPE)) THEN
       A%YSCALE_LTYPE=YSCALE_LTYPE
    ELSE
       A%YSCALE_LTYPE=9
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
       A%FSCALE_SIZE=FSCALE_SIZE*GFACTOR
    ELSE
       A%FSCALE_SIZE=0.2*GFACTOR
    END IF

    IF(PRESENT(FSCALE_TYPE)) THEN
       A%FSCALE_TYPE=FSCALE_TYPE
    ELSE
       A%FSCALE_TYPE=9
    END IF

    IF(PRESENT(FSCALE_LTYPE)) THEN
       A%FSCALE_LTYPE=FSCALE_LTYPE
    ELSE
       A%FSCALE_LTYPE=9
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

    RETURN
  END SUBROUTINE GRF_CONVERT
END MODULE grfconvert
