MODULE grftype

  TYPE grf_attr_type
     REAL:: GPXMIN,GPXMAX,GPYMIN,GPYMAX ! graph size
     REAL:: XMIN,XMAX,XORG    ! X axis range, orgin of scale and value
     REAL:: YMIN,YMAX,YORG    ! Y axis range, orgin of scale and value
     REAL:: FMIN,FMAX,FORG    ! Z axis range, orgin of scale and value
     REAL:: XSPACE_FACTOR,YSPACE_FACTOR ! Axis range reduction
     INTEGER:: NXMIN,NYMIN   ! Top of plot data
     INTEGER:: NXSTEP,NYSTEP ! Interval of plot data
     INTEGER:: NLMAX         ! Number of line attribute
     INTEGER:: NPMAX         ! Number of paint attribute
     REAL,DIMENSION(:),ALLOCATABLE:: LINE_VALUE   ! Line value (NLMAX)
     REAL,DIMENSION(:,:),ALLOCATABLE :: LINE_RGB  ! Line color array (NLMAX)
     REAL,DIMENSION(:),ALLOCATABLE:: LINE_THICKNESS ! Line thickness (NLMAX)
     INTEGER,DIMENSION(:),ALLOCATABLE:: LINE_PAT       ! Line pattern (NLMAX)
     INTEGER,DIMENSION(:),ALLOCATABLE:: LINE_MARK      ! Mark type (NLMAX)
     INTEGER,DIMENSION(:),ALLOCATABLE:: LINE_MARK_STEP ! Mark interval (NLMAX)
     REAL,DIMENSION(:),ALLOCATABLE:: LINE_MARK_SIZE ! Mark size (NLMAX)
     REAL,DIMENSION(:),ALLOCATABLE:: PAINT_VALUE ! Paint value (NPMAX)
     REAL,DIMENSION(:,:),ALLOCATABLE :: PAINT_RGB! Paint color array (NPMAX)
     REAL:: ASPECT       ! Aspect ratio of the graph 
!                                  Default 0.75
!                                  Square for 1.0
!                                  0.0 for (YMAX-YMIN)/(XMAX-XMIN)
     REAL:: BEV_XLEN     ! Bird's eye view : length of x-axis [15 cm]
     REAL:: BEV_YLEN     ! Bird's eye view : length of y-axis [30 cm]
     REAL:: BEV_ZLEN     ! Bird's eye view : length of z-axis [15 cm]
     REAL:: BEV_PHI      ! Bird's eye view : angle in xy plane [-60 deg.]
     REAL:: BEV_CHI      ! Bird's eye view : angle from xy plane [65 deg.]
     REAL:: BEV_DISTANCE ! Bird's eye view : distance from origin [100 cm]

     INTEGER:: BEV_TYPE     ! Bird's eye view : type of distance from origin
     CHARACTER(LEN=80):: TITLE,XTITLE,YTITLE       ! Title contents
     INTEGER:: TITLE_LEN,XTITLE_LEN,YTITLE_LEN     ! Number of chars in title
     REAL:: TITLE_SIZE,XTITLE_SIZE,YTITLE_SIZE       ! Title font size
     INTEGER:: TITLE_FONT,XTITLE_FONT,YTITLE_FONT       ! Title font type
     REAL:: TITLE_RGB(3),XTITLE_RGB(3),YTITLE_RGB(3) ! Title font color
     REAL:: TITLE_SEP,XTITLE_SEP,YTITLE_SEP  ! Title separation from graph
     INTEGER:: TITLE_POS ! Title horizontal position
!                           0 : Flush left (default)
!                           1 : Centering
!                           2 : Flush right
!                           Always centering for X/Y/Z titles
     REAL:: FRAME_RGB(3)     ! Frame and scale color
     REAL:: FRAME_THICKNESS  ! Frame thickness
     INTEGER:: MODE_LS   ! Scale type
!                           0 : X:LINEAR  Y:LINEAR
!                           1 : X:LOG     Y:LINEAR
!                           2 : X:LINEAR  Y:LOG
!                           3 : X:LOG     Y:LOG
     REAL:: SCALE_RGB(3)          ! Scale color
     REAL:: SCALE_Zero_RGB(3)     ! Zero scale color
     REAL:: SCALE_THICKNESS       ! Scale thickness
     REAL:: SCALE_ZERO_THICKNESS  ! Zero Scale thickness
     REAL:: XSCALE_STEP,YSCALE_STEP,FSCALE_STEP ! Scale step size
     REAL:: XSCALE_SIZE,YSCALE_SIZE,FSCALE_SIZE ! Scale length or grid pat.
     INTEGER:: XSCALE_TYPE,YSCALE_TYPE,FSCALE_TYPE ! Scale type (0 for grid)
     INTEGER:: XSCALE_LTYPE,YSCALE_LTYPE,FSCALE_LTYPE ! Scale type for log plot
     INTEGER:: XSCALE_ZERO,YSCALE_ZERO,FSCALE_ZERO ! /=0 for origin grid
     INTEGER:: XSCALE_ZERO_PAT,YSCALE_ZERO_PAT,FSCALE_ZERO_PAT ! org grid pat
     REAL:: XGRID_STEP,YGRID_STEP,FGRID_STEP ! Grid step size (/=0 for grid)
     INTEGER:: XGRID_LTYPE,YGRID_LTYPE,FGRID_LTYPE ! Scale type for log plot
     INTEGER:: VALUE_FONT   ! Value font type
     REAL:: VALUE_SIZE   ! Value font size
     REAL:: VALUE_RGB(3) ! Value font color
     INTEGER:: XVALUE_STEP,YVALUE_STEP,FVALUE_STEP ! Value step size
     INTEGER:: XVALUE_POS,YVALUE_POS,FVALUE_POS    ! Value pos for 3D
     INTEGER:: XVALUE_TYPE,YVALUE_TYPE,FVALUE_TYPE ! Value type
     INTEGER:: XVALUE_LTYPE,YVALUE_LTYPE,FVALUE_LTYPE ! Value type for log plot
     INTEGER:: MODE_2D   ! Graph tupe
!                           0 : 1D plot (multi lines)  GRF1D
!                           1 : 2D contour lines       CONTG
!                           2 : 2D equi-value paint    CONTP
!                           3 : 2D contour lines and equi-value paint
!                           4 : 2D contour lines       CONTQ fixed line prop.
!                          11 : Bird's eye view contour lines    GRF2DC
!                          12 : Bird's eye view equi-value paint GRF2DD
!                          13 : Bird's eye view contour and paint
!                          14 : Bird's eye view contour on Z=0 plane
!                          15 : Bird's eye view contour on X=0 plane
!                          16 : Bird's eye view contour on Y=0 plane
!                          17 : Bird's eye view X mesh
!                          18 : Bird's eye view Y mesh
!                          19 : Bird's eye view X and Y mesh
!                          21 : X-Y plot
     INTEGER:: MODE_XY   ! Coordinate tupe
!                           0 : Rectangular coordinates
!                           1 : Circle
!                           2 : Half circle
     INTEGER:: MODE_PRD  ! Periodicity mode
!                           0 : no periodicity
!                           1 : periodic in X
!                           2 : periodic in Y
!                           3 : periodic in both X and Y
     INTEGER:: MODE_SPL     ! Spline tupe for 2D contour
!                           0 : No spline
!                           1 : Same point number for spline
!                           n : n times point number for spline
     INTEGER:: FRAME_TYPE   ! Frame type
!                           0 : Rectangular frame (Default)
!                           1 : XY axis only
     INTEGER:: NOTITLE                      ! /= 0 for no title
     INTEGER:: NOFRAME                      ! /= 0 for no frame
     INTEGER:: NOINFO                       ! /= 0 for no info (min,max,step)
     INTEGER:: NOXSCALE,NOYSCALE,NOFSCALE   ! /= 0 for no scale
     INTEGER:: NOXVALUE,NOYVALUE,NOFVALUE   ! /= 0 for no value
  END TYPE grf_attr_type

END MODULE grftype
