MODULE grfxyexec

  PRIVATE
  PUBLIC grfxy_exec

CONTAINS

! ***** DRAW 2D trajectory *****

  SUBROUTINE GRFXY_EXEC(GD,NDIM,NTM,NTMAX,NDMAX,A)

    USE grftype, ONLY: grf_attr_type
    USE grfutils, ONLY: grf_title, grf_frame2d, grf_frame3d ,grf_info,setrgba
    IMPLICIT NONE

    INTEGER,INTENT(IN):: NDIM,NTM,NDMAX,NTMAX(NDMAX)
    REAL(4),INTENT(IN):: GD(NDIM,NTM,NDMAX)
    REAL(4):: line_width,x1,x2,x3
    TYPE(grf_attr_type),INTENT(IN):: A

    REAL(4),ALLOCATABLE:: GX(:),GY(:)
    INTEGER(4):: ND,NT,NL,NTMAX_ALL

    CALL INQLNW(line_width)

    SELECT CASE(A%MODE_2D)
    CASE(21,22) ! 2D trajectory

       CALL GDEFIN(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX, &
                   A%XMIN,A%XMAX,A%YMIN,A%YMAX)
       CALL SETCLP(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX)

       NTMAX_ALL=NTMAX(1)
       DO ND=2,NDMAX
          NTMAX_ALL=MAX(NTMAX_ALL,NTMAX(ND))
       END DO
       ALLOCATE(GX(NTMAX_ALL),GY(NTMAX_ALL))
       DO ND=1,NDMAX
          NL=MOD(ND-1,A%NLMAX)+1
          CALL SETLNW(A%LINE_THICKNESS(NL))
          CALL SETMKS(0,A%LINE_MARK_SIZE(NL))
          SELECT CASE(A%MODE_2D)
          CASE(21)
             CALL SETRGBA(A%LINE_RGB(1:3,NL))
             DO NT=1,NTMAX(ND)
                GX(NT)=GD(1,NT,ND)
                GY(NT)=GD(2,NT,ND)
             END DO
             CALL GPLOTP(GX,GY,A%NXMIN,NTMAX(ND),A%NXSTEP, &
                         A%LINE_MARK(NL),A%LINE_MARK_STEP(NL),A%LINE_PAT(NL))
          CASE(22)
             CALL SETRGBA(A%LINE_RGB(1:3,NL))
             DO NT=1,NTMAX(ND)
                GX(NT)=GD(1,NT,ND)
                GY(NT)=GD(2,NT,ND)
             END DO
             CALL GPLOTPG(GX,GY,A%NXMIN,NTMAX(ND),A%NXSTEP, &
                          A%LINE_MARK(NL),A%LINE_MARK_STEP(NL), &
                          A%LINE_PAT(NL),1)
          END SELECT
       END DO
       CALL OFFCLP

       IF(A%NOFRAME.EQ.0) CALL GRF_FRAME2D(A)
       IF(A%NOTITLE.EQ.0) CALL GRF_TITLE(A)
       IF(A%NOINFO.EQ.0)  CALL GRF_INFO(A)

    END SELECT
    CALL SETRGB(0.0,0.0,0.0)
    CALL SETLNW(line_width)
    RETURN
  END SUBROUTINE GRFXY_EXEC

END MODULE grfxyexec
