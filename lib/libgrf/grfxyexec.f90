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
    REAL,INTENT(IN):: GD(NDIM,NTM,NDMAX)
    INTEGER:: NL1,NL2
    REAL:: line_width,line_rgb(3),RNL,FACTOR1,FACTOR2
    TYPE(grf_attr_type),INTENT(IN):: A

    REAL,ALLOCATABLE:: GX(:),GY(:)
    INTEGER:: ND,NT,NL,NTMAX_ALL

    CALL INQLNW(line_width)

    SELECT CASE(A%MODE_2D)
    CASE(21,22,23) ! 2D trajectory

       CALL GDEFIN(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX, &
                   A%XMIN,A%XMAX,A%YMIN,A%YMAX)
       CALL SETCLP(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX)

       NTMAX_ALL=NTMAX(1)
       DO ND=2,NDMAX
          NTMAX_ALL=MAX(NTMAX_ALL,NTMAX(ND))
       END DO
       ALLOCATE(GX(NTMAX_ALL),GY(NTMAX_ALL))
       DO ND=1,NDMAX
          SELECT CASE(A%MODE_2D)
          CASE(21,22)
             NL=MOD(ND-1,A%NLMAX)+1
             CALL SETLNW(A%LINE_THICKNESS(NL))
             CALL SETMKS(0,A%LINE_MARK_SIZE(NL))
          END SELECT
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
          CASE(23)
             DO NT=A%NXMIN,NTMAX(ND)-1,A%NXSTEP
                GX(1)=GD(1,NT  ,ND)
                GY(1)=GD(2,NT  ,ND)
                GX(2)=GD(1,NT+1,ND)
                GY(2)=GD(2,NT+1,ND)
                RNL=(NTMAX(ND)-REAL(NT))/(REAL(NTMAX(ND))-A%NXMIN)
                NL1=MAX(1,INT((A%NLMAX-1)*RNL)+1)
                NL2=MIN(A%NLMAX,INT((A%NLMAX-1)*RNL)+2)
                FACTOR1=NL1-RNL*(A%NLMAX-1)
                FACTOR2=1.0-FACTOR1
!                WRITE(6,'(A,4I8,3ES12.4)') 'NT:',NT,NL1,NL2, &
!                     A%NLMAX,RNL,FACTOR1,FACTOR2
                line_rgb(1)=factor1*A%LINE_RGB(1,NL1)+factor2*A%LINE_RGB(1,NL2)
                line_rgb(2)=factor1*A%LINE_RGB(2,NL1)+factor2*A%LINE_RGB(2,NL2)
                line_rgb(3)=factor1*A%LINE_RGB(3,NL1)+factor2*A%LINE_RGB(3,NL2)
                CALL SETRGBA(line_rgb)
                IF(MOD(NT,A%LINE_MARK_STEP(1)).EQ.0) THEN
                   CALL GPLOTP(GX,GY,1,2,1,A%LINE_MARK(1),1,A%LINE_PAT(1))
                ELSE
                   CALL GPLOTP(GX,GY,1,2,1,0,1,A%LINE_PAT(1))
                END IF
             END DO
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
