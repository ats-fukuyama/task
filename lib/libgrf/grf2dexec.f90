MODULE grf2dexec

  PRIVATE
  PUBLIC grf2d_exec

CONTAINS

! ***** DRAW 2D GRAPH *****

  SUBROUTINE GRF2D_EXEC(GX,GY,GF,NXM,NXMAX,NYMAX,A, &
                        LINE_RGB_SUB,PAINT_RGB_SUB)

    USE grftype, ONLY: grf_attr_type
    USE grfutils, ONLY: grf_title, grf_frame2d, grf_frame3d ,grf_info,setrgba
    IMPLICIT NONE

    REAL(4),INTENT(IN):: GX(NXMAX),GY(NYMAX),GF(NXM,NYMAX)
    INTEGER,INTENT(IN):: NXM,NXMAX,NYMAX
    REAL(4):: line_width
    TYPE(grf_attr_type),INTENT(IN):: A
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

    INTEGER(4),DIMENSION(:,:),ALLOCATABLE:: KA
    REAL(4):: FSTEP
    INTEGER(4):: NL
    INTERFACE
       SUBROUTINE R2W2B(VALUE,RGB)
         REAL(4),INTENT(IN):: VALUE
         REAL(4),DIMENSION(3),INTENT(OUT)::RGB
       END SUBROUTINE R2W2B
       SUBROUTINE W2G2B(VALUE,RGB)
         REAL(4),INTENT(IN):: VALUE
         REAL(4),DIMENSION(3),INTENT(OUT):: RGB
       END SUBROUTINE W2G2B
       SUBROUTINE R2Y2W(VALUE,RGB)
         REAL(4),INTENT(IN):: VALUE
         REAL(4),DIMENSION(3),INTENT(OUT):: RGB
       END SUBROUTINE R2Y2W
    END INTERFACE

    CALL INQLNW(line_width)
    SELECT CASE(A%MODE_2D)
    CASE(1) ! 2D contour line

       CALL GDEFIN(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX, &
                   A%XMIN,A%XMAX,A%YMIN,A%YMAX)
       CALL SETCLP(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX)

       ALLOCATE(KA(2,NXMAX*NYMAX))
       SELECT CASE(A%MODE_XY)
       CASE(0)
          CALL CONTG2(GF,GX,GY,NXM,NXMAX,NYMAX, &
                      A%LINE_VALUE,A%LINE_RGB,A%LINE_PAT,A%LINE_THICKNESS, &
                      A%NLMAX,A%MODE_SPL,A%MODE_PRD,KA)
       CASE(1)
          CALL CONTG3(GF,GX,NXM,NXMAX,NYMAX, &
                      A%LINE_VALUE,A%LINE_RGB,A%LINE_PAT,A%LINE_THICKNESS, &
                      A%NLMAX,A%MODE_SPL,KA)
       CASE(2)
          CALL CONTG4(GF,GX,GY,NXM,NXMAX,NYMAX, &
                      A%LINE_VALUE,A%LINE_RGB,A%LINE_PAT,A%LINE_THICKNESS, &
                      A%NLMAX,A%MODE_SPL,KA)
       END SELECT

       DEALLOCATE(KA)

       IF(A%NOFRAME.EQ.0) CALL GRF_FRAME2D(A)
       IF(A%NOTITLE.EQ.0) CALL GRF_TITLE(A)
       IF(A%NOINFO.EQ.0)  CALL GRF_INFO(A)

    CASE(2) ! 2D contour paint
       CALL GDEFIN(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX, &
                   A%XMIN,A%XMAX,A%YMIN,A%YMAX)
       CALL SETCLP(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX)

       CALL SETLNW(A%LINE_THICKNESS(1))
       
       SELECT CASE(A%MODE_XY)
       CASE(0)
          CALL CONTF2(GF,GX,GY,NXM,NXMAX,NYMAX, &
                      A%PAINT_VALUE,A%PAINT_RGB,A%NPMAX-1,A%MODE_PRD)
       CASE(1)
          CALL CONTF3(GF,GX,NXM,NXMAX,NYMAX, &
                      A%PAINT_VALUE,A%PAINT_RGB,A%NPMAX-1)
       CASE(2)
          CALL CONTF4(GF,GX,GY,NXM,NXMAX,NYMAX, &
                      A%PAINT_VALUE,A%PAINT_RGB,A%NPMAX-1)
       END SELECT
       IF(A%NOFRAME.EQ.0) CALL GRF_FRAME2D(A)
       IF(A%NOTITLE.EQ.0) CALL GRF_TITLE(A)
       IF(A%NOINFO.EQ.0)  CALL GRF_INFO(A)

    CASE(3) ! 2D contour paint
       CALL GDEFIN(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX, &
                   A%XMIN,A%XMAX,A%YMIN,A%YMAX)
       CALL SETCLP(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX)

       CALL SETLNW(A%LINE_THICKNESS(1))
       CALL CONTF2(GF,GX,GY,NXM,NXMAX,NYMAX, &
                   A%PAINT_VALUE,A%PAINT_RGB,A%NPMAX-1,A%MODE_PRD)

       ALLOCATE(KA(2,NXMAX*NYMAX))
       CALL CONTG2(GF,GX,GY,NXM,NXMAX,NYMAX, &
                   A%LINE_VALUE,A%LINE_RGB,A%LINE_PAT,A%LINE_THICKNESS, &
                   A%NLMAX,A%MODE_SPL,A%MODE_PRD,KA)
       DEALLOCATE(KA)

       IF(A%NOFRAME.EQ.0) CALL GRF_FRAME2D(A)
       IF(A%NOTITLE.EQ.0) CALL GRF_TITLE(A)
       IF(A%NOINFO.EQ.0)  CALL GRF_INFO(A)

    CASE(4) ! 2D contour line

       CALL GDEFIN(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX, &
                   A%XMIN,A%XMAX,A%YMIN,A%YMAX)
       CALL SETCLP(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX)

       ALLOCATE(KA(8,NXMAX*NYMAX))
       FSTEP=A%FMAX-A%FMIN

       DO NL=1,A%NLMAX
          CALL SETRGB(A%LINE_RGB(1,NL),A%LINE_RGB(2,NL),A%LINE_RGB(3,NL))
          CALL SETLNW(A%LINE_THICKNESS(NL))
          SELECT CASE(A%MODE_XY)
          CASE(0)
             CALL CONTQ2(GF,GX,GY,NXM,NXMAX,NYMAX, &
                         A%LINE_VALUE(NL),FSTEP,1,A%MODE_PRD,A%LINE_PAT(NL),KA)
          CASE(1)
             CALL CONTQ3(GF,GX,NXM,NXMAX,NYMAX, &
                         A%LINE_VALUE(NL),FSTEP,1,A%LINE_PAT(NL),KA)
          CASE(2)
             CALL CONTQ4(GF,GX,GY,NXM,NXMAX,NYMAX, &
                         A%LINE_VALUE(NL),FSTEP,1,A%MODE_PRD,A%LINE_PAT(NL),KA)
          END SELECT
       END DO

       DEALLOCATE(KA)

       IF(A%NOFRAME.EQ.0) CALL GRF_FRAME2D(A)
       IF(A%NOTITLE.EQ.0) CALL GRF_TITLE(A)
       IF(A%NOINFO.EQ.0)  CALL GRF_INFO(A)

    CASE(11:12) ! 3D Bird's eye view
       CALL GDEFIN3D(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX, &
                     A%BEV_XLEN,A%BEV_YLEN,A%BEV_ZLEN)
       CALL GVIEW3D(A%BEV_PHI,A%BEV_CHI,A%BEV_DISTANCE,1.0,1, &
                    0.5*(A%XMIN+A%XMAX),0.5*(A%YMIN+A%YMAX), &
                    0.5*(A%FMIN+A%FMAX))
       CALL GDATA3D1(GF,NXM,NXMAX,NYMAX, &
                     A%XMIN,A%XMAX,A%YMIN,A%YMAX,A%FMIN,A%FMAX)

       SELECT CASE(A%MODE_2D)
       CASE(11)
          IF(PRESENT(PAINT_RGB_SUB)) THEN
             CALL CPLOT3D1(A%BEV_TYPE,PAINT_RGB_SUB)
          ELSE
             IF(A%FMIN*A%FMAX.LT.0.0) THEN
                CALL CPLOT3D1(A%BEV_TYPE,R2W2B)
             ELSEIF(A%FMIN.LT.0.0) THEN
                CALL CPLOT3D1(A%BEV_TYPE,W2G2B)
             ELSE
                CALL CPLOT3D1(A%BEV_TYPE,R2Y2W)
             ENDIF
          END IF
       CASE(12)
          IF(PRESENT(PAINT_RGB_SUB)) THEN
             CALL CPLOT3D1(MOD(A%BEV_TYPE,10),PAINT_RGB_SUB)
          ELSE
             IF(A%FMIN*A%FMAX.LT.0.0) THEN
                CALL CPLOT3D1(MOD(A%BEV_TYPE,10),R2W2B)
             ELSEIF(A%FMIN.LT.0.0) THEN
                CALL CPLOT3D1(MOD(A%BEV_TYPE,10),W2G2B)
             ELSE
                CALL CPLOT3D1(MOD(A%BEV_TYPE,10),R2Y2W)
             ENDIF
          END IF
          CALL SETLNW(A%LINE_THICKNESS(1))
          ALLOCATE(KA(8,NXMAX*NYMAX))
          IF(PRESENT(LINE_RGB_SUB)) THEN
             CALL CONTQ3D1(A%FORG,A%LINE_VALUE(2)-A%LINE_VALUE(1), &
                           A%NLMAX,A%MODE_PRD, &
                           A%LINE_PAT,KA,LINE_RGB_SUB,A%BEV_TYPE/10)
          ELSE
             IF(A%FMIN*A%FMAX.LT.0.0) THEN
                CALL CONTQ3D1(A%FORG,A%LINE_VALUE(2)-A%LINE_VALUE(1), &
                              A%NLMAX,A%MODE_PRD, &
                              A%LINE_PAT,KA,R2W2B,A%BEV_TYPE/10)
                CALL CONTQ3D1(A%FORG,A%LINE_VALUE(1)-A%LINE_VALUE(2), &
                              A%NLMAX,A%MODE_PRD, &
                              A%LINE_PAT,KA,R2W2B,A%BEV_TYPE/10)
             ELSEIF(A%FMIN.LT.0.0) THEN
                CALL CONTQ3D1(A%FORG,A%LINE_VALUE(1)-A%LINE_VALUE(2), &
                              A%NLMAX,A%MODE_PRD, &
                              A%LINE_PAT,KA,W2G2B,A%BEV_TYPE/10)
             ELSE
                CALL CONTQ3D1(A%FORG,A%LINE_VALUE(2)-A%LINE_VALUE(1), &
                              A%NLMAX,A%MODE_PRD, &
                              A%LINE_PAT,KA,R2Y2W,A%BEV_TYPE/10)
             ENDIF
          END IF
          DEALLOCATE(KA)
       END SELECT

       CALL GRF_FRAME3D(A)

    END SELECT
    CALL SETRGB(0.0,0.0,0.0)
    CALL SETLNW(line_width)
    RETURN
  END SUBROUTINE GRF2D_EXEC

END MODULE grf2dexec
