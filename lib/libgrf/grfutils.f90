MODULE grfutils

  PRIVATE
  PUBLIC grf_title,grf_frame,grfut1,grfut2,grfut3,grfut4,setrgba, &
         grf2da,grf2db,grf2dc

CONTAINS

! ***** DRAW GRAPH TITLES *****

  SUBROUTINE GRF_TITLE(A)

    USE grftype, ONLY: grf_attr_type
    IMPLICIT NONE

    TYPE(grf_attr_type),INTENT(IN):: A
    REAL(4):: SIZE
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

  SUBROUTINE GRF_FRAME(A)

    USE grftype, ONLY: grf_attr_type
    IMPLICIT NONE

    TYPE(grf_attr_type),INTENT(IN):: A
    REAL(4):: SIZE
    INTEGER:: NG,NL

    CALL GDEFIN(A%GPXMIN,A%GPXMAX,A%GPYMIN,A%GPYMAX, &
                A%XMIN,A%XMAX,A%FMIN,A%FMAX)

    CALL SETRGBA(A%FRAME_RGB(1:3))
    CALL SETLNW(A%FRAME_THICKNESS)
    CALL GFRAME

    CALL SETRGBA(A%SCALE_RGB(1:3))
    CALL SETLNW(A%SCALE_THICKNESS)
    CALL SETCHS(A%VALUE_SIZE,0.0)
    CALL SETFNT(A%VALUE_FONT)

    IF(A%MODE_LS.EQ.0.OR.A%MODE_LS.EQ.2) THEN
       IF(A%XSCALE_ZERO.NE.0) &
            CALL GSCALE(0.0,2*(A%XMAX-A%XMIN),0.0,0.0, &
                        FLOAT(A%XSCALE_ZERO_PAT),0)
       CALL GSCALE(A%XORG,A%XSCALE_STEP,0.0,0.0,A%XSCALE_SIZE,A%XSCALE_TYPE)
       CALL SETRGBA(A%VALUE_RGB(1:3))
       CALL GVALUE(A%XORG,A%XVALUE_STEP*A%XSCALE_STEP,0.0,0.0,A%XVALUE_TYPE)
    ELSE
       CALL GSCALL(A%XORG,A%XSCALE_LTYPE,0.0,0,A%XSCALE_SIZE,A%XSCALE_TYPE)
       CALL SETRGBA(A%VALUE_RGB(1:3))
       CALL GVALUL(A%XORG,A%XVALUE_LTYPE,0.0,0,A%XVALUE_TYPE)
    ENDIF

    IF(A%MODE_LS.EQ.0.OR.A%MODE_LS.EQ.1) THEN
       IF(A%FSCALE_ZERO.NE.0) &
            CALL GSCALE(0.0,0.0,0.0,2*(A%FMAX-A%FMIN), &
                        FLOAT(A%FSCALE_ZERO_PAT),0)
       CALL GSCALE(0.0,0.0,A%FORG,A%FSCALE_STEP,A%FSCALE_SIZE,A%FSCALE_TYPE)
       CALL SETRGBA(A%VALUE_RGB(1:3))
       CALL GVALUE(0.0,0.0,A%FORG,A%FVALUE_STEP*A%FSCALE_STEP,A%FVALUE_TYPE)
    ELSE
       CALL GSCALL(0.0,0,A%FORG,A%FSCALE_LTYPE,A%FSCALE_SIZE,A%FSCALE_TYPE)
       CALL SETRGBA(A%VALUE_RGB(1:3))
       CALL GVALUL(0.0,0,A%FORG,A%FVALUE_LTYPE,A%FVALUE_TYPE)
    ENDIF

  END SUBROUTINE GRF_FRAME


!     ***** EVALUATE GXMIN,GXMAX,GYMIN,GYMAX *****

  SUBROUTINE GRFUT1(GX,NXMAX,GXMIN,GXMAX,NXMIN,NXSTEP)

    IMPLICIT NONE

    REAL(4),INTENT(IN):: GX(NXMAX)
    INTEGER,INTENT(IN):: NXMAX
    REAL(4),INTENT(OUT):: GXMIN,GXMAX
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
    IF(ABS(GXMAX-GXMIN).LT.1.D-6) THEN
       GXMIN=GXMIN-0.999E-6
       GXMAX=GXMAX+1.000E-6
    ENDIF

    RETURN
  END SUBROUTINE GRFUT1

!     ***** EVALUATE GXMIN,GXMAX,GYMIN,GYMAX *****

  SUBROUTINE GRFUT2(GY,NXM,NXMAX,NYMAX,GYMIN,GYMAX, &
                       NXMIN,NXSTEP,NYMIN,NYSTEP)

    IMPLICIT NONE

    REAL(4),INTENT(IN):: GY(NXM,NYMAX)
    INTEGER,INTENT(IN):: NXM,NXMAX,NYMAX
    REAL(4),INTENT(OUT):: GYMIN,GYMAX
    INTEGER,INTENT(IN),OPTIONAL:: NXMIN,NXSTEP,NYMIN,NYSTEP
    INTEGER:: NXMINL,NXSTEPL,NYMINL,NYSTEPL
    REAL(4):: GYNORM

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
    IF(ABS(GYMAX-GYMIN).LT.1.D-6*GYNORM) THEN
       GYMIN=GYMIN-0.999E-6*GYNORM
       GYMAX=GYMAX+1.000E-6*GYNORM
    ENDIF

    RETURN
  END SUBROUTINE GRFUT2

! ***** EVALUATE GSMIN,GSMAX,GSTEP,GORG *****

  SUBROUTINE GRFUT3(GMIN,GMAX,GSMIN,GSMAX,GSTEP,GORG)

    IMPLICIT NONE
    REAL(4),INTENT(IN):: GMIN,GMAX
    REAL(4),INTENT(OUT):: GSMIN,GSMAX,GSTEP,GORG

    CALL GQSCAL(GMIN,GMAX,GSMIN,GSMAX,GSTEP)

    IF(GMIN*GMAX.LE.0.0) THEN
       GORG=0.0
    ELSE
       GORG=(NINT(GSMIN/(2*GSTEP))+1)*2*GSTEP
    ENDIF
    RETURN
  END SUBROUTINE GRFUT3

! ***** SET GP FOR NGP *****

  SUBROUTINE GRFUT4(NGP,GP)

    IMPLICIT NONE

    INTEGER NGP
    REAL*4 GP(4)
    INTEGER IX,IY
    REAL*4 X0,XL,Y0,YL

    IF(NGP.EQ.0) THEN
       GP(1)= 3.0
       GP(2)=23.0
       GP(3)= 2.0
       GP(4)=17.0
    ELSE
       IF(NGP.LE.4) THEN
          X0= 3.0
          XL=12.0
          Y0= 1.75
          YL= 8.5
          IX=MOD(NGP-1,2)
          IY=1-(NGP-1)/2
       ELSEIF(NGP.LE.13) THEN
          X0= 3.0*0.666
          XL=12.0*0.666
          Y0= 1.75*0.666
          YL= 8.5*0.666
          IX=MOD(NGP-5,3)
          IY=2-(NGP-5)/3
       ELSEIF(NGP.LE.29) THEN
          X0= 3.0*0.5
          XL=12.0*0.5
          Y0= 1.75*0.5
          YL= 8.5*0.5
          IX=MOD(NGP-14,4)
          IY=3-(NGP-14)/4
       ENDIF
       GP(1)=IX*XL+X0
       GP(2)=IX*XL+XL
       GP(3)=IY*YL+Y0
       GP(4)=IY*YL+YL
    ENDIF
    RETURN
  END SUBROUTINE GRFUT4

! ***** SET RGB BY ARRAY *****

  SUBROUTINE SETRGBA(RGB)

    IMPLICIT NONE
    REAL(4),DIMENSION(3),INTENT(IN):: RGB

    CALL SETRGB(RGB(1),RGB(2),RGB(3))
    RETURN
  END SUBROUTINE SETRGBA

!     ***** DRAW 2D GRAPH *****

      SUBROUTINE GRF2DA(NGP,GX,GY,GZ,NXM,NXMAX,NYMAX,STR,MODE_LS,IPRD)

      IMPLICIT REAL*8 (A-F,H,O-Z)

      DIMENSION GP(4),GX(NXMAX),GY(NYMAX),GZ(NXM,NYMAX)
      CHARACTER STR*(*)

      IF(NGP.LT.0.OR.NGP.GT.29) THEN
         WRITE(6,*) 'XX GRF2DC: INVALID NGP'
         RETURN
      ENDIF
      CALL GRFUT4(NGP,GP)

      CALL GRFUT1(GX,NXMAX,GXMIN,GXMAX)
      CALL GRFUT1(GY,NYMAX,GYMIN,GYMAX)
      CALL GRFUT2(GZ,NXM,NXMAX,NYMAX,GFMIN,GFMAX)

      CALL GRFUT3(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSTEP,GXORG)
      CALL GRFUT3(GYMIN,GYMAX,GSYMIN,GSYMAX,GYSTEP,GYORG)
      CALL GRFUT3(GFMIN,GFMAX,GSFMIN,GSFMAX,GFSTEP,GFORG)
      GFSTEP=0.5*GFSTEP

      CALL GRF2DAX(GP, &
                   GXMIN,GXMAX,GXSTEP,GXORG, &
                   GYMIN,GYMAX,GYSTEP,GYORG, &
                   GFMIN,GFMAX,GFSTEP,GFORG, &
                   GX,GY,GZ,NXM,NXMAX,NYMAX,STR,MODE_LS,IPRD)

      RETURN
    END SUBROUTINE GRF2DA

!     ***** DRAW 2D GRAPH *****

      SUBROUTINE GRF2DAX(GP, &
                         GXMIN,GXMAX,GXSTEP,GXORG, &
                         GYMIN,GYMAX,GYSTEP,GYORG, &
                         GFMIN,GFMAX,GFSTEP,GFORG, &
                         GX,GY,GZ,NXM,NXMAX,NYMAX, &
                         STR,MODE_LS,IPRD)

      IMPLICIT REAL*8 (A-F,H,O-Z)

      DIMENSION GP(4),GX(NXMAX),GY(NYMAX),GZ(NXM,NYMAX)
      CHARACTER STR*(*),KT*80,KDL*1
      DIMENSION GZL(101),GRGB(3,101)

      GL=GP(2)-GP(1)
      IF(GL.GT.7.9) then
         GFACTOR=1.0
      ELSEIF(GL.GT.7.9*0.666) then
         GFACTOR=0.666
      ELSE
         GFACTOR=0.5
      ENDIF

      CALL SETCHS(0.3*GFACTOR,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.035*GFACTOR)
      CALL SETRGB(0.0,0.0,0.0)

      KDL=STR(1:1)
      I=2
    1 IF(STR(I:I).EQ.KDL.OR.I.EQ.LEN(STR)) GOTO 2
         KT(I-1:I-1)=STR(I:I)
         I=I+1
      GOTO 1

    2 CALL MOVE(GP(1),GP(4)+0.3)
      CALL TEXT(KT,I-2)

      CALL GDEFIN(GP(1),GP(2),GP(3),GP(4), &
                  GXMIN,GXMAX,GYMIN,GYMAX)

      NSTEP=NINT(ABS((GFMAX-GFMIN)/GFSTEP))+1
      IF(NSTEP.GT.101) NSTEP=101
      GFSTEPL=(GFMAX-GFMIN)/(NSTEP-1)
      DO I=1,NSTEP
         GZL(I)=GFMIN+GFSTEPL*(I-1)
      ENDDO

      IF(GFMIN*GFMAX.GE.0.D0) THEN
         IF(GFORG.GE.0.D0) THEN
            DO I=1,NSTEP
               GFACTOR=FLOAT(I-1)/FLOAT(NSTEP-1)
               CALL R2Y2W(GFACTOR,GRGB(1,I))
            ENDDO
            CALL CONTF2(GZ,GX,GY,NXM,NXMAX,NYMAX,GZL,GRGB,NSTEP,IPRD)
         ELSE
            DO I=1,NSTEP
               GFACTOR=FLOAT(I-1)/FLOAT(NSTEP-1)
               CALL W2G2B(GFACTOR,GRGB(1,I))
            ENDDO
            CALL CONTF2(GZ,GX,GY,NXM,NXMAX,NYMAX,GZL,GRBG,NSTEP,IPRD)
         ENDIF
      ELSE
         DO I=1,NSTEP
            GFACTOR=FLOAT(I-1)/FLOAT(NSTEP-1)
            CALL R2W2B(GFACTOR,GRGB(1,I))
         ENDDO
         CALL CONTF2(GZ,GX,GY,NXM,NXMAX,NYMAX,GZL,GRGB,NSTEP,IPRD)
      ENDIF

      CALL SETRGB(0.0,0.0,0.0)
      CALL GFRAME
      IF(MODE_LS.EQ.0.OR.MODE_LS.EQ.2) THEN
         CALL GSCALE(0.0,2*(GXMAX-GXMIN),0.0,0.0,0.0,0)
         CALL GSCALE(GXORG,GXSTEP,0.0,0.0,0.1,9)
         CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
      ELSE
         CALL GSCALL(GXORG,9,0.0,0,0.1,9)
         CALL GVALUL(GXORG,1,0.0,0,NGULEN(2*GXSTEP))
      ENDIF
      IF(MODE_LS.EQ.0.OR.MODE_LS.EQ.1) THEN
         CALL GSCALE(0.0,0.0,0.0,2*(GYMAX-GYMIN),0.0,0)
         CALL GSCALE(0.0,0.0,GYORG,GYSTEP,0.1,9)
         CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))
      ELSE
         CALL GSCALL(0.0,0,GYORG,9,0.1,9)
         CALL GVALUL(0.0,0,GYORG,1,NGULEN(2*GYSTEP))
      ENDIF

      CALL MOVE(0.5*(GP(1)+GP(2))-0.2*32*GFACTOR,GP(3)-1.5*GFACTOR)
      CALL TEXT('MIN  =',6)
      CALL NUMBR(GFMIN ,'(1PE12.4)',12)
      CALL TEXT('     MAX  =',11)
      CALL NUMBR(GFMAX ,'(1PE12.4)',12)
      CALL TEXT('     STEP =',11)
      CALL NUMBR(GFSTEP,'(1PE12.4)',12)

      CALL SETRGB(0.0,0.0,0.0)
      RETURN
    END SUBROUTINE GRF2DAX

!     ***** DRAW 2D GRAPH *****

      SUBROUTINE GRF2DB(NGP,GX,GY,GZ,NXM,NXMAX,NYMAX,STR)

      IMPLICIT REAL*8 (A-F,H,O-Z)

      DIMENSION GP(4),GX(NXMAX),GY(NYMAX),GZ(NXM,NYMAX)
      DIMENSION GP3D(6)
      CHARACTER STR*(*)

      IF(NGP.LT.0.OR.NGP.GT.29) THEN
         WRITE(6,*) 'XX GRF2DB: INVALID NGP'
         RETURN
      ENDIF
      CALL GRFUT4(NGP,GP)

      GP3D(1)=10.0*1.5
      GP3D(2)=20.0*1.5
      GP3D(3)=10.0*1.5
      GP3D(4)=-60.0
      GP3D(5)=65.0
      GP3D(6)=100.0

      CALL GRFUT1(GX,NXMAX,GXMIN,GXMAX)
      CALL GRFUT1(GY,NYMAX,GYMIN,GYMAX)
      CALL GRFUT2(GZ,NXM,NXMAX,NYMAX,GFMIN,GFMAX)

      CALL GRFUT3(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSTEP,GXORG)
      CALL GRFUT3(GYMIN,GYMAX,GSYMIN,GSYMAX,GYSTEP,GYORG)
      CALL GRFUT3(GFMIN,GFMAX,GSFMIN,GSFMAX,GFSTEP,GFORG)

      CALL GRF2DBX(GP, &
                   GXMIN,GXMAX,GXSTEP,GXORG, &
                   GYMIN,GYMAX,GYSTEP,GYORG, &
                   GFMIN,GFMAX,GFSTEP,GFORG, &
                   GX,GY,GZ,NXM,NXMAX,NYMAX, &
                   STR,GP3D)
      RETURN
    END SUBROUTINE GRF2DB

!     ***** DRAW 2D GRAPH *****

      SUBROUTINE GRF2DBX(GP, &
                         GXMIN,GXMAX,GXSTEP,GXORG, &
                         GYMIN,GYMAX,GYSTEP,GYORG, &
                         GFMIN,GFMAX,GFSTEP,GFORG, &
                         GX,GY,GZ,NXM,NXMAX,NYMAX, &
                         STR,GP3D)

      IMPLICIT REAL*8 (A-F,H,O-Z)

      DIMENSION GP(4),GX(NXMAX),GY(NYMAX),GZ(NXM,NYMAX)
      DIMENSION GP3D(6)
      CHARACTER STR*(*),KT*80,KDL*1
      EXTERNAL R2W2B,W2G2B,R2Y2W

      GL=GP(2)-GP(1)
      IF(GL.GT.7.9) then
         GFACTOR=1.0
      ELSEIF(GL.GT.7.9*0.666) then
         GFACTOR=0.666
      ELSE
         GFACTOR=0.5
      ENDIF

      CALL SETCHS(0.3*GFACTOR,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.035*GFACTOR)
      CALL SETRGB(0.0,0.0,0.0)

      KDL=STR(1:1)
      I=2
    1 IF(STR(I:I).EQ.KDL.OR.I.EQ.LEN(STR)) GOTO 2
         KT(I-1:I-1)=STR(I:I)
         I=I+1
      GOTO 1

    2 CALL MOVE(GP(1),GP(4)+0.3)
      CALL TEXT(KT,I-2)

      GXL    =GP3D(1)
      GYL    =GP3D(2)
      GZL    =GP3D(3)
      GPHI   =GP3D(4)
      GTHETA =GP3D(5)
      GRADIUS=GP3D(6)
      GOX=0.5*(GXMIN+GXMAX)
      GOY=0.5*(GYMIN+GYMAX)
      GOZ=0.5*(GFMIN+GFMAX)

      CALL GDEFIN3D(GP(1),GP(2),GP(3),GP(4),GXL,GYL,GZL)
      CALL GVIEW3D(GPHI,GTHETA,GRADIUS,1.0,1,GOX,GOY,GOZ)
      CALL GDATA3D1(GZ,NXM,NXMAX,NYMAX, &
                    GXMIN,GXMAX,GYMIN,GYMAX,GFMIN,GFMAX)

      CALL SETCHS(0.2,0.0)
      CALL SETLIN(0,0,7)

      IF(NXMAX.GT.100.OR.NYMAX.GT.100) THEN
         ID=4
      ELSE
         ID=7
      ENDIF

      IF(GFMIN*GFMAX.LT.0.0) THEN
         CALL CPLOT3D1(ID,R2W2B)
      ELSEIF(GFMIN.LT.0.0) THEN
         CALL CPLOT3D1(ID,W2G2B)
      ELSE
         CALL CPLOT3D1(ID,R2Y2W)
      ENDIF
      CALL GAXIS3D(0)

      CALL SETLIN(0,0,7)
      CALL GSCALE3DX(GXORG,GXSTEP,0.3,0)
      CALL GSCALE3DY(GYORG,GYSTEP,0.3,0)
      CALL GSCALE3DZ(GFORG,GFSTEP,0.3,0)
      CALL GVALUE3DX(GXORG,2*GXSTEP,1,NGULEN(2*GXSTEP))
      CALL GVALUE3DY(GYORG,2*GYSTEP,1,NGULEN(2*GYSTEP))
      CALL GVALUE3DZ(GFORG,2*GFSTEP,2,NGULEN(2*GFSTEP))
      RETURN
    END SUBROUTINE GRF2DBX

!     ***** DRAW 2D GRAPH *****

      SUBROUTINE GRF2DC(NGP,GX,GY,GZ,NXM,NXMAX,NYMAX,STR,MODE_LS)

      IMPLICIT REAL*8 (A-F,H,O-Z)

      DIMENSION GP(4),GX(NXMAX),GY(NYMAX),GZ(NXM,NYMAX)
      CHARACTER STR*(*)

      IF(NGP.LT.0.OR.NGP.GT.29) THEN
         WRITE(6,*) 'XX GRF2DC: INVALID NGP'
         RETURN
      ENDIF
      CALL GRFUT4(NGP,GP)

      CALL GRFUT1(GX,NXMAX,GXMIN,GXMAX)
      CALL GRFUT1(GY,NYMAX,GYMIN,GYMAX)
      CALL GRFUT2(GZ,NXM,NXMAX,NYMAX,GFMIN,GFMAX)

      CALL GRFUT3(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSTEP,GXORG)
      CALL GRFUT3(GYMIN,GYMAX,GSYMIN,GSYMAX,GYSTEP,GYORG)
      CALL GRFUT3(GFMIN,GFMAX,GSFMIN,GSFMAX,GFSTEP,GFORG)
      GFSTEP=0.2*GFSTEP

      CALL GRF2DCX(GP, &
                   GXMIN,GXMAX,GXSTEP,GXORG, &
                   GYMIN,GYMAX,GYSTEP,GYORG, &
                   GFMIN,GFMAX,GFSTEP,GFORG, &
                   GX,GY,GZ,NXM,NXMAX,NYMAX, &
                   STR,MODE_LS)

      RETURN
    END SUBROUTINE GRF2DC

!     ***** DRAW 2D GRAPH *****

      SUBROUTINE GRF2DCX(GP, &
                         GXMIN,GXMAX,GXSTEP,GXORG, &
                         GYMIN,GYMAX,GYSTEP,GYORG, &
                         GFMIN,GFMAX,GFSTEP,GFORG, &
                         GX,GY,GZ,NXM,NXMAX,NYMAX, &
                         STR,MODE_LS)

      IMPLICIT REAL*8 (A-F,H,O-Z)

      DIMENSION GP(4),GX(NXMAX),GY(NYMAX),GZ(NXM,NYMAX)
      DIMENSION KA(8,NXM,NYMAX)
      CHARACTER STR*(*),KT*80,KDL*1

      GL=GP(2)-GP(1)
      IF(GL.GT.7.9) then
         GFACTOR=1.0
      ELSEIF(GL.GT.7.9*0.666) then
         GFACTOR=0.666
      ELSE
         GFACTOR=0.5
      ENDIF

      CALL SETCHS(0.3*GFACTOR,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.035*GFACTOR)
      CALL SETRGB(0.0,0.0,0.0)

      KDL=STR(1:1)
      I=2
    1 IF(STR(I:I).EQ.KDL.OR.I.EQ.LEN(STR)) GOTO 2
         KT(I-1:I-1)=STR(I:I)
         I=I+1
      GOTO 1

    2 CALL MOVE(GP(1),GP(4)+0.3)
      CALL TEXT(KT,I-2)

      CALL GDEFIN(GP(1),GP(2),GP(3),GP(4), &
                  GXMIN,GXMAX,GYMIN,GYMAX)

      IF(GFMIN*GFMAX.GE.0.D0) THEN
         NSTEP=NINT(ABS((GFMAX-GFMIN)/GFSTEP))+1
         IF(GFORG.GE.0.D0) THEN
            CALL SETRGB(1.0,0.0,0.0)
            CALL CONTQ2(GZ,GX,GY,NXM,NXMAX,NYMAX, &
                        GFORG,GFSTEP,NSTEP,0,0,KA)
         ELSE
            CALL SETRGB(0.0,0.0,1.0)
            CALL CONTQ2(GZ,GX,GY,NXM,NXMAX,NYMAX, &
                        GFORG,GFSTEP,NSTEP,0,3,KA)
         ENDIF
      ELSE
         CALL SETRGB(1.0,0.0,0.0)
         NSTEP=NINT(ABS(GFMAX/GFSTEP))+1
         CALL CONTQ2(GZ,GX,GY,NXM,NXMAX,NYMAX, &
                     GFSTEP,GFSTEP,NSTEP,0,0,KA)
         CALL SETRGB(0.0,1.0,0.0)
         NSTEP=NINT(ABS(GFMAX/GFSTEP))+1
         CALL CONTQ2(GZ,GX,GY,NXM,NXMAX,NYMAX, &
                     0.0,GFSTEP,1,0,0,KA)
         CALL SETRGB(0.0,0.0,1.0)
         NSTEP=NINT(ABS(GFMIN/GFSTEP))+1
         CALL CONTQ2(GZ,GX,GY,NXM,NXMAX,NYMAX, &
                     -GFSTEP,-GFSTEP,NSTEP,0,3,KA)
      ENDIF

      CALL SETRGB(0.0,0.0,0.0)
      CALL GFRAME

      IF(MODE_LS.EQ.0.OR.MODE_LS.EQ.2) THEN
         CALL GSCALE(0.0,2*(GXMAX-GXMIN),0.0,0.0,0.0,0)
         CALL GSCALE(GXORG,GXSTEP,0.0,0.0,0.1,9)
         CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
      ELSE
         CALL GSCALL(GXORG,9,0.0,0,0.1,9)
         CALL GVALUL(GXORG,1,0.0,0,NGULEN(2*GXSTEP))
      ENDIF
      IF(MODE_LS.EQ.0.OR.MODE_LS.EQ.1) THEN
         CALL GSCALE(0.0,0.0,0.0,2*(GYMAX-GYMIN),0.0,0)
         CALL GSCALE(0.0,0.0,GYORG,GYSTEP,0.1,9)
         CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))
      ELSE
         CALL GSCALL(0.0,0,GYORG,9,0.1,9)
         CALL GVALUL(0.0,0,GYORG,1,NGULEN(2*GYSTEP))
      ENDIF

      CALL MOVE(0.5*(GP(1)+GP(2))-0.2*32*GFACTOR,GP(3)-1.5*GFACTOR)
      CALL TEXT('MIN  =',6)
      CALL NUMBR(GFMIN ,'(1PE12.4)',12)
      CALL TEXT('     MAX  =',11)
      CALL NUMBR(GFMAX ,'(1PE12.4)',12)
      CALL TEXT('     STEP =',11)
      CALL NUMBR(GFSTEP,'(1PE12.4)',12)

      CALL SETRGB(0.0,0.0,0.0)
      RETURN
    END SUBROUTINE GRF2DCX

END MODULE grfutils
