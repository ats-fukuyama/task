C     $Id$
C  
C     **************************************************************
C
C           GRAPHIC 3D : UNIVERSAL ROUTINE
C
C     **************************************************************
C
      SUBROUTINE TRGRUR(GVD,STR,KV,INQ)
C
      INCLUDE 'trcomm.inc'
      CHARACTER STR*80,KV*80
      DIMENSION GVD(NRM,NTM)
C
      GX1=3.0
      GX2=18.0
      GY1=2.0
      GY2=17.0
C
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      CALL PAGES
      CALL GSGLENABLELIGHTING
      CALL TRGR3D(GX1,GX2,GY1,GY2,GRM,GT,GVD,NRM,NRMAX,NGT,
     &            STR,KV,2+INQ)
      CALL PAGEE
      IF(RHOA.NE.1.D0) NRMAX=NRAMAX
C
      RETURN
      END
C
C     ***********************************************************
C
C           SUBPROGRAM FOR 3D PROFILE
C
C     ***********************************************************
C
      SUBROUTINE TRGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,
     &     NXM,NXMAX,NYMAX,STR,KV,MODE)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
      DIMENSION GX(NXMAX),GY(NYMAX),GZ(NXM,NYMAX)
      CHARACTER STR*80,KT*80,KDL*1,KV*80
      EXTERNAL R2G2B
C
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(0,0,7)
      KDL=STR(1:1)
      I=2
    1 IF(STR(I:I).EQ.KDL.OR.I.EQ.80) GOTO 2
         KT(I-1:I-1)=STR(I:I)
         I=I+1
      GOTO 1
C
    2 CALL MOVE(GX1,GY2+0.2)
      CALL TEXT(KT,I-2)
C
      CALL GMNMX2(GZ,NXM,1,NXMAX,1,1,NYMAX,1,GZMIN,GZMAX)
      IF(ABS(GZMAX-GZMIN).LT.1.D-6) THEN
         GZMIN=GZMIN-0.999E-6
         GZMAX=GZMAX+1.000E-6
      ENDIF
C
      IF(MOD(MODE,2).EQ.0) THEN
         IF(GZMIN.GE.0.0) THEN
            GZMIN=0.0
         ELSEIF(GZMAX.LE.0.0) THEN
            GZMAX=0.0
         ENDIF
      ENDIF
C
      CALL GMNMX1(GX,1,NXMAX,1,GXMIN,GXMAX)
      CALL GMNMX1(GY,1,NYMAX,1,GYMIN,GYMAX)
      IF(ABS(GXMAX-GXMIN).LT.1.D-6) THEN
         GXMIN=GXMIN-0.999E-6
         GXMAX=GXMAX+1.000E-6
      ENDIF
C
C      GXMIN=GX(1)
C      GXMAX=GX(NXMAX)
C
      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSTEPY)
      CALL GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GSTEPZ)
C
      IF(MOD(MODE,2).EQ.0) THEN
         IF(GZMIN.GE.0.0) THEN
            GSZMIN=0.0
         ELSEIF(GZMAX.LE.0.0) THEN
            GSZMAX=0.0
         ENDIF
      ENDIF
      GZMIN=GSZMIN
      GZMAX=GSZMAX
      IF(MOD(MODE/4,2).EQ.1) THEN
         CALL CHMODE
         WRITE(6,*) '## TRGR : XMIN,XMAX,YMIN,YMAX = ',
     &              GXMIN,GXMAX,GYMIN,GYMAX
         READ(5,*) GXMIN,GXMAX,GYMIN,GYMAX
         CALL GRMODE
      ENDIF
      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSTEPY)
      CALL GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GSTEPZ)
C
C      IF(GXMIN*GXMAX.LE.0.0) THEN
C         GXORG=0.0
C      ELSE
C         GXORG=GSXMIN
C      ENDIF
C      IF(GYMIN*GYMAX.LE.0.0) THEN
C         GYORG=0.0
C      ELSE
C         GYORG=GSYMIN
C      ENDIF
C
      CALL GDEFIN(GX1,GX2,GY1,GY2,
     &            GSXMIN,GSXMAX,GSYMIN,GSYMAX)
      GXL=10.0*1.5
      GYL=20.0*1.5
      GZL=10.0*1.5
      CALL GDEFIN3D(GXL,GYL,GZL,GSXMIN,GSXMAX,GYMIN,GYMAX,
     &     GSZMIN,GSZMAX)
      GPHI=-20.0
      GTHETA=60.0
      GRADIUS=15.0
      GOX=0.5*(GSXMIN+GSXMAX)
      GOY=0.5*(GYMIN+GYMAX)
      GOZ=0.5*(GSZMIN+GSZMAX)
      CALL GVIEW3D(GPHI,GTHETA,GRADIUS,GOX,GOY,GOZ)
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(0,0,4)
C
C      CALL GSCALE3DX(GSXMIN,GSTEPX,0.3,0)
C      CALL GSCALE3DY(GT(1),GSTEPT,0.3,0)
C      CALL GSCALE3DZ(GSYMIN,GSTEPY,0.3,10)
C      CALL GVALUE3DX(GSXMIN,GSTEPX,1,1)
C      CALL GVALUE3DY(GT(1),GSTEPT,1,1)
C      CALL GVALUE3DZ(GSYMIN,GSTEPY,11,-2)
      CALL GSCALE3DX(GSXMIN,GSTEPX,0.3,0)
      CALL GSCALE3DY(GSYMIN,GSTEPY,0.3,0)
      CALL GSCALE3DZ(GSZMIN,GSTEPZ,0.3,10)
      CALL GVALUE3DX(GSXMIN,GSTEPX,1,1)
      CALL GVALUE3DY(GSYMIN,GSTEPY,1,1)
      CALL GVALUE3DZ(GSZMIN,GSTEPZ,11,-2)
C
      CALL Set3DTextBaseLine(0.0, 1.0, 0.0, -1.0, 0.0, 0.0)
      CALL GTEXTX3D(GSXMAX+0.15*(GSXMAX-GSXMIN),
     &              0.5*(GY(1)+GY(NYMAX)),
     &              GSYMIN,
     &              '@TIME (sec)@',
     &              2)
      CALL Set3DTextBaseLine(0.0, 1.0, 0.0, -1.0, 0.0, 0.0)
      CALL GTEXTX3D(0.5*(GSXMIN+GSXMAX),
     &              GY(1)+0.1*(GY(1)-GY(NYMAX)),
     &              GSYMIN,
     &              '@RHO@',
     &              2)
      CALL Set3DTextBaseLine(0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
      CALL GTEXTX3D(GSXMIN,
     &              GY(1)+0.05*(GY(1)-GY(NYMAX)),
     &              GSZMAX+0.1*(GSZMAX-GSZMIN),
     &              KV,
     &              2)
C
      CALL PERS3D1(GZ,NXM,NXMAX,NYMAX,-27,R2G2B)
      CALL GAxis3D(0)
      CALL GDrawBack3D(0.5, 0.5, 0.5)
C
      CALL SETLIN(0,0,4)
      RETURN
      END
