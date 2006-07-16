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
      DIMENSION GVDX(NRM,NTM)
      DIMENSION GTX(NTM)
      DIMENSION KA(NRM,NTM)
C
      NG=NGT
      NSTEP=(NG-1)/100+1
      NTLMAX=NGT/NSTEP+1
      DO NTL=1,NTLMAX-1
         NGTL=NSTEP*(NTL-1)+1
         GTX(NTL)=GT(NGTL)
         DO NR=1,NRMAX
            GVDX(NR,NTL)=GVD(NR,NGTL)
         ENDDO
      ENDDO
      NTL=NTLMAX
      NGTL=NGT
         GTX(NTL)=GT(NGTL)
         DO NR=1,NRMAX
            GVDX(NR,NTL)=GVD(NR,NGTL)
         ENDDO
C
      GX1=3.0
      GX2=20.0
      GY1=2.0
      GY2=17.0
C
      CALL PAGES
      CALL TRGR3D(GX1,GX2,GY1,GY2,GRM,GTX,GVDX,NRM,NRMAX,NTLMAX,KA,
     &            STR,KV,2+INQ)
      CALL PAGEE
C
      RETURN
      END
C
C     ***********************************************************
C
C           SUBPROGRAM FOR 2D PROFILE
C
C     ***********************************************************
C
      SUBROUTINE TRGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,
     &     NXM,NXMAX,NYMAX,KA,STR,KV,MODE)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
      DIMENSION GX(NXMAX),GY(NYMAX),GZ(NXM,NYMAX)
      DIMENSION KA(8,NXM,NYMAX)
      CHARACTER STR*80,KT*80,KDL*1,KV*80
      EXTERNAL R2W2B
C      EXTERNAL R2G2B
C
      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.035)
      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLIN(-1,-1,7)
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
      GZVAL=0.5*(ABS(GZMIN)+ABS(GZMAX))
      IF((GZVAL.LE.1.D-6).OR.
     &   (ABS(GZMAX-GZMIN)/GZVAL.LT.1.D-6)) THEN
         WRITE(6,*) GZMIN,GZMAX
         CALL GUFLSH
         RETURN
      ENDIF
         WRITE(6,*) GZMIN,GZMAX
         CALL GUFLSH
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
C
      IF(MOD(MODE/4,2).EQ.1) THEN
         CALL CHMODE
         WRITE(6,*) '## TRGR : XMIN,XMAX,YMIN,YMAX = ',
     &              GXMIN,GXMAX,GYMIN,GYMAX
         READ(5,*) GXMIN,GXMAX,GYMIN,GYMAX
         CALL GRMODE
      ENDIF
C
      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSTEPY)
      CALL GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GSTEPZ)
         WRITE(6,*) GSZMIN,GSZMAX,GSTEPZ
         CALL GUFLSH
C
      IF(GXMIN*GXMAX.LE.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GSXMIN
      ENDIF
      IF(GYMIN*GYMAX.LE.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=GSYMIN
      ENDIF
C
      GXL=10.0*1.5
      GYL=20.0*1.5
      GZL=10.0*1.5
      CALL GDEFIN3D(GX1,GX2,GY1,GY2,GXL,GYL,GZL)
      GPHI=-60.0
      GTHETA=65.0
      GRADIUS=100.0
      GOX=0.5*(GSXMIN+GSXMAX)
      GOY=0.5*(GYMIN+GYMAX)
      GOZ=0.5*(GSZMIN+GSZMAX)
      CALL GVIEW3D(GPHI,GTHETA,GRADIUS,1.0,1,GOX,GOY,GOZ)
      CALL GDATA3D1(GZ,NXM,NXMAX,NYMAX,
     &              GSXMIN,GSXMAX,GSYMIN,GSYMAX,GSZMIN,GSZMAX)
      CALL SETCHS(0.2,0.0)
      CALL SETLIN(0,0,7)
C
      CALL GSCALE3DX(GSXMIN,GSTEPX,0.3,0)
      CALL GSCALE3DY(GSYMIN,GSTEPY,0.3,0)
      CALL GSCALE3DZ(GSZMIN,GSTEPZ,0.3,0)
      CALL GVALUE3DX(GSXMIN,GSTEPX,1,1)
      CALL GVALUE3DY(GSYMIN,GSTEPY,1,1)
      CALL GVALUE3DZ(GSZMIN,GSTEPZ,2,-2)
C
C      CALL GTEXTX(GX1-0.3,0.5*(GY1+GY2),
C     &              '@TIME (sec)@',
C     &              2)
C      CALL GTEXTX(0.5*(GX1+GX2),GY1-0.3,
C     &              '@R@',
C     &              2)
C      CALL GTEXTX(GX1,GY2+0.1,
C     &              KV,
C     &              0)
C
      CALL CPLOT3D1(7,R2W2B)
C      CALL CONTQ3D1(GZMIN,0.1*(GZMAX-GZMIN),11,0,0,KA,R2W2B,2)
C      CALL PERSE3D(3,1)
      CALL GAXIS3D(0)
C      CALL GDrawBack3D(0.5, 0.5, 0.5)
C
      CALL SETLIN(0,0,7)
      RETURN
      END
