C     $Id$
C
C     ****** CALCULATE RANGE OF WINDOW ******
C
      SUBROUTINE WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
C
      INCLUDE 'wfcomm.inc'
C
      DXLEN= XNDMAX-XNDMIN
      DYLEN= YNDMAX-YNDMIN
      DRATIO=DYLEN/DXLEN
C
      IF(NWXMAX.EQ.0) THEN
         IF(DRATIO.GT.1.25D0) THEN
            NWW=6
         ELSEIF(DRATIO.GT.0.5D0) THEN
            IF(NWMAX.LE.3) then
               NWW=3
            ELSE IF(NWMAX.LE.4) THEN
               NWW=2
            ELSE IF(NWMAX.LE.6) THEN
               NWW=3
            ELSE IF(NWMAX.LE.8) THEN
               NWW=4
            ELSE
               NWW=5
            ENDIF
            NWW=4
         ELSE
            NWW=2
         ENDIF
      ELSE
         NWW=NWXMAX
      ENDIF
C
      PXMIN=0.0D0
      PXMAX=25.6D0
      PYMIN=0.0D0
      PYMAX=14.0D0
C
      NWYMAX=(NWMAX-1)/NWW+1
      PXLEN=(PXMAX-PXMIN)/MIN(NWW,NWMAX)
      PYLEN=(PYMAX-PYMIN)/NWYMAX
      NWX=MOD(NW-1,NWW)+1
      PXMIN=PXMIN+(NWX-1)*PXLEN
      PXMAX=PXMIN+PXLEN
      NWY=NWYMAX-(NW-1)/NWW
      PYMIN=PYMIN+(NWY-1)*PYLEN
      PYMAX=PYMIN+PYLEN
      RETURN
      END
C
C     ****** PLOT 2D STRUCTURE OF FIELD AND POWER ******
C
      SUBROUTINE WFGPPC(NW,NWMAX,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      CHARACTER KWD*(NCHM),KID*1
      DIMENSION KA(4,NGXM,NGYM)
      DIMENSION GAX(NGXM),GAY(NGYM)
C
      IF(NWMAX.LE.0) RETURN
C
      CALL WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
C
      KID=KWD(4:4)
      IF(KID.EQ.'X') THEN
         XMIN = YNDMIN
         XMAX = YNDMAX
         YMIN = ZNDMIN
         YMAX = ZNDMAX
      ELSE IF(KID.EQ.'Y') THEN
         XMIN = XNDMIN
         XMAX = XNDMAX
         YMIN = ZNDMIN
         YMAX = ZNDMAX
      ELSE IF(KID.EQ.'Z') THEN
         XMIN = XNDMIN
         XMAX = XNDMAX
         YMIN = YNDMIN
         YMAX = YNDMAX
      ELSE IF(KID.EQ.'A') THEN
         XMIN = 0.D0
         XMAX = 1.D0
         YMIN = 0.D0
         YMAX = 1.D0
      ENDIF
C
      DXLEN= XMAX-XMIN
      DYLEN= YMAX-YMIN
      DRATIO=DYLEN/DXLEN
      PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
      GXMIN=GUCLIP(XMIN)
      GXMAX=GUCLIP(XMAX)
      GYMIN=GUCLIP(YMIN)
      GYMAX=GUCLIP(YMAX)
      IF(DRATIO.GT.PRATIO) THEN
         GPYMIN=GUCLIP(PYMIN)+1.5
         GPYMAX=GUCLIP(PYMAX)
         PXMID=0.5D0*(PXMIN+PXMAX)
         PYLEN=PYMAX-PYMIN-1.5D0
         PXLEN=PYLEN/DRATIO
         GPXMIN=GUCLIP(PXMID-0.5D0*PXLEN)
         GPXMAX=GUCLIP(PXMID+0.5D0*PXLEN)
      ELSE
         GPXMIN=GUCLIP(PXMIN)+0.3
         GPXMAX=GUCLIP(PXMAX)-0.3
         PYMID=0.5D0*(PYMIN+PYMAX+1.5D0)
         PXLEN=PXMAX-PXMIN-0.6D0
         PYLEN=PXLEN*DRATIO
         GPYMIN=GUCLIP(PYMID-0.5D0*PYLEN)
         GPYMAX=GUCLIP(PYMID+0.5D0*PYLEN)
      ENDIF
C
      DX=(XMAX-XMIN)/(NGXMAX-1)
      DO NGX=1,NGXMAX
        GAX(NGX)=GUCLIP(XMIN+(NGX-1)*DX)
      ENDDO
      DY=(YMAX-YMIN)/(NGYMAX-1)
      DO NGY=1,NGYMAX
        GAY(NGY)=GUCLIP(YMIN+(NGY-1)*DY)
      ENDDO
C
      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSCAL)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSXMAX,GYSCAL)
      IF(GXMIN*GXMAX.LT.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GSXMIN+GXSCAL
      ENDIF
      IF(GYMIN*GYMAX.LT.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=GSYMIN+GYSCAL
      ENDIF
C
C      IF(KID.EQ.'Z') THEN
C         CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
C     &               GXMIN,GXMAX,GYMIN,GYMAX)
C
C         CALL SETLIN(0,0,7)
C         CALL WFGBDY
C
C         CALL SETLIN(0,0,4)
C         CALL WFGPLA
C
C         CALL SETLIN(0,0,4)
C         CALL WFGANT
C
C         CALL OFFVEW
C      ENDIF
C
      CALL GMNMX2(GZ,NGXM,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
      CALL GQSCAL(GZMIN,GZMAX,GQZMIN,GQZMAX,GZSCAL)
      GZDEL=0.5*GZSCAL
      IF(GZDEL.EQ.0.0) GOTO 1000
      ISTEP=INT((GZMAX-GZMIN)/GZDEL)
C
      CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
C
      CALL SETLIN(0,0,7)
      CALL GFRAME
      CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.1,9)
      CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.1,9)
C      IF(NWMAX.LE.2) THEN
C         CALL GVALUE(GXORG,4*GXSCAL,0.0,0.0,NGULEN(GXSCAL))
C         CALL GVALUE(0.0,0.0,GYORG,4*GYSCAL,NGULEN(GYSCAL))
C      ENDIF
C
      CALL SETLIN(0,0,7)
      IF(GZMIN*GZMAX.GT.0.0) THEN
         GZORG=GQZMIN-GZDEL
      ELSE
         GZORG=0.5*GZDEL
      ENDIF
C
C         GZORG=GQZMIN
C         CALL SETLIN(0,0,6)
C         CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,
C     &               GZORG,GZDEL,ISTEP,0,0,KA)
C
         CALL SETLIN(0,0,6)
         CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,
     &               GZORG,GZDEL,ISTEP,0,0,KA)
         CALL SETLIN(0,0,5)
         CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,
     &               -GZORG,-GZDEL,ISTEP,0,2,KA)
C      ENDIF
C
 1000 CALL SETLIN(0,0,7)
      CALL SETCHS(0.2,0.0)
      GXPOS=0.5*(GPXMIN+GPXMAX-17*0.2)
      GYPOS=GPYMIN-0.4
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT(KWD(1:3),3)
         IF(KWD(4:4).EQ.'X') THEN
            CALL TEXT('(YZ) ',5)
         ELSE IF(KWD(4:4).EQ.'Y') THEN
            CALL TEXT('(XZ) ',5)
         ELSE IF(KWD(4:4).EQ.'Z') THEN
            CALL TEXT('(XY) ',5)
         ENDIF
         CALL TEXT(KWD(4:4),1)
         CALL TEXT('=',1)
         CALL TEXT(KWD(5:NCHM),NCHM-4)
      GYPOS=GYPOS-0.3
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT('MAX=',4)
         CALL NUMBR(GZMAX,'(1PE9.2)',9)
      GYPOS=GYPOS-0.3
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT('MIN=',4)
         CALL NUMBR(GZMIN,'(1PE9.2)',9)
      GYPOS=GYPOS-0.3
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT('STP=',4)
         CALL NUMBR(GZDEL,'(1PE9.2)',9)
      RETURN
      END
C
C     ****** PLOT 2D STRUCTURE OF FIELD AND POWER ******
C
      SUBROUTINE WFGPFC(NW,NWMAX,KWD)
C
      INCLUDE 'wfcomm.inc'
      PARAMETER(NSTEPM=101)
      PARAMETER(NRGBA=5,NRGBB=7)
C
      CHARACTER KWD*(NCHM),KID*1
      DIMENSION GAX(NGXM),GAY(NGYM)
      DIMENSION GZL(NSTEPM),GRGBL(3,0:NSTEPM)
      DIMENSION GRGBA(3,NRGBA),GLA(NRGBA)
      DIMENSION GRGBB(3,NRGBB),GLB(NRGBB)
      DATA GRGBA/0.0,0.0,1.0,
     &           0.0,1.0,1.0,
     &           1.0,1.0,1.0,
     &           1.0,1.0,0.0,
     &           1.0,0.0,0.0/
      DATA GLA/0.0,0.40,0.5,0.60,1.0/
      DATA GRGBB/1.0,1.0,1.0,
     &           0.0,1.0,1.0,
     &           0.0,0.0,1.0,
     &           0.0,1.0,0.0,
     &           1.0,0.0,0.0,
     &           1.0,1.0,0.0,
     &           1.0,1.0,1.0/
      DATA GLB/0.0,0.15,0.3,0.5,0.7,0.85,1.0/
C
      IF(NWMAX.LE.0) RETURN
C
      CALL WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
C
      KID=KWD(4:4)
      IF(KID.EQ.'X') THEN
         XMIN = YNDMIN
         XMAX = YNDMAX
         YMIN = ZNDMIN
         YMAX = ZNDMAX
      ELSE IF(KID.EQ.'Y') THEN
         XMIN = XNDMIN
         XMAX = XNDMAX
         YMIN = ZNDMIN
         YMAX = ZNDMAX
      ELSE IF(KID.EQ.'Z') THEN
         XMIN = XNDMIN
         XMAX = XNDMAX
         YMIN = YNDMIN
         YMAX = YNDMAX
      ELSE IF(KID.EQ.'A') THEN
         XMIN = 0.D0
         XMAX = 1.D0
         YMIN = 0.D0
         YMAX = 1.D0
      ENDIF
C
      DXLEN= XMAX-XMIN
      DYLEN= YMAX-YMIN
      DRATIO=DYLEN/DXLEN
      PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
      GXMIN=GUCLIP(XMIN)
      GXMAX=GUCLIP(XMAX)
      GYMIN=GUCLIP(YMIN)
      GYMAX=GUCLIP(YMAX)
      IF(DRATIO.GT.PRATIO) THEN
         GPYMIN=GUCLIP(PYMIN)+1.5
         GPYMAX=GUCLIP(PYMAX)
         PXMID=0.5D0*(PXMIN+PXMAX)
         PYLEN=PYMAX-PYMIN-1.5D0
         PXLEN=PYLEN/DRATIO
         GPXMIN=GUCLIP(PXMID-0.5D0*PXLEN)
         GPXMAX=GUCLIP(PXMID+0.5D0*PXLEN)
      ELSE
         GPXMIN=GUCLIP(PXMIN)+0.3
         GPXMAX=GUCLIP(PXMAX)-0.3
         PYMID=0.5D0*(PYMIN+PYMAX+1.5D0)
         PXLEN=PXMAX-PXMIN-0.6D0
         PYLEN=PXLEN*DRATIO
         GPYMIN=GUCLIP(PYMID-0.5D0*PYLEN)
         GPYMAX=GUCLIP(PYMID+0.5D0*PYLEN)
      ENDIF
C
      DX=(XMAX-XMIN)/(NGXMAX-1)
      DO NGX=1,NGXMAX
        GAX(NGX)=GUCLIP(XMIN+(NGX-1)*DX)
      ENDDO
      DY=(YMAX-YMIN)/(NGYMAX-1)
      DO NGY=1,NGYMAX
        GAY(NGY)=GUCLIP(YMIN+(NGY-1)*DY)
      ENDDO
C
      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSCAL)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSXMAX,GYSCAL)
      IF(GXMIN*GXMAX.LT.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GSXMIN+GXSCAL
      ENDIF
      IF(GYMIN*GYMAX.LT.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=GSYMIN+GYSCAL
      ENDIF
C
      ISTEP=50
C
      CALL GMNMX2(GZ,NGXM,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
         IF(GZMIN*GZMAX.LT.0.D0) THEN
            GZA=MAX(ABS(GZMAX),ABS(GZMIN))
            GDZ=2*GZA/ISTEP
            DO I=1,ISTEP
               GZL(I)=GDZ*(I-0.5)-GZA
            ENDDO
C      
            DO I=0,ISTEP
               GFACT=REAL(I)/REAL(ISTEP)
               CALL GUSRGB(GFACT,GRGBL(1,I),NRGBA,GLA,GRGBA)
            ENDDO
         ELSE
            GZA=GZMIN
            GDZ=(GZMAX-GZMIN)/ISTEP
            DO I=1,ISTEP
               GZL(I)=GDZ*I+GZA
            ENDDO
            DO I=0,ISTEP
               GFACT=REAL(I)/REAL(ISTEP)
               CALL GUSRGB(GFACT,GRGBL(1,I),NRGBB,GLB,GRGBB)
            ENDDO
         ENDIF
C
      CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
C
         CALL CONTF2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,GZL,GRGBL,ISTEP,0)
         CALL SETLIN(0,0,7)
         CALL GFRAME
         CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.1,9)
         CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.1,9)
C
         CALL RGBBAR(GPXMAX+0.2,GPXMAX+0.5,GPYMIN,GPYMAX,
     &               GRGBL,ISTEP+1,1)
C
C      IF(KID.EQ.'Z') THEN
C         CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
C     &               GXMIN,GXMAX,GYMIN,GYMAX)
C
C         CALL SETLIN(0,0,7)
C         CALL WFGBDY
C
C         CALL SETLIN(0,0,4)
C         CALL WFGPLA
C
C         CALL SETLIN(0,0,4)
C         CALL WFGANT
C
C         CALL OFFVEW
C      ENDIF
C
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.2,0.0)
      GXPOS=0.5*(GPXMIN+GPXMAX-17*0.2)
      GYPOS=GPYMIN-0.4
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT(KWD(1:3),3)
         IF(KWD(4:4).EQ.'X') THEN
            CALL TEXT('(YZ) ',5)
         ELSE IF(KWD(4:4).EQ.'Y') THEN
            CALL TEXT('(XZ) ',5)
         ELSE IF(KWD(4:4).EQ.'Z') THEN
            CALL TEXT('(XY) ',5)
         ENDIF
         CALL TEXT(KWD(4:4),1)
         CALL TEXT('=',1)
         CALL TEXT(KWD(5:NCHM),NCHM-4)
      GYPOS=GYPOS-0.3
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT('MAX=',4)
         CALL NUMBR(GZMAX,'(1PE9.2)',9)
      GYPOS=GYPOS-0.3
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT('MIN=',4)
         CALL NUMBR(GZMIN,'(1PE9.2)',9)
      GYPOS=GYPOS-0.3
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT('TOP=',4)
         CALL NUMBR(GZA,'(1PE9.2)',9)
      RETURN
      END
C
C     ****** PLOT 2D STRUCTURE OF FIELD AND POWER ******
C
      SUBROUTINE WFGPBC(NW,NWMAX,KWD)
C
      INCLUDE 'wfcomm.inc'
      DIMENSION KA(8,NGXM,NGYM)
      EXTERNAL R2W2B
C
      CHARACTER KWD*(NCHM),KID*1
      DIMENSION GAX(NGXM),GAY(NGYM)
C
      IF(NWMAX.LE.0) RETURN
C
      CALL WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
C
      KID=KWD(4:4)
      IF(KID.EQ.'X') THEN
         XMIN = YNDMIN
         XMAX = YNDMAX
         YMIN = ZNDMIN
         YMAX = ZNDMAX
      ELSE IF(KID.EQ.'Y') THEN
         XMIN = XNDMIN
         XMAX = XNDMAX
         YMIN = ZNDMIN
         YMAX = ZNDMAX
      ELSE IF(KID.EQ.'Z') THEN
         XMIN = XNDMIN
         XMAX = XNDMAX
         YMIN = YNDMIN
         YMAX = YNDMAX
      ELSE IF(KID.EQ.'A') THEN
         XMIN = 0.D0
         XMAX = 1.D0
         YMIN = 0.D0
         YMAX = 1.D0
      ENDIF
C
      DXLEN= XMAX-XMIN
      DYLEN= YMAX-YMIN
      DRATIO=DYLEN/DXLEN
      PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
      GXMIN=GUCLIP(XMIN)
      GXMAX=GUCLIP(XMAX)
      GYMIN=GUCLIP(YMIN)
      GYMAX=GUCLIP(YMAX)
      IF(DRATIO.GT.PRATIO) THEN
         GPYMIN=GUCLIP(PYMIN)+1.5
         GPYMAX=GUCLIP(PYMAX)
         PXMID=0.5D0*(PXMIN+PXMAX)
         PYLEN=PYMAX-PYMIN-1.5D0
         PXLEN=PYLEN/DRATIO
         GPXMIN=GUCLIP(PXMID-0.5D0*PXLEN)
         GPXMAX=GUCLIP(PXMID+0.5D0*PXLEN)
      ELSE
         GPXMIN=GUCLIP(PXMIN)+0.3
         GPXMAX=GUCLIP(PXMAX)-0.3
         PYMID=0.5D0*(PYMIN+PYMAX+1.5D0)
         PXLEN=PXMAX-PXMIN-0.6D0
         PYLEN=PXLEN*DRATIO
         GPYMIN=GUCLIP(PYMID-0.5D0*PYLEN)
         GPYMAX=GUCLIP(PYMID+0.5D0*PYLEN)
      ENDIF
C
      DX=(XMAX-XMIN)/(NGXMAX-1)
      DO NGX=1,NGXMAX
        GAX(NGX)=GUCLIP(XMIN+(NGX-1)*DX)
      ENDDO
      DY=(YMAX-YMIN)/(NGYMAX-1)
      DO NGY=1,NGYMAX
        GAY(NGY)=GUCLIP(YMIN+(NGY-1)*DY)
      ENDDO
C
      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSCAL)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSXMAX,GYSCAL)
      IF(GXMIN*GXMAX.LT.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=(GSXMIN/(2*GXSCAL)+1)*2*GXSCAL
      ENDIF
      IF(GYMIN*GYMAX.LT.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=(GSYMIN/(2*GYSCAL)+1)*2*GYSCAL
      ENDIF
C
      CALL GMNMX2(GZ,NGXM,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
      GZMAX=MAX(ABS(GZMIN),ABS(GZMAX))
      GZMIN=-GZMAX
      CALL GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GZSCAL)
      GZORG=0.0
C
      CALL SETLIN(0,0,7)
      GXL=(GPXMAX-GPXMIN)
      GYL=GXL
      GZL=0.5*(GPYMAX-GPYMIN)
      IF(NGRAPH.EQ.3) THEN
         GPHI=-80.0
      ELSEIF(NGRAPH.EQ.4) THEN
         GPHI=190.0
      ELSEIF(NGRAPH.EQ.5) THEN
         GPHI=100.0
      ELSEIF(NGRAPH.EQ.6) THEN
         GPHI=10.0
      ENDIF
      GTHETA=40.0
      GRADIUS=100.0
C
      CALL GDEFIN3D(GPXMIN,GPXMAX,GPYMIN,GPYMAX,GXL,GYL,GZL)
      CALL GVIEW3D(GPHI,GTHETA,GRADIUS,0.9,1,
     &             0.5*(GXMIN+GXMAX),0.5*(GYMIN+GYMAX),0.0)
      CALL GDATA3D1(GZ,NGXM,NGXMAX,NGYMAX,
     &              GXMIN,GXMAX,GYMIN,GYMAX,GZMIN,GZMAX)
C
      CALL CPLOT3D1(1,R2W2B)
      CALL CONTQ3D1(GZMIN,0.1*(GZMAX-GZMIN),11,0,0,KA,R2W2B,0)
C
      CALL GAXIS3D(0)
      CALL GSCALE3DX(GXORG,GXSCAL,0.2,2)
      CALL GSCALE3DY(GYORG,GYSCAL,0.2,2)
      CALL GSCALE3DZ(GZORG,GZSCAL,0.2,2)
      CALL SETCHS(0.2,0.0)
      CALL GVALUE3DX(GXORG,2.*GXSCAL,-6,2)
      CALL GVALUE3DY(GYORG,2.*GYSCAL,-6,2)
      CALL GVALUE3DZ(GZORG,2.*GZSCAL,-2,0)
C
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.2,0.0)
      GXPOS=0.5*(GPXMIN+GPXMAX-17*0.2)
      GYPOS=GPYMIN-0.4
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT(KWD(1:3),3)
         IF(KWD(4:4).EQ.'X') THEN
            CALL TEXT('(YZ) ',5)
         ELSE IF(KWD(4:4).EQ.'Y') THEN
            CALL TEXT('(XZ) ',5)
         ELSE IF(KWD(4:4).EQ.'Z') THEN
            CALL TEXT('(XY) ',5)
         ENDIF
         CALL TEXT(KWD(4:4),1)
         CALL TEXT('=',1)
         CALL TEXT(KWD(5:NCHM),NCHM-4)
      GYPOS=GYPOS-0.3
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT('MAX=',4)
         CALL NUMBR(GZMAX,'(1PE9.2)',9)
      GYPOS=GYPOS-0.3
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT('MIN=',4)
         CALL NUMBR(GZMIN,'(1PE9.2)',9)
      GYPOS=GYPOS-0.3
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT('TOP=',4)
         CALL NUMBR(GZA,'(1PE9.2)',9)
      RETURN
      END
C
C     ****** DRAW 1D PROFILE ******
C
      SUBROUTINE WFGPFR(NW,NWMAX,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION IPAT(3)
      CHARACTER KWD*(NCHM)
      DATA IPAT/0,2,4/
C
      IF(NWMAX.LE.0) RETURN
C
      CALL WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
C
C     ------ DO NOT ALLOW LANDSCAPE GRAPH ------
C
      PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-2.5D0)
      IF(PRATIO.LT.1.D0) THEN
         PYLEN=PYMAX-PYMIN-1.5D0
         PXLEN=PYLEN+2.5D0
         PXMID=0.5D0*(PXMIN+PXMAX)
         PXMIN=PXMID-0.5D0*PXLEN
         PXMAX=PXMID+0.5D0*PXLEN
      ENDIF
C      
      GPXMIN=GUCLIP(PXMIN)+2.2
      GPXMAX=GUCLIP(PXMAX)-0.3
      GPYMIN=GUCLIP(PYMIN)+1.5
      GPYMAX=GUCLIP(PYMAX)
C
      IF(KWD(1:1).EQ.'E'.OR.
     &   KWD(1:1).EQ.'D'.OR.
     &   KWD(1:1).EQ.'B'.OR.
     &   KWD(1:1).EQ.'A') THEN
         NGMAX=3
      ELSE
         NGMAX=1
      ENDIF
C
      CALL GMNMX1(GX,1,NGVMAX,1,GXMIN1,GXMAX1)
      CALL GMNMX2(GV,NGVM,1,NGVMAX,1,1,NGMAX,1,GYMIN1,GYMAX1)
      CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP)
      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      IF(GXMIN*GXMAX.LT.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GXMIN
      ENDIF
      IF(GYMIN*GYMAX.LT.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=GYMIN
      ENDIF
C
      CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
C
      CALL SETLIN(0,0,4)
      CALL SETCHS(0.25,0.0)
      CALL GFRAME
      CALL GSCALE(GXORG,GXSTEP,0.0,0.0,0.2,9)
      CALL GSCALE(0.0,0.0,GYORG,GYSTEP,0.2,9)
      IF(GYMIN*GYMAX.LT.0.0) THEN
         CALL GSCALE(0.0,0.0,GYORG,100*GYSTEP,0.0,0)
      ENDIF
      CALL GVALUE(GXORG,4*GXSTEP,0.0,0.0,NGULEN(4*GXSTEP))
      CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))
C
      DO NG=1,NGMAX
         CALL SETLIN(0,0,8-NG)
         CALL GPLOTP(GX,GV(1,NG),1,NGVMAX,1,0,0,IPAT(NG))
      ENDDO
C
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.2,0.0)
      GXPOS=0.5*(GPXMIN+GPXMAX-20*0.2)
      GYPOS=GPYMIN-1.0
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT(KWD(1:2),2)
         IF(KWD(3:3).EQ.'X') THEN
            CALL TEXT('(X): Y,Z=',9)
         ELSE IF(KWD(3:3).EQ.'Y') THEN
            CALL TEXT('(Y): X,Z=',9)
         ELSE IF(KWD(3:3).EQ.'Z') THEN
            CALL TEXT('(Z): X,Y=',9)
         ENDIF
         CALL TEXT(KWD(4:NCHM),NCHM-3)
      RETURN
      END
C
C     ****** DRAW PARAMETER ON GRAPHIC SCREEN ******
C
      SUBROUTINE WFGPRM
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION REST(NAM),REAT(NAM)
      DIMENSION SRFR(NMDM,NBM),SRFI(NMDM,NBM),SRFL(NMDM,NBM)
C
      DO NA=1,NAMAX
         REST(NA)=DBLE(CIMP(NA))
         REAT(NA)=IMAG(CIMP(NA))
      ENDDO
      DO NB=1,NBMAX
         DO L=1,NMBDY(NB)
            SRFR(L,NB)=DBLE(CRFL(L,NB))
            SRFI(L,NB)=IMAG(CRFL(L,NB))
            SRFL(L,NB)=ABS(CRFL(L,NB))**2
         ENDDO
      ENDDO
C
      GXMIN=0.
      GYMAX=18.2
      GRCHH=0.30
      CALL SETCHS(GRCHH,0.)
      GDX=15.*GRCHH
      GDY=-1.5*GRCHH
      GXL=GXMIN
      GYL=GYMAX
C
      GYL=GYL+GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('RF  =',5)
      IF(RF.GE.1.D5) THEN
         CALL NUMBD(RF  ,'(F7.0)',7)
      ELSEIF(RF.GE.1.D4) THEN
         CALL NUMBD(RF  ,'(F7.1)',7)
      ELSEIF(RF.GE.1.D3) THEN
         CALL NUMBD(RF  ,'(F7.2)',7)
      ELSE
         CALL NUMBD(RF  ,'(F7.3)',7)
      ENDIF
C
      GXL=GXL+GDX
      CALL MOVE(GXL,GYL)
      CALL TEXT('BB  =',5)
      CALL NUMBD(BB  ,'(F7.3)',7)
C
      GXL=GXL+GDX
      CALL MOVE(GXL,GYL)
      CALL TEXT('M=',2)
      CALL NUMBI(MODELG,'(I2)',2)
      CALL NUMBI(MODELB,'(I2)',2)
      CALL NUMBI(MODELD,'(I2)',2)
      CALL NUMBI(MODELP,'(I2)',2)
      CALL NUMBI(MODELS,'(I2)',2)
      CALL NUMBI(MODELN,'(I2)',2)
C
      GXL=GXMIN
      GYL=GYL+GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('NNMAX=',6)
      CALL NUMBI(NNMAX,'(I6)',6)
C
      GXL=GXL+GDX
      CALL MOVE(GXL,GYL)
      CALL TEXT('XYZ MAX=',8)
      CALL NUMBD(XNDMAX,'(F7.3)',7)
      CALL NUMBD(YNDMAX,'(F7.3)',7)
      CALL NUMBD(ZNDMAX,'(F7.3)',7)
C
      GXL=GXMIN
      GYL=GYL+GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('NEMAX=',6)
      CALL NUMBI(NEMAX,'(I6)',6)
C
      GXL=GXL+GDX
      CALL MOVE(GXL,GYL)
      CALL TEXT('XYZ MIN=',8)
      CALL NUMBD(XNDMIN,'(F7.3)',7)
      CALL NUMBD(YNDMIN,'(F7.3)',7)
      CALL NUMBD(ZNDMIN,'(F7.3)',7)
C
C      GXL=GXL+GDX
C      CALL MOVE(GXL,GYL)
C      CALL TEXT('P=',2)
C      CALL NUMBD(TSPWR,'(1PE10.3)',10)
C
      GXL=GXMIN
      GYL=GYL+GDY
C
      CALL MOVE(GXL,GYL)
      CALL TEXT(' NK',3)
      CALL TEXT(' NM',3)
      CALL TEXT(' PABS     ',10)
      CALL TEXT(' NK',3)
      CALL TEXT(' NM',3)
      CALL TEXT(' PABS     ',10)
C     
      I=0
      DO NK=1,NKMAX
         NM=NMKA(NK)
         IF(MOD(I,2).EQ.0) THEN
            GXL=GXMIN
            GYL=GYL+GDY
         ELSE
            GXL=GXMIN+16.*GRCHH
         ENDIF
         I=I+1
         CALL MOVE(GXL,GYL)
         CALL NUMBI(NK,'(I3)',3)
         CALL NUMBI(NM,'(I3)',3)
         CALL NUMBD(PABSK(NK),'(1PE10.2)',10)
      ENDDO
C
      IF(NSMAX.GT.0) THEN
         GXL=GXMIN
         GYL=GYL+GDY
         CALL MOVE(GXL,GYL)
         CALL TEXT(' NS',3)
         CALL TEXT('    PA    ',10)
         CALL TEXT(' PZ ',4)
         CALL TEXT('    PN    ',10)
         CALL TEXT(' PZCL ',6)
         CALL TEXT('   PABS   ',10)
C     
         DO NS=1,NSMAX
            GXL=GXMIN
            GYL=GYL+GDY
            CALL MOVE(GXL,GYL)
            CALL NUMBI(NS,'(I3)',3)
            CALL NUMBD(PA(NS),'(1PE10.2)',10)
            CALL NUMBD(PZ(NS),'(F4.0)',4)
            CALL NUMBD(PN(NS),'(1PE10.2)',10)
            CALL NUMBD(PZCL(NS),'(F6.2)',6)
            CALL NUMBD(PABSS(NS),'(1PE10.2)',10)
         ENDDO
      ENDIF
C
      GXL=GXMIN+45.*GRCHH
      GYL=GYMAX+GDY
      IF(NAMAX.GT.0) THEN
         CALL MOVE(GXL,GYL)
         CALL TEXT('IJ',2)
         CALL TEXT('  AJ ',5)
         CALL TEXT(' PHASE ',7)
         CALL TEXT('     R     ',11)
         CALL TEXT('     X     ',11)
C
         DO NA=1,NAMAX
            GXL=GXMIN+45.*GRCHH
            GYL=GYL+GDY
            CALL MOVE(GXL,GYL)
            CALL NUMBI(NA,'(I2)',2)
            CALL NUMBD(AJ(NA),'(F5.1)',5)
            CALL NUMBD(APH(NA),'(F7.1)',7)
            CALL NUMBD(REST(NA),'(1PE11.3)',11)
            CALL NUMBD(REAT(NA),'(1PE11.3)',11)
         ENDDO
      ELSEIF(NBMAX.GT.0) THEN
         CALL MOVE(GXL,GYL)
         CALL TEXT('NB',2)
         CALL TEXT(' MD',3)
         CALL TEXT('  Real   ',10)
         CALL TEXT('  Imag   ',10)
         CALL TEXT('  REFL     ',11)
C
         DO NB=1,NBMAX
            IF(KABDY(NB).GE.8) THEN
               DO L=1,NMBDY(NB)
                  GXL=GXMIN+45.*GRCHH
                  GYL=GYL+GDY
                  CALL MOVE(GXL,GYL)
                  CALL NUMBI(NB,'(I2)',2)
                  CALL NUMBI(L,'(I3)',3)
                  CALL NUMBD(SRFR(L,NB),'(1PE10.2)',10)
                  CALL NUMBD(SRFI(L,NB),'(1PE10.2)',10)
                  CALL NUMBD(SRFL(L,NB),'(1PE11.3)',11)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
C
C     ****** Draw Vessel Boundary ******
C
      SUBROUTINE WFGBDY
C
      INCLUDE 'wfcomm.inc'
C
      RETURN
      END
C
C     ****** Draw Plasma Boundary ******
C
      SUBROUTINE WFGPLA
C
      INCLUDE 'wfcomm.inc'
C
      NPMAX=100
      DTH=2.D0*PI/NPMAX
      GXL=GUCLIP(RA)
      GYL=0.0
      CALL MOVE(GXL,GYL)
      DO 10 I=1,NPMAX
         THETA=DTH*I
         GXL=GUCLIP(RA*COS(THETA))
         GYL=GUCLIP(RA*SIN(THETA))
         CALL DRAW(GXL,GYL)
   10 CONTINUE
      RETURN
      END
C
C     ****** Draw Antenna Path ******
C
      SUBROUTINE WFGANT
C
      INCLUDE 'wfcomm.inc'
C
      IF(NDRAWA.EQ.0) THEN
         DO 100 NA=1,NAMAX
            GXL=GUCLIP(XJ0(1,NA))
            GYL=GUCLIP(YJ0(1,NA))
            CALL MOVE(GXL,GYL)
            DO 90 I=2,JNUM0(NA)
               GXL=GUCLIP(XJ0(I,NA))
               GYL=GUCLIP(YJ0(I,NA))
               CALL DRAW(GXL,GYL)
   90       CONTINUE
  100    CONTINUE
      ELSE
         DO 200 NA=1,NAMAX
            GXL=GUCLIP(XJ(1,NA))
            GYL=GUCLIP(YJ(1,NA))
            CALL MOVE(GXL,GYL)
            DO 190 I=2,JNUM(NA)
               GXL=GUCLIP(XJ(I,NA))
               GYL=GUCLIP(YJ(I,NA))
               CALL DRAW(GXL,GYL)
  190       CONTINUE
  200    CONTINUE
      ENDIF
      RETURN
      END
C
C     ****** Draw Element Data ******
C
      SUBROUTINE WFGELM
C
      INCLUDE 'wfcomm.inc'
      DIMENSION I1(4),I2(4),I3(4)
      DATA I1,I2,I3/1,2,3,4, 2,3,4,1, 3,4,1,2/
      DATA EPS/1.D-6/
C
      CALL SETCHR(0.2,0.15,0.2,0.,-30.)
C
      DO IE=1,NEMAX
         DO K=1,4
         IN1=NDELM(I1(K),IE)
         IN2=NDELM(I2(K),IE)
         IN3=NDELM(I3(K),IE)
         IF(ZND(IN1).GT.ZNDMAX-EPS.AND.
     &      ZND(IN2).GT.ZNDMAX-EPS.AND.
     &      ZND(IN3).GT.ZNDMAX-EPS) THEN
            GX1=GUCLIP(XND(IN1))
            GY1=GUCLIP(YND(IN1))
            GX2=GUCLIP(XND(IN2))
            GY2=GUCLIP(YND(IN2))
            GX3=GUCLIP(XND(IN3))
            GY3=GUCLIP(YND(IN3))
            IF (GX1.GT.GX2) THEN
               CALL MOVE(GX2,GY2)
               CALL DRAW(GX1,GY1)
            ELSE
               CALL MOVE(GX1,GY1)
               CALL DRAW(GX2,GY2)
            END IF
            IF (GX2.GT.GX3) THEN
               CALL MOVE(GX3,GY3)
               CALL DRAW(GX2,GY2)
            ELSE
               CALL MOVE(GX2,GY2)
               CALL DRAW(GX3,GY3)
            END IF
            IF (GX3.GT.GX1) THEN
               CALL MOVE(GX1,GY1)
               CALL DRAW(GX3,GY3)
            ELSE
               CALL MOVE(GX3,GY3)
               CALL DRAW(GX1,GY1)
            END IF
C
            IF(NDRAWD.GE.2) THEN
               GXC=(GX1+GX2+GX3)/3.
               GYC=(GY1+GY2+GY3)/3.
               IEL=IE
               CALL GNUMBI(GXC,GYC,IEL,2)
            ENDIF
         ENDIF
         ENDDO
      ENDDO
C
      IF(NDRAWD.GE.3) THEN
         CALL SETCHS(0.2,0.)
         DO 200 IN=1,NNMAX
            GX1=GUCLIP(XND(IN))
            GY1=GUCLIP(YND(IN))
            INL=IN
            CALL GNUMBI(GX1,GY1,INL,0)
  200    CONTINUE
      ENDIF
C
      RETURN
      END
C
C     ****** Draw Element Data ******
C
      SUBROUTINE WFGNAS(ID)
C
      INCLUDE 'wfcomm.inc'
C
      CALL WFGWIN(1,1,PXMIN,PXMAX,PYMIN,PYMAX)
      GPXMIN=GUCLIP(PXMIN)
      GPXMAX=GUCLIP(PXMAX)
      GPYMIN=GUCLIP(PYMIN)+1.0
      GPYMAX=GUCLIP(PYMAX)+1.0
C
      DXLEN= XNDMAX-XNDMIN+0.5D0*(YNDMAX-YNDMIN)
      DYLEN= ZNDMAX-ZNDMIN+0.5D0*(YNDMAX-YNDMIN)
      DRATIO=DYLEN/DXLEN
      PRATIO=(PYMAX-PYMIN)/(PXMAX-PXMIN)
      IF(DRATIO.GT.PRATIO) THEN
         GYMIN=GUCLIP(ZNDMIN+0.5D0*YNDMIN)
         GYMAX=GUCLIP(ZNDMAX+0.5D0*YNDMAX)
         XMID=0.5D0*(XNDMIN+XNDMAX+0.5D0*(YNDMAX+YNDMIN))
         XLEN=DYLEN/PRATIO
         GXMIN=GUCLIP(XMID-0.5D0*XLEN)
         GXMAX=GUCLIP(XMID+0.5D0*XLEN)
      ELSE
         GXMIN=GUCLIP(XNDMIN+0.5D0*YNDMIN)
         GXMAX=GUCLIP(XNDMAX+0.5D0*YNDMIN)
         YMID=0.5D0*(ZNDMIN+ZNDMAX+0.5D0*(YNDMAX+YNDMIN))
         YLEN=DXLEN*PRATIO
         GYMIN=GUCLIP(YMID-0.5D0*YLEN)
         GYMAX=GUCLIP(YMID+0.5D0*YLEN)
      ENDIF
C
      CALL PAGES
      CALL WFPRME
      CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
C     &            GXMIN,GXMAX,GYMIN,GYMAX)
C     &     -0.03,-0.01,-0.005,0.005)
C     &     -0.03,0.03,-0.02,0.02)
C
      CALL SETLIN(0,0,7)
      CALL WFGELM3(ID)
C
      CALL PAGEE
      RETURN
      END
C
C     ****** Draw 3D Element Data ******
C
      SUBROUTINE WFGELM3(ID)
C
      INCLUDE 'wfcomm.inc'
C
C     ***** DRAW EDGE *****
C
      IF(ID.EQ.0) THEN
         DO IE=1,NEMAX
            IN1=NDELM(1,IE)
            IN2=NDELM(2,IE)
            IN3=NDELM(3,IE)
            IN4=NDELM(4,IE)
            GX1=GUCLIP(XND(IN1))
            GY1=GUCLIP(ZND(IN1)+0.5D0*YND(IN1))
            GX2=GUCLIP(XND(IN2))
            GY2=GUCLIP(ZND(IN2)+0.5D0*YND(IN2))
            GX3=GUCLIP(XND(IN3))
            GY3=GUCLIP(ZND(IN3)+0.5D0*YND(IN3))
            GX4=GUCLIP(XND(IN4))
            GY4=GUCLIP(ZND(IN4)+0.5D0*YND(IN4))
            CALL MOVE(GX1,GY1)
            CALL DRAW(GX2,GY2)
            CALL DRAW(GX3,GY3)
            CALL DRAW(GX4,GY4)
            CALL DRAW(GX2,GY2)
            CALL MOVE(GX4,GY4)
            CALL DRAW(GX1,GY1)
            CALL DRAW(GX3,GY3)
C     
            IF(NDRAWD.GE.2) THEN
               GXC=(GX1+GX2+GX3)/3.
               GYC=(GY1+GY2+GY3)/3.
               IEL=IE
               CALL GNUMBI(GXC,GYC,IEL,2)
            ENDIF
         ENDDO
C
         IF(NDRAWD.GE.3) THEN
            CALL SETCHS(0.2,0.)
            DO IN=1,NNMAX
               GX1=GUCLIP(XND(IN))
               GY1=GUCLIP(YND(IN))
               INL=IN
               CALL GNUMBI(GX1,GY1,INL,0)
            ENDDo
         ENDIF
C
C     ***** DRAW SURFACE EDGE *****
C
      ELSEIF(ID.EQ.1) THEN
         CALL SETRGB(0.0,0.0,1.0)
         DO NSF=1,NSFMAX
            ND1=NDSRF(1,NSF)
            ND2=NDSRF(2,NSF)
            ND3=NDSRF(3,NSF)
            GX1=GUCLIP(XND(ND1))
            GY1=GUCLIP(ZND(ND1)+0.5D0*YND(ND1))
            GX2=GUCLIP(XND(ND2))
            GY2=GUCLIP(ZND(ND2)+0.5D0*YND(ND2))
            GX3=GUCLIP(XND(ND3))
            GY3=GUCLIP(ZND(ND3)+0.5D0*YND(ND3))
            CALl MOVE(GX1,GY1)
            CALl DRAW(GX2,GY2)
            CALl DRAW(GX3,GY3)
            CALl DRAW(GX1,GY1)
         ENDDO
C
C     ***** DRAW NODES *****
C
      ELSEIF(ID.EQ.2) THEN
         CALL SETMKS(1,0.5)
         DO NN=1,NNMAX
            GX1=GUCLIP(XND(NN))
            GY1=GUCLIP(ZND(NN)+0.5D0*YND(NN))
            KA=KANOD(NN)
            IF(KA.EQ.0) THEN
               CALL SETRGB(0.0,1.0,0.0)
            ELSEIF(KA.EQ.1) THEN
               CALL SETRGB(0.0,0.0,1.0)
            ELSEIF(KA.EQ.2) THEN
               CALL SETRGB(1.0,0.0,0.0)
            ELSE
               CALL SETRGB(0.0,1.0,0.0)
            ENDIF
            CALL MARK(GX1,GY1)
         ENDDO
C
      ELSEIF(ID.EQ.3) THEN
         CALL SETMKS(1,0.5)
         DO NN=1,NNMAX
            GX1=GUCLIP(XND(NN))
            GY1=GUCLIP(ZND(NN)+0.5D0*YND(NN))
            KA=KANOD(NN)
            IF(KA.GE.1.OR.KA.LT.0) THEN
               IF(KA.EQ.1) THEN
                  CALL SETRGB(0.0,0.0,1.0)
               ELSEIF(KA.EQ.2) THEN
                  CALL SETRGB(1.0,0.0,0.0)
               ELSE
                  CALL SETRGB(0.0,1.0,0.0)
               ENDIF
               CALL MARK(GX1,GY1)
            ENDIF
         ENDDO
C
      ELSEIF(ID.EQ.4) THEN
         CALL SETMKS(1,0.5)
         DO NN=1,NNMAX
            GX1=GUCLIP(XND(NN))
            GY1=GUCLIP(ZND(NN)+0.5D0*YND(NN))
            KA=KANOD(NN)
            IF(KA.GE.2.OR.KA.LT.0) THEN
               IF(KA.EQ.2) THEN
                  CALL SETRGB(1.0,0.0,0.0)
               ELSEIF(KA.EQ.-4) THEN
                  CALL SETRGB(1.0,0.0,1.0)
               ELSE
                  CALL SETRGB(0.0,1.0,0.0)
               ENDIF
               CALL MARK(GX1,GY1)
            ENDIF
         ENDDO
C
      ELSEIF(ID.EQ.5) THEN
         CALL SETMKS(1,0.5)
         DO NN=1,NNMAX
            GX1=GUCLIP(XND(NN))
            GY1=GUCLIP(ZND(NN)+0.5D0*YND(NN))
            KA=KANOD(NN)
            IF(KA.EQ.2) THEN
               CALL SETRGB(1.0,0.0,0.0)
               CALL MARK(GX1,GY1)
            ENDIF
         ENDDO
      ENDIF
C
      RETURN
      END
C
C     ****** Draw Element Data ******
C
      SUBROUTINE WFGDIV
C
      INCLUDE 'wfcomm.inc'
C
      CALL WFGWIN(1,1,PXMIN,PXMAX,PYMIN,PYMAX)
      GPXMIN=GUCLIP(PXMIN)
      GPXMAX=GUCLIP(PXMAX)
      GPYMIN=GUCLIP(PYMIN)+1.0
      GPYMAX=GUCLIP(PYMAX)+1.0
C
      DXLEN= XNDMAX-XNDMIN
      DYLEN= YNDMAX-YNDMIN
      DRATIO=DYLEN/DXLEN
      PRATIO=(PYMAX-PYMIN)/(PXMAX-PXMIN)
      IF(DRATIO.GT.PRATIO) THEN
         GYMIN=GUCLIP(YNDMIN)
         GYMAX=GUCLIP(YNDMAX)
         XMID=0.5D0*(XNDMIN+XNDMAX)
         XLEN=DYLEN/PRATIO
         GXMIN=GUCLIP(XMID-0.5D0*XLEN)
         GXMAX=GUCLIP(XMID+0.5D0*XLEN)
      ELSE
         GXMIN=GUCLIP(XNDMIN)
         GXMAX=GUCLIP(XNDMAX)
         YMID=0.5D0*(YNDMIN+YNDMAX)
         YLEN=DXLEN*PRATIO
         GYMIN=GUCLIP(YMID-0.5D0*YLEN)
         GYMAX=GUCLIP(YMID+0.5D0*YLEN)
      ENDIF
C
      CALL PAGES
      CALL WFPRME
      CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
C
      CALL SETLIN(0,0,7)
      IF(NDRAWD.EQ.0) THEN
         CALL WFGBDY
      ELSE
         CALL WFGELM
      ENDIF
C
      CALL PAGEE
      RETURN
      END
C
C     ****** Draw Element Paramters ******
C
      SUBROUTINE WFPRME
C
      INCLUDE 'wfcomm.inc'
C
      GXMIN=20.
      GYMAX=17.
      CALL SETCHS(0.3,0.)
      GDY=1.5*0.3
C
      GXL=GXMIN
      GYL=GYMAX
      CALL MOVE(GXL,GYL)
      CALL TEXT('NNMAX=',6)
      CALL NUMBI(NNMAX,'(I5)',5)
C
      GYL=GYL-GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('NEMAX=',6)
      CALL NUMBI(NEMAX,'(I5)',5)
C
      GYL=GYL-GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('NBMAX=',6)
      CALL NUMBI(NBMAX,'(I5)',5)
C
      GYL=GYL-GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('MBND =',6)
      CALL NUMBI(MBND,'(I5)',5)
C
      GYL=GYL-GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('MLEN =',6)
      CALL NUMBI(MLEN,'(I5)',5)
      RETURN
      END
C
C     ****** Draw Antenna ******
C
      SUBROUTINE WFPLTA
C
      INCLUDE 'wfcomm.inc'
C
      CALL WFGWIN(1,1,PXMIN,PXMAX,PYMIN,PYMAX)
      GPXMIN=GUCLIP(PXMIN)
      GPXMAX=GUCLIP(PXMAX)
      GPYMIN=GUCLIP(PYMIN)+1.0
      GPYMAX=GUCLIP(PYMAX)+1.0
C
      DXLEN= XNDMAX-XNDMIN
      DYLEN= YNDMAX-YNDMIN
      DRATIO=DYLEN/DXLEN
      PRATIO=(PYMAX-PYMIN)/(PXMAX-PXMIN)
      IF(DRATIO.GT.PRATIO) THEN
         GYMIN=GUCLIP(YNDMIN)
         GYMAX=GUCLIP(YNDMAX)
         XMID=0.5D0*(XNDMIN+XNDMAX)
         XLEN=DYLEN/PRATIO
         GXMIN=GUCLIP(XMID-0.5D0*XLEN)
         GXMAX=GUCLIP(XMID+0.5D0*XLEN)
      ELSE
         GXMIN=GUCLIP(XNDMIN)
         GXMAX=GUCLIP(XNDMAX)
         YMID=0.5D0*(YNDMIN+YNDMAX)
         YLEN=DXLEN*PRATIO
         GYMIN=GUCLIP(YMID-0.5D0*YLEN)
         GYMAX=GUCLIP(YMID+0.5D0*YLEN)
      ENDIF
C
      CALL PAGES
      CALL WFPRMJ
      CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
C
      CALL SETLIN(0,0,4)
      IF(NDRAWA.LE.1) THEN
         CALL WFGBDY
      ELSE
         NTEMP=NDRAWD
         NDRAWD=NDRAWA-1
         CALL WFGELM
         NDRAWD=NTEMP
      ENDIF
C
      CALL SETLIN(0,0,5)
      CALL WFGPLA
C
      CALL SETLIN(0,0,6)
      CALL WFGANT
C
      CALL PAGEE
      RETURN
      END
C
C     ****** Draw Antenna Paramters ******
C
      SUBROUTINE WFPRMJ
C
      INCLUDE 'wfcomm.inc'
C
      GXMIN=20.
      GYMAX=17.
      CALL SETCHS(0.3,0.)
      GDY=1.5*0.3
C
      GXL=GXMIN
      GYL=GYMAX
      CALL MOVE(GXL,GYL)
      CALL TEXT('NNMAX=',6)
      CALL NUMBI(NNMAX,'(I5)',5)
C
      GYL=GYL-GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('NEMAX=',6)
      CALL NUMBI(NEMAX,'(I5)',5)
C
      GYL=GYL-GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('NBMAX=',6)
      CALL NUMBI(NBMAX,'(I5)',5)
C
      GYL=GYL-GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('MBND =',6)
      CALL NUMBI(MBND,'(I5)',5)
C
      GYL=GYL-GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('MLEN =',6)
      CALL NUMBI(MLEN,'(I5)',5)
C
      GYL=GYL-GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('JNUM =',6)
      GXL=GXL+6*0.3
      DO 10 NA=1,NAMAX
         CALL MOVE(GXL,GYL)
         CALL NUMBI(JNUM(NA),'(I5)',5)
         GYL=GYL-GDY
   10 CONTINUE
      RETURN
      END
