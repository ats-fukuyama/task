C     $Id$
C
C     ****** CALCULATE RANGE OF WINDOW ******
C
      SUBROUTINE WFGINI
C
      INCLUDE 'wfcomm.inc'
C
      XDMIN=XD(1)
      XDMAX=XD(1)
      YDMIN=YD(1)
      YDMAX=YD(1)
C
      DO IN=2,NNOD
         XDMIN=MIN(XDMIN,XD(IN))
         XDMAX=MAX(XDMAX,XD(IN))
         YDMIN=MIN(YDMIN,YD(IN))
         YDMAX=MAX(YDMAX,YD(IN))
      ENDDO
C
      IF(NGXORG.EQ.0) XDMIN=0.D0
C
      XDMID=0.5D0*(XDMIN+XDMAX)
      XDMIN=XDMID+(XDMIN-XDMID)*0.999999D0
      XDMAX=XDMID+(XDMAX-XDMID)*0.999999D0
      YDMID=0.5D0*(YDMIN+YDMAX)
      YDMIN=YDMID+(YDMIN-YDMID)*0.999999D0
      YDMAX=YDMID+(YDMAX-YDMID)*0.999999D0
C
      DO IN=1,NNOD
         GXD(IN)=GCLIP(XD(IN))
         GYD(IN)=GCLIP(YD(IN))
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE RANGE OF WINDOW ******
C
      SUBROUTINE WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
C
      INCLUDE 'wfcomm.inc'
C
      DXLEN= XDMAX-XDMIN
      DYLEN= YDMAX-YDMIN
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
C     ***** SET PAGE-RECT AND VAR-RECT *****
C
      SUBROUTINE WFGVEW(NW,NWMAX,GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &                           GXMIN,GXMAX,GYMIN,GYMAX)
C
      INCLUDE 'wfcomm.inc'
C
      CALL WFGWIN(1,1,PXMIN,PXMAX,PYMIN,PYMAX)
      GPXMIN=GCLIP(PXMIN)
      GPXMAX=GCLIP(PXMAX)
      GPYMIN=GCLIP(PYMIN)+1.0
      GPYMAX=GCLIP(PYMAX)+1.0
C
      DXLEN= XDMAX-XDMIN
      DYLEN= YDMAX-YDMIN
      DRATIO=DYLEN/DXLEN
      PRATIO=FRATIO*(PYMAX-PYMIN)/(PXMAX-PXMIN)
      IF(DRATIO.GT.PRATIO) THEN
         GYMIN=GCLIP(YDMIN)
         GYMAX=GCLIP(YDMAX)
         XMID=0.5D0*(XDMIN+XDMAX)
         XLEN=DYLEN/PRATIO
         GXMIN=GCLIP(XMID-0.5D0*XLEN)
         GXMAX=GCLIP(XMID+0.5D0*XLEN)
      ELSE
         GXMIN=GCLIP(XDMIN)
         GXMAX=GCLIP(XDMAX)
         YMID=0.5D0*(YDMIN+YDMAX)
         YLEN=DXLEN*PRATIO
         GYMIN=GCLIP(YMID-0.5D0*YLEN)
         GYMAX=GCLIP(YMID+0.5D0*YLEN)
      ENDIF
      RETURN
      END
C
C     ****** PAINT 2D STRUCTURE OF FIELD AND POWER ******
C
      SUBROUTINE WFGPFX(NW,NWMAX,GV,KWD,IDP)
C
      INCLUDE 'wfcomm.inc'
C
      PARAMETER(NGXM=101,NGYM=101)
      PARAMETER(NSTEPM=101)
      PARAMETER(NRGBA=5,NRGBB=7)
      DIMENSION GV(NNODM),GZ(NGXM,NGYM)
      DIMENSION GAX(NGXM),GAY(NGYM)
      DIMENSION GZL(NSTEPM),GRGBL(3,0:NSTEPM)
      DIMENSION GRGBA(3,NRGBA),GLA(NRGBA)
      DIMENSION GRGBB(3,NRGBB),GLB(NRGBB)
      DIMENSION KAF(2,NELMM),KA(4,NGXM,NGYM)
      REAL*8 A(3),B(3),C(3)
      CHARACTER KWD*(NCHM)
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
      IF(NGRAPH.EQ.1) THEN
         GPXMIN=GCLIP(PXMIN)+0.3
         GPXMAX=GCLIP(PXMAX)-0.3
         GPYMIN=GCLIP(PYMIN)+1.5
         GPYMAX=GCLIP(PYMAX)
C
         DXLEN= XDMAX-XDMIN
         DYLEN= YDMAX-YDMIN
         DRATIO=DYLEN/DXLEN
         PRATIO=FRATIO*(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
         IF(DRATIO.GT.PRATIO) THEN
            GYMIN=GCLIP(YDMIN)
            GYMAX=GCLIP(YDMAX)
            XMID=0.5D0*(XDMIN+XDMAX)
            XLEN=DYLEN/PRATIO
            GXMIN=GCLIP(XMID-0.5D0*XLEN)
            GXMAX=GCLIP(XMID+0.5D0*XLEN)
         ELSE
            GXMIN=GCLIP(XDMIN)
            GXMAX=GCLIP(XDMAX)
            YMID=0.5D0*(YDMIN+YDMAX)
            YLEN=DXLEN*PRATIO
            GYMIN=GCLIP(YMID-0.5D0*YLEN)
            GYMAX=GCLIP(YMID+0.5D0*YLEN)
         ENDIF
      ELSEIF(NGRAPH.EQ.3) THEN
         DXLEN= DBLE(GXN2-GXN1)
         DYLEN= DBLE(GYN2-GYN1)
         DRATIO=DYLEN/DXLEN
         PRATIO=FRATIO*(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
         IF(DRATIO.GT.PRATIO) THEN
            GPYMIN=GCLIP(PYMIN)+1.5
            GPYMAX=GCLIP(PYMAX)
            PXMID=0.5D0*(PXMIN+PXMAX)
            PYLEN=PYMAX-PYMIN-1.5D0
            PXLEN=PYLEN/DRATIO
            GPXMIN=GCLIP(PXMID-0.5D0*PXLEN)
            GPXMAX=GCLIP(PXMID+0.5D0*PXLEN)
         ELSE
            GPXMIN=GCLIP(PXMIN)+0.3
            GPXMAX=GCLIP(PXMAX)-0.3
            PYMID=0.5D0*(PYMIN+PYMAX+1.5D0)
            PXLEN=PXMAX-PXMIN-0.6D0
            PYLEN=PXLEN*DRATIO
            GPYMIN=GCLIP(PYMID-0.5D0*PYLEN)
            GPYMAX=GCLIP(PYMID+0.5D0*PYLEN)
         ENDIF
      ELSE
         DXLEN= XDMAX-XDMIN
         DYLEN= YDMAX-YDMIN
         DRATIO=DYLEN/DXLEN
         PRATIO=FRATIO*(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
         GXMIN=GCLIP(XDMIN)
         GXMAX=GCLIP(XDMAX)
         GYMIN=GCLIP(YDMIN)
         GYMAX=GCLIP(YDMAX)
         IF(DRATIO.GT.PRATIO) THEN
            GPYMIN=GCLIP(PYMIN)+1.5
            GPYMAX=GCLIP(PYMAX)
            PXMID=0.5D0*(PXMIN+PXMAX)
            PYLEN=PYMAX-PYMIN-1.5D0
            PXLEN=PYLEN/DRATIO*FRATIO
            GPXMIN=GCLIP(PXMID-0.5D0*PXLEN)
            GPXMAX=GCLIP(PXMID+0.5D0*PXLEN)
         ELSE
            GPXMIN=GCLIP(PXMIN)+0.3
            GPXMAX=GCLIP(PXMAX)-0.3
            PYMID=0.5D0*(PYMIN+PYMAX+1.5D0)
            PXLEN=PXMAX-PXMIN-0.6D0
            PYLEN=PXLEN*DRATIO/FRATIO
            GPYMIN=GCLIP(PYMID-0.5D0*PYLEN)
            GPYMAX=GCLIP(PYMID+0.5D0*PYLEN)
         ENDIF
         IF(NGRAPH.EQ.2) GPXMAX=GPXMAX-0.5
      ENDIF
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
      IF(NGRAPH.NE.3) THEN
         CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &               GXMIN,GXMAX,GYMIN,GYMAX)
C
         CALL SETLIN(0,0,7)
         CALL WFGBDY
C
         CALL SETLIN(0,0,4)
         CALL WFGANT
C
         CALL OFFVEW
      ENDIF
C
      IF(NGRAPH.EQ.1) THEN
         GZMIN= 1.E30
         GZMAX=-1.E30
         IF(IDP.EQ.0) THEN
            DO IN=1,NNOD
               GZMIN=MIN(GZMIN,GV(IN))
               GZMAX=MAX(GZMAX,GV(IN))
            ENDDO
         ELSE
            DO IN=1,NNOD
               IF(KNODP(IN).LE.9) THEN
                  GZMIN=MIN(GZMIN,GV(IN))
                  GZMAX=MAX(GZMAX,GV(IN))
               ENDIF
            ENDDO
         ENDIF
         CALL GQSCAL(GZMIN,GZMAX,GQZMIN,GQZMAX,GZSCAL)
         GZDEL=0.5*GZSCAL
         IF(KWD(1:3).EQ.'EZA') GZDEL=0.05*GZSCAL
         ISTEP=INT((GQZMAX-GQZMIN)/GZDEL)
         IF(IDP.EQ.2) THEN
            GZDEL=1.0
            ISTEP=10
         ENDIF
C
         IF(GZMIN*GZMAX.GE.0.0) THEN
            IF(GZMAX.GT.0.D0) THEN
               GZORG=GQZMIN-1.E-5*GZDEL
            ELSE
               GZORG=GQZMAX+1.E-5*GZDEL
            ENDIF
         ELSE
            IF(IDP.EQ.2) THEN
               GZORG=GZDEL
            ELSE
               GZORG=0.5*GZDEL
            ENDIF
         ENDIF
      ENDIF
C
      IF(NGRAPH.NE.1) THEN
         NGXMAX=41
         NGYMAX=101
         DX=(XDMAX-XDMIN)/(NGXMAX-1)
         DY=(YDMAX-YDMIN)/(NGYMAX-1)
         DO NGY=1,NGYMAX
            GAY(NGY)=GCLIP(YDMIN+DY*(NGY-1))
         ENDDO
         DO NGX=1,NGXMAX
            GAX(NGX)=GCLIP(XDMIN+DX*(NGX-1))
         ENDDO
C
         DO NGY=1,NGYMAX
            Y=YDMIN+DY*(NGY-1)
            DO NGX=1,NGXMAX
               X=XDMIN+DX*(NGX-1)
               CALL WFFEP(X,Y,IE)
               IF(IE.LE.0) THEN
                  GZ(NGX,NGY)=0.0
               ELSE
                  CALL WFABC(IE,A,B,C,S)
                  GVL=0.0
                  DO N=1,3
                     IN=IELM(N,IE)
                     WEIGHT=A(N)+B(N)*X+C(N)*Y
                     GVL=GVL+SNGL(WEIGHT)*GV(IN)
                  ENDDO
                  GZ(NGX,NGY)=GVL
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C
      IF(NGRAPH.EQ.4) THEN
         CALL GMNMX2(GZ,NGXM,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
         GMAX=MAX(ABS(GZMIN),ABS(GZMAX))
         IF(GMAX.LE.0.D0) GOTO 1000
         GZDEL=GZMAX-GZMIN
         IF(ABS(GZDEL)/GMAX.LE.1.D-6) THEN
            GZDEL=0.1*GMAX
            GZORG=GZMIN
            ISTEP=0
         ELSE
            CALL GQSCAL(GZMIN,GZMAX,GQZMIN,GQZMAX,GZSCAL)
            GZDEL=0.5*GZSCAL
            IF(KWD(1:3).EQ.'EZA') GZDEL=0.05*GZSCAL
            ISTEP=INT((GQZMAX-GQZMIN)/GZDEL)
            IF(GZMIN*GZMAX.GE.0.0) THEN
               IF(GZMAX.GT.0.D0) THEN
                  GZORG=GQZMIN
               ELSE
                  GZORG=GQZMAX
               ENDIF
            ELSE
               GZORG=0.5*GZDEL
            ENDIF
         ENDIF
      ENDIF
C
      IF(NGRAPH.EQ.2) THEN
         ISTEP=50
      ENDIF
C
      IF(MODIFY.EQ.1) THEN
         CALL CHMODE
         WRITE(6,'(A,1P2E12.4)') '#MOD# GZMIN,GZMAX=',GZMIN,GZMAX
         WRITE(6,'(A,1P2E12.4,I5)') '#MOD# GZORG,GZDEL,ISTEP=',
     &                              GZORG,GZDEL,ISTEP
         WRITE(6,'(A)')          '#MOD# GZORG,GZDEL,ISTEP ?'
         READ(5,*,ERR=9,END=9) GZORG,GZDEL,ISTEP
    9    CONTINUE
         CALL GRMODE
      ENDIF
C
      IF(NGRAPH.EQ.1) THEN
         CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &               GXMIN,GXMAX,GYMIN,GYMAX)
C
         CALL SETLIN(0,0,7)
         CALL GFRAME
         CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.1,9)
         CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.1,9)
C
         IF(MODIFY.EQ.1) THEN
            GVL=GZORG
            DO I=1,ISTEP
               IF(GVL.GT.0.D0) THEN
                  CALL SETLIN(0,0,6)
                  IF(IDP.EQ.0) THEN
                     CALL CONTPF(GV,GXD,GYD,NNOD,IELM,NELM,
     &                           GZL,GZDEL,1,0,KAF)
                  ELSE
                     CALL CONTPF(GV,GXD,GYD,NNOD,IELMP,NELMP,
     &                           GZL,GZDEL,1,0,KAF)
                  ENDIF
               ELSEIF(GVL.LT.0.D0) THEN
                  CALL SETLIN(0,0,5)
                  IF(IDP.EQ.0) THEN
                     CALL CONTPF(GV,GXD,GYD,NNOD,IELM,NELM,
     &                           GVL,-GZDEL,1,2,KAF)
                  ELSE
                     CALL CONTPF(GV,GXD,GYD,NNOD,IELMP,NELMP,
     &                           GVL,-GZDEL,1,2,KAF)
                  ENDIF
               ELSE
                  CALL SETLIN(0,0,7)
                  IF(IDP.EQ.0) THEN
                     CALL CONTPF(GV,GXD,GYD,NNOD,IELM,NELM,
     &                           GVL,GZDEL,1,4,KAF)
                  ELSE
                     CALL CONTPF(GV,GXD,GYD,NNOD,IELMP,NELMP,
     &                           GVL,GZDEL,1,4,KAF)
                  ENDIF
               ENDIF
               GVL=GVL+GZDEL
            ENDDO
         ELSE
            IF(GZMIN*GZMAX.GE.0.0) THEN
               IF(GZMAX.GT.0.D0) THEN
                  CALL SETLIN(0,0,6)
                  IF(IDP.EQ.0) THEN
                     CALL CONTPF(GV,GXD,GYD,NNOD,IELM,NELM,
     &                           GZORG,GZDEL,ISTEP,0,KAF)
                  ELSEIF(IDP.LE.2) THEN
                     CALL CONTPF(GV,GXD,GYD,NNOD,IELMP,NELMP,
     &                           GZORG,GZDEL,ISTEP,0,KAF)
                  ENDIF
               ELSE
                  CALL SETLIN(0,0,5)
                  IF(IDP.EQ.0) THEN
                     CALL CONTPF(GV,GXD,GYD,NNOD,IELM,NELM,
     &                           GZORG,-GZDEL,ISTEP,2,KAF)
                  ELSEIF(IDP.LE.2) THEN
                     CALL CONTPF(GV,GXD,GYD,NNOD,IELMP,NELMP,
     &                           GZORG,-GZDEL,ISTEP,2,KAF)
                  ENDIF
               ENDIF
            ELSE
               CALL SETLIN(0,0,6)
               IF(IDP.EQ.0) THEN
                  CALL CONTPF(GV,GXD,GYD,NNOD,IELM,NELM,
     &                        GZORG,GZDEL,ISTEP,0,KAF)
               ELSEIF(IDP.LE.2) THEN
                  CALL CONTPF(GV,GXD,GYD,NNOD,IELMP,NELMP,
     &                        GZORG,GZDEL,ISTEP,0,KAF)
               ENDIF
               IF(IDP.EQ.2) THEN
                  CALL SETLIN(0,0,7)
                  CALL CONTPF(GV,GXD,GYD,NNOD,IELM,NELM,
     &                        0.0,GZDEL,1,4,KAF)
               ENDIF
               IF(IDP.EQ.3) THEN
                  CALL SETLIN(0,0,7)
                  CALL CONTPF(GV,GXD,GYD,NNOD,IELM,NELM,
     &                        0.0,GZDEL,1,0,KAF)
               ENDIF
               CALL SETLIN(0,0,5)
               IF(IDP.EQ.0) THEN
                  CALL CONTPF(GV,GXD,GYD,NNOD,IELM,NELM,
     &                        -GZORG,-GZDEL,ISTEP,2,KAF)
               ELSEIF(IDP.LE.2) THEN
                  CALL CONTPF(GV,GXD,GYD,NNOD,IELMP,NELMP,
     &                        -GZORG,-GZDEL,ISTEP,2,KAF)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
C
      IF(NGRAPH.EQ.2) THEN
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
     &               GXMIN,GXMAX,GYMIN,GYMAX)
         CALL CONTF2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,GZL,GRGBL,ISTEP,0)
         CALL SETLIN(0,0,7)
         CALL GFRAME
         CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.1,9)
         CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.1,9)
C
         CALL RGBBAR(GPXMAX+0.2,GPXMAX+0.5,GPYMIN,GPYMAX,
     &               GRGBL,ISTEP+1,1)
C
         CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &               GXMIN,GXMAX,GYMIN,GYMAX)
         CALL SETLIN(0,0,7)
         CALL WFGBDY
         CALL OFFVEW
      ENDIF
C
      IF(NGRAPH.EQ.3) THEN
         CALL GMNMX2(GZ,NGXM,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
         GMAX=MAX(ABS(GZMIN),ABS(GZMAX))
         IF(GMAX.LE.0.D0) GOTO 1000
         CALL GQSCAL(GZMIN,GZMAX,GQZMIN,GQZMAX,GZSCAL)
C
         CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &               GXN1,GXN2,GYN1,GYN2)
C
         CALL SETLIN(0,0,7)
         CALL PERSE1(GZ,NGXM,NGXMAX,NGYMAX,GQZMIN,GQZMAX,IXY,IDN,
     &               GXN,GYN,GZN,GA,GB,GC,GD,GE)
      ENDIF
C
      IF(NGRAPH.EQ.4) THEN
         CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &               GXMIN,GXMAX,GYMIN,GYMAX)
         CALL GFRAME
         CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.1,9)
         CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.1,9)
C
         IF(MODIFY.EQ.1) THEN
            GVL=GZORG
            DO I=1,ISTEP
               IF(GVL.GT.0.D0) THEN
                  CALL SETLIN(0,0,6)
                  CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,
     &                        GVL,GZDEL,1,0,0,KA)
               ELSEIF(GVL.LT.0.D0) THEN
                  CALL SETLIN(0,0,5)
                  CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,
     &                        GVL,GZDEL,1,0,2,KA)
               ELSE
                  CALL SETLIN(0,0,7)
                  CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,
     &                        GVL,GZDEL,1,0,4,KA)
               ENDIF
               GVL=GVL+GZDEL
            ENDDO
         ELSE
            IF(GZMIN*GZMAX.GE.0.0) THEN
               IF(GZMAX.GT.0.D0) THEN
                  CALL SETLIN(0,0,6)
                  CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,
     &                        GZORG,GZDEL,ISTEP,0,0,KA)
               ELSE
                  CALL SETLIN(0,0,5)
                  CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,
     &                        GZORG,-GZDEL,ISTEP,0,2,KA)
               ENDIF
            ELSE
               CALL SETLIN(0,0,6)
               CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,
     &                     GZORG,GZDEL,ISTEP,0,0,KA)
               CALL SETLIN(0,0,5)
               CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,
     &                     -GZORG,-GZDEL,ISTEP,0,2,KA)
            ENDIF
         ENDIF
      ENDIF
C
 1000 CALL SETLIN(0,0,7)
      CALL SETCHS(0.25,0.0)
      GXPOS=0.5*(GPXMIN+GPXMAX)-6.5*0.25
      GYPOS=GPYMIN-0.35
      CALL MOVE(GXPOS,GYPOS)
      IF(KWD(1:1).EQ.'P') THEN
         CALL TEXT(KWD(1:3),3)
         CALL TEXT(' ',1)
         CALL TEXT(KWD(4:NCHM),NCHM-3)
      ELSE
         CALL TEXT(KWD(1:2),2)
         CALL TEXT(' ',1)
         CALL TEXT(KWD(3:NCHM),NCHM-2)
      ENDIF
      GYPOS=GYPOS-0.35
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT('MAX=',4)
         CALL NUMBR(GZMAX,'(1PE9.2)',9)
      GYPOS=GYPOS-0.35
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT('MIN=',4)
         CALL NUMBR(GZMIN,'(1PE9.2)',9)
      GYPOS=GYPOS-0.35
      CALL MOVE(GXPOS,GYPOS)
         CALL TEXT('STP=',4)
         CALL NUMBR(GZDEL,'(1PE9.2)',9)
      RETURN
      END
C
C     ****** DRAW 1D PROFILE ******
C
      SUBROUTINE WFGPPR(NW,NWMAX,GZ,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION GZ(NNODM)
      DIMENSION GX(NRM),GY(NRM,3)
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
      GPXMIN=GCLIP(PXMIN)+2.2
      GPXMAX=GCLIP(PXMAX)-0.3
      GPYMIN=GCLIP(PYMIN)+1.5
      GPYMAX=GCLIP(PYMAX)
C
      IF(KWD(1:1).EQ.'E'.OR.
     &   KWD(1:1).EQ.'B'.OR.
     &   KWD(1:1).EQ.'A') THEN
         IGMAX=3
      ELSE
         IGMAX=1
      ENDIF
      DO NR=1,NRMAX
         GX(NR)  =GZ(        NR)
         DO IG=1,IGMAX
            GY(NR,IG)=GZ(IG*NRMAX+NR)
         ENDDO
      ENDDO
      CALL GMNMX1(GX,1,NRMAX,1,GXMIN1,GXMAX1)
      CALL GMNMX2(GY,NRM,1,NRMAX,1,1,IGMAX,1,GYMIN1,GYMAX1)
      CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP)
      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
C
      IF(MODIFY.EQ.1) THEN
         CALL CHMODE
         WRITE(6,'(A,1P3E12.4)') 
     &        '#MOD# GXMIN,GXMAX,GXSTEP=',GXMIN,GXMAX,GXSTEP
         READ(5,*,ERR=8,END=8) GXMIN,GXMAX,GXSTEP
    8    CONTINUE
         WRITE(6,'(A,1P3E12.4)') 
     &        '#MOD# GYMIN,GYMAX,GYSTEP=',GYMIN,GYMAX,GYSTEP
         READ(5,*,ERR=9,END=9) GYMIN,GYMAX,GYSTEP
    9    CONTINUE
         CALL GRMODE
      ENDIF
C
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
      CALL GVALUE(GXORG,4*GXSTEP,0.0,0.0,NGVLEN(4*GXSTEP))
      CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGVLEN(2*GYSTEP))
C
      DO IG=1,IGMAX
         CALL SETLIN(0,0,8-IG)
         CALL GPLOTP(GX,GY(1,IG),1,NRMAX,1,0,0,IPAT(IG))
      ENDDO
C
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.25,0.0)
      GXPOS=0.5*(GPXMIN+GPXMAX)-7.5*0.25
      GYPOS=GPYMIN-1.0
      CALL MOVE(GXPOS,GYPOS)
      IF(KWD(1:1).EQ.'P') THEN
         CALL TEXT(KWD(1:3),3)
         IF(KWD(4:4).EQ.'X') THEN
            CALL TEXT('(X): Y=',7)
         ELSE
            CALL TEXT('(Y): X=',7)
         ENDIF
         CALL TEXT(KWD(5:NCHM),NCHM-4)
      ELSE
         CALL TEXT(KWD(1:2),2)
         IF(KWD(3:3).EQ.'X') THEN
            CALL TEXT('(X): Y=',7)
         ELSE
            CALL TEXT('(Y): X=',7)
         ENDIF
         CALL TEXT(KWD(4:NCHM),NCHM-3)
      ENDIF
      RETURN
      END
C
C     ****** DRAW 1D PROFILE ******
C
      SUBROUTINE WFGPPRX(NW,NWMAX,GZ1,GZ2,GZ3,GZ4,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION GZ1(NNODM),GZ2(NNODM),GZ3(NNODM),GZ4(NNODM)
      DIMENSION GX(NRM),GY(NRM,4)
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
      GPXMIN=GCLIP(PXMIN)+2.2
      GPXMAX=GCLIP(PXMAX)-0.3
      GPYMIN=GCLIP(PYMIN)+1.5
      GPYMAX=GCLIP(PYMAX)
C
      IF(KWD(1:1).EQ.'E'.OR.
     &   KWD(1:1).EQ.'B'.OR.
     &   KWD(1:1).EQ.'A') THEN
         IGMAX=3
      ELSE
         IGMAX=4
      ENDIF
      DO NR=1,NRMAX
         GX(NR)  =GZ1(      NR)
         GY(NR,1)=GZ1(NR+NRMAX)
         GY(NR,2)=GZ2(NR+NRMAX)
         GY(NR,3)=GZ3(NR+NRMAX)
         GY(NR,4)=GZ4(NR+NRMAX)
      ENDDO
      CALL GMNMX1(GX,1,NRMAX,1,GXMIN1,GXMAX1)
      CALL GMNMX2(GY,NRM,1,NRMAX,1,1,IGMAX,1,GYMIN1,GYMAX1)
      CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP)
      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
C
      IF(MODIFY.EQ.1) THEN
         CALL CHMODE
         WRITE(6,'(A,1P3E12.4)') 
     &        '#MOD# GXMIN,GXMAX,GXSTEP=',GXMIN,GXMAX,GXSTEP
         READ(5,*,ERR=8,END=8) GXMIN,GXMAX,GXSTEP
    8    CONTINUE
         WRITE(6,'(A,1P3E12.4)') 
     &        '#MOD# GYMIN,GYMAX,GYSTEP=',GYMIN,GYMAX,GYSTEP
         READ(5,*,ERR=9,END=9) GYMIN,GYMAX,GYSTEP
    9    CONTINUE
         CALL GRMODE
      ENDIF
C
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
      CALL GVALUE(GXORG,4*GXSTEP,0.0,0.0,NGVLEN(4*GXSTEP))
      CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGVLEN(2*GYSTEP))
C
      DO IG=1,IGMAX
         CALL SETLIN(0,0,8-IG)
         CALL GPLOTP(GX,GY(1,IG),1,NRMAX,1,0,0,IPAT(IG))
      ENDDO
C
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.25,0.0)
      GXPOS=0.5*(GPXMIN+GPXMAX)-7.5*0.25
      GYPOS=GPYMIN-1.0
      CALL MOVE(GXPOS,GYPOS)
      IF(KWD(1:1).EQ.'P') THEN
         CALL TEXT(KWD(1:3),3)
         IF(KWD(4:4).EQ.'X') THEN
            CALL TEXT('(X): Y=',7)
         ELSE
            CALL TEXT('(Y): X=',7)
         ENDIF
         CALL TEXT(KWD(5:NCHM),NCHM-4)
      ELSE
         CALL TEXT(KWD(1:2),2)
         IF(KWD(3:3).EQ.'X') THEN
            CALL TEXT('(X): Y=',7)
         ELSE
            CALL TEXT('(Y): X=',7)
         ENDIF
         CALL TEXT(KWD(4:NCHM),NCHM-3)
      ENDIF
      RETURN
      END
C
C     ****** DRAW 1D PROFILE ******
C
      SUBROUTINE WFGPFR(NW,NWMAX,GZ,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION GZ(NNODM)
      DIMENSION GX(NRM),GY(NRM,3)
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
      GPXMIN=GCLIP(PXMIN)+2.2
      GPXMAX=GCLIP(PXMAX)-0.3
      GPYMIN=GCLIP(PYMIN)+1.5
      GPYMAX=GCLIP(PYMAX)
C
      IGMAX=2
      DO NR=1,NRMAX
         GX(NR)  =GZ(        NR)
         DO IG=1,IGMAX
            GY(NR,IG)=GZ(IG*NRMAX+NR)
         ENDDO
      ENDDO
      CALL GMNMX1(GX,1,NRMAX,1,GXMIN1,GXMAX1)
      CALL GMNMX2(GY,NRM,1,NRMAX,1,1,IGMAX,1,GYMIN1,GYMAX1)
      GMAX=MAX(ABS(GYMIN1),ABS(GYMAX1))
      CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP)
      CALL GQSCAL(-GMAX,GMAX,GYMIN,GYMAX,GYSTEP)
C
      IF(MODIFY.EQ.1) THEN
         CALL CHMODE
         WRITE(6,'(A,1P3E12.4)') 
     &        '#MOD# GXMIN,GXMAX,GXSTEP=',GXMIN,GXMAX,GXSTEP
         READ(5,*,ERR=8,END=8) GXMIN,GXMAX,GXSTEP
    8    CONTINUE
         WRITE(6,'(A,1P3E12.4)') 
     &        '#MOD# GYMIN,GYMAX,GYSTEP=',GYMIN,GYMAX,GYSTEP
         READ(5,*,ERR=9,END=9) GYMIN,GYMAX,GYSTEP
    9    CONTINUE
         CALL GRMODE
      ENDIF
C
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
      CALL GVALUE(GXORG,4*GXSTEP,0.0,0.0,NGVLEN(4*GXSTEP))
      CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGVLEN(2*GYSTEP))
C
      DO IG=1,IGMAX
         CALL SETLIN(0,0,8-IG)
         CALL GPLOTP(GX,GY(1,IG),1,NRMAX,1,0,0,IPAT(IG))
      ENDDO
C
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.25,0.0)
      GXPOS=0.5*(GPXMIN+GPXMAX)-7.5*0.25
      GYPOS=GPYMIN-1.0
      CALL MOVE(GXPOS,GYPOS)
      IF(KWD(1:1).EQ.'P') THEN
         CALL TEXT(KWD(1:3),3)
         IF(KWD(4:4).EQ.'X') THEN
            CALL TEXT('(X): Y=',7)
         ELSE
            CALL TEXT('(Y): X=',7)
         ENDIF
         CALL TEXT(KWD(5:NCHM),NCHM-4)
      ELSE
         CALL TEXT(KWD(1:2),2)
         IF(KWD(3:3).EQ.'X') THEN
            CALL TEXT('(X): Y=',7)
         ELSE
            CALL TEXT('(Y): X=',7)
         ENDIF
         CALL TEXT(KWD(4:NCHM),NCHM-3)
      ENDIF
      RETURN
      END
C
C     ****** DRAW 1D PROFILE (AMP AND PHASE) ******
C
      SUBROUTINE WFGPAR(NW,NWMAX,GZ,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION GZ(NNODM)
      DIMENSION GX(NRM),GY(NRM,2)
      DIMENSION IPAT(2)
      CHARACTER KWD*(NCHM)
      DATA IPAT/0,2/
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
      GPXMIN=GCLIP(PXMIN)+2.2
      GPXMAX=GCLIP(PXMAX)-0.3
      GPYMIN=GCLIP(PYMIN)+1.5
      GPYMAX=GCLIP(PYMAX)
C
      DO NR=1,NRMAX
         GX(NR)  =GZ(        NR)
         GY(NR,1)=GZ(3*NRMAX+NR)
         GY(NR,2)=ATAN2(GZ(2*NRMAX+NR),GZ(NRMAX+NR))
      ENDDO
C
      CALL GMNMX1(GX,     1,NRMAX,1,GXMIN1,GXMAX1)
      CALL GMNMX1(GY(1,1),1,NRMAX,1,GYMIN1,GYMAX1)
      GYMIN1=0.0
      CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP)
      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
      GYMIN=0.0
C
      IF(MODIFY.EQ.1) THEN
         CALL CHMODE
         WRITE(6,'(A,1P3E12.4)') 
     &        '#MOD# GXMIN,GXMAX,GXSTEP=',GXMIN,GXMAX,GXSTEP
         READ(5,*,ERR=8,END=8) GXMIN,GXMAX,GXSTEP
    8    CONTINUE
         WRITE(6,'(A,1P3E12.4)') 
     &        '#MOD# GYMIN,GYMAX,GYSTEP=',GYMIN,GYMAX,GYSTEP
         READ(5,*,ERR=9,END=9) GYMIN,GYMAX,GYSTEP
    9    CONTINUE
         CALL GRMODE
      ENDIF
C
      IF(GXMIN*GXMAX.LT.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GXMIN
      ENDIF
      GYORG=0.0
C
      CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
C
      CALL SETLIN(0,0,4)
      CALL SETCHS(0.25,0.0)
      CALL GFRAME
      CALL GSCALE(GXORG,GXSTEP,0.0,0.0,0.2,9)
      CALL GSCALE(0.0,0.0,GYORG,GYSTEP,0.2,1)
      CALL GVALUE(GXORG,4*GXSTEP,0.0,0.0,NGVLEN(4*GXSTEP))
      CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGVLEN(2*GYSTEP))
C
      CALL SETLIN(0,0,6)
      CALL GPLOTP(GX,GY(1,1),1,NRMAX,1,0,0,IPAT(1))
C
      GPI=SNGL(PI)
      CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,-GPI,GPI)
      CALL SETLIN(0,0,4)
      CALL GSCALE(0.0,0.0,0.0,GPI/4,0.2,5)
      CALL GVALUE(0.0,0.0,0.0,GPI/2,2+400)
C
      CALL SETLIN(0,0,5)
      IMARK=2
      CALL GPLOTP(GX,GY(1,2),1,NRMAX,1,IMARK,1,IPAT(2))
C
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.25,0.0)
      GXPOS=0.5*(GPXMIN+GPXMAX)-7.5*0.25
      GYPOS=GPYMIN-1.0
      CALL MOVE(GXPOS,GYPOS)
      IF(KWD(1:1).EQ.'P') THEN
         CALL TEXT(KWD(1:3),3)
         IF(KWD(4:4).EQ.'X') THEN
            CALL TEXT('(X): Y=',7)
         ELSE
            CALL TEXT('(Y): X=',7)
         ENDIF
         CALL TEXT(KWD(5:NCHM),NCHM-4)
      ELSE
         CALL TEXT(KWD(1:2),2)
         IF(KWD(3:3).EQ.'X') THEN
            CALL TEXT('(X): Y=',7)
         ELSE
            CALL TEXT('(Y): X=',7)
         ENDIF
         CALL TEXT(KWD(4:NCHM),NCHM-3)
      ENDIF
      RETURN
      END
C
C     ****** DRAW 1D HISTORY ******
C
      SUBROUTINE WFGPFT(NW,NWMAX,GY,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION GX(NGTM),GY(NGTM,3)
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
C      PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-2.5D0)
C      IF(PRATIO.LT.1.D0) THEN
C         PYLEN=PYMAX-PYMIN-1.5D0
C         PXLEN=PYLEN+2.5D0
C         PXMID=0.5D0*(PXMIN+PXMAX)
C         PXMIN=PXMID-0.5D0*PXLEN
C         PXMAX=PXMID+0.5D0*PXLEN
C      ENDIF
C      
      GPXMIN=GCLIP(PXMIN)+2.2
      GPXMAX=GCLIP(PXMAX)-0.3
      GPYMIN=GCLIP(PYMIN)+1.5
      GPYMAX=GCLIP(PYMAX)
C
      IGMAX=3
      DO NGT=1,NGTMAX
         GX(NGT)  =GCLIP(TG(NGT))
      ENDDO
      CALL GMNMX1(GX,1,NGTMAX,1,GXMIN1,GXMAX1)
      CALL GMNMX2(GY,NGTM,1,NGTMAX,1,1,IGMAX,1,GYMIN1,GYMAX1)
      CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP)
      CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
C
      IF(MODIFY.EQ.1) THEN
         CALL CHMODE
         WRITE(6,'(A,1P3E12.4)') 
     &        '#MOD# GXMIN,GXMAX,GXSTEP=',GXMIN,GXMAX,GXSTEP
         READ(5,*,ERR=8,END=8) GXMIN,GXMAX,GXSTEP
    8    CONTINUE
         WRITE(6,'(A,1P3E12.4)') 
     &        '#MOD# GYMIN,GYMAX,GYSTEP=',GYMIN,GYMAX,GYSTEP
         READ(5,*,ERR=9,END=9) GYMIN,GYMAX,GYSTEP
    9    CONTINUE
         CALL GRMODE
      ENDIF
C
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
      CALL GVALUE(GXORG,4*GXSTEP,0.0,0.0,NGVLEN(4*GXSTEP))
      CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGVLEN(2*GYSTEP))
C
      DO IG=1,IGMAX
         CALL SETLIN(0,0,8-IG)
         CALL GPLOTP(GX,GY(1,IG),1,NGTMAX,1,0,0,IPAT(IG))
      ENDDO
C
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.25,0.0)
      GXPOS=0.5*(GPXMIN+GPXMAX)-7.5*0.25
      GYPOS=GPYMIN-1.0
      CALL MOVE(GXPOS,GYPOS)
      CALL TEXT(KWD(1:3),3)
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
C
      DO 10 NA=1,NAMAX
         REST(NA)=DBLE(CIMP(NA))
         REAT(NA)=DIMAG(CIMP(NA))
   10 CONTINUE
C
      GXMIN=0.
      GYMAX=18.0
      GRCHH=0.20
      CALL SETCHS(GRCHH,0.)
      GDX=15.*GRCHH
      GDY=1.5*GRCHH
      GX=GXMIN
      GY=GYMAX
C
      GY=GY-GDY
      CALL MOVE(GX,GY)
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
      GX=GX+GDX
      CALL MOVE(GX,GY)
      CALL TEXT('BB  =',5)
      CALL NUMBD(BB  ,'(F7.3)',7)
C
      GX=GX+GDX
      CALL MOVE(GX,GY)
      IF(MODELS.EQ.2) THEN
         CALL TEXT('RR  =',5)
         CALL NUMBD(RR   ,'(F7.3)',7)
      ELSEIF(MODELS.EQ.1) THEN
         CALL TEXT('RMIR=',5)
         CALL NUMBD(RMIR ,'(F7.3)',7)
      ELSEIF(MODELS.EQ.3) THEN
         CALL TEXT('RRC =',5)
         CALL NUMBD(RRC  ,'(F7.3)',7)
      ELSE
         CALL TEXT('MODELB= ',8)
         CALL NUMBI(MODELB,'(I4)',4)
      ENDIF
C
      GY=GY-GDY
      GX=GXMIN
      CALL MOVE(GX,GY)
      IF(NZMAX.EQ.1) THEN
         IF(MODELS.EQ.2.OR.MODELS.EQ.1) THEN
            CALL TEXT('NPHI=',5)
            CALL NUMBI(NPHI,'(I7)',7)
         ELSE
            CALL TEXT('RKZ =',5)
            CALL NUMBD(RKZ,'(F7.3)',7)
         ENDIF
      ELSE
         IF(MODELS.EQ.2.OR.MODELS.EQ.1) THEN
            CALL TEXT('NPHI SUM',8)
         ELSE
            CALL TEXT('RK  =',5)
            CALL NUMBD(RZ,'(F7.3)',7)
         ENDIF
      ENDIF
C
      GX=GX+GDX
      CALL MOVE(GX,GY)
      CALL TEXT('RA  =',5)
      CALL NUMBD(RA,'(F7.3)',7)
C
      GX=GX+GDX
      CALL MOVE(GX,GY)
      IF(MODELS.EQ.2) THEN
         CALL TEXT('QA  =',5)
         CALL NUMBD(QA  ,'(F7.3)',7)
      ELSEIF(MODELS.EQ.1) THEN
         CALL TEXT('ZBB =',5)
         CALL NUMBD(ZBB ,'(F7.3)',7)
      ELSEIF(MODELS.EQ.3) THEN
         CALL TEXT('H1  =',5)
         CALL NUMBD(H1  ,'(F7.3)',7)
      ENDIF
C
      GY=GY-GDY
      GX=GXMIN
      CALL MOVE(GX,GY)
      CALL TEXT(' ',1)
      CALL NUMBI(NNOD,'(I5)',5)
      CALL TEXT('/',1)
      CALL NUMBI(NELM,'(I5)',5)
C
      GX=GX+GDX
      CALL MOVE(GX,GY)
      CALL TEXT('T=',2)
      CALL NUMBD(T,'(1PE10.3)',10)
C
      GX=GX+GDX
      CALL MOVE(GX,GY)
      CALL TEXT('P=',2)
      CALL NUMBD(PWR,'(1PE10.3)',10)
C
      IF(NEVOL.EQ.0) THEN
C
      GX=GXMIN
      GY=GY-4.0*GDY
C
      CALL MOVE(GX,GY)
      CALL TEXT('IS',2)
      CALL TEXT('     PA    ',11)
      CALL TEXT('  PZ ',5)
      CALL TEXT('     PN    ',11)
      CALL TEXT('     PNS   ',11)
      CALL TEXT('  PZCL ',7)
C
      CALL TEXT('   PTPR  ',9)
      CALL TEXT('   PTPP  ',9)
      CALL TEXT('   PTS   ',9)
      CALL TEXT('    PABS   ',11)
C
      DO 20 IS=1,NSMAX
         GX=GXMIN
         GY=GY-GDY
         CALL MOVE(GX,GY)
         IST=IS
         CALL NUMBI(IST,'(I2)',2)
         CALL NUMBD(PA(IS),'(1PE11.3)',11)
         CALL NUMBD(PZ(IS),'(F5.0)',5)
         CALL NUMBD(PN(IS),'(1PE11.3)',11)
         CALL NUMBD(PNS(IS),'(1PE11.3)',11)
         CALL NUMBD(PZCL(IS),'(F7.3)',7)
C
         CALL NUMBD(PTPR(IS),'(F9.3)',9)
         CALL NUMBD(PTPP(IS),'(F9.3)',9)
         CALL NUMBD(PTS(IS),'(F9.3)',9)
         CALL NUMBD(PWRS(IS),'(1PE11.3)',11)
   20 CONTINUE
C
      ELSEIF(NEVOL.EQ.1.OR.NEVOL.EQ.2) THEN
C
      GX=GXMIN
      GY=GY-1.3*GDY
C
      CALL MOVE(GX,GY)
      CALL TEXT('MODELS,B,D,P,W,T=',17)
      CALL NUMBI(MODELS,'(I1)',1)
      CALL TEXT(',',1)
      CALL NUMBI(MODELB,'(I1)',1)
      CALL TEXT(',',1)
      CALL NUMBI(MODELD,'(I1)',1)
      CALL TEXT(',',1)
      CALL NUMBI(MODELP,'(I1)',1)
      CALL TEXT(',',1)
      CALL NUMBI(MODELW,'(I1)',1)
      CALL TEXT(',',1)
      CALL NUMBI(MODELT,'(I1)',1)
C
      GX=GX+3.0*GDX
      CALL MOVE(GX,GY)
      CALL TEXT('      PIN  =',12)
      CALL NUMBD(PIN,'(F7.1)',7)
C
      GX=GXMIN
      GY=GY-GDY
C
      CALL MOVE(GX,GY)
      CALL TEXT('PPN0 =',6)
      CALL NUMBD(PPN0,'(1PE10.3)',10)
      CALL TEXT(' PNE0 =',7)
      CALL NUMBD(PNE0,'(1PE10.3)',10)
      CALL TEXT(' PTE0 =',7)
      CALL NUMBD(PTE0,'(1PE10.3)',10)
      CALL TEXT(' PTI0 =',7)
      CALL NUMBD(PTI0,'(1PE10.3)',10)
C
      GX=GXMIN
      GY=GY-GDY
C
      CALL MOVE(GX,GY)
      CALL TEXT('PTN0 =',6)
      CALL NUMBD(PTN0,'(1PE10.3)',10)
      CALL TEXT(' PNES =',7)
      CALL NUMBD(PNES,'(1PE10.3)',10)
      CALL TEXT(' PTES =',7)
      CALL NUMBD(PTES,'(1PE10.3)',10)
      CALL TEXT(' PTIS =',7)
      CALL NUMBD(PTIS,'(1PE10.3)',10)
C
      GX=GXMIN
      GY=GY-GDY
C
      CALL MOVE(GX,GY)
      CALL TEXT('NEMAX=',6)
      CALL NUMBD(PNT(NGTMAX,2,1),'(1PE10.3)',10)
      CALL TEXT(' NIMAX=',7)
      CALL NUMBD(PNT(NGTMAX,2,2),'(1PE10.3)',10)
      CALL TEXT(' TEMAX=',7)
      CALL NUMBD(PTT(NGTMAX,2,1),'(1PE10.3)',10)
      CALL TEXT(' TIMAX=',7)
      CALL NUMBD(PTT(NGTMAX,2,2),'(1PE10.3)',10)
      CALL TEXT(' PHIMAX=',8)
      CALL NUMBD(PHIT(NGTMAX,2),'(1PE10.3)',10)
C
      GX=GXMIN
      GY=GY-GDY
C
      CALL MOVE(GX,GY)
      CALL TEXT('NEAVE=',6)
      CALL NUMBD(PNT(NGTMAX,1,1),'(1PE10.3)',10)
      CALL TEXT(' NIAVE=',7)
      CALL NUMBD(PNT(NGTMAX,1,2),'(1PE10.3)',10)
      CALL TEXT(' TEAVE=',7)
      CALL NUMBD(PTT(NGTMAX,1,1),'(1PE10.3)',10)
      CALL TEXT(' TIAVE=',7)
      CALL NUMBD(PTT(NGTMAX,1,2),'(1PE10.3)',10)
      CALL TEXT(' PHIAVE=',8)
      CALL NUMBD(PHIT(NGTMAX,1),'(1PE10.3)',10)
C
      GX=GXMIN+4.8*GDX
      GY=GYMAX-6.3*GDY
      CALL MOVE(GX,GY)
      CALL TEXT('IS  ',4)
      CALL TEXT('   PABS  ',9)
C
      DO 30 IS=1,NSMAX
         GX=GXMIN+4.8*GDX
         GY=GY-GDY
         CALL MOVE(GX,GY)
         IST=IS
         CALL NUMBI(IST,'(I2)',2)
         CALL NUMBD(PWRS(IS),'(1PE11.3)',11)
   30 CONTINUE
C
      ENDIF
C
      GX=GXMIN+45.*GRCHH
      GY=GYMAX-GDY
      CALL MOVE(GX,GY)
      CALL TEXT('IJ',2)
      CALL TEXT('   AJ  ',7)
      CALL TEXT('  PHASE ',8)
      CALL TEXT('     R     ',11)
      CALL TEXT('     X     ',11)
C
      DO 40 IJ=1,NAMAX
         GX=GXMIN+45.*GRCHH
         GY=GY-GDY
         CALL MOVE(GX,GY)
         IJT=IJ
         CALL NUMBI(IJT,'(I2)',2)
         CALL NUMBD(AJ(IJ),'(F7.1)',7)
         CALL NUMBD(APH(IJ),'(F7.1)',7)
         CALL NUMBD(REST(IJ),'(1PE11.3)',11)
         CALL NUMBD(REAT(IJ),'(1PE11.3)',11)
   40 CONTINUE
C
      RETURN
      END
C
C     ****** Draw Boundary Node Points ******
C
      SUBROUTINE WFGNOD(ID)
C
      INCLUDE 'wfcomm.inc'
C
      CALL SETLIN(0,0,2)
      IF(NBDYP.GT.0) THEN
         IN=IBDYP(NBDYP)
         GX=GCLIP(XD(IN))
         GY=GCLIP(YD(IN))
         CALL MOVE(GX,GY)
         DO IB=1,NBDYP
            IN=IBDYP(IB)
            GX=GCLIP(XD(IN))
            GY=GCLIP(YD(IN))
            CALL DRAW(GX,GY)
         ENDDO
         DO IB=1,NBDYP
            IN=IBDYP(IB)
            GX=GCLIP(XD(IN))
            GY=GCLIP(YD(IN))
            IF(ID.EQ.0) THEN
               KN=KNODW(IN)
            ELSE
               KN=KNODP(IN)
            ENDIF
            IF(KN.GE.1) THEN
               IF(KN.GE.1.AND.KN.LE.4) THEN
                  CALL SETLIN(0,0,7)
               ELSEIF(KN.GE.5.AND.KN.LE.8) THEN
                  CALL SETLIN(0,0,6)
               ELSEIF(KN.GE.9.AND.KN.LE.12) THEN
                  CALL SETLIN(0,0,5)
               ELSEIF(KN.GE.13.AND.KN.LE.16) THEN
                  CALL SETLIN(0,0,4)
               ENDIF
               CALL SETMKS(MOD(KN-1,4)+1,0.2)
               CALL MARK(GX,GY)
            ENDIF
         ENDDO
      ENDIF
C
      CALL SETLIN(0,0,1)
      IF(NBDYW.GT.0) THEN
         IN=IBDYW(NBDYW)
         GX=GCLIP(XD(IN))
         GY=GCLIP(YD(IN))
         CALL MOVE(GX,GY)
         DO IB=1,NBDYW
            IN=IBDYW(IB)
            GX=GCLIP(XD(IN))
            GY=GCLIP(YD(IN))
            CALL DRAW(GX,GY)
         ENDDO
         DO IB=1,NBDYW
            IN=IBDYW(IB)
            GX=GCLIP(XD(IN))
            GY=GCLIP(YD(IN))
            IF(ID.EQ.0) THEN
               KN=KNODW(IN)
            ELSE
               KN=KNODP(IN)
            ENDIF
            IF(KN.GE.1) THEN
               IF(KN.GE.1.AND.KN.LE.4) THEN
                  CALL SETLIN(0,0,7)
               ELSEIF(KN.GE.5.AND.KN.LE.8) THEN
                  CALL SETLIN(0,0,6)
               ELSEIF(KN.GE.9.AND.KN.LE.12) THEN
                  CALL SETLIN(0,0,5)
               ELSEIF(KN.GE.13.AND.KN.LE.16) THEN
                  CALL SETLIN(0,0,4)
               ENDIF
               CALL SETMKS(MOD(KN-1,4)+1,0.2)
               CALL MARK(GX,GY)
            ENDIF
         ENDDO
      ENDIF
C
      CALL INQVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
      CALL OFFVEW
C
      DO KN=1,13
         GX=0.2
         GY=16.0-KN*0.5
            IF(KN.GE.1) THEN
               IF(KN.GE.1.AND.KN.LE.4) THEN
                  CALL SETLIN(0,0,7)
               ELSEIF(KN.GE.5.AND.KN.LE.8) THEN
                  CALL SETLIN(0,0,6)
               ELSEIF(KN.GE.9.AND.KN.LE.12) THEN
                  CALL SETLIN(0,0,5)
               ELSEIF(KN.GE.13.AND.KN.LE.16) THEN
                  CALL SETLIN(0,0,4)
               ENDIF
               CALL SETMKS(MOD(KN-1,4)+1,0.2)
               CALL MARK(GX,GY)
            ENDIF
         CALL SETCHS(0.15,0.0)
         CALL MOVE(GX+0.5,GY)
         CALL NUMBI(KN,'(I2)',2)
      ENDDO
C
      CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
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
      IF(NBDYP.GT.0) THEN
         GX=GCLIP(XD(IBDYP(NBDYP)))
         GY=GCLIP(YD(IBDYP(NBDYP)))
         CALL MOVE(GX,GY)
         DO IB=1,NBDYP
            GX=GCLIP(XD(IBDYP(IB)))
            GY=GCLIP(YD(IBDYP(IB)))
            CALL DRAW(GX,GY)
         ENDDO
      ENDIF
C
      IF(NBDYW.GT.0) THEN
         GX=GCLIP(XD(IBDYW(NBDYW)))
         GY=GCLIP(YD(IBDYW(NBDYW)))
         CALL MOVE(GX,GY)
         DO IB=1,NBDYW
            GX=GCLIP(XD(IBDYW(IB)))
            GY=GCLIP(YD(IBDYW(IB)))
            CALL DRAW(GX,GY)
         ENDDO
      ENDIF
      RETURN
      END
C
C     ****** Draw Plasma Boundary ******
C
      SUBROUTINE WFGPLA
C
      INCLUDE 'wfcomm.inc'
C
      IF(MODELS.EQ.1) RETURN
C
      NPMAX=100
      DTH=2.D0*PI/NPMAX
      GX=GCLIP(RA)
      GY=0.0
      CALL MOVE(GX,GY)
      DO 10 I=1,NPMAX
         THETA=DTH*I
         GX=GCLIP(RA*COS(THETA))
         GY=GCLIP(RA*SIN(THETA)*RKAP)
         CALL DRAW(GX,GY)
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
         DO NANT=1,NAMAX
            GX=GCLIP(XJ0(1,NANT))
            GY=GCLIP(YJ0(1,NANT))
            IF(JNUM0(NANT).EQ.1) THEN
               CALL MOVE(GX-0.005,GY)
               CALL DRAW(GX+0.005,GY)
               CALL MOVE(GX,GY-0.005)
               CALL DRAW(GX,GY+0.005)
            ELSE
               CALL MOVE(GX,GY)
               DO I=2,JNUM0(NANT)
                  GX=GCLIP(XJ0(I,NANT))
                  GY=GCLIP(YJ0(I,NANT))
                  CALL DRAW(GX,GY)
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO NANT=1,NAMAX
            GX=GCLIP(XJ(1,NANT))
            GY=GCLIP(YJ(1,NANT))
            IF(JNUM(NANT).EQ.1) THEN
               CALL MOVE(GX-0.005,GY)
               CALL DRAW(GX+0.005,GY)
               CALL MOVE(GX,GY-0.005)
               CALL DRAW(GX,GY+0.005)
            ELSE
               CALL MOVE(GX,GY)
               DO I=2,JNUM(NANT)
                  GX=GCLIP(XJ(I,NANT))
                  GY=GCLIP(YJ(I,NANT))
                  CALL DRAW(GX,GY)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      RETURN
      END
C
C     ****** Draw Element Data ******
C
      SUBROUTINE WFGELM
C
      INCLUDE 'wfcomm.inc'
      DIMENSION GX(4),GY(4)
C
      CALL INQRGB(GRGB1,GRGB2,GRGB3)
      CALL SETRGB(0.0,0.5,0.0)
      DO IE=1,NELM
         IN1=ABS(IELM(1,IE))
         IN2=ABS(IELM(2,IE))
         IN3=ABS(IELM(3,IE))
         GX(1)=GCLIP(XD(IN1))
         GY(1)=GCLIP(YD(IN1))
         GX(2)=GCLIP(XD(IN2))
         GY(2)=GCLIP(YD(IN2))
         GX(3)=GCLIP(XD(IN3))
         GY(3)=GCLIP(YD(IN3))
         GX(4)=GX(1)
         GY(4)=GY(1)
         KE=KELM(IE)
         IF(KE.EQ.0) THEN
            CALL SETRGB(1.0,1.0,1.0)
         ELSEIF(KE.EQ.1) THEN
            CALL SETRGB(0.8,1.0,0.8)
            CALL POLY(GX,GY,4)
         ELSEIF(KE.EQ.2) THEN
            CALL SETRGB(1.0,0.8,0.8)
            CALL POLY(GX,GY,4)
         ELSEIF(KE.EQ.3) THEN
            CALL SETRGB(0.8,0.8,1.0)
            CALL POLY(GX,GY,4)
         ENDIF
         CALL SETRGB(0.0,1.0,0.0)
         CALL LINES(GX,GY,4)
      ENDDO
C
      IF(NDRAWD.GE.2) THEN
         CALL SETRGB(1.0,0.0,0.0)
         CALL SETCHR(0.2,0.15,0.2,0.,-30.)
         DO IE=1,NELM
            IN1=ABS(IELM(1,IE))
            IN2=ABS(IELM(2,IE))
            IN3=ABS(IELM(3,IE))
            XC=(XD(IN1)+XD(IN2)+XD(IN3))/3.D0
            YC=(YD(IN1)+YD(IN2)+YD(IN3))/3.D0
            GXC=GCLIP(XC)
            GYC=GCLIP(YC)
            IEL=IE
            CALL GNUMBI(GXC,GYC,IEL,2)
         ENDDO
      ENDIF
C
      IF(NDRAWD.GE.3) THEN
         CALL SETRGB(0.0,0.0,0.5)
         CALL SETCHS(0.2,0.)
         DO IN=1,NNOD
            GX1=GCLIP(XD(IN))
            GY1=GCLIP(YD(IN))
            INL=IN
            CALL GNUMBI(GX1,GY1,INL,0)
         ENDDO
      ENDIF
      CALL SETRGB(GRGB1,GRGB2,GRGB3)
C
      RETURN
      END
C
C     *********************
C
C     CLIPPING FOR GRAPHICS
C
C     *********************
C
      FUNCTION GCLIP(D)
      REAL*8 D
      IF(ABS(D).GT.1.D-15) THEN
         GCLIP=REAL(D)
      ELSE
         GCLIP=0.0
      ENDIF
      RETURN
      END
C
C     *****************************
C
C     OPTIMUM NUM LENGTH FOR GVALUE
C
C     *****************************
C
      FUNCTION NGVLEN(GSTEP)
C
      NGX = -INT(LOG10(DBLE(GSTEP*0.11)))
      IF(NGX.LT.-5)THEN
         NGX=-1
      ELSE
         IF(NGX.LT.0) NGX=0
         IF(NGX.GT.5) NGX=-1
      ENDIF
      NGVLEN=NGX
      RETURN
      END
