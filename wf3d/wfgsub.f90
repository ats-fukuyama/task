!     $Id$

!     ****** CALCULATE RANGE OF WINDOW ******

SUBROUTINE WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)

  use wfcomm
  implicit none
  integer :: NWW,NWMAX,NWYMAX,MIN,NW,NWX,NWY
  real(8) :: DXLEN,DYLEN,DRATIO,PXMIN,PXMAX,PYMIN,PYMAX,PXLEN,PYLEN

  DXLEN= XNDMAX-XNDMIN
  DYLEN= YNDMAX-YNDMIN
  DRATIO=DYLEN/DXLEN
  
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
  
  PXMIN=0.0D0
  PXMAX=25.6D0
  PYMIN=0.0D0
  PYMAX=14.0D0
  
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
END SUBROUTINE WFGWIN

!
! ----- Add. By YOKOYAMA Mar./04/2013 ----
!
!     ****** CALCULATE RANGE OF WINDOW (Y-Z PLANE) ******
!
      SUBROUTINE WFGWIN_YZ(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
!
      use wfcomm
      IMPLICIT NONE
      REAL(8):: DXLEN,DYLEN,DRATIO,PXMIN,PXMAX,PYMIN,PYMAX,PXLEN,PYLEN
      INTEGER:: NWW,NWMAX,NWYMAX,NWX,NWY,NW
!
!     XNDMAX,XNDMIN & YNDMAX,YNDMINは， WFINDX/WFVLIM で設定されている
!     XNDMAX: 節点のX座標の最大値
!     XNDMIN: 節点のX座標の最小値
!     DXLEN: 分割領域のX方向の長さ

!      DXLEN= XNDMAX-XNDMIN
      DXLEN= ZNDMAX-ZNDMIN
      DYLEN= YNDMAX-YNDMIN
      DRATIO=DYLEN/DXLEN
!
!     NWXMAXは，WFGOUTで値が代入されている．
!     要素分割を図示する場合には，WFGOUTよりも先にWFGWINがCALLされているように見えるが...
!     確認したところ，要素分割の図示の際には，NWXMAX=0であった．
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
!
      PXMIN=0.0D0
      PXMAX=25.6D0
      PYMIN=0.0D0
      PYMAX=14.0D0
!
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
    END SUBROUTINE WFGWIN_YZ
!
! ----- Mar./04/2013 -----
!
!     ****** PLOT 2D STRUCTURE OF FIELD AND POWER ******

SUBROUTINE WFGPPC(NW,NWMAX,KWD)
  
  use wfcomm
  USE libgrf
  implicit none
  integer :: KA(4,NGXM,NGYM),NWMAX,NW,NGX,NGY,ISTEP
  real(4) :: GAX(NGXM),GAY(NGYM)
  real(4) :: GXMIN,GXMAX,GYMIN,GYMAX,GPYMIN,GPYMAX
  real(4) :: GPXMIN,GPXMAX,GSXMIN,GSXMAX,GXSCAL,GSYMIN,GYSCAL
  real(4) :: GXORG,GYORG,GZMIN,GZMAX,GQZMIN,GQZMAX,GZSCAL,GZDEL,GZORG,GXPOS,GYPOS
  real(8) :: XMIN,XMAX,YMIN,YMAX,PXMIN,PXMAX,PYMIN,PYMAX
  real(8) :: DXLEN,DYLEN,DRATIO,PRATIO
  real(8) :: PXMID,PYLEN,PXLEN,PYMID,DX,DY
  CHARACTER KWD*(NCHM),KID*1
  
  IF(NWMAX.LE.0) RETURN
  
  CALL WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
  
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
!
! ----- Add. By YOKOYAMA Mar./04/2013 ----
!
      IF((KWD(4:4).EQ.'X').OR.(KWD(4:4).EQ.'Y')) THEN
         DXLEN = FRATIO*(XMAX-XMIN)
      ELSE
         DXLEN = XMAX-XMIN
      ENDIF
!
! ----- Mar./04/2013 -----
!
  DYLEN= YMAX-YMIN
  DRATIO=DYLEN/DXLEN
!
!  ---- Add. By YOKOYAMA Mar./04/2013 ----
!
      PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
! org.  PRATIO=FRATIO*(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
! ----- Mar./04/2013 -----
!
!     FUNCTION gdclip(XXX)
!        XXXの絶対値が1E-15よりも大きければ，REAL(XXX)を返す
!        XXXの絶対値が1E-15以下ならば，0.0を返す
  GXMIN=gdclip(XMIN)
  GXMAX=gdclip(XMAX)
  GYMIN=gdclip(YMIN)
  GYMAX=gdclip(YMAX)
  IF(DRATIO.GT.PRATIO) THEN
     GPYMIN=gdclip(PYMIN)+1.5
     GPYMAX=gdclip(PYMAX)
     PXMID=0.5D0*(PXMIN+PXMAX)
     PYLEN=PYMAX-PYMIN-1.5D0
     PXLEN=PYLEN/DRATIO
     GPXMIN=gdclip(PXMID-0.5D0*PXLEN)
     GPXMAX=gdclip(PXMID+0.5D0*PXLEN)
  ELSE
     GPXMIN=gdclip(PXMIN)+0.3
     GPXMAX=gdclip(PXMAX)-0.3
     PYMID=0.5D0*(PYMIN+PYMAX+1.5D0)
     PXLEN=PXMAX-PXMIN-0.6D0
     PYLEN=PXLEN*DRATIO
     GPYMIN=gdclip(PYMID-0.5D0*PYLEN)
     GPYMAX=gdclip(PYMID+0.5D0*PYLEN)
  ENDIF
  
  DX=(XMAX-XMIN)/(NGXMAX-1)
  DO NGX=1,NGXMAX
     GAX(NGX)=gdclip(XMIN+(NGX-1)*DX)
  ENDDO
  DY=(YMAX-YMIN)/(NGYMAX-1)
  DO NGY=1,NGYMAX
     GAY(NGY)=gdclip(YMIN+(NGY-1)*DY)
  ENDDO
  
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

!      IF(KID.EQ.'Z') THEN
!         CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,&
!     &               GXMIN,GXMAX,GYMIN,GYMAX)

!         CALL SETLIN(0,0,7)
!         CALL WFGBDY

!         CALL SETLIN(0,0,4)
!         CALL WFGPLA

!         CALL SETLIN(0,0,4)
!         CALL WFGANT

!         CALL OFFVEW
!      ENDIF

  CALL GMNMX2(GZ,NGXM,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
  CALL GQSCAL(GZMIN,GZMAX,GQZMIN,GQZMAX,GZSCAL)
  GZDEL=0.5*GZSCAL
  IF(GZDEL.EQ.0.0) GOTO 1000
  ISTEP=INT((GZMAX-GZMIN)/GZDEL)
  
  CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,GXMIN,GXMAX,GYMIN,GYMAX)
  
  CALL SETLIN(0,0,7)
  CALL GFRAME
  CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.1,9)
  CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.1,9)

  CALL SETLIN(0,0,7)
  IF(GZMIN*GZMAX.GT.0.0) THEN
     GZORG=GQZMIN-GZDEL
  ELSE
     GZORG=0.5*GZDEL
  ENDIF
  
!         GZORG=GQZMIN
!         CALL SETLIN(0,0,6)
!         CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,&
!     &               GZORG,GZDEL,ISTEP,0,0,KA)

  CALL SETLIN(0,0,6)
  CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX, GZORG, GZDEL,ISTEP,0,0,KA)
  CALL SETLIN(0,0,5)
  CALL CONTP2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,-GZORG,-GZDEL,ISTEP,0,2,KA)
!      ENDIF

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
END SUBROUTINE WFGPPC

!     ****** PLOT 2D STRUCTURE OF FIELD AND POWER ******

SUBROUTINE WFGPFC(NW,NWMAX,KWD)

  use wfcomm
  USE libgrf
  implicit none
  integer,parameter :: NSTEPM=101
  integer,parameter :: NRGBA=5
  integer,parameter :: NRGBB=7
  integer :: NWMAX,NW,NGX,NGY,ISTEP,I
  real(8) :: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
  real(8) :: DXLEN,DYLEN,DRATIO,PRATIO
  real(8) :: PXMID,PYLEN,PXLEN,PYMID,DX,DY
  real(4) :: GXMIN,GXMAX,GYMIN,GYMAX,GPYMIN,GPYMAX,GPXMIN,GPXMAX
  real(4) :: GSXMAX,GSXMIN,GXSCAL,GSYMIN,GYSCAL,GXORG,GYORG,GZMIN,GZMAX,GZA
  real(4) :: GDZ,GFACT,GXPOS,GYPOS
  real(4) :: GAX(NGXM),GAY(NGYM),GZL(NSTEPM),GRGBL(3,0:NSTEPM)
  real(4) :: GRGBA(3,NRGBA),GLA(NRGBA),GRGBB(3,NRGBB),GLB(NRGBB)
  CHARACTER KWD*(NCHM),KID*1
  DATA GRGBA/0.0,0.0,1.0,&
       &     0.0,1.0,1.0,&
       &     1.0,1.0,1.0,&
       &     1.0,1.0,0.0,&
       &     1.0,0.0,0.0/
  DATA GLA/0.0,0.40,0.5,0.60,1.0/
  DATA GRGBB/1.0,1.0,1.0,&
       &     0.0,1.0,1.0,&
       &     0.0,0.0,1.0,&
       &     0.0,1.0,0.0,&
       &     1.0,0.0,0.0,&
       &     1.0,1.0,0.0,&
       &     1.0,1.0,1.0/
  DATA GLB/0.0,0.15,0.3,0.5,0.7,0.85,1.0/
  
  IF(NWMAX.LE.0) RETURN
  
  CALL WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
  
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
  
  DXLEN= XMAX-XMIN
  DYLEN= YMAX-YMIN
  DRATIO=DYLEN/DXLEN
  PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
  GXMIN=gdclip(XMIN)
  GXMAX=gdclip(XMAX)
  GYMIN=gdclip(YMIN)
  GYMAX=gdclip(YMAX)
  IF(DRATIO.GT.PRATIO) THEN
     GPYMIN=gdclip(PYMIN)+1.5
     GPYMAX=gdclip(PYMAX)
     PXMID=0.5D0*(PXMIN+PXMAX)
     PYLEN=PYMAX-PYMIN-1.5D0
     PXLEN=PYLEN/DRATIO
     GPXMIN=gdclip(PXMID-0.5D0*PXLEN)
     GPXMAX=gdclip(PXMID+0.5D0*PXLEN)
  ELSE
     GPXMIN=gdclip(PXMIN)+0.3
     GPXMAX=gdclip(PXMAX)-0.3
     PYMID=0.5D0*(PYMIN+PYMAX+1.5D0)
     PXLEN=PXMAX-PXMIN-0.6D0
     PYLEN=PXLEN*DRATIO
     GPYMIN=gdclip(PYMID-0.5D0*PYLEN)
     GPYMAX=gdclip(PYMID+0.5D0*PYLEN)
  ENDIF
  
  DX=(XMAX-XMIN)/(NGXMAX-1)
  DO NGX=1,NGXMAX
     GAX(NGX)=gdclip(XMIN+(NGX-1)*DX)
  ENDDO
  DY=(YMAX-YMIN)/(NGYMAX-1)
  DO NGY=1,NGYMAX
     GAY(NGY)=gdclip(YMIN+(NGY-1)*DY)
  ENDDO
  
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
  
  ISTEP=50
  
  CALL GMNMX2(GZ,NGXM,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
  IF(GZMIN*GZMAX.LT.0.D0) THEN
     GZA=MAX(ABS(GZMAX),ABS(GZMIN))
     GDZ=2*GZA/ISTEP
     DO I=1,ISTEP
        GZL(I)=GDZ*(I-0.5)-GZA
     ENDDO
     
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
  
  CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,GXMIN,GXMAX,GYMIN,GYMAX)
  
  CALL CONTF2(GZ,GAX,GAY,NGXM,NGXMAX,NGYMAX,GZL,GRGBL,ISTEP,0)
  CALL SETLIN(0,0,7)
  CALL GFRAME
  CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.1,9)
  CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.1,9)
  
  CALL RGBBAR(GPXMAX+0.2,GPXMAX+0.5,GPYMIN,GPYMAX,GRGBL,ISTEP+1,1)
  
!      IF(KID.EQ.'Z') THEN
!         CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,&
!     &               GXMIN,GXMAX,GYMIN,GYMAX)

!         CALL SETLIN(0,0,7)
!         CALL WFGBDY

!         CALL SETLIN(0,0,4)
!         CALL WFGPLA

!         CALL SETLIN(0,0,4)
!         CALL WFGANT

!         CALL OFFVEW
!      ENDIF

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
END SUBROUTINE WFGPFC

!     ****** PLOT 2D STRUCTURE OF FIELD AND POWER ******

SUBROUTINE WFGPBC(NW,NWMAX,KWD)
  
  use wfcomm
  USE libgrf
  implicit none
  integer :: KA(8,NGXM,NGYM),NWMAX,NW,NGX,NGY
  real(4) :: GAX(NGXM),GAY(NGYM),GXMIN,GXMAX,GYMIN,GYMAX
  real(4) :: GPYMIN,GPYMAX,GPXMIN,GPXMAX
  real(8) :: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
  real(8) :: DXLEN,DYLEN,DRATIO,PRATIO
  real(8) :: PXMID,PYLEN,PXLEN,PYMID,DX,DY
  real(4) :: GSXMIN,GSXMAX,GXSCAL,GSYMIN,GYSCAL,GXORG,GYORG,GZMIN,GZMAX,GSZMIN
  real(4) :: GSZMAX,GZSCAL,GZORG,GXL,GYL,GZL,GPHI,GTHETA,GRADIUS,GXPOS,GYPOS,GZA
  EXTERNAL R2W2B
  CHARACTER KWD*(NCHM),KID*1
  
  IF(NWMAX.LE.0) RETURN
  
  CALL WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
  
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
  
  DXLEN= XMAX-XMIN
  DYLEN= YMAX-YMIN
  DRATIO=DYLEN/DXLEN
  PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
  GXMIN=gdclip(XMIN)
  GXMAX=gdclip(XMAX)
  GYMIN=gdclip(YMIN)
  GYMAX=gdclip(YMAX)
  IF(DRATIO.GT.PRATIO) THEN
     GPYMIN=gdclip(PYMIN)+1.5
     GPYMAX=gdclip(PYMAX)
     PXMID=0.5D0*(PXMIN+PXMAX)
     PYLEN=PYMAX-PYMIN-1.5D0
     PXLEN=PYLEN/DRATIO
     GPXMIN=gdclip(PXMID-0.5D0*PXLEN)
     GPXMAX=gdclip(PXMID+0.5D0*PXLEN)
  ELSE
     GPXMIN=gdclip(PXMIN)+0.3
     GPXMAX=gdclip(PXMAX)-0.3
     PYMID=0.5D0*(PYMIN+PYMAX+1.5D0)
     PXLEN=PXMAX-PXMIN-0.6D0
     PYLEN=PXLEN*DRATIO
     GPYMIN=gdclip(PYMID-0.5D0*PYLEN)
     GPYMAX=gdclip(PYMID+0.5D0*PYLEN)
  ENDIF
  
  DX=(XMAX-XMIN)/(NGXMAX-1)
  DO NGX=1,NGXMAX
     GAX(NGX)=gdclip(XMIN+(NGX-1)*DX)
  ENDDO
  DY=(YMAX-YMIN)/(NGYMAX-1)
  DO NGY=1,NGYMAX
     GAY(NGY)=gdclip(YMIN+(NGY-1)*DY)
  ENDDO
  
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
  
  CALL GMNMX2(GZ,NGXM,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
  GZMAX=MAX(ABS(GZMIN),ABS(GZMAX))
  GZMIN=-GZMAX
  CALL GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GZSCAL)
  GZORG=0.0
  
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
  
  CALL GDEFIN3D(GPXMIN,GPXMAX,GPYMIN,GPYMAX,GXL,GYL,GZL)
  CALL GVIEW3D(GPHI,GTHETA,GRADIUS,0.9,1,0.5*(GXMIN+GXMAX),0.5*(GYMIN+GYMAX),0.0)
  CALL GDATA3D1(GZ,NGXM,NGXMAX,NGYMAX,GXMIN,GXMAX,GYMIN,GYMAX,GZMIN,GZMAX)
  
  CALL CPLOT3D1(1,R2W2B)
  CALL CONTQ3D1(GZMIN,0.1*(GZMAX-GZMIN),11,0,0,KA,R2W2B,0)
  
  CALL GAXIS3D(0)
  CALL GSCALE3DX(GXORG,GXSCAL,0.2,2)
  CALL GSCALE3DY(GYORG,GYSCAL,0.2,2)
  CALL GSCALE3DZ(GZORG,GZSCAL,0.2,2)
  CALL SETCHS(0.2,0.0)
  CALL GVALUE3DX(GXORG,2.*GXSCAL,-6,2)
  CALL GVALUE3DY(GYORG,2.*GYSCAL,-6,2)
  CALL GVALUE3DZ(GZORG,2.*GZSCAL,-2,0)
  
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
END SUBROUTINE WFGPBC

!     ****** DRAW 1D PROFILE ******

SUBROUTINE WFGPFR(NW,NWMAX,KWD)
  
  use wfcomm
  USE libgrf
  implicit none
  integer :: IPAT(3),NWMAX,NW,NGMAX,NG
  real(8) :: PXMIN,PYMIN,PXMAX,PYMAX,PRATIO,PXLEN,PYLEN,PXMID
  real(4) :: GPXMIN
  real(4) :: GXMIN1,GXMAX1,GYMIN1,GYMAX1,GXMIN,GXMAX,GXSTEP,GYMIN,GYMAX,GYSTEP
  real(4) :: GXORG,GYORG,GPXMAX,GPYMIN,GPYMAX,GXPOS,GYPOS
  CHARACTER KWD*(NCHM)
  DATA IPAT/0,2,4/
  
  IF(NWMAX.LE.0) RETURN
  
  CALL WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)

  !     ------ DO NOT ALLOW LANDSCAPE GRAPH ------
  
  PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-2.5D0)
  IF(PRATIO.LT.1.D0) THEN
     PYLEN=PYMAX-PYMIN-1.5D0
     PXLEN=PYLEN+2.5D0
     PXMID=0.5D0*(PXMIN+PXMAX)
     PXMIN=PXMID-0.5D0*PXLEN
     PXMAX=PXMID+0.5D0*PXLEN
  ENDIF
  
  GPXMIN=gdclip(PXMIN)+2.2
  GPXMAX=gdclip(PXMAX)-0.3
  GPYMIN=gdclip(PYMIN)+1.5
  GPYMAX=gdclip(PYMAX)
  
  IF(KWD(1:1).EQ.'E'.OR.&
  &  KWD(1:1).EQ.'D'.OR.&
  &  KWD(1:1).EQ.'B'.OR.&
  &  KWD(1:1).EQ.'A'    ) THEN
     NGMAX=3
  ELSE
     NGMAX=1
  ENDIF
  
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
  
  CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,&
       &      GXMIN ,GXMAX ,GYMIN ,GYMAX  )
  
  CALL SETLIN(0,0,4)
  CALL SETCHS(0.25,0.0)
  CALL GFRAME
  CALL GSCALE(GXORG,GXSTEP,0.0,0.0,0.2,9)
  CALL GSCALE(0.0,0.0,GYORG,GYSTEP,0.2,9)
  IF(GYMIN*GYMAX.LT.0.0) THEN
     CALL GSCALE(0.0,0.0,GYORG,100*GYSTEP,0.0,0)
  ENDIF
  CALL GVALUE(GXORG,4*GXSTEP,0.0,0.0,NGSLEN(4*GXSTEP))
  CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGSLEN(2*GYSTEP))
  
  DO NG=1,NGMAX
     CALL SETLIN(0,0,8-NG)
     CALL GPLOTP(GX,GV(1,NG),1,NGVMAX,1,0,0,IPAT(NG))
  ENDDO
  
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
END SUBROUTINE WFGPFR

!     ****** DRAW PARAMETER ON GRAPHIC SCREEN ******

SUBROUTINE WFGPRM
  
  use wfcomm
  implicit none
  integer :: NA,NB,L,I,NK,NM,NS
  real(8) :: REST(NAM),REAT(NAM)
  real(4) :: GXMIN,GYMAX,GRCHH,GDX,GDY,GXL,GYL
  real(8) :: SRFR(NMDM,NBM),SRFI(NMDM,NBM),SRFL(NMDM,NBM)
  
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
  
  GXMIN=0.0
  GYMAX=18.2
  GRCHH=0.30
  CALL SETCHS(GRCHH,0.0)
  GDX=15.*GRCHH
  GDY=-1.5*GRCHH
  GXL=GXMIN
  GYL=GYMAX
  
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
  
  GXL=GXL+GDX
  CALL MOVE(GXL,GYL)
  CALL TEXT('BB  =',5)
  CALL NUMBD(BB  ,'(F7.3)',7)
  
  GXL=GXL+GDX
  CALL MOVE(GXL,GYL)
  CALL TEXT('M=',2)
  CALL NUMBI(MODELG,'(I2)',2)
  CALL NUMBI(MODELB,'(I2)',2)
  CALL NUMBI(MODELD,'(I2)',2)
  CALL NUMBI(MODELP,'(I2)',2)
  CALL NUMBI(MODELS,'(I2)',2)
  CALL NUMBI(MODELN,'(I2)',2)
  
  GXL=GXMIN
  GYL=GYL+GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('NNMAX=',6)
  CALL NUMBI(NNMAX,'(I6)',6)
  
  GXL=GXL+GDX
  CALL MOVE(GXL,GYL)
  CALL TEXT('XYZ MAX=',8)
  CALL NUMBD(XNDMAX,'(F7.3)',7)
  CALL NUMBD(YNDMAX,'(F7.3)',7)
  CALL NUMBD(ZNDMAX,'(F7.3)',7)
  
  GXL=GXMIN
  GYL=GYL+GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('NEMAX=',6)
  CALL NUMBI(NEMAX,'(I6)',6)
  
  GXL=GXL+GDX
  CALL MOVE(GXL,GYL)
  CALL TEXT('XYZ MIN=',8)
  CALL NUMBD(XNDMIN,'(F7.3)',7)
  CALL NUMBD(YNDMIN,'(F7.3)',7)
  CALL NUMBD(ZNDMIN,'(F7.3)',7)

!      GXL=GXL+GDX
!      CALL MOVE(GXL,GYL)
!      CALL TEXT('P=',2)
!      CALL NUMBD(TSPWR,'(1PE10.3)',10)

  GXL=GXMIN
  GYL=GYL+GDY
  
  CALL MOVE(GXL,GYL)
  CALL TEXT(' NK',3)
  CALL TEXT(' NM',3)
  CALL TEXT(' PABS     ',10)
  CALL TEXT(' NK',3)
  CALL TEXT(' NM',3)
  CALL TEXT(' PABS     ',10)
  
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
  
  IF(NSMAX.GT.0) THEN
     GXL=GXMIN
     GYL=GYL+GDY
     CALL MOVE(GXL,GYL)
     CALL TEXT(' NS',3)
     CALL TEXT('    PA    ',10)
     CALL TEXT('    PZ    ',10)
     CALL TEXT('    PN    ',10)
     CALL TEXT(' PZCL ',6)
     CALL TEXT('   PABS   ',10)
     
     DO NS=1,NSMAX
        GXL=GXMIN
        GYL=GYL+GDY
        CALL MOVE(GXL,GYL)
        CALL NUMBI(NS,'(I3)',3)
        CALL NUMBD(PA(NS),'(1PE10.2)',10)
        CALL NUMBD(PZ(NS),'(1PE10.2)',10)
        CALL NUMBD(PN(NS),'(1PE10.2)',10)
        CALL NUMBD(PZCL(NS),'(F6.2)',6)
        CALL NUMBD(PABSS(NS),'(1PE10.2)',10)
     ENDDO
  ENDIF
  
  GXL=GXMIN+45.*GRCHH
  GYL=GYMAX+GDY
  IF(NAMAX.GT.0) THEN
     CALL MOVE(GXL,GYL)
     CALL TEXT('IJ',2)
     CALL TEXT('  AJ ',5)
     CALL TEXT(' PHASE ',7)
     CALL TEXT('     R     ',11)
     CALL TEXT('     X     ',11)
     
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
  
  RETURN
END SUBROUTINE WFGPRM

!     ****** Draw Vessel Boundary ******

SUBROUTINE WFGBDY

  RETURN
END SUBROUTINE WFGBDY

!     ****** Draw Plasma Boundary ******

SUBROUTINE WFGPLA
  
  use wfcomm
  USE libgrf
  implicit none
  integer :: NPMAX,I
  real(8) :: DTH,THETA
  real(4) :: GXL,GYL
  
  NPMAX=100
  DTH=2.D0*PI/NPMAX
  GXL=gdclip(RA)
  GYL=0.0
  CALL MOVE(GXL,GYL)
  DO I=1,NPMAX
     THETA=DTH*I
     GXL=gdclip(RA*COS(THETA))
     GYL=gdclip(RA*SIN(THETA))
     CALL DRAW(GXL,GYL)
  END DO
  RETURN
END SUBROUTINE WFGPLA

!     ****** Draw Antenna Path ******

SUBROUTINE WFGANT

  use wfcomm
  USE libgrf
  implicit none
  integer :: NA,I
  real(4) :: GXL,GYL

  IF(NDRAWA.EQ.0) THEN
     DO NA=1,NAMAX
        GXL=gdclip(XJ0(1,NA))
        GYL=gdclip(YJ0(1,NA))
        CALL MOVE(GXL,GYL)
        DO I=2,JNUM0(NA)
           GXL=gdclip(XJ0(I,NA))
           GYL=gdclip(YJ0(I,NA))
           CALL DRAW(GXL,GYL)
        END DO
     END DO
  ELSE
     DO NA=1,NAMAX
        GXL=gdclip(XJ(1,NA))
        GYL=gdclip(YJ(1,NA))
        CALL MOVE(GXL,GYL)
        DO I=2,JNUM(NA)
           GXL=gdclip(XJ(I,NA))
           GYL=gdclip(YJ(I,NA))
           CALL DRAW(GXL,GYL)
        END DO
     END DO
  ENDIF
  RETURN
END SUBROUTINE WFGANT

!     ****** Draw Element Data ******

SUBROUTINE WFGELM

  use wfcomm
  USE libgrf
  implicit none
  integer :: I1(4),I2(4),I3(4),IE,K,IN1,IN2,IN3,IEL,IN,INL,NN
  real(4) :: GX1,GX2,GX3,GY1,GY2,GY3,GXC,GYC
  real(8) :: EPS,ANTZP,DELTZ,DZMIN,ZAXIS
  real(8) :: XT(3),YT(3),ZT(3)!test
  DATA I1,I2,I3/1,2,3,4, 2,3,4,1, 3,4,1,2/
  DATA EPS/1.D-6/
  
  CALL SETCHR(0.2,0.15,0.2,0.,-30.)

  ANTZP=0.0D0
  DZMIN=5.D0
  DO NN=1,NNMAX
     DELTZ=ABS(ZND(NN)-ANTZP)
     IF(DELTZ.LT.DZMIN) THEN
        DZMIN=DELTZ
        ZAXIS=ZND(NN)
     ENDIF
  ENDDO
  
  DO IE=1,NEMAX
     DO K=1,4
        IN1=NDELM(I1(K),IE)
        IN2=NDELM(I2(K),IE)
        IN3=NDELM(I3(K),IE)
!
! ----- Add. By YOKOYAMA Mar./05/2013 ----
!
!        if(MODELG.eq.1) then
!           CALL RCTORT(XND(IN1),YND(IN1),ZND(IN1),XT(1),YT(1),ZT(1))
!           CALL RCTORT(XND(IN2),YND(IN2),ZND(IN2),XT(2),YT(2),ZT(2))
!           CALL RCTORT(XND(IN3),YND(IN3),ZND(IN3),XT(3),YT(3),ZT(3))
!           IF(ZT(1).GT.ZNDMAX-EPS.AND.&
!            & ZT(2).GT.ZNDMAX-EPS.AND.&
!            & ZT(3).GT.ZNDMAX-EPS) THEN
!
         IF(ABS(ZND(IN1)-ZAXIS).LT.EPS.AND.&
     &      ABS(ZND(IN2)-ZAXIS).LT.EPS.AND.&
     &      ABS(ZND(IN3)-ZAXIS).LT.EPS) THEN
!
!  ---- Mar./05/2013 -----
              GX1=gdclip(XT(1))
              GY1=gdclip(YT(1))
              GX2=gdclip(XT(2))
              GY2=gdclip(YT(2))
              GX3=gdclip(XT(3))
              GY3=gdclip(YT(3))
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
              
              IF(NDRAWD.GE.2) THEN
                 GXC=(GX1+GX2+GX3)/3.
                 GYC=(GY1+GY2+GY3)/3.
                 IEL=IE
                 CALL GNUMBI(GXC,GYC,IEL,2)
              ENDIF
        else
           IF(ZND(IN1).GT.ZNDMAX-EPS.AND.&
            & ZND(IN2).GT.ZNDMAX-EPS.AND.&
            & ZND(IN3).GT.ZNDMAX-EPS) THEN
              GX1=gdclip(XND(IN1))
              GY1=gdclip(YND(IN1))
              GX2=gdclip(XND(IN2))
              GY2=gdclip(YND(IN2))
              GX3=gdclip(XND(IN3))
              GY3=gdclip(YND(IN3))
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
           END IF
           IF(NDRAWD.GE.2) THEN
              GXC=(GX1+GX2+GX3)/3.
              GYC=(GY1+GY2+GY3)/3.
              IEL=IE
              CALL GNUMBI(GXC,GYC,IEL,2)
           ENDIF
        end if
     ENDDO
  ENDDO
  
  IF(NDRAWD.GE.3) THEN
     CALL SETCHS(0.2,0.)
     DO IN=1,NNMAX
        GX1=gdclip(XND(IN))
        GY1=gdclip(YND(IN))
        INL=IN
        CALL GNUMBI(GX1,GY1,INL,0)
     END DO
  ENDIF
  
  RETURN
END SUBROUTINE WFGELM

!     ****** Draw Element Data ******

SUBROUTINE WFGNAS(ID)
  
  use wfcomm
  USE libgrf
  implicit none
  integer :: ID
  real(8) :: PXMIN,PXMAX,PYMIN,PYMAX
  real(4) :: GPXMIN,GPXMAX,GPYMIN,GPYMAX,GYMIN,GYMAX,GXMIN,GXMAX
  real(8) :: DXLEN,DYLEN,DRATIO,PRATIO,XMID,XLEN,YMID,YLEN
  
  CALL WFGWIN(1,1,PXMIN,PXMAX,PYMIN,PYMAX)
  GPXMIN=gdclip(PXMIN)
  GPXMAX=gdclip(PXMAX)
  GPYMIN=gdclip(PYMIN)+1.0
  GPYMAX=gdclip(PYMAX)+1.0
  
  DXLEN= XNDMAX-XNDMIN+0.5D0*(YNDMAX-YNDMIN)
  DYLEN= ZNDMAX-ZNDMIN+0.5D0*(YNDMAX-YNDMIN)
  DRATIO=DYLEN/DXLEN
  PRATIO=(PYMAX-PYMIN)/(PXMAX-PXMIN)
  IF(DRATIO.GT.PRATIO) THEN
     GYMIN=gdclip(ZNDMIN+0.5D0*YNDMIN)
     GYMAX=gdclip(ZNDMAX+0.5D0*YNDMAX)
     XMID=0.5D0*(XNDMIN+XNDMAX+0.5D0*(YNDMAX+YNDMIN))
     XLEN=DYLEN/PRATIO
     GXMIN=gdclip(XMID-0.5D0*XLEN)
     GXMAX=gdclip(XMID+0.5D0*XLEN)
  ELSE
     GXMIN=gdclip(XNDMIN+0.5D0*YNDMIN)
     GXMAX=gdclip(XNDMAX+0.5D0*YNDMIN)
     YMID=0.5D0*(ZNDMIN+ZNDMAX+0.5D0*(YNDMAX+YNDMIN))
     YLEN=DXLEN*PRATIO
     GYMIN=gdclip(YMID-0.5D0*YLEN)
     GYMAX=gdclip(YMID+0.5D0*YLEN)
  ENDIF
  
  CALL PAGES
  CALL WFPRME
  CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,&
       &      GXMIN ,GXMAX ,GYMIN ,GYMAX )
  !     &            GXMIN,GXMAX,GYMIN,GYMAX)
  !     &     -0.03,-0.01,-0.005,0.005)
  !     &     -0.03,0.03,-0.02,0.02)
  
  CALL SETLIN(0,0,7)
  CALL WFGELM3(ID)
  
  CALL PAGEE
  RETURN
END SUBROUTINE WFGNAS

!     ****** Draw 3D Element Data ******

SUBROUTINE WFGELM3(ID)

  use wfcomm
  USE libgrf
  implicit none
  integer :: ID,IE,IN1,IN2,IN3,IN4,IN,INL,NSF,ND1,ND2,ND3,NN,KA,IEL
  real(4) :: GX1,GY1,GX2,GY2,GX3,GY3,GX4,GY4,GXC,GYC

!     ***** DRAW EDGE *****

  IF(ID.EQ.0) THEN
     DO IE=1,NEMAX
        IN1=NDELM(1,IE)
        IN2=NDELM(2,IE)
        IN3=NDELM(3,IE)
        IN4=NDELM(4,IE)
        GX1=gdclip(XND(IN1))
        GY1=gdclip(ZND(IN1)+0.5D0*YND(IN1))
        GX2=gdclip(XND(IN2))
        GY2=gdclip(ZND(IN2)+0.5D0*YND(IN2))
        GX3=gdclip(XND(IN3))
        GY3=gdclip(ZND(IN3)+0.5D0*YND(IN3))
        GX4=gdclip(XND(IN4))
        GY4=gdclip(ZND(IN4)+0.5D0*YND(IN4))
        CALL MOVE(GX1,GY1)
        CALL DRAW(GX2,GY2)
        CALL DRAW(GX3,GY3)
        CALL DRAW(GX4,GY4)
        CALL DRAW(GX2,GY2)
        CALL MOVE(GX4,GY4)
        CALL DRAW(GX1,GY1)
        CALL DRAW(GX3,GY3)
        
        IF(NDRAWD.GE.2) THEN
           GXC=(GX1+GX2+GX3)/3.
           GYC=(GY1+GY2+GY3)/3.
           IEL=IE
           CALL GNUMBI(GXC,GYC,IEL,2)
        ENDIF
     ENDDO
     
     IF(NDRAWD.GE.3) THEN
        CALL SETCHS(0.2,0.)
        DO IN=1,NNMAX
           GX1=gdclip(XND(IN))
           GY1=gdclip(YND(IN))
           INL=IN
           CALL GNUMBI(GX1,GY1,INL,0)
        ENDDO
     ENDIF

!     ***** DRAW SURFACE EDGE *****

  ELSEIF(ID.EQ.1) THEN
     CALL SETRGB(0.0,0.0,1.0)
     DO NSF=1,NSFMAX
        ND1=NDSRF(1,NSF)
        ND2=NDSRF(2,NSF)
        ND3=NDSRF(3,NSF)
        GX1=gdclip(XND(ND1))
        GY1=gdclip(ZND(ND1)+0.5D0*YND(ND1))
        GX2=gdclip(XND(ND2))
        GY2=gdclip(ZND(ND2)+0.5D0*YND(ND2))
        GX3=gdclip(XND(ND3))
        GY3=gdclip(ZND(ND3)+0.5D0*YND(ND3))
        CALl MOVE(GX1,GY1)
        CALl DRAW(GX2,GY2)
        CALl DRAW(GX3,GY3)
        CALl DRAW(GX1,GY1)
     ENDDO

!     ***** DRAW NODES *****

  ELSEIF(ID.EQ.2) THEN
     CALL SETMKS(1,0.5)
     DO NN=1,NNMAX
        GX1=gdclip(XND(NN))
        GY1=gdclip(ZND(NN)+0.5D0*YND(NN))
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
     
  ELSEIF(ID.EQ.3) THEN
     CALL SETMKS(1,0.5)
     DO NN=1,NNMAX
        GX1=gdclip(XND(NN))
        GY1=gdclip(ZND(NN)+0.5D0*YND(NN))
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
     
  ELSEIF(ID.EQ.4) THEN
     CALL SETMKS(1,0.5)
     DO NN=1,NNMAX
        GX1=gdclip(XND(NN))
        GY1=gdclip(ZND(NN)+0.5D0*YND(NN))
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
     
  ELSEIF(ID.EQ.5) THEN
     CALL SETMKS(1,0.5)
     DO NN=1,NNMAX
        GX1=gdclip(XND(NN))
        GY1=gdclip(ZND(NN)+0.5D0*YND(NN))
        KA=KANOD(NN)
        IF(KA.EQ.2) THEN
           CALL SETRGB(1.0,0.0,0.0)
           CALL MARK(GX1,GY1)
        ENDIF
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE WFGELM3


!     ---- Draw Element Data (ADD. BY YAMA 30/SEP./2008) ----

      SUBROUTINE WFGELMA

      use wfcomm
      USE libgrf  
      implicit none
      INTEGER,DIMENSION(4):: I1=(/1,2,3,4/)
      INTEGER,DIMENSION(4):: I2=(/2,3,4,1/)
      INTEGER,DIMENSION(4):: I3=(/3,4,1,2/)
      REAL(8),PARAMETER:: EPS=1.D-6
      REAL(8):: ANTZP,DZMIN,ZAXIS,DELTZ
      REAL(4):: GX1,GY1,GX2,GY2,GX3,GY3
      REAL(4):: GXC,GYC
      INTEGER:: NN,IE,K,IN1,IN2,IN3,IEL,IN,INL

      CALL SETCHR(0.2,0.15,0.2,0.,-30.)

!      ANTZP = 1.715D0
!      ANTZP = 2.055D0
!      ANTZP = 5.6D0
      ANTZP = ZJ0(1,1)
      DZMIN = 1.D0
      DO NN=1,NNMAX
         DELTZ = ABS(ZND(NN) - ANTZP)
         IF(DELTZ.LT.DZMIN) THEN
            DZMIN = DELTZ
            ZAXIS = ZND(NN)
         ENDIF
      ENDDO

      DO IE=1,NEMAX
         DO K=1,4
         IN1=NDELM(I1(K),IE)
         IN2=NDELM(I2(K),IE)
         IN3=NDELM(I3(K),IE)
         IF(ABS(ZND(IN1)-ZAXIS).LT.EPS.AND. &
            ABS(ZND(IN2)-ZAXIS).LT.EPS.AND. &
            ABS(ZND(IN3)-ZAXIS).LT.EPS) THEN
            GX1=gdclip(XND(IN1))
            GY1=gdclip(YND(IN1))
            GX2=gdclip(XND(IN2))
            GY2=gdclip(YND(IN2))
            GX3=gdclip(XND(IN3))
            GY3=gdclip(YND(IN3))
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

            IF(NDRAWD.GE.2) THEN
               GXC=(GX1+GX2+GX3)/3.
               GYC=(GY1+GY2+GY3)/3.
               IEL=IE
               CALL GNUMBI(GXC,GYC,IEL,2)
            ENDIF
         ENDIF
         ENDDO
      ENDDO

      IF(NDRAWD.GE.3) THEN
         CALL SETCHS(0.2,0.)
         DO IN=1,NNMAX
            GX1=gdclip(XND(IN))
            GY1=gdclip(YND(IN))
            INL=IN
            CALL GNUMBI(GX1,GY1,INL,0)
         END DO
      ENDIF

      RETURN
    END SUBROUTINE WFGELMA

!
!
! ----- Add. By YOKOYAMA Mar./05/2013 ----
!
!     ****** Draw Element Data (ON Y-Z PLANE) ******

      SUBROUTINE WFGELM_YZ

      use wfcomm
      USE libgrf  
      implicit none
      INTEGER,DIMENSION(4),PARAMETER:: I1=(/1,2,3,4/)
      INTEGER,DIMENSION(4),PARAMETER:: I2=(/2,3,4,1/)
      INTEGER,DIMENSION(4),PARAMETER:: I3=(/3,4,1,2/)
      REAL(8),PARAMETER:: EPS=1.D-6
      INTEGER:: IE,KK,IIN1,IIN2,IIN3,IN1,IN2,IN3,K,IZLN,IEL,IN,INL
      REAL(8):: XNDS,YNDS,ZAVG
      REAL(4):: GX1,GY1,GX2,GY2,GX3,GY3,GXC,GYC

      CALL SETCHR(0.2,0.15,0.2,0.,-30.)

      DO IE=1,NEMAX

         XNDS=0.D0
         YNDS=0.D0
         DO KK=1,4
            IIN1=NDELM(I1(KK),IE)
            IIN2=NDELM(I2(KK),IE)
            IIN3=NDELM(I3(KK),IE)
            XNDS=XND(IIN1)+XND(IIN2)+XND(IIN3)
            YNDS=YND(IIN1)+YND(IIN2)+YND(IIN3)
         ENDDO

!         IF(((YNDS.GT.0.D0).AND.(XNDS.GT.0.D0)).OR.
!     &      ((YNDS.LT.0.D0).AND.(XNDS.LT.0.D0))) THEN
         IF(XNDS.GE.0.D0) THEN

         DO K=1,4
         IN1=NDELM(I1(K),IE)
         IN2=NDELM(I2(K),IE)
         IN3=NDELM(I3(K),IE)

!  --  重ね描きの回避 --
         IZLN=0
         ZAVG=(ZND(IN1)+ZND(IN2)+ZND(IN3))/3.D0
         IF(ZND(IN1).LT.ZAVG) IZLN=IZLN+1
         IF(ZND(IN2).LT.ZAVG) IZLN=IZLN+1
         IF(ZND(IN3).LT.ZAVG) IZLN=IZLN+1

         IF(ABS(XND(IN1)).LT.EPS.AND. &
            ABS(XND(IN2)).LT.EPS.AND. &
            ABS(XND(IN3)).LT.EPS.AND. &
            IZLN.EQ.1) THEN


            GX1=gdclip(ZND(IN1))
            GY1=gdclip(YND(IN1))
            GX2=gdclip(ZND(IN2))
            GY2=gdclip(YND(IN2))
            GX3=gdclip(ZND(IN3))
            GY3=gdclip(YND(IN3))
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

            IF(NDRAWD.GE.2) THEN
               GXC=(GX1+GX2+GX3)/3.
               GYC=(GY1+GY2+GY3)/3.
               IEL=IE
               CALL GNUMBI(GXC,GYC,IEL,2)
            ENDIF
         ENDIF
         ENDDO

         ENDIF

      ENDDO

      IF(NDRAWD.GE.3) THEN
         CALL SETCHS(0.2,0.)
         DO IN=1,NNMAX
            GX1=gdclip(ZND(IN))
            GY1=gdclip(YND(IN))
            INL=IN
            CALL GNUMBI(GX1,GY1,INL,0)
         END DO
      ENDIF

      RETURN
    END SUBROUTINE WFGELM_YZ

!
!     ****** Draw Element Data (ON X-Z PLANE) ******
!
      SUBROUTINE WFGELM_XZ
!
      use wfcomm
      USE libgrf  
      implicit none
      INTEGER,DIMENSION(4),PARAMETER:: I1=(/1,2,3,4/)
      INTEGER,DIMENSION(4),PARAMETER:: I2=(/2,3,4,1/)
      INTEGER,DIMENSION(4),PARAMETER:: I3=(/3,4,1,2/)
      REAL(8),PARAMETER:: EPS=1.D-6
      INTEGER:: IE,KK,IIN1,IIN2,IIN3,IN1,IN2,IN3,K,IZLN,IEL,IN,INL
      REAL(8):: XNDS,YNDS,ZAVG
      REAL(4):: GX1,GY1,GX2,GY2,GX3,GY3,GXC,GYC
!
      CALL SETCHR(0.2,0.15,0.2,0.,-30.)
!
!      DO IE=1,NEMAX,2
      DO IE=1,NEMAX

         YNDS=0.D0
         DO KK=1,4
            IIN1=NDELM(I1(KK),IE)
            IIN2=NDELM(I2(KK),IE)
            IIN3=NDELM(I3(KK),IE)
            YNDS=YND(IIN1)+YND(IIN2)+YND(IIN3)
         ENDDO

         IF(YNDS.GT.0.D0) THEN

         DO K=1,4
         IN1=NDELM(I1(K),IE)
         IN2=NDELM(I2(K),IE)
         IN3=NDELM(I3(K),IE)

!  --  重ね描きの回避 --
         IZLN=0
         ZAVG=(ZND(IN1)+ZND(IN2)+ZND(IN3))/3.D0
         IF(ZND(IN1).LT.ZAVG) IZLN=IZLN+1
         IF(ZND(IN2).LT.ZAVG) IZLN=IZLN+1
         IF(ZND(IN3).LT.ZAVG) IZLN=IZLN+1

         IF(ABS(YND(IN1)).LT.EPS.AND.&
     &      ABS(YND(IN2)).LT.EPS.AND.&
     &      ABS(YND(IN3)).LT.EPS.AND.&
     &      IZLN.EQ.2) THEN


!            GX1=gdclip(XND(IN1))
            GX1=gdclip(ZND(IN1))
            GY1=gdclip(XND(IN1))
!            GX2=gdclip(XND(IN2))
            GX2=gdclip(ZND(IN2))
            GY2=gdclip(XND(IN2))
!            GX3=gdclip(XND(IN3))
            GX3=gdclip(ZND(IN3))
            GY3=gdclip(XND(IN3))
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
!
            IF(NDRAWD.GE.2) THEN
               GXC=(GX1+GX2+GX3)/3.
               GYC=(GY1+GY2+GY3)/3.
               IEL=IE
               CALL GNUMBI(GXC,GYC,IEL,2)
            ENDIF
         ENDIF
         ENDDO

         ENDIF

      ENDDO
!
      IF(NDRAWD.GE.3) THEN
         CALL SETCHS(0.2,0.)
         DO 200 IN=1,NNMAX
!            GX1=gdclip(XND(IN))
            GX1=gdclip(ZND(IN))
            GY1=gdclip(YND(IN))
            INL=IN
            CALL GNUMBI(GX1,GY1,INL,0)
  200    CONTINUE
      ENDIF
!
      RETURN
      END SUBROUTINE WFGELM_XZ

!
! ----- Mar./05/2013 -----
!
!     ****** Draw Element Data ******

SUBROUTINE WFGDIV

  use wfcomm
  USE libgrf
  implicit none
  real(8) :: PXMIN,PXMAX,PYMIN,PYMAX
  real(4) :: GPXMIN,GPXMAX,GPYMIN,GPYMAX,GYMIN,GYMAX,GXMIN,GXMAX
  real(8) :: DXLEN,DYLEN,DRATIO,PRATIO,XMID,XLEN,YMID,YLEN
  
  CALL WFGWIN(1,1,PXMIN,PXMAX,PYMIN,PYMAX)
  GPXMIN=gdclip(PXMIN)
  GPXMAX=gdclip(PXMAX)
  GPYMIN=gdclip(PYMIN)+1.0
  GPYMAX=gdclip(PYMAX)+1.0
  
  DXLEN= XNDMAX-XNDMIN
  DYLEN= YNDMAX-YNDMIN
  DRATIO=DYLEN/DXLEN
  PRATIO=(PYMAX-PYMIN)/(PXMAX-PXMIN)

  IF(DRATIO.GT.PRATIO) THEN
     GYMIN=gdclip(YNDMIN)
     GYMAX=gdclip(YNDMAX)
     XMID=0.5D0*(XNDMIN+XNDMAX)
     XLEN=DYLEN/PRATIO
     GXMIN=gdclip(XMID-0.5D0*XLEN)
     GXMAX=gdclip(XMID+0.5D0*XLEN)
  ELSE
     GXMIN=gdclip(XNDMIN)
     GXMAX=gdclip(XNDMAX)
     YMID=0.5D0*(YNDMIN+YNDMAX)
     YLEN=DXLEN*PRATIO
     GYMIN=gdclip(YMID-0.5D0*YLEN)
     GYMAX=gdclip(YMID+0.5D0*YLEN)
  ENDIF
  
  CALL PAGES
  CALL WFPRME
  CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,GXMIN,GXMAX,GYMIN,GYMAX)
  
  CALL SETLIN(0,0,7)
  IF(NDRAWD.EQ.0) THEN
     CALL WFGBDY
  ELSE
     CALL WFGELM
  ENDIF
  
  CALL PAGEE
  RETURN
END SUBROUTINE WFGDIV

!
!
!
! ----- Add. By YOKOYAMA Mar./05/2013 ----
!
!     ****** Draw Element Data (Y-Z PLANE) ******
!
      SUBROUTINE WFGDIV_YZ
!
      use wfcomm
      USE libgrf  
      implicit none
      REAL(8):: PXMIN,PXMAX,PYMIN,PYMAX,DXLEN,DYLEN,DRATIO,PRATIO
      REAL(8):: XMID,XLEN,YMID,YLEN
      REAL(4):: GPXMIN,GPXMAX,GPYMIN,GPYMAX,GYMIN,GYMAX,GXMIN,GXMAX
!
!      CALL WFGWIN(1,1,PXMIN,PXMAX,PYMIN,PYMAX)
      CALL WFGWIN_YZ(1,1,PXMIN,PXMAX,PYMIN,PYMAX)
      GPXMIN=gdclip(PXMIN)
      GPXMAX=gdclip(PXMAX)
      GPYMIN=gdclip(PYMIN)+1.0
      GPYMAX=gdclip(PYMAX)+1.0
!
!      DXLEN= XNDMAX-XNDMIN
      DXLEN= ZNDMAX-ZNDMIN
      DYLEN= YNDMAX-YNDMIN
      DRATIO=DYLEN/DXLEN
      PRATIO=(PYMAX-PYMIN)/(PXMAX-PXMIN)
      IF(DRATIO.GT.PRATIO) THEN
         GYMIN=gdclip(YNDMIN)
         GYMAX=gdclip(YNDMAX)
!         XMID=0.5D0*(XNDMIN+XNDMAX)
         XMID=0.5D0*(ZNDMIN+ZNDMAX)
         XLEN=DYLEN/PRATIO
         GXMIN=gdclip(XMID-0.5D0*XLEN)
         GXMAX=gdclip(XMID+0.5D0*XLEN)
      ELSE
!         GXMIN=gdclip(XNDMIN)
!         GXMAX=gdclip(XNDMAX)
         GXMIN=gdclip(ZNDMIN)
         GXMAX=gdclip(ZNDMAX)
         YMID=0.5D0*(YNDMIN+YNDMAX)
         YLEN=DXLEN*PRATIO
         GYMIN=gdclip(YMID-0.5D0*YLEN)
         GYMAX=gdclip(YMID+0.5D0*YLEN)
      ENDIF
!
      CALL PAGES
      CALL WFPRME
      CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,&
     &            GXMIN,GXMAX,GYMIN,GYMAX)
!
      CALL SETLIN(0,0,7)
      IF(NDRAWD.EQ.0) THEN
         CALL WFGBDY
      ELSE
!         CALL WFGELM
         CALL WFGELM_YZ
      ENDIF
!
      CALL PAGEE
      RETURN
      END SUBROUTINE WFGDIV_YZ

!     ****** Draw Element Data (X-Z PLANE) ******
!
      SUBROUTINE WFGDIV_XZ
!
      use wfcomm
      USE libgrf  
      implicit none
      REAL(8):: PXMIN,PXMAX,PYMIN,PYMAX,DXLEN,DYLEN,DRATIO,PRATIO
      REAL(8):: XMID,XLEN,YMID,YLEN
      REAL(4):: GPXMIN,GPXMAX,GPYMIN,GPYMAX,GYMIN,GYMAX,GXMIN,GXMAX
!
!      CALL WFGWIN(1,1,PXMIN,PXMAX,PYMIN,PYMAX)
      CALL WFGWIN_YZ(1,1,PXMIN,PXMAX,PYMIN,PYMAX)
      GPXMIN=gdclip(PXMIN)
      GPXMAX=gdclip(PXMAX)
      GPYMIN=gdclip(PYMIN)+1.0
      GPYMAX=gdclip(PYMAX)+1.0
!
!      DXLEN= XNDMAX-XNDMIN
      DXLEN= ZNDMAX-ZNDMIN
      DYLEN= YNDMAX-YNDMIN
      DRATIO=DYLEN/DXLEN
      PRATIO=(PYMAX-PYMIN)/(PXMAX-PXMIN)
      IF(DRATIO.GT.PRATIO) THEN
         GYMIN=gdclip(YNDMIN)
         GYMAX=gdclip(YNDMAX)
!         XMID=0.5D0*(XNDMIN+XNDMAX)
         XMID=0.5D0*(ZNDMIN+ZNDMAX)
         XLEN=DYLEN/PRATIO
         GXMIN=gdclip(XMID-0.5D0*XLEN)
         GXMAX=gdclip(XMID+0.5D0*XLEN)
      ELSE
!         GXMIN=gdclip(XNDMIN)
!         GXMAX=gdclip(XNDMAX)
         GXMIN=gdclip(ZNDMIN)
         GXMAX=gdclip(ZNDMAX)
         YMID=0.5D0*(YNDMIN+YNDMAX)
         YLEN=DXLEN*PRATIO
         GYMIN=gdclip(YMID-0.5D0*YLEN)
         GYMAX=gdclip(YMID+0.5D0*YLEN)
      ENDIF
!
      CALL PAGES
      CALL WFPRME
      CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,&
     &            GXMIN,GXMAX,GYMIN,GYMAX)
!
      CALL SETLIN(0,0,7)
      IF(NDRAWD.EQ.0) THEN
         CALL WFGBDY
      ELSE
!         CALL WFGELM
         CALL WFGELM_XZ
!         CALL WFGELM3(0)
      ENDIF
      CALL PAGEE
      RETURN
      END SUBROUTINE WFGDIV_XZ
!
! ----- Mar./05/2013 -----
!
!     ****** Draw Element Paramters ******

SUBROUTINE WFPRME

  use wfcomm
  implicit none
  real(4) :: GXMIN,GYMAX,GDY,GXL,GYL
  
  GXMIN=20.
  GYMAX=17.
  CALL SETCHS(0.3,0.)
  GDY=1.5*0.3
  
  GXL=GXMIN
  GYL=GYMAX
  CALL MOVE(GXL,GYL)
  CALL TEXT('NNMAX=',6)
  CALL NUMBI(NNMAX,'(I5)',5)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('NEMAX=',6)
  CALL NUMBI(NEMAX,'(I5)',5)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('NBMAX=',6)
  CALL NUMBI(NBMAX,'(I5)',5)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('MBND =',6)
  CALL NUMBI(MBND,'(I5)',5)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('MLEN =',6)
  CALL NUMBI(MLEN,'(I5)',5)
  RETURN
END SUBROUTINE WFPRME

!     ****** Draw Antenna ******

SUBROUTINE WFPLTA

  use wfcomm
  USE libgrf
  implicit none
  integer :: NTEMP
  real(8) :: PXMIN,PXMAX,PYMIN,PYMAX
  real(4) :: GYMIN,GYMAX,GXMIN,GXMAX,GPXMIN,GPXMAX,GPYMIN,GPYMAX
  real(8) :: DXLEN,DYLEN,DRATIO,PRATIO,XMID,XLEN,YMID,YLEN

  CALL WFGWIN(1,1,PXMIN,PXMAX,PYMIN,PYMAX)
  GPXMIN=gdclip(PXMIN)
  GPXMAX=gdclip(PXMAX)
  GPYMIN=gdclip(PYMIN)+1.0
  GPYMAX=gdclip(PYMAX)+1.0
  
  DXLEN= XNDMAX-XNDMIN
  DYLEN= YNDMAX-YNDMIN
  DRATIO=DYLEN/DXLEN
  PRATIO=(PYMAX-PYMIN)/(PXMAX-PXMIN)
  IF(DRATIO.GT.PRATIO) THEN
     GYMIN=gdclip(YNDMIN)
     GYMAX=gdclip(YNDMAX)
     XMID=0.5D0*(XNDMIN+XNDMAX)
     XLEN=DYLEN/PRATIO
     GXMIN=gdclip(XMID-0.5D0*XLEN)
     GXMAX=gdclip(XMID+0.5D0*XLEN)
  ELSE
     GXMIN=gdclip(XNDMIN)
     GXMAX=gdclip(XNDMAX)
     YMID=0.5D0*(YNDMIN+YNDMAX)
     YLEN=DXLEN*PRATIO
     GYMIN=gdclip(YMID-0.5D0*YLEN)
     GYMAX=gdclip(YMID+0.5D0*YLEN)
  ENDIF
  
  CALL PAGES
  CALL WFPRMJ
  CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,GXMIN,GXMAX,GYMIN,GYMAX)
  
  CALL SETLIN(0,0,4)
  IF(NDRAWA.LE.1) THEN
     CALL WFGBDY
  ELSE
     NTEMP=NDRAWD
     NDRAWD=NDRAWA-1
!
! ----- Add. By YOKOYAMA Mar./05/2013 ----
!         CALL WFGELM
! ----- Mar./05/2013 -----
         CALL WFGELMA
     CALL WFGELM
     NDRAWD=NTEMP
  ENDIF
  
  CALL SETLIN(0,0,5)
  CALL WFGPLA
  
  CALL SETLIN(0,0,6)
  CALL WFGANT
  
  CALL PAGEE
  RETURN
END SUBROUTINE WFPLTA

!     ****** Draw Antenna Paramters ******

SUBROUTINE WFPRMJ

  use wfcomm
  implicit none
  integer :: NA
  real(4) :: GXMIN,GYMAX,GDY,GXL,GYL

  GXMIN=20.
  GYMAX=17.
  CALL SETCHS(0.3,0.)
  GDY=1.5*0.3
  
  GXL=GXMIN
  GYL=GYMAX
  CALL MOVE(GXL,GYL)
  CALL TEXT('NNMAX=',6)
  CALL NUMBI(NNMAX,'(I5)',5)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('NEMAX=',6)
  CALL NUMBI(NEMAX,'(I5)',5)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('NBMAX=',6)
  CALL NUMBI(NBMAX,'(I5)',5)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('MBND =',6)
  CALL NUMBI(MBND,'(I5)',5)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('MLEN =',6)
  CALL NUMBI(MLEN,'(I5)',5)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('JNUM =',6)
  GXL=GXL+6*0.3
  DO NA=1,NAMAX
     CALL MOVE(GXL,GYL)
     CALL NUMBI(JNUM(NA),'(I5)',5)
     GYL=GYL-GDY
  END DO
  RETURN
END SUBROUTINE WFPRMJ
