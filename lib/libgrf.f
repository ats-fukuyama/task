C     $Id$
C
C     ***** DRAW 1D GRAPH *****
C
      SUBROUTINE GRF1D(NGP,GX,GY,NXM,NXMAX,NGMAX,STR,MODE)
C
C     +++++ MODE=0 X:LINEAR  Y:LINEAR
C     +++++ MODE=1 X:LOG     Y:LINEAR
C     +++++ MODE=2 X:LINEAR  Y:LOG
C     +++++ MODE=3 X:LOG     Y:LOG
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GP(4,0:4),GX(NXMAX),GY(NXM,NGMAX)
      DIMENSION GRGB(3,5),IPAT(5)
      CHARACTER STR*(*)
C
      DATA GP/ 3.0,23.0, 2.0,17.0,
     &         3.0,11.0,11.0,17.0,
     &         3.0,11.0, 2.0, 8.0,
     &        15.0,23.0,11.0,17.0,
     &        15.0,23.0, 2.0, 8.0/
      DATA IPAT/0,2,3,4,6/
      DATA GRGB/0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,0.0,1.0,
     &          0.0,1.0,0.0, 1.0,0.0,1.0/
C
      IF(NGP.LT.0.OR.NGP.GT.4) THEN
         WRITE(6,*) 'XX GRF1D: INVALID NGP'
         RETURN
      ENDIF
C
      CALL GRFUT1(GX,NXMAX,GXMIN,GXMAX)
      CALL GRFUT2(GY,NXM,NXMAX,NGMAX,GYMIN,GYMAX)
C
      CALL GRFUT3(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSTEP,GXORG)
      CALL GRFUT3(GYMIN,GYMAX,GSYMIN,GSYMAX,GYSTEP,GYORG)
C
      NLMAX=5
      CALL GRF1DX(GP(1,NGP),
     &            GSXMIN,GSXMAX,GXSTEP,GXORG,
     &            GSYMIN,GSYMAX,GYSTEP,GYORG,
     &            GX,GY,NXM,NXMAX,NGMAX,
     &            GRGB,IPAT,NLMAX,STR,MODE)
      RETURN
      END
C
C     ***** DRAW 1D GRAPH *****
C
      SUBROUTINE GRF1DX(GP,
     &                  GSXMIN,GSXMAX,GXSTEP,GXORG,
     &                  GSYMIN,GSYMAX,GYSTEP,GYORG,
     &                  GX,GY,NXM,NXMAX,NGMAX,
     &                  GRGB,IPAT,NLMAX,STR,MODE)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GP(4),GX(NXMAX),GY(NXM,NGMAX)
      DIMENSION GRGB(3,NLMAX),IPAT(NLMAX)
      CHARACTER STR*(*),KT*80,KDL*1
C
      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.035)
      CALL SETRGB(0.0,0.0,0.0)
C
      KDL=STR(1:1)
      I=2
    1 IF(STR(I:I).EQ.KDL.OR.I.EQ.LEN(STR)) GOTO 2
         KT(I-1:I-1)=STR(I:I)
         I=I+1
      GOTO 1
C
    2 CALL MOVE(GP(1),GP(4)+0.3)
      CALL TEXT(KT,I-2)
C
      CALL GDEFIN(GP(1),GP(2),GP(3),GP(4),
     &            GSXMIN,GSXMAX,GSYMIN,GSYMAX)
      CALL GFRAME
C
      IF(MODE.EQ.0.OR.MODE.EQ.2) THEN
         CALL GSCALE(0.0,2*(GSXMAX-GSXMIN),0.0,0.0,0.0,0)
         CALL GSCALE(GXORG,GXSTEP,0.0,0.0,0.1,9)
         CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
      ELSE
         CALL GSCALL(GXORG,9,0.0,0,0.1,9)
         CALL GVALUL(GXORG,1,0.0,0,NGULEN(2*GXSTEP))
      ENDIF
      IF(MODE.EQ.0.OR.MODE.EQ.1) THEN
         CALL GSCALE(0.0,0.0,0.0,2*(GSYMAX-GSYMIN),0.0,0)
         CALL GSCALE(0.0,0.0,GYORG,GYSTEP,0.1,9)
         CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))
      ELSE
         CALL GSCALL(0.0,0,GYORG,9,0.1,9)
         CALL GVALUL(0.0,0,GYORG,1,NGULEN(2*GYSTEP))
      ENDIF
C
      DO NG=1,NGMAX
         NL=MOD(NG-1,NLMAX)+1
         CALL SETRGB(GRGB(1,NL),GRGB(2,NL),GRGB(3,NL))
         CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,IPAT(NL))
      ENDDO
C
      CALL SETRGB(0.0,0.0,0.0)
      RETURN
      END
C
C     ***** DRAW 2D GRAPH *****
C
      SUBROUTINE GRF2DC(NGP,GX,GY,GZ,NXM,NXMAX,NYMAX,STR,KA,MODE)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GP(4,0:4),GX(NXMAX),GY(NYMAX),GZ(NXM,NYMAX)
      DIMENSION KA(8,NXM,NYMAX)
      CHARACTER STR*(*)
      DATA GP/ 3.0,23.0, 2.0,17.0,
     &         3.0,11.0,11.0,17.0,
     &         3.0,11.0, 2.0, 8.0,
     &        15.0,23.0,11.0,17.0,
     &        15.0,23.0, 2.0, 8.0/
C
      IF(NGP.LT.0.OR.NGP.GT.4) THEN
         WRITE(6,*) 'XX GRF2DC: INVALID NGP'
         RETURN
      ENDIF
C
      CALL GRFUT1(GX,NXMAX,GXMIN,GXMAX)
      CALL GRFUT1(GY,NYMAX,GYMIN,GYMAX)
      CALL GRFUT2(GZ,NXM,NXMAX,NYMAX,GZMIN,GZMAX)
C
      CALL GRFUT3(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSTEP,GXORG)
      CALL GRFUT3(GYMIN,GYMAX,GSYMIN,GSYMAX,GYSTEP,GYORG)
      CALL GRFUT3(GZMIN,GZMAX,GSZMIN,GSZMAX,GZSTEP,GZORG)
      GZSTEP=0.5*GZSTEP
C
      CALL GRF2DCX(GP,
     &             GXMIN,GXMAX,GXSTEP,GXORG,
     &             GYMIN,GYMAX,GYSTEP,GYORG,
     &             GZMIN,GZMAX,GZSTEP,GZORG,
     &             GZ,NXM,NXMAX,NYMAX,KA,MODE,STR)
C
      RETURN
      END
C
C     ***** DRAW 2D GRAPH *****
C
      SUBROUTINE GRF2DCX(GP,
     &                   GXMIN,GXMAX,GXSTEP,GXORG,
     &                   GYMIN,GYMAX,GYSTEP,GYORG,
     &                   GZMIN,GZMAX,GZSTEP,GZORG,
     &                   GZ,NXM,NXMAX,NYMAX,
     &                   KA,MODE,STR)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GP(4),GZ(NXM,NYMAX)
      DIMENSION KA(8,NXM,NYMAX)
      CHARACTER STR*(*),KT*80,KDL*1
C
      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.035)
      CALL SETRGB(0.0,0.0,0.0)
C
      KDL=STR(1:1)
      I=2
    1 IF(STR(I:I).EQ.KDL.OR.I.EQ.LEN(STR)) GOTO 2
         KT(I-1:I-1)=STR(I:I)
         I=I+1
      GOTO 1
C
    2 CALL MOVE(GP(1),GP(4)+0.3)
      CALL TEXT(KT,I-2)
C
      CALL GDEFIN(GP(1),GP(2),GP(3),GP(4),
     &            GXMIN,GXMAX,GYMIN,GYMAX)
      CALL GFRAME
C
      IF(MODE.EQ.0.OR.MODE.EQ.2) THEN
         CALL GSCALE(0.0,2*(GXMAX-GXMIN),0.0,0.0,0.0,0)
         CALL GSCALE(GXORG,GXSTEP,0.0,0.0,0.1,9)
         CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
      ELSE
         CALL GSCALL(GXORG,9,0.0,0,0.1,9)
         CALL GVALUL(GXORG,1,0.0,0,NGULEN(2*GXSTEP))
      ENDIF
      IF(MODE.EQ.0.OR.MODE.EQ.1) THEN
         CALL GSCALE(0.0,0.0,0.0,2*(GYMAX-GYMIN),0.0,0)
         CALL GSCALE(0.0,0.0,GYORG,GYSTEP,0.1,9)
         CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGULEN(2*GYSTEP))
      ELSE
         CALL GSCALL(0.0,0,GYORG,9,0.1,9)
         CALL GVALUL(0.0,0,GYORG,1,NGULEN(2*GYSTEP))
      ENDIF
C
      IF(GZMIN*GZMAX.GE.0.D0) THEN
         NSTEP=NINT(ABS((GZMAX-GZMIN)/GZSTEP))+1
         IF(GZORG.GE.0.D0) THEN
            CALL SETRGB(1.0,0.0,0.0)
            CALL CONTQ1(GZ,NXM,NXMAX,NYMAX,GZORG,GZSTEP,NSTEP,0,0,KA)
         ELSE
            CALL SETRGB(0.0,0.0,1.0)
            CALL CONTQ1(GZ,NXM,NXMAX,NYMAX,GZORG,GZSTEP,NSTEP,0,3,KA)
         ENDIF
      ELSE
         CALL SETRGB(1.0,0.0,0.0)
         NSTEP=NINT(ABS(GZMAX/GZSTEP))+1
         CALL CONTQ1(GZ,NXM,NXMAX,NYMAX, 0.5*GZSTEP, GZSTEP,
     &               NSTEP,0,0,KA)
C         CALL SETRGB(0.0,1.0,0.0)
C         CALL CONTQ1(GZ,NXM,NXMAX,NYMAX, 0.0,    GZSTEP,
C     &               1,0,4,KA)
         CALL SETRGB(0.0,0.0,1.0)
         NSTEP=NINT(ABS(GZMIN/GZSTEP))+1
         CALL CONTQ1(GZ,NXM,NXMAX,NYMAX,-0.5*GZSTEP,-GZSTEP,
     &               NSTEP,0,3,KA)
      ENDIF
      CALL SETRGB(0.0,0.0,0.0)
      CALL MOVE(0.5*(GP(1)+GP(2))-0.2*8,GP(3)-1.5)
      CALL TEXT('STEP =',6)
      CALL NUMBR(GZSTEP,'(1PE10.2)',10)
C
      CALL SETRGB(0.0,0.0,0.0)
      RETURN
      END
C
C     ***** DRAW 2D GRAPH *****
C
      SUBROUTINE GRF2DB(NGP,GX,GY,GZ,NXM,NXMAX,NYMAX,STR)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GP(4,0:4),GX(NXMAX),GY(NYMAX),GZ(NXM,NYMAX)
      DIMENSION GP3D(6)
      CHARACTER STR*(*)
      DATA GP/ 3.0,23.0, 2.0,17.0,
     &         3.0,11.0,11.0,17.0,
     &         3.0,11.0, 2.0, 8.0,
     &        15.0,23.0,11.0,17.0,
     &        15.0,23.0, 2.0, 8.0/
C
      IF(NGP.LT.0.OR.NGP.GT.4) THEN
         WRITE(6,*) 'XX GRF2DB: INVALID NGP'
         RETURN
      ENDIF
C
      GP3D(1)=10.0*1.5
      GP3D(2)=20.0*1.5
      GP3D(3)=10.0*1.5
      GP3D(4)=-60.0
      GP3D(5)=65.0
      GP3D(6)=100.0
C
      CALL GRFUT1(GX,NXMAX,GXMIN,GXMAX)
      CALL GRFUT1(GY,NYMAX,GYMIN,GYMAX)
      CALL GRFUT2(GZ,NXM,NXMAX,NYMAX,GZMIN,GZMAX)
C
      CALL GRFUT3(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSTEP,GXORG)
      CALL GRFUT3(GYMIN,GYMAX,GSYMIN,GSYMAX,GYSTEP,GYORG)
      CALL GRFUT3(GZMIN,GZMAX,GSZMIN,GSZMAX,GZSTEP,GZORG)
C
      CALL GRF2DBX(GP,
     &             GXMIN,GXMAX,GXSTEP,GXORG,
     &             GYMIN,GYMAX,GYSTEP,GYORG,
     &             GZMIN,GZMAX,GZSTEP,GZORG,
     &             GZ,NXM,NXMAX,NYMAX,GP3D,STR)
      RETURN
      END
C
C     ***** DRAW 2D GRAPH *****
C
      SUBROUTINE GRF2DBX(GP,
     &                   GXMIN,GXMAX,GXSTEP,GXORG,
     &                   GYMIN,GYMAX,GYSTEP,GYORG,
     &                   GZMIN,GZMAX,GZSTEP,GZORG,
     &                   GZ,NXM,NXMAX,NYMAX,
     &                   GP3D,STR)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GP(4),GZ(NXM,NYMAX)
      DIMENSION GP3D(6)
      CHARACTER STR*(*),KT*80,KDL*1
C
      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.035)
      CALL SETRGB(0.0,0.0,0.0)
C
      KDL=STR(1:1)
      I=2
    1 IF(STR(I:I).EQ.KDL.OR.I.EQ.LEN(STR)) GOTO 2
         KT(I-1:I-1)=STR(I:I)
         I=I+1
      GOTO 1
C
    2 CALL MOVE(GP(1),GP(4)+0.3)
      CALL TEXT(KT,I-2)
C
      GXL    =GP3D(1)
      GYL    =GP3D(2)
      GZL    =GP3D(3)
      GPHI   =GP3D(4)
      GTHETA =GP3D(5)
      GRADIUS=GP3D(6)
      GOX=0.5*(GXMIN+GXMAX)
      GOY=0.5*(GYMIN+GYMAX)
      GOZ=0.5*(GZMIN+GZMAX)
      CALL GDEFIN3D(GP(1),GP(2),GP(3),GP(4),GXL,GYL,GZL)
      CALL GVIEW3D(GPHI,GTHETA,GRADIUS,1.0,1,GOX,GOY,GOZ)
      CALL GDATA3D1(GZ,NXM,NXMAX,NYMAX,
     &              GXMIN,GXMAX,GYMIN,GYMAX,GZMIN,GZMAX)
C
      CALL SETCHS(0.2,0.0)
      CALL GSCALE3DX(GXORG,GXSTEP,0.3,0)
      CALL GSCALE3DY(GYORG,GYSTEP,0.3,0)
      CALL GSCALE3DZ(GZORG,GZSTEP,0.3,0)
      CALL GVALUE3DX(GXORG,2*GXSTEP,1,NGULEN(2*GXSTEP))
      CALL GVALUE3DY(GYORG,2*GYSTEP,1,NGULEN(2*GYSTEP))
      CALL GVALUE3DZ(GZORG,2*GZSTEP,2,NGULEN(2*GZSTEP))
C
      CALL PERSE3D(3,1)
      CALL GAXIS3D(0)
C
      CALL SETRGB(0.0,0.0,0.0)
      RETURN
      END
C
C     ***** EVALUATE GXMIN,GXMAX,GYMIN,GYMAX *****
C
      SUBROUTINE GRFUT1(GX,NXMAX,GXMIN,GXMAX)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GX(NXMAX)
C
      CALL GMNMX1(GX,1,NXMAX,1,GXMIN,GXMAX)
      IF(ABS(GXMAX-GXMIN).LT.1.D-6) THEN
         GXMIN=GXMIN-0.999E-6
         GXMAX=GXMAX+1.000E-6
      ENDIF
C
      RETURN
      END
C
C     ***** EVALUATE GXMIN,GXMAX,GYMIN,GYMAX *****
C
      SUBROUTINE GRFUT2(GY,NXM,NXMAX,NYMAX,GYMIN,GYMAX)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GY(NXM,NYMAX)
C
      CALL GMNMX2(GY,NXM,1,NXMAX,1,1,NYMAX,1,GYMIN,GYMAX)
      IF(ABS(GYMAX-GYMIN).LT.1.D-6) THEN
         GYMIN=GYMIN-0.999E-6
         GYMAX=GYMAX+1.000E-6
      ENDIF
C
      RETURN
      END
C
C     ***** EVALUATE GSMIN,GSMAX,GSTEP,GORG *****
C
      SUBROUTINE GRFUT3(GMIN,GMAX,GSMIN,GSMAX,GSTEP,GORG)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      CALL GQSCAL(GMIN,GMAX,GSMIN,GSMAX,GSTEP)
C
      IF(GMIN*GMAX.LE.0.0) THEN
         GORG=0.0
      ELSE
         GORG=(NINT(GSMIN/(2*GSTEP))+1)*2*GSTEP
      ENDIF
      RETURN
      END
