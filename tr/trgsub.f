C     $Id$
C
C     ***********************************************************
C
C           SUBPROGRAM FOR 1D PROFILE
C
C                   MODE =  0  : Y=0 INCLUDED
C                          +1  : USING YMIN/YMAX
C                          +2  : LINE PATTERN CHANGE
C                          +4  : INPUT YMIN/YMAX
C                          +8  : LOG SCALE
C
C     ***********************************************************
C
      SUBROUTINE TRGR1D(GX1,GX2,GY1,GY2,GX,GY,NXM,NXMAX,NGMAX,STR,MODE)
C
C      IMPLICIT REAL (KIND=8) (A-F,H,O-Z)
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GX(NXM),GY(NXM,NGMAX),IPAT(5)
      CHARACTER STR*(*),KT*80,KDL*1
      DATA IPAT/0,2,3,4,6/
C
      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.035)
      CALL SETRGB(0.0,0.0,0.0)
      CALL SETLIN(-1,-1,7)
      KDL=STR(1:1)
      I=2
    1 IF(STR(I:I).EQ.KDL.OR.I.EQ.LEN(STR)) GOTO 2
         KT(I-1:I-1)=STR(I:I)
         I=I+1
      GOTO 1
C
    2 CALL MOVE(GX1,GY2+0.3)
      CALL TEXT(KT,I-2)
C
      CALL GMNMX2(GY,NXM,1,NXMAX,1,1,NGMAX,1,GYMIN,GYMAX)
      IF(ABS(GYMAX-GYMIN).LT.1.D-6) THEN
         GYMIN=GYMIN-0.999E-6
         GYMAX=GYMAX+1.000E-6
      ENDIF
C
      IF(MOD(MODE,2).EQ.0) THEN
         IF(GYMIN.GE.0.0) THEN
            GYMIN=0.0
         ELSEIF(GYMAX.LE.0.0) THEN
            GYMAX=0.0
         ENDIF
      ENDIF
C
      CALL GMNMX1(GX,1,NXMAX,1,GXMIN,GXMAX)
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
C
      IF(MOD(MODE,2).EQ.0) THEN
         IF(GYMIN.GE.0.0) THEN
            GSYMIN=0.0
         ELSEIF(GYMAX.LE.0.0) THEN
            GSYMAX=0.0
         ENDIF
      ENDIF
      GYMIN=GSYMIN
      GYMAX=GSYMAX
      IF(MOD(MODE/4,2).EQ.1) THEN
         CALL CHMODE
  701    WRITE(6,*) '## TRGR : XMIN,XMAX,YMIN,YMAX = ',
     &              GXMIN,GXMAX,GYMIN,GYMAX
         READ(5,*,ERR=701,END=900) GXMIN,GXMAX,GYMIN,GYMAX
         CALL GRMODE
      ENDIF
      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSTEPY)
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
      CALL GDEFIN(GX1,GX2,GY1,GY2,
     &            GSXMIN,GSXMAX,GSYMIN,GSYMAX)
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(-1,-1,7)
      CALL GFRAME
      CALL GSCALE(GXORG,GSTEPX,0.0,0.0,0.1,9)
      CALL GSCALE(0.0,0.0,0.0,2*(GSYMAX-GSYMIN),0.0,0)
      CALL GVALUE(GXORG,2*GSTEPX,0.0,0.0,NGVLEN(2*GSTEPX))
      IF(MOD(MODE/8,2).EQ.1) THEN
         CALL GSCALL(0.0,0,GYORG,9,0.1,9)
         CALL GVALUL(0.0,0,GYORG,1,NGVLEN(2*GSTEPY))
      ELSE
         CALL GSCALE(0.0,0.0,GYORG,GSTEPY,0.1,9)
         CALL GVALUE(0.0,0.0,GYORG,2*GSTEPY,NGVLEN(2*GSTEPY))
      ENDIF
C
      DO 10 NG=1,NGMAX
         CALL SETLIN(-1,-1,7-MOD(NG-1,5))
         IF(MOD(MODE/2,2).EQ.0) THEN
            CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,0)
C            IF(NG.EQ.5) THEN
C            CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,2)
C            ELSE
C            CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,0)
C            ENDIF
         ELSE
            CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,IPAT(MOD(NG-1,5)+1))
         ENDIF
   10 CONTINUE
C
      CALL SETLIN(-1,-1,7)
  900 RETURN
      END
C
C     ***********************************************************
C
C           SUBPROGRAM FOR 1D PROFILE
C
C                   MODE =  0  : Y=0 INCLUDED
C                          +1  : USING YMIN/YMAX
C                          +2  : LINE PATTERN CHANGE
C                          +4  : INPUT YMIN/YMAX
C                          +8  : LOG SCALE
C
C     ***********************************************************
C
      SUBROUTINE TRGR1DX(GX1,GX2,GY1,GY2,GX,GY,NXM,NXMAX,NGMAX,STR,MODE)
C
C      IMPLICIT REAL (KIND=8) (A-F,H,O-Z)
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GX(NXM,NGMAX),GY(NXM,NGMAX),IPAT(5)
      CHARACTER STR*(*),KT*80,KDL*1
      DATA IPAT/0,2,3,4,6/
C
      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.035)
      CALL SETLIN(-1,-1,7)
      KDL=STR(1:1)
      I=2
    1 IF(STR(I:I).EQ.KDL.OR.I.EQ.LEN(STR)) GOTO 2
         KT(I-1:I-1)=STR(I:I)
         I=I+1
      GOTO 1
C
    2 CALL MOVE(GX1,GY2+0.3)
      CALL TEXT(KT,I-2)
C
      CALL GMNMX2(GY,NXM,1,NXMAX,1,1,NGMAX,1,GYMIN,GYMAX)
      IF(ABS(GYMAX-GYMIN).LT.1.D-6) THEN
         GYMIN=GYMIN-0.999E-6
         GYMAX=GYMAX+1.000E-6
      ENDIF
C
      CALL GMNMX2(GX,NXM,1,NXMAX,1,1,NGMAX,1,GXMIN,GXMAX)
      IF(ABS(GXMAX-GXMIN).LT.1.D-6) THEN
         GXMIN=GXMIN-0.999E-6
         GXMAX=GXMAX+1.000E-6
      ENDIF
C
      IF(MOD(MODE,2).EQ.0) THEN
         IF(GYMIN.GE.0.0) THEN
            GYMIN=0.0
         ELSEIF(GYMAX.LE.0.0) THEN
            GYMAX=0.0
         ENDIF
      ENDIF
C
      IF(MOD(MODE,2).EQ.0) THEN
         IF(GXMIN.GE.0.0) THEN
            GXMIN=0.0
         ELSEIF(GXMAX.LE.0.0) THEN
            GXMAX=0.0
         ENDIF
      ENDIF
C
      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSTEPY)
C
      IF(MOD(MODE,2).EQ.0) THEN
         IF(GYMIN.GE.0.0) THEN
            GSYMIN=0.0
         ELSEIF(GYMAX.LE.0.0) THEN
            GSYMAX=0.0
         ENDIF
      ENDIF
C
      IF(MOD(MODE,2).EQ.0) THEN
         IF(GXMIN.GE.0.0) THEN
            GSXMIN=0.0
         ELSEIF(GXMAX.LE.0.0) THEN
            GSXMAX=0.0
         ENDIF
      ENDIF
C
      GYMIN=GSYMIN
      GYMAX=GSYMAX
      GXMIN=GSXMIN
      GXMAX=GSXMAX
      IF(MOD(MODE/4,2).EQ.1) THEN
         CALL CHMODE
  701    WRITE(6,*) '## TRGR : XMIN,XMAX,YMIN,YMAX = ',
     &              GXMIN,GXMAX,GYMIN,GYMAX
         READ(5,*,ERR=701,END=900) GXMIN,GXMAX,GYMIN,GYMAX
         CALL GRMODE
      ENDIF
      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSTEPY)
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
      CALL GDEFIN(GX1,GX2,GY1,GY2,
     &            GSXMIN,GSXMAX,GSYMIN,GSYMAX)
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(-1,-1,7)
      CALL GFRAME
      CALL GSCALE(GXORG,GSTEPX,0.0,0.0,0.1,9)
      CALL GSCALE(0.0,0.0,0.0,2*(GSYMAX-GSYMIN),0.0,0)
      CALL GVALUE(GXORG,2*GSTEPX,0.0,0.0,NGVLEN(2*GSTEPX))
      IF(MOD(MODE/8,2).EQ.1) THEN
         CALL GSCALL(0.0,0,GYORG,9,0.1,9)
         CALL GVALUL(0.0,0,GYORG,1,NGVLEN(2*GSTEPY))
      ELSE
         CALL GSCALE(0.0,0.0,GYORG,GSTEPY,0.1,9)
         CALL GVALUE(0.0,0.0,GYORG,2*GSTEPY,NGVLEN(2*GSTEPY))
      ENDIF
C
      DO 10 NG=1,NGMAX
         CALL SETLIN(-1,-1,7-MOD(NG-1,5))
         IF(MOD(MODE/2,2).EQ.0) THEN
            CALL GPLOTP(GX(1,NG),GY(1,NG),1,NXMAX,1,0,0,0)
C            IF(NG.EQ.5) THEN
C            CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,2)
C            ELSE
C            CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,0)
C            ENDIF
         ELSE
            CALL GPLOTP(GX(1,NG),GY(1,NG),1,NXMAX,1,
     &           0,0,IPAT(MOD(NG-1,5)+1))
         ENDIF
   10 CONTINUE
C
      CALL SETLIN(-1,-1,7)
  900 RETURN
      END
C
C     ***********************************************************
C
C           SUBPROGRAM FOR 2D CONTOUR
C
C     ***********************************************************
C
      SUBROUTINE TRGR2C(GX1,GX2,GY1,GY2,GXMIN,GXMAX,GYMIN,GYMAX,
     &                  GZ,NXM,NXMAX,NYMAX,KA,STR,MODE)
C
      DIMENSION GZ(NXM,NYMAX),KA(4,NXM,NYMAX)
      CHARACTER STR*(*),KT*80,KDL*1
C
      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.035)
      CALL SETLIN(-1,-1,7)
      KDL=STR(1:1)
      I=2
    1 IF(STR(I:I).EQ.KDL.OR.I.EQ.LEN(STR)) GOTO 2
         KT(I-1:I-1)=STR(I:I)
         I=I+1
      GOTO 1
C
    2 CALL MOVE(GX1,GY2+0.3)
      CALL TEXT(KT,I-2)
C
      CALL GMNMX2(GZ,NXM,1,NXMAX,1,1,NYMAX,1,GZMIN,GZMAX)
      IF(GZMAX.EQ.GZMIN) THEN
         GZMAX=GZMIN+1.0
      ENDIF
      CALL GQSCAL(GZMIN,GZMAX,GGMIN,GGMAX,GZSTEP)
      GZMIN=GGMIN+GZSTEP
      GZMAX=GGMAX-GZSTEP
      NSTEP=INT((GZMAX-GZMIN)/GZSTEP)+3
      IF(MOD(MODE/4,2).EQ.1) THEN
         CALL CHMODE
  701    WRITE(6,*) '## TRGR2D : ZMIN,ZMAX = ',GZMIN,GZMAX
         WRITE(6,*) '## INPUT  : ZORG,ZSTEP,NSTEP'
         READ(5,*,ERR=701,END=900) GZORG,GZSTEP,NSTEP
         CALL GRMODE
      ENDIF
C
      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSTEP)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GYSTEP)
C
      IF(GXMIN*GXMAX.LE.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GSXMIN+GXSTEP
      ENDIF
      IF(GYMIN*GYMAX.LE.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=GSYMIN+GYSTEP
      ENDIF      
C
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(-1,-1,7)
      CALL GDEFIN(GX1,GX2,GY1,GY2,GXMIN,GXMAX,GYMIN,GYMAX)
      CALL GFRAME
      CALL GSCALE(GXORG,GXSTEP,GYORG,GYSTEP,0.1,9)
      CALL GSCALE(0.0,100*GXSTEP,0.0,100*GYSTEP,0.0,0)
      CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGVLEN(2*GXSTEP))
      CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGVLEN(2*GYSTEP))
      CALL SETLIN(-1,-1,7)
      IF(GZMIN*GZMAX.GE.0.D0) THEN
         IF(MOD(MODE,2).EQ.0) THEN
            GZORG=GZMIN
         ELSE
            GZORG=GZMAX
            GZSTEP=-GZSTEP
         ENDIF
         IF(GZORG.GE.0.D0) THEN
            CALL CONTP1(GZ,NXM,NXMAX,NYMAX,GZORG,GZSTEP,NSTEP,0,0,KA)
         ELSE
            CALL CONTP1(GZ,NXM,NXMAX,NYMAX,GZORG,GZSTEP,NSTEP,0,3,KA)
         ENDIF
         CALL MOVE(GX2-10.0,GY2+0.3)
         CALL TEXT('ZORG =',6)
         CALL NUMBR(GZORG,'(1PE10.2)',10)
      ELSE
         CALL CONTP1(GZ,NXM,NXMAX,NYMAX, GZSTEP, GZSTEP,
     &               NSTEP,0,0,KA)
         CALL CONTP1(GZ,NXM,NXMAX,NYMAX, 0.0,    GZSTEP,
     &               1,0,4,KA)
         CALL CONTP1(GZ,NXM,NXMAX,NYMAX,-GZSTEP,-GZSTEP,
     &               NSTEP,0,3,KA)
      ENDIF
      CALL MOVE(GX2-5.0,GY2+0.3)
      CALL TEXT('STEP =',6)
      CALL NUMBR(GZSTEP,'(1PE10.2)',10)
C
  900 CALL SETLIN(-1,-1,7)
      RETURN
      END
C
C     ***********************************************************
C
C           SUBPROGRAM FOR 2D PROFILE
C
C     ***********************************************************
C
      SUBROUTINE TRGR2D(GX1,GX2,GY1,GY2,GZ,NXM,NXMAX,NYMAX,STR,MODE)
C
      DIMENSION GZ(NXM,NYMAX)
      CHARACTER STR*(*),KT*80,KDL*1
C
      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
      CALL SETLNW(0.035)
      CALL SETLIN(-1,-1,7)
      KDL=STR(1:1)
      I=2
    1 IF(STR(I:I).EQ.KDL.OR.I.EQ.LEN(STR)) GOTO 2
         KT(I-1:I-1)=STR(I:I)
         I=I+1
      GOTO 1
C
    2 CALL MOVE(GX1,GY2+0.3)
      CALL TEXT(KT,I-2)
C
      IXY= 3
      IND=-3
      CALL GMNMX2(GZ,NXM,1,NXMAX,1,1,NYMAX,1,GZMIN,GZMAX)
      IF(MOD(MODE,2).EQ.0) THEN
         GZMIN=0.0
      ENDIF
      IF(GZMAX.EQ.GZMIN) THEN
         GZMAX=GZMIN+1.0
      ENDIF
      IF(MOD(MODE/4,2).EQ.1) THEN
         CALL CHMODE
  701    WRITE(6,*) '## TRGR2D : ZMIN,ZMAX = ',GZMIN,GZMAX
         READ(5,*,ERR=701,END=900) GZMIN,GZMAX
         CALL GRMODE
      ENDIF
C
      XL=10.0
      YL=10.0
      ZL= 5.0
      A=-30.0
      B=0.0
      C=-30.0
      D=0.0
      E=1000.0
      X1=-10.0
      X2= 10.0
      Y1= 0.0
      Y2= 15.0
C
C   1 WRITE(6,*) 'INPUT : IND,ZMIN,ZMAX,HX,HY,HZ,A,B,C,D,E,X1,X2,Y1,Y2 '
C     READ(5,*,END=9999) IND,ZMIN,ZMAX,HX,HY,HZ,A,B,C,D,E,X1,X2,Y1,Y2
C
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(-1,-1,7)
      CALL GDEFIN(GX1,GX2,GY1,GY2,X1,X2,Y1,Y2)
      CALL GFRAME
      CALL SETLIN(-1,-1,7)
      CALL PERSE1(GZ,NXM,NXMAX,NYMAX,GZMIN,GZMAX,IXY,IND,
     &            XL,YL,ZL,A,B,C,D,E)
C
  900 CALL SETLIN(-1,-1,7)
      RETURN
      END
C
C
C     ***********************************************************
C
C           CEILING FUNCTION FOR LOG10 PLOT
C
C     ***********************************************************
C
      FUNCTION GLOG(X,XMIN,XMAX)
C
C      REAL (KIND=8) X,XMIN,XMAX,PLOG
      REAL*8 X,XMIN,XMAX,PLOG
C
      GLOG=GCLIP(PLOG(X,XMIN,XMAX))
C
      RETURN
      END
C
C     *****************************
C
C     WRITE TIME ON FIGURE
C
C     *****************************
C
      SUBROUTINE TRGRTM
C
      INCLUDE 'trcomm.h'
C
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.3,0.0)
      CALL SETFNT(32)
      CALL MOVE(11.8,18.0)
      CALL TEXT('T=',2)
      CALL NUMBD(T,'(1F7.3)',7)
      CALL TEXT(' SEC',4)
      RETURN
      END
