C     $Id$
C
C     ****** DRAW COUNTOUR IN MAGNETIC SURFACE COORDINATES ******
C
      SUBROUTINE EQGR2D(GF,GR,GZ,GRS,GZS,NRM,NRMAX,NTHMAX,KA,TXT)
C
      DIMENSION GF(NRM,NTHMAX),GR(NRM,NTHMAX),GZ(NRM,NTHMAX)
      DIMENSION GRS(NTHMAX+1),GZS(NTHMAX+1)
      DIMENSION KA(4,NRM,NTHMAX)
      CHARACTER TXT*72
C
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(0,2,7)
C
      CALL GMNMX2(GF,NRM,1,NRMAX,1,1,NTHMAX,1,GFMIN,GFMAX)
      CALL GQSCAL(GFMIN,GFMAX,GGFMIN,GGFMAX,GGFSTP)
      GGFSTP=0.5*GGFSTP
      NSTEP=INT((GGFMAX-GGFMIN)/GGFSTP)+1
C
      CALL GMNMX1(GRS,1,NTHMAX+1,1,GRMIN,GRMAX)
      CALL GMNMX1(GZS,1,NTHMAX+1,1,GZMIN,GZMAX)
C
      GRLEN=GRMAX-GRMIN
      GZLEN=GZMAX-GZMIN
      IF(GRLEN.GT.GZLEN) THEN
         GPR=15.0
         GPZ=15.0*GZLEN/GRLEN
      ELSE
         GPR=15.0*GRLEN/GZLEN
         GPZ=15.0
      ENDIF
      CALL GQSCAL(GRMIN,GRMAX,GGRMIN,GGRMAX,GGRSTP)
      CALL GQSCAL(GZMIN,GZMAX,GGZMIN,GGZMAX,GGZSTP)
C
      CALL GDEFIN(2.0,2.0+GPR,2.0,2.0+GPZ,GRMIN,GRMAX,GZMIN,GZMAX)
      CALL SETLIN(0,2,7)
      CALL GFRAME
      CALL GSCALE(GGRMIN+GGRSTP,GGRSTP,0.0,0.0,0.1,9)
      CALL GVALUE(GGRMIN+GGRSTP,GGRSTP*2,0.0,0.0,1)
      CALL GSCALE(0.0,0.0,0.0,GGZSTP,0.1,9)
      CALL GVALUE(0.0,0.0,0.0,GGZSTP*2,1)
C
      CALL SETLIN(0,-1,7)
      IF(GFMIN*GFMAX.GT.0.) THEN
         CALL CONTP5(GF,GR,GZ,NRM,NRMAX,NTHMAX,
     &               GGFMIN,GGFSTP,NSTEP,2,0,KA)
      ELSE
         CALL CONTP5(GF,GR,GZ,NRM,NRMAX,NTHMAX,
     &                0.5*GGFSTP, GGFSTP,NSTEP,2,0,KA)
         CALL CONTP5(GF,GR,GZ,NRM,NRMAX,NTHMAX,
     &               -0.5*GGFSTP,-GGFSTP,NSTEP,2,2,KA)
      ENDIF
C
      CALL SETLIN(0,-1,5)
      CALL GPLOTP(GRS,GZS,1,NTHMAX+1,1,0,0,0)
C
      CALL SETLIN(0,-1,7)
      CALL MOVE(20.0,17.0)
      CALL TEXT('MAX :',5)
      CALL NUMBR(GFMAX,'(1PE12.4)',12)
      CALL MOVE(20.0,16.5)
      CALL TEXT('MIN :',5)
      CALL NUMBR(GFMIN,'(1PE12.4)',12)
      CALL MOVE(20.0,16.0)
      CALL TEXT('STEP:',5)
      CALL NUMBR(GGFSTP,'(1PE12.4)',12)
      CALL MOVE(2.0,17.3)
      CALL TEXTX(TXT)
C
      RETURN
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
      SUBROUTINE EQGR1D(GX1,GX2,GY1,GY2,GX,GY,NXM,NXMAX,NGMAX,STR,MODE)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GX(NXM),GY(NXM,NGMAX),IPAT(5)
      CHARACTER STR*80,KT*80,KDL*1
      DATA IPAT/0,2,3,4,6/
C
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(0,2,7)
      KDL=STR(1:1)
      I=2
    1 IF(STR(I:I).EQ.KDL.OR.I.EQ.80) GOTO 2
         KT(I-1:I-1)=STR(I:I)
         I=I+1
      GOTO 1
C
    2 CALL MOVE(GX1,GY2+0.1)
      CALL TEXT(KT,I-2)
C
      CALL GMNMX2(GY,NXM,1,NXMAX,1,1,NGMAX,1,GYMIN,GYMAX)
      IF(ABS(GYMAX-GYMIN).LT.1.E-6) THEN
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
      IF(ABS(GXMAX-GXMIN).LT.1.E-6) THEN
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
         WRITE(6,*) '## WMGR : XMIN,XMAX,YMIN,YMAX = ',
     &              GXMIN,GXMAX,GYMIN,GYMAX
         READ(5,*) GXMIN,GXMAX,GYMIN,GYMAX
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
      CALL SETLIN(0,2,7)
      CALL GFRAME
      CALL GSCALE(GXORG,GSTEPX,0.0,0.0,0.1,9)
      CALL GSCALE(0.0,0.0,0.0,2*(GSYMAX-GSYMIN),0.0,9)
      CALL GVALUE(GXORG,2*GSTEPX,0.0,0.0,NGULEN(GSTEPX))
      IF(MOD(MODE/8,2).EQ.1) THEN
         CALL GSCALL(0.0,0,GYORG,9,0.1,9)
         CALL GVALUL(0.0,0,GYORG,1,NGULEN(GSTEPY))
      ELSE
         CALL GSCALE(0.0,0.0,GYORG,GSTEPY,0.1,9)
         CALL GVALUE(0.0,0.0,GYORG,2*GSTEPY,NGULEN(GSTEPY))
      ENDIF
C
      DO 10 NG=1,NGMAX
         CALL SETLIN(0,-1,7-MOD(NG-1,5))
         IF(MOD(MODE/2,2).EQ.0) THEN
            CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,0)
         ELSE
            CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,IPAT(MOD(NG-1,5)+1))
         ENDIF
   10 CONTINUE
C
      CALL SETLIN(0,-1,4)
      RETURN
      END

































