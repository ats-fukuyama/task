C
C     ***** READ EQ VARIABLE *****
C
      SUBROUTINE TRXREQ(KFID,T,R,Z,F2,NRM,NZM,NRXMAX,NZXMAX,MODE)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(NRM),Z(NZM),F2(NRM,NZM)
      COMMON /TRKID1/ KDIRR1,KDIRR2
      CHARACTER*80 KDIRR1,KDIRR2
      CHARACTER*80 KFID,KFILE,KFILEB
      CHARACTER*80 KMATCH1,KMATCH2,KMATCH3,KMATCH4,KMATCH5
      CHARACTER*80 KKLINE,KKLINE1,KKLINE2
      LOGICAL KMATCH
C
      IMATCH=0
      KMATCH1='# OF X PTS'
      KMATCH2='# OF Y PTS'
      KMATCH3='# OF PTS'
      KMATCH4='# OF T PTS'
      KMATCH5='# OF Z PTS'
C
      CALL KTRIM(KDIRR2,KL2)
      KFILE=KDIRR2(1:KL2)//KFID
      CALL KTRIM(KFILE,KL2)
      KFILEB=KFILE(1:KL2)//'.bin'
      CALL KTRIM(KFILEB,KL3)
C
      IF(MODE.EQ.1) THEN
         OPEN(15,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',
     &        STATUS='OLD',ERR=8)
         READ(15) T,NRXMAX,NZXMAX
         READ(15) (R(NRX),NRX=1,NRXMAX)
         READ(15) (Z(NZX),NZX=1,NZXMAX)
         READ(15) ((F2(NRX,NZX),NRX=1,NRXMAX),NZX=1,NZXMAX)
         CLOSE(15)
         WRITE(6,'(A14,A,A2,I6,A1,I6,A1)')
     &        ' - READ FILE :',KFILEB(1:KL3),' (',NRXMAX,',',NZXMAX,')'
         GOTO 9000
      ENDIF
C
    8 OPEN(15,FILE=KFILE,IOSTAT=IST,FORM='FORMATTED',
     &     STATUS='OLD',ERR=10)
      GOTO 100
C
   10 CONTINUE
      WRITE(6,*) 'XX NOT FOUND :',KFILE(1:60)
      STOP
C
  100 CONTINUE
         READ(15,'(A80)',ERR=8000,END=9000) KKLINE
         CALL KTRIM(KKLINE,KL)
         IF(KL.EQ.0) GOTO 100
         INCL=KKINDEX(KKLINE,';')
         IF(INCL.EQ.0) THEN
            IF(IMATCH.GE.2) GOTO 1000
            GOTO 100
         ENDIF
         CALL KSPLIT(KKLINE,';',KKLINE1,KKLINE2)
         CALL KTRIM(KKLINE1,KL)
         CALL KEXTR(KKLINE2)
         IF(KMATCH(KKLINE2,KMATCH1)) THEN
            READ(KKLINE1,*,ERR=8000) NRXMAX
            IMATCH=IMATCH+1
         ELSE IF(KMATCH(KKLINE2,KMATCH2)) THEN
            READ(KKLINE1,*,ERR=8000) NZXMAX
            IMATCH=IMATCH+1
         ELSE IF(KMATCH(KKLINE2,KMATCH3)) THEN
            READ(KKLINE1,*,ERR=8000) NTXMAX
            IMATCH=IMATCH+1
         ELSE IF(KMATCH(KKLINE2,KMATCH4)) THEN
            READ(KKLINE1,*,ERR=8000) NTXMAX
            IMATCH=IMATCH+1
         ELSE IF(KMATCH(KKLINE2,KMATCH5)) THEN
            READ(KKLINE1,*,ERR=8000) NFXMAX
            IMATCH=IMATCH+1
         ENDIF
      GOTO 100
C
 1000 CONTINUE
C
      BACKSPACE 15
C
      CALL NDINIT
      DO NTX=1,NTXMAX
         CALL NDREAD(15,T,IERR)
         IF(IERR.NE.0) GOTO 8000
      ENDDO
      DO 1030 NRX=1,NRXMAX
         CALL NDREAD(15,R(NRX),IERR)
         IF(IERR.NE.0) GOTO 8000
 1030 CONTINUE
      DO 1040 NZX=1,NZXMAX
         CALL NDREAD(15,Z(NZX),IERR)
         IF(IERR.NE.0) GOTO 8000
 1040 CONTINUE
      DO 1050 NZX=1,NZXMAX
      DO 1050 NRX=1,NRXMAX
         CALL NDREAD(15,F2(NRX,NZX),IERR)
         IF(IERR.NE.0) GOTO 8000
 1050 CONTINUE
C
      CLOSE(15)
      WRITE(6,'(A14,A,A2,I6,A1,I6,A1)') ' - READ FILE :',KFILE(1:KL2),
     &           ' (',NRXMAX,',',NZXMAX,')'
C
      IF(MODE.EQ.1) THEN
         OPEN(15,FILE=KFILEB,IOSTAT=IST,FORM='UNFORMATTED',
     &        STATUS='NEW',ERR=8010)
         WRITE(15) T,NRXMAX,NZXMAX
         WRITE(15) (R(NRX),NRX=1,NRXMAX)
         WRITE(15) (Z(NZX),NZX=1,NZXMAX)
         WRITE(15) ((F2(NRX,NZX),NRX=1,NRXMAX),NZX=1,NZXMAX)
         WRITE(6,*) ' - WRITE FILE:',KFILEB(1:60)
         CLOSE(15)
      ENDIF
C
 9000 IF(KFID(1:6).EQ.'VOLUME'.OR.
     &   KFID(1:4).EQ.'SURF'.OR.
     &   KFID(1:6).EQ.'RMINOR') THEN
      IF(R(1).NE.0.D0) THEN
         NRXMAX=NRXMAX+1
         DO NRX=NRXMAX,2,-1
            R(NRX)=R(NRX-1)
         ENDDO
         R(1)=0.D0
         DO NZX=1,NZXMAX
            DO NRX=NRXMAX,2,-1
               F2(NRX,NZX)=F2(NRX-1,NZX)
            ENDDO
            F2(1,NZX)=0.D0
         ENDDO
      ENDIF
      ENDIF
      RETURN
C
 8000 WRITE(6,*) 'XX FILE READ ERROR'
      RETURN
 8010 WRITE(6,*) 'XX FILE OPEN ERROR: STATSU = ',IST
      RETURN
      END
C
C     ****** INTERPOLATE 1D VARIABLE ******
C
      SUBROUTINE TRXT1D(T0,F0,T,F1,NTM,NTXMAX)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION T(NTM),F1(NTM)
C
      IF(NTXMAX.EQ.0) THEN
         RETURN
      ENDIF
C
      IF(NTXMAX.EQ.1) THEN
         F0=F1(1)
         RETURN
      ENDIF
C
      TMIN=T(1)
      TMAX=T(NTXMAX)
      NT=INT((NTXMAX-1)*(T0-TMIN)/(TMAX-TMIN))+1
    1 IF(NT.LT.1) THEN
         F0=F1(1)
         GOTO 2000
      ELSE IF(NT.GE.NTXMAX) THEN
         F0=F1(NTXMAX)
         GOTO 2000
      ELSE
         IF(T0.LT.T(NT)) THEN
            NT=NT-1
         ELSE 
            IF(T0.GT.T(NT+1)) THEN
               NT=NT+1
            ELSE
               GOTO 1000
            ENDIF
         ENDIF
      ENDIF
      GOTO 1
C
 1000 CONTINUE
      F0=F1(NT)+(F1(NT+1)-F1(NT))*(T0-T(NT))/(T(NT+1)-T(NT))
 2000 RETURN
      END
C
C     ****** INTERPOLATE EXP ******
C
      SUBROUTINE TRXT2D(T0,R0,F0,T,R,F2,NRM,NTM,NRXMAX,NTXMAX)
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION T(NTM),R(NRM),F2(NRM,NTM)
C
      IF(NTXMAX.EQ.0) THEN
         WRITE(6,*) 'XX WRONG CALL OF TRXT2D'
         RETURN
      ENDIF
C
      IF(NTXMAX.EQ.1) THEN
         RMIN=R(1)
         RMAX=R(NRXMAX)
         NR=INT(NRXMAX*(R0-RMIN)/(RMAX-RMIN))+1
 3000    IF(NR.LT.1) THEN
            NR=0
            GOTO 3010
         ELSE IF(NR.GE.NRXMAX) THEN
            NR=NRXMAX
            GOTO 3010
         ELSE
            IF(R0.LT.R(NR)) THEN
               NR=NR-1
            ELSE 
               IF(R0.GT.R(NR+1)) THEN
                  NR=NR+1
               ELSE
                  GOTO 3010
               ENDIF
            ENDIF
         ENDIF
         GOTO 3000
C
 3010    CONTINUE
         IF(NR.EQ.0) THEN
            NR1=1
            NR2=1
            FACTORR=0.D0
         ELSE IF(NR.EQ.NRXMAX) THEN
            NR1=NRXMAX-1
            NR2=NRXMAX
            FACTORR=(R0-R(NR1))/(R(NR2)-R(NR1))
         ELSE
            NR1=NR
            NR2=NR+1
            FACTORR=(R0-R(NR1))/(R(NR2)-R(NR1))
         ENDIF
         FR2=F2(NR2,1)
         FR1=F2(NR1,1)
         F0=FR1+FACTORR*(FR2-FR1)
         RETURN
      ENDIF
C
      TMIN=T(1)
      TMAX=T(NTXMAX)
      NT=INT(NTXMAX*(T0-TMIN)/(TMAX-TMIN))+1
    1 IF(NT.LT.1) THEN
         NT=0
         GOTO 1000
      ELSE IF(NT.GE.NTXMAX) THEN
         NT=NTXMAX
         GOTO 1000
      ELSE
         IF(T0.LT.T(NT)) THEN
            NT=NT-1
         ELSE 
            IF(T0.GT.T(NT+1)) THEN
               NT=NT+1
            ELSE
               GOTO 1000
            ENDIF
         ENDIF
      ENDIF
      GOTO 1
C
 1000 RMIN=R(1)
      RMAX=R(NRXMAX)
      NR=INT(NRXMAX*(R0-RMIN)/(RMAX-RMIN))+1
    2 IF(NR.LT.1) THEN
         NR=0
         GOTO 2000
      ELSE IF(NR.GE.NRXMAX) THEN
         NR=NRXMAX
         GOTO 2000
      ELSE
         IF(R0.LT.R(NR)) THEN
            NR=NR-1
         ELSE 
            IF(R0.GT.R(NR+1)) THEN
               NR=NR+1
            ELSE
               GOTO 2000
            ENDIF
         ENDIF
      ENDIF
      GOTO 2
C
 2000 CONTINUE
      IF(NT.EQ.0) THEN
         NT1=1
         NT2=1
         FACTORT=0.D0
      ELSE IF(NT.EQ.NTXMAX) THEN
         NT1=NTXMAX
         NT2=NTXMAX
         FACTORT=0.D0
      ELSE
         NT1=NT
         NT2=NT+1
         FACTORT=(T0-T(NT1))/(T(NT2)-T(NT1))
      ENDIF
      IF(NR.EQ.0) THEN
         NR1=1
         NR2=2
      ELSE IF(NR.EQ.NRXMAX) THEN
         NR1=NRXMAX-1
         NR2=NRXMAX
      ELSE
         NR1=NR
         NR2=NR+1
      ENDIF
      FACTORR=(R0-R(NR1))/(R(NR2)-R(NR1))
      FR2=F2(NR2,NT1)+FACTORT*(F2(NR2,NT2)-F2(NR2,NT1))
      FR1=F2(NR1,NT1)+FACTORT*(F2(NR1,NT2)-F2(NR1,NT1))
      F0=FR1+FACTORR*(FR2-FR1)
C
C      WRITE(6,*) NT1,NT2,NR1,NR2
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
      SUBROUTINE TRGR1D(GX1,GX2,GY1,GY2,GX,GY,NXM,NXMAX,NGMAX,STR,MODE)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GX(NXM),GY(NXM,NGMAX),IPAT(5)
      CHARACTER STR*80,KT*80,KDL*1
      DATA IPAT/0,2,3,4,6/
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
         WRITE(6,*) '## TRGR : XMIN,XMAX,YMIN,YMAX = ',
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
      CALL SETLIN(0,0,4)
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
         CALL SETLIN(0,0,7-MOD(NG-1,5))
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
      CALL SETLIN(0,0,7)
      RETURN
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
      CHARACTER STR*80,KT*80,KDL*1
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
      IF(GZMAX.EQ.GZMIN) THEN
         GZMAX=GZMIN+1.0
      ENDIF
      CALL GQSCAL(GZMIN,GZMAX,GGMIN,GGMAX,GZSTEP)
      GZMIN=GGMIN+GZSTEP
      GZMAX=GGMAX-GZSTEP
      NSTEP=(GZMAX-GZMIN)/GZSTEP+3
      IF(MOD(MODE/4,2).EQ.1) THEN
         CALL CHMODE
         WRITE(6,*) '## TRGR2D : ZMIN,ZMAX = ',GZMIN,GZMAX
         WRITE(6,*) '## INPUT  : ZORG,ZSTEP,NSTEP'
         READ(5,*) GZORG,GZSTEP,NSTEP
         CALL GRMODE
      ENDIF
C
      CALL GQSCAL(GXMIN,GXMAX,GGMIN,GGMAX,GXSTEP)
      CALL GQSCAL(GYMIN,GYMAX,GGMIN,GGMAX,GYSTEP)
C
      IF(GXMIN*GXMAX.LE.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GXMIN
      ENDIF
      IF(GYMIN*GYMAX.LE.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=GYMIN
      ENDIF      
C
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(0,0,4)
      CALL GDEFIN(GX1,GX2,GY1,GY2,GXMIN,GXMAX,GYMIN,GYMAX)
      CALL GFRAME
      CALL GSCALE(GXORG,GXSTEP,GYORG,GYSTEP,0.1,9)
      CALL GSCALE(0.0,100*GXSTEP,0.0,100*GYSTEP,0.0,0)
      CALL GVALUE(GXORG,2*GXSTEP,0.0,0.0,NGVLEN(2*GXSTEP))
      CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGVLEN(2*GYSTEP))
      CALL SETLIN(0,0,7)
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
         CALL MOVE(GX2-10.0,GY2+0.2)
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
      CALL MOVE(GX2-5.0,GY2+0.2)
      CALL TEXT('STEP =',6)
      CALL NUMBR(GZSTEP,'(1PE10.2)',10)
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
      SUBROUTINE TRGR2D(GX1,GX2,GY1,GY2,GZ,NXM,NXMAX,NYMAX,STR,MODE)
C
      DIMENSION GZ(NXM,NYMAX)
      CHARACTER STR*80,KT*80,KDL*1
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
         WRITE(6,*) '## TRGR2D : ZMIN,ZMAX = ',GZMIN,GZMAX
         READ(5,*) GZMIN,GZMAX
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
      CALL SETLIN(0,0,4)
      CALL GDEFIN(GX1,GX2,GY1,GY2,X1,X2,Y1,Y2)
      CALL GFRAME
      CALL SETLIN(0,0,7)
      CALL PERSE1(GZ,NXM,NXMAX,NYMAX,GZMIN,GZMAX,IXY,IND,
     &            XL,YL,ZL,A,B,C,D,E)
C
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
      NGX = -LOG10(DBLE(GSTEP*0.11))
      IF(NGX.LT.-5)THEN
         NGX=-1
      ELSE
         IF(NGX.LT.0) NGX=0
         IF(NGX.GT.5) NGX=-1
      ENDIF
      NGVLEN=NGX
      RETURN
      END
