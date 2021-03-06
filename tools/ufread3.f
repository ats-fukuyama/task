      PROGRAM UFILE_READER_WITH_GSGL
      USE libitp
      USE libchar
      IMPLICIT REAL*8(A-F,H,O-Z)
      PARAMETER (NRM=100,NTM=2001)
      DIMENSION T(NTM),R(NRM),F1(NTM),F2(NRM,NTM),RL(NRM),FL(NRM)
      DIMENSION GT(NTM),GR(NRM),GF1(NTM),GF2(NRM,NTM)
      COMMON /TRKID1/ KDIRR1,KDIRR2
      CHARACTER KDEV*80,KDCG*80
      CHARACTER KDIRR1*80,KDIRR2*80,KDIR*80,KVAL*80
      CHARACTER KDIRX*80
      CHARACTER KFID*80,KFIDX*80,KVAR*80
      CHARACTER KFIDCK*90,KFILE*110
      CHARACTER KTIME*12
      REAL STIME
      LOGICAL LEX
      COMMON /DIR/ KDIR
C
      CALL GSOPEN
C
      KDEV='jt60u'
      KDCG='29728'
      KDIR='data'
      CALL parameter_read
      NDIM=2
C
    1 WRITE(6,*) '# DEVICE NAME ?'
      READ(5,'(A40)',ERR=1,END=9000) KDEV
    2 WRITE(6,*) '# DISCHARGE NUMBER ?'
      READ(5,'(A40)',ERR=2,END=1) KDCG
C
      CALL KTRIM(KDEV,IKDEV)
      CALL KTRIM(KDCG,IKDCG)
      CALL KTRIM(KDIR,IKDIR )
      CALL TOLOWER(KDEV)
C
      KDIRX=KDIR(1:IKDIR)//'/'//KDEV(1:IKDEV)//'/'
     &    //KDCG(1:IKDCG)//'/in/'
      CALL KTRIM(KDIRX,IKDIRX)
      KFILE=KDIRX(1:IKDIRX)//KDEV(1:IKDEV)//'2d'//KDCG(1:IKDCG)//
     &     '.NE'
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
      IF(LEX.EQV..FALSE.) THEN
         WRITE(6,600) 'XX: DIRECTORY DOES NOT EXIST! ',KDIRX
         GOTO 1
      ENDIF
 600  FORMAT(' ',A31,A)
C
      KDIRR1=KDIRX(1:IKDIRX)//KDEV(1:IKDEV)
     &       //'1d'//KDCG(1:IKDCG)//'.'
      KDIRR2=KDIRX(1:IKDIRX)//KDEV(1:IKDEV)
     &       //'2d'//KDCG(1:IKDCG)//'.'
C
    5 WRITE(6,*) 'INPUT DIM'
      READ(5,'(I1)',ERR=5,END=9000) NDIM
      IF(NDIM.LE.0) GOTO 9000
    6 WRITE(6,*) 'INPUT FILEID'
      READ(5,'(A80)',ERR=6,END=5) KFID
      CALL KTRIM(KDIRR1,IKDIRR1)
      CALL KTRIM(KDIRR2,IKDIRR2)
      CALL KTRIM(KFID,IKFID)
C     converting lower case to upper case
      DO I=1,IKFID
         CALL toupper(KFID(I:I))
      ENDDO
C
      IF(NDIM.EQ.1) THEN
         KFIDCK=KDIRR1(1:IKDIRR1)//KFID(1:IKFID)
      ELSEIF(NDIM.EQ.2) THEN
         KFIDCK=KDIRR2(1:IKDIRR2)//KFID(1:IKFID)
      ELSE
         GOTO 100
      ENDIF
      INQUIRE(FILE=KFIDCK,EXIST=LEX,ERR=9000)
      IF(LEX.EQV..FALSE.) THEN
         WRITE(6,610) 'XX: FILE DOES NOT EXIST! ',KFIDCK
         GOTO 5
      ENDIF
 610  FORMAT (' ',A25,A)
C
 100  IF(NDIM.EQ.1) THEN
         CALL TRXR1D(KDIRR1,KFID,T,F1,NTM,NTXMAX,1)
      ELSE IF(NDIM.EQ.2) THEN
         CALL TRXR2D(KDIRR2,KFID,T,R,F2,NRM,NTM,NRXMAX,NTXMAX,1)
      ELSE IF(NDIM.EQ.3) THEN
         CALL TRXR1D(KDIRR1,KFID,T,F1,NTM,NTXMAX,1)
    7    WRITE(6,*) '# INPUT T'
         READ(5,*,ERR=7,END=5) TL
         IF(TL.LT.0.0) GOTO 5
         CALL TRXT1D(TL,FLL,T,F1,NTM,NTXMAX)
         WRITE(6,'(1PE12.4)') FLL
         GOTO 7
      ELSE IF(NDIM.EQ.4) THEN
         CALL TRXR2D(KDIRR2,KFID,T,R,F2,NRM,NTM,NRXMAX,NTXMAX,1)
    8    WRITE(6,*) '# INPUT T,R'
         READ(5,*,ERR=8,END=5) TL,RLL
         IF(TL.LT.0.0) GOTO 5
         CALL TRXT2D(TL,RLL,FLL,T,R,F2,NRM,NTM,NRXMAX,NTXMAX)
         WRITE(6,'(1PE12.4)') FLL
         GOTO 8
      ENDIF
C      write(6,*) NDIM,NRM,NTM,NRXMAX,NTXMAX
C     
      GX1=2.0
      GX2=17.0
      GY1=2.0
      GY2=17.0
C
C      WRITE(6,*) NRXMAX,NTXMAX
C      WRITE(6,'(1P6E12.4)') (R(NRX),NRX=1,NRXMAX)
C      WRITE(6,'(1P6E12.4)') (T(NTX),NTX=1,NTXMAX)
C      WRITE(6,'(1P6E12.4)') ((F2(NRX,NTX),NRX=1,NRXMAX),NTX=1,NTXMAX)
C
      CALL PAGES
      IF(NDIM.EQ.1) THEN
         IF(NTXMAX.EQ.1) THEN
            WRITE(KVAL,'(A13)') '@** VALUE **@'
            CALL GTEXTX(18.0,15.0,KVAL,0)
            WRITE(KVAL,'(A9,1PE9.3,A1)') '@TIME  = ',REAL(T(1)),'@'
            CALL GTEXTX(18.5,14.3,KVAL,0)
            WRITE(KVAL,'(A9,1PE9.3,A1)') '@VALUE = ',REAL(F1(1)),'@'
            CALL GTEXTX(18.5,13.7,KVAL,0)
         ELSE
            DO NTX=1,NTXMAX
               GT(NTX)=GUCLIP(T(NTX))
               GF1(NTX)=GUCLIP(F1(NTX))
            ENDDO
            CALL KTRIM(KFID,KL)
            KFIDX='@'//KDEV(1:IKDEV)//'/'//KDCG(1:IKDCG)//'/'
     &           //KFID(1:KL)//'@'
            CALL TRGR1D(GX1,GX2,GY1,GY2,
     &                  GT,GF1,NTM,NTXMAX,1,KFIDX,2)
         ENDIF
      ELSE
C
         CALL label_of_slice_time(STIME,T,NTM,NTXMAX,NTSL,IERR)
         IF(IERR.NE.0) STOP
C
         DO NTX=1,NTXMAX
            GT(NTX)=GUCLIP(T(NTX))
         ENDDO
C
         NTXMAX1=NTXMAX
         IF(KFID(1:1).EQ.'T') THEN
            FACT=1.D-3
         ELSEIF(KFID(1:1).EQ.'N') THEN
            FACT=1.D-19
         ELSE
            FACT=1.D0
         ENDIF
         IF(NTSL.EQ.0) THEN
C     *** F-t-r graph ***
            DO NRX=1,NRXMAX
               GR(NRX)=GUCLIP(R(NRX))
               DO NTX=1,NTXMAX
                  GF2(NRX,NTX)=GUCLIP(F2(NRX,NTX)*FACT)
               ENDDO
            ENDDO
         ELSE
C     *** F-r graph ***
            NRLMAX=NRXMAX
            DO NRX=1,NRXMAX
               NTXMAX=1
               RL(NRX)=R(NRX)
               FL(NRX)=F2(NRX,NTSL)
               GR(NRX)=GUCLIP(R(NRX))
               GF2(NRX,1)=GUCLIP(F2(NRX,NTSL)*FACT)
            ENDDO
            IF(RL(1).NE.0.D0) THEN
               DO NRL=NRLMAX,1,-1
                  RL(NRL+1)=RL(NRL)
                  FL(NRL+1)=FL(NRL)*FACT
               ENDDO
               RL(1)=0.D0
               FL(1)=FCTR(RL(2),RL(3),FL(2),FL(3))
               NRLMAX=NRLMAX+1
            ENDIF
            IF(RL(NRLMAX).NE.1.D0) THEN
               RL(NRLMAX+1)=1.D0
               FL(NRLMAX+1)=AITKEN2P(1.D0,FL(NRLMAX),FL(NRLMAX-1),
     &                               FL(NRLMAX-2),RL(NRLMAX),
     &                               RL(NRLMAX-1),RL(NRLMAX-2))*FACT
               NRLMAX=NRLMAX+1
            ENDIF
         ENDIF
C
c$$$         CALL KTRIM(KFID,KL)
c$$$         KFIDX='@'//KDEV(1:IKDEV)//'/'//KDCG(1:IKDCG)//'/'
c$$$     &            //KFID(1:KL)//'@'
c$$$         IF (NTXMAX.EQ.1) THEN
c$$$            CALL TRGR1D(GX1,GX2,GY1,GY2,
c$$$     &           GR,GF2,NRM,NRXMAX,NTXMAX,KFIDX,2)
c$$$         ELSE
c$$$            CALL GSGLENABLELIGHTING
c$$$            KVAR='@'//KFID(1:KL)//'@'
c$$$            CALL GRAPH3(GX1,GX2,GY1,GY2,GT,
c$$$     &           GR,GF2,NTM,NRM,NRXMAX,NTXMAX,KFIDX,KVAR,2)
c$$$         ENDIF
         CALL KTRIM(KDEV,IKDEV)
         CALL KTRIM(KDCG,IKDCG)
         CALL KTRIM(KFID,IKFID)
C         KFIDX='@'//KDEV(1:IKDEV)//'/'//KDCG(1:IKDCG)//'/'
C     &            //KFID(1:IKFID)//'@'
         KFIDX=' '
         IF (NTXMAX.EQ.1) THEN
            CALL TRGR1D(GX1,GX2,GY1,GY2,
     &                  GR,GF2,NRM,NRXMAX,NTXMAX,KFIDX,2)
            WRITE(KVAL,'(A26)') '@** INTERPOLATED VALUE **@'
            CALL GTEXTX(18.0,15.0,KVAL,0)
            WRITE(KVAL,'(A16)') '@RHO      VALUE@'
            CALL GTEXTX(18.5,14.3,KVAL,0)
            DO NL=1,11
               RP=(NL-1)*0.1D0
               CALL AITKEN(RP,FP,RL,FL,2,NRLMAX)
               WRITE(KVAL,'(A1,F6.4,A3,1PE9.3,A1)') '@',RP,'   ',FP,'@'
               CALL GTEXTX(18.0,14.0-0.6*REAL(NL),KVAL,0)
            ENDDO
            IF(NTXMAX1.NE.1.AND.NTSL.NE.1) THEN
               WRITE(KTIME,'(A5,F6.3,A1)') '@t = ',STIME,'@'
               CALL GTEXTX(21.2, 2.0, KTIME, 0)
            ENDIF
         ELSE IF(NTXMAX.GT.1) THEN
            CALL GSGLENABLELIGHTING
            KVAR='@'//KFID(1:IKFID)//'@'
            CALL GRAPH3(GX1,GX2,GY1,GY2,GT,
     &           GR,GF2,NTM,NRM,NRXMAX,NTXMAX,KFIDX,KVAR,2)
         ENDIF
      ENDIF
      CALL PAGEE
      GOTO 5
C
 9000 CALL GSCLOS
      STOP
      END
C
C     *** RETRIEVING LABEL NUMBER OF SLICE TIME ***
C
      SUBROUTINE label_of_slice_time(STIME,TL,NTM,NTXMAX,NTSL,IERR)
      IMPLICIT NONE
      INTEGER NTX, NTM, NTXMAX, NTSL, IERR
      REAL STIME
      REAL*8 TL_MIN, TL_MIN_OLD
      REAL*8 TL(NTM)
      CHARACTER KERR*100
C
      IF(NTXMAX.EQ.1) THEN
         NTSL=NTXMAX
         RETURN
      ENDIF
      IF(STIME.EQ.0.0) THEN
         NTSL=0
         RETURN
      ENDIF
      IF(DBLE(STIME).LT.TL(1).OR.DBLE(STIME).GT.TL(NTXMAX)) GOTO 100
C
      DO NTX=1,NTXMAX
         IF(ABS(TL(NTX)-DBLE(STIME)).LE.1.D-5) THEN
            NTSL=NTX
            RETURN
         ENDIF
      ENDDO
C
      TL_MIN=TL(NTXMAX)
      DO NTX=1,NTXMAX
         TL_MIN_OLD=TL_MIN
         TL_MIN=MIN(ABS(TL(NTX)-DBLE(STIME)),TL_MIN)
         IF(TL_MIN_OLD.EQ.TL_MIN) THEN
            NTSL=NTX-1
            STIME=SNGL(TL(NTSL))
            RETURN
         ENDIF
      ENDDO
C
 100  WRITE(KERR,'(A,2(F7.4,A))') 
     &     '@XX: SLICE TIME IS OUT OF THE RANGE OF ',
     &     SNGL(TL(1)),' - ',SNGL(TL(NTXMAX)),'@'
      CALL GTEXTX(10.0,10.0,KERR,2)
      IERR=1
C
      RETURN
      END
C
C     *** READING PARAMETER FILE ***
C
      SUBROUTINE parameter_read
C
      USE libchar
      IMPLICIT NONE
      INTEGER IDOPEN, IST1, IST2, KL
      LOGICAL LEX
      CHARACTER*80 LINE, DIR
      COMMON /DIR/ DIR
      NAMELIST /UFPARM/ DIR
C
      IDOPEN=25
      LINE='ufparm'
C
      INQUIRE(FILE=LINE,EXIST=LEX,ERR=9100)
      IF(.NOT.LEX) THEN
         WRITE(6,'(A)') '## NO INPUT FILE EXISTS.'
         RETURN
      ENDIF
      OPEN(IDOPEN,FILE=LINE,IOSTAT=IST1,STATUS='OLD',ERR=9100)
      READ(IDOPEN,UFPARM,IOSTAT=IST2,ERR=9200,END=9900)
      CLOSE(IDOPEN)
C
      CALL KTRIM(LINE,KL)
      WRITE(6,'(A,A,A)') 
     &     '## FILE (',LINE(1:KL),') IS ASSIGNED FOR PARM INPUT'
      GOTO 9900
C
 9100 WRITE(6,'(A,I6)') 'XX: FAILED TO OPEN PARM FILE : IOSTAT = ', IST1
      STOP
 9200 WRITE(6,'(A,I6)') 'XX: FAILED TO READ PARM FILE : IOSTAT = ', IST2
      STOP
 9900 RETURN
      END
C
C     *** 3D GRAPHICS SUBROUTINE ***
C
      SUBROUTINE GRAPH3(GX1,GX2,GY1,GY2,GT,GX,GY,
     &     NTM,NXM,NXMAX,NGMAX,STR,KV,MODE)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      DIMENSION GX(NXM),GY(NXM,NGMAX),GT(NTM)
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
C      DO NG=1, NGMAX
C         GT(NG)=REAL(NG)
C      ENDDO
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
      CALL GQSCAL(GT(1),GT(NGMAX),GSTMIN,GSTMAX,GSTEPT)
C      write(6,*) gt(1),gt(ngmax),gstept
C
      GXL=10.0*1.5
      GYL=20.0*1.5
      GZL=10.0*1.5
      CALL GDEFIN3D(GXL,GYL,GZL,GSXMIN,GSXMAX,GT(1),GT(NGMAX),
     &     GSYMIN,GSYMAX)
      GPHI=-20.0
      GTHETA=60.0
      GRADIUS=15.0
      GOX=0.5*(GSXMIN+GSXMAX)
      GOY=0.5*(GT(1)+GT(NGMAX))
      GOZ=0.5*(GSYMIN+GSYMAX)
      CALL GVIEW3D(GPHI,GTHETA,GRADIUS,GOX,GOY,GOZ)
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(0,0,4)
C
      CALL GSCALE3DX(GSXMIN,GSTEPX,0.3,0)
      CALL GSCALE3DY(GT(1),GSTEPT,0.3,0)
      CALL GSCALE3DZ(GSYMIN,GSTEPY,0.3,10)
      CALL GVALUE3DX(GSXMIN,GSTEPX,1,1)
      CALL GVALUE3DY(GT(1),GSTEPT,1,1)
      CALL GVALUE3DZ(GSYMIN,GSTEPY,11,-2)
C
      CALL Set3DTextBaseLine(0.0, 1.0, 0.0, -1.0, 0.0, 0.0)
      CALL GTEXTX3D(GSXMAX+0.15*(GSXMAX-GSXMIN),
     &              0.5*(GT(1)+GT(NGMAX)),
     &              GSYMIN,
     &              '@TIME (sec)@',
     &              2)
      CALL Set3DTextBaseLine(0.0, 1.0, 0.0, -1.0, 0.0, 0.0)
      CALL GTEXTX3D(0.5*(GSXMIN+GSXMAX),
     &              GT(1)+0.1*(GT(1)-GT(NGMAX)),
     &              GSYMIN,
     &              '@RHO@',
     &              2)
      CALL Set3DTextBaseLine(0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
      CALL GTEXTX3D(GSXMIN,
     &              GT(1)+0.05*(GT(1)-GT(NGMAX)),
     &              GSYMAX+0.1*(GSYMAX-GSYMIN),
     &              KV,
     &              2)
C
      CALL PERS3D1(GY,NXM,NXMAX,NGMAX,-27,R2G2B)
      CALL GAxis3D(0)
      CALL GDrawBack3D(0.5, 0.5, 0.5)
C      DO NG=1,NGMAX
C         DO NX=1,NXMAX-1
C            CALL MOVE3D(GX(NX),GT(NG),GY(NX,NG))
C            CALL DRAW3D(GX(NX+1),GT(NG),GY(NX+1,NG))
C         ENDDO
C      ENDDO
C
      CALL SETLIN(0,0,4)
      RETURN
      END
C
C     ***************************************************************
C
C        Convert Strings to Lower Case
C
C     ***************************************************************
C
      SUBROUTINE TOLOWER(KTEXT)
C
      IMPLICIT NONE
      CHARACTER KTEXT*(*)
C
      INTEGER NCHAR, I, ID
      INTEGER IASC(256)
C
      NCHAR = LEN(KTEXT)
      CALL CHRASC(KTEXT, IASC, NCHAR)
      DO I = 1, NCHAR
         ID = IASC(I)
         IF(ID .GE. 65 .AND. ID .LE. 90) IASC(I) = ID + 32
      END DO
      CALL ASCCHR(IASC, KTEXT, NCHAR)
C
      RETURN
      END
