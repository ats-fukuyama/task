C     ***********************************************
C     *                                             *
C     *    GUI-based UFILE reader                   *
C     *         Programed by HONDA Mitsuru          *
C     *                                             *
C     ***********************************************
C
      PROGRAM GUIREAD
      IMPLICIT NONE
C
      CALL GSOPEN
C
      CALL PAGES
      CALL MAIN
      CALL PAGEE
C
      CALL GSCLOS
C
      STOP
      END
C
C     *** READING PARAMETER FILE ***
C
      SUBROUTINE parameter_read
C
      IMPLICIT NONE
      INTEGER IDOPEN, IST1, IST2, DIM, KL
      REAL STIME
      LOGICAL LEX
      CHARACTER*80 LINE, DIR, DEV, SHOT, ID
      NAMELIST /GUIPARM/ DIR, DEV, SHOT, DIM, ID, STIME
      COMMON /DATA/ DIR,DEV,SHOT,DIM,ID,STIME
C
      IDOPEN=25
      LINE='guiparm'
C
      INQUIRE(FILE=LINE,EXIST=LEX,ERR=9100)
      IF(.NOT.LEX) THEN
         WRITE(6,'(A)') '## NO INPUT FILE EXISTS.'
         RETURN
      ENDIF
      OPEN(IDOPEN,FILE=LINE,IOSTAT=IST1,STATUS='OLD',ERR=9100)
      READ(IDOPEN,GUIPARM,IOSTAT=IST2,ERR=9200,END=9900)
      CLOSE(IDOPEN)
C
      WRITE(6,'(A,A,A)') 
     &     '## FILE (',TRIM(LINE),') IS ASSIGNED FOR PARM INPUT'
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
      IMPLICIT NONE
      INTEGER NTM,NXM,NXMAX,NGMAX,MODE,I
      REAL GX1,GX2,GY1,GY2,GXMIN,GXMAX,GYMIN,GYMAX
      REAL GSXMIN,GSXMAX,GSYMIN,GSYMAX,GSTEPX,GSTEPY
      REAL GSTMIN,GSTMAX,GSTEPT,GXL,GYL,GZL,GPHI,GTHETA,GRADIUS
      REAL GOX,GOY,GOZ
      REAL*8 R2G2B
      REAL GX(NXM),GY(NXM,NGMAX),GT(NTM)
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
C     *** MAIN SUBROUTINE ***
C
      SUBROUTINE MAIN
      USE libchar
      IMPLICIT NONE
      COMMON /MDSPLUS/ MDS
      COMMON /DATA/ KDIR,KDEV,KNUM,NDIM,KFID,STIME
      INTEGER MDS,IFNT,IKDEV,IKFID,IKNUM
      INTEGER NDIM
      REAL W, H, STIME
      CHARACTER CHOICE(4)*256
      CHARACTER*80 KDIR,KDEV,KFID,KNUM
      EXTERNAL CLICKED,SETDEV,SETNUM,SETDIM,SETID,SETMODE,SETSTIME,
     &         END_MAIN
C
      CHOICE(1) = '@UFILE@'
      CHOICE(2) = '@MDSPLUS pr98@'
      CHOICE(3) = '@MDSPLUS wrk@'
      CHOICE(4) = '@MDSPLUS itb@'
C
      MDS  = 0
      IFNT = 2
      KDEV = 'jet'
      KNUM = '38285'
      NDIM  = 2
      KFID = 'te'
      KDIR = 'data'
C
      CALL parameter_read
C
      CALL KTRIM(KDEV,IKDEV)
      CALL KTRIM(KNUM,IKNUM)
      CALL KTRIM(KFID,IKFID)
C
      CALL SETFNT (IFNT)
      CALL SETCHR (0.35, 0.18, 0.35, 0.0, 0.0)
C
      CALL GUIEXTERNAL (0)
C
      CALL INQLABELSIZE ('@DEVICE = @', W, H)
      CALL GUILABEL (0.5, 17.0, W, 1.0, '@DEVICE = @', 0)
      CALL GUILINEEDIT (0.6+W, 17.0, 2.0, 1.0, 10, KDEV, IKDEV, SETDEV)
C
      CALL INQLABELSIZE ('@SHOT = @', W, H)
      CALL GUILABEL (5.0, 17.0, W, 1.0, '@SHOT = @', 0)
      CALL GUILINEEDIT (5.1+W, 17.0, 3.0, 1.0, 10, KNUM, IKNUM, SETNUM)
C
      CALL INQLABELSIZE ('@DIM. = @', W, H)
      CALL GUILABEL (10.2, 17.0, W, 1.0, '@DIM. = @', 0)
      CALL GUISPINBUTTONI (10.3+W, 17.0, 1.0, 1.0, NDIM, 1, 2,
     &     SETDIM)
C 
      CALL INQLABELSIZE ('@ID = @', W, H)
      CALL GUILABEL (13.1, 17.0, W, 1.0, '@ID = @', 0)
      CALL GUILINEEDIT (13.2+W, 17.0, 2.5, 1.0, 6, KFID, IKFID, SETID)
C
      CALL GUICOMBOBOX (17.2, 17.0, 3.5, 1.0, 4, CHOICE, 1, SETMODE)
C
      CALL INQLABELSIZE ('@SLICE = @', W, H)
      CALL GUILABEL (21.0, 0.5, W, 1.0, '@SLICE = @', 0)
      CALL GUISPINBUTTONF (21.1+W, 0.5, 2.5, 1.0, STIME, 0.0, 100.0, 3,
     &                     SETSTIME)
C
      CALL GUIBUTTON (21.0, 17.0, 2.0, 1.0, '@DRAW@', CLICKED)
      CALL GUIBUTTON (23.0, 17.0, 2.0, 1.0, '@END@', END_MAIN)
C
      RETURN
      END
C
C     *** ONCLICK PROCESS ON 'DRAW' BUTTON ***
C
      SUBROUTINE CLICKED
      USE libchar
      IMPLICIT NONE
      COMMON /MDSPLUS/ MDS
      COMMON /DATA/ KDIR,KDEV,KNUM,NDIM,KFID,STIME
      CHARACTER*80 KDIR,KDEV,KDEV1,KFID,KFID1,KNUM,KNUM1
      INTEGER MDS,NDIM,NDIM1
      INTEGER I,IKFID,IERR
      REAL STIME,STIME1
C
      CALL ERAS
C
      KDEV1=KDEV
      KNUM1=KNUM
      NDIM1=NDIM
      CALL KTRIM(KFID,IKFID)
      DO I=1,IKFID
         CALL toupper(KFID(I:I))
      ENDDO
      KFID1=KFID
      STIME1=STIME
C
      CALL SETDIR(KDIR,KDEV1,KNUM1,NDIM1,KFID1,IERR)
      IF(IERR.NE.0) RETURN
      CALL VIEW_GRAPH (NDIM1,STIME1)
C
      RETURN
      END
C
C     *** ONCLICK PROCESS ON 'END' BUTTON ***
C
      SUBROUTINE END_MAIN
      IMPLICIT NONE
      INTEGER MDS
      COMMON /MDSPLUS/ MDS
C
      IF(MDS.NE.0) CALL IPDB_CLOSE
      CALL GUIEXITLOOP
C
      RETURN
      END
C
C     *** DETERMINATING READING DIRECTORY ***
C
      SUBROUTINE SETDIR(KDIR,KDEV1,KDCG1,IDM,KFID1,IERR)
C
      USE libchar
      IMPLICIT NONE
      INTEGER MDS,IERR,IKDIR,IKDEV,IKDCG,IKDIRX,IDM,IKFID,IKDIRR1,
     &        IKDIRR2,IKFILE
      COMMON /MDSPLUS/ MDS
      COMMON /KID1/ KDEV,KDCG,KFID
      COMMON /TRKID1/ KDIRR1,KDIRR2
      CHARACTER KDEV*80,KDCG*80,KFID*80
      CHARACTER KDEV1*80,KDCG1*80,KFID1*80
      CHARACTER KDIRR1*80,KDIRR2*80
      CHARACTER KDIR*80,KDIRX*90,KFIDCK*90,KFILE*110
      CHARACTER KERR*100,KADD*10
      LOGICAL LEX
      INTEGER NSTR,NSTR1
C
      IERR=0
      KDEV=KDEV1
      KDCG=KDCG1
      KFID=KFID1
C
      CALL KTRIM(KDIR,IKDIR)
      CALL KTRIM(KDEV,IKDEV)
      CALL KTRIM(KDCG,IKDCG)
      IF(MDS.NE.0) THEN
         IF(MDS.EQ.1) THEN
            KDEV='pr98_'
         ELSEIF(MDS.EQ.2) THEN
            KDEV=''
         ELSEIF(MDS.EQ.3) THEN
            KDEV='itb_'
         ENDIF
         CALL KTRIM(KDEV,NSTR)
         CALL KTRIM(KDEV1,NSTR1)
         CALL APSTOS(KDEV,NSTR,KDEV1,NSTR1)
C
         CALL IPDB_OPEN(KDEV, KDCG)
         RETURN
      ENDIF
C
      KDIRX=KDIR(1:IKDIR)//'/'//KDEV(1:IKDEV)//'/'//KDCG(1:IKDCG)
     &     //'/in/'
      CALL KTRIM(KDIRX,IKDIRX)
      KFILE=KDIRX(1:IKDIRX)//KDEV(1:IKDEV)//'2d'//KDCG(1:IKDCG)//
     &     '.NE'
      INQUIRE(FILE=KFILE,EXIST=LEX,ERR=9000)
      IF(LEX.EQV..FALSE.) THEN
         CALL KTRIM(KDIRX,IKDIRX)
         KERR='@XX: NONEXISTENT DIRECTORY = '//KDIRX(1:IKDIRX)//'@'
         CALL GTEXTX(10.0,10.0,KERR,2)
         IERR=1
         GOTO 9000
      ENDIF
C
      KDIRR1=KDIRX(1:IKDIRX)//KDEV(1:IKDEV)
     &       //'1d'//KDCG(1:IKDCG)//'.'
      KDIRR2=KDIRX(1:IKDIRX)//KDEV(1:IKDEV)
     &       //'2d'//KDCG(1:IKDCG)//'.'
C
      CALL KTRIM(KDIRR1,IKDIRR1)
      CALL KTRIM(KDIRR2,IKDIRR2)
      CALL KTRIM(KFID,IKFID)
      IF(IDM.EQ.1) THEN
         KFIDCK=KDIRR1(1:IKDIRR1)//KFID(1:IKFID)
         KFILE =KDEV(1:IKDEV)//'1d'//KDCG(1:IKDCG)//'.'//KFID(1:IKFID)
      ELSEIF(IDM.EQ.2) THEN
         KFIDCK=KDIRR2(1:IKDIRR2)//KFID(1:IKFID)
         KFILE =KDEV(1:IKDEV)//'2d'//KDCG(1:IKDCG)//'.'//KFID(1:IKFID)
      ELSE
         GOTO 9000
      ENDIF
      INQUIRE(FILE=KFIDCK,EXIST=LEX,ERR=9000)
      IF(LEX.EQV..FALSE.) THEN
         CALL KTRIM(KFILE,IKFILE)
         KERR='@XX: NONEXISTENT FILE = '//KFILE(1:IKFILE)//'@'
         CALL GTEXTX(10.0,10.0,KERR,2)
         IERR=1
      ENDIF
C
 9000 RETURN
      END
C
C     *** SHOWING GRAPH ***
C
      SUBROUTINE VIEW_GRAPH(IDM,STIME)
C
      USE libitp
      USE libchar
      IMPLICIT NONE
      INTEGER NRM,NTM
      PARAMETER (NRM=100,NTM=2001)
      INTEGER MDS,NTX,NTXMAX,NRX,NRXMAX,NTSL,IERR,NTXMAX1,NRL,NRLMAX,NL
      INTEGER IDM,IKDEV,IKDCG,IKFID
      REAL GX1,GX2,GY1,GY2,GUCLIP,STIME
      REAL*8 FACT
      CHARACTER KDEV*80,KDCG*80,KFID*80,KTIME*12,KVAL*80
      CHARACTER KFIDX*80,KVAR*80
      CHARACTER KDIRR1*80,KDIRR2*80
      CHARACTER KERR*100,KDEVL*80
      REAL*8 T(NTM),R(NRM),F1(NTM),F2(NRM,NTM),F3(NTM,NRM)
      REAL*8 RL(NRM),FL(NRM),RP,FP
      REAL GT(NTM),GR(NRM),GF1(NTM),GF2(NRM,NTM)
      COMMON /MDSPLUS/ MDS
      COMMON /KID1/ KDEV,KDCG,KFID
      COMMON /TRKID1/ KDIRR1,KDIRR2
C
      IERR=0
C
      IF(IDM.EQ.1) THEN
         IF(MDS.EQ.0) THEN
            CALL TRXR1D(KDIRR1,KFID,T,F1,NTM,NTXMAX,1)
            IF(NTXMAX.EQ.0) IERR=1
         ELSE
            CALL IPDB_MDS1(KDEV,KDCG,KFID,NTM,T,F1,NTXMAX,IERR)
         ENDIF
      ELSE IF(IDM.EQ.2) THEN
         IF(MDS.EQ.0) THEN
            CALL TRXR2D(KDIRR2,KFID,T,R,F2,NRM,NTM,NRXMAX,NTXMAX,1)
            IF(NTXMAX.EQ.0) IERR=1
         ELSE
            CALL IPDB_MDS2(KDEV,KDCG,KFID,NRM,NTM,
     &                     R,T,F3,NRXMAX,NTXMAX,IERR)
            DO NTX=1,NTXMAX
               DO NRX=1,NRXMAX
                 F2(NRX,NTX)=F3(NTX,NRX) 
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      IF(IERR.NE.0) THEN
         KERR='@XX: DEVICE, NONEXISTENT DISCHARGE OR FILE ID@'
         CALL GTEXTX(10.0,10.0,KERR,2)
         GOTO 9000
      ENDIF
C
      GX1=2.5
      GX2=17.0
      GY1=1.5
      GY2=16.5
C
C     *** F-t graph ***
      IF(IDM.EQ.1) THEN
         IF(NTXMAX.EQ.1) THEN
            WRITE(KVAL,'(A13)') '@** VALUE **@'
            CALL GTEXTX(18.0,15.0,KVAL,0)
            WRITE(KVAL,'(A9,F15.7,A1)') '@TIME  = ',REAL(T(1)),'@'
            CALL GTEXTX(18.5,14.3,KVAL,0)
            WRITE(KVAL,'(A9,F15.7,A1)') '@VALUE = ',REAL(F1(1)),'@'
            CALL GTEXTX(18.5,13.7,KVAL,0)
         ELSE
            DO NTX=1,NTXMAX
               GT(NTX)=GUCLIP(T(NTX))
               GF1(NTX)=GUCLIP(F1(NTX))
            ENDDO
            CALL KTRIM(KDEV,IKDEV)
            CALL KTRIM(KDCG,IKDCG)
            CALL KTRIM(KFID,IKFID)
C            KFIDX='@'//KDEV(1:IKDEV)//'/'//KDCG(1:IKDCG)//'/'
C     &           //KFID(1:IKFID)//'@'
            KFIDX=' '
            CALL TRGR1D(GX1,GX2,GY1,GY2,
     &                  GT,GF1,NTM,NTXMAX,1,KFIDX,2)
         ENDIF
      ELSE
C
         CALL label_of_slice_time(STIME,T,NTM,NTXMAX,NTSL,IERR)
         IF(IERR.NE.0) RETURN
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
            WRITE(KVAL,'(A21)') '@RHO           VALUE@'
            CALL GTEXTX(18.5,14.3,KVAL,0)
            DO NL=1,11
               RP=(NL-1)*0.1D0
               CALL AITKEN(RP,FP,RL,FL,2,NRLMAX)
               WRITE(KVAL,'(A1,F6.4,A3,F15.7,A1)') '@',RP,'   ',FP,'@'
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
C
 9000 RETURN
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
C     *** SUBSTITUTING PARAMETERS INTO COMMON BLOCK ***
C
      SUBROUTINE SETDEV(KDEV1,L)
      IMPLICIT NONE
      COMMON /DATA/ KDIR,KDEV,KNUM,NDIM,KFID,STIME
      CHARACTER*80 KDIR,KDEV,KDEV1,KNUM,KFID
      INTEGER NDIM
      INTEGER I,L
      REAL STIME
C
      KDEV = ' '
      DO I = 1, L
         KDEV(I:I) = KDEV1(I:I) 
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE SETNUM(KNUM1,L)
      IMPLICIT NONE
      COMMON /DATA/ KDIR,KDEV,KNUM,NDIM,KFID,STIME
      CHARACTER*80 KDIR,KDEV,KNUM,KNUM1,KFID
      INTEGER NDIM
      INTEGER I,L
      REAL STIME
C
      KNUM = ' '
      DO I = 1, L
         KNUM(I:I) = KNUM1(I:I) 
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE SETDIM(NDIM1)
      IMPLICIT NONE
      COMMON /DATA/ KDIR,KDEV,KNUM,NDIM,KFID,STIME
      CHARACTER*80 KDIR,KDEV,KNUM,KFID
      INTEGER NDIM,NDIM1
      REAL STIME
C
      NDIM=NDIM1
C
      RETURN
      END
C
      SUBROUTINE SETID(KFID1,L)
      IMPLICIT NONE
      COMMON /DATA/ KDIR,KDEV,KNUM,NDIM,KFID,STIME
      CHARACTER*80 KDIR,KDEV,KNUM,KFID,KFID1
      INTEGER NDIM
      INTEGER I,L
      REAL STIME
C
      KFID = ' '
      DO I = 1, L
         KFID(I:I) = KFID1(I:I) 
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE SETMODE(MDS1)
      IMPLICIT NONE
      COMMON /MDSPLUS/ MDS
      INTEGER MDS, MDS1
C
      MDS=MDS1-1
C
      RETURN
      END
C
      SUBROUTINE SETSTIME(STIME1)
      IMPLICIT NONE
      COMMON /DATA/ KDIR,KDEV,KNUM,NDIM,KFID,STIME
      CHARACTER*80 KDIR,KDEV,KNUM,KFID
      INTEGER NDIM
      REAL STIME,STIME1
C
      STIME=STIME1
C
      RETURN
      END
C
C     *****************************************************************
C
C     SUBROUTINE APpend Strings TO Strings
C        INPUT  : NSTR, INSTR, NINSTR
C                 NSTR : Number of STR. First, NSTR = 0.
C                 NINSTR : Number of INSTR
C        OUTPUT : STR, NSTR
C                 STR(NSTR(original+1):NSTR(return))
C
C     *****************************************************************
C
      SUBROUTINE APSTOS(STR, NSTR, INSTR, NINSTR)
      IMPLICIT REAL*8(A-F, H, O-Z)
      CHARACTER STR*(*), INSTR*(*)
C
      STR(NSTR+1:NSTR+NINSTR) = INSTR(1:NINSTR)
      NSTR = NSTR + NINSTR
C
      RETURN
      END
