C     $Id$
C  
C     ***********************************************************
C
C           GRAPHIC 3D : CONTROL ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TRGRD0(K2,K3,INQ)
C
      INCLUDE 'trcomm.h'
      CHARACTER K2*1,K3*3,KK*4
      CHARACTER STRL*80,KVL*80
C
      READ(K2,'(I1)',ERR=600) I2
      READ(K3,'(I1)') I3
      IF (K3.EQ.' ') THEN
         NMB=I2
         IF (K2.EQ.' ') THEN
            CALL VIEW3DLIST
            GOTO 100
         ENDIF
      ELSE
         NMB=I2*10+I3
      ENDIF
      GOTO 200
C
 600  KK=K2//K3
      CALL CHNBPR(KK,NMB,IERR)
      IF (IERR.EQ.1) GOTO 100
C
 200  CALL TRGRTD(STRL,KVL,NMB)
      CALL TRGRUR(G3D(1,1,NMB),STRL,KVL,INQ)
C
 100  RETURN
      END
C
C     **************************************************************
C
C           GRAPHIC 3D : UNIVERSAL ROUTINE
C
C     **************************************************************
C
      SUBROUTINE TRGRUR(GVD,STR,KV,INQ)
C
      INCLUDE 'trcomm.h'
      CHARACTER STR*80,KV*80
      DIMENSION GVD(NRM,NTM)
C
      GX1=3.0
      GX2=18.0
      GY1=2.0
      GY2=17.0
C
      CALL PAGES
      CALL TRGR3D(GX1,GX2,GY1,GY2,GRM,GT,GVD,NRM,NRMAX,NGT,
     &            STR,KV,2+INQ)
      CALL PAGEE
C
      RETURN
      END
C
C     **************************************************************
C
C           GRAPHIC 3D : GRAPH TITLE DATA
C
C     **************************************************************
C
      SUBROUTINE TRGRTD(STRL,KVL,NMB)
C
      INCLUDE 'trcomm.h'
      CHARACTER STRL*80,KVL*80
      CHARACTER STR0(29)*80,KV0(29)*80
C
      DATA STR0/'@TE [keV] vs t@','@TD [keV] vs t@',
     &          '@TT [keV] vs t@','@TA [keV] vs t@',
     &          '@NE [10^20/m^3] vs t@','@ND [10^20/m^3] vs t@',
     &          '@NT [10^20/m^3] vs t@','@NA [10^20/m^3] vs t@',
     &          '@AJ [MA/m^2] vs t@','@AJOH [MA/m^2] vs t@',
     &          '@AJNB [MA/m^2] vs t@','@AJRF [MA/m^2] vs t@',
     &          '@AJBS [MA/m^2] vs t@',
     &          '@PIN [MW/m^3] vs t@','@POH [MW/m^3] vs t@',
     &          '@PNB [MW/m^3] vs t@','@PNF [MW/m^3] vs t@',
     &          '@PRFE [MW/m^3] vs t@','@PRFD [MW/m^3] vs t@',
     &          '@PRFT [MW/m^3] vs t@','@PRFA [MW/m^3] vs t@',
     &          '@PRL [MW/m^3] vs t@','@PCX [MW/m^3] vs t@',
     &          '@PIE [MW/m^3] vs t@',
     &          '@QP vs t@','@EZOH [V/m] vs t@',
     &          '@BETA vs t@','@BETAP vs t@','@VLOOP [V] vs t@'/
C
      DATA KV0 /'@TE@','@TD@','@TT@','@TA@',
     &          '@NE@','@ND@','@NT@','@NA@',
     &          '@AJ@','@AJOH@','@AJNB@','@AJRF@','@AJBS@',
     &          '@PIN@','@POH@','@PNB@','@PNF@',
     &          '@PRFE@','@PRFD@','@PRFT@','@PRFA@',
     &          '@PRL@','@PCX@','@PIE@',
     &          '@QP@','@EZOH@','@BETA@','@BETAP@','@VLOOP@'/
      SAVE STR0, KV0
C
      STRL=STR0(NMB)
      KVL =KV0 (NMB)
C
      RETURN
      END
C
C     **************************************************************
C
C           GRAPHIC 3D : PARAMETER STRING FOR NP
C
C     **************************************************************
C
      SUBROUTINE CHKVPL(KVPL,NP)
C
      CHARACTER KVPL*4,KVP(29)*4
C
      DATA KVP /'TE  ','TD  ','TT  ','TA  ',
     &          'NE  ','ND  ','NT  ','NA  ',
     &          'AJ  ','AJOH','AJNB','AJRF','AJBS',
     &          'PIN ','POH ','PNB ','PNF ',
     &          'PRFE','PRFD','PRFT','PRFA',
     &          'PRL ','PCX ','PIE ',
     &          'QP  ','EZOH','BETA','BETP','VLOP'/
C
      KVPL=KVP(NP)
      RETURN
      END
C
C     **************************************************************
C
C           GRAPHIC 3D : CHECK A NUMBER CORRESPONDING A PARAMETER
C
C     **************************************************************
C
      SUBROUTINE CHNBPR(KK,NMB,IERR)
C
      CHARACTER KK*4,KVPL*4
C
      IERR=0
      DO NP=1,29
         CALL CHKVPL(KVPL,NP)
         IF (KK.EQ.KVPL) THEN
            NMB=NP
            GOTO 1000
         ENDIF
      ENDDO
      IERR=1
C
 1000 RETURN
      END
C
C     **************************************************************
C
C           GRAPHIC 3D : VIEW 3D GRAPHICS LIST 
C
C     **************************************************************
C
      SUBROUTINE VIEW3DLIST
C
      CHARACTER KK*4,KVP(88)*4
      COMMON KVP
C
      WRITE(6,700)
      WRITE(6,710) (KVP(I),I= 1, 5)
      WRITE(6,720) (KVP(I),I= 6,10)
      WRITE(6,730) (KVP(I),I=11,15)
      WRITE(6,740) (KVP(I),I=16,20)
      WRITE(6,750) (KVP(I),I=21,25)
      WRITE(6,760) (KVP(I),I=26,29)
C
 700  FORMAT(' ',
     &       '# 3D GRAPHICS LIST (EACH ONE IS CONSIST OF 4 LETTERS)')
 710  FORMAT(' ',' 1- 5: ',4(A4,', '),A4)
 720  FORMAT(' ',' 6-10: ',4(A4,', '),A4)
 730  FORMAT(' ','11-15: ',4(A4,', '),A4)
 740  FORMAT(' ','16-20: ',4(A4,', '),A4)
 750  FORMAT(' ','21-25: ',4(A4,', '),A4)
 760  FORMAT(' ','26-29: ',3(A4,', '),A4)
C
      RETURN
      END
C
C     ***********************************************************
C
C           SUBPROGRAM FOR 3D PROFILE
C
C     ***********************************************************
C
      SUBROUTINE TRGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,
     &     NXM,NXMAX,NYMAX,STR,KV,MODE)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
      DIMENSION GX(NXMAX),GY(NYMAX),GZ(NXM,NYMAX)
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
      CALL GMNMX2(GZ,NXM,1,NXMAX,1,1,NYMAX,1,GZMIN,GZMAX)
      IF(ABS(GZMAX-GZMIN).LT.1.D-6) THEN
         GZMIN=GZMIN-0.999E-6
         GZMAX=GZMAX+1.000E-6
      ENDIF
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
      CALL GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GSTEPZ)
C
      IF(MOD(MODE,2).EQ.0) THEN
         IF(GZMIN.GE.0.0) THEN
            GSZMIN=0.0
         ELSEIF(GZMAX.LE.0.0) THEN
            GSZMAX=0.0
         ENDIF
      ENDIF
      GZMIN=GSZMIN
      GZMAX=GSZMAX
      IF(MOD(MODE/4,2).EQ.1) THEN
         CALL CHMODE
         WRITE(6,*) '## TRGR : XMIN,XMAX,YMIN,YMAX = ',
     &              GXMIN,GXMAX,GYMIN,GYMAX
         READ(5,*) GXMIN,GXMAX,GYMIN,GYMAX
         CALL GRMODE
      ENDIF
      CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
      CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSTEPY)
      CALL GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GSTEPZ)
C
C      IF(GXMIN*GXMAX.LE.0.0) THEN
C         GXORG=0.0
C      ELSE
C         GXORG=GSXMIN
C      ENDIF
C      IF(GYMIN*GYMAX.LE.0.0) THEN
C         GYORG=0.0
C      ELSE
C         GYORG=GSYMIN
C      ENDIF
C
      CALL GDEFIN(GX1,GX2,GY1,GY2,
     &            GSXMIN,GSXMAX,GSYMIN,GSYMAX)
      GXL=10.0*1.5
      GYL=20.0*1.5
      GZL=10.0*1.5
      CALL GDEFIN3D(GXL,GYL,GZL,GSXMIN,GSXMAX,GYMIN,GYMAX,
     &     GSZMIN,GSZMAX)
      GPHI=-20.0
      GTHETA=60.0
      GRADIUS=15.0
      GOX=0.5*(GSXMIN+GSXMAX)
      GOY=0.5*(GYMIN+GYMAX)
      GOZ=0.5*(GSZMIN+GSZMAX)
      CALL GVIEW3D(GPHI,GTHETA,GRADIUS,GOX,GOY,GOZ)
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(0,0,4)
C
C      CALL GSCALE3DX(GSXMIN,GSTEPX,0.3,0)
C      CALL GSCALE3DY(GT(1),GSTEPT,0.3,0)
C      CALL GSCALE3DZ(GSYMIN,GSTEPY,0.3,10)
C      CALL GVALUE3DX(GSXMIN,GSTEPX,1,1)
C      CALL GVALUE3DY(GT(1),GSTEPT,1,1)
C      CALL GVALUE3DZ(GSYMIN,GSTEPY,11,-2)
      CALL GSCALE3DX(GSXMIN,GSTEPX,0.3,0)
      CALL GSCALE3DY(GSYMIN,GSTEPY,0.3,0)
      CALL GSCALE3DZ(GSZMIN,GSTEPZ,0.3,10)
      CALL GVALUE3DX(GSXMIN,GSTEPX,1,1)
      CALL GVALUE3DY(GSYMIN,GSTEPY,1,1)
      CALL GVALUE3DZ(GSZMIN,GSTEPZ,11,-2)
C
      CALL Set3DTextBaseLine(0.0, 1.0, 0.0, -1.0, 0.0, 0.0)
      CALL GTEXTX3D(GSXMAX+0.15*(GSXMAX-GSXMIN),
     &              0.5*(GY(1)+GY(NYMAX)),
     &              GSYMIN,
     &              '@TIME (sec)@',
     &              2)
      CALL Set3DTextBaseLine(0.0, 1.0, 0.0, -1.0, 0.0, 0.0)
      CALL GTEXTX3D(0.5*(GSXMIN+GSXMAX),
     &              GY(1)+0.1*(GY(1)-GY(NYMAX)),
     &              GSYMIN,
     &              '@R@',
     &              2)
      CALL Set3DTextBaseLine(0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
      CALL GTEXTX3D(GSXMIN,
     &              GY(1)+0.05*(GY(1)-GY(NYMAX)),
     &              GSZMAX+0.1*(GSZMAX-GSZMIN),
     &              KV,
     &              2)
C
      CALL PERS3D1(GZ,NXM,NXMAX,NYMAX,-27,R2G2B)
      CALL GAxis3D(0)
      CALL GDrawBack3D(0.5, 0.5, 0.5)
C
      CALL SETLIN(0,0,4)
      RETURN
      END
