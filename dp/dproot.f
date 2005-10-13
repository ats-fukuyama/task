C     $Id$
C
C     *************************
C           PLOT 1D GRAPH      
C     *************************
C
      SUBROUTINE DPGRP1
C
      INCLUDE 'dpcomm.inc'
      DIMENSION GX(NXM),GY(NXM)
C
    1 WRITE(6,*)'## SELECT : (1)RF (2)RFI (3)KX (4)KY (5)KZ (0)END'
      READ(5,*,ERR=1,END=9000) ISW
      IF(ISW.EQ.0) GOTO 9000
C
  100 CONTINUE
C
      IF(ISW.EQ.1) THEN
    2    WRITE(6,*) ' INPUT : RF1,RF2,NXMAX'
         NXMAX1=NXMAX
         READ(5,*,ERR=2,END=1) RF1,RF2,NXMAX1
         IF(NXMAX1.GT.NXM) GOTO 2
         IF(NXMAX1.EQ.0) GOTO 1
         NXMAX=NXMAX1
         DRF=(RF2-RF1)/NXMAX
         DO 10 NX=1,NXMAX
            CRF=DCMPLX(RF1,RFI0)+DRF*(NX-1)
            CKX=RKX0
            CKY=RKY0
            CKZ=RKZ0
            CD=CFDISP(CRF,CKX,CKY,CKZ,RX0,RY0,RZ0)
            GX(NX)=GUCLIP(DBLE(CRF))
            GY(NX)=GUCLIP(LOG10(ABS(CD)))
   10    CONTINUE
      ELSEIF(ISW.EQ.2) THEN
  102    WRITE(6,*) ' INPUT : RFI1,RFI2,NXMAX'
         NXMAX1=NXMAX
         READ(5,*,ERR=102,END=1) RFI1,RFI2,NXMAX1
         IF(NXMAX1.GT.NXM) GOTO 102
         IF(NXMAX1.EQ.0) GOTO 1
         NXMAX=NXMAX1
         DRFI=(RFI2-RFI1)/NXMAX
         DO 110 NX=1,NXMAX
            CRF=DCMPLX(RF0,RFI1)+DCMPLX(0.D0,DRFI)*(NX-1)
            CKX=RKX0
            CKY=RKY0
            CKZ=RKZ0
            CD=CFDISP(CRF,CKX,CKY,CKZ,RX0,RY0,RZ0)
            GX(NX)=GUCLIP(DIMAG(CRF))
            GY(NX)=GUCLIP(LOG10(ABS(CD)))
  110    CONTINUE
      ELSEIF(ISW.EQ.3) THEN
  202    WRITE(6,*) ' INPUT : RKX1,RKX2,NXMAX,RF0,RX0,RNPHII'
         NXMAX1=NXMAX
         READ(5,*,ERR=202,END=1) RKX1,RKX2,NXMAX1,RF0,RX0,RNPHII
         IF(NXMAX1.GT.NXM) GOTO 202
         IF(NXMAX1.EQ.0) GOTO 1
         NXMAX=NXMAX1
         DKX=(RKX2-RKX1)/NXMAX
         DO 210 NX=1,NXMAX
            CRF=DCMPLX(RF0,RFI0)
            CKX=RKX1+DKX*(NX-1)
            RKY0=2.D6*PI*RF0*RNPHII/VC
            CKY=RKY0
            CKZ=RKZ0
            CD=CFDISP(CRF,CKX,CKY,CKZ,RX0,RY0,RZ0)
C            write(6,*) 
C            write(6,*)crf,ckx,cky,ckz,rx0,ry0,rz0,CD
            GX(NX)=GUCLIP(DBLE(CKX))
            GY(NX)=GUCLIP(LOG10(ABS(CD)))
  210    CONTINUE
      ELSEIF(ISW.EQ.4) THEN
  302    WRITE(6,*) ' INPUT : RKY1,RKY2,NXMAX'
         NXMAX1=NXMAX
         READ(5,*,ERR=302,END=1) RKY1,RKY2,NXMAX1
         IF(NXMAX1.GT.NXM) GOTO 302
         IF(NXMAX1.EQ.0) GOTO 1
         NXMAX=NXMAX1
         DKY=(RKY2-RKY1)/NXMAX
         DO 310 NX=1,NXMAX
            CRF=DCMPLX(RF0,RFI0)
            CKX=RKX0
            CKY=RKY1+DKY*(NX-1)
            CKZ=RKZ0
            CD=CFDISP(CRF,CKX,CKY,CKZ,RX0,RY0,RZ0)
            GX(NX)=GUCLIP(DBLE(CKY))
            GY(NX)=GUCLIP(LOG10(ABS(CD)))
  310    CONTINUE
      ELSEIF(ISW.EQ.5) THEN
  402    WRITE(6,*) ' INPUT : RKZ1,RKZ2,NXMAX'
         NXMAX1=NXMAX
         READ(5,*,ERR=402,END=1) RKZ1,RKZ2,NXMAX1
         IF(NXMAX1.GT.NXM) GOTO 402
         IF(NXMAX1.EQ.0) GOTO 1
         NXMAX=NXMAX1
         DKZ=(RKZ2-RKZ1)/NXMAX
         DO 410 NX=1,NXMAX
            CRF=DCMPLX(RF0,RFI0)
            CKX=RKX0
            CKY=RKY0
            CKZ=RKZ1+DKZ*(NX-1)
            CD=CFDISP(CRF,CKX,CKY,CKZ,RX0,RY0,RZ0)
            GX(NX)=GUCLIP(DBLE(CKZ))
            GY(NX)=GUCLIP(LOG10(ABS(CD)))
  410    CONTINUE
      ELSE
         GOTO 1
      ENDIF
C
      CALL PAGES
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.3,0.)
      CALL GMNMX1(GX,1,NXMAX,1,GXMIN,GXMAX)
      CALL GMNMX1(GY,1,NXMAX,1,GYMIN,GYMAX)
      CALL GQSCAL(GXMIN,GXMAX,GXSMN,GXSMX,GSCALX)
      CALL GQSCAL(GYMIN,GYMAX,GYSMN,GYSMX,GSCALY)
      CALL GDEFIN(3.,18.,2.,17.,GXSMN,GXSMX,GYSMN,GYSMX)
      CALL GFRAME
      CALL GSCALE(GXSMN,GSCALX,0.0,0.0,0.2,9)
      CALL GSCALL(0.0,0,GYSMN,2,0.2,9)
      CALL GVALUE(GXSMN,2*GSCALX,0.0,0.0,-2)
      CALL GVALUL(0.0,0,GYSMN,1,0)
      CALL GPLOTP(GX,GY,1,NXMAX,1,0,0,0)
C
      CALL MOVE(9.0,18.0)
      IF (ISW.EQ.1) THEN
       CALL TEXT (' RF ',4)
      ELSEIF (ISW.EQ.2) THEN
       CALL TEXT (' RFI',4) 
      ELSEIF (ISW.EQ.3) THEN
       CALL TEXT (' KPP',4) 
      ELSEIF (ISW.EQ.4) THEN
       CALL TEXT (' KPR',4) 
      ENDIF
      CALL MOVE(20.0,16.5)
      CALL TEXT('RF  =',5)
      CALL NUMBD(RF0,  '(1PE11.3)',11)
      CALL MOVE(20.0,16.0)
      CALL TEXT('RFI0=',5)
      CALL NUMBD(RFI0, '(1PE11.3)',11)
      CALL MOVE(20.0,15.5)
      CALL TEXT('KX  =',5)
      CALL NUMBD(RKX0,'(1PE11.3)',11)
      CALL MOVE(20.0,15.0)
      CALL TEXT('KY  =',5)
      CALL NUMBD(RKY0,'(1PE11.3)',11)
      CALL MOVE(20.0,14.5)
      CALL TEXT('KZ  =',5)
      CALL NUMBD(RKZ0,'(1PE11.3)',11)
      CALL PAGEE
      GOTO 100
C
 9000 RETURN
      END
C
C     *************************
C             FIND ROOT
C     *************************
C
      SUBROUTINE DPROOT
C
      INCLUDE 'dpcomm.inc'
C
      PARAMETER (IGM=100)
      DIMENSION GX(NXM),GY(NXM),GY2(NXM)
      DIMENSION GXA(NXM,IGM),GYA(NXM,IGM),GY2A(NXM,IGM),NXMAXA(IGM)
C
    1 WRITE(6,*)'## SELECT : CRF : (1)ROOT (2)KX (3)KY (4)KZ (5) X'
      WRITE(6,*)'  (0) END : CKX : (6)ROOT (7)X1 (8)X2'
      READ(5,*,ERR=1,END=9000) ISW
      IF(ISW.EQ.0) GOTO 9000
      RFX=RF0
      RFIX=RFI0
C
      IF(ISW.EQ.1) THEN
    2    WRITE(6,*) ' INPUT : RFX,RFIX'
         READ(5,*,ERR=2,END=1) RFX,RFIX
            ILIST=1
            CKX0=RKX0
            CKY0=RKY0
            CKZ0=RKZ0
            XPOS0=RX0
            YPOS0=RY0
            ZPOS0=RZ0
            CRF=DCMPLX(RFX,RFIX)
            CALL DPBRENT1(CRF)
            WRITE(6,*) CRF
            GOTO 2
      ELSEIF(ISW.EQ.2) THEN
  202    WRITE(6,*) ' INPUT : RKX1,RKX2,NXMAX'
         READ(5,*,ERR=202,END=1) RKX1,RKX2,NXMAX
         IF(NXMAX.GT.NXM) GOTO 202
         DKX=(RKX2-RKX1)/(NXMAX-1)
         CRF=DCMPLX(RF0,RFI0)
         CKY0=RKY0
         CKZ0=RKZ0
         XPOS0=RX0
         YPOS0=RY0
         ZPOS0=RZ0
         DO 210 NX=1,NXMAX
            ILIST=0
            CKX0=RKX1+DKX*(NX-1)
            WRITE(6,*) CRF
            CALL DPBRENT1(CRF)
            GX(NX)=GUCLIP(DBLE(CKX0))
            GY(NX)=GUCLIP(DBLE(CRF))
            GY2(NX)=GUCLIP(DIMAG(-CRF))
  210    CONTINUE
         CALL DPGR1A(GX,GY,GY2,NXMAX,ISW)
      ELSEIF(ISW.EQ.3) THEN
  302    WRITE(6,*) ' INPUT : RKY1,RKY2,NXMAX'
         READ(5,*,ERR=302,END=1) RKY1,RKY2,NXMAX
         IF(NXMAX.GT.NXM) GOTO 302
         DKY=(RKY2-RKY1)/(NXMAX-1)
         CRF=DCMPLX(RF0,RFI0)
         CKX0=RKX0
         CKZ0=RKZ0
         XPOS0=RX0
         YPOS0=RY0
         ZPOS0=RZ0
         DO 310 NX=1,NXMAX
            ILIST=0
            CKY0=RKY1+DKY*(NX-1)
            CALL DPBRENT1(CRF)
            WRITE(6,*) CRF
            GX(NX)=GUCLIP(DBLE(CKY0))
            GY(NX)=GUCLIP(DBLE(CRF))
            GY2(NX)=GUCLIP(DIMAG(-CRF))
  310    CONTINUE
         CALL DPGR1A(GX,GY,GY2,NXMAX,ISW)
      ELSEIF(ISW.EQ.4) THEN
  402    WRITE(6,*) ' INPUT : RKZ1,RKZ2,NXMAX'
         READ(5,*,ERR=402,END=1) RKZ1,RKZ2,NXMAX
         IF(NXMAX.GT.NXM) GOTO 402
         DKZ=(RKZ2-RKZ1)/(NXMAX-1)
         CRF=DCMPLX(RF0,RFI0)
         CKX0=RKX0
         CKY0=RKY0
         XPOS0=RX0
         YPOS0=RY0
         ZPOS0=RZ0
         DO 410 NX=1,NXMAX
            ILIST=0
            CKZ0=RKZ1+DKZ*(NX-1)
            CALL DPBRENT1(CRF)
            WRITE(6,*) CRF
            GX(NX)=GUCLIP(DBLE(CKZ0))
            GY(NX)=GUCLIP(DBLE(CRF))
            GY2(NX)=GUCLIP(DIMAG(-CRF))
  410    CONTINUE
         CALL DPGR1A(GX,GY,GY2,NXMAX,ISW)
      ELSEIF(ISW.EQ.5) THEN
  502    WRITE(6,*) ' INPUT : RX1,RX2,NXMAX'
         READ(5,*,ERR=502,END=1) RX1,RX2,NXMAX
         IF(NXMAX.GT.NXM) GOTO 502
         DXPOS=(RX2-RX1)/(NXMAX-1)
         CRF=DCMPLX(RF0,RFI0)
         CKX0=RKX0
         CKY0=RKY0
         CKZ0=RKZ0
         YPOS0=RY0
         ZPOS0=RZ0
         DO 510 NX=1,NXMAX
            ILIST=0
            XPOS0=RX1+DXPOS*(NX-1)
            CALL DPBRENT1(CRF)
            WRITE(6,*) CRF
            GX(NX)=GUCLIP(XPOS0)
            GY(NX)=GUCLIP(DBLE(CRF))
            GY2(NX)=GUCLIP(DIMAG(-CRF))
  510    CONTINUE
         CALL DPGR1A(GX,GY,GY2,NXMAX,ISW)
      ELSEIF(ISW.EQ.6) THEN
  602    WRITE(6,*) ' INPUT : RKX,RKXI'
         READ(5,*,ERR=602,END=1) RKXX,RKXIX
            ILIST=1
            CRF0=DCMPLX(RF0,0.D0)
            CKY0=RKY0
            CKZ0=RKZ0
            XPOS0=RX0
            YPOS0=RY0
            ZPOS0=RZ0
            CKX=DCMPLX(RKXX,RKXIX)
            CALL DPBRENT2(CKX)
            WRITE(6,*) CKX
            GOTO 602
      ELSEIF(ISW.EQ.7) THEN
  702    WRITE(6,*) ' INPUT : RX1,RX2,NXMAX'
         READ(5,*,ERR=702,END=1) RX1,RX2,NXMAX
         IF(NXMAX.GT.NXM) GOTO 702
         DXPOS=(RX2-RX1)/(NXMAX-1)
         CRF0=DCMPLX(RF0,0.D0)
         CKX0=RKX0
         CKY0=RKY0
         CKZ0=RKZ0
         YPOS0=RY0
         ZPOS0=RZ0
         DO 710 NX=1,NXMAX
            ILIST=0
            XPOS0=RX1+DXPOS*(NX-1)
            CALL DPBRENT2(CKX)
            GX(NX)=GUCLIP(DBLE(XPOS0))
            GY(NX)=GUCLIP(DBLE(CKX**2))
            GY2(NX)=GUCLIP(DIMAG(CKX**2))
            WRITE(6,*) CKX
  710    CONTINUE
         CALL DPGR1A(GX,GY,GY2,NXMAX,ISW)
      ELSEIF(ISW.EQ.8) THEN
         I=0
  802    WRITE(6,*) ' INPUT : RX1,RX2,NXMAX,RKX0'
         READ(5,*,ERR=802,END=8000) RX1,RX2,NXMAX,RKX0
         IF(RX1.EQ.RX2)   GO TO 8000
         IF(NXMAX.GT.NXM) GOTO 802
         I=I+1
         DXPOS=(RX2-RX1)/(NXMAX-1)
         CRF0=DCMPLX(RF0,0.D0)
         CKX=DCMPLX(RKX0,0.D0)
         CKY0=RKY0
         CKZ0=RKZ0
         YPOS0=RY0
         ZPOS0=RZ0
         NXMAXA(I)=NXMAX
         DO 810 NX=1,NXMAX
            ILIST=0
            XPOS0=RX1+DXPOS*(NX-1)
            CALL DPBRENT2(CKX)
            GXA(NX,I)=GUCLIP(XPOS0)
              IF(DBLE(CKX**2).GT.0) THEN
                GYA(NX,I)=GUCLIP(SQRT(DBLE(CKX**2)))
              ELSEIF(DBLE(CKX**2).LT.0) THEN         
                GYA(NX,I)=GUCLIP(-SQRT(-DBLE(CKX**2)))
              ENDIF
              IF(DIMAG(CKX**2).GT.0) THEN
                GY2A(NX,I)=GUCLIP(SQRT(DIMAG(CKX**2)))
              ELSEIF(DIMAG(CKX**2).LT.0) THEN         
                GY2A(NX,I)=GUCLIP(-SQRT(DIMAG(CKX**2)))
              ENDIF
            WRITE(6,*) CKX
  810    CONTINUE
         CALL DPGR1A(GXA(1,I),GYA(1,I),GY2A(1,I),NXMAX,ISW)
         GOTO 802
      ENDIF
C
      GOTO 1
C
 8000 CALL DPGR1B(GXA,GYA,GY2A,NXM,NXMAXA,I,ISW)
      GOTO 1
C
 9000 RETURN
      END
C
C     *************************
C           BRENT METHOD 1
C     *************************
C
      SUBROUTINE DPBRENT1(CX)
C
      INCLUDE 'dpcomm.inc'
C
      EXTERNAL DPFUNC1
C
      PARAMETER (NBRM=2,LWA=NBRM*(NBRM+3))
      DIMENSION X(NBRM),F(NBRM),WA(LWA)
C
      M = LMAXRT
      N = NBRM
      X(1)=DBLE(CX)
      X(2)=DIMAG(CX)
C
      CALL FBRENTN(DPFUNC1,N,X,F,EPSRT,INFO,WA,LWA)
      IF(INFO.GE.1.AND.INFO.LE.3) THEN
         X(1)=RF0
         X(2)=RFI0
      ENDIF
C
      CX=DCMPLX(X(1),X(2))
      RETURN
      END
C
C     *******************************
C           Determinant Function
C     *******************************
C
      SUBROUTINE DPFUNC1(N,X,F,I)
C
      INCLUDE 'dpcomm.inc'
C
      DIMENSION X(N),F(N)
C
      CRF=DCMPLX(X(1),X(2))
      CFUNC=CFDISP(CRF,CKX0,CKY0,CKZ0,XPOS0,YPOS0,ZPOS0)
      IF(I.EQ.1)THEN
         F(1)=DBLE(CFUNC)
      ELSE
         F(2)=DIMAG(CFUNC)
      END IF
      IF(ILIST.NE.0) WRITE(6,'(1P3E12.4,I5)') X(1),X(2),F(I),I
      RETURN
      END
C
C     *******************************
C           BRENT METHOD 2
C     *******************************
C
      SUBROUTINE DPBRENT2(CX)
C
      INCLUDE 'dpcomm.inc'
C
      EXTERNAL DPFUNC2
C
      PARAMETER (NBRM=2,LWA=NBRM*(NBRM+3))
      DIMENSION X(NBRM),F(NBRM),WA(LWA)
C
      EPS  = 1.D2
      M    = 100
      N    = 2
      X(1)=DBLE(CX)
      X(2)=DIMAG(CX)
C
      CALL FBRENTN(DPFUNC2,N,X,F,EPS,INFO,WA,LWA)
      IF(INFO.GE.1.AND.INFO.LE.3) THEN
         X(1)=RKX0
         X(2)=0.D0
      ENDIF
C
      CX=DCMPLX(X(1),X(2))
      RETURN
      END
C
C     *******************************
C           Determinant Function
C     *******************************
C
      SUBROUTINE DPFUNC2(N,X,F,I)
C
      INCLUDE 'dpcomm.inc'
C
      DIMENSION X(N),F(N)
C
      CKX=DCMPLX(X(1),X(2))
      CFUNC2=CFDISP(CRF0,CKX,CKY0,CKZ0,XPOS0,YPOS0,ZPOS0)
      IF(I.EQ.1)THEN
         F(1)=DBLE(CFUNC2)
      ELSE
         F(2)=DIMAG(CFUNC2)
      END IF
      IF(ILIST.NE.0) WRITE(6,'(1P3E12.4,I5)') X(1),X(2),F(I),I
      RETURN
      END
C
C     *******************************
C           BRENT METHOD 3
C     *******************************
C
      SUBROUTINE DPBRENT3(CX)
C
      INCLUDE 'dpcomm.inc'
C
      EXTERNAL DPFUNC3
C
      PARAMETER (NBRM=2,LWA=NBRM*(NBRM+3))
      DIMENSION X(NBRM),F(NBRM),WA(LWA)
C
      EPS  = 1.D2
      M    = 100
      N    = 2
      X(1)=DBLE(SQRT(CX))
      X(2)=DIMAG(SQRT(CX))
C
      CALL FBRENTN(DPFUNC3,N,X,F,EPS,INFO,WA,LWA)
      IF(INFO.GE.1.AND.INFO.LE.3) THEN
         X(1)=RKX0**2
         X(2)=0.D0
      ENDIF
C
      CX=DCMPLX(X(1),X(2))**2
      RETURN
      END
C
C     *******************************
C           Determinant Function
C     *******************************
C
      SUBROUTINE DPFUNC3(N,X,F,I)
C
      INCLUDE 'dpcomm.inc'
C
      DIMENSION X(N),F(N)
C
      CKX=DCMPLX(X(1),X(2))
      CFUNC3=CFDISP(CRF0,CKX,CKY0,CKZ0,XPOS0,YPOS0,ZPOS0)
      IF(I.EQ.1)THEN
         F(1)=DBLE(CFUNC3)
      ELSE
         F(2)=DIMAG(CFUNC3)
      END IF
      IF(ILIST.NE.0) WRITE(6,'(1P3E12.4,I5)') X(1),X(2),F(I),I
      RETURN
      END
C
C     *******************************
C           1D plot
C     *******************************
C
      SUBROUTINE DPGR1A(GX,GY,GY2,NXMAX,ISW)
C
      DIMENSION GX(NXMAX),GY(NXMAX),GY2(NXMAX)
C
      CALL PAGES
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.3,0.)
      CALL GMNMX1(GX,1,NXMAX,1,GXMIN,GXMAX)
      CALL GMNMX1(GY,1,NXMAX,1,GYMIN,GYMAX)
      CALL GQSCAL(GXMIN,GXMAX,GXSMN,GXSMX,GSCALX)
      CALL SETLIN(0,0,5)
      CALL GQSCAL(GYMIN,GYMAX,GYSMN,GYSMX,GSCALY)
      IF(GXMIN*GXMAX.LT.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GXSMN
      ENDIF
      IF(GYMIN*GYMAX.LT.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=GYSMN
      ENDIF
      CALL GDEFIN(4.,17.,2.,17.,GXSMN,GXSMX,GYSMN,GYSMX)
      CALL SETLIN(0,0,7)
      CALL GFRAME
      CALL GSCALE(GXORG,GSCALX,0.0,0.0,0.2,9)
      CALL SETLIN(0,0,5)
      CALL GSCALE(0.0,0.0,GYORG,GSCALY,0.2,1)
      CALL GSCALE(0.0,0.0,0.0,10*(GYMAX-GYMIN),0.0,0)
      CALL SETLIN(0,0,7)
      CALL GVALUE(GXORG,2*GSCALX,0.0,0.0,2)
      CALL SETLIN(0,0,5)
      CALL GVALUE(0.0,0.0,GYORG,2*GSCALY,-2)
      CALL GPLOTP(GX,GY,1,NXMAX,1,0,0,0)
C
      CALL SETLIN(0,0,2)
      CALL GMNMX1(GY2,1,NXMAX,1,GYMIN,GYMAX)
      IF(GYMAX-GYMIN.GT.1.D-15) THEN
         CALL GQSCAL(GYMIN,GYMAX,GYSMN,GYSMX,GSCALY)
         IF(GYMIN*GYMAX.LT.0.0) THEN
            GYORG=0.0
         ELSE
            GYORG=GYSMN
         ENDIF
         CALL GDEFIN(4.,17.,2.,17.,GXSMN,GXSMX,GYSMN,GYSMX)
         CALL GSCALE(0.0,0.0,GYORG,GSCALY,0.2,5)
         CALL GVALUE(0.0,0.0,GYORG,2*GSCALY,-402)
         CALL GVALUE(0.0,0.0,GYORG,2*GSCALY,-402)
         CALL GPLOTP(GX,GY2,1,NXMAX,1,0,0,2)
      ENDIF
C
      CALL DPGPRM(ISW)
      CALL PAGEE
      RETURN
      END
C
C     *******************************
C           1D plot
C     *******************************
C
      SUBROUTINE DPGR1B(GXA,GYA,GY2A,NXM,NXMAXA,IGMAX,ISW)
C
      DIMENSION GXA(NXM,IGMAX),GYA(NXM,IGMAX),GY2A(NXM,IGMAX)
      DIMENSION NXMAXA(IGMAX)
C
      CALL PAGES
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.3,0.)
      GXMIN  = 1.E20
      GXMAX  =-1.E20
      GYMIN  = 1.E20
      GYMAX  =-1.E20
      GYMIN3 = 1.E20
      GYMAX3 =-1.E20
      DO 100 IG=1,IGMAX
         CALL GMNMX1(GXA(1,IG),1,NXMAXA(IG),1,GXMIN1,GXMAX1)
         CALL GMNMX1(GYA(1,IG),1,NXMAXA(IG),1,GYMIN1,GYMAX1)
         CALL GMNMX1(GY2A(1,IG),1,NXMAXA(IG),1,GYMIN2,GYMAX2)
         GXMIN  = MIN(GXMIN,GXMIN1)
         GXMAX  = MAX(GXMAX,GXMAX1)
         GYMIN  = MIN(GYMIN,GYMIN1)
         GYMAX  = MAX(GYMAX,GYMAX1)
         GYMIN3 = MIN(GYMIN3,GYMIN2)
         GYMAX3 = MAX(GYMAX3,GYMAX2)
C         WRITE(6,*) I,GXMIN,GXMAX,GYMIN,GYMAX  
  100 CONTINUE
         GYMIN  = MIN(GYMIN,GYMIN3)
         GYMAX  = MAX(GYMAX,GYMAX3)
      CALL GQSCAL(GXMIN,GXMAX,GXSMN,GXSMX,GSCALX)
      CALL GQSCAL(GYMIN,GYMAX,GYSMN,GYSMX,GSCALY)
      IF(GXMIN*GXMAX.LT.0.0) THEN
         GXORG=0.0
      ELSE
         GXORG=GXSMN
      ENDIF
      IF(GYMIN*GYMAX.LT.0.0) THEN
         GYORG=0.0
      ELSE
         GYORG=GYSMN
      ENDIF
      CALL GDEFIN(4.,17.,2.,17.,GXSMN,GXSMX,GYSMN,GYSMX)
      CALL GFRAME
      CALL GSCALE(GXORG,GSCALX,0.0,0.0,0.2,9)
      CALL GSCALE(0.0,0.0,GYORG,GSCALY,0.2,9)
      CALL GSCALE(0.0,0.0,0.0,10*(GYMAX-GYMIN),0.0,0)
      CALL GVALUE(GXORG,2*GSCALX,0.0,0.0,2)
      CALL GVALUE(0.0,0.0,GYORG,2*GSCALY,-2)
      DO 200 IG=1,IGMAX
      CALL SETLIN(0,0,4)
      	 CALL GPLOTP(GXA(1,IG),GYA(1,IG),1,NXMAXA(IG),1,0,0,0)
      CALL SETLIN(0,0,6)
	 CALL GPLOTP(GXA(1,IG),GY2A(1,IG),1,NXMAXA(IG),1,0,0,1)
  200 CONTINUE
C
C      GYMIN= 1.E20
C      GYMAX=-1.E20
C      DO 300 IG=1,IGMAX
C         CALL GMNMX1(GY2A(1,IG),1,NXMAXA(IG),1,GYMIN1,GYMAX1)
C         GYMIN=MIN(GYMIN,GYMIN1)
C         GYMAX=MAX(GYMAX,GYMAX1)
C  300 CONTINUE
C      IF(GYMAX-GYMIN.GT.1.D-15) THEN
C         CALL GQSCAL(GYMIN,GYMAX,GYSMN,GYSMX,GSCALY)
C         IF(GYMIN*GYMAX.LT.0.0) THEN
C            GYORG=0.0
C         ELSE
C            GYORG=GYSMN
C         ENDIF
C         CALL GDEFIN(4.,17.,2.,17.,GXSMN,GXSMX,GYSMN,GYSMX)
C         CALL GSCALE(0.0,0.0,GYORG,GSCALY,0.2,5)
C         CALL GVALUE(0.0,0.0,GYORG,2*GSCALY,-402)
C         DO 400 IG=1,IGMAX
C      	    CALL GPLOTP(GXA(1,IG),GY2A(1,IG),1,NXMAXA(IG),1,0,0,2)
C  400    CONTINUE
C      ENDIF
C
      CALL DPGPRM(ISW)
      CALL PAGEE
      RETURN
      END
C
C     *******************************
C           DRAW PARAMETERS
C     *******************************
C
      SUBROUTINE DPGPRM(ISW)
C
      INCLUDE 'dpcomm.inc'
C
      CALL SETLIN(0,0,7)
      IF (ISW.EQ.1) THEN
        CALL MOVE(10.0,0.5)
         CALL TEXT (' RF ',4)
        CALL MOVE(0.5,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' RFI ',4)
         CALL SETCHS(0.3,0.0)
      ELSEIF (ISW.EQ.2) THEN
        CALL MOVE(10.0,0.5)
         CALL TEXT (' KX ',4) 
        CALL MOVE(0.5,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' RF ',4) 
         CALL SETCHS(0.3,0.0)
      ELSEIF (ISW.EQ.3) THEN
        CALL MOVE(10.0,0.5)
         CALL TEXT (' KY ',4) 
        CALL MOVE(0.5,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' RF ',4) 
         CALL SETCHS(0.3,0.0)
      ELSEIF (ISW.EQ.4) THEN
        CALL MOVE(10.0,0.5)
         CALL TEXT (' KZ ',4) 
        CALL MOVE(0.5,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' RF ',4) 
         CALL SETCHS(0.3,0.0)
      ELSEIF (ISW.EQ.5) THEN
        CALL MOVE(10.0,0.5)
         CALL TEXT (' X  ',4) 
        CALL MOVE(0.5,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' RF ',4) 
         CALL SETCHS(0.3,0.0)
      ELSEIF (ISW.EQ.6) THEN
        CALL MOVE(10.0,0.5)
         CALL TEXT (' KX ',4) 
        CALL MOVE(0.5,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' KXI ',4) 
         CALL SETCHS(0.3,0.0)
      ELSEIF (ISW.EQ.7) THEN
        CALL MOVE(10.0,0.5)
         CALL TEXT (' X  ',4) 
        CALL MOVE(0.5,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' KX*KX ',7) 
         CALL SETCHS(0.3,0.0)
      ELSEIF (ISW.EQ.8) THEN
        CALL MOVE(10.0,0.5)
         CALL TEXT (' X  ',4) 
        CALL MOVE(0.5,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' KX*KX ',7) 
         CALL SETCHS(0.3,0.0)
      ENDIF
C
      CALL MOVE(20.5,16.5)
      CALL TEXT('RF  =',5)
      CALL NUMBD(RF0,  '(1PE11.3)',11)
      CALL MOVE(20.5,16.0)
      CALL TEXT('RFI0=',5)
      CALL NUMBD(RFI0, '(1PE11.3)',11)
      CALL MOVE(20.5,15.5)
      CALL TEXT('KX  =',5)
      CALL NUMBD(RKX0,'(1PE11.3)',11)
      CALL MOVE(20.5,15.0)
      CALL TEXT('KY  =',5)
      CALL NUMBD(RKY0,'(1PE11.3)',11)
      CALL MOVE(20.5,14.5)
      CALL TEXT('KZ  =',5)
      CALL NUMBD(RKZ0,'(1PE11.3)',11)
      CALL MOVE(20.5,14.0)
      CALL TEXT('BB  =',5)
      CALL NUMBD(BB,'(1PE11.3)',11)
      CALL MOVE(20.5,13.5)
      CALL TEXT('RR  =',5)
      CALL NUMBD(RR,'(1PE11.3)',11)
      CALL MOVE(20.5,13.0)
      CALL TEXT('RA  =',5)
      CALL NUMBD(RA,'(1PE11.3)',11)
      CALL MOVE(20.5,12.5)
      CALL TEXT('PA1 =',5)
      CALL NUMBD(PA(1),'(1PE11.3)',11)
      CALL MOVE(20.5,12.0)
      CALL TEXT('PA2 =',5)
      CALL NUMBD(PA(2),'(1PE11.3)',11)
      CALL MOVE(20.5,11.5)
      CALL TEXT('PA3 =',5)
      CALL NUMBD(PA(3),'(1PE11.3)',11)
      CALL MOVE(20.5,11.0)
      CALL TEXT('PZ1 =',5)
      CALL NUMBD(PZ(1),'(1PE11.3)',11)
      CALL MOVE(20.5,10.5)
      CALL TEXT('PZ2 =',5)
      CALL NUMBD(PZ(2),'(1PE11.3)',11)
      CALL MOVE(20.5,10.0)
      CALL TEXT('PZ3 =',5)
      CALL NUMBD(PZ(3),'(1PE11.3)',11)
      CALL MOVE(20.5,9.5)
      CALL TEXT('PN1 =',5)
      CALL NUMBD(PN(1),'(1PE11.3)',11)
      CALL MOVE(20.5,9.0)
      CALL TEXT('PN2 =',5)
      CALL NUMBD(PN(2),'(1PE11.3)',11)
      CALL MOVE(20.5,8.5)
      CALL TEXT('PN3 =',5)
      CALL NUMBD(PN(3),'(1PE11.3)',11)
      CALL MOVE(20.5,8.0)
      CALL TEXT('PNS1=',5)
      CALL NUMBD(PNS(1),'(1PE11.3)',11)
      CALL MOVE(20.5,7.5)
      CALL TEXT('PNS2=',5)
      CALL NUMBD(PNS(2),'(1PE11.3)',11)
      CALL MOVE(20.5,7.0)
      CALL TEXT('PNS3=',5)
      CALL NUMBD(PNS(3),'(1PE11.3)',11)
      CALL MOVE(20.5,6.5)
      CALL TEXT('PTPR(e)=',8)
      CALL NUMBD(PTPR(1),'(1PE9.3)',9)
      CALL MOVE(20.5,6.0)
      CALL TEXT('PTPP(e)=',8)
      CALL NUMBD(PTPP(1),'(1PE9.3)',9)
      CALL MOVE(20.5,5.5)
      CALL TEXT('PTPR(i)=',8)
      CALL NUMBD(PTPR(2),'(1PE9.3)',9)
      CALL MOVE(20.5,5.0)
      CALL TEXT('PTPP(i)=',8)
      CALL NUMBD(PTPP(2),'(1PE9.3)',9)
      CALL MOVE(20.5,4.5)
      CALL TEXT('PTS    =',8)
      CALL NUMBD(PTS(1),'(1PE9.3)',9)
      CALL MOVE(20.5,3.5)
      CALL TEXT('MODELP =',8)
      CALL NUMBI(MODELP(1),'(I1)',1)
C
      RETURN
      END
