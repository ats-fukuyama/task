C     $Id$
C
C     ***** PLOT CONTOUR GRAPH *************************************
C
      SUBROUTINE DPCONT
C
      USE plcomm
      USE pllocal
      USE plparm,ONLY: pl_parm,pl_view
      INCLUDE 'dpcomm.inc'
      DIMENSION GX(NGXM),GY(NGYM),GZ(NGXM,NGYM)
      DIMENSION Z(NGXM,NGYM)
      DIMENSION KA(8,NGXM,NGYM)
      CHARACTER KID*2,KID1*1,KID2*1
C
      DATA INIT/0/
C
      IF(INIT.EQ.0) THEN
         KID='YW'
         NXMAX=201
         NYMAX=201
         INIT=1
      ENDIF
C
    1 WRITE(6,*)'## SELECT : XY : W/RF X/KX Y/KY Z/KZ : ',
     &          'P,D,V/PARM 0/EXIT'
      READ(5,'(A2)',ERR=1,END=9000) KID
      KID1=KID(1:1)
      KID2=KID(2:2)
      CALL GUCPTL(KID1)
      CALL GUCPTL(KID2)
      IF(KID1.EQ.'0') THEN
         GOTO 9000
      ELSEIF(KID1.EQ.'P') THEN
         CALL PL_PARM(0,'PL',IERR)
         GOTO 1
      ELSEIF(KID1.EQ.'D') THEN
         CALL DPPARM(0,'DP',IERR)
         GOTO 1
      ELSEIF(KID1.EQ.'V') THEN
         CALL PL_VIEW
         CALL DPVIEW
         GOTO 1
      ENDIF
C
    2 IF(KID1.EQ.'W') THEN
         WRITE(6,'(A,1P2E12.4,I5)') 
     &        ' INPUT X-AXIS: RF1,RF2,NXMAX',RF1,RF2,NXMAX
         XMIN=RF1
         XMAX=RF2
      ELSE IF(KID1.EQ.'X') THEN
         WRITE(6,'(A,1P2E12.4,I5)') 
     &        ' INPUT X-AXIS: RKX1,RKX2,NXMAX',RKX1,RKX2,NXMAX
         XMIN=RKX1
         XMAX=RKX2
      ELSE IF(KID1.EQ.'Y') THEN
         WRITE(6,'(A,1P2E12.4,I5)') 
     &        ' INPUT X-AXIS: RKY1,RKY2,NXMAX',RKY1,RKY2,NXMAX
         XMIN=RKY1
         XMAX=RKY2
      ELSE IF(KID1.EQ.'Z') THEN
         WRITE(6,'(A,1P2E12.4,I5)') 
     &        ' INPUT X-AXIS: RKZ1,RKZ2,NXMAX',RKZ1,RKZ2,NXMAX
         XMIN=RKZ1
         XMAX=RKZ2
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID1'
         GOTO 1
      ENDIF
      NXMAX_SAVE=NXMAX
      READ(5,*,ERR=2,END=101) XMIN,XMAX,NXMAX
      IF(NXMAX.EQ.0) GO TO 101
      IF(NXMAX.GT.NGXM) THEN
         WRITE(6,*) 'XX NXMAX EXCEEDS NGXM =',NGXM
         GOTO 2
      ENDIF
      IF(KID1.EQ.'W') THEN
         RF1=XMIN
         RF2=XMAX
      ELSE IF(KID1.EQ.'X') THEN
         RKX1=XMIN
         RKX2=XMAX
      ELSE IF(KID1.EQ.'Y') THEN
         RKY1=XMIN
         RKY2=XMAX
      ELSE IF(KID1.EQ.'Z') THEN
         RKZ1=XMIN
         RKZ2=XMAX
      ENDIF
C
    3 IF(KID2.EQ.'W') THEN
         WRITE(6,'(A,1P2E12.4,I5)') 
     &        ' INPUT Y-AXIS: RF1,RF2,NYMAX',RF1,RF2,NYMAX
         YMIN=RF1
         YMAX=RF2
      ELSE IF(KID2.EQ.'X') THEN
         WRITE(6,'(A,1P2E12.4,I5)') 
     &        ' INPUT Y-AXIS: RKX1,RKX2,NYMAX',RKX1,RKX2,NYMAX
         YMIN=RKX1
         YMAX=RKX2
      ELSE IF(KID2.EQ.'Y') THEN
         WRITE(6,'(A,1P2E12.4,I5)') 
     &        ' INPUT Y-AXIS: RKY1,RKY2,NYMAX',RKY1,RKY2,NYMAX
         YMIN=RKY1
         YMAX=RKY2
      ELSE IF(KID2.EQ.'Z') THEN
         WRITE(6,'(A,1P2E12.4,I5)') 
     &        ' INPUT Y-AXIS: RKZ1,RKZ2,NYMAX',RKZ1,RKZ2,NYMAX
         YMIN=RKZ1
         YMAX=RKZ2
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID2'
         GOTO 1
      ENDIF
      NYMAX_SAVE=NYMAX
      READ(5,*,ERR=3,END=102) YMIN,YMAX,NYMAX
      IF(NYMAX.EQ.0) GO TO 102
      IF(NYMAX.GT.NGYM) THEN
         WRITE(6,*) 'XX NYMAX EXCEEDS NGYM =',NGYM
         GOTO 3
      ENDIF
      IF(KID2.EQ.'W') THEN
         RF1=YMIN
         RF2=YMAX
      ELSE IF(KID2.EQ.'X') THEN
         RKX1=YMIN
         RKX2=YMAX
      ELSE IF(KID2.EQ.'Y') THEN
         RKY1=YMIN
         RKY2=YMAX
      ELSE IF(KID2.EQ.'Z') THEN
         RKZ1=YMIN
         RKZ2=YMAX
      ENDIF
C
      DY=(YMAX-YMIN)/(NYMAX-1)
      DX=(XMAX-XMIN)/(NXMAX-1)
C
      CRF=DCMPLX(RF0,0.D0)
      CKX=RKX0
      CKY=RKY0
      CKZ=RKZ0
C
      DO NY=1,NYMAX
         Y=YMIN+DY*(NY-1)
         IF(KID2.EQ.'W') THEN
            CRF=DCMPLX(Y,0.D0)
         ELSEIF(KID2.EQ.'X') THEN
            CKX=Y
         ELSEIF(KID2.EQ.'Y') THEN
            CKY=Y
         ELSEIF(KID2.EQ.'Z') THEN
            CKZ=Y
         ENDIF
         DO NX=1,NXMAX
            X=XMIN+DX*(NX-1)
            IF(KID1.EQ.'W') THEN
               CRF=DCMPLX(X,0.D0)
            ELSEIF(KID1.EQ.'X') THEN
               CKX=X
            ELSEIF(KID1.EQ.'Y') THEN
               CKY=X
            ELSEIF(KID1.EQ.'Z') THEN
               CKZ=X
            ENDIF
            CD=CFDISP(CRF,CKX,CKY,CKZ,RX0,RY0,RZ0)

C *** remove sign change due to cyclotron resonance ***
            
            CW=2.D0*PI*1.D6*CRF
            DO NS=1,NSMAX
               CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
               CD=CD*(1.D0-CWC)
            ENDDO
C
            GX(NX)=GUCLIP(X)
            GY(NY)=GUCLIP(Y)
            Z(NX,NY)=DBLE(CD)
         ENDDO
      ENDDO
C
      DO NY=1,NYMAX
         IF(NY.EQ.1) THEN
            NY1=NY
            NY2=NY+1
            FACTY=1.D0
         ELSEIF(NY.EQ.NYMAX) THEN
            NY1=NY-1
            NY2=NY
            FACTY=1.D0
         ELSE
            NY1=NY-1
            NY2=NY+1
            FACTY=0.5D0
         ENDIF
      DO NX=1,NXMAX
         IF(NX.EQ.1) THEN
            NX1=NX
            NX2=NX+1
            FACTX=1.D0
         ELSEIF(NX.EQ.NXMAX) THEN
            NX1=NX-1
            NX2=NX
            FACTX=1.D0
         ELSE
            NX1=NX-1
            NX2=NX+1
            FACTX=0.5D0
         ENDIF
         DZY=FACTY*(Z(NX,NY2)-Z(NX,NY1))
         DZX=FACTX*(Z(NX2,NY)-Z(NX1,NY))
         DZ2=MAX(SQRT(DZX**2+DZY**2),1.D0)
         GZ(NX,NY)=GUCLIP(Z(NX,NY)/DZ2)
      ENDDO
      ENDDO
C
C      DO NY=1,NYMAX
C      DO NX=1,NXMAX
C            IF(MOD(NX-1,10).EQ.0.AND.MOD(NY-1,10).EQ.0) THEN
C               WRITE(6,'(2I5,1P4E15.7)') 
C     &              NX,NY,GX(NX),GY(NY),Z(NX,NY),GZ(NX,NY)
C            ENDIF
C      ENDDO
C      ENDDO
C
      CALL PAGES
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.3,0.)
      CALL GMNMX1(GX,1,NXMAX,1,GXMIN,GXMAX)
      CALL GMNMX1(GY,1,NYMAX,1,GYMIN,GYMAX)
      CALL GQSCAL(GXMIN,GXMAX,GXSMN,GXSMX,GSCALX)
      CALL GQSCAL(GYMIN,GYMAX,GYSMN,GYSMX,GSCALY)
      CALL GDEFIN(3.,18.,2.,17.,GXSMN,GXSMX,GYSMN,GYSMX)
      CALL GFRAME
      CALL GSCALE(GXSMN,GSCALX,0.0,0.0,0.2,9)
      CALL GSCALE(0.0,0.0,GYSMN,GSCALY,0.2,9)
      CALL GVALUE(GXSMN,2*GSCALX,0.0,0.0,NGULEN(2*GSCALX))
      CALL GVALUE(0.0,0.0,GYSMN,2*GSCALY,NGULEN(2*GSCALY))
      CALL CONTQ2(GZ,GX,GY,NGXM,NXMAX,NYMAX,0.0,2.0,1,0,0,KA)
C
      CALL MOVE(3.0,17.5)
      IF(KID1.EQ.'W') THEN
         CALL TEXT ('RF',2)
      ELSEIF(KID1.EQ.'X') THEN
         CALL TEXT ('KX',2)
      ELSEIF(KID1.EQ.'Y') THEN
         CALL TEXT ('KY',2)
      ELSEIF(KID1.EQ.'Z') THEN
         CALL TEXT ('KZ',2)
      ENDIF
      CALL TEXT('-',1)
      IF(KID2.EQ.'W') THEN
         CALL TEXT ('RF',2)
      ELSEIF(KID2.EQ.'X') THEN
         CALL TEXT ('KX',2)
      ELSEIF(KID2.EQ.'Y') THEN
         CALL TEXT ('KY',2)
      ELSEIF(KID2.EQ.'Z') THEN
         CALL TEXT ('KZ',2)
      ENDIF
C
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
      GOTO 2
C
 9000 RETURN
  101 NXMAX=NXMAX_SAVE
      GO TO 1
  102 NYMAX=NYMAX_SAVE
      GO TO 2
      END
C
C     ***** PLOT CONTOUR GRAPH *****
C     ***** DENSITY SCAN *****
C
      SUBROUTINE DPCONTX
C
      USE plcomm
      USE pllocal
      INCLUDE 'dpcomm.inc'
      PARAMETER (NGPM=21)
      DIMENSION GX(NGXM),GY(NGYM),GZ(NGXM,NGYM,NGPM)
      DIMENSION Z(NGXM,NGYM,NGPM)
      DIMENSION PNSAVE(NSM),PNNGP(NGPM)
      DIMENSION KA(8,NGXM,NGYM)
      CHARACTER KID*2,KID1*1,KID2*1
C
      DATA INIT/0/
C
      IF(INIT.EQ.0) THEN
         KID='YW'
         XMIN=1.D0
         XMAX=100.D0
         NXMAX=51
         YMIN=1.D0
         YMAX=100.D0
         NYMAX=51
         INIT=1
      ENDIF
C
    1 WRITE(6,*)'## SELECT : XY : W/RF X/KX Y/KY Z/KZ : P,V/PARM Q/END'
      READ(5,'(A2)',ERR=1,END=9000) KID
      KID1=KID(1:1)
      KID2=KID(2:2)
      CALL GUCPTL(KID1)
      CALL GUCPTL(KID2)
      IF(KID1.EQ.'Q') THEN
         GOTO 9000
      ELSEIF(KID1.EQ.'P') THEN
         CALL DPPARM(0,'DP',IERR)
         GOTO 1
      ELSEIF(KID1.EQ.'V') THEN
         CALL DPVIEW
         GOTO 1
      ENDIF
C
    2 IF(KID1.EQ.'W') THEN
         WRITE(6,*) ' INPUT X-AXIS: RF1,RF2,NXMAX'
      ELSE IF(KID1.EQ.'X') THEN
         WRITE(6,*) ' INPUT X-AXIS: RKX1,RKX2,NXMAX'
      ELSE IF(KID1.EQ.'Y') THEN
         WRITE(6,*) ' INPUT X-AXIS: RKY1,RKY2,NXMAX'
      ELSE IF(KID1.EQ.'Z') THEN
         WRITE(6,*) ' INPUT X-AXIS: RKZ1,RKZ2,NXMAX'
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID1'
         GOTO 1
      ENDIF
      READ(5,*,ERR=2,END=1) XMIN,XMAX,NXMAX
      IF(NXMAX.GT.NGXM) THEN
         WRITE(6,*) 'XX NXMAX EXCEEDS NGXM =',NGXM
         GOTO 2
      ENDIF
C
    3 IF(KID2.EQ.'W') THEN
         WRITE(6,*) ' INPUT Y-AXIS: RF1,RF2,NYMAX'
      ELSE IF(KID2.EQ.'X') THEN
         WRITE(6,*) ' INPUT Y-AXIS: RKX1,RKX2,NYMAX'
      ELSE IF(KID2.EQ.'Y') THEN
         WRITE(6,*) ' INPUT Y-AXIS: RKY1,RKY2,NYMAX'
      ELSE IF(KID2.EQ.'Z') THEN
         WRITE(6,*) ' INPUT Y-AXIS: RKZ1,RKZ2,NYMAX'
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID2'
         GOTO 1
      ENDIF
      READ(5,*,ERR=3,END=2) YMIN,YMAX,NYMAX
      IF(NYMAX.GT.NGYM) THEN
         WRITE(6,*) 'XX NYMAX EXCEEDS NGYM =',NGYM
         GOTO 3
      ENDIF
C
      DY=(YMAX-YMIN)/(NYMAX-1)
      DX=(XMAX-XMIN)/(NXMAX-1)
C
      CRF=DCMPLX(RF0,0.D0)
      CKX=RKX0
      CKY=RKY0
      CKZ=RKZ0
C
    4 WRITE(6,*) ' INPUT : NUMBER OF DENSITY ?'
      READ(5,*,ERR=4,END=3) NGPMAX
      IF(NGPMAX.LT.1.OR.NGPMAx.GT.NGPM) THEN
         WRITE(6,*) 'XX NGPMAX IS ZERO OR EXCEEDS NGPM =',NGPM
         GOTO 4
      ENDIF
C
      DO NS=1,NSMAX
         PNSAVE(NS)=PN(NS)
      ENDDO
C
      DO NGP=1,NGPMAX
    5    WRITE(6,*) ' INPUT : DENSITY FOR NGP =',NGP
         READ(5,*,ERR=5,END=4) PNNGP(NGP)
         PN(1)=PNNGP(NGP)
         PN(2)=PNNGP(NGP)/PZ(2)
C
      DO NY=1,NYMAX
         Y=YMIN+DY*(NY-1)
         IF(KID2.EQ.'W') THEN
            CRF=DCMPLX(Y,0.D0)
         ELSEIF(KID2.EQ.'X') THEN
            CKX=Y
         ELSEIF(KID2.EQ.'Y') THEN
            CKY=Y
         ELSEIF(KID2.EQ.'Z') THEN
            CKZ=Y
         ENDIF
         DO NX=1,NXMAX
            X=XMIN+DX*(NX-1)
            IF(KID1.EQ.'W') THEN
               CRF=DCMPLX(X,0.D0)
            ELSEIF(KID1.EQ.'X') THEN
               CKX=X
            ELSEIF(KID1.EQ.'Y') THEN
               CKY=X
            ELSEIF(KID1.EQ.'Z') THEN
               CKZ=X
            ENDIF
            CD=CFDISP(CRF,CKX,CKY,CKZ,RX0,RY0,RZ0)
            CW=2.D0*PI*1.D6*CRF
            DO NS=1,NSMAX
               CWC=BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
               CD=CD*(1.D0-CWC)
            ENDDO
C
            GX(NX)=GUCLIP(X)
            GY(NY)=GUCLIP(Y)
            Z(NX,NY,NGP)=DBLE(CD)
         ENDDO
      ENDDO
      ENDDO
C
      DO NGP=1,NGPMAX
      DO NY=1,NYMAX
         IF(NY.EQ.1) THEN
            NY1=NY
            NY2=NY+1
            FACTY=1.D0
         ELSEIF(NY.EQ.NYMAX) THEN
            NY1=NY-1
            NY2=NY
            FACTY=1.D0
         ELSE
            NY1=NY-1
            NY2=NY+1
            FACTY=0.5D0
         ENDIF
      DO NX=1,NXMAX
         IF(NX.EQ.1) THEN
            NX1=NX
            NX2=NX+1
            FACTX=1.D0
         ELSEIF(NX.EQ.NXMAX) THEN
            NX1=NX-1
            NX2=NX
            FACTX=1.D0
         ELSE
            NX1=NX-1
            NX2=NX+1
            FACTX=0.5D0
         ENDIF
         DZY=FACTY*(Z(NX,NY2,NGP)-Z(NX,NY1,NGP))
         DZX=FACTX*(Z(NX2,NY,NGP)-Z(NX1,NY,NGP))
         DZ2=MAX(SQRT(DZX**2+DZY**2),1.D0)
         GZ(NX,NY,NGP)=GUCLIP(Z(NX,NY,NGP)/DZ2)
      ENDDO
      ENDDO
      ENDDO
C
      DO NS=1,NSMAX
         PN(NS)=PNSAVE(NS)
      ENDDO
C
      CALL PAGES
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.3,0.)
      CALL GMNMX1(GX,1,NXMAX,1,GXMIN,GXMAX)
      CALL GMNMX1(GY,1,NYMAX,1,GYMIN,GYMAX)
      CALL GQSCAL(GXMIN,GXMAX,GXSMN,GXSMX,GSCALX)
      CALL GQSCAL(GYMIN,GYMAX,GYSMN,GYSMX,GSCALY)
      CALL GDEFIN(3.,18.,2.,17.,GXSMN,GXSMX,GYSMN,GYSMX)
      CALL GFRAME
      CALL GSCALE(GXSMN,GSCALX,0.0,0.0,0.2,9)
      CALL GSCALE(0.0,0.0,GYSMN,GSCALY,0.2,9)
      CALL GVALUE(GXSMN,2*GSCALX,0.0,0.0,NGULEN(2*GSCALX))
      CALL GVALUE(0.0,0.0,GYSMN,2*GSCALY,NGULEN(2*GSCALY))
      DO NGP=1,NGPMAX
         CALL SETLIN(0,0,7-MOD(NGP-1,5))
         CALL CONTQ2(GZ(1,1,NGP),GX,GY,NGXM,NXMAX,NYMAX,
     &               0.0,2.0,1,0,0,KA)
      ENDDO
      CALL SETLIN(0,0,7)
C
      CALL MOVE(3.0,17.5)
      IF(KID1.EQ.'W') THEN
         CALL TEXT ('RF',2)
      ELSEIF(KID1.EQ.'X') THEN
         CALL TEXT ('KX',2)
      ELSEIF(KID1.EQ.'Y') THEN
         CALL TEXT ('KY',2)
      ELSEIF(KID1.EQ.'Z') THEN
         CALL TEXT ('KZ',2)
      ENDIF
      CALL TEXT('-',1)
      IF(KID2.EQ.'W') THEN
         CALL TEXT ('RF',2)
      ELSEIF(KID2.EQ.'X') THEN
         CALL TEXT ('KX',2)
      ELSEIF(KID2.EQ.'Y') THEN
         CALL TEXT ('KY',2)
      ELSEIF(KID2.EQ.'Z') THEN
         CALL TEXT ('KZ',2)
      ENDIF
C
      CALL MOVE(20.0,16.5)
      CALL TEXT('BB  =',5)
      CALL NUMBD(BB,'(1PE11.3)',11)
      CALL MOVE(20.0,16.0)
      CALL TEXT('RF  =',5)
      CALL NUMBD(RF0,  '(1PE11.3)',11)
      CALL MOVE(20.0,15.5)
      CALL TEXT('RFI0=',5)
      CALL NUMBD(RFI0, '(1PE11.3)',11)
      CALL MOVE(20.0,15.0)
      CALL TEXT('KX  =',5)
      CALL NUMBD(RKX0,'(1PE11.3)',11)
      CALL MOVE(20.0,14.5)
      CALL TEXT('KY  =',5)
      CALL NUMBD(RKY0,'(1PE11.3)',11)
      CALL MOVE(20.0,14.0)
      CALL TEXT('KZ  =',5)
      CALL NUMBD(RKZ0,'(1PE11.3)',11)
      CALL MOVE(20.0,13.5)
      CALL TEXT('BB  =',5)
      CALL NUMBD(BB,'(1PE11.3)',11)
C
      DO NGP=1,NGPMAX
         CALL MOVE(20.0,12.5-NGP*0.5)
         CALL SETLIN(0,0,7-MOD(NGP-1,5))
         CALL TEXT('PN  =',5)
         CALL NUMBD(PNNGP(NGP),'(1PE11.3)',11)
      ENDDO
      CALL SETLIN(0,0,7)
      CALL PAGEE
      GOTO 1
C
 9000 RETURN
      END
