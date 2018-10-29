MODULE dpcont

CONTAINS

!     ***** PLOT CONTOUR GRAPH *************************************

  SUBROUTINE DP_CONT

    USE dpcomm
    USE pllocal
    USE PLPARM,ONLY: pl_parm,pl_view
    USE dpparm
    USE dpdisp
    IMPLICIT NONE
    REAL(4),ALLOCATABLE:: GX(:),GY(:),GZ(:,:)
    REAL(rkind),ALLOCATABLE:: Z(:,:)
    INTEGER,ALLOCATABLE:: KA(:,:,:)
    CHARACTER(LEN=2):: KID
    CHARACTER(LEN=1):: KID1,KID2
    INTEGER,SAVE:: INIT=0
    INTEGER:: NX,NY,NS,NY1,NY2,NX1,NX2,NGULEN,IERR
    REAL(rkind):: XMIN,XMAX,YMIN,YMAX,DX,DY,Y,X,FACTY,FACTX,DZY,DZX,DZ2
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CD,CW,CWC
    REAL(4):: GUCLIP
    REAL(4):: GXMIN,GXMAX,GYMIN,GYMAX,GXSMN,GXSMX,GSCALX,GYSMN,GYSMX,GSCALY

      IF(INIT.EQ.0) THEN
         KID='YW'
         XMIN=0.D0
         XMAX=100.D0
         NGXMAX=51
         YMIN=0.D0
         YMAX=100.D0
         NGYMAX=51
         INIT=1
      ENDIF

    1 WRITE(6,*)'## SELECT : XY : W/RF X/KX Y/KY Z/KZ : ', &
                'P,D,V/PARM X/EXIT'
      READ(5,'(A2)',ERR=1,END=9000) KID
      KID1=KID(1:1)
      KID2=KID(2:2)
      CALL GUCPTL(KID1)
      CALL GUCPTL(KID2)
      IF(KID1.EQ.'X'.AND.KID2.EQ.' ') THEN
         GOTO 9000
      ELSEIF(KID1.EQ.'P') THEN
         CALL PL_PARM(0,'PL',IERR)
         GOTO 1
      ELSEIF(KID1.EQ.'D') THEN
         CALL DP_PARM(0,'DP',IERR)
         GOTO 1
      ELSEIF(KID1.EQ.'V') THEN
         CALL PL_VIEW
         CALL DP_VIEW
         GOTO 1
      ENDIF

    2 IF(KID1.EQ.'W') THEN
         WRITE(6,*) ' INPUT X-AXIS: RF1,RF2,NGXMAX'
      ELSE IF(KID1.EQ.'X') THEN
         WRITE(6,*) ' INPUT X-AXIS: RKX1,RKX2,NGXMAX'
      ELSE IF(KID1.EQ.'Y') THEN
         WRITE(6,*) ' INPUT X-AXIS: RKY1,RKY2,NGXMAX'
      ELSE IF(KID1.EQ.'Z') THEN
         WRITE(6,*) ' INPUT X-AXIS: RKZ1,RKZ2,NGXMAX'
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID1'
         GOTO 1
      ENDIF
      READ(5,*,ERR=2,END=1) XMIN,XMAX,NGXMAX
      IF(NGXMAX.LE.0) GOTO 1

    3 IF(KID2.EQ.'W') THEN
         WRITE(6,*) ' INPUT Y-AXIS: RF1,RF2,NGYMAX'
      ELSE IF(KID2.EQ.'X') THEN
         WRITE(6,*) ' INPUT Y-AXIS: RKX1,RKX2,NGYMAX'
      ELSE IF(KID2.EQ.'Y') THEN
         WRITE(6,*) ' INPUT Y-AXIS: RKY1,RKY2,NGYMAX'
      ELSE IF(KID2.EQ.'Z') THEN
         WRITE(6,*) ' INPUT Y-AXIS: RKZ1,RKZ2,NGYMAX'
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID2'
         GOTO 1
      ENDIF
      READ(5,*,ERR=3,END=2) YMIN,YMAX,NGYMAX
      IF(NGYMAX.LE.0) GOTO 2

      ALLOCATE(GX(NGXMAX),GY(NGYMAX),GZ(NGXMAX,NGYMAX))
      ALLOCATE(Z(NGXMAX,NGYMAX),KA(8,NGXMAX,NGYMAX))

      DY=(YMAX-YMIN)/(NGYMAX-1)
      DX=(XMAX-XMIN)/(NGXMAX-1)

      CRF=DCMPLX(RF0,0.D0)
      CKX=RKX0
      CKY=RKY0
      CKZ=RKZ0

      DO NY=1,NGYMAX
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
         DO NX=1,NGXMAX
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

            GX(NX)=GUCLIP(X)
            GY(NY)=GUCLIP(Y)
            Z(NX,NY)=DBLE(CD)
         ENDDO
      ENDDO

      DO NY=1,NGYMAX
         IF(NY.EQ.1) THEN
            NY1=NY
            NY2=NY+1
            FACTY=1.D0
         ELSEIF(NY.EQ.NGYMAX) THEN
            NY1=NY-1
            NY2=NY
            FACTY=1.D0
         ELSE
            NY1=NY-1
            NY2=NY+1
            FACTY=0.5D0
         ENDIF
      DO NX=1,NGXMAX
         IF(NX.EQ.1) THEN
            NX1=NX
            NX2=NX+1
            FACTX=1.D0
         ELSEIF(NX.EQ.NGXMAX) THEN
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

      DO NY=1,NGYMAX
         DO NX=1,NGXMAX
            IF(MOD(NX-1,10).EQ.0.AND.MOD(NY-1,10).EQ.0) THEN
               WRITE(6,'(2I5,1P4E15.7)')  &
                    NX,NY,GX(NX),GY(NY),Z(NX,NY),GZ(NX,NY)
            ENDIF
         ENDDO
      ENDDO

      CALL PAGES
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.3,0.)
      CALL GMNMX1(GX,1,NGXMAX,1,GXMIN,GXMAX)
      CALL GMNMX1(GY,1,NGYMAX,1,GYMIN,GYMAX)
      CALL GQSCAL(GXMIN,GXMAX,GXSMN,GXSMX,GSCALX)
      CALL GQSCAL(GYMIN,GYMAX,GYSMN,GYSMX,GSCALY)
      CALL GDEFIN(3.,18.,2.,17.,GXSMN,GXSMX,GYSMN,GYSMX)
      CALL GFRAME
      CALL GSCALE(GXSMN,GSCALX,0.0,0.0,0.2,9)
      CALL GSCALE(0.0,0.0,GYSMN,GSCALY,0.2,9)
      CALL GVALUE(GXSMN,2*GSCALX,0.0,0.0,NGULEN(2*GSCALX))
      CALL GVALUE(0.0,0.0,GYSMN,2*GSCALY,NGULEN(2*GSCALY))
      CALL CONTQ2(GZ,GX,GY,NGXMAX,NGXMAX,NGYMAX,0.0,2.0,1,0,0,KA)

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
      DEALLOCATE(GX,GY,GZ,Z,KA)
      GOTO 2

 9000 RETURN
  END SUBROUTINE DP_CONT

!     ***** PLOT CONTOUR GRAPH *****
!     ***** DENSITY SCAN *****

  SUBROUTINE DP_CONTX

    USE dpcomm
    USE pllocal
    USE plinit
    USE dpparm
    USE dpdisp
    IMPLICIT NONE
    REAL(4),ALLOCATABLE:: GX(:),GY(:),GZ(:,:,:)
    REAL(rkind),ALLOCATABLE:: Z(:,:,:),PNSAVE(:),PNNGP(:)
    INTEGER,ALLOCATABLE:: KA(:,:,:)
    CHARACTER(LEN=2):: KID
    CHARACTER(LEN=1):: KID1,KID2
    INTEGER,SAVE:: INIT=0
    INTEGER:: NS,NGP,NY,NX,NY1,NY2,NX1,NX2,IERR,NGULEN
    REAL(rkind):: XMIN,XMAX,YMIN,YMAX,DX,DY,Y,X,FACTY,FACTX,DZY,DZX,DZ2
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CD,CW,CWC
    REAL(4):: GUCLIP
    REAL(4):: GXMIN,GXMAX,GYMIN,GYMAX,GXSMN,GXSMX,GSCALX,GYSMN,GYSMX,GSCALY

      IF(INIT.EQ.0) THEN
         KID='YW'
         XMIN=1.D0
         XMAX=100.D0
         NGXMAX=51
         YMIN=1.D0
         YMAX=100.D0
         NGYMAX=51
         INIT=1
      ENDIF

    1 WRITE(6,*)'## SELECT : XY : W/RF X/KX Y/KY Z/KZ : P,V/PARM Q/END'
      READ(5,'(A2)',ERR=1,END=9000) KID
      KID1=KID(1:1)
      KID2=KID(2:2)
      CALL GUCPTL(KID1)
      CALL GUCPTL(KID2)
      IF(KID1.EQ.'Q') THEN
         GOTO 9000
      ELSEIF(KID1.EQ.'P') THEN
         CALL DP_PARM(0,'DP',IERR)
         GOTO 1
      ELSEIF(KID1.EQ.'V') THEN
         CALL DP_VIEW
         GOTO 1
      ENDIF

    2 IF(KID1.EQ.'W') THEN
         WRITE(6,*) ' INPUT X-AXIS: RF1,RF2,NGXMAX'
      ELSE IF(KID1.EQ.'X') THEN
         WRITE(6,*) ' INPUT X-AXIS: RKX1,RKX2,NGXMAX'
      ELSE IF(KID1.EQ.'Y') THEN
         WRITE(6,*) ' INPUT X-AXIS: RKY1,RKY2,NGXMAX'
      ELSE IF(KID1.EQ.'Z') THEN
         WRITE(6,*) ' INPUT X-AXIS: RKZ1,RKZ2,NGXMAX'
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID1'
         GOTO 1
      ENDIF
      READ(5,*,ERR=2,END=1) XMIN,XMAX,NGXMAX
      IF(NGXMAX.LE.0) GOTO 1

    3 IF(KID2.EQ.'W') THEN
         WRITE(6,*) ' INPUT Y-AXIS: RF1,RF2,NGYMAX'
      ELSE IF(KID2.EQ.'X') THEN
         WRITE(6,*) ' INPUT Y-AXIS: RKX1,RKX2,NGYMAX'
      ELSE IF(KID2.EQ.'Y') THEN
         WRITE(6,*) ' INPUT Y-AXIS: RKY1,RKY2,NGYMAX'
      ELSE IF(KID2.EQ.'Z') THEN
         WRITE(6,*) ' INPUT Y-AXIS: RKZ1,RKZ2,NGYMAX'
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID2'
         GOTO 1
      ENDIF
      READ(5,*,ERR=3,END=2) YMIN,YMAX,NGYMAX
      IF(NGYMAX.LE.0) GOTO 2

      DY=(YMAX-YMIN)/(NGYMAX-1)
      DX=(XMAX-XMIN)/(NGXMAX-1)

      CRF=DCMPLX(RF0,0.D0)
      CKX=RKX0
      CKY=RKY0
      CKZ=RKZ0

    4 WRITE(6,*) ' INPUT : NUMBER OF DENSITY ?'
      READ(5,*,ERR=4,END=3) NGPMAX
      IF(NGPMAX.LE.0) GOTO 3

      ALLOCATE(GX(NGXMAX),GY(NGYMAX),GZ(NGXMAX,NGYMAX,NGPMAX))
      ALLOCATE(Z(NGXMAX,NGYMAX,NGPMAX),KA(8,NGXMAX,NGYMAX))
      ALLOCATE(PNSAVE(NSMAX),PNNGP(NGPMAX))

      DO NS=1,NSMAX
         PNSAVE(NS)=PN(NS)
      ENDDO

      DO NGP=1,NGPMAX
    5    WRITE(6,*) ' INPUT : DENSITY FOR NGP =',NGP
         READ(5,*,ERR=5,END=4) PNNGP(NGP)
         PN(1)=PNNGP(NGP)
         PN(2)=PNNGP(NGP)/PZ(2)

      DO NY=1,NGYMAX
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
         DO NX=1,NGXMAX
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

            GX(NX)=GUCLIP(X)
            GY(NY)=GUCLIP(Y)
            Z(NX,NY,NGP)=DBLE(CD)
         ENDDO
      ENDDO
      ENDDO

      DO NGP=1,NGPMAX
      DO NY=1,NGYMAX
         IF(NY.EQ.1) THEN
            NY1=NY
            NY2=NY+1
            FACTY=1.D0
         ELSEIF(NY.EQ.NGYMAX) THEN
            NY1=NY-1
            NY2=NY
            FACTY=1.D0
         ELSE
            NY1=NY-1
            NY2=NY+1
            FACTY=0.5D0
         ENDIF
      DO NX=1,NGXMAX
         IF(NX.EQ.1) THEN
            NX1=NX
            NX2=NX+1
            FACTX=1.D0
         ELSEIF(NX.EQ.NGXMAX) THEN
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

      DO NS=1,NSMAX
         PN(NS)=PNSAVE(NS)
      ENDDO

      CALL PAGES
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.3,0.)
      CALL GMNMX1(GX,1,NGXMAX,1,GXMIN,GXMAX)
      CALL GMNMX1(GY,1,NGYMAX,1,GYMIN,GYMAX)
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
         CALL CONTQ2(GZ(1:NGXMAX,1:NGYMAX,NGP),GX,GY,NGXMAX,NGXMAX,NGYMAX, &
                     0.0,2.0,1,0,0,KA)
      ENDDO
      CALL SETLIN(0,0,7)

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

      DO NGP=1,NGPMAX
         CALL MOVE(20.0,12.5-NGP*0.5)
         CALL SETLIN(0,0,7-MOD(NGP-1,5))
         CALL TEXT('PN  =',5)
         CALL NUMBD(PNNGP(NGP),'(1PE11.3)',11)
      ENDDO
      CALL SETLIN(0,0,7)
      CALL PAGEE
      DEALLOCATE(GX,GY,GZ,Z,KA,PNSAVE,PNNGP)
      GOTO 1

 9000 RETURN
  END SUBROUTINE DP_CONTX
END MODULE dpcont
