MODULE dpcont

  PRIVATE
  PUBLIC dp_cont2,dp_cont3

CONTAINS

!     ***** PLOT CONTOUR GRAPH *************************************

  SUBROUTINE DP_CONT2

    USE dpcomm_local
    USE plprof
    USE PLPARM,ONLY: pl_parm,pl_view
    USE dpparm
    USE dpdisp
    IMPLICIT NONE
    REAL(4),ALLOCATABLE:: GX(:),GY(:),GZ(:,:)
    REAL(rkind),ALLOCATABLE:: Z(:,:)
    INTEGER,ALLOCATABLE:: KA(:,:,:)
    TYPE(pl_mag_type):: mag
    CHARACTER(LEN=2):: KID
    CHARACTER(LEN=1):: KID1,KID2
    INTEGER,SAVE:: INIT=0
    INTEGER:: NX,NY,NS,NY1,NY2,NX1,NX2,NGULEN,IERR
    INTEGER:: NGXMAX_SAVE,NGYMAX_SAVE
    REAL(rkind):: XMIN,XMAX,YMIN,YMAX,DX,DY,Y,X,FACTY,FACTX,DZY,DZX,DZ2
    REAL(rkind):: RX,RY,RZ
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CD,CW,CWC
    REAL(4):: GUCLIP
    REAL(4):: GXMIN,GXMAX,GYMIN,GYMAX,GXSMN,GXSMX,GSCALX,GYSMN,GYSMX,GSCALY

      IF(INIT.EQ.0) THEN
         KID='YW'
         NGXMAX=51
         NGYMAX=51
         INIT=1
      ENDIF

    1 WRITE(6,*)'## SELECT : XY : W/RF X/KX Y/KY Z/KZ 1/X 2/Y 3/Z: ', &
                'P,D,V/PARM 0/EXIT'
      READ(5,*,ERR=1,END=9000) KID
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
         CALL DP_PARM(0,'DP',IERR)
         GOTO 1
      ELSEIF(KID1.EQ.'V') THEN
         CALL PL_VIEW
         CALL DP_VIEW
         GOTO 1
      ENDIF

    2 IF(KID1.EQ.'W') THEN
         XMIN=RF1
         XMAX=RF2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RF1,RF2,NGXMAX',RF1,RF2,NGXMAX
      ELSE IF(KID1.EQ.'X') THEN
         XMIN=RKX1
         XMAX=RKX2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RKX1,RKX2,NGXMAX',RKX1,RKX2,NGXMAX
      ELSE IF(KID1.EQ.'Y') THEN
         XMIN=RKY1
         XMAX=RKY2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RKY1,RKY2,NGXMAX',RKY1,RKY2,NGXMAX
      ELSE IF(KID1.EQ.'Z') THEN
         XMIN=RKZ1
         XMAX=RKZ2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RKZ1,RKZ2,NGXMAX',RKZ1,RKZ2,NGXMAX
      ELSE IF(KID1.EQ.'1') THEN
         XMIN=RX1
         XMAX=RX2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RX1,RX2,NGXMAX',RX1,RX2,NGXMAX
      ELSE IF(KID1.EQ.'2') THEN
         XMIN=RY1
         XMAX=RY2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RY1,RY2,NGXMAX',RY1,RY2,NGXMAX
      ELSE IF(KID1.EQ.'3') THEN
         XMIN=RX1
         XMAX=RX2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RZ1,RZ2,NGXMAX',RZ1,RZ2,NGXMAX
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID1'
         GOTO 1
      ENDIF
      NGXMAX_SAVE=NGXMAX
      READ(5,*,ERR=2,END=1) XMIN,XMAX,NGXMAX
      IF(NGXMAX.LE.0) THEN
         NGXMAX=NGXMAX_SAVE
         GOTO 1
      END IF
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
      ELSE IF(KID1.EQ.'1') THEN
         RX1=XMIN
         RX2=XMAX
      ELSE IF(KID1.EQ.'2') THEN
         RY1=XMIN
         RY2=XMAX
      ELSE IF(KID1.EQ.'3') THEN
         RZ1=XMIN
         RZ2=XMAX
      END IF
      
    3 IF(KID2.EQ.'W') THEN
         YMIN=RF1
         YMAX=RF2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT Y-AXIS: RF1,RF2,NGYMAX',RF1,RF2,NGYMAX
      ELSE IF(KID2.EQ.'X') THEN
         YMIN=RKX1
         YMAX=RKX2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT Y-AXIS: RKX1,RKX2,NGYMAX',RKX1,RKX2,NGYMAX
      ELSE IF(KID2.EQ.'Y') THEN
         YMIN=RKY1
         YMAX=RKY2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT Y-AXIS: RKY1,RKY2,NGYMAX',RKY1,RKY2,NGYMAX
      ELSE IF(KID2.EQ.'Z') THEN
         YMIN=RKZ1
         YMAX=RKZ2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT Y-AXIS: RKZ1,RKZ2,NGYMAX',RKZ1,RKZ2,NGYMAX
      ELSE IF(KID2.EQ.'1') THEN
         YMIN=RX1
         YMAX=RX2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT Y-AXIS: RX1,RX2,NGYMAX',RX1,RX2,NGYMAX
      ELSE IF(KID2.EQ.'2') THEN
         YMIN=RY1
         YMAX=RY2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT Y-AXIS: RY1,RY2,NGYMAX',RY1,RY2,NGYMAX
      ELSE IF(KID2.EQ.'3') THEN
         YMIN=RZ1
         YMAX=RZ2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT Y-AXIS: RZ1,RZ2,NGYMAX',RZ1,RZ2,NGYMAX
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID2'
         GOTO 1
      ENDIF
      NGYMAX_SAVE=NGYMAX
      READ(5,*,ERR=3,END=2) YMIN,YMAX,NGYMAX
      IF(NGYMAX.LE.0) THEN
         NGYMAX=NGYMAX_SAVE
         GOTO 2
      END IF
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
      ELSE IF(KID2.EQ.'1') THEN
         RX1=YMIN
         RX2=YMAX
      ELSE IF(KID2.EQ.'2') THEN
         RY1=YMIN
         RY2=YMAX
      ELSE IF(KID2.EQ.'3') THEN
         RZ1=YMIN
         RZ2=YMAX
      END IF

      ALLOCATE(GX(NGXMAX),GY(NGYMAX),GZ(NGXMAX,NGYMAX))
      ALLOCATE(Z(NGXMAX,NGYMAX),KA(8,NGXMAX,NGYMAX))

      DY=(YMAX-YMIN)/(NGYMAX-1)
      DX=(XMAX-XMIN)/(NGXMAX-1)

      CRF=DCMPLX(RF0,RFI0)
      CKX=RKX0
      CKY=RKY0
      CKZ=RKZ0
      RX=RX0
      RY=RY0
      RZ=RZ0

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
         ELSEIF(KID2.EQ.'1') THEN
            RX=Y
         ELSEIF(KID2.EQ.'2') THEN
            RY=Y
         ELSEIF(KID2.EQ.'3') THEN
            RZ=Y
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
            ELSEIF(KID1.EQ.'1') THEN
               RX=X
            ELSEIF(KID1.EQ.'2') THEN
               RY=X
            ELSEIF(KID1.EQ.'3') THEN
               RZ=X
            ENDIF
            CD=CFDISP(CRF,CKX,CKY,CKZ,RX,RY,RZ)

! --- remove sign change due to cyc

            CW=2.D0*PI*1.D6*CRF
            IF(ABS(CW).LE.1.D-8) CW=(1.D-8,0.D0)
            CALL pl_mag(RX,RY,RZ,mag)
            DO NS=1,NSMAX
               CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
!               CD=CD*(1.D0-CWC)
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

!      DO NY=1,NGYMAX
!         DO NX=1,NGXMAX
!            IF(MOD(NX-1,10).EQ.0.AND.MOD(NY-1,10).EQ.0) THEN
!               WRITE(6,'(2I5,1P4E15.7)')  &
!                    NX,NY,GX(NX),GY(NY),Z(NX,NY),GZ(NX,NY)
!            ENDIF
!         ENDDO
!      ENDDO

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
      ELSEIF(KID1.EQ.'1') THEN
         CALL TEXT ('X',1)
      ELSEIF(KID1.EQ.'2') THEN
         CALL TEXT ('Y',1)
      ELSEIF(KID1.EQ.'Z') THEN
         CALL TEXT ('Z',1)
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
      ELSEIF(KID2.EQ.'1') THEN
         CALL TEXT ('X',1)
      ELSEIF(KID2.EQ.'2') THEN
         CALL TEXT ('Y',1)
      ELSEIF(KID2.EQ.'3') THEN
         CALL TEXT ('Z',1)
      ENDIF

      CALL MOVE(20.0,16.5)
      CALL TEXT('RF  =',5)
      CALL NUMBD(RF0,  '(1PE11.3)',11)
      CALL MOVE(20.0,16.0)
      CALL TEXT('RFI=',5)
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
      CALL MOVE(20.0,14.0)
      CALL TEXT('X   =',5)
      CALL NUMBD(RX0,'(1PE11.3)',11)
      CALL MOVE(20.0,13.5)
      CALL TEXT('Y   =',5)
      CALL NUMBD(RY0,'(1PE11.3)',11)
      CALL MOVE(20.0,13.0)
      CALL TEXT('Z   =',5)
      CALL NUMBD(RZ0,'(1PE11.3)',11)
      CALL PAGEE
      DEALLOCATE(GX,GY,GZ,Z,KA)
      GOTO 2

 9000 RETURN
  END SUBROUTINE DP_CONT2

!     ***** PLOT CONTOUR GRAPH *****
!     ***** DENSITY SCAN *****

  SUBROUTINE DP_CONT3

    USE dpcomm_local
    USE plprof
    USE plinit
    USE dpparm
    USE dpdisp
    IMPLICIT NONE
    REAL(4),ALLOCATABLE:: GX(:),GY(:),GZ(:,:,:)
    REAL(rkind),ALLOCATABLE:: Z(:,:,:),PNSAVE(:),PNNGP(:)
    INTEGER,ALLOCATABLE:: KA(:,:,:)
    TYPE(pl_mag_type):: mag
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
      READ(5,*,ERR=1,END=9000) KID
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
            CALL pl_mag(RX0,RY0,RZ0,mag)
            CW=2.D0*PI*1.D6*CRF
            DO NS=1,NSMAX
               CWC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS)*CW)
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
  END SUBROUTINE DP_CONT3
END MODULE dpcont
