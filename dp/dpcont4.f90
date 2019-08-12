MODULE dpcont4

  PRIVATE
  PUBLIC dp_cont4

CONTAINS

!     ***** PLOT CONTOUR GRAPH *************************************

  SUBROUTINE DP_CONT4(NID)

    USE dpcomm_local
    USE plprof
    USE PLPARM,ONLY: pl_parm,pl_view
    USE dpparm
    USE dpdisp
    USE dpglib
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NID
    REAL(4),ALLOCATABLE:: GX(:),GY(:),GZ(:,:)
    REAL(rkind),ALLOCATABLE:: RFNORM(:),RKNORM(:)
    REAL(rkind),ALLOCATABLE:: Z(:,:)
    INTEGER,ALLOCATABLE:: KA(:,:,:)
    TYPE(pl_mag_type):: mag
    TYPE(pl_plf_type),DIMENSION(nsmax):: plf
    CHARACTER(LEN=1):: KID
    INTEGER,SAVE:: INIT=0
    INTEGER:: NX,NY,NS,NY1,NY2,NX1,NX2,NGULEN,IERR
    INTEGER:: NGXMAX_SAVE,NGYMAX_SAVE,INFO
    REAL(rkind):: XMIN,XMAX,YMIN,YMAX,DX,DY,Y,X,FACTY,FACTX,DZY,DZX,DZ2
    REAL(rkind):: RX,RY,RZ,Y1,Y2,VAL1,VAL2,VAL3,Y0,RF,RFI,V,VL,RFPREV
    REAL(rkind):: RHON,WC,WCI,WP2,VA,VT,XN,YN,YNI
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CD,CW,CWC,CD0,CD1
    REAL(4):: GUCLIP,F
    REAL(4):: GXMIN,GXMAX,GYMIN,GYMAX,GXSMN,GXSMX,GSCALX,GYSMN,GYSMX,GSCALY

      IF(INIT.EQ.0) THEN
         KID='1'
         INIT=1
      ENDIF

    1 WRITE(6,*)'## SELECT : X-var : 1/KX 2/KY 3/KZ 4/X 5/Y 6/Z 7/K: ', &
                'P,D,V/PARM X/EXIT'
      READ(5,*,ERR=1,END=9000) KID
      CALL GUCPTL(KID)
      IF(KID.EQ.'X') THEN
         GOTO 9000
      ELSEIF(KID.EQ.'P') THEN
         CALL DP_PARM(0,'DP',IERR)
         GOTO 1
      ELSEIF(KID.EQ.'V') THEN
         CALL PL_VIEW
         CALL DP_VIEW
         GOTO 1
      ENDIF

    2 IF(KID.EQ.'1') THEN
         XMIN=RKX1
         XMAX=RKX2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RKX1,RKX2,NGXMAX',RKX1,RKX2,NGXMAX
      ELSE IF(KID.EQ.'2') THEN
         XMIN=RKY1
         XMAX=RKY2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RKY1,RKY2,NGXMAX',RKY1,RKY2,NGXMAX
      ELSE IF(KID.EQ.'3') THEN
         XMIN=RKZ1
         XMAX=RKZ2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RKZ1,RKZ2,NGXMAX',RKZ1,RKZ2,NGXMAX
      ELSE IF(KID.EQ.'4') THEN
         XMIN=RX1
         XMAX=RX2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RX1,RX2,NGXMAX',RX1,RX2,NGXMAX
      ELSE IF(KID.EQ.'5') THEN
         XMIN=RY1
         XMAX=RY2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RY1,RY2,NGXMAX',RY1,RY2,NGXMAX
      ELSE IF(KID.EQ.'6') THEN
         XMIN=RX1
         XMAX=RX2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RZ1,RZ2,NGXMAX',RZ1,RZ2,NGXMAX
      ELSE IF(KID.EQ.'7') THEN
         XMIN=RK1
         XMAX=RK2
         WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT X-AXIS: RK1,RK2,NGXMAX',RK1,RK2,NGXMAX
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
      IF(KID.EQ.'1') THEN
         RKX1=XMIN
         RKX2=XMAX
      ELSE IF(KID.EQ.'2') THEN
         RKY1=XMIN
         RKY2=XMAX
      ELSE IF(KID.EQ.'3') THEN
         RKZ1=XMIN
         RKZ2=XMAX
      ELSE IF(KID.EQ.'4') THEN
         RX1=XMIN
         RX2=XMAX
      ELSE IF(KID.EQ.'5') THEN
         RY1=XMIN
         RY2=XMAX
      ELSE IF(KID.EQ.'6') THEN
         RZ1=XMIN
         RZ2=XMAX
      ELSE IF(KID.EQ.'7') THEN
         RK1=XMIN
         RK2=XMAX
      END IF
      
    3 YMIN=RF1
      YMAX=RF2
      WRITE(6,'(A,1P2E12.4,I8)') &
              ' INPUT Y-AXIS: RF1,RF2,NGYMAX',RF1,RF2,NGYMAX
      NGYMAX_SAVE=NGYMAX
      READ(5,*,ERR=3,END=2) YMIN,YMAX,NGYMAX
      IF(NGYMAX.LE.0) THEN
         NGYMAX=NGYMAX_SAVE
         GOTO 2
      END IF
      RF1=YMIN
      RF2=YMAX

      ALLOCATE(GX(NGXMAX),GY(NGYMAX),GZ(NGXMAX,NGYMAX))
      ALLOCATE(RFNORM(NGXMAX),RKNORM(NGXMAX))
      ALLOCATE(Z(NGXMAX,NGYMAX),KA(8,NGXMAX,NGYMAX))

      DY=(YMAX-YMIN)/(NGYMAX-1)
      DX=(XMAX-XMIN)/(NGXMAX-1)

      CKX=RKX0
      CKY=RKY0
      CKZ=RKZ0
      RX=RX0
      RY=RY0
      RZ=RZ0

      DO NX=1,NGXMAX
         X=XMIN+DX*(NX-1)
         IF(KID.EQ.'1') THEN
            CKX=X
         ELSEIF(KID.EQ.'2') THEN
            CKY=X
         ELSEIF(KID.EQ.'3') THEN
            CKZ=X
         ELSEIF(KID.EQ.'4') THEN
            RX=X
         ELSEIF(KID.EQ.'5') THEN
            RY=X
         ELSEIF(KID.EQ.'6') THEN
            RZ=X
         ELSEIF(KID.EQ.'7') THEN
            CKX=X*SIN(RKANG0*PI/180.D0)
            CKY=X*COS(RKANG0*PI/180.D0)
         ENDIF

         IF(NORMF.NE.0.OR.NORMK.NE.0) THEN
            CALL PL_MAG(RX,RY,RZ,mag)
            RHON=mag%rhon
            IF(MODELG.LE.1.OR.MODELG.GT.10) THEN
               CALL PL_PROF3D(RX,RY,RZ,plf)
            ELSE
               CALL PL_PROF(RHON,plf)
            END IF
            IF(NORMK.LT.0) THEN
               VA=1.D0
               DO NS=1,NSMAX
                  WC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS))
                  WP2=plf(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE &
                       /(EPS0*AMP*PA(NS))
                  VA=VA+WP2/WC**2
               END DO
               VA=VC/SQRT(VA)
            END IF
            IF(NORMF.GE.1.AND.NORMF.LE.NSMAX) THEN
               NS=NORMF
               WC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS))
               RFNORM(NX)=2.D0*PI*1.D6/WC
            ELSE
               RFNORM(NX)=1.D0
            END IF
            IF(NORMK.GE.1.AND.NORMK.LE.NSMAX) THEN
               NS=NORMK
               WC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS))
               VT=SQRT(plf(NS)%RTPP*AEE*1.D3/(AMP*PA(NS)))
               RKNORM(NX)=VT/WC
            ELSEIF(-NORMK.GE.1.AND.-NORMK.LE.NSMAX) THEN
               NS=-NORMK
               WC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS))
               RKNORM(NX)=VA/WC
            ELSE
               RKNORM(NX)=1.D0
            END IF
         ELSE
            RFNORM(NX)=1.D0
            RKNORM(NX)=1.D0
         END IF

         DO NY=1,NGYMAX
            Y=YMIN+DY*(NY-1)
            CRF=DCMPLX(Y,0.D0)
            CD=CFDISP(CRF,CKX,CKY,CKZ,RX,RY,RZ)

            SELECT CASE(KID)
            CASE('1','2','3','7')
               GX(NX)=GUCLIP(X*RKNORM(NX))
               GY(NY)=GUCLIP(Y*RFNORM(NX))
            CASE('4','5','6')
               GX(NX)=GUCLIP(X)
               GY(NY)=GUCLIP(Y*RFNORM(NX))
            END SELECT
            Z(NX,NY)=DBLE(CD)
         ENDDO
      ENDDO

      DO NY=1,NGYMAX
      DO NX=1,NGXMAX
         GZ(NX,NY)=GUCLIP(Z(NX,NY))
      ENDDO
      ENDDO

      IF(IDEBUG.GE.1) THEN
         DO NY=1,NGYMAX
            DO NX=1,NGXMAX
               IF(MOD(NX-1,10).EQ.0.AND.MOD(NY-1,10).EQ.0) THEN
                  WRITE(6,'(2I5,1P4E15.7)')  &
                       NX,NY,GX(NX),GY(NY),Z(NX,NY),GZ(NX,NY)
               ENDIF
            ENDDO
         ENDDO
      END IF

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
      IF(NID.EQ.4) THEN
         CALL CONTQ2(GZ,GX,GY,NGXMAX,NGXMAX,NGYMAX,0.0,2.0,1,0,0,KA)
      END IF

      CALL SETMKS(3,0.1)
      RFPREV=0.D0
      DO NX=1,NGXMAX
         VAL1=Z(NX,1)
         DO NY=2,NGYMAX
            VAL2=Z(NX,NY)
            IF(VAL1*VAL2.LT.0.D0) THEN
               X=XMIN+DX*(NX-1)
               Y1=YMIN+DY*(NY-1)
               Y2=YMIN+DY* NY
               Y0=Y1-(Y2-Y1)*VAL1/(VAL2-VAL1)
               CKX0=RKX0
               CKY0=RKY0
               CKZ0=RKZ0
               XPOS0=RX0
               YPOS0=RY0
               ZPOS0=RZ0
               IF(KID.EQ.'1') THEN
                  CKX0=X
               ELSEIF(KID.EQ.'2') THEN
                  CKY0=X
               ELSEIF(KID.EQ.'3') THEN
                  CKZ0=X
               ELSEIF(KID.EQ.'4') THEN
                  XPOS0=X
               ELSEIF(KID.EQ.'5') THEN
                  YPOS0=X
               ELSEIF(KID.EQ.'6') THEN
                  ZPOS0=X
               ELSEIF(KID.EQ.'7') THEN
                  CKX0=X*SIN(RKANG0*PI/180.D0)
                  CKY0=X*COS(RKANG0*PI/180.D0)
               ENDIF

               CRF=DCMPLX(Y0,0.D0)
               CD0=CFDISP(CRF,CKX0,CKY0,CKZ0,XPOS0,YPOS0,ZPOS0)
               CALL DPBRENTX(CRF,INFO)
               CD1=CFDISP(CRF,CKX0,CKY0,CKZ0,XPOS0,YPOS0,ZPOS0)
               IF(IDEBUG.GE.2) THEN
                  WRITE(22,'(A,2I5,1P5E11.3)')    'RT:',NX,NY,X,Y0,0.D0,CD0
                  WRITE(22,'(A,2I5,1P5E11.3,I5)') '   ',NX,NY,X,CRF,    CD1,INFO
               END IF
               RF=DBLE(CRF)
               RFI=AIMAG(CRF)
               IF(INFO.GE.1.AND.INFO.LE.3.AND. &
                    RF.GE.YMIN.AND.RF.LE.YMAX.AND. &
                    ABS(CD1).LE.EPSDP.AND. &
                    ABS(RF-RFPREV).GT.EPSRF) THEN
                  V=RFI/RF
                  IF(V.GT.1.D-12) THEN
                     VL=MAX(MIN(0.D0,LOG10(V)),-12.D0)+12.D0  ! 0<VL<12
                     f=GUCLIP(VL/12.D0)
                     CALL SETRGB(f,1.0-f,0.0)
                  ELSE IF(V.LT.-1.D-12) THEN
                     VL=MAX(MIN(0.D0,LOG10(-V)),-12.D0)+12.D0  ! 0<VL<12
                     f=GUCLIP(VL/12.D0)
                     CALL SETRGB(0.0,1.0-f,f)
                  ELSE
                     VL=0.D0
                     f=0.0
                     CALL SETRGB(0.0,1.0,0.0)
                  END IF
                  
                  SELECT CASE(KID)
                  CASE('1','2','3','7')
                     XN=X*RKNORM(NX)
                     YN=RF*RFNORM(NX)
                     YNI=RFI*RFNORM(NX)
                  CASE('4','5','6')
                     XN=X
                     YN=RF*RFNORM(NX)
                     YNI=RFI*RFNORM(NX)
                  END SELECT
                  
                  CALL MARK2D(GUCLIP(XN),GUCLIP(YN))
                  IF(NFLOUT.EQ.21) THEN
                     WRITE(21,'(1P6E15.7,I5,1P2E15.7)') &
                          X,RF,RFI,XN,YN,YNI,INFO,CD1
                    RFPREV=RF
                 END IF
               END IF
               RFPREV=RF
            END IF
            VAL1=VAL2
         END DO
      END DO
      CALL SETRGB(0.0,0.0,0.0)

      CALL MOVE(3.0,17.5)
      IF(KID.EQ.'1') THEN
         CALL TEXT ('KX',2)
      ELSEIF(KID.EQ.'2') THEN
         CALL TEXT ('KY',2)
      ELSEIF(KID.EQ.'3') THEN
         CALL TEXT ('KZ',2)
      ELSEIF(KID.EQ.'4') THEN
         CALL TEXT ('X',1)
      ELSEIF(KID.EQ.'5') THEN
         CALL TEXT ('Y',1)
      ELSEIF(KID.EQ.'6') THEN
         CALL TEXT ('Z',1)
      ELSEIF(KID.EQ.'7') THEN
         CALL TEXT ('K',1)
      ENDIF
      CALL TEXT('-',1)
      CALL TEXT ('RF',2)

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
      CALL MOVE(20.0,12.5)
      CALL TEXT('ANG =',5)
      CALL NUMBD(RKANG0,'(1PE11.3)',11)
      CALL PAGEE
      DEALLOCATE(GX,GY,GZ,Z,KA,RFNORM,RKNORM)
      GOTO 2

 9000 RETURN
 END SUBROUTINE DP_CONT4

!     *************************
!           BRENT METHOD 2
!     *************************

  SUBROUTINE DPBRENTX(CX,INFO)

    USE dpcomm_local
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(INOUT):: CX
    INTEGER,INTENT(OUT):: INFO
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: X,F,WA
    INTEGER:: M,NBRMAX,LWA

    M = LMAXRT
    NBRMAX = 2
    LWA=NBRMAX*(NBRMAX+3)
    ALLOCATE(X(NBRMAX),F(NBRMAX),WA(LWA))
    X(1)=DBLE(CX)
    X(2)=AIMAG(CX)
    CALL FBRENTN(DPFUNCX,NBRMAX,X,F,EPSRT,INFO,WA,LWA)
    CX=DCMPLX(X(1),X(2))
    DEALLOCATE(X,F,WA)
    RETURN
 END SUBROUTINE DPBRENTX

  SUBROUTINE DPFUNCX(N,X,F,ID)

    USE dpcomm_local
    USE dpdisp
    IMPLICIT NONE
    INTEGER,INTENT(IN):: N,ID
    REAL(rkind),INTENT(IN):: X(N)
    REAL(rkind),INTENT(INOUT):: F(N)
    COMPLEX(rkind):: CRF,CFUNC

    CRF=DCMPLX(X(1),X(2))
    CFUNC=CFDISP(CRF,CKX0,CKY0,CKZ0,XPOS0,YPOS0,ZPOS0)
    IF(ID.EQ.1)THEN
       F(1)=DBLE(CFUNC)
    ELSE
       F(2)=DIMAG(CFUNC)
    END IF
    IF(ILIST.NE.0) WRITE(6,'(1P3E12.4,I5)') X(1),X(2),F(ID),ID
    RETURN
  END SUBROUTINE DPFUNCX

END MODULE dpcont4
