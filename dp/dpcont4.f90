MODULE dpcont4

  PRIVATE
  PUBLIC dp_cont4

CONTAINS

!     ***** PLOT CONTOUR GRAPH *************************************

  SUBROUTINE DP_CONT4(NID_)

    USE dpcomm_local
    USE plprof
    USE plprofw
    USE PLPARM,ONLY: pl_parm,pl_view
    USE dpparm
    USE dpdisp
    USE dpglib
    USE libchar
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NID_
    REAL,ALLOCATABLE:: GX(:),GY(:),GZ(:,:)
    REAL(rkind),ALLOCATABLE:: RFNORM(:),RKNORM(:)
    REAL(rkind),ALLOCATABLE:: Z(:,:)
    INTEGER,ALLOCATABLE:: KA(:,:,:)
    INTEGER,ALLOCATABLE:: NXA(:)
    REAL(rkind),ALLOCATABLE:: XA(:),RFA(:),RFIA(:)
    TYPE(pl_mag_type):: mag
    TYPE(pl_prfw_type),DIMENSION(nsmax):: plfw
    CHARACTER(LEN=1):: KID
    INTEGER,SAVE:: INIT=0
    INTEGER:: NID
    INTEGER:: NX,NY,NS,NGULEN,IERR
    INTEGER:: NGXMAX_SAVE,NGYMAX_SAVE,INFO
    REAL(rkind):: XMIN,XMAX,YMIN,YMAX,DX,DY,Y,X
    REAL(rkind):: RX,RY,RZ,Y1,Y2,VAL1,VAL2,Y0,RF,RFI,V
    REAL(rkind):: RHON,WC,WP2,VA,VT,XN,YN,YNI
    COMPLEX(rkind):: CRF,CKX,CKY,CKZ,CD,CD0,CD1
    REAL:: GUCLIP,F
    REAL:: GXMIN,GXMAX,GYMIN,GYMAX,GXSMN,GXSMX,GSCALX,GYSMN,GYSMX,GSCALY
    INTEGER:: ncount,ncount_max
    REAL(rkind):: vmin,vmax
    EXTERNAL PAGES,SETLIN,SETCHS,GMNMX1,GQSCAL,MOVE,TEXT,NUMBD,PAGEE
    EXTERNAL GDEFIN,GFRAME,GSCALE,GVALUE,CONTQ2,SETMKS,SETRGB,MARK2D

      IF(INIT.EQ.0) THEN
         KID='1'
         INIT=1
      ENDIF
      NID=NID_

    1 WRITE(6,*)'## SELECT: X-var: 1/KX 2/KY 3/KZ 4/X 5/Y 6/Z 7/K: ', &
                'P,D,V/parm A,B/C/type X/EXIT'
      READ(5,*,ERR=1,END=9000) KID
      CALL toupper(KID)
      IF(KID.EQ.'X') THEN
         GOTO 9000
      ELSEIF(KID.EQ.'P') THEN
         CALL DP_PARM(0,'DP',IERR)
         GOTO 1
      ELSEIF(KID.EQ.'A') THEN
         NID=4
         GOTO 1
      ELSEIF(KID.EQ.'B') THEN
         NID=5
         GOTO 1
      ELSEIF(KID.EQ.'C') THEN
         NID=6
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
               CALL PL_PROFW3D(RX,RY,RZ,plfw)
            ELSE
               CALL PL_PROFW(RHON,plfw)
            END IF
            IF(NORMK.LT.0) THEN
               VA=1.D0
               DO NS=1,NSMAX
                  WC=mag%BABS*PZ(NS)*AEE/(AMP*PA(NS))
                  WP2=plfw(NS)%RN*1.D20*PZ(NS)*PZ(NS)*AEE*AEE &
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
               VT=SQRT(plfw(NS)%RTPP*AEE*1.D3/(AMP*PA(NS)))
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

         ! --- count fliping points ---

      ncount=0
      DO NX=1,NGXMAX
         VAL1=Z(NX,1)
         DO NY=2,NGYMAX
            VAL2=Z(NX,NY)
            IF(VAL1*VAL2.LT.0.D0) ncount=ncount+1
            VAL1=VAL2
         END DO
      END DO
      ncount_max=ncount
      ALLOCATE(NXA(ncount_max),XA(ncount_max),RFA(ncount_max),RFIA(ncount_max))

      ncount=0
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

               CRF=CMPLX(Y0,0.D0)
               CD0=CFDISP(CRF,CKX0,CKY0,CKZ0,XPOS0,YPOS0,ZPOS0)
               CALL DPBRENTX(CRF,INFO)
               CD1=CFDISP(CRF,CKX0,CKY0,CKZ0,XPOS0,YPOS0,ZPOS0)
!               IF(ABS(CKX0).GT.1900.D0.AND.ABS(CKX0).LT.2000.D0.AND. &
!                  REAL(CRF).GT.210.D0.AND.REAL(CRF).LT.220.D0) THEN
!                  WRITE(6,'(A,5ES12.4)') &
!                       'CD0:',Y0,0.D0,REAL(CKX0),CD0
!                  WRITE(6,'(A,5ES12.4)') &
!                       'CD1:',CRF,REAL(CKX0),CD1
!               END IF
               IF(IDEBUG.GE.2) THEN
                  WRITE(22,'(A,2I5,1P5E11.3)')    'RT:',NX,NY,X,Y0,0.D0,CD0
                  WRITE(22,'(A,2I5,1P5E11.3,I5)') '   ',NX,NY,X,CRF,    CD1,INFO
               END IF
               RF=DBLE(CRF)
               RFI=AIMAG(CRF)
               IF(INFO.GE.1.AND.INFO.LE.3.AND. &     ! converged
                    RF.GE.YMIN.AND.RF.LE.YMAX.AND. & ! not out-of-range
                    ABS(CD1).LE.EPSDP.AND. &         ! small residue
                    RFI*RFNORM(NX).GT.-EPSRF) THEN   ! not strong damping
                  ncount=ncount+1                    ! these data are saved
                  NXA(ncount)=NX
                  XA(ncount)=X
                  RFA(ncount)=RF
                  RFIA(ncount)=RFI
               END IF
            END IF
            VAL1=VAL2
         END DO
      END DO
      ncount_max=ncount

      ! --- evaluate min and max of rfi/rf ---
      
      vmin=0.D0
      vmax=0.D0
      DO ncount=1,ncount_max
         v=RFIA(ncount)/RFA(ncount)
         vmin=MIN(v,vmin)
         vmax=MAX(v,vmax)
      END DO
      WRITE(6,'(A,I6,2ES12.4)') 'rfi/rf min/max:',ncount_max,vmin,vmax
      IF(vmax.LT. 1.D-8) vmax= 1.D-8
      IF(vmin.GT.-1.D-8) vmin=-1.D-8

      ! --- File output ---
      
      IF(NFLOUT.NE.0) THEN
         DO ncount=1,ncount_max
            NX=NXA(ncount)
            X=XA(ncount)
            RF=RFA(ncount)
            RFI=RFIA(ncount)
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
                  
            WRITE(NFLOUT,'(1P6E15.7,I5,1P2E15.7)') &
                 X,RF,RFI,XN,YN,YNI,INFO,CD1
         END DO
      END IF

      SELECT CASE(NID)
      CASE(4,5)
         CALL PAGES
         CALL SETLIN(0,0,7)
         CALL SETCHS(0.3,0.)
         CALL SETMKS(3,0.1)
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

         DO ncount=1,ncount_max
            NX=NXA(ncount)
            X=XA(ncount)
            RF=RFA(ncount)
            RFI=RFIA(ncount)
!            IF(RFI.LE.-0.01D0*ABS(RF)) CYCLE
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

            v=RFI/RF
            IF(v.GT.1.D-8) THEN
               f=GUCLIP(MIN(v/vmax,1.D0))
               CALL SETRGB(1.0,0.8*(1.0-f),0.0) ! 0 to 1: yellow to red
            ELSE IF(v.LT.-1.D-8) THEN
               f=GUCLIP(MIN(-v/MIN(vmax,EPSRF),1.D0))
               CALL SETRGB(0.0,0.5*(1.0-f),1.0)   ! 0 to -1: green to blue
            ELSE
               CALL SETRGB(0.8,0.8,0.8)
            END IF

            CALL MARK2D(GUCLIP(XN),GUCLIP(YN))
         END DO
      
         CALL DP_CONT4_PARM(KID,vmin,vmax,RFNORM(1),RKNORM(1))
         CALL PAGEE

      CASE(6)
         CALL PAGES
         CALL SETLIN(0,0,7)
         CALL SETCHS(0.3,0.)
         CALL SETMKS(3,0.1)
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

         DO ncount=1,ncount_max
            NX=NXA(ncount)
            X=XA(ncount)
            RF=RFA(ncount)
            RFI=RFIA(ncount)
            SELECT CASE(KID)
            CASE('1','2','3','7')
               XN=X*RKNORM(NX)
               YN=RF*RFNORM(NX)
            CASE('4','5','6')
               XN=X
               YN=RF*RFNORM(NX)
            END SELECT
                  
            CALL SETRGB(0.8,0.8,0.8)
            CALL MARK2D(GUCLIP(XN),GUCLIP(YN))

            v=RFI/RF
            IF(v.GT.1.D-8) THEN
               f=GUCLIP(MIN(v/vmax,1.D0))
               CALL SETRGB(1.0,0.8*(1.0-f),0.0) ! 0 to 1: yellow to red
               CALL MARK2D(GUCLIP(XN),GUCLIP(YN+RFI*RFNORM(NX)))
            ELSE IF(v.LT.-1.D-8) THEN
               f=GUCLIP(MIN(-v/MIN(vmax,EPSRF),1.D0))
               CALL SETRGB(0.0,0.5*(1.0-f),1.0)   ! 0 to -1: green to blue
               CALL MARK2D(GUCLIP(XN),GUCLIP(YN+RFI*RFNORM(NX)))
            END IF

         END DO
         CALL DP_CONT4_PARM(KID,vmin,vmax,RFNORM(1),RKNORM(1))
         CALL PAGEE
      END SELECT

      DEALLOCATE(NXA,XA,RFA,RFIA)
      DEALLOCATE(GX,GY,GZ,Z,KA,RFNORM,RKNORM)
      GOTO 2

 9000 RETURN
    END SUBROUTINE DP_CONT4

 
    SUBROUTINE DP_CONT4_PARM(KID,vmin,vmax,RFNORM,RKNORM)
      USE dpcomm_local
      IMPLICIT NONE
      CHARACTER(LEN=1),INTENT(IN):: KID
      REAL(rkind),INTENT(IN):: vmin,vmax
      REAL(rkind):: RFNORM,RKNORM
      EXTERNAL SETRGB,MOVE,TEXT,NUMBD,NUMBI

      CALL SETRGB(0.0,0.0,0.0)
      
      CALL MOVE(3.0,17.5)
      CALL TEXT ('RF',2)
      CALL TEXT(' vs ',4)
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
      IF(NORMF.GE.1.AND.NORMF.LE.NSMAX) THEN
         CALL TEXT('  fc=',5)
         CALL NUMBD(1.D0/RFNORM,'(ES11.3)',11)
      END IF
      IF(NORMK.GE.1.AND.NORMK.LE.NSMAX) THEN
         CALL TEXT('  Lc=',5)
         CALL NUMBD(RKNORM,'(ES11.3)',11)
      ELSEIF(NORMK.LE.-1.AND.NORMK.GE.-NSMAX) THEN
         CALL TEXT('  LA=',5)
         CALL NUMBD(RKNORM,'(ES11.3)',11)
      END IF

      CALL MOVE(20.0,17.5)
      CALL TEXT('BB  =',5)
      CALL NUMBD(BB,'(1PE11.3)',11)

      CALL MOVE(20.0,16.5)
      CALL TEXT('PN1 =',5)
      CALL NUMBD(PN(1),'(1PE11.3)',11)
      CALL MOVE(20.0,16.0)
      CALL TEXT('PTPR=',5)
      CALL NUMBD(PTPR(1),'(1PE11.3)',11)
      CALL MOVE(20.0,15.5)
      CALL TEXT('PTPP=',5)
      CALL NUMBD(PTPP(1),'(1PE11.3)',11)
      CALL MOVE(20.0,15.0)
      CALL TEXT('PUPR=',5)
      CALL NUMBD(PUPR(1),'(1PE11.3)',11)
      CALL MOVE(20.0,14.5)
      CALL TEXT('PUPP=',5)
      CALL NUMBD(PUPP(1),'(1PE11.3)',11)
      CALL MOVE(20.0,14.0)
      CALL TEXT('MODELP(1)=',10)
      CALL NUMBI(MODELP(1),'(I6)',6)

      CALL MOVE(20.0,13.0)
      CALL TEXT('PN2 =',5)
      CALL NUMBD(PN(2),'(1PE11.3)',11)
      CALL MOVE(20.0,12.5)
      CALL TEXT('PTPR=',5)
      CALL NUMBD(PTPR(2),'(1PE11.3)',11)
      CALL MOVE(20.0,12.0)
      CALL TEXT('PTPP=',5)
      CALL NUMBD(PTPP(2),'(1PE11.3)',11)
      CALL MOVE(20.0,11.5)
      CALL TEXT('PUPR=',5)
      CALL NUMBD(PUPR(2),'(1PE11.3)',11)
      CALL MOVE(20.0,11.0)
      CALL TEXT('PUPP=',5)
      CALL NUMBD(PUPP(2),'(1PE11.3)',11)
      CALL MOVE(20.0,10.5)
      CALL TEXT('MODELP(2)=',10)
      CALL NUMBI(MODELP(2),'(I6)',6)

      CALL MOVE(20.0, 9.5)
      CALL TEXT('PN3 =',5)
      CALL NUMBD(PN(3),'(1PE11.3)',11)
      CALL MOVE(20.0, 9.0)
      CALL TEXT('PTPR=',5)
      CALL NUMBD(PTPR(3),'(1PE11.3)',11)
      CALL MOVE(20.0, 8.5)
      CALL TEXT('PTPP=',5)
      CALL NUMBD(PTPP(3),'(1PE11.3)',11)
      CALL MOVE(20.0, 8.0)
      CALL TEXT('PUPR=',5)
      CALL NUMBD(PUPR(3),'(1PE11.3)',11)
      CALL MOVE(20.0, 7.5)
      CALL TEXT('PUPP=',5)
      CALL NUMBD(PUPP(3),'(1PE11.3)',11)
      CALL MOVE(20.0, 7.0)
      CALL TEXT('MODELP(3)=',10)
      CALL NUMBI(MODELP(3),'(I6)',6)

      SELECT CASE(KID)
      CASE('1','2','3','7')
         SELECT CASE(KID)
         CASE('1')
            CALL MOVE(20.0, 6.0)
            CALL TEXT('KY  =',5)
            CALL NUMBD(RKY0,'(1PE11.3)',11)
            CALL MOVE(20.0, 5.5)
            CALL TEXT('KZ  =',5)
            CALL NUMBD(RKZ0,'(1PE11.3)',11)
         CASE('2')
            CALL MOVE(20.0, 6.0)
            CALL TEXT('KX  =',5)
            CALL NUMBD(RKX0,'(1PE11.3)',11)
            CALL MOVE(20.0, 5.5)
            CALL TEXT('KZ  =',5)
            CALL NUMBD(RKZ0,'(1PE11.3)',11)
         CASE('3')
            CALL MOVE(20.0, 6.0)
            CALL TEXT('KX  =',5)
            CALL NUMBD(RKX0,'(1PE11.3)',11)
            CALL MOVE(20.0, 5.5)
            CALL TEXT('KY  =',5)
            CALL NUMBD(RKY0,'(1PE11.3)',11)
         CASE('7')
            CALL MOVE(20.0, 6.0)
            CALL TEXT('ANG =',5)
            CALL NUMBD(RKANG0,'(1PE11.3)',11)
            CALL MOVE(20.0, 5.5)
            CALL TEXT('KZ  =',5)
            CALL NUMBD(RKZ0,'(1PE11.3)',11)
         END SELECT
         CALL MOVE(20.0, 5.0)
         CALL TEXT('X   =',5)
         CALL NUMBD(RX0,'(1PE11.3)',11)
         CALL MOVE(20.0, 4.5)
         CALL TEXT('Y   =',5)
         CALL NUMBD(RY0,'(1PE11.3)',11)
         CALL MOVE(20.0, 4.0)
         CALL TEXT('Z   =',5)
         CALL NUMBD(RZ0,'(1PE11.3)',11)
      CASE('4','5','6')
         CALL MOVE(20.0, 6.0)
         CALL TEXT('KX  =',5)
         CALL NUMBD(RKX0,'(1PE11.3)',11)
         CALL MOVE(20.0, 5.5)
         CALL TEXT('KY  =',5)
         CALL NUMBD(RKY0,'(1PE11.3)',11)
         CALL MOVE(20.0, 5.0)
         CALL TEXT('KZ  =',5)
         CALL NUMBD(RKZ0,'(1PE11.3)',11)
         SELECT CASE(KID)
         CASE('4')
            CALL MOVE(20.0, 4.5)
            CALL TEXT('Y   =',5)
            CALL NUMBD(RY0,'(1PE11.3)',11)
            CALL MOVE(20.0, 4.0)
            CALL TEXT('Z   =',5)
            CALL NUMBD(RZ0,'(1PE11.3)',11)
         CASE('5')
            CALL MOVE(20.0, 4.5)
            CALL TEXT('X   =',5)
            CALL NUMBD(RX0,'(1PE11.3)',11)
            CALL MOVE(20.0, 4.0)
            CALL TEXT('Z   =',5)
            CALL NUMBD(RZ0,'(1PE11.3)',11)
         CASE('6')
            CALL MOVE(20.0, 4.5)
            CALL TEXT('X   =',5)
            CALL NUMBD(RX0,'(1PE11.3)',11)
            CALL MOVE(20.0, 4.0)
            CALL TEXT('Y   =',5)
            CALL NUMBD(RY0,'(1PE11.3)',11)
         END SELECT
      END SELECT
         
      CALL MOVE(20.0, 3.0)
      CALL TEXT('MODEL_ES =',10)
      CALL NUMBI(MODEL_ES,'(I6)',6)
      CALL MOVE(20.0, 2.5)
      CALL TEXT('Imax=',5)
      CALL NUMBD(vmax,'(ES11.3)',11)
      CALL MOVE(20.0, 2.0)
      CALL TEXT('Imin=',5)
      CALL NUMBD(vmin,'(ES11.3)',11)
      RETURN
    END SUBROUTINE DP_CONT4_PARM


 !     *************************
!           BRENT METHOD 2
!     *************************

  SUBROUTINE DPBRENTX(CX,INFO)

    USE dpcomm_local
    USE libbrent
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
