!   pl module for loading profile data

MODULE plp2D
  USE bpsd_kinds
  INTEGER(ikind):: NXMAX,NYMAX
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: XD,YD
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: VA
  REAL(rkind),DIMENSION(:,:,:,:,:),ALLOCATABLE:: UA
END MODULE plp2D

MODULE pl_trdata
  USE bpsd_kinds
  INTEGER(ikind):: NRMAX_TR,NSMAX_TR,NFMAX_TR
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: RM_TR,RG_TR,DERIV
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: RN_TR,RT_TR,RW_TR,RNF_TR,RTF_TR
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: URN_TR,URT_TR
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: URW_TR,URNF_TR,URTF_TR
END MODULE pl_trdata

MODULE plload

  private
  public pl_load,pl_read_xprf,pl_read_p2D,pl_read_p2Dmag, &
         pl_load_trdata,pl_read_trdata

CONTAINS

  SUBROUTINE pl_load(ierr)

    USE plcomm,ONLY: modeln,modelg
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
  
    SELECT CASE(modeln)
    CASE(8)
       CALL pl_load_xprf(ierr)
    CASE(12)
       CALL pl_load_p2D(ierr)
    CASE(21)
       CALL pl_load_trdata(0,ierr)
    END SELECT

    SELECT CASE(modelg)
    CASE(12)
       CALL pl_load_p2D(ierr)
    END SELECT

    RETURN
  END SUBROUTINE pl_load

!     ***** LOAD PROFILE DATA from TOPICS *****

    SUBROUTINE pl_load_xprf(ierr)

      USE plcomm,ONLY: NSMAX,PZ,KNAMPF
      USE plxprf
      USE libspl1d
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: ierr
      REAL(rkind),DIMENSION(NXPRF,NXSPC):: PRFN,PRFT
      INTEGER(ikind):: ifno,nr,ns,irc
      REAL(rkind):: val
      CHARACTER(LEN=80),SAVE:: KNAMPF_SAVE=' '

      ierr = 0
      IF(KNAMPF.EQ.KNAMPF_SAVE) RETURN
      KNAMPF_SAVE=KNAMPF

!----  Open profile data file and read
!----  PRFNE, PRFTE is data at the point divided equally by rho 
!----    defined by toroidal magnetic flux

      IFNO=22
      OPEN ( IFNO, FILE=KNAMPF, ERR=9995 )
      READ ( IFNO, '(I3)', END=9996, ERR=9996 ) NPRFMAX
      DO NR=1,NPRFMAX
         READ ( IFNO, '(13E14.7)', END=9996, ERR=9996 ) &
             PRFRHO(NR), (PRFN(NR,NS), NS=1,NXSPC), &
                         (PRFT(NR,NS), NS=1,NXSPC)
      ENDDO

!----  Modification for charge neutrality

      DO NR=1,NPRFMAX
         VAL=0.D0
         DO NS=2,NSMAX
            VAL=VAL+PZ(NS)*PRFN(NR,NS)
         ENDDO
         PRFN(NR,1)=VAL
      ENDDO

!----  Set coefficient for spline

      DO NS=1,NSMAX
         CALL SPL1D(PRFRHO,PRFN(1,NS),DERIV,UPRFN(1,1,NS),NPRFMAX,0,IRC)
         IF (IRC.NE.0) GO TO 9997
         CALL SPL1D(PRFRHO,PRFT(1,NS),DERIV,UPRFT(1,1,NS),NPRFMAX,0,IRC)
         IF (IRC.NE.0) GO TO 9997
      ENDDO

      GO TO 9999

 9995 WRITE(6,*) '==========  pl_load_xprf FILE OPEN ERROR  =========='
      GO TO 9999
 9996 WRITE(6,*) '==========  pl_load_xprf FILE READ ERROR  =========='
      GO TO 9999
 9997 WRITE(6,*) '==========  pl_load_xprf SPL1D ERROR  =========='

 9999 CLOSE( IFNO )
      RETURN
    END SUBROUTINE pl_load_xprf

!     ***** Interpolation of profile at a given point *****

    SUBROUTINE pl_read_xprf(Rhol,NS,PNL,PTL)

      USE plcomm,ONLY: PNS,PTS,modeln
      USE plxprf
      USE libspl1d
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: rhol   ! Normalized radius
      INTEGER(ikind),INTENT(IN):: NS  ! Particle species
      REAL(rkind),INTENT(OUT):: PNL   ! Density at Rhol
      REAL(rkind),INTENT(OUT):: PTL   ! Temperature at Rhol
      REAL(rkind):: PPL
      INTEGER(ikind):: IERR

      IF (Rhol.GT.1.0D0) THEN
         IF(modeln.EQ.1) THEN
            PNL = PNS(NS)
         ELSE
            PNL = 0.D0
         END IF
         PTL = PTS(NS)
      ELSE
         CALL SPL1DF(Rhol,PPL,PRFRHO,UPRFN(1,1,NS),NPRFMAX,IERR)
         PNL=PPL
         CALL SPL1DF(Rhol,PPL,PRFRHO,UPRFT(1,1,NS),NPRFMAX,IERR)
         PTL=PPL
      ENDIF

      RETURN
    END SUBROUTINE pl_read_xprf

!     ***** LOAD PROFILE DATA from trdata *****

    SUBROUTINE pl_load_trdata(nid,ierr)

      USE pl_trdata
      USE libspl1d
      USE libfio
      IMPLICIT NONE
      INTEGER,INTENT(IN):: nid
      INTEGER,INTENT(OUT):: ierr
      INTEGER:: NFL,NR,NS,NF
      CHARACTER(LEN=80):: KFLNAM

      ierr = 0
      IF(ALLOCATED(RM_TR)) THEN
         DEALLOCATE(RM_TR,RG_TR,DERIV)
         DEALLOCATE(RN_TR,RT_TR,RW_TR,RNF_TR,RTF_TR)
         DEALLOCATE(URN_TR,URT_TR,URW_TR,URNF_TR,URTF_TR)
      END IF

!----  Open profile data file and read

      NFL=25
      SELECT CASE(nid)
      CASE(0)
         KFLNAM='trdata'
      CASE(1)
         KFLNAM='trdata1'
      CASE(2)
         KFLNAM='trdata2'
      CASE(3)
         KFLNAM='trdata3'
      CASE(4)
         KFLNAM='trdata4'
      CASE(5)
         KFLNAM='trdata5'
      CASE(6)
         KFLNAM='trdata6'
      END SELECT
      WRITE(6,*) 'KFLNAM=',KFLNAM
      CALL FROPEN(NFL,KFLNAM,0,0,'trdata',IERR)
      IF(IERR.NE.0) GO TO 9991
      READ(NFL) NRMAX_TR,NSMAX_TR,NFMAX_TR
      WRITE(6,'(A,3I5)') &
           'NRMAX_TR,NSMAX_TR,NFMAX_TR=', &
           NRMAX_TR,NSMAX_TR,NFMAX_TR

      ALLOCATE(RM_TR(NRMAX_TR),RG_TR(NRMAX_TR),DERIV(NRMAX_TR))
      ALLOCATE(RN_TR(NRMAX_TR,NSMAX_TR),RT_TR(NRMAX_TR,NSMAX_TR))
      ALLOCATE(RW_TR(NRMAX_TR,NSMAX_TR))
      ALLOCATE(RNF_TR(NRMAX_TR,NFMAX_TR),RTF_TR(NRMAX_TR,NFMAX_TR))
      ALLOCATE(URN_TR(4,NRMAX_TR,NSMAX_TR),URT_TR(4,NRMAX_TR,NSMAX_TR))
      ALLOCATE(URW_TR(4,NRMAX_TR,NSMAX_TR))
      ALLOCATE(URNF_TR(4,NRMAX_TR,NFMAX_TR),URTF_TR(4,NRMAX_TR,NFMAX_TR))

      READ(NFL) (RM_TR(NR),RG_TR(NR),NR=1,NRMAX_TR)
      READ(NFL) ((RN_TR(NR,NS),RT_TR(NR,NS),NR=1,NRMAX_TR),NS=1,NSMAX_TR)
      READ(NFL) ((RW_TR(NR,NF),RNF_TR(NR,NF),RTF_TR(NR,NF),NR=1,NRMAX_TR), &
                                                           NF=1,NFMAX_TR)
      CLOSE(NFL)
      WRITE(6,'(1P5E12.4)') (RM_TR(NR),RN_TR(NR,1),RT_TR(NR,1),&
                                       RN_TR(NR,2),RT_TR(NR,2),NR=1,NRMAX_TR)
      WRITE(6,*) '## Data loaded frm trdata'

!----  Set coefficient for spline

      DERIV(1:NRMAX_TR)=0.D0
      DO NS=1,NSMAX_TR
         CALL SPL1D(RM_TR,RN_TR(1:NRMAX_TR,NS),DERIV, &
                    URN_TR(1:4,1:NRMAX_TR,NS),NRMAX_TR,0,IERR)
         IF (IERR.NE.0) GO TO 9992
         CALL SPL1D(RM_TR,RT_TR(1:NRMAX_TR,NS),DERIV, &
                    URT_TR(1:4,1:NRMAX_TR,NS),NRMAX_TR,0,IERR)
         IF (IERR.NE.0) GO TO 9992
      ENDDO
      DO NF=1,NFMAX_TR
         CALL SPL1D(RM_TR,RW_TR(1:NRMAX_TR,NF),DERIV, &
                    URW_TR(1:4,1:NRMAX_TR,NF),NRMAX_TR,0,IERR)
         IF (IERR.NE.0) GO TO 9992
         CALL SPL1D(RM_TR,RNF_TR(1:NRMAX_TR,NF),DERIV, &
                    URNF_TR(1:4,1:NRMAX_TR,NF),NRMAX_TR,0,IERR)
         IF (IERR.NE.0) GO TO 9992
         CALL SPL1D(RM_TR,RTF_TR(1:NRMAX_TR,NF),DERIV, &
                    URTF_TR(1:4,1:NRMAX_TR,NF),NRMAX_TR,0,IERR)
         IF (IERR.NE.0) GO TO 9992
      ENDDO
      RETURN

9991  WRITE(6,*) '==========  pl_load_trdata FILE OPEN ERROR  =========='
      RETURN
9992  WRITE(6,*) '==========  pl_load_trdata SPLINE ERROR  =========='
      RETURN

    END SUBROUTINE pl_load_trdata

    SUBROUTINE pl_read_trdata(rho,NS,PNL,PTL)

      USE pl_trdata
      USE libspl1d
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: rho    ! Normalized radius
      INTEGER(ikind),INTENT(IN):: NS  ! Particle species
      REAL(rkind),INTENT(OUT):: PNL   ! Density at rho
      REAL(rkind),INTENT(OUT):: PTL   ! Temperature at rho
      REAL(rkind):: rhol
      INTEGER(ikind):: IERR

      rhol=MIN(MAX(RM_TR(1),rho),RM_TR(NRMAX_TR))
      IF(NS.LE.NSMAX_TR) THEN
         CALL SPL1DF(rhol,PNL,RM_TR,URN_TR(1:4,1:NRMAX_TR,NS),NRMAX_TR,IERR)
         CALL SPL1DF(rhol,PTL,RM_TR,URT_TR(1:4,1:NRMAX_TR,NS),NRMAX_TR,IERR)
      ELSE
         PNL=0.D0
         PTL=0.003D0
      END IF
!      WRITE(6,'(A,I5,1P3E12.4)') 'NZ,rho,PNL,PTL=',NS,rho,PNL,PTL

      RETURN
    END SUBROUTINE pl_read_trdata

!     ***** LOAD 2D PROFILE DATA *****

    SUBROUTINE pl_load_p2D(ierr)

      USE plcomm,ONLY: ikind,rkind,KNAMPF
      USE plp2D
      USE libspl2d
      USE libfio
      USE libgrf
      USE libmpi
      USE commpi
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: ierr
      INTEGER:: NFL,NX,NY,NV,IERSPL
      REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: FX,FY,FXY
      REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: VATEMP
      CHARACTER(LEN=80),SAVE:: KNAMPF_SAVE=' '
      INTEGER:: istatus

      ierr = 0
      istatus=0

      IF(nrank.EQ.0) THEN
         IF(KNAMPF.NE.KNAMPF_SAVE) THEN
            KNAMPF_SAVE=KNAMPF
            istatus=1

! --- open 2D profile data faile ---
            
            NFL=22
            CALL FROPEN(NFL,KNAMPF,1,0,'PF',ierr) 
                        ! open to read formatted without pronpt
            IF(ierr.NE.0) GOTO 9990
            READ(NFL,*,END=9991) ! skip 2 lines: rf,kx,ky,kz
            READ(NFL,*,END=9991)
            READ(NFL,*,ERR=9992,END=9991) NXMAX
            READ(NFL,*,END=9991)
            READ(NFL,*,ERR=9993,END=9991) NYMAX

            ALLOCATE(VATEMP(NXMAX,NYMAX,10))

            READ(NFL,*,END=9991) ! skip 1 line: variable name
            DO NX=1,NXMAX
               DO NY=1,NYMAX
                  READ(NFL,'(10E22.12)',ERR=9994,END=9991) &
                       (VATEMP(NX,NY,NV),NV=1,10)
               END DO
            END DO
            
            WRITE(6,'(A)') '## profile data successfully loaded'
            WRITE(6,'(A)') '   profile data file: '
            WRITE(6,'(A)') KNAMPF
            GO TO 9995

9990        WRITE(6,*) '==========  pl_load_p2D OPEN ERROR  ========'
            WRITE(6,*) '   IERR = ',IERR
            GO TO 9995
9991        WRITE(6,*) '==========  pl_load_p2D ABNORMAL END of FILE  ========'
            GO TO 9995
9992        WRITE(6,*) '==========  pl_load_p2D FILE READ ERROR: NXMAX  ======'
            GO TO 9995
9993        WRITE(6,*) '==========  pl_load_p2D FILE READ ERROR: NYMAX  ======'
            GO TO 9995
9994        WRITE(6,*) '==========  pl_load_p2D FILE READ ERROR: VA  ========='
            GO TO 9995
      
9995        REWIND NFL
            CLOSE(NFL)
            
            IF(ALLOCATED(XD)) DEALLOCATE(XD,YD,VA,UA)
            ALLOCATE(XD(NXMAX),YD(NYMAX),VA(NXMAX,NYMAX,8))
            ALLOCATE(UA(4,4,NXMAX,NYMAX,8))
            ALLOCATE(FX(NXMAX,NYMAX),FY(NXMAX,NYMAX),FXY(NXMAX,NYMAX))

            DO NX=1,NXMAX
               XD(NX)=VATEMP(NX,1,1)
            END DO
            DO NY=1,NYMAX
               YD(NY)=VATEMP(1,NY,2)
            END DO

            DO NV=1,8
               DO NY=1,NYMAX
                  DO NX=1,NXMAX
                     VA(NX,NY,NV)=VATEMP(NX,NY,NV+2)
                  END DO
               END DO
            END DO

            DEALLOCATE(VATEMP)
            
            WRITE(6,'(A,1P4E12.4)') 'XMIN,XMAX,YMIN,YMAX =', &
                 XD(1),XD(NXMAX),YD(1),YD(NYMAX)
         END IF
      END IF

      CALL mtx_broadcast1_integer(istatus)

      IF(istatus.EQ.1) THEN
         CALL mtx_broadcast1_integer(nxmax)
         CALL mtx_broadcast1_integer(nymax)

         IF(nrank.NE.0) THEN
            IF(ALLOCATED(XD)) DEALLOCATE(XD,YD,VA,UA,FX,FY,FXY)
            ALLOCATE(XD(NXMAX),YD(NYMAX),VA(NXMAX,NYMAX,8))
            ALLOCATE(UA(4,4,NXMAX,NYMAX,8))
            ALLOCATE(FX(NXMAX,NYMAX),FY(NXMAX,NYMAX),FXY(NXMAX,NYMAX))
         END IF

         CALL mtx_broadcast_real8(xd,nxmax)
         CALL mtx_broadcast_real8(yd,nymax)
         DO nv=1,8
            CALL mtx_broadcast2D_real8(va(1:nxmax,1:nymax,nv), &
                 nxmax,nxmax,nymax)
         END DO

!         IF(nrank.EQ.0) THEN
!            IF(IDEBUG.EQ.1) THEN
!               CALL PAGES
!               CALL GRD2D(14,XD,YD,VA(1:NXMAX,1:NYMAX,1), &
!                    NXMAX,NXMAX,NYMAX,'@BX@',MODE_2D=2)
!               CALL GRD2D(15,XD,YD,VA(1:NXMAX,1:NYMAX,2), &
!                    NXMAX,NXMAX,NYMAX,'@BY@',MODE_2D=2)
!               CALL GRD2D(16,XD,YD,VA(1:NXMAX,1:NYMAX,3), &
!                    NXMAX,NXMAX,NYMAX,'@BZ@',MODE_2D=2)
!               CALL GRD2D(17,XD,YD,VA(1:NXMAX,1:NYMAX,4), &
!                    NXMAX,NXMAX,NYMAX,'@BTOT@',MODE_2D=2)
!               CALL GRD2D(18,XD,YD,VA(1:NXMAX,1:NYMAX,5), &
!                    NXMAX,NXMAX,NYMAX,'@TE@',MODE_2D=2)
!               CALL GRD2D(19,XD,YD,VA(1:NXMAX,1:NYMAX,6), &
!                    NXMAX,NXMAX,NYMAX,'@NE@',MODE_2D=2)
!               CALL GRD2D(20,XD,YD,VA(1:NXMAX,1:NYMAX,7), &
!                    NXMAX,NXMAX,NYMAX,'@TI@',MODE_2D=2)
!               CALL GRD2D(21,XD,YD,VA(1:NXMAX,1:NYMAX,8), &
!                    NXMAX,NXMAX,NYMAX,'@NI@',MODE_2D=2)
!               CALL PAGEE
!            END IF
!         END IF

!----  Set coefficient for spline

         DO NV=1,8
            CALL SPL2D(XD,YD,VA(1,1,NV),FX,FY,FXY,UA(1,1,1,1,NV), &
                       NXMAX,NXMAX,NYMAX,0,0,IERSPL)
            IF (IERSPL.NE.0) THEN
               WRITE(6,*) 'XX pl_load_p2D: SPL2D: NV=', NV
               GO TO 9997
            END IF
         END DO
         GO TO 9999

9997     WRITE(6,*) '==========  pl_load_p2D SPL2D ERROR  =========='
         WRITE(6,*) '   IERSPL = ',IERSPL

9999     CONTINUE
      END IF
      RETURN
    END SUBROUTINE pl_load_p2D

!     ***** 2D magnetic field profile *****

    SUBROUTINE pl_read_p2Dmag(X,Y,BX,BY,BZ,IERR)

      USE plcomm,ONLY: rkind,ikind
      USE plp2d
      USE libspl2d
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,Y       ! Position
      REAL(rkind),INTENT(OUT):: BX,BY,BZ ! magnetic field
      INTEGER(ikind),INTENT(OUT):: IERR  ! ERROR Indicator 
      REAL(rkind):: XL,YL
      INTEGER(ikind):: IERL

      XL=X
      IF(XL.LT.XD(1))     XL=XD(1)
      IF(XL.GT.XD(NXMAX)) XL=XD(NXMAX)
      YL=Y
      IF(YL.LT.YD(1))     YL=YD(1)
      IF(YL.GT.YD(NYMAX)) YL=YD(NYMAX)

      IERR=0
      CALL SPL2DF(XL,YL,BX,XD,YD,UA(1,1,1,1,1),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8001
      CALL SPL2DF(XL,YL,BY,XD,YD,UA(1,1,1,1,2),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8002
      CALL SPL2DF(XL,YL,BZ,XD,YD,UA(1,1,1,1,3),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8003

      RETURN
    END SUBROUTINE pl_read_p2Dmag

!     ***** 2D density and temperature profile *****

    SUBROUTINE pl_read_p2D(X,Y,RN,RTPR,RTPP,RU,IERR)

      USE plcomm,ONLY: rkind,ikind,NSMAX
      USE plp2d
      USE libspl2d
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,Y    ! Position
      REAL(rkind),DIMENSION(NSMAX),INTENT(OUT):: &
           RN,    &! Density [10^{20} m^{-3}]
           RTPR,  &! Parallel Temperature [keV]
           RTPP,  &! Parallel Temperature [keV]
           RU      ! Flow velosity [m/s]
      INTEGER(ikind),INTENT(OUT):: &
           IERR    ! ERROR Indicator 
      REAL(rkind):: XL,YL,RN_PL,RT_PL
      INTEGER(ikind):: IERL

      XL=X
      IF(XL.LT.XD(1))     XL=XD(1)
      IF(XL.GT.XD(NXMAX)) XL=XD(NXMAX)
      YL=Y
      IF(YL.LT.YD(1))     YL=YD(1)
      IF(YL.GT.YD(NYMAX)) YL=YD(NYMAX)

      IERR=0
      CALL SPL2DF(XL,YL,RN_PL,XD,YD,UA(1,1,1,1, 6),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8001
      CALL SPL2DF(XL,YL,RT_PL,XD,YD,UA(1,1,1,1, 5),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8002
      RN(1)=RN_PL*1.D-20
      RTPR(1)=RT_PL*1.D-3
      RTPP(1)=RT_PL*1.D-3
      RU(1)=0.D0
      CALL SPL2DF(XL,YL,RN_PL,XD,YD,UA(1,1,1,1, 8),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8003
      CALL SPL2DF(XL,YL,RT_PL,XD,YD,UA(1,1,1,1, 7),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8004
      RN(2)=RN_PL*1.D-20
      RTPR(2)=RT_PL*1.D-3
      RTPP(2)=RT_PL*1.D-3
      RU(2)=0.D0
         
      RETURN
    END SUBROUTINE pl_read_p2D

  END MODULE plload
