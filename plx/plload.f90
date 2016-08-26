!   pl module for loading profile data

MODULE plp2D
  USE bpsd_kinds
  INTEGER(ikind):: NXMAX,NYMAX
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: XD,YD
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: VA
  REAL(rkind),DIMENSION(:,:,:,:,:),ALLOCATABLE:: UA
END MODULE plp2D

MODULE plload

  private
  public pl_load,pl_read_xprf,pl_read_p2D,pl_read_p2Dmag

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

!----  Debug write

 8000 FORMAT(' N ',3X,'PRFRHO',6X,'PRFNE',6X,'PRFNI1',5X,'PRFNI2', &
                   5X,'PRFNI3',5X,'PRFNI4')
 8010 FORMAT(' N ',3X,'PRFRHO',6X,'PRFTE',6X,'PRFTI1',5X,'PRFTI2', &
                   5X,'PRFTI3',5X,'PRFTI4')
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

!     ***** LOAD 2D PROFILE DATA *****

    SUBROUTINE pl_load_p2D(ierr)

      USE plcomm,ONLY: ikind,rkind,idebug,KNAMPF
      USE plp2D
      USE libgrf
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: ierr
      INTEGER:: NFL,NX,NY,NV,IX,IY,IERSPL
      REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: FX,FY,FXY
      CHARACTER(LEN=80),SAVE:: KNAMPF_SAVE=' '

      ierr = 0
      IF(KNAMPF.EQ.KNAMPF_SAVE) RETURN
      KNAMPF_SAVE=KNAMPF

!----    open defined by toroidal magnetic flux

      NFL=22
      CALL FROPEN(NFL,KNAMPF,1,0,'PF',ierr) 
                                       !open to read formatted without pronpt
      IF(ierr.NE.0) GOTO 9990
      READ(NFL,*,END=9991) ! skip 4 lines: rf,kx,ky,kz
      READ(NFL,*,END=9991)
      READ(NFL,*,END=9991)
      READ(NFL,*,END=9991)
      READ(NFL,'(8X,I4)',ERR=9992,END=9991) NXMAX
      READ(NFL,'(8X,I4)',ERR=9993,END=9991) NYMAX

      IF(ALLOCATED(XD)) DEALLOCATE(XD,YD,VA,UA)
      ALLOCATE(XD(NXMAX),YD(NYMAX),VA(NXMAX,NYMAX,13))
      ALLOCATE(UA(4,4,NXMAX,NYMAX,13))
      ALLOCATE(FX(NXMAX,NYMAX),FY(NXMAX,NYMAX),FXY(NXMAX,NYMAX))

      READ(NFL,*,END=9991) ! skip 1 line: variable name
      DO NX=1,NXMAX
         DO NY=1,NYMAX
            READ(NFL,'(I3,I4,13E18.7)',ERR=9994,END=9991) &
                 IX,IY,(VA(NX,NY,NV),NV=1,13)
!            WRITE(6,'(I3,I4,13E18.7)') IX,IY,(VA(NX,NY,NV),NV=1,13)
            IF(IX.NE.NX) THEN
               WRITE(6,'(A,2I10)') 'IX,NX=',IX,NX
               GO TO 9995
            END IF
            IF(IY.NE.NY) THEN
               WRITE(6,'(A,2I10)') 'IY,NY=',IY,NY
               GO TO 9996
            END IF
         END DO
         READ(NFL,*,END=9991) ! skip 1 line: next NX
      END DO

      DO NX=1,NXMAX
         XD(NX)=VA(NX,1,1)
      END DO
      DO NY=1,NYMAX
         YD(NY)=VA(1,NY,2)
      END DO

      WRITE(6,'(A)') '## profile data successfully loaded'
      WRITE(6,'(A)') '   profile data file: '
      WRITE(6,'(A)') KNAMPF

      WRITE(6,'(A,1P4E12.4)') 'XMIN,XMAX,YMIN,YMAX =', &
           XD(1),XD(NXMAX),YD(1),YD(NYMAX)

      IF(IDEBUG.EQ.1) THEN
         CALL PAGES
         CALL GRD2D(14,XD,YD,VA(1:NXMAX,1:NYMAX,1),NXMAX,NXMAX,NYMAX,'@XD@')
         CALL GRD2D(15,XD,YD,VA(1:NXMAX,1:NYMAX,2),NXMAX,NXMAX,NYMAX,'@YD@')
         CALL GRD2D(16,XD,YD,VA(1:NXMAX,1:NYMAX,3),NXMAX,NXMAX,NYMAX,'@BX@')
         CALL GRD2D(17,XD,YD,VA(1:NXMAX,1:NYMAX,4),NXMAX,NXMAX,NYMAX,'@BY@')
         CALL GRD2D(18,XD,YD,VA(1:NXMAX,1:NYMAX,5),NXMAX,NXMAX,NYMAX,'@BZ@')
         CALL GRD2D(19,XD,YD,VA(1:NXMAX,1:NYMAX,6),NXMAX,NXMAX,NYMAX,'@BTOT@')
         CALL GRD2D(20,XD,YD,VA(1:NXMAX,1:NYMAX,7),NXMAX,NXMAX,NYMAX,'@TE@')
         CALL GRD2D(21,XD,YD,VA(1:NXMAX,1:NYMAX,8),NXMAX,NXMAX,NYMAX,'@NE@')
         CALL GRD2D(22,XD,YD,VA(1:NXMAX,1:NYMAX,9),NXMAX,NXMAX,NYMAX,'@TI@')
         CALL GRD2D(23,XD,YD,VA(1:NXMAX,1:NYMAX,10),NXMAX,NXMAX,NYMAX,'@NI@')
         CALL GRD2D(24,XD,YD,VA(1:NXMAX,1:NYMAX,11),NXMAX,NXMAX,NYMAX,'@X@')
         CALL GRD2D(25,XD,YD,VA(1:NXMAX,1:NYMAX,12),NXMAX,NXMAX,NYMAX,'@Y@')
         CALL GRD2D(26,XD,YD,VA(1:NXMAX,1:NYMAX,13),NXMAX,NXMAX,NYMAX,'@Z@')
         CALL PAGEE
      END IF

!----  Set coefficient for spline

      DO NV=3,11
         CALL SPL2D(XD,YD,VA(1,1,NV),FX,FY,FXY,UA(1,1,1,1,NV), &
                    NXMAX,NXMAX,NYMAX,0,0,IERSPL)
         IF (IERSPL.NE.0) THEN
            WRITE(6,*) 'XX pl_load_p2D: SPL2D: NV=', NV
            GO TO 9997
         END IF
      END DO
      GO TO 9999

 9990 WRITE(6,*) '==========  pl_load_p2D OPEN ERROR  =========='
      WRITE(6,*) '   IERR = ',IERR
      GO TO 9999
 9991 WRITE(6,*) '==========  pl_load_p2D ABNORMAL END of FILE  =========='
      GO TO 9999
 9992 WRITE(6,*) '==========  pl_load_p2D FILE READ ERROR: NXMAX  =========='
      GO TO 9999
 9993 WRITE(6,*) '==========  pl_load_p2D FILE READ ERROR: NYMAX  =========='
      GO TO 9999
 9994 WRITE(6,*) '==========  pl_load_p2D FILE READ ERROR: VA  =========='
      GO TO 9999
 9995 WRITE(6,*) '==========  pl_load_p2D INCONSISTENT DATA: IX != NX'
      GO TO 9999
 9996 WRITE(6,*) '==========  pl_load_p2D INCONSISTENT DATA: IY != NY'
      GO TO 9999
 9997 WRITE(6,*) '==========  pl_load_p2D SPL2D ERROR  =========='
      WRITE(6,*) '   IERSPL = ',IERSPL
      GO TO 9999

 9999 REWIND NFL
      CLOSE(NFL)
      IF(ALLOCATED(FX)) DEALLOCATE(FX,FY,FXY)
      RETURN
    END SUBROUTINE pl_load_p2D

!     ***** 2D magnetic field profile *****

    SUBROUTINE pl_read_p2Dmag(X,Y,BX,BY,BZ,IERR)

      USE plcomm,ONLY: rkind,ikind,NSMAX
      USE plp2d
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
      CALL SPL2DF(XL,YL,BX,XD,YD,UA(1,1,1,1,3),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8001
      CALL SPL2DF(XL,YL,BY,XD,YD,UA(1,1,1,1,4),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8002
      CALL SPL2DF(XL,YL,BZ,XD,YD,UA(1,1,1,1,5),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8003

      RETURN
    END SUBROUTINE pl_read_p2Dmag

!     ***** 2D density and temperature profile *****

    SUBROUTINE pl_read_p2D(X,Y,RNPL,RTPL,RUPL,NSMAXL,IERR)

      USE plcomm,ONLY: rkind,ikind,NSMAX
      USE plp2d
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X,Y    ! Position
      REAL(rkind),DIMENSION(NSMAX),INTENT(OUT):: &
           RNPL,  &! Density [10^{20} m^{-3}]
           RTPL,  &! Temperature [keV]
           RUPL    ! Flow velosity [m/s]
      INTEGER(ikind),INTENT(OUT):: &
           NSMAXL,&! Number of provided particle species
           IERR    ! ERROR Indicator 
      REAL(rkind):: XL,YL
      INTEGER(ikind):: IERL

      XL=X
      IF(XL.LT.XD(1))     XL=XD(1)
      IF(XL.GT.XD(NXMAX)) XL=XD(NXMAX)
      YL=Y
      IF(YL.LT.YD(1))     YL=YD(1)
      IF(YL.GT.YD(NYMAX)) YL=YD(NYMAX)

      NSMAXL=2

      IERR=0
      CALL SPL2DF(XL,YL,RNPL(1),XD,YD,UA(1,1,1,1, 8),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8001
      RNPL(1)=RNPL(1)*1.D-20
      CALL SPL2DF(XL,YL,RTPL(1),XD,YD,UA(1,1,1,1, 7),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8002
      RTPL(1)=RTPL(1)*1.D-3
      RUPL(1)=0.D0
      CALL SPL2DF(XL,YL,RNPL(2),XD,YD,UA(1,1,1,1,10),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8003
      RNPL(2)=RNPL(2)*1.D-20
      CALL SPL2DF(XL,YL,RTPL(2),XD,YD,UA(1,1,1,1, 9),NXMAX,NXMAX,NYMAX,IERL)
      IF(IERL.NE.0) IERR=8004
      RTPL(2)=RTPL(2)*1.D-3
      RUPL(2)=0.D0
         
      RETURN
    END SUBROUTINE pl_read_p2D

  END MODULE plload
