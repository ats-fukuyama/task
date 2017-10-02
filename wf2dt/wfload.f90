! ***** FILE IO *****

! ***** Module for WG E filed data *****

MODULE wfwg1D

  USE bpsd_kinds
  INTEGER(ikind):: NYMAX
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: YD
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: VA
  REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: UA

END MODULE wfwg1D

MODULE wfload

  PRIVATE
  PUBLIC wf_load_wg,wf_read_wg

CONTAINS

!     ***** LOAD WAVEGUIDE DATA *****

  SUBROUTINE wf_load_wg(ierr)

      USE bpsd_kinds
      USE plcomm,ONLY: idebug
      USE wfcomm,ONLY: KNAMWG
      USE wfwg1D
      USE libgrf
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: ierr
      INTEGER:: NFL,NY,NV,IY,IERSPL,NYEND
      INTEGER,PARAMETER:: NYM=1000
      REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: VATEMP
      REAL(rkind),DIMENSION(:),ALLOCATABLE:: FX
      CHARACTER(LEN=80),SAVE:: KNAMWG_SAVE=' '

      ierr = 0
      IF(KNAMWG.EQ.KNAMWG_SAVE) RETURN
      KNAMWG_SAVE=KNAMWG

!----    open defined by toroidal magnetic flux

      NFL=22
      CALL FROPEN(NFL,KNAMWG,1,0,'WG',ierr) 
                                       !open to read formatted without pronpt
      IF(ierr.NE.0) GOTO 9990

      READ(NFL,*,END=9991) ! skip 2 lines: title,XYZ
      READ(NFL,*,END=9991)
      READ(NFL,*,ERR=9994,END=9991) NYMAX

      ALLOCATE(VATEMP(NYMAX,8))

      READ(NFL,*,END=9991) ! skip 1 line: variable name
      DO NY=1,NYMAX
         READ(NFL,'(8E22.12)',ERR=9994,END=9991) &
              (VATEMP(NY,NV),NV=1,8)
      END DO

      IF(ALLOCATED(YD)) DEALLOCATE(YD,VA,UA,FX)
      ALLOCATE(YD(NYMAX),VA(NYMAX,7))
      ALLOCATE(UA(4,4,NYMAX,7))
      ALLOCATE(FX(NYMAX))

      DO NY=1,NYMAX
         YD(NY)=VATEMP(NY,1)
         DO NV=1,7
            VA(NY,NV)=VATEMP(NY,NV+1)
         END DO
      END DO

      DEALLOCATE(VATEMP)

      WRITE(6,'(A)') '## WG 1D data successfully loaded'
      WRITE(6,'(A)') '   WG 1D data file: '
      WRITE(6,'(A)') KNAMWG

      WRITE(6,'(A,1P2E12.4)') 'YMIN,YMAX =',YD(1),YD(NYMAX)

      IF(IDEBUG.EQ.1) THEN
         CALL PAGES
         CALL GRD1D( 1,YD,VA(1,1),NYMAX,NYMAX,2,'@EX@')
         CALL GRD1D( 2,YD,VA(1,3),NYMAX,NYMAX,2,'@EY@')
         CALL GRD1D( 3,YD,VA(1,5),NYMAX,NYMAX,2,'@EZ@')
         CALL GRD1D( 4,YD,VA(1,7),NYMAX,NYMAX,1,'@EABS@')
         CALL PAGEE
      END IF

!----  Set coefficient for spline

      DO NV=1,6
         CALL SPL1D(YD,VA(1,NV),FX,UA(1,1,1,NV), &
                    NYMAX,0,IERSPL)
         IF (IERSPL.NE.0) THEN
            WRITE(6,*) 'XX wf_load_wg: SPL1D: NV=', NV
            GO TO 9997
         END IF
      END DO
      GO TO 9999

 9990 WRITE(6,*) '==========  wf_load_wg OPEN ERROR  =========='
      WRITE(6,*) '   IERR = ',IERR
      GO TO 9999
 9991 WRITE(6,*) '==========  wf_load_wg ABNORMAL END of FILE  =========='
      GO TO 9999
 9994 WRITE(6,*) '==========  wf_load_wg FILE READ ERROR: VA  =========='
      GO TO 9999
 9997 WRITE(6,*) '==========  wf_load_wg SPL1D ERROR  =========='
      WRITE(6,*) '   IERSPL = ',IERSPL
      GO TO 9999

 9999 REWIND NFL
      CLOSE(NFL)
      IF(ALLOCATED(FX)) DEALLOCATE(FX)
      RETURN
    END SUBROUTINE wf_load_wg

!     ***** 2D density and temperature profile *****

    SUBROUTINE wf_read_wg(Y,CEX,CEY,CEZ,IERR)

      USE bpsd_kinds
      USE wfcomm,ONLY: MODELWG
      USE wfwg1D
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: Y    ! Position
      COMPLEX(rkind),INTENT(OUT):: CEX,CEY,CEZ ! WG mouth electric field 
      INTEGER(ikind),INTENT(OUT):: IERR    ! ERROR Indicator 
      REAL(rkind):: YL,ER,EI
      INTEGER(ikind):: IERL

      YL=Y
      IF(YL.LT.YD(1))     YL=YD(1)
      IF(YL.GT.YD(NYMAX)) YL=YD(NYMAX)

      IERR=0
      CALL SPL1DF(YL,ER,YD,UA(1,1,1,1),NYMAX,IERL)
      IF(IERL.NE.0) IERR=8001
      CALL SPL1DF(YL,EI,YD,UA(1,1,1,2),NYMAX,IERL)
      IF(IERL.NE.0) IERR=8002
      CEX=CMPLX(ER,EI)
      CALL SPL1DF(YL,ER,YD,UA(1,1,1,3),NYMAX,IERL)
      IF(IERL.NE.0) IERR=8003
      CALL SPL1DF(YL,EI,YD,UA(1,1,1,4),NYMAX,IERL)
      IF(IERL.NE.0) IERR=8004
      CEY=CMPLX(ER,EI)
      CALL SPL1DF(YL,ER,YD,UA(1,1,1,5),NYMAX,IERL)
      IF(IERL.NE.0) IERR=8005
      CALL SPL1DF(YL,EI,YD,UA(1,1,1,6),NYMAX,IERL)
      IF(IERL.NE.0) IERR=8006
      CEZ=CMPLX(ER,EI)
         
      RETURN
    END SUBROUTINE wf_read_wg


  END MODULE wfload
