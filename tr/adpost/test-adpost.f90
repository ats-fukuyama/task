!
! Test program to read ADPOST file
!
PROGRAM test_adpost

  USE ADPOST
  USE libgrf
  IMPLICIT NONE
  integer:: ID,IZ0,NXMAX,NX,IERR
  REAL(rkind):: PTMIN,PTMAX,DPT,PT
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: XDATA
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: FDATA

  ID   =    1
  IZ0  =   74
  PTMIN=-3.d0
  PTMAX= 2.d0
  NXMAX=  200

  CALL GSOPEN
  CALL READ_ADPOST(IERR)
  IF(IERR.NE.0) THEN
     WRITE(6,'(A,I4)') 'XX test_adpost: READ_ADPOST: IERR =',IERR
     STOP
  END IF

1 CONTINUE
  WRITE(6,'(A)') '## Input Atomic Number IZ0 (0 for quit):'
  READ(5,*,ERR=1,END=9000) IZ0
  IF(IZ0.LE.0) GOTO 9000
  IF(.NOT.AVAILABLE_IZ0(IZ0)) GOTO 1

2 CONTINUE
  WRITE(6,'(A)') '## Input Data Type,LOG10_PTMIN,LO10_PTMAX,NXMAX:'
  WRITE(6,'(A)') '             1:Z_AV, 2:Z^2_av, 3:I_rad, 0:end'
  WRITE(6,'(A,I5,1P2E12.4,I5)') 'ID,PTMIN,PTMAX,NXMAX=', &
                                 ID,PTMIN,PTMAX,NXMAX
  READ(5,*,ERR=2,END=1) ID,PTMIN,PTMAX,NXMAX
  IF(ID.LE.0) GOTO 1
  IF(NXMAX.LE.1) GOTO 2

  ALLOCATE(XDATA(NXMAX),FDATA(NXMAX,1))

  DPT=(PTMAX-PTMIN)/NXMAX
  DO NX=1,NXMAX
     PT=PTMIN+DPT*(NX-1)
     XDATA(NX)=PT
     FDATA(NX,1)=FUNC_ADPOST(IZ0,ID,10.D0**PT)
!     WRITE(6,'(A,I5,1P3E12.4)') 'NX,PT,10.D0**PT,FDATA(NX,1)=', &
!                                 NX,PT,10.D0**PT,FDATA(NX,1)
  END DO

  CALL PAGES
  SELECt CASE(ID)
  CASE(1)
     CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,1,'@Z_av vs log_PT@',1)
  CASE(2)
     CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,1,'@Z^2_av vs log_PT@',1)
  CASE(3)
     CALL GRD1D(0,XDATA,FDATA,NXMAX,NXMAX,1,'@log_I_rad vs log_PT@',1)
  END SELECt
  CALL PAGEE
  DEALLOCATE(XDATA,FDATA)
  GO TO 2

9000 CONTINUE
  STOP
END PROGRAM test_adpost
