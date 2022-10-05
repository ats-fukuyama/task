Module wimparm

  PRIVATE

  PUBLIC wim_parm, wim_view

CONTAINS

!     ****** PARAMETER INPUT ******

  SUBROUTINE wim_parm(MODE,KIN,IERR)

!     MODE=0 : standard namelinst input
!     MODE=1 : namelist file input
!     MODE=2 : namelist line input

!     IERR=0 : normal end
!     IERR=1 : namelist standard input error
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE
!     IERR=10X : input parameter out of range

    USE libkio
    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode
    CHARACTER(LEN=*),INTENT(IN)::  kin
    INTEGER,INTENT(OUT):: ierr

1   CALL TASK_PARM(MODE,'WIM',KIN,wim_nlin,wim_plst,IERR)
    IF(IERR.NE.0 .AND. IERR.NE.2) RETURN

    CALl wim_check(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
    IF(IERR.NE.0) IERR=IERR+100

    RETURN
  END SUBROUTINE wim_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE wim_nlin(NID,IST,IERR)

    USE wimcomm

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NID
    INTEGER,INTENT(OUT) :: IST,IERR
    INTEGER,PARAMETER:: NSM=100

    NAMELIST /WIM/ NZMAX,NWMAX,ZMIN,ZMAX,PN0,DBDZ,ANX,BETA,NTMAX,TMAX, &
                   CER1,CEL1,CER2,CEL2,DZMAX,DZWID,PZCL,MODELW,MODELP, &
                   IDEBUG

    READ(NID,WIM,IOSTAT=IST,ERR=9800,END=9900)

    IERR=0
    RETURN

9800 IERR=8
    RETURN
9900 IERR=9
    RETURN
  END SUBROUTINE wim_nlin

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE wim_plst

    IMPLICIT NONE
    WRITE(6,'(A/)') &
    '# &WIM : NZMAX,NWMAX,ZMIN,ZMAX,PN0,DBDZ,ANX,BETA,NTMAX,IDEBUG,',&
    '         TMAX,CER1,CEL1,CER2,CEL2,DZMAX,DZWID,PZCL,MODELW,MODELP'
    RETURN

  END SUBROUTINE wim_plst

!     ****** CHECK INPUT PARAMETER ******

  SUBROUTINE wim_check(IERR)

    USE wimcomm
    IMPLICIT NONE
    INTEGER:: IERR

    IERR=0

    IF(NZMAX .LE. 0) THEN
       WRITE(6,'(A,I8)') 'WIM wim_check: INVALID NZMAX: NZMAX=',NZMAX
       IERR=1
    ENDIF
    RETURN
  END SUBROUTINE wim_check

!     ****** SHOW PARAMETERS ******

  SUBROUTINE wim_view

    use wimcomm
    implicit none

    WRITE(6,601) 'NZMAX   ',NZMAX   ,'NWMAX   ',NWMAX   , &
                 'NTMAX   ',NTMAX   ,'IDEBUG  ',IDEBUG
    WRITE(6,601) 'MODELW  ',MODELW  ,'MODELP  ',MODELP
    WRITE(6,602) 'ZMIN    ',ZMIN    ,'ZMAX    ',ZMAX    , &
                 'PN0     ',PN0     ,'DBDZ    ',DBDZ    , &
                 'ANX     ',ANX     ,'BETA    ',BETA    , &
                 'TMAX    ',TMAX    ,'DZMAX   ',DZMAX   , &
                 'DZWID   ',DZWID   ,'PZCL    ',PZCL
    WRITE(6,603) 'CER1    ',CER1    ,'CEL1    ',CEL1    ,&
                 'CER2    ',CER2    ,'CEL2    ',CEL2   
    RETURN

601 FORMAT(' ',A8,'=',I8:2X,A8,'=',I8: &
            2X,A8,'=',I8:2X,A8,'=',I8)
602 FORMAT(' ',A8,'=',1PE12.4:2X,A8,'=',1PE12.4:2X,A8,'=',1PE12.4)
603 FORMAT(' ',A8,'=',1P2E12.4:4X,A8,'=',1P2E12.4)
  END SUBROUTINE wim_view

END MODULE wimparm
  
