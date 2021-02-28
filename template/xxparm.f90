Module xxparm

  PRIVATE
  PUBLIC xx_parm, xx_view

CONTAINS

!     ****** PARAMETER INPUT ******

  SUBROUTINE xx_parm(MODE,KIN,IERR)

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

1   CALL TASK_PARM(MODE,'XX',KIN,xx_nlin,xx_plst,IERR)
    IF(IERR.NE.0 .AND. IERR.NE.2) RETURN

    CALl xx_check(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
    IF(IERR.NE.0) IERR=IERR+100

    RETURN
  END SUBROUTINE xx_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE xx_nlin(NID,IST,IERR)

    USE xxcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NID
    INTEGER,INTENT(OUT) :: IST,IERR
    INTEGER,PARAMETER:: NSM=100

    NAMELIST /XX/ nxmax,x0,dx,am

    READ(NID,XX,IOSTAT=IST,ERR=9800,END=9900)

    IERR=0
    RETURN

9800 IERR=8
    RETURN
9900 IERR=9
    RETURN
  END SUBROUTINE xx_nlin

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE xx_plst

    IMPLICIT NONE
    WRITE(6,'(A)') '# &XX : nxmax,x0,dx,am'
    RETURN

  END SUBROUTINE xx_plst

!     ****** CHECK INPUT PARAMETER ******

  SUBROUTINE xx_check(IERR)

    USE xxcomm_parm
    IMPLICIT NONE
    INTEGER:: IERR

    IERR=0

    IF(dx <= 0.D0) THEN
       WRITE(6,'(A,1PE12.4)') 'XX xx_check: INVALID dx: dx=',dx
       IERR=1
    ENDIF
    RETURN
  END SUBROUTINE xx_check

!     ****** SHOW PARAMETERS ******

  SUBROUTINE xx_view

    use xxcomm_parm
    implicit none

    WRITE(6,601) 'x0    ',x0    ,'dx    ',dx    , &
                 'am    ',am

    WRITE(6,602) 'nxmax ',nxmax
    RETURN

601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
            2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
602 FORMAT(' ',A6,'=',I7,4X  :2X,A6,'=',I7,4X  : &
            2X,A6,'=',I7,4X  :2X,A6,'=',I7)
  END SUBROUTINE xx_view

END MODULE xxparm
