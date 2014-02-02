!     $Id$

Module wiparm

  PRIVATE

  PUBLIC wi_parm, wi_view

CONTAINS

!     ****** PARAMETER INPUT ******

  SUBROUTINE wi_parm(MODE,KIN,IERR)

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

    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode
    CHARACTER(LEN=*),INTENT(IN)::  kin
    INTEGER,INTENT(OUT):: ierr

1   CALL TASK_PARM(MODE,'WI',KIN,wi_nlin,wi_plst,IERR)
    IF(IERR.NE.0 .AND. IERR.NE.2) RETURN

    CALl wi_check(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
    IF(IERR.NE.0) IERR=IERR+100

    RETURN
  END SUBROUTINE wi_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE wi_nlin(NID,IST,IERR)

    USE wicomm, ONLY: nxmax,nwmax,modewi,xmax,pn0,alfa,aky,beta,cfyn, &
                      ntaumax,taumin,taumax

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NID
    INTEGER,INTENT(OUT) :: IST,IERR

    NAMELIST /WI/ nxmax,nwmax,modewi,xmax,pn0,alfa,aky,beta,cfyn, &
                  ntaumax,taumin,taumax

    READ(NID,WI,IOSTAT=IST,ERR=9800,END=9900)

    IERR=0
    RETURN

9800 IERR=8
    RETURN
9900 IERR=9
    RETURN
  END SUBROUTINE wi_nlin

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE wi_plst

    IMPLICIT NONE
    WRITE(6,'(A)') '# &WI : nxmax,nwmax,modewi,xmax,pn0,alfa,aky,beta,cfyn,'
    WRITE(6,'(A)') '        ntaumax,taumin,taumax'
    RETURN

  END SUBROUTINE wi_plst

!     ****** CHECK INPUT PARAMETER ******

  SUBROUTINE wi_check(IERR)

    USE wicomm,ONLY: xmax,taumin,taumax
    IMPLICIT NONE
    INTEGER:: IERR

    IERR=0

    IF(xmax <= 0.D0) THEN
       WRITE(6,'(A,1PE12.4)') 'XX wi_check: INVALID xmax: xmax=',xmax
       IERR=1
    ENDIF
    IF(taumax <= taumin) THEN
       WRITE(6,'(A,A,1P2E12.4)') 'XX wi_check: INVALID taumin,taumax: ',&
                                 'taumin,taumax =',taumin,taumax
       IERR=1
    ENDIF
    RETURN
  END SUBROUTINE wi_check

!     ****** SHOW PARAMETERS ******

  SUBROUTINE wi_view

    use wicomm
    implicit none

    WRITE(6,602) 'nxmax   ',nxmax, 'nwmax   ',nwmax , &
                 'modewi  ',modewi,'ntaumax ',ntaumax

    WRITE(6,601) 'xmax  ',xmax  ,'pn0   ',pn0   , &
                 'alfa  ',alfa  ,'aky   ',aky
    WRITE(6,601) 'beta  ',beta  ,'taumin',taumin, &
                 'taumax',taumax

    WRITE(6,603) 'cfyn  ',cfyn
    RETURN

601 FORMAT(' ',A6,'=',1PE11.3:2X,A6,'=',1PE11.3: &
            2X,A6,'=',1PE11.3:2X,A6,'=',1PE11.3)
602 FORMAT(' ',A8,'=',I7,2X  :2X,A8,'=',I7,2X  : &
            2X,A8,'=',I7,2X  :2X,A8,'=',I7)
603 FORMAT(' ',A6,'=',1P2E11.3:2X,A6,'=',1P2E11.3)

  END SUBROUTINE wi_view

END MODULE wiparm
  
