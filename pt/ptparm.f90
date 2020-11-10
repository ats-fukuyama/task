! ptparm.f90

MODULE ptparm

  PRIVATE
  PUBLIC pt_parm

CONTAINS

!     ****** INPUT PARAMETERS ******

  SUBROUTINE pt_parm(mode,kin,ierr)

!     mode=0 : standard namelinst input
!     mode=1 : namelist file input
!     mode=2 : namelist line input

!     ierr=0 : normal end
!     ierr=1 : namelist standard input error
!     ierr=2 : namelist file does not exist
!     ierr=3 : namelist file open error
!     ierr=4 : namelist file read error
!     ierr=5 : namelist file abormal end of file
!     ierr=6 : namelist line input error
!     ierr=7 : unknown MODE
!     ierr=10X : input parameter out of range

    IMPLICIT NONE
    INTEGER,INTENT(IN):: mode
    CHARACTER(LEN=*),INTENT(IN):: kin
    INTEGER,INTENT(OUT):: ierr

    ierr=0

1   CALL task_parm(mode,'PT',kin,pt_nlin,pt_plst,ierr)
    IF(ierr.NE.0) RETURN

    CALl pt_chek(ierr)
    IF(mode.EQ.0.AND.ierr.NE.0) GOTO 1

    RETURN
  END SUBROUTINE pt_parm

!     ****** INPUT NAMELIST ******

  SUBROUTINE pt_nlin(nid,ist,ierr)

    USE ptcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nid
    INTEGER,INTENT(OUT):: ist,ierr
    INTEGER:: NS

    NAMELIST /pt/ &
         bb,rr,ra,rkap,rdlt,ngxmax,ngymax,nthmax

    READ(nid,pt,IOSTAT=ist,ERR=9800,END=9900)
    
    ierr=0
    RETURN

9800 CONTINUE
    ierr=8
    RETURN
9900 CONTINUE
    ierr=9
    RETURN
  END SUBROUTINE pt_nlin

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE pt_plst

    WRITE(6,'(A)') &
         '# &pt : bb,rr,ra,rkap,rdlt,ngxmax,ngymax,nthmax'
    RETURN
  END SUBROUTINE pt_plst

!     ***** CHECK INPUT PARAMETERS *****

  SUBROUTINE pt_chek(ierr)

    USE ptcomm_parm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr

    ierr=0

    RETURN
  END SUBROUTINE pt_chek
END MODULE ptparm

