! ptmenu.f90

! *** pt menu ***

MODULE ptmenu

  PRIVATE
  PUBLIC pt_menu

CONTAINS

  SUBROUTINE pt_menu

    USE ptcomm
    USE ptinit
    USE ptparm
    USE ptview
    USE ptexec1
    USE libmpi
      
    IMPLICIT NONE
    INTEGER(ikind) :: IERR, MODE
    INTEGER(ikind), SAVE :: INIT=0
    CHARACTER(LEN=1) :: KID
    CHARACTER(LEN=80):: LINE

!     ------ SELECTION OF TASK TYPE ------

    IERR=0

1   CONTINUE
    WRITE(6,601)
601 FORMAT('## pt menu: P,V/parm  R/exec  Q/quit')
    CALL TASK_KLIN(LINE,KID,MODE,pt_parm)

    CALL mtx_broadcast1_character(KID)
    CALL mtx_broadcast1_integer(MODE)

    IF(MODE.NE.1) GOTO 1

    SELECT CASE(KID)
    CASE('P')                     ! parameter input
       CALL pt_parm(0,'PT',IERR)
    CASE('V')                     ! view input parameters
       CALL pt_view
    CASE('R')                     ! load bulk component profile from TR
       CALL pt_exec1
    CASE('Q')                     ! quit 
       GOTO 9000
    CASE('X','#','!')
       CONTINUE
    CASE default
       WRITE(6,*) 'XX ptmenu: UNKNOWN KID'
    END SELECT
    
    GOTO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE pt_menu
END MODULE ptmenu
