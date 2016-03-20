!    $Id$

Module xxgout
  PRIVATE
  PUBLIC xx_gout
 
CONTAINS

  SUBROUTINE xx_gout

    USE xxcomm,ONLY: rkind,nxmax,x,y
    USE libgrf
    IMPLICIT NONE

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') &
         '#### XX GOUT: A X/exit'
    CALL TASK_KLIN(line,kid,mode,xx_parm)
    IF(mode == 2 .OR. mode == 3) GOTO 1

    ICH=ICHAR(LINE(1:1))
    IF(ICH.GE.97.AND.ICH.LE.122) ICH=ICH-32
    KCH=CHAR(ICH)

    SELECT CASE(kch)
    CASE('A') 
       CALL PAGES
       CALL GRD1D(0,x,y,nxmax,nxmax,1,'xx')
       CALL PAGEE
    CASE('X') GO TO 9000
    END SELECT
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE xx_gout
END Module xxgout
