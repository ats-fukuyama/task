!     $Id$
!
!               ############# TASK/XX #############
!
!                  Template for TASK components
!
!-----------------------------------------------------------------------

PROGRAM xx_main
  USE xxcomm
  USE xxinit
  USE xxparm
  USE xxmenu

  IMPLICIT none
  INTEGER(ikind)  :: ierr

  CALL GSOPEN
  WRITE(6,*) '## TASK/XX 2013/11/15'
  OPEN(7,STATUS='SCRATCH')
  
  CALL xx_init
  CALL xx_parm(1,'xxparm',ierr)
  IF(ierr /= 0 .AND. ierr /= 2) THEN
     WRITE(6,*) 'XX Error during reading the namelist file: xxparm'
     WRITE(6,*) '     ierr = ',ierr
     STOP
  END IF

  CALL xx_menu

  CLOSE(7)
  CALL GSCLOS
  STOP
END PROGRAM xx_main
