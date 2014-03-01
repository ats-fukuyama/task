!     $Id$
!
!               ############# TASK/WI #############
!
!                  Template for TASK components
!
!-----------------------------------------------------------------------

PROGRAM wi_main
  USE wicomm
  USE wiparm
  USE wimenu

  IMPLICIT none
  INTEGER(ikind)  :: ierr

  CALL GSOPEN
  WRITE(6,'(A)') '## TASK/WI 2014/03/01'
  OPEN(7,STATUS='SCRATCH')
  
  CALL wi_parm(1,'wiparm.nl',ierr)
  IF(ierr /= 0 .AND. ierr /= 2) THEN
     WRITE(6,*) 'XX Error during reading the namelist file: wiparm.nl'
     WRITE(6,*) '     ierr = ',ierr
     STOP
  END IF

  CALL wi_menu

  CLOSE(7)
  CALL GSCLOS
  STOP
END PROGRAM wi_main
