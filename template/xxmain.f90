!
!               ############# TASK/XX #############
!
!                  Template for TASK components
!
!-----------------------------------------------------------------------

PROGRAM xx_main
  USE xxinit,ONLY: xx_init
  USE xxparm,ONLY: xx_parm
  USE xxmenu,ONLY: xx_menu

  IMPLICIT none
  INTEGER :: ierr

  CALL GSOPEN
  WRITE(6,*) '## TASK/XX 2015/05/16'
  OPEN(7,STATUS='SCRATCH')
  
  CALL xx_init
  CALL xx_parm(1,'xxparm.nl',ierr)
  IF(ierr /= 0 .AND. ierr /= 2) THEN
     WRITE(6,*) 'XX Error during reading the namelist file: xxparm.nl'
     WRITE(6,*) '     ierr = ',ierr
     STOP
  END IF

  CALL xx_menu

  CLOSE(7)
  CALL GSCLOS
  STOP
END PROGRAM xx_main
