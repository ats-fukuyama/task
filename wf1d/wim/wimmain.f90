!
!               ############# TASK/WIM #############
!
!             1D wave solver with integral formulation
!                    magnetized plasma
!
!-----------------------------------------------------------------------

PROGRAM wim_main
  USE wimcomm
  USE wiminit
  USE wimparm
  USE wimmenu

  IMPLICIT none
  INTEGER(ikind)  :: ierr

  CALL GSOPEN
  WRITE(6,*) '## TASK/WIM 2016/03/19'
  OPEN(7,STATUS='SCRATCH')
  
  CALL wim_init
  CALL wim_parm(1,'wimparm',ierr)
  IF(ierr /= 0 .AND. ierr /= 2) THEN
     WRITE(6,*) 'WIM Error during reading the namelist file: wimparm'
     WRITE(6,*) '     ierr = ',ierr
     STOP
  END IF

  CALL wim_menu

  CLOSE(7)
  CALL GSCLOS
  STOP
END PROGRAM wim_main
