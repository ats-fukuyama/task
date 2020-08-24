! cvmain.f90

!               ############# TASK/CV #############

!                      COVID-19 data processing

!                           A. Fukuyama

!                      V1.00  : 2020 AUGUST 23
!     --------------------------------------------------------

PROGRAM cv

  USE cvcomm
  USE cvinit,ONLY: cv_init
  USE cvparm,ONLY: cv_parm
  USE cvmenu,ONLY: cv_menu
  USE cvlib
  IMPLICIT NONE
  INTEGER:: ierr

  WRITE(6,*) '## TASK/CV 2020/08/23'
  CALL GSOPEN

  CALL cv_init
  OPEN(7,STATUS='SCRATCH')
  CALL cv_parm(1,'cvparm',IERR)

  CALL cv_menu

  CLOSE(7)
  CALL GSCLOS

  STOP
END PROGRAM cv
