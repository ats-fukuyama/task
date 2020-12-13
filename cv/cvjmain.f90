! cvjmain.f90

!               ############# TASK/CV #############

!                      COVID-19 data processing for japan

!                           A. Fukuyama

!                      V1.00  : 2020 DECEMBER 06
!     --------------------------------------------------------

PROGRAM cvj

  USE cvcomm
  USE cvinit,ONLY: cv_init
  USE cvparm,ONLY: cv_parm
  USE cvjmenu,ONLY: cvj_menu
  USE cvlib
  IMPLICIT NONE
  INTEGER:: ierr
  EXTERNAL GSOPEN,GSCLOS

  WRITE(6,*) '## TASK/CVJ 2020/12/06'
  CALL GSOPEN

  CALL cv_init
  OPEN(7,STATUS='SCRATCH')
  CALL cv_parm(1,'cvjparm',IERR)

  CALL cvj_menu

  CLOSE(7)
  CALL GSCLOS

  STOP
END PROGRAM cvj
