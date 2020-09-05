! cvview.f90

MODULE cvview

CONTAINS

  !     ****** SHOW PARAMETERS ******

  SUBROUTINE cv_view

    USE cvcomm_parm
    IMPLICIT NONE

    WRITE(6,'(A,A)') 'knam_csv_in           =',knam_csv_in
    WRITE(6,'(A,A)') 'knam_csv_global_cases =',knam_csv_global_cases
    WRITE(6,'(A,A)') 'knam_csv_global_deaths=',knam_csv_global_deaths
    WRITE(6,'(A,A)') 'knam_csv_region_cases =',knam_csv_region_cases
    WRITE(6,'(A,A)') 'knam_csv_region_deaths=',knam_csv_region_deaths
    WRITE(6,'(A,A)') 'knam_csv_select_cases =',knam_csv_select_cases
    WRITE(6,'(A,A)') 'knam_csv_select_deaths=',knam_csv_select_deaths
    WRITE(6,'(A,A)') 'knam_cv_select        =',knam_cv_select
    WRITE(6,'(A,A)') 'knam_cv_population    =',knam_cv_population
    WRITE(6,'(A,I8)') 'ndate_step_global          =',ndate_step_global
    WRITE(6,'(A,I8)') 'ndate_step_region          =',ndate_step_region
    WRITE(6,'(A,I8)') 'ndate_step_select          =',ndate_step_select
    WRITE(6,'(A,I8)') 'ndate_start                =',ndate_start
    RETURN
  END SUBROUTINE cv_view
END MODULE cvview
