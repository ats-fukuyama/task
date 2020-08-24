! cvview.f90

MODULE cvview

CONTAINS

  !     ****** SHOW PARAMETERS ******

  SUBROUTINE cv_view

    USE cvcomm_parm
    IMPLICIT NONE

    WRITE(6,'(A,A)') 'knam_csv_in        =',knam_csv_in
    WRITE(6,'(A,A)') 'knam_csv_out_cases =',knam_csv_out_cases
    WRITE(6,'(A,A)') 'knam_csv_out_deaths=',knam_csv_out_deaths
    RETURN
  END SUBROUTINE cv_view
END MODULE cvview
