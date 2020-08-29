! cvinit.f90

MODULE cvinit

  PRIVATE
  PUBLIC cv_init

CONTAINS

!     ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE cv_init

    USE cvcomm_parm
    IMPLICIT NONE

    knam_csv_in='WHO-COVID-19-global-data-200821.csv'
    knam_csv_out_cases='cv-data-cases.csv'
    knam_csv_out_deaths='cv-data-deaths.csv'
    knam_csv_region_cases='cv-region-cases.csv'
    knam_csv_region_deaths='cv-region-deaths.csv'
    knam_csv_select_cases='cv-select-cases.csv'
    knam_csv_select_deaths='cv-select-deaths.csv'
    knam_cv_select='cv-select.data'
    
    RETURN
  END SUBROUTINE cv_init
END MODULE cvinit
