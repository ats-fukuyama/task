! cvinit.f90

MODULE cvinit

  PRIVATE
  PUBLIC cv_init

CONTAINS

!     ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE cv_init

    USE cvcomm_parm
    IMPLICIT NONE

    knam_csv_in='WHO-COVID-19-global-data-200821.csv' ! input data
    knam_csv_global_cases='cv-global-cases.csv'
    knam_csv_global_deaths='cv-global-deaths.csv'
    knam_csv_region_cases='cv-region-cases.csv'
    knam_csv_region_deaths='cv-region-deaths.csv'
    knam_csv_select_cases='cv-select-cases.csv'
    knam_csv_select_deaths='cv-select-deaths.csv'
    knam_cv_select='cv-select.data' ! list of country_id for select
    knam_cv_population='cv-population.data' ! list of population and area
    ndate_step_global=1  ! ndate step of global data: 1: all, 7: a week step
    ndate_step_region=1  ! ndate step of region data: 1: all, 7: a week step
    ndate_step_select=1  ! ndata step of select data: 1: all, 7: a week step
    ndate_start=1        ! initial ndate: 1 for 2020-01-04 (Sat)
    
    RETURN
  END SUBROUTINE cv_init
END MODULE cvinit
