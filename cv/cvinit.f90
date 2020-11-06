! cvinit.f90

MODULE cvinit

  PRIVATE
  PUBLIC cv_init

CONTAINS

!     ****** INITIALIZE INPUT PARAMETERS ******

  SUBROUTINE cv_init

    USE cvcomm_parm
    IMPLICIT NONE

    knam_csv_in='WHO-COVID-19-global-data.csv' ! input data
    knam_csv_out_global='cv-global.csv'
    knam_csv_out_region='cv-region.csv'
    knam_csv_out_select='cv-select.csv'
    knam_cv_select='cv-select.data' ! list of country_id for select
    knam_cv_population='cv-population.data' ! list of population and area
    ndate_step_global=1  ! ndate step of global data: 1: all, 7: week interval
    ndate_step_region=1  ! ndate step of region data: 1: all, 7: week interval
    ndate_step_select=7  ! ndata step of select data: 1: all, 7: week interval
    ndate_start_global=1 ! initial ndate of global data:
    ndate_start_region=1 ! initial ndate of region data:
    ndate_start_select=7 ! initial ndate of select data:
                         !         1 for 2020-01-04 (Sat), 7 for Friday 
    ndays_ave=7          ! range of day averageing
    nrank_max=12         ! number of items in ranking lists
    ndate_min_g=1        ! minimum ndate for graphics
    ndate_max_g=1        ! maximum ndate for graphics (be updated after load)
    population_min_rank=10.D0 ! minimum population per million in ranking lists

    cases_number_log_min=10.D0    ! minimum number for ncases in log
    deaths_number_log_min=1.D0    ! minimum number for ndeaths in log
    cases_rate_log_min=0.1D0      ! minimum rate for ncases in log
    deaths_rate_log_min=0.01D0    ! minimum rate for ndeaths in log
    ratio_new_total_log_min=0.1D0 ! ratio of log minimum between new and total

    
    RETURN
  END SUBROUTINE cv_init
END MODULE cvinit
