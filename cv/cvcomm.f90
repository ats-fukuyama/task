! cvcomm.f90

MODULE cvcomm_parm
  INTEGER,PARAMETER :: dp = selected_real_kind(15) !double precision
  CHARACTER(LEN=256):: knam_csv_in,knam_csv_out_cases,knam_csv_out_deaths
  CHARACTER(LEN=256):: knam_csv_region_cases,knam_csv_region_deaths
  CHARACTER(LEN=256):: knam_csv_select_cases,knam_csv_select_deaths
  CHARACTER(LEN=256):: knam_cv_select
  CHARACTER(LEN=256):: knam_cv_population

END MODULE cvcomm_parm

MODULE cvcomm
  USE cvcomm_parm

  INTEGER:: ncountry_max,ndate_max,nregion_max,nselect_max
  CHARACTER(LEN=2),ALLOCATABLE:: country_id_ncountry(:)
  CHARACTER(LEN=80),ALLOCATABLE:: country_name_ncountry(:)
  CHARACTER(LEN=4),ALLOCATABLE:: region_id_ncountry(:)
  CHARACTER(LEN=10),ALLOCATABLE:: date_id_ndate(:)
  CHARACTER(LEN=4),ALLOCATABLE:: region_id_nregion(:)
  CHARACTER(LEN=80),ALLOCATABLE:: region_name_nregion(:)
  CHARACTER(LEN=2),ALLOCATABLE:: country_id_nselect(:)
  CHARACTER(LEN=80),ALLOCATABLE:: country_name_nselect(:)
  INTEGER,DIMENSION(:,:),ALLOCATABLE:: &
       ncases_new_ndate_ncountry,ncases_total_ndate_ncountry, &
       ndeaths_new_ndate_ncountry,ndeaths_total_ndate_ncountry
  INTEGER,DIMENSION(:,:),ALLOCATABLE:: &
       ncases_new_ndate_nregion,ncases_total_ndate_nregion, &
       ndeaths_new_ndate_nregion,ndeaths_total_ndate_nregion
  INTEGER,DIMENSION(:,:),ALLOCATABLE:: &
       ncases_new_ndate_nselect,ncases_total_ndate_nselect, &
       ndeaths_new_ndate_nselect,ndeaths_total_ndate_nselect
  INTEGER,DIMENSION(:),ALLOCATABLE:: &
       nregion_ncountry,ncount_nregion
  INTEGER,DIMENSION(:),ALLOCATABLE:: &
       ncountry_nselect
  REAL(dp),DIMENSION(:),ALLOCATABLE:: &
       population_ncountry(:),area_ncountry(:)
END MODULE cvcomm
