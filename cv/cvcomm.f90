! cvcomm.f90

MODULE cvcomm_parm
  CHARACTER(LEN=256):: knam_csv_in,knam_csv_out_cases,knam_csv_out_deaths

END MODULE cvcomm_parm

MODULE cvcomm
  USE cvcomm_parm

  INTEGER:: nlist_max,ndate_max
  CHARACTER(LEN=2),ALLOCATABLE:: country_id_nlist(:)
  CHARACTER(LEN=80),ALLOCATABLE:: country_name_nlist(:)
  CHARACTER(LEN=4),ALLOCATABLE:: region_name_nlist(:)
  CHARACTER(LEN=10),ALLOCATABLE:: date_id_ndate(:)
  INTEGER,DIMENSION(:,:),ALLOCATABLE:: &
       ncases_new_ndate_nlist,ncases_total_ndate_nlist, &
       ndeaths_new_ndate_nlist,ndeaths_total_ndate_nlist
END MODULE cvcomm
