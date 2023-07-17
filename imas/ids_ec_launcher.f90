! ids_ec_launcher.f90

MODULE ids_ec_launchers

  USE task_kinds
  USE task_constants
  INTEGER:: nbeam
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: &
       r_lc,z_lc,phi_lc
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: &
       ang_t_lc,ang_p_lc
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: &
       freq_lc,power_lc
  INTEGER,ALLOCATABLE,DIMENSION(:):: &
       mode_lc
  
  PRIVATE
  PUBLIC get_ids_ec_launchers
  PUBLIC clear_ids_ec_launchers

CONTAINS

  SUBROUTINE get_ids_ec_launchers

    USE ids_schemas,ONLY: ids_ec_launchers
    USE ids_routines,ONLY: &
         imas_open_env,imas_create_env,imas_close,ids_get
    IMPLICIT NONE
    type(ids_ec_launchers):: ec_launchers_ids
    character(len=200):: user,machine
    integer:: idx,err_flag
    character(len=:), pointer:: err_msg
    INTEGER:: nb,nt

    ! DEFINE LOCAL DATABASE
    call getenv('USER',user)
    machine = 'ITER'

    ! OPEN INPUT DATAFILE FROM OFFICIAL ITER DB
    call imas_open_env('ids' ,131024,1,idx, &
         'public',trim(machine),'3')
    call ids_get(idx,'ec_launchers',ec_launchers_ids)
    call imas_close(idx)
    
    ! --- set number of beam ---

    if (associated(ec_launchers_ids%launcher)) then
       nbeam = size(ec_launchers_ids%launcher)
    else
       nbeam = 0
    endif

    ! --- allocate ec beam initial values ---

    IF(ALLOCATED(r_lc)) DEALLOCATE(r_lc)
    IF(ALLOCATED(z_lc)) DEALLOCATE(z_lc)
    IF(ALLOCATED(phi_lc)) DEALLOCATE(phi_lc)
    IF(ALLOCATED(ang_t_lc)) DEALLOCATE(ang_t_lc)
    IF(ALLOCATED(ang_p_lc)) DEALLOCATE(ang_p_lc)
    IF(ALLOCATED(freq_lc)) DEALLOCATE(freq_lc)
    IF(ALLOCATED(power_lc)) DEALLOCATE(power_lc)
    IF(ALLOCATED(mode_lc)) DEALLOCATE(mode_lc)

    ALLOCATE(r_lc(nbeam))
    ALLOCATE(z_lc(nbeam))
    ALLOCATE(phi_lc(nbeam))
    ALLOCATE(ang_t_lc(nbeam))
    ALLOCATE(ang_p_lc(nbeam))
    ALLOCATE(freq_lc(nbeam))
    ALLOCATE(power_lc(nbeam))
    ALLOCATE(mode_lc(nbeam))

    nt=1
    DO nb=1,nbeam
       r_lc(nb)= ec_launchers_ids%launcher(nb)%launching_position%r(nt)
       z_lc(nb)= ec_launchers_ids%launcher(nb)%launching_position%z(nt)
       phi_lc(nb)= ec_launchers_ids%launcher(nb)%launching_position%phi(nt)
       ang_t_lc(nb)= ec_launchers_ids%launcher(nb)%steering_angle_tor%data(nt)
       ang_p_lc(nb)= ec_launchers_ids%launcher(nb)%steering_angle_pol%data(nt)
       freq_lc(nb)= ec_launchers_ids%launcher(nb)%frequency%data(nt)
       power_lc(nb)= ec_launchers_ids%launcher(nb)%power_launched%data(nt)
       mode_lc(nb)= ec_launchers_ids%launcher(nb)%mode%data(nt)
    END DO

    RETURN

  END SUBROUTINE get_ids_ec_launchers

  SUBROUTINE clean_ids_ec_launchers
    IMPLICIT NONE
    deallocate(r_lc,z_lc,phi_lc)
    deallocate(ang_t_lc,ang_p_lc)
    deallocate(freq_lc,power_lc)
    deallocate(mode_lc)
  END SUBROUTINE clean_ids_ec_launchers
END MODULE ids_ec_launchers
