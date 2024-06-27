! ids_ec_launchers.f90

MODULE ids_ec_launchers

  USE wrcomm
  USE task_kinds
  USE task_constants
  
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
    call imas_open_env('ids' ,131024,1,idx,'public',trim(machine),'3')
    call ids_get(idx,'ec_launchers',ec_launchers_ids)
    call imas_close(idx)
    
    ! --- set number of beam ---

    if (associated(ec_launchers_ids%beam)) then
       nraymax = size(ec_launchers_ids%beam)
    else
       nraymax = 0
    endif

    ! --- allocate ec beam initial values ---

    IF(ALLOCATED(RPIN)) DEALLOCATE(RPIN)
    IF(ALLOCATED(ZPIN)) DEALLOCATE(ZPIN)
    IF(ALLOCATED(PHIIN)) DEALLOCATE(PHIIN)
    IF(ALLOCATED(ANGTIN)) DEALLOCATE(ANGTIN)
    IF(ALLOCATED(ANGPIN)) DEALLOCATE(ANGPIN)
    IF(ALLOCATED(RFIN)) DEALLOCATE(RFIN)
    IF(ALLOCATED(UUIN)) DEALLOCATE(UUIN)
    IF(ALLOCATED(MODEWIN)) DEALLOCATE(MODEWIN)

    ALLOCATE(RPIN(nraymax))
    ALLOCATE(ZPIN(nraymax))
    ALLOCATE(PHIIN(nraymax))
    ALLOCATE(ANGTIN(nraymax))
    ALLOCATE(ANGPIN(nraymax))
    ALLOCATE(RFIN(nraymax))
    ALLOCATE(UUIN(nraymax))
    ALLOCATE(MODEWIN(nraymax))

    nt=1
    DO nray=1,nraymax
       RPIN(nray)= ec_launchers_ids%beam(nray)%launching_position%r(nt)
       ZPIN(nray)= ec_launchers_ids%beam(nray)%launching_position%z(nt)
       PHIIN(nray)= ec_launchers_ids%beam(nray)%launching_position%phi(nt)
       ANGTIN(nray)= ec_launchers_ids%beam(nray)%steering_angle_tor(nt)
       ANGPIN(nray)= ec_launchers_ids%beam(nray)%steering_angle_pol(nt)
       RFIN(nray)= ec_launchers_ids%beam(nray)%frequency%data(nt)
       UUIN(nray)= ec_launchers_ids%beam(nray)%power_launched%data(nt)
!       MODEWIN(nray)= ec_launchers_ids%beam(nray)%mode%data(nt)
    END DO

    RETURN

  END SUBROUTINE get_ids_ec_launchers

  SUBROUTINE clean_ids_ec_launchers
    IMPLICIT NONE
    deallocate(RPIN,ZPIN,PHIIN)
    deallocate(ANGTIN,ANGPIN)
    deallocate(RFIN,UUIN)
    deallocate(MODEWIN)
  END SUBROUTINE clean_ids_ec_launchers
END MODULE ids_ec_launchers
