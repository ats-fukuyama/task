! ids_equilibrium.f90

MODULE ids_equilibrium

  USE task_kinds
  INTEGER:: nrgmax,nzgmax
  REAL(rkind):: b0,rip,psip0,psipa
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: rg,zg
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: br,bt,bz,psip2d
  
  PRIVATE
  PUBLIC get_ids_equilibrium
  PUBLIC clear_ids_equilibrium

CONTAINS

  SUBROUTINE get_ids_equilibrium

    USE ids_schemas,ONLY: ids_equilibrium
    USE ids_routines,ONLY: &
         imas_open_env,imas_create_env,imas_close,ids_get
    IMPLICIT NONE
    type(ids_equilibrium):: equilibrium_ids
    character(len=200):: user,machine
    integer:: idx,err_flag
    character(len=:), pointer:: err_msg
    INTEGER:: nrgmax,nzgmax

    ! DEFINE LOCAL DATABASE
    call getenv('USER',user)
    machine = 'ITER'

    ! OPEN INPUT DATAFILE FROM OFFICIAL ITER DB
    call imas_open_env('ids' ,131024,1,idx, &
         'public',trim(machine),'3')
    call ids_get(idx,'equilibrium',equilibrium_ids)
    call imas_close(idx)
    
    ! --- set size of R and Z grid ---

    nrgmax = size(equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim1)
    nzgmax = size(equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim2)

    IF(ALLOCATED(rg)) DEALLOCATE(rg)
    IF(ALLOCATED(zg)) DEALLOCATE(zg)
    IF(ALLOCATED(br)) DEALLOCATE(br)
    IF(ALLOCATED(bt)) DEALLOCATE(bt)
    IF(ALLOCATED(bt)) DEALLOCATE(bt)
    IF(ALLOCATED(psip2d)) DEALLOCATE(psip2d)
    ALLOCATE(RG(nrgmax))
    ALLOCATE(ZG(nzgmax))
    allocate(br(nrgmax,nzgmax))
    allocate(bt(nrgmax,nzgmax))
    allocate(bz(nrgmax,nzgmax))
    allocate(psip2d(nrgmax,nzgmax))
    rg    = equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim1
    zg    = equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim2
    br    = equilibrium_ids%time_slice(1)%profiles_2d(1)%b_field_r
    bt    = equilibrium_ids%time_slice(1)%profiles_2d(1)%b_field_tor
    bz    = equilibrium_ids%time_slice(1)%profiles_2d(1)%b_field_z
    psip2d= equilibrium_ids%time_slice(1)%profiles_2d(1)%psi
    psipa = equilibrium_ids%time_slice(1)%global_quantities%psi_boundary
    psip0 = equilibrium_ids%time_slice(1)%global_quantities%psi_axis
    b0    = -equilibrium_ids%time_slice(1)%global_quantities%magnetic_axis%b_field_tor
    rip   = equilibrium_ids%time_slice(1)%global_quantities%ip
  END SUBROUTINE get_ids_equilibrium

  SUBROUTINE clean_ids_equilibrium
    IMPLICIT NONE
    deallocate(psip2d)
    deallocate(rg,zg)
    deallocate(br)
    deallocate(bt)
    deallocate(bz)
  END SUBROUTINE clean_ids_equilibrium
END MODULE ids_equilibrium
