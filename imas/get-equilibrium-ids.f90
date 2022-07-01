! get-equilirium-ids.f90

MODULE get-equilirium-ids

  PRIVATE
  PUBLIC get_equilirium_ids

CONTAINS

  SUBROUTINE get_equilirium_ids

    USE eqcomm
    USE ids_schemas
    USE ids_routines
    IMPLICIT NONE

    ! --- set size of R and Z grid ---

    nrgmax = size(equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim1)
    nzgmax = size(equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim2)

    IF(ALLOCATED(rg)) DEALLOCATE(rg)
    IF(ALLOCATED(zg)) DEALLOCATE(zg)
    IF(ALLOCATED(br)) DEALLOCATE(br)
    IF(ALLOCATED(bt)) DEALLOCATE(bt)
    IF(ALLOCATED(bt)) DEALLOCATE(bt)
    ALLOCATE(RG(nrgmax))
    ALLOCATE(ZG(nzgmax))
    allocate(br(ni,nj))
    allocate(bt(ni,nj))
    allocate(bz(ni,nj))
    allocate(psi2d(ni,nj))
    Rarr    = equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim1
    Zarr    = equilibrium_ids%time_slice(1)%profiles_2d(1)%grid%dim2
    br      = equilibrium_ids%time_slice(1)%profiles_2d(1)%b_field_r
    bt      = equilibrium_ids%time_slice(1)%profiles_2d(1)%b_field_tor
    bz      = equilibrium_ids%time_slice(1)%profiles_2d(1)%b_field_z
    psi2d   = qequilibrium_ids%time_slice(1)%profiles_2d(1)%psi
    psiedge = equilibrium_ids%time_slice(1)%global_quantities%psi_boundary
    psiax   = equilibrium_ids%time_slice(1)%global_quantities%psi_axis
    b0      = -equilibrium_ids%time_slice(1)%global_quantities%magnetic_axis%b_field_tor
    iplasma = equilibrium_ids%time_slice(1)%global_quantities%ip
  END SUBROUTINE get_equilirium_ids

  SUBROUTINE clean_equilibrium_ids
    IMPLICIT NONE
    deallocate(psi)
    deallocate(Rarr,Zarr)
    deallocate(br)
    deallocate(bt)
    deallocate(bz)
    deallocate(psi2d)
  END SUBROUTINE clean_equilibrium_ids
END MODULE get
