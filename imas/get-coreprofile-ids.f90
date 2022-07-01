! get-coreprofile-ids.f90

MODULE get-coreprofile-ids

  PRIVATE
  PUBLIC get_coreprofile_ids
  PUBLIC clean_coreprofile_ids

CONTAINS

  SUBROUTINE get_coreprofile_ids

    USE eqcomm
    USE ids_schemas
    USE ids_routines
    IMPLICIT NONE

    ! --- set psi grid ---
    
    npsi = size(core_profiles_ids%profiles_1d(1)%grid%psi)
    
    IF(Allocated(te)) DEALLOCATE(te))
    IF(Allocated(ne)) DEALLOCATE(ne))
    IF(Allocated(zeff)) DEALLOCATE(zeff))

    allocate(te(npsi),ne(npsi),Zeff(npsi))
    te    = core_profiles_ids%profiles_1d(1)%electrons%temperature
    ne    = core_profiles_ids%profiles_1d(1)%electrons%density
    Zeff  = core_profiles_ids%profiles_1d(1)%zeff

    RETURN
  END SUBROUTINE get_coreprofile_ids

  SUBROUTINE clean_core_profile_ids


