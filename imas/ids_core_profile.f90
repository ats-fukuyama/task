! ids_core_profiles.f90

MODULE ids_core_profiles

  USE task_kinds
  USE task_constants
  INTEGER:: nrmax,nsmax
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: &
       rhotn
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: &
       pa,pz
  REAL(rkind),ALLOCATABLE,DIMENSION(:,:):: &
       rn,ru,rtpr,rtpp
  REAL(rkind),ALLOCATABLE,DIMENSION(:):: &
       zeff
  
  PRIVATE
  PUBLIC get_ids_core_profiles
  PUBLIC clear_ids_core_profiles

CONTAINS

  SUBROUTINE get_ids_core_profiles

    USE ids_schemas,ONLY: ids_core_profiles
    USE ids_routines,ONLY: &
         imas_open_env,imas_create_env,imas_close,ids_get
    IMPLICIT NONE
    type(ids_core_profiles):: core_profiles_ids
    character(len=200):: user,machine
    integer:: idx,err_flag
    character(len=:), pointer:: err_msg
    INTEGER:: ns,nr

    ! DEFINE LOCAL DATABASE
    call getenv('USER',user)
    machine = 'ITER'

    ! OPEN INPUT DATAFILE FROM OFFICIAL ITER DB
    call imas_open_env('ids' ,131024,1,idx,'public',trim(machine),'3')
    call ids_get(idx,'core_profiles',core_profiles_ids)
    call imas_close(idx)
    
    ! --- set size of R and Z grid ---

    nrmax = size(core_profiles_ids%profiles_1d(1)%grid%rho_tor_norm)
    nsmax = size(core_profiles_ids%profiles_1d(1)%ion)+1

    IF(ALLOCATED(rhotn)) DEALLOCATE(rhotn)
    IF(ALLOCATED(pa)) DEALLOCATE(pa)
    IF(ALLOCATED(pz)) DEALLOCATE(pz)
    IF(ALLOCATED(rn)) DEALLOCATE(rn)
    IF(ALLOCATED(ru)) DEALLOCATE(ru)
    IF(ALLOCATED(rtpr)) DEALLOCATE(rtpr)
    IF(ALLOCATED(rtpp)) DEALLOCATE(rtpp)
    IF(ALLOCATED(zeff)) DEALLOCATE(zeff)

    ALLOCATE(rhotn(nrmax))
    ALLOCATE(pa(nsmax))
    ALLOCATE(pz(nsmax))
    ALLOCATE(rn(nrmax,nsmax))
    ALLOCATE(ru(nrmax,nsmax))
    ALLOCATE(rtpr(nrmax,nsmax))
    ALLOCATE(rtpp(nrmax,nsmax))
    ALLOCATE(zeff(nrmax))

    rhotn = core_profiles_ids%profiles_1d(1)%grid%rho_tor_norm
    pa(1)=AME/AMP
    pz(1)=-1.D0
    rn(1:nrmax,1)= core_profiles_ids%profiles_1d(1)%electrons%density(1:nrmax)
    ru(1:nrmax,1)= 0.D0
    rtpr(1:nrmax,1) &
         = core_profiles_ids%profiles_1d(1)%electrons%temperature(1:nrmax)
    rtpp(1:nrmax,1) &
         = core_profiles_ids%profiles_1d(1)%electrons%temperature(1:nrmax)
    DO ns=2,nsmax
       pa(ns)=core_profiles_ids%profiles_1d(1)%ion(ns)%element(1)%a
       pz(ns)=core_profiles_ids%profiles_1d(1)%ion(ns)%element(1)%z_n
       rn(1:nrmax,ns) &
            = core_profiles_ids%profiles_1d(1)%ion(ns)%density(1:nrmax)
       ru(1:nrmax,ns) &
          = core_profiles_ids%profiles_1d(1)%ion(ns)%velocity%toroidal(1:nrmax)
       rtpr(1:nrmax,ns) &
            = core_profiles_ids%profiles_1d(1)%ion(ns)%temperature(1:nrmax)
       rtpp(1:nrmax,1) &
            = core_profiles_ids%profiles_1d(1)%ion(ns)%temperature(1:nrmax)
       Zeff  = core_profiles_ids%profiles_1d(1)%zeff
    END DO
    RETURN

  END SUBROUTINE get_ids_core_profiles

  SUBROUTINE clean_ids_core_profiles
    IMPLICIT NONE
    deallocate(rhotn)
    deallocate(pa,pz)
    deallocate(rn,ru,rtpr,rtpp)
    deallocate(zeff)
  END SUBROUTINE clean_ids_core_profiles
END MODULE ids_core_profiles
