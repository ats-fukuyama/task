! ids_ec_launchers.f90

MODULE ids_ec_launchers

  USE task_kinds
  USE task_constants

  INTEGER:: nbeam_max,ntmax_rf,ntmax_pw,ntmax_ec
  CHARACTER(LEN=256),ALLOCATABLE:: &
       BNAME_L(:),BIDNAME_L(:)
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       RPINMIN_L,RPINMAX_L
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       RFIN_L,RFIN_T,PECIN_L,PECIN_T,TIME_L,RPIN_L,ZPIN_L,PHIIN_L, &
       ANGTIN_L,ANGPIN_L,OXFRAC_L,PCURV1_L,PCURV2_L,PANGLE_L, &
       BSIZE1_L,BSIZE2_L,BANGLE_L

  PRIVATE
  PUBLIC get_ids_ec_launchers
  PUBLIC clear_ids_ec_launchers

CONTAINS

  SUBROUTINE get_ids_ec_launchers(ierr)

    USE wrcomm
    USE ids_schemas,ONLY: ids_ec_launchers
    USE ids_routines,ONLY: &
         imas_open_env,imas_create_env,imas_close,ids_get
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    type(ids_ec_launchers):: ec_launchers_ids
    character(len=200):: user,machine
    INTEGER:: nbeam,nt_rf,nt_pw,nt_ec
    integer:: idx,err_flag
    character(len=:), pointer:: err_msg

    !     beta  = -ec_launchers_ids%launcher(ibeam)%steering_angle_tor%data(1)
    !     alpha =  ec_launchers_ids%launcher(ibeam)%steering_angle_pol%data(1)
    !     xpoldeg = asin(cos(beta)*sin(alpha))*180./pi
    !     xtordeg = atan(tan(beta)/cos(alpha))*180./pi

    ierr=0

    IF(ALLOCATED(BNAME_L)) CALL clear_ids_ec_launchers
    
    ! DEFINE LOCAL DATABASE
    call getenv('USER',user)
    machine = 'ITER'

    call imas_open_env('ids' ,131024,1,idx,'public',trim(machine),'3')
    call ids_get(idx,'ec_launchers',ec_launchers_ids)
    call imas_close(idx)
    
    ! --- set number of beam ---

    if (associated(ec_launchers_ids%beam)) then
       nbeam_max = size(ec_launchers_ids%beam)
    else
       nbeam_max = 0
    endif

    IF(nbeam_max.EQ.0) THEN
       WRITE(6,*) 'XX nbeam_max.EQ.0: no EC launcher data' 
       ierr=1001
       RETURN
    END IF

    ntmax_rf=0
    ntmax_pw=0
    ntmax_ec=0
    DO nbeam=1,nbeam_max
       ntmax_rf=MAX(ntmax_rf, &
                    SIZE(ec_launchers_ids%beam(nbeam)%frequency%data))
       ntmax_pw=MAX(ntmax_pw, &
                    SIZE(ec_launchers_ids%beam(nbeam)%power_launched%data))
       ntmax_ec=MAX(ntmax_ec, &
                    SIZE(ec_launchers_ids%beam(nbeam)%time))
    END DO

    ! --- allocate ec beam initial values ---

    
    ALLOCATE(BNAME_L(nbeam_max))    ! beam name
    ALLOCATE(BIDNAME_L(nbeam_max))  ! beam identifier
    ALLOCATE(RPINMIN_L(nbeam_max))          ! position: min of major radius [m]
    ALLOCATE(RPINMAX_L(nbeam_max))          ! position: max of major radius [m]
    ALLOCATE(RFIN_L(ntmax_rf,nbeam_max))    ! frequency [Hz]
    ALLOCATE(RFIN_T(ntmax_rf,nbeam_max))    ! timing of frequency [s]
    ALLOCATE(PECIN_L(ntmax_pw,nbeam_max))   ! input power [W]
    ALLOCATE(PECIN_T(ntmax_pw,nbeam_max))   ! timing of input power [W]
    ALLOCATE(TIME_L(ntmax_ec,nbeam_max))   ! prescribed time [s]
    ALLOCATE(RPIN_L(ntmax_ec,nbeam_max))   ! position: major radius [m]
    ALLOCATE(ZPIN_L(ntmax_ec,nbeam_max))   ! position: height [m]
    ALLOCATE(PHIIN_L(ntmax_ec,nbeam_max))  ! position: toroidal angle [rad]
    ALLOCATE(ANGTIN_L(ntmax_ec,nbeam_max)) ! toroidal ang: Asin(k_phi,k) [rad]
    ALLOCATE(ANGPIN_L(ntmax_ec,nbeam_max)) ! poloidal ang: Atan2(-kz,-kR)[rad]
    ALLOCATE(OXFRAC_L(ntmax_ec,nbeam_max)) ! O-mode power fraction: 0 for X
    ALLOCATE(PCURV1_L(ntmax_ec,nbeam_max)) ! phase horizontal curv in [m-1]
    ALLOCATE(PCURV2_L(ntmax_ec,nbeam_max)) ! phase vertical curv in [m-1]
    ALLOCATE(PANGLE_L(ntmax_ec,nbeam_max)) ! rotation angle of ph curv [rad] 
    ALLOCATE(BSIZE1_L(ntmax_ec,nbeam_max)) ! horizontal size of spot [m]
    ALLOCATE(BSIZE2_L(ntmax_ec,nbeam_max)) ! vertical size of spot [m]
    ALLOCATE(BANGLE_L(ntmax_ec,nbeam_max)) ! rotation angle of spot rad]

    DO nbeam=1,nbeam_max
       BNAME_L(nbeam)=ec_launchers_ids%beam(nbeam)%name(1)
       BIDNAME_L(nbeam)=ec_launchers_ids%beam(nbeam)%identifier(1)
       DO nt_rf=1,ntmax_rf
          RFIN_L(nt_rf,nbeam)= ec_launchers_ids%beam(nbeam)% &
               frequency%data(nt_rf)
          RFIN_T(nt_rf,nbeam)= ec_launchers_ids%beam(nbeam)% &
               frequency%time(nt_rf)
       END DO
       DO nt_pw=1,ntmax_pw
          PECIN_L(nt_pw,nbeam)= ec_launchers_ids%beam(nbeam)% &
               power_launched%data(nt_pw)
          PECIN_L(nt_pw,nbeam)= ec_launchers_ids%beam(nbeam)% &
               power_launched%time(nt_pw)
       END DO
       DO nt_ec=1,ntmax_ec
          TIME_L(nt_ec,nbeam)=ec_launchers_ids%beam(nbeam)%time(nt_ec)
          RPIN_L(nt_ec,nbeam)= ec_launchers_ids%beam(nbeam)% &
               launching_position%r(nt_ec)
          ZPIN_L(nt_ec,nbeam)= ec_launchers_ids%beam(nbeam)% &
               launching_position%z(nt_ec)
          PHIIN_L(nt_ec,nbeam)= ec_launchers_ids%beam(nbeam)% &
               launching_position%phi(nt_ec)
          ANGTIN_L(nt_ec,nbeam)= ec_launchers_ids%beam(nbeam)% &
               steering_angle_tor(nt_ec)
          ANGPIN_L(nt_ec,nbeam)= ec_launchers_ids%beam(nbeam)% &
               steering_angle_pol(nt_ec)
          OXFRAC_L(nt_ec,nbeam)= ec_launchers_ids%beam(nbeam)% &
               o_mode_fraction(nt_ec)
          PCURV1_L(nt_ec,nbeam)=ec_launchers_ids%beam(nbeam)% &
               phase%curvature(1,nt_ec)
          PCURV2_L(nt_ec,nbeam)=ec_launchers_ids%beam(nbeam)% &
               phase%curvature(2,nt_ec)
          PCURV2_L(nt_ec,nbeam)=ec_launchers_ids%beam(nbeam)% &
               phase%curvature(2,nt_ec)
          PANGLE_L(nt_ec,nbeam)=ec_launchers_ids%beam(nbeam)% &
               spot%angle(nt_ec)
          BSIZE1_L(nt_ec,nbeam)=ec_launchers_ids%beam(nbeam)% &
               spot%size(1,nt_ec)
          BSIZE2_L(nt_ec,nbeam)=ec_launchers_ids%beam(nbeam)% &
               spot%size(2,nt_ec)
          BANGLE_L(nt_ec,nbeam)=ec_launchers_ids%beam(nbeam)% &
               spot%angle(nt_ec)
       END DO
    END DO

    RETURN

  END SUBROUTINE get_ids_ec_launchers

  SUBROUTINE clear_ids_ec_launchers
    IMPLICIT NONE
    DEALLOCATE(BNAME_L,BIDNAME_L,RPINMIN_L,RPINMAX_L)
    DEALLOCATE(RFIN_L,RFIN_T,PECIN_L,PECIN_T,TIME_L)
    DEALLOCATE(RPIN_L,ZPIN_L,PHIIN_L,ANGTIN_L,ANGPIN_L,OXFRAC_L)
    DEALLOCATE(PCURV1_L,PCURV2_L,PANGLE_L,BSIZE1_L,BSIZE2_L,BANGLE_L)
    RETURN
  END SUBROUTINE clear_ids_ec_launchers
  
END MODULE ids_ec_launchers
