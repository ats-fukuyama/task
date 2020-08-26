! obcomm.f90

MODULE obcomm_parm
  USE commpi
  USE plcomm_parm
  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER:: &
       nobt_m=100             ! maximum number of orbits

! --- input parameters ---

  INTEGER:: &
       nobt_max, &            ! number of orbits
       nstp_max, &            ! maximum number of orbit steps
       ns_ob, &               ! id of particle species
       lmax_nw                ! maximum number of iteration (initial condition)
  INTEGER:: &
       mdlobp, &              ! model id of equation of motion
       mdlobi, &              ! model id of input scheme of initial parameters
       mdlobq, &              ! model id of ODE solver
       mdlobt, &              ! model id of time normalization
       mdlobc, &              ! model id of cycle (0: not cycle, 1: one cycle)
       mdlobw, &              ! model id of output interval
       mdlobg                 ! model id of graphics
  REAL(rkind) &
       tmax, &                ! maximum of orbit following time in omega_bounce
       delt, &                ! time step size in onega_bounce
       eps_obt, &             ! convergence criterion of orbit solution
       del_obt, &             ! step size of iteration (initial condition)
       eps_nw                 ! convergence criterion of iteration (initial c.)

  REAL(rkind),DIMENSION(nobt_m):: &
       penergy_in, &          ! initial particle energy (mdlobi=0,1) [keV]
       pcangle_in, &          ! initial cosine of pitch angle (mdlobi=0,1)
       zeta_in, &             ! initial toroidal angle (mdlobi=0,1) [deg]
       psipn_in, &            ! initial normalized flux (mdlobi=0)
       theta_in, &            ! initial poloidal angle (mdlobi=0) [deg]
       rr_in, &               ! initial major radius (mdlobi=1) [m]
       zz_in                  ! initial vertical position (mdlobi=1) [m]

  INTEGER:: nrmax_ob, &       ! number of radial mesh
            nthmax_ob, &      ! number of pooidal mesh
            nsumax_ob         ! number plasma surface points (same as wall)

CONTAINS

  SUBROUTINE open_obcomm_parm
    RETURN
  END SUBROUTINE open_obcomm_parm

END MODULE obcomm_parm

MODULE obcomm
  USE obcomm_parm
  USE plcomm
  USE commpi
  IMPLICIT NONE

  INTEGER,PARAMETER:: neq_max=4

  INTEGER,DIMENSION(:),ALLOCATABLE:: &
       nstp_max_nobt           ! number of steps for nobt
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       obt_in                 ! initial condition of orbit equation
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       obts                   ! solution of orbit equation
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       time_ob, &             ! time
       zetab_ob, &            ! toroidal boozer angle zeta translated from obts
       thetab_ob, &           ! poloidal boozer angle theta
       psip_ob, &             ! poloidal magnetic flux
       rhopara_ob, &          ! parallel velocity devidede by cyclotron freq.
       pzeta_ob, &            ! toroidal momentum pzeta
       ptheta_ob, &           ! poloidal momentum ptheta
       babs_ob, &             ! absolute value of magnetic field
       phi_ob, &              ! electrostatic potential
       vpara_ob, &            ! parallel velocity
       vperp_ob, &            ! perpendicular velocity
       psit_ob, &             ! toroidal magnetic flux
       zeta_ob, &             ! toroidal angle
       rr_ob, &               ! major radius
       zz_ob, &               ! vertical positon
       rs_ob, &               ! minor radius
       theta_ob               ! poloidal angle

CONTAINS
       
  SUBROUTINE ob_allocate
    IMPLICIT NONE
    INTEGER,SAVE:: init=0
    INTEGER,SAVE:: nobt_max_save=0
    INTEGER,SAVE:: nstp_max_save=0

    IF(init.EQ.0) THEN
       init=1
    ELSE
       IF((nobt_max.EQ.nobt_max_save).AND. &
          (nstp_max.EQ.nstp_max_save)) RETURN
       CALL ob_deallocate
    END IF

    ALLOCATE(obt_in(neq_max,nobt_max))
    ALLOCATE(nstp_max_nobt(nobt_max))
    ALLOCATE(obts(0:neq_max,0:nstp_max,nobt_max))

    ALLOCATE(time_ob(0:nstp_max,nobt_max))
    ALLOCATE(zetab_ob(0:nstp_max,nobt_max))
    ALLOCATE(thetab_ob(0:nstp_max,nobt_max))
    ALLOCATE(psip_ob(0:nstp_max,nobt_max))
    ALLOCATE(rhopara_ob(0:nstp_max,nobt_max))
    ALLOCATE(pzeta_ob(0:nstp_max,nobt_max))
    ALLOCATE(ptheta_ob(0:nstp_max,nobt_max))
    ALLOCATE(babs_ob(0:nstp_max,nobt_max))
    ALLOCATE(phi_ob(0:nstp_max,nobt_max))
    ALLOCATE(vpara_ob(0:nstp_max,nobt_max))
    ALLOCATE(vperp_ob(0:nstp_max,nobt_max))
    ALLOCATE(psit_ob(0:nstp_max,nobt_max))
    ALLOCATE(zeta_ob(0:nstp_max,nobt_max))
    ALLOCATE(rr_ob(0:nstp_max,nobt_max))
    ALLOCATE(zz_ob(0:nstp_max,nobt_max))
    ALLOCATE(rs_ob(0:nstp_max,nobt_max))
    ALLOCATE(theta_ob(0:nstp_max,nobt_max))
  END SUBROUTINE ob_allocate

  SUBROUTINE ob_deallocate
    IMPLICIT NONE

    DEALLOCATE(obt_in)
    DEALLOCATE(nstp_max_nobt)
    DEALLOCATE(obts)

    DEALLOCATE(time_ob)
    DEALLOCATE(zetab_ob)
    DEALLOCATE(thetab_ob)
    DEALLOCATE(psip_ob)
    DEALLOCATE(rhopara_ob)
    DEALLOCATE(pzeta_ob)
    DEALLOCATE(ptheta_ob)
    DEALLOCATE(babs_ob)
    DEALLOCATE(phi_ob)
    DEALLOCATE(vpara_ob)
    DEALLOCATE(vperp_ob)
    DEALLOCATE(psit_ob)
    DEALLOCATE(zeta_ob)
    DEALLOCATE(rr_ob)
    DEALLOCATE(zz_ob)
    DEALLOCATE(rs_ob)
    DEALLOCATE(theta_ob)
  END SUBROUTINE ob_deallocate
END MODULE obcomm
