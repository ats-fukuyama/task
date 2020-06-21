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
       mdlobw, &              ! model id of output interval
       mdlobg                 ! model id of graphics
  REAL(rkind) &
       smax, &                ! maximum of orbit length
       dels, &                ! step size of orbit length
       eps_obt, &             ! convergence criterion of orbit solution
       del_obt, &             ! step size of iteration (initial condition)
       eps_nw                 ! convergence criterion of iteration (initial c.)

  REAL(rkind),DIMENSION(nobt_m):: &
       penergy_in, &          ! initial particle energy (mdlobi=0,1) [keV]
       pangle_in, &           ! initial sine of pitch angle (mdlobi=0,1) [mu/E]
       zeta_in, &             ! initial toroidal angle (mdlobi=0,1) [deg]
       pzeta_in, &            ! initial toroidal momentum (mdlobi=0) [P/E]
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
       zetab_ob, &            ! toroidal boozer angle zeta translated from obts
       thetab_ob, &           ! poloidal boozer angle theta
       pzeta_ob, &            ! toroidal momentum pzeta
       ptheta_ob, &           ! poloidal momentum ptheta
       psip_ob, &             ! poloidal magnetic flux
       rhopara_ob, &          ! parallel velocity devidede by cyclotron freq.
       babs_ob, &             ! absolute value of magnetic field
       phi_ob, &              ! electrostatic potential
       penergy_ob, &          ! particle energy
       pangle_ob, &           ! cos(pitch angle) (1:para,0:perp,-1:anti-para)
       rr_ob, &               ! major radius
       zz_ob, &               ! vertical positon
       zeta_ob, &             ! toroidal angle
       rs_ob, &               ! minor radius
       theta_ob               ! poloidal angle

  REAL(rkind):: &
       rr_axis,zz_axis        ! position of magnetic axis

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

    ALLOCATE(zetab_ob(0:nstp_max,nobt_max))
    ALLOCATE(thetab_ob(0:nstp_max,nobt_max))
    ALLOCATE(pzeta_ob(0:nstp_max,nobt_max))
    ALLOCATE(ptheta_ob(0:nstp_max,nobt_max))
    ALLOCATE(psip_ob(0:nstp_max,nobt_max))
    ALLOCATE(rhopara_ob(0:nstp_max,nobt_max))
    ALLOCATE(babs_ob(0:nstp_max,nobt_max))
    ALLOCATE(phi_ob(0:nstp_max,nobt_max))
    ALLOCATE(penergy_ob(0:nstp_max,nobt_max))
    ALLOCATE(pangle_ob(0:nstp_max,nobt_max))
    ALLOCATE(rr_ob(0:nstp_max,nobt_max))
    ALLOCATE(zz_ob(0:nstp_max,nobt_max))
    ALLOCATE(zeta_ob(0:nstp_max,nobt_max))
    ALLOCATE(rs_ob(0:nstp_max,nobt_max))
    ALLOCATE(theta_ob(0:nstp_max,nobt_max))
  END SUBROUTINE ob_allocate

  SUBROUTINE ob_deallocate
    IMPLICIT NONE

    DEALLOCATE(obt_in)
    DEALLOCATE(nstp_max_nobt)
    DEALLOCATE(obts)

    DEALLOCATE(zetab_ob)
    DEALLOCATE(thetab_ob)
    DEALLOCATE(pzeta_ob)
    DEALLOCATE(ptheta_ob)
    DEALLOCATE(psip_ob)
    DEALLOCATE(rhopara_ob)
    DEALLOCATE(babs_ob)
    DEALLOCATE(phi_ob)
    DEALLOCATE(penergy_ob)
    DEALLOCATE(pangle_ob)
    DEALLOCATE(rr_ob)
    DEALLOCATE(zz_ob)
    DEALLOCATE(zeta_ob)
    DEALLOCATE(rs_ob)
    DEALLOCATE(theta_ob)
  END SUBROUTINE ob_deallocate
END MODULE obcomm
