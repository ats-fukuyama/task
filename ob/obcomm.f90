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
       nobtmax, &             ! number of orbits
       nstpmax, &             ! maximum number of orbits
       ns_ob, &               ! id of particle species
       lmaxnw                 ! maximum number of iteration (initial condition)
  INTEGER:: &
       mdlobi, &              ! model id of initial input parameters
       mdlobg                 ! model id of graphics
  REAL(rkind) &
       smax, &                ! maximum of orbit length
       dels, &                ! step size of orbit length
       epsobt, &              ! convergence criterion of orbit solution
       delobt, &              ! step size of iteration (initial condition)
       epsnw                  ! convergence criterion of iteration (initial c.)

  REAL(rkind),DIMENSION(nobt_m):: &
       zetab_in, &            ! initial toroidal boozer angle
       thetab_in, &           ! initial poloidal boozer angle
       pzeta_in, &            ! initial toroidal momentum
       ptheta_in, &           ! initial poloidal momentum
       psip_in, &             ! initial poloidal magnetic flux
       rhopara_in, &          ! initial parallel velocity divided by 
       penerg_in, &           ! initial particle energy (mdlobi=1)
       pangle_in, &           ! initial sine of pitch angle (mdlobj=1)
       rr_in, &               ! initial major radius (mdlobj=1)
       zz_in, &               ! initial vertical position (mdlobj=1)
       zeta_in                ! initial toroidal angle (mdlobj=1)

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

  INTEGER,PARAMETER:: neq=6

  INTEGER,DIMENSION(:),ALLOCATABLE:: &
       nstpmax_nobt           ! number of steps for nobt
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
       penerg_ob, &           ! particle energy
       pangle_ob, &           ! pitch angle
       rr_ob, &               ! major radius
       zz_ob, &               ! vertical positon
       zeta_ob, &             ! toroidal angle
       rs_ob, &               ! minor radius
       theta_ob               ! poloidal angle

CONTAINS
       
  SUBROUTINE ob_allocate
    IMPLICIT NONE
    INTEGER,SAVE:: init=0
    INTEGER,SAVE:: nobtmax_save=0
    INTEGER,SAVE:: nstpmax_save=0

    IF(init.EQ.0) THEN
       init=1
    ELSE
       IF((nobtmax.EQ.nobtmax_save).AND. &
          (nstpmax.EQ.nstpmax_save)) RETURN
       CALL ob_deallocate
    END IF

    ALLOCATE(obt_in(neq,nobtmax))
    ALLOCATE(nstpmax_nobt(nobtmax))
    ALLOCATE(obts(0:neq,0:nstpmax,nobtmax))

    ALLOCATE(zetab_ob(0:nstpmax,nobtmax))
    ALLOCATE(thetab_ob(0:nstpmax,nobtmax))
    ALLOCATE(pzeta_ob(0:nstpmax,nobtmax))
    ALLOCATE(ptheta_ob(0:nstpmax,nobtmax))
    ALLOCATE(psip_ob(0:nstpmax,nobtmax))
    ALLOCATE(rhopara_ob(0:nstpmax,nobtmax))
    ALLOCATE(babs_ob(0:nstpmax,nobtmax))
    ALLOCATE(phi_ob(0:nstpmax,nobtmax))
    ALLOCATE(penerg_ob(0:nstpmax,nobtmax))
    ALLOCATE(pangle_ob(0:nstpmax,nobtmax))
    ALLOCATE(rr_ob(0:nstpmax,nobtmax))
    ALLOCATE(zz_ob(0:nstpmax,nobtmax))
    ALLOCATE(zeta_ob(0:nstpmax,nobtmax))
    ALLOCATE(rs_ob(0:nstpmax,nobtmax))
    ALLOCATE(theta_ob(0:nstpmax,nobtmax))
  END SUBROUTINE ob_allocate

  SUBROUTINE ob_deallocate
    IMPLICIT NONE

    DEALLOCATE(obt_in)
    DEALLOCATE(nstpmax_nobt)
    DEALLOCATE(obts)

    DEALLOCATE(zetab_ob)
    DEALLOCATE(thetab_ob)
    DEALLOCATE(pzeta_ob)
    DEALLOCATE(ptheta_ob)
    DEALLOCATE(psip_ob)
    DEALLOCATE(rhopara_ob)
    DEALLOCATE(babs_ob)
    DEALLOCATE(phi_ob)
    DEALLOCATE(penerg_ob)
    DEALLOCATE(pangle_ob)
    DEALLOCATE(rr_ob)
    DEALLOCATE(zz_ob)
    DEALLOCATE(zeta_ob)
    DEALLOCATE(rs_ob)
    DEALLOCATE(theta_ob)
  END SUBROUTINE ob_deallocate
END MODULE obcomm
