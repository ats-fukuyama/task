MODULE trcomm

  USE plcomm
  IMPLICIT NONE

! ----- contral parameters -----

  INTEGER(ikind):: nrmax       ! number of radial step (except rg=0)
  INTEGER(ikind):: ntmax       ! number of time step
  REAL(rkind):: dt             ! size of time step
  REAL(rkind),DIMENSION(3,0:NSM):: &
       rg_fixed                ! minimum radius of fixed profile

  INTEGER(ikind):: nsamax      ! number of active particle species
  INTEGER(ikind),DIMENSION(nsm):: ns_nsa
                               ! table of NS for NSA
  INTEGER(ikind):: lmaxtr      ! maximum number of iterations
  REAL(rkind)::    epsltr      ! tolerance of iteration

  INTEGER(ikind):: mdltr_nc    ! model type of neoclassical transport coef.
  INTEGER(ikind):: mdltr_tb    ! model type of turbulent transport coef.
  INTEGER(ikind):: mdltr_prv   ! model type of Pereverzev method
  REAL(rkind)::    d0          ! lower diffusion coefficient
  REAL(rkind)::    d1          ! upper diffusion coefficient
  REAL(rkind)::    ltcr        ! critical scale length
  REAL(rkind)::    ph0         ! heating power density at r = 0
  REAL(rkind)::    phs         ! heating power density at r = a
  REAL(rkind)::    dprv1       ! enhanced diffusion coefficient
  REAL(rkind)::    dprv2       ! diffusion enhancement factor
  REAL(rkind)::    cdtrn       ! factor for particle diffusion
  REAL(rkind)::    cdtru       ! factor for toroidal flow viscosity
  REAL(rkind)::    cdtrt       ! factor for thermal diffusivity

  INTEGER(ikind):: ntstep      ! number of time step for status report
  INTEGER(ikind):: ngtmax      ! number of saved data
  INTEGER(ikind):: ngtstp      ! number of time step for data save

! ----- global variables -----

  REAL(rkind)::                t       ! time [s]
  REAL(rkind)::                dr      ! size of radial step
  INTEGER(ikind)::             neqmax  ! number of equations
  INTEGER(ikind)::             neqrmax ! number of active equations
  INTEGER(ikind)::             nvmax   ! number of variables
  INTEGER(ikind)::             nvrmax  ! number of active variables

! ----- plasma variables -----

  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       rg,      &! radial mesh position [m]
       qp,      &! safety factor
       qp_prev   ! previous,safety factor
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       rn,      &! particle density [10^{20] m^{-3}]
       ru,      &! toroidal velocity [m/s]
       rt,      &! temperature [keV]
       rn_prev, &! previous particle density [10^{20] m^{-3}]
       ru_prev, &! previous toroidal velocity [m/s]
       rt_prev   ! previous temperature [keV]
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       dtr,      &! diffusion coefficient [m^2/s]
       vtr,      &! convection velocity [m^2/s]
       ctr        ! exchange freuency [1/s]
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       str,      &! source density [MW/m^3]
       htr        ! current density [MA/m^2]
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       dtr_tb,   &! turbulent diffusion coefficient [m^2/s]
       vtr_tb,   &! turbulent convection velocity [m^2/s]
       dtr_nc,   &! neoclassical diffusion coefficient [m^2/s]
       vtr_nc,   &! neoclassical convection velocity [m^2/s]
       ctr_ex     ! exchange frequency [1/s]
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE::&
                  ! variables for Pereverzev method
       dtr_prv,  &! additional diffusion coefficient [m^2/s]
       vtr_prv,  &! additional convection velocity [m^2/s]
       add_prv    ! numerically additional term
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       eta_nc,   &! neoclassical resistivity [ohm m]
       jbs_nc,   &! bootstrap current by neoclassical effect [A/m^2]
       jcd_nb,   &! current driven by NBI [A/m^2]
       jcd_ec,   &! current driven by ECRF waves [A/m^2]
       jcd_lh,   &! current driven by LHRF waves [A/m^2]
       jcd_ic     ! current driven by ICRF waves [A/m^2]
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       lt_save    ! temperature scale length [m]

! ----- computation variables -----

  INTEGER(ikind),DIMENSION(:),ALLOCATABLE:: &
       nsa_neq,  &! particle species of equation neq
       nva_neq,  &! variable type of equation neq
       id_neq     ! boundary condition type of equation neq
  INTEGER(ikind),DIMENSION(:,:),ALLOCATABLE:: &
       id_neqnr   ! variable type of (neq,nr)

! ----- matrix equaton variables -----

  INTEGER(ikind),DIMENSION(:),ALLOCATABLE:: &
       neqr_neq,neq_neqr
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       rhv,xv_new,xv,xv_prev
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       elmtx,lhmtx
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       limtx, rjimtx, rsimtx
  REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: &
       r1imtx, r2imtx, r3imtx, rimtx

! ----- diagnostics -----

  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       error_it
  INTEGER(ikind):: &
       nitmax

! ----- save data parameters -----
  INTEGER(ikind):: ngt
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: gvt
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: gvts
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: gvrt
  REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: gvrts

! ----- species id -----
  INTEGER(ikind),DIMENSION(:),ALLOCATABLE:: idnsa
  CHARACTER(LEN=1),DIMENSION(:),ALLOCATABLE:: kidnsa

CONTAINS

! ----- allocation for radial profile -----

  SUBROUTINE tr_nr_allocate

    INTEGER(ikind),SAVE:: nrmax_save=0
    INTEGER(ikind),SAVE:: nsamax_save=0
    INTEGER(ikind)::      ierr

    IF(nrmax /= nrmax_save .OR. &
       nsamax /= nsamax_save ) THEN

       IF(nrmax_save /= 0) CALL tr_nr_deallocate

       ALLOCATE(rg(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(qp(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(qp_prev(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(rn(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(ru(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(rt(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(rn_prev(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(ru_prev(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(rt_prev(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(dtr(3*neqmax,3*neqmax,0:nrmax),STAT=ierr)
          IF(ierr /= 0) GOTO 9000
       ALLOCATE(vtr(3*neqmax,3*neqmax,0:nrmax),STAT=ierr) 
          IF(ierr /= 0) GOTO 9000
       ALLOCATE(ctr(3*neqmax,3*neqmax,0:nrmax),STAT=ierr)
          IF(ierr /= 0) GOTO 9000
       ALLOCATE(htr(3*neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(str(3*neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000

       ALLOCATE(dtr_tb(3*neqmax,3*neqmax,0:nrmax),STAT=ierr)
          IF(ierr /= 0) GOTO 9000
       ALLOCATE(vtr_tb(3*neqmax,3*neqmax,0:nrmax),STAT=ierr)
          IF(ierr /= 0) GOTO 9000
       ALLOCATE(dtr_nc(3*neqmax,3*neqmax,0:nrmax),STAT=ierr)
          IF(ierr /= 0) GOTO 9000
       ALLOCATE(vtr_nc(3*neqmax,3*neqmax,0:nrmax),STAT=ierr)
          IF(ierr /= 0) GOTO 9000
       ALLOCATE(ctr_ex(3*neqmax,3*neqmax,0:nrmax),STAT=ierr)
          IF(ierr /= 0) GOTO 9000
       ALLOCATE(eta_nc(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(jbs_nc(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(jcd_nb(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(jcd_ec(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(jcd_lh(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(jcd_ic(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000

       ALLOCATE(lt_save(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000

       ! for Pereverzev method
       ALLOCATE(dtr_prv(3*neqmax,0:nrmax),STAT=ierr)
          IF(ierr /= 0) GOTO 9000
       ALLOCATE(vtr_prv(3*neqmax,0:nrmax),STAT=ierr)
          IF(ierr /= 0) GOTO 9000
       ALLOCATE(add_prv(3*neqmax,0:nrmax),STAT=ierr)
          IF(ierr /= 0) GOTO 9000

       nsamax_save=nsamax
       nrmax_save=nrmax
    END IF

    RETURN
9000 WRITE(6,*) 'XX tr_nr_allocate: allocation error: ierr=',ierr
    STOP
  END SUBROUTINE tr_nr_allocate

  SUBROUTINE tr_nr_deallocate

    IF(ALLOCATED(rg)) DEALLOCATE(rg)
    IF(ALLOCATED(qp)) DEALLOCATE(qp)
    IF(ALLOCATED(qp_prev)) DEALLOCATE(qp_prev)
    IF(ALLOCATED(rn)) DEALLOCATE(rn)
    IF(ALLOCATED(ru)) DEALLOCATE(ru)
    IF(ALLOCATED(rt)) DEALLOCATE(rt)
    IF(ALLOCATED(rn_prev)) DEALLOCATE(rn_prev)
    IF(ALLOCATED(ru_prev)) DEALLOCATE(ru_prev)
    IF(ALLOCATED(rt_prev)) DEALLOCATE(rt_prev)
    IF(ALLOCATED(dtr)) DEALLOCATE(dtr)
    IF(ALLOCATED(vtr)) DEALLOCATE(vtr)
    IF(ALLOCATED(ctr)) DEALLOCATE(ctr)
    IF(ALLOCATED(htr)) DEALLOCATE(htr)
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(dtr_nc)) DEALLOCATE(dtr_nc)
    IF(ALLOCATED(dtr_tb)) DEALLOCATE(dtr_tb)
    IF(ALLOCATED(vtr_nc)) DEALLOCATE(vtr_nc)
    IF(ALLOCATED(vtr_tb)) DEALLOCATE(vtr_tb)
    IF(ALLOCATED(ctr_ex)) DEALLOCATE(ctr_ex)
    IF(ALLOCATED(eta_nc)) DEALLOCATE(eta_nc)
    IF(ALLOCATED(jbs_nc)) DEALLOCATE(jbs_nc)
    IF(ALLOCATED(jcd_nb)) DEALLOCATE(jcd_nb)
    IF(ALLOCATED(jcd_ec)) DEALLOCATE(jcd_ec)
    IF(ALLOCATED(jcd_lh)) DEALLOCATE(jcd_lh)
    IF(ALLOCATED(jcd_ic)) DEALLOCATE(jcd_ic)
    IF(ALLOCATED(lt_save)) DEALLOCATE(lt_save)

    ! for Pereverzev method
    IF(ALLOCATED(dtr_prv)) DEALLOCATE(dtr_prv)
    IF(ALLOCATED(vtr_prv)) DEALLOCATE(vtr_prv)
    IF(ALLOCATED(add_prv)) DEALLOCATE(add_prv)

    RETURN
  END SUBROUTINE tr_nr_deallocate

! ----- allocation for coefficient calculation -----

  SUBROUTINE tr_neq_allocate

    INTEGER(ikind),SAVE:: neqmax_save=0
    INTEGER(ikind),SAVE:: nrmax_save=0
    INTEGER(ikind),SAVE:: nvmax_save=0
    INTEGER(ikind)::      ierr

    IF(nrmax == nrmax_save .AND. &
       neqmax == neqmax_save .AND. &
       nvmax == nvmax_save) RETURN

    IF(neqmax_save /= 0) CALL tr_neq_deallocate

    ALLOCATE(nsa_neq(neqmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(nva_neq(neqmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(id_neq(neqmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(id_neqnr(neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000

    ALLOCATE(xv(nvmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(xv_new(nvmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(xv_prev(nvmax),STAT=ierr); IF(ierr /= 0) GOTO 9000

    ALLOCATE(neqr_neq(neqmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(limtx(2,2,neqmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(rjimtx(2,2,neqmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(rsimtx(2,2,neqmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(rimtx(2,2,neqmax,neqmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(r1imtx(2,2,neqmax,neqmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(r2imtx(2,2,neqmax,neqmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(r3imtx(2,2,neqmax,neqmax),STAT=ierr); IF(ierr /= 0) GOTO 9000

    nrmax_save=nrmax
    neqmax_save=neqmax
    nvmax_save=nvmax

    RETURN
9000 WRITE(6,*) 'XX tr_neq_allocate: allocation error: ierr=',ierr
    STOP
  END SUBROUTINE tr_neq_allocate

  SUBROUTINE tr_neq_deallocate

    IF(ALLOCATED(nsa_neq))  DEALLOCATE(nsa_neq)
    IF(ALLOCATED(nva_neq))  DEALLOCATE(nva_neq)
    IF(ALLOCATED(id_neq))   DEALLOCATE(id_neq)
    IF(ALLOCATED(id_neqnr)) DEALLOCATE(id_neqnr)

    IF(ALLOCATED(xv)) DEALLOCATE(xv)
    IF(ALLOCATED(xv_new)) DEALLOCATE(xv_new)
    IF(ALLOCATED(xv_prev)) DEALLOCATE(xv_prev)

    IF(ALLOCATED(neqr_neq)) DEALLOCATE(neqr_neq)
    IF(ALLOCATED(limtx)) DEALLOCATE(limtx)
    IF(ALLOCATED(rjimtx)) DEALLOCATE(rjimtx)
    IF(ALLOCATED(rsimtx)) DEALLOCATE(rsimtx)
    IF(ALLOCATED(rimtx)) DEALLOCATE(rimtx)
    IF(ALLOCATED(r1imtx)) DEALLOCATE(r1imtx)
    IF(ALLOCATED(r2imtx)) DEALLOCATE(r2imtx)
    IF(ALLOCATED(r3imtx)) DEALLOCATE(r3imtx)

    RETURN
  END SUBROUTINE tr_neq_deallocate

! ----- allocation for solver -----

  SUBROUTINE tr_neqr_allocate

    INTEGER(ikind),SAVE:: neqrmax_save=0
    INTEGER(ikind),SAVE:: nvrmax_save=0
    INTEGER(ikind)::      ierr

    IF(neqrmax == neqrmax_save .AND. &
       nvrmax == nvrmax_save) RETURN

    IF(neqrmax_save /= 0) CALL tr_neqr_deallocate

    ALLOCATE(neq_neqr(neqrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(rhv(nvrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(elmtx(2*neqrmax,2*neqrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
    ALLOCATE(lhmtx(4*neqrmax-1,nvrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000

    neqrmax_save=neqrmax
    nvrmax_save=nvrmax

    RETURN
9000 WRITE(6,*) 'XX tr_neqr_allocate: allocation error: ierr=',ierr
    STOP
  END SUBROUTINE tr_neqr_allocate

  SUBROUTINE tr_neqr_deallocate

    IF(ALLOCATED(neq_neqr)) DEALLOCATE(neq_neqr)
    IF(ALLOCATED(rhv)) DEALLOCATE(rhv)
    IF(ALLOCATED(elmtx)) DEALLOCATE(elmtx)
    IF(ALLOCATED(lhmtx)) DEALLOCATE(lhmtx)

    RETURN
  END SUBROUTINE tr_neqr_deallocate

! ----- allocation for diagnostics -----

  SUBROUTINE tr_nit_allocate

    INTEGER(ikind),SAVE:: lmaxtr_save=0
    INTEGER(ikind)::      ierr

    IF(lmaxtr == lmaxtr_save) RETURN

    IF(lmaxtr_save /= 0) CALL tr_nit_deallocate

    ALLOCATE(error_it(lmaxtr),STAT=ierr); IF(ierr /= 0) GOTO 9000

    lmaxtr_save=lmaxtr

    RETURN
9000 WRITE(6,*) 'XX tr_nit_allocate: allocation error: ierr=',ierr
    STOP
  END SUBROUTINE tr_nit_allocate

  SUBROUTINE tr_nit_deallocate

    IF(ALLOCATED(error_it)) DEALLOCATE(error_it)

    RETURN
  END SUBROUTINE tr_nit_deallocate

! ----- allocation for data save -----

  SUBROUTINE tr_ngt_allocate

    INTEGER(ikind),SAVE:: ngtmax_save=0
    INTEGER(ikind),SAVE:: nsamax_save=0
    INTEGER(ikind),SAVE:: nrmax_save=0
    INTEGER(ikind)::      ierr

    IF(ngtmax /= ngtmax_save .OR. &
       nsamax /= nsamax_save .OR. &
       nrmax  /= nrmax_save ) THEN

       IF(ngtmax_save /= 0) CALL tr_ngt_deallocate

       ALLOCATE(gvt(0:ngtmax,0:2),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(gvts(0:ngtmax,nsamax,3),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(gvrt(0:nrmax,0:ngtmax,1),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(gvrts(0:nrmax,0:ngtmax,nsamax,3),STAT=ierr)
            IF(ierr /= 0) GOTO 9000

       ngtmax_save=ngtmax
       nsamax_save=nsamax
       nrmax_save=nrmax
    END IF

    RETURN
9000 WRITE(6,*) 'XX tr_ngt_allocate: allocation error: ierr=',ierr
    STOP
  END SUBROUTINE tr_ngt_allocate

  SUBROUTINE tr_ngt_deallocate

    IF(ALLOCATED(gvt)) DEALLOCATE(gvt)
    IF(ALLOCATED(gvts)) DEALLOCATE(gvts)
    IF(ALLOCATED(gvrt)) DEALLOCATE(gvrt)
    IF(ALLOCATED(gvrts)) DEALLOCATE(gvrts)

    RETURN
  END SUBROUTINE tr_ngt_deallocate

! ----- allocation for species id -----

  SUBROUTINE tr_nsa_allocate

    INTEGER(ikind),SAVE:: nsamax_save=0
    INTEGER(ikind)::      ierr

    IF(nsamax /= nsamax_save) THEN

       IF(nsamax_save /= 0) CALL tr_nsa_deallocate

       ALLOCATE(idnsa(nsamax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(kidnsa(nsamax),STAT=ierr); IF(ierr /= 0) GOTO 9000

       nsamax_save=nsamax
    END IF

    RETURN
9000 WRITE(6,*) 'XX tr_nsa_allocate: allocation error: ierr=',ierr
    STOP
  END SUBROUTINE tr_nsa_allocate

  SUBROUTINE tr_nsa_deallocate

    IF(ALLOCATED(kidnsa)) DEALLOCATE(kidnsa)

    RETURN
  END SUBROUTINE tr_nsa_deallocate

END MODULE trcomm
