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

  INTEGER(ikind):: mdld        ! model type of diffusion coefficient
  INTEGER(ikind):: mdlprv      ! model type of Pereverzev method
  REAL(rkind)::    d0          ! lower diffusion coefficient
  REAL(rkind)::    d1          ! upper diffusion coefficient
  REAL(rkind)::    ltcr        ! critical scale length
  REAL(rkind)::    ph0         ! heating power density at r = 0
  REAL(rkind)::    phs         ! heating power density at r = a
  REAL(rkind)::    dprv1       ! enhanced diffusion coefficient
  REAL(rkind)::    dprv2       ! diffusion enhancement factor

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
       DFa,      &! diffusion coefficient [m^2/s]
       VCa,      &! convection velocity [m^2/s]
       EXa        ! heat exchange freuency [1/s]
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       PHa        ! heating power density [MW/m^3]

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
       RHV,XV_NEW,XV,XV_PREV
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       ELMTX,LHMTX
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       LIMTX,RPIMTX
  REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: &
       R1IMTX, R2IMTX, R3IMTX, RIMTX

! ----- diagnostics -----

  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       error_it,lt_save
  INTEGER(ikind):: &
       nitmax

! ----- save data parameters -----
  INTEGER(ikind):: ngt
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: gvt
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: gvts
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: gvrt
  REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: gvrts

! ----- species id -----
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
       ALLOCATE(DFa(neqmax,neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(VCa(neqmax,neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(EXa(neqmax,neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(PHa(neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(lt_save(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000

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
    IF(ALLOCATED(DFa)) DEALLOCATE(DFa)
    IF(ALLOCATED(VCa)) DEALLOCATE(VCa)
    IF(ALLOCATED(EXa)) DEALLOCATE(EXa)
    IF(ALLOCATED(PHa)) DEALLOCATE(PHa)
    IF(ALLOCATED(lt_save)) DEALLOCATE(lt_save)

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
    ALLOCATE(rpimtx(2,2,neqmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
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
    IF(ALLOCATED(rpimtx)) DEALLOCATE(rpimtx)
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
