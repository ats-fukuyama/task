MODULE trcomm

  USE plcomm
  IMPLICIT NONE

!  --- module bpsd_constants ---
!  real(rkind),parameter :: pi   = 3.14159265358979323846_dp
!  real(rkind),parameter :: aee  = 1.602176487E-19_dp ! elementary charge
!  real(rkind),parameter :: ame  = 9.10938215E-31_dp  ! electron mass
!  real(rkind),parameter :: amp  = 1.672621637E-27_dp ! proton mass [kg]
!  real(rkind),parameter :: vc   = 2.99792458E8_dp    ! speed of light [m/s]
!  real(rkind),parameter :: rmu0 = 4.E-7_dp*PI        ! permeability
!  real(rkind),parameter :: eps0 = ONE/(VC*VC*RMU0)   ! permittivity

  REAL(rkind),PARAMETER :: rkev = aee*1.d3 ! the factor ([keV] -> [J])

! ----- contral parameters -----

  INTEGER(ikind):: nrmax       ! number of radial step (except rg=0)
  INTEGER(ikind):: ntmax       ! number of time step
  REAL(rkind)::    dt          ! size of time step
  REAL(rkind),DIMENSION(3,0:NSM):: &
                   rg_fixed    ! minimum radius of fixed profile

  INTEGER(ikind):: nsamax      ! number of active particle species
  INTEGER(ikind),DIMENSION(nsm)::  &
                   ns_nsa      ! table of NS for NSA
  INTEGER(ikind):: lmaxtr      ! maximum number of iterations
  REAL(rkind)::    epsltr      ! tolerance of iteration

  INTEGER(ikind):: mdltr_nc    ! model type of neoclassical transport coef.
  INTEGER(ikind):: mdltr_tb    ! model type of turbulent transport coef.
  INTEGER(ikind):: mdltr_prv   ! model type of Pereverzev method
  REAL(rkind)::    dtr0        ! lower diffusion coefficient in simple model
  REAL(rkind)::    dtr1        ! upper diffusion coefficient in simple model
  REAL(rkind)::    ltcr        ! critical scale length in simple model
  REAL(rkind)::    ph0         ! heating power density at r = 0
  REAL(rkind)::    phs         ! heating power density at r = a
  REAL(rkind)::    dprv1       ! enhanced diffusion coefficient in Prv method
  REAL(rkind)::    dprv2       ! diffusion enhancement factor in Prv method
  REAL(rkind)::    cdtrn       ! factor for particle diffusion
  REAL(rkind)::    cdtru       ! factor for toroidal flow viscosity
  REAL(rkind)::    cdtrt       ! factor for thermal diffusivity

  INTEGER(ikind):: ntstep      ! number of time step for status report
  INTEGER(ikind):: ngtmax      ! number of saved data
  INTEGER(ikind):: ngtstp      ! number of time step for data save

! ----- global variables -----

  REAL(rkind)::        t       ! time [s]
  REAL(rkind)::        dr      ! size of radial step
  INTEGER(ikind)::     neqmax  ! number of equations
  INTEGER(ikind)::     neqrmax ! number of active equations
  INTEGER(ikind)::     nvmax   ! number of variables
  INTEGER(ikind)::     nvrmax  ! number of active variables

! ----- plasma variables -----

  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       rg,      &! radial mesh position [m]
       rm,      &! radial half mesh position [m]
       qp,      &! safety factor
       qp_prev   ! previous safety factor
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
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
                  ! variables for Pereverzev method
       dtr_prv,  &! additional diffusion coefficient [m^2/s]
       vtr_prv    ! additional convection velocity [m^2/s]
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       eta_nc,   &! neoclassical resistivity [ohm m]
       jbs_nc,   &! bootstrap current by neoclassical effect [A/m^2]
       jcd_nb,   &! current driven by NBI [A/m^2]
       jcd_ec,   &! current driven by ECRF waves [A/m^2]
       jcd_lh,   &! current driven by LHRF waves [A/m^2]
       jcd_ic     ! current driven by ICRF waves [A/m^2]
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       lt_save    ! temperature scale length [m]

! ----- profile variables -----

  REAL(rkind),DIMENSION(:),ALLOCATABLE::&
       vtor,     &! toroidal rotation velocity [m/s]
       vpol       ! poloidal rotation velocity [m/s]

! --- equilibrium interface variables -----
! ------ interface variables fot bpsd_equ1D
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       psitrho,  &! Toroidal magnetic flux [Wb] ~ pi*r^2*B
       psiprho,  &! Poloidal magnetic flux [Wb] ~ pi*R*r*Bp
       ppprho,   &! Plasma pressure [Pa]
       piqrho,   &! Inverse of safety factor, iota
       ttrho,    &! R*B
       pirho      ! Toroidal current [A] ~ 2*pi*r*Bp/mu_0

! ------ interface variables for bpsd_metric1D
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       pvolrho,  &! Plasma volume
       psurrho,  &! Plasma surface
       dvrho,    &! dV/drho
       rdpvrho,  &! dpsi/dV
       arrho,    &! <1/R^2> ; aver2i
       abb2rho,  &! <B^2>
       aib2rho,  &! <1/B^2>
       abvrho,   &! <|grad V|^2/R^2>
       ar1rho,   &! <|grad rho|>
       ar2rho,   &! <|grad rho|^2>
       abrho,    &! <|grad rho|^2/R^2> ; avegrr2
       rmjrho,   &! local R
       rmnrho,   &! local r
       rkprho,   &! local kappa
       abb1rho,  &! <B>
       arhbrho,  &! not used
       epsrho,   &! r/R
!
       rmnrhom,  &! local r (half-grid)
       rkprhom    ! local kappa (half-grid)

! ----- normalized variables -----
  REAL(rkind) :: &
       rhoa       ! normalized minor radius
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       dpdrho,   &! dpsi/drho
       rhog,     &! normalized minor radius mesh position
       rhom,     &! normalized minor radius half-mesh position
       rjcb       ! 1/rho : rho ~ kappa * r : effective minor radius ?

! ----- unclassified -----
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       bp,       &! poloidal magnetic field [T]
       er         ! radial electric field [V/m]

! ----- switch variables -----
  INTEGER(ikind) :: &
       mdluf,    &! model type of UFILE
       mdler,    &! model type of radial electric field
       mdleqn     !
         
!       modelg,   &! control plasma geometry model
!       nteqit     ! step interval of EQ calculation

! ----- daignostic variables -----
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       nrd1,     &! a diagnostic array for radial grid
       nrd2,     &! a diagnostic array for radial grid
       nrd3,     &! a diagnostic array for radial grid
       nrd4       ! a diagnostic array for radial grid
  
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
  REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: gvrts,gparts

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

       IF(nrmax_save /= 0 ) CALL tr_nr_deallocate

       ALLOCATE(rg(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(rm(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
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

       ! profile variables
       ALLOCATE(vtor(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
       ALLOCATE(vpol(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000

       ! for Pereverzev method
       ALLOCATE(dtr_prv(3*neqmax,0:nrmax),STAT=ierr)
          IF(ierr /= 0) GOTO 9000
       ALLOCATE(vtr_prv(3*neqmax,0:nrmax),STAT=ierr)
          IF(ierr /= 0) GOTO 9000

      ! geometric factors
      ! +-- interface variables for bpsd_equ1D
      ALLOCATE(psitrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(psiprho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(ppprho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(piqrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(ttrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(pirho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000

      ! +-- interface variables for bpsd_metric1D
      ALLOCATE(pvolrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(psurrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(dvrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000 
      ALLOCATE(rdpvrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(arrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000 
      ALLOCATE(abb2rho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(aib2rho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(abvrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000 
      ALLOCATE(ar1rho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(ar2rho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000 
      ALLOCATE(abrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000 
      ALLOCATE(rmjrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000 
      ALLOCATE(rmnrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000 
      ALLOCATE(rkprho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000 
      ALLOCATE(abb1rho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000 
      ALLOCATE(arhbrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000 
      ALLOCATE(epsrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000 

      ALLOCATE(rmnrhom(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(rkprhom(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000  

      !    normalized variables
      ALLOCATE(dpdrho(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(rhog(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000   
      ALLOCATE(rhom(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000   
      ALLOCATE(rjcb(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000   

      ! unclassified
      ALLOCATE(bp(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(er(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000

      ! for diagnostic
      ALLOCATE(nrd1(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(nrd2(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(nrd3(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      ALLOCATE(nrd4(0:nrmax),STAT=ierr); IF(ierr /= 0) GOTO 9000
      nrd1 = 0.d0
      nrd2 = 0.d0
      nrd3 = 0.d0
      nrd4 = 0.d0
     
      nrmax_save  = nrmax
      nsamax_save = nsamax
    END IF

    RETURN
9000 WRITE(6,*) 'XX tr_nr_allocate: allocation error: ierr=',ierr
    STOP
  END SUBROUTINE tr_nr_allocate

  SUBROUTINE tr_nr_deallocate

    IF(ALLOCATED(rg)) DEALLOCATE(rg)
    IF(ALLOCATED(rm)) DEALLOCATE(rm)
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

    ! profile variables
    IF(ALLOCATED(vtor)) DEALLOCATE(vtor)
    IF(ALLOCATED(vpol)) DEALLOCATE(vpol)

    ! for Pereverzev method
    IF(ALLOCATED(dtr_prv)) DEALLOCATE(dtr_prv)
    IF(ALLOCATED(vtr_prv)) DEALLOCATE(vtr_prv)

    ! interface variables for bpsd_equ1D
    IF(ALLOCATED(psitrho)) DEALLOCATE(psitrho)   
    IF(ALLOCATED(psiprho)) DEALLOCATE(psiprho)   
    IF(ALLOCATED(ppprho)) DEALLOCATE(ppprho)   
    IF(ALLOCATED(piqrho)) DEALLOCATE(piqrho)   
    IF(ALLOCATED(ttrho)) DEALLOCATE(ttrho)   
    IF(ALLOCATED(pirho)) DEALLOCATE(pirho)   

    ! interface variables for bpsd_metric1D
    IF(ALLOCATED(pvolrho)) DEALLOCATE(pvolrho)
    IF(ALLOCATED(psurrho)) DEALLOCATE(psurrho)
    IF(ALLOCATED(dvrho)) DEALLOCATE(dvrho)
    IF(ALLOCATED(rdpvrho)) DEALLOCATE(rdpvrho)
    IF(ALLOCATED(arrho)) DEALLOCATE(arrho)
    IF(ALLOCATED(abb2rho)) DEALLOCATE(abb2rho)
    IF(ALLOCATED(aib2rho)) DEALLOCATE(aib2rho)
    IF(ALLOCATED(abvrho)) DEALLOCATE(abvrho)
    IF(ALLOCATED(ar1rho)) DEALLOCATE(ar1rho)
    IF(ALLOCATED(ar2rho)) DEALLOCATE(ar2rho)
    IF(ALLOCATED(abrho)) DEALLOCATE(abrho)
    IF(ALLOCATED(rmjrho)) DEALLOCATE(rmjrho)
    IF(ALLOCATED(rmnrho)) DEALLOCATE(rmnrho)
    IF(ALLOCATED(rkprho)) DEALLOCATE(rkprho)
    IF(ALLOCATED(abb1rho)) DEALLOCATE(abb1rho)
    IF(ALLOCATED(arhbrho)) DEALLOCATE(arhbrho)
    IF(ALLOCATED(epsrho)) DEALLOCATE(epsrho)

    IF(ALLOCATED(rmnrhom)) DEALLOCATE(rmnrhom)
    IF(ALLOCATED(rkprhom)) DEALLOCATE(rkprhom)

    ! normalized variables
    IF(ALLOCATED(dpdrho)) DEALLOCATE(dpdrho)
    IF(ALLOCATED(rhog)) DEALLOCATE(rhog)
    IF(ALLOCATED(rhom)) DEALLOCATE(rhom)
    IF(ALLOCATED(rjcb)) DEALLOCATE(rjcb)

    ! unclassified
    IF(ALLOCATED(bp)) DEALLOCATE(bp)
    IF(ALLOCATED(er)) DEALLOCATE(er)

    ! for diagnostic
    IF(ALLOCATED(nrd1)) DEALLOCATE(nrd1)
    IF(ALLOCATED(nrd2)) DEALLOCATE(nrd2)
    IF(ALLOCATED(nrd3)) DEALLOCATE(nrd3)
    IF(ALLOCATED(nrd4)) DEALLOCATE(nrd4)
       
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

    nrmax_save  = nrmax
    neqmax_save = neqmax
    nvmax_save  = nvmax

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

    neqrmax_save = neqrmax
    nvrmax_save  = nvrmax

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

! ----- allocation for diagnostics of iteration -----

  SUBROUTINE tr_nit_allocate

    INTEGER(ikind),SAVE:: lmaxtr_save=0
    INTEGER(ikind)::      ierr

    IF(lmaxtr == lmaxtr_save) RETURN

    IF(lmaxtr_save /= 0) CALL tr_nit_deallocate

    ALLOCATE(error_it(lmaxtr),STAT=ierr); IF(ierr /= 0) GOTO 9000

    lmaxtr_save = lmaxtr

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

       ! for Pereverzev method
       ALLOCATE(gparts(0:nrmax,0:ngtmax,nsamax,3),STAT=ierr)
            IF(ierr /= 0) GOTO 9000

       nrmax_save  = nrmax
       ngtmax_save = ngtmax
       nsamax_save = nsamax
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
    IF(ALLOCATED(gparts)) DEALLOCATE(gparts)

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

       nsamax_save = nsamax
    END IF

    RETURN
9000 WRITE(6,*) 'XX tr_nsa_allocate: allocation error: ierr=',ierr
    STOP
  END SUBROUTINE tr_nsa_allocate

  SUBROUTINE tr_nsa_deallocate

    IF(ALLOCATED(idnsa)) DEALLOCATE(idnsa)
    IF(ALLOCATED(kidnsa)) DEALLOCATE(kidnsa)

    RETURN
  END SUBROUTINE tr_nsa_deallocate

END MODULE trcomm
