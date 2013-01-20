MODULE trcomm

  USE plcomm
  IMPLICIT NONE

!!$                   *** module bpsd_constants ***
!!$
!!$  real(rkind),parameter :: pi   = 3.14159265358979323846_dp
!!$  real(rkind),parameter :: aee  = 1.602176487E-19_dp ! elementary charge
!!$  real(rkind),parameter :: ame  = 9.10938215E-31_dp  ! electron mass
!!$  real(rkind),parameter :: amp  = 1.672621637E-27_dp ! proton mass [kg]
!!$  real(rkind),parameter :: vc   = 2.99792458E8_dp    ! speed of light [m/s]
!!$  real(rkind),parameter :: rmu0 = 4.E-7_dp*PI        ! permeability
!!$  real(rkind),parameter :: eps0 = ONE/(VC*VC*RMU0)   ! permittivity

!!$               *** Variables defined in plcomm.f90 ***
!!$
!!$  INTEGER,PARAMETER:: NSM=100 ! Maximum number of particle species
!!$
!!$  INTEGER:: NSMAX,MODELG,MODELN,MODELQ,IDEBUG,MODEFR,MODEFW
!!$
!!$  REAL(rkind):: RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ
!!$  REAL(rkind):: PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2
!!$  REAL(rkind):: RHOMIN,QMIN,RHOEDG,RHOITB,RHOGMN,RHOGMX
!!$
!!$  REAL(rkind),DIMENSION(NSM)::                       &
!!$       &      PA,PZ,PZ0,PN,PNS,PTPR,PTPP,PTS,PU,PUS, &
!!$       &      PNITB,PTITB,PUITB
!!$
!!$  CHARACTER(len=80):: KNAMEQ,KNAMWR,KNAMFP,KNAMWM,KNAMPF,KNAMFO,KNAMTR
!!$

  REAL(rkind),PARAMETER :: rkev = aee*1.d3 ! the factor ([keV] -> [J])

  ! *** the maximum size of UFILE ***
  INTEGER(ikind),PARAMETER :: ntum = 1001
  INTEGER(ikind),PARAMETER :: nrum = 201
  INTEGER(ikind),PARAMETER :: nsum = 10   ! number of species: NMn + electron

! ----- contral parameters -----
  INTEGER(ikind):: nrmax       ! number of radial step (except rg=0)
  INTEGER(ikind):: ntmax       ! number of time step
  REAL(rkind)::    dt          ! size of time step
  REAL(rkind),DIMENSION(3,0:NSM):: &
                   rg_fixed    ! minimum radius of fixed profile

  INTEGER(ikind):: nsamax      ! number of active particle species
  INTEGER(ikind):: nsafmax     ! number of active (fast ion) particle species
  INTEGER(ikind):: nsabmax     ! number of active bulk particle species
                               !  (including electron)
  INTEGER(ikind):: nsanmax     ! number of active (neutral) particle species
  INTEGER(ikind),DIMENSION(nsm)::  &
                   idnsa,     &! species id (-1:e, 0:n, 1:bulk, 2:fast)
                   ns_nsa,    &! conversion table of NS for NSA
                   nsab_nsa,  &!
                   nsaf_nsa,  &!
                   nsan_nsa,  &!
                   nsab_nsaf   ! conversion table of [fast ion > bulk ion]
  REAL(rkind)::    pa_mion     ! atomic number of main hydrogenic ion
  REAL(rkind)::    pz_mion     ! charge number of main hydrogenic ion
  REAL(rkind)::    pa_mimp     ! atomic number of main impurity ion
  REAL(rkind)::    pz_mimp     ! charge number of main impurity ion
  INTEGER(ikind):: lmaxtr      ! maximum number of iterations
  REAL(rkind)::    epsltr      ! tolerance of iteration

  INTEGER(ikind):: mdltr_nc    ! model type of neoclassical transport coef.
  INTEGER(ikind):: mdltr_tb    ! model type of turbulent transport coef.
  INTEGER(ikind):: mdltr_prv   ! model type of Pereverzev method
  REAL(rkind)::    dtr0        ! lower diffusion coefficient in simple model
  REAL(rkind)::    dtr1        ! upper diffusion coefficient in simple model
  REAL(rkind)::    ltcr        ! critical scale length in simple model
  REAL(rkind)::    ph0         ! heating power density 
                               !  at r = 0 in simple model [MW/m^3]
  REAL(rkind)::    phs         ! heating power density 
                               !  at r = a in simple model [MW/m^3]
  REAL(rkind)::    dprv1       ! enhanced diffusion coefficient in Prv method
  REAL(rkind)::    dprv2       ! diffusion enhancement factor in Prv method
  REAL(rkind)::    rhog_prv    ! enhanced diffusion region in Prv method
  REAL(rkind)::    cdtrn       ! factor for particle diffusion
  REAL(rkind)::    cdtru       ! factor for toroidal flow viscosity
  REAL(rkind)::    cdtrt       ! factor for thermal diffusivity

  INTEGER(ikind):: ntstep      ! number of time step for status report
  INTEGER(ikind):: ngtmax      ! number of saved data
  INTEGER(ikind):: ngtstp      ! number of time step for data save

  INTEGER(ikind):: mdlnb,mdlec,mdllh,mdlic
  REAL(rkind)::    pnbtot     ! total power of NBI heating [MW]
  REAL(rkind)::    pnbrw      ! radial width of deposition [m]
  REAL(rkind)::    pnbcd      ! current drive factor
  REAL(rkind)::    pnbr0      ! radial position of deposition
  REAL(rkind)::    pnbeng     ! [keV]
  REAL(rkind)::    pectot,pictot,plhtot ! total power of RF heating [MW]
  REAL(rkind)::    pecrw, picrw, plhrw  ! radial width of deposition [m]
  REAL(rkind)::    peccd, piccd, plhcd  ! current drice factor
  REAL(rkind)::    pecr0, picr0, plhr0  ! radial position of deposition [m]
  REAL(rkind)::    pectoe,pictoe,plhtoe ! power partition to electron
  REAL(rkind)::    pecnpr,picnpr,plhnpr ! parallel refractive index


! ----- switch variables -----
  INTEGER(ikind) :: &
       mdleqb,   &! EQs. selection parameter for B_theta(d psi/d rho)
       mdleqn,   &! EQs. selection parameter for density
       mdlequ,   &! EQs. selection parameter for rotation
       mdleqt,   &! EQs. selection parameter for heat
       mdler,    &! model type of radial electric field
       nteqit     ! step interval of EQ calculation

  INTEGER(ikind) :: &
       mdlglb,   &!
       mdlgmt,   &!
       mdlsrc,   &!
       mdlijq


! ----- global variables -----
  REAL(rkind)::      t       ! time [s]
  REAL(rkind)::      t_prev  ! time at the previous time step
  REAL(rkind)::      dr      ! size of radial step
  INTEGER(ikind)::   neqmax  ! number of equations
  INTEGER(ikind)::   neqrmax ! number of active equations
  INTEGER(ikind)::   nvmax   ! number of variables
  INTEGER(ikind)::   nvrmax  ! number of active variables

  REAL(rkind)::      wp_t    ! total stored energy
  REAL(rkind)::      wp_th   ! total stored energy of thermal particles
  REAL(rkind)::      wp_inc  ! incremental stored energy of thermal particles
  REAL(rkind)::      wpu_inc ! incremental stored energy of thermal particles
  REAL(rkind)::      rw      ! the deviation of Wp_inc
  REAL(rkind)::      taue1   ! energy confinment time (steady state)
  REAL(rkind)::      taue2   ! energy confinment time (transient)
  REAL(rkind)::      taue3   ! energy confinment time (thermal, transient)
  REAL(rkind)::      taue89  ! ITER89-P L-mode scaling
  REAL(rkind)::      taue98  ! IPB98(y,2) H-mode scaling with ELMs
  REAL(rkind)::      h89     ! confinement enhancement factor for ITER89-P
  REAL(rkind)::      h98y2   ! confinement enhancement factor for IPB98(y,2)
  REAL(rkind)::      betan   ! normalized toroidal beta

  REAL(rkind)::      pin_t   ! [MW/m^3]
  REAL(rkind)::      poh_t   !
  REAL(rkind)::      pnb_t   !
  REAL(rkind)::      pec_t   !
  REAL(rkind)::      pibw_t  !
  REAL(rkind)::      pic_t   !
  REAL(rkind)::      plh_t   !
  REAL(rkind)::      pnf_t   !
  REAL(rkind)::      prl_t   !
  REAL(rkind)::      pout_t  !

  REAL(rkind),DIMENSION(:),ALLOCATABLE:: & ! 0:nrmax
       beta,     &! toroidal beta
       beta_va,  &! volume-averaged toroidal beta
       betap,    &! poloidal beta
       betap_va, &! volume-averaged poloidal beta
       betaq      ! toroidal beta for reaction rate

  REAL(rkind),DIMENSION(:),ALLOCATABLE:: & ! 1:nsamax
       rns_va,   &! volume-averaged density
       rts_va,   &! volume-averaged temperature
       ws_t,     &! stored energy of each species
       stdrt,    &! the deviations of each temperature profiles
       offrt      ! the offsets of each temperature profiles

! ----- plasma variables -----
  REAL(rkind) :: phia   ! total toroidal flux enclosed by the plasma [Wb]
  REAL(rkind) :: rips   ! toroidal current at the beginning [MA]
  REAL(rkind) :: ripe   ! toroidal current at the end       [MA]
  REAL(rkind) :: vloop  ! loop valtage

  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       rg,      &! radial mesh position [m]
       rm,      &! radial half mesh position [m]
       qp,      &! safety factor
       qp_prev, &! previous safety factor
       dpdrho,   &! d psi/d rho
       dpdrho_prev! d psi/d rho
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       rn,      &! particle density [10^{20] m^{-3}]
       ru,      &! toroidal velocity [m/s]
       rt,      &! particle temperature [keV]
       rn_prev, &! previous particle density [10^{20] m^{-3}]
       ru_prev, &! previous toroidal velocity [m/s]
       rt_prev   ! previous temperature [keV]
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       dtr,      &! diffusion coefficient [m^2/s]
       vtr,      &! convection velocity [m^2/s]
       ctr        ! charge exchange freuency [1/s]
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       str,      &! source density [keV/(10^{20}m^3 s)]
       htr        ! current density [A/m^2]
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       ! ***_tb, ***_nc are variables except for magnetic equation.
       dtr_tb,   &! turbulent diffusion coefficient [m^2/s]
       vtr_tb,   &! turbulent convection velocity [m/s]
       dtr_nc,   &! neoclassical diffusion coefficient [m^2/s]
       vtr_nc,   &! neoclassical convection velocity [m/s]
       ctr_ex     ! charge (energy) exchange frequency [1/s]
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       fluxtb,   &! fluxes by turbulent transport [***/(m^2 s)]
       fluxnc,   &! fluxes by neoclassical transport [***/(m^2 s)]
       grdpf,    &! gradients of the variables to be solved
       dtr_nl     ! assumed non-linear turbulent transport coefficient [m^2/s]
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
                  ! variables for Pereverzev method
       dtr_prv,  &! additional diffusion coefficient [m^2/s]
       vtr_prv    ! additional convection velocity [m/s]
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       htr_simple !  simple model of external driven current density [A/m^2]
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       eta   ,   &! pararell resistivity [ohm m]
       jtot  ,   &! (parallel) total current density [A/m^2]
       joh   ,   &! (parallel) ohmic current density [A/m^2]
       jtor  ,   &! toroidal current density [A/m^2]
       eta_nc,   &! neoclassical resistivity [ohm m]
       etam_nc,  &! neoclassical resistivity (half grid) [ohm m]
       jbs_nc,   &! bootstrap current density by neoclassical effect [A/m^2]
       jex_nc,   &! externally driven current density (NCLASS) [A/m^2]
       jcd_nb,   &! current density driven by NBI [A/m^2]
       jcd_ec,   &! current density driven by ECRF waves [A/m^2]
       jcd_lh,   &! current density driven by LHRF waves [A/m^2]
       jcd_ic     ! current density driven by ICRF waves [A/m^2]
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       str_simple ! simple model of source density [MW/m^3]
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       ptot,    & ! total heating power density [W/m^3]
       poh,     & ! ohmic heating power density [W/m^3]
       pnb,     & ! NBI heating power density [W/m^3]
       pec,     & ! EC heating power density [W/m^3]
       pibw,    & ! IBW heating power density [W/m^3]
       pic,     & ! IC heating power density [W/m^3]
       plh,     & ! LH heating power density [W/m^3]
       pnf,     & ! heating power density due to fusion alpha [W/m^3]
       prl,     & !  [W/m^3]
       pwl,     & !  [W/m^3]
!
       snb,     & ! [10^{20] m^{-3} s^{-1}]
       spl,     & ! [10^{20] m^{-3} s^{-1}]
       swl        ! [10^{20] m^{-3} s^{-1}]

! ----- profile variables -----
  REAL(rkind) :: profj1, profj2 ! current density profile factors

  REAL(rkind),DIMENSION(:),ALLOCATABLE::&
       vtor,     &! toroidal rotation velocity [m/s]
       vpol,     &! poloidal rotation velocity [m/s]
       vpar,     &! parallel roration velocity [m/s]
       vprp,     &! perpendicular rotation velocity [m/s]
       wrot       ! toroidal angular speed [rad/s]

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
       rmjrho,   &! local R [m]
       rmnrho,   &! local r [m]
       rkprho,   &! local kappa
       abb1rho,  &! <B>
       epsrho     ! r/R

! ----- normalized variables -----
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       rhog,     &! normalized minor radius mesh position
       rhom,     &! normalized minor radius half-mesh position
       rjcb       ! 1/rho : rho ~ kappa * r : effective minor radius ?

! ----- derivatives of the quantities -----
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE ::&
       rp,       &!the pressure of each species (nT) [Pa]
       rp_d       !the deriv. of pressure of each species (dnT/dr) 
       
  REAL(rkind),DIMENSION(:),ALLOCATABLE ::&
       rp_tot,   &! the total pressure
       rp_totd,  &! the deriv. of total pressure
       rp_add,   &! the additional pressure
       rp_beam,  &! the beam pressure
!
       rt_e,     &! the electron temperature
       rt_ed,    &! the deriv. of electron temperature
       rt_ecl,   &! the scale length of electron temperature 
       rt_i,     &! the effective hydrogenic ion temperature
       rt_id,    &! the deriv. of effective hydrogenic ion temperature 
       rt_icl,   &! the scale length of hydrogenic ion temperature 
!
       rn_e,     &! the electron density
       rn_ed,    &! the deriv. of electron density 
       rn_ecl,   &! the scale length of electron density 
       rn_i,     &! the sum of hydrogenic ion density
       rn_id,    &! the deriv. of hydrogenic ion density 
       rn_icl,   &! the scale length of hydrogenic ion density 
       qp_d,     &! the deriv. of safety factor
!
       ai_ave     ! mean atomic mass of thermal ions [AMU]

  REAL(rkind),DIMENSION(:),ALLOCATABLE ::&
       mshear,   &! magnetic shear            r/q * (dq/dr)
       mshear_cl,&! magnetic shear length  R*q**2/(r*dq/dr)
       mcurv,    &! magnetic curvature
       vexbp,    &! ExBp velocity [m/s]
       dvexbpdr, &! ExBp velocity gradient [1/s]
       wexbp,    &! ExBp shearing rate [rad/s]
!       v_alf,    &! Alfven wave velocity
       v_se,     &! speed of sound for electron
       alpha      ! MHD alpha

  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       z_eff       ! Z_eff: effective charge

! ----- unclassified -----
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       bp,       &! poloidal magnetic field [T]
       er,       &! radial electric field [V/m]
       ezoh       !

! ----- daignostic variables for debug -----
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
       elmtx,lhmtx,elmtx_ofd
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: &
       limtx, rjimtx, rsimtx
  REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: &
       r1imtx, r2imtx, r3imtx, rimtx, r1imtx_ofd

! ----- diagnostics -----
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       error_it
  INTEGER(ikind):: &
       nitmax

! ----- save data parameters -----
  INTEGER(ikind):: ngt
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: gvt,gvtu
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: gvts
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: gvrt
  REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: gvrts

  INTEGER(ikind),DIMENSION(5):: &
       unitid     ! unit identifier for input and output
                  ! 1: csv output of radial profiles
                  ! 2: csv output of time evolutions
                  ! 3: UFILE output of radial profiles
                  ! 4: UFILE output of time evolutions
                  ! 5: csv input
                  ! 6: UFILE output

  INTEGER(ikind) :: &
       mdlwrt,   &! switch for writing down output
       nwrstp     ! the interval number of time step for writing down output

  CHARACTER(LEN=30) :: &
       kwpnam,   &! the filename for writing down output of radial profiles
       kwtnam     ! the filename for writing down output of time evolutions

  INTEGER(ikind) :: &  
       mdlxp,    &! Select UFILE or MDSplus
       mdluf,    &! Model type of UFILE
       mdlugt,   &! Select the way to set the time of snap shot for graphic
       mdlni,    &! Select how to determine main ion density, impurity density
                  !  or effective charge number
                  ! 1 : complete n_i and n_imp  from Zeff, n_e (and n_bulk)
                  ! 2 : complete n_imp and Zeff from n_e, n_i (and n_bulk)
                  ! 3 : complete n_i and Zeff   from n_e, n_imp (and n_bulk)
                  ! 8 : complete n_i and Zeff using the assumption n_e = n_i (+ n_fast)
                  ! 9 : complete n_i and Zeff using the assumption n_e = n_i
       ufid_bin   ! Parameter which determines how to handle exp. files.
                  ! 0 : Binary files are loaded if available, or ASCII files
                  !      are loaded and aftermath binary files are created.
                  ! 1 : Only binary files are loaded.
                  ! 2 : Only ASCII files are loaded and binary files are NOT
                  !      created.

  ! ----- UFILE control -----
  CHARACTER(80) :: kuf_dir  ! UFILE database directory path
  CHARACTER(80) :: kuf_dev  ! device name
  CHARACTER(80) :: kuf_dcg  ! discharge number
  CHARACTER(80) :: kdirx    ! directory containig a set of UFILE
  REAL(rkind)   :: time_slc ! time slicing point for steady state simulation
  REAL(rkind)   :: time_snap! time slicing point for graphic

  REAL(rkind)   :: uf_tinit !
  REAL(rkind)   :: uf_tdura !

  ! ----- Stored variables for UFILE -----
  INTEGER(ikind) :: ndmax  ! number of 0D experimental (UFILE) data
  INTEGER(ikind) :: ntxmax ! number of experimental time data
  INTEGER(ikind) :: ntlmax ! number of time step of experimental data
                           ! * Step size 'dt' is that of simulation.
  REAL(rkind)    :: tlmax  ! end of time of experimental data

  INTEGER(ikind),DIMENSION(1:nsum) :: nsa_nsu ! conversion table [nsu  -> nsa]
  INTEGER(ikind),DIMENSION(1:nsum) :: nsa_nsfu! conversion table [nsuf -> nsa]

  REAL(rkind),DIMENSION(1:ntum) :: & 
       ! global variables
       tmu,   &! time point data vector [s]
       rru,   &! plasma major radius [m]
       rau,   &! plasma minor radius
       phiau, &! total toroidal flux enclosed by the plasma [Wb]
       bbu,   &! vaccume toroidal field at geometric axis [T]
       rkapu, &! plasma elongation
       rdltu, &! mean triangularity of the plasma boundary
       ripu,  &! plasma current [A]
       wthu,  &! thermal plasma energy content [J]
       wtotu, &! total plasma energy content [J]
       zeffu, &! line averaged effective charge
!
       ! source
       pnbu,  &! total injected NBI power (minus shine through) [W]
       pecu,  &! coupled ECH power [W]
       pibwu, &! coupled IBW power [W]
       picu,  &! coupled ICRH power [W]
       plhu,  &! coupled LH power [W]
       pohmu, &! total ohmic power [W]
       pradu   ! total radiated power [W]

  REAL(rkind),DIMENSION(1:ntum,1:nrum) :: & ! profile variables
       qpu,    &! safety factor
       bpu,    &! surface averaged poloidal magnetic field [T]
       zeffru, &! plasma effective charge radial profile
       wrotu,  &! toroidal angular speed [rad/s]
!
       jtotu,  &! total current density [A/m^2]
       jnbu,   &! NBI driven current profile [A/m^2]
       jecu,   &! ECH driven current profile [A/m^2]
       jicu,   &! ICRH driven current profile [A/m^2]
       jlhu,   &! LH driven current profile [A/m^2]
       jbsu     ! bootstrap current profile [A/m^2]
!
  REAL(rkind),DIMENSION(2,1:ntum,1:nrum) :: & ! electron and ion
       qnbu,   &! power deposition profile on thermal particles by NBI [W/m^3]
       qecu,   &! power deposition profile on thermal particles by ECH [W/m^3]
       qibwu,  &! power deposition profile on thermal particles by IBW [W/m^3]
       qicu,   &! power deposition profile on thermal particles by ICRH[W/m^3]
       qlhu,   &! power deposition profile on thermal particles by LH  [W/m^3]
       qwallu, &! heat loss due to ionisation of wall neutrals [W/m^3]
       qfusu,  &! heating density due to DT fusion [W/m^2]
!
       snbu     ! source of thermal particles from NBI [/(m^3 s)]
!               ! in the case of ions, it's due to thermalization and 
!               !  charge exchange process [/(m^3 s)]

  REAL(rkind),DIMENSION(1:ntum,1:nrum) :: & ! profile variables
       qohmu,  &! ohmic power density [W/m^3]
       qradu,  &! total radiated power density [W/m^3]
!
       swallu, &! main thermal ion particle source from 
!               !  due to ionisation recycling wall neutral [/(m^3 s)]
!
       pvolu,  &! volume enclosed by the magnetic surface [m^3]
       psuru,  &! surface area of the magnetic surface [m^2]
       rmjrhou,&! geometrical major radius [m]
       rmnrhou,&! geometrical minor radius [m]
       ar1rhou,&! metric quantity: <|nabla rho|>
       ar2rhou,&! metric quantity: <|nabla rho|^2>
       rkprhou,&! average elongation of the magnetic surface
       dvrhou, &! d V /d rho
       arrhou, &! <1/R^2>
       abrhou, &! <|grad rho|^2/R^2>
       ttrhou   ! R*B

  REAL(rkind),DIMENSION(1:nsum,1:ntum,1:nrum) :: &
       rnu,    &! particle density [10^{20} m^{-3}]
       rnfu,   &! fast particle density [10^{20} m^{-3}]
       rtu,    &! particle temperature [keV]
       rpu      ! particle pressure [Pa]
  REAL(rkind),DIMENSION(1:nsum) :: &
       pau,    &! atomic number
       pzu,    &! charge number
       pafu,   &! atomic number of fast ions (nsu=1 is dummy)
       pzfu     ! charge number of fast ions (nsu=1 is dummy)
  
CONTAINS

  SUBROUTINE tr_nr_allocate
  ! ----- allocation for radial profile -----

    INTEGER(ikind),SAVE:: nrmax_save=0
    INTEGER(ikind),SAVE:: nsamax_save=0
    INTEGER(ikind)::      ierr

    IF(nrmax /= nrmax_save .OR. &
       nsamax /= nsamax_save ) THEN

       IF(nrmax_save /= 0 ) CALL tr_nr_deallocate

       nr_allocation : DO
          ALLOCATE(rg(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(rm(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(qp(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(qp_prev(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(dpdrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT 
          ALLOCATE(dpdrho_prev(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT 
          ALLOCATE(rn(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(ru(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(rt(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(rn_prev(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(ru_prev(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(rt_prev(nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(dtr(neqmax,neqmax,0:nrmax),STAT=ierr)
            IF(ierr /= 0) EXIT
          ALLOCATE(vtr(neqmax,neqmax,0:nrmax),STAT=ierr) 
            IF(ierr /= 0) EXIT
          ALLOCATE(ctr(neqmax,neqmax,0:nrmax),STAT=ierr)
            IF(ierr /= 0) EXIT
          ALLOCATE(str(neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(htr(neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
            
          ALLOCATE(dtr_tb(neqmax,neqmax,0:nrmax),STAT=ierr)
            IF(ierr /= 0) EXIT
          ALLOCATE(vtr_tb(neqmax,neqmax,0:nrmax),STAT=ierr)
            IF(ierr /= 0) EXIT
          ALLOCATE(dtr_nc(neqmax,neqmax,0:nrmax),STAT=ierr)
            IF(ierr /= 0) EXIT
          ALLOCATE(vtr_nc(neqmax,neqmax,0:nrmax),STAT=ierr)
            IF(ierr /= 0) EXIT
          ALLOCATE(ctr_ex(neqmax,neqmax,0:nrmax),STAT=ierr)
            IF(ierr /= 0) EXIT

          ALLOCATE(jtot(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(joh(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(jtor(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(eta_nc(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(etam_nc(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(jbs_nc(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(jex_nc(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(jcd_nb(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(jcd_ec(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(jcd_lh(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(jcd_ic(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(eta(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(htr_simple(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(str_simple(neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(ptot(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(poh(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(pnb(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(pec(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(pibw(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(pic(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(plh(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(prl(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(pwl(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(pnf(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          
          ALLOCATE(snb(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(spl(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(swl(1:nsamax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          
          ! profile variables
          ALLOCATE(vtor(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vpol(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vpar(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vprp(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(wrot(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          
          ! for Pereverzev method
          ALLOCATE(dtr_prv(neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(vtr_prv(neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(fluxtb(neqmax,0:nrmax),STAT=ierr);  IF(ierr /= 0) EXIT
          ALLOCATE(fluxnc(neqmax,0:nrmax),STAT=ierr);  IF(ierr /= 0) EXIT
          ALLOCATE(grdpf(neqmax,0:nrmax),STAT=ierr);   IF(ierr /= 0) EXIT
          ALLOCATE(dtr_nl(neqmax,0:nrmax),STAT=ierr);  IF(ierr /= 0) EXIT

          ! geometric factors
          ! +-- interface variables for bpsd_equ1D
          ALLOCATE(psitrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(psiprho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(ppprho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(piqrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(ttrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(pirho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          
          ! +-- interface variables for bpsd_metric1D
          ALLOCATE(pvolrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(psurrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(dvrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT 
          ALLOCATE(rdpvrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(arrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT 
          ALLOCATE(abb2rho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(aib2rho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(abvrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT 
          ALLOCATE(ar1rho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(ar2rho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT 
          ALLOCATE(abrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT 
          ALLOCATE(rmjrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT 
          ALLOCATE(rmnrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT 
          ALLOCATE(rkprho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT 
          ALLOCATE(abb1rho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT 
          ALLOCATE(epsrho(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT 
                    
          !    normalized variables
          ALLOCATE(rhog(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT   
          ALLOCATE(rhom(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT   
          ALLOCATE(rjcb(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT   
          
          ! unclassified
          ALLOCATE(bp(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(er(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(ezoh(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          
          ! global variables
          ALLOCATE(beta(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(beta_va(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(betap(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(betap_va(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(betaq(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          
          ! for diagnostic
          ALLOCATE(nrd1(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(nrd2(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(nrd3(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(nrd4(0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT

          ALLOCATE(rp(nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(rp_d(nsamax,0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          
          ALLOCATE(rp_tot(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT   
          ALLOCATE(rp_totd(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT  
          ALLOCATE(rp_add(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT   
          ALLOCATE(rp_beam(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT  
          
          ALLOCATE(rt_e(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(rt_ed(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(rt_ecl(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT   
          ALLOCATE(rt_i(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(rt_id(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(rt_icl(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT   
          
          ALLOCATE(rn_e(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(rn_ed(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(rn_ecl(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(rn_i(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT     
          ALLOCATE(rn_id(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(rn_icl(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT        
          ALLOCATE(qp_d(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT      

          ALLOCATE(ai_ave(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT      
          
          ALLOCATE(mshear(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT   
          ALLOCATE(mcurv(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(vexbp(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(dvexbpdr(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT 
          ALLOCATE(wexbp(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT       
!          ALLOCATE(v_alf(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT    
          ALLOCATE(v_se(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT
          ALLOCATE(alpha(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT

          ALLOCATE(z_eff(0:nrmax),STAT=ierr); IF(ierr /=0 ) EXIT

          nrmax_save  = nrmax
          nsamax_save = nsamax
          RETURN
       END DO nr_allocation
       
       WRITE(6,*) 'XX tr_nr_allocate: allocation error: ierr=',ierr
       STOP
    END IF

    RETURN
  END SUBROUTINE tr_nr_allocate

  SUBROUTINE tr_nr_deallocate

    IF(ALLOCATED(rg)) DEALLOCATE(rg)
    IF(ALLOCATED(rm)) DEALLOCATE(rm)
    IF(ALLOCATED(qp)) DEALLOCATE(qp)
    IF(ALLOCATED(qp_prev)) DEALLOCATE(qp_prev)
    IF(ALLOCATED(dpdrho)) DEALLOCATE(dpdrho)
    IF(ALLOCATED(dpdrho_prev)) DEALLOCATE(dpdrho_prev)
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
    IF(ALLOCATED(eta)) DEALLOCATE(eta)
    IF(ALLOCATED(jtot)) DEALLOCATE(jtot)
    IF(ALLOCATED(joh)) DEALLOCATE(joh)
    IF(ALLOCATED(jtor)) DEALLOCATE(jtor)
    IF(ALLOCATED(eta_nc)) DEALLOCATE(eta_nc)
    IF(ALLOCATED(etam_nc)) DEALLOCATE(etam_nc)
    IF(ALLOCATED(jbs_nc)) DEALLOCATE(jbs_nc)
    IF(ALLOCATED(jex_nc)) DEALLOCATE(jex_nc)
    IF(ALLOCATED(jcd_nb)) DEALLOCATE(jcd_nb)
    IF(ALLOCATED(jcd_ec)) DEALLOCATE(jcd_ec)
    IF(ALLOCATED(jcd_lh)) DEALLOCATE(jcd_lh)
    IF(ALLOCATED(jcd_ic)) DEALLOCATE(jcd_ic)

    IF(ALLOCATED(htr_simple)) DEALLOCATE(htr_simple)
    IF(ALLOCATED(str_simple)) DEALLOCATE(str_simple)
    IF(ALLOCATED(ptot)) DEALLOCATE(ptot)
    IF(ALLOCATED(poh)) DEALLOCATE(poh)
    IF(ALLOCATED(pnb)) DEALLOCATE(pnb)
    IF(ALLOCATED(pec)) DEALLOCATE(pec)
    IF(ALLOCATED(pibw)) DEALLOCATE(pibw)
    IF(ALLOCATED(pic)) DEALLOCATE(pic)
    IF(ALLOCATED(plh)) DEALLOCATE(plh)
    IF(ALLOCATED(prl)) DEALLOCATE(prl)
    IF(ALLOCATED(pwl)) DEALLOCATE(pwl)
    IF(ALLOCATED(pnf)) DEALLOCATE(pnf)

    IF(ALLOCATED(snb)) DEALLOCATE(snb)
    IF(ALLOCATED(spl)) DEALLOCATE(spl)
    IF(ALLOCATED(swl)) DEALLOCATE(swl)

    ! profile variables
    IF(ALLOCATED(vtor)) DEALLOCATE(vtor)
    IF(ALLOCATED(vpol)) DEALLOCATE(vpol)
    IF(ALLOCATED(vpar)) DEALLOCATE(vpar)
    IF(ALLOCATED(vprp)) DEALLOCATE(vprp)
    IF(ALLOCATED(wrot)) DEALLOCATE(wrot)

    ! for Pereverzev method
    IF(ALLOCATED(dtr_prv)) DEALLOCATE(dtr_prv)
    IF(ALLOCATED(vtr_prv)) DEALLOCATE(vtr_prv)
    IF(ALLOCATED(fluxtb)) DEALLOCATE(fluxtb)
    IF(ALLOCATED(fluxnc)) DEALLOCATE(fluxnc)
    IF(ALLOCATED(grdpf)) DEALLOCATE(grdpf)
    IF(ALLOCATED(dtr_nl)) DEALLOCATE(dtr_nl)

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
    IF(ALLOCATED(epsrho)) DEALLOCATE(epsrho)

    ! normalized variables
    IF(ALLOCATED(dpdrho)) DEALLOCATE(dpdrho)
    IF(ALLOCATED(rhog)) DEALLOCATE(rhog)
    IF(ALLOCATED(rhom)) DEALLOCATE(rhom)
    IF(ALLOCATED(rjcb)) DEALLOCATE(rjcb)

    ! unclassified
    IF(ALLOCATED(bp)) DEALLOCATE(bp)
    IF(ALLOCATED(er)) DEALLOCATE(er)
    IF(ALLOCATED(ezoh)) DEALLOCATE(ezoh)

    IF(ALLOCATED(beta)) DEALLOCATE(beta)
    IF(ALLOCATED(beta_va)) DEALLOCATE(beta_va)
    IF(ALLOCATED(betap)) DEALLOCATE(betap)
    IF(ALLOCATED(betap_va)) DEALLOCATE(betap_va)
    IF(ALLOCATED(betaq)) DEALLOCATE(betaq)

    ! for diagnostic
    IF(ALLOCATED(nrd1)) DEALLOCATE(nrd1)
    IF(ALLOCATED(nrd2)) DEALLOCATE(nrd2)
    IF(ALLOCATED(nrd3)) DEALLOCATE(nrd3)
    IF(ALLOCATED(nrd4)) DEALLOCATE(nrd4)

! ***

    IF(ALLOCATED(rp)) DEALLOCATE(rp)
    IF(ALLOCATED(rp_d)) DEALLOCATE(rp_d)
!    
    IF(ALLOCATED(rp_tot)) DEALLOCATE(rp_tot)
    IF(ALLOCATED(rp_totd)) DEALLOCATE(rp_totd)
    IF(ALLOCATED(rp_add)) DEALLOCATE(rp_add)
    IF(ALLOCATED(rp_beam)) DEALLOCATE(rp_beam)
!
    IF(ALLOCATED(rt_e)) DEALLOCATE(rt_e)
    IF(ALLOCATED(rt_ed)) DEALLOCATE(rt_ed)
    IF(ALLOCATED(rt_ecl)) DEALLOCATE(rt_ecl)
    IF(ALLOCATED(rt_i)) DEALLOCATE(rt_i)
    IF(ALLOCATED(rt_id)) DEALLOCATE(rt_id)
    IF(ALLOCATED(rt_icl)) DEALLOCATE(rt_icl)
!
    IF(ALLOCATED(rn_e)) DEALLOCATE(rn_e)
    IF(ALLOCATED(rn_ed)) DEALLOCATE(rn_ed)
    IF(ALLOCATED(rn_ecl)) DEALLOCATE(rn_ecl)
    IF(ALLOCATED(rn_i)) DEALLOCATE(rn_i)
    IF(ALLOCATED(rn_id)) DEALLOCATE(rn_id)
    IF(ALLOCATED(rn_icl)) DEALLOCATE(rn_icl)
    IF(ALLOCATED(qp_d)) DEALLOCATE(qp_d)
!
    IF(ALLOCATED(ai_ave)) DEALLOCATE(ai_ave)
!       
    IF(ALLOCATED(mshear)) DEALLOCATE(mshear)
    IF(ALLOCATED(mcurv)) DEALLOCATE(mcurv)
    IF(ALLOCATED(vexbp)) DEALLOCATE(vexbp)
    IF(ALLOCATED(dvexbpdr)) DEALLOCATE(dvexbpdr)
    IF(ALLOCATED(wexbp)) DEALLOCATE(wexbp)
!    IF(ALLOCATED(v_alf)) DEALLOCATE(v_alf)
    IF(ALLOCATED(v_se)) DEALLOCATE(v_se)
    IF(ALLOCATED(alpha)) DEALLOCATE(alpha)

    IF(ALLOCATED(z_eff)) DEALLOCATE(z_eff)

    RETURN
  END SUBROUTINE tr_nr_deallocate


  SUBROUTINE tr_neq_allocate
  ! ----- allocation for coefficient calculation -----

    INTEGER(ikind),SAVE:: neqmax_save=0
    INTEGER(ikind),SAVE:: nrmax_save=0
    INTEGER(ikind),SAVE:: nvmax_save=0
    INTEGER(ikind)::      ierr

    IF(nrmax  /= nrmax_save .OR. &
       neqmax /= neqmax_save .OR. &
       nvmax  /= nvmax_save) THEN

       IF(neqmax_save /= 0) CALL tr_neq_deallocate

    neq_allocation: DO
       ALLOCATE(nsa_neq(neqmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(nva_neq(neqmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(id_neq(neqmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(id_neqnr(neqmax,0:nrmax),STAT=ierr); IF(ierr /= 0) EXIT
       
       ALLOCATE(xv(nvmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(xv_new(nvmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(xv_prev(nvmax),STAT=ierr); IF(ierr /= 0) EXIT
       
       ALLOCATE(neqr_neq(neqmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(limtx(2,2,neqmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(rjimtx(2,2,neqmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(rsimtx(2,2,neqmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(rimtx(2,2,neqmax,neqmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(r1imtx(2,2,neqmax,neqmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(r2imtx(2,2,neqmax,neqmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(r3imtx(2,2,neqmax,neqmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(r1imtx_ofd(2,2,neqmax,neqmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(elmtx_ofd(2*neqmax,2*neqmax),STAT=ierr); IF(ierr /= 0) EXIT

       nrmax_save  = nrmax
       neqmax_save = neqmax
       nvmax_save  = nvmax
       RETURN
    END DO neq_allocation

    WRITE(6,*) 'XX tr_neq_allocate: allocation error: ierr=',ierr
    STOP

    END IF
    RETURN
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
    IF(ALLOCATED(r1imtx_ofd)) DEALLOCATE(r1imtx_ofd)
    IF(ALLOCATED(elmtx_ofd)) DEALLOCATE(elmtx_ofd)

    RETURN
  END SUBROUTINE tr_neq_deallocate


  SUBROUTINE tr_neqr_allocate
  ! ----- allocation for solver -----

    INTEGER(ikind),SAVE:: neqrmax_save=0
    INTEGER(ikind),SAVE:: nvrmax_save=0
    INTEGER(ikind)::      ierr

    IF(neqrmax /= neqrmax_save .OR. &
       nvrmax  /= nvrmax_save) THEN

       IF(neqrmax_save /= 0) CALL tr_neqr_deallocate

    neqr_allocation: DO
       ALLOCATE(neq_neqr(neqrmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(rhv(nvrmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(elmtx(2*neqrmax,2*neqrmax),STAT=ierr); IF(ierr /= 0) EXIT
       ALLOCATE(lhmtx(4*neqrmax-1,nvrmax),STAT=ierr); IF(ierr /= 0) EXIT

       neqrmax_save = neqrmax
       nvrmax_save  = nvrmax
       RETURN
    END DO neqr_allocation

    WRITE(6,*) 'XX tr_neqr_allocate: allocation error: ierr=',ierr
    STOP

    END IF
    RETURN
  END SUBROUTINE tr_neqr_allocate

  SUBROUTINE tr_neqr_deallocate

    IF(ALLOCATED(neq_neqr)) DEALLOCATE(neq_neqr)
    IF(ALLOCATED(rhv)) DEALLOCATE(rhv)
    IF(ALLOCATED(elmtx)) DEALLOCATE(elmtx)
    IF(ALLOCATED(lhmtx)) DEALLOCATE(lhmtx)

    RETURN
  END SUBROUTINE tr_neqr_deallocate


  SUBROUTINE tr_nit_allocate
  ! ----- allocation for diagnostics of iteration -----

    INTEGER(ikind),SAVE:: lmaxtr_save=0
    INTEGER(ikind)::      ierr

    ierr = 0
    IF(lmaxtr /= lmaxtr_save) THEN
       IF(lmaxtr_save /= 0) CALL tr_nit_deallocate

       nit_allocation: DO
          ALLOCATE(error_it(lmaxtr),STAT=ierr); IF(ierr /= 0) EXIT

          lmaxtr_save = lmaxtr
          RETURN
       END DO nit_allocation

       WRITE(6,*) 'XX tr_nit_allocate: allocation error: ierr=',ierr
       STOP

    END IF
    RETURN
  END SUBROUTINE tr_nit_allocate

  SUBROUTINE tr_nit_deallocate

    IF(ALLOCATED(error_it)) DEALLOCATE(error_it)

    RETURN
  END SUBROUTINE tr_nit_deallocate


  SUBROUTINE tr_ngt_allocate
  ! ----- allocation for data save -----

    INTEGER(ikind),SAVE:: ngtmax_save=0
    INTEGER(ikind),SAVE:: nsamax_save=0
    INTEGER(ikind),SAVE:: nrmax_save=0
    INTEGER(ikind)::      ierr

    IF(ngtmax /= ngtmax_save .OR. &
       nsamax /= nsamax_save .OR. &
       nrmax  /= nrmax_save ) THEN

       IF(ngtmax_save /= 0) CALL tr_ngt_deallocate

       ngt_allocation: DO
          ALLOCATE(gvt(0:ngtmax,0:50),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gvtu(0:ngtmax,0:10),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gvts(0:ngtmax,nsamax,20),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gvrt(0:nrmax,0:ngtmax,10),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(gvrts(0:nrmax,0:ngtmax,nsamax,10),STAT=ierr);IF(ierr /= 0) EXIT

          nrmax_save  = nrmax
          ngtmax_save = ngtmax
          nsamax_save = nsamax
          RETURN
       END DO ngt_allocation

    WRITE(6,*) 'XX tr_ngt_allocate: allocation error: ierr=',ierr
    STOP

    END IF
    RETURN
  END SUBROUTINE tr_ngt_allocate

  SUBROUTINE tr_ngt_deallocate

    IF(ALLOCATED(gvt)) DEALLOCATE(gvt)
    IF(ALLOCATED(gvtu)) DEALLOCATE(gvtu)
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

       nsa_allocation: DO
          ALLOCATE(ws_t(nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(stdrt(nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(offrt(nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(rns_va(nsamax),STAT=ierr); IF(ierr /= 0) EXIT
          ALLOCATE(rts_va(nsamax),STAT=ierr); IF(ierr /= 0) EXIT

          nsamax_save = nsamax
          RETURN
       END DO nsa_allocation

    WRITE(6,*) 'XX tr_nsa_allocate: allocation error: ierr=',ierr
    STOP

    END IF
    RETURN
  END SUBROUTINE tr_nsa_allocate

  SUBROUTINE tr_nsa_deallocate

    IF(ALLOCATED(ws_t)) DEALLOCATE(ws_t)
    IF(ALLOCATED(stdrt)) DEALLOCATE(stdrt)
    IF(ALLOCATED(offrt)) DEALLOCATE(offrt)
    IF(ALLOCATED(rns_va)) DEALLOCATE(rns_va)
    IF(ALLOCATED(rts_va)) DEALLOCATE(rts_va)

    RETURN
  END SUBROUTINE tr_nsa_deallocate

END MODULE trcomm
