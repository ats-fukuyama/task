module tx_commons
  implicit none
  public

  integer(4), parameter :: NQM=40, NCM=34
  integer(4), parameter :: NSM=3, NFM=2
  integer(4), parameter :: LQm1=1,  LQm2=2,  LQm3=3,  LQm4=4,  LQm5=5,&
       &                   LQe1=6,  LQe2=7,  LQe3=8,  LQe4=9,  LQe5=10,&
       &                   LQe6=11, LQe7=12, LQe8=13, &
       &                   LQi1=14, LQi2=15, LQi3=16, LQi4=17, LQi5=18,&
       &                   LQi6=19, LQi7=20, LQi8=21, &
       &                   LQz1=22, LQz2=23, LQz3=24, LQz4=25, LQz5=26,&
       &                   LQz6=27, LQz7=28, LQz8=29, &
       &                   LQb1=30, LQb2=31, LQb3=32, LQb4=33, &
       &                            LQb7=34, LQb8=35, &
       &                   LQn1=36, LQn2=37, LQn3=38, LQnz=39,&
       &                   LQr1=40
  integer(4), parameter :: nmax_file = 100

  real(4) :: gkilo = 1e3, gmega = 1e6

  !**********************************************!
  ! Set physical constants, based on CODATA 2010 !
  !    http://physics.nist.gov/cuu/Constants/    !
  !**********************************************!

  !   Electron charge (C)
  real(8), parameter :: AEE  = 1.602176565D-19

  !   Electron mass (kg)
  real(8), parameter :: AME  = 9.10938291D-31

  !   Proton mass (kg)
  real(8), parameter :: AMP  = 1.672621777D-27

  !   Speed of light in vacuum (m/s)
  real(8), parameter :: VC   = 2.99792458D8

  !   Pi (circle ratio)
  real(8), parameter :: Pi   = 3.14159265358979323846D0

  !   mu0 (H/m)
  real(8), parameter :: rMU0 = 4.D0 * PI * 1.D-7
  
  !   epsilon0 (F/m)
  real(8), parameter :: EPS0 = 1.D0 / (rMU0 * VC**2)

  !   Hydrogen ionization energy (eV)
  real(8), parameter :: EION = 13.6D0

  !   Electron-proton mass ratio (-)
  real(8), parameter :: AEP = AME / AMP

  !   Mass-to-charge ratio for a proton (kg/C)
  real(8), parameter :: AMQP = AMP / AEE

  !*******************!
  !   Set constants   !
  !*******************!


  !   Conversion factor from keV to joule
  real(8), parameter :: rKeV = 1.D3 * AEE

  !   Conversion factor from keV to eV
  real(8), parameter :: rKilo = 1.D3
  
  !   Square of pi
  real(8), parameter :: Pisq = pi * pi

  !   Square root permittivity for LQm1
  !     for the sake of acceleration of convergence
  real(8), parameter :: sqeps0 = sqrt(EPS0)

  !**********************!
  !   INPUT PARAMETERS   !
  !**********************!

  ! Configuration parameters
  real(8) :: RA, RB, ravl, rbvl, RR, BB, rbvt, rIPs, rIPe
  real(8) :: rhob, rhoaccum

  ! Species
  real(8) :: Zeffin
  real(8) :: amb, achgb
  real(8), dimension(NSM) :: amas, achg

  ! Profile shape
  real(8) :: PN0, PNa, PTe0, PTea, PTi0, PTia, PTz0, PTza
  real(8) :: PROFJ, PROFN1, PROFN2, PROFT1, PROFT2, Uiph0

  ! Diffusivities and viscosities
  real(8), dimension(NSM) :: Dfs0, rMus0, Chis0
  real(8) :: ChiNC, VWpch0, WPM0

  ! Amplitude parameters for transport
  real(8) :: FSCBAL, FSCBKP, FSCBEL, FSCBSH, FSBOHM, FSVAHL
  real(8) :: PROFD, PROFD1, PROFD2, PROFDB, PROFM, PROFM1, PROFMB, PROFC, PROFC1, PROFCB
  real(8) :: FSCX, FSLC, FSLP, FSLPB, FSION, FSNCOL, FSD01, FSD02, FSD03, FSD0z, rG1, FSRP, FSNF

  !   PROTOTYPE, not for production.
  !   MDLNZD : radial transport model for the C0-C5 superstage LQnz.
  !     0 : the hydrogen slow-neutral diffusivity D01 (current code).
  !     1 : the charged-impurity diffusivity Dfs(:,1), i.e. the bundle is
  !         transported as the magnetized C1+-C5+ population it is.
  !     2 : charge-state-weighted,  D = F_0 D_C0 + (1-F_0) Dfs(:,1),  with the
  !         effective pinch -(D_C0 - Dfs) dF_0/dV that the weighting implies.
  !     3 : as 2, with F_0 from the local sputtering-source balance.
  !     4 : as 3, plus the parallel loss of the charged part in the SOL.
  !     5 : as 4, but with the C0 density n_C0 = min(S_wall/(S_0 n_e), n_b)
  !         carried as a known field, so the flux -D_c dn_b/dV
  !         - (D_C0-D_c) dn_C0/dV and the loss are linear in n_b.
  !   MDLN0FL : Larsen flux limiter on the slow-neutral diffusivity D01,
  !     0 off, 1 on (D/(1+R) lagged), 2 on (linearized: D/(1+R)^2 plus a
  !     convective term, same flux at convergence), 3 on (D/(1+R) with R
  !     from the previous time level, held fixed during the Picard sweeps),
  !     4 on (1/(1+R) averaged with its previous-sweep value; same fixed
  !     point as 1).  D -> D / (1 + D |grad n| / (FLlim n v_th)).
  integer(4) :: MDLNZD, MDLN0FL
  !   MDLC0FL : 1 applies the Larsen limiter to D_C0 of MDLNZD = 5 even when
  !     MDLN0FL = 0 (no limiter on the hydrogen slow neutrals).
  integer(4) :: MDLC0FL
  !   RLXN0  : Picard under-relaxation weight on LQn1-LQn3 only (0 = off).
  !   FLRAMP : sweeps over which the D01 limiter is ramped to full strength.
  real(8) :: RLXN0
  integer(4) :: FLRAMP
  real(8) :: FLlim
  !   MDLSOL : width over which the parallel loss turns on outside the separatrix.
  !     0 : the historical DBW, built from half the on-axis ion temperature,
  !         Q(NRA) and B_phi (current code).
  !     1 : the ion banana width at the separatrix, sqrt(eps) rho_theta, from the
  !         local edge T_i and B_theta.  Parameter-free.
  integer(4) :: MDLSOL
  !   MDLWBC : density boundary condition at the wall rho = rhob.
  !     0 : zero gradient imposed only weakly (current code).
  !     1 : n(NRMAX) = n(NRMAX-1), the zero gradient imposed explicitly
  !         (Dirichlet value from the previous Picard sweep).
  !     2 : the same condition as a row of the system, solved implicitly.
  !     3 : as 2, with the wall row of the continuity equation first added
  !         to row NRMAX-1 (constrained-DOF Galerkin form).
  integer(4) :: MDLWBC
  real(8) :: FSADV, FSADVB, FSUG, FSG0iz
  real(8), dimension(NSM) :: FSMPCH, FSTPTM, FSPARV, FSLTs
  real(8), dimension(1:2) :: FSNC, FSNCB
  real(8), dimension(1:3) :: FSDFIX, FSANOM, RhoETB, FSPCL
  integer(4) :: MDLC, MDANOM

  ! SOL parameters
  real(8) :: rLn, rLT  ! Scale lengths in SOL
  real(8), dimension(NSM) :: PNsSUP, PTsSUP

  ! Heat sources
  real(8) :: Ebmax, RNBP, RNBP0, RNBT1, RNBT2, RNBT10, RNBT20, &
       &     PNBH, PNBHP, PNBHT1, PNBHT2, PNBHex, PNBCD, PNBMPD, PNBPTC, &
       &     rNRFe, RRFew, RRFe0, PRFHe, rNRFi, RRFiw, RRFi0, PRFHi, Tqt0, Tqp0
  real(8), dimension(1:3) :: esps
  integer(4) :: MDLNBD, MDBMCX

  ! Neutral parameters
  real(8) :: rGamm0_in, rGamm0zz_in, rGamm0iz_in
  real(8) :: PN0s, V0, rGamm0, rGASPF
  real(8) :: PN0zs, V0z, rGamm0zz, rGamm0iz, rGASPFz

  ! Ripple parameters
  real(8) :: DltRPn, kappa
  integer(4) :: NTCOIL, m_pol, n_tor

  ! Numerical parameters
  real(8) :: DT, EPS, ADV, tiny_cap, CMESH0, WMESH0, CMESH, WMESH

  !   Nonlinear-iteration index up to which the lagged quantities that feed the
  !   transport coefficients -- PNsV_FIX, PTsV_FIX, ErV_FIX, SCX and ddPhidpsi --
  !   are refreshed.  Beyond it they stay frozen at the iterate they had, which
  !   is an O(DT) distance from the converged state and therefore costs the time
  !   integration its second order however good the scheme is.  Freezing is there
  !   to steady the nonlinear iteration and to keep an expensive transport model
  !   from being re-evaluated every sweep; raise it to ICMAX where accuracy in
  !   time matters and the model is cheap.
  integer(4) :: ICCALC
  integer(4) :: ICMAX

  ! Mesh parameters
  integer(4) :: NRMAX, NEMAX, NRA, NRC

  ! Equilibrium parameters
  integer(4) :: ieqread, irktrc, ipbtdir

  ! Parameters for density perturbation experiment
  real(8) :: DelRho, DelN

  ! Helical parameters
  real(8) :: EpsH, FSHL, Q0, QA
  integer(4) :: NCph, NCth
  !!  multiple Fourier mode      miki_m 10-08-10
!  integer(4) :: NHFMmx
!  integer(4), dimension(:,:), allocatable :: HPN   ! helical pitch numbers
!  real(8), dimension(:,:), allocatable :: EpsHM   ! helical amplitude
  integer(4), parameter :: NHFMmx = 20   !! max value of helical Fourier modes
  integer(4), dimension(1:NHFMmx, 1:2) :: HPN   ! helical pitch numbers
  real(8), dimension(1:NHFMmx, 0:3) :: EpsHM   ! helical amplitudes for rho^0:rho^3

  ! Magnetic braiding parameters
  real(8) :: DMAG0, RMAGMN, RMAGMX

  ! Diagnostics
  integer(4) :: IDIAG, NTSTEP, midbg(2)

  ! LAPACK
  integer(4) :: MDLPCK

  ! Numerical parameteres
  !   Implicit Euler or Gear's backward differential formula
  integer(4) :: IGBDF
  !   Scheme to evaluate dPhi/dV and dPs/dpsi for suppressing wiggles near the axis
  integer(4) :: ISMTHD

  ! Choice of equations
  integer(4) :: MDFIXT, MDBEAM

  ! Transport model
  integer(4) :: MDOSQZ, MDOSQZN, MDLETA, MDLETB, MDLNEO, MDBSETA, imodel_neo(6)

  ! Initial condition
  integer(4) :: MDITSN, MDITST, MDINTN, MDINTT, MDINTC

  ! Poynting flux calculation
  integer(4) :: iPoyntpol, iPoynttor

  !*********************************!
  !   INTERNAL CONTROL PARAMETERS   !
  !*********************************!

  ! Initializer of arrays
  real(8), dimension(:),   allocatable :: array_init_NR
  real(8), dimension(:,:), allocatable :: array_init_NRNS

  ! Time control parameters
  integer(4) :: NT, NTCUM, NTMAX
  real(8) :: T_TX = 0.D0, TMAX

  ! Configuration parameters
  integer(4) :: NQMAX, IERR, ICONT, IRPIN, irestart

  !   Set once standard input has reached end of file.  After the first EOF a
  !   read may report a positive iostat instead of a negative one (gfortran
  !   returns 5001, "read past end of file"), and every interactive loop here
  !   treats a positive iostat as a typo and asks again -- which spins forever
  !   on a redirected or closed stdin.  Each loop consults this first.
  logical :: stdin_ended = .false.
  real(4) :: AVE_IC
  real(8) :: rIP
  real(8) :: UHth, UHph
  real(8) :: Rax, Zax, surflcfs
  real(8), dimension(:), allocatable :: vv, rho
  real(8), dimension(:), allocatable :: hv

  ! Convergence accelerator
  real(8) :: rMUb1, rMUb2

  ! Convergence parameter
  integer(4) :: MODECV, iprestab
  real(8) :: oldmix
  ! Adaptive under-relaxation factor, engaged only on time steps where the
  ! Picard sweep is measured not to contract.  See TXLOOP.  0 disables it.
  real(8) :: oldmixa

  ! Diagnostic only: project the near-axis profiles onto the regular (even in rho)
  ! subspace q0 + q2 rho^2 + q4 rho^4 after each time step.  Default off.
  integer(4) :: iaxproj
  real(8) :: rhoaxproj

  ! Diagnostic only: a conservative particle diffusion added to the continuity
  ! equations, as a fraction of each species' own momentum diffusivity.  This is
  ! a regularity-selection experiment, not a model term.  0 = off (default).
  real(8) :: Dnfac

  ! SUPG parameter
  integer(4) :: iSUPG3, iSUPG6, iSUPG8
  real(8) :: SUPGstab

  ! Variables
  real(8), dimension(:), allocatable :: &
       & ErV,    BEpol,   Etor,    BthV,  BphV, &
       & PNbV,   UbrV,    UbrVV,   BUbparV, RUbphV, PTbV, UbphVR, UbphV, BVbdiag, &
       & PN01V,  PN02V,   PN03V,   PN0zV,&
       & PsiV,   PsitV,   PhiV,  &
       & PT01V,  PT02V,   PT03V,   PNbrpV, &
       & PsidotV,PsitdotV,dPhidV,  dPhidrho, bbt, bthco, Zeff
!!rp_conv       &, PNbrpLV

  real(8), dimension(:), allocatable :: ErV_FIX

  ! Species variables
  type species_var
     real(8) :: n      = 0.1d0 ! [10^{20}m^{-3}]
     real(8) :: T      = 1.d0  ! [keV]
     real(8) :: p      = 0.1d0 ! [10^{20}m^{-3}*keV]
     real(8) :: UrV    = 0.d0  ! [1/s]    : n <u.grad V>
     real(8) :: BUpar  = 0.d0  ! [Tm/s]
     real(8) :: Bqpar  = 0.d0  ! [Tm/s]   : <Bq//>/(5p)
     real(8) :: Uthhat = 0.d0  ! [m/s/T]  : <u.grad theta>/<B.grad theta>
     real(8) :: UphR   = 0.d0  ! [/s]
     real(8) :: RUph   = 0.d0  ! [m^2/s]
     real(8) :: Ur     = 0.d0  ! [m/s]
     real(8) :: Uth    = 0.d0  ! [m/s]
     real(8) :: Uph    = 0.d0  ! [m/s]
     real(8) :: BV1    = 0.d0  ! [Tm/s]
  end type species_var
  type(species_var), dimension(:,:), allocatable :: Var

  real(8), dimension(:,:), allocatable :: PNsV_FIX, PTsV_FIX

  real(8), dimension(:,:,:), allocatable :: BVsdiag, BusparNCL, UsthNCL, QsthNCL

  ! Coefficients
  real(8), dimension(:), allocatable :: &
       & rNuION, rNu0b, rNuL, rNuiCX, rNubCX, rNuiCXT, &
       & rNuei, rNuii, rNuiz, rNuTei, rNuTez, rNuTiz, rNuD, rNuOL, &
       & WPM, FVpch, FQLcoef, FQLcoef1, FQLcoef2, FQCcoef, &
       & rNuB, rNuLB, ft, VWpch, D01, D02, D03, D0z, Dnzb, Vnzb, FracC0, nC0z, DC0z, W01, &
       & SiVizA, SiVcxA, SiVa6A, SiVsefA, wexb, &
       & UgV, PNbVinv, Vbpara
  real(8), dimension(:), allocatable :: & ! Ripple transport
       & Ubrp, RUbrp, Dbrp, DltRP, DltRP_mid, rNubL, rip_rat, &
       & rNubrp1, rNubrp2
  real(8), dimension(:), allocatable :: & ! Quantities related to Helicals 
       & rNueHLthth,rNueHLthph, rNueHLphth, rNueHLphph, &
       & rNuiHLthth,rNuiHLthph, rNuiHLphth, rNuiHLphph, &
       & DMAG, DMAGe, DMAGi
  real(8), dimension(:,:), allocatable :: Dfs, rMus, Chis, ChiNCp, ChiNCt
  real(8), dimension(:,:), allocatable :: rNu0s, rNuLTs, rNuAss, Vhps, Vmps, PiRess
  real(8), dimension(:,:), allocatable :: gamITG

  ! Coefficients related to neoclassical transport
  integer(4), dimension(:),allocatable :: mxneo
  real(8), dimension(:),   allocatable :: gamneo
  real(8), dimension(:,:), allocatable :: fmneo, xmuf, BnablaPi, cL31
  real(8), dimension(:,:,:), allocatable :: lff, gflux
  real(8), dimension(:,:,:,:), allocatable :: xmu, laf, lfb
  real(8), dimension(:,:,:,:,:), allocatable :: lab

  ! Equilibrium metrics 
  real(8), dimension(:), allocatable :: &
       & aat, rrt, ckt, suft, sst, vro, vlt, rhov, art, epst, ait, elip, trig, rtt, rpt, lpt, drhodr
  real(8), dimension(:), allocatable :: fipol, Bpsq, qhatsq, Fqhatsq, BEpara, bri
  real(8), dimension(:), allocatable :: sdt, hdt
  real(8), dimension(:), allocatable :: bit, bbrt
  real(8), dimension(:), allocatable :: gtti

  ! Multiple helical Fourier modes    miki_m 10-08-06~
  integer(4) :: UHphSwitch           !  To support (m=0, n>0) Fourier component
  real(8), dimension(:,:), allocatable :: rNueHLM, rNuiHLM, &
       & rNueHLththM, rNueHLthphM, rNueHLphthM, rNueHLphphM,&
       & rNuiHLththM, rNuiHLthphM, rNuiHLphthM, rNuiHLphphM 

  ! CDBM
  real(8), dimension(:), allocatable :: rG1h2, FCDBM, S, Alpha, rKappa

  ! For numerical stability
  real(8), dimension(:), allocatable :: pres0, ErV0

  ! CDIM
  ! 09/06/17~ miki_m
  real(8), dimension(:), allocatable :: rG1h2IM, FCDIM, RAQPR

  ! Sources and sinks
  real(8), dimension(:), allocatable :: PNB, PNBTG, PNBPD, PNBcol_e, PNBcol_i, &
       &                                SNB, SNBe, SNBi, SNBb, SNBPDi, SNBTGi, &
       &                                PNBe, PNBi,PNBz, MNB, PRFe, PRFi, PRFz, &
       &                                POH, PEQei, PEQez, PEQiz, &
       &                                SiLC, SiLCB, SiLCph, PALFe, PALFi, PALFz, &
       &                                BSmb, Tqt, Tqp
  real(8), dimension(:), allocatable :: PIE, PCX, SIE, SCX, PBr, RatCX, ratecxb
  real(8), dimension(:,:), allocatable :: POHs
  real(8) :: Eb, Vb

  ! Safety factor, currents, resistivity
  real(8), dimension(:), allocatable :: Q, AJ, BJPARA, AJOH, BJOH, AJRF, BJRF, AJNB, BJNB, &
       &                                AJBS, BJBS, ETA, ETAS, AJV
  real(8), dimension(:,:), allocatable :: BJBSvar, ETAvar

  ! Derivatives
  real(8), dimension(:,:), allocatable :: dPsdpsi, dTsdpsi

  ! For display
  real(8), dimension(:), allocatable :: ANS0, TS0, ANSAV, TSAV, WST
  real(8), dimension(:), allocatable :: ANF0, TF0, ANFAV, TFAV, WFT
  real(8), dimension(:), allocatable :: PBCLT, PFCLT, PLT, SPET, SLT
  real(8), dimension(:), allocatable :: PoyntS, PoyntI, PoyntR, CPsi, CPsi_old, VPoynt
  real(8), dimension(:), allocatable :: thrp, qneut
  real(8), dimension(:), allocatable :: PnumN0
  real(8), dimension(:,:), allocatable :: Deff
  real(8) :: WBULKT, WTAILT, WPT
  real(8) :: AJT, AJOHT, AJNBT, AJRFT, AJBST
  real(8) :: PINT, POHT, PNBT, PRFT, PRFTe, PRFTi, PNFT
  real(8) :: PBINT, PFINT, POUT, PCXT, PIET, PRLT, SINT, SIET
  real(8) :: SNBT, SNFT, SOUT, TNBcol, TTqt
  real(8) :: VLOOP, VLOOPpol, ALI, RQ1, RPE, ZEFF0, QF
  real(8) :: WPDOT, TAUE1, TAUE2, TAUEP, TAUEH, TAUP, TAUPA
  real(8) :: BETAP0, BETAPA, BETA0, BETAA, BETAQ0, BETAN
  real(8) :: TPRE, WPPRE, totmnRV, totjxB
  real(8) :: VOLAVN, Gamma_a, CEjima, PnumN0z

  ! Internal variables for transport matrix
  real(8),    dimension(:,:,:), allocatable :: ALC, BLC, CLC, PLC
  integer(4), dimension(:,:),   allocatable :: NLC
  integer(4), dimension(:,:,:), allocatable :: NLCR
  integer(4), dimension(:),     allocatable :: NLCMAX
  real(8),    dimension(:,:),   allocatable :: X, XOLD

  ! .false. while the two-step BDF has no genuine X(t-dt) to lean on: at a cold
  ! start, and after a LOAD because XOLD does not travel in the save file.  A
  ! CONT does have one and keeps it.  See TXCALB.
  logical :: bdf_history_ready = .false.
  ! Step size the stored XOLD belongs to.  The two-step BDF weights are the
  ! fixed-step ones, so the history is only usable while DT is unchanged.
  real(8) :: bdf_history_dt = 0.d0

  ! Diagnostics
  integer(4) :: MODEAV

  ! I/O
  character(len=20) :: SLID

  ! Input data from ascii files
  type infiles_type
     character(len=6) :: name
     integer(4) :: nol ! number of lines, i.e. upper bound of the array
     integer(4) :: ncol_mesh, ncol_data ! indices indicating the columns where the mesh
                                        !   and the data are stored
     real(8) :: totS, totP ! Total number of ions per second and total power
     real(8), dimension(1:nmax_file) :: r, data, vb
!!     real(8), dimension(:), pointer :: r, data, vb
  end type infiles_type
  type(infiles_type), allocatable :: infiles(:)
  integer(4) :: n_infiles = 0 ! number of data which are read from the file
  integer(4) :: iflag_file = 0 ! Flag indicating whether infiles is null or not.
                               ! Default "null"
  character(len=6), dimension(1:5) :: datatype
  data datatype /'PNBP', 'PNBT1', 'PNBT2', 'PRF', 'LQe4'/

contains

  subroutine allocate_txcomm(ier, icont_in)

    integer(4), intent(out) :: ier
    integer(4), intent(in), optional :: icont_in
    integer(4) :: iflag, NS, NF
    integer(4), dimension(1:20) :: ierl

    ierl(:) = 0
    if(nrmax <= 1) then
      write(6,*) "XXX ALLOCATE_TXCOMM : ILLEGAL PARAMETER    NRMAX=", nrmax
      ier = 1
      return
    endif

    NEMAX = NRMAX
    NS    = NSM
    NF    = NFM

    ! allocation check
    if(allocated(array_init_NR)) then
       if(present(icont_in) .and. icont_in /= 0) then
          call deallocate_txcomm
       else
          return
       end if
    end if

    do
       ! Array for numerical purpose
       allocate(array_init_NR(0:NRMAX),         source=0.d0,  stat=ierl(1))
       allocate(array_init_NRNS(0:NRMAX,1:NSM), source=0.d0,  stat=ierl(2))

       ! Mesh-related arrays
       allocate(vv, rho,                source=array_init_NR, stat=ierl(3))
       allocate(hv(0:NEMAX),                                  stat=ierl(4))
       ier = sum(ierl) ; iflag = 1
       if (ier == 0) then ; ierl(:) = 0 ; else ; exit ; end if

       ! Arrays associated with dependent variables
       allocate(ErV,     BEpol,    Etor,    BthV,   BphV,                 source=array_init_NR, stat=ierl(1))
       allocate(PNbV,    UbrV,     UbrVV,   BUbparV, RUbphV, PTbV,  UbphVR, UbphV, BVbdiag, source=array_init_NR, stat=ierl(2))
       allocate(PN01V,   PN02V,    PN03V,   PN0zV,                        source=array_init_NR, stat=ierl(3))
       allocate(PsiV,    PsitV,    PhiV,    bbt,    bthco, Zeff,          source=array_init_NR, stat=ierl(4))
       allocate(PT01V,   PT02V,    PT03V,   PNbrpV,                       source=array_init_NR, stat=ierl(5))
       allocate(PsidotV, PsitdotV, dPhidV,  dPhidrho,                     source=array_init_NR, stat=ierl(6))
       ! Derived-type dependent variables
       allocate(Var(0:NRMAX,NS),                                                                stat=ierl(7))
       ier = sum(ierl) ; iflag = 2
       if (ier == 0) then ; ierl(:) = 0 ; else ; exit ; end if

       ! For numerical stabilization
       allocate(ErV_FIX,  pres0,    ErV0,                               source=array_init_NR,   stat=ierl(1))
       allocate(PNsV_FIX, PTsV_FIX,                                     source=array_init_NRNS, stat=ierl(2))
       ier = sum(ierl) ; iflag = 3
       if (ier == 0) then ; ierl(:) = 0 ; else ; exit ; end if

       ! Physical coefficients
       allocate(rNuION,  rNu0b,   rNuL,    rNuiCX,  rNubCX,  rNuiCXT,   source=array_init_NR,   stat=ierl(1))
       allocate(rNuei,   rNuii,   rNuiz,   rNuTei,  rNuTez,  rNuTiz,    source=array_init_NR,   stat=ierl(2))
       allocate(rNuD,    rNuOL,                                         source=array_init_NR,   stat=ierl(3))
       allocate(WPM,     FVpch,   FQLcoef, FQLcoef1, FQLcoef2, FQCcoef, source=array_init_NR,   stat=ierl(4))
       allocate(rNuB,    rNuLB,   ft,      VWpch,                       source=array_init_NR,   stat=ierl(5))
       allocate(D01,     D02,     D03,     D0z,                         source=array_init_NR,   stat=ierl(6))
       allocate(Dnzb,    Vnzb,    FracC0,  nC0z, DC0z, W01,                             source=array_init_NR,   stat=ierl(6))
       allocate(SiVizA,  SiVcxA,  SiVa6A,  SiVsefA,  wexb,              source=array_init_NR,   stat=ierl(7))
       allocate(UgV,     PNbVinv, Vbpara,                               source=array_init_NR,   stat=ierl(8))
       allocate(Ubrp,    RUbrp,   Dbrp,    DltRP,    DltRP_mid,         source=array_init_NR,   stat=ierl(9))
       allocate(rNubL,   rip_rat, rNubrp1, rNubrp2,                     source=array_init_NR,   stat=ierl(10))
       allocate(Dfs,     rMus,    Chis,    ChiNCp,   ChiNCt,            source=array_init_NRNS, stat=ierl(11))
       allocate(rNu0s,   rNuLTs,  rNuAss,  Vhps,     Vmps,     PiRess,  source=array_init_NRNS, stat=ierl(12))
       allocate(gamITG(0:NRMAX,1:3),                                                            stat=ierl(13))
       ier = sum(ierl) ; iflag = 4
       if (ier == 0) then ; ierl(:) = 0 ; else ; exit ; end if

       ! Neoclassical transport quantities
       allocate(xmu(0:NRMAX,NS,3,3), lab(0:NRMAX,NS,NS,3,3),                           stat=ierl(1))
       allocate(laf(0:NRMAX,NS,2,2), lfb(0:NRMAX,NS,2,2), lff(0:NRMAX,2,2),            stat=ierl(2))
       allocate(xmuf(0:NRMAX,3),                                                       stat=ierl(3))
       allocate(BnablaPi, cL31,                               source=array_init_NRNS,  stat=ierl(4))
       allocate(gflux(0:NRMAX,NS,3),                                                   stat=ierl(5))
       allocate(mxneo(0:NRMAX),      fmneo(1:10,0:NRMAX), gamneo(0:NRMAX),             stat=ierl(6))
       allocate(BVsdiag(0:NRMAX,NS,2), BusparNCL(0:NRMAX,NS,2), UsthNCL(0:NRMAX,NS,2), QsthNCL(0:NRMAX,NS,2), stat=ierl(7))
       ier = sum(ierl) ; iflag = 5
       if (ier == 0) then ; ierl(:) = 0 ; else ; exit ; end if

       ! Equilibrium quantities
       allocate(aat, rrt, ckt, suft, sst, vro, vlt, rhov, art, &
            &   epst, ait, elip, trig, rtt, rpt, lpt, drhodr,   source=array_init_NR, stat=ierl(1))
       allocate(fipol, Bpsq, qhatsq, Fqhatsq, BEpara, bri,      source=array_init_NR, stat=ierl(2))
       allocate(sdt, hdt, bit, bbrt, gtti,                      source=array_init_NR, stat=ierl(3))
       ier = sum(ierl) ; iflag = 6
       if (ier == 0) then ; ierl(:) = 0 ; else ; exit ; end if

       ! CDBM and CDIM
       allocate(rG1h2,   FCDBM,  S,     Alpha, rKappa,          source=array_init_NR, stat=ierl(1))
       allocate(rG1h2IM, FCDIM,  RAQPR,                         source=array_init_NR, stat=ierl(2)) !09/06/17 miki_m
       ier = sum(ierl) ; iflag = 7
       if (ier == 0) then ; ierl(:) = 0 ; else ; exit ; end if

       ! Sources and sinks
       allocate(PNB,    PNBTG,   PNBPD,  PNBcol_e, PNBcol_i,        source=array_init_NR,   stat=ierl(1))
       allocate(SNB,    SNBe,    SNBi,   SNBb,                      source=array_init_NR,   stat=ierl(2))
       allocate(SNBPDi, SNBTGi,  MNB,                               source=array_init_NR,   stat=ierl(3))
       allocate(PNBe,   PNBi,    PNBz,   PRFe,     PRFi,     PRFz,  source=array_init_NR,   stat=ierl(4))
       allocate(POH,    PEQei,   PEQez,  PEQiz,                     source=array_init_NR,   stat=ierl(5))
       allocate(SiLC,   SiLCB,   SiLCph, PALFe,    PALFi,    PALFz, source=array_init_NR,   stat=ierl(6))
       allocate(PIE,    PCX,     SIE,    SCX,      PBr,             source=array_init_NR,   stat=ierl(7))
       allocate(BSmb,   Tqt,     Tqp,                               source=array_init_NR,   stat=ierl(8))
       allocate(RatCX,  ratecxb,                                    source=array_init_NR,   stat=ierl(9))
       allocate(POHs,                                               source=array_init_NRNS, stat=ierl(10))
       ier = sum(ierl) ; iflag = 8
       if (ier == 0) then ; ierl(:) = 0 ; else ; exit ; end if

       ! Currents and the safety factor
       allocate(Q,    AJ,   BJPARA, AJOH, BJOH, AJRF, BJRF,     source=array_init_NR, stat=ierl(1))
       allocate(AJNB, BJNB, AJBS,   BJBS,                       source=array_init_NR, stat=ierl(2))
       allocate(ETA,  ETAS, AJV,                                source=array_init_NR, stat=ierl(3))
       allocate(BJBSvar(0:NRMAX,0:3), ETAvar(0:NRMAX,0:4),                            stat=ierl(4))
       ier = sum(ierl) ; iflag = 9
       if (ier == 0) then ; ierl(:) = 0 ; else ; exit ; end if

       ! Derivatives with respect to psi
       allocate(dPsdpsi, dTsdpsi,                             source=array_init_NRNS, stat=ierl(1))
       ier = sum(ierl) ; iflag = 10
       if (ier == 0) then ; ierl(:) = 0 ; else ; exit ; end if

       ! Global quantities mainly for statistics and graphics
       allocate(ANS0(1:NS), TS0(1:NS),   ANSAV(1:NS), TSAV(1:NS), WST(1:NS),          stat=ierl(1))
       allocate(ANF0(1:NF), TF0(1:NF),   ANFAV(1:NF), TFAV(1:NF), WFT(1:NF),          stat=ierl(2))
       allocate(PBCLT(1:NS),PFCLT(1:NS), PLT(1:NS),   SPET(1:NS), SLT(1:NS),          stat=ierl(3))
       allocate(PoyntS(1:2),PoyntI(1:4), PoyntR(1:2), &
            &   CPsi(0:3),  CPsi_old(0:3),VPoynt(0:3),                                stat=ierl(4))
       allocate(thrp(1:2*NRMAX),qneut(0:NRMAX),                                       stat=ierl(5))
       allocate(PnumN0(0:3),                                                          stat=ierl(6))
       allocate(Deff,                                         source=array_init_NRNS, stat=ierl(7))
       ier = sum(ierl) ; iflag = 11
       if (ier == 0) then ; ierl(:) = 0 ; else ; exit ; end if

       ! Internal arrays for numerics
       allocate(ALC(0:NRMAX,0:NCM,1:NQMAX),                                           stat=ierl(1))
       allocate(BLC, CLC,                                                source=ALC,  stat=ierl(2))
       allocate(PLC(0:NRMAX,1:NCM,1:NQMAX),                                           stat=ierl(3))
       allocate(NLC(0:NCM,1:NQMAX), NLCR(0:NCM,0:NQMAX,0:1), NLCMAX(1:NQMAX),         stat=ierl(4))
       allocate(X(0:NRMAX,1:NQMAX),                                      source=0.d0, stat=ierl(5))
       allocate(XOLD,                                                    source=X,    stat=ierl(6))
       ier = sum(ierl) ; iflag = 12
       if (ier == 0) then ; ierl(:) = 0 ; else ; exit ; end if

       ! Coefficients associated with non-axisymmetric transport
!       allocate(rNueHL, rNuiHL, rNueHLthth, rNuiHLthth, rNueHLthph, &
!            &   rNuiHLthph, rNueHLphth, rNuiHLphth, rNueHLphph, &
!            &   rNuiHLphph,                             source=array_init_NR, stat=ierl(1))
!       allocate(HPN(1:2,1:NHFMmx),EpsHM(1:2,1:NHFMmx),                        stat=ierl(2))
       allocate(rNueHLthth, rNuiHLthth, rNueHLthph, rNuiHLthph, rNueHLphth, rNuiHLphth, &
            &   rNueHLphph, rNuiHLphph,                         source=array_init_NR, stat=ierl(1))
       allocate(rNueHLM(1:NHFMmx,0:NRMAX),                      source=0.d0,          stat=ierl(2))
       allocate(rNuiHLM, rNueHLththM, rNueHLthphM, rNueHLphthM, rNueHLphphM, rNuiHLththM, &
            &   rNuiHLthphM, rNuiHLphthM, rNuiHLphphM,          source=rNueHLM,       stat=ierl(3))
       allocate(DMAG,    DMAGe,   DMAGi,                        source=array_init_NR, stat=ierl(4))
       ier = sum(ierl) ; iflag = 13
       if (ier /= 0) exit

       ! normal end
       iflag = 0
       exit
    end do

    ! All the memories allocated above are clear if some errors occur.
    if(iflag /= 0) then
       write(6,*) "XX Allocation error category = ",iflag
       call deallocate_txcomm
    end if

  end subroutine allocate_txcomm

  subroutine deallocate_txcomm

    if( .not. allocated(array_init_NR) ) return

    deallocate(array_init_NR)
    deallocate(vv,      rho)
    deallocate(hv)

    deallocate(ErV,    BEpol,   Etor,   BthV,  BphV)
    deallocate(PNbV,   UbrV,    UbrVV,  BUbparV,RUbphV, PTbV, UbphVR, UbphV, BVbdiag)
    deallocate(PN01V,  PN02V,   PN03V,  PN0zV)
    deallocate(PsiV,   PsitV,   PhiV)
    deallocate(PT01V,  PT02V,   PT03V,  PNbrpV)
    deallocate(PsidotV,PsitdotV,dPhidV, dPhidrho, bbt, bthco, Zeff)
    deallocate(Var)

    deallocate(ErV_FIX, pres0,  ErV0)
    deallocate(PNsV_FIX, PTsV_FIX)

    deallocate(rNuION, rNu0b,  rNuL,   rNuiCX,  rNubCX,  rNuiCXT)
    deallocate(rNuei,  rNuii,  rNuiz,  rNuTei,  rNuTez,  rNuTiz)
    deallocate(rNuOL,  rNuD)
    deallocate(WPM,    FVpch,  FQLcoef,FQLcoef1,FQLcoef2,FQCcoef)
    deallocate(rNuB,   rNuLB,  ft,     VWpch)
    deallocate(D01,    D02,    D03,    D0z)
    deallocate(Dnzb,   Vnzb,   FracC0, nC0z, DC0z, W01)
    deallocate(UgV,    PNbVinv,Vbpara)
    deallocate(Ubrp,   RUbrp,  Dbrp,   DltRP,   DltRP_mid)
    deallocate(rNubL,  rip_rat,rNubrp1,rNubrp2)
    deallocate(Dfs,    rMus,   Chis,   ChiNCp,  ChiNCt) 
    deallocate(rNu0s,  rNuLTs, rNuAss, Vhps,    Vmps,   PiRess)
    deallocate(gamITG)
    deallocate(SiVizA, SiVcxA, SiVa6A, SiVsefA, wexb)

    deallocate(xmu, xmuf, lab, laf, lfb, lff, BnablaPi, gflux, cL31, mxneo, fmneo, gamneo)
    deallocate(BVsdiag, BusparNCL, UsthNCL, QsthNCL)

    deallocate(aat, rrt, ckt, suft, sst, vro, vlt, rhov, art, epst, ait, elip, trig, rtt, rpt, lpt, drhodr)
    deallocate(fipol, Bpsq, qhatsq, Fqhatsq, BEpara, bri)
    deallocate(sdt, hdt, bit, bbrt, gtti)

    deallocate(rG1h2,  FCDBM,  S,     Alpha, rKappa)
    deallocate(rG1h2IM, FCDIM, RAQPR)  !***miki_m 09/06/17~

    deallocate(PNB,   PNBTG, PNBPD, PNBcol_e, PNBcol_i)
    deallocate(SNB,   SNBe,  SNBi,  SNBb,   SNBPDi, SNBTGi,  MNB)
    deallocate(PNBe,  PNBi,  PNBz,  PRFe,   PRFi,   PRFz)
    deallocate(POH,   PEQei, PEQez, PEQiz)
    deallocate(SiLC,  SiLCB, SiLCph,PALFe,  PALFi,  PALFz)
    deallocate(PIE,   PCX,   SIE,   SCX,    PBr)
    deallocate(BSmb,  Tqt,   Tqp)
    deallocate(RatCX, ratecxb)
    deallocate(POHs)

    deallocate(Q, AJ, BJPARA, AJOH, BJOH, AJRF, BJRF)
    deallocate(AJNB, BJNB, AJBS, BJBS)
    deallocate(ETA, ETAS, AJV)
    deallocate(BJBSvar, ETAvar)

    deallocate(dPsdpsi, dTsdpsi)

    deallocate(ANS0, TS0, ANSAV, TSAV, WST)
    deallocate(ANF0, TF0, ANFAV, TFAV, WFT)
    deallocate(PBCLT, PFCLT, PLT, SPET, SLT)
    deallocate(PoyntS, PoyntI, PoyntR, CPsi, CPsi_old, VPoynt)
    deallocate(thrp, qneut)
    deallocate(PnumN0)
    deallocate(Deff)

    deallocate(ALC, BLC, CLC)
    deallocate(PLC)
    deallocate(NLC,NLCR,NLCMAX)
    deallocate(X, XOLD)

!    deallocate(rNueHL, rNuiHL, rNueHLthth, rNuiHLthth, rNueHLthph, &
!         &     rNuiHLthph, rNueHLphth, rNuiHLphth, rNueHLphph, &
!         &     rNuiHLphph)
!    deallocate(HPN, EpsHM) ! miki_m 10-08-06
    deallocate(rNueHLM, rNuiHLM, &
         &     rNueHLthth, rNuiHLthth, rNueHLthph, rNuiHLthph, &
         &     rNueHLphth, rNuiHLphth, rNueHLphph, rNuiHLphph, &
         &     rNueHLththM, rNueHLthphM, rNueHLphthM, rNueHLphphM, &
         &     rNuiHLththM, rNuiHLthphM, rNuiHLphthM, rNuiHLphphM) ! miki_m 10-08-06
    deallocate(DMAG, DMAGe, DMAGi)  !***AF (2008-06-08)

    if(allocated(infiles)) deallocate(infiles)

  end subroutine deallocate_txcomm

end module tx_commons
