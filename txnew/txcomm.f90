module tx_commons
  implicit none
  public

  integer(4), parameter :: NQM=28, NCM=28
  integer(4), parameter :: NSM=2, NFM=2
  integer(4), parameter :: LQm1=1,  LQm2=2,  LQm3=3,  LQm4=4,  LQm5=5,&
       &                   LQe1=6,  LQe2=7,  LQe3=8,  LQe4=9,  LQe5=10,&
       &                   LQe6=11, LQe7=12, &
       &                   LQi1=13, LQi2=14, LQi3=15, LQi4=16, LQi5=17,&
       &                   LQi6=18, LQi7=19, &
       &                   LQb1=20, LQb2=21, LQb3=22, LQb4=23, &
       &                            LQb7=24, &
       &                   LQn1=25, LQn2=26, LQn3=27, &
       &                   LQr1=28
  integer(4), parameter :: nmax_file = 100

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
  real(8), parameter :: PI   = 3.14159265358979323846D0

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
  real(8) :: RA, RB, RR, BB, rIPs, rIPe
  real(8) :: rhob, rhoaccum

  ! Species
  real(8) :: Zeff
  real(8) :: amb, achgb
  real(8), dimension(NSM) :: amas, achg

  ! Profiles
  real(8) :: PN0, PNa, PTe0, PTea, PTi0, PTia, PROFJ, PROFN1, PROFN2, PROFT1, PROFT2, Uiph0

  ! Diffusivities and viscosities
  real(8) :: De0, Di0, rMue0, rMui0, Chie0, Chii0, ChiNC, VWpch0, WPM0

  ! Amplitude parameters for transport
  real(8) :: FSCBAL, FSCBKP, FSCBEL, FSCBSH, FSBOHM, FSPCLD, FSPCLM, FSPCLC, FSVAHL
  real(8) :: PROFD, PROFD1, PROFD2, PROFDB, PROFM, PROFM1, PROFMB, PROFC, PROFC1, PROFCB
  real(8) :: FSCX, FSLC, FSNC, FSNCB, FSLP, FSLTE, FSLTI, FSION, FSD01, FSD02, FSD03, rG1, FSRP, FSNF
  real(8) :: FSADV, FSADVB, FSUG
  real(8), dimension(1:NSM) :: FSMPCH, FSPARV
  real(8), dimension(1:3) :: FSDFIX, FSANOM, RhoETB
  integer(4) :: MDLC, MDANOM

  ! Scale lengths in SOL
  real(8) :: rLn, rLT

  ! Heat sources
  real(8) :: Ebmax, RNBP, RNBP0, RNBT1, RNBT2, RNBT10, RNBT20, &
       &     PNBH, PNBHP, PNBHT1, PNBHT2, PNBHex, PNBCD, PNBMPD, PNBPTC, &
       &     rNRFe, RRFew, RRFe0, PRFHe, rNRFi, RRFiw, RRFi0, PRFHi, Tqt0, Tqp0
  real(8), dimension(1:3) :: esps
  integer(4) :: MDLNBD

  ! Neutral parameters
  real(8) :: PN0s, V0, rGamm0, rGASPF, PNeDIV, PNiDIV, PTeDIV, PTiDIV

  ! Ripple parameters
  real(8) :: DltRPn, kappa
  integer(4) :: NTCOIL, m_pol, n_tor

  ! Numerical parameters
  real(8) :: DT, EPS, ADV, tiny_cap, CMESH0, WMESH0, CMESH, WMESH
  integer(4) :: ICMAX

  ! Mesh parameters
  integer(4) :: NRMAX, NEMAX, NRA, NRC

  ! Equilibrium parameters
  integer(4) :: ieqread

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
  integer(4) :: IDIAG, NTSTEP

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
  integer(4) :: MDOSQZ, MDLETA, MDLETB, MDLNEO, MDBSETA

  ! Initial condition
  integer(4) :: MDITSN, MDITST, MDINTN, MDINTT, MDINTC

  !**********************************!
  !   INTERNAL CONTOROL PARAMETERS   !
  !**********************************!

  ! Time control parameters
  integer(4) :: NT, NTCUM, NTMAX
  real(8) :: T_TX = 0.D0, TMAX

  ! Configuration parameters
  integer(4) :: NQMAX, IERR, ICONT, IRPIN
  real(4) :: AVE_IC
  real(8) :: rIP
  real(8) :: UHth, UHph
  real(8) :: Rax, Zax, perimlcfs
  real(8), dimension(:), allocatable :: R, vv, vvn, rho, rhosq
  real(8), dimension(:), allocatable :: hv

  ! Convergence accelerator
  real(8) :: rMUb1, rMUb2

  ! Convergence parameter
  integer(4) :: MODECV
  real(8) :: oldmix

  ! SUPG parameter
  integer(4) :: iSUPG2, iSUPG3, iSUPG6
  real(8) :: SUPGstab

  ! Variables
  real(8), dimension(:), allocatable :: &
       & ErV,    BEpol,   Etor,    BthV,  BphV, &
       & PNbV,   UbrVV,   BUbparV, RUbphV, PTbV, UbphVR, UbphV, &
       & PN01V,  PN02V,   PN03V, &
       & PsiV,   PsitV,   PhiV,  &
       & PT01V,  PT02V,   PT03V,   PNbrpV, &
       & PsidotV,PsitdotV,dPhiV,   bbt,   bthco
!!rp_conv       &, PNbrpLV

  real(8), dimension(:), allocatable :: ErV_FIX

  ! Species variables
  type species_var
     real(8) :: n, T, p, UrV, BUpar, Bqpar, Uthhat, UphR, RUph, Ur, Uth, Uph
  end type species_var
  type(species_var), dimension(:,:), allocatable :: Var

  real(8), dimension(:,:), allocatable :: PNsV_FIX, PTsV_FIX

  real(8), dimension(:,:,:), allocatable :: BVsdiag, UsparNCL, UsthNCL

  ! Coefficients
  real(8), dimension(:), allocatable :: &
       & rNuION, rNu0e, rNu0i, rNu0b, rNuL, rNuiCX, rNubCX, rNuiCXT, &
       & rNuei, rNuii, rNuTei, rNuD, &
       & rNubrp1, rNubrp2, &
       & rNuLTe, rNuLTi, &
       & rNuAse, rNuAsi, &
!       & rNueHL, rNuiHL, &
       & rNueHLthth,rNueHLthph, rNueHLphth, rNueHLphph, &
       & rNuiHLthth,rNuiHLthph, rNuiHLphth, rNuiHLphph, &
       & WPM, FVpch, FQLcoef, FQLcoef1, FQLcoef2, &
       & rMue, rMui, rNuB, rNuLB, ft, &
       & Chie, Chii, De, Di, VWpch, D01, D02, D03, &
       & ChiNCpe, ChiNCte, ChiNCpi, ChiNCti, &
       & DMAG, DMAGe, DMAGi, &
       & Ubrp, RUbrp, Dbrp, DltRP, DltRP_mid, rNubL, rip_rat, rNuOL, Vbpara, &
       & SiVizA, SiVcxA, wexb, &
       & UgV, PNbVinv
  real(8), dimension(:,:), allocatable :: gamITG
  real(8), dimension(:,:), allocatable :: Vhps, Vmps, PiRess

  ! Coefficients related to neoclassical transport
  real(8), dimension(:,:), allocatable :: xmuf, BnablaPi, gflux
  real(8), dimension(:,:,:), allocatable :: lff
  real(8), dimension(:,:,:,:), allocatable :: xmu, laf, lfb
  real(8), dimension(:,:,:,:,:), allocatable :: lab

  ! Equilibrium metrics 
  real(8), dimension(:), allocatable :: &
       & aat, rrt, ckt, suft, sst, vro, vlt, rhov, art, epst, d_rrr, elip, trig
  real(8), dimension(:), allocatable :: fipol, Bpsq, qhatsq, Fqhatsq, BEpara, bri
  real(8), dimension(:), allocatable :: sdt, hdt
  real(8), dimension(:), allocatable :: bit, bbrt

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
       &                                PNBe, PNBi, MNB, PRFe, PRFi, &
       &                                POH, POHe, POHi, PEQe, PEQi, &
       &                                SiLC, SiLCth, SiLCph, PALFe, PALFi, &
       &                                BSmb, Tqt, Tqp
  real(8), dimension(:), allocatable :: PIE, PCX, SIE, SCX, PBr
  real(8) :: Eb, Vb, RatCX

  ! Safety factor, currents, resistivity
  real(8), dimension(:), allocatable :: Q, AJ, BJPARA, AJOH, BJOH, AJV, AJRF, AJNB, BJNB, &
       &                                BJBS, AJBS, ETA, ETAS 
  real(8), dimension(:,:), allocatable :: BJBSvar, ETAvar

  ! Derivatives
  real(8), dimension(:,:), allocatable :: dPsdpsi, dTsdpsi

  ! For display
  real(8), dimension(:), allocatable :: ANS0, TS0, ANSAV, TSAV, WST
  real(8), dimension(:), allocatable :: ANF0, TF0, ANFAV, TFAV, WFT
  real(8), dimension(:), allocatable :: PBCLT, PFCLT, PLT, SPET, SLT
  real(8), dimension(:), allocatable :: Deff, thrp
  real(8) :: WBULKT, WTAILT, WPT
  real(8) :: AJT, AJOHT, AJNBT, AJRFT, AJBST
  real(8) :: PINT, POHT, PNBT, PRFT, PRFTe, PRFTi, PNFT
  real(8) :: PBINT, PFINT, POUT, PCXT, PIET, PRLT, SINT, SIET
  real(8) :: SNBT, SNFT, SOUT, TNBcol, TTqt
  real(8) :: VLOOP, ALI, RQ1, RPE, ZEFF0, QF
  real(8) :: WPDOT, TAUE1, TAUE2, TAUEP, TAUEH, TAUP, TAUPA
  real(8) :: BETAP0, BETAPA, BETA0, BETAA, BETAQ0, BETAN
  real(8) :: TPRE, WPPRE
  real(8) :: VOLAVN, Gamma_a

  ! Internal variables for transport matrix
  real(8),    dimension(:,:,:), allocatable :: ALC, BLC, CLC, PLC
  integer(4), dimension(:,:),   allocatable :: NLC
  integer(4), dimension(:,:,:), allocatable :: NLCR
  integer(4), dimension(:),     allocatable :: NLCMAX
  real(8),    dimension(:,:),   allocatable :: X, XOLD

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
    integer(4) :: iflag, N, NS, NF
    integer(4), dimension(1:20) :: ierl

    ierl(:) = 0
    if(nrmax <= 1) then
      write(6,*) "XXX ALLOCATE_TXCOMM : ILLEGAL PARAMETER    NRMAX=", nrmax
      ier = 1
      return
    endif

    N     = NRMAX
    NEMAX = NRMAX
    NS    = NSM
    NF    = NFM

    ! allocation check
    if(allocated(R)) then
       if(present(icont_in) .and. icont_in /= 0) then
          call deallocate_txcomm
       else
          return
       end if
    end if

    do
       allocate(R(0:N),      vv(0:N),    vvn(0:N),   rho(0:N),   rhosq(0:N),  stat = ierl(1))
       allocate(hv(0:NEMAX),                                                  stat = ierl(2))
       ier = sum(ierl) ; iflag = 1
       if (ier /= 0) exit

       allocate(ErV(0:N),    BEpol(0:N), Etor(0:N),  BthV(0:N),  BphV(0:N),    stat = ierl(1))
       allocate(PNbV(0:N),   UbrVV(0:N), BUbparV(0:N), RUbphV(0:N), PTbV(0:N), UbphVR(0:N), UbphV(0:N), stat = ierl(2))
       allocate(PN01V(0:N),  PN02V(0:N), PN03V(0:N),                           stat = ierl(3))
       allocate(PsiV(0:N),   PsitV(0:N), PhiV(0:N),  bbt(0:N), bthco(0:N),     stat = ierl(4))
       allocate(PT01V(0:N),  PT02V(0:N), PT03V(0:N), PNbrpV(0:N),              stat = ierl(5))
       allocate(PsidotV(0:N),PsitdotV(0:N),dPhiV(0:N),                         stat = ierl(6))
       allocate(Var(0:N,NS),                                                   stat = ierl(7))
       ier = sum(ierl) ; iflag = 2
       if (ier /= 0) exit

       allocate(ErV_FIX(0:N),                                                  stat = ierl(1))
       allocate(PNsV_FIX(0:N,NS), PTsV_FIX(0:N,NS),                            stat = ierl(2))
       allocate(BVsdiag(0:N,NS,2), UsparNCL(0:N,NS,2), UsthNCL(0:N,NS,2),      stat = ierl(3))
       ier = sum(ierl) ; iflag = 3
       if (ier /= 0) exit

       allocate(rNuION(0:N), rNu0e(0:N),  rNu0i(0:N), rNu0b(0:N), rNuL(0:N),  stat = ierl(1))
       allocate(rNuiCX(0:N), rNubCX(0:N), rNuei(0:N), rNuiCXT(0:N),           stat = ierl(2))
       allocate(rNuii(0:N),  rNuTei(0:N), rNuD(0:N),                          stat = ierl(3))
       allocate(rNubrp1(0:N),rNubrp2(0:N),                                    stat = ierl(4))
       allocate(rNuLTe(0:N),rNuLTi(0:N),                                      stat = ierl(5))
       allocate(rNuAse(0:N),rNuAsi(0:N),                                      stat = ierl(6))
!       allocate(rNueHL(0:N), rNuiHL(0:N), rNueHLthth(0:N), rNuiHLthth(0:N), rNueHLthph(0:N), &
!            &   rNuiHLthph(0:N), rNueHLphth(0:N), rNuiHLphth(0:N),rNueHLphph(0:N), &
!            &   rNuiHLphph(0:N),                                              stat = ierl(7))
!       allocate(HPN(1:2,1:NHFMmx),EpsHM(1:2,1:NHFMmx),                        stat = ierl(7))
       allocate(rNueHLM(1:NHFMmx,0:N), rNuiHLM(1:NHFMmx,0:N), &
            &   rNueHLthth(0:N), rNuiHLthth(0:N), rNueHLthph(0:N), rNuiHLthph(0:N), &
            &   rNueHLphth(0:N), rNuiHLphth(0:N), rNueHLphph(0:N), rNuiHLphph(0:N), &
            &   rNueHLththM(1:NHFMmx,0:N), rNueHLthphM(1:NHFMmx,0:N), &
            &   rNueHLphthM(1:NHFMmx,0:N), rNueHLphphM(1:NHFMmx,0:N), &
            &   rNuiHLththM(1:NHFMmx,0:N), rNuiHLthphM(1:NHFMmx,0:N), &
            &   rNuiHLphthM(1:NHFMmx,0:N), rNuiHLphphM(1:NHFMmx,0:N),         stat = ierl(7))
       ! miki_m 10-08-06
       allocate(WPM(0:N),    FVpch(0:N),  FQLcoef(0:N), FQLcoef1(0:N), FQLcoef2(0:N),stat = ierl(8))
       allocate(rMue(0:N),   rMui(0:N),                                       stat = ierl(9))
       allocate(rNuB(0:N),   rNuLB(0:N),  ft(0:N),    Chie(0:N),  Chii(0:N),  stat = ierl(10))
       allocate(De(0:N),     Di(0:N),     VWpch(0:N), D01(0:N),    D02(0:N),  &
            &   D03(0:N),                                                     stat = ierl(11))
       allocate(ChiNCpe(0:N),ChiNCte(0:N),ChiNCpi(0:N),ChiNCti(0:N),          stat = ierl(12))
       allocate(Ubrp(0:N),   RUbrp(0:N),  Dbrp(0:N),  DltRP(0:N), DltRP_mid(0:N), &
            &   rNubL(0:N),                                                   stat = ierl(13))
       allocate(rip_rat(0:N),rNuOL(0:N),  Vbpara(0:N),                        stat = ierl(14))
       allocate(gamITG(0:N,1:3),                                              stat = ierl(15))
       allocate(DMAG(0:N),   DMAGe(0:N),  DMAGi(0:N),                         stat = ierl(16))
       allocate(SiVizA(0:N), SiVcxA(0:N), wexb(0:N),                          stat = ierl(17))
       allocate(UgV(0:N),    PNbVinv(0:N),                                    stat = ierl(18))
       allocate(Vhps(0:N,NS),Vmps(0:N,NS), PiRess(0:N,NS),                    stat = ierl(19))
       ier = sum(ierl) ; iflag = 4
       if (ier /= 0) exit

       allocate(xmu(0:N,NS,3,3), lab(0:N,NS,NS,3,3),                          stat = ierl(1))
       allocate(laf(0:N,NS,2,2), lfb(0:N,NS,2,2), lff(0:N,2,2), xmuf(0:N,3),  stat = ierl(2))
       allocate(BnablaPi(0:N,NS),gflux(0:N,NS),                               stat = ierl(3))
       ier = sum(ierl) ; iflag = 5
       if (ier /= 0) exit

       allocate(aat(0:N), rrt(0:N), ckt(0:N), suft(0:N), sst(0:N), vro(0:N), vlt(0:N), rhov(0:N), &
            &   art(0:N), epst(0:N), d_rrr(0:N), elip(0:N), trig(0:N), stat = ierl(1))
       allocate(fipol(0:N), Bpsq(0:N), qhatsq(0:N), Fqhatsq(0:N), BEpara(0:N), bri(0:N),stat = ierl(2))
       allocate(sdt(0:N), hdt(0:N),                                           stat = ierl(3))
       allocate(bit(0:N), bbrt(0:N),                                          stat = ierl(4))
       ier = sum(ierl) ; iflag = 6
       if (ier /= 0) exit

       allocate(rG1h2(0:N),  FCDBM(0:N),  S(0:N),     Alpha(0:N), rKappa(0:N),stat = ierl(1))
       allocate(rG1h2IM(0:N),  FCDIM(0:N),  RAQPR(0:N),                       stat = ierl(2)) !09/06/17 miki_m
       allocate(pres0(0:N),  ErV0(0:N),                                       stat = ierl(3))
       ier = sum(ierl) ; iflag = 7
       if (ier /= 0) exit

       allocate(PNB(0:N),   PNBTG(0:N),PNBPD(0:N),PNBcol_e(0:N),PNBcol_i(0:N),stat = ierl(1))
       allocate(SNB(0:N),   SNBe(0:N),   SNBi(0:N),   SNBb(0:N),              stat = ierl(2))
       allocate(SNBPDi(0:N),SNBTGi(0:N),                                      stat = ierl(3))
       allocate(PNBe(0:N),  PNBi(0:N),   MNB(0:N),    PRFe(0:N),  PRFi(0:N),  stat = ierl(4))
       allocate(POH(0:N),   POHe(0:N),   POHi(0:N),   PEQe(0:N),  PEQi(0:N),  stat = ierl(5))
       allocate(SiLC(0:N),  SiLCth(0:N), SiLCph(0:N), PALFe(0:N), PALFi(0:N), stat = ierl(6))
       allocate(PIE(0:N),   PCX(0:N),    SIE(0:N),    SCX(0:N),   PBr(0:N),   stat = ierl(7))
       allocate(BSmb(0:N),  Tqt(0:N),    Tqp(0:N),                            stat = ierl(8))
       ier = sum(ierl) ; iflag = 8
       if (ier /= 0) exit

       allocate(Q(0:N), AJ(0:N), BJPARA(0:N), AJOH(0:N), BJOH(0:N), AJV(0:N), stat = ierl(1))
       allocate(AJRF(0:N), AJNB(0:N), BJNB(0:N), BJBS(0:N), AJBS(0:N),        stat = ierl(2))
       allocate(ETA(0:N), ETAS(0:N),                                          stat = ierl(3))
       allocate(BJBSvar(0:N,0:3), ETAvar(0:N,0:4),                            stat = ierl(4))
       ier = sum(ierl) ; iflag = 9
       if (ier /= 0) exit

       allocate(dPsdpsi(0:NRMAX,NSM), dTsdpsi(0:NRMAX,NSM),                   stat = ierl(1))
       ier = sum(ierl) ; iflag = 10
       if (ier /= 0) exit

       allocate(ANS0(1:NS), TS0(1:NS),   ANSAV(1:NS), TSAV(1:NS), WST(1:NS),  stat = ierl(1))
       allocate(ANF0(1:NF), TF0(1:NF),   ANFAV(1:NF), TFAV(1:NF), WFT(1:NF),  stat = ierl(2))
       allocate(PBCLT(1:NS),PFCLT(1:NS), PLT(1:NS),   SPET(1:NS), SLT(1:NS),  stat = ierl(3))
       allocate(Deff(0:N),  thrp(1:2*N),                                      stat = ierl(4))
       ier = sum(ierl) ; iflag = 11
       if (ier /= 0) exit

       allocate(ALC(0:N,0:NCM,1:NQMAX),BLC(0:N,0:NCM,1:NQMAX),CLC(0:N,0:NCM,1:NQMAX),stat = ierl(1))
       allocate(PLC(0:N,1:NCM,1:NQMAX),                                              stat = ierl(2))
       allocate(NLC(0:NCM,1:NQMAX),NLCR(0:NCM,0:NQMAX,0:1),NLCMAX(1:NQMAX),          stat = ierl(3))
       allocate(X(0:N,1:NQMAX), XOLD(0:N,1:NQMAX),                                   stat = ierl(4))
       ier = sum(ierl) ; iflag = 12
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

    deallocate(R,      vv,      vvn,    rho,   rhosq)
    deallocate(hv)

    deallocate(ErV,    BEpol,   Etor,   BthV,  BphV)
    deallocate(PNbV,   UbrVV,   BUbparV,RUbphV, PTbV, UbphVR, UbphV)
    deallocate(PN01V,  PN02V,   PN03V)
    deallocate(PsiV,   PsitV,   PhiV)
    deallocate(PT01V,  PT02V,   PT03V,  PNbrpV)
    deallocate(PsidotV,PsitdotV,dPhiV,  bbt,   bthco)
    deallocate(Var)

    deallocate(ErV_FIX)
    deallocate(PNsV_FIX, PTsV_FIX)
    deallocate(BVsdiag, UsparNCL, UsthNCL)

    deallocate(rNuION, rNu0e,  rNu0i, rNu0b, rNuL)
    deallocate(rNuiCX, rNubCX, rNuei, rNuiCXT)
    deallocate(rNuii,  rNuTei, rNuD)
    deallocate(rNubrp1,rNubrp2)
    deallocate(rNuLTe, rNuLTi)
    deallocate(rNuAse, rNuAsi)
!    deallocate(rNueHL, rNuiHL, rNueHLthth, rNuiHLthth, rNueHLthph, &
!         &     rNuiHLthph, rNueHLphth, rNuiHLphth, rNueHLphph, &
!         &     rNuiHLphph)
!    deallocate(HPN, EpsHM) ! miki_m 10-08-06
    deallocate(rNueHLM, rNuiHLM, &
         &     rNueHLthth, rNuiHLthth, rNueHLthph, rNuiHLthph, &
         &     rNueHLphth, rNuiHLphth, rNueHLphph, rNuiHLphph, &
         &     rNueHLththM, rNueHLthphM, rNueHLphthM, rNueHLphphM, &
         &     rNuiHLththM, rNuiHLthphM, rNuiHLphthM, rNuiHLphphM) ! miki_m 10-08-06
    deallocate(WPM,   FVpch, FQLcoef, FQLcoef1, FQLcoef2)
    deallocate(rMue,  rMui)
    deallocate(rNuB,   rNuLB,  ft,    Chie,  Chii)
    deallocate(De,     Di,     VWpch, D01,   D02,  D03)
    deallocate(ChiNCpe,ChiNCte,ChiNCpi,ChiNCti)
    deallocate(Ubrp,   RUbrp,  Dbrp,  DltRP, DltRP_mid, rNubL)
    deallocate(rip_rat,rNuOL,  Vbpara)
    deallocate(gamITG)
    deallocate(DMAG,   DMAGe,  DMAGi)  !***AF (2008-06-08)
    deallocate(SiVizA, SiVcxA, wexb)
    deallocate(UgV,    PNbVinv)
    deallocate(Vhps,   Vmps,   PiRess)

    deallocate(xmu, xmuf, lab, laf, lfb, lff, BnablaPi, gflux)

    deallocate(aat, rrt, ckt, suft, sst, vro, vlt, rhov, art, epst, d_rrr, elip, trig)
    deallocate(fipol, Bpsq, qhatsq, Fqhatsq, BEpara, bri)
    deallocate(sdt, hdt)
    deallocate(bit, bbrt)

    deallocate(rG1h2,  FCDBM,  S,     Alpha, rKappa)
    deallocate(pres0,  ErV0)

    deallocate(rG1h2IM, FCDIM, RAQPR)  !***miki_m 09/06/17~

    deallocate(PNB,    PNBTG, PNBPD, PNBcol_e, PNBcol_i)
    deallocate(SNB,    SNBe,  SNBi,  SNBb,   SNBPDi, SNBTGi)
    deallocate(PNBe,   PNBi,  MNB,   PRFe,   PRFi)
    deallocate(POH,    POHe,  POHi,  PEQe,   PEQi)
    deallocate(SiLC,   SiLCth,SiLCph,PALFe,  PALFi)
    deallocate(PIE,    PCX,   SIE,   SCX,    PBr)
    deallocate(BSmb,   Tqt,   Tqp)

    deallocate(Q, AJ, BJPARA, AJOH, BJOH, AJV)
    deallocate(AJRF, AJNB, BJNB, BJBS, AJBS)
    deallocate(ETA, ETAS)
    deallocate(BJBSvar, ETAvar)

    deallocate(dPsdpsi, dTsdpsi)

    deallocate(ANS0, TS0, ANSAV, TSAV, WST)
    deallocate(ANF0, TF0, ANFAV, TFAV, WFT)
    deallocate(PBCLT, PFCLT, PLT, SPET, SLT)
    deallocate(Deff, thrp)

    deallocate(ALC, BLC, CLC)
    deallocate(PLC)
    deallocate(NLC,NLCR,NLCMAX)
    deallocate(X, XOLD)

    if(allocated(infiles)) deallocate(infiles)

  end subroutine deallocate_txcomm

end module tx_commons
