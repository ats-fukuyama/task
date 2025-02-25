module tx_commons
  implicit none
  public

  integer(4), parameter :: NQM=22, NCM=32, NGRM=20, &
       &                   NGTM=5000, NGVM=5000, NGYRM=143, NGYTM=51, &
       &                   NGYVM=56, NGPRM=21, NGPTM=8, NGPVM=15, &
       &                   NMNQM=446, M_POL_M=64
!09/07/13 miki_m original : NGYRM=136
  integer(4), parameter :: NSM=2, NFM=2
  integer(4), parameter :: LQm1=1,  LQm2=2,  LQm3=3,  LQm4=4,  LQm5=5,&
       &                   LQe1=6,  LQe2=7,  LQe3=8,  LQe4=9,  LQe5=10,&
       &                   LQi1=11, LQi2=12, LQi3=13, LQi4=14, LQi5=15,&
       &                   LQb1=16,          LQb3=17, LQb4=18,&
       &                   LQn1=19, LQn2=20, LQn3=21, &
       &                   LQr1=22
  integer(4), parameter :: nmax_file = 100

  !**********************************************!
  ! Set physical constants, based on CODATA 2010 !
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

  !   Conversion factor from keV to joule
  real(8), parameter :: rKeV = 1.D3 * AEE
  
  !   Hydrogen ionization energy
  real(8), parameter :: EION = 13.6D0

  !   Electron-proton mass ratio
  real(8), parameter :: AEP = AME / AMP

  !**********************!
  !   INPUT PARAMETERS   !
  !**********************!

  ! Configuration parameters
  real(8) :: RA, RB, RACCUM, RR, BB, rIPs, rIPe

  ! Species
  real(8) :: PA, PZ, Zeff

  ! Profiles
  real(8) :: PN0, PNa, PTe0, PTea, PTi0, PTia, PROFJ, PROFN1, PROFN2, PROFT1, PROFT2

  ! Diffusivities and viscosities
  real(8) :: De0, Di0, rMue0, rMui0, Chie0, Chii0, ChiNC, VWpch0, WPM0

  ! Amplitude parameters for transport
  real(8) :: FSCBKP, FSCBEL, FSCBSH, FSBOHM, FSPCLD, FSPCLM, FSPCLC, FSVAHL
  real(8) :: PROFD, PROFD1, PROFD2, PROFM, PROFM1, PROFC, PROFC1
  real(8) :: FSCX, FSLC, FSNC, FSLP, FSLTE, FSLTI, FSION
  real(8) :: FSD01, FSD02, FSD03, FSNCPL, rG1, FSRP, FSNF, FSLIM
  real(8), dimension(1:3) :: FSDFIX, FSANOM, RhoETB
  integer(4) :: MDLC, MDANOM

  ! Scale lengths in SOL
  real(8) :: rLn, rLT

  ! Heat sources
  real(8) :: Eb, RNBP, RNBP0, RNBT1, RNBT2, RNBT10, RNBT20, PNBH, PNBHP, PNBHT1, PNBHT2, PNBHex, &
       &     PNBCD, PNBMPD, &
       &     rNRFe, RRFew, RRFe0, PRFHe, rNRFi, RRFiw, RRFi0, PRFHi, Tqt0, Tqp0
  integer(4) :: MDLNBD, MDLMOM

  ! Neutral parameters
  real(8) :: PN0s, V0, rGamm0, rGASPF, PNeDIV, PTeDIV, PTiDIV

  ! Ripple parameters
  real(8) :: DltRPn, kappa
  integer(4) :: NTCOIL, m_pol, n_tor

  ! Numerical parameters
  real(8) :: DT, EPS, ADV, tiny_cap, CMESH0, WMESH0, CMESH, WMESH
  integer(4) :: ICMAX

  ! Mesh parameters
  integer(4) :: NRMAX, NEMAX, NTMAX, NTSTEP, NGRSTP, NGTSTP, NGVSTP, NRA, NRC

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

  ! Diagnostic message
  integer(4) :: IDIAG

  ! Graphic patameters --- 4 byte real
  real(4) :: GRMIN,GRMAX

  ! LAPACK
  integer(4) :: MDLPCK

  ! Implicit Euler or Gear's backward differential formula
  integer(4) :: IGBDF

  ! Nonlinear iteration method
  integer(4) :: MDSOLV

  !  Transport model
  integer(4) :: MDOSQZ, MDLETA, MDFIXT, MDVAHL, MDLETB

  !  Initial condition
  integer(4) :: MDITSN, MDITST, MDINTN, MDINTT, MDINTC, MDINIT

  !**********************************!
  !   INTERNAL CONTOROL PARAMETERS   !
  !**********************************!

  ! Configuration parameters
  integer(4) :: NT, NQMAX, IERR, ICONT, IRPIN
  real(4) :: AVE_IC
  real(8) :: T_TX = 0.D0, TMAX
  real(8) :: AMI, AMB, Vb, sqeps0
  real(8) :: rIP, Bthb
  real(8) :: UHth, UHph
  real(8), dimension(:), allocatable :: R, PSI, RHO
  real(8), dimension(:), allocatable :: H, HPSI

  ! Convergence accelerator
  real(8) :: rMUb1, rMUb2

  ! Convergenve parameter
  integer(4) :: MODECV
  
  ! Variables
  real(8), dimension(:), allocatable :: &
       & ErV,    EthV,   EphV,  BthV,  BphV, &
       & PNeV,   UerV,   UethV, UephV, PTeV, &
       & PNiV,   UirV,   UithV, UiphV, PTiV, &
       & PNbV,   UbthV,  UbphV, PN01V, PN02V, PN03V, &
       & AphV,   PhiV,   RAthV, PeV,   PiV,  &
       & RUethV, RUithV, PT01V, PT02V, PT03V, PNbrpV
!!rp_conv       &, PNbrpLV

  real(8), dimension(:), allocatable :: PNeV_FIX, PTeV_FIX, PNiV_FIX, PTiV_FIX, ErV_FIX

  ! Coefficients
  real(8), dimension(:), allocatable :: &
       & rNuION, rNu0e, rNu0i, rNu0b, rNuL, rNuiCX, rNubCX, rNuiCXT, &
       & rNuee, rNuei, rNuie, rNuii, rNuTei, rNube, rNubi, rNuD, &
       & rNubrp1, rNubrp2, rNuei1, rNuei2, rNuei3, rNuei2Bth, &
       & rNube1, rNube2, rNube3, rNube2Bth, rNuLTe, rNuLTi, &
       & rNueNC, rNuiNC, rNue2NC, rNui2NC, rNuAse, rNuAsi, &
!       & rNueHL, rNuiHL, &
       & rNueHLthth,rNueHLthph, rNueHLphth, rNueHLphph, &
       & rNuiHLthth,rNuiHLthph, rNuiHLphth, rNuiHLphph, &
       & FWe, FWi, FWthe, FWthi, FWthphe, FWthphi, &
       & WPM, FVpch, &
       & rMue, rMui, rNuB, rNuLB, ft, &
       & Chie, Chii, De, Di, VWpch, D01, D02, D03, &
       & ChiNCpe, ChiNCte, ChiNCpi, ChiNCti, &
       & DMAG, DMAGe, DMAGi, &
       & Ubrp, RUbrp, Dbrp, DltRP, DltRP_mid, rNubL, rip_rat, rNuOL, &
       & rNuNTV, UastNC, Vbpara, &
       & SiVizA, SiVcxA, wexb, Ys, FQeth, FQith
  real(8), dimension(:,:), allocatable :: deltam, gamITG

  ! Multiple helical Fourier modes    miki_m 10-08-06~
  integer(4) :: UHphSwitch           !  To support (m=0, n>0) Fourier component
  real(8), dimension(:,:), allocatable :: rNueHLM, rNuiHLM, &
       & rNueHLththM, rNueHLthphM, rNueHLphthM, rNueHLphphM,&
       & rNuiHLththM, rNuiHLthphM, rNuiHLphthM, rNuiHLphphM 

  ! Read Wnm table
  real(8), dimension(:),   allocatable :: Fmnq, Wnm
  real(8), dimension(:,:), allocatable :: Umnq
 
  ! CDBM
  real(8), dimension(:), allocatable :: rG1h2, FCDBM, S, Alpha, rKappa

  ! For numerical stability
  real(8), dimension(:), allocatable :: pres0, ErV0

  ! CDIM
  ! 09/06/17~ miki_m
  real(8), dimension(:), allocatable :: rG1h2IM, FCDIM, RAQPR

  ! Sources and sinks
  real(8), dimension(:), allocatable :: PNB, PNBTG, PNBPD, PNBcol_e, PNBcol_i,  &
       &                                SNB, SNBe, SNBi, SNBb, SNBPDi, SNBTGi, Tqt, Tqp, &
       &                                PNBe, PNBi, MNB, PRFe, PRFi, &
       &                                POH, POHe, POHi, PEQe, PEQi, &
       &                                SiLC, SiLCth, SiLCph, PALFe, PALFi
  real(8), dimension(:), allocatable :: PIE, PCX, SIE, SCX, PBr
  real(8) :: RatCX

  ! Safety factor, currents, resistivity
  real(8), dimension(:), allocatable :: Q, AJ, AJOH, AJV, AJRF, AJNB, AJPARA, &
       &                                AJBS, AJBS1, AJBS2, AJBS3, AJBS4, &
       &                                ETA, ETAS, ETA1, ETA2, ETA3, ETA4

  ! For display
  real(8), dimension(:), allocatable :: ANS0, TS0, ANSAV, TSAV, WST
  real(8), dimension(:), allocatable :: ANF0, TF0, ANFAV, TFAV, WFT
  real(8), dimension(:), allocatable :: PBCLT, PFCLT, PLT, SPET, SLT
  real(8), dimension(:), allocatable :: Deff, thrp
  real(8), dimension(:,:), allocatable :: Ueth_NC, Uith_NC
  real(8) :: WBULKT, WTAILT, WPT
  real(8) :: AJT, AJOHT, AJNBT, AJRFT, AJBST
  real(8) :: PINT, POHT, PNBT, PRFT, PRFTe, PRFTi, PNFT
  real(8) :: PBINT, PFINT, POUT, PCXT, PIET, PRLT, SINT, SIET
  real(8) :: SNBT, SNFT, SOUT
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
  
  ! Graphics
  integer(4) :: NGR, NGT, NGVV
  integer(4) :: MODEG, MODEAV, MODEGL
  real(4), dimension(:),     allocatable :: GX
  real(4), dimension(:,:,:), allocatable :: GY, GQY, GYT
  real(4), dimension(0:NGRM) :: GT
  real(4), dimension(0:NGTM) :: GTX
  real(4), dimension(0:NGTM,1:NGYTM) :: GTY
  real(4), dimension(1:NGYRM) :: gDIV
  real(4), dimension(0:NGVM) :: GVX
  real(4), dimension(0:NGVM,1:NGYVM) :: GVY

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
       allocate(R(0:N),      PSI(0:N),    RHO(0:N),                           stat = ierl(1))
       allocate(H(0:NEMAX),  HPSI(0:NEMAX),                                   stat = ierl(2))
       ier = sum(ierl) ; iflag = 1
       if (ier /= 0) exit

       allocate(ErV(0:N),    EthV(0:N),   EphV(0:N),  BthV(0:N),  BphV(0:N),  stat = ierl(1))
       allocate(PNeV(0:N),   UerV(0:N),   UethV(0:N), UephV(0:N), PTeV(0:N),  stat = ierl(2))
       allocate(PNiV(0:N),   UirV(0:N),   UithV(0:N), UiphV(0:N), PTiV(0:N),  stat = ierl(3))
       allocate(PNbV(0:N),   UbthV(0:N),  UbphV(0:N), PN01V(0:N), PN02V(0:N), &
            &   PN03V(0:N),                                                   stat = ierl(4))
       allocate(AphV(0:N),   PhiV(0:N),   RAthV(0:N), PeV(0:N),   PiV(0:N),   stat = ierl(5))
       allocate(RUethV(0:N), RUithV(0:N), PT01V(0:N), PT02V(0:N), PT03V(0:N), &
            &   PNbrpV(0:N),                                                  stat = ierl(6))
       ier = sum(ierl) ; iflag = 2
       if (ier /= 0) exit

       allocate(PNeV_FIX(0:N), PTeV_FIX(0:N), PNiV_FIX(0:N), PTiV_FIX(0:N),   &
            &   ErV_FIX(0:N),                                                 stat = ierl(1))
       ier = sum(ierl) ; iflag = 3
       if (ier /= 0) exit

       allocate(rNuION(0:N), rNu0e(0:N),  rNu0i(0:N), rNu0b(0:N), rNuL(0:N),  stat = ierl(1))
       allocate(rNuiCX(0:N), rNubCX(0:N), rNuee(0:N), rNuei(0:N), rNuie(0:N), &
            &   rNuiCXT(0:N),                                                 stat = ierl(2))
       allocate(rNuii(0:N),  rNuTei(0:N), rNube(0:N), rNubi(0:N), rNuD(0:N),  stat = ierl(3))
       allocate(rNubrp1(0:N),rNubrp2(0:N),rNuei1(0:N),rNuei2(0:N),rNuei3(0:N), &
            &   rNuei2Bth(0:N),                                               stat = ierl(4))
       allocate(rNube1(0:N), rNube2(0:N), rNube3(0:N),rNuLTe(0:N),rNuLTi(0:N), &
            &   rNube2Bth(0:N),                                               stat = ierl(5))
       allocate(rNueNC(0:N), rNuiNC(0:N), rNuAse(0:N),rNuAsi(0:N),            stat = ierl(6))
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
       allocate(FWe(0:N),    FWi(0:N),    FWthe(0:N), FWthi(0:N), FWthphe(0:N), &
            &   FWthphi(0:N),                                                 stat = ierl(8))
       allocate(WPM(0:N),    FVpch(0:N),                                      stat = ierl(9))
       allocate(rMue(0:N),   rMui(0:N),                                       stat = ierl(10))
       allocate(rNuB(0:N),   rNuLB(0:N),  ft(0:N),    Chie(0:N),  Chii(0:N),  stat = ierl(11))
       allocate(De(0:N),     Di(0:N),     VWpch(0:N), D01(0:N),    D02(0:N),  &
            &   D03(0:N),                                                     stat = ierl(12))
       allocate(ChiNCpe(0:N),ChiNCte(0:N),ChiNCpi(0:N),ChiNCti(0:N),          stat = ierl(13))
       allocate(Ubrp(0:N),   RUbrp(0:N),  Dbrp(0:N),  DltRP(0:N), DltRP_mid(0:N), &
            &   rNubL(0:N),                                                   stat = ierl(14))
       allocate(rip_rat(0:N),rNuOL(0:N),  rNuNTV(0:N),UastNC(0:N),Vbpara(0:N),stat = ierl(15))
       allocate(Fmnq(1:NMNQM), Wnm(1:NMNQM), Umnq(1:4,1:NMNQM),               stat = ierl(16))
       allocate(deltam(0:N,0:M_POL_M),gamITG(0:N,1:3),                        stat = ierl(17))
       allocate(DMAG(0:N),   DMAGe(0:N),  DMAGi(0:N),                         stat = ierl(18))
       allocate(SiVizA(0:N), SiVcxA(0:N), wexb(0:N),  Ys(0:N),                stat = ierl(19))
       allocate(rNue2NC(0:N),rNui2NC(0:N),FQeth(0:N), FQith(0:N),             stat = ierl(20))
       ier = sum(ierl) ; iflag = 4
       if (ier /= 0) exit

       allocate(rG1h2(0:N),  FCDBM(0:N),  S(0:N),     Alpha(0:N), rKappa(0:N),stat = ierl(1))
       allocate(rG1h2IM(0:N),  FCDIM(0:N),  RAQPR(0:N),                       stat = ierl(2)) !09/06/17 miki_m
       allocate(pres0(0:N),  ErV0(0:N),                                       stat = ierl(3))
       ier = sum(ierl) ; iflag = 5
       if (ier /= 0) exit

       allocate(PNB(0:N),   PNBTG(0:N),PNBPD(0:N),PNBcol_e(0:N),PNBcol_i(0:N),stat = ierl(1))
       allocate(SNB(0:N),   SNBe(0:N),   SNBi(0:N),   SNBb(0:N),              stat = ierl(2))
       allocate(SNBPDi(0:N),SNBTGi(0:N), Tqt(0:N),    Tqp(0:N),               stat = ierl(3))
       allocate(PNBe(0:N),  PNBi(0:N),   MNB(0:N),    PRFe(0:N),  PRFi(0:N),  stat = ierl(4))
       allocate(POH(0:N),   POHe(0:N),   POHi(0:N),   PEQe(0:N),  PEQi(0:N),  stat = ierl(5))
       allocate(SiLC(0:N),  SiLCth(0:N), SiLCph(0:N), PALFe(0:N), PALFi(0:N), stat = ierl(6))
       allocate(PIE(0:N),   PCX(0:N),    SIE(0:N),    SCX(0:N),   PBr(0:N),   stat = ierl(7))
       ier = sum(ierl) ; iflag = 6
       if (ier /= 0) exit

       allocate(Q(0:N), AJ(0:N), AJOH(0:N), AJV(0:N), AJRF(0:N),  AJNB(0:N),  stat = ierl(1))
       allocate(AJBS(0:N),  AJBS1(0:N),  AJBS2(0:N),  AJBS3(0:N), AJBS4(0:N), stat = ierl(2))
       allocate(ETA(0:N),ETAS(0:N),ETA1(0:N),ETA2(0:N),ETA3(0:N),ETA4(0:N),   stat = ierl(3))
       allocate(AJPARA(0:N),                                                  stat = ierl(4))
       ier = sum(ierl) ; iflag = 7
       if (ier /= 0) exit

       allocate(ANS0(1:NS), TS0(1:NS),   ANSAV(1:NS), TSAV(1:NS), WST(1:NS),  stat = ierl(1))
       allocate(ANF0(1:NF), TF0(1:NF),   ANFAV(1:NF), TFAV(1:NF), WFT(1:NF),  stat = ierl(2))
       allocate(PBCLT(1:NS),PFCLT(1:NS), PLT(1:NS),   SPET(1:NS), SLT(1:NS),  stat = ierl(3))
       allocate(Deff(0:N),  thrp(1:2*N), Ueth_NC(0:N,1:2),Uith_NC(0:N,1:2),   stat = ierl(4))
       ier = sum(ierl) ; iflag = 8
       if (ier /= 0) exit

       allocate(ALC(0:NCM,1:NQM,0:N),BLC(0:NCM,1:NQM,0:N),CLC(0:NCM,1:NQM,0:N),stat = ierl(1))
       allocate(PLC(1:NCM,1:NQM,0:N),                                          stat = ierl(2))
       allocate(NLC(0:NCM,1:NQM),NLCR(0:NCM,1:NQM,0:N),NLCMAX(1:NQM),          stat = ierl(3))
       allocate(X(1:NQM,0:N), XOLD(1:NQM,0:N),                                 stat = ierl(4))
       ier = sum(ierl) ; iflag = 9
       if (ier /= 0) exit

       allocate(GX(0:N), GY(0:N,0:NGRM,1:NGYRM), GQY(0:N,1:NCM,1:NQM),         stat = ierl(1))
       allocate(GYT(0:N,0:NGTM,1:NGYRM),                                       stat = ierl(2))
       ier = sum(ierl) ; iflag = 10
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

    deallocate(R,      PSI,    RHO)
    deallocate(H,      HPSI)

    deallocate(ErV,    EthV,   EphV,  BthV,  BphV)
    deallocate(PNeV,   UerV,   UethV, UephV, PTeV)
    deallocate(PNiV,   UirV,   UithV, UiphV, PTiV)
    deallocate(PNbV,   UbthV,  UbphV, PN01V, PN02V, PN03V)
    deallocate(AphV,   PhiV,   RAthV, PeV,   PiV)
    deallocate(RUethV, RUithV, PT01V, PT02V, PT03V, PNbrpV)

    deallocate(PNeV_FIX, PTeV_FIX, PNiV_FIX, PTiV_FIX, ErV_FIX)

    deallocate(rNuION, rNu0e,  rNu0i, rNu0b, rNuL)
    deallocate(rNuiCX, rNubCX, rNuee, rNuei, rNuie, rNuiCXT)
    deallocate(rNuii,  rNuTei, rNube, rNubi, rNuD)
    deallocate(rNubrp1,rNubrp2,rNuei1,rNuei2,rNuei3,rNuei2Bth)
    deallocate(rNube1, rNube2, rNube3,rNuLTe,rNuLTi,rNube2Bth)
    deallocate(rNueNC, rNuiNC, rNuAse,rNuAsi)
!    deallocate(rNueHL, rNuiHL, rNueHLthth, rNuiHLthth, rNueHLthph, &
!         &     rNuiHLthph, rNueHLphth, rNuiHLphth, rNueHLphph, &
!         &     rNuiHLphph)
!    deallocate(HPN, EpsHM) ! miki_m 10-08-06
    deallocate(rNueHLM, rNuiHLM, &
         &     rNueHLthth, rNuiHLthth, rNueHLthph, rNuiHLthph, &
         &     rNueHLphth, rNuiHLphth, rNueHLphph, rNuiHLphph, &
         &     rNueHLththM, rNueHLthphM, rNueHLphthM, rNueHLphphM, &
         &     rNuiHLththM, rNuiHLthphM, rNuiHLphthM, rNuiHLphphM) ! miki_m 10-08-06
    deallocate(FWe,    FWi,    FWthe,  FWthi,  WPM,   FVpch)
    deallocate(rMue,  rMui)
    deallocate(rNuB,   rNuLB,  ft,    Chie,  Chii)
    deallocate(De,     Di,     VWpch, D01,   D02,  D03)
    deallocate(ChiNCpe,ChiNCte,ChiNCpi,ChiNCti)
    deallocate(FWthphe,FWthphi)
    deallocate(Ubrp,   RUbrp,  Dbrp,  DltRP, DltRP_mid, rNubL)
    deallocate(rip_rat,rNuOL,  rNuNTV,UastNC,Vbpara)
    deallocate(Fmnq,   Wnm,    Umnq)
    deallocate(deltam, gamITG)
    deallocate(DMAG,   DMAGe,  DMAGi)  !***AF (2008-06-08)
    deallocate(SiVizA, SiVcxA, wexb,  Ys)
    deallocate(rNue2NC,rNui2NC,FQeth, FQith)

    deallocate(rG1h2,  FCDBM,  S,     Alpha, rKappa)
    deallocate(pres0,  ErV0)

    deallocate(rG1h2IM, FCDIM, RAQPR)  !***miki_m 09/06/17~

    deallocate(PNB,    PNBTG, PNBPD, PNBcol_e, PNBcol_i)
    deallocate(SNB,    SNBe,  SNBi,  SNBb)
    deallocate(SNBPDi, SNBTGi,Tqt,   Tqp)
    deallocate(PNBe,   PNBi,  MNB,   PRFe,   PRFi)
    deallocate(POH,    POHe,  POHi,  PEQe,   PEQi)
    deallocate(SiLC,   SiLCth,SiLCph,PALFe,  PALFi)
    deallocate(PIE,    PCX,   SIE,   SCX,    PBr)

    deallocate(Q, AJ, AJOH, AJV, AJRF, AJNB, AJPARA)
    deallocate(AJBS, AJBS1, AJBS2, AJBS3, AJBS4)
    deallocate(ETA, ETAS, ETA1, ETA2, ETA3, ETA4)

    deallocate(ANS0, TS0, ANSAV, TSAV, WST)
    deallocate(ANF0, TF0, ANFAV, TFAV, WFT)
    deallocate(PBCLT, PFCLT, PLT, SPET, SLT)
    deallocate(Deff, thrp, Ueth_NC, Uith_NC)

    deallocate(ALC, BLC, CLC)
    deallocate(PLC)
    deallocate(NLC,NLCR,NLCMAX)
    deallocate(X, XOLD)

    deallocate(GX, GY, GQY)
    deallocate(GYT)

    if(allocated(infiles)) deallocate(infiles)

  end subroutine deallocate_txcomm

end module tx_commons
