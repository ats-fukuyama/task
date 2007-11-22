module commons
  implicit none
  public

  integer(4), parameter :: NRM=101, NEM=NRM, NQM=21, NCM=29, NGRM=20, &
       &                   NGTM=5000, NGVM=5000, NGYRM=120, NGYTM=46, &
       &                   NGYVM=49, NGPRM=18, NGPTM=8, NGPVM=15
  integer(4), parameter :: NSM=2, NFM=2
  integer(4), parameter :: LQm1=1,  LQm2=2,  LQm3=3,  LQm4=4,  LQm5=5,&
       &                   LQe1=6,  LQe2=7,  LQe3=8,  LQe4=9,  LQe5=10,&
       &                   LQi1=11, LQi2=12, LQi3=13, LQi4=14, LQi5=15,&
       &                   LQb1=16,          LQb3=17, LQb4=18,&
       &                   LQn1=19, LQn2=20, &
       &                   LQr1=21

  !**********************!
  !   INPUT PARAMETERS   !
  !**********************!

  ! Configuration parameters
  real(8) :: RA, RB, RC, RR, BB, rIPs, rIPe

  ! Species
  real(8) :: PA, PZ, Zeff

  ! Profiles
  real(8) :: PN0, PNa, PTe0, PTea, PTi0, PTia, PROFJ

  ! Diffusivities and viscosities
  real(8) :: De0, Di0, rMue0, rMui0, Chie0, Chii0, WPM0, WPE0, WPI0

  ! Amplitude parameters for transport
  real(8) :: FSDFIX, FSCDBM, FSBOHM, FSPSCL, PROFD
  real(8) :: FSCX, FSLC, FSNC, FSLP, FSLTE, FSLTI, FSION, FSD0, rG1, FSRP
  integer(4) :: MDLC

  ! Scale lengths in SOL
  real(8) :: rLn, rLT

  ! Heat sources
  real(8) :: Eb, RNBP, RNBP0, RNBT1, RNBT2, RNBT10, RNBT20, PNBH, PNBHP, PNBHT1, PNBHT2, &
       &     PNBCD, rNRF, RRF, RRF0, PRFH

  ! Neutral parameters
  real(8) :: PN0s, V0, rGamm0, rGASPF, PNeDIV, PTeDIV, PTiDIV

  ! Ripple parameters
  real(8) :: DIN, DltRP0
  integer(4) :: NTCOIL

  ! Numerical parameters
  real(8) :: DT, EPS, ADV, tiny_cap, CMESH0, WMESH0, CMESH, WMESH
  integer(4) :: ICMAX

  ! Mesh parameters
  integer(4) :: NRMAX, NEMAX, NTMAX, NTSTEP, NGRSTP, NGTSTP, NGVSTP, NRA, NRC

  ! Parameters for parameter survey
  real(8) :: DelR, DelN

  ! Helical parameters
  real(8) :: EpsH, FSHL, Q0, QA
  integer(4) :: NCPHI

  ! Diagnostic message
  integer(4) :: IDIAG

  ! LAPACK
  integer(4) :: MDLPCK

  ! Implicit Euler or Gear's backward differential formula
  integer(4) :: IGBDF

  !  Transport model
  integer(4) :: MDLWTB, MDLETA, MDFIXT

  !**********************************!
  !   INTERNAL CONTOROL PARAMETERS   !
  !**********************************!

  ! Configuration parameters
  integer(4) :: NT, NQMAX, IERR
  real(8) :: T_TX, TMAX
  real(8) :: AMI, AMB, Vb, sqeps0
  real(8) :: rIP, Bthb
  real(8) :: UHth, UHph
  real(8), dimension(:), allocatable :: R, PSI
  real(8), dimension(:), allocatable :: H, HPSI

  ! Convergence accelerator
  real(8) :: rMUb1, rMUb2, AMPe4 

  ! Variables
  real(8), dimension(:), allocatable :: &
       & ErV,    EthV,   EphV,  BthV,  BphV, &
       & PNeV,   UerV,   UethV, UephV, PTeV, &
       & PNiV,   UirV,   UithV, UiphV, PTiV, &
       & PNbV,   UbthV,  UbphV, PN01V, PN02V, &
       & AphV,   PhiV,    RAthV, PeV,   PiV,  &
       & RUethV, RUithV, PT01V, PT02V, PNbrpV
!!rp_conv       &, PNbrpLV

  real(8), dimension(:), allocatable :: PNeV_FIX, PTeV_FIX, dPNeV_FIX, dPNiV_FIX

  ! Coefficients
  real(8), dimension(:), allocatable :: &
       & rNuION, rNu0e, rNu0i, rNu0b, rNuL, rNuiCX, rNubCX, &
       & rNuee, rNuei, rNuie, rNuii, rNuTei, rNube, rNubi, rNuD, &
       & rNubrp1, rNubrp2, rNuei1, rNuei2, rNuei3, &
       & rNube1, rNube2, rNube3, rNuLTe, rNuLTi, &
       & rNueNC, rNuiNC, rNuAse, rNuAsi, rNueHL, rNuiHL, &
       & FWthe, FWthi, WPM, rMue, rMui, rNuB, rNuLB, ft, &
       & Chie, Chii, De, Di, D01, D02, &
       & WNthe, WEMthe, WWthe, WT1the, WT2the, &
       & WNthi, WEMthi, WWthi, WT1thi, WT2thi, &
       & FWthphe, FWthphi, rlnLe, rlnLi, &
       & Ubrp, RUbrp, Dbrp, DltRP, rNubL, RATIO, rNuOL
 
  ! CDBM
  real(8), dimension(:), allocatable :: rG1h2, FCDBM, S, Alpha, rKappa

  ! Sources and sinks
  real(8), dimension(:), allocatable :: PNB, PNBTG, PNBPD, PNBcol_e, PNBcol_i,  &
       &                                SNB, SNBTG, SNBPD, MNB, &
       &                                PNBe, PNBi, PRFe, PRFi, &
       &                                POH, POHe, POHi, PEQe, PEQi, &
       &                                SiLC, SiLCth, SiLCph
  real(8), dimension(:), allocatable :: PIE, PCX, SIE, PBr
  real(8) :: RatCX

  ! Safety factor, currents, resistivity
  real(8), dimension(:), allocatable :: Q, AJ, AJOH, AJV, AJRF, AJNB, &
       &                                AJBS, AJBS1, AJBS2, AJBS3, AJBS4, &
       &                                ETA, ETAS, ETA1, ETA2, ETA3, ETA4

  ! Global parameters for display
  real(8), dimension(:), allocatable :: ANS0, TS0, ANSAV, TSAV, WST
  real(8), dimension(:), allocatable :: ANF0, TF0, ANFAV, TFAV, WFT
  real(8), dimension(:), allocatable :: PBCLT, PFCLT, PLT, SPET, SLT
  real(8), dimension(:), allocatable :: Deff
  real(8) :: WBULKT, WTAILT, WPT
  real(8) :: AJT, AJOHT, AJNBT, AJRFT, AJBST
  real(8) :: PINT, POHT, PNBT, PRFT, PNFT
  real(8) :: PBINT, PFINT, POUT, PCXT, PIET, PRLT, SINT, SIET
  real(8) :: SNBT, SNFT, SOUT
  real(8) :: VLOOP, ALI, RQ1, RPE, ZEFF0, QF
  real(8) :: WPDOT, TAUE1, TAUE2, TAUEP, TAUEH
  real(8) :: BETAP0, BETAPA, BETA0, BETAA, BETAQ0, BETAN
  real(8) :: TPRE, WPPRE
  real(8) :: VOLAVN

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

contains

  subroutine allocate_txcomm(ier, icont)

    integer(4), intent(out) :: ier
    integer(4), intent(in), optional :: icont
    integer(4) :: iflag, N, NS, NF
    integer(4), dimension(1:15) :: ierl

    ierl(1:15) = 0
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
       if(present(icont) .and. icont /= 0) then
          call deallocate_txcomm
       else
          return
       end if
    end if

    do
       allocate(R(0:N),      PSI(0:N),                                        stat = ierl(1))
       allocate(H(0:NEMAX),  HPSI(0:NEMAX),                                   stat = ierl(2))
       ier = sum(ierl) ; iflag = 1
       if (ier /= 0) exit

       allocate(ErV(0:N),    EthV(0:N),   EphV(0:N),  BthV(0:N),  BphV(0:N),  stat = ierl(1))
       allocate(PNeV(0:N),   UerV(0:N),   UethV(0:N), UephV(0:N), PTeV(0:N),  stat = ierl(2))
       allocate(PNiV(0:N),   UirV(0:N),   UithV(0:N), UiphV(0:N), PTiV(0:N),  stat = ierl(3))
       allocate(PNbV(0:N),   UbthV(0:N),  UbphV(0:N), PN01V(0:N), PN02V(0:N), stat = ierl(4))
       allocate(AphV(0:N),   PhiV(0:N),   RAthV(0:N), PeV(0:N),   PiV(0:N),   stat = ierl(5))
       allocate(RUethV(0:N), RUithV(0:N), PT01V(0:N), PT02V(0:N), PNbrpV(0:N),stat = ierl(6))
       ier = sum(ierl) ; iflag = 2
       if (ier /= 0) exit

       allocate(PNeV_FIX(0:N), PTeV_FIX(0:N), dPNeV_FIX(0:N), dPNiV_FIX(0:N), stat = ierl(1))
       ier = sum(ierl) ; iflag = 3
       if (ier /= 0) exit

       allocate(rNuION(0:N), rNu0e(0:N),  rNu0i(0:N), rNu0b(0:N), rNuL(0:N),  stat = ierl(1))
       allocate(rNuiCX(0:N), rNubCX(0:N), rNuee(0:N), rNuei(0:N), rNuie(0:N), stat = ierl(2))
       allocate(rNuii(0:N),  rNuTei(0:N), rNube(0:N), rNubi(0:N), rNuD(0:N),  stat = ierl(3))
       allocate(rNubrp1(0:N),rNubrp2(0:N),rNuei1(0:N),rNuei2(0:N),rNuei3(0:N),stat = ierl(4))
       allocate(rNube1(0:N), rNube2(0:N), rNube3(0:N),rNuLTe(0:N),rNuLTi(0:N),stat = ierl(5))
       allocate(rNueNC(0:N), rNuiNC(0:N), rNuAse(0:N),rNuAsi(0:N),            stat = ierl(6))
       allocate(rNueHL(0:N), rNuiHL(0:N),                                     stat = ierl(7))
       allocate(FWthe(0:N),  FWthi(0:N),  WPM(0:N),   rMue(0:N),  rMui(0:N),  stat = ierl(8))
       allocate(rNuB(0:N),   rNuLB(0:N),  ft(0:N),    Chie(0:N),  Chii(0:N),  stat = ierl(9))
       allocate(De(0:N),     Di(0:N),     D01(0:N),   D02(0:N),               stat = ierl(10))
       allocate(WNthe(0:N),  WEMthe(0:N), WWthe(0:N), WT1the(0:N),WT2the(0:N),stat = ierl(11))
       allocate(WNthi(0:N),  WEMthi(0:N), WWthi(0:N), WT1thi(0:N),WT2thi(0:N),stat = ierl(12))
       allocate(FWthphe(0:N),FWthphi(0:N),rlnLe(0:N), rlnLi(0:N),             stat = ierl(13))
       allocate(Ubrp(0:N),   RUbrp(0:N),  Dbrp(0:N),  DltRP(0:N), rNubL(0:N), stat = ierl(14))
       allocate(RATIO(0:N),  rNuOL(0:N),                                      stat = ierl(15))
       ier = sum(ierl) ; iflag = 4
       if (ier /= 0) exit

       allocate(rG1h2(0:N),  FCDBM(0:N),  S(0:N),     Alpha(0:N), rKappa(0:N),stat = ierl(1))
       ier = sum(ierl) ; iflag = 5
       if (ier /= 0) exit

       allocate(PNB(0:N),   PNBTG(0:N),PNBPD(0:N),PNBcol_e(0:N),PNBcol_i(0:N),stat = ierl(1))
       allocate(SNB(0:N),   SNBTG(0:N),  SNBPD(0:N),  MNB(0:N),               stat = ierl(2))
       allocate(PNBe(0:N),  PNBi(0:N),   PRFe(0:N),   PRFi(0:N),              stat = ierl(3))
       allocate(POH(0:N),   POHe(0:N),   POHi(0:N),   PEQe(0:N),  PEQi(0:N),  stat = ierl(4))
       allocate(SiLC(0:N),  SiLCth(0:N), SiLCph(0:N),                         stat = ierl(5))
       allocate(PIE(0:N),   PCX(0:N),    SIE(0:N),    PBr(0:N),               stat = ierl(6))
       ier = sum(ierl) ; iflag = 6
       if (ier /= 0) exit

       allocate(Q(0:N), AJ(0:N), AJOH(0:N), AJV(0:N), AJRF(0:N),  AJNB(0:N),  stat = ierl(1))
       allocate(AJBS(0:N),  AJBS1(0:N),  AJBS2(0:N),  AJBS3(0:N), AJBS4(0:N), stat = ierl(2))
       allocate(ETA(0:N),ETAS(0:N),ETA1(0:N),ETA2(0:N),ETA3(0:N),ETA4(0:N),   stat = ierl(3))
       ier = sum(ierl) ; iflag = 7
       if (ier /= 0) exit

       allocate(ANS0(1:NS), TS0(1:NS),   ANSAV(1:NS), TSAV(1:NS), WST(1:NS),  stat = ierl(1))
       allocate(ANF0(1:NF), TF0(1:NF),   ANFAV(1:NF), TFAV(1:NF), WFT(1:NF),  stat = ierl(2))
       allocate(PBCLT(1:NS),PFCLT(1:NS), PLT(1:NS),   SPET(1:NS), SLT(1:NS),  stat = ierl(3))
       allocate(Deff(0:N),                                                    stat = ierl(4))
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
       call deallocate_txcomm
       write(6,*) "XX Allocation error category = ",iflag
    end if

  end subroutine allocate_txcomm

  subroutine deallocate_txcomm

    deallocate(R,      PSI)
    deallocate(H,      HPSI)

    deallocate(ErV,    EthV,   EphV,  BthV,  BphV)
    deallocate(PNeV,   UerV,   UethV, UephV, PTeV)
    deallocate(PNiV,   UirV,   UithV, UiphV, PTiV)
    deallocate(PNbV,   UbthV,  UbphV, PN01V, PN02V)
    deallocate(AphV,   PhiV,   RAthV, PeV,   PiV)
    deallocate(RUethV, RUithV, PT01V, PT02V, PNbrpV)

    deallocate(PNeV_FIX, PTeV_FIX)

    deallocate(rNuION, rNu0e,  rNu0i, rNu0b, rNuL)
    deallocate(rNuiCX, rNubCX, rNuee, rNuei, rNuie)
    deallocate(rNuii,  rNuTei, rNube, rNubi, rNuD)
    deallocate(rNubrp1,rNubrp2,rNuei1,rNuei2,rNuei3)
    deallocate(rNube1, rNube2, rNube3,rNuLTe,rNuLTi)
    deallocate(rNueNC, rNuiNC, rNuAse,rNuAsi)
    deallocate(rNueHL, rNuiHL)
    deallocate(FWthe,  FWthi,  WPM,   rMue,  rMui)
    deallocate(rNuB,   rNuLB,  ft,    Chie,  Chii)
    deallocate(De,     Di,     D01,   D02)
    deallocate(WNthe,  WEMthe, WWthe, WT1the,WT2the)
    deallocate(WNthi,  WEMthi, WWthi, WT1thi,WT2thi)
    deallocate(FWthphe,FWthphi,rlnLe, rlnLi)
    deallocate(Ubrp,   RUbrp,  Dbrp,  DltRP, rNubL)
    deallocate(RATIO,  rNuOL)

    deallocate(rG1h2,  FCDBM,  S,     Alpha, rKappa)

    deallocate(PNB,    PNBTG, PNBPD, PNBcol_e, PNBcol_i)
    deallocate(SNB,    SNBTG, SNBPD, MNB)
    deallocate(PNBe,   PNBi,  PRFe,  PRFi)
    deallocate(POH,    POHe,  POHi,  PEQe,   PEQi)
    deallocate(SiLC,   SiLCth,SiLCph)
    deallocate(PIE,    PCX,   SIE,   PBr)

    deallocate(Q, AJ, AJOH, AJV, AJRF, AJNB)
    deallocate(AJBS, AJBS1, AJBS2, AJBS3, AJBS4)
    deallocate(ETA, ETAS, ETA1, ETA2, ETA3, ETA4)

    deallocate(ANS0, TS0, ANSAV, TSAV, WST)
    deallocate(ANF0, TF0, ANFAV, TFAV, WFT)
    deallocate(PBCLT, PFCLT, PLT, SPET, SLT)
    deallocate(Deff)

    deallocate(ALC, BLC, CLC)
    deallocate(PLC)
    deallocate(NLC,NLCR,NLCMAX)
    deallocate(X, XOLD)

    deallocate(GX, GY, GQY)
    deallocate(GYT)

  end subroutine deallocate_txcomm

end module commons
