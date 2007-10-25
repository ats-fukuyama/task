module commons
  implicit none
  public

  integer, parameter :: NRM=101, NEM=NRM, NQM=21, NCM=29, NGRM=20, &
       &                NGTM=5000, NGVM=5000, NGYRM=118, NGYTM=46, &
       &                NGYVM=49, NGPRM=18, NGPTM=8, NGPVM=15
  integer, parameter :: NSM=2, NFM=2
  integer, parameter :: LQm1=1,  LQm2=2,  LQm3=3,  LQm4=4,  LQm5=5,&
       &                LQe1=6,  LQe2=7,  LQe3=8,  LQe4=9,  LQe5=10,&
       &                LQi1=11, LQi2=12, LQi3=13, LQi4=14, LQi5=15,&
       &                LQb1=16,          LQb3=17, LQb4=18,&
       &                LQn1=19, LQn2=20, &
       &                LQr1=21

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
  integer :: MDLC

  ! Scale lengths in SOL
  real(8) :: rLn, rLT

  ! Heat sources
  real(8) :: Eb, RNBP, RNBP0, RNBT1, RNBT2, RNBT10, RNBT20, PNBH, PNBHP, PNBHT1, PNBHT2, &
       &     PNBCD, rNRF, RRF, RRF0, PRFH

  ! Neutral parameters
  real(8) :: PN0s, V0, rGamm0, rGASPF, PNeDIV, PTeDIV, PTiDIV

  ! Ripple parameters
  real(8) :: DIN, DltRP0
  integer :: NTCOIL

  ! Numerical parameters
  real(8) :: DT, EPS, ADV, tiny_cap, CMESH, WMESH
  integer :: ICMAX

  ! Mesh parameters
  integer :: NRMAX, NEMAX, NTMAX, NTSTEP, NGRSTP, NGTSTP, NGVSTP, NRA, NRC

  ! Parameters for parameter survey
  real(8) :: DelR, DelN

  ! Helical parameters
  real(8) :: EpsH, FSHL, Q0, QA
  integer :: NCPHI

  ! Diagnostic message
  integer :: IDIAG

  ! LAPACK
  integer :: MDLPCK

  ! Implicit Euler or Gear's backward differential formula
  integer :: IGBDF

  !  Transport model
  integer :: MDLWTB, MDLETA, MDFIXT

  !**********************************!
  !   INTERNAL CONTOROL PARAMETERS   !
  !**********************************!

  ! Configuration parameters
  integer :: NT, NQMAX, IERR
  real(8) :: T_TX, TMAX
  real(8) :: AMI, AMB, Vb, sqeps0
  real(8) :: rIP, Bthb
  real(8) :: UHth, UHph
  real(8), dimension(0:NRM) :: R, PSI
  real(8), dimension(0:NEM) :: H, HPSI

  ! Convergence accelerator
  real(8) :: rMUb1, rMUb2, AMPe4 

  ! Variables
  real(8), dimension(0:NRM) :: ErV,    EthV,   EphV,  BthV,  BphV, &
       &                       PNeV,   UerV,   UethV, UephV, PTeV, &
       &                       PNiV,   UirV,   UithV, UiphV, PTiV, &
       &                       PNbV,   UbthV,  UbphV, PN01V, PN02V, &
       &                       AphV,   Phi,    RAthV, PeV,   PiV,  &
       &                       RUethV, RUithV, PT01V, PT02V, PNbrpV, PNbrpLV

  real(8), dimension(0:NRM) :: PNeV_FIX, PTeV_FIX, dPNeV_FIX, dPNiV_FIX

  ! Coefficients
  real(8), dimension(0:NRM) :: rNuION, rNu0e, rNu0i, rNu0b, rNuL, rNuiCX, rNubCX, &
       &                       rNuee, rNuei, rNuie, rNuii, rNuTei, rNube, rNubi, &
       &                       rNubrp1, rNubrp2, rNuD, rNuei1, rNuei2, rNuei3, &
       &                       rNube1, rNube2, rNube3, &
       &                       rNueNC, rNuiNC, rNuAse, rNuAsi, rNueHL, rNuiHL, &
       &                       FWthe, FWthi, WPM, rMue, rMui, rNuB, rNuLB, ft, &
       &                       Chie, Chii, De, Di, D01, D02, rNuLTe, rNuLTi, &
       &                       WNthe, WEMthe, WWthe, WT1the, WT2the, &
       &                       WNthi, WEMthi, WWthi, WT1thi, WT2thi, &
       &                       FWthphe, FWthphi, rlnLe, rlnLi, &
       &                       Ubrp, Dbrp, DltRP, rNubL, RATIO, rNuOL
  real(8) :: FWthea, FWthia
 
  ! CDBM
  real(8), dimension(0:NRM) :: rG1h2, FCDBM, S, Alpha, rKappa

  ! Sources and sinks
  real(8), dimension(0:NRM) :: PNB, PNBPD, PNBcol_e, PNBcol_i, PNBTG, SNB, SNBTG, SNBPD, &
       &                       MNB, PNBe, PNBi, PRFe, PRFi, &
       &                       POH, POHe, POHi, PEQe, PEQi, SiLC, SiLCth, SiLCph, RUiloss
  real(8), dimension(0:NRM) :: PIE, PCX, SIE, PBr
  real(8) :: RatCX

  ! Safety factor, currents, resistivity
  real(8), dimension(0:NRM) :: Q, AJ, AJOH, AJV, AJRF, AJNB, AJBS, AJBS1, AJBS2, AJBS3, &
       &                       AJBS4, ETA, ETAS, ETA1, ETA2, ETA3, ETA4

  ! Global parameters for display
  real(8), dimension(1:NSM) :: ANS0, TS0, ANSAV, TSAV, WST
  real(8), dimension(1:NFM) :: ANF0, TF0, ANFAV, TFAV, WFT
  real(8), dimension(1:NSM) :: PBCLT, PFCLT, PLT, SPET, SLT
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
  real(8), dimension(0:NCM,1:NQM,0:NRM) :: ALC, BLC, CLC
  real(8), dimension(1:NCM,1:NQM,0:NRM) :: PLC
  integer, dimension(0:NCM,1:NQM)       :: NLC
  integer, dimension(0:NCM,1:NQM,0:NRM) :: NLCR
  integer, dimension(1:NQM)             :: NLCMAX
  real(8), dimension(1:NQM,0:NRM) :: X
  
  ! Graphics
  integer :: NGR, NGT, NGVV
  integer :: MODEG, MODEAV, MODEGL
  real, dimension(0:NRM) :: GX
  real, dimension(0:NRM,0:NGRM,1:NGYRM) :: GY
  real, dimension(0:NGRM) :: GT
  real, dimension(0:NGTM) :: GTX
  real, dimension(0:NGTM,1:NGYTM) :: GTY
  real, dimension(1:NGYRM) :: gDIV
  real, dimension(0:NGVM) :: GVX
  real, dimension(0:NGVM,1:NGYVM) :: GVY
  real, dimension(0:NRM,1:NCM,1:NQM) :: GQY

  ! I/O
  character(len=20) :: SLID

end module commons
