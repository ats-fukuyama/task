module commons
  implicit none
  public

  integer, parameter :: NRM=101, NEM=NRM, NQM=20, NCM=21, NGRM=20, &
       &                NGTM=1000, NGVM=1000, NGYRM=80, NGYTM=39, &
       &                NGYVM=43, NGPRM=14, NGPTM=13, NGPVM=7
  integer, parameter :: NSM=2, NFM=2
  integer, parameter :: LQm1=1,  LQm2=2,  LQm3=3,  LQm4=4,  LQm5=5,&
       &                LQe1=6,  LQe2=7,  LQe3=8,  LQe4=9,  LQe5=10,&
       &                LQi1=11, LQi2=12, LQi3=13, LQi4=14, LQi5=15,&
       &                LQb1=16,          LQb3=17, LQb4=18,&
       &                LQn1=19, LQn2=20

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
  real(8) :: De0, Di0, rMue0, rMui0, WPM0, Chie0, Chii0

  ! Amplitude parameters for transport
  real(8) :: FSDFIX, FSCDBM, FSBOHM, FSPSCL, PROFD
  real(8) :: FSCX, FSLC, FSNC, FSLP, FSION, FSD0, rG1

  ! Scale lengths in SOL
  real(8) :: rLn, rLT

  ! Heat sources
  real(8) :: Eb, RNB, PNBH, PNBCD, rNRF, RRF, PRFH

  ! Neutral parameters
  real(8) :: PN0s, V0, rGamm0, rGASPF, PNeDIV, PTeDIV, PTiDIV

  ! Numerical parameters
  real(8) :: DLT, DT, EPS
  integer :: ICMAX

  ! Mesh parameters
  integer :: NRMAX, NEMAX, NTMAX, NTSTEP, NGRSTP, NGTSTP, NGVSTP, NRA

  ! Parameters for parameter survey
  real(8) :: DelR, DelN

  ! Helical parameters
  real(8) :: EpsH, FSHL, Q0, QA
  integer :: NCPHI

  !**********************************!
  !   INTERNAL CONTOROL PARAMETERS   !
  !**********************************!

  ! Configuration parameters
  integer :: NT, NQMAX, IERR
  real(8) :: T_TX, TMAX
  real(8) :: AMI, AMB, Vb
  real(8) :: rIP, Bthb
  real(8) :: UHth, UHph
  real(8), dimension(0:NRM) :: R
  real(8), dimension(0:NEM) :: H
  
  ! Variables
  real(8), dimension(0:NRM) :: ErV, EthV, EphV, BthV, BphV, &
       &                       PNeV, UerV, UethV, UephV, PTeV, &
       &                       PNiV, UirV, UithV, UiphV, PTiV, &
       &                       PNbV, UbthV, UbphV, PN01V, PN02V

  ! Derivative
  real(8), dimension(0:NRM) :: RdBthV

  ! Coefficients
  real(8), dimension(0:NRM) :: rNuION, rNu0e, rNu0i, rNuL, rNuiCX, &
       &                       rNuei, rNuii, rNuTei, rNube, rNubi, &
       &                       rNueNC, rNuiNC, rNueHL, rNuiHL, &
       &                       FWthe, FWthi, WPM, rMue, rMui, rNuB, &
       &                       Chie, Chii, De, Di, D01, D02

  ! CDBM
  real(8), dimension(0:NRM) :: rG1h2, FCDBM, S, Alpha, rKappa

  ! Sources and sinks
  real(8), dimension(0:NRM) :: PNB, SNB, PRFe, PRFi, POH, SiLC, SiLCth, SiLCph
  real(8), dimension(0:NRM) :: PIE, PCX, SIE

  ! Safety factor, currents
  real(8), dimension(0:NRM) :: Q, AJ, AJOH, AJV, AJRF, AJNB, AJBS

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
  real(8) :: WPDOT, TAUE1, TAUE2, TAUEP
  real(8) :: BETAP0, BETAPA, BETA0, BETAA, BETAQ0
  real(8) :: TPRE, WPPRE

  ! Internal variables for transport matrix
  real(8), dimension(0:NCM,1:NQM,0:NRM) :: ALC, BLC, CLC
  real(8), dimension(1:NCM,1:NQM,0:NRM) :: PLC
  integer, dimension(0:NCM,1:NQM)       :: NLC
  integer, dimension(0:NCM,1:NQM,0:NRM) :: NLCR
  integer, dimension(1:NQM)             :: NLCMAX
  real(8), dimension(1:4*NQM-1,1:NQM*(NRM+1)) :: BA
  real(8), dimension(1:NQM*(NRM+1)) :: BX
  real(8), dimension(1:NQM,0:NRM) :: X
  
  ! Graphics
  integer :: NGR, NGT, NGVV
  integer :: MODEG, MODEAV, MODEL
  real, dimension(0:NRM) :: GX
  real, dimension(0:NRM,0:NGRM,1:NGYRM) :: GY
  real, dimension(0:NGRM) :: GT
  real, dimension(0:NGTM) :: GTX
  real, dimension(0:NGTM,1:NGYTM) :: GTY
  real, dimension(1:NGYRM) :: gDIV
  real, dimension(0:NGVM) :: GVX
  real, dimension(0:NGVM,1:NGYVM) :: GVY

  ! I/O
  character(len=20) :: SLID

end module commons
