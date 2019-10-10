MODULE TRCOMM
  USE TRCOM0
  IMPLICIT NONE

!     ****** CONSTANTS, based on CODATA 2006 ******
! TRCNS
!        PI    : Pi
!        AEE   : Elementaty charge
!        AME   : Electron mass
!        AMM   : Proton mass
!        VC    : Speed of light in vacuum
!        RMU0  : Permeability of free space
!        EPS0  : Permittivity of free space
!        RKEV  : Factor ([keV] -> [J])
  REAL(8), PARAMETER :: PI   = 3.14159265358979323846D0
  REAL(8), PARAMETER :: AEE  = 1.602176487D-19
  REAL(8), PARAMETER :: AME  = 9.10938215D-31
  REAL(8), PARAMETER :: AMM  = 1.672621637D-27
  REAL(8), PARAMETER :: VC   = 2.99792458D8
  REAL(8), PARAMETER :: RMU0 = 4.D0*PI*1.D-7
  REAL(8), PARAMETER :: EPS0 = 1.D0/(VC*VC*RMU0)
  REAL(8), PARAMETER :: RKEV = AEE*1.D3

  INTEGER, PARAMETER :: NPSCM = 10 ! Maximum number of particle source

!     ****** PARAMETERS ******
! TRPRM
  REAL(8)   :: &
       RR, RA, RKAP, RDLT, BB, RIPS, RIPE, RIPSS, PHIA, PNC, PNFE, PNNU, &
       PNNUS, PROFN1, PROFN2, PROFT1, PROFT2, PROFU1, PROFU2, &
       PROFNU1, PROFNU2, PROFJ1, &
       PROFJ2, AD0, AV0, CNP, CNH, CDP, CDH, CNN, CWEB, DT, EPSLTR, &
       CHP, CK0, CK1, CKALFA, CKBETA, CKGUMA, CNB, CALF, CSPRS, TSST, &
       SYNCABS, SYNCSELF
  REAL(8), DIMENSION(3) :: &
       ALP
  REAL(8), DIMENSION(8) :: &
       CDW
  REAL(8), DIMENSION(NSMM) :: &
       PA,PZ,PN,PNS,PT,PTS,PU,PUS
  INTEGER, DIMENSION(NSMM) :: &
       NPA
  INTEGER(4):: &
       LMAXTR, MDLKAI, MDLETA, MDLAD, MDLAVK, MDLKNC, MDLTPF, &
       NTMAX, NTSTEP, NGTSTP, NGRSTP, NGPST, MODELG,MDLDSK,MDTC, &
       MDLPR

!     ****** MODEL PARAMETERS ******
! TRMDL
  REAL(8)   :: &
       TPRST, PBSCD, &
       PNBTOT, PNBR0, PNBRW, PNBCD, PNBVY, PNBVW, PNBENG, PNBRTG, &
       PECTOT, PECR0, PECRW, PECCD, PECTOE, PECNPR, &
       PLHTOT, PLHR0, PLHRW, PLHCD, PLHTOE, PLHNPR, &
       PICTOT, PICR0, PICRW, PICCD, PICTOE, PICNPR, &
       PELTOT, PELR0, PELRW, PELRAD, PELVEL, PELTIM, &
       ELMWID, ELMDUR
  REAL(8), DIMENSION(NSMM) :: &
       PELPAT, ELMNRD, ELMTRD, ELMENH
  REAL(8), DIMENSION(NPSCM) :: &
       PSCTOT,PSCR0,PSCRW
  INTEGER, DIMENSION(NPSCM) :: &
       NSPSC
  INTEGER(4):: &
       MDLNB, MDLEC, MDLLH, MDLIC, MDLCD, MDLPEL, MDLJBS, MDLST, MDLNF, &
       IZERO, MDLELM, MDLPSC, NPSCMAX

!     ****** CONTROL VARIABLES ******
!TRCTL
  REAL(4)   :: &
       GTCPU1
  REAL(8)   :: &
       T, TST, TPRE, WPPRE, RIP, DR, FKAP, RHOA, RIPA, VSEC, RDPS,DIPDT
  REAL(8), DIMENSION(:),POINTER :: & ! (NSTM)
       PNSS
  INTEGER(4):: &
       NT, NRAMAX, NROMAX, NREDGE, NTMAX_SAVE, IREAD

!     ****** MATRIX VARIBALES ******
! TRMTX
  REAL(8), DIMENSION(:,:),POINTER :: & ! (NVM,NRM)
       XV
  REAL(8), DIMENSION(:,:),POINTER :: & ! (NFM,NRM)
       YV, AY, Y
  REAL(8), DIMENSION(:,:),POINTER :: & ! (NSM,NRM)
       ZV, AZ, Z
  REAL(8), DIMENSION(:,:),POINTER :: & ! (LDAB,MLM)
       AX
  REAL(8), DIMENSION(:)  ,POINTER :: & ! (MLM)
       X
!
!     ****** FUNDAMENTAL VARIABLES ******
! TRVAR
  REAL(8), DIMENSION(:)  ,POINTER :: & ! (NRM)
       RG, RM, RHOM, RHOG, BP, RDP, RPSI
  REAL(8), DIMENSION(:,:),POINTER :: & ! (NRM,NSTM)
       RN, RT, RU
  REAL(8), DIMENSION(:,:),POINTER :: & ! (NRM,NFM)
       RW

!     ****** PROFILE VARIABLES ******
! TRPFL
  REAL(8), DIMENSION(:,:), POINTER :: & ! (NRM,NFM)
       RNF, RTF
  REAL(8), DIMENSION(:)  , POINTER :: & ! (NRM)
       ANC, ANFE, ANNU, ZEFF, PZC, PZFE, BETA, BETAP, BETAL, BETAPL, &
       BETAQ, PBM, PADD, VTOR, VPAR, VPRP, VPOL, WROT, ER, VEXB, WEXB, AGMP, &
       VEXBP, WEXBP

!     ****** SOURCE VARIABLES ******
! TRSRC
  REAL(8), DIMENSION(10)      :: &
       RTG
  REAL(8), DIMENSION(:),     POINTER :: & ! (NRM)
       AJ, AJOH, EZOH, QP, AJTOR, AJNB, AJRF, AJBS, QPINV, PNB, SNB, &
       PBIN, PNF, SNF, PFIN, POH, PRB, PRC, PRL, PRSUM, &
       PCX, PIE, SIE, SCX, TSIE, TSCX
  REAL(8), DIMENSION(:,:),   POINTER :: & ! (NRM,NSTM)
       PIN, SSIN, PBCL, PFCL, PRF, SPE
  REAL(8), DIMENSION(:,:),   POINTER :: & ! (NRM,NSM)
       RGFLX,SPSC
  REAL(8), DIMENSION(:,:),   POINTER :: & ! (NRM,3)
       AJRFV
  REAL(8), DIMENSION(:,:,:), POINTER :: & ! (NRM,NSTM,3)
       PRFV

!     ****** COEFFICIENT VARIABLES ******
! TRCEF
  REAL(8), DIMENSION(:)  ,POINTER :: & ! (NRM)
       ETA, S, ALPHA, RKCV, TAUB, TAUF, TAUK
  REAL(8), DIMENSION(:,:),POINTER :: & ! (NRM,NSTM)
       AK, AVK, AD, AV, AKNC, AKDW, ADNC, ADDW, AVNC, AVDW, AVKNC, AVKDW
  REAL(8), DIMENSION(:,:),POINTER :: & ! (NRM,4)
       VGR1, VGR2, VGR3, VGR4
  CHARACTER(LEN=40)       :: &
       KGR1, KGR2, KGR3, KGR4

!     ****** GLOBAL VARIABLES ******
! TRGLB
  REAL(8)  ::  &
       WBULKT, WTAILT, WPT, AJT, AJOHT, AJNBT, AJRFT, AJBST, AJTTOR, &
       PINT, POHT, PNBT, PNFT, PBINT, PFINT, POUT, PCXT, PIET, &
       PRBT, PRCT, PRLT, PRSUMT, &
       PEXST, PRFST, SINT, SIET, SNBT, SNFT, SOUT, VLOOP, ALI, RQ1, &
       RPE, Q0, ZEFF0, QF, WPDOT, TAUE1, TAUE2, TAUE89, TAUE98, H98Y2, &
       BETAP0, BETAPA, BETA0, BETAA, BETAQ0, BETAN
  REAL(8), DIMENSION(:)  ,POINTER :: &  ! (NSM)
       SPSCT
  REAL(8), DIMENSION(3)   :: &
       AJRFVT
  REAL(8), DIMENSION(NFM) :: &
       ANF0, TF0, ANFAV, TFAV, WFT
  REAL(8), DIMENSION(:)  ,POINTER :: &  ! (NSTM)
       ANS0, TS0, ANSAV, ANLAV, TSAV, WST, PRFT, PBCLT, PFCLT, PLT, SPET, SLT
  REAL(8), DIMENSION(:,:),POINTER :: &  ! (NSTM,3)
       PRFVT

!     ****** GRAPHIC DATA VARIABLES ******
! TRPLR
  REAL(4), DIMENSION(NTM) :: GT, GTS
  REAL(4), DIMENSION(NGM) :: GTR
  REAL(4), DIMENSION(NLM) :: GBL
  REAL(4), DIMENSION(NTM,8) :: GYT
  REAL(4), DIMENSION(NTM,NCTM) :: GVT
  REAL(4), DIMENSION(NLM,10) :: GBR, GBRH, GBP1, GBAN
  REAL(4), DIMENSION(:)    ,POINTER :: GRM                       ! (NRM)
  REAL(4), DIMENSION(:)    ,POINTER :: GRG                       ! (NRMP)
  REAL(4), DIMENSION(:,:)  ,POINTER :: GJB, GAD                  ! (NRMP,4)
  REAL(4), DIMENSION(:,:)  ,POINTER :: GET                       ! (NRMP,5)
  REAL(4), DIMENSION(:,:)  ,POINTER :: GAK                       ! (NRMP,6)
  REAL(4), DIMENSION(:,:)  ,POINTER :: GYR, GER                  ! (NRMP,8)
  REAL(4), DIMENSION(:,:)  ,POINTER :: GPNB                      ! (4*NRM,10)
  REAL(4), DIMENSION(:,:,:),POINTER :: GVR                       ! (NRMP,NGM,NCGM)
  REAL(4), DIMENSION(:,:,:),POINTER :: GVRT                       ! (NRM,NGM,NCRTM)
  INTEGER(4)              :: NGR, NGT, NGST
  INTEGER(4),DIMENSION(10):: NLMAX

!     ****** EQUILIBRIUM INTERFACE VARIABLES ******
! TREQV1
  REAL(8)                 :: BPSEQ,RIPEQ
  REAL(8), DIMENSION(:),POINTER :: &  ! (NRM)
       RHOTR, PRHO, HJRHO, VTRHO, TRHO, TTRHO, DVRHO, RKPRHO, ABRHO, ABVRHO, &
       ARRHO, AR1RHO, AR2RHO,  RJCB, EPSRHO, BPRHO, RMJRHO, RMNRHO, &
       TTRHOG, DVRHOG, RKPRHOG, ABRHOG, ABVRHOG, ARRHOG, AR1RHOG, AR2RHOG, &
       ABB2RHOG, AIB2RHOG, ARHBRHOG, RMJRHOG, RMNRHOG, RDPVRHOG, &
       PSITRHO, PSIPRHO, PPPRHO, PIQRHO, PIRHO, FACTQ, &
       PVOLRHOG, PSURRHOG, ABB1RHO
  REAL(8), DIMENSION(:),POINTER :: &  ! (NRMP)
       QRHO

!     ****** LOG FILE NAME ******
! TRNAM
  CHARACTER(LEN=80) :: KNAMEQ,KNAMEQ2,KNAMTR,KFNLOG
! TRCOM2
  CHARACTER(LEN=80) :: KXNDEV,KXNDCG,KXNID
  CHARACTER(LEN=80) :: KDIRW1,KDIRW2

!     ****** ADDITIONAL PART BY HONDA MITSURU ******
! TRADD
  REAL(8), PARAMETER          :: VOID = 0.D0
  REAL(8), DIMENSION(3,3)     :: FA, FB, FC
  REAL(8), DIMENSION(:)      ,POINTER :: & ! (NSTM)
       RTM, AMZ, PEXT, PNSSA, PNSA, PTSA
  REAL(8), DIMENSION(:,:)    ,POINTER :: & ! (NRM,NSTM)
       PEX, SEX
  REAL(8), DIMENSION(:,:,:)  ,POINTER :: & ! (NRM,NSM,NSM)
       AKDWD, AKDWP, ADDWD, ADDWP
  REAL(8), DIMENSION(:,:,:,:),POINTER :: & ! (NVM,NVM,4,3)
       VV,DD
  REAL(8), DIMENSION(:,:,:,:),POINTER :: & ! (NVM,NVM,2,3)
       VI,DI
  INTEGER(4)                      :: NTEQIT

!     ****** MODEL SELECTION VARIABLES ******
! TRMDS
  REAL(8)     :: SUMPBM
  INTEGER(4)  :: &
       NEQMAX, MDLEQB, MDLEQN, MDLEQT, MDLEQU, MDLEQZ, MDLEQ0, MDLEQE, &
       MDLEOI, NSCMAX, NSTMAX, MDLWLD, MDDIAG, MDDW, MDLFLX, MDLER, MDCD05, &
       MDEDGE
  INTEGER(4), DIMENSION(:)  ,POINTER :: & ! (NEQM)
       NSS, NSV, NNS, NST
  INTEGER(4), DIMENSION(:,:),POINTER :: & ! (0:NSTM,0:3)
       NEA

!     ****** NCLASS ******
! TRNCL
  REAL(8), DIMENSION(:)    ,POINTER :: & !  (NRM)
       AJBSNC, ETANC, AJEXNC
  REAL(8), DIMENSION(:,:)  ,POINTER :: & !  (NRM,NSM)
       CJBSP, CJBST, ADNCG, AVNCG
  REAL(8), DIMENSION(:,:,:),POINTER :: & !  (NRM,NSM,NSM)
       AKNCP, AKNCT, ADNCP, ADNCT, AKLP, AKLD, ADLP, ADLD
  REAL(8), DIMENSION(:,:,:),POINTER :: & !  (NRM,5,NSM)
       RGFLS,   RQFLS
  INTEGER(4) :: MDNCLS, NSLMAX

!     ****** STORED VARIABLES FOR UFILE ******
! TRSVU
  REAL(8), DIMENSION(NTUM) :: &
       PNFU, PNFUA, RRU, RAU, PHIAU, VOLAU, RIPU, BBU, RKAPU, PNBIU
  REAL(8), DIMENSION(:,:)  ,POINTER :: & ! (NTUM,NRMP)
       QPU,  AJU, AJNBU, BPU, PRLU, PECU, POHU, DVRHOU, AJBSU, ZEFFU, PBMU, &
       RNFU, WROTU, ZEFFU_ORG, RKPRHOU, RMJRHOU, RMNRHOU, ARRHOU, AR1RHOU, &
       AR2RHOU, ABRHOU, TTRHOU, SWLU
  REAL(8), DIMENSION(:,:)  ,POINTER :: & ! (NTUM,NSTM)
       PTSU, PNSU, PTSUA, PNSUA
  REAL(8), DIMENSION(:,:,:),POINTER :: & ! (NTUM,NRMP,NSM)
       RNU, RTU, PNBU, PICU, SNBU, RNU_ORG

!     ****** UFILE CONTROL ******
! TRUFL
  REAL(8)    :: TIME_INT
  INTEGER(4) :: MDLXP,MDLUF, MODEP,MDNI,MDCURT,MDNM1,MDLJQ,MDPHIA,NTS
  CHARACTER(LEN=80) :: KUFDIR,KUFDEV,KUFDCG
!     ****** LAPACK ******
! TRLPCK
  INTEGER(4) :: MDLPCK

!     ************
!TRPROF
  REAL(8), DIMENSION(:),POINTER :: & ! (NSTM)
       PNSSO,PTSO,PNSSAO,PTSAO

!     ******************************************************
  CONTAINS

  SUBROUTINE ALLOCATE_TRCOMM(ierr)

    use trcom1, ONLY : ALLOCATE_TRCOM1
    integer(4), intent(out):: ierr

    ierr = 0
    if(nsmax<1) then
      write(6,*) "XXX ALLOCATE_TRCOMM : ILLEGAL PARAMETER    NSMAX=",nsmax
      ierr = 1
      return
    endif
    if(nrmax<1) then
      write(6,*) "XXX ALLOCATE_TRCOMM : ILLEGAL PARAMETER    NRMAX=",nrmax
      ierr = 1
      return
    endif

    if(nrmax_old==nrmax .and. nsmax_old==nsmax .and. nszmax_old==nszmax .and. nsnmax_old==nsnmax .and. associated(PNSS)) return
    if(associated(PNSS)) call DEALLOCATE_TRCOMM

    NSTMAX  = NSMAX+NSZMAX+NSNMAX
    NEQMAXM = 3*NSTMAX+1
    NVM     = NEQMAXM
    MWM     = 4*NEQMAXM-1
    MLM     = NEQMAXM*NRMAX
    NRMP    = NRMAX+1
    NGLF    = NRMAX
    LDAB    = 6*NEQMAXM


    ALLOCATE(PNSS(NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(XV(NEQMAXM,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(YV(NFM,NRMAX),AY(NFM,NRMAX),Y(NFM,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(ZV(NSMAX,NRMAX),AZ(NSMAX,NRMAX),Z(NSMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AX(6*NEQMAXM,NEQMAXM*NRMAX),X(NEQMAXM*NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(RG(NRMAX),RM(NRMAX),RHOM(NRMAX),RHOG(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(BP(NRMAX),RDP(NRMAX),RPSI(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(RN(NRMAX,NSTM),RT(NRMAX,NSTM),RU(NRMAX,NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(RW(NRMAX,NFM),RNF(NRMAX,NFM),RTF(NRMAX,NFM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(ANC(NRMAX),ANFE(NRMAX),ANNU(NRMAX),ZEFF(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PZC(NRMAX),PZFE(NRMAX),BETA(NRMAX),BETAP(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(BETAL(NRMAX),BETAPL(NRMAX),BETAQ(NRMAX),PBM(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PADD(NRMAX),VTOR(NRMAX),VPAR(NRMAX),VPRP(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(VPOL(NRMAX),WROT(NRMAX),ER(NRMAX),VEXB(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(VEXBP(NRMAX),WEXBP(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(WEXB(NRMAX),AGMP(NRMAX),AJ(NRMAX),AJOH(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(EZOH(NRMAX),QP(NRMAX),AJTOR(NRMAX),AJNB(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(QPINV(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AJRF(NRMAX),AJBS(NRMAX),PNB(NRMAX),SNB(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SPSC(NRMAX,NSM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PBIN(NRMAX),PNF(NRMAX),SNF(NRMAX),PFIN(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(POH(NRMAX),PRB(NRMAX),PRC(NRMAX),PRSUM(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(POH(NRMAX),PRL(NRMAX),PCX(NRMAX),PIE(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SIE(NRMAX),SCX(NRMAX),TSIE(NRMAX),TSCX(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PIN(NRMAX,NSTM),SSIN(NRMAX,NSTM),PBCL(NRMAX,NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SPE(NRMAX,NSTM),PFCL(NRMAX,NSTM),PRF(NRMAX,NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PRFV(NRMAX,NSTM,3),AJRFV(NRMAX,3),RGFLX(NRMAX,NSM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SPSCT(NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(ETA(NRMAX),S(NRMAX),ALPHA(NRMAX),RKCV(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(TAUB(NRMAX),TAUF(NRMAX),TAUK(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AK(NRMAX,NSTM),AVK(NRMAX,NSTM),AD(NRMAX,NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AV(NRMAX,NSTM),AKNC(NRMAX,NSTM),AKDW(NRMAX,NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(ADNC(NRMAX,NSTM),ADDW(NRMAX,NSTM),AVNC(NRMAX,NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AVDW(NRMAX,NSTM),AVKNC(NRMAX,NSTM),AVKDW(NRMAX,NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(VGR1(NRMAX,4),VGR2(NRMAX,4),VGR3(NRMAX,4),VGR4(NRMAX,4),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(ANS0(NSTM),TS0(NSTM),ANSAV(NSTM),ANLAV(NSTM),TSAV(NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(WST(NSTM),PRFT(NSTM),PBCLT(NSTM),PFCLT(NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PLT(NSTM),SPET(NSTM),SLT(NSTM),PRFVT(NSTM,3),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(GRM(NRMAX),GRG(NRMP),GJB(NRMP,4),GAD(NRMP,4),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(GET(NRMP,5),GAK(NRMP,6),GYR(NRMP,8),GER(NRMP,8),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(GPNB(4*NRMAX,10),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(GVR(NRMP,NGM,NCGM),GVRT(NRMAX,NTM,NCRTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(RHOTR(NRMAX),PRHO(NRMAX),HJRHO(NRMAX),VTRHO(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(TRHO(NRMAX),TTRHO(NRMAX),DVRHO(NRMAX),RKPRHO(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(ABRHO(NRMAX),ABVRHO(NRMAX),PSITRHO(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PSIPRHO(NRMAX),PPPRHO(NRMAX),PIQRHO(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PIRHO(NRMAX),ARRHO(NRMAX),AR1RHO(NRMAX),AR2RHO(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(RJCB(NRMAX),EPSRHO(NRMAX),BPRHO(NRMAX),RMJRHO(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(RMNRHO(NRMAX) ,TTRHOG(NRMAX),DVRHOG(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(RKPRHOG(NRMAX),ABRHOG(NRMAX),ABVRHOG(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(ARRHOG(NRMAX),AR1RHOG(NRMAX),AR2RHOG(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(ABB2RHOG(NRMAX),AIB2RHOG(NRMAX),ARHBRHOG(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PVOLRHOG(NRMAX),PSURRHOG(NRMAX),ABB1RHO(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(RMJRHOG(NRMAX),RMNRHOG(NRMAX),RDPVRHOG(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(FACTQ(NRMAX),QRHO(NRMP),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(RTM(NSTM),AMZ(NSTM),PEXT(NSTM),PNSSA(NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PNSA(NSTM),PTSA(NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PEX(NRMAX,NSTM),SEX(NRMAX,NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AKDWD(NRMAX,NSMAX,NSMAX),AKDWP(NRMAX,NSMAX,NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(ADDWD(NRMAX,NSMAX,NSMAX),ADDWP(NRMAX,NSMAX,NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(VV(NEQMAXM,NEQMAXM,4,3),DD(NEQMAXM,NEQMAXM,4,3),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(VI(NEQMAXM,NEQMAXM,2,3),DI(NEQMAXM,NEQMAXM,2,3),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(NSS(NEQMAXM),NSV(NEQMAXM),NNS(NEQMAXM),NST(NEQMAXM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(NEA(0:NSTM,0:3),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(AJBSNC(NRMAX),ETANC(NRMAX),AJEXNC(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(CJBSP(NRMAX,NSM),CJBST(NRMAX,NSM),ADNCG(NRMAX,NSM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AVNCG(NRMAX,NSM) ,STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AKNCP(NRMAX,NSM,NSM),AKNCT(NRMAX,NSM,NSM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(ADNCP(NRMAX,NSM,NSM),ADNCT(NRMAX,NSM,NSM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AKLP(NRMAX,NSM,NSM),AKLD(NRMAX,NSM,NSM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(ADLP(NRMAX,NSM,NSM),ADLD(NRMAX,NSM,NSM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(RGFLS(NRMAX,5,NSM),RQFLS(NRMAX,5,NSM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(QPU(NTUM,NRMP),AJU(NTUM,NRMP),AJNBU(NTUM,NRMP),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(BPU(NTUM,NRMP),PRLU(NTUM,NRMP),PECU(NTUM,NRMP),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(POHU(NTUM,NRMP),DVRHOU(NTUM,NRMP),AJBSU(NTUM,NRMP),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(ZEFFU(NTUM,NRMP),PBMU(NTUM,NRMP),RNFU(NTUM,NRMP),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(WROTU(NTUM,NRMP),ZEFFU_ORG(NTUM,NRMP),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(RKPRHOU(NTUM,NRMP),RMJRHOU(NTUM,NRMP),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(RMNRHOU(NTUM,NRMP),ARRHOU(NTUM,NRMP),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AR1RHOU(NTUM,NRMP),AR2RHOU(NTUM,NRMP),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(ABRHOU(NTUM,NRMP),TTRHOU(NTUM,NRMP),SWLU(NTUM,NRMP),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PTSU(NTUM,NSTM),PNSU(NTUM,NSTM),PTSUA(NTUM,NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PNSUA(NTUM,NSTM) ,STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(RNU(NTUM,NRMP,NSM),RTU(NTUM,NRMP,NSM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PNBU(NTUM,NRMP,NSMAX),PICU(NTUM,NRMP,NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNBU(NTUM,NRMP,NSMAX),RNU_ORG(NTUM,NRMP,NSM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(PNSSO(NSTM),PTSO(NSTM),PNSSAO(NSTM),PTSAO(NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900


    CALL ALLOCATE_TRCOM1(IERR)
      IF(IERR.NE.0) GOTO 900

    nrmax_old  = nrmax
    nsmax_old  = nsmax
    nszmax_old = nszmax
    nsnmax_old = nsnmax
    return

 900 continue
    write(6,*) "XX  TRCOMM ALLOCATION ERROR IERR=",ierr
    call DEALLOCATE_ERR_TRCOMM
    return

  END SUBROUTINE ALLOCATE_TRCOMM


  SUBROUTINE DEALLOCATE_TRCOMM
    use trcom1, ONLY : DEALLOCATE_TRCOM1

    DEALLOCATE(PNSS)
    DEALLOCATE(XV,YV, AY, Y,ZV, AZ, Z,AX,X)
    DEALLOCATE(RG,RM,RHOM,RHOG,BP,RDP,RPSI,RN,RT,RU,RW)
    DEALLOCATE(RNF,RTF,ANC,ANFE,ANNU,ZEFF,PZC,PZFE,BETA,BETAP,BETAL,BETAPL)
    DEALLOCATE(BETAQ,PBM,PADD,VTOR,VPAR,VPRP,VPOL,WROT,ER,VEXB,WEXB,AGMP)
    DEALLOCATE(VEXBP,WEXBP)
    DEALLOCATE(AJ,AJOH, EZOH,QP,AJTOR,AJNB,AJRF,AJBS,QPINV)
    DEALLOCATE(PNB,SNB,PBIN,PNF,SNF,PFIN,POH,PRB,PRC,PRL,PRSUM,PCX,PIE)
    DEALLOCATE(SIE,SCX,TSIE,TSCX,PIN,SSIN,PBCL, SPE)
    DEALLOCATE(PFCL,PRF,PRFV,AJRFV,RGFLX)
    DEALLOCATE(ETA,S,ALPHA,RKCV,TAUB,TAUF,TAUK,AK,AVK,AD,AV)
    DEALLOCATE(AKNC,AKDW,ADNC,ADDW,AVNC,AVDW,AVKNC,AVKDW)
    DEALLOCATE(VGR1,VGR2,VGR3,VGR4)
    DEALLOCATE(ANS0,TS0,ANSAV,ANLAV,TSAV,WST,PRFT,PBCLT,PFCLT)
    DEALLOCATE(PLT,SPET,SLT,PRFVT)
    DEALLOCATE(GRM,GRG,GJB,GAD,GET,GAK,GYR,GER,GPNB,GVR,GVRT)
    DEALLOCATE(RHOTR,PRHO,HJRHO,VTRHO,TRHO,TTRHO,DVRHO,RKPRHO,ABRHO,ABVRHO)
    DEALLOCATE(PSITRHO,PSIPRHO,PPPRHO,PIQRHO,PIRHO)
    DEALLOCATE(ARRHO,AR1RHO,AR2RHO,RJCB,EPSRHO,BPRHO,RMJRHO,RMNRHO)
    DEALLOCATE(TTRHOG,DVRHOG,RKPRHOG,ABRHOG,ABVRHOG,ARRHOG,AR1RHOG,AR2RHOG)
    DEALLOCATE(ABB2RHOG,AIB2RHOG,ARHBRHOG,RMJRHOG,RMNRHOG,RDPVRHOG,FACTQ,QRHO)
    DEALLOCATE(PVOLRHOG,PSURRHOG,ABB1RHO)
    DEALLOCATE(RTM,AMZ,PEXT,PNSSA,PNSA,PTSA,PEX, SEX,AKDWD)
    DEALLOCATE(AKDWP,ADDWD,ADDWP,VV,DD,VI,DI)
    DEALLOCATE(NSS,NSV,NNS,NST,NEA)
    DEALLOCATE(AJBSNC, ETANC, AJEXNC,CJBSP,CJBST,ADNCG,AVNCG)
    DEALLOCATE(AKNCP,AKNCT,ADNCP,ADNCT,AKLP,AKLD,ADLP,ADLD,RGFLS,RQFLS)
    DEALLOCATE(QPU,AJU,AJNBU,BPU,PRLU,PECU,POHU,DVRHOU)
    DEALLOCATE(AJBSU,ZEFFU,PBMU,RNFU,WROTU,ZEFFU_ORG,RKPRHOU)
    DEALLOCATE(RMJRHOU,RMNRHOU,ARRHOU,AR1RHOU,AR2RHOU,ABRHOU)
    DEALLOCATE(TTRHOU,SWLU,PTSU,PNSU,PTSUA,PNSUA)
    DEALLOCATE(RNU,RTU,PNBU,PICU,SNBU,RNU_ORG)
    DEALLOCATE(PNSSO,PTSO,PNSSAO,PTSAO)

    CALL  DEALLOCATE_TRCOM1

    return

  END SUBROUTINE DEALLOCATE_TRCOMM


  SUBROUTINE DEALLOCATE_ERR_TRCOMM
    use trcom1, ONLY : DEALLOCATE_ERR_TRCOM1

    IF(ASSOCIATED(PNSS     ))     DEALLOCATE(PNSS     )
    IF(ASSOCIATED(XV       ))     DEALLOCATE(XV       )
    IF(ASSOCIATED(YV       ))     DEALLOCATE(YV       )
    IF(ASSOCIATED(AY       ))     DEALLOCATE(AY       )
    IF(ASSOCIATED(Y        ))     DEALLOCATE(Y        )
    IF(ASSOCIATED(ZV       ))     DEALLOCATE(ZV       )
    IF(ASSOCIATED(AZ       ))     DEALLOCATE(AZ       )
    IF(ASSOCIATED(Z        ))     DEALLOCATE(Z        )
    IF(ASSOCIATED(AX       ))     DEALLOCATE(AX       )
    IF(ASSOCIATED(X        ))     DEALLOCATE(X        )
    IF(ASSOCIATED(RG       ))     DEALLOCATE(RG       )
    IF(ASSOCIATED(RM       ))     DEALLOCATE(RM       )
    IF(ASSOCIATED(RHOM     ))     DEALLOCATE(RHOM     )
    IF(ASSOCIATED(RHOG     ))     DEALLOCATE(RHOG     )
    IF(ASSOCIATED(BP       ))     DEALLOCATE(BP       )
    IF(ASSOCIATED(RDP      ))     DEALLOCATE(RDP      )
    IF(ASSOCIATED(RPSI     ))     DEALLOCATE(RPSI     )
    IF(ASSOCIATED(RN       ))     DEALLOCATE(RN       )
    IF(ASSOCIATED(RT       ))     DEALLOCATE(RT       )
    IF(ASSOCIATED(RU       ))     DEALLOCATE(RU       )
    IF(ASSOCIATED(RW       ))     DEALLOCATE(RW       )
    IF(ASSOCIATED(RNF      ))     DEALLOCATE(RNF      )
    IF(ASSOCIATED(RTF      ))     DEALLOCATE(RTF      )
    IF(ASSOCIATED(ANC      ))     DEALLOCATE(ANC      )
    IF(ASSOCIATED(ANFE     ))     DEALLOCATE(ANFE     )
    IF(ASSOCIATED(ANNU     ))     DEALLOCATE(ANNU     )
    IF(ASSOCIATED(ZEFF     ))     DEALLOCATE(ZEFF     )
    IF(ASSOCIATED(PZC      ))     DEALLOCATE(PZC      )
    IF(ASSOCIATED(PZFE     ))     DEALLOCATE(PZFE     )
    IF(ASSOCIATED(BETA     ))     DEALLOCATE(BETA     )
    IF(ASSOCIATED(BETAP    ))     DEALLOCATE(BETAP    )
    IF(ASSOCIATED(BETAL    ))     DEALLOCATE(BETAL    )
    IF(ASSOCIATED(BETAPL   ))     DEALLOCATE(BETAPL   )
    IF(ASSOCIATED(BETAQ    ))     DEALLOCATE(BETAQ    )
    IF(ASSOCIATED(PBM      ))     DEALLOCATE(PBM      )
    IF(ASSOCIATED(PADD     ))     DEALLOCATE(PADD     )
    IF(ASSOCIATED(VTOR     ))     DEALLOCATE(VTOR     )
    IF(ASSOCIATED(VPAR     ))     DEALLOCATE(VPAR     )
    IF(ASSOCIATED(VPRP     ))     DEALLOCATE(VPRP     )
    IF(ASSOCIATED(VPOL     ))     DEALLOCATE(VPOL     )
    IF(ASSOCIATED(WROT     ))     DEALLOCATE(WROT     )
    IF(ASSOCIATED(ER       ))     DEALLOCATE(ER       )
    IF(ASSOCIATED(VEXB     ))     DEALLOCATE(VEXB     )
    IF(ASSOCIATED(WEXB     ))     DEALLOCATE(WEXB     )
    IF(ASSOCIATED(AGMP     ))     DEALLOCATE(AGMP     )
    IF(ASSOCIATED(AJ       ))     DEALLOCATE(AJ       )
    IF(ASSOCIATED(AJOH     ))     DEALLOCATE(AJOH     )
    IF(ASSOCIATED(EZOH     ))     DEALLOCATE(EZOH     )
    IF(ASSOCIATED(QP       ))     DEALLOCATE(QP       )
    IF(ASSOCIATED(AJTOR    ))     DEALLOCATE(AJTOR    )
    IF(ASSOCIATED(AJNB     ))     DEALLOCATE(AJNB     )
    IF(ASSOCIATED(AJRF     ))     DEALLOCATE(AJRF     )
    IF(ASSOCIATED(AJBS     ))     DEALLOCATE(AJBS     )
    IF(ASSOCIATED(PNB      ))     DEALLOCATE(PNB      )
    IF(ASSOCIATED(SNB      ))     DEALLOCATE(SNB      )
    IF(ASSOCIATED(PBIN     ))     DEALLOCATE(PBIN     )
    IF(ASSOCIATED(PNF      ))     DEALLOCATE(PNF      )
    IF(ASSOCIATED(SNF      ))     DEALLOCATE(SNF      )
    IF(ASSOCIATED(PFIN     ))     DEALLOCATE(PFIN     )
    IF(ASSOCIATED(POH      ))     DEALLOCATE(POH      )
    IF(ASSOCIATED(PRB      ))     DEALLOCATE(PRB      )
    IF(ASSOCIATED(PRC      ))     DEALLOCATE(PRC      )
    IF(ASSOCIATED(PRL      ))     DEALLOCATE(PRL      )
    IF(ASSOCIATED(PRSUM    ))     DEALLOCATE(PRSUM    )
    IF(ASSOCIATED(PCX      ))     DEALLOCATE(PCX      )
    IF(ASSOCIATED(PIE      ))     DEALLOCATE(PIE      )
    IF(ASSOCIATED(SIE      ))     DEALLOCATE(SIE      )
    IF(ASSOCIATED(SCX      ))     DEALLOCATE(SCX      )
    IF(ASSOCIATED(TSIE     ))     DEALLOCATE(TSIE     )
    IF(ASSOCIATED(TSCX     ))     DEALLOCATE(TSCX     )
    IF(ASSOCIATED(PIN      ))     DEALLOCATE(PIN      )
    IF(ASSOCIATED(SSIN     ))     DEALLOCATE(SSIN     )
    IF(ASSOCIATED(PBCL     ))     DEALLOCATE(PBCL     )
    IF(ASSOCIATED(SPE      ))     DEALLOCATE(SPE      )
    IF(ASSOCIATED(PFCL     ))     DEALLOCATE(PFCL     )
    IF(ASSOCIATED(PRF      ))     DEALLOCATE(PRF      )
    IF(ASSOCIATED(PRFV     ))     DEALLOCATE(PRFV     )
    IF(ASSOCIATED(AJRFV    ))     DEALLOCATE(AJRFV    )
    IF(ASSOCIATED(RGFLX    ))     DEALLOCATE(RGFLX    )
    IF(ASSOCIATED(ETA      ))     DEALLOCATE(ETA      )
    IF(ASSOCIATED(S        ))     DEALLOCATE(S        )
    IF(ASSOCIATED(ALPHA    ))     DEALLOCATE(ALPHA    )
    IF(ASSOCIATED(RKCV     ))     DEALLOCATE(RKCV     )
    IF(ASSOCIATED(TAUB     ))     DEALLOCATE(TAUB     )
    IF(ASSOCIATED(TAUF     ))     DEALLOCATE(TAUF     )
    IF(ASSOCIATED(TAUK     ))     DEALLOCATE(TAUK     )
    IF(ASSOCIATED(AK       ))     DEALLOCATE(AK       )
    IF(ASSOCIATED(AVK      ))     DEALLOCATE(AVK      )
    IF(ASSOCIATED(AD       ))     DEALLOCATE(AD       )
    IF(ASSOCIATED(AV       ))     DEALLOCATE(AV       )
    IF(ASSOCIATED(AKNC     ))     DEALLOCATE(AKNC     )
    IF(ASSOCIATED(AKDW     ))     DEALLOCATE(AKDW     )
    IF(ASSOCIATED(ADNC     ))     DEALLOCATE(ADNC     )
    IF(ASSOCIATED(ADDW     ))     DEALLOCATE(ADDW     )
    IF(ASSOCIATED(AVNC     ))     DEALLOCATE(AVNC     )
    IF(ASSOCIATED(AVDW     ))     DEALLOCATE(AVDW     )
    IF(ASSOCIATED(AVKNC    ))     DEALLOCATE(AVKNC    )
    IF(ASSOCIATED(AVKDW    ))     DEALLOCATE(AVKDW    )
    IF(ASSOCIATED(VGR1     ))     DEALLOCATE(VGR1     )
    IF(ASSOCIATED(VGR2     ))     DEALLOCATE(VGR2     )
    IF(ASSOCIATED(VGR3     ))     DEALLOCATE(VGR3     )
    IF(ASSOCIATED(VGR4     ))     DEALLOCATE(VGR4     )
    IF(ASSOCIATED(ANS0     ))     DEALLOCATE(ANS0     )
    IF(ASSOCIATED(TS0      ))     DEALLOCATE(TS0      )
    IF(ASSOCIATED(ANSAV    ))     DEALLOCATE(ANSAV    )
    IF(ASSOCIATED(ANLAV    ))     DEALLOCATE(ANLAV    )
    IF(ASSOCIATED(TSAV     ))     DEALLOCATE(TSAV     )
    IF(ASSOCIATED(WST      ))     DEALLOCATE(WST      )
    IF(ASSOCIATED(PRFT     ))     DEALLOCATE(PRFT     )
    IF(ASSOCIATED(PBCLT    ))     DEALLOCATE(PBCLT    )
    IF(ASSOCIATED(PFCLT    ))     DEALLOCATE(PFCLT    )
    IF(ASSOCIATED(PLT      ))     DEALLOCATE(PLT      )
    IF(ASSOCIATED(SPET     ))     DEALLOCATE(SPET     )
    IF(ASSOCIATED(SLT      ))     DEALLOCATE(SLT      )
    IF(ASSOCIATED(PRFVT    ))     DEALLOCATE(PRFVT    )
    IF(ASSOCIATED(GRM      ))     DEALLOCATE(GRM      )
    IF(ASSOCIATED(GRG      ))     DEALLOCATE(GRG      )
    IF(ASSOCIATED(GJB      ))     DEALLOCATE(GJB      )
    IF(ASSOCIATED(GAD      ))     DEALLOCATE(GAD      )
    IF(ASSOCIATED(GET      ))     DEALLOCATE(GET      )
    IF(ASSOCIATED(GAK      ))     DEALLOCATE(GAK      )
    IF(ASSOCIATED(GYR      ))     DEALLOCATE(GYR      )
    IF(ASSOCIATED(GER      ))     DEALLOCATE(GER      )
    IF(ASSOCIATED(GPNB     ))     DEALLOCATE(GPNB     )
    IF(ASSOCIATED(GVR      ))     DEALLOCATE(GVR      )
    IF(ASSOCIATED(GVRT     ))     DEALLOCATE(GVRT     )
    IF(ASSOCIATED(RHOTR    ))     DEALLOCATE(RHOTR    )
    IF(ASSOCIATED(PRHO     ))     DEALLOCATE(PRHO     )
    IF(ASSOCIATED(HJRHO    ))     DEALLOCATE(HJRHO    )
    IF(ASSOCIATED(VTRHO    ))     DEALLOCATE(VTRHO    )
    IF(ASSOCIATED(TRHO     ))     DEALLOCATE(TRHO     )
    IF(ASSOCIATED(TTRHO    ))     DEALLOCATE(TTRHO    )
    IF(ASSOCIATED(PSITRHO  ))     DEALLOCATE(PSITRHO  )
    IF(ASSOCIATED(PSIPRHO  ))     DEALLOCATE(PSIPRHO  )
    IF(ASSOCIATED(PPPRHO   ))     DEALLOCATE(PPPRHO   )
    IF(ASSOCIATED(PIQRHO   ))     DEALLOCATE(PIQRHO   )
    IF(ASSOCIATED(PIRHO    ))     DEALLOCATE(PIRHO    )
    IF(ASSOCIATED(DVRHO    ))     DEALLOCATE(DVRHO    )
    IF(ASSOCIATED(RKPRHO   ))     DEALLOCATE(RKPRHO   )
    IF(ASSOCIATED(ABRHO    ))     DEALLOCATE(ABRHO    )
    IF(ASSOCIATED(ABVRHO   ))     DEALLOCATE(ABVRHO   )
    IF(ASSOCIATED(ARRHO    ))     DEALLOCATE(ARRHO    )
    IF(ASSOCIATED(AR1RHO   ))     DEALLOCATE(AR1RHO   )
    IF(ASSOCIATED(AR2RHO   ))     DEALLOCATE(AR2RHO   )
    IF(ASSOCIATED(RJCB     ))     DEALLOCATE(RJCB     )
    IF(ASSOCIATED(EPSRHO   ))     DEALLOCATE(EPSRHO   )
    IF(ASSOCIATED(BPRHO    ))     DEALLOCATE(BPRHO    )
    IF(ASSOCIATED(RMJRHO   ))     DEALLOCATE(RMJRHO   )
    IF(ASSOCIATED(RMNRHO   ))     DEALLOCATE(RMNRHO   )
    IF(ASSOCIATED(TTRHOG   ))     DEALLOCATE(TTRHOG   )
    IF(ASSOCIATED(DVRHOG   ))     DEALLOCATE(DVRHOG   )
    IF(ASSOCIATED(RKPRHOG  ))     DEALLOCATE(RKPRHOG  )
    IF(ASSOCIATED(ABRHOG   ))     DEALLOCATE(ABRHOG   )
    IF(ASSOCIATED(ABVRHOG  ))     DEALLOCATE(ABVRHOG  )
    IF(ASSOCIATED(ARRHOG   ))     DEALLOCATE(ARRHOG   )
    IF(ASSOCIATED(AR1RHOG  ))     DEALLOCATE(AR1RHOG  )
    IF(ASSOCIATED(AR2RHOG  ))     DEALLOCATE(AR2RHOG  )
    IF(ASSOCIATED(ABB2RHOG ))     DEALLOCATE(ABB2RHOG )
    IF(ASSOCIATED(AIB2RHOG ))     DEALLOCATE(AIB2RHOG )
    IF(ASSOCIATED(ARHBRHOG ))     DEALLOCATE(ARHBRHOG )
    IF(ASSOCIATED(RMJRHOG  ))     DEALLOCATE(RMJRHOG  )
    IF(ASSOCIATED(RMNRHOG  ))     DEALLOCATE(RMNRHOG  )
    IF(ASSOCIATED(PVOLRHOG ))     DEALLOCATE(PVOLRHOG )
    IF(ASSOCIATED(PSURRHOG ))     DEALLOCATE(PSURRHOG )
    IF(ASSOCIATED(ABB1RHO  ))     DEALLOCATE(ABB1RHO  )
    IF(ASSOCIATED(RDPVRHOG ))     DEALLOCATE(RDPVRHOG )
    IF(ASSOCIATED(FACTQ    ))     DEALLOCATE(FACTQ    )
    IF(ASSOCIATED(QRHO     ))     DEALLOCATE(QRHO     )
    IF(ASSOCIATED(RTM      ))     DEALLOCATE(RTM      )
    IF(ASSOCIATED(AMZ      ))     DEALLOCATE(AMZ      )
    IF(ASSOCIATED(PEXT     ))     DEALLOCATE(PEXT     )
    IF(ASSOCIATED(PNSSA    ))     DEALLOCATE(PNSSA    )
    IF(ASSOCIATED(PNSA     ))     DEALLOCATE(PNSA     )
    IF(ASSOCIATED(PTSA     ))     DEALLOCATE(PTSA     )
    IF(ASSOCIATED(PEX      ))     DEALLOCATE(PEX      )
    IF(ASSOCIATED(SEX      ))     DEALLOCATE(SEX      )
    IF(ASSOCIATED(AKDWD    ))     DEALLOCATE(AKDWD    )
    IF(ASSOCIATED(AKDWP    ))     DEALLOCATE(AKDWP    )
    IF(ASSOCIATED(ADDWD    ))     DEALLOCATE(ADDWD    )
    IF(ASSOCIATED(ADDWP    ))     DEALLOCATE(ADDWP    )
    IF(ASSOCIATED(VV       ))     DEALLOCATE(VV       )
    IF(ASSOCIATED(DD       ))     DEALLOCATE(DD       )
    IF(ASSOCIATED(VI       ))     DEALLOCATE(VI       )
    IF(ASSOCIATED(DI       ))     DEALLOCATE(DI       )
    IF(ASSOCIATED(NSS      ))     DEALLOCATE(NSS      )
    IF(ASSOCIATED(NSV      ))     DEALLOCATE(NSV      )
    IF(ASSOCIATED(NNS      ))     DEALLOCATE(NNS      )
    IF(ASSOCIATED(NST      ))     DEALLOCATE(NST      )
    IF(ASSOCIATED(NEA      ))     DEALLOCATE(NEA      )
    IF(ASSOCIATED(AJBSNC   ))     DEALLOCATE(AJBSNC   )
    IF(ASSOCIATED(ETANC    ))     DEALLOCATE(ETANC    )
    IF(ASSOCIATED(AJEXNC   ))     DEALLOCATE(AJEXNC   )
    IF(ASSOCIATED(CJBSP    ))     DEALLOCATE(CJBSP    )
    IF(ASSOCIATED(CJBST    ))     DEALLOCATE(CJBST    )
    IF(ASSOCIATED(ADNCG    ))     DEALLOCATE(ADNCG    )
    IF(ASSOCIATED(AVNCG    ))     DEALLOCATE(AVNCG    )
    IF(ASSOCIATED(AKNCP    ))     DEALLOCATE(AKNCP    )
    IF(ASSOCIATED(AKNCT    ))     DEALLOCATE(AKNCT    )
    IF(ASSOCIATED(ADNCP    ))     DEALLOCATE(ADNCP    )
    IF(ASSOCIATED(ADNCT    ))     DEALLOCATE(ADNCT    )
    IF(ASSOCIATED(AKLP     ))     DEALLOCATE(AKLP     )
    IF(ASSOCIATED(AKLD     ))     DEALLOCATE(AKLD     )
    IF(ASSOCIATED(ADLP     ))     DEALLOCATE(ADLP     )
    IF(ASSOCIATED(ADLD     ))     DEALLOCATE(ADLD     )
    IF(ASSOCIATED(RGFLS    ))     DEALLOCATE(RGFLS    )
    IF(ASSOCIATED(RQFLS    ))     DEALLOCATE(RQFLS    )
    IF(ASSOCIATED(QPU      ))     DEALLOCATE(QPU      )
    IF(ASSOCIATED(AJU      ))     DEALLOCATE(AJU      )
    IF(ASSOCIATED(AJNBU    ))     DEALLOCATE(AJNBU    )
    IF(ASSOCIATED(BPU      ))     DEALLOCATE(BPU      )
    IF(ASSOCIATED(PRLU     ))     DEALLOCATE(PRLU     )
    IF(ASSOCIATED(PECU     ))     DEALLOCATE(PECU     )
    IF(ASSOCIATED(POHU     ))     DEALLOCATE(POHU     )
    IF(ASSOCIATED(DVRHOU   ))     DEALLOCATE(DVRHOU   )
    IF(ASSOCIATED(AJBSU    ))     DEALLOCATE(AJBSU    )
    IF(ASSOCIATED(ZEFFU    ))     DEALLOCATE(ZEFFU    )
    IF(ASSOCIATED(PBMU     ))     DEALLOCATE(PBMU     )
    IF(ASSOCIATED(RNFU     ))     DEALLOCATE(RNFU     )
    IF(ASSOCIATED(WROTU    ))     DEALLOCATE(WROTU    )
    IF(ASSOCIATED(ZEFFU_ORG))     DEALLOCATE(ZEFFU_ORG)
    IF(ASSOCIATED(RKPRHOU  ))     DEALLOCATE(RKPRHOU  )
    IF(ASSOCIATED(RMJRHOU  ))     DEALLOCATE(RMJRHOU  )
    IF(ASSOCIATED(RMNRHOU  ))     DEALLOCATE(RMNRHOU  )
    IF(ASSOCIATED(ARRHOU   ))     DEALLOCATE(ARRHOU   )
    IF(ASSOCIATED(AR1RHOU  ))     DEALLOCATE(AR1RHOU  )
    IF(ASSOCIATED(AR2RHOU  ))     DEALLOCATE(AR2RHOU  )
    IF(ASSOCIATED(ABRHOU   ))     DEALLOCATE(ABRHOU   )
    IF(ASSOCIATED(TTRHOU   ))     DEALLOCATE(TTRHOU   )
    IF(ASSOCIATED(SWLU     ))     DEALLOCATE(SWLU     )
    IF(ASSOCIATED(PTSU     ))     DEALLOCATE(PTSU     )
    IF(ASSOCIATED(PNSU     ))     DEALLOCATE(PNSU     )
    IF(ASSOCIATED(PTSUA    ))     DEALLOCATE(PTSUA    )
    IF(ASSOCIATED(PNSUA    ))     DEALLOCATE(PNSUA    )
    IF(ASSOCIATED(RNU      ))     DEALLOCATE(RNU      )
    IF(ASSOCIATED(RTU      ))     DEALLOCATE(RTU      )
    IF(ASSOCIATED(PNBU     ))     DEALLOCATE(PNBU     )
    IF(ASSOCIATED(PICU     ))     DEALLOCATE(PICU     )
    IF(ASSOCIATED(SNBU     ))     DEALLOCATE(SNBU     )
    IF(ASSOCIATED(RNU_ORG  ))     DEALLOCATE(RNU_ORG  )
    IF(ASSOCIATED(PNSSO    ))     DEALLOCATE(PNSSO    )
    IF(ASSOCIATED(PTSO     ))     DEALLOCATE(PTSO     )
    IF(ASSOCIATED(PNSSAO   ))     DEALLOCATE(PNSSAO   )
    IF(ASSOCIATED(PTSAO    ))     DEALLOCATE(PTSAO    )

    CALL DEALLOCATE_ERR_TRCOM1

    return

  END SUBROUTINE DEALLOCATE_ERR_TRCOMM

  SUBROUTINE open_trcomm
    RETURN
  END SUBROUTINE open_trcomm

END MODULE TRCOMM
