! trcomm.f90

MODULE trcomm_parm

  USE plcomm,pm=>pa
  USE commpi
  IMPLICIT NONE

! IMPORTED FROM plcomm
!     RR,RA,RKAP,RDLT,BB,RIP,
!     NSM,NSMAX,
!     NPA(NSM),PM(NSM),PZ(NSM),PN(NSM),PNS(NSM),PTPR(NSM),PTPP(NSM),PTS(NSM),
!     PU(NSM),PUS(NSM),ID_NS(NSM),KID_NS(NSM)
!     PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2
!     MODELG,MODELN,MODELQ,MODEL_NPROF
!     KNAMEQ,KNAMEQ2,KNAMTR
!     MODEFR,MODEFW,IDEBUG
!     RHOA,Q0,QA (not an input parameter for tr)

  INTEGER, PARAMETER :: NSAM=6     ! Max number of bulk species
  INTEGER, PARAMETER :: NSFM=6     ! Max number of fast ion species
  INTEGER, PARAMETER :: NSZM=2     ! Max number of impurity ions species
  INTEGER, PARAMETER :: NSNM=2     ! Max number of neutral species

  INTEGER, PARAMETER :: NTM=10001  ! Max number of time steps to save globals
  INTEGER, PARAMETER :: NGM=1001   ! Max number of time steps to save profiles
  INTEGER, PARAMETER :: NCTM=110   ! Max number of save variables in globals
  INTEGER, PARAMETER :: NCGM=30    ! Max number of save variables in profiles
  INTEGER, PARAMETER :: NCRTM=75   ! Max number of save variables in profiles
  
  INTEGER, PARAMETER :: NSTM=NSM+NSZM+NSNM ! number of thermal species
  INTEGER, PARAMETER :: NEQM=3*NSTM+1      ! number of equations
  INTEGER, PARAMETER :: NLM =1001  ! Max number of NBI path sections
  INTEGER, PARAMETER :: NPSCM = 4  ! Maximum number of particle source
  INTEGER, PARAMETER :: NNBM = 4   ! Maximum number of NBI source
  INTEGER, PARAMETER :: NECM = 4   ! Maximum number of ECRF source
  INTEGER, PARAMETER :: NLHM = 4   ! Maximum number of LHRF source
  INTEGER, PARAMETER :: NICM = 4   ! Maximum number of ICRF source
  INTEGER, PARAMETER :: NNFM = 4   ! Maximum number of fusion source
  INTEGER, PARAMETER :: NPELM = 4  ! Maximum number of PELLET source

  REAL(rkind), PARAMETER :: RKEV=AEE*1.D3

  ! ****** INPUT PARAMETERS ******
  
  ! === configuration parameters ===

  INTEGER:: NTEQIT

  ! === plasma parameters ===

  REAL(rkind),DIMENSION(NSM):: PT
  REAL(rkind):: RIPS,RIPE
  INTEGER:: NSAMAX,NSFMAX,NSZMAX,NSNMAX
  INTEGER:: &
       NS_e,NS_D,NS_T,NS_A,NS_H,NS_He3,NS_C,NS_Fe
  
  ! === profile parameters ===

  CHARACTER(LEN=128):: knam_prof
  REAL(rkind):: &
       PROFNU1,PROFNU2,PROFJ1,PROFJ2
  REAL(rkind):: ALP(6)
  INTEGER:: model_nfixed,model_tfixed
  CHARACTER(LEN=128):: knam_nfixed,knam_tfixed

  ! === impurity and neutral parameters ===

  INTEGER:: MDLIMP,MDLNI
  REAL(rkind):: PNC,PNFE,PNNU,PNNUS
  
  ! === transport model parameters ===

  INTEGER:: MDLKAI
  INTEGER:: MDLETA,MDLAD,MDLAVK,MDLJBS,MDLKNC,MDLTPF

  ! === NCLASS switch ===
  
  INTEGER:: MDLNCL,NSLMAX

  ! === Turbulence model switch ===
  
  INTEGER:: MDLWLD,MDLCD05,MDLDW
  INTEGER:: MDLTC

  ! === Transport parameters ===
  
  REAL(rkind):: AD0,AV0,CNP,CNH,CDP,CDH,CNN,CDW(8)
  REAL(rkind):: CHP,CK0,CK1,CWEB,CALF,CKALFA,CKBETA,CKGUMA

  ! ==- radial electric field model parameter ===

  INTEGER:: MDLER

  ! === Edge model parameter ===

  INTEGER:: MDLEDGE,NREDGE
  REAL(rkind):: CSPRS,RHOA

  ! == ELM MODEL parameters ===

  INTEGER:: MDLElM
  REAL(rkind)   :: ELMWID,ELMDUR  
  REAL(rkind),DIMENSION(NSM) :: ELMNRD,ELMTRD,ELMENH

  ! === radiation parameter ===

  INTEGER:: MDLPR
  REAL(rkind):: SYNCABS,SYNCSELF

  ! === Source parameters: NB,EC,LH,IC,PEL,PSC ===

  INTEGER:: NNBMAX
  INTEGER,DIMENSION(nnbm):: &
       model_nnb,ns_nnb,nrmax_nnb
  REAL(rkind),DIMENSION(nnbm):: &
       PNBIN,PNBR0,PNBRW,PNBCD,PNBVY,PNBVW,PNBENG,PNBRTG

  INTEGER:: NECMAX,MDLEC(necm)
  REAL(rkind),DIMENSION(necm):: &
       PECIN,PECR0,PECRW,PECCD,PECTOE,PECNPR
  
  INTEGER:: NLHMAX,MDLLH(nlhm)
  REAL(rkind),DIMENSION(nlhm):: &
       PLHIN,PLHR0,PLHRW,PLHCD,PLHTOE,PLHNPR
  
  INTEGER:: NICMAX,MDLIC(nicm)
  REAL(rkind),DIMENSION(nicm):: &
       PICIN,PICR0,PICRW,PICCD,PICTOE,PICNPR

  INTEGER:: NPELMAX,MDLPEL(npelm),number_of_pellet_repeat(npelm)
  REAL(rkind),DIMENSION(npelm):: &
       PELIN,PELR0,PELRW,PELRAD,PELVEL,PELTIM, &
       pellet_time_start,pellet_time_interval
  REAL(rkind),DIMENSION(nsm,npelm):: &
       PELPAT
!  INTEGER:: NPELMAX,MDLPEL(npelm),number_of_pellet_repeat(npelm)
!  REAL(rkind),DIMENSION(npelm):: &
!       PELTOT,PELR0,PELRW,PELRAD,PELVEL,PELTIM, &
!       pellet_time_start,pellet_time_interval
!  REAL(rkind),DIMENSION(nsmm,npelm):: &
!       PELPAT
  
  INTEGER:: NPSCMAX,MDLPSC(npscm),NSPSC(npscm)
  REAL(rkind), DIMENSION(npscm) :: &
       PSCIN,PSCR0,PSCRW

  INTEGER:: NNFMAX,model_nnf(nnfm),ns_nnf(nnfm)

  ! === current drive parameters ===

  INTEGER:: MDLCD
  REAL(rkind):: PBSCD

  ! === sawtooth parameters ===

  INTEGER:: MDLST,IZERO,NGPST
  REAL(rkind):: TPRST,TSST

  ! === experimental data parameters ===

  INTEGER:: MODEP

  ! === current profile switch ===

  INTEGER:: MDLJQ

  ! === flux switch ===

  INTEGER:: MDLFLX

  ! === simulation parameter ===

  REAL(rkind):: DT
  INTEGER:: NRMAX,NTMAX,NTSTEP,NGTSTP,NGRSTP
  REAL(rkind):: EPSLTR
  INTEGER:: LMAXTR

  ! === equation selection parameter ===

  INTEGER:: MDLEQB,MDLEQN,MDLEQT,MDLEQU,MDLEQZ,MDLEQ0,MDLEQE,MDLEOI

  ! === LAPACK parameter ===

  INTEGER:: MDLPCK

  ! === DATA file name ===

  CHARACTER(LEN=80) :: KFNLOG,KFNTXT,KFNCVS
  
END MODULE trcomm_parm

MODULE trcomx ! common for tr solver core
  USE trcomm_parm,ONLY: rkind
  IMPLICIT NONE
  REAL(rkind), DIMENSION(:,:,:),ALLOCATABLE :: A, B, C
  REAL(rkind), DIMENSION(:,:)  ,ALLOCATABLE :: D
  REAL(rkind), DIMENSION(:,:)  ,ALLOCATABLE :: RD, PPA, PPB, PPC

CONTAINS

  SUBROUTINE allocate_trcomx(NEQMAXM,NRMAX)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NEQMAXM,NRMAX
    INTEGER:: IERR

    IF(ALLOCATED(A)) CALL deallocate_trcomx
    IERR=0
    ALLOCATE(A(NEQMAXM,NEQMAXM,NRMAX),B(NEQMAXM,NEQMAXM,NRMAX),STAT=IERR)
    IF(IERR.NE.0) GOTO 900
    ALLOCATE(C(NEQMAXM,NEQMAXM,NRMAX),D(NEQMAXM,NRMAX),STAT=IERR)
    IF(IERR.NE.0) GOTO 900
    ALLOCATE(RD(NEQMAXM,NRMAX),STAT=IERR)
    IF(IERR.NE.0) GOTO 900
    ALLOCATE(PPA(NEQMAXM,NRMAX),PPB(NEQMAXM,NRMAX),STAT=IERR)
    IF(IERR.NE.0) GOTO 900
    ALLOCATE(PPC(NEQMAXM,NRMAX),STAT=IERR)
    IF(IERR.NE.0) GOTO 900
    RETURN

900 WRITE(6,'(A,I8)') 'XX allocate_trcomx: ierr=',IERR
    IF(ALLOCATED(A)) DEALLOCATE(A)
    IF(ALLOCATED(B)) DEALLOCATE(B)
    IF(ALLOCATED(C)) DEALLOCATE(C)
    IF(ALLOCATED(D)) DEALLOCATE(D)
    IF(ALLOCATED(RD)) DEALLOCATE(RD)
    IF(ALLOCATED(PPA)) DEALLOCATE(PPA)
    IF(ALLOCATED(PPB)) DEALLOCATE(PPB)
    IF(ALLOCATED(PPC)) DEALLOCATE(PPC)
    STOP
  END SUBROUTINE allocate_trcomx

  SUBROUTINE deallocate_trcomx

    IMPLICIT NONE

    DEALLOCATE(A,B,C,D,RD,PPA,PPB,PPC)
    RETURN
  END SUBROUTINE deallocate_trcomx
END MODULE trcomx

MODULE trcomm
  
  USE trcomm_parm
  IMPLICIT NONE

!     ****** CONTROL VARIABLES ******
!TRCTL
  REAL   :: &
       GTCPU1
  REAL(rkind)   :: &
       T, TST, TPRE, WPPRE, DR, FKAP, RIPA, VSEC, RDPS,DIPDT
  REAL(rkind), DIMENSION(:),ALLOCATABLE :: & ! (NSTM)
       PNSS
  INTEGER:: &
       NT, NRAMAX, NROMAX, NTMAX_SAVE, IREAD
  INTEGER:: &
       NFMAX      ! Max number of fast ion species (NNBMAX+NNFMAX)
  INTEGER:: &
       NEQMAXM, NVM, MWM, MLM, NRMP, NGLF, LDAB
  INTEGER:: &
       icount_of_pellet(npelm) ! 0 : before start
                        ! positive: number of pellet from t_start

!     ****** MATRIX VARIBALES ******
! TRMTX
  REAL(rkind), DIMENSION(:,:),ALLOCATABLE :: & ! (NVM,NRM)
       XV
  REAL(rkind), DIMENSION(:,:),ALLOCATABLE :: & ! (NFM,NRM)
       YV, AY, Y
  REAL(rkind), DIMENSION(:,:),ALLOCATABLE :: & ! (NSM,NRM)
       ZV, AZ, Z
  REAL(rkind), DIMENSION(:,:),ALLOCATABLE :: & ! (LDAB,MLM)
       AX
  REAL(rkind), DIMENSION(:)  ,ALLOCATABLE :: & ! (MLM)
       X
!
!     ****** FUNDAMENTAL VARIABLES ******
! TRVAR
  REAL(rkind), DIMENSION(:)  ,ALLOCATABLE :: & ! (NRM)
       RG, RM, RHOM, RHOG, BP, RDP, RPSI
  REAL(rkind), DIMENSION(:,:),ALLOCATABLE :: & ! (NRM,NSTM)
       RN, RT, RU
  REAL(rkind), DIMENSION(:,:),ALLOCATABLE :: & ! (NRM,NFM)
       RW

!     ****** PROFILE VARIABLES ******
! TRPFL
  REAL(rkind), DIMENSION(:,:), ALLOCATABLE :: & ! (NRM,NFM)
       RNF, RTF
  REAL(rkind), DIMENSION(:)  , ALLOCATABLE :: & ! (NRM)
       ANC, ANFE, ANNU, ZEFF, PZC, PZFE, BETA, BETAP, BETAL, BETAPL, &
       BETAQ, PBM, PADD, VTOR, VPAR, VPRP, VPOL, WROT, ER, VEXB, WEXB, AGMP, &
       VEXBP, WEXBP

!     ****** SOURCE VARIABLES ******
! TRSRC
  REAL(rkind), DIMENSION(10)      :: &
       RTG
  REAL(rkind), DIMENSION(:),     ALLOCATABLE :: & ! (NRM)
       AJ, AJOH, EZOH, QP, AJTOR, AJNB, AJRF, AJBS, QPINV, &
       POH, PRB, PRC, PRL, PRSUM, &
       PCX, PIE, SIE, SCX, TSIE, TSCX
  REAL(rkind), DIMENSION(:,:),   ALLOCATABLE :: & ! (NSCM,NRM)
       AJNB_NNBNR,PEC_NEC,PLH_NLH,PIC_NIC
  REAL(rkind), DIMENSION(:,:),   ALLOCATABLE :: & ! (NSCM,NRM)
       SNB_NNBNR,PNB_NNBNR
  REAL(rkind), DIMENSION(:,:,:),     ALLOCATABLE :: & ! (NSM,NRM)
       SNB_NSNNBNR,PNB_NSNNBNR
  REAL(rkind), DIMENSION(:,:),     ALLOCATABLE :: & ! (NSM,NRM)
       SNB_NSNR,PNB_NSNR
  REAL(rkind), DIMENSION(:),     ALLOCATABLE :: & ! (NRM)
       SNB_NR,PNB_NR
  REAL(rkind), DIMENSION(:),     ALLOCATABLE :: & ! (NRM)
       SNB_NS,PNB_NS
  REAL(rkind), DIMENSION(:),     ALLOCATABLE :: & ! (NNBM)
       SNB_NNB,PNB_NNB
  REAL(rkind), DIMENSION(:,:),   ALLOCATABLE :: & ! (NSCM,NRM)
       SNF_NNFNR,PNF_NNFNR
  REAL(rkind), DIMENSION(:,:,:),   ALLOCATABLE :: & ! (NSM,NRM)
       SNF_NSNNFNR,PNF_NSNNFNR
  REAL(rkind), DIMENSION(:,:),     ALLOCATABLE :: & ! (NSM,NRM)
       SNF_NSNR,PNF_NSNR
  REAL(rkind), DIMENSION(:),     ALLOCATABLE :: & ! (NRM)
       SNF_NR,PNF_NR
  REAL(rkind), DIMENSION(:),     ALLOCATABLE :: & ! (NRM)
       SNF_NS,PNF_NS
  REAL(rkind), DIMENSION(:),     ALLOCATABLE :: & ! (NNFM)
       SNF_NNF,PNF_NNF,ENF_NNF
  REAL(rkind), DIMENSION(:,:),     ALLOCATABLE :: & ! (NXM,NRM)
       PNBIN_NNBNR,PNFIN_NNFNR
  REAL(rkind), DIMENSION(:,:),     ALLOCATABLE :: & ! (NSM,NXM,NRM)
       PNBCL_NSNR,PNFCL_NSNR
  REAL(rkind), DIMENSION(:,:,:),     ALLOCATABLE :: & ! (NSM,NXM,NRM)
       PNBCL_NSNNBNR,PNFCL_NSNNFNR
  REAL(rkind), DIMENSION(:,:),   ALLOCATABLE :: & ! (NRM,NSTM)
       PIN, SSIN, PBCL, PFCL, PRF, SPE
  REAL(rkind), DIMENSION(:,:,:),   ALLOCATABLE :: & ! (NSM,NPELM,NRM))
       SPE_NSNPELNR
  REAL(rkind), DIMENSION(:,:),   ALLOCATABLE :: & ! (NRM,NSM)
       RGFLX,SPSC
  REAL(rkind), DIMENSION(:,:),   ALLOCATABLE :: & ! (NRM,3)
       AJRFV
  REAL(rkind), DIMENSION(:,:,:), ALLOCATABLE :: & ! (NRM,NSTM,3)
       PRFV
  REAL(rkind):: &
       PECTOT,PLHTOT,PICTOT,PNBTOT,PELTOT,PSCTOT

!     ****** COEFFICIENT VARIABLES ******
! TRCEF
  REAL(rkind), DIMENSION(:)  ,ALLOCATABLE :: & ! (NRM)
       ETA, S, ALPHA, RKCV, TAUK
  REAL(rkind), DIMENSION(:,:)  ,ALLOCATABLE :: & ! (NRM)
       TAUB, TAUF
  REAL(rkind), DIMENSION(:,:),ALLOCATABLE :: & ! (NRM,NSTM)
       AK, AVK, AD, AV, AKNC, AKDW, ADNC, ADDW, AVNC, AVDW, AVKNC, AVKDW
  REAL(rkind), DIMENSION(:,:),ALLOCATABLE :: & ! (NRM,4)
       VGR1, VGR2, VGR3, VGR4
  CHARACTER(LEN=46)       :: &
       KGR1, KGR2, KGR3, KGR4

!     ****** GLOBAL VARIABLES ******
! TRGLB
  REAL(rkind)  ::  &
       WBULKT, WTAILT, WPT, AJT, AJOHT, AJNBT, AJRFT, AJBST, AJTTOR, &
       PINT, POHT, PNBT, PNFT, PNBINT, PNFINT, POUT, PCXT, PIET, &
       PRBT, PRCT, PRLT, PRSUMT, &
       PEXST, PRFST, SINT, SIET, SNBT, SNFT, SOUT, VLOOP, ALI, RQ1, &
       RPE, ZEFF0, QF, WPDOT, TAUE1, TAUE2, TAUE89, TAUE98, H98Y2, &
       BETAP0, BETAPA, BETA0, BETAA, BETAQ0, BETAN
  REAL(rkind), DIMENSION(:), ALLOCATABLE :: &  ! (NSM)
       SPSCT,PNBCLT,PNFCLT
  REAL(rkind), DIMENSION(3)   :: &
       AJRFVT
  REAL(rkind), DIMENSION(:), ALLOCATABLE :: &
       ANF0, TF0, ANFAV, TFAV, WFT ! (NFM)
  REAL(rkind), DIMENSION(:)  ,ALLOCATABLE :: &  ! (NSTM)
       ANS0, TS0, ANSAV, ANLAV, TSAV, WST, PRFT, PBCLT, PFCLT, PLT, SPET, SLT
  REAL(rkind), DIMENSION(:,:),ALLOCATABLE :: &  ! (NSTM,3)
       PRFVT

!     ****** GRAPHIC DATA VARIABLES ******
! TRPLR
  REAL, DIMENSION(NTM) :: GT, GTS
  REAL, DIMENSION(NGM) :: GTR
  REAL, DIMENSION(NLM) :: GBL
  REAL, DIMENSION(NTM,8) :: GYT
  REAL, DIMENSION(NTM,NCTM) :: GVT
  REAL, DIMENSION(NLM,10) :: GBR, GBRH, GBP1, GBAN
  REAL, DIMENSION(:)    ,ALLOCATABLE :: GRM                       ! (NRM)
  REAL, DIMENSION(:)    ,ALLOCATABLE :: GRG                       ! (NRMP)
  REAL, DIMENSION(:,:)  ,ALLOCATABLE :: GJB, GAD                  ! (NRMP,4)
  REAL, DIMENSION(:,:)  ,ALLOCATABLE :: GET                       ! (NRMP,5)
  REAL, DIMENSION(:,:)  ,ALLOCATABLE :: GAK                       ! (NRMP,6)
  REAL, DIMENSION(:,:)  ,ALLOCATABLE :: GYR, GER                  ! (NRMP,8)
  REAL, DIMENSION(:,:)  ,ALLOCATABLE :: GPNB                      ! (4*NRM,10)
  REAL, DIMENSION(:,:,:),ALLOCATABLE :: GVR                       ! (NRMP,NGM,NCGM)
  REAL, DIMENSION(:,:,:),ALLOCATABLE :: GVRT                       ! (NRM,NGM,NCRTM)
  INTEGER              :: NGR, NGT, NGST
  INTEGER,DIMENSION(10):: NLMAX
  CHARACTER(10),DIMENSION(NCTM):: KVT
  CHARACTER(10),DIMENSION(NCRTM):: KVRT

!     ****** EQUILIBRIUM INTERFACE VARIABLES ******
! TREQV1
  REAL(rkind)                 :: BPSEQ,RIPEQ
  REAL(rkind), DIMENSION(:),ALLOCATABLE :: &  ! (NRM)
       RHOTR, PRHO, HJRHO, VTRHO, TRHO, TTRHO, DVRHO, RKPRHO, ABRHO, ABVRHO, &
       ARRHO, AR1RHO, AR2RHO,  RJCB, EPSRHO, BPRHO, RMJRHO, RMNRHO, &
       TTRHOG, DVRHOG, RKPRHOG, ABRHOG, ABVRHOG, ARRHOG, AR1RHOG, AR2RHOG, &
       ABB2RHOG, AIB2RHOG, ARHBRHOG, RMJRHOG, RMNRHOG, RDPVRHOG, &
       PSITRHO, PSIPRHO, PPPRHO, PIQRHO, PIRHO, FACTQ, &
       PVOLRHOG, PSURRHOG, ABB1RHO
  REAL(rkind), DIMENSION(:),ALLOCATABLE :: &  ! (NRMP)
       QRHO

!     ****** LOG FILE NAME ******
! TRCOM2
  CHARACTER(LEN=80) :: KXNDEV,KXNDCG,KXNID
  CHARACTER(LEN=80) :: KDIRW1,KDIRW2

!     ****** ADDITIONAL PART BY HONDA MITSURU ******
! TRADD
  REAL(rkind), PARAMETER          :: VOID = 0.D0
  REAL(rkind), DIMENSION(3,3)     :: FA, FB, FC
  REAL(rkind), DIMENSION(:)      ,ALLOCATABLE :: & ! (NSTM)
       RTM, AMZ, PEXT, PNSSA, PNSA, PTSA
  REAL(rkind), DIMENSION(:,:)    ,ALLOCATABLE :: & ! (NRM,NSTM)
       PEX, SEX
  REAL(rkind), DIMENSION(:,:,:)  ,ALLOCATABLE :: & ! (NRM,NSM,NSM)
       AKDWD, AKDWP, ADDWD, ADDWP
  REAL(rkind), DIMENSION(:,:,:,:),ALLOCATABLE :: & ! (NVM,NVM,4,3)
       VV,DD
  REAL(rkind), DIMENSION(:,:,:,:),ALLOCATABLE :: & ! (NVM,NVM,2,3)
       VI,DI

!     ****** MODEL SELECTION VARIABLES ******
! TRMDS
  REAL(rkind)     :: SUMPBM
  INTEGER:: &
       NEQMAX,NSCMAX,NSTMAX,MDDIAG
  INTEGER, DIMENSION(:)  ,ALLOCATABLE :: & ! (NEQM)
       NSS, NSV, NNS, NST
  INTEGER, DIMENSION(:,:),ALLOCATABLE :: & ! (0:NSTM,0:3)
       NEA

!     ****** NCLASS ******
! TRNCL
  REAL(rkind), DIMENSION(:)    ,ALLOCATABLE :: & !  (NRM)
       AJBSNC, ETANC, AJEXNC
  REAL(rkind), DIMENSION(:,:)  ,ALLOCATABLE :: & !  (NRM,NSM)
       CJBSP, CJBST, ADNCG, AVNCG
  REAL(rkind), DIMENSION(:,:,:),ALLOCATABLE :: & !  (NRM,NSM,NSM)
       AKNCP, AKNCT, ADNCP, ADNCT, AKLP, AKLD, ADLP, ADLD
  REAL(rkind), DIMENSION(:,:,:),ALLOCATABLE :: & !  (NRM,5,NSM)
       RGFLS,   RQFLS

! fixed profile variables
  REAL(rkind):: time_initial_nfixed,time_initial_tfixed

! *** trcom1 ***
  
! TMSLC
  INTEGER              :: NTAMAX, NTXMAX, NTXMAX1
  REAL(rkind)                 :: PNBI

! TRERU
  REAL(rkind), DIMENSION(:,:),ALLOCATABLE :: RTEXU, RTIXU, RNEXU, RTEXEU, RTIXEU, RNEXEU   ! (NTUM,NRUM)

! TRINS
  INTEGER:: INS

  CONTAINS

  SUBROUTINE allocate_trcomm(ierr)

    USE trcomx,ONLY: allocate_trcomx
    IMPLICIT NONE
    integer, intent(out):: ierr
    INTEGER,SAVE:: nrmax_save=0
    INTEGER,SAVE:: nsmax_save=0
    INTEGER,SAVE:: nszmax_save=0
    INTEGER,SAVE:: nsnmax_save=0
    INTEGER,SAVE:: nnbmax_save=0
    INTEGER,SAVE:: nnfmax_save=0
    INTEGER,SAVE:: necmax_save=0
    INTEGER,SAVE:: nlhmax_save=0
    INTEGER,SAVE:: nicmax_save=0

    ierr = 0

    if(nrmax==nrmax_save .and. &
       nsmax==nsmax_save .and. &
       nszmax==nszmax_save .and. &
       nsnmax==nsnmax_save .and. &
       nnbmax==nnbmax_save .and. &
       nnfmax==nnfmax_save .and. &
       necmax==necmax_save .and. &
       nlhmax==nlhmax_save .and. &
       nicmax==nicmax_save .and. &
       ALLOCATED(PNSS)) return
    if(ALLOCATED(PNSS)) call DEALLOCATE_TRCOMM

    NSTMAX  = NSMAX+NSZMAX+NSNMAX
    NEQMAXM = 3*NSTMAX+1
    NVM     = NEQMAXM
    MWM     = 4*NEQMAXM-1
    MLM     = NEQMAXM*NRMAX
    NRMP    = NRMAX+1
    NGLF    = NRMAX
    LDAB    = 6*NEQMAXM
    NFMAX   = NNBMAX+NNFMAX


    ALLOCATE(PNSS(NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(XV(NEQMAXM,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    IF(NFMAX.GT.0) &
    ALLOCATE(YV(NFMAX,NRMAX),AY(NFMAX,NRMAX),Y(NFMAX,NRMAX),STAT=IERR)
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
    ALLOCATE(RW(NRMAX,NFMAX),RNF(NRMAX,NFMAX),RTF(NRMAX,NFMAX),STAT=IERR)
    IF(IERR.NE.0) GOTO 900
    IF(NFMAX.GT.0) &
    ALLOCATE(ANF0(NFMAX),TF0(NFMAX),ANFAV(NFMAX),TFAV(NFMAX),STAT=IERR)
    IF(IERR.NE.0) GOTO 900
    IF(NFMAX.GT.0) &
    ALLOCATE(WFT(NFMAX),STAT=IERR)
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
    ALLOCATE(AJRF(NRMAX),AJBS(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNB_NNBNR(NNBMAX,NRMAX),PNB_NNBNR(NNBMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNB_NSNNBNR(NSMAX,NNBMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PNB_NSNNBNR(NSMAX,NNBMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNB_NSNR(NSMAX,NRMAX),PNB_NSNR(NSMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNB_NR(NRMAX),PNB_NR(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNB_NS(NSMAX),PNB_NS(NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNB_NNB(NNBMAX),PNB_NNB(NNBMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNF_NNFNR(NNFMAX,NRMAX),PNF_NNFNR(NNFMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNF_NSNNFNR(NSMAX,NNFMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PNF_NSNNFNR(NSMAX,NNFMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNF_NSNR(NSMAX,NRMAX),PNF_NSNR(NSMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNF_NR(NRMAX),PNF_NR(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNF_NS(NSMAX),PNF_NS(NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SNF_NNF(NNFMAX),PNF_NNF(NNFMAX),ENF_NNF(NNFMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AJNB_NNBNR(NNBMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PEC_NEC(NECMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PLH_NLH(NLHMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PIC_NIC(NICMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    
    ALLOCATE(SPSC(NRMAX,NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(POH(NRMAX),PRB(NRMAX),PRC(NRMAX),PRSUM(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PRL(NRMAX),PCX(NRMAX),PIE(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SIE(NRMAX),SCX(NRMAX),TSIE(NRMAX),TSCX(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PIN(NRMAX,NSTM),SSIN(NRMAX,NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SPE(NRMAX,NSTM),PRF(NRMAX,NSTM),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SPE_NSNPELNR(NSMAX,NPELMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PNBIN_NNBNR(NNBMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PNBCL_NSNR(NSMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PNBCL_NSNNBNR(NSMAX,NNBMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PNFIN_NNFNR(NNFMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PNFCL_NSNR(NSMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PNFCL_NSNNFNR(NSMAX,NNFMAX,NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PRFV(NRMAX,NSTM,3),AJRFV(NRMAX,3),RGFLX(NRMAX,NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(PNBCLT(NSMAX),PNFCLT(NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(SPSCT(NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    ALLOCATE(ETA(NRMAX),S(NRMAX),ALPHA(NRMAX),RKCV(NRMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(TAUB(NNBMAX,NRMAX),TAUF(NNFMAX,NRMAX),TAUK(NRMAX),STAT=IERR)
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
    ALLOCATE(ADDWD(NRMAX,NSM,NSMAX),ADDWP(NRMAX,NSMAX,NSMAX),STAT=IERR)
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
    ALLOCATE(CJBSP(NRMAX,NSMAX),CJBST(NRMAX,NSMAX),ADNCG(NRMAX,NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AVNCG(NRMAX,NSMAX) ,STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AKNCP(NRMAX,NSMAX,NSMAX),AKNCT(NRMAX,NSMAX,NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(ADNCP(NRMAX,NSMAX,NSMAX),ADNCT(NRMAX,NSMAX,NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(AKLP(NRMAX,NSMAX,NSMAX),AKLD(NRMAX,NSMAX,NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(ADLP(NRMAX,NSMAX,NSMAX),ADLD(NRMAX,NSMAX,NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900
    ALLOCATE(RGFLS(NRMAX,5,NSMAX),RQFLS(NRMAX,5,NSMAX),STAT=IERR)
      IF(IERR.NE.0) GOTO 900

    CALL allocate_trcomx(NEQMAXM,NRMAX)

    nrmax_save = nrmax
    nsmax_save = nsmax
    nszmax_save = nszmax
    nsnmax_save = nsnmax
    nnbmax_save = nnbmax
    necmax_save = necmax
    nlhmax_save = nlhmax
    nicmax_save = nicmax
    return

 900 continue
    write(6,*) "XX  TRCOMM ALLOCATION ERROR IERR=",ierr
    call DEALLOCATE_ERR_TRCOMM
    return

  END SUBROUTINE ALLOCATE_TRCOMM


  SUBROUTINE DEALLOCATE_TRCOMM

    USE trcomx,ONLY: deallocate_trcomx
    IMPLICIT NONE
    DEALLOCATE(PNSS)
    DEALLOCATE(XV,YV, AY, Y,ZV, AZ, Z,AX,X)
    DEALLOCATE(RG,RM,RHOM,RHOG,BP,RDP,RPSI,RN,RT,RU,RW)
    DEALLOCATE(RNF,RTF,ANC,ANFE,ANNU,ZEFF,PZC,PZFE,BETA,BETAP,BETAL,BETAPL)
    DEALLOCATE(BETAQ,PBM,PADD,VTOR,VPAR,VPRP,VPOL,WROT,ER,VEXB,WEXB,AGMP)
    DEALLOCATE(VEXBP,WEXBP)
    DEALLOCATE(AJ,AJOH, EZOH,QP,AJTOR,AJNB,AJRF,AJBS,QPINV)
    DEALLOCATE(SNB_NNBNR,PNB_NNBNR)
    DEALLOCATE(SNB_NSNNBNR,PNB_NSNNBNR)
    DEALLOCATE(SNB_NSNR,PNB_NSNR,SNB_NR,PNB_NR,SNB_NS,PNB_NS)
    DEALLOCATE(SNF_NNFNR,PNF_NNFNR,ENF_NNF)
    DEALLOCATE(SNF_NSNNFNR,PNF_NSNNFNR)
    DEALLOCATE(SNF_NSNR,PNF_NSNR,SNF_NR,PNF_NR,SNF_NS,PNF_NS)
    DEALLOCATE(PNBIN_NNBNR,PNFIN_NNFNR)
    DEALLOCATE(PNBCl_NSNR,PNFCL_NSNR)
    DEALLOCATE(PNBCl_NSNNBNR,PNFCL_NSNNFNR)
    DEALLOCATE(POH,PRB,PRC,PRL,PRSUM,PCX,PIE)
    DEALLOCATE(AJNB_NNBNR,PEC_NEC,PLH_NLH,PIC_NIC)
    DEALLOCATE(SIE,SCX,TSIE,TSCX,PIN,SSIN,PBCL,SPE,SPE_NSNPELNR)
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

    CALL deallocate_trcomx

    return

  END SUBROUTINE DEALLOCATE_TRCOMM


  SUBROUTINE DEALLOCATE_ERR_TRCOMM

    IF(ALLOCATED(PNSS     ))     DEALLOCATE(PNSS     )
    IF(ALLOCATED(XV       ))     DEALLOCATE(XV       )
    IF(ALLOCATED(YV       ))     DEALLOCATE(YV       )
    IF(ALLOCATED(AY       ))     DEALLOCATE(AY       )
    IF(ALLOCATED(Y        ))     DEALLOCATE(Y        )
    IF(ALLOCATED(ZV       ))     DEALLOCATE(ZV       )
    IF(ALLOCATED(AZ       ))     DEALLOCATE(AZ       )
    IF(ALLOCATED(Z        ))     DEALLOCATE(Z        )
    IF(ALLOCATED(AX       ))     DEALLOCATE(AX       )
    IF(ALLOCATED(X        ))     DEALLOCATE(X        )
    IF(ALLOCATED(RG       ))     DEALLOCATE(RG       )
    IF(ALLOCATED(RM       ))     DEALLOCATE(RM       )
    IF(ALLOCATED(RHOM     ))     DEALLOCATE(RHOM     )
    IF(ALLOCATED(RHOG     ))     DEALLOCATE(RHOG     )
    IF(ALLOCATED(BP       ))     DEALLOCATE(BP       )
    IF(ALLOCATED(RDP      ))     DEALLOCATE(RDP      )
    IF(ALLOCATED(RPSI     ))     DEALLOCATE(RPSI     )
    IF(ALLOCATED(RN       ))     DEALLOCATE(RN       )
    IF(ALLOCATED(RT       ))     DEALLOCATE(RT       )
    IF(ALLOCATED(RU       ))     DEALLOCATE(RU       )
    IF(ALLOCATED(RW       ))     DEALLOCATE(RW       )
    IF(ALLOCATED(RNF      ))     DEALLOCATE(RNF      )
    IF(ALLOCATED(RTF      ))     DEALLOCATE(RTF      )
    IF(ALLOCATED(ANC      ))     DEALLOCATE(ANC      )
    IF(ALLOCATED(ANFE     ))     DEALLOCATE(ANFE     )
    IF(ALLOCATED(ANNU     ))     DEALLOCATE(ANNU     )
    IF(ALLOCATED(ZEFF     ))     DEALLOCATE(ZEFF     )
    IF(ALLOCATED(PZC      ))     DEALLOCATE(PZC      )
    IF(ALLOCATED(PZFE     ))     DEALLOCATE(PZFE     )
    IF(ALLOCATED(BETA     ))     DEALLOCATE(BETA     )
    IF(ALLOCATED(BETAP    ))     DEALLOCATE(BETAP    )
    IF(ALLOCATED(BETAL    ))     DEALLOCATE(BETAL    )
    IF(ALLOCATED(BETAPL   ))     DEALLOCATE(BETAPL   )
    IF(ALLOCATED(BETAQ    ))     DEALLOCATE(BETAQ    )
    IF(ALLOCATED(PBM      ))     DEALLOCATE(PBM      )
    IF(ALLOCATED(PADD     ))     DEALLOCATE(PADD     )
    IF(ALLOCATED(VTOR     ))     DEALLOCATE(VTOR     )
    IF(ALLOCATED(VPAR     ))     DEALLOCATE(VPAR     )
    IF(ALLOCATED(VPRP     ))     DEALLOCATE(VPRP     )
    IF(ALLOCATED(VPOL     ))     DEALLOCATE(VPOL     )
    IF(ALLOCATED(WROT     ))     DEALLOCATE(WROT     )
    IF(ALLOCATED(ER       ))     DEALLOCATE(ER       )
    IF(ALLOCATED(VEXB     ))     DEALLOCATE(VEXB     )
    IF(ALLOCATED(WEXB     ))     DEALLOCATE(WEXB     )
    IF(ALLOCATED(AGMP     ))     DEALLOCATE(AGMP     )
    IF(ALLOCATED(AJ       ))     DEALLOCATE(AJ       )
    IF(ALLOCATED(AJOH     ))     DEALLOCATE(AJOH     )
    IF(ALLOCATED(EZOH     ))     DEALLOCATE(EZOH     )
    IF(ALLOCATED(QP       ))     DEALLOCATE(QP       )
    IF(ALLOCATED(AJTOR    ))     DEALLOCATE(AJTOR    )
    IF(ALLOCATED(AJNB     ))     DEALLOCATE(AJNB     )
    IF(ALLOCATED(AJRF     ))     DEALLOCATE(AJRF     )
    IF(ALLOCATED(AJBS     ))     DEALLOCATE(AJBS     )
    IF(ALLOCATED(POH      ))     DEALLOCATE(POH      )
    IF(ALLOCATED(PRB      ))     DEALLOCATE(PRB      )
    IF(ALLOCATED(PRC      ))     DEALLOCATE(PRC      )
    IF(ALLOCATED(PRL      ))     DEALLOCATE(PRL      )
    IF(ALLOCATED(PRSUM    ))     DEALLOCATE(PRSUM    )
    IF(ALLOCATED(PCX      ))     DEALLOCATE(PCX      )
    IF(ALLOCATED(PIE      ))     DEALLOCATE(PIE      )
    IF(ALLOCATED(SIE      ))     DEALLOCATE(SIE      )
    IF(ALLOCATED(SCX      ))     DEALLOCATE(SCX      )
    IF(ALLOCATED(TSIE     ))     DEALLOCATE(TSIE     )
    IF(ALLOCATED(TSCX     ))     DEALLOCATE(TSCX     )
    IF(ALLOCATED(PIN      ))     DEALLOCATE(PIN      )
    IF(ALLOCATED(SSIN     ))     DEALLOCATE(SSIN     )
    IF(ALLOCATED(PBCL     ))     DEALLOCATE(PBCL     )
    IF(ALLOCATED(SPE      ))     DEALLOCATE(SPE      )
    IF(ALLOCATED(SPE_NSNPELNR))  DEALLOCATE(SPE_NSNPELNR)
    IF(ALLOCATED(PFCL     ))     DEALLOCATE(PFCL     )
    IF(ALLOCATED(PRF      ))     DEALLOCATE(PRF      )
    IF(ALLOCATED(PRFV     ))     DEALLOCATE(PRFV     )
    IF(ALLOCATED(AJRFV    ))     DEALLOCATE(AJRFV    )
    IF(ALLOCATED(RGFLX    ))     DEALLOCATE(RGFLX    )
    IF(ALLOCATED(ETA      ))     DEALLOCATE(ETA      )
    IF(ALLOCATED(S        ))     DEALLOCATE(S        )
    IF(ALLOCATED(ALPHA    ))     DEALLOCATE(ALPHA    )
    IF(ALLOCATED(RKCV     ))     DEALLOCATE(RKCV     )
    IF(ALLOCATED(TAUB     ))     DEALLOCATE(TAUB     )
    IF(ALLOCATED(TAUF     ))     DEALLOCATE(TAUF     )
    IF(ALLOCATED(TAUK     ))     DEALLOCATE(TAUK     )
    IF(ALLOCATED(AK       ))     DEALLOCATE(AK       )
    IF(ALLOCATED(AVK      ))     DEALLOCATE(AVK      )
    IF(ALLOCATED(AD       ))     DEALLOCATE(AD       )
    IF(ALLOCATED(AV       ))     DEALLOCATE(AV       )
    IF(ALLOCATED(AKNC     ))     DEALLOCATE(AKNC     )
    IF(ALLOCATED(AKDW     ))     DEALLOCATE(AKDW     )
    IF(ALLOCATED(ADNC     ))     DEALLOCATE(ADNC     )
    IF(ALLOCATED(ADDW     ))     DEALLOCATE(ADDW     )
    IF(ALLOCATED(AVNC     ))     DEALLOCATE(AVNC     )
    IF(ALLOCATED(AVDW     ))     DEALLOCATE(AVDW     )
    IF(ALLOCATED(AVKNC    ))     DEALLOCATE(AVKNC    )
    IF(ALLOCATED(AVKDW    ))     DEALLOCATE(AVKDW    )
    IF(ALLOCATED(VGR1     ))     DEALLOCATE(VGR1     )
    IF(ALLOCATED(VGR2     ))     DEALLOCATE(VGR2     )
    IF(ALLOCATED(VGR3     ))     DEALLOCATE(VGR3     )
    IF(ALLOCATED(VGR4     ))     DEALLOCATE(VGR4     )
    IF(ALLOCATED(ANS0     ))     DEALLOCATE(ANS0     )
    IF(ALLOCATED(TS0      ))     DEALLOCATE(TS0      )
    IF(ALLOCATED(ANSAV    ))     DEALLOCATE(ANSAV    )
    IF(ALLOCATED(ANLAV    ))     DEALLOCATE(ANLAV    )
    IF(ALLOCATED(TSAV     ))     DEALLOCATE(TSAV     )
    IF(ALLOCATED(WST      ))     DEALLOCATE(WST      )
    IF(ALLOCATED(PRFT     ))     DEALLOCATE(PRFT     )
    IF(ALLOCATED(PBCLT    ))     DEALLOCATE(PBCLT    )
    IF(ALLOCATED(PFCLT    ))     DEALLOCATE(PFCLT    )
    IF(ALLOCATED(PLT      ))     DEALLOCATE(PLT      )
    IF(ALLOCATED(SPET     ))     DEALLOCATE(SPET     )
    IF(ALLOCATED(SLT      ))     DEALLOCATE(SLT      )
    IF(ALLOCATED(PRFVT    ))     DEALLOCATE(PRFVT    )
    IF(ALLOCATED(GRM      ))     DEALLOCATE(GRM      )
    IF(ALLOCATED(GRG      ))     DEALLOCATE(GRG      )
    IF(ALLOCATED(GJB      ))     DEALLOCATE(GJB      )
    IF(ALLOCATED(GAD      ))     DEALLOCATE(GAD      )
    IF(ALLOCATED(GET      ))     DEALLOCATE(GET      )
    IF(ALLOCATED(GAK      ))     DEALLOCATE(GAK      )
    IF(ALLOCATED(GYR      ))     DEALLOCATE(GYR      )
    IF(ALLOCATED(GER      ))     DEALLOCATE(GER      )
    IF(ALLOCATED(GPNB     ))     DEALLOCATE(GPNB     )
    IF(ALLOCATED(GVR      ))     DEALLOCATE(GVR      )
    IF(ALLOCATED(GVRT     ))     DEALLOCATE(GVRT     )
    IF(ALLOCATED(RHOTR    ))     DEALLOCATE(RHOTR    )
    IF(ALLOCATED(PRHO     ))     DEALLOCATE(PRHO     )
    IF(ALLOCATED(HJRHO    ))     DEALLOCATE(HJRHO    )
    IF(ALLOCATED(VTRHO    ))     DEALLOCATE(VTRHO    )
    IF(ALLOCATED(TRHO     ))     DEALLOCATE(TRHO     )
    IF(ALLOCATED(TTRHO    ))     DEALLOCATE(TTRHO    )
    IF(ALLOCATED(PSITRHO  ))     DEALLOCATE(PSITRHO  )
    IF(ALLOCATED(PSIPRHO  ))     DEALLOCATE(PSIPRHO  )
    IF(ALLOCATED(PPPRHO   ))     DEALLOCATE(PPPRHO   )
    IF(ALLOCATED(PIQRHO   ))     DEALLOCATE(PIQRHO   )
    IF(ALLOCATED(PIRHO    ))     DEALLOCATE(PIRHO    )
    IF(ALLOCATED(DVRHO    ))     DEALLOCATE(DVRHO    )
    IF(ALLOCATED(RKPRHO   ))     DEALLOCATE(RKPRHO   )
    IF(ALLOCATED(ABRHO    ))     DEALLOCATE(ABRHO    )
    IF(ALLOCATED(ABVRHO   ))     DEALLOCATE(ABVRHO   )
    IF(ALLOCATED(ARRHO    ))     DEALLOCATE(ARRHO    )
    IF(ALLOCATED(AR1RHO   ))     DEALLOCATE(AR1RHO   )
    IF(ALLOCATED(AR2RHO   ))     DEALLOCATE(AR2RHO   )
    IF(ALLOCATED(RJCB     ))     DEALLOCATE(RJCB     )
    IF(ALLOCATED(EPSRHO   ))     DEALLOCATE(EPSRHO   )
    IF(ALLOCATED(BPRHO    ))     DEALLOCATE(BPRHO    )
    IF(ALLOCATED(RMJRHO   ))     DEALLOCATE(RMJRHO   )
    IF(ALLOCATED(RMNRHO   ))     DEALLOCATE(RMNRHO   )
    IF(ALLOCATED(TTRHOG   ))     DEALLOCATE(TTRHOG   )
    IF(ALLOCATED(DVRHOG   ))     DEALLOCATE(DVRHOG   )
    IF(ALLOCATED(RKPRHOG  ))     DEALLOCATE(RKPRHOG  )
    IF(ALLOCATED(ABRHOG   ))     DEALLOCATE(ABRHOG   )
    IF(ALLOCATED(ABVRHOG  ))     DEALLOCATE(ABVRHOG  )
    IF(ALLOCATED(ARRHOG   ))     DEALLOCATE(ARRHOG   )
    IF(ALLOCATED(AR1RHOG  ))     DEALLOCATE(AR1RHOG  )
    IF(ALLOCATED(AR2RHOG  ))     DEALLOCATE(AR2RHOG  )
    IF(ALLOCATED(ABB2RHOG ))     DEALLOCATE(ABB2RHOG )
    IF(ALLOCATED(AIB2RHOG ))     DEALLOCATE(AIB2RHOG )
    IF(ALLOCATED(ARHBRHOG ))     DEALLOCATE(ARHBRHOG )
    IF(ALLOCATED(RMJRHOG  ))     DEALLOCATE(RMJRHOG  )
    IF(ALLOCATED(RMNRHOG  ))     DEALLOCATE(RMNRHOG  )
    IF(ALLOCATED(PVOLRHOG ))     DEALLOCATE(PVOLRHOG )
    IF(ALLOCATED(PSURRHOG ))     DEALLOCATE(PSURRHOG )
    IF(ALLOCATED(ABB1RHO  ))     DEALLOCATE(ABB1RHO  )
    IF(ALLOCATED(RDPVRHOG ))     DEALLOCATE(RDPVRHOG )
    IF(ALLOCATED(FACTQ    ))     DEALLOCATE(FACTQ    )
    IF(ALLOCATED(QRHO     ))     DEALLOCATE(QRHO     )
    IF(ALLOCATED(RTM      ))     DEALLOCATE(RTM      )
    IF(ALLOCATED(AMZ      ))     DEALLOCATE(AMZ      )
    IF(ALLOCATED(PEXT     ))     DEALLOCATE(PEXT     )
    IF(ALLOCATED(PNSSA    ))     DEALLOCATE(PNSSA    )
    IF(ALLOCATED(PNSA     ))     DEALLOCATE(PNSA     )
    IF(ALLOCATED(PTSA     ))     DEALLOCATE(PTSA     )
    IF(ALLOCATED(PEX      ))     DEALLOCATE(PEX      )
    IF(ALLOCATED(SEX      ))     DEALLOCATE(SEX      )
    IF(ALLOCATED(AKDWD    ))     DEALLOCATE(AKDWD    )
    IF(ALLOCATED(AKDWP    ))     DEALLOCATE(AKDWP    )
    IF(ALLOCATED(ADDWD    ))     DEALLOCATE(ADDWD    )
    IF(ALLOCATED(ADDWP    ))     DEALLOCATE(ADDWP    )
    IF(ALLOCATED(VV       ))     DEALLOCATE(VV       )
    IF(ALLOCATED(DD       ))     DEALLOCATE(DD       )
    IF(ALLOCATED(VI       ))     DEALLOCATE(VI       )
    IF(ALLOCATED(DI       ))     DEALLOCATE(DI       )
    IF(ALLOCATED(NSS      ))     DEALLOCATE(NSS      )
    IF(ALLOCATED(NSV      ))     DEALLOCATE(NSV      )
    IF(ALLOCATED(NNS      ))     DEALLOCATE(NNS      )
    IF(ALLOCATED(NST      ))     DEALLOCATE(NST      )
    IF(ALLOCATED(NEA      ))     DEALLOCATE(NEA      )
    IF(ALLOCATED(AJBSNC   ))     DEALLOCATE(AJBSNC   )
    IF(ALLOCATED(ETANC    ))     DEALLOCATE(ETANC    )
    IF(ALLOCATED(AJEXNC   ))     DEALLOCATE(AJEXNC   )
    IF(ALLOCATED(CJBSP    ))     DEALLOCATE(CJBSP    )
    IF(ALLOCATED(CJBST    ))     DEALLOCATE(CJBST    )
    IF(ALLOCATED(ADNCG    ))     DEALLOCATE(ADNCG    )
    IF(ALLOCATED(AVNCG    ))     DEALLOCATE(AVNCG    )
    IF(ALLOCATED(AKNCP    ))     DEALLOCATE(AKNCP    )
    IF(ALLOCATED(AKNCT    ))     DEALLOCATE(AKNCT    )
    IF(ALLOCATED(ADNCP    ))     DEALLOCATE(ADNCP    )
    IF(ALLOCATED(ADNCT    ))     DEALLOCATE(ADNCT    )
    IF(ALLOCATED(AKLP     ))     DEALLOCATE(AKLP     )
    IF(ALLOCATED(AKLD     ))     DEALLOCATE(AKLD     )
    IF(ALLOCATED(ADLP     ))     DEALLOCATE(ADLP     )
    IF(ALLOCATED(ADLD     ))     DEALLOCATE(ADLD     )
    IF(ALLOCATED(RGFLS    ))     DEALLOCATE(RGFLS    )
    IF(ALLOCATED(RQFLS    ))     DEALLOCATE(RQFLS    )

    RETURN

  END SUBROUTINE DEALLOCATE_ERR_TRCOMM

END MODULE trcomm
