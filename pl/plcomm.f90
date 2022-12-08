MODULE plcomm_type
  USE bpsd_kinds

  TYPE pl_mag_type
     real(rkind):: BABS,BNX,BNY,BNZ,RHON
                           ! BABS: magnetic field strength [T]
                           ! BNX:  normalized X component of B (phi=0 deg)
                           ! BNY:  normalized Y component of B (phi=90 deg)
                           ! BNZ:  normalized Z component of B (vertical)
                           ! RHON: normalized minor radius
  END TYPE pl_mag_type

  TYPE pl_prf_type       ! local plasma parameter
     real(rkind):: RN,RTPR,RTPP,RU,RUPL,RNUC
                           ! RN:   number density [10^{20}m^{-3}]
                           ! RTPR: parallel temperature [keV]
                           ! RTPP: perpendicular temperature [keV]
                           ! RU:   toroidal fluid velocity [m/s]
                           ! RUPL: poloidal fluid velocity [m/s]
                           ! RNUC: collision frequency [1/s]
  END TYPE pl_prf_type

  TYPE pl_grd_type
     real(rkind):: grdn,grdtpr,grdtpp,grdu,grdupl
                           ! GRDNN:  density gradient [1/m]
                           ! GRDTPR: parallel temperature gradient [1/m]
                           ! GRDTPP: perpendicular temperature gradient [1/m]
                           ! GRDU:   toroidal fluid velocity gradient [1/m]
                           ! GRDUPL: poloidal fluid velocity gradient [1/m]
  END TYPE pl_grd_type
END MODULE plcomm_type

MODULE plcomm_parm

      USE bpsd_kinds
      USE bpsd_constants

      INTEGER,PARAMETER:: NSM=100   ! Maximum number of particle species
      INTEGER,PARAMETER:: NCOILM=30 ! Maximum number of mirror coils

      INTEGER:: NSMAX,NCOILMAX
      INTEGER:: MODELG,MODELB,MODELN,MODELQ,IDEBUG,MODEFR,MODEFW
      INTEGER:: mdlplw
      INTEGER:: MODEL_PROF,MODEL_NPROF,model_coll

      REAL(rkind):: RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ
      REAL(rkind):: RMIR,ZBB,Hpitch1,Hpitch2,RRCH
      REAL(rkind),DIMENSION(NCOILM):: RCOIL,ZCOIL,BCOIL
      REAL(rkind):: RHOMIN,QMIN,RHOEDG,RHOGMN,RHOGMX
      REAL(rkind):: PPN0,PTN0,RF_PL
      REAL(rkind):: BAXIS_SCALED
      REAL(rkind),DIMENSION(3):: r_corner,z_corner
      REAL(rkind),DIMENSION(3):: br_corner,bz_corner,bt_corner
      REAL(rkind),DIMENSION(3,NSM):: pn_corner,ptpr_corner,ptpp_corner
      REAL(rkind):: profn_travis_g,profn_travis_h,profn_travis_p, &
           profn_travis_q,profn_travis_w,proft_travis_g,proft_travis_h, &
           proft_travis_p,proft_travis_q,proft_travis_w

      REAL(rkind),DIMENSION(NSM):: & 
           PA,PZ,PN,PNS,PTPR,PTPP,PTS, &
           PU,PUS,PUPR,PUPP, &
           RHOITB,PNITB,PTITB,PUITB, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           PZCL,PNUC
      INTEGER,DIMENSION(NSM)::NPA,ID_NS
      CHARACTER(len=2),DIMENSION(NSM)::KID_NS

      CHARACTER(len=80):: KNAMEQ,KNAMWR,KNAMFP,KNAMWM,KNAMPF,KNAMFO,KNAMTR
      CHARACTER(len=80):: KNAMEQ2

END MODULE plcomm_parm

MODULE plcomm

  USE plcomm_parm
  USE plcomm_type

  REAL(rkind):: HA1

CONTAINS

  SUBROUTINE pl_allocate_ns
    ! DUMMY SUBROUTINE 
    NSMAX=1
    return
  end subroutine pl_allocate_ns

END MODULE plcomm

MODULE plxprf

  USE bpsd_kinds

!     NXPRF : Maximum number of spatial points read from external file
!     NXSPC : Maximum number of species read from external file
    
  INTEGER(ikind),PARAMETER:: NXPRF=101,NXSPC=6

  INTEGER(ikind):: NPRFMAX
  REAL(rkind),DIMENSION(NXPRF):: PRFRHO,DERIV
  REAL(rkind),DIMENSION(4,NXPRF,NXSPC):: UPRFN,UPRFT

END MODULE plxprf

