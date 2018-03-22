MODULE plcomm

      USE bpsd_kinds
      USE bpsd_constants

      INTEGER,PARAMETER:: NSM=100 ! Maximum number of particle species

      INTEGER:: NSMAX,MODELG,MODELN,MODELQ,IDEBUG,MODEFR,MODEFW
      INTEGER:: MODEL_PROF,MODEL_NPROF

      REAL(rkind):: RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ
      REAL(rkind):: RHOMIN,QMIN,RHOEDG,RHOGMN,RHOGMX
      REAL(rkind):: PPN0,PTN0,RF_PL
      REAL(rkind),DIMENSION(3):: r_corner,z_corner
      REAL(rkind),DIMENSION(3):: br_corner,bz_corner,bt_corner
      REAL(rkind),DIMENSION(3,NSM):: pn_corner,ptpr_corner,ptpp_corner

      REAL(rkind),DIMENSION(NSM):: & 
           PA,PZ,PZ0,IDION,PN,PNS,PTPR,PTPP,PTS,PU,PUS, &
           RHOITB,PNITB,PTITB,PUITB, &
           PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
           PZCL
      INTEGER,DIMENSION(NSM)::ID_NS
      CHARACTER(len=2),DIMENSION(NSM)::KID_NS

      CHARACTER(len=80):: KNAMEQ,KNAMWR,KNAMFP,KNAMWM,KNAMPF,KNAMFO,KNAMTR
      CHARACTER(len=80):: KNAMEQ2

      CONTAINS

        SUBROUTINE pl_allocate_ns
          ! DUMMY SUBROUTINE 
          NSMAX=1
          return
        end subroutine pl_allocate_ns

      end module plcomm

MODULE plxprf

  USE bpsd_kinds

!     NXPRF : Maximum number of spatial points read from external file
!     NXSPC : Maximum number of species read from external file
    
  INTEGER(ikind),PARAMETER:: NXPRF=101,NXSPC=6

  INTEGER(ikind):: NPRFMAX
  REAL(rkind),DIMENSION(NXPRF):: PRFRHO,DERIV
  REAL(rkind),DIMENSION(4,NXPRF,NXSPC):: UPRFN,UPRFT

END MODULE plxprf

