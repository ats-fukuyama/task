MODULE plcomm

      USE bpsd_kinds
      USE bpsd_constants

      INTEGER,PARAMETER:: NSM=100 ! Maximum number of particle species

      INTEGER:: NSMAX,MODELG,MODELN,MODELQ,IDEBUG,MODEFR,MODEFW

      REAL(rkind):: RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ
      REAL(rkind):: PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2
      REAL(rkind):: RHOMIN,QMIN,RHOEDG,RHOITB,RHOGMN,RHOGMX

      REAL(rkind),DIMENSION(NSM):: & 
           PA,PZ,PZ0,IDION,PN,PNS,PTPR,PTPP,PTS,PU,PUS, &
           PNITB,PTITB,PUITB,PZCL

      CHARACTER(len=80):: KNAMEQ,KNAMWR,KNAMFP,KNAMWM,KNAMPF,KNAMFO,KNAMTR
      CHARACTER(len=80):: KNAMEQ2
      CHARACTER(len=1),DIMENSION(NSM)::KIDNS

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

