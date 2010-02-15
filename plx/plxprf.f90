!     $Id$

  MODULE plxprf

    USE bpsd_kinds

!     NXPRF : Maximum number of spatial points read from external file
!     NXSPC : Maximum number of species read from external file
    
    INTEGER(ikind),PARAMETER:: NXPRF=101,NXSPC=6

    INTEGER(ikind):: NPRF
    REAL(rkind),DIMENSION(NXPRF):: PRFRHO,DERIV
    REAL(rkind),DIMENSION(4,NXPRF,NXSPC):: UPRFN,UPRFT

  END MODULE plxprf
