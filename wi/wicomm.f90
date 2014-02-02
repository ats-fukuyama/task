!     $Id$

MODULE wicomm

  USE bpsd_kinds
  USE bpsd_constants

! --- Input parameters ---

  INTEGER(ikind):: nxmax,nwmax,modewi,ntaumax
  REAL(rkind)::    xmax,pn0,alfa,aky,beta,taumin,taumax
  COMPLEX(rkind):: cfyn

! --- Global variables ---

  INTEGER(ikind):: MLEN,MWID
  REAL(rkind):: PTOT
  REAL(rkind):: D0(0:1,0:1),D1(0:1,0:1),D2(0:1,0:1),D3(0:1,0:1,0:1)
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: &
       CFY,CSO,CWP,CWE,CPOWER
  COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       CU,CK

CONTAINS

  SUBROUTINE wi_allocate

    IMPLICIT NONE
    INTEGER(ikind),SAVE:: nxmax_save = 0
    INTEGER(ikind),SAVE:: nwmax_save = 0

    IF(nxmax == nxmax_save .AND. &
       nwmax == nwmax_save) RETURN

    CALL wi_deallocate

    mlen=nxmax*2+3
    mwid=4*nwmax+3

    ALLOCATE(CFY(mlen))
    ALLOCATE(CU(2,-nxmax:nxmax))
    ALLOCATE(CK(mwid,mlen),CSO(mlen))
    ALLOCATE(CWP(0:nxmax),CWE(0:nxmax))
    ALLOCATE(CPOWER(0:nxmax))

    nxmax_save=nxmax
    nwmax_save=nwmax

    RETURN
  END SUBROUTINE wi_allocate

  SUBROUTINE wi_deallocate

    IMPLICIT NONE

    IF(ALLOCATED(cfy)) DEALLOCATE(cfy)
    IF(ALLOCATED(cu)) DEALLOCATE(cu)
    IF(ALLOCATED(ck)) DEALLOCATE(ck)
    IF(ALLOCATED(cSO)) DEALLOCATE(cSO)
    IF(ALLOCATED(cwp)) DEALLOCATE(cwp)
    IF(ALLOCATED(cwe)) DEALLOCATE(cwe)
    IF(ALLOCATED(cpower)) DEALLOCATE(cpower)

    RETURN
  END SUBROUTINE wi_deallocate
END MODULE wicomm

MODULE wigcom
  USE bpsd_kinds
  INTEGER(ikind):: N1
  COMPLEX(rkind):: G1,G2,G3,G4,G5
END MODULE wigcom

