MODULE wicomm

  USE bpsd_kinds
  USE bpsd_constants

! --- Input parameters ---

  INTEGER(ikind):: modelg = 0      ! calculation mode: 0:unmag
  INTEGER(ikind):: modelp = 2      ! plasma mode: 0:cold, 1:warm, 2:hot
  REAL(rkind)::    xmin   = -20.D0   ! minimum value of x
  REAL(rkind)::    xmax   = 100.D0 ! maximum value of x
  REAL(rkind)::    pn0    = 1.D0   ! normalized plasma density at x=0.D0
  REAL(rkind)::    alfa   = 0.01D0 ! normalized density gradient 
                                   !    [ alfa = vte / L omega ]
                                   !    [ n=n_0 exp(-x/L) ]
                                   !    [ alfa = c / L omega ]
  REAL(rkind)::    any    = 0.2D0  ! refractive index in the y direction
                                   !    [ any = k_y c / omega ]
  REAL(rkind)::    beta   = 0.1D0  ! ratio of thermal velocity to light veloc.
                                   !    [ beta = vte / c ]
  REAL(rkind)::    pnu    = 0.003D0! normalized collision frequency
                                   !    [ pnu = nu / omega ]

  INTEGER(ikind):: ntaumax= 51     ! number of TAU scan points
  REAL(rkind)::    taumin = 0.D0   ! minimum of TAU scan
  REAL(rkind)::    taumax = 1.65D0 ! maximum of TAU scan
  INTEGER(ikind):: nalfamax= 31    ! number of ALFA scan points
  REAL(rkind)::    alfamin=  0.1D0 ! minimum of ALFA scan (in log step)
  REAL(rkind)::    alfamax= 100.D0 ! maximum of ALFA scan (in log step)
  REAL(rkind)::    xwint  = 10.D0  ! range of kernel integral in vte
  REAL(rkind)::    dx0    = 0.05D0 ! default grid size
  REAL(rkind)::    dxmin  = 0.D0   ! minimum grid size at omegape = omega
  REAL(rkind)::    xwmin  = 1.D0   ! range of reduction near omegape = omega
  COMPLEX(rkind):: cfyn = (1.D0,0.D0) ! E field of incident wave at nx=nxmax
  INTEGER(ikind):: modela = 0      ! modea =0: acceleration included
                                   !        1: acceleration neglected
  INTEGER(ikind):: idebug = 0      ! debug option index
  CHARACTER(len=80):: kfscan=''    ! filename to save scan data

! --- Global variables ---

  INTEGER(ikind):: nxmax ! number of grid points in the x direction
  INTEGER(ikind):: nwmax ! number of grid points for kernel integral
  INTEGER(ikind):: mlmax ! matrix length
  INTEGER(ikind):: mwmax ! matrix width
  INTEGER(ikind):: MLEN  ! coefficient matrix length
  INTEGER(ikind):: MWID  ! coefficient matrix width
  REAL(rkind):: PTOT     ! total absorption power
  REAL(rkind):: D0(0:1,0:1),D1(0:1,0:1),D2(0:1,0:1),D3(0:1,0:1,0:1)
                         ! FEM integral coefficients
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       XGRID
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

