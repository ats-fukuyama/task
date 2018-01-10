
MODULE dpcomm_parm

  USE plcomm
  IMPLICIT NONE

  PUBLIC

! --- input paramters ---

  INTEGER:: MODELP(NSM),MODELV(NSM),NDISP1(NSM),NDISP2(NSM)
  REAL(rkind):: PMAX(NSM)

  INTEGER:: NPMAX,NTHMAX,NRMAX
  INTEGER:: NGXMAX,NGYMAX,NGPMAX,NSAMAX,NCHMAX,MODEFA
  inTEGER:: NHMAX

  REAL(rkind):: RF0,RFI0,RKX0,RKY0,RKZ0,RX0,RY0,RZ0
  REAL(rkind):: RF1,RFI1,RKX1,RKY1,RKZ1,RX1
  REAL(rkind):: RF2,RFI2,RKX2,RKY2,RKZ2,RX2
  REAL(rkind):: RMIN,RMAX,EPSRT
  INTEGER:: NXMAX,NYMAX,LMAXRT

END MODULE dpcomm_parm

MODULE dpcomm
  USE dpcomm_parm
  IMPLICIT NONE

  PUBLIC

  COMPLEX(rkind):: CRF0,CKX0,CKY0,CKZ0
  REAL(rkind):: XPOS0,YPOS0,ZPOS0
  INTEGER:: ILIST

  REAL(rkind):: RHON_MIN,RHON_MAX

  REAL(rkind):: DELTH,DELR
  REAL(rkind):: PN0,PT0,PTH0
  REAL(rkind),ALLOCATABLE:: FP(:,:,:) ! FP(NTHMAX,NPMAX,NRMAX)
  INTEGER,ALLOCATABLE:: NS_NSA(:) ! NS_NSA(NSAMAX)
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       AEFP,AMFP,RNFP0,RTFP0,DELP ! (NSAMAX)
  REAL(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: FNS
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       FM,DFP,DFT ! (NPMAX,NTHMAX)
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       THM,THG,TSNM,TSNG,TCSM,TCSG,TTNM,TTNG ! (NTHMAX)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       PM,PG ! (NPMAX,NSAMAX)
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       RGMM,RGMG ! (NPMAX)
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       RM ! (NRMAX)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       DGP1,DGP2,DGT1,DGT2 ! (NPM,NTHM)
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       ADJ,ADJD ! (0:NHMAX)
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: &
       CHI,CCHI,SCHI ! (NCHMAX)

CONTAINS

  SUBROUTINE dpfp_allocate
    IMPLICIT NONE
    INTEGER,SAVE:: init=0
    INTEGER,SAVE:: NTHMAX_SAVE=0
    INTEGER,SAVE:: NPMAX_SAVE=0
    INTEGER,SAVE:: NRMAX_SAVE=0
    INTEGER,SAVE:: NSAMAX_SAVE=0
    INTEGER,SAVE:: NCHMAX_SAVE=0

    IF(init.EQ.0) THEN
       init=1
    ELSE
       IF((NTHMAX.EQ.NTHMAX_SAVE).AND. &
          (NPMAX.EQ.NPMAX_SAVE).AND. &
          (NRMAX.EQ.NRMAX_SAVE).AND. &
          (NSAMAX.EQ.NSAMAX_SAVE).AND. &
          (NCHMAX.EQ.NCHMAX_SAVE)) RETURN
       CALL dpfp_deallocate
    ENDIF

    ALLOCATE(NS_NSA(NSAMAX))
    ALLOCATE(AEFP(NSAMAX),AMFP(NSAMAX),RNFP0(NSAMAX), &
             RTFP0(NSAMAX),DELP(NSAMAX))
    ALLOCATE(FNS(NTHMAX,NPMAX,NRMAX,NSAMAX))
    ALLOCATE(FP(NTHMAX,NPMAX,NRMAX))
    ALLOCATE(FM(NPMAX,NTHMAX),DFP(NPMAX,NTHMAX),DFT(NPMAX,NTHMAX))
    ALLOCATE(THM(NTHMAX),THG(NTHMAX),TSNM(NTHMAX),TSNG(NTHMAX), &
             TCSM(NTHMAX),TCSG(NTHMAX),TTNM(NTHMAX),TTNG(NTHMAX))
    ALLOCATE(PM(NPMAX,NSAMAX),PG(NPMAX,NSAMAX))
    ALLOCATE(RGMM(NPMAX),RGMG(NPMAX))
    ALLOCATE(RM(NRMAX))
    ALLOCATE(DGP1(NPMAX,NTHMAX),DGP2(NPMAX,NTHMAX), &
             DGT1(NPMAX,NTHMAX),DGT2(NPMAX,NTHMAX))
    ALLOCATE(CHI(NCHMAX),CCHI(NCHMAX),SCHI(NCHMAX))
    RETURN
  END SUBROUTINE dpfp_allocate

  SUBROUTINE dpfp_deallocate
    DEALLOCATE(NS_NSA)
    DEALLOCATE(AEFP,AMFP,RNFP0,RTFP0,DELP)
    DEALLOCATE(FNS)
    DEALLOCATE(FP)
    DEALLOCATE(FM,DFP,DFT)
    DEALLOCATE(THM,THG,TSNM,TSNG,TCSM,TCSG,TTNM,TTNG)
    DEALLOCATE(PM,PG)
    DEALLOCATE(RGMM,RGMG)
    DEALLOCATE(RM)
    DEALLOCATE(DGP1,DGP2,DGT1,DGT2)
    DEALLOCATE(CHI,CCHI,SCHI)
    RETURN
  END SUBROUTINE dpfp_deallocate

END MODULE dpcomm
