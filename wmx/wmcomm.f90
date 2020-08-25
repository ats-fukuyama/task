! wmcomm.f90

MODULE wmcomm_parm

!  USE bpsd_kinds
!  USE bpsd_constants
  USE commpi                    ! ncomm,nrank,nsize

! --- plcomm_parm

  USE plcomm_parm

!,ONLY: NSM,NSMAX,MODELG,MODELB,MODELN,MODELQ, &
!       IDEBUG,MODEFR,MODEFW,MODEL_PROF,MODEL_NPROF, &
!       RR,RA,RB,RKAP,RDLT,BB,Q0,QA,RIP,PROFJ, &
!       RHOMIN,QMIN,RHOEDG,PPN0,PTN0, &
!       PA,PZ,PN,PNS,PTPR,PTPP,PTS,PU,PUS,PZCL,NPA,ID_NS,KID_NS, &
!       RHOITB,PNITB,PTITB,PUITB, &
!       PROFN1,PROFN2,PROFT1,PROFT2,PROFU1,PROFU2, &
!       KNAMEQ,KNAMWR,KNAMFP,KNAMWM,KNAMPF,KNAMFO,KNAMTR,KNAMEQ2

! --- dpcomm_parm

  USE dpcomm_parm

!  MODELP(NSM),MODELV(NSM),NCMIN(NSM),NCMAX(NSM)
!  RF0,RFI0,RKX0,RKY0,RKZ0,RX0,RY0,RZ0,RK0,RKANG0
!  EPSRT,LMAXRT
!  NS_NSA_DP(NSM)
!  PMAX(NSM),EMAX(NSM),RHON_MIN(NSM),RHON_MAX(NSM)
!  NPMAX_DP,NTHMAX_DP,NRMAX_DP,NSAMAX_DP

  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER:: NAM=8       ! maximum number of antenna

! --- wm specific input parameters ---  

  INTEGER:: NRMAX       ! number of radial mesh (element) in plcomm.parm
  INTEGER:: NTHMAX      ! number of poloidal mesh         in plcomm.parm

  INTEGER:: NHHMAX      ! number of helical mesh (0 for axisymmetric)
  INTEGER:: NPHMAX      ! number of toroidal mesh (0 for single toroidal mode)
                        !   NPHMAX should be >= NHHMAX*NHC for helical sym.

  INTEGER:: NSUMAX          ! Number of plasma surface plot points
  INTEGER:: NSWMAX          ! Number of wall surface plot points

  REAL(rkind):: B0_FACT     ! B factor for equilibrium data

  REAL(rkind):: RF      ! Real part of wave frequency [MHz]
  REAL(rkind):: RFI     ! Imaginary part of wave frequency [MHz]
  REAL(rkind):: RD      ! Antenna minor radius [m]
  REAL(rkind):: PRFIN   ! Input power (0 for given antenna current) [W]:

  INTEGER:: NTH0        ! Central valude of poloidal mode number
  INTEGER:: NPH0        ! Central valude of toroidal mode number
  INTEGER:: NHC         ! Number of helical coils 

  REAL(rkind):: factor_nth ! Factor of nthmax expansion for couping tensor 
  REAL(rkind):: factor_nhh ! Factor of nhhmax expansion for couping tensor 
  
  INTEGER:: NAMAX           ! number of antenna
  REAL(rkind):: AJ(NAM)     ! Antenna current density [A/m]
  REAL(rkind):: AEWGT(NAM)  ! Poloidal waveguide electric field [V/m]
  REAL(rkind):: AEWGZ(NAM)  ! Toloidal waveguide electric field [V/m]
  REAL(rkind):: APH(NAM)    ! Antenna phase [degree]
  REAL(rkind):: THJ1(NAM)   ! Start poloidal angle of antenna
  REAL(rkind):: THJ2(NAM)   ! End poloidal angle of antenna
  REAL(rkind):: PHJ1(NAM)   ! Start toroidal angle of antenna
  REAL(rkind):: PHJ2(NAM)   ! End toroidal angle of antenna
  REAL(rkind):: ANTANG(NAM) ! Antenna angle to vertical
  REAL(rkind):: BETAJ(NAM)  ! Antenna current profile parameter

  INTEGER:: NPRINT          ! Print output control
  INTEGER:: NGRAPH          ! Graph output contral
  INTEGER:: MODELJ          ! Antenna current model parameter
  INTEGER:: MODELA          ! Alpha particle contribution model parameter
  INTEGER:: MODELM          ! Matrix solver parameter
  INTEGER:: MDLWMK          ! k_paralle toroidal effect
  
  REAL(rkind):: PNA         ! Alpha denisty [10^20 m^-3]
  REAL(rkind):: PNAL        ! Density scale length [m]
  REAL(rkind):: PTA         ! Effective emperature
  REAL(rkind):: ZEFF        ! Effective Z 

  REAL(rkind):: FRMIN       ! Minimum real part of frequency in amplitude scan
  REAL(rkind):: FRMAX       ! Maximum real part of frequency in amplitude scan
  REAL(rkind):: FIMIN       ! Minimum imag part of frequency in amplitude scan
  REAL(rkind):: FIMAX       ! Maximum imag part of frequency in amplitude scan
  REAL(rkind):: FI0         ! Imag part of frequency in 1D amplitude scan

  INTEGER:: NGFMAX          ! Number of real freq mesh in 1D amplitude scan
  INTEGER:: NGXMAX          ! Number of real f. mesh in 2D scan in plcomm_parm
  INTEGER:: NGYMAX          ! Number of imag f. mesh in 2D scan in plasma_parm

  REAL(rkind):: SCMIN       ! Minimum value in parameter scan
  REAL(rkind):: SCMAX       ! Maximum value in parameter scan
  INTEGER:: NSCMAX          ! Number of mesh in parameter scan
  INTEGER:: LISTEG          ! Listing parameter in parameter scan
 
  REAL(rkind):: FRINI       ! Initial real part of frequency in Newton method
  REAL(rkind):: FIINI       ! Initial imag part of frequency in Newton method

  REAL(rkind):: DLTNW       ! Step size in evaluating derivs in Newton method
  REAL(rkind):: EPSNW       ! Convergence criterion in Newton method
  INTEGER:: LMAXNW          ! Maximum iteration count in Newton method
  INTEGER:: LISTNW          ! Listing parameter in Newton method
  INTEGER:: MODENW          ! Type of Newton method

  INTEGER:: NCONT           ! Number of contour lines
  INTEGER:: ILN1            ! Line type of lower contours
  INTEGER:: IBL1            ! Line boldness of lower contours
  INTEGER:: ICL1            ! Line color of lower contours
  INTEGER:: ILN2            ! Line type of higher contours
  INTEGER:: IBL2            ! Line boldness of higher contours
  INTEGER:: ICL2            ! Line color of higher contours

  REAL(rkind):: WAEMIN      ! minimum frequency range in Alfven frequency
  REAL(rkind):: WAEMAX      ! maximum frequency range in Alfven frequency

  INTEGER:: NTHGMAX         ! number of poloidal mesh in graphics

END MODULE wmcomm_parm

MODULE wmcomm

  USE wmcomm_parm
  USE commpi
  IMPLICIT NONE

  REAL(rkind),ALLOCATABLE:: XR(:),XRHO(:)
  INTEGER:: NTHMAX_F    ! number of extended poloidal mesh =NTHMAX*factor_nth
  INTEGER:: NHHMAX_F    ! number of extended helical mesh  =NHHMAX*factor_nhh
  INTEGER:: MDSIZ,MDMIN,MDMAX,LDSIZ,LDMIN,LDMAX
  INTEGER:: MDSIZ_F,MDMIN_F,MDMAX_F,LDSIZ_F,LDMIN_F,LDMAX_F
  INTEGER:: NDSIZ,NDMIN,NDMAX,KDSIZ,KDMIN,KDMAX
  INTEGER:: NDSIZ_F,NDMIN_F,NDMAX_F,KDSIZ_F,KDMIN_F,KDMAX_F
  INTEGER:: MODEWG,MWGMAX

  COMPLEX(rkind),ALLOCATABLE:: CJANT(:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CEWALL(:,:,:)
  INTEGER:: MLEN,MBND
  INTEGER:: istart,iend,nr_start,nr_end
  INTEGER:: MODELK  ! control wave number spectrum model
  INTEGER:: MODEEG  ! indicate eigen mode calculation
  COMPLEX(rkind),ALLOCATABLE:: CFVG(:)
  COMPLEX(rkind),ALLOCATABLE:: CTNSR(:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CGD(:,:,:,:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CPSF(:,:,:,:,:,:,:)

  REAL(rkind),ALLOCATABLE:: RG11(:,:,:),RG12(:,:,:),RG13(:,:,:)
  REAL(rkind),ALLOCATABLE:: RG22(:,:,:),RG23(:,:,:),RG33(:,:,:)
  REAL(rkind),ALLOCATABLE:: RGI11(:,:,:),RGI12(:,:,:),RGI13(:,:,:)
  REAL(rkind),ALLOCATABLE:: RGI22(:,:,:),RGI23(:,:,:),RGI33(:,:,:)
  REAL(rkind),ALLOCATABLE:: RJ(:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CGF11(:,:,:),CGF12(:,:,:),CGF13(:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CGF22(:,:,:),CGF23(:,:,:),CGF33(:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CMAF(:,:,:,:,:),CRMAF(:,:,:,:,:)

  COMPLEX(rkind),ALLOCATABLE:: CEFLD(:,:,:,:),CEFLDK(:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CBFLD(:,:,:,:),CBFLDK(:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CEFLD3D(:,:,:,:),CEFLD3DK(:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CBFLD3D(:,:,:,:),CBFLD3DK(:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CEN(:,:,:,:),CEP(:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CBN(:,:,:,:),CBP(:,:,:,:)

  REAL(rkind),ALLOCATABLE:: PABS(:,:,:,:),PABSK(:,:,:,:),PABSR(:,:,:)
  REAL(rkind),ALLOCATABLE:: PABSRT(:,:),PABST(:,:),PABSTT(:),PABSTS(:)
  COMPLEX(rkind),ALLOCATABLE:: CPABS(:,:,:,:),CPABSK(:,:,:,:),CPABSR(:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CPABS3D(:,:,:,:),CPABS3DK(:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CPABS3DR(:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CFLX(:,:,:),CFLXK(:,:,:),CFLXR(:,:)
  COMPLEX(rkind),ALLOCATABLE:: CFLX3D(:,:,:),CFLX3DK(:,:,:),CFLX3DR(:,:)
  COMPLEX(rkind),ALLOCATABLE:: CFLXRR(:),CFLXT(:)
  COMPLEX(rkind),ALLOCATABLE:: CRADK(:,:,:),CRADKT(:,:),CRADT(:)
  COMPLEX(rkind):: CRADTT
  REAL(rkind),ALLOCATABLE:: PCUR(:,:,:),PCURR(:)
  REAL(rkind):: PCURT

! --- allocatable for FFT
  COMPLEX(rkind),ALLOCATABLE:: CFFT(:)
  REAL(rkind),ALLOCATABLE:: RFFT(:)
  INTEGER,ALLOCATABLE:: LFFT(:)

! --- equilibrium

  REAL(rkind),ALLOCATABLE:: RPST(:,:,:),ZPST(:,:,:)
  REAL(rkind),ALLOCATABLE:: PPST(:,:,:),BPST(:,:,:)
  REAL(rkind),ALLOCATABLE:: BFLD(:,:,:,:)
  REAL(rkind)::RGMIN,RGMAX,ZGMIN,ZGMAX
  REAL(rkind),ALLOCATABLE:: RSU(:,:),ZSU(:,:),RSW(:,:),ZSW(:,:)

  REAL(rkind),ALLOCATABLE:: RHOT(:),PSIP(:)
  REAL(rkind):: PSIPA,PSITA,RAXIS,ZAXIS
  REAL(rkind),ALLOCATABLE:: RPS(:,:),ZPS(:,:)
  REAL(rkind),ALLOCATABLE:: DRPSI(:,:),DZPSI(:,:)
  REAL(rkind),ALLOCATABLE:: DRCHI(:,:),DZCHI(:,:)
  REAL(rkind),ALLOCATABLE:: PPS(:),QPS(:),RBPS(:),VPS(:),RLEN(:)
  REAL(rkind),ALLOCATABLE:: BPR(:,:),BPZ(:,:),BPT(:,:),BTP(:,:)
  REAL(rkind),ALLOCATABLE:: SIOTA(:)
  REAL(rkind),ALLOCATABLE:: RPSG(:,:),ZPSG(:,:)

CONTAINS

  SUBROUTINE wm_allocate
    IMPLICIT NONE
    INTEGER,SAVE:: nrmax_save,nthmax_save,nhhmax_save,nphmax_save
    INTEGER,SAVE:: nsumax_save,nswmax_save,nthgmax_save
    INTEGER,SAVE:: INIT=0
    
    IF(INIT.EQ.0) THEN
       INIT=1
    ELSE
       IF(nrmax.EQ.nrmax_save.AND. &
          nthmax.EQ.nthmax_save.AND. &
          nhhmax.EQ.nhhmax_save.AND. &
          nphmax.EQ.nphmax_save.AND. &
          nsumax.EQ.nsumax_save.AND. &
          nswmax.EQ.nswmax_save.AND. &
          nthgmax.EQ.nthgmax_save) RETURN
       CALL wm_deallocate
    END IF

    IF(modewg.EQ.0) THEN
       mlen=3*nrmax*nthmax*nhhmax
    ELSE
       mlen=3*nrmax*nthmax*nhhmax+mwgmax*namax
    END IF
    mbnd=12*nthmax*nhhmax-1

    ALLOCATE(XR(nrmax),XRHO(nrmax))
    ALLOCATE(CJANT(3,nthmax,nhhmax))
    ALLOCATE(CEWALL(3,nthmax,nhhmax))
    ALLOCATE(CFVG(mlen))
    ALLOCATE(CTNSR(3,3,nthmax_f,nhhmax_f))
    ALLOCATE(CGD(3,3,nthmax_f,nthmax,nhhmax_f,nhhmax,3))
    ALLOCATE(CPSF(3,3,nthmax_f,nthmax,nhhmax_f,nhhmax,3))
    ALLOCATE(RG11(nthmax_f,nhhmax_f,nrmax),RG12(nthmax_f,nhhmax_f,nrmax))
    ALLOCATE(RG13(nthmax_f,nhhmax_f,nrmax),RG22(nthmax_f,nhhmax_f,nrmax))
    ALLOCATE(RG23(nthmax_f,nhhmax_f,nrmax),RG33(nthmax_f,nhhmax_f,nrmax))
    ALLOCATE(RGI11(nthmax_f,nhhmax_f,nrmax),RGI12(nthmax_f,nhhmax_f,nrmax))
    ALLOCATE(RGI13(nthmax_f,nhhmax_f,nrmax),RGI22(nthmax_f,nhhmax_f,nrmax))
    ALLOCATE(RGI23(nthmax_f,nhhmax_f,nrmax),RGI33(nthmax_f,nhhmax_f,nrmax))

    ALLOCATE(RJ(nthmax_f,nhhmax_f,nrmax))
    ALLOCATE(CGF11(nthmax_f,nhhmax_f,3),CGF12(nthmax_f,nhhmax_f,3))
    ALLOCATE(CGF13(nthmax_f,nhhmax_f,3),CGF22(nthmax_f,nhhmax_f,3))
    ALLOCATE(CGF23(nthmax_f,nhhmax_f,3),CGF33(nthmax_f,nhhmax_f,3))
    ALLOCATE(CMAF(3,3,nthmax_f,nhhmax_f,3))
    ALLOCATE(CRMAF(3,3,nthmax_f,nhhmax_f,3))

    ALLOCATE(CEFLD(3,nthmax,nhhmax,nrmax),CEFLDK(3,nthmax,nhhmax,nrmax))
    ALLOCATE(CEFLD3D(3,nthmax,nphmax,nrmax),CEFLD3DK(3,nthmax,nphmax,nrmax))
    ALLOCATE(CBFLD(3,nthmax,nhhmax,nrmax),CBFLDK(3,nthmax,nhhmax,nrmax))
    ALLOCATE(CBFLD3D(3,nthmax,nphmax,nrmax),CBFLD3DK(3,nthmax,nphmax,nrmax))
    ALLOCATE(CEN(3,nthmax,nhhmax,nrmax),CEP(3,nthmax,nhhmax,nrmax))
    ALLOCATE(CBN(3,nthmax,nhhmax,nrmax),CBP(3,nthmax,nhhmax,nrmax))

    ALLOCATE(PABS(nthmax,nhhmax,nrmax,nsmax))
    ALLOCATE(PABSK(nthmax,nhhmax,nrmax,nsmax))
    ALLOCATE(PABSR(nhhmax,nrmax,nsmax))
    ALLOCATE(PABSRT(nrmax,nsmax))
    ALLOCATE(PABST(nphmax,nsmax))
    ALLOCATE(PABSTT(nphmax))
    ALLOCATE(PABSTS(nsmax))
    ALLOCATE(CPABS(nthmax,nhhmax,nrmax,nsmax))
    ALLOCATE(CPABSK(nthmax,nhhmax,nrmax,nsmax))
    ALLOCATE(CPABSR(nhhmax,nrmax,nsmax))
    ALLOCATE(CPABS3D(nthmax,nphmax,nrmax,nsmax))
    ALLOCATE(CPABS3DK(nthmax,nphmax,nrmax,nsmax))
    ALLOCATE(CPABS3DR(nphmax,nrmax,nsmax))
    ALLOCATE(CFLX(nthmax,nhhmax,nrmax))
    ALLOCATE(CFLXK(nthmax,nhhmax,nrmax))
    ALLOCATE(CFLXR(nhhmax,nrmax))
    ALLOCATE(CFLX3D(nthmax,nphmax,nrmax))
    ALLOCATE(CFLX3DK(nthmax,nphmax,nrmax))
    ALLOCATE(CFLX3DR(nphmax,nrmax))
    ALLOCATE(CFLXRR(nrmax))
    ALLOCATE(CFLXT(nphmax))
    ALLOCATE(CRADK(nthmax,nhhmax,namax),CRADKT(nthmax,nhhmax),CRADT(namax))
    ALLOCATE(PCUR(nthmax,nhhmax,nrmax),PCURR(nrmax))

    ALLOCATE(RPST(nthmax_f,nhhmax_f,nrmax))
    ALLOCATE(ZPST(nthmax_f,nhhmax_f,nrmax))
    ALLOCATE(PPST(nthmax_f,nhhmax_f,nrmax))
    ALLOCATE(BPST(nthmax_f,nhhmax_f,nrmax))
    ALLOCATE(BFLD(2:3,nthmax_f,nhhmax_f,nrmax))

    ALLOCATE(RSU(nsumax,nhhmax),ZSU(nsumax,nhhmax))
    ALLOCATE(RSW(nswmax,nhhmax),ZSU(nswmax,nhhmax))

    ALLOCATE(RHOT(nrmax),PSIP(nrmax))
    ALLOCATE(RPS(nthmax_f,nrmax),ZPS(nthmax_f,nrmax))
    ALLOCATE(DRPSI(nthmax_f,nrmax),DZPSI(nthmax_f,nrmax))
    ALLOCATE(DRCHI(nthmax_f,nrmax),DZCHI(nthmax_f,nrmax))
    ALLOCATE(PPS(nrmax),QPS(nrmax),RBPS(nrmax),VPS(nrmax),RLEN(nrmax))
    ALLOCATE(BPR(nthmax_f,nrmax),BPZ(nthmax_f,nrmax))
    ALLOCATE(BPT(nthmax_f,nrmax),BTP(nthmax_f,nrmax))
    ALLOCATE(SIOTA(nrmax))

    ALLOCATE(RPSG(nthgmax,nrmax),ZPSG(nthgmax,nrmax))

    nrmax_save=nrmax
    nthmax_save=nthmax
    nhhmax_save=nhhmax
    nphmax_save=nphmax
    nsumax_save=nsumax
    nswmax_save=nswmax
    nthgmax_save=nthgmax
  END SUBROUTINE wm_allocate

  SUBROUTINE wm_deallocate
    IMPLICIT NONE

    DEALLOCATE(XR,XRHO)
    DEALLOCATE(CJANT)
    DEALLOCATE(CEWALL)
    DEALLOCATE(CFVG)
    DEALLOCATE(CTNSR)
    DEALLOCATE(CGD)
    DEALLOCATE(CPSF)
    DEALLOCATE(RG11,RG12,RG13,RG22,RG23,RG33)
    DEALLOCATE(RGI11,RGI12,RGI13,RGI22,RGI23,RGI33)
    DEALLOCATE(RJ)
    DEALLOCATE(CGF11,CGF12,CGF13,CGF22,CGF23,CGF33)
    DEALLOCATE(CMAF,CRMAF)

    DEALLOCATE(CEFLD,CEFLDK,CEFLD3D,CEFLD3DK)
    DEALLOCATE(CBFLD,CBFLDK,CBFLD3D,CBFLD3DK)
    DEALLOCATE(CEN,CEP,CBN,CBP)

    DEALLOCATE(PABS,PABSK,PABSR)
    DEALLOCATE(CPABS,CPABSK,CPABSR)
    DEALLOCATE(CPABS3D,CPABS3DK,CPABS3DR)
    DEALLOCATE(PABSR,PABST,PABSTT,PABSTS)
    DEALLOCATE(CFLX,CFLXK,CFLXR)
    DEALLOCATE(CFLX3D,CFLX3DK,CFLX3DR,CFLXR,CFLXT)
    DEALLOCATE(CRADK,CRADKT,CRADT)
    DEALLOCATE(PCUR,PCURR)

    DEALLOCATE(RSU,ZSU)
    DEALLOCATE(RSW,ZSW)

    DEALLOCATE(RHOT,PSIP)
    DEALLOCATE(RPST,ZPST,PPST,BPST,BFLD)
    DEALLOCATE(RPS,ZPS,DRPSI,DZPSI,DRCHI,DZCHI)
    DEALLOCATE(PPS,QPS,RBPS,VPS,RLEN)
    DEALLOCATE(BPR,BPZ,BPT,BTP)
    DEALLOCATE(SIOTA)

    DEALLOCATE(RPSG,ZPSG)
  END SUBROUTINE wm_deallocate
END MODULE wmcomm
