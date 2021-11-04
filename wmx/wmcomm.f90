! wmcomm.f90

MODULE wmcomm_parm

  USE commpi                    ! ncomm,nrank,nsize

! --- plcomm_parm
!  USE bpsd_kinds
!  USE bpsd_constants

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

! --- fixed parameter necessary for parameter input
  
  INTEGER,PARAMETER:: NAM=8         ! maximum number of antenna
  INTEGER,PARAMETER:: idebug_max=99 ! maximum number of idebuga

! --- wm specific input parameters ---  

  INTEGER:: NRMAX       ! number of radial mesh (element) in plcomm.parm
  INTEGER:: NTHMAX      ! number of poloidal mesh         in plcomm.parm

  INTEGER:: NHHMAX      ! number of helical mode (1 for tokamak)
  INTEGER:: NPHMAX      ! number of toroidal mode (valid for tokamak)
                        !                          1: for single mode
  INTEGER:: NPPMAX      ! number of toroidal mode group (valid for helical)
                        !                          1: for single mode group
                        !                        nhc: for full helical modes
                        !    (nphmax will be to set to nhc*nhhmax for helical)

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

  INTEGER:: factor_nth  ! Factor of nthmax expansion for couping tensor 
  INTEGER:: factor_nhh  ! Factor of nhhmax expansion for couping tensor 
  INTEGER:: factor_nph  ! Factor of nphmax expansion for couping tensor 
  
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
  INTEGER:: MDLWMX          ! model id: 0:wm, 1:wm_seki, 2:wmx

  
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

  INTEGER:: nthmax_g        ! number of poloidal mesh for graphics

  INTEGER:: idebuga(idebug_max) ! control of debug info
  CHARACTER(LEN=80):: knam_dump ! matrix and vector dump :idebuga(61)

END MODULE wmcomm_parm

MODULE wmcomm

  USE wmcomm_parm
  USE commpi
  IMPLICIT NONE

  REAL(rkind),ALLOCATABLE:: XR(:),XRHO(:),XTH(:),XTHF(:)
  INTEGER:: NTHMAX_F    ! number of extended poloidal mesh =NTHMAX*factor_nth
  INTEGER:: NHHMAX_F    ! number of extended helical mesh  =NHHMAX*factor_nhh
  INTEGER:: NPHMAX_F    ! number of extended toroidal mesh =NPHMAX*factor_nph
  INTEGER:: NPHTOT      ! number of toroidal mode number for FFT 2**n
  INTEGER:: MDSIZ,MDMIN,MDMAX,LDSIZ,LDMIN,LDMAX
  INTEGER:: MDSIZ_F,MDMIN_F,MDMAX_F,LDSIZ_F,LDMIN_F,LDMAX_F
  INTEGER:: NDSIZ,NDMIN,NDMAX,KDSIZ,KDMIN,KDMAX
  INTEGER:: NDSIZ_F,NDMIN_F,NDMAX_F,KDSIZ_F,KDMIN_F,KDMAX_F
  INTEGER:: MODEWG,MWGMAX
  INTEGER:: NR_S,NBST,NBED,NBMODE  ! for mdlwmx=1

  INTEGER,ALLOCATABLE:: NPH_LOOP(:)
  
  COMPLEX(rkind),ALLOCATABLE:: CJANT(:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CEWALL(:,:,:)
  INTEGER:: MLEN  ! matrix height
  INTEGER:: MBND  ! matrix width
  INTEGER:: MCENT ! matrix center (diagonal position in width)
  INTEGER:: mblock_size ! matrix block size
  INTEGER:: istart,iend ! matrix line range
  INTEGER:: nr_start,nr_end ! radial point range
  INTEGER:: MODELK  ! control wave number spectrum model
  INTEGER:: MODEEG  ! indicate eigen mode calculation
  COMPLEX(rkind),ALLOCATABLE:: CFVG(:)
  COMPLEX(rkind),ALLOCATABLE:: CTNSR(:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CGD(:,:,:,:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CGDD(:,:,:,:,:,:,:)
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
  COMPLEX(rkind),ALLOCATABLE:: CEFLD3D(:,:,:,:),CEFLDK3D(:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CBFLD(:,:,:,:),CBFLDK(:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CBFLD3D(:,:,:,:),CBFLDK3D(:,:,:,:)
  COMPLEX(rkind),ALLOCATABLE:: CEN(:,:,:,:),CEP(:,:,:,:)

  COMPLEX(rkind),ALLOCATABLE:: CPABS(:,:,:,:),CPABSK(:,:,:,:)
  REAL(rkind),ALLOCATABLE:: PABS(:,:,:,:),PABSK(:,:,:,:),PABSR(:,:)
  REAL(rkind),ALLOCATABLE:: PABST(:)

  COMPLEX(rkind),ALLOCATABLE:: CPABS3D(:,:,:,:),CPABSK3D(:,:,:,:)
  
  REAL(rkind),ALLOCATABLE:: PABS3D(:,:,:,:),PABS2D(:,:,:)
  REAL(rkind),ALLOCATABLE:: PABSR3D(:,:,:),PABSKT(:,:,:)
  REAL(rkind),ALLOCATABLE:: PABST3D(:,:),PABSTT3D(:)
  REAL(rkind):: PABSTT

!  REAL(rkind),ALLOCATABLE:: FLUX3DR(:,:,:),FLUX3DTH(:,:,:),FLUX3DPH(:,:,:)
!  REAL(rkind),ALLOCATABLE:: FLUX2DR(:,:),FLUX2DTH(:,:),FLUXR(:)
!  REAL(rkind):: FLUXT
  
  COMPLEX(rkind),ALLOCATABLE:: CPRADK(:,:)
  COMPLEX(rkind):: CPRAD
  REAL(rkind),ALLOCATABLE:: PCUR(:,:,:),PCURR(:)
  REAL(rkind):: PCURT

  REAL(rkind):: AMP_EIGEN

! --- allocatable for FFT
  COMPLEX(rkind),ALLOCATABLE:: CFFT(:)
  REAL(rkind),ALLOCATABLE:: RFFT(:)
  INTEGER,ALLOCATABLE:: LFFT(:)

! --- equilibrium

  REAL(rkind),ALLOCATABLE:: RPST(:,:,:),ZPST(:,:,:)
  REAL(rkind),ALLOCATABLE:: PPST(:,:,:),BPST(:,:,:)
  REAL(rkind),ALLOCATABLE:: BFLD(:,:,:,:)
  REAL(rkind),ALLOCATABLE:: BTHOBN(:)
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

! --- wmxprf

  REAL(rkind),ALLOCATABLE:: PTSSV(:),PNSSV(:)
  REAL(rkind),ALLOCATABLE:: PN60(:,:),PT60(:,:)
  REAL(rkind),ALLOCATABLE:: RNPRF(:,:),RTPRF(:,:)
  CHARACTER(LEN=80):: KNAMEQ_SAVE

! --- graphics

  INTEGER,ALLOCATABLE:: KACONT(:,:,:)

! --- MPI

  INTEGER,ALLOCATABLE:: nr_start_nrank(:),nr_end_nrank(:)
  INTEGER,ALLOCATABLE:: nr_pos_nrank(:),nr_len_nrank(:)
  INTEGER,ALLOCATABLE:: istart_nrank(:),iend_nrank(:)

! --- Interface

  INTERFACE
     FUNCTION GUCLIP(X)
       USE bpsd_kinds
       REAL(rkind):: X
       REAL:: GUCLIP
     END FUNCTION GUCLIP
     FUNCTION NGULEN(X)
       USE bpsd_kinds
       REAL:: X
       INTEGER:: NGULEN
     END FUNCTION NGULEN
  END INTERFACE

CONTAINS

  SUBROUTINE wm_allocate
    IMPLICIT NONE
    INTEGER,SAVE:: nrmax_save,nthmax_save,nhhmax_save,nphmax_save
    INTEGER,SAVE:: nsumax_save,nswmax_save
    INTEGER,SAVE:: INIT=0
    
    mblock_size=3*nthmax*nhhmax
    mbnd= 4*mblock_size-1
    mcent=2*mblock_size
    SELECT CASE(mdlwmx)
    CASE(0,2)
       mlen=mblock_size*nrmax
    CASE(1)
       mlen=mblock_size*(nrmax+2)
    END SELECT

    IF(nthmax.EQ.1) THEN
       nthmax_f=1
    ELSE
       nthmax_f=nthmax*factor_nth
    END IF
    IF(nhhmax.EQ.1) THEN
       nhhmax_f=1
    ELSE
       nhhmax_f=nhhmax*factor_nhh
    END IF
    IF(nphmax.EQ.1) THEN
       nphmax_f=1
    ELSE
       nphmax_f=nphmax*factor_nph
    END IF

    IF(INIT.EQ.0) THEN
       INIT=1
    ELSE
       IF(nrmax.EQ.nrmax_save.AND. &
          nthmax.EQ.nthmax_save.AND. &
          nhhmax.EQ.nhhmax_save.AND. &
          nphmax.EQ.nphmax_save.AND. &
          nsumax.EQ.nsumax_save.AND. &
          nswmax.EQ.nswmax_save) RETURN
       CALL wm_deallocate
    END IF

    ALLOCATE(XR(nrmax+1),XRHO(nrmax+1),XTH(nthmax+1),XTHF(nthmax_f+1))
    ALLOCATE(NPH_LOOP(nphmax))
    ALLOCATE(CJANT(3,nthmax,nhhmax))
    ALLOCATE(CEWALL(3,nthmax,nhhmax))
    ALLOCATE(CFVG(mlen))
    ALLOCATE(CTNSR(3,3,nthmax_f,nhhmax_f))
    ALLOCATE(CGD(3,3,nthmax_f,nthmax,nhhmax_f,nhhmax,3))
    ALLOCATE(CGDD(3,3,nthmax_f,nthmax,nhhmax_f,nhhmax,3))
    ALLOCATE(CPSF(3,3,nthmax_f,nthmax,nhhmax_f,nhhmax,3))
    ALLOCATE(RG11(nthmax_f,nhhmax_f,nrmax+1),RG12(nthmax_f,nhhmax_f,nrmax+1))
    ALLOCATE(RG13(nthmax_f,nhhmax_f,nrmax+1),RG22(nthmax_f,nhhmax_f,nrmax+1))
    ALLOCATE(RG23(nthmax_f,nhhmax_f,nrmax+1),RG33(nthmax_f,nhhmax_f,nrmax+1))
    ALLOCATE(RGI11(nthmax_f,nhhmax_f,nrmax+1),RGI12(nthmax_f,nhhmax_f,nrmax+1))
    ALLOCATE(RGI13(nthmax_f,nhhmax_f,nrmax+1),RGI22(nthmax_f,nhhmax_f,nrmax+1))
    ALLOCATE(RGI23(nthmax_f,nhhmax_f,nrmax+1),RGI33(nthmax_f,nhhmax_f,nrmax+1))

    ALLOCATE(RJ(nthmax_f,nhhmax_f,nrmax+1))
    ALLOCATE(CGF11(nthmax_f,nhhmax_f,3),CGF12(nthmax_f,nhhmax_f,3))
    ALLOCATE(CGF13(nthmax_f,nhhmax_f,3),CGF22(nthmax_f,nhhmax_f,3))
    ALLOCATE(CGF23(nthmax_f,nhhmax_f,3),CGF33(nthmax_f,nhhmax_f,3))
    ALLOCATE(CMAF(3,3,nthmax_f,nhhmax_f,3))
    ALLOCATE(CRMAF(3,3,nthmax_f,nhhmax_f,3))

    ALLOCATE(CEFLDK(3,nthmax,nhhmax,nrmax+1))   ! wm_efield E_rpppr(m,nh,rho)
    ALLOCATE(CEFLD(3,nthmax,nhhmax,nrmax+1))    ! wm_efield E_rpppr(th,hh,rho)
    ALLOCATE(CEN(3,nthmax,nhhmax,nrmax+1))      ! wm_efield E_rpppr(th,hh,rho)
    ALLOCATE(CEP(3,nthmax,nhhmax,nrmax+1))      ! wm_efield E_+-pr(th,hh,rho)
    ALLOCATE(CBFLDK(3,nthmax,nhhmax,nrmax+1))   ! wm_bfield B_rpppr(m,nh,rho)
    ALLOCATE(CBFLD(3,nthmax,nhhmax,nrmax+1))    ! wm_bfield B_rpppr(th,hh,rho)

    ALLOCATE(CPABSK(nthmax_f,nhhmax_f,nrmax+1,nsmax)) ! wm_pabs P(m,nh,rho,ns)
    ALLOCATE(CPABS(nthmax_f,nhhmax_f,nrmax+1,nsmax))  ! wm_pabs P(th,hh,rho,ns)
    ALLOCATE(PCUR(nthmax,nhhmax,nrmax+1))             ! wm_pabs J(th,hh,rho)

    ALLOCATE(PABS(nthmax_f,nhhmax_f,nrmax+1,nsmax))   ! wm_pabs P(th,hh,rho,ns)
    ALLOCATE(PABSK(nthmax_f,nhhmax_f,nrmax+1,nsmax))  ! wm_pabs P(m,nh,rho,ns)
    
    ALLOCATE(PABSR(nrmax+1,nsmax))                    ! wm_pout_sum P(rho,ns)
    ALLOCATE(PABST(nsmax))                            ! wm_pout_sum P(ns)
    ALLOCATE(PCURR(nrmax))                            ! wm_pout_sum J(rho)
    
    ALLOCATE(CEFLD3D(3,nthmax,nphmax,nrmax+1))  ! wm_loop E_rpppr(th,ph,rho)
    ALLOCATE(CEFLDK3D(3,nthmax,nphmax,nrmax+1)) ! wm_loop E_rpppr(m,n,,rho)
    ALLOCATE(CBFLD3D(3,nthmax,nphmax,nrmax+1))  ! wm_loop B_rpppr(th,ph,rho)
    ALLOCATE(CBFLDK3D(3,nthmax,nphmax,nrmax+1)) ! wm_loop B_rpppr(m,n,rho)

    ALLOCATE(CPABS3D(nthmax_f,nphmax_f,nrmax+1,nsmax))  ! wm_loop P(th,ph,r,ns)
    ALLOCATE(CPABSK3D(nthmax_f,nphmax_f,nrmax+1,nsmax)) ! wm_loop P(m,n,rho,ns)

    ALLOCATE(PABS3D(nthmax_f,nphmax_f,nrmax+1,nsmax)) ! wm_loop P(th,ph,rho,ns)
    ALLOCATE(PABS2D(nthmax_f,nrmax+1,nsmax))          ! wm_loop P(th,rho,ns)
    ALLOCATE(PABSR3D(nphmax_f,nrmax+1,nsmax))         ! wm_loop P(n,rho,ns)
    ALLOCATE(PABSKT(nthmax_f,nhhmax_f,nsmax))         ! wm_loop P(m,nh,rho,ns)
    ALLOCATE(PABST3D(nphmax_f,nsmax))                 ! wm_loop P(n,ns)
    ALLOCATE(PABSTT3D(nphmax_f))                      ! wm_loop P(n)

!    ALLOCATE(FLUX3DR(nthmax_f,nphmax_f,nrmax+1))
!    ALLOCATE(FLUX3DTH(nthmax_f,nphmax_f,nrmax+1))
!    ALLOCATE(FLUX3DPH(nthmax_f,nphmax_f,nrmax+1))
!    ALLOCATE(FLUX2DR(nthmax_f,nrmax+1))
!    ALLOCATE(FLUX2DTH(nthmax_f,nrmax+1))
!    ALLOCATE(FLUXR(nrmax+1))

    ALLOCATE(CPRADK(nthmax,nhhmax))                  ! wm_pwrant

    ALLOCATE(RPST(nthmax_f,nhhmax_f,nrmax+1))
    ALLOCATE(ZPST(nthmax_f,nhhmax_f,nrmax+1))
    ALLOCATE(PPST(nthmax_f,nhhmax_f,nrmax+1))
    ALLOCATE(BPST(nthmax_f,nhhmax_f,nrmax+1))
    ALLOCATE(BFLD(2:3,nthmax_f,nhhmax_f,nrmax+1))
    ALLOCATE(BTHOBN(nrmax+1))

    ALLOCATE(RSU(nsumax,nhhmax),ZSU(nsumax,nhhmax))
    ALLOCATE(RSW(nswmax,nhhmax),ZSW(nswmax,nhhmax))

    ALLOCATE(RHOT(nrmax+1),PSIP(nrmax+1))
    ALLOCATE(RPS(nthmax_f,nrmax+1),ZPS(nthmax_f,nrmax+1))
    ALLOCATE(DRPSI(nthmax_f,nrmax+1),DZPSI(nthmax_f,nrmax+1))
    ALLOCATE(DRCHI(nthmax_f,nrmax+1),DZCHI(nthmax_f,nrmax+1))
    ALLOCATE(PPS(nrmax+1),QPS(nrmax+1),RBPS(nrmax+1),VPS(nrmax+1),RLEN(nrmax+1))
    ALLOCATE(BPR(nthmax_f,nrmax+1),BPZ(nthmax_f,nrmax+1))
    ALLOCATE(BPT(nthmax_f,nrmax+1),BTP(nthmax_f,nrmax+1))
    ALLOCATE(SIOTA(nrmax+1))
    ALLOCATE(RPSG(nthmax_g,nrmax+1),ZPSG(nthmax_g,nrmax+1))

    ALLOCATE(PTSSV(nsmax),PNSSV(nsmax))
    PTSSV(1:nsmax)=0.D0
    PNSSV(1:nsmax)=0.D0
    ALLOCATE(PN60(nrmax+1,nsmax),PT60(nrmax+1,nsmax))
    ALLOCATE(RNPRF(nrmax+1,nsmax),RTPRF(nrmax+1,nsmax))

    ALLOCATE(KACONT(8,nrmax+1,nthmax_f))

    nrmax_save=nrmax
    nthmax_save=nthmax
    nhhmax_save=nhhmax
    nphmax_save=nphmax
    nsumax_save=nsumax
    nswmax_save=nswmax
  END SUBROUTINE wm_allocate

  SUBROUTINE wm_deallocate
    IMPLICIT NONE

    DEALLOCATE(XR,XRHO,XTH,XTHF)
    DEALLOCATE(NPH_LOOP)
    DEALLOCATE(CJANT)
    DEALLOCATE(CEWALL)
    DEALLOCATE(CFVG)
    DEALLOCATE(CTNSR)
    DEALLOCATE(CGD,CGDD)
    DEALLOCATE(CPSF)
    DEALLOCATE(RG11,RG12,RG13,RG22,RG23,RG33)
    DEALLOCATE(RGI11,RGI12,RGI13,RGI22,RGI23,RGI33)
    DEALLOCATE(RJ)
    DEALLOCATE(CGF11,CGF12,CGF13,CGF22,CGF23,CGF33)
    DEALLOCATE(CMAF,CRMAF)

    DEALLOCATE(CEFLD,CEFLDK,CEFLD3D,CEFLDK3D)
    DEALLOCATE(CBFLD,CBFLDK,CBFLD3D,CBFLDK3D)
    DEALLOCATE(CEN,CEP)

    DEALLOCATE(CPABS,CPABSK)
    DEALLOCATE(PABS,PABSK,PABSR,PABST)

    DEALLOCATE(CPABS3D,CPABSK3D)
    DEALLOCATE(PABS3D,PABS2D,PABSR3D,PABSKT,PABST3D,PABSTT3D)
!    DEALLOCATE(FLUX3DR,FLUX3DTH,FLUX3DPH,FLUX2DR,FLUX2DTH,FLUXR)
    DEALLOCATE(CPRADK)
    DEALLOCATE(PCUR,PCURR)

    DEALLOCATE(RPST,ZPST,PPST,BPST)
    DEALLOCATE(BFLD,BTHOBN)

    DEALLOCATE(RSU,ZSU)
    DEALLOCATE(RSW,ZSW)

    DEALLOCATE(RHOT,PSIP)
    DEALLOCATE(RPS,ZPS,DRPSI,DZPSI,DRCHI,DZCHI)
    DEALLOCATE(PPS,QPS,RBPS,VPS,RLEN)
    DEALLOCATE(BPR,BPZ,BPT,BTP)
    DEALLOCATE(SIOTA)
    DEALLOCATE(RPSG,ZPSG)
    DEALLOCATE(PNSSV,PTSSV)
    DEALLOCATE(PN60,PT60)
    DEALLOCATE(RNPRF,RTPRF)
    DEALLOCATE(KACONT)

  END SUBROUTINE wm_deallocate
END MODULE wmcomm
