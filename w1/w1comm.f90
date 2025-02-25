MODULE w1comm_parm

  USE bpsd_kinds
  USE bpsd_constants
  IMPLICIT NONE
  
  INTEGER,PARAMETER:: NSM=100  ! maximum number of particle species
  INTEGER,PARAMETER:: NAM=20   ! maximum number of antenna


  INTEGER:: NXMAX,NZMAX,NSMAX,NAMAX
  REAL(rkind):: BB,RR,RA,RD,RB,WALLR,APRFPN,APRFTR,APRFTP
  REAL(rkind):: RF,RKZ,DRF,DRKZ,DXFACT,DXWDTH
  REAL(rkind),DIMENSION(NSM):: PA,PZ,PN,PTPP,PTPR,PU,PNS,PTS,PZCL
  INTEGER,DIMENSION(NSM):: IELEC,IHARM
  REAL(rkind),DIMENSION(NAM):: AJYH,AJZH,APYH,APZH,ALZH,APHH
  REAL(rkind),DIMENSION(NAM):: AJYL,AJZL,APYL,APZL,ALZL,APHL
  REAL(rkind):: RZ,DZ,EPSH,XDMAX,ZEFF,WVYSIZ
  INTEGER:: NPRINT,NFILE,NGRAPH,NLOOP,NSYM
  INTEGER:: NMODEL,NALPHA,NDMAX,NSYS,NGDSP,NCDTYP,NXABS
  INTEGER:: MODELN
  INTEGER:: MDLWG,MDLWGS
  REAL(rkind):: WGZ1,WGZ2,WGAMP,WGNZ
  
END MODULE w1comm_parm

MODULE w1comm

  USE w1comm_parm
  IMPLICIT NONE

  INTEGER:: NXPMAX,NXVMAX,NCMAX,NHMAX,NXANT1,NXANT2
  COMPLEX(rkind):: CGIN(4,5),CGOT(4,5),CFJY1,CFJY2,CFJZ1,CFJZ2
  COMPLEX(rkind):: CFWG1,CFWG2,CFWG3,CFWG4
  COMPLEX(rkind),DIMENSION(:,:,:),ALLOCATABLE:: CE2DA
  COMPLEX(rkind),DIMENSION(:,:,:,:),ALLOCATABLE:: CJ2DA,CJ2DB
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CJ1,CJ2,CJ3,CJ4
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CWG1,CWG2,CWG3,CWG4
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: ZA,AKZ
  INTEGER,DIMENSION(:),ALLOCATABLE:: NZANTYH,NZANTZH,NZANTLH
  INTEGER,DIMENSION(:),ALLOCATABLE:: NZANTYL,NZANTZL,NZANTLL
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: XA,XAM
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: AJCDX,AJCDK
  REAL(rkind):: AJCDT
  INTEGER,DIMENSION(:),ALLOCATABLE:: NXPRNT

  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: PABS2D  ! (NZPM,NXPM,ISM)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: PABS,PABSX  ! (NXPM,ISM)
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: FLUX,FLUXX  ! (NXTM)
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: PABSXZ  ! (ISM)
  COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: CPABS  ! (NXPM,ISM)
  REAL(rkind):: PABSTT,PANT1,PANT2,PANT
  REAL(rkind):: PWALL1,PWALL2,PWALL,PIBW1,PIBW2,PIBW,PERR
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: RANT1,RANT2,XANT1,XANT2 ! (IAM)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: PAK  ! (NZPM,ISM)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: PAKT  ! (NZPM,3)
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CPANTK  ! (NZPM)
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CPANTK1,CPANTK2  ! (NZPM)

  REAL(rkind),DIMENSION(:),ALLOCATABLE:: PROFB ! (NXPM)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       PROFPN,PROFPU,PROFTR,PROFTP ! (NXPM,ISM)

  INTEGER:: MLEN,MWID,MCEN
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CA
  COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: CF
  COMPLEX(rkind),DIMENSION(:,:,:),ALLOCATABLE:: CFS
  COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: CD0,CD1,CD2
  COMPLEX(rkind),DIMENSION(:,:,:),ALLOCATABLE:: CM0,CM1,CM2

  COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: CEF,CBF,CAF,CED
  COMPLEX(rkind),DIMENSION(:,:,:),ALLOCATABLE:: CJD
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: RHL,AHLT
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: FHL,AHL
  COMPLEX(rkind),DIMENSION(:,:,:),ALLOCATABLE:: CJF
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CSKX,CSPX,CSPZ

  REAL,DIMENSION(:),ALLOCATABLE:: GDATA1,GDATA2,GDATA3,GX,GXM,GZ,GDATAZ

  INTEGER,SAVE:: nxmax_save = 0
  INTEGER,SAVE:: nxpmax_save = 0
  INTEGER,SAVE:: nxvmax_save = 0
  INTEGER,SAVE:: nzmax_save = 0
  INTEGER,SAVE:: nsmax_save = 0
  INTEGER,SAVE:: namax_save = 0

CONTAINS

!  PARAMETER (NXPM=5000,NXVM=200,NZLM=5,ISM=5,IAM=8)
!  PARAMETER (NHARMM=30,MATLM=70,NCLM=140*NXPM,NDM=101)

!  PARAMETER (NZPM=2**NZLM)
!  PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
!  PARAMETER (NHM=2*NHARMM+1,NHM1=NHARMM+1)
!  INTEGER:: NXM,NDM,NHM,NHM1,MATLM,NCLM

  SUBROUTINE w1_allocate

    IMPLICIT NONE

    NXPMAX=NXMAX*RA/RB
    NXVMAX=(NXMAX-NXPMAX)/2
    NXPMAX=NXMAX-2*NXVMAX
    
    IF(nxmax == nxmax_save   .AND. &
       nxpmax == nxpmax_save .AND. &
       nxvmax == nxvmax_save .AND. &
       nzmax  == nzmax_save  .AND. &
       nsmax  == nsmax_save  .AND. &
       namax  == namax_save) RETURN

    IF(nxpmax_save.NE.0) CALL w1_deallocate

    ALLOCATE(CE2DA(NZMAX,NXMAX+1,3))
    ALLOCATE(CJ2DA(NZMAX,NXMAX+1,3,NSMAX),CJ2DB(NZMAX,NXMAX+1,3,NSMAX))
    ALLOCATE(CJ1(NZMAX),CJ2(NZMAX),CJ3(NZMAX),CJ4(NZMAX))
    ALLOCATE(CWG1(NZMAX),CWG2(NZMAX),CWG3(NZMAX),CWG4(NZMAX))
    ALLOCATE(ZA(NZMAX),AKZ(NZMAX))
    ALLOCATE(NZANTYH(NAMAX),NZANTZH(NAMAX),NZANTLH(NAMAX))
    ALLOCATE(NZANTYL(NAMAX),NZANTZL(NAMAX),NZANTLL(NAMAX))
    ALLOCATE(XA(NXMAX+1),XAM(NXMAX))
    ALLOCATE(AJCDX(NXMAX),AJCDK(NZMAX))
    ALLOCATE(NXPRNT(NXMAX))

    ALLOCATE(PABS(NXMAX,NSMAX),PABSX(NXMAX,NSMAX),PABSXZ(NSMAX))
    ALLOCATE(CPABS(NXMAX,NSMAX))
    ALLOCATE(PABS2D(NZMAX,NXMAX,NSMAX))
    ALLOCATE(PAK(NZMAX,NSMAX))
    ALLOCATE(PAKT(NZMAX,3))
    ALLOCATE(CPANTK(NZMAX),CPANTK1(NZMAX),CPANTK2(NZMAX))
    ALLOCATE(FLUX(NXMAX+1),FLUXX(NXMAX+1))
    ALLOCATE(RANT1(NAMAX),RANT2(NAMAX),XANT1(NAMAX),XANT2(NAMAX))

    ALLOCATE(PROFB(NXMAX+1))
    ALLOCATE(PROFPN(NXMAX+1,NSMAX),PROFPU(NXMAX+1,NSMAX))
    ALLOCATE(PROFTR(NXMAX+1,NSMAX),PROFTP(NXMAX+1,NSMAX))

    ALLOCATE(CA(6*NXMAX+10))
    ALLOCATE(CD0(4,NXMAX),CD1(2,NXMAX),CD2(4,NXMAX))
    ALLOCATE(CM0(4,NXMAX,NSMAX),CM1(2,NXMAX,NSMAX),CM2(4,NXMAX,NSMAX))

    ALLOCATE(CEF(NXMAX,3),CBF(NXMAX,3),CAF(NXMAX,4))
    ALLOCATE(CED(NXMAX,3),CJD(NXMAX,NSMAX,4))
    ALLOCATE(RHL(NXMAX,7),CJF(NXMAX,NSMAX,7),FHL(NXMAX,NSMAX,2))
    ALLOCATE(AHLT(NSMAX,4),AHL(NXMAX,NSMAX,4))
    ALLOCATE(CSKX(NXMAX*3),CSPX(NXMAX*3),CSPZ(NXMAX*3))

    ALLOCATE(GDATA1(NXMAX),GDATA2(NXMAX),GDATA3(NXMAX))
    ALLOCATE(GX(NXMAX),GXM(NXMAX),GZ(NZMAX+1),GDATAZ(NZMAX+1))


    nxmax_save = nxmax
    nxpmax_save = nxpmax
    nxvmax_save = nxvmax
    nzmax_save = nzmax
    nsmax_save = nsmax
    namax_save = namax

  END SUBROUTINE w1_allocate

  SUBROUTINE w1_deallocate
    
    IMPLICIT NONE

    IF(.NOT.ALLOCATED(CE2DA)) RETURN

    DEALLOCATE(CE2DA,CJ2DA,CJ2DB)
    DEALLOCATE(CJ1,CJ2,CJ3,CJ4)
    DEALLOCATE(CWG1,CWG2,CWG3,CWG4)
    DEALLOCATE(ZA,AKZ)
    DEALLOCATE(NZANTYH,NZANTZH,NZANTLH)
    DEALLOCATE(NZANTYL,NZANTZL,NZANTLL)
    DEALLOCATE(XA,XAM)
    DEALLOCATE(AJCDX,AJCDK)
    DEALLOCATE(NXPRNT)
    
    DEALLOCATE(PABS,PABSX,PABSXZ)
    DEALLOCATE(PABS2D)
    DEALLOCATE(PAK)
    DEALLOCATE(PAKT)
    DEALLOCATE(CPANTK,CPANTK1,CPANTK2)
    DEALLOCATE(FLUX,FLUXX)
    DEALLOCATE(RANT1,RANT2,XANT1,XANT2)

    DEALLOCATE(PROFB)
    DEALLOCATE(PROFPN,PROFPU,PROFTR,PROFTP)

    DEALLOCATE(CA)
    DEALLOCATE(CD0,CD1,CD2)
    DEALLOCATE(CM0,CM1,CM2)

    DEALLOCATE(CEF,CBF,CAF)
    DEALLOCATE(CED,CJD)
    DEALLOCATE(RHL,CJF,FHL)
    DEALLOCATE(AHLT,AHL)
    DEALLOCATE(CSKX,CSPX,CSPZ)

    DEALLOCATE(GDATA1,GDATA2,GDATA3)
    DEALLOCATE(GX,GXM,GZ,GDATAZ)

    nxpmax_save = 0
    nxvmax_save = 0
    nzmax_save = 0
    nsmax_save = 0
    namax_save = 0

  END SUBROUTINE w1_deallocate

END MODULE w1comm

