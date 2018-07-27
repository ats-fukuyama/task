MODULE w1comm_parm

  USE bpsd_kinds
  USE bpsd_constants
  IMPLICIT NONE
  
  INTEGER,PARAMETER:: NSM=100  ! maximum number of particle species
  INTEGER,PARAMETER:: NAM=20   ! maximum number of antenna


  INTEGER:: NXPMAX,NXVMAX,NZPMAX,NSMAX,NAMAX
  REAL(rkind):: BB,RR,RA,RD,RB,WALLR,APRFPN,APRFTR,APRFTP
  REAL(rkind):: RF,RKZ,DRF,DRKZ,DXFACT,DXWDTH
  REAL(rkind),DIMENSION(NSM):: PA,PZ,PN,PTPP,PTPR,PU,PNS,PTS,PZCL
  INTEGER,DIMENSION(NSM):: IELEC,IHARM
  REAL(rkind),DIMENSION(NAM):: AJYH,AJZH,APYH,APZH,ALZH,APHH
  REAL(rkind),DIMENSION(NAM):: AJYL,AJZL,APYL,APZL,ALZL,APHL
  REAL(rkind):: RZ,DZ,EPSH,XDMAX,ZEFF,WVYSIZ
  INTEGER:: NPRINT,NFILE,NGRAPH,NLOOP,NSYM
  INTEGER:: NMODEL,NALPHA,NDMAX,NSYS,NDISP,NCDTYP,NXABS
  
END MODULE w1comm_parm

MODULE w1comm

  USE w1comm_parm
  IMPLICIT NONE

  INTEGER:: NXTMAX,NCMAX,NHMAX
  COMPLEX(rkind):: CGIN(3,5),CGOT(3,5),CFJY1,CFJY2,CFJZ1,CFJZ2
  COMPLEX(rkind),DIMENSION(:,:,:),ALLOCATABLE:: CE2DA
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CJ1,CJ2,CJ3,CJ4
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
  REAL(rkind):: PABSTT,PANT1,PANT2,PANT
  REAL(rkind):: PWALL1,PWALL2,PWALL,PIBW1,PIBW2,PIBW,PERR
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: RANT1,RANT2,XANT1,XANT2 ! (IAM)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: PAK  ! (NZPM,ISM)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: PAKT  ! (NZPM,3)
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: PANTK  ! (NZPM)

  REAL(rkind),DIMENSION(:),ALLOCATABLE:: PROFB ! (NXPM)
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: &
       PROFPN,PROFPU,PROFTR,PROFTP ! (NXPM,ISM)

  INTEGER:: MLEN,MWID
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CA
  COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: CF
  COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: CD0,CD1,CD2
  COMPLEX(rkind),DIMENSION(:,:,:),ALLOCATABLE:: CM0,CM1,CM2

  COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: CEF,CBF,CAF,CED
  COMPLEX(rkind),DIMENSION(:,:,:),ALLOCATABLE:: CJD
  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: RHL,AHLT
  REAL(rkind),DIMENSION(:,:,:),ALLOCATABLE:: FHL,AHL
  COMPLEX(rkind),DIMENSION(:,:,:),ALLOCATABLE:: CJF
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CSKX,CSPX,CSPZ

  REAL(4),DIMENSION(:),ALLOCATABLE:: GDATA1,GDATA2,GDATA3,GX,GXM,GZ,GDATAZ

  INTEGER,SAVE:: nxpmax_save = 0
  INTEGER,SAVE:: nxvmax_save = 0
  INTEGER,SAVE:: nzpmax_save = 0
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

    IF(nxpmax == nxpmax_save .AND. &
       nxvmax == nxvmax_save .AND. &
       nzpmax == nzpmax_save .AND. &
       nsmax  == nsmax_save  .AND. &
       namax  == namax_save) RETURN

    IF(nxpmax_save.NE.0) CALL w1_deallocate

    NXTMAX=NXPMAX+2*NXVMAX

    ALLOCATE(CE2DA(NZPMAX,NXTMAX,3))
    ALLOCATE(CJ1(NZPMAX),CJ2(NZPMAX),CJ3(NZPMAX),CJ4(NZPMAX))
    ALLOCATE(ZA(NZPMAX),AKZ(NZPMAX))
    ALLOCATE(NZANTYH(NAMAX),NZANTZH(NAMAX),NZANTLH(NAMAX))
    ALLOCATE(NZANTYL(NAMAX),NZANTZL(NAMAX),NZANTLL(NAMAX))
    ALLOCATE(XA(NXTMAX+1),XAM(NXTMAX))
    ALLOCATE(AJCDX(NXPMAX),AJCDK(NZPMAX))
    ALLOCATE(NXPRNT(NXTMAX))

    ALLOCATE(PABS(NXTMAX,NSMAX),PABSX(NXTMAX,NSMAX),PABSXZ(NSMAX))
    ALLOCATE(PABS2D(NZPMAX,NXTMAX,NSMAX))
    ALLOCATE(PAK(NZPMAX,NSMAX))
    ALLOCATE(PAKT(NZPMAX,3))
    ALLOCATE(PANTK(NZPMAX))
    ALLOCATE(FLUX(NXTMAX),FLUXX(NXTMAX))
    ALLOCATE(RANT1(NAMAX),RANT2(NAMAX),XANT1(NAMAX),XANT2(NAMAX))

    ALLOCATE(PROFB(NXTMAX))
    ALLOCATE(PROFPN(NXTMAX,NSMAX),PROFPU(NXTMAX,NSMAX))
    ALLOCATE(PROFTR(NXTMAX,NSMAX),PROFTP(NXTMAX,NSMAX))

    ALLOCATE(CA(6*NXPMAX+10))
    ALLOCATE(CD0(4,NXTMAX),CD1(2,NXTMAX),CD2(4,NXTMAX))
    ALLOCATE(CM0(4,NXTMAX,NSMAX),CM1(2,NXTMAX,NSMAX),CM2(4,NXTMAX,NSMAX))

    ALLOCATE(CEF(NXTMAX,3),CBF(NXTMAX,3),CAF(NXTMAX,4))
    ALLOCATE(CED(NXTMAX,3),CJD(NXTMAX,NSMAX,4))
    ALLOCATE(RHL(NXTMAX,7),CJF(NXTMAX,NSMAX,7),FHL(NXTMAX,NSMAX,2))
    ALLOCATE(AHLT(NSMAX,4),AHL(NXTMAX,NSMAX,4))
    ALLOCATE(CSKX(NXTMAX*3),CSPX(NXTMAX*3),CSPZ(NXTMAX*3))

    ALLOCATE(GDATA1(NXTMAX),GDATA2(NXTMAX),GDATA3(NXTMAX))
    ALLOCATE(GX(NXTMAX),GXM(NXTMAX),GZ(NZPMAX+1),GDATAZ(NZPMAX+1))


    nxpmax_save = nxpmax
    nxvmax_save = nxvmax
    nzpmax_save = nzpmax
    nsmax_save = nsmax
    namax_save = namax

  END SUBROUTINE w1_allocate

  SUBROUTINE w1_deallocate
    
    IMPLICIT NONE

    IF(.NOT.ALLOCATED(CE2DA)) RETURN

    DEALLOCATE(CE2DA)
    DEALLOCATE(CJ1,CJ2,CJ3,CJ4)
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
    DEALLOCATE(PANTK)
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
    nzpmax_save = 0
    nsmax_save = 0
    namax_save = 0

  END SUBROUTINE w1_deallocate

END MODULE w1comm

