! wfcomm.f90

MODULE wfcomm

  !   Define and allocate global variables

  !   nelm_max=>NEMAX
  !   node_max=>NNMAX
  !   node_nside_nelm=>NDELM
  !   xnode=>RNODE
  !   ynode=>ZNODE
  !   xnode_min=>BDRMIN
  !   xnode_max=>BDRMAX
  !   ynode_min=>BDZMIN
  !   ynode_max=>BDZMAX

  !   nseg_max=>NSDMAX
  !   node_nseg=>NDSID
  !   nelm_nseg=>NESID
  !   nside_nseg=>INSID
  !   nseg_nside_nelm=>NSDELM
  !   mode_nseg=>KASID
  !   mode_node=>KANOD
  !   nelm1_nside_nelm=>KNELM
  !   nside1_nside_nelm=>KSELM

  use bpsd
  use plcomm
  use dpcomm_parm
  use commpi
  implicit none

  public

!       --- for parallel computing ---
!  integer :: nrank,nprocs
  integer :: istart,iend
  real(rkind) :: tolerance

!       --- input parameters ---

  integer :: NRM,NZM
  integer :: NZMH

!      ----- graphics -----

  integer :: NBPM
  integer,parameter::  NMDM = 10
  integer,parameter::  NCNM = 12+NMDM

  integer,parameter::   NBM = 100  ! Maximum number of boundary type
  integer,parameter::   NKM = 100  ! Maximum number of material type
  integer,parameter::   NMM = 100  ! Maximum number of medium type
  integer,parameter::   NAM =   8  ! Maximum number of antenna
  integer,parameter::   NJM =2000  ! Maximum number of antenna current position
  integer,parameter::   NCM =  30  ! Maximum number of magnectic coil position
  
!  integer,parameter::  NGXM = 301  ! resolution in x direction
!  integer,parameter::  NGYM = 301  ! resolution in y direction
!  integer,parameter::  NGVM = 301  ! resolution used by 1D plot
  integer,parameter::   NGM =   3
  
  integer,parameter::  NWDM = 12
  integer,parameter::  NCHM = 80

!       --- common variables ---
  complex(rkind):: CII

!       /WFPRM/
  real(rkind):: RF,RKZ
  integer(ikind):: NAMAX,NPH
!  moved to plx
!  real(rkind):: PPN0,PTN0
  real(rkind):: PIN
  integer(ikind):: NPRINT,NDRAWD,NDRAWA,NDRAWE,NGRAPH,NDRAWV
  integer(ikind):: MODELI
  integer(ikind):: MODELD
  REAL(rkind):: sort_weight_x,sort_weight_y
  REAL(rkind):: PSIA
  real(rkind):: R1WG,Z1WG,R2WG,Z2WG,PH1WG,PH2WG,AMPWG,ANGWG,ELPWG,DPHWG
  REAL(rkind):: th_wg_min,th_wg_max,phase_wg_min,phase_wg_cen,phase_wg_max
  REAL(rkind):: gauss_wg
  integer(ikind):: MODELWG
  CHARACTER(LEN=80):: KNAMWG
  real(rkind):: gfactor
  integer(ikind):: MODELWF

  INTEGER(ikind):: model_coll_enhance
  INTEGER(ikind):: model_interpolation
  REAL(rkind):: factor_coll_enhance
  REAL(rkind):: xpos_coll_enhance,xwidth_coll_enhance
  REAL(rkind):: ypos_coll_enhance,ywidth_coll_enhance

!       /WFPRK/
  character(len=32) :: KFNAME,KFNAMA,KFNAMF,KFNAMN
  character(len=32) :: KFNAMB

!       /WFPRD/
  real(rkind):: BDRMIN,BDRMAX,BDZMIN,BDZMAX
  real(rkind):: DELR,DELZ
  real(rkind):: RD,THETJ1,THETJ2
  integer(ikind):: NJMAX
  real(rkind),dimension(NAM):: AJ,APH,APOS,AWD
  integer(ikind):: MDAMP
  real(rkind):: WDAMP,FDAMP,rdamp_min,rdamp_max,zdamp_min,zdamp_max
  real(rkind):: thdamp_min,thdamp_max

  Integer:: idebuga(100)

!       /WFDIV/
  integer(ikind):: iddiv
         
!       /WFELM/
  integer(ikind):: NNMAX,NEMAX,NBNOD,NBSID
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: RNODE,ZNODE !(NNMAX)
                                                ! poisition of node
  integer(ikind),dimension(:)  ,ALLOCATABLE :: KANOD,KBNOD !(NNMAX)
                                                ! if boundary
                                                !  KANOD=1
                                                !  KBNOD=boundary node number
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: SELM        !(NEMAX)
                                                ! area of element
  integer(ikind),dimension(:)  ,ALLOCATABLE :: KAELM !(NEMAX)
                                                ! KAELM=dielectric id number
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: REMIN,ZEMIN !(NEMAX)
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: REMAX,ZEMAX !(NEMAX)
                                                ! range of element area
  integer(ikind),dimension(:,:),ALLOCATABLE :: NDELM       !(3,NEMAX)
                                                ! node number of element
  integer(ikind),dimension(:,:),ALLOCATABLE :: KNELM       !(3,NEMAX)
                                                ! element number of adjascent
  integer(ikind),dimension(:,:),ALLOCATABLE :: KSELM       !(3,NEMAX)
                                                ! sude number of adjascent
  integer(ikind),dimension(:,:),ALLOCATABLE :: NSDELM     !(3,NEMAX)
                                                ! side number of element
  integer(ikind),dimension(:)  ,ALLOCATABLE :: NVNN        !(NNMAX)
                                                ! variable number of node E
        
!       /WFSID/
  integer(ikind):: NSDMAX
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: LSID       !(NSDMAX)
                                                ! length of side
  integer(ikind),dimension(:,:),ALLOCATABLE :: NDSID      !(2,NSDMAX)
                                                ! node number of side
  integer(ikind),dimension(:,:)  ,ALLOCATABLE :: INSID,NESID!(2,NSDMAX) 
                                                ! side number in a element
                                                ! element number of a side
  integer(ikind),dimension(:)  ,ALLOCATABLE :: KASID,KBSID!(NSDMAX) 
                                                ! if boundary
                                                !   KASID=1
                                                !   KBSID: boundary side number
  integer(ikind),dimension(:)  ,ALLOCATABLE :: NVNSD      !(NSDMAX)
                                                ! variable number of side E
        
!       /WFSRT/
  real(rkind)   ,dimension(:),ALLOCATABLE :: SINDEX       !(NEMAX)
  integer(ikind),dimension(:),ALLOCATABLE :: IVELM,IWELM  !(NEMAX)  
  integer(ikind),dimension(:),ALLOCATABLE :: IDELM        !(NEMAX)
  real(rkind)   ,dimension(:),ALLOCATABLE :: SINDEX_MIN   !(NEMAX)
  real(rkind)   ,dimension(:),ALLOCATABLE :: SINDEX_MAX   !(NEMAX)
        
!       /WFMED/
  integer(ikind):: NMMAX,NKMAX
  real(rkind)   ,dimension(:),ALLOCATABLE :: EPSDM,AMUDM,SIGDM !(NMMAX)
  integer(ikind),dimension(:),ALLOCATABLE :: NMKA              !(NKMAX)

  INTEGER:: nbdy_max
  INTEGER,ALLOCATABLE:: nseg_nbdy(:),nbdy_nseg(:)
  
!       /WFBDY/
  integer(ikind):: NBMAX
  integer(ikind),dimension(:)  ,ALLOCATABLE :: KABDY             !(NBMAX)
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: PHIBDY,RESBDY     !(NBMAX)
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: PWRBDY,PHABDY     !(NBMAX)
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: XGBDY,YGBDY,ZGBDY !(NBMAX)
  real(rkind)   ,dimension(:,:),ALLOCATABLE :: XNBDY,YNBDY,ZNBDY !(3,NBMAX)
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: XPBDY,YPBDY,ZPBDY !(NBMAX)
  real(rkind)   ,dimension(:,:),ALLOCATABLE :: SZBDY             !(2,NBMAX)
  integer(ikind),dimension(:)  ,ALLOCATABLE :: NDBDY,NMBDY,NBPMAX!(NBMAX)
  integer(ikind),dimension(:,:),ALLOCATABLE :: NENBP,NDNBP       !(NBPM,NBMAX)

  INTEGER(ikind),DIMENSION(:),ALLOCATABLE:: NSDBS,NNDBS
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CEBSD,CEBND
        
!       /WFSLV/
  integer(ikind):: MLEN,NNBMAX
  integer(ikind):: NMDMAX
  complex(rkind):: CM(6,6)
  complex(rkind),dimension(:)  ,ALLOCATABLE :: CSV    !(MLEN)
  complex(rkind),dimension(:,:),ALLOCATABLE :: CVTOT  !(6,NEMAX)
  complex(rkind),dimension(:)  ,ALLOCATABLE :: CQQ    !(MBND)

!       /WFAIF/
  real(rkind),dimension(3,3,3):: AIF3,AIE3
  real(rkind),dimension(3,3)  :: AIF2,AIE2
  real(rkind),dimension(3)    :: AIF1,AIE1
        
!       /WFFLD/
  complex(rkind),dimension(:)  ,ALLOCATABLE :: CESD   !(NSDMAX)
  complex(rkind),dimension(:)  ,ALLOCATABLE :: CEND   !(NDMAX)
  complex(rkind),dimension(:,:),ALLOCATABLE :: CBF,CBP!(3,NNM)
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: EMAX   !(4)
  real(rkind)::ETMAX,PNMAX
  complex(rkind),dimension(:,:),ALLOCATABLE :: CRFL!(NMDM,NBM)
       
!       /WFPWR/
  real(rkind):: PABSTT
  real(rkind),dimension(NSM):: PABST
!  real(rkind),dimension(:)  ,ALLOCATABLE :: PABSS !(NSM)
!  real(rkind),dimension(:)  ,ALLOCATABLE :: PABSK !(NKMAX)
!  real(rkind),dimension(:)  ,ALLOCATABLE :: PABSTN!(NNMAX)
!  real(rkind),dimension(:,:),ALLOCATABLE :: PABSSN!(NNMAX,NSM)
!  real(rkind),dimension(:,:),ALLOCATABLE :: PABSKN!(NNMAX,NKMAX)
!  real(rkind),dimension(:,:),ALLOCATABLE :: PFV   !(NNMAX,3)
        
!       /WFANT/
  integer(ikind),dimension(NAM)    :: JNUM0
  integer(ikind),dimension(NJM,NAM):: JELMT
  integer(ikind),dimension(NAM)    :: JNUM
  real(rkind),dimension(NJM,NAM)   :: RJ,ZJ
  real(rkind),dimension(NJM,NAM)   :: RJ0,ZJ0
  complex(rkind),dimension(NAM)    :: CIMP
  complex(rkind)                   :: CTIMP
        
!       /WFNAS/
!  real(rkind):: FACT_LEN
!  integer(ikind),dimension(:),ALLOCATABLE :: IDND !(NNMAX)
!  real(rkind)   ,dimension(:),ALLOCATABLE :: EX1WG,EY1WG,EZ1WG!(NBMAX)
!  real(rkind)   ,dimension(:),ALLOCATABLE :: EX2WG,EY2WG,EZ2WG!(NBMAX)
!  real(rkind)   ,dimension(:),ALLOCATABLE :: PWRWG,PHAWG!(NBMAX)
!  integer(ikind),dimension(:),ALLOCATABLE :: IDKA !(NKMAX)
!  integer(ikind),dimension(:),ALLOCATABLE :: IDMAT!(NMMAX)
!  integer(ikind),dimension(:),ALLOCATABLE :: IDBDY!(NBMAX)
!  CHARACTER     ,dimension(:),ALLOCATABLE :: KDKA*25 !(NKMAX)
!  CHARACTER     ,dimension(:),ALLOCATABLE :: KDMAT*25!(NMMAX)
!  CHARACTER     ,dimension(:),ALLOCATABLE :: KDBDY*25!(NBMAX)
!  integer(ikind):: IDNMIN,IDNMAX,IDEMIN,IDEMAX
        
!       /WFWIN/
  real(rkind):: RNDMIN,RNDMAX,ZNDMIN,ZNDMAX
  real(rkind):: LNDMIN,LNDMAX
  integer(ikind):: NFOPEN
  integer(ikind):: NWXMAX
  real,dimension(:,:),ALLOCATABLE :: GZ,GZ_temp !(NGXM,NGYM)
  INTEGER,dimension(:,:),ALLOCATABLE :: IEGZ !(NXVM,NGYM)
  real,dimension(:)  ,ALLOCATABLE :: G2X!(NGXM)
  real,dimension(:)  ,ALLOCATABLE :: G2Y!(NGYM)
  real,dimension(:,:),ALLOCATABLE :: GV !(NGVM,NGM)
  real,dimension(:)  ,ALLOCATABLE :: GX !(NGVM)
  integer(ikind):: NGXMAX,NGYMAX,NGVMAX
        
!       /WFDBG/
  integer(ikind):: NDFILE

! --- allocation flag ---
  integer(ikind):: divinit,elminit,sidinit,srtinit,medinit
  integer(ikind):: srfinit,slvinit,fldinit,pwrinit,nasinit,wininit

! --- fem module variables ---
  INTEGER(ikind):: nxzone_max,nyzone_max

! -----------------------------------------------------------------

contains

!  subroutine wfdiv_allocate

!    use libmtx
!    implicit none
!    integer,save:: NRM_save,NZM_save
    
!    if((NRM.eq.NRM_save).and.&
!       (NZM.eq.NZM_save)) then
!       divinit=0
!       return
!    end if
 
!    allocate(RLEFT(NZM),RRIGHT(NZM),NRA(NZM))
!    allocate(RNDA(NRM,NZM),ZNDA(NRM,NZM))
!    allocate(NDA(NRM,NZM))
    
!    NRM_save = NRM
!    NZM_save = NZM
!    divinit=1

!    return
!  end subroutine wfdiv_allocate
! ----
!  subroutine wfdiv_deallocate
!    implicit none

!    if(divinit.eq.1) then
!       deallocate(RLEFT,RRIGHT,NRA,RNDA,ZNDA,NDA)
!    end if
!    return
!  end subroutine wfdiv_deallocate
! ----
  subroutine wfelm_allocate
    implicit none
    integer,save :: NNMAX_save,NEMAX_save

    if((elminit.eq.0).or.&
      &(elminit.eq.2)) then

       if(elminit.eq.2) then
          if((NNMAX.eq.NNMAX_save).and.&
            &(NEMAX.eq.NEMAX_save)) then
             return
          else
             call wfelm_deallocate
          end if
       end if

       allocate(RNODE(NNMAX),ZNODE(NNMAX),KANOD(NNMAX),KBNOD(NNMAX))
       ALLOCATE(NVNN(NNMAX))
       elminit = 1
       NNMAX_save = NNMAX

    else if(elminit.eq.1) then
       allocate(SELM(NEMAX),KAELM(NEMAX))
       allocate(REMIN(NEMAX),ZEMIN(NEMAX))
       allocate(REMAX(NEMAX),ZEMAX(NEMAX))
       allocate(NDELM(3,NEMAX),KNELM(3,NEMAX),KSELM(3,NEMAX))
       elminit = 2
       NEMAX_save = NEMAX
    end if

    return
  end subroutine wfelm_allocate
!----
  subroutine wfelm_deallocate
    implicit none

    deallocate(RNODE,ZNODE,KANOD,KBNOD,SELM,KAELM)
    deallocate(REMIN,ZEMIN,REMAX,ZEMAX,NDELM,KNELM,KSELM,NVNN)

    return
  end subroutine wfelm_deallocate
!----
  subroutine wfsid_allocate
    implicit none
    integer,save :: NSDMAX_save,NEMAX_save

    if(sidinit.eq.1) then
       if((NSDMAX.eq.NSDMAX_save).and.&
         &( NEMAX.eq. NEMAX_save)) then
          return
       else
          call wfsid_deallocate
       end if
    end if

    allocate(NDSID(2,NSDMAX),KASID(NSDMAX),KBSID(NSDMAX),INSID(2,NSDMAX))
    allocate(NESID(2,NSDMAX),NSDELM(3,NEMAX),LSID(NSDMAX),NVNSD(NSDMAX))

    NSDMAX_save = NSDMAX
    sidinit = 1
    
    return
  end subroutine wfsid_allocate
!-----
  subroutine wfsid_deallocate
    implicit none

    deallocate(NDSID,KASID,KBSID,NSDELM,INSID,NESID,LSID,NVNSD)

    return
  end subroutine wfsid_deallocate
!-----
  subroutine wfsrt_allocate
    implicit none
    integer,save :: NEMAX_save
    
    if(srtinit.eq.1) then
       if(NEMAX.eq.NEMAX_save) then
          return
       else
          call wfsrt_deallocate
       end if
    end if

    allocate(SINDEX(NEMAX),IVELM(NEMAX),IWELM(NEMAX),IDELM(NEMAX))
    allocate(SINDEX_MIN(NEMAX),SINDEX_MAX(NEMAX))

    NEMAX_save = NEMAX
    srtinit = 1
    
    return
  end subroutine wfsrt_allocate
!-----
  subroutine wfsrt_deallocate
    implicit none

    deallocate(SINDEX,IVELM,IWELM,IDELM,SINDEX_MIN,SINDEX_MAX)

    return
  end subroutine wfsrt_deallocate
!-----
  subroutine wfmed_allocate
    implicit none
    
    allocate(EPSDM(NMM),AMUDM(NMM),SIGDM(NMM))
    allocate(NMKA(NKM))

    return
  end subroutine wfmed_allocate
!-----
  subroutine wfmed_deallocate
    implicit none

    IF(ALLOCATED(epsdm)) DEALLOCATE(epsdm)
    IF(ALLOCATED(amudm)) DEALLOCATE(amudm)
    IF(ALLOCATED(sigdm)) DEALLOCATE(sigdm)
    IF(ALLOCATED(nmka)) DEALLOCATE(nmka)

    return
  end subroutine wfmed_deallocate
!-----
!   subroutine wfbdy_allocate
!     implicit none
!     integer,save :: NBMAX_save
!     integer,save :: bdyinit = 0

!     if(bdyinit.eq.0) then
!        NBMAX = NBM
!        NBPM  = NEMAX
!     end if

!     if(bdyinit.eq.1) then
!        if((NBMAX.eq.NBMAX_save)) then
!           return
!        else
!           call wfbdy_deallocate
!        end if
!     end if
    
!     allocate(KABDY(NBMAX),PHIBDY(NBMAX),RESBDY(NBMAX),PWRBDY(NBMAX))
!     allocate(PHABDY(NBMAX),XGBDY(NBMAX),YGBDY(NBMAX),ZGBDY(NBMAX))
!     allocate(XNBDY(3,NBMAX),YNBDY(3,NBMAX),ZNBDY(3,NBMAX))
!     allocate(XPBDY(NBMAX),YPBDY(NBMAX),ZPBDY(NBMAX))
!     allocate(SZBDY(2,NBMAX),NDBDY(NBMAX),NMBDY(NBMAX),NBPMAX(NBMAX))
!     allocate(NENBP(NBPM,NBMAX),NDNBP(NBPM,NBMAX))
  
!     NBMAX_save = NBMAX
!     bdyinit = 1

!     return
!   end subroutine wfbdy_allocate
! !-----
!   subroutine wfbdy_deallocate
!     implicit none

!     deallocate(KABDY,PHIBDY,RESBDY,PWRBDY,PHABDY,XGBDY,YGBDY,ZGBDY)
!     deallocate(XNBDY,YNBDY,ZNBDY,XPBDY,YPBDY,ZPBDY,SZBDY)
!     deallocate(NDBDY,NMBDY,NBPMAX,NENBP,NDNBP)

!     return
!   end subroutine wfbdy_deallocate
!-----
  subroutine wfslv_allocate
    implicit none
    integer,save :: NEMAX_save,MLEN_save
    integer,save :: slvinit = 0

    if(slvinit.eq.1) then
       if((NEMAX .eq. NEMAX_save).and.&
         &(MLEN  .eq.  MLEN_save)) then
          return
       else
          call wfslv_deallocate
       end if
    end if

    allocate(CSV(MLEN),CVTOT(6,NEMAX))
    slvinit = 1
    NEMAX_save = NEMAX
    MLEN_save = MLEN

    return
  end subroutine wfslv_allocate
!-----
  subroutine wfslv_deallocate
    implicit none

    deallocate(CSV,CVTOT)

    return
  end subroutine wfslv_deallocate
!-----
  subroutine wffld_allocate
    implicit none
    integer,save :: NSDMAX_save,NNMAX_save,NMDMAX_save,NBMAX_save

    if(fldinit.eq.1) then
       if((NSDMAX.eq.NSDMAX_save).and.&
         &( NNMAX.eq. NNMAX_save).and.&
         &(NMDMAX.eq.NMDMAX_save).and.&
         &( NBMAX.eq. NBMAX_save)) then
          return
       else
          call wffld_deallocate
       end if
    end if
    allocate(CESD(NSDMAX),CEND(NNMAX))!,CEP(3,NNMAX))
    allocate(CBF(3,NNMAX),CBP(3,NNMAX),EMAX(4),CRFL(NMDMAX,NBMAX))

    NSDMAX_save = NSDMAX
     NNMAX_save = NNMAX
    NMDMAX_save = NMDMAX
     NBMAX_save = NBMAX
    fldinit = 1

    return
  end subroutine wffld_allocate
!-----
  subroutine wffld_deallocate
    implicit none

    deallocate(CESD,CEND,CBF,CBP,EMAX,CRFL)!,CEP)

    return
  end subroutine wffld_deallocate
!-----
!  subroutine wfpwr_allocate
!    implicit none
!    integer,save :: NKMAX_save,NNMAX_save

!    if(pwrinit.eq.1)then
!       if((NKMAX.eq.NKMAX_save).and.&
!         &(NNMAX.eq.NNMAX_save)) then
!          return
!       else
!          call wfpwr_deallocate
!       end if
!    end if
    
!    allocate(PABSS(NSM),PABSK(NKMAX),PABSTN(NNMAX))
!    allocate(PABSSN(NNMAX,NSM),PABSKN(NNMAX,NKMAX),PFV(NNMAX,3))

!    NKMAX_save = NKMAX
!    NNMAX_save = NNMAX
!    pwrinit = 1

!    return
!  end subroutine wfpwr_allocate
!-----
!  subroutine wfpwr_deallocate
!    implicit none

!    deallocate(PABSS,PABSK,PABSTN,PABSSN,PABSKN,PFV)

!    return
!  end subroutine wfpwr_deallocate
!-----
!   subroutine wfnas_allocate
!     implicit none
!     integer,save :: NNMAX_save,NBMAX_save,NKMAX_save,NMMAX_save

!     if(nasinit.eq.0) NBMAX = NBM

!     if(nasinit.eq.1) then
!        if((NNMAX.eq.NNMAX_save).and.&
!          &(NBMAX.eq.NBMAX_save).and.&
!          &(NKMAX.eq.NKMAX_save).and.&
!          &(NMMAX.eq.NMMAX_save)) then
!           return
!        else
!           call wfnas_deallocate
!        end if
!     end if

!     allocate(IDND(NNMAX))
!     allocate(EX1WG(NBMAX),EY1WG(NBMAX),EZ1WG(NBMAX))
!     allocate(EX2WG(NBMAX),EY2WG(NBMAX),EZ2WG(NBMAX))
!     allocate(PWRWG(NBMAX),PHAWG(NBMAX),IDKA(NKMAX))
!     allocate(IDMAT(NMMAX),IDBDY(NBMAX))
!     allocate(KDKA(NKMAX),KDMAT(NMMAX),KDBDY(NBMAX))

!     NNMAX_save = NNMAX
!     NBMAX_save = NBMAX
!     NKMAX_save = NKMAX
!     NMMAX_save = NMMAX
!     nasinit =1

!     return
!   end subroutine wfnas_allocate
! !-----
!   subroutine wfnas_deallocate
!     implicit none
!     deallocate(IDND,EX1WG,EY1WG,EZ1WG,EX2WG,EY2WG,EZ2WG)
!     deallocate(PWRWG,PHAWG,IDKA,IDMAT,IDBDY,KDKA,KDMAT,KDBDY)
!     return
!   end subroutine wfnas_deallocate
!-----
  subroutine wfwin_allocate
    implicit none
    integer,save :: NGXMAX_save,NGYMAX_save,NGVMAX_save

    if(wininit.eq.1) then
       if((NGXMAX.eq.NGXMAX_save).and.&
         &(NGYMAX.eq.NGYMAX_save).and.&
         &(NGVMAX.eq.NGVMAX_save)) then
          return
       else
          call wfwin_deallocate
       end if
    end if

    allocate(GZ(NGXMAX,NGYMAX),GZ_temp(NGXMAX,NGYMAX),G2X(NGXMAX),G2Y(NGYMAX))
    allocate(IEGZ(NGXMAX,NGYMAX))
    allocate(GV(NGVMAX,3),GX(NGVMAX))

    NGXMAX_save = NGXMAX
    NGYMAX_save = NGYMAX
    NGVMAX_save = NGVMAX
    wininit = 1
    
    return
  end subroutine wfwin_allocate
!-----
  subroutine wfwin_deallocate
    implicit none

    deallocate(GZ,GZ_temp,G2X,G2Y,GV,GX,IEGZ)

    return
  end subroutine wfwin_deallocate

end module wfcomm
