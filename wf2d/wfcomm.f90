! wfcomm.f90

module wfcomm

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
  !
  !   nbdy_max=NBMAX


  !   Define and allocate global variables

  use bpsd
  use plcomm
  USE dpcomm_parm
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

  INTEGER,PARAMETER::  NLM = 11    ! maximum number of layers in divider
  
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

  REAL(rkind):: rmin_div,rmax_div,thmin_div,thmax_div
  REAL(rkind):: delr_div,delth_div
  
  INTEGER:: nlayer_max
  character(LEN=1):: ch_layer_mode
  REAL(rkind):: posl_nlayer(NLM+1),thickness_nlayer(NLM)
  REAL(rkind):: pos_min,pos_max,step_size
  
  real(rkind):: RD,THETJ1,THETJ2
  integer(ikind):: NJMAX
  real(rkind),dimension(NAM):: AJ,APH,APOS,AWD
  integer(ikind):: MDAMP
  real(rkind):: WDAMP,FDAMP,rdamp_min,rdamp_max,zdamp_min,zdamp_max
  real(rkind):: thdamp_min,thdamp_max

  Integer:: idebug_wf(100)

!       /WFDIV/
  integer(ikind):: mode_div
         
!       /WFELM/
  integer(ikind):: node_max,nelm_max,NBNOD,NBSID
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: RNODE,ZNODE !(node_max)
                                                ! poisition of node
  integer(ikind),dimension(:)  ,ALLOCATABLE :: KANOD,KBNOD !(node_max)
                                                ! if boundary
                                                !  KANOD=1
                                                !  KBNOD=boundary node number
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: SELM        !(nelm_max)
                                                ! area of element
  integer(ikind),dimension(:)  ,ALLOCATABLE :: KAELM !(nelm_max)
                                                ! KAELM=dielectric id number
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: REMIN,ZEMIN !(nelm_max)
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: REMAX,ZEMAX !(nelm_max)
                                                ! range of element area
  integer(ikind),dimension(:,:),ALLOCATABLE :: NDELM       !(3,nelm_max)
                                                ! node number of element
  integer(ikind),dimension(:,:),ALLOCATABLE :: KNELM       !(3,nelm_max)
                                                ! element number of adjascent
  integer(ikind),dimension(:,:),ALLOCATABLE :: KSELM       !(3,nelm_max)
                                                ! sude number of adjascent
  integer(ikind),dimension(:,:),ALLOCATABLE :: NSDELM     !(3,nelm_max)
                                                ! side number of element
  integer(ikind),dimension(:)  ,ALLOCATABLE :: NVNN        !(node_max)
                                                ! variable number of node E
        
!       /WFSID/
  integer(ikind):: NSDMAX
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: LSID       !(NSDMAX)
                                                ! length of side
  integer(ikind),dimension(:,:),ALLOCATABLE :: NDSID      !(2,NSDMAX)
                                                ! node number of side
  integer(ikind),dimension(:)  ,ALLOCATABLE :: INSID,NESID!(NSDMAX) 
                                                ! side number in a element
                                                ! element number of a side
  integer(ikind),dimension(:)  ,ALLOCATABLE :: KASID,KBSID!(NSDMAX) 
                                                ! if boundary
                                                !   KASID=1
                                                !   KBSID: boundary side number
  integer(ikind),dimension(:)  ,ALLOCATABLE :: NVNSD      !(NSDMAX)
                                                ! variable number of side E
        
!       /WFSRT/
  real(rkind)   ,dimension(:),ALLOCATABLE :: SINDEX       !(nelm_max)
  integer(ikind),dimension(:),ALLOCATABLE :: IVELM,IWELM  !(nelm_max)  
  integer(ikind),dimension(:),ALLOCATABLE :: IDELM        !(nelm_max)
  real(rkind)   ,dimension(:),ALLOCATABLE :: SINDEX_MIN   !(nelm_max)
  real(rkind)   ,dimension(:),ALLOCATABLE :: SINDEX_MAX   !(nelm_max)
        
!       /WFMED/
  integer(ikind):: NMMAX,NKMAX
  real(rkind)   ,dimension(:),ALLOCATABLE :: EPSDM,AMUDM,SIGDM !(NMMAX)
  integer(ikind),dimension(:),ALLOCATABLE :: NMKA              !(NKMAX)
        
!       /WFBDY/
!  integer(ikind):: NBMAX
!  integer(ikind),dimension(:)  ,ALLOCATABLE :: KABDY             !(NBMAX)
!  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: PHIBDY,RESBDY     !(NBMAX)
!  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: PWRBDY,PHABDY     !(NBMAX)
!  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: XGBDY,YGBDY,ZGBDY !(NBMAX)
!  real(rkind)   ,dimension(:,:),ALLOCATABLE :: XNBDY,YNBDY,ZNBDY !(3,NBMAX)
!  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: XPBDY,YPBDY,ZPBDY !(NBMAX)
!  real(rkind)   ,dimension(:,:),ALLOCATABLE :: SZBDY             !(2,NBMAX)
!  integer(ikind),dimension(:)  ,ALLOCATABLE :: NDBDY,NMBDY,NBPMAX!(NBMAX)
!  integer(ikind),dimension(:,:),ALLOCATABLE :: NENBP,NDNBP       !(NBPM,NBMAX)

  INTEGER(ikind),DIMENSION(:),ALLOCATABLE:: NSDBS,NNDBS
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CEBSD,CEBND
        
!       /WFSLV/
  integer(ikind):: MLEN !,NNBMAX
  integer(ikind):: NMDMAX
  complex(rkind):: CM(6,6)
  complex(rkind),dimension(:)  ,ALLOCATABLE :: CSV    !(MLEN)
  complex(rkind),dimension(:,:),ALLOCATABLE :: CVTOT  !(6,nelm_max)
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
!  complex(rkind),dimension(:,:),ALLOCATABLE :: CRFL!(NMDM,NBM)
       
!       /WFPWR/
  real(rkind):: PABSTT
  real(rkind),dimension(NSM):: PABST
!  real(rkind),dimension(:)  ,ALLOCATABLE :: PABSS !(NSM)
!  real(rkind),dimension(:)  ,ALLOCATABLE :: PABSK !(NKMAX)
!  real(rkind),dimension(:)  ,ALLOCATABLE :: PABSTN!(node_max)
!  real(rkind),dimension(:,:),ALLOCATABLE :: PABSSN!(node_max,NSM)
!  real(rkind),dimension(:,:),ALLOCATABLE :: PABSKN!(node_max,NKMAX)
!  real(rkind),dimension(:,:),ALLOCATABLE :: PFV   !(node_max,3)
        
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
!  integer(ikind),dimension(:),ALLOCATABLE :: IDND !(node_max)
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
  INTEGER(ikind):: nxzone_max,nyzone_max,ncount_zone_max
  REAL(rkind):: xlen_zone,ylen_zone ! length of rectangular zone in x or y 

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


  SUBROUTINE wf_node_allocate
    IMPLICIT NONE
    INTEGER,SAVE :: node_max_save=0

    IF(node_max.EQ.node_max_save) THEN
       RETURN
    ELSE
       CALL wf_node_deallocate
       node_max_save=0
    END IF

    ALLOCATE(RNODE(node_max),ZNODE(node_max),KANOD(node_max),KBNOD(node_max))
    ALLOCATE(NVNN(node_max))
    node_max_save = node_max
    RETURN
  END SUBROUTINE wf_node_allocate
    
  SUBROUTINE wf_node_deallocate
    IMPLICIT NONE
    IF(.NOT.ALLOCATED(RNODE)) RETURN
    DEALLOCATE(RNODE,ZNODE,KANOD,KBNOD,NVNN)
    RETURN
  END SUBROUTINE wf_node_deallocate

  SUBROUTINE wf_elm_allocate
    IMPLICIT NONE
    INTEGER,SAVE :: nelm_max_save=0

    IF(nelm_max.EQ.nelm_max_save) THEN
       RETURN
    ELSE
       CALL wf_elm_deallocate
       nelm_max_save=0
    END IF

    ALLOCATE(SELM(nelm_max),KAELM(nelm_max))
    ALLOCATE(REMIN(nelm_max),ZEMIN(nelm_max))
    ALLOCATE(REMAX(nelm_max),ZEMAX(nelm_max))
    ALLOCATE(NDELM(3,nelm_max),KNELM(3,nelm_max))
    nelm_max_save = nelm_max
    RETURN
  END SUBROUTINE wf_elm_allocate

  SUBROUTINE wf_elm_deallocate
    IMPLICIT NONE
    IF(.NOT.ALLOCATED(KAELM)) RETURN
    DEALLOCATE(KAELM,REMIN,ZEMIN,REMAX,ZEMAX,NDELM,KNELM)
    RETURN
  END SUBROUTINE wf_elm_deallocate
  
!----
  subroutine wfsid_allocate
    implicit none
    integer,save :: NSDMAX_save=0
    integer,save :: nelm_max_save=0

    if(sidinit.eq.1) then
       if((NSDMAX.eq.NSDMAX_save).and.&
         &( nelm_max.eq. nelm_max_save)) then
          return
       else
          call wfsid_deallocate
          NSDMAX_save=0
          nelm_max_save=0
       end if
    end if

    allocate(NDSID(2,NSDMAX),KASID(NSDMAX),KBSID(NSDMAX),INSID(NSDMAX))
    allocate(NESID(NSDMAX),NSDELM(3,nelm_max),LSID(NSDMAX),NVNSD(NSDMAX))

    NSDMAX_save = NSDMAX
    nelm_max_save = nelm_max
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
    integer,save :: nelm_max_save=0
    
    if(srtinit.eq.1) then
       if(nelm_max.eq.nelm_max_save) then
          return
       else
          call wfsrt_deallocate
          nelm_max_save=0
       end if
    end if

    allocate(SINDEX(nelm_max),IVELM(nelm_max),IWELM(nelm_max),IDELM(nelm_max))
    allocate(SINDEX_MIN(nelm_max),SINDEX_MAX(nelm_max))

    nelm_max_save = nelm_max
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
!        NBPM  = nelm_max
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
    integer,save :: nelm_max_save=0
    integer,save :: MLEN_save=0
    integer,save :: slvinit = 0

    if(slvinit.eq.1) then
       if((nelm_max .eq. nelm_max_save).and.&
         &(MLEN  .eq.  MLEN_save)) then
          return
       else
          call wfslv_deallocate
          nelm_max_save=0
          MLEN_save=0
       end if
    end if

    allocate(CSV(MLEN),CVTOT(6,nelm_max))
    slvinit = 1
    nelm_max_save = nelm_max
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
    integer,save :: NSDMAX_save=0
    integer,save :: node_max_save=0
    integer,save :: NMDMAX_save=0

    if(fldinit.eq.1) then
       if((NSDMAX.eq.NSDMAX_save).and.&
          ( node_max.eq. node_max_save).and.&
          (NMDMAX.eq.NMDMAX_save)) THEN
          return
       else
          call wffld_deallocate
          NSDMAX_save=0
          node_max_save=0
          NMDMAX_save=0
       end if
    end if
    allocate(CESD(NSDMAX),CEND(node_max))!,CEP(3,node_max))
    allocate(CBF(3,node_max),CBP(3,node_max),EMAX(4)) !,CRFL(NMDMAX,NBMAX))

    NSDMAX_save = NSDMAX
    node_max_save  = node_max
    NMDMAX_save = NMDMAX
    fldinit = 1

    return
  end subroutine wffld_allocate
!-----
  subroutine wffld_deallocate
    implicit none

    deallocate(CESD,CEND,CBF,CBP,EMAX)! ,CRFL,CEP)

    return
  end subroutine wffld_deallocate
!-----
!  subroutine wfpwr_allocate
!    implicit none
!    integer,save :: NKMAX_save,node_max_save

!    if(pwrinit.eq.1)then
!       if((NKMAX.eq.NKMAX_save).and.&
!         &(node_max.eq.node_max_save)) then
!          return
!       else
!          call wfpwr_deallocate
!       end if
!    end if
    
!    allocate(PABSS(NSM),PABSK(NKMAX),PABSTN(node_max))
!    allocate(PABSSN(node_max,NSM),PABSKN(node_max,NKMAX),PFV(node_max,3))

!    NKMAX_save = NKMAX
!    node_max_save = node_max
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
!     integer,save :: node_max_save,NBMAX_save,NKMAX_save,NMMAX_save

!     if(nasinit.eq.0) NBMAX = NBM

!     if(nasinit.eq.1) then
!        if((node_max.eq.node_max_save).and.&
!          &(NBMAX.eq.NBMAX_save).and.&
!          &(NKMAX.eq.NKMAX_save).and.&
!          &(NMMAX.eq.NMMAX_save)) then
!           return
!        else
!           call wfnas_deallocate
!        end if
!     end if

!     allocate(IDND(node_max))
!     allocate(EX1WG(NBMAX),EY1WG(NBMAX),EZ1WG(NBMAX))
!     allocate(EX2WG(NBMAX),EY2WG(NBMAX),EZ2WG(NBMAX))
!     allocate(PWRWG(NBMAX),PHAWG(NBMAX),IDKA(NKMAX))
!     allocate(IDMAT(NMMAX),IDBDY(NBMAX))
!     allocate(KDKA(NKMAX),KDMAT(NMMAX),KDBDY(NBMAX))

!     node_max_save = node_max
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
    integer,save :: NGXMAX_save=0
    integer,save :: NGYMAX_save=0
    integer,save :: NGVMAX_save=0

    if(wininit.eq.1) then
       if((NGXMAX.eq.NGXMAX_save).and.&
         &(NGYMAX.eq.NGYMAX_save).and.&
         &(NGVMAX.eq.NGVMAX_save)) then
          return
       else
          call wfwin_deallocate
          NGXMAX_save=0
          NGYMAX_save=0
          NGVMAX_save=0
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
