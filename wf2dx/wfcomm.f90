! wfcomm.f90

!   Define and allocate global variables

MODULE wfcomm_parm

  USE plcomm_parm
  USE dpcomm_parm

  IMPLICIT NONE

  ! === input parameter array size ====
  
  INTEGER,PARAMETER:: NAM =   8 ! maximum number of antenna
  INTEGER,PARAMETER:: NJM =2000 ! maximum number of antenna path points
  INTEGER,PARAMETER:: NMM =  10 ! maximum number of medium type
  INTEGER,PARAMETER:: NKM =  10 ! Maximum number of material type
!  INTEGER,PARAMETER:: NBM =  10 ! Maximum number of boundary type
!  INTEGER,PARAMETER:: NBPM = 10 ! Maximum number of boundary type
  INTEGER,PARAMETER:: idebuga_max=100 ! maximum number of debug parameter array

  ! === controll parameters ====
  
  INTEGER:: MODELI
  CHARACTER(LEN=80):: KFNAME,KFNAMA,KFNAMF,KFNAMN,KFNAMB
  CHARACTER(LEN=80):: KNAMWG

  ! === RF parameters ===
  
  REAL(rkind):: RF,RKZ
  INTEGER:: NPH

  ! === antenna paramters ===

  INTEGER:: nant_max
  REAL(rkind):: AJ(NAM),APH(NAM),AWD(NAM),APOS(NAM)
  REAL(rkind):: PIN(NAM),RD(NAM),THETJ1(NAM),THETJ2(NAM)

  ! === waveguide parameters ===

  INTEGER:: modelwg,modelwf
  REAL(rkind):: x1wg,y1wg,x2wg,y2wg,ph1WG,ph2wg,ampwg,angwg,elpwg,dphwg

  ! === edge dump parameters ===

  INTEGER:: MDAMP
  REAL(rkind):: WDAMP,FDAMP,rdamp_min,rdamp_max,zdamp_min,zdamp_max

  ! === collsion inhancement parameters ===

  INTEGER:: model_coll_enhance
  REAL(rkind):: &
       factor_coll_enhance, &
       xpos_coll_enhance,xwidth_coll_enhance, &
       ypos_coll_enhance,ywidth_coll_enhance

  ! === interpolation level ===

  INTEGER:: model_interpolation

  ! === medium parameters ===

  INTEGER:: nmmax,nkmax
  REAL(rkind):: epsdm(nmm),amudm(nmm),sigdm(nmm)
  INTEGER:: nmka(nkm)

!  ! === boundary type parameters ===

!  INTEGER:: nbmax,nbpmax(nbm)
!  integer(ikind):: KABDY(nbm)
!  real(rkind):: PHIBDY(nbm),RESBDY(nbm)
!  real(rkind):: PWRBDY(nbm),PHABDY(nbm)
!  real(rkind):: XGBDY(nbm),YGBDY(nbm),ZGBDY(nbm)
!  real(rkind):: XNBDY(3,nbm),YNBDY(3,nbm),ZNBDY(3,nbm)
!  real(rkind):: XPBDY(nbm),YPBDY(nbm),ZPBDY(nbm)
!  real(rkind):: SZBDY(2,nbm)
!  integer:: NDBDY(nbm),NMBDY(nbm)
!  integer:: NENBP(nbpm,nbm),NDNBP(nbpm,nbm)

  ! === output parameters ===

  INTEGER:: nprint,ngraph,ndrawd,ndraws,ndrawa,ndrawe,ndrawv

  ! === sort weight parameters ===

  REAL(rkind):: sort_weight_x,sort_weight_y

  ! === divider parameters ===

  REAL(rkind):: delr,delz,bdrmin,bdrmax,bdzmin,bdzmax

  ! === graphic paramters ===

  INTEGER:: ngxmax,ngymax,ngvmax
  REAL(rkind):: gfactor

  ! === division control parameter ===

  INTEGER:: nxzone_max,nyzone_max

  ! === computation control parameters ===

  REAL(rkind):: tolerance

  ! === debug controll parameters ===

  INTEGER:: idebuga(idebuga_max)

END MODULE wfcomm_parm

MODULE wfcomm

  USE wfcomm_parm

  IMPLICIT NONE
  
  !       --- for parallel computing ---
!  integer :: nrank,nprocs
  integer :: istart,iend

  INTEGER :: node_max_save=0
  INTEGER :: nelm_max_save=0
  INTEGER :: nseg_max_save=0
  INTEGER :: nelm_max_sort_save=0
  INTEGER :: nelm_max_solve_save=0
  INTEGER :: mtx_len_solve_save=0
  INTEGER :: nseg_max_field_save=0
  INTEGER :: node_max_field_save=0
  INTEGER :: ngxmax_save=0
  INTEGER :: ngymax_save=0
  INTEGER :: ngvmax_save=0
  
!      ----- graphics -----

!  integer :: NBPM
!  integer,parameter::  NMDM = 10
!  integer,parameter::  NCNM = 12+NMDM

!  integer,parameter::   NGM =   3
  
  integer,parameter::  NWDM = 12
  integer,parameter::  NCHM = 80

!       --- common variables ---

  complex(rkind):: CII

  REAL(rkind):: PSIA

!       /WFDIV/
  integer:: itype_mesh=0
  integer:: NRM,NZM,NZMH
         
!       /WFELM/
  integer(ikind):: node_max,nelm_max,node_bdy_max,nseg_bdy_max
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: xnode,ynode !(node_max)
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
  integer(ikind),dimension(:,:),ALLOCATABLE :: nseg_nside_nelm       !(3,nelm_max)
                                                ! node number of element
  integer(ikind),dimension(:,:),ALLOCATABLE :: KNELM       !(3,nelm_max)
                                                ! element number of adjascent
  integer(ikind),dimension(:,:),ALLOCATABLE :: node_nside_nelm! NSDELM     !(3,nelm_max)
                                                ! side number of element
  integer(ikind),dimension(:)  ,ALLOCATABLE :: NVNN        !(node_max)
                                                ! variable number of node E
        
!       /WFSID/
  integer:: nseg_max
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: LSID       !(nseg_max)
                                                ! length of side
  integer(ikind),dimension(:,:),ALLOCATABLE :: node_nseg ! NDSID  !(2,nseg_max)
                                                ! node number of side
  integer(ikind),dimension(:)  ,ALLOCATABLE :: INSID,NESID!(nseg_max) 
                                                ! side number in a element
                                                ! element number of a side
  integer(ikind),dimension(:)  ,ALLOCATABLE :: KASID,KBSID!(nseg_max) 
                                                ! if boundary
                                                !   KASID=1
                                                !   KBSID: boundary side number
  integer(ikind),dimension(:)  ,ALLOCATABLE :: NVNSD      !(nseg_max)
                                                ! variable number of side E
        
!       /WFSRT/
  real(rkind)   ,dimension(:),ALLOCATABLE :: SINDEX       !(nelm_max)
  integer(ikind),dimension(:),ALLOCATABLE :: IVELM,IWELM  !(nelm_max)  
  integer(ikind),dimension(:),ALLOCATABLE :: IDELM        !(nelm_max)
  real(rkind)   ,dimension(:),ALLOCATABLE :: SINDEX_MIN   !(nelm_max)
  real(rkind)   ,dimension(:),ALLOCATABLE :: SINDEX_MAX   !(nelm_max)
        
  INTEGER(ikind),DIMENSION(:),ALLOCATABLE:: NSDBS,NNDBS
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CEBSD,CEBND
        
!       /WFSLV/
  integer(ikind):: mtx_len,NNBMAX
  complex(rkind):: CM(6,6)
  complex(rkind),dimension(:)  ,ALLOCATABLE :: CSV    !(MLEN)
  complex(rkind),dimension(:,:),ALLOCATABLE :: CVTOT  !(6,nelm_max)
  complex(rkind),dimension(:)  ,ALLOCATABLE :: CQQ    !(MBND)

!       /WFAIF/
  real(rkind),dimension(3,3,3):: AIF3,AIE3
  real(rkind),dimension(3,3)  :: AIF2,AIE2
  real(rkind),dimension(3)    :: AIF1,AIE1
        
!       /WFFLD/
  complex(rkind),dimension(:)  ,ALLOCATABLE :: CESD   !(nseg_max)
  complex(rkind),dimension(:)  ,ALLOCATABLE :: CEND   !(NDMAX)
  complex(rkind),dimension(:,:),ALLOCATABLE :: CBF,CBP!(3,NNM)
  real(rkind):: EMAX(4)
  real(rkind):: ETMAX,PNMAX
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
  complex(rkind),dimension(NAM)    :: CIMP
  complex(rkind)                   :: CTIMP
  INTEGER:: npoint_max_nant(NAM)
  real(rkind),dimension(NJM,NAM)   :: RJ,ZJ
  real(rkind),dimension(NJM,NAM)   :: RJ0,ZJ0
        
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
        
!       /WFDBG/
  integer(ikind):: NDFILE

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


  ! === allocate node array ===

  SUBROUTINE wf_node_allocate
    IMPLICIT NONE
    
    IF(node_max.EQ.node_max_save) RETURN
    IF(node_max_save.NE.0) CALL wf_node_deallocate

    ALLOCATE(xnode(node_max),ynode(node_max))
    ALLOCATE(KANOD(node_max),KBNOD(node_max))
    ALLOCATE(NVNN(node_max))
    node_max_save = node_max

  END SUBROUTINE wf_node_allocate

  SUBROUTINE wf_node_deallocate
    IMPLICIT NONE

    DEALLOCATE(xnode,ynode,KANOD,KBNOD,NVNN)
    node_max_save=0
    
    RETURN
  END SUBROUTINE wf_node_deallocate

  ! === allocate nelm array ===

  SUBROUTINE wf_nelm_allocate
    IMPLICIT NONE

    IF(nelm_max.eq.nelm_max_save) RETURN
    IF(nelm_max_save.NE.0) CALL wf_nelm_deallocate

    ALLOCATE(SELM(nelm_max),KAELM(nelm_max))
    ALLOCATE(REMIN(nelm_max),ZEMIN(nelm_max))
    ALLOCATE(REMAX(nelm_max),ZEMAX(nelm_max))
    ALLOCATE(node_nside_nelm(3,nelm_max))
    ALLOCATE(nseg_nside_nelm(3,nelm_max))
    ALLOCATE(KNELM(3,nelm_max))
    nelm_max_save = nelm_max

  END SUBROUTINE wf_nelm_allocate

  SUBROUTINE wf_nelm_deallocate
    IMPLICIT NONE

    DEALLOCATE(SELM,KAELM,REMIN,ZEMIN,REMAX,ZEMAX,node_nside_nelm,KNELM)
    nelm_max_save = 0

  END SUBROUTINE wf_nelm_deallocate
  
  ! === allocate nseg array ===
  
  SUBROUTINE wf_nseg_allocate
    IMPLICIT NONE

    IF(nseg_max.EQ.nseg_max_save) RETURN
    !    IF(nseg_max_save.NE.0) CALL wf_nseg_deallocate
    IF(ALLOCATED(node_nseg)) CALL wf_nseg_deallocate

    ALLOCATE(node_nseg(2,nseg_max))
    ALLOCATE(KASID(nseg_max),KBSID(nseg_max),INSID(nseg_max))
    ALLOCATE(NESID(nseg_max),LSID(nseg_max),NVNSD(nseg_max))
    nseg_max_save = nseg_max
    
    RETURN
  END SUBROUTINE wf_nseg_allocate

  SUBROUTINE wf_nseg_deallocate
    IMPLICIT NONE

    DEALLOCATE(node_nseg,KASID,KBSID,INSID,NESID,LSID,NVNSD)
    nseg_max_save=0

    RETURN
  END SUBROUTINE wf_nseg_deallocate

  ! === allocate array for sort ===
  
  SUBROUTINE wf_sort_allocate
    IMPLICIT NONE
    
    IF(nelm_max.eq.nelm_max_sort_save) RETURN
    IF(nelm_max_sort_save.NE.0) CALL wf_sort_deallocate

    ALLOCATE(SINDEX(nelm_max),IVELM(nelm_max),IWELM(nelm_max),IDELM(nelm_max))
    ALLOCATE(SINDEX_MIN(nelm_max),SINDEX_MAX(nelm_max))
    nelm_max_sort_save = nelm_max
    
    RETURN
  END SUBROUTINE wf_sort_allocate
  
  SUBROUTINE wf_sort_deallocate
    IMPLICIT NONE

    DEALLOCATE(SINDEX,IVELM,IWELM,IDELM,SINDEX_MIN,SINDEX_MAX)
    nelm_max_sort_save = 0

    RETURN
  END SUBROUTINE wf_sort_deallocate

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

  ! === allocate array for solve ===
  
  SUBROUTINE wf_solve_allocate
    IMPLICIT NONE

    IF((mtx_len .EQ.mtx_len_solve_save).AND. &
       (nelm_max.EQ.nelm_max_solve_save)) RETURN
    IF(mtx_len_solve_save.NE.0) CALL wf_solve_deallocate

    ALLOCATE(CSV(mtx_len),CVTOT(6,nelm_max))
    mtx_len_solve_save=mtx_len
    nelm_max_solve_save=nelm_max
    
    RETURN
  END SUBROUTINE wf_solve_allocate

  SUBROUTINE wf_solve_deallocate
    IMPLICIT NONE

    DEALLOCATE(CSV,CVTOT)
    mtx_len_solve_save=0
    nelm_max_solve_save=0

    RETURN
  END SUBROUTINE wf_solve_deallocate

  ! === allocate array for field ===
  
  SUBROUTINE wf_field_allocate
    IMPLICIT NONE

    IF((nseg_max.eq.nseg_max_field_save).AND.&
       (node_max.eq.node_max_field_save)) RETURN
    IF(nseg_max_field_save.NE.0) CALL wf_field_deallocate

    ALLOCATE(CESD(nseg_max),CEND(node_max))
    ALLOCATE(CBF(3,node_max),CBP(3,node_max))
    nseg_max_field_save = nseg_max
    node_max_field_save = node_max

  END SUBROUTINE wf_field_allocate

  SUBROUTINE wf_field_deallocate
    IMPLICIT NONE

    DEALLOCATE(CESD,CEND,CBF,CBP)
    nseg_max_field_save = 0
    node_max_field_save = 0

  END SUBROUTINE wf_field_deallocate
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

  SUBROUTINE wf_win_allocate
    IMPLICIT NONE
    INTEGER,SAVE :: NGXMAX_save,NGYMAX_save,NGVMAX_save

    IF((ngxmax.EQ.ngxmax_save).AND. &
       (ngymax.EQ.ngymax_save).AND. &
       (ngvmax.EQ.ngvmax_save)) RETURN
    IF(ngxmax_save.NE.0) CALL wf_win_deallocate

    ALLOCATE(gz(ngxmax,ngymax),gz_temp(ngxmax,ngymax),g2x(ngxmax),g2y(ngymax))
    ALLOCATE(iegz(ngxmax,ngymax))
    ALLOCATE(gv(ngvmax,3),gx(ngvmax))
    ngxmax_save = ngxmax
    ngymax_save = ngymax
    ngvmax_save = ngvmax
    
    RETURN
  END SUBROUTINE wf_win_allocate

  SUBROUTINE wf_win_deallocate
    IMPLICIT none

    DEALLOCATE(GZ,GZ_temp,G2X,G2Y,GV,GX,IEGZ)
    ngxmax_save = 0
    ngymax_save = 0
    ngvmax_save = 0

    RETURN
  END SUBROUTINE wf_win_deallocate

END MODULE wfcomm
