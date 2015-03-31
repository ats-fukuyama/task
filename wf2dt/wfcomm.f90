!     $Id: wfcomm.f90,v 1.23 2012/03/05 06:29:02 maruyama Exp $

module wfcomm

  use bpsd
  use plcomm
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
  integer,parameter::   NJM = 800  ! Maximum number of antenna current
  
  integer,parameter::  NGXM = 301  ! resolution in x direction
  integer,parameter::  NGYM = 301  ! resolution in y direction
  integer,parameter::  NGVM = 301  ! resolution used by 1D plot
  integer,parameter::   NGM =   3
  
  integer,parameter::  NWDM = 12
  integer,parameter::  NCHM = 80

!       --- common variables ---
  complex(rkind):: CII

!       /WFPRM/
  real(rkind):: RF
  integer(ikind):: NAMAX,NPH
  real(rkind):: PPN0,PTN0,PIN
  integer(ikind):: NPRINT,NDRAWD,NDRAWA,NDRAWE,NGRAPH,NDRAWV
  integer(ikind):: MODELI,MODELB
  integer(ikind):: MODELD,MODELP
  real(rkind),DIMENSION(3):: r_corner,z_corner
  real(rkind),DIMENSION(3):: br_corner,bz_corner,bt_corner
  real(rkind),DIMENSION(3,NSM):: pn_corner,ptpr_corner,ptpp_corner
  real(rkind):: R1WG,Z1WG,R2WG,Z2WG,PH1WG,PH2WG,AMPWG,ANGWG

!       /WFPRK/
  character(len=32) :: KFNAME,KFNAMA,KFNAMF,KFNAMN
  character(len=32) :: KFNAMB

!       /WFPRD/
  real(rkind):: BRMIN,BRMAX,BZMIN,BZMAX
  real(rkind):: DELR,DELZ
  real(rkind):: RD,THETJ1,THETJ2
  integer(ikind):: NJMAX
  real(rkind),dimension(NAM):: AJ,APH,APOS,AWD
  real(rkind):: WDUMP

!       /WFDIV/
  integer(ikind):: iddiv
         
!       /WFELM/
  integer(ikind):: NNMAX,NEMAX
  real(rkind)   ,dimension(:)  ,pointer :: RNODE,ZNODE !(NNMAX)
  integer(ikind),dimension(:)  ,pointer :: KANOD,KBNOD !(NNMAX)
  real(rkind)   ,dimension(:)  ,pointer :: SELM        !(NEMAX)
  integer(ikind),dimension(:)  ,pointer :: KAELM,NBELM !(NEMAX)
  real(rkind)   ,dimension(:)  ,pointer :: REMIN,ZEMIN !(NEMAX)
  real(rkind)   ,dimension(:)  ,pointer :: REMAX,ZEMAX !(NEMAX)
  integer(ikind),dimension(:,:),pointer :: NDELM       !(5,NEMAX)
  integer(ikind),dimension(:,:),pointer :: KNELM       !(3,NEMAX)
  integer(ikind),dimension(:)  ,pointer :: NVNN        !(NNMAX)
        
!       /WFSID/
  integer(ikind):: NSDMAX
  real(ikind)   ,dimension(:)  ,pointer :: LSID       !(NSDMAX)
  integer(ikind),dimension(:,:),pointer :: NDSID      !(2,NSDMAX)
  integer(ikind),dimension(:)  ,pointer :: INSID,NESID!(NSDMAX) 
  integer(ikind),dimension(:)  ,pointer :: KASID,KBSID!(NSDMAX) 
  integer(ikind),dimension(:,:),pointer :: NSDELM     !(3,NEM)
  integer(ikind),dimension(:)  ,pointer :: NVNSD      !(NSDMAX)
        
!       /WFSRT/
  real(rkind)   ,dimension(:),pointer :: SINDEX       !(NEMAX)
  integer(ikind),dimension(:),pointer :: IVELM,IWELM  !(NEMAX)  
  integer(ikind),dimension(:),pointer :: IDELM        !(NEMAX)
  real(rkind)   ,dimension(:),pointer :: SINDEX_MIN   !(NEMAX)
  real(rkind)   ,dimension(:),pointer :: SINDEX_MAX   !(NEMAX)
        
!       /WFMED/
  integer(ikind):: NMMAX,NKMAX
  real(rkind)   ,dimension(:),pointer :: EPSDM,AMUDM,SIGDM !(NMMAX)
  integer(ikind),dimension(:),pointer :: NMKA              !(NKMAX)
        
!       /WFBDY/
  integer(ikind):: NBMAX
  integer(ikind),dimension(:)  ,pointer :: KABDY             !(NBMAX)
  real(rkind)   ,dimension(:)  ,pointer :: PHIBDY,RESBDY     !(NBMAX)
  real(rkind)   ,dimension(:)  ,pointer :: PWRBDY,PHABDY     !(NBMAX)
  real(rkind)   ,dimension(:)  ,pointer :: XGBDY,YGBDY,ZGBDY !(NBMAX)
  real(rkind)   ,dimension(:,:),pointer :: XNBDY,YNBDY,ZNBDY !(3,NBMAX)
  real(rkind)   ,dimension(:)  ,pointer :: XPBDY,YPBDY,ZPBDY !(NBMAX)
  real(rkind)   ,dimension(:,:),pointer :: SZBDY             !(2,NBMAX)
  integer(ikind),dimension(:)  ,pointer :: NDBDY,NMBDY,NBPMAX!(NBMAX)
  integer(ikind),dimension(:,:),pointer :: NENBP,NDNBP       !(NBPM,NBMAX)

  INTEGER(ikind),DIMENSION(:),ALLOCATABLE:: NSDBS,NNDBS
  COMPLEX(rkind),DIMENSION(:),ALLOCATABLE:: CEBSD,CEBND
        
!       /WFSLV/
  integer(ikind):: MLEN,NNBMAX
  integer(ikind):: NMDMAX
  complex(rkind):: CM(6,6)
  complex(rkind),dimension(:)  ,pointer :: CSV    !(MLEN)
  complex(rkind),dimension(:,:),pointer :: CVTOT  !(6,NEMAX)
  complex(rkind),dimension(:)  ,pointer :: CQQ    !(MBND)

!       /WFAIF/
  real(rkind),dimension(3,3,3):: AIF3
  real(rkind),dimension(3,3)  :: AIF2
  real(rkind),dimension(3)    :: AIF1
        
!       /WFFLD/
  complex(rkind),dimension(:)  ,pointer :: CESD   !(NSDMAX)
  complex(rkind),dimension(:)  ,pointer :: CEND   !(NDMAX)
  complex(rkind),dimension(:,:),pointer :: CBF,CBP!(3,NNM)
  real(rkind)   ,dimension(:)  ,pointer :: EMAX   !(4)
  real(rkind)::ETMAX,PNMAX
  complex(rkind),dimension(:,:),pointer :: CRFL!(NMDM,NBM)
       
!       /WFPWR/
  real(rkind),dimension(NSM):: PABST
!  real(rkind),dimension(:)  ,pointer :: PABSS !(NSM)
!  real(rkind),dimension(:)  ,pointer :: PABSK !(NKMAX)
!  real(rkind),dimension(:)  ,pointer :: PABSTN!(NNMAX)
!  real(rkind),dimension(:,:),pointer :: PABSSN!(NNMAX,NSM)
!  real(rkind),dimension(:,:),pointer :: PABSKN!(NNMAX,NKMAX)
!  real(rkind),dimension(:,:),pointer :: PFV   !(NNMAX,3)
        
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
!  integer(ikind),dimension(:),pointer :: IDND !(NNMAX)
!  real(rkind)   ,dimension(:),pointer :: EX1WG,EY1WG,EZ1WG!(NBMAX)
!  real(rkind)   ,dimension(:),pointer :: EX2WG,EY2WG,EZ2WG!(NBMAX)
!  real(rkind)   ,dimension(:),pointer :: PWRWG,PHAWG!(NBMAX)
!  integer(ikind),dimension(:),pointer :: IDKA !(NKMAX)
!  integer(ikind),dimension(:),pointer :: IDMAT!(NMMAX)
!  integer(ikind),dimension(:),pointer :: IDBDY!(NBMAX)
!  CHARACTER     ,dimension(:),pointer :: KDKA*25 !(NKMAX)
!  CHARACTER     ,dimension(:),pointer :: KDMAT*25!(NMMAX)
!  CHARACTER     ,dimension(:),pointer :: KDBDY*25!(NBMAX)
!  integer(ikind):: IDNMIN,IDNMAX,IDEMIN,IDEMAX
        
!       /WFWIN/
  real(rkind):: RNDMIN,RNDMAX,ZNDMIN,ZNDMAX
  real(rkind):: LNDMIN,LNDMAX
  integer(ikind):: NFOPEN
  integer(ikind):: NWXMAX
  real(4),dimension(:,:),pointer :: GZ,GZ_temp !(NGXM,NGYM)
  real(4),dimension(:)  ,pointer :: G2X!(NGXM)
  real(4),dimension(:)  ,pointer :: G2Y!(NGYM)
  real(4),dimension(:,:),pointer :: GV !(NGVM,NGM)
  real(4),dimension(:)  ,pointer :: GX !(NGVM)
  integer(ikind):: NGXMAX,NGYMAX,NGVMAX
        
!       /WFDBG/
  integer(ikind):: NDFILE

! --- allocation flag ---
  integer(rkind):: divinit,elminit,sidinit,srtinit,medinit
  integer(rkind):: srfinit,slvinit,fldinit,pwrinit,nasinit,wininit

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
       allocate(SELM(NEMAX),KAELM(NEMAX),NBELM(NEMAX))
       allocate(REMIN(NEMAX),ZEMIN(NEMAX))
       allocate(REMAX(NEMAX),ZEMAX(NEMAX))
       allocate(NDELM(3,NEMAX),KNELM(3,NEMAX))
       elminit = 2
       NEMAX_save = NEMAX
    end if

    return
  end subroutine wfelm_allocate
!----
  subroutine wfelm_deallocate
    implicit none

    deallocate(RNODE,ZNODE,KANOD,KBNOD,SELM,KAELM,NBELM)
    deallocate(REMIN,ZEMIN,REMAX,ZEMAX,NDELM,KNELM,NVNN)

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

    allocate(NDSID(2,NSDMAX),KASID(NSDMAX),KBSID(NSDMAX),INSID(NSDMAX))
    allocate(NESID(NSDMAX),NSDELM(3,NEMAX),LSID(NSDMAX),NVNSD(NSDMAX))

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

    IF(ASSOCIATED(epsdm)) DEALLOCATE(epsdm)
    IF(ASSOCIATED(amudm)) DEALLOCATE(amudm)
    IF(ASSOCIATED(sigdm)) DEALLOCATE(sigdm)
    IF(ASSOCIATED(nmka)) DEALLOCATE(nmka)

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
    integer,save :: NGXM_save,NGYM_save,NGVM_save

    if(wininit.eq.1) then
       if((NGXM.eq.NGXM_save).and.&
         &(NGYM.eq.NGYM_save).and.&
         &(NGVM.eq.NGVM_save)) then
          return
       else
          call wfwin_deallocate
       end if
    end if

    allocate(GZ(NGXM,NGYM),GZ_temp(NGXM,NGYM),G2X(NGXM),G2Y(NGYM))
    allocate(GV(NGVM,NGM),GX(NGVM))

    NGXM_save = NGXM
    NGYM_save = NGYM
    NGVM_save = NGVM
    wininit = 1
    
    return
  end subroutine wfwin_allocate
!-----
  subroutine wfwin_deallocate
    implicit none

    deallocate(GZ,GZ_temp,G2X,G2Y,GV,GX)

    return
  end subroutine wfwin_deallocate

end module wfcomm
