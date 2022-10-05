!     $Id$

module wfcomm

  use bpsd
  use plcomm
  use commpi
  implicit none

!  public
!       --- for parallel computing ---
!  integer(ikind) :: nrank,nsize
  integer(ikind) :: istart,iend

!       --- input parameters ---

!  integer(ikind),parameter:: NXM = 201
!  integer(ikind),parameter:: NYM = 201
!  integer(ikind),parameter:: NZM = 101

  integer(ikind) :: NXM,NYM,NZM
  integer(ikind) :: NYMH

!      ----- graphics -----

  integer(ikind) :: NBPM
  integer(ikind),parameter::  NMDM = 10
  integer(ikind),parameter::  NCNM = 12+NMDM

  integer(ikind),parameter::   NBM = 100
  integer(ikind),parameter::   NKM = 100   
  integer(ikind),parameter::   NMM = 100   
  integer(ikind),parameter::   NAM =  25

  integer(ikind),parameter::   NJM = 800
  
  integer(ikind),parameter::  NGXM = 101
  integer(ikind),parameter::  NGYM = 101
  integer(ikind),parameter::  NGVM = 101
  integer(ikind),parameter::   NGM =   3
  
  integer(ikind),parameter::  NWDM = 12
  integer(ikind),parameter::  NCHM = 80

!       --- common variables ---
  complex(rkind):: CII

!       /WFPRM/
  real(rkind):: RF
  integer(ikind):: NAMAX
!  real(rkind),dimension(NSM):: PZCL
  real(rkind):: ZPMIN,ZPMAX,ZPLEN
!  real(rkind):: PPN0,PTN0,PIN
  real(rkind):: PIN
  integer(ikind):: NPRINT,NDRAWD,NDRAWA,NGRAPH
  integer(ikind):: MODELI
  integer(ikind):: MODELD,MODELP,MODELS,MODELX,MODELA
  real(rkind):: POSRES,POSABS,EPSABS,DLTABS

!       /WFPRK/
  character(len=32) :: KFNAME,KFNAMA,KFNAMF,KFNAMN
  character(len=32) :: KFNAMB

!       /WFPRD/
  real(rkind):: BXMIN,BXMAX,BYMIN,BYMAX,RBAX
  real(rkind):: BZMIN,BZMIDL,BZMIDH,BZMAX,DELX,DELY,DELZ,DelZM
  real(rkind):: RD,THETJ1,THETJ2
  integer(ikind):: NJMAX
  real(rkind),dimension(NAM):: AJ,APH,APOS,AWD
  real(rkind):: THETS1,THETS2,RD1,RD2,ZANT,ZWALL
 
!       /WFDIV/
  real(rkind)   ,dimension(:)    ,ALLOCATABLE:: XL,XR!(NYM)
  integer(ikind),dimension(:)    ,ALLOCATABLE:: NXA  !(NYM)
  integer(ikind):: NYMAX
  integer(ikind):: IDDIV
  real(rkind)   ,dimension(:,:)  ,ALLOCATABLE:: XNDA !(NXM,NYM)
  real(rkind)   ,dimension(:,:)  ,ALLOCATABLE:: YNDA !(NXM,NYM)
  integer(ikind),dimension(:,:,:),ALLOCATABLE:: NDA  !(NXM,NYM,NZM)
  real(rkind):: XF,YF
  real(rkind)   ,dimension(:)    ,ALLOCATABLE:: ZNDA !(NZM)
  integer(ikind):: NZMAX
        
!       /WFELM/
  integer(ikind):: NNMAX,NEMAX
  real(rkind)   :: VTOT
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: XND,YND,ZND,VNOD !(NNMAX)
  integer(ikind),dimension(:)  ,ALLOCATABLE :: KANOD            !(NNMAX)
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: VELM             !(NEMAX)
  integer(ikind),dimension(:)  ,ALLOCATABLE :: KAELM,NBELM      !(NEMAX)
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: XEMIN,YEMIN,ZEMIN!(NEMAX)
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: XEMAX,YEMAX,ZEMAX!(NEMAX)
  integer(ikind),dimension(:,:),ALLOCATABLE :: NDELM            !(5,NEMAX)
  integer(ikind),dimension(:,:),ALLOCATABLE :: KNELM            !(4,NEMAX)
        
!       /WFSID/
  integer(ikind):: NSDMAX
  integer(ikind),dimension(:,:),ALLOCATABLE :: NDSID  !(2,NSDM)
  integer(ikind),dimension(:)  ,ALLOCATABLE :: KASID  !(NSDM) 
  integer(ikind),dimension(:,:),ALLOCATABLE :: NSDELM !(7,NEM)
  integer(ikind),dimension(:,:),ALLOCATABLE :: NSDSRF !(3,NSFM)
        
!       /WFSRT/
  real(rkind)   ,dimension(:),ALLOCATABLE :: SINDEX               !(NEMAX)
  integer(ikind),dimension(:),ALLOCATABLE :: IVELM,IWELM          !(NEMAX)  
  integer(ikind),dimension(:),ALLOCATABLE :: IDELM!(NEMAX)
  real(rkind)   ,dimension(:),ALLOCATABLE :: SINDEX_MIN,SINDEX_MAX!(NEMAX)
        
!       /WFMED/
  integer(ikind):: NMMAX,NKMAX
  real(rkind)   ,dimension(:),ALLOCATABLE :: EPSDM,AMUDM,SIGDM !(NMMAX)
  integer(ikind),dimension(:),ALLOCATABLE :: NMKA              !(NKMAX)
        
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
        
!       /WFSRF/
  integer(ikind):: NSFMAX
  integer(ikind),dimension(:)  ,ALLOCATABLE :: INSRF,NESRF !(NSFMAX)
  integer(ikind),dimension(:,:),ALLOCATABLE :: NDSRF,KNSRF !(3,NSFMAX)
  integer(ikind),dimension(:,:),ALLOCATABLE :: NSFELM      !(4,NEMAX)
        
!       /WFSLV/
  integer(ikind):: MBND,MLEN,NNBMAX
  integer(ikind):: NMDMAX
  complex(rkind):: CM(NMDM,NMDM,7,7),CV(NMDM,7)
  integer(ikind),dimension(:,:),ALLOCATABLE :: ISDELM        !(7,NEMAX)
  integer(ikind),dimension(:)  ,ALLOCATABLE :: IMLEN,INLEN   !(NSDMAX)
!  complex(rkind),dimension(:)  ,ALLOCATABLE :: CRV           !(MLENM)
  complex(rkind),dimension(:)  ,ALLOCATABLE :: CSV           !(MLENM)
  complex(rkind),dimension(:,:),ALLOCATABLE :: CVTOT         !(7,NEMAX)
  integer(ikind),dimension(:)  ,ALLOCATABLE :: LDEST,NODEK   !(NCNM)
  integer(ikind),dimension(:)  ,ALLOCATABLE :: NFLG          !(NSDMAX)
!  complex(rkind),dimension(:,:),ALLOCATABLE :: CEQ           !(MBND,MBND)
  complex(rkind),dimension(:)  ,ALLOCATABLE :: CQQ           !(MBND)
  integer(ikind),dimension(:)  ,ALLOCATABLE :: LPIV          !(MBND)
  integer(ikind),dimension(:)  ,ALLOCATABLE :: LHED          !(MLEN)
  integer(ikind),dimension(:)  ,ALLOCATABLE :: MLCO,MCOL,MPOS!(MLEN)
  integer(ikind),dimension(:)  ,ALLOCATABLE :: NSDNV         !(MLEN)
        
!       /WFAIF/
  real(rkind),dimension(3,3,3):: AIF3
  real(rkind),dimension(3,3)  :: AIF2
  real(rkind),dimension(3)    :: AIF1
  real(rkind),dimension(4,4,4):: AIG3
  real(rkind),dimension(4,4)  :: AIG2
  real(rkind),dimension(4)    :: AIG1
        
!       /WFFLD/
  complex(rkind),dimension(:)  ,ALLOCATABLE :: CESD   !(NSDM)
  complex(rkind),dimension(:,:),ALLOCATABLE :: CEF,CEP!(3,NNM)
  complex(rkind),dimension(:,:),ALLOCATABLE :: CBF,CBP!(3,NNM)
  real(rkind)   ,dimension(:)  ,ALLOCATABLE :: EMAX   !(4)
  real(rkind)::ETMAX,PNMAX
  complex(rkind),dimension(:,:),ALLOCATABLE :: CRFL!(NMDM,NBM)
       
!       /WFPWR/
  real(rkind):: PABST
  real(rkind),dimension(:)  ,ALLOCATABLE :: PABSS !(NSM)
  real(rkind),dimension(:)  ,ALLOCATABLE :: PABSK !(NKMAX)
  real(rkind),dimension(:)  ,ALLOCATABLE :: PABSTN!(NNMAX)
  real(rkind),dimension(:,:),ALLOCATABLE :: PABSSN!(NNMAX,NSM)
  real(rkind),dimension(:,:),ALLOCATABLE :: PABSKN!(NNMAX,NKMAX)
  real(rkind),dimension(:,:),ALLOCATABLE :: PFV   !(NNMAX,3)
        
!       /WFANT/
  complex(rkind),dimension(NAM)  :: CIMP
  complex(rkind):: CTIMP
  real(rkind),dimension(NJM,NAM)   :: XJ,YJ,ZJ
  integer(ikind),dimension(NJM,NAM):: JELMT
  integer(ikind),dimension(NAM)    :: JNUM
  real(rkind),dimension(NJM,NAM)   :: XJ0,YJ0,ZJ0
  integer(ikind),dimension(NAM)    :: JNUM0
  complex(rkind),dimension(NJM,NAM):: CEIMP
        
!       /WFNAS/
  real(rkind):: FACT_LEN
  integer(ikind),dimension(:),ALLOCATABLE :: IDND !(NNMAX)
  real(rkind)   ,dimension(:),ALLOCATABLE :: EX1WG,EY1WG,EZ1WG!(NBMAX)
  real(rkind)   ,dimension(:),ALLOCATABLE :: EX2WG,EY2WG,EZ2WG!(NBMAX)
  real(rkind)   ,dimension(:),ALLOCATABLE :: PWRWG,PHAWG!(NBMAX)
  integer(ikind),dimension(:),ALLOCATABLE :: IDKA !(NKMAX)
  integer(ikind),dimension(:),ALLOCATABLE :: IDMAT!(NMMAX)
  integer(ikind),dimension(:),ALLOCATABLE :: IDBDY!(NBMAX)
  CHARACTER     ,dimension(:),ALLOCATABLE :: KDKA*25 !(NKMAX)
  CHARACTER     ,dimension(:),ALLOCATABLE :: KDMAT*25!(NMMAX)
  CHARACTER     ,dimension(:),ALLOCATABLE :: KDBDY*25!(NBMAX)
  integer(ikind):: IDNMIN,IDNMAX,IDEMIN,IDEMAX
        
!       /WFWIN/
  real(rkind):: XNDMIN,XNDMAX,YNDMIN,YNDMAX,ZNDMIN,ZNDMAX
  real(rkind):: RNDMIN,RNDMAX
  integer(ikind):: NFOPEN
  integer(ikind):: NWXMAX
  real(4),dimension(:,:),ALLOCATABLE :: GZ !(NGXM,NGYM)
  real(4),dimension(:)  ,ALLOCATABLE :: G2X!(NGXM)
  real(4),dimension(:)  ,ALLOCATABLE :: G2Y!(NGYM)
  real(4),dimension(:,:),ALLOCATABLE :: GV !(NGVM,NGM)
  real(4),dimension(:)  ,ALLOCATABLE :: GX !(NGVM)
  integer(ikind):: NGXMAX,NGYMAX,NGVMAX
  CHARACTER,dimension(0:9) :: KGINX*80,KGINV*80

! ----- Add. By YOKOYAMA Mar./05/2013 ----
!     (wfgout.f) Magnetic Field and Plasma density Profile
   real(rkind),dimension(:),ALLOCATABLE:: YBABS,YDEN,YDENI,YPSI
   real(rkind):: DUMMY1
   real(rkind):: DUMMY2(NSM),DUMMY3(NSM),DUMMY4(NSM)
!
!     (wfdiv.f) Elements
!     /YAMA03/YAMA04/
   real(rkind):: RBIN,RBOUT,DXIN,DYIN,DXOUT,DYOUT
   integer(ikind):: INOD,NYTEMP,NZTEMP
!
!     (wfgsub.f) FRATIO
!     /YAMA05/
   real(rkind):: FRATIO
!
!     (wffile.f) B-FIELD LINE
!    /YAMA06/07/08/
    integer(ikind):: NGFLIN
    real(rkind),dimension(:),ALLOCATABLE:: FLZ,FLX,FLY,XYR
    real(rkind):: FACTC, FACTA
!
!     (wfwave.f) Loaging Impedance
!    /YAMA09/
    complex(rkind),dimension(NAM):: CLOAD
!
!     (wfwave.f) E-r, E-theta
!    /YAMA10/
	real(rkind),dimension(:),ALLOCATABLE::ANGLE
    complex(rkind),dimension(:,:),ALLOCATABLE:: CERT,CBRT
!
!     (mbant.f) Rotarion Angle [deg.]
!    /YAMA11/
     real(rkind):: RTDEG
!
!     (wfwave.f) COLLISION FREQ.
!    /YAMA12/
      real(rkind),dimension(:,:),ALLOCATABLE:: RZCO
! ----- Mar./05/2013 -----
!        
!       /WFDBG/
  integer(ikind):: NDFILE

!------
contains
  subroutine wfdiv_initialize

    use libmpi
    use libmtx
    implicit none
    integer(ikind) :: idata(3)

    if (nrank.eq.0) then
       WRITE(6,*) "## DIV: NXM, NYM, NZM ?"
       READ(*,*)           NXM, NYM, NZM
       idata(1)=NXM
       idata(2)=NYM
       idata(3)=NZM
    end if

    call mtx_barrier
    call mtx_broadcast_integer(idata,3)
    NXM=idata(1)
    NYM=idata(2)
    NZM=idata(3)
    
    NYMH = NYM/2

    return
  end subroutine wfdiv_initialize
! ----------------------------------  
  subroutine wfdiv_allocate

    use libmtx
    implicit none
    integer,save:: NXM_save=0,NYM_save=0,NZM_save=0
    integer,save:: div_init=0
    
    if((NXM.eq.NXM_save).and.&
       (NYM.eq.NYM_save).and.&
       (NZM.eq.NZM_save)) then
       return
    else if(div_init.ne.0) then
       call wfdiv_deallocate
    end if
 
    allocate(XL(NYM),XR(NYM),NXA(NYM))
    allocate(XNDA(NXM,NYM),YNDA(NXM,NYM))
    allocate(NDA(NXM,NYM,NZM))
    allocate(ZNDA(NZM))
    
    NXM_save = NXM
    NYM_save = NYM
    NZM_save = NZM
    div_init = 1

    return
  end subroutine wfdiv_allocate
!----
  subroutine wfdiv_deallocate
    implicit none

    deallocate(XL,XR,NXA,XNDA,YNDA,ZNDA,NDA)
    
    return
  end subroutine wfdiv_deallocate
!----
  subroutine wfelm_allocate
    implicit none
    integer,save :: NNMAX_save,NEMAX_save
    integer,save :: elminit = 0

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

       allocate(XND(NNMAX),YND(NNMAX),ZND(NNMAX),VNOD(NNMAX),KANOD(NNMAX))
       elminit = 1
       NNMAX_save = NNMAX

    else if(elminit.eq.1) then
       allocate(VELM(NEMAX),KAELM(NEMAX),NBELM(NEMAX))
       allocate(XEMIN(NEMAX),YEMIN(NEMAX),ZEMIN(NEMAX))
       allocate(XEMAX(NEMAX),YEMAX(NEMAX),ZEMAX(NEMAX))
       allocate(NDELM(5,NEMAX),KNELM(4,NEMAX))
       elminit = 2
       NEMAX_save = NEMAX
    end if

    return
  end subroutine wfelm_allocate
!----
  subroutine wfelm_deallocate
    implicit none

    IF(.NOT.ALLOCATED(XND)) RETURN
    deallocate(XND,YND,ZND,VNOD,KANOD,VELM,KAELM,NBELM)
    deallocate(XEMIN,YEMIN,ZEMIN,XEMAX,YEMAX,ZEMAX,NDELM,KNELM)

    return
  end subroutine wfelm_deallocate
!----
  subroutine wfsid_allocate
    implicit none
    integer,save :: NSDMAX_save,NSFMAX_save,NEMAX_save
    integer,save :: sidinit = 0

    if(sidinit.eq.1) then
       if((NSDMAX.eq.NSDMAX_save).and.&
         &(NSFMAX.eq.NSFMAX_save).and.&
         &( NEMAX.eq. NEMAX_save)) then
          return
       else
          call wfsid_deallocate
       end if
    end if

    allocate(NDSID(2,NSDMAX),KASID(NSDMAX))
    allocate(NSDSRF(3,NSFMAX),NSDELM(7,NEMAX))

    NSDMAX_save = NSDMAX
    NSFMAX_save = NSFMAX
    sidinit = 1
    
    return
  end subroutine wfsid_allocate
!-----
  subroutine wfsid_deallocate
    implicit none

    IF(.NOT.ALLOCATED(NDSID)) RETURN
    deallocate(NDSID,KASID,NSDELM,NSDSRF)

    return
  end subroutine wfsid_deallocate
!-----
  subroutine wfsrt_allocate
    implicit none
    integer,save :: NEMAX_save
    integer,save :: srtinit = 0
    
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

    IF(.NOT.ALLOCATED(SINDEX)) RETURN
    deallocate(SINDEX,IVELM,IWELM,IDELM,SINDEX_MIN,SINDEX_MAX)

    return
  end subroutine wfsrt_deallocate
!-----
  subroutine wfmed_allocate
    implicit none
    integer,save :: NMMAX_save,NKMAX_save
    integer,save :: medinit = 0

    if(medinit.eq.0) then
       NMMAX = NMM
       NKMAX = NKM
    end if

    if(medinit.eq.1) then
       if((NMMAX.eq.NMMAX_save).and.&
         &(NKMAX.eq.NKMAX_save)) then
          return
       else
          call wfmed_deallocate
       end if
    end if
    
    allocate(EPSDM(NMMAX),AMUDM(NMMAX),SIGDM(NMMAX))
    allocate(NMKA(NKMAX))

    NMMAX_save = NMMAX
    NKMAX_save = NKMAX
    medinit = 1
    return
  end subroutine wfmed_allocate
!-----
  subroutine wfmed_deallocate
    implicit none
    
    IF(.NOT.ALLOCATED(EPSDM)) RETURN
    deallocate(EPSDM,AMUDM,SIGDM,NMKA)

    return
  end subroutine wfmed_deallocate
!-----
  subroutine wfbdy_allocate
    implicit none
    integer,save :: NBMAX_save
    integer,save :: bdyinit = 0

    if(bdyinit.eq.0) then
       NBMAX = NBM
       NBPM  = NEMAX
    end if

    if(bdyinit.eq.1) then
       if((NBMAX.eq.NBMAX_save)) then
          return
       else
          call wfbdy_deallocate
       end if
    end if
    
    allocate(KABDY(NBMAX),PHIBDY(NBMAX),RESBDY(NBMAX),PWRBDY(NBMAX))
    allocate(PHABDY(NBMAX),XGBDY(NBMAX),YGBDY(NBMAX),ZGBDY(NBMAX))
    allocate(XNBDY(3,NBMAX),YNBDY(3,NBMAX),ZNBDY(3,NBMAX))
    allocate(XPBDY(NBMAX),YPBDY(NBMAX),ZPBDY(NBMAX))
    allocate(SZBDY(2,NBMAX),NDBDY(NBMAX),NMBDY(NBMAX),NBPMAX(NBMAX))
    allocate(NENBP(NBPM,NBMAX),NDNBP(NBPM,NBMAX))
  
    NBMAX_save = NBMAX
    bdyinit = 1

    return
  end subroutine wfbdy_allocate
!-----
  subroutine wfbdy_deallocate
    implicit none

    IF(.NOT.ALLOCATED(KABDY)) RETURN
    deallocate(KABDY,PHIBDY,RESBDY,PWRBDY,PHABDY,XGBDY,YGBDY,ZGBDY)
    deallocate(XNBDY,YNBDY,ZNBDY,XPBDY,YPBDY,ZPBDY,SZBDY)
    deallocate(NDBDY,NMBDY,NBPMAX,NENBP,NDNBP)

    return
  end subroutine wfbdy_deallocate
!-----
  subroutine wfsrf_allocate
    implicit none
    integer,save :: NSFMAX_save,NEMAX_save
    integer,save :: srfinit = 0

    if(srfinit.eq.1) then
       if((NSFMAX.eq.NSFMAX_save).and.&
          ( NEMAX.eq. NEMAX_save)) then
          return
       else
          call wfsrf_deallocate
       end if
    end if
    
    allocate(NSFELM(4,NEMAX))
    allocate(INSRF(NSFMAX),NESRF(NSFMAX))
    allocate(NDSRF(3,NSFMAX),KNSRF(3,NSFMAX))

    NSFMAX_save = NSFMAX
     NEMAX_save = NEMAX
    srfinit = 1
    
    return
  end subroutine wfsrf_allocate
!-----
  subroutine wfsrf_deallocate
    implicit none

    IF(.NOT.ALLOCATED(INSRF)) RETURN
    deallocate(INSRF,NESRF,NDSRF,KNSRF,NSFELM)

    return
  end subroutine wfsrf_deallocate
!-----
  subroutine wfslv_allocate
    implicit none
    integer,save :: NSDMAX_save,NEMAX_save,MLEN_save,MBND_save
    integer,save :: slvinit = 0

    if(slvinit.eq.1) then
       if((NSDMAX.eq.NSDMAX_save).and.&
         &(NEMAX .eq. NEMAX_save).and.&
         &(MLEN  .eq.  MLEN_save).and.&
         &(MBND  .eq.  MBND_save)) then
          return
       else
          call wfslv_deallocate
       end if
    end if

    allocate(LDEST(NCNM),NODEK(NCNM))
    allocate(IMLEN(NSDMAX+1),INLEN(NSDMAX))
    allocate(NFLG(NSDMAX),ISDELM(7,NEMAX))
    allocate(LHED(MLEN),MLCO(MLEN),MCOL(MLEN),MPOS(MLEN),NSDNV(MLEN))
    allocate(CSV(MLEN),CVTOT(7,NEMAX))
    slvinit = 1
    NSDMAX_save = NSDMAX
    NEMAX_save = NEMAX
    MLEN_save = MLEN
!       allocate(CQQ(MLEN))

    return
  end subroutine wfslv_allocate
!-----
  subroutine wfslv_deallocate
    implicit none

    !deallocate(CRV,CEQ)
    IF(.NOT.ALLOCATED(ISDELM)) RETURN
    deallocate(ISDELM,IMLEN,INLEN,CSV,CVTOT)
    deallocate(LDEST,NODEK,NFLG)
!    deallocate(CQQ)
    deallocate(LHED,MLCO,MCOL,MPOS,NSDNV)

    return
  end subroutine wfslv_deallocate
!-----
  subroutine wffld_allocate
    implicit none
    integer,save :: NSDMAX_save,NNMAX_save,NMDMAX_save,NBMAX_save
    integer,save :: fldinit = 0

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
    
    allocate(CESD(NSDMAX),CEF(3,NNMAX),CEP(3,NNMAX))
    allocate(CBF(3,NNMAX),CBP(3,NNMAX),EMAX(4),CRFL(NMDMAX,NBMAX))
    allocate(ANGLE(NNMAX))

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

    IF(.NOT.ALLOCATED(CESD)) RETURN
    deallocate(CESD,CEF,CEP,CBF,CBP,EMAX,CRFL)

    return
  end subroutine wffld_deallocate
!-----
  subroutine wfpwr_allocate
    implicit none
    integer,save :: NKMAX_save,NNMAX_save
    integer,save :: pwrinit = 0

    if(pwrinit.eq.1)then
       if((NKMAX.eq.NKMAX_save).and.&
         &(NNMAX.eq.NNMAX_save)) then
          return
       else
          call wfpwr_deallocate
       end if
    end if
    
    allocate(PABSS(NSM),PABSK(NKMAX),PABSTN(NNMAX))
    allocate(PABSSN(NNMAX,NSM),PABSKN(NNMAX,NKMAX),PFV(NNMAX,3))
    allocate(RZCO(NSM,NNMAX))

    NKMAX_save = NKMAX
    NNMAX_save = NNMAX
    pwrinit = 1

    return
  end subroutine wfpwr_allocate
!-----
  subroutine wfpwr_deallocate
    implicit none

    IF(.NOT.ALLOCATED(PABSS)) RETURN
    deallocate(PABSS,PABSK,PABSTN,PABSSN,PABSKN,PFV,RZCO)

    return
  end subroutine wfpwr_deallocate
!-----
  subroutine wfnas_allocate
    implicit none
    integer,save :: NNMAX_save,NBMAX_save,NKMAX_save,NMMAX_save
    integer,save :: nasinit = 0

    if(nasinit.eq.0) NBMAX = NBM

    if(nasinit.eq.1) then
       if((NNMAX.eq.NNMAX_save).and.&
         &(NBMAX.eq.NBMAX_save).and.&
         &(NKMAX.eq.NKMAX_save).and.&
         &(NMMAX.eq.NMMAX_save)) then
          return
       else
          call wfnas_deallocate
       end if
    end if

    allocate(IDND(NNMAX))
    allocate(EX1WG(NBMAX),EY1WG(NBMAX),EZ1WG(NBMAX))
    allocate(EX2WG(NBMAX),EY2WG(NBMAX),EZ2WG(NBMAX))
    allocate(PWRWG(NBMAX),PHAWG(NBMAX),IDKA(NKMAX))
    allocate(IDMAT(NMMAX),IDBDY(NBMAX))
    allocate(KDKA(NKMAX),KDMAT(NMMAX),KDBDY(NBMAX))

    NNMAX_save = NNMAX
    NBMAX_save = NBMAX
    NKMAX_save = NKMAX
    NMMAX_save = NMMAX
    nasinit =1

    return
  end subroutine wfnas_allocate
!-----
  subroutine wfnas_deallocate
    implicit none
    IF(.NOT.ALLOCATED(IDND)) RETURN
    deallocate(IDND,EX1WG,EY1WG,EZ1WG,EX2WG,EY2WG,EZ2WG)
    deallocate(PWRWG,PHAWG,IDKA,IDMAT,IDBDY,KDKA,KDMAT,KDBDY)
    return
  end subroutine wfnas_deallocate
!-----
  subroutine wfwin_allocate
    implicit none
    integer,save :: NGXM_save,NGYM_save,NGVM_save
    integer,save :: wininit = 0

    if(wininit.eq.1) then
       if((NGXM.eq.NGXM_save).and.&
         &(NGYM.eq.NGYM_save).and.&
         &(NGVM.eq.NGVM_save)) then
          return
       else
          call wfwin_deallocate
       end if
    end if

    allocate(GZ(NGXM,NGYM),G2X(NGXM),G2Y(NGYM))
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

    IF(.NOT.ALLOCATED(GZ)) RETURN
    deallocate(GZ,G2X,G2Y,GV,GX)

    return
  end subroutine wfwin_deallocate

end module wfcomm
