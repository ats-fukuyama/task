! wfcomm.f90

!   Define and allocate global variables

MODULE wfcomm_parm

  USE plcomm_parm
  USE dpcomm_parm

  IMPLICIT NONE

  ! === input parameter array size ====
  
  INTEGER,PARAMETER:: nantm =       8 ! maximum number of antenna
  INTEGER,PARAMETER:: npointm =  2000 ! maximum number of antenna path points
  INTEGER,PARAMETER:: nmedm =      10 ! maximum number of medium type
  INTEGER,PARAMETER:: idebuga_max=100 ! maximum number of debug parameter array

  ! === wfdiv parameters ===

  INTEGER:: model_config,model_shape
  
  ! --- model_config=1 : rectangular ---
  
  REAL(rkind):: xdiv_min,xdiv_max,ydiv_min,ydiv_max
  REAL(rkind):: delx,dely
  REAL(rkind):: rdiv_min,rdiv_max,thdiv_min,thdiv_max

  ! === RF parameters ===
  
  REAL(rkind):: RF,RKZ
  INTEGER:: NPH

  ! === antenna paramters ===

  INTEGER:: nant_max
  REAL(rkind):: AJ(nantm),APH(nantm),AWD(nantm),APOS(nantm)
  REAL(rkind):: PIN(nantm),RD(nantm),THETJ1(nantm),THETJ2(nantm)

  ! === waveguide parameters ===

  INTEGER:: model_wg
  REAL(rkind):: xwg_min,xwg_max,ywg_min,ywg_max,thwg_min,thwg_max
  REAL(rkind):: phase_wg_min,phase_wg_cen,phase_wg_max
  REAL(rkind):: amp_wg,angle_wg,ellip_wg

  ! === edge dump parameters ===

  INTEGER:: model_damp
  REAL(rkind):: xdamp_min,xdamp_max,ydamp_min,ydamp_max
  REAL(rkind):: thdamp_min,thdamp_max,width_damp,factor_damp

  ! === collsion inhancement parameters ===

  INTEGER:: model_coll_enhance
  REAL(rkind):: &
       factor_coll_enhance, &
       xpos_coll_enhance,xwidth_coll_enhance, &
       ypos_coll_enhance,ywidth_coll_enhance

  ! === interpolation level ===

  INTEGER:: model_interpolation

  ! === seg field parity ===

  INTEGER:: model_wf    ! 0: positive, 1: alternative

  ! === medium parameters ===

  INTEGER:: nmed_max
  INTEGER:: model_nmed(nmedm)
  REAL(rkind):: epsilon_nmed(nmedm),amu_nmed(nmedm),sigma_nmed(nmedm)
  REAL(rkind):: xmin_nmed(nmedm),xmax_nmed(nmedm)
  REAL(rkind):: ymin_nmed(nmedm),ymax_nmed(nmedm)
  REAL(rkind):: rmin_nmed(nmedm),rmax_nmed(nmedm)
  REAL(rkind):: thmin_nmed(nmedm),thmax_nmed(nmedm)

  ! === output parameters ===

  INTEGER:: nprint,ngraph,ndrawd,ndraws,ndrawa,ndrawe,ndrawv

  ! === sort weight parameters ===

  REAL(rkind):: sort_weight_x,sort_weight_y

  ! === graphic paramters ===

  INTEGER:: ngxmax,ngymax,ngvmax
  REAL(rkind):: gaspect

  ! === division control parameter ===

  INTEGER:: nxzone_max,nyzone_max

  ! === computation control parameters ===

  REAL(rkind):: tolerance

  ! === controll parameters ====
  
  INTEGER:: MODELI

  ! === input file name ===

  CHARACTER(LEN=80):: KFNAME,KFNAMA,KFNAMF,KFNAMB,KNAMWG

  ! === debug controll parameters ===

  INTEGER:: idebuga(idebuga_max)

END MODULE wfcomm_parm

MODULE wfcomm

  USE wfcomm_parm
  USE femcomm
  USE commpi
  IMPLICIT NONE
  
  ! === for parallel processing ===

  INTEGER :: istart,iend

  ! === 

!      ----- graphics -----

!  integer :: NBPM
!  integer,parameter::  NMDM = 10
!  integer,parameter::  NCNM = 12+NMDM

!  integer,parameter::   NGM =   3
  
  integer,parameter::  NWDM = 12
  integer,parameter::  NCHM = 80

!       --- common variables ---

  COMPLEX(rkind):: CII

  REAL(rkind):: PSIA

  ! --- magnetic field relation between nodes and elements ---

  REAL(rkind),ALLOCATABLE:: bx_nelm(:) ! x component of static magnetic field 
  REAL(rkind),ALLOCATABLE:: by_nelm(:) ! y component of static magnetic field 
  REAL(rkind),ALLOCATABLE:: bz_nelm(:) ! z component of static magnetic field 
  REAL(rkind),ALLOCATABLE:: babs_nelm(:) ! babs=SQRT(bx**2+by**2+bz**2)

  ! === field variables ===

  COMPLEX(rkind),ALLOCATABLE:: CESD_nseg(:)
  COMPLEX(rkind),ALLOCATABLE:: CEND_node(:)
  COMPLEX(rkind),ALLOCATABLE:: CESD_nbdy(:)
  COMPLEX(rkind),ALLOCATABLE:: CEND_nbdy(:)

  COMPLEX(rkind),ALLOCATABLE:: CBF_node(:,:)
  COMPLEX(rkind),ALLOCATABLE:: CBP_node(:,:)
        
  REAL(rkind),ALLOCATABLE:: pabs_ns_nelm(:,:)
  REAL(rkind),ALLOCATABLE:: pabs_ns(:)
  REAL(rkind),ALLOCATABLE:: pabs_tot
        
  ! === matrix variables ===
  
  INTEGER:: mtx_len
  INTEGER:: NNBMAX
  COMPLEX(rkind):: CM(6,6)
  COMPLEX(rkind),ALLOCATABLE :: CSV(:)      !(MLEN)
  COMPLEX(rkind),ALLOCATABLE :: CQQ(:)      !(MBND)
  COMPLEX(rkind),ALLOCATABLE :: CVTOT(:,:)  !(6,nelm_max)


  ! === antenna variables ===
  
  INTEGER :: np0_max_nant(nantm)
  INTEGER :: nelm_np0_nant(npointm,nantm)
  INTEGER :: np_max_nant(nantm)
  INTEGER :: nelm_np_nant(npointm,nantm)
  COMPLEX(rkind) :: cimp_nant(nantm)
  COMPLEX(rkind) :: cimp_tot
  REAL(rkind),DIMENSION(npointm,nantm):: x_np0_nant,y_np0_nant
  REAL(rkind),DIMENSION(npointm,nantm):: x_np_nant,y_np_nant
        
!       /WFWIN/
  real(rkind):: len_min,len_max
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

  LOGICAL wf_solve_allocated
  LOGICAL wf_win_allocated

! -----------------------------------------------------------------

contains

  ! === allocate arrays for solve ===
  
  SUBROUTINE wf_solve_allocate
    IMPLICIT NONE
    INTEGER,SAVE:: node_max_save=0
    INTEGER,SAVE:: nelm_max_save=0
    INTEGER,SAVE:: nseg_max_save=0
    INTEGER,SAVE:: nbdy_max_save=0
    INTEGER,SAVE:: nsmax_save=0
    INTEGER,SAVE:: mtx_len_save=0

    IF(wf_solve_allocated) THEN
       IF((node_max.EQ.node_max_save).AND. &
          (nelm_max.EQ.nelm_max_save).AND. &
          (nseg_max.EQ.nseg_max_save).AND. &
          (nbdy_max.EQ.nbdy_max_save).AND. &
          (nsmax   .EQ.nsmax_save).AND. &
          (mtx_len .EQ.mtx_len_save )) RETURN
       CALL wf_solve_deallocate
    END IF

    ALLOCATE(CSV(mtx_len),CVTOT(6,nelm_max))
    ALLOCATE(CESD_nseg(nseg_max),CEND_node(node_max))
    ALLOCATE(CESD_nbdy(nbdy_max),CEND_nbdy(nbdy_max))
    ALLOCATE(CBF_node(3,node_max),CBP_node(3,node_max))
    ALLOCATE(pabs_ns_nelm(nsmax,nelm_max))
    ALLOCATE(pabs_ns(nsmax))

    node_max_save=node_max
    nelm_max_save=nelm_max
    nseg_max_save=nseg_max
    nsmax_save=nsmax
    mtx_len_save =mtx_len

    wf_solve_allocated=.TRUE.
    
    RETURN
  END SUBROUTINE wf_solve_allocate

  SUBROUTINE wf_solve_deallocate
    IMPLICIT NONE

    DEALLOCATE(CSV,CVTOT)
    DEALLOCATE(CESD_nseg,CEND_node)
    DEALLOCATE(CESD_nbdy,CEND_nbdy)
    DEALLOCATE(CBF_node,CBP_node)
    DEALLOCATE(pabs_ns_nelm,pabs_ns)

    wf_solve_allocated=.FALSE.

    RETURN
  END SUBROUTINE wf_solve_deallocate

  ! === allocate arrays for graphics ===

  SUBROUTINE wf_win_allocate
    IMPLICIT NONE
    INTEGER,SAVE :: ngxmax_save=0
    INTEGER,SAVE :: ngymax_save=0
    INTEGER,SAVE :: ngvmax_save=0

    IF(wf_win_allocated) THEN
       IF ((ngxmax.EQ.ngxmax_save).AND. &
           (ngymax.EQ.ngymax_save).AND. &
           (ngvmax.EQ.ngvmax_save)) RETURN
       CALL wf_win_deallocate
    END IF

    ALLOCATE(gz(ngxmax,ngymax),gz_temp(ngxmax,ngymax),g2x(ngxmax),g2y(ngymax))
    ALLOCATE(iegz(ngxmax,ngymax))
    ALLOCATE(gv(ngvmax,3),gx(ngvmax))

    ngxmax_save = ngxmax
    ngymax_save = ngymax
    ngvmax_save = ngvmax

    wf_win_allocated=.TRUE.
    
    RETURN
  END SUBROUTINE wf_win_allocate

  SUBROUTINE wf_win_deallocate
    IMPLICIT none

    DEALLOCATE(GZ,GZ_temp,G2X,G2Y,GV,GX,IEGZ)

    wf_win_allocated=.FALSE.

    RETURN
  END SUBROUTINE wf_win_deallocate

END MODULE wfcomm
