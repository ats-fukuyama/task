module equ_params
  implicit none
  integer(4), parameter :: irdm=2049, izdm=1025, ivdm=2049
  integer(4), parameter :: izdm2=2*izdm-1
  integer(4), parameter :: isrzdm=2*(irdm+izdm2)
  integer(4), parameter :: nmax = 1001 ! division number of a flux surface
  real(8) :: tol = 1.d-8

!  integer(4), parameter :: icvdm = 15, nsfix = 30, icvdm1 = icvdm + 1
  integer(4), parameter :: icvdm = 19, nsfix = 19, icvdm1 = icvdm + 1
  integer(4) :: iudsym
  character(len=80) :: eqfile = 'eqdata'

  integer(4) :: iraxis,izaxis
  real(8) :: dr,dz
  real(8) :: btv,raxis,zaxis,saxis
  real(8) :: dsep,psep,rspmx,rspmn
  real(8) :: ccpl,qsurf,ell,trg,zzlp,zzli
  real(8) :: q95,qqj,el95,prfac,sigcu
  real(8) :: psign
  real(8), dimension(:), allocatable :: psi,rbp
  real(8), dimension(:), allocatable :: csu,rsu,zsu
  ! -------------------------------------------------
  real(8), dimension(:,:,:,:), allocatable :: upsi

  ! --- Miller parameters ----------------------------
  type miller
     real(8) :: Rgeo = 0.d0, Zgeo = 0.d0, rminor = 0.d0
     real(8) :: Rmax = 0.d0, Rmin = 0.d0, zpla = 0.d0
     real(8) :: elip = 0.d0, trig = 0.d0, zeta = 0.d0
  end type miller
  type(miller), dimension(:), allocatable :: vmiller
  ! --------------------------------------------------

  integer(4) :: nv,nsu
  real(8), dimension(:), allocatable :: pds,fds,vlv,qqv,prv
  real(8), dimension(:), allocatable :: hiv,siv,siw,sdw,ckv,ssv,aav,rrv, &
       &                      rbv,arv,bbv,biv,r2b2v,shv,grbm2v, &
       &                      rov,aiv,brv,epsv,elipv,trigv,ftv,rtv,rpv
  real(8), dimension(:), allocatable :: gttiv

  integer(4) ::nsr,nsz
  real(8), dimension(:), allocatable :: rg, zg

  integer(4), dimension(:), allocatable :: icp
  real(8) :: rmaj,rpla,elip,trig,elipup,trigup,elipdw,trigdw,qaxi,qsur
  real(8), dimension(:), allocatable :: cp

  ! -----

  real(8), dimension(:), allocatable :: cvac,dsr,dsz

end module equ_params
 
