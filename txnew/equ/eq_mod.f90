module equ_params
  integer(4), parameter :: irdm=2049, izdm=1025, ivdm=2049
  integer(4), parameter :: izdm2=2*izdm-1
  integer(4), parameter :: irzdm=irdm*izdm
  integer(4), parameter :: irzdm2=irdm*izdm2, isrzdm=2*(irdm+izdm2)

!  integer(4), parameter :: icvdm = 15, nsfix = 30, icvdm1 = icvdm + 1
  integer(4), parameter :: icvdm = 19, nsfix = 19, icvdm1 = icvdm + 1
  integer(4) :: iudsym
  character(len=80) :: eqfile = 'eqdata'

  integer(4) :: ieqrd = 16
  integer(4) :: jeqrd
  integer(4), dimension(:), allocatable :: ieqout,ieqerr

  integer(4) :: iaxis,iraxis,izaxis,nsu,isep
  real(8) :: btv,cpl,raxis,zaxis,saxis,qaxis
  real(8) :: dsep,psep,rspmx,rspmn
  real(8) :: bets,beta,betj,ccpl,ttcu,qsurf,ttpr,ell,trg,zzlp,zzli
  real(8) :: zsdw,zaav,zlen,q95,qqj,el95,prfac,sigcu
  real(8), dimension(:), allocatable :: psi,rbp
  real(8), dimension(:), allocatable :: csu,rsu,zsu

  integer(4) :: nv,nvm
  real(8), dimension(:), allocatable :: pds,fds,vlv,qqv,prv
  real(8), dimension(:), allocatable :: hiv,siv,siw,sdw,ckv,ssv,aav,rrv, &
       &                      rbv,arv,bbv,biv,r2b2v,shv,grbm2v, &
       &                      rov,aiv,brv,epsv,elipv,trigv,ftv

  integer(4) ::nr,nz,nrm,nsr,nsz,nszm
  real(8) :: dr,dz,dr2i,dz2i,ddri,ddzi
  real(8), dimension(:), allocatable :: rg, zg

  integer(4), dimension(:), allocatable :: icp
  real(8) :: rmaj,rpla,elip,trig,elipup,trigup,elipdw,trigdw,qaxi,qsur
  real(8), dimension(:), allocatable :: cp

  integer(4) :: ilimt
  integer(4), dimension(:), allocatable :: ivac,ncoil
  real(8),    dimension(:), allocatable :: cvac,rvac,zvac
  real(8), dimension(:,:),  allocatable :: rcoil,zcoil,ccoil
  real(8), dimension(:),    allocatable :: rlimt,zlimt

end module equ_params
 
