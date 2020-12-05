module equ_params
  USE bpsd_kinds,ONLY: dp
  integer, parameter :: irdm=2049, izdm=1025, ivdm=2049
  integer, parameter :: izdm2=2*izdm-1
  integer, parameter :: irzdm=irdm*izdm
  integer, parameter :: irzdm2=irdm*izdm2, isrzdm=2*(irdm+izdm2)

!  integer, parameter :: icvdm = 15, nsfix = 30, icvdm1 = icvdm + 1
  integer, parameter :: icvdm = 19, nsfix = 19, icvdm1 = icvdm + 1
  integer :: iudsym
  character(len=80) :: eqfile = 'eqdb.dat'

!  integer :: ieqrd = 16
  integer :: ieqrd = 17
  integer :: jeqrd
  integer, dimension(:), allocatable :: ieqout,ieqerr

  integer :: iaxis,iraxis,izaxis,nsu,isep
  real(dp) :: btv,cpl,raxis,zaxis,saxis,qaxis
  real(dp) :: dsep,psep,rspmx,rspmn
  real(dp) :: bets,beta,betj,ccpl,ttcu,qsurf,ttpr,ell,trg,zzlp,zzli
  real(dp) :: zsdw,zaav,zlen,q95,qqj,el95,prfac,sigcu
  real(dp), dimension(:), allocatable :: psi,rbp
  real(dp), dimension(:), allocatable :: csu,rsu,zsu

  integer :: nv,nvm
  real(dp), dimension(:), allocatable :: pds,fds,vlv,qqv,prv
  real(dp), dimension(:), allocatable :: hiv,siv,siw,sdw,ckv,ssv,aav,rrv, &
       &                      rbv,arv,bbv,biv,r2b2v,shv,grbm2v, &
       &                      rov,aiv,brv,epsv,elipv,trigv,ftv

  integer ::nr,nz,nrm,nsr,nsz,nszm
  real(dp) :: dr,dz,dr2i,dz2i,ddri,ddzi
  real(dp), dimension(:), allocatable :: rg, zg

  integer, dimension(:), allocatable :: icp
  real(dp) :: rmaj,rpla,elip,trig,elipup,trigup,elipdw,trigdw,qaxi,qsur
  real(dp), dimension(:), allocatable :: cp

  integer :: ilimt
  integer, dimension(:), allocatable :: ivac,ncoil
  real(dp),    dimension(:), allocatable :: cvac,rvac,zvac
  real(dp), dimension(:,:),  allocatable :: rcoil,zcoil,ccoil
  real(dp), dimension(:),    allocatable :: rlimt,zlimt

end module equ_params
 
