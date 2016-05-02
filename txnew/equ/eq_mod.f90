module equ_params
  integer(4), parameter :: irdm=2049, izdm=1025, ivdm=2049
  integer(4), parameter :: izdm2=2*izdm-1
  integer(4), parameter :: irzdm=irdm*izdm
  integer(4), parameter :: irzdm2=irdm*izdm2, isrzdm=2*(irdm+izdm2)

!  integer(4), parameter :: icvdm = 15, nsfix = 30, icvdm1 = icvdm + 1
  integer(4), parameter :: icvdm = 19, nsfix = 19, icvdm1 = icvdm + 1
  integer(4) :: iudsym

  integer(4) :: ieqrd = 16
  integer(4) :: jeqrd
  integer(4), dimension(10) :: ieqout,ieqerr

  integer(4) :: iaxis,iraxis,izaxis,nsu,isep
  real(8) :: btv,cpl,raxis,zaxis,saxis,qaxis
  real(8) :: dsep,psep,rspmx,rspmn
  real(8) :: bets,beta,betj,ccpl,ttcu,qsurf,ttpr,ell,trg,zzlp,zzli
  real(8) :: zsdw,zaav,zlen,q95,qqj,el95,prfac,sigcu
  real(8), dimension(irzdm2) :: psi,rbp
  real(8), dimension(isrzdm) :: csu,rsu,zsu

  integer(4) :: nv,nvm
  real(8), dimension(ivdm) :: hiv,siv,siw,sdw,vlv,ckv,ssv,aav,rrv,pds,fds, &
       &                      prv,rbv,qqv,arv,bbv,biv,r2b2v,shv,grbm2v,rov, &
       &                      aiv,brv,epsv,elipv,trigv

  integer(4) ::nr,nz,nrm,nsr,nsz,nszm
  real(8) :: dr,dz,dr2i,dz2i,ddri,ddzi
  real(8), dimension(irdm)   :: rg
  real(8), dimension(izdm2)  :: zg

  integer(4), dimension(10) :: icp
  real(8) :: rmaj,rpla,elip,trig,elipup,trigup,elipdw,trigdw,qaxi,qsur
  real(8), dimension(10) :: cp

  integer(4) :: ilimt
  integer(4), dimension(0:nsfix) :: ivac,ncoil
  real(8),    dimension(0:nsfix) :: cvac,rvac,zvac
  real(8), dimension(100,icvdm)     :: rcoil,zcoil,ccoil
  real(8), dimension(200)           :: rlimt,zlimt

end module equ_params
 
