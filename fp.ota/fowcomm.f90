module fowcomm

  USE commpi
  USE plcomm_parm
  USE fpcomm,ONLY:rkind,nrmax,nthmax,npmax,nsamax

  implicit none

  public
  integer :: nrgmax&                                         ! number of major radius grid points
            ,nzgmax&                                         ! number of vertical coordinate grid points
            ,nthpmax                                         ! number of poloidal angle grid points

  real(rkind) :: psi0

  real(rkind),allocatable,dimension(:) :: psim,&                 ! maximum poloidal magnetic flux in an orbit
                                          Fpsi,&                 ! poloidai current at psim
                                          Bout,&                 ! abs(magnetic field) at psim when psim is outside the axis
                                          Bin,&                  ! abs(magnetic field) at psim when psim is inside the axis
                                          psimg,&                ! psim for grid points
                                          Fpsig,&                ! Fpsi for grid points
                                          Boutg,&                ! Bout for grid points
                                          Bing                   ! Bin for grid points

  real(rkind),allocatable,dimension(:) :: xi,&                 ! cos(pitch angle) at psim
                                          xig                  ! xi for grid points

  type :: orbit                                              ! orbit information culcurated by TASK/OB
    integer :: nstp_max
    real(rkind),allocatable,dimension(:) :: time, &              ! time
                                            psip, &              ! poloidal magnetic flux
                                            babs, &              ! absolute value of magnetic field
                                            phi, &               ! electrostatic potential
                                            vpara, &             ! parallel velocity
                                            vperp, &             ! perpendicular velocity
                                            psit, &              ! toroidal magnetic flux
                                            zeta, &              ! toroidal angle
                                            rr, &                ! major radius
                                            zz, &                ! vertical positon
                                            rs, &                ! minor radius
                                            theta                ! poloidal angle
  end type orbit

  type(orbit),allocatable,dimension(:,:,:,:) :: orbit_p,&    ! (np,nth,nr)=(-0.5,  1 ,  1 ) to (npmax+0.5,nthmax,nrmax)
                                                orbit_th,&   ! (np,nth,nr)=(  1 ,-0.5,  1 ) to (npmax,nthmax+0.5,nrmax)
                                                orbit_r,&   ! (np,nth,nr)=(  1 ,  1 ,-0.5) to (npmax,nthmax,nrmax+0.5)
                                                orbit_m       ! (np,nth,nr)=(1,1,1) to (npmax,nthmax,nrmax)

  real(rkind),allocatable,dimension(:,:) :: Babs             ! abs(mabnetic field) in poloidal cross-section, (th,psi) plane

  real(rkind),allocatable,dimension(:,:) :: Brz,&                ! abs(magnetic field) in (R,Z) plane
                                            psirz,&              ! poloidal magnetic flux in (R,Z) plane
                                            Frz,&                ! poloidai current in (R,Z) plane
                                            BBm                  ! BBm(nth,nr) = Bout  for cos(pitch angle)>0
                                                                 !             = Bin   for cos(pitch angle)<0
  real(rkind),allocatable,dimension(:,:,:,:) :: UBspl,UFspl
                                                
contains

  subroutine fow_allocate
    use fpcomm, only:npmax,nthmax,nrmax,nsamax
    allocate(xi(nthmax))
    allocate(xig(nthmax+1))
    allocate(psim(nrmax),Fpsi(nrmax),Bout(nrmax),Bin(nrmax))
    allocate(psimg(nrmax+1),Fpsig(nrmax+1),Boutg(nrmax+1),Bing(nrmax+1))
    allocate(BBm(nthmax,nrmax))
    allocate(Babs(nthpmax,nrmax+1))
    allocate(Brz(nrgmax,nzgmax),psirz(nrgmax,nzgmax),Frz(nrgmax,nzgmax))
    allocate(UBspl(4,4,nrgmax,nzgmax),UFspl(4,4,nrgmax,nzgmax))
    allocate(orbit_p(npmax+1,nthmax,nrmax,nsamax),orbit_th(npmax,nthmax+1,nrmax,nsamax)&
            ,orbit_r(npmax,nthmax,nrmax+1,nsamax),orbit_m(npmax,nthmax,nrmax,nsamax))
  end subroutine fow_allocate

  subroutine fow_deallocate
    deallocate(xi,xig)
    deallocate(psim,Fpsi,Bout,Bin)
    deallocate(psimg,Fpsig,Boutg,Bing)
    deallocate(BBm)
    deallocate(Babs)
    deallocate(Brz,psirz,Frz)
    deallocate(UBspl,UFspl)
    deallocate(orbit_p,orbit_th,orbit_r,orbit_m)
  end subroutine fow_deallocate

  subroutine fow_read_namelist
    namelist /fow/nrgmax,nzgmax,nthpmax

    open(11,file="fpparm",status='old',action='read')
    read(11,nml=fow)
    close(11)
  end subroutine fow_read_namelist

end module fowcomm
