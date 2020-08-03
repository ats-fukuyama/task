module fowcomm

  USE commpi
  USE plcomm_parm
  USE fpcomm,only:rkind

  implicit none

  public
  integer :: npsmax&                                         ! number of psip grid points
            ,nzemax&                                         ! number of zeta(=cos(theta_m)) grid points
            ,nrgmax&                                         ! number of major radius grid points
            ,nzgmax&                                         ! number of vertical coordinate grid points
            ,nthpmax                                         ! number of poloidal angle grid points

  real(8) :: psi0

  real(8),allocatable,dimension(:) :: psim,&                 ! maximum poloidal magnetic flux in an orbit
                                      Fpsi,&                 ! poloidai current at psim
                                      Bout,&                 ! abs(magnetic field) at psim when psim is outside the axis
                                      Bin,&                  ! abs(magnetic field) at psim when psim is inside the axis
                                      psimg,&                ! psim for grid points
                                      Fpsig,&                ! Fpsi for grid points
                                      Boutg,&                ! Bout for grid points
                                      Bing                   ! Bin for grid points

  real(8),allocatable,dimension(:) :: zeta,&                 ! cos(pitch angle) at psim
                                      zetag                  ! zeta for grid points

  type :: orbit                                            ! orbit information culcurated by TASK/OB
    integer :: nstp
    real(8),allocatable,dimension(:) :: time, &              ! time
                                        psip, &              ! poloidal magnetic flux
                                        babs, &              ! absolute value of magnetic field
                                        phi, &               ! electrostatic potential
                                        vpara, &             ! parallel velocity
                                        vperp, &             ! perpendicular velocity
                                        psit, &              ! toroidal magnetic flux
                                        toloidal_angle, &    ! toroidal angle
                                        rr, &                ! major radius
                                        zz, &                ! vertical positon
                                        rs, &                ! minor radius
                                        theta                ! poloidal angle
  end type orbit

  type(orbit),allocatable,dimension(:,:,:,:) :: orbit_p,&    ! (np,nze,nps)=(-0.5,  1 ,  1 ) to (npmax+0.5,nzemax,npsmax)
                                                orbit_ze,&   ! (np,nze,nps)=(  1 ,-0.5,  1 ) to (npmax,nzemax+0.5,npsmax)
                                                orbit_ps,&   ! (np,nze,nps)=(  1 ,  1 ,-0.5) to (npmax,nzemax,npsmax+0.5)
                                                orbitm       ! (np,nze,nps)=(1,1,1) to (npmax,nzemax,npsmax)

  real(rkind),allocatable,dimension(:,:) :: Babs             ! abs(mabnetic field) in poloidal cross-section, (th,psi) plane

  real(8),allocatable,dimension(:,:) :: Brz,&                ! abs(magnetic field) in (R,Z) plane
                                        psirz,&              ! poloidal magnetic flux in (R,Z) plane
                                        Frz,&                ! poloidai current in (R,Z) plane
                                        BBm                  ! BBm(nze,nps) = Bout  for cos(pitch angle)>0
                                                             !              = Bin   for cos(pitch angle)<0
  real(8),allocatable,dimension(:,:,:,:) :: UBspl,UFspl
                                                
contains

  subroutine fow_allocate
    use fpcomm, only:npmax,nsamax
    allocate(zeta(nzemax))
    allocate(zetag(nzemax+1))
    allocate(psim(npsmax),Fpsi(npsmax),Bout(npsmax),Bin(npsmax))
    allocate(psimg(npsmax+1),Fpsig(npsmax+1),Boutg(npsmax+1),Bing(npsmax+1))
    allocate(BBm(nzemax,npsmax))
    allocate(Babs(nthpmax,npsmax+1))
    allocate(Brz(nrgmax,nzgmax),psirz(nrgmax,nzgmax),Frz(nrgmax,nzgmax))
    allocate(UBspl(4,4,nrgmax,nzgmax),UFspl(4,4,nrgmax,nzgmax))
    allocate(orbit_p(npmax+1,nzemax,npsmax,nsamax),orbit_ze(npmax,nzemax+1,npsmax,nsamax)&
            ,orbit_ps(npmax,nzemax,npsmax+1,nsamax),orbitm(npmax,nzemax,npsmax,nsamax))
  end subroutine fow_allocate

  subroutine fow_deallocate
    deallocate(zeta,zetag)
    deallocate(psim,Fpsi,Bout,Bin)
    deallocate(psimg,Fpsig,Boutg,Bing)
    deallocate(BBm)
    deallocate(Babs)
    deallocate(Brz,psirz,Frz)
    deallocate(UBspl,UFspl)
    deallocate(orbit_p,orbit_ze,orbit_ps,orbitm)
  end subroutine fow_deallocate

  subroutine fow_read_namelist
    namelist /fow/npsmax,nzemax,nrgmax,nzgmax,nthpmax

    open(11,file="fpparm",status='old',action='read')
    read(11,nml=fow)
    close(11)
  end subroutine fow_read_namelist

  function construct_orbit(n,input) result(ret)
    type(orbit) :: ret
    integer, intent(in) :: n
    real(8), intent(in) :: input(:,:)
    integer :: i

    allocate(ret%time(n),ret%psip(n),ret%babs(n),ret%phi(n),ret%vpara(n),ret%vperp(n)&
            ,ret%psit(n),ret%toloidal_angle(n),ret%rr(n),ret%zz(n),ret%rs(n),ret%theta(n))

    ret%nstp=n
    do i=1,n
      ret%time(i)=input(1,i)
      ret%psip(i)=input(2,i)
      ret%babs(i)=input(3,i)
      ret%phi(i)=input(4,i)
      ret%vpara(i)=input(5,i)
      ret%vperp(i)=input(6,i)
      ret%psit(i)=input(7,i)
      ret%toloidal_angle(i)=input(8,i)
      ret%rr(i)=input(9,i)
      ret%zz(i)=input(10,i)
      ret%rs(i)=input(11,i)
      ret%theta(i)=input(12,i)
    end do

  end function construct_orbit

  subroutine destruct_orbit(orbit_in)
    type(orbit) :: orbit_in
    deallocate(orbit_in%time,orbit_in%psip,orbit_in%babs,orbit_in%phi,orbit_in%vpara,orbit_in%vperp&
            ,orbit_in%psit,orbit_in%toloidal_angle,orbit_in%rr,orbit_in%zz,orbit_in%rs,orbit_in%theta)
  end subroutine destruct_orbit

end module fowcomm
