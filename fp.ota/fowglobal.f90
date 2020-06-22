module fow_global

  implicit none

  integer :: npsmax,nzemax,nrgmax,nzgmax
  real(8) :: psi0

  ! nzemax
  real(8),allocatable,dimension(:) :: zeta
  ! npsmax
  real(8),allocatable,dimension(:) :: psi,Fpsi,Bout,Bin
  ! nrgmax * nzgmax
  real(8),allocatable,dimension(:,:) :: Brz,psirz,Frz
  ! 4 * 4 * nrgmax * nzgmax
  real(8),allocatable,dimension(:,:,:,:) :: UBspl,UFspl

contains

  subroutine fow_allocate
    allocate(zeta(nzemax))
    allocate(psi(npsmax),Fpsi(npsmax),Bout(npsmax),Bin(npsmax))
    allocate(Brz(nrgmax,nzgmax),psirz(nrgmax,nzgmax),Frz(nrgmax,nzgmax))
    allocate(UBspl(4,4,nrgmax,nzgmax),UFspl(4,4,nrgmax,nzgmax))
  end subroutine fow_allocate

  subroutine fow_deallocate
    deallocate(zeta)
    deallocate(psi,Fpsi,Bout,Bin)
    deallocate(Brz,psirz,Frz)
    deallocate(UBspl,UFspl)
  end subroutine fow_deallocate

  subroutine fow_read_namelist
    namelist /fow/npsmax,nzemax,nrgmax,nzgmax

    open(11,file="fpparm",status='old',action='read')
    read(11,nml=fow)
    close(11)
  end subroutine fow_read_namelist

end module fow_global
