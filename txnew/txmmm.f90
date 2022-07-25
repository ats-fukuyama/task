subroutine txmmm95(dNsdrho,dTsdrho,dQdrho,cexb,gamma)
  use tx_commons, only : NRMAX, rpt, RR, achg, PNbV, Zeff, Var, Q, &
       &                 BphV, amas, Chis, Dfs, rMus, FSANOM, &
       &                 FSPCL, rKeV, AEE, NRA, BthV, wexb
!  use tx_interface, only : dfdx
  implicit none

  ! 0:NRMAX => 1:NRMAX+1 if dimension(*) used
  real(8), intent(in) :: cexb
  real(8), dimension(:,:), intent(in) :: dNsdrho, dTsdrho
  real(8), dimension(:),   intent(in) :: dQdrho
  real(8), dimension(0:NRMAX), intent(out), optional :: gamma

  integer(4), parameter :: mxmode = 12 ! at least 5
  integer(4) :: N, NR, igamma, i
  integer(4) :: matdim, npoints, nprout, lprint, nerr, lsuper, lreset
  integer(4), dimension(:), allocatable :: lswitch
  real(4), dimension(:), allocatable :: &
       & cswitch, fig, frb, fkb, &
       & zrmajor, zelong, zdense, zdensh, zdensimp, zdensfe, zxzeff, zavezimp, zmassimp, zmasshyd, &
       & zaimass, zwexbs, zgrdne, zgrdni, zgrdnh, zgrdnz, zgrdte, zgrdti, zgrdq, zthiig, zthdig, &
       & ztheig, zthzig, zthirb, zthdrb, ztherb, zthzrb, zthikb, zthdkb, zthekb, zthzkb
  real(4), dimension(:,:), allocatable :: zgamma, zomega, zvelthi, zvflux
  real(4), dimension(:,:,:), allocatable :: zdifthi
  real(8) :: factor_bohm, DeL, cap_val

  N = NRMAX ! used for abbreviation

  if(present(gamma)) then
     igamma = 1
  else
     igamma = 0
  end if

  allocate(lswitch(1:8))
  allocate(cswitch(1:25),fig(1:4),frb(1:4),fkb(1:4))
  allocate(zrmajor(0:N), source=0.0)
  allocate(zelong,zdense, zdensh,zdensimp,zdensfe,zxzeff,zavezimp, &
       &   zmassimp,zmasshyd,zaimass,zwexbs,zgrdne, &
       &   zgrdni,zgrdnh,zgrdnz,zgrdte,zgrdti,zgrdq, &
       &   zthiig,zthdig,ztheig,zthzig,zthirb,zthdrb, &
       &   ztherb,zthzrb,zthikb,zthdkb,zthekb,zthzkb, source=zrmajor)
  allocate(zgamma(1:mxmode,0:N),zomega(1:mxmode,0:N),zdifthi(1:mxmode,1:mxmode,0:N), &
       &   zvelthi(1:mxmode,0:N),zvflux(1:mxmode,0:N))

  !     *** Dummy high Z impurity ***
  !         No impurities induce no ITG growth rate.

  zrmajor(0:N)  = real(RR)
  zelong(0:N)   = 1.0

  zdense(0:N)   = real(Var(0:N,1)%n) * 1.e20
  zdensh(0:N)   = real(Var(0:N,2)%n) * 1.e20
  zdensimp(0:N) = real(Var(0:N,3)%n) * 1.e20
  zdensfe(0:N)  = real(PNbV(0:N))    * 1.e20

  zxzeff(0:N)   = real(Zeff(0:N))

  zavezimp(0:N) = real(achg(3))
  zmassimp(0:N) = real(amas(3))
  zmasshyd(0:N) = real(amas(2))
  zaimass(0:N)  = ( zmasshyd(0:N) * zdensh(0:N) + zmassimp(0:N) * zdensimp(0:N) ) &
       &        / ( zdensh(0:N) + zdensimp(0:N) )

  if(cexb == 0.d0) then
     zwexbs(0:N)   = 0.0
  else
     zwexbs(0:N)   = real(abs(cexb) * wexb(0:N))
  end if

  ! Electron density gradient
  zgrdne(0:N)   = - zrmajor(0:N) * real((dNsdrho(1:N+1,1)) / Var(0:N,1)%n)
  ! Hydrogenic ion density gradient (excluding impurities)
  zgrdnh(0:N)   = - zrmajor(0:N) * real((dNsdrho(1:N+1,2)) / Var(0:N,2)%n)
  ! Impurity ion density gradient
  zgrdnz(0:N)   = - zrmajor(0:N) * real((dNsdrho(1:N+1,3)) / Var(0:N,3)%n)
  ! Thermal ion density gradient (including all thermal ions)
  zgrdni(0:N)   = zgrdnh(0:N) + zgrdnz(0:N)
  ! Electron temperature gradient
  zgrdte(0:N)   = - zrmajor(0:N) * real((dTsdrho(1:N+1,1)) / Var(0:N,1)%T)
  ! Ion temperature gradient
  zgrdti(0:N)   = - zrmajor(0:N) * real((dTsdrho(1:N+1,2)) / Var(0:N,2)%T)
  ! Safety factor gradient (~ magnetic shear)
  zgrdq(0:N)    =   zrmajor(0:N) * real((dQdrho(1:N+1) ) / Q(0:N))

  matdim  = mxmode
  npoints = N + 1
  nprout  = 6       ! output unit number for long printout

  lprint  = 0       ! controls amount of printout
  lsuper  = 0       ! used for supershots (reducing kinetic ballooning mode significantly)
  lreset  = 0       ! use default values for cswitch and lswitch

  call mmm95( &
     ! Input:
     &   real(rpt),  zrmajor,   zelong &
     & , zdense,   zdensh,    zdensimp,  zdensfe &
     & , zxzeff,   real(Var(0:N,1)%T),real(Var(0:N,2)%T),real(Q),real(BphV) &
     & , zavezimp, zmassimp,  zmasshyd,  zaimass,  zwexbs &
     & , zgrdne,   zgrdni,    zgrdnh,    zgrdnz,   zgrdte,   zgrdti,  zgrdq &
     ! Output:
     & , zthiig,   zthdig,    ztheig,    zthzig &
     & , zthirb,   zthdrb,    ztherb,    zthzrb &
     & , zthikb,   zthdkb,    zthekb,    zthzkb &
     & , zgamma,   zomega,    zdifthi,   zvelthi,  zvflux &
     ! Input integers and switches, output:
     & , matdim,  npoints,  nprout,   lprint,  nerr &
     ! Internal control variables:
     & , lsuper,  lreset,   lswitch,  cswitch, fig,    frb,     fkb)
  if(nerr .ne. 0) write(6,*) 'mmm95 outputs some error. nr,nerr=', nr, nerr

  if(FSANOM(3) /= 0.d0 .and. igamma == 0) then
     !  Ion thermal diffusivity
     Chis(0:N,2) = FSANOM(3) * (  real(zthiig(0:N),8) &  ! ITG mode
          &                     + real(zthirb(0:N),8) &  ! Resistive Ballooning mode
          &                     + real(zthikb(0:N),8))   ! Kinetic Ballooning mode
!!$     do nr=0,n
!!$        write(6,*) rho(nr),zthiig(nr),zthirb(nr),zthikb(nr)
!!$     end do

     !  Electron thermal diffusivity
     Chis(0:N,1) = FSANOM(3) * (  real(ztheig(0:N),8) &
          &                     + real(ztherb(0:N),8) &
          &                     + real(zthekb(0:N),8))
  end if

  !  Momentum viscosity regarded as the same as thermal diffusivity
  if(FSANOM(2) /= 0.d0 .and. igamma == 0) then
     rMus(0:N,1) = FSANOM(2) * Chis(0:N,1)
     rMus(0:N,2) = FSANOM(2) * Chis(0:N,2)
  end if

  !  Particle diffusivity
  if(FSANOM(1) /= 0.d0 .and. igamma == 0) then
     Dfs(0:N,1) = FSANOM(1) * (  real(zthdig(0:N),8) &
          &                    + real(zthdrb(0:N),8) &
          &                    + real(zthdkb(0:N),8))
  end if

  if(igamma == 1) then
     ! Most unstable mode at each radial location
     do NR = 0, NRMAX
        gamma(NR) = real(maxval(zgamma(1:mxmode,NR)),8)
     end do
  end if

  deallocate(lswitch,cswitch,fig,frb,fkb, &
       &     zrmajor,zelong,zdense,zdensh,zdensimp,zdensfe,zxzeff,zavezimp,zmassimp,zmasshyd, &
       &     zaimass,zwexbs,zgrdne,zgrdni,zgrdnh,zgrdnz,zgrdte,zgrdti,zgrdq,zthiig,zthdig, &
       &     ztheig,zthzig,zthirb,zthdrb,ztherb,zthzrb,zthikb,zthdkb,zthekb,zthzkb, &
       &     zgamma,zomega,zdifthi,zvelthi,zvflux)

  ! Outside the separatrix

  if(FSPCL(1) == 0.d0) then
     do NR = NRA, NRMAX
        factor_bohm = Dfs(NR,1) / (Var(NRA,1)%T * rKeV &
             &      / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
        Dfs(NR,1) = factor_bohm * Var(NR,1)%T * rKeV / (16.D0 * AEE * SQRT(BphV(NR)**2+BthV(NR)**2))
     end do
  else
     do NR = NRA + 1, NRMAX
        Dfs(NR,1) = Dfs(NRA,1)
     end do
  end if
     
  if(FSPCL(3) == 0.d0) then
     do NR = NRA, NRMAX
        factor_bohm = Chis(NR,1) / (Var(NRA,1)%T * rKeV &
             &      / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
        DeL = factor_bohm * Var(NR,1)%T * rKeV / (16.D0 * AEE * SQRT(BphV(NR)**2+BthV(NR)**2))
        Chis(NR,1) = DeL
        rMus(NR,1) = DeL
        factor_bohm = Chis(NR,2) / (Var(NRA,1)%T * rKeV &
             &      / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
        DeL = factor_bohm * Var(NR,1)%T * rKeV / (16.D0 * AEE * SQRT(BphV(NR)**2+BthV(NR)**2))
        Chis(NR,2) = DeL
        rMus(NR,2) = DeL
     end do
  else
     do NR = NRA + 1, NRMAX
        do i = 1, 2
           Chis(NR,i) = Chis(NRA,i)
           rMus(NR,i) = rMus(NRA,i)
        end do
     end do
  end if

  ! Limitation of diffusivities

  cap_val = 5.d1
  do i = 1, 2
     where(Chis(:,i) > cap_val) Chis(:,i) = cap_val
     where(rMus(:,i) > cap_val) rMus(:,i) = cap_val
  end do

  ! Minimum diffusivities for numerical purpose

  cap_val = 0.1d0
  do i = 1, 2
     where(Chis(:,i) <= 0.d0) Chis(:,i) = cap_val
     where(rMus(:,i) <= 0.d0) rMus(:,i) = cap_val
  end do

end subroutine txmmm95

!=======================================================================

subroutine ITG_growthrate
  use tx_commons, only : NRMAX, NSM, Q, Var, Rho, vv, wexb, gamITG, RR, S, vro, array_init_NRNS
  use tx_interface, only : dfdx, txmmm95
  implicit none

  integer(4) :: N, NR, i
  real(8)    :: epsN
  real(8), dimension(:,:), allocatable :: dTsdrho, dNsdrho
  real(8), dimension(:),   allocatable :: dQdrho, gamma, gamma_exb

  N = NRMAX

  allocate(dTsdrho,dNsdrho,mold=array_init_NRNS)
  allocate(dQdrho,gamma,gamma_exb,mold=vv)

  dQdrho (0:N) = vro(0:N) * dfdx(vv,Q           ,N,0)
  do i = 1, NSM
     dTsdrho(0:N,i) = vro(0:N) * dfdx(vv,Var(0:N,i)%T,N,0)
     dNsdrho(0:N,i) = vro(0:N) * dfdx(vv,Var(0:N,i)%n,N,0)
  end do

  ! With ExB shearing effect
  call txmmm95(dNsdrho,dTsdrho,dQdrho,1.d0,gamma_exb)

  ! Without ExB shearing effect
  call txmmm95(dNsdrho,dTsdrho,dQdrho,0.d0,gamma)

  write(6,*) "MMM95: Growth rates by Weiland14 model"
  do NR = 0, NRMAX
     write(6,'(F15.7,3ES15.7)') Rho(NR),gamma_exb(NR),gamma(NR),wexb(NR)
  end do

  write(6,*)

  write(6,*) "Linear growth rate (Newman, Rogister, Candy)"
  do NR = 1, NRMAX
     epsN = Var(NR,2)%n / abs(dNsdrho(NR,2)) / RR
     write(6,'(F9.6,5ES14.6)') Rho(NR),epsN*S(NR)/Q(NR),gamITG(NR,1),gamITG(NR,2),gamITG(NR,3),wexb(NR)
  end do

  deallocate(dQdrho,dTsdrho,dNsdrho,gamma,gamma_exb)

end subroutine ITG_growthrate
