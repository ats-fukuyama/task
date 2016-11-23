subroutine txmmm95(dNedrho,dNidrho,dTedrho,dTidrho,dQdrho,gamma)
  use tx_commons, only : NRMAX, R, RR, achg, PNbV, Zeff, Var, Q, &
       &                 BphV, amas, Chie, Chii, De, rMue, rMui, FSANOM, &
       &                 FSCBSH, FSPCLD, FSPCLC, rKeV, AEE, NRA, BthV, wexb
!  use tx_interface, only : dfdx
  implicit none

  ! 0:NRMAX => 1:NRMAX+1 if dimension(*) used
  real(8), dimension(:), intent(in) :: dNedrho, dNidrho, dTedrho, dTidrho, dQdrho
  real(8), dimension(0:NRMAX), intent(out), optional :: gamma

  integer(4), parameter :: mxmode = 12 ! at least 5
  integer(4) :: N, NR, igamma
  integer(4) :: matdim, npoints, nprout, lprint, nerr, lsuper, lreset
  integer(4), dimension(:), allocatable :: lswitch
  real(4), dimension(:), allocatable :: &
       & cswitch, fig, frb, fkb, &
       & zrmajor, zelong, zdense, zdensh, zdensimp, zdensfe, zxzeff, zavezimp, zmassimp, zmasshyd, &
       & zaimass, zwexbs, zgrdne, zgrdni, zgrdnh, zgrdnz, zgrdte, zgrdti, zgrdq, zthiig, zthdig, &
       & ztheig, zthzig, zthirb, zthdrb, ztherb, zthzrb, zthikb, zthdkb, zthekb, zthzkb
  real(4), dimension(:,:), allocatable :: zgamma, zomega, zvelthi, zvflux
  real(4), dimension(:,:,:), allocatable :: zdifthi
  real(8) :: factor_bohm, DeL, cap_val, PAL, PZL

  N = NRMAX ! used for abbreviation

  if(present(gamma)) then
     igamma = 1
  else
     igamma = 0
  end if

  allocate(lswitch(1:8))
  allocate(cswitch(1:25),fig(1:4),frb(1:4),fkb(1:4))
  allocate(zrmajor(0:N),zelong(0:N),zdense(0:N), zdensh(0:N),zdensimp(0:N),zdensfe(0:N),&
       &   zxzeff(0:N),zavezimp(0:N), &
       &   zmassimp(0:N),zmasshyd(0:N),zaimass(0:N),zwexbs(0:N),zgrdne(0:N), &
       &   zgrdni(0:N),zgrdnh(0:N),zgrdnz(0:N),zgrdte(0:N),zgrdti(0:N),zgrdq(0:N), &
       &   zthiig(0:N),zthdig(0:N),ztheig(0:N),zthzig(0:N),zthirb(0:N),zthdrb(0:N), &
       &   ztherb(0:N),zthzrb(0:N),zthikb(0:N),zthdkb(0:N),zthekb(0:N),zthzkb(0:N))
  allocate(zgamma(1:mxmode,0:N),zomega(1:mxmode,0:N),zdifthi(1:mxmode,1:mxmode,0:N), &
       &   zvelthi(1:mxmode,0:N),zvflux(1:mxmode,0:N))

  !     *** Dummy high Z impurity ***
  !         No impurities induce no ITG growth rate.

  PAL = 56.D0 ! Fe
  PZL = 18.D0 ! Fe

  zrmajor(0:N)  = real(RR)
  zelong(0:N)   = 1.0

  zdense(0:N)   = real(Var(0:N,1)%n) * 1.e20
  zdensh(0:N)   = real((PZL*achg(2)-Zeff)/(achg(2)*(PZL-achg(2)))*Var(0:N,2)%n) * 1.e20
  zdensimp(0:N) = real((Zeff-achg(2)**2)/(PZL*(PZL-achg(2)))*Var(0:N,2)%n) * 1.e20
  zdensfe(0:N)  = real(achg(2)*PNbV(0:N)) * 1.e20

  zxzeff(0:N)   = real(Zeff)

  zavezimp(0:N) = real(PZL)
  zmassimp(0:N) = real(PAL)
  zmasshyd(0:N) = real(amas(2))
  zaimass(0:N)  = ( zmasshyd(0:N) * zdensh(0:N) + zmassimp(0:N) * zdensimp(0:N) ) &
       &        / ( zdensh(0:N) + zdensimp(0:N) )

  if(FSCBSH == 0.d0) then
     zwexbs(0:N)   = 0.0
  else
     zwexbs(0:N)   = real(FSCBSH * wexb(0:N))
  end if

  ! Electron density gradient
  zgrdne(0:N)   = - zrmajor(0:N) * real((dNedrho(1:N+1)) / Var(0:N,1)%n)
  ! Thermal ion density gradient (including all thermal ions)
  zgrdni(0:N)   = - zrmajor(0:N) * real((dNidrho(1:N+1)) / Var(0:N,2)%n) ! Same scale length
  ! Hydrogenic ion density gradient (excluding impurities)
  zgrdnh(0:N)   = - zrmajor(0:N) * real((dNidrho(1:N+1)) / Var(0:N,2)%n) ! Same scale length
  ! Impurity ion density gradient
  zgrdnz(0:N)   = - zrmajor(0:N) * real((dNidrho(1:N+1)) / Var(0:N,2)%n) ! Same scale length
  ! Electron temperature gradient
  zgrdte(0:N)   = - zrmajor(0:N) * real((dTedrho(1:N+1)) / Var(0:N,1)%T)
  ! Ion temperature gradient
  zgrdti(0:N)   = - zrmajor(0:N) * real((dTidrho(1:N+1)) / Var(0:N,2)%T)
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
     &   real(R),  zrmajor,   zelong &
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
     Chii(0:N) = FSANOM(3) * (  dble(zthiig(0:N)) &  ! ITG mode
          &                   + dble(zthirb(0:N)) &  ! Resistive Ballooning mode
          &                   + dble(zthikb(0:N)))   ! Kinetic Ballooning mode
!!$     do nr=0,n
!!$        write(6,*) rho(nr),zthiig(nr),zthirb(nr),zthikb(nr)
!!$     end do

     !  Electron thermal diffusivity
     Chie(0:N) = FSANOM(3) * (  dble(ztheig(0:N)) &
          &                   + dble(ztherb(0:N)) &
          &                   + dble(zthekb(0:N)))
  end if

  !  Momentum viscosity regarded as the same as thermal diffusivity
  if(FSANOM(2) /= 0.d0 .and. igamma == 0) then
     rMue(0:N) = FSANOM(2) * Chie(0:N)
     rMui(0:N) = FSANOM(2) * Chii(0:N)
  end if

  !  Particle diffusivity
  if(FSANOM(1) /= 0.d0 .and. igamma == 0) then
     De  (0:N) = FSANOM(1) * (  dble(zthdig(0:N)) &
          &                   + dble(zthdrb(0:N)) &
          &                   + dble(zthdkb(0:N)))
  end if

  if(igamma == 1) then
     ! Most unstable mode at each radial location
     do NR = 0, NRMAX
        gamma(NR) = dble(maxval(zgamma(1:mxmode,NR)))
     end do
  end if

  deallocate(lswitch,cswitch,fig,frb,fkb, &
       &     zrmajor,zelong,zdense,zdensh,zdensimp,zdensfe,zxzeff,zavezimp,zmassimp,zmasshyd, &
       &     zaimass,zwexbs,zgrdne,zgrdni,zgrdnh,zgrdnz,zgrdte,zgrdti,zgrdq,zthiig,zthdig, &
       &     ztheig,zthzig,zthirb,zthdrb,ztherb,zthzrb,zthikb,zthdkb,zthekb,zthzkb, &
       &     zgamma,zomega,zdifthi,zvelthi,zvflux)

  ! Outside the separatrix

  if(FSPCLD == 0.d0) then
     do NR = NRA, NRMAX
        factor_bohm = De(NR) / (Var(NRA,1)%T * rKeV &
             &      / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
        De(NR) = factor_bohm * Var(NR,1)%T * rKeV / (16.D0 * AEE * SQRT(BphV(NR)**2+BthV(NR)**2))
     end do
  else
     do NR = NRA + 1, NRMAX
        De(NR) = De(NRA)
     end do
  end if
     
  if(FSPCLC == 0.d0) then
     do NR = NRA, NRMAX
        factor_bohm = Chie(NR) / (Var(NRA,1)%T * rKeV &
             &      / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
        DeL = factor_bohm * Var(NR,1)%T * rKeV / (16.D0 * AEE * SQRT(BphV(NR)**2+BthV(NR)**2))
        Chie(NR) = DeL
        rMue(NR) = DeL
        factor_bohm = Chii(NR) / (Var(NRA,1)%T * rKeV &
             &      / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
        DeL = factor_bohm * Var(NR,1)%T * rKeV / (16.D0 * AEE * SQRT(BphV(NR)**2+BthV(NR)**2))
        Chii(NR) = DeL
        rMui(NR) = DeL
     end do
  else
     do NR = NRA + 1, NRMAX
        Chie(NR) = Chie(NRA)
        rMue(NR) = rMue(NRA)
        Chii(NR) = Chii(NRA)
        rMui(NR) = rMui(NRA)
     end do
  end if

  ! Limitation of diffusivities

  cap_val = 5.d1
  where(Chie > cap_val) Chie = cap_val
  where(rMue > cap_val) rMue = cap_val
  where(Chii > cap_val) Chie = cap_val
  where(rMui > cap_val) rMui = cap_val

  ! Minimum diffusivities for numerical purpose

  cap_val = 0.1d0
  where(Chie <= 0.d0) Chie = cap_val
  where(rMue <= 0.d0) rMue = cap_val
  where(Chii <= 0.d0) Chie = cap_val
  where(rMui <= 0.d0) rMui = cap_val

end subroutine txmmm95


subroutine ITG_growthrate
  use tx_commons, only : NRMAX, Q, Var, Rho, vv, FSCBSH, &
       &                 wexb, gamITG, RR, S, vro
  use tx_interface, only : dfdx, txmmm95
  implicit none

  integer(4) :: N, NR
  real(8)    :: FSCBSH_STR, epsN
  real(8), dimension(:), allocatable :: dQdrho, dTedrho, dTidrho, dNedrho, dNidrho
  real(8), dimension(:), allocatable :: gamma, gamma_exb

  N = NRMAX

  allocate(dQdrho(0:N),dTedrho(0:N),dTidrho(0:N),dNedrho(0:N),dNidrho(0:N))
  allocate(gamma(0:N),gamma_exb(0:N))

  dQdrho (0:N) = vro(0:N) * dfdx(vv,Q          ,N,0)
  dTedrho(0:N) = vro(0:N) * dfdx(vv,Var(0:N,1)%T,N,0)
  dTidrho(0:N) = vro(0:N) * dfdx(vv,Var(0:N,2)%T,N,0)
  dNedrho(0:N) = vro(0:N) * dfdx(vv,Var(0:N,1)%n,N,0)
  dNidrho(0:N) = vro(0:N) * dfdx(vv,Var(0:N,2)%n,N,0)

  FSCBSH_STR = FSCBSH
  ! With ExB shearing effect
  FSCBSH = 1.d0
  call txmmm95(dNedrho,dNidrho,dTedrho,dTidrho,dQdrho,gamma_exb)

  ! Without ExB shearing effect
  FSCBSH = 0.d0
  call txmmm95(dNedrho,dNidrho,dTedrho,dTidrho,dQdrho,gamma)

  write(6,*) "MMM95: Growth rates by Weiland14 model"
  do NR = 0, NRMAX
     write(6,'(F15.7,1P3E15.7)') Rho(NR),gamma_exb(NR),gamma(NR),wexb(NR)
  end do

  write(6,*)

  write(6,*) "Linear growth rate (Newman, Rogister, Candy)"
  do NR = 1, NRMAX
     epsN = Var(NR,2)%n / abs(dNidrho(NR)) / RR
     write(6,'(F9.6,1P5E14.6)') Rho(NR),epsN*S(NR)/Q(NR),gamITG(NR,1),gamITG(NR,2),gamITG(NR,3),wexb(NR)
  end do

  FSCBSH = FSCBSH_STR

  deallocate(dQdrho,dTedrho,dTidrho,dNedrho,dNidrho,gamma,gamma_exb)

end subroutine ITG_growthrate
