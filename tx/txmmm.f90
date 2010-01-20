subroutine txmmm95(dNedr,dNidr,dTedr,dTidr,dQdr,gamma)
  use tx_commons, only : NRMAX, R, RR, PNeV, PNiV, PZ, PNbV, Zeff, PTeV, PTiV, Q, &
       &                 BphV, PA, Chie, Chii, De, rMue, rMui, FSANOM, &
       &                 FSCBSH, FSPCLD, FSPCLC, rKeV, AEE, NRA, BthV, wexb, Rho
!  use tx_interface, only : dfdx
  implicit none

  ! 0:NRMAX => 1:NRMAX+1 if dimension(*) used
  real(8), dimension(*), intent(in) :: dNedr, dNidr, dTedr, dTidr, dQdr
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

  zdense(0:N)   = real(PNeV(0:N)) * 1.e20
  zdensh(0:N)   = real((PZL*PZ-Zeff)/(PZ*(PZL-PZ))*PNiV(0:N)) * 1.e20
  zdensimp(0:N) = real((Zeff-PZ**2)/(PZL*(PZL-PZ))*PNiV(0:N)) * 1.e20
  zdensfe(0:N)  = real(PZ*PNbV(0:N)) * 1.e20

  zxzeff(0:N)   = real(Zeff)

  zavezimp(0:N) = real(PZL)
  zmassimp(0:N) = real(PAL)
  zmasshyd(0:N) = real(PA)
  zaimass(0:N)  = ( zmasshyd(0:N) * zdensh(0:N) + zmassimp(0:N) * zdensimp(0:N) ) &
       &        / ( zdensh(0:N) + zdensimp(0:N) )

  if(FSCBSH == 0.d0) then
     zwexbs(0:N)   = 0.0
  else
     zwexbs(0:N)   = real(FSCBSH * wexb(0:N))
  end if

  ! Electron density gradient
  zgrdne(0:N)   = - zrmajor(0:N) * real(dNedr(1:N+1) / PNeV(0:N))
  ! Thermal ion density gradient (including all thermal ions)
  zgrdni(0:N)   = - zrmajor(0:N) * real(dNidr(1:N+1) / PNiV(0:N)) ! Same scale length
  ! Hydrogenic ion density gradient (excluding impurities)
  zgrdnh(0:N)   = - zrmajor(0:N) * real(dNidr(1:N+1) / PNiV(0:N)) ! Same scale length
  ! Impurity ion density gradient
  zgrdnz(0:N)   = - zrmajor(0:N) * real(dNidr(1:N+1) / PNiV(0:N)) ! Same scale length
  ! Electron temperature gradient
  zgrdte(0:N)   = - zrmajor(0:N) * real(dTedr(1:N+1) / PTeV(0:N))
  ! Ion temperature gradient
  zgrdti(0:N)   = - zrmajor(0:N) * real(dTidr(1:N+1) / PTiV(0:N))
  ! Safety factor gradient (~ magnetic shear)
  zgrdq(0:N)    =   zrmajor(0:N) * real(dQdr(1:N+1)  / Q(0:N))

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
     & , zxzeff,   real(PTeV),real(PTiV),real(Q),real(BphV) &
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
        factor_bohm = De(NR) / (PTeV(NRA) * rKeV &
             &      / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
        De(NR) = factor_bohm * PTeV(NR) * rKeV / (16.D0 * AEE * SQRT(BphV(NR)**2+BthV(NR)**2))
     end do
  else
     do NR = NRA + 1, NRMAX
        De(NR) = De(NRA)
     end do
  end if
     
  if(FSPCLC == 0.d0) then
     do NR = NRA, NRMAX
        factor_bohm = Chie(NR) / (PTeV(NRA) * rKeV &
             &      / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
        DeL = factor_bohm * PTeV(NR) * rKeV / (16.D0 * AEE * SQRT(BphV(NR)**2+BthV(NR)**2))
        Chie(NR) = DeL
        rMue(NR) = DeL
        factor_bohm = Chii(NR) / (PTeV(NRA) * rKeV &
             &      / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
        DeL = factor_bohm * PTeV(NR) * rKeV / (16.D0 * AEE * SQRT(BphV(NR)**2+BthV(NR)**2))
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
  use tx_commons, only : NRMAX, Q, PTeV, PTiV, PNeV, PNiV, Rho, R, PSI, FSCBSH, &
       &                 wexb, Ys, gamITG, RR, S
  use tx_interface, only : dfdx, txmmm95
  implicit none

  integer(4) :: N, NR
  real(8)    :: FSCBSH_STR, epsN
  real(8), dimension(:), allocatable :: dQdr, dTedr, dTidr, dNedr, dNidr
  real(8), dimension(:), allocatable :: gamma, gamma_exb

  N = NRMAX

  allocate(dQdr(0:N),dTedr(0:N),dTidr(0:N),dNedr(0:N),dNidr(0:N))
  allocate(gamma(0:N),gamma_exb(0:N))

  dQdr (0:N) = 2.D0 * R(0:N) * dfdx(PSI,Q    ,N,0)
  dTedr(0:N) = 2.D0 * R(0:N) * dfdx(PSI,PTeV ,N,0)
  dTidr(0:N) = 2.D0 * R(0:N) * dfdx(PSI,PTiV ,N,0)
  dNedr(0:N) = 2.D0 * R(0:N) * dfdx(PSI,PNeV ,N,0)
  dNidr(0:N) = 2.D0 * R(0:N) * dfdx(PSI,PNiV ,N,0)

  FSCBSH_STR = FSCBSH
  ! With ExB shearing effect
  FSCBSH = 1.d0
  call txmmm95(dNedr,dNidr,dTedr,dTidr,dQdr,gamma_exb)

  ! Without ExB shearing effect
  FSCBSH = 0.d0
  call txmmm95(dNedr,dNidr,dTedr,dTidr,dQdr,gamma)

  write(6,*) "MMM95: Growth rates by Weiland14 model"
  do NR = 0, NRMAX
     write(6,'(F15.7,1P3E15.7)') Rho(NR),gamma_exb(NR),gamma(NR),wexb(NR)
  end do

  write(6,*)

  write(6,*) "Linear stability theory parameter (valid almost 0.2 < rho < 0.8)"
  do NR = 0, NRMAX
     write(6,'(F15.7,1PE15.7)') Rho(NR),Ys(NR)
  end do

  write(6,*) 

  write(6,*) "Linear growth rate (Newman, Rogister, Candy)"
  do NR = 1, NRMAX
     epsN = PNiV(NR) / abs(dNidr(NR)) / RR
     write(6,'(F9.6,1P5E14.6)') Rho(NR),epsN*S(NR)/Q(NR),gamITG(NR,1),gamITG(NR,2),gamITG(NR,3),wexb(NR)
  end do

  FSCBSH = FSCBSH_STR

  deallocate(dQdr,dTedr,dTidr,dNedr,dNidr,gamma,gamma_exb)

end subroutine ITG_growthrate
