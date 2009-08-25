subroutine txmmm95(dNedr,dNidr,dTedr,dTidr,dQdr)
  use tx_commons, only : NRMAX, R, RR, PNeV, PNiV, PZ, PNbV, Zeff, PTeV, PTiV, Q, &
       &                 BphV, PA, Chie, Chii, De, rMue, rMui, FSANOM, ErV, BphV, &
       &                 PSI, FSCBSH, FSPCLD, FSPCLC, rKeV, AEE, NRA, BthV, wexb
  use tx_interface, only : dfdx
  implicit none

  ! 0:NRMAX => 1:NRMAX+1 if dimension(*) used
  real(8), dimension(*), intent(in) :: dNedr, dNidr, dTedr, dTidr, dQdr

  integer(4), parameter :: mxmode = 12 ! at least 5
  integer(4) :: N, matdim, npoints, nprout, lprint, nerr, lsuper, lreset, NR
  integer(4), dimension(:), allocatable :: lswitch, cswitch
  real(4), dimension(:), allocatable :: &
       & fig, frb, fkb, &
       & zrmajor, zelong, zdensimp, zxzeff, zavezimp, zmassimp, zmasshyd, zaimass, &
       & zwexbs, zgrdne, zgrdni, zgrdnh, zgrdnz, zgrdte, zgrdti, zgrdq, zthiig, zthdig, &
       & ztheig, zthzig, zthirb, zthdrb, ztherb, zthzrb, zthikb, zthdkb, zthekb, zthzkb
  real(4), dimension(:,:), allocatable :: zgamma, zomega, zvelthi, zvflux
  real(4), dimension(:,:,:), allocatable :: zdifthi
  real(8), dimension(:), allocatable :: ErBph, dErBph, qr, dqr
  real(8) :: factor_bohm, DeL

  N = NRMAX ! used for abbreviation

  allocate(lswitch(1:8),cswitch(1:25))
  allocate(fig(1:4),frb(1:4),fkb(1:4))
  allocate(zrmajor(0:N),zelong(0:N),zdensimp(0:N),zxzeff(0:N),zavezimp(0:N), &
       &   zmassimp(0:N),zmasshyd(0:N),zaimass(0:N),zwexbs(0:N),zgrdne(0:N), &
       &   zgrdni(0:N),zgrdnh(0:N),zgrdnz(0:N),zgrdte(0:N),zgrdti(0:N),zgrdq(0:N), &
       &   zthiig(0:N),zthdig(0:N),ztheig(0:N),zthzig(0:N),zthirb(0:N),zthdrb(0:N), &
       &   ztherb(0:N),zthzrb(0:N),zthikb(0:N),zthdkb(0:N),zthekb(0:N),zthzkb(0:N))
  allocate(zgamma(1:mxmode,0:N),zomega(1:mxmode,0:N),zdifthi(1:mxmode,1:mxmode,0:N), &
       &   zvelthi(1:mxmode,0:N),zvflux(1:mxmode,0:N))
  allocate(ErBph(0:N),dErBph(0:N),qr(0:N),dqr(0:N))

  zrmajor(0:N)  = real(RR)
  zelong(0:N)   = 1.0
  zdensimp(0:N) = 0.0
  zxzeff(0:N)   = real(Zeff)
  zavezimp(0:N) = 0.0
  zmassimp(0:N) = 0.0
  zmasshyd(0:N) = real(PA)
  zaimass(0:N)  = real(PA)

  ! ExB shearing rate = r/q d/dr(q v_E/r)
  !                   = r/q [q/r d/dr(v_E) + v_E d/dr(q/r)], where v_E = - ErV / BphV
  if(FSCBSH == 0.d0) then
     wexb(0:N)     = 0.D0
     zwexbs(0:N)   = 0.0
  else
     ErBph(0:N)    = ErV(0:N) / BphV(0:N)
     qr(0)         = 0.d0 ! owing to l'Hopital's rule
     qr(1:N)       = Q(1:N) / R(1:N)
     dErBph(0:N)   = 2.D0 * R(0:NRMAX) * dfdx(PSI,ErBph,NRMAX,0)
     dqr(0:N)      = 2.D0 * R(0:NRMAX) * dfdx(PSI,qr,NRMAX,0)
     wexb(0:N)     = FSCBSH * (- dErBph(0:N) - R(0:N)*ErV(0:N)/(Q(0:N)*BphV(0:N))*dqr(0:N))
     zwexbs(0:N)   = real(wexb)
  end if

  zgrdne(0:N)   = - zrmajor(0:N) * real(dNedr(1:N+1) / PNeV(0:N))
  zgrdni(0:N)   = - zrmajor(0:N) * real(dNidr(1:N+1) / PNiV(0:N))
  zgrdnh(0:N)   = zgrdni(0:N)
  zgrdnz(0:N)   = 0.d0
  zgrdte(0:N)   = - zrmajor(0:N) * real(dTedr(1:N+1) / PTeV(0:N))
  zgrdti(0:N)   = - zrmajor(0:N) * real(dTidr(1:N+1) / PTiV(0:N))
  zgrdq(0:N)    = - zrmajor(0:N) * real(dQdr(1:N+1)  / Q(0:N))

  matdim  = mxmode
  npoints = N + 1
  nprout  = 6       ! output unit number for long printout

  lprint  = 0       ! controls amount of printout
  lsuper  = 0       ! used for supershots
  lreset  = 0       ! use default values for cswitch and lswitch

  call mmm95( &
     ! Input:
     &   real(R),  zrmajor,   zelong &
     & , real(PNeV)*1.e20,   real(PNiV)*1.e20,    zdensimp,  real(PZ*PNbV) &
     & , zxzeff,   real(PTeV),    real(PTiV),    real(Q),       real(BphV) &
     & , zavezimp, zmassimp, zmasshyd, zaimass,  zwexbs &
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

  if(FSANOM(3) /= 0.d0) then
     !  Ion thermal diffusivity
     Chii(0:N) = FSANOM(3) * (  dble(zthiig(0:N)) &  ! ITG mode
          &                   + dble(zthirb(0:N)) &  ! Resistive Ballooning mode
          &                   + dble(zthikb(0:N)))   ! Kinetic Ballooning mode

     !  Electron thermal diffusivity
     Chie(0:N) = FSANOM(3) * (  dble(ztheig(0:N)) &
          &                   + dble(ztherb(0:N)) &
          &                   + dble(zthekb(0:N)))
  end if

  !  Momentum viscosity regarded as the same as thermal diffusivity
  if(FSANOM(2) /= 0.d0) then
     rMue(0:N) = FSANOM(2) * Chie(0:N)
     rMui(0:N) = FSANOM(2) * Chii(0:N)
  end if

  !  Particle diffusivity
  if(FSANOM(1) /= 0.d0) then
     De  (0:N) = FSANOM(1) * (  dble(zthdig(0:N)) &
          &                   + dble(zthdrb(0:N)) &
          &                   + dble(zthdkb(0:N)))
  end if

  deallocate(lswitch,cswitch,ErBph,qr,dErBph,dqr,fig,frb,fkb, &
       &     zrmajor,zelong,zdensimp,zxzeff,zavezimp,zmassimp,zmasshyd,zaimass,zwexbs, &
       &     zgrdne,zgrdni,zgrdnh,zgrdnz,zgrdte,zgrdti,zgrdq,zthiig,zthdig,ztheig,zthzig, &
       &     zthirb,zthdrb,ztherb,zthzrb,zthikb,zthdkb,zthekb,zthzkb, &
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

end subroutine txmmm95
