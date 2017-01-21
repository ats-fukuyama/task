!mmm95.f
c
c  Revision History
c  ----------------
c       date            Description
c
c   18-Dec-2001         Resistive ballooning diffusion = abs ( ... )
c
c   13-Aug-2001         Fixed normalization of gamma and omega
c
c   26-Mar-1999         Added statement at beginning of do 300 loop
c                       to avoid computation of fluxes at
c                       rminor(jz) .lt. 1.e-4 * rmajor(jz)
c
c   09-Mar-1999         Revamped comments and included changes as
c                       suggested by D. McCune, module reviewer for
c                       NTCC.
c

c@mmm95.tex
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine mmm95 (
     &   rminor,  rmajor,   elong
     & , dense,   densh,    densimp,  densfe
     & , xzeff,   tekev,    tikev,    q,       btor
     & , avezimp, amassimp, amasshyd, aimass,  wexbs
     & , grdne,   grdni,    grdnh,    grdnz,   grdte,   grdti,  grdq
     & , thiig,   thdig,    theig,    thzig
     & , thirb,   thdrb,    therb,    thzrb
     & , thikb,   thdkb,    thekb,    thzkb
     & , gamma,   omega,    difthi,   velthi,  vflux
     & , matdim,  npoints,  nprout,   lprint,  nerr
     & , lsuper,  lreset,   lswitch,  cswitch, fig,    frb,     fkb)
c
c
c    All the following 1-D arrays are assumed to be defined on flux
c    surfaces called zone boundaries where the transport fluxes are
c    to be computed.  The number of flux surfaces is given by npoints
c    (see below).  For example, if you want to compute the transport
c    on only one flux surface, set npoints = 1.
c
c    Note, if the mmm95 module is used in another code, the minimum
c    dimension required for cswitch is 23, the minimum dimension for
c    lswitch is 5 and fig, frb, and fkb must all be dimensioned to a
c    minimum of 4.

c  Input arrays:
c  -------------
c
c  rminor(jz)   = minor radius (half-width) of zone boundary [m]
c  rmajor(jz)   = major radius to geometric center of zone bndry [m]
c  elong(jz)    = local elongation of zone boundary
c
c  dense(jz)    = electron density [m^-3]
c  densh(jz)    = sum over thermal hydrogenic ion densities [m^-3]
c  densimp(jz)  = sum over impurity ion densities [m^-3]
c  densfe(jz)   = electron density from fast (non-thermal) ions [m^-3]
c
c  xzeff(jz)    = Z_eff
c  tekev(jz)    = T_e (electron temperature) [keV]
c  tikev(jz)    = T_i (temperature of thermal ions) [keV]
c  q(jz)        = magnetic q-value
c  btor(jz)     = ( R B_tor ) / rmajor(jz)  [tesla]
c
c  avezimp(jz)  = average density weighted charge of impurities
c               = ( sum_imp n_imp Z_imp ) / ( sum_imp n_imp ) where
c                 sum_imp = sum over impurity ions with charge state Z_imp
c
c  amassimp(jz) = average density weighted atomic mass of impurities
c               = ( sum_imp n_imp M_imp ) / ( sum_imp n_imp ) where
c                 sum_imp = sum over impurity ions, each with mass M_imp
c
c  amasshyd(jz) = average density weighted atomic mass of hydrogen ions
c               = ( sum_hyd n_hyd M_hyd ) / ( sum_hyd n_hyd ) where
c                 sum_hyd = sum over hydrogenic ions, each with mass M_hyd
c
c  aimass(jz)   = mean atomic mass of thermal ions [AMU]
c               = ( sum_i n_i M_i ) / ( sum_i n_i ) where
c                 sum_i = sum over all ions, each with mass M_i
c
c  wexbs(jz)    = ExB shearing rate in [rad/s].  See  K.H. Burrell,
c                 "Effects of {ExB} velocity shear and magnetic shear
c                 on turbulence and transport in magnetic confinement
c                 devices", Phys. of Plasmas, 4, 1499 (1997).
c
c    All of the following normalized gradients are at zone boundaries.
c    r = half-width, R = major radius to center of flux surface
c
c  grdne(jz) = -R ( d n_e / d r ) / n_e
c  grdni(jz) = -R ( d n_i / d r ) / n_i
c  grdnh(jz) = -R ( d n_h / d r ) / n_h
c  grdnz(jz) = -R ( d Z n_Z / d r ) / ( Z n_Z )
c  grdte(jz) = -R ( d T_e / d r ) / T_e
c  grdti(jz) = -R ( d T_i / d r ) / T_i
c  grdq (jz) =  R ( d q   / d r ) / q    related to magnetic shear
c
c  where:
c    n_i     = thermal ion density (sum over hydrogenic and impurity)
c    n_h     = thermal hydrogenic density (sum over hydrogenic species)
c    n_Z     = thermal impurity density,  Z = average impurity charge
c                      sumed over all impurities
c
c  Output:
c  -------
c
c    The following effective diffusivities represent contributions
c    to the total diffusivity matrix (difthi and velthi given below)
c    from each of the models that contribute to the Multi-Mode model.
c    Generally, these arrays are used for diagnostic output only.
c
c  thiig(jz) = ion thermal diffusivity from the Weiland model
c  thdig(jz) = hydrogenic ion diffusivity from the Weiland model
c  theig(jz) = elelctron thermal diffusivity from the Weiland model
c  thzig(jz) = impurity ion diffusivity from the Weiland model
c
c  thirb(jz) = ion thermal diffusivity from resistive ballooning modes
c  thdrb(jz) = hydrogenic ion diffusivity from resistive ballooning modes
c  therb(jz) = elelctron thermal diffusivity from resistive ballooning modes
c  thzrb(jz) = impurity ion diffusivity from resistive ballooning modes
c
c  thikb(jz) = ion thermal diffusivity from kinetic ballooning modes
c  thdkb(jz) = hydrogenic ion diffusivity from kinetic ballooning modes
c  thekb(jz) = elelctron thermal diffusivity from kinetic ballooning modes
c  thzkb(jz) = impurity ion diffusivity from kinetic ballooning modes
c
c    The following are growth rates and mode frequencies from the
c    Weiland model for drift modes such as ITG and TEM.
c    These arrays are intended for diagnostic output.
c
c  gamma(jm,jz) = growth rate for mode jm at point jz ( 1/sec )
c  omega(jm,jz) = frequency for mode jm at point jz ( radians/sec )
c
c    All of the transport coefficients are given in the following two
c    matricies for diffusion difthi and convection velthi in MKS units.
c    See the LaTeX documentation for difthi and velthi just below.
c
c    NOTE:  difthi and velthi include all of the anomalous transport.
c    There are no additional contributions to the heat fluxs from
c    charged particle convection.
c
c  difthi(j1,j2,jz) = full matrix of anomalous transport diffusivities
c  velthi(j1,jz)    = convective velocities
c  vflux(j1,jz)     = flux matrix
c
c  Input integers:
c  ---------------
c
c  matdim  = first and second dimension of transport matricies
c            difthi(j1,j2,jz) and the first dimension of
c            velthi(j1,jz), vflux(j1,jz), gamma(j1,jz), and omega(j1,jz).
c            matdim must be at least 5
c
c  npoints = number of values of jz in all of the above arrays
c
c  nprout  = output unit number for long printout
c
c
c  Input switches
c  --------------
c
c  lprint      controls the amount of printout (0 => no printout)
c              higher values yield more diagnostic output
c
c  lsuper   = 0 for simulations of all other discharges
c           > 0 for supershot simulations; substantially reduces
c               contribution from kinetic ballooning mode
c
c  lreset  = 0 to use only internal settings for lswitch, cswitch
c              and for the coefficients fig, frb, and fkb that control
c              the contributions form the various instability modes
c
c    Note that when lreset = 0, the values of the switches and
c    coefficients in the argument list are ignored and all the
c    switches and coefficients are set internally.
c
c    WARNING:  use lreset > 0 only if you want to pass all the switches
c              lswitch, cswitch, fig, frb, and fkb through the
c              argument list.
c
c    WARNING:  NTCC users desiring to use anything other than lreset = 0
c              should consult with the mmm95 code authors first.
c
c
c  Output Integer
c  --------------
c
c  nerr        status code returned; 0 = OK; .ne. 0 indicates error
c
c
c  Internal control variables:
c  ---------------------------
c
c  lswitch(j), j=1,8   integer control variables:
c
c  cswitch(j), j=1,25   general control variables:
c
c  lswitch(1)  controls which version of the Weiland model is used
c                  Default lswitch(1) = 10
c             = 2  2 eqn  Weiland model Hydrogen \eta_i mode only
c             = 4  4 eqn  Weiland model with Hydrogen and trapped electrons
c             = 5  5 eqn  Weiland model with trapped electrons, FLR effects,
c                         and parallel ion motion
c             = 6  6 eqn  Weiland model Hydrogen, trapped electrons,
c                    and one impurity species
c             = 7  7 eqn   Weiland model Hydrogen, trapped electrons,
c                  one impurity species, and collisions
c             = 8  8 eqn  Weiland model Hydrogen, trapped electrons,
c                  one impurity species, collisions, and parallel
c                  ion (hydrogenic) motion
c             = 9  9 eqn  Weiland model Hydrogen, trapped electrons,
c                  one impurity species, collisions, and finite beta
c             = 10 10 eqn Weiland model Hydrogen, trapped electrons,
c                  one impurity species, collisions, parallel
c                  ion (hydrogenic) motion, and finite beta
c             = 11 11 eqn Weiland model Hydrogen, trapped electrons,
c                  one impurity species, collisions, parallel
c                  ion (hydrogenic, impurity) motion, and finite beta
c
c  lswitch(2) = 0  full matrix representation for difthi and velthi
c                  Default lswitch(2) = 2
c             = 1  set diagonal matrix elements difthi and velthi
c             = 2  set diagonal matrix elements = effective diffusivities
c
c  lswitch(3)  controls \kappa scaling
c                  Default lswitch(3) = 0
c             = 0  use \kappa scaling raised to
c                  exponents (cswitch(3) - cswitch(5))
c             = 1  use (1+\kappa^2)/2 instead of \kappa scaling
c
c  lswitch(4) > 0  to replace negative diffusivity with velocity
c                  Default lswitch(4) = 1
c
c  lswitch(5) = 1  to limit magnitude of all normalized gradients
c                  to ( major radius ) / ( ion Larmor radius )
c                  Default lswitch(5) = 1
c
c  cswitch(1)   0.5  minimum value of shear
c  cswitch(2)   3.5  coeff in exponential (fbeta-th) of kinetic ballooning model
c  cswitch(3)  -4.0  exponent of local elongation multiplying drift waves
c  cswitch(4)  -4.0  exponent of local elongation multiplying resistive
c                     ballooning modes
c  cswitch(5)  -4.0  exponent of local elongation multiplying
c                     kinetic balllooning modes
c  cswitch(6)   0.0  k_y \rho_s (= 0.316 if abs(cswitch(6)) < zepslon)
c  cswitch(8)   1.0  coeff of beta_prime_1 in kinetic ballooning mode
c  cswitch(9)  0.15  alpha in diamagnetic stabil. in kinetic ballooning model
c  cswitch(10)  0.0  rel fract of ion thermal diffusivity given to convection
c  cswitch(11)  0.0  rel fract of hydrogen particle diffusivity given to convection
c  cswitch(12)  0.0  rel fract of el thermal diffusivity given to convection
c  cswitch(13)  0.0  rel fract of impurity particle diffusivity given to convection
c  cswitch(14)  1.0  coef of finite beta effect in weiland14 = cetain(20)
c  cswitch(15)  0.0  min value of impurity charge state zimpz
c  cswitch(16)  0.0  coef of fast particle fraction (superthermal ions)
c                    in weiland model -- coef of densfe
c  cswitch(17)  1.0  coeff of k_\parallel (parallel ion motion) in weiland14
c                    = cetain(10)
c  cswitch(18)  0.0  coeff of nuhat (effect of collisions) in weiland14
c                    = cetain(15)
c  cswitch(19)  0.0  coeff for including v_parallel in strong ballooning limit
c                    = cetain(12); cswitch(19) = 1 for inclusion of v_par effect
c  cswitch(20)  0.0  trapping fraction used in weiland14 (when > 0.0)
c                    multiplies electron trapping fraction when < 0.0
c                    no effect when cswitch(20) = 0.0
c  cswitch(21)  1.0  multiplier for wexbs (flow shear rate) in Weiland model
c  cswitch(22)  0.0  ranges from 0.0 to 1.0 adds impurity heat flow to total
c                    ionic heat flow for the weiland model
c  cswitch(23)  0.0  controls finite diff to construct the zgm matrix
c                    = cetain(30)
c
c     contributions to vfluxes and interchanges:
c
c  fig(1)   hydrogen particle transport from ITG (eta_i) mode
c  fig(2)   electron thermal  transport from ITG (eta_i) mode
c  fig(3)   ion      thermal  transport from ITG (eta_i) mode
c  fig(4)   impurity particle transport from ITG (eta_i) mode
c
c  frb(1)   hydrogen particle transport from resistive ballooning mode
c  frb(2)   electron thermal  transport from resistive ballooning mode
c  frb(3)   ion      thermal  transport from resistive ballooning mode
c  frb(4)   impurity particle transport from resistive ballooning mode
c
c  fkb(1)   hydrogen particle transport from kinetic ballooning mode
c  fkb(2)   electron thermal  transport from kinetic ballooning mode
c  fkb(3)   ion      thermal  transport from kinetic ballooning mode
c  fkb(4)   impurity particle transport from kinetic ballooning mode
c
c
c***********************************************************************
c
c-----------------------------------------------------------------------
c
c  Compile this routine and routines that it calls with a compiler
c  option, such as -r8, to convert real to double precision when used on
c  workstations.
c
c-----------------------------------------------------------------------
c
c  External dependencies:
c
c  Call tree: MMM95 calls the following routines
c
c  WEILAND14       - Computes diffusion matrix and convect velocities
c                        for the Weiland transport model
c    WEILAND14FLUX - Calculates fluxes and effective diffusivities
c      TOMSQZ      - Wrapper for QZ algorithm solving Ax = lambda Bx
c         CQZHES   - First step in QZ algorithm
c         CQZVAL   - Second and third step in QZ algorithm
c         CQZVEC   - Fourth step in QZ algorithm
c
c-----------------------------------------------------------------------

      implicit none
c
      integer km, klswitch, kcswitch
c
      integer  matdim,  npoints, nprout,  lprint,   nerr
c
      integer  lsuper,  lreset,  lswitch(*)
c
      parameter ( km = 12, klswitch = 8, kcswitch = 25 )
c
      real
     &   rminor(*),  rmajor(*),   elong(*)
     & , dense(*),   densh(*),    densimp(*),  densfe(*)
     & , xzeff(*),   tekev(*),    tikev(*),    q(*),       btor(*)
     & , avezimp(*), amassimp(*), amasshyd(*), aimass(*),  wexbs(*)
     & , grdne(*),   grdni(*),    grdnh(*),    grdnz(*)
     & , grdte(*),   grdti(*),    grdq(*)
c
      real
     &   thiig(*),   thdig(*),    theig(*),    thzig(*)
     & , thirb(*),   thdrb(*),    therb(*),    thzrb(*)
     & , thikb(*),   thdkb(*),    thekb(*),    thzkb(*)
     & , omega(matdim,*),         gamma(matdim,*)
     & , difthi(matdim,matdim,*), velthi(matdim,*)
     & , vflux(matdim,*)
c
      real     cswitch(*)
c
      real     fig(*),  fkb(*),  frb(*)
c
c..physical constants
c
      real zpi,  zcc,  zcmu0,  zceps0,  zckb,  zcme,  zcmp,  zce
c
c  zpi     = pi
c  zcc     = speed of light                  [m/sec]
c  zcmu0   = vacuum magnetic permeability    [henrys/m]
c  zceps0  = vacuum electrical permittivity  [farads/m]
c  zckb    = energy conversion factor        [Joule/keV]
c  zcme    = electron mass                   [kg]
c  zcmp    = proton mass                     [kg]
c  zce     = electron charge                 [Coulomb]
c
c..computer constants
c
      real  zepslon, zlgeps
c
c  zepslon = machine epsilon [smallest number so that 1.0+zepslon>1.0]
c  zlgeps  = ln ( zepslon )
c
c
c..local variables
c
      integer  jz, j1, j2, jm

      real  zelong,  zelonf,  zai,     zne,     zni,    zte,    zti
     & ,    zq,      zeff,    zgne,    zgni,    zgnh,   zgnz,   zgte
     & ,    zgth,    zgtz,    zshear,  zrmin,   zrmaj,  zbtor,  zep
     & ,    zgyrfi,  zbeta,   zvthe,   zvthi,   zsound, zlog,   zcf
     & ,    znuei,   znueff,  zlari,   zlarpo,  zrhos,  zwn,    znude
     & ,    znuhat,  zgpr,    zscyl,   zsmin,   zshat,  zgmax
c
c.. variables for Weiland model
c
c  iletai(j1) and cetain(j1) are control variables
c
      integer        iletai(32)
c
      real  cetain(32), zomega(km), zgamma(km), zchieff(km)
c
      real           zdfthi(km,km),    zvlthi(km),     zflux(km)
c
      integer        idim,    ieq,     imodes
c
      real  zthte,   zbetae,  znz,     zmass,  zimpz,  ztzte
     & ,    zfnzne,  zmzmh,   zfnsne,  zftrap, zkyrho, zomegde
     & ,    zwexb,   znormd,  znormv
c
c  zexb    = local copy of ExB shearing rate
c  znormd  = factor to convert normalized diffusivities
c  znormv  = factor to convert normalized convective velocities
c
c..local variables for resistive ballooning modes
c
      real  zgyrfe, zlare,   zgddia, zgdp
c
c..local variables for kinetic ballooning modes
c
      real  zbprim, zbcoef1, zbc1,   zelfkb,  zfbthn, zdk
c
c-----------------------------------------------------------------------
c
c..initialize imodes
        imodes  = 0
c..physical constants
c
        zpi     = atan2 ( 0.0, -1.0 )
        zcc     = 2.997925e+8
        zcmu0   = 4.0e-7 * zpi
        zceps0  = 1.0 / ( zcc**2 * zcmu0 )
        zckb    = 1.60210e-16
        zcme    = 9.1091e-31
        zcmp    = 1.67252e-27
        zce     = 1.60210e-19
c
c..computer constants
c
        zepslon = 1.0e-34
        zlgeps  = log ( zepslon )
c
c
c..initialize arrays
c
      do jz = 1, npoints
        thiig(jz)  = 0.
        thdig(jz)  = 0.
        theig(jz)  = 0.
        thzig(jz)  = 0.
        therb(jz)  = 0.
        thirb(jz)  = 0.
        thdkb(jz)  = 0.
        thekb(jz)  = 0.
        thzkb(jz)  = 0.
        thikb(jz)  = 0.
        thdkb(jz)  = 0.
        thekb(jz)  = 0.
        thzkb(jz)  = 0.
      enddo
c
      do jz = 1, npoints
        do j1 = 1, matdim
          velthi(j1,jz) = 0.0
          vflux(j1,jz) = 0.0
          gamma(j1,jz) = 0.0
          omega(j1,jz) = 0.0
          do j2 = 1, matdim
            difthi(j1,j2,jz) = 0.0
          enddo
        enddo
      enddo
c
      nerr = 0
c
c..if lreset < 1, use internal settings for switches and coefficients
c  otherwise, use values passed through the argument list above
c
      if ( lreset .lt. 1 ) then
c
c..initialize switches
c
      do j1=1,kcswitch
        cswitch(j1) = 0.0
      enddo
c
      do j1=1,klswitch
        lswitch(j1) = 0
      enddo
c
c
c  Multi Mode Model in sbrtn THEORY version MMM95
c  for use in the BALDUR transport code
c
      lswitch(1) = 10 ! Weiland ITG model weiland14 (10 eqns, no collisions)
      lswitch(2) = 2  ! use effective diffusivities
      lswitch(3) = 0  ! use kappa instead of (1+\kappa^2)/2
      lswitch(4) = 1  ! replace -ve diffusivity with convective velocity
      lswitch(5) = 1  ! limit gradients by major radius / ion Larmor radius
c
c  misc. parameters for subroutine mmm95
c
      cswitch(1)  =  0.5  ! minimum value of shear
      cswitch(2)  =  3.5  ! coef in exponential (fbeta-th) in kinetic ballooning
      cswitch(3)  = -4.0  ! elongation scaling for drift wave mode
      cswitch(4)  = -4.0  ! elongation scaling for resistive ballooning mode
      cswitch(5)  = -4.0  ! elongation scaling for kinetic ballooning mode
      cswitch(6)  =  0.0  ! k_y \rho_s (= 0.316 if abs(cswitch(6)) < zepslon)
      cswitch(8)  =  1.0  ! coeff of beta_prime_1 in kinetic ballooning mode
      cswitch(9)  = 0.15  ! alpha in diamagnetic stabil. in kinetic balloon model
      cswitch(10) =  0.0  ! relative fraction of ion thermal diffusivity
                          ! given to convection
      cswitch(11) =  0.0  ! relative fract of hydrogen particle diffusivity
                          ! given to convection
      cswitch(12) =  0.0  ! relative fraction of el thermal diffusivity
                          ! given to convection
      cswitch(13) =  0.0  ! relative fract of impurity particle diffusivity
                          ! given to convection
      cswitch(14) =  1.0  ! coef of finite beta effect in weiland14 = cetain(20)
      cswitch(15) =  0.0  ! min value of impurity charge state zimpz
      cswitch(16) =  1.0  ! coef of fast particle fraction (superthermal ions)
                          ! in weiland model -- coef of densfe
      cswitch(17) =  1.0  ! coeff of k_\parallel (parallel ion motion) in
                          ! weiland14 = cetain(10)
      cswitch(18) =  0.0  ! coeff of nuhat in weiland14 = cetain(15)
      cswitch(19) =  0.0  ! coeff for including v_parallel in strong ballooning
                          ! limit = cetain(12); cswitch(19) = 1 for inclusion
                          ! of v_par effect
      cswitch(20) =  0.0  ! trapping fraction used in weiland14 (when > 0.0)
                          ! multiplies electron trapping fraction when < 0.0
                          ! no effect when cswitch(20) = 0.0
      cswitch(21) =  1.0  ! multiplier for wexbs (flow shear rate)
                          ! in Weiland model
      cswitch(22) =  0.0  ! multiplier to impurity heat flux
      cswitch(23) =  0.0  ! controls finite diff to construct the
                          ! zgm matrix = cetain(30)

c  contributions to hydrogenic particle, elec-energy, ion-energy,
c    and impurity ion fluxes

        fig(1) = 0.80
        fig(2) = 0.80
        fig(3) = 0.80
        fig(4) = 0.80
c
        fkb(1) = 1.00
        fkb(2) = 0.65
        fkb(3) = 0.65
        fkb(4) = 1.00
c
        if ( lsuper .gt. 0 ) then
          fkb(1) = 0.045
          fkb(2) = 0.010
          fkb(3) = 0.010
          fkb(4) = 0.045
        endif
c
        frb(1) = 1.00
        frb(2) = 1.00
        frb(3) = 1.00
        frb(4) = 1.00
c
      endif

c
c.. start the main do-loop over the radial index "jz"..........
c
c
      do 300 jz = 1, npoints

c avoid computation of fluxes at
c rminor(jz) .lt. 1.e-4 * rmajor(jz)
c
      if (rminor(jz) .lt. (1.e-4 * rmajor(jz))) go to 300
c
c  transfer common to local variables to compact the notation
c
      zelong = max (zepslon,elong(jz))
      if ( lswitch(3) .eq. 1 ) then
        zelonf = ( 1. + zelong**2 ) / 2.
      else
        zelonf = zelong
      endif
c
      zai    = aimass(jz)
      zne    = dense(jz)
      zni    = densh(jz) + densimp(jz)
      znz    = densimp(jz)
      zte    = tekev(jz)
      zti    = tikev(jz)
      zq     = q(jz)
      zeff   = xzeff(jz)
c
c  normalized gradients
c
      zgne   = grdne(jz)
      zgni   = grdni(jz)
      zgnh   = grdnh(jz)
      zgnz   = grdnz(jz)
      zgte   = grdte(jz)
      zgth   = grdti(jz)
      zgtz   = grdti(jz)

      zrmin  = max( rminor(jz), zepslon )
      zrmaj  = rmajor(jz)
      zshear = grdq(jz) * zrmin / zrmaj
      zbtor  = btor(jz)
c
c  compute inverse aspect ratio
c
      zep    = max( zrmin/zrmaj, zepslon )
c
c
      zgyrfi = zce * zbtor / (zcmp * zai)
      zbeta  = (2. * zcmu0 * zckb / zbtor**2) * (zne * zte + zni * zti)
      zvthe  = sqrt(2. * zckb * zte / zcme)
      zvthi  = sqrt(2. * zckb * zti / (zcmp * zai))
      zsound = sqrt(zckb * zte / (zcmp * zai))
      zlog   = 37.8-log(sqrt(zne) / zte)
      zcf    = (4. * sqrt(zpi) / 3.)
      zcf    = zcf * (zce / (4. * zpi * zceps0))**2
      zcf    = zcf * (zce / zckb) * sqrt( (zce/zcme) * (zce/zckb) )
      znuei  = zcf * sqrt(2.) * zne * zlog * zeff / (zte * sqrt(zte))
c
      znueff = znuei / zep
      zlari  = zvthi / zgyrfi
      zlarpo = max(zlari * zq / zep, zepslon)
      zrhos  = zsound / zgyrfi
      zwn    = 0.3 / zrhos
      znude  = 2 * zwn * zrhos * zsound / zrmaj
      znuhat = znueff / znude
c
c..if lswitch(5) = 1, limit magnitude of normalized gradients
c                    to ( major radius ) / ( ion Larmor radius )
c
      zgmax = zrmaj / zlarpo
c
      if ( lswitch(5) .eq. 1 ) then
c
        zgne = sign ( min ( abs ( zgne ), zgmax ), zgne )
        zgni = sign ( min ( abs ( zgni ), zgmax ), zgni )
        zgnh = sign ( min ( abs ( zgnh ), zgmax ), zgnh )
        zgnz = sign ( min ( abs ( zgnz ), zgmax ), zgnz )
        zgte = sign ( min ( abs ( zgte ), zgmax ), zgte )
        zgth = sign ( min ( abs ( zgth ), zgmax ), zgth )
        zgtz = sign ( min ( abs ( zgtz ), zgmax ), zgtz )
c
      endif
c
c  zgpr = -R ( d p   / d r ) / p    for thermal pressure
c
c  Compute the pressure scale length using smoothed and bounded
c  density and temperature
c
      zgpr = ( zne * zte * ( zgne + zgte )
     &         + zni * zti * ( zgni + zgth ) )
     &         / ( zne * zte + zni * zti )
c
      if ( lswitch(5) .eq. 1 )
     &  zgpr = sign ( min ( abs ( zgpr ), zgmax ), zgpr )
c
c
c
      zscyl=max(abs(zshear),zepslon)
      zsmin=max(cswitch(1),zepslon)
      zshat=max(zsmin,zscyl)
c
c
        do j1=1,32
          iletai(j1) = 0
          cetain(j1) = 0.0
        enddo
c
        thiig(jz) = 0.0
        theig(jz) = 0.0
        thdig(jz) = 0.0
        thzig(jz) = 0.0
c
c..set the number of equations to use in the Weiland model
c
        if ( (lswitch(1) .lt. 2) .or. (lswitch(1) .gt. 11 )) then
          nerr = -10
          return
        elseif (lswitch(1) .eq. 3) then
          nerr = -10
          return
        else
          ieq = lswitch(1)
        endif
c
        cetain(11) = 1.0
c
c.. coefficient of k_parallel for parallel ion motion
c.. cswitch(19) for v_parallel in strong ballooning limit
c.. in 9 eqn model
c
        cetain(10) = cswitch(17)
        cetain(12) = cswitch(19)
        cetain(15) = cswitch(18)
        cetain(20) = cswitch(14)
c
        iletai(10) = 0
c
        idim   = km
c
c  Hydrogen species
c
        zthte  = zti / zte

        zbetae = 2. * zcmu0 * zckb * zne * zte / zbtor**2
c
c  Impurity species (use only impurity species 1 for now)
c  assume T_Z = T_H throughout the plasma here
c
        znz    = densimp(jz)
        zmass  = amassimp(jz)
        zimpz  = avezimp(jz)
        zimpz  = max ( zimpz, cswitch(15) )
c
        ztzte  = zti / zte
        zfnzne = znz / zne
        zmzmh  = zmass / amasshyd(jz)
c
c  superthermal ions
c
c  zfnsne = ratio of superthermal ions to electrons
c  L_ns   = gradient length of superthermal ions
c
        zfnsne = max ( cswitch(16) * densfe(jz) / dense(jz) , 0.0 )
c
        zftrap = sqrt ( 2. * zrmin / ( zrmaj * ( 1. + zrmin / zrmaj )))
        if ( cswitch(20) .gt. zepslon ) zftrap = cswitch(20)
        if ( cswitch(20) .lt. -zepslon )
     &       zftrap = abs(cswitch(20))*zftrap
c
        if ( abs(cswitch(6)) .lt. zepslon ) then
          zkyrho = 0.316
        else
          zkyrho = cswitch(6)
        endif
c
c
c...Define a local copy of normalized ExB shearing rate : pis
c
        zomegde = 2.0 * zkyrho * zsound / zrmaj
c
        zwexb = cswitch(21) * wexbs(jz) / zomegde
c
c
        cetain(30) = cswitch(23)
        iletai(6)  = 0
        if ( lswitch(2) .lt. 1 ) iletai(7) = 1
c
c  if lswitch(2) .lt. 1, compute only the effective diffusivities
c
        iletai(9) = lswitch(2)
c
        call weiland14 (
     &   iletai,   cetain,   lprint,   ieq,      nprout,   zgne
     & , zgnh,     zgnz,     zgte,     zgth,     zgtz,     zthte
     & , ztzte,    zfnzne,   zimpz,    zmzmh,    zfnsne,   zbetae
     & , zftrap,   znuhat,   zq,       zshat,    zkyrho,   zwexb
     & , idim,     zomega,   zgamma,   zdfthi,   zvlthi,   zchieff
     & , zflux,    imodes,   nerr )
c
c  If nerr not equal to 0 an error has occured
c
        if (nerr .ne. 0) return
c
c
c  Growth rates for diagnostic output
c    Note that all frequencies are normalized by \omega_{De}
c      consequently, trapped electron modes rotate in the positive
c      direction (zomega > 0) while eta_i modes have zomega < 0.
c
        jm = 0
        do j1=1,imodes
          if ( zgamma(j1) .gt. zepslon ) then
            jm = jm + 1
            gamma(jm,jz) = zgamma(j1) * zomegde
            omega(jm,jz) = zomega(j1) * zomegde
          endif
        enddo
c
c  conversion factors for diffusion and convective velocity
c
        znormd = zelonf**cswitch(3) *
     &    2.0 * zsound * zrhos**2 / ( zrmaj * zkyrho )
        znormv = zelonf**cswitch(3) *
     &    2.0 * zsound * zrhos**2 / ( zrmaj**2 * zkyrho )
c
c  compute effective diffusivites for diagnostic purposes only
c
        thdig(jz) = fig(1) * znormd * zchieff(2)
        theig(jz) = fig(2) * znormd * zchieff(3)
        thiig(jz) = fig(3) * znormd * zchieff(1)
     &  + fig(3) * znormd * zchieff(5) * cswitch(22) * znz / zni
        thzig(jz) = fig(4) * znormd * zchieff(4)
c
c  start computing the fluxes
c
        vflux(1,jz) = vflux(1,jz) + thiig(jz) * zgth / zrmaj
        vflux(2,jz) = vflux(2,jz) + thdig(jz) * zgnh / zrmaj
        vflux(3,jz) = vflux(3,jz) + theig(jz) * zgte / zrmaj
        vflux(4,jz) = vflux(4,jz) + thzig(jz) * zgnz / zrmaj
c
c  compute diffusivity matrix
c
        do j1=1,matdim
          velthi(j1,jz) = 0.0
          vflux(j1,jz) = 0.0
          do j2=1,matdim
            difthi(j1,j2,jz) = 0.0
          enddo
        enddo
c
c..set difthi and velthi
c
        if ( lswitch(2) .gt. 1 ) then
c
c  diagonal elements of matrix = effective diffusivities
c
          difthi(1,1,jz) = difthi(1,1,jz) + thiig(jz)
          difthi(2,2,jz) = difthi(2,2,jz) + thdig(jz)
          difthi(3,3,jz) = difthi(3,3,jz) + theig(jz)
          difthi(4,4,jz) = difthi(4,4,jz) + thzig(jz)
c
        else
c
c..full matrix form of model
c
          if ( ieq .eq. 2 ) then
            difthi(1,1,jz) = difthi(1,1,jz) +
     &        fig(3) * znormd * zdfthi(1,1)
            velthi(1,jz)   = velthi(1,jz) +
     &        fig(3) * znormv * zvlthi(1)
          elseif ( ieq .eq. 4 ) then
            do j2=1,3
              difthi(1,j2,jz) = difthi(1,j2,jz) +
     &          fig(3) * znormd * zdfthi(1,j2)
              difthi(2,j2,jz) = difthi(2,j2,jz) +
     &          fig(1) * znormd * zdfthi(2,j2)
              difthi(3,j2,jz) = difthi(3,j2,jz) +
     &          fig(2) * znormd * zdfthi(3,j2)
              difthi(4,j2,jz) = difthi(4,j2,jz) +
     &          fig(4) * znormd * zdfthi(2,j2)
            enddo
              velthi(1,jz)    = velthi(1,jz) +
     &          fig(3) * znormv * zvlthi(1)
              velthi(2,jz)    = velthi(2,jz) +
     &          fig(1) * znormv * zvlthi(2)
              velthi(3,jz)    = velthi(3,jz) +
     &          fig(2) * znormv * zvlthi(3)
              velthi(4,jz)    = velthi(4,jz) +
     &          fig(4) * znormv * zvlthi(2)
          else
            do j2=1,4
              difthi(1,j2,jz) = difthi(1,j2,jz) +
     &          fig(3) * znormd * zdfthi(1,j2)
              difthi(2,j2,jz) = difthi(2,j2,jz) +
     &          fig(1) * znormd * zdfthi(2,j2)
              difthi(3,j2,jz) = difthi(3,j2,jz) +
     &          fig(2) * znormd * zdfthi(3,j2)
              difthi(4,j2,jz) = difthi(4,j2,jz) +
     &          fig(4) * znormd * zdfthi(4,j2)
            enddo
              velthi(1,jz)    = velthi(1,jz) +
     &          fig(3) * znormv * zvlthi(1)
              velthi(2,jz)    = velthi(2,jz) +
     &          fig(1) * znormv * zvlthi(2)
              velthi(3,jz)    = velthi(3,jz) +
     &          fig(2) * znormv * zvlthi(3)
              velthi(4,jz)    = velthi(4,jz) +
     &          fig(4) * znormv * zvlthi(4)
          endif
c
        endif
c
c
c..transfer from diffusivity to convective velocity
c
        if ( lswitch(4) .gt. 0 ) then
c
          if ( thiig(jz) .lt. 0.0 ) then
            velthi(1,jz) = velthi(1,jz) - thiig(jz) * zgth / zrmaj
            thiig(jz) = 0.0
            do j2=1,4
              difthi(1,j2,jz) = 0.0
            enddo
          endif
c
          if ( thdig(jz) .lt. 0.0 ) then
            velthi(2,jz) = velthi(2,jz) - thdig(jz) * zgnh / zrmaj
            thdig(jz) = 0.0
            do j2=1,4
              difthi(2,j2,jz) = 0.0
            enddo
          endif
c
          if ( theig(jz) .lt. 0.0 ) then
            velthi(3,jz) = velthi(3,jz) - theig(jz) * zgte / zrmaj
            theig(jz) = 0.0
            do j2=1,4
              difthi(3,j2,jz) = 0.0
            enddo
          endif
c
          if ( thzig(jz) .lt. 0.0 ) then
            velthi(4,jz) = velthi(4,jz) - thzig(jz) * zgnz / zrmaj
            thzig(jz) = 0.0
            do j2=1,4
              difthi(4,j2,jz) = 0.0
            enddo
          endif
c
        else
c
c..shift from diffusion to convective velocity
c
          if ( abs(cswitch(10)) + abs(cswitch(11)) + abs(cswitch(12))
     &       + abs(cswitch(13)) .gt. zepslon ) then
c
            velthi(1,jz) = velthi(1,jz)
     &       + cswitch(10) * thiig(jz) * zgth / zrmaj
            velthi(2,jz) = velthi(2,jz)
     &       + cswitch(11) * thdig(jz) * zgnh / zrmaj
            velthi(3,jz) = velthi(3,jz)
     &       + cswitch(12) * theig(jz) * zgte / zrmaj
            velthi(4,jz) = velthi(4,jz)
     &       + cswitch(13) * thzig(jz) * zgnz / zrmaj
c
c..alter the effective diffusivities
c  if they are used for more than diagnostic purposes
c
            thiig(jz) = ( 1.0 - cswitch(10) ) * thiig(jz)
            thdig(jz) = ( 1.0 - cswitch(11) ) * thdig(jz)
            theig(jz) = ( 1.0 - cswitch(12) ) * theig(jz)
            thzig(jz) = ( 1.0 - cswitch(13) ) * thzig(jz)
c
            do j2=1,4
              difthi(1,j2,jz) = ( 1.0 - cswitch(10) ) * difthi(1,j2,jz)
              difthi(2,j2,jz) = ( 1.0 - cswitch(11) ) * difthi(2,j2,jz)
              difthi(3,j2,jz) = ( 1.0 - cswitch(12) ) * difthi(3,j2,jz)
              difthi(4,j2,jz) = ( 1.0 - cswitch(13) ) * difthi(4,j2,jz)
            enddo
c
          endif
c
        endif
c
c..end of Weiland model
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c
c..Guzdar-Drake theory (Phys Fluids B 5 (1993) 3712
c..L_p used instead of L_n
c
        zgyrfe = zce * zbtor / zcme  ! electron plasma frequency
        zlare = zvthe / zgyrfe    ! electron Larmor radius
c
c..   Diamagnetic stabilization
c
          zgddia = cswitch(9)
c
c..   Diffusivities
c
        zgdp = abs ( 2. * zpi * ((zq * zlare)**2.) * znuei
     &    * zgpr * 100. * zgddia )

        thdrb(jz) = frb(1) * zgdp * zelonf**cswitch(4)
        therb(jz) = frb(2) * zgdp * zelonf**cswitch(4)
        thirb(jz) = frb(3) * zgdp * zelonf**cswitch(4)
        thzrb(jz) = frb(4) * zgdp * zelonf**cswitch(4)
c
c  add to the fluxes
c
        vflux(1,jz) = vflux(1,jz) + thirb(jz) * zgth / zrmaj
        vflux(2,jz) = vflux(2,jz) + thdrb(jz) * zgnh / zrmaj
        vflux(3,jz) = vflux(3,jz) + therb(jz) * zgte / zrmaj
        vflux(4,jz) = vflux(4,jz) + thzrb(jz) * zgnz / zrmaj
c
        difthi(1,1,jz) = difthi(1,1,jz) + thirb(jz)
        difthi(2,2,jz) = difthi(2,2,jz) + thdrb(jz)
        difthi(3,3,jz) = difthi(3,3,jz) + therb(jz)
        difthi(4,4,jz) = difthi(4,4,jz) + thzrb(jz)
c

c ..................................
c .  the kinetic ballooning model  .
c ..................................
c
c       zbprim and zbc1 computed above under drift model
c
      if (  abs(cswitch(2)) .gt. zepslon
     &   .and.  zgpr .gt. 0.0  ) then
c
      zbprim = abs( zbeta * zgpr / zrmaj )
      zbcoef1 = 1.0
      if ( abs(cswitch(8)) .gt. zepslon ) zbcoef1 = cswitch(8)
      zbc1   = zbcoef1 * abs(zshat)/(1.7*zq**2*zrmaj)
      zelfkb = zelonf**cswitch(5)
c
        zfbthn = exp( min(abs(zlgeps),
     &     max(-abs(zlgeps),cswitch(2)*(zbprim/zbc1 - 1.))) )
c
        zdk = abs( zsound * zrhos**2 * zfbthn * zgpr / zrmaj )
c
        thdkb(jz) = zdk*fkb(1)*zelfkb
        thekb(jz) = zdk*fkb(2)*zelfkb
        thikb(jz) = zdk*fkb(3)*zelfkb
        thzkb(jz) = zdk*fkb(4)*zelfkb
c
c  add to the fluxes
c
        vflux(1,jz) = vflux(1,jz) + thikb(jz) * zgth / zrmaj
        vflux(2,jz) = vflux(2,jz) + thdkb(jz) * zgnh / zrmaj
        vflux(3,jz) = vflux(3,jz) + thekb(jz) * zgte / zrmaj
        vflux(4,jz) = vflux(4,jz) + thzkb(jz) * zgnz / zrmaj
c
        difthi(1,1,jz) = difthi(1,1,jz) + thikb(jz)
        difthi(2,2,jz) = difthi(2,2,jz) + thdkb(jz)
        difthi(3,3,jz) = difthi(3,3,jz) + thekb(jz)
        difthi(4,4,jz) = difthi(4,4,jz) + thzkb(jz)
c
      endif
c
 300  continue
c
c
c   end of the main do-loop over the radial index, "jz"----------
c
      return
      end
