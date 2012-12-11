Module modmmm7_1

!                 >>>>>   N O T E   <<<<<
! {{{ and }}} are VIM editor folding markers. DO NOT REMOVE THEM.

!----- DESCRIPTION -------------------------------------------------{{{
!
! This Fortran 90 module contains a subroutine called mmm7_1 and its
! supporting codes. Together these codes provide the capability for
! computing plasma transport coefficients using the Multi-Mode transport
! model (MMM7.1). The current version of MMM includes contributions from
! four transport models: Weiland (W20), Drift-resistive-inertial
! Ballooning Mode (DRIBM) and Horton ETG (ETG).
!
! This module contains two public subroutines:
!
!    * Subroutine mmm7_1
!    * Subroutine set_mmm_switches
!
! Detailed description can be found in each subroutine.
!
! Revision History
! ----------------
! Mar 1, 2011   Lixiang Luo
!               V7.1 Original Release
!
!-------------------------------------------------------------------}}}

!----- MODULE SPECIFICATIONS ---------------------------------------{{{
Implicit None

Private ! All definitions are assumed private to avoid naming
        ! conflicts with external codes.

!.. Define a better Real*8 type
Integer, Parameter :: R8 = SELECTED_REAL_KIND(12,100)

!.. Physical constants

Real(R8), Parameter :: &
   zpi    = 3.1415926535898_R8  ,&! Pi
   zcc    = 299792458._R8       ,&! Speed of light                 [m/sec]
   zcmu0  = 4E-7_R8 * zpi       ,&! Vacuum magnetic permeability   [henrys/m]
   zceps0 = 1._R8/(zcc**2*zcmu0),&! Vacuum electrical permittivity [farads/m]
   zckb   = 1.602176487E-16_R8  ,&! Energy conversion factor       [Joule/keV]
   zcme   = 9.10938215E-31_R8   ,&! Electron mass                  [kg]
   zcmp   = 1.672621638E-27_R8  ,&! Proton mass                    [kg]
   zce    = 1.602176487E-19_R8    ! Electron charge                [Coulomb]

!.. Machine Epsilon (smallest number so that 1.0+zepslon>1.0)
Real(R8), Parameter :: &
   zepslon = 2._R8**(-53) ,&! As defined by IEEE 754-2008 standard
   zepsqrt = 2._R8**(-26)   ! Square root of Epsilon, approximated

!----- PUBLIC ENTRIES ------------------------------------------
Integer, Parameter, PUBLIC :: &
   MMM_NCH   = 6, &! Maximum number of transport channels
   MMM_NMODE = 3, &! Number of internal models
   MAXNOPT   = 10  ! Maximum number of internal switches

!.. Identifiers for internal models
Integer, Parameter, PUBLIC :: &
   KW20 = 1 ,&! Weiland 20
   KDBM = 2 ,&! DRIBM
   KETG = 3   ! Horton ETG

Integer, Parameter, PUBLIC :: &
! Raised if gelong is not specified and the number of radial points is
! less than 3
   MMM_ERR_ETG_NOT_ENOUGH_ZONE = -1

Public :: mmm7_1
Public :: set_mmm7_1_switches

!!! END OF MODULE SPECIFICATION }}}

Contains

Subroutine mmm7_1 ( &
   !.. Physical inputs
   rmin,  rmaj,   rmaj0, elong,                            &
   ne,    nh,     nz,    nf,                               &
   zeff,  te,     ti,    q,      btor,                     &
   zimp,  aimp,   ahyd,  aimass, wexbs,                    &
   gne,   gni,    gnh,   gnz,    gte,    gti,   gq,        &
   vtor,  vpol,   vpar,  gvtor,  gvpol,  gvpar, gelong,    &
   !.. Physical outputs
   xti,    xdi,    xte,    xdz,    xvt,    xvp,            &
   xtiW20, xdiW20, xteW20, xtiDBM, xdiDBM, xteDBM, xteETG, &
   gammaW20, omegaW20, gammaDBM, omegaDBM, vconv,  vflux,  &
   !.. Feature controls
   npoints, lprint, nprout, nerr, cmodel, cswitch, lswitch )

!----- DESCRIPTION -------------------------------------------------{{{
!  
!  This subroutine calculates coefficients for anomalous transport
!  in tokamak plasmas, using Weiland, ETG and DRIBM models. Detailed
!  description of the subroutine can be found in the documentation.
!
!  Internal switches
!
!    They are passed by LSWITCH and CSWITCH. Their minimum size of the
!    first dimension required is 3 in the current version of MMM. However,
!    this is subjected to change in the future versions. Each entry in the
!    following list is identified by a component name (W20...), a type
!    (L/C), an index and a default value. The corresponding parameter
!    should be located at
!
!       <TYPE LETTER>SWITCH(K<MODEL>,<INDEX>).
!
!    For example, the scaling factor for electromagnetic realm in the ETG
!    component is located in CSWITCH(KETG,2), with a default value of
!    0.06. Detailed instructions on using optional arguments are given in
!    the documentation.
!
!  List
!  -------------------
!  W20   C1   1.0   ExB shearing multiplier (0-OFF, 1-ON)
!        C2   1.0   Momentum pinch scaling factor
!        C3   1E-4  Lower bound of xteW20 [m^2/s]
!        C4   1E+2  Upper bound of xteW20 [m^2/s]
!        C5   1E-4  Lower bound of xtiW20 [m^2/s]
!        C6   1E+2  Upper bound of xtiW20 [m^2/s]
!
!  DBM   C1   0.0   ExB shearing multiplier (0-OFF, 1-ON)
!
!  ETG   L1    1    Jenko threshold selection
!                     0: Disabled, use original Horton model instead
!                     1: Applied to ES but not EM
!                     2: Applied to both ES and EM
!        C1   0.06  scaling factor for electrostatic realm
!        C2   0.06  scaling factor for electromagnetic realm
!
!----- END OF DESCRIPTION ------------------------------------------}}}

Use w20mod

!----- SUBROUTINE SPECIFICATIONS -----------------------------------{{{
Implicit None

! All the following 1-D profiles are assumed to be defined on flux
! surfaces called zone boundaries where the transport fluxes are
! to be computed.  The number of flux surfaces is given by npoints
! (see below).  For example, if you want to compute the transport
! on only one flux surface, set npoints = 1.

Real(R8), Intent(In), Dimension(1:) :: &
   rmin, rmaj, elong, ne, nh, nz, nf, zeff, &
   te, ti, q, btor, zimp, aimp, ahyd, aimass, wexbs, &
   gne, gni, gnh, gnz, gte, gti, gq
!
! rmin(jz)   = half-width of the magnetic surface [m]
! rmaj(jz)   = major radius to geometric center of the magnetic surface [m]
! elong(jz)  = local elongation
! ne(jz)     = electron density [m^-3]
! nh(jz)     = hydrogenic thermal particle density [m^-3]
! nz(jz)     = impurity ion density [m^-3]
! nf(jz)     = density from fast (non-thermal) ions [m^-3]
! zeff(jz)   = mean charge, Z_effective
! te(jz)     = T_e (electron temperature) [keV]
! ti(jz)     = T_i (temperature of thermal ions) [keV]
! q(jz)      = magnetic q-value
! btor(jz)   = ( R B_tor ) / rmaj(jz)  [Tesla]
!
! zimp(jz)   = mean charge of impurities
!            = ( sum_imp n_imp Z_imp ) / ( sum_imp n_imp ) where
!              sum_imp = sum over impurity ions with charge state Z_imp
!
! aimp(jz)   = mean atomic mass of impurities
!            = ( sum_imp n_imp M_imp ) / ( sum_imp n_imp ) where
!              sum_imp = sum over impurity ions, each with mass M_imp
!
! ahyd(jz)   = mean atomic mass of hydrogen ions
!            = ( sum_hyd n_hyd M_hyd ) / ( sum_hyd n_hyd ) where
!              sum_hyd = sum over hydrogenic ions, each with mass M_hyd
!
! aimass(jz) = mean atomic mass of thermal ions
!            = ( sum_i n_i M_i ) / ( sum_i n_i ) where
!              sum_i = sum over all ions, each with mass M_i
!
! wexbs(jz)  = ExB shearing rate [rad/s].
!              See K.H. Burrell, Phys. of Plasmas, 4, 1499 (1997).
!
!   All of the following normalized gradients are at radial points.
!   r = half-width, R = major radius to center of flux surface
!
! gne(jz) = -R ( d n_e / d r ) / n_e
! gni(jz) = -R ( d n_i / d r ) / n_i
! gnh(jz) = -R ( d n_h / d r ) / n_h
! gnz(jz) = -R ( d Z n_Z / d r ) / ( Z n_Z )
! gte(jz) = -R ( d T_e / d r ) / T_e
! gti(jz) = -R ( d T_i / d r ) / T_i
! gq(jz)  =  R ( d q   / d r ) / q
!
! where:
!   n_i  = thermal ion density (sum over hydrogenic and impurity)
!   n_h  = thermal hydrogenic density (sum over hydrogenic species)
!   n_Z  = thermal impurity density,  Z = average impurity charge
!                     summed over all impurities

Real(R8), Intent(In), Optional :: &
   rmaj0  ! Major radius at plasma axis [m]

!.. Profiles related to momentum transport in 2006 Weiland model
Real(R8), Intent(In), Optional, Dimension(1:) :: &
   gvtor, vtor, gvpol, vpol, gvpar, vpar
!
! gvtor(jz)  = Normalized toroidal velocity gradient (r/v_tor)*(dv_tor/dr)
! vtor(jz) = Toroidal velocity [m/s]
! gvpol(jz)  = Normalized poloidal velocity gradient (r/v_pol)*(dv_pol/dr)
! vpol(jz) = Poloidal velocity [m/s]
! gvpar(jz)= Normalized parallel velocity gradient (r/v_par)*(dv_par/dr)
! vpar(jz) = Parallel velocity [m/s]

Real(R8), Intent(In), Optional, Dimension(1:) :: &
   gelong ! Elongation gradient
!
!  gelong = d elong / d rho
!  where elong is local elongation and rho=r/R
!
! > If not associated:
!   The elongation gradient will be calculated internally using other profiles.
! > If associated:
!   The elongation gradient will be assigned by gelong.

! Some output arguments are optional. The calling program is not required to
! use them. When they are associated, the actual arguments have to be
! allocated with enough space to store all return values (npoints is the
! minimal dimension) in advance.

Real(R8), Intent(Out), Dimension(1:) :: & ![m^2/s]
   xti, xdi, xte, xdz, xvt, xvp
!
! Diffusivity profiles, which are the sum of all of the component
! diffusivities in the mmm7.1 model:
!
!  xti(jz) = Effective ion thermal diffusivity
!  xdi(jz) = Effective hydrogenic ion diffusivity
!  xte(jz) = Effective electron thermal diffusivity
!
!  xdz(jz) = Impurity ion diffusivity from the Weiland model
!  xvt(jz) = Toroidal momentum diffusivity from the Weiland model
!  xvp(jz) = Poloidal momentum diffusivity from the Weiland model

! The following component output arrays give the separate contribution from
! each internal model. Note that the momentum diffusivities are only provided
! by the Weiland model. Generally, these arrays are used for diagnostic
! output only.

Real(R8), Intent(Out), Optional, Dimension(1:) :: & ![m^2/s]
   xtiW20, xdiW20, xteW20, xtiDBM, xdiDBM, xteDBM, xteETG
!
! The component output profiles are optional. They should be used for the
! diagnostic output purpose.
!
! xtiW20(jz) = Ion thermal diffusivity from the Weiland (W20) component
! xdiW20(jz) = Hydrogenic ion particle diffusivity from the Weiland component
! xteW20(jz) = Electron thermal diffusivity from the Weiland component
! xtiDBM(jz) = Ion thermal diffusivity from the DRIBM component
! xdiDBM(jz) = Hydrogenic ion diffusivity from the DRIBM component
! xteDBM(jz) = Electron thermal diffusivity from the DRIBM component
! xteETG(jz) = Electron thermal diffusivity from the Horton ETG component

! The following are growth rates and mode frequencies from the
! Weiland model for drift modes such as ITG and TEM.
! These arrays are intended for diagnostic output.
!
!  gamma(jz) = growth rate for the dominating mode at point jz ( 1/sec )
!  omega(jz) = frequency for dominating mode jm at point jz ( radians/sec )
!
Real(R8), Intent(Out), Optional, Dimension(1:,1:) :: &
   omegaW20, gammaW20
Real(R8), Intent(Out), Optional, Dimension(1:) :: &
   omegaDBM, gammaDBM

Real(R8), Intent(Out), Dimension(1:,1:) :: &
   vconv, &! Convective velocities [m/s]
   vflux    ! Flux matrix

Integer, Intent(In) :: &
   npoints ! Number of values in all of the 1-D arrays

Integer, Intent(In) :: &
   lprint, &! Verbose level
   nprout   ! Output unit number for long printout

Integer, Intent(Out) :: &
   nerr ! Error return value

Real(R8), Intent(In), Optional, Dimension(1:) :: &
   cmodel ! Weights of internal models
!
! This is an optional argument. Only the first four elements are used.
! > If associated:
!   cmodel(1)~cmodel(3) are assigned as the weights for Weiland20, DRIBM
!   and ETG models, respectively.
! > If not associated:
!   Equivalent to cmodel=(/1.0,1.0,1.0/)

Real(R8), Intent(In), Optional, Dimension(1:,1:) :: &
   cswitch ! Real internal parameters
!
! This is an optional argument. The second dimension is corresponding to
! the index of an internal model, and the first dimension to the index
! of a real internal parameter for that particular model. For example,
! the first real internal parameter for the Weiland model would be
! stored in cswitch(1,KW20) or cswitch(1,1).
!
! > If associated:
!   The internal parameters are assigned to the given values of the
!   actual argument.
! > If not associated:
!   all actual values default to the internally set values.
!
! An up-to-date list of parameters can be found at the header of
! this subroutine, along with the default values.

Integer, Intent(In), Optional, Dimension(1:,1:) :: &
   lswitch ! Integral internal parameters
!
! This is an optional argument. The second dimension is corresponding to
! the index of an internal model, and the first dimension to the index
! of a integral internal parameter for that particular model. For example,
! the second integral internal parameter for the DRIBM model would be
! stored in cswitch(2,KDBM) or cswitch(2,2).
!
! > If associated:
!   The internal switches are assigned to the given values of the
!   actual argument.
! > If not associated:
!   all actual values default to the internally set values.
!
! An up-to-date list of parameters can be found at the header of
! this subroutine, along with the default values.

!----- LOCAL VARIABLES  -----------------------------------------------

Real(R8) :: &
   cscal(MMM_NMODE)! Weights for individual models

Real(R8) :: &
   csw(MAXNOPT,MMM_NMODE) ! Real internal parameters

Integer :: &
   lsw(MAXNOPT,MMM_NMODE) ! Integral internal parameters

Integer :: &
   jz ! Loop counter

Real(R8) :: &
   amin,       &! Minor radius of plasma
   zep,        &! Inverse aspect ratio
   zgyrfe,     &! Electron gyrofrequency
   zgyrfi,     &! Ion gyrofrequency
   zbeta,      &! Thermal/magnetic energy ratio
   zvthe,      &! Thermal velocity of electrons
   zvthi,      &! Thermal velocity of thermal ions
   zsound,     &! Speed of sound
   zsound_axis,&! Speed of sound at the magnetic center
   zlog,       &! Coulomb logarithm
   zcf,        &! A factor in collision frequency
   zlare,      &! Electron Larmor radius
   zlari,      &! Ion Larmor radius
   zlarpo,     &! Poloidal ion Larmor radius
   zrhos,      &! Ion larmor radius with Te
   zwn,        &! 0.3/rhos
   zgmax,      &! Upper bound a gradients
   znuei,      &! Electron collision frequency
   znude,      &! Electron magnetic drift frequency
                ! Also used as the normalization factor in W20 and DBM
   znuhat       ! Normalized electron collision frequency

Real(R8), Dimension(1:4) :: &
   gamma, omega

!.. Weiland variables{{{
Real(R8), Dimension(MMM_NCH) :: &
   zxteff, &! Diffusivitiy return values
   zvconv    ! Effective convective velocity return values

Real(R8) :: &
   zthte,  &! Ti/Te
   zbetae, &! Thermal/magnetic energy ratio for electrons
   ztzte,  &! Tz/Te
   zfnzne, &! n_z/n_e
   zfnsne, &! n_f/n_e, Nf is the density of superthermal ions
   zftrap, &! Fraction of trapped particles
   zkyrho, &! k_y*rho
   znormd, &! Factor to convert normalized diffusivities
   znormv, &! Factor to convert normalized convective velocities
   shear    ! Magnetic shear

Real(R8) :: &
   wexb_w20 ! ExB shearing rate for W20

Real(R8) :: & ! Copies of gne, gnh, gte and gth
   zgne, zgnh, zgte, zgth !}}}

!.. Horton ETG variables {{{
Real(R8) :: &
   cees,  &! Scaling factor for electrostatic case
   ceem,  &! Scaling factor for electromagnetic case
   cl,    &! 3/2*sqrt(pi/2)
   zwpe,  &! Electron plasma frequency
   zlce,  &! Mixing length of ETG mode
   zdeltae ! Collisionless skin depth

Real(R8) :: &
   zgtec,           &! Critical electron temperature gradient
   zgte_thr_horton, &! in Horton's original model
   zgte_thr_jenko    ! in Jenko's modified model

Real(R8), Dimension(npoints) :: &
   zepa,   &! rmin/rmaj
   dk_deps,&! rate of change of elongation wrt epsilon
   theeg,  &! Electron thermal diffusivity
   thees,  &! Electron thermal diffusivity for electrostatic case
   theem    ! Electron thermal diffusivity for electromagnetic case
!}}}

!.. DRIBM variables {{{
Integer :: nmodes ! Number of DRIBM unstable modes

Real(R8) :: &
   wexb_drbm,&! ExB shearing rate for DRIBM
   difftemp   ! Temporary storage
!}}}

!----- END OF SUBROUTINE SPECIFICATIONS ----------------------------}}}


!----- EXECUTION BODY -------------------------------------------------

!.. Initialize variables {{{

xti = 0._R8; xdi = 0._R8; xte = 0._R8
xdz = 0._R8; xvt = 0._R8; xvp = 0._R8
theeg = 0._R8; thees = 0._R8; theem = 0._R8
vconv = 0._R8; vflux = 0._R8; zxteff= 0._R8
zvconv= 0._R8
If ( Present(gammaW20) ) gammaW20 = 0._R8
If ( Present(omegaW20) ) omegaW20 = 0._R8
If ( Present(gammaDBM) ) gammaDBM = 0._R8
If ( Present(omegaDBM) ) omegaDBM = 0._R8

nerr = 0

!.. Assign internal parameters accoding to optional arguments or
!   default values if the optional argument is omitted
If ( Present( cmodel ) ) Then
   cscal(1:MMM_NMODE) = cmodel(1:MMM_NMODE)
Else
   cscal(1:MMM_NMODE) = 1._R8
End If
!
csw = 0._R8
If ( Present( cswitch ) ) Then
   jz = min( MAXNOPT, size(cswitch,1) )
   csw(1:jz,1:MMM_NMODE) = cswitch(1:jz,1:MMM_NMODE)
Else
   call set_mmm7_1_switches(cmmm=csw)
End If
!
lsw = 0
If ( Present( lswitch ) ) Then
   jz = min( MAXNOPT, size(lswitch,1) )
   lsw(1:jz,1:MMM_NMODE) = lswitch(1:jz,1:MMM_NMODE)
Else
   call set_mmm7_1_switches(lmmm=lsw)
End If

!.. Diffusivities components
If ( Present( xtiW20 ) ) xtiW20 = 0._R8
If ( Present( xdiW20 ) ) xdiW20 = 0._R8
If ( Present( xteW20 ) ) xteW20 = 0._R8
If ( Present( xtiDBM ) ) xtiDBM = 0._R8
If ( Present( xdiDBM ) ) xdiDBM = 0._R8
If ( Present( xteDBM ) ) xteDBM = 0._R8
If ( Present( xteETG ) ) xteETG = 0._R8

!.. Calculate some ETG variables only if ETG is included
dk_deps = 0._R8
If ( Present( gelong ) ) Then
   dk_deps(1:npoints) = gelong(1:npoints)
Else
   If ( cscal(KETG) > 0._R8 ) Then
      If ( npoints >= 3 ) Then
         zepa(1:npoints) = rmin(1:npoints) / rmaj0
         do jz = 2, npoints-1
            dk_deps(jz) = ( elong(jz+1) - elong(jz-1) ) / &
               ( 1E-6_R8 + zepa(jz+1) - zepa(jz-1) )
         end do
         dk_deps(1)=dk_deps(2)
         dk_deps(npoints)=dk_deps(npoints-1)
      Else
         nerr = MMM_ERR_ETG_NOT_ENOUGH_ZONE
         return
      End If
   End If
End If

!}}}

!.. Start the main do-loop over the radial index "jz"
MAINLOOP: &
Do jz = 1, npoints

   ! Calculation of physical quantities {{{

   shear = gq(jz) * rmin(jz) / rmaj(jz)

   !  compute inverse aspect ratio
   zep    = max( rmin(jz)/rmaj(jz), zepslon )

   zgyrfe = zce * btor(jz) / zcme  ! electron plasma frequency
   zgyrfi = zce * btor(jz) / (zcmp * aimass(jz))
   zbeta  = (2._R8* zcmu0 * zckb / btor(jz)**2) * (ne(jz) * te(jz) + nh(jz) * ti(jz))
   zvthe  = sqrt(2._R8* zckb * te(jz) / zcme)
   zvthi  = sqrt(2._R8* zckb * ti(jz) / (zcmp * aimass(jz)))
   zsound = sqrt(zckb * te(jz) / (zcmp * aimass(jz)))

   ! The values of zsound_axis and amin does not affect results.
   ! It is kept only for compatibility with original W20 module.
   zsound_axis = 1._R8
   amin = 1._R8

   zrhos  = zsound / zgyrfi

   zlog   = 37.8_R8-log(sqrt(ne(jz)) / te(jz))
   zcf    = (4._R8* sqrt(zpi) / 3._R8) * &
             (zce / (4._R8* zpi * zceps0))**2 * &
             (zce / zckb) * sqrt( (zce/zcme) * (zce/zckb) )
   znuei  = zcf*sqrt(2._R8)* ne(jz) * zlog * zeff(jz) / (te(jz) * sqrt(te(jz)))

   zlari  = zvthi / zgyrfi
   zlare  = zvthe / zgyrfe
   zlarpo = max(zlari * q(jz) / zep, zepslon)

   zwn    = 0.3_R8 / zrhos
   znude  = 2._R8 * zwn * zrhos * zsound / rmaj(jz)

   zgmax = rmaj(jz) / zlarpo

   !.. Hydrogen species
   zthte  = ti(jz) / te(jz)

   zbetae = 2._R8 * zcmu0 * zckb * ne(jz) * te(jz) / btor(jz)**2

   !.. Impurity species (use only impurity species 1 for now)
   !  assume T_Z = T_H throughout the plasma here

   !Ratio of impurity temperature to electron temperature
   ztzte  = ti(jz) / te(jz)

   !Ratio of impurity density to electron density
   zfnzne = nz(jz) / ne(jz)

   !Fraction of superthermal ions (beam, RF minority) to electrons
   zfnsne = nf(jz) / ne(jz)

   !Fraction of trapped particles
   zftrap = sqrt( 2._R8* rmin(jz) / ( rmaj(jz) * ( 1._R8+ rmin(jz) / rmaj(jz) )))

   !Default k_y rho for Weiland model
   zkyrho = 0.316_R8

   ! End calculation of physical quantities!}}}

   ! Weiland20 {{{
   If ( cscal(KW20) > 0._R8 ) Then

      If ( rmin(jz)<1E-4_R8*rmaj(jz) ) Then
      ! Avoid computation of fluxes at rmin(jz) < 1e-4 * rmaj(jz)
      ! Too close to the plasma center. Assume minimal transportation.
         xte(jz)=csw(3,KW20)
      Else

         wexb_w20 = csw(1,KW20)*wexbs(jz) ! ExB effect included

         Call w20main ( lprint, nprout, &
              te(jz), ne(jz), vtor(jz), vpol(jz), vpar(jz), &
              btor(jz), amin, rmin(jz), rmaj(jz), amin/rmaj0, &
              aimp(jz), ahyd(jz), zimp(jz), &
              gte(jz), gti(jz), gti(jz), &
              gne(jz), gni(jz), gnz(jz), &
              gvtor(jz), gvpol(jz), gvpar(jz), &
              zkyrho, zthte, ztzte, zfnzne, zfnsne, zftrap, &
              q(jz), shear, elong(jz), wexb_w20, &
              zsound_axis, &
              zxteff, zvconv, &
              gamma(1:4), omega(1:4)  )

         !.. If nerr not equal to 0 an error has occured
         If (nerr /= 0) Return

         If (Present( gammaW20 )) gammaW20(1:4,jz) = gamma(1:4) * znude
         If (Present( omegaW20 )) omegaW20(1:4,jz) = omega(1:4) * znude

         !  Conversion factors for diffusion and convective velocity
         znormd = 2._R8 * zsound * zrhos**2 / ( rmaj(jz) * zkyrho )
         znormv = 2._R8 * zsound * zrhos**2 / ( rmaj(jz)**2 * zkyrho )

         !  Store effective diffusivities from Weiland19 model
         xti(jz) = max(csw(5,KW20),min(csw(6,KW20), cscal(KW20)*zxteff(1) ))
         xdi(jz) = cscal(KW20) * zxteff(2)
         xte(jz) = max(csw(3,KW20),min(csw(4,KW20), cscal(KW20)*zxteff(3) ))
         xdz(jz) = cscal(KW20) * zxteff(4)
         xvt(jz) = cscal(KW20) * zxteff(5)
         xvp(jz)= cscal(KW20) * zxteff(6)

         !.. Weiland's own diffusivities
         If ( Present( xtiW20 ) ) xtiW20(jz) = xti(jz)
         If ( Present( xdiW20 ) ) xdiW20(jz) = xdi(jz)
         If ( Present( xteW20 ) ) xteW20(jz) = xte(jz)

         !  start computing the fluxes
         vflux(1,jz) = cscal(KW20) * xti(jz) * gti(jz) / rmaj(jz)
         vflux(2,jz) = cscal(KW20) * xdi(jz) * gnh(jz) / rmaj(jz)
         vflux(3,jz) = cscal(KW20) * xte(jz) * gte(jz) / rmaj(jz)
         vflux(4,jz) = cscal(KW20) * xdz(jz) * gnz(jz) / rmaj(jz)

         !  Add momentum transport convective velocities from Reynold's source term
         !  in Weiland model and other pinches
         vconv(1,jz) = cscal(KW20)*zvconv(1)
         vconv(2,jz) = cscal(KW20)*zvconv(2)
         vconv(3,jz) = cscal(KW20)*zvconv(3)
         vconv(4,jz) = cscal(KW20)*csw(2,KW20)*zvconv(4)
         vconv(5,jz) = cscal(KW20)*csw(2,KW20)*zvconv(5)
         vconv(6,jz) = cscal(KW20)*csw(2,KW20)*zvconv(6)

         vconv(5,jz) = vconv(4,jz) + vconv(5,jz)

      End If
   End If ! End of Weiland model!}}}

   ! DRIBM {{{
   If ( cscal(KDBM) > 0._R8 ) Then

      !.. An upper bound for gradients
      zgmax = rmaj(jz) / zlarpo
      zgne = sign( min( abs( gne(jz) ), zgmax ), gne(jz) )
      zgnh = sign( min( abs( gnh(jz) ), zgmax ), gnh(jz) )
      zgte = sign( min( abs( gte(jz) ), zgmax ), gte(jz) )
      zgth = sign( min( abs( gti(jz) ), zgmax ), gti(jz) )

      zkyrho = 0.2_R8
      zwn    = 0.2_R8 / zrhos
      znude  = 2._R8 * zwn * zrhos * zsound / rmaj(jz)
      znuhat = znuei / znude

      wexb_drbm = csw(1,KDBM) * wexbs(jz)/znude

      Call drbm ( zgne,   zgnh,   zgte,   zgth, &
         zthte,   zbetae, znuhat, q(jz),  zkyrho, wexb_drbm, &
         zxteff, nmodes, gamma(1), omega(1), nerr )

      If (Present( gammaDBM )) gammaDBM(jz) = gamma(1) * znude
      If (Present( omegaDBM )) omegaDBM(jz) = omega(1) * znude

      !.. Conversion factors for diffusion and convective velocity
      znormd = 2._R8 * zsound * zrhos**2 / ( rmaj(jz) * zkyrho )
      znormv = 2._R8 * zsound * zrhos**2 / ( rmaj(jz)**2 * zkyrho )

      difftemp = cscal(KDBM)*znormd* max(0._R8,zxteff(1))
      difftemp = min(25._R8,difftemp)
      If ( Present( xtiDBM ) ) xtiDBM(jz) = difftemp
      xti(jz) = xti(jz) + difftemp

      difftemp = cscal(KDBM)*znormd* max(0._R8,zxteff(2))
      difftemp = min(25._R8,difftemp)
      If ( Present( xdiDBM ) ) xdiDBM(jz) = difftemp
      xdi(jz) = xdi(jz) + difftemp

      difftemp = cscal(KDBM)*znormd* max(0._R8,zxteff(3))
      difftemp = min(25._R8,difftemp)
      If ( Present( xteDBM ) ) xteDBM(jz) = difftemp
      xte(jz) = xte(jz) + difftemp

      !.. Start computing the fluxes
      vflux(1,jz) = vflux(1,jz) + xti(jz) * zgth / rmaj(jz)
      vflux(2,jz) = vflux(2,jz) + xdi(jz) * zgnh / rmaj(jz)
      vflux(3,jz) = vflux(3,jz) + xte(jz) * zgte / rmaj(jz)

   End If!}}}

   ! Horton ETG{{{
   If ( cscal(KETG) > 0._R8 ) Then

      cees  = csw(1,KETG)
      ceem  = csw(2,KETG)

      cl = 1.88_R8

      !Thermal velocity of electrons [ m / s ]
      zvthe = sqrt( zckb * te(jz) / zcme )

      !Electron gyrofrequency [ 1 / s ]
      zgyrfe = zce * btor(jz) / zcme

      !Electron gyroradius [ m ]
      zlare = zvthe / zgyrfe

      !Electron plasma frequency [ 1 / s ]
      zwpe = sqrt( (ne(jz)*zce**2) / (zcme * zceps0) )

      !Mixing length of ETG mode
      zlce = max( 0._R8, q(jz) * zlare * gte(jz))

      !Collisionless skin depth
      zdeltae = zcc / zwpe

      !Critical electron temperature gradient in Horton's paper
      zgte_thr_horton = cl * abs( shear ) / q(jz) * &
           (1._R8 + zeff(jz) * te(jz) / ti(jz) )

      !Critical electron temperature gradient by Jenko ETG model
      zgte_thr_jenko = max( ( 1._R8 + zeff(jz) * te(jz) / ti(jz) ) * &
           ( 1.33_R8 + 1.91_R8*abs( shear ) / q(jz) ) * &
           ( 1._R8 - 1.5_R8 * zep ) * ( 1._R8 + 0.3_R8 * zep * dk_deps(jz) ), &
           0.8_R8 * gne(jz) )

      If ( zlce > zdeltae ) then ! Electrostatic regime

         If ( lsw(1,KETG) == 1 .or. lsw(1,KETG) == 2  ) Then
            zgtec = zgte_thr_jenko
         Else
            zgtec = zgte_thr_horton
         End If

         theeg(jz) = cees * (q(jz) * zlare)**2 * zvthe &
            * gte(jz)* sqrt( abs(gte(jz)) ) / rmaj(jz) * &
            max( 0._R8, ( gte(jz) - zgtec ) )

      Else if ( zlce <= zdeltae ) Then ! electro-magnetic regime

         theeg(jz) = ceem * zdeltae**2 * zvthe &
                 * sqrt(max( 1E-6_R8, gte(jz) )) / rmaj(jz)

         If ( lsw(1,KETG) == 2 ) Then ! Apply Jenko threshold to EM regime
            If ( gte(jz) < zgte_thr_jenko ) Then
               theeg(jz) = 0._R8
            Else
               theeg(jz) = tanh( gte(jz) - zgte_thr_jenko ) * theeg(jz)
            End If
         End If

      End If

      !Apply bound on diffusivity
      theeg(jz) = min( theeg(jz), 1E2_R8 )

      !.. ETG's own diffusivity
      difftemp = cscal(KETG) * theeg(jz)

      If ( Present( xteETG ) ) xteETG(jz) = difftemp

      xte(jz) = xte(jz) + difftemp

      !thr = zgte_thr_jenko
   End If!}}}

End Do MAINLOOP ! End of the main do-loop over the radial index, "jz"

End Subroutine mmm7_1


Subroutine set_mmm7_1_switches & !{{{
   ( cmmm, lmmm,                                         &
     KW20_C_EXB,     KW20_C_MOM_PINCH_SCALE,             &
     KW20_C_XTE_MIN, KW20_C_XTE_MAX,                     &
     KW20_C_XTI_MIN, KW20_C_XTI_MAX,                     &
     KDBM_C_EXB,                                         &
     KETG_C_CEES_SCALE, KETG_C_CEEM_SCALE, KETG_L_NLTHR )

!----- DESCRIPTION -------------------------------------------------{{{
!  
!  This subroutine set up the argument arrays (lswitch and cswitch)
!  for subroutine mmm7_1. Users can use a <keyword>=<value> approach
!  to set up internal parameters. Detailed discussion on this
!  subroutine can be found in the documentation.
!
!  *** Input ***
!     <keyword>=<value> pairs
!  *** Output ***
!     cmmm and lmmm, which are passed to cswitch and lswtich for mmm7_1
!
!----- END OF DESCRIPTION ------------------------------------------}}}

Implicit None

Real(R8), Optional, Intent(Out), Dimension(:,:) :: cmmm

Integer, Optional, Intent(Out), Dimension(:,:) :: lmmm

Real(R8), Optional, Intent(In) :: &
   KW20_C_EXB,             &! W20 C1
   KW20_C_MOM_PINCH_SCALE, &! W20 C2
   KW20_C_XTE_MIN,         &! W20 C3
   KW20_C_XTE_MAX,         &! W20 C4
   KW20_C_XTI_MIN,         &! W20 C5
   KW20_C_XTI_MAX,         &! W20 C6
   KDBM_C_EXB,             &! DBM C1
   KETG_C_CEES_SCALE,      &! ETG C1
   KETG_C_CEEM_SCALE        ! ETG C1

Integer, Optional, Intent(In) :: &
   KETG_L_NLTHR             ! ETG L1

!------- Execution Body --------------------------------------------------

If ( Present( cmmm ) ) Then
   cmmm = 0._R8

   cmmm(1,KW20) = 1.0_R8
   If (Present(KW20_C_EXB)) cmmm(6,KW20) = KW20_C_EXB

   cmmm(2,KW20) = 1._R8
   If (Present(KW20_C_MOM_PINCH_SCALE)) cmmm(1,KW20) = KW20_C_MOM_PINCH_SCALE

   cmmm(3,KW20) = 1E-4_R8
   If (Present(KW20_C_XTE_MIN)) cmmm(2,KW20) = KW20_C_XTE_MIN

   cmmm(4,KW20) = 1E+2_R8
   If (Present(KW20_C_XTE_MAX)) cmmm(3,KW20) = KW20_C_XTE_MAX

   cmmm(5,KW20) = 1E-4_R8
   If (Present(KW20_C_XTI_MIN)) cmmm(4,KW20) = KW20_C_XTI_MIN

   cmmm(6,KW20) = 1E+2_R8
   If (Present(KW20_C_XTI_MAX)) cmmm(5,KW20) = KW20_C_XTI_MAX

   cmmm(1,KDBM) = 0.0_R8
   If (Present(KDBM_C_EXB)) cmmm(1,KDBM) = KDBM_C_EXB

   cmmm(1,KETG) = 6E-2_R8
   If (Present(KETG_C_CEES_SCALE)) cmmm(1,KETG) = KETG_C_CEES_SCALE

   cmmm(2,KETG) = 6E-2_R8
   If (Present(KETG_C_CEEM_SCALE)) cmmm(2,KETG) = KETG_C_CEEM_SCALE
End If

If ( Present( lmmm ) ) Then
   lmmm = 0

   lmmm(1,KETG) = 1
   If (Present(KETG_L_NLTHR)) lmmm(1,KETG) = KETG_L_NLTHR
End If

End Subroutine set_mmm7_1_switches!}}}


Subroutine drbm & !{{{
   (  gnein,   gnhin,   gtein,   gthin,               &
      tauhin,  betaein, vefin,   q,  ekyrhoin,  wexb, &
      xteff,  nmodes,  gamma_max,  omega_max,  nerr  )

!----- DESCRIPTION -------------------------------------------------{{{
!  
!  This private module subroutine calculates DRIBM transport
!  coefficients for a single radial point. It is not user-callable.
!
!----- END OF DESCRIPTION ------------------------------------------}}}

! SUBROUTINE SPECIFICATIONS {{{
Use w20mod
Implicit None

!----- CONSTANTS ---------------------------------------------------

Integer, Parameter :: &
   neq = 6 ! Number of equations

!----- ARGUMENTS ---------------------------------------------------

Real(R8), Intent(In) :: &
   gnein,  gnhin,  gtein, gthin, &
   tauhin, betaein, q, vefin, ekyrhoin, wexb

Integer, Intent(Out) :: &
   nmodes, & ! nmodes = number of unstable modes
   nerr      ! Error value

Real(R8), Intent(Out) :: &
   gamma_max, omega_max

Real(R8), Intent(Out) :: &
   xteff(:)

!----- LOCAL VARIABLES  --------------------------------------------

Integer :: j1, j

Real(R8) :: &
   vef, zflh, zflxph, zflxnh, &
   zflxpe, zphsph, zphsnh, zphspe, ztemp1, &
   zreal, zimag

!.. Arrays for eigenvalue problem
! ( zamr(i,j), zami(i,j) ) = matrix A
! ( zbmr(i,j), zbmi(i,j) ) = matrix B
!
! Note that the eigenvalues are
!
!   zomega(j) = zalfr(j) / zbeta(j)
!   zgamma(j) = zalfi(j) / zbeta(j)
!
Real(R8), Dimension(1:neq) :: zgamma, zomega
!
! and the eigenvectors are
!
!   zevec(j) = cmplx( zvr(j), zvi(j) )
! Here, zbeta(j) will be 0.0 in the case of an infinite eigenvalue
!
! Here, zbeta(j) will be 0.0 in the case of an infinite eigenvalue
Real(R8), Dimension(neq,neq) :: &
   zamr, zami, zbmr, zbmi, &
   zvr, zvi

Real(R8), Dimension(neq) :: &
   zalfr, zalfi, zbeta

Complex*16 :: &
   zevec(neq,neq) = 0._R8

!  These are the thermal and particle fluxes and effective diffusivities
!  Normalized transport fluxes                  Eff. Diffusivities
!       zflxm(1)   Ion thermal flux               zchim(1)    chi_i
!       zflxm(2)   Hydrogenic flux                zchim(2)    D_i
!       zflxm(3)   Electron thermal flux          zchim(3)    chi_e
Real(R8) :: zflxm(3), zchim(3)

Real(R8) :: &
   zkpsh, zalp, zalf

Integer :: ifail

! END SUBROUTINE SPECIFICATIONS }}}

!.. Initialize variables  {{{
!
! CAUTION: Initial values during definition may not have the expected
! effects, because all local variables have an implicit "SAVE"
! attribute in Fortran 90.
!
zomega = 0._R8
zgamma = 0._R8
zalfr  = 0._R8
zalfi  = 0._R8
zbeta  = 0._R8
zflxm  = 0._R8
zchim  = 0._R8
zami   = 0._R8
zbmr   = 0._R8
zevec  = ( 0._R8, 0._R8 )
zamr   = 0._R8
zami   = 0._R8
zbmr   = 0._R8
zbmi   = 0._R8
zvr    = 0._R8
zvi    = 0._R8
xteff  = 0._R8
nmodes = 0
nerr   = 0

vef = 1._R8 * vefin ! Possible scaling

!.. End of initialization

!.. compute the rest of the dimensionless variables needed
!     zflh  = k_y^2 \rho_{sH}^2

zflh = ekyrhoin**2

If ( gnhin + gthin < zepsqrt ) then
   zalp  = 0._R8
   zalf  = 0._R8
   zkpsh = 0._R8
Else
   zalp = 1._R8
   zkpsh = 1._R8 * 0.5_R8 * sqrt( zalp / zflh ) / q ![eric]scaling
End If
!}}}

!.. Equations for e \phi / T_e, T_H, n_H, T_e, A// and v// {{{

!.. Hydrogen density
zamr(1,1) = - 1._R8 + 0.5_R8 * (gnhin - zflh * tauhin * (gnhin+gthin))
zamr(1,2) = - tauhin
zamr(1,3) = - tauhin
zamr(1,6) = zkpsh
zbmr(1,1) = zflh
zbmr(1,3) = 1._R8

!.. Total momentum equation
zamr(2,2) = zkpsh*tauhin
zamr(2,3) = zkpsh*(1._R8+tauhin)
zamr(2,4) = zkpsh
zamr(2,5) = -0.5_R8 * (tauhin*(gthin+gnhin) + gtein+gnein)
zbmr(2,6) = 1._R8

!.. Hydrogen energy
zamr(3,1) = -0.5_R8 * ( gthin - (2._R8/3._R8)*gnhin )
zamr(3,2) = tauhin * (5._R8/3._R8)
zbmr(3,2) = - 1._R8
zbmr(3,3) =  (2._R8/3._R8)

!.. Electron energy
zamr(4,1) =  0.5_R8 * ( gtein - (2._R8/3._R8)*gnein )
zamr(4,4) =  (5._R8/3._R8)
zamr(4,5) = - zkpsh*0.96D0*zflh / betaein
zami(4,4) = - 3600._R8*zkpsh**2 / vef
zami(4,5) =   918._R8*zkpsh*gtein / vef
zbmr(4,3) = - (2._R8/3._R8)
zbmr(4,4) = 1._R8

!.. Vorticity equation
zamr(5,1) = - 0.5_R8*zflh*tauhin*(gnhin+gthin)
zamr(5,2) = - tauhin
zamr(5,3) = - (1._R8 + tauhin)
zamr(5,4) = - 1._R8
zamr(5,5) = 2._R8*zkpsh*zflh / betaein
zbmr(5,1) = zflh

!.. Ohms law
zamr(6,1) = zkpsh
zamr(6,3) = - zkpsh
zamr(6,4) = - 1.71_R8*zkpsh
zamr(6,5) = 0.5_R8*(gnein + 1.71_R8*gtein)
zami(6,5) = - 0.0054_R8 * vef * zflh / betaein
zbmr(6,5) = 1._R8 + 2._R8*0.0054_R8*zflh / betaein
zbmr(6,6) = - 0.0054_R8

!.. Find the eigenvalues and eigenvectors using ACM/TOMS routine 535
Ifail = -1

!.. Solution of the generalized eigenvalue problem
Call r8tomsqz_mmm(neq,neq,zamr,zami,zbmr,zbmi,zalfr,zalfi,zbeta,zvr,zvi,ifail)
!
If ( ifail /= 0 ) Then
   nerr = ifail
   Return
End If

!.. Storing the results -  eigenvalues and eigenvectors
gamma_max = 0._R8
omega_max = 0._R8
Do j=1,neq
   ztemp1 = zbeta(j)
   If ( abs(zbeta(j)) < zepsqrt ) ztemp1 = zepsqrt

   zomega(j) = zalfr(j)  / ztemp1
   zgamma(j) = zalfi(j)  / ztemp1
   
   If ( zgamma(j) > gamma_max ) Then ! To find the most significant mode
      gamma_max = zgamma(j)
      omega_max = zomega(j)
   End If

   Do j1=1,neq
      zevec(j1,j) = DCMPLX ( zvr(j1,j), zvi(j1,j))
   End Do
End Do
!}}}

!.. Calculation of diffusivities and fluxes if transport exists {{{
If_HAVETXP: &
If ( abs(gamma_max) < 1E-8_R8 ) Then ! Yes: No transport, go to finalization

   ! This is to avoid artifacts
   gamma_max = 0._R8
   omega_max = 0._R8

Else

! Real and imaginary parts
   nmodes = 0

! Loop over number of possible modes
   LOOP_ALLEQ: &
   Do j=1,neq

      !.. Compute effective diffusivities directly from eigenvalues
      ! Assuming eigenvectors are arranged in the order of
      !    e\Phi/Te, Ti, nH, Te, A//, v//

      !.. Calculate norm of solution in terms of (e\phi/T_e) **2
      ztemp1 = DBLE(zevec(1,j)) *  DBLE(zevec(1,j)) &
            + DIMAG(zevec(1,j)) * DIMAG(zevec(1,j))

      !.. Check if current mode j is unstable
      If ( zgamma(j)>zepsqrt .and. ztemp1>zepsqrt ) Then

         nmodes = nmodes + 1

         !.. Define fluxes : Thermal hydrogen flux
         zreal =  DBLE(zevec(2,j)) +  DBLE(zevec(3,j))
         zimag = DIMAG(zevec(2,j)) + DIMAG(zevec(3,j))
         !
         !.. phase difference
         zphsph = - ( zimag * DBLE(zevec(1,j) ) &
                  -   zreal *DIMAG(zevec(1,j) ) ) / ztemp1
         !
         !.. flux from j:th mode
         zflxph = 2._R8* zphsph * max(0._R8,zgamma(j)-abs(wexb))**2
         !
         !.. flux summed over all unstable modes
         zflxm(1) = zflxm(1) + zflxph
         !
         !.. Define hydrogen density flux - phase shift
         zphsnh = - ( DIMAG(zevec(3,j)) * DBLE(zevec(1,j)) &
                - DBLE(zevec(3,j)) * DIMAG(zevec(1,j)) ) / ztemp1
         !
         !.. Flux from j:th mode
         zflxnh = 2._R8* zphsnh * max(0._R8,zgamma(j)-abs(wexb))**2
         !
         !.. Flux summed over all unstable modes
         zflxm(2) = zflxm(2) + zflxnh
         !
         !.. Define thermal electron flux
         zreal =  DBLE(zevec(4,j)) +  DBLE(zevec(3,j))
         zimag = DIMAG(zevec(4,j)) + DIMAG(zevec(3,j))
         !
         zphspe = - ( zimag * DBLE(zevec(1,j)) &
                  -  zreal * DIMAG(zevec(1,j)) ) / ztemp1
         !
         zflxpe = 2._R8* zphspe * max(0._R8,zgamma(j)-abs(wexb))**2
         zflxm(3) = zflxm(3) + zflxpe

      End If

   End Do LOOP_ALLEQ !.. End of flux and diffusivity definitions loop

   !..compute effective total diffusivities
   zchim(1) = zflxm(1) / sign( max( abs( gthin ), zepsqrt ),  gthin )
   zchim(2) = zflxm(2) / sign( max( abs( gnhin ), zepsqrt ),  gnhin )
   zchim(3) = zflxm(3) / sign( max( abs( gtein ), zepsqrt ),  gtein )

   !..save effective diffusivities and fluxes
   xteff(1:3) = zchim(1:3)

   ! Effective growth rate
   gamma_max = max( 0._R8, gamma_max - abs(wexb) )

End If IF_HAVETXP! No transport }}}

End Subroutine drbm!}}}


End Module modmmm7_1
