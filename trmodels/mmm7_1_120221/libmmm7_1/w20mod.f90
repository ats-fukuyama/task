module w20mod

!----- DESCRIPTION -------------------------------------------------{{{
!
! This Fortran 90 module contains the Weiland20 anomalous transport
! subroutine w20main and its supporting codes. 
!
! This module contains a single public subroutines:
!
!    * Subroutine w20main
!
! Detailed description can be found in each subroutine.
!
! Revision History
! ----------------
! Apr 7, 2011   Lixiang Luo
!               Improved diagnostic output handling
!
! Mar 1, 2011   Lixiang Luo
!               Repackaging as a Fortran 90 module
!               V7.1 Original Release
!
! 2008          Federico Halpern
!               Original codes of MMM08
!
!-------------------------------------------------------------------}}}


! Module specification {{{
  implicit none

  private

  !Precision constants
  INTEGER, PARAMETER :: RKIND = SELECTED_REAL_KIND(12,100)
  INTEGER, PARAMETER :: CKIND = RKIND

  !Order of the dispersion relation
  integer, parameter :: neq = 9

  !Flags used for searching ion or electron drift modes
  integer :: searchmode
  integer, parameter :: S_ION = 2
  integer, parameter :: S_ELC = 3
  integer, parameter :: S_ANY = 0
  integer, parameter :: N_NONE = 0
  real, parameter :: grdmin = 0.001e0_rkind

  !Printing flag
  integer :: lprint
  integer :: hprint      ! Unit number for diagnostic output

  !Species mass
  real(rkind) :: aimp    ! Atomic mass a_z (a.m.u.)
  real(rkind) :: ahyd    ! Atomic mass a_i (a.m.u.)

  !Species charge
  real(rkind) :: zimp    ! Impurity charge Z_z
  real(rkind) :: zhyd    ! Hydrogenic charge Z_i
  real(rkind) :: zfion   ! Fast ion charge Z_f
  real(rkind) :: zeff    ! Effective charge

  !Density and temperature gradients
  !Defined as -(R/U)(dU/dr), derivative respect to minor radius r (not rho)
  real(rkind) :: gte     ! Electron temperature gradient 
  real(rkind) :: gti     ! Ion temperature gradient
  real(rkind) :: gtz     ! Impurity temperature gradient
  real(rkind) :: gne     ! Electron density gradient
  real(rkind) :: gni     ! Ion density gradient
  real(rkind) :: gnz     ! Impurity density gradient
  real(rkind) :: gvt     ! Toroidal velocity gradient
  real(rkind) :: gvp     ! Poloidal velocity gradient
  real(rkind) :: gvpar   ! Parallel velocity gradient

  !Correlation lengths
  real(rkind) :: kyrho   ! k_y rho = 0.316, ( k_y rho )^2 = 0.1
  real(rkind) :: flh     ! Hydrogenic modes correlation length ( k_y rho )**2~0.10
  real(rkind) :: fle     ! Electron modes correlation length
  real(rkind) :: zflh    ! Atomic mass * Correlation length / hydrogenic charge squared
  real(rkind) :: zflz    ! Same, for impurities
  real(rkind) :: epsncorr ! epsilon_n correlation factor ( \epsilon_n = 2 L_n / L_B = 2 / gne )

  !Electron temperature and temperature ratios
  real(rkind) :: te      ! Te
  real(rkind) :: th      ! Th (hydrogenic ion temperature)
  real(rkind) :: tz      ! Tz (impurity temperature)
  real(rkind) :: tau     ! Te/Ti
  real(rkind) :: tau_inv ! Ti/Te
  real(rkind) :: tauh    ! Ti/Te
  real(rkind) :: tauz    ! Tz/Te

  real(rkind) :: ztauh   ! Z_i Ti/Te
  real(rkind) :: ztauz   ! Z_z Tz/Te

  !Electron density and density ratios
  real(rkind) :: ne      ! Electron density n_e
  real(rkind) :: nh      ! Main hydrogenic ion density n_h
  real(rkind) :: nz      ! Impurity ion density n_z
  real(rkind) :: ns      ! Fast particle density n_s
  real(rkind) :: znz_ne  ! Z_imp nz/ne
  real(rkind) :: zns_ne  ! Z_(fast ion) ns/ne
  real(rkind) :: ne_nh   ! ne/nh
  real(rkind) :: nh_ne   ! ni/ne (fraction of main ions)
  real(rkind) :: nz_ne   ! nz/ne (fraction of impurity ions)
  real(rkind) :: ns_ne   ! ns/ne (fraction of fast ions)
  real(rkind) :: fte     ! n_e (trapped) / n_e

  !Fractions
  real(rkind) :: fft
  real(rkind) :: fzft
  real(rkind) :: ftrt
  real(rkind) :: gm
  real(rkind) :: bta
  real(rkind) :: xt

  !Collision frequency related
  real(rkind) :: vei     ! Electron-ion collision rate
  real(rkind) :: vef     ! nu effective
  real(rkind) :: bt
  real(rkind) :: bt1

  !Geometry and equilibrium quantities
  real(rkind) :: rhos    ! Gyroradius
  real(rkind) :: rmaj    ! Local major radius
  real(rkind) :: rmaj0   ! Major radius at magnetic axis
  real(rkind) :: rmin    ! Local minor radius
  real(rkind) :: amin    ! Plasma minor radius a
  real(rkind) :: eps0    ! Global aspect ratio a/R_0
  real(rkind) :: eps     ! Local aspect ratio r/R
  real(rkind) :: btor    ! Local toroidal magnetic field
  real(rkind) :: q       ! Magnetic safety factor
  real(rkind) :: shear   ! Magnetic shear (r/q)(dq/dr)
  real(rkind) :: shat    ! Effective magnetic shear modified for elongation
  real(rkind) :: kappa   ! Plasma elongation
  real(rkind) :: alpha_MHD   ! MHD alpha

  !Plasma beta
  real(rkind) :: beta    ! Plasma beta
  real(rkind) :: betae   ! Electron beta

  !ExB shear rates
  real(rkind) :: wexb    ! Shear rate w_ExB

  !Saved internal variables used in dispersion relation

  real(rkind) :: str = 7.0e0_rkind/3.0e0_rkind
  real(rkind) :: ftr = 5.0e0_rkind/3.0e0_rkind   
  real(rkind) :: tvr = 2.0e0_rkind/3.0e0_rkind
  real(rkind) :: em2 = 1.0e0_rkind
  real(rkind) :: em1 = 1.0e0_rkind
  real(rkind) :: em  = 1.0e0_rkind
  real(rkind) :: kpc = 1.0e0_rkind

  !Plasma velocity normalized by the local ion soundspeed Cs = sqrt ( Te/mi )
  real(rkind) :: vtor       ! Toroidal velocity
  real(rkind) :: vpol       ! Poloidal velocity
  real(rkind) :: vpar       ! Parallel velocity

  !Some other useful quantities
  real(rkind) :: wde        ! 2.0 * kyrho**2 * csound / R
  real(rkind) :: Csound     ! speed of sound
  real(rkind) :: Csound0    ! speed of sound at magnetic axis
  real(rkind) :: diffBohm   ! Bohm diffusion coefficient

public :: w20main!}}}

contains

Subroutine w20main &!{{{
   ( i_print, h_print, &
     z_te,   z_ne,    z_vtor, z_vpol, z_vpar, &
     z_btor, z_amin,  z_rmin, z_rmaj, z_eps0, &
     z_aimp, z_ahyd,  z_zimp, &
     z_gte,  z_gti,   z_gtz , &
     z_gne,  z_gni,   z_gnz, &
     z_gvt,  z_gvp,   z_gvpar, &
     z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
     z_q,    z_shear, z_kappa, z_wexb, z_Csound0, &
     diffs,  vconv,   gamma,  omega)
!For definition of input parameters look at file w20_README.txt
!or at the variable definitions in w20data.f90

  implicit none

  !Variable declation
  integer, intent(in) :: i_print, h_print
  
  real(rkind), intent(in) :: &
       z_te,   z_ne,    z_vtor,  z_vpol,  z_vpar, &
       z_btor, z_amin,  z_rmin,  z_rmaj,  z_eps0, &
       z_aimp, z_ahyd,  z_zimp,  &
       z_gte,  z_gti,   z_gtz ,  &
       z_gne,  z_gni,   z_gnz,   &
       z_gvt,  z_gvp,   z_gvpar, & 
       z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
       z_q,    z_shear, z_kappa, z_wexb,  z_Csound0

  integer :: ierr, iiu, ieu, iimax, iemax, eiu, eeu, eimax, eemax, j

  !Resulting fluxes
  real(rkind), intent(inout), dimension(6) :: diffs, vconv

  ! Growth rates and frequency are given in this order:
  ! 1 - ii - max ion transport due to the most significant ion mode
  ! 2 - ie - max electron transport due to the most significant ion mode
  ! 3 - ei - max ion transport due to the most significant electron mode
  ! 4 - ee - max electron transport due to the most significant electron mode
  Real(rkind), Intent(Out), Dimension(1:) :: &
     gamma, &! Growth rate
     omega   ! Frequency

  !Individual fluxes from ion and electron call
  real(rkind), dimension(6) :: idiffs, ivconv, ediffs, evconv

  !Eigenvectors resulting from solved system of equations
  real(rkind), dimension(neq,neq) :: izvr, izvi, ezvr, ezvi, zvr, zvi

  !Complex frequency of the mode 
  complex(ckind), dimension(neq) :: izz, ezz, zz

  !Geometric average factor multiplying W_de, used to compute momentum thermoelectric pinch
  real(rkind) :: G_ave_i, G_ave_e

  !Parallel wavenumber of drift mode, used to compute momentum pinches
  real(rkind) :: kps_i, kps_e

  !Logical variable decides if separate correlation length needed
  logical :: l_elc

  !Initialize all internal variables
  call w20init( i_print, h_print, &
     z_te,   z_ne,    z_vtor,  z_vpol,  z_vpar, &
     z_btor, z_amin,  z_rmin,  z_rmaj,  z_eps0, &
     z_aimp, z_ahyd,  z_zimp, &
     z_gte,  z_gti,   z_gtz,  &
     z_gne,  z_gni,   z_gnz,  &
     z_gvt,  z_gvp,   z_gvpar,&
     z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
     z_q,    z_shear, z_kappa, z_wexb,  z_Csound0 )


  ! Algorithm starts here

  !Initialize geometric factor to 1.0 (circular)
  G_ave_e = 1.0
  G_ave_i = 1.0

  !Search for most unstable ion mode
  searchmode = S_ION

  !Variable correlation length for ion modes

  call w20solv( izz, izvr, izvi, iiu, ieu, iimax, iemax, G_ave_i, kps_i, ierr )

  ! Record the most signifiant modes
  gamma = 0D0
  omega = 0D0
  if ( iimax > 0 ) then
     gamma(1) = dimag( izz(iimax) )
     omega(1) = dreal( izz(iimax) )
  end if
  if ( iemax > 0 ) then
     gamma(2) = dimag( izz(iemax) )
     omega(2) = dreal( izz(iemax) )
  end if

  !Search for most unstable electron mode using varying correlation length
  searchmode = S_ELC

  !Different correlation length for electrons
  zflh = fle
!  zflh = ahyd * fle / zhyd**2.0
  zflz = aimp * fle / zimp**2.0
  kyrho = sqrt(fle)
  wde  = 2.0 * kyrho * csound / rmaj

  call w20solv( ezz, ezvr, ezvi, eiu, eeu, eimax, eemax, G_ave_e, kps_e, ierr )

  ! Record the most signifiant modes
  if ( eimax > 0 ) then
     gamma(3) = dimag( ezz(eimax) )
     omega(3) = dreal( ezz(eimax) )
  end if
  if ( eemax > 0 ) then
     gamma(4) = dimag( ezz(eemax) )
     omega(4) = dreal( ezz(eemax) )
  end if

  do j=1,4
     if ( omega(j).lt.0.0 ) then !Ion mode
        gamma(j) = max( 0D0, gamma(j) - abs(wexb) )
     else
        gamma(j) = max( 0D0, gamma(j)**2-0.25*wexb**2 )
        if ( gamma(j) > 0D0 ) gamma(j) = sqrt(gamma(j))
     end if
  end do


  !Compute transport coefficients

  diffs = 0.0
  vconv = 0.0

  idiffs = 0.0
  ediffs = 0.0

  ivconv = 0.0
  evconv = 0.0

  !G_ave_i = 1.0

  !No unstable modes, then no transport
  if ( iiu .eq. 0 .and. ieu .eq. 0 .and. eiu .eq. 0 .and. eeu .eq. 0 ) then
     return
  end if

  !Consider whether a separate correlation length for electrons is needed
  !Conditions for yes:
  ! Fastest growing mode in electron w20solv call gives more transport than
  ! fastest growing mode in ion w20solv call
  ! (note that unstable modes must exist for the comparison to occur)
  ! No electron modes in ion w20solv call
  !Conditions for no
  ! All others
  if ( eeu > 0 ) then

     if ( ieu > 0 ) then
        if ( dimag(izz(iemax))/sqrt(flh) .lt. dimag(ezz(eemax))/sqrt(fle) ) then
           l_elc=.true.
        else
           l_elc=.false.
        end if
     else
        l_elc=.true.
     end if

  else if ( eeu .eq. 0 ) then
     l_elc=.false.
  end if

  if ( l_elc ) then
    
     !Set up call to flux computation using only
     ! -Eigenvector associated with fastest growing electron mode in electron w20solv call
     ! -Eigenvalue associated with fastest growing electron mode in electron w20solv call
     zz = ezz
     ezz= (0.0,0.0)
     ezz(1) = zz(eemax)
     
     ezvr(1:neq,1) = ezvr(1:neq,eemax) 
     ezvi(1:neq,1) = ezvi(1:neq,eemax) 

     !Use different correlation length for electron modes
     zflh = fle
!     zflh = ahyd * fle / zhyd**2.0
     zflz = aimp * fle / zimp**2.0
     kyrho = sqrt(fle)
     wde  = 2.0 * kyrho * csound / rmaj

     !Compute fluxes
     call w20diff( ezz, ezvr, ezvi, G_ave_i, kps_e, ediffs, evconv )

     !Set up call to flux computation with ion modes
     ! -Remove fastest growing electron mode in ion w20solv call
     if (iemax>0) izz(iemax)=(0.0,0.0)

  end if

  !Use correlation length for ion modes
  kyrho = (dsqrt(flh)-0.1*fte)*epsncorr
  zflh = kyrho**2.0 / zhyd**2.0
!  zflh = ahyd * flh / zhyd**2.0
  zflz = aimp * kyrho**2.0 / zimp**2.0

  wde  = 2.0 * kyrho * csound / rmaj

  !Compute fluxes
  call w20diff( izz, izvr, izvi, G_ave_i, kps_i, idiffs, ivconv )

  !Add transport from two calls
  diffs = idiffs + ediffs
  vconv = ivconv + evconv

End Subroutine w20main!}}}

    subroutine w20init &!{{{
      (  i_print, h_print, &
         z_te,   z_ne,    z_vtor, z_vpol,   z_vpar, &
         z_btor, z_amin,  z_rmin, z_rmaj,   z_eps0, &
         z_aimp, z_ahyd,  z_zimp, &
         z_gte,  z_gti,   z_gtz , &
         z_gne,  z_gni,   z_gnz, &
         z_gvt,  z_gvp,   z_gvpar, &
         z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
         z_q,    z_shear, z_kappa, z_wexb , z_Csound0 )

      implicit none

      real(rkind), intent(in) :: &
           z_te,   z_ne,   z_vtor, z_vpol, z_vpar,  &
           z_btor, z_amin, z_rmin, z_rmaj, z_eps0,  &
           z_aimp, z_ahyd, z_zimp, &
           z_gte,  z_gti,  z_gtz , &
           z_gne,  z_gni,  z_gnz,  &
           z_gvt,  z_gvp,  z_gvpar, &
           z_kyrho,z_tauh,  z_tauz,  z_nz_ne, z_ns_ne, z_fte, &
           z_q,    z_shear, z_kappa, z_wexb, z_Csound0

      integer, intent(in) :: i_print, h_print
      
      lprint = i_print
      hprint = h_print
      te     = z_te
      ne     = z_ne*1.0E-19_rkind
      btor   = z_btor
      amin   = z_amin
      rmin   = z_rmin
      rmaj   = z_rmaj
      eps0   = z_eps0
      aimp   = z_aimp
      ahyd   = z_ahyd
      zimp   = z_zimp
      kyrho  = z_kyrho
      tauh   = z_tauh
      tauz   = z_tauz
      nz_ne  = z_nz_ne
      ns_ne  = z_ns_ne
      fte    = z_fte
      q      = z_q
      shear  = z_shear
      kappa  = z_kappa
      Csound0= z_Csound0 
      
      !Initialize some additional variables
      zhyd  = 1.0e0_rkind
      zfion = 1.0e0_rkind

      !Aspect ratio and major radius at axis

      rmaj0 = amin / eps0
      eps   = max( 0.01e0_rkind, rmin / rmaj0 )
!      rmaj  = rmaj0 * ( 1.0e0_rkind+1.084e0_rkind*eps )

      !Gradients
      gte    = sign( max( abs(z_gte), grdmin ), z_gte ) * rmaj / rmaj0
      gti    = sign( max( abs(z_gti), grdmin ), z_gti ) * rmaj / rmaj0
      gtz    = sign( max( abs(z_gtz), grdmin ), z_gtz ) * rmaj / rmaj0
      gne    = sign( max( abs(z_gne), grdmin ), z_gne ) * rmaj / rmaj0
      gni    = sign( max( abs(z_gni), grdmin ), z_gni ) * rmaj / rmaj0
      gnz    = sign( max( abs(z_gnz), grdmin ), z_gnz ) * rmaj / rmaj0
      gvp    = sign( max( abs(z_gvp), grdmin ), z_gvp ) * rmaj / rmaj0
      gvt    = sign( max( abs(z_gvt), grdmin ), z_gvt ) * rmaj / rmaj0
      gvpar  = sign( max( abs(z_gvpar), grdmin ), z_gvpar ) * rmaj / rmaj0

!      gte    = max(z_gte, grdmin ) * rmaj / rmaj0
!      gti    = max(z_gti, grdmin ) * rmaj / rmaj0
!      gtz    = max(z_gtz, grdmin ) * rmaj / rmaj0
!      gne    = max(z_gne, grdmin ) * rmaj / rmaj0
!      gni    = max(z_gni, grdmin ) * rmaj / rmaj0
!      gnz    = max(z_gnz, grdmin ) * rmaj / rmaj0
!      gvp    = max(z_gvp, grdmin ) * rmaj / rmaj0
!      gvt    = max(z_gvt, grdmin ) * rmaj / rmaj0
!      gvpar  = max(z_gvpar, grdmin ) * rmaj / rmaj0


      !Effective magnetic shear
      shat = sqrt(max(2*shear-1.0e0_rkind+(kappa*(shear-1))**2,0.0))

      !Temperatures and temperature ratios
      th = te*tauh
      tz = te*tauz

      tau      = te  / max( th, 0.001e0_rkind )    ! Te/Ti
      tau_inv  = th  / max( te, 0.001e0_rkind )    ! Te/Ti

      ztauh = tauh / zhyd               ! T_h / ( Z_h Te )
      ztauz = tauz / zimp               ! T_z / ( Z_z Te )

      !Correlation lengths (k_y rho )**2
      epsncorr = min( 2.0e0_rkind / abs(gne), 4.0e0_rkind )**0.17e0_rkind
      flh   = (0.7e0_rkind+2.4e0_rkind/(7.14e0_rkind*q*shat+0.1e0_rkind)) &
              *kyrho**2.0e0_rkind*2.0e0_rkind/(1.0e0_rkind+tau_inv)
!      zflh  = flh
      kyrho = (sqrt(flh)-0.1e0_rkind*fte)*epsncorr
      zflh  = kyrho**2.0e0_rkind / zhyd**2.0e0_rkind
      zflz  = aimp * kyrho**2.0e0_rkind / zimp**2.0e0_rkind    ! A_z * ( k_y rho ) / Z_z**2.0

      znz_ne = z_zimp*nz_ne
      zns_ne = zfion*ns_ne

      !Densities and density ratios
      nh_ne = (1.0 - znz_ne - zns_ne )
      ne_nh = 1.0 / nh_ne
      nh    = ne * nh_ne
      nz    = ne * nz_ne
      ns    = ne * ns_ne

      !Effective charge Zeff
      zeff = (nh_ne*zhyd**2 + nz_ne*zimp**2 + ns_ne*zfion**2 ) 

      fft   = (1.0e0_rkind - znz_ne - zns_ne ) / ( 1.0e0_rkind - fte )
      fzft  = znz_ne / ( 1.0e0_rkind - fte )
      ftrt  = ( 5.0e0_rkind/3.0e0_rkind ) * tauh

      !Collision related quantities
      bt    = 1.5e0_rkind
      bt1   = bt - 2.5e0_rkind

      !Plasma beta
      beta = 4.027E-22_rkind * ( z_ne*te*(1.0e0_rkind + nh_ne*tauh + nz_ne*tauz ) ) / (btor**2)

      !Electron beta
!      betae = 0.01D0*0.4092D0*ne*te/(btor**2)
      betae = 4.027E-22_rkind * z_ne*te / (btor**2)      

      !MHD alpha
      alpha_MHD = em*(q**2.0e0_rkind)*betae*( gne + gte + tauh*(gni + gti) )

      !Speed of sound
      Csound = sqrt( te * 1000.0 * 1.602E-19_rkind / ( 1.6725E-27_rkind * ahyd ) )
!      Csound = 0.311*sqrt(te)*10**6.0/sqrt(ahyd)

      !Ion gyroradius
      rhos = Csound / ( 1.e8_rkind*0.957e0_rkind*btor )
!      rhos = (4.57E-3 * sqrt(ahyd * te*tauh)) / btor 

      !Temporary: use Csound at axis for momentum transport stuff
!      Csound = Csound0

      !Toroidal velocity normalized by speed of sound
      vtor = z_vtor / Csound0

      !Poloidal velocity normalized by speed of sound
      vpol = z_vpol / Csound0

      !Parallel velocity normalized by speed of sound
      vpar = z_vpar / Csound0

      !Magnetic drift frequency with fixed k_y rho
      wde = 2.0e0_rkind * sqrt(flh) * Csound / rmaj
!      wde = 2.0 * kyrho * Csound / rmaj
!      wde = 2.0 * kyrho * Csound / (rmaj*(1.0+1.084*eps))

      !ExB shear rates
      wexb   = abs(z_wexb / wde)

      !Bohm diffusion coefficient
      diffBohm = rhos * Csound

      !Electron ion collision rate (main ions only)
      vei = 0.09e0_rkind*1.0E4_rkind * nh * ( 15.95e0_rkind - log( sqrt(ne)/te) ) / te**1.5
!      vei = 0.09174*1.0E4 * nh * (15.2 - 0.5*dlog(ne/10) + dlog(te)) / te**1.5
!      vei = (4.0/3.0*sqrt(2.0*3.14159265))*nh*(37.8-dlog(dsqrt(ne)/te))/te**1.5

      !Effective collision rate
      vef = vei / ( eps * wde )
!      vef = vei * q / ( eps * wde )

      !Electron correlation length
      fle = max(0.005,0.024-0.0064*(vei/10000.0-0.61)+0.4*betae)

      gm   = 1.0e0_rkind / (1.0e0_rkind-fte)
      bta  = ftrt*(gm+tauh)
      xt   = 1.0e0_rkind / (1.0e0_rkind+tauh)

!      lprint = 9
      lprint = i_print

    end subroutine w20init!}}}

subroutine w20diff &!{{{
   ( zz, zvr, zvi, G_ave, kps, diffs, conv )
  !use w20data
  implicit none

  integer, parameter :: ichn = 6

  complex*16, intent(in), dimension(neq) :: zz
  real*8, intent(in), dimension(neq,neq) :: zvr, zvi
  real*8, intent(inout), dimension(6) :: diffs, conv
  real*8, intent(in) :: G_ave, kps

  !Internal variables
  integer :: i,j                        !Looping indeces

  !Printing flag
  integer :: iw = 0

  !Pointer to axis point (temporary, must undo all this)
  integer :: ik = 1

  !Unstable mode
  complex*16 :: wu                    !Unstable mode complex frequency
  real*8 :: wu_sq                     !Magnitude of wu squared
  real*8 :: re_wu, im_wu              !Real and imaginary part of unstable mode
  real*8 :: im_zz                     !Imaginary part of unstable mode without ExB shear stabilization
  real*8 :: im_wu_inv                 !Inverse of the growth rate
  real*8 :: wu_res                    !Residual part of omega unstable?
  integer :: iunst                    !Number of unstable modes

  !Some physics quantities
  real*8 :: d1

  !Scaling factor for geometric correction (1/gdro2 in Jetto code)
  real*8 :: gk = 1.0

  !Fractions
  real*8 :: thrd2
  real*8 :: thrd = 0.3333333333

  !Random internal variables of unknown origin, destination, or meaning
  real*8 :: hx

  !Local quantities copied from module and bounded
  real*8 :: sn

  !eta_i threshold
  real*8 :: eith

  !Various quantities related to diffusivities
  real*8, dimension(ichn) :: xi, xe, xd, xq, xdq, xtv, hpt, dm, dmt, &
       chi, che, d, chq, dhq
  real*8 :: hp, shp, shpe, dmi, hpe

  !Local diffusivities and accumulators for diffusivities
  real*8 :: chic, schef, def, dqeff, chqeff, schq, ceft, deft
  real*8 :: schi, gci, sche, gce, sd, gd, sdm, gcq, sdq, gnq, gei
  real*8 :: xih, xqh, xeh, xdh, xhh

  real*8 :: pe, pq, pn, pnq,tiq, pi, piq

  !Possibly, non-diagonal elements and pinches from EM effects
  real*8 :: dms1, dms2, dms21, dms22, dmi21, dmi22, smp, smpt, bv, dmip, adisp, bdisp, &
       dt, dmit, dms, dmst, dmd, dmdn, dmdt, dmdn2, dmdt2, dmi1, dmi2, dadk, dbdk

  !Variables for output of poloidal, toroidal, perpendicular, thermoelectric
  !components of the pinch velocities
  real*8::v_p_tor,v_p_pol,v_p_per,v_p_thr

  !Collision related parameters
  complex*16 :: ga, gb, gm2
  real*8 :: re_ga, im_ga, re_gb, im_gb
  real*8 :: hr, xde, xdi, yda

  !Related to momentum transport
  real*8 :: dmh, dmef, dmef1, dmef2, dmef21, dmef22, dmdiag, dmdiagt, dmeftest
  real*8 :: dmeff, dmtef, dmteff
  real*8 :: sp, smef1, smef21, smef22, sre
  real*8 :: tsour
  real*8 :: kxi
  real*8 :: ainf, imf, tc, dp1, gp1, gp2
  real*8 :: elfs

  !Lengthscales related to momentum transport not yet purified
  real*8 :: lvf, lvft

  !Variables related to the parallel wavenumber in momentum transoprt computation
  real*8 :: Kkap, Hm, kap1, acn, RDexp, ks

  !Magnitude square, real and imaginary parts of various perturbations
  real*8 :: nehat_sq, re_nehat, im_nehat, nzhat_sq, re_nzhat, im_nzhat
  real*8 :: re_eln, im_eln
  real*8 :: dni, dne, dnz, di, de, dz, dnin
  real*8 :: dznt, dzn1, dzn2
  real*8 :: dtimp
  real*8 :: dqn, dqt, dtiq
  real*8 :: phs

  !More unknown variables
  !Random pieces of equations which need variables
  real*8 :: tt, et, ht, h, k, t, ts, t1, t2, a3
  real*8 :: c1, c2, e1, e2, h1, h2, at, ct, a2, a1, a, b, c, e, dh, f
  real*8 :: stf, eqh

  !Divergences of stuff???
  real*8 :: divgat, divat, devgat, devat, diva3, devga1, devga2, devga3, deva1, deva2, deva3
  real*8 :: divga1, divga2, divga3, diva1, diva2, deva, devb
  real*8 :: divga, divgb, diva, divb, devga, devgb

  !Sources of stuff????
  real*8 :: svat, sva1, sva2, sva3, sva, svb

  !Effective collisionality normalized by electron perturbation magnitude square
  real*8 :: vefn

  !Complex variables from Weiland's code
  COMPLEX*16 :: Fm,hf,EMP,EMPT,VPL,Vg,WRES
  COMPLEX*16 :: FIH,NEF,NRAT,ELN,TEF,FRAT,AV,NIF,MRT,HQ,CDexp
  COMPLEX*16 :: RNI,TII,RTI,TICF,CN,AVRAT,MRT1,MRT2,NII
  COMPLEX*16 :: nehat,nzhat,BT2
  COMPLEX*16 :: IU,HC,WZ,WZP
  COMPLEX*16 :: DMSP,KPF,kpf1,kpf2,KPX,ELMS,KPDX,TCONT,NCONT,EMCONT

  wde = 2.0 * csound * kyrho / rmaj0

  !Lengthscales of velocity as defined by J. Weiland
  !Note: these are replaced by gradients (gvt, gvp) in most places

  lvf   = rmaj / ( amin * gvp )
  lvft  = rmaj / ( amin * gvt )

  !Initialize some variables
  diffs = 0.0
!  vconv = 0.0
  iunst = 0
  wu    = ( 0.0, 0.0 ) 

  THRD2=THRD*THRD
  IU=(0.D0,1.D0)
  SHP=0.D0
  CHIC=0.
  SHPE=0.
  SCHEF=0.
  DEF=0.
  EQH=(gtz/gnz)-STR+FTR*(2.0/gnz)
  d1 = 6.616*dsqrt(ahyd)/(rmaj*btor**2)

  !Two missing variables
  RDexp = 1.0

  !Initialize diffusion matrix elements 
  chi(1:6)=0.0 ; che(1:6)=0.0 ; d(1:6)  =0.0 ; chq(1:6)=0.0 
  dm(1:6) =0.0 ; dhq(1:6)=0.0 ; dmt(1:6)=0.0 ; hpt(1:6)=0.0

  !Initialize effective diffusivity
  SCHI=0.D0 ; SCHE=0.D0 ; SD  =0.D0 ; SDM =0.D0 ; SDQ =0.D0 ; XHH =0.D0

  !Main loop over unstable modes
  !(further down is second loop for momentum transport)
  !Thermal and particle diffusivities are computed here

  main_modes: do j=1,neq

     re_wu=DREAL(zz(j))
     im_wu=DIMAG(zz(j))

     !Cycle when growth rate is negative (stable mode)
     if (im_wu.le.0.0) cycle main_modes

     !Notice different substraction of ExB shear for
     !ion modes and electron modes
     if (re_wu.lt.0.0) then
        im_wu = im_wu - dabs(wexb)
     else if ( dreal(zz(j)).gt.0.0) then
        wu_res = im_wu**2-0.25*wexb**2
        if (wu_res .lt. 0.0) cycle
        im_wu = sqrt(wu_res)
     end if

     wu=re_wu+(0.0,1.0)*im_wu

     if(im_wu .le. 1.0E-3) then
        cycle main_modes
     end if

     !Define complex frequency of mode after ExB shear effect
     wu=CMPLX(re_wu,im_wu)
     WRES=wu+FTRT
     gm2=1.0+0.5*gte/(wu-1.0+IU*vef)
     BT2=BT-2.5*GM2
     HC=wu-FTR+TVR*BT1
     nehat=wu*wu-2.D0*FTR*wu+FTR+IU*VEF*HC
     wu_sq=re_wu*re_wu+im_wu*im_wu
     im_wu_inv=1.0/im_wu
     !
     !   ******   COLLISION PARAMETERS  *******
     !
     GA=wu-FTR+TVR*BT2
     GB=(wu-FTR)/(wu-1.D0+IU*VEF)
     re_ga=DREAL(GA)
     im_ga=DIMAG(GA)
     re_gb=DREAL(GB)
     im_gb=DIMAG(GB)

     HR=re_wu-FTR+TVR*BT1
     XDE=wu_sq-FTR*re_wu
     XDI=wu_sq+FTR*tau_inv*re_wu
     YDA=re_wu*(1.D0-(2.0/gne))+(gte/gne)-STR
     !   ***************************************
     !

     re_nehat=DREAL(nehat)
     im_nehat=DIMAG(nehat)
     nehat_sq=(re_nehat)**2+(im_nehat)**2
     !
     !      IF(VEF.EQ.0.) GO TO 25
     DIVGA=FTRT*((2.0/gne)*(re_ga*im_nehat-im_ga*re_nehat) &
          -YDA*im_wu+im_wu*HR*(1.-(2.0/gne)))

     DIVGB=FTRT*(re_gb*im_nehat-im_gb*re_nehat-im_wu)

     DIVA=XDI*((2.0/gne)*(re_ga*re_nehat+im_ga*im_nehat)-im_wu*im_wu*(1.-(2.0/gne))-YDA*HR)

     DIVB=XDI*(re_gb*re_nehat+im_gb*im_nehat-HR)

     DEVGA=BT1*((YDA-VEF*im_ga*(2.0/gne))*im_nehat    &
          -(im_wu*(1.-(2.0/gne))+VEF*re_ga*(2.0/gne)) &
          *re_nehat) &
          +FTR*(im_wu*YDA+(2.0/gne)*(im_ga*re_nehat-re_ga*im_nehat)-im_wu*HR*(1.-(2.0/gne)))

     DEVGB=BT1*((1.-VEF*im_gb)*im_nehat-VEF*re_gb*re_nehat) &
          +FTR*(im_wu+im_gb*re_nehat-re_gb*im_nehat)

     DEVA=XDE*((2.0/gne)*(re_ga*re_nehat+im_ga*im_nehat)-im_wu*im_wu*(1.-(2.0/gne))-YDA*HR) &
          -BT1*(re_wu-FTR)*((YDA-VEF*im_ga*(2.0/gne))*re_nehat+(im_wu*(1.-(2.0/gne)) &
          +VEF*re_ga*(2.0/gne))*im_nehat)

     DEVB=XDE*(re_gb*re_nehat+im_gb*im_nehat-HR)-BT1*(re_wu-FTR)*((1.- &
          VEF*im_gb)*re_nehat+VEF*re_gb*im_nehat)

     SVA=(2.0/gne)*(re_ga*re_nehat+im_ga*im_nehat)-im_wu*im_wu*(1.-(2.0/gne))-YDA*HR

     SVB=re_gb*re_nehat+im_gb*im_nehat-HR

     VEFN=VEF/nehat_sq

     DIVGA=DIVGA*VEFN
     DIVGB=DIVGB*VEFN
     DEVGA=DEVGA*VEFN
     DEVGB=DEVGB*VEFN
     DIVA=DIVA*VEFN
     DIVB=DIVB*VEFN
     DEVA=DEVA*VEFN
     DEVB=DEVB*VEFN
     SVA=SVA*VEFN
     SVB=SVB*VEFN
     !   25 CONTINUE
     !
     STF=STR-FTR*(2.0/gne)
     A1=wu_sq*((2.0/gne)-1.D0)+2.D0*re_wu*STF
     A=wu_sq*(A1+FTR*(STR*(2.0/gne)-11.D0*THRD-tau_inv*(1.D0- &
          FTR*(2.0/gne))))+FTR*FTR*tau_inv*(2.D0*re_wu*(1.D0-(2.0/gne))-STF)
     A=A/nehat_sq
     B=(wu_sq*(2.D0*(re_wu-FTR)+FTR*tau_inv)- &
          FTR*FTR*tau_inv)/nehat_sq
     C=wu_sq*(A1+TVR*FTR*((2.0/gne)-4.D0))+ &
          FTR*FTR*(2.D0*re_wu*((2.0/gne)-1.D0)+STF)
     C=C/nehat_sq
     DH=(wu_sq*(2.D0*re_wu-5.D0)+FTR*FTR)/nehat_sq
     E=wu_sq*(1.D0-(2.0/gne))-2.D0*re_wu*STF+FTR*(11.D0*THRD &
          -STR*(2.0/gne))
     E=E/nehat_sq
     F=2.D0*(-re_wu+FTR)/nehat_sq
     !
     re_nzhat=re_wu**2-im_wu**2+2.D0*FTR*ztauz*re_wu+FTR*ztauz*ztauz
     im_nzhat=2.D0*im_wu*(re_wu+FTR*ztauz)
     nzhat_sq=re_nzhat**2+im_nzhat**2
     !
     !   ****   SPLITTING IN  EN  AND EE  ************
     !
     !
     !  Ion thermal conductivity
     !
     DIVGA1=FTRT*((STR-re_wu)*im_wu+im_wu*HR)*VEFN
     DIVGA2=FTRT*(re_ga*im_nehat-im_ga*re_nehat+re_wu*im_wu-im_wu*HR)*VEFN
     DIVGA3=-FTRT*im_wu*VEFN

     DIVA1=XDI*(-im_wu*im_wu+(STR-re_wu)*HR)*VEFN
     DIVA2=XDI*(re_ga*re_nehat+im_ga*im_nehat+im_wu*im_wu+re_wu*HR)*VEFN
     DIVA3=-XDI*HR*VEFN

     !
     !   Electron thermal conductivity
     !
     DEVGA1=(BT1*((re_wu-STR)*im_nehat-im_wu*re_nehat)+FTR*(im_wu*(re_wu-STR)-im_wu*HR))&
          *VEFN
     DEVGA2=(BT1*((-re_wu-VEF*im_ga)*im_nehat+(im_wu-VEF*re_ga)*re_nehat)&
          +FTR*(-im_wu*re_wu+im_ga*re_nehat-re_ga*im_nehat+im_wu*HR))*VEFN
     DEVGA3=(BT1*im_nehat+FTR*im_wu)*VEFN

     DEVA1=(XDE*(-im_wu*im_wu-(re_wu-STR)*HR)-BT1*(re_wu-FTR)*((re_wu-STR)*re_nehat&
          +im_wu*im_nehat))*VEFN
     DEVA2=(XDE*(re_ga*re_nehat+im_ga*im_nehat+im_wu*im_wu+re_wu*HR)-BT1*(re_wu-FTR)&
          *((-re_wu-VEF*im_ga)*re_nehat-(im_wu-VEF*re_ga)*im_nehat))*VEFN
     DEVA3=(-XDE*HR-BT1*(re_wu-FTR)*re_nehat)*VEFN

     SVA1=(-im_wu*im_wu-(re_wu-STR)*HR)*VEFN
     SVA2=(re_ga*re_nehat+im_ga*im_nehat+im_wu*im_wu+re_wu*HR)*VEFN
     SVA3=-HR*VEFN

     !------------------------------------------------------------------------ 
     ! Contributions to the ion conductivity 
     !------------------------------------------------------------------------ 

     A2=(wu_sq*wu_sq+FTR*((STR+FTR*tau_inv-2.D0*re_wu)*wu_sq+FTR*tau_inv*(FTR-&
          2.D0*re_wu)))/nehat_sq
     A3=(wu_sq*(-wu_sq+2.D0*STR*re_wu-FTR*(11.D0*THRD+tau_inv))+FTR*FTR&
          *tau_inv*(2.D0*re_wu-STR))/nehat_sq

     !------------------------------------------------------------------------ 
     ! Contributions to the electron conductivity 
     !------------------------------------------------------------------------ 

     C1=(nehat_sq-wu_sq*wu_sq+14.D0*THRD*wu_sq*re_wu-40.D0*THRD2*wu_sq&
          -50.D0*THRD2*re_wu+175.D0/27.D0)/nehat_sq
     C2=(wu_sq*wu_sq-10.D0*THRD*wu_sq*re_wu+10.D0*THRD2*wu_sq&
          +50.D0*THRD2*re_wu-125.D0/27.D0)/nehat_sq

     !------------------------------------------------------------------------ 
     !   Contributions to the electron diffusivity 
     !------------------------------------------------------------------------ 

     E1=(wu_sq+11.D0*THRD*FTR-2.D0*re_wu*STR)/nehat_sq
     E2=-(-2.D0*re_wu*FTR+(wu_sq+STR*FTR))/nehat_sq

     !------------------------------------------------------------------------ 
     ! Contributions to the main ion Conductivity 
     !------------------------------------------------------------------------ 

     H1=(wu_sq*(-wu_sq-2.D0*ztauz*re_wu*STR&
          +FTR*ztauz*ztauz*(-11.D0*THRD)+FTR*tau_inv*ztauz)&
          +FTR*FTR*tau_inv*ztauz*ztauz*(2.D0*re_wu+STR*ztauz))/nzhat_sq

     H2=(wu_sq*(wu_sq+2.D0*ZTAUZ*re_wu*FTR+FTR*ZTAUZ*ZTAUZ*STR-FTR*tau_inv*ZTAUZ*FTR)&
          -FTR*FTR*tau_inv*ZTAUZ*ZTAUZ*(2.D0*re_wu+FTR*ZTAUZ))/nzhat_sq

     !------------------------------------------------------------------------ 
     ! Contributions to the impurity conductivity 
     !------------------------------------------------------------------------ 

     T1=(wu_sq*(-wu_sq-2.D0*ZTAUZ*re_wu*STR-FTR*ZTAUZ*ZTAUZ*8.D0*THRD)&
          +FTR*FTR*ZTAUZ**3*(2.D0*re_wu+ZTAUZ*STR))/nzhat_sq
     T2=(wu_sq*(wu_sq+2.D0*ZTAUZ*re_wu*FTR+FTR*ZTAUZ*ZTAUZ*TVR)&
          -FTR*FTR*ZTAUZ**3*(2.d0*re_wu+ZTAUZ*FTR))/nzhat_sq

     !------------------------------------------------------------------------ 
     ! Contributions to the impurity diffusivity 
     !------------------------------------------------------------------------ 

     DZN1=(-re_nzhat+2.D0*(re_wu+ZTAUZ*STR)*(re_wu+FTR*ZTAUZ))/nzhat_sq
     DZN2=(re_nzhat-2.D0*(re_wu+ZTAUZ*FTR)*(re_wu+FTR*ZTAUZ))/nzhat_sq
     !
     DIVGAT=DIVGA1+(2.0/gne)*DIVGA2+(gte/gne)*DIVGA3
     DIVAT=DIVA1+(2.0/gne)*DIVA2+(gte/gne)*DIVA3
     DEVGAT=DEVGA1+(2.0/gne)*DEVGA2+(gte/gne)*DEVGA3
     DEVAT=DEVA1+(2.0/gne)*DEVA2+(gte/gne)*DEVA3
     SVAT=SVA1+(2.0/gne)*SVA2+(gte/gne)*SVA3
     AT=(2.0/gne)*A2+A3
     CT=C1-1.+(2.0/gne)*C2
     ET=E1+(2.0/gne)*E2
     HT=H1+(2.0/gnz)*H2
     TT=T1+(2.0/gnz)*T2
     DZNT=DZN1+(2.0/gne)*DZN2
     !
     ! **** IMPURITIES *****
     !
     nzhat_sq=(re_wu*(re_wu+2.*FTR*ZTAUZ)-im_wu*im_wu+FTR*ZTAUZ*ZTAUZ)**2 &
          +4.*im_wu*im_wu*(re_wu+FTR*ZTAUZ)**2
     DTIMP=(wu_sq*(wu_sq*((2.0/gnz)-1.)+2.*ZTAUZ*re_wu*EQH+FTR*ZTAUZ &
          *ZTAUZ*(2.*(gtz/gnz)-11./3.+STR*(2.0/gnz))+FTR*tau_inv*ZTAUZ*(1. &
          +(gtz/gnz)-FTR*(2.0/gnz)))+FTR*tau_inv*(2.*FTR*ZTAUZ*ZTAUZ*re_wu*(1.-(2.0/gnz)) &
          -FTR*ZTAUZ**3*EQH))/nzhat_sq
     !  *************
     !
     DNI=(re_wu+FTR*tau_inv)**2+im_wu*im_wu
     DNE=(re_wu-FTR)**2+im_wu*im_wu
     DNZ=(re_wu+FTR*ZTAUZ)**2+im_wu*im_wu
     DI=(gti/gni)*kyrho*DNI
     DE=(gte/gne)*kyrho*DNE
     DZ=(gtz/gnz)*kyrho*DNZ
     !
     ! **** IMPURITIES *****
     !
     nzhat_sq=re_nzhat**2+im_nzhat**2
     !
     H=(wu_sq*(wu_sq*((2.0/gnz)-1.)-2.*ZTAUZ*re_wu*(STR-FTR*(2.0/gnz)) &
          +FTR*ZTAUZ*ZTAUZ*(-11.*THRD+STR*(2.0/gnz))+FTR*tau_inv*ZTAUZ*(1.-FTR &
          *(2.0/gnz)))+FTR*FTR*tau_inv*ZTAUZ*ZTAUZ*(2.*re_wu*(1.-(2.0/gnz))+(STR-FTR &
          *(2.0/gnz))*ZTAUZ))/nzhat_sq
     K=ZTAUZ*(FTR*FTR*tau_inv*ZTAUZ*ZTAUZ-wu_sq*(2.*re_wu &
          +FTR*(2.*ZTAUZ+tau_inv)))/nzhat_sq
     !
     !  *************
     !
     T=(wu_sq*(wu_sq*((2.0/gnz)-1.)-2.*ZTAUZ*re_wu*(STR-FTR*(2.0/gnz)) &
          +FTR*ZTAUZ*ZTAUZ*(-8.*THRD+TVR*(2.0/gnz)))+FTR*FTR*(ZTAUZ)**3 &
          *(2.*re_wu*(1.-(2.0/gnz))+ZTAUZ*(STR-FTR*(2.0/gnz))))/nzhat_sq
     TS=ZTAUZ*(FTR*FTR*(ZTAUZ)**3-wu_sq*(2.*re_wu+5.*ZTAUZ))/nzhat_sq
     !
     T=(GNZ/GNE)*T
     TS=(GNZ/GNE)*TS
     DQN=((2.0/gnz)-1.)*re_nzhat+2.*((1.-(2.0/gnz))*re_wu+ZTAUZ*(STR-FTR*(2.0/gnz)) &
          )*(re_wu+FTR*ZTAUZ)
     DQT=2.*ZTAUZ*(re_wu+FTR*ZTAUZ)
     DQN=DQN/nzhat_sq
     DQT=DQT/nzhat_sq
     !
     DTIQ=H-(gtz/gnz)*K
     !      DTQQ=T-EQ*TS
     !
     PHS=(re_wu-FTR)*DREAL(BT2)+im_wu*DIMAG(BT2)
     XDH=D1*(TE**1.5D0)*im_wu**3/kyrho
     XIH=XDH/DNI
     XEH=XDH/DNE
     XQH=XDH/DNZ

     !----------------------
     ! Ion transport channel
     !----------------------
     ! Diffusive part
     !
     XI(1)=XIH
     XI(2)=tvr*fte*xih*ne_nh*(B-DIVGA3-DIVGB-im_wu_inv*(DIVA3+DIVB))
     XI(3)=-tvr*xih*(gni/gne+fte*ne_nh*(A3+DIVGA1+DIVA1*im_wu_inv))
     !      XI(4)=-XIH*TVR*nz_ne*Z*ne_nh*K*TAUZ*tau_inv
     XI(4)=0.
     XI(5)=tvr*xih*zimp*h1/nh
     XI(6)=0.d0
     !
     ! Convective part
     !
     PI=fte*TVR*XIH*ne_nh*(A2+DIVGA2+im_wu_inv*DIVA2)*amin/rmaj0
     PIQ=-TVR*XIH*zimp*H2*nz_ne*ne_nh*amin/rmaj0
     HP=XIH*ne_nh*tau_inv*TVR*FTR*(1.D0-fte)*amin/rmaj0
     !
     ! Effective diffusivity
     !
     GCI = XI(1) &
          + XI(2)*gte/gti &
          + XI(3)*gne/gti &
          + XI(5)*nz*gne/gti &
          - (HP+PI+PIQ)*2.0*gk*rmaj/amin/gti

     !-----------------------------------
     ! Electron thermal transport channel
     !-----------------------------------
     ! Diffusive part
     !
     XE(1)=0.D0
     XE(2)=fte*XEH*(1.D0+TVR*(DH-DEVGB-DEVGA3-im_wu_inv*(DEVB+DEVA3)))
     XE(3)=-TVR*fte*XEH*(C1+DEVGA1+im_wu_inv*DEVA1)
     XE(4)=0.D0
     XE(5)=0.D0
     XE(6)=0.D0
     !
     ! Convective part
     !
     PE=fte*TVR*XEH*(C2+DEVGA2+im_wu_inv*(DEVA2+VEF*PHS))*amin/rmaj0
     !
     ! Effective diffusivity
     !
     GCE = XE(2) &
          + XE(3)*gne/gte &
          - PE*gk*2.0*rmaj/amin/gte

     !------------------------------------
     ! Electron particle transport channel
     !------------------------------------
     ! Diffusive part
     !
     XD(1)=0.D0
     XD(2)=-fte*XDH*(ne/te)*(F+im_wu_inv*(SVB+SVA3))
     XD(3)=fte*XDH*(E1-im_wu_inv*SVA1)
     XD(4)=0.D0
     XD(5)=0.D0
     XD(6)=0.D0
     !
     ! Convective part
     !
     PN=-fte*XDH*(E2-im_wu_inv*SVA2)*amin/rmaj0
     !
     ! Effective diffusivity
     !
     GD = XD(2)*te/ne*gte/gne &
          + XD(3) &
          -  PN*gk*2.0*rmaj/amin/gne

     !-----------------------------------
     ! Impurity thermal transport channel
     !-----------------------------------
     ! Diffusive part
     !
     XQ(1)=0.D0
     XQ(2)=0.D0
     XQ(3)=0.D0
     XQ(4)=1.D0
     XQ(5)=0.D0
     XQ(6)=0.D0
     !
     ! Convective part
     !
     PQ=TVR*XQH*T2*2/rmaj
     !
     ! Effective diffusivity
     !
     !      GCQ=XQ(4)+XQ(5)*NQ/(TZ*(gtz/gnz))-PQ*GK*2.*LTZ

     GCQ = -XDH*DQT*NZ/TZ !XDQ(4)
     !------------------------------------
     ! Impurity particle transport channel
     !------------------------------------
     ! Diffusive part
     !
     XDQ(1)=0.D0
     XDQ(2)=0.D0
     XDQ(3)=0.D0
     XDQ(4)=-XDH*DQT*NZ/TZ
     XDQ(5)=XDH*DZN1
     XDQ(6)=0.D0
     !
     ! Convective part
     !
     PNQ=-XDH*DZN2*2/rmaj
     !
     ! Effective diffusivity
     !
     GNQ = XDQ(4)*TZ*(gtz/gnz)/NZ &
          + XDQ(5) &
          - PNQ*GK*rmaj/gnz

     XTV(1)=0.D0
     XTV(2)=0.D0
     XTV(3)=0.D0
     XTV(4)=0.D0
     XTV(5)=0.D0
     XTV(6)=1.D0
     !
     TIQ=1.D0/(TAUZ*(GNI/GNZ))

     SHP=SHP+HP
     !
     HPT(1)=HPT(1)+HP+PI+PIQ
     HPT(2)=HPT(2)+PE
     HPT(3)=HPT(3)+PN
     !      HPT(4)=HPT(4)+PQ
     HPT(5)=HPT(5)+PNQ
     !
     SCHI=SCHI+GCI
     SCHE=SCHE+GCE
     SD=SD+GD
     SDM=SDM+GCQ
     SDQ=SDQ+GNQ
     !
     !
     !      CHQEFF=D1*TE**1.5*WI**3*((gtz/gnz)-TVR-TVR*DTQQ)/DQ
     chqeff=0.0
     schq  =0.0
     DQEFF=XDH*(DQN-(gtz/gnz)*DQT)

     XHH=XHH+XIH

     Do i=1,ichn
        CHI(I)=CHI(I)+XI(I)
        CHE(I)=CHE(I)+XE(I)
        D(I)=D(I)+XD(I)
        !      CHQ(I)=CHQ(I)+XQ(I)
        DHQ(I)=DHQ(I)+XDQ(I)
     end do

  End Do main_modes
  ! 
  !  IF electromagnetic effects or collisions on free electrons are included
  !  the transport coefficients are corrected for this.
  !------------------------------------------------------------
  DT=D1*TE**1.5
  SHPE=0.D0
  SCHEF=0.D0
  DEF=0.D0
  DMI=0.D0
  TSOUR=0.D0
  DMIT=0.D0
  DMS=0.D0
  DMST=0.D0
  DMSP=(0.D0,0.D0)
  DMD=0.D0
  DMDT=0.D0
  DMDT2=0.D0
  DMI1=0.D0
  DMI2=0.D0
  DMS1=0.D0
  DMS2=0.D0
  DMI21=0.D0
  DMI22=0.D0
  SMP=0.D0
  SMPT=0.D0
  BV=0.D0
  DMIP=0.D0    !!  Used for summing up curvature pinch terms
  v_p_pol=0.0
  v_p_tor=0.0
  v_p_per=0.0
  v_p_thr=0.0
  adisp=((2.0/gne)-1.D0-flh*(FTRT*(2.0/gne)-tau_inv*(1.D0+(gti/gni))))/(1.D0+flh)
  bdisp=tau_inv*(2.0/gne)/(1.D0+fte)
  dadk=-2.D0*kyrho*rhos*(FTRT*(2.0/gne)-tau_inv*(1.D0+(gti/gni))+((2.0/gne)-1.D0)/(1.D0+flh))
  dadk=dadk/(1.D0+flh)
  dbdk=-2.D0*(2.0/gne)*tau_inv*kyrho*rhos/(1.D0+flh)**2
  EITH=TVR+10.D0*(2.0/gne)*tau_inv/9.D0
  KXI=2.D0*q*kyrho/(eps*rmaj0)
  !------------------------------------------- loop

  main_momentum: DO J=1,NEQ

     re_wu=DREAL(zz(j))
     im_wu=DIMAG(zz(j))

     if ( im_wu.lt.0.0 ) cycle

     if ( re_wu.lt.0.0 ) then !Ion mode
        im_wu = im_wu - dabs(wexb)
     else if ( re_wu .gt.0.0 ) then
        wu_res =im_wu**2-0.25*wexb**2
        if (wu_res.lt.0.0) cycle
        im_wu = dsqrt(wu_res)
     end if

     wu=re_wu+(0.0,1.0)*im_wu

     if(im_wu .le. 1.0E-2) then
        cycle
     end if

     wu=CMPLX(re_wu,im_wu)

     im_zz = max( 0.01, dimag(zz(j)) )

     ! --- contr. to chii from em free electr. ----
     !
     DNI=(re_wu+FTR*tau_inv)**2+im_wu*im_wu
     DNIN=(re_wu+FTR*tau_inv)**2+im_zz*im_zz
     wu_sq=re_wu*re_wu+im_wu*im_wu
     XDH=DT*im_wu**3/kyrho
     XIH=XDH/DNI
     FIH=DCMPLX(ZVR(1,J),ZVI(1,J))
     ! F. Halpern: place a lower bound on the absolute value of FIH 
     IF ( abs(fih) .lt. 0.001 ) then
        fih = (0.001,0.0)
     end if

     !-- Here NEF for disp 9 is defined --
     AV=DCMPLX(ZVR(8,J),ZVI(8,J))
     NEF=FIH-(ZZ(J)-0.5*gne)*AV/KPC
     !--------------------------------------------

     IF(ABS(FIH).LT.0.0001) FIH=(0.0001,0.)
     NRAT=NEF/FIH
     ELN=NRAT-1.
     re_eln=DREAL(ELN)
     im_eln=DIMAG(ELN)
     AINF=TVR*(FTR*tau_inv*re_eln+im_eln*(im_wu*im_wu+re_wu*(re_wu+FTR*tau_inv))/im_zz)
     HPE=XIH*ne_nh*(1.-fte)*eps0*AINF

     SHPE=SHPE+HPE
     ! ****  Free electron heat flux *********************
     !
     !-----------------------------------------------------
     IF(NEQ.EQ.11) GO TO 10099
     IF(NEQ.EQ.9) GO TO 10098
     !--- Here TEF for disp10 is defined ---
     TEF=DCMPLX(ZVR(6,J),ZVI(6,J))
     GO TO 10099
10098 CONTINUE
     !-- Here TEF for disp9 is defined ---
     !      TEF=(gte/gne)*0.5*gne*AV/KPC
     TEF = 0.5*gte*AV/KPC
     !-------------------------------------------------------------
10099 CONTINUE
     FRAT=TEF/FIH
     IMF=-IMAG(FRAT)*DNIN/DNI
     CEFT=(1.-fte)*IMF*DT*(2.0/gte)*im_wu**3/(kyrho*im_zz)
     SCHEF=SCHEF+CEFT
     !**********************************************************
     ! ----Free electron particle flux -----------
     GEI=-IMAG(NRAT)/im_wu
     DEFT=(1.D0-fte)*GEI*(2.0/gne)*XDH
     DEF=DEF+DEFT
     ! ----Momentum transport ------
     TII=DCMPLX(ZVR(2,J),ZVI(2,J))
     NIF=DCMPLX(ZVR(3,J),ZVI(3,J))
     AVRAT=AV/FIH
     MRT=(TII+NIF)/FIH
     MRT1=0.5*(gti*gne/gni-gne*tvr)/(wu+ftrt*G_ave)
     TICF=TII/FIH
     NII=NIF/FIH
     MRT2=NIF*(1.+TVR*wu/(wu+FTRT*G_ave))/FIH
     DMS21=im_wu*im_wu*tau_inv*DREAL(MRT1)
     DMS22=im_wu*im_wu*tau_inv*DREAL(MRT2)
     DMS1=im_wu*im_wu
     DMS2=im_wu*im_wu*tau_inv*DREAL(MRT)
     DMS=DMS1+DMS2
     DMI1=DMI1+DMS1
     DMI2=DMI2+DMS2
     DMI21=DMI21+DMS21    !! part of temp. pert 
     DMI22=DMI22+DMS22    !! density pert including part from temp. pert.
     DMI=DMI+DMS
     DMD=DMD+im_wu**3/wu_sq
     !Comparison between old and new diagonal element
     DMDN=(re_wu+2.*tau_inv*G_ave)**2+im_wu*im_wu
     DMDT=DMDT+im_wu**3/DMDN  !! Diagonal element for toroidal momentum transport
     TC=(re_wu+2.*tau_inv*G_ave)/DMDN  !! Correl. time for TEP and Termoel. mom pinch (Hahm)
     DP1=2.*DT*im_wu*im_wu/kyrho

! By F.Halpern: 18-Feb-2009: Introduce gp1 and gp2 as stated in Jan Weiland's MOMCODE from 25-Jan-2009

! By A Kritz: 31-Mar-2009: Moved line defining ELMS before it is used in computing gp2
     ELMS = (wu+0.5*tau_inv*(gne+gti*gne/gni))*AVRAT/KPC ! From Weiland's momcode

     gp1 = - im_wu * dp1 * vtor * G_ave / dmdn
     gp2 = dreal( IU * ( MRT - tau*em2*ELMS )/( wu + 2.0*tau_inv*G_ave ) ) * dp1 * vtor * G_ave

!     GP1=0.5*TC*DP1*vtor  !!  TEP  momentum pinch
!     GP2=DREAL(TICF/(wu+2.*tau_inv*G_ave))*DP1*vtor  !!  Termoel. momentum pinch


     ELFS = (1.D0+0.5*tau_inv*(gne+gti*gne/gni)/wu)*AVRAT/KPC ! From Weiland's momcode

! By F. Halpern: Add factor of average curvature as per 25-Jan-2009 MOMCODE
     DMSP=im_wu*im_wu*(IU/(wu+2.*tau_inv*G_ave))*(1+tau_inv*MRT-em2*ELMS) !! due to v_par incl TC
!     DMSP=im_wu*im_wu*(IU/(wu+2.*tau_inv))*(1+tau_inv*MRT-em2*ELMS) !! due to v_par incl TC

     !----- Computation of parallel wavenumber k_\parallel

     !K in Weiland's 2006 EPS paper
     !F. Halpern: replaced toroidal velocity gradient with parallel gradient
     !      Kkap = rmaj0/rmaj*q*kyrho*vtor*gvt
     Kkap = rmaj0/rmaj*q*kyrho*vpar*gvpar

     !\kappa_1 in Weiland's 2006 EPS paper
     !F. Halpern: replaced toroidal velocity gradient with parallel gradient

     kap1 = vpol*gvp/(gne*sign(max(abs(shear),0.1),shear)*kyrho)
     !      kap1 = vpol*gvpar/(gne*sign(max(abs(shear),0.1),shear)*kyrho*dsqrt(kappa))
     kap1 = sign(min(abs(kap1),1.0),kap1)

     !ks in Weiland's Jan-2009 MOMCODE
     ks   = -2.0*vtor*kyrho*q/(tau_inv*tau_inv)

     Fm=wu*(1.+FTRT)+tau_inv*0.5*(gti*gne/gni-TVR*gne)
     Hm=FTRT*(1.+tau_inv)
     hf=4.0*flh*q**2*wu*0.5*gne

     !Flux average of k_\parallel normalized, plus separate expressions 
     !poloidal and toroidal driving terms
     KPF=-(0.5*(wu+ftr)*(Kkap+ks) + IU*hf*im_wu*KAP1)/((Fm+Hm)*Q*Rmaj0) !!This is dimk_par
     kpf1=-(0.5*(wu+ftr)*(Kkap+ks))/((Fm+Hm)*Q*Rmaj0)
     kpf2=-(IU*hf*im_wu*KAP1)/((Fm+Hm)*Q*Rmaj0)

!     KPF=-(0.5*(wu+ftr)*Kkap+IU*hf*im_wu*KAP1)/((Fm+Hm)*Q*Rmaj0) !!This is dimk_par
!     kpf1=-(0.5*(wu+ftr)*Kkap)/((Fm+Hm)*Q*Rmaj0)

     !Normalization is removed here
     KPF=KPF*CSound/Wde
     kpf1=kpf1*CSound/Wde
     kpf2=kpf2*CSound/Wde

     !Compute the fluxes here
     !Split the different pinches due to velocity shear
     ! v_p_per ---> perpendicular part
     ! v_p_pol ---> poloidal part of parallel velocity pinch
     ! v_p_tor ---> toroidal part of parallel velocity pinch

     v_p_per = v_p_per -2.D0*d1*te**1.5*(eps/q)*dms
     v_p_pol = v_p_pol +2.D0*d1*te**1.5*dreal(DMSP*kpf2)/kyrho
     v_p_tor = v_p_tor +2.D0*d1*te**1.5*dreal(DMSP*kpf1)/kyrho

     !Toroidal momentum pinch term
     DMST=-(eps/q)*dms +dreal(DMSP*kpf)/kyrho
     DMIT=DMIT+DMST

     !---- Summation of curvature pinch term fluxes here
     DMIP=DMIP+GP1+GP2     !!  Summation of curvature pinch terms

     v_p_thr = dmip

     !---- Quantities related to poloidal momentum transport

     !---- Obtained from Jan Weiland's momcode dated from Feb 4 2008
     !---- Not debugged or tested

     Vg=-((wu+ftrt)*dadk+0.5*(gti*gne/gni-EITH*gne)*dbdk) & 
          /(2.0*(wu+ftrt)+0.5*gne*adisp)
     RDexp=DREAL( 0.25*(gne*gti-eith*gne**2)*Vg/ (wu+ftrt)**2 )

     TSOUR=TSOUR+4.*RDexp*im_wu*im_wu*diffBohm*KXI*KPS/rmaj0**2
     BV=BV+im_wu*im_wu*DREAL(IU*NII*CONJG(NII-TICF))

     CN=wu+FTRT
     ACN=ABS(CN)

     ! End of looping through unstable modes
     ! -------
  End Do main_momentum

  !---- Final computation of momentum fluxes
  !---- First, some unknown computations

  CHIC=-SHPE*GK*2.*(rmaj/gti)
  HPT(1)=HPT(1)+EM*SHPE
  CHE(2)=CHE(2)+EM*SCHEF
  D(3)=D(3)+EM1*DEF
  !
  DMH=-2.D0*D1*TE**1.5*eps0*(Csound/CSOUND0)/vpol
  DMEF=DMH*DMI*LVF
  DMEF1=DMH*DMI1*LVF
  DMEF2=DMH*DMI2*LVF
  DMEF21=DMH*DMI21*LVF  !! splitting of DMEF2
  DMEF22=DMH*DMI22*LVF   !! splitting of DMEF2

  !----- Toroidal momentum transport

  !Diagonal element of momentum diffusivity
  DMDIAGT=D1*TE**1.5*DMDT/kyrho
  DMT(6)=DMDIAGT

  !Momentum pinches expressed as an effective diffusivity
  DMTEF = rmaj / rmaj0 * &
       ( DMIP  &
       + 2.0*d1*te**1.5*DMIT ) &
       /sign( max(abs(gvt*vtor),0.01) , gvt*vtor )

  !Total momentum pinch expressed as a convection velocity in m/s
  !Diagnostic only, not output in interface
  HPT(6) = ( dmip + 2.0*d1*te**1.5*dmit ) / sign( max( abs(vtor) ,0.01), vtor )

  !Parts of momentum pinches expressed as a convective velocity
  !These are output through the convective velocity array

  !Sum of thermoelectric and turbulent equipartition pinches
  v_p_thr = DMIP / sign( max(abs(rmaj0*vtor), 0.01), rmaj0*vtor )
!  v_p_thr = DMIP / sign( max(abs(vtor), 0.01), vtor )

  !Sum of parallel velocity and Reynolds stress pinches
  v_p_tor = 2.0*d1*te**1.5*DMIT / sign( max(abs(rmaj0*vtor), 0.01), rmaj0*vtor )

  !Total flux expressed as an effective diffusivity
  !WARNING: it is better to use diagonal element + convection
  !         Momentum pinches can produce large negative diffusion
  DMTEFF = DMT(6) + DMTEF

  !----- Poloidal momentum diffusivity
  ! WARNING: These computations have not been tested
  !          (There is no poloidal momentum eq. in PTRANSP,
  !           where the code was implemented and tested)

  !Diagonal element of poloidal momentum diffusivity
  !This should be fine
  DMDIAG=D1*TE**1.5*DMD/kyrho
  DM(4)=DMDIAG

  !Poloidal momentum velocity pinches
  !This looks like the Reynolds source terms
  HPT(4)=0.5D0*DMH*(DMI1+DMI21+DMI22)

  !Poloidal momentum pinches expressed as a diffusivity
  DMEFTEST=-2.*HPT(4)*LVF

  !Total flux expressed as a diffusivity 
  DMEFF=DM(4) -HPT(4)*2.*LVF

  !Unintelligible stuff. Quite likely, the Reynolds stress term
  !is being expressed as a source term here
  SP=1.
  IF(vpol.LT.0.) SP=-1.
  SMEF1=SP*D1*TE**1.5*DMI1/rmaj
  SMEF21=SP*D1*TE**1.5*DMI21/rmaj   !! splitting of DMEF2
  SMEF22=SP*D1*TE**1.5*DMI22/rmaj   !! splitting of DMEF2
!  Vconv(IK)=2.*(SMEF1+SMEF21+SMEF22)*(Csound/CSOUND0)
!  SRE=Vconv(IK)/(rmin*LVF)

  !End of momentum transport computations
  !-----------------------------------------------------------
  !Write effective diffusion and convection to returned arrays
  !

  diffs(1) = SCHI   !Ion heat diffusivity
  diffs(2) = SD     !Electron particle diffusivity?
  diffs(3) = SCHE   !Electron heat diffusivity
  diffs(4) = SDQ    !Impurity particle diffusivity
  diffs(5) = DMT(6) !Toroidal momentum diffusivity
  diffs(6) = DM(4)  !Poloidal momentum diffusivity
                    !Note that with the curvature effect omitted
                    !DM(4) is also the toroidal momentum diffusivity
  conv(1) = 0.0 
  conv(2) = 0.0 
  conv(3) = 0.0 
  conv(4) = v_p_tor ! Sum of parellel velocity + Reynolds stress pinch
  conv(5) = v_p_thr ! Thermoelectric pinch
  conv(6) = hpt(4) ! Sum of poloidal momentum related convection

  return
end subroutine w20diff!}}}

subroutine w20disp &!{{{
   ( zamr, zami, zbmr, zbmi, wz, G_ave, kps )
  !use w20data
  implicit none

  !Equation matrices
  real*8, dimension(neq,neq) :: zamr, zami, zbmr, zbmi

  !Input frequency
  complex*16, intent(in) :: wz

  REAL*8 RAV, G_ave
  REAL*8 H1,XH,R,HQR,HQI,WM

  COMPLEX*16 ALPC,ALPHA,HQ,IU,H2
  COMPLEX*16 A_lpc, A_lpk

  REAL*8 wr,wi,H

  real*8 :: alp, alf, kps
!  REAL*8 alp,kps,alf

  real*8 zdenom

  integer :: j1, j2

!
!    variables i=1,6: e Phi/Te, dTi/Ti, dni/ni, dTe/Te, dnq/nq, dTq/Tq
!    variables j=1,6 same as for i
!

  IU=(0.,1.)
!
!..print header
!
  if ( lprint .gt. 2 ) then
     write(hprint,*)
     write(hprint,*) &
          'Weiland-Nordman eigenvalue equations, subroutine etawn6'
     write(hprint,*) '(all frequencies normalized by omega_{De})'
     write(hprint,*) '(all diffusivities normalized by ' &
          ,'omega_{De} / k_y^2'
  endif


!     zflz   = aimp * zflh / zimp**2
!
!..diagnostic output
!

  if ( lprint .gt. 6 ) then
     write(hprint,*)
     write(hprint,*) '--------------------------------------'
     write(hprint,*)
     write(hprint,*)
     write(hprint,*) gne,' = gne'
     write(hprint,*) gni,' = gni'
     write(hprint,*) gnz,' = gnz'
     write(hprint,*) gte,' = gte'
     write(hprint,*) gti,' = gti'
     write(hprint,*) gtz,' = gtz'
     write(hprint,*)
     write(hprint,*) vef,' = vef'
     write(hprint,*) betae, '=betae'
     write(hprint,*) zimp,' = zimp'
     write(hprint,*) aimp,' = aimp'
     write(hprint,*) 
     write(hprint,*) zflh,' = zflh'
     write(hprint,*) zflz,' = zflz'
     write(hprint,*)
     write(hprint,*) bt,' =bt'
     write(hprint,*) bt1,'  =bt1'
     write(hprint,*) kappa,'  =kappa'
     write(hprint,*) wz,'   =wz'
  endif

  !Reset system of equations to zero
  zamr(1:neq,1:neq) = 0.0
  zami(1:neq,1:neq) = 0.0
  zbmr(1:neq,1:neq) = 0.0
  zbmi(1:neq,1:neq) = 0.0

!
!..Nine  equations with impurities, trapped electrons, parallel ion
!  motion, collisions,  FLR , finite beta and parallel motion of impurities
!
!  Equations for e phi/Te, Ti, ni, Tet, n_z, T_z, F, Av, K
!  Here F = omega*gamma (collisions), K = omega*Av and Av is the parallel
!   magnetic vector potential.
!

      H1=4.*q*q*zflh

      H=0.5*ABS(shat)/q
      H2=IU*H

      A_lpk=0.5*shat*SQRT(H1*XT*zflh*(1.0+TAUH*(gni+gti)/(2*WZ)) )
!      A_lpk=0.5*shat*SQRT(H1*XT*zflh*(1.0+TAUH*(gni+gti))/(2*WZ) )

      IF ( DREAL(A_lpk) .LT. 0) then
         A_lpk=-A_lpk
      END IF

      if ( abs(dreal(A_lpk)).lt.0.01 ) then
         A_lpk = A_lpk - dreal(A_lpk) + 0.01
      end if

      ALPC=-IU*A_lpk
      ALPHA=-IU*ABS(shat)*q*zflh
      XH=ABS(ALPHA/ALPC)
      ALPC=XH*ALPC
      R=2.*ABS(DREAL(WZ*ALPC))
      IF(R.LT.0.001) R=0.001    !! NEW 01.03.8
      HQ=2.*ALPC/H1

      !Geometric average multiplying \omega_{De}
      G_ave=(1.0+0.5*shear/R)*EXP(-0.25/R)-0.5*alpha_MHD*(1.-EXP(-1/R))
      G_ave=max(G_ave,1.0E-2)
!      G_ave = max(G_ave,1.0E-4)
! Temporary investigation
!      G_ave = 1.0

      alp=max(0.1,0.5*R)

      zdenom = (2.D0*zflh*q*q*betae*(1.D0 - fte))
      alf = (alp+1.0E-6) / (zdenom+1.0E-6)
      kps = 0.5*SQRT(alp/zflh)/q
      kpc = 1.D0
      RAV=1.+0.25*shat**2.0/alp

!
!  *********
      HQR=DREAL(HQ)
      HQI=DIMAG(HQ)
!  *********


      if ( lprint > 6 ) then
         write(hprint,10002) alp,shat**2,RAV
10002    FORMAT(2X,'alp=',ES15.6,'  SH2=',ES15.6,' RAV=',ES15.6)
         write(hprint,10003) XH,G_ave,alf
10003    FORMAT(2X,'XH=',ES15.6,' G_ave=',ES15.6,' alf=',ES15.6)
         write(hprint,10005) alpha_MHD
10005    FORMAT(2X,'MHD Alpha=',ES15.6)
         write(hprint,10006) DREAL(A_lpk),R
10006    FORMAT(2X,' A_lpk=',ES15.6,' R='ES15.6)
         write(hprint,*)"HQR=",hqr,"   HQI=",hqi
      end if

!
!--- WE NOW DEFINE THE MATRIX ELEMENTS --------------
! hydrogen density
!
!
!...Equations rewritten in terms of the profile gradients
!
        zamr(1,1) = - G_ave+HQR + 0.5 * ( gni - zflh * ztauh * ( gni + gti ) )
        zami(1,1) = HQI
        zamr(1,2) = (HQR-G_ave)*ztauh
        zami(1,2) = ztauh*HQI
        zamr(1,3) = (HQR-G_ave)*ztauh
        zami(1,3) = ztauh*HQI
        zamr(1,8) = -em*ztauh*HQR*0.5*(gni+gti)/kpc
        zami(1,8) = -em*ztauh*HQI*0.5*(gni+gti)/kpc
        zamr(1,9) = -em*HQR/kpc
        zami(1,9) = -em*HQI/kpc
!
        zbmr(1,1) = zflh
        zbmr(1,3) = 1.
!
!  hydrogen energy
!
        zamr(2,1) = 0.5*(gti - tvr * gni)
        zamr(2,2) = - ztauh * ftr
!
        zbmr(2,2) = 1.
        zbmr(2,3) = - tvr
!
!  total electron density expressed in ion density and imp density
!
        zamr(3,1) = -1.0 + 0.5*fte*gne
       zami(3,1) = vef*(1.-fte)
       zamr(3,3) = 1. - znz_ne -zns_ne
       zami(3,3) = -vef*(1. - znz_ne - zns_ne)
       zamr(3,4) = fte
       zamr(3,5) = znz_ne
       zami(3,5) = -vef*znz_ne
       zami(3,7) = vef*fte
       zamr(3,8) = -em*0.5*gne*(1. - fte)/kpc
       zami(3,8) =  em*0.5*gne*(1. - fte)*vef/kpc
       zamr(3,9) =  em* (1. - fte)* (1.+0.5*gne) / kpc
       zami(3,9) = -em* (1. - fte)* vef / kpc

      zbmr(3,1) = fte - 1.
      zbmr(3,3) = 1. - znz_ne - zns_ne
      zbmr(3,5) = znz_ne
      zbmr(3,9) = em*(1. - fte)/kpc
!
!  trapped electron energy
!
      zamr(4,1) = fte * 0.5 * ( gte - tvr*gne )
      zami(4,1) = vef*tvr*(bt-2.5*(1.-fte))
      zami(4,3) = -vef*tvr*bt1*(1.-znz_ne -zns_ne)
      zamr(4,4) = fte*ftr
      zami(4,5) = -vef*tvr*bt1*znz_ne
      zami(4,7) = -ftr*vef*fte
!
      zbmr(4,1) = (1. - fte)*tvr
      zbmr(4,3) = -(1. - znz_ne -zns_ne)*tvr
      zbmr(4,4) = fte
      zbmr(4,5) = -znz_ne*tvr
!
!  impurity density
!
      zamr(5,1) = - G_ave + zimp*HQR/aimp        &
           +0.5*( gnz - zflz*ztauz*(gnz+gtz) )
      zami(5,1) = zimp*HQI/aimp

      zamr(5,5) = (HQR*zimp/aimp-G_ave)*ztauz
      zami(5,5) = zimp*ztauz*HQI/aimp

      zamr(5,6) = (HQR*zimp/aimp-G_ave)*ztauz
      zami(5,6) = zimp*ztauz*HQI/aimp

      zamr(5,8) = -em*HQR*zimp*ztauz*0.5*(gnz+gtz)/(kpc*aimp)
      zami(5,8) = -em*HQI*zimp*ztauz*0.5*(gnz+gtz)/(kpc*aimp)

      zamr(5,9) = -em*HQR*zimp/(kpc*aimp)
      zami(5,9) = -em*HQI*zimp/(kpc*aimp)
 
      zbmr(5,1) = zflz
      zbmr(5,5) = 1.
!
!  impurity energy
!
      zamr(6,1) = 0.5*(gtz - tvr * gnz)
      zamr(6,6) = -ztauz*ftr
!
      zbmr(6,5) = -tvr
      zbmr(6,6) = 1.
!
!  variable F
!
      zamr(7,1) = 0.5*gte - 1.0
      zami(7,1) = vef
      ZAMR(7,7) = 1.
      zami(7,7) = -vef
!
      zbmr(7,1) = -1.
      zbmr(7,7) = 1.
!
!
!  electromagnetic parallel vectorpotential Av = e A_par/Te
!
      zamr(8,1) = em1*kpc*(0.5*gne+HQR*(fft+zimp*fzft/aimp))
      zami(8,1) = em1*HQI*(fft+zimp*fzft/aimp)*kpc

      zamr(8,2) = em1*HQR*ztauh*fft*kpc
      zami(8,2) = em1*HQI*ztauh*fft*kpc

      zamr(8,3) = em1*HQR*ztauh*fft*kpc
      zami(8,3) = em1*HQI*ztauh*fft*kpc

      zamr(8,5) = em1*HQR*zimp*ztauz*fzft*kpc/aimp
      zami(8,5) = em1*HQI*zimp*ztauz*fzft*kpc/aimp

      zamr(8,6) = em1*HQR*zimp*ztauz*fzft*kpc/aimp
      zami(8,6) = em1*HQI*zimp*ztauz*fzft*kpc/aimp

      zamr(8,8) = em1*(0.5*(gne + gte) - alf*zflh*RAV ) &
           - em1*HQR*(fft*ztauh*0.5*(gni+gti)           &
           + zimp*fzft*ztauz*0.5*(gtz+gnz)/aimp)


      zami(8,8) = -em1*HQI*(fft*ztauh*0.5*(gni + gti)   &
           +zimp*fzft*ztauz*0.5*(gnz + gtz)/aimp)


      zamr(8,9) = -em1*(0.5*gne+HQR*(fft+zimp*fzft/aimp))
      zami(8,9) = -em1*HQI*(fft+zimp*fzft/aimp)

      zbmr(8,1) = em1*kpc
      zbmr(8,8) = em1
      zbmr(8,9)= -em1
!
!     K = omega*Av
!
      zamr(9,9) = em1
!
      zbmr(9,8) = em1
!
      if ( lprint .gt. 6 ) then
        write(hprint,*)
        write(hprint,*) ' zamr(j1,j2)  j2 ->'
        do j1=1,9
          write(hprint,192) (zamr(j1,j2),j2=1,9)
        enddo

        write(hprint,*)

        write(hprint,*) ' zami(j1,j2)  j2 ->'
        do j1=1,9
          write(hprint,192) (zami(j1,j2),j2=1,9)
        enddo

       write(hprint,*)
        write(hprint,*) ' zbmr(j1,j2)  j2->'
        do j1=1,9
          write(hprint,192) (zbmr(j1,j2),j2=1,9)
        enddo

       write(hprint,*)
       write(hprint,*) ' zbmi(j1,j2)  j2->'
       do j1=1,9
         write(hprint,192) (zbmi(j1,j2),j2=1,9)
       enddo

 192  format (1p10e12.4)
      endif

    End subroutine w20disp!}}}

Subroutine w20solv &!{{{
   ( zz, zvr, zvi, iu, eu, imax, emax, G_ave, kps, ierr )
  !use w20data
  implicit none

  !Maximum number of iterations
  integer, parameter :: nitmax = 50

  !Eigenvectors resulting from solved system of equations
  real*8, intent(inout), dimension(neq,neq) :: zvr, zvi

  !Complex frequency of the mode 
  complex*16, intent(inout), dimension(neq) :: zz

  !Pointing indices (im-->unstable mode used, iu-->number of unstable modes)
  integer, intent(inout) :: iu, eu, imax, emax

  !External error checking flag
  integer, intent(inout) :: ierr

  !Geometric average needed to compute momentum pinches
  real*8 :: G_ave, G_ave_tmp, kps, kps_tmp

  !Internal error checking flag for eigenvalue solver
  integer :: ifail, kflag, it_err, it_conv
  
  !Internal counters for number of unstable modes and unstable mode used
  integer :: im, nu

  !Current number of iterations
  integer :: niter, niterm

  !Iteration indices
  integer :: i, j, nmod

  !Skip mode flag
  integer :: iskip, foundmode

  !Complex frequencies of fastest growing eigenmode
  complex*16, dimension(nitmax+1) :: wz, xbest, wzfnd

  !Input frequency for dispersion relation
  complex*16 :: wzin, dw1, dw

  !Temporary storage for growth rates and eigenmodes
  real*8                     :: wimax, wemax, wamax
  real*8, dimension(neq,neq) :: zvrtmp, zvitmp
  complex*16, dimension(neq) :: zztmp

!  complex*16, external :: cmplex

  !Mask array for active modes
!  integer, dimension(neq) :: mask_mod
  real*8,  dimension(neq) :: mask

  !Variables dealing with error control
  real*8 ::  zepsilon, ztol
  real*8 ::  delWZ, delZmin, delZ

  !Pointer to mode being used
  integer :: nmod_ptr

  !Initial count of unstable modes
  integer :: unst

  !Final test for largest growing mode in the other direction
  real :: womax, woomax, wiimax
  integer :: omax, ou, au, amax

  !Matrices and vectors for eigenvalue system
  real*8, dimension(neq) :: omega, gamma, zbeta, zalfi, zalfr, somega, sgamma
  real*8, dimension(neq,neq) :: zamr, zami, zbmr, zbmi
  
  real*8 :: zb

!  OPEN (UNIT=10,FILE='gamma.out')
!  OPEN (UNIT=11,FILE='gamma2.out')

  au = 0
  iu = 0
  eu = 0
  ou = 0
  unst = 0
  amax = 0
  imax = 0
  emax = 0
  omax = 0
  ierr = 0
  zepsilon = 1.0E-6
  ztol = 1.0E-3
  zz(1:neq)=(0.0,0.0)
  mask(1:neq)=0.0
!  mask_mod(1:neq)=0
  wamax = 0.0
  wimax = 0.0
  wemax = 0.0
  womax = 0.0
  wiimax = 0.0
  woomax = 0.0
  it_err = 0
  it_conv = 0
!  write(hprint,*) "Enter dispersion relation solver:"

  modes:  do nmod=1,neq

!     write(hprint,*) "Processing mode number",nmod
!     write(hprint,*) "======================"

     !Initialization of some parameters
     nmod_ptr = 0
     niter = 0
     ierr  = 0
     wz(1:nitmax+1)=(0.0,0.0)
     wzin=(0.0,0.0)
     wzfnd(1:nitmax+1)=(0.0,0.0)
     zztmp(1:neq)=(0.0,0.0)
     xbest(1:nitmax+1)=(0.0,0.0)
     delZ = 0.0
     delWZ = 0.0
     dw = (0.0,0.0)
     dw1 = (0.0,0.0)

     call w20wguess( wzin )

        if ( searchmode .eq. S_ELC .and. dreal(wzin).lt.0.0 ) then
           wzin=wzin+2*abs(dreal(wzin))
        end if

     wz(1)    = wzin
     xbest(1) = wzin     

     iskip = 0

     !Algorithm to find converged frequency begins here
     freq: Do
        !Initialize complex frequencies

        niter = niter + 1

        omega(1:neq)=0.0
        gamma(1:neq)=0.0


        !Maximum number of iterations 
        If ( niter .gt. nitmax ) Then
           it_err = 1
           If ( 0.1*abs(delWZ) .lt. ztol ) then
              exit freq
           !Hard error if tolerance not found within an order of magnitude
           Else If ( 0.1*abs(delWZ) .ge. ztol ) then
              ierr = 2
              it_err = 2
              Cycle modes
           End If
        End If

        !Obtain dispersion relation
        call w20disp(zamr, zami, zbmr, zbmi, wzin, G_ave_tmp, kps_tmp )

        !Diagonalize equation system
        call r8tomsqz(neq,neq,zamr,zami,zbmr,zbmi,zalfr,zalfi,zbeta,zvrtmp,zvitmp,ifail)

        !Check for errors
        If ( ifail .gt. 0 ) Then
           write(hprint,*) "In w20solv: eror in eigenvalue solver"
           ierr = 1
           exit modes
        End If

        zztmp = w20omg ( zalfr, zalfi, zbeta, neq )
        omega = dreal(zztmp)
        gamma = dimag(zztmp)

        !Begin carry out a series out tests to improve convergence and reliability

        !Test 1: if all modes are unstable, then return from this call
        unst = w20nunst ( zztmp, neq, S_ANY )

        if ( unst .eq. 0 .and. niter .ge. 3) then 
           it_conv = 1
           cycle modes
        end if

        !Test 2: lock onto the mode being looked at
        !        It should be nmod, but sometimes the solver rotates eigenvalues
        !        Thus find the mode that is most like the mode last used
        !        In the first iteration, choose the mode to be nmod
        if ( niter .eq. 1 ) then
           nmod_ptr = nmod
        else
           !For niter > 1, lock onto the mode last used
!           nmod_ptr = w20wmatch ( wzin, zztmp, neq )
           nmod_ptr = w20wmatch ( wzfnd(niter-1), zztmp, neq )
        end if

        !Save the matched mode for the next iteration
        wzfnd(niter)=zztmp(nmod_ptr)

        !Test 3: if the mode is stable, go to the next mode
        If ( gamma(nmod_ptr) .lt. 0.0 .and. niter .ge. 3 ) then
           it_conv = 1
           cycle modes
        end if

        foundmode = w20s_fnd ( zztmp(nmod_ptr) )

        !Test 4: if the mode is in the wrong direction, go to the next mode
        if ( foundmode .ne. searchmode .and. niter .gt. 3 ) then
           it_conv = 1
           cycle modes
        end if
        !Choose the next trial frequency as the average of the old and the new mode
        if ( niter .le. 1 ) then
           wz(niter+1)   = zztmp(nmod_ptr)
        else
           wz(niter+1)   = 0.5*(wzin+(0.D0,1.D0)*gamma(nmod_ptr)+omega(nmod_ptr))
        end if


        wzin          = wz(niter+1)
        Xbest(niter+1)= wzin

!        write(hprint,*) "Next mode is", wz(niter+1)

!        If ( niter .gt. 1 ) CALL AITKEN(wz,Xbest,niter+1,nitmax+1)       
        
        !Compute estimate of error
        delWZ = w20delwz( Xbest(niter+1), Xbest(niter) )

!        write(hprint,*) "Computed error is", delWz

        !Exit inner loop (freq) when convergence is found
        if ( delWZ .lt. ztol .and. niter .ge. 2) then
!           write(hprint,*) "Mode converged, modes found are:"
!           do j=1,niter
!              write(hprint,*) j, wzfnd(j)
!           end do
           it_conv = 1
           Exit freq
        end if

     End Do freq

     !If the unstable mode is in the direction we are looking for
     ! then evaluate whether this mode is faster growing than the 
     ! previously stored mode. If so, carry out the replacement

     foundmode = w20s_fnd ( zztmp(nmod_ptr) )
     
     if ( foundmode .eq. searchmode ) then
        if ( w20gamma( zztmp, neq, nmod_ptr ) .gt. wamax ) then       
           wamax = w20gamma( zztmp, neq, nmod_ptr )
           amax  = nmod_ptr
           G_ave = G_ave_tmp
           kps   = kps_tmp

           do j=1,neq
              foundmode = w20s_fnd( zztmp(j) )
              if ( j .ne. amax .and. foundmode .eq. searchmode .and. w20gamma( zztmp, neq, j ) .gt. wamax ) then
!                 zztmp(j)=zztmp(j)-(0.0,1.0)*zztmp(j)
              end if
           end do

           imax  = w20wmunst ( zztmp, neq, S_ION ) !nmod_ptr
           wimax = w20gamma ( zztmp, neq, imax ) !wamax
           emax  = w20wmunst ( zztmp, neq, S_ELC )
           wemax = w20gamma ( zztmp, neq, emax )

           eu  = w20nunst( zztmp, neq, S_ELC )
           iu  = w20nunst( zztmp, neq, S_ION )

           zz  = zztmp
           zvr = zvrtmp
           zvi = zvitmp
        end if
              
     end if

  End Do modes

  if ( it_conv .ne. 1 ) then
     write(hprint,*) "In w20solv: Excessive iterations in Weiland model"    

     if ( it_err .gt. 1 ) then
!        write(hprint,*)"Mode not converged, modes found are:"

!        do j=1,niter
!           write(hprint,*) j, wzfnd(j)
!        end do

        write(hprint,*) "mode=",nmod,"  ptr=",nmod_ptr
        write(hprint,*) "delWZ=",delWZ
        write(hprint,*) "GRDNE=",gne
        write(hprint,*) "GRDTE=",gte
        write(hprint,*) "GRDNI=",gni
        write(hprint,*) "GRDTI=",gti
        write(hprint,*) "GRDNZ=",gnz
        write(hprint,*) "GRDTZ=",gtz
        write(hprint,*) "NE=",ne
        write(hprint,*) "TE=",te
        write(hprint,*) "NH=",nh
        write(hprint,*) "TI=",th
        write(hprint,*) "NZ=",nz
        write(hprint,*) "NS=",ns_ne*ne
        write(hprint,*) "Q=",q
        write(hprint,*) "Shear=",shear
        write(hprint,*) "Elong=",kappa
        write(hprint,*) "Aimp=",aimp
        write(hprint,*) "Zimp=",zimp
        write(hprint,*) "Ahyd=",ahyd
        write(hprint,*) "Rmaj=",rmaj
        write(hprint,*) "Rmin=",rmin
        write(hprint,*) "btor=",btor
        !              write(hprint)') (wz(j),wzfnd(j),w20delwz(wz(j),wzfnd(j) ),j=1,niter)
        !              stop
     end if
  end if

 ! if (searchmode .eq. S_ION) then
 !    write(15,*)betae,gne,beta,wimax,dreal(zz(imax))
 !    write(hprint,*)"Vector found is:"
 !    do j=1,9
 !       write(hprint,*) zz(j)
 !    end do
 ! else if (searchmode .eq. S_ELC) then
 !    write(hprint,*)betae,gne,beta,wimax,dreal(zz(emax))
 !    write(hprint,*)"Vector found is:"
 !    do j=1,9
 !       write(hprint,*) zz(j)
 !    end do
 ! end if

!  stop

End Subroutine w20solv!}}}

function w20wmunst &!{{{
  ( zz, nmod, S_SRCH ) result ( nmod_ptr )
!This function returns the position of the mode with
!the largest growthrate
!
!INPUT: zz(nmod): eigenvalues of drift mode dispersion relation )
!       nmode   : dimension of the vector
!       S_SRCH  : direction S_ION or S_ELC or S_ANY
!
!OUTPUT: nmod_ptr : position of the largest growing mode in
!                   direction S_SRCH in array zz
!                   will return 0 if no 
implicit none

integer, intent(in) :: nmod
complex(ckind), intent(in), dimension(nmod) :: zz
integer, intent(in) :: S_SRCH

integer :: nmod_ptr, i, S_FND 
real(rkind) :: gamma_max

gamma_max = -10.0e3_rkind
S_FND = 0
nmod_ptr = N_NONE

do i=1,nmod

   if ( S_SRCH .ne. S_ANY ) then
      S_FND = w20s_fnd( zz(i) )
   end if
   
   if ( aimag(zz(i)) > gamma_max .and. S_FND .eq. S_SRCH ) then
      gamma_max = aimag(zz(i))
      nmod_ptr = i
   end if
end do

end function w20wmunst!}}}

    function w20s_fnd &!{{{
      ( zz ) result ( S_FND )
    !This function takes in a complex frequency and
    !evaluates whether the mode is an electron or ion mode
    !
    ! INPUT: zz = complex frequency ( w_r, gamma )
    ! OUTPUT: S_FND = S_ELC or S_ION ( integer flag for electron or ion mode)
    !
      implicit none

      complex(ckind), intent(in) :: zz

      integer :: S_FND

      if ( real(zz,kind=rkind) > 0.0 ) then
         S_FND = S_ELC
      else
         S_FND = S_ION
      end if

    end function w20s_fnd!}}}

function w20nunst &!{{{
( zz, nmod, S_SRCH ) result ( iunst )
!This function takes in a set of complex frequencies
!and counts the number of unstable modes in any direction
!
! INPUT: zz(nmod) = complex frequency
!        nmod     = dimension of zz vector
!        S_SRCH   = direction of search S_ANY, S_ION, S_ELC
!
! OUTPUT: iunst   = number of unstable modes in S_SRCH direction

implicit none

integer, intent(in) :: nmod, S_SRCH
complex(ckind), intent(in), dimension(nmod) :: zz

integer :: iunst, i, S_FND

iunst = N_NONE
S_FND = S_ANY

do i=1,nmod
   
   if ( S_SRCH .ne. S_ANY ) then
      S_FND = w20s_fnd( zz(i) )
   end if
   
   if ( aimag(zz(i)) > 0.0 .and. S_SRCH .eq. S_FND ) then
      iunst = iunst + 1
   end if
   
end do

end function w20nunst!}}}

function w20gamma &!{{{
  ( zz, nmod, nmod_ptr ) result ( gamma )
!This function returns the growthrate of the drift mode
!at the position nmod_ptr. If nmod_ptr is not between
!1 and nmod, then the result is zero. This is useful
!when no unstable modes are found in one of the functions
!above
!
!INPUT: zz(nmod) = complex eigenfrequencies
!       nmode    = dimension of zz vector
!       nmod_ptr = position of the growthrate
!
implicit none

integer, intent(in) :: nmod, nmod_ptr
complex(ckind), intent(in), dimension(nmod) :: zz

real(rkind) :: gamma
integer :: i

if ( nmod_ptr .eq. N_NONE ) then
   gamma = 0.0
else if ( nmod_ptr .ge. 1 .and. nmod_ptr .le. nmod ) then
   gamma = aimag(zz(nmod_ptr))
else
   write (*,*) "Invalid mode number: ",nmod_ptr
end if

end function w20gamma!}}}

function w20omg &!{{{
( zomega, zgamma, zbeta, nmod ) result ( zz )
!This function takes in the real and imaginary parts of the
!eigenvalues and folds them into a complex frequency
!The frequencies are then normalized

implicit none

integer, intent(in) :: nmod
real(rkind) , intent(in), dimension(nmod) :: zomega, zgamma, zbeta

complex(ckind), dimension(nmod) :: zz
real, dimension(nmod) :: zb

zz = (0.0,0.0)
zb = max(1.0E-4_rkind,zbeta)

zz(1:nmod) = ( zomega(1:nmod)+(0.0,1.0e0_rkind)*zgamma(1:nmod) ) / zb(1:nmod)

end function w20omg!}}}

function w20delwz &!{{{
( wz1, wz2 ) result ( delwz )
!This function computes the relative difference between two modes

implicit none

complex(ckind) :: wz1, wz2

real(rkind) :: delwz,denom
real(rkind) :: eps = 0.001e0_rkind
real(rkind) :: wr1, wr2, wi1, wi2

eps = 0.001e0_rkind
wr1 = real(wz1,kind=rkind)
wr1 = sign( max( abs(wr1), eps ), wr1)
wi1 = aimag(wz1)
wi1 = sign( max( abs(wi1), eps ), wi1)

wr2 = real(wz2,kind=rkind)
wi2 = aimag(wz2)

denom = ( max( abs(wz1), eps ) )**2
denom = 0.5e0_rkind * ( abs(wz1) + abs(wz2) )

delwz = sqrt( ( wr1-wr2 )**2 + ( wi1-wi2 )**2 ) / denom
!      delwz = ( (1.0e0_rkind-wr2/wr1)**2 + (1.0e0_rkind-wi2/wi1)**2 ) ! / denom 



end function w20delwz!}}}

function w20wmatch &!{{{
( wz, zz, nmod ) result ( nmod_ptr )
!This function returns the mode closest to the input mode
implicit none

integer, intent(in) :: nmod
complex(ckind), intent(in) :: wz
complex(ckind), intent(in), dimension(nmod):: zz

integer :: nmod_ptr,i
real :: delwz, delwzmin

nmod_ptr = N_NONE
delwzmin = 10.0e6_rkind
delwz = 0.0

do i=1,nmod
   delwz = w20delwz( wz, zz(i) )
   if (  delwz .lt. delwzmin ) then
      delwzmin = delwz
      nmod_ptr = i
   end if
end do

end function w20wmatch!}}}

subroutine w20wguess &!{{{
( wz )
!use w20data
implicit none

COMPLEX*16, intent(inout) :: wz

COMPLEX*16 :: wz1, wz2, H_q, E1, E2
COMPLEX*16 :: iu
REAL*8     :: G_ave

G_ave = 1.0
IU = (0.0, 1.0)

H_q = cmplx(0.0,0.5*shat/q)

E1 = ftrt*(1.0+zflh)-0.5*gne+0.5*zflh*tauh*(gne+gti) &
  + (gm+ftrt)*G_ave+H_q*(1.0+ftr*tauh)

E1 = E1 * 0.5 / ( 1.0 + zflh )


E2 =  (  ( 0.5 * tauh * gm * ( gti - tvr*gne ) + bta )*( G_ave + H_q ) &
    - 0.5 * tauh * ftr * ( gni - zflh * tauh * (gti + gni))  ) &
    / ( 1.0 + zflh )

wz1=-E1+SQRT(E1*E1-E2)
wz2=-E1-SQRT(E1*E1-E2)

wz=wz1

IF(DIMAG(wz2).GT.DIMAG(wz1)) wz=wz2

if ( dimag(wz) .le. 0.01 ) then
  wz = wz - iu*dimag(wz) + iu*0.01
end if

end subroutine w20wguess!}}}

end module w20mod
