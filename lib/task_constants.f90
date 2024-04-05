! task_constants.f90

module task_constants
  use task_kinds
  implicit none

  real(rkind),parameter :: ZERO = 0.0_dp
  real(rkind),parameter :: HALF = 0.5_dp
  real(rkind),parameter :: ONE  = 1.0_dp
  real(rkind),parameter :: TWO  = 2.0_dp
  real(rkind),parameter :: ONEHALF = 1.5_dp

  real(rkind),parameter :: PI   = 3.14159265358979323846_dp
  real(rkind),parameter :: TWOPI= PI+PI

  ! ***** Physical constants based on CODATA 2006 *****

!  real(rkind),parameter :: AEE  = 1.602176487E-19_dp ! elementary charge
!  real(rkind),parameter :: AME  = 9.10938215E-31_dp  ! electron mass
!  real(rkind),parameter :: AMP  = 1.672621637E-27_dp ! proton mass
!  real(rkind),parameter :: VC   = 2.99792458E8_dp    ! speed of light
!  real(rkind),parameter :: RMU0 = 4.E-7_dp*PI        ! permeability
!  real(rkind),parameter :: EPS0 = ONE/(VC*VC*RMU0)   ! permittivity

  ! ***** Physical constants based on CODATA 2010 *****

!  real(rkind),parameter :: AEE  = 1.602176565E-19_dp ! elementary charge [C]
!  real(rkind),parameter :: AME  = 9.10938291E-31_dp  ! electron mass     [kg]
!  real(rkind),parameter :: AMP  = 1.672621777E-27_dp ! proton mass       [kg]
!  real(rkind),parameter :: VC   = 2.99792458E8_dp    ! speed of light    [m/s]
!  real(rkind),parameter :: RMU0 = 4.E-7_dp*PI        ! permeability      [H/m]
!  real(rkind),parameter :: EPS0 = ONE/(VC*VC*RMU0)   ! permittivity

  ! ***** Physical constants based on CODATA 2018 *****

  real(rkind),parameter :: AEE  = 1.602176634E-19_dp  ! elementary charge [C]
  real(rkind),parameter :: AME  = 9.1093837015E-31_dp ! electron mass     [kg]
  real(rkind),parameter :: AMDa = 1.66053906660E-27_dp! Dalton (amu)      [kg]
  real(rkind),parameter :: AMP  = 1.67262192369E-27_dp! proton mass       [kg]
  real(rkind),parameter :: AMN  = 1.67492749804E-27_dp! neutron mass      [kg]
  real(rkind),parameter :: AMD  = 3.3435837724E-27_dp ! deuteron mass     [kg]
  real(rkind),parameter :: AMT  = 5.0073567446E-27_dp ! triton mass       [kg]
  real(rkind),parameter :: AMA  = 6.6446573357E-27_dp ! alpha mass        [kg]
  real(rkind),parameter :: AMH1 = 1.00782503207*AMDa  ! Hydrogen mass     [kg]
  real(rkind),parameter :: AMH2 = 2.0141017778*AMDa   ! deuterium mass    [kg]
  real(rkind),parameter :: AMH3 = 3.0160492777*AMDa   ! tritium mass      [kg]
  real(rkind),parameter :: AMH  = 1.00794*AMDa        ! Hydrogen mass     [kg]
  real(rkind),parameter :: AMHe3= 3.0160293191*AMDa   ! helium3 mass      [kg]
  real(rkind),parameter :: AMHe4= 4.00260325415*AMDa  ! helium4 mass      [kg]
  real(rkind),parameter :: AMHe = 4.002602*AMDa       ! helium mass       [kg]
  real(rkind),parameter :: VC   = 2.99792458E8_dp     ! speed of light    [m/s]
  real(rkind),parameter :: RMU0 = 4.E-7_dp*PI         ! permeability      [H/m]
  real(rkind),parameter :: EPS0 = ONE/(VC*VC*RMU0)    ! permittivity
  real(rkind),parameter :: HPLK = 6.62607015E-34_dp   ! Planck constant   [Js]
  real(rkind),parameter :: RKBL = 1.380649E-23_dp     ! Boltzmann constant[J/K]
  real(rkind),parameter :: RNAV = 6.02214076E23_dp    ! Avogadro constant[/mol]

  complex(rkind),parameter :: CI   = (0.0_dp,1.0_dp)
end module task_constants
