module physical_constants
  implicit none
  public

  !------------------------!
  ! Set physical constants !
  !------------------------!

  !   Electron charge (C)
  real(8), parameter :: AEE  = 1.60217733D-19

  !   Electron mass (kg)
  real(8), parameter :: AME  = 9.1093897D-31

  !   Proton mass (kg)
  real(8), parameter :: AMP  = 1.6726231D-27

  !   Light velocity (m/s)
  real(8), parameter :: VC   = 2.99792458D8

  !   Pi (circle ratio)
  real(8), parameter :: PI   = 3.14159265358979323846D0

  !   mu0 (H/m)
  real(8), parameter :: rMU0 = 4.D0 * PI * 1.D-7
  
  !   epsilon0 (F/m)
  real(8), parameter :: EPS0 = 1.D0 / (rMU0 * VC**2)

  !   Conversion factor from keV to joule
  real(8), parameter :: rKeV = 1.D3 * AEE
  
end module physical_constants
