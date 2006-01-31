MODULE physical_constants
  IMPLICIT NONE

  !------------------------!
  ! Set physical constants !
  !------------------------!

  !   Electron charge (C)
  REAL(8), PARAMETER :: AEE  = 1.60217733D-19

  !   Electron mass (kg)
  REAL(8), PARAMETER :: AME  = 9.1093897D-31

  !   Proton mass (kg)
  REAL(8), PARAMETER :: AMP  = 1.6726231D-27

  !   Light velocity (m/s)
  REAL(8), PARAMETER :: VC   = 2.99792458D8

  !   Pi (circle ratio)
  REAL(8), PARAMETER :: PI   = 3.14159265358979323846D0

  !   mu0 (H/m)
  REAL(8), PARAMETER :: rMU0 = 4.D0 * PI * 1.D-7
  
  !   epsilon0 (F/m)
  REAL(8), PARAMETER :: EPS0 = 1.D0 / (rMU0 * VC**2)

  !   Conversion factor from keV to joule
  REAL(8), PARAMETER :: rKEV = 1.D3 * AEE
  
END MODULE physical_constants

