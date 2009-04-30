MODULE tr_cytran_mod

!-------------------------------------
! Interface of TASK/TR with CYTRAN_MOD
!-------------------------------------

USE bpsd_kinds
use cytran_mod
IMPLICIT NONE

CONTAINS

SUBROUTINE tr_cytran(nrmax,den_rm,te_rm,bavg_rm,dvol_rm,area_rg, &
                     fabs,fself,psync_rm)

!------------------------------------------------------------------------------
!Declaration of interface variables

INTEGER, INTENT(IN) :: &
  nrmax                  !no. of radial plasma nodes [-]
REAL(KIND=rkind), INTENT(IN) :: &
  den_rm(nrmax),       & !electron density in cell [/m::3]
  te_rm(nrmax),        & !electron temperature in cell [keV]
  bavg_rm(nrmax),      & !<B> in cell [T]
  dvol_rm(nrmax),      & !volume of cell [m**3]
  area_rg(nrmax)         !surface area of inner boundary of cell [m**2]
REAL(KIND=rkind), INTENT(IN) :: &
  fabs,                & !fraction of cyclotron radiation by walls [-]
  fself                  !fraction of x (o) mode reflected as x (o) mode [-]
REAL(KIND=rkind), INTENT(OUT) :: &
  psync_rm(nrmax)

!Declaration of internal variables

INTEGER :: &
  k_cyt_res              !option for resolution of frequency interval
                         !=1    gives default of 100 intervals
                         !>1    gives enhancement 100*k_cyt_res
                         !      10 recommended to keep graininess < a couple %
REAL(KIND=rkind) :: &
  ree,                 & !refl of X mode from incident X mode [0-1]
  reo,                 & !refl of X mode from incident O mode [0-1]
  roe,                 & !refl of O mode from incident X mode [0-1]
  roo                    !refl of O mode from incident O mode [0-1]
  
!Radial grid
!-----------------------------------------------------------------------------
!Schematic:
!  rho_g            |     |     |     | ... |     |
!    Label          1     2     3     4    n-1    n
!  rho_m               |     |     |           |     |
!    Label             1     2     3     4    n-1    n
!                   |                             |
!                   Magnetic axis (rho=0)         Plasma boundary (rho=1)

!------------------------------------------------------------------------------
!Set reflection factors
!------------------------------------------------------------------------------
!ree-fraction fself of the wall reflection coefficient assigned to the
!    reflection of x mode relative to incident x mode intensity
ree=(1.0-fabs)*fself

!reo-fraction 1-fself of the wall reflection coefficient assigned to the
!    reflection of x mode relative to incident o mode intensity
reo=(1.0-fabs)*(1.0-fself)

!roo-fraction fself of the wall reflection coefficient assigned to the
!    reflection of o mode relative to incident o mode intensity
roo=(1.0-fabs)*fself

!roe-fraction 1-fself of the wall reflection coefficient assigned to the
!    reflection of o mode relative to incident x mode intensity
roe=(1.0-fabs)*(1.0-fself)

!------------------------------------------------------------------------------
!Get the synchrotron radiation power profile
!------------------------------------------------------------------------------
k_cyt_res=10

CALL CYTRAN(ree,reo,roo,roe,nrmax-1,bavg_rm,den_rm,te_rm,area_rg, &
            dvol_rm,psync_rm, &
            K_CYT_RES=k_cyt_res)

!Set outside ghost point value same as last node inside plasma
psync_rm(nrmax)=psync_rm(nrmax-1)

!------------------------------------------------------------------------------
!Cleanup and exit
!------------------------------------------------------------------------------
RETURN

END SUBROUTINE tr_cytran

END MODULE tr_cytran_mod
      
