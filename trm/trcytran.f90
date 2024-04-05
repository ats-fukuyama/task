MODULE tr_cytran_mod

!-------------------------------------
! Interface of TASK/TR with CYTRAN_MOD
!-------------------------------------

CONTAINS

SUBROUTINE tr_cytran

USE bpsd_kinds
use cytran_mod

USE trcomm, ONLY: &
     nrmax,rn,rt,abb1rho,pvolrhog,psurrhog,prc,syncabs,syncself,aee

!------------------------------------------------------------------------------
!Declaration of local variables

REAL(KIND=rkind) :: &
  den_rm(nrmax),       & !electron density in cell [/m**3]
  te_rm(nrmax),        & !electron temperature in cell [keV]
  bavg_rm(nrmax),      & !<B> in cell [T]
  dvol_rm(nrmax),      & !volume of cell [m**3]
  area_rg(nrmax+1),    & !surface area of inner boundary of cell [m**2]
  psync_rm(nrmax)        !net power source/loss in cell i [keV/m**3/s]
REAL(KIND=rkind) :: &
  fabs,                & !fraction of cyclotron radiation absorption by walls
  fself                  !fraction of x (o) mode reflected as x (o) mode

INTEGER :: &
  k_cyt_res              !option for resolution of frequency interval
                         !=1    gives default of 100 intervals
                         !>1    gives enhancement 100*k_cyt_res
                         !      10 recommended to keep graininess < a couple %
INTEGER :: &
     nr
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

do nr=1,nrmax
   den_rm(nr)=rn(nr,1)*1.D20
   te_rm(nr)=rt(nr,1)
   bavg_rm(nr)=abb1rho(nr)
   area_rg(nr+1)=psurrhog(nr)
enddo
area_rg(1)=0.D0
dvol_rm(1)=pvolrhog(1)
do nr=2,nrmax
   dvol_rm(nr)=pvolrhog(nr)-pvolrhog(nr-1)
enddo
fabs=syncabs
fself=syncself

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

CALL CYTRAN(ree,reo,roo,roe,nrmax,bavg_rm,den_rm,te_rm,area_rg, &
            dvol_rm,psync_rm, &
            K_CYT_RES=k_cyt_res)

!Set outside ghost point value same as last node inside plasma
do nr=1,nrmax
   prc(nr)=-psync_rm(nr)*aee*1.D3
end do

!do nr=1,nrmax
!   write(6,'(A,I5,1P6E12.4)') '@ ',nr,bavg_rm(nr),den_rm(nr),te_rm(nr), &
!        area_rg(nr),dvol_rm(nr),prc(nr)
!end do
!write(6,'(A,1P4E12.4)') '@ ',ree,reo,roo,roe
!------------------------------------------------------------------------------
!Cleanup and exit
!------------------------------------------------------------------------------
RETURN

END SUBROUTINE tr_cytran

END MODULE tr_cytran_mod
      
