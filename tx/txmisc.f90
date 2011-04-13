!!! Miscellaneous libraries related to the physics or physical model

!***************************************************************
!
!   Coefficient function of CDBM model
!
!***************************************************************

pure REAL(8) FUNCTION TRCOFS(S,ALFA,RKCV)

  implicit none
  real(8), intent(in) :: S, ALFA, RKCV
  real(8) :: SA, FS1, FS2

  IF(ALFA > 0.D0) THEN
     SA = S - ALFA
     IF(SA > 0.D0) THEN
        FS1 = (1.D0 + 9.0D0 * SQRT(2.D0) * SA**2.5D0) &
             &  / (SQRT(2.D0) * (1.D0 - 2.D0 * SA + 3.D0 * SA**2 + 2.0D0 * SA**3))
     ELSE
        FS1 = 1.D0 / SQRT(2.D0 * (1.D0 - 2.D0 * SA) * (1.D0 - 2.D0 * SA + 3.D0 * SA**2))
     ENDIF
     IF(RKCV > 0.D0) THEN
        FS2 = SQRT(RKCV)**3 / S**2
     ELSE
        FS2 = 0.D0
     ENDIF
  ELSE
     SA = ALFA - S
     IF(SA > 0.D0) THEN
        FS1 = (1.D0 + 9.0D0 * SQRT(2.D0) * SA**2.5D0) &
             &  / (SQRT(2.D0) * (1.D0 - 2.D0 * SA + 3.D0 * SA**2 + 2.0D0 * SA**3))
     ELSE
        FS1 = 1.D0 / SQRT(2.D0 * (1.D0 - 2.D0 * SA) * (1.D0 - 2.D0 * SA + 3.D0 * SA**2))
     ENDIF
     IF(RKCV < 0.D0) THEN
        FS2 = SQRT(-RKCV)**3 / S**2
     ELSE
        FS2 = 0.D0
     ENDIF
  ENDIF
  TRCOFS = MAX(FS1,FS2)

END FUNCTION TRCOFS

!***************************************************************
!
!   Coulomb logarithm
!
!     imodel = 1 : electron - electron
!              2 : electron - ion
!              3 : ion - ion
!
!***************************************************************

pure real(8) function coulog(imodel, Ne, Te, Ni, Ti, PA, PZ) result(f)

  implicit none
  integer(4), intent(in) :: imodel
  real(8), intent(in) :: Ne, Te
  real(8), intent(in), optional :: Ni, Ti, PA, PZ
  integer(4), parameter :: memp = 5.446170221d-4 ! me/mp
  real(8) :: Ne_m3, Ni_m3, Te_eV, Ti_eV, rat_mass

  !     (NRL Plasma Formulary p34,35 (2007))
  Ne_m3 = Ne * 1.d20 ; Te_eV = Te * 1.d3

  if (imodel == 1) then
     f = 30.4d0 - log(sqrt(Ne_m3)/(Te_eV**1.25d0)) &
          &     - sqrt(1.d-5+(log(Te_eV)-2.d0)**2/16.d0)
  else
     if(present(Ni) .and. present(Ti) .and. present(PA) .and. present(PZ)) then
        Ni_m3 = Ni * 1.d20 ; Ti_eV = Ti * 1.d3
        if (imodel == 2) then
           rat_mass = memp / PA
           if(10.d0*PZ**2 > Ti_eV*rat_mass .and. 10.d0*PZ**2 < Te_eV) then ! usual case
              f = 30.9d0 - log(sqrt(Ne_m3)/Te_eV)
           else if(Te_eV > Ti_eV*rat_mass .and. Te_eV < 10.d0*PZ**2) then ! rare case
              f = 29.9d0 - log(sqrt(Ne_m3)*PZ/Te_eV**1.5d0)
           else if(Te_eV < Ti_eV*PZ*rat_mass) then ! very rare case
              f = 36.9d0 - log(sqrt(Ni_m3)/Ti_eV**1.5d0*PZ**2/PA)
           else ! assumption
              f = 30.9d0 - log(sqrt(Ne_m3)/Te_eV)
           end if
!tokamaks       f = 37.8d0 - LOG(SQRT(Ne_m3)/(Te))
        else ! imodel = 3
           f = 29.9d0 - log(PZ**2/Ti_eV*sqrt(2.d0*Ni_m3*PZ**2/Ti_eV))
!tokamaks       f = 40.3d0 - LOG(PZ**2/Ti*SQRT(2.D0*Ni_m3*PZ**2/Ti))
        end if
     else
        
     end if
  end if

end function coulog

!***************************************************************
!
!   Correction factor for resistivity
!     (Hirshman and Sigmar, (1981), Eq. (7.36))
!
!***************************************************************

pure REAL(8) FUNCTION CORR(X)
  ! X is the effective charge number
  real(8), intent(in) :: X

  CORR = (1.D0 + (1.198D0 + 0.222D0 * X) * X) * X &
  &    / (1.D0 + (2.966D0 + 0.753D0 * X) * X)

END FUNCTION CORR
