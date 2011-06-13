!*****************************************************************************
!   Coulomb logarithm Formulae given by Mitsuru HONDA (2011/06/13)
!     References: M. Honda, in preparation
!                 D.V. Sivukhin, Rev. Plasma Phys. Vol. 1 (1966)
!
! Note: Beam temperature defined by Tb = (2/3)Eb, where Eb: injection energy
!*****************************************************************************

!     Note: Coulog assumes that the ions share a common temperature, ti.
!           Set "2" to beam ions if "1" is the ion species.
function coulog( zeff, ne, te, ti, A1, Z1, A2, Z2, tb ) result( lambda )
  implicit none
    
  ! arguments
  real(8), intent(in) :: zeff, ne, te, ti, A1, Z1, A2, Z2
  real(8), intent(in), optional :: tb
  ! zeff   : effective charge
  ! ne     : electron density in m^{-3} or 10^20 m^{-3}
  ! te, ti : electron and ion temperatures in eV or keV
  ! A1, A2 : mass number
  ! Z1, Z2 : ABSOLUTE charge number, e.g. Ze = 1
  ! tb     : mean temperature of beam ions defined by Tb = (2/3)Eb in eV or keV
  real(8), parameter :: Ae = 5.446170219d-4 ! from CODATA 2010
  real(8), parameter :: eps_max = 1.d3 ! margin
  real(8) :: lambda, CD, t1, t2, fac_n, fac_t, eps
  real(8) :: coef, coef1, coef2, thres_l, thres_r

  eps = eps_max * epsilon(1.d0)

  if (ne > 1.d10) then
     ! density in m^{-3}, temperature in eV
     fac_n = 1.d0  ; fac_t = 1.d0
  else
     ! density in 10^20 m^{-3}, temperature in keV
     fac_n = 1.d20 ; fac_t = 1.d3
  end if

  CD = ne * ( 1.d0 / te + zeff / ti ) * fac_n / fac_t

  if ( present( tb ) ) then
     if ( ( A1 - Ae ) < eps ) then
        t1 = te ; t2 = tb
        if ( ( A2 - Ae ) < eps ) stop 'Error!'
     else if ( ( A2 - Ae ) < eps ) then
        t1 = tb ; t2 = te
        if ( ( A1 - Ae ) < eps ) stop 'Error!'
     else
        t1 = ti ; t2 = tb
     end if
  else
     t1 = te ; t2 = te
     if ( ( A1 - Ae ) > eps ) t1 = ti
     if ( ( A2 - Ae ) > eps ) t2 = ti
  end if

  coef1 = ( A1 + A2 ) / ( A2 * t1 + A1 * t2 ) / fac_t
  coef2 = A1 * A2 / ( A1 + A2 )
  coef  = log( coef1 * sqrt(CD) )

  thres_l = 1.d0 / coef1
  thres_r = 2.45d4 * ( Z1 * Z2 )**2 * coef2
  if ( thres_l <= thres_r ) then
     !     Classical formula
     lambda = 30.4d0 - log( Z1 * Z2 ) - coef
  else
     !     Quantum-mechanical formula
     lambda = 35.4d0 - coef + 0.5d0 * log( coef1 * coef2 )
  end if

  return
end function coulog

!

pure function coulog_gen( ne, te, CDi, A1, Z1, t1, A2, Z2, t2 ) result( lambda )

  ! arguments
  real(8), intent(in) :: ne, te, CDi, A1, Z1, t1, A2, Z2, t2
  ! ne     : electron density in m^{-3} or 10^20 m^{-3}
  ! te     : electron temperature in eV or keV
  ! CDi    : sum_j (Z_j^2 n_j/T_j)
  ! A1, A2 : mass number
  ! Z1, Z2 : ABSOLUTE charge number, e.g. Ze = 1
  ! t1, t2 : temperature in eV or keV
  real(8) :: lambda, CD, fac_n, fac_t
  real(8) :: coef, coef1, coef2, thres_l, thres_r

  if (ne > 1.d10) then
     ! density in m^{-3}, temperature in eV
     fac_n = 1.d0  ; fac_t = 1.d0
  else
     ! density in 10^20 m^{-3}, temperature in keV
     fac_n = 1.d20 ; fac_t = 1.d3
  end if

  CD = ( ne / te + CDi ) * fac_n / fac_t

  coef1 = ( A1 + A2 ) / ( A2 * t1 + A1 * t2 ) / fac_t
  coef2 = A1 * A2 / ( A1 + A2 )
  coef  = log( coef1 * sqrt(CD) )

  thres_l = 1.d0 / coef1
  thres_r = 2.45d4 * ( Z1 * Z2 )**2 * coef2
  if ( thres_l <= thres_r ) then
     !     Classical formula
     lambda = 30.4d0 - log( Z1 * Z2 ) - coef
  else
     !     Quantum-mechanical formula
     lambda = 35.4d0 - coef + 0.5d0 * log( coef1 * coef2 )
  end if

  return
end function coulog_gen

!***************************************************************
!
!   Coulomb logarithm
!     Reference: NRL Formulary (2009)
!
!     imodel = 1 : electron - electron
!              2 : electron - ion
!              3 : ion - ion
!
!***************************************************************

pure real(8) function coulog_NRL(imodel, Ne, Te, Ni, Ti, PA, PZ) result(f)

  implicit none
  integer(4), intent(in) :: imodel
  real(8), intent(in) :: Ne, Te
  real(8), intent(in), optional :: Ni, Ti, PA, PZ
  integer(4), parameter :: memp = 5.446170219d-4 ! me/mp
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

end function coulog_NRL
