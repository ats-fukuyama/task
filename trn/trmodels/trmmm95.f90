MODULE trmmm95

  PRIVATE
  PUBLIC tr_mmm95

CONTAINS

! **********************************************************************
!
!  The interface between TASK/TR(trcoeftb) and MMM95 model(mmm95)
!
! **********************************************************************

  SUBROUTINE tr_mmm95

    USE trcomm, ONLY: &
         rkind,ikind,nrmax,nsamax,neqmax,idnsa,ns_nsa,pa,pz, &
         rmnrho,rmjrho,rkprho,abb1rho,BB,RR,rn,qp,           &
         cdtrn,cdtru,cdtrt,dtr_tb,vtr_tb,                    &
         rt_e,rt_ecl,rt_i,rt_icl,rn_e,rn_ecl,rn_i,rn_icl,    &
         qp_d,ai_ave,mshear,wexbp,z_eff
         !nrd1,nrd2,nrd3

    IMPLICIT NONE

    ! Inputs array for mmm95 (subroutine): the dimension number is 'npoints'
    REAL(rkind),DIMENSION(1) :: &
         rminor,rmajor,elong,dense,densh,densimp,densfe,xzeff,      &
         tekev,tikev,q,btor,avezimp,amassimp,amasshyd,aimass,wexbs, &
         grdne,grdni,grdnh,grdnz,grdte,grdti,grdq

    ! Outputs array for mmm95 (subroutine): 
    ! **  The dimension number of radial grid is 'npoints'
    ! **  The dimension number of transport matrix is 'matdim = 5'
    REAL(rkind),DIMENSION(1) :: &
         thiig,thdig,theig,thzig,thirb,thdrb,therb,thzrb,  &
         thikb,thdkb,thekb,thzkb
    REAL(rkind),DIMENSION(1:5,1) :: &
         gamma,omega
    REAL(rkind),DIMENSION(1:5,1) :: &
         velthi,vflux
    REAL(rkind),DIMENSION(1:5,1:5,1) :: &
         difthi

    ! Inputs integer for mmm95 (subroutine)
    INTEGER(ikind) :: &
         npoints,nprout,matdim
    ! Inputs switches for mmm95 (subroutine)
    INTEGER(ikind) :: &
         lprint,lsuper,lreset
    ! Outputs integer for mmm95 (subroutine)
    INTEGER(ikind) :: nerr
    
    ! Outputs control variables for mmm95 (subroutine)
    ! +-- These variables are determined in the mmm95 subroutine.
    INTEGER(ikind),DIMENSION(4) :: fig,frb,fkb
    INTEGER(ikind),DIMENSION(5) :: lswitch
    INTEGER(ikind),DIMENSION(23) :: cswitch


    ! Internal varialbes
    INTEGER(ikind) :: nr,ns,nsa
    REAL(rkind) :: sum_pan,pa_ave
    REAL(rkind):: &
         rmnrhom,rmjrhom,rkprhom,rn_em,rn_im,rt_em,rt_im,z_effm, &
         qpm,abb1rhom
    REAL(rkind),DIMENSION(0:nrmax) :: &
       mmm_diff,mmm_chie,mmm_chii,mmm_vele,mmm_veli
    REAL(rkind),DIMENSION(1:nrmax) :: &
       mmm_diffm,mmm_chiem,mmm_chiim,mmm_velem,mmm_velim


    ! Initilization
    dtr_tb(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    vtr_tb(1:neqmax,1:neqmax,0:nrmax) = 0.d0

    !  << Input integers >>
    matdim  = 5  ! the dimension of transport matricies ( j1 and j2 )
                 !   matdim must be at least 5
    npoints = 1  ! number of values of jz in all of the arrays
    nprout  = 95 ! output unit number for long printout
                 !  : for error and diagnostic
    
    !  +++ Input switches +++
    lprint  = 0  ! controls the amount of printout (0 => no printout)
                 !   higher values yield more diagnostic output
    lsuper  = 0  ! = 0 for simulations of all other discharges
                 ! > 0 for supershot simulations; substantially reduces 
                 !       contribution from kinetic ballooning mode
    lreset  = 0  ! = 0 to use only internal settings for lswitch, cswitch

    lswitch = 0
    cswitch = 0
    fig     = 0
    frb     = 0
    fkb     = 0

    DO nr = 1, nrmax
       rmnrhom  = 0.5d0*(rmnrho(nr)  +  rmnrho(nr-1))
       rmjrhom  = 0.5d0*(rmjrho(nr)  +  rmjrho(nr-1))
       rkprhom  = 0.5d0*(rkprho(nr)  +  rkprho(nr-1))
       rn_em    = 0.5d0*(rn_e(nr)    +    rn_e(nr-1))
       rn_im    = 0.5d0*(rn_i(nr)    +    rn_i(nr-1))
       rt_em    = 0.5d0*(rt_e(nr)    +    rt_e(nr-1))
       rt_im    = 0.5d0*(rt_i(nr)    +    rt_i(nr-1))
       z_effm   = 0.5d0*(z_eff(nr)   +   z_eff(nr-1))
       qpm      = 0.5d0*(qp(nr)      +      qp(nr-1))
       abb1rhom = 0.5d0*(abb1rho(nr) + abb1rho(nr-1))

       rminor(1) = rmnrhom
       rmajor(1) = rmjrhom
       elong(1)  = rkprhom

       dense(1) = rn_em*1.d20 ! electron density [m^-3]
       densh(1) = rn_im*1.d20 ! sum over thermal hyd. ion densities [m^-3]

       ! --- for the time being ---
       densimp(1) = 0.d0! sum over impurity ion densities [m^-3]
       densfe(1)  = 0.d0! electron density from fast (non-thermal) ions [m^-3]

       xzeff(1) = z_effm
       tekev(1) = rt_em
       tikev(1) = rt_im
       q(1)     = qpm
!       btor(1)  = BB
       btor(1)  = abb1rhom

       ! --- for the time being ---
       ! ------ calculate in sbrtn. trcalv ( Non-zero values are required.)
       ! average density weighted charge of impurities
       avezimp(1)  = PZ(3)
       ! average density weighted atomic mass of impurities
       amassimp(1) = PA(3)
       ! average density weighted atomic mass of hydrogen ions
       amasshyd(1) = PA(2)

       ! average ion mass [AMU]   
       aimass(1)   = ai_ave(nr)
       ! -------------------------------
       
       wexbs(1) = wexbp(nr)

       grdne(1) = - rmjrhom * rn_ecl(nr) ! -R ( d n_e / d r ) / n_e
       grdni(1) = - rmjrhom * rn_icl(nr) ! -R ( d n_i / d r ) / n_i
       grdnh(1) = - rmjrhom * rn_icl(nr) ! -R ( d n_h / d r ) / n_h
       grdnz(1) = - rmjrhom * rn_icl(nr) ! -R ( d Z n_Z / d r ) / ( Z n_Z )
       grdte(1) = - rmjrhom * rt_ecl(nr) ! -R ( d T_e / d r ) / T_e
       grdti(1) = - rmjrhom * rt_icl(nr) ! -R ( d T_i / d r ) / T_i
       !  R ( d q   / d r ) / q    related to magnetic shear
       grdq (1) =   rmjrhom * qp_d(nr) / qpm

       CALL mmm95( &
! << Input arrays >> (jz: nr,  jm: mode)
  rminor,   &! minor radius (half-width) of zone boundary [m]
  rmajor,   &! major radius to geometric center of zone bndry [m]
  elong,    &! local elongation of zone boundary
!
  dense,    &! electron density [m^-3]
  densh,    &! sum over thermal hydrogenic ion densities [m^-3]
  densimp,  &! sum over impurity ion densities [m^-3]
  densfe,   &! electron density from fast (non-thermal) ions [m^-3]
!
  xzeff,    &! Z_eff
  tekev,    &! T_e (electron temperature) [keV] 
  tikev,    &! T_i (temperature of thermal ions) [keV]
  q,        &! magnetic q-value
  btor,     &! ( R B_tor ) / rmajor(jz) [tesla]
!
  avezimp,  &! average density weighted charge of impurities
!               = ( sum_imp n_imp Z_imp ) / ( sum_imp n_imp ) 
!               where sum_imp = sum over impurity ions with charge state Z_imp
!
  amassimp, &! average density weighted atomic mass of impurities
!               = ( sum_imp n_imp M_imp ) / ( sum_imp n_imp )
!               where sum_imp = sum over impurity ions, each with mass M_imp
!
  amasshyd, &! average density weighted atomic mass of hydrogen ions
!               = ( sum_hyd n_hyd M_hyd ) / ( sum_hyd n_hyd )
!               where sum_hyd = sum over hydrogenic ions, each with mass M_hyd
!
  aimass,   &! mean atomic mass of thermal ions [AMU]
!               = ( sum_i n_i M_i ) / ( sum_i n_i )
!               where sum_i = sum over all ions, each with mass M_i
!
  wexbs,    &! ExB shearing rate in [rad/s].  See  K.H. Burrell,
!            "Effects of {ExB} velocity shear and magnetic shear 
!             on turbulence and transport in magnetic confinement devices",
!             Phys. of Plasmas, 4, 1499 (1997).
!
!    All of the following normalized gradients are at zone boundaries.
!    r = half-width, R = major radius to center of flux surface
!
  grdne,    &! -R ( d n_e / d r ) / n_e
  grdni,    &! -R ( d n_i / d r ) / n_i
  grdnh,    &! -R ( d n_h / d r ) / n_h
  grdnz,    &! -R ( d Z n_Z / d r ) / ( Z n_Z )
  grdte,    &! -R ( d T_e / d r ) / T_e
  grdti,    &! -R ( d T_i / d r ) / T_i
  grdq ,    &!  R ( d q   / d r ) / q    related to magnetic shear
!
!  where:
!    n_i   : thermal ion density (sum over hydrogenic and impurity)
!    n_h   : thermal hydrogenic density (sum over hydrogenic species)
!    n_Z   : thermal impurity density,  Z = average impurity charge
!                      sumed over all impurities
!
! ----------------------------------------------------------------------------
! << Output array>>
!
!    The following effective diffusivities represent contributions
!    to the total diffusivity matrix (difthi and velthi given below)
!    from each of the models that contribute to the Multi-Mode model.
!    Generally, these arrays are used for diagnostic output only.
!
  thiig,    &! ion thermal diffusivity from the Weiland model
  thdig,    &! hydrogenic ion diffusivity from the Weiland model
  theig,    &! elelctron thermal diffusivity from the Weiland model
  thzig,    &! impurity ion diffusivity from the Weiland model
!	    
  thirb,    &! ion thermal diffusivity from resistive ballooning modes
  thdrb,    &! hydrogenic ion diffusivity from resistive ballooning modes
  therb,    &! elelctron thermal diffusivity from resistive ballooning modes
  thzrb,    &! impurity ion diffusivity from resistive ballooning modes
!	    
  thikb,    &! ion thermal diffusivity from kinetic ballooning modes
  thdkb,    &! hydrogenic ion diffusivity from kinetic ballooning modes
  thekb,    &! elelctron thermal diffusivity from kinetic ballooning modes
  thzkb,    &! impurity ion diffusivity from kinetic ballooning modes
!
!    The following are growth rates and mode frequencies from the
!    Weiland model for drift modes such as ITG and TEM.
!    These arrays are intended for diagnostic output.
!
  gamma,    &! gamma(jm,jz) growth rate for mode jm at point jz ( 1/sec )
  omega,    &! omega(jm,jz) frequency for mode jm at point jz ( rad/sec )
!
!    All of the transport coefficients are given in the following two
!    matricies for diffusion difthi and convection velthi in MKS units.
!    See the LaTeX documentation for difthi and velthi just below.
!
!    NOTE:  difthi and velthi include all of the anomalous transport.
!    There are no additional contributions to the heat fluxs from
!    charged particle convection.
!
  difthi,   &! (j1,j2,jz) full matrix of anomalous transport diffusivities
  velthi,   &! (j1,jz) convective velocities
  vflux,    &! (j1,jz) flux matrix
!
! ****************************************************************************
!  << Input integers >>
!
  matdim,   &! = first and second dimension of transport matricies
!            difthi(j1,j2,jz) and the first dimension of 
!            velthi(j1,jz), vflux(j1,jz), gamma(j1,jz), and omega(j1,jz).
!            matdim must be at least 5
!
  npoints,  &!= number of values of jz in all of the above arrays
!
  nprout,   &! = output unit number for long printout
!
!
!    +++ Input switches +++
!
  lprint,   &! controls the amount of printout (0 => no printout)
!              higher values yield more diagnostic output
!
  lsuper,   &!  = 0 for simulations of all other discharges
!           > 0 for supershot simulations; substantially reduces 
!               contribution from kinetic ballooning mode
!
  lreset,   &! = 0 to use only internal settings for lswitch, cswitch
!              and for the coefficients fig, frb, and fkb that control
!              the contributions form the various instability modes
!
!    Note that when lreset = 0, the values of the switches and
!    coefficients in the argument list are ignored and all the 
!    switches and coefficients are set internally.
!
!    WARNING:  use lreset > 0 only if you want to pass all the switches
!              lswitch, cswitch, fig, frb, and fkb through the 
!              argument list.
!
!    WARNING:  NTCC users desiring to use anything other than lreset = 0
!              should consult with the mmm95 code authors first.
!
! ----------------------------------------------------------------------------
!  << Output Integer >>
!
  nerr,      &!       status code returned; 0 = OK; .ne. 0 indicates error
!
!    +++ Internal control variables +++
!
  lswitch,   &!j=1,8   integer control variables: 
!
  cswitch,   &!j=1,25   general control variables:
!
!  lswitch(1)  controls which version of the Weiland model is used
!              Default lswitch(1) = 10
!             = 2  2 eqn  Weiland model Hydrogen \eta_i mode only
!             = 4  4 eqn  Weiland model with Hydrogen and trapped electrons
!             = 5  5 eqn  Weiland model with trapped electrons, FLR effects, 
!                         and parallel ion motion
!             = 6  6 eqn  Weiland model Hydrogen, trapped electrons,
!                    and one impurity species
!             = 7  7 eqn   Weiland model Hydrogen, trapped electrons,
!                  one impurity species, and collisions
!             = 8  8 eqn  Weiland model Hydrogen, trapped electrons,
!                  one impurity species, collisions, and parallel
!                  ion (hydrogenic) motion
!             = 9  9 eqn  Weiland model Hydrogen, trapped electrons,
!                  one impurity species, collisions, and finite beta
!             = 10 10 eqn Weiland model Hydrogen, trapped electrons,
!                  one impurity species, collisions, parallel
!                  ion (hydrogenic) motion, and finite beta
!             = 11 11 eqn Weiland model Hydrogen, trapped electrons,
!                  one impurity species, collisions, parallel
!                  ion (hydrogenic, impurity) motion, and finite beta
!
!  lswitch(2) = 0  full matrix representation for difthi and velthi
!                  Default lswitch(2) = 2
!             = 1  set diagonal matrix elements difthi and velthi
!             = 2  set diagonal matrix elements = effective diffusivities
!
!  lswitch(3)  controls \kappa scaling
!                  Default lswitch(3) = 0
!             = 0  use \kappa scaling raised to
!                  exponents (cswitch(3) - cswitch(5))
!             = 1  use (1+\kappa^2)/2 instead of \kappa scaling
!               
!  lswitch(4) > 0  to replace negative diffusivity with velocity
!                  Default lswitch(4) = 1
!
!  lswitch(5) = 1  to limit magnitude of all normalized gradients
!                  to ( major radius ) / ( ion Larmor radius )
!                  Default lswitch(5) = 1
!
!  cswitch(1)   0.5  minimum value of shear
!  cswitch(2)   3.5  coeff in exponential (fbeta-th) of kinetic ballooning model
!  cswitch(3)  -4.0  exponent of local elongation multiplying drift waves
!  cswitch(4)  -4.0  exponent of local elongation multiplying resistive
!                     ballooning modes
!  cswitch(5)  -4.0  exponent of local elongation multiplying
!                     kinetic balllooning modes
!  cswitch(6)   0.0  k_y \rho_s (= 0.316 if abs(cswitch(6)) < zepslon)
!  cswitch(8)   1.0  coeff of beta_prime_1 in kinetic ballooning mode
!  cswitch(9)  0.15  alpha in diamagnetic stabil. in kinetic ballooning model
!  cswitch(10)  0.0  rel fract of ion thermal diffusivity given to convection 
!  cswitch(11)  0.0  rel fract of hydrogen particle diffusivity given to convection 
!  cswitch(12)  0.0  rel fract of el thermal diffusivity given to convection 
!  cswitch(13)  0.0  rel fract of impurity particle diffusivity given to convection 
!  cswitch(14)  1.0  coef of finite beta effect in weiland14 = cetain(20) 
!  cswitch(15)  0.0  min value of impurity charge state zimpz
!  cswitch(16)  0.0  coef of fast particle fraction (superthermal ions) 
!                    in weiland model -- coef of densfe
!  cswitch(17)  1.0  coeff of k_\parallel (parallel ion motion) in weiland14 
!                    = cetain(10)
!  cswitch(18)  0.0  coeff of nuhat (effect of collisions) in weiland14 
!                    = cetain(15)
!  cswitch(19)  0.0  coeff for including v_parallel in strong ballooning limit
!                    = cetain(12); cswitch(19) = 1 for inclusion of v_par effect
!  cswitch(20)  0.0  trapping fraction used in weiland14 (when > 0.0)
!                    multiplies electron trapping fraction when < 0.0
!                    no effect when cswitch(20) = 0.0
!  cswitch(21)  1.0  multiplier for wexbs (flow shear rate) in Weiland model
!  cswitch(22)  0.0  ranges from 0.0 to 1.0 adds impurity heat flow to total 
!                    ionic heat flow for the weiland model
!  cswitch(23)  0.0  controls finite diff to construct the zgm matrix 
!                    = cetain(30)
!
!     contributions to vfluxes and interchanges: 
  fig,      &!
  frb,      &!
  fkb )      !
!
!  fig(1)   hydrogen particle transport from ITG (eta_i) mode
!  fig(2)   electron thermal  transport from ITG (eta_i) mode
!  fig(3)   ion      thermal  transport from ITG (eta_i) mode
!  fig(4)   impurity particle transport from ITG (eta_i) mode
!
!  frb(1)   hydrogen particle transport from resistive ballooning mode
!  frb(2)   electron thermal  transport from resistive ballooning mode
!  frb(3)   ion      thermal  transport from resistive ballooning mode
!  frb(4)   impurity particle transport from resistive ballooning mode
!
!  fkb(1)   hydrogen particle transport from kinetic ballooning mode
!  fkb(2)   electron thermal  transport from kinetic ballooning mode
!  fkb(3)   ion      thermal  transport from kinetic ballooning mode
!  fkb(4)   impurity particle transport from kinetic ballooning mode

       IF(nerr /= 0)THEN
          WRITE(*,*) 'error in subroutine "mmm95": nerr =',nerr
          STOP
       END IF

! --- for diagnostic ---
!       nrd1(nr) = theig(1)
!       nrd2(nr) = therb(1)
!       nrd3(nr) = thekb(1)

       ! in the case of isotropic (effective) diffusivities
       !  ** lswitch(2) = 2 **
       mmm_diffm(nr) = difthi(2,2,1)
       mmm_chiem(nr) = difthi(3,3,1)
       mmm_chiim(nr) = difthi(1,1,1)
       mmm_velem(nr) = velthi(3,1)
       mmm_velim(nr) = velthi(1,1)
       ! << NOTE >>
       ! 'difthi' and 'velthi' include all of the anomalous transport.
       ! There are no additional contributions to the heat flux 
       !  from charged particle convection.
    END DO

    DO nsa = 1, nsamax
       IF(idnsa(nsa) == -1)THEN ! electron
          dtr_tb(1+3*nsa-2,1+3*nsa-2,1:nrmax) = cdtrn * mmm_diffm(1:nrmax)
!          dtr_tb(1+3*nsa-1,1+3*nsa-1,1:nrmax) = cdtru *
          dtr_tb(1+3*nsa  ,1+3*nsa  ,1:nrmax) = cdtrt * mmm_chiem(1:nrmax)
          
          vtr_tb(1+3*nsa-2,1+3*nsa-2,1:nrmax) = cdtrn * 0.d0
!          vtr_tb(1+3*nsa-1,1+3*nsa-1,1:nrmax) = cdtru *
          vtr_tb(1+3*nsa  ,1+3*nsa  ,1:nrmax) = cdtrt * mmm_velem(1:nrmax)
       ELSE IF(idnsa(nsa) /= 0)THEN ! (hydrogenic) ion
          dtr_tb(1+3*nsa-2,1+3*nsa-2,1:nrmax) = cdtrn * mmm_diffm(1:nrmax)
!          dtr_tb(1+3*nsa-1,1+3*nsa-1,1:nrmax) = cdtru *
          dtr_tb(1+3*nsa  ,1+3*nsa  ,1:nrmax) = cdtrt * mmm_chiim(1:nrmax)

          vtr_tb(1+3*nsa-2,1+3*nsa-2,1:nrmax) = cdtrn * 0.d0
!          vtr_tb(1+3*nsa-1,1+3*nsa-1,1:nrmax) = cdtru *
          vtr_tb(1+3*nsa  ,1+3*nsa  ,1:nrmax) = cdtrt * mmm_velim(1:nrmax)
       END IF
    END DO

  END SUBROUTINE tr_mmm95

END MODULE trmmm95
