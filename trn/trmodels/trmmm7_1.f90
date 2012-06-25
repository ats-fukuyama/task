MODULE trmmm7_1

  PUBLIC tr_mmm7_1
  PRIVATE

CONTAINS

! *************************************************************************
!
!  The interface between TASK/TR(trcoeftb) and MMM7_1 model (mmm7_1)
!
! *************************************************************************

  SUBROUTINE tr_mmm7_1
    USE modmmm7_1
    USE trcomm, ONLY: &
         ikind,rkind,nrmax,nsamax,neqmax,idnsa,ns_nsa,pa,pz, &
         rmnrho,rmjrho,rkprho,abb1rho,BB,RR,rn,qp,           &
         vtor,vpol,vpar,cdtrn,cdtru,cdtrt,dtr_tb,vtr_tb

    USE trcalv, ONLY: &
         rt_e,rt_ecl,rt_i,rt_icl,rn_e,rn_ecl,rn_i,rn_icl, &
         qp_d,mshear,wexbp,z_eff

    IMPLICIT NONE

    ! Inputs array for mmm7_1 (subroutine): the dimension number is 'npoints'
    REAL(rkind),DIMENSION(1) :: &
         rmin, rmaj, elong, ne, nh, nz, nf, zeff, &
         te, ti, q, btor, zimp, aimp, ahyd, aimass, wexbs, &
         gne, gni, gnh, gnz, gte, gti, gq
    REAL(rkind),DIMENSION(1) :: &
         gvtor, mvtor, gvpol, mvpol, gvpar, mvpar
    REAL(rkind),DIMENSION(1) :: &
         gelong

    ! Outputs array for mmm7_1 (subroutine): the dimension number is 'npoints'
    REAL(rkind),DIMENSION(1) :: & ![m^2/s]
         xti, xdi, xte, xdz, xvt, xvp

    ! *** following argumets are for diagnostic outputs ***
    REAL(rkind),DIMENSION(1) :: & ![m^2/s]
         xtiW20, xdiW20, xteW20, xtiDBM, xdiDBM, xteDBM, xteETG
    REAL(rkind),DIMENSION(4,1) :: &
         omegaW20, gammaW20
    REAL(rkind),DIMENSION(1) :: &
         omegaDBM, gammaDBM

    REAL(rkind),DIMENSION(6,1) :: vconv
    REAL(rkind),DIMENSION(4,1) :: vflux

    ! input
    INTEGER(ikind) :: npoints, lprint, nprout
    REAL(rkind) :: rmaj0

    REAL(rkind),DIMENSION(3) :: cmodel

    REAL(rkind),DIMENSION(6,3)    :: cswitch
    INTEGER(ikind),DIMENSION(1,3) :: lswitch

    ! output
    INTEGER(ikind) :: nerr


    !  *** Internal variables for tr_mmm7_1 ***
    INTEGER(ikind) :: nr,ns,nsa
    REAL(rkind) :: sum_pan
    REAL(rkind) :: deriv3 ! in TASK/lib
    REAL(rkind),DIMENSION(0:nrmax) :: &
         mmm7_diff,mmm7_chie,mmm7_chii,mmm7_difz,mmm7_dtm,mmm7_dpm, &
         mmm7_veli,mmm7_vele,mmm7_velh
    REAL(rkind),DIMENSION(1:nrmax) :: &
         mmm7_diffm,mmm7_chiem,mmm7_chiim, &
         mmm7_velim,mmm7_velem,mmm7_velhm


    ! Initilization
    dtr_tb(1:neqmax,1:neqmax,0:nrmax) = 0.d0
    vtr_tb(1:neqmax,1:neqmax,0:nrmax) = 0.d0

    mmm7_diff(0:nrmax) = 0.d0
    mmm7_chie(0:nrmax) = 0.d0
    mmm7_chii(0:nrmax) = 0.d0
    mmm7_velim(1:nrmax) = 0.d0
    mmm7_velhm(1:nrmax) = 0.d0
    mmm7_velem(1:nrmax) = 0.d0
    mmm7_difz(0:nrmax) = 0.d0
    mmm7_dtm(0:nrmax)  = 0.d0
    mmm7_dpm(0:nrmax)  = 0.d0

    
    ! << Input integers >>
    npoints = 1  ! number of values of jz in all of the arrays

    nprout  = 71 ! output unit number for long printout
                 !  : for error and diagnostic
    lprint  = 50  ! controls the amount of printout (0 => no printout)
                 !  higher values yield more diagnostic ouput

    OPEN(nprout,FILE='mmm7_1_err_output.txt')

    ! +++ Input switches +++
    ! cmodel(1)~cmodel(3) are assigned as the weights 
    !  for Weiland20, DRIBM and ETG models, respectively.
    cmodel(1) = 1.d0
    cmodel(2) = 1.d0
    cmodel(3) = 1.d0

    
    DO nr = 1, nrmax

       rmin(1)  = rmnrho(nr)
       rmaj(1)  = rmjrho(nr)
       elong(1) = rkprho(nr)
       rmaj0    = RR

       ne(1) = rn_e(nr)*1.d20
       nh(1) = rn_i(nr)*1.d20

       ! --- for the time being ---
       nz(1) = 1.d10 ! impurity ion density [m^-3]
       nf(1) = 1.d5 ! density from fast (non-thermal) ions [m^-3]

       zeff(1) = z_eff(nr)
       te(1)   = rt_e(nr)
       ti(1)   = rt_i(nr)
       q(1)    = qp(nr)
       btor(1) = abb1rho(nr)

       ! --- for the time being ---
       ! ------ calculate in sbrtn. trcalv ( Non-zero values are required.) 
       ! mean charge of impurities
       zimp(1) = pz(3)
       ! mean atomic mass of impurities
       aimp(1) = pa(3)
       ! mean atomic mass of hydrogenic ions
       ahyd(1) = pa(2)


       sum_pan = 0.d0
       DO nsa = 1, nsamax
          IF(idnsa(nsa)==1)THEN
             ns = ns_nsa(nsa)
             sum_pan = sum_pan+pa(ns)*rn(nsa,nr)
          END IF
       END DO
       ! mean atomic mass of thermal ions
       aimass(1) = sum_pan / rn_i(nr)

       ! w_exb shearing rate [rad/s] [Phys. of Plasmas, 4, 1499 (1997)]
       wexbs(1) = wexbp(nr)

       gne(1) = - rmjrho(nr) * rn_ecl(nr) ! -R ( d n_e / d r) / n_e
       gni(1) = - rmjrho(nr) * rn_icl(nr) ! -R ( d n_i / d r) / n_i
       gnh(1) = - rmjrho(nr) * rn_icl(nr) ! -R ( d n_h / d r) / n_h
       gnz(1) = - rmjrho(nr) * rn_icl(nr) ! -R ( d Z n_z / d r) / ( Z n_z )
       gte(1) = - rmjrho(nr) * rt_ecl(nr) ! -R ( d T_e / d r) / T_e
       gti(1) = - rmjrho(nr) * rt_icl(nr) ! -R ( d T_i / d r) / T_i
       ! R ( d q / d r) / q : related to magnetic shear
       gq(1)  = rmjrho(nr) * qp_d(nr) / qp(nr)


       ! profiles related to momentum transport in 2006 Weiland model
       ! optional arguments of subrouitne mmm7_1
       !  --- for the time being ---
!!$       gvtor(1) = 0.d0!rmnrho(nr)/vtor(nr) * deriv3(nr,rmnrho,vtor,nrmax,0)
!!$       mvtor(1) = 0.d0!vtor(nr)
!!$       gvpol(1) = rmnrho(nr)/vpol(nr) * deriv3(nr,rmnrho,vpol,nrmax,0)
!!$       mvpol(1) = vpol(nr)
!!$       gvpar(1) = rmnrho(nr)/vpar(nr) * deriv3(nr,rmnrho,vpar,nrmax,0)
!!$       mvpar(1) = vpar(nr)
       gvtor(1) = 0.d0
       mvtor(1) = 0.d0
       gvpol(1) = 0.d0
       mvpol(1) = 0.d0
       gvpar(1) = 0.d0
       mvpar(1) = 0.d0
       
       ! Elongation gradient (optional argument)
       !  --- for the time being ---
       gelong(1) = 0.d0


       ! set optional arguments 'cswitch' and 'lswitch'
       !  for now, setting the defalut values
       CALL set_mmm7_1_switches( &
            cmmm = cswitch, lmmm = lswitch,   &
            KW20_C_EXB             = 1.d0,    &
            KW20_C_MOM_PINCH_SCALE = 1.d0,    &
            KW20_C_XTE_MIN         = 1.d-4,   &
            KW20_C_XTE_MAX         = 100.d0,  &
            KW20_C_XTI_MIN         = 1.d-4,   &
            KW20_C_XTI_MAX         = 100.d0,  &
            KDBM_C_EXB             = 0.d0,    &
            KETG_C_CEES_SCALE      = 0.06d0,  &
            KETG_C_CEEM_SCALE      = 0.06d0,  &
            KETG_L_NLTHR           = 1        &
            )

  CALL mmm7_1 ( &
       ! << Physical inputs >>
rmin=rmin,  &!: minor radius (half-width) of zone boundary [m]
rmaj=rmaj,  &!: major radius to geometric center of zone boundary [m]
rmaj0=rmaj0, &!: major radius at plasma axis (scalar) [m]
elong=elong, &!: local elongation of zone boundary
ne=ne,    &!: electron density [m^-3]
nh=nh,    &!: sum over thermal hydrogenic ion densities [m^-3]
nz=nz,    &!: sum over impurity ion densities [m^-3]
nf=nf,    &!: electron density from fast (non-thermal) ions [m^-3]
zeff=zeff,  &!: Z_eff
te=te,    &!: T_e (electron temperature) [keV]
ti=ti,    &!: T_i (temperature of thermal ions) [keV]
q=q,     &!: magnetic q-value
btor=btor,  &!: ( R B_tor ) / rmaj(jz)  [tesla]
zimp=zimp,  &!: average density weighted charge of impurities
               !   ( sum_imp n_imp Z_imp ) / ( sum_imp n_imp ) where
               !   sum_imp = sum over impurity ions with charge state Z_imp
aimp=aimp,  &!: average density weighted atomic mass of impurities
               !   ( sum_imp n_imp M_imp ) / ( sum_imp n_imp ) where
               !   sum_imp = sum over impurity ions, each with mass M_imp
ahyd=ahyd,  &!: average density weighted atomic mass of hydrogen ions
               !   ( sum_hyd n_hyd M_hyd ) / ( sum_hyd n_hyd ) where
               !   sum_hyd = sum over hydrogenic ions, each with mass M_hyd
aimass=aimass,&!: mean atomic mass of thermal ions [AMU]
               !   ( sum_i n_i M_i ) / ( sum_i n_i ) where
               !   sum_i = sum over all ions, each with mass M_i
wexbs=wexbs, &!: ExB shearing rate in [rad/s].  See  K.H. Burrell,
               !  "Effects of {ExB} velocity shear and magnetic shear
               !   on turbulence and transport in magnetic confinement
               !   devices", Phys. of Plasmas, 4, 1499 (1997).
gne=gne,   &!: -R ( d n_e / d r ) / n_e
gni=gni,   &!: -R ( d n_i / d r ) / n_i
gnh=gnh,   &!: -R ( d n_h / d r ) / n_h
gnz=gnz,   &!: -R ( d Z n_Z / d r ) / ( Z n_Z )
gte=gte,   &!: -R ( d T_e / d r ) / T_e
gti=gti,   &!: -R ( d T_i / d r ) / T_i
gq=gq,    &!:  R ( d q   / d r ) / q   related to magnetic shear
             ! where:
             ! n_i: thermal ion density (sum over hydrogenic and impurity)
             ! n_h: thermal hydrogenic density (sum over hydrogenic species)
             ! n_Z: thermal impurity density,  Z = average impurity charge
             !       summed over all impurities
vtor=mvtor,  &!: Toroidal velocity [m/s]
vpol=mvpol,  &!: Poloidal velocity [m/s]
vpar=mvpar,  &!: Parallel velocity [m/s]                 
gvtor=gvtor, &!: Normalized toroidal velocity grad. (r/v_tor)*(dv_tor/dr)
gvpol=gvpol, &!: Normalized poloidal velocity grad. (r/v_pol)*(dv_pol/dr)
gvpar=gvpar, &!: Normalized parallel velocity grad. (r/v_par)*(dv_par/dr)
gelong=gelong,&!: d elong / d rho  (elong is local elongation and rho=r/R)
       ! << Physical outputs >>
xti=xti,   &!: Effective ion thermal diffusivity
xdi=xdi,   &!: Effective hydrogenic ion diffusivity
xte=xte,   &!: Effective electron thermal diffusivity
xdz=xdz,   &!: Impurity ion diffusivity from the Weiland model
xvt=xvt,   &!: Toroidal momentum diffusivity from the Weiland model
xvp=xvp,   &!: Poloidal momentum diffusivity from the Weiland model
xtiW20=xtiW20,&!: Ion thermal diffusivity from the Weiland (W20) component
xdiW20=xdiW20,&!: Hydrogenic ion particle diffusivity from the Weiland component
xteW20=xteW20,&!: Electron thermal diffusivity from the Weiland component
xtiDBM=xtiDBM,&!: Ion thermal diffusivity from the DRIBM component
xdiDBM=xdiDBM,&!: Hydrogenic ion diffusivity from the DRIBM component
xteDBM=xteDBM,&!: Electron thermal diffusivity from the DRIBM component
xteETG=xteETG,&!: Electron thermal diffusivity from the Horton ETG component
gammaW20=gammaW20,&!: growth rate for the dominating mode at point jz (1/sec)
                 !   in Weiland
omegaW20=omegaW20,&!: frequency for dominating mode jm at point jz (rad/sec)
                 !   in Weiland
gammaDBM=gammaDBM,&!: growth rate for the dominating mode at point jz (1/sec)
                 !   of DRIBM mode
omegaDBM=omegaDBM,&!: frequency for dominating mode jm at point jz (rad/sec)
                 !   of DRIBM mode
vconv=vconv, &!: Convection velocities [m/s]
vflux=vflux, &!: Flux matrix
       ! << Feature controls >>
npoints=npoints,&!: Number of values in all of the 1-D arrays
lprint=lprint, &!: Vercose level
nprout=nprout, &!: Output unit number for long printout
nerr=nerr,   &!: Error return value
cmodel=cmodel, &!: Weights of internal models
                !: cmodel(1)~cmodel(3) are assigned as the weights 
                ! for Weiland20, DRIBM, ETG and Paleoclassical, respectively
cswitch=cswitch,&!: Real internal parameters
lswitch=lswitch )!: Integral internal parameters


      mmm7_chii(nr) = xti(1)
      mmm7_diff(nr) = xdi(1)
      mmm7_chie(nr) = xte(1)

      mmm7_veli(nr) = vconv(1,1) ! ion thermal convective velocity
      mmm7_velh(nr) = vconv(2,1) ! hydrogenic ion particle conv. vel.
      mmm7_vele(nr) = vconv(3,1) ! electron thermal convective velocity

      mmm7_difz(nr) = xdz(1)
      mmm7_dtm(nr)  = xvt(1)
      mmm7_dpm(nr)  = xvp(1)

   END DO ! End of nr loop

   ! on grid -> on half grid
    mmm7_chiim(1:nrmax) = 0.5d0*(mmm7_chii(0:nrmax-1)+mmm7_chii(1:nrmax))
    mmm7_diffm(1:nrmax) = 0.5d0*(mmm7_diff(0:nrmax-1)+mmm7_diff(1:nrmax))
    mmm7_chiem(1:nrmax) = 0.5d0*(mmm7_chie(0:nrmax-1)+mmm7_chie(1:nrmax))
    mmm7_velim(1:nrmax) = 0.5d0*(mmm7_veli(0:nrmax-1)+mmm7_veli(1:nrmax))
    mmm7_velhm(1:nrmax) = 0.5d0*(mmm7_velh(0:nrmax-1)+mmm7_velh(1:nrmax))
    mmm7_velem(1:nrmax) = 0.5d0*(mmm7_vele(0:nrmax-1)+mmm7_vele(1:nrmax))

    DO nsa = 1, nsamax
       IF(idnsa(nsa) == -1)THEN ! electron                                
          dtr_tb(1+3*nsa-2,1+3*nsa-2,1:nrmax) = cdtrn * mmm7_diffm(1:nrmax)
!          dtr_tb(1+3*nsa-1,1+3*nsa-1,1:nrmax) = cdtru *
          dtr_tb(1+3*nsa  ,1+3*nsa  ,1:nrmax) = cdtrt * mmm7_chiem(1:nrmax)

          vtr_tb(1+3*nsa-2,1+3*nsa-2,1:nrmax) = cdtrn * mmm7_velhm(1:nrmax)
!          vtr_tb(1+3*nsa-1,1+3*nsa-1,1:nrmax) = cdtru *
          vtr_tb(1+3*nsa  ,1+3*nsa  ,1:nrmax) = cdtrt * mmm7_velem(1:nrmax)
       ELSE IF(idnsa(nsa) /= 0)THEN ! (hydrogenic) ion                    
          dtr_tb(1+3*nsa-2,1+3*nsa-2,1:nrmax) = cdtrn * mmm7_diffm(1:nrmax)
!          dtr_tb(1+3*nsa-1,1+3*nsa-1,1:nrmax) = cdtru *
          dtr_tb(1+3*nsa  ,1+3*nsa  ,1:nrmax) = cdtrt * mmm7_chiim(1:nrmax)

          vtr_tb(1+3*nsa-2,1+3*nsa-2,1:nrmax) = cdtrn * mmm7_velhm(1:nrmax)
!          vtr_tb(1+3*nsa-1,1+3*nsa-1,1:nrmax) = cdtru *
          vtr_tb(1+3*nsa  ,1+3*nsa  ,1:nrmax) = cdtrt * mmm7_velim(1:nrmax)
       END IF
    END DO

  END SUBROUTINE tr_mmm7_1

END MODULE trmmm7_1
