!*************************** Sauter model ******************************

module sauter_mod
  implicit none
  real(8), parameter :: AEE  = 1.602176487D-19 ! elementary charge
  private
  public :: SAUTER

contains

!***********************************************************************
!
!       Sauter model
!
!***********************************************************************
  
  SUBROUTINE SAUTER(NE,TE,DTE,DPE_IN,NI,TI,DTI,DPI_IN, &
       &            Q,BB,DPSI,IPSI,EPS,RR,ZI,ZEFF,FT,rlnLei_IN,rlnLii_IN,NUE_IN,NUI_IN, &
       &            JBS,ETA,JBSB,ETAS)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Programmed by HONDA Mitsuru (2007/05/24)
!     
!     Note that all the parameters are on the GRID point.
!
!     *** INPUT VARIABLES ***
!     NE    : electron density [10^20/m^3]
!     TE    : electron temperature [keV]
!     DTE   : partial derivative of electron temperature to rho [keV/rho]
!     DPE   : partial derivative of electron pressure to rho [10^20 kev/m^3/rho]
!     NI    : ion density [10^20/m^3]
!     TI    : ion temperature [keV]
!     DTI   : partial derivative of ion temperature to rho [keV/rho]
!     DPI   : partial derivative of ion pressure to rho [10^20 kev/m^3/rho]
!     Q     : safety factor
!     BB    : toloidal magnetic field [T]
!     DPSI  : partial derivative of poloidal flux to rho [Wb/rho]
!     IPSI  : magnetic flux function, i.e. the product of major radius
!             and toroidal magnetic field (RR*BT) [mT]
!     EPS   : local inverse aspect ratio
!     RR    : major radius [m]
!     ZI    : charge number of bulk ion
!     ZEFF  : effective charge number
!     FT    : trapped particle fraction
!     rlnLei: Coulomb logarithm for electron collisions, optional
!     rlnLii: Coulomb logarithm for ion collisions,      optional
!     NUE   : Normalized collisionality for electrons,   optional
!     NUI   : Normalized collisionality for ions,        optional
!
!     *** OUTPUT VARIABLES ***
!     JBS   : bootstrap current [A/m^2]
!     JBSB  : parallel current <J . B> [AT/m^2],         optional
!     ETA   : neoclassical resistivity [Ohm m]
!     ETAS  : classical resistivity [Ohm m],             optional
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    real(8), intent(in) :: NE,TE,DTE,DPE_IN,NI,TI,DTI,DPI_IN,Q,BB,DPSI,IPSI,EPS,RR, &
    &                      ZI,ZEFF,FT
    real(8), intent(in), optional :: rlnLei_IN,rlnLii_IN,NUE_IN,NUI_IN
    real(8), intent(out) :: JBS, ETA
    real(8), intent(out), optional :: JBSB, ETAS

    real(8) :: EPSS,LnLame,NUE,LnLamii,NUI,RPE,ALFA0,ALFA,L31,L32,L34,RNZ,SGMSPTZ,QABS,&
         &     PE,PI,DPE,DPI
    real(8) :: F31TEFF,F32EETEFF,F32EITEFF,F34TEFF,F33TEFF

    ! On magnetic axis, the neoclassical resistivity reduces to the classical resistivity
    ! and the bootstrap current vanishes.
    IF(EPS == 0.D0) THEN
       RNZ = 0.58D0+0.74D0/(0.76D0+ZEFF)
       LnLame=31.3D0-LOG(SQRT(NE*1.D20)/ABS(TE*1.D3))
       IF(PRESENT(rlnLei_IN)) LnLame = rlnLei_IN

       JBS = 0.D0
       ETA = 1.D0 / (1.9012D4*(TE*1.D3)**1.5/(ZEFF*RNZ*LnLame))
       IF(PRESENT(JBSB)) JBSB = 0.D0
       IF(PRESENT(ETAS)) ETAS = ETA
       RETURN
    END IF

    EPSS = SQRT(EPS)**3
    QABS = ABS(Q)

    PE  = NE * TE * 1.D20 * AEE * 1.D3
    PI  = NI * TI * 1.D20 * AEE * 1.D3
    DPE = DPE_IN  * 1.D20 * AEE * 1.D3
    DPI = DPI_IN  * 1.D20 * AEE * 1.D3

!     LnLam : coulomb logarithm
!     NU    : collisional frequency [/s]

    LnLame=31.3D0-LOG(SQRT(NE*1.D20)/ABS(TE*1.D3))
    IF(PRESENT(rlnLei_IN)) LnLame = rlnLei_IN
    NUE=6.921D-18*QABS*RR*NE*1.D20*ZEFF*LnLame/(ABS(TE*1.D3)**2*EPSS)
    IF(PRESENT(NUE_IN))   NUE    = NUE_IN

    LnLamii=30.D0-LOG(ZI**3*SQRT(NI*1.D20)/(ABS(TI*1.D3)**1.5D0))
    IF(PRESENT(rlnLii_IN)) LnLamii = rlnLii_IN
    NUI = 4.90D-18*QABS*RR*NI*1.D20*ZI**4*LnLamii/(ABS(TI*1.D3)**2*EPSS)
    IF(PRESENT(NUI_IN))   NUI     = NUI_IN

!     RPE   : ratio of electron pressure to total pressure 

    RPE = PE/(PE+PI)

    F31TEFF   = FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(NUE)+0.5D0*(1.D0-FT)*NUE/ZEFF)
    F32EETEFF = FT/(1.D0+0.26D0*(1.D0-FT)*SQRT(NUE)+0.18D0*(1.D0-0.37D0*FT)*NUE/SQRT(ZEFF))
    F32EITEFF = FT/(1.D0+(1.D0+0.6D0*FT)*SQRT(NUE)+0.85D0*(1.D0-0.37D0*FT)*NUE*(1.D0+ZEFF))
    F34TEFF   = FT/(1.D0+(1.D0-0.1D0*FT)*SQRT(NUE)+0.5D0*(1.D0-0.5D0*FT)*NUE/ZEFF)

    ALFA0 = -1.17D0*(1.D0-FT)/(1.D0-0.22D0*FT-0.19D0*FT**2)
    ALFA  = ((ALFA0+0.25D0*(1.D0-FT**2)*SQRT(NUI)) &
         &     /(1.D0+0.5D0*SQRT(NUI))+0.315D0*NUI**2*FT**6)/(1.D0+0.15D0*NUI**2*FT**6)

    L31 = F31(F31TEFF,ZEFF)
    L32 = F32EE(F32EETEFF,ZEFF)+F32EI(F32EITEFF,ZEFF)
    L34 = F31(F34TEFF,ZEFF)

!     *** Bootstrap Current, JBS ***

    JBS = -IPSI*PE/DPSI/BB &
         & *( L31*(DPE/PE+DPI/PE)+L32*DTE/TE+L34*ALFA*(1.D0-RPE)/RPE*DTI/TI)
    if(present(jbsb)) JBSB = -IPSI*PE/DPSI &
         & *( L31*(DPE/PE+DPI/PE)+L32*DTE/TE+L34*ALFA*(1.D0-RPE)/RPE*DTI/TI)

!     *** Spitzer and Neoclassical Resistivity, ETAS, ETA ***

    F33TEFF = FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(NUE)+0.45D0*(1.D0-FT)*NUE/ZEFF**1.5)
    RNZ     = 0.58D0+0.74D0/(0.76D0+ZEFF)
    SGMSPTZ = 1.9012D4*(TE*1.D3)**1.5/(ZEFF*RNZ*LnLame)
    ETA  = 1.D0 / (SGMSPTZ * F33(F33TEFF,ZEFF))
    if(present(etas)) ETAS = 1.D0 / SGMSPTZ

  END SUBROUTINE SAUTER

!     *********************
!     *  Fitting Function *
!     *********************

  pure real(8) FUNCTION F33(X,Z)
    
    real(8), intent(in) :: X,Z

    F33 = 1.D0-(1.D0+0.36D0/Z)*X+0.59D0/Z*X**2-0.23D0/Z*X**3

  END FUNCTION F33

  pure real(8) FUNCTION F31(X,Z)

    real(8), intent(in) :: X,Z

    F31 = (1.D0+1.4D0/(Z+1.D0))*X-1.9D0/(Z+1.D0)*X**2 &
         &     +0.3D0/(Z+1.D0)*X**3+0.2D0/(Z+1.D0)*X**4

  END FUNCTION F31

  pure real(8) FUNCTION F32EE(X,Z)

    real(8), intent(in) :: X,Z

    F32EE = (0.05D0+0.62D0*Z)/(Z*(1.D0+0.44D0*Z))*(X-X**4) &
         &       +1.D0/(1.D0+0.22D0*Z)*(X**2-X**4-1.2D0*(X**3-X**4)) &
         &       +1.2D0/(1.D0+0.5D0*Z)*X**4 

  END FUNCTION F32EE

  pure real(8) FUNCTION F32EI(X,Z)

    real(8), intent(in) :: X,Z

    F32EI =-(0.56D0+1.93D0*Z)/(Z*(1.D0+0.44D0*Z))*(X-X**4) &
         &       +4.95D0/(1.D0+2.48D0*Z)*(X**2-X**4-0.55D0*(X**3-X**4)) &
         &       -1.2D0/(1.D0+0.5D0*Z)*X**4

  END FUNCTION F32EI

end module sauter_mod


module tx_nclass_mod
  implicit none
  private
  public :: TX_NCLASS

contains

!***********************************************************
!
!       NCLASS
!
!***********************************************************

  SUBROUTINE TX_NCLASS(NR,NueNC,NuiNC,Nue2NC,Nui2NC,ETAout,JBSout,ChiNCpel,ChiNCtel,ChiNCpil,ChiNCtil, &
       &               dTedr,dTidr,dPedr,dPidr,IER,dErdr,dBthdr,dErdr0,dBthdr0,p_gr2phi_in)
!****************************************************************************
!
!  Input : NR,dErdr,dBthdr,dTedr,dTidr,dPedr,dPidr
!          (optional) dErdr0,dBthdr0
!  Output: NueNC,NuiNC,Nue2NC,Nui2NC,ETAout,JBSout,ChiNCpel,ChiNCtel,ChiNCpil,ChiNCtil,IER
!
!****************************************************************************
!
!TX_NCLASS calculates various parameters and arrays for NCLASS.
!Please note that type declarations of all variables except "INTEGER" 
!  in NCLASS subroutine are "REAL(*4)" or "SINGLE" but not "REAL*8" 
!  or "DOUBLE".
!Input:
!  k_order - order of v moments to be solved (-)
!          =2 u and q
!          =3 u, q, and u2
!          =else error
!  k_potato - option to include potato orbits (-)
!           =0 off
!           =else on
!  m_i - number of isotopes (1<m_i<mx_mi+1)
!  m_z - highest charge state of all species (0<m_z<mx_mz+1)
!  c_den - density cutoff below which species is ignored (/m**3)
!  c_potb - kappa(0)*Bt(0)/[2*q(0)**2] (T)
!  c_potl - q(0)*R(0) (m)
!  p_b2 - <B**2> (T**2)
!  p_bm2 - <1/B**2> (/T**2)
!  p_eb - <E.B> (V*T/m)
!  p_fhat - mu_0*F/(dPsi/dr) (rho/m)
!  p_fm(3) - poloidal moments of geometric factor for PS viscosity (1/m**2)
!  p_ft - trapped particle fraction (-)
!  p_grbm2 - <grad(rho)**2/BB**2> (rho**2/m**2/T**2)
!  p_grphi - radial electric field Phi' (V/rho)
!  p_gr2phi - radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
!  p_ngrth - <n.grad(Theta)> (1/m)
!  amu_i(i) - atomic mass number of i (-)
!  grt_i(i) - temperature gradient of i (keV/rho)
!  temp_i(i) - temperature of i (keV)
!  den_iz(i,z) - density of i,z (/m**3)
!  fex_iz(3,i,z) - moments of external parallel force on i,z (T*j/m**3)
!  grp_iz(i,z) - pressure gradient of i,z (keV/m**3/rho)
!Using Output:
!  p_bsjb - <J_bs.B> (A*T/m**2)
!  p_etap - parallel electrical resistivity (Ohm*m)
!  p_exjb - <J_ex.B> current response to fex_iz (A*T/m**2)
!  bsjbp_s(s) - <J_bs.B> driven by unit p'/p of s (A*T*rho/m**3)
!  bsjbt_s(s) - <J_bs.B> driven by unit T'/T of s (A*T*rho/m**3)
!  gfl_s(m,s) - radial particle flux comps of s (rho/m**3/s)
!             m=1, banana-plateau, p' and T'
!             m=2, Pfirsch-Schluter
!             m=3, classical
!             m=4, banana-plateau, <E.B>
!             m=5, banana-plateau, external parallel force fex_iz
!  qfl_s(m,s) - radial heat conduction flux comps of s (W*rho/m**3)
!             m=1, banana-plateau, p' and T'
!             m=2, Pfirsch-Schluter
!             m=3, classical
!             m=4, banana-plateau, <E.B>
!             m=5, banana-plateau, external parallel force fex_iz
!  veb_s(s) - <E.B> particle convection velocity of s (rho/s)
!  qeb_s(s) - <E.B> heat convection velocity of s (rho/s)
!  ymu_s(s) - normalized viscosity for s (kg/m**3/s)
!  chip_ss(s1,s2) - heat cond coefficient of s2 on p'/p of s1 (rho**2/s)
!  chit_ss(s1,s2) - heat cond coefficient of s2 on T'/T of s1 (rho**2/s)
!  dp_ss(s1,s2) - diffusion coefficient of s2 on p'/p of s1 (rho**2/s)
!  dt_ss(s1,s2) - diffusion coefficient of s2 on T'/T of s1 (rho**2/s)
!  iflag - warning and error flag
!        =-4 warning: no viscosity
!        =-3 warning: no banana viscosity
!        =-2 warning: no Pfirsch-Schluter viscosity
!        =-1 warning: no potato orbit viscosity
!        = 0 no warnings or errors
!        = 1 error: order of v moments to be solved must be 2 or 3
!        = 2 error: number of species must be 1<m_i<mx_mi+1
!        = 3 error: number of species must be 0<m_z<mx_mz+1
!        = 4 error: number of species must be 1<m_s<mx_ms+1
!        = 5 error: inversion of flow matrix failed
!        = 6 error: trapped fraction must be 0.0.le.p_ft.le.1.0
!***********************************************************************
    use tx_commons
    use sauter_mod
    INCLUDE 'nclass/pamx_mi.inc'
    INCLUDE 'nclass/pamx_ms.inc'
    INCLUDE 'nclass/pamx_mz.inc'
    INCLUDE 'txncls.inc'
    INTEGER(4), INTENT(IN)  :: NR
    INTEGER(4), INTENT(OUT) :: IER
    real(8), intent(in)  :: dTedr,dTidr,dPedr,dPidr
    real(8), intent(in), optional :: dErdr,dBthdr,dErdr0,dBthdr0
    real(4), intent(in), optional :: p_gr2phi_in
    REAL(8), INTENT(OUT) :: NueNC, NuiNC, Nue2NC, Nui2NC, ETAout, JBSout
    INTEGER(4) :: i, k_out, k_v, ier_check
    REAL(4) :: a0, bt0, e0, p_eps, p_q, q0l, r0
    REAL(8) :: EpsL, BBL, PZMAX,  &
         &     PAL, PZL, RKAP, &
         &     ChiNCpel, ChiNCtel, ChiNCpil, ChiNCtil
    real(8) :: RL, BphVL, BthVL, EphVL, EthVL, QL, ErVL, PTeVL, PTiVL, &
         &     PNeVL, PNiVL, PeVL, PiVL, &
         &     dErdrL, dBthdrL, dTedrL, dTidrL, dPedrL, dPidrL
!!    real(8) :: p_fhat1, p_fhat2, p_fhat3, btot, uthai, VPOL(0:NRMAX), ppr, AJBSL, ETAL
!!    REAL(8) :: AITKEN2P

    !     *** Ellipticity on axis ***

    RKAP = 1.D0

    !     *** Dummy impurity for using high Zeff ***

!!$    PAL = 12.D0
!!$    PZL = 6.D0
    PAL = 56.D0 ! Fe
    PZL = 18.D0 ! Fe

    !     *** Initialization ***

    amu_i (1:mx_mi) = 0.0
    grt_i (1:mx_mi) = 0.0
    temp_i(1:mx_mi) = 0.0
    den_iz(1:mx_mi,1:mx_mz) = 0.0
    grp_iz(1:mx_mi,1:mx_mz) = 0.0
    fex_iz(1:3,1:mx_mi,1:mx_mz) = 0.0

    IER = 0

    !  k_out-option for output to nout (-)
    !       = 0 no display
    !        +1 warnings
    !        +2 errors
    !        +4 results
    !  k_v-option for neoclassical v_tor,v_pol,v_para,v_perp
    !       =1 output
    !       =else no output
    k_out    = 0
    k_v      = 1
    k_order  = 2
    k_potato = 1

    m_i      = 2
    IF(Zeff > 1.D0) m_i   = 3
    PZMAX    = PZ
    IF(Zeff > 1.D0) PZMAX = PZL
    m_z      = INT(PZMAX)
    c_den    = 1.E10
    !  *** Potate orbit factors ****************
    c_potb   = REAL(RKAP*BphV(0)/(2.D0*Q(0)**2))
    c_potl   = REAL(Q(0)*RR)
    !  *****************************************

    amu_i(1) = REAL(AME/AMP)
    amu_i(2) = REAL(PA)
    IF(Zeff > 1.D0) amu_i(3) = REAL(PAL)

    !***** !! IMPORTANT !! **********************************************!
    !  When NR == 0, the values are evaluated at the position            !
    !    intermediate between the magnetic axis and the adjacent mesh.   !
    !********************************************************************!
    
    IF (NR /= 0) THEN
       RL    = R(NR)
       BphVL = BphV(NR) ; BthVL = BthV(NR)
       EphVL = EphV(NR) ; EthVL = EthV(NR)
       QL    = Q   (NR) ; ErVL  = ErV (NR)
       PTeVL = PTeV(NR) ; PTiVL = PTiV(NR)
       PNeVL = PNeV(NR) ; PNiVL = PNiV(NR)
       PeVL  = PeV (NR) ; PiVL  = PiV (NR)
       dTedrL = dTedr   ; dTidrL  = dTidr
       dPedrL = dPedr   ; dPidrL  = dPidr
       if(present(dErdr)) then
          dErdrL = dErdr
       else
          dErdrL = 0.D0
       end if
       if(present(dBthdr)) then
          dBthdrL = dBthdr
       else
          dBthdrL = 0.D0
       end if
    ELSE
       RL    = 0.5D0 * R(NR+1)
       BphVL = 0.5D0 *(BphV(NR) + BphV(NR+1)) ; BthVL = 0.5D0 * BthV(NR+1)
       EphVL = 0.5D0 *(EphV(NR) + EphV(NR+1)) ; EthVL = 0.5D0 * EthV(NR+1)
       QL    = 0.5D0 *(Q   (NR) + Q   (NR+1)) ; ErVL  = 0.5D0 * ErV (NR+1)
       PTeVL = 0.5D0 *(PTeV(NR) + PTeV(NR+1)) ; PTiVL = 0.5D0 *(PTiV(NR) + PTiV(NR+1))
       PNeVL = 0.5D0 *(PNeV(NR) + PNeV(NR+1)) ; PNiVL = 0.5D0 *(PNiV(NR) + PNiV(NR+1))
       PeVL  = 0.5D0 *(PeV (NR) + PeV (NR+1)) ; PiVL  = 0.5D0 *(PiV (NR) + PiV (NR+1))
       dTedrL = 0.5D0 * dTedr   ; dTidrL  = 0.5D0 * dTidr
       dPedrL = 0.5D0 * dPedr   ; dPidrL  = 0.5D0 * dPidr
       if(present(dErdr)) then
          dErdrL = 0.5D0 *(dErdr + dErdr0)
       else
          dErdrL = 0.D0
       end if
       if(present(dBthdr)) then
          dBthdrL = 0.5D0 *(dBthdr + dBthdr0)
       else
          dBthdrL = 0.D0
       end if
    END IF

    BBL   = SQRT(BphVL**2 + BthVL**2)
    EpsL  = RL / RR
    p_b2  = REAL(BBL**2)
    p_bm2 = REAL(1.D0 / BBL**2)

    p_eb  = REAL(EphVL*BphVL + EthVL*BthVL)
!    p_eb  = REAL(EphV*BphV)
!!$    rlnLei(NR) = 37.8d0 - LOG(SQRT(PNeVL*1.D20)/(PTeVL))
!!$    rlnLii(NR) = 40.3d0 - LOG(PZ**2/PTiVL*SQRT(2.D0*PNiVL*1.D20*PZ**2/PTiVL))
!!$    CALL SAUTER(PNeVL,PTeVL,dTedrL,dPedrL,PNiVL,PTiVL,dTidrL,dPidrL, &
!!$         &      QL,BphVL,RR*RA*BthVL,RR*BphVL,EpsL,RR,PZ,Zeff,ft(nr), &
!!$         &      rlnLei_IN=rlnLei(NR),rlnLii_IN=rlnLii(NR),JBS=AJBSL,ETA=ETAL)
!!$    IF(NR == 0) AJBSL = 0.D0
!!$    p_eb  = REAL(ETAL*(( (-   AEE*PNeVL*UephV(NR) &
!!$         &                +PZ*AEE*PNiVL*UiphV(NR) &
!!$         &                +PZ*AEE*PNbV(NR)*UbphV(NR))*BphVL &
!!$         &              +(-   AEE*PNeVL*UethV(NR) &
!!$         &                +PZ*AEE*PNiVL*UithV(NR) &
!!$         &                +PZ*AEE*PNbV(NR)*UbthV(NR))*BthVL)*1.D20-AJBSL*BBL))
!!$    p_eb = REAL(ETAL*( (-   AEE*PNeVL*UephV(NR) &
!!$         &              +PZ*AEE*PNiVL*UiphV(NR) &
!!$         &              +PZ*AEE*PNbV(NR)*UbphV(NR))*BphVL &
!!$         &            +(-   AEE*PNeVL*UethV(NR) &
!!$         &              +PZ*AEE*PNiVL*UithV(NR) &
!!$         &              +PZ*AEE*PNbVL*UbthV(NR))*BthVL)*1.D20)
    ! No ohmic current causes infinite resistivity in the NCLASS module.
    IF(p_eb == 0.D0) p_eb = 1.e-10

!!$    IF(NR == 0) THEN
!!$       p_fhat1 = BphV(NR+1)/(RA*BthV(NR+1))
!!$       p_fhat2 = BphV(NR+2)/(RA*BthV(NR+2))
!!$       p_fhat3 = BphV(NR+3)/(RA*BthV(NR+3))
!!$       p_fhat  = REAL(AITKEN2P(R(0),p_fhat1,p_fhat2,p_fhat3,R(1),R(2),R(3)))
!!$    ELSE
       p_fhat  = REAL(BphVL/(RA*BthVL))
!!$    END IF

    !  Approximation inverse aspect ratio at the magnetix axis
!!$    IF(NR == 0) EpsL = 0.01D0*R(NR+1)/RR
    IF(REAL(EpsL) > 0.0) THEN
       ! poloidal moments of geometric factor for PS viscosity
       DO i=1,3
          p_fm(i)=REAL(DBLE(i)*( (1.D0-SQRT(1.D0-EpsL**2))/EpsL)**(2*i) &
               &                *(1.D0+DBLE(i)*SQRT(1.D0-EpsL**2))/((1.D0-EpsL**2)**1.5D0 &
               &                *(QL*RR)**2))
       ENDDO
    ENDIF
!!    p_fm(1:3) = 0.0 ! No Pfirsch-Schulter viscosity
    p_ft=REAL(1.46D0 * SQRT(EpsL) - 0.46 * EpsL * SQRT(EpsL))

    p_grbm2   = REAL(1.D0/RA**2)
    p_grphi   = REAL(-RA*ErVL)
!    p_gr2phi  = REAL(-RA**2*dErdrL) ! Orbit squeezing
    ! For orbit squeezing (Houlberg, PoP, 1997, Eq. (B2))
    if(present(p_gr2phi_in)) then
       p_gr2phi = p_gr2phi_in
    else
!!$       if(nr == 0) then
!!$          p_gr2phi = 0.0 ! Any value is OK because the value at nr=0 is discarded.
!!$       else
          p_gr2phi  = REAL(-RA**2*dErdrL+RA**2*ErVL*dBthdrL/BthVL)
!!$          if(nt==ntmax)write(6,*) Rho(NR),-RA**2*dErdrL,RA**2*ErVL*dBthdrL/BthVL
!!$          if(nr == 1) write(6,*) T_TX,p_gr2phi
!!$       end if
    end if
    p_ngrth   = REAL(BphVL/(RR*QL*BBL))
    temp_i(1) = REAL(PTeVL)
    temp_i(2) = REAL(PTiVL)
    grt_i(1)  = REAL(RA * dTedrL)
    grt_i(2)  = REAL(RA * dTidrL)
    den_iz(1,1)       = REAL(PNeVL) * 1.E20
    den_iz(2,INT(PZ)) = REAL(PNiVL) * 1.E20
    grp_iz(1,1)       = REAL(RA * dPedrL) * 1.E20
    grp_iz(2,INT(PZ)) = REAL(RA * dPidrL) * 1.E20
    IF (Zeff > 1.D0) THEN
       temp_i(3) = temp_i(2)
       grt_i(3)  = grt_i(2)
       den_iz(2,INT(PZ))  = REAL((PZL*PZ-Zeff)/(PZ*(PZL-PZ))*PNiVL) * 1.E20
       den_iz(3,INT(PZL)) = REAL((Zeff-PZ**2)/(PZL*(PZL-PZ))*PNiVL) * 1.E20
       grp_iz(2,INT(PZ))  = REAL(RA * dPidrL * (PZL*PZ-Zeff)/(PZ*(PZL-PZ))) * 1.E20
       grp_iz(3,INT(PZL)) = REAL(RA * dPidrL * (Zeff-PZ**2)/(PZL*(PZL-PZ))) * 1.E20
    END IF

    ! Even when NBI is activated, any external parallel force concerning NBI
    ! never acts on electrons and bulk ions. Force term due to NBI only appear
    ! in LQb4.
    fex_iz(1:3,1:mx_mi,1:mx_mz) = 0.0

    CALL NCLASS( &
!     Input
         &      k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2, &
         &      p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi, &
         &      p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz, &
!     Output
         &      m_s,jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii, &
         &      capm_ii,capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s, &
         &      sqz_s,upar_s,utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s, &
         &      chip_ss,chit_ss,dp_ss,dt_ss,iflag)

    IF(k_out >0 .and.k_out <= 7) THEN
       p_eps = REAL(RL/RR)
       p_q   = REAL(QL)
       r0    = REAL(RR)
       a0    = REAL(RA)
       e0    = REAL(1.D0)
       bt0   = REAL(BphV(0))
       q0l   = REAL(Q(0))
       CALL NCLASS_CHECK(6,NR,k_out,k_order,m_i,m_z, &
            &        p_fhat,p_grphi,amu_i,grt_i,temp_i,den_iz,grp_iz,p_eps,bt0, &
            &        m_s,jm_s,jz_s,p_bsjb, &
            &        p_etap,p_exjb,calm_i,caln_ii, &
            &        capm_ii,capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s, &
            &        upar_s,utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s, &
            &        chip_ss,chit_ss,dp_ss,dt_ss,iflag,ier_check)
       IER = ier_check
    ENDIF

!    write(6,'(4F18.7)') r(nr)/ra,upar_s(1,1,1)/BBL,upar_s(1,2,1)/BBL,sum(upar_s(1,1:3,1))/BBL
!    write(6,'(4F18.7)') r(nr)/ra,upar_s(1,1,2)/BBL,upar_s(1,2,2)/BBL,sum(upar_s(1,1:3,2))/BBL

    !     *** Takeover Parameters ***

    IF(k_potato == 0) THEN
       IF(iflag /= -1 .and. k_out > 0) THEN
          WRITE(6,*) "XX iflag=",iflag
          IER=1
          RETURN
       ENDIF
    ELSE
       IF(iflag /= 0 .and. k_out > 0) WRITE(6,*) "XX iflag=",iflag
    ENDIF

    !   Bootstrap current
    IF(Zeff == 1.D0) THEN
       JBSout =-(  DBLE(bsjbt_s(1)) *(RA * dTedrL/ PTeVL) &
              &  + DBLE(bsjbp_s(1)) *(RA * dPedrL/ PeVL) &
              &  + DBLE(bsjbt_s(2)) *(RA * dTidrL/ PTiVL) &
              &  + DBLE(bsjbp_s(2)) *(RA * dPidrL/ PiVL)) / BBL
    ELSE
       JBSout =-(  DBLE(bsjbt_s(1)) *(RA * dTedrL/ PTeVL) &
              &  + DBLE(bsjbp_s(1)) *(RA * dPedrL/ PeVL) &
              &  + DBLE(bsjbt_s(2)) *(RA * dTidrL/ PTiVL) &
              &  + DBLE(bsjbp_s(2)) *(RA * dPidrL/ PiVL * (PZL*PZ-Zeff)/(PZ*(PZL-PZ))) &
              &  + DBLE(bsjbt_s(3)) *(RA * dTidrL/ PTiVL) &
              &  + DBLE(bsjbp_s(3)) *(RA * dPidrL/ PiVL * (Zeff-PZ**2)/(PZL*(PZL-PZ)))) / BBL
    END IF
    IF(k_potato == 0 .and. NR == 0) JBSout = 0.D0

    !   Neoclassical resistivity
    ETAout = DBLE(p_etap)

    !   Neoclassical viscosity (so-called "Heuristic closure")
    !     T.A. Gianakon, S.E. Kruger, C.C. Hegna, PoP 9 (2002) 536, Eq.(12)
    !     D.D. Schnack, et al., PoP 13 (2006) 058103, Eqs.(11),(88) and (89)
    NueNC  = FSNC * DBLE(p_b2 * ymu_s(1,1,1)) / (PNeVL * 1.D20 * AME * BthVL**2)  ! [/s]
    NuiNC  = FSNC * DBLE(p_b2 * ymu_s(1,1,2)) / (PNiVL * 1.D20 * AMI * BthVL**2)  ! [/s]
    Nue2NC = FSNC * DBLE(ymu_s(1,2,1)) * BphVL / (PNeVL * 1.D20 * AME * BthVL**2) ! [/Ts]
    Nui2NC = FSNC * DBLE(ymu_s(1,2,2)) * BphVL / (PNiVL * 1.D20 * AMI * BthVL**2) ! [/Ts]

    !   Neoclassical thermal diffusivity (Diagonal effect only)
    ChiNCpel = ChiNC * DBLE(chip_ss(1,1))
    ChiNCtel = ChiNC * DBLE(chit_ss(1,1))
    ChiNCpil = ChiNC * DBLE(chip_ss(2,2))
    ChiNCtil = ChiNC * DBLE(chit_ss(2,2))
    if(ChiNCpel < 0.d0) ChiNCpel = 0.d0
    if(ChiNCtel < 0.d0) ChiNCtel = 0.d0
    if(ChiNCpil < 0.d0) ChiNCpil = 0.d0
    if(ChiNCtil < 0.d0) ChiNCtil = 0.d0

!!$    !   Now not using
!!$       DO i=1,m_s
!!$          uthai=DBLE(utheta_s(1,1,i)+utheta_s(1,2,i)+utheta_s(1,3,i))
!!$          IF(DBLE(amu_i(jm_s(i))) == PA .and. jz_s(i) == INT(PZ)) then
!!$             !     Poloidal
!!$             VPOL(NR)=uthai*BthVL
!!$!             !     Parallel
!!$!             VPAR(NR)=DBLE(BthVL/btot*VPOL(NR)+BphVL/btot*VTOR(NR))
!!$!             !     Perpendicular
!!$!             VPRP(NR)=DBLE(BphVL/btot*VPOL(NR)-BthVL/btot*VTOR(NR))
!!$          ENDIF
!!$       ENDDO
!!$!       IF(MOD(NT,1000)==0) write(6,'(5F15.7)') RL/RA,UithV(NR),VPOL(NR),utheta_s(1,1,2)*BthVL,utheta_s(1,2,2)*BthVL

    RETURN
  end SUBROUTINE TX_NCLASS

!***********************************************************
!
!    PARAMETER AND CONSISTENCY CHECK
!
!***********************************************************

  SUBROUTINE NCLASS_CHECK(nout,nr,k_out,k_order,m_i,m_z, &
       &            p_fhat,p_grphi,amu_i,grt_i,temp_i,den_iz,grp_iz,p_eps,bt0, &
       &            m_s,jm_s,jz_s,p_bsjb, &
       &            p_etap,p_exjb,calm_i,caln_ii, &
       &            capm_ii,capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s, &
       &            upar_s,utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s, &
       &            chip_ss,chit_ss,dp_ss,dt_ss,iflag,ier)

    INCLUDE 'nclass/pamx_mi.inc'
    INCLUDE 'nclass/pamx_ms.inc'
    INCLUDE 'nclass/pamx_mz.inc'
    INCLUDE 'txncls.inc'
    !Declaration of local variables
    CHARACTER :: label*120
    INTEGER(4), intent(in)  :: nr, nout, k_out
    integer(4), intent(out) :: ier
    INTEGER(4) :: i, im, iz, iza, j, jm, jza, k, l
    INTEGER(4) :: idum(8)
    REAL(4)    :: bt0, bpol, btor, btot, p_eps, ppr, uthai
    REAL(4)    :: dq_s(mx_ms), vq_s(mx_ms)
    REAL(4)    :: z_coulomb, z_electronmass, z_j7kv, z_mu0, z_pi, z_protonmass
    REAL(4)    :: dum(8), edum(8), rdum(8)
    !Declaration of functions
    REAL(4)    ::  RARRAY_SUM

    !Physical and conversion constants
    z_coulomb=1.6022e-19
    z_electronmass=9.1095e-31
    z_j7kv=1.6022e-16
    z_mu0=1.2566e-06
    z_pi=ACOS(-1.0)
    z_protonmass=1.6726e-27

    ier = 0

    IF(mod(k_out,2) == 1) THEN
       IF(iflag /= 0) WRITE(nout,'(A3,I3)') "NR=",NR
       !Check warning flags
       IF(iflag == -1) THEN
          label='WARNING:NCLASS-no potato orbit viscosity'
          CALL WRITE_LINE(nout,label,0,0)
       ELSEIF(iflag == -2) THEN
          label='WARNING:NCLASS-Pfirsch-Schluter viscosity'
          CALL WRITE_LINE(nout,label,0,0)
       ELSEIF(iflag == -3) THEN
          label='WARNING:NCLASS-no banana viscosity'
          CALL WRITE_LINE(nout,label,0,0)
       ELSEIF(iflag == -4) THEN
          label='WARNING:NCLASS-no viscosity'
          CALL WRITE_LINE(nout,label,0,0)
       ENDIF
    END IF
    !Check error flags
    IF(iflag > 0) THEN
       IF(int(k_out/2) == 1 .or. int(k_out/2) == 3) THEN
          IF(iflag == 1) THEN
             label='ERROR:NCLASS-k_order must be 2 or 3, k_order='
             idum(1)=k_order
             CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag == 2) THEN
             label='ERROR:NCLASS-require 1<m_i<mx_mi, m_i='
             idum(1)=m_i
             CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag == 3) THEN
             label='ERROR:NCLASS-require 0<m_z<mx_mz, m_z='
             idum(1)=m_z
             CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag == 4) THEN
             label='ERROR:NCLASS-require 0<m_s<mx_ms, m_s='
             idum(1)=m_s
             CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          ELSEIF(iflag == 5) THEN
             label='ERROR:NCLASS-inversion of flow matrix failed'
             CALL WRITE_LINE(nout,label,0,0)
          ELSEIF(iflag == 6) THEN
             label='ERROR:NCLASS-trapped fraction must be 0.0<p_ft<1.0'
             CALL WRITE_LINE(nout,label,0,0)
          ENDIF
       ENDIF
       WRITE(6,*) 'XX NCLASS_CHECK: non zero iflag: ',iflag
       ier = 1
       return
    ENDIF
    !Check for optional output
    IF(k_out >= 4) THEN
       !  Species identification
       label='     *** Species Identification ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='     Isotope     Species      Charge        Mass'// &
            &'     Density Temperature  Chg Factor'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -           -     Coulomb         AMU'// &
            &'       /m**3         keV           -'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          im=jm_s(i)
          iz=jz_s(i)
          iza=IABS(iz)
          idum(1)=im
          idum(2)=i
          idum(3)=iz
          rdum(1)=amu_i(im)
          rdum(2)=den_iz(im,iza)
          rdum(3)=temp_i(im)
          rdum(4)=xi_s(i)
          CALL WRITE_IR(nout,3,idum,4,rdum,2)
       ENDDO
       !  Friction coefficients
       label='     *** Friction Coefficients ***'
       CALL WRITE_LINE(nout,label,2,0)
       DO im=1,m_i
          DO jm=1,m_i
             CALL WRITE_LINE(nout,' ',0,0)            
             label='  Isotopes ='
             idum(1)=im
             idum(2)=jm
             CALL WRITE_LINE_IR(nout,label,2,idum,0,rdum,0)           
             !           Mkl
             label='  Mkl(-) ='
             CALL WRITE_LINE(nout,label,0,0)           
             DO k=1,k_order
                DO l=1,k_order
                   rdum(l)=capm_ii(k,l,im,jm)
                ENDDO
                CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
             ENDDO
             !           Nkl
             label='  Nkl(-) ='
             CALL WRITE_LINE(nout,label,0,0)           
             DO k=1,k_order
                DO l=1,k_order
                   rdum(l)=capn_ii(k,l,im,jm)
                ENDDO
                CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
             ENDDO
          ENDDO
       ENDDO
       !  Reduced friction coefficients
       label='     *** Reduced Friction Coefficients ***'
       CALL WRITE_LINE(nout,label,2,0)
       !       cal(M)
       DO im=1,m_i
          CALL WRITE_LINE(nout,' ',0,0)
          label='  calMij (kg/m**3/s) for Isotope ='
          idum(1)=im
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)           
          DO j=1,k_order
             DO k=1,k_order
                rdum(k)=calm_i(j,k,im)
             ENDDO
             CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
          ENDDO
       ENDDO
       !       cal(N)
       DO im=1,m_i
          CALL RARRAY_ZERO(3,dum)
          DO jm=1,m_i
             CALL WRITE_LINE(nout,' ',0,0)
             label='  calNij (kg/m**3/s) for Isotopes ='
             idum(1)=im
             idum(2)=jm
             CALL WRITE_LINE_IR(nout,label,2,idum,0,rdum,0)           
             DO k=1,k_order
                dum(k)=dum(k)-caln_ii(k,1,im,jm)
                DO l=1,k_order
                   rdum(l)=caln_ii(k,l,im,jm)
                ENDDO
                CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
             ENDDO
          ENDDO
          !         Momentum check
          CALL WRITE_LINE(nout,' ',0,0)
          label='Momentum check for Isotope ='
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          label='  sum_b[-calNk1] (kg/m**3/s) ='
          CALL WRITE_LINE_IR(nout,label,0,idum,k_order,dum,2)           
          label='  calMk1 (kg/m**3/s)         ='
          CALL WRITE_LINE_IR(nout,label,0,idum,k_order,calm_i(1,1,im),2)           
          label='  Difference ='
          DO k=1,k_order
             dum(k+3)=dum(k)-calm_i(k,1,im)
          ENDDO
          CALL WRITE_LINE_IR(nout,label,0,idum,k_order,dum(4),2)           
       ENDDO
       !  Normalized viscosities
       label='     *** Normalized Viscosities ***'
       CALL WRITE_LINE(nout,label,2,0)
       DO i=1,m_s
          CALL WRITE_LINE(nout,' ',0,0)
          label='  muij (kg/m**3/s) for Species ='
          idum(1)=i
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          DO k=1,k_order
             DO l=1,k_order
                rdum(l)=ymu_s(k,l,i)
             ENDDO
             CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
          ENDDO
       ENDDO
       !  Normalized parallel flows
       label='     *** Normalized Parallel Flows ***'
       CALL WRITE_LINE(nout,label,2,0)
       DO i=1,m_s
          CALL WRITE_LINE(nout,' ',0,0)
          label='  upar_ij (T*m/s) for Species ='
          idum(1)=i
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          DO j=1,k_order
             DO k=1,k_order
                rdum(k)=upar_s(j,k,i)
             ENDDO
             CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
          ENDDO
       ENDDO
       !  Normalized poloidal flows
       label='     *** Normalized Poloidal Flows ***'
       CALL WRITE_LINE(nout,label,2,0)
       DO i=1,m_s
          CALL WRITE_LINE(nout,' ',0,0)
          label='  utheta_ij (m/s/T) for Species ='
          idum(1)=i
          CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
          DO j=1,k_order
             DO k=1,k_order
                rdum(k)=utheta_s(j,k,i)
             ENDDO
             CALL WRITE_IR(nout,0,idum,k_order,rdum,2)           
          ENDDO
       ENDDO
       !  Radial particle fluxes
       label='     *** Radial Particle Fluxes ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species          BP          PS          CL'// &
            &'       <E.B>         src       total'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -     /m**2/s     /m**2/s     /m**2/s'// &
            &'     /m**2/s     /m**2/s     /m**2/s'
       CALL WRITE_LINE(nout,label,0,0)
       CALL RARRAY_ZERO(6,dum)
       CALL RARRAY_ZERO(6,edum)
       !       Load into rdum the five flux components and total
       !       Load into dum the z-weighted + charge (ion) components
       !       Load into edum the z-weighted - charge (electron) components 
       DO i=1,m_s
          iz=jz_s(i)
          idum(1)=i
          CALL RARRAY_COPY(5,gfl_s(1,i),1,rdum,1)
          rdum(6)=RARRAY_SUM(5,rdum,1)
          CALL WRITE_IR(nout,1,idum,6,rdum,2)           
          DO k=1,5
             IF(iz > 0) THEN
                dum(k)=dum(k)+iz*gfl_s(k,i)
             ELSE
                edum(k)=edum(k)+iz*gfl_s(k,i)
             ENDIF
          ENDDO
          IF(iz > 0) THEN
             dum(6)=dum(6)+iz*rdum(6)
          ELSE
             edum(6)=edum(6)+iz*rdum(6)
          ENDIF
       ENDDO
       !  Ambipolarity check
       label='     *** Ambipolarity Check ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='  Sum_i Z_i*Gamma_i ='
       CALL WRITE_LINE(nout,label,0,0)
       idum(1)=0
       CALL WRITE_IR(nout,1,idum,6,dum,2)           
       label='  Sum_e Z_e*Gamma_e ='
       CALL WRITE_LINE(nout,label,0,0)
       CALL WRITE_IR(nout,1,idum,6,edum,2)           
       label='  Difference ='
       CALL WRITE_LINE(nout,label,0,0)
       DO k=1,6
          rdum(k)=edum(k)+dum(k)
       ENDDO
       CALL WRITE_IR(nout,1,idum,6,rdum,2)
       !  Particle transport is returned in three forms
       !  Particle flux consistency check
       label='     *** Particle Flux Consistency Check ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species       gfl_s   dn_s,vn_s   dp_s,dt_s'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -     /m**2/s     /m**2/s     /m**2/s'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          im=jm_s(i)
          iza=IABS(jz_s(i))
          idum(1)=i
          !    Gamma total is sum of components
          rdum(1)=RARRAY_SUM(4,gfl_s(1,i),1)
          !    Gamma total is diffusive plus pinch
          rdum(2)=-dn_s(i)*(grp_iz(im,iza)-den_iz(im,iza)*grt_i(im)) &
               & /temp_i(im)+den_iz(im,iza)*(vn_s(i)+veb_s(i))
          !    Gamma total is sum over T' and p' of all species
          rdum(3)=0.0
          DO j=1,m_s
             jm=jm_s(j)
             jza=IABS(jz_s(j))
             rdum(3)=rdum(3)-dt_ss(j,i)*grt_i(jm)/temp_i(jm) &
                  & -dp_ss(j,i)*grp_iz(jm,jza)/den_iz(jm,jza)/temp_i(jm)
          ENDDO
          rdum(3)=(rdum(3)+veb_s(i))*den_iz(im,iza)
          CALL WRITE_IR(nout,1,idum,3,rdum,2)
       ENDDO
       !  Particle diffusion, velocity
       label='     *** Particle Diffusion, Velocity ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species         D_s         V_s'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      m**2/s         m/s'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          im=jm_s(i)
          iza=IABS(jz_s(i))
          idum(1)=i
          rdum(1)=dn_s(i)
          rdum(2)=vn_s(i)+veb_s(i)+gfl_s(5,i)/den_iz(im,iza)
          CALL WRITE_IR(nout,1,idum,2,rdum,2)
       ENDDO
       !  Particle diffusivity matrices
       !    On p'/p           
       label='     *** Particle Diffusivity Matrices ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species    on dp/dr ... by Species'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      m**2/s ...'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          idum(1)=i
          CALL WRITE_IR(nout,1,idum,m_s,dp_ss(1,i),2)
       ENDDO
       !    On T'/T
       label='     Species    on dT/dr ... by Species'
       CALL WRITE_LINE(nout,label,1,0)
       label='           -      m**2/s ...'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          idum(1)=i
          CALL WRITE_IR(nout,1,idum,m_s,dt_ss(1,i),2)
       ENDDO
       !  Radial conduction fluxes
       label='     *** Radial Conduction Fluxes ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species          BP          PS          CL'// &
            &'       <E.B>         src       total'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      w/m**2      w/m**2      w/m**2'// &
            &'      w/m**2      w/m**2      w/m**2'
       CALL WRITE_LINE(nout,label,0,0)
       CALL RARRAY_ZERO(6,dum)
       CALL RARRAY_ZERO(6,edum)
       DO i=1,m_s
          iz=jz_s(i)
          idum(1)=i
          CALL RARRAY_COPY(5,qfl_s(1,i),1,rdum,1)
          rdum(6)=RARRAY_SUM(5,qfl_s(1,i),1)
          CALL WRITE_IR(nout,1,idum,6,rdum,2)           
          DO k=1,5
             IF(iz > 0) THEN
                dum(k)=dum(k)+qfl_s(k,i)
             ELSE
                edum(k)=edum(k)+qfl_s(k,i)
             ENDIF
          ENDDO
          IF(iz > 0) THEN
             dum(6)=dum(6)+rdum(6)
          ELSE
             edum(6)=edum(6)+rdum(6)
          ENDIF
       ENDDO
       !     Total electron radial conduction flux
       label='  Sum of electron conduction fluxes ='
       CALL WRITE_LINE(nout,label,1,0)
       idum(1)=0
       CALL WRITE_IR(nout,1,idum,6,edum,2)
       !     Total ion radial conduction flux
       label='  Sum of ion conduction fluxes ='
       CALL WRITE_LINE(nout,label,1,0)
       idum(1)=0
       CALL WRITE_IR(nout,1,idum,6,dum,2)
       !  Heat conduction is returned in two forms, add diagonal form
       !  Heat flux consistency check
       label='     *** Heat Flux Consistency Check ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species       qfl_s   dq_s,vq_s,  cp_s,ct_s'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      w/m**2      w/m**2      w/m**2'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          im=jm_s(i)
          iza=IABS(jz_s(i))
          idum(1)=i
          !    Conduction total is sum of components
          rdum(1)=RARRAY_SUM(5,qfl_s(1,i),1)
          !    Diagonal conductivity plus convective velocity
          dq_s(i)=chit_ss(i,i)+chip_ss(i,i)
          vq_s(i)=rdum(1)/(den_iz(im,iza)*z_j7kv*temp_i(im)) &
               &              +dq_s(i)*grt_i(im)/temp_i(im)
          rdum(2)=den_iz(im,iza)*z_j7kv*(-dq_s(i)*grt_i(im) &
               &              +vq_s(i)*temp_i(im))
          !    Conduction total is sum over T' and p' of all species
          rdum(3)=0.0
          DO j=1,m_s
             jm=jm_s(j)
             jza=IABS(jz_s(j))
             rdum(3)=rdum(3)-chit_ss(j,i)*grt_i(jm)/temp_i(jm) &
                  &         -chip_ss(j,i)*grp_iz(jm,jza) &
                  &                      /den_iz(jm,jza)/temp_i(jm)
          ENDDO
          rdum(3)=(rdum(3)+qeb_s(i))*den_iz(im,iza)*temp_i(im)*z_j7kv
          CALL WRITE_IR(nout,1,idum,3,rdum,2)
       ENDDO
       !  Heat conduction, velocity
       label='     *** Heat Conduction, Velocity ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species         X_s         V_s'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      m**2/s         m/s'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          idum(1)=i
          rdum(1)=dq_s(i)
          rdum(2)=vq_s(i)
          CALL WRITE_IR(nout,1,idum,2,rdum,2)
       ENDDO
       !     Effective ion heat conduction, velocity
       label='  Effective ion conduction, velocity ='
       CALL WRITE_LINE(nout,label,1,0)
       idum(1)=0
       CALL RARRAY_ZERO(4,rdum)
       DO i=1,m_s
          im=jm_s(i)
          iz=jz_s(i)
          iza=IABS(iz)
          IF(iz > 0) THEN
             rdum(1)=rdum(1)+(chit_ss(i,i)+chip_ss(i,i))*den_iz(im,iza)
             rdum(2)=rdum(2)+RARRAY_SUM(5,qfl_s(1,i),1)/temp_i(im)/z_j7kv
             rdum(3)=rdum(1)+den_iz(im,iza)
             rdum(4)=rdum(4)+(chit_ss(i,i)+chip_ss(i,i))*den_iz(im,iza) &
                  &                      *grt_i(im)/temp_i(im)
          ENDIF
       ENDDO
       rdum(1)=rdum(1)/rdum(3)
       rdum(2)=(rdum(2)+rdum(4))/rdum(3)
       CALL WRITE_IR(nout,1,idum,2,rdum,2)
       !  Thermal diffusivity matrices
       !    On p'/p           
       label='     *** Thermal Diffusivity Matrices ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species    on dp/dr ... by Species'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      m**2/s ...'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          idum(1)=i
          CALL WRITE_IR(nout,1,idum,m_s,chip_ss(1,i),2)
       ENDDO
       !    On T'/T
       label='     Species    on dT/dr ... by Species'
       CALL WRITE_LINE(nout,label,1,0)
       label='           -      m**2/s ...'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          idum(1)=i
          CALL WRITE_IR(nout,1,idum,m_s,chit_ss(1,i),2)
       ENDDO
       !  Radial energy (conduction+convection) fluxes
       label='     *** Radial Energy (Cond+Conv) Fluxes ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species    ban-plat          PS   classical'// &
            &'       <E.B>  extern src       total'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -      w/m**2      w/m**2      w/m**2'// &
            &'      w/m**2      w/m**2      w/m**2'
       CALL WRITE_LINE(nout,label,0,0)
       CALL RARRAY_ZERO(6,dum)
       DO i=1,m_s
          idum(1)=i
          im=jm_s(i)
          iz=jz_s(i)
          DO k=1,5
             rdum(k)=qfl_s(k,i)+2.5*gfl_s(k,i)*temp_i(im)*z_j7kv
          ENDDO
          rdum(6)=RARRAY_SUM(5,rdum,1)
          CALL WRITE_IR(nout,1,idum,6,rdum,2)           
          IF(iz > 0) THEN
             DO k=1,5
                dum(k)=dum(k)+qfl_s(k,i)+2.5*gfl_s(k,i)*temp_i(im)*z_j7kv
             ENDDO
             dum(6)=dum(6)+rdum(6)
          ENDIF
       ENDDO
       !     Total ion radial energy flux
       label='  Sum of ion energy fluxes ='
       CALL WRITE_LINE(nout,label,1,0)
       idum(1)=0
       CALL WRITE_IR(nout,1,idum,6,dum,2)           
       !  Bootstrap current is returned in two forms
       !  Bootstrap current consistency check
       label='     *** Bootstrap Current Consistency Check ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='      p_bsjb   dp_s,dt_s'
       CALL WRITE_LINE(nout,label,0,0)
       label='    A*T/m**2    A*T/m**2'
       CALL WRITE_LINE(nout,label,0,0)
       !    Total bootstrap current
       rdum(1)=p_bsjb
       !    Bootstrap current summed over T' and p' components
       rdum(2)=0.0
       DO i=1,m_s
          im=jm_s(i)
          iza=IABS(jz_s(i))
          rdum(2)=rdum(2)-bsjbt_s(i)*grt_i(im)/temp_i(im) &
               & -bsjbp_s(i)*grp_iz(im,iza) &
               &               /den_iz(im,iza)/temp_i(im)
       ENDDO
       CALL WRITE_IR(nout,0,idum,2,rdum,2)
       !  Bootstrap current arrays
       !    On p'/p           
       label='     *** Bootstrap Current Arrays ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species    on dp/dr    on dT/dr'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -    A*T/m**2    A*T/m**2'
       CALL WRITE_LINE(nout,label,0,0)
       DO i=1,m_s
          idum(1)=i
          rdum(1)=bsjbp_s(i)
          rdum(2)=bsjbt_s(i)
          CALL WRITE_IR(nout,1,idum,2,rdum,2)
       ENDDO
       !  Current response to external source
       label='     *** Bootstrap and External Source Current ***'           
       CALL WRITE_LINE(nout,label,2,1)
       label='   Bootstrap    External       Total'
       CALL WRITE_LINE(nout,label,0,0)
       label='    A*T/m**2    A*T/m**2    A*T/m**2'
       CALL WRITE_LINE(nout,label,0,0)
       rdum(1)=p_bsjb
       rdum(2)=p_exjb
       rdum(3)=rdum(1)+rdum(2)
       CALL WRITE_IR(nout,0,idum,3,rdum,2)
       !  Flow velocities on outside midplane
       label='     *** Flow Velocities on Outside Midplane ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='     Species       v-tor       v-pol       v-par'// &
            &'      v-perp'
       CALL WRITE_LINE(nout,label,0,0)
       label='           -         m/s         m/s         m/s'// &
            &'         m/s'
       CALL WRITE_LINE(nout,label,0,0)                        
       btor=bt0/(1.0+p_eps)
       bpol=btor/p_fhat
       btot=SQRT(btor**2+bpol**2)*btor/ABS(btor)
       DO i=1,m_s
          idum(1)=i
          im=jm_s(i)
          iz=jz_s(i)
          iza=IABS(iz)
          ppr=p_fhat*grp_iz(im,iza)*z_j7kv &
               & /(z_coulomb*iz*den_iz(im,iza))+p_fhat*p_grphi
          uthai=utheta_s(1,1,i)+utheta_s(1,2,i)+utheta_s(1,3,i)
          !         Toroidal
          rdum(1)=uthai*btor-ppr/btor
          !         Poloidal
          rdum(2)=uthai*bpol
          !         Parallel
          rdum(3)=uthai*btot-ppr/btot
          !         Perpendicular
          rdum(4)=ppr*bpol/btot/btor
          CALL WRITE_IR(nout,1,idum,4,rdum,2)           
       ENDDO
       !  Miscellaneous parameters
       label='     *** Miscellaneous Parameters ***'
       CALL WRITE_LINE(nout,label,2,1)
       label='<J_bs.B>/Bt0 (A/m**2) ='
       rdum(1)=p_bsjb/bt0
       CALL WRITE_LINE_IR(nout,label,0,idum,1,rdum,2)
       label='<J_ex.B>/Bt0 (A/m**2) ='
       rdum(1)=p_exjb/bt0
       CALL WRITE_LINE_IR(nout,label,0,idum,1,rdum,2)
       label='eta parallel (Ohm*m) ='
       rdum(1)=p_etap
       CALL WRITE_LINE_IR(nout,label,0,idum,1,rdum,2)
    ENDIF

    RETURN
  end SUBROUTINE NCLASS_CHECK

end module tx_nclass_mod

