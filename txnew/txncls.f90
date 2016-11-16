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
       &            Q,DPSI,IPSI,EPS,RR,ZI,ZEFF,FT,rlnLei_IN,rlnLii_IN,NUE_IN,NUI_IN, &
       &            BJBS,ETA,ETAS)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Programmed by HONDA Mitsuru (2007/05/24, 2010/04/07 modified)
!     
!     Note that all the parameters are on the GRID point.
!     "rho" denotes the radial coordinate and whether or not it is dimensionless
!        does not matter. It only requires that derivatives are consistent to DPSI.
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
!     DPSI  : partial derivative of poloidal flux to rho [Wb/rho]
!     IPSI  : magnetic flux function, i.e. the product of major radius
!             and toroidal magnetic field (RR*BT) [mT]
!     EPS   : local inverse aspect ratio
!     RR    : major radius [m]
!     ZI    : charge number of bulk ion
!     ZEFF  : effective charge number
!     FT    : trapped particle fraction
!     rlnLei: Coulomb logarithm for electron collisions,   optional
!     rlnLii: Coulomb logarithm for ion collisions,        optional
!     NUE   : Normalized collisionality for electrons,     optional
!     NUI   : Normalized collisionality for ions,          optional
!
!     *** OUTPUT VARIABLES ***
!     BJBS  : bootstrap parallel current <J . B> [AT/m^2]
!     ETA   : neoclassical resistivity [Ohm m]
!     ETAS  : classical resistivity [Ohm m],               optional
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    real(8), intent(in) :: NE,TE,DTE,DPE_IN,NI,TI,DTI,DPI_IN,Q,DPSI,IPSI,EPS,RR, &
    &                      ZI,ZEFF,FT
    real(8), intent(in), optional :: rlnLei_IN,rlnLii_IN,NUE_IN,NUI_IN
    real(8), intent(out) :: BJBS, ETA
    real(8), intent(out), optional :: ETAS

    real(8) :: EPSS,LnLame,NUE,LnLamii,NUI,RPE,ALFA0,ALFA,L31,L32,L34,RNZ,SGMSPTZ,QABS,&
         &     PE,PI,DPE,DPI,TEeV,TIeV
    real(8) :: F31TEFF,F32EETEFF,F32EITEFF,F34TEFF,F33TEFF

    TEeV = TE * 1.D3
    TIeV = TI * 1.D3

    ! On magnetic axis, the neoclassical resistivity reduces to the classical resistivity
    ! and the bootstrap current vanishes.
    IF(EPS == 0.D0) THEN
       RNZ = 0.58D0+0.74D0/(0.76D0+ZEFF)
       LnLame=31.3D0-LOG(SQRT(NE*1.D20)/ABS(TEeV))
       IF(PRESENT(rlnLei_IN)) LnLame = rlnLei_IN

       BJBS = 0.D0
       ETA = 1.D0 / (1.9012D4*TEeV**1.5D0/(ZEFF*RNZ*LnLame))
       IF(PRESENT(ETAS)) ETAS = ETA
       RETURN
    END IF

    EPSS = EPS*SQRT(EPS)
    QABS = ABS(Q)

    PE  = NE * TEeV * 1.D20 * AEE
    PI  = NI * TIeV * 1.D20 * AEE
    DPE = DPE_IN    * 1.D20 * AEE * 1.D3
    DPI = DPI_IN    * 1.D20 * AEE * 1.D3

!     LnLam : coulomb logarithm
!     NU    : collisional frequency [/s]

    LnLame=31.3D0-LOG(SQRT(NE*1.D20)/ABS(TEeV))
    IF(PRESENT(rlnLei_IN)) LnLame = rlnLei_IN
    NUE=6.921D-18*QABS*RR*NE*1.D20*ZEFF*LnLame/(ABS(TEeV)**2*EPSS)
    IF(PRESENT(NUE_IN))   NUE    = NUE_IN

    LnLamii=30.D0-LOG(ZI**3*SQRT(NI*1.D20)/(ABS(TIeV)**1.5D0))
    IF(PRESENT(rlnLii_IN)) LnLamii = rlnLii_IN
    NUI = 4.90D-18*QABS*RR*NI*1.D20*ZI**4*LnLamii/(ABS(TIeV)**2*EPSS)
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

!     *** Bootstrap Current, <B J_BS> ***

    BJBS = -IPSI*PE/DPSI &
         & *( L31*(DPE+DPI)/PE+L32*DTE/TE+L34*ALFA*(1.D0-RPE)/RPE*DTI/TI)

!     *** Spitzer and Neoclassical Resistivity, ETAS, ETA ***

    F33TEFF = FT/(1.D0+(0.55D0-0.1D0*FT)*SQRT(NUE)+0.45D0*(1.D0-FT)*NUE/ZEFF**1.5D0)
    RNZ     = 0.58D0+0.74D0/(0.76D0+ZEFF)
    SGMSPTZ = 1.9012D4*TEeV**1.5D0/(ZEFF*RNZ*LnLame)
    ETA  = 1.D0 / (SGMSPTZ * F33(F33TEFF,ZEFF))
    if(present(etas)) ETAS = 1.D0 / SGMSPTZ

  END SUBROUTINE SAUTER

!     *********************
!     *  Fitting Function *
!     *********************

  pure real(8) FUNCTION F33(X,Z)
    
    real(8), intent(in) :: X,Z

    F33 = 1.D0+(-(1.D0+0.36D0/Z)+(0.59D0-0.23D0*X)*X/Z)*X

  END FUNCTION F33

  pure real(8) FUNCTION F31(X,Z)

    real(8), intent(in) :: X,Z

    F31 = ((1.D0+1.4D0/(Z+1.D0))+(-1.9D0+(0.3D0+0.2D0*X)*X)*X/(Z+1.D0))*X

  END FUNCTION F31

  pure real(8) FUNCTION F32EE(X,Z)

    real(8), intent(in) :: X,Z

    F32EE = ((0.05D0+0.62D0*Z)/(Z*(1.D0+0.44D0*Z))*(1.D0-X**3) &
         &       +( 1.D0/(1.D0+0.22D0*Z)*(1.D0-1.2D0*X+0.2D0*X**2) &
         &         +1.2D0/(1.D0+0.5D0*Z)*X**2)*X)*X

  END FUNCTION F32EE

  pure real(8) FUNCTION F32EI(X,Z)

    real(8), intent(in) :: X,Z

    F32EI =(-(0.56D0+1.93D0*Z)/(Z*(1.D0+0.44D0*Z))*(1.D0-X**3) &
         &       +( 4.95D0/(1.D0+2.48D0*Z)*(1.D0-0.55D0*X-0.45D0*X**2) &
         &         -1.2D0/(1.D0+0.5D0*Z)*X**2)*X)*X

  END FUNCTION F32EI

end module sauter_mod

!------------------------------------------------------------------------

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

  SUBROUTINE TX_NCLASS(NR,ETAout,BJBSout, &
       &               ChiNCpel,ChiNCtel,ChiNCpil,ChiNCtil, &
       &               dTsdrho,dPsdrho,p_gr2phi_in,IER)
!****************************************************************************
!
!  Input : NR,dTsdrho,dPsdrho
!  Output: ETAout,BJBSout,
!          ChiNCpel,ChiNCtel,ChiNCpil,ChiNCtil,IER
!
!****************************************************************************
!
!TX_NCLASS calculates various parameters and arrays for NCLASS.
!Please note that type declarations of all variables except "INTEGER" 
!  in NCLASS subroutine are "REAL(*4)" or "SINGLE" but not "REAL*8" 
!  or "DOUBLE".
!
!The radial coordinate in NCLASS is arbitrary, though it is expressed as 'rho'.
!Considering the compatibility with TASK/TX, we choose the volume as 'rho'.
!
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
!  upar_s(3,m,s)-parallel flow of s from force m (T*m/s)
!                m=1, p', T', Phi'
!                m=2, <E.B>
!                m=3, fex_iz
!  utheta_s(3,m,s)-poloidal flow of s from force m (m/s/T)
!                  m=1, p', T'
!                  m=2, <E.B>
!                  m=3, fex_iz
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
    INCLUDE 'nclass/pamx_mi.inc'
    INCLUDE 'nclass/pamx_ms.inc'
    INCLUDE 'nclass/pamx_mz.inc'
    INCLUDE 'txncls.inc'
    INTEGER(4), INTENT(IN)  :: NR
    INTEGER(4), INTENT(OUT) :: IER
    real(8), intent(in)  :: dTsdrho(NSM),dPsdrho(NSM)
    real(8), intent(in) :: p_gr2phi_in
    REAL(8), INTENT(OUT) :: ETAout, BJBSout
    integer(4) :: i, k_out, k_v, ier_check, j, k, l
    real    :: a0, bt0, e0, p_eps, p_q, q0l, r0, ds
    REAL(8) :: EpsL, PAL, PZL, &
         &     ChiNCpel, ChiNCtel, ChiNCpil, ChiNCtil
    real(8) :: RAL, sgn_xi_s
    real(8), dimension(NSM) :: PTsVL, PNsVL, PsVL
    real(8) :: zp43, Vte, Vti, Wc3, Vc3, Vbc3
!    real(8) :: PZMAX
    real    :: smallvalue = 1.e-5
    real    :: tau_ss(mx_ms,mx_ms) ! Added argument from NCLASS

    !     *** Internal minor radius ***

    RAL = RA

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
    k_potato = 0

    c_den    = 1.E10
    !  *** Potate orbit factors ****************
    c_potb   = REAL(elip(NRA)*fipol(0)/(2.D0*Q(0)**2*rr))
    c_potl   = REAL(Q(0)*RR)
    !  *****************************************

    amu_i(1) = real(amas(1))
    amu_i(2) = real(amas(2))
    IF(Zeff > 1.D0) amu_i(3) = REAL(PAL)

    do i = 1, NSM
       PTsVL(i) = Var(NR,i)%T
       PNsVL(i) = Var(NR,i)%n
       PsVL(i)  = Var(NR,i)%p

       den_iz(i,int(abs(achg(i)))) = REAL(PNsVL(i)) * 1.e20
    end do

!  m_i-number of isotopes (1<mi<mx_mi+1)
!  m_z-highest charge state of all species (0<mz<mx_mz+1)

    m_i = 0
    m_z = 0
    DO i = 1,mx_mi
       ds = 0.0
       DO j = 1,mx_mz
          IF(den_iz(i,j) > c_den) THEN
             ds = den_iz(i,j)
             IF(j > m_z) m_z = j
          ENDIF
       ENDDO
       IF(ds > c_den) m_i = m_i+1
    ENDDO

!!$    m_i      = 2
!!$    IF(Zeff > 1.D0) m_i   = 3
!!$    PZMAX    = achg(2)
!!$    IF(Zeff > 1.D0) PZMAX = PZL
!!$    m_z      = INT(PZMAX)
    
    EpsL  = epst(NR)
    p_b2  = real(bbt(NR)) ! <B^2>
    p_bm2 = real(bit(NR)) ! <B^-2>

    p_eb  = real(BEpara(NR))
    ! No ohmic current causes infinite resistivity in the NCLASS module.
    IF(p_eb == 0.D0) p_eb = 1.e-10

    if( EpsL > 0.d0 ) then
       p_ft    = real(ft(NR))
       p_fhat  = real(fipol(NR)/sdt(NR))
       ! poloidal moments of geometric factor for PS viscosity
       DO i=1,3
          p_fm(i)=REAL(real(i,8)*( (1.D0-SQRT(1.D0-EpsL**2))/EpsL)**(2*i) &
               &                *(1.D0+real(i,8)*SQRT(1.D0-EpsL**2))/((1.D0-EpsL**2)**1.5D0 &
               &                *(Q(NR)*RR)**2))
       ENDDO
    else
       p_ft      = smallvalue ! how to choose this value significantly affect results
       p_fhat    = real(fipol(1)/sdt(1))
       p_fm(1:3) = 0.0
    end if
!!    p_fm(1:3) = 0.0 ! No Pfirsch-Schulter viscosity

    ! p_grbm2 : <|nabla rho|^2/B^2>
    !     Note: One may think that p_grbm2 = (<R^2>/I^2-<B^-2>)(I/dV/drho dpsi/dV)^2.
    !           If the large aspect ratio appoximation is used with fipol=R_0 B_0 assumed,
    !           however, <R^2>/I^2-<B^-2> is identically zero.
    !           Calculating <|nabla rho|^2/B^2> directly in an equilibrium solver enables us
    !           to precisely calculate the value, but p_grbm2 is used in NCLASS solely for
    !           estimating the classical flux. In this sense, for the TASK/TX purpose, 
    !           how p_grbm2 is does not matter.
    p_grbm2   = real(8.d0*(Pi*rr)**4*epsl**2*(2.d0+3.d0*epsl**2)/bbt(0)) ! Large aspect ratio approx.
    IF(ISMTHD == 0) THEN
       p_grphi   = real(dPhiV(NR)) ! [V/m^3]
    ELSE
       IF(NR == 0) THEN
          p_grphi   = real(dPhiV(NR+1)/vro(NR+1)) ! [V/m^3], assump.
       ELSE
          p_grphi   = real(dPhiV(NR)/vro(NR)) ! [V/m^3]
       END IF
    END IF
    ! For orbit squeezing (Houlberg, PoP, 1997, Eq. (B2))
    p_gr2phi  = real(p_gr2phi_in)

    ! p_ngrth : n.grad(Theta) = <B^theta>/<B>
    p_ngrth   = real(bthco(NR)/bbrt(NR)) ! [1/m]
    do i = 1, m_i
       temp_i(i) = real(PTsVL(i))   ! [keV]
       grt_i(i)  = real(dTsdrho(i)) ! dTs/dV [keV/m^3]
       do j = 1, m_z
          grp_iz(i,j) = real(dPsdrho(i)) * 1.e20 ! dPs/dV [keV/m^6]
       end do
    end do
    IF (Zeff > 1.D0) THEN
       temp_i(3) = temp_i(2)
       grt_i(3)  = grt_i(2)
       den_iz(2,INT(achg(2)))  = real((PZL*achg(2)-Zeff)/(achg(2)*(PZL-achg(2)))*PNsVL(2)) * 1.E20
       den_iz(3,INT(PZL))      = real((Zeff-achg(2)**2)/(PZL*(PZL-achg(2)))*PNsVL(2)) * 1.E20
       grp_iz(2,INT(achg(2)))  = real(dPsdrho(2) * (PZL*achg(2)-Zeff)/(achg(2)*(PZL-achg(2)))) * 1.e20
       grp_iz(3,INT(PZL))      = real(dPsdrho(2) * (Zeff-achg(2)**2)/(PZL*(PZL-achg(2)))) * 1.e20
    END IF

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
         &      chip_ss,chit_ss,dp_ss,dt_ss,iflag,tau_ss)

    IF(k_out >0 .and.k_out <= 7) THEN
       p_eps = REAL(epst(NR))
       p_q   = REAL(Q(NR))
       r0    = REAL(RR)
       a0    = REAL(RAL)
       e0    = REAL(1.D0)
       bt0   = REAL(fipol(0)/rr)
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

!    Check particle flux
!    if(NR/=0.and.NT==500) write(6,*) rho(NR),(sum(gfl_s(1:2,1))+sum(gfl_s(4:5,1)))/vro(NR) ! w/o classical flux
!    if(NR/=0.and.NT==500) write(6,*) rho(NR),sum(gfl_s(1:5,1))/vro(NR) ! total
!    if(NR/=0.and.NT==500) write(6,*) rho(NR),gfl_s(4,1)/vro(NR) ! Ware
!    if(NR/=0.and.NT==500) write(6,'(F10.7,1P5E15.7)') rho(NR),(gfl_s(i,1)/vro(NR),i=1,5) ! breakdown
!    write(6,*) rho(NR),gfl_s(4,1)*sdt(NR)*1.d-20,gfl_s(4,2)*sdt(NR)*1.d-20 ! Ware flux
!    write(6,*) rho(NR),sum(gfl_s(1:5,1))*sdt(NR)*1.d-20,sum(gfl_s(1:5,2))*sdt(NR)*1.d-20 ! Particle flux

    !     *** Parameters to pass ***

    IF(k_potato == 0) THEN
       IF(iflag /= -1 .and. k_out > 0) THEN
          WRITE(6,'(A,I4,2X,A,I4)') "XX iflag=",iflag,"NR=",NR
          IER=1
          RETURN
       ENDIF
    ELSE
       IF(iflag /= 0 .and. k_out > 0) WRITE(6,'(A,I4,2X,A,I4)') "XX iflag=",iflag,"NR=",NR
    ENDIF

    !---   Bootstrap current
    BJBSout =-sum(  real(bsjbt_s(1:NSM),8) *(dTsdrho(1:NSM) / PTsVL(1:NSM)) &
         &        + real(bsjbp_s(1:NSM),8) *(dPsdrho(1:NSM) / PsVL (1:NSM)))
    IF(Zeff > 1.D0) THEN
       BJBSout =-( real(bsjbt_s(1),8) *(dTsdrho(1) / PTsVL(1)) &
              &  + real(bsjbp_s(1),8) *(dPsdrho(1) / PsVL(1)) &
              &  + real(bsjbt_s(2),8) *(dTsdrho(2) / PTsVL(2)) &
              &  + real(bsjbp_s(2),8) *(dPsdrho(2) / PsVL(2) * (PZL*achg(2)-Zeff)/(achg(2)*(PZL-achg(2)))) &
              &  + real(bsjbt_s(3),8) *(dTsdrho(2) / PTsVL(2)) &
              &  + real(bsjbp_s(3),8) *(dPsdrho(2) / PsVL(2) * (Zeff-achg(2)**2)/(PZL*(PZL-achg(2)))))
    END IF
    IF(k_potato == 0 .and. NR == 0) BJBSout = 0.D0

    !---   Neoclassical resistivity
    ETAout = real(p_etap,8)

    !---   Neoclassical thermal diffusivity (Diagonal effect only)
    !     chi*_ss = <|\nabla \rho|^2>*chi_NC, where chi_NC is in unit of [m^2/s].
    !     chi*_ss = <|\nabla V|^2>*chi_NC if rho = V is assumed in NCLASS.
    if( NR /= 0 ) then
       ChiNCpel = ChiNC * real(chip_ss(1,1),8) / sst(NR)
       ChiNCtel = ChiNC * real(chit_ss(1,1),8) / sst(NR)
       ChiNCpil = ChiNC * real(chip_ss(2,2),8) / sst(NR)
       ChiNCtil = ChiNC * real(chit_ss(2,2),8) / sst(NR)
    else
       ChiNCpel = 0.d0
       ChiNCtel = 0.d0
       ChiNCpil = 0.d0
       ChiNCtil = 0.d0
    end if
!!$    if(ChiNCpel < 0.d0) ChiNCpel = 0.d0
!!$    if(ChiNCtel < 0.d0) ChiNCtel = 0.d0
!!$    if(ChiNCpil < 0.d0) ChiNCpil = 0.d0
!!$    if(ChiNCtil < 0.d0) ChiNCtil = 0.d0

    !---   Poloidal flows by NCLASS for graphics: uthhat(moment, from what, species)
    do i = 1, NSM
       UsthNCL(NR,i,1) = real(utheta_s(1,1,i),8) ! (particle flow, p' T')
       UsthNCL(NR,i,2) = real(utheta_s(1,2,i),8) ! (particle flow, <E.B>)
    end do

!    write(6,*) rho(nr),utheta_s(1,1,2)-upar_s(1,1,2)/p_b2,-BVsdiag(NR,2,1)/p_b2
!    write(6,*) rho(nr),utheta_s(1,1,2)+utheta_s(1,2,2),Var(NR,2)%Uthhat
!    write(6,*) rho(nr),upar_s(1,1,2)+upar_s(1,2,2),Var(NR,2)%BUpar

    !---   Neoclassical viscosities

!    if(NR /= 0) then
       do i = 1, NSM ! species
          do j = 1, k_order
             do k = 1, k_order
                xmu(NR,i,j,k) = FSNC * real(ymu_s(j,k,i),8)
             end do
          end do
       end do
!    else
!       xmu(NR,:,:,:) = 0.d0
!    end if

    !---   Friction coefficients normalized by m_a n_a : lab(NR,i,j,k,l)/(m_a n_a)

    do i = 1, NSM ! species
       do j = 1, NSM ! species
          do l = 1, k_order
             do k = 1, k_order
                if( l == k ) then ! l + k = even
                   sgn_xi_s = xi_s(i)
                else              ! l + k = odd
                   sgn_xi_s =-xi_s(i)
                end if
                if( i == j ) then
                   lab(NR,i,j,k,l) = sgn_xi_s * ( calm_i(k,l,jm_s(i)) + xi_s(j) * caln_ii(k,l,i,j) )
                else
                   lab(NR,i,j,k,l) = sgn_xi_s * (                       xi_s(j) * caln_ii(k,l,i,j) )
                end if
                ! normalization: divided by m_i n_i
                lab(NR,i,j,k,l) = lab(NR,i,j,k,l) / (amas(i) * amp * Var(NR,i)%n * 1.D20)
!                write(6,*) i,k,k,l,jm_s(i),real(lab(NR,i,j,k,l))
             end do
          end do
       end do
    end do
!!$    !---   Sanity check
!!$    write(6,*) lab(NR,1,1,1,1)*tau_ss(1,1)/(amas(1)*amp*Var(nr,1)%n*1.d20),-zeff ! l_11^ee
!!$    write(6,*) lab(NR,1,2,1,1)*tau_ss(1,1)/(amas(1)*amp*Var(nr,1)%n*1.d20),Var(nr,2)%n/Var(nr,1)%n ! l_11^ei
!!$    write(6,*) lab(NR,1,1,1,2),lab(NR,1,1,2,1),1.5d0*lab(NR,1,1,1,1) ! l_12^ee
!!$    write(6,*) lab(NR,1,2,1,2),0.d0 ! l_12^ei
!!$    write(6,*) lab(NR,1,2,2,1),1.5d0*lab(NR,1,2,1,1) ! l_21^ei
!!$    write(6,*) lab(NR,1,1,2,2)*tau_ss(1,1)/(amas(1)*amp*Var(nr,1)%n*1.d20),-(sqrt(2.d0)+13.d0/4.d0*zeff) ! l_22^ee
!!$    write(6,*) lab(NR,1,2,2,2),0.d0 ! l_22^ei
!!$    stop
!!$    !---   Correction for exact symmetry of friction coefficients
!!$    lab(NR,2,1,1,1) = lab(NR,1,2,1,1)*(amas(1)*Var(NR,1)%n)/(amas(2)*Var(NR,2)%n)
!!$    lab(NR,2,1,2,1) = lab(NR,1,2,1,2)*(amas(1)*Var(NR,1)%n)/(amas(2)*Var(NR,2)%n)
!!$    lab(NR,2,1,2,2) = lab(NR,1,2,2,2)*(amas(1)*Var(NR,1)%n)/(amas(2)*Var(NR,2)%n)

    !---   Beam neoclassical viscosity

    do j = 1, k_order
       xmuf(NR,j) = FSNCB * 0.d0
    end do

    !---   Beam friction coefficients
    !     laf : thermal and fast,   normalized by m_a n_a
    !     lfb : fast and thermal,   normalized by m_e n_e
    !     lff : fast and fast,      normalized by m_e n_e

    laf(NR,1:NSM,1:2,1:2) = 0.d0
    lfb(NR,1:NSM,1:2,1:2) = 0.d0
    lff(NR,1:2,1:2) = 0.d0

    zp43 = 0.75d0 * sqrt(pi)
    Vte = sqrt(2.d0 * abs(Var(NR,1)%T) * rKilo / (amas(1) * amqp))
    Vti = sqrt(2.d0 * abs(Var(NR,2)%T) * rKilo / (amas(2) * amqp))
    Wc3 = zp43 * Vte*Vte*Vte * amas(1)/amb * (achg(2)*achg(2) * Var(NR,2)%n / Var(NR,1)%n)
    Vc3 = Wc3 * (achg(2)*achg(2) * Var(NR,2)%n / amas(2)) / (achg(2)*achg(2) * Var(NR,2)%n) * amb
    Vbc3 = Vb*Vb*Vb+Vc3
    !   --- electron and beam ---
    laf(NR,1,1,1)   = achgb*achgb * PNbV(NR) / Var(NR,1)%n / tau_ss(1,1)        ! l_11^ef
    lab(NR,1,1,1,1) = lab(NR,1,1,1,1) - laf(NR,1,1,1)                          ! l_11^ee corrected
    lfb(NR,1,1,1)   = laf(NR,1,1,1)                                            ! l_11^fe
    laf(NR,1,2,1)   = 1.5d0 * laf(NR,1,1,1)                                    ! l_21^ef
    lab(NR,1,1,2,1) = lab(NR,1,1,2,1) - laf(NR,1,2,1)                          ! l_21^ee corrected
    !   --- ion and beam ---
    if( Vc3 /= 0.d0 ) then                                                     ! l_11^if
       laf(NR,2,1,1) = zp43*(achgb*achgb*PNbV(NR)/(achg(2)*achg(2)*Var(NR,2)%n)) &
            & *(1.d0+amas(2)/amb)/log(Vbc3/Vc3)*(Vti*Vti*Vti)/Vc3 / tau_ss(2,2)
    else
       laf(NR,2,1,1) = 0.d0
    end if
    lfb(NR,2,1,1)   = laf(NR,2,1,1) * (amas(2)*Var(NR,2)%n*tau_ss(1,1) &
         &                           /(amas(1)*Var(NR,1)%n*tau_ss(2,2)))         ! l_11^fi
    lab(NR,2,2,1,1) = lab(NR,2,2,1,1) - laf(NR,2,1,1)                           ! l_11^ii corrected
    !   --- beam and beam ---
    lff(NR,1,1) = - lfb(NR,1,1,1) - lfb(NR,2,1,1)                               ! l_11^ff

    !---   Parallel heat flows: 2<B q_s,para>/(5p_s)

    do j = 1, 2
       do i = 1, NSM
          UsparNCL(NR,i,j) = real(sum(upar_s(j,1:2,i)),8)
       end do
    end do
    
    !---   <B.nabla.Pi> for viscous heating

    do i = 1, NSM
       BnablaPi(NR,i) = p_b2 * sum(ymu_s(1,1:2,i)*(utheta_s(1:2,1,i)+utheta_s(1:2,2,i)))
    end do

    !---   Particle flux : < Gamma . nabla psi >

    do i = 1, NSM
       gflux(NR,i,3) = sum(gfl_s(1:5,i)) * sdt(NR)
    end do

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
