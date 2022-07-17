module tx_variables
  implicit none
  public

contains

!***************************************************************
!
!        Calculate variables
!
!***************************************************************

  subroutine TXCALV(XL,ID)

    use tx_commons
    use tx_interface, only : dfdx, replace_interpolate_value
    use tx_core_module, only : intg_vol_p
    use libitp, only: aitken2p
!    use aux_system, only : Vbabsmax
    real(8), DIMENSION(0:NRMAX,1:NQMAX), INTENT(IN) :: XL
    integer(4), intent(in), optional :: ID
    integer(4) :: NR, i, JSMTHD, j
    real(8), parameter :: fourPisq = 4.d0 * Pi * Pi
    real(8) :: sdtvac, dPsitVdVvac, sum_int, fricsum(NSM)!, BBL

    JSMTHD = ISMTHD / 10

    if(present(ID)) then
      if( ID == 0 .or. ID == 1 ) then
         ! The pres0 and ErV0 are the values evaluated at the previous time step
         !   for numerical stability when using turbulent transport models.
         if(MDFIXT == 0) then
            pres0(:) = (  XL(:,LQe5) + XL(:,LQi5) + XL(:,LQz5)) * 1.d20 * rKeV
         else
            pres0(:) = (  XL(:,LQe1) * XL(:,LQe5) &
                 &      + XL(:,LQi1) * XL(:,LQi5) &
                 &      + XL(:,LQz1) * XL(:,LQz5)) * 1.d20 * rKeV
         end if
         ErV0 (0)    =   0.d0
         dPhidV(:)   = dfdx( vv,XL(:,LQm1),NRMAX,0)           ! dPhi/dV
         dPhidrho(:) = dfdx(rho,XL(:,LQm1),NRMAX,0,daxs=0.d0) ! dPhi/drho
         if(ISMTHD == 0) then
            ErV0(1:NRMAX) = - sst(1:NRMAX) / vro(1:NRMAX) / drhodr(1:NRMAX) * dPhidV(1:NRMAX)
         else
            if(JSMTHD == 1) call replace_interpolate_value(dPhidrho(1),1,rho,dPhidrho)
            ErV0(1:NRMAX) = - sst(1:NRMAX) / vro(1:NRMAX)**2 / drhodr(1:NRMAX) * dPhidrho(1:NRMAX)
         end if
         if(ID /= 0) then
            return
         else
            !     Only right after restart file load

            !     Averaged beam injection energy
            Eb =  (esps(1) + 0.5d0 * esps(2) + esps(3) / 3.d0) * Ebmax
         end if
      end if
    end if

    PhiV(:)     = XL(:,LQm1)

    PsitdotV(:) = XL(:,LQm2)

    PsidotV(:)  = XL(:,LQm3)

    ! Etor: <E_t> = <1/R>dpsi/dt
    Etor(:) = PsidotV(:) * ait(:)

    PsiV(:) = XL(:,LQm4) * rMUb2
    ! sdt: dpsi/dV
    !    sdt(NRMAX) =2.d0*Pi*rMU0*rIp*1.d6/ckt(NRMAX) is a boundary condition.
    sdtvac = 2.d0 * Pi * rMU0 * rIp * 1.d6 / ckt(NRMAX)
    sdt(:) = dfdx(vv,PsiV,NRMAX,0,dbnd=sdtvac)

    ErV(0)      = 0.d0
    dPhidV(:)   = dfdx( vv,PhiV,NRMAX,0)           ! dPhi/dV
    dPhidrho(:) = dfdx(rho,PhiV,NRMAX,0,daxs=0.d0) ! dPhi/drho
    if(ISMTHD == 0) then
       ErV(1:NRMAX) = - sst(1:NRMAX) / vro(1:NRMAX) / drhodr(1:NRMAX) * dPhidV(1:NRMAX)
    else
       ! Replace dPhidrho(1) by the interpolated value
       if(JSMTHD == 1) call replace_interpolate_value(dPhidrho(1),1,rho,dPhidrho)
       ErV(1:NRMAX) = - sst(1:NRMAX) / vro(1:NRMAX)**2 / drhodr(1:NRMAX) * dPhidrho(1:NRMAX)
    end if

!    BthV (:) =   fourPisq * rho(:)  / drhodr(:) * sdt(:) ! NOT FSA quantity
    BthV (:) =   sqrt(ckt(:)) * sdt(:)
    PsitV(:) =   XL(:,LQm5) * rMU0
!    ! bbt = <B^2>
!    bbt  (:) =   fourPisq*fourPisq / aat(:) * hdt(:)*hdt(:) &
!         &     + ckt(:) * sdt(:)*sdt(:)

    ! Kinetic profiles
    !  -- electrons
    Var(:,1)%n = XL(:,LQe1)
    if(MDFIXT == 0) then
       Var(:,1)%p = XL(:,LQe5)
       Var(:,1)%T = XL(:,LQe5) / Var(:,1)%n
    else
       Var(:,1)%p = XL(:,LQe5) * Var(:,1)%n
       Var(:,1)%T = XL(:,LQe5)
    end if
    !  -- ions
    Var(:,2)%n = XL(:,LQi1)
    if(MDFIXT == 0) then
       Var(:,2)%p = XL(:,LQi5)
       Var(:,2)%T = XL(:,LQi5) / Var(:,2)%n
    else
       Var(:,2)%p = XL(:,LQi5) * Var(:,2)%n
       Var(:,2)%T = XL(:,LQi5)
    end if
    !  -- impurities
    Var(:,3)%n = XL(:,LQz1)
    if(MDFIXT == 0) then
       Var(:,3)%p = XL(:,LQz5)
       Var(:,3)%T = XL(:,LQz5) / Var(:,3)%n
    else
       Var(:,3)%p = XL(:,LQz5) * Var(:,3)%n
       Var(:,3)%T = XL(:,LQz5)
    end if

    ! Radial flows
    Var(:,1)%UrV = XL(:,LQe2) / Var(:,1)%n
    Var(:,2)%UrV = XL(:,LQi2) / Var(:,2)%n
    Var(:,3)%UrV = XL(:,LQz2) / Var(:,3)%n
    do i = 1, NSM
       Var(0,i)%Ur       = 0.d0
       Var(1:NRMAX,i)%Ur = Var(1:NRMAX,i)%UrV / suft(1:NRMAX) ! <u.nabla V>/<|nabla V|>
    end do

    ! Parallel flows
    Var(:,1)%BUpar = XL(:,LQe3)
    Var(:,2)%BUpar = XL(:,LQi3)
    Var(:,3)%BUpar = XL(:,LQz3)
    Var(:,1)%Bqpar = XL(:,LQe6)
    Var(:,2)%Bqpar = XL(:,LQi6)
    Var(:,3)%Bqpar = XL(:,LQz6)

    ! Toroidal flows
    Var(:,1)%RUph = XL(:,LQe4) / Var(:,1)%n
    Var(:,2)%RUph = XL(:,LQi4) / Var(:,2)%n
    Var(:,3)%RUph = XL(:,LQz4) / Var(:,3)%n
    Var(:,1)%UphR = XL(:,LQe7) / Var(:,1)%n
    Var(:,2)%UphR = XL(:,LQi7) / Var(:,2)%n
    Var(:,3)%UphR = XL(:,LQz7) / Var(:,3)%n
    Var(:,1)%Uph  = Var(:,1)%UphR / ait(:) ! Uph = <Uph/R>/<1/R>
    Var(:,2)%Uph  = Var(:,2)%UphR / ait(:)
    Var(:,3)%Uph  = Var(:,3)%UphR / ait(:)

    ! Diamagnetic particle flows
    Var(:,1)%BV1 = XL(:,LQe8)
    Var(:,2)%BV1 = XL(:,LQi8)
    Var(:,3)%BV1 = XL(:,LQz8)

    ! Poloidal current function: RB_t
    !    fipol(NRMAX) = rbvt is a boundary condition.
    dPsitVdVvac = rbvt / fourPisq * aat(NRMAX)

    fipol(:) = fourPisq / aat(:) * dfdx(vv,PsitV,NRMAX,0,dbnd=dPsitVdVvac)
    BphV (:) = fipol(:) / rr ! NOT FSA quantity
    ! hdt: dpsit/dV
    hdt  (:) =   fipol(:) * aat(:) / fourPisq
    ! BEpol: <B_p E_p> = -<R^-2> I dPsit/dt
    BEpol(:) = - fipol(:) * aat(:) * PsitdotV(:)

    ! --- Beam ions ---
    PNbV (:) =   XL(:,LQb1)
    PTbV (:) = 2.d0 / 3.d0 * Eb ! Mean temperature of beam ions
    PNbVinv(:) = 0.d0
    forall (NR = 0:NRMAX, PNbV(NR) /= 0.d0) PNbVinv(NR) = 1.d0 / PNbV(NR)
    do NR = 0, NRMAX
       UbrVV(NR)   = XL(NR,LQb2) * PNbVinv(NR) * mod(MDBEAM,2)
       BUbparV(NR) = XL(NR,LQb3) * PNbVinv(NR)
       RUbphV(NR)  = XL(NR,LQb4) * PNbVinv(NR)
       UbphVR(NR)  = XL(NR,LQb7) * PNbVinv(NR)
       BVbdiag(NR) = XL(NR,LQb8) * PNbVinv(NR)
!       UbphV(NR) = BUbparV(NR) * fipol(NR) / bbt(NR) / rr
       UbphV(NR)   = UbphVR(NR) / ait(NR)
    end do
    UbrV(0) = 0.d0
    UbrV(1:NRMAX) = UbrVV(1:NRMAX) / suft(1:NRMAX) ! <u.nabla V>/<|nabla V|>

    if(MDBEAM /= 1 .and. abs(FSRP) == 0.d0) then ! No beam ions in the SOL
       sum_int = 0.d0
       do NR = NRA, NRMAX
          sum_int = sum_int + intg_vol_p(PNbV,nr)
          PNbV(NR)    = 0.d0
          PNbVinv(NR) = 0.d0
       end do
       PNbV(NRA-1) = PNbV(NRA-1) + sum_int / (vv(NRA) - vv(NRA-1))
       if(PNbV(NRA-1) /= 0.d0) PNbVinv(NRA-1) = 1.d0 / PNbV(NRA-1)
       UbrVV(NRA:NRMAX)   = 0.d0
       BUbparV(NRA:NRMAX) = 0.d0
       UbphV(NRA:NRMAX)   = 0.d0
       RUbphV(NRA:NRMAX)  = 0.d0
       UbphVR(NRA:NRMAX)  = 0.d0
    end if

    PNbRPV(:)=   XL(:,LQr1)
    ! -----------------

    ! Effective charge
    do NR = 0, NRMAX
       Zeff(NR) = ( sum(achg(2:)**2 * Var(NR,2:)%n) + achgb**2 * (PNbV(NR) + PNbRPV(NR)) ) &
            &   / Var(NR,1)%n
    end do

    PN01V(:) =   XL(:,LQn1)
    PN02V(:) =   XL(:,LQn2)
    PN03V(:) =   XL(:,LQn3)

    PN0zV(:) =   XL(:,LQnz)

    if(ISMTHD == 0) then
       do i = 1, NSM
          dTsdpsi(:,i) = dfdx(PsiV,Var(:,i)%T,NRMAX,0) ! dT/dpsi
          dPsdpsi(:,i) = dfdx(PsiV,Var(:,i)%p,NRMAX,0) ! dp/dpsi
       end do
    else
       do i = 1, NSM
          dTsdpsi(:,i) = dfdx(rho,Var(:,i)%T,NRMAX,0,daxs=0.d0) ! dT/drho
          dPsdpsi(:,i) = dfdx(rho,Var(:,i)%p,NRMAX,0,daxs=0.d0) ! dp/drho
          ! Replace dPsdpsi(1) by the interpolated value
          if(JSMTHD == 1) then
             call replace_interpolate_value(dTsdpsi(1,i),1,rho,dTsdpsi(:,i))
             call replace_interpolate_value(dPsdpsi(1,i),1,rho,dPsdpsi(:,i))
          end if
          dTsdpsi(1:NRMAX,i) = dTsdpsi(1:NRMAX,i) / (sdt(1:NRMAX) * vro(1:NRMAX))! dTs/dpsi
          dTsdpsi(0,i) = AITKEN2P(vv(0),dTsdpsi(1,i),dTsdpsi(2,i),dTsdpsi(3,i) &
               &                 ,vv(1),vv(2),vv(3))
          dPsdpsi(1:NRMAX,i) = dPsdpsi(1:NRMAX,i) / (sdt(1:NRMAX) * vro(1:NRMAX))! dps/dpsi
          dPsdpsi(0,i) = AITKEN2P(vv(0),dPsdpsi(1,i),dPsdpsi(2,i),dPsdpsi(3,i) &
               &                 ,vv(1),vv(2),vv(3))
       end do
   end if

    ! Safety factor: q = dpsit/dpsi = I<R^-2>/(4 Pi^2)*dV/dpsi
    !   Note: The latter definition is preferable because the components I(=fipol) and
    !         dV/dpsi(=1/sdt) have been already computed with constraints by B.C.
    Q(:) = fipol(:) * aat(:) / ( fourPisq * sdt(:) )

    ! Square of the poloidal magnetic field: <B_p^2> = ckt * (dpsi/dV)^2
    Bpsq(:)  = ckt(:) * sdt(:)*sdt(:)

    do NR = 1, NRMAX
       ! Grid velocity
       UgV(NR) = - FSUG * 0.5d0 * PsitdotV(NR) * vro(NR) / ( rho(NR) * PsitV(NRA) )
       ! qhat square: q^^2 = I^2/(2<B_p^2>)(<1/R^2>-1/<R^2>)
       qhatsq(NR) = 0.5d0 * fipol(NR)**2 / Bpsq(NR) * (aat(NR) - 1.d0 / rrt(NR))
       ! Metric coefficient: <B^2><R^2>-I^2
!       bri(NR) = bbt(NR) * rrt(NR) - fipol(NR)**2 ! => causing instability through LQe2, LQi2 and LQz2
       bri(NR) = rrt(NR) * Bpsq(NR) * (1.d0 + 2.d0 * qhatsq(NR))
!       if(abs(bri(NR)) < 1.d-10) bri(NR) = 0.d0
    end do
    NR = 0
    bri(NR) = 0.d0
    UgV(NR) = 0.d0
    if(ieqread == 0) then ! Perfect cylinder (Pfirsch-Schluter contribution nil)
       qhatsq(NR) = 0.d0
    else
       qhatsq(NR) = Q(NR)**2
    end if

    do NR = 0, NRMAX
       ! Coefficients regarding Pfirsch-Schluter contribution: 2q^^2/(1+2q^^2)
       Fqhatsq(NR) = 2.d0 * qhatsq(NR) / (1.d0 + 2.d0 * qhatsq(NR))
       ! Parallel Electric field: <BE_//>
       BEpara(NR) = fipol(NR) * aat(NR) * (- PsitdotV(NR) / Q(NR) + PsidotV(NR))
       ! <B^2>
       bbt(NR) = fourPisq**2 * hdt(NR)**2 / aat(NR) + ckt(NR) * sdt(NR)**2
       ! <B^theta> = 4 pi^2 dpsi/dV
       bthco(NR) = fourPisq * sdt(NR)
    end do

    !  Diamagnetic particle and heat flows: B V_1s, B V_2s
    forall (i = 1:NSM, NR = 0:NRMAX) 
       ! Diamag. particle flow
       !       BVsdiag(NR,i,1) = - fipol(NR) &
       !            & * ( dPsdpsi(NR,i) * rKilo / ( achg(i) * Var(NR,i)%n ) + dPhidpsiL(NR) )
       BVsdiag(NR,i,1) = Var(NR,i)%BV1
       ! Diamag. heat flow
       BVsdiag(NR,i,2) = - fipol(NR) / achg(i) * dTsdpsi(NR,i) * rKilo

       ! Var%Uthhat = U.nabla theta/B.nabla theta [m/s/T]
       !    Var%Uthhat = ( <B U_{s//}> - BV_{1s} ) / <B^2> or
       !               = ( <R U_{s zeta}> - <R^2>/I BV_{1s} ) / I, alternatively
       Var(NR,i)%Uthhat = ( Var(NR,i)%BUpar - Var(NR,i)%BV1 ) / bbt(NR)
       ! Poloidal rotation for graphics: u_theta = uthhatV * BthV [m/s]
       Var(NR,i)%Uth = Var(NR,i)%Uthhat * BthV(NR)
    end forall

!!$    ! For JSPF meeting 2020
!!$    do NR = 0, NRMAX
!!$       fricsum(:) = 0.d0
!!$       do i = 1, NSM
!!$          do j = 1, NSM
!!$             fricsum(i) = fricsum(i) + amas(i)*amp*Var(NR,i)%n*1.d20 &
!!$                  & *( lab(NR,i,j,1,1)*Var(NR,j)%RUph &
!!$                  &   -lab(NR,i,j,1,2)*( fipol(NR)/bbt(NR)*Var(NR,j)%Bqpar &
!!$                  &                     +  bri(NR)/bbt(NR)/fipol(NR)*BVsdiag(NR,j,2)))
!!$          end do
!!$       end do
!!$
!!$       write(221,'(2F8.5,8ES15.7)') t_tx,rho(NR) &
!!$            & , -sum(achg(:)*Var(NR,:)%n)*aee*1.d20*PsidotV(NR) &
!!$            & , -sum(fricsum(:))
!!$    end do

  end subroutine TXCALV

!***************************************************************
!
!   Calculate coefficients
!
!***************************************************************

  subroutine TXCALC(IC)

    use tx_commons
    use tx_interface, only : dfdx, txmmm95, &
         &                   moving_average, coulog, CORR, coll_freq
    use tx_core_module, only : sub_intg_vol
    use tx_nclass_mod
    use sauter_mod
    use aux_system
    use matrix_inversion, only : tx_matrix_inversion
    use tx_ripple
    use cdbm_mod
    use mod_eqneo, only : wrap_eqneo
    use mod_cross_section, only : RateCoef_C_RC_A6, RateCoef_C_IZ_S5, Fraction_C0toC5
    use libitp, only: aitken2p
!    use tx_ntv, only : NTVcalc, rNuNTV, UastNC

    integer(4), intent(in) :: IC
    integer(4), save :: NRB = 1
    integer(4) :: NR, NR1, IER, i, MDANOMabs, model_cdbm, izvpch, MDLNEOL
    integer(4) :: NHFM ! miki_m 10-08-06
    real(8) :: Sigma0, Vte, Vti, Vtz, Vtb, Wte, Wti, Wtz, EpsL, &
         &     rNuAsI_inv, BBL, Va, Wpe2, PN0tot, &
         &     PROFML, PROFCL, Dturb, DeL, &
         &     Cs, Lc, RhoIT, ExpArg, AiP, DISTAN, UbparaL, &
         &     rNuOLL, SiLCL, SiLCBL, SiLCphL, RL, DBW, PTiVA, &
         &     Chicl, factor_bohm, rNustar, &
         &     RLOSS, sqz, rNuDL, Ln, LT, etai_chk, kthrhos, &
         &     RhoSOL, V0ave, Viave, DturbA, rLmean, Sitot, &
         &     rGCIM, rGIM, rHIM, OMEGAPR !09/06/17~ miki_m
    real(8), dimension(0:NRMAX) :: gr2phi
    real(8), dimension(1:NHFMmx) :: EpsLM 
!    real(8) :: rLmeanL, QL, rNuAsE_inv 
    real(8), save :: Fcoef = 1.d0
    real(8) :: PTiVav, N02INT, RatSCX, sum1, sum2, zvpch
    real(8) :: Frdc, Dcoef
    real(8) :: omegaer, omegaere, omegaeri
    real(8) :: zrNuee, zxi, znusteinv, zLamEZ, zCRZ, zsigma0, zsigmaneo
    real(8) :: rhoni, dvexbdr ! CDBM
    real(8) :: xb, fp, fdp, gfun
    real(8), dimension(0:NRMAX) :: pres, ddPhidpsi, tmp
    ! Mainly for derivatives
    real(8), dimension(:), allocatable :: dErdr, dpdr, dErdrS, ErVlc
    real(8), dimension(:), allocatable :: dQdrho, dlnNedrhov
    real(8), dimension(:,:), allocatable :: dTsdV, dPsdV, dNsdV, dNsdrho, dTsdrho

    MDANOMabs = abs(MDANOM)
    MDLNEOL   = mod(MDLNEO,10)
    izvpch = 0

    !     *** Constants ***

    !     Neutral cross section (Characteristic atomic collision cross section)
    !     (NRL Plasma Formulary p53 Eq. (1) (2019))

    Sigma0 = 8.8d-21

    ! ************* Auxiliary heating part *************

    call txauxs

    ! ************** For turbulent transport **************
    !   The reason that we keep the pressure throughout iteration is to stabilize numerical
    !      instability caused by the CDBM model, strongly dependent on the pressure gradient.
    !   In order to suppress oscillation of the pressure in the direction of time, 
    !      we take the average between pres and pres0, evaluated at the previous time step,
    !      when differentiating the pressure with respect to vv.
    !   In addition, during iteration, pres is fixed.
    !   This holds true with the radial electric field, ErV, if one prefers.
    pres(:) = sum(Var(:,:)%p, dim=2) * 1.d20 * rKeV
    select case(iprestab)
    case(1,3)
       do NR = 0, NRMAX
          pres0(NR)  = sum(PNsV_FIX(NR,:)*PTsV_FIX(NR,:)) * 1.d20 * rKeV
       end do
    case(2)
       ! pres0 estimated in txcalv is used; hence doing nothing here.
    case default
       pres0(:)  = pres(:)
    end select
    if(MDANOM > 0 .and. maxval(FSANOM) > 0.d0) pres(:)  = 0.5d0 * (pres(:) + pres0(:))

    allocate(ErVlc, mold=ErV0)
    ErVlc(:) = 0.5d0 * (ErV_FIX(:) + ErV0(:))

    if(PROFM == 0.d0 .and. FSDFIX(2) /= 0.d0) then
       PROFML = (Var(NRA,1)%T * rKeV / (16.d0 * AEE * sqrt(bbt(NRA)))) / FSDFIX(2)
    else
       PROFML = PROFM
    end if

    if(PROFC == 0.d0 .and. FSDFIX(3) /= 0.d0) then
       PROFCL = (Var(NRA,1)%T * rKeV / (16.d0 * AEE * sqrt(bbt(NRA)))) / FSDFIX(3)
    else
       PROFCL = PROFC
    end if

    ! ************** Turbulent transport end **************

    ! Banana width
!    Wbane = (Q(0) * sqrt(RR * amas(1) * Var(0,1)%T * rKeV * amqp) / BphV(0))**(2.d0/3.d0)
!    Wbani = (Q(0) * sqrt(RR * amas(2) * Var(0,2)%T * rKeV * amqp) / (achg(2) * BphV(0)))**(2.d0/3.d0)

    ! Banana width at separatrix
!    PTiVA = Var(NRA,2)%T
    PTiVA = 0.5d0 * Var(0,2)%T
    DBW = 3.d0 * sqrt(PTiVA * rKeV * (amas(2)*amp)) * Q(NRA) / (achg(2) * AEE * BphV(NRA)) &
         & / sqrt(ra / RR)

    !     *** Calculate derivatives in advance ***
    !     !!! Caution !!!
    !        The r-derivatives of variables, or near-variables (ex. temperature) should be
    !          estimated by their vv-derivatives multiplied by dV/dr because they are
    !          evaluated on the vv-abscissa. On the other hand, those of the other
    !          parameters (ex. radial electric field, poloidal magnetic field) should be
    !          directly calculated.

    allocate(dErdr, dErdrS, dpdr, mold=array_init_NR)
    dErdr(:) = dfdx(rpt,ErVlc,NRMAX,0)
    dpdr (:) = dfdx(rpt,pres ,NRMAX,0)
!    dpdr (:) = vro(:) / ravl * dfdx(vv ,pres ,NRMAX,0)

    allocate(dQdrho, dlnNedrhov, mold=vv)
    allocate(dTsdV, dPsdV, dNsdV, mold=array_init_NRNS)
    do i = 1, NSM
       dTsdV(:,i) = dfdx(vv,Var(:,i)%T,NRMAX,0) ! [keV/m^3]
       dPsdV(:,i) = dfdx(vv,Var(:,i)%p,NRMAX,0) ! [10^{20}keV/m^6]
       dNsdV(:,i) = (dPsdV(:,i) - Var(:,i)%n * dTsdV(:,i)) / Var(:,i)%T ! [10^{20}/m^6]
    end do
    
    dQdrho(:) = vro(:) * dfdx(vv,Q,NRMAX,0)
    ! (1/ne)dne/drho_v
    dlnNedrhov(0) = 0.d0
    do NR = 1, NRMAX
       dlnNedrhov(NR) = 2.d0 * vlt(NR) / ( rhov(NR) * Var(NR,1)%n ) * dNsdV(NR,1)
    end do

!!D02    write(6,'(F8.5,I4,2F11.6)') T_TX,NRB,Rho(NRB),PT02V(NR)

    !  Smoothing Er gradient for numerical stability
    do NR = 0, NRMAX
       dErdrS(NR) = moving_average(NR,dErdr,NRMAX,NRA)
    end do

    ! *** Temperatures for neutrals ***

    PT01V(:) = 0.5d0 * amas(2) * V0**2 * amqp / rKilo

    !  --- For thermal neutrals originating from slow neutrals ---

    if(IC == 1) then
       ! SCX : source of PN02V
       do NR = 0, NRMAX
          SCX(NR) = Var(NR,2)%n * SiVcx(Var(NR,2)%T) * PN01V(NR) * 1.D40
       end do

       ! NRB : Boundary of N02 source
       NRB = NRA
       do NR = NRMAX - 1, 0, -1
          RatSCX = SCX(NR) / SCX(NRMAX)
          if(RatSCX < 1.D-8) then
             NRB = NR
             exit
          else if (RatSCX < tiny_cap) then
             NRB = NR - 1
          end if
       end do

       ! Neoclassical quantities associated with equilibrium 
       call wrap_eqneo
    end if

    ! PTiVav : Particle (or density) averaged temperature across the N02 source region
    sum1 = 0.d0 ; sum2 = 0.d0
    do nr = nrb, nrmax
       call sub_intg_vol(PN02V,NR,N02int,NR)
       sum1 = sum1 + N02int
       sum2 = sum2 + N02int * (0.5d0 * (Var(NR-1,2)%T + Var(NR,2)%T))
    end do
    if(sum1 > epsilon(1.d0)) then
       PTiVav = sum2 / sum1
    else
       PTiVav = Var(NRA,2)%T
    end if
    do nr = 0, nrmax
       PT02V(NR) = PTiVav
    end do

    PT03V(:) =   Var(:,2)%T

    ! *********************************

    ! Set flag for CDBM model
    model_cdbm = 0
    select case(int(FSCBSH))
    case(1)
       model_cdbm = model_cdbm + 2
    case(2)
       model_cdbm = model_cdbm + 4
    case(3)
       model_cdbm = model_cdbm + 6
    end select
    if( FSCBEL /= 0.d0 ) model_cdbm = model_cdbm + 1

    ! Calculate CDIM coefficient
    RAQPR(:) = vro(:) * drhodr(:) * dfdx(vv, rho**4 / Q, NRMAX , 0) ! cf. txcalv.f90 L501

    ! Orbit squeezing effect for neoclassical solvers
    !   gr2phi = psi'(Phi'/psi')' = (dpsi/dV dV/drho)^2 d/dpsi(dPhi/dpsi), 
    !      where rho is an arbitrary radial coordinate.
    !                      ^^^^^^^^^
    !   gr2phi = (dpsi/dV)^2 d/dpsi(dPhi/dpsi) if rho = V is assumed.
    if( mod(MDOSQZ,10) /= 1 .or. IC == 1 ) THEN
       ddPhidpsi(0) = 0.d0 ! d/dpsi(dPhi/dpsi)
       do NR = 1, NRMAX-1
          ddPhidpsi(NR) = (  ( PhiV(NR+1) - PhiV(NR) ) / ( PsiV(NR+1) - PsiV(NR) ) &
               &     - ( PhiV(NR) - PhiV(NR-1) ) / ( PsiV(NR) - PsiV(NR-1) ) ) &
               &     / ( PsiV(NR+1) - PsiV(NR-1) ) * 2.d0
       end do
       ! Phi(NRMAX+1)=Phi(NRMAX) is assumed.
       NR = NRMAX
       ddPhidpsi(NR) = (- ( PhiV(NR) - PhiV(NR-1) ) / ( PsiV(NR) - PsiV(NR-1) ) ) &
            &     / ( PsiV(NR) - PsiV(NR-1) )
       
       if(MDOSQZ > 10) then
          do NR = 0, NRMAX
             tmp(NR) = moving_average(NR,ddPhidpsi,NRMAX)
          end do
          ddPhidpsi(:) = tmp(:)
       end if
       ! For NCLASS
       gr2phi(:) = sdt(:)*sdt(:) * ddPhidpsi(:) * MDOSQZN
    end if

    !     *** ExB shearing rate (Hahm & Burrel PoP 1995) ***
    ! omega_ExB = |- <R^2B_p^2>/B d/dpsi(dPhi/dpsi)|
    wexb(:)  = abs(- sst(:) * sdt(:) / BB * ddPhidpsi(:))

    !  Coefficients

    L_NR: do NR = 0, NRMAX

       Vte = sqrt(2.d0 * abs(Var(NR,1)%T) * rKilo / (amas(1) * amqp))
       Vti = sqrt(2.d0 * abs(Var(NR,2)%T) * rKilo / (amas(2) * amqp))
       Vtz = sqrt(2.d0 * abs(Var(NR,3)%T) * rKilo / (amas(3) * amqp))
       Vtb = sqrt(2.d0 * abs(Var(NR,2)%T) * rKilo / (amb     * amqp)) ! ??

       ! For neutrals
       PN0tot = (PN01V(NR) + PN02V(NR) + PN03V(NR)) * 1.d20
       SiVizA(NR)  = SiViz(Var(NR,1)%T,Var(NR,2)%T/amas(2))
       SiVcxA(NR)  = SiVcx(Var(NR,2)%T)
       ! Maxwellian rate coefficient of effective recombination
       SiVa6A(NR)  = RateCoef_C_RC_A6(Var(NR,1)%T,Var(NR,1)%n)
       ! Maxwellian rate coefficient of effective ionization
       SiVsefA(NR) = RateCoef_C_IZ_S5(Var(NR,1)%T,Var(NR,1)%n) &
            &      * Fraction_C0toC5(Var(NR,1)%T,Var(NR,1)%n)

       !     *** Ionization rate ***

!old       !     (NRL Plasma Formulary p55 Eq. (12) (2019))
!old       XXX = MAX(Var(NR,1)%T * 1.D3 / EION, 1.D-2)
!old       SiV = 1.D-11 * sqrt(XXX) * exp(- 1.d0 / XXX) &
!old            &              / (EION**1.5d0 * (6.d0 + XXX))
!old       rNuION(NR) = FSION * SiV * (PN01V(NR) + PN02V(NR)) * 1.d20
       rNuION(NR) = FSION * SiVizA(NR) * PN0tot

       !     *** Slow neutral diffusion coefficient ***
       !  For example,
       !    E.L. Vold et al., NF 32 (1992) 1433

       !  Maxwellian velocity for slow neutrals
       V0ave = sqrt(4.d0 * V0*V0 / Pi)

       !  Total Maxwellian rate coefficients
!!$       if(nr <= NRB) then
!!$          Sitot = (SiVcxA(NRB) * Var(NRB,2)%n + SiVizA(NRB) * Var(NRB,1)%n) *1.d20
!!$       else
          Sitot = (SiVcxA(NR) * Var(NR,2)%n + SiVizA(NR) * Var(NR,1)%n) *1.d20
!!$       end if

       !  Diffusion coefficient for slow neutrals (short m.f.p.)
!old       D01(NR) = FSD01 * V0**2 &
!old            &   / (Sigma0 * (PN01V(NR) * V0  + (Var(NR,2)%n + PN02V(NR)) &
!old            &      * Vti) * 1.d20)
!       D01(NR) = FSD01 * V0ave**2 &
!            &   / (3.d0 * SiVcxA(NR) * Var(NR,2)%n * 1.d20)
       D01(NR) = FSD01 * V0ave*V0ave / (3.d0 * Sitot)

       !     *** Thermal neutral diffusion coefficient ***

       !  Maxwellian thermal velocity at the separatrix
!       Viave = sqrt(8.d0 * Var(NR,2)%T * rKilo / (Pi * amas(2) * amqp))
       Viave = sqrt(8.d0 * PT02V(NR) * rKilo / (Pi * amas(2) * amqp))

       !  Mean free path for fast neutrals
       rLmean = Viave / Sitot
!!D02       !  Locally determined mean free path for fast neutrals
!!D02       rLmeanL = sqrt(8.d0 * Var(NR,2)%T * rKilo / (Pi * amas(2) * amqp)) &
!!D02            &  / ((SiVcxA(NR) * Var(NR,2)%n + SiVizA(NR) * Var(NR,1)%n) *1.d20)

       !  Diffusion coefficient for fast neutrals (short to long m.f.p.)
!old       D02(NR) = FSD02 * Vti**2 &
!old            &   / (Sigma0 * (Var(NR,2)%n + PN01V(NR) + PN02V(NR)) &
!old            &      * Vti * 1.d20)
!       D02(NR) = FSD02 * Viave**2 &
!            &   / (3.d0 * SiVcxA(NR) * Var(NR,2)%n * 1.d20)
!!       if(rho(nr) < 0.9d0) then
!!          D02(NR) = D02(NR) * (-0.33333d0*rho(nr)+0.35d0)
!!       else if (rho(nr) >= 0.9d0 .and. rho(nr) <= 1.0d0) then
!!          D02(NR) = D02(NR) * (3.5055d0*rho(nr)**3-2.5055d0)
!!       end if
       D02(NR) = FSD02 * Viave*Viave / (3.d0 * Sitot)
!!!!!!!!!!       if(nr >= NRB .and. nr <= NRA) D02(NR) = D02(NR) * reduce_D02(NR,rLmean)
!!D02       if(nt >= ntmax-1) write(6,'(I3,F10.6,F11.3,3F10.6,2ES11.3)') nr,r(nr),d02(nr),rLmean,rLmeanL,Var(NR,2)%T,SCX(NR)

       !     *** Halo neutral diffusion coefficient ***

       Viave = sqrt(8.d0 * Var(NR,2)%T * rKilo / (Pi * amas(2) * amqp))
       D03(NR) = FSD03 * Viave*Viave / (3.d0 * Sitot)

       !     *** Impurity neutral diffusion coefficient (slow) ***

       !  Maxwellian velocity for slow impurity neutrals
       V0ave = sqrt(4.d0 * V0*V0 / Pi)

       !  Total Maxwellian rate coefficients
       Sitot = (SiVsefA(NR) + SiVa6A(NR)) * Var(NR,1)%n *1.d20

       !  Diffusion coefficient for slow impurity neutrals (short m.f.p.)
       D0z(NR) = FSD0z * V0ave*V0ave / (3.d0 * Sitot)

       !     *** Charge exchange rate ***
       !  For thermal ions (assuming that energy of deuterium
       !                    is equivalent to that of proton)

!old       !     (Riviere, NF 11 (1971) 363, Eq.(4))
!old       Scxi = 6.937D-19 * (1.d0 - 0.155d0 * LOG10(Var(NR,2)%T*1.D3))**2 &
!old            & / (1.d0 + 0.1112D-14 * (Var(NR,2)%T*1.D3)**3.3d0) ! in m^2
!old       Vave = sqrt(8.d0 * Var(NR,2)%T * rKilo / (PI * amas(2) * amqp))
!old       rNuiCX(NR) = FSCX * Scxi * Vave * (PN01V(NR) + PN02V(NR)) * 1.d20
       rNuiCX(NR)  = FSCX * SiVcxA(NR) * PN0tot
       !  For thermal loss by charge exchange
       rNuiCXT(NR) = FSCX * SiVcxA(NR) * PN01V(NR) * 1.d20

       !  For beam ions
!old       !     (C.E. Singer et al., Comp. Phys. Comm. 49 (1988) 275, p.323, B163)
!oldold       Vave = sqrt(8.d0 * Eb * rKilo / (PI * amas(2) * amqp))
!oldold       rNubCX(NR) = FSCX * Scxb * Vave * (PN01V(NR) + PN02V(NR)) * 1.d20
!old       rNubCX(NR) = FSCX * Scxb * Vb * PN0tot
       rNubCX(NR) = FSCX * ratecxb(NR) * PN0tot

       !     *** Collision frequency (with neutral) ***

       rNu0s(NR,1) = FSNCOL * PN0tot * Sigma0 * Vte
       rNu0s(NR,2) = FSNCOL * PN0tot * Sigma0 * Vti
       rNu0s(NR,3) = FSNCOL * PN0tot * Sigma0 * Vtz
       rNu0b(NR)   = FSNCOL * PN0tot * Sigma0 * Vtb

       !     *** Collision frequency (90 degree deflection) ***
       !         (see  Helander & Sigmar textbook (2002), p.5, p.79,
       !               Hirshman & Sigmar, p.1105, around (4.7)      )
       !     c.f. Braginskii's collision time: rNue = rNuei, rNui = rNuii / sqrt(2)

       rNuei(NR) = coll_freq(NR,1,2)
       rNuii(NR) = coll_freq(NR,2,2)
       rNuiz(NR) = coll_freq(NR,2,3)

       !     *** Collision frequency (energy equipartition) ***
       !     (D. V. Sivukhin, Rev. Plasma Phys. Vol. 4, Consultants Bureau (1966))
       !         Energy relaxation time between two kinds of particles
!approximate formula (amas(1) << amas(2))      rNuTei(NR) = rNuei(NR) * (2.d0 * amas(1) / amas(2))
       rNuTei(NR) = rNuTeq(NR,1,2) ! e-i
       rNuTez(NR) = rNuTeq(NR,1,3) ! e-z
       rNuTiz(NR) = rNuTeq(NR,2,3) ! i-z

       BBL = sqrt(bbt(NR))

       !     *** Transit frequency ***
       !     Omega_{Ta} = v_{Ta} / L_c^*
       !     (see Hirshman & Sigmar, p.1115, below (4.65))
       Wte = Vte / (Q(NR) * RR) ! Omega_te; transit frequency for electrons
       Wti = Vti / (Q(NR) * RR) ! Omega_ti; transit frequency for ions
       Wtz = Vtz / (Q(NR) * RR) ! Omega_tz; transit frequency for impurities
       
       !     *** Collisionality parameter ***
       !     nu_{*a} = 1 / (Omega_{Ta} * tau_a * eps**1.5)
       !     (see Hirshman & Sigmar, p.1112, above (4.45))
       if(NR /= 0) then
!!$          rNuAss(NR,1) = sqrt(2.d0) * rNuei(NR) / ( Wte * epst(NR)**1.5d0 )
!!$          rNuAss(NR,2) = sqrt(2.d0) * rNuii(NR) / ( Wti * epst(NR)**1.5d0 )
          rNuAss(NR,1) = rNuei(NR) / ( Wte * epst(NR)**1.5d0 )
          rNuAss(NR,2) = rNuii(NR) / ( Wti * epst(NR)**1.5d0 )
          rNuAss(NR,3) = rNuiz(NR) / ( Wtz * epst(NR)**1.5d0 )
       end if

       if( MDLNEOL == 3 .or. MDLNEO > 10 ) then
          !     *** Neoclassical coefficients by NCLASS ***
          !     (W. A. Houlberg, et al., Phys. Plasmas 4 (1997) 3230)

          call TX_NCLASS(NR,ETAvar(NR,2),BJBSvar(NR,2), &
               &         ChiNCp(NR,:),ChiNCt(NR,:), &
               &         dTsdV(NR,:),dPsdV(NR,:),gr2phi(NR),IER)
          if(IER /= 0) IERR = IER
!          write(6,'(A,F8.5,6ES15.7)') 'NCLASS',rho(NR),chincp(nr,:),chinct(nr,:)
       end if
       if( MDLNEOL <= 2 .or. MDLNEO > 10 ) then
          !     *** Neoclassical coefficients by Matrix Inversion ***
          !     (M. Kikuchi and M. Azumi, Plasma Phys. Control. Fusion 37 (1995) 1215)

          call tx_matrix_inversion(NR,ETAvar(NR,1),BJBSvar(NR,1), &
               &         ChiNCp(NR,:),ChiNCt(NR,:), &
               &         ddPhidpsi(NR)*MDOSQZN,MDLNEOL)
!          write(6,'(A,F8.5,6ES15.7)') 'MATINV',rho(NR),chincp(nr,:),chinct(nr,:)
!!$          !     Resistivity should be the classical (Spitzer) one at axis.
!!$          if( NR == 0 ) ETAvar(NR,1) = CORR(1.d0) * amas(1) * amqp * rNuei(NR) &
!!$               &                     / (Var(NR,1)%n * 1.d20 * AEE)
       end if

       !     *** Helical neoclassical viscosity ***

!       if(abs(FSHL) > 0.d0 .and. NR > 0) then
       if(abs(FSHL) > 0.d0 ) then
!          IF(int(FSHL) .EQ. 1 ) THEN !! for FSHL = 1 (single Fourier mode)
!!$            Wte = Vte * NCph / RR
!!$            Wti = Vti * NCph / RR
!!$            EpsL = EpsH * rho(NR)**2
!!$            rNuAsE_inv = EpsL**1.5d0 * Wte / (sqrt(2.d0) * rNuei(NR))
!!$            rNuAsI_inv = EpsL**1.5d0 * Wti / (sqrt(2.d0) * rNuii(NR))
!!$            IF(NR.EQ.0) THEN
!!$               BLinv=0.d0
!!$               omegaer=0.d0
!!$            ELSE
!!$!             QL=(Q0-QA)*(1.d0-rho(NR)**2)+QA
!!$!             Bthl = BB*rpt(NR)/(QL*RR)
!!$!             BLinv=BB/Bthl
!!$               BBL=sqrt(BphV(NR)**2 + BthV(NR)**2)
!!$               BLinv=BBL/BthV(NR)
!!$               omegaer=ErVlc(NR)/(BBL*rpt(NR))
!!$            ENDIF
!!$            omegaere=EpsL*rpt(NR) / RR * omegaer**2 / rNuei(NR)**2
!!$!            rNueHL(NR) = FSHL * Wte * BLinv * rNuAsE_inv &
!!$            rNueHL(NR) =        Wte * BLinv * rNuAsE_inv &
!!$            &            /(3.d0+1.67d0*omegaere)
!!$
!!$            omegaeri=EpsL*rpt(NR) / RR * omegaer**2 / rNuii(NR)**2
!!$!            rNuiHL(NR) = FSHL * Wti * BLinv * rNuAsI_inv &
!!$            rNuiHL(NR) =        Wti * BLinv * rNuAsI_inv &
!!$            &            /(3.d0+1.67d0*omegaeri)
!!$
!!$!          UHth=(RR/NCph)/sqrt((RR/NCph)**2+(rpt(NR)/NCth)**2)
!!$!          UHph=(rpt(NR)/NCth)/sqrt((RR/NCph)**2+(rpt(NR)/NCth)**2)
!!$!          UHth  = real(NCth,8) / NCph
!!$!          UHph  = 1.d0
!!$!    09/02/11 mm
!!$
!!$            UHth=(RR*NCth)/sqrt((RR*NCth)**2+(rpt(NR)*NCph)**2)
!!$!---- 09/11/26 AF
!!$!          UHph=-(rpt(NR)*Ncph)/sqrt((RR*NCth)**2+(rpt(NR)*NCph)**2)
!!$            UHph=-Ncph/sqrt((RR*NCth)**2+(rpt(NR)*NCph)**2)
!!$
!!$            rNueHLthth(NR)=UHth*UHth*rNueHL(NR) ! [s^-1]
!!$            rNueHLthph(NR)=UHth*UHph*rNueHL(NR) ! [m^-1 s^-1]
!!$            rNueHLphth(NR)=UHth*UHph*rNueHL(NR) ! [m^-1 s^-1]
!!$            rNueHLphph(NR)=UHph*UHph*rNueHL(NR) ! [m^-2 s^-1]
!!$            rNuiHLthth(NR)=UHth*UHth*rNuiHL(NR)
!!$            rNuiHLthph(NR)=UHth*UHph*rNuiHL(NR)
!!$            rNuiHLphth(NR)=UHth*UHph*rNuiHL(NR)
!!$            rNuiHLphph(NR)=UHph*UHph*rNuiHL(NR)
!!$
!          ELSE IF (abs(FSHL) .GE. 2.d0) THEN  !! for FSHL = 2 (multiple Fourier modes)
!!!kokokara
          if(NR == 0) then
             omegaer=0.d0
          else
             omegaer=ErVlc(NR)/(BBL*rpt(NR))
          end if
!         
          do NHFM = 1, NHFMmx 
             if (HPN(NHFM,1) == 0 .and. HPN(NHFM,2) == 0) then
                rNueHLththM(NHFM,NR) = 0.d0
                rNueHLthphM(NHFM,NR) = 0.d0
                rNueHLphthM(NHFM,NR) = 0.d0
                rNueHLphphM(NHFM,NR) = 0.d0
                rNuiHLththM(NHFM,NR) = 0.d0
                rNuiHLthphM(NHFM,NR) = 0.d0
                rNuiHLphthM(NHFM,NR) = 0.d0
                rNuiHLphphM(NHFM,NR) = 0.d0
                cycle
             endif

             EpsLM(NHFM) = EpsHM(NHFM,0)              + EpsHM(NHFM,1) * RHO(NR)        &
                  &      + EpsHM(NHFM,2) * RHO(NR)**2 + EpsHM(NHFM,3) * RHO(NR)**3
!
             omegaere   = abs(EpsLM(NHFM)) * rpt(NR) / RR * omegaer**2 / rNuei(NR)**2
             rNueHLM(NHFM, NR) = FSHL * abs(EpsLM(NHFM))**1.5d0 * Var(NR,1)%T * rKilo          &
                  &            / (amas(1) * amqp * RR**2 * rNuei(NR) * (3.d0 + 1.67d0 * omegaere)) ! [s^-1]
!
             omegaeri   = abs(EpsLM(NHFM)) * rpt(NR) / RR * omegaer**2 / rNuii(NR)**2
             rNuiHLM(NHFM, NR) = FSHL * abs(EpsLM(NHFM))**1.5d0 * Var(NR,2)%T * rKilo          &
                  &            / (amas(2) * amqp * RR**2 * rNuii(NR) * (3.d0 + 1.67d0 * omegaeri)) ! [s^-1]
!
             if (UHphSwitch == 0) then    ! For the case which (m=0, n>0) component DOES NOT exist
                UHth =    RR * HPN(NHFM,1) / sqrt((RR * HPN(NHFM,1))**2 + (rpt(NR) * HPN(NHFM,2))**2) ! [nondimensional]
                UHph =         HPN(NHFM,2) / sqrt((RR * HPN(NHFM,1))**2 + (rpt(NR) * HPN(NHFM,2))**2) ! [m^-1]
             else if (HPN(NHFM,1) == 0 .and. HPN(NHFM,2) > 0) then ! to avoid NaN and Infty for (m=0, n>0) component
                UHth = 0.d0 ! [nondimensional]
                UHph = 1.d0 ! [nondimensional]
             else   ! For the case which (m=0, n>0) component exist
                UHth =      RR * HPN(NHFM,1) / sqrt((RR * HPN(NHFM,1))**2 + (rpt(NR) * HPN(NHFM,2))**2) ! [nondimensional]
                UHph = rpt(NR) * HPN(NHFM,2) / sqrt((RR * HPN(NHFM,1))**2 + (rpt(NR) * HPN(NHFM,2))**2) ! [nondimensional]
             endif

             ! Dimension of thph, phth, phph will change according to the value of UHphSwitch
             rNueHLththM(NHFM,NR)=UHth*UHth*rNueHLM(NHFM,NR)
             rNueHLthphM(NHFM,NR)=UHth*UHph*rNueHLM(NHFM,NR)
             rNueHLphthM(NHFM,NR)=UHth*UHph*rNueHLM(NHFM,NR)
             rNueHLphphM(NHFM,NR)=UHph*UHph*rNueHLM(NHFM,NR)
             rNuiHLththM(NHFM,NR)=UHth*UHth*rNuiHLM(NHFM,NR)
             rNuiHLthphM(NHFM,NR)=UHth*UHph*rNuiHLM(NHFM,NR)
             rNuiHLphthM(NHFM,NR)=UHth*UHph*rNuiHLM(NHFM,NR)
             rNuiHLphphM(NHFM,NR)=UHph*UHph*rNuiHLM(NHFM,NR)
             
          end do
          rNueHLthth(NR) = sum(rNueHLththM(1:NHFMmx,NR))
          rNueHLthph(NR) = sum(rNueHLthphM(1:NHFMmx,NR))
          rNueHLphth(NR) = sum(rNueHLphthM(1:NHFMmx,NR))
          rNueHLphph(NR) = sum(rNueHLphphM(1:NHFMmx,NR))
          rNuiHLthth(NR) = sum(rNuiHLththM(1:NHFMmx,NR))
          rNuiHLthph(NR) = sum(rNuiHLthphM(1:NHFMmx,NR))
          rNuiHLphth(NR) = sum(rNuiHLphthM(1:NHFMmx,NR))
          rNuiHLphph(NR) = sum(rNuiHLphphM(1:NHFMmx,NR))
!          write(*,*) ' rNueHLthth=', rNueHLthth(NR)
!          write(*,*) ' rNueHLthph=', rNueHLthph(NR)
!          write(*,*) ' rNueHLphth=', rNueHLphth(NR)
!          write(*,*) ' rNueHLphph=', rNueHLphph(NR)
!          write(*,*) ' rNuiHLthth=', rNuiHLthth(NR)
!          write(*,*) ' rNuiHLthph=', rNuiHLthph(NR)
!          write(*,*) ' rNuiHLphth=', rNuiHLphth(NR)
!          write(*,*) ' rNuiHLphph=', rNuiHLphph(NR)
!!!kokomade 10-08-06
       else
          rNueHLthth(:) = 0.d0
          rNueHLthph(:) = 0.d0
          rNueHLphth(:) = 0.d0
          rNueHLphph(:) = 0.d0
          rNuiHLthth(:) = 0.d0
          rNuiHLthph(:) = 0.d0
          rNuiHLphth(:) = 0.d0
          rNuiHLphph(:) = 0.d0
       end if
!!$       if (ic == 0 .or. ic == 1) then
!!$          if (nr == 0 .or. nr == 1) then
!!$             write(*,*) ' rNueHLthth(',NR,')=', rNueHLthth(NR)
!!$             write(*,*) ' rNueHLthph(',NR,')=', rNueHLthph(NR)
!!$             write(*,*) ' rNueHLphth(',NR,')=', rNueHLphth(NR)
!!$             write(*,*) ' rNueHLphph(',NR,')=', rNueHLphph(NR)
!!$             write(*,*) ' rNuiHLthth(',NR,')=', rNuiHLthth(NR)
!!$             write(*,*) ' rNuiHLthph(',NR,')=', rNuiHLthph(NR)
!!$             write(*,*) ' rNuiHLphth(',NR,')=', rNuiHLphth(NR)
!!$             write(*,*) ' rNuiHLphph(',NR,')=', rNuiHLphph(NR)
!!$          endif
!!$       endif

       !  Magnetic shear, normalized pressure gradient
       S(NR) = rho(NR) / Q(NR) * dQdrho(NR)
       Alpha(NR) = - Q(NR)*Q(NR) * RR * dpdr(NR) * 2.d0 * rMU0 / bbt(NR)

       !   ***** Heat pinch *****

       Vhps(NR,1:NSM) = 0.d0

       !   ***** Thermal diffusivity *****

       !   *** CDBM model ***

       if (maxval(FSANOM) > 0.d0) then

          !   *** CDBM model ***
          select case(MDANOMabs)
          case(1)
             rhoni = Var(NR,2)%n * 1.d20 * (amas(2)*amp)
             dvexbdr = dErdrS(NR) / bbrt(NR)
!             dvexbdr = dErdr(NR) / bbrt(NR)

             call cdbm(BBL,rr,rpt(NR),elip(NR),Q(NR),S(NR),Var(NR,1)%n*1.d20,rhoni,dpdr(NR), &
                  &    dvexbdr,FSCBAL,FSCBKP,rG1,model_cdbm,Dturb, &
                  &    FCDBM(NR),rKappa(NR),rG1h2(NR),wexb(NR))

          !   *** CDIM model ***
          case(2)
             ! Alfven velocity
             Va = sqrt(BBL*BBL / (rMU0 * Var(NR,2)%n * 1.d20 * (amas(2)*amp)))
             ! Squared plasma frequency
             Wpe2 = Var(NR,1)%n * 1.d20 * AEE / (amas(1) * amqp * EPS0)

             ! Arbitrary coefficient for CDIM model
             rGCIM = 10.d0
             OMEGAPR = (ra / rr)**2 * (real(NCph,8) / NCth) * RAQPR(NR)
         
             if(NR == 0) then  ! for s=0
                FCDIM(NR) = 0
             else
                FCDIM(NR) = 3.d0 * (0.5d0 * OMEGAPR)**1.5d0 * (rr / ra)**1.5d0 / (Q(NR) * S(NR)**2)
             end if

             ! ExB rotational shear
             if(NR == 0) then
                rGIM = 0.d0
                rHIM = 0.d0
             else
!                rGIM = rG1
                rGIM = 1.04d0 * rpt(NR)**2 / ( dpdr(NR) * 2.d0 * rMU0 / bbt(NR) &
                     &                         * OMEGAPR * ra**2 * rr**2 )
                rHIM = ra * sqrt( rMU0 * (amas(2)*amp) * Var(NR,2)%n * 1.d20 ) / BthV(NR) * dErdrS(NR) / BBL
             end if

             ! Turbulence suppression by ExB shear for CDIM mode
             rG1h2IM(NR) = 1.d0 / (1.d0 + rGIM * rHIM**2)
             ! Turbulent transport coefficient calculated by CDIM model
             Dturb = rGCIM * FCDIM(NR) * rG1h2IM(NR) * abs(Alpha(NR))**1.5d0 &
                  &              * VC*VC / Wpe2 * Va / (Q(NR) * RR)

!-----------------memo:beta'=dpdr(NR) * 2.d0 * rMU0 / bbt(NR)-----------

!          IF(NR == 0 .OR. NR == 1) & 
!                         write(6,*) ',NR=',NR,'RAQPR=',RAQPR(NR),'OMEGAPR=',OMEGAPR, &
!               &     'Q=',Q(NR),'S=',S(NR),'FCDIM=',FCDIM(NR),'Dturb=',Dturb

          end select
          if(Rho(NR) == 1.d0) DturbA = Dturb 
       else
          rG1h2(NR)   = 0.d0
          FCDBM(NR)   = 0.d0
          rG1h2IM(NR) = 0.d0
          FCDIM(NR)   = 0.d0
          Dturb       = 0.d0
       end if

       !     *** Turbulent transport of particles ***
       !     ***     Wave-particle interaction    ***

       RhoSOL = 1.d0
       if (RHO(NR) < RhoSOL) then
!          DeL = diff_prof(RHO(NR),FSDFIX(1),PROFD,PROFD1,PROFDB) + FSANOM(1) * Dturb
          DeL = diff_prof(RHO(NR),FSDFIX(1),PROFD,PROFD1,PROFDB,PROFD2,0.99d0,0.07d0) &
!          DeL = diff_prof(RHO(NR),FSDFIX(1),PROFD,PROFD1,PROFDB,PROFD2,0.99d0,0.125d0) &
               & + FSANOM(1) * Dturb
       else
          if(FSPCL(1) == 0.d0) then
             ! Bohm-like diffusivity
             factor_bohm = (FSDFIX(1) * PROFD + PROFDB + FSANOM(1) * Dturb) &
                  &  / (Var(NRA,1)%T * rKeV / (16.d0 * AEE * sqrt(bbt(NRA))))
             DeL = factor_bohm * Var(NR,1)%T * rKeV / (16.d0 * AEE * BBL)
          else
             if(FSANOM(1) == 0.d0) then
                ! Fixed value and fixed profile
                DeL = FSPCL(1) * diff_prof(RhoSOL,FSDFIX(1),PROFD,PROFD1,PROFDB,0.d0)
             else
                ! Theory-based anomalous diffusivity
                DeL = FSANOM(1) * DturbA
             end if
          end if
!!$          DeL = FSDFIX(1) * PROFD + FSANOM(1) * Dturb
       end if
       ! Particle diffusivity
!!$       if(rho(nr) > 0.7d0) then
!!$          DeL = DeL * (-5.d0/3.d0*Rho(NR)+13.d0/6.d0)
!!$!          DeL = DeL * (50.d0/9.d0*(Rho(NR)-1.d0)**2+0.5d0)
!!$       end if
       Dfs(NR,:) = Dfs0(:) * DeL

       ! Turbulent pinch term [1/Wb]
       VWpch(NR) = VWpch0 * 0.d0

       !     *** Turbulent transport of momentum ***

       RhoSOL = 1.d0

       if (RHO(NR) < RhoSOL) then
          DeL = diff_prof(RHO(NR),FSDFIX(2),PROFML,PROFM1,PROFMB,0.d0) + FSANOM(2) * Dturb
!pedestal          if(rho(nr) > 0.9d0) DeL = DeL * exp(-120.d0*(rho(nr)-0.9d0)**2)
       else
          if(FSPCL(2) == 0.d0) then
             factor_bohm = (FSDFIX(2) * PROFML + PROFMB + FSANOM(2) * Dturb) &
                  &  / (Var(NRA,1)%T * rKeV / (16.d0 * AEE * sqrt(bbt(NRA))))
!bohm_model2             DeL =  (1.d0 - MOD(FSBOHM,2.d0)) * FSDFIX(2) * PROFML &
!bohm_model2                  &+ FSBOHM * Var(NR,1)%T * rKeV / (16.d0 * AEE * BBL)
             DeL = factor_bohm * Var(NR,1)%T * rKeV / (16.d0 * AEE * BBL)
          else
             if(FSANOM(2) == 0.d0) then
                DeL = FSPCL(2) * diff_prof(RhoSOL,FSDFIX(2),PROFML,PROFM1,PROFMB,0.d0)
             else
                DeL = FSANOM(2) * DturbA
             end if
!pedestal             DeL = FSPCL(2) * FSDFIX(2) * PROFML * exp(-120.d0*(rho(nra)-0.9d0)**2)
          end if
       end if
       ! Viscosity
       rMus(NR,:) = rMus0(:) * DeL

       !   ***** Momentum pinch *****
       !   -R Vmpi/rMui = 1.1 R/Lne + 1.0
       !   [T. Tala et al., 2012, Proc of 24th IAEA FEC (San Diego) ITR/P1-1]

       if(NR /= 0) then
!          zvpch = - rMus(NR,2) / rr * (1.1d0 * rr * moving_average(NR,dlnNedrhov,NRMAX) + 1.d0) ! smoothing 1/Lne
          zvpch = rMus(NR,2) / rr * (1.1d0 * rr * moving_average(NR,dlnNedrhov,NRMAX) + 1.d0) ! smoothing 1/Lne
!          ! +++ reducing Vpch near the separatrix +++
!          zvpch = 0.5d0 * ( tanh(-30.d0*(rho(NR) - 0.98d0)) + 1.d0 ) * zvpch
!          zvpch = (1.d0 - 0.8d0 * exp(- 0.5d0 * ((rho(NR) - 1.d0) / 0.05d0)**2)) * zvpch
          ! +++++++++++++++++++++++++++++++++++++++++
          if(izvpch == 0 .and. rho(NR) >= 0.1d0) then
             izvpch = 1
             do NR1 = 1, NR
                ! replace Vpch by the linearly interpolated value inside rho~0.1
                Vmps(NR1,1:NSM) = FSMPCH(1:NSM) * (zvpch / rho(NR) * (rho(NR1) - rho(NR)) + zvpch)
             end do
          else
             Vmps(NR,1:NSM) = FSMPCH(1:NSM) * zvpch
          end if
       else
          Vmps(NR,1:NSM) = 0.d0
       end if

       !   ***** Residual stress *****

       PiRess(NR,1:NSM) = 0.d0

       !     *** Turbulent transport of heat ***

       if (RHO(NR) < RhoSOL) then
          DeL = diff_prof(RHO(NR),FSDFIX(3),PROFCL,PROFC1,PROFCB,0.d0) + FSANOM(3) * Dturb
!pedestal          if(rho(nr) > 0.9d0) DeL = DeL * exp(-120.d0*(rho(nr)-0.9d0)**2)
       else
          if(FSPCL(3) == 0.d0) then
             factor_bohm = (FSDFIX(3) * PROFCL + PROFCB + FSANOM(3) * Dturb) &
                  &  / (Var(NRA,1)%T * rKeV / (16.d0 * AEE * sqrt(bbt(NRA))))
!bohm_model2             DeL =  (1.d0 - MOD(FSBOHM,2.d0)) * FSDFIX(3) * PROFCL &
!bohm_model2                  &+ FSBOHM * Var(NR,1)%T * rKeV / (16.d0 * AEE * BBL)
             DeL = factor_bohm * Var(NR,1)%T * rKeV / (16.d0 * AEE * BBL)
          else
             if(FSANOM(3) == 0.d0) then
                DeL = FSPCL(3) * diff_prof(RhoSOL,FSDFIX(3),PROFCL,PROFC1,PROFCB,0.d0)
             else
                DeL = FSANOM(3) * DturbA
             end if
!pedestal             DeL = FSPCL(3) * FSDFIX(3) * PROFCL * exp(-120.d0*(rho(nra)-0.9d0)**2)
          end if
       end if
!       DeL = 3.d0
       ! Thermal diffusivity
       Chis(NR,:) = Chis0(:) * DeL

       ! <omega/m>
!!$       WPM(NR) = WPM0 * Var(NR,1)%T * rKeV / (ra**2 * AEE * BphV(NR))
       WPM(NR) = 0.d0
       ! Ad hoc turbulent pinch velocity
       if(NR /= 0) then
          FVpch(NR) = Var(NR,1)%T * rKilo * sdt(NR) * VWpch(NR)
       else
          FVpch(NR) = 0.d0
       end if

       ! Coefficient of the force causing the quasilinear flux, induced by drift wave
       FQLcoef (NR) = sst(NR) * sdt(NR) * Dfs(NR,1) / (Var(NR,1)%T * rKilo)
       FQLcoef1(NR) = FQLcoef(NR) * sdt(NR) / fipol(NR)
       if(NR /= 0) then
          FQLcoef2(NR) = FQLcoef(NR) * sdt(NR) * fipol(NR) / bri(NR)
       else
          FQLcoef2(0) = 0.d0
       end if

       !     *** Loss to divertor ***

!       IF (rpt(NR) + DBW > ra) THEN
       if (rho(NR) > 1.d0) then
!          Cs = sqrt(2.d0 * Var(NR,1)%T * rKilo / amas(2) / amqp)
          Cs = sqrt((achg(2) * Var(NR,1)%T + 3.d0 * Var(NR,2)%T) * rKilo / (amas(2) * amqp))
          Lc = 2.d0 * PI * Q(NR) * RR ! Connection length to the divertor
          RL = (rpt(NR) - ra) / DBW! / 2.d0
          rNuL  (NR) = FSLP  * Cs / Lc &
               &             * RL*RL / (1.d0 + RL*RL)
          ! Classical heat conduction [s**4/(kg**2.5*m**6)]
          ! (C S Pitcher and P C Stangeby, PPCF 39 (1997) 779)
          Chicl = (4.d0*PI*EPS0)**2 &
               & /(  sqrt(amas(1)*amp)*AEE**4*Zeff(NR) &
               &   * coulog(Zeff(NR),Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amas(1),achg(1),amas(2),achg(2)))

          ! When calculating rNuLTs(1), we fix Var(:,1)%n and Var(:,1)%T constant during iteration
          !   to obain good convergence.
          rNuLTs(NR,1) = FSLTs(1) * Chicl * (PTsV_FIX(NR,1)*rKeV)**2.5d0 &
               &                  /(Lc**2 * PNsV_FIX(NR,1)*1.d20) &
               &                  * RL*RL / (1.d0 + RL*RL)
          do i = 2, NSM
             rNuLTs(NR,i) = FSLTs(i) * Cs / Lc &
                  &                  * RL*RL / (1.d0 + RL*RL)
          end do
!!$          IF(abs(FSRP) > 0.d0) THEN
             UbparaL = BUbparV(NR) / BBL
!             IF(NR == NRMAX) Ubpara(NR) = AITKEN2P(rpt(NRMAX), &
!                  & Ubpara(NRMAX-1),Ubpara(NRMAX-2),Ubpara(NRMAX-3),&
!                  & rpt(NRMAX-1),rpt(NRMAX-2),rpt(NRMAX-3))
             UbparaL = max(UbparaL, FSLP*Cs)
             rNuLB(NR) = FSLPB * UbparaL / Lc &
                  &          * RL*RL / (1.d0 + RL*RL)
!!$          END IF
       else
          rNuL(NR)   = 0.d0
          rNuLTs(NR,:) = 0.d0
          rNuLB(NR)  = 0.d0
       end if

       !     *** Current density profiles ***

       ! Parallel beam current density
       BJNB(NR)  = achgb * aee * PNbV(NR) * BUbparV(NR) * 1.d20
       ! Parallel current density <Bj//>
       BJPARA(NR)= aee * sum(achg(:)*Var(NR,:)%n*Var(NR,:)%BUpar) * 1.d20 &
            &    + BJNB(NR)
       ! Toroidal beam current density : j_{b,phi} = <j_{b,zeta}/R>/<1/R>
       AJNB(NR)  = achgb * aee * PNbV(NR) * UbphVR(NR) * 1.d20 / ait(NR)
       ! Total toroidal current density : j_phi = <j_zeta/R>/<1/R>
       AJ(NR)    = aee * sum(achg(:)*Var(NR,:)%n*Var(NR,:)%UphR) * 1.d20 / ait(NR) &
            &    + AJNB(NR)

       !     *** Equipartition power for graphics ***

       PEQei(NR) = 1.5d0 * rNuTei(NR) * Var(NR,1)%n * 1.d20 * (Var(NR,1)%T - Var(NR,2)%T) * rKeV
       PEQez(NR) = 1.5d0 * rNuTez(NR) * Var(NR,1)%n * 1.d20 * (Var(NR,1)%T - Var(NR,3)%T) * rKeV
       PEQiz(NR) = 1.5d0 * rNuTiz(NR) * Var(NR,2)%n * 1.d20 * (Var(NR,2)%T - Var(NR,3)%T) * rKeV

       !     *** Ohmic heating power ***
       ! POH = <j.E> = -PsitdotV * sum_s (e_s n_s Var%Uthhat) + PsidotV * sum_s (e_s n_s Var%UphR)
       do i = 1, NSM
          POHs(NR,i) = achg(i) * aee * Var(NR,i)%n * 1.d20 &
               &   * ( - PsitdotV(NR) * Var(NR,i)%Uthhat + PsidotV(NR) * Var(NR,i)%UphR )
       end do
       POH(NR) = sum(POHs(NR,:))

       !     *** Bremsstraulung loss ***
       !     (NRL Plasma Formulary p58 Eq. (30) (2019))
!       PBr(NR) = 1.69d-38 * Var(NR,1)%n * sqrt(Var(NR,1)%T * rKilo) * sum(achg(2:NSM)**2 * Var(NR,2:NSM)%n) * 1.d40
       !     (J. Wesson, Tokamaks 4th ed., p.229, Eq. 4.25.5)
       PBr(NR) = 5.35D-37 * Var(NR,1)%n * sum(achg(2:NSM)**2 * Var(NR,2:NSM)%n) * 1.d40 * sqrt(Var(NR,1)%T)

       !     *** Particle diffusion due to magnetic braiding ***AF 2008-06-08

       if (rpt(NR) > RMAGMN .and. rpt(NR) < RMAGMX) then
          DMAG(NR)=DMAG0*16.d0*(rpt(NR)-RMAGMN)**2*(RMAGMX-rpt(NR))**2 &
          &                   /(RMAGMX-RMAGMN)**4
          DMAGe(NR)=DMAG(NR)*Vte
          DMAGi(NR)=DMAG(NR)*Vti
       else
          DMAG(NR)=0.d0
          DMAGe(NR)=0.d0
          DMAGi(NR)=0.d0
       end if

    end do L_NR

    rNuAss(0,1) = AITKEN2P(rpt(0),rNuAss(1,1),rNuAss(2,1),rNuAss(3,1),rpt(1),rpt(2),rpt(3))
    rNuAss(0,2) = AITKEN2P(rpt(0),rNuAss(1,2),rNuAss(2,2),rNuAss(3,2),rpt(1),rpt(2),rpt(3))
    rNuAss(0,3) = AITKEN2P(rpt(0),rNuAss(1,3),rNuAss(2,3),rNuAss(3,3),rpt(1),rpt(2),rpt(3))

    !  Linear extrapolation

!!$    ChiNCp(0,:) = 2.d0 * ChiNCp(0,:) - ChiNCp(1,:)
!!$    ChiNCt(0,:) = 2.d0 * ChiNCt(0,:) - ChiNCt(1,:)
    ChiNCp(0,:) = ChiNCp(1,:)
    ChiNCt(0,:) = ChiNCt(1,:)
!    ETAvar(0,2)  = 2.d0 * ETAvar(0,2)  - ETAvar(1,2)
!    BJBSvar(0,2) = 2.d0 * BJBSvar(0,2) - BJBSvar(1,2)
!!$    ! For Neumann condition, finite viscosity is required at the magnetic axis.
!!$    Dfs(0,:)  = Dfs(1,:)
!!$    rMus(0,:) = rMus(1,:)

    !     *** Linear growth rate for toroidal gamma_etai branch of the ITG mode ***
    !        (F.Crisanti et al, NF 41 (2001) 883)
    allocate(dNsdrho, dTsdrho, mold=array_init_NRNS)
    do i = 1, NSM
       dNsdrho(:,i) = vro(:) * dNsdV(:,i)
       dTsdrho(:,i) = vro(:) * dTsdV(:,i)
    end do
    gamITG(:,1) = 0.1d0 * sqrt(Var(:,1)%T*rKilo/(amas(2)*amqp))/ra * sqrt(ra/RR) &
         &              * sqrt(abs(dNsdrho(:,2))/Var(:,2)%n &
         &                   + abs(dTsdrho(:,2))/Var(:,2)%T) &
         &              * sqrt(Var(:,2)%T/Var(:,1)%T)

    gamITG(0,2:3) = 0.d0
    i = 0
    do NR = 1, NRMAX
       if(dNsdrho(NR,2) /= 0.d0) then
          Ln = Var(NR,2)%n / abs(dNsdrho(NR,2))
       else
          i = i + 1
       end if
       if(dTsdrho(NR,2) /= 0.d0) then
          LT = Var(NR,2)%T / abs(dTsdrho(NR,2))
       else
          i = i + 1
       end if
       if(i /= 0) then
          gamITG(NR,2) = 0.d0
          gamITG(NR,3) = 0.d0
       else

          !     *** Linear growth rate valid for low abs(S) ***
          !        (A.L.Rogister, NF 41 (2001) 1101)
          !        (B.Esposito et al, PPCF 45 (2003) 933)
          etai_chk =  Ln / LT - 2.d0 / 3.d0
          if(etai_chk < 0.d0) then
             gamITG(NR,2) = 0.d0 ! marginally stable
          else
             gamITG(NR,2) = sqrt(Ln / LT - 2.d0 / 3.d0) * abs(S(NR)) &
                  &       * sqrt(Var(NR,2)%T*rKilo/(amas(2)*amqp)) / (Q(NR) * RR)
          end if

          !     *** Linear growth rate for q>2 and s=0 ***
          !        (J.Candy, PoP 11 (2004) 1879)
          kthrhos = sqrt(0.1d0)
          gamITG(NR,3) = kthrhos * sqrt(Var(NR,1)%T*rKilo/(amas(2)*amqp)) / ra &
               &       * sqrt(2.d0 * ra / RR * (1.d0 / Ln + 1.d0 / LT))
       end if
    end do

    if(MDANOMabs == 3) &
         & call txmmm95(dNsdrho,dTsdrho,dQdrho,FSCBSH*rG1)

    !     *** ETB model ***

    if(MDLETB /= 0) then ! ETB on
       Frdc = 0.1d0
       if(Fcoef > Frdc) then
          Dcoef = (1.d0 - Frdc) / NTMAX
          Fcoef = 1.d0 - NT * Dcoef
       end if
       do NR = 0, NRA
          if(RhoETB(1) /= 0.d0) then
             if(Rho(NR) > RhoETB(1)) then
                Dfs(NR,:) = Dfs(NR,:) * Fcoef
             end if
          end if
          if(RhoETB(2) /= 0.d0) then
             if(Rho(NR) > RhoETB(2)) then
                rMus(NR,:) = rMus(NR,:) * Fcoef
             end if
          end if
          if(RhoETB(3) /= 0.d0) then
             if(Rho(NR) > RhoETB(3)) then
                Chis(NR,:) = Chis(NR,:) * Fcoef
             end if
          end if
       end do
    end if

    !     *** Resistivity and current density ***

    do NR = 0, NRMAX
       !     *** Resistivity by Hirshman, Hawryluk and Birge model ***
       !     (S.P. Hirshman, R.J. Hawryluk and B. Birge, Nucl. Fusion 17 (1977) 3)
       zrNuee = coll_freq(NR,1,1)
       Vte = sqrt(2.d0 * Var(NR,1)%T * rKilo / amas(1) * amqp)
       zxi = 0.58d0 + 0.2d0 * Zeff(NR)
       znusteinv = epst(NR)**1.5d0 * Vte / (sqrt(2.d0) * RR * Q(NR) * zrNuee)
       zLamEZ = 3.4d0 / Zeff(NR) * (1.13d0 + Zeff(NR)) / (2.67d0 + Zeff(NR))
       zCRZ = 0.56d0 / Zeff(NR) * (3.d0 - Zeff(NR)) / (3.d0 + Zeff(NR))
       zsigma0 = Var(NR,1)%n * 1.d20 * aee / (amas(1) * amqp * zrNuee)
       zsigmaneo = zsigma0 * zLamEZ * (1.d0 -        ft(NR) * znusteinv / (znusteinv + zxi)) &
            &                       * (1.d0 - zCRZ * ft(NR) * znusteinv / (znusteinv + zxi))
       ETAvar(NR,4) = 1.d0 / zsigmaneo
       ! Spitzer resistivity for hydrogen plasma (parallel direction)
       ETAS(NR) = CORR(1.d0) / zsigma0
    end do

    if( MDLNEO > 10 ) then
       !     *** Bootstrap current and resistivity by Sauter model ***
       !     (O. Sauter, et al.,  Phys. Plasmas 6 (1999) 2834, ibid. 9 (2002) 5140)
       NR = 0
       ETAvar(NR,3) = ETAS(NR) * (CORR(Zeff(NR)) / CORR(1.d0))
       BJBSvar(NR,3) = 0.d0

       do NR = 1, NRMAX
          call sauter(Var(NR,1)%n,Var(NR,1)%T,dTsdV(NR,1),dPsdV(NR,1), &
               &      Var(NR,2)%n,Var(NR,2)%T,dTsdV(NR,2),dPsdV(NR,2), &
               &      Q(NR),sdt(NR),fipol(NR),epst(NR),RR,achg(2),Zeff(NR),ft(NR), &
               &      rlnLei_IN=coulog(Zeff(NR),Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T, &
               &                       amas(1),achg(1),amas(2),achg(2)), &
               &      rlnLii_IN=coulog(Zeff(NR),Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T, &
               &                       amas(2),achg(2),amas(2),achg(2)), &
               &      BJBS=BJBSvar(NR,3),ETA=ETAvar(NR,3))
       end do
    end if

    if(FSNC(1) /= 0.d0) then
       select case(MDLETA)
       case default
          ETA(:) = ETAvar(:,MDLNEOL) ! depending upon MDLNEO
       case(2)
          if(MDLNEO > 10) then
             ETA(:) = ETAvar(:,3) ! Sauter
          else
             ETA(:) = ETAvar(:,MDLNEOL) ! depending upon MDLNEO
          end if
       case(3)
          ETA(:) = ETAvar(:,4) ! HHB
       end select
    else
       ! Spitzer resistivity when no neoclassical effects
       ETA(:) = ETAS(:)
    end if

    if( MDBSETA == 0 ) then
       do NR = 0, NRMAX
          ! Bootstrap current density estimated by the neoclassical transport model
          select case(MDLNEOL)
          case(1:2) ! MI
             BJBS(NR) = BJBSvar(NR,1)
          case(3)   ! NCLASS
             BJBS(NR) = BJBSvar(NR,2)
          case default
             BJBS(NR) = BJBSvar(NR,3)
          end select

          ! Parallel Ohmic current density
          BJOH(NR) = BJPARA(NR) - BJBS(NR) - BJNB(NR)
       end do
    else
       do NR = 0, NRMAX
          ! Parallel Ohmic current density using the resistivity estimation 
          BJOH(NR) = BEpara(NR) / ETA(NR)

          ! Rough estimate of the bootstrap current
          !    This way of estimation strongly depends upon how ETA(NR) in BJOH(NR) is
          !    , which means that this estimate may be no use.
          BJBS(NR) = BJPARA(NR) - BJOH(NR) - BJNB(NR)
       end do
    end if
    do NR = 0, NRMAX
       ! Toroidal bootstrap current density : 
       !   j_{BS,phi} = <j_{BS,zeta}/R>/<1/R> = (<BJ_BS><1/R^2>I/<B^2>)/<1/R>
       AJBS(NR) = BJBS(NR) * aat(NR) * fipol(NR) / (bbt(NR) * ait(NR))

       ! Toroidal Ohmic current density :
       !   j_{OH,phi} = <j_{OH,zeta}/R>/<1/R> = (<BJ_OH><1/R^2>I/<B^2>)/<1/R>
       AJOH(NR) = BJOH(NR) * aat(NR) * fipol(NR) / (bbt(NR) * ait(NR))

       ! Rough estimate of the resistivity using the bootstrap estimation
       !    ETA is less sensitive than the bootstrap current.
!       ETAvar(NR,0) = BEpara(NR) / (BJPARA(NR) - BJBSvar(NR,MDLNEOL) - BJNB(NR))
       if( BJOH(NR) /= 0.d0 ) then
          ETAvar(NR,0) = BEpara(NR) / BJOH(NR)
       else
          ! BJOH = 0.d0 sometimes happen in the SOL region.
          ETAvar(NR,0) = 0.d0
       end if
!       write(6,'(F8.5,3ES15.7)') rho(nr),BEpara(NR),BJOH(NR),ETAvar(NR,0)
    end do
!    write(6,*)

    !     ***** Ion Orbit Loss *****

    if (abs(FSLC) > 0.d0) then
       SiLC  (0) = 0.d0
       SiLCB (0) = 0.d0
       SiLCph(0) = 0.d0
       rNuOL (0) = 0.d0

       if(MDLC == 1) then
          ! Ref. [K. C. Shaing, Phys. Fluids B 4 (1992) 3310]
          do NR = 1, NRMAX
             Vti = sqrt(2.d0 * Var(NR,2)%T * rKilo / (amas(2) * amqp))
             ! Orbit squeezing factor, given just below Eq. (3) of Ref.
             sqz = 1.d0 + (fipol(NR) / bb)**2 * amqp * amas(2) / achg(2) * abs(ddPhidpsi(NR))

             EpsL = epst(NR)
             BBL = sqrt(bbt(NR))

             ! rNuDL : deflection collisional frequency at V = Vti, see [Kikuchi PPCF 1995 p.1236]
             rNuDL = 0.d0
             do i = 1, NSM
                xb    = Vti / sqrt(2.d0 * Var(NR,i)%T * rKilo / (amas(i) * amqp))
                fp    = erf(xb)
                fdp   = 2.d0/sqrt(pi)*exp(-xb**2) ! derivative of erf
                gfun  = (fp-xb*fdp)/(2.d0*xb**2)  ! Chandrasekhar function
                rNuDL = rNuDL + 0.75d0 * sqrt(pi) &
                     & * (Var(NR,i)%n * achg(i)**2) / (Var(NR,2)%n * achg(2)**2) &
                     & * (fp - gfun) ! Actually (fp-gfun)/xa**3, but xa at V=Vti is equal to unity.
             end do
             rNuDL = rNuDL * coll_freq(NR,2,2)

             rNustar = rNuDL * RR * Q(NR) / (Vti * (abs(sqz) * EpsL)**1.5d0)

             rNuOLL  = 2.25d0 * rNuDL / (sqrt(PI) * sqrt(2.d0 * abs(sqz) * EpsL)) &
                  &  * exp(-(rNustar**0.25d0 + (achg(2) * BBL / (amas(2) * amqp)) &
                  &  * sqrt(abs(sqz)) / (fipol(NR) * Vti) &
                  &  * abs(PsiV(NR) - PsiV(NRA)) / sqrt(2.d0 * EpsL))**2)

             if(FSLC <= 1.d0) then
                SiLC  (NR) = 0.d0
                SiLCB (NR) = 0.d0
                SiLCph(NR) = 0.d0
                rNuOL (NR) = FSLC * rNuOLL
             else
                SiLC  (NR) = - mod(FSLC,1.d1) * Var(NR,2)%n * rNuOLL
                SiLCB (NR) = SiLC(NR) * Var(NR,2)%BUpar
                SiLCph(NR) = SiLC(NR) * Var(NR,2)%RUph
                rNuOL (NR) = 0.d0
             end if
          end do

          do NR = 0, NRMAX
             tmp(NR) = moving_average(NR,rNuOL,NRMAX)
          end do
          rNuOL(:) = tmp(:)

       else if(MDLC == 2) then
          !     S. -I. Itoh and K. Itoh, Nucl. Fusion 29 (1989) 1031
          if(FSLC <= 1.d0) then
             ! RLOSS : Numerical coefficient proportional to the relative number of ions
             !         in the loss cone in velocity space
             RLOSS = 0.1d0
             rNuOL(0) = 0.d0
             do NR = 1, NRMAX
                EpsL = epst(NR)
                Vti = sqrt(Var(NR,2)%T * rKilo / (amas(2) * amqp))
                RhoIT = Vti * amas(2) * amqp / (achg(2) * BthV(NR))
                RL = (rpt(NR) - (ra - 1.5d0 * RhoIT)) / DBW ! Alleviation factor
                if(rpt(NR) > (ra - RhoIT)) then
!                if(abs(ra - rpt(NR)) <= RhoIT .AND. RHO(NR) < 1.d0) then
                   ExpArg = -2.d0 * EpsL * (ErVlc(NR) / BthV(NR))**2 / Vti**2
                   ExpArg = ExpArg * (rpt(NR) / ra)**2
                   rNuOL(NR) = FSLC * RLOSS * rNuii(NR) / sqrt(EpsL) * exp(ExpArg) &
                        &    * RL**2 / (1.d0 + RL**2)
                else
                   rNuOL(NR) = 0.d0
                end if
             end do

             SiLC  (:) = 0.d0
             SiLCB (:) = 0.d0
             SiLCph(:) = 0.d0

          else
             do NR = 1, NRA
                EpsL = epst(NR)
                Vti = sqrt(Var(NR,2)%T * rKilo / (amas(2) * amqp))
                RhoIT = Vti * amas(2) * amqp / (achg(2) * BthV(NR))
                RhoIT = min(RhoIT,0.1d0)
                Wti = Vti / (Q(NR) * RR)
                rNuAsI_inv = EpsL**1.5d0 * Wti / (sqrt(2.d0) * rNuii(NR))
                ExpArg = 2.d0 * EpsL / Vti**2 * (ErVlc(NR) / BthV(NR))**2
                AiP = rNuii(NR) * sqrt(EpsL) * rNuAsI_inv / (1.d0 + rNuAsI_inv) &
                     & * exp(- ExpArg)
                do NR1 = NRA, NRMAX
                   DISTAN = (rpt(NR1) - rpt(NR)) / RhoIT
                   SiLCL = AiP * exp( - DISTAN**2) * Var(NR,2)%n
                   SiLC(NR) = SiLC(NR) - SiLCL
                   SiLC(NR1) = SiLC(NR1) + SiLCL * rpt(NR) / rpt(NR1)

                   SiLCBL = SiLCL * (amas(2)*amp) * Var(NR,2)%Uth * rpt(NR)
                   SiLCB(NR) = SiLCB(NR) - SiLCBL
                   SiLCB(NR1) = SiLCB(NR1) + SiLCBL * rpt(NR) / rpt(NR1)

                   SiLCphL = SiLCL * (amas(2)*amp) * Var(NR,2)%Uph
                   SiLCph(NR) = SiLCph(NR) - SiLCphL
                   SiLCph(NR1) = SiLCph(NR1) + SiLCphL * rpt(NR) / rpt(NR1)
                end do
             end do

             SiLC  (:) = FSLC * SiLC  (:)
             SiLCB (:) = FSLC * SiLCB (:)
             SiLCph(:) = FSLC * SiLCph(:)
             rNuOL (:) = 0.d0
          end if
       end if
    end if

    !     ***** Toroidal ripple effect *****
    call ripple_effect(dQdrho)

!    !     ***** Neoclassical toroidal viscosity (NTV) *****
!    !      "rNuNTV" and "UastNC" are obtained from NTVcalc
!    
!    CALL NTVcalc
!    rNuNTV(:) = 0.d0
!    UastNC(:) = 0.d0

    deallocate(dErdr,dpdr,dErdrS,ErVlc)
    deallocate(dQdrho,dlnNedrhov)
    deallocate(dTsdV,dTsdrho,dPsdV,dNsdrho)

  contains

!***************************************************************
!
!   Equipartition frequency
!
!***************************************************************

    real(8) function rNuTeq(NR,i,j)
      use tx_interface, only : coulog
      integer(4), intent(in) :: NR, i, j
      real(8) :: tau_ij, vcoulog

      !   *** Collision frequency (energy equipartition) ***
      !   [D. V. Sivukhin, Rev. Plasma Phys. Vol. 4, Consultants Bureau (1966) Sec. 9]
      !       Energy relaxation time between two kinds of particles
      !   See also [J.F. Artaud et al., Nucl Fusion 50 (2010) 043001]

      !   tau_ij : momentum exchange time between i and j species
      !   taueq  : energy exchange time between i and j species
      !      taueq = m_j / (2 m_i) * tau_ij
      !   rNuTeq = 1 / taueq
 
      vcoulog = coulog(Zeff(NR),Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amas(i),achg(i),amas(j),achg(j))

      tau_ij  = 6.d0 * Pi * sqrt(2.d0 * Pi) * eps0**2 / (aee**4 * vcoulog) &
           & * amas(i)**2 * sqrt(amp) / ( Var(NR,j)%n * 1.d20 * ( achg(i) * achg(j) )**2 ) &
           & * ( ( Var(NR,i)%T / amas(i) + Var(NR,j)%T / amas(j) ) * rKeV )**1.5d0

      rNuTeq = 2.d0 * amas(i) / ( amas(j) * tau_ij )

    end function rNuTeq

!***************************************************************
!
!   Given diffusion coefficient profile
!     Input : factor : FSDFIX
!             profd  : shape factor
!             rho    : normalized radius
!             npower : power of rho
!             base   : pedestal
!             fgfact : switch for superimpose of Gaussian profile
!          <optional> (valid when fgfact /= 0)
!             mu     : average
!             sigma  : standard deviation
!
!     diff_prof = factor         at rho=0
!                 factor * profd at rho=1
!
!***************************************************************

    real(8) function diff_prof(rho,factor,profd,power,base,fgfact,mu,sigma)
      use tx_interface, only : fgaussian
      real(8), intent(in) :: rho, factor, profd, power, base, fgfact
      real(8), intent(in), optional :: mu, sigma
      real(8) :: fmod, fmodmax

      ! Gaussian profile modified by parabolic profile
      if(fgfact /= 0.d0) then
         fmod    = fgaussian(rho,mu,sigma) * (- 4.d0 * (rho - 0.5d0)**2 + 1.d0)
         fmodmax = fgaussian(mu, mu,sigma) * (- 4.d0 * (mu  - 0.5d0)**2 + 1.d0)
         fmod = fgfact * fmod / fmodmax + base
      else
         fmod = base
      end if

      ! Sum  jof usual and modified Gaussian profiles
      diff_prof = factor * (1.d0 + (profd - 1.d0) * rho**power) + fmod

    end function diff_prof

!***************************************************************
!
!   Rate coefficients for electron impact hydrogen (e+H) ionization process
!     valid for 1eV   => 10^5eV
!   Rate coefficients for proton   impact hydrogen (p+H) ionization process
!     valid for 100eV => 10^7eV
! 
!   Cross-section data:
!     (R.L. Freeman and E.M. Jones, Culham Laboratory Report, CLM-R-137 (May 1974))
!   Example of use:
!     (C.E. Singer et al., Comp. Phys. Comm. 49 (1988) 275, Eq. (2.9.5l))
!
!     Inputs (real*8): tekev : Electron temperature [keV]
!            (real*8): tikev : Proton   temperature [keV]
!     Output (real*8): SiViz : Ionization maxwellian rate coefficient [m^3/s]
!
!     Internal (real*8): a, b : Table of Maxwellian rate coefficients in TABLE 3
!                                        ^^^^^^^^^^
!***************************************************************

    real(8) function SiViz(tekev,tikev)

      real(8), intent(in) :: tekev, tikev
      real(8) :: x, tekev_temp, SiVize, tikev_temp, SiVizp
      real(8), dimension(0:6) :: a
      real(8), dimension(0:8) :: b
      data a /-0.3173850d02, 0.1143818d02, -0.3833998d01,  0.7046692d0, &
           &  -0.7431486d-1, 0.4153749d-2, -0.9486967d-4/
      data b /-0.1490861d03, 0.7592575d02, -0.2209281d02,  0.3909709d01, &
           &  -0.4402168d00, 0.3209047d-1, -0.1493409d-2,  0.4094151d-4, -0.5069777d-6/

      ! Electron ionization

      if(tekev < 1.d-3) then
         write(6,'(A,ES12.4)') &
              'Function SiViz: Out of energy range. tekev=', tekev
         tekev_temp=1.D-3
      else if(tekev > 1.d2) then
!!$         write(6,'(A,ES12.4)') &
!!$              'Function SiViz: Out of energy range. tekev=', tekev
         tekev_temp=1.D2
      else
         tekev_temp=tekev
      endif
      x = log(tekev_temp * 1.d3)
      SiVize = exp(a(0)+x*(a(1)+x*(a(2)+x*(a(3)+x*(a(4)+x*(a(5)+x*a(6)))))))*1.D-6

      SiViz = SiVize

      ! Above 6 keV, proton ionization process exceeds electron's.
      ! To save computation, the following evaluation is not done below 6 keV.

      if( tekev > 6.d0 .and. tikev > 6.d0 ) then

         ! Proton ionization

!!$         if(tikev < 1.d-1) then
!!$            write(6,'(A,ES12.4)') &
!!$                 'Function SiViz: Out of energy range. tikev=', tikev
!!$            tikev_temp=1.D-3
!!$         else if(tikev > 1.d4) then
         if(tikev > 1.d4) then
            write(6,'(A,ES12.4)') &
                 'Function SiViz: Out of energy range. tikev=', tikev
            tikev_temp=1.D4
         else
            tikev_temp=tikev
         endif
         x = log(tikev_temp * 1.d3)
         SiVizp = exp( b(0)+x*(b(1)+x*(b(2)+x*(b(3)+x*(b(4) &
              &      +x*(b(5)+x*(b(6)+x*(b(7)+x*b(8)))))))))*1.D-6

         SiViz = max(SiVize, SiVizp)

      end if

    end function SiViz

!***************************************************************
!
!   Rate coefficients for charge exchange cross-section of protons on atomic hydrogen
!     valid for 1eV => 10^5eV
! 
!   Cross-section data:
!     (R.L. Freeman and E.M. Jones, Culham Laboratory Report, CLM-R-137 (May 1974))
!
!     Inputs (real*8): tikev : Ion temperature [keV]
!     Output (real*8): SiVcx : Charge exchange maxwellian rate coefficient [m^3/s]
!
!     Internal (real*8): a : Table of Maxwellian rate coefficients in TABLE 3
!                                     ^^^^^^^^^^
!***************************************************************

    real(8) function SiVcx(tikev)

      real(8), intent(in) :: tikev
      real(8) :: x, tikev_temp
      real(8), dimension(0:8) :: a
      data a /-0.1841757d02, 0.5282950d0, -0.2200477d0,   0.9750192d-1, &
           &  -0.1749183d-1, 0.4954298d-3, 0.2174910d-3, -0.2530205d-4, 0.8230751d-6/

      if(tikev < 1.d-3) then
         write(6,'(A,ES12.4)') &
              'Function SiVcx: Out of energy range. tikev=', tikev
         tikev_temp=1.D-3
      else if(tikev > 1.d2) then
         write(6,'(A,ES12.4)') &
              'Function SiVcx: Out of energy range. tikev=', tikev
         tikev_temp=1.D2
      else
         tikev_temp=tikev
      endif
      x = log(tikev_temp * 1.d3)
      SiVcx = exp( a(0)+x*(a(1)+x*(a(2)+x*(a(3)+x*(a(4) &
           &      +x*(a(5)+x*(a(6)+x*(a(7)+x*a(8)))))))))*1.D-6

    end function SiVcx

!***************************************************************
!
!   Reduce diffusivity for thermal neutrals LQn2 due to cylindrical geometry
!
!     Inputs (integer*4): NRctr  : interest radial grid number
!            (real*8)   : rLmean : mean free path of LQn2 at NRctr [m]
!     Output (real*8)   : reduce_D02 : Reduction factor of D02 [*]
!
!***************************************************************

    real(8) function reduce_D02(NRctr,rLmean)
      use tx_commons, only : Pi, NRMAX, rpt
      integer(4), intent(in) :: NRctr
      real(8), intent(in) :: rLmean

      integer(4) :: nr, nr0, idebug = 0
      real(8) :: Rctr, costh0, theta0, frac, DltL, costh, theta, theta1, rLmean_eff, rLmean_av

      Rctr = rpt(NRctr)
      costh0 = 0.5d0 * rLmean / Rctr
      if(abs(costh0) > 1.d0) then
         reduce_D02 = 1.d0
         return
      end if
      theta0 = 2.d0 * acos(costh0)
      frac = theta0 / Pi
      if(idebug /= 0) write(6,*) "frac=",frac,"rLmean=",rLmean
      
      DltL = abs(Rctr - rLmean)
      do nr = 0, nrmax
         if(rpt(nr) > DltL) then
            nr0 = nr
            exit
         end if
      end do

      rLmean_av = 0.d0
      theta     = 0.d0
      do nr = nr0, nrctr
         costh  = (Rctr**2 + rLmean**2 - rpt(nr)**2) / (2.d0 * Rctr * rLmean)
         theta1 = 2.d0 * acos(costh)
         theta  = theta1 - theta
!         rLmean_eff = rLmean * costh
         rLmean_eff = Rctr - rpt(nr)
         rLmean_av  = rLmean_av + rLmean_eff * (theta / theta0)
         if(idebug /= 0) write(6,'(I3,5F15.7)') nr,theta1,theta,theta/theta0,rLmean_eff,rLmean_av
         theta  = theta1
      end do
 
      reduce_D02 = frac * (rLmean_av / rLmean)
      if(idebug /= 0) write(6,*) "reduce_D02=",reduce_D02

    end function reduce_D02

  end subroutine TXCALC

end module tx_variables
