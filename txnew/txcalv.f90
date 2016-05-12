!     $Id: txcalv.f90,v 1.112 2011/06/13 07:53:20 honda Exp $
module tx_variables
  implicit none
  public

contains

!***************************************************************
!
!        Calculate variables
!
!***************************************************************

  SUBROUTINE TXCALV(XL,ID)

    use tx_commons
    use tx_interface, only : dfdx, replace_interpolate_value
    use tx_core_module, only : intg_vol_p
!    use aux_system, only : Vbabsmax
    REAL(8), DIMENSION(0:NRMAX,1:NQMAX), INTENT(IN) :: XL
    integer(4), intent(in), optional :: ID
    INTEGER(4) :: NR, i, JSMTHD, n1, n2
    real(8), parameter :: fourPisq = 4.d0 * Pi * Pi
    real(8) :: sdtvac, dPsitVdVvac, sum_int!, BBL
    real(8), dimension(:), allocatable :: dPhidpsiL
    real(8) :: AITKEN2P ! function

    JSMTHD = ISMTHD / 10

    IF(present(ID)) THEN
       ! The pres0 and ErV0 are the values evaluated at the previous time step
       !   for numerical stability when using turbulent transport models.
       IF(MDFIXT == 0) THEN
          pres0(:) = (  XL(:,LQe5) + XL(:,LQi5)) * 1.D20 * rKeV
       ELSE
          pres0(:) = (  XL(:,LQe1) * XL(:,LQe5) &
               &      + XL(:,LQi1) * XL(:,LQi5)) * 1.D20 * rKeV
       END IF
       ErV0 (0)       =   0.d0
       if(ISMTHD == 0) then
          dPhiV(:) =   dfdx(vv,XL(:,LQm1),NRMAX,0) ! dPhiV/dV
          ErV0 (1:NRMAX) = - sst(1:NRMAX) / vro(1:NRMAX) * rbvl * dPhiV(1:NRMAX)
       else
          dPhiV(:) =   dfdx(rho,XL(:,LQm1),NRMAX,0,daxs=0.d0) ! dPhiV/drho
          if(JSMTHD == 1) call replace_interpolate_value(dPhiV(1),1,rho,dPhiV)
          ErV0 (1:NRMAX) = - sst(1:NRMAX) / vro(1:NRMAX)**2 * rbvl * dPhiV(1:NRMAX)
       end if
       IF(ID /= 0) return
    END IF

    allocate(dPhidpsiL(0:NRMAX))

!$omp parallel
!$omp workshare
    PhiV (:) =   XL(:,LQm1)

    PsitdotV(:)= XL(:,LQm2)

    PsidotV(:) = XL(:,LQm3)

    ! Etor: <E_t> = <1/R>dpsi/dt
    Etor (:) =   PsidotV(:) * ait(:)

    PsiV (:) =   XL(:,LQm4) * rMUb2
!$omp end workshare
    ! sdt: dpsi/dV
    !    sdt(NRMAX) =2.D0*Pi*rMU0*rIp*1.D6/ckt(NRMAX) is a boundary condition.
    sdtvac = 2.d0 * Pi * rMU0 * rIp * 1.D6 / ckt(NRMAX)
    sdt  (:) =   dfdx(vv,PsiV,NRMAX,0,dbnd=sdtvac)

    ErV  (0)       =   0.d0
    if(ISMTHD == 0) then
       dPhiV(:) =   dfdx(vv,PhiV,NRMAX,0) ! dPhiV/dV
       ErV  (1:NRMAX) = - sst(1:NRMAX) / vro(1:NRMAX) * rbvl * dPhiV(1:NRMAX)
       dPhidpsiL(:) = dPhiV(:) / sdt(:) ! dPhi/dpsi
    else
       dPhiV(:) =   dfdx(rho,PhiV,NRMAX,0,daxs=0.d0) ! dPhiV/drho
       ! Replace dPhiV(1) by the interpolated value
       if(JSMTHD == 1) call replace_interpolate_value(dPhiV(1),1,rho,dPhiV)
       ErV  (1:NRMAX) = - sst(1:NRMAX) / vro(1:NRMAX)**2 * rbvl * dPhiV(1:NRMAX)
       dPhidpsiL(1:NRMAX) = dPhiV(1:NRMAX) / (sdt(1:NRMAX) * vro(1:NRMAX))! dPhi/dpsi
!       dPhidpsiL(0) = AITKEN2P(rho(0),dPhidpsiL(1),dPhidpsiL(2),dPhidpsiL(3) &
!            &                 ,rho(1),rho(2),rho(3))
       dPhidpsiL(0) = AITKEN2P(vv(0),dPhidpsiL(1),dPhidpsiL(2),dPhidpsiL(3) &
            &                 ,vv(1),vv(2),vv(3))
    end if

!$omp workshare
!    BthV (:) =   fourPisq * rho(:) * rbvl * sdt(:) ! NOT FSA quantity
    BthV (:) =   sqrt(ckt(:)) * sdt(:)
    PsitV(:) =   XL(:,LQm5) * rMU0
!    ! bbt = <B^2>
!    bbt  (:) =   fourPisq*fourPisq / aat(:) * hdt(:)*hdt(:) &
!         &     + ckt(:) * sdt(:)*sdt(:)

    Var(:,1)%n = XL(:,LQe1)
    Var(:,1)%UrV = XL(:,LQe2) / Var(:,1)%n
    Var(0,1)%Ur       = 0.d0
    Var(1:NRMAX,1)%Ur = Var(1:NRMAX,1)%UrV / suft(1:NRMAX) ! <u.nabla V>/<|nabla V|>

    Var(:,1)%RUph  = XL(:,LQe4) / Var(:,1)%n
!$omp end workshare
    IF(MDFIXT == 0) THEN
       Var(:,1)%p = XL(:,LQe5)
       Var(:,1)%T = XL(:,LQe5) / Var(:,1)%n
    ELSE
       Var(:,1)%p = XL(:,LQe5) * Var(:,1)%n
       Var(:,1)%T = XL(:,LQe5)
    END IF
!$omp workshare
    Var(:,2)%n = XL(:,LQi1)
    Var(:,2)%UrV = XL(:,LQi2) / Var(:,2)%n
    Var(0,2)%Ur        = 0.D0
    Var(1:NRMAX,2)%Ur  = Var(1:NRMAX,2)%UrV / suft(1:NRMAX) ! <u.nabla V>/<|nabla V|>

    Var(:,2)%RUph  = XL(:,LQi4) / Var(:,2)%n
!$omp end workshare
    IF(MDFIXT == 0) THEN
       Var(:,2)%p = XL(:,LQi5)
       Var(:,2)%T = XL(:,LQi5) / Var(:,2)%n
    ELSE
       Var(:,2)%p = XL(:,LQi5) * Var(:,2)%n
       Var(:,2)%T = XL(:,LQi5)
    END IF
!$omp workshare
    ! Parallel flows
    Var(:,1)%BUpar = XL(:,LQe3)
    Var(:,2)%BUpar = XL(:,LQi3)
    Var(:,1)%Bqpar = XL(:,LQe6)
    Var(:,2)%Bqpar = XL(:,LQi6)

    Var(:,1)%UphR = XL(:,LQe7) / Var(:,1)%n
    Var(:,2)%UphR = XL(:,LQi7) / Var(:,2)%n
    Var(:,1)%Uph = Var(:,1)%UphR / ait(:)
    Var(:,2)%Uph = Var(:,2)%UphR / ait(:)
!$omp end workshare
!$omp end parallel

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
    forall (NR=0:NRMAX, PNbV(NR) /= 0.d0) PNbVinv(NR) = 1.d0 / PNbV(NR)
    do NR = 0, NRMAX
       UbrVV(NR)   = XL(NR,LQb2) * PNbVinv(NR) * MDBEAM
       BUbparV(NR) = XL(NR,LQb3) * PNbVinv(NR)
       RUbphV(NR)  = XL(NR,LQb4) * PNbVinv(NR)
       UbphVR(NR)  = XL(NR,LQb7) * PNbVinv(NR)
!       UbphV(NR) = BUbparV(NR) * fipol(NR) / bbt(NR) / rr
       UbphV(NR)   = UbphVR(NR) / ait(NR)
    end do

    if(MDBEAM == 0 .and. abs(FSRP) == 0.d0) then ! No beam ions in the SOL
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
       UbphV(NRA:NRMAX)   = 0.D0
       RUbphV(NRA:NRMAX)  = 0.D0
       UbphVR(NRA:NRMAX)  = 0.D0
    end if

    PNbRPV(:)=   XL(:,LQr1)
    ! -----------------

    PN01V(:) =   XL(:,LQn1)
    PN02V(:) =   XL(:,LQn2)
    PN03V(:) =   XL(:,LQn3)

    do i = 1, NSM
       dTsdpsi(:,i) = dfdx(PsiV,Var(:,i)%T,NRMAX,0) ! dT/dpsi
    end do
    if(ISMTHD == 0) then
       do i = 1, NSM
          dPsdpsi(:,i) = dfdx(PsiV,Var(:,i)%p,NRMAX,0) ! dp/dpsi
       end do
    else
       do i = 1, NSM
          dPsdpsi(:,i) = dfdx(rho,Var(:,i)%p,NRMAX,0,daxs=0.d0) ! dp/drho
          ! Replace dPsdpsi(1) by the interpolated value
          if(JSMTHD == 1) call replace_interpolate_value(dPsdpsi(1,i),1,rho,dPsdpsi(:,i))
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
       qhatsq(NR) = 0.5d0 * fipol(NR)*fipol(NR) / Bpsq(NR) * (aat(NR) - 1.d0 / rrt(NR))
       ! Metric coefficient: <B^2><R^2>-I^2
!       bri(NR) = bbt(NR) * rrt(NR) - fipol(NR)**2 ! => causing instability through LQe2 and LQi2
       bri(NR) = rrt(NR) * Bpsq(NR) * (1.d0 + 2.d0 * qhatsq(NR))
!       if(abs(bri(NR)) < 1.d-10) bri(NR) = 0.d0
    end do
    NR=0
    bri(NR) = 0.d0
    UgV(NR) = 0.d0
    if(ieqread == 0) then ! Perfect cylinder (Pfirsch-Schluter contribution nil)
       qhatsq(NR) = 0.d0
    else if(ieqread == 1) then ! Large-aspect-ratio tokamak
       qhatsq(NR) = Q(NR)*Q(NR)
    else ! General tokamak equilibria
       qhatsq(NR) = AITKEN2P(vv(0),qhatsq(1),qhatsq(2),qhatsq(3),vv(1),vv(2),vv(3))
    end if

    do NR = 0, NRMAX
       ! Coefficients regarding Pfirsch-Schluter contribution: 2q^^2/(1+2q^^2)
       Fqhatsq(NR) = 2.d0 * qhatsq(NR) / (1.d0 + 2.d0 * qhatsq(NR))
       ! Parallel Electric field: <BE_//>
       BEpara(NR) = fipol(NR) * aat(NR) * (- PsitdotV(NR) / Q(NR) + PsidotV(NR))
       ! <B^2>
       bbt(NR) = fourPisq*fourPisq * hdt(NR)*hdt(NR) / aat(NR) + ckt(NR) * sdt(NR)*sdt(NR)
       ! <B^theta> = 4 pi^2 dpsi/dV
       bthco(NR) = fourPisq * sdt(NR)
    end do

    !  Diamagnetic particle and heat flows: B V_1s, B V_2s
    forall (i = 1:NSM, NR = 0:NRMAX) 
       ! Diamag. particle flow
       BVsdiag(NR,i,1) = - fipol(NR) &
            & * ( dPsdpsi(NR,i) * rKilo / ( achg(i) * Var(NR,i)%n ) + dPhidpsiL(NR) )
       ! Diamag. heat flow
       BVsdiag(NR,i,2) = - fipol(NR) / achg(i) * dTsdpsi(NR,i) * rKilo

       ! Var%Uthhat = U.nabla theta/B.nabla theta [m/s/T]
       !    Var%Uthhat = ( <B U_{s//}> - BV_{1s} ) / <B^2> or
       !               = ( <R U_{s zeta}> - <R^2>/I BV_{1s} ) / I, alternatively
       Var(NR,i)%Uthhat = ( Var(NR,i)%BUpar - BVsdiag(NR,i,1) ) / bbt(NR)
       ! Poloidal rotation for graphics: u_theta = uthhatV * BthV [m/s]
       Var(NR,i)%Uth = Var(NR,i)%Uthhat * BthV(NR)
    end forall

    deallocate(dPhidpsiL)

    RETURN
  END SUBROUTINE TXCALV

!***************************************************************
!
!   Calculate coefficients
!
!***************************************************************

  SUBROUTINE TXCALC(IC)

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
!    use tx_ntv, only : NTVcalc, rNuNTV, UastNC

    integer(4), intent(in) :: IC
    integer(4), save :: NRB = 1
    INTEGER(4) :: NR, NR1, IER, i, MDANOMabs, model_cdbm, izvpch, MDLNEOL
    INTEGER(4) :: NHFM ! miki_m 10-08-06
    REAL(8) :: Sigma0, Vte, Vti, Vtb, Wte, Wti, EpsL, &
         &     rNuAsE_inv, rNuAsI_inv, BBL, Va, Wpe2, PN0tot, &
         &     PROFML, PROFCL, Dturb, DeL, &
         &     Cs, Lc, RhoIT, ExpArg, AiP, DISTAN, UbparaL, &
         &     rNuOLL, SiLCL, SiLCBL, SiLCphL, RL, DBW, PTiVA, &
         &     Chicl, factor_bohm, rNustar, &
         &     RLOSS, sqz, rNuDL, Ln, LT, etai_chk, kthrhos, &
         &     RhoSOL, V0ave, Viave, DturbA, rLmean, Sitot, &
         &     rGCIM, rGIM, rHIM, OMEGAPR !09/06/17~ miki_m
    real(8), dimension(0:NRMAX) :: gr2phi
    real(8), dimension(1:NHFMmx) :: EpsLM 
!    real(8) :: rLmeanL, QL
    real(8), save :: Fcoef = 1.d0
    real(8) :: PTiVav, N02INT, RatSCX, sum1, sum2, zvpch
    real(8) :: Frdc, Dcoef
    real(8) :: omegaer, omegaere, omegaeri
    real(8) :: EFT, CR
    real(8) :: rhoni, dvexbdr ! CDBM
    real(8) :: xb, fp, fdp, gfun
    real(8) :: AITKEN2P ! function
    real(8), dimension(0:NRMAX) :: pres, ddPhidpsi, tmp
    ! Mainly for derivatives
    real(8), dimension(:), allocatable :: dErdr, dpdr, dErdrS, ErVlc
    real(8), dimension(:), allocatable :: dQdrho, dlnNedrhov
    real(8), dimension(:,:), allocatable :: dTsdV, dPsdV, dNsdrho, dTsdrho

    MDANOMabs = abs(MDANOM)
    MDLNEOL   = mod(MDLNEO,10)
    izvpch = 0

    !     *** Constants ***

    !     Neutral cross section
    !     (NRL Plasma Formulary p52 Eq. (1) (2002))

    Sigma0 = 8.8D-21

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

    pres(:)  = ( PNsV_FIX(:,1)*PTsV_FIX(:,1) &
         &      +PNsV_FIX(:,2)*PTsV_FIX(:,2)) * 1.D20 * rKeV
    IF(MDANOM > 0 .and. maxval(FSANOM) > 0.D0) pres(:)  = 0.5d0 * (pres(:) + pres0(:))

    allocate(ErVlc(0:NRMAX))
    ErVlc(:) = 0.5d0 * (ErV_FIX(:) + ErV0(:))

    IF(PROFM == 0.D0 .AND. FSDFIX(2) /= 0.D0) THEN
       PROFML = (Var(NRA,1)%T * rKeV / (16.D0 * AEE * SQRT(bbt(NRA)))) / FSDFIX(2)
    ELSE
       PROFML = PROFM
    END IF

    IF(PROFC == 0.D0 .AND. FSDFIX(3) /= 0.D0) THEN
       PROFCL = (Var(NRA,1)%T * rKeV / (16.D0 * AEE * SQRT(bbt(NRA)))) / FSDFIX(3)
    ELSE
       PROFCL = PROFC
    END IF

    ! ************** Turbulent transport end **************

    ! Banana width
!    Wbane = (Q(0) * SQRT(RR * amas(1) * Var(0,1)%T * rKeV * amqp) / BphV(0))**(2.D0/3.D0)
!    Wbani = (Q(0) * SQRT(RR * amas(2) * Var(0,2)%T * rKeV * amqp) / (achg(2) * BphV(0)))**(2.D0/3.D0)

    ! Banana width at separatrix
!    PTiVA = Var(NRA,2)%T
    PTiVA = 0.5D0 * Var(0,2)%T
    DBW = 3.D0 * SQRT(PTiVA * rKeV * (amas(2)*amp)) * Q(NRA) / (achg(2) * AEE * BphV(NRA)) &
         & / SQRT(R(NRA) / RR)

    !     *** Calculate derivatives in advance ***
    !     !!! Caution !!!
    !        The r-derivatives of variables, or near-variables (ex. temperature) should be
    !          estimated by their vv-derivatives multiplied by 2*r because they are
    !          evaluated on the vv-abscissa. On the other hand, those of the other
    !          parameters (ex. radial electric field, poloidal magnetic field) should be
    !          directly calculated.

    allocate(dErdr(0:NRMAX),dErdrS(0:NRMAX))
    allocate(dpdr(0:NRMAX))
    dErdr (:) =                 dfdx(R  ,ErVlc,NRMAX,0)
    dpdr  (:) = vro(:) / ravl * dfdx(vv ,pres ,NRMAX,0)

    allocate(dQdrho(0:NRMAX), dlnNedrhov(0:NRMAX))
    allocate(dTsdV(0:NRMAX,NSM), dPsdV(0:NRMAX,NSM))
    do i = 1, NSM
       dTsdV(:,i) = dfdx(vv  ,Var(:,i)%T,NRMAX,0)
       dPsdV(:,i) = dfdx(vv  ,Var(:,i)%p,NRMAX,0)
    end do
    dQdrho    (:) = vro(:) * dfdx(vv  ,Q   ,NRMAX,0)
    dlnNedrhov(:) = dfdx(vv,Var(:,1)%n,NRMAX,0) ! dne/dV, temporarily
    dlnNedrhov(0) = 0.d0                        ! (1/ne)dne/drho_v at rho=0
    do NR = 1, NRMAX
       dlnNedrhov(NR) = 2.d0 * vlt(NR) / ( rhov(NR) * Var(NR,1)%n ) * dlnNedrhov(NR) ! (1/ne)dne/drho_v
    end do

!!D02    write(6,'(F8.5,I4,2F11.6)') T_TX,NRB,Rho(NRB),PT02V(NR)

    !  Smoothing Er gradient for numerical stability
    do NR = 0, NRMAX
       dErdrS(NR) = moving_average(NR,dErdr,NRMAX,NRA)
    end do

    ! *** Temperatures for neutrals ***

    PT01V(:) =   0.5D0 * amas(2) * V0**2 * amqp / rKilo

    !  --- For thermal neutrals originating from slow neutrals ---

    IF(IC == 1) THEN
       ! SCX : source of PN02V
       DO NR = 0, NRMAX
          SCX(NR) = Var(NR,2)%n * SiVcx(Var(NR,2)%T) * PN01V(NR) * 1.D40
       END DO

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
    END IF

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

!    PT02V(:) =   Var(:,2)%T
    PT03V(:) =   Var(:,2)%T

    ! *********************************

    ! Set flag for CDBM model
    model_cdbm = 0
    if( FSCBSH > 0.d0 ) then
       model_cdbm = model_cdbm + 2
    else if( FSCBSH < 0.d0 ) then
       model_cdbm = model_cdbm + 4
    end if
    if( FSCBEL /= 0.d0 ) model_cdbm = model_cdbm + 1

    ! Calculate CDIM coefficient
    RAQPR(:) = vro(:) / ravl * dfdx (vv, rho**4 / Q, NRMAX , 0) ! cf. txcalv.f90 L501

    ! Orbit squeezing effect for neoclassical solvers
    !   gr2phi = psi'(Phi'/psi')' = (dpsi/dV dV/drho)^2 d/dpsi(dPhi/dpsi), 
    !      where rho is an arbitrary radial coordinate.
    !                      ^^^^^^^^^
    !   gr2phi = (dpsi/dV)^2 d/dpsi(dPhi/dpsi) if rho = V is assumed.
    IF( mod(MDOSQZ,10) /= 1 .or. IC == 1 ) THEN
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
    END IF

    !  Coefficients

    L_NR: DO NR = 0, NRMAX

       Vte = SQRT(2.D0 * ABS(Var(NR,1)%T) * rKilo / (amas(1) * amqp))
       Vti = SQRT(2.D0 * ABS(Var(NR,2)%T) * rKilo / (amas(2) * amqp))
       Vtb = SQRT(2.D0 * ABS(Var(NR,2)%T) * rKilo / (amb     * amqp)) ! ??

       PN0tot = (PN01V(NR) + PN02V(NR) + PN03V(NR)) * 1.D20
       SiVizA(NR) = SiViz(Var(NR,1)%T)
       SiVcxA(NR) = SiVcx(Var(NR,2)%T)

       !     *** Ionization rate ***

!old       !     (NRL Plasma Formulary p54 Eq. (12) (2002))
!old       XXX = MAX(Var(NR,1)%T * 1.D3 / EION, 1.D-2)
!old       SiV = 1.D-11 * SQRT(XXX) * EXP(- 1.D0 / XXX) &
!old            &              / (EION**1.5D0 * (6.D0 + XXX))
!old       rNuION(NR) = FSION * SiV * (PN01V(NR) + PN02V(NR)) * 1.D20
       rNuION(NR) = FSION * SiVizA(NR) * PN0tot

       !     *** Slow neutral diffusion coefficient ***
       !  For example,
       !    E.L. Vold et al., NF 32 (1992) 1433

       !  Maxwellian velocity for slow neutrals
       V0ave = sqrt(4.D0 * V0*V0 / Pi)

       !  Total Maxwellian rate coefficients
!!$       if(nr <= NRB) then
!!$          Sitot = (SiVcxA(NRB) * Var(NRB,2)%n + SiVizA(NRB) * Var(NRB,1)%n) *1.D20
!!$       else
          Sitot = (SiVcxA(NR) * Var(NR,2)%n + SiVizA(NR) * Var(NR,1)%n) *1.D20
!!$       end if

       !  Diffusion coefficient for slow neutrals (short m.f.p.)
!old       D01(NR) = FSD01 * V0**2 &
!old            &   / (Sigma0 * (PN01V(NR) * V0  + (Var(NR,2)%n + PN02V(NR)) &
!old            &      * Vti) * 1.D20)
!       D01(NR) = FSD01 * V0ave**2 &
!            &   / (3.D0 * SiVcxA(NR) * Var(NR,2)%n * 1.D20)
       D01(NR) = FSD01 * V0ave*V0ave &
            &  / (3.D0 * (SiVcxA(NR) * Var(NR,2)%n + SiVizA(NR) * Var(NR,1)%n) *1.D20)

       !     *** Thermal neutral diffusion coefficient ***

       !  Maxwellian thermal velocity at the separatrix
!       Viave = sqrt(8.D0 * Var(NR,2)%T * rKilo / (Pi * amas(2) * amqp))
       Viave = sqrt(8.D0 * PT02V(NR) * rKilo / (Pi * amas(2) * amqp))

       !  Mean free path for fast neutrals
       rLmean = Viave / Sitot
!!D02       !  Locally determined mean free path for fast neutrals
!!D02       rLmeanL = sqrt(8.D0 * Var(NR,2)%T * rKilo / (Pi * amas(2) * amqp)) &
!!D02            &  / ((SiVcxA(NR) * Var(NR,2)%n + SiVizA(NR) * Var(NR,1)%n) *1.D20)

       !  Diffusion coefficient for fast neutrals (short to long m.f.p.)
!old       D02(NR) = FSD02 * Vti**2 &
!old            &   / (Sigma0 * (Var(NR,2)%n + PN01V(NR) + PN02V(NR)) &
!old            &      * Vti * 1.D20)
!       D02(NR) = FSD02 * Viave**2 &
!            &   / (3.D0 * SiVcxA(NR) * Var(NR,2)%n * 1.D20)
!!       if(rho(nr) < 0.9d0) then
!!          D02(NR) = D02(NR) * (-0.33333d0*rho(nr)+0.35d0)
!!       else if (rho(nr) >= 0.9d0 .and. rho(nr) <= 1.0d0) then
!!          D02(NR) = D02(NR) * (3.5055d0*rho(nr)**3-2.5055d0)
!!       end if
       D02(NR) = FSD02 * Viave*Viave / (3.D0 * Sitot)
!!!!!!!!!!       if(nr >= NRB .and. nr <= NRA) D02(NR) = D02(NR) * reduce_D02(NR,rLmean)
!!D02       if(nt >= ntmax-1) write(6,'(I3,F10.6,F11.3,3F10.6,1P2E11.3)') nr,r(nr),d02(nr),rLmean,rLmeanL,Var(NR,2)%T,SCX(NR)

       !     *** Halo neutral diffusion coefficient ***

       Viave = sqrt(8.D0 * Var(NR,2)%T * rKilo / (Pi * amas(2) * amqp))
       D03(NR) = FSD03 * Viave*Viave / (3.D0 * Sitot)

       !     *** Charge exchange rate ***
       !  For thermal ions (assuming that energy of deuterium
       !                    is equivalent to that of proton)

!old       !     (Riviere, NF 11 (1971) 363, Eq.(4))
!old       Scxi = 6.937D-19 * (1.D0 - 0.155D0 * LOG10(Var(NR,2)%T*1.D3))**2 &
!old            & / (1.D0 + 0.1112D-14 * (Var(NR,2)%T*1.D3)**3.3d0) ! in m^2
!old       Vave = SQRT(8.D0 * Var(NR,2)%T * rKilo / (PI * amas(2) * amqp))
!old       rNuiCX(NR) = FSCX * Scxi * Vave * (PN01V(NR) + PN02V(NR)) * 1.D20
       rNuiCX(NR)  = FSCX * SiVcxA(NR) * PN0tot
       !  For thermal loss by charge exchange
       rNuiCXT(NR) = FSCX * SiVcxA(NR) * PN01V(NR) * 1.D20

       !  For beam ions
       !     (C.E. Singer et al., Comp. Phys. Comm. 49 (1988) 275, p.323, B163)
!old       Vave = SQRT(8.D0 * Eb * rKilo / (PI * amas(2) * amqp))
!old       rNubCX(NR) = FSCX * Scxb * Vave * (PN01V(NR) + PN02V(NR)) * 1.D20
       rNubCX(NR) = FSCX * Scxb * Vb * PN0tot

       !     *** Collision frequency (with neutral) ***

       rNu0e(NR) = PN0tot * Sigma0 * Vte
       rNu0i(NR) = PN0tot * Sigma0 * Vti
       rNu0b(NR) = PN0tot * Sigma0 * Vtb

       !     *** Collision frequency (90 degree deflection) ***
       !         (see  Helander & Sigmar textbook (2002), p.5, p.79,
       !               Hirshman & Sigmar, p.1105, around (4.7)      )
       !     c.f. Braginskii's collision time: rNue = rNuei, rNui = rNuii / sqrt(2)

       rNuei(NR) = coll_freq(NR,1,2)
       rNuii(NR) = coll_freq(NR,2,2)

       !     *** Collision frequency (energy equipartition) ***
       !     (D. V. Sivukhin, Rev. Plasma Phys. Vol. 4, Consultants Bureau (1966))
       !         Energy relaxation time between two kinds of particles
!approximate formula (amas(1) << amas(2))      rNuTei(NR) = rNuei(NR) * (2.D0 * amas(1) / amas(2))
       rNuTei(NR) = Var(NR,2)%n * 1.D20 * achg(1)*achg(1) * achg(2)*achg(2) * AEE*AEE*AEE*AEE &
            &     * coulog(1.d0,Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amas(1),achg(1),amas(2),achg(2)) &
            &     / (3.D0 * SQRT(2.D0 * PI) * PI * EPS0*EPS0 * amas(1) * amas(2) * amp*amp &
            &     * (  ABS(Var(NR,1)%T)*rKilo / (amas(1) * amqp) + ABS(Var(NR,2)%T)*rKilo / (amas(2) * amqp))**1.5D0)

       BBL = sqrt(bbt(NR))

       Wte = Vte / (Q(NR) * RR) ! Omega_te; transit frequency for electrons
       Wti = Vti / (Q(NR) * RR) ! Omega_ti; transit frequency for ions
       
       EpsL = epst(NR)          ! Inverse aspect ratio
       rNuAsE_inv = EpsL*sqrt(EpsL) * Wte / (SQRT(2.D0) * rNuei(NR))
       rNuAsI_inv = EpsL*sqrt(EpsL) * Wti / (SQRT(2.D0) * rNuii(NR))
       IF(NR /= 0) THEN
          rNuAse(NR) = 1.D0 / rNuAsE_inv
          rNuAsi(NR) = 1.D0 / rNuAsI_inv
       END IF

       if( MDLNEO == 2 .or. MDLNEO > 10 ) then
          !     *** Neoclassical coefficients by NCLASS ***
          !     (W. A. Houlberg, et al., Phys. Plasmas 4 (1997) 3230)

          CALL TX_NCLASS(NR,ETAvar(NR,2),BJBSvar(NR,2), &
               &         ChiNCpe(NR),ChiNCte(NR),ChiNCpi(NR),ChiNCti(NR), &
               &         dTsdV(NR,:),dPsdV(NR,:),gr2phi(NR),IER)
          IF(IER /= 0) IERR = IER
!          write(6,*) NR,chincpe(nr),chincte(nr)
!          write(6,*) NR,chincpi(nr),chincti(nr)
       end if
       if( MDLNEO == 1 .or. MDLNEO > 10 ) then
          !     *** Neoclassical coefficients by Matrix Inversion ***
          !     (M. Kikuchi and M. Azumi, Plasma Phys. Control. Fusion 37 (1995) 1215)

          call tx_matrix_inversion(NR,ETAvar(NR,1),BJBSvar(NR,1), &
               &         ChiNCpe(NR),ChiNCte(NR),ChiNCpi(NR),ChiNCti(NR), &
               &         ddPhidpsi(NR)*MDOSQZN)
!          write(6,*) NR,chincpe(nr),chincte(nr)
!          write(6,*) NR,chincpi(nr),chincti(nr)
       end if

       !     *** Helical neoclassical viscosity ***

!       IF(ABS(FSHL) > 0.D0 .AND. NR > 0) THEN
       IF(ABS(FSHL) > 0.D0 ) THEN
!          IF(int(FSHL) .EQ. 1 ) THEN !! for FSHL = 1 (single Fourier mode)
!!$            Wte = Vte * NCph / RR
!!$            Wti = Vti * NCph / RR
!!$            EpsL = EpsH * rho(NR)**2
!!$            rNuAsE_inv = EpsL**1.5D0 * Wte / (SQRT(2.D0) * rNuei(NR))
!!$            rNuAsI_inv = EpsL**1.5D0 * Wti / (SQRT(2.D0) * rNuii(NR))
!!$            IF(NR.EQ.0) THEN
!!$               BLinv=0.d0
!!$               omegaer=0.d0
!!$            ELSE
!!$!             QL=(Q0-QA)*(1.D0-rho(NR)**2)+QA
!!$!             Bthl = BB*R(NR)/(QL*RR)
!!$!             BLinv=BB/Bthl
!!$               BBL=SQRT(BphV(NR)**2 + BthV(NR)**2)
!!$               BLinv=BBL/BthV(NR)
!!$               omegaer=ErVlc(NR)/(BBL*R(NR))
!!$            ENDIF
!!$            omegaere=EpsL*R(NR) / RR * omegaer**2 / rNuei(NR)**2
!!$!            rNueHL(NR) = FSHL * Wte * BLinv * rNuAsE_inv &
!!$            rNueHL(NR) =        Wte * BLinv * rNuAsE_inv &
!!$            &            /(3.D0+1.67D0*omegaere)
!!$
!!$            omegaeri=EpsL*R(NR) / RR * omegaer**2 / rNuii(NR)**2
!!$!            rNuiHL(NR) = FSHL * Wti * BLinv * rNuAsI_inv &
!!$            rNuiHL(NR) =        Wti * BLinv * rNuAsI_inv &
!!$            &            /(3.D0+1.67D0*omegaeri)
!!$
!!$!          UHth=(RR/NCph)/SQRT((RR/NCph)**2+(R(NR)/NCth)**2)
!!$!          UHph=(R(NR)/NCth)/SQRT((RR/NCph)**2+(R(NR)/NCth)**2)
!!$!          UHth  = real(NCth,8) / NCph
!!$!          UHph  = 1.D0
!!$!    09/02/11 mm
!!$
!!$            UHth=(RR*NCth)/SQRT((RR*NCth)**2+(R(NR)*NCph)**2)
!!$!---- 09/11/26 AF
!!$!          UHph=-(R(NR)*Ncph)/SQRT((RR*NCth)**2+(R(NR)*NCph)**2)
!!$            UHph=-Ncph/SQRT((RR*NCth)**2+(R(NR)*NCph)**2)
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
          IF(NR == 0) THEN
             omegaer=0.d0
          ELSE
             omegaer=ErVlc(NR)/(BBL*R(NR))
          ENDIF
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
             omegaere   = abs(EpsLM(NHFM)) * R(NR) / RR * omegaer**2 / rNuei(NR)**2
             rNueHLM(NHFM, NR) = FSHL * abs(EpsLM(NHFM))**1.5D0 * Var(NR,1)%T * rKilo          &
                  &            / (amas(1) * amqp * RR**2 * rNuei(NR) * (3.D0 + 1.67D0 * omegaere)) ! [s^-1]
!
             omegaeri   = abs(EpsLM(NHFM)) * R(NR) / RR * omegaer**2 / rNuii(NR)**2
             rNuiHLM(NHFM, NR) = FSHL * abs(EpsLM(NHFM))**1.5D0 * Var(NR,2)%T * rKilo          &
                  &            / (amas(2) * amqp * RR**2 * rNuii(NR) * (3.D0 + 1.67D0 * omegaeri)) ! [s^-1]
!
             if (UHphSwitch == 0) then    ! For the case which (m=0, n>0) component DOES NOT exist
                UHth =    RR * HPN(NHFM,1) / SQRT((RR * HPN(NHFM,1))**2 + (R(NR) * HPN(NHFM,2))**2) ! [nondimensional]
                UHph =         HPN(NHFM,2) / SQRT((RR * HPN(NHFM,1))**2 + (R(NR) * HPN(NHFM,2))**2) ! [m^-1]
             else if (HPN(NHFM,1) == 0 .and. HPN(NHFM,2) > 0) then ! to avoid NaN and Infty for (m=0, n>0) component
                UHth = 0.d0 ! [nondimensional]
                UHph = 1.d0 ! [nondimensional]
             else   ! For the case which (m=0, n>0) component exist
                UHth =    RR * HPN(NHFM,1) / SQRT((RR * HPN(NHFM,1))**2 + (R(NR) * HPN(NHFM,2))**2) ! [nondimensional]
                UHph = R(NR) * HPN(NHFM,2) / SQRT((RR * HPN(NHFM,1))**2 + (R(NR) * HPN(NHFM,2))**2) ! [nondimensional]
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
             
          ENDDO
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
       ELSE
          rNueHLthth(:) = 0.D0
          rNueHLthph(:) = 0.D0
          rNueHLphth(:) = 0.D0
          rNueHLphph(:) = 0.D0
          rNuiHLthth(:) = 0.D0
          rNuiHLthph(:) = 0.D0
          rNuiHLphth(:) = 0.D0
          rNuiHLphph(:) = 0.D0
       END IF
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
       Alpha(NR) = - Q(NR)*Q(NR) * RR * dpdr(NR) * 2.D0 * rMU0 / bbt(NR)

       !   ***** Heat pinch *****

       Vhps(NR,1:NSM) = 0.d0

       !   ***** Thermal diffusivity *****

       !   *** CDBM model ***

       IF (maxval(FSANOM) > 0.D0) THEN

          !   *** CDBM model ***
          IF (MDANOMabs == 1) THEN
             rhoni = Var(NR,2)%n * 1.D20 * (amas(2)*amp)
             dvexbdr = dErdrS(NR) / bbrt(NR)
!             dvexbdr = dErdr(NR) / bbrt(NR)

             ! Magnetic curvature
             rKappa(NR) = FSCBKP * (- epst(NR) * (1.D0 - 1.D0 / Q(NR)**2))

             call cdbm(BBL,rr,r(NR),elip(NR),Q(NR),S(NR),Var(NR,1)%n*1.d20,rhoni,dpdr(NR), &
                  &    dvexbdr,FSCBAL,FSCBKP,abs(FSCBSH)*rG1,model_cdbm,Dturb, &
                  &    FCDBM(NR),rKappa(NR),rG1h2(NR))

          !   *** CDIM model ***
          ELSE IF (MDANOMabs == 2) THEN
             ! Alfven velocity
             Va = SQRT(BBL*BBL / (rMU0 * Var(NR,2)%n * 1.D20 * (amas(2)*amp)))
             ! Squared plasma frequency
             Wpe2 = Var(NR,1)%n * 1.D20 * AEE / (amas(1) * amqp * EPS0)

             ! Arbitrary coefficient for CDIM model
             rGCIM = 10.D0
             OMEGAPR = (RAVL / RR)**2 * (real(NCph,8) / NCth) * RAQPR(NR)
         
             IF(NR == 0) THEN  ! for s=0
                FCDIM(NR) = 0
             ELSE
                FCDIM(NR) = 3.D0 * (0.5D0 * OMEGAPR)**1.5D0 * (RR / RAVL)**1.5D0 / (Q(NR) * S(NR)**2)
             END IF

             ! ExB rotational shear
             IF(NR == 0) THEN
                rGIM = 0.D0
                rHIM = 0.D0
             ELSE
!                rGIM = rG1
                rGIM = 1.04D0 * R(NR)*R(NR) / ( dpdr(NR) * 2.D0 * rMU0 / bbt(NR) &
                     &                         * OMEGAPR * RAVL*RAVL * RR*RR )
                rHIM = RAVL * SQRT( rMU0 * (amas(2)*amp) * Var(NR,2)%n * 1.D20 ) / BthV(NR) * dErdrS(NR) / BBL
             END IF

             ! Turbulence suppression by ExB shear for CDIM mode
             rG1h2IM(NR) = 1.D0 / (1.D0 + abs(FSCBSH) * (rGIM * rHIM**2))
             ! Turbulent transport coefficient calculated by CDIM model
             Dturb = rGCIM * FCDIM(NR) * rG1h2IM(NR) * ABS(Alpha(NR))**1.5D0 &
                  &              * VC*VC / Wpe2 * Va / (Q(NR) * RR)

!-----------------memo:beta'=dpdr(NR) * 2.D0 * rMU0 / bbt(NR)-----------

!          IF(NR == 0 .OR. NR == 1) & 
!                         write(6,*) ',NR=',NR,'RAQPR=',RAQPR(NR),'OMEGAPR=',OMEGAPR, &
!               &     'Q=',Q(NR),'S=',S(NR),'FCDIM=',FCDIM(NR),'Dturb=',Dturb

          END IF
          IF(Rho(NR) == 1.D0) DturbA = Dturb 
       ELSE
          rG1h2(NR)   = 0.D0
          FCDBM(NR)   = 0.D0
          rG1h2IM(NR) = 0.D0
          FCDIM(NR)   = 0.D0
          Dturb       = 0.D0
       END IF

       !     *** Turbulent transport of particles ***
       !     ***     Wave-particle interaction    ***

       RhoSOL = 1.D0
       IF (RHO(NR) < RhoSOL) THEN
!          DeL = diff_prof(RHO(NR),FSDFIX(1),PROFD,PROFD1,PROFDB) + FSANOM(1) * Dturb
          DeL = diff_prof(RHO(NR),FSDFIX(1),PROFD,PROFD1,PROFDB,PROFD2,0.99d0,0.07d0) &
!          DeL = diff_prof(RHO(NR),FSDFIX(1),PROFD,PROFD1,PROFDB,PROFD2,0.99d0,0.125d0) &
               & + FSANOM(1) * Dturb
       ELSE
          IF(FSPCLD == 0.D0) THEN
             ! Bohm-like diffusivity
             factor_bohm = (FSDFIX(1) * PROFD + PROFDB + FSANOM(1) * Dturb) &
                  &  / (Var(NRA,1)%T * rKeV / (16.D0 * AEE * SQRT(bbt(NRA))))
             DeL = factor_bohm * Var(NR,1)%T * rKeV / (16.D0 * AEE * BBL)
          ELSE
             IF(FSANOM(1) == 0.D0) THEN
                ! Fixed value and fixed profile
                DeL = FSPCLD * diff_prof(RhoSOL,FSDFIX(1),PROFD,PROFD1,PROFDB,0.d0)
             ELSE
                ! Theory-based anomalous diffusivity
                DeL = FSANOM(1) * DturbA
             END IF
          END IF
!!$          DeL = FSDFIX(1) * PROFD + FSANOM(1) * Dturb
       END IF
       ! Particle diffusivity
!!$       if(rho(nr) > 0.7d0) then
!!$          DeL = DeL * (-5.d0/3.d0*Rho(NR)+13.d0/6.d0)
!!$!          DeL = DeL * (50.d0/9.d0*(Rho(NR)-1.d0)**2+0.5d0)
!!$       end if
       De(NR)   = De0   * DeL
       Di(NR)   = Di0   * DeL

       ! Turbulent pinch term [1/Wb]
       VWpch(NR) = VWpch0 * 0.d0

       !     *** Turbulent transport of momentum ***

       RhoSOL = 1.D0

       IF (RHO(NR) < RhoSOL) THEN
          DeL = diff_prof(RHO(NR),FSDFIX(2),PROFML,PROFM1,PROFMB,0.d0) + FSANOM(2) * Dturb
!pedestal          if(rho(nr) > 0.9d0) DeL = DeL * exp(-120.d0*(rho(nr)-0.9d0)**2)
       ELSE
          IF(FSPCLM == 0.D0) THEN
             factor_bohm = (FSDFIX(2) * PROFML + PROFMB + FSANOM(2) * Dturb) &
                  &  / (Var(NRA,1)%T * rKeV / (16.D0 * AEE * SQRT(bbt(NRA))))
!bohm_model2             DeL =  (1.D0 - MOD(FSBOHM,2.D0)) * FSDFIX(2) * PROFML &
!bohm_model2                  &+ FSBOHM * Var(NR,1)%T * rKeV / (16.D0 * AEE * BBL)
             DeL = factor_bohm * Var(NR,1)%T * rKeV / (16.D0 * AEE * BBL)
          ELSE
             IF(FSANOM(2) == 0.D0) THEN
                DeL = FSPCLM * diff_prof(RhoSOL,FSDFIX(2),PROFML,PROFM1,PROFMB,0.d0)
             ELSE
                DeL = FSANOM(2) * DturbA
             END IF
!pedestal             DeL = FSPCLM * FSDFIX(2) * PROFML * exp(-120.d0*(rho(nra)-0.9d0)**2)
          END IF
       END IF
       ! Viscosity
       rMue(NR) = rMue0 * DeL
       rMui(NR) = rMui0 * DeL

       !   ***** Momentum pinch *****
       !   -R Vmpi/rMui = 1.1 R/Lne + 1.0
       !   [T. Tala et al., 2012, Proc of 24th IAEA FEC (San Diego) ITR/P1-1]

       if(NR /= 0) then
!          zvpch = - rMui(NR) / rr * (1.1d0 * rr * moving_average(NR,dlnNedrhov,NRMAX) + 1.d0) ! smoothing 1/Lne
          zvpch = rMui(NR) / rr * (1.1d0 * rr * moving_average(NR,dlnNedrhov,NRMAX) + 1.d0) ! smoothing 1/Lne
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

       IF (RHO(NR) < RhoSOL) THEN
          DeL = diff_prof(RHO(NR),FSDFIX(3),PROFCL,PROFC1,PROFCB,0.d0) + FSANOM(3) * Dturb
!pedestal          if(rho(nr) > 0.9d0) DeL = DeL * exp(-120.d0*(rho(nr)-0.9d0)**2)
       ELSE
          IF(FSPCLC == 0.D0) THEN
             factor_bohm = (FSDFIX(3) * PROFCL + PROFCB + FSANOM(3) * Dturb) &
                  &  / (Var(NRA,1)%T * rKeV / (16.D0 * AEE * SQRT(bbt(NRA))))
!bohm_model2             DeL =  (1.D0 - MOD(FSBOHM,2.D0)) * FSDFIX(3) * PROFCL &
!bohm_model2                  &+ FSBOHM * Var(NR,1)%T * rKeV / (16.D0 * AEE * BBL)
             DeL = factor_bohm * Var(NR,1)%T * rKeV / (16.D0 * AEE * BBL)
          ELSE
             IF(FSANOM(3) == 0.D0) THEN
                DeL = FSPCLC * diff_prof(RhoSOL,FSDFIX(3),PROFCL,PROFC1,PROFCB,0.d0)
             ELSE
                DeL = FSANOM(3) * DturbA
             END IF
!pedestal             DeL = FSPCLC * FSDFIX(3) * PROFCL * exp(-120.d0*(rho(nra)-0.9d0)**2)
          END IF
       END IF
!       DeL = 3.d0
       ! Thermal diffusivity
       Chie(NR) = Chie0 * DeL
       Chii(NR) = Chii0 * DeL

!!$       ! <omega/m>
!!$       WPM(NR) = WPM0 * Var(NR,1)%T * rKeV / (RAVL**2 * AEE * BphV(NR))
       ! Ad hoc turbulent pinch velocity
       IF(NR == 0) THEN
          FVpch(NR) = 0.D0
       ELSE
          FVpch(NR) = Var(NR,1)%T * rKilo * sdt(NR) * VWpch(NR)
       END IF

       ! Coefficient of the force causing the quasilinear flux, induced by drift wave
       FQLcoef(NR) = sst(NR) * sdt(NR) * De(NR) / (Var(NR,1)%T * rKilo)
       if(NR /= 0) then
          FQLcoef1(NR) = FQLcoef(NR) * sdt(NR) * bbt  (NR) / bri(NR)
          FQLcoef2(NR) = FQLcoef(NR) * sdt(NR) * fipol(NR) / bri(NR)
       else
          FQLcoef1(0) = bbt(0)   * De(0) / (Var(0,1)%T * rKilo)
          FQLcoef2(0) = fipol(0) * De(0) / (Var(0,1)%T * rKilo)
       end if

       !     *** Loss to divertor ***

!       IF (R(NR) + DBW > RAVL) THEN
       IF (rho(NR) > 1.d0) THEN
!          Cs = SQRT(2.D0 * Var(NR,1)%T * rKilo / amas(2) / amqp)
          Cs = SQRT((achg(2) * Var(NR,1)%T + 3.D0 * Var(NR,2)%T) * rKilo / (amas(2) * amqp))
          Lc = 2.D0 * PI * Q(NR) * RR ! Connection length to the divertor
          RL = (R(NR) - RAVL) / DBW! / 2.D0
          rNuL  (NR) = FSLP  * Cs / Lc &
               &             * RL*RL / (1.D0 + RL*RL)
          ! Classical heat conduction [s**4/(kg**2.5*m**6)]
          ! (C S Pitcher and P C Stangeby, PPCF 39 (1997) 779)
          Chicl = (4.D0*PI*EPS0)**2 &
               & /(SQRT(amas(1)*amp)*coulog(1.d0,Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amas(1),achg(1),amas(2),achg(2))*AEE**4*Zeff)

          ! When calculating rNuLTe, we fix Var(:,1)%n and Var(:,1)%T constant during iteration
          !   to obain good convergence.
          rNuLTe(NR) = FSLTE * Chicl * (PTsV_FIX(NR,1)*rKeV)**2.5D0 &
               &                  /(Lc**2 * PNsV_FIX(NR,1)*1.D20) &
               &             * RL*RL / (1.D0 + RL*RL)
          rNuLTi(NR) = FSLTI * Cs / Lc &
               &             * RL*RL / (1.D0 + RL*RL)
!!$          IF(ABS(FSRP) > 0.D0) THEN
             UbparaL = BUbparV(NR) / BBL
!             IF(NR == NRMAX) Ubpara(NR) = AITKEN2P(R(NRMAX), &
!                  & Ubpara(NRMAX-1),Ubpara(NRMAX-2),Ubpara(NRMAX-3),&
!                  & R(NRMAX-1),R(NRMAX-2),R(NRMAX-3))
             UbparaL = max(UbparaL, FSLP*Cs)
             rNuLB(NR) = UbparaL / Lc &
                  &          * RL*RL / (1.D0 + RL*RL)
!!$          END IF
       ELSE
          rNuL(NR) = 0.D0
          rNuLTe(NR) = 0.D0
          rNuLTi(NR) = 0.D0
          rNuLB(NR) = 0.D0
       END IF

       !     *** Current density profile ***

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

       PEQe(NR)  = - 1.5D0 * rNuTei(NR) * Var(NR,1)%n * 1.D20 * (Var(NR,1)%T - Var(NR,2)%T) * rKeV
       PEQi(NR)  = - 1.5D0 * rNuTei(NR) * Var(NR,1)%n * 1.D20 * (Var(NR,2)%T - Var(NR,1)%T) * rKeV

       !     *** Ohmic heating power ***
       ! POH = <j.E> = -PsitdotV * sum_s (e_s n_s Var%Uthhat) + PsidotV * sum_s (e_s n_s Var%UphR)
       POHe(NR) = achg(1) * aee * Var(NR,1)%n * 1.d20 &
            &   * ( - PsitdotV(NR) * Var(NR,1)%Uthhat + PsidotV(NR) * Var(NR,1)%UphR )
       POHi(NR) = achg(2) * aee * Var(NR,2)%n * 1.d20 &
            &   * ( - PsitdotV(NR) * Var(NR,2)%Uthhat + PsidotV(NR) * Var(NR,2)%UphR )
       POH(NR)  = POHe(NR) + POHi(NR)

       !     *** Bremsstraulung loss ***
       !     (NRL Plasma Formulary p57 Eq. (30) (2002))

       PBr(NR) = 5.35D-37 * achg(2)*achg(2) * Var(NR,1)%n * Var(NR,2)%n * 1.D40 * SQRT(Var(NR,1)%T)

       !     *** Particle diffusion due to magnetic braiding ***AF 2008-06-08

       IF (R(NR) > RMAGMN .AND. R(NR) < RMAGMX) THEN
          DMAG(NR)=DMAG0*16.D0*(R(NR)-RMAGMN)**2*(RMAGMX-R(NR))**2 &
          &                   /(RMAGMX-RMAGMN)**4
          DMAGe(NR)=DMAG(NR)*Vte
          DMAGi(NR)=DMAG(NR)*Vti
       ELSE
          DMAG(NR)=0.D0
          DMAGe(NR)=0.D0
          DMAGi(NR)=0.D0
       ENDIF

    END DO L_NR

    rNuAse(0) = AITKEN2P(R(0),rNuAse(1),rNuAse(2),rNuAse(3),R(1),R(2),R(3))
    rNuAsi(0) = AITKEN2P(R(0),rNuAsi(1),rNuAsi(2),rNuAsi(3),R(1),R(2),R(3))

    !  Linear extrapolation

!!    ChiNCpe(0) = 2.D0 * ChiNCpe(0) - ChiNCpe(1)
!!    ChiNCte(0) = 2.D0 * ChiNCte(0) - ChiNCte(1)
!!    ChiNCpi(0) = 2.D0 * ChiNCpi(0) - ChiNCpi(1)
!!    ChiNCti(0) = 2.D0 * ChiNCti(0) - ChiNCti(1)
    ChiNCpe(0) = ChiNCpe(1)
    ChiNCte(0) = ChiNCte(1)
    ChiNCpi(0) = ChiNCpi(1)
    ChiNCti(0) = ChiNCti(1)
!    ETAvar(0,2)  = 2.D0 * ETAvar(0,2)  - ETAvar(1,2)
!    BJBSvar(0,2) = 2.D0 * BJBSvar(0,2) - BJBSvar(1,2)
!!$    ! For Neumann condition, finite viscosity is required at the magnetic axis.
!!$    De(0)    = De(1)
!!$    Di(0)    = Di(1)
!!$    rMue(0)  = rMue(1)
!!$    rMui(0)  = rMui(1)

    !     *** ExB shearing rate (Hahm & Burrel PoP 1995) ***
    ! omega_ExB = |- <R^2B_p^2>/B d/dpsi(dPhi/dpsi)|
    wexb(:)  = abs(- sst(:) * sdt(:) / BB * ddPhidpsi(:))

    !     *** Linear growth rate for toroidal gamma_etai branch of the ITG mode ***
    !        (F.Crisanti et al, NF 41 (2001) 883)
    allocate(dNsdrho(0:NRMAX,NSM), dTsdrho(0:NRMAX,NSM))
    do i = 1, NSM
       dNsdrho(:,i) = vro(:) * dfdx(vv  ,Var(:,i)%n,NRMAX,0)
       dTsdrho(:,i) = vro(:) * dTsdV(:,i)
    end do
    gamITG(:,1) = 0.1d0 * sqrt(Var(:,1)%T*rKilo/(amas(2)*amqp))/RAVL * sqrt(RAVL/RR) &
         &            * sqrt(abs(dNsdrho(:,2))/Var(:,2)%n &
         &                 + abs(dTsdrho(:,2))/Var(:,2)%T) &
         &            * sqrt(Var(:,2)%T/Var(:,1)%T)

    gamITG(0,2:3) = 0.d0
    i = 0
    do nr = 1, nrmax
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
          gamITG(NR,3) = kthrhos * sqrt(Var(NR,1)%T*rKilo/(amas(2)*amqp)) / RAVL &
               &       * sqrt(2.d0 * RAVL / RR * (1.d0 / Ln + 1.d0 / LT))
       end if
    end do

    if(MDANOMabs == 3) then
       call txmmm95(dNsdrho(:,1),dNsdrho(:,2),dTsdrho(:,1),dTsdrho(:,1),dQdrho,FSCBSH)
    end if

    !     *** ETB model ***

    IF(MDLETB /= 0) THEN ! ETB on
       Frdc = 0.1d0
       if(Fcoef > Frdc) then
          Dcoef = (1.d0 - Frdc) / NTMAX
          Fcoef = 1.d0 - NT * Dcoef
       end if
       DO NR = 0, NRA
          IF(RhoETB(1) /= 0.D0) THEN
             IF(Rho(NR) > RhoETB(1)) THEN
                De(NR) = De(NR) * Fcoef
                Di(NR) = Di(NR) * Fcoef
             END IF
          END IF
          IF(RhoETB(2) /= 0.D0) THEN
             IF(Rho(NR) > RhoETB(2)) THEN
                rMue(NR) = rMue(NR) * Fcoef
                rMui(NR) = rMui(NR) * Fcoef
             END IF
          END IF
          IF(RhoETB(3) /= 0.D0) THEN
             IF(Rho(NR) > RhoETB(3)) THEN
                Chie(NR) = Chie(NR) * Fcoef
                Chii(NR) = Chii(NR) * Fcoef
             END IF
          END IF
       END DO
    END IF

    !     *** Resistivity and current density ***

    do NR = 0, NRMAX
        ! +++ Hirshman, Hawryluk and Birge model +++
       Vte = SQRT(2.D0 * ABS(Var(NR,1)%T) * rKilo / amas(1) * amqp)
       Wte = Vte / (Q(NR) * RR) ! Omega_te; transit frequency for electrons
       rNuAsE_inv = epst(NR)*sqrt(epst(NR)) * Wte / (SQRT(2.D0) * rNuei(NR))
       EFT  = ft(NR) * rNuAsE_inv / (rNuAsE_inv + (0.58D0 + 0.2D0 * Zeff))
       CR   = 0.56D0 * (3.D0 - Zeff) / ((3.D0 + Zeff) * Zeff)
       ! Spitzer resistivity for hydrogen plasma (parallel direction)
       ETAS(NR) = CORR(1.D0) * amas(1) * amqp * rNuei(NR) / (Var(NR,1)%n * 1.D20 * AEE)
       ETAvar(NR,4) = ETAS(NR) * Zeff * (1.D0 + 0.27D0 * (Zeff - 1.D0)) &
            &   /((1.D0 - EFT) * (1.D0 - CR * EFT) * (1.D0 + 0.47D0 * (Zeff - 1.D0)))
    end do

    if( MDLNEO > 10 ) then
       !     *** Bootstrap current and resistivity by Sauter model ***
       !     (O. Sauter, et al.,  Phys. Plasmas 6 (1999) 2834, ibid. 9 (2002) 5140)
       NR = 0
       ETAvar(NR,3) = ETAS(NR) * (CORR(Zeff) / CORR(1.D0))
       BJBSvar(NR,3) = 0.D0

       do NR = 1, NRMAX
          CALL SAUTER(Var(NR,1)%n,Var(NR,1)%T,dTsdV(NR,1),dPsdV(NR,1), &
               &      Var(NR,2)%n,Var(NR,2)%T,dTsdV(NR,2),dPsdV(NR,2), &
               &      Q(NR),sdt(NR),fipol(NR),epst(NR),RR,achg(2),Zeff,ft(nr), &
               &      rlnLei_IN=coulog(1.d0,Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T, &
               &                       amas(1),achg(1),amas(2),achg(2)), &
               &      rlnLii_IN=coulog(1.d0,Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T, &
               &                       amas(2),achg(2),amas(2),achg(2)), &
               &      BJBS=BJBSvar(NR,3),ETA=ETAvar(NR,3))
       end do
    end if

    IF(FSNC /= 0) THEN
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
    ELSE
       ! Spitzer resistivity when no neoclassical effects
       ETA(:) = ETAS(:)
    END IF

    if( MDBSETA == 0 ) then
       do NR = 0, NRMAX
          ! Bootstrap current density estimated by the neoclassical transport model
          BJBS(NR) = BJBSvar(NR,MDLNEOL)

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
       ETAvar(NR,0) = BEpara(NR) / (BJPARA(NR) - BJBSvar(NR,MDLNEOL) - BJNB(NR))
    end do

    !     ***** Ion Orbit Loss *****

    SiLC  (:) = 0.D0
    SiLCB (:) = 0.D0
    SiLCph(:) = 0.D0
    rNuOL (:) = 0.D0
    IF (ABS(FSLC) > 0.D0) THEN
       IF(MDLC == 1) THEN
          ! Ref. [K. C. Shaing, Phys. Fluids B 4 (1992) 3310]
          do NR = 1, NRMAX
             Vti = SQRT(2.D0 * Var(NR,2)%T * rKilo / (amas(2) * amqp))
             ! Orbit squeezing factor, given just below Eq. (3) of Ref.
             sqz = 1.d0 + (fipol(NR) / bb)**2 * amqp * amas(2) / achg(2) * abs(ddPhidpsi(NR))

             EpsL = epst(NR)
             BBL = sqrt(bbt(NR))

             ! rNuDL : deflection collisional frequency at V = Vti, see [Kikuchi PPCF 1995 p.1236]
             rNuDL = 0.d0
             do i = 1, NSM
                xb    = Vti / SQRT(2.D0 * Var(NR,i)%T * rKilo / (amas(i) * amqp))
                fp    = erf(xb)
                fdp   = 2.d0/sqrt(pi)*exp(-xb**2) ! derivative of erf
                gfun  = (fp-xb*fdp)/(2.d0*xb**2)  ! Chandrasekhar function
                rNuDL = rNuDL + 0.75d0 * sqrt(pi) &
                     & * (Var(NR,i)%n * achg(i)**2) / (Var(NR,2)%n * achg(2)**2) &
                     & * (fp - gfun) ! Actually (fp-gfun)/xa**3, but xa at V=Vti is equal to unity.
             end do
             rNuDL = rNuDL * coll_freq(NR,2,2)

             rNustar = rNuDL * RR * Q(NR) / (Vti * (abs(sqz) * EpsL)**1.5D0)

             rNuOLL  = 2.25D0 * rNuDL / (sqrt(PI) * sqrt(2.D0 * abs(sqz) * EpsL)) &
                  &  * EXP(-(rNustar**0.25D0 + (achg(2) * BBL / (amas(2) * amqp)) &
                  &  * sqrt(abs(sqz)) / (fipol(NR) * Vti) &
                  &  * ABS(PsiV(NR) - PsiV(NRA)) / sqrt(2.D0 * EpsL))**2)
             IF(FSLC <= 1.D0) THEN
                rNuOL(NR) = FSLC * rNuOLL
             ELSE
                SiLC  (NR) = - mod(FSLC,1.d1) * Var(NR,2)%n * rNuOLL
                SiLCB (NR) = SiLC(NR) * Var(NR,2)%BUpar
                SiLCph(NR) = SiLC(NR) * Var(NR,2)%RUph
             END IF
          end do

          do NR = 0, NRMAX
             tmp(NR) = moving_average(NR,rNuOL,NRMAX)
          end do
          rNuOL(:) = tmp(:)

       ELSEIF(MDLC == 2) THEN
          !     S. -I. Itoh and K. Itoh, Nucl. Fusion 29 (1989) 1031
          IF(FSLC <= 1.D0) THEN
             ! RLOSS : Numerical coefficient proportional to the relative number of ions
             !         in the loss cone in velocity space
             RLOSS = 0.1D0
             rNuOL(0) = 0.D0
             DO NR = 1, NRMAX
                EpsL = epst(NR)
                Vti = SQRT(Var(NR,2)%T * rKilo / (amas(2) * amqp))
                RhoIT = Vti * amas(2) * amqp / (achg(2) * BthV(NR))
                RL = (R(NR) - (RAVL - 1.5D0 * RhoIT)) / DBW ! Alleviation factor
                IF(R(NR) > (RAVL - RhoIT)) THEN
!                IF(ABS(RAVL - R(NR)) <= RhoIT .AND. RHO(NR) < 1.D0) THEN
                   ExpArg = -2.D0 * EpsL * (ErVlc(NR) / BthV(NR))**2 / Vti**2
                   ExpArg = ExpArg * (R(NR) / RAVL)**2
                   rNuOL(NR) = FSLC * RLOSS * rNuii(NR) / SQRT(EpsL) * EXP(ExpArg) &
                        &    * RL**2 / (1.D0 + RL**2)
                ELSE
                   rNuOL(NR) = 0.D0
                END IF
             END DO
          ELSE
             DO NR = 1, NRA
                EpsL = epst(NR)
                Vti = SQRT(Var(NR,2)%T * rKilo / (amas(2) * amqp))
                RhoIT = Vti * amas(2) * amqp / (achg(2) * BthV(NR))
                RhoIT = MIN(RhoIT,0.1D0)
                Wti = Vti / (Q(NR) * RR)
                rNuAsI_inv = EpsL**1.5D0 * Wti / (SQRT(2.D0) * rNuii(NR))
                ExpArg = 2.D0 * EpsL / Vti**2 * (ErVlc(NR) / BthV(NR))**2
                AiP = rNuii(NR) * SQRT(EpsL) * rNuAsI_inv / (1.D0 + rNuAsI_inv) &
                     & * EXP(- ExpArg)
                DO NR1 = NRA, NRMAX
                   DISTAN = (R(NR1) - R(NR)) / RhoIT
                   SiLCL = AiP * EXP( - DISTAN**2) * Var(NR,2)%n
                   SiLC(NR) = SiLC(NR) - SiLCL
                   SiLC(NR1) = SiLC(NR1) + SiLCL * R(NR) / R(NR1)

                   SiLCBL = SiLCL * (amas(2)*amp) * Var(NR,2)%Uth * R(NR)
                   SiLCB(NR) = SiLCB(NR) - SiLCBL
                   SiLCB(NR1) = SiLCB(NR1) + SiLCBL * R(NR) / R(NR1)

                   SiLCphL = SiLCL * (amas(2)*amp) * Var(NR,2)%Uph
                   SiLCph(NR) = SiLCph(NR) - SiLCphL
                   SiLCph(NR1) = SiLCph(NR1) + SiLCphL * R(NR) / R(NR1)
                END DO
             END DO

             SiLC  (:) = FSLC * SiLC  (:)
             SiLCB (:) = FSLC * SiLCB (:)
             SiLCph(:) = FSLC * SiLCph(:)
          END IF
       END IF
    END IF

    !     ***** Toroidal ripple effect *****
    call ripple_effect(dQdrho)

!    !     ***** Neoclassical toroidal viscosity (NTV) *****
!    !      "rNuNTV" and "UastNC" are obtained from NTVcalc
!    
!    CALL NTVcalc
!    rNuNTV(:) = 0.D0
!    UastNC(:) = 0.D0

    deallocate(dErdr,dpdr,dErdrS,ErVlc)
    deallocate(dQdrho,dlnNedrhov)
    deallocate(dTsdV,dTsdrho,dPsdV,dNsdrho)

    RETURN
  END SUBROUTINE TXCALC

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
!     valid for 1eV => 10^5eV
! 
!     (C.E. Singer et al., Comp. Phys. Comm. 49 (1988) 275, Eq. (2.9.5l))
!     (R.L. Freeman and E.M. Jones, Culham Laboratory Report, CLM-R-137 (May 1974)
!
!     Inputs (real*8): tekev : Electron temperature [keV]
!     Output (real*8): SiViz : Ionization maxwellian rate coefficient [m^3/s]
!
!     Internal (real*8): a : Table of Maxwellian rate coefficients in TABLE 3
!                                     ^^^^^^^^^^
!***************************************************************

  real(8) function SiViz(tekev)

    real(8), intent(in) :: tekev
    real(8) :: x, tekev_temp
    real(8), dimension(0:6) :: a
    data a /-0.3173850D02, 0.1143818D02, -0.3833998D01, 0.7046692D0, &
         &  -0.7431486D-1, 0.4153749D-2, -0.9486967D-4/

    if(tekev < 1.d-3) then
       write(6,'(A,1pE12.4)') &
            'Function SiViz: Out of energy range. tekev=', tekev
       tekev_temp=1.D-3
    else if(tekev > 1.d2) then
       write(6,'(A,1pE12.4)') &
            'Function SiViz: Out of energy range. tekev=', tekev
       tekev_temp=1.D2
    else
       tekev_temp=tekev
    endif
    x = log(tekev_temp * 1.d3)
    SiViz = exp(a(0)+x*(a(1)+x*(a(2)+x*(a(3)+x*(a(4)+x*(a(5)+x*a(6)))))))*1.D-6

  end function SiViz

!***************************************************************
!
!   Rate coefficients for charge exchange cross-section protons on atomic hydrogen
!     valid for 1eV => 10^5eV
! 
!     (R.L. Freeman and E.M. Jones, Culham Laboratory Report, CLM-R-137 (May 1974)
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
    data a /-0.1841757D02, 0.5282950D0, -0.2200477D0,   0.9750192D-1, &
         &  -0.1749183D-1, 0.4954298D-3, 0.2174910D-3, -0.2530205D-4, 0.8230751D-6/

    if(tikev < 1.d-3) then
       write(6,'(A,1pE12.4)') &
            'Function SiVcx: Out of energy range. tikev=', tikev
       tikev_temp=1.D-3
    else if(tikev > 1.d2) then
       write(6,'(A,1pE12.4)') &
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
!     Internal (real*8): a : Table of Maxwellian rate coefficients in TABLE 3
!                                     ^^^^^^^^^^
!***************************************************************

  real(8) function reduce_D02(NRctr,rLmean)
    use tx_commons, only : Pi, NRMAX, R
    integer(4), intent(in) :: NRctr
    real(8), intent(in) :: rLmean

    integer(4) :: nr, nr0, idebug = 0
    real(8) :: Rctr, costh0, theta0, frac, DltL, costh, theta, theta1, rLmean_eff, rLmean_av

    Rctr = R(NRctr)
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
       if(R(nr) > DltL) then
          nr0 = nr
          exit
       end if
    end do

    rLmean_av = 0.d0
    theta     = 0.d0
    do nr = nr0, nrctr
       costh  = (Rctr**2 + rLmean**2 - R(nr)**2) / (2.d0 * Rctr * rLmean)
       theta1 = 2.d0 * acos(costh)
       theta  = theta1 - theta
!       rLmean_eff = rLmean * costh
       rLmean_eff = Rctr - R(nr)
       rLmean_av  = rLmean_av + rLmean_eff * (theta / theta0)
       if(idebug /= 0) write(6,'(I3,5F15.7)') nr,theta1,theta,theta/theta0,rLmean_eff,rLmean_av
       theta  = theta1
    end do
 
    reduce_D02 = frac * (rLmean_av / rLmean)
    if(idebug /= 0) write(6,*) "reduce_D02=",reduce_D02

  end function reduce_D02

end module tx_variables
