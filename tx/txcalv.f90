!     $Id$
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
    use tx_interface, only : dfdx
    use aux_system, only : Vbabsmax
    REAL(8), DIMENSION(1:NQM,0:NRMAX), INTENT(IN) :: XL
    integer(4), intent(in), optional :: ID
    INTEGER(4) :: NR
    real(8) :: BBL
    real(8) :: FCTR

    IF(present(ID)) THEN
       ! The pres0 and ErV0 are the values evaluated at the previous time step
       !   for numerical stability when using turbulent transport models.
       IF(MDFIXT == 0) THEN
          pres0(0:NRMAX) = (XL(LQe5,0:NRMAX) + XL(LQi5,0:NRMAX)) * 1.D20 * rKeV
       ELSE
          pres0(0:NRMAX) = (  XL(LQe1,0:NRMAX) * XL(LQe5,0:NRMAX) &
               &            + XL(LQi1,0:NRMAX) * XL(LQi5,0:NRMAX)) * 1.D20 * rKeV
       END IF
       ErV0 (0:NRMAX) = - 2.D0 * R(0:NRMAX) * dfdx(PSI,XL(LQm1,0:NRMAX),NRMAX,0)
       IF(ID /= 0) return
    END IF

    PhiV (0:NRMAX) =   XL(LQm1,0:NRMAX)
    ErV  (0:NRMAX) = - 2.D0 * R(0:NRMAX) * dfdx(PSI,PhiV,NRMAX,0)
    EthV (0)       =   0.D0
    EthV (1:NRMAX) = - XL(LQm2,1:NRMAX) / R(1:NRMAX)
    EphV (0:NRMAX) = - XL(LQm3,0:NRMAX)
    AphV (0:NRMAX) =   XL(LQm4,0:NRMAX) * rMUb2
    BthV (0:NRMAX) = - 2.D0 * R(0:NRMAX) * dfdx(PSI,AphV ,NRMAX,0)
    RAthV(0:NRMAX) =   XL(LQm5,0:NRMAX) * rMU0
    BphV (0:NRMAX) =   2.D0              * dfdx(PSI,RAthV,NRMAX,0)

    PNeV (0:NRMAX) =   XL(LQe1,0:NRMAX)
    UerV (0)       =   0.D0
    UerV (1:NRMAX) =   XL(LQe2,1:NRMAX )/ PNeV(1:NRMAX) / R(1:NRMAX)
    RUethV(0:NRMAX)=   XL(LQe3,0:NRMAX) / PNeV(0:NRMAX)
    UethV(0)       =   0.D0
    UethV(1:NRMAX) =   XL(LQe3,1:NRMAX) / PNeV(1:NRMAX) / R(1:NRMAX)
    UephV(0:NRMAX) =   XL(LQe4,0:NRMAX) / PNeV(0:NRMAX)
    IF(MDFIXT == 0) THEN
       PeV  (0:NRMAX) =   XL(LQe5,0:NRMAX)
       PTeV (0:NRMAX) =   XL(LQe5,0:NRMAX) / PNeV(0:NRMAX)
    ELSE
       PeV  (0:NRMAX) =   XL(LQe5,0:NRMAX) * PNeV(0:NRMAX)
       PTeV (0:NRMAX) =   XL(LQe5,0:NRMAX)
    END IF
    PNiV (0:NRMAX) =   XL(LQi1,0:NRMAX)
    UirV (0)       =   0.D0
    UirV (1:NRMAX) =   XL(LQi2,1:NRMAX) / PNiV(1:NRMAX) / R(1:NRMAX)
    RUithV(0:NRMAX)=   XL(LQi3,0:NRMAX) / PNiV(0:NRMAX)
    UithV(0)       =   0.D0
    UithV(1:NRMAX) =   XL(LQi3,1:NRMAX) / PNiV(1:NRMAX) / R(1:NRMAX)
    UiphV(0:NRMAX) =   XL(LQi4,0:NRMAX) / PNiV(0:NRMAX)
    IF(MDFIXT == 0) THEN
       PiV  (0:NRMAX) =   XL(LQi5,0:NRMAX)
       PTiV (0:NRMAX) =   XL(LQi5,0:NRMAX) / PNiV(0:NRMAX)
    ELSE
       PiV  (0:NRMAX) =   XL(LQi5,0:NRMAX) * PNiV(0:NRMAX)
       PTiV (0:NRMAX) =   XL(LQi5,0:NRMAX)
    END IF
    PNbV (0:NRMAX) =   XL(LQb1,0:NRMAX)
    IF(ABS(FSRP) > 0.D0) THEN
       DO NR = 0, NRMAX
          BBL = SQRT(BphV(NR)**2 + BthV(NR)**2)
          IF(PNbV(NR) == 0.D0) THEN ! The region without beam particles
             UbthV(NR) = 0.D0
             UbphV(NR) = 0.D0
          ELSE ! The region with beam particles
             IF(NR == 0) THEN ! On axis, poloidal beam velocity is assumed to be zero.
                UbthV(NR) = 0.D0
                UbphV(NR) = XL(LQb4,NR) / PNbV(NR)
             ELSE ! The region except the magnetic axis
                UbthV(NR) = XL(LQb3,NR) / PNbV(NR) / R(NR)
                UbphV(NR) = XL(LQb4,NR) / PNbV(NR)
                if(abs(UbthV(NR)) > Vbabsmax * (BthV(NR) / BBL)) UbthV(NR) = 0.D0
                if(abs(UbphV(NR)) > Vbabsmax * (BphV(NR) / BBL)) UbphV(NR) = 0.D0
             END IF
          END IF
       END DO
    ELSE
       DO NR = 0, NRA-1
          IF(PNbV(NR) == 0.D0) THEN
             UbthV(NR) = 0.D0
             UbphV(NR) = 0.D0
          ELSE
             IF(PNBH == 0.D0 .AND. PNbV(NR) < 1.D-8) THEN
                UbthV(NR) = 0.D0
                UbphV(NR) = 0.D0
             ELSE
                IF(NR == 0) THEN
                   UbthV(NR) = 0.D0
                   UbphV(NR) = XL(LQb4,NR) / PNbV(NR)
                ELSE
                   UbthV(NR) = XL(LQb3,NR) / PNbV(NR) / R(NR)
                   UbphV(NR) = XL(LQb4,NR) / PNbV(NR)
                END IF
             END IF
          END IF
       END DO
       UbthV(NRA:NRMAX) = 0.D0
       UbphV(NRA:NRMAX) = 0.D0
    END IF

!!$    DO NR = 0, NRMAX
!!$       IF (ABS(PNbV(NR)) < 1.D-6) THEN
!!$          UbthV(NR) = 0.D0
!!$          UbphV(NR) = 0.D0
!!$       ELSE
!!$          IF (NR == 0) THEN
!!$             UbthV(NR) = 0.D0
!!$          ELSE
!!$             UbthV(NR) = XL(LQb3,NR) / PNbV(NR) / R(NR)
!!$          END IF
!!$          UbphV(NR) =    XL(LQb4,NR) / PNbV(NR)
!!$       END IF
!!$    END DO
    PN01V(0:NRMAX) =   XL(LQn1,0:NRMAX)
    PN02V(0:NRMAX) =   XL(LQn2,0:NRMAX)
    PN03V(0:NRMAX) =   XL(LQn3,0:NRMAX)

    PNbRPV(0:NRMAX)=   XL(LQr1,0:NRMAX)

    Q(1:NRMAX) = ABS(R(1:NRMAX) * BphV(1:NRMAX) / (RR * BthV(1:NRMAX)))
    Q(0) = FCTR(R(1),R(2),Q(1),Q(2))

!!$    pres0(0:NRMAX) = (XL(LQe5,0:NRMAX) + XL(LQi5,0:NRMAX)) * 1.D20 * rKeV
!!$    ErV0 (0:NRMAX) = - 2.D0 * R(0:NRMAX) * dfdx(PSI,PhiV,NRMAX,0)

    RETURN
  END SUBROUTINE TXCALV

!***************************************************************
!
!   Calculate coefficients
!
!***************************************************************

  SUBROUTINE TXCALC(IC)

    use tx_commons
    use tx_interface, only : EXPV, VALINT_SUB, TRCOFS, dfdx, txmmm95, &
         &                   moving_average, coulog, CORR, INTG_F
    use tx_nclass_mod
    use sauter_mod
    use aux_system
    use tx_ripple

    integer(4) :: IC
    integer(4), save :: NRB = 1
    INTEGER(4) :: NR, NR1, IER, i, MDANOMabs
    REAL(8) :: Sigma0, QL, &
         &     Vte, Vti, Vtb, Wte, Wti, EpsL, rNuPara, &
         &     rNuAsE_inv, rNuAsI_inv, BBL, Va, Wpe2, rGC, PN0tot, &
         &     rH, Smod, PROFML, PROFCL, Dturb, fk, DeL, AJPH, AJTH, EPARA, &
         &     Cs, RhoIT, ExpArg, AiP, DISTAN, UbparaL, &
         &     SiLCL, SiLCthL, SiLCphL, RL, DBW, PTiVA, &
         &     Chicl, factor_bohm, rNuAsIL, ALFA, &
         &     RLOSS, SQZ, rNuDL, ETASL, Ln, LT, etai_chk, kthrhos, &
         &     RhoSOL, V0ave, Viave, DturbA, rLmean, Sitot, &
         &     rGCIM, rGIM, rHIM, OMEGAPR, &  !09/06/17~ miki_m
         &     rLMFL ! 11/06/15 AF
    INTEGER(4) :: NHFM ! miki_m 10-08-06
    real(8), dimension(1:NHFMmx) :: EpsLM 
!    real(8) :: rLmeanL
!!    real(8) :: XXX, SiV, ScxV, Vcr, Wbane, Wbani, ALFA, cap_val, Scxi, Vave, bthl
    real(8), save :: Fcoef = 1.d0
    real(8) :: PTiVav, N02INT, RatSCX, sum1, sum2
    real(8) :: Frdc, Dcoef
    real(8) :: omegaer, omegaere, omegaeri
    real(8) :: EFT, CR
    real(8) :: AITKEN2P
    real(4), dimension(0:NRMAX) :: p_gr2phi
    real(8), dimension(0:NRMAX) :: pres, Ubpara, ETA_coef, JBS_coef
!!    real(8), dimension(0:NRMAX) :: Vexbr
    ! Mainly for derivatives
    real(8), dimension(:), allocatable :: dQdr, dVebdr, dErdr, dBthdr, dTedr, dTidr, &
         &                                dPedr, dPidr, dpdr, dNedr, dNidr, dErdrS, ErVlc, &
         &                                qr, dqr, dBph

    MDANOMabs = abs(MDANOM)

    !     *** Constants ***

    !     Neutral cross section
    !     (NRL Plasma Formulary p52 Eq. (1) (2002))

    Sigma0 = 8.8D-21

    !     Poloidal magnetic field on wall

    IF(FSHL == 0.D0) THEN
       Bthb = rMU0 * rIP * 1.D6 / (2.D0 * PI * RB)
    ELSE
       QL = (Q0 - QA) * (1.D0 - (RB / RA)**2) + QA
       Bthb = BB * RB / (QL * RR)
    END IF

    !     *** Trapped particle fraction ***
    !     (Y. B. Kim, et al., Phys. Fluids B 3 (1990) 2050)
    DO NR = 0, NRMAX
       EpsL = R(NR) / RR        ! Inverse aspect ratio
       ft(NR) = 1.46D0 * SQRT(EpsL) - 0.46 * EpsL * SQRT(EpsL)
    END DO

    ! ************* Auxiliary heating part *************
    call txauxs

    ! ************** For turbulent transport **************
    !   The reason that we keep the pressure throughout iteration is to stabilize numerical
    !      instability caused by the CDBM model, strongly dependent on the pressure gradient.
    !   In order to suppress oscillation of the pressure in the direction of time, 
    !      we take the average between pres and pres0, evaluated at the previous time step,
    !      when differentiating the pressure with respect to psi.
    !   In addition, during iteration, pres is fixed.
    !   This holds true with the radial electric field, ErV, if one prefers.

    pres(0:NRMAX)  = ( PNeV_FIX(0:NRMAX)*PTeV_FIX(0:NRMAX) &
         &            +PNiV_FIX(0:NRMAX)*PTiV_FIX(0:NRMAX)) * 1.D20 * rKeV
    IF(MDANOM > 0 .and. maxval(FSANOM) > 0.D0) pres(0:NRMAX)  = 0.5d0 * (pres(0:NRMAX) + pres0(0:NRMAX))

    allocate(ErVlc(0:NRMAX))
    ErVlc(0:NRMAX) = 0.5d0 * (ErV_FIX(0:NRMAX) + ErV0(0:NRMAX))

!    Vexbr(0)       = 0.d0
!    Vexbr(1:NRMAX) = ErVlc(1:NRMAX) &
!         &         / (R(1:NRMAX) * SQRT(BphV(1:NRMAX)**2 + BthV(1:NRMAX)**2))

    IF(PROFM == 0.D0 .AND. FSDFIX(2) /= 0.D0) THEN
       PROFML = (PTeV(NRA) * rKeV / (16.D0 * AEE * &
            &    SQRT(BphV(NRA)**2 + BthV(NRA)**2))) / FSDFIX(2)
    ELSE
       PROFML = PROFM
    END IF

    IF(PROFC == 0.D0 .AND. FSDFIX(3) /= 0.D0) THEN
       PROFCL = (PTeV(NRA) * rKeV / (16.D0 * AEE * &
            &    SQRT(BphV(NRA)**2 + BthV(NRA)**2))) / FSDFIX(3)
    ELSE
       PROFCL = PROFC
    END IF

    ! ************** Turbulent transport end **************

    ! Banana width
!    Wbane = (Q(0) * SQRT(RR * AME * PTeV(0) * rKeV) / (AEE * BphV(0)))**(2.D0/3.D0)
!    Wbani = (Q(0) * SQRT(RR * AMI * PTiV(0) * rKeV) / (PZ * AEE * BphV(0)))**(2.D0/3.D0)

    ! Banana width at separatrix
!    PTiVA = PTiV(NRA)
    PTiVA = 0.5D0 * PTiV(0)
    DBW = 3.D0 * SQRT(PTiVA * rKEV * AMI) * Q(NRA) / (PZ * AEE * BphV(NRA)) &
         & / SQRT(R(NRA) / RR)

    !     *** Calculate derivatives in advance ***
    !     !!! Caution !!!
    !        The r-derivatives of variables, or near-variables (ex. temperature) should be
    !          estimated by their psi-derivatives multiplied by 2*r because they are
    !          evaluated on the psi-abscissa. On the other hand, those of the other
    !          parameters (ex. radial electric field, poloidal magnetic field) should be
    !          directly calculated.

    allocate(dQdr(0:NRMAX),dVebdr(0:NRMAX),dErdr(0:NRMAX),dBthdr(0:NRMAX),dErdrS(0:NRMAX))
    allocate(dTedr(0:NRMAX),dTidr(0:NRMAX),dPedr(0:NRMAX),dPidr(0:NRMAX),dpdr(0:NRMAX), &
         &   dNedr(0:NRMAX),dNidr(0:NRMAX))
    dQdr  (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,Q    ,NRMAX,0)
!    dVebdr(0:NRMAX) =                     dfdx(R  ,Vexbr,NRMAX,0)
!    dErdr (0:NRMAX) =                     dfdx(R  ,ErV  ,NRMAX,0)
    dErdr (0:NRMAX) =                     dfdx(R  ,ErVlc,NRMAX,0)
    dBthdr(0:NRMAX) =                     dfdx(R  ,BthV ,NRMAX,0)
    dTedr (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PTeV ,NRMAX,0)
    dTidr (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PTiV ,NRMAX,0)
    dPedr (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PeV  ,NRMAX,0)
    dPidr (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PiV  ,NRMAX,0)
    dpdr  (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,pres ,NRMAX,0)
    dNedr (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PNeV ,NRMAX,0)
    dNidr (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PNiV ,NRMAX,0)

!!D02    write(6,'(F8.5,I4,2F11.6)') T_TX,NRB,Rho(NRB),PT02V(NR)

    !  Smoothing Er gradient for numerical stability
    do NR = 0, NRMAX
       dErdrS(NR) = moving_average(NR,dErdr,NRMAX,NRA)
!       if(nr < 9) write(6,'(I4,5F15.7)') NT,Rho(nr),ErV(NR),ErVlc(NR),dErdrS(NR)
    end do
!!$    !  Double smoothing
!!$    do NR = 0, NRMAX
!!$       dErdr(NR) = moving_average(NR,dErdrS,NRMAX,NRA)
!!$    end do
!!$    dErdrS(0:NRMAX) = dErdr(0:NRMAX)

    ! *** Temperatures for neutrals ***

    PT01V(0:NRMAX) =   0.5D0 * AMI * V0**2 / rKeV

    !  --- For thermal neutrals originating from slow neutrals ---

    IF(IC == 1) THEN
       ! SCX : source of PN02V
       DO NR = 0, NRMAX
          SCX(NR) = PNiV(NR) * SiVcx(PTiV(NR)) * PN01V(NR) * 1.D40
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
       call VALINT_SUB(PN02V,NR,N02int,NR)
       sum1 = sum1 + N02int
       sum2 = sum2 + N02int * (0.5d0 * (PTiV(NR-1) + PTiV(NR)))
    end do
    if(sum1 > epsilon(1.d0)) then
       PTiVav = sum2 / sum1
    else
       PTiVav = PTiV(NRA)
    end if
    do nr = 0, nrmax
       PT02V(NR) = PTiVav
    end do

!    PT02V(0:NRMAX) =   PTiV(0:NRMAX)
    PT03V(0:NRMAX) =   PTiV(0:NRMAX)

    ! *********************************

    ! Calculate CDIM coefficient
    RAQPR(0:NRMAX) = 2.D0 * R(0:NRMAX) * &  ! cf. txcalv.f90 L501
         &     dfdx (PSI(0:NRMAX) , RHO(0:NRMAX)**4 / Q(0:NRMAX) , NRMAX , 0)

    ! Orbit squeezing effect for NCLASS
!    IF(IC == 1) THEN
       p_gr2phi(0)  = 0.0
       IF(MDOSQZ == 2) THEN
          DO NR = 1, NRMAX
             p_gr2phi(NR)  = REAL(-RA**2*dErdrS(NR)+RA**2*ErVlc(NR)*dBthdr(NR)/BthV(NR))
          END DO
       ELSE
          p_gr2phi(1:NRMAX) = 0.0
       END IF
!    END IF

    !  Coefficients

    L_NR: DO NR = 0, NRMAX

       Vte = SQRT(2.D0 * ABS(PTeV(NR)) * rKeV / AME)
       Vti = SQRT(2.D0 * ABS(PTiV(NR)) * rKeV / AMI)
       Vtb = SQRT(2.D0 * ABS(PTiV(NR)) * rKeV / AMB)

       PN0tot = (PN01V(NR) + PN02V(NR) + PN03V(NR)) * 1.D20
       SiVizA(NR) = SiViz(PTeV(NR))
       SiVcxA(NR) = SiVcx(PTiV(NR))

       !     *** Ionization rate ***

!old       !     (NRL Plasma Formulary p54 Eq. (12) (2002))
!old       XXX = MAX(PTeV(NR) * 1.D3 / EION, 1.D-2)
!old       SiV = 1.D-11 * SQRT(XXX) * EXP(- 1.D0 / XXX) &
!old            &              / (EION**1.5D0 * (6.D0 + XXX))
!old       rNuION(NR) = FSION * SiV * (PN01V(NR) + PN02V(NR)) * 1.D20
       rNuION(NR) = FSION * SiVizA(NR) * PN0tot

       !     *** Slow neutral diffusion coefficient ***
       !  For example,
       !    E.L. Vold et al., NF 32 (1992) 1433

       !  Maxwellian velocity for slow neutrals
       V0ave = sqrt(4.D0 * V0**2 / Pi)

       !  Total Maxwellian rate coefficients
!!$       if(nr <= NRB) then
!!$          Sitot = (SiVcxA(NRB) * PNiV(NRB) + SiVizA(NRB) * PNeV(NRB)) *1.D20
!!$       else
          Sitot = (SiVcxA(NR) * PNiV(NR) + SiVizA(NR) * PNeV(NR)) *1.D20
!!$       end if

       !  Diffusion coefficient for slow neutrals (short m.f.p.)
!old       D01(NR) = FSD01 * V0**2 &
!old            &   / (Sigma0 * (PN01V(NR) * V0  + (PNiV(NR) + PN02V(NR)) &
!old            &      * Vti) * 1.D20)
!       D01(NR) = FSD01 * V0ave**2 &
!            &   / (3.D0 * SiVcxA(NR) * PNiV(NR) * 1.D20)
       D01(NR) = FSD01 * V0ave**2 &
            &  / (3.D0 * (SiVcxA(NR) * PNiV(NR) + SiVizA(NR) * PNeV(NR)) *1.D20)

       !     *** Thermal neutral diffusion coefficient ***

       !  Maxwellian thermal velocity at the separatrix
!       Viave = sqrt(8.D0 * PTiV(NR) * rKeV / (Pi * AMi))
       Viave = sqrt(8.D0 * PT02V(NR) * rKeV / (Pi * AMi))

       !  Mean free path for fast neutrals
       rLmean = Viave / Sitot
!!D02       !  Locally determined mean free path for fast neutrals
!!D02       rLmeanL = sqrt(8.D0 * PTiV(NR) * rKeV / (Pi * AMi)) &
!!D02            &  / ((SiVcxA(NR) * PNiV(NR) + SiVizA(NR) * PNeV(NR)) *1.D20)

       !  Diffusion coefficient for fast neutrals (short to long m.f.p.)
!old       D02(NR) = FSD02 * Vti**2 &
!old            &   / (Sigma0 * (PNiV(NR) + PN01V(NR) + PN02V(NR)) &
!old            &      * Vti * 1.D20)
!       D02(NR) = FSD02 * Viave**2 &
!            &   / (3.D0 * SiVcxA(NR) * PNiV(NR) * 1.D20)
!!       if(rho(nr) < 0.9d0) then
!!          D02(NR) = D02(NR) * (-0.33333d0*rho(nr)+0.35d0)
!!       else if (rho(nr) >= 0.9d0 .and. rho(nr) <= 1.0d0) then
!!          D02(NR) = D02(NR) * (3.5055d0*rho(nr)**3-2.5055d0)
!!       end if
       D02(NR) = FSD02 * Viave**2 / (3.D0 * Sitot)
       if(nr >= NRB .and. nr <= NRA) D02(NR) = D02(NR) * reduce_D02(NR,rLmean)
!!D02       if(nt >= ntmax-1) write(6,'(I3,F10.6,F11.3,3F10.6,1P2E11.3)') nr,r(nr),d02(nr),rLmean,rLmeanL,PTiV(NR),SCX(NR)

       !     *** Halo neutral diffusion coefficient ***

       Viave = sqrt(8.D0 * PTiV(NR) * rKeV / (Pi * AMi))
       D03(NR) = FSD03 * Viave**2 / (3.D0 * Sitot)

       !     *** Charge exchange rate ***
       !  For thermal ions (assuming that energy of deuterium
       !                    is equivalent to that of proton)

!old       !     (Riviere, NF 11 (1971) 363, Eq.(4))
!old       Scxi = 6.937D-19 * (1.D0 - 0.155D0 * LOG10(PTiV(NR)*1.D3))**2 &
!old            & / (1.D0 + 0.1112D-14 * (PTiV(NR)*1.D3)**3.3d0) ! in m^2
!old       Vave = SQRT(8.D0 * PTiV(NR) * rKeV / (PI * AMI))
!old       rNuiCX(NR) = FSCX * Scxi * Vave * (PN01V(NR) + PN02V(NR)) * 1.D20
       rNuiCX(NR)  = FSCX * SiVcxA(NR) * PN0tot
       !  For thermal loss by charge exchange
       rNuiCXT(NR) = FSCX * SiVcxA(NR) * PN01V(NR) * 1.D20

       !  For beam ions
       !     (C.E. Singer et al., Comp. Phys. Comm. 49 (1988) 275, p.323, B163)
!old       Vave = SQRT(8.D0 * Eb * rKeV / (PI * AMI))
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

       rNuee(NR) = PNeV(NR) * 1.D20 *         AEE**4 &
            &     * coulog(1.d0,PNeV(NR),PTeV(NR),PTiV(NR),AEP,1.d0,AEP,1.d0) &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AME) &
            &     * (ABS(PTeV(NR)) * rKeV)**1.5D0)
       rNuei(NR) = PNiV(NR) * 1.D20 * PZ**2 * AEE**4 &
            &     * coulog(1.d0,PNeV(NR),PTeV(NR),PTiV(NR),AEP,1.d0,PA,PZ) &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AME) &
            &     * (ABS(PTeV(NR)) * rKeV)**1.5D0)
       rNuie(NR) = (AME / AMI) * rNuei(NR)
       rNuii(NR) = PNiV(NR) * 1.D20 * PZ**4 * AEE**4 &
            &     * coulog(1.d0,PNeV(NR),PTeV(NR),PTiV(NR),PA,PZ,PA,PZ) &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AMI) &
            &     * (ABS(PTiV(NR)) * rKeV)**1.5D0)

       !     *** Collision frequency (energy equipartition) ***
       !     (D. V. Sivukhin, Rev. Plasma Phys. Vol. 4, Consultants Bureau (1966))
       !         Energy relaxation time between two kinds of particles
!approximate formula (AME << AMI)      rNuTei(NR) = rNuei(NR) * (2.D0 * AME / AMI)
       rNuTei(NR) = PNiV(NR) * 1.D20 * PZ**2 * AEE**4 &
            &     * coulog(1.d0,PNeV(NR),PTeV(NR),PTiV(NR),AEP,1.d0,PA,PZ) &
            &     / (3.D0 * SQRT(2.D0 * PI) * PI * EPS0**2 * AME * AMI &
            &     * (  ABS(PTeV(NR))*rKeV / AME + ABS(PTiV(NR))*rKeV / AMI)**1.5D0)

       BBL = SQRT(BphV(NR)**2 + BthV(NR)**2)

       rNuPara = CORR(Zeff) * rNuei(NR)
       rNuei1(NR)    =(BthV(NR)**2 * rNuPara + BphV(NR)**2 * rNuei(NR)) / BBL**2
       rNuei2(NR)    = BphV(NR) * BthV(NR) / BBL**2 * (rNuPara - rNuei(NR))
       rNuei3(NR)    =(BphV(NR)**2 * rNuPara + BthV(NR)**2 * rNuei(NR)) / BBL**2
       rNuei2Bth(NR) = BphV(NR) / BBL**2 * (rNuPara - rNuei(NR))
!!$       rNuei1(NR)    = rNuei(NR)
!!$       rNuei2(NR)    = 0.D0
!!$       rNuei3(NR)    = rNuei(NR)
!!$       rNuei2Bth(NR) = 0.D0
!!$       rNuei1(NR)    = rNuPara
!!$       rNuei2(NR)    = 0.D0
!!$       rNuei3(NR)    = rNuPara
!!$       rNuei2Bth(NR) = 0.D0

       Wte = Vte / (Q(NR) * RR) ! Omega_te; transit frequency for electrons
       Wti = Vti / (Q(NR) * RR) ! Omega_ti; transit frequency for ions
       
       EpsL = R(NR) / RR        ! Inverse aspect ratio
       rNuAsE_inv = EpsL**1.5D0 * Wte / (SQRT(2.D0) * rNuei(NR))
       rNuAsI_inv = EpsL**1.5D0 * Wti / (SQRT(2.D0) * rNuii(NR))
       IF(NR /= 0) THEN
          rNuAse(NR) = 1.D0 / rNuAsE_inv
          rNuAsi(NR) = 1.D0 / rNuAsI_inv
       END IF

       IF(ABS(FSHL)==0.D0) THEN
       !     *** Parallel neoclassical viscosity ***
       !     (W. A. Houlberg, et al., Phys. Plasmas 4 (1997) 3230)

          IF(MDOSQZ == 0 .OR. MDOSQZ == 2) THEN
             IF(NR == 0) THEN
                NR1 = 1
             ELSE
                NR1 = NR
             END IF
             CALL TX_NCLASS(NR,rNueNC(NR),rNuiNC(NR),rNue2NC(NR),rNui2NC(NR), &
                  &         FQeth(NR),FQith(NR),ETA_coef(NR),JBS_coef(NR),ETA2(NR),AJBS2(NR), &
                  &         ChiNCpe(NR),ChiNCte(NR),ChiNCpi(NR),ChiNCti(NR), &
                  &         dTedr(NR1),dTidr(NR1),dPedr(NR1),dPidr(NR1),IER, &
                  &         p_gr2phi_in=p_gr2phi(NR))
          ELSE IF(MDOSQZ == 1) THEN
             IF(NR == 0) THEN
                CALL TX_NCLASS(NR,rNueNC(NR),rNuiNC(NR),rNue2NC(NR),rNui2NC(NR), &
                     &         FQeth(NR),FQith(NR),ETA_coef(NR),JBS_coef(NR),ETA2(NR),AJBS2(NR), &
                     &         ChiNCpe(NR),ChiNCte(NR),ChiNCpi(NR),ChiNCti(NR), &
                     &         dTedr(NR+1),dTidr(NR+1),dPedr(NR+1),dPidr(NR+1),IER, &
                     &         dErdrS(NR+1),dBthdr(NR+1),dErdrS(0),dBthdr(0))
             ELSE
                CALL TX_NCLASS(NR,rNueNC(NR),rNuiNC(NR),rNue2NC(NR),rNui2NC(NR), &
                     &         FQeth(NR),FQith(NR),ETA_coef(NR),JBS_coef(NR),ETA2(NR),AJBS2(NR), &
                     &         ChiNCpe(NR),ChiNCte(NR),ChiNCpi(NR),ChiNCti(NR), &
                     &         dTedr(NR),dTidr(NR),dPedr(NR),dPidr(NR),IER, &
                     &         dErdrS(NR),dBthdr(NR))
             END IF
          END IF
!!$          if(nr == 0) then
!!$             write(6,*) 0.5d0*Rho(nr+1),rNueNC(NR),ChiNCpe(NR)
!!$          else
!!$             write(6,*) Rho(nr),rNueNC(NR),ChiNCte(NR)
!!$          end if
          IF(IER /= 0) IERR = IER
!       if(mod(nt,10)==0) write(6,*) rho(nr),((BthV(NR)/BphV(NR))**2-1.d0)*rNuPara+rNueNC(NR)*(BthV(NR)/BphV(NR))**2,UerV(NR)

       ELSE
!!$       IF(RHO(NR) < 1.D0) THEN
!!$          rNueNC(NR) = FSNC * (BphV(0) / BphV(NR))**2 / (2.D0*(1.D0-EpsL**2)**1.5D0) &
!!$               &     * SQRT(PI) * Q(NR)**2 * Wte &
!!$               &     * (1.53D0 / (rNuAsE_inv + SQRT(3.48D0 * rNuAsE_inv) + 1.52D0)) &
!!$               &     / (1.D0 + SQRT(0.37D0 * SQRT(2.D0) * rNuei(NR) / Wte) &
!!$               &                  + 1.25D0 * SQRT(2.D0) * rNuei(NR) / Wte)
!!$          rNuiNC(NR) = FSNC * (BphV(0) / BphV(NR))**2 / (2.D0*(1.D0-EpsL**2)**1.5D0) &
!!$               &     * SQRT(PI) * Q(NR)**2 * Wti &
!!$               &     * (0.53D0 / (rNuAsI_inv + SQRT(0.52D0 * rNuAsI_inv) + 0.56D0)) &
!!$               &     / (1.D0 + SQRT(0.14D0 * SQRT(2.D0) * rNuii(NR) / Wti) &
!!$               &                  + 0.70D0 * SQRT(2.D0) * rNuii(NR) / Wti)
!!$       ELSE
!!$          RL = (R(NR) - RA) / DBW
!!$          rNueNC(NR) = FSNC * (BphV(0) / BphV(NR))**2 / (2.D0*(1.D0-EpsL**2)**1.5D0) &
!!$               &     * SQRT(PI) * Q(NR)**2 * Wte &
!!$               &     * (1.53D0 / (rNuAsE_inv + SQRT(3.48D0 * rNuAsE_inv) + 1.52D0)) &
!!$               &     / (1.D0 + SQRT(0.37D0 * SQRT(2.D0) * rNuei(NR) / Wte) &
!!$               &                  + 1.25D0 * SQRT(2.D0) * rNuei(NR) / Wte) &
!!$               &     * 1.D0 / (1.D0 + RL**2)
!!$          rNuiNC(NR) = FSNC * (BphV(0) / BphV(NR))**2 / (2.D0*(1.D0-EpsL**2)**1.5D0) &
!!$               &     * SQRT(PI) * Q(NR)**2 * Wti &
!!$               &     * (0.53D0 / (rNuAsI_inv + SQRT(0.52D0 * rNuAsI_inv) + 0.56D0)) &
!!$               &     / (1.D0 + SQRT(0.14D0 * SQRT(2.D0) * rNuii(NR) / Wti) &
!!$               &                  + 0.70D0 * SQRT(2.D0) * rNuii(NR) / Wti) &
!!$               &     * 1.D0 / (1.D0 + RL**2)
!!$       END IF

          IF(nr == 0) THEN
             rNueNC(NR) = 0.d0
             rNuiNC(NR) = 0.d0
          ELSE
             rNueNC(NR) = FSNC * SQRT(PI) * Q(NR)**2 &
                  &     * Wte * 1.78D0 / (rNuAsE_inv + 1.78D0)
             rNuiNC(NR) = FSNC * SQRT(PI) * Q(NR)**2 &
                  &     * Wti * 1.78D0 / (rNuAsI_inv + 1.78D0)
          ENDIF
       !     &               * (1 + EpsL**1.5D0 * rNuAsI)
       !     &               / (1 + 1.44D0
       !     &                      * ((EpsL**1.5D0 * rNuAsI)**2
       !     &                         + ( ErVlc(NR)
       !     &                             / ( Vti * BthV(NR)) )**2))
       !!     &                         + ( ErVlc(NR) * BBL
       !!     &                             / ( Vti * BthV(NR)**2) )**2))
          ETA_coef(NR) = 0.d0
          JBS_coef(NR) = 0.d0

       ENDIF

       !     *** Helical neoclassical viscosity ***

!       IF(ABS(FSHL) > 0.D0 .AND. NR > 0) THEN
       IF(ABS(FSHL) > 0.D0 ) THEN
!          IF(int(FSHL) .EQ. 1 ) THEN !! for FSHL = 1 (single Fourier mode)
!!$            Wte = Vte * NCph / RR
!!$            Wti = Vti * NCph / RR
!!$            EpsL = EpsH * (R(NR) / RA)**2
!!$            rNuAsE_inv = EpsL**1.5D0 * Wte / (SQRT(2.D0) * rNuei(NR))
!!$            rNuAsI_inv = EpsL**1.5D0 * Wti / (SQRT(2.D0) * rNuii(NR))
!!$            IF(NR.EQ.0) THEN
!!$               BLinv=0.d0
!!$               omegaer=0.d0
!!$            ELSE
!!$!             QL=(Q0-QA)*(1.D0-(R(NR)/RA)**2)+QA
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
!!$!          UHth  = DBLE(NCth) / DBLE(NCph)
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
             rNueHLM(NHFM, NR) = FSHL * abs(EpsLM(NHFM))**1.5D0 * PTeV(NR) * rKeV          &
                  &            / (AME * RR**2 * rNuei(NR) * (3.D0 + 1.67D0 * omegaere)) ! [s^-1]
!
             omegaeri   = abs(EpsLM(NHFM)) * R(NR) / RR * omegaer**2 / rNuii(NR)**2
             rNuiHLM(NHFM, NR) = FSHL * abs(EpsLM(NHFM))**1.5D0 * PTiV(NR) * rKeV          &
                  &            / (AMP * RR**2 * rNuii(NR) * (3.D0 + 1.67D0 * omegaeri)) ! [s^-1]
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
             
!!$             if (ic == 0 .or. ic == 1) then
!!$                write(*,*) 'IC=',IC, 'NR=', NR, 'R/RA=', R(NR)/RA
!!$                write(*,*) ' NHFM=',NHFM, ',  EpsLM=', EpsLM(NHFM) 
!!$!                write(*,*) '  rNuei=', rNuei(NR)
!!$!                write(*,*) '  omegaere=', omegaere
!!$                write(*,*) '  rNueHLM=', rNueHLM(NHFM, NR) 
!!$             
!!$!                write(*,*) '  rNuii=', rNuii(NR)
!!$!                write(*,*) '  omegaeri=', omegaeri
!!$                write(*,*) '  rNuiHLM=', rNuiHLM(NHFM, NR) 
!!$                write(*,*) ' UHth=',UHth, ',  UHph=', UHph 
!!$             endif

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
          rNueHLthth(0:NRMAX) = 0.D0
          rNueHLthph(0:NRMAX) = 0.D0
          rNueHLphth(0:NRMAX) = 0.D0
          rNueHLphph(0:NRMAX) = 0.D0
          rNuiHLthth(0:NRMAX) = 0.D0
          rNuiHLthph(0:NRMAX) = 0.D0
          rNuiHLphth(0:NRMAX) = 0.D0
          rNuiHLphph(0:NRMAX) = 0.D0
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

       !  Derivatives (beta, safety factor, mock ExB velocity)
       S(NR) = R(NR) / Q(NR) * dQdr(NR)
       Alpha(NR) = - Q(NR)**2 * RR * dpdr(NR) * 2.D0 * rMU0 / (BphV(NR)**2 + BthV(NR)**2)

       !   ***** Thermal diffusivity *****

       !   *** CDBM model ***

       IF (maxval(FSANOM) > 0.D0) THEN
          ! Alfven velocity
          Va = SQRT(BBL**2 / (rMU0 * PNiV(NR) * 1.D20 * AMI))
          ! Squared plasma frequency
          Wpe2 = PNeV(NR) * 1.D20 * AEE**2 / (AME * EPS0)
          ! Magnetic curvature
          rKappa(NR) = FSCBKP * (- R(NR) / RR * (1.D0 - 1.D0 / Q(NR)**2))

          !   *** CDBM model ***
          IF (MDANOMabs == 1) THEN
             ! Arbitrary coefficient for CDBM model
!             rGC = 8.D0
             rGC = 12.D0
             ! Calculate CDBM coefficient
             FCDBM(NR) = TRCOFS(S(NR),Alpha(NR),rKappa(NR))
             ! ExB rotational shear 
             !  e.g. [A.Fukuyama et al PPCF 38 (1996) 1319]
             IF(NR == 0) THEN
                rH = 0.D0
                rG1h2(NR) = 1.D0
             ELSE
                IF(S(NR) /= 0.D0) THEN
                   Smod = sqrt(S(NR)**2 + 0.1d0**2) ! To avoid that S(NR) becomes zero.
                   rH = (Q(NR) * RR / Va) / (Smod * BBL) * dErdrS(NR)
                   ! Turbulence suppression by ExB shear
                   rG1h2(NR) = 1.D0 / (1.D0 + FSCBSH * (rG1 * rH**2))
                ELSE
                   rG1h2(NR) = 0.D0
                END IF
             END IF
             ! Elongation effect (M.Honda and A.Fukuyama, NF 46 (2006) 580)
             fk = (2.D0*SQRT(FSCBEL)/(1.D0+FSCBEL**2))**1.5D0

             ! Turbulent transport coefficient calculated by CDBM model
             Dturb = fk * rGC * FCDBM(NR) * rG1h2(NR) * ABS(Alpha(NR))**1.5D0 &
                  &              * VC**2 / Wpe2 * Va / (Q(NR) * RR)
             !Dturb = MAX(Dturb,1.D-05)

          !   *** CDIM model ***
          ELSE IF (MDANOMabs == 2) THEN
             ! Arbitrary coefficient for CDIM model
             rGCIM = 10.D0
             OMEGAPR = (RA / RR)**2 * (dble(NCph) / dble(NCth)) * RAQPR(NR)
         
             IF(NR == 0) THEN  ! for s=0
                FCDIM(NR) = 0
             ELSE
                FCDIM(NR) = 3.D0 * (0.5D0 * OMEGAPR)**1.5D0 * (RR / RA)**1.5D0 / (Q(NR) * S(NR)**2)
             END IF

             ! ExB rotational shear
             IF(NR == 0) THEN
                rGIM = 0.D0
                rHIM = 0.D0
             ELSE
!                rGIM = rG1
                rGIM = 1.04D0 * R(NR)**2 / ( dpdr(NR) * 2.D0 * rMU0 / (BphV(NR)**2 + BthV(NR)**2) &
                     &                         * OMEGAPR * RA**2 * RR**2 )
!                rHIM = Q(NR) * RR * R(NR)**2 * dVebdr(NR) / (Va * RA)
                rHIM = RA * SQRT( rMU0 * AMI * PNiV(NR) * 1.D20 ) / BthV(NR) * dErdrS(NR) / BBL
             END IF

             ! Turbulence suppression by ExB shear for CDIM mode
             rG1h2IM(NR) = 1.D0 / (1.D0 + FSCBSH * (rGIM * rHIM**2))
             ! Turbulent transport coefficient calculated by CDIM model
             Dturb = rGCIM * FCDIM(NR) * rG1h2IM(NR) * ABS(Alpha(NR))**1.5D0 &
                  &              * VC**2 / Wpe2 * Va / (Q(NR) * RR)

!-----------------memo:beta'=dpdr(NR) * 2.D0 * rMU0 / (BphV(NR)**2 + BthV(NR)**2)-----------

!          IF(NR == 0 .OR. NR == 1) & 
!                         write(6,*) ',NR=',NR,'RAQPR=',RAQPR(NR),'OMEGAPR=',OMEGAPR, &
!               &     'Q=',Q(NR),'S=',S(NR),'FCDIM=',FCDIM(NR),'Dturb=',Dturb

          END IF
          IF(Rho(NR) == 1.D0) DturbA = Dturb 
       ELSE
          rG1h2(NR)  = 0.D0
          FCDBM(NR)  = 0.D0
          rG1h2IM(NR) = 0.D0
          FCDIM(NR)   = 0.D0
          Dturb      = 0.D0
       END IF

       !     *** Turbulent transport of particles ***
       !     ***     Wave-particle interaction    ***

!parail       PROFD1 = 4
       RhoSOL = 1.D0
       IF (RHO(NR) < RhoSOL) THEN
!          DeL = diff_prof(RHO(NR),FSDFIX(1),PROFD,PROFD1) + FSANOM(1) * Dturb
          DeL = diff_prof(RHO(NR),FSDFIX(1),PROFD,PROFD1,PROFD2,0.99d0,0.07d0) &
!          DeL = diff_prof(RHO(NR),FSDFIX(1),PROFD,PROFD1,PROFD2,0.99d0,0.125d0) &
               & + FSANOM(1) * Dturb
       ELSE
          IF(FSPCLD == 0.D0) THEN
             ! Bohm-like diffusivity
             factor_bohm = (FSDFIX(1) * PROFD + FSANOM(1) * Dturb) &
                  &  / (PTeV(NRA) * rKeV / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
             DeL = factor_bohm * PTeV(NR) * rKeV / (16.D0 * AEE * BBL)
          ELSE
             IF(FSANOM(1) == 0.D0) THEN
                ! Fixed value and fixed profile
                DeL = FSPCLD * diff_prof(RhoSOL,FSDFIX(1),PROFD,PROFD1,0.d0)
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

       ! Turbulent pinch term
       VWpch(NR) = VWpch0 * RHO(NR)

       !     *** Turbulent transport of momentum ***

       RhoSOL = 1.D0
!parail       RhoSOL = 0.93D0

       IF (RHO(NR) < RhoSOL) THEN
          DeL = diff_prof(RHO(NR),FSDFIX(2),PROFML,PROFM1,0.d0) + FSANOM(2) * Dturb
!pedestal          if(rho(nr) > 0.9d0) DeL = DeL * exp(-120.d0*(rho(nr)-0.9d0)**2)
       ELSE
          IF(FSPCLM == 0.D0) THEN
             factor_bohm = (FSDFIX(2) * PROFML + FSANOM(2) * Dturb) &
                  &  / (PTeV(NRA) * rKeV / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
!bohm_model2             DeL =  (1.D0 - MOD(FSBOHM,2.D0)) * FSDFIX(2) * PROFML &
!bohm_model2                  &+ FSBOHM * PTeV(NR) * rKeV / (16.D0 * AEE * BBL)
             DeL = factor_bohm * PTeV(NR) * rKeV / (16.D0 * AEE * BBL)
          ELSE
             IF(FSANOM(2) == 0.D0) THEN
                DeL = FSPCLM * diff_prof(RhoSOL,FSDFIX(2),PROFML,PROFM1,0.d0)
             ELSE
                DeL = FSANOM(2) * DturbA
             END IF
!pedestal             DeL = FSPCLM * FSDFIX(2) * PROFML * exp(-120.d0*(rho(nra)-0.9d0)**2)
          END IF
       END IF
       ! Viscosity
       rMue(NR) = rMue0 * DeL
       rMui(NR) = rMui0 * DeL

       !     *** Turbulent transport of heat ***

       IF (RHO(NR) < RhoSOL) THEN
          DeL = diff_prof(RHO(NR),FSDFIX(3),PROFCL,PROFC1,0.d0) + FSANOM(3) * Dturb
!pedestal          if(rho(nr) > 0.9d0) DeL = DeL * exp(-120.d0*(rho(nr)-0.9d0)**2)
       ELSE
          IF(FSPCLC == 0.D0) THEN
             factor_bohm = (FSDFIX(3) * PROFCL + FSANOM(3) * Dturb) &
                  &  / (PTeV(NRA) * rKeV / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
!bohm_model2             DeL =  (1.D0 - MOD(FSBOHM,2.D0)) * FSDFIX(3) * PROFCL &
!bohm_model2                  &+ FSBOHM * PTeV(NR) * rKeV / (16.D0 * AEE * BBL)
             DeL = factor_bohm * PTeV(NR) * rKeV / (16.D0 * AEE * BBL)
          ELSE
             IF(FSANOM(3) == 0.D0) THEN
                DeL = FSPCLC * diff_prof(RhoSOL,FSDFIX(3),PROFCL,PROFC1,0.d0)
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

       ! <omega/m>
       WPM(NR) = WPM0 * PTeV(NR) * rKeV / (RA**2 * AEE * BphV(NR))
       ! Force induced by drift wave (e.q.(8),(13))
       FWe(NR)     = AEE**2         * De(NR) / (PTeV(NR) * rKeV)
       FWi(NR)     = AEE**2 * PZ**2 * De(NR) / (PTeV(NR) * rKeV)
       FWthphe(NR) = FWe(NR)     * BphV(NR)
       FWthphi(NR) = FWi(NR)     * BphV(NR)
       FWthe(NR)   = FWthphe(NR) * BphV(NR)
       FWthi(NR)   = FWthphi(NR) * BphV(NR)
       ! Ad hoc turbulent pinch velocity
       IF(NR == 0) THEN
          FVpch(NR) = 0.D0
       ELSE
          FVpch(NR) = AEE * BphV(NR) * VWpch(NR) / R(NR)
!!$          FVpch(NR) = AEE * BphV(NR) * VWpch(NR)
       END IF

       !     *** Loss to divertor ***

       IF (R(NR) > RA) THEN
          Cs = SQRT((PZ * PTeV(NR) + 3.D0 * PTiV(NR)) * rKeV / AMI)
!          rLMFL = 2.D0 * PI * Q(NR) * RR   ! Length of field line to divertor
          rLMFL = 2.D0 * PI * Q(NR) * RR  &
                  * (1.D0 + LOG(1.D0 + 0.03D0 / (R(NR) - RA)))

          RL = (R(NR) - RA) / DBW! / 2.D0
          rNuL  (NR) = FSLP  * Cs / rLMFL &
               &             * RL**2 / (1.D0 + RL**2)

! Classical heat conduction [s**4/(kg**2.5*m**6)]
! (C S Pitcher and P C Stangeby, PPCF 39 (1997) 779)

          Chicl = (4.D0*PI*EPS0)**2 &
                /(SQRT(AME)*coulog(1.d0,PNeV(NR),PTeV(NR),PTiV(NR), &
                                   AEP,1.d0,PA,PZ)*AEE**4*Zeff)

! When calculating rNuLTe, we fix PNeV and PTeV constant during iteration
!   to obain good convergence.

          rNuLTe(NR) = FSLTE * Chicl * (PTeV_FIX(NR)*rKeV)**2.5D0 &
                             /(rLMFL**2 &
                               * PNeV_FIX(NR)*1.D20) &
                             * RL**2 / (1.D0 + RL**2)
          rNuLTi(NR) = FSLTI * Cs / rLMFL &
                             * RL**2 / (1.D0 + RL**2)
          IF(ABS(FSRP) > 0.D0) THEN
             Ubpara(NR) = (BphV(NR) * UbphV(NR) + BthV(NR) * UbthV(NR)) / BBL
             IF(NR == NRMAX) Ubpara(NR) = AITKEN2P(R(NRMAX), &
                  & Ubpara(NRMAX-1),Ubpara(NRMAX-2),Ubpara(NRMAX-3),&
                  & R(NRMAX-1),R(NRMAX-2),R(NRMAX-3))
             UbparaL = max(Ubpara(NR), FSLP*Cs)
             rNuLB(NR) = UbparaL / rLMFL &
                       * RL**2 / (1.D0 + RL**2)
          END IF
       ELSE
          rNuL(NR) = 0.D0
          rNuLTe(NR) = 0.D0
          rNuLTi(NR) = 0.D0
          rNuLB(NR) = 0.D0
       END IF

       !     *** Current density profile ***

       ! Toloidal current density
       AJPH  = -      AEE * PNeV(NR) * 1.D20 * UephV(NR) &
            &  + PZ * AEE * PNiV(NR) * 1.D20 * UiphV(NR) &
            &  + PZ * AEE * PNbV(NR) * 1.D20 * UbphV(NR)
       ! Poroidal current density
       AJTH  = -      AEE * PNeV(NR) * 1.D20 * UethV(NR) &
            &  + PZ * AEE * PNiV(NR) * 1.D20 * UithV(NR) &
            &  + PZ * AEE * PNbV(NR) * 1.D20 * UbthV(NR)

       ! Parallel current density
       AJPARA(NR)=(BthV(NR)*AJTH     + BphV(NR)*AJPH    )/BBL
       ! Parallel electric field
       EPARA =(BthV(NR)*EthV(NR) + BphV(NR)*EphV(NR))/BBL
       ! Total current density = toroidal current density
       AJ(NR)   = AJPH
       ! Beam current density (parallel)
       AJNB(NR) = (  (PZ * AEE * PNbV(NR) * 1.D20 * UbphV(NR)) * BphV(NR) &
            &      + (PZ * AEE * PNbV(NR) * 1.D20 * UbthV(NR)) * BthV(NR))/BBL! &
!            &    *(1.D0 - PZ / Zeff)
       ! Ohmic heating power
       POH(NR)  = EPARA*(AJPARA(NR)-AJNB(NR)) ! This form neglects BS current!!
!       POH(NR)  = EthV(NR)*AJTH + EphV(NR)*AJPH
!       POH(NR)  = EphV(NR) * AJPH

       !     *** Equipartition power for graphics ***

       PEQe(NR)  = - 1.5D0 * rNuTei(NR) * PNeV(NR) * 1.D20 * (PTeV(NR) - PTiV(NR)) * rKeV
       PEQi(NR)  = - 1.5D0 * rNuTei(NR) * PNeV(NR) * 1.D20 * (PTiV(NR) - PTeV(NR)) * rKeV

       !     *** Ohmic power from Equations for graphics ***

       POHe(NR) = - AEE * EthV(NR) * PNeV(NR) * 1.D20 * UethV(NR) &
            &     - AEE * EphV(NR) * PNeV(NR) * 1.D20 * UephV(NR)
       POHi(NR) =   AEE * EthV(NR) * PNiV(NR) * 1.D20 * UithV(NR) &
            &     + AEE * EphV(NR) * PNiV(NR) * 1.D20 * UiphV(NR)

       !     *** Bremsstraulung loss ***
       !     (NRL Plasma Formulary p57 Eq. (30) (2002))

       PBr(NR) = 5.35D-37 * PZ**2 * PNeV(NR) * PNiV(NR) * 1.D40 * SQRT(PTeV(NR))

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

!!Parail_SOL       if(nr>=nra) write(6,*) nr,RB-RA,sqrt(rmui(nr)*(RR*q(nr)/Vti))

    END DO L_NR

    rNuAse(0) = AITKEN2P(R(0),rNuAse(1),rNuAse(2),rNuAse(3),R(1),R(2),R(3))
    rNuAsi(0) = AITKEN2P(R(0),rNuAsi(1),rNuAsi(2),rNuAsi(3),R(1),R(2),R(3))

    !  Linear extrapolation

!    rNueNC(0) = 2.D0 * rNueNC(0) - rNueNC(1)
!    rNuiNC(0) = 2.D0 * rNuiNC(0) - rNuiNC(1)
!    if(rNueNC(0) < 0.d0) rNueNC(0) = 0.d0
!    if(rNuiNC(0) < 0.d0) rNuiNC(0) = 0.d0
!!    rNueNC(0)  = 0.D0
!!    rNuiNC(0)  = 0.D0
!    rNue2NC(0) = 0.D0
!    rNui2NC(0) = 0.D0
    rNueNC(0) = rNueNC(1)
    rNuiNC(0) = rNuiNC(1)
    FQeth(0) = 0.d0 ! because it's proportional to poloidal rotation.
    FQith(0) = 0.d0 ! because it's proportional to poloidal rotation. 

!!    ChiNCpe(0) = 2.D0 * ChiNCpe(0) - ChiNCpe(1)
!!    ChiNCte(0) = 2.D0 * ChiNCte(0) - ChiNCte(1)
!!    ChiNCpi(0) = 2.D0 * ChiNCpi(0) - ChiNCpi(1)
!!    ChiNCti(0) = 2.D0 * ChiNCti(0) - ChiNCti(1)
    ChiNCpe(0) = ChiNCpe(1)
    ChiNCte(0) = ChiNCte(1)
    ChiNCpi(0) = ChiNCpi(1)
    ChiNCti(0) = ChiNCti(1)
!!$    cap_val = 20.d0
!!$    where(ChiNCpe > cap_val) ChiNCpe = cap_val
!!$    where(ChiNCte > cap_val) ChiNCte = cap_val
!!$    where(ChiNCpi > cap_val) ChiNCpi = cap_val
!!$    where(ChiNCti > cap_val) ChiNCti = cap_val
    ETA2(0)  = 2.D0 * ETA2(0)  - ETA(1)
    AJBS2(0) = 2.D0 * AJBS2(0) - AJBS2(1)
!!$    ! For Neumann condition, finite viscosity is required at the magnetic axis.
!!$    De(0)    = De(1)
!!$    Di(0)    = Di(1)
!!$    rMue(0)  = rMue(1)
!!$    rMui(0)  = rMui(1)

!    write(6,*) INTG_F(SNBe),INTG_F(SNBi)
!    write(6,*) (INTG_F(SNBi*PNBcol_i))*Eb*(1.D20 * rKeV)*2.D0*PI*RR*2.D0*PI/1.D6

    !     *** ExB shearing rate (Hahm & Burrel PoP 1995) ***
    ! omega_ExB = r/q d/dr(q v_E/r)
    !           = r/q [q/r d/dr(v_E) + v_E d/dr(q/r)], where v_E = - ErV / BphV

    allocate(qr(0:NRMAX),dqr(0:NRMAX),dBph(0:NRMAX))
    qr(0)          = 0.d0 ! owing to l'Hopital's rule
    qr(1:NRMAX)    = Q(1:NRMAX) / R(1:NRMAX)
    dqr(0:NRMAX)   = dfdx(R,qr,NRMAX,0)
    dBph(0:NRMAX)  = 2.d0 * R(0:NRMAX) * dfdx(PSI,BphV,NRMAX,0)
    wexb(0:NRMAX)  = abs(- dErdrS(0:NRMAX)/BphV(0:NRMAX) + ErVlc(0:NRMAX)/BphV(0:NRMAX)**2*dBph(0:NRMAX) &
         &               - R(0:NRMAX)*ErVlc(0:NRMAX)/(Q(0:NRMAX)*BphV(0:NRMAX))*dqr(0:NRMAX))

    !     *** Linear stability theory parameter (Zhu, Horton, Sugama PoP 6 (1999) 2503) ***
    ! Ys = sqrt(m_i/Te) abs[R d/dr(Er/(R Bth)) / d/dr(ln q)]
    !      (point: d/dr(Er/(R Bth)) = d/dr (q/r Er/Bph) )
    Ys(0) = 0.d0
    do nr = 1, nrmax
       if(s(nr) /= 0.d0) then
          Ys(nr) = sqrt(AMi/(PTeV(NR)*rKeV)) &
            & *abs(RR*R(NR)*( Q(NR)/R(NR)*(dErdrS(NR)*BphV(NR) &
            &                              -ErVlc(NR)*dBph(NR))/BphV(NR)**2 &
            &                +ErVlc(NR)/BphV(NR)*dqr(NR))/S(NR))
       else
          Ys(nr) = 0.d0
       end if
    end do
    deallocate(qr,dqr,dBph)

    !     *** Linear growth rate for toroidal gamma_etai branch of the ITG mode ***
    !        (F.Crisanti et al, NF 41 (2001) 883)
    gamITG(0:NRMAX,1) = 0.1d0 * sqrt(PTeV(0:NRMAX)*rKeV/AMi)/RA * sqrt(RA/RR) &
         &            * sqrt(RA*abs(dNidr(0:NRMAX))/PNiV(0:NRMAX) &
         &                 + RA*abs(dTidr(0:NRMAX))/PTiV(0:NRMAX)) &
         &            * sqrt(PTiV(0:NRMAX)/PTeV(0:NRMAX))

    gamITG(0,2:3) = 0.d0
    i = 0
    do nr = 1, nrmax
       if(dNidr(NR) /= 0.d0) then
          Ln = PNiV(NR) / abs(dNidr(NR))
       else
          i = i + 1
       end if
       if(dTidr(NR) /= 0.d0) then
          LT = PTiV(NR) / abs(dTidr(NR))
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
                  &       * sqrt(PTiV(NR)*rKeV/AMi) / (Q(NR) * RR)
          end if

          !     *** Linear growth rate for q>2 and s=0 ***
          !        (J.Candy, PoP 11 (2004) 1879)
          kthrhos = sqrt(0.1d0)
          gamITG(NR,3) = kthrhos * sqrt(PTeV(NR)*rKeV/Ami) / RA &
               &       * sqrt(2.d0 * RA / RR * (RA / Ln + RA / LT))
       end if
    end do

    if(MDANOMabs == 3) call txmmm95(dNedr,dNidr,dTedr,dTidr,dQdr)

    !     *** ETB model ***

    IF(MDLETB /= 0) THEN ! ETB on
       Frdc = 0.1d0
       if(Fcoef > Frdc) then
          Dcoef = (1.d0 - Frdc) / DBLE(NTMAX)
          Fcoef = 1.d0 - DBLE(NT) * Dcoef
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

    !     *** Resistivity ***

    DO NR = 0, NRMAX
       ! +++ Original model +++
       EpsL = R(NR) / RR
       BBL = SQRT(BphV(NR)**2 + BthV(NR)**2)
       rNuPara = CORR(Zeff) * rNuei(NR)
       ETASL = AME * rNuPara / (PNeV(NR) * 1.D20 * AEE**2)
       ETA1(NR) = ETASL * (1.D0+(BthV(NR)**2/BBL**2)*rNueNC(NR)/(CORR(Zeff)*rNuei(NR))) &
            &           / (1.D0 + ETA_coef(NR))

       ! +++ Original model +++
!       ALFA = (rNuei1(NR)+rNueNC(NR))/rNuei3(NR)*(BthV(NR)/BphV(NR))**2 &
!            & + 2.D0*rNuei2(NR)/rNuei3(NR)*BthV(NR)/BphV(NR)
!       AJBS1(NR) =- 1.D0 / ((1.D0 + ALFA) * rNuei3(NR) * BphV(NR)) &
!            &    * (  BthV(NR) / BBL * rNueNC(NR) * dpdr(NR) - JBS_coef(NR))
       ALFA = (BBL**2 * rNuPara + BthV(NR)**2 * rNueNC(NR)) / BphV(NR)**2
       AJBS1(NR) =- 1.D0 / (ALFA * BphV(NR)) &
            &    * (  BthV(NR) / BBL * rNueNC(NR) * dpdr(NR) - JBS_coef(NR))
!!$       AJBS1(NR) =- 1.D0 / (ALFA * BphV(NR)) &
!!$            &    * (  BthV(NR) / BBL * rNueNC(NR) * dpdr(NR))

       ! +++ Sauter model +++
       ! Inverse aspect ratio
       IF(NR == 0) THEN
          AJBS3(NR) = 0.D0
       ELSE
          CALL SAUTER(PNeV(NR),PTeV(NR),dTedr(NR),dPedr(NR),PNiV(NR),PTiV(NR),dTidr(NR), &
               &      dPidr(NR),Q(NR),BphV(NR),RR*BthV(NR),RR*BphV(NR),EpsL,RR,PZ,Zeff, &
               &      ft(nr),rlnLei_IN=coulog(1.d0,PNeV(NR),PTeV(NR),PTiV(NR),AEP,1.d0,PA,PZ), &
               &             rlnLii_IN=coulog(1.d0,PNeV(NR),PTeV(NR),PTiV(NR),PA,PZ,PA,PZ), &
               &      JBS=AJBS3(NR),ETA=ETA3(NR))
       END IF

        ! +++ Hirshman, Hawryluk and Birge model +++
       Vte = SQRT(2.D0 * ABS(PTeV(NR)) * rKeV / AME)
       Wte = Vte / (Q(NR) * RR) ! Omega_te; transit frequency for electrons
       rNuAsE_inv = EpsL**1.5D0 * Wte / (SQRT(2.D0) * rNuei(NR))
       EFT  = ft(NR) * rNuAsE_inv / (rNuAsE_inv + (0.58D0 + 0.2D0 * Zeff))
       CR   = 0.56D0 * (3.D0 - Zeff) / ((3.D0 + Zeff) * Zeff)
       ! Spitzer resistivity for hydrogen plasma (parallel direction)
       ETAS(NR) = CORR(1.D0) * AME * rNuei(NR) / (PNeV(NR) * 1.D20 * AEE**2)
       IF(NR == 0) ETA3(NR) = ETAS(NR) * (CORR(Zeff) / CORR(1.D0))
       ETA4(NR) = ETAS(NR) * Zeff * (1.D0 + 0.27D0 * (Zeff - 1.D0)) &
            &   /((1.D0 - EFT) * (1.D0 - CR * EFT) * (1.D0 + 0.47D0 * (Zeff - 1.D0)))
       IF(FSNC /= 0) THEN
          select case(MDLETA)
          case(0)
             ETA(NR) = ETA1(NR)
          case(1)
             ETA(NR) = ETA2(NR)
          case(2)
             ETA(NR) = ETA3(NR)
          case(3)
             ETA(NR) = ETA4(NR)
          case default
             ETA(NR) = ETA2(NR)
          end select
       ELSE
          ! Spitzer resistivity when no neoclassical effects
          ETA(NR) = ETAS(NR)
       END IF

       ! Ohmic current density (parallel)
!       AJOH(NR) = EphV(NR) / ETA(NR)
       EPARA =(BthV(NR)*EthV(NR) + BphV(NR)*EphV(NR))/BBL

       AJOH(NR) = EPARA / ETA(NR)! - AJNB(NR) ! ????????????
    END DO

    !     ***** Ion Orbit Loss *****

    SiLC  (0:NRMAX) = 0.D0
    SiLCth(0:NRMAX) = 0.D0
    SiLCph(0:NRMAX) = 0.D0
    rNuOL (0:NRMAX) = 0.D0
    IF (ABS(FSLC) > 0.D0) THEN
       IF(MDLC == 1) THEN
          ! K. C. Shaing, Phys. Fluids B 4 (1992) 3310
          do nr=1,nrmax
             Vti = SQRT(2.D0 * PTiV(NR) * rKeV / AMI)
             ! Orbit squeezing factor (K.C.Shaing, et al., Phys. Plasmas 1 (1994) 3365)
!?? JCP version             SQZ = 1.D0 - AMI / (PZ * AEE) / BthV(NR)**2 * dVebdr(NR)
             SQZ = 1.D0 - AMI / (PZ * AEE) / BthV(NR)**2 * dErdrS(NR) / Vti

             EpsL = R(NR) / RR
             BBL = sqrt(BphV(NR)**2 + BthV(NR)**2)
             ! rNuDL : deflection collisional frequency at V = Vti
             rNuDL = PNiV(NR) *1.D20 * PZ**2 * PZ**2 * AEE**4 & 
                  &   * coulog(1.d0,PNeV(NR),PTeV(NR),PTiV(NR),PA,PZ,PA,PZ) &
                  &   / (2.D0 * PI * EPS0**2 * AMI**2 * Vti**3) &
                  &   * 0.6289d0 ! <-- Numerical value of (Phi - G) at V = Vti
             rNuAsIL = rNuDL * RR * Q(NR) / (Vti * (abs(SQZ) * EpsL)**1.5D0)
             IF(FSLC == 1.D0) THEN
                rNuOL(NR) = 2.25D0 * rNuDL / (sqrt(PI) * sqrt(2.D0 * abs(SQZ) * EpsL)) &
                     &   * EXP(-(rNuAsIL**0.25D0 + (PZ * AEE * BBL / AMI) &
                     &   * sqrt(abs(SQZ)) / (BphV(NR) * Vti) &
                     &   * ABS(- AphV(NR) + AphV(NRA)) / sqrt(2.D0 * EpsL))**2)
             ELSE
                SiLC(NR) = - 2.25D0 * PNiV(NR) * rNuii(NR) / (sqrt(PI) * sqrt(2.D0 * EpsL)) &
                     &   * EXP(-(rNuAsIL**0.25D0 + AEE * BBL / (BphV(NR) * Vti * AMI) &
                     &         * ABS(- AphV(NR) + AphV(NRA)) / sqrt(2.D0 * EpsL))**2)
                SiLCth(NR) = SiLC(NR) * AMI * UithV(NR) * R(NR)
                SiLCph(NR) = SiLC(NR) * AMI * UiphV(NR)
             END IF
          end do
!          stop
       ELSEIF(MDLC == 2) THEN
          !     S. -I. Itoh and K. Itoh, Nucl. Fusion 29 (1989) 1031
          IF(FSLC == 1.D0) THEN
             ! RLOSS : Numerical coefficient proportional to the relative number of ions
             !         in the loss cone in velocity space
             RLOSS = 0.1D0
             rNuOL(0) = 0.D0
             DO NR = 1, NRMAX
                EpsL = R(NR) / RR
                Vti = SQRT(PTiV(NR) * rKeV / AMI)
                RhoIT = Vti * AMI / (PZ * AEE * BthV(NR))
                RL = (R(NR) - (RA - 1.5D0 * RhoIT)) / DBW ! Alleviation factor
                IF(R(NR) > (RA - RhoIT)) THEN
!                IF(ABS(RA - R(NR)) <= RhoIT .AND. RHO(NR) < 1.D0) THEN
                   ExpArg = -2.D0 * EpsL * (ErVlc(NR) / BthV(NR))**2 / Vti**2
                   ExpArg = ExpArg * (R(NR) / RA)**2
                   rNuOL(NR) = RLOSS * rNuii(NR) / SQRT(EpsL) * EXP(ExpArg) &
                        &    * RL**2 / (1.D0 + RL**2)
                ELSE
                   rNuOL(NR) = 0.D0
                END IF
             END DO
          ELSE
             DO NR = 1, NRA
                EpsL = R(NR) / RR
                Vti = SQRT(PTiV(NR) * rKeV / AMI)
                RhoIT = Vti * AMI / (PZ * AEE * BthV(NR))
                RhoIT = MIN(RhoIT,0.1D0)
                Wti = Vti / (Q(NR) * RR)
                rNuAsI_inv = EpsL**1.5D0 * Wti / (SQRT(2.D0) * rNuii(NR))
                ExpArg = 2.D0 * EpsL / Vti**2 * (ErVlc(NR) / BthV(NR))**2
                AiP = rNuii(NR) * SQRT(EpsL) * rNuAsI_inv / (1.D0 + rNuAsI_inv) &
                     & * EXPV(- ExpArg)
                DO NR1 = NRA, NRMAX
                   DISTAN = (R(NR1) - R(NR)) / RhoIT
                   SiLCL = AiP * EXPV( - DISTAN**2) * PNiV(NR)
                   SiLC(NR) = SiLC(NR) - SiLCL
                   SiLC(NR1) = SiLC(NR1) + SiLCL * R(NR) / R(NR1)
                   SiLCthL = SiLCL * AMI * UithV(NR) * R(NR)
                   SiLCth(NR) = SiLCth(NR) - SiLCthL
                   SiLCth(NR1) = SiLCth(NR1) + SiLCthL * R(NR) / R(NR1)
                   SiLCphL = SiLCL * AMI * UiphV(NR)
                   SiLCph(NR) = SiLCph(NR) - SiLCphL
                   SiLCph(NR1) = SiLCph(NR1) + SiLCphL * R(NR) / R(NR1)
                END DO
             END DO

             SiLC  (0:NRMAX) = FSLC * SiLC  (0:NRMAX)
             SiLCth(0:NRMAX) = FSLC * SiLCth(0:NRMAX)
             SiLCph(0:NRMAX) = FSLC * SiLCph(0:NRMAX)
          END IF
       END IF
    END IF

    !     ***** Toroidal ripple effect *****
    call ripple_effect(dQdr)

    !     ***** Neoclassical toroidal viscosity (NTV) *****
    !      "rNuNTV" and "UastNC" are obtained from NTVcalc
    
!    CALL NTVcalc
    rNuNTV(0:NRMAX) = 0.D0
    UastNC(0:NRMAX) = 0.D0

    deallocate(dQdr,dVebdr,dErdr,dBthdr,dTedr,dTidr,dPedr,dPidr,dpdr,dNedr,dNidr,dErdrS)
    deallocate(ErVlc)

    RETURN
  END SUBROUTINE TXCALC

!***************************************************************
!
!   Given diffusion coefficient profile
!     Input : factor : FSDFIX
!             profd  : shape factor
!             rho    : normalized radius
!             npower : power of rho
!             fgfact : switch for superimpose of Gaussian profile
!          <optional> (valid when fgfact /= 0)
!             mu     : average
!             sigma  : standard deviation
!
!     diff_prof = factor         at rho=0
!                 factor * profd at rho=1
!
!***************************************************************

  real(8) function diff_prof(rho,factor,profd,power,fgfact,mu,sigma)
    use tx_interface, only : fgaussian
    real(8), intent(in) :: rho, factor, profd, power, fgfact
    real(8), intent(in), optional :: mu, sigma
    real(8) :: fmod, fmodmax

    ! Gaussian profile modified by parabolic profile
    if(fgfact /= 0.d0) then
       fmod    = fgaussian(rho,mu,sigma) * (- 4.d0 * (rho - 0.5d0)**2 + 1.d0)
       fmodmax = fgaussian(mu, mu,sigma) * (- 4.d0 * (mu  - 0.5d0)**2 + 1.d0)
       fmod = fgfact * fmod / fmodmax
    else
       fmod = 0.d0
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
