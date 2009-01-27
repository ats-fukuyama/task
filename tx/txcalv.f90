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
    REAL(8), DIMENSION(1:NQM,0:NRMAX), INTENT(IN) :: XL
    integer(4), intent(in), optional :: ID
    INTEGER(4) :: NR
    real(8) :: FCTR

    IF(present(ID)) THEN
       ! The pres0 is the pressure evaluated at the previous time step.
       pres0(0:NRMAX) = (XL(LQe5,0:NRMAX) + XL(LQi5,0:NRMAX)) * 1.D20 * rKeV
       return
    END IF

    PhiV (0:NRMAX) =   XL(LQm1,0:NRMAX)
    ErV  (0:NRMAX) = - 2.D0 * R(0:NRMAX) * dfdx(PSI,PhiV,NRMAX,0)
    EthV (0)       =   0.D0
    EthV (1:NRMAX) = - XL(LQm2,1:NRMAX) / R(1:NRMAX)
    EphV (0:NRMAX) = - XL(LQm3,0:NRMAX)
    AphV (0:NRMAX) =   XL(LQm4,0:NRMAX) * rMUb2
    BthV (0:NRMAX) = - 2.D0 * R(0:NRMAX) * dfdx(PSI,AphV ,NRMAX,0)
    RAthV(0:NRMAX) =   XL(LQm5,0:NRMAX) * rMU0 * AMPm5
    BphV (0:NRMAX) =   2.D0              * dfdx(PSI,RAthV,NRMAX,0)

    PNeV (0:NRMAX) =   XL(LQe1,0:NRMAX)
    UerV (0)       =   0.D0
    UerV (1:NRMAX) =   XL(LQe2,1:NRMAX )/ PNeV(1:NRMAX) / R(1:NRMAX)
    RUethV(0:NRMAX)=   XL(LQe3,0:NRMAX) / PNeV(0:NRMAX)
    UethV(0)       =   0.D0
    UethV(1:NRMAX) =   XL(LQe3,1:NRMAX) / PNeV(1:NRMAX) / R(1:NRMAX)
    UephV(0:NRMAX) =   XL(LQe4,0:NRMAX) / PNeV(0:NRMAX) * AMPe4
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

    PNbRPV(0:NRMAX)=   XL(LQr1,0:NRMAX)

    PT01V(0:NRMAX) =   0.5D0 * AMI * V0**2 / rKeV
    PT02V(0:NRMAX) =   PTiV(0:NRMAX)

    Q(1:NRMAX) = ABS(R(1:NRMAX) * BphV(1:NRMAX) / (RR * BthV(1:NRMAX)))
    Q(0) = FCTR(R(1),R(2),Q(1),Q(2))

    pres0(0:NRMAX) = (XL(LQe5,0:NRMAX) + XL(LQi5,0:NRMAX)) * 1.D20 * rKeV

    RETURN
  END SUBROUTINE TXCALV

!***************************************************************
!
!   Calculate coefficients
!
!***************************************************************

  SUBROUTINE TXCALC

    use tx_commons
    use tx_interface, only : EXPV, VALINT_SUB, TRCOFS, INTG_F, inexpolate, dfdx
    use tx_core_module, only : inv_int
    use tx_nclass_mod
    use sauter_mod

    real(8), parameter :: PAHe = 4.D0, & ! Atomic mass number of He
         &                Enf  = 3.5D3   ! in keV, equal to 3.5 MeV

    INTEGER(4) :: NR, NP, NR1, IER, i, imax, nrl, ist, irip, nr_potato
    REAL(8) :: Sigma0, QL, SL, SLT1, SLT2, PNBP0, PNBT10, PNBT20, SNBPi_INTG, &
         &     PNBPi0, PNBTi10, PNBTi20, SNBTG, SNBPD, PRFe0, PRFi0, &
         &     Vte, Vti, Vtb, XXX, SiV, ScxV, Wte, Wti, EpsL, rNuPara, rNubes, &
         &     rNuAsE_inv, rNuAsI_inv, BBL, Va, Wpe2, rGC, SP, rGBM, &
         &     Ne_m3, Ni_m3, Te_eV, Ti_eV, rat_mass, &
         &     rGIC, rH, PROFCL, PALFL, DCDBM, DeL, AJPH, AJTH, EPARA, Vcr, &
         &     Cs, RhoIT, ExpArg, AiP, DISTAN, UbparaL, &
         &     SiLCL, SiLCthL, SiLCphL, Wbane, Wbani, RL, ALFA, DBW, PTiVA, &
         &     Chicl, rNuBAR, Ecr, factor_bohm, rNuAsIL, &
         &     rhob, rNueff, rNubnc, DCB, DRP, Dltcr, Dlteff, DltR, Vdrift, &
         &     theta1, theta2, thetab, sinthb, dlt, width0, width1, ARC, &
         &     DltRP_rim, theta_rim, diff_min, theta_min, sum_rp, DltRP_ave, &
         &     EbL, logEbL, Scx, Scxb, Vave, Sion, Left, Right, RV0, tmp, &
         &     RLOSS, SQZ, rNuDL, xl, alpha_l, facST, ellE, ellK, Rpotato, ETASL, Tqi0L!, &
    real(8) :: omegaer, omegaere, omegaeri, blinv, bthl
    real(8) :: FCL, EFT, CR, dPTeV, dPTiV, dPPe, dPPi
    real(8) :: DERIV3, AITKEN2P, ELLFC, ELLEC, deriv4
    real(8), dimension(0:NRMAX) :: pres, Vexbr, SNBP, SNBT1, SNBT2, &
         &                         SNBPi, SNBTi1, SNBTi2, &
         &                         SRFe, SRFi, th1, th2, Ubpara
!!rp_conv         &                         ,PNbrpL, DERIV
!!rp_conv    real(8), dimension(1:4,0:NRMAX) :: U
    ! For derivatives
    real(8), dimension(:), allocatable :: dQdr, dVebdr, dErdr, dBthdr, dTedr, dTidr, &
         &                                dPedr, dPidr, dpdr

    !     *** Constants ***

    !     Neutral cross section
    !     (NRL Plasma Formulary p52 Eq. (1) (2002))

    Sigma0 = 8.8D-21

    !     NBI beam speed

    Vb =  SQRT(2.D0 * Eb * rKeV / AMB)
    Vbpara(0:NRMAX) = Vb

    !     Poloidal magnetic field on wall

    IF(FSHL == 0.D0) THEN
       Bthb = rMU0 * rIP * 1.D6 / (2.D0 * PI * RB)
    ELSE
       QL=(Q0-QA)*(1.D0-(RB/RA)**2)+QA
       Bthb = BB*RB/(QL*RR)
    END IF

    !     *** Trapped particle fraction ***
    !     (Y. B. Kim, et al., Phys. Fluids B 3 (1990) 2050)
    DO NR = 0, NRMAX
       EpsL = R(NR) / RR        ! Inverse aspect ratio
       ft(NR) = 1.46D0 * SQRT(EpsL) - 0.46 * EpsL * SQRT(EpsL)
    END DO

    ! **************** Heating part ****************

    !     *** Normalization factor for heating profile ***
    !
    !    SL is a normalization factor for a given heating profile.
    !    It is assumed that the heating profile has a shape of 
    !    exp(-r^2/r_NB^2)*(1-(r/a)^4), in order to renormalize the
    !    heating profile we therefore need to integrate that profile with 
    !    maximum value of unity at axis in advance and calculate the 
    !    normalized factor (i.e. PNBP0, PNBT0) which allows an integration
    !    value of renormalized profile to equal designated value (i.e. PNBH).
    !    The factor (1-(r/a)^4), which is an arbitrary function, is
    !    required to damp the profile of exp(-r^2/r_NB^2) rapidly
    !    because if this factor is not imposed unrealistically large current
    !    driven by NBI is generated near the edge where electron and bulk ion
    !    density are dilute although in this case fast ion density from 
    !    NBI is relatively large.

    !  For NBI heating
    !  *** Perpendicular
    CALL deposition_profile(SNBP,SL,RNBP0,RNBP,'NB')
    PNBP0 = ABS(PNBHP) * 1.D6 / (2.D0 * Pi * RR * SL)

    IF(MDLNBD /= 0) THEN ! Finite Delta effect
       ! For ions
       IF(PNBMPD == 0.D0) THEN ! Exact perpendicular NBI
          CALL deposition_profile(SNBPi,SL,RNBP0,RNBP,'NB_TRAP',0.D0)
       ELSE ! Near perpendicular NBI
          CALL deposition_profile(SNBPi,SL,RNBP0,RNBP,'NB_TRAP',SIGN(1.D0,PNBMPD))
       END IF
       PNBPi0 = ABS(PNBHP) * 1.D6 / (2.D0 * Pi * RR * SL)
    END IF

    !  *** Tangential
    CALL deposition_profile(SNBT1,SLT1,RNBT10,RNBT1,'NB')
    PNBT10 = ABS(PNBHT1) * 1.D6 / (2.D0 * Pi * RR * SLT1)

    IF(MDLNBD > 1) THEN
       ! For ions
       CALL deposition_profile(SNBTi1,SLT1,RNBT10,RNBT1,'NB_PASS',SIGN(1.D0,PNBCD))
       PNBTi10 = ABS(PNBHT1) * 1.D6 / (2.D0 * Pi * RR * SLT1)
    END IF

    CALL deposition_profile(SNBT2,SLT2,RNBT20,RNBT2,'NB')
    PNBT20 = ABS(PNBHT2) * 1.D6 / (2.D0 * Pi * RR * SLT2)

    IF(MDLNBD > 1) THEN
       ! For ions
       CALL deposition_profile(SNBTi2,SLT2,RNBT20,RNBT2,'NB_PASS',SIGN(1.D0,PNBCD))
       PNBTi20 = ABS(PNBHT2) * 1.D6 / (2.D0 * Pi * RR * SLT2)
    END IF

    !  For RF heating (equally heating for electrons and ions)
    CALL deposition_profile(SRFe,SL,RRFe0,RRFew,'RFe')
    PRFe0 = PRFHe * 1.D6 / (2.D0 * Pi * RR * SL)

    CALL deposition_profile(SRFi,SL,RRFi0,RRFiw,'RFi')
    PRFi0 = PRFHi * 1.D6 / (2.D0 * Pi * RR * SL)

    ! Deposition profiles are loaded from the file

    ! *** Defined input *********************************************
    if(iflag_file == 1) then

       ! Birth TOTAL (SNB for heating profiles)
       i = 1
       call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNB)
       ! Calibration by using total power of all ions
       SL = 2.D0 * Pi * INTG_F(SNB)
       PNBHP = infiles(i)%totP * 1.D-6
       SNB(0:NRMAX) = SNB(0:NRMAX) * 1.D-20 &
            &       * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))
       PNBP0 = Eb * rKeV * 1.D20
       ! "or"= infiles(i)%totP / (2.D0 * Pi * RR * (2.D0 * Pi * INTG_F(SNB)))
       ! Birth profile for electrons
       SNBe(0:NRMAX) = SNB(0:NRMAX)
       ! Birth profile for thermal ions
       SNBi(0:NRMAX) = SNB(0:NRMAX)

       ! Birth Trapped (SNBP for graphic of PNBPD)
       i = 2
       call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,2,SNBP)
       ! Calibration by using total power of all ions
       SL = 2.D0 * Pi * INTG_F(SNBP)
       SNBP(0:NRMAX) = SNBP(0:NRMAX) * 1.D-20 &
            &        * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))

       ! *** No orbit effect for all ions ***
       if(MDLNBD == 0) then
          ! Birth Trapped (SNBPDi for trapped beam ions)
          SNBPDi(0:NRMAX) = SNBP(0:NRMAX)
          PNBPi0 = Eb * rKeV * 1.D20
          ! "or" = infiles(2)%totP / (2.D0 * Pi * RR * (2.D0 * Pi * INTG_F(SNBPDi)))

          ! Birth Passing (SNBTGi for passing beam ions)
          i = 3
          call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBTGi)
          ! Calibration by using total power of trapped ions
          SL = 2.D0 * Pi * INTG_F(SNBTGi)
          SNBTGi(0:NRMAX) = SNBTGi(0:NRMAX) * 1.D-20 &
               &          * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))
          PNBT10 = Eb * rKeV * 1.D20
          ! "or" = infiles(i)%totP / (2.D0 * Pi * RR * (2.D0 * Pi * INTG_F(SNBTGi)))

          ! Birth TOTAL (SNBb for beam ions)
          SNBb(0:NRMAX) = SNB(0:NRMAX)

       ! *** Orbit effect ***
       else
          ! Orbit Trapped (SNBPDi for trapped beam ions)
          i = 5
          call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBPDi)
          ! Calibration by using total power of trapped ions
          SL = 2.D0 * Pi * INTG_F(SNBPDi)
          SNBPDi(0:NRMAX) = SNBPDi(0:NRMAX) * 1.D-20 &
               &          * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))
          PNBPi0 = Eb * rKeV * 1.D20
          ! "or" = infiles(i)%totP / (2.D0 * Pi * RR * (2.D0 * Pi * INTG_F(SNBPDi)))

          ! *** Orbit effect for banana ions only ***
          if(MDLNBD == 1) then
             ! Birth Passing (SNBTGi for passing beam ions)
             i = 3
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBTGi)
             ! Calibration by using total power of trapped ions
             SL = 2.D0 * Pi * INTG_F(SNBTGi)
             SNBTGi(0:NRMAX) = SNBTGi(0:NRMAX) * 1.D-20 &
                  &          * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))
             PNBT10 = Eb * rKeV * 1.D20
             ! "or" = infiles(i)%totP / (2.D0 * Pi * RR * (2.D0 * Pi * INTG_F(SNBTGi)))

             ! Birth Passing + Orbit TOTAL (SNBb for beam ions)
             SNBb(0:NRMAX) = SNBPDi(0:NRMAX) + SNBTGi(0:NRMAX)

             ! *** Orbit effect for all ions ***
          else if(MDLNBD == 2) then
             ! Orbit TOTAL (SNBb for beam ions)
             i = 4
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBb)
             ! Calibration by using total power of trapped ions
             SL = 2.D0 * Pi * INTG_F(SNBb)
             SNBb(0:NRMAX) = SNBb(0:NRMAX) * 1.D-20 &
                  &        * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))

             ! Orbit Passing (SNBTGi for passing beam ions)
             i = 6
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBTGi)
             ! Calibration by using total power of trapped ions
             SL = 2.D0 * Pi * INTG_F(SNBTGi)
             SNBTGi(0:NRMAX) = SNBTGi(0:NRMAX) * 1.D-20 &
                  &          * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))
             PNBT10 = Eb * rKeV * 1.D20
             ! "or" = infiles(i)%totP / (2.D0 * Pi * RR * (2.D0 * Pi * INTG_F(SNBTGi)))
          end if
       end if
       ! Torque injection part
       MNB(0:NRMAX)  = PNBCD * SNBTGi(0:NRMAX) * PNBMPD

       ! Local parallel velocity at birth for passing ions
       i = 6
       call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%vb,NRMAX,RHO,5,Vbpara)

    ! *** Arbitrary input *********************************************
    else if(iflag_file == 2) then
       do i = 1, n_infiles
          if(infiles(i)%name == datatype(1)) then ! Perp NB
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBP)
             SL   = 2.D0 * PI * INTG_F(SNBP)
             PNBP0  = ABS(PNBHP) * 1.D6 / (2.D0 * Pi * RR * SL)
             ! Banana orbit effect
             if(MDLNBD /= 0) then! .and. PNBMPD /= 0.D0) then
                do NR = 0, NRMAX
                   ! Passing particles generated by Perp NB
                   SNBTi2(NR) = (1.d0 - ft(NR)) * SNBP(NR) ! tangential part of ions
                   SNBT2(NR)  = SNBTi2(NR)                 ! tangential part of electrons
                   ! Trapped particles generated by Perp NB
                   SNBPi(NR) = ft(NR) * SNBP(NR)           ! perpendicular part of ions
                   SNBP(NR)  = SNBPi(NR)                   ! perpendicular part of electrons
                end do
                PNBTi20 = PNBP0 ! tangential part of ions
                PNBT20  = PNBP0 ! tangential part of electrons
                SNBPi_INTG = INTG_F(SNBPi)
                call shift_prof(SNBPi, 'TRAP',SIGN(1.D0,PNBCD))
                ! calibration of Perp NB amplitude of ions
                PNBPi0 = PNBP0 * (SNBPi_INTG / INTG_F(SNBPi))
!!                call shift_prof(SNBTi2,'PASS',SIGN(1.D0,PNBCD))
             end if
          else if(infiles(i)%name == datatype(2)) then ! Tang NB 1
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,9,SNBT1)
             SLT1 = 2.D0 * PI * INTG_F(SNBT1)
             PNBT10 = ABS(PNBHT1) * 1.D6 / (2.D0 * Pi * RR * SLT1)
          else if(infiles(i)%name == datatype(3)) then ! Tang NB 2
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,9,SNBT2)
             SLT2 = 2.D0 * PI * INTG_F(SNBT2)
             PNBT20 = ABS(PNBHT2) * 1.D6 / (2.D0 * Pi * RR * SLT2)
          else if(infiles(i)%name == datatype(4)) then ! RF
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,0,SRFe)
             SL   = 2.D0 * PI * INTG_F(SRFe)
             PRFe0 = ABS(PRFHe) * 1.D6 / (2.D0 * Pi * RR * SL)
             SL   = 2.D0 * PI * INTG_F(SRFe)
             PRFi0 = ABS(PRFHi) * 1.D6 / (2.D0 * Pi * RR * SL)
          end if
       end do
    end if

    !   NBI total input power (MW)
    PNBH = PNBHP + PNBHT1 + PNBHT2

    !   Ratio of CX deposition rate to IZ deposition rate
    !     (Riviere, NF 11 (1971) 363)

    EbL = Eb * 1.D3 / PA
    logEbL = log10(EbL)

    Scxb = 6.937D-19 * (1.D0 - 0.155D0 * logEbL)**2 / (1.D0 + 1.112D-15 * EbL**3.3D0)

    IF(PNBH == 0.D0) THEN
       RatCX = 0.D0
    ELSE
       IF(Eb > 150.D0) THEN
          Sion = 3.6D-16 / EbL * (- 0.7783D0 + logEbL)
       ELSE
          Sion = 10.D0**(-0.8712D0 * logEbL**2 + 8.156D0 * logEbL - 38.833D0)
       END IF
       RatCX = Scxb / (Scxb + Sion)
    END IF

    !   Alpha heating

    do nr = 0, nrmax
       PALFL = FSNF * bosch_fusion(PTiV(NR),0.5D0*PNiV(NR),0.5D0*PNiV(NR))
       Ecr = 14.8D0 * (PAHe / PA**(2.D0/3.D0)) * PTeV(NR) ! in keV
       PALFi(NR) = PALFL * rate_to_ion(Enf/Ecr)
       PALFe(NR) = PALFL - PALFi(NR)
    end do

    !   Virtual torque input to LQi4

    CALL deposition_profile(Tqi,SL,0.d0,0.d0,'Virtual')
    Tqi0L = Tqi0 / (2.D0 * Pi * RR * SL)
    Tqi(0:NRMAX) = Tqi0L * Tqi(0:NRMAX)

    ! ************** Heating part end **************

    ! ************** For turbulent transport **************
    !   The reason that we keep the pressure throughout iteration is to stabilize numerical
    !      instability caused by the CDBM model, strongly dependent on the pressure gradient.
    !   In order to suppress oscillation of the pressure in the direction of time, 
    !      we take the average between pres and pres0, evaluated at the previous time step,
    !      when differentiating the pressure with respect to psi.

    pres(0:NRMAX)  = ( PNeV_FIX(0:NRMAX)*PTeV_FIX(0:NRMAX) &
         &            +PNiV_FIX(0:NRMAX)*PTiV_FIX(0:NRMAX)) * 1.D20 * rKeV
    IF(ABS(FSCDBM) > 0.D0) pres(0:NRMAX)  = 0.5d0 * (pres(0:NRMAX) + pres0(0:NRMAX))
    Vexbr(0)       = 0.d0
    Vexbr(1:NRMAX) = ErV(1:NRMAX) &
         &         / (R(1:NRMAX) * SQRT(BphV(1:NRMAX)**2 + BthV(1:NRMAX)**2))

    IF(PROFC == 0.D0 .AND. FSDFIX /= 0.D0) THEN
       PROFCL = (PTeV(NRA) * rKeV / (16.D0 * AEE * &
            &    SQRT(BphV(NRA)**2 + BthV(NRA)**2))) / FSDFIX
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

    ! Ripple amplitude
    thetab = 0.5D0 * PI ! pitch angle of a typical banana particle
    sinthb = sin(thetab)
    DO NR = 0, NRMAX
       RL = R(NR) * (1.D0 + (kappa - 1.D0) * sin(thetab))
       DltRP(NR) = ripple(RL,thetab,FSRP) ! Ripple amplitude at banana tip point
                                          ! Mainly use for estimation of diffusive processes
    END DO

    ellK = ELLFC(sin(0.5d0*thetab),IER) ! first kind of complete elliptic function 
    ellE = ELLEC(sin(0.5d0*thetab),IER) ! second kind of complete elliptic function 

!!rp_conv    PNbrpLV(0:NRMAX) = 0.D0

    !     *** Calculate derivatives in advance ***
    !     !!! Caution !!!
    !        The r-derivatives of variables, or near-variables (ex. temperature) should be
    !          estimated by their psi-derivatives multiplied by 2*r because they are
    !          evaluated on the psi-abscissa. On the other hand, those of the other
    !          parameters (ex. radial electric field, poloidal magnetic field) should be
    !          directly calculated.

    allocate(dQdr(0:NRMAX),dVebdr(0:NRMAX),dErdr(0:NRMAX),dBthdr(0:NRMAX))
    allocate(dTedr(0:NRMAX),dTidr(0:NRMAX),dPedr(0:NRMAX),dPidr(0:NRMAX),dpdr(0:NRMAX))
    dQdr  (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,Q    ,NRMAX,0)
    dVebdr(0:NRMAX) =                     dfdx(R  ,Vexbr,NRMAX,0)
    dErdr (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,ErV  ,NRMAX,0)
    dBthdr(0:NRMAX) =                     dfdx(R  ,BthV ,NRMAX,0)
    dTedr (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PTeV ,NRMAX,0)
    dTidr (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PTiV ,NRMAX,0)
    dPedr (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PeV  ,NRMAX,0)
    dPidr (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,PiV  ,NRMAX,0)
    dpdr  (0:NRMAX) = 2.D0 * R(0:NRMAX) * dfdx(PSI,pres ,NRMAX,0)

    !  Coefficients

    L_NR: DO NR = 0, NRMAX

       Vte = SQRT(2.D0 * ABS(PTeV(NR)) * rKeV / AME)
       Vti = SQRT(2.D0 * ABS(PTiV(NR)) * rKeV / AMI)
       Vtb = SQRT(2.D0 * ABS(PTiV(NR)) * rKeV / AMB)

       !     *** Coulomb logarithms ***
       !     (NRL Plasma Formulary p34,35 (2007))
       Ne_m3 = PNeV(NR) * 1.d20 ; Ni_m3 = PNiV(NR) * 1.d20
       Te_eV = PTeV(NR) * 1.d3  ; Ti_eV = PTiV(NR) * 1.d3
       rlnLee(NR) = 30.4d0 - log(sqrt(Ne_m3)/(Te_eV**1.25d0)) &
            &              - sqrt(1.d-5+(log(Te_eV)-2.d0)**2/16.d0)
       rat_mass = AME / AMI
       if(10.d0*PZ**2 > Ti_eV*rat_mass .and. 10.d0*PZ**2 < Te_eV) then ! usual case
          rlnLei(NR) = 30.9d0 - log(sqrt(Ne_m3)/Te_eV)
       else if(Te_eV > Ti_eV*rat_mass .and. Te_eV < 10.d0*PZ**2) then ! rare case
          rlnLei(NR) = 29.9d0 - log(sqrt(Ne_m3)*PZ/Te_eV**1.5d0)
       else if(Te_eV < Ti_eV*PZ*rat_mass) then ! very rare case
          rlnLei(NR) = 36.9d0 - log(sqrt(Ni_m3)/Ti_eV**1.5d0*PZ**2/PA)
       else ! assumption
          rlnLei(NR) = 30.9d0 - log(sqrt(Ne_m3)/Te_eV)
       end if
       rlnLii(NR) = 29.9d0 - log(PZ**2/Ti_eV*sqrt(2.d0*Ni_m3*PZ**2/Ti_eV))
!tokamaks       rlnLei(NR) = 37.8d0 - LOG(SQRT(PNeV(NR)*1.D20)/(PTeV(NR)))
!tokamaks       rlnLii(NR) = 40.3d0 - LOG(PZ**2/PTiV(NR)*SQRT(2.D0*PNiV(NR)*1.D20*PZ**2/PTiV(NR)))

       !     *** Ionization rate ***

!old       !     (NRL Plasma Formulary p54 Eq. (12) (2002))
!old       XXX = MAX(PTeV(NR) * 1.D3 / EION, 1.D-2)
!old       SiV = 1.D-11 * SQRT(XXX) * EXP(- 1.D0 / XXX) &
!old            &              / (EION**1.5D0 * (6.D0 + XXX))
!old       rNuION(NR) = FSION * SiV * (PN01V(NR) + PN02V(NR)) * 1.D20
       rNuION(NR) = FSION * SiViz(PTeV(NR)) * (PN01V(NR) + PN02V(NR)) * 1.D20

       !     *** Slow neutral diffusion coefficient ***

       D01(NR) = FSD0 * V0**2 &
            &   / (Sigma0 * (PN01V(NR) * V0  + (PNiV(NR) + PN02V(NR)) &
            &      * Vti) * 1.D20)

       !     *** Fast neutral diffusion coefficient ***

       D02(NR) = FSD0 * Vti**2 &
            &   / (Sigma0 * (PNiV(NR) + PN01V(NR) + PN02V(NR)) &
            &      * Vti * 1.D20)

       !     *** Charge exchange rate ***
       !  For thermal ions (assuming that energy of deuterium
       !                    is equivalent to that of proton)

!old       !     (Riviere, NF 11 (1971) 363, Eq.(4))
!old       Scx = 6.937D-19 * (1.D0 - 0.155D0 * LOG10(PTiV(NR)*1.D3))**2 &
!old            & / (1.D0 + 0.1112D-14 * (PTiV(NR)*1.D3)**3.3d0) ! in m^2
!old       Vave = SQRT(8.D0 * PTiV(NR) * rKeV / (PI * AMI))
!old       rNuiCX(NR) = FSCX * Scx * Vave * (PN01V(NR) + PN02V(NR)) * 1.D20
       rNuiCX(NR) = FSCX * SiVcx(PTiV(NR)) * (PN01V(NR) + PN02V(NR)) * 1.D20

       !  For beam ions
       !     (C.E. Singer et al., Comp. Phys. Comm. 49 (1988) 275, p.323, B163)
!old       Vave = SQRT(8.D0 * Eb * rKeV / (PI * AMI))
!old       rNubCX(NR) = FSCX * Scxb * Vave * (PN01V(NR) + PN02V(NR)) * 1.D20
       rNubCX(NR) = FSCX * Scxb * Vb * (PN01V(NR) + PN02V(NR)) * 1.D20

       !     *** Collision frequency (with neutral) ***

       rNu0e(NR) = (PN01V(NR) + PN02V(NR)) * 1.D20 * Sigma0 * Vte
       rNu0i(NR) = (PN01V(NR) + PN02V(NR)) * 1.D20 * Sigma0 * Vti
       rNu0b(NR) = (PN01V(NR) + PN02V(NR)) * 1.D20 * Sigma0 * Vtb

       !     *** Collision frequency (momentum transfer, Braginskii's formula) ***

       ! Braginskii's collision time
       rNuee(NR) = PNeV(NR) * 1.D20 *         AEE**4 * rlnLee(NR) &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AME) &
            &     * (ABS(PTeV(NR)) * rKeV)**1.5D0)
       rNuei(NR) = PNiV(NR) * 1.D20 * PZ**2 * AEE**4 * rlnLei(NR) &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AME) &
            &     * (ABS(PTeV(NR)) * rKeV)**1.5D0)
       rNuie(NR) = (AME / AMI) * rNuei(NR)
       ! Caution!: rNuii is not the same as rNui, rNui = sqrt(2) * rNuii
       rNuii(NR) = PNiV(NR) * 1.D20 * PZ**4 * AEE**4 * rlnLii(NR) &
            &     / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * SQRT(AMI) &
            &     * (ABS(PTiV(NR)) * rKeV)**1.5D0)

       !     *** Collision frequency (energy equipartition) ***
       !     (D. V. Sivukhin, Rev. Plasma Phys. Vol. 4, Consultants Bureau (1966))
       !         Energy relaxation time between two kinds of particles
!approximate formula (AME << AMI)      rNuTei(NR) = rNuei(NR) * (2.D0 * AME / AMI)
       rNuTei(NR) = PNiV(NR) * 1.D20 * PZ**2 * AEE**4 * rlnLei(NR) &
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
!       write(6,*) r(nr)/ra,rNueNC(NR)*(BthV(NR)/BBL)**2/(rNuei(NR)+rNueNC(NR)*(BthV(NR)/BBL)**2)*BphV(NR)/BBL

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

       CALL TX_NCLASS(NR,rNueNC(NR),rNuiNC(NR),ETA2(NR),AJBS2(NR), &
            &         ChiNCpe(NR),ChiNCte(NR),ChiNCpi(NR),ChiNCti(NR), &
            &         dErdr(NR),dBthdr(NR),dTedr(NR),dTidr(NR),dPedr(NR),dPidr(NR),IER)
       IF(IER /= 0) IERR = IER

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
       !     &                         + ( ErV(NR)
       !     &                             / ( Vti * BthV(NR)) )**2))
       !!     &                         + ( ErV(NR) * BBL
       !!     &                             / ( Vti * BthV(NR)**2) )**2))

          ENDIF

       !     *** Beam slowing down time (momentum transfer with beam) ***
       ! reference : memo (92/04/02, 92/04/21)
       !             Tokamaks 3rd pp.246 - 252

!!!       Vcr = (3.D0 * SQRT(PI / 2.D0) * PNiV(NR) * PZ**2 / PNeV(NR) * AME / AMI &
!!!            &   * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)**(1.D0/3.D0)
!       Ecr = (9.D0 * PI / 16.D0 * AMI / AME)**(1.D0/3.D0) * AMB / AMI * PTeV(NR) ! in keV
       Ecr = 14.8D0 * (PA / PA**(2.D0/3.D0)) * PTeV(NR) ! in keV
       PNBcol_i(NR) = rate_to_ion(Eb/Ecr)
       PNBcol_e(NR) = 1.d0 - PNBcol_i(NR)
       IF(PNBH == 0.D0 .AND. PNbV(NR) < 1.D-8) THEN
          rNube(NR) = 0.D0
          rNubi(NR) = 0.D0
          rNuB (NR) = 0.D0
       ELSE
          rNube(NR) = PNeV(NR) * 1.D20 * PZ**2 * AEE**4 * rlnLei(NR) &
               &     / (3.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * AMB * AME &
               &             * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)
          rNubi(NR) = PNiV(NR) * 1.D20 * PZ**2 * PZ**2 * AEE**4 * rlnLii(NR) &
               &     / (4.D0 * PI * EPS0**2 * AMB) &
               &     * (1.D0 / AMB + 1.D0 / AMI) &
               &     * 1.D0 / (Vb**3 + 0.75D0 * SQRT(PI) * Vti**3)

          rNuPara = CORR(Zeff) * rNube(NR)
          rNube1(NR)    =(BthV(NR)**2 * rNuPara + BphV(NR)**2 * rNube(NR)) / BBL**2
          rNube2(NR)    = BphV(NR) * BthV(NR) / BBL**2 * (rNuPara - rNube(NR))
          rNube3(NR)    =(BphV(NR)**2 * rNuPara + BthV(NR)**2 * rNube(NR)) / BBL**2
          rNube2Bth(NR) = BphV(NR) / BBL**2 * (rNuPara - rNube(NR))

          !     *** deflection time of beam ions against bulk ions ***
          !     (Takamura (3.26) + Tokamaks 3rd p64)
          !     ** rNuD is similar to rNubi. **
          xl = Vb / Vti
          rNuD(NR) = PNiV(NR) *1.D20 * PZ**2 * PZ**2 * AEE**4 * rlnLii(NR) &
               &   / (2.D0 * PI * EPS0**2 * AMB**2 * Vb**3) &
               &   * (  SQRT(1.D0 - EXP(- 4.D0 * xl**2 / PI)) &
               &      - 2.D0 * xl / (4.D0 * xl**3 + 3.D0 * SQRT(PI)))

          !     *** Beam ion-electron slowing-down time
          !     (Tokamaks 3rd, A_D from below (2.13.2) + (2.14.1), assuming AMB >> AME
          !      and Vb/(SQRT(2)*Vte) << 1)
          !     (originally from the book by Spitzer (1962))
          !     Note that this is a half value of rNube(NR).
          rNubes = PNeV(NR) * 1.D20 * PZ**2 * AEE**4 * rlnLei(NR) &
               & / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * AMB * AME &
               & * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)

          ! The definition given below is a "energy slowing down time" and is not
          ! "particle slowing down time". At present, we assume the former is the
          ! same as the latter.
          rNuB(NR) = rNubes * 3.D0 / LOG(1.D0 + (Eb / Ecr)**1.5D0)
!!          rNuB(NR) = rNube(NR) * 1.5D0 / LOG(1.D0 + (Eb / Ecr)**1.5D0)
!!          rNuB(NR) = rNube(NR) + rNubi(NR)
!       write(6,'(1P4E15.7)') rho(nr),2.d0*pi*q(nr)*RR/Vbpara(nr),1.d0/rnube(nr),1.d0/rnubi(nr)!1.d0/rnub(nr)
       END IF

       !     *** Helical neoclassical viscosity ***

       IF(ABS(FSHL) > 0.D0 .AND. NR > 0) THEN
          Wte = Vte * NCph / RR
          Wti = Vti * NCph / RR
          EpsL = EpsH * (R(NR) / RA)**2
          rNuAsE_inv = EpsL**1.5D0 * Wte / (SQRT(2.D0) * rNuei(NR))
          rNuAsI_inv = EpsL**1.5D0 * Wti / (SQRT(2.D0) * rNuii(NR))
          IF(NR.EQ.0) THEN
             omegaer=0.d0
             BLinv=0.d0
          ELSE
             QL=(Q0-QA)*(1.D0-(R(NR)/RA)**2)+QA
             Bthl = BB*R(NR)/(QL*RR)
             BLinv=BB/Bthl
             omegaer=ErV(NR)/(BB*R(NR))
          ENDIF
          omegaere=EpsL*R(NR) / RR * omegaer**2 / rNuei(NR)**2
          omegaeri=EpsL*R(NR) / RR * omegaer**2 / rNuii(NR)**2
          rNueHL(NR) = FSHL * Wte * BLinv * rNuAsE_inv &
          &            /(3.D0+1.67*omegaere)
          rNuiHL(NR) = FSHL * Wti * BLinv * rNuAsI_inv &
          &            /(3.D0+1.67*omegaeri)

          UHth=(RR/Ncph)/SQRT((RR/NCph)**2+(R(NR)/NCth)**2)
          UHph=(R(NR)/Ncth)/SQRT((RR/NCph)**2+(R(NR)/NCth)**2)
          UHth  = DBLE(NCth) / DBLE(NCph)
          UHph  = 1.D0

          rNueHLthth(NR)=UHth*UHth*rNueHL(NR)
          rNueHLthph(NR)=UHth*UHph*rNueHL(NR)
          rNueHLphth(NR)=UHth*UHph*rNueHL(NR)
          rNueHLphph(NR)=UHph*UHph*rNueHL(NR)
          rNuiHLthth(NR)=UHth*UHth*rNuiHL(NR)
          rNuiHLthph(NR)=UHth*UHph*rNuiHL(NR)
          rNuiHLphth(NR)=UHth*UHph*rNuiHL(NR)
          rNuiHLphph(NR)=UHph*UHph*rNuiHL(NR)

!          rNueHL(NR) = FSHL * SQRT(PI) &
!               &     * Wte * 1.78D0 / (rNuAsE_inv + 1.78D0)
!          rNuiHL(NR) = FSHL * SQRT(PI) &
!               &     * Wti * 1.78D0 / (rNuAsI_inv + 1.78D0)
       ELSE
          rNueHL(0:NRMAX) = 0.D0
          rNuiHL(0:NRMAX) = 0.D0
       END IF

       !  Derivatives (beta, safety factor, mock ExB velocity)
       S(NR) = R(NR) / Q(NR) * dQdr(NR)
       Alpha(NR) = - Q(NR)**2 * RR * dpdr(NR) * 2.D0 * rMU0 / (BphV(NR)**2 + BthV(NR)**2)

       !     *** Wave-particle interaction ***

       IF (ABS(FSCDBM) > 0.D0) THEN
          ! Alfven velocity
          Va = SQRT(BBL**2 / (rMU0 * PNiV(NR) * 1.D20 * AMI))
          ! Squared plasma frequency
          Wpe2 = PNeV(NR) * 1.D20 * AEE**2 / (AME * EPS0)
          ! Arbitrary coefficient for CDBM model
          rGC = 8.D0
          ! Magnetic curvature
          rKappa(NR) = - R(NR) / RR * (1.D0 - 1.D0 / Q(NR)**2)
          ! Calculate CDBM coefficient
          FCDBM(NR) = TRCOFS(S(NR),Alpha(NR),rKappa(NR))
          ! ExB rotational shear
          IF(NR == 0) THEN
             rH=0.D0
          ELSE
             rH = Q(NR) * RR * R(NR) * dVebdr(NR) / (Va * S(NR))
          END IF
          
          ! Turbulence suppression by ExB shear
          rG1h2(NR) = 1.D0 / (1.D0 + rG1 * rH**2)
          ! Turbulent transport coefficient calculated by CDBM model
          DCDBM = rGC * FCDBM(NR) * rG1h2(NR) * ABS(Alpha(NR))**1.5D0 &
               &              * VC**2 / Wpe2 * Va / (Q(NR) * RR)
          !DCDBM = MAX(DCDBM,1.D-05)
       ELSE
          rG1h2(NR)  = 0.D0
          FCDBM(NR)  = 0.D0
          DCDBM      = 0.D0
       END IF

       !     *** Turbulent transport of particles ***

!       PROFD = 8.D0
!       PROFD = 3.D0
!       PROFD = 2.D0
       IF (RHO(NR) < 1.D0) THEN
!parail          DeL = FSDFIX * (1.D0 + (PROFD - 1.D0) * (R(NR) / RA)**4) + FSCDBM * DCDBM
          DeL = FSDFIX * (1.D0 + (PROFD - 1.D0) * (R(NR) / RA)**3) + FSCDBM * DCDBM
!          DeL = FSDFIX * (1.D0 + (PROFD - 1.D0) * (R(NR) / RA)**2) + FSCDBM * DCDBM
       ELSE
          IF(FSPCLD == 0.D0) THEN
             factor_bohm = (FSDFIX * PROFD + FSCDBM * DCDBM) &
                  &  / (PTeV(NRA) * rKeV / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
             DeL = factor_bohm * PTeV(NR) * rKeV / (16.D0 * AEE * BBL)
          ELSE
             DeL = FSPCLD * FSDFIX * PROFD
          END IF
!!$          DeL = FSDFIX * PROFD + FSCDBM * DCDBM
       END IF
       ! Particle diffusivity
       De(NR)   = De0   * DeL
       Di(NR)   = Di0   * DeL

       !     *** Turbulent transport of momentum and heat ***

       IF (RHO(NR) < 1.D0) THEN
!parail       IF (RHO(NR) < 0.93D0) THEN
          DeL = FSDFIX * (1.D0 + (PROFCL - 1.D0) * (R(NR) / RA)**2) + FSCDBM * DCDBM
!pedestal          if(rho(nr) > 0.9d0) DeL = DeL * exp(-120.d0*(rho(nr)-0.9d0)**2)
       ELSE
          IF(FSPCLC == 0.D0) THEN
             factor_bohm = (FSDFIX * PROFCL + FSCDBM * DCDBM) &
                  &  / (PTeV(NRA) * rKeV / (16.D0 * AEE * SQRT(BphV(NRA)**2 + BthV(NRA)**2)))
!bohm_model2             DeL =  (1.D0 - MOD(FSBOHM,2.D0)) * FSDFIX * PROFCL &
!bohm_model2                  &+ FSBOHM * PTeV(NR) * rKeV / (16.D0 * AEE * BBL)
             DeL = factor_bohm * PTeV(NR) * rKeV / (16.D0 * AEE * BBL)
          ELSE
             DeL = FSPCLC * FSDFIX * PROFCL
!pedestal             DeL = FSPCLC * FSDFIX * PROFCL * exp(-120.d0*(rho(nra)-0.9d0)**2)
          END IF
       END IF
!       DeL = 3.d0
       ! Viscosity
       rMue(NR) = rMue0 * DeL
       rMui(NR) = rMui0 * DeL
       ! Thermal diffusivity
       Chie(NR) = Chie0 * DeL
       Chii(NR) = Chii0 * DeL

       ! <omega/m>
       WPM(NR) = WPM0 * PTeV(NR) * rKeV / (RA**2 * AEE * BphV(NR))
       ! Force induced by drift wave (e.q.(8),(13))
       FWthe(NR)   = AEE**2         * BphV(NR)**2 * De(NR) / (PTeV(NR) * rKeV)
       FWthi(NR)   = AEE**2 * PZ**2 * BphV(NR)**2 * Di(NR) / (PTiV(NR) * rKeV)
       FWthphe(NR) = AEE**2         * BphV(NR)    * De(NR) / (PTeV(NR) * rKeV)
       FWthphi(NR) = AEE**2 * PZ**2 * BphV(NR)    * Di(NR) / (PTiV(NR) * rKeV)

       ! Work induced by drift wave
       IF(MDLWTB == 1) THEN
          WNthe (NR) = WPE0 *         AEE    * BphV(NR)    * De(NR)
          WEMthe(NR) = WPE0 *         AEE**2 * BphV(NR)    * De(NR) &
               &                       / (PTeV(NR) * rKeV)
          WWthe (NR) = WPE0 *         AEE**2 * BphV(NR)**2 * De(NR) &
               &            * WPM (NR) / (PTeV(NR) * rKeV)
          WT1the(NR) = WPE0 *         AEE    * BphV(NR)    *(rMue(NR) - 0.5D0 * De(NR)) &
               &                       / (PTeV(NR) * rKeV)
          WT2the(NR) = WPE0 *         AEE    * BphV(NR)    *(rMue(NR) - 0.5D0 * De(NR))
          WNthi (NR) = WPI0 * PZ    * AEE    * BphV(NR)    * Di(NR)
          WEMthi(NR) = WPI0 * PZ**2 * AEE**2 * BphV(NR)    * Di(NR) &
               &                       / (PTiV(NR) * rKeV)
          WWthi (NR) = WPI0 * PZ**2 * AEE**2 * BphV(NR)**2 * Di(NR) &
               &            * WPM (NR) / (PTiV(NR) * rKeV)
          WT1thi(NR) = WPI0 * PZ    * AEE    * BphV(NR)    *(rMui(NR) - 0.5D0 * Di(NR)) &
               &            / (PTiV(NR) * rKeV)
          WT2thi(NR) = WPI0 * PZ    * AEE    * BphV(NR)    *(rMui(NR) - 0.5D0 * Di(NR))
       END IF

       !     *** Heating profile ***

       ! For graphic
       IF(iflag_file == 1) THEN
          PNBPD(NR) = PNBP0 * SNBP(NR) ! Power of trapped ions
          PNBTG(NR) = PNBT10 * SNBTGi(NR) ! Power of passing ions
          PNB(NR)   = PNBPD(NR) + PNBTG(NR)
       ELSE
          PNBPD(NR) = PNBP0 * SNBP(NR)
          PNBTG(NR) = PNBT10 * SNBT1(NR) + PNBT20 * SNBT2(NR)
          PNB(NR)   = PNBPD(NR) + PNBTG(NR)
       END IF
       ! For graphic and calculation
       IF(MDLNBD == 0) THEN
          SNBTG     = PNBTG(NR) / (Eb * rKeV * 1.D20)
          SNBPD     = PNBP0 * SNBP(NR) / (Eb * rKeV * 1.D20)
          SNBTGi(NR)= SNBTG
          SNBPDi(NR)= SNBPD
          SNB(NR)   = SNBTG + SNBPD
          SNBe(NR)  = SNB(NR)
          SNBi(NR)  = SNB(NR)
          SNBb(NR)  = SNB(NR)
          MNB(NR)   = PNBCD * SNBTG
       ELSE
          IF(iflag_file /= 1) THEN
             ! Source profile for passing ions in temporal and graphic use
             IF(MDLNBD == 1) THEN
                SNBTGi(NR)=(PNBT10 * SNBT1(NR) + PNBT20 * SNBT2(NR)) / (Eb * rKeV * 1.D20)
             ELSE
                SNBTGi(NR)=(PNBTi10 * SNBTi1(NR) + PNBTi20 * SNBTi2(NR)) / (Eb * rKeV * 1.D20)
             END IF
             ! Source profile for trapped ions with banana orbit effect
             !   in temporal and graphic use
             SNBPDi(NR)= PNBPi0 * SNBPi(NR) / (Eb * rKeV * 1.D20)
             ! Birth profiles for heating power
             SNB(NR)   = PNB(NR) / (Eb * rKeV * 1.D20)
             ! Birth profiles for electrons and thermal ions
             SNBe(NR)  = SNB(NR)
             SNBi(NR)  = SNB(NR)
!!old fashion             SNBi(NR)  = SNBPDi(NR) + SNBTGi(NR)
             ! Source profiles for beam ions with banana orbit effect
             SNBb(NR)  = SNBPDi(NR) + SNBTGi(NR)
             ! Torque injection part
             MNB(NR)   = PNBCD * SNBTGi(NR) * PNBMPD
          END IF
          ! Note: in case of iflag_file == 1, these terms have been already defined above.
       END IF
       PRFe(NR)  = PRFe0 * SRFe(NR)
       PRFi(NR)  = PRFi0 * SRFi(NR)

       !     *** Loss to divertor ***

!       IF (R(NR) + DBW > RA) THEN
       IF (R(NR) > RA) THEN
!          Cs = SQRT(2.D0 * PTeV(NR) * rKeV / AMI)
          Cs = SQRT((PZ * PTeV(NR) + 3.D0 * PTiV(NR)) * rKeV / AMI)
          RL = (R(NR) - RA) / DBW! / 2.D0
          rNuL  (NR) = FSLP  * Cs / (2.D0 * PI * Q(NR) * RR) &
               &             * RL**2 / (1.D0 + RL**2)
          Chicl = (4.D0*PI*EPS0)**2/(SQRT(AME)*rlnLei(NR)*AEE**4*Zeff)*AEE**2.5D0

          ! When calculating rNuLTe, we fix PNeV and PTeV constant during iteration
          !   to obain good convergence.
          rNuLTe(NR) = FSLTE * Chicl * (PTeV_FIX(NR)*1.D3)**2.5D0 &
               &                  /((2.D0 * PI * Q(NR) * RR)**2 * PNeV_FIX(NR)*1.D20) &
               &             * RL**2 / (1.D0 + RL**2)
          rNuLTi(NR) = FSLTI * Cs / (2.D0 * PI * Q(NR) * RR) &
               &             * RL**2 / (1.D0 + RL**2)
          IF(ABS(FSRP) > 0.D0) THEN
             Ubpara(NR) = (BphV(NR) * UbphV(NR) + BthV(NR) * UbthV(NR)) / BBL
             IF(NR == NRMAX) Ubpara(NR) = AITKEN2P(R(NRMAX), &
                  & Ubpara(NRMAX-1),Ubpara(NRMAX-2),Ubpara(NRMAX-3),&
                  & R(NRMAX-1),R(NRMAX-2),R(NRMAX-3))
             UbparaL = max(Ubpara(NR), FSLP*Cs)
             rNuLB(NR) = UbparaL / (2.D0 * PI * Q(NR) * RR) * RL**2 / (1.D0 + RL**2)
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

       !     *** NBI power deposition for graphics ***

       PNBe(NR) = Eb * SNB(NR) * PNBcol_e(NR) * (1.D20 * rKeV)
       PNBi(NR) = Eb * SNB(NR) * PNBcol_i(NR) * (1.D20 * rKeV) &
            &   + AMb * Vb * MNB(NR) * (BthV(NR)*UithV(NR)+BphV(NR)*UiphV(NR))/BBL * 1.D20

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
    rNueNC(0) = AITKEN2P(R(0),rNueNC(1),rNueNC(2),rNueNC(3),R(1),R(2),R(3))
    rNuiNC(0) = AITKEN2P(R(0),rNuiNC(1),rNuiNC(2),rNuiNC(3),R(1),R(2),R(3))
    if(rNueNC(0) < 0.d0) rNueNC(0) = 0.d0
    if(rNuiNC(0) < 0.d0) rNuiNC(0) = 0.d0
    ChiNCpe(0) = AITKEN2P(R(0),ChiNCpe(1),ChiNCpe(2),ChiNCpe(3),R(1),R(2),R(3))
    ChiNCte(0) = AITKEN2P(R(0),ChiNCte(1),ChiNCte(2),ChiNCte(3),R(1),R(2),R(3))
    ChiNCpi(0) = AITKEN2P(R(0),ChiNCpi(1),ChiNCpi(2),ChiNCpi(3),R(1),R(2),R(3))
    ChiNCti(0) = AITKEN2P(R(0),ChiNCti(1),ChiNCti(2),ChiNCti(3),R(1),R(2),R(3))

!    write(6,*) INTG_F(SNBe),INTG_F(SNBi)
!    write(6,*) (INTG_F(SNBi*PNBcol_i))*Eb*(1.D20 * rKeV)*2.D0*PI*RR*2.D0*PI/1.D6

    !     *** Resistivity ***

    DO NR = 0, NRMAX
       ! +++ Original model +++
       EpsL = R(NR) / RR
       ETASL = CORR(Zeff) * AME * rNuei(NR) / (PNeV(NR) * 1.D20 * AEE**2)
       ETA1(NR) = ETASL * (1.D0+(BthV(NR)**2/(BthV(NR)**2+BphV(NR)**2))*rNueNC(NR)/(CORR(Zeff)*rNuei(NR)))
!!$       ETA1(NR) = ETASL * (1.D0+(BthV(NR)**2/(BthV(NR)**2+BphV(NR)**2))*rNueNC(NR)/(2.D0*CORR(Zeff)*rNuei(NR)))
!!$       ETA1(NR) = AME / (AEE**2 * PNeV(NR)*1.D20) &
!!$            &   * (CORR(Zeff)*rNuei(NR) + BthV(NR)**2/(BphV(NR)**2+BthV(NR)**2)*rNueNC(NR))

       ! +++ Sauter model +++
       ! Inverse aspect ratio
       dPTeV = dTedr(NR) * RA
       dPTiV = dTidr(NR) * RA
       dPPe  = dPedr(NR) * RA
       dPPi  = dPidr(NR) * RA
       CALL SAUTER(PNeV(NR),PTeV(NR),dPTeV,dPPe,PNiV(NR),PTiV(NR),dPTiV,dPPi, &
            &      Q(NR),BphV(NR),RR*RA*BthV(NR),RR*BphV(NR),EpsL,RR,PZ,Zeff,ft(nr), &
            &      rlnLei_IN=rlnLei(NR),rlnLii_IN=rlnLii(NR),&
            &      JBS=AJBS3(NR),ETA=ETA3(NR))
       IF(NR == 0) AJBS3(NR) = 0.D0

        ! +++ Hirshman, Hawryluk and Birge model +++
       Vte = SQRT(2.D0 * ABS(PTeV(NR)) * rKeV / AME)
       Wte = Vte / (Q(NR) * RR) ! Omega_te; transit frequency for electrons
       rNuAsE_inv = EpsL**1.5D0 * Wte / (SQRT(2.D0) * rNuei(NR))
       EFT  = ft(NR) * rNuAsE_inv / (rNuAsE_inv + (0.58D0 + 0.2D0 * Zeff))
       CR   = 0.56D0 * (3.D0 - Zeff) / ((3.D0 + Zeff) * Zeff)
       ! Spitzer resistivity for hydrogen plasma (parallel direction)
       ETAS(NR) = CORR(1.D0) * AME * rNuei(NR) / (PNeV(NR) * 1.D20 * AEE**2)
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

!!$       EpsL  = R(NR) / RR
!!$       BBL = SQRT(BphV(NR)**2 + BthV(NR)**2)
!!$       ALFA = (rNuei1(NR)+rNueNC(NR))/rNuei3(NR)*(BthV(NR)/BphV(NR))**2 &
!!$            & + 2.D0*rNuei2(NR)/rNuei3(NR)*BthV(NR)/BphV(NR)
!!$       AJBS1(NR) = -1.D0 / (1.D0 + ALFA) * BthV(NR) / (BBL * BphV(NR)) * rNueNC(NR) / rNuei3(NR) * (dPPe(NR) + dPPi(NR)) * 1.D20 * rKeV

!       write(6,*) r(nr)/ra,epara/eta(nr),aj(nr)-ajbs1(nr)-ajnb(nr)
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
             SQZ = 1.D0 - AMI / (PZ * AEE) / BthV(NR)**2 * dErdr(NR) / Vti

             EpsL = R(NR) / RR
             BBL = sqrt(BphV(NR)**2 + BthV(NR)**2)
             ! rNuDL : deflection collisional frequency at V = Vti
             rNuDL = PNiV(NR) *1.D20 * PZ**2 * PZ**2 * AEE**4 * rlnLii(NR) &
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
                   ExpArg = -2.D0 * EpsL * (ErV(NR) / BthV(NR))**2 / Vti**2
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
                ExpArg = 2.D0 * EpsL / Vti**2 * (ErV(NR) / BthV(NR))**2
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

    !     ***** Ripple loss transport *****

    rip_rat(0:NRMAX) = 0.D0
    IF(ABS(FSRP) > 0.D0) THEN
       ! +++ Convective loss +++

       ! *** Physical aspect of trapped particles in local ripple wells ***
       ! As pointed out in Goldston and Towner (1981) in section 2.5, banana particles
       ! which have small perpendicular component of their velocity, i.e. those
       ! residing near the banana tip point, are apt to be trapped in ripple wells,
       ! provided that they lie in ripple well region with alpha < 1.
       ! In this code, however, we cannot calculate banana tip point for each particle.
       ! Then we assume that all the ripple trapped particles are generated from the
       ! banana particles at the rim of the ripple well region, i.e. alpha = 1.
       ! We readily find two rims in a upper-half plane of a flux surface: LFS and HFS.
       ! Since the contribution from the HFS is usually negligible of their smallness,
       ! we only take account of the rim at the LFS.
       ! The banana particles have an opportunity to be trapped in local ripple wells
       ! when their banana tips lie in ripple well region. As will be discussed in 
       ! "Diffusive loss", we assume the poloidal angle of banana tip is 90 degrees on
       ! behalf of all the banana particles. However, on most of magnetic surfaces,
       ! the banana tips of 90 degrees do not lie in ripple well region, hence we cannot
       ! assume it for ripple well trapping. In this case we assume that the banana
       ! particles have an maxwellian distribution over the poloidal plane whose peak
       ! is taken at theta=90 degree, which is taken to be consistent with the assumption
       ! for diffusive loss case. Therefore, when we consider the scattering at the rim
       ! of ripple well region, we estimate the banana particle density which are going
       ! to be trapped as {sqrt(delta(theta_rim)) * fmaxwell(theta_rim)} * nb.

       ! Ripple well region
       DO NR = 1, NRMAX
          RL = R(NR)
          EpsL = RL / RR
          ! For LFS
          theta1  = 0.d0
          i = 0
          imax = 101
          dlt = 1.d0 / (imax - 1)
          irip = 0
          do
             i = i + 1
             irip = irip + 1
             if(i == imax) then
!                write(6,'(A,I3)') "LFS rim of ripple well region not detected at NR = ",NR
                theta1 = PI
                exit
             end if
             theta1 = theta1 + PI * dlt
             RL = R(NR) * (1.D0 + (kappa - 1.D0) * sin(theta1))
             width0 = ripple(RL,theta1,FSRP)
             width1 = EpsL * sin(theta1) / (NTCOIL * Q(NR))
             ! Poloidal angle at which the difference between width0 and width1 is minimized.
             if(i == 1) then
                theta_min = theta1
                diff_min = abs(width0 - width1)
             else
                if(abs(width0 - width1) < diff_min) then
                   diff_min = abs(width0 - width1)
                   theta_min = theta1
                end if
             end if
             ! Rim of ripple well region detected
             if(abs(width0 - width1) < 1.d-6) exit
             ! Overreached a rim of ripple well. Go back and use finer step size.
             if(width0 < width1) then
                theta1  = theta1 - PI * dlt
                dlt = 0.1d0 * dlt
                i = 0
                irip = irip - 1
                cycle
             end if
             ! Ripple amplitude at the rim of the ripple well region
             DltRP_rim = width0
             theta_rim = theta1
          end do
          ARC = 2.d0 * theta1
          th1(nr) = theta_min

          ! For HFS
          theta2 = PI
          i = 0
          imax = 101
          dlt = 1.d0 / (imax - 1)
          irip = 0
          do 
             i = i + 1
             irip = irip + 1
             if(i == imax) then
!                write(6,'(A,I3)') "HFS rim of ripple well region not detected at NR = ",NR
                theta2 = PI
                exit
             end if
             theta2 = theta2 - PI * dlt
             RL = R(NR) * (1.D0 + (kappa - 1.D0) * sin(theta2))
             width0 = ripple(RL,theta2,FSRP)
             width1 = EpsL * sin(theta2) / (NTCOIL * Q(NR))
             if(i == 1) then
                theta_min = theta2
                diff_min = abs(width0 - width1)
             else
                if(abs(width0 - width1) < diff_min) then
                   diff_min = abs(width0 - width1)
                   theta_min = theta2
                end if
             end if
             if(abs(width0 - width1) < 1.d-6) exit
             if(width0 < width1) then
                theta2  = theta2 + PI * dlt
                dlt = 0.1d0 * dlt
                i = 0
                irip = irip - 1
                cycle
             end if
          end do
          ARC = ARC + 2.d0 * (PI - theta2)
          ! Ratio of ripple well region in a certain flux surface
          rip_rat(NR) = ARC / (2.d0 * PI)
          th2(nr) = theta_min

!!rpl_ave          sum_rp = 0.d0
!!rpl_ave          imax = 51
!!rpl_ave          dlt = 1.d0 / (imax - 1)
!!rpl_ave          do i = 1, imax
!!rpl_ave             theta = (i - 1) * PI * dlt
!!rpl_ave             RL = R(NR) * (1.D0 + (kappa - 1.D0) * sin(theta))
!!rpl_ave             sum_rp = sum_rp + ripple(RL,theta,FSRP)**2
!!rpl_ave          end do
!!rpl_ave          DltRP_ave = sqrt(sum_rp/imax)

!!$          ! alpha_l : ripple well parameter
!!$          alpha_l = EpsL * sin(thetab) / (NTCOIL * Q(NR) * DltRP(NR))
!!$          ! Dlteff : effective depth of well along the magnetic field line
!!$          Dlteff = 2.D0*DltRP(NR)*(SQRT(1.D0-alpha_l**2)-alpha_l*acos(alpha_l))

          ! effective time of detrapping
          ! (Yushmanov NF (1982), Stringer NF (1972) 689, Takamura (5.31))
          rNubrp1(NR) = rNuD(NR) / DltRP_rim
          ! See the description of "Convective loss"
          rNubrp2(NR) = rNubrp1(NR) * SQRT(DltRP_rim) * fmaxwell(theta_rim,0.5D0*PI,0.85D0)
!!rpl_ave          rNubrp2(NR) = rNubrp1(NR) * SQRT(DltRP_ave)

          ! Convectitve loss (vertical grad B drift velocity)
          Vdrift = 0.5D0 * AMb * Vb**2 / (PZ * AEE * RR * SQRT(BphV(NR)**2 + BthV(NR)**2))
!          RUbrp(NR)=(NTCOIL*Q(NR)*RR*DltRP(NR))*Vdrift
!          if(nr/=0) Ubrp(NR)=(NTCOIL*Q(NR)*RR*DltRP(NR))/R(NR)*Vdrift
          RUbrp(NR)=(NTCOIL*Q(NR)*RR*DltRP_rim)*Vdrift
          if(nr/=0) Ubrp(NR)=(NTCOIL*Q(NR)*RR*DltRP_rim)/R(NR)*Vdrift
!!$          Ubrp(NR) = 0.5D0 * Vdrift
!!$            &  * (theta1*sin(theta1) + (PI - theta2)*sin(theta2)) / (PI + theta1 - theta2))
!!$          IF(NR == NRMAX) THEN
!!$             rNubL(NR) = rNubL(NR-1)
!!$          ELSE
!!$             rNubL(NR) = Ubrp(NR) / SQRT(RB**2 - R(NR)**2)
!!$          END IF
!!$          rNubL(NR) = Ubrp(NR) / (R(NR) * sin(theta1))
!          if(nr >=5) stop
       END DO
!!$       RV0 = AITKEN2P(R(0),r(1)*(pi-th2(1)),r(1)*th1(1),r(2)*th1(2),-R(1),R(1),R(2))
!!$       rNubL(0) = Ubrp(0) / RV0
       Ubrp(0) = AITKEN2P(R(0),Ubrp(1),Ubrp(2),Ubrp(3),R(1),R(2),R(3))
       ! On the axis ripple amplitude is uniquely defined because of no poloidal variation.
       rNubrp1(0) = rNuD(0) / DltRP(0)
       rNubrp2(0) = rNubrp1(0) * SQRT(DltRP(0))

       ! Save for graphic
       thrp(1:nrmax) = th2(nrmax:1:-1)
       thrp(nrmax+1:2*nrmax) = th1(1:nrmax)

!!rp_conv       CALL SPL1D(R,PNbrpV,DERIV,U,NRMAX+1,0,IER)
!!rp_conv       do nr = 0, nrmax
!!rp_conv          if(R(NR) <= RB * cos(th1(nr))) then
!!rp_conv             tmp = r(nr)*cos(th1(nr))
!!rp_conv             call wherenr(r,tmp,nrl,Left,Right)
!!rp_conv             CALL SPL1DF(tmp,PNbrpL(NR),R,U,NRMAX+1,IER)!
!!rp_conv             PNbrpLV(NRL-1) = PNbrpLV(NRL-1) + Left  * PNbrpL(NR)
!!rp_conv             PNbrpLV(NRL)   = PNbrpLV(NRL)   + Right * PNbrpL(NR)
!!rp_conv             write(6,*) nrl,sngl(tmp),sngl(r(nrl-1)),sngl(r(nrl)),sngl(PNbrpL(NR)),sngl(Left  * PNbrpL(NR)),sngl(Right * PNbrpL(NR)),sngl(PNbrpLV(NRL-1)),sngl(PNbrpLV(NRL))
!!rp_conv          end if
!!rp_conv       end do

       !  +++ Diffusive loss +++

       ! *** Physical aspect of banana particles suffered from diffusive loss ***
       ! Banana particles are not affected by local ripple wells regardless whether
       ! the wells exist or not. All the diffusion processes for them occur at the
       ! banana tip point "thetab", hence we only consider the ripple amplitude at
       ! the banana tip point, DltRP(NR).

       !  -- Collisional diffusion of trapped fast particles --
       IF(PNBH == 0.D0) THEN
          Dbrp(0:NRMAX) = 0.D0
       ELSE
       do nr = 1, nrmax
          RL = R(NR)
          EpsL = RL / RR

          ! rhob : Larmor radius of beam ions
          rhob = AMb * Vb / (PZ * AEE * SQRT(BphV(NR)**2 + BthV(NR)**2))
          ! DltR : Step size of banana particles
!          DltR = SQRT(PI * NTCOIL * (Q(NR) / EpsL)**3 * DltRP(NR)**2 * rhob**2)
          DltR = rhob*DltRP(NR)*SQRT(PI*NTCOIL*Q(nr)**3)*ABS(sinthb) &
               & /((NTCOIL*Q(NR)*DltRP(NR))**1.5D0+(EpsL*ABS(sinthb))**1.5D0)

          ! effective collisional frequency (G. Park et al., PoP 10 (2003) 4004, eq.(6))
!          rNueff = 1.82d0 * Q(NR)**2 * NTCOIL**2 / EpsL * rNuD(NR) ! for thetab=0.5*PI
          rNueff = 8.D0 * Q(NR)**2 * NTCOIL**2 * rNuD(NR) / (EpsL * sinthb**2) &
               & * (ellE/ellK - cos(0.5d0*thetab)**2)
          ! rNubnc : bounce frequency of beam ions (Helander and Sigmar, p132 eq.(7.27))
!          rNubnc = SQRT(EpsL) * Vb / (10.5D0 * Q(NR) * (RR + R(NR))) ! for thetab=0.5*PI
          rNubnc = SQRT(EpsL) * Vb / (4.D0 * SQRT(2.D0) * ellK * Q(NR) * RR)
!!$          ! DCB : confined banana diffusion coefficient
!!$          ! (V. Ya Goloborod'ko, et al., Physica Scripta T16 (1987) 46)
!!$          DCB = NTCOIL**2.25D0*Q(NR)**3.25D0*RR*rhob*DltRP(NR)**1.5d0*rNuD(NR)/EpsL**2.5D0
          ! DRP : ripple-plateau diffusion coefficient
          DRP = rNubnc * DltR**2

          ! Collisional ripple well diffusion
!!$          if (rNueff < rNubnc) then
!!$             Dbrp(NR) = DCB
!!$          else
!!$             Dbrp(NR) = DRP
!!$          end if
!!$!          Dbrp(NR) = DCB * DRP / (DCB + DRP)
          Dbrp(NR) = DRP/SQRT(1.D0+(rNubnc/rNueff*DltR*NTCOIL*dQdr(NR))**2)

          ! Dltcr : GWB criterion of stochastic diffusion at banana tip point
!          Dltcr = (EpsL / (PI * NTCOIL * Q(NR)))**1.5D0 / (rhob * dQdr(NR))
          Dltcr = DltRP(NR) / (DltR * NTCOIL * (2.D0 * thetab * dQdr(NR) &
          &     + 2.D0 * Q(NR) / R(NR) * cos(thetab) / sinthb))
!!$          ! Fraction of stochastic region occupied in a flux surface
!!$          theta1 = 0.d0
!!$          i = 0
!!$          dlt = 1.d0 / (imax - 1)
!!$          ist = -1
!!$          do
!!$             i = i + 1
!!$             RL = R(NR) * (1.D0 + (kappa - 1.D0) * sin(theta1))
!!$             if(ripple(RL,theta1,fsrp) > Dltcr) ist = ist + 1
!!$             theta1 = theta1 + PI * dlt
!!$             if(i == imax) exit
!!$          end do
!!$          ! Collisionless stochastic (ergodic) diffusion (whose value is 
!!$          ! equivalent to that of ripple-plateau diffusion)
!!$          if (DltRP(NR) > Dltcr) then
!!$             ! facST : Fraction of stochastic region in a flux surface
!!$             facST = dble(ist) / dble(imax - 1)
!!$             Dbrp(NR) = DRP * facST + Dbrp(NR) * (1.d0 - facST)
!!$          end if
          if (DltRP(NR) > Dltcr) Dbrp(NR) = DRP

!!$          ! Old version for display and comparison
!!$          DRP = rNubnc * PI * NTCOIL * (Q(NR) / EpsL)**3 * DltRP(NR)**2 * rhob**2
!!$          DCB = NTCOIL**2.25D0*Q(NR)**3.25D0*RR*rhob*DltRP(NR)**1.5d0*rNuD(NR)/EpsL**2.5D0
!!$          Dbrp(NR) = DCB * DRP / (DCB + DRP)
!!$          if(DltRP(NR) > Dltcr) Dbrp(NR) = DRP
!!$          write(6,*) r(nr)/ra,ft(NR)*rip_rat(NR)*Dbrp(NR)
       end do
       Dbrp(0) = AITKEN2P(R(0),Dbrp(1),Dbrp(2),Dbrp(3),R(1),R(2),R(3))

       ! Potato orbit effect
       !   Approaching to the magnetic axis, we reach the point where the banana
       !   particles changes themselves to the potato particles.
       !   Inside the point, the particle orbit is no longer changed.
       do nr = 0, nrmax
          ! potato width (Helander and Sigmar, p133)
          Rpotato = (Q(NR)**2*(AMb * Vb / (PZ * AEE * BphV(NR)))**2*RR)**(1.D0/3.D0)
          if(r(nr) > Rpotato) then
             nr_potato = nr - 1
             exit
          end if
       end do
       do nr = 0, nr_potato
          Dbrp(nr)  = Dbrp(nr_potato+1)
       end do
       END IF

       !     ***** Neoclassical toroidal viscosity (NTV) *****
       !      "rNuNTV" and "UastNC" are obtained from NTVcalc

!       CALL NTVcalc

    ELSE
       rNubrp1(0:NRMAX) = 0.D0
       rNubrp2(0:NRMAX) = 0.D0
       Ubrp(0:NRMAX) = 0.D0
       Dbrp(0:NRMAX) = 0.D0
       rNuNTV(0:NRMAX) = 0.D0
       UastNC(0:NRMAX) = 0.D0
    END IF

    deallocate(dQdr,dVebdr,dErdr,dBthdr,dTedr,dTidr,dPedr,dPidr,dpdr)

    RETURN
  END SUBROUTINE TXCALC

!***************************************************************
!
!   Heating deposition profile
!     Input : R0    : Deposition center (m)
!             RW    : Deposition width (m)
!             CHR   : Trapped NB, passing NB, RF or Virtual torque
!             PNBCD : Injection direction, optional
!     Output : S(0:NRMAX) : Deposition profile
!              SINT       : Normalization factor
!
!***************************************************************

  subroutine deposition_profile(S,SINT,R0,RW,CHR,PNBCD)
    use tx_commons, only : NRMAX, NRA, FSRP, R, RA, RB, PI, RR, &
         &                 AMb, Vb, PZ, AEE, BthV, Q, BphV
    use tx_interface, only : INTG_F
    real(8), intent(in)  :: R0, RW
    real(8), intent(in), optional :: PNBCD
    character(len=*), intent(in) :: CHR
    real(8), intent(out), dimension(0:NRMAX) :: S
    real(8), intent(out) :: SINT
    integer(4) :: nr
    real(8) :: EpsL, Rshift, Rpotato, rhop
    real(8) :: AITKEN2P

    if(CHR == 'Virtual') then
       S(0:NRA) = 1.D0 - (R(0:NRA) / RA)**2
       S(NRA+1:NRMAX) = 0.D0
    else if(CHR == 'RFe' .OR. CHR == 'RFi') then
       S(0:NRMAX) = EXP(- ((R(0:NRMAX) - R0) / RW)**2) * (1.D0 - (R(0:NRMAX) / RB)**4)
    else if(CHR == 'NB') then
       if(abs(FSRP) > 0.D0) then
          S(0:NRMAX) = EXP(- ((R(0:NRMAX) - R0) / RW)**2) * (1.D0 - (R(0:NRMAX) / RB)**4)
       else
          S(0:NRA) = EXP(- ((R(0:NRA) - R0) / RW)**2) * (1.D0 - (R(0:NRA) / RA)**4)
          S(NRA+1:NRMAX) = 0.D0
       end if
    else
       if(present(PNBCD) .EQV. .FALSE.) stop 'deposition_profile: input error!'
       do nr = 0, nrmax
          EpsL = R(NR) / RR
          if(nr /= 0) rhop = AMb * Vb / (PZ * AEE * BthV(NR)) ! poloidal Larmor radius
          if(CHR == 'NB_TRAP') then ! trapped particle
             ! potato width
             Rpotato = (Q(NR)**2*(AMb * Vb / (PZ * AEE * BphV(NR)))**2*RR)**(1.D0/3.D0)
             if(nr == 0) then
                Rshift = PNBCD * Rpotato ! potato particle
             else
                Rshift = PNBCD * MIN(SQRT(EpsL) * rhop, Rpotato) ! potato or banana particle
             end if
          else
             if (nr == 0) then ! passing particle
                Rshift = PNBCD * (AMb * Vb * Q(NR) / (PZ * AEE * BphV(NR)))
             else
                Rshift = PNBCD * (     EpsL  * rhop)
             end if
          end if

          if(abs(FSRP) > 0.D0) then
             S(NR) = EXP(- (((R(NR) + Rshift) - R0) / RW)**2) * (1.D0 - (R(NR) / RB)**4)
          else
             if(nr <= nra) then
                S(NR) = EXP(- (((R(NR) + Rshift) - R0) / RW)**2) * (1.D0 - (R(NR) / RA)**4)
             else
                S(NR) = 0.D0
             end if
          end if
       end do
!!$       if(CHR == 'NB_TRAP') then
!!$          if(PNBCD > 0.D0) then
!!$             S(0) = AITKEN2P(R(0),S(1),S(2),S(3),R(1),R(2),R(3))
!!$          else if(PNBCD < 0.D0) then
!!$             S(0) = 0.D0
!!$          end if
!!$       end if
    end if

    SINT = 2.D0 * PI * INTG_F(S)

  end subroutine deposition_profile

!***************************************************************
!
!   Shifting deposition profile loaded from files
!     Input  : kchar      : 'TRAP' or 'PASS'
!              direct     : direction of NBI, -1 or 1
!     In/Out : f(0:NRMAX) : Deposition center (m)
!
!***************************************************************

  subroutine shift_prof(f,kchar,direct)
    use tx_commons, only : NRMAX, R, RR, AMb, Vb, PZ, AEE, BthV, Q, BphV, PNBCD
    real(8), dimension(0:NRMAX), intent(inout) :: f
    character(len=4), intent(in) :: kchar
    real(8), intent(in) :: direct
    integer(4) :: nr, nrl
    real(8) :: EpsL, rhop, Rpotato, Rshift
    real(8), dimension(0:nrmax) :: r_shift, fl1, fl2
    real(8), dimension(:), allocatable :: r_alloc, f_alloc

    ! Shift the horizontal axis
    do nr = 0, nrmax
       EpsL = R(NR) / RR
       if(nr /= 0) rhop = AMb * Vb / (PZ * AEE * BthV(NR)) ! poloidal Larmor radius
       if(kchar == 'TRAP') then
          ! potato width
          Rpotato = (Q(NR)**2*(AMb * Vb / (PZ * AEE * BphV(NR)))**2*RR)**(1.D0/3.D0)
          if(nr == 0) then
             Rshift = direct * Rpotato ! potato particle
          else
             Rshift = direct * MIN(SQRT(EpsL) * rhop, Rpotato) ! potato or banana particle
          end if
       else if(kchar == 'PASS') then
          if (nr == 0) then ! passing particle
             Rshift = direct * (AMb * Vb * Q(NR) / (PZ * AEE * BphV(NR)))
          else
             Rshift = direct * (     EpsL  * rhop)
          end if
       else
          stop 'shift_prof: input error!'
       end if

       r_shift(nr) = r(nr) - Rshift
    end do

    fl1(0:nrmax) = 0.d0
    fl2(0:nrmax) = 0.d0

    ! Fold back at the magnetic axis
    if(r_shift(0) < 0.d0) then
       do nr = 1, nrmax
          if(r_shift(nr) > 0.d0) then
             nrl = nr - 1
             exit
          end if
       end do

       allocate(r_alloc(0:nrl),f_alloc(0:nrl))
       r_alloc(0:nrl) = abs(r_shift(nrl:0:-1))
       f_alloc(0:nrl) = f(nrl:0:-1)
       do nr = 0, nrl
          call aitken(r(nr),fl1(nr),r_alloc,f_alloc,2,nrl+1)
       end do
       deallocate(r_alloc,f_alloc)
    end if

    ! Interpolate
    do nr = 0, nrmax
       call aitken(r(nr),fl2(nr),r_shift,f,2,nrmax+1)
    end do

    ! Sum of both contributions 
    f(0:nrmax) = fl1(0:nrmax) + fl2(0:nrmax)
    
  end subroutine shift_prof

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

!***************************************************************
!
!   Ion-electron heating fraction
!     (Tokamaks 3rd, p250)
!
!   (input)
!     x     : fraction of energy
!
!***************************************************************

  pure real(8) function rate_to_ion(x) result(f)
    use tx_commons, only : PI
    real(8), intent(in) :: x

    if (x == 0.d0) then
       f = 1.d0
    else
       f = 1.d0 / x * (  1.d0 / 3.d0 * log((1.d0 - sqrt(x) + x) / (1.d0 + sqrt(x))**2) &
            &          + 2.d0 / sqrt(3.d0) * (atan((2.d0 * sqrt(x) - 1.d0) / sqrt(3.d0)) &
            &          + PI / 6.d0))
    end if

  end function Rate_to_ion

!***************************************************************
!
!   Ripple amplitude function
!     (Yushmanov review pp.122, 123)
!
!***************************************************************

  real(8) function ripple(RL,theta,FSRP) result(f)
    use tx_commons, only : RR, NTCOIL, DltRPn, RA
    real(8), intent(in) :: RL, theta, FSRP
    real(8) :: a, L0, Rmag0 = 2.4D0 ! specific value for JT-60U
    real(8) :: BESIN

    if(FSRP /= 0.D0) then
       L0 = RR - Rmag0
       a = sqrt((RL**2+L0**2+2.D0*RL*L0*cos(theta))*(RR-L0)/(RR+RL*cos(theta)))
       f = DltRPn * BESIN(0,NTCOIL/(RR-L0)*a)
    else
       f = 0.d0
    end if

  end function ripple

!***************************************************************
!
!   Maxwellian distribution
!
!   (input)
!     x     : position
!     mu    : average
!     sigma : standard deviation
!
!***************************************************************

  real(8) function fmaxwell(x,mu,sigma) result(f)
    real(8), intent(in) :: x, mu, sigma

    f = exp(- (x - mu)**2 / (2.d0 * sigma**2))

  end function fmaxwell

!***************************************************************
!
!   Alpha heating power
! 
!     T(d,n)4He : D + T -> 4He + n + 17.6 MeV
!     valid for 0.2 keV <= Ti <= 100 keV
!
!     << Bosch-Hale fusion reactivity model >>
!     (H.-S. Bosch and G.M. Hale, Nucl. Fusion 32 (1992) 611)
!     written by HONDA Mitsuru based on ITPA SSO Plans (2007/10/19)
!
!     Inputs (real*8): tikev  : Ion temperature [keV]
!                      den_D  : Deuterium density [10^20/m^3]
!                      den_T  : Tritium density [10^20/m^3]
!     Output (real*8): bosch_fusion  : Fusion reaction power [W/m^3]
!
!***************************************************************

  real(8) function bosch_fusion(tikev,den_D,den_T)
    real(8), intent(in) :: tikev, den_D, den_T
    real(8) :: c1,c2,c3,c4,c5,c6,c7,bg,mrcsq
    real(8) :: denDcgs,denTcgs,theta,sk,svdt

    data c1,c2,c3,c4,c5,c6,c7/1.17302d-9,1.51361d-2,7.51886d-2, &
         &     4.60643d-3,1.35000d-2,-1.06750d-4,1.36600d-5/
    data bg,mrcsq/34.3827d0,1.124656d6/

    denDcgs = den_D * 1.d-6 * 1.d20 ! nD [/cm^3]
    denTcgs = den_T * 1.d-6 * 1.d20 ! nT [/cm^3]
    !...Bosch-Hale formulation
    !     theta : Eq.(13)
    theta = tikev/(1.d0-((tikev*(c2+(tikev*(c4+tikev*c6)))) &
         &      /(1.d0+tikev*(c3+tikev*(c5+tikev*c7)))))
    !     sk : Eq.(14)
    sk    = (bg**2/(4.d0*theta))**0.333d0
    !     svdt : Eq.(12), [cm^3/s]
    svdt  = c1*theta*sqrt(sk/(mrcsq*tikev**3))*exp(-3.d0*sk)
    !     Fusion reaction power for D-T reaction
    !        Alpha particle energy : 3.5 [MeV] = 5.6d-13 [J]
    bosch_fusion = 5.6d-13*denDcgs*denTcgs*svdt*1.d6
    
  end function bosch_fusion

!***************************************************************
!
!   Rate coefficients for electron impact hydrogen (e+H) ionization process
!     valid for 1eV => 10^5eV
! 
!     (C.E. Singer et al., Comp. Phys. Comm. 49 (1988) 275, Eq. (2.9.5l))
!     (R.L. Freeman and E.M. Jones, Culham Laboratory Report, CLM-R-137 (May 1974)
!
!     Inputs (real*8): tekev : Electron temperature [keV]
!     Output (real*8): SiViz : Ionization maxwellian rate coefficient (m^3/s)
!
!     Internal (real*8): a : Table of Maxwellian rate coefficients in TABLE 3
!                                     ^^^^^^^^^^
!***************************************************************

  real(8) function SiViz(tekev)

    implicit none
    real(8), intent(in) :: tekev
    real(8) :: x
    real(8), dimension(0:6) :: a
    data a /-0.3173850D02, 0.1143818D02, -0.3833998D01, 0.7046692D0, &
         &  -0.7431486D-1, 0.4153749D-2, -0.9486967D-4/

    if(tekev < 1.d-3 .or. tekev > 1.d2) stop 'Function SiViz: Out of energy range.'
    x = log(tekev * 1.d3)
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
!     Output (real*8): SiVcx : Charge exchange maxwellian rate coefficient (m^3/s)
!
!     Internal (real*8): a : Table of Maxwellian rate coefficients in TABLE 3
!                                     ^^^^^^^^^^
!***************************************************************

  real(8) function SiVcx(tikev)

    implicit none
    real(8), intent(in) :: tikev
    real(8) :: x
    real(8), dimension(0:8) :: a
    data a /-0.1841757D02, 0.5282950D0, -0.2200477D0,   0.9750192D-1, &
         &  -0.1749183D-1, 0.4954298D-3, 0.2174910D-3, -0.2530205D-4, 0.8230751D-6/

    if(tikev < 1.d-3 .or. tikev > 1.d2) stop 'Function SiVcx: Out of energy range.'
    x = log(tikev * 1.d3)
    SiVcx = exp( a(0)+x*(a(1)+x*(a(2)+x*(a(3)+x*(a(4) &
         &      +x*(a(5)+x*(a(6)+x*(a(7)+x*a(8)))))))))*1.D-6

  end function SiVcx

!!$  real(8) function ripple(NR,theta,FSRP) result(f)
!!$    use tx_commons, only : RR, R, RA, NTCOIL
!!$    integer(4), intent(in) :: NR
!!$    real(8), intent(in) :: theta, FSRP
!!$    real(8) :: DIN = 0.2D0, DltRP0 = 0.015D0
!!$
!!$    if(FSRP /= 0.D0) then
!!$       f = DltRP0 * (       ((RR + R(NR) * cos(theta)) / (RR + RA))**(NTCOIL-1) &
!!$            &        + DIN *((RR - RA) / (RR + R(NR) * cos(theta)))**(NTCOIL+1))
!!$    else
!!$       f = 0.D0
!!$    end if
!!$       
!!$  end function ripple

!!rp_conv  ! Search minimum radial number NR satisfying R(NR) > X.
!!rp_conv
!!rp_conv  subroutine wherenr(R,X,NR,Left,Right)
!!rp_conv    real(8), dimension(0:NRMAX), intent(in) :: R
!!rp_conv    real(8), intent(in) :: X
!!rp_conv    integer, intent(out) :: NR
!!rp_conv    real(8), intent(out) :: Left, Right
!!rp_conv    integer :: NRL
!!rp_conv
!!rp_conv    if(X == 0.d0) then
!!rp_conv       NR = 1
!!rp_conv       Left  = 0.d0
!!rp_conv       Right = 0.d0
!!rp_conv       return
!!rp_conv    end if
!!rp_conv
!!rp_conv    do nrl = 1, nrmax
!!rp_conv       if(r(nrl) > x) then
!!rp_conv          NR = nrl
!!rp_conv          Right = (x - r(nr-1)) / (r(nr) - r(nr-1))
!!rp_conv          Left  = 1.d0 - Right
!!rp_conv          exit
!!rp_conv       end if
!!rp_conv    end do
!!rp_conv
!!rp_conv  end subroutine wherenr
end module tx_variables
