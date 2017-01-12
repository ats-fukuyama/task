!***************************************************************
!
!   Auxiliary heating system
!
!***************************************************************

module aux_system
  implicit none
  private
  real(8), save :: Vbabsmax, Scxb
  public :: txauxs, Vbabsmax, Scxb

contains

  subroutine txauxs

    ! SNB*(NR) : The birth rate of beam ions per unit volume
    ! MNB(NR)  : In essence MNB is similar to SNB, but MNB distinguishes itself
    !            from SNB in that MNB takes into account the tangential NB direction,
    !            the pitch between the NB chord and the field line and the contribution
    !            of the perpendicular NBI to the momentum input.

    use tx_commons
    use tx_interface, only : inexpolate, fgaussian, coulog
    use tx_core_module, only : intg_vol

    INTEGER(4) :: NR, i, ideriv = 1, nrbound
    REAL(8) :: SL, SLT1, SLT2, PNBP0, PNBT10, PNBT20, PNBex0, SNBPDi_INTG, &
         &     PNBPi0, PNBTi10, PNBTi20, PRFe0, PRFi0, SL1, SL2, &
         &     Vti, rNubes, BBL, PALFL, Ecr, &
         &     EbL, logEbL, Sion, xl, Tqt0L
    real(8), parameter :: PAHe = 4.D0, & ! Atomic mass number of He
         &                Enf  = 3.5D3   ! in keV, equal to 3.5 MeV
    real(8), dimension(0:NRMAX) :: SNBP, SNBT1, SNBT2, SNBTi1, SNBTi2, SRFe, SRFi
    ! Mainly for derivatives
    real(8), dimension(:), allocatable :: Tqt_tmp, Tqp_tmp

    !     *** Constants ***

    !     Averaged beam injection energy

    Eb =  (esps(1) + 0.5d0 * esps(2) + esps(3) / 3.d0) * Ebmax

    !     NBI beam speed

    Vb =  SQRT(2.D0 * Eb * rKilo / (amb * amqp) )
    Vbpara(0:NRMAX) = Vb
    Vbabsmax = Vb


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
    IF(PNBHP /= 0.D0) THEN
       CALL deposition_profile(SNBP,SL,RNBP0,RNBP,'NB')
       PNBP0 = PNBHP * 1.D6 / SL
       PNBPi0 = PNBP0
       SNBPDi(0:NRMAX) = SNBP(0:NRMAX)

       IF(MDLNBD /= 0) THEN ! Orbit effect
          ! For ions
          IF(PNBMPD == 0.D0) THEN ! Exact perpendicular NBI
             CALL deposition_profile(SNBPDi,SL,RNBP0,RNBP,'NB_TRAP',0.D0)
          ELSE ! Near perpendicular NBI
             CALL deposition_profile(SNBPDi,SL,RNBP0,RNBP,'NB_TRAP',SIGN(1.D0,PNBMPD))
          END IF
          PNBPi0 = PNBHP * 1.D6 / SL
       END IF
    ELSE
       PNBP0  = 0.D0
       PNBPi0 = 0.D0
       SNBP(0:NRMAX)   = 0.D0
       SNBPDi(0:NRMAX) = 0.D0
    END IF

    !  *** Tangential
    IF(PNBHT1 /= 0.D0) THEN
       CALL deposition_profile(SNBT1,SLT1,RNBT10,RNBT1,'NB')
       PNBT10 = PNBHT1 * 1.D6 / SLT1

       IF(MDLNBD > 1) THEN ! Orbit effect
          ! For ions
          CALL deposition_profile(SNBTi1,SLT1,RNBT10,RNBT1,'NB_PASS',SIGN(1.D0,PNBCD))
          PNBTi10 = PNBHT1 * 1.D6 / SLT1
       END IF
    ELSE
       PNBT10 = 0.D0
       SNBT1(0:NRMAX)  = 0.D0
       SNBTi1(0:NRMAX) = 0.D0
    END IF

    IF(PNBHT2 /= 0.D0) THEN
       CALL deposition_profile(SNBT2,SLT2,RNBT20,RNBT2,'NB')
       PNBT20 = PNBHT2 * 1.D6 / SLT2

       IF(MDLNBD > 1) THEN
          ! For ions
          CALL deposition_profile(SNBTi2,SLT2,RNBT20,RNBT2,'NB_PASS',SIGN(1.D0,PNBCD))
          PNBTi20 = PNBHT2 * 1.D6 / SLT2
       END IF
    ELSE
       PNBT20 = 0.D0
       SNBT2(0:NRMAX)  = 0.D0
       SNBTi2(0:NRMAX) = 0.D0
    END IF

    !  For RF heating
    IF(PRFHe /= 0.D0) THEN
       CALL deposition_profile(SRFe,SL,RRFe0,RRFew,'RFe')
       PRFe0 = PRFHe * 1.D6 / SL
    ELSE
       PRFe0 = 0.D0
       SRFe(0:NRMAX) = 0.D0
    END IF

    IF(PRFHi /= 0.D0) THEN
       CALL deposition_profile(SRFi,SL,RRFi0,RRFiw,'RFi')
       PRFi0 = PRFHi * 1.D6 / SL
    ELSE
       PRFi0 = 0.D0
       SRFi(0:NRMAX) = 0.D0
    END IF

    ! Deposition profiles are loaded from the file

     !   In case of "NBI input from OFMC (1)", sequence of data is already
     !   defined as follows:
     !     (1) S_birth_ele,  (2) S_birth_tot, (3) S_birth_trap, (4) S_birth_pass,
     !     (5) S_birth_loss, (6) S_orbit_tot, (7) S_orbit_trap, (8) S_orbit_pass

     !     (1) S_birth_ele,   (2) S_birth_trap, (3) S_birth_pass
     !     (4) S_orbit_total, (5) S_orbit_trap, (6) S_orbit_pass

    ! *** Pre-defined input (OFMC: OrbitEffectDist.dat) ***************
    if(iflag_file == 1) then

       ! (1) Birth electrons
       i = 1
       call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBe,ideriv)
       ! Calibration by using total power of electrons
       SL = intg_vol(SNBe)
       if(SL /= 0.d0) SNBe(0:NRMAX) = SNBe(0:NRMAX) * 1.D-20 &
            &        * (infiles(i)%totP / (Eb * rKeV * SL))

       ! (2) Birth TOTAL ions (SNB for heating profiles)
       i = 2
       call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBi,ideriv)
       ! Calibration by using total power of all ions
       SL = intg_vol(SNBi)
       PNBHex = infiles(i)%totP * 1.D-6
       if(SL /= 0.d0) SNBi(0:NRMAX) = SNBi(0:NRMAX) * 1.D-20 &
            &        * (infiles(i)%totP / (Eb * rKeV * SL))
       PNBex0 = Eb * rKeV * 1.D20
       !  "or"= infiles(i)%totP / intg_vol(SNBi)
       ! Birth profiles for heating power
       !   (plasma is usually heated by beam ions, not beam electrons)
       SNB(0:NRMAX) = SNBi(0:NRMAX)

       ! *** No orbit effect for all ions ***
       if(MDLNBD == 0) then
          ! (3) Birth Trapped (SNBPDi for trapped beam ions)
          i = 3
          call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,2,SNBPDi)
          ! Calibration by using total power of all ions
          SL = intg_vol(SNBPDi)
          if(SL /= 0.d0) SNBPDi(0:NRMAX) = SNBPDi(0:NRMAX) * 1.D-20 &
               &          * (infiles(i)%totP / (Eb * rKeV * SL))
          PNBPi0 = Eb * rKeV * 1.D20
          ! "or" = infiles(2)%totP / intg_vol(SNBPDi)

          ! (4) Birth Passing (SNBTGi for passing beam ions)
          i = 4
          call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBTGi,ideriv)
          ! Calibration by using total power of trapped ions
          SL = intg_vol(SNBTGi)
          if(SL /= 0.d0) SNBTGi(0:NRMAX) = SNBTGi(0:NRMAX) * 1.D-20 &
               &          * (infiles(i)%totP / (Eb * rKeV * SL))
          PNBT10 = Eb * rKeV * 1.D20
          ! "or" = infiles(i)%totP / intg_vol(SNBTGi)

          ! (2) Birth TOTAL (SNBb for beam ions)
          SNBb(0:NRMAX) = SNB(0:NRMAX)

       ! *** Orbit effect ***
       else
          ! (7) Orbit Trapped (SNBPDi for trapped beam ions)
          i = 7
          call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBPDi,ideriv)
          ! Calibration by using total power of trapped ions
          SL = intg_vol(SNBPDi)
          if(SL /= 0.d0) SNBPDi(0:NRMAX) = SNBPDi(0:NRMAX) * 1.D-20 &
               &          * (infiles(i)%totP / (Eb * rKeV * SL))
          PNBPi0 = Eb * rKeV * 1.D20
          ! "or" = infiles(i)%totP / intg_vol(SNBPDi)

       ! *** Orbit effect for banana ions only ***
          if(MDLNBD == 1) then
             ! (4) Birth Passing (SNBTGi for passing beam ions)
             i = 4
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBTGi,ideriv)
             ! Calibration by using total power of trapped ions
             SL = intg_vol(SNBTGi)
             if(SL /= 0.d0) SNBTGi(0:NRMAX) = SNBTGi(0:NRMAX) * 1.D-20 &
                  &          * (infiles(i)%totP / (Eb * rKeV * SL))
             PNBT10 = Eb * rKeV * 1.D20
             ! "or" = infiles(i)%totP / intg_vol(SNBTGi)

             ! (4) Birth Passing + (7) Orbit Trapped (SNBb for beam ions)
             SNBb(0:NRMAX) = SNBTGi(0:NRMAX) + SNBPDi(0:NRMAX)

       ! *** Orbit effect for all ions ***
          else if(MDLNBD == 2) then
             ! (6) Orbit TOTAL (SNBb for beam ions)
             i = 6
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBb,ideriv)
             ! Calibration by using total power of trapped ions
             SL = intg_vol(SNBb)
             if(SL /= 0.d0) SNBb(0:NRMAX) = SNBb(0:NRMAX) * 1.D-20 &
                  &        * (infiles(i)%totP / (Eb * rKeV * SL))

             ! (8) Orbit Passing (SNBTGi for passing beam ions)
             i = 8
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBTGi,ideriv)
             ! Calibration by using total power of trapped ions
             SL = intg_vol(SNBTGi)
             if(SL /= 0.d0) SNBTGi(0:NRMAX) = SNBTGi(0:NRMAX) * 1.D-20 &
                  &          * (infiles(i)%totP / (Eb * rKeV * SL))
             PNBT10 = Eb * rKeV * 1.D20
             ! "or" = infiles(i)%totP / intg_vol(SNBTGi)
          end if
       end if
       ! Collisional torque injection part
       MNB(0:NRMAX)  = PNBCD * SNBTGi(0:NRMAX) * PNBPTC + SNBPDi(NR) * PNBMPD

       ! Local parallel velocity at birth for passing ions
       i = 8
       call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%vb,NRMAX,RHO,5,Vbpara,nrbound=nrbound,idx=0)
       do nr = nrbound+1, nrmax
          Vbpara(nr) = Vbpara(nrbound)
       end do
       Vbabsmax = maxval(abs(Vbpara))

    ! *** Pre-defined input (OFMC: Torque.txt) ************************
    else if(iflag_file == 2) then

       allocate(Tqt_tmp(0:NRMAX))
       ! Total torque input
       i = 1
       call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,Tqt_tmp,ideriv,idx=0)
       SL = intg_vol(Tqt_tmp)
       if(SL /= 0.d0) Tqt_tmp(0:NRMAX) = Tqt_tmp(0:NRMAX) * (infiles(i)%totS / SL)

    ! *** Arbitrary input *********************************************
    else if(iflag_file == 3) then
       do i = 1, n_infiles
          if(infiles(i)%name == datatype(1)) then ! Perp NB
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBP)
             SL   = intg_vol(SNBP)
             PNBP0  = PNBHP * 1.D6 / SL
             ! Banana orbit effect
             if(MDLNBD /= 0) then! .and. PNBMPD /= 0.D0) then
                do NR = 0, NRMAX
                   ! Passing particles generated by Perp NB
                   SNBTi2(NR) = (1.d0 - ft(NR)) * SNBP(NR) ! tangential part of ions
                   SNBT2(NR)  = SNBTi2(NR)                 ! tangential part of electrons
                   ! Trapped particles generated by Perp NB
                   SNBPDi(NR) = ft(NR) * SNBP(NR)          ! perpendicular part of ions
                   SNBP(NR)   = SNBPDi(NR)                 ! perpendicular part of electrons
                end do
                PNBTi20 = PNBP0 ! tangential part of ions
                PNBT20  = PNBP0 ! tangential part of electrons
                SNBPDi_INTG = intg_vol(SNBPDi)
                call shift_prof(SNBPDi, 'TRAP',SIGN(1.D0,PNBCD))
                ! calibration of Perp NB amplitude of ions
                PNBPi0 = PNBP0 * (SNBPDi_INTG / intg_vol(SNBPDi))
!!                call shift_prof(SNBTi2,'PASS',SIGN(1.D0,PNBCD))
             end if
          else if(infiles(i)%name == datatype(2)) then ! Tang NB 1
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,9,SNBT1)
             SLT1 = intg_vol(SNBT1)
             PNBT10 = PNBHT1 * 1.D6 / SLT1
          else if(infiles(i)%name == datatype(3)) then ! Tang NB 2
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,9,SNBT2)
             SLT2 = intg_vol(SNBT2)
             PNBT20 = PNBHT2 * 1.D6 / SLT2
          else if(infiles(i)%name == datatype(4)) then ! RF
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,0,SRFe)
             SL   = intg_vol(SRFe)
             PRFe0 = PRFHe * 1.D6 / SL
             SL   = intg_vol(SRFe)
             PRFi0 = PRFHi * 1.D6 / SL
          end if
       end do
    end if

    !   NBI total input power (MW)
    PNBH = PNBHP + PNBHT1 + PNBHT2 + PNBHex

    !   Ratio of CX deposition rate to IZ deposition rate
    !     (Riviere, NF 11 (1971) 363)

    EbL = Eb * 1.D3 / amb
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
       PALFL = FSNF * bosch_fusion(Var(NR,2)%T,0.5D0*Var(NR,2)%n,0.5D0*Var(NR,2)%n)
       Ecr = 14.8D0 * (PAHe / amb**(2.D0/3.D0)) * Var(NR,1)%T ! in keV
       PALFi(NR) = PALFL * rate_to_ion(Enf/Ecr)
       PALFe(NR) = PALFL - PALFi(NR)
    end do

    !   Additional torque input to LQi3 without net torque

    if(Tqp0 /= 0.d0) then
       allocate(Tqp_tmp(0:NRMAX))
       do nr = 0, nrmax
          Tqp(nr) = Tqp0 * fgaussian(Rho(nr),0.6d0,0.02d0)
       end do
       SL1 = intg_vol(Tqp)
       do nr = 0, nrmax
          Tqp_tmp(nr) = Tqp0 * fgaussian(Rho(nr),0.7d0,0.02d0)
       end do
       SL2 = intg_vol(Tqp_tmp)
       do nr = 0, nrmax
          Tqp(nr) = Tqp(nr) - (SL1 / SL2) * Tqp_tmp(nr)
       end do
       deallocate(Tqp_tmp)
    else
       Tqp(0:NRMAX) = 0.d0
    end if

    !   Additional torque input to LQi4

!    CALL deposition_profile(Tqt,SL,0.d0,0.d0,'Additional')
    CALL deposition_profile(Tqt,SL,RNBT10,RNBT1,'Additional')
    Tqt0L = Tqt0 / SL ! Tqt0 [N m], Tqt0L [-]
    Tqt(0:NRMAX) = Tqt0L * Tqt(0:NRMAX) ! <R.Source> [Nm/m^3]
!    Tqt(0:NRMAX) = Tqt0L * Tqt(0:NRMAX) * vlt(0:NRMAX) / vlt(NRMAX)!test1
!    Tqt(0:NRA) = Tqt0 * EXP(- ((vlt(0:NRA)/vlt(NRA) - RNBT10) / RNBT1)**2) * (1.D0 - (vlt(0:NRA)/vlt(NRA))**4)!test2
!    Tqt(NRA+1:NRMAX) = 0.d0
    if(iflag_file == 2) then
       Tqt(0:NRMAX) = Tqt(0:NRMAX) + Tqt_tmp(0:NRMAX)
       deallocate(Tqt_tmp)
    end if

    ! ************** Heating part end **************

    !     *** Beam slowing down time (momentum transfer with beam) ***
    ! reference : memo (92/04/02, 92/04/21)
    !             Tokamaks 3rd pp.246 - 252

    do NR = 0, NRMAX
       BBL = SQRT(BphV(NR)**2 + BthV(NR)**2)

!!!       Vcr = (3.D0 * SQRT(PI / 2.D0) * Var(NR,2)%n * achgb**2 / Var(NR,1)%n * amas(1) / amas(2) &
!!!            &   * (ABS(Var(NR,1)%T) * rKilo / (amas(1) * amqp))**1.5D0)**(1.D0/3.D0)
!       Ecr = (9.D0 * PI / 16.D0 * amas(2) / amas(1) )**(1.D0/3.D0) * amb / amas(2) * Var(NR,1)%T ! in keV
       Ecr = 14.8D0 * (amb / amb**(2.D0/3.D0)) * Var(NR,1)%T ! in keV
       Vti = SQRT(2.D0 * ABS(Var(NR,2)%T) * rKilo / (amas(2) * amqp))
       PNBcol_i(NR) = rate_to_ion(Eb/Ecr)
       PNBcol_e(NR) = 1.d0 - PNBcol_i(NR)
       IF(PNBH == 0.D0 .AND. PNbV(NR) < 1.D-8) THEN
          rNuD (NR) = 0.D0
          rNuB (NR) = 0.D0
       ELSE
          !     *** deflection time of beam ions against bulk ions ***
          !     (Takamura (3.26) + Tokamaks 3rd p64)
          xl = Vb / Vti
          rNuD(NR) = Var(NR,2)%n *1.D20 * achg(2)**2 * achgb**2 * AEE**4 &
               &   * coulog(1.d0,Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amas(2),achg(2),amb,achgb,PTbV(NR)) &
               &   / (2.D0 * PI * EPS0**2 * (amb*amp)**2 * Vb**3) &
               &   * (  SQRT(1.D0 - EXP(- 4.D0 * xl**2 / PI)) &
               &      - 2.D0 * xl / (4.D0 * xl**3 + 3.D0 * SQRT(PI)))

          !     *** Beam ion-electron slowing-down time
          !     (Tokamaks 3rd, A_D from below (2.13.2) + (2.14.1), assuming AMB >> amas(1)
          !      and Vb/(SQRT(2)*Vte) << 1)
          !     (originally from the book by Spitzer (1962))
          rNubes = Var(NR,1)%n * 1.D20 * achg(1)**2 * achgb**2 * AEE**4 &
               & * coulog(1.d0,Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amb,achgb,amas(1),achg(1)) &
!               & * coulog(1.d0,Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amb,achgb,amb,achgb) & ! original, but could be wrong
               & / (6.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * amb * amas(1) * amp**2 &
               & * (ABS(Var(NR,1)%T) * rKilo / (amas(1) * amqp))**1.5D0)

          ! The definition given below is a "energy slowing down time" and is not
          ! "particle slowing down time". At present, we assume the former is the
          ! same as the latter.
          rNuB(NR) = rNubes * 3.D0 / LOG(1.D0 + (Eb / Ecr)**1.5D0)
       END IF

       !     *** Heating profile ***

       ! For graphic
       if(iflag_file == 1) then
          PNBPD(NR) = PNBPi0 * SNBPDi(NR)  ! Power of trapped ions
          PNBTG(NR) = PNBT10 * SNBTGi(NR)  ! Power of passing ions
       else
          PNBPD(NR) = PNBP0 * SNBP(NR)     ! Power of perpendicular NBI
          PNBTG(NR) = PNBT10 * SNBT1(NR) + PNBT20 * SNBT2(NR)  ! Power of tangential NBIs
       end if
       PNB(NR)   = PNBPD(NR) + PNBTG(NR)

       !     *** For graphic and calculation ***

       !   Note: in case of iflag_file == 1, following terms have been already defined above.
       if(iflag_file /= 1) then
          ! Source profile for passing ions in temporal and graphic use
          if(MDLNBD == 2) then
             SNBTGi(NR)=(PNBTi10 * SNBTi1(NR) + PNBTi20 * SNBTi2(NR)) / (Eb * rKeV * 1.D20)
          else
             SNBTGi(NR)=(PNBT10  * SNBT1(NR)  + PNBT20  * SNBT2(NR))  / (Eb * rKeV * 1.D20)
          end if
          ! Source profile for trapped ions with banana orbit effect
          !   in temporal and graphic use
          SNBPDi(NR)= PNBPi0 * SNBPDi(NR) / (Eb * rKeV * 1.D20)
          ! Birth profiles for heating power
          SNB(NR)   = PNB(NR) / (Eb * rKeV * 1.D20)
          ! Birth profiles for electrons and thermal ions
          SNBe(NR)  = SNB(NR)
          SNBi(NR)  = SNB(NR)
!!old fashion          SNBi(NR)  = SNBPDi(NR) + SNBTGi(NR)
          ! Source profiles for beam ions with banana orbit effect
          SNBb(NR)  = SNBPDi(NR) + SNBTGi(NR)
          ! Torque injection part
          MNB(NR)   = PNBCD * SNBTGi(NR) * PNBPTC + SNBPDi(NR) * PNBMPD
       end if

       PRFe(NR)  = PRFe0 * SRFe(NR)
       PRFi(NR)  = PRFi0 * SRFi(NR)

       !     *** NBI power deposition for graphics ***
       !     see eq. (9.7) of [Hirshman and Sigmar, NF 22 (1981) 1079] 

       PNBe(NR) = Eb * SNB(NR) * PNBcol_e(NR) * (1.D20 * rKeV)
       PNBi(NR) = Eb * SNB(NR) * PNBcol_i(NR) * (1.D20 * rKeV) &
            &   + amb*amp * MNB(NR) * Var(NR,2)%BUpar / bbt(NR) * BUbparV(NR) * 1.D20

       !     *** Parallel torque injection ***
       !     BSmb = m_b dot{n}_b <B v_//0>
       BSmb(NR) = (amb * amp) * MNB(NR) * (sqrt(bbt(NR)) * Vbpara(NR))
    end do

  end subroutine txauxs

!***************************************************************
!
!   Heating deposition profile
!     Input : R0    : Deposition center (-)
!             RW    : Deposition width (-)
!             CHR   : Trapped NB, passing NB, RF or Additional torque
!             PNBCD : Injection direction, optional
!     Output : S(0:NRMAX) : Deposition profile
!              SINT       : Normalization factor (volume integrated)
!
!***************************************************************

  subroutine deposition_profile(S,SINT,R0,RW,CHR,PNBCD)
    use tx_commons, only : NRMAX, NRA, FSRP, RA, PI, RR, &
         &                 amb, amqp, Vb, achgb, BthV, Q, BphV, rho, epst
    use tx_core_module, only : intg_vol
    real(8), intent(in)  :: R0, RW
    real(8), intent(in), optional :: PNBCD
    character(len=*), intent(in) :: CHR
    real(8), intent(out), dimension(0:NRMAX) :: S
    real(8), intent(out) :: SINT
    integer(4) :: nr
    real(8) :: Rshift, Rpotato, rhop
!!    real(8) :: AITKEN2P

    if(CHR == 'Additional') then
!       S(0:NRA) = 1.D0 - rho(0:NRA)**2
       S(0:NRA) = EXP(- ((rho(0:NRA) - R0) / RW)**2) * (1.D0 - rho(0:NRA)**4)
       S(NRA+1:NRMAX) = 0.D0
    else if(CHR == 'RFe' .OR. CHR == 'RFi') then
       S(0:NRMAX) = EXP(- ((rho(0:NRMAX) - R0) / RW)**2) * (1.D0 - (rho(0:NRMAX)/rho(nrmax))**4)
    else if(CHR == 'NB') then
       if(abs(FSRP) > 0.D0) then
          S(0:NRMAX) = EXP(- ((rho(0:NRMAX) - R0) / RW)**2) * (1.D0 - (rho(0:NRMAX)/rho(nrmax))**4)
       else
          S(0:NRA) = EXP(- ((rho(0:NRA) - R0) / RW)**2) * (1.D0 - rho(0:NRA)**4)
          S(NRA+1:NRMAX) = 0.D0
       end if
    else
       if(present(PNBCD) .EQV. .FALSE.) stop 'deposition_profile: input error!'
       do nr = 0, nrmax
          if(nr /= 0) rhop = amb * Vb * amqp / (achgb * BthV(NR)) ! poloidal Larmor radius
          if(CHR == 'NB_TRAP') then ! trapped particle
             ! potato width
             Rpotato = (Q(NR)**2*(amb * Vb * amqp / (achgb * BphV(NR)))**2*RR)**(1.D0/3.D0)
             if(nr == 0) then
                Rshift = PNBCD * Rpotato ! potato particle
             else
                Rshift = PNBCD * MIN(SQRT(epst(nr)) * rhop, Rpotato) ! potato or banana particle
             end if
          else
             if (nr == 0) then ! passing particle
                Rshift = PNBCD * (amb * Vb * Q(NR) * amqp / (achgb * BphV(NR)))
             else
                Rshift = PNBCD * (     epst(nr)  * rhop)
             end if
          end if

          if(abs(FSRP) > 0.D0) then
             S(NR) = EXP(- (((rho(nr) + Rshift/RA) - R0) / RW)**2) * (1.D0 - (rho(nr)/rho(nrmax))**4)
          else
             if(nr <= nra) then
                S(NR) = EXP(- (((rho(nr) + Rshift/RA) - R0) / RW)**2) * (1.D0 - rho(NR)**4)
             else
                S(NR) = 0.D0
             end if
          end if
       end do
!!$       if(CHR == 'NB_TRAP') then
!!$          if(PNBCD > 0.D0) then
!!$             S(0) = AITKEN2P(rho(0),S(1),S(2),S(3),rho(1),rho(2),rho(3))
!!$          else if(PNBCD < 0.D0) then
!!$             S(0) = 0.D0
!!$          end if
!!$       end if
    end if

    ! Modify S(0) and S(1) so that S' becomes zero at the axis
    call sctr(rho(1),rho(2),rho(3),rho(4),S(2),S(3),S(4),S(0),S(1))

    SINT = intg_vol(S)

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
    use tx_commons, only : NRMAX, RR, amb, amqp, Vb, achgb, BthV, Q, BphV, epst, r
    real(8), dimension(0:NRMAX), intent(inout) :: f
    character(len=4), intent(in) :: kchar
    real(8), intent(in) :: direct
    integer(4) :: nr, nrl
    real(8) :: rhop, Rpotato, Rshift
    real(8), dimension(0:nrmax) :: r_shift, fl1, fl2
    real(8), dimension(:), allocatable :: r_alloc, f_alloc

    ! Shift the horizontal axis
    do nr = 0, nrmax
       if(nr /= 0) rhop = amb * Vb * amqp / (achgb * BthV(NR)) ! poloidal Larmor radius
       if(kchar == 'TRAP') then
          ! potato width
          Rpotato = (Q(NR)**2*(amb * Vb * amqp / (achgb * BphV(NR)))**2*RR)**(1.D0/3.D0)
          if(nr == 0) then
             Rshift = direct * Rpotato ! potato particle
          else
             Rshift = direct * MIN(SQRT(epst(nr)) * rhop, Rpotato) ! potato or banana particle
          end if
       else if(kchar == 'PASS') then
          if (nr == 0) then ! passing particle
             Rshift = direct * (amb * Vb * Q(NR) * amqp / (achgb * BphV(NR)))
          else
             Rshift = direct * (     epst(nr)  * rhop)
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

end module aux_system
