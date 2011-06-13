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

    use tx_commons
    use tx_interface, only : INTG_F, inexpolate, fgaussian, coulog, CORR

    INTEGER(4) :: NR, i, ideriv = 1, nrbound
    REAL(8) :: SL, SLT1, SLT2, PNBP0, PNBT10, PNBT20, PNBex0, SNBPDi_INTG, &
         &     PNBPi0, PNBTi10, PNBTi20, PRFe0, PRFi0, SL1, SL2, &
         &     Vti, rNuPara, rNubes, BBL, PALFL, Ecr, &
         &     EbL, logEbL, Sion, xl, Tqt0L, Tb
    real(8), parameter :: PAHe = 4.D0, & ! Atomic mass number of He
         &                Enf  = 3.5D3   ! in keV, equal to 3.5 MeV
    real(8), dimension(0:NRMAX) :: SNBP, SNBT1, SNBT2, SNBTi1, SNBTi2, SRFe, SRFi
    ! Mainly for derivatives
    real(8), dimension(:), allocatable :: Tqt_tmp, Tqp_tmp

    !     *** Constants ***

    !     NBI beam speed

    Vb =  SQRT(2.D0 * Eb * rKeV / AMB)
    Tb =  2.d0 / 3.d0 * Eb ! Mean temperature of beam ions
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
       PNBP0 = PNBHP * 1.D6 / (2.D0 * Pi * RR * SL)
       PNBPi0 = PNBP0
       SNBPDi(0:NRMAX) = SNBP(0:NRMAX)

       IF(MDLNBD /= 0) THEN ! Orbit effect
          ! For ions
          IF(PNBMPD == 0.D0) THEN ! Exact perpendicular NBI
             CALL deposition_profile(SNBPDi,SL,RNBP0,RNBP,'NB_TRAP',0.D0)
          ELSE ! Near perpendicular NBI
             CALL deposition_profile(SNBPDi,SL,RNBP0,RNBP,'NB_TRAP',SIGN(1.D0,PNBMPD))
          END IF
          PNBPi0 = PNBHP * 1.D6 / (2.D0 * Pi * RR * SL)
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
       PNBT10 = PNBHT1 * 1.D6 / (2.D0 * Pi * RR * SLT1)

       IF(MDLNBD > 1) THEN ! Orbit effect
          ! For ions
          CALL deposition_profile(SNBTi1,SLT1,RNBT10,RNBT1,'NB_PASS',SIGN(1.D0,PNBCD))
          PNBTi10 = PNBHT1 * 1.D6 / (2.D0 * Pi * RR * SLT1)
       END IF
    ELSE
       PNBT10 = 0.D0
       SNBT1(0:NRMAX)  = 0.D0
       SNBTi1(0:NRMAX) = 0.D0
    END IF

    IF(PNBHT2 /= 0.D0) THEN
       CALL deposition_profile(SNBT2,SLT2,RNBT20,RNBT2,'NB')
       PNBT20 = PNBHT2 * 1.D6 / (2.D0 * Pi * RR * SLT2)

       IF(MDLNBD > 1) THEN
          ! For ions
          CALL deposition_profile(SNBTi2,SLT2,RNBT20,RNBT2,'NB_PASS',SIGN(1.D0,PNBCD))
          PNBTi20 = PNBHT2 * 1.D6 / (2.D0 * Pi * RR * SLT2)
       END IF
    ELSE
       PNBT20 = 0.D0
       SNBT2(0:NRMAX)  = 0.D0
       SNBTi2(0:NRMAX) = 0.D0
    END IF

    !  For RF heating
    IF(PRFHe /= 0.D0) THEN
       CALL deposition_profile(SRFe,SL,RRFe0,RRFew,'RFe')
       PRFe0 = PRFHe * 1.D6 / (2.D0 * Pi * RR * SL)
    ELSE
       PRFe0 = 0.D0
       SRFe(0:NRMAX) = 0.D0
    END IF

    IF(PRFHi /= 0.D0) THEN
       CALL deposition_profile(SRFi,SL,RRFi0,RRFiw,'RFi')
       PRFi0 = PRFHi * 1.D6 / (2.D0 * Pi * RR * SL)
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
       SL = 2.D0 * Pi * INTG_F(SNBe)
       if(SL /= 0.d0) SNBe(0:NRMAX) = SNBe(0:NRMAX) * 1.D-20 &
            &        * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))

       ! (2) Birth TOTAL ions (SNB for heating profiles)
       i = 2
       call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBi,ideriv)
       ! Calibration by using total power of all ions
       SL = 2.D0 * Pi * INTG_F(SNBi)
       PNBHex = infiles(i)%totP * 1.D-6
       if(SL /= 0.d0) SNBi(0:NRMAX) = SNBi(0:NRMAX) * 1.D-20 &
            &        * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))
       PNBex0 = Eb * rKeV * 1.D20
       !  "or"= infiles(i)%totP / (2.D0 * Pi * RR * (2.D0 * Pi * INTG_F(SNBi)))
       ! Birth profiles for heating power
       !   (plasma is usually heated by beam ions, not beam electrons)
       SNB(0:NRMAX) = SNBi(0:NRMAX)

       ! *** No orbit effect for all ions ***
       if(MDLNBD == 0) then
          ! (3) Birth Trapped (SNBPDi for trapped beam ions)
          i = 3
          call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,2,SNBPDi)
          ! Calibration by using total power of all ions
          SL = 2.D0 * Pi * INTG_F(SNBPDi)
          if(SL /= 0.d0) SNBPDi(0:NRMAX) = SNBPDi(0:NRMAX) * 1.D-20 &
               &          * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))
          PNBPi0 = Eb * rKeV * 1.D20
          ! "or" = infiles(2)%totP / (2.D0 * Pi * RR * (2.D0 * Pi * INTG_F(SNBPDi)))

          ! (4) Birth Passing (SNBTGi for passing beam ions)
          i = 4
          call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBTGi,ideriv)
          ! Calibration by using total power of trapped ions
          SL = 2.D0 * Pi * INTG_F(SNBTGi)
          if(SL /= 0.d0) SNBTGi(0:NRMAX) = SNBTGi(0:NRMAX) * 1.D-20 &
               &          * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))
          PNBT10 = Eb * rKeV * 1.D20
          ! "or" = infiles(i)%totP / (2.D0 * Pi * RR * (2.D0 * Pi * INTG_F(SNBTGi)))

          ! (2) Birth TOTAL (SNBb for beam ions)
          SNBb(0:NRMAX) = SNB(0:NRMAX)

       ! *** Orbit effect ***
       else
          ! (7) Orbit Trapped (SNBPDi for trapped beam ions)
          i = 7
          call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBPDi,ideriv)
          ! Calibration by using total power of trapped ions
          SL = 2.D0 * Pi * INTG_F(SNBPDi)
          if(SL /= 0.d0) SNBPDi(0:NRMAX) = SNBPDi(0:NRMAX) * 1.D-20 &
               &          * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))
          PNBPi0 = Eb * rKeV * 1.D20
          ! "or" = infiles(i)%totP / (2.D0 * Pi * RR * (2.D0 * Pi * INTG_F(SNBPDi)))

       ! *** Orbit effect for banana ions only ***
          if(MDLNBD == 1) then
             ! (4) Birth Passing (SNBTGi for passing beam ions)
             i = 4
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBTGi,ideriv)
             ! Calibration by using total power of trapped ions
             SL = 2.D0 * Pi * INTG_F(SNBTGi)
             if(SL /= 0.d0) SNBTGi(0:NRMAX) = SNBTGi(0:NRMAX) * 1.D-20 &
                  &          * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))
             PNBT10 = Eb * rKeV * 1.D20
             ! "or" = infiles(i)%totP / (2.D0 * Pi * RR * (2.D0 * Pi * INTG_F(SNBTGi)))

             ! (4) Birth Passing + (7) Orbit Trapped (SNBb for beam ions)
             SNBb(0:NRMAX) = SNBTGi(0:NRMAX) + SNBPDi(0:NRMAX)

       ! *** Orbit effect for all ions ***
          else if(MDLNBD == 2) then
             ! (6) Orbit TOTAL (SNBb for beam ions)
             i = 6
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBb,ideriv)
             ! Calibration by using total power of trapped ions
             SL = 2.D0 * Pi * INTG_F(SNBb)
             if(SL /= 0.d0) SNBb(0:NRMAX) = SNBb(0:NRMAX) * 1.D-20 &
                  &        * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))

             ! (8) Orbit Passing (SNBTGi for passing beam ions)
             i = 8
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBTGi,ideriv)
             ! Calibration by using total power of trapped ions
             SL = 2.D0 * Pi * INTG_F(SNBTGi)
             if(SL /= 0.d0) SNBTGi(0:NRMAX) = SNBTGi(0:NRMAX) * 1.D-20 &
                  &          * (infiles(i)%totP / (Eb * rKeV * (2.D0 * Pi * RR * SL)))
             PNBT10 = Eb * rKeV * 1.D20
             ! "or" = infiles(i)%totP / (2.D0 * Pi * RR * (2.D0 * Pi * INTG_F(SNBTGi)))
          end if
       end if
       ! Collisional torque injection part
       MNB(0:NRMAX)  = PNBCD * SNBTGi(0:NRMAX) * PNBMPD

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
       SL = 4.D0 * Pi**2 * RR * INTG_F(Tqt_tmp)
       if(SL /= 0.d0) Tqt_tmp(0:NRMAX) = Tqt_tmp(0:NRMAX) * (infiles(i)%totS / SL)

    ! *** Arbitrary input *********************************************
    else if(iflag_file == 3) then
       do i = 1, n_infiles
          if(infiles(i)%name == datatype(1)) then ! Perp NB
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBP)
             SL   = 2.D0 * PI * INTG_F(SNBP)
             PNBP0  = PNBHP * 1.D6 / (2.D0 * Pi * RR * SL)
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
                SNBPDi_INTG = INTG_F(SNBPDi)
                call shift_prof(SNBPDi, 'TRAP',SIGN(1.D0,PNBCD))
                ! calibration of Perp NB amplitude of ions
                PNBPi0 = PNBP0 * (SNBPDi_INTG / INTG_F(SNBPDi))
!!                call shift_prof(SNBTi2,'PASS',SIGN(1.D0,PNBCD))
             end if
          else if(infiles(i)%name == datatype(2)) then ! Tang NB 1
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,9,SNBT1)
             SLT1 = 2.D0 * PI * INTG_F(SNBT1)
             PNBT10 = PNBHT1 * 1.D6 / (2.D0 * Pi * RR * SLT1)
          else if(infiles(i)%name == datatype(3)) then ! Tang NB 2
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,9,SNBT2)
             SLT2 = 2.D0 * PI * INTG_F(SNBT2)
             PNBT20 = PNBHT2 * 1.D6 / (2.D0 * Pi * RR * SLT2)
          else if(infiles(i)%name == datatype(4)) then ! RF
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,0,SRFe)
             SL   = 2.D0 * PI * INTG_F(SRFe)
             PRFe0 = PRFHe * 1.D6 / (2.D0 * Pi * RR * SL)
             SL   = 2.D0 * PI * INTG_F(SRFe)
             PRFi0 = PRFHi * 1.D6 / (2.D0 * Pi * RR * SL)
          end if
       end do
    end if

    !   NBI total input power (MW)
    PNBH = PNBHP + PNBHT1 + PNBHT2 + PNBHex

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

    !   Additional torque input to LQi3 without net torque

    if(Tqp0 /= 0.d0) then
       allocate(Tqp_tmp(0:NRMAX))
       do nr = 0, nrmax
          Tqp(nr) = Tqp0 * fgaussian(Rho(nr),0.6d0,0.02d0)
       end do
       SL1 = 2.D0 * PI * INTG_F(Tqp)
       do nr = 0, nrmax
          Tqp_tmp(nr) = Tqp0 * fgaussian(Rho(nr),0.7d0,0.02d0)
       end do
       SL2 = 2.D0 * PI * INTG_F(Tqp_tmp)
       do nr = 0, nrmax
          Tqp(nr) = Tqp(nr) - (SL1 / SL2) * Tqp_tmp(nr)
       end do
!!for_check       SL = 2.D0 * PI * INTG_F(Tqp)
       deallocate(Tqp_tmp)
    else
       Tqp(0:NRMAX) = 0.d0
    end if

    !   Additional torque input to LQi4

    CALL deposition_profile(Tqt,SL,0.d0,0.d0,'Additional')
    Tqt0L = Tqt0 / (2.D0*Pi*RR*2.D0*Pi*SL) ! Tqt0 [N m], Tqt0L [N/m^2]
    Tqt(0:NRMAX) = Tqt0L * Tqt(0:NRMAX)
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

!!!       Vcr = (3.D0 * SQRT(PI / 2.D0) * PNiV(NR) * PZ**2 / PNeV(NR) * AME / AMI &
!!!            &   * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)**(1.D0/3.D0)
!       Ecr = (9.D0 * PI / 16.D0 * AMI / AME)**(1.D0/3.D0) * AMB / AMI * PTeV(NR) ! in keV
       Ecr = 14.8D0 * (PA / PA**(2.D0/3.D0)) * PTeV(NR) ! in keV
       Vti = SQRT(2.D0 * ABS(PTiV(NR)) * rKeV / AMI)
       PNBcol_i(NR) = rate_to_ion(Eb/Ecr)
       PNBcol_e(NR) = 1.d0 - PNBcol_i(NR)
       IF(PNBH == 0.D0 .AND. PNbV(NR) < 1.D-8) THEN
          rNube(NR) = 0.D0
          rNubi(NR) = 0.D0
          rNuB (NR) = 0.D0
       ELSE
          rNube(NR) = PNeV(NR) * 1.D20 * PZ**2 * AEE**4 &
               &     * coulog(1.d0,PNeV(NR),PTeV(NR),PTiV(NR),AEP,1.d0,PA,PZ,Tb) &
               &     / (3.D0 * PI * SQRT(2.D0 * PI) * EPS0**2 * AMB * AME &
               &             * (ABS(PTeV(NR)) * rKeV / AME)**1.5D0)
          rNubi(NR) = PNiV(NR) * 1.D20 * PZ**2 * PZ**2 * AEE**4 &
               &     * coulog(1.d0,PNeV(NR),PTeV(NR),PTiV(NR),PA,PZ,PA,PZ,Tb) &
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
          rNuD(NR) = PNiV(NR) *1.D20 * PZ**2 * PZ**2 * AEE**4 &
               &   * coulog(1.d0,PNeV(NR),PTeV(NR),PTiV(NR),PA,PZ,PA,PZ,Tb) &
               &   / (2.D0 * PI * EPS0**2 * AMB**2 * Vb**3) &
               &   * (  SQRT(1.D0 - EXP(- 4.D0 * xl**2 / PI)) &
               &      - 2.D0 * xl / (4.D0 * xl**3 + 3.D0 * SQRT(PI)))

          !     *** Beam ion-electron slowing-down time
          !     (Tokamaks 3rd, A_D from below (2.13.2) + (2.14.1), assuming AMB >> AME
          !      and Vb/(SQRT(2)*Vte) << 1)
          !     (originally from the book by Spitzer (1962))
          !     Note that this is a half value of rNube(NR).
          rNubes = PNeV(NR) * 1.D20 * PZ**2 * AEE**4 &
               & * coulog(1.d0,PNeV(NR),PTeV(NR),PTiV(NR),PA,PZ,PA,PZ) &
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
          MNB(NR)   = PNBCD * SNBTGi(NR) * PNBMPD
       end if

       PRFe(NR)  = PRFe0 * SRFe(NR)
       PRFi(NR)  = PRFi0 * SRFi(NR)

       !     *** NBI power deposition for graphics ***

       PNBe(NR) = Eb * SNB(NR) * PNBcol_e(NR) * (1.D20 * rKeV)
       PNBi(NR) = Eb * SNB(NR) * PNBcol_i(NR) * (1.D20 * rKeV) &
            &   + AMb * Vb * MNB(NR) * (BthV(NR)*UithV(NR)+BphV(NR)*UiphV(NR))/BBL * 1.D20

    end do

  end subroutine txauxs

!***************************************************************
!
!   Heating deposition profile
!     Input : R0    : Deposition center (m)
!             RW    : Deposition width (m)
!             CHR   : Trapped NB, passing NB, RF or Additional torque
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
!!    real(8) :: AITKEN2P

    if(CHR == 'Additional') then
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

    ! Modify S(0) and S(1) so that S' becomes zero at the axis
    call sctr(R(1),R(2),R(3),R(4),S(2),S(3),S(4),S(0),S(1))

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
