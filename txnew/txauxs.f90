!!$!***********************************************************************
!!$!
!!$!   Rate coefficients with Chebyshev polynomials fit
!!$!     valid for background Ti = 1eV => 2*10^4eV
!!$!               beam energy   = 10keV/amu => 500keV/amu
!!$! 
!!$!     (C.F. Barnett, "ATOMIC DATA FOR FUSION VOLUME 1", ORNL-6086 (1990))
!!$!
!!$!     Inputs (real*8): TieV    : Hydrogenic ion temperature [eV]
!!$!                      EbeVamu : Beam energy [eV]
!!$!     Output (real*8): SiVcxb  : Charge exchange maxwellian rate coefficient with beam [m^3/s]
!!$!                      SiVizb  : Ionization maxwellian rate coefficient with beam [m^3/s]
!!$!
!!$!***********************************************************************
!!$
!!$module mod_chebyshev
!!$  implicit none
!!$  integer(4), parameter :: iemax = 1001, jmax = 7
!!$  real(8), dimension(jmax) :: EBdata
!!$  real(8), dimension(:,:,:,:), allocatable :: ucxb, uizb
!!$  real(8), dimension(:)      , allocatable :: Earray
!!$
!!$contains
!!$
!!$!***********************************************************************
!!$!
!!$!   Spline Table of rate coefficients for Chebyshev polynomials fit
!!$!
!!$!     Common (real*8): Earray : Hydrogenic ion temperature range from Emin to Emax
!!$!                      Ebdata : Discrete beam energy
!!$!                      ucxb   : 2D spline coefficients for charge exchange
!!$!                      uizb   : 2D spline coefficients for ionization by protons
!!$!     Internal (real*8): acxb : Table of Maxwellian rate coefficient in A-27
!!$!                                                   ^^^^^^^^^^^^^^^^
!!$!                        aizb : Table of Maxwellian rate coefficient in D-10
!!$!                                                   ^^^^^^^^^^^^^^^^
!!$!***********************************************************************
!!$
!!$
!!$  subroutine spline_table_beam_rate_coef
!!$    implicit none
!!$    real(8), dimension(0:6,jmax) :: acx, aiz
!!$    real(8), dimension(:,:), allocatable :: ratedata, dx, dy, dxy
!!$    integer(4) :: i, j, ierr
!!$    real(8) :: Emax, Emin, dElin, Elin
!!$
!!$    !  H + H^+ -> H + H^+
!!$    !  Charge exchange (A-27)
!!$    !  Beam - Maxwellian rate coefficients
!!$
!!$    acx(0:6,1:jmax) = &
!!$         & reshape((/ -3.23588d+1, -2.28428d-1, -1.47539d-1, -7.73418d-2  &
!!$                & ,   -3.29002d-2, -1.09141d-2, -2.93158d-3  & ! 10keV/amu
!!$                & ,   -3.27829d+1, -2.26775d-1, -1.38750d-1, -7.04432d-2  &
!!$                & ,   -3.03096d-2, -1.00013d-2, -1.47228d-3  & ! 20keV/amu
!!$                & ,   -3.39547d+1, -1.36550d-1, -7.15786d-2, -3.65682d-2  &
!!$                & ,   -2.10115d-2, -1.11355d-2, -4.14062d-3  & ! 40keV/amu
!!$                & ,   -3.60736d+1,  1.09084d-1,  7.35296d-2,  2.31264d-2  &
!!$                & ,   -7.74992d-3, -1.38165d-2, -8.49906d-3  & ! 70keV/amu
!!$                & ,   -3.81082d+1,  3.09107d-1,  1.95996d-1,  8.14408d-2  &
!!$                & ,    1.35787d-2, -9.21346d-3, -1.02081d-2  & ! 100keV/amu
!!$                & ,   -4.29608d+1,  3.96726d-1,  2.99000d-1,  1.81207d-1  &
!!$                & ,    8.74182d-2,  3.26701d-2,  7.33350d-3  & ! 200keV/amu
!!$                & ,   -5.17564d+1,  3.37283d-1,  2.43298d-1,  1.48098d-1  &
!!$                & ,    7.72749d-2,  3.54105d-2,  1.56884d-2  & ! 500keV/amu
!!$                &   /), shape(acx) )
!!$
!!$    !  H + H^+ -> H^+ + H^ + e^-
!!$    !  Ionization by protons (D-11)
!!$    !  Beam - Maxwellian rate coefficients
!!$
!!$    aiz(0:6,1:jmax) = &
!!$         & reshape((/ -3.82968d+1,  1.12152d00,  6.44681d-1,  1.00717d-1  &
!!$                & ,   -1.00166d-1, -3.95957d-2, -2.65964d-3  & ! 10keV/amu
!!$                & ,   -3.57290d+1,  3.56274d-1,  2.02092d-1,  6.20743d-2  &
!!$                & ,   -2.76894d-3, -1.15918d-2, -4.10154d-3  & ! 20keV/amu
!!$                & ,   -3.41679d+1, -1.62401d-2,  3.246363d-4, 7.54856d-3  &
!!$                & ,    8.28055d-3,  4.20183d-3,  1.24748d-3  & ! 40keV/amu
!!$                & ,   -3.37606d+1, -4.00418d-2, -2.58461d-2, -1.20904d-2  &
!!$                & ,   -3.36001d-3,  2.548592d-4, 9.663075d-4 & ! 70keV/amu
!!$                & ,   -3.37204d+1, -2.56605d-2, -1.80094d-2, -1.07132d-2  &
!!$                & ,   -5.35053d-3, -2.08862d-3, -4.645604d-4 & ! 100keV/amu
!!$                & ,   -3.38801d+1, -1.46171d-2, -9.02218d-3, -4.77874d-3  &
!!$                & ,   -2.29897d-3, -9.828762d-4,-3.658987d-4 & ! 200keV/amu
!!$                & ,   -3.44565d+1, -7.913823d-4,-5.710595d-4,-3.486254d-4 &
!!$                & ,   -1.820028d-4,-8.154704d-5,-3.512063d-5 & ! 500keV/amu
!!$                &   /), shape(aiz) )
!!$
!!$    EBdata(1:jmax) = (/ 1.d4, 2.d4, 4.d4, 7.d4, 1.d5, 2.d5, 5.d5 /)
!!$
!!$    Emin = 1.d0 ! [eV]
!!$    Emax = 2.d4 ! [eV]
!!$
!!$    allocate(Earray(iemax))
!!$    allocate(ucxb(4,4,iemax,jmax),uizb(4,4,iemax,jmax))
!!$
!!$    allocate(ratedata(iemax,jmax))
!!$    allocate(dx, dy, dxy, mold=ratedata)
!!$
!!$    ! --- Equally-spaced in logspace ---
!!$    dElin = (log10(Emax) - log10(Emin)) / (iemax - 1)
!!$
!!$    Elin      = log10(Emin)
!!$    Earray(1) = 1.d0**Elin
!!$    do i = 2, iemax
!!$       Elin      = Elin + dElin
!!$       Earray(i) = 1.d1**Elin
!!$    end do
!!$
!!$    ! --- Charge exchange ---
!!$
!!$    do j = 1, jmax
!!$       do i = 1, iemax
!!$          ratedata(i,j) = func_ratecoef(Earray(i),Emax,Emin,acx(:,j))
!!$       end do
!!$    end do
!!$
!!$    !   Compute spline coefficients
!!$
!!$    call spl2d(Earray,EBdata,ratedata,dx,dy,dxy,ucxb,iemax,iemax,jmax,0,0,ierr)
!!$    if( ierr /= 0 ) stop 'spl2d error in spline_table_beam_rate_coef; cx.'
!!$
!!$    ! --- Ionization ---
!!$
!!$    do j = 1, jmax
!!$       do i = 1, iemax
!!$          ratedata(i,j) = func_ratecoef(Earray(i),Emax,Emin,aiz(:,j))
!!$       end do
!!$    end do
!!$
!!$    !   Compute spline coefficients
!!$
!!$    call spl2d(Earray,EBdata,ratedata,dx,dy,dxy,uizb,iemax,iemax,jmax,0,0,ierr)
!!$    if( ierr /= 0 ) stop 'spl2d error in spline_table_beam_rate_coef; iz.'
!!$
!!$    deallocate(ratedata,dx,dy,dxy)
!!$
!!$  end subroutine spline_table_beam_rate_coef
!!$
!!$  subroutine deallocate_spline_table_beam_rate_coef
!!$    implicit none
!!$
!!$    if(allocated(Earray)) deallocate(Earray,ucxb,uizb)
!!$
!!$  end subroutine deallocate_spline_table_beam_rate_coef
!!$
!!$!***********************************************************************
!!$!
!!$!   Charge exchange rate coefficients for beam H + Maxwellian H^+ -> H + H^+
!!$!
!!$!     valid for background Ti = 1eV => 2*10^4eV
!!$!               beam energy   = 10keV/amu => 500keV/amu
!!$! 
!!$!     Inputs (real*8): TieV    : Hydrogenic ion temperature [eV]
!!$!                    : EbeVamu : Beam energy [eV]
!!$!     Output (real*8): SiVcxb  : Charge exchange maxwellian rate coefficient with beam [m^3/s]
!!$!
!!$!***********************************************************************
!!$
!!$  function SiVcxB(TieV, EbeVamu)
!!$    implicit none
!!$    ! input
!!$    real(8), intent(in) :: TieV, EbeVamu ! [eV]
!!$    ! output
!!$    real(8) :: SiVcxB
!!$    ! local
!!$    integer(4) :: ierr
!!$
!!$    if( TieV < Earray(1) .or. TieV > Earray(iemax) ) then
!!$       write(6,*) 'Function SiVcxB: Out of energy range. Ti [eV]=', TieV
!!$    else if( EbeVamu < EBdata(1) .or. EbeVamu > EBdata(jmax) ) then
!!$       write(6,*) 'Function SiVcxB: Out of energy range. Eb [eV/amu]=', EbeVamu
!!$    end if
!!$
!!$    call spl2df(TieV,EbeVamu,SiVcxB,Earray,EBdata,ucxb,iemax,iemax,jmax,ierr)
!!$
!!$  end function SiVcxB
!!$
!!$!***********************************************************************
!!$!
!!$!   Ionization rate coefficients for beam H + Maxwellian H^+ -> H^+ + H^+ + e^-
!!$!
!!$!     valid for background Ti = 1eV => 2*10^4eV
!!$!               beam energy   = 10keV/amu => 500keV/amu
!!$! 
!!$!     Inputs (real*8): TieV    : Hydrogenic ion temperature [eV]
!!$!                    : EbeVamu : Beam energy [eV]
!!$!     Output (real*8): SiVizb  : Ionization maxwellian rate coefficient with beam [m^3/s]
!!$!
!!$!***********************************************************************
!!$
!!$  function SiVizB(TieV, EbeVamu)
!!$    implicit none
!!$    ! input
!!$    real(8), intent(in) :: TieV, EbeVamu ! [eV]
!!$    ! output
!!$    real(8) :: SiVizB
!!$    ! local
!!$    integer(4) :: ierr
!!$
!!$    if( TieV < Earray(1) .or. TieV > Earray(iemax) ) then
!!$       write(6,*) 'Function SiVizB: Out of energy range. Ti [eV]=', TieV
!!$    else if( EbeVamu < EBdata(1) .or. EbeVamu > EBdata(jmax) ) then
!!$       write(6,*) 'Function SiVizB: Out of energy range. Eb [eV/amu]=', EbeVamu
!!$    end if
!!$
!!$    call spl2df(TieV,EbeVamu,SiVizB,Earray,EBdata,uizb,iemax,iemax,jmax,ierr)
!!$
!!$  end function SiVizB
!!$
!!$! ---
!!$
!!$  function func_ratecoef(E,Emax,Emin,a)
!!$    implicit none
!!$    ! input
!!$    real(8), intent(in) :: E, Emax, Emin
!!$    real(8), dimension(0:), intent(in) :: a
!!$    ! output
!!$    real(8) :: func_ratecoef
!!$    ! local
!!$    real(8) :: lnsigE, x
!!$    real(8) :: cm3tom3 = 1.d-6
!!$
!!$    x = ( ( log(E) - log(Emin) ) - ( log(Emax) - log(E) ) ) / ( log(Emax) - log(Emin) )
!!$
!!$    lnsigE = 0.5d0 * a(0) + chebypoly(x,a)
!!$
!!$    func_ratecoef = exp(lnsigE) * cm3tom3
!!$
!!$  end function func_ratecoef
!!$
!!$! ---
!!$
!!$  function chebypoly(x,a) result(AT)
!!$    implicit none
!!$    ! input
!!$    real(8), intent(in) :: x
!!$    real(8), dimension(0:), intent(in) :: a
!!$    ! output
!!$    real(8) :: AT
!!$    ! local
!!$    integer(4) :: i, j, imax
!!$    real(8) :: T(0:8)
!!$
!!$    imax = size(a) - 2
!!$
!!$    T(0) = 1.d0
!!$    T(1) = x
!!$    AT = a(1) * T(1)
!!$
!!$    do i = 1, imax
!!$       j = i + 1
!!$       ! Chebyshev polynomials
!!$       T(j) = 2.d0 * x * T(j-1) - T(j-2)
!!$       ! Sum A_jT_j
!!$       AT = AT + a(j) * T(j)
!!$    end do
!!$
!!$  end function chebypoly
!!$
!!$end module mod_chebyshev
!!$
! ----------------------------------------------------------------------

!***************************************************************
!
!   Auxiliary heating system
!
!***************************************************************

module aux_system
  implicit none
  private
  real(8), save :: Vbabsmax!, Scxb
  public :: txauxs, Vbabsmax!, Scxb

contains

  subroutine txauxs

    ! SNB*(NR) : The birth rate of beam ions per unit volume
    ! MNB(NR)  : In essence MNB is similar to SNB, but MNB distinguishes itself
    !            from SNB in that MNB takes into account the tangential NB direction,
    !            the pitch between the NB chord and the field line and the contribution
    !            of the perpendicular NBI to the momentum input.

    use tx_commons
    use tx_interface,      only : inexpolate, fgaussian, coulog
    use tx_core_module,    only : intg_vol
    use mod_cross_section, only : SiVcxb, SiVizb

    integer(4) :: NR, i, ideriv = 1, nrbound!, izmodel = 1
    real(8) :: SL, SLT1, SLT2, PNBP0, PNBT10, PNBT20, PNBex0, SNBPDi_INTG, &
         &     PNBPi0, PNBTi10, PNBTi20, PRFe0, PRFi0, SL1, SL2, &
         &     Vti, rNubes, BBL, PALFL, Ecr, &
         &     zEbeV, zEbkev, xl, Tqt0L, TieV, rateizb!, Sion, logEbeV
    real(8), parameter :: PAHe = 4.d0, & ! Atomic mass number of He
         &                Enf  = 3.5d3   ! in keV, equal to 3.5 MeV
    real(8), dimension(0:NRMAX) :: SNBP, SNBT1, SNBT2, SNBTi1, SNBTi2, SRFe, SRFi
    ! Mainly for derivatives
    real(8), dimension(:), allocatable :: zTqt, zTqp

    !     *** Constants ***

    !     Averaged beam injection energy

    Eb =  (esps(1) + 0.5d0 * esps(2) + esps(3) / 3.d0) * Ebmax

    !     NBI beam speed

    Vb =  sqrt(2.d0 * Eb * rKilo / (amb * amqp) )
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
    if(PNBHP /= 0.d0) then
       call deposition_profile(SNBP,SL,RNBP0,RNBP,'NB')
       PNBP0 = PNBHP * 1.d6 / SL
       PNBPi0 = PNBP0
       SNBPDi(0:NRMAX) = SNBP(0:NRMAX)

       if(MDLNBD /= 0) then ! Orbit effect
!!$          ! For ions
!!$          if(PNBMPD == 0.d0) then ! Exact perpendicular NBI
!!$             call deposition_profile(SNBPDi,SL,RNBP0,RNBP,'NB_TRAP',0.d0)
!!$          else ! Near perpendicular NBI
!!$             call deposition_profile(SNBPDi,SL,RNBP0,RNBP,'NB_TRAP',sign(1.d0,PNBMPD))
!!$          end if
          ! For ions
          if(PNBCD == 0.d0) then ! Exact perpendicular NBI
             call deposition_profile(SNBPDi,SL,RNBP0,RNBP,'NB_TRAP',0.d0)
          else ! Near perpendicular NBI
             if(PNBMPD == 0.d0) then ! No collisional torque
                call deposition_profile(SNBPDi,SL,RNBP0,RNBP,'NB_TRAP',PNBCD)
             else ! W/ collisional torque
                call deposition_profile(SNBPDi,SL,RNBP0,RNBP,'NB_TRAP',sign(1.d0,PNBMPD))
             end if
          end if
          PNBPi0 = PNBHP * 1.d6 / SL
       end if
    else
       PNBP0  = 0.d0
       PNBPi0 = 0.d0
       SNBP(0:NRMAX)   = 0.d0
       SNBPDi(0:NRMAX) = 0.d0
    end if

    !  *** Tangential
    if(PNBHT1 /= 0.d0) then
       call deposition_profile(SNBT1,SLT1,RNBT10,RNBT1,'NB')
       PNBT10 = PNBHT1 * 1.d6 / SLT1

       if(MDLNBD > 1) then ! Orbit effect
          ! For ions
          call deposition_profile(SNBTi1,SLT1,RNBT10,RNBT1,'NB_PASS',sign(1.d0,PNBCD))
          PNBTi10 = PNBHT1 * 1.d6 / SLT1
       end if
    else
       PNBT10 = 0.d0
       SNBT1(0:NRMAX)  = 0.d0
       SNBTi1(0:NRMAX) = 0.d0
    end if

    if(PNBHT2 /= 0.d0) then
       call deposition_profile(SNBT2,SLT2,RNBT20,RNBT2,'NB')
       PNBT20 = PNBHT2 * 1.d6 / SLT2

       if(MDLNBD > 1) then
          ! For ions
          call deposition_profile(SNBTi2,SLT2,RNBT20,RNBT2,'NB_PASS',sign(1.d0,PNBCD))
          PNBTi20 = PNBHT2 * 1.d6 / SLT2
       end if
    else
       PNBT20 = 0.d0
       SNBT2(0:NRMAX)  = 0.d0
       SNBTi2(0:NRMAX) = 0.d0
    end if

    !  For RF heating
    if(PRFHe /= 0.d0) then
       call deposition_profile(SRFe,SL,RRFe0,RRFew,'RFe')
       PRFe0 = PRFHe * 1.d6 / SL
    else
       PRFe0 = 0.d0
       SRFe(0:NRMAX) = 0.d0
    end if

    if(PRFHi /= 0.d0) then
       call deposition_profile(SRFi,SL,RRFi0,RRFiw,'RFi')
       PRFi0 = PRFHi * 1.d6 / SL
    else
       PRFi0 = 0.d0
       SRFi(0:NRMAX) = 0.d0
    end if

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
       PNBex0 = Eb * rKeV * 1.d20
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
          PNBPi0 = Eb * rKeV * 1.d20
          ! "or" = infiles(2)%totP / intg_vol(SNBPDi)

          ! (4) Birth Passing (SNBTGi for passing beam ions)
          i = 4
          call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBTGi,ideriv)
          ! Calibration by using total power of trapped ions
          SL = intg_vol(SNBTGi)
          if(SL /= 0.d0) SNBTGi(0:NRMAX) = SNBTGi(0:NRMAX) * 1.D-20 &
               &          * (infiles(i)%totP / (Eb * rKeV * SL))
          PNBT10 = Eb * rKeV * 1.d20
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
          PNBPi0 = Eb * rKeV * 1.d20
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
             PNBT10 = Eb * rKeV * 1.d20
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
             PNBT10 = Eb * rKeV * 1.d20
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

       allocate(zTqt, mold=rho)
       ! Total torque input
       i = 1
       call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,zTqt,ideriv,idx=0)
       SL = intg_vol(zTqt)
       if(SL /= 0.d0) zTqt(0:NRMAX) = zTqt(0:NRMAX) * (infiles(i)%totS / SL)

    ! *** Arbitrary input *********************************************
    else if(iflag_file == 3) then
       do i = 1, n_infiles
          if(infiles(i)%name == datatype(1)) then ! Perp NB
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,5,SNBP)
             SL   = intg_vol(SNBP)
             PNBP0  = PNBHP * 1.d6 / SL
             ! Banana orbit effect
             if(MDLNBD /= 0) then! .and. PNBMPD /= 0.d0) then
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
                call shift_prof(SNBPDi, 'TRAP',sign(1.d0,PNBCD))
                ! calibration of Perp NB amplitude of ions
                PNBPi0 = PNBP0 * (SNBPDi_INTG / intg_vol(SNBPDi))
!!                call shift_prof(SNBTi2,'PASS',sign(1.d0,PNBCD))
             end if
          else if(infiles(i)%name == datatype(2)) then ! Tang NB 1
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,9,SNBT1)
             SLT1 = intg_vol(SNBT1)
             PNBT10 = PNBHT1 * 1.d6 / SLT1
          else if(infiles(i)%name == datatype(3)) then ! Tang NB 2
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,9,SNBT2)
             SLT2 = intg_vol(SNBT2)
             PNBT20 = PNBHT2 * 1.d6 / SLT2
          else if(infiles(i)%name == datatype(4)) then ! RF
             call inexpolate(infiles(i)%nol,infiles(i)%r,infiles(i)%data,NRMAX,RHO,0,SRFe)
             SL   = intg_vol(SRFe)
             PRFe0 = PRFHe * 1.d6 / SL
             SL   = intg_vol(SRFe)
             PRFi0 = PRFHi * 1.d6 / SL
          end if
       end do
    end if

    !   NBI total input power (MW)
    PNBH = PNBHP + PNBHT1 + PNBHT2 + PNBHex

    !   Ratio of CX deposition rate to IZ deposition rate

    RatCX(:) = 0.d0

    zEbkeV  = Eb / amb
    zEbeV   = zEbkeV * 1.d3
!!$    logEbeV = log10(zEbeV)
!!$
!!$    !     CX cross-section
!!$    !     (Riviere, NF 11 (1971) 363)
!!$    Scxb = 6.937D-19 * (1.d0 - 0.155d0 * logEbeV)**2 / (1.d0 + 1.112D-15 * zEbeV**3.3d0)
!!$
!!$    IF(PNBH /= 0.d0) THEN
!!$       !     Ionization cross-section
!!$       if( izmodel == 0 ) then
!!$          !     (Riviere, NF 11 (1971) 363)
!!$          IF(zEbeV > 150.d0) THEN
!!$             Sion = 3.6D-16 / zEbeV * (- 0.7783d0 + logEbeV)
!!$          ELSE
!!$             Sion = 10.d0**(-0.8712d0 * logEbeV**2 + 8.156d0 * logEbeV - 38.833d0)
!!$          END IF
!!$       else
!!$          !     (R.L. Freeman and E.M. Jones, Culham Laboratory Report, CLM-R-137 (May 1974))
!!$          Sion = SigIzPro(zEbkeV)
!!$       end if
!!$
!!$       RatCX(0:NRMAX) = Scxb / (Scxb + Sion)
!!$    END IF
    if( PNBH /= 0.d0 ) then
       do nr = 0, nrmax
          TieV        = Var(NR,2)%T * rKilo
          rateizb     = SiVizB(TieV, zEbeV)
          ratecxb(NR) = SiVcxB(TieV, zEbeV)
!          RatCX(NR)   = 0.d0 ! Ionization only
!          RatCX(NR)   = 1.d0 ! Charge Exchange only
          RatCX(NR)   = ratecxb(NR) / (ratecxb(NR) + rateizb)
       end do
    end if

    !   Alpha heating

    do nr = 0, nrmax
       PALFL = FSNF * bosch_fusion(Var(NR,2)%T,0.5d0*Var(NR,2)%n,0.5d0*Var(NR,2)%n)
       Ecr = 14.8d0 * (PAHe / amb**(2.d0/3.d0)) * Var(NR,1)%T ! in keV
       PALFi(NR) = PALFL * rate_to_ion(Enf/Ecr)
       PALFe(NR) = PALFL - PALFi(NR)
    end do

    !   Additional torque input to LQi3 without net torque

    if(Tqp0 /= 0.d0) then
       allocate(zTqp, mold=rho)
       do nr = 0, nrmax
          Tqp(nr) = Tqp0 * fgaussian(Rho(nr),0.6d0,0.02d0)
       end do
       SL1 = intg_vol(Tqp)
       do nr = 0, nrmax
          zTqp(nr) = Tqp0 * fgaussian(Rho(nr),0.7d0,0.02d0)
       end do
       SL2 = intg_vol(zTqp)
       do nr = 0, nrmax
          Tqp(nr) = Tqp(nr) - (SL1 / SL2) * zTqp(nr)
       end do
       deallocate(zTqp)
    else
       Tqp(0:NRMAX) = 0.d0
    end if

    !   Additional torque input to LQi4

!    call deposition_profile(Tqt,SL,0.d0,0.d0,'Additional')
    call deposition_profile(Tqt,SL,RNBT10,RNBT1,'Additional')
    Tqt0L = Tqt0 / SL ! Tqt0 [N m], Tqt0L [-]
    Tqt(0:NRMAX) = Tqt0L * Tqt(0:NRMAX) ! <R.Source> [Nm/m^3]
!    Tqt(0:NRMAX) = Tqt0L * Tqt(0:NRMAX) * vlt(0:NRMAX) / vlt(NRMAX)!test1
!    Tqt(0:NRA) = Tqt0 * exp(- ((vlt(0:NRA)/vlt(NRA) - RNBT10) / RNBT1)**2) * (1.d0 - (vlt(0:NRA)/vlt(NRA))**4)!test2
!    Tqt(NRA+1:NRMAX) = 0.d0
    if(iflag_file == 2) then
       Tqt(0:NRMAX) = Tqt(0:NRMAX) + zTqt(0:NRMAX)
       deallocate(zTqt)
    end if

    ! ************** Heating part end **************

    !     *** Beam slowing down time (momentum transfer with beam) ***
    ! reference : memo (92/04/02, 92/04/21)
    !             Tokamaks 3rd pp.246 - 252

    do NR = 0, NRMAX
       BBL = sqrt(BphV(NR)**2 + BthV(NR)**2)

!!!       Vcr = (3.d0 * sqrt(Pi / 2.d0) * Var(NR,2)%n * achgb**2 / Var(NR,1)%n * amas(1) / amas(2) &
!!!            &   * (abs(Var(NR,1)%T) * rKilo / (amas(1) * amqp))**1.5d0)**(1.d0/3.d0)
!       Ecr = (9.d0 * Pi / 16.d0 * amas(2) / amas(1) )**(1.d0/3.d0) * amb / amas(2) * Var(NR,1)%T ! in keV
       Ecr = 14.8d0 * (amb / amb**(2.d0/3.d0)) * Var(NR,1)%T ! in keV
       Vti = sqrt(2.d0 * abs(Var(NR,2)%T) * rKilo / (amas(2) * amqp))
       PNBcol_i(NR) = rate_to_ion(Eb/Ecr)
       PNBcol_e(NR) = 1.d0 - PNBcol_i(NR)
       if(PNBH == 0.d0 .and. PNbV(NR) < 1.D-8) then
          rNuD (NR) = 0.d0
          rNuB (NR) = 0.d0
       else
          !     *** deflection time of beam ions against bulk ions ***
          !     (Takamura (3.26) + Tokamaks 3rd p64)
          xl = Vb / Vti
          rNuD(NR) = Var(NR,2)%n *1.d20 * achg(2)**2 * achgb**2 * AEE**4 &
               &   * coulog(Zeff(NR),Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amas(2),achg(2),amb,achgb,PTbV(NR)) &
               &   / (2.d0 * Pi * EPS0**2 * (amb*amp)**2 * Vb**3) &
               &   * (  sqrt(1.d0 - exp(- 4.d0 * xl**2 / Pi)) &
               &      - 2.d0 * xl / (4.d0 * xl**3 + 3.d0 * sqrt(Pi)))

          !     *** Beam ion-electron slowing-down time
          !     (Tokamaks 3rd, A_D from below (2.13.2) + (2.14.1), assuming AMB >> amas(1)
          !      and Vb/(sqrt(2)*Vte) << 1)
          !     (originally from the book by Spitzer (1962))
          rNubes = Var(NR,1)%n * 1.d20 * achg(1)**2 * achgb**2 * AEE**4 &
               & * coulog(Zeff(NR),Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amb,achgb,amas(1),achg(1)) &
!               & * coulog(Zeff(NR),Var(NR,1)%n,Var(NR,1)%T,Var(NR,2)%T,amb,achgb,amb,achgb) & ! original, but could be wrong
               & / (6.d0 * Pi * sqrt(2.d0 * Pi) * EPS0**2 * amb * amas(1) * amp**2 &
               & * (abs(Var(NR,1)%T) * rKilo / (amas(1) * amqp))**1.5d0)

          ! The definition given below is a "energy slowing down time" and is not
          ! "particle slowing down time". At present, we assume the former is the
          ! same as the latter.
          rNuB(NR) = rNubes * 3.d0 / log(1.d0 + (Eb / Ecr)**1.5d0)
       end if

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
             SNBTGi(NR)=(PNBTi10 * SNBTi1(NR) + PNBTi20 * SNBTi2(NR)) / (Eb * rKeV * 1.d20)
          else
             SNBTGi(NR)=(PNBT10  * SNBT1(NR)  + PNBT20  * SNBT2(NR))  / (Eb * rKeV * 1.d20)
          end if
          ! Source profile for trapped ions with banana orbit effect
          !   in temporal and graphic use
          SNBPDi(NR)= PNBPi0 * SNBPDi(NR) / (Eb * rKeV * 1.d20)
          ! Birth profiles for heating power
          SNB(NR)   = PNB(NR) / (Eb * rKeV * 1.d20)
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

       PNBe(NR) = Eb * SNB(NR) * PNBcol_e(NR) * (1.d20 * rKeV)
       PNBi(NR) = Eb * SNB(NR) * PNBcol_i(NR) * (1.d20 * rKeV) &
            &   + amb*amp * MNB(NR) * Var(NR,2)%BUpar / bbt(NR) * BUbparV(NR) * 1.d20

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
    use tx_commons, only : NRMAX, NRA, FSRP, RA, Pi, RR, &
         &                 amb, amqp, Vb, achgb, BthV, Q, BphV, rho, epst
    use tx_core_module, only : intg_vol
    use libitp, only : sctr
    real(8), intent(in)  :: R0, RW
    real(8), intent(in), optional :: PNBCD
    character(len=*), intent(in) :: CHR
    real(8), intent(out), dimension(0:NRMAX) :: S
    real(8), intent(out) :: SINT
    integer(4) :: nr
    real(8) :: Rshift, Rpotato, rhop

    if(CHR == 'Additional') then
!       S(0:NRA) = 1.d0 - rho(0:NRA)**2
       S(0:NRA) = exp(- ((rho(0:NRA) - R0) / RW)**2) * (1.d0 - rho(0:NRA)**4)
       S(NRA+1:NRMAX) = 0.d0
    else if(CHR == 'RFe' .or. CHR == 'RFi') then
       S(0:NRMAX) = exp(- ((rho(0:NRMAX) - R0) / RW)**2) * (1.d0 - (rho(0:NRMAX)/rho(nrmax))**4)
    else if(CHR == 'NB') then
       if(abs(FSRP) > 0.d0) then
          S(0:NRMAX) = exp(- ((rho(0:NRMAX) - R0) / RW)**2) * (1.d0 - (rho(0:NRMAX)/rho(nrmax))**4)
       else
          S(0:NRA) = exp(- ((rho(0:NRA) - R0) / RW)**2) * (1.d0 - rho(0:NRA)**4)
          S(NRA+1:NRMAX) = 0.d0
       end if
    else
       if(present(PNBCD) .eqv. .false.) stop 'deposition_profile: input error!'
       do nr = 0, nrmax
          if(nr /= 0) rhop = amb * Vb * amqp / (achgb * BthV(NR)) ! poloidal Larmor radius
          if(CHR == 'NB_TRAP') then ! trapped particle
             ! potato width
             Rpotato = (Q(NR)**2*(amb * Vb * amqp / (achgb * BphV(NR)))**2*RR)**(1.d0/3.d0)
             if(nr == 0) then
                Rshift = PNBCD * Rpotato ! potato particle
             else
                Rshift = PNBCD * min(sqrt(epst(nr)) * rhop, Rpotato) ! potato or banana particle
             end if
          else
             if (nr == 0) then ! passing particle
                Rshift = PNBCD * (amb * Vb * Q(NR) * amqp / (achgb * BphV(NR)))
             else
                Rshift = PNBCD * (     epst(nr)  * rhop)
             end if
          end if

          if(abs(FSRP) > 0.d0) then
             S(NR) = exp(- (((rho(nr) + Rshift/RA) - R0) / RW)**2) * (1.d0 - (rho(nr)/rho(nrmax))**4)
          else
             if(nr <= nra) then
                S(NR) = exp(- (((rho(nr) + Rshift/RA) - R0) / RW)**2) * (1.d0 - rho(NR)**4)
             else
                S(NR) = 0.d0
             end if
          end if
       end do
!!$       if(CHR == 'NB_TRAP') then
!!$          if(PNBCD > 0.d0) then
!!$             S(0) = AITKEN2P(rho(0),S(1),S(2),S(3),rho(1),rho(2),rho(3))
!!$          else if(PNBCD < 0.d0) then
!!$             S(0) = 0.d0
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
    use tx_commons, only : NRMAX, RR, amb, amqp, Vb, achgb, BthV, Q, BphV, epst, rpt
    use libitp, only: aitken
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
          Rpotato = (Q(NR)**2*(amb * Vb * amqp / (achgb * BphV(NR)))**2*RR)**(1.d0/3.d0)
          if(nr == 0) then
             Rshift = direct * Rpotato ! potato particle
          else
             Rshift = direct * min(sqrt(epst(nr)) * rhop, Rpotato) ! potato or banana particle
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

       r_shift(nr) = rpt(nr) - Rshift
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
          call aitken(rpt(nr),fl1(nr),r_alloc,f_alloc,2,nrl+1)
       end do
       deallocate(r_alloc,f_alloc)
    end if

    ! Interpolate
    do nr = 0, nrmax
       call aitken(rpt(nr),fl2(nr),r_shift,f,2,nrmax+1)
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
    use tx_commons, only : Pi
    real(8), intent(in) :: x

    if (x == 0.d0) then
       f = 1.d0
    else
       f = 1.d0 / x * (  1.d0 / 3.d0 * log((1.d0 - sqrt(x) + x) / (1.d0 + sqrt(x))**2) &
            &          + 2.d0 / sqrt(3.d0) * (atan((2.d0 * sqrt(x) - 1.d0) / sqrt(3.d0)) &
            &          + Pi / 6.d0))
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

!***************************************************************
!
!   Ionization of atomic hydrogen by protons
!     valid for 100eV => 5*10^6eV
! 
!     (R.L. Freeman and E.M. Jones, Culham Laboratory Report, CLM-R-137 (May 1974))
!
!     Inputs (real*8): Ekev     : Energy [keV]
!     Output (real*8): SigIzPro : Ionization cross-section by protons [m^2]
!
!     Internal (real*8): a : Table of cross-section in TABLE 2
!                                     ^^^^^^^^^^^^^
!***************************************************************

  real(8) function SigIzPro(Ekev)

    real(8), intent(in) :: Ekev
    real(8) :: x, Ekev_temp
    real(8), dimension(0:6) :: a
    data a /-0.4203309d02, 0.3557321d01, -0.1045134d01, 0.3139238d0, &
         &  -0.7454475d-1, 0.8459113d-2, -0.3495444d-3/

    if(Ekev < 1.d-1) then
       write(6,'(A,ES12.4)') &
            'Function SigIzPro: Out of energy range. Ekev=', Ekev
       Ekev_temp=1.d-1
    else if(Ekev > 5.d3) then
       write(6,'(A,ES12.4)') &
            'Function SigIzPro: Out of energy range. Ekev=', Ekev
       Ekev_temp=5.d3
    else
       Ekev_temp=Ekev
    endif
    x = log(Ekev_temp)
    SigIzPro = exp(a(0)+x*(a(1)+x*(a(2)+x*(a(3)+x*(a(4)+x*(a(5)+x*a(6)))))))*1.d-4

  end function SigIzPro

end module aux_system

