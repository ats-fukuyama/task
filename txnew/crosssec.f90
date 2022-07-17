!***********************************************************************
!
!   Rate coefficients
!
!***********************************************************************

module mod_cross_section
  implicit none
  ! Carbon crosssections
  integer :: ntemp, ndens
  real(8), dimension(2) :: templim, denslim
  real(8), dimension(:),       allocatable :: ztemp, zdens
  real(8), dimension(:,:,:,:), allocatable :: urca6, urcs5, urcf5
  

  ! Beam crosssections
  integer(4), parameter :: iemax = 1001, jmax = 7
  real(8), dimension(jmax) :: EBdata, EHdata
  real(8), dimension(:,:,:,:), allocatable :: ucxb, uizb, ucxc
  real(8), dimension(:)      , allocatable :: Earray, ECarray

  !  e^- + C^q+ -> C^(q+1)+ + 2e^-
  !  Electron-impact ionization cross sections ([2] 4-20)

  real(8), dimension(0:5,jmax) :: &
       & aizc = reshape((/  1.090d-7,  2.992d-8, -4.748d-8, -9.354d-9  &
       &               ,    1.802d-8, -1.348d-9, -5.845d-9  & ! C
       &               ,    3.044d-8,  8.641d-9, -1.213d-8, -1.607d-9  &
       &               ,    4.621d-9, -1.031d-9, -1.710d-9  & ! C^+
       &               ,    8.804d-9,  2.795d-9, -3.188d-9, -6.010d-10 &
       &               ,    1.180d-9, -2.497d-10,-4.485d-10 & ! C^2+
       &               ,    2.741d-9,  1.113d-9, -7.997d-10,-4.238d-10 &
       &               ,    2.343d-10, 5.672d-11,-5.945d-11 & ! C^3+
       &               ,    4.783d-10, 2.513d-10,-7.327d-11,-9.965d-11 &
       &               ,    1.234d-11, 3.381d-11, 7.284d-12 & ! C^4+
       &               ,    1.594d-10, 8.604d-11,-2.042d-11,-3.194d-11 &
       &               ,    2.997d-12, 1.052d-11, 2.367d-12 & ! C^5+
       &                 /), shape(aizc))

contains

!***********************************************************************
!
!   Spline Table of Maxwellian rate coefficients of superstaged carbon
!
!     Data produced with OpenADAS via "paramsurvey" program
!
!     valid for Te = 1 eV => 100 keV
!               Ne = 1e17 m^{-3} => 1e21 m^{-3}
!     * A valid range of data can vary by modifying "Input.data" of paramsurvey.
! 
!     Common (real*8): Earray : Hydrogenic ion temperature range from Emin to Emax
!                      Ebdata : Discrete beam energy
!                      ucxb   : 2D spline coefficients for charge exchange
!                      uizb   : 2D spline coefficients for ionization by protons
!     Internal (real*8): acxb : Table of Maxwellian rate coefficient in A-27
!                                                   ^^^^^^^^^^^^^^^^
!                        aizb : Table of Maxwellian rate coefficient in D-10
!                                                   ^^^^^^^^^^^^^^^^
!***********************************************************************

  subroutine spline_table_carbon_rate_coef_adas
    use libspl2d, only : spl2d
    integer :: iocrca, ist, ierr
    integer :: i, it, in
    character(len=:), allocatable :: fcrca 
    character(len=80) :: readline
    real(8) :: zdum1, zdum2
    real(8), dimension(:,:),   allocatable :: dx, dy, dxy
    real(8), dimension(:,:,:), allocatable :: ratedata

    fcrca = 'data/Carbon-A6S5F6_20220606.txt'

    open(newunit=iocrca,file=fcrca,iostat=ist,form='formatted',action='read')
    if(ist /= 0) stop 'file open error.'

    deallocate(fcrca)

    !   Read header
    read(iocrca,*)
    read(iocrca,'(A)') readline
    read(iocrca,*)

    read(readline(3:),*) ntemp, ndens

    allocate(ztemp(ntemp),zdens(ndens),ratedata(ntemp,ndens,3))
    allocate(dx(ntemp,ndens))
    allocate(dy, dxy, mold=dx)
    allocate(urca6(4,4,ntemp,ndens),urcs5(4,4,ntemp,ndens),urcf5(4,4,ntemp,ndens))
 
    !   Reac main body of data
    do in = 1, ndens
       do it = 1, ntemp
          read(iocrca,*) zdum1, zdum2, (ratedata(it,in,i), i=1,3)
          if( in == 1 ) ztemp(it) = zdum1
       end do
       zdens(in) = zdum2
    end do

    templim(1) = minval(ztemp)
    templim(2) = maxval(ztemp)
    denslim(1) = minval(zdens)
    denslim(2) = maxval(zdens)

    !   Compute spline coefficients

    ! A_6; effective recombination coefficients for iz=6
    call spl2d(ztemp,zdens,ratedata(:,:,1),dx,dy,dxy,urca6,ntemp,ntemp,ndens,0,0,ierr)
    if( ierr /= 0 ) stop 'spl2d error in spline_table_carbon_rate_coef_adas; A_6.'
    ! S_5; effective ionization coefficients for iz=5
    call spl2d(ztemp,zdens,ratedata(:,:,2),dx,dy,dxy,urcs5,ntemp,ntemp,ndens,0,0,ierr)
    if( ierr /= 0 ) stop 'spl2d error in spline_table_carbon_rate_coef_adas; S_5.'
    ! F_5; fraction of carbon with the charge state of iz=0 to 5
    call spl2d(ztemp,zdens,ratedata(:,:,3),dx,dy,dxy,urcf5,ntemp,ntemp,ndens,0,0,ierr)
    if( ierr /= 0 ) stop 'spl2d error in spline_table_carbon_rate_coef_adas; F_5.'

    deallocate(ratedata,dx,dy,dxy)

  end subroutine spline_table_carbon_rate_coef_adas

! ---

  subroutine deallocate_spline_table_carbon_rate_coef_adas

    deallocate(ztemp,zdens,urca6,urcs5,urcf5)

  end subroutine deallocate_spline_table_carbon_rate_coef_adas

!***********************************************************************
!
!   Effective recombination rate coefficients of carbon
!     from hexavalent to pentavalent, taken from OpenADAS
! 
!     Inputs (real*8): TekeV   : electron temperature [keV]
!                    : dene    : Hydrogen neutral temperature [eV]
!     Output (real*8): RateCoef_C_RC_A6 : Effective recombination maxwellian
!                                         rate coefficient [m^3/s]
!
!***********************************************************************

  function RateCoef_C_RC_A6(TekeV, dene20)
    use libspl2d, only : spl2df
    ! input
    real(8), intent(in) :: TekeV & ! [keV]
         &               , dene20  ! [10^{20} m^{-3}]
    ! output
    real(8) :: RateCoef_C_RC_A6
    ! local
    integer(4) :: ierr
    real(8) :: TeeV, dene

    TeeV = TekeV  * 1.d3
    dene = dene20 * 1.d20

    if( TeeV < templim(1) .or. TeeV > templim(2) ) then  
       write(6,*) 'Function RateCoef_RC_A6: Out of e.energy range. Te [eV]=', TeeV
    else if( dene < denslim(1) .or. dene > denslim(2) ) then
       write(6,*) 'Function RateCoef_RC_A6: Out of e.density range. Ne [m^-3]=', dene
    end if

    call spl2df(TeeV,dene,RateCoef_C_RC_A6,ztemp,zdens,urca6,ntemp,ntemp,ndens,ierr)

  end function RateCoef_C_RC_A6

!***********************************************************************
!
!   Effective ionization rate coefficients of carbon
!     from pentavalent to hexavalent, taken from OpenADAS
! 
!     Inputs (real*8): TekeV   : electron temperature [keV]
!                    : dene    : Hydrogen neutral temperature [eV]
!     Output (real*8): RateCoef_C_IZ_S5 : Effective ionization maxwellian
!                                         rate coefficient [m^3/s]
!
!***********************************************************************

  function RateCoef_C_IZ_S5(TekeV, dene20)
    use libspl2d, only : spl2df
    ! input
    real(8), intent(in) :: TekeV & ! [keV]
         &               , dene20  ! [10^{20} m^{-3}]
    ! output
    real(8) :: RateCoef_C_IZ_S5
    ! local
    integer(4) :: ierr
    real(8) :: TeeV, dene

    TeeV = TekeV  * 1.d3
    dene = dene20 * 1.d20

    if( TeeV < templim(1) .or. TeeV > templim(2) ) then  
       write(6,*) 'Function RateCoef_IZ_S5: Out of e.energy range. Te [eV]=', TeeV
    else if( dene < denslim(1) .or. dene > denslim(2) ) then
       write(6,*) 'Function RateCoef_IZ_S5: Out of e.density range. Ne [m^-3]=', dene
    end if

    call spl2df(TeeV,dene,RateCoef_C_IZ_S5,ztemp,zdens,urcs5,ntemp,ntemp,ndens,ierr)

  end function RateCoef_C_IZ_S5

!***********************************************************************
!
!   Fraction of carbon with the charge state of iz=0 to 5
!     based on the superstaged strategy, taken from OpenADAS
! 
!     Inputs (real*8): TekeV   : electron temperature [keV]
!                    : dene    : Hydrogen neutral temperature [eV]
!     Output (real*8): RateCoef_RC_A6 : Effective ionization maxwellian
!                                       rate coefficient [m^3/s]
!
!***********************************************************************

  function Fraction_C0toC5(TekeV, dene20)
    use libspl2d, only : spl2df
    ! input
    real(8), intent(in) :: TekeV & ! [keV]
         &               , dene20  ! [10^{20} m^{-3}]
    ! output
    real(8) :: Fraction_C0toC5
    ! local
    integer(4) :: ierr
    real(8) :: TeeV, dene

    TeeV = TekeV  * 1.d3
    dene = dene20 * 1.d20

    if( TeeV < templim(1) .or. TeeV > templim(2) ) then  
       write(6,*) 'Function RateCoef_IZ_S5: Out of e.energy range. Te [eV]=', TeeV
    else if( dene < denslim(1) .or. dene > denslim(2) ) then
       write(6,*) 'Function RateCoef_IZ_S5: Out of e.density range. Ne [m^-3]=', dene
    end if

    call spl2df(TeeV,dene,Fraction_C0toC5,ztemp,zdens,urcf5,ntemp,ntemp,ndens,ierr)

  end function Fraction_C0toC5

!***********************************************************************
!
!   Spline Table of rate coefficients of beam ionization and charge exchange
!     with Chebyshev polynomials fit
!
!     valid for background Ti = 1eV => 2*10^4eV
!               beam energy   = 10keV/amu => 500keV/amu
! 
!     [1] C.F. Barnett, "ATOMIC DATA FOR FUSION VOLUME 1", ORNL-6086 (1990)
!
!     Common (real*8): Earray : Hydrogenic ion temperature range from Emin to Emax
!                      Ebdata : Discrete beam energy
!                      ucxb   : 2D spline coefficients for charge exchange
!                      uizb   : 2D spline coefficients for ionization by protons
!     Internal (real*8): acxb : Table of Maxwellian rate coefficient in A-27
!                                                   ^^^^^^^^^^^^^^^^
!                        aizb : Table of Maxwellian rate coefficient in D-10
!                                                   ^^^^^^^^^^^^^^^^
!***********************************************************************

  subroutine spline_table_beam_rate_coef
    use libspl2d, only : spl2d
    real(8), dimension(0:6,jmax) :: acxb, aizb
    real(8), dimension(:,:), allocatable :: ratedata, dx, dy, dxy
    integer(4) :: i, j, ierr
    real(8) :: Emax, Emin, dElin, Elin

    !  H + H^+ -> H + H^+
    !  Charge exchange (A-27)
    !  Beam - Maxwellian rate coefficients

    acxb(0:6,1:jmax) = &
         & reshape((/ -3.23588d+1, -2.28428d-1, -1.47539d-1, -7.73418d-2  &
                & ,   -3.29002d-2, -1.09141d-2, -2.93158d-3  & ! 10keV/amu
                & ,   -3.27829d+1, -2.26775d-1, -1.38750d-1, -7.04432d-2  &
                & ,   -3.03096d-2, -1.00013d-2, -1.47228d-3  & ! 20keV/amu
                & ,   -3.39547d+1, -1.36550d-1, -7.15786d-2, -3.65682d-2  &
                & ,   -2.10115d-2, -1.11355d-2, -4.14062d-3  & ! 40keV/amu
                & ,   -3.60736d+1,  1.09084d-1,  7.35296d-2,  2.31264d-2  &
                & ,   -7.74992d-3, -1.38165d-2, -8.49906d-3  & ! 70keV/amu
                & ,   -3.81082d+1,  3.09107d-1,  1.95996d-1,  8.14408d-2  &
                & ,    1.35787d-2, -9.21346d-3, -1.02081d-2  & ! 100keV/amu
                & ,   -4.29608d+1,  3.96726d-1,  2.99000d-1,  1.81207d-1  &
                & ,    8.74182d-2,  3.26701d-2,  7.33350d-3  & ! 200keV/amu
                & ,   -5.17564d+1,  3.37283d-1,  2.43298d-1,  1.48098d-1  &
                & ,    7.72749d-2,  3.54105d-2,  1.56884d-2  & ! 500keV/amu
                &   /), shape(acxb) )

    !  H + H^+ -> H^+ + H^ + e^-
    !  Ionization by protons (D-11)
    !  Beam - Maxwellian rate coefficients

    aizb(0:6,1:jmax) = &
         & reshape((/ -3.82968d+1,  1.12152d00,  6.44681d-1,  1.00717d-1  &
                & ,   -1.00166d-1, -3.95957d-2, -2.65964d-3  & ! 10keV/amu
                & ,   -3.57290d+1,  3.56274d-1,  2.02092d-1,  6.20743d-2  &
                & ,   -2.76894d-3, -1.15918d-2, -4.10154d-3  & ! 20keV/amu
                & ,   -3.41679d+1, -1.62401d-2,  3.246363d-4, 7.54856d-3  &
                & ,    8.28055d-3,  4.20183d-3,  1.24748d-3  & ! 40keV/amu
                & ,   -3.37606d+1, -4.00418d-2, -2.58461d-2, -1.20904d-2  &
                & ,   -3.36001d-3,  2.548592d-4, 9.663075d-4 & ! 70keV/amu
                & ,   -3.37204d+1, -2.56605d-2, -1.80094d-2, -1.07132d-2  &
                & ,   -5.35053d-3, -2.08862d-3, -4.645604d-4 & ! 100keV/amu
                & ,   -3.38801d+1, -1.46171d-2, -9.02218d-3, -4.77874d-3  &
                & ,   -2.29897d-3, -9.828762d-4,-3.658987d-4 & ! 200keV/amu
                & ,   -3.44565d+1, -7.913823d-4,-5.710595d-4,-3.486254d-4 &
                & ,   -1.820028d-4,-8.154704d-5,-3.512063d-5 & ! 500keV/amu
                &   /), shape(aizb) )

    EBdata(1:jmax) = (/ 1.d4, 2.d4, 4.d4, 7.d4, 1.d5, 2.d5, 5.d5 /)

    Emin = 1.d0 ! [eV]
    Emax = 2.d4 ! [eV]

    allocate(Earray(iemax))
    allocate(ucxb(4,4,iemax,jmax),uizb(4,4,iemax,jmax))

    allocate(ratedata(iemax,jmax))
    allocate(dx, dy, dxy, mold=ratedata)

    ! --- Equally-spaced in logspace ---
    dElin = (log10(Emax) - log10(Emin)) / (iemax - 1)

    Elin      = log10(Emin)
    Earray(1) = 1.d0**Elin
    do i = 2, iemax
       Elin      = Elin + dElin
       Earray(i) = 1.d1**Elin
    end do

    ! --- Charge exchange ---

    do j = 1, jmax
       do i = 1, iemax
          ratedata(i,j) = calc_ratecoef(Earray(i),Emax,Emin,acxb(:,j))
       end do
    end do

    !   Compute spline coefficients

    call spl2d(Earray,EBdata,ratedata,dx,dy,dxy,ucxb,iemax,iemax,jmax,0,0,ierr)
    if( ierr /= 0 ) stop 'spl2d error in spline_table_beam_rate_coef; cx.'

    ! --- Ionization ---

    do j = 1, jmax
       do i = 1, iemax
          ratedata(i,j) = calc_ratecoef(Earray(i),Emax,Emin,aizb(:,j))
       end do
    end do

    !   Compute spline coefficients

    call spl2d(Earray,EBdata,ratedata,dx,dy,dxy,uizb,iemax,iemax,jmax,0,0,ierr)
    if( ierr /= 0 ) stop 'spl2d error in spline_table_beam_rate_coef; iz.'

    deallocate(ratedata,dx,dy,dxy)

  end subroutine spline_table_beam_rate_coef

! ---

  subroutine deallocate_spline_table_beam_rate_coef

    if(allocated(Earray)) deallocate(Earray,ucxb,uizb)

  end subroutine deallocate_spline_table_beam_rate_coef

!***********************************************************************
!
!   Spline Table of rate coefficients of charge exchange of singly ionized carbon
!     with hydrogen neutral by Chebyshev polynomials fit
!
!     valid for background Ti =  1eV => 20keV
!               C^+ Tz        = 48eV => 20keV
! 
!     [2] R.A. Phaneuf et al, "ATOMIC DATA FOR FUSION VOLUME 5", ORNL-6090 (1987)
!
!     Common (real*8): ECarray : Carbon ion temperature range from Emin to Emax
!                      EHdata  : Discrete hydrogen temperature
!                      ucxc    : 2D spline coefficients for charge exchange
!     Internal (real*8): acxc  : Table of Maxwellian rate coefficient in 1-14
!                                                    ^^^^^^^^^^^^^^^^
!***********************************************************************

  subroutine spline_table_carbon_rate_coef
    use libspl2d, only : spl2d
    real(8), dimension(0:6,jmax) :: acxc
    real(8), dimension(:,:), allocatable :: ratedata, dx, dy, dxy
    integer(4) :: i, j, ierr
    real(8) :: Emax, Emin, dElin, Elin

    !  C^+ + H -> C + H^+
    !  Charge exchange (1-14)
    !  Maxwellian - Maxwellian rate coefficients

    acxc(0:6,1:jmax) = &
       & reshape((/  2.542d-8,  2.076d-8,  1.143d-8,  3.961d-9  &
       &         ,   4.393d-10,-3.370d-10,-1.523d-10 & ! 1eV
       &         ,   2.589d-8,  2.075d-8,  1.130d-8,  3.776d-9  &
       &         ,   2.705d-10,-4.344d-10,-1.885d-10 & ! 10eV
       &         ,   3.090d-8,  2.072d-8,  1.048d-8,  2.995d-9  &
       &         ,  -1.066d-11,-4.362d-10,-1.755d-10 & ! 100eV
       &         ,   8.273d-8,  1.082d-8,  5.220d-9,  1.524d-9  &
       &         ,   1.318d-10,-1.100d-10,-6.651d-11 & ! 1000eV
       &         ,   1.409d-7,  1.055d-9,  5.203d-10, 1.559d-10 &
       &         ,   1.439d-11,-1.384d-11,-9.706d-12 & ! 5000eV
       &         ,   1.438d-7, -4.313d-10,-2.618d-10,-1.280d-10 &
       &         ,  -5.361d-11,-2.080d-11,-7.884d-12 & ! 10000eV
       &         ,   1.239d-7, -7.127d-10,-4.050d-10,-1.734d-10 &
       &         ,  -5.861d-11,-1.623d-11,-3.746d-12 & ! 20000eV
       &           /), shape(acxc))

    EHdata(1:jmax) = (/ 1.d0, 1.d1, 1.d2, 1.d3, 5.d3, 1.d4, 2.d4 /)

    Emin = 4.8d0 ! [eV]
    Emax = 2.0d4 ! [eV]

    allocate(ECarray(iemax))
    allocate(ucxc(4,4,iemax,jmax))

    allocate(ratedata(iemax,jmax))
    allocate(dx, dy, dxy, mold=ratedata)

    ! --- Equally-spaced in logspace ---
    dElin = (log10(Emax) - log10(Emin)) / (iemax - 1)

    Elin      = log10(Emin)
    ECarray(1) = 1.d0**Elin
    do i = 2, iemax
       Elin      = Elin + dElin
       ECarray(i) = 1.d1**Elin
    end do

    ! --- Charge exchange ---

    do j = 1, jmax
       do i = 1, iemax
          ratedata(i,j) = calc_ratecoef(ECarray(i),Emax,Emin,acxc(:,j))
       end do
    end do

    !   Compute spline coefficients

    call spl2d(ECarray,EHdata,ratedata,dx,dy,dxy,ucxc,iemax,iemax,jmax,0,0,ierr)
    if( ierr /= 0 ) stop 'spl2d error in spline_table_carbon_rate_coef; cx.'

    deallocate(ratedata,dx,dy,dxy)

  end subroutine spline_table_carbon_rate_coef

! ---

  subroutine deallocate_spline_table_carbon_rate_coef

    if(allocated(ECarray)) deallocate(ECarray,ucxc)

  end subroutine deallocate_spline_table_carbon_rate_coef

! ----------------------------------------------------------------------

  function calc_ratecoef(E,Emax,Emin,a)
    ! input
    real(8), intent(in) :: E, Emax, Emin
    real(8), dimension(0:), intent(in) :: a
    ! output
    real(8) :: calc_ratecoef
    ! local
    real(8) :: lnsigE, x
    real(8) :: cm3tom3 = 1.d-6

    x = ( ( log(E) - log(Emin) ) - ( log(Emax) - log(E) ) ) / ( log(Emax) - log(Emin) )

    lnsigE = 0.5d0 * a(0) + chebypoly(x,a)

    calc_ratecoef = exp(lnsigE) * cm3tom3

  end function calc_ratecoef

! ---

  function chebypoly(x,a) result(AT)
    ! input
    real(8), intent(in) :: x
    real(8), dimension(0:), intent(in) :: a
    ! output
    real(8) :: AT
    ! local
    integer(4) :: i, j, imax
    real(8) :: T(0:8)

    imax = size(a) - 2

    T(0) = 1.d0
    T(1) = x
    AT = a(1) * T(1)

    do i = 1, imax
       j = i + 1
       ! Chebyshev polynomials
       T(j) = 2.d0 * x * T(j-1) - T(j-2)
       ! Sum A_jT_j
       AT = AT + a(j) * T(j)
    end do

  end function chebypoly

!***********************************************************************
!
!   Charge exchange rate coefficients for beam H + Maxwellian H^+ -> H + H^+
!
!     valid for background Ti = 1eV => 2*10^4eV
!               beam energy   = 10keV/amu => 500keV/amu
! 
!     Inputs (real*8): TieV    : Hydrogenic ion temperature [eV]
!                    : EbeVamu : Beam energy [eV]
!     Output (real*8): SiVcxb  : Charge exchange maxwellian rate coefficient with beam [m^3/s]
!
!***********************************************************************

  function SiVcxB(TieV, EbeVamu)
    use libspl2d, only : spl2df
    ! input
    real(8), intent(in) :: TieV, EbeVamu ! [eV]
    ! output
    real(8) :: SiVcxB
    ! local
    integer(4) :: ierr

    if( TieV < Earray(1) .or. TieV > Earray(iemax) ) then  
       write(6,*) 'Function SiVcxB: Out of energy range. Ti [eV]=', TieV
    else if( EbeVamu < EBdata(1) .or. EbeVamu > EBdata(jmax) ) then
       write(6,*) 'Function SiVcxB: Out of energy range. Eb [eV/amu]=', EbeVamu
    end if

    call spl2df(TieV,EbeVamu,SiVcxB,Earray,EBdata,ucxb,iemax,iemax,jmax,ierr)

  end function SiVcxB

!***********************************************************************
!
!   Ionization rate coefficients for beam H + Maxwellian H^+ -> H^+ + H^+ + e^-
!
!     valid for background Ti = 1eV => 2*10^4eV
!               beam energy   = 10keV/amu => 500keV/amu
! 
!     Inputs (real*8): TieV    : Hydrogenic ion temperature [eV]
!                    : EbeVamu : Beam energy [eV]
!     Output (real*8): SiVizb  : Ionization maxwellian rate coefficient with beam [m^3/s]
!
!***********************************************************************

  function SiVizB(TieV, EbeVamu)
    use libspl2d, only : spl2df
    ! input
    real(8), intent(in) :: TieV, EbeVamu ! [eV]
    ! output
    real(8) :: SiVizB
    ! local
    integer(4) :: ierr

    if( TieV < Earray(1) .or. TieV > Earray(iemax) ) then
       write(6,*) 'Function SiVizB: Out of energy range. Ti [eV]=', TieV
    else if( EbeVamu < EBdata(1) .or. EbeVamu > EBdata(jmax) ) then
       write(6,*) 'Function SiVizB: Out of energy range. Eb [eV/amu]=', EbeVamu
    end if

    call spl2df(TieV,EbeVamu,SiVizB,Earray,EBdata,uizb,iemax,iemax,jmax,ierr)

  end function SiVizB

!***********************************************************************
!
!   Charge exchange rate coefficients of singly ionized carbon with hydrogen neutral
!     C^+ + H -> C + H^+
!
!     valid for background Ti =  1eV => 20keV
!               C^+ Tz        = 48eV => 20keV
! 
!     Inputs (real*8): TzeV    : carbon ion temperature [eV]
!                    : TieV    : Hydrogen neutral temperature [eV]
!     Output (real*8): SiVcxC  : Charge exchange maxwellian rate coefficient [m^3/s]
!
!***********************************************************************

  function SiVcxC(TzeV, TieV)
    use libspl2d, only : spl2df
    ! input
    real(8), intent(in) :: TzeV, TieV ! [eV]
    ! output
    real(8) :: SiVcxC
    ! local
    integer(4) :: ierr
    real(8) :: Emax = 2.d4, Emin = 48.d0
    real(8), dimension(7) :: &
         & a = (/ 7.115d8, 3.946d8, -2.767d-9, -1.080d-8, -2.262d-9, 1.163d-9, 3.938d-10/)

    if( TzeV < ECarray(1) .or. TzeV > ECarray(iemax) ) then  
       write(6,*) 'Function SiVcxC: Out of energy range. Tz [eV]=', TzeV
    else if( TieV < EHdata(1) .or. TieV > EHdata(jmax) ) then
       write(6,*) 'Function SiVcxC: Out of energy range. Ti [eV]=', TieV
    end if

    if( TzeV /= TieV ) then
       call spl2df(TzeV,TieV,SiVcxC,ECarray,EHdata,ucxc,iemax,iemax,jmax,ierr)
    else
       SiVcxC = calc_ratecoef(TzeV,Emax,Emin,a)
    end if

  end function SiVcxC

!***********************************************************************
!
!   Ionization rate coefficients by electron impact of e^- + C^q+ -> C^(q+1)+ + 2e^-
!     with Chebyshev polynomials fit
!
!     valid for background Te = 2eV - 40eV => 2*10^4eV
! 
!     [2] R.A. Phaneuf et al, "ATOMIC DATA FOR FUSION VOLUME 5", ORNL-6090 (1987)
!
!     Inputs (real*8): ZC    : Carbon charge state
!                      TeeV  : Electron temperature [eV]
!     Output (real*8): SiVizC: Maxwellian rate coefficient of carbon ionization
!                                         ^^^^^^^^^^^^^^^^
!***********************************************************************

  function SiVizC(ZC, TeeV)
    ! input
    real(8), intent(in) :: ZC, TeeV ! [eV]
    ! output
    real(8) :: SiVizC
    ! local
    integer(4) :: ierr, intZC
    real(8) :: Emax, Emin

    Emax = 2.d4
    if( TeeV > Emax ) then
       write(6,'(A,ES12.4,A)') 'Function SiVizC: Out of energy range. Te [eV]=', TeeV &
            &                , ' > 20.0 [keV]' 
    end if

    intZC = int(ZC + 1.d-6)
    select case(intZC)
    case(0)
       Emin = 2.d0
       if( TeeV < Emin ) &
            & write(6,'(A,ES12.4,A)') 'Function SiVizC: Out of energy range. Te [eV]=', TeeV &
            &                       , ' < 2.0 [eV]'
       
    case(1)
       Emin = 4.d0
       if( TeeV < Emin ) &
            & write(6,'(A,ES12.4,A)') 'Function SiVizC: Out of energy range. Te [eV]=', TeeV &
            &                       , ' < 4.0 [eV]'
       
    case(2:3)
       Emin = 7.d0
       if( TeeV < Emin ) &
            & write(6,'(A,ES12.4,A)') 'Function SiVizC: Out of energy range. Te [eV]=', TeeV &
            &                       , ' < 7.0 [eV]'
       
    case(4:5)
       Emin = 4.d1
       if( TeeV < Emin ) &
            & write(6,'(A,ES12.4,A)') 'Function SiVizC: Out of energy range. Te [eV]=', TeeV &
            &                       , ' < 40.0 [eV]'
       
    end select

    SiVizC = calc_ratecoef(TeeV,Emax,Emin,aizc(:,intZC))

  end function SiVizC

end module mod_cross_section
