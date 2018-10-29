Program testmmm

!                 >>>>>   N O T E   <<<<<
! {{{ and }}} are VIM folding markers. DO NOT REMOVE THEM.

Use modmmm7_1

Implicit None

! SPECIFICtiON {{{

!------- Constants -------------------------------------------------

Integer, Parameter :: R8 = SELECTED_REAL_KIND( 12, 100 )

!.. Machine Epsilon
Real(R8), Parameter :: &
   zepslon = 2._R8**(-53) ! As defined by IEEE 754-2008 standard

Integer, Parameter :: &
   NPMAX  = 200, &! Maximum data points allowed in arrays
   hfIn   = 34,  &! Input file handle
   hfOut  = 35,  &! Output file handle
   hfDebug= 36    ! Diagnostic output

!.. Constants for missing input dectection
Real(R8), Parameter :: BADREAL = -1E10_R8
Integer, Parameter :: BADINT = -1000000

!------- Input Variables -------------------------------------------
! The content of the input variables follows their corresponding arguments
! of the mmm7_1 subroutine. See modmmm7_1.f90 for more details.
Integer :: lprint, npoints, nerr

Real(R8), Dimension(3) :: &
   cmodel  ! Internal model weights

Real(R8), Dimension(MAXNOPT) :: &
   cW20, cDBM, cETG ! To be passed to cswitch

Integer, Dimension(MAXNOPT) :: &
   lW20, lDBM, lETG ! To be passed to lswitch

!.. Input profile variables of the first kind
Real(R8), Dimension(NPMAX) :: &
   rmin   ,&! half-width of the mgnetic surface [m]
   rmaj   ,&! major radius to geometric center of the mgnetic surface [m]
   elong  ,&! local elongtion
   ne     ,&! electron density [m^-3]
   nh     ,&! thermal hydrogenic ion densitie [m^-3]
   nz     ,&! impurity ion densitie [m^-3]
   nf     ,&! density from fast (non-thermal) ions [m^-3]
   zeff   ,&! Z_eff
   te     ,&! T_e (electron temperature) [keV]
   ti     ,&! T_i (temperature of thermal ions) [keV]
   q      ,&! mgnetic q-value
   btor   ,&! ( R B_tor ) / rmaj(jz)  [Tesla]
   zimp   ,&! mean charge of impurities
   aimp   ,&! mean atomic mass of impurities
   ahyd   ,&! mean atomic mass of hydrogen ions
   aimass ,&! mean atomic mass of thermal ions
   wexbs  ,&! ExB shearing rate in [rad/s]
   gne    ,&! -R ( d n_e / d r ) / n_e
   gni    ,&! -R ( d n_i / d r ) / n_i
   gnh    ,&! -R ( d n_h / d r ) / n_h
   gnz    ,&! -R ( d Z n_Z / d r ) / ( Z n_Z )
   gte    ,&! -R ( d T_e / d r ) / T_e
   gti    ,&! -R ( d T_i / d r ) / T_i
   gq     ,&!  R ( d q   / d r ) / q
   gvtor  ,&! Normalized toroidal velocity gradient (r/v_tor)*(dv_tor/dr)
   vtor   ,&! Toroidal velocity [m/s]
   gvpol  ,&! Normalized poloidal velocity gradient (r/v_pol)*(dv_pol/dr)
   vpol   ,&! Poloidal velocity [m/s]
   gvpar  ,&! Normalized parallel velocity gradient (r/v_par)*(dv_par/dr)
   vpar     ! Parallel velocity [m/s]

namelist /testmmm_input_1stkind/ &
   cmodel, npoints, lprint,                         & 
   cW20, cDBM, cETG, lW20, lDBM, lETG,              &
   rmin, rmaj, elong, ne, nh, nz, nf, zeff, te, ti, &
   q, btor, zimp, aimp, ahyd, aimass, wexbs,        &
   gne, gni, gnh, gnz, gte, gti, gq,                &
   gvtor, vtor, gvpol, vpol, gvpar, vpar

!.. Input variables of the second kind
Real(R8) :: &
   k_rminor ,&! Minor radius of plasma [m]
   k_rmajor ,&! Major radius at plasma axis [m]
   k_elong  ,&! Elongtion
   k_denmin ,&! Minimal electron density [m^-3]
   k_temin  ,&! Minimal electron temperature [keV]
   k_btor   ,&! Toroidal mgnetic field [Tesla]
   k_amassh ,&! Mean atomic mass of hydrogenic ions
   k_amassz   ! Mean atomic mass of impurities
!
!.. Parameters for generting profiles
!   See documenttion for more explantion
Real(R8) :: &
   denhaxis, denhedge, denhexp, &! Ion density
   denzaxis, denzedge, denzexp, &! Average impurity density
   chrzaxis, chrzedge, chrzexp, &! Average impurity charge
   denfaxis, denfedge, denfexp, &! Fast ion density
   chrfaxis, chrfedge, chrfexp, &! Average charge of super-thermal ions
   teaxis,   teedge,   teexp,   &! Electron temperature
   tiaxis,   tiedge,   tiexp,   &! Ion temperature
   qaxis,    qedge,    qexp,    &! Mgnetic q
   wexbmax,  xwexbinn, xwexbout  ! ExB flow shear rate

namelist /testmmm_input_2ndkind/ &
   cmodel, npoints, lprint,               &
   cW20, cDBM, cETG, lW20, lDBM, lETG,    &
   k_rminor, k_rmajor, k_elong, k_denmin, &
   k_temin,  k_amassh, k_btor,  k_amassz, &
   denhaxis, denhedge, denhexp,           &
   denzaxis, denzedge, denzexp,           &
   denfaxis, denfedge, denfexp,           & 
   chrzaxis, chrzedge, chrzexp,           &
   chrfaxis, chrfedge, chrfexp,           &
   teaxis,   teedge,   teexp,             &
   tiaxis,   tiedge,   tiexp,             &
   qaxis,    qedge,    qexp,              &
   wexbmax,  xwexbinn, xwexbout 

!------- Output Variables ------------------------------------------
! The content of the output variables follows their corresponding arguments
! of the mmm7_1 subroutine. See modmmm7_1.f90 for more details.

Real(R8), Dimension(NPMAX) :: &
   xti, xdi, xte, xdz, xvt, xvp, gammaDBM, omegaDBM,      &
   xtiW20, xdiW20, xteW20, xtiDBM, xdiDBM, xteDBM, xteETG

Real(R8), Dimension(4,NPMAX) :: &
   gammaW20, omegaW20

Real(R8), Dimension(6,NPMAX) :: &
   vconv, vflux

!------- Local Variables ------------------------------------------

!.. Internal variables
Real(R8), Dimension(NPMAX) :: &
   zchrgz, zchrgf, zxb, zdensf, zdensi, zgrdnf

Real(R8) :: &
   zdx, zxfact

Real(R8) :: cswitch(MAXNOPT,4)

Integer :: lswitch(MAXNOPT,4)

Integer :: jr

Integer :: input_kind

Integer :: errcount

Character(len=255) :: strbuf
!}}}

!------- EXECUTION BODY --------------------------------------------

Open( hfIn, file='input', &
   form='formatted', status='old', iostat=nerr)
if ( nerr /= 0 ) stop "ERROR: Cannot open input file!"

Read( hfIn,'(A)') strbuf
Rewind( hfIn )

lprint = 0
npoints = 0
cmodel = BADREAL
cW20 = BADREAL; cDBM = BADREAL; cETG = BADREAL
lW20 = BADINT;  lDBM = BADINT;  lETG = BADINT

If ( index( strbuf, '1stkind') > 0 ) Then

   ! Input of the 1st kind {{{
   input_kind = 1

   Write(*,'("Input of the first kind (values) is detected. Processing...")')

   rmin  = BADREAL; rmaj  = BADREAL; elong = BADREAL; ne    = BADREAL
   nh    = BADREAL; nz    = BADREAL; nf    = BADREAL; zeff  = BADREAL
   te    = BADREAL; ti    = BADREAL; q     = BADREAL; btor  = BADREAL
   zimp  = BADREAL; aimp  = BADREAL; ahyd  = BADREAL; aimass= BADREAL
   wexbs = BADREAL; gne   = BADREAL; gni   = BADREAL; gnh   = BADREAL
   gnz   = BADREAL; gte   = BADREAL; gti   = BADREAL; gq    = BADREAL
   gvtor = BADREAL; vtor  = BADREAL; gvpol = BADREAL; vpol  = BADREAL
   gvpar = BADREAL; vpar  = BADREAL

   Read( hfIn, NML=testmmm_input_1stkind )

   errcount = 0

   Call check_input( rmin(1)  , "rmin"  )
   Call check_input( rmaj(1)  , "rmaj"  )
   Call check_input( elong(1) , "elong" )
   Call check_input( ne(1)    , "ne"    )
   Call check_input( nh(1)    , "nh"    )
   Call check_input( nz(1)    , "nz"    )
   Call check_input( nf(1)    , "nf"    )
   Call check_input( zeff(1)  , "zeff"  )
   Call check_input( te(1)    , "te"    )
   Call check_input( ti(1)    , "ti"    )
   Call check_input( q(1)     , "q"     )
   Call check_input( btor(1)  , "btor"  )
   Call check_input( zimp(1)  , "zimp"  )
   Call check_input( aimp(1)  , "aimp"  )
   Call check_input( ahyd(1)  , "ahyd"  )
   Call check_input( aimass(1), "aimass")
   Call check_input( wexbs(1) , "wexbs" )
   Call check_input( gq(1)    , "gq"    )
   Call check_input( gvtor(1) , "gvtor" )
   Call check_input( vtor(1)  , "vtor"  )
   Call check_input( gvpol(1) , "gvpol" )
   Call check_input( vpol(1)  , "vpol"  )
   Call check_input( gvpar(1) , "gvpar" )
   Call check_input( vpar(1)  , "vpar"  )

   If ( errcount > 0 ) Stop

   !.. Calcultes the gradients if they < -100
100 Format ("** ERROR ** Too few radial points to calculte gradients!")

   If ( gte(1) < -1E2_R8 ) Then

      If ( npoints < 2) Then
         Write(*,100)
         Stop
      End If

      gte(1) = 0._R8
      Do jr=2,npoints-1
         gte(jr) = - rmaj(jr)/te(jr)* &
            ( te(jr+1)-te(jr-1) )/( rmin(jr+1)-rmin(jr-1) )
      End Do
      gte(npoints) = - rmaj(npoints)/te(npoints)* &
         ( te(npoints)-te(npoints-1) )/&
         ( rmin(npoints)-rmin(npoints-1) )
   End If

   If ( gti(1) < -1E2_R8 ) Then

      If ( npoints < 2) Then
         Write(*,100)
         Stop
      End If

      gti(1) = 0._R8
      Do jr=2,npoints-1
         gti(jr) = - rmaj(jr)/ti(jr)* &
            ( ti(jr+1)-ti(jr-1) )/( rmin(jr+1)-rmin(jr-1) )
      End Do
      gti(npoints) = - rmaj(npoints)/ti(npoints)* &
         ( ti(npoints)-ti(npoints-1) )/&
         ( rmin(npoints)-rmin(npoints-1) )
   End If

   If ( gne(1) < -1E2_R8 ) Then

      If ( npoints < 2) Then
         Write(*,100)
         Stop
      End If

      gne(1) = 0._R8
      Do jr=2,npoints-1
         gne(jr) = - rmaj(jr)/ne(jr)* &
            ( ne(jr+1)-ne(jr-1) )/( rmin(jr+1)-rmin(jr-1) )
      End Do
      gne(npoints) = - rmaj(npoints)/ne(npoints)* &
         ( ne(npoints)-ne(npoints-1) )/&
         ( rmin(npoints)-rmin(npoints-1) )
   End If

   If ( gnh(1) < -1E2_R8 ) Then

      If ( npoints < 2) Then
         Write(*,100)
         Stop
      End If

      gnh(1) = 0._R8
      Do jr=2,npoints-1
         gnh(jr) = - rmaj(jr)/nh(jr)* &
            ( nh(jr+1)-nh(jr-1) )/( rmin(jr+1)-rmin(jr-1) )
      End Do
      gnh(npoints) = - rmaj(npoints)/nh(npoints)* &
         ( nh(npoints)-nh(npoints-1) )/&
         ( rmin(npoints)-rmin(npoints-1) )
   End If

   If ( gnz(1) < -1E2_R8 ) Then

      If ( npoints < 2) Then
         Write(*,100)
         Stop
      End If

      gnz(1) = 0._R8
      Do jr=2,npoints-1
         gnz(jr) = - rmaj(jr)/nz(jr)* &
            ( nz(jr+1)-nz(jr-1) )/( rmin(jr+1)-rmin(jr-1) )
      End Do
      gnz(npoints) = - rmaj(npoints)/nz(npoints)* &
         ( nz(npoints)-nz(npoints-1) )/&
         ( rmin(npoints)-rmin(npoints-1) )
   End If

   If ( gni(1) < -1E2_R8 ) Then
      gni=gnh+gnz
   End If
!}}}

Elseif ( index( strbuf, '2ndkind') > 0 ) Then

   ! Input of the 2nd kind {{{

   input_kind = 2

   Write(*,'("Input of the second kind (shape functions) is detected. Processing...")')

   k_rminor = BADREAL; k_rmajor = BADREAL; k_elong  = BADREAL
   k_denmin = BADREAL; k_temin  = BADREAL; k_btor   = BADREAL
   k_amassh = BADREAL; k_amassz = BADREAL
   denhaxis = BADREAL; denhedge = BADREAL; denhexp  = BADREAL
   denzaxis = BADREAL; denzedge = BADREAL; denzexp  = BADREAL
   denfaxis = BADREAL; denfedge = BADREAL; denfexp  = BADREAL
   chrzaxis = BADREAL; chrzedge = BADREAL; chrzexp  = BADREAL
   chrfaxis = BADREAL; chrfedge = BADREAL; chrfexp  = BADREAL
   teaxis   = BADREAL; teedge   = BADREAL; teexp    = BADREAL
   tiaxis   = BADREAL; tiedge   = BADREAL; tiexp    = BADREAL
   qaxis    = BADREAL; qedge    = BADREAL; qexp     = BADREAL
   wexbmax  = BADREAL; xwexbinn = BADREAL; xwexbout = BADREAL

   read( hfIn, NML=testmmm_input_2ndkind )

   errcount = 0

   Call check_input( k_rminor , "k_rminor" )
   Call check_input( k_rmajor , "k_rmajor" )
   Call check_input( k_elong  , "k_elong " )
   Call check_input( k_denmin , "k_denmin" )
   Call check_input( k_temin  , "k_temin " )
   Call check_input( k_btor   , "k_btor  " )
   Call check_input( k_amassh , "k_amassh" )
   Call check_input( k_amassz , "k_amassz" )
   Call check_input( denhaxis , "denhaxis" )
   Call check_input( denhedge , "denhedge" )
   Call check_input( denhexp  , "denhexp " )
   Call check_input( denzaxis , "denzaxis" )
   Call check_input( denzedge , "denzedge" )
   Call check_input( denzexp  , "denzexp " )
   Call check_input( denfaxis , "denfaxis" )
   Call check_input( denfedge , "denfedge" )
   Call check_input( denfexp  , "denfexp " )
   Call check_input( chrzaxis , "chrzaxis" )
   Call check_input( chrzedge , "chrzedge" )
   Call check_input( chrzexp  , "chrzexp " )
   Call check_input( chrfaxis , "chrfaxis" )
   Call check_input( chrfedge , "chrfedge" )
   Call check_input( chrfexp  , "chrfexp " )
   Call check_input( teaxis   , "teaxis  " )
   Call check_input( teedge   , "teedge  " )
   Call check_input( teexp    , "teexp   " )
   Call check_input( tiaxis   , "tiaxis  " )
   Call check_input( tiedge   , "tiedge  " )
   Call check_input( tiexp    , "tiexp   " )
   Call check_input( qaxis    , "qaxis   " )
   Call check_input( qedge    , "qedge   " )
   Call check_input( qexp     , "qexp    " )
   Call check_input( wexbmax  , "wexbmax " )
   Call check_input( xwexbinn , "xwexbinn" )
   Call check_input( xwexbout , "xwexbout" )

   If ( errcount > 0 ) Stop

   !..check to see if input makes sense
   If ( npoints < 5 .or. npoints > NPMAX ) then
     write(*,'("** ERROR ** npoints < 5  .or.  npoints > NPMAX ")')
     stop
   End If

   If ( k_denmin < zepslon ) then
     write(*,'("# k_denmin changed from "E14.6" to 1e7")') k_denmin
     k_denmin = 1E7_R8
   End If

   If ( k_temin < zepslon ) then
     write(*,'("# tenmin changed from "E14.6" to 1e-6")') k_temin
     k_temin = 1E-6_R8
   End If

   !..set up argument list for subroutine theory

   !..equilibrium on zone boundaries (radial points)
   zdx = 1._R8 / (npoints-1)
   Do jr=1,npoints
      zxb(jr)     = ( jr - 1 ) * zdx
      rmin(jr) = ( jr - 1 ) * zdx * k_rminor
      rmaj(jr) = k_rmajor ! Ignore Shafranov Shift
      elong(jr)  = k_elong
   End Do

   !..construct ion density profiles on zone boundaries (radial points)
   Do jr=1,npoints-1
      nh(jr) = denhedge + (denhaxis-denhedge) * (1._R8 - zxb(jr)**2)**denhexp
      nz(jr) = denzedge + (denzaxis-denzedge) * (1._R8 - zxb(jr)**2)**denzexp
      zdensf(jr) = denfedge + (denfaxis-denfedge) * (1._R8 - zxb(jr)**2)**denfexp
      te(jr) = max ( k_temin, teedge + (teaxis-teedge) * (1._R8 - zxb(jr)**2)**teexp )
      ti(jr) = max ( k_temin, tiedge + (tiaxis-tiedge) * (1._R8 - zxb(jr)**2)**tiexp )
      zchrgz(jr) = max (1._R8, chrzedge + (chrzaxis-chrzedge) * (1._R8 - zxb(jr)**2)**chrzexp )
      zchrgf(jr) = max (1._R8, chrfedge + (chrfaxis-chrfedge) * (1._R8 - zxb(jr)**2)**chrfexp )
   End Do

   nh(npoints) = denhedge
   nz(npoints) = denzedge
   zdensf(npoints) = denfedge
   te(npoints) = max ( k_temin, teedge )
   ti(npoints) = max ( k_temin, tiedge )
   zchrgz(npoints) = max ( 1._R8, chrzedge )
   zchrgf(npoints) = max ( 1._R8, chrfedge )

   !..electron and thermal ion densities
   Do jr=1,npoints
     ne(jr) = max ( k_denmin,  nh(jr) &
        + zchrgz(jr) * nz(jr) + zchrgf(jr) * zdensf(jr) )
     zdensi(jr) = max ( k_denmin, nh(jr) + nz(jr) )
     nf(jr) = zchrgf(jr) * zdensf(jr)
   End Do

   !..normalized gradients
   Do jr=1,npoints

     zxfact = max ( zepslon, 1._R8 - zxb(jr)**2 )

     gnh(jr) = rmaj(jr) * 2._R8 * zxb(jr) * denhexp * &
        (denhaxis-denhedge) * zxfact**(denhexp-1._R8) &
        / ( max( k_denmin, nh(jr) ) * rmin(npoints) )

     gnz(jr) = rmaj(jr) * 2._R8 * zxb(jr) * ( &
        (denzaxis-denzedge) * zxfact**(denzexp-1._R8) * &
          denzexp / ( max( k_denmin, nz(jr) ) * rmin(npoints)))

     zgrdnf(jr) = rmaj(jr) * 2._R8 * zxb(jr) * ( &
        (denfaxis-denfedge) * zxfact**(denfexp-1._R8) * &
          denfexp / ( max( k_denmin, zdensf(jr) ) * rmin(npoints)))

     gte(jr) = rmaj(jr) * 2._R8 * zxb(jr) * teexp * &
        (teaxis-teedge) * zxfact**(teexp-1._R8) &
          / ( te(jr) * rmin(npoints) )

     gti(jr) = rmaj(jr) * 2._R8 * zxb(jr) * tiexp * &
        (tiaxis-tiedge) * zxfact**(tiexp-1._R8) &
          / ( ti(jr) * rmin(npoints) )

     gni(jr) = 2._R8 * zxb(jr) * rmaj(jr) * ( &
        + denhexp*(denhaxis-denhedge) * zxfact**(denhexp-1._R8) &
        + denzexp*(denzaxis-denzedge) * zxfact**(denzexp-1._R8) &
        ) / ( rmin(npoints) * zdensi(jr) )

   End Do

   !..normalized electron density gradient
   Do jr=1,npoints
     gne(jr) = ( nh(jr) * gnh(jr) &
        + zchrgz(jr) * nz(jr)*gnz(jr) &
        + zchrgf(jr) * zdensf(jr) * zgrdnf(jr)  ) / ne(jr)
   End Do

   !..average mass and charge
   do jr=1,npoints
      zimp(jr) = zchrgz(jr)
      aimp(jr) = k_amassz
      ahyd(jr) = k_amassh
      aimass(jr) = &
         ( ahyd(jr) * nh(jr) + aimp(jr) * nz(jr) ) &
         / ( nh(jr) + nz(jr) )
      zeff(jr) =  &
         ( nh(jr) + zchrgz(jr)**2 * nz(jr) &
         + zchrgf(jr)**2 * zdensf(jr) ) &
         / ne(jr)
   End Do

   !.. Mgnetics profiles
   do jr=1,npoints
     q(jr) = qaxis + (qedge-qaxis) * abs(zxb(jr))**qexp
     gq(jr) = rmaj(jr) * &
        (qedge-qaxis) * (abs(zxb(jr)))**(qexp-1._R8) &
          / ( q(jr) * rmin(npoints) )
     !qprime = (qedge-qaxis) * qexp * abs(zxb(jr))**(qexp-1D0)
   End Do

   do jr=1,npoints
     btor(jr)  = k_btor
   End Do

   !   Note:  The flow shear rate is given as a fourth-order
   !     polynomial with maximum value wexbmax between
   !     xwexbinn < r/a < xwexbout.
   !
   ! wexbmax    maximum flow shear rate [sec^-1]
   ! xwexbinn   inner r/a cutoff for the flow shear rate
   ! xwexbout   outer r/a cutoff for the flow shear rate
   do jr=1,npoints
     wexbs(jr) = 0._R8
     if ( zxb(jr) > xwexbinn .and.  zxb(jr) < xwexbout ) then
       wexbs(jr) = wexbmax * 16._R8 &
          * ( ( zxb(jr) - xwexbinn ) / ( xwexbout - xwexbinn ) )**2 &
          * ( ( xwexbout - zxb(jr) ) / ( xwexbout - xwexbinn ) )**2
     End If
   End Do
   
   gvtor = 0._R8
   vtor  = 0._R8
   gvpol = 0._R8
   vpol  = 0._R8
   gvpar = 0._R8
   vpar  = 0._R8

!}}}

Else ! Unsupported kind

   write(*,'("** ERROR ** Unsupported input!")')
   Close( unit=hfIn )
   Stop

End If

Close( unit=hfIn )

errcount = 0
Call check_input( cmodel(1), "cmodel" )
If ( errcount > 0 ) Stop

!.. Fill parameter arrays with default values
Call set_mmm7_1_switches( cmmm = cswitch, lmmm = lswitch )

!.. Assign user specified parameters
If ( abs( cW20(1) - BADREAL ) > 1E-6_R8 ) cswitch(1:MAXNOPT,1)=cW20
If ( abs( cDBM(1) - BADREAL ) > 1E-6_R8 ) cswitch(1:MAXNOPT,2)=cDBM
If ( abs( cETG(1) - BADREAL ) > 1E-6_R8 ) cswitch(1:MAXNOPT,3)=cETG

If ( lW20(1) /= BADINT ) lswitch(1:MAXNOPT,1)=lW20
If ( lDBM(1) /= BADINT ) lswitch(1:MAXNOPT,2)=lDBM
If ( lETG(1) /= BADINT ) lswitch(1:MAXNOPT,3)=lETG

Open( hfOut, file='output', &
   form='formatted', status='replace', iostat=nerr)
if ( nerr /= 0 ) stop "ERROR: Cannot creat output file!"

Open( hfDebug, file='testmmm_debug', &
   form='formatted', status='replace', iostat=nerr)
if ( nerr /= 0 ) stop "ERROR: Cannot creat debug file!"

!.. Write input profiles
Write(hfOut,'("#.. Each column from 1 to 7 below represents&
   & an input profile, which was used in the file named ""input""")')
write(hfOut,'("# "I6,30I11)') (jr,jr=1,30)
write(hfOut,'("#Unit"A6,30A11)') "m", "m", &
   " ", "m^-3", "m^-3", "m^-3", "m^-3", " ", &
   "keV", "keV", " ", "Tesla", " ",          &
   " ", " ", " ", "rad/s", " ",              &
   " ", " ", " ", " ", " ", " ",             &
   " ", "m/s", " ", "m/s", " ", "m/s"
write(hfOut,'("# "A9,30A11)') "rmin","rmaj", &
   "elong","ne", "nh","nz","nf","zeff",      &
   "te","ti", "q","btor","zimp",             &
   "aimp","ahyd","aimass","wexbs","gne",     &
   "gni","gnh","gnz","gte","gti","gq",       &
   "gvtor","vtor","gvpol","vpol","gvpar","vpar"

Do jr=1, npoints
   write(hfOut,'(0P2F11.6,1P29E11.3)') &
      rmin(jr), rmaj(jr), elong(jr),       &
      ne(jr), nh(jr), nz(jr),              &
      nf(jr), zeff(jr), te(jr), ti(jr),    &
      q(jr), btor(jr), zimp(jr), aimp(jr), &
      ahyd(jr), aimass(jr), wexbs(jr),     &
      gne(jr), gni(jr), gnh(jr),           &
      gnz(jr), gte(jr), gti(jr), gq(jr),   &
      gvtor(jr), vtor(jr), gvpol(jr),    &
      vpol(jr), gvpar(jr), vpar(jr)
End Do

!.. Initializtion
xti    = 0._R8; xdi    = 0._R8; xte    = 0._R8
xdz    = 0._R8; xvt    = 0._R8; xvp    = 0._R8
xtiW20 = 0._R8; xdiW20 = 0._R8; xteW20 = 0._R8
xtiDBM = 0._R8; xdiDBM = 0._R8; xteDBM = 0._R8; xteETG = 0._R8
gammaW20 = 0._R8; omegaW20 = 0._R8
gammaDBM = 0._R8; omegaDBM = 0._R8

call mmm7_1( &
   rmin   = rmin,   rmaj   = rmaj,   rmaj0  = rmaj(1),    &
   elong  = elong,  ne     = ne,     nh     = nh,         &
   nz     = nz,     nf     = nf,     zeff   = zeff,       &
   te     = te,     ti     = ti,     q      = q,          &
   btor   = btor,   zimp   = zimp,   aimp   = aimp,       &
   ahyd   = ahyd,   aimass = aimass, wexbs  = wexbs,      &
   gne    = gne,    gni    = gni,    gnh    = gnh,        &
   gnz    = gnz,    gte    = gte,    gti    = gti,        &
   gq     = gq,                                           &
   gvtor  = gvtor,  vtor   = vtor,   gvpol  = gvpol,      &
   vpol   = vpol,   gvpar  = gvpar,  vpar   = vpar,       &
   xti    = xti,    xdi    = xdi,    xte    = xte,        &
   xdz    = xdz,    xvt    = xvt,    xvp    = xvp,        &
   xtiW20 = xtiW20, xdiW20 = xdiW20, xteW20 = xteW20,     &
   xtiDBM = xtiDBM, xdiDBM = xdiDBM, xteDBM = xteDBM,     &
   xteETG = xteETG,                                       &
   gammaW20 = gammaW20, omegaW20 = omegaW20,              &
   gammaDBM = gammaDBM, omegaDBM = omegaDBM,              &
   npoints = npoints,                                     &
   lprint  = lprint, nprout  = hfDebug, nerr    = nerr,   &
   vconv   = vconv,  vflux   = vflux ,                    &
   cmodel  = cmodel, cswitch = cswitch, lswitch = lswitch)

!.. Check for errors by mmm7_1
If ( nerr /= 0 ) Then
   write(*,'("MMM7.1 finished with an error "I3)') nerr
Else
   write(*,'("MMM7.1 finished sucessfully!")')
End If

!.. Write return values 
Write(hfOut,'("#.. Output profiles: electron thermal diffusivity&
   & and electron thermal flux")')
write(hfOut,'("#Unit"A6,23A11)') "m",              &
   "m^2/s","m^2/s","m^2/s","m^2/s","m^2/s","m^2/s",&
   "m^2/s","m^2/s","m^2/s",                        &
   "m^2/s","m^2/s","m^2/s","m^2/s",                &
   "s^-1","rad/s","s^-1","rad/s",                  &
   "s^-1","rad/s","s^-1","rad/s",                  &
   "s^-1","rad/s"
write(hfOut,'("# "A9,23A11)') "rmin",              &
   "xti", "xdi", "xte", "xdz", "xvt", "xvp",       &
   "xtiW20", "xdiW20", "xteW20",                   &
   "xtiDBM", "xdiDBM", "xteDBM", "xteETG",         &
   "gmaW20ii", "omgW20ii", "gmaW20ie", "omgW20ie", &
   "gmaW20ei", "omgW20ei", "gmaW20ee", "omgW20ee", &
   "gmaDBM",   "omgDBM"

Do jr=1, npoints
   write(hfOut,'(0PF11.6,23ES11.3)') rmin(jr), &
      xti(jr), xdi(jr), xte(jr),  &
      xdz(jr), xvt(jr), xvp(jr), &
      xtiW20(jr), xdiW20(jr), xteW20(jr),  &
      xtiDBM(jr), xdiDBM(jr), xteDBM(jr),  &
      xteETG(jr),                          &
      gammaW20(1,jr), omegaW20(1,jr), &
      gammaW20(2,jr), omegaW20(2,jr), &
      gammaW20(3,jr), omegaW20(3,jr), &
      gammaW20(4,jr), omegaW20(4,jr), &
      gammaDBM(jr),   omegaDBM(jr)
End Do

Close(unit=hfOut)
Close(unit=hfDebug)

Contains

Subroutine check_input( vval, vname ) !{{{
! This internal subroutine checks for variable with BADREAL
! If so, print an error message and increase the error count
Implicit None

Real(R8), Intent(In) :: vval ! the value to be checked
Character(len=*), Intent(In) :: vname ! the corresponding variable name
Character(2) :: vname_len ! length of the variable name

If ( abs( vval - dble(BADREAL) ) < 1e-6 ) Then
   Write(vname_len,'(I2)') len_trim(vname)
   Write(*,'("** ERROR ** """A'//trim(adjustl(vname_len))//&
      '""" is not found in input")') vname
   errcount = errcount + 1
End If

End Subroutine check_input !}}}

End Program testmmm
