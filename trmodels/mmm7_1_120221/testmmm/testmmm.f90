Program testmmm

!                 >>>>>   N O T E   <<<<<
! {{{ and }}} are VIM folding markers. DON NOT REMOVE THEM.

Use modmmm7_1

Implicit None

! SPECIFICATION {{{

!------- Constants -------------------------------------------------

Integer, Parameter :: R8 = KIND( 0D0 )

!.. Machine Epsilon
Real(R8), Parameter :: &
   zepslon = 2D0**(-53) ,&! As defined by IEEE 754-2008
   zepsqrt = 2D0**(-26)   ! Square root of Epsilon, approximated

Integer, Parameter :: &
   NPMAX  = 200, &! Maximum data points allowed in arrays
   hfIn   = 34,  &! Input file handle
   hfOut  = 35,  &! Output file handle
   hfDebug= 36    ! Diagnostic output

Integer, Parameter :: &
   MAXNOPT   = 10  ! Maximum number of internal switches

!------- Input Variables -------------------------------------------
! The content of the input variables follows their corresponding arguments
! of the mmm subroutine. See modmmm.f90 for more details.
Integer :: lprint, npoints, nerr

Real(R8), Dimension(4) :: cmodel

Real(R8), Dimension(MAXNOPT) :: &
   cW20=0D0, cDBM=0D0, cETG=0D0, cPLC=0D0 ! To be passed to cswitch

Integer, Dimension(MAXNOPT) :: &
   lW20=0, lDBM=0, lETG=0, lPLC=0 ! To be passed to lswitch

Real(R8), Dimension(NPMAX) :: &
   zrminor, zrmajor, zelong,                              &
   zdense, zdensh, zdensimp, zdensfe,                     &
   zxzeff, ztekev, ztikev, zq, zbtor, etanc,              &
   zavezimp, zamassimp, zamasshyd, zaimass, zwexbs,       &
   zgrdne, zgrdni, zgrdnh, zgrdnz, zgte, zgti, zshear,    &
   zgrdvphi, zvtorin, zgrdvtht, zvpolin, zgrdvpar, zvpar, &
   zchrgz, zchrgf, zxb, zdensf, zdensi, zgrdnf,           &
   zqprime

Real(R8) :: &
   rminor,   rmajor,   elong,    denmin,   temin,             &
   denhaxis, denhedge, denhexp,  amassh,   denzaxis, denzedge,&
   denzexp,  amassz,   chrzaxis, chrzedge, chrzexp,  denfaxis,&
   denfedge, denfexp,  chrfaxis, chrfedge, chrfexp,           &
   teaxis,   teedge,   teexp,    tiaxis,   tiedge,   tiexp,   &
   wexbmax,  xwexbinn, xwexbout,                              &
   qaxis,    qedge,    qexp,     btor, zdx, zxfact,           &
   etancis,  etancedge,etancexp, zwidth, rmajbnd

namelist /testmmm_input_1stkind/ &
   cmodel, npoints, lprint,                            &
   cW20, cDBM, cETG, cPLC,                             &
   lW20, lDBM, lETG, lPLC,                             &
   zrminor, zrmajor, zelong,                           &
   zdense, zdensh, zdensimp, zdensfe,                  &
   zxzeff, ztekev, ztikev, zq, zbtor, etanc,           &
   zavezimp, zamassimp, zamasshyd, zaimass, zwexbs,    &
   zgrdne, zgrdni, zgrdnh, zgrdnz, zgte, zgti, zshear, &
   zgrdvphi, zvtorin, zgrdvtht, zvpolin, zgrdvpar,     &
   zvpar, zqprime, zwidth, rmajbnd

namelist /testmmm_input_2ndkind/ &
   cmodel, npoints, lprint, &
   cW20, cDBM, cETG, cPLC,  &
   lW20, lDBM, lETG, lPLC,  &
   rminor,   rmajor,        &
   elong,    denmin,   temin,    denhaxis, denhedge, denhexp, &
   amassh,   denzaxis, denzedge, denzexp,  amassz,   denfaxis,&
   denfedge, denfexp,  chrzaxis, chrzedge, chrzexp,           &
   chrfaxis, chrfedge, chrfexp,  teaxis,   teedge,   teexp,   &
   tiaxis,   tiedge,   tiexp,    qaxis,    qedge,    qexp,    &
   wexbmax,  xwexbinn, xwexbout, btor, &
   etancis,  etancedge,etancexp

!------- Output Variables ------------------------------------------
! The content of the output variables follows their corresponding arguments
! of the mmm7_1 subroutine. See modmmm7_1.f90 for more details.

Real(R8), Dimension(NPMAX) :: &
   zthiig, zthdig, ztheig, zthzig, zthtig, zthttig, &
   xkiW20, xdiW20, xkeW20,                          &
   xkiDBM, xkhDBM, xkeDBM, gammaDBM, omegaDBM,      &
   xkeETG, xkePLC

Real(R8), Dimension(4,NPMAX) :: &
   gammaW20, omegaW20

Real(R8), Dimension(6,NPMAX) :: &
   zvelthi, zvflux

!------- Local Variables ------------------------------------------

Real(R8) :: cmmm(MAXNOPT,4)

Integer :: lmmm(MAXNOPT,4)

Integer :: jr

Integer :: input_kind

Character(len=255) :: strbuf
!}}}

!------- EXECUTION BODY --------------------------------------------

Open( hfIn, file='input', &
   form='formatted', status='old', iostat=nerr)
if ( nerr /= 0 ) stop "ERROR: Cannot open input file!"

Open( hfOut, file='output', &
   form='formatted', status='replace', iostat=nerr)
if ( nerr /= 0 ) stop "ERROR: Cannot creat output file!"

Open( hfDebug, file='testmmm_debug', &
   form='formatted', status='replace', iostat=nerr)
if ( nerr /= 0 ) stop "ERROR: Cannot creat debug file!"

Read( hfIn,'(A)') strbuf
Rewind( hfIn )

If ( index( strbuf, '1stkind') > 0 ) Then

! 1st kind input {{{
   input_kind = 1

   Write(*,'("Input of the first kind (values) is detected. Processing...")')

   Read( hfIn, NML=testmmm_input_1stkind )

   If ( zgte(1) < -1D2 ) Then
      zgte(1) = 0D0
      Do jr=2,npoints-1
         zgte(jr) = - zrmajor(jr)/ztekev(jr)* &
            ( ztekev(jr+1)-ztekev(jr-1) )/( zrminor(jr+1)-zrminor(jr-1) )
      End Do
      zgte(npoints) = - zrmajor(npoints)/ztekev(npoints)* &
         ( ztekev(npoints)-ztekev(npoints-1) )/&
         ( zrminor(npoints)-zrminor(npoints-1) )
   End If

   If ( zgti(1) < -1D2 ) Then
      zgti(1) = 0D0
      Do jr=2,npoints-1
         zgti(jr) = - zrmajor(jr)/ztikev(jr)* &
            ( ztikev(jr+1)-ztikev(jr-1) )/( zrminor(jr+1)-zrminor(jr-1) )
      End Do
      zgti(npoints) = - zrmajor(npoints)/ztikev(npoints)* &
         ( ztikev(npoints)-ztikev(npoints-1) )/&
         ( zrminor(npoints)-zrminor(npoints-1) )
   End If

   If ( zgrdne(1) < -1D2 ) Then
      zgrdne(1) = 0D0
      Do jr=2,npoints-1
         zgrdne(jr) = - zrmajor(jr)/zdense(jr)* &
            ( zdense(jr+1)-zdense(jr-1) )/( zrminor(jr+1)-zrminor(jr-1) )
      End Do
      zgrdne(npoints) = - zrmajor(npoints)/zdense(npoints)* &
         ( zdense(npoints)-zdense(npoints-1) )/&
         ( zrminor(npoints)-zrminor(npoints-1) )
   End If

   If ( zgrdnh(1) < -1D2 ) Then
      zgrdnh(1) = 0D0
      Do jr=2,npoints-1
         zgrdnh(jr) = - zrmajor(jr)/zdensh(jr)* &
            ( zdensh(jr+1)-zdensh(jr-1) )/( zrminor(jr+1)-zrminor(jr-1) )
      End Do
      zgrdnh(npoints) = - zrmajor(npoints)/zdensh(npoints)* &
         ( zdensh(npoints)-zdensh(npoints-1) )/&
         ( zrminor(npoints)-zrminor(npoints-1) )
   End If

   If ( zgrdnz(1) < -1D2 ) Then
      zgrdnz(1) = 0D0
      Do jr=2,npoints-1
         zgrdnz(jr) = - zrmajor(jr)/zdensimp(jr)* &
            ( zdensimp(jr+1)-zdensimp(jr-1) )/( zrminor(jr+1)-zrminor(jr-1) )
      End Do
      zgrdnz(npoints) = - zrmajor(npoints)/zdensimp(npoints)* &
         ( zdensimp(npoints)-zdensimp(npoints-1) )/&
         ( zrminor(npoints)-zrminor(npoints-1) )
   End If

   If ( zgrdni(1) < -1D2 ) Then
      zgrdni=zgrdnh+zgrdnz
   End If
!}}}

Elseif ( index( strbuf, '2ndkind') > 0 ) Then

! Input of the 2nd Kind{{{

   input_kind = 2

   Write(*,'("Input of the second kind (shape functions) is detected. Processing...")')

   read( hfIn, NML=testmmm_input_2ndkind )

   !..check to see if input makes sense
   If ( npoints < 5  .or.  npoints > NPMAX ) then
     write(*,'("Abort: npoints < 5  .or.  npoints > NPMAX ")')
     stop
   End If

   If ( denmin < 1.0 ) then
     write(hfOut,'("# denmin changed from "E14.6" to 1e7")') denmin
     denmin = 1.e7
   End If

   If ( temin < 1.0 ) then
     write(hfOut,'("# tenmin changed from "E14.6" to 1e7")') temin
     temin = 1.e-6
   End If

   !..set up argument list for subroutine theory

   !..equilibrium on zone boundaries
   zdx = 1D0 / (npoints-1)
   Do jr=1,npoints
      zxb(jr)     = ( jr - 1 ) * zdx
      zrminor(jr) = ( jr - 1 ) * zdx * rminor
      zrmajor(jr) = rmajor + zrminor(jr)
      zelong(jr)  = elong
   End Do

   zwidth = zdx * rminor
   rmajbnd = zrmajor(npoints)

   !..construct ion density profiles on zone boundaries
   Do jr=1,npoints-1
      zdensh(jr) = denhedge + (denhaxis-denhedge) * (1.0 - zxb(jr)**2)**denhexp
      zdensimp(jr) = denzedge + (denzaxis-denzedge) * (1.0 - zxb(jr)**2)**denzexp
      zdensf(jr) = denfedge + (denfaxis-denfedge) * (1.0 - zxb(jr)**2)**denfexp
      ztekev(jr) = max ( temin, teedge + (teaxis-teedge) * (1.0 - zxb(jr)**2)**teexp )
      ztikev(jr) = max ( temin, tiedge + (tiaxis-tiedge) * (1.0 - zxb(jr)**2)**tiexp )
      zchrgz(jr) = max (1.0, chrzedge + (chrzaxis-chrzedge) * (1.0 - zxb(jr)**2)**chrzexp )
      zchrgf(jr) = max (1.0, chrfedge + (chrfaxis-chrfedge) * (1.0 - zxb(jr)**2)**chrfexp )
      etanc(jr)  = (etancedge-etancis)* zxb(jr)**2**etancexp
   End Do

   zdensh(npoints) = denhedge
   zdensimp(npoints) = denzedge
   zdensf(npoints) = denfedge
   ztekev(npoints) = max ( temin, teedge )
   ztikev(npoints) = max ( temin, tiedge )
   zchrgz(npoints) = max ( 1.0, chrzedge )
   zchrgf(npoints) = max ( 1.0, chrfedge )
   etanc(npoints) = etancedge

   !..electron and thermal ion densities
   Do jr=1,npoints
     zdense(jr) = max ( denmin,  zdensh(jr) &
        + zchrgz(jr) * zdensimp(jr) + zchrgf(jr) * zdensf(jr) )
     zdensi(jr) = max ( denmin, zdensh(jr) + zdensimp(jr) )
     zdensfe(jr) = zchrgf(jr) * zdensf(jr)
   End Do

   !..normalized gradients
   Do jr=1,npoints

     zxfact = max ( zepslon, 1.0 - zxb(jr)**2 )

     zgrdnh(jr) = zrmajor(jr) * 2.0 * zxb(jr) * denhexp * &
        (denhaxis-denhedge) * zxfact**(denhexp-1.0) &
        / ( max( denmin, zdensh(jr) ) * zrminor(npoints) )

     zgrdnz(jr) = zrmajor(jr) * 2.0 * zxb(jr) * ( &
        (denzaxis-denzedge) * zxfact**(denzexp-1.0) * &
          denzexp / ( max( denmin, zdensimp(jr) ) * zrminor(npoints) ) &
        + (chrzaxis-chrzedge) * zxfact**(chrzexp-1.0) * &
          chrzexp / ( zchrgz(jr) * zrminor(npoints) ) )

     zgrdnf(jr) = zrmajor(jr) * 2.0 * zxb(jr) * ( &
        (denfaxis-denfedge) * zxfact**(denfexp-1.0) * &
          denfexp / ( max( denmin, zdensf(jr) ) * zrminor(npoints) ) &
        + (chrfaxis-chrfedge) * zxfact**(chrfexp-1.0) * &
          chrfexp / ( zchrgf(jr) * zrminor(npoints) ) )

     zgte(jr) = zrmajor(jr) * 2.0 * zxb(jr) * teexp * &
        (teaxis-teedge) * zxfact**(teexp-1.0) &
          / ( ztekev(jr) * zrminor(npoints) )

     zgti(jr) = zrmajor(jr) * 2.0 * zxb(jr) * tiexp * &
        (tiaxis-tiedge) * zxfact**(tiexp-1.0) &
          / ( ztikev(jr) * zrminor(npoints) )

     zgrdni(jr) = 2.0 * zxb(jr) * zrmajor(jr) * ( &
        + denhexp*(denhaxis-denhedge) * zxfact**(denhexp-1.0) &
        + denzexp*(denzaxis-denzedge) * zxfact**(denzexp-1.0) &
        ) / ( zrminor(npoints) * zdensi(jr) )

   End Do

   !..normalized electron density gradient
   Do jr=1,npoints
     zgrdne(jr) = ( zdensh(jr) * zgrdnh(jr) &
        + zchrgz(jr) * zdensimp(jr)*zgrdnz(jr) &
        + zchrgf(jr) * zdensf(jr) * zgrdnf(jr)  ) / zdense(jr)
   End Do

   !..average mass and charge
   do jr=1,npoints
      zavezimp(jr) = zchrgz(jr)
      zamassimp(jr) = amassz
      zamasshyd(jr) = amassh
      zaimass(jr) = &
         ( zamasshyd(jr) * zdensh(jr) + zamassimp(jr) * zdensimp(jr) ) &
         / ( zdensh(jr) + zdensimp(jr) )
      zxzeff(jr) =  &
         ( zdensh(jr) + zchrgz(jr)**2 * zdensimp(jr) &
         + zchrgf(jr)**2 * zdensf(jr) ) &
         / zdense(jr)
   End Do

   !..simple magnetics profiles
   do jr=1,npoints
     zq(jr) = qaxis + (qedge-qaxis) * abs(zxb(jr))**qexp
     zshear(jr) = zrmajor(jr) * &
        (qedge-qaxis) * (abs(zxb(jr)))**(qexp-1.0) &
          / ( zq(jr) * zrminor(npoints) )
     zqprime = (qedge-qaxis) * qexp * abs(zxb(jr))**(qexp-1D0)
   End Do

   do jr=1,npoints
     zbtor(jr)  = btor
   End Do

   !   Note:  The flow shear rate is given as a fourth-order
   !     polynomial with maximum value wexbmax between
   !     xwexbinn < r/a < xwexbout.
   !
   ! wexbmax    maximum flow shear rate [sec^-1]
   ! xwexbinn   inner r/a cutoff for the flow shear rate
   ! xwexbout   outer r/a cutoff for the flow shear rate
   do jr=1,npoints
     zwexbs(jr) = 0.0
     if ( zxb(jr) > xwexbinn .and.  zxb(jr) < xwexbout ) then
       zwexbs(jr) = wexbmax * 16.0 &
          * ( ( zxb(jr) - xwexbinn ) / ( xwexbout - xwexbinn ) )**2 &
          * ( ( xwexbout - zxb(jr) ) / ( xwexbout - xwexbinn ) )**2
     End If
   End Do
   
   zgrdvphi = 0D0
   zvtorin  = 0D0
   zgrdvtht = 0D0
   zvpolin  = 0D0
   zgrdvpar = 0D0
   zvpar    = 0D0

!}}}

Else ! Unsupported kind

   write(*,'("Unsupported input! Stop.")')
   Close( unit=hfIn )
   stop

End If

Close( unit=hfIn )

cmmm(1:MAXNOPT,1)=cW20
cmmm(1:MAXNOPT,2)=cDBM
cmmm(1:MAXNOPT,3)=cETG

lmmm(1:MAXNOPT,1)=lW20
lmmm(1:MAXNOPT,2)=lDBM
lmmm(1:MAXNOPT,3)=lETG

!.. Write profiles
write(hfOut,'("#I"I6,31I11)') (jr,jr=1,31)
write(hfOut,'("#I"A9,31A11)') "rmin","rmaj", &
   "elong","dense", "densh","densimp","densfe","xzeff", &
   "tekev","tikev", "q","btor","etanc","avezimp",       &
   "amassimp","amasshyd","aimass","wexbs","grdne",      &
   "grdni","grdnh","grdnz","gte","gti","shear",         &
   "grdvphi","vtorin","grdvtht","vpolin","grdvpar","vpar"

Do jr=1, npoints
   write(hfOut,'(0P2F11.6,1P29E11.3)') &
      zrminor(jr), zrmajor(jr),   &
      zelong(jr), zdense(jr),   &
      zdensh(jr), zdensimp(jr), zdensfe(jr),     &
      zxzeff(jr), ztekev(jr), ztikev(jr),        &
      zq(jr), zbtor(jr), etanc(jr), zavezimp(jr),     &
      zamassimp(jr), zamasshyd(jr), zaimass(jr), &
      zwexbs(jr), zgrdne(jr), zgrdni(jr), zgrdnh(jr), &
      zgrdnz(jr), zgte(jr), zgti(jr), zshear(jr),     &
      zgrdvphi(jr), zvtorin(jr), zgrdvtht(jr),   &
      zvpolin(jr), zgrdvpar(jr), zvpar(jr)
End Do

zthiig = 0D0
zthdig = 0D0
ztheig = 0D0

If ( input_kind == 1 ) Then

   call mmm7_1( &
      rmin  = zrminor,  rmaj  = zrmajor,  elong = zelong,      &
      ne    = zdense,   ni    = zdensh,   nz    = zdensimp,    &
      nf    = zdensfe,  xzeff = zxzeff,   te    = ztekev,      &
      ti    = ztikev,   q     = zq,       btor  = zbtor,       &
      zimp  = zavezimp, aimp  = zamassimp, ahyd = zamasshyd,   &
      aimass= zaimass,  wexbs = zwexbs,                        &
      gne   = zgrdne,   gni   = zgrdni,    gnh  = zgrdnh,      &
      gnz   = zgrdnz,   gte   = zgte,      gti  = zgti,        &
      gq    = zshear,                                          &
      gvrin  = zgrdvphi, vtorin = zvtorin,  gvpin  = zgrdvtht, &
      vpolin = zvpolin, gvparin = zgrdvpar, vparin = zvpar,    &
      thiig  = zthiig, thdig  = zthdig, theig  = ztheig,       &
      thzig  = zthzig, thtig  = zthtig, thttig = zthttig,      &
      xkiW20 = xkiW20, xdiW20 = xdiW20, xkeW20 = xkeW20,       &
      xkiDBM = xkiDBM, xkhDBM = xkhDBM, xkeDBM = xkeDBM,       &
      xkeETG = xkeETG, &
      gammaW20 = gammaW20, omegaW20 = omegaW20,                &
      gammaDBM = gammaDBM, omegaDBM = omegaDBM,                &
      npoints= npoints,                                        &
      lprint = lprint, nprout = hfDebug, nerr = nerr,          &
      velthi = zvelthi, vflux = zvflux ,                       &
      cmodel = cmodel, cswitch = cmmm, lswitch = lmmm)

Else

   !print *,"csnd0=",1D-4
   call mmm7_1( &
      rmin  = zrminor,  rmaj  = zrmajor,  elong = zelong,      &
      ne    = zdense,   ni    = zdensh,   nz    = zdensimp,    &
      nf    = zdensfe,  xzeff = zxzeff,   te    = ztekev,      &
      ti    = ztikev,   q     = zq,       btor  = zbtor,       &
      zimp  = zavezimp, aimp  = zamassimp, ahyd = zamasshyd,   &
      aimass= zaimass,  wexbs = zwexbs,                        &
      gne   = zgrdne,   gni   = zgrdni,    gnh  = zgrdnh,      &
      gnz   = zgrdnz,   gte   = zgte,      gti  = zgti,        &
      gq    = zshear,                                          &
      gvrin  = zgrdvphi, vtorin = zvtorin,  gvpin  = zgrdvtht, &
      vpolin = zvpolin, gvparin = zgrdvpar, vparin = zvpar,    &
      thiig  = zthiig, thdig  = zthdig, theig  = ztheig,       &
      thzig  = zthzig, thtig  = zthtig, thttig = zthttig,      &
      xkiW20 = xkiW20, xdiW20 = xdiW20, xkeW20 = xkeW20,       &
      xkiDBM = xkiDBM, xkhDBM = xkhDBM, xkeDBM = xkeDBM,       &
      xkeETG = xkeETG, &
      gammaW20 = gammaW20, omegaW20 = omegaW20,                &
      gammaDBM = gammaDBM, omegaDBM = omegaDBM,                &
      npoints = npoints,                                       &
      lprint = lprint, nprout = hfDebug, nerr=nerr,            &
      velthi = zvelthi, vflux = zvflux ,                       &
      cmodel = cmodel, cswitch = cmmm, lswitch = lmmm)
End If

If ( nerr /= 0 ) Then
   write(*,'("MMM finished with an error "I3)') nerr
Else
   write(*,'("MMM finished sucessfully!")')
End If

!.. Write return values 
write(hfOut,'("#O"A9,23A11)') "rmin",&
   "thiig",  "thdig",  "theig",  &
   "xkiW20", "xdiW20", "xkeW20", &
   "xkiDBM", "xkhDBM", "xkeDBM", &
   "xkeETG",                     &
   "thzig",  "thtig",  "thttig", &
   "gmaW20ii", "omgW20ii", "gmaW20ie", "omgW20ie", &
   "gmaW20ei", "omgW20ei", "gmaW20ee", "omgW20ee", &
   "gmaDBM",   "omgDBM"

Do jr=1, npoints
   write(hfOut,'(0PF11.6,23ES11.3)') zrminor(jr), &
      zthiig(jr), zthdig(jr), ztheig(jr),  &
      xkiW20(jr), xdiW20(jr), xkeW20(jr),  &
      xkiDBM(jr), xkhDBM(jr), xkeDBM(jr),  &
      xkeETG(jr),                          &
      zthzig(jr), zthtig(jr), zthttig(jr), &
      gammaW20(1,jr), omegaW20(1,jr), &
      gammaW20(2,jr), omegaW20(2,jr), &
      gammaW20(3,jr), omegaW20(3,jr), &
      gammaW20(4,jr), omegaW20(4,jr), &
      gammaDBM(jr),   omegaDBM(jr)
End Do

Close(unit=hfOut)
Close(unit=hfDebug)

End Program testmmm
