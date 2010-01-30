!| %
!| % test_mmm95.f  Glenn Bateman  Lehigh University
!| %---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!| \documentstyle{article}
!| \headheight 0pt \headsep 0pt  \topmargin 0pt  \textheight 9.0in
!| \oddsidemargin 0pt \textwidth 6.5in
!| \begin{document}
!| 
!| \begin{center} 
!| {\bf {\tt test-mmm95.tex} \\
!| Stand-alone Code to Test Subroutine MMM95 } \\
!| \vspace{1pc}  Glenn Bateman and Arnold Kritz \\
!| Lehigh University \\
!| Version 1.0$\quad$ 9 October 1998
!| \end{center}
!| 
!| This is a stand-alone program to test subroutine mmm95.
!| 
!| For the purposes of this test program, only three ion species are
!| considered: a thermal hydrogen species, a thermal impurity species,
!| and a fast ion (super-thermal) ion species.
!| 
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c-----------------------------------------------------------------------
c
c  Compile this routine and routines that it calls with a compiler 
c  option, such as -r8, to convert real to double precision when used on 
c  workstations.
c
c-----------------------------------------------------------------------
c
c  External dependencies:
c
c  Call tree: TEST_MMM95 calls the following routines
c
c  MMM95	     - Computes Multi-Mode Model transport coefficients
c    WEILAND14       - Computes diffusion matrix and convect velocities
c      WEILAND14FLUX - Calculates fluxes and effective diffusivities
c        TOMSQZ      - Wrapper for QZ algorithm solving Ax = lambda Bx
c          CQZHES    - First step in QZ algorithm 
c          CQZVAL    - Second and third step in QZ algorithm
c          CQZVEC    - Fourth step in QZ algorithm
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer      kr,       km,       klswitch, kcswitch
      parameter ( kr = 55, km = 12, klswitch = 8, kcswitch = 25 )
c
c  kr = max number of radial points
c  km = first dimension of diffusivity matrix
c  klswitch = number of elements in lswitch array
c  kcswitch = number of elements in cswitch array
c
      integer lswitch(8)
     & , nin,      nout,     ntemp,    lprint,   maxis,     medge
     & , matdim,   npoints,  lsuper,   lreset,   nerr
c
      real cswitch(25)
     & , rminor,   rmajor,   elong,    denmin,   temin
     & , denhaxis, denhedge, denhexp,  amassh,   denzaxis, denzedge
     & , denzexp,  amassz,   chrzaxis, chrzedge, chrzexp,  denfaxis
     & , denfedge, denfexp,  chrfaxis, chrfedge, chrfexp
     & , teaxis,   teedge,   teexp,    tiaxis,   tiedge,   tiexp
     & , wexbmax,  xwexbinn, xwexbout
     & , qaxis,    qedge,    qexp,     btor,     fig(4),   fkb(4)   
     & , frb(4)
c
c..Namelist input:
c
c  lswitch(jc), j=1,8    integer control variables
c    (see documentation in file mmm95.tex)
c  cswitch(jc), jc=1,25  general control variables
c    (see documentation in file mmm95.tex)
c
c  maxis   = zone boundary at magnetic axis
c  medge   = zone boundary at edge of plasma
c
c  rminor = minor radius (half-width) at edge of plasma
c  rmajor = major radius to geometric center at edge of plasma
c  elong  = elongation at edge of plasma
c
c  denmin = minimum value allowed for n_e
c  temin  = minimum value allowed for T_e, T_i, and T_f
c
c    In this test program, the profiles are all taken to be
c  parabolas raised to powers:   For example
c  T_i = (tiaxis-tiedge)*(1.-x^2)^tiexp + tiedge
c  where x = sqrt normalized toroidal flux at zone boundaries
c
c  The namelist input parameters are:
c
c maxis      index of magnetic axis
c medge      index of edge of plasma
c
c rminor     minor radius [m]
c rmajor     major radius [m]
c elong      elongation
c
c denmin     minimum density [keV]
c temin      minimum temperature  [keV]
c
c denhaxis   hydrogen density at magnetic axis [m^-3]
c denhedge   hydrogen density at plasma edge [m^-3]
c denhexp    exponential of parabola for profile
c amassh     isotopic mass of hydrogenic species
c
c denzaxis   impurity density at magnetic axis [m^-3]
c denzedge   impurity density at plasma edge [m^-3]
c denzexp    exponential of parabola for profile
c amassz     isotopic mass of impurity species
c
c denfaxis   fast ion density at magnetic axis [m^-3]
c denfedge   fast ion density at plasma edge [m^-3]
c denfexp    exponential of parabola for profile
c
c chrzaxis   impurity charge at magnetic axis
c chrzedge   impurity charge at plasma edge
c chrzexp    exponential of parabola for profile
c chrfaxis   fast ion charge at magnetic axis
c chrfedge   fast ion charge at plasma edge
c chrfexp    exponential of parabola for profile
c
c teaxis     electron temperature at magnetic axis [keV]
c teedge     electron temperature at plasma edge [keV]
c teexp      exponential of parabola for profile
c
c tiaxis     ion temperature at magnetic axis [keV]
c tiedge     ion temperature at plasma edge [keV]
c tiexp      exponential of parabola for profile
c
c qaxis      magnetic q-value at magnetic axis
c qedge      magnetic q-value at plasma edge
c qexp       q(x) = qaxis + (qedge-qaxis) * x**qexp
c
c btor       toroidal field at rmajor [Tesla]
c
c   Note:  The flow shear rate is given as a fourth-order
c     polynomial with maximum value wexbmax between
c     xwexbinn < r/a < xwexbout.
c
c wexbmax    maximum flow shear rate [sec^-1]
c xwexbinn   inner r/a cutoff for the flow shear rate
c xwexbout   outer r/a cutoff for the flow shear rate
c 
c
c
c lprint     controls amount of printout
c lsuper     used for supershots
c lreset     use default values for cswitch and lswitch
c nerr       error if nerr returns not equal to 0
c
c    The magnetic q profile is given by
c  q(x) = qaxis + ( qedge - qaxis ) * x**(qexp)
c
c  fig(j), ..., fkb(j), j=1,4
c    coefficients for the various modes of turbulence
c    (see documentation in file mmm95.tex)
c
c
      integer        jr
c
      real           zdx,          zxfact,       zepslon
c
      real 
     &   zxb(kr),    zrminor(kr),  zrmajor(kr),  zelong(kr), zdense(kr)
     & , zdensi(kr), zdensh(kr),   zdensz(kr),   zdensf(kr), zdensfe(kr)
     & , zxzeff(kr), ztekev(kr),   ztikev(kr),   zq(kr),     zbtor(kr)
     & , zwexb(kr),  zavezimp(kr), zmassimp(kr), zmasshyd(kr)
     & , zchrgz(kr), zchrgf(kr),   zaimass(kr)
c
c
c  matdim  = first and second dimension of transport matricies
c              difthi(j1,j2,jz) and velthi(j1,jz)
c
c    All the following 1-D arrays are assumed to be on zone boundaries
c    including the densities (in m^-3) and temperatures (in keV).
c
c  zrminor(jz) = minor radius (half-width) of zone boundary [m]
c  zrmajor(jz) = major radius to geometric center of zone bndry [m]
c  zelong(jz)  = local elongation of zone boundary
c
c  zaimass(jz) = mean atomic mass of thermal ions [AMU]
c  zdense(jz)  = electron density [m^-3]
c  zdensi(jz)  = thermal ion density [m^-3]
c  zdensh(jz)  = thermal hydrogen ion density [m^-3]
c  zdensz(jz)  = sum over impurity ion densities [m^-3]
c  zdensf(jz)  = fast (non-thermal) ion density [m^-3]
c  zdensfe(jz) = electron density from fast (non-thermal) ions [m^-3]
c  zxzeff(jz)  = Z_eff
c  ztekev(jz)  = T_e [keV]  (electron temperature)
c  ztikev(jz)  = T_i [keV]  (temperature of thermal ions)
c  zq(jz)      = magnetic q-value
c  zbtor(jz)   = ( R B_tor ) / rmajor(jz)  [tesla]
c  zresist(jz) = plasma resistivity [Ohm-m]
c
c  zavezimp(jz) = average density weighted charge of impurities
c  zmassimp(jz) = average density weighted atomic mass of impurities
c  zmasshyd(jz) = average density weighted atomic mass of hydrogen ions
c
c  zchrgz(jz)   = ionic charge of the thermal impurity ions
c  zchrgf(jz)   = ionic charge of the super-thermal ions
c
c..normalized gradients:
c
      real 
     &   zgrdne(kr), zgrdni(kr),   zgrdnh(kr),   zgrdnz(kr), zgrdnf(kr)
     & , zgrdte(kr), zgrdti(kr),   zgrdq(kr)
c
c    All of the following normalized gradients are at zone boundaries.
c    r = half-width, R = major radius to center of flux surface
c    carry out any smoothing before calling sbrtn theory
c
c  zgrdne(jz) = - R ( d n_e / d r ) / n_e
c  zgrdni(jz) = - R ( d n_i / d r ) / n_i
c     n_i = thermal ion density (sum over hydrogenic and impurity)
c  zgrdnh(jz) = - R ( d n_h / d r ) / n_h
c     n_h = thermal hydrogenic density (sum over hydrogenic species)
c  zgrdnz(jz) = - R ( d Z n_Z / d r ) / ( Z n_Z )
c     n_Z = thermal impurity density,  Z = average impurity charge
c           sumed over all impurities
c  zgrdnf(jz) = - R ( d Z_f n_f / d r ) / ( Z_f n_f )
c     n_f = superthermal ion density,
c     Z_f = average superthermal ion charge at gridpoint jz
c  zgrdte(jz) = - R ( d T_e / d r ) / T_e
c  zgrdti(jz) = - R ( d T_i / d r ) / T_i
c  zgrdq (jz) =   R ( d q   / d r ) / q    related to magnetic shear
c
c..transport output:
c
      real 
     &   zdifthi(km,km,kr),       zvelthi(km,kr),zflux(km,kr) 
     & , zgamma(km,kr),           zomega(km,kr), zthekb(kr), zthzkb(kr)
     & , zthiig(kr), zthdig(kr),  ztheig(kr),    zthzig(kr), zthirb(kr)
     & , zthdrb(kr), ztherb(kr),  zthzrb(kr),    zthikb(kr), zthdkb(kr)
     & 
c
c  zdifthi(j1,j2,jz) = full matrix of anomalous transport diffusivities
c                     [m^2/sec]
c  zvelthi(j1,jz)    = convective velocities [m/sec]
c  zflux(j1,jz)      = total vflux [m/sec]
c
c    j1 = 1  for ion thermal transport
c    j1 = 2  for Hydrogenic particle transport
c    j1 = 3  for electron thermal transport
c    j1 = 4  for impurity particle transport
c
c  zgamma(jm,jr) = growth rate for mode jm [1/sec]
c  zomega(jm,jr) = frequency for mode jm [radians/sec]
c
c  zthiig(jr) = ion thermal diffusivity from the Weiland model
c
c
c..namelist input:
c
      namelist /nt/  
     &   lswitch,  cswitch,  maxis,    medge,    rminor,   rmajor
     & , elong,    denmin,   temin,    denhaxis, denhedge, denhexp
     & , amassh,   denzaxis, denzedge, denzexp,  amassz,   denfaxis
     & , denfedge, denfexp,  chrzaxis, chrzedge, chrzexp
     & , chrfaxis, chrfedge, chrfexp,  teaxis,   teedge,   teexp
     & , tiaxis,   tiedge,   tiexp,    qaxis,    qedge,    qexp
     & , wexbmax,  xwexbinn, xwexbout
     & , btor,     fig,      fkb,      frb,      lprint,   lsuper
     & , lreset
c
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c..open input and output files
c  and a temporary file for the namelist with comments stripped off
c
      open (2,file='temp')
      open (5,file='input')
      open (6,file='output')
c
c..defaults
c
      ntemp   =  2
      nin     =  5
      nout    =  6
c
      matdim  = km
c
      nerr    =  0
c
c  for subroutine mmm95
c
      do jr=1,kcswitch
        cswitch(jr) = 0.0
      enddo
c
      do jr=1,klswitch
        lswitch(jr) = 0
      enddo
c
c..defaults for the contributions to the fluxes and interchange
c
      do jr=1,4
        fig(jr)=0.0
        fkb(jr)=0.0
        frb(jr)=0.0
      enddo
c..default input values
c
      maxis   = 0         ! index of magnetic axis
      medge   = 20        ! index of edge of plasma

      rminor  = 0.9       ! minor radius [m]
      rmajor  = 2.5       ! major radius [m]
      elong   = 1.0       ! elongation
      
      denmin  = 0.0       ! minimum density = 1.e7 [keV]
      temin   = 0.0       ! minimum temperature = 1.e-6 [keV]

      denhaxis = 0.3e19   ! hydrogen density at magnetic axis [m^-3]
      denhedge = 1.0e18   ! hydrogen density at plasma edge [m^-3]
      denhexp  = 1.0      ! exponential of parabola for profile [m^-3]
      amassh   = 2.0      ! isotopic mass of hydrogenic species

      denzaxis = 0.7e18   ! impurity density at magnetic axis [m^-3]
      denzedge = 3.0e17   ! impurity density at plasma edge [m^-3]
      denzexp  = 1.0      ! exponential of parabola for profile [m^-3]
      amassz   = 12.0     ! isotopic mass of impurity species
      
      denfaxis = 0.0e17   ! fast ion density at magnetic axis [m^-3]
      denfedge = 0.0e16   ! fast ion density at plasma edge [m^-3]
      denfexp  = 1.0      ! exponential of parabola for profile
      
      chrzaxis = 6.0      ! impurity charge at magnetic axis
      chrzedge = 6.0      ! impurity charge at plasma edge
      chrzexp  = 1.0      ! exponential of parabola for profile
      chrfaxis = 1.0      ! fast ion charge at magnetic axis
      chrfedge = 1.0      ! fast ion charge at plasma edge
      chrfexp  = 1.0      ! exponential of parabola for profile
      
      teaxis   = 2.0      ! electron temperature at magnetic axis [keV]
      teedge   = 0.2      ! electron temperature at plasma edge [keV]
      teexp    = 1.0      ! exponential of parabola for profile
      
      tiaxis   = 2.0      ! ion temperature at magnetic axis [keV]
      tiedge   = 0.2      ! ion temperature at plasma edge [keV]
      tiexp    = 1.0      ! exponential of parabola for profile
      
      qaxis    = 0.8      ! magnetic q-value at magnetic axis
      qedge    = 4.4      ! magnetic q-value at plasma edge
      qexp     = 4.0      ! q(x) = qaxis + (qedge-qaxis) * x**qexp
      
      btor     = 5.0      ! toroidal field at rmajor [Tesla]

      wexbmax  = 0.0      ! maximum flow shear rate [sec^-1]
      xwexbinn = 0.6      ! inner r/a cutoff for the flow shear rate
      xwexbout = 0.8      ! outer r/a cutoff for the flow shear rate
      
      lprint   = 0        ! controls amount of printout
      lsuper   = 0        ! used for supershots
      lreset   = 0        ! use default values for cswitch and lswitch
c
c..strip comments off namelist input
c
      call stripx ( nin, ntemp, nout )
c
c..read namelist
c
      read (ntemp, nt)
c
c
c..check to see if input makes sense
c
      if ( medge .lt. 5  .or.  medge .gt. kr ) then
        write (nout,*) 'Abort: medge .lt. 5  .or.  medge .gt. kr'
        stop
      endif
c
      if ( denmin .lt. 1.0 ) then
        write (nout,*) 'denmin changed from ',denmin,' to 1.e7'
        denmin = 1.e7
      endif
c
      if ( temin .lt. 1.0 ) then
        write (nout,*) 'temin changed from ',temin,' to 1.e-6'
        temin = 1.e-6
      endif
c
c
c..set up argument list for subroutine theory
c

c
c..equilibrium on zone boundaries
c
      zdx = 1.0 / real (medge - maxis)
      do jr=1,medge
        zxb(jr)     = ( jr - maxis ) * zdx
        zrminor(jr) = ( jr - maxis ) * zdx * rminor
        zrmajor(jr) = rmajor
        zelong(jr)  = elong
      enddo
c
c..construct ion density profiles on zone boundaries
c
      do jr=1,medge-1
        zdensh(jr) = denhedge
     &    + (denhaxis-denhedge) * (1.0 - zxb(jr)**2)**denhexp
        zdensz(jr) = denzedge
     &    + (denzaxis-denzedge) * (1.0 - zxb(jr)**2)**denzexp
        zdensf(jr) = denfedge
     &    + (denfaxis-denfedge) * (1.0 - zxb(jr)**2)**denfexp
        ztekev(jr) = max ( temin, teedge
     &    + (teaxis-teedge) * (1.0 - zxb(jr)**2)**teexp )
        ztikev(jr) = max ( temin, tiedge
     &    + (tiaxis-tiedge) * (1.0 - zxb(jr)**2)**tiexp )
        zchrgz(jr) = max (1.0, chrzedge
     &    + (chrzaxis-chrzedge) * (1.0 - zxb(jr)**2)**chrzexp )
        zchrgf(jr) = max (1.0, chrfedge
     &    + (chrfaxis-chrfedge) * (1.0 - zxb(jr)**2)**chrfexp )
      enddo
c
        zdensh(medge) = denhedge
        zdensz(medge) = denzedge
        zdensf(medge) = denfedge
        ztekev(medge) = max ( temin, teedge )
        ztikev(medge) = max ( temin, tiedge )
        zchrgz(medge) = max ( 1.0, chrzedge )
        zchrgf(medge) = max ( 1.0, chrfedge )
c
c..electron and thermal ion densities
c
      do jr=1,medge
        zdense(jr) = max ( denmin,  zdensh(jr)
     &    + zchrgz(jr) * zdensz(jr) + zchrgf(jr) * zdensf(jr) )
        zdensi(jr) = max ( denmin, zdensh(jr) + zdensz(jr) )
        zdensfe(jr) = zchrgf(jr) * zdensf(jr)
      enddo
c
c..normalized gradients
c
      zepslon = 1.e-9
c
      do jr=1,medge
c
        zxfact = max ( zepslon, 1.0 - zxb(jr)**2 )
c
        zgrdnh(jr) = zrmajor(jr) * 2.0 * zxb(jr) * denhexp *
     &    (denhaxis-denhedge) * zxfact**(denhexp-1.0)
     &      / ( max( denmin, zdensh(jr) ) * zrminor(medge) )
c
        zgrdnz(jr) = zrmajor(jr) * 2.0 * zxb(jr) * (
     &    (denzaxis-denzedge) * zxfact**(denzexp-1.0) *
     &      denzexp / ( max( denmin, zdensz(jr) ) * zrminor(medge) )
     &    + (chrzaxis-chrzedge) * zxfact**(chrzexp-1.0) *
     &      chrzexp / ( zchrgz(jr) * zrminor(medge) ) )
c
        zgrdnf(jr) = zrmajor(jr) * 2.0 * zxb(jr) * (
     &    (denfaxis-denfedge) * zxfact**(denfexp-1.0) *
     &      denfexp / ( max( denmin, zdensf(jr) ) * zrminor(medge) )
     &    + (chrfaxis-chrfedge) * zxfact**(chrfexp-1.0) *
     &      chrfexp / ( zchrgf(jr) * zrminor(medge) ) )
c
        zgrdte(jr) = zrmajor(jr) * 2.0 * zxb(jr) * teexp *
     &    (teaxis-teedge) * zxfact**(teexp-1.0)
     &      / ( ztekev(jr) * zrminor(medge) )
c
        zgrdti(jr) = zrmajor(jr) * 2.0 * zxb(jr) * tiexp *
     &    (tiaxis-tiedge) * zxfact**(tiexp-1.0)
     &      / ( ztikev(jr) * zrminor(medge) )
c
        zgrdni(jr) = 2.0 * zxb(jr) * zrmajor(jr) * (
     &    + denhexp*(denhaxis-denhedge) * zxfact**(denhexp-1.0)
     &    + denzexp*(denzaxis-denzedge) * zxfact**(denzexp-1.0)
     &    ) / ( zrminor(medge) * zdensi(jr) )
c
      enddo
c
c..normalized electron density gradient
c
      do jr=1,medge
        zgrdne(jr) = (   zdensh(jr) * zgrdnh(jr) 
     &    + zchrgz(jr) * zdensz(jr) * zgrdnz(jr)
     &    + zchrgf(jr) * zdensf(jr) * zgrdnf(jr) ) / zdense(jr)
      enddo
c
c..average mass and charge
c
      do jr=1,medge
        zavezimp(jr) = zchrgz(jr)
        zmassimp(jr) = amassz
        zmasshyd(jr) = amassh
        zaimass(jr) = 
     &    ( zmasshyd(jr) * zdensh(jr) + zmassimp(jr) * zdensz(jr) )
     &    / ( zdensh(jr) + zdensz(jr) )
        zxzeff(jr) = 
     &    ( zdensh(jr) + zchrgz(jr)**2 * zdensz(jr)
     &    + zchrgf(jr)**2 * zdensf(jr) )
     &    / zdense(jr)
       enddo
c
c
c..simple magnetics profiles
c
      do jr=maxis+1,medge
        zq(jr) = qaxis + (qedge-qaxis) * abs(zxb(jr))**qexp
        zgrdq(jr) = zrmajor(jr) *
     &    (qedge-qaxis) * (abs(zxb(jr)))**(qexp-1.0)
     &      / ( zq(jr) * zrminor(medge) )
      enddo
c
c
      if ( maxis .gt. 0 ) then
        zq(maxis) = qaxis
        zgrdq(maxis) = 0.0
      endif
      if ( maxis .gt. 1 ) then
        zq(maxis-1) = zq(maxis+1)
        zgrdq(maxis-1) = - zgrdq(maxis+1)
      endif
c
      do jr=1,medge
        zbtor(jr)  = btor
      enddo
c
c   Note:  The flow shear rate is given as a fourth-order
c     polynomial with maximum value wexbmax between
c     xwexbinn < r/a < xwexbout.
c
c wexbmax    maximum flow shear rate [sec^-1]
c xwexbinn   inner r/a cutoff for the flow shear rate
c xwexbout   outer r/a cutoff for the flow shear rate
c
      do jr=1,medge
        zwexb(jr) = 0.0
        if ( zxb(jr) .gt. xwexbinn
     &       .and.  zxb(jr) .lt. xwexbout ) then
          zwexb(jr) = wexbmax * 16.0
     &      * ( ( zxb(jr) - xwexbinn ) / ( xwexbout - xwexbinn ) )**2
     &      * ( ( xwexbout - zxb(jr) ) / ( xwexbout - xwexbinn ) )**2
        endif
      enddo
c
      npoints = medge
c
c
c..compute transport coefficients
c
      call mmm95( 
     &   zrminor,  zrmajor,  zelong
     & , zdense,   zdensh,   zdensz,   zdensfe
     & , zxzeff,   ztekev,   ztikev,   zq,       zbtor
     & , zavezimp, zmassimp, zmasshyd, zaimass,  zwexb
     & , zgrdne,   zgrdni,   zgrdnh,   zgrdnz,   zgrdte,  zgrdti, zgrdq
     & , zthiig,   zthdig,   ztheig,   zthzig
     & , zthirb,   zthdrb,   ztherb,   zthzrb
     & , zthikb,   zthdkb,   zthekb,   zthzkb
     & , zgamma,   zomega,   zdifthi,  zvelthi,  zflux
     & , matdim,   npoints,  nout,     lprint,   nerr
     & , lsuper,   lreset,   lswitch,  cswitch,  fig,     frb,   fkb)
c
c..long printouts
c
      write (nout,*)
      write (nout,*) '  OUTPUT FROM TEST_MMM95.TEX'
      write (nout,*)
c..print out namelist
c
      write (nout, 210) lswitch(1),  lswitch(2),  lswitch(3), lswitch(4) 
      write (nout, 220) lswitch(5),  lswitch(6),  lswitch(7), lswitch(8) 
      write (nout, 230) cswitch(1),  cswitch(2),  cswitch(3)
      write (nout, 240) cswitch(4),  cswitch(5),  cswitch(6)
      write (nout, 250) cswitch(7),  cswitch(8),  cswitch(9)
      write (nout, 260) cswitch(10), cswitch(11), cswitch(12)
      write (nout, 270) cswitch(13), cswitch(14), cswitch(15)
      write (nout, 280) cswitch(16), cswitch(17), cswitch(18)
      write (nout, 290) cswitch(19), cswitch(20), cswitch(21)
      write (nout, 300) cswitch(22), cswitch(23), cswitch(24)
      write (nout, 310) cswitch(25)

 210  format(//,2x, 'lswitch(1)=', I3,6x,'lswitch(2)=', I3,6x,
     &       'lswitch(3)=', I3,6x,'lswitch(4)=', I3)

 220  format(/,2x, 'lswitch(5)=', I3,6x,'lswitch(6)=', I3,6x,
     &       'lswitch(7)=', I3,6x,'lswitch(8)=', I3)

 230  format(//,2x, 'cswitch(1)=',1pe10.3,6x,'cswitch(2)=',1pe10.3,6x,
     &       'cswitch(3)=',1pe10.3)

 240  format(/,2x, 'cswitch(4)=', 1pe10.3,6x,'cswitch(5)=',1pe10.3,6x,
     &       'cswitch(6)=',1pe10.3)
 250  format(/,2x, 'cswitch(7)=', 1pe10.3,6x,'cswitch(8)=',1pe10.3,6x,
     &       'cswitch(9)=',1pe10.3)
 260  format(/,2x,'cswitch(10)=',1pe10.3,5x,'cswitch(11)=',1pe10.3,5x,
     &       'cswitch(12)=',1pe10.3)
 270  format(/,2x,'cswitch(13)=',1pe10.3,5x,'cswitch(14)=',1pe10.3,5x,
     &       'cswitch(15)=',1pe10.3)
 280  format(/,2x,'cswitch(16)=',1pe10.3,5x,'cswitch(17)=',1pe10.3,5x,
     &       'cswitch(18)=',1pe10.3)
 290  format(/,2x,'cswitch(19)=',1pe10.3,5x,'cswitch(20)=',1pe10.3,5x,
     &       'cswitch(21)=',1pe10.3)
 300  format(/,2x,'cswitch(22)=',1pe10.3,5x,'cswitch(23)=',1pe10.3,5x,
     &       'cswitch(24)=',1pe10.3)
 310  format(/,2x,'cswitch(25)=',1pe10.3,1pe10.3)
c
      write(nout, 320)  maxis,      medge,      rminor
      write(nout, 330)  rmajor,     elong,      denmin
      write(nout, 340)  temin,      denhaxis,   denhedge
      write(nout, 350)  denhexp,    amassh,     denzaxis
      write(nout, 360)  denzedge,   denzexp,    amassz
      write(nout, 370)  denfaxis,   denfedge,   denfexp
      write(nout, 380)  chrzaxis,   chrzedge,   chrzexp
      write(nout, 390)  chrfaxis,   chrfedge,   chrfexp
      write(nout, 400)  teaxis,     teedge,     teexp
      write(nout, 410)  tiaxis,     tiedge,     tiexp
      write(nout, 420)  qaxis,      qedge,      qexp
      write(nout, 430)  btor
      write(nout, 440)  fig(1),     fig(2),     fig(3),     fig(4)
      write(nout, 450)  frb(1),     frb(2),     frb(3),     frb(4)
      write(nout, 460)  fkb(1),     fkb(2),     fkb(3),     fkb(4)
      write(nout, 470)  lprint,     lsuper,     lreset
c
c
 320  format(/t2,'maxis=',i3,t27,'medge=',i3,t54,'rminor[m]=',1p1e9.2)
 330  format(/t2,'rmajor[m]=',f7.3,t27,'elong=',f7.3,t54,
     &       'denmin[m^-3]=',1p1e9.2)
 340  format(/t2,'temin[keV]=',1p1e9.2,t27,'denhaxis[m^-3]=',
     &        1p1e9.2,t54,'denhedge[m^-3]=',1p1e9.2)
 350  format(/t2,'denhexp=',1p1e9.2,t27,'amassh=',
     &        1p1e9.2,t54,'denzaxis[m^-3]=',1p1e9.2)
 360  format(/t2,'denzedge[m^-3]=',1p1e9.2,t27,'denzexp=',
     &        1p1e9.2,t54,'amassz=',1p1e9.2)
 370  format(/t2,'denfaxis[m^-3]=',1p1e9.2,t27,'denfedge[m^-3]=',
     &        1p1e9.2,t54,'denfexp=',1p1e9.2)
 380  format(/t2,'chrzaxis=',1p1e9.2,t27,'chrzedge=',
     &        1p1e9.2,t54,'chrzexp=',1p1e9.2)
 390  format(/t2,'chrfaxis=',1p1e9.2,t27,'chrfedge=',
     &        1p1e9.2,t54,'chrfexp=',1p1e9.2)
 400  format(/t2,'teaxis[keV]=',1p1e9.2,t27,'teedge[keV]=',
     &        1p1e9.2,t54,'teexp=',1p1e9.2)
 410  format(/t2,'tiaxis[keV]=',1p1e9.2,t27,'tiedge[keV]=',
     &        1p1e9.2,t54,'tiexp=',1p1e9.2)
 420  format(/t2,'qaxis=',1p1e9.2,t27,'qedge=',
     &        1p1e9.2,t54,'qexp=',1p1e9.2)
 430  format(/t2,'btor[T]=',1p1e9.2)
 440  format(/t2,'fig(1)=',0p1f5.3,t22,'fig(2)=',0p1f5.3,t42,
     &        'fig(3)=',0p1f5.3,t62,'fig(4)=',0p1f5.3)
 450  format(/t2,'frb(1)=',0p1f5.3,t22,'frb(2)=',0p1f5.3,t42,
     &        'frb(3)=',0p1f5.3,t62,'frb(4)=',0p1f5.3)
 460  format(/t2,'fkb(1)=',0p1f5.3,t22,'fkb(2)=',0p1f5.3,t42,
     &        'fkb(3)=',0p1f5.3,t62,'fkb(4)=',0p1f5.3)
 470  format(/t2,'lprint=',i3,t22,'lsuper=',i3,t42,'lreset=',i3)
c
	if (nerr .eq. -10) then 
	  write (nout,*) 'ABORT: Number of eqs. in Weiland model'
          write (nout,*) 'are incorrect. Only 2, 4 and 6-11 allowed.'
          write (nout,*) lswitch(1),' = lswitch(1) '
          stop
        elseif  (nerr .eq. -7) then
          write (nout,*) 'ABORT: Error on return from weiland14flux'
          stop
	elseif (nerr .eq. -6) then
          write(nout,*) 'ABORT:  ndim .gt. idim in sbrtn weiland14flux'
          stop
        endif
c
c..Print Header
c
        write (nout,*)
        write (nout,*)
     &   'Weiland-Nordman eigenvalue equations, sbrtn weiland14flux'
        write (nout,*)
c
          write (nout,*)
     &        '  Eigenvalues computed using ACM/TOMS routine 535'
c
        if ( lswitch(1) .eq. 11) then
          write (nout,*)
          write (nout,*)
     &     'Eleven eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z,'
          write (nout,*)
     &     '            F, Vp_H, Av, K, and Vp_Z'
        elseif ( lswitch(1) .eq. 10) then
          write (nout,*)
          write (nout,*)
     &     'Ten eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z,' 
          write (nout,*)
     &     '            F, Vp, Av, and K'
        elseif ( lswitch(1) .eq. 9) then
          write (nout,*)
          write (nout,*)
     &     'Nine eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z,'
          write (nout,*)
     &     '          F, Av, K, and Vp in strong ballooning limit'
        elseif ( lswitch(1) .eq. 8) then
          write (nout,*)
          write (nout,*)
     &     'Eight eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, and Vp'
        elseif ( lswitch(1) .eq. 7) then
          write (nout,*)
          write (nout,*)
     &     'Seven eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z, and F'
        elseif ( lswitch(1) .eq. 6 ) then
          write (nout,*)
          write (nout,*) 
     &     'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
        elseif ( lswitch(1) .eq. 5 ) then
          write (nout,*)
          write (nout,*)
     &     ' Five eqns for e phi/T_e, T_H, n_H, T_e, and Vp'
        elseif ( lswitch(1) .eq. 4 ) then
          write (nout,*)
          write (nout,*) ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
        elseif ( lswitch(1) .eq. 2 ) then
           write (nout,*)
           write (nout,*) ' Two eqns for e phi/T_e and T_H'
        endif
      if ( nerr .eq. 5 ) then
        write(nout,*) 
     &   'ifail .gt. 0 after call tomsqz in sbrtn weiland14flux'
      endif

      write (nout,110)
 110  format (/t7,'xb',t15,'rminor[m]',t29,'rmajor[m]',t45,'elong')
      do jr=1,medge
        write (nout,101) zxb(jr), zrminor(jr), zrmajor(jr)
     &    , zelong(jr)
      enddo
c
      write (nout,111)
 111  format (/t2,'rminor',t11,'zdense',t23,'zdensi'
     &  ,t35,'zdensh',t47,'zdensz',t59,'zdensf',t71,'zdensfe',
     &   /t4,'[m]',t11,'[m^-3]',t23,'[m^-3]',t35,'[m^-3]',
     &   t47,'[m^-3]',t59,'[m^-3]',t72,'[m^-3]')
      do jr=1,medge
        write (nout,102) zrminor(jr), zdense(jr)
     &    , zdensi(jr), zdensh(jr), zdensz(jr), zdensf(jr), zdensfe(jr)
      enddo
c
      write (nout,112)
 112  format (/t2,'rminor',t11,'ztekev',t23,'ztikev'
     &  ,t35,'zxzeff',t47,'zavezimp',t59,'zmassimp',t71,'zmasshyd',
     &   /t4,'[m]',t11,'[keV]',t23,'[keV]',t59,'[AMU]',t72,'[AMU]')
      do jr=1,medge
        write (nout,102) zrminor(jr)
     &    , ztekev(jr),   ztikev(jr),   zxzeff(jr)
     &    , zavezimp(jr), zmassimp(jr), zmasshyd(jr)
      enddo
c
      write (nout,114)
 114  format (/t2,'rminor',t11,'zq',t23,'zbtor'
     &  ,t35,'zwexb',
     &   /t4,'[m]',t11,' ',t23,'[tesla]',t35,'[sec^-1]')
      do jr=1,medge
        write (nout,102) zrminor(jr)
     &    , zq(jr),   zbtor(jr),   zwexb(jr)
      enddo
c
      write (nout,113)
 113  format (/t2,'rminor'
     &  ,t11,'zgrdne',t23,'zgrdni',t35,'zgrdnh'
     &  ,t47,'zgrdnz',t59,'zgrdte',t71,'zgrdti')
      do jr=1,medge
        write (nout,102) zrminor(jr)
     &    , zgrdne(jr), zgrdni(jr), zgrdnh(jr)
     &    , zgrdnz(jr), zgrdte(jr), zgrdti(jr)
      enddo
c
c
      write (nout,*)
      write (nout,*) ' Ion thermal diffusivities'
      write (nout,122)
 122  format (/t6,'rminor'
     &  ,t17,'chi_i_ig',t32,'chi_i_rb',t47,'chi_i_kb'
     &  ,t62,'chi_i_tot',/t7,'[m]',t17,'[m^2/s]',t32,'[m^2/s]',
     &   t47,'[m^2/s]',t62,'[m^2/s]')
 
      do jr=1,medge
        write (nout,103) zrminor(jr)
     &    , zthiig(jr), zthirb(jr), zthikb(jr)
     &    , zthiig(jr) + zthirb(jr) + zthikb(jr)
      enddo
c
      write (nout,*)
      write (nout,*) ' Electron thermal diffusivities'
      write (nout,124)
 124  format (/t6,'rminor'
     &  ,t17,'chi_e_ig',t32,'chi_e_rb',t47,'chi_e_kb'
     &  ,t62,'chi_e_tot',/t7,'[m]',t17,'[m^2/s]',t32,'[m^2/s]',
     &   t47,'[m^2/s]',t62,'[m^2/s]')

      do jr=1,medge
        write (nout,103) zrminor(jr)
     &    , ztheig(jr), ztherb(jr), zthekb(jr)
     &    , ztheig(jr) + ztherb(jr) + zthekb(jr)
      enddo
c
      write (nout,*)
      write (nout,*) ' Hydrogenic particle diffusivities'
      write (nout,126)
 126  format (/t6,'rminor'
     &  ,t18,'D_H_ig',t33,'D_H_rb',t48,'D_H_kb'
     &  ,t62,'D_H_tot',/t7,'[m]',t17,'[m^2/s]',t32,'[m^2/s]',
     &   t47,'[m^2/s]',t62,'[m^2/s]')
c
      do jr=1,medge
        write (nout,103) zrminor(jr)
     &    , zthdig(jr), zthdrb(jr), zthdkb(jr)
     &    , zthdig(jr) + zthdrb(jr) + zthdkb(jr)
      enddo
c
      write (nout,*)
      write (nout,*) ' Impurity particle diffusivities'
      write (nout,128)
 128  format (/t6,'rminor'
     &  ,t18,'D_Z_ig',t33,'D_Z_rb',t48,'D_Z_kb'
     &  ,t62,'D_Z_tot',/t7,'[m]',t17,'[m^2/s]',t32,'[m^2/s]',
     &   t47,'[m^2/s]',t62,'[m^2/s]')
      do jr=1,medge
        write (nout,103)zrminor(jr)
     &    , zthzig(jr), zthzrb(jr), zthzkb(jr)
     &    , zthzig(jr) + zthzrb(jr) + zthzkb(jr)
      enddo
c
      write (nout,*)
      write (nout,*) ' Diffusion matrix:'
      write (nout,141)
 141  format (/t2,'rminor'
     &  ,t10,'vel(1,jr)',t21,'dif(1,1,jr)',t33,'dif(1,2,jr)'
     &  ,t45,'dif(1,3,jr)',t57,'dif(1,4,jr)',t69,'dif(1,5,jr)',
     &   /t4,'[m]',t12,'[m/s]',t23,'[m^2/s]',t35,'[m^2/s]',
     &   t47,'[m^2/s]',t59,'[m^2/s]',t71,'[m^2/s]')

      do jr=1,medge
        write (nout,102) zrminor(jr)
     &    , zvelthi(1,jr),   zdifthi(1,1,jr), zdifthi(1,2,jr)
     &    , zdifthi(1,3,jr), zdifthi(1,4,jr), zdifthi(1,5,jr)
      enddo
c
      write (nout,*)
      write (nout,142)
 142  format (/t2,'rminor'
     &  ,t10,'vel(2,jr)',t21,'dif(2,1,jr)',t33,'dif(2,2,jr)'
     &  ,t45,'dif(2,3,jr)',t57,'dif(2,4,jr)',t69,'dif(2,5,jr)',
     &   /t4,'[m]',t12,'[m/s]',t23,'[m^2/s]',t35,'[m^2/s]',
     &   t47,'[m^2/s]',t59,'[m^2/s]',t71,'[m^2/s]')

      do jr=1,medge
        write (nout,102) zrminor(jr)
     &    , zvelthi(2,jr),   zdifthi(2,1,jr), zdifthi(2,2,jr)
     &    , zdifthi(2,3,jr), zdifthi(2,4,jr), zdifthi(2,5,jr)
      enddo
c
      write (nout,*)
      write (nout,143)
 143  format (/t2,'rminor'
     &  ,t10,'vel(3,jr)',t21,'dif(3,1,jr)',t33,'dif(3,2,jr)'
     &  ,t45,'dif(3,3,jr)',t57,'dif(3,4,jr)',t69,'dif(3,5,jr)',
     &   /t4,'[m]',t12,'[m/s]',t23,'[m^2/s]',t35,'[m^2/s]',
     &   t47,'[m^2/s]',t59,'[m^2/s]',t71,'[m^2/s]')

      do jr=1,medge
        write (nout,102) zrminor(jr)
     &    , zvelthi(3,jr),   zdifthi(3,1,jr), zdifthi(3,2,jr)
     &    , zdifthi(3,3,jr), zdifthi(3,4,jr), zdifthi(3,5,jr)
      enddo
c
      write (nout,*)
      write (nout,144)
 144  format (/t2,'rminor'
     &  ,t10,'vel(4,jr)',t21,'dif(4,1,jr)',t33,'dif(4,2,jr)'
     &  ,t45,'dif(4,3,jr)',t57,'dif(4,4,jr)',t69,'dif(4,5,jr)',
     &   /t4,'[m]',t12,'[m/s]',t23,'[m^2/s]',t35,'[m^2/s]',
     &   t47,'[m^2/s]',t59,'[m^2/s]',t71,'[m^2/s]')

      do jr=1,medge
        write (nout,102) zrminor(jr)
     &    , zvelthi(4,jr),   zdifthi(4,1,jr), zdifthi(4,2,jr)
     &    , zdifthi(4,3,jr), zdifthi(4,4,jr), zdifthi(4,5,jr)
      enddo
c
      write (nout,*)
      write (nout,145)
 145  format (/t2,'rminor'
     &  ,t10,'vel(5,jr)',t21,'dif(5,1,jr)',t33,'dif(5,2,jr)'
     &  ,t45,'dif(5,3,jr)',t57,'dif(5,4,jr)',t69,'dif(5,5,jr)',
     &   /t4,'[m]',t12,'[m/s]',t23,'[m^2/s]',t35,'[m^2/s]',
     &   t47,'[m^2/s]',t59,'[m^2/s]',t71,'[m^2/s]')

      do jr=1,medge
        write (nout,102) zrminor(jr)
     &    , zvelthi(5,jr),   zdifthi(5,1,jr), zdifthi(5,2,jr)
     &    , zdifthi(5,3,jr), zdifthi(5,4,jr), zdifthi(5,5,jr)
      enddo
c
      write (nout,*)
      write (nout,*) ' Total thermal transport:'
      write (nout,130)
 130  format (/t6,'rminor',t18,'chi_i',t33,'vel_i',t48,'chi_e'
     &    ,t62,'vel_e',/t7,'[m]',t17,'[m^2/s]',t33,'[m/s]',
     &   t47,'[m^2/s]',t62,'[m/s]')

      do jr=1,medge
        write (nout,103) zrminor(jr)
     &    , zdifthi(1,1,jr), zvelthi(1,jr)
     &    , zdifthi(3,3,jr), zvelthi(3,jr)
      enddo
c
      write (nout,*)
      write (nout,*) ' Total particle transport:'
      write (nout,131)
 131  format (/t6,'rminor',t18,'D_H',t33,'vel_H',t48,'D_Z',t62,'v_Z'
     & ,/t7,'[m]',t17,'[m^2/s]',t33,'[m/s]',t47,'[m^2/s]',t62,'[m/s]')
      do jr=1,medge
        write (nout,103) zrminor(jr)
     &    , zdifthi(2,2,jr), zvelthi(2,jr)
     &    , zdifthi(4,4,jr), zvelthi(4,jr)
      enddo
c
c
      write (nout,*)
      write (nout,*) ' Total vFlux:'
      write (nout,132)
 132  format (/t6,'rminor',t18,'vFlux_i',t33,'vFlux_H',t48,
     &      'vFlux_e',t62,'vFlux_Z',/t7,'[m]',t19,'[m/s]',t34,'[m/s]',
     &   t49,'[m/s]',t63,'[m/s]')
      do jr=1,medge
        write (nout,103) zrminor(jr)
     &    , zflux(1,jr), zflux(2,jr),  zflux(3,jr),  zflux(4,jr)
      enddo
c
c
c
 101  format (0p2f10.4,1p6e17.4)
 102  format (0p1f6.3,1p6e12.3)
 103  format (0p1f10.4,1p6e15.4)
c
      stop
      end
!| 
!| \end{document}
