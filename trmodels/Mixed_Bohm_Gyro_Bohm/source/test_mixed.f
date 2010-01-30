!| %
!| %  This is a LaTeX ASCII file.  To typeset document type:
!| % latex test_mixed.tex
!| %  To extract the fortran source code, type:
!| % python ./s2tex.py test_mixed.f
!| %
!| % test_mixed.f  Glenn Bateman  Lehigh University
!| %---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!| \documentstyle{article}
!| \headheight 0pt \headsep 0pt  \topmargin 0pt  \textheight 9.0in
!| \oddsidemargin 0pt \textwidth 6.5in
!| \begin{document}
!| 
!| \begin{center}
!| {\bf {\tt test-mixed.tex} \\
!| Stand-alone Code to Test the Mixed Bohm/gyroBohm model } \\
!| \vspace{1pc}  Glenn Bateman and Arnold Kritz \\
!| Lehigh University \\
!| Version 1.0$\quad$ 22 August 2003
!| \end{center}
!| 
!| This is a stand-alone program to test Mixed Bohm/gyroBohm model.
!| 
!| For the purposes of this test program, only two ion species are
!| considered: a thermal hydrogen species and thermal impurity species
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
c  Call tree: TEST_MIXED calls the following routines
c
c  mixed_model       - Computes transport coefficients
c    mixed           - Computes diffusivity from the Mixed Bohm/gyro-Bohm Model
c
c-----------------------------------------------------------------------
c
      use mixed_Bohm_gyro_Bohm
c
      implicit none
c
      integer      kr,       km
      parameter ( kr = 55, km = 12)
c
c  kr = max number of radial points
c  km = first dimension of diffusivity matrix
c
      integer
     &   nin,      nout,     ntemp,    maxis,     medge
     & , npoints,  lflowshear, nerr
c
      real
     &   rminor,   rmajor,   denmin,   temin,  amassh, charge
     & , deneaxis, deneedge, deneexp
     & , teaxis,   teedge,   teexp,    tiaxis,   tiedge,   tiexp
     & , wexbmax,  xwexbinn, xwexbout
     & , qaxis,    qedge,    qexp,     btor
c
c..Namelist input:
c
c  maxis   = zone boundary at magnetic axis
c  medge   = zone boundary at edge of plasma
c
c  rminor = minor radius (half-width) at edge of plasma
c  rmajor = major radius to geometric center at edge of plasma
c
c  denmin = minimum value allowed for n_e
c  temin  = minimum value allowed for T_e, T_i, and T_f
c
c    In this test program, the profiles are all taken to be
c  parabolas raised to powers:   For example
c  T_i = (tiaxis-tiedge)*(1.-x^2)^tiexp + tiedge
c  where x = sqrt normalized toroidal flux at zone boundaries
c
c deneaxis   electron density at magnetic axis [m^-3]
c deneedge   electron density at plasma edge [m^-3]
c deneexp    exponential of parabola for profile
c
c amassh     isotopic mass of hydrogenic species
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
c lflowshear = 0 for no flow shear stabilization
c            = 1 for the flow shear stabilizatioin model of
c                 T.J. Tala et al Plasma Phys. Controlled Fusion 44
c                 (2002) A495]
c
c nerr       error if nerr returns not equal to 0
c
c    The magnetic q profile is given by
c  q(x) = qaxis + ( qedge - qaxis ) * x**(qexp)
c
c
      integer        jr
c
      real           zdx,          zxfact,       zepslon
     & , zt_e_kev_edge  ,          zrminor_edge
c
      real
     &   zxb(kr),     zrminor(kr),   zrmajor(kr),  zdense(kr)
     & , ztekev(kr),  ztikev(kr),    zq(kr),       zbtor(kr)
     & , zwexb(kr),   zmass_hyd(kr), zcharge_hyd(kr)
c
c
c    All the following 1-D arrays are assumed to be on zone boundaries
c    including the densities (in m^-3) and temperatures (in keV).
c
c  zxb(jr)     = normalized minor radius ( r / a ) at zone boundaries
c
c  zrminor(jz) = minor radius (half-width) of zone boundary [m]
c  zrmajor(jz) = major radius to geometric center of zone bndry [m]
c
c  zdense(jz)  = electron density [m^-3]
c  ztekev(jz)  = T_e [keV]  (electron temperature)
c  ztikev(jz)  = T_i [keV]  (temperature of thermal ions)
c  zq(jz)      = magnetic q-value
c  zbtor(jz)   = ( R B_tor ) / rmajor(jz)  [tesla]
c
c  zwexb(jz)   = flow shear rate [radians/sec]
c
c  zmass_hyd(jz)   = mean atomic mass of main thermal ions [AMU]
c  zcharge_hyd(jz) = ionic charge of main ions ( = 1.0 for hydrogenic ions)
c
c..normalized gradients:
c
      real   zgrdne(kr), zgrdte(kr),   zshear(kr)
c
c    All of the following normalized gradients are at zone boundaries.
c    r = half-width, R = major radius to center of flux surface
c    carry out any smoothing before calling sbrtn theory
c
c  zgrdne(jz) = - R ( d n_e / d r ) / n_e
c  zgrdte(jz) = - R ( d T_e / d r ) / T_e
c  zshear(jz) =   r ( d q   / d r ) / q     magnetic shear
c
c..transport output:
c
      real
     &   chi_i_gyro_bohm(kr),   chi_e_gyro_bohm(kr)
     & , chi_i_bohm(kr), chi_e_bohm(kr)
     & , chi_i_mix(kr),  chi_e_mix(kr),   d_hyd_mix(kr)
c
c  effective diffusivites:
c
c  chi_i_mix = ion thermal effective diffusivity [m^2/sec]
c  chi_e_mix = electron thermal effective diffusivity [m^2/sec]
c  d_hyd_mix = ion particle effective diffusivity [m^2/sec]
c
c  values of gyro-Bohm and Bohm terms for diagnostic purposes only:
c
c  chi_i_gyro_bohm   = gyro-Bohm contribution to chi_i_mix [m^2/sec]
c  chi_e_gyro_bohm   = gyro-Bohm contribution to chi_e_mix [m^2/sec]
c  chi_i_bohm = Bohm contribution to chi_i_mix [m^2/sec]
c  chi_e_bohm = Bohm contribution to chi_e_mix [m^2/sec]
c
c..namelist input:
c
      namelist /nt/
     &   maxis,    medge,    rminor,   rmajor
     & , denmin,   temin,    amassh,   charge
     & , deneaxis, deneedge, deneexp
     & , teaxis,   teedge,   teexp
     & , tiaxis,   tiedge,   tiexp
     & , qaxis,    qedge,    qexp
     & , wexbmax,  xwexbinn, xwexbout
     & , btor,     lflowshear
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
      nerr    =  0
c
c  for subroutine mixed_model
c..default input values
c
      maxis   = 0         ! index of magnetic axis
      medge   = 20        ! index of edge of plasma

      rminor  = 0.9       ! minor radius [m]
      rmajor  = 2.5       ! major radius [m]

      denmin  = 0.0       ! minimum density = 1.e7 [keV]
      temin   = 0.0       ! minimum temperature = 1.e-6 [keV]

      deneaxis = 0.3e19   ! hydrogen density at magnetic axis [m^-3]
      deneedge = 1.0e18   ! hydrogen density at plasma edge [m^-3]
      deneexp  = 1.0      ! exponential of parabola for profile [m^-3]

      amassh   = 2.0141   ! Mass of Deuterium in amu units
                          ! isotopic mass of hydrogenic species
      charge   = 1.0      ! hydrogenic charge

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
c
      lflowshear = 0      ! no flow shear stabilization
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
      enddo
c
c..construct electron density and temperature profiles on zone boundaries
c
      do jr=1,medge-1
        zdense(jr) = deneedge
     &    + (deneaxis-deneedge) * (1.0 - zxb(jr)**2)**deneexp
        ztekev(jr) = max ( temin, teedge
     &    + (teaxis-teedge) * (1.0 - zxb(jr)**2)**teexp )
        ztikev(jr) = max ( temin, tiedge
     &    + (tiaxis-tiedge) * (1.0 - zxb(jr)**2)**tiexp )
      enddo
c
        zdense(medge) = deneedge
        ztekev(medge) = teedge
        ztikev(medge) = tiedge
c
c..normalized gradients
c
      zepslon = 1.e-9
c
      do jr=1,medge
c
        zxfact = max ( zepslon, 1.0 - zxb(jr)**2 )
c
        zgrdne(jr) = zrmajor(jr) * 2.0 * zxb(jr) * deneexp *
     &    (deneaxis-deneedge) * zxfact**(deneexp-1.0)
     &      / ( max( denmin, zdense(jr) ) * zrminor(medge) )
c
        zgrdte(jr) = zrmajor(jr) * 2.0 * zxb(jr) * teexp *
     &    (teaxis-teedge) * zxfact**(teexp-1.0)
     &      / ( ztekev(jr) * zrminor(medge) )
c
      enddo
c
c..average mass and charge
c
      do jr=1,medge
        zcharge_hyd(jr)  = charge
        zmass_hyd(jr)    = 2.0141 ! Mass of Deuterium in amu units
      enddo
c
c
c..simple magnetics profiles
c
      do jr=maxis+1,medge
        zq(jr) = qaxis + (qedge-qaxis) * abs(zxb(jr))**qexp
        zshear(jr) = zrminor(jr) *
     &    (qedge-qaxis) * (abs(zxb(jr)))**(qexp-1.0)
     &      / ( zq(jr) * zrminor(medge) )
      enddo
c
c
      if ( maxis .gt. 0 ) then
        zq(maxis) = qaxis
        zshear(maxis) = 0.0
      endif
      if ( maxis .gt. 1 ) then
        zq(maxis-1) = zq(maxis+1)
        zshear(maxis-1) = - zshear(maxis+1)
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
c..edge values
c
      zt_e_kev_edge = ztekev(medge)
      zrminor_edge  = zrminor(medge)
c
      npoints = medge
c
c
c..compute transport coefficients
c
      call mixed_model (
     &   zrminor,  zrmajor
     & , ztekev,   ztikev,    zq
     & , zbtor,    zmass_hyd, zcharge_hyd,  zwexb
     & , zgrdte,   zgrdne,    zshear
     & , zt_e_kev_edge, zrminor_edge
     & , npoints
     & , chi_i_mix,  chi_e_mix, d_hyd_mix
     & , chi_i_gyro_bohm,   chi_e_gyro_bohm
     & , chi_i_bohm, chi_e_bohm
     & , nerr
     & , lflowshear=lflowshear)
c
c..long printouts
c
      write (nout,*)
      write (nout,*) '  OUTPUT FROM TEST_MIXED.F'
      write (nout,*)
c..print out namelist
c
c
      write(nout, 320)  maxis,      medge,      rminor
      write(nout, 330)  rmajor,     denmin,     temin
      write(nout, 340)  deneaxis,   deneedge,   deneexp
      write(nout, 350)  amassh,     charge
      write(nout, 400)  teaxis,     teedge,     teexp
      write(nout, 410)  tiaxis,     tiedge,     tiexp
      write(nout, 420)  qaxis,      qedge,      qexp
      write(nout, 430)  btor
c
c
 320  format(/t2,'maxis=',i3,t27,'medge=',i3,t54,'rminor[m]=',1p1e9.2)
 330  format(/t2,'rmajor[m]=',f7.3,t27,'denmin[m^-3]=',1p1e9.2
     &  ,t54,'temin[keV]=',1p1e9.2)
 340  format(/t2,'deneaxis[m^-3]=',1p1e9.2
     &        ,t27,'deneedge[m^-3]=',1p1e9.2
     &        ,t54,'deneexp=',1p1e9.2)
 350  format(/t2,'amassh=',1p1e9.2,t27,'charge=',1p1e9.2)
 400  format(/t2,'teaxis[keV]=',1p1e9.2,t27,'teedge[keV]=',
     &        1p1e9.2,t54,'teexp=',1p1e9.2)
 410  format(/t2,'tiaxis[keV]=',1p1e9.2,t27,'tiedge[keV]=',
     &        1p1e9.2,t54,'tiexp=',1p1e9.2)
 420  format(/t2,'qaxis=',1p1e9.2,t27,'qedge=',
     &        1p1e9.2,t54,'qexp=',1p1e9.2)
 430  format(/t2,'btor[T]=',1p1e9.2)
c
      select case (nerr)
        case (1)
         write(nout,*) 'test_mixed input error> Te < 0'
         stop 'test_mixed input error> Te < 0'
        case (2)
         write(nout,*) 'test_mixed input error> Ti < 0'
         stop 'test_mixed input error> Ti < 0'
        case (3)
         write(nout,*) 'test_mixed input error> Te_p8 < 0'
         stop 'test_mixed input error> Te_p8 < 0'
        case (4)
         write(nout,*) 'test_mixed input error> Te_edge < 0'
         stop
        case (5)
         write(nout,*) 'test_mixed input error> rmaj < 0'
         stop 'test_mixed input error> rmaj < 0'
      end select
c
c..Print Header
c
      write (nout,110)
 110  format (/t7,'xb',t15,'rminor[m]',t29,'rmajor[m]')
      do jr=1,medge
        write (nout,101) zxb(jr), zrminor(jr), zrmajor(jr)
      enddo
c
      write (nout,111)
 111  format (/t2,'rminor',t11,'zdense',t23,'ztekev',t35,'ztikev'
     &  ,t47,'zmass_hyd'
     &   /t4,'[m]',t11,'[m^-3]',t23,'[keV]',t35,'[keV]',
     &   t47,'AMU')
      do jr=1,medge
        write (nout,102) zrminor(jr), zdense(jr)
     &   , ztekev(jr),   ztikev(jr) , zmass_hyd(jr)
      enddo
c
      write (nout,114)
 114  format (/t2,'rminor',t11,'zbtor',t23,'zq',t35,'shear'
     &  ,t47,'zwexb',
     &   /t4,'[m]',t11,'[tesla]',t23,' ',t35,' ',t47,'[sec^-1]')
      do jr=1,medge
        write (nout,102) zrminor(jr)
     &    ,   zbtor(jr), zq(jr), zshear(jr),   zwexb(jr)
      enddo
c
      write (nout,113)
 113  format (/t2,'rminor'
     &  ,t11,'zgrdte',t23,'zgrdne')
      do jr=1,medge
        write (nout,102) zrminor(jr)
     &       , zgrdte(jr), zgrdne(jr)
      enddo
c
c
      write (nout,*)
      write (nout,*)
     &   ' Ion thermal diffusivities: Bohm and gyro-Bohm terms'
      write (nout,122)
 122  format (/t6,'rminor'
     &  ,t17,'chi_i_gB',t32,'chi_e_gB',t47,'chi_i_B', t60,'chi_e_B'
     &  ,/t7,'[m]',t17,'[m^2/s]',t32,'[m^2/s]',
     &   t47,'[m^2/s]',t60,'[m^2/s]')

      do jr=1,medge
        write (nout,103) zrminor(jr)
     & , chi_i_gyro_bohm(jr),   chi_e_gyro_bohm(jr)
     & , chi_i_bohm(jr), chi_e_bohm(jr)
      enddo
c
      write (nout,*)
      write (nout,*)
     & ' Total thermal and hydrogenic particle diffusivities'
      write (nout,124)
 124  format (/t6,'rminor'
     &  ,t17,'chi_i',t32,'chi_e',t47,'D_H',
     &  /t7,'[m]',t17,'[m^2/s]',t32,'[m^2/s]',t47,'[m^2/s]')

      do jr=1,medge
        write (nout,103) zrminor(jr)
     &    , chi_i_mix(jr),  chi_e_mix(jr),   d_hyd_mix(jr)
      enddo
c
 101  format (0p2f10.4,1p6e17.4)
 102  format (0p1f6.3,1p6e12.3)
 103  format (0p1f10.4,1p6e15.4)
c
      stop
      end
!| %
!| % \newpage
!| %
!| \begin{center}
!| \Large {\bf Subroutine STRIPX} \\
!| \vspace{1pc} \normalsize
!| Glenn Bateman \\
!|  Lehigh University, Physics Department \\
!| 16 Memorial Drive East, Bethlehem PA 18015 \\
!| bateman@plasma.physics.lehigh.edu \\
!| \end{center}
!| 
!| This file documents a subroutine sbrtn that reads input data file
!| logical unit number nin one line at a time up to 80 characters long,
!| prints out the input verbatum to output file logical unit number nout
!| then searches for the first appearance of the character ! on each line
!| and prints out everything to the left of ! to output unit number ntrans
!| In this way, everything to the right of ! is treated as a commentwhich computes
!| plasma transport coefficients using the Mixed transport model.
!| 
      subroutine stripx (nin,ntrans,nout)
c
      parameter (kc = 132)
c
      character line * 132
c
  10  read (nin,100,end=20) line
 100  format (a132)
c
c..find number of characters before spaces
c
      ilength = 0
      do j=1,kc
        if ( line(j:j) .ne. ' ' ) ilength = j
      enddo
c
c..echo input data on output file
c
      if ( ilength .gt. 0 ) then
        write (nout,110) line(1:ilength)
      else
        write (nout,*)
      endif
c
c..ixlen = number of nontrivial characters before "!"
c
      ixlen = 0
      do j=1,kc
        if ( line(j:j) .eq. '!' ) go to 14
        if ( line(j:j) .ne. ' ' ) ixlen = j
      enddo
  14  continue
c
c..echo input data up to "!" on temporary file
c
      if ( ixlen .gt. 0 ) write (ntrans,110) line(1:ixlen)
c
 110  format (a)
c
c
      go to 10
c
  20  continue
      rewind ntrans
c
      return
      end
!| 
!| \end{document}
