       TRANSPORT MODEL (Weiland- Nordman model for reactive drift waves)
C      modified version October  1997 
C -----------------------------------------------------------------
      This is a transport model for the good confinement region
      describing transport of ion temperature, electron temperature,
      electron density, auxiliary ion temperature and auxiliary ion
      density. A 5x5 transport matrix is calculated with the dynamic 
      quantities in this order. 
      The electrons are free (Boltzmann) and trapped with the fraction FT.
      Both electron species are assumed to have the same temperature
      and Xe is the electron thermal conductivity seen as corresponding
      to all electrons. The auxiliary ions are treated with the same
      accuracy as the main ions. The transport is due to toroidal eta-i modes
      for main and auxiliary ions and trapped electron modes using the same
      two-pole reactive fluid model for all species. This leads do a six
      degree dispersion relation for the collisionless case (NDISP=6) and
      to a seventh degree dispersion relation with collisions (NDISP=7).
      The code can also be run with only one ion species (NDISP=4) or with
      no electron trapping and one ion species (NDISP=2). The dispersion
      relation is solved in DISP7P(NDISP,ZZ) where ZZ is a complex vector
      containing the roots. This programme was in the collisionless version
      without parallel ion motion written by Glenn Bateman. DISP7P uses a NAG 
      routine for solving systems of coupled equations that are linear in
      omega. Parallel ion motion is included by solving the eigenvalue problem 
      analytically in the strong ballooning limit. This gives a dependencs on q
      and S but does not give the improvement for negative shear. Parallel ion
      motion is only included in the seven equation set. The transport involves
      Ti,Te,Ne,Tq and Nq, using a 5 x 5 transport matrix represented in
      vector form CHI(I), CHE(I), D(I), CHQ(I) and DQ(I), I = 1 to 5 .
      The ion temperature flux, however, also requires a nondiffusive
      part as explained in the beginning of the subroutine DIFF.
      The main ion density Ni is calculated from quasineutrality.
      The diffusion matrix is calculated in the subroutine DIFF.
      Also the effective diffusivities are calculated. These are
      denoted SCHI, SCHE, SD, SCHQ and SDQ. The routine DIFF sums
      all transport coefficients over the unstable modes.
C
      Because of numerical difficulties with impurities, the impurity
      density is usually kept fixed or at a given fraction of the electron
      density. As it turns out the latter choice may give considerably
      better agreement for the ion temperature. The impurity temperature
      is integrated self consistently. For numerical stability reasons
      it has, however, been necessary to divide the heating power between
      main ions and impurity ions so that they remain approximately the same.
C
     In the most recent version of the code, the flux surface geometry has been      included. Here also the dependence of magnetic field and curvature
     radius on the small radius have been included.
C
     In the present version the convective fluxes have been extracted and
     put in the commonlist COMMON/TP/. They should be added according to the
     transport equations on page 4 in the document "Standard list of variables 
     and file format". This holds exactly for the density while the energy
     equation here is replaced by the temperature equations. This means that a
     factor 2/3 is absorbed into the thermal conductivities and the factor 3/2 
     in front of the convective heatflux is absent. It also means that the
     convective parts, corresponding to Gamma, do not contain the density, ie
     they are pure temperature fluxes. In COMMON/TP/, EA is the inverse
     aspect ratio a/R, HPS is the convective ion temperature flux, PES
     is the convective electron temperature flux, PNS is the convective 
     electron density flux, PQS is the convective impurity temperature flux and
     PNQS is the convective impurity density flux.
C
C  ******************************************************
C As an example we here give the input data used in the TOR code 
C at Chalmers University of Technology
C
      The input data in TRANSTOR.START are
C
      E                                  a/R
      N0         initial peak electron density
      T0I        edge Ti (fixed)
      T0E        edge Te (fixed)
      TCI        initial peak Ti
      TCE        initial peak Te
      TM         maximum time in seconds
      PHI        aux heating power to ions (MW)
      PHE        aux heating power to electrons (MW)
      SAN        strength of interior electron source
      W          exponent for initial elect. density prof.
      RIW        code for printout (1 gives full printout)
      TAU        redundant
      FL         flr parameter (k*ro)**2
      XA         location of outer boundary (X = 1 - r/a)
      XB         location of inner boundary
      XSI        location of inner ion heat source
      XSE        location of inner electron heat source.
      XSN        location of inner particle source (electr-ion)
      RIT        number of space points
      CX         amplitude of external X from XL to the axis  (sawteeth)
      TWR        time after which RIW=1 (at stationary state)
      SEI        relative strength of edge ion heat source
      SEE        relative strength of edge electron heat source
      SEN        edge particle source strength
      KT         initial time step in seconds
      DX         normalisation factor for transport accord. to mashin par.
      PR         horisontal plasma radius ( a ) in meters
      ZE         Z value for auxiliary ions
      KNC        factor multiplying neoclassical diffusivity
      COL        factor multiplying collisionality 
      RIB        code for internal boundary cond.
      KP         elongation (b/a)
      TEQ        factor multiplying the resistive temp. equilibration term
      TPL        time in seconds between global printouts and plots
      BF         initial fraction of auxiliary ions
      EV         loop voltage
      T0Q        edge  Tq (fixed)
      TCQ        initial peak temp of aux ions
      RNEQ       number of equations (NDISP)
      SEQ        strength of auxiliary ion edge source
      XL         x-point where the artificial diffusion starts (x=1. -r/a)
      RAD        factor that multiplies the radiation sink.
      BTOR       toroidal vacuum magnetic field in tesla
      ROT        factor that multiplies the EXB shearing rate.
C  --------------------------------------------------------------
C     These input data are particular to the TOR code but may give
C     some hints to how the model has been implemented

      When input data are taken from the ITER database the equilibrium
      profiles are usually taken as initial profiles. Then the parameters
      related to initial profiles and heat or particle sources are not
      used. An exception is the edge electron source (SEN) which is kept
      as a free parameter. When the impurity diffusion is included also
      the sources related to impurities are used
C  *********************************************************************
C
     Most notations are standard. Note, however, the definition
     EN = omegad/omega* = 2*Ln/Lb (epsilon-n). The index q (Q) is used
     for auxiliary ions.

C    --------------------------------------------------------
C   The subroutine disp7p(neq,ZZ) needs the following input.
C
C    Argument list:  neq = number of equations linear in omega
C    which define the linear dispersion relation, neq< or = 7. 
C
C    ZZ(7) (complex) Contains the solution of the dispersion
C    relation. Frequencies are normalized by the electron 
C    magnetic drift frequency omega_de.
C
C    COMMON/GRAD/ EN,ENI,EI,EE          
C    EN = 2*Ln/Lb , ENH = 2*Lni/LB, EI = Lni/Lti, EE = Ln/Lte
C
C    Here Ln,Lt,Lb etc denote scale lengths of density, temperature,
C    magnetic field etc. and i stands for main ions, e for electrons
C    and q for auxiliary ions. Index n stands for electron density NE.
C    Snce ENI is used for the inverse of EN, ENH is used for the EN of main 
C    ions. In the code capitals are used also for indices.
C
C    COMMON/IMP/ BQ,G,SI,Z,TAUZ,ENQ,EQ,KIQ,KXQ
C    BQ = NQ/NE, G = 1 - Z*BQ, SI=0 , TAUZ = TE/TQ, ENQ =2*LNQ/LB
C    EQ = LNQ/LTQ, KIQ = LNQ/LNI, KXQ = LN/LNQ, Z is the charge
C    number for the auxiliary ions.
C
     COMMON/ETAWN6/ CETAIN,ETE,ETI,ETQ,TAUH,AZ,FT,RFLT,LPRINTIN,NDIM
C    CETAIN(32) control vector. Present setting: first 31 elements zero
C    CETAIN(32) = 0.001 which controlles the accuracy in the NAG routine
C    f02bje. ETE = 2*LTE/LB, ETI = 2*LTI/LB, ETQ = 2*LTQ/LB, TAUH = TE/TI
C    AZ = mass number of aux. ions, FT = fraction of trapped electrons,
C    FFLT = (Kperp*rhos)**2 = 0.1 (FLR par.), LPRINTIN = control variable
C    for printout. Larger values generate more output, NDIM = 5 (dimension
C    of transport matrix).
C
C    COMMON/COLL/ BTA,VEF
C    BTA = 1.5 , VEF = effevtive electron ion collision frequency for
C    trapped electrons, normalized by omega_de.
C    BTA is a parameter used in the collision model
C
C    COMMON/BETAE/ betae, q, S, Cs, em
C    This common is used also for an electromagnetic version of the code
C    which is now beeing tested so betae and em are not used here.
C    q is the safety factor, S is the magnetic shear parameter and Cs is
C    the ion acoustic velocity for hydrogen in meters/s.
C
C    OUTPUT VARIABLE 
C    COMMON/IRET/ IRET
C    IRET nonzero indicates that the routine f02bje did not find the roots.
C
C  -----------------------------------------------------
C    INPUT NEEDED FOR
C    SUBROUTINE DIFF(RP,IR,TAUI,FT,U,CHI,CHE,D,CHQ,DHQ)
C
C    ARGUMENT LIST:
C    RP(7) contains the unstable roots found by disp7.
C    IR is the number of unstable roots. TAUI=TI/TE, FT has been defined above.
C    U(J,IX) contains the radial profiles of the dynamic variables.
C    Here             J = 1  corresponds to main ion temperature
C                     J = 2      "           electron temperature
C                     J = 3      "           electron density
C                     J = 4      "           aux. ion temperature
C                     J = 5      "           aux. ion density
C
      IX is the space point index, starting from the edge.
C     The remaining variables in the list are the output vectors
C     that define the transport matrix as discussed above and in the 
C     beginning of DIFF.
C
C     COMMON BLOCKS
C     COMMON/GRAD/  same as above
C     COMMON/COLL/ same as above
C     COMMON/IMP/ same as above
C     COMMON/PAR/ THRD,FTR,STR,RFL,D1,XIH,IX,IW
C     THRD = 1/3, FTR = 5/3, STR = 7/3, RFL = 0.1 (FLR par), 
C     D1 = (2.*rhos**2*Cs/LB)/TE**1.5, where rhos is the ion
C     gyroradius at the electron temp. and Cs is the ion acoustic velocity,
C     D1 is a machin dependent parameter in SI units. As an example, for JET
C     with B = 3.5 Tesla and LB = R = 3m we obtain D1 = 0.175.
C     XIH is the kernel of the transport matrix components sometimes
C     needed in discretizing the transport equations (output variable)
C     IX = current space point number, IW = 1 gives printout
C
C     COMMON/ETAWN6/  same as above. Note that the trapped fraction
C     was included before in the argument list. Here we use another
C     variable name so that we can keep the same common block as in other 
C     subroutines.
C
      COMMON/WROT/ WEXB,ROT
C     WEXB is the EXB shearing rate according to Hahm and Burrell. It is
      subtracted from the growthrate before the transport is calculated.
      ROT is a factor multiplying WEXB for turning it off.
C
      COMMON/SCALE/ NLX
      NLX is a denominator in a scale-length. It has been given a small
      value with the correct sign if its magnitude is too small.
C
      COMMON/GRKVOT/ GK
      GK is the ratio GRHO1/(a*GRHO2)
      where a is the small radius. GK=1 for cylindrical flux surfaces.
      It is an input parameter needed to calculate the effective transport
      coefficients in COMMON/EFF DIFF/.
C
C     COMMON/EFF DIFF/ SCHI,SCHE,SD,SCHQ,SDQ
C     OUTPUT VARIABLES DEFINED ABOVE
C
      COMMON/TP/ EA,HPS,PES,PNS,PQS,PNQS

      EA is the inverse aspect ratio a/R. The remaining variables are the
      convective fluxes for Ti,Te,Ne,Tq and Nq. It would be preferrable to put
      these in a vector. The interpretation is given above.

C  ------------------------------------------------------
      subroutine disp7P(neq,ZZ)
c
c  ***************************************************
c  THIS ROUTINE IS A MODIFICATION OF THE LINEAR PART OF ETAWN6
c  WRITTEN BY GLENN BATEMAN. THE MODIFICATIONS CONSIST OF THE INCLUSION
c  OF IMPURITIES IN THE DILUTION APPROX. IN THE SYSTEM WITH
c  4 EQUATIONS AND THE INCLUSION OF COLLISIONS ON TRAPPED
c  ELECTRONS IN THE FULL SYSTEM. ALSO PARALLEL ION MOTION IS INCLUDED
c  IN THE SYSTEM WITH 7 EQUATIONS.
c  MOST PARAMETERS ARE TRANSFERRED THROUGH COMMON BLOCKS LIKE GRAD,
c  IMP AND ETAWN6. NOTE THE INVERSE DEFINITION OF ZTAUH AND ZTAUZ !
c  ********************************************************
      implicit none
c
      INTEGER idp
      PARAMETER ( idp = 11 )
c
      LOGICAL inital, lmatv
      data inital /.true./, lmatv /.false./
c
      REAL*8  cetain(32)
     &  , omega(idp), gamma(idp)
      COMPLEX*16 ZZ(*)
c
      REAL*8 epsnhin, epsnzin, epstein, epsthin, epstzin, tauhin, tauzin
     & , fnzin, czin, azin, ftrapein, epsnin
     & , ekyrhoin, g, gzn, epsrat, etai, etae
     & , zb, si, eq, kiq, kxq, bt, bt1, vef, tvr, ftr
c
      INTEGER lprintin, lprint, neq, idim, ndim 
     & , ieq, j1, j2, j, iret
c
c ndim  = first dimension of the 2-D array difthi
c           and the maximum number of unstable modes allowed
c ieq   = number of equations
c
      REAL*8 zamr(idp,idp),zami(idp,idp),zbmr(idp,idp),zbmi(idp,idp)
     &  ,zamrt(idp,idp), zbmrt(idp,idp)
     &  ,zalfr(idp),zalfi(idp),zbeta(idp),zvr(idp,idp),zvi(idp,idp),ztol
c
      INTEGER iter(idp), ifail
c
c zamr(i,j) = matrix A
c zbmr(i,j) = matrix B
c   Note that the eigenvalues are
c omega(j) = ( zalfr(j) + i zalfi(j) ) / zbeta(j)
c where beta(j) will be 0.0 in the case of an infinite eigenvalue
c zvr(j) = eigenvector real part
c zvi(j) = eigenvector imag part
c
      REAL*8 wr, wi, wsq, H 
c
      REAL*8 zepsilon, zepsmach, zepsqrt
     & , zetah, zetaz, zetae
     & , zepsne, zepsnh, zepste, zepsth
     & , ztauh, ztauz, zep2nh, zep2nz, zep2ne, zft
     & , zimp, zfnz, zmass, zflh, zflz
      REAL*8 q,S,Cs,betae,em
c
c
      COMMON/GRAD/ epsnin,epsnhin,etai,etae
      COMMON/IMP/ fnzin,g,si,czin,tauzin,epsnzin,eq,kiq,kxq
      COMMON/ETAWN6/ cetain,epstein,
     &epsthin,epstzin,tauhin,azin,
     &ftrapein,ekyrhoin,lprintin,ndim
       COMMON/IRET/ iret
       COMMON/COLL/ bt,vef
       COMMON/BETAE/ betae, q, S, Cs,em
c
c ...  STRUCTURE OF MATRIX EQUATION ...
c
c ...  omega*zbmr(i,j) = zamr(i,j)
c
c    variables i=1,6: efi/Te, dTi/Ti, dni/ni, dTe/Te, dnq/nq, dTq/Tq
c    variables j=1,6 same as for i
c
c  ...................................................
c
c
c..initialize variables
c
      if ( inital ) then
c
        idim = idp
c
        tvr = 2./3.
        ftr = 5./3.
c
        zepsilon = 1.e-4
c
        zepsmach = 0.5
  2     if ( 0.5 * zepsmach + 1.0 .gt. 1.0 ) then
          zepsmach = 0.5 * zepsmach
          go to 2
        endif
c
        zepsqrt = sqrt ( zepsmach )
c
        inital = .false.
c
      endif
c
      lprint = lprintin
c
      ieq = max ( 2, neq )
c
      bt1 = bt - 2.5
c
      iret = 0
c..print header
c
      if ( lprint .gt. 2 ) then
c
        write (6,*)
        write (6,*)
     & 'Weiland-Nordman eigenvalue equations, subroutine etawn6'
        write (6,*) '(all frequencies normalized by omega_{De})'
        write (6,*) '(all diffusivities normalized by '
     &    ,'omega_{De} / k_y^2'
c
      endif
c
c..check validity of input data
c
      if ( neq .lt. 2 ) call abortb (6
     & ,' neq .lt. 2 in sbrtn etawn6')
c
c
        if ( czin .lt. 1.0 ) call abortb (6
     &   ,' czin .lt. 1.0 in sbrtn etawn6')
c
c      endif
c
      do 10 j1=1,ndim
        omega(j1) = 0.0
        gamma(j1) = 0.0
  10  continue
c
c
c..compute the rest of the dimensionless variables needed
c
      zetah  = etai
      zetae  = etae
c
c  *******  NOTE THE INVERSE DEFINITION OF ZTAUH ! ******
      ztauh  =1./ tauhin
c
      zep2nh=epsnhin
      zep2ne=epsnin
      zft    = ftrapein
      zflh   = ekyrhoin**2
c
      zimp   = czin
      zmass  = azin
c
        zfnz   = fnzin * zimp
        zetaz=epsnzin/epstzin
c
c  ******  NOTE THE INVERSE DEFINITION OF ZTAUZ ! ******
        ztauz = 1./(czin*tauzin)
        zep2nz=epsnzin
        zflz   = zmass * zflh / zimp**2
c
c
c..Note:
c
c  ztauz = T_Z / ( Z T_e )
c  zfnz  = Z n_Z / n_e
c  zimp  = Z
c  zmass = m_Z / m_H
c  zflh  = k_y^2 \rho_{sH}^2
c  zflz  = k_y^2 \rho_{sZ}^2
c
c..diagnostic output
c
      if ( lprint .gt. 6 ) then
        write (6,*)
        write (6,*) '--------------------------------------'
        write (6,*)
        write (6,*)
        write (6,*) zetah,' = zetah'
        write (6,*) zetaz,' = zetaz'
        write (6,*) zetae,' = zetae'
        write (6,*) ztauh,' = ztauh'
        write (6,*) ztauz,' = ztauz'
        write (6,*) zep2nh,' = zep2nh'
        write (6,*) zep2nz,' = zep2nz'
        write (6,*) zep2ne,' = zep2ne'
        write (6,*) zft,' = zft'
        write (6,*) zimp,' = zimp'
        write (6,*) zmass,' = zmass'
        write (6,*) zfnz,' = zfnz'
        write (6,*) zflh,' = zflh'
        write (6,*) zflz,' = zflz'
        write (6,*) zepsqrt,' = zepsqrt'
        write (6,*) zepsmach,' = zepsmach'
      endif
c
c..set matricies for eigenvalue equation
c
      do j1=1,idim
        zalfr(j1) = 0.0
        zalfi(j1) = 0.0
        zbeta(j1) = 0.0
        do j2=1,idim
          zamr(j1,j2) = 0.0
          zami(j1,j2) = 0.0
          zbmr(j1,j2) = 0.0
          zbmi(j1,j2) = 0.0
          zvr(j1,j2) = 0.0
          zvi(j1,j2) = 0.0
        enddo
      enddo
c
      if ( ieq .eq. 2 ) then
c
c..two equations when trapped particles and FLR effects omitted
c
        zamr(1,1) = ( 1.0 / zep2nh ) - ztauh - 1.0
        zamr(1,2) = - ztauh
        zamr(2,1) = ( zetah - tvr ) / zep2nh
        zamr(2,2) = - ztauh * ftr
c
        zbmr(1,1) = 1.0
        zbmr(1,2) = 0.0
        zbmr(2,1) = - tvr
        zbmr(2,2) = 1.0
c
      elseif ( ieq .eq. 4 ) then
c
c..4 equations with trapped electrons and FLR effects
c
c  equations for e phi/T_e, T_H, n_i, and T_e
c
      if  ( lprint .gt. 5 ) then
        write (6,*)
        write (6,*) ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
      endif
c
c       ion continuity
c
        zamr(1,1) = 1.0 - zep2nh - zflh * ztauh * ( 1.0 + zetah )
        zamr(1,2) = - zep2nh * ztauh
        zamr(1,3) = - zep2nh * ztauh
c
        zbmr(1,1) = zflh * zep2nh
        zbmr(1,3) = zep2nh
c
c  ion energy
c
        zamr(2,1) = zetah - tvr
        zamr(2,2) = - zep2nh * ztauh * ftr
c
        zbmr(2,2) =   zep2nh
        zbmr(2,3) = - zep2nh * tvr
c
c  trapped electron continuity
c
c   Calculates the total electron density perturbation and replaces it
c   by the ion density perturbation. The dilution factor 1-zfnz has now been 
c   added.
c
        zamr(3,1) = zft - zep2ne
        zamr(3,3) = zep2ne * (1. - zfnz)
        zamr(3,4) = zft * zep2ne
c
        zbmr(3,1) = ( zft - 1.0 ) * zep2ne
        zbmr(3,3) = zep2ne * (1. - zfnz)
c
c  trapped electron energy
c
        zamr(4,1) = zft * ( zetae - tvr )
        zamr(4,4) = zft * zep2ne * ftr
c
        zbmr(4,1) = ( 1.0 - zft ) * zep2ne * tvr
        zbmr(4,3) = - zep2ne * tvr
        zbmr(4,4) = zft * zep2ne
c
      else if ( ieq .eq. 6 ) then
c
c..Six equations with impurities, trapped electrons, and FLR
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, and T_Z
      if ( lprint .gt. 5 ) then
      write (6,*)
      write (6,*)
     & 'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
      endif
c
c  hydrogen density
c
        zamr(1,1) = - 1.0
     &   + ( 1.0 - zflh * ztauh * ( 1.0 + zetah ) ) / zep2nh
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = 1.0
c
c  hydrogen energy
c
        zamr(2,1) = ( zetah - tvr ) / zep2nh
        zamr(2,2) = - ztauh * ftr
c
        zbmr(2,2) = 1.0
        zbmr(2,3) = - tvr
c
c  trapped electron density
c
        zamr(3,1) = - 1.0 + zft / zep2ne
        zamr(3,3) = 1.0 - zfnz
        zamr(3,4) = zft
        zamr(3,5) = zfnz
c
        zbmr(3,1) = zft - 1.0
        zbmr(3,3) = 1.0 - zfnz
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft * ( zetae - tvr ) / zep2ne
        zamr(4,4) = zft * ftr
c
        zbmr(4,1) = ( 1.0 - zft ) * tvr
        zbmr(4,3) = - ( 1.0 - zfnz ) * tvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * tvr
c
c  impurity density
c
        zamr(5,1) = - 1.0
     &    + ( 1.0 - zflz * ztauz * ( 1.0 + zetaz ) ) / zep2nz
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
c
        zbmr(5,1) = zflz
        zbmr(5,5) = 1.0
c
c  impurity energy
c
        zamr(6,1) = ( zetaz - tvr ) / zep2nz
        zamr(6,6) = - ztauz * ftr
c
        zbmr(6,5) = - tvr
        zbmr(6,6) = 1.0
c
      else if ( ieq .eq. 7 ) then
c
c..Seven equations with impurities, trapped electrons, collisions and FLR
c  Parallel ion motion in strong ballooning approx.
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z and F
c  Here F is defined as F = GM*e phi/T_e where GM=1+etae/(epsn*(omega-1+i*vef))
c
c  ****
      H = 0.5*DABS(S)/q
c  ***
      if ( lprint .gt. 5 ) then
      write (6,*)
      write (6,*)
     & 'Seven eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z and F'
      endif
c
c  hydrogen density
c
        zamr(1,1) = - 1.0 
     &   + ( 1.0 - zflh * ztauh * ( 1.0 + zetah ) ) / zep2nh
        zami(1,1) = -H
        zamr(1,2) = - ztauh
        zami(1,2) = -ztauh*H
        zamr(1,3) = - ztauh
        zami(1,3) = -ztauh*H
c
        zbmr(1,1) = zflh
        zbmr(1,3) = 1.0
c
c  hydrogen energy
c
        zamr(2,1) = ( zetah - tvr ) / zep2nh
        zamr(2,2) = - ztauh * ftr
c
        zbmr(2,2) = 1.0
        zbmr(2,3) = - tvr
c
c  trapped electron density
c
        zamr(3,1) = - 1.0 + zft / zep2ne 
        zami(3,1) = vef*(1.-zft)
        zamr(3,3) = 1.0 - zfnz
        zami(3,3) = -vef*(1. - zfnz)
        zamr(3,4) = zft
        zamr(3,5) = zfnz
        zami(3,5) = -vef*zfnz
        zami(3,7) = vef*zft
c
        zbmr(3,1) = zft - 1.0
        zbmr(3,3) = 1.0 - zfnz
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft*(zetae - tvr)/zep2ne 
        zami(4,1) = vef*tvr*(bt-2.5*(1.-zft))
        zami(4,3) = -vef*tvr*bt1*(1.-zfnz)
        zamr(4,4) = zft * ftr
        zami(4,5) = -vef*tvr*bt1*zfnz
        zami(4,7) = -ftr*vef*zft
c
        zbmr(4,1) = ( 1.0 - zft ) *tvr
        zbmr(4,3) = - ( 1.0 - zfnz ) *tvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * tvr
c
c  impurity density
c
        zamr(5,1) = - 1.0
     &    + ( 1.0 - zflz * ztauz * ( 1.0 + zetaz ) ) / zep2nz
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
c
        zbmr(5,1) = zflz
        zbmr(5,5) = 1.0
c
c  impurity energy
c
        zamr(6,1) = ( zetaz - tvr ) / zep2nz
        zamr(6,6) = - ztauz * 5.0 / 3.0
c
        zbmr(6,5) = - tvr
        zbmr(6,6) = 1.0
c
c  variable F
c
         zamr(7,1) = zetae/zep2ne - 1. 
         zami(7,1) = vef
         zamr(7,7) = 1. 
         zami(7,7) = -vef
c
         zbmr(7,1) = -1.
         zbmr(7,7) = 1.
c
c
      else
c
        write (6,*)
        write (6,*) ieq,' = ieq in sbrtn etawn6'
        call abortb(6
     & ,'the value of ieq is wrong in sbrtn etawn6')
c
      endif
c
c..save copy of matrix which is over-written by sbrtn f02bjf
c
      do j2=1,ieq
        do j1=1,ieq
          zamrt(j1,j2) = zamr(j1,j2)
          zbmrt(j1,j2) = zbmr(j1,j2)
        enddo
      enddo
c
c..diagnostic output
c
      if ( lprint .gt. 6 ) then
        write (6,*)
        write (6,*) ' zamr(j1,j2)  j2 ->'
        do j1=1,ieq
          write (6,192) (zamr(j1,j2),j2=1,ieq)
        enddo
c
        write (6,*)
        write (6,*) ' zbmr(j1,j2)  j2->'
        do j1=1,ieq
          write (6,192) (zbmr(j1,j2),j2=1,ieq)
        enddo
 192  format (1p6e13.5)
      endif
c
c..find the eigenvalues using NAG14 routine f02bjf or f02gjf
c
      ztol = max ( 0.0, cetain(32) )
      ifail = 1
      lmatv = .false.
c
      if( ieq .eq. 7 ) go to 201
c
      call f02bjf ( ieq, zamr, idim, zbmr, idim, ztol
     & , zalfr, zalfi, zbeta, lmatv, zvr, idim, iter, ifail )
c
      go to 202
c
  201 continue
c
      call f02gjf ( ieq,zamr,idim,zami,idim,zbmr,idim,zbmi,idim,ztol
     & ,zalfr,zalfi,zbeta,lmatv,zvr,idim,zvi,idim,iter,ifail )
c
  202 continue
c
c      if ( ifail .gt. 0 ) call abortb ( 6
c     & ,'ifail .gt. 0 after call f02bjf in sbrtn etawn6')
c
       if( ifail .le. 0 ) goto 210
       iret = 1
c      write(*,205) zetah,zetaz,zetae,ztauh,ztauz,zep2nh,zep2nz
c     &,zep2ne,zft,zfnz
c  205 format(x,'zetah',G10.4,' zetaz',G10.4,' zetae',G10.4,' ztauh',
c     1G10.4,/,' ztauz',G10.4,' zep2nh',G10.4,' zep2nz',G10.4,' zep2ne'
c     1,G10.4,/,' zft',G10.4,' zfnz',G10.4)
      go to 215
c
  210 continue
c..compute the complex eigenvalues
c
      do j=1,ieq
c
        zb = zbeta(j)
        if ( abs(zbeta(j)) .lt. zepsilon ) zb = zepsilon
        omega(j) = zalfr(j) / zb
        gamma(j) = zalfi(j) / zb
        ZZ(j)=omega(j)+(0.D0,1.D0)*gamma(j)
c
       enddo
c
  220 format(2X,'wr=',G11.4,' wi=',G11.4)
c
  215 continue
c
        if ( lprint .gt. 6 ) then
        do j=1,ieq
        write(*,220) omega(j),gamma(j)
      enddo
      endif
c
      return
      end
C     File DIFFTD.F
C     Assumes maximum 7 unstable roots
C
      SUBROUTINE DIFF(RP,IR,TAUI,FT,U,CHI,CHE,D,CHQ,DHQ)
C
      IMPLICIT NONE
      INTEGER I,IR,IX,J,IW,LPRINTIN,NDIM
      REAL*8 WR,WI,EN,ENI,ENH,EI,EE,TAUI,FT,TE,N
      REAL*8 XI(5),XE(5),XD(5),CHI(5),CHE(5),D(5)
      REAL*8 CHQ(5),DHQ(5),XQ(5),XDQ(5)
      REAL*8 U(5,100)
      REAL*8 DNI,DNE,WSQ,NN
      REAL*8 A,B,C,E,F,A1,DH,WR1,WI1
      REAL*8 SCHI,SCHE,SD,HP,SHP,NR,NI
      REAL*8 THRD,TVR,FTR,STR,FTRT,RFL,D1,XHH,THRD2
      REAL*8 XIH,XEH,XDH,STF,PHS,NLX
      REAL*8 XDE,XDI,VEF,BTA,GAR,GAI,GBR,GBI,YDA
      REAL*8 HR,DIVGA,DIVGB,DIVA,DIVB,DEVGA,DEVGB,DEVA
      REAL*8 DEVB,SVA,SVB,WII,VEFN,BT1,EIH,EEH,EQH,IMP,TZ
      REAL*8 SCHQ,SDQ,KXQ,CHQEFF,DQEFF
      REAL*8 DTIQ,DTQQ,XQH,NQ,TQ,TIQ,NG,LNEH,LNHE
      REAL*8 BE,G,GI,SI,Z,TAUZ,EQ,ENQ,DTIMP,NIMP,H
      REAL*8 KIQ,DNQ,DI,DE,DQ,NQR,NQI,K,T,S,DQT,DQN
      REAL*8 WEXB,ROT
      REAL*8 DIVGA1,DIVGA2,DIVGA3,DIVA1,DIVA2,DIVA3
      REAL*8 DEVGA1,DEVGA2,DEVGA3,DEVA1,DEVA2,DEVA3
      REAL*8 DIVGAT,DIVAT,DEVGAT,DEVAT,SVAT,AT,CT,ET,HT,TT,DQNT
      REAL*8 CETAIN(32),ETE,ETI,ETQ,TAUH,AZ,FX,RFLT
      REAL*8 H1,H2,A2,A3,C1,C2,E1,E2,T1,T2,DQN1,DQN2
      REAL*8 SVA1,SVA2,SVA3
      REAL*8 PI,PE,PN,PQ,PNQ,PIQ,EA
      REAL*8 HPS,PES,PNS,PQS,PNQS
      REAL*8 GCI,GCE,GD,GCQ,GNQ,GK
      COMPLEX*16 RP(7),W,NE,BT2
      COMPLEX*16 GA,GB,GM,IU,HC
      COMMON/GRAD/ EN,ENH,EI,EE
      COMMON/PAR/ THRD,FTR,STR,RFL,D1,XIH,IX,IW
      COMMON/EFF DIFF/ SCHI,SCHE,SD,SCHQ,SDQ
      COMMON/COLL/ BTA,VEF
      COMMON/IMP/ BE,G,SI,Z,TAUZ,ENQ,EQ,KIQ,KXQ
      COMMON/WROT/ WEXB,ROT
      COMMON/ETAWN6/ CETAIN,ETE,ETI,ETQ,TAUH,AZ,FX,RFLT,
     1LPRINTIN,NDIM
      COMMON/SCALE/ NLX
      COMMON/TP/ EA,HPS,PES,PNS,PQS,PNQS
      COMMON/GRKVOT/ GK
C
      TVR=2.D0*THRD
      THRD2=THRD*THRD
      FTRT=FTR*TAUI
      IU=(0.D0,1.D0)
      SHP=0.D0
      ENI=1./EN
      EIH=EI-STR+FTR*ENH
      EEH=EE-STR+FTR*EN
      EQH=EQ-STR+FTR*ENQ
      LNEH=EN/ENH
      LNHE=ENH/EN
      G=1.-Z*BE
      NG=MAX(G,1.D-10)
      GI=1./NG
      TZ=Z*TAUZ
      IF(TZ.LT.0.0001) GO TO 00005
      IMP=1./TZ
00005 CONTINUE
C
      DO 00010 I=1,5
      CHI(I)=0.D0
      CHE(I)=0.D0
      D(I)=0.D0
      CHQ(I)=0.D0
00010 DHQ(I)=0.D0
C
      SCHI=0.D0
      SCHE=0.D0
      SD=0.D0
      SCHQ=0.D0
      SDQ=0.D0
      XHH=0.D0
C
C
      HPS=0.D0
      PES=0.D0
      PNS=0.D0
      PQS=0.
      PNQS=0.D0
C
      N=U(3,IX)
      TE=U(2,IX)
      TQ=U(4,IX)
      NQ=U(5,IX)
C
C  IR IS THE NUMBER OF UNSTABLE MODES
C  A SUMMATION IS PERFORMED OVER THE DIFFUSION
C  DUE TO ALL UNSTABLE MODES
C
      IF(IR.EQ.0) GOTO 00110
      DO 01100 J=1,IR
      W=RP(J)-IU*ROT*DABS(WEXB)
      WR=DREAL(W)
      WI=DIMAG(W)
      IF(WI.GT.0.001) GO TO 11 
      WI=0.001
      GO TO 01100
   11 CONTINUE
      W=CMPLX(WR,WI)
      GM=1.+EE*ENI/(W-1.D0+IU*VEF)
      BT1=BTA-2.5
      BT2=BTA-2.5*GM
      HC=W-FTR+TVR*BT1
      NE=W*W-2.D0*FTR*W+FTR+IU*VEF*HC
      WSQ=WR*WR+WI*WI
      WII=1./WI
C
C
C   ******   COLLISION PARAMETERS  *******
C
      GA=W-FTR+TVR*BT2
      GB=(W-FTR)/(W-1.D0+IU*VEF)
      GAR=DREAL(GA)
      GAI=DIMAG(GA)
      GBR=DREAL(GB)
      GBI=DIMAG(GB)
      HR=WR-FTR+TVR*BT1
      XDE=WSQ-FTR*WR
      XDI=WSQ+FTR*TAUI*WR
      YDA=WR*(1.D0-EN)+EE-STR
C   ***************************************
C
      IF(IW.NE.1) GOTO 00021
      WR1=DABS(EN)*WR
      WI1=DABS(EN)*WI
      WRITE(*,00020) WR1,WI1,J
00020 FORMAT(2X,'WR=',D11.5,' WI=',D11.5,' J=',I5)
00021 CONTINUE
C
      NR=DREAL(NE)
      NI=DIMAG(NE)
      NN=(NR)**2+(NI)**2
C
      IF(VEF.EQ.0.) GO TO 25
      DIVGA=FTRT*(EN*(GAR*NI-GAI*NR)-YDA*WI+WI*HR*(1.-
     1EN))
      DIVGB=FTRT*(GBR*NI-GBI*NR-WI)
      DIVA=XDI*(EN*(GAR*NR+GAI*NI)-WI*WI*(1.-EN)-YDA*HR)
      DIVB=XDI*(GBR*NR+GBI*NI-HR)
      DEVGA=BT1*((YDA-VEF*GAI*EN)*NI-(WI*(1.-EN)+VEF*GAR*EN
     1)*NR)+FTR*(WI*YDA+EN*(GAI*NR-GAR*NI)-WI*HR*(1.-EN))
      DEVGB=BT1*((1.-VEF*GBI)*NI-VEF*GBR*NR)+FTR*(WI
     1+GBI*NR-GBR*NI)
      DEVA=XDE*(EN*(GAR*NR+GAI*NI)-WI*WI*(1.-EN)-YDA*HR)
     1-BT1*(WR-FTR)*((YDA-VEF*GAI*EN)*NR+(WI*(1.-EN)
     1+VEF*GAR*EN)*NI)
      DEVB=XDE*(GBR*NR+GBI*NI-HR)-BT1*(WR-FTR)*((1.-
     1VEF*GBI)*NR+VEF*GBR*NI)
      SVA=EN*(GAR*NR+GAI*NI)-WI*WI*(1.-EN)-YDA*HR
      SVB=GBR*NR+GBI*NI-HR
      VEFN=VEF/NN
      DIVGA=DIVGA*VEFN
      DIVGB=DIVGB*VEFN
      DEVGA=DEVGA*VEFN
      DEVGB=DEVGB*VEFN
      DIVA=DIVA*VEFN
      DIVB=DIVB*VEFN
      DEVA=DEVA*VEFN
      DEVB=DEVB*VEFN
      SVA=SVA*VEFN
      SVB=SVB*VEFN
   25 CONTINUE
C
      STF=STR-FTR*EN
      A1=WSQ*(EN-1.D0)+2.D0*WR*STF
      A=WSQ*(A1+FTR*(STR*EN-11.D0*THRD-TAUI*(1.D0-
     1FTR*EN)))+FTR*FTR*TAUI*(2.D0*WR*(1.D0-EN)-STF)
      A=A/NN
      B=(WSQ*(2.D0*(WR-FTR)+FTR*TAUI)-
     1FTR*FTR*TAUI)/NN
      C=WSQ*(A1+TVR*FTR*(EN-4.D0))+
     1FTR*FTR*(2.D0*WR*(EN-1.D0)+STF)
      C=C/NN
      DH=(WSQ*(2.D0*WR-5.D0)+FTR*FTR)/NN
      E=WSQ*(1.D0-EN)-2.D0*WR*STF+FTR*(11.D0*THRD
     1-STR*EN)
      E=E/NN
      F=2.D0*(-WR+FTR)/NN
C
      NQR=WR**2-WI**2+2.D0*FTR*IMP*WR+FTR*IMP*IMP
      NQI=2.D0*WI*(WR+FTR*IMP)
      NIMP=NQR**2+NQI**2
C
C   ****   SPLITTING IN  EN  AND EE  ************
C
      DIVGA1=FTRT*((STR-WR)*WI+WI*HR)*VEFN
      DIVGA2=FTRT*(GAR*NI-GAI*NR+WR*WI-WI*HR)*VEFN
      DIVGA3=-FTRT*WI*VEFN
      DIVA1=XDI*(-WI*WI+(STR-WR)*HR)*VEFN
      DIVA2=XDI*(GAR*NR+GAI*NI+WI*WI+WR*HR)*VEFN
      DIVA3=-XDI*HR*VEFN
      DEVGA1=(BT1*((WR-STR)*NI-WI*NR)+FTR*(WI*(WR-STR)-WI*HR))
     &*VEFN
      DEVGA2=(BT1*((-WR-VEF*GAI)*NI+(WI-VEF*GAR)*NR)
     &+FTR*(-WI*WR+GAI*NR-GAR*NI+WI*HR))*VEFN
       DEVGA3=(BT1*NI+FTR*WI)*VEFN
      DEVA1=(XDE*(-WI*WI-(WR-STR)*HR)-BT1*(WR-FTR)*((WR-STR)*NR
     &+WI*NI))*VEFN
      DEVA2=(XDE*(GAR*NR+GAI*NI+WI*WI+WR*HR)-BT1*(WR-FTR)
     &*((-WR-VEF*GAI)*NR-(WI-VEF*GAR)*NI))*VEFN
      DEVA3=(-XDE*HR-BT1*(WR-FTR)*NR)*VEFN
      SVA1=(-WI*WI-(WR-STR)*HR)*VEFN
      SVA2=(GAR*NR+GAI*NI+WI*WI+WR*HR)*VEFN
      SVA3=-HR*VEFN
      A2=(WSQ*WSQ+FTR*((STR+FTR*TAUI-2.D0*WR)*WSQ+FTR*TAUI*(FTR-
     &2.D0*WR)))/NN
      A3=(WSQ*(-WSQ+2.D0*STR*WR-FTR*(11.D0*THRD+TAUI))+FTR*FTR
     &*TAUI*(2.D0*WR-STR))/NN
      C1=(NN-WSQ*WSQ+14.D0*THRD*WSQ*WR-40.D0*THRD2*WSQ
     &-50.D0*THRD2*WR+175.D0/27.D0)/NN
      C2=(WSQ*WSQ-10.D0*THRD*WSQ*WR+10.D0*THRD2*WSQ
     &+50.D0*THRD2*WR-125.D0/27.D0)/NN
      E1=(WSQ+11.D0*THRD*FTR-2.D0*WR*STR)/NN
      E2=-(-2.D0*WR*FTR+(WSQ+STR*FTR))/NN
      H1=(WSQ*(-WSQ-2.D0*IMP*WR*STR
     &+FTR*IMP*IMP*(-11.D0*THRD)+FTR*TAUI*IMP)
     &+FTR*FTR*TAUI*IMP*IMP*(2.D0*WR+STR*IMP))/NIMP
      H2=(WSQ*(WSQ+2.D0*IMP*WR*FTR+FTR*IMP*IMP*STR-FTR*TAUI*IMP*FTR)
     &-FTR*FTR*TAUI*IMP*IMP*(2.D0*WR+FTR*IMP))/NIMP
      T1=KXQ*(WSQ*(-WSQ-2.D0*IMP*WR*STR-FTR*IMP*IMP*8.D0*THRD)
     &+FTR*FTR*IMP**3*(2.D0*WR+IMP*STR))/NIMP
      T2=KXQ*(WSQ*(WSQ+2.D0*IMP*WR*FTR+FTR*IMP*IMP*TVR)
     &-FTR*FTR*IMP**3*(2.d0*WR+IMP*FTR))/NIMP
      DQN1=(-NQR+2.D0*(WR+IMP*STR)*(WR+FTR*IMP))/NIMP
      DQN2=(NQR-2.D0*(WR+IMP*FTR)*(WR+FTR*IMP))/NIMP
C
      DIVGAT=DIVGA1+EN*DIVGA2+EE*DIVGA3
      DIVAT=DIVA1+EN*DIVA2+EE*DIVA3
      DEVGAT=DEVGA1+EN*DEVGA2+EE*DEVGA3
      DEVAT=DEVA1+EN*DEVA2+EE*DEVA3
      SVAT=SVA1+EN*SVA2+EE*SVA3
      AT=EN*A2+A3
      CT=C1-1.+EN*C2
      ET=E1+EN*E2
      HT=H1+ENQ*H2
      TT=T1+ENQ*T2
      DQNT=DQN1+EN*DQN2
C
      IF(IW.NE.1) GOTO 00024
C      WRITE(*,22) NN,A,B,C,DH,E,F
00022 FORMAT(2X,7G11.3)
00024 CONTINUE
C
C **** IMPURITIES *****
C
      NIMP=(WR*(WR+2.*FTR*IMP)-WI*WI+FTR*IMP*IMP)**2
     1+4.*WI*WI*(WR+FTR*IMP)**2
      DTIMP=(WSQ*(WSQ*(ENQ-1.)+2.*IMP*WR*EQH+FTR*IMP
     1*IMP*(2.*EQ-11./3.+STR*ENQ)+FTR*TAUI*IMP*(1.+
     1+EQ-FTR*ENQ))+FTR*TAUI*(2.*FTR*IMP*IMP*WR*(1.-ENQ)
     1-FTR*IMP**3*EQH))/NIMP
C  *************
C
      DNI=(WR+FTR*TAUI)**2+WI*WI
      DNE=(WR-FTR)**2+WI*WI
      DNQ=(WR+FTR*IMP)**2+WI*WI
      DI=EI*RFL*DNI
      DE=EE*RFL*DNE
      DQ=EQ*RFL*DNQ
C
C **** IMPURITIES *****
C
      NIMP=NQR**2+NQI**2
C
      H=(WSQ*(WSQ*(ENQ-1.)-2.*IMP*WR*(STR-FTR*ENQ)
     1+FTR*IMP*IMP*(-11.*THRD+STR*ENQ)+FTR*TAUI*IMP*(1.-FTR
     1*ENQ))+FTR*FTR*TAUI*IMP*IMP*(2.*WR*(1.-ENQ)+(STR-FTR
     1*ENQ)*IMP))/NIMP
      K=IMP*(FTR*FTR*TAUI*IMP*IMP-WSQ*(2.*WR
     1+FTR*(2.*IMP+TAUI)))/NIMP
C
C  *************
      T=(WSQ*(WSQ*(ENQ-1.)-2.*IMP*WR*(STR-FTR*ENQ)
     1+FTR*IMP*IMP*(-8.*THRD+TVR*ENQ))+FTR*FTR*(IMP)**3
     1*(2.*WR*(1.-ENQ)+IMP*(STR-FTR*ENQ)))/NIMP
      S=IMP*(FTR*FTR*(IMP)**3-WSQ*(2.*WR+5.*IMP))/NIMP
C
      T=KXQ*T
      S=KXQ*S
      DQN=(ENQ-1.)*NQR+2.*((1.-ENQ)*WR+IMP*(STR-FTR*ENQ)
     1)*(WR+FTR*IMP)
      DQT=2.*IMP*(WR+FTR*IMP)
      DQN=DQN/NIMP
      DQT=DQT/NIMP
C
      DTIQ=H-EQ*K
      DTQQ=T-EQ*S
C
      PHS=(WR-FTR)*DREAL(BT2)+WI*DIMAG(BT2)
      XDH=D1*(TE**1.5D0)*WI**3/RFL
      XIH=XDH/DNI
      XEH=XDH/DNE
      XQH=XDH/DNQ
C
C
      XI(1)=XIH
      XI(2)=TVR*FT*XIH*TAUI*GI*(B-DIVGA3-DIVGB-WII*(DIVA3+DIVB))
      XI(3)=-TVR*XIH*(TE/N)*TAUI*(LNEH+FT*GI*(A3+DIVGA1+DIVA1*WII))
      XI(4)=-XIH*TVR*BE*Z*GI*K*TAUZ*TAUI
      XI(5)=TVR*XIH*(TE/N)*TAUI*Z*H1*GI
C
      PI=FT*TVR*XIH*GI*(A2+DIVGA2+WII*DIVA2)*EA
      PIQ=-TVR*XIH*Z*H2*BE*GI*EA
C
C
      XE(1)=0.D0
      XE(2)=FT*XEH*(1.D0+TVR*(DH-DEVGB-DEVGA3-WII*(DEVB+DEVA3)))
      XE(3)=-TVR*FT*XEH*(C1+DEVGA1+WII*DEVA1)*TE/N
      XE(4)=0.D0
      XE(5)=0.D0
C
      PE=FT*TVR*XEH*(C2+DEVGA2+WII*(DEVA2+VEF*PHS))*EA
C
C
      XD(1)=0.D0
      XD(2)=-FT*XDH*(N/TE)*(F+WII*(SVB+SVA3))
      XD(3)=FT*XDH*(E1-WII*SVA1)
      XD(4)=0.D0
      XD(5)=0.D0
C
      PN=-FT*XDH*(E2-WII*SVA2)*EA
C
C
      XQ(1)=0.D0
      XQ(2)=0.D0
      XQ(3)=0.D0
      XQ(4)=(1.D0+TVR*S)*XQH
      XQ(5)=-TVR*TQ*(1.D0+T1)*XQH/NQ
C
      PQ=TVR*XQH*T2*EA
C
C
      XDQ(1)=0.
      XDQ(2)=0.
      XDQ(3)=0.
      XDQ(4)=-XDH*DQT*NQ/TQ
      XDQ(5)=XDH*DQN1
C
      PNQ=-XDH*DQN2*EA
C
      TIQ=1.D0/(TAUZ*KIQ)
C
      HP=XIH*GI*TAUI*TVR*FTR*(1.D0-FT)*EA
      GCI=XI(1)+(LNHE*(XI(2)*EE+XI(3)*N/TE)+XI(4)*EQ*TIQ
     1+XI(5)*NQ/(TE*KIQ))/(EI*TAUI)-(HP+PI+PIQ)*GK*ETI/EA
      GCE=XE(2)+XE(3)*N/(TE*EE)-PE*GK*ETE/EA
      GD=XD(2)*TE*EE/N+XD(3)-PN*GK*EN/EA
      GCQ=XQ(4)+XQ(5)*NQ/(TQ*EQ)-PQ*GK*ETQ/EA
      GNQ=XDQ(4)*TQ*EQ/NQ+XDQ(5)-PNQ*GK*ENQ/EA
C
C      WRITE(*,00092) HP
00092 FORMAT(2X,'HP=',F11.5)
C
      SHP=SHP+HP
C
      HPS=HPS+HP+PI+PIQ
      PES=PES+PE
      PNS=PNS+PN
      PQS=PQS+PQ
      PNQS=PNQS+PNQ
C
      SCHI=SCHI+GCI
      SCHE=SCHE+GCE
      SD=SD+GD
      SCHQ=SCHQ+GCQ
      SDQ=SDQ+GNQ
C
      IF(IW.NE.1) GOTO 00111
      WRITE(*,00094) SHP
00094 FORMAT(2X,'SHP=',G11.3)
      WRITE(*,00096) XDH,XIH,XEH,XQH
00096 FORMAT(2X,'XDH=',G11.3,' XIH=',G11.3,' XEH=',G11.3,' XQH=',G11.3)
00111 CONTINUE
C
      CHQEFF=D1*TE**1.5*WI**3*(EQ-TVR-TVR*DTQQ)/DQ
      DQEFF=XDH*(DQN-EQ*DQT)
00222 FORMAT(2X,'CHQEFF=',G11.4,' DQEFF=',G11.4)
C
      IF(IW.NE.1) GO TO 00093
      WRITE(*,00222) CHQEFF,DQEFF
      WRITE(*,00099) SCHI,SCHE,SD,SCHQ,SDQ,J
00099 FORMAT(2X,'SCHI=',F11.5,' SCHE=',F11.5,' SD=',
     1F11.5,' SCHQ=',F11.5,' SDQ=',F11.5,' J=',I5)
00093 CONTINUE
C
      XHH=XHH+XIH
C
      DO 00100 I=1,5
      CHI(I)=CHI(I)+XI(I)
      CHE(I)=CHE(I)+XE(I)
      D(I)=D(I)+XD(I)
      CHQ(I)=CHQ(I)+XQ(I)
      DHQ(I)=DHQ(I)+XDQ(I)
C
00100 CONTINUE
01100 CONTINUE
C
C
C      TPX=TP*XIH*GI
C      HPH=HPS/(XIH*GI)
C
      IF(IW.NE.1) GO TO 00110
C
C      WRITE(*,211) TPX,HPS,TP,HPH
C  211 FORMAT(/,2X,'difftd',' TPX=',G12.4,' HPH=',G12.4,' TP=',G12.4
C     &,' HPH=',G12.4)
00110 CONTINUE
C
      RETURN
      END
E
      0.378
N0
      4.76
T0I
      0.56
T0E
      0.41
TCI
      4.83
TCE
      4.86
TM
      0.5
PHI
      7.4
PHE
      3.6
SAN
      0.02
W      
      3.71
RIW
      15
TAU
      1
FL
      0.1
XA
      0.1
XB
      0.9
XSI
      0.75
XSE
      0.75
XSN
      0.75
RIT
      31.
CX
      1.
TWR
      0.499999999995
SEI
      0.001
SEE
      0.177
SEN
      5.
KT
      0.00001D0
DX
      1.
PR
      1.17
ZE
      6.
KNC
      1.
COL
      1.
RIB
      2.
KP
      1.4
TEQ
      1.
TPL
      0.01
BF
      0.038
EV
      0.005
T0Q
      0.56
TCQ
      4.83
RNEQ
      7.
SEQ
      4.
XL
      0.7
RAD
      1.
BTOR
      3.1
ROT
      0.
