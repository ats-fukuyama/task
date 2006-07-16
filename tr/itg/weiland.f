C***********************************************************************
C  Linear solver of dispersion equation in matrix form
C**********************************************************************
      subroutine disp9t(neq,ZZ)
c
c  ***************************************************
c  THIS ROUTINE IS A MODIFICATION OF THE LINEAR PART OF ETAWN6
c  WRITTEN BY GLENN BATEMAN. THE MODIFICATIONS CONSIST OF THE INCLUSION
c  OF IMPURITIES IN THE DILUTION APPROX. IN THE SYSTEM WITH
c  4 EQUATIONS AND THE INCLUSION OF COLLISIONS ON TRAPPED
c  ELECTRONS IN THE FULL SYSTEM. THIS SYSTEM THEN CONSISTS
c  OF 7 EQUATIONS. WHEN PARALLEL ION MOTION IS INCLUDED THERE
c  ARE 8 EQUATIONS WITH COLLISIONS. IN disp10 ALSO ELECTROMAGNETIC 
c  EFFECTS ARE INCLUDED.  
c  MOST PARAMETERS ARE TRANSFERRED THROUGH COMMON BLOCKS LIKE GRAD,
c  IMP AND ETAWN6. NOTE THE INVERSE DEFINITION OF ZTAUH AND ZTAUZ !
c  ********************************************************
      implicit none
c      SAVE
c
      INTEGER idp
c------------------------------------------------------------------
c  Note  that idp gives the dimensions of zvr and zvi which also have to be
c  defined in the main programme and in difftd.f. These declarations thus
c  have to be  the same.
      PARAMETER ( idp = 10)
c-----------------------------------------------------------------------
c
      LOGICAL inital, lmatv
      data inital /.true./, lmatv /.true./
c
      REAL*8  cetain(32)
     &  , omega(idp), gamma(idp)
      COMPLEX*16 ZZ(*), ZZS(10)
c
      REAL*8 epsnhin, epsnzin, epstein, epsthin, epstzin, tauhin, tauzin
     & , fnzin, czin, azin, ftrapein, epsnin
     & , ekyrhoin, g,   etai, etae
     & , zb, si, eq, kiq, kxq, bt, bt1, vef, tvr, ftr, ftrt
c
      REAL*8 KAPPA, RAV, GAV,ALA,SH,SH2,WZR,WZI,WZIMAX
      REAL*8 H1,XH,EN,EI,EE,TAU,TAUI,FL,FT,GM,BTA,XT,R,HQR,HQI,WM
      COMPLEX*16 ALPC,ALPK,ALPHA,HQ,WZ,WZ1,WZ2,IU,H2,E,E1,WZJ(100),WZH
      COMPLEX*16 WZP,WZPP,DWN,DEW,DEWP,D2EW,W1,W2,AM1,BM1,DR
      COMPLEX*16 DWN1,WZA,WS,WSF,WSFF,WZO,WZE
      REAL*8 DW1,DW2,CTEST
c
      SAVE WZJ,idim,tvr,ftr,em1,A,zepsilon,zepsmach
      INTEGER lprintin, lprint, neq, idim, ndim 
     & , ieq, j1, j2, j, iret, IK, IST, IM, ISB,IMX
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
      REAL*8 wr,wi,H,fft,fzft
c
      REAL*8 zepsilon, zepsmach, zepsqrt
     & , zetah, zetaz, zetae 
     & , zepsnh, zepste, zepsth
     & , ztauh, ztauz, zep2nh, zep2nz, zep2ne, zft
     & , zimp, zfnz, zmass, zflh, zflz,zfs,A,em,em1
      REAL*8 q,S,Cs,alp,k1,k2,kps,betae,eni,alf,kpc
      INTEGER ITC,ITL,ITS,ITERA,IMET,ISEARCH,SEARCHMODE
      REAL*8 TOL
c
c
      COMMON/GRAD/ epsnin,epsnhin,etai,etae
      COMMON/IMP/ fnzin,g,si,czin,tauzin,epsnzin,eq,kiq,kxq,zfs
      COMMON/ETAWN6/ cetain,epstein,
     &epsthin,epstzin,tauhin,azin,
     &ftrapein,ekyrhoin,lprintin,ndim
      COMMON/IRET/ iret
      COMMON/COLL/ bt,vef
      COMMON/BETAE/ betae, q, S, Cs, em
      COMMON/TEST/ alp,alf,kpc,kps
      COMMON/HQ/ H1,HQ,GAV,WZ,WZP,WZIMAX
      COMMON/SHAFRS/ ALA
      COMMON/KAPPA/ KAPPA,RAV
      COMMON/IK/ IK,IST,ITC,ITL,ITS,ITERA,TOL
      COMMON/ZV/ zvr,zvi
C      COMMON/FZ/ fft,faft,ftt,ftf,fzf,fh
      COMMON/ISB/ ISB
      COMMON/SEARCHMODE/ SEARCHMODE
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
c  ---------------------------------------------------
c  Iteration control variables used only when neq=9.
c
c   ITC=1  means  iterate ITC=2 means use only average
c   ITS=1  means  iterations converged
c   ITL    means  maximum number of iterations
c   IST=1 means  start iterations from analytical approx.
c   ITER   means  number of iterations performed
c   TOL    means  relative error allowed for convergence
c   IF IST.NE.1 the iterations will start from whatever value stored in 
c   WZJ(IK). This will during the simulation be the root from the previous time
c   step with the largest growthrate corresponding to propagation in the ion 
c   drift direction.
c   It is necessary that IST=1 in the first time step
c  --------------------------------------------------------
c
c..initialize variables
c
      if ( inital ) then
c
        idim = idp
c
      if(neq.LE.idp) go to 40001
      write(*,40002)
40002 format(' The dimension idp is not sufficient') 
      go to 00085
40001 continue
        tvr = 2./3.
        ftr = 5./3.
        ftrt = ftr/tauhin
        em1=em
c        em1=1.D0
c
        A = azin
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
      ftrt=ftr/tauhin
      A = azin
      em1 = em
      ieq = max ( 2, neq )
c
      IU=(0.,1.)
c
      ftrt = ftr/tauhin
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
c      if ( abs(epstein) .lt. zepsqrt ) call abortb (6
c     & ,' abs(epstein) .lt. zepsqrt in sbrtn etawn6')
c
c      if ( abs(epsthin) .lt. zepsqrt ) call abortb (6
c     & ,' abs(epsthin) .lt. zepsqrt in sbrtn etawn6')
c
c      if ( neq .gt. 4 ) then
c
c        if ( abs(epstzin) .lt. zepsqrt ) call abortb (6
c     &   ,' abs(epstzin) .lt. zepsqrt in sbrtn etawn6')
c
c        if ( abs(epsnzin) .lt. zepsqrt ) call abortb (6
c     &   ,' abs(epsnzin) .lt. zepsqrt in sbrtn etawn6')
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
      zepsth = epsthin
      zepsnh = epsnhin 
      zepste = epstein 
c
c..compute the rest of the dimensionless variables needed
c
      zetah  = zepsnh / zepsth
*USCH      zetae  = zepsne / zepste
c
c  *******  NOTE THE INVERSE DEFINITION OF ZTAUH ! ******
      ztauh  =1./tauhin
c
      zep2nh=epsnhin
      zep2ne=epsnin
      zetae = zep2ne/zepste
      zft    = ftrapein
      zflh   = ekyrhoin**2
      eni = 1./zep2ne
c
      zimp   = czin
      zmass  = azin
c
        zfnz   = fnzin * zimp
c        zetaz  = zepsnz / zepstz
        zetaz=epsnzin/epstzin
c
c  ******  NOTE THE INVERSE DEFINITION OF ZTAUZ ! ******
        ztauz = 1./(czin*tauzin)
        zep2nz=epsnzin
c        zep2nz = 2.0 * zepsnz
        zflz   = zmass * zflh / zimp**2
c
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
c------------------------------------------------
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
        zamr(3,3) = zep2ne * (1. - zfnz - zfs)
        zamr(3,4) = zft * zep2ne
c
        zbmr(3,1) = ( zft - 1.0 ) * zep2ne
        zbmr(3,3) = zep2ne * (1. - zfnz -zfs)
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
        zamr(3,3) = 1.0 - zfnz -zfs
        zamr(3,4) = zft
        zamr(3,5) = zfnz
c
        zbmr(3,1) = zft - 1.0
        zbmr(3,3) = 1.0 - zfnz -zfs
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft * ( zetae - tvr ) / zep2ne
        zamr(4,4) = zft * ftr
c
        zbmr(4,1) = ( 1.0 - zft ) * tvr
        zbmr(4,3) = - ( 1.0 - zfnz -zfs) * tvr
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
c..Seven equations with impurities, trapped electrons, parallel ion motion
c and FLR
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z and F
c  Here F is defined as F = GM*e phi/T_e where GM=1+etae/(epsn*(omega-1+i*vef))
c
      if ( lprint .gt. 5 ) then
      write (6,*)
      write (6,*)
     & 'Seven eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z and Vp'
      endif
c
c********
      H=0.5*DABS(S)/q
c********
c
c  hydrogen density
c
      zamr(1,1) = -1.
     & + (1. - zflh * ztauh * (1. + zetah ))/zep2nh 
      zami(1,1) = -H
      zamr(1,2) = - ztauh
      zami(1,2) = -ztauh*H
      zamr(1,3) = - ztauh
      zami(1,3) = -ztauh*H
C
      zbmr(1,1) = zflh
      zbmr(1,3) = 1.
c
c hydrogen energy
c
      zamr(2,1) = (zetah - tvr )/zep2nh
      zamr(2,2) = - ztauh * ftr
c
      zbmr(2,2) = 1.
      zbmr(2,3) = -tvr
c
c trapped electron density
c
      zamr(3,1) = -1. + zft/zep2ne
      zami(3,1) = vef*(1. - zft)
      zamr(3,3) = 1. - zfnz -zfs
      zami(3,3) = -vef*(1. -zfnz - zfs)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = -vef*zfnz
      zami(3,7) = vef*zft
c
      zbmr(3,1) = zft - 1.
      zbmr(3,3) = 1. - zfnz -zfs
      zbmr(3,5) = zfnz
c
c trapped electron energy
c
      zamr(4,1) = zft * (zetae - tvr )/zep2ne
      zami(4,1) = vef*tvr*(bt-2.5*(1.-zft))
      zami(4,3) = -vef*tvr*bt1*(1. -zfnz - zfs)
      zamr(4,4) = zft * ftr
      zami(4,5) = -vef*tvr*bt1*zfnz
      zami(4,7) = -ftr*vef*zft
c
      zbmr(4,1) = (1. - zft) * tvr
      zbmr(4,3) = - (1. - zfnz -zfs) * tvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * tvr
c
c impurity density
c
      zamr(5,1) = -1.
     & + ( 1. - zflz * ztauz * (1. + zetaz ))/zep2nz
      zami(5,1)=-zimp*H/A
      zamr(5,5) = - ztauz
      zami(5,5) = -zimp*ztauz*H/A
      zamr(5,6) = - ztauz
      zami(5,6) = -zimp*ztauz*H/A
c
      zbmr(5,1) = zflz
      zbmr(5,5) = 1.
c
c impurity energy
c
      zamr(6,1) = (zetaz - tvr)/zep2nz
      zamr(6,6) = - ztauz * ftr
c
      zbmr(6,5) = -tvr
      zbmr(6,6) = 1.
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
      else if ( ieq .eq. 8 ) then
c
c..Eight equations with impurities, trapped electrons, parallel ion
c  motion collisions and FLR 
c  Equations for e phi/Te, Ti, ni, Tet, n_z, T_z, F, Vp 
c  Here F = omega*gamma (collisions).
c
      k1 = 0.25*q*q*zflh*SQRT(DABS((1.+zetah)*ztauh/((1.-zft)*zep2nh)))
      k2 = q*q*zflh*zflh*(1.+zetah)*ztauh/((1.-zft)*zep2nh)
      alp = 0.5*(k1+SQRT(k1+S*S*k2))
      alf = alp/(2.D0*zflh*q*q*betae*(1.D0 - zft))
      kps = 0.5*SQRT(alp/zflh)/q
      kpc = 1.D0
c
c
c hydrogen density
c
        zamr(1,1) = - 1.0
     &   + ( 1.0 - zflh * ztauh * ( 1.0 + zetah ) ) / zep2nh
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
        zamr(1,8) = kps
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
c  total electron density expressed in ion and imp densities
c
        zamr(3,1) = - 1.0 + zft / zep2ne 
        zami(3,1) = vef*(1.-zft)
        zamr(3,3) = 1.0 - zfnz -zfs
        zami(3,3) = -vef*(1. - zfnz -zfs)
        zamr(3,4) = zft
        zamr(3,5) = zfnz
        zami(3,5) = -vef*zfnz
        zami(3,7) = vef*zft
c
        zbmr(3,1) = zft - 1.0
        zbmr(3,3) = 1.0 - zfnz -zfs
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft*(zetae - tvr)/zep2ne 
        zami(4,1) = vef*tvr*(bt-2.5*(1.-zft))
        zami(4,3) = -vef*tvr*bt1*(1.-zfnz -zfs)
        zamr(4,4) = zft * ftr
        zami(4,5) = -vef*tvr*bt1*zfnz
        zami(4,7) = -ftr*vef*zft
c
        zbmr(4,1) = ( 1.0 - zft ) *tvr
        zbmr(4,3) = - ( 1.0 - zfnz -zfs ) *tvr
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
c     Parallel ion motion Vpi/Cs
c
         zamr(8,1) = kps
         zamr(8,2) = kps*ztauh
         zamr(8,3) = kps*ztauh
c
         zbmr(8,8) = 1.
c
      ELSE
      GO TO 08886
      ENDIF
      GOTO 08888
08886 CONTINUE
c
c--------------------------------------------------------------------------------
C      else if ( ieq .eq. 9 ) then
      IF(IEQ.NE.9) GO TO 08887
c
c..Nine  equations with impurities, trapped electrons, parallel ion
c  motion, collisions,  FLR , finite beta and parallel motion of impurities
c
c  Equations for e phi/Te, Ti, ni, Tet, n_z, T_z, F, Av, K
c  Here F = omega*gamma (collisions), K = omega*Av and Av is the parallel
c   magnetic vector potential.
c
c      WRITE(*,20061) WZJ(9),IK
20061 FORMAT(' Enter 9 Eq:s WZJ(9)=',2G11.3,' IK=',I13)

      EN=zep2nh
      ENI=1./EN
      EI=zetah
      EE=zetae
      TAU=1./ztauh
      TAUI=ztauh
      FL=zflh
      FT=zft
c
      GM=1./(1.-FT)
      BTA=FTRT*(GM+TAUI)
      H1=4.*q*q*FL
      XT=1./(1.+TAUI)
c
      ALA=em*q*q*betae*(1.+EE+TAUI*(1.+EI))/EN
      SH2=2.*S-1.+(KAPPA*(S-1.))**2
      SH=SQRT(SH2)
      H=0.5*ABS(SH)/q
      H2=IU*H
      GAV=1.
      ITS=0
      ITERA=1
c
      ISEARCH=1
c      SEARCHMODE=2
c  If ISB=1 we use the strong ballooning approx (GAV=1)
c     IST is 1 at the first step. Then we use an analytical approx.
      IF(ISB.EQ.1) GO TO 20001
      IF(IST.NE.1) GO TO 800
c********************************************************************
c    Here the analytical approximation for WZ is calculated 
c    and stored in WZJ(IK)
c
      E1=FTRT*(1.+FL)-ENI+FL*TAUI*ENI*(1.+EI)+(GM+FTRT)*GAV
     &+H2*(1.+FTRT)
      E1=0.5*E1/(1.+FL)
      E=(TAUI*ENI*GM*(EI-TVR)+BTA)*(GAV+H2)-FTRT*ENI*(1.-FL*TAUI*(1.
     &+EI))
      E=E/(1.+FL)
      WZ1=-E1+SQRT(E1*E1-E)
      WZ2=-E1-SQRT(E1*E1-E)
      WZ=WZ1
      IF(DIMAG(WZ2).GT.DIMAG(WZ1)) WZ=WZ2
      WZI=DIMAG(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI)
C      WRITE(*,10021)EI,E,EN,ENI,GM,BTA,GAV,FTRT,H2
10021 FORMAT(' EI=',G11.4,' E=',2G11.4,' EN=',G11.4,' ENI=',G11.4,
     &' GM=',G11.4,' BTA=',G11.4,' GAV=',G11.4,' FTRT=',G11.4,
     &' H2=',2G11.4)
c
      WZJ(IK)=WZ
      WZH=EN*WZ
      IF(LPRINT.EQ.2)   WRITE(*,10001) WZH
10001 FORMAT(2X,' DISP9T WZ=',2G11.3)
c*********************************************************************
C--THE ITERATIVE PROCEEDURE FOR WZ STARTS HERE --
c
      IMET=0
  800 CONTINUE
      WZ=WZJ(IK)
c----------------------------------------------------
      IF(lprint.NE.2) GO TO 30006
      WRITE(*,00798) WZ,ITERA,ISEARCH,IM
00798 FORMAT(' WZ=',2G11.3,' ITERA=',I5,' ISEARCH=',I5,'   IM=',I5)
      IF(IMET.EQ.1) GO TO 30001
      IF(IMET.EQ.2) GO TO 30002
      WRITE(*,30003)
30003 FORMAT(' Average WZ')
      GO TO 30006
30001 WRITE(*,30004) 
30004 FORMAT(' Newton-Raps')
      GO TO 30006
30002 WRITE(*,30005)
30005 FORMAT(' Muller')
c
c--------------------------------------------------
30006 CONTINUE
c
      WZP=WZ
      ALPK=0.5*SH*SQRT(H1*XT*FL*(1.+TAUI*(1.+EI)/(EN*WZ)))
      IF(DREAL(ALPK).GE.0.) GOTO 801
      ALPK=-ALPK
  801 CONTINUE
      ALPC=-IU*ALPK
      ALPHA=-IU*ABS(SH)*q*FL
      XH=ABS(ALPHA/ALPC)
      ALPC=XH*ALPC
      R=2.*ABS(DREAL(WZ*ALPC))
      IF(R.LT.0.001) R=0.001    !! NEW 01.03.8
      HQ=2.*ALPC/H1
  802 GAV=(1.+0.5*S/R)*EXP(-0.25/R)
      GAV=GAV-0.5*ALA*(1.-EXP(-1./R))
      IF(GAV.LT.0.01) GAV=0.01
c
20001 CONTINUE
      IF(ISB.NE.1) GO TO 20002
      k1=0.25*q*q*zflh*DSQRT(DABS((1.+zetah)*ztauh/((1.-zft)*zep2nh)))
      k2=q*q*zflh*zflh*(1.+zetah)*ztauh/((1.-zft)*zep2nh)
      R=k1+SQRT(ABS(k1+S*S*k2))
      H=0.5*DABS(SH)/q
      HQ=-IU*H
      ITC=2
20002 CONTINUE
      alp=0.5*R
      IF(alp.LT.0.1) alp=0.1
      alf = alp/(2.D0*zflh*q*q*betae*(1.D0 - zft))
      kps = 0.5*SQRT(alp/zflh)/q
      kpc = 1.D0
      RAV=1.+0.25*SH2/alp
c      WRITE(*,10002) alp,SH2,RAV
10002 FORMAT(2X,'alp=',G11.3,' SH2=',G11.3,' RAV=',G11.3)
c      WRITE(*,10003) XH,GAV,alf
10003 FORMAT(2X,'XH=',G11.3,' GAV=',G11.3,' alf=',G11.3)
c
c
c  *********
      HQR=DREAL(HQ)
      HQI=DIMAG(HQ)
c  *********
c--- WE NOW DEFINE THE MATRIX ELEMENTS --------------
c hydrogen density
c
        zamr(1,1) = - GAV+HQR
     &   + ( 1.0 - zflh * ztauh * ( 1.0 + zetah ) ) / zep2nh
        zami(1,1) = HQI
        zamr(1,2) = (HQR-GAV)*ztauh
        zami(1,2) = ztauh*HQI
        zamr(1,3) = (HQR-GAV)*ztauh
        zami(1,3) = ztauh*HQI
        zamr(1,8) = -em*ztauh*HQR*(1.+zetah)/(kpc*zep2nh)
        zami(1,8) = -em*ztauh*HQI*(1.+zetah)/(kpc*zep2nh)
        zamr(1,9) = -em*HQR/kpc
        zami(1,9) = -em*HQI/kpc
c
        zbmr(1,1) = zflh
        zbmr(1,3) = 1.
c
c  hydrogen energy
c
        zamr(2,1) = ( zetah - tvr )/zep2nh
        zamr(2,2) = - ztauh * ftr
c
        zbmr(2,2) = 1.
        zbmr(2,3) = - tvr
c
c  total electron density expressed in ion density and imp density
c
       zamr(3,1) = -1. + zft/zep2ne
       zami(3,1) = vef*(1.-zft)
       zamr(3,3) = 1. - zfnz -zfs
       zami(3,3) = -vef*(1. - zfnz - zfs)
       zamr(3,4) = zft
       zamr(3,5) = zfnz
       zami(3,5) = -vef*zfnz
       zami(3,7) = vef*zft
       zamr(3,8) = -em*(1. - zft)/(kpc*zep2ne)
       zami(3,8) = em*(1.-zft)*vef/(kpc*zep2ne)
       zamr(3,9) = em*(1. - zft)*(1.+eni)/kpc
       zami(3,9) = -em*(1.-zft)*vef/kpc
c
      zbmr(3,1) = zft - 1.
      zbmr(3,3) = 1. - zfnz - zfs
      zbmr(3,5) = zfnz
      zbmr(3,9) = em*(1. - zft)/kpc
c
c  trapped electron energy
c
      zamr(4,1) = zft*(zetae - tvr)/zep2ne
      zami(4,1) = vef*tvr*(bt-2.5*(1.-zft))
      zami(4,3) = -vef*tvr*bt1*(1.-zfnz -zfs)
      zamr(4,4) = zft*ftr
      zami(4,5) = -vef*tvr*bt1*zfnz
      zami(4,7) = -ftr*vef*zft
c
      zbmr(4,1) = (1. - zft)*tvr
      zbmr(4,3) = -(1. - zfnz -zfs)*tvr
      zbmr(4,4) = zft
      zbmr(4,5) = -zfnz*tvr
c
c  impurity density
c
      zamr(5,1) = - GAV +zimp*HQR/A
     & +(1. -zflz*ztauz*(1.+zetaz))/zep2nz
      zami(5,1) = zimp*HQI/A
      zamr(5,5) = (HQR*zimp/A-GAV)*ztauz
      zami(5,5) = zimp*ztauz*HQI/A
      zamr(5,6) = (HQR*zimp/A-GAV)*ztauz
      zami(5,6) = zimp*ztauz*HQI/A
      zamr(5,8) = -em*HQR*zimp*ztauz*(1.+zetaz)/(kpc*zep2nz*A)
      zami(5,8) = -em*HQI*zimp*ztauz*(1.+zetaz)/(kpc*zep2nz*A)
      zamr(5,9) = -em*HQR*zimp/(kpc*A)
      zami(5,9) = -em*HQI*zimp/(kpc*A)
c
      zbmr(5,1) = zflz
      zbmr(5,5) = 1.
c
c  impurity energy
c
      zamr(6,1) = (zetaz - tvr)/zep2nz
      zamr(6,6) = -ztauz*ftr
c
      zbmr(6,5) = -tvr
      zbmr(6,6) = 1.
c
c  variable F
c
      zamr(7,1) = zetae/zep2ne - 1.
      zami(7,1) = vef
      ZAMR(7,7) = 1.
      zami(7,7) = -vef
c
      zbmr(7,1) = -1.
      zbmr(7,7) = 1.
c
c
c  electromagnetic parallel vectorpotential Av = e A_par/Te
c
      fft=(1.-zfnz-zfs)/(1.-zft)
      fzft=zfnz/(1.-zft)
      zamr(8,1) = em1*kpc*(eni+HQR*(fft+zimp*fzft/A))
      zami(8,1) = em1*HQI*(fft+zimp*fzft/A)*kpc
      zamr(8,2) = em1*HQR*ztauh*fft*kpc
      zami(8,2) = em1*HQI*ztauh*fft*kpc
      zamr(8,3) = em1*HQR*ztauh*fft*kpc
      zami(8,3) = em1*HQI*ztauh*fft*kpc
      zamr(8,5) = em1*HQR*zimp*ztauz*fzft*kpc/A
      zami(8,5) = em1*HQI*zimp*ztauz*fzft*kpc/A
      zamr(8,6) = em1*HQR*zimp*ztauz*fzft*kpc/A
      zami(8,6) = em1*HQI*zimp*ztauz*fzft*kpc/A
      zamr(8,8) = em1*((1.+zetae)*eni - alf*zflh*RAV)
     &-em1*HQR*(fft*ztauh*(1.+zetah)/zep2nh
     &+zimp*fzft*ztauz*(1.+zetaz)/(zep2nz*A))
      zami(8,8) = -em1*HQI*(fft*ztauh*(1.+zetah)/zep2nh
     &+zimp*fzft*ztauz*(1.+zetaz)/(zep2nz*A))
      zamr(8,9)= -em1*(eni+HQR*(fft+zimp*fzft/A))
      zami(8,9) = -em1*HQI*(fft+zimp*fzft/A)
c
      zbmr(8,1) = em1*kpc
      zbmr(8,8) = em1
      zbmr(8,9)= -em1
c
c     K = omega*Av
c
      zamr(9,9) = em1
c
      zbmr(9,8) = em1
c
      GO TO 08888
c      else
08887 CONTINUE
c
      write(6,*)
      write(6,*) ieq,' = ieq in sbrtn disp9t'
      call abortb(6
     &,'the value of ieq is wrong in sbrtn disp9t')
c
c       endif
c
08888 CONTINUE
c---------------------------------------------------------------
c  This part is common to all values of ieq
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
      if ( lprint .gt. 9 ) then
        write (6,*)
        write (6,*) ' zamr(j1,j2)  j2 ->'
        do j1=1,ieq
          write (6,192) (zamr(j1,j2),j2=1,ieq)
        enddo
c
        write (6,*)
c
        write (6,*) ' zami(j1,j2)  j2 ->'
        do j1=1,ieq
          write(6,192) (zami(j1,j2),j2=1,ieq)
        enddo
c
       write (6,*)
        write (6,*) ' zbmr(j1,j2)  j2->'
        do j1=1,ieq
          write (6,192) (zbmr(j1,j2),j2=1,ieq)
        enddo
c
       write (6,*)
       write(6,*) ' zbmi(j1,j2)  j2->'
       do j1=1,ieq
         write(6,192) (zbmi(j1,j2),j2=1,ieq) 
       enddo
c
 192  format (1p5e12.4)
      endif
c
c      write(*,193) zep2ne,zep2nh,zep2nz
  193 format(2X,'zep2n=',G12.5,' zep2nh=',G12.5,' zep2nz=',G12.5)
c      zgne=2./zep2ne
c      zgnh=2./zep2nh
c      zgnz=2./zep2nz
c
c      write(*,194) zgne,zgnh,zgnz
  194 format(2X,'zgne=',G12.5,' zgnh=',G12.5,' zgnz=',G12.5)
c

c..find the eigenvalues using NAG14 routine  f02gjf
c
      ztol = max ( 0.D0, cetain(32) )
      ifail = 1
c
  201 continue
c
      call r8tomsqz( idim,ieq,zamr,zami,zbmr,zbmi,
     &zalfr,zalfi,zbeta,zvr,zvi,ifail )
c
  202 continue
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
C        write(6,*) j,omega(j),gamma(j)
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
c--------------------------------------------------------------------
      IF(ieq.le.8) go to 00086 
      IF(ISB.EQ.1) go to 00086
      IF(ISEARCH.EQ.3) go to 10083
c**************************************************************************
c**  HERE THE PROCEEDURE TO FIND THE NEW WZ STARTS  ******************
c
c  Find the fastest growing mode, IM. If ISEARCH=1 only modes propagating
c  in the ion drift direction are considered.
      WM=0.001
      IM=0
      DO 00082 j=1,ieq
      IF(ISEARCH.EQ.2) GO TO 20021
      IF(omega(j).GE.0.) GO TO 00082
20021 IF(gamma(j).LT.WM) GO TO 00082
      WM=gamma(j)
      IM=j
00082 CONTINUE
C
      IMX=MAX(1,IM)
      WS=ZZ(IMX)
      IF(ITERA.EQ.1.AND.ISEARCH.EQ.1) WSF=WZP
c --  DEFINE THE NEW WZ ----
c -- If there is no unstable mode only average WZ is used--
      IMET=0
      IF(IM.GE.1) GO TO 10025
c      WZ=(WZ+WZP)/2.
      IF(ITERA.GT.1) GO TO 20019
      WZI=DIMAG(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI)
      WZA=WZ
      WZ=(WZ+WZP)/2.
      IF(lPRINT.EQ.2) WRITE(*,10028) WZA
10028 FORMAT(' WZA=',2G11.3)
      GO TO 00083
20019 CONTINUE
      WZ=(WZ+WZP)/2.
      GO TO 00083
10025 CONTINUE
      WZ=WS
c----------------------------------
c A new WZ corresponding to an unstable mode has been found --
c ----------------------------------
10026 CONTINUE
c  CHECK FOR CONVERGENCE 
c Independently of method used the average is always taken first
      WZ=(WZ+WZP)/2.
      IF(ITC.EQ.0) GO TO 00083
c      CTEST=ABS((WZ-WZP)/WZP)
      CTEST=ABS((WS-WSF)/WSF)
c      CTEST=ABS((WS-WZ)/WS)
      IF(CTEST.LE.TOL) ITS=1
      IF(ITS.EQ.1) GO TO 00083
c-----
      WZI=DIMAG(WZ)
      IF(ITERA.GT.1) GO TO 10011
      WZI=DIMAG(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI)
      WZA=WZ
c -- In the first iteration we only use average
       GO TO 00083
c------------------------------------
10011 CONTINUE
c-- If ITC=2 we only use average ---------
      IF(ITC.EQ.2) GO TO 00083
c---------------------------------
c-- In the second and higher iteration we expand W=W(WK)
c-- In the second iteration we use Newton-Raphsons method--
      DWN=WZP-WZPP
      IF(ABS(DWN).LT.0.001) GO TO 00083
      DEW=(WZ-WZP)/DWN
c-- In the third and higher iterations we use Mullers method--
      IF(ITERA.GE.3) GO TO 10027  !! We will use Mullers method
c--------------------------------
20003 CONTINUE
      IF(ABS(1.-DEW).LT.0.001) GO TO 00083
      WZ=(WZ-WZP*DEW)/(1.-DEW)  !! This is the Newton Raphson result
      IMET=1
      GO TO 00083
10027 CONTINUE
c---- We will here use Mullers method -----------
      DWN1=WZ-WZP
      IF(ABS(DWN1).LT.0.001) GO TO 00083
c      DEW=(WSF-WSFF)/DWN
      D2EW=(DEW-DEWP)/DWN1
      IF(ABS(D2EW).LT.0.001) GO TO 20003
      AM1=(1.+D2EW*WZP-DEW)/D2EW
      BM1=2.*(WZP*DEW-WSF)/D2EW-WZP*WZP
      DR=SQRT(AM1**2+BM1)
      W1=AM1+DR
      W2=AM1-DR
      DW1=ABS(W1-WZP)
      DW2=ABS(W2-WZP)
      WZ=W1
      IF(DW1.GT.DW2) WZ=W2
c     We have now obtained WZ by Mullers method 
      IMET=2
c------------------------------------------------------------
00083 CONTINUE
c -- THE NEW WZ HAS BEEN OBTAINED
      WZI=DIMAG(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI) !! WZI always larger than 0.01
      WZI=DIMAG(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI) !! May  sometimes  be necessary  
      IF(ABS(WZ).GT.1000.) WZ=WZ*1000./ABS(WZ)  !! NEW 01.03.08
      WZJ(IK)=WZ
c      IF((ABS((WZ-WZP)/WZP)).LE.TOL) ITS=1
      IF(IM.EQ.0) ITS=0
c
C      WRITE(*,20032) IM,ITS,WZJ(IK),IK
20032 FORMAT(' IM=',I5,' ITS=',I5,' WZJ=',2G11.3,' IK=',I5)
c ---  Different conditions for continued iteration are tested ----
      IF(ITC.EQ.0) GO TO 00085
      IF(ITS.EQ.1) GO TO 00085
      IF(ITERA.GE.ITL) GO TO 00084
      IF(IM.EQ.0.AND.ITERA.GE.3) GO TO 00084
c------------------------------------------------------------
c  *** A new iteration will be made ***
20022 CONTINUE
      WZPP=WZP
      DEWP=DEW
      WZP=WZ
      WSFF=WSF
      WSF=WS
C
      DO j1=1,idim
      zalfr(j1)=0.
      zalfi(j1)=0.
      zbeta(j1)=0.
c
      DO j2=1,idim
      zamr(j1,j2)=0.
      zami(j1,j2)=0.
      zbmr(j1,j2)=0.
      zbmi(j1,j2)=0.
      zvr(j1,j2)=0.
      zvi(j1,j2)=0.
       enddo
      enddo
c
      ITERA=ITERA+1
      GO TO 800
c
00084 CONTINUE
      IF(IM.GE.1) GO TO 20024
      IF(ISEARCH.EQ.SEARCHMODE) GO TO 20024
      IF(ISEARCH.EQ.2) GO TO 20024
      ISEARCH=ISEARCH+1  !! Here the search is made with general WZ 
      WZJ(IK)=WZA
      ITERA=0
      GO TO 20022
20024 CONTINUE
C----WHEN THE ITERATIONS DO NOT CONVERGE WE USE THE FIRST AVERAGE--
      WZ=WZA
      WZR=DREAL(WZ)
      WZI=DIMAG(WZ)
      IF(WZR.LT.-10.) WZ=WZ-WZR-10.
      IF(WZI.GT.10.) WZ=WZ+IU*(10.-WZI)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI)
      WZI=DIMAG(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI)
C
      WZJ(IK)=WZ
      IF(lprint.NE.2) GO TO 00085
      WRITE(*,50001) WZ
50001 FORMAT(' NONCONVERGENT WZ=',2G11.4)
c----------------------------------------------------------------------
00085 CONTINUE
c New part using electron WZ for electron root 
      IF(SEARCHMODE.NE.3) GO TO 10085
C
      WRITE(*,90075)
90075 FORMAT(' SEARCHMODE=3')
      IF(ISEARCH.EQ.3) GO TO 10083
      IF(DREAL(WZ).GT.0.) GO TO 10085
      ISEARCH=3
      IM=0
      WM=0.
      DO 10082 j=1,ieq
      WR=DREAL(ZZ(j))
      WI=DIMAG(ZZ(j))
      ZZS(j)=ZZ(j)
      IF(WR.LT.0.) GO TO 10082
      IF(WI.LE.WM) GO TO 10082
      IM=j
      WM=WI
10082 CONTINUE
c
      IF(WM.LT.0.001) GO TO 10085
      WZO=WZ
      WZ=ZZ(IM)
c
      DO j1=1,idim
      zalfr(j1)=0.
      zalfi(j1)=0.
      zbeta(j1)=0.
c
      DO j2=1,idim
      zamr(j1,j2)=0.
      zami(j1,j2)=0.
      zbmr(j1,j2)=0.
      zbmi(j1,j2)=0.
      zvr(j1,j2)=0.
      zvi(j1,j2)=0.
       enddo
      enddo
c
c      WRITE(*,10075) ISEARCH, WZ
10075 FORMAT(' ELECTRON MODE ISEARCH=',I5,' WZ=',2G11.3)
c
      GO TO 30006
10083 CONTINUE  !!! Here ISEARCH=3 and a new results has been obtained
c
      DO 10084 j=1,ieq
      IF(DREAL(ZZ(j)).GT.0.) GO TO 10084
      ZZ(j)=ZZS(j)
10084 CONTINUE
      WZE=WZ
      WZ=WZO
C----------------------------------------------------
10085 CONTINUE
      IF(lprint.NE.2) GO TO 00086
      WRITE(*,50002) EN*WZ,WZJ(IK),IK
50002 FORMAT(' EXIT DISP9T  WZ=',2G11.3,' WZJ(IK)=',2G11.3,'  IK=',I5)
      IF(SEARCHMODE.NE.3) GO TO 00086
      WRITE(*,10087) EN*WZE
10087 FORMAT(' SEARCHMODE=3, Electron WZ=',2G11.3)
00086 CONTINUE
      return
      end
C***************************************************************************
C    DIFFTD calculates transport coefficients
C***************************************************************************
C     File DIFFTD.F
C     Assumes maximum 10 unstable roots
C
      SUBROUTINE DIFF(RP,IR,TAUI,FT,U,CHI,CHE,D,CHQ,DHQ)
C
      IMPLICIT NONE
      SAVE
      INTEGER I,IR,IX,J,IW,LPRINTIN,NDIM,NEQ
      REAL*8 WR,WI,EN,ENI,ENH,EI,EE,TAUI,FT,TE,N,WIN,WDE
      REAL*8 XI(5),XE(5),XD(5),CHI(5),CHE(5),D(5)
      REAL*8 CHQ(5),DHQ(5),XQ(5),XDQ(5)
      REAL*8 U(5,100)
      REAL*8 DNI,DNE,WSQ,NN,DNIN,WRED,fre
      REAL*8 A,B,C,E,F,A1,DH,WR1,WI1,WR2,WI2
      REAL*8 SCHI,SCHE,SD,HP,SHP,NR,NI
      REAL*8 THRD,TVR,FTR,STR,FTRT,RFL,D1,XHH,DT,THRD2
      REAL*8 XIH,XEH,XDH,STF,PHS,KPC
      REAL*8 XDE,XDI,VEF,BTA,GAR,GAI,GBR,GBI,YDA
      REAL*8 HR,DIVGA,DIVGB,DIVA,DIVB,DEVGA,DEVGB,DEVA
      REAL*8 DEVB,SVB,WII,VEFN,BT1,EIH,EEH,EQH,IMP,TZ
      REAL*8 SCHQ,SDQ,KXQ,CHQEFF,DQEFF
      REAL*8 DTIQ,DTQQ,XQH,NQ,TQ,TIQ,NG,LNEH,LNHE
      REAL*8 BE,G,GI,SI,Z,TAUZ,EQ,ENQ,DTIMP,NIMP,H,ZFS
      REAL*8 KIQ,DNQ,DI,DE,DQ,NQR,NQI,K,T,DQT,DQN,TS
      REAL*8 WEXB,ROT,WROR,WIOR,WJR,WJI
      REAL*8 RVA,RVAID,RVAC,RVAC1,RVAC2,RVAC3
      REAL*8 RVAID1,RVAID2,RVAID3,RVB
      REAL*8 SVAID,SVAC,SVAC1,SVAC2,SVAC3
      REAL*8 SVAID1,SVAID2,SVAID3
      REAL*8 NIER,NIEI,RN
      REAL*8 DIVGA1,DIVGA2,DIVGA3,DIVA1,DIVA2,DIVA3
      REAL*8 DEVGA1,DEVGA2,DEVGA3,DEVA1,DEVA2,DEVA3
      REAL*8 DIVGAT,DIVAT,DEVGAT,DEVAT,SVAT,AT,CT,ET,HT,TT,DQNT
      REAL*8 CETAIN(32),ETE,ETI,ETQ,TAUH,AZ,FX,RFLT
      REAL*8 H1,H2,A2,A3,C1,C2,E1,E2,T1,T2,DQN1,DQN2
      REAL*8 HPT(5)
      REAL*8 PI,PE,PN,PQ,PNQ,PIQ,EA
      REAL*8 GCI,GCE,GD,GCQ,GNQ,GK
      REAL*8 LTI,LTE,LN,LTQ,LNQ
      REAL*8 ZVR(10,10),ZVI(10,10)
      COMPLEX*16 FIH,NEF,NRAT,ELN,TEF,FRAT,AV
      COMPLEX*16 NIE
      REAL*8 ELNR,ELNI,AINF,HPE,SHPE,CHIC,IMF,CEFT,SCHEF,GEI,DEFT,DEF
      REAL*8 betae,q,S,Cs,EM,EM1
      COMMON/LT/ LTI,LTE,LN,LTQ,LNQ
      COMPLEX*16 ZZ(10),RP(10),W,WJ,NE,BT2
      COMPLEX*16 GA,GB,GM,IU,HC
      COMMON/GRAD/ EN,ENH,EI,EE
      COMMON/PAR/ THRD,FTR,STR,RFL,D1,XHH,IX,IW
      COMMON/EFF DIFF/ SCHI,SCHE,SD,SCHQ,SDQ
      COMMON/COLL/ BTA,VEF
      COMMON/IMP/ BE,G,SI,Z,TAUZ,ENQ,EQ,KIQ,KXQ,ZFS
      COMMON/WROT/ WEXB,ROT
      COMMON/ETAWN6/ CETAIN,ETE,ETI,ETQ,TAUH,AZ,FX,RFLT,
     1LPRINTIN,NDIM
      COMMON/TP/ EA,HPT
      COMMON/GRKVOT/ GK
      COMMON/ZV/ ZVR,ZVI
      COMMON/ZZ/ ZZ
      COMMON/NEQ/ NEQ
      COMMON/EM/ SHPE,CHIC,SCHEF,DEF
      COMMON/BETAE/ betae,q,S,Cs,EM
      COMMON/WDE/ WDE
C
C--------------------------------------------------------------
      EM1=EM
C      EM1=0.    !!! DEF not used if this line is active !
C-------------------------------------------------------------
C
      TVR=2.D0*THRD
      THRD2=THRD*THRD
      FTRT=FTR*TAUI
      IU=(0.D0,1.D0)
      SHP=0.D0
      CHIC=0.
      SHPE=0.
      SCHEF=0.
      DEF=0.
      ENI=1./EN
      EIH=EI-STR+FTR*ENH
      EEH=EE-STR+FTR*EN
      EQH=EQ-STR+FTR*ENQ
      LNEH=EN/ENH
      LNHE=ENH/EN
      G=1.-Z*BE-ZFS
      NG=MAX(G,1.D-10)
      GI=1./NG
      TZ=Z*TAUZ
      IF(TZ.LT.0.0001) GO TO 00005
      IMP=1./TZ
      KPC=1.
00005 CONTINUE
C
c      WRITE(*,00998) BE,G,SI,Z,TAUZ,ENQ,EQ,KIQ,KXQ,ZFS
c00998 FORMAT(' BQ=',G11.3,' G=',G11.3,' SI=',G11.3,' Z=',G11.3,/,
c     &' TAUZ=',G11.3,' ENQ=',G11.3,' EQ=',G11.3,' KIQ=',G11.3,/,
c     &' KXQ=',G11.3,' ZFS=',G11.3)
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
      DIVGA1=0.D0
      DIVGA2=0.D0
      DIVGA3=0.D0
      DIVGB=0.D0
      DIVA1=0.D0
      DIVA2=0.D0
      DIVA3=0.D0
      DIVB=0.D0
      DEVGA1=0.D0
      DEVGA2=0.D0
      DEVGA3=0.D0
      DEVGB=0.D0
      DEVA1=0.D0
      DEVA2=0.D0
      DEVA3=0.D0
      DEVB=0.D0
      RVA=0.D0
      RVAID=0.D0
      RVAID1=0.D0
      RVAID2=0.D0
      RVAID3=0.D0
      RVAC=0.D0
      RVAC1=0.D0
      RVAC2=0.D0
      RVAC3=0.D0
      SVAID=0.D0
      SVAID1=0.D0
      SVAID2=0.D0
      SVAID3=0.D0
      SVAC=0.D0
      SVAC1=0.D0
      SVAC2=0.D0
      SVAC3=0.D0
      SVB=0.D0
      PHS=0.D0
      RN=0.D0
C
      IF(NEQ.LT.7) VEF=0.D0
c
      DO 34 J=1,5
   34 HPT(J)=0.
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
c      WRITE(*,20011) EA
20011 FORMAT(/,' DIFFTD  EA=',G11.3)
      fre=0.25
c      fre=0.5
      IF(IR.EQ.0) GOTO 00110
c  ---------------------------------------------------
c  Main Loop *******************
c
      DO 01100 J=1,IR
      WJ=RP(J)
      WJR=DREAL(WJ)
      WJI=DIMAG(WJ)
      IF(WJI.LT.0.001) GO TO 01100
c
      IF(IW.NE.1) GO TO 10022
      WROR=WJR*WDE
      WIOR=WJI*WDE
      WRITE(*,10021) WROR,WIOR,J
10021 FORMAT(' Orig W, WR=',G11.3,'/s WI=', G11.3,'/s  J=',I5)
10022 CONTINUE
c
      W=RP(J)-IU*ROT*DABS(WEXB)
      WR=DREAL(W)
      IF(WR.LE.0.) GO TO 10011
      WRED=(DIMAG(RP(J)))**2-fre*ROT*WEXB**2
      IF(WRED.LT.0.D0) WRED=0.D0
      W=WR+IU*DSQRT(WRED)
C      W=RP(J)    !!!!!!!  NOTE TEMPORARY CHANGE  NO STAB OF TE MODE
10011 CONTINUE
      WI=DIMAG(W)
      IF(WI.GT.0.001) GO TO 11 
      WI=0.001
      GO TO 01100
   11 CONTINUE
      W=DCMPLX(WR,WI)
      GM=1.D0+EE*ENI/(W-1.D0+IU*VEF)
      BT1=BTA-2.5D0
      BT2=BTA-2.5D0*GM
      HC=W-FTR+TVR*BT1
      NE=W*W-2.D0*FTR*W+FTR+IU*VEF*HC
      NIE=W*W-2.D0*FTR*W+FTR
      NIER=DREAL(NIE)
      NIEI=DIMAG(NIE)
      WSQ=WR*WR+WI*WI
      WII=1.D0/WI
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
      YDA=WR*(1.D0-EN)+EE-STR+FTR*EN
c Linear trapped electron density response dn?n = FX ephi/Te where **
c **  FX = ENI*(YDA+IU*WI*(1-EN)+IU*VEF*(EN*GA+EE*GB))/NE **
C ** where NE=NER+IU*NEI; NER=NIER-WI*VEF, NEI=NIEI+VEF*HR ****
C   ***************************************
C
      IF(IW.NE.1) GOTO 00021
c      WR1=DABS(EN)*WR
c      WI1=DABS(EN)*WI
      WR2=WDE*WR
      WI2=WDE*WI
      WR1=WR
      WI1=WI
      WRITE(*,00020) WR1,WI1,J
00020 FORMAT(2X,'WR=',D11.5,' WI=',D11.5,' J=',I5)
      WRITE(*,10020) WR2,WI2,J
10020 FORMAT('  WR in sec**-1 ',G11.3,' WI in sec**-1 ',G11.3,' J=',I5)
00021 CONTINUE
C
      NR=DREAL(NE)
      NI=DIMAG(NE)
      NN=(NR)**2+(NI)**2
C
      IF(VEF.EQ.0.D0) GO TO 25
c
c    We write FX = KK*(KR + IU*KI) where KK=1/(NN*EN),
c    KR = RVA + EE*RVB and KI = SVA + EE*SVB
c    We divide into ideal and collisional parts as:
c    RVA = RVAID + RVAC, SVA = SVAID + SVAC 
c
      RVAID=NIER*YDA+NIEI*W*(1.D0-EN)
      SVAID=NIER*WI*(1.D0-EN)-NIEI*YDA
c
      RVAC=VEF*(EN*(NI*GAR-NR*GAI)-WI*YDA+HR*WI*(1.D0-EN))
      RVB=VEF*(NI*GBR-NR*GBI)
c
      SVAC=VEF*(EN*(NR*GAR+NI*GAI)-(1.D0-EN)*WI**2-HR*YDA)
      SVB=VEF*(NR*GBR+NI*GBI)
c
c     These parts are now divided into diagonal, off diagonal and convective 
c     parts as e.g.  RVA = RVA1 + EN*RVA2 + EE*RVA3 etc ...
c
c ** The following parts, due to the ideal part of the density responce
c    enter only for electron thermal transport -----------------------
c
      RVAID1=NIER*(WR-STR)+WI*NIEI
      RVAID2=NIER*(FTR-WR)-WI*NIEI
      RVAID3=NIER
c
      SVAID1=WI*NIER-NIEI*(WR-STR)
      SVAID2=-WI*NIER+NIEI*(WR-FTR)
      SVAID3=-NIEI
c----------------------------------------------------------------------
c
      RVAC1=VEF*WI*(HR-WR+STR)
      RVAC2=VEF*(GAR*NI-GAI*NR+WI*(WR-FTR-HR))
      RVAC3=-VEF*WI
c
      SVAC1=VEF*(HR*(STR-WR)-WI*WI)
      SVAC2=VEF*(GAR*NR+GAI*NI+WI*WI+HR*(WR-FTR))
      SVAC3=-VEF*HR
c
c     The ideal parts RVAID etc .. will only enter into the electron thermal 
c     transport. For the particle transport we need only SVAC and SVB
c
c  Ion thermal conductivity
c
      DIVGA1=FTRT*RVAC1
      DIVGA2=FTRT*RVAC2
      DIVGA3=FTRT*RVAC3
c
      DIVA1=XDI*SVAC1
      DIVA2=XDI*SVAC2
      DIVA3=XDI*SVAC3
c
      DIVGB=FTRT*RVB
      DIVB=XDI*SVB
c
c   Electron thermal conductivity
c
      DEVGA1=-FTR*RVAC1-BT1*VEF*(SVAID1+SVAC1)
      DEVGA2=-FTR*RVAC2-BT1*VEF*(SVAID2+SVAC2)
      DEVGA3=-FTR*RVAC3-BT1*VEF*(SVAID3+SVAC3)
c
      DEVA1=XDE*SVAC1-BT1*VEF*(WR-FTR)*(RVAID1+RVAC1)
      DEVA2=XDE*SVAC2-BT1*VEF*(WR-FTR)*(RVAID2+RVAC2)
      DEVA3=XDE*SVAC3-BT1*VEF*(WR-FTR)*(RVAID3+RVAC3)
c
      DEVGB=-(FTR*RVB+BT1*VEF*SVB)
      DEVB=XDE*SVB-BT1*VEF*(WR-FTR)*RVB
c
      PHS=(WR-FTR)*DREAL(BT2)+WI*DIMAG(BT2) !! independent of density resp.
c
      RN=1.D0/NN

      DIVGA1=DIVGA1*RN
      DIVGA2=DIVGA2*RN
      DIVGA3=DIVGA3*RN
      DIVGB=DIVGB*RN
      DEVGA1=DEVGA1*RN
      DEVGA2=DEVGA2*RN
      DEVGA3=DEVGA3*RN
      DEVGB=DEVGB*RN
      DIVA1=DIVA1*RN
      DIVA2=DIVA2*RN
      DIVA3=DIVA3*RN
      DIVB=DIVB*RN
      DEVA1=DEVA1*RN
      DEVA2=DEVA2*RN
      DEVA3=DEVA3*RN
      DEVB=DEVB*RN
      SVAC1=SVAC1*RN
      SVAC2=SVAC2*RN
      SVAC3=SVAC3*RN
      SVB=SVB*RN
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
     &-FTR*FTR*IMP**3*(2.D0*WR+IMP*FTR))/NIMP
      DQN1=(-NQR+2.D0*(WR+IMP*STR)*(WR+FTR*IMP))/NIMP
      DQN2=(NQR-2.D0*(WR+IMP*FTR)*(WR+FTR*IMP))/NIMP
C
      DIVGAT=DIVGA1+EN*DIVGA2+EE*DIVGA3
      DIVAT=DIVA1+EN*DIVA2+EE*DIVA3
      DEVGAT=DEVGA1+EN*DEVGA2+EE*DEVGA3
      DEVAT=DEVA1+EN*DEVA2+EE*DEVA3
      SVAT=SVAC1+EN*SVAC2+EE*SVAC3
      AT=EN*A2+A3
      CT=C1-1.D0+EN*C2
      ET=E1+EN*E2
      HT=H1+ENQ*H2
      TT=T1+ENQ*T2
      DQNT=DQN1+EN*DQN2
C
C **** IMPURITIES *****
C
      NIMP=(WR*(WR+2.*FTR*IMP)-WI*WI+FTR*IMP*IMP)**2
     1+4.D0*WI*WI*(WR+FTR*IMP)**2
      DTIMP=(WSQ*(WSQ*(ENQ-1.D0)+2.*IMP*WR*EQH+FTR*IMP
     1*IMP*(2.D0*EQ-11.D0/3.D0+STR*ENQ)+FTR*TAUI*IMP*(1.D0
     1+EQ-FTR*ENQ))+FTR*TAUI*(2.D0*FTR*IMP*IMP*WR*(1.D0-ENQ)
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
      H=(WSQ*(WSQ*(ENQ-1.D0)-2.D0*IMP*WR*(STR-FTR*ENQ)
     1+FTR*IMP*IMP*(-11.D0*THRD+STR*ENQ)+FTR*TAUI*IMP*(1.D0-FTR
     1*ENQ))+FTR*FTR*TAUI*IMP*IMP*(2.D0*WR*(1.D0-ENQ)+(STR-FTR
     1*ENQ)*IMP))/NIMP
      K=IMP*(FTR*FTR*TAUI*IMP*IMP-WSQ*(2.D0*WR
     1+FTR*(2.D0*IMP+TAUI)))/NIMP
C
C  *************
      T=(WSQ*(WSQ*(ENQ-1.D0)-2.D0*IMP*WR*(STR-FTR*ENQ)
     1+FTR*IMP*IMP*(-8.D0*THRD+TVR*ENQ))+FTR*FTR*(IMP)**3
     1*(2.D0*WR*(1.D0-ENQ)+IMP*(STR-FTR*ENQ)))/NIMP
      TS=IMP*(FTR*FTR*(IMP)**3-WSQ*(2.D0*WR+5.D0*IMP))/NIMP
C
      T=KXQ*T
      TS=KXQ*TS
      DQN=(ENQ-1.D0)*NQR+2.D0*((1.D0-ENQ)*WR+IMP*(STR-FTR*ENQ)
     1)*(WR+FTR*IMP)
      DQT=2.D0*IMP*(WR+FTR*IMP)
      DQN=DQN/NIMP
      DQT=DQT/NIMP
C
      DTIQ=H-EQ*K
      DTQQ=T-EQ*TS
C
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
C      writE(6,'(I3,5F15.7)') IX,DH,-DEVGB,-DEVGA3,-WII*(DEVB+DEVA3),
C     &     DH-DEVGB-DEVGA3-WII*(DEVB+DEVA3)
C      write(6,'(I3,3F15.7)') IX,WII,DEVB
      XE(4)=0.D0
      XE(5)=0.D0
C
      PE=FT*TVR*XEH*(C2+DEVGA2+WII*(DEVA2+VEF*PHS))*EA
C
C
      XD(1)=0.D0
      XD(2)=-FT*XDH*(N/TE)*(F+WII*(SVB+SVAC3))
      XD(3)=FT*XDH*(E1-WII*SVAC1)
      XD(4)=0.D0
      XD(5)=0.D0
C
      PN=-FT*XDH*(E2-WII*SVAC2)*EA
C
C
      XQ(1)=0.D0
      XQ(2)=0.D0
      XQ(3)=0.D0
      XQ(4)=(1.D0+TVR*TS)*XQH
      XQ(5)=-TVR*TQ*(1.D0+T1)*XQH/NQ
C      write(6,*) XQ(5),TVR*TQ*(1.D0+T1)*XQH,NQ
C
      PQ=TVR*XQH*T2*EA
C
C
      XDQ(1)=0.D0
      XDQ(2)=0.D0
      XDQ(3)=0.D0
      XDQ(4)=-XDH*DQT*NQ/TQ
      XDQ(5)=XDH*DQN1
C
      PNQ=-XDH*DQN2*EA
C
      TIQ=1.D0/(TAUZ*KIQ)
C
      HP=XIH*GI*TAUI*TVR*FTR*(1.D0-FT)*EA
      GCI=XI(1)+(LNHE*(XI(2)*EE+XI(3)*N/TE)+XI(4)*EQ*TIQ
     1+XI(5)*NQ/(TE*KIQ))/(EI*TAUI)-(HP+PI+PIQ)*GK*2.D0*LTI
      GCE=XE(2)+XE(3)*N/(TE*EE)-PE*GK*2.D0*LTE
      GD=XD(2)*TE*EE/N+XD(3)-PN*GK*2.D0*LN
      GCQ=XQ(4)+XQ(5)*NQ/(TQ*EQ)-PQ*GK*2.D0*LTQ
C      write(6,*) XQ(4),XQ(5),XQ(5)*NQ/(TQ*EQ)
      GNQ=XDQ(4)*TQ*EQ/NQ+XDQ(5)-PNQ*GK*2.D0*LNQ
C
C
      SHP=SHP+HP
C
      HPT(1)=HPT(1)+HP+PI+PIQ
      HPT(2)=HPT(2)+PE
      HPT(3)=HPT(3)+PN
      HPT(4)=HPT(4)+PQ
      HPT(5)=HPT(5)+PNQ
C
      SCHI=SCHI+GCI
      SCHE=SCHE+GCE
      SD=SD+GD
      SCHQ=SCHQ+GCQ
      SDQ=SDQ+GNQ
C
C
      CHQEFF=D1*TE**1.5*WI**3*(EQ-TVR-TVR*DTQQ)/DQ
      DQEFF=XDH*(DQN-EQ*DQT)
00222 FORMAT(2X,'CHQEFF=',G11.4,' DQEFF=',G11.4,' XIH=',G11.3,
     &' XEH=',G11.3)
C
      IF(IW.NE.1) GO TO 00093
      WRITE(*,00222) CHQEFF,DQEFF,XIH,XEH
      WRITE(*,00099) SCHI,SCHE,SD,SCHQ,SDQ,J
00099 FORMAT(2X,'SCHI=',F11.5,' SCHE=',F11.5,' SD=',
     1F11.5,' SCHQ=',F11.5,' SDQ=',F11.5,' J=',I5)
      WRITE(*,00223) XDH,DNE,WI
00223 FORMAT(' XDH=',G11.3,' DNE=',G11.3,' WI=',G11.3)
C      WRITE(*,00229) XEH,DH,DEVGB,DEVGA3,WII,DEVB,DEVA3,XE(2)
00229 FORMAT(' XEH=',G11.3,' DH=',G11.3,' DEVGB=',G11.3,' DEVGA3=',G11.3,
     &/,' WII=',G11.3,' DEVB=',G11.3,' DEVA3=',G11.3,' XE(2)=',G11.3)
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
00100 CONTINUE
C
C
01100 CONTINUE
c-----------------------------------------------  MAIN LOOP
C
C--------------------------------------------------------------
      IF(NEQ.LE.8) GO TO 00110
C 
      IF(EM.EQ.0.AND.EM1.EQ.0.) GO TO 00110
C  IF electromagnetic effects or collisions on free electrons are included
C  the transport coefficients are corrected for this.
C------------------------------------------------------------
      DT=D1*TE**1.5D0
      SHPE=0.D0
      SCHEF=0.D0
      DEF=0.D0
      DO 00200 J=1,NEQ
      WR=DREAL(ZZ(J))
      IF(DIMAG(ZZ(J)).LE.0.01D0) GO TO 00200
      WI=DIMAG(ZZ(J))-ROT*DABS(WEXB)
      W=WR+IU*WI
      IF(WR.LE.0.D0) GO TO 10012
      WRED=(DIMAG(ZZ(J)))**2-fre*ROT*WEXB**2
      IF(WRED.LT.0.D0) WRED=0.D0
      W=WR+IU*DSQRT(WRED)
C      W=ZZ(J)  !!!!  NOTE TEMPORARY CHANGE NO STAB OF TE MODE
10012 CONTINUE
      WI=DIMAG(W)
      WIN=DIMAG(ZZ(J))
      IF(WIN.LE.1.D-3) WIN=1.D-3
      IF(WI.LT.1.D-3) GO TO 00200
c --- contr. to chii from em free electr. ----
c
      DNI=(WR+FTR*TAUI)**2+WI*WI
      DNIN=(WR+FTR*TAUI)**2+WIN*WIN
      XDH=DT*WI**3/RFL
      XIH=XDH/DNI
      FIH=DCMPLX(ZVR(1,J),ZVI(1,J))
      IF(NEQ.EQ.11) GO TO 197
C----------------------------------------------
      IF(NEQ.EQ.9) GO TO 00095
c-- Here NEF for disp10 is defined ---
      NEF=DCMPLX(ZVR(4,J),ZVI(4,J))
      GO TO 00097
C----------------------------------------
00095 CONTINUE
C-- Here NEF for disp 9 is defined --
      AV=DCMPLX(ZVR(8,J),ZVI(8,J))
      NEF=FIH-(ZZ(J)-ENI)*AV/KPC
C--------------------------------------------
      GO TO 00097
  197 CONTINUE
      AV=DCMPLX(ZVR(9,J),ZVI(9,J))
      NEF=FIH-(ZZ(J)-ENI)*AV/KPC
      TEF=EE*ENI*AV/KPC
00097 CONTINUE
      IF(CDABS(FIH).LT.0.0001) FIH=(0.0001,0.)
      NRAT=NEF/FIH
      ELN=NRAT-1.D0
      ELNR=DREAL(ELN)
      ELNI=DIMAG(ELN)
      AINF=TVR*(FTR*TAUI*ELNR+ELNI*(WI*WI+WR*(WR+FTR*TAUI))/WIN)
      HPE=XIH*GI*(1.D0-FT)*EA*AINF
      IF(IW.NE.1) GO TO 20035
C      WRITE(*,20021) ELNR,ELNI,WR,WI,WIN,AINF,HPE
20021 FORMAT(' ELNR=',G11.3,' ELNI=',G11.3,' WR=',G11.3,' WI=',G11.3,
     &' WIN=',G11.3,' AINF=',G11.3,' HPE=',G11.3)
20035 CONTINUE
      SHPE=SHPE+HPE
c
c ****  Free electron heat flux *********************
c
C-----------------------------------------------------
      IF(NEQ.EQ.11) GO TO 10099
      IF(NEQ.EQ.9) GO TO 10098
c--- Here TEF for disp10 is defined ---
      TEF=DCMPLX(ZVR(6,J),ZVI(6,J))
      GO TO 10099
10098 CONTINUE
c-- Here TEF for disp9 is defined ---
      TEF=EE*ENI*AV/KPC
C-------------------------------------------------------------
10099 CONTINUE
      FRAT=TEF/FIH
      IMF=-DIMAG(FRAT)*DNIN/DNI
      CEFT=(1.-FT)*IMF*DT*ETE*WI**3/(RFL*WIN)
      SCHEF=SCHEF+CEFT
      IF(IW.NE.1) GO TO 20036
C      WRITE(*,20022) FRAT,IMF,CEFT,SCHEF
20022 FORMAT(' FRAT=',2G11.3,' IMF=',G11.3,' CEFT=',G11.3,
     &' SCHEF=',G11.3)
20036 CONTINUE
C**********************************************************
c ----Free electron particle flux -----------
      GEI=-DIMAG(NRAT)/WI
      DEFT=(1.D0-FT)*GEI*EN*XDH
      DEF=DEF+DEFT
c -------
00200 CONTINUE
c
      CHIC=-SHPE*GK*2.D0*LTI
      HPT(1)=HPT(1)+EM*SHPE
      CHE(2)=CHE(2)+EM*SCHEF
      D(3)=D(3)+EM1*DEF
c
C-------------------------------------------------------
00110 CONTINUE
c
c      WRITE(*,21001) HPT(1),HPT(2),HPT(3),HPT(4),SHPE
21001 FORMAT(' HP1',G11.3,' HP2',G11.3,' HP3',G11.3,' HP4',G11.3,
     &' SHPE=',G11.3)
C
      RETURN
      END
c$$$C******************************************************************
c$$$C  r8tmosqz solver of linear system of equations
c$$$C******************************************************************
c$$$!dmc -- real*8 version generated using `fgtok'.
c$$$!  names changed:  tomsqz -> r8tomsqz, cqzhes -> r8cqzhes, etc.
c$$$!  all constants converted to "D" expontent form
c$$$!  all declarations remapped to REAL*8 / COMPLEX*16
c$$$!  cpp for standardizing REAL/COMPLEX intrinsics:
c$$$!
c$$$!#include "f77_dcomplx.h"
c$$$!
c$$$      SUBROUTINE R8TOMSQZ(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI,IFAIL)
c$$$C-----------------------------------------------------------------------
c$$$C TOMSQZ  written by P. Strand 27-apr-98,       elfps@elmagn.chalmers.se
c$$$C-----------------------------------------------------------------------
c$$$C CQZHES, CQZVEC, and CQZVAL: Fortran subroutines implementing the QZ
c$$$C algorithm for solving the generalized eigenvalue problem for complex
c$$$C matrices. (See B.S.C Garbow, ACM TOMS 4 (1978) pp. 404-410.).
c$$$C-----------------------------------------------------------------------
c$$$C
c$$$C ON INPUT THE GENERALIZED EIGENVALUE PROBLEM IS DEFINED THROUGH THE
c$$$C COMPLEX MATRICES
c$$$C
c$$$C       A = cmplx (AR, AI)  AND   B = cmplx (BR, BI)
c$$$C
c$$$C WHERE LEADING DIMENSION N IS AS DEFINED IN THE CALLING ROUTINE AND
c$$$C WHERE NA IS THE ROW  RANGE IN THE CURRENT PROBLEM. THE EIGENVALUE
c$$$C PROBLEM IS THEN DEFINED THROUGH
c$$$C
c$$$C       A x = w B x
c$$$C
c$$$C WHERE  THE COMPLEX EIGENVECTORS
c$$$C
c$$$C       x = cmplx (ZVR, ZVI)
c$$$C
c$$$C TOGETHER WITH THE COMPLEX EIGENVALUE
c$$$C
c$$$C        w = cmplx(alfr, alfi)/beta
c$$$C
c$$$C IS OUTPUT FROM THE ROUTINE
c$$$C
c$$$C IFAIL WILL BE NONZERO IF CONVERGENCE HAS NOT BEEN REACH WITHIN 50
c$$$C ITERATIONS
c$$$C-----------------------------------------------------------------------
c$$$C DECLARATIONS FOR INPUT VARIABLES
c$$$C-----------------------------------------------------------------------
c$$$ 
c$$$      IMPLICIT NONE
c$$$      INTEGER N, NA
c$$$      REAL*8 AR(N,NA),AI(N,NA),BR(N,NA),BI(N,NA)
c$$$ 
c$$$C-----------------------------------------------------------------------
c$$$C DECALRATIONS FOR OUTPUT VARIABLES
c$$$C-----------------------------------------------------------------------
c$$$ 
c$$$      REAL*8 ALFR(N),ALFI(N),BETA(N)
c$$$      REAL*8 ZVR(N,NA), ZVI(N,NA)
c$$$      INTEGER IFAIL
c$$$ 
c$$$C-----------------------------------------------------------------------
c$$$C LOCAL VARIABLES
c$$$C-----------------------------------------------------------------------
c$$$ 
c$$$      LOGICAL WANTX
c$$$      REAL*8 EPS1
c$$$      REAL*8 ZERO,ZONE
c$$$ 
c$$$C-----------------------------------------------------------------------
c$$$C START OF ACTUAL CODING
c$$$C-----------------------------------------------------------------------
c$$$ 
c$$$      WANTX = .TRUE.
c$$$      EPS1  = -0.0D0
c$$$ 
c$$$      CALL R8CQZHES(N,NA,AR,AI,BR,BI,WANTX,ZVR,ZVI)
c$$$      CALL R8CQZVAL(N,NA,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,WANTX,
c$$$     &            ZVR,ZVI,IFAIL)
c$$$      CALL R8CQZVEC(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI)
c$$$      RETURN
c$$$      END
c$$$ 
c$$$C     ------------------------------------------------------------------
c$$$C
c$$$      SUBROUTINE R8CQZHES(NM,N,AR,AI,BR,BI,MATZ,ZR,ZI)
c$$$C
c$$$      IMPLICIT NONE
c$$$      INTEGER I,J,K,L,N,K1,LB,L1,NM,NK1,NM1
c$$$      REAL*8 AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ZR(NM,N),ZI(NM,N)
c$$$      REAL*8 R,S,T,TI,U1,U2,XI,XR,YI,YR,RHO,U1I
c$$$      LOGICAL MATZ
c$$$      REAL*8 ZERO,ZONE
c$$$C
c$$$C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FIRST STEP OF THE
c$$$C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
c$$$C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
c$$$C
c$$$C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX GENERAL MATRICES AND
c$$$C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM WITH REAL (AND NON-
c$$$C     NEGATIVE) SUBDIAGONAL ELEMENTS AND THE OTHER TO UPPER TRIANGULAR
c$$$C     FORM USING UNITARY TRANSFORMATIONS.  IT IS USUALLY FOLLOWED BY
c$$$C     CQZVAL  AND POSSIBLY  CQZVEC.
c$$$C
c$$$C     ON INPUT-
c$$$C
c$$$C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
c$$$C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
c$$$C          DIMENSION STATEMENT,
c$$$C
c$$$C        N IS THE ORDER OF THE MATRICES,
c$$$C
c$$$C        A=(AR,AI) CONTAINS A COMPLEX GENERAL MATRIX,
c$$$C
c$$$C        B=(BR,BI) CONTAINS A COMPLEX GENERAL MATRIX,
c$$$C
c$$$C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
c$$$C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
c$$$C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
c$$$C
c$$$C     ON OUTPUT-
c$$$C
c$$$C        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
c$$$C          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO, AND THE
c$$$C          SUBDIAGONAL ELEMENTS HAVE BEEN MADE REAL (AND NON-NEGATIVE),
c$$$C
c$$$C        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
c$$$C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
c$$$C
c$$$C        Z=(ZR,ZI) CONTAINS THE PRODUCT OF THE RIGHT HAND
c$$$C          TRANSFORMATIONS IF MATZ HAS BEEN SET TO .TRUE.
c$$$C          OTHERWISE, Z IS NOT REFERENCED.
c$$$C
c$$$C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
c$$$C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
c$$$C
c$$$C     ------------------------------------------------------------------
c$$$C
c$$$C     ********** INITIALIZE Z **********
c$$$ 
c$$$      ZERO = 0.0D0
c$$$      ZONE = 1.0D0
c$$$      IF (.NOT. MATZ) GO TO 10
c$$$C
c$$$      DO 3 I = 1, N
c$$$C
c$$$         DO 2 J = 1, N
c$$$            ZR(I,J) = 0.0D0
c$$$            ZI(I,J) = 0.0D0
c$$$    2    CONTINUE
c$$$C
c$$$         ZR(I,I) = 1.0D0
c$$$    3 CONTINUE
c$$$C     ********** REDUCE B TO UPPER TRIANGULAR FORM WITH
c$$$C                TEMPORARILY REAL DIAGONAL ELEMENTS **********
c$$$   10 IF (N .LE. 1) GO TO 170
c$$$      NM1 = N - 1
c$$$C
c$$$      DO 100 L = 1, NM1
c$$$         L1 = L + 1
c$$$         S = 0.0D0
c$$$C
c$$$         DO 20 I = L, N
c$$$            S = S + ABS(BR(I,L)) + ABS(BI(I,L))
c$$$   20    CONTINUE
c$$$C
c$$$         IF (S .EQ. ZERO) GO TO 100
c$$$         RHO = 0.0D0
c$$$C
c$$$         DO 25 I = L, N
c$$$            BR(I,L) = BR(I,L) / S
c$$$            BI(I,L) = BI(I,L) / S
c$$$            RHO = RHO + BR(I,L)**2 + BI(I,L)**2
c$$$   25    CONTINUE
c$$$C
c$$$         R = SQRT(RHO)
c$$$         XR = abs(CMPLX(BR(L,L),BI(L,L)))
c$$$         IF (XR .EQ. ZERO) GO TO 27
c$$$         RHO = RHO + XR * R
c$$$         U1 = -BR(L,L) / XR
c$$$         U1I = -BI(L,L) / XR
c$$$         YR = R / XR + 1.0D0
c$$$         BR(L,L) = YR * BR(L,L)
c$$$         BI(L,L) = YR * BI(L,L)
c$$$         GO TO 28
c$$$C
c$$$   27    BR(L,L) = R
c$$$         U1 = -1.0D0
c$$$         U1I = 0.0D0
c$$$C
c$$$   28    DO 50 J = L1, N
c$$$            T = 0.0D0
c$$$            TI = 0.0D0
c$$$C
c$$$            DO 30 I = L, N
c$$$               T = T + BR(I,L) * BR(I,J) + BI(I,L) * BI(I,J)
c$$$               TI = TI + BR(I,L) * BI(I,J) - BI(I,L) * BR(I,J)
c$$$   30       CONTINUE
c$$$C
c$$$            T = T / RHO
c$$$            TI = TI / RHO
c$$$C
c$$$            DO 40 I = L, N
c$$$               BR(I,J) = BR(I,J) - T * BR(I,L) + TI * BI(I,L)
c$$$               BI(I,J) = BI(I,J) - T * BI(I,L) - TI * BR(I,L)
c$$$   40       CONTINUE
c$$$C
c$$$            XI = U1 * BI(L,J) - U1I * BR(L,J)
c$$$            BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
c$$$            BI(L,J) = XI
c$$$   50    CONTINUE
c$$$C
c$$$         DO 80 J = 1, N
c$$$            T = 0.0D0
c$$$            TI = 0.0D0
c$$$C
c$$$            DO 60 I = L, N
c$$$               T = T + BR(I,L) * AR(I,J) + BI(I,L) * AI(I,J)
c$$$               TI = TI + BR(I,L) * AI(I,J) - BI(I,L) * AR(I,J)
c$$$   60       CONTINUE
c$$$C
c$$$            T = T / RHO
c$$$            TI = TI / RHO
c$$$C
c$$$            DO 70 I = L, N
c$$$               AR(I,J) = AR(I,J) - T * BR(I,L) + TI * BI(I,L)
c$$$               AI(I,J) = AI(I,J) - T * BI(I,L) - TI * BR(I,L)
c$$$   70       CONTINUE
c$$$C
c$$$            XI = U1 * AI(L,J) - U1I * AR(L,J)
c$$$            AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
c$$$            AI(L,J) = XI
c$$$   80    CONTINUE
c$$$C
c$$$         BR(L,L) = R * S
c$$$         BI(L,L) = 0.0D0
c$$$C
c$$$         DO 90 I = L1, N
c$$$            BR(I,L) = 0.0D0
c$$$            BI(I,L) = 0.0D0
c$$$   90    CONTINUE
c$$$C
c$$$  100 CONTINUE
c$$$C     ********** REDUCE A TO UPPER HESSENBERG FORM WITH REAL SUBDIAGONAL
c$$$C                ELEMENTS, WHILE KEEPING B TRIANGULAR **********
c$$$      DO 160 K = 1, NM1
c$$$         K1 = K + 1
c$$$C     ********** SET BOTTOM ELEMENT IN K-TH COLUMN OF A REAL **********
c$$$         IF (AI(N,K) .EQ. ZERO) GO TO 105
c$$$         R = abs(CMPLX(AR(N,K),AI(N,K)))
c$$$         U1 = AR(N,K) / R
c$$$         U1I = AI(N,K) / R
c$$$         AR(N,K) = R
c$$$         AI(N,K) = 0.0D0
c$$$C
c$$$         DO 103 J = K1, N
c$$$            XI = U1 * AI(N,J) - U1I * AR(N,J)
c$$$            AR(N,J) = U1 * AR(N,J) + U1I * AI(N,J)
c$$$            AI(N,J) = XI
c$$$  103    CONTINUE
c$$$C
c$$$         XI = U1 * BI(N,N) - U1I * BR(N,N)
c$$$         BR(N,N) = U1 * BR(N,N) + U1I * BI(N,N)
c$$$         BI(N,N) = XI
c$$$  105    IF (K .EQ. NM1) GO TO 170
c$$$         NK1 = NM1 - K
c$$$C     ********** FOR L=N-1 STEP -1 UNTIL K+1 DO -- **********
c$$$         DO 150 LB = 1, NK1
c$$$            L = N - LB
c$$$            L1 = L + 1
c$$$C     ********** ZERO A(L+1,K) **********
c$$$            S = ABS(AR(L,K)) + ABS(AI(L,K)) + AR(L1,K)
c$$$            IF (S .EQ. ZERO) GO TO 150
c$$$            U1 = AR(L,K) / S
c$$$            U1I = AI(L,K) / S
c$$$            U2 = AR(L1,K) / S
c$$$            R = SQRT(U1*U1+U1I*U1I+U2*U2)
c$$$            U1 = U1 / R
c$$$            U1I = U1I / R
c$$$            U2 = U2 / R
c$$$            AR(L,K) = R * S
c$$$            AI(L,K) = 0.0D0
c$$$            AR(L1,K) = 0.0D0
c$$$C
c$$$            DO 110 J = K1, N
c$$$               XR = AR(L,J)
c$$$               XI = AI(L,J)
c$$$               YR = AR(L1,J)
c$$$               YI = AI(L1,J)
c$$$               AR(L,J) = U1 * XR + U1I * XI + U2 * YR
c$$$               AI(L,J) = U1 * XI - U1I * XR + U2 * YI
c$$$               AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
c$$$               AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
c$$$  110       CONTINUE
c$$$C
c$$$            XR = BR(L,L)
c$$$            BR(L,L) = U1 * XR
c$$$            BI(L,L) = -U1I * XR
c$$$            BR(L1,L) = -U2 * XR
c$$$C
c$$$            DO 120 J = L1, N
c$$$               XR = BR(L,J)
c$$$               XI = BI(L,J)
c$$$               YR = BR(L1,J)
c$$$               YI = BI(L1,J)
c$$$               BR(L,J) = U1 * XR + U1I * XI + U2 * YR
c$$$               BI(L,J) = U1 * XI - U1I * XR + U2 * YI
c$$$               BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
c$$$               BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
c$$$  120       CONTINUE
c$$$C     ********** ZERO B(L+1,L) **********
c$$$            S = ABS(BR(L1,L1)) + ABS(BI(L1,L1)) + ABS(BR(L1,L))
c$$$            IF (S .EQ. ZERO) GO TO 150
c$$$            U1 = BR(L1,L1) / S
c$$$            U1I = BI(L1,L1) / S
c$$$            U2 = BR(L1,L) / S
c$$$            R = SQRT(U1*U1+U1I*U1I+U2*U2)
c$$$            U1 = U1 / R
c$$$            U1I = U1I / R
c$$$            U2 = U2 / R
c$$$            BR(L1,L1) = R * S
c$$$            BI(L1,L1) = 0.0D0
c$$$            BR(L1,L) = 0.0D0
c$$$C
c$$$            DO 130 I = 1, L
c$$$               XR = BR(I,L1)
c$$$               XI = BI(I,L1)
c$$$               YR = BR(I,L)
c$$$               YI = BI(I,L)
c$$$               BR(I,L1) = U1 * XR + U1I * XI + U2 * YR
c$$$               BI(I,L1) = U1 * XI - U1I * XR + U2 * YI
c$$$               BR(I,L) = U1 * YR - U1I * YI - U2 * XR
c$$$               BI(I,L) = U1 * YI + U1I * YR - U2 * XI
c$$$  130       CONTINUE
c$$$C
c$$$            DO 140 I = 1, N
c$$$               XR = AR(I,L1)
c$$$               XI = AI(I,L1)
c$$$               YR = AR(I,L)
c$$$               YI = AI(I,L)
c$$$               AR(I,L1) = U1 * XR + U1I * XI + U2 * YR
c$$$               AI(I,L1) = U1 * XI - U1I * XR + U2 * YI
c$$$               AR(I,L) = U1 * YR - U1I * YI - U2 * XR
c$$$               AI(I,L) = U1 * YI + U1I * YR - U2 * XI
c$$$  140       CONTINUE
c$$$C
c$$$            IF (.NOT. MATZ) GO TO 150
c$$$C
c$$$            DO 145 I = 1, N
c$$$               XR = ZR(I,L1)
c$$$               XI = ZI(I,L1)
c$$$               YR = ZR(I,L)
c$$$               YI = ZI(I,L)
c$$$               ZR(I,L1) = U1 * XR + U1I * XI + U2 * YR
c$$$               ZI(I,L1) = U1 * XI - U1I * XR + U2 * YI
c$$$               ZR(I,L) = U1 * YR - U1I * YI - U2 * XR
c$$$               ZI(I,L) = U1 * YI + U1I * YR - U2 * XI
c$$$  145       CONTINUE
c$$$C
c$$$  150    CONTINUE
c$$$C
c$$$  160 CONTINUE
c$$$C
c$$$  170 RETURN
c$$$C     ********** LAST CARD OF CQZHES **********
c$$$      END
c$$$C
c$$$C     ------------------------------------------------------------------
c$$$C
c$$$      SUBROUTINE R8CQZVAL(NM,N,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,
c$$$     X                                       MATZ,ZR,ZI,IERR)
c$$$C
c$$$      IMPLICIT NONE
c$$$      INTEGER I,J,K,L,N,EN,K1,K2,LL,L1,NA,NM,ITS,KM1,LM1,
c$$$     X        ENM2,IERR,LOR1,ENORN
c$$$      REAL*8 AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),ALFI(N),
c$$$     X       BETA(N),ZR(NM,N),ZI(NM,N)
c$$$      REAL*8 R,S,A1,A2,EP,SH,U1,U2,XI,XR,YI,YR,ANI,A1I,A33,A34,A43,A44,
c$$$     X       BNI,B11,B33,B44,SHI,U1I,A33I,A34I,A43I,A44I,B33I,B44I,
c$$$     X       EPSA,EPSB,EPS1,ANORM,BNORM,B3344,B3344I
c$$$      INTEGER max
c$$$      LOGICAL MATZ
c$$$      COMPLEX*16 Z3
c$$$C
c$$$      REAL*8 ZERO,ZONE,ZTWO
c$$$C
c$$$C
c$$$C
c$$$C
c$$$C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF STEPS 2 AND 3 OF THE
c$$$C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
c$$$C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,
c$$$C     AS MODIFIED IN TECHNICAL NOTE NASA TN E-7305(1973) BY WARD.
c$$$C
c$$$C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES, ONE OF THEM
c$$$C     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM,
c$$$C     THE HESSENBERG MATRIX MUST FURTHER HAVE REAL SUBDIAGONAL ELEMENTS.
c$$$C     IT REDUCES THE HESSENBERG MATRIX TO TRIANGULAR FORM USING
c$$$C     UNITARY TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM
c$$$C     OF THE OTHER MATRIX AND FURTHER MAKING ITS DIAGONAL ELEMENTS
c$$$C     REAL AND NON-NEGATIVE.  IT THEN RETURNS QUANTITIES WHOSE RATIOS
c$$$C     GIVE THE GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY
c$$$C     CQZHES  AND POSSIBLY FOLLOWED BY  CQZVEC.
c$$$C
c$$$C     ON INPUT-
c$$$C
c$$$C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
c$$$C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
c$$$C          DIMENSION STATEMENT,
c$$$C
c$$$C        N IS THE ORDER OF THE MATRICES,
c$$$C
c$$$C        A=(AR,AI) CONTAINS A COMPLEX UPPER HESSENBERG MATRIX
c$$$C          WITH REAL SUBDIAGONAL ELEMENTS,
c$$$C
c$$$C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
c$$$C
c$$$C        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
c$$$C          EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
c$$$C          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
c$$$C          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
c$$$C          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
c$$$C          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
c$$$C          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
c$$$C          BUT LESS ACCURATE RESULTS,
c$$$C
c$$$C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
c$$$C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
c$$$C          EIGENVECTORS, AND TO .FALSE. OTHERWISE,
c$$$C
c$$$C        Z=(ZR,ZI) CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
c$$$C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
c$$$C          BY  CQZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
c$$$C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
c$$$C
c$$$C     ON OUTPUT-
c$$$C
c$$$C        A HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
c$$$C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
c$$$C
c$$$C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
c$$$C          HAVE BEEN ALTERED.  IN PARTICULAR, ITS DIAGONAL HAS BEEN SET
c$$$C          REAL AND NON-NEGATIVE.  THE LOCATION BR(N,1) IS USED TO
c$$$C          STORE EPS1 TIMES THE NORM OF B FOR LATER USE BY  CQZVEC,
c$$$C
c$$$C        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE
c$$$C          DIAGONAL ELEMENTS OF THE TRIANGULARIZED A MATRIX,
c$$$C
c$$$C        BETA CONTAINS THE REAL NON-NEGATIVE DIAGONAL ELEMENTS OF THE
c$$$C          CORRESPONDING B.  THE GENERALIZED EIGENVALUES ARE THEN
c$$$C          THE RATIOS ((ALFR+I*ALFI)/BETA),
c$$$C
c$$$C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
c$$$C          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE.,
c$$$C
c$$$C        IERR IS SET TO
c$$$C          ZERO       FOR NORMAL RETURN,
c$$$C          J          IF AR(J,J-1) HAS NOT BECOME
c$$$C                     ZERO AFTER 50 ITERATIONS.
c$$$C
c$$$C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
c$$$C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
c$$$C
c$$$C     ------------------------------------------------------------------
c$$$C
c$$$      ZTWO = 2.0D0
c$$$      ZONE = 1.0D0
c$$$      ZERO = 0.0D0
c$$$ 
c$$$      IERR = 0
c$$$C     ********** COMPUTE EPSA,EPSB **********
c$$$      ANORM = 0.0D0
c$$$      BNORM = 0.0D0
c$$$C
c$$$      DO 30 I = 1, N
c$$$         ANI = 0.0D0
c$$$         IF (I .NE. 1) ANI = ABS(AR(I,I-1))
c$$$         BNI = 0.0D0
c$$$C
c$$$         DO 20 J = I, N
c$$$            ANI = ANI + ABS(AR(I,J)) + ABS(AI(I,J))
c$$$            BNI = BNI + ABS(BR(I,J)) + ABS(BI(I,J))
c$$$   20    CONTINUE
c$$$C
c$$$         IF (ANI .GT. ANORM) ANORM = ANI
c$$$         IF (BNI .GT. BNORM) BNORM = BNI
c$$$   30 CONTINUE
c$$$C
c$$$      IF (ANORM .EQ. ZERO) ANORM = 1.0D0
c$$$      IF (BNORM .EQ. ZERO) BNORM = 1.0D0
c$$$      EP = EPS1
c$$$      IF (EP .GT. ZERO) GO TO 50
c$$$C     ********** COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO **********
c$$$      EP = 1.0D0
c$$$   40 EP = EP / ZTWO
c$$$      IF (ZONE + EP .GT. ZONE) GO TO 40
c$$$   50 EPSA = EP * ANORM
c$$$      EPSB = EP * BNORM
c$$$C     ********** REDUCE A TO TRIANGULAR FORM, WHILE
c$$$C                KEEPING B TRIANGULAR **********
c$$$      LOR1 = 1
c$$$      ENORN = N
c$$$      EN = N
c$$$C     ********** BEGIN QZ STEP **********
c$$$   60 IF (EN .EQ. 0) GO TO 1001
c$$$      IF (.NOT. MATZ) ENORN = EN
c$$$      ITS = 0
c$$$      NA = EN - 1
c$$$      ENM2 = NA - 1
c$$$C     ********** CHECK FOR CONVERGENCE OR REDUCIBILITY.
c$$$C                FOR L=EN STEP -1 UNTIL 1 DO -- **********
c$$$   70 DO 80 LL = 1, EN
c$$$         LM1 = EN - LL
c$$$         L = LM1 + 1
c$$$         IF (L .EQ. 1) GO TO 95
c$$$         IF (ABS(AR(L,LM1)) .LE. EPSA) GO TO 90
c$$$   80 CONTINUE
c$$$C
c$$$   90 AR(L,LM1) = 0.0D0
c$$$C     ********** SET DIAGONAL ELEMENT AT TOP OF B REAL **********
c$$$   95 B11 = abs(CMPLX(BR(L,L),BI(L,L)))
c$$$      IF (B11     .EQ. ZERO) GO TO 98
c$$$      U1 = BR(L,L) / B11
c$$$      U1I = BI(L,L) / B11
c$$$C
c$$$      DO 97 J = L, ENORN
c$$$         XI = U1 * AI(L,J) - U1I * AR(L,J)
c$$$         AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
c$$$         AI(L,J) = XI
c$$$         XI = U1 * BI(L,J) - U1I * BR(L,J)
c$$$         BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
c$$$         BI(L,J) = XI
c$$$   97 CONTINUE
c$$$C
c$$$      BI(L,L) = 0.0D0
c$$$   98 IF (L .NE. EN) GO TO 100
c$$$C     ********** 1-BY-1 BLOCK ISOLATED **********
c$$$      ALFR(EN) = AR(EN,EN)
c$$$      ALFI(EN) = AI(EN,EN)
c$$$      BETA(EN) = B11
c$$$      EN = NA
c$$$      GO TO 60
c$$$C     ********** CHECK FOR SMALL TOP OF B **********
c$$$  100 L1 = L + 1
c$$$      IF (B11 .GT. EPSB) GO TO 120
c$$$      BR(L,L) = 0.0D0
c$$$      S = ABS(AR(L,L)) + ABS(AI(L,L)) + ABS(AR(L1,L))
c$$$      U1 = AR(L,L) / S
c$$$      U1I = AI(L,L) / S
c$$$      U2 = AR(L1,L) / S
c$$$      R = SQRT(U1*U1+U1I*U1I+U2*U2)
c$$$      U1 = U1 / R
c$$$      U1I = U1I / R
c$$$      U2 = U2 / R
c$$$      AR(L,L) = R * S
c$$$      AI(L,L) = 0.0D0
c$$$C
c$$$      DO 110 J = L1, ENORN
c$$$         XR = AR(L,J)
c$$$         XI = AI(L,J)
c$$$         YR = AR(L1,J)
c$$$         YI = AI(L1,J)
c$$$         AR(L,J) = U1 * XR + U1I * XI + U2 * YR
c$$$         AI(L,J) = U1 * XI - U1I * XR + U2 * YI
c$$$         AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
c$$$         AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
c$$$         XR = BR(L,J)
c$$$         XI = BI(L,J)
c$$$         YR = BR(L1,J)
c$$$         YI = BI(L1,J)
c$$$         BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
c$$$         BR(L,J) = U1 * XR + U1I * XI + U2 * YR
c$$$         BI(L,J) = U1 * XI - U1I * XR + U2 * YI
c$$$         BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
c$$$  110 CONTINUE
c$$$C
c$$$      LM1 = L
c$$$      L = L1
c$$$      GO TO 90
c$$$C     ********** ITERATION STRATEGY **********
c$$$  120 IF (ITS .EQ. 50) GO TO 1000
c$$$      IF (ITS .EQ. 10) GO TO 135
c$$$C     ********** DETERMINE SHIFT **********
c$$$      B33 = BR(NA,NA)
c$$$      B33I = BI(NA,NA)
c$$$      IF (abs(CMPLX(B33,B33I)) .GE. EPSB) GO TO 122
c$$$      B33 = EPSB
c$$$      B33I = 0.0D0
c$$$  122 B44 = BR(EN,EN)
c$$$      B44I = BI(EN,EN)
c$$$      IF (abs(CMPLX(B44,B44I)) .GE. EPSB) GO TO 124
c$$$      B44 = EPSB
c$$$      B44I = 0.0D0
c$$$  124 B3344 = B33 * B44 - B33I * B44I
c$$$      B3344I = B33 * B44I + B33I * B44
c$$$      A33 = AR(NA,NA) * B44 - AI(NA,NA) * B44I
c$$$      A33I = AR(NA,NA) * B44I + AI(NA,NA) * B44
c$$$      A34 = AR(NA,EN) * B33 - AI(NA,EN) * B33I
c$$$     X    - AR(NA,NA) * BR(NA,EN) + AI(NA,NA) * BI(NA,EN)
c$$$      A34I = AR(NA,EN) * B33I + AI(NA,EN) * B33
c$$$     X     - AR(NA,NA) * BI(NA,EN) - AI(NA,NA) * BR(NA,EN)
c$$$      A43 = AR(EN,NA) * B44
c$$$      A43I = AR(EN,NA) * B44I
c$$$      A44 = AR(EN,EN) * B33 - AI(EN,EN) * B33I - AR(EN,NA) * BR(NA,EN)
c$$$      A44I = AR(EN,EN) * B33I + AI(EN,EN) * B33 - AR(EN,NA) * BI(NA,EN)
c$$$      SH = A44
c$$$      SHI = A44I
c$$$      XR = A34 * A43 - A34I * A43I
c$$$      XI = A34 * A43I + A34I * A43
c$$$      IF (XR .EQ. ZERO .AND. XI .EQ. ZERO) GO TO 140
c$$$      YR = (A33 - SH) / 2.0D0
c$$$      YI = (A33I - SHI) / 2.0D0
c$$$      Z3 = sqrt(CMPLX(YR**2-YI**2+XR,2.0D0*YR*YI+XI))
c$$$      U1 = dble(Z3)
c$$$      U1I = dimag(Z3)
c$$$      IF (YR * U1 + YI * U1I .GE. ZERO) GO TO 125
c$$$      U1 = -U1
c$$$      U1I = -U1I
c$$$  125 Z3 = (CMPLX(SH,SHI) - CMPLX(XR,XI) / CMPLX(YR+U1,YI+U1I))
c$$$     X   / CMPLX(B3344,B3344I)
c$$$      SH = dble(Z3)
c$$$      SHI = dimag(Z3)
c$$$      GO TO 140
c$$$C     ********** AD HOC SHIFT **********
c$$$  135 SH = AR(EN,NA) + AR(NA,ENM2)
c$$$      SHI = 0.0D0
c$$$C     ********** DETERMINE ZEROTH COLUMN OF A **********
c$$$  140 A1 = AR(L,L) / B11 - SH
c$$$      A1I = AI(L,L) / B11 - SHI
c$$$      A2 = AR(L1,L) / B11
c$$$      ITS = ITS + 1
c$$$      IF (.NOT. MATZ) LOR1 = L
c$$$C     ********** MAIN LOOP **********
c$$$      DO 260 K = L, NA
c$$$         K1 = K + 1
c$$$         K2 = K + 2
c$$$         KM1 = max(K-1,L)
c$$$C     ********** ZERO A(K+1,K-1) **********
c$$$         IF (K .EQ. L) GO TO 170
c$$$         A1 = AR(K,KM1)
c$$$         A1I = AI(K,KM1)
c$$$         A2 = AR(K1,KM1)
c$$$  170    S = ABS(A1) + ABS(A1I) + ABS(A2)
c$$$         U1 = A1 / S
c$$$         U1I = A1I / S
c$$$         U2 = A2 / S
c$$$         R = SQRT(U1*U1+U1I*U1I+U2*U2)
c$$$         U1 = U1 / R
c$$$         U1I = U1I / R
c$$$         U2 = U2 / R
c$$$C
c$$$         DO 180 J = KM1, ENORN
c$$$            XR = AR(K,J)
c$$$            XI = AI(K,J)
c$$$            YR = AR(K1,J)
c$$$            YI = AI(K1,J)
c$$$            AR(K,J) = U1 * XR + U1I * XI + U2 * YR
c$$$            AI(K,J) = U1 * XI - U1I * XR + U2 * YI
c$$$            AR(K1,J) = U1 * YR - U1I * YI - U2 * XR
c$$$            AI(K1,J) = U1 * YI + U1I * YR - U2 * XI
c$$$            XR = BR(K,J)
c$$$            XI = BI(K,J)
c$$$            YR = BR(K1,J)
c$$$            YI = BI(K1,J)
c$$$            BR(K,J) = U1 * XR + U1I * XI + U2 * YR
c$$$            BI(K,J) = U1 * XI - U1I * XR + U2 * YI
c$$$            BR(K1,J) = U1 * YR - U1I * YI - U2 * XR
c$$$            BI(K1,J) = U1 * YI + U1I * YR - U2 * XI
c$$$  180    CONTINUE
c$$$C
c$$$         IF (K .EQ. L) GO TO 240
c$$$         AI(K,KM1) = 0.0D0
c$$$         AR(K1,KM1) = 0.0D0
c$$$         AI(K1,KM1) = 0.0D0
c$$$C     ********** ZERO B(K+1,K) **********
c$$$  240    S = ABS(BR(K1,K1)) + ABS(BI(K1,K1)) + ABS(BR(K1,K))
c$$$         U1 = BR(K1,K1) / S
c$$$         U1I = BI(K1,K1) / S
c$$$         U2 = BR(K1,K) / S
c$$$         R = SQRT(U1*U1+U1I*U1I+U2*U2)
c$$$         U1 = U1 / R
c$$$         U1I = U1I / R
c$$$         U2 = U2 / R
c$$$         IF (K .EQ. NA) GO TO 245
c$$$         XR = AR(K2,K1)
c$$$         AR(K2,K1) = U1 * XR
c$$$         AI(K2,K1) = -U1I * XR
c$$$         AR(K2,K) = -U2 * XR
c$$$C
c$$$  245    DO 250 I = LOR1, K1
c$$$            XR = AR(I,K1)
c$$$            XI = AI(I,K1)
c$$$            YR = AR(I,K)
c$$$            YI = AI(I,K)
c$$$            AR(I,K1) = U1 * XR + U1I * XI + U2 * YR
c$$$            AI(I,K1) = U1 * XI - U1I * XR + U2 * YI
c$$$            AR(I,K) = U1 * YR - U1I * YI - U2 * XR
c$$$            AI(I,K) = U1 * YI + U1I * YR - U2 * XI
c$$$            XR = BR(I,K1)
c$$$            XI = BI(I,K1)
c$$$            YR = BR(I,K)
c$$$            YI = BI(I,K)
c$$$            BR(I,K1) = U1 * XR + U1I * XI + U2 * YR
c$$$            BI(I,K1) = U1 * XI - U1I * XR + U2 * YI
c$$$            BR(I,K) = U1 * YR - U1I * YI - U2 * XR
c$$$            BI(I,K) = U1 * YI + U1I * YR - U2 * XI
c$$$  250    CONTINUE
c$$$C
c$$$         BI(K1,K1) = 0.0D0
c$$$         BR(K1,K) = 0.0D0
c$$$         BI(K1,K) = 0.0D0
c$$$         IF (.NOT. MATZ) GO TO 260
c$$$C
c$$$         DO 255 I = 1, N
c$$$            XR = ZR(I,K1)
c$$$            XI = ZI(I,K1)
c$$$            YR = ZR(I,K)
c$$$            YI = ZI(I,K)
c$$$            ZR(I,K1) = U1 * XR + U1I * XI + U2 * YR
c$$$            ZI(I,K1) = U1 * XI - U1I * XR + U2 * YI
c$$$            ZR(I,K) = U1 * YR - U1I * YI - U2 * XR
c$$$            ZI(I,K) = U1 * YI + U1I * YR - U2 * XI
c$$$  255    CONTINUE
c$$$C
c$$$  260 CONTINUE
c$$$C     ********** SET LAST A SUBDIAGONAL REAL AND END QZ STEP **********
c$$$      IF (AI(EN,NA) .EQ. ZERO) GO TO 70
c$$$      R = abs(CMPLX(AR(EN,NA),AI(EN,NA)))
c$$$      U1 = AR(EN,NA) / R
c$$$      U1I = AI(EN,NA) / R
c$$$      AR(EN,NA) = R
c$$$      AI(EN,NA) = 0.0D0
c$$$C
c$$$      DO 270 J = EN, ENORN
c$$$         XI = U1 * AI(EN,J) - U1I * AR(EN,J)
c$$$         AR(EN,J) = U1 * AR(EN,J) + U1I * AI(EN,J)
c$$$         AI(EN,J) = XI
c$$$         XI = U1 * BI(EN,J) - U1I * BR(EN,J)
c$$$         BR(EN,J) = U1 * BR(EN,J) + U1I * BI(EN,J)
c$$$         BI(EN,J) = XI
c$$$  270 CONTINUE
c$$$C
c$$$      GO TO 70
c$$$C     ********** SET ERROR -- BOTTOM SUBDIAGONAL ELEMENT HAS NOT
c$$$C                BECOME NEGLIGIBLE AFTER 50 ITERATIONS **********
c$$$ 1000 IERR = EN
c$$$C     ********** SAVE EPSB FOR USE BY CQZVEC **********
c$$$ 1001 IF (N .GT. 1) BR(N,1) = EPSB
c$$$      RETURN
c$$$C     ********** LAST CARD OF CQZVAL **********
c$$$      END
c$$$C
c$$$C     ------------------------------------------------------------------
c$$$C
c$$$      SUBROUTINE R8CQZVEC(NM,N,AR,AI,BR,BI,ALFR,ALFI,BETA,ZR,ZI)
c$$$C
c$$$      IMPLICIT NONE
c$$$      INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN
c$$$      REAL*8 AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),ALFI(N),
c$$$     X       BETA(N),ZR(NM,N),ZI(NM,N)
c$$$      REAL*8 R,T,RI,TI,XI,ALMI,ALMR,BETM,EPSB
c$$$      COMPLEX*16 Z3
c$$$      REAL*8 ZERO
c$$$C
c$$$C
c$$$C
c$$$C
c$$$C
c$$$C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FOURTH STEP OF THE
c$$$C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
c$$$C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
c$$$C
c$$$C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES IN UPPER
c$$$C     TRIANGULAR FORM, WHERE ONE OF THEM FURTHER MUST HAVE REAL DIAGONAL
c$$$C     ELEMENTS.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM
c$$$C     AND TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
c$$$C     IT IS USUALLY PRECEDED BY  CQZHES  AND  CQZVAL.
c$$$C
c$$$C     ON INPUT-
c$$$C
c$$$C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
c$$$C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
c$$$C          DIMENSION STATEMENT,
c$$$C
c$$$C        N IS THE ORDER OF THE MATRICES,
c$$$C
c$$$C        A=(AR,AI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
c$$$C
c$$$C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX WITH REAL
c$$$C          DIAGONAL ELEMENTS.  IN ADDITION, LOCATION BR(N,1) CONTAINS
c$$$C          THE TOLERANCE QUANTITY (EPSB) COMPUTED AND SAVED IN  CQZVAL,
c$$$C
c$$$C        ALFR, ALFI, AND BETA ARE VECTORS WITH COMPONENTS WHOSE
c$$$C          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
c$$$C          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  CQZVAL,
c$$$C
c$$$C        Z=(ZR,ZI) CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
c$$$C          REDUCTIONS BY  CQZHES  AND  CQZVAL, IF PERFORMED.
c$$$C          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
c$$$C          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
c$$$C
c$$$C     ON OUTPUT-
c$$$C
c$$$C        A IS UNALTERED,
c$$$C
c$$$C        B HAS BEEN DESTROYED,
c$$$C
c$$$C        ALFR, ALFI, AND BETA ARE UNALTERED,
c$$$C
c$$$C        Z CONTAINS THE EIGENVECTORS.  EACH EIGENVECTOR IS NORMALIZED
c$$$C          SO THAT THE MODULUS OF ITS LARGEST COMPONENT IS 1.0 .
c$$$C
c$$$C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
c$$$C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
c$$$C
c$$$C     ------------------------------------------------------------------
c$$$C
c$$$      ZERO = 0.0D0
c$$$      IF (N .LE. 1) GO TO 1001
c$$$      EPSB = BR(N,1)
c$$$C     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********
c$$$      DO 800 NN = 2, N
c$$$         EN = N + 2 - NN
c$$$         NA = EN - 1
c$$$         ALMR = ALFR(EN)
c$$$         ALMI = ALFI(EN)
c$$$         BETM = BETA(EN)
c$$$C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
c$$$         DO 700 II = 1, NA
c$$$            I = EN - II
c$$$            R = 0.0D0
c$$$            RI = 0.0D0
c$$$            M = I + 1
c$$$C
c$$$            DO 610 J = M, EN
c$$$               T = BETM * AR(I,J) - ALMR * BR(I,J) + ALMI * BI(I,J)
c$$$               TI = BETM * AI(I,J) - ALMR * BI(I,J) - ALMI * BR(I,J)
c$$$               IF (J .EQ. EN) GO TO 605
c$$$               XI = T * BI(J,EN) + TI * BR(J,EN)
c$$$               T = T * BR(J,EN) - TI * BI(J,EN)
c$$$               TI = XI
c$$$  605          R = R + T
c$$$               RI = RI + TI
c$$$  610       CONTINUE
c$$$C
c$$$            T = ALMR * BETA(I) - BETM * ALFR(I)
c$$$            TI = ALMI * BETA(I) - BETM * ALFI(I)
c$$$            IF (T .EQ. ZERO .AND. TI .EQ. ZERO) T = EPSB
c$$$            Z3 = CMPLX(R,RI) / CMPLX(T,TI)
c$$$            BR(I,EN) = dble(Z3)
c$$$            BI(I,EN) = dimag(Z3)
c$$$  700    CONTINUE
c$$$C
c$$$  800 CONTINUE
c$$$C     ********** END BACK SUBSTITUTION.
c$$$C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
c$$$C                FOR J=N STEP -1 UNTIL 2 DO -- **********
c$$$      DO 880 JJ = 2, N
c$$$         J = N + 2 - JJ
c$$$         M = J - 1
c$$$C
c$$$         DO 880 I = 1, N
c$$$C
c$$$            DO 860 K = 1, M
c$$$               ZR(I,J) = ZR(I,J) + ZR(I,K) * BR(K,J) - ZI(I,K) * BI(K,J)
c$$$               ZI(I,J) = ZI(I,J) + ZR(I,K) * BI(K,J) + ZI(I,K) * BR(K,J)
c$$$  860       CONTINUE
c$$$C
c$$$  880 CONTINUE
c$$$C     ********** NORMALIZE SO THAT MODULUS OF LARGEST
c$$$C                COMPONENT OF EACH VECTOR IS 1 **********
c$$$      DO 950 J = 1, N
c$$$         T = 0.0D0
c$$$C
c$$$         DO 930 I = 1, N
c$$$            R = abs(CMPLX(ZR(I,J),ZI(I,J)))
c$$$            IF (R .GT. T) T = R
c$$$  930    CONTINUE
c$$$C
c$$$         DO 940 I = 1, N
c$$$            ZR(I,J) = ZR(I,J) / T
c$$$            ZI(I,J) = ZI(I,J) / T
c$$$  940    CONTINUE
c$$$C
c$$$  950 CONTINUE
c$$$C
c$$$ 1001 RETURN
c$$$C     ********** LAST CARD OF CQZVEC **********
c$$$      END
      SUBROUTINE ABORTB(NID,KSTR)
      CHARACTER KSTR*(*)
      WRITE(NID,'(A)') KSTR
      STOP
      END
