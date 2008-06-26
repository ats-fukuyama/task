c     
c==================================================c
c==================================================c
      program ftrd20
c
c==================================================c
c
c                         coded by S. Murakami
c
c                       last modified 2001/08/20
c==================================================c
c   05/03/23    n, T limit is added.
c==================================================c
c
      implicit real*8 (a-h,o-y)
c
      parameter (maxdv=10000)
      parameter	(MAXP=30)
      character *256	DBFILE
      data DBFILE/'/home/murakami/fort/nbidb/dataindex.in'/
c     data DBFILE/'/home/nakajima/NEC2/NBI/data/dataindex.in'/
      real*8 TE, TI, N, ALPHA
      real*8 r(MAXP), P1(MAXP), P2(MAXP), rra(maxp), cvol(maxp)
      dimension qe(maxp),qi(maxp),cb(maxp),cpperp(maxp),cppara(maxp)
     &        , cjb(maxp),cmmnte(maxp),cmmnti(maxp),pts(maxp)
      dimension anexm(10), atixm(10), atexm(10)
      integer NP, STATUS
      integer MODE
c
c
      namelist /nmbgpl/ te,ti,n,alpha,zeff,z1,cn0,cnew,tew,tiw,rmaj,rmin 
      namelist /nmbeam/ enbi,azb,amb,pnbi
      namelist /nmtime/ time
      namelist /prof/  anexm, atixm, atexm
c
      print *, 'start'
c
      open(20,file='prof.out20',status='unknown')
      open(10,file='prof.out10',status='unknown')
      open(11,file='prof.out11',status='unknown')
c
      np=maxp-1
c       
c     pi
      pi  = atan2(0.d0, -1.d0)
      pi2 = 2.d0*pi
c     charge
      echg=1.6021892d-19
c
c     elementaly charge in cgs unit
c      echg=4.803242d-10
c
c     speed of light
      clight=2.99792458d8
c
c     electron mass
      cme=9.109534d-31
c     proton mass
      cmp=1.6726485d-27
c
      TE    = 1000.0
      TI    = 1000.0
      N     = 10.0
      ALPHA = 2.0
      z1    = 1.d0
      zeff  = 1.d0
c
      cn0    = 0.d0
      cnew    = 1.d18
      tew    = 100.d0
      tiw    = 100.d0
c
      azb   = 1.d0 
      amb   = 1.d0
c
c...enbi [eV]
      enbi   = 180.d3
c
      read(5,nmbgpl)
      write(6,nmbgpl)
c
      read(5,nmbeam)
      write(6,nmbeam)
c
      read(5,nmtime)
      write(6,nmtime)
c
      read(5,prof)
      write(6,prof)
c
      pnorm = pnbi/1.d6
c
      print *, " "
      print *, "----------------------- "
      print *, " Parameters"
      print *, " "
      print *, '   pi          =',pi
      print *, '   e           =',echg
      print *, '   c [m/sec]  =',clight
      print *, '   me [Kg]      =',cme
      print *, '   mi [Kg]      =',cmp
c
c
      enbikv = enbi/1.d3
c
      vf = sqrt(2.d0*echg*enbi/(cmp*amb))
c
      print *, ' Zf           =',azb
      print *, ' Mf/Mp        =',amb
      print *, ' E_nbi [eV]   =',enbi
      print *, ' v_f [m/sec]  =',vf
      print *, ' Z_eff        =',zeff
      print *, ' Z_1          =',z1
      print *, ' pnorm        =',pnorm
c     
c      cn0    = 1.0d17
c
      itr = 0
      ti0 = ti
c     
10000 continue
c
      itr = itr+1
      qet = 0.d0
      qit = 0.d0
      ppt = 0.d0
      plt = 0.d0
      cjt = 0.d0
      ptst = 0.d0
c
      vmcm = 1.d2
      delr = 1.d0/dfloat(maxp-1)
c
c     NP    = MAXP
c     MODE  = 0
c     call CALC_PROFILE(TE, TI, N, ALPHA, R, P1, P2, NP, MODE, STATUS)
c      
      if (itr.eq.1) then
         write(6,6020) TE, TI, N, ALPHA
         do i=1, NP
         read(20,11000) iir, r(i), p1(i), p2(i), rra(i), cvol(i)      
         if (p1(i).le.0.d0) then
             p1(i)=1.d-20
         end if
         if (p2(i).le.0.d0) then
             p2(i)=1.d-20
         end if
         write(6,11000) i, R(i), P1(i), P2(i), rra(i), cvol(i)
         end do
      end if
c
11000 format(1h ,i4,5e13.4)
6000	format(I3,')',F7.4,5x,2F10.4)
6010	format('[GET_PROFILE] Te=',f6.0,' Ti=',f6.0,' N=',f6.0,
     &	       ' Alpha=',f6.0)
6020	format('[CALC_PROFILE] Te=',f6.0,' Ti=',f6.0,' N=',f6.0,
     &	       ' Alpha=',f6.0,' MODE=',i2)
c
c
      do 100 ir=1,np
c
      zeta0  = sqrt(p1(ir)/(p1(ir)+p2(ir)))
c
      cne = 0.d0
      teev = 0.d0
      tiev = 0.d0
c
       do 535 jj=1,10
        cne  = cne  +anexm(jj)*r(ir)**(dfloat(jj-1))
        teev = teev +atexm(jj)*r(ir)**(dfloat(jj-1))
        tiev = tiev +atixm(jj)*r(ir)**(dfloat(jj-1))
 535  continue
c
c...n,T low limit
      if (cne.lt.1.d11) then
         cne=1.d11
         write(6,*) " **** cne < 10.d11 ****"
      end if
c
      if (teev.lt.10.d0) then
         teev=10.d0
         write(6,*) " **** teev < 10.d0 ****"
      end if
c
      if (tiev.lt.10.d0) then
         tiev=10.d0
         write(6,*) " **** tiev < 10.d0 ****"
      end if
c
c...transform to CGS to MKS
      cne    = cne*1.d6
      cnn    = cn0*exp(-(1.d0-r(ir))/0.1d0)
c     
      tempi0= echg*tiev
      vti2 = 2.d0*tempi0/cmp
      vti  = sqrt(vti2)
c 
c
      ecev  = teev*(9.d0*pi/16.d0*cmp/cme)**(1.d0/3.d0)
     &        *(amb)**(1.d0/3.d0)*(amb*z1)**(2.d0/3.d0)
c
      vc = sqrt(2.d0*ecev*echg/(amb*cmp))
c
      taus = 0.12*(teev/1.d3)**(1.5d0)/
     &         ((cne/1.d19)*azb**2)*(amb)
c
      ttf = taus/3.d0*log((1.e0+(vf/vc)**3)/(1.e0+(vti/vc)**3)) 
c
      if (time.gt.ttf) then 
         vlst = vti
      else
         tex  = exp(-3*time/taus) 
         vlst = vf*( tex -(vc/vf)**3*(1.e0-tex) )**(1.e0/3.e0)
      end if
c
      print *, ' tau_f        =',ttf
      print *, ' time         =',time
      print *, ' vlst         =',vlst
      print *, ' exp(-3t/ts)  =',tex
      print *, ' v_Ti         =',vti
      print *, ' zeta0        =',zeta0
c
c--------------------------------------------c
c--------------------------------------------c
      a1cx = 3.2345d0
      a2cx = 235.88d0
      a3cx = 0.038371d0
      a4cx = 3.8068d-6
      a5cx = 1.1832d-10
      a6cx = 2.3713d0
c
      sgmcx= 1.d-16*a1cx*dlog(a2cx/enbikv +a6cx)
     &      /(1.d0+a3cx*enbikv +a4cx*enbikv**3.5d0 
     &                      +a5cx*enbikv**5.4d0)
c
      if (cnn.gt.0.d0) then
         taucx = 1.d0/(cnn*vf*sgmcx*1.d-4)
      else
         taucx = 1.d20
      end if
c
c     
      deltv = (vf-vlst)/dfloat(maxdv-1)
      cbden = 0.d0
      ge    = 0.d0
      gi    = 0.d0
      pperp = 0.d0
      ppara = 0.d0
c     
      do 1000 iv=1,maxdv
c
         vv = vlst +deltv*dfloat(iv-1)
         pcx = ((vf**3+vc**3)/(vv**3+vc**3))**(-taus/(3.d0*taucx))
         bb = ((vf**3 +vc**3)/(vv**3+vc**3)
     &          *(vv**3/vf**3))**(1.d0/(3.d0*amb)*(zeff/z1))
c
         cbden = cbden +deltv*(vv**2/(vv**3+vc**3)*pcx)
         ge    = ge    +deltv*(vv**4/(vv**3+vc**3)*pcx)
         gi    = gi    +deltv*(vv*vc**3/(vv**3+vc**3)*pcx)
         cke   = cke   +deltv*(vv**3/(vv**3+vc**3)*pcx*bb)
         cki   = cki   +deltv*(vc**3/(vv**3+vc**3)*pcx*bb)
     &                   *(1.d0 +zeff/(amb*z1))
         pperp = pperp +deltv*(vv**4/(vv**3+vc**3)*pcx
     &                   *(1.d0-0.5d0*(3.d0*zeta0**2-1.d0)*bb**3))
         ppara = ppara +deltv*(vv**4/(vv**3+vc**3)*pcx
     &                   *(1.d0+(3.d0*zeta0**2-1.d0)*bb**3))
 1000 continue
c
      ge = 2.d0/vf**2*ge
      gi = 2.d0/vf**2*gi
      pperp = 2.d0/(3.d0*vf**2)*pperp
      ppara = 2.d0/(3.d0*vf**2)*ppara
      cke = cke/vf   
      cki = cki/vf   

c
      qe(ir) = pnorm*(p1(ir)+p2(ir))*ge
      qi(ir) = pnorm*(p1(ir)+p2(ir))*gi
c
      cb(ir) = pnorm*(p1(ir)+p2(ir))/(echg*enbi)*taus*cbden
c
      cpperp(ir) = pnorm*(p1(ir)+p2(ir))*taus*pperp
      cppara(ir) = pnorm*(p1(ir)+p2(ir))*taus*ppara
c
      cjb(ir) = pnorm*(p1(ir)+p2(ir))/(echg*enbi)
     &         *azb*echg*vf*zeta0*taus*cke*1.d6 
      cmmnte(ir) = pnorm*(p1(ir)+p2(ir))/(echg*enbi)
     &            *amb*cmp*vf*zeta0*cke
      cmmnti(ir) = pnorm*(p1(ir)+p2(ir))/(echg*enbi)
     &            *amb*cmp*vf*zeta0*cki
      pts(ir) = pnorm*(p1(ir)+p2(ir))/(echg*enbi)
c
c--------------------------------------------c
c--------------------------------------------c
c
      if (itr.eq.1) then
         write(6,*) 
         write(6,*) '--------------------------'
         write(6,*) '  ir =',ir, '  r/a=',r(ir)
         write(6,*) 
         
c--------------------------------------------c
c     
         print *, ' n0  [m^-3]   =',cnn
         print *, ' Te0 [eV]     =',teev
         print *, ' Ti0 [eV]     =',tiev
         print *, ' Vthi [m/sec] =',vti
         print *, ' E_c [eV]     =',ecev
         print *, ' v_c  [m/sec] =',vc
         print *, ' tau_s [sec]  =',taus
         print *, ' sgmcx [cm^2] =',sgmcx
         print *, ' tau_cx [sec] =',taucx
c     
         print *, ' cbden            =',cbden
         print *, ' ge               =',ge
         print *, ' gi               =',gi
         print *, ' cke               =',cke
         print *, ' cki               =',cki
         print *, ' pperp            =',pperp
         print *, ' ppara            =',ppara
         print *, ' q_e  [W/cm^-3]   =',qe(ir)
         print *, ' q_i  [W/cm^-3]   =',qi(ir)
         print *, ' n_b  [m^-3]      =',cb(ir)
         print *, ' jb   [A/cm^-3]   =',cjb(ir)
         print *, ' p_perp           =',cpperp(ir)
         print *, ' p_para           =',cppara(ir)
c     
c
      qet  = qet+qe(ir)*cvol(ir)
      qit  = qit+qi(ir)*cvol(ir)
      qtot = qet+qit
      ppt = ppt +cpperp(ir)*cvol(ir)
      plt = plt +cppara(ir)*cvol(ir)
c      cjt = cjt +cjb(ir)*cvol(ir)
       cjt = cjt +cjb(ir)*(cvol(ir)*1.d-6)/(2.d0*pi*rmaj)
       ptst = ptst +pts(ir)*cvol(ir)
c
      qtei= qe(ir)+qi(ir) 
c
       write(10,7000) ir, r(ir),qtei,qe(ir),qi(ir),qtot,qet,qit
     &                 ,cb(ir),cpperp(ir),cppara(ir),ppt,plt
     &                 ,cjb(ir),cjt,pts(ir),ptst
      end if
c
c     
 100  continue
c
      qtot = qet+qit
      qtot = qet+qit
         print *, ' cbden            =',cbden
         print *, ' ge               =',ge
         print *, ' gi               =',gi
c     
 7000 format(i4,17e15.5)
c
c================================c
c     p_b iteration      
c================================c
c
c
      pbm = (2.d0*cpperp(1) +cppara(1))/3.d0
      pbgi = n*1.d13*echg*(ti0)
c
      prate = pbm/pbgi
c
      tinew = ti0*(1+prate)
      write(6,*) '  itr =',itr,'  prate ,tinew=',prate,tinew
c
      deltti = (tinew-ti)/tinew
c
      if (deltti.gt.1.d-2) then
         ti = tinew
         go to 10000
      else
         do 200 ir=1,np
         write(11,7000) ir, r(ir),qe(ir),qi(ir)
     &                    ,cb(ir),cpperp(ir),cppara(ir)
 200     continue
      end if
c     
      stop
      end
