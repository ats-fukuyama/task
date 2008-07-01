!
!==================================================c
!==================================================c
!      PROGRAM ftrd20
       SUBROUTINE fit(infile)
!==================================================c

!                         coded by S. Murakami

!                       last modified 2001/08/20
!==================================================c
!   05/03/23    n, T limit is added.
!   07/12       to Fortran90, subroutine
!==================================================c

  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: infile
  INTEGER(4),PARAMETER:: maxdv=10000, maxp=30
  INTEGER(4)::i, iir, ir, itr, iv, jj, np
  REAL(4):: z1, zeff, zeta0
  REAL(8):: a1cx, a2cx, a3cx, a4cx, a5cx, a6cx, alpha, amb, azb, bb, cbden, cjt, cke,&
    cki, clight, cme, cmp, cn0, cne, cnew, cnn, deltti, deltv, ecev, echg, enbi,     &
    enbikv, ge, gi, n, pbgi, pbm, pcx, pi, plt, pnbi, pnorm, ppara, pperp, ppt, prate,&
    ptst, qet, qit, qtei, qtot, rmaj, rmin, sgmcx, taucx, taus, te, teev, tempi0,    &
    tew, tex, ti, ti0, tiev, time, tinew, tiw, ttf, vc, vf, vlst, vti, vti2, vv,     &
    anexm(10), atexm(10), atixm(10)
  REAL(8), ALLOCATABLE::                                                             &
    cb(:), cjb(:), cmmnte(:), cmmnti(:), cppara(:), cpperp(:), cvol(:), p1(:), p2(:),&
    pts(:), qe(:), qi(:), r(:), rra(:)         ! (maxp)

      NAMELIST /nmbgpl/ te,ti,n,alpha,zeff,z1,cn0,cnew,tew,tiw,rmaj,rmin
      NAMELIST /nmbeam/ enbi,azb,amb,pnbi
      NAMELIST /nmtime/ time
      NAMELIST /prof/  anexm, atixm, atexm


      ALLOCATE(cb(maxp), cjb(maxp), cmmnte(maxp), cmmnti(maxp), cppara(maxp), cpperp(maxp))
      ALLOCATE(cvol(maxp),p1(maxp),p2(maxp),pts(maxp),qe(maxp), qi(maxp),r(maxp),rra(maxp))

      WRITE(6,*) 'start'

!!      open(20,file='prof.out20',status='unknown')
!!      open(10,file='prof.out10',status='unknown')
!!      open(11,file='prof.out11',status='unknown')

      np=maxp-1
!
!     pi
      pi  = atan2(0.d0, -1.d0)
!!      pi2 = 2.d0*pi    ! 200710
!     charge
      echg=1.6021892d-19

!     elementaly charge in cgs unit
!      echg=4.803242d-10

!     speed of light
      clight=2.99792458d8

!     electron mass
      cme=9.109534d-31
!     proton mass
      cmp=1.6726485d-27

      TE    = 1000.0
      TI    = 1000.0
      N     = 10.0
      ALPHA = 2.0
      z1    = 1.d0
      zeff  = 1.d0

      cn0    = 0.d0
      cnew    = 1.d18
      tew    = 100.d0
      tiw    = 100.d0

      azb   = 1.d0
      amb   = 1.d0

!...enbi [eV]
      enbi   = 180.d3

      READ(5,nmbgpl)
      WRITE(6,nmbgpl)

      READ(5,nmbeam)
      WRITE(6,nmbeam)

      READ(5,nmtime)
      WRITE(6,nmtime)

      READ(5,prof)
      WRITE(6,prof)

      pnorm = pnbi/1.d6

      WRITE(6,*) " "
      WRITE(6,*) "----------------------- "
      WRITE(6,*) " Parameters"
      WRITE(6,*) " "
      WRITE(6,*) '   pi          =',pi
      WRITE(6,*) '   e           =',echg
      WRITE(6,*) '   c [m/sec]  =',clight
      WRITE(6,*) '   me [Kg]      =',cme
      WRITE(6,*) '   mi [Kg]      =',cmp


      enbikv = enbi/1.d3

      vf = sqrt(2.d0*echg*enbi/(cmp*amb))

      WRITE(6,*) ' Zf           =',azb
      WRITE(6,*) ' Mf/Mp        =',amb
      WRITE(6,*) ' E_nbi [eV]   =',enbi
      WRITE(6,*) ' v_f [m/sec]  =',vf
      WRITE(6,*) ' Z_eff        =',zeff
      WRITE(6,*) ' Z_1          =',z1
      WRITE(6,*) ' pnorm        =',pnorm
!
!      cn0    = 1.0d17

      itr = 0
      ti0 = ti
!
10000 CONTINUE

      itr = itr+1
      qet = 0.d0
      qit = 0.d0
      ppt = 0.d0
      plt = 0.d0
      cjt = 0.d0
      ptst = 0.d0

!!      vmcm = 1.d2                    ! 200710
!!      delr = 1.d0/dfloat(maxp-1)     ! 200710

!     NP    = MAXP
!     MODE  = 0
!     call CALC_PROFILE(TE, TI, N, ALPHA, R, P1, P2, NP, MODE, STATUS)
!
      IF (itr.eq.1) THEN
        WRITE(6,6020) TE, TI, N, ALPHA
        DO i=1, NP
!!          READ(20,11000) iir, r(i), p1(i), p2(i), rra(i), cvol(i)
          READ(infile,11000) iir, r(i), p1(i), p2(i), rra(i), cvol(i)
          IF (p1(i).le.0.d0) THEN
             p1(i)=1.d-20
          END IF
          IF (p2(i).le.0.d0) THEN
             p2(i)=1.d-20
          END IF
          WRITE(6,11000) i, R(i), P1(i), P2(i), rra(i), cvol(i)
        END DO
      END IF

11000 FORMAT(1h ,i4,5e13.4)
6000  FORMAT(I3,')',F7.4,5x,2F10.4)
6010  FORMAT('[GET_PROFILE] Te=',f6.0,' Ti=',f6.0,' N=',f6.0, ' Alpha=',f6.0)
6020  FORMAT('[CALC_PROFILE] Te=',f6.0,' Ti=',f6.0,' N=',f6.0,' Alpha=',f6.0,' MODE=',i2)


      DO ir=1,np

      zeta0  = sqrt(p1(ir)/(p1(ir)+p2(ir)))

      cne = 0.d0
      teev = 0.d0
      tiev = 0.d0

       DO jj=1,10
        cne  = cne  +anexm(jj)*r(ir)**(dfloat(jj-1))
        teev = teev +atexm(jj)*r(ir)**(dfloat(jj-1))
        tiev = tiev +atixm(jj)*r(ir)**(dfloat(jj-1))
       END DO

!...n,T low limit
      IF (cne.lt.1.d11) THEN
        cne=1.d11
        WRITE(6,*) " **** cne < 10.d11 ****"
      END IF

      IF (teev.lt.10.d0) THEN
        teev=10.d0
        WRITE(6,*) " **** teev < 10.d0 ****"
      END IF

      IF (tiev.lt.10.d0) THEN
        tiev=10.d0
        WRITE(6,*) " **** tiev < 10.d0 ****"
      END IF

!...transform to CGS to MKS
      cne    = cne*1.d6
      cnn    = cn0*exp(-(1.d0-r(ir))/0.1d0)
!
      tempi0= echg*tiev
      vti2 = 2.d0*tempi0/cmp
      vti  = sqrt(vti2)
!

      ecev  = teev*(9.d0*pi/16.d0*cmp/cme)**(1.d0/3.d0) &
              *(amb)**(1.d0/3.d0)*(amb*z1)**(2.d0/3.d0)

      vc = sqrt(2.d0*ecev*echg/(amb*cmp))

      taus = 0.12*(teev/1.d3)**(1.5d0)/((cne/1.d19)*azb**2)*(amb)

      ttf = taus/3.d0*log((1.e0+(vf/vc)**3)/(1.e0+(vti/vc)**3))

      IF (time.gt.ttf) THEN
        tex=0.d0
        vlst = vti
      ELSE
        tex  = exp(-3*time/taus)
        vlst = vf*( tex -(vc/vf)**3*(1.e0-tex) )**(1.e0/3.e0)
      END IF

      WRITE(6,*), ' tau_f        =',ttf
      WRITE(6,*) ' time         =',time
      WRITE(6,*) ' vlst         =',vlst
      WRITE(6,*) ' exp(-3t/ts)  =',tex
      WRITE(6,*) ' v_Ti         =',vti
      WRITE(6,*) ' zeta0        =',zeta0

!--------------------------------------------c
!--------------------------------------------c
      a1cx = 3.2345d0
      a2cx = 235.88d0
      a3cx = 0.038371d0
      a4cx = 3.8068d-6
      a5cx = 1.1832d-10
      a6cx = 2.3713d0

      sgmcx= 1.d-16*a1cx*dlog(a2cx/enbikv +a6cx)           &
            /(1.d0+a3cx*enbikv +a4cx*enbikv**3.5d0+a5cx*enbikv**5.4d0)

      IF (cnn.gt.0.d0) THEN
        taucx = 1.d0/(cnn*vf*sgmcx*1.d-4)
      ELSE
        taucx = 1.d20
      END IF

!
      deltv = (vf-vlst)/dfloat(maxdv-1)
      cbden = 0.d0
      ge    = 0.d0
      gi    = 0.d0
      cke   = 0.d0
      cki   = 0.d0
      pperp = 0.d0
      ppara = 0.d0
!
      DO iv=1,maxdv

        vv = vlst +deltv*dfloat(iv-1)
        pcx = ((vf**3+vc**3)/(vv**3+vc**3))**(-taus/(3.d0*taucx))
        bb = ((vf**3 +vc**3)/(vv**3+vc**3)*(vv**3/vf**3))**(1.d0/(3.d0*amb)*(zeff/z1))

        cbden = cbden +deltv*(vv**2/(vv**3+vc**3)*pcx)
        ge    = ge    +deltv*(vv**4/(vv**3+vc**3)*pcx)
        gi    = gi    +deltv*(vv*vc**3/(vv**3+vc**3)*pcx)
        cke   = cke   +deltv*(vv**3/(vv**3+vc**3)*pcx*bb)
        cki   = cki   +deltv*(vc**3/(vv**3+vc**3)*pcx*bb)*(1.d0 +zeff/(amb*z1))
        pperp = pperp +deltv*(vv**4/(vv**3+vc**3)*pcx                  &
                                      *(1.d0-0.5d0*(3.d0*zeta0**2-1.d0)*bb**3))
        ppara = ppara +deltv*(vv**4/(vv**3+vc**3)*pcx                  &
                                             *(1.d0+(3.d0*zeta0**2-1.d0)*bb**3))
      END DO

      ge = 2.d0/vf**2*ge
      gi = 2.d0/vf**2*gi
      pperp = 2.d0/(3.d0*vf**2)*pperp
      ppara = 2.d0/(3.d0*vf**2)*ppara
      cke = cke/vf
      cki = cki/vf


      qe(ir) = pnorm*(p1(ir)+p2(ir))*ge
      qi(ir) = pnorm*(p1(ir)+p2(ir))*gi

      cb(ir) = pnorm*(p1(ir)+p2(ir))/(echg*enbi)*taus*cbden

      cpperp(ir) = pnorm*(p1(ir)+p2(ir))*taus*pperp
      cppara(ir) = pnorm*(p1(ir)+p2(ir))*taus*ppara

      cjb(ir) =    pnorm*(p1(ir)+p2(ir))/(echg*enbi)*azb*echg*vf*zeta0*taus*cke*1.d6
      cmmnte(ir) = pnorm*(p1(ir)+p2(ir))/(echg*enbi)*amb*cmp*vf*zeta0*cke
      cmmnti(ir) = pnorm*(p1(ir)+p2(ir))/(echg*enbi)*amb*cmp*vf*zeta0*cki
      pts(ir) =    pnorm*(p1(ir)+p2(ir))/(echg*enbi)

!--------------------------------------------c
!--------------------------------------------c

      IF (itr.eq.1) THEN
        WRITE(6,*)
        WRITE(6,*) '--------------------------'
        WRITE(6,*) '  ir =',ir, '  r/a=',r(ir)
        WRITE(6,*)

!--------------------------------------------c
!
        WRITE(6,*) ' n0  [m^-3]   =',cnn
        WRITE(6,*) ' Te0 [eV]     =',teev
        WRITE(6,*) ' Ti0 [eV]     =',tiev
        WRITE(6,*) ' Vthi [m/sec] =',vti
        WRITE(6,*) ' E_c [eV]     =',ecev
        WRITE(6,*) ' v_c  [m/sec] =',vc
        WRITE(6,*) ' tau_s [sec]  =',taus
        WRITE(6,*) ' sgmcx [cm^2] =',sgmcx
        WRITE(6,*) ' tau_cx [sec] =',taucx
!
        WRITE(6,*) ' cbden            =',cbden
        WRITE(6,*) ' ge               =',ge
        WRITE(6,*) ' gi               =',gi
        WRITE(6,*) ' cke               =',cke
        WRITE(6,*) ' cki               =',cki
        WRITE(6,*) ' pperp            =',pperp
        WRITE(6,*) ' ppara            =',ppara
        WRITE(6,*) ' q_e  [W/cm^-3]   =',qe(ir)
        WRITE(6,*) ' q_i  [W/cm^-3]   =',qi(ir)
        WRITE(6,*) ' n_b  [m^-3]      =',cb(ir)
        WRITE(6,*) ' jb   [A/cm^-3]   =',cjb(ir)
        WRITE(6,*) ' p_perp           =',cpperp(ir)
        WRITE(6,*) ' p_para           =',cppara(ir)
!

        qet  = qet+qe(ir)*cvol(ir)
        qit  = qit+qi(ir)*cvol(ir)
        qtot = qet+qit
        ppt = ppt +cpperp(ir)*cvol(ir)
        plt = plt +cppara(ir)*cvol(ir)
!       cjt = cjt +cjb(ir)*cvol(ir)
        cjt = cjt +cjb(ir)*(cvol(ir)*1.d-6)/(2.d0*pi*rmaj)
        ptst = ptst +pts(ir)*cvol(ir)

        qtei= qe(ir)+qi(ir)

!!        WRITE(10,7000) ir, r(ir),qtei,qe(ir),qi(ir),qtot,qet,qit      &
        WRITE(infile+40,7000) ir, r(ir),qtei,qe(ir),qi(ir),qtot,qet,qit      &
                     ,cb(ir),cpperp(ir),cppara(ir),ppt,plt,cjb(ir),cjt,pts(ir),ptst
      END IF

!
      END DO          ! do ir=1,np

      qtot = qet+qit
      qtot = qet+qit
         WRITE(6,*) ' cbden            =',cbden
         WRITE(6,*) ' ge               =',ge
         WRITE(6,*) ' gi               =',gi
!
 7000 FORMAT(i4,17e15.5)

!================================c
!     p_b iteration
!================================c


      pbm = (2.d0*cpperp(1) +cppara(1))/3.d0
      pbgi = n*1.d13*echg*(ti0)

      prate = pbm/pbgi

      tinew = ti0*(1+prate)
      WRITE(6,*) '  itr =',itr,'  prate ,tinew=',prate,tinew

      deltti = (tinew-ti)/tinew

      IF (deltti.gt.1.d-2) THEN
         ti = tinew
         GO TO 10000
      ELSE
         DO ir=1,np
!!           WRITE(11,7000) ir, r(ir),qe(ir),qi(ir),cb(ir),cpperp(ir),cppara(ir)
           WRITE(infile+41,7000) ir, r(ir),qe(ir),qi(ir),cb(ir),cpperp(ir),cppara(ir)
         END DO
      END IF
!
      DEALLOCATE(cb,cjb,cmmnte,cmmnti,cppara,cpperp,cvol,p1,p2,pts,qe,qi,r,rra)

!      STOP
      RETURN
      END SUBROUTINE
