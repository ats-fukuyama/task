c@testglf.f
c 12-mar-03 version 1.61
c stand-alone driver for the GLF23 model
c written by Jon Kinsey, General Atomics
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
c
c...declare variables
c
      character cdate*24, line*132
c     character ctime*8  ! used if calling routine clock w/ +U77 option
      integer jpd
 
      double precision epsilon
      parameter ( jpd=50 )
c
      double precision te_m(0:jpd), ti_m(0:jpd)
     & , ne_m(0:jpd), ni_m(0:jpd), ns_m(0:jpd)
     & , zpte_m(0:jpd), zpti_m(0:jpd), zpne_m(0:jpd), zpni_m(0:jpd)
     & , angrotp_exp(0:jpd), egamma_exp(0:jpd), gamma_p_exp(0:jpd)
     & , vphi_m(0:jpd), vpar_m(0:jpd), vper_m(0:jpd)
     & , zeff_exp(0:jpd), bt_exp, bteff_exp(0:jpd), rho(0:jpd), arho_exp
     & , gradrho_exp(0:jpd), gradrhosq_exp(0:jpd)
     & , rmin_exp(0:jpd), rmaj_exp(0:jpd), rmajor_exp
     & , q_exp(0:jpd), shat_exp(0:jpd), alpha_exp(0:jpd)
     & , elong_exp(0:jpd), zimp_exp, amassimp_exp, amassgas_exp
     & , alpha_e, x_alpha
c
      double precision zpte_in, zpti_in, zpne_in, zpni_in, drho
c
      double precision diffnem, chietem, chiitim
     & , etaphim, etaparm, etaperm, exchm
     & , diff_m(0:jpd), chie_m(0:jpd), chii_m(0:jpd), etaphi_m(0:jpd)
     & , etapar_m(0:jpd), etaper_m(0:jpd), exch_m(0:jpd)
     & , egamma_m(0:jpd), egamma_d(0:jpd,10), gamma_p_m(0:jpd)
     & , anrate_m(0:jpd), anrate2_m(0:jpd)
     & , anfreq_m(0:jpd), anfreq2_m(0:jpd)
c
      integer lprint, nroot, jshoot, jmm, jmaxm, itport_pt(1:5)
     & , igrad, idengrad, i_delay, j, k, leigen, irotstab, bt_flag, iglf
c
      namelist /nlglf/ leigen, lprint, nroot, iglf, jshoot, jmm, jmaxm
     & , itport_pt, irotstab, te_m, ti_m, ne_m, ni_m, ns_m
     & , igrad, idengrad, zpte_in, zpti_in, zpne_in, zpni_in
     & , angrotp_exp, egamma_exp, gamma_p_exp, vphi_m, vpar_m, vper_m
     & , zeff_exp, bt_exp, bt_flag, rho, arho_exp
     & , gradrho_exp, gradrhosq_exp
     & , rmin_exp, rmaj_exp, rmajor_exp, q_exp, shat_exp
     & , alpha_exp, elong_exp, zimp_exp, amassimp_exp, amassgas_exp
     & , alpha_e, x_alpha
c
 
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      open (4,file='temp')
      open (5,file='in')
      open (6,file='out')
c
      call stripx (5,4,6)
c
      cdate = ' '
c     ctime = ' '
c
      call c9date (cdate)
cray      call clock (ctime)
cibm      call clock_ (ctime)
c
      write (6,*)
      write (6,*) ' GLF23 stand-alone code by Kinsey, GA  ',cdate
c
c..default inputs
c
      epsilon  = 1.e-10
      leigen   = 0   ! for cgg eigenvalue solver
      nroot    = 8   ! number of roots in eigenvalue solver
      iglf     = 0   ! original GLF23 normalization
      jshoot   = 0   ! for time-dependent code
      jmaxm    = 2
      igrad    = 0   ! compute gradients
      idengrad = 2   ! simple dilution
      i_delay  = 0
      itport_pt(1) = 0
      itport_pt(2) = 0
      itport_pt(3) = 0
      itport_pt(4) = 0
      itport_pt(5) = 0
      irotstab     = 1    ! use internally computed ExB shear, 0 for prescribed
      bt_exp       = 1.0
      bt_flag      = 0    ! do not use effective B-field
      bteff_exp    = 1.0  ! effective B-field (used when bt_flag > 0)
      rmajor_exp   = 1.0
      amassgas_exp = 1.0
      zimp_exp     = 6.0
      amassimp_exp = 12.0
      arho_exp     = 1.0
      alpha_e      = 0.   ! ExB shear stabilization
      x_alpha      = 0.   ! alpha stabilization
      zpte_in      = 0.
      zpti_in      = 0.
      zpne_in      = 0.
      zpni_in      = 0.
c
      do j=0,jpd
        te_m(j)   = 0.0
        ti_m(j)   = 0.0
        ne_m(j)   = 0.0
        ni_m(j)   = 0.0
        ns_m(j)   = 0.0
c
        zpte_m(j) = 0.0
        zpti_m(j) = 0.0
        zpne_m(j) = 0.0
        zpni_m(j) = 0.0
c
        angrotp_exp(j)   = 0.0
        egamma_exp(j)    = 0.0
        gamma_p_exp(j)   = 0.0
        vphi_m(j)        = 0.0
        vpar_m(j)        = 0.0
        vper_m(j)        = 0.0
c
        zeff_exp(j)   = 1.0
        rho(j)        = 0.0
        gradrho_exp(j)   = 1.0
        gradrhosq_exp(j) = 1.0
        rmin_exp(j)   = 0.0
        rmaj_exp(j)   = 0.0
        q_exp(j)      = 1.0
        shat_exp(j)   = 0.0
        alpha_exp(j)  = 0.0
        elong_exp(j)  = 1.0
      enddo
c
c..default outputs
c
      diffnem      = 0
      chietem      = 0
      chiitim      = 0
      etaphim      = 0
      etaparm      = 0
      etaperm      = 0
      exchm        = 0
c
      do j=0,jpd
        diff_m(j)   = 0.0
        chie_m(j)   = 0.0
        chii_m(j)   = 0.0
        etaphi_m(j) = 0.0
        etapar_m(j) = 0.0
        etaper_m(j) = 0.0
        exch_m(j)   = 0.0
        egamma_m(j) = 0.0
        gamma_p_m(j)= 0.0
        anrate_m(j) = 0.0
        anrate2_m(j) = 0.0
        anfreq_m(j) = 0.0
        anfreq2_m(j) = 0.0
        do k=1,10
          egamma_d(j,k) = 0.0
        enddo
      enddo
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c..read input file
c
  10  continue
c
      read  (4,20,end=900,err=900) line
  20  format (a)
c
      if ( index ( line, '$nlglf' ) .gt. 0
     &    .or. index ( line, '&nlglf' ) .gt. 0 ) then
c
c  read namelist input
c
        backspace 4
        read  (4,nlglf)
      else
        go to 10
      endif
c
      if ( leigen .gt. 0 ) then
        write (6,*) ' tomsqz eigenvalue solver '
      else
        write (6,*) ' cgg eigenvalue solver '
      endif
c
      if ( lprint .gt. 100 ) write (6,nlglf)
c
      write (6,*)
c
        call callglf2d( leigen, nroot, iglf
     & , jshoot, jmm, jmaxm, itport_pt
     & , irotstab, te_m, ti_m, ne_m, ni_m, ns_m
     & , igrad, idengrad, zpte_in, zpti_in, zpne_in, zpni_in
     & , angrotp_exp, egamma_exp, gamma_p_exp, vphi_m, vpar_m, vper_m
     & , zeff_exp, bt_exp, bt_flag, rho
     & , arho_exp, gradrho_exp, gradrhosq_exp
     & , rmin_exp, rmaj_exp, rmajor_exp, zimp_exp, amassimp_exp
     & , q_exp, shat_exp, alpha_exp, elong_exp, amassgas_exp
     & , alpha_e, x_alpha, i_delay
     & , diffnem, chietem, chiitim, etaphim, etaparm, etaperm
     & , exchm, diff_m, chie_m, chii_m, etaphi_m, etapar_m, etaper_m
     & , exch_m, egamma_m, egamma_d, gamma_p_m
     & , anrate_m, anrate2_m, anfreq_m, anfreq2_m )
c
      do j=1,jmaxm
        drho=rho(j-1)-rho(j)+epsilon
        zpte_m(j)=-(log(te_m(j-1))-log(te_m(j)))/drho
        zpti_m(j)=-(log(ti_m(j-1))-log(ti_m(j)))/drho
        zpne_m(j)=-(log(ne_m(j-1))-log(ne_m(j)))/drho
        zpni_m(j)=-(log(ni_m(j-1))-log(ni_m(j)))/drho
      enddo
c
c..printout
c
 
      if ( lprint .gt. 0 ) then
        write(6,100)
        do j=0,jmaxm
          write (6,110) rho(j), te_m(j), ti_m(j), ne_m(j), ni_m(j)
     & ,        zeff_exp(j), q_exp(j), shat_exp(j)
        enddo
c
        write(6,120)
        do j=0,jmaxm
          write(6,110) rho(j), zpte_m(j), zpti_m(j)
     & ,       zpne_m(j), zpni_m(j)
        enddo
c
        write(6,130)
        do j=0,jmaxm
          write (6,110) rho(j), diff_m(j), chie_m(j), chii_m(j)
     & ,        etaphi_m(j), etapar_m(j), etaper_m(j), exch_m(j)
        enddo
c
        write(6,140)
        do j=0,jmaxm
          write(6,110) rho(j), egamma_m(j), gamma_p_m(j)
     & ,       anrate_m(j),anrate2_m(j), anfreq_m(j), anfreq2_m(j)
        enddo
      endif
c
 100    format(t5,'rho',t13,'Te',t21,'Ti',t29,'ne',t37,'ni'
     & ,      t44,'Zeff',t53,'q',t60,'shear',t68,'#prof')
 110    format (11(0pf8.4))
 120    format(/,t5,'rho',t12,'zpte',t20,'zpti',t28,'zpne'
     & ,      t36,'zpni',t68,'#log-grad')
 130    format(/,t5,'rho',t12,'diff',t20,'chie',t28,'chii'
     & ,      t35,'etaphi',t43,'etapar',t51,'etaper',t60
     & ,      'exch',t68,'#chi')
 140    format(/,t5,'rho',t10,'egamma_m',t19, 'gamma_p',t27
     & ,      'anrate', t35,'anrate2',t43,'anfreq',t51,'anfreq2'
     & ,      t68,'#gamma')
c
 900  continue
c
      stop
      end
