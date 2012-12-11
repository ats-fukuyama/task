c@callglf2d.f
c 12-mar-10 clf, use "kind"
c 12-feb-03 Kinsey, v1.61 retuned GLF23 model
c 11-apr-01 Kinsey, updated for v1.50
c 29-aug-00 Kinsey, added ave_ve for 3pt smoothing of ve (0 for none)
c 05-nov-99 Kinsey, precall routine for glf2d.f
c added i_dengrad switch for dilution
************************************************************************
       subroutine callglf2d(
     >                 !INPUTS
     > leigen,         ! eigenvalue solver
     >                 ! 0 for cgg (default), 1 for tomsqz, 2 for zgeev
     > nroot,          ! no. roots,8 for default, 12 for impurity dynamics
     > iglf,           ! 0 for original GLF23, 1 for retuned version
     > jshoot,         ! jshoot=0 time-dep code;jshoot=1 shooting code
     > jmm,            ! grid number;jmm=0 does full grid jm=1 to jmaxm-1
     > jmaxm,          ! profile grids 0 to jmaxm
     > itport_pt,      ! 1:5 transport flags
     > irotstab,       ! 0 to use egamma_exp; 1 use egamma_m
     > te_m,           ! 0:jmaxm te Kev           itport_pt(2)=1 transport
     > ti_m,           ! 0:jmaxm ti Kev           itport_pt(3)=1 transport
     > ne_m,           ! 0:jmaxm ne 10**19 1/m**3
     > ni_m,           ! 0:jmaxm ni 10**19 1/m**3 itport_pt(1)=1 transport
     > ns_m,           ! 0:jmaxm ns 10**19 1/m**3 
     > idengrad,       ! default 2, for simple dilution
     > zpte_in,        ! externally provided log gradient te w.r.t rho (i_grad=1)
     > zpti_in,        ! externally provided log gradient ti w.r.t rho
     > zpne_in,        ! externally provided log gradient ne w.r.t rho
     > zpni_in,        ! externally provided log gradient ni w.r.t rho
     > angrotp_exp,    ! 0:jmaxm exp plasma toroidal angular velocity 1/sec
     >                 ! if itport_pt(4)=0 itport_pt(5)=0
     > egamma_exp,     ! 0:jmaxm exp exb shear rate in units of csda_exp
     >                 ! if itport_pt(4)=-1 itport_pt(5)=0
     > gamma_p_exp,    ! 0:jmaxm exp par. vel. shear rate in units of csda_exp
     >                 ! if itport_pt(4)=-1 itport_pt(5)=0
     > vphi_m,         ! 0:jmaxm toroidal velocity m/sec
     >                 ! if itport_pt(4)=1 itport_pt(5)=0 otherwise output
     > vpar_m,         ! 0:jmaxm parallel velocity m/sec
     >                 ! if itport_pt(4)=1 itport_pt(5)=1 otherwise output
     > vper_m,         ! 0:jmaxm perp. velocity m/sec
     >                 ! if itport_pt(4)=1 itport_pt(5)=1 otherwise output
     > zeff_exp,       ! 0:jmaxm ne in 10**19 1/m**3
     > bt_exp,         ! vaccuum axis toroidal field in tesla
     > bt_flag,        ! switch for effective toroidal field use in rhosda
     > rho,            ! 0:jmaxm 0 < rho < 1 normalized toroidal flux (rho=rho/rho(a))
     > arho_exp,       ! rho(a), toroidal flux at last closed flux surface (LCFS)
     >                 !   toroidal flux= B0*rho_phys**2/2 (m)
     >                 !   B0=bt_exp, arho_exp=rho_phys_LCFS
     > gradrho_exp,    ! 0:jmaxm dimensionless <|grad rho_phys |**2>
     > gradrhosq_exp,  ! 0:jmaxm dimensionless <|grad rho_phys |>
     >                 !NOTE:can set arho_exp=1.,if gradrho_exp=<|grad rho |>
     >                 !                 and gradrhosq_exp = <|grad rho |**2>
     > rmin_exp,       ! 0:jmaxm minor radius in meters
     > rmaj_exp,       ! 0:jmaxm major radius in meters
     > rmajor_exp,     ! axis major radius
     > zimp_exp,       ! effective Z of impurity
     > amassimp_exp,   ! effective A of impurity
     > q_exp,          ! 0:jmaxm safety factor
     > shat_exp,       ! 0:jmaxm magnetic shear, d (ln q_exp)/ d (ln rho)
     > alpha_exp,      ! 0:jmaxm MHD alpha from experiment
     > elong_exp,      ! 0:jmaxm elongation
     > amassgas_exp,   !  atomic number working hydrogen gas
     > alpha_e,        ! 1 full (0 no) no ExB shear stab
     > x_alpha,        ! 1 full (0 no) alpha stabilization  with alpha_exp
     >                 !-1 full (0 no) self consistent alpha_m stab.
     > i_delay,        !i_delay time delay for ExB shear should be non-zero only
     >                 ! once per step and is less than or equal 10
     >                 !OUTPUTS
     > diffnem,        ! ion plasma diffusivity in m**2/sec
     > chietem,        ! electron ENERGY diffuivity in m**2/sec
     > chiitim,        ! ion      ENERGY diffuivity in m**2/sec
     > etaphim,        ! toroidal velocity diffusivity in m**2/sec
     > etaparm,        ! parallel velocity diffusivity in m**2/sec
     > etaperm,        ! perpendicular velocity diffusivity in m**2/sec
     > exchm,          ! turbulent electron to ion ENERGY exchange in MW/m**3
     >                 ! 0:jmaxm values
     >  diff_m,
     >  chie_m,
     >  chii_m,
     >  etaphi_m,
     >  etapar_m,
     >  etaper_m,
     >  exch_m,
     >
     >  egamma_m,      !0:jmaxm exb shear rate in units of local csda_m
     >  egamma_d,      !0:jmaxm exb shear rate delayed by i_delay steps
     >  gamma_p_m,     !0:jmaxm par. vel. shear rate in units of local csda_m
     >  anrate_m,      !0:jmaxm leading mode rate in unints of local csda_m
     >  anrate2_m,     !0:jmaxm 2nd mode rate in units of local csda_m
     >  anfreq_m,      !0:jmaxm leading mode frequency
     >  anfreq2_m      !0:jmaxm 2nd mode frequency
     > )
c************************************************************************
 
c see glf2d documentation for use of diffusivities in transport equations
c ITER definitions of diffusivity used.
c must add neoclassical diffusion
 
c.......begin common block ....................
c      only common used for communication with glf2d....designation by xxxxx_gf

      use glf23_data_mod

      implicit none
 
c.......end common block....................
 
c.......begin dimensions.................
 
      REAL(KIND=rspec) epsilon, zeps, zpi
      parameter ( epsilon = 1.e-34_rspec, zeps = 1.e-6_rspec )
 
c external arrays
 
      integer jmaxm, jshoot, jmm, i_grad, idengrad, itport_pt(1:5),
     >  i_delay, j, jin, jout, jm, irotstab, iglf,
     >  jpt, jptf, jptt, jj, ii, ik, bt_flag, leigen, nroot
      REAL(KIND=rspec) alpha_e, x_alpha, xalpha,
     >  diffnem, chietem, chiitim, etaphim,
     >  etaparm, etaperm, exchm,
     >  rmajor_exp, zimp_exp, amassimp_exp, 
     >  bt_exp, arho_exp, amassgas_exp,
     >  cbetae, cxnu, relx, cmodel, drho, zeff_e,
     >  zpmte, zpmti, zpmne, zpmni, vstar_sign,
     >  egeo_local, pgeo_local, rdrho_local, rdrho_local_p1, fc,
     >  akappa1, alpha_neo, alpha_neo_hold,
     >  zpmnimp, gfac
      REAL(KIND=rspec) te_m(0:jmaxm),ti_m(0:jmaxm),
     >  ne_m(0:jmaxm),ni_m(0:jmaxm), ns_m(0:jmaxm),
     >  vphi_m(0:jmaxm),angrotp_exp(0:jmaxm),
     >  egamma_exp(0:jmaxm),gamma_p_exp(0:jmaxm),
     >  vpar_m(0:jmaxm),vper_m(0:jmaxm),
     >  rho(0:jmaxm),rmin_exp(0:jmaxm),rmaj_exp(0:jmaxm),
     >  gradrho_exp(0:jmaxm),gradrhosq_exp(0:jmaxm),
     >  zeff_exp(0:jmaxm),q_exp(0:jmaxm),shat_exp(0:jmaxm),
     >  bteff_exp(0:jmaxm)
      REAL(KIND=rspec) alpha_exp(0:jmaxm),elong_exp(0:jmaxm),
     >  diff_m(0:jmaxm),chie_m(0:jmaxm),chii_m(0:jmaxm),
     >  etaphi_m(0:jmaxm),
     >  etapar_m(0:jmaxm),etaper_m(0:jmaxm), exch_m(0:jmaxm),
     >  egamma_m(0:jmaxm),egamma_d(0:jmaxm,10),gamma_p_m(0:jmaxm),
     >  anrate_m(0:jmaxm), anrate2_m(0:jmaxm),
     >  anfreq_m(0:jmaxm), anfreq2_m(0:jmaxm)
 
c internal arrays (which can be converted to externals)
 
      REAL(KIND=rspec) zpte_in, zpti_in, zpne_in, zpni_in,
     >  zpte_m(0:jmaxm),zpti_m(0:jmaxm),
     >  zpne_m(0:jmaxm),zpni_m(0:jmaxm),
     >  drhodr(0:jmaxm),drhodrrrho(0:jmaxm),geofac(0:jmaxm),
     >  rhosda_m(0:jmaxm),csda_m(0:jmaxm),cgyrobohm_m(0:jmaxm),
     >  betae_m(0:jmaxm),xnu_m(0:jmaxm),
     >  alpha_m(0:jmaxm),vstarp_m(0:jmaxm)
 
c working arrays and variables
 
      REAL(KIND=rspec) ve(0:jmaxm),vpar(0:jmaxm)
 
c some internals

      REAL(KIND=rspec) tem,tim,nem,nim,nsm,zeffm, aiwt_jp1, 
     >       xnimp_jp1, xnimp, vnewk3x
 
c.......end   dimensions.................
 
c    jm is local grid 0 < jin_m < j < jout_m < jmaxm
c    jm must be greater than 0 and less than jmaxm
c
c
c...constants
      zpi = atan2 ( 0.0e0_rspec, -1.0e0_rspec )
c
c...initialize variables
c
      do j=0,jmaxm
        zpte_m(j) = 0.e0_rspec
        zpti_m(j) = 0.e0_rspec
        zpne_m(j) = 0.e0_rspec
        zpni_m(j) = 0.e0_rspec
 
        betae_m(j)     = 0.e0_rspec
        xnu_m(j)       = 0.e0_rspec
        cgyrobohm_m(j) = 0.e0_rspec
        rhosda_m(j)    = 0.e0_rspec
        csda_m(j)      = 0.e0_rspec
 
        geofac(j)      = 0.e0_rspec
        drhodr(j)      = 0.e0_rspec
        drhodrrrho(j)  = 0.e0_rspec
 
        gamma_p_m(j)   = 0.e0_rspec
        egamma_m(j)    = 0.e0_rspec
        ve(j)          = 0.e0_rspec
        vper_m(j)      = 0.e0_rspec
        vpar(j)        = 0.e0_rspec
        vstarp_m(j)    = 0.e0_rspec
        alpha_m(j)     = 0.e0_rspec
 
        anrate_m(j)    = 0.e0_rspec
        anrate2_m(j)   = 0.e0_rspec
        anfreq_m(j)    = 0.e0_rspec
        anfreq2_m(j)   = 0.e0_rspec
 
        exch_m(j)      = 0.e0_rspec
        diff_m(j)      = 0.e0_rspec
        chie_m(j)      = 0.e0_rspec
        chii_m(j)      = 0.e0_rspec
        etaphi_m(j)    = 0.e0_rspec
        etapar_m(j)    = 0.e0_rspec
        etaper_m(j)    = 0.e0_rspec
      enddo


      eigen_gf = leigen
      nroot_gf=nroot
      iflagin_gf(1)=0
      iflagin_gf(2)=1
      iflagin_gf(3)=1
      iflagin_gf(4)=0
      iflagin_gf(5)=3
 
      xparam_gf(1)=0.e0_rspec
      xparam_gf(2)=0
      xparam_gf(3)=.7e0_rspec
      xparam_gf(4)=0.e0_rspec
      xparam_gf(6)=0.e0_rspec
      xparam_gf(7)=1.e0_rspec
      xparam_gf(8)=0.e0_rspec
      xparam_gf(9)=1.e0_rspec
      xparam_gf(10)=0.e0_rspec
      xparam_gf(11)=0.e0_rspec
      xparam_gf(12)=0.e0_rspec
      xparam_gf(13)=0.2e0_rspec
      xparam_gf(14)=1.e0_rspec
      xparam_gf(15)=-0.1e0_rspec
      xparam_gf(16)=0.e0_rspec
      xparam_gf(17)=0.1e0_rspec
      xparam_gf(18)=.0e0_rspec
      xparam_gf(19)=0.e0_rspec
      xparam_gf(20)=0.e0_rspec
      xparam_gf(21)=0.e0_rspec
      xparam_gf(22)=0.e0_rspec
      xparam_gf(23)=1.e0_rspec
      xparam_gf(24)=0.e0_rspec
      xparam_gf(25)=0.e0_rspec
      xparam_gf(26)=0.e0_rspec
      xparam_gf(27)=0.e0_rspec
      xparam_gf(28)=0.e0_rspec
      xparam_gf(29)=0.e0_rspec
      xparam_gf(30)=0.e0_rspec
c
      xky0_gf= .2e0_rspec
      rms_theta_gf=zpi/3.e0_rspec
      park_gf  =0.7e0_rspec
      ghat_gf  =1.e0_rspec
      gchat_gf =1.e0_rspec
 
      adamp_gf=.50e0_rspec
      alpha_star_gf  =0.e0_rspec
      alpha_mode_gf=0.e0_rspec
      gamma_e_gf  =0.e0_rspec
ctemp
      gamma_e_gf  =-.000000000001e0_rspec
      xkdamp_gf     =0.e0_rspec
 
      alpha_p_gf=0.50e0_rspec
 
c   cbetae=1 is full electromagetic
      cbetae=1.D-6
c      cbetae=1.e0_rspec
c   full collisionality
      cxnu=1.e0_rspec
 
      cnorm_gf=100.e0_rspec
 
      ikymax_gf=10
      xkymin_gf=.02e0_rspec
      xkymax_gf=.5e0_rspec
 
c# non glf23 paramerter
 
       cmodel=1.e0_rspec
       xalpha=x_alpha
 
       bteff_exp = 0.
 
c......begin important optional settings
 
c turn on self-consistant alpha-stabilization
c       ialphastab=1
 
c       turn on EXB shear stabilization
c       alpha_e_gf=1. full on ExB shear
        alpha_e_gf=alpha_e
 
c       turn on self consistant EXB stabilization
c       irotstab=1
 
c itport_pt(1)=1 plasma transport on; itport_pt(2:3)=1; electron and ion heat on
c itport_pt(4)=1 ;itport_pt(5)=0  transport vphi with neoclassical determining vtheta
c itport_pt(4)=1 ;itport_pt(5)=1  transport vphi and vtheta with fast time scale
c       neoclassical drag built into vphi and vtheta transport equations...
c       consult G.M. Staebler
 
c if only vphi_exp is available itport_pt(4)=0
 
c      grid-centering in computing EXB shear  span jm-jptf to jm+jpt
 
c      turn on high-k eta-e modes
        xparam_gf(10)=1.e0_rspec
c
c    relaxation turned off relx=0.
c    relaxation can be turned on for one call per step
       relx=0.e0_rspec
c
c settings for retuned GLF23 model
c
       if (iglf.eq.1) then      ! retuned model
         cnorm_gf=50.e0_rspec         ! ITG normalization (via GYRO runs)
         xparam_gf(10)=12.e0_rspec    ! ETG normalization (cnorm*xparam(10))
         iflagin_gf(5)=5        ! rms theta fit formula
         xparam_gf(13)=0.15     ! rms_theta q-dependence
         xparam_gf(16)=0.15     ! rms_theta shat dependence
         xparam_gf(17)=0.25     ! rms_theta shat dependence
         xparam_gf(19)=1.0      ! rms_theta alpha dependence
         adamp_gf=.70e0_rspec         ! radial mode damping exponent
         alpha_p_gf=0.35e0_rspec      ! parallel velocity shear fit
         park_gf=0.8e0_rspec          ! parallel ion motion fit
         bt_flag=1              ! use real geometry ExB shear
       endif
c
c.......end important optional settings
c*********************************************************************************
c GEOMETRY FACTORS NEEDED FOR SETUP
c external
c      rho(jm)    :toroidal flux co-ordinate 0:50 grids 0->1
c      gradrhosq_exp(jm) : <|grad rho_phys |**2> toroidal flux= B0*rho_phys**2/2
c                 rho_phys=rho*arho_exp
c                 hence   gradrhosq_exp ->1 for a circle
c      gradrho_exp(jm)  : <|grad rho_phys |>
c internal
c      drhodr(jm)
c      drhodrrrho(jm)
c      geofac(jm)
c
c        geofac(j)=gradrho_exp(j)*(rho(j+1)-rho(j))*arho_exp
c     >   /(rmin_exp(j+1)-rmin_exp(j))/gradrhosq_exp(j)
c
c        drhodr(j)=(rho(j+1)-rho(j))*arho_exp/
c     >   (rmin_exp(j+1)-rmin_exp(j))
c
c        drhodrrrho(j)=drhodr(j)*rmin_exp(j)/arho_exp/rho(j)
 
c  surface factor for power flow
c        sfactor(j)=2.*pi_m*arho_exp*rho(j)*h_exp(j)*2.*pi_m*rmajor_exp
c        h_exp(j-1)=hcap_d(j) hcap in ONETWO
 
c******************************************************************************

      if(jmm.gt.0) then
       jin=jmm
       jout=jmm
      endif
      if(jmm.eq.0) then
       jin=1
       jout=jmaxm-1
      endif
      do jm=jin,jout
 
c time dependent codes jshoot=0
c diffusion coefficients and gradients between jm+1 and jm
 
c outside to inside shooting codes jshoot=1 diffusion coefficient at jm
c gradient is forward jm to jm-1
c backward gradient is jm to jm+1
c shear is between forward and backward gradient.
c forward is implicit and backward is already updated

       tem   = te_m(jm)
       tim   = ti_m(jm)
       nem   = ne_m(jm)
       nim   = ni_m(jm)
       nsm   = ns_m(jm)
       zeffm = zeff_exp(jm)
 
       betae_m(jm) = 400.e0_rspec*nem*tem/(1.D5*bt_exp**2)
 
crew    gks collisionality (xnu/w_star_i)*(ky*rho_i)
       vnewk3x=
     >   0.117e0_rspec*nem*tem**(-1.5e0_rspec)/(tim**0.5e0_rspec)
     >   *arho_exp*(amassgas_exp/2.e0_rspec)**0.5e0_rspec

crew   as used in gks multiply by 1/2 and takout any 1/2 factor in solfp
crew          vnewk3x=vnewk3x/2.

       xnu_m(jm) =vnewk3x/(2.e0_rspec*tem/tim)**0.5e0_rspec

crew  10/25/95 fixed zeff+1 factor: zeff col with ions;1 col with elecs.
       zeff_e=0.e0_rspec
       xnu_m(jm) = xnu_m(jm)*(zeff_exp(jm)+zeff_e)

 
      if(jm.eq.0.or.jm.eq.jmaxm) then
       write(6,*) 'can not call callglf2d for this jm'
      endif
 
      jptf=0
       if(jm.eq.1) jptf=0
      jpt=1
      jptt=jpt
       if(jptt.lt.1) jptt=1
       if(jm.eq.jmaxm-1) jptt=1
      jpt=jptt
 
      do j=jm-jptf,jm+jpt
c
c... some geometric factors
c
        geofac(j)=gradrho_exp(j)*(rho(j)-rho(j-1))*arho_exp
     >   /(rmin_exp(j)-rmin_exp(j-1)+epsilon)/gradrhosq_exp(j)
 
        drhodr(j)=(rho(j)-rho(j-1))*arho_exp/
     >   (rmin_exp(j)-rmin_exp(j-1)+epsilon)
 
        drhodrrrho(j)=drhodr(j)*rmin_exp(j)/
     >   arho_exp/(rho(j)+epsilon)
c
c... local rate unit
c
       csda_m(j)=9.79D5*(te_m(j)*1.D3)**.5e0_rspec/
     >    (arho_exp*100.e0_rspec)/amassgas_exp**0.5e0_rspec
c
c... local rho_star
c Note: use effective B-field if bt_flag > 0
c
       if (bt_flag .gt. 0) then
         bteff_exp(j)=bt_exp*rho(j)*arho_exp/
     >         rmin_exp(j)*drhodr(j)
         rhosda_m(j)=((1.02D2*(te_m(j)*1.D3)**.5e0_rspec)/bteff_exp(j)
     >         /1.D4)*amassgas_exp**.5e0_rspec/(arho_exp*100.e0_rspec)
       else
         rhosda_m(j)=((1.02D2*(te_m(j)*1.D3)**.5e0_rspec)/bt_exp/1.D4)
     >         *amassgas_exp**.5e0_rspec/(arho_exp*100.e0_rspec)
       endif
      enddo
c
c   local gyrobohm unit of diffusion
c
       cgyrobohm_m(jm)=1.D-4*
     >  9.79D5*(tem*1.D3)**.5e0_rspec/(arho_exp*100.e0_rspec)
     >  *(1.02D2*(tem*1.D3)**.5e0_rspec/bt_exp/1.D4)**2*amassgas_exp**.5e0_rspec
c 
       zpmte=zpte_in
       zpmti=zpti_in
       zpmne=zpne_in
       zpmni=zpni_in
 
c MHD alpha parameter
 
         alpha_m(jm)=drhodr(jm)*
     >    q_exp(jm)**2*rmaj_exp(jm)/arho_exp*
     >    betae_m(jm)*((tim/tem*nim/nem)*
     >   (zpmni+zpmti)
     >    +zpmne+zpmte)
 
c vstarp_m is diamagnetic part of egamma_m (doppler shear rate)
c vstar_sign is negative for negative vstar_i. Thus for co-injection or positive
c angrot toroidal rotation cancels the diamgnetic rotation
 
       vstar_sign=-1.e0_rspec
 
        j=jm
        rho(j-jptf)=rho(j-jptf)+epsilon
        rho(j+jpt)=rho(j+jpt)+epsilon
        egeo_local=1.e0_rspec
        pgeo_local=drhodr(j)
        rdrho_local=rmin_exp(j-jptf)/arho_exp/rho(j-jptf)
        rdrho_local_p1=rmin_exp(j+jpt)/arho_exp/rho(j+jpt)
 
        vstarp_m(jm)=(
     >     +egeo_local*vstar_sign*
     >     (rho(j+jpt)*rdrho_local_p1+rho(j-jptf)*rdrho_local)
     >     /2.e0_rspec*(
     >     (ti_m(j+jpt)/te_m(j+jpt))*csda_m(j+jpt)*
     >     (zpti_m(j+jpt)+zpni_m(j+jpt))
     >     *pgeo_local*rhosda_m(j+jpt)/rho(j+jpt)/rdrho_local_p1
     >     -(ti_m(j-jptf)/te_m(j-jptf))*csda_m(j-jptf)*
     >     (zpti_m(j-jptf)+zpni_m(j-jptf))
     >     *pgeo_local*rhosda_m(j-jptf)/rho(j-jptf)/rdrho_local
     >     )/(rho(j+jpt)-rho(j-jptf)+epsilon)/csda_m(j)
     >     )
c
       do jj=1,2
 
        if(jj.eq.1) j=jm-jptf
 
        if(jj.eq.2) j=jm+jpt

c banana regime ie collisionless limit formulas
 
        fc=1-1.46e0_rspec*(rmin_exp(j)/rmaj_exp(j))**0.5e0_rspec+
     >      0.46e0_rspec*(rmin_exp(j)/rmaj_exp(j))**1.5e0_rspec
        akappa1=0.8839e0_rspec*fc/(0.3477e0_rspec+0.4058e0_rspec*fc)
        alpha_neo=-akappa1+1.e0_rspec
        alpha_neo_hold=alpha_neo
 
        egeo_local=1.e0_rspec
        pgeo_local=drhodr(j)
        rdrho_local=rmin_exp(j)/arho_exp/(rho(j)+epsilon)
 
        ve(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     > (zpni_m(j)+alpha_neo*zpti_m(j))*vstar_sign*pgeo_local
     > -rho(j)*rdrho_local*
     >  arho_exp/rmajor_exp/q_exp(j)*rmajor_exp*angrotp_exp(j)

 
        vpar(j)=rmajor_exp*angrotp_exp(j)-vstar_sign*
     > (ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*pgeo_local*
     >((alpha_neo-1.e0_rspec)*zpti_m(j))*rho(j)*rdrho_local*
     > arho_exp/rmajor_exp/q_exp(j)
 
 
       if(itport_pt(4).eq.0.and.itport_pt(5).eq.0) then
         vphi_m(j)=rmajor_exp*angrotp_exp(j)
 
         vpar_m(j)=vpar(j)
 
         vper_m(j)=ve(j)
     >  +(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     >   (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
       endif
 
       if(abs(itport_pt(4)).eq.1.and.itport_pt(5).eq.1) then
         vpar(j)=vpar_m(j)
       endif
 
       if(abs(itport_pt(4)).eq.1.and.itport_pt(5).eq.0) then
c this option vpar is vphi vexb from neo+vphi
         ve(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     >  (zpni_m(j)+alpha_neo*zpti_m(j))*vstar_sign*pgeo_local
     >  -rho(j)*rdrho_local*
     >   arho_exp/rmajor_exp/q_exp(j)*vphi_m(j)
 
         vpar(j)=vphi_m(j)-vstar_sign*
     >   (ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*pgeo_local*
     >   ((alpha_neo-1.e0_rspec)*zpti_m(j))*rho(j)*rdrho_local*
     >   arho_exp/rmajor_exp/q_exp(j)
 
         vpar_m(j)=vpar(j)
 
         vper_m(j)=ve(j)
     >   +(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     >   (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
       endif
 
       if(itport_pt(5).eq.1) then
c this option vexb from vper and vpar with neo dampng built into
c vpar and vper transport equations
         ve(j)=vper_m(j)
     >   -(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     >   (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
       endif
 
      enddo
c 
c compute shears from outside minus inside
c 
        j=jm
        rho(j-jptf)=rho(j-jptf)+epsilon
        rho(j+jpt)=rho(j+jpt)+epsilon
        egeo_local=1.e0_rspec
        pgeo_local=drhodr(j)
        rdrho_local=rmin_exp(j-jptf)/arho_exp/rho(j-jptf)
        rdrho_local_p1=rmin_exp(j+jpt)/arho_exp/rho(j+jpt)
c
        egamma_m(jm)=relx*egamma_m(jm)+(1.e0_rspec-relx)*(
     >  egeo_local*drhodrrrho(j)*
     >  (rho(j-jptf)+rho(j+jpt))/(q_exp(j-jptf)+q_exp(j+jpt))*
     >  (ve(j+jpt)*q_exp(j+jpt)/rho(j+jpt)/rdrho_local_p1-
     >  ve(j-jptf)*q_exp(j-jptf)/rho(j-jptf)/rdrho_local)/
     >  (rho(j+jpt)-rho(j-jptf)+epsilon)/arho_exp/csda_m(j)
     >  )

         
        gamma_p_m(jm)=relx*gamma_p_m(jm)+(1.e0_rspec-relx)*(
     >   -drhodr(j)*
     >   (vpar(j+jpt)-vpar(j-jptf))
     >   /(rho(j+jpt)-rho(j-jptf)+epsilon)/arho_exp/csda_m(j)
     >  )
 

       if (jm.eq.1.and.gamma_p_m(jm).gt.10)
     >    gamma_p_m(jm)=10.e0_rspec
        alpha_neo=alpha_neo_hold
 
c.......end   setups...........
 
c***********************************************************************
cvv      subroutine model
************************************************************************
 
c   units:
c     diffusion (m**2/sec) note: computed in derived modprofiles
c     density (10**13 cm**-3 or 10**19 m**-3)
c     arho_exp and rmajor_exp (m)
c     power (MW)
c     flow (MW/kev=kA)
c
c   kev/sec per MW
c   kevdsecpmw=1.6022e-19*1.0e3*1.e-6
c
c       cgyrobohm_m(jm)=1.e-4*
c    >  9.79e5*(tem*1.e3)**.5/(arho_exp*100.)
c    >  *(1.02e2*(tem*1.e3)**.5/bt_exp/1.e4)**2*(amassgas_exp)**.5
c
c   sfactor(j)=
c     >      2.*pi_m*arho_exp*rho(j)*h_exp(j)*2.*pi_m*rmajor_exp
c   units of meter**2
c
c   9000 format (1i6,7e10.3)
 
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cmnt                         the model
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cmnt
cmnt supply model for chietem,chietim,chienem
cmnt                  chiitem,chiitim,chiinem
cmnt                  difftem,difftim,diffnem
cmnt
cmnt    chi_s and diff_s must be in meters**2/sec units
cmnt     and recall chi_s refer to total energy flow
cmnt
cmnt    if the model chi_s refer to "heat conduction" flow
cmnt    then a convection term xconv*3./2.*t_m*flow_exp is added.
cmnt    normally input xconv=0. otherwise xconv=1. or 5./3.
cmnt
cmnt    it is also possible to build convection into the model
cmnt    with "aconv".  aconv and xconv should not be double counted.
cmnt
cmnt note: can use diagonal forms with off diagonal dependence on
cmnt zpmte,zpmti,zpmne intrinsic to the diagonal elements as in sample
cmnt normall models are written in diagonal for with dependence on off
cmnt diagonal gradients implicit
cmnt
cmnt note: when flow is large anomalous e-i exchange should be added
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cxx      if (imodel.eq.8) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 2DGLF quasilinear model  GLF23
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
       rmaj_gf  = rmaj_exp(jm)/arho_exp
       rmin_gf  = rmin_exp(jm)/arho_exp
       q_gf     = q_exp(jm)
       betae_gf = dmax1(cbetae*betae_m(jm), 1.D-6)
       shat_gf  = shat_exp(jm)*drhodrrrho(jm)
 
       if(xalpha.lt.0.) alpha_gf=-xalpha*alpha_m(jm)
       if(xalpha.gt.0.) alpha_gf=xalpha*alpha_exp(jm)
       if(alpha_gf.gt.4.e0_rspec) alpha_gf=4.e0_rspec
 
       elong_gf=elong_exp(jm)
 
       xnu_gf=cxnu*xnu_m(jm)
       taui_gf=tim/tem
       amassgas_gf=amassgas_exp
 
       apwt_gf=1.e0_rspec
c impurity dynamics not turned on by default 
c and simple dilution included (idengrad=2, dil_gf=1-nim/nem)
c to turn on impurity dynamics need to change number of roots
c supply zimp_exp, amassimp_exp, and fractional density weights
c apwt_gf and aiwt_gf
       dil_gf      = 0.e0_rspec
       aiwt_gf     = 0.e0_rspec
       zimp_gf     = zimp_exp
       amassimp_gf = amassimp_exp
       rlnimp_gf   = 1.e0_rspec
       zpmnimp     = 1.e0_rspec
       if(idengrad.eq.2) dil_gf=1.e0_rspec-nim/nem
       if(idengrad.eq.3) then
         apwt_gf=nim/nem
         aiwt_jp1=(zeffm*nem-ni_m(jm+1)
     >            -ns_m(jm+1))/(zimp_gf**2*nem)
         xnimp_jp1=aiwt_jp1*ne_m(jm+1)
         aiwt_gf=(zeffm*nem
     >           -nim-ns_m(jm))/(zimp_gf**2*ne_m(jm))
         xnimp=aiwt_gf*ne_m(jm)
         zpmnimp=-(dlog(xnimp_jp1)-dlog(xnimp))/
     >           (rho(jm+1)-rho(jm))
         rlnimp_gf=zpmnimp*elong_exp(jm)**0.5
       endif
 
       rlte_gf   = zpmte   
       rlti_gf   = zpmti   
       rlne_gf   = zpmne 
       rlni_gf   = zpmni 
       rlnimp_gf = rmaj_exp(jm)*zpmnimp*drhodr(jm)
 
       gamma_star_gf = vstarp_m(jm)
       gamma_e_gf    = egamma_m(jm)
       gamma_p_gf    = gamma_p_m(jm)
       if(itport_pt(4).eq.-1) then
         gamma_e_gf = egamma_exp(jm)
         gamma_p_gf = gamma_p_exp(jm)
       endif

c..jek 8/15/00
        if (irotstab.eq.0) then
          gamma_e_gf=egamma_exp(jm)
          gamma_p_gf=gamma_p_exp(jm)
          egamma_m(jm)=egamma_exp(jm)
        endif
        gamma_mode_gf=0.0e0_rspec
 
        if(i_delay.ne.0) then
 
         if(i_delay.gt.1) then
          do ii=1,i_delay-1
           egamma_d(jm,ii)=egamma_d(jm,ii+1)
          enddo
 
          egamma_d(jm,i_delay)=egamma_m(jm)
         endif
          gamma_e_gf=egamma_d(jm,1)
        endif
 
 
c.......THE  CALL TO GLF23
 
        call glf2d(iglf)
 
c.......POST CALL TO GLF23
 
       anrate_m(jm)=gamma_gf(1)
       anrate2_m(jm)=gamma_gf(2)
       anfreq_m(jm)=freq_gf(1)
       anfreq2_m(j)=freq_gf(2)

       gfac=geofac(jm)
c
c exch_m in MW/m**3
c   kev/sec per MW
ck       kevdsecpmw=1.6022e-19*1.0e3*1.e-6
 
      exch_m(jm)=1.D19*1.6022D-19*1.0D3*1.D-6*
     > nem*tem*csda_m(jm)*rhosda_m(jm)**2*exch_gf*cmodel
 
      exchm=exch_m(jm)
 
      chietem=cmodel*gfac*chie_gf*cgyrobohm_m(jm)
      chiitim=cmodel*gfac*chii_gf*cgyrobohm_m(jm)
      diffnem=cmodel*gfac*diff_gf*cgyrobohm_m(jm)
 
      etaphim=cmodel*gfac*eta_phi_gf*cgyrobohm_m(jm)
      etaparm=cmodel*gfac*eta_par_gf*cgyrobohm_m(jm)
      etaperm=cmodel*gfac*eta_per_gf*cgyrobohm_m(jm)
 
      if(itport_pt(1) .eq. 0 ) diffnem = 0.e0_rspec
      if(itport_pt(2) .eq. 0 ) chietem = 0.e0_rspec
      if(itport_pt(3) .eq. 0 ) chiitim = 0.e0_rspec
      if((itport_pt(4) .eq. 0) .and. (itport_pt(5) .eq. 0) ) then
         etaphim=0.e0_rspec
         etaparm=0.e0_rspec
         etaperm=0.e0_rspec
      endif

      diff_m(jm)=diffnem
      chie_m(jm)=chietem
      chii_m(jm)=chiitim
 
      etaphi_m(jm)=etaphim
      etapar_m(jm)=etaparm
      etaper_m(jm)=etaperm
 
        enddo
 
       return
       end
