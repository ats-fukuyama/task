c@callglf2d.f
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
     > i_grad,         ! default 0, for D-V method use i_grad=1 to input gradients
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
 
      implicit none

      include 'glf.m'
 
c.......end common block....................
 
c.......begin dimensions.................
 
      double precision epsilon, zeps, zpi
      parameter ( epsilon = 1.D-34, zeps = 1.D-6 )
 
c external arrays
 
      integer jmaxm, jshoot, jmm, i_grad, idengrad, itport_pt(1:5),
     >  i_delay, j, jin, jout, jm, irotstab, iglf,
     >  jpt, jptf, jptt, jj, ii, ik, bt_flag, leigen, nroot
      double precision alpha_e, x_alpha, xalpha,
     >  diffnem, chietem, chiitim, etaphim,
     >  etaparm, etaperm, exchm,
     >  rmajor_exp, zimp_exp, amassimp_exp, 
     >  bt_exp, arho_exp, amassgas_exp,
     >  cbetae, cxnu, relx, cmodel, drho, zeff_e,
     >  zpmte, zpmti, zpmne, zpmni, vstar_sign,
     >  egeo_local, pgeo_local, rdrho_local, rdrho_local_p1, fc,
     >  akappa1, alpha_neo, alpha_neo_hold,
     >  zpmne_q, zpmni_q, zpmnimp, gfac
      double precision te_m(0:jmaxm),ti_m(0:jmaxm),
     >  ne_m(0:jmaxm),ni_m(0:jmaxm), ns_m(0:jmaxm),
     >  vphi_m(0:jmaxm),angrotp_exp(0:jmaxm),
     >  egamma_exp(0:jmaxm),gamma_p_exp(0:jmaxm),
     >  vpar_m(0:jmaxm),vper_m(0:jmaxm),
     >  rho(0:jmaxm),rmin_exp(0:jmaxm),rmaj_exp(0:jmaxm),
     >  gradrho_exp(0:jmaxm),gradrhosq_exp(0:jmaxm),
     >  zeff_exp(0:jmaxm),q_exp(0:jmaxm),shat_exp(0:jmaxm),
     >  bteff_exp(0:jmaxm)
      double precision alpha_exp(0:jmaxm),elong_exp(0:jmaxm),
     >  diff_m(0:jmaxm),chie_m(0:jmaxm),chii_m(0:jmaxm),
     >  etaphi_m(0:jmaxm),
     >  etapar_m(0:jmaxm),etaper_m(0:jmaxm), exch_m(0:jmaxm),
     >  egamma_m(0:jmaxm),egamma_d(0:jmaxm,10),gamma_p_m(0:jmaxm),
     >  anrate_m(0:jmaxm), anrate2_m(0:jmaxm),
     >  anfreq_m(0:jmaxm), anfreq2_m(0:jmaxm)
 
c internal arrays (which can be converted to externals)
 
      double precision zpte_in, zpti_in, zpne_in, zpni_in,
     >  zpte_m(0:jmaxm),zpti_m(0:jmaxm),
     >  zpne_m(0:jmaxm),zpni_m(0:jmaxm),
     >  drhodr(0:jmaxm),drhodrrrho(0:jmaxm),geofac(0:jmaxm),
     >  rhosda_m(0:jmaxm),csda_m(0:jmaxm),cgyrobohm_m(0:jmaxm),
     >  betae_m(0:jmaxm),xnu_m(0:jmaxm),
     >  alpha_m(0:jmaxm),vstarp_m(0:jmaxm)
 
c working arrays and variables
 
      double precision ve(0:jmaxm),vpar(0:jmaxm)
c     real*8 vmode(0:jmaxm)
c     real*8 kevdsecpmw
 
c diagnostic arrays (if used)
 
c     real*8 vstar_m(0:jmaxm),vexb_m(0:jmaxm)
c     real*8 vmode_m(0:jmaxm)
c     real*8 gamma_mode_m(0:jmaxm)
c     real*8 gamma_k_j(20,0:jmaxm),freq_k_j(20,0:jmaxm)
c     real*8 chie_k_j(20,0:jmaxm),chii_k_j(20,0:jmaxm)
c     real*8 vnewstare_m(0:jmaxm),vnewstari_m(0:jmaxm)
c     real*8 ky_j(0:jmaxm)
c     real*8 gamma_j(0:jmaxm,1:4),freq_j(0:jmaxm,1:4)
c     real*8 phi_norm_j(0:jmaxm,1:4)
c     real*8 dnrate_m(0:jmaxm), dtnrate_m(0:jmaxm)
c     real*8 dnfreq_m(0:jmaxm)
 
c some internals

      double precision tem,tim,nem,nim,nsm,zeffm, aiwt_jp1, 
     >       xnimp_jp1, xnimp, vnewk3x
 
c.......end   dimensions.................
 
c    jm is local grid 0 < jin_m < j < jout_m < jmaxm
c    jm must be greater than 0 and less than jmaxm
c
c
c...constants
      zpi = atan2 ( 0.0D0, -1.0D0 )
c
cdmc      write(6,7701) zpi
cdmc 7701 format(' zpi = ',1pe17.10)
cdmc      write(6,7702) jshoot,jmm,jmaxm,(itport_pt(j),j=1,5)
cdmc 7702 format(/' jshoot,jmm,jmaxm = ',3(1x,i6)/' itport_pt = ',5(1x,i3))
cdmc      call echo('te_m',te_m,jmaxm)
cdmc      call echo('ti_m',ti_m,jmaxm)
cdmc      call echo('ne_m',ne_m,jmaxm)
cdmc      call echo('ni_m',ni_m,jmaxm)
c
cdmc      write(6,7703) i_grad,zpte_in,zpti_in,zpne_in,zpni_in
cdmc 7703 format(/' igz: ',i1,1x,4(1pe17.10))
c
cdmc      call echo('angrotp_exp',angrotp_exp,jmaxm)
cdmc      call echo('egamma_exp',egamma_exp,jmaxm)
cdmc      call echo('gamma_p_exp',gamma_p_exp,jmaxm)
cdmc      call echo('vphi_m',vphi_m,jmaxm)
cdmc      call echo('vpar_m',vpar_m,jmaxm)
cdmc      call echo('vper_m',vper_m,jmaxm)
cdmc      call echo('zeff_exp',zeff_exp,jmaxm)
c
cdmc      write(6,7705) bt_exp
cdmc 7705 format(/' bt_exp = ',1pe17.10)
c
cdmc      call echo('rho',rho,jmaxm)
c
cdmc      write(6,7706) arho_exp
cdmc 7706 format(/' arho_exp = ',1pe17.10)
c
cdmc      call echo('gradrho_exp',gradrho_exp,jmaxm)
cdmc      call echo('gradrhosq_exp',gradrhosq_exp,jmaxm)
cdmc      call echo('rmin_exp',rmin_exp,jmaxm)
cdmc      call echo('rmaj_exp',rmaj_exp,jmaxm)
c
cdmc      write(6,7707) rmajor_exp
cdmc 7707 format(/' rmajor_exp = ',1pe17.10)
c
cdmc      call echo('q_exp',q_exp,jmaxm)
cdmc      call echo('shat_exp',shat_exp,jmaxm)
cdmc      call echo('alpha_exp',alpha_exp,jmaxm)
cdmc      call echo('elong_exp',elong_exp,jmaxm)
c
cdmc      write(6,7708) amassgas_exp
cdmc 7708 format(/' atomic no.:  ',1pe17.10)
c
cdmc      write(6,7709) alpha_e,x_alpha
cdmc 7709 format(/' alpha_e, x_alpha = ',2(1x,1pe17.10))
c
cdmc      write(6,7710) i_delay
cdmc 7710 format(/' i_delay = ',i5)
c
c...initialize variables
c
      do j=0,jmaxm
        zpte_m(j) = 0.D0
        zpti_m(j) = 0.D0
        zpne_m(j) = 0.D0
        zpni_m(j) = 0.D0
 
        betae_m(j)     = 0.D0
        xnu_m(j)       = 0.D0
        cgyrobohm_m(j) = 0.D0
        rhosda_m(j)    = 0.D0
        csda_m(j)      = 0.D0
 
        geofac(j)      = 0.D0
        drhodr(j)      = 0.D0
        drhodrrrho(j)  = 0.D0
 
        gamma_p_m(j)   = 0.D0
        egamma_m(j)    = 0.D0
        ve(j)          = 0.D0
        vper_m(j)      = 0.D0
        vpar(j)        = 0.D0
c        vphi_m(j)      = 0.D0
        vstarp_m(j)    = 0.D0
        alpha_m(j)     = 0.D0
 
        anrate_m(j)    = 0.D0
        anrate2_m(j)   = 0.D0
        anfreq_m(j)    = 0.D0
        anfreq2_m(j)   = 0.D0
 
        exch_m(j)      = 0.D0
        diff_m(j)      = 0.D0
        chie_m(j)      = 0.D0
        chii_m(j)      = 0.D0
        etaphi_m(j)    = 0.D0
        etapar_m(j)    = 0.D0
        etaper_m(j)    = 0.D0
      enddo
c
c diagnostic arrays (not used)
c
c     do j=0,jmaxm
c       vstar_m(j)     = 0.
c       vexb_m(j)      = 0.
c       dnrate_m(j)    = 0.
c       dtnrate_m(j)   = 0.
c       dnfreq_m(j)    = 0.
c      do k=1,4
c       gamma_j(j,k)    = 0.
c       freq_j(j,k)     = 0.
c       phi_norm_j(j,k) = 0.
c      enddo
c     enddo
c
c     do j=0,jmaxm
c      do ik=1,20
c       gamma_k_j(ik,j) = 0.D0
c       freq_k_j(ik,j)  = 0.D0
c       chie_k_j(ik,j)  = 0.D0
c       chii_k_j(ik,j)  = 0.D0
c      enddo
c      vnewstare_m(j) = 0.
c      vnewstari_m(j) = 0.
c     enddo
c
c
************************************************************************
cmnt   profiles of quantities derived from model profiles of te,ti,ne
cmnt   revised derivedmodlocal to center between jm-jptf and jm+ptf
************************************************************************
 
c.......begin switches and settings......
 
c**********************************************************************
c     glf23 parameters
cxx      settings same as function glf23_v4_1_10
 
ck      write(6,*)  'jmm=',jmm,'going in'

      eigen_gf = leigen
c      nroot_gf=8     ! 8 for pure plasma, 12 for full impurity dynamics
      nroot_gf=nroot
      iflagin_gf(1)=0
      iflagin_gf(2)=1
      iflagin_gf(3)=1
      iflagin_gf(4)=0
      iflagin_gf(5)=3
 
      xparam_gf(1)=0.D0
      xparam_gf(2)=0
      xparam_gf(3)=.7D0
      xparam_gf(4)=0.D0
      xparam_gf(6)=0.D0
      xparam_gf(7)=1.D0
      xparam_gf(8)=0.D0
      xparam_gf(9)=1.D0
      xparam_gf(10)=0.D0
      xparam_gf(11)=0.D0
      xparam_gf(12)=0.D0
      xparam_gf(13)=0.2D0
      xparam_gf(14)=1.D0
      xparam_gf(15)=-0.1D0
      xparam_gf(16)=0.D0
      xparam_gf(17)=0.1D0
      xparam_gf(18)=.0D0
      xparam_gf(19)=0.D0
      xparam_gf(20)=0.D0
      xparam_gf(21)=0.D0
      xparam_gf(22)=0.D0
      xparam_gf(23)=1.D0
      xparam_gf(24)=0.D0
      xparam_gf(25)=0.D0
      xparam_gf(26)=0.D0
      xparam_gf(27)=0.D0
      xparam_gf(28)=0.D0
      xparam_gf(29)=0.D0
      xparam_gf(30)=0.D0
c
      xky0_gf= .2D0
      rms_theta_gf=zpi/3.D0
      park_gf  =0.7D0
      ghat_gf  =1.D0
      gchat_gf =1.D0
 
      adamp_gf=.50D0
      alpha_star_gf  =0.D0
      alpha_mode_gf=0.D0
      gamma_e_gf  =0.D0
ctemp
      gamma_e_gf  =-.000000000001D0
      xkdamp_gf     =0.D0
 
      alpha_p_gf=0.50D0
 
c   cbetae=1 is full electromagetic
      cbetae=1.D-6
c      cbetae=1.D0
c   full collisionality
      cxnu=1.D0
 
      cnorm_gf=100.D0
 
      ikymax_gf=10
      xkymin_gf=.02D0
      xkymax_gf=.5D0
 
c# non glf23 paramerter
 
       cmodel=1.D0
       xalpha=x_alpha
 
c      ialphastab=1
c      ineo=-2
 
c      iexch_m=1
c      iexp_exch=-1
 
c      i_dengrad=2
c      iexp_imp=1
c      igeo_m=3
 
c      irotstab=1
c      irot1=1
c      irot2=1
c
cxx      endf
 
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
        xparam_gf(10)=1.D0
c
c    relaxation turned off relx=0.
c    relaxation can be turned on for one call per step
       relx=0.D0
c
c settings for retuned GLF23 model
c
       if (iglf.eq.1) then      ! retuned model
         cnorm_gf=50.D0         ! ITG normalization (via GYRO runs)
         xparam_gf(10)=12.D0    ! ETG normalization (cnorm*xparam(10))
         iflagin_gf(5)=5        ! rms theta fit formula
         xparam_gf(13)=0.15     ! rms_theta q-dependence
         xparam_gf(16)=0.15     ! rms_theta shat dependence
         xparam_gf(17)=0.25     ! rms_theta shat dependence
         xparam_gf(19)=1.0      ! rms_theta alpha dependence
         adamp_gf=.70D0         ! radial mode damping exponent
         alpha_p_gf=0.35D0      ! parallel velocity shear fit
         park_gf=0.8D0          ! parallel ion motion fit
         bt_flag=1              ! use real geometry ExB shear
       endif
c
c.......end important optional settings
 
c.......end switches and settings......
 
c.......start setups.....................
c
 
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
 
       tem=(te_m(jm+1-jshoot)+te_m(jm))/2.D0
       tim=(ti_m(jm+1-jshoot)+ti_m(jm))/2.D0
       nem=(ne_m(jm+1-jshoot)+ne_m(jm))/2.D0
       nim=(ni_m(jm+1-jshoot)+ni_m(jm))/2.D0
       nsm=(ns_m(jm+1-jshoot)+ns_m(jm))/2.D0
       zeffm=(zeff_exp(jm+1-jshoot)+zeff_exp(jm))/2.D0
       betae_m(jm) = 400.D0*nem*tem/(1.D5*bt_exp**2)
c      betai_m(jm) = 400.*nim*tim/(1.e5*bt_exp**2)
 
 
crew    gks collisionality (xnu/w_star_i)*(ky*rho_i)
       vnewk3x=
     >   0.117D0*nem*tem**(-1.5D0)/(tim**0.5D0)*arho_exp*
     >   (amassgas_exp/2.D0)**0.5D0
crew   as used in gks multiply by 1/2 and takout any 1/2 factor in solfp
crew          vnewk3x=vnewk3x/2.
       xnu_m(jm) =vnewk3x/(2.D0*tem/tim)**0.5D0
crew  10/25/95 fixed zeff+1 factor: zeff col with ions;1 col with elecs.
       zeff_e=0.D0
       xnu_m(jm) = xnu_m(jm)*(zeff_exp(jm)+zeff_e)
 
c      vnewstare_m(jm)=zeff_exp(jm) *2.91e-6*nem*1.e13*15./
c    >  (tem*1.e3)**2*rmaj_exp(jm)*100.*q_exp(jm)
c    >     /(rmin_exp(jm)/rmaj_exp(jm)+epsilon)**1.5/4.19e7
c
c      vnewstari_m(jm)=4.78e-8*nem*1.e13*15./
c    >  (tim*1.e3)**2*rmaj_exp(jm)*100.*q_exp(jm)
c    >     /(rmin_exp(jm)/rmaj_exp(jm)+epsilon)**1.5/9.79e5
 
c
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
c  note: dependence of shear's  on zpxx_m(jm-1) and zpxx_m(jm+1)
c  and alpha(jm) depends on zpxx_m(jm+1)
 200  format(i2,2x,0p1f4.2,0p5f10.5)
c
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
       csda_m(j)=9.79D5*(te_m(j)*1.D3)**.5D0/
     >    (arho_exp*100.D0)/amassgas_exp**0.5D0
c
c... local rho_star
c Note: use effective B-field if bt_flag > 0
c
       if (bt_flag .gt. 0) then
         bteff_exp(j)=bt_exp*rho(j)*arho_exp/
     >         rmin_exp(j)*drhodr(j)
         rhosda_m(j)=((1.02D2*(te_m(j)*1.D3)**.5D0)/bteff_exp(j)
     >         /1.D4)*amassgas_exp**.5D0/(arho_exp*100.D0)
       else
         rhosda_m(j)=((1.02D2*(te_m(j)*1.D3)**.5D0)/bt_exp/1.D4)
     >         *amassgas_exp**.5D0/(arho_exp*100.D0)
       endif
      enddo
c
c   local gyrobohm unit of diffusion
c
       cgyrobohm_m(jm)=1.D-4*
     >  9.79D5*(tem*1.D3)**.5D0/(arho_exp*100.D0)
     >  *(1.02D2*(tem*1.D3)**.5D0/bt_exp/1.D4)**2*amassgas_exp**.5D0
c 
      do j=jm-jptf, jm+jpt
        drho=rho(j-1)-rho(j)+epsilon
        zpte_m(j)=-(dlog(te_m(j-1))-dlog(te_m(j)))/drho
        zpti_m(j)=-(dlog(ti_m(j-1))-dlog(ti_m(j)))/drho
        zpne_m(j)=-(dlog(ne_m(j-1))-dlog(ne_m(j)))/drho
        zpni_m(j)=-(dlog(ni_m(j-1))-dlog(ni_m(j)))/drho
      enddo

        zpmte=zpte_m(jm+1-jshoot)
        zpmti=zpti_m(jm+1-jshoot)
        zpmne=zpne_m(jm+1-jshoot)
        zpmni=zpni_m(jm+1-jshoot)
c
c... check on zero norm gradients:
c
        if (DABS(zpmti).lt.zeps) zpmti=zeps
        if (DABS(zpmte).lt.zeps) zpmte=zeps
        if (DABS(zpmne).lt.zeps) zpmne=zeps
        if (DABS(zpmni).lt.zeps) zpmni=zeps
c
        if(i_grad.eq.1) then
         zpmte=zpte_in
         zpmti=zpti_in
         zpmne=zpne_in
         zpmni=zpni_in
        endif
 
c MHD alpha parameter
 
         alpha_m(jm)=drhodr(jm)*
     >    q_exp(jm)**2*rmaj_exp(jm)/arho_exp*
     >    betae_m(jm)*((tim/tem*nim/nem)*
     >   (zpmni+zpmti)
     >    +zpmne+zpmte)
 
c vstarp_m is diamagnetic part of egamma_m (doppler shear rate)
c vstar_sign is negative for negative vstar_i. Thus for co-injection or positive
c angrot toroidal rotation cancels the diamgnetic rotation
 
       vstar_sign=-1.D0
 
        j=jm
          rho(j-jptf)=rho(j-jptf)+epsilon
          rho(j+jpt)=rho(j+jpt)+epsilon
          egeo_local=1.D0
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j-jptf)/arho_exp/rho(j-jptf)
          rdrho_local_p1=rmin_exp(j+jpt)/arho_exp/rho(j+jpt)
 
        vstarp_m(jm)=(
     > +egeo_local*vstar_sign*
     > (rho(j+jpt)*rdrho_local_p1+rho(j-jptf)*rdrho_local)/2.D0*(
     >(ti_m(j+jpt)/te_m(j+jpt))*csda_m(j+jpt)*
     > (zpti_m(j+jpt)+zpni_m(j+jpt))
     >  *pgeo_local*rhosda_m(j+jpt)/rho(j+jpt)/rdrho_local_p1
     >-(ti_m(j-jptf)/te_m(j-jptf))*csda_m(j-jptf)*
     > (zpti_m(j-jptf)+zpni_m(j-jptf))
     >  *pgeo_local*rhosda_m(j-jptf)/rho(j-jptf)/rdrho_local
     >  )/(rho(j+jpt)-rho(j-jptf)+epsilon)/csda_m(j)
     >  )
c
       do jj=1,2
 
        if(jj.eq.1) j=jm-jptf
 
        if(jj.eq.2) j=jm+jpt

c banana regime ie collisionless limit formulas
 
        fc=1-1.46D0*(rmin_exp(j)/rmaj_exp(j))**0.5D0+
     >      0.46D0*(rmin_exp(j)/rmaj_exp(j))**1.5D0
        akappa1=0.8839D0*fc/(0.3477D0+0.4058D0*fc)
        alpha_neo=-akappa1+1.D0
        alpha_neo_hold=alpha_neo
 
cxc angrotp is plasma rotation
cxc angrot is impurity rotation from experiment
cxc if angrotp_exp is not suppied must insert
cx
cx        akappa2=(1.-fc)/(1.+1.1671*fc)
cxc trace impurity limit to go from impurity rotation angrot_exp to
cxc plasma rotation angrotp_exp
cx        angrotp_exp(j)=angrot_exp(j)+
cx     >     akappa2*3./2.*csda_exp(j)*zpti_exp(j)*(ti_exp(j)/te_exp(j))
cx     >     /rho(j)*q_exp(j)*rhosda_exp(j)*pgeo_local/rdrho_local
cx        if(angrot_exp(j).eq.0.) angrotp_exp(j)=0.
cx        angrotp_exp(j)=corot*angrotp_exp(j)
 
          egeo_local=1.D0
          pgeo_local=drhodr(j)
          rdrho_local=rmin_exp(j)/arho_exp/(rho(j)+epsilon)
 
        ve(j)=-(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
     > (zpni_m(j)+alpha_neo*zpti_m(j))*vstar_sign*pgeo_local
     > -rho(j)*rdrho_local*
     >  arho_exp/rmajor_exp/q_exp(j)*rmajor_exp*angrotp_exp(j)
 
        vpar(j)=rmajor_exp*angrotp_exp(j)-vstar_sign*
     > (ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*pgeo_local*
     >((alpha_neo-1.D0)*zpti_m(j))*rho(j)*rdrho_local*
     > arho_exp/rmajor_exp/q_exp(j)
 
c        vmode(j)=anfreq_m(j)/(ky_j(j)+epsilon)*
c     >    csda_m(j)*arho_exp*rhosda_m(j)
 
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
     >   ((alpha_neo-1.D0)*zpti_m(j))*rho(j)*rdrho_local*
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
 
c       vexb_m(j)=ve(j)
c       vmode_m(j)=vmode(j)
c       vstar_m(j)=(ti_m(j)/te_m(j))*csda_m(j)*arho_exp*rhosda_m(j)*
c    > (zpni_m(j)+zpti_m(j))*vstar_sign*pgeo_local
 
      enddo
c 
c compute shears from outside minus inside
c 
        j=jm
        rho(j-jptf)=rho(j-jptf)+epsilon
        rho(j+jpt)=rho(j+jpt)+epsilon
        egeo_local=1.D0
        pgeo_local=drhodr(j)
        rdrho_local=rmin_exp(j-jptf)/arho_exp/rho(j-jptf)
        rdrho_local_p1=rmin_exp(j+jpt)/arho_exp/rho(j+jpt)
c
        egamma_m(jm)=relx*egamma_m(jm)+(1.D0-relx)*(
     >  egeo_local*drhodrrrho(j)*
     >  (rho(j-jptf)+rho(j+jpt))/(q_exp(j-jptf)+q_exp(j+jpt))*
     >  (ve(j+jpt)*q_exp(j+jpt)/rho(j+jpt)/rdrho_local_p1-
     >  ve(j-jptf)*q_exp(j-jptf)/rho(j-jptf)/rdrho_local)/
     >  (rho(j+jpt)-rho(j-jptf)+epsilon)/arho_exp/csda_m(j)
     >  )
c
c       write(*,*) jm, rho(jm), egamma_m(jm), ' egamma'
c
c       gamma_mode_m(jm)=relx*gamma_mode_m(jm)+(1.-relx)*(
c    >  egeo_local*drhodrrrho(j)*
c    >  (rho(j-jptf)+rho(j+jpt))/2.*
c    >  (vmode(j+jpt)/rho(j+jpt)/rdrho_local_p1-
c    >  vmode(j-jptf)/rho(j-jptf)/rdrho_local)/
c    >  (rho(j+jpt)-rho(j-jptf)+epsilon)/arho_exp/csda_m(j)
c    >  )
 
        gamma_p_m(jm)=relx*gamma_p_m(jm)+(1.D0-relx)*(
     >   -drhodr(j)*
     >   (vpar(j+jpt)-vpar(j-jptf))
     >   /(rho(j+jpt)-rho(j-jptf)+epsilon)/arho_exp/csda_m(j)
     >  )
 
       if (jm.eq.1.and.gamma_p_m(jm).gt.10)
     >    gamma_p_m(jm)=10.D0
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
 
c  test effect of canonical density gradients
c note pzmn_sol only in zpmne_q and zpmni_q which go into
c glf2d drivers rlne_gf and rlni_gf
 
c       zpmne_q=(1.+xparam_gf(20))*zpmne*pzmn_sol(jm)-xparam_gf(19)*
c     > (alog(q_exp(jm-1))-alog(q_exp(jm)))/(rho(jm-1)-rho(jm))
c       zpmni_q=(1.+xparam_gf(20))*zpmni*pzmn_sol(jm)-xparam_gf(19)*
c     > (alog(q_exp(jm-1))-alog(q_exp(jm)))/(rho(jm-1)-rho(jm))
 
       zpmne_q= zpmne
       zpmni_q= zpmni
 
       rmaj_gf=rmaj_exp(jm)/arho_exp
       rmin_gf=rmin_exp(jm)/arho_exp
       q_gf=q_exp(jm)
       betae_gf=dmax1(cbetae*betae_m(jm), 1.D-6)
       shat_gf=shat_exp(jm)*drhodrrrho(jm)
 
       if(xalpha.lt.0.) alpha_gf=-xalpha*alpha_m(jm)
       if(xalpha.gt.0.) alpha_gf=xalpha*alpha_exp(jm)
       if(alpha_gf.gt.4.D0) alpha_gf=4.D0
 
       elong_gf=elong_exp(jm)
 
       xnu_gf=cxnu*xnu_m(jm)
       taui_gf=tim/tem
       amassgas_gf=amassgas_exp
 
       apwt_gf=1.D0
c impurity dynamics not turned on by default 
c and simple dilution included (idengrad=2, dil_gf=1-nim/nem)
c to turn on impurity dynamics need to change number of roots
c supply zimp_exp, amassimp_exp, and fractional density weights
c apwt_gf and aiwt_gf
       dil_gf=0.D0
       aiwt_gf=0.D0
c      zimp_gf=6.D0
c      amassimp_gf=12.D0
       zimp_gf=zimp_exp
       amassimp_gf=amassimp_exp
       rlnimp_gf=1.D0
       zpmnimp=1.D0
       if (idengrad.eq.2) dil_gf=1.D0-nim/nem
       if (idengrad.eq.3) then
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
c       write(*,66) rho(jm),nem,nim,rlnimp_gf,zpmnimp,amassimp_gf,
c     >             xnimp,zimp_gf,dil_gf
c 66    format(2x,0p1f4.2,0p9f10.5)
 
         rlte_gf=zpmte*drhodr(jm)
         rlti_gf=zpmti*drhodr(jm)
         rlne_gf=zpmne_q*drhodr(jm)
         rlni_gf=zpmni_q*drhodr(jm)
         rlnimp_gf=zpmnimp*drhodr(jm)
c        write(*,200) jm, rho(jm), rlti_gf, rlte_gf, rlne_gf, rlni_gf
 
        gamma_star_gf=vstarp_m(jm)
        gamma_e_gf=egamma_m(jm)
        gamma_p_gf=gamma_p_m(jm)
        if(itport_pt(4).eq.-1) then
          gamma_e_gf=egamma_exp(jm)
          gamma_p_gf=gamma_p_exp(jm)
        endif
c..jek 8/15/00
        if (irotstab.eq.0) then
          gamma_e_gf=egamma_exp(jm)
          gamma_p_gf=gamma_p_exp(jm)
          egamma_m(jm)=egamma_exp(jm)
        endif
c       gamma_mode_gf=gamma_mode_m(jm)
        gamma_mode_gf=0.0D0
 
        if(i_delay.ne.0) then
c   i_delay should be negative on any intermediate step
 
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
 
c...diagnostic arrays (not presently used)
c
c      ky_j(jm)=xkyf_gf
c      gamma_j(jm,1)=gamma_gf(1)
c      gamma_j(jm,2)=gamma_gf(2)
c      gamma_j(jm,3)=gamma_gf(3)
c      gamma_j(jm,4)=gamma_gf(4)
c
c      freq_j(jm,1)=freq_gf(1)
c      freq_j(jm,2)=freq_gf(2)
c      freq_j(jm,3)=freq_gf(3)
c      freq_j(jm,4)=freq_gf(4)
c
c      phi_norm_j(jm,1)=phi_norm_gf(1)
c      phi_norm_j(jm,2)=phi_norm_gf(2)
c      phi_norm_j(jm,3)=phi_norm_gf(3)
c      phi_norm_j(jm,4)=phi_norm_gf(4)
c
c      do ik=1,ikymax_gf
c       gamma_k_j(ik,jm)= gamma_k_gf(1,ik)
c       freq_k_j(ik,jm) = freq_k_gf(1,ik)
c       chie_k_j(ik,jm) = chie_k_gf(ik)
c       chii_k_j(ik,jm) = chii_k_gf(ik)
c      enddo
 
       anrate_m(jm)=gamma_gf(1)
       anrate2_m(jm)=gamma_gf(2)
c      dnrate_m(jm)=0.
c      dtnrate_m(jm)=0.
       anfreq_m(jm)=freq_gf(1)
       anfreq2_m(j)=freq_gf(2)
c      dnfreq_m(jm)=0.
c       xkymax_m(jm)=xky_gf(1)
c       xkymax2_m(jm)=xky_gf(2)

       gfac=geofac(jm)
c
c exch_m in MW/m**3
c   kev/sec per MW
ck       kevdsecpmw=1.6022e-19*1.0e3*1.e-6
 
ck      exch_m(jm)=1.e19*
ck     > kevdsecpmw*nem*tem*csda_m(jm)*rhosda_m(jm)**2*exch_gf*cmodel
      exch_m(jm)=1.D19*1.6022D-19*1.0D3*1.D-6*
     > nem*tem*csda_m(jm)*rhosda_m(jm)**2*exch_gf*cmodel
 
      exchm=exch_m(jm)
 
c exch_m is directly related to the flow
c for a single mode branch exch_gf=-(-freq_gf(1)/xkyf_gf)*diff_gf*rln_gf.
c we can  not expect to find exch_m without knowing flow_exp as input.
c and solving self consistant flow eq. flown=flow_exp for density
c density profile.
c
c however, knowing freq_gf(1) from the gf model we can compute exch_exp
c from flow_exp using
c       flowm=kevdsecpmw*1.*nem*1.e19/arho_exp*gradrhosq_exp(jm)*
c     >       sfactor(jm)*(difftem*zpmte+difftim*zpmti+diffnem*zpmne)
c we have:
 
c       diffgb_local=flow_exp(jm)/
c     > (kevdsecpmw*1.*nem*1.e19/arho_exp*gradrhosq_exp(jm)*sfactor(jm)*
c     > zpmne_q)/cgyrobohm_m(jm)
 
c       exchgb_local=-(-freq_gf(1)/xkyf_gf)*diffgb_local*rlni_gf
 
c       exch_exp(jm)=1.e19*
c     > kevdsecpmw*nem*tem*csda_m(jm)*rhosda_m(jm)**2*exchgb_local
 
c       exch_exp(jm)=flow_exp(jm)*tem*(-1.)*(-freq_gf(1)/xkyf_gf)*
c     > sqrt(elong_exp(jm))*arho_exp/gradrhosq_exp(jm)/sfactor(jm)
c     >/arho_exp(jm)**2
 
 
c   note electron(ion) wave freq > 0(<0) cool(heat) electrons
c (-1) denotes electron to ion
 
c to emphasize, we can not know exch_exp better than we know flow_exp
 
c chietem, chiitim, diffen in m**2/sec
c ITER definition of "chi" assumes will be proceeded with gradrhosq factor
 
      chietem=cmodel*gfac*chie_gf*cgyrobohm_m(jm)
      chiitim=cmodel*gfac*chii_gf*cgyrobohm_m(jm)
      diffnem=cmodel*gfac*diff_gf*cgyrobohm_m(jm)
 
      etaphim=cmodel*gfac*eta_phi_gf*cgyrobohm_m(jm)
      etaparm=cmodel*gfac*eta_par_gf*cgyrobohm_m(jm)
      etaperm=cmodel*gfac*eta_per_gf*cgyrobohm_m(jm)
 
cxx endif
c
      if ( itport_pt(1) .eq. 0 ) diffnem = 0.D0
      if ( itport_pt(2) .eq. 0 ) chietem = 0.D0
      if ( itport_pt(3) .eq. 0 ) chiitim = 0.D0
      if ( (itport_pt(4) .eq. 0) .and. (itport_pt(5) .eq. 0) ) then
         etaphim=0.D0
         etaparm=0.D0
         etaperm=0.D0
      endif
c
      diff_m(jm)=diffnem
      chie_m(jm)=chietem
      chii_m(jm)=chiitim
 
      etaphi_m(jm)=etaphim
      etapar_m(jm)=etaparm
      etaper_m(jm)=etaperm
c
c     write(6,*) jm,rho(jm),zpmte, zpmti, zpmne, zpmni
 
        enddo
 
c      if(jmm.eq.0)  write(6,*) 'jmm=', jmm,'going out'
c
cdmc        write(6,7801) diffnem,chietem,chiitim,etaphim,etaparm,etaperm,
cdmc     >     exchm
cdmc 7801   format(//' OUTPUTS...'/
cdmc     >'  diffnem = ',1pe17.10,' chietem = ',1pe17.10/
cdmc     >'  chiitim = ',1pe17.10,' etaphim = ',1pe17.10/
cdmc     >'  etaparm = ',1pe17.10,' etaperm = ',1pe17.10/
cdmc     >'  exchm = ',1pe17.10)
c
cdmc        call echo('diff_m output',diff_m,jmaxm)
cdmc        call echo('chie_m output',chie_m,jmaxm)
cdmc        call echo('chii_m output',chii_m,jmaxm)
cdmc        call echo('etaphi_m output',etaphi_m,jmaxm)
cdmc        call echo('etapar_m output',etapar_m,jmaxm)
cdmc        call echo('etaper_m output',etaper_m,jmaxm)
cdmc        call echo('exch_m output',exch_m,jmaxm)
c
cdmc        call echo('egamma_m output',egamma_m,jmaxm)
cdmc        call echo('egamma_d output',egamma_d,jmaxm)
cdmc        call echo('gamma_p_m output',gamma_p_m,jmaxm)
cdmc        call echo('anrate_m output',anrate_m,jmaxm)
cdmc        call echo('anrate2_m output',anrate2_m,jmaxm)
cdmc        call echo('anfreq_m output',anfreq_m,jmaxm)
cdmc        call echo('anfreq2_m output',anfreq2_m,jmaxm)
c
       return
       end
c
cdmc      subroutine echo(lbl,arr,jmaxm)
cdmc      character*(*) lbl
cdmc      real*8 arr(0:jmaxm)
c
cdmc      write(6,1001) lbl
cdmc 1001 format(/5x,a)
cdmc      write(6,1002) (arr(j),j=0,jmaxm)
cdmc 1002 format(4(1x,1pe17.10))
c
cdmc      return
cdmc      end
