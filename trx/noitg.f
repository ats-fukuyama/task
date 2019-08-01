c
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

c      include 'glf.m'
 
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
     >  zpte_m(0:100),zpti_m(0:100),
     >  zpne_m(0:100),zpni_m(0:100),
     >  drhodr(0:100),drhodrrrho(0:100),geofac(0:100),
     >  rhosda_m(0:100),csda_m(0:100),cgyrobohm_m(0:100),
     >  betae_m(0:100),xnu_m(0:100),
     >  alpha_m(0:100),vstarp_m(0:100)
 
c working arrays and variables
 
      double precision ve(0:100),vpar(0:100)
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
 
      return
      end

      subroutine disp9t(neq,ZZ)
      implicit none
      INTEGER neq
      COMPLEX*16 ZZ(*)
      return
      end
      SUBROUTINE DIFF(RP,IR,TAUI,FT,U,CHI,CHE,D,CHQ,DHQ)
C
      IMPLICIT NONE
      COMPLEX*16 RP(10)
      INTEGER IR
      REAL*8 TAUI,FT
      REAL*8 CHI(5),CHE(5),D(5)
      REAL*8 CHQ(5),DHQ(5)
      REAL*8 U(5,100)
      return
      end

      subroutine IFSPPPL( znine, zncne, znbne, zrlt, zrln,
     &                    zq, zshat, zeps, zkappa, omegaexb,
     &                    ne19, tekev, tikev, rmajor, btesla,
     &                    switches, grhoi, gvti, gnu, gtau,
     &                    chii, chie,
     &                    zchiicyc, zchii1, zchii2, zchie1, zchie2,
     &                    zrlt1, zrlt2, ierr )
c
      IMPLICIT NONE
c
c
c
      INTEGER switches(32), ierr
      REAL znine, zncne, znbne, zrlt, zrln, zq, zshat, zeps,
     &       ne19, tekev, tikev, rmajor, rhoi, vthermi, omegaexb,
     &       ztau, znu, zzeffstr, zrlnstr, ztaub, zkappa,
     &       zchii, zchii1, zchii2, zscriptg1, zscriptg2, zx,
     &       zscriptz, zf, zg, zh, zscriptd, zscripte,
     &       zrlt1, zrlt2, zchie, zchie1, zchie2, zscriptj,
     &       btesla, chie, chii, omegaci, zcorrect, zkapfact,
     &       grhoi, gvti, gtau, gnu, zscripty, gamma,
     &       zc1, zc2, zc3, zscriptg, zchiicyc, zscripts
      CHARACTER*35 calcflag
      return
      end

