MODULE trglf23

  PRIVATE
  PUBLIC tr_glf23

CONTAINS

! ******************************************************************
!
!  The interface between TASK/TR(trcoeftb) and GLF23 model(callglf2d)
!
! ******************************************************************

  SUBROUTINE tr_glf23

! ******************************************************************
!  In case of jshoot=0 and sometimes jmm=0, zeroth arguments of
!   te_m, ti_m, ne_m and ni_m  require finite values
!   typically the same ones as first arguments, and those of rho and 
!   rmin_exp require zero avoiding numerical error due to the absence
!   of the value. 
!  However, the value calculated by using zeroth value of rmin_exp 
!   are not used in this case.
! ******************************************************************

    USE trcomm, ONLY : &
         rkind,ikind,pi,rkev,BB,RR,ra,rkap,mdltr_tb,   &
         idnsa,mdluf,pz_mion,pa_mion,pz_mimp,pa_mimp,  &
         nrmax,nsamax,neqmax,pa,pz,pz0,rhog,rhom,phia, &
         qp,rn,rt,rmu0,ns_nsa,dtr_tb,vtr_tb,           &
         abb1rho,ar1rho,ar2rho,rkprho,rmjrho,rmnrho,   &
         cdtrn,cdtrt,vtor,vpar,vprp,wrot,              &
         rt_e,    &! the electron temperature
         rt_i,    &! the ion temperature
         rn_e,    &! the electron density
         rn_i,    &! the ion density
         rt_ecl,  &! the scale length of electron temperature
         rt_icl,  &! the scale length of ion temperature
         rn_ecl,  &! the scale length of electron density
         rn_icl,  &! the scale length of ion density
         rp_totd, &! the derivative of total pressure w.r.t. 'r'
         ai_ave,  &! mean atomic mass of thermal ions [AMU]
         mshear,  &! the magnetic shear
         alpha,   &! MHD alpha
         z_eff,   &! effective charge
         wexbp

    IMPLICIT NONE
    ! Inputs for callglf2d
    INTEGER(ikind),DIMENSION(5) :: itport_pt
    REAL(rkind),DIMENSION(0:nrmax) :: &
         te_m, ti_m, ne_m, ni_m, ns_m, angrotp_exp, egamma_exp, &
         gamma_p_exp, vphi_m, vpar_m, vper_m, zeff_exp, rho, &
         gradrho_exp, gradrhosq_exp, rmin_exp, rmaj_exp, q_exp, &
         shat_exp, alpha_exp, elong_exp

    ! Outputs for callglf2d
    REAL(rkind),DIMENSION(0:nrmax) :: &
         diff_m, chie_m, chii_m, etaphi_m, etapar_m, etaper_m, exch_m, &
         egamma_m, rgamma_p_m, anrate_m, anrate2_m, anfreq_m, anfreq2_m
    REAL(rkind),DIMENSION(0:nrmax,10) :: egamma_d
    
    ! Auxiliaries
    INTEGER(ikind) :: &
         i_delay, idengrad, iglf, i_grad, irotstab, jmaxm, jmm, &
         jshoot, leigen, bt_flag, nr, nroot, ns, ns1
    REAL(rkind) :: &
         alpha_e, amassgas_exp, amassimp_exp, arho_exp, bt_exp, chietem, &
         chiitim, diffnem, etaparm, etaperm, etaphim, exchm, &
         rmajor_exp, x_alpha, zimp_exp, zpne_in, &
         zpni_in, zpte_in, zpti_in
    
    ! Internal variables and others
    REAL(rkind),DIMENSION(0:nrmax)::&
         diff_jm,chie_jm,chii_jm,diff_jmg,chie_jmg,chii_jmg, &
         vtorm,vprpm,vparm
    INTEGER(ikind) :: mdleqn,mdleqt,mdleoi,nsa,nbase

    ! TASK/lib
    REAL(rkind) :: FCTR
    
    ! initilization
    dtr_tb(1:neqmax,1:neqmax,1:nrmax)=0.D0
    vtr_tb(1:neqmax,1:neqmax,1:nrmax)=0.D0

    diff_jmg(0:nrmax) = 0.d0
    chie_jmg(0:nrmax) = 0.d0
    chii_jmg(0:nrmax) = 0.d0
    diff_jm(0:nrmax)  = 0.d0
    chie_jm(0:nrmax)  = 0.d0
    chii_jm(0:nrmax)  = 0.d0

    jmaxm = nrmax     ! maximum num. of grid points

    leigen=1    ! 1 for tomsqz, 0 for cgg solver
    IF(mdluf /= 0 .AND. nsamax > 2) THEN
       nroot = 12 ! num. of equations, 8 for default, 12 for impurity dynamics
       idengrad = 3 ! actual dilution
    ELSE
       nroot = 8
       idengrad = 2 ! compute simple dilution
    ENDIF
    iglf   = 1  ! 0 for original model, 1 for retuned model

    jshoot = 0  ! 0 for time-dep code, 1 for shooting code
    ! In case of jshoot=0, maximum argument of array is important;
    !  values of zeroth argument are not important, but avoiding Inf error
    !  due to logarithm calculation some finite values need to be stored.


    ! -----------------------------------------------------------
    ! itport_pt(1) : desnity trasnport
    ! itport_pt(2) : electron trasnport
    ! itport_pt(3) : ion trasnport
    ! itport_pt(4) : v_phi trasnport   ( = -1 : use egamma_exp)
    ! itport_pt(5) : v_theta trasnport ( = -1 : use gamma_p_exp)
    !
    !               = 1 : calculate,  = 0 : not calculate
    !
    ! --- switches and variables from TASK/TR ---
    ! ***  These substitution is a interim way. ***
    mdleqn= 0 ! density transport
    mdleqt= 1 ! energy transport
    mdleoi= 0 ! electron or ion

    itport_pt(1) = MDLEQN
    IF(MDLEQT /= 0) THEN
       IF(nsamax == 1) THEN
          IF(MDLEOI == 1) THEN
             itport_pt(2) = MDLEQT
             itport_pt(3) = 0
          ELSEIF(MDLEOI == 2) THEN
             itport_pt(2) = 0
             itport_pt(3) = MDLEQT
          ELSE ! dafalut
             itport_pt(2) = MDLEQT
             itport_pt(3) = MDLEQT
          ENDIF
       ELSE
          itport_pt(2) = MDLEQT
          itport_pt(3) = MDLEQT
       ENDIF
    ELSE
       itport_pt(2) = MDLEQT
       itport_pt(3) = MDLEQT
    ENDIF
    itport_pt(4) = 0  ! v_phi transport
    itport_pt(5) = 0  ! v_theta transport
    
    irotstab=1    ! 1 uses internally computed ExB shear, 0 for prescribed
        
    ! variables from experimental data ---------------------
    IF(     itport_pt(4) ==  0 .AND. itport_pt(5) ==  0)THEN
       ! experimental toroidal angular velocity
       angrotp_exp(0:nrmax) = wrot(0:nrmax) ! taking over 'VROT' in exp. data

    ELSE IF(itport_pt(4) == -1 .AND. itport_pt(5) ==  0)THEN
       ! experimental ExB shearing rate
       egamma_exp(0:nrmax) = 0.d0

    ELSE IF(itport_pt(4) ==  0 .AND. itport_pt(5) == -1)THEN
       ! experimental  parallel velocity shearing rate
       gamma_p_exp(0:nrmax) = 0.d0

    ELSE IF(itport_pt(4) ==  1 .AND. itport_pt(5) ==  0)THEN
       vphi_m(0:nrmax) = vtor(0:nrmax)

    ELSE IF(itport_pt(4) ==  1 .AND. itport_pt(5) ==  1)THEN
       vpar_m(0:nrmax) = vpar(0:nrmax)

    ELSE IF(itport_pt(4) ==  1 .AND. itport_pt(5) ==  1)THEN
       vper_m(0:nrmax) = vprp(0:nrmax)
    END IF
    ! ------------------------------------------------------
    
    alpha_e = 1.D0     ! ExB shear stabilization (0=off,>0=on)
    x_alpha = 1.D0     ! alpha stabilization (0=off,>0=on)
    i_delay = 0.d0     ! default (usually 0 is recommended for stabilityXX)
    
    bt_exp     = BB    ! vaccume axis toroidal field [T]
    bt_flag    = 1     ! >0 for Beff, Bt otherwise
    rmajor_exp = RR    ! geometrical major radius of magnetic axis [m]

    zimp_exp     = pz_mimp  ! Zimp; finite data is necessary
    amassimp_exp = pa_mimp  ! Aimp; finite data is necessary
    amassgas_exp = pa_mion  ! atomic num. of working hydrogen gas


    ! --- geometiric factor
    IF(phia.EQ.0.D0) THEN
       arho_exp = 1.d0 ! rho at last closed flux surface [m]
    ELSE ! in the case experimental data has been read.
       arho_exp=SQRT(phia/(pi*BB))
    END IF
    rho(0:nrmax) = rhog(0:nrmax)
    
    rmin_exp(0:nrmax) = rmnrho(0:nrmax)    ! local minor radius [m]
    rmaj_exp(0:nrmax) = rmjrho(0:nrmax)    ! local major radius [m]
    
    gradrho_exp(0:nrmax)   = ar1rho(0:nrmax)
    gradrhosq_exp(0:nrmax) = ar2rho(0:nrmax)
    
    ! local elongation
    elong_exp(0:nrmax) = rkprho(0:nrmax)

    ! Values related to physical gradients
    !                          are evaluated by forward difference.
    ! magnetic shear
    shat_exp(0)         = FCTR(rhom(1),rhom(2),mshear(1),mshear(2))
    shat_exp(1:nrmax-1) = 0.5d0*(mshear(1:nrmax-1)+mshear(2:nrmax))
    shat_exp(nrmax)     = mshear(nrmax)
    ! MHD alpha: - q^2 R (d beta/d r)
    alpha_exp(0)         = FCTR(rhom(1),rhom(2),alpha(1),alpha(2))
    alpha_exp(1:nrmax-1) = 0.5d0*(alpha(1:nrmax-1)+alpha(2:nrmax))
    alpha_exp(nrmax)     = alpha(nrmax)
    
    te_m(0:nrmax) = rt_e(0:nrmax)       ! Te [keV]
    ti_m(0:nrmax) = rt_i(0:nrmax)       ! Ti [keV]
    ne_m(0:nrmax) = rn_e(0:nrmax)*10.d0 ! Ne [10^19 /m^3]
    ni_m(0:nrmax) = rn_i(0:nrmax)*10.d0 ! Ni [10^19 /m^3]

    ns_m(0:nrmax) = 0.d0
    DO nsa = 1, nsamax
       IF(idnsa(nsa) == 2)THEN
          ! Fast ion density [10^19 /m^3]
          ns_m(0:nrmax) = ns_m(0:nrmax) + rn(nsa,0:nrmax)*10.d0          
       END IF
    END DO
    
    q_exp(0:nrmax)    = qp(0:nrmax)     ! safety factor
    zeff_exp(0:nrmax) = z_eff(0:nrmax)  ! effective charge
    
    ! --- NR LOOP -----------------------------------------
    ! [nr=1, nrmax] corresponds in TASK/TR to [jm=1, jmaxm] in callglf2d.f.
    DO nr = 1, nrmax-1 ! on grid points
       jmm   = nr      ! jmm=0 does full grid from jm=1,jmaxm-1
       
       IF(mdltr_tb.EQ.60) THEN
          !  +++ Normal type +++
          
          ! compute gradients (1=input gradients)
          i_grad=0
          ! the derivatives below variables are with respect to 'rho'
          ! ( see callglf2d.f l.592~l.603 )
          ! 1/Lte (necessary if i_grad and jmm != 0)
          zpte_in = 0.d0
          ! 1/Lne (necessary if i_grad and jmm != 0)
          zpti_in = 0.d0
          ! 1/Lti (necessary if i_grad and jmm != 0)
          zpne_in = 0.d0
          ! 1/Lni (necessary if i_grad and jmm != 0)
          zpni_in = 0.d0
          !            write(*,*) zpte_in,zpti_in,zpne_in,zpni_in
       ENDIF


       call callglf2d( &
!          << INPUTS >>
            leigen,    &! eigenvalue solver
!                       ! 0 for cgg (default), 1 for tomsqz, 2 for zgeev
            nroot,     &! no. roots, 8 for default, 12 for impurity dynamics
            iglf,      &! 0 for original GLF23, 1 for retuned version
            jshoot,    &! jshoot=0 time-dep code; jshoot=1 shooting code
            jmm,       &! grid number; jmm=0 does full grid jm=1 to jmaxm-1
            jmaxm,     &! profile grids 0 to jmaxm
            itport_pt, &! transport flags (array[1:5])
            irotstab,  &! 0 to use egamma_exp; 1 use egamma_m
            te_m,      &! 0:jmaxm te [keV]        : itport_pt(2)=1 transport
            ti_m,      &! 0:jmaxm ti [keV]        : itport_pt(3)=1 transport
            ne_m,      &! 0:jmaxm ne [10^19 /m^3]
            ni_m,      &! 0:jmaxm ni [10^19 /m^3] : itport_pt(1)=1 transport
            ns_m,      &! 0:jmaxm ns [10^19 /m^3]
            i_grad,    &! default 0, 
!                       ! for D-V method use i_grad=1 to input gradients
            idengrad,  &! default 2, for simple dilution
!                       ! ** below zpXX's are for i_grad=1
            zpte_in,   &! externally provided log gradient te w.r.t rho
            zpti_in,   &! externally provided log gradient ti w.r.t rho
            zpne_in,   &! externally provided log gradient ne w.r.t rho
            zpni_in,   &! externally provided log gradient ni w.r.t rho
            angrotp_exp,&!0:jmaxm exp plasma toroidal angular velocity 1/sec
                         ! if itport_pt(4)=0 itport_pt(5)=0 
            egamma_exp, &!0:jmaxm exp exb shear rate in units of csda_exp
                         ! if itport_pt(4)=-1 itport_pt(5)=0 
            gamma_p_exp,&!0:jmaxm exp par.vel.shear rate in units of csda_exp
                         ! if itport_pt(4)=0  itport_pt(5)=-1
            vphi_m,    &! 0:jmaxm toroidal velocity m/sec
!                       !  if itport_pt(4)=1 itport_pt(5)=0 otherwise output
            vpar_m,    &! 0:jmaxm parallel velocity m/sec
!                       !  if itport_pt(4)=1 itport_pt(5)=1 otherwise output
            vper_m,    &! 0:jmaxm perp. velocity m/sec
!                       !  if itport_pt(4)=1 itport_pt(5)=1 otherwise output 
            zeff_exp,  &! effective charge
            bt_exp,    &! vaccume axis toroidal field in Tesla
            bt_flag,   &! switch for effective toroidal field use in rhosda
            rho,       &! 0:jmaxm 0<rho<1 normalized toloidal flux
!                       !  (rho = rho/rho(a))
            arho_exp,  &! rho(a); toroidal flux at last closed flux surface
!                       ! (LCFS)   toroidal flux = B0*rho_phys^2/2 [m]
!                       !          B0 = bt_exp,  arho_exp = rho_phys_LCFS
            gradrho_exp,&  ! 0:jmaxm dimensionless <|grad rho_phys|^2>
            gradrhosq_exp,&! 0:jmaxm dimensionless <|grad rho_phys|>
            rmin_exp,  &! 0:jmaxm minor radius in meters
            rmaj_exp,  &! 0:jmaxm major radius in meters
            rmajor_exp,&! axis major radius
            zimp_exp,  &  ! effective Z of impurity 
            amassimp_exp,&! effective A of impurity
            q_exp,     &! 0:jmaxm safety factor   
            shat_exp,  &! 0:jmaxm magnetic shear, d (ln q_exp)/ d (ln rho)
            alpha_exp, &! 0:jmaxm MHD alpha from experiment   
            elong_exp, &! 0:jmaxm elongation 
            amassgas_exp,&! atomic number working hydrogen gas
            alpha_e,   &! 1 full (0 no) no ExB shear stab
            x_alpha,   &! 1 full (0 no) alpha stabilization  with alpha_exp
!                       !-1 full (0 no) self consistent alpha_m stab.
            i_delay,   &! i_delay time delay for ExB shear 
!                       !  should be non-zero only once per step 
!                          and is less than or equal 10.
!
!        << OUTPUTS >>
            diffnem,   &! ion plasma diffusivity             [m^2/sec]
            chietem,   &! electron ENERGY diffuivity         [m^2/sec]
            chiitim,   &! ion      ENERGY diffusivity        [m^2/sec]
            etaphim,   &! toroidal velocity diffusivity      [m^2/sec]
            etaparm,   &! parallel velocity diffusivity      [m^2/sec]
            etaperm,   &! perpendicular velocity diffusivity [m^2_sec]
            exchm,     &! turbulent electron to ion ENERGY exchange [MW/m^3]
            diff_m,    &! **
            chie_m,    &! ** 
            chii_m,    &! ** 
            etaphi_m,  &! **  These arguments are arrays of above arguments.
            etapar_m,  &! **
            etaper_m,  &! **
            exch_m,    &! **
            egamma_m,  &!0:jmaxm exb shear rate in units of local csda_m
            egamma_d,  &!0:jmaxm exb shear rate delayed by i_delay steps
            rgamma_p_m,&!0:jmaxm par.vel.shear rate in units of local csda_m
            anrate_m,  &!0:jmaxm leading mode rate in unints of local csda_m
            anrate2_m, &!0:jmaxm 2nd mode rate in units of local csda_m
            anfreq_m,  &!0:jmaxm leading mode frequency
            anfreq2_m ) !0:jmaxm 2nd mode frequency

       diff_jmg(nr) = MAX(diffnem,0.d0) 
       chie_jmg(nr) = MAX(chietem,0.d0)
       chii_jmg(nr) = MAX(chiitim,0.d0)

    ENDDO

    diff_jm(1:nrmax-1) = 0.5d0*(diff_jmg(0:nrmax-2)+diff_jmg(1:nrmax-1))
    chie_jm(1:nrmax-1) = 0.5d0*(chie_jmg(0:nrmax-2)+chie_jmg(1:nrmax-1))
    chii_jm(1:nrmax-1) = 0.5d0*(chii_jmg(0:nrmax-2)+chii_jmg(1:nrmax-1))
    diff_jm(nrmax) = diff_jmg(nrmax-1)
    chie_jm(nrmax) = chie_jmg(nrmax-1)
    chii_jm(nrmax) = chii_jmg(nrmax-1)

    vtorm(1:nrmax-1) = 0.5d0*(vphi_m(0:nrmax-2)+vphi_m(1:nrmax-1))
    vparm(1:nrmax-1) = 0.5d0*(vpar_m(0:nrmax-2)+vpar_m(1:nrmax-1))
    vprpm(1:nrmax-1) = 0.5d0*(vper_m(0:nrmax-2)+vper_m(1:nrmax-1))
    vtorm(nrmax) = vphi_m(nrmax-1)
    vparm(nrmax) = vpar_m(nrmax-1)
    vprpm(nrmax) = vper_m(nrmax-1)

    DO nsa=1,nsamax
       ns = ns_nsa(nsa)
       nbase = 3*(nsa-1)+1
       IF(pz0(ns) < 0.D0) THEN ! for electron
          dtr_tb(nbase+1,nbase+1,1:nrmax)=cdtrn*MAX(diff_jm(1:nrmax),0.d0)
          dtr_tb(nbase+3,nbase+3,1:nrmax)=cdtrt*MAX(chie_jm(1:nrmax),0.d0)
!          write(*,*) chietem,diffnem
       ELSE 
          IF(pz(ns) /= 0.d0) THEN ! for ion
             dtr_tb(nbase+1,nbase+1,1:nrmax)= &
                                          cdtrn*MAX(diff_jm(1:nrmax),0.d0)
             dtr_tb(nbase+3,nbase+3,1:nrmax)= &
                                          cdtrt*MAX(chii_jm(1:nrmax),0.d0)
          END IF
       END IF
    END DO

    RETURN
  END SUBROUTINE tr_glf23
END MODULE trglf23
