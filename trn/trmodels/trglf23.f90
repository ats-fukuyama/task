MODULE trglf23

  PRIVATE
  PUBLIC tr_glf23

CONTAINS

! ******************************************************************
!
!  The interface between TASK/TR(trcalv) and GLF23 model(callglf2d)
!
! ******************************************************************

  SUBROUTINE tr_glf23

! ******************************************************************
!  In case of jshoot=0 and sometimes jmm=0, zeroth arguments of
!   te_m, ti_m, ne_m and ni_m require finite values
!   typically the same ones as first arguments, and those of rho and 
!   rmin_exp require zero avoiding numerical error due to the absence
!   of the value. 
!  However, the value calculated by using zeroth value of rmin_exp 
!   are not used in this case.

! ** Input and output profile are on grid points (nr=1, nrmax-1) **

! ******************************************************************

    USE trcomm, ONLY : &
         rkind,ikind,BB, abb1rho, mdltr_tb, rkev, &
         nrmax, nsamax, pa, pi, pz, pz0, qp, RR,ra,rhog,rkap, &
         rn,rt,rmu0,ns_nsa,dtr_tb,vtr_tb,  &
         ar1rho,ar2rho,rkprho,rmjrho,rmnrho,    &
         vtor
    USE trlib, ONLY: mesh_convert_mtog
    
    USE trcalv, ONLY: &
         rt_e,    &! the electron temperature
         rt_i,    &! the ion temperature
         rn_e,    &! the electron density
         rn_i,    &! the ion density
         rt_ecl,  &! the scale length of electron temperature
         rt_icl,  &! the scale length of ion temperature
         rn_ecl,  &! the scale length of electron density
         rn_icl,  &! the scale length of ion density
         rp_totd, &! the derivative of total pressure w.r.t. 'r'
         mshear,  &! the magnetic shear
         z_eff     ! effective charge

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
         Wexb_exp,Wrot,Vpar,Vprp,Vpar_shear,rp_totdg, &
         diff_jm,chie_jm,chii_jm,diff_jmg,chie_jmg,chii_jmg
    REAL(rkind) :: phia
    INTEGER(ikind) :: nsmax,mdluf,mdleqn,mdleqt,mdleoi,nsa,nbase
    
    ! initilization
    dtr_tb(1:3*nsamax,1:3*nsamax,1:nrmax)=0.D0
    vtr_tb(1:3*nsamax,1:3*nsamax,1:nrmax)=0.D0


    phia  = 0.d0
    nsmax = nsamax
    jmaxm = nrmax     ! maximum num. of grid points

    ! --- switches and variables from TASK/TR ---
    mdluf = 0
    mdleqn= 0 ! density transport
    mdleqt= 1 ! energy transport
    mdleoi= 0 ! electron or ion

    Vtor = 0.d0
    Vpar = 0.d0
    Vprp = 0.d0

    ! from experimental data
    Wrot       = 0.d0
    Wexb_exp   = 0.d0
    Vpar_shear = 0.d0
    !---------------------------------------------


    leigen=1    ! 1 for tomsqz, 0 for cgg solver
    IF(MDLUF.NE.0.AND.NSMAX.GT.2) THEN
       nroot=12 ! num. of equations, 8 for default, 12 for impurity dynamics
    ELSE
       nroot=8
    ENDIF
    iglf   = 1  ! 0 for original model, 1 for retuned model
    jshoot = 0  ! 0 for time-dep code, 1 for shooting code

! *** In case of jshoot=0, maximum argument of array is important;
! ***  values of zeroth argument are not important, but avoiding Inf error 
! ***  due to logarithm calculation some finite values need to be stored.
! *** nr=1, nrmax corresponds in TASK/TR to jm=1, jmaxm in callglf2d.f.


    itport_pt(1) = MDLEQN          ! density transport
    IF(MDLEQT /= 0) THEN
       IF(nsamax == 1) THEN
          IF(MDLEOI == 1) THEN
             itport_pt(2) = MDLEQT ! electron transport
             itport_pt(3) = 0      ! ion transport
          ELSEIF(MDLEOI == 2) THEN
             itport_pt(2) = 0
             itport_pt(3) = MDLEQT
          ELSE
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

    IF(MDLUF.NE.0.AND.NSMAX.GT.2) THEN
       idengrad=3 ! actual dilution
    ELSE
       idengrad=2 ! compute simple dilution
    ENDIF
    

    bt_exp     = BB   ! vaccume axis toroidal field [T]
    bt_flag    = 1    ! >0 for Beff, Bt otherwise
    rmajor_exp = RR   ! geometrical major radius of magnetic axis [m]
    
    q_exp(0:nrmax)    = qp(0:nrmax)     ! safety factor
    shat_exp(0:nrmax) = mshear(0:nrmax) ! magnetic shear
    
    
    ! --- geometiric factor
    ! ------ for the time being ---
    ! ------ arho_exp(& phia) should be calculated trprof/trsetup ---
    phia = 0.d0
    IF(phia.EQ.0.D0) THEN
       arho_exp=SQRT(rkap)*ra ! rho at last closed flux surface [m]
    ELSE
       arho_exp=SQRT(phia/(pi*BB))
    END IF
    rho(0:nrmax) = rhog(0:nrmax)
    
    rmin_exp(0:nrmax) = rmnrho(0:nrmax) ! local minor radius [m]
    rmaj_exp(0:nrmax) = rmjrho(0:nrmax) ! local major radius [m]
    
    gradrho_exp(0:nrmax)   = ar1rho(0:nrmax)
    gradrhosq_exp(0:nrmax) = ar2rho(0:nrmax)
    
    elong_exp(0:nrmax) = rkprho(0:nrmax)
    
    
    te_m(0:nrmax) = rt_e(0:nrmax)       ! Te [keV]
    ti_m(0:nrmax) = rt_i(0:nrmax)       ! Ti [keV]
    ne_m(0:nrmax) = rn_e(0:nrmax)*10.d0 ! Ne [10^19 /m^3]
    ni_m(0:nrmax) = rn_i(0:nrmax)*10.d0 ! Ni [10^19 /m^3]
    ns_m(0:nrmax) = 0.d0*10.d0          ! Fast ion density [10^19 /m^3]
    
    
    ! variables form experimental data ---------------------
    zeff_exp(0:nrmax) = z_eff(0:nrmax)
    
    IF(     itport_pt(4) ==  0 .AND. itport_pt(5) ==  0)THEN
       angrotp_exp(0:nrmax) = Wrot(0:nrmax)
       
    ELSE IF(itport_pt(4) == -1 .AND. itport_pt(5) ==  0)THEN
       egamma_exp(0:nrmax) = Wexb_exp(0:nrmax)
       
    ELSE IF(itport_pt(4) ==  0 .AND. itport_pt(5) == -1)THEN
       gamma_p_exp(0:nrmax) = Vpar_shear(0:nrmax)
       
       
    ELSE IF(itport_pt(4) ==  1 .AND. itport_pt(5) ==  0)THEN
       vphi_m(0:nrmax) = Vtor(0:nrmax)
       
    ELSE IF(itport_pt(4) ==  1 .AND. itport_pt(5) ==  1)THEN
       vpar_m(0:nrmax) = Vpar(0:nrmax)
       
    ELSE IF(itport_pt(4) ==  1 .AND. itport_pt(5) ==  1)THEN
       vper_m(0:nrmax) = Vprp(0:nrmax)
    END IF
    ! ------------------------------------------------------
    
    alpha_exp(0:nrmax) = -2.d0*rmu0*qp(0:nrmax)**2*rmjrho(0:nrmax)  &
                         /abb1rho(0:nrmax)**2*rp_totd(0:nrmax)

    alpha_e = 1.D0        ! ExB shear stabilization (0=off,>0=on)
    x_alpha = 1.D0        ! alpha stabilization (0=off,>0=on)
    i_delay = 0.d0        ! default(usually recommended)
    
    zimp_exp     = PZ(3)  ! Zimp; finite data is necessary
    amassimp_exp = PA(3)  ! Aimp; finite data is necessary
      
    amassgas_exp = PA(2)  ! atomic num. of working hydrogen gas
    
    !--- NR LOOP -----------------------------------------
    DO nr = 1, nrmax-1
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

       diff_jm(nr) = MAX(diffnem,0.d0) 
       chie_jm(nr) = MAX(chietem,0.d0)
       chii_jm(nr) = MAX(chiitim,0.d0)

    ENDDO

    ! *** < nr=1, nrmax-1 on grid >  to  < nr=1,nrmax on half-grid > ***
    ! this subroutine is used for the time being
    call mesh_convert_mtog(diff_jm(1:nrmax-1),diff_jmg(1:nrmax),nrmax-1)
    call mesh_convert_mtog(chie_jm(1:nrmax-1),chie_jmg(1:nrmax),nrmax-1)
    call mesh_convert_mtog(chii_jm(1:nrmax-1),chii_jmg(1:nrmax),nrmax-1)

    DO nsa=1,nsamax
       ns=ns_nsa(nsa)
       nbase=3*(nsa-1)
       IF(pz0(ns) < 0.D0) THEN ! for electron
!          dtr_tb(nbase+1,nbase+1,1:nrmax)=MAX(diff_jmg(1:nrmax),0.d0)
          dtr_tb(nbase+3,nbase+3,1:nrmax)=MAX(chie_jmg(1:nrmax),0.d0)
!          write(*,*) chietem,diffnem
       ELSE 
          IF(pz(ns) /= 0.d0) THEN ! for ion
             dtr_tb(nbase+1,nbase+1,1:nrmax)=MAX(diff_jmg(1:nrmax),0.d0)
             dtr_tb(nbase+3,nbase+3,1:nrmax)=MAX(chii_jmg(1:nrmax),0.d0)    
          END IF
       END IF
    END DO

    RETURN
  END SUBROUTINE tr_glf23
END MODULE trglf23
    
