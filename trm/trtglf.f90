! trtglf.f90

!-------------------------------------------------------------------------
! based on tglf_interface.f90
!
! PURPOSE:
!  Setup input parameters for tglf_run!
!-------------------------------------------------------------------------

MODULE trtglf

  PRIVATE
  PUBLIC tr_tglf
  PUBLIC tr_tglf_test

CONTAINS

  SUBROUTINE tr_tglf_test
    USE trcomm
    USE trprep
    IMPLICIT NONE
    INTEGER:: ierr,nr

    CALL tr_prep(ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,*) 'XX tr_prep ERROR: ierr=',ierr
       STOP
    END IF
    
    CALL tr_tglf

    WRITE(6,'(A,A)') &
         '    NR  pflux:elec  eflux:elec  mflux_elec', &
               '  pflux,ion   eflux:ion   mflux_ion'
    DO nr=1,nrmax
       WRITE(6,'(I6,6ES12.4)') &
            nr,tglf_flux(1,1,nr),tglf_flux(1,2,nr),tglf_flux(1,3,nr), &
               tglf_flux(2,1,nr),tglf_flux(2,2,nr),tglf_flux(2,3,nr)
    END DO
    RETURN
  END SUBROUTINE tr_tglf_test

  SUBROUTINE tr_tglf

    USE trcomm
    IMPLICIT NONE
    INTEGER:: nr,ierr

    DO nr=1,nrmax
       CALL tr_tglf_input(nr)
       CALL tglf_run
       CALL tr_tglf_output(nr,ierr)
       IF(Ierr.NE.0) THEN
          WRITE(6,'(A,2I6)') 'XX tr_tglf: ERROR: nr,ierr=',nr,ierr
          STOP
       END IF
    END DO
    RETURN

  END SUBROUTINE tr_tglf

  SUBROUTINE tr_tglf_input(NR)

    USE trcomm
    USE tglf_interface

    IMPLICIT none

    INTEGER,INTENT(IN):: nr
    INTEGER:: ns

    ! *** following parameters are initialized in module tglf_interface ***
    
    ! === CONTROL PARAMETERS
    ! tglf_path_in       = ''
    ! file_dump_local = 'out.tglf.localdump'
    ! tglf_dump_flag_in  = .false.
    ! tglf_quiet_flag_in = .true.
    ! tglf_test_flag_in  = 0

    ! === units switch
    ! tglf_units_in = 'GYRO'

    ! === INPUT PARAMETERS
    ! tglf_use_transport_model_in = .true.
    ! tglf_geometry_flag_in = 1
    ! tglf_write_wavefunction_flag_in=0

    ! === Data passed to: put_signs
    ! tglf_sign_bt_in        = 1.0
    ! tglf_sign_it_in        = 1.0

    ! === Data passed to: put_rare_switches
    ! tglf_theta_trapped_in  = 0.7
    ! tglf_wdia_trapped_in   = 0.0
    ! tglf_park_in           = 1.0
    ! tglf_ghat_in           = 1.0
    ! tglf_gchat_in          = 1.0
    ! tglf_wd_zero_in        = 0.1
    ! tglf_linsker_factor_in = 0.0
    ! tglf_gradb_factor_in   = 0.0
    ! tglf_filter_in         = 2.0
    ! tglf_damp_psi_in       = 0.0
    ! tglf_damp_sig_in       = 0.0

    ! === Data passed to: put_switches
    ! tglf_iflux_in          = .true.
    ! tglf_use_bper_in       = .false.
    ! tglf_use_bpar_in       = .false.
    ! tglf_use_mhd_rule_in   = .true.
    ! tglf_use_bisection_in  = .true.
    ! tglf_use_inboard_detrapped_in = .false.
    ! tglf_use_ave_ion_grid_in = .false.
    ! tglf_ibranch_in        = -1
    ! tglf_nmodes_in         = 2
    ! tglf_nbasis_max_in     = 4
    ! tglf_nbasis_min_in     = 2
    ! tglf_nxgrid_in         = 16
    ! tglf_nky_in            = 12

    ! === Data passed to: put_model_parameters
    ! tglf_adiabatic_elec_in = .false.
    ! tglf_alpha_mach_in     = 0.0
    ! tglf_alpha_e_in        = 1.0
    ! tglf_alpha_p_in        = 1.0
    ! tglf_alpha_quench_in   = 0.0
    ! tglf_alpha_zf_in       = 1.0
    ! tglf_xnu_factor_in     = 1.0
    ! tglf_debye_factor_in   = 1.0
    ! tglf_etg_factor_in     = 1.25
    ! tglf_rlnp_cutoff_in     = 18.0
    ! tglf_sat_rule_in       = 0
    ! tglf_kygrid_model_in   = 1
    ! tglf_xnu_model_in      = 2
    ! tglf_vpar_model_in     = 0
    ! tglf_vpar_shear_model_in = 1

    ! === Data passed to: put_species
    ! tglf_ns_in             = 2
    ! tglf_mass_in(nsm)      = 0.0
    ! tglf_zs_in(nsm)        = 0.0

    ! ===  Data passed to: put_kys
    ! tglf_ky_in             = 0.3

    ! === Data passed to: put_gaussian_width
    ! tglf_width_in          = 1.65
    ! tglf_width_min_in      = 0.3
    ! tglf_nwidth_in         = 21
    ! tglf_find_width_in     = .true.

    ! === Data passed to: put_gradients
    ! tglf_rlns_in(nsm)      = 0.0
    ! tglf_rlts_in(nsm)      = 0.0
    ! tglf_vpar_shear_in(nsm)= 0.0
    ! tglf_vexb_shear_in     = 0.0

    ! === Data passed to: put_profile_shear
    ! tglf_vns_shear_in(nsm) = 0.0
    ! tglf_vts_shear_in(nsm) = 0.0

    ! === Data passed to: put_averages
    ! tglf_taus_in(nsm)    = 0.0
    ! tglf_as_in(nsm)      = 0.0
    ! tglf_vpar_in(nsm)    = 0.0
    ! tglf_vexb_in         = 0.0
    ! tglf_betae_in        = 0.0
    ! tglf_xnue_in         = 0.0
    ! tglf_zeff_in         = 1.0
    ! tglf_debye_in        = 0.0

    ! === Data passed to: put_eikonal
    ! tglf_new_eikonal_in    = .true.

    ! === Data passed to: put_s_alpha_geometry
    ! tglf_rmin_sa_in        = 0.5
    ! tglf_rmaj_sa_in        = 3.0
    ! tglf_q_sa_in           = 2.0
    ! tglf_shat_sa_in        = 1.0
    ! tglf_alpha_sa_in       = 0.0
    ! tglf_xwell_sa_in       = 0.0
    ! tglf_theta0_sa_in      = 0.0
    ! tglf_b_model_sa_in     = 1
    ! tglf_ft_model_sa_in    = 1

    ! === Data passed to: put_Miller_geometry
    ! tglf_rmin_loc_in       = 0.5
    ! tglf_rmaj_loc_in       = 3.0
    ! tglf_zmaj_loc_in       = 0.0
    ! tglf_drmindx_loc_in    = 1.0
    ! tglf_drmajdx_loc_in    = 0.0
    ! tglf_dzmajdx_loc_in    = 0.0
    ! tglf_kappa_loc_in      = 1.0
    ! tglf_s_kappa_loc_in    = 0.0
    ! tglf_delta_loc_in      = 0.0
    ! tglf_s_delta_loc_in    = 0.0
    ! tglf_zeta_loc_in       = 0.0
    ! tglf_s_zeta_loc_in     = 0.0
    ! tglf_q_loc_in          = 2.0
    ! tglf_q_prime_loc_in    = 16.0
    ! tglf_p_prime_loc_in    = 0.0
    ! tglf_beta_loc_in       = 0.0
    ! tglf_kx0_loc_in        = 0.0

    ! === Data passed to Fourier_geometry
    ! tglf_nfourier_in        = 16
    ! tglf_q_fourier_in          = 2.0
    ! tglf_q_prime_fourier_in    = 16.0
    ! tglf_p_prime_fourier_in    = 0.0
    ! tglf_fourier_in(8,0:max_fourier)=0.0

    ! === Data passed to: put_ELITE_geometry
    ! tglf_n_elite_in   = 700
    ! tglf_q_elite_in   = 2.0
    ! tglf_q_prime_elite_in = 16.0
    ! tglf_R_elite_in(max_ELITE) 
    ! tglf_Z_elite_in(max_ELITE) 
    ! tglf_Bp_elite_in(max_ELITE) 

    ! *** Data passed to: put_signs ***
    tglf_sign_bt_in        = 1.0
    tglf_sign_it_in        = 1.0
    tglf_sat_rule_in       = 0

    ! *** Data passed to: put_species ***
    tglf_ns_in             = nsmax
    tglf_vexb_shear_in     = 0.0        ! ExB velocity shear
    DO ns=1,nsmax
       tglf_mass_in(ns)       = 0.5D0*PA(ns)
       tglf_zs_in(ns)         = PZ(ns)
       tglf_as_in(ns)         = RN(NR,NS)/RN(NR,1) ! n_s/n_e
       tglf_taus_in(ns)       = RT(NR,NS)/RT(NR,1) ! T_s/T_e
       tglf_vpar_in(ns)       = 0.0  ! toroidal velocity R_0 V/R c_s
       tglf_rlns_in(ns)       = 0.0  ! density gradient
       tglf_rlts_in(ns)       = 0.0  ! temperature gradient
       IF(NR.GT.1) THEN
          IF(NR.LT.NRMAX) THEN
             tglf_rlns_in(ns)=RA*(RN(NR+1,NS)-RN(NR,NS))/(RM(NR+1)-RM(NR)) &
                  /(0.5D0*(RN(NR+1,NS)+RN(NR,NS)))
             tglf_rlts_in(ns)=RA*(RT(NR+1,NS)-RT(NR,NS))/(RM(NR+1)-RM(NR)) &
                  /(0.5D0*(RT(NR+1,NS)+RT(NR,NS)))
          ELSE
             tglf_rlns_in(ns)=RA*(RN(NR,NS)-RN(NR-1,NS))/(RM(NR)-RM(NR-1)) &
                  /(0.5D0*(RN(NR,NS)+RN(NR-1,NS)))
             tglf_rlts_in(ns)=RA*(RT(NR,NS)-RT(NR-1,NS))/(RM(NR)-RM(NR-1)) &
                  /(0.5D0*(RT(NR,NS)+RT(NR-1,NS)))
          END IF
       END IF
       tglf_vpar_shear_in(ns) = 0.0  ! toroidal velocity shear R_0 V/R c_s
    END DO
    tglf_vexb_in         = 0.0       ! ExB velocity
    tglf_vexb_shear_in   = 0.0       ! ExB velocity shear

    tglf_betae_in        = 0.0
    tglf_xnue_in         = 0.0
    tglf_zeff_in         = 1.0
    tglf_debye_in        = 0.0

    ! tglf_geometry_flag_in  = 0       ! s-a geometry
    tglf_geometry_flag_in  = 1       ! Miller-a geometry

    ! *** Data passed to: put_s_alpha_geometry ***
    tglf_rmin_sa_in        = RG(NR)
    tglf_rmaj_sa_in        = RR/RA
    tglf_q_sa_in           = ABS(QP(NR))
    tglf_shat_sa_in        = S(NR)
    tglf_alpha_sa_in       = ALPHA(NR)
    tglf_xwell_sa_in       = 0.0
    tglf_theta0_sa_in      = 0.0
    tglf_b_model_sa_in     = 1
    tglf_ft_model_sa_in    = 1

    ! *** Data passed to: put_Miller_geometry ***
    tglf_rmin_loc_in       = RG(NR)
    tglf_rmaj_loc_in       = RR/RA
    tglf_zmaj_loc_in       = 0.0
    tglf_drmindx_loc_in    = 1.0
    tglf_drmajdx_loc_in    = 0.0
    tglf_dzmajdx_loc_in    = 0.0
    tglf_kappa_loc_in      = 1.0
    tglf_s_kappa_loc_in    = 0.0
    tglf_delta_loc_in      = 0.0
    tglf_s_delta_loc_in    = 0.0
    tglf_zeta_loc_in       = 0.0
    tglf_s_zeta_loc_in     = 0.0
    tglf_q_loc_in          = ABS(QP(NR))
    tglf_q_prime_loc_in    = QP(NR)**2*RA**2*S(NR)/RM(NR)**2
    tglf_p_prime_loc_in    = 0.0
    tglf_beta_loc_in       = 0.0
    tglf_kx0_loc_in        = 0.0

    RETURN
  END SUBROUTINE tr_tglf_input

  ! TRANSPORT OUTPUT PARAMETERS

  SUBROUTINE tr_tglf_output(nr,ierr)

    USE trcomm
    USE tglf_interface

    IMPLICIT none
    INTEGER,INTENT(IN):: nr
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: ns

    tglf_flux(1,1,NR)=tglf_elec_pflux_out
    tglf_flux(1,2,NR)=tglf_elec_eflux_out
    tglf_flux(1,3,NR)=tglf_elec_mflux_out
!    tglf_flux(1,4)=tglf_elec_eflux_low_out
!    tglf_flux(1,5)=tglf_elec_expwd_out

    DO ns=2,NSMAX
       tglf_flux(ns,1,NR)=tglf_ion_pflux_out(ns-1)
       tglf_flux(ns,2,NR)=tglf_ion_eflux_out(ns-1)
       tglf_flux(ns,3,NR)=tglf_ion_mflux_out(ns-1)
    END DO

    ! ERROR OUTPUT
    tglf_error_message='null'
    IERR=tglf_error_status
    IF(IERR.NE.0) WRITE(6,*) 'XX '//TRIM(tglf_error_message)

  END SUBROUTINE tr_tglf_output
END MODULE trtglf

