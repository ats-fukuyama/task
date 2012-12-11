!c@glf.m 11-Apr-01 J. Kinsey, General Atomics
! 03mar2007 pankin converted into F90 module
! 05-mar-01 changed 20 to nmode
! 23-aug-00 aligned common block, added xky_gf
! 14-june-00 added ngrow_k_gf
! 13-june-00 added ipert_gf
! 03-aug-99 added jeigen
!---------------------------------------------------------------------
MODULE glf23_data_mod
!*FD This module contains data structure necessary for glf2d.F
  IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!For REAL variables:
!  Use p=6, r=35   for single precision on 32 bit machine
!  Use p=12,r=100  for double precision on 32 bit machine
!                  for single precision on 64 bit machine
!Parameters for SELECTED_REAL_KIND:
!  p                   -number of digits
!  r                   -range of exponent from -r to r
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER ::                         &
      rspec = SELECTED_REAL_KIND(p=12,r=100),       &
      r8 = SELECTED_REAL_KIND(p=12,r=100)  ,        &
      cspec = SELECTED_REAL_KIND(p=12,r=100),       &
      ispec = SELECTED_INT_KIND(r=6)

  INTEGER nmode
  PARAMETER (nmode=20)

  INTEGER ::iflagin_gf(30)=0, ngrow_k_gf(0:nmode)=0

  INTEGER :: nroot_gf=0

  INTEGER :: jeigen=0           &
             ,lprint_gf=0       &
             ,ikymax_gf=0       &
             ,eigen_gf=0        &
             ,i_err=0           &
             ,first_order_gf=0  &
             ,ipert_gf=0
  !***NOTE: necessary for have quantities with no spaces and comma leading the
  !line for f90doc.py to work.
  !
  !        DOUBLE PRECISION :: yparam_k_gf(nmode,nmode)=0. &
  REAL(r8) :: yparam_k_gf(nmode,nmode)=0. &
             ,gamma_k_gf(1:4,nmode)=0.    &
             ,freq_k_gf(1:4,nmode)=0.     &
             ,phi_norm_k_gf(1:4,nmode)=0. &
             ,xparam_gf(30)=0.            &
             ,xkyf_k_gf(nmode)=0.         &
             ,diff_k_gf(nmode)=0.         &
             ,diff_im_k_gf(nmode)=0.      &
             ,chii_k_gf(nmode)=0.         &
             ,chie_k_gf(nmode)=0.         &
             ,exch_k_gf(nmode)=0.         &
             ,eta_par_k_gf(nmode)=0.      &
             ,eta_per_k_gf(nmode)=0.      &
             ,eta_phi_k_gf(nmode)=0.      &
             ,chie_e_k_gf(nmode)=0.       &
             ,yparam_gf(nmode)=0.         &
             ,gamma_gf(1:4)=0.            &
             ,freq_gf(1:4)=0.             &
             ,phi_norm_gf(1:4)=0.         &
             ,xky_gf(1:4)=0.              &
             ,xky0_gf=0.                  &
             ,rms_theta_gf=0.             &
             ,rlti_gf=0.                  &
             ,rlte_gf=0.                  &
             ,rlne_gf=0.                  &
             ,rlni_gf=0.                  &
             ,rlnimp_gf=0.                &
             ,dil_gf=0.                   &
             ,apwt_gf=0.                  &
             ,aiwt_gf=0.                  &
             ,taui_gf=0.                  &
             ,rmin_gf=0.                  &
             ,rmaj_gf=0.                  &
             ,q_gf=0.                     &
             ,xnu_gf=0.                   &
             ,betae_gf=1.e-6              &
             ,shat_gf=0.                  &
             ,alpha_gf=0.                 &
             ,elong_gf=0.                 &
             ,xwell_gf=0.                 &
             ,park_gf=0.                  &
             ,ghat_gf=0.                  &
             ,gchat_gf=0.                 &
             ,adamp_gf=0.                 &
             ,alpha_star_gf=0.            &
             ,gamma_star_gf=0.            &
             ,alpha_e_gf=0.               &
             ,gamma_e_gf=0.               &
             ,alpha_mode_gf=0.            &
             ,gamma_mode_gf=0.            &
             ,alpha_p_gf=0.               &
             ,gamma_p_gf=0.               &
             ,xkdamp_gf=0.                &
             ,xkyf_gf=0.                  &
             ,diff_gf=0.                  &
             ,diff_im_gf=0.               &
             ,chii_gf=0.                  &
             ,chie_gf=0.                  &
             ,exch_gf=0.                  &
             ,eta_par_gf=0.               &
             ,eta_per_gf=0.               &
             ,eta_phi_gf=0.               &
             ,chie_e_gf=0.                &
             ,cnorm_gf=0.                 &
             ,xkymin_gf=0.                &
             ,xkymax_gf=0.                &
             ,amassgas_gf=0.              &
             ,amassimp_gf=0.              &
             ,zimp_gf=0. 

  COMPLEX(r8) :: zevec_k_gf(nmode,12,12)=0. &
                 ,zomega_k_gf(nmode,12)=0.

CONTAINS

  SUBROUTINE glf_print_input(ionum)
    implicit none
    integer, intent(in), OPTIONAL :: ionum
    INTEGER                       :: filenum

    IF (PRESENT(ionum)) THEN
        filenum=ionum    
    ELSE
        filenum=6
    ENDIF

    write(filenum,*) '==================================================='
    write(filenum,*) 'ipert_gf = ', ipert_gf, ' ngrow_k_gf = ', ngrow_k_gf
    write(filenum,*) 'iflagin_gf = ', iflagin_gf
    write(filenum,*) 'xparam_gf = ', xparam_gf
    write(filenum,*) 'ikymax_gf = ', ikymax_gf
    write(filenum,*) 'nroot_gf = ', nroot_gf
    write(filenum,*)  'xky0_gf = ', xky0_gf
    write(filenum,*)  'rms_theta_gf = ', rms_theta_gf
    write(filenum,*)  'rlti_gf = ', rlti_gf, ' rlte_gf = ', rlte_gf
    write(filenum,*)  'rlne_gf = ', rlne_gf, ' rlni_gf = ', rlni_gf
    write(filenum,*)  'rlnimp_gf = ', rlnimp_gf
    write(filenum,*)  'dil_gf = ', dil_gf
    write(filenum,*)  'apwt_gf = ', apwt_gf, ' aiwt_gf = ', aiwt_gf
    write(filenum,*)  'taui_gf = ', taui_gf
    write(filenum,*)  'rmin_gf = ',  rmin_gf, ' rmaj_gf = ', rmaj_gf
    write(filenum,*)  'q_gf = ', q_gf
    write(filenum,*)  'xnu_gf = ', xnu_gf
    write(filenum,*)  'betae_gf = ', betae_gf
    write(filenum,*)  'shat_gf = ', shat_gf
    write(filenum,*)  'alpha_gf = ', alpha_gf
    write(filenum,*)  'elong_gf = ', elong_gf
    write(filenum,*)  'xwell_gf = ', xwell_gf
    write(filenum,*)  'park_gf = ', park_gf
    write(filenum,*)  'ghat_gf = ', ghat_gf
    write(filenum,*)  'gchat_gf = ', gchat_gf
    write(filenum,*)  'adamp_gf = ', adamp_gf
    write(filenum,*)  'alpha_star_gf = ', alpha_star_gf
    write(filenum,*)  'gamma_star_gf = ', gamma_star_gf
    write(filenum,*)  'alpha_e_gf = ', alpha_e_gf
    write(filenum,*)  'gamma_e_gf = ', gamma_e_gf
    write(filenum,*)  'alpha_mode_gf = ', alpha_mode_gf
    write(filenum,*)  'gamma_mode_gf = ', gamma_mode_gf
    write(filenum,*)  'alpha_p_gf = ', alpha_p_gf
    write(filenum,*)  'gamma_p_gf = ', gamma_p_gf
    write(filenum,*)  'kdamp_gf = ', xkdamp_gf
    write(filenum,*)  'amassgas_gf = ', amassgas_gf
    write(filenum,*)  'amassimp_gf = ', amassimp_gf
    write(filenum,*)  'zimp_gf = ', zimp_gf
    write(filenum,*) '==================================================='

  END SUBROUTINE glf_print_input
  SUBROUTINE glf_print_output(ionum)
    implicit none
    integer, intent(in), OPTIONAL :: ionum
    INTEGER                       :: filenum

    IF (PRESENT(ionum)) THEN
        filenum=ionum    
    ELSE
        filenum=6
    ENDIF


    write(filenum,*) '------OUTPUT---------------------------------------'
    write(filenum,*) 'diff_gf = ', diff_gf
    write(filenum,*) 'diff_im_gf=', diff_im_gf
    write(filenum,*) 'chii_gf=', chii_gf
    write(filenum,*) 'chie_gf=', chie_gf
    write(filenum,*) 'exch_gf=', exch_gf
    write(filenum,*) 'eta_par_gf=', eta_par_gf
    write(filenum,*) 'eta_per_gf=', eta_per_gf
    write(filenum,*) 'eta_phi_gf=', eta_phi_gf
    write(filenum,*) 'chie_e_gf=', chie_e_gf
    write(filenum,*) 'gamma_gf(j1)=', gamma_gf
    write(filenum,*) 'freq_gf(j1)=', freq_gf
    write(filenum,*) 'xky_gf(j1)=', xky_gf
    write(filenum,*) '-----------------------------------------------------'      

  END SUBROUTINE glf_print_output


END MODULE glf23_data_mod
