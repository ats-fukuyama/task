c@glf.m 11-Apr-01 J. Kinsey, General Atomics
c 05-mar-01 changed 20 to nmode
c 23-aug-00 aligned common block, added xky_gf
c 14-june-00 added ngrow_k_gf
c 13-june-00 added ipert_gf
c 03-aug-99 added jeigen
c---------------------------------------------------------------------
      integer nmode
      parameter (nmode=20)
      integer iflagin_gf(30), ngrow_k_gf(0:nmode)
      integer nroot_gf, jeigen
     & , lprint_gf, ikymax_gf, eigen_gf
     & , i_err, first_order_gf, ipert_gf

      double precision yparam_k_gf(nmode,nmode)
     & , gamma_k_gf(1:4,nmode),freq_k_gf(1:4,nmode)
     & , phi_norm_k_gf(1:4,nmode)
     & , xparam_gf(30)
     & , xkyf_k_gf(nmode),diff_k_gf(nmode)
     & , diff_im_k_gf(nmode),chii_k_gf(nmode)
     & , chie_k_gf(nmode),exch_k_gf(nmode)
     & , eta_par_k_gf(nmode),eta_per_k_gf(nmode),eta_phi_k_gf(nmode)
     & , chie_e_k_gf(nmode), yparam_gf(nmode)
     & , gamma_gf(1:4),freq_gf(1:4),phi_norm_gf(1:4),xky_gf(1:4)
     & , xky0_gf,rms_theta_gf,rlti_gf
     & , rlte_gf,rlne_gf,rlni_gf,rlnimp_gf,dil_gf,apwt_gf
     & , aiwt_gf,taui_gf,rmin_gf,rmaj_gf,q_gf,xnu_gf,betae_gf
     & , shat_gf,alpha_gf,elong_gf,xwell_gf,park_gf,ghat_gf
     & , gchat_gf,adamp_gf,alpha_star_gf,gamma_star_gf
     & , alpha_e_gf,gamma_e_gf,alpha_mode_gf,gamma_mode_gf
     & , alpha_p_gf,gamma_p_gf,xkdamp_gf,xkyf_gf
     & , diff_gf,diff_im_gf,chii_gf,chie_gf,exch_gf,eta_par_gf
     & , eta_per_gf,eta_phi_gf,chie_e_gf
     & , cnorm_gf,xkymin_gf,xkymax_gf,amassgas_gf,amassimp_gf,zimp_gf

C      double precision complex zevec_k_gf(nmode,12,12)
C     & , zomega_k_gf(nmode,12)

C      common /glf/ zevec_k_gf, zomega_k_gf
      common /glf/ 
     &   yparam_k_gf,gamma_k_gf,freq_k_gf,phi_norm_k_gf
     & , xparam_gf,xkyf_k_gf,diff_k_gf
     & , diff_im_k_gf,chii_k_gf,chie_k_gf,exch_k_gf
     & , eta_par_k_gf,eta_per_k_gf,eta_phi_k_gf
     & , chie_e_k_gf,yparam_gf
     & , gamma_gf,freq_gf,phi_norm_gf,xky_gf
     & , xky0_gf,rms_theta_gf
     & , rlti_gf,rlte_gf,rlne_gf,rlni_gf,rlnimp_gf,dil_gf
     & , apwt_gf,aiwt_gf,taui_gf,rmin_gf,rmaj_gf,q_gf,xnu_gf
     & , betae_gf,shat_gf,alpha_gf,elong_gf,xwell_gf,park_gf
     & , ghat_gf,gchat_gf,adamp_gf,alpha_star_gf
     & , gamma_star_gf,alpha_e_gf,gamma_e_gf,alpha_mode_gf
     & , gamma_mode_gf,alpha_p_gf,gamma_p_gf,xkdamp_gf
     & , xkyf_gf,diff_gf,diff_im_gf,chii_gf,chie_gf,exch_gf
     & , eta_par_gf,eta_per_gf,eta_phi_gf,chie_e_gf
     & , cnorm_gf,xkymin_gf,xkymax_gf
     & , amassgas_gf,amassimp_gf,zimp_gf
     & , iflagin_gf,ngrow_k_gf
     & , nroot_gf,lprint_gf,ikymax_gf
     & , i_err,jeigen, eigen_gf
     & , first_order_gf,ipert_gf
