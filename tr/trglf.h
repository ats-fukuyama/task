C     Inputs
      DIMENSION itport_pt(5)
      DIMENSION te_m(0:NGLF),ti_m(0:NGLF),rne_m(0:NGLF),rni_m(0:NGLF),
     &          rns_m(0:NGLF)
      DIMENSION angrotp_exp(0:NGLF),egamma_exp(0:NGLF),
     &          rgamma_p_exp(0:NGLF)
      DIMENSION vphi_m(0:NGLF),vpar_m(0:NGLF),vper_m(0:NGLF)
      DIMENSION zeff_exp(0:NGLF)
      DIMENSION rho(0:NGLF),rgradrho_exp(0:NGLF),rgradrhosq_exp(0:NGLF)
      DIMENSION rmin_exp(0:NGLF),rmaj_exp(0:NGLF)
      DIMENSION q_exp(0:NGLF),shat_exp(0:NGLF),alpha_exp(0:NGLF),
     &          elong_exp(0:NGLF)
C
C     Outputs
      DIMENSION diff_m(0:NGLF),chie_m(0:NGLF),chii_m(0:NGLF)
      DIMENSION etaphi_m(0:NGLF),etapar_m(0:NGLF),etaper_m(0:NGLF)
      DIMENSION exch_m(0:NGLF)
      DIMENSION egamma_m(0:NGLF),egamma_d(0:NGLF,10),rgamma_p_m(0:NGLF)
      DIMENSION anrate_m(0:NGLF),anrate2_m(0:NGLF)
      DIMENSION anfreq_m(0:NGLF),anfreq2_m(0:NGLF)
C
C     Auxiliary
      DIMENSION S_AR(NRM),ALFA_AR(NRM),RKCV_AR(NRM)
      DIMENSION zpte_m(0:NGLF),zpti_m(0:NGLF),zpne_m(0:NGLF),
     &          zpni_m(0:NGLF)
      DIMENSION diff_st(0:NGLF),chie_st(0:NGLF),chii_st(0:NGLF)
      DIMENSION diff_mn(0:NGLF),chie_mn(0:NGLF),chii_mn(0:NGLF)
