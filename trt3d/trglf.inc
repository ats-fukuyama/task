!     Inputs
      INTEGER(4),DIMENSION(5)  :: itport_pt
      REAL(8),DIMENSION(0:NGLF):: te_m, ti_m, rne_m, rni_m, rns_m, angrotp_exp, egamma_exp, rgamma_p_exp, vphi_m, vpar_m, &
     &                  vper_m, zeff_exp, rho, rgradrho_exp, rgradrhosq_exp, rmin_exp, rmaj_exp, q_exp, shat_exp, &
     &                  alpha_exp, elong_exp

!     Outputs
      REAL(8),DIMENSION(0:NGLF)   :: diff_m, chie_m, chii_m, etaphi_m, etapar_m, etaper_m, exch_m, egamma_m, &
     &                  rgamma_p_m, anrate_m, anrate2_m, anfreq_m, anfreq2_m
      REAL(8),DIMENSION(0:NGLF,10)::egamma_d

!     Auxiliaries
      REAL(8),DIMENSION(0:NGLF):: zpte_m, zpti_m, zpne_m, zpni_m, qe0, qi0, qn0, ddnn, ddne, ddni, chien, chiee, chiei, &
     &                  chiin,chiie,chiii
