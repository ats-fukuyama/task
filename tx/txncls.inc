!Declaration of input to NCLASS
  INTEGER :: k_order, k_potato
  INTEGER :: m_i, m_z
  REAL    :: c_den, c_potb, c_potl
  REAL    :: p_b2, p_bm2, p_eb, p_fhat, p_fm(3), p_ft, p_grbm2, p_ngrth, p_grphi, p_gr2phi
  REAL    :: amu_i(mx_mi), grt_i(mx_mi), temp_i(mx_mi)
  REAL    :: den_iz(mx_mi,mx_mz), fex_iz(3,mx_mi,mx_mz), grp_iz(mx_mi,mx_mz)

!Declaration of output from NCLASS
  INTEGER :: iflag, m_s
  INTEGER :: jm_s(mx_ms), jz_s(mx_ms)
  REAL    :: p_bsjb, p_etap, p_exjb
  REAL    :: calm_i(3,3,mx_mi)
  REAL    :: caln_ii(3,3,mx_mi,mx_mi), capm_ii(3,3,mx_mi,mx_mi), capn_ii(3,3,mx_mi,mx_mi)
  REAL    :: bsjbp_s(mx_ms), bsjbt_s(mx_ms), dn_s(mx_ms), gfl_s(5,mx_ms), qfl_s(5,mx_ms),&
       &     sqz_s(mx_ms), upar_s(3,3,mx_ms), utheta_s(3,3,mx_ms), vn_s(mx_ms), &
       &     veb_s(mx_ms), qeb_s(mx_ms), xi_s(mx_ms), ymu_s(3,3,mx_ms)
  REAL    :: chip_ss(mx_ms,mx_ms), chit_ss(mx_ms,mx_ms), &
       &     dp_ss(mx_ms,mx_ms),      dt_ss(mx_ms,mx_ms)
