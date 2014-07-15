!Declaration of input to NCLASS
      INTEGER(4)::                        k_order,  k_potato, m_i,      m_z
      REAL(4)::                           c_den,    c_potb,   c_potl,   p_b2,     p_bm2,    p_eb,   &
     &                                    p_fhat,   p_ft,     p_grbm2,  p_ngrth,  p_grphi,  p_gr2phi
      REAL(4),DIMENSION(3)            ::  p_fm
      REAL(4),DIMENSION(mx_mi)        ::  amu_i,    grt_i,    temp_i
      REAL(4),DIMENSION(mx_mi,mx_mz)  ::  den_iz,   grp_iz
      REAL(4),DIMENSION(3,mx_mi,mx_mz)::  fex_iz

!Declaration of output from NCLASS
      INTEGER(4)::                        iflag,    m_s
      INTEGER(4),DIMENSION(mx_ms)::       jm_s,     jz_s
      REAL(4)::                           p_bsjb,   p_etap,    p_exjb
      REAL(4),DIMENSION(3,3,mx_mi)::      calm_i
      REAL(4),DIMENSION(3,3,mx_mi,mx_mi)::caln_ii,  capm_ii,  capn_ii
      REAL(4),DIMENSION(mx_ms)::          bsjbp_s,  bsjbt_s,  dn_s,     sqz_s,    vn_s,     veb_s,    qeb_s,    xi_s
      REAL(4),DIMENSION(5,mx_ms)::        gfl_s,    qfl_s
      REAL(4),DIMENSION(3,3,mx_ms)::      upar_s,   utheta_s, ymu_s
      REAL(4),DIMENSION(mx_ms,mx_ms)::    chip_ss,  chit_ss,  dp_ss,    dt_ss
