! ptcomm.f90

MODULE ptcomm_parm
  USE bpsd_kinds
  USE bpsd_constants

  REAL(rkind):: bb,rr,ra,rkap,rdlt
  INTEGER:: ngxmax,ngymax,nthmax
END MODULE ptcomm_parm

MODULE ptcomm
  USE ptcomm_parm

END MODULE ptcomm
