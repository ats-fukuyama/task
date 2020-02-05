MODULE plprof_travis

PRIVATE

PUBLIC plprof_travis_n,plprof_travis_t

CONTAINS

!     ****** density profile : n=n_0*prof ******

  FUNCTION plprof_travis_n(rhon)
    USE plcomm,ONLY: rkind, &
         profn_travis_g,profn_travis_h,profn_travis_p,profn_travis_q, &
         profn_travis_w
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: rhon
    REAL(rkind):: plprof_travis_n
    
    plprof_travis_n &
         =profn_travis_g-profn_travis_h &
         +(1.D0-profn_travis_g+profn_travis_h) &
         *(1.D0-rhon**profn_travis_p)**profn_travis_q &
         +profn_travis_h*(1.D0-EXP(-rhon**2/profn_travis_w**2))
    RETURN
  END FUNCTION plprof_travis_n

!     ****** temperature profile : T=T_0*prof ******

  FUNCTION plprof_travis_t(rhon)
    USE plcomm,ONLY: rkind, &
         proft_travis_g,proft_travis_h,proft_travis_p,proft_travis_q, &
         proft_travis_w
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: rhon
    REAL(rkind):: plprof_travis_t
    
    plprof_travis_t &
         =proft_travis_g-proft_travis_h &
         +(1.D0-proft_travis_g+proft_travis_h) &
         *(1.D0-rhon**proft_travis_p)**proft_travis_q &
         +proft_travis_h*(1.D0-EXP(-rhon**2/proft_travis_w**2))
    RETURN
  END FUNCTION plprof_travis_t
END MODULE plprof_travis

