PROGRAM test

  USE grd2dp_mod,ONLY: GRD2DP
  IMPLICIT NONE
  REAL(8):: x(1000),y(1000),z(1000)
  REAL(4):: R
  INTEGER:: nfmax,nf

  CALL SEED(1952)
  nfmax=100
  DO nf=1,nfmax
     CALL RANDOM(R)
     x(nf)=R
     CALL RANDOM(R)
     y(nf)=R
     z(nf)=SQRT(x(nf)**2+y(nf)**2)
  END DO

  CALL GSOPEN
  CALL PAGES
  CALL GRD2DP(0,x,y,z,nfmax,'@x-y@',NMMAX=1)
  CALL PAGEE

  STOP
END PROGRAM test
