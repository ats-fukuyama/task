MODULE T2CNST 
   
  IMPLICIT NONE

  !C
  INTEGER, PARAMETER :: i0rkind=8
  INTEGER, PARAMETER :: i0ikind=4
  INTEGER, PARAMETER :: i0lmaxm=256
  INTEGER, PARAMETER :: i0spcsm=100
  !C PHYSICAL CONSTANTS
  !C                    based on CODATA 2006
  !C D0PI   : PI
  !C D0AEE  : ELEMENTARY CHARGE 
  !C D0AME  : ELECTRON MASS
  !C D0AMP  : PROTON MASS
  !C D0VC   : SPEED OF LIGHT IN VACCUM
  !C D0RMU0 : PERMEABILITY OF FREE SPACE
  !C D0EPS0 : PERMITTIVITY OF FREE SPACE
  REAL(8), PARAMETER :: d0pi   = 3.14159265358979323846D0
  REAL(8), PARAMETER :: d0aee  = 1.602176487D-19
  REAL(8), PARAMETER :: d0ame  = 9.10938215D-31
  REAL(8), PARAMETER :: d0amp  = 1.672621637D-27
  REAL(8), PARAMETER :: d0vc   = 2.99792458D8
  REAL(8), PARAMETER :: d0rmu0 = 4.D0*d0pi*1.D-7
  REAL(8), PARAMETER :: d0eps0 = 1.D0/(d0vc*d0vc*d0rmu0)
  REAL(8), PARAMETER :: d0vci2  = 1.D0/(d0vc*d0vc)
  !C NUMERICAL COEFFICIENTS FOR EACH VARIABLE 
  !C 
  !C 
  !C
  !C
  !C
!  REAL(8), PARAMETER :: d0bpcst = 1.D0 
!  REAL(8), PARAMETER :: d0btcst = 1.D0 
!  REAL(8), PARAMETER :: d0ercst = 1.D0 
!  REAL(8), PARAMETER :: d0epcst = 1.D0 
!  REAL(8), PARAMETER :: d0etcst = 1.D0 
!  REAL(8), PARAMETER :: d0nncst = 1.D20 
!  REAL(8), PARAMETER :: d0frcst = 1.D20 
!  REAL(8), PARAMETER :: d0fbcst = 1.D20
!  REAL(8), PARAMETER :: d0ftcst = 1.D20
!  REAL(8), PARAMETER :: d0ppcst = d0aee*1.D23 
!  REAL(8), PARAMETER :: d0qrcst = d0aee*1.D23 
!  REAL(8), PARAMETER :: d0qbcst = d0aee*1.D23 
!  REAL(8), PARAMETER :: d0qtcst = d0aee*1.D23 
  REAL(8), PARAMETER :: d0de   = 4.D0/(3.D0*SQRT(d0pi))

 
  !C ABSCISSAS AND WEIGHT FACTOR FOR GAUSIAN INTEGRATION (32 points) 
  !C
  !C FROM HAND BOOK OF MATHEMATICAL FUNCTIONS pp.917
  !C ISBN: 0-486-61272-4
  !C
  !C D0ABSCXX: ABSCISSAS FOR GAUSSIAN INTEGRATION
  !C D0WFCTXX: WEIGHT FACTOR FOR GAUSSIAN INTEGRATION
  REAL(8),DIMENSION(16),PARAMETER::&
       d1absc32(1:16)=(/&
       0.048307665687738D0, 0.144471961582796D0,& 
       0.239287362252137D0, 0.331868602282127D0,&
       0.421351276130635D0, 0.506899908932229D0,&
       0.587715757240762D0, 0.663044266930215D0,&
       0.732182118740289D0, 0.794483795967942D0,&
       0.849367613732569D0, 0.896321155766052D0,&
       0.934906075937739D0, 0.964762255587506D0,&
       0.985611511545268D0, 0.997263861849481D0/),& 
       d1wfct32(1:16)=(/&
       0.096540088514727D0, 0.095638720079274D0,&
       0.093844399080804D0, 0.091173878695763D0,&
       0.087652093004403D0, 0.083311924226946D0,&
       0.078193895787070D0, 0.072345794108848D0,&
       0.065822222776361D0, 0.058684093478535D0,&
       0.050998059262376D0, 0.042835898022226D0,&
       0.034273862913021D0, 0.025392065309262D0,&
       0.016274394730905D0, 0.007018610009470D0/)
END MODULE T2CNST
