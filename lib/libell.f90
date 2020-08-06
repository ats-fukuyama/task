!     $Id$
      MODULE libell

      private ELLRF, ELLRD
      public ELLFC, ELLEC, ELLF, ELLE

      contains
!
!     ***** ELLIPTIC INTEGRAL LIBRARY *****
!
!     *** Complete elliptic integral of the first kind ***
!
      DOUBLE PRECISION FUNCTION ELLFC(AK,IERR)
!
      IMPLICIT NONE
      REAL(8):: AK
      INTEGER:: IERR
!U    USES ELLRF
      REAL(8):: AAK!, ELLRF
!
      AAK=ABS(AK)
      IF(AAK.GT.1.D0) THEN
         ELLFC=0.D0
         IERR=1
      ELSEIF(AAK.EQ.1.D0) THEN
         ELLFC=3.8718E+01
         IERR=1
      ELSE
         ELLFC=ELLRF(0.D0,(1.D0-AK)*(1.D0+AK),1.D0,IERR)
      ENDIF
      RETURN
      END FUNCTION ELLFC
!
!     *** Complete elliptic integral of the second kind ***
!
      DOUBLE PRECISION FUNCTION ELLEC(AK,IERR)
!
      IMPLICIT NONE
      REAL(8):: AK
      INTEGER:: IERR,IERR1,IERR2
!U    USES ELLRD,ELLRF
      REAL(8):: Q,AAK!,ELLRD,ELLRF
!
      AAK=ABS(AK)
      IF(AAK.GT.1.D0) THEN
         ELLEC=1.D0
         IERR=1
      ELSEIF(AAK.EQ.1.D0) THEN
         ELLEC=1.D0
         IERR=0
      ELSE
         Q=(1.D0-AK)*(1.D0+AK)
         ELLEC=ELLRF(0.D0,Q,1.D0,IERR1) &
              -(AK**2)*ELLRD(0.D0,Q,1.D0,IERR2)/3.D0
         IERR=IERR1+IERR2
      ENDIF
      RETURN
      END FUNCTION ELLEC
!
      DOUBLE PRECISION FUNCTION ELLF(PHI,AK,IERR)
!
      IMPLICIT NONE
      REAL(8):: AK,PHI
      INTEGER:: IERR
!U    USES ELLRF
      REAL(8):: S,Q!,ELLRF
!
      IF(ABS(AK).GT.1.D0) THEN
         ELLF=0.D0
         IERR=1
      ELSE
         S=SIN(PHI)
         Q=(1.D0-AK)*(1.D0+AK)
         ELLF=S*ELLRF(COS(PHI)**2,Q,1.D0,IERR)
      ENDIF
      RETURN
      END FUNCTION ELLF
!
      DOUBLE PRECISION FUNCTION ELLE(PHI,AK,IERR)
      REAL(8):: AK,PHI
      INTEGER:: IERR,IERR1,IERR2
!U    USES ELLRD,ELLRF
      REAL(8):: CC,Q,S!,ELLRD,ELLRF
!
      IF(ABS(AK).GT.1.D0) THEN
         ELLE=0.D0
         IERR=1
      ELSE
         S=SIN(PHI)
         CC=COS(PHI)**2
         Q=(1.D0-S*AK)*(1.D0+S*AK)
         ELLE=S*(ELLRF(CC,Q,1.D0,IERR1) &
              -((S*AK)**2)*ELLRD(CC,Q,1.D0,IERR2)/3.D0)
         IERR=IERR1+IERR2
      ENDIF
      RETURN
      END FUNCTION ELLE
!
!     ------
!
      DOUBLE PRECISION FUNCTION ELLRF (X, Y, Z, IER)
!
!***BEGIN PROLOGUE  DRF
!***PURPOSE  Compute the incomplete or complete elliptic integral of the
!            1st kind.  For X, Y, and Z non-negative and at most one of
!            them zero, RF(X,Y,Z) = Integral from zero to infinity of
!                                -1/2     -1/2     -1/2
!                      (1/2)(t+X)    (t+Y)    (t+Z)    dt.
!            If X, Y or Z is zero, the integral is complete.
!***LIBRARY   SLATEC
!***CATEGORY  C14
!***TYPE      DOUBLE PRECISION (RF-S, DRF-D)
!***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
!             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE FIRST KIND,
!             TAYLOR SERIES
!***AUTHOR  Carlson, B. C.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Notis, E. M.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Pexton, R. L.
!             Lawrence Livermore National Laboratory
!             Livermore, CA  94550
!***DESCRIPTION
!
!   1.     DRF
!          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
!          of the first kind
!          Standard FORTRAN function routine
!          Double precision version
!          The routine calculates an approximation result to
!          DRF(X,Y,Z) = Integral from zero to infinity of
!
!                               -1/2     -1/2     -1/2
!                     (1/2)(t+X)    (t+Y)    (t+Z)    dt,
!
!          where X, Y, and Z are nonnegative and at most one of them
!          is zero.  If one of them  is zero, the integral is COMPLETE.
!          The duplication theorem is iterated until the variables are
!          nearly equal, and the function is then expanded in Taylor
!          series to fifth order.
!
!   2.     Calling sequence
!          DRF( X, Y, Z, IER )
!
!          Parameters On entry
!          Values assigned by the calling routine
!
!          X      - Double precision, nonnegative variable
!
!          Y      - Double precision, nonnegative variable
!
!          Z      - Double precision, nonnegative variable
!
!
!
!          On Return    (values assigned by the DRF routine)
!
!          DRF     - Double precision approximation to the integral
!
!          IER    - Integer
!
!                   IER = 0 Normal and reliable termination of the
!                           routine. It is assumed that the requested
!                           accuracy has been achieved.
!
!                   IER >  0 Abnormal termination of the routine
!
!          X, Y, Z are unaltered.
!
!
!   3.    Error Messages
!
!
!         Value of IER assigned by the DRF routine
!
!                  Value assigned         Error Message Printed
!                  IER = 1                MIN(X,Y,Z) .LT. 0.0D0
!                      = 2                MIN(X+Y,X+Z,Y+Z) .LT. LOLIM
!                      = 3                MAX(X,Y,Z) .GT. UPLIM
!
!
!
!   4.     Control Parameters
!
!                  Values of LOLIM, UPLIM, and ERRTOL are set by the
!                  routine.
!
!          LOLIM and UPLIM determine the valid range of X, Y and Z
!
!          LOLIM  - Lower limit of valid arguments
!
!                   Not less than 5 * (machine minimum).
!
!          UPLIM  - Upper limit of valid arguments
!
!                   Not greater than (machine maximum) / 5.
!
!
!                     Acceptable values for:   LOLIM      UPLIM
!                     IBM 360/370 SERIES   :   3.0D-78     1.0D+75
!                     CDC 6000/7000 SERIES :   1.0D-292    1.0D+321
!                     UNIVAC 1100 SERIES   :   1.0D-307    1.0D+307
!                     CRAY                 :   2.3D-2466   1.09D+2465
!                     VAX 11 SERIES        :   1.5D-38     3.0D+37
!
!
!
!          ERRTOL determines the accuracy of the answer
!
!                 The value assigned by the routine will result
!                 in solution precision within 1-2 decimals of
!                 "machine precision".
!
!
!
!          ERRTOL - Relative error due to truncation is less than
!                   ERRTOL ** 6 / (4 * (1-ERRTOL)  .
!
!
!
!        The accuracy of the computed approximation to the integral
!        can be controlled by choosing the value of ERRTOL.
!        Truncation of a Taylor series after terms of fifth order
!        introduces an error less than the amount shown in the
!        second column of the following table for each value of
!        ERRTOL in the first column.  In addition to the truncation
!        error there will be round-off error, but in practice the
!        total error from both sources is usually less than the
!        amount given in the table.
!
!
!
!
!
!          Sample choices:  ERRTOL   Relative Truncation
!                                    error less than
!                           1.0D-3    3.0D-19
!                           3.0D-3    2.0D-16
!                           1.0D-2    3.0D-13
!                           3.0D-2    2.0D-10
!                           1.0D-1    3.0D-7
!
!
!                    Decreasing ERRTOL by a factor of 10 yields six more
!                    decimal digits of accuracy at the expense of one or
!                    two more iterations of the duplication theorem.
!
! *Long Description:
!
!   DRF Special Comments
!
!
!
!          Check by addition theorem: DRF(X,X+Z,X+W) + DRF(Y,Y+Z,Y+W)
!          = DRF(0,Z,W), where X,Y,Z,W are positive and X * Y = Z * W.
!
!
!          On Input:
!
!          X, Y, and Z are the variables in the integral DRF(X,Y,Z).
!
!
!          On Output:
!
!
!          X, Y, Z are unaltered.
!
!
!
!          ********************************************************
!
!          WARNING: Changes in the program may improve speed at the
!                   expense of robustness.
!
!
!
!   Special double precision functions via DRF
!
!
!
!
!                  Legendre form of ELLIPTIC INTEGRAL of 1st kind
!
!                  -----------------------------------------
!
!
!
!                                             2         2   2
!                  F(PHI,K) = SIN(PHI) DRF(COS (PHI),1-K SIN (PHI),1)
!
!
!                                  2
!                  K(K) = DRF(0,1-K ,1)
!
!
!                         PI/2     2   2      -1/2
!                       = INT  (1-K SIN (PHI) )   D PHI
!                          0
!
!
!
!                  Bulirsch form of ELLIPTIC INTEGRAL of 1st kind
!
!                  -----------------------------------------
!
!
!                                          2 2    2
!                  EL1(X,KC) = X DRF(1,1+KC X ,1+X )
!
!
!                  Lemniscate constant A
!
!                  -----------------------------------------
!
!
!                       1      4 -1/2
!                  A = INT (1-S )    DS = DRF(0,1,2) = DRF(0,2,1)
!                       0
!
!
!
!    -------------------------------------------------------------------
!
!***REFERENCES  B. C. Carlson and E. M. Notis, Algorithms for incomplete
!                 elliptic integrals, ACM Transactions on Mathematical
!                 Software 7, 3 (September 1981), pp. 398-403.
!               B. C. Carlson, Computing elliptic integrals by
!                 duplication, Numerische Mathematik 33, (1979),
!                 pp. 1-16.
!               B. C. Carlson, Elliptic integrals of the first kind,
!                 SIAM Journal of Mathematical Analysis 8, (1977),
!                 pp. 231-242.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900510  Changed calls to XERMSG to standard form, and some
!           editorial changes.  (RWC))
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DRF
      IMPLICIT NONE
      INTEGER IER
      DOUBLE PRECISION LOLIM, UPLIM, EPSLON, ERRTOL, D1MACH
      DOUBLE PRECISION C1, C2, C3, E2, E3, LAMDA
      DOUBLE PRECISION MU, S, X, XN, XNDEV
      DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV, &
           ZNROOT
      LOGICAL FIRST
      SAVE ERRTOL,LOLIM,UPLIM,C1,C2,C3,FIRST
      DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  DRF
!
      IF (FIRST) THEN
!         ERRTOL = (4.0D0*D1MACH(3))**(1.0D0/6.0D0)
!         LOLIM  = 5.0D0 * D1MACH(1)
!         UPLIM  = D1MACH(2)/5.0D0
!
         ERRTOL=1.D-2
         LOLIM=3.0D-78
         UPLIM=1.0D+75
!
         C1 = 1.0D0/24.0D0
         C2 = 3.0D0/44.0D0
         C3 = 1.0D0/14.0D0
      ENDIF
      FIRST = .FALSE.
!
!         CALL ERROR HANDLER IF NECESSARY.
!
      ELLRF = 0.0D0
      IF (MIN(X,Y,Z).LT.0.0D0) THEN
         IER = 1
         WRITE(6,'(A,3ES12.4)')  &
             'XX ELLRF: MIN(X,Y,Z).LT.0.D0: X,Y,Z=',X,Y,Z
         RETURN
      ENDIF
!
      IF (MAX(X,Y,Z).GT.UPLIM) THEN
         IER = 3
         WRITE(6,'(A,3ES12.4)')  &
             'XX ELLRF: MAN(X,Y,Z).GT.UPLIM: X,Y,Z=',X,Y,Z
         RETURN
      ENDIF
!
      IF (MIN(X+Y,X+Z,Y+Z).LT.LOLIM) THEN
         IER = 2
         WRITE(6,'(A,3ES12.4)')  &
             'XX ELLRF: MIN(X+Y,Y+Z,Z+X).LT.LOLIM: X,Y,Z=',X,Y,Z
         RETURN
      ENDIF
!
      IER = 0
      XN = X
      YN = Y
      ZN = Z
!
   30 MU = (XN+YN+ZN)/3.0D0
      XNDEV = 2.0D0 - (MU+XN)/MU
      YNDEV = 2.0D0 - (MU+YN)/MU
      ZNDEV = 2.0D0 - (MU+ZN)/MU
      EPSLON = MAX(ABS(XNDEV),ABS(YNDEV),ABS(ZNDEV))
      IF (EPSLON.LT.ERRTOL) GO TO 40
      XNROOT = SQRT(XN)
      YNROOT = SQRT(YN)
      ZNROOT = SQRT(ZN)
      LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
      XN = (XN+LAMDA)*0.250D0
      YN = (YN+LAMDA)*0.250D0
      ZN = (ZN+LAMDA)*0.250D0
      GO TO 30
!
   40 E2 = XNDEV*YNDEV - ZNDEV*ZNDEV
      E3 = XNDEV*YNDEV*ZNDEV
      S  = 1.0D0 + (C1*E2-0.10D0-C2*E3)*E2 + C3*E3
      ELLRF = S/SQRT(MU)
!
      RETURN
      END FUNCTION ELLRF
!
      DOUBLE PRECISION FUNCTION ELLRD (X, Y, Z, IER)
!
!***BEGIN PROLOGUE  DRD
!***PURPOSE  Compute the incomplete or complete elliptic integral of
!            the 2nd kind. For X and Y nonnegative, X+Y and Z positive,
!            DRD(X,Y,Z) = Integral from zero to infinity of
!                                -1/2     -1/2     -3/2
!                      (3/2)(t+X)    (t+Y)    (t+Z)    dt.
!            If X or Y is zero, the integral is complete.
!***LIBRARY   SLATEC
!***CATEGORY  C14
!***TYPE      DOUBLE PRECISION (RD-S, DRD-D)
!***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
!             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE SECOND KIND,
!             TAYLOR SERIES
!***AUTHOR  Carlson, B. C.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Notis, E. M.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Pexton, R. L.
!             Lawrence Livermore National Laboratory
!             Livermore, CA  94550
!***DESCRIPTION
!
!   1.     DRD
!          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
!          of the second kind
!          Standard FORTRAN function routine
!          Double precision version
!          The routine calculates an approximation result to
!          DRD(X,Y,Z) = Integral from zero to infinity of
!                              -1/2     -1/2     -3/2
!                    (3/2)(t+X)    (t+Y)    (t+Z)    dt,
!          where X and Y are nonnegative, X + Y is positive, and Z is
!          positive.  If X or Y is zero, the integral is COMPLETE.
!          The duplication theorem is iterated until the variables are
!          nearly equal, and the function is then expanded in Taylor
!          series to fifth order.
!
!   2.     Calling Sequence
!
!          DRD( X, Y, Z, IER )
!
!          Parameters On Entry
!          Values assigned by the calling routine
!
!          X      - Double precision, nonnegative variable
!
!          Y      - Double precision, nonnegative variable
!
!                   X + Y is positive
!
!          Z      - Double precision, positive variable
!
!
!
!          On Return    (values assigned by the DRD routine)
!
!          DRD     - Double precision approximation to the integral
!
!
!          IER    - Integer
!
!                   IER = 0 Normal and reliable termination of the
!                           routine. It is assumed that the requested
!                           accuracy has been achieved.
!
!                   IER >  0 Abnormal termination of the routine
!
!
!          X, Y, Z are unaltered.
!
!   3.    Error Messages
!
!         Value of IER assigned by the DRD routine
!
!                  Value assigned         Error message printed
!                  IER = 1                MIN(X,Y) .LT. 0.0D0
!                      = 2                MIN(X + Y, Z ) .LT. LOLIM
!                      = 3                MAX(X,Y,Z) .GT. UPLIM
!
!
!   4.     Control Parameters
!
!                  Values of LOLIM, UPLIM, and ERRTOL are set by the
!                  routine.
!
!          LOLIM and UPLIM determine the valid range of X, Y, and Z
!
!          LOLIM  - Lower limit of valid arguments
!
!                    Not less  than 2 / (machine maximum) ** (2/3).
!
!          UPLIM  - Upper limit of valid arguments
!
!                 Not greater than (0.1D0 * ERRTOL / machine
!                 minimum) ** (2/3), where ERRTOL is described below.
!                 In the following table it is assumed that ERRTOL will
!                 never be chosen smaller than 1.0D-5.
!
!
!                    Acceptable values for:   LOLIM      UPLIM
!                    IBM 360/370 SERIES   :   6.0D-51     1.0D+48
!                    CDC 6000/7000 SERIES :   5.0D-215    2.0D+191
!                    UNIVAC 1100 SERIES   :   1.0D-205    2.0D+201
!                    CRAY                 :   3.0D-1644   1.69D+1640
!                    VAX 11 SERIES        :   1.0D-25     4.5D+21
!
!
!          ERRTOL determines the accuracy of the answer
!
!                 The value assigned by the routine will result
!                 in solution precision within 1-2 decimals of
!                 "machine precision".
!
!          ERRTOL    Relative error due to truncation is less than
!                    3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2.
!
!
!
!        The accuracy of the computed approximation to the integral
!        can be controlled by choosing the value of ERRTOL.
!        Truncation of a Taylor series after terms of fifth order
!        introduces an error less than the amount shown in the
!        second column of the following table for each value of
!        ERRTOL in the first column.  In addition to the truncation
!        error there will be round-off error, but in practice the
!        total error from both sources is usually less than the
!        amount given in the table.
!
!
!
!
!          Sample choices:  ERRTOL   Relative truncation
!                                    error less than
!                           1.0D-3    4.0D-18
!                           3.0D-3    3.0D-15
!                           1.0D-2    4.0D-12
!                           3.0D-2    3.0D-9
!                           1.0D-1    4.0D-6
!
!
!                    Decreasing ERRTOL by a factor of 10 yields six more
!                    decimal digits of accuracy at the expense of one or
!                    two more iterations of the duplication theorem.
!
! *Long Description:
!
!   DRD Special Comments
!
!
!
!          Check: DRD(X,Y,Z) + DRD(Y,Z,X) + DRD(Z,X,Y)
!          = 3 / SQRT(X * Y * Z), where X, Y, and Z are positive.
!
!
!          On Input:
!
!          X, Y, and Z are the variables in the integral DRD(X,Y,Z).
!
!
!          On Output:
!
!
!          X, Y, Z are unaltered.
!
!
!
!          ********************************************************
!
!          WARNING: Changes in the program may improve speed at the
!                   expense of robustness.
!
!
!
!    -------------------------------------------------------------------
!
!
!   Special double precision functions via DRD and DRF
!
!
!                  Legendre form of ELLIPTIC INTEGRAL of 2nd kind
!
!                  -----------------------------------------
!
!
!                                             2         2   2
!                  E(PHI,K) = SIN(PHI) DRF(COS (PHI),1-K SIN (PHI),1) -
!
!                     2      3             2         2   2
!                  -(K/3) SIN (PHI) DRD(COS (PHI),1-K SIN (PHI),1)
!
!
!                                  2        2            2
!                  E(K) = DRF(0,1-K ,1) - (K/3) DRD(0,1-K ,1)
!
!                         PI/2     2   2      1/2
!                       = INT  (1-K SIN (PHI) )  D PHI
!                          0
!
!                  Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind
!
!                  -----------------------------------------
!
!                                               2 2    2
!                  EL2(X,KC,A,B) = AX DRF(1,1+KC X ,1+X ) +
!
!                                              3          2 2    2
!                                 +(1/3)(B-A) X DRD(1,1+KC X ,1+X )
!
!
!
!
!                  Legendre form of alternative ELLIPTIC INTEGRAL
!                  of 2nd kind
!
!                  -----------------------------------------
!
!
!
!                            Q     2       2   2  -1/2
!                  D(Q,K) = INT SIN P  (1-K SIN P)     DP
!                            0
!
!
!
!                                     3          2     2   2
!                  D(Q,K) = (1/3) (SIN Q) DRD(COS Q,1-K SIN Q,1)
!
!
!
!
!                  Lemniscate constant  B
!
!                  -----------------------------------------
!
!
!
!
!                       1    2    4 -1/2
!                  B = INT  S (1-S )    DS
!                       0
!
!
!                  B = (1/3) DRD (0,2,1)
!
!
!                  Heuman's LAMBDA function
!
!                  -----------------------------------------
!
!
!
!                  (PI/2) LAMBDA0(A,B) =
!
!                                    2                2
!                 = SIN(B) (DRF(0,COS (A),1)-(1/3) SIN (A) *
!
!                            2               2         2       2
!                  *DRD(0,COS (A),1)) DRF(COS (B),1-COS (A) SIN (B),1)
!
!                            2       3             2
!                  -(1/3) COS (A) SIN (B) DRF(0,COS (A),1) *
!
!                           2         2       2
!                   *DRD(COS (B),1-COS (A) SIN (B),1)
!
!
!
!                  Jacobi ZETA function
!
!                  -----------------------------------------
!
!                             2                 2       2   2
!                  Z(B,K) = (K/3) SIN(B) DRF(COS (B),1-K SIN (B),1)
!
!
!                                       2             2
!                             *DRD(0,1-K ,1)/DRF(0,1-K ,1)
!
!                               2       3           2       2   2
!                            -(K /3) SIN (B) DRD(COS (B),1-K SIN (B),1)
!
!
! ---------------------------------------------------------------------
!
!***REFERENCES  B. C. Carlson and E. M. Notis, Algorithms for incomplete
!                 elliptic integrals, ACM Transactions on Mathematical
!                 Software 7, 3 (September 1981), pp. 398-403.
!               B. C. Carlson, Computing elliptic integrals by
!                 duplication, Numerische Mathematik 33, (1979),
!                 pp. 1-16.
!               B. C. Carlson, Elliptic integrals of the first kind,
!                 SIAM Journal of Mathematical Analysis 8, (1977),
!                 pp. 231-242.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900510  Modify calls to XERMSG to put in standard form.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DRD
!
      IMPLICIT NONE
      INTEGER IER
      REAL*8 LOLIM, TUPLIM, UPLIM, EPSLON, ERRTOL
      REAL*8 C1, C2, C3, C4, EA, EB, EC, ED, EF, LAMDA
      REAL*8 MU, POWER4, SIGMA, S1, S2, X, XN, XNDEV
      REAL*8 XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV, ZNROOT
      LOGICAL FIRST
      SAVE ERRTOL, LOLIM, UPLIM, C1, C2, C3, C4, FIRST
      DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  DRD
      IF (FIRST) THEN
!         ERRTOL = (D1MACH(3)/3.0D0)**(1.0D0/6.0D0)
!         LOLIM  = 2.0D0/(D1MACH(2))**(2.0D0/3.0D0)
!         TUPLIM = D1MACH(1)**(1.0E0/3.0E0)
!         TUPLIM = (0.10D0*ERRTOL)**(1.0E0/3.0E0)/TUPLIM
!         UPLIM  = TUPLIM**2.0D0
!
         ERRTOL=1.D-2
         LOLIM=3.0D-78
         UPLIM=1.0D+75
!
         C1 = 3.0D0/14.0D0
         C2 = 1.0D0/6.0D0
         C3 = 9.0D0/22.0D0
         C4 = 3.0D0/26.0D0
      ENDIF
      FIRST = .FALSE.
!
!         CALL ERROR HANDLER IF NECESSARY.
!
      ELLRD = 0.0D0
      IF( MIN(X,Y).LT.0.0D0) THEN
         IER = 1
         WRITE(6,'(A,2ES12.4)')  &
             'XX ELLRD: MIN(X,Y).LT.0: X,Y=',X,Y
         RETURN
      ENDIF
!
      IF (MAX(X,Y,Z).GT.UPLIM) THEN
         IER = 3
         WRITE(6,'(A,3ES12.4)')  &
             'XX ELLRD: MAN(X,Y,Z).GT.UPLIM: X,Y,Z=',X,Y,Z
         RETURN
      ENDIF
!
      IF (MIN(X+Y,Z).LT.LOLIM) THEN
         IER = 2
         WRITE(6,'(A,3ES12.4)')  &
             'XX ELLRD: MIN(X+Y,Z).LT.LOLIM: X,Y,Z=',X,Y,Z
         RETURN
      ENDIF
!
      IER = 0
      XN = X
      YN = Y
      ZN = Z
      SIGMA = 0.0D0
      POWER4 = 1.0D0
!
   30 MU = (XN+YN+3.0D0*ZN)*0.20D0
      XNDEV = (MU-XN)/MU
      YNDEV = (MU-YN)/MU
      ZNDEV = (MU-ZN)/MU
      EPSLON = MAX(ABS(XNDEV), ABS(YNDEV), ABS(ZNDEV))
      IF (EPSLON.LT.ERRTOL) GO TO 40
      XNROOT = SQRT(XN)
      YNROOT = SQRT(YN)
      ZNROOT = SQRT(ZN)
      LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
      SIGMA = SIGMA + POWER4/(ZNROOT*(ZN+LAMDA))
      POWER4 = POWER4*0.250D0
      XN = (XN+LAMDA)*0.250D0
      YN = (YN+LAMDA)*0.250D0
      ZN = (ZN+LAMDA)*0.250D0
      GO TO 30
!
   40 EA = XNDEV*YNDEV
      EB = ZNDEV*ZNDEV
      EC = EA - EB
      ED = EA - 6.0D0*EB
      EF = ED + EC + EC
      S1 = ED*(-C1+0.250D0*C3*ED-1.50D0*C4*ZNDEV*EF)
      S2 = ZNDEV*(C2*EF+ZNDEV*(-C3*EC+ZNDEV*C4*EA))
      ELLRD = 3.0D0*SIGMA + POWER4*(1.0D0+S1+S2)/(MU*SQRT(MU))
!
      RETURN
      END FUNCTION ELLRD


      END MODULE libell
