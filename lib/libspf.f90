MODULE libspf

  PRIVATE
  PUBLIC DGAMMA,DLGAMA,DERF,DERFC,DERFCX,ERF0,ERF1,DPLEG

CONTAINS
!
!S    REAL FUNCTION GAMMA(X)
      DOUBLE PRECISION FUNCTION DGAMMA(X)
!----------------------------------------------------------------------
!
! This routine calculates the GAMMA function for a real argument X.
!   Computation is based on an algorithm outlined in reference 1.
!   The program uses rational functions that approximate the GAMMA
!   function to at least 20 significant decimal digits.  Coefficients
!   for the approximation over the interval (1,2) are unpublished.
!   Those for the approximation for X .GE. 12 are from reference 2.
!   The accuracy achieved depends on the arithmetic system, the
!   compiler, the intrinsic functions, and proper selection of the
!   machine-dependent constants.
!
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
! beta   - radix for the floating-point representation
! maxexp - the smallest positive power of beta that overflows
! XBIG   - the largest argument for which GAMMA(X) is representable
!          in the machine, i.e., the solution to the equation
!                  GAMMA(XBIG) = beta**maxexp
! XINF   - the largest machine representable floating-point number;
!          approximately beta**maxexp
! EPS    - the smallest positive floating-point number such that
!          1.0+EPS .GT. 1.0
! XMININ - the smallest positive floating-point number such that
!          1/XMININ is machine representable
!
!     Approximate values for some important machines are:
!
!                            beta       maxexp        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! Cyber 180/855
!   under NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-Format   (D.P.)        2          127        34.844
! VAX G-Format   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! Cyber 180/855
!   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns the value XINF for singularities or
!     when overflow would occur.  The computation is believed
!     to be free of underflow and overflow.
!
!
!  Intrinsic functions required are:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! References: "An Overview of Software Development for Special
!              Functions", W. J. Cody, Lecture Notes in Mathematics,
!              506, Numerical Analysis Dundee, 1975, G. A. Watson
!              (ed.), Springer Verlag, Berlin, 1976.
!
!              Computer Approximations, Hart, Et. Al., Wiley and
!              sons, New York, 1968.
!
!  Latest modification: October 12, 1989
!
!  Authors: W. J. Cody and L. Stoltz
!           Applied Mathematics Division
!           Argonne National Laboratory
!           Argonne, IL 60439
!
!----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
!S    REAL 
      DOUBLE PRECISION &
           C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE, &
           TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
!----------------------------------------------------------------------
!  Mathematical constants
!----------------------------------------------------------------------
!S    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,
!S   1     SQRTPI/0.9189385332046727417803297E0/,
!S   2     PI/3.1415926535897932384626434E0/
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/, &
           SQRTPI/0.9189385332046727417803297D0/, &
           PI/3.1415926535897932384626434D0/
!----------------------------------------------------------------------
!  Machine dependent parameters
!----------------------------------------------------------------------
!S    DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,
!S   1     XINF/3.4E38/
      DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/, &
           XINF/1.79D308/
!----------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minimax
!     approximation over (1,2).
!----------------------------------------------------------------------
!S    DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,
!S   1       -3.79804256470945635097577E+2,6.29331155312818442661052E+2,
!S   2       8.66966202790413211295064E+2,-3.14512729688483675254357E+4,
!S   3       -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
!S    DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,
!S   1      -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,
!S   2        2.25381184209801510330112E+4,4.75584627752788110767815E+3,
!S   3      -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
      DATA P/-1.71618513886549492533811D+0, 2.47656508055759199108314D+1, &
             -3.79804256470945635097577D+2, 6.29331155312818442661052D+2, &
              8.66966202790413211295064D+2,-3.14512729688483675254357D+4, &
             -3.61444134186911729807069D+4, 6.64561438202405440627855D+4/
      DATA Q/-3.08402300119738975254353D+1, 3.15350626979604161529144D+2, &
             -1.01515636749021914166146D+3,-3.10777167157231109440444D+3, &
              2.25381184209801510330112D+4, 4.75584627752788110767815D+3, &
             -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
!----------------------------------------------------------------------
!  Coefficients for minimax approximation over (12, INF).
!----------------------------------------------------------------------
!S    DATA C/-1.910444077728E-03,8.4171387781295E-04,
!S   1     -5.952379913043012E-04,7.93650793500350248E-04,
!S   2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
!S   3      5.7083835261E-03/
      DATA C/-1.910444077728D-03,8.4171387781295D-04, &
           -5.952379913043012D-04,7.93650793500350248D-04, &
          -2.777777777777681622553D-03,8.333333333333333331554247D-02, &
            5.7083835261D-03/
!----------------------------------------------------------------------
!  Statement functions for conversion between integer and float
!----------------------------------------------------------------------
!S    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .LE. ZERO) THEN
!----------------------------------------------------------------------
!  Argument is negative
!----------------------------------------------------------------------
         Y = -X
         Y1 = AINT(Y)
         RES = Y - Y1
         IF (RES .NE. ZERO) THEN
            IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
            FACT = -PI / SIN(PI*RES)
            Y = Y + ONE
         ELSE
            RES = XINF
            GO TO 900
         END IF
      END IF
!----------------------------------------------------------------------
!  Argument is positive
!----------------------------------------------------------------------
      IF (Y .LT. EPS) THEN
!----------------------------------------------------------------------
!  Argument .LT. EPS
!----------------------------------------------------------------------
         IF (Y .GE. XMININ) THEN
            RES = ONE / Y
         ELSE
            RES = XINF
            GO TO 900
         END IF
      ELSE IF (Y .LT. TWELVE) THEN
         Y1 = Y
         IF (Y .LT. ONE) THEN
!----------------------------------------------------------------------
!  0.0 .LT. argument .LT. 1.0
!----------------------------------------------------------------------
            Z = Y
            Y = Y + ONE
         ELSE
!----------------------------------------------------------------------
!  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
!----------------------------------------------------------------------
            N = INT(Y) - 1
            Y = Y - CONV(N)
            Z = Y - ONE
         END IF
!----------------------------------------------------------------------
!  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
!----------------------------------------------------------------------
         XNUM = ZERO
         XDEN = ONE
         DO I = 1, 8
            XNUM = (XNUM + P(I)) * Z
            XDEN = XDEN * Z + Q(I)
         END DO
         RES = XNUM / XDEN + ONE
         IF (Y1 .LT. Y) THEN
!----------------------------------------------------------------------
!  Adjust result for case  0.0 .LT. argument .LT. 1.0
!----------------------------------------------------------------------
            RES = RES / Y1
         ELSE IF (Y1 .GT. Y) THEN
!----------------------------------------------------------------------
!  Adjust result for case  2.0 .LT. argument .LT. 12.0
!----------------------------------------------------------------------
            DO I = 1, N
               RES = RES * Y
               Y = Y + ONE
            END DO
         END IF
      ELSE
!----------------------------------------------------------------------
!  Evaluate for argument .GE. 12.0,
!----------------------------------------------------------------------
         IF (Y .LE. XBIG) THEN
            YSQ = Y * Y
            SUM = C(7)
            DO I = 1, 6
               SUM = SUM / YSQ + C(I)
               END DO
               SUM = SUM/Y - Y + SQRTPI
               SUM = SUM + (Y-HALF)*LOG(Y)
               RES = EXP(SUM)
            ELSE
               RES = XINF
               GO TO 900
            END IF
         END IF
!----------------------------------------------------------------------
!  Final adjustments and return
!----------------------------------------------------------------------
      IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
!S900 GAMMA = RES
  900 DGAMMA = RES
      RETURN
! ---------- Last line of GAMMA ----------
    END FUNCTION DGAMMA
!
!S    REAL FUNCTION ALGAMA(X)
    DOUBLE PRECISION FUNCTION DLGAMA(X)
!----------------------------------------------------------------------
!
! This routine calculates the LOG(GAMMA) function for a positive real
!   argument X.  Computation is based on an algorithm outlined in
!   references 1 and 2.  The program uses rational functions that
!   theoretically approximate LOG(GAMMA) to at least 18 significant
!   decimal digits.  The approximation for X > 12 is from reference
!   3, while approximations for X < 12.0 are similar to those in
!   reference 1, but are unpublished.  The accuracy achieved depends
!   on the arithmetic system, the compiler, the intrinsic functions,
!   and proper selection of the machine-dependent constants.
!
!
!*********************************************************************
!*********************************************************************
!
! Explanation of machine-dependent constants
!
! beta   - radix for the floating-point representation
! maxexp - the smallest positive power of beta that overflows
! XBIG   - largest argument for which LN(GAMMA(X)) is representable
!          in the machine, i.e., the solution to the equation
!                  LN(GAMMA(XBIG)) = beta**maxexp
! XINF   - largest machine representable floating-point number;
!          approximately beta**maxexp.
! EPS    - The smallest positive floating-point number such that
!          1.0+EPS .GT. 1.0
! FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!     Approximate values for some important machines are:
!
!                           beta      maxexp         XBIG
!
! CRAY-1        (S.P.)        2        8191       9.62E+2461
! Cyber 180/855
!   under NOS   (S.P.)        2        1070       1.72E+319
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)        2         128       4.08E+36
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)        2        1024       2.55D+305
! IBM 3033      (D.P.)       16          63       4.29D+73
! VAX D-Format  (D.P.)        2         127       2.05D+36
! VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                           XINF        EPS        FRTBIG
!
! CRAY-1        (S.P.)   5.45E+2465   7.11E-15    3.13E+615
! Cyber 180/855
!   under NOS   (S.P.)   1.26E+322    3.55E-15    6.44E+79
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)   3.40E+38     1.19E-7     1.42E+9
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)   1.79D+308    2.22D-16    2.25D+76
! IBM 3033      (D.P.)   7.23D+75     2.22D-16    2.56D+18
! VAX D-Format  (D.P.)   1.70D+38     1.39D-17    1.20D+9
! VAX G-Format  (D.P.)   8.98D+307    1.11D-16    1.89D+76
!
!**************************************************************
!**************************************************************
!
! Error returns
!
!  The program returns the value XINF for X .LE. 0.0 or when
!     overflow would occur.  The computation is believed to 
!     be free of underflow and overflow.
!
!
! Intrinsic functions required are:
!
!      LOG
!
!
! References:
!
!  1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for
!     the Natural Logarithm of the Gamma Function,' Math. Comp. 21,
!     1967, pp. 198-203.
!
!  2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May,
!     1969.
! 
!  3) Hart, Et. Al., Computer Approximations, Wiley and sons, New
!     York, 1968.
!
!
!  Authors: W. J. Cody and L. Stoltz
!           Argonne National Laboratory
!
!  Latest modification: June 16, 1988
!
!----------------------------------------------------------------------
      INTEGER I
!S    REAL      
      DOUBLE PRECISION &
           C,CORR,D1,D2,D4,EPS,FRTBIG,FOUR,HALF,ONE,PNT68,P1,P2,P4, &
           Q1,Q2,Q4,RES,SQRTPI,THRHAL,TWELVE,TWO,X,XBIG,XDEN,XINF, &
           XM1,XM2,XM4,XNUM,Y,YSQ,ZERO
      DIMENSION C(7),P1(8),P2(8),P4(8),Q1(8),Q2(8),Q4(8)
!----------------------------------------------------------------------
!  Mathematical constants
!----------------------------------------------------------------------
!S    DATA ONE,HALF,TWELVE,ZERO/1.0E0,0.5E0,12.0E0,0.0E0/,
!S   1     FOUR,THRHAL,TWO,PNT68/4.0E0,1.5E0,2.0E0,0.6796875E0/,
!S   2     SQRTPI/0.9189385332046727417803297E0/
      DATA ONE,HALF,TWELVE,ZERO/1.0D0,0.5D0,12.0D0,0.0D0/, &
           FOUR,THRHAL,TWO,PNT68/4.0D0,1.5D0,2.0D0,0.6796875D0/, &
           SQRTPI/0.9189385332046727417803297D0/
!----------------------------------------------------------------------
!  Machine dependent parameters
!----------------------------------------------------------------------
!S    DATA XBIG,XINF,EPS,FRTBIG/4.08E36,3.401E38,1.19E-7,1.42E9/
      DATA XBIG,XINF,EPS,FRTBIG/2.55D305,1.79D308,2.22D-16,2.25D76/
!----------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minimax
!     approximation over (0.5,1.5).
!----------------------------------------------------------------------
!S    DATA D1/-5.772156649015328605195174E-1/
!S    DATA P1/4.945235359296727046734888E0,2.018112620856775083915565E2,
!S   1        2.290838373831346393026739E3,1.131967205903380828685045E4,
!S   2        2.855724635671635335736389E4,3.848496228443793359990269E4,
!S   3        2.637748787624195437963534E4,7.225813979700288197698961E3/
!S    DATA Q1/6.748212550303777196073036E1,1.113332393857199323513008E3,
!S   1        7.738757056935398733233834E3,2.763987074403340708898585E4,
!S   2        5.499310206226157329794414E4,6.161122180066002127833352E4,
!S   3        3.635127591501940507276287E4,8.785536302431013170870835E3/
      DATA D1/-5.772156649015328605195174D-1/
      DATA P1/4.945235359296727046734888D0,2.018112620856775083915565D2, &
              2.290838373831346393026739D3,1.131967205903380828685045D4, &
              2.855724635671635335736389D4,3.848496228443793359990269D4, &
              2.637748787624195437963534D4,7.225813979700288197698961D3/
      DATA Q1/6.748212550303777196073036D1,1.113332393857199323513008D3, &
              7.738757056935398733233834D3,2.763987074403340708898585D4, &
              5.499310206226157329794414D4,6.161122180066002127833352D4, &
              3.635127591501940507276287D4,8.785536302431013170870835D3/
!----------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minimax
!     Approximation over (1.5,4.0).
!----------------------------------------------------------------------
!S    DATA D2/4.227843350984671393993777E-1/
!S    DATA P2/4.974607845568932035012064E0,5.424138599891070494101986E2,
!S   1        1.550693864978364947665077E4,1.847932904445632425417223E5,
!S   2        1.088204769468828767498470E6,3.338152967987029735917223E6,
!S   3        5.106661678927352456275255E6,3.074109054850539556250927E6/
!S    DATA Q2/1.830328399370592604055942E2,7.765049321445005871323047E3,
!S   1        1.331903827966074194402448E5,1.136705821321969608938755E6,
!S   2        5.267964117437946917577538E6,1.346701454311101692290052E7,
!S   3        1.782736530353274213975932E7,9.533095591844353613395747E6/
      DATA D2/4.227843350984671393993777D-1/
      DATA P2/4.974607845568932035012064D0,5.424138599891070494101986D2, &
              1.550693864978364947665077D4,1.847932904445632425417223D5, &
              1.088204769468828767498470D6,3.338152967987029735917223D6, &
              5.106661678927352456275255D6,3.074109054850539556250927D6/
      DATA Q2/1.830328399370592604055942D2,7.765049321445005871323047D3, &
              1.331903827966074194402448D5,1.136705821321969608938755D6, &
              5.267964117437946917577538D6,1.346701454311101692290052D7, &
              1.782736530353274213975932D7,9.533095591844353613395747D6/
!----------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minimax
!     Approximation over (4.0,12.0).
!----------------------------------------------------------------------
!S    DATA D4/1.791759469228055000094023E0/
!S    DATA P4/1.474502166059939948905062E4,2.426813369486704502836312E6,
!S   1        1.214755574045093227939592E8,2.663432449630976949898078E9,
!S   2      2.940378956634553899906876E10,1.702665737765398868392998E11,
!S   3      4.926125793377430887588120E11,5.606251856223951465078242E11/
!S    DATA Q4/2.690530175870899333379843E3,6.393885654300092398984238E5,
!S   2        4.135599930241388052042842E7,1.120872109616147941376570E9,
!S   3      1.488613728678813811542398E10,1.016803586272438228077304E11,
!S   4      3.417476345507377132798597E11,4.463158187419713286462081E11/
      DATA D4/1.791759469228055000094023D0/
      DATA P4/1.474502166059939948905062D4,2.426813369486704502836312D6, &
              1.214755574045093227939592D8,2.663432449630976949898078D9, &
              2.940378956634553899906876D10,1.702665737765398868392998D11, &
              4.926125793377430887588120D11,5.606251856223951465078242D11/
      DATA Q4/2.690530175870899333379843D3,6.393885654300092398984238D5, &
              4.135599930241388052042842D7,1.120872109616147941376570D9, &
              1.488613728678813811542398D10,1.016803586272438228077304D11, &
              3.417476345507377132798597D11,4.463158187419713286462081D11/
!----------------------------------------------------------------------
!  Coefficients for minimax approximation over (12, INF).
!----------------------------------------------------------------------
!S    DATA C/-1.910444077728E-03,8.4171387781295E-04,
!S   1     -5.952379913043012E-04,7.93650793500350248E-04,
!S   2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
!S   3      5.7083835261E-03/
      DATA C/-1.910444077728D-03,8.4171387781295D-04, &
           -5.952379913043012D-04,7.93650793500350248D-04, &
           -2.777777777777681622553D-03,8.333333333333333331554247D-02, &
            5.7083835261D-03/
!----------------------------------------------------------------------
      Y = X
      IF ((Y .GT. ZERO) .AND. (Y .LE. XBIG)) THEN
         IF (Y .LE. EPS) THEN
            RES = -LOG(Y)
         ELSE IF (Y .LE. THRHAL) THEN
!----------------------------------------------------------------------
!  EPS .LT. X .LE. 1.5
!----------------------------------------------------------------------
            IF (Y .LT. PNT68) THEN
               CORR = -LOG(Y)
               XM1 = Y
            ELSE
               CORR = ZERO
               XM1 = (Y - HALF) - HALF
            END IF
            IF ((Y .LE. HALF) .OR. (Y .GE. PNT68)) THEN
               XDEN = ONE
               XNUM = ZERO
               DO I = 1, 8
                  XNUM = XNUM*XM1 + P1(I)
                  XDEN = XDEN*XM1 + Q1(I)
               END DO
               RES = CORR + (XM1 * (D1 + XM1*(XNUM/XDEN)))
            ELSE
               XM2 = (Y - HALF) - HALF
               XDEN = ONE
               XNUM = ZERO
               DO I = 1, 8
                  XNUM = XNUM*XM2 + P2(I)
                  XDEN = XDEN*XM2 + Q2(I)
               END DO
               RES = CORR + XM2 * (D2 + XM2*(XNUM/XDEN))
            END IF
         ELSE IF (Y .LE. FOUR) THEN
!----------------------------------------------------------------------
!  1.5 .LT. X .LE. 4.0
!----------------------------------------------------------------------
            XM2 = Y - TWO
            XDEN = ONE
            XNUM = ZERO
            DO I = 1, 8
               XNUM = XNUM*XM2 + P2(I)
               XDEN = XDEN*XM2 + Q2(I)
            END DO
            RES = XM2 * (D2 + XM2*(XNUM/XDEN))
         ELSE IF (Y .LE. TWELVE) THEN
!----------------------------------------------------------------------
!  4.0 .LT. X .LE. 12.0
!----------------------------------------------------------------------
            XM4 = Y - FOUR
            XDEN = -ONE
            XNUM = ZERO
            DO I = 1, 8
               XNUM = XNUM*XM4 + P4(I)
               XDEN = XDEN*XM4 + Q4(I)
            END DO
            RES = D4 + XM4*(XNUM/XDEN)
         ELSE 
!----------------------------------------------------------------------
!  Evaluate for argument .GE. 12.0,
!----------------------------------------------------------------------
            RES = ZERO
            IF (Y .LE. FRTBIG) THEN
               RES = C(7)
               YSQ = Y * Y
               DO I = 1, 6
                  RES = RES / YSQ + C(I)
               END DO
            END IF
            RES = RES/Y
            CORR = LOG(Y)
            RES = RES + SQRTPI - HALF*CORR
            RES = RES + Y*(CORR-ONE)
         END IF
      ELSE
!----------------------------------------------------------------------
!  Return for bad arguments
!----------------------------------------------------------------------
         RES = XINF
      END IF
!----------------------------------------------------------------------
!  Final adjustments and return
!----------------------------------------------------------------------
!S    ALGAMA = RES
      DLGAMA = RES
      RETURN
! ---------- Last line of DLGAMA ----------
    END FUNCTION DLGAMA
!
    SUBROUTINE CALERF(ARG,RESULT,ID)
!------------------------------------------------------------------
!
! This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
!   for a real argument  x.  It contains three FUNCTION type
!   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
!   and one SUBROUTINE type subprogram, CALERF.  The calling
!   statements for the primary entries are:
!
!                   Y=ERF(X)     (or   Y=DERF(X)),
!
!                   Y=ERFC(X)    (or   Y=DERFC(X)),
!   and
!                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
!
!   The routine  CALERF  is intended for internal packet use only,
!   all computations within the packet being concentrated in this
!   routine.  The function subprograms invoke  CALERF  with the
!   statement
!
!          CALL CALERF(ARG,RESULT,ID)
!
!   where the parameter usage is as follows
!
!      Function                     Parameters for CALERF
!       call              ARG                  Result          ID
!
!     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
!     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
!     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
!
!   The main computation evaluates near-minimax approximations
!   from "Rational Chebyshev approximations for the error function"
!   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
!   transportable program uses rational functions that theoretically
!   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!   decimal digits.  The accuracy achieved depends on the arithmetic
!   system, the compiler, the intrinsic functions, and proper
!   selection of the machine-dependent constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   XMIN   = the smallest positive floating-point number.
!   XINF   = the largest positive finite floating-point number.
!   XNEG   = the largest negative argument acceptable to ERFCX;
!            the negative of the solution to the equation
!            2*exp(x*x) = XINF.
!   XSMALL = argument below which erf(x) may be represented by
!            2*x/sqrt(pi)  and above which  x*x  will not underflow.
!            A conservative value is the largest machine number X
!            such that   1.0 + X = 1.0   to machine precision.
!   XBIG   = largest argument acceptable to ERFC;  solution to
!            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
!            W(x) = exp(-x*x)/[x*sqrt(pi)].
!   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
!            machine precision.  A conservative value is
!            1/[2*sqrt(XSMALL)]
!   XMAX   = largest acceptable argument to ERFCX; the minimum
!            of XINF and 1/[sqrt(pi)*XMIN].
!
!   Approximate values for some important machines are:
!
!                          XMIN       XINF        XNEG     XSMALL
!
!  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
!  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
!  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
!  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
!  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
!  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
!
!
!                          XBIG       XHUGE       XMAX
!
!  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
!  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
!  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
!  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
!  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
!  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns  ERFC = 0      for  ARG .GE. XBIG;
!
!                       ERFCX = XINF  for  ARG .LT. XNEG;
!      and
!                       ERFCX = 0     for  ARG .GE. XMAX.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, EXP
!
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: March 19, 1990
!
!------------------------------------------------------------------
      INTEGER I,ID
!S    REAL
      DOUBLE PRECISION &
           A,ARG,B,C,D,DEL,FOUR,HALF,P,ONE,Q,RESULT,SIXTEN,SQRPI, &
           TWO,THRESH,X,XBIG,XDEN,XHUGE,XINF,XMAX,XNEG,XNUM,XSMALL, &
           Y,YSQ,ZERO
      DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)
!------------------------------------------------------------------
!  Mathematical constants
!------------------------------------------------------------------
!S    DATA FOUR,ONE,HALF,TWO,ZERO/4.0E0,1.0E0,0.5E0,2.0E0,0.0E0/,
!S   1     SQRPI/5.6418958354775628695E-1/,THRESH/0.46875E0/,
!S   2     SIXTEN/16.0E0/
      DATA FOUR,ONE,HALF,TWO,ZERO/4.0D0,1.0D0,0.5D0,2.0D0,0.0D0/, &
           SQRPI/5.6418958354775628695D-1/,THRESH/0.46875D0/, &
           SIXTEN/16.0D0/
!------------------------------------------------------------------
!  Machine-dependent constants
!------------------------------------------------------------------
!S    DATA XINF,XNEG,XSMALL/3.40E+38,-9.382E0,5.96E-8/,
!S   1     XBIG,XHUGE,XMAX/9.194E0,2.90E3,4.79E37/
      DATA XINF,XNEG,XSMALL/1.79D308,-26.628D0,1.11D-16/, &
           XBIG,XHUGE,XMAX/26.543D0,6.71D7,2.53D307/
!------------------------------------------------------------------
!  Coefficients for approximation to  erf  in first interval
!------------------------------------------------------------------
!S    DATA A/3.16112374387056560E00,1.13864154151050156E02,
!S   1       3.77485237685302021E02,3.20937758913846947E03,
!S   2       1.85777706184603153E-1/
!S    DATA B/2.36012909523441209E01,2.44024637934444173E02,
!S   1       1.28261652607737228E03,2.84423683343917062E03/
      DATA A/3.16112374387056560D00,1.13864154151050156D02, &
             3.77485237685302021D02,3.20937758913846947D03, &
             1.85777706184603153D-1/
      DATA B/2.36012909523441209D01,2.44024637934444173D02, &
             1.28261652607737228D03,2.84423683343917062D03/
!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in second interval
!------------------------------------------------------------------
!S    DATA C/5.64188496988670089E-1,8.88314979438837594E0,
!S   1       6.61191906371416295E01,2.98635138197400131E02,
!S   2       8.81952221241769090E02,1.71204761263407058E03,
!S   3       2.05107837782607147E03,1.23033935479799725E03,
!S   4       2.15311535474403846E-8/
!S    DATA D/1.57449261107098347E01,1.17693950891312499E02,
!S   1       5.37181101862009858E02,1.62138957456669019E03,
!S   2       3.29079923573345963E03,4.36261909014324716E03,
!S   3       3.43936767414372164E03,1.23033935480374942E03/
      DATA C/5.64188496988670089D-1,8.88314979438837594D0, &
             6.61191906371416295D01,2.98635138197400131D02, &
             8.81952221241769090D02,1.71204761263407058D03, &
             2.05107837782607147D03,1.23033935479799725D03, &
             2.15311535474403846D-8/
      DATA D/1.57449261107098347D01,1.17693950891312499D02, &
             5.37181101862009858D02,1.62138957456669019D03, &
             3.29079923573345963D03,4.36261909014324716D03, &
             3.43936767414372164D03,1.23033935480374942D03/
!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in third interval
!------------------------------------------------------------------
!S    DATA P/3.05326634961232344E-1,3.60344899949804439E-1,
!S   1       1.25781726111229246E-1,1.60837851487422766E-2,
!S   2       6.58749161529837803E-4,1.63153871373020978E-2/
!S    DATA Q/2.56852019228982242E00,1.87295284992346047E00,
!S   1       5.27905102951428412E-1,6.05183413124413191E-2,
!S   2       2.33520497626869185E-3/
      DATA P/3.05326634961232344D-1,3.60344899949804439D-1, &
             1.25781726111229246D-1,1.60837851487422766D-2, &
             6.58749161529837803D-4,1.63153871373020978D-2/
      DATA Q/2.56852019228982242D00,1.87295284992346047D00, &
             5.27905102951428412D-1,6.05183413124413191D-2, &
             2.33520497626869185D-3/
!------------------------------------------------------------------
      X = ARG
      Y = ABS(X)
      IF (Y .LE. THRESH) THEN
!------------------------------------------------------------------
!  Evaluate  erf  for  |X| <= 0.46875
!------------------------------------------------------------------
         YSQ = ZERO
         IF (Y .GT. XSMALL) YSQ = Y * Y
         XNUM = A(5)*YSQ
         XDEN = YSQ
         DO I = 1, 3
            XNUM = (XNUM + A(I)) * YSQ
            XDEN = (XDEN + B(I)) * YSQ
         END DO
         RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
         IF (ID .NE. 0) RESULT = ONE - RESULT
         IF (ID .EQ. 2) RESULT = EXP(YSQ) * RESULT
         GO TO 800
!------------------------------------------------------------------
!  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
!------------------------------------------------------------------
      ELSE IF (Y .LE. FOUR) THEN
         XNUM = C(9)*Y
         XDEN = Y
         DO I = 1, 7
            XNUM = (XNUM + C(I)) * Y
            XDEN = (XDEN + D(I)) * Y
         END DO
         RESULT = (XNUM + C(8)) / (XDEN + D(8))
         IF (ID .NE. 2) THEN
            YSQ = AINT(Y*SIXTEN)/SIXTEN
            DEL = (Y-YSQ)*(Y+YSQ)
            RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
         END IF
         !------------------------------------------------------------------
         !  Evaluate  erfc  for |X| > 4.0
!------------------------------------------------------------------
      ELSE
         RESULT = ZERO
         IF (Y .GE. XBIG) THEN
            IF ((ID .NE. 2) .OR. (Y .GE. XMAX)) GO TO 300
            IF (Y .GE. XHUGE) THEN
               RESULT = SQRPI / Y
               GO TO 300
            END IF
         END IF
         YSQ = ONE / (Y * Y)
         XNUM = P(6)*YSQ
         XDEN = YSQ
         DO I = 1, 4
            XNUM = (XNUM + P(I)) * YSQ
            XDEN = (XDEN + Q(I)) * YSQ
         END DO
         RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
         RESULT = (SQRPI -  RESULT) / Y
         IF (ID .NE. 2) THEN
            YSQ = AINT(Y*SIXTEN)/SIXTEN
            DEL = (Y-YSQ)*(Y+YSQ)
            RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
            END IF
      END IF
!------------------------------------------------------------------
!  Fix up for negative argument, erf, etc.
!------------------------------------------------------------------
  300 IF (ID .EQ. 0) THEN
         RESULT = (HALF - RESULT) + HALF
         IF (X .LT. ZERO) RESULT = -RESULT
      ELSE IF (ID .EQ. 1) THEN
         IF (X .LT. ZERO) RESULT = TWO - RESULT
      ELSE
         IF (X .LT. ZERO) THEN
            IF (X .LT. XNEG) THEN
               RESULT = XINF
            ELSE
               YSQ = AINT(X*SIXTEN)/SIXTEN
               DEL = (X-YSQ)*(X+YSQ)
               Y = EXP(YSQ*YSQ) * EXP(DEL)
               RESULT = (Y+Y) - RESULT
            END IF
         END IF
      END IF
  800 RETURN
!---------- Last card of CALERF ----------
    END SUBROUTINE CALERF
!
!S    REAL FUNCTION ERF(X)
    DOUBLE PRECISION FUNCTION DERF(X)
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for erf(x).
!   (see comments heading CALERF).
!
!   Author/date: W. J. Cody, January 8, 1985
!
!--------------------------------------------------------------------
      INTEGER ID
!S    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
!------------------------------------------------------------------
      ID = 0
      CALL CALERF(X,RESULT,ID)
!S    ERF = RESULT
      DERF = RESULT
      RETURN
!---------- Last card of DERF ----------
    END FUNCTION DERF
!
!S    REAL FUNCTION ERFC(X)
    DOUBLE PRECISION FUNCTION DERFC(X)
!--------------------------------------------------------------------
!
! This subprogram computes approximate values for erfc(x).
!   (see comments heading CALERF).
!
!   Author/date: W. J. Cody, January 8, 1985
!
!--------------------------------------------------------------------
      INTEGER ID
!S    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
!------------------------------------------------------------------
      ID = 1
      CALL CALERF(X,RESULT,ID)
!S    ERFC = RESULT
      DERFC = RESULT
      RETURN
!---------- Last card of DERFC ----------
    END FUNCTION DERFC
!
!S    REAL FUNCTION ERFCX(X)
    DOUBLE PRECISION FUNCTION DERFCX(X)
!------------------------------------------------------------------
!
! This subprogram computes approximate values for exp(x*x) * erfc(x).
!   (see comments heading CALERF).
!
!   Author/date: W. J. Cody, March 30, 1987
!
!------------------------------------------------------------------
      INTEGER ID
!S    REAL             X, RESULT
      DOUBLE PRECISION X, RESULT
!------------------------------------------------------------------
      ID = 2
      CALL CALERF(X,RESULT,ID)
!S    ERFCX = RESULT
      DERFCX = RESULT
      RETURN
!---------- Last card of DERFCX ----------
    END FUNCTION DERFCX

! ***********************************************************

!     ERROR FUNCTION

! ***********************************************************

  FUNCTION ERF0(X)
    IMPLICIT NONE
    REAL(8):: ERF0
    REAL(8),INTENT(IN):: X
    REAL(8):: DERF

    ERF0=DERF(X)
    RETURN
  END FUNCTION ERF0

! ***********************************************************

!     DERIVATIVE OF ERROR FUNCTION

! ***********************************************************

  FUNCTION ERF1(U)
    IMPLICIT NONE
    REAL(8):: ERF1
    REAL(8),INTENT(IN):: U
    REAL(8):: U2
    REAL(8),PARAMETER:: PI=3.14159265358979323846D0

    U2=U**2
    IF (U2.GT.100.D0)THEN
       ERF1=0.D0
    ELSE
       ERF1=2/SQRT(PI)*EXP(-U2)
    END IF

    RETURN
  END FUNCTION ERF1

! ***************************************************************

!     LEGENDRE POLYNOMIALS

! ***************************************************************

  SUBROUTINE DPLEG(X,N,Y,IERR)

    IMPLICIT NONE
    REAL(8),INTENT(IN):: X
    INTEGER,INTENT(IN):: N
    INTEGER,INTENT(OUT):: IERR
    REAL(8),INTENT(OUT),DIMENSION(N+1):: Y
    INTEGER:: I
    REAL(8):: W,WY

    IERR = 0
    IF(N .LT. 0) THEN
       WRITE(6,*) 'XX DPLEG: N.LT.0: N=',N
       IERR = 1
       RETURN
    ENDIF

    Y(1) = 1.D0
    IF(N .EQ. 0) RETURN
    Y(2) = X
    IF(N .EQ. 1) RETURN
    DO I=2,N
       W  = X * Y(I)
       WY = W - Y(I-1)
       Y(I+1) = WY + W - WY / DBLE(I)
    ENDDO
    RETURN
  END SUBROUTINE DPLEG
END MODULE libspf
