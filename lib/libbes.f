C     $Id$
C
C     ****** LIBRARY OF BESSEL FUNCTIONS *****
C
      DOUBLE PRECISION FUNCTION BESJN(N,X)
      REAL*8 X,BJ(101),ALPHA,BESJ0,BESJ1
      DATA ALPHA/0.D0/
C
      NN=ABS(N)
      IF(NN.EQ.0) THEN
         BESJN=BESJ0(X)
      ELSEIF(NN.EQ.1) THEN
         BESJN=BESJ1(X)
      ELSEIF(NN.LE.100) THEN
         CALL DJBESL(X,ALPHA,N+1,BJ,NCALC)
         IF(NCALC.LT.N+1) WRITE(6,*) 'XX BESJN: NCALC=',NCALC
         BESJN=BJ(N+1)
      ELSE
         WRITE(6,*) 'XX BESJN: ABS(N).GT.100: N=',N
         BESJN=0.D0
      ENDIF
      IF(N.LT.0.D0.AND.MOD(NN,2).EQ.1) BESJN=-BESJN
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION BESYN(N,X)
      REAL*8 X,BY(101),ALPHA,BESY0,BESY1
      DATA ALPHA/0.D0/
C
      NN=ABS(N)
      IF(NN.EQ.0) THEN
         BESYN=BESY0(X)
      ELSEIF(NN.EQ.1) THEN
         BESYN=BESY1(X)
      ELSEIF(NN.LE.100) THEN
         CALL DYBESL(X,ALPHA,N+1,BY,NCALC)
         IF(NCALC.LT.N+1) WRITE(6,*) 'XX BESYN: NCALC=',NCALC
         BESYN=BY(N+1)
      ELSE
         WRITE(6,*) 'XX BESYN: ABS(N).GT.100: N=',N
         BESYN=0.D0
      ENDIF
      IF(N.LT.0.D0.AND.MOD(NN,2).EQ.1) BESYN=-BESYN
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION BESIN(N,X)
      REAL*8 X,BI(101),ALPHA,BESI0,BESI1
      DATA ALPHA/0.D0/
C
      NN=ABS(N)
      IF(NN.EQ.0) THEN
         BESIN=BESI0(X)
      ELSEIF(NN.EQ.1) THEN
         BESIN=BESI1(X)
      ELSEIF(NN.LE.100) THEN
         CALL DIBESL(X,ALPHA,N+1,1,BI,NCALC)
         IF(NCALC.LT.N+1) WRITE(6,*) 'XX BESIN: NCALC=',NCALC
         BESIN=BI(N+1)
      ELSE
         WRITE(6,*) 'XX BESIN: ABS(N).GT.100: N=',N
         BESIN=0.D0
      ENDIF
      IF(N.LT.0.D0.AND.MOD(NN,2).EQ.1) BESIN=-BESIN
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION BESKN(N,X)
      REAL*8 X,BK(101),ALPHA,BESK0,BESK1
      DATA ALPHA/0.D0/
C
      NN=ABS(N)
      IF(NN.EQ.0) THEN
         BESKN=BESK0(X)
      ELSEIF(NN.EQ.1) THEN
         BESKN=BESK1(X)
      ELSEIF(NN.LE.100) THEN
         CALL DKBESL(X,ALPHA,N+1,1,BK,NCALC)
         IF(NCALC.LT.N+1) WRITE(6,*) 'XX BESKN: NCALC=',NCALC
         BESKN=BK(N+1)
      ELSE
         WRITE(6,*) 'XX BESKN: ABS(N).GT.100: N=',N
         BESKN=0.D0
      ENDIF
      IF(N.LT.0.D0.AND.MOD(NN,2).EQ.1) BESKN=-BESKN
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION BESEIN(N,X)
      REAL*8 X,BI(101),ALPHA,BESEI0,BESEI1
      DATA ALPHA/0.D0/
C
      NN=ABS(N)
      IF(NN.EQ.0) THEN
         BESEIN=BESEI0(X)
      ELSEIF(NN.EQ.1) THEN
         BESEIN=BESEI1(X)
      ELSEIF(NN.LE.100) THEN
         CALL DIBESL(X,ALPHA,N+1,2,BI,NCALC)
         IF(NCALC.LT.N+1) WRITE(6,*) 'XX BESEIN: NCALC=',NCALC
         BESEIN=BI(N+1)
      ELSE
         WRITE(6,*) 'XX BESEIN: ABS(N).GT.100: N=',N
         BESEIN=0.D0
      ENDIF
      IF(N.LT.0.D0.AND.MOD(NN,2).EQ.1) BESEIN=-BESEIN
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION BESEKN(N,X)
      REAL*8 X,BK(101),ALPHA,BESEK0,BESEK1
      DATA ALPHA/0.D0/
C
      NN=ABS(N)
      IF(NN.EQ.0) THEN
         BESEKN=BESEK0(X)
      ELSEIF(NN.EQ.1) THEN
         BESEKN=BESEK1(X)
      ELSEIF(NN.LE.100) THEN
         CALL DKBESL(X,ALPHA,N+1,2,BK,NCALC)
         IF(NCALC.LT.N+1) WRITE(6,*) 'XX BESEKN: NCALC=',NCALC
         BESEKN=BK(N+1)
      ELSE
         WRITE(6,*) 'XX BESEKN: ABS(N).GT.100: N=',N
         BESEKN=0.D0
      ENDIF
      IF(N.LT.0.D0.AND.MOD(NN,2).EQ.1) BESEKN=-BESEKN
      RETURN
      END
C
      SUBROUTINE BESJNV(N,X,V,IERR)
      REAL*8 X,V(0:N),ALPHA
      DATA ALPHA/0.D0/
C
      IF(N.LT.0) THEN
         IERR=1
         RETURN
      ENDIF
C         
      CALL DJBESL(X,ALPHA,N+1,V,NCALC)
      IF(NCALC.LT.N+1) THEN
         WRITE(6,*) 'XX BESJNV: NCALC,N=',NCALC,N
         IERR=10+NCALC
      ELSE
         IERR=0
      ENDIF
      RETURN
      END
C
      SUBROUTINE BESYNV(N,X,V,IERR)
      REAL*8 X,V(0:N),ALPHA
      DATA ALPHA/0.D0/
C
      IF(N.LT.0) THEN
         IERR=1
         RETURN
      ENDIF
C         
      CALL DYBESL(X,ALPHA,N+1,V,NCALC)
      IF(NCALC.LT.N+1) THEN
         WRITE(6,*) 'XX BESYNV: NCALC,N=',NCALC,N
         IERR=10+NCALC
      ELSE
         IERR=0
      ENDIF
      RETURN
      END
C
      SUBROUTINE BESINV(N,X,V,IERR)
      REAL*8 X,V(0:N),ALPHA
      DATA ALPHA/0.D0/
C
      IF(N.LT.0) THEN
         IERR=1
         RETURN
      ENDIF
C         
      CALL DIBESL(X,ALPHA,N+1,1,V,NCALC)
      IF(NCALC.LT.N+1) THEN
         WRITE(6,*) 'XX BESINV: NCALC,N=',NCALC,N
         IERR=10+NCALC
      ELSE
         IERR=0
      ENDIF
      RETURN
      END
C
      SUBROUTINE BESKNV(N,X,V,IERR)
      REAL*8 X,V(0:N),ALPHA
      DATA ALPHA/0.D0/
C
      IF(N.LT.0) THEN
         IERR=1
         RETURN
      ENDIF
C         
      CALL DKBESL(X,ALPHA,N+1,1,V,NCALC)
      IF(NCALC.LT.N+1) THEN
         WRITE(6,*) 'XX BESKNV: NCALC,N=',NCALC,N
         IERR=10+NCALC
      ELSE
         IERR=0
      ENDIF
      RETURN
      END
C
      SUBROUTINE BESEINV(N,X,V,IERR)
      REAL*8 X,V(0:N),ALPHA
      DATA ALPHA/0.D0/
C
      IF(N.LT.0) THEN
         IERR=1
         RETURN
      ENDIF
C         
      CALL DIBESL(X,ALPHA,N+1,2,V,NCALC)
      IF(NCALC.LT.N+1) THEN
         WRITE(6,*) 'XX BESEINV: NCALC,N=',NCALC,N
         IERR=10+NCALC
      ELSE
         IERR=0
      ENDIF
      RETURN
      END
C
      SUBROUTINE BESEKNV(N,X,V,IERR)
      REAL*8 X,V(0:N),ALPHA
      DATA ALPHA/0.D0/
C
      IF(N.LT.0) THEN
         IERR=1
         RETURN
      ENDIF
C         
      CALL DKBESL(X,ALPHA,N+1,2,V,NCALC)
      IF(NCALC.LT.N+1) THEN
         WRITE(6,*) 'XX BESEKNV: NCALC,N=',NCALC,N
         IERR=10+NCALC
      ELSE
         IERR=0
      ENDIF
      RETURN
      END
C
C     ***** netlib-specfun *****
C
      SUBROUTINE DJBESL(X, ALPHA, NB, B, NCALC)
C---------------------------------------------------------------------
C This routine calculates Bessel functions J sub(N+ALPHA) (X)
C   for non-negative argument X, and non-negative order N+ALPHA.
C
C
C  Explanation of variables in the calling sequence.
C
C   X     - working precision non-negative real argument for which
C           J's are to be calculated.
C   ALPHA - working precision fractional part of order for which
C           J's or exponentially scaled J'r (J*exp(X)) are
C           to be calculated.  0 <= ALPHA < 1.0.
C   NB  - integer number of functions to be calculated, NB > 0.
C           The first function calculated is of order ALPHA, and the
C           last is of order (NB - 1 + ALPHA).
C   B  - working precision output vector of length NB.  If RJBESL
C           terminates normally (NCALC=NB), the vector B contains the
C           functions J/ALPHA/(X) through J/NB-1+ALPHA/(X), or the
C           corresponding exponentially scaled functions.
C   NCALC - integer output variable indicating possible errors.
C           Before using the vector B, the user should check that
C           NCALC=NB, i.e., all orders have been calculated to
C           the desired accuracy.  See Error Returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C  Explanation of machine-dependent constants
C
C   it     = Number of bits in the mantissa of a working precision
C            variable
C   NSIG   = Decimal significance desired.  Should be set to
C            INT(LOG10(2)*it+1).  Setting NSIG lower will result
C            in decreased accuracy while setting NSIG higher will
C            increase CPU time without increasing accuracy.  The
C            truncation error is limited to a relative error of
C            T=.5*10**(-NSIG).
C   ENTEN  = 10.0 ** K, where K is the largest integer such that
C            ENTEN is machine-representable in working precision
C   ENSIG  = 10.0 ** NSIG
C   RTNSIG = 10.0 ** (-K) for the smallest integer K such that
C            K .GE. NSIG/4
C   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
C   XLARGE = Upper limit on the magnitude of X.  If ABS(X)=N,
C            then at least N iterations of the backward recursion
C            will be executed.  The value of 10.0 ** 4 is used on
C            every machine.
C
C
C     Approximate values for some important machines are:
C
C
C                            it    NSIG    ENTEN       ENSIG
C
C   CRAY-1        (S.P.)     48     15    1.0E+2465   1.0E+15
C   Cyber 180/855
C     under NOS   (S.P.)     48     15    1.0E+322    1.0E+15
C   IEEE (IBM/XT,
C     SUN, etc.)  (S.P.)     24      8    1.0E+38     1.0E+8
C   IEEE (IBM/XT,
C     SUN, etc.)  (D.P.)     53     16    1.0D+308    1.0D+16
C   IBM 3033      (D.P.)     14      5    1.0D+75     1.0D+5
C   VAX           (S.P.)     24      8    1.0E+38     1.0E+8
C   VAX D-Format  (D.P.)     56     17    1.0D+38     1.0D+17
C   VAX G-Format  (D.P.)     53     16    1.0D+307    1.0D+16
C
C
C                           RTNSIG      ENMTEN      XLARGE
C
C   CRAY-1        (S.P.)    1.0E-4    1.84E-2466   1.0E+4
C   Cyber 180/855
C     under NOS   (S.P.)    1.0E-4    1.25E-293    1.0E+4
C   IEEE (IBM/XT,
C     SUN, etc.)  (S.P.)    1.0E-2    4.70E-38     1.0E+4
C   IEEE (IBM/XT,
C     SUN, etc.)  (D.P.)    1.0E-4    8.90D-308    1.0D+4
C   IBM 3033      (D.P.)    1.0E-2    2.16D-78     1.0D+4
C   VAX           (S.P.)    1.0E-2    1.17E-38     1.0E+4
C   VAX D-Format  (D.P.)    1.0E-5    1.17D-38     1.0D+4
C   VAX G-Format  (D.P.)    1.0E-4    2.22D-308    1.0D+4
C
C*******************************************************************
C*******************************************************************
C
C  Error returns
C
C    In case of an error,  NCALC .NE. NB, and not all J's are
C    calculated to the desired accuracy.
C
C    NCALC .LT. 0:  An argument is out of range. For example,
C       NBES .LE. 0, ALPHA .LT. 0 or .GT. 1, or X is too large.
C       In this case, B(1) is set to zero, the remainder of the
C       B-vector is not calculated, and NCALC is set to
C       MIN(NB,0)-1 so that NCALC .NE. NB.
C
C    NB .GT. NCALC .GT. 0: Not all requested function values could
C       be calculated accurately.  This usually occurs because NB is
C       much larger than ABS(X).  In this case, B(N) is calculated
C       to the desired accuracy for N .LE. NCALC, but precision
C       is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
C       for N .GT. NCALC (because it is too small to be represented),
C       and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
C       significant figures of B(N) can be trusted.
C
C
C  Intrinsic and other functions required are:
C
C     ABS, AINT, COS, DBLE, GAMMA (or DGAMMA), INT, MAX, MIN,
C
C     REAL, SIN, SQRT
C
C
C  Acknowledgement
C
C   This program is based on a program written by David J. Sookne
C   (2) that computes values of the Bessel functions J or I of real
C   argument and integer order.  Modifications include the restriction
C   of the computation to the J Bessel function of non-negative real
C   argument, the extension of the computation to arbitrary positive
C   order, and the elimination of most underflow.
C
C  References: "A Note on Backward Recurrence Algorithms," Olver,
C               F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
C               pp 941-947.
C
C              "Bessel Functions of Real Argument and Integer Order,"
C               Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
C               125-132.
C
C  Latest modification: March 19, 1990
C
C  Author: W. J. Cody
C          Applied Mathematics Division
C          Argonne National Laboratory
C          Argonne, IL  60439
C
C---------------------------------------------------------------------
      INTEGER I,J,K,L,M,MAGX,N,NB,NBMX,NCALC,NEND,NSTART
CS    REAL               GAMMA,
      DOUBLE PRECISION  DGAMMA,
     1 ALPHA,ALPEM,ALP2EM,B,CAPP,CAPQ,CONV,EIGHTH,EM,EN,ENMTEN,ENSIG,
     2 ENTEN,FACT,FOUR,FUNC,GNU,HALF,HALFX,ONE,ONE30,P,PI2,PLAST,
     3 POLD,PSAVE,PSAVEL,RTNSIG,S,SUM,T,T1,TEMPA,TEMPB,TEMPC,TEST,
     4 THREE,THREE5,TOVER,TWO,TWOFIV,TWOPI1,TWOPI2,X,XC,XIN,XK,XLARGE,
     5 XM,VCOS,VSIN,Z,ZERO
      DIMENSION B(NB), FACT(25)
C---------------------------------------------------------------------
C  Mathematical constants
C
C   PI2    - 2 / PI
C   TWOPI1 - first few significant digits of 2 * PI
C   TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,
C            TWOPI1 + TWOPI2 = 2 * PI to extra precision.
C---------------------------------------------------------------------
CS    DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535E0,6.28125E0,
CS   1 1.935307179586476925286767E-3/
CS    DATA ZERO, EIGHTH, HALF, ONE /0.0E0,0.125E0,0.5E0,1.0E0/
CS    DATA TWO, THREE, FOUR, TWOFIV /2.0E0,3.0E0,4.0E0,25.0E0/
CS    DATA ONE30, THREE5 /130.0E0,35.0E0/
      DATA PI2, TWOPI1, TWOPI2 /0.636619772367581343075535D0,6.28125D0,
     1 1.935307179586476925286767D-3/
      DATA ZERO, EIGHTH, HALF, ONE /0.0D0,0.125D0,0.5D0,1.0D0/
      DATA TWO, THREE, FOUR, TWOFIV /2.0D0,3.0D0,4.0D0,25.0D0/
      DATA ONE30, THREE5 /130.0D0,35.0D0/
C---------------------------------------------------------------------
C  Machine-dependent parameters
C---------------------------------------------------------------------
CS    DATA ENTEN, ENSIG, RTNSIG /1.0E38,1.0E8,1.0E-2/
CS    DATA ENMTEN, XLARGE /1.2E-37,1.0E4/
      DATA ENTEN, ENSIG, RTNSIG /1.0D38,1.0D17,1.0D-4/
      DATA ENMTEN, XLARGE /1.2D-37,1.0D4/
C---------------------------------------------------------------------
C     Factorial(N)
C---------------------------------------------------------------------
CS    DATA FACT /1.0E0,1.0E0,2.0E0,6.0E0,24.0E0,1.2E2,7.2E2,5.04E3,
CS   1 4.032E4,3.6288E5,3.6288E6,3.99168E7,4.790016E8,6.2270208E9,
CS   2 8.71782912E10,1.307674368E12,2.0922789888E13,3.55687428096E14,
CS   3 6.402373705728E15,1.21645100408832E17,2.43290200817664E18,
CS   4 5.109094217170944E19,1.12400072777760768E21,
CS   5 2.585201673888497664E22,6.2044840173323943936E23/
      DATA FACT /1.0D0,1.0D0,2.0D0,6.0D0,24.0D0,1.2D2,7.2D2,5.04D3,
     1 4.032D4,3.6288D5,3.6288D6,3.99168D7,4.790016D8,6.2270208D9,
     2 8.71782912D10,1.307674368D12,2.0922789888D13,3.55687428096D14,
     3 6.402373705728D15,1.21645100408832D17,2.43290200817664D18,
     4 5.109094217170944D19,1.12400072777760768D21,
     5 2.585201673888497664D22,6.2044840173323943936D23/
C---------------------------------------------------------------------
C Statement functions for conversion and the gamma function.
C---------------------------------------------------------------------
CS    CONV(I) = REAL(I)
CS    FUNC(X) = GAMMA(X)
      CONV(I) = DBLE(I)
      FUNC(X) = DGAMMA(X)
C---------------------------------------------------------------------
C Check for out of range arguments.
C---------------------------------------------------------------------
      MAGX = INT(X)
      IF ((NB.GT.0) .AND. (X.GE.ZERO) .AND. (X.LE.XLARGE) 
     1       .AND. (ALPHA.GE.ZERO) .AND. (ALPHA.LT.ONE))  
     2   THEN
C---------------------------------------------------------------------
C Initialize result array to zero.
C---------------------------------------------------------------------
            NCALC = NB
            DO 20 I=1,NB
              B(I) = ZERO
   20       CONTINUE
C---------------------------------------------------------------------
C Branch to use 2-term ascending series for small X and asymptotic
C form for large X when NB is not too large.
C---------------------------------------------------------------------
            IF (X.LT.RTNSIG) THEN
C---------------------------------------------------------------------
C Two-term ascending series for small X.
C---------------------------------------------------------------------
               TEMPA = ONE
               ALPEM = ONE + ALPHA
               HALFX = ZERO
               IF (X.GT.ENMTEN) HALFX = HALF*X
               IF (ALPHA.NE.ZERO)
     1            TEMPA = HALFX**ALPHA/(ALPHA*FUNC(ALPHA))
               TEMPB = ZERO
               IF ((X+ONE).GT.ONE) TEMPB = -HALFX*HALFX
               B(1) = TEMPA + TEMPA*TEMPB/ALPEM
               IF ((X.NE.ZERO) .AND. (B(1).EQ.ZERO)) NCALC = 0
               IF (NB .NE. 1) THEN
                  IF (X .LE. ZERO) THEN
                        DO 30 N=2,NB
                          B(N) = ZERO
   30                   CONTINUE
                     ELSE
C---------------------------------------------------------------------
C Calculate higher order functions.
C---------------------------------------------------------------------
                        TEMPC = HALFX
                        TOVER = (ENMTEN+ENMTEN)/X
                        IF (TEMPB.NE.ZERO) TOVER = ENMTEN/TEMPB
                        DO 50 N=2,NB
                          TEMPA = TEMPA/ALPEM
                          ALPEM = ALPEM + ONE
                          TEMPA = TEMPA*TEMPC
                          IF (TEMPA.LE.TOVER*ALPEM) TEMPA = ZERO
                          B(N) = TEMPA + TEMPA*TEMPB/ALPEM
                          IF ((B(N).EQ.ZERO) .AND. (NCALC.GT.N))
     1                       NCALC = N-1
   50                   CONTINUE
                  END IF
               END IF
            ELSE IF ((X.GT.TWOFIV) .AND. (NB.LE.MAGX+1)) THEN
C---------------------------------------------------------------------
C Asymptotic series for X .GT. 21.0.
C---------------------------------------------------------------------
               XC = SQRT(PI2/X)
               XIN = (EIGHTH/X)**2
               M = 11
               IF (X.GE.THREE5) M = 8
               IF (X.GE.ONE30) M = 4
               XM = FOUR*CONV(M)
C---------------------------------------------------------------------
C Argument reduction for SIN and COS routines.
C---------------------------------------------------------------------
               T = AINT(X/(TWOPI1+TWOPI2)+HALF)
               Z = ((X-T*TWOPI1)-T*TWOPI2) - (ALPHA+HALF)/PI2
               VSIN = SIN(Z)
               VCOS = COS(Z)
               GNU = ALPHA + ALPHA
               DO 80 I=1,2
                 S = ((XM-ONE)-GNU)*((XM-ONE)+GNU)*XIN*HALF
                 T = (GNU-(XM-THREE))*(GNU+(XM-THREE))
                 CAPP = S*T/FACT(2*M+1)
                 T1 = (GNU-(XM+ONE))*(GNU+(XM+ONE))
                 CAPQ = S*T1/FACT(2*M+2)
                 XK = XM
                 K = M + M
                 T1 = T
                 DO 70 J=2,M
                   XK = XK - FOUR
                   S = ((XK-ONE)-GNU)*((XK-ONE)+GNU)
                   T = (GNU-(XK-THREE))*(GNU+(XK-THREE))
                   CAPP = (CAPP+ONE/FACT(K-1))*S*T*XIN
                   CAPQ = (CAPQ+ONE/FACT(K))*S*T1*XIN
                   K = K - 2
                   T1 = T
   70            CONTINUE
                 CAPP = CAPP + ONE
                 CAPQ = (CAPQ+ONE)*(GNU*GNU-ONE)*(EIGHTH/X)
                 B(I) = XC*(CAPP*VCOS-CAPQ*VSIN)
                 IF (NB.EQ.1) GO TO 300
                 T = VSIN
                 VSIN = -VCOS
                 VCOS = T
                 GNU = GNU + TWO
   80         CONTINUE
C---------------------------------------------------------------------
C If  NB .GT. 2, compute J(X,ORDER+I)  I = 2, NB-1
C---------------------------------------------------------------------
               IF (NB .GT. 2) THEN
                  GNU = ALPHA + ALPHA + TWO
                  DO 90 J=3,NB
                    B(J) = GNU*B(J-1)/X - B(J-2)
                    GNU = GNU + TWO
   90             CONTINUE
               END IF
C---------------------------------------------------------------------
C Use recurrence to generate results.  First initialize the
C calculation of P*S.
C---------------------------------------------------------------------
            ELSE
               NBMX = NB - MAGX
               N = MAGX + 1
               EN = CONV(N+N) + (ALPHA+ALPHA)
               PLAST = ONE
               P = EN/X
C---------------------------------------------------------------------
C Calculate general significance test.
C---------------------------------------------------------------------
               TEST = ENSIG + ENSIG
               IF (NBMX .GE. 3) THEN
C---------------------------------------------------------------------
C Calculate P*S until N = NB-1.  Check for possible overflow.
C---------------------------------------------------------------------
                  TOVER = ENTEN/ENSIG
                  NSTART = MAGX + 2
                  NEND = NB - 1
                  EN = CONV(NSTART+NSTART) - TWO + (ALPHA+ALPHA)
                  DO 130 K=NSTART,NEND
                     N = K
                     EN = EN + TWO
                     POLD = PLAST
                     PLAST = P
                     P = EN*PLAST/X - POLD
                     IF (P.GT.TOVER) THEN
C---------------------------------------------------------------------
C To avoid overflow, divide P*S by TOVER.  Calculate P*S until
C ABS(P) .GT. 1.
C---------------------------------------------------------------------
                        TOVER = ENTEN
                        P = P/TOVER
                        PLAST = PLAST/TOVER
                        PSAVE = P
                        PSAVEL = PLAST
                        NSTART = N + 1
  100                   N = N + 1
                           EN = EN + TWO
                           POLD = PLAST
                           PLAST = P
                           P = EN*PLAST/X - POLD
                        IF (P.LE.ONE) GO TO 100
                        TEMPB = EN/X
C---------------------------------------------------------------------
C Calculate backward test and find NCALC, the highest N such that
C the test is passed.
C---------------------------------------------------------------------
                        TEST = POLD*PLAST*(HALF-HALF/(TEMPB*TEMPB))
                        TEST = TEST/ENSIG
                        P = PLAST*TOVER
                        N = N - 1
                        EN = EN - TWO
                        NEND = MIN(NB,N)
                        DO 110 L=NSTART,NEND
                           POLD = PSAVEL
                           PSAVEL = PSAVE
                           PSAVE = EN*PSAVEL/X - POLD
                           IF (PSAVE*PSAVEL.GT.TEST) THEN
                              NCALC = L - 1
                              GO TO 190
                           END IF
  110                   CONTINUE
                        NCALC = NEND
                        GO TO 190
                     END IF
  130             CONTINUE
                  N = NEND
                  EN = CONV(N+N) + (ALPHA+ALPHA)
C---------------------------------------------------------------------
C Calculate special significance test for NBMX .GT. 2.
C---------------------------------------------------------------------
                  TEST = MAX(TEST,SQRT(PLAST*ENSIG)*SQRT(P+P))
               END IF
C---------------------------------------------------------------------
C Calculate P*S until significance test passes.
C---------------------------------------------------------------------
  140          N = N + 1
                  EN = EN + TWO
                  POLD = PLAST
                  PLAST = P
                  P = EN*PLAST/X - POLD
               IF (P.LT.TEST) GO TO 140
C---------------------------------------------------------------------
C Initialize the backward recursion and the normalization sum.
C---------------------------------------------------------------------
  190          N = N + 1
               EN = EN + TWO
               TEMPB = ZERO
               TEMPA = ONE/P
               M = 2*N - 4*(N/2)
               SUM = ZERO
               EM = CONV(N/2)
               ALPEM = (EM-ONE) + ALPHA
               ALP2EM = (EM+EM) + ALPHA
               IF (M .NE. 0) SUM = TEMPA*ALPEM*ALP2EM/EM
               NEND = N - NB
               IF (NEND .GT. 0) THEN
C---------------------------------------------------------------------
C Recur backward via difference equation, calculating (but not
C storing) B(N), until N = NB.
C---------------------------------------------------------------------
                  DO 200 L=1,NEND
                     N = N - 1
                     EN = EN - TWO
                     TEMPC = TEMPB
                     TEMPB = TEMPA
                     TEMPA = (EN*TEMPB)/X - TEMPC
                     M = 2 - M
                     IF (M .NE. 0) THEN
                        EM = EM - ONE
                        ALP2EM = (EM+EM) + ALPHA
                        IF (N.EQ.1) GO TO 210
                        ALPEM = (EM-ONE) + ALPHA
                        IF (ALPEM.EQ.ZERO) ALPEM = ONE
                        SUM = (SUM+TEMPA*ALP2EM)*ALPEM/EM
                     END IF
  200             CONTINUE
               END IF
C---------------------------------------------------------------------
C Store B(NB).
C---------------------------------------------------------------------
  210          B(N) = TEMPA
               IF (NEND .GE. 0) THEN
                  IF (NB .LE. 1) THEN
                        ALP2EM = ALPHA
                        IF ((ALPHA+ONE).EQ.ONE) ALP2EM = ONE
                        SUM = SUM + B(1)*ALP2EM
                        GO TO 250
                     ELSE
C---------------------------------------------------------------------
C Calculate and store B(NB-1).
C---------------------------------------------------------------------
                        N = N - 1
                        EN = EN - TWO
                        B(N) = (EN*TEMPA)/X - TEMPB
                        IF (N.EQ.1) GO TO 240
                        M = 2 - M
                        IF (M .NE. 0) THEN
                           EM = EM - ONE
                           ALP2EM = (EM+EM) + ALPHA
                           ALPEM = (EM-ONE) + ALPHA
                           IF (ALPEM.EQ.ZERO) ALPEM = ONE
                           SUM = (SUM+B(N)*ALP2EM)*ALPEM/EM
                        END IF
                  END IF
               END IF
               NEND = N - 2
               IF (NEND .NE. 0) THEN
C---------------------------------------------------------------------
C Calculate via difference equation and store B(N), until N = 2.
C---------------------------------------------------------------------
                  DO 230 L=1,NEND
                     N = N - 1
                     EN = EN - TWO
                     B(N) = (EN*B(N+1))/X - B(N+2)
                     M = 2 - M
                     IF (M .NE. 0) THEN
                        EM = EM - ONE
                        ALP2EM = (EM+EM) + ALPHA
                        ALPEM = (EM-ONE) + ALPHA
                        IF (ALPEM.EQ.ZERO) ALPEM = ONE
                        SUM = (SUM+B(N)*ALP2EM)*ALPEM/EM
                     END IF
  230             CONTINUE
               END IF
C---------------------------------------------------------------------
C Calculate B(1).
C---------------------------------------------------------------------
               B(1) = TWO*(ALPHA+ONE)*B(2)/X - B(3)
  240          EM = EM - ONE
               ALP2EM = (EM+EM) + ALPHA
               IF (ALP2EM.EQ.ZERO) ALP2EM = ONE
               SUM = SUM + B(1)*ALP2EM
C---------------------------------------------------------------------
C Normalize.  Divide all B(N) by sum.
C---------------------------------------------------------------------
  250          IF ((ALPHA+ONE).NE.ONE)
     1              SUM = SUM*FUNC(ALPHA)*(X*HALF)**(-ALPHA)
               TEMPA = ENMTEN
               IF (SUM.GT.ONE) TEMPA = TEMPA*SUM
               DO 260 N=1,NB
                 IF (ABS(B(N)).LT.TEMPA) B(N) = ZERO
                 B(N) = B(N)/SUM
  260          CONTINUE
            END IF
C---------------------------------------------------------------------
C Error return -- X, NB, or ALPHA is out of range.
C---------------------------------------------------------------------
         ELSE
            B(1) = ZERO
            NCALC = MIN(NB,0) - 1
      END IF
C---------------------------------------------------------------------
C Exit
C---------------------------------------------------------------------
  300 RETURN
C ---------- Last line of RJBESL ----------
      END
C
      SUBROUTINE DYBESL(X,ALPHA,NB,BY,NCALC)
C----------------------------------------------------------------------
C
C  This routine calculates Bessel functions Y SUB(N+ALPHA) (X)
C  for non-negative argument X, and non-negative order N+ALPHA.
C
C
C Explanation of variables in the calling sequence
C
C X     - Working precision non-negative real argument for which
C         Y's are to be calculated.
C ALPHA - Working precision fractional part of order for which
C         Y's are to be calculated.  0 .LE. ALPHA .LT. 1.0.
C NB    - Integer number of functions to be calculated, NB .GT. 0.
C         The first function calculated is of order ALPHA, and the 
C         last is of order (NB - 1 + ALPHA).
C BY    - Working precision output vector of length NB.  If the
C         routine terminates normally (NCALC=NB), the vector BY
C         contains the functions Y(ALPHA,X), ... , Y(NB-1+ALPHA,X),
C         If (0 .LT. NCALC .LT. NB), BY(I) contains correct function
C         values for I .LE. NCALC, and contains the ratios
C         Y(ALPHA+I-1,X)/Y(ALPHA+I-2,X) for the rest of the array.
C NCALC - Integer output variable indicating possible errors.
C         Before using the vector BY, the user should check that 
C         NCALC=NB, i.e., all orders have been calculated to
C         the desired accuracy.  See error returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta   = Radix for the floating-point system
C   p      = Number of significant base-beta digits in the
C            significand of a floating-point number
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C   EPS    = beta ** (-p)
C   DEL    = Machine number below which sin(x)/x = 1; approximately
C            SQRT(EPS).
C   XMIN   = Smallest acceptable argument for RBESY; approximately
C            max(2*beta**minexp,2/XINF), rounded up
C   XINF   = Largest positive machine number; approximately
C            beta**maxexp
C   THRESH = Lower bound for use of the asymptotic form; approximately
C            AINT(-LOG10(EPS/2.0))+1.0
C   XLARGE = Upper bound on X; approximately 1/DEL, because the sine
C            and cosine functions have lost about half of their 
C            precision at that point.
C
C
C     Approximate values for some important machines are:
C
C                        beta    p     minexp      maxexp      EPS
C
C  CRAY-1        (S.P.)    2    48     -8193        8191    3.55E-15
C  Cyber 180/185 
C    under NOS   (S.P.)    2    48      -975        1070    3.55E-15
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)    2    24      -126         128    5.96E-8
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)    2    53     -1022        1024    1.11D-16
C  IBM 3033      (D.P.)   16    14       -65          63    1.39D-17
C  VAX           (S.P.)    2    24      -128         127    5.96E-8
C  VAX D-Format  (D.P.)    2    56      -128         127    1.39D-17
C  VAX G-Format  (D.P.)    2    53     -1024        1023    1.11D-16
C
C
C                         DEL      XMIN      XINF     THRESH  XLARGE
C
C CRAY-1        (S.P.)  5.0E-8  3.67E-2466 5.45E+2465  15.0E0  2.0E7
C Cyber 180/855
C   under NOS   (S.P.)  5.0E-8  6.28E-294  1.26E+322   15.0E0  2.0E7
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)  1.0E-4  2.36E-38   3.40E+38     8.0E0  1.0E4
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)  1.0D-8  4.46D-308  1.79D+308   16.0D0  1.0D8
C IBM 3033      (D.P.)  1.0D-8  2.77D-76   7.23D+75    17.0D0  1.0D8
C VAX           (S.P.)  1.0E-4  1.18E-38   1.70E+38     8.0E0  1.0E4
C VAX D-Format  (D.P.)  1.0D-9  1.18D-38   1.70D+38    17.0D0  1.0D9
C VAX G-Format  (D.P.)  1.0D-8  2.23D-308  8.98D+307   16.0D0  1.0D8
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  In case of an error, NCALC .NE. NB, and not all Y's are
C  calculated to the desired accuracy.
C
C  NCALC .LT. -1:  An argument is out of range. For example,
C       NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
C       XMAX.  In this case, BY(1) = 0.0, the remainder of the
C       BY-vector is not calculated, and NCALC is set to
C       MIN0(NB,0)-2  so that NCALC .NE. NB.
C  NCALC = -1:  Y(ALPHA,X) .GE. XINF.  The requested function
C       values are set to 0.0.
C  1 .LT. NCALC .LT. NB: Not all requested function values could
C       be calculated accurately.  BY(I) contains correct function
C       values for I .LE. NCALC, and and the remaining NB-NCALC
C       array elements contain 0.0.
C
C
C Intrinsic functions required are:
C
C     DBLE, EXP, INT, MAX, MIN, REAL, SQRT
C
C
C Acknowledgement
C
C  This program draws heavily on Temme's Algol program for Y(a,x)
C  and Y(a+1,x) and on Campbell's programs for Y_nu(x).  Temme's
C  scheme is used for  x < THRESH, and Campbell's scheme is used
C  in the asymptotic region.  Segments of code from both sources
C  have been translated into Fortran 77, merged, and heavily modified.
C  Modifications include parameterization of machine dependencies,
C  use of a new approximation for ln(gamma(x)), and built-in
C  protection against over/underflow.
C
C References: "Bessel functions J_nu(x) and Y_nu(x) of real
C              order and real argument," Campbell, J. B.,
C              Comp. Phy. Comm. 18, 1979, pp. 133-142.
C
C             "On the numerical evaluation of the ordinary
C              Bessel function of the second kind," Temme,
C              N. M., J. Comput. Phys. 21, 1976, pp. 343-350.
C
C  Latest modification: March 19, 1990
C
C  Modified by: W. J. Cody
C               Applied Mathematics Division
C               Argonne National Laboratory
C               Argonne, IL  60439
C
C----------------------------------------------------------------------
      INTEGER I,K,NA,NB,NCALC
CS    REAL
      DOUBLE PRECISION
     1  ALFA,ALPHA,AYE,B,BY,C,CH,COSMU,D,DEL,DEN,DDIV,DIV,DMU,D1,D2,
     2  E,EIGHT,EN,ENU,EN1,EPS,EVEN,EX,F,FIVPI,G,GAMMA,H,HALF,ODD,
     3  ONBPI,ONE,ONE5,P,PA,PA1,PI,PIBY2,PIM5,Q,QA,QA1,Q0,R,S,SINMU,
     4  SQ2BPI,TEN9,TERM,THREE,THRESH,TWO,TWOBYX,X,XINF,XLARGE,XMIN,
     5  XNA,X2,YA,YA1,ZERO
      DIMENSION BY(NB),CH(21)
C----------------------------------------------------------------------
C  Mathematical constants
C    FIVPI = 5*PI
C    PIM5 = 5*PI - 15
C    ONBPI = 1/PI
C    PIBY2 = PI/2
C    SQ2BPI = SQUARE ROOT OF 2/PI
C----------------------------------------------------------------------
CS    DATA ZERO,HALF,ONE,TWO,THREE/0.0E0,0.5E0,1.0E0,2.0E0,3.0E0/
CS    DATA EIGHT,ONE5,TEN9/8.0E0,15.0E0,1.9E1/
CS    DATA FIVPI,PIBY2/1.5707963267948966192E1,1.5707963267948966192E0/
CS    DATA PI,SQ2BPI/3.1415926535897932385E0,7.9788456080286535588E-1/
CS    DATA PIM5,ONBPI/7.0796326794896619231E-1,3.1830988618379067154E-1/
      DATA ZERO,HALF,ONE,TWO,THREE/0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/
      DATA EIGHT,ONE5,TEN9/8.0D0,15.0D0,1.9D1/
      DATA FIVPI,PIBY2/1.5707963267948966192D1,1.5707963267948966192D0/
      DATA PI,SQ2BPI/3.1415926535897932385D0,7.9788456080286535588D-1/
      DATA PIM5,ONBPI/7.0796326794896619231D-1,3.1830988618379067154D-1/
C----------------------------------------------------------------------
C  Machine-dependent constants
C----------------------------------------------------------------------
CS    DATA DEL,XMIN,XINF,EPS/1.0E-4,2.36E-38,3.40E38,5.96E-8/
CS    DATA THRESH,XLARGE/8.0E0,1.0E4/
      DATA DEL,XMIN,XINF,EPS/1.0D-8,4.46D-308,1.79D308,1.11D-16/
      DATA THRESH,XLARGE/16.0D0,1.0D8/
C----------------------------------------------------------------------
C  Coefficients for Chebyshev polynomial expansion of 
C         1/gamma(1-x), abs(x) .le. .5
C----------------------------------------------------------------------
CS    DATA CH/-0.67735241822398840964E-23,-0.61455180116049879894E-22,
CS   1         0.29017595056104745456E-20, 0.13639417919073099464E-18,
CS   2         0.23826220476859635824E-17,-0.90642907957550702534E-17,
CS   3        -0.14943667065169001769E-14,-0.33919078305362211264E-13,
CS   4        -0.17023776642512729175E-12, 0.91609750938768647911E-11,
CS   5         0.24230957900482704055E-09, 0.17451364971382984243E-08,
CS   6        -0.33126119768180852711E-07,-0.86592079961391259661E-06,
CS   7        -0.49717367041957398581E-05, 0.76309597585908126618E-04,
CS   8         0.12719271366545622927E-02, 0.17063050710955562222E-02,
CS   9        -0.76852840844786673690E-01,-0.28387654227602353814E+00,
CS   A         0.92187029365045265648E+00/
      DATA CH/-0.67735241822398840964D-23,-0.61455180116049879894D-22,
     1         0.29017595056104745456D-20, 0.13639417919073099464D-18,
     2         0.23826220476859635824D-17,-0.90642907957550702534D-17,
     3        -0.14943667065169001769D-14,-0.33919078305362211264D-13,
     4        -0.17023776642512729175D-12, 0.91609750938768647911D-11,
     5         0.24230957900482704055D-09, 0.17451364971382984243D-08,
     6        -0.33126119768180852711D-07,-0.86592079961391259661D-06,
     7        -0.49717367041957398581D-05, 0.76309597585908126618D-04,
     8         0.12719271366545622927D-02, 0.17063050710955562222D-02,
     9        -0.76852840844786673690D-01,-0.28387654227602353814D+00,
     A         0.92187029365045265648D+00/
C----------------------------------------------------------------------
      EX = X
      ENU = ALPHA
      IF ((NB .GT. 0) .AND. (X .GE. XMIN) .AND. (EX .LT. XLARGE)
     1       .AND. (ENU .GE. ZERO) .AND. (ENU .LT. ONE))  THEN
            XNA = AINT(ENU+HALF)
            NA = INT(XNA)
            IF (NA .EQ. 1) ENU = ENU - XNA
            IF (ENU .EQ. -HALF) THEN
                  P = SQ2BPI/SQRT(EX)
                  YA = P * SIN(EX)
                  YA1 = -P * COS(EX)
               ELSE IF (EX .LT. THREE) THEN
C----------------------------------------------------------------------
C  Use Temme's scheme for small X
C----------------------------------------------------------------------
                  B = EX * HALF
                  D = -LOG(B)
                  F = ENU * D
                  E = B**(-ENU)
                  IF (ABS(ENU) .LT. DEL) THEN
                        C = ONBPI
                     ELSE
                        C = ENU / SIN(ENU*PI)
                  END IF
C----------------------------------------------------------------------
C  Computation of sinh(f)/f
C----------------------------------------------------------------------
                  IF (ABS(F) .LT. ONE) THEN
                        X2 = F*F
                        EN = TEN9
                        S = ONE
                        DO 80 I = 1, 9
                           S = S*X2/EN/(EN-ONE)+ONE
                           EN = EN - TWO
   80                   CONTINUE
                     ELSE 
                        S = (E - ONE/E) * HALF / F
                  END IF
C----------------------------------------------------------------------
C  Computation of 1/gamma(1-a) using Chebyshev polynomials
C----------------------------------------------------------------------
                  X2 = ENU*ENU*EIGHT
                  AYE = CH(1)
                  EVEN = ZERO
                  ALFA = CH(2)
                  ODD = ZERO
                  DO 40 I = 3, 19, 2
                     EVEN = -(AYE+AYE+EVEN)
                     AYE = -EVEN*X2 - AYE + CH(I)
                     ODD = -(ALFA+ALFA+ODD)
                     ALFA = -ODD*X2 - ALFA + CH(I+1)
   40             CONTINUE
                  EVEN = (EVEN*HALF+AYE)*X2 - AYE + CH(21)
                  ODD = (ODD+ALFA)*TWO
                  GAMMA = ODD*ENU + EVEN
C----------------------------------------------------------------------
C  End of computation of 1/gamma(1-a)
C----------------------------------------------------------------------
                  G = E * GAMMA
                  E = (E + ONE/E) * HALF
                  F = TWO*C*(ODD*E+EVEN*S*D)
                  E = ENU*ENU
                  P = G*C
                  Q = ONBPI / G
                  C = ENU*PIBY2
                  IF (ABS(C) .LT. DEL) THEN
                        R = ONE
                     ELSE 
                        R = SIN(C)/C
                  END IF
                  R = PI*C*R*R
                  C = ONE
                  D = - B*B
                  H = ZERO
                  YA = F + R*Q
                  YA1 = P
                  EN = ZERO
  100             EN = EN + ONE
                  IF (ABS(G/(ONE+ABS(YA)))
     1                      + ABS(H/(ONE+ABS(YA1))) .GT. EPS) THEN
                        F = (F*EN+P+Q)/(EN*EN-E)
                        C = C * D/EN
                        P = P/(EN-ENU)
                        Q = Q/(EN+ENU)
                        G = C*(F+R*Q)
                        H = C*P - EN*G
                        YA = YA + G
                        YA1 = YA1+H
                        GO TO 100
                  END IF
                  YA = -YA
                  YA1 = -YA1/B
               ELSE IF (EX .LT. THRESH) THEN
C----------------------------------------------------------------------
C  Use Temme's scheme for moderate X
C----------------------------------------------------------------------
                  C = (HALF-ENU)*(HALF+ENU)
                  B = EX + EX
                  E = (EX*ONBPI*COS(ENU*PI)/EPS)
                  E = E*E
                  P = ONE
                  Q = -EX
                  R = ONE + EX*EX
                  S = R
                  EN = TWO
  200             IF (R*EN*EN .LT. E) THEN
                        EN1 = EN+ONE
                        D = (EN-ONE+C/EN)/S
                        P = (EN+EN-P*D)/EN1
                        Q = (-B+Q*D)/EN1
                        S = P*P + Q*Q
                        R = R*S
                        EN = EN1
                        GO TO 200
                  END IF
                  F = P/S
                  P = F
                  G = -Q/S
                  Q = G
  220             EN = EN - ONE  
                  IF (EN .GT. ZERO) THEN
                        R = EN1*(TWO-P)-TWO
                        S = B + EN1*Q
                        D = (EN-ONE+C/EN)/(R*R+S*S)
                        P = D*R
                        Q = D*S
                        E = F + ONE
                        F = P*E - G*Q
                        G = Q*E + P*G
                        EN1 = EN
                        GO TO 220
                  END IF
                  F = ONE + F
                  D = F*F + G*G
                  PA = F/D
                  QA = -G/D
                  D = ENU + HALF -P
                  Q = Q + EX
                  PA1 = (PA*Q-QA*D)/EX
                  QA1 = (QA*Q+PA*D)/EX
                  B = EX - PIBY2*(ENU+HALF)
                  C = COS(B)
                  S = SIN(B)
                  D = SQ2BPI/SQRT(EX)
                  YA = D*(PA*S+QA*C)
                  YA1 = D*(QA1*S-PA1*C)
               ELSE
C----------------------------------------------------------------------
C  Use Campbell's asymptotic scheme.
C----------------------------------------------------------------------
                  NA = 0
                  D1 = AINT(EX/FIVPI)
                  I = INT(D1)
                  DMU = ((EX-ONE5*D1)-D1*PIM5)-(ALPHA+HALF)*PIBY2
                  IF (I-2*(I/2) .EQ. 0) THEN
                        COSMU = COS(DMU)
                        SINMU = SIN(DMU)
                     ELSE
                        COSMU = -COS(DMU)
                        SINMU = -SIN(DMU)
                  END IF
                  DDIV = EIGHT * EX
                  DMU = ALPHA
                  DEN = SQRT(EX)
                  DO 350 K = 1, 2
                     P = COSMU
                     COSMU = SINMU
                     SINMU = -P
                     D1 = (TWO*DMU-ONE)*(TWO*DMU+ONE)
                     D2 = ZERO
                     DIV = DDIV
                     P = ZERO
                     Q = ZERO
                     Q0 = D1/DIV
                     TERM = Q0
                     DO 310 I = 2, 20
                        D2 = D2 + EIGHT
                        D1 = D1 - D2
                        DIV = DIV + DDIV
                        TERM = -TERM*D1/DIV
                        P = P + TERM
                        D2 = D2 + EIGHT
                        D1 = D1 - D2
                        DIV = DIV + DDIV
                        TERM = TERM*D1/DIV
                        Q = Q + TERM
                        IF (ABS(TERM) .LE. EPS) GO TO 320
  310                CONTINUE
  320                P = P + ONE
                     Q = Q + Q0
                     IF (K .EQ. 1) THEN
                           YA = SQ2BPI * (P*COSMU-Q*SINMU) / DEN
                        ELSE
                           YA1 = SQ2BPI * (P*COSMU-Q*SINMU) / DEN
                     END IF
                     DMU = DMU + ONE
  350             CONTINUE
            END IF
            IF (NA .EQ. 1) THEN
               H = TWO*(ENU+ONE)/EX
               IF (H .GT. ONE) THEN
                  IF (ABS(YA1) .GT. XINF/H) THEN
                     H = ZERO
                     YA = ZERO
                  END IF
               END IF
               H = H*YA1 - YA
               YA = YA1
               YA1 = H
            END IF
C----------------------------------------------------------------------
C  Now have first one or two Y's
C----------------------------------------------------------------------
            BY(1) = YA
            BY(2) = YA1
            IF (YA1 .EQ. ZERO) THEN
                  NCALC = 1
               ELSE
                  AYE = ONE + ALPHA
                  TWOBYX = TWO/EX
                  NCALC = 2
                  DO 400 I = 3, NB
                     IF (TWOBYX .LT. ONE) THEN
                           IF (ABS(BY(I-1))*TWOBYX .GE. XINF/AYE)
     1                                                     GO TO 450
                        ELSE
                           IF (ABS(BY(I-1)) .GE. XINF/AYE/TWOBYX )
     1                                                     GO TO 450
                     END IF
                     BY(I) = TWOBYX*AYE*BY(I-1) - BY(I-2) 
                     AYE = AYE + ONE
                     NCALC = NCALC + 1
  400             CONTINUE
            END IF
  450       DO 460 I = NCALC+1, NB
               BY(I) = ZERO
  460       CONTINUE
         ELSE
            BY(1) = ZERO
            NCALC = MIN(NB,0) - 1
      END IF
      RETURN
C---------- Last line of RYBESL ----------
      END
C
      SUBROUTINE DIBESL(X,ALPHA,NB,IZE,B,NCALC)
C-------------------------------------------------------------------
C
C  This routine calculates Bessel functions I SUB(N+ALPHA) (X)
C  for non-negative argument X, and non-negative order N+ALPHA,
C  with or without exponential scaling.
C
C
C Explanation of variables in the calling sequence
C
C X     - Working precision non-negative real argument for which
C         I's or exponentially scaled I's (I*EXP(-X))
C         are to be calculated.  If I's are to be calculated,
C         X must be less than EXPARG (see below).
C ALPHA - Working precision fractional part of order for which
C         I's or exponentially scaled I's (I*EXP(-X)) are
C         to be calculated.  0 .LE. ALPHA .LT. 1.0.
C NB    - Integer number of functions to be calculated, NB .GT. 0.
C         The first function calculated is of order ALPHA, and the 
C         last is of order (NB - 1 + ALPHA).
C IZE   - Integer type.  IZE = 1 if unscaled I's are to calculated,
C         and 2 if exponentially scaled I's are to be calculated.
C B     - Working precision output vector of length NB.  If the routine
C         terminates normally (NCALC=NB), the vector B contains the 
C         functions I(ALPHA,X) through I(NB-1+ALPHA,X), or the
C         corresponding exponentially scaled functions.
C NCALC - Integer output variable indicating possible errors.
C         Before using the vector B, the user should check that
C         NCALC=NB, i.e., all orders have been calculated to
C         the desired accuracy.  See error returns below.
C 
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta   = Radix for the floating-point system
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C   it     = Number of bits in the mantissa of a working precision
C            variable
C   NSIG   = Decimal significance desired.  Should be set to
C            INT(LOG10(2)*it+1).  Setting NSIG lower will result
C            in decreased accuracy while setting NSIG higher will
C            increase CPU time without increasing accuracy.  The
C            truncation error is limited to a relative error of
C            T=.5*10**(-NSIG).
C   ENTEN  = 10.0 ** K, where K is the largest integer such that
C            ENTEN is machine-representable in working precision
C   ENSIG  = 10.0 ** NSIG
C   RTNSIG = 10.0 ** (-K) for the smallest integer K such that
C            K .GE. NSIG/4
C   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
C   XLARGE = Upper limit on the magnitude of X when IZE=2.  Bear
C            in mind that if ABS(X)=N, then at least N iterations
C            of the backward recursion will be executed.  The value
C            of 10.0 ** 4 is used on every machine.
C   EXPARG = Largest working precision argument that the library
C            EXP routine can handle and upper limit on the
C            magnitude of X when IZE=1; approximately 
C            LOG(beta**maxexp)
C
C
C     Approximate values for some important machines are:
C
C                        beta       minexp      maxexp       it
C
C  CRAY-1        (S.P.)    2        -8193        8191        48
C  Cyber 180/855
C    under NOS   (S.P.)    2         -975        1070        48
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)    2         -126         128        24
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)    2        -1022        1024        53
C  IBM 3033      (D.P.)   16          -65          63        14
C  VAX           (S.P.)    2         -128         127        24
C  VAX D-Format  (D.P.)    2         -128         127        56
C  VAX G-Format  (D.P.)    2        -1024        1023        53
C
C
C                        NSIG       ENTEN       ENSIG      RTNSIG
C
C CRAY-1        (S.P.)    15       1.0E+2465   1.0E+15     1.0E-4
C Cyber 180/855
C   under NOS   (S.P.)    15       1.0E+322    1.0E+15     1.0E-4
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)    16       1.0D+308    1.0D+16     1.0D-4
C IBM 3033      (D.P.)     5       1.0D+75     1.0D+5      1.0D-2
C VAX           (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
C VAX D-Format  (D.P.)    17       1.0D+38     1.0D+17     1.0D-5
C VAX G-Format  (D.P.)    16       1.0D+307    1.0D+16     1.0D-4
C
C
C                         ENMTEN      XLARGE   EXPARG 
C
C CRAY-1        (S.P.)   1.84E-2466   1.0E+4    5677 
C Cyber 180/855
C   under NOS   (S.P.)   1.25E-293    1.0E+4     741
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)   4.70E-38     1.0E+4      88  
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)   8.90D-308    1.0D+4     709
C IBM 3033      (D.P.)   2.16D-78     1.0D+4     174
C VAX           (S.P.)   1.17E-38     1.0E+4      88
C VAX D-Format  (D.P.)   1.17D-38     1.0D+4      88
C VAX G-Format  (D.P.)   2.22D-308    1.0D+4     709
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  In case of an error,  NCALC .NE. NB, and not all I's are
C  calculated to the desired accuracy.
C
C  NCALC .LT. 0:  An argument is out of range. For example,
C     NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE. EXPARG.
C     In this case, the B-vector is not calculated, and NCALC is
C     set to MIN0(NB,0)-1 so that NCALC .NE. NB.
C
C  NB .GT. NCALC .GT. 0: Not all requested function values could
C     be calculated accurately.  This usually occurs because NB is
C     much larger than ABS(X).  In this case, B(N) is calculated
C     to the desired accuracy for N .LE. NCALC, but precision
C     is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
C     for N .GT. NCALC (because it is too small to be represented),
C     and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
C     significant figures of B(N) can be trusted.
C
C
C Intrinsic functions required are:
C
C     DBLE, EXP, DGAMMA, GAMMA, INT, MAX, MIN, REAL, SQRT
C
C
C Acknowledgement
C
C  This program is based on a program written by David J.
C  Sookne (2) that computes values of the Bessel functions J or
C  I of real argument and integer order.  Modifications include
C  the restriction of the computation to the I Bessel function
C  of non-negative real argument, the extension of the computation
C  to arbitrary positive order, the inclusion of optional
C  exponential scaling, and the elimination of most underflow.
C  An earlier version was published in (3).
C
C References: "A Note on Backward Recurrence Algorithms," Olver,
C              F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
C              pp 941-947.
C
C             "Bessel Functions of Real Argument and Integer Order,"
C              Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
C              125-132.
C
C             "ALGORITHM 597, Sequence of Modified Bessel Functions
C              of the First Kind," Cody, W. J., Trans. Math. Soft.,
C              1983, pp. 242-245.
C
C  Latest modification: May 30, 1989
C
C  Modified by: W. J. Cody and L. Stoltz
C               Applied Mathematics Division
C               Argonne National Laboratory
C               Argonne, IL  60439
C
C-------------------------------------------------------------------
      INTEGER IZE,K,L,MAGX,N,NB,NBMX,NCALC,NEND,NSIG,NSTART
CS    REAL              GAMMA,
      DOUBLE PRECISION DGAMMA,             
     1 ALPHA,B,CONST,CONV,EM,EMPAL,EMP2AL,EN,ENMTEN,ENSIG,
     2 ENTEN,EXPARG,FUNC,HALF,HALFX,ONE,P,PLAST,POLD,PSAVE,PSAVEL,
     3 RTNSIG,SUM,TEMPA,TEMPB,TEMPC,TEST,TOVER,TWO,X,XLARGE,ZERO
      DIMENSION B(NB)
C-------------------------------------------------------------------
C  Mathematical constants
C-------------------------------------------------------------------
CS    DATA ONE,TWO,ZERO,HALF,CONST/1.0E0,2.0E0,0.0E0,0.5E0,1.585E0/
      DATA ONE,TWO,ZERO,HALF,CONST/1.0D0,2.0D0,0.0D0,0.5D0,1.585D0/
C-------------------------------------------------------------------
C  Machine-dependent parameters
C-------------------------------------------------------------------
CS    DATA NSIG,XLARGE,EXPARG /8,1.0E4,88.0E0/
CS    DATA ENTEN,ENSIG,RTNSIG/1.0E38,1.0E8,1.0E-2/
CS    DATA ENMTEN/4.7E-38/
      DATA NSIG,XLARGE,EXPARG /16,1.0D4,709.0D0/
      DATA ENTEN,ENSIG,RTNSIG/1.0D308,1.0D16,1.0D-4/
      DATA ENMTEN/8.9D-308/
C-------------------------------------------------------------------
C  Statement functions for conversion
C-------------------------------------------------------------------
CS    CONV(N) = REAL(N)
CS    FUNC(X) = GAMMA(X)
      CONV(N) = DBLE(N)
      FUNC(X) = DGAMMA(X)
C-------------------------------------------------------------------
C Check for X, NB, OR IZE out of range.
C-------------------------------------------------------------------
      IF ((NB.GT.0) .AND. (X .GE. ZERO) .AND.
     1    (ALPHA .GE. ZERO) .AND. (ALPHA .LT. ONE) .AND.
     2    (((IZE .EQ. 1) .AND. (X .LE. EXPARG)) .OR.
     3     ((IZE .EQ. 2) .AND. (X .LE. XLARGE)))) THEN
C-------------------------------------------------------------------
C Use 2-term ascending series for small X
C-------------------------------------------------------------------
            NCALC = NB
            MAGX = INT(X)
            IF (X .GE. RTNSIG) THEN
C-------------------------------------------------------------------
C Initialize the forward sweep, the P-sequence of Olver
C-------------------------------------------------------------------
                  NBMX = NB-MAGX
                  N = MAGX+1
                  EN = CONV(N+N) + (ALPHA+ALPHA)
                  PLAST = ONE
                  P = EN / X
C-------------------------------------------------------------------
C Calculate general significance test
C-------------------------------------------------------------------
                  TEST = ENSIG + ENSIG
                  IF (2*MAGX .GT. 5*NSIG) THEN
                        TEST = SQRT(TEST*P)
                     ELSE
                        TEST = TEST / CONST**MAGX
                  END IF
                  IF (NBMX .GE. 3) THEN
C-------------------------------------------------------------------
C Calculate P-sequence until N = NB-1.  Check for possible overflow.
C-------------------------------------------------------------------
                     TOVER = ENTEN / ENSIG
                     NSTART = MAGX+2
                     NEND = NB - 1
                     DO 100 K = NSTART, NEND
                        N = K
                        EN = EN + TWO
                        POLD = PLAST
                        PLAST = P
                        P = EN * PLAST/X + POLD
                        IF (P .GT. TOVER) THEN
C-------------------------------------------------------------------
C To avoid overflow, divide P-sequence by TOVER.  Calculate
C P-sequence until ABS(P) .GT. 1.
C-------------------------------------------------------------------
                           TOVER = ENTEN
                           P = P / TOVER
                           PLAST = PLAST / TOVER
                           PSAVE = P
                           PSAVEL = PLAST
                           NSTART = N + 1
   60                      N = N + 1
                              EN = EN + TWO
                              POLD = PLAST
                              PLAST = P
                              P = EN * PLAST/X + POLD
                           IF (P .LE. ONE) GO TO 60
                           TEMPB = EN / X
C-------------------------------------------------------------------
C Calculate backward test, and find NCALC, the highest N
C such that the test is passed.
C-------------------------------------------------------------------
                           TEST = POLD*PLAST / ENSIG
                           TEST = TEST*(HALF-HALF/(TEMPB*TEMPB))
                           P = PLAST * TOVER
                           N = N - 1
                           EN = EN - TWO
                           NEND = MIN0(NB,N)
                           DO 80 L = NSTART, NEND
                              NCALC = L
                              POLD = PSAVEL
                              PSAVEL = PSAVE
                              PSAVE = EN * PSAVEL/X + POLD
                              IF (PSAVE*PSAVEL .GT. TEST) GO TO 90
   80                      CONTINUE
                           NCALC = NEND + 1
   90                      NCALC = NCALC - 1
                           GO TO 120
                        END IF
  100                CONTINUE
                     N = NEND
                     EN = CONV(N+N) + (ALPHA+ALPHA)
C-------------------------------------------------------------------
C Calculate special significance test for NBMX .GT. 2.
C-------------------------------------------------------------------
                     TEST = MAX(TEST,SQRT(PLAST*ENSIG)*SQRT(P+P))
                  END IF
C-------------------------------------------------------------------
C Calculate P-sequence until significance test passed.
C-------------------------------------------------------------------
  110             N = N + 1
                     EN = EN + TWO
                     POLD = PLAST
                     PLAST = P
                     P = EN * PLAST/X + POLD
                  IF (P .LT. TEST) GO TO 110
C-------------------------------------------------------------------
C Initialize the backward recursion and the normalization sum.
C-------------------------------------------------------------------
  120             N = N + 1
                  EN = EN + TWO
                  TEMPB = ZERO
                  TEMPA = ONE / P
                  EM = CONV(N) - ONE
                  EMPAL = EM + ALPHA
                  EMP2AL = (EM - ONE) + (ALPHA + ALPHA)
                  SUM = TEMPA * EMPAL * EMP2AL / EM
                  NEND = N - NB
                  IF (NEND .LT. 0) THEN
C-------------------------------------------------------------------
C N .LT. NB, so store B(N) and set higher orders to zero.
C-------------------------------------------------------------------
                        B(N) = TEMPA
                        NEND = -NEND
                        DO 130 L = 1, NEND
  130                      B(N+L) = ZERO
                     ELSE
                        IF (NEND .GT. 0) THEN
C-------------------------------------------------------------------
C Recur backward via difference equation, calculating (but
C not storing) B(N), until N = NB.
C-------------------------------------------------------------------
                           DO 140 L = 1, NEND
                              N = N - 1
                              EN = EN - TWO
                              TEMPC = TEMPB
                              TEMPB = TEMPA
                              TEMPA = (EN*TEMPB) / X + TEMPC
                              EM = EM - ONE
                              EMP2AL = EMP2AL - ONE
                              IF (N .EQ. 1) GO TO 150
                              IF (N .EQ. 2) EMP2AL = ONE
                              EMPAL = EMPAL - ONE
                              SUM = (SUM + TEMPA*EMPAL) * EMP2AL / EM
  140                      CONTINUE
                        END IF
C-------------------------------------------------------------------
C Store B(NB)
C-------------------------------------------------------------------
  150                   B(N) = TEMPA
                        IF (NB .LE. 1) THEN
                           SUM = (SUM + SUM) + TEMPA
                           GO TO 230
                        END IF
C-------------------------------------------------------------------
C Calculate and Store B(NB-1)
C-------------------------------------------------------------------
                        N = N - 1
                        EN = EN - TWO
                        B(N)  = (EN*TEMPA) / X + TEMPB
                        IF (N .EQ. 1) GO TO 220
                        EM = EM - ONE
                        EMP2AL = EMP2AL - ONE
                        IF (N .EQ. 2) EMP2AL = ONE
                        EMPAL = EMPAL - ONE
                        SUM = (SUM + B(N)*EMPAL) * EMP2AL / EM
                  END IF
                  NEND = N - 2
                  IF (NEND .GT. 0) THEN
C-------------------------------------------------------------------
C Calculate via difference equation and store B(N), until N = 2.
C-------------------------------------------------------------------
                     DO 200 L = 1, NEND
                        N = N - 1
                        EN = EN - TWO
                        B(N) = (EN*B(N+1)) / X +B(N+2)
                        EM = EM - ONE
                        EMP2AL = EMP2AL - ONE
                        IF (N .EQ. 2) EMP2AL = ONE
                        EMPAL = EMPAL - ONE
                        SUM = (SUM + B(N)*EMPAL) * EMP2AL / EM
  200                CONTINUE
                  END IF
C-------------------------------------------------------------------
C Calculate B(1)
C-------------------------------------------------------------------
                  B(1) = TWO*EMPAL*B(2) / X + B(3)
  220             SUM = (SUM + SUM) + B(1)
C-------------------------------------------------------------------
C Normalize.  Divide all B(N) by sum.
C-------------------------------------------------------------------
  230             IF (ALPHA .NE. ZERO)
     1               SUM = SUM * FUNC(ONE+ALPHA) * (X*HALF)**(-ALPHA)
                  IF (IZE .EQ. 1) SUM = SUM * EXP(-X)
                  TEMPA = ENMTEN
                  IF (SUM .GT. ONE) TEMPA = TEMPA * SUM
                  DO 260 N = 1, NB
                     IF (B(N) .LT. TEMPA) B(N) = ZERO
                     B(N) = B(N) / SUM
  260             CONTINUE
                  RETURN
C-------------------------------------------------------------------
C Two-term ascending series for small X.
C-------------------------------------------------------------------
               ELSE
                  TEMPA = ONE
                  EMPAL = ONE + ALPHA
                  HALFX = ZERO
                  IF (X .GT. ENMTEN) HALFX = HALF * X
                  IF (ALPHA .NE. ZERO) TEMPA = HALFX**ALPHA /FUNC(EMPAL)
                  IF (IZE .EQ. 2) TEMPA = TEMPA * EXP(-X)
                  TEMPB = ZERO
                  IF ((X+ONE) .GT. ONE) TEMPB = HALFX * HALFX
                  B(1) = TEMPA + TEMPA*TEMPB / EMPAL
                  IF ((X .NE. ZERO) .AND. (B(1) .EQ. ZERO)) NCALC = 0
                  IF (NB .GT. 1) THEN
                     IF (X .EQ. ZERO) THEN
                           DO 310 N = 2, NB
                              B(N) = ZERO
  310                      CONTINUE
                        ELSE
C-------------------------------------------------------------------
C Calculate higher-order functions.
C-------------------------------------------------------------------
                           TEMPC = HALFX
                           TOVER = (ENMTEN + ENMTEN) / X
                           IF (TEMPB .NE. ZERO) TOVER = ENMTEN / TEMPB
                           DO 340 N = 2, NB
                              TEMPA = TEMPA / EMPAL
                              EMPAL = EMPAL + ONE
                              TEMPA = TEMPA * TEMPC
                              IF (TEMPA .LE. TOVER*EMPAL) TEMPA = ZERO
                              B(N) = TEMPA + TEMPA*TEMPB / EMPAL
                              IF ((B(N) .EQ. ZERO) .AND. (NCALC .GT. N))
     1                             NCALC = N-1
  340                      CONTINUE
                     END IF
                  END IF
            END IF
         ELSE
            NCALC = MIN0(NB,0)-1
      END IF
      RETURN
C---------- Last line of RIBESL ----------
      END
C
      SUBROUTINE DKBESL(X,ALPHA,NB,IZE,BK,NCALC)
C-------------------------------------------------------------------
C
C  This FORTRAN 77 routine calculates modified Bessel functions
C  of the second kind, K SUB(N+ALPHA) (X), for non-negative
C  argument X, and non-negative order N+ALPHA, with or without
C  exponential scaling.
C
C  Explanation of variables in the calling sequence
C
C  Description of output values ..
C
C X     - Working precision non-negative real argument for which
C         K's or exponentially scaled K's (K*EXP(X))
C         are to be calculated.  If K's are to be calculated,
C         X must not be greater than XMAX (see below).
C ALPHA - Working precision fractional part of order for which 
C         K's or exponentially scaled K's (K*EXP(X)) are
C         to be calculated.  0 .LE. ALPHA .LT. 1.0.
C NB    - Integer number of functions to be calculated, NB .GT. 0.
C         The first function calculated is of order ALPHA, and the 
C         last is of order (NB - 1 + ALPHA).
C IZE   - Integer type.  IZE = 1 if unscaled K's are to be calculated,
C         and 2 if exponentially scaled K's are to be calculated.
C BK    - Working precision output vector of length NB.  If the
C         routine terminates normally (NCALC=NB), the vector BK
C         contains the functions K(ALPHA,X), ... , K(NB-1+ALPHA,X),
C         or the corresponding exponentially scaled functions.
C         If (0 .LT. NCALC .LT. NB), BK(I) contains correct function
C         values for I .LE. NCALC, and contains the ratios
C         K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
C NCALC - Integer output variable indicating possible errors.
C         Before using the vector BK, the user should check that 
C         NCALC=NB, i.e., all orders have been calculated to
C         the desired accuracy.  See error returns below.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta   = Radix for the floating-point system
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C   EPS    = The smallest positive floating-point number such that 
C            1.0+EPS .GT. 1.0
C   XMAX   = Upper limit on the magnitude of X when IZE=1;  Solution 
C            to equation:
C               W(X) * (1-1/8X+9/128X**2) = beta**minexp
C            where  W(X) = EXP(-X)*SQRT(PI/2X)
C   SQXMIN = Square root of beta**minexp
C   XINF   = Largest positive machine number; approximately
C            beta**maxexp
C   XMIN   = Smallest positive machine number; approximately
C            beta**minexp
C
C
C     Approximate values for some important machines are:
C
C                          beta       minexp      maxexp      EPS
C
C  CRAY-1        (S.P.)      2        -8193        8191    7.11E-15
C  Cyber 180/185 
C    under NOS   (S.P.)      2         -975        1070    3.55E-15
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)      2         -126         128    1.19E-7
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)      2        -1022        1024    2.22D-16
C  IBM 3033      (D.P.)     16          -65          63    2.22D-16
C  VAX           (S.P.)      2         -128         127    5.96E-8
C  VAX D-Format  (D.P.)      2         -128         127    1.39D-17
C  VAX G-Format  (D.P.)      2        -1024        1023    1.11D-16
C
C
C                         SQXMIN       XINF        XMIN      XMAX
C
C CRAY-1        (S.P.)  6.77E-1234  5.45E+2465  4.59E-2467 5674.858
C Cyber 180/855
C   under NOS   (S.P.)  1.77E-147   1.26E+322   3.14E-294   672.788
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)  1.08E-19    3.40E+38    1.18E-38     85.337
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)  1.49D-154   1.79D+308   2.23D-308   705.342
C IBM 3033      (D.P.)  7.35D-40    7.23D+75    5.40D-79    177.852
C VAX           (S.P.)  5.42E-20    1.70E+38    2.94E-39     86.715
C VAX D-Format  (D.P.)  5.42D-20    1.70D+38    2.94D-39     86.715
C VAX G-Format  (D.P.)  7.46D-155   8.98D+307   5.57D-309   706.728
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  In case of an error, NCALC .NE. NB, and not all K's are
C  calculated to the desired accuracy.
C
C  NCALC .LT. -1:  An argument is out of range. For example,
C       NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
C       XMAX.  In this case, the B-vector is not calculated,
C       and NCALC is set to MIN0(NB,0)-2  so that NCALC .NE. NB.
C  NCALC = -1:  Either  K(ALPHA,X) .GE. XINF  or
C       K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) .GE. XINF.  In this case,
C       the B-vector is not calculated.  Note that again 
C       NCALC .NE. NB.
C
C  0 .LT. NCALC .LT. NB: Not all requested function values could
C       be calculated accurately.  BK(I) contains correct function
C       values for I .LE. NCALC, and contains the ratios
C       K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, EXP, INT, LOG, MAX, MIN, SINH, SQRT
C
C
C Acknowledgement
C
C  This program is based on a program written by J. B. Campbell
C  (2) that computes values of the Bessel functions K of real
C  argument and real order.  Modifications include the addition
C  of non-scaled functions, parameterization of machine
C  dependencies, and the use of more accurate approximations
C  for SINH and SIN.
C
C References: "On Temme's Algorithm for the Modified Bessel
C              Functions of the Third Kind," Campbell, J. B.,
C              TOMS 6(4), Dec. 1980, pp. 581-586.
C
C             "A FORTRAN IV Subroutine for the Modified Bessel
C              Functions of the Third Kind of Real Order and Real
C              Argument," Campbell, J. B., Report NRC/ERB-925,
C              National Research Council, Canada.
C
C  Latest modification: May 30, 1989
C
C  Modified by: W. J. Cody and L. Stoltz
C               Applied Mathematics Division
C               Argonne National Laboratory
C               Argonne, IL  60439
C
C-------------------------------------------------------------------
      INTEGER I,IEND,ITEMP,IZE,J,K,M,MPLUS1,NB,NCALC
CS    REAL
      DOUBLE PRECISION  
     1    A,ALPHA,BLPHA,BK,BK1,BK2,C,D,DM,D1,D2,D3,ENU,EPS,ESTF,ESTM,
     2    EX,FOUR,F0,F1,F2,HALF,ONE,P,P0,Q,Q0,R,RATIO,S,SQXMIN,T,TINYX,
     3    TWO,TWONU,TWOX,T1,T2,WMINF,X,XINF,XMAX,XMIN,X2BY4,ZERO
      DIMENSION BK(1),P(8),Q(7),R(5),S(4),T(6),ESTM(6),ESTF(7)
C---------------------------------------------------------------------
C  Mathematical constants
C    A = LOG(2.D0) - Euler's constant
C    D = SQRT(2.D0/PI)
C---------------------------------------------------------------------
CS    DATA HALF,ONE,TWO,ZERO/0.5E0,1.0E0,2.0E0,0.0E0/
CS    DATA FOUR,TINYX/4.0E0,1.0E-10/
CS    DATA A/ 0.11593151565841244881E0/,D/0.797884560802865364E0/
      DATA HALF,ONE,TWO,ZERO/0.5D0,1.0D0,2.0D0,0.0D0/
      DATA FOUR,TINYX/4.0D0,1.0D-10/
      DATA A/ 0.11593151565841244881D0/,D/0.797884560802865364D0/
C---------------------------------------------------------------------
C  Machine dependent parameters
C---------------------------------------------------------------------
CS    DATA EPS/1.19E-7/,SQXMIN/1.08E-19/,XINF/3.40E+38/
CS    DATA XMIN/1.18E-38/,XMAX/85.337E0/
      DATA EPS/2.22D-16/,SQXMIN/1.49D-154/,XINF/1.79D+308/
      DATA XMIN/2.23D-308/,XMAX/705.342D0/
C---------------------------------------------------------------------
C  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA
C                                         + Euler's constant
C         Coefficients converted from hex to decimal and modified
C         by W. J. Cody, 2/26/82
C  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA)
C  T    - Approximation for SINH(Y)/Y
C---------------------------------------------------------------------
CS    DATA P/ 0.805629875690432845E00,    0.204045500205365151E02,
CS   1        0.157705605106676174E03,    0.536671116469207504E03,
CS   2        0.900382759291288778E03,    0.730923886650660393E03,
CS   3        0.229299301509425145E03,    0.822467033424113231E00/
CS    DATA Q/ 0.294601986247850434E02,    0.277577868510221208E03,
CS   1        0.120670325591027438E04,    0.276291444159791519E04,
CS   2        0.344374050506564618E04,    0.221063190113378647E04,
CS   3        0.572267338359892221E03/
CS    DATA R/-0.48672575865218401848E+0,  0.13079485869097804016E+2,
CS   1       -0.10196490580880537526E+3,  0.34765409106507813131E+3,
CS   2        0.34958981245219347820E-3/
CS    DATA S/-0.25579105509976461286E+2,  0.21257260432226544008E+3,
CS   1       -0.61069018684944109624E+3,  0.42269668805777760407E+3/
CS    DATA T/ 0.16125990452916363814E-9, 0.25051878502858255354E-7,
CS   1        0.27557319615147964774E-5, 0.19841269840928373686E-3,
CS   2        0.83333333333334751799E-2, 0.16666666666666666446E+0/
CS    DATA ESTM/5.20583E1, 5.7607E0, 2.7782E0, 1.44303E1, 1.853004E2,
CS   1          9.3715E0/
CS    DATA ESTF/4.18341E1, 7.1075E0, 6.4306E0, 4.25110E1, 1.35633E0,
CS   1          8.45096E1, 2.0E1/
      DATA P/ 0.805629875690432845D00,    0.204045500205365151D02,
     1        0.157705605106676174D03,    0.536671116469207504D03,
     2        0.900382759291288778D03,    0.730923886650660393D03,
     3        0.229299301509425145D03,    0.822467033424113231D00/
      DATA Q/ 0.294601986247850434D02,    0.277577868510221208D03,
     1        0.120670325591027438D04,    0.276291444159791519D04,
     2        0.344374050506564618D04,    0.221063190113378647D04,
     3        0.572267338359892221D03/
      DATA R/-0.48672575865218401848D+0,  0.13079485869097804016D+2,
     1       -0.10196490580880537526D+3,  0.34765409106507813131D+3,
     2        0.34958981245219347820D-3/
      DATA S/-0.25579105509976461286D+2,  0.21257260432226544008D+3,
     1       -0.61069018684944109624D+3,  0.42269668805777760407D+3/
      DATA T/ 0.16125990452916363814D-9, 0.25051878502858255354D-7,
     1        0.27557319615147964774D-5, 0.19841269840928373686D-3,
     2        0.83333333333334751799D-2, 0.16666666666666666446D+0/
      DATA ESTM/5.20583D1, 5.7607D0, 2.7782D0, 1.44303D1, 1.853004D2,
     1          9.3715D0/
      DATA ESTF/4.18341D1, 7.1075D0, 6.4306D0, 4.25110D1, 1.35633D0,
     1          8.45096D1, 2.0D1/
C---------------------------------------------------------------------
      EX = X
      ENU = ALPHA
      NCALC = MIN(NB,0)-2
      IF ((NB .GT. 0) .AND. ((ENU .GE. ZERO) .AND. (ENU .LT. ONE))
     1     .AND. ((IZE .GE. 1) .AND. (IZE .LE. 2)) .AND.
     2     ((IZE .NE. 1) .OR. (EX .LE. XMAX)) .AND.
     3     (EX .GT. ZERO))  THEN
            K = 0
            IF (ENU .LT. SQXMIN) ENU = ZERO
            IF (ENU .GT. HALF) THEN
                  K = 1
                  ENU = ENU - ONE
            END IF
            TWONU = ENU+ENU
            IEND = NB+K-1
            C = ENU*ENU
            D3 = -C
            IF (EX .LE. ONE) THEN
C---------------------------------------------------------------------
C  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA
C                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA
C---------------------------------------------------------------------
                  D1 = ZERO
                  D2 = P(1)
                  T1 = ONE
                  T2 = Q(1)
                  DO 10 I = 2,7,2
                     D1 = C*D1+P(I)
                     D2 = C*D2+P(I+1)
                     T1 = C*T1+Q(I)
                     T2 = C*T2+Q(I+1)
   10             CONTINUE
                  D1 = ENU*D1
                  T1 = ENU*T1
                  F1 = LOG(EX)
                  F0 = A+ENU*(P(8)-ENU*(D1+D2)/(T1+T2))-F1
                  Q0 = EXP(-ENU*(A-ENU*(P(8)+ENU*(D1-D2)/(T1-T2))-F1))
                  F1 = ENU*F0
                  P0 = EXP(F1)
C---------------------------------------------------------------------
C  Calculation of F0 = 
C---------------------------------------------------------------------
                  D1 = R(5)
                  T1 = ONE
                  DO 20 I = 1,4
                     D1 = C*D1+R(I)
                     T1 = C*T1+S(I)
   20             CONTINUE
                  IF (ABS(F1) .LE. HALF) THEN
                        F1 = F1*F1
                        D2 = ZERO
                        DO 30 I = 1,6
                           D2 = F1*D2+T(I)
   30                   CONTINUE
                        D2 = F0+F0*F1*D2
                     ELSE
                        D2 = SINH(F1)/ENU
                  END IF
                  F0 = D2-ENU*D1/(T1*P0)
                  IF (EX .LE. TINYX) THEN
C--------------------------------------------------------------------
C  X.LE.1.0E-10
C  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X)
C--------------------------------------------------------------------
                        BK(1) = F0+EX*F0
                        IF (IZE .EQ. 1) BK(1) = BK(1)-EX*BK(1)
                        RATIO = P0/F0
                        C = EX*XINF
                        IF (K .NE. 0) THEN
C--------------------------------------------------------------------
C  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X),
C  ALPHA .GE. 1/2
C--------------------------------------------------------------------
                              NCALC = -1
                              IF (BK(1) .GE. C/RATIO) GO TO 500
                              BK(1) = RATIO*BK(1)/EX
                              TWONU = TWONU+TWO
                              RATIO = TWONU
                        END IF
                        NCALC = 1
                        IF (NB .EQ. 1) GO TO 500
C--------------------------------------------------------------------
C  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1, 2, ... , NB-1
C--------------------------------------------------------------------
                        NCALC = -1
                        DO 80 I = 2,NB
                           IF (RATIO .GE. C) GO TO 500
                           BK(I) = RATIO/EX
                           TWONU = TWONU+TWO
                           RATIO = TWONU
   80                   CONTINUE
                        NCALC = 1
                        GO TO 420
                     ELSE
C--------------------------------------------------------------------
C  1.0E-10 .LT. X .LE. 1.0
C--------------------------------------------------------------------
                        C = ONE
                        X2BY4 = EX*EX/FOUR
                        P0 = HALF*P0
                        Q0 = HALF*Q0
                        D1 = -ONE
                        D2 = ZERO
                        BK1 = ZERO
                        BK2 = ZERO
                        F1 = F0
                        F2 = P0
  100                   D1 = D1+TWO
                        D2 = D2+ONE
                        D3 = D1+D3
                        C = X2BY4*C/D2
                        F0 = (D2*F0+P0+Q0)/D3
                        P0 = P0/(D2-ENU)
                        Q0 = Q0/(D2+ENU)
                        T1 = C*F0
                        T2 = C*(P0-D2*F0)
                        BK1 = BK1+T1
                        BK2 = BK2+T2
                        IF ((ABS(T1/(F1+BK1)) .GT. EPS) .OR.
     1                     (ABS(T2/(F2+BK2)) .GT. EPS))  GO TO 100
                        BK1 = F1+BK1
                        BK2 = TWO*(F2+BK2)/EX
                        IF (IZE .EQ. 2) THEN
                              D1 = EXP(EX)
                              BK1 = BK1*D1
                              BK2 = BK2*D1
                        END IF
                        WMINF = ESTF(1)*EX+ESTF(2)
                  END IF
               ELSE IF (EPS*EX .GT. ONE) THEN
C--------------------------------------------------------------------
C  X .GT. ONE/EPS
C--------------------------------------------------------------------
                  NCALC = NB
                  BK1 = ONE / (D*SQRT(EX))
                  DO 110 I = 1, NB
                     BK(I) = BK1
  110             CONTINUE
                  GO TO 500
               ELSE
C--------------------------------------------------------------------
C  X .GT. 1.0
C--------------------------------------------------------------------
                  TWOX = EX+EX
                  BLPHA = ZERO
                  RATIO = ZERO
                  IF (EX .LE. FOUR) THEN
C--------------------------------------------------------------------
C  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 .LE. X .LE. 4.0
C--------------------------------------------------------------------
                        D2 = AINT(ESTM(1)/EX+ESTM(2))
                        M = INT(D2)
                        D1 = D2+D2
                        D2 = D2-HALF
                        D2 = D2*D2
                        DO 120 I = 2,M
                           D1 = D1-TWO
                           D2 = D2-D1
                           RATIO = (D3+D2)/(TWOX+D1-RATIO)
  120                   CONTINUE
C--------------------------------------------------------------------
C  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
C    recurrence and K(ALPHA,X) from the wronskian
C--------------------------------------------------------------------
                        D2 = AINT(ESTM(3)*EX+ESTM(4))
                        M = INT(D2)
                        C = ABS(ENU)
                        D3 = C+C
                        D1 = D3-ONE
                        F1 = XMIN
                        F0 = (TWO*(C+D2)/EX+HALF*EX/(C+D2+ONE))*XMIN
                        DO 130 I = 3,M
                           D2 = D2-ONE
                           F2 = (D3+D2+D2)*F0
                           BLPHA = (ONE+D1/D2)*(F2+BLPHA)
                           F2 = F2/EX+F1
                           F1 = F0
                           F0 = F2
  130                   CONTINUE
                        F1 = (D3+TWO)*F0/EX+F1
                        D1 = ZERO
                        T1 = ONE
                        DO 140 I = 1,7
                           D1 = C*D1+P(I)
                           T1 = C*T1+Q(I)
  140                   CONTINUE
                        P0 = EXP(C*(A+C*(P(8)-C*D1/T1)-LOG(EX)))/EX
                        F2 = (C+HALF-RATIO)*F1/EX
                        BK1 = P0+(D3*F0-F2+F0+BLPHA)/(F2+F1+F0)*P0
                        IF (IZE .EQ. 1) BK1 = BK1*EXP(-EX)
                        WMINF = ESTF(3)*EX+ESTF(4)
                     ELSE
C--------------------------------------------------------------------
C  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
C  recurrence, for  X .GT. 4.0
C--------------------------------------------------------------------
                        DM = AINT(ESTM(5)/EX+ESTM(6))
                        M = INT(DM)
                        D2 = DM-HALF
                        D2 = D2*D2
                        D1 = DM+DM
                        DO 160 I = 2,M
                           DM = DM-ONE
                           D1 = D1-TWO
                           D2 = D2-D1
                           RATIO = (D3+D2)/(TWOX+D1-RATIO)
                           BLPHA = (RATIO+RATIO*BLPHA)/DM
  160                   CONTINUE
                        BK1 = ONE/((D+D*BLPHA)*SQRT(EX))
                        IF (IZE .EQ. 1) BK1 = BK1*EXP(-EX)
                        WMINF = ESTF(5)*(EX-ABS(EX-ESTF(7)))+ESTF(6)
                  END IF
C--------------------------------------------------------------------
C  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and
C    K(ALPHA+1,X)/K(ALPHA,X)
C--------------------------------------------------------------------
                  BK2 = BK1+BK1*(ENU+HALF-RATIO)/EX
            END IF
C--------------------------------------------------------------------
C  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1,
C  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1
C--------------------------------------------------------------------
            NCALC = NB
            BK(1) = BK1
            IF (IEND .EQ. 0) GO TO 500
            J = 2-K
            IF (J .GT. 0) BK(J) = BK2
            IF (IEND .EQ. 1) GO TO 500
            M = MIN(INT(WMINF-ENU),IEND)
            DO 190 I = 2,M
               T1 = BK1
               BK1 = BK2
               TWONU = TWONU+TWO
               IF (EX .LT. ONE) THEN
                     IF (BK1 .GE. (XINF/TWONU)*EX) GO TO 195
                     GO TO 187
                  ELSE 
                     IF (BK1/EX .GE. XINF/TWONU) GO TO 195
               END IF
  187          CONTINUE
               BK2 = TWONU/EX*BK1+T1
               ITEMP = I
               J = J+1
               IF (J .GT. 0) BK(J) = BK2
  190       CONTINUE
  195       M = ITEMP
            IF (M .EQ. IEND) GO TO 500
            RATIO = BK2/BK1
            MPLUS1 = M+1
            NCALC = -1
            DO 410 I = MPLUS1,IEND
               TWONU = TWONU+TWO
               RATIO = TWONU/EX+ONE/RATIO
               J = J+1
               IF (J .GT. 1) THEN
                     BK(J) = RATIO
                  ELSE
                     IF (BK2 .GE. XINF/RATIO) GO TO 500
                     BK2 = RATIO*BK2
               END IF
  410       CONTINUE
            NCALC = MAX(MPLUS1-K,1)
            IF (NCALC .EQ. 1) BK(1) = BK2
            IF (NB .EQ. 1) GO TO 500
  420       J = NCALC+1
            DO 430 I = J,NB
               IF (BK(NCALC) .GE. XINF/BK(I)) GO TO 500
               BK(I) = BK(NCALC)*BK(I)
               NCALC = I
  430       CONTINUE
      END IF
  500 RETURN
C---------- Last line of RKBESL ----------
      END
C
      SUBROUTINE CALJY0(ARG,RESULT,JINT)
C---------------------------------------------------------------------
C
C This packet computes zero-order Bessel functions of the first and
C   second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
C   for Y0, and |X| <= XMAX for J0.  It contains two function-type
C   subprograms,  BESJ0  and  BESY0,  and one subroutine-type
C   subprogram,  CALJY0.  The calling statements for the primary
C   entries are:
C
C           Y = BESJ0(X)
C   and
C           Y = BESY0(X),
C
C   where the entry points correspond to the functions J0(X) and Y0(X),
C   respectively.  The routine  CALJY0  is intended for internal packet
C   use only, all computations within the packet being concentrated in
C   this one routine.  The function subprograms invoke  CALJY0  with
C   the statement
C           CALL CALJY0(ARG,RESULT,JINT),
C   where the parameter usage is as follows:
C
C      Function                  Parameters for CALJY0
C       call              ARG             RESULT          JINT
C
C     BESJ0(ARG)     |ARG| .LE. XMAX       J0(ARG)          0
C     BESY0(ARG)   0 .LT. ARG .LE. XMAX    Y0(ARG)          1
C
C   The main computation uses unpublished minimax rational
C   approximations for X .LE. 8.0, and an approximation from the 
C   book  Computer Approximations  by Hart, et. al., Wiley and Sons, 
C   New York, 1968, for arguments larger than 8.0   Part of this
C   transportable packet is patterned after the machine-dependent
C   FUNPACK program BESJ0(X), but cannot match that version for
C   efficiency or accuracy.  This version uses rational functions
C   that are theoretically accurate to at least 18 significant decimal
C   digits for X <= 8, and at least 18 decimal places for X > 8.  The
C   accuracy achieved depends on the arithmetic system, the compiler,
C   the intrinsic functions, and proper selection of the machine-
C   dependent constants.
C
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   XINF   = largest positive machine number
C   XMAX   = largest acceptable argument.  The functions AINT, SIN
C            and COS must perform properly for  ABS(X) .LE. XMAX.
C            We recommend that XMAX be a small integer multiple of
C            sqrt(1/eps), where eps is the smallest positive number
C            such that  1+eps > 1. 
C   XSMALL = positive argument such that  1.0-(X/2)**2 = 1.0
C            to machine precision for all  ABS(X) .LE. XSMALL.
C            We recommend that  XSMALL < sqrt(eps)/beta, where beta
C            is the floating-point radix (usually 2 or 16).
C
C     Approximate values for some important machines are
C
C                          eps      XMAX     XSMALL      XINF  
C
C  CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
C  CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
C  IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
C  IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
C  IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
C  UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
C  VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38
C
C*******************************************************************
C*******************************************************************
C
C Error Returns
C
C  The program returns the value zero for  X .GT. XMAX, and returns
C    -XINF when BESLY0 is called with a negative or zero argument.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, COS, LOG, SIN, SQRT
C
C
C  Latest modification: June 2, 1989
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division 
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C--------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL
      DOUBLE PRECISION
     1       ARG,AX,CONS,DOWN,EIGHT,FIVE5,FOUR,ONE,ONEOV8,PI2,PJ0,
     2       PJ1,PLG,PROD,PY0,PY1,PY2,P0,P1,P17,QJ0,QJ1,QLG,QY0,QY1,
     3       QY2,Q0,Q1,RESJ,RESULT,R0,R1,SIXTY4,THREE,TWOPI,TWOPI1,
     4       TWOPI2,TWO56,UP,W,WSQ,XDEN,XINF,XMAX,XNUM,XSMALL,XJ0,
     5       XJ1,XJ01,XJ02,XJ11,XJ12,XY,XY0,XY01,XY02,XY1,XY11,XY12,
     6       XY2,XY21,XY22,Z,ZERO,ZSQ
      DIMENSION PJ0(7),PJ1(8),PLG(4),PY0(6),PY1(7),PY2(8),P0(6),P1(6),
     1          QJ0(5),QJ1(7),QLG(4),QY0(5),QY1(6),QY2(7),Q0(5),Q1(5)
C-------------------------------------------------------------------
C  Mathematical constants
C    CONS = ln(.5) + Euler's gamma
C-------------------------------------------------------------------
CS    DATA ZERO,ONE,THREE,FOUR,EIGHT/0.0E0,1.0E0,3.0E0,4.0E0,8.0E0/,
CS   1     FIVE5,SIXTY4,ONEOV8,P17/5.5E0,64.0E0,0.125E0,1.716E-1/,
CS   2     TWO56,CONS/256.0E0,-1.1593151565841244881E-1/,
CS   3     PI2,TWOPI/6.3661977236758134308E-1,6.2831853071795864769E0/,
CS   4     TWOPI1,TWOPI2/6.28125E0,1.9353071795864769253E-3/
      DATA ZERO,ONE,THREE,FOUR,EIGHT/0.0D0,1.0D0,3.0D0,4.0D0,8.0D0/,
     1     FIVE5,SIXTY4,ONEOV8,P17/5.5D0,64.0D0,0.125D0,1.716D-1/,
     2     TWO56,CONS/256.0D0,-1.1593151565841244881D-1/,
     3     PI2,TWOPI/6.3661977236758134308D-1,6.2831853071795864769D0/,
     4     TWOPI1,TWOPI2/6.28125D0,1.9353071795864769253D-3/
C-------------------------------------------------------------------
C  Machine-dependent constants
C-------------------------------------------------------------------
CS    DATA XMAX/8.19E+03/,XSMALL/1.22E-09/,XINF/1.7E+38/
      DATA XMAX/1.07D+09/,XSMALL/9.31D-10/,XINF/1.7D+38/
C-------------------------------------------------------------------
C  Zeroes of Bessel functions
C-------------------------------------------------------------------
CS    DATA XJ0/2.4048255576957727686E+0/,XJ1/5.5200781102863106496E+0/,
CS   1     XY0/8.9357696627916752158E-1/,XY1/3.9576784193148578684E+0/,
CS   2     XY2/7.0860510603017726976E+0/,
CS   3     XJ01/ 616.0E+0/, XJ02/-1.4244423042272313784E-03/,
CS   4     XJ11/1413.0E+0/, XJ12/ 5.4686028631064959660E-04/,
CS   5     XY01/ 228.0E+0/, XY02/ 2.9519662791675215849E-03/,
CS   6     XY11/1013.0E+0/, XY12/ 6.4716931485786837568E-04/,
CS   7     XY21/1814.0E+0/, XY22/ 1.1356030177269762362E-04/
      DATA XJ0/2.4048255576957727686D+0/,XJ1/5.5200781102863106496D+0/,
     1     XY0/8.9357696627916752158D-1/,XY1/3.9576784193148578684D+0/,
     2     XY2/7.0860510603017726976D+0/,
     3     XJ01/ 616.0D+0/, XJ02/-1.4244423042272313784D-03/,
     4     XJ11/1413.0D+0/, XJ12/ 5.4686028631064959660D-04/,
     5     XY01/ 228.0D+0/, XY02/ 2.9519662791675215849D-03/,
     6     XY11/1013.0D+0/, XY12/ 6.4716931485786837568D-04/,
     7     XY21/1814.0D+0/, XY22/ 1.1356030177269762362D-04/
C-------------------------------------------------------------------
C  Coefficients for rational approximation to ln(x/a)
C--------------------------------------------------------------------
CS    DATA PLG/-2.4562334077563243311E+01,2.3642701335621505212E+02,
CS   1         -5.4989956895857911039E+02,3.5687548468071500413E+02/
CS    DATA QLG/-3.5553900764052419184E+01,1.9400230218539473193E+02,
CS   1         -3.3442903192607538956E+02,1.7843774234035750207E+02/
      DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02,
     1         -5.4989956895857911039D+02,3.5687548468071500413D+02/
      DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02,
     1         -3.3442903192607538956D+02,1.7843774234035750207D+02/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C  J0(X) / (X**2 - XJ0**2),  XSMALL  <  |X|  <=  4.0
C--------------------------------------------------------------------
CS    DATA PJ0/6.6302997904833794242E+06,-6.2140700423540120665E+08,
CS   1         2.7282507878605942706E+10,-4.1298668500990866786E+11,
CS   2        -1.2117036164593528341E-01, 1.0344222815443188943E+02,
CS   3        -3.6629814655107086448E+04/
CS    DATA QJ0/4.5612696224219938200E+05, 1.3985097372263433271E+08,
CS   1         2.6328198300859648632E+10, 2.3883787996332290397E+12,
CS   2         9.3614022392337710626E+02/
      DATA PJ0/6.6302997904833794242D+06,-6.2140700423540120665D+08,
     1         2.7282507878605942706D+10,-4.1298668500990866786D+11,
     2        -1.2117036164593528341D-01, 1.0344222815443188943D+02,
     3        -3.6629814655107086448D+04/
      DATA QJ0/4.5612696224219938200D+05, 1.3985097372263433271D+08,
     1         2.6328198300859648632D+10, 2.3883787996332290397D+12,
     2         9.3614022392337710626D+02/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C  J0(X) / (X**2 - XJ1**2),  4.0  <  |X|  <=  8.0
C-------------------------------------------------------------------
CS    DATA PJ1/4.4176707025325087628E+03, 1.1725046279757103576E+04,
CS   1         1.0341910641583726701E+04,-7.2879702464464618998E+03,
CS   2        -1.2254078161378989535E+04,-1.8319397969392084011E+03,
CS   3         4.8591703355916499363E+01, 7.4321196680624245801E+02/
CS    DATA QJ1/3.3307310774649071172E+02,-2.9458766545509337327E+03,
CS   1         1.8680990008359188352E+04,-8.4055062591169562211E+04,
CS   2         2.4599102262586308984E+05,-3.5783478026152301072E+05,
CS   3        -2.5258076240801555057E+01/
      DATA PJ1/4.4176707025325087628D+03, 1.1725046279757103576D+04,
     1         1.0341910641583726701D+04,-7.2879702464464618998D+03,
     2        -1.2254078161378989535D+04,-1.8319397969392084011D+03,
     3         4.8591703355916499363D+01, 7.4321196680624245801D+02/
      DATA QJ1/3.3307310774649071172D+02,-2.9458766545509337327D+03,
     1         1.8680990008359188352D+04,-8.4055062591169562211D+04,
     2         2.4599102262586308984D+05,-3.5783478026152301072D+05,
     3        -2.5258076240801555057D+01/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y0(X) - 2 LN(X/XY0) J0(X)) / (X**2 - XY0**2),
C        XSMALL  <  |X|  <=  3.0
C--------------------------------------------------------------------
CS    DATA PY0/1.0102532948020907590E+04,-2.1287548474401797963E+06,
CS   1         2.0422274357376619816E+08,-8.3716255451260504098E+09,
CS   2         1.0723538782003176831E+11,-1.8402381979244993524E+01/
CS    DATA QY0/6.6475986689240190091E+02, 2.3889393209447253406E+05,
CS   1         5.5662956624278251596E+07, 8.1617187777290363573E+09,
CS   2         5.8873865738997033405E+11/
      DATA PY0/1.0102532948020907590D+04,-2.1287548474401797963D+06,
     1         2.0422274357376619816D+08,-8.3716255451260504098D+09,
     2         1.0723538782003176831D+11,-1.8402381979244993524D+01/
      DATA QY0/6.6475986689240190091D+02, 2.3889393209447253406D+05,
     1         5.5662956624278251596D+07, 8.1617187777290363573D+09,
     2         5.8873865738997033405D+11/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y0(X) - 2 LN(X/XY1) J0(X)) / (X**2 - XY1**2),
C        3.0  <  |X|  <=  5.5
C--------------------------------------------------------------------
CS    DATA PY1/-1.4566865832663635920E+04, 4.6905288611678631510E+06,
CS   1         -6.9590439394619619534E+08, 4.3600098638603061642E+10,
CS   2         -5.5107435206722644429E+11,-2.2213976967566192242E+13,
CS   3          1.7427031242901594547E+01/
CS    DATA QY1/ 8.3030857612070288823E+02, 4.0669982352539552018E+05,
CS   1          1.3960202770986831075E+08, 3.4015103849971240096E+10,
CS   2          5.4266824419412347550E+12, 4.3386146580707264428E+14/
      DATA PY1/-1.4566865832663635920D+04, 4.6905288611678631510D+06,
     1         -6.9590439394619619534D+08, 4.3600098638603061642D+10,
     2         -5.5107435206722644429D+11,-2.2213976967566192242D+13,
     3          1.7427031242901594547D+01/
      DATA QY1/ 8.3030857612070288823D+02, 4.0669982352539552018D+05,
     1          1.3960202770986831075D+08, 3.4015103849971240096D+10,
     2          5.4266824419412347550D+12, 4.3386146580707264428D+14/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y0(X) - 2 LN(X/XY2) J0(X)) / (X**2 - XY2**2),
C        5.5  <  |X|  <=  8.0
C--------------------------------------------------------------------
CS    DATA PY2/ 2.1363534169313901632E+04,-1.0085539923498211426E+07,
CS   1          2.1958827170518100757E+09,-1.9363051266772083678E+11,
CS   2         -1.2829912364088687306E+11, 6.7016641869173237784E+14,
CS   3         -8.0728726905150210443E+15,-1.7439661319197499338E+01/
CS    DATA QY2/ 8.7903362168128450017E+02, 5.3924739209768057030E+05,
CS   1          2.4727219475672302327E+08, 8.6926121104209825246E+10,
CS   2          2.2598377924042897629E+13, 3.9272425569640309819E+15,
CS   3          3.4563724628846457519E+17/
      DATA PY2/ 2.1363534169313901632D+04,-1.0085539923498211426D+07,
     1          2.1958827170518100757D+09,-1.9363051266772083678D+11,
     2         -1.2829912364088687306D+11, 6.7016641869173237784D+14,
     3         -8.0728726905150210443D+15,-1.7439661319197499338D+01/
      DATA QY2/ 8.7903362168128450017D+02, 5.3924739209768057030D+05,
     1          2.4727219475672302327D+08, 8.6926121104209825246D+10,
     2          2.2598377924042897629D+13, 3.9272425569640309819D+15,
     3          3.4563724628846457519D+17/
C-------------------------------------------------------------------
C  Coefficients for Hart,s approximation,  |X| > 8.0
C-------------------------------------------------------------------
CS    DATA P0/3.4806486443249270347E+03, 2.1170523380864944322E+04,
CS   1        4.1345386639580765797E+04, 2.2779090197304684302E+04,
CS   2        8.8961548424210455236E-01, 1.5376201909008354296E+02/
CS    DATA Q0/3.5028735138235608207E+03, 2.1215350561880115730E+04,
CS   1        4.1370412495510416640E+04, 2.2779090197304684318E+04,
CS   2        1.5711159858080893649E+02/
CS    DATA P1/-2.2300261666214198472E+01,-1.1183429920482737611E+02,
CS   1        -1.8591953644342993800E+02,-8.9226600200800094098E+01,
CS   2        -8.8033303048680751817E-03,-1.2441026745835638459E+00/
CS    DATA Q1/1.4887231232283756582E+03, 7.2642780169211018836E+03,
CS   1        1.1951131543434613647E+04, 5.7105024128512061905E+03,
CS   2        9.0593769594993125859E+01/
      DATA P0/3.4806486443249270347D+03, 2.1170523380864944322D+04,
     1        4.1345386639580765797D+04, 2.2779090197304684302D+04,
     2        8.8961548424210455236D-01, 1.5376201909008354296D+02/
      DATA Q0/3.5028735138235608207D+03, 2.1215350561880115730D+04,
     1        4.1370412495510416640D+04, 2.2779090197304684318D+04,
     2        1.5711159858080893649D+02/
      DATA P1/-2.2300261666214198472D+01,-1.1183429920482737611D+02,
     1        -1.8591953644342993800D+02,-8.9226600200800094098D+01,
     2        -8.8033303048680751817D-03,-1.2441026745835638459D+00/
      DATA Q1/1.4887231232283756582D+03, 7.2642780169211018836D+03,
     1        1.1951131543434613647D+04, 5.7105024128512061905D+03,
     2        9.0593769594993125859D+01/
C-------------------------------------------------------------------
C  Check for error conditions
C-------------------------------------------------------------------
      AX = ABS(ARG)
      IF ((JINT .EQ. 1) .AND. (ARG .LE. ZERO)) THEN
            RESULT = -XINF
            GO TO 2000
         ELSE IF (AX .GT. XMAX) THEN
            RESULT = ZERO
            GO TO 2000
      END IF
      IF (AX .GT. EIGHT) GO TO 800
      IF (AX .LE. XSMALL) THEN
         IF (JINT .EQ. 0) THEN
               RESULT = ONE
            ELSE
               RESULT = PI2 * (LOG(AX) + CONS)
         END IF
         GO TO 2000
      END IF
C-------------------------------------------------------------------
C  Calculate J0 for appropriate interval, preserving
C     accuracy near the zero of J0
C-------------------------------------------------------------------
      ZSQ = AX * AX
      IF (AX .LE. FOUR) THEN
            XNUM = (PJ0(5) * ZSQ + PJ0(6)) * ZSQ + PJ0(7)
            XDEN = ZSQ + QJ0(5)
            DO 50 I = 1, 4
               XNUM = XNUM * ZSQ + PJ0(I)
               XDEN = XDEN * ZSQ + QJ0(I)
   50       CONTINUE
            PROD = ((AX - XJ01/TWO56) - XJ02) * (AX + XJ0)
         ELSE
            WSQ = ONE - ZSQ / SIXTY4
            XNUM = PJ1(7) * WSQ + PJ1(8)
            XDEN = WSQ + QJ1(7)
            DO 220 I = 1, 6
               XNUM = XNUM * WSQ + PJ1(I)
               XDEN = XDEN * WSQ + QJ1(I)
  220       CONTINUE
            PROD = (AX + XJ1) * ((AX - XJ11/TWO56) - XJ12)
      END IF
      RESULT = PROD * XNUM / XDEN
      IF (JINT .EQ. 0) GO TO 2000
C-------------------------------------------------------------------
C  Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
C    where xn is a zero of Y0
C-------------------------------------------------------------------
      IF (AX .LE. THREE) THEN
            UP = (AX-XY01/TWO56)-XY02
            XY = XY0
         ELSE IF (AX .LE. FIVE5) THEN
            UP = (AX-XY11/TWO56)-XY12
            XY = XY1
         ELSE
            UP = (AX-XY21/TWO56)-XY22
            XY = XY2
      END IF
      DOWN = AX + XY
      IF (ABS(UP) .LT. P17*DOWN) THEN
            W = UP/DOWN
            WSQ = W*W
            XNUM = PLG(1)
            XDEN = WSQ + QLG(1)
            DO 320 I = 2, 4
               XNUM = XNUM*WSQ + PLG(I)
               XDEN = XDEN*WSQ + QLG(I)
  320       CONTINUE
            RESJ = PI2 * RESULT * W * XNUM/XDEN
         ELSE
            RESJ = PI2 * RESULT * LOG(AX/XY)
      END IF
C-------------------------------------------------------------------
C  Now calculate Y0 for appropriate interval, preserving
C     accuracy near the zero of Y0
C-------------------------------------------------------------------
      IF (AX .LE. THREE) THEN
            XNUM = PY0(6) * ZSQ + PY0(1)
            XDEN = ZSQ + QY0(1)
            DO 340 I = 2, 5
               XNUM = XNUM * ZSQ + PY0(I)
               XDEN = XDEN * ZSQ + QY0(I)
  340       CONTINUE
         ELSE IF (AX .LE. FIVE5) THEN
            XNUM = PY1(7) * ZSQ + PY1(1)
            XDEN = ZSQ + QY1(1)
            DO 360 I = 2, 6
               XNUM = XNUM * ZSQ + PY1(I)
               XDEN = XDEN * ZSQ + QY1(I)
  360       CONTINUE
         ELSE
            XNUM = PY2(8) * ZSQ + PY2(1)
            XDEN = ZSQ + QY2(1)
            DO 380 I = 2, 7
               XNUM = XNUM * ZSQ + PY2(I)
               XDEN = XDEN * ZSQ + QY2(I)
  380       CONTINUE
      END IF
      RESULT = RESJ + UP * DOWN * XNUM / XDEN
      GO TO 2000
C-------------------------------------------------------------------
C  Calculate J0 or Y0 for |ARG|  >  8.0
C-------------------------------------------------------------------
  800 Z = EIGHT / AX
      W = AX / TWOPI
      W = AINT(W) + ONEOV8
      W = (AX - W * TWOPI1) - W * TWOPI2
      ZSQ = Z * Z
      XNUM = P0(5) * ZSQ + P0(6)
      XDEN = ZSQ + Q0(5)
      UP = P1(5) * ZSQ + P1(6)
      DOWN = ZSQ + Q1(5)
      DO 850 I = 1, 4
         XNUM = XNUM * ZSQ + P0(I)
         XDEN = XDEN * ZSQ + Q0(I)
         UP = UP * ZSQ + P1(I)
         DOWN = DOWN * ZSQ + Q1(I)
  850 CONTINUE
      R0 = XNUM / XDEN
      R1 = UP / DOWN
      IF (JINT .EQ. 0) THEN
            RESULT = SQRT(PI2/AX) * (R0*COS(W) - Z*R1*SIN(W))
         ELSE
            RESULT = SQRT(PI2/AX) * (R0*SIN(W) + Z*R1*COS(W))
      END IF
 2000 RETURN
C---------- Last line of CALJY0 ----------
      END
C
      DOUBLE PRECISION FUNCTION BESJ0(X)
CS    REAL FUNCTION BESJ0(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for Bessel functions
C   of the first kind of order zero for arguments  |X| <= XMAX
C   (see comments heading CALJY0).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL  X, RESULT
      DOUBLE PRECISION  X, RESULT
C--------------------------------------------------------------------
      JINT=0
      CALL CALJY0(X,RESULT,JINT)
      BESJ0 = RESULT
      RETURN
C---------- Last line of BESJ0 ----------
      END
C
      DOUBLE PRECISION FUNCTION BESY0(X)
CS    REAL FUNCTION BESY0(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for Bessel functions
C   of the second kind of order zero for arguments 0 < X <= XMAX
C   (see comments heading CALJY0).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL  X, RESULT
      DOUBLE PRECISION  X, RESULT
C--------------------------------------------------------------------
      JINT=1
      CALL CALJY0(X,RESULT,JINT)
      BESY0 = RESULT
      RETURN
C---------- Last line of BESY0 ----------
      END
C
      SUBROUTINE CALJY1(ARG,RESULT,JINT)
C---------------------------------------------------------------------
C
C This packet computes first-order Bessel functions of the first and
C   second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
C   for Y1, and |X| <= XMAX for J1.  It contains two function-type
C   subprograms,  BESJ1  and  BESY1,  and one subroutine-type
C   subprogram,  CALJY1.  The calling statements for the primary
C   entries are:
C
C           Y = BESJ1(X)
C   and
C           Y = BESY1(X),
C
C   where the entry points correspond to the functions J1(X) and Y1(X),
C   respectively.  The routine  CALJY1  is intended for internal packet
C   use only, all computations within the packet being concentrated in
C   this one routine.  The function subprograms invoke  CALJY1  with
C   the statement
C           CALL CALJY1(ARG,RESULT,JINT),
C   where the parameter usage is as follows:
C
C      Function                  Parameters for CALJY1
C       call              ARG             RESULT          JINT
C
C     BESJ1(ARG)     |ARG| .LE. XMAX       J1(ARG)          0
C     BESY1(ARG)   0 .LT. ARG .LE. XMAX    Y1(ARG)          1
C
C   The main computation uses unpublished minimax rational
C   approximations for X .LE. 8.0, and an approximation from the 
C   book  Computer Approximations  by Hart, et. al., Wiley and Sons, 
C   New York, 1968, for arguments larger than 8.0   Part of this
C   transportable packet is patterned after the machine-dependent
C   FUNPACK program BESJ1(X), but cannot match that version for
C   efficiency or accuracy.  This version uses rational functions
C   that are theoretically accurate to at least 18 significant decimal
C   digits for X <= 8, and at least 18 decimal places for X > 8.  The
C   accuracy achieved depends on the arithmetic system, the compiler,
C   the intrinsic functions, and proper selection of the machine-
C   dependent constants.
C
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   XINF   = largest positive machine number
C   XMAX   = largest acceptable argument.  The functions AINT, SIN
C            and COS must perform properly for  ABS(X) .LE. XMAX.
C            We recommend that XMAX be a small integer multiple of
C            sqrt(1/eps), where eps is the smallest positive number
C            such that  1+eps > 1. 
C   XSMALL = positive argument such that  1.0-(1/2)(X/2)**2 = 1.0
C            to machine precision for all  ABS(X) .LE. XSMALL.
C            We recommend that  XSMALL < sqrt(eps)/beta, where beta
C            is the floating-point radix (usually 2 or 16).
C
C     Approximate values for some important machines are
C
C                          eps      XMAX     XSMALL      XINF  
C
C  CDC 7600      (S.P.)  7.11E-15  1.34E+08  2.98E-08  1.26E+322
C  CRAY-1        (S.P.)  7.11E-15  1.34E+08  2.98E-08  5.45E+2465
C  IBM PC (8087) (S.P.)  5.96E-08  8.19E+03  1.22E-04  3.40E+38
C  IBM PC (8087) (D.P.)  1.11D-16  2.68D+08  3.72D-09  1.79D+308
C  IBM 195       (D.P.)  2.22D-16  6.87D+09  9.09D-13  7.23D+75
C  UNIVAC 1108   (D.P.)  1.73D-18  4.30D+09  2.33D-10  8.98D+307
C  VAX 11/780    (D.P.)  1.39D-17  1.07D+09  9.31D-10  1.70D+38
C
C*******************************************************************
C*******************************************************************
C
C Error Returns
C
C  The program returns the value zero for  X .GT. XMAX, and returns
C    -XINF when BESLY1 is called with a negative or zero argument.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, COS, LOG, SIN, SQRT
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division 
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: November 10, 1987
C
C--------------------------------------------------------------------
      INTEGER I,JINT
      DIMENSION PJ0(7),PJ1(8),PLG(4),PY0(7),PY1(9),P0(6),P1(6),
     1          QJ0(5),QJ1(7),QLG(4),QY0(6),QY1(8),Q0(6),Q1(6)
CS    REAL
      DOUBLE PRECISION
     1   ARG,AX,DOWN,EIGHT,FOUR,HALF,PI2,PJ0,PJ1,PLG,PROD,PY0,
     2   PY1,P0,P1,P17,QJ0,QJ1,QLG,QY0,QY1,Q0,Q1,RESJ,RESULT,
     3   RTPI2,R0,R1,THROV8,TWOPI,TWOPI1,TWOPI2,TWO56,UP,W,WSQ,
     4   XDEN,XINF,XMAX,XNUM,XSMALL,XJ0,XJ1,XJ01,XJ02,XJ11,XJ12,
     5   XY,XY0,XY01,XY02,XY1,XY11,XY12,Z,ZERO,ZSQ
C-------------------------------------------------------------------
C  Mathematical constants
C-------------------------------------------------------------------
CS    DATA EIGHT/8.0E0/,
CS   1     FOUR/4.0E0/,HALF/0.5E0/,THROV8/0.375E0/,
CS   2     PI2/6.3661977236758134308E-1/,P17/1.716E-1/
CS   3     TWOPI/6.2831853071795864769E+0/,ZERO/0.0E0/,
CS   4     TWOPI1/6.28125E0/,TWOPI2/1.9353071795864769253E-03/
CS   5     TWO56/256.0E+0/,RTPI2/7.9788456080286535588E-1/
      DATA EIGHT/8.0D0/,
     1     FOUR/4.0D0/,HALF/0.5D0/,THROV8/0.375D0/,
     2     PI2/6.3661977236758134308D-1/,P17/1.716D-1/
     3     TWOPI/6.2831853071795864769D+0/,ZERO/0.0D0/,
     4     TWOPI1/6.28125D0/,TWOPI2/1.9353071795864769253D-03/
     5     TWO56/256.0D+0/,RTPI2/7.9788456080286535588D-1/
C-------------------------------------------------------------------
C  Machine-dependent constants
C-------------------------------------------------------------------
CS    DATA XMAX/8.19E+03/,XSMALL/1.22E-09/,XINF/1.7E+38/
      DATA XMAX/1.07D+09/,XSMALL/9.31D-10/,XINF/1.7D+38/
C-------------------------------------------------------------------
C  Zeroes of Bessel functions
C-------------------------------------------------------------------
CS    DATA XJ0/3.8317059702075123156E+0/,XJ1/7.0155866698156187535E+0/,
CS   1     XY0/2.1971413260310170351E+0/,XY1/5.4296810407941351328E+0/,
CS   2     XJ01/ 981.0E+0/, XJ02/-3.2527979248768438556E-04/,
CS   3     XJ11/1796.0E+0/, XJ12/-3.8330184381246462950E-05/,
CS   4     XY01/ 562.0E+0/, XY02/ 1.8288260310170351490E-03/,
CS   5     XY11/1390.0E+0/, XY12/-6.4592058648672279948E-06/
      DATA XJ0/3.8317059702075123156D+0/,XJ1/7.0155866698156187535D+0/,
     1     XY0/2.1971413260310170351D+0/,XY1/5.4296810407941351328D+0/,
     2     XJ01/ 981.0D+0/, XJ02/-3.2527979248768438556D-04/,
     3     XJ11/1796.0D+0/, XJ12/-3.8330184381246462950D-05/,
     4     XY01/ 562.0D+0/, XY02/ 1.8288260310170351490D-03/,
     5     XY11/1390.0D+0/, XY12/-6.4592058648672279948D-06/
C-------------------------------------------------------------------
C  Coefficients for rational approximation to ln(x/a)
C--------------------------------------------------------------------
CS    DATA PLG/-2.4562334077563243311E+01,2.3642701335621505212E+02,
CS   1         -5.4989956895857911039E+02,3.5687548468071500413E+02/
CS    DATA QLG/-3.5553900764052419184E+01,1.9400230218539473193E+02,
CS   1         -3.3442903192607538956E+02,1.7843774234035750207E+02/
      DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02,
     1         -5.4989956895857911039D+02,3.5687548468071500413D+02/
      DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02,
     1         -3.3442903192607538956D+02,1.7843774234035750207D+02/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C  J1(X) / (X * (X**2 - XJ0**2)),  XSMALL  <  |X|  <=  4.0
C--------------------------------------------------------------------
CS    DATA PJ0/9.8062904098958257677E+05,-1.1548696764841276794E+08,
CS   1       6.6781041261492395835E+09,-1.4258509801366645672E+11,
CS   2      -4.4615792982775076130E+03, 1.0650724020080236441E+01,
CS   3      -1.0767857011487300348E-02/
CS    DATA QJ0/5.9117614494174794095E+05, 2.0228375140097033958E+08,
CS   1       4.2091902282580133541E+10, 4.1868604460820175290E+12,
CS   2       1.0742272239517380498E+03/
      DATA PJ0/9.8062904098958257677D+05,-1.1548696764841276794D+08,
     1       6.6781041261492395835D+09,-1.4258509801366645672D+11,
     2      -4.4615792982775076130D+03, 1.0650724020080236441D+01,
     3      -1.0767857011487300348D-02/
      DATA QJ0/5.9117614494174794095D+05, 2.0228375140097033958D+08,
     1       4.2091902282580133541D+10, 4.1868604460820175290D+12,
     2       1.0742272239517380498D+03/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C  J1(X) / (X * (X**2 - XJ1**2)),  4.0  <  |X|  <=  8.0
C-------------------------------------------------------------------
CS    DATA PJ1/4.6179191852758252280E+00,-7.1329006872560947377E+03,
CS   1       4.5039658105749078904E+06,-1.4437717718363239107E+09,
CS   2       2.3569285397217157313E+11,-1.6324168293282543629E+13,
CS   3       1.1357022719979468624E+14, 1.0051899717115285432E+15/
CS    DATA QJ1/1.1267125065029138050E+06, 6.4872502899596389593E+08,
CS   1       2.7622777286244082666E+11, 8.4899346165481429307E+13,
CS   2       1.7128800897135812012E+16, 1.7253905888447681194E+18,
CS   3       1.3886978985861357615E+03/
      DATA PJ1/4.6179191852758252280D+00,-7.1329006872560947377D+03,
     1       4.5039658105749078904D+06,-1.4437717718363239107D+09,
     2       2.3569285397217157313D+11,-1.6324168293282543629D+13,
     3       1.1357022719979468624D+14, 1.0051899717115285432D+15/
      DATA QJ1/1.1267125065029138050D+06, 6.4872502899596389593D+08,
     1       2.7622777286244082666D+11, 8.4899346165481429307D+13,
     2       1.7128800897135812012D+16, 1.7253905888447681194D+18,
     3       1.3886978985861357615D+03/
C-------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y1(X) - 2 LN(X/XY0) J1(X)) / (X**2 - XY0**2),
C        XSMALL  <  |X|  <=  4.0
C--------------------------------------------------------------------
CS    DATA PY0/2.2157953222280260820E+05,-5.9157479997408395984E+07,
CS   1         7.2144548214502560419E+09,-3.7595974497819597599E+11,
CS   2         5.4708611716525426053E+12, 4.0535726612579544093E+13,
CS   3        -3.1714424660046133456E+02/
CS    DATA QY0/8.2079908168393867438E+02, 3.8136470753052572164E+05,
CS   1         1.2250435122182963220E+08, 2.7800352738690585613E+10,
CS   2         4.1272286200406461981E+12, 3.0737873921079286084E+14/
      DATA PY0/2.2157953222280260820D+05,-5.9157479997408395984D+07,
     1         7.2144548214502560419D+09,-3.7595974497819597599D+11,
     2         5.4708611716525426053D+12, 4.0535726612579544093D+13,
     3        -3.1714424660046133456D+02/
      DATA QY0/8.2079908168393867438D+02, 3.8136470753052572164D+05,
     1         1.2250435122182963220D+08, 2.7800352738690585613D+10,
     2         4.1272286200406461981D+12, 3.0737873921079286084D+14/
C--------------------------------------------------------------------
C  Coefficients for rational approximation of
C    (Y1(X) - 2 LN(X/XY1) J1(X)) / (X**2 - XY1**2),
C        4.0  <  |X|  <=  8.0
C--------------------------------------------------------------------
CS    DATA PY1/ 1.9153806858264202986E+06,-1.1957961912070617006E+09,
CS   1          3.7453673962438488783E+11,-5.9530713129741981618E+13,
CS   2          4.0686275289804744814E+15,-2.3638408497043134724E+16,
CS   3         -5.6808094574724204577E+18, 1.1514276357909013326E+19,
CS   4         -1.2337180442012953128E+03/
CS    DATA QY1/ 1.2855164849321609336E+03, 1.0453748201934079734E+06,
CS   1          6.3550318087088919566E+08, 3.0221766852960403645E+11,
CS   2          1.1187010065856971027E+14, 3.0837179548112881950E+16,
CS   3          5.6968198822857178911E+18, 5.3321844313316185697E+20/
      DATA PY1/ 1.9153806858264202986D+06,-1.1957961912070617006D+09,
     1          3.7453673962438488783D+11,-5.9530713129741981618D+13,
     2          4.0686275289804744814D+15,-2.3638408497043134724D+16,
     3         -5.6808094574724204577D+18, 1.1514276357909013326D+19,
     4         -1.2337180442012953128D+03/
      DATA QY1/ 1.2855164849321609336D+03, 1.0453748201934079734D+06,
     1          6.3550318087088919566D+08, 3.0221766852960403645D+11,
     2          1.1187010065856971027D+14, 3.0837179548112881950D+16,
     3          5.6968198822857178911D+18, 5.3321844313316185697D+20/
C-------------------------------------------------------------------
C  Coefficients for Hart,s approximation,  |X| > 8.0
C-------------------------------------------------------------------
CS    DATA P0/-1.0982405543459346727E+05,-1.5235293511811373833E+06,
CS   1         -6.6033732483649391093E+06,-9.9422465050776411957E+06,
CS   2         -4.4357578167941278571E+06,-1.6116166443246101165E+03/
CS    DATA Q0/-1.0726385991103820119E+05,-1.5118095066341608816E+06,
CS   1         -6.5853394797230870728E+06,-9.9341243899345856590E+06,
CS   2         -4.4357578167941278568E+06,-1.4550094401904961825E+03/
CS    DATA P1/ 1.7063754290207680021E+03, 1.8494262873223866797E+04,
CS   1          6.6178836581270835179E+04, 8.5145160675335701966E+04,
CS   2          3.3220913409857223519E+04, 3.5265133846636032186E+01/
CS    DATA Q1/ 3.7890229745772202641E+04, 4.0029443582266975117E+05,
CS   1          1.4194606696037208929E+06, 1.8194580422439972989E+06,
CS   2          7.0871281941028743574E+05, 8.6383677696049909675E+02/
      DATA P0/-1.0982405543459346727D+05,-1.5235293511811373833D+06,
     1         -6.6033732483649391093D+06,-9.9422465050776411957D+06,
     2         -4.4357578167941278571D+06,-1.6116166443246101165D+03/
      DATA Q0/-1.0726385991103820119D+05,-1.5118095066341608816D+06,
     1         -6.5853394797230870728D+06,-9.9341243899345856590D+06,
     2         -4.4357578167941278568D+06,-1.4550094401904961825D+03/
      DATA P1/ 1.7063754290207680021D+03, 1.8494262873223866797D+04,
     1          6.6178836581270835179D+04, 8.5145160675335701966D+04,
     2          3.3220913409857223519D+04, 3.5265133846636032186D+01/
      DATA Q1/ 3.7890229745772202641D+04, 4.0029443582266975117D+05,
     1          1.4194606696037208929D+06, 1.8194580422439972989D+06,
     2          7.0871281941028743574D+05, 8.6383677696049909675D+02/
C-------------------------------------------------------------------
C  Check for error conditions
C-------------------------------------------------------------------
      AX = ABS(ARG)
      IF ((JINT .EQ. 1) .AND. ((ARG .LE. ZERO) .OR.
     1   ((ARG .LT. HALF) .AND. (AX*XINF .LT. PI2)))) THEN
            RESULT = -XINF
            GO TO 2000
         ELSE IF (AX .GT. XMAX) THEN
            RESULT = ZERO
            GO TO 2000
      END IF
      IF (AX .GT. EIGHT) THEN
            GO TO 800
         ELSE IF (AX .LE. XSMALL) THEN
            IF (JINT .EQ. 0) THEN
                  RESULT = ARG * HALF
               ELSE
                  RESULT = -PI2 / AX
            END IF
            GO TO 2000
      END IF
C-------------------------------------------------------------------
C  Calculate J1 for appropriate interval, preserving
C     accuracy near the zero of J1
C-------------------------------------------------------------------
      ZSQ = AX * AX
      IF (AX .LE. FOUR) THEN
            XNUM = (PJ0(7) * ZSQ + PJ0(6)) * ZSQ + PJ0(5)
            XDEN = ZSQ + QJ0(5)
            DO 50 I = 1, 4
               XNUM = XNUM * ZSQ + PJ0(I)
               XDEN = XDEN * ZSQ + QJ0(I)
   50       CONTINUE
            PROD = ARG * ((AX - XJ01/TWO56) - XJ02) * (AX + XJ0)
         ELSE
            XNUM = PJ1(1)
            XDEN = (ZSQ + QJ1(7)) * ZSQ + QJ1(1)
            DO 220 I = 2, 6
               XNUM = XNUM * ZSQ + PJ1(I)
               XDEN = XDEN * ZSQ + QJ1(I)
  220       CONTINUE
            XNUM = XNUM * (AX - EIGHT) * (AX + EIGHT) + PJ1(7)
            XNUM = XNUM * (AX - FOUR) * (AX + FOUR) + PJ1(8)
            PROD = ARG * ((AX - XJ11/TWO56) - XJ12) * (AX + XJ1)
      END IF
      RESULT = PROD * (XNUM / XDEN)
      IF (JINT .EQ. 0) GO TO 2000
C-------------------------------------------------------------------
C  Calculate Y1.  First find  RESJ = pi/2 ln(x/xn) J1(x),
C    where xn is a zero of Y1
C-------------------------------------------------------------------
      IF (AX .LE. FOUR) THEN
            UP = (AX-XY01/TWO56)-XY02
            XY = XY0
         ELSE
            UP = (AX-XY11/TWO56)-XY12
            XY = XY1
      END IF
      DOWN = AX + XY
      IF (ABS(UP) .LT. P17*DOWN) THEN
            W = UP/DOWN
            WSQ = W*W
            XNUM = PLG(1)
            XDEN = WSQ + QLG(1)
            DO 320 I = 2, 4
               XNUM = XNUM*WSQ + PLG(I)
               XDEN = XDEN*WSQ + QLG(I)
  320       CONTINUE
            RESJ = PI2 * RESULT * W * XNUM/XDEN
         ELSE
            RESJ = PI2 * RESULT * LOG(AX/XY)
      END IF
C-------------------------------------------------------------------
C  Now calculate Y1 for appropriate interval, preserving
C     accuracy near the zero of Y1
C-------------------------------------------------------------------
      IF (AX .LE. FOUR) THEN
            XNUM = PY0(7) * ZSQ + PY0(1)
            XDEN = ZSQ + QY0(1)
            DO 340 I = 2, 6
               XNUM = XNUM * ZSQ + PY0(I)
               XDEN = XDEN * ZSQ + QY0(I)
  340       CONTINUE
         ELSE
            XNUM = PY1(9) * ZSQ + PY1(1)
            XDEN = ZSQ + QY1(1)
            DO 360 I = 2, 8
               XNUM = XNUM * ZSQ + PY1(I)
               XDEN = XDEN * ZSQ + QY1(I)
  360       CONTINUE
      END IF
      RESULT = RESJ + (UP*DOWN/AX) * XNUM / XDEN
      GO TO 2000
C-------------------------------------------------------------------
C  Calculate J1 or Y1 for |ARG|  >  8.0
C-------------------------------------------------------------------
  800 Z = EIGHT / AX
      W = AINT(AX/TWOPI) + THROV8
      W = (AX - W * TWOPI1) - W * TWOPI2
      ZSQ = Z * Z
      XNUM = P0(6)
      XDEN = ZSQ + Q0(6)
      UP = P1(6)
      DOWN = ZSQ + Q1(6)
      DO 850 I = 1, 5
         XNUM = XNUM * ZSQ + P0(I)
         XDEN = XDEN * ZSQ + Q0(I)
         UP = UP * ZSQ + P1(I)
         DOWN = DOWN * ZSQ + Q1(I)
  850 CONTINUE
      R0 = XNUM / XDEN
      R1 = UP / DOWN
      IF (JINT .EQ. 0) THEN
            RESULT = (RTPI2/SQRT(AX)) * (R0*COS(W) - Z*R1*SIN(W))
         ELSE
            RESULT = (RTPI2/SQRT(AX)) * (R0*SIN(W) + Z*R1*COS(W))
      END IF
      IF ((JINT .EQ. 0) .AND. (ARG .LT. ZERO)) RESULT = -RESULT
 2000 RETURN
C---------- Last card of CALJY1 ----------
      END
C
      FUNCTION BESJ1(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for Bessel functions
C   of the first kind of order zero for arguments  |X| <= XMAX
C   (see comments heading CALJY1).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL
      DOUBLE PRECISION
     1   BESJ1,RESULT,X
C--------------------------------------------------------------------
      JINT=0
      CALL CALJY1(X,RESULT,JINT)
      BESJ1 = RESULT
      RETURN
C---------- Last card of BESJ1 ----------
      END
C
      FUNCTION BESY1(X)
C--------------------------------------------------------------------
C
C This subprogram computes approximate values for Bessel functions
C   of the second kind of order zero for arguments 0 < X <= XMAX
C   (see comments heading CALJY1).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL
      DOUBLE PRECISION
     1   BESY1,RESULT,X
C--------------------------------------------------------------------
      JINT=1
      CALL CALJY1(X,RESULT,JINT)
      BESY1 = RESULT
      RETURN
C---------- Last card of BESY1 ----------
      END
C
      SUBROUTINE CALCI0(ARG,RESULT,JINT)
C--------------------------------------------------------------------
C
C This packet computes modified Bessel functions of the first kind
C   and order zero, I0(X) and EXP(-ABS(X))*I0(X), for real
C   arguments X.  It contains two function type subprograms, BESI0
C   and BESEI0, and one subroutine type subprogram, CALCI0.
C   The calling statements for the primary entries are 
C
C                   Y=BESI0(X)
C   and
C                   Y=BESEI0(X)
C
C   where the entry points correspond to the functions I0(X) and
C   EXP(-ABS(X))*I0(X), respectively.  The routine CALCI0 is
C   intended for internal packet use only, all computations within
C   the packet being concentrated in this routine.  The function
C   subprograms invoke CALCI0 with the statement
C          CALL CALCI0(ARG,RESULT,JINT)
C   where the parameter usage is as follows
C
C      Function                     Parameters for CALCI0
C       Call              ARG                  RESULT          JINT
C
C     BESI0(ARG)    ABS(ARG) .LE. XMAX        I0(ARG)           1
C     BESEI0(ARG)    any real ARG        EXP(-ABS(ARG))*I0(ARG) 2
C
C   The main computation evaluates slightly modified forms of
C   minimax approximations generated by Blair and Edwards, Chalk
C   River (Atomic Energy of Canada Limited) Report AECL-4928,
C   October, 1974.  This transportable program is patterned after 
C   the machine-dependent FUNPACK packet NATSI0, but cannot match
C   that version for efficiency or accuracy.  This version uses
C   rational functions that theoretically approximate I-SUB-0(X)
C   to at least 18 significant decimal digits.  The accuracy
C   achieved depends on the arithmetic system, the compiler, the
C   intrinsic functions, and proper selection of the machine-
C   dependent constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta   = Radix for the floating-point system
C   maxexp = Smallest power of beta that overflows
C   XSMALL = Positive argument such that 1.0 - X = 1.0 to
C            machine precision for all ABS(X) .LE. XSMALL.
C   XINF =   Largest positive machine number; approximately
C            beta**maxexp
C   XMAX =   Largest argument acceptable to BESI0;  Solution to
C            equation:  
C               W(X) * (1+1/(8*X)+9/(128*X**2) = beta**maxexp
C            where  W(X) = EXP(X)/SQRT(2*PI*X)
C
C
C     Approximate values for some important machines are:
C
C                          beta       maxexp       XSMALL
C
C CRAY-1        (S.P.)       2         8191       3.55E-15
C Cyber 180/855
C   under NOS   (S.P.)       2         1070       3.55E-15
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)       2          128       2.98E-8
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)       2         1024       5.55D-17
C IBM 3033      (D.P.)      16           63       6.95D-18
C VAX           (S.P.)       2          127       2.98E-8
C VAX D-Format  (D.P.)       2          127       6.95D-18
C VAX G-Format  (D.P.)       2         1023       5.55D-17
C
C
C                               XINF          XMAX
C
C CRAY-1        (S.P.)       5.45E+2465     5682.810
C Cyber 180/855
C   under NOS   (S.P.)       1.26E+322       745.893
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)       3.40E+38         91.900
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)       1.79D+308       713.986
C IBM 3033      (D.P.)       7.23D+75        178.182
C VAX           (S.P.)       1.70D+38         91.203
C VAX D-Format  (D.P.)       1.70D+38         91.203
C VAX G-Format  (D.P.)       8.98D+307       713.293
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns XINF for BESI0 for ABS(ARG) .GT. XMAX.
C
C
C  Intrinsic functions required are:
C
C     ABS, SQRT, EXP
C
C
C  Authors: W. J. Cody and L. Stoltz
C           Mathematics and Computer Science Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C  Latest modification: June 7, 1988
C
C--------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL
      DOUBLE PRECISION
     1       A,ARG,B,EXP40,FORTY,ONE,ONE5,P,PP,Q,QQ,RESULT,
     2       REC15,SUMP,SUMQ,TWO25,X,XINF,XMAX,XSMALL,XX
      DIMENSION P(15),PP(8),Q(5),QQ(7)
C--------------------------------------------------------------------
C  Mathematical constants
C--------------------------------------------------------------------
CS    DATA ONE/1.0E0/,ONE5/15.0E0/,EXP40/2.353852668370199854E17/,
CS   1     FORTY/40.0E0/,REC15/6.6666666666666666666E-2/,
CS   2     TWO25/225.0E0/
      DATA ONE/1.0D0/,ONE5/15.0D0/,EXP40/2.353852668370199854D17/,
     1     FORTY/40.0D0/,REC15/6.6666666666666666666D-2/,
     2     TWO25/225.0D0/
C--------------------------------------------------------------------
C  Machine-dependent constants
C--------------------------------------------------------------------
CS    DATA XSMALL/2.98E-8/,XINF/3.40E38/,XMAX/91.9E0/
      DATA XSMALL/5.55D-17/,XINF/1.79D308/,XMAX/713.986D0/
C--------------------------------------------------------------------
C  Coefficients for XSMALL .LE. ABS(ARG) .LT. 15.0
C--------------------------------------------------------------------
CS    DATA  P/-5.2487866627945699800E-18,-1.5982226675653184646E-14,
CS   1        -2.6843448573468483278E-11,-3.0517226450451067446E-08,
CS   2        -2.5172644670688975051E-05,-1.5453977791786851041E-02,
CS   3        -7.0935347449210549190E+00,-2.4125195876041896775E+03,
CS   4        -5.9545626019847898221E+05,-1.0313066708737980747E+08,
CS   5        -1.1912746104985237192E+10,-8.4925101247114157499E+11,
CS   6        -3.2940087627407749166E+13,-5.5050369673018427753E+14,
CS   7        -2.2335582639474375249E+15/
CS    DATA  Q/-3.7277560179962773046E+03, 6.5158506418655165707E+06,
CS   1        -6.5626560740833869295E+09, 3.7604188704092954661E+12,
CS   2        -9.7087946179594019126E+14/
      DATA  P/-5.2487866627945699800D-18,-1.5982226675653184646D-14,
     1        -2.6843448573468483278D-11,-3.0517226450451067446D-08,
     2        -2.5172644670688975051D-05,-1.5453977791786851041D-02,
     3        -7.0935347449210549190D+00,-2.4125195876041896775D+03,
     4        -5.9545626019847898221D+05,-1.0313066708737980747D+08,
     5        -1.1912746104985237192D+10,-8.4925101247114157499D+11,
     6        -3.2940087627407749166D+13,-5.5050369673018427753D+14,
     7        -2.2335582639474375249D+15/
      DATA  Q/-3.7277560179962773046D+03, 6.5158506418655165707D+06,
     1        -6.5626560740833869295D+09, 3.7604188704092954661D+12,
     2        -9.7087946179594019126D+14/
C--------------------------------------------------------------------
C  Coefficients for 15.0 .LE. ABS(ARG)
C--------------------------------------------------------------------
CS    DATA PP/-3.9843750000000000000E-01, 2.9205384596336793945E+00,
CS   1        -2.4708469169133954315E+00, 4.7914889422856814203E-01,
CS   2        -3.7384991926068969150E-03,-2.6801520353328635310E-03,
CS   3         9.9168777670983678974E-05,-2.1877128189032726730E-06/
CS    DATA QQ/-3.1446690275135491500E+01, 8.5539563258012929600E+01,
CS   1        -6.0228002066743340583E+01, 1.3982595353892851542E+01,
CS   2        -1.1151759188741312645E+00, 3.2547697594819615062E-02,
CS   3        -5.5194330231005480228E-04/
      DATA PP/-3.9843750000000000000D-01, 2.9205384596336793945D+00,
     1        -2.4708469169133954315D+00, 4.7914889422856814203D-01,
     2        -3.7384991926068969150D-03,-2.6801520353328635310D-03,
     3         9.9168777670983678974D-05,-2.1877128189032726730D-06/
      DATA QQ/-3.1446690275135491500D+01, 8.5539563258012929600D+01,
     1        -6.0228002066743340583D+01, 1.3982595353892851542D+01,
     2        -1.1151759188741312645D+00, 3.2547697594819615062D-02,
     3        -5.5194330231005480228D-04/
C--------------------------------------------------------------------
      X = ABS(ARG)
      IF (X .LT. XSMALL) THEN
            RESULT = ONE
         ELSE IF (X .LT. ONE5) THEN
C--------------------------------------------------------------------
C  XSMALL .LE.  ABS(ARG)  .LT. 15.0
C--------------------------------------------------------------------
            XX = X * X
            SUMP = P(1)
            DO 50 I = 2, 15
              SUMP = SUMP * XX + P(I)
   50       CONTINUE
            XX = XX - TWO25
            SUMQ = ((((XX+Q(1))*XX+Q(2))*XX+Q(3))*XX+Q(4))*XX+Q(5)
            RESULT = SUMP / SUMQ
            IF (JINT .EQ. 2) RESULT = RESULT * EXP(-X)
         ELSE IF (X .GE. ONE5) THEN
            IF ((JINT .EQ. 1) .AND. (X .GT. XMAX)) THEN
                  RESULT = XINF
               ELSE
C--------------------------------------------------------------------
C  15.0  .LE.  ABS(ARG)
C--------------------------------------------------------------------
                  XX = ONE / X - REC15
                  SUMP = ((((((PP(1)*XX+PP(2))*XX+PP(3))*XX+PP(4))*XX+
     1                   PP(5))*XX+PP(6))*XX+PP(7))*XX+PP(8)
                  SUMQ = ((((((XX+QQ(1))*XX+QQ(2))*XX+QQ(3))*XX+
     1                   QQ(4))*XX+QQ(5))*XX+QQ(6))*XX+QQ(7)
                  RESULT = SUMP / SUMQ
                  IF (JINT .EQ. 2) THEN
                        RESULT = (RESULT - PP(1)) / SQRT(X)
                     ELSE
C--------------------------------------------------------------------
C  Calculation reformulated to avoid premature overflow
C--------------------------------------------------------------------
                        IF (X .LE.(XMAX-ONE5)) THEN
                              A = EXP(X)
                              B = ONE
                           ELSE
                              A = EXP(X-FORTY)
                              B = EXP40
                        END IF
                        RESULT = ((RESULT*A-PP(1)*A)/SQRT(X))*B
                  END IF
            END IF
      END IF
C--------------------------------------------------------------------
C  Return for ABS(ARG) .LT. XSMALL
C--------------------------------------------------------------------
      RETURN
C----------- Last line of CALCI0 -----------
      END
C
CS    REAL 
      DOUBLE PRECISION 
     1    FUNCTION BESI0(X)
C--------------------------------------------------------------------
C
C This long precision subprogram computes approximate values for
C   modified Bessel functions of the first kind of order zero for
C   arguments ABS(ARG) .LE. XMAX  (see comments heading CALCI0).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL  
      DOUBLE PRECISION  
     1    X, RESULT
C--------------------------------------------------------------------
      JINT=1
      CALL CALCI0(X,RESULT,JINT)
      BESI0=RESULT
      RETURN
C---------- Last line of BESI0 ----------
      END
C
CS    REAL 
      DOUBLE PRECISION 
     1    FUNCTION BESEI0(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   modified Bessel function of the first kind of order zero
C   multiplied by EXP(-ABS(X)), where EXP is the
C   exponential function, ABS is the absolute value, and X
C   is any argument.
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL  
      DOUBLE PRECISION  
     1    X, RESULT
C--------------------------------------------------------------------
      JINT=2
      CALL CALCI0(X,RESULT,JINT)
      BESEI0=RESULT
      RETURN
C---------- Last line of BESEI0 ----------
      END
C
      SUBROUTINE CALCI1(ARG,RESULT,JINT)
C--------------------------------------------------------------------
C
C This packet computes modified Bessel functioons of the first kind
C   and order one, I1(X) and EXP(-ABS(X))*I1(X), for real
C   arguments X.  It contains two function type subprograms, BESI1
C   and BESEI1, and one subroutine type subprogram, CALCI1.
C   The calling statements for the primary entries are
C
C                   Y=BESI1(X)
C   and
C                   Y=BESEI1(X)
C
C   where the entry points correspond to the functions I1(X) and
C   EXP(-ABS(X))*I1(X), respectively.  The routine CALCI1 is
C   intended for internal packet use only, all computations within
C   the packet being concentrated in this routine.  The function
C   subprograms invoke CALCI1 with the statement
C          CALL CALCI1(ARG,RESULT,JINT)
C   where the parameter usage is as follows
C
C      Function                     Parameters for CALCI1
C       Call              ARG                  RESULT          JINT
C
C     BESI1(ARG)    ABS(ARG) .LE. XMAX        I1(ARG)           1
C     BESEI1(ARG)    any real ARG        EXP(-ABS(ARG))*I1(ARG) 2
C
C   The main computation evaluates slightly modified forms of
C   minimax approximations generated by Blair and Edwards, Chalk
C   River (Atomic Energy of Canada Limited) Report AECL-4928,
C   October, 1974.  This transportable program is patterned after
C   the machine-dependent FUNPACK packet NATSI1, but cannot match
C   that version for efficiency or accuracy.  This version uses
C   rational functions that theoretically approximate I-SUB-1(X)
C   to at least 18 significant decimal digits.  The accuracy
C   achieved depends on the arithmetic system, the compiler, the
C   intrinsic functions, and proper selection of the machine-
C   dependent constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta   = Radix for the floating-point system
C   maxexp = Smallest power of beta that overflows
C   XSMALL = Positive argument such that 1.0 - X = 1.0 to
C            machine precision for all ABS(X) .LE. XSMALL.
C   XINF =   Largest positive machine number; approximately
C            beta**maxexp
C   XMAX =   Largest argument acceptable to BESI1;  Solution to
C            equation: 
C               EXP(X) * (1-3/(8*X)) / SQRT(2*PI*X) = beta**maxexp
C
C
C     Approximate values for some important machines are:
C
C                          beta       maxexp       XSMALL
C
C CRAY-1        (S.P.)       2         8191       3.55E-15
C Cyber 180/855
C   under NOS   (S.P.)       2         1070       3.55E-15
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)       2          128       2.98E-8
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)       2         1024       5.55D-17
C IBM 3033      (D.P.)      16           63       6.95D-18
C VAX           (S.P.)       2          127       2.98E-8
C VAX D-Format  (D.P.)       2          127       6.95D-18
C VAX G-Format  (D.P.)       2         1023       5.55D-17
C
C
C                               XINF          XMAX
C
C CRAY-1        (S.P.)       5.45E+2465     5682.810
C Cyber 180/855
C   under NOS   (S.P.)       1.26E+322       745.894
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)       3.40E+38         91.906
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)       1.79D+308       713.987
C IBM 3033      (D.P.)       7.23D+75        178.185
C VAX           (S.P.)       1.70D+38         91.209
C VAX D-Format  (D.P.)       1.70D+38         91.209
C VAX G-Format  (D.P.)       8.98D+307       713.293
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns the value XINF for ABS(ARG) .GT. XMAX.
C
C
C Intrinsic functions required are:
C
C     ABS, SQRT, EXP
C
C
C  Authors: W. J. Cody and L. Stoltz
C           Mathematics and Computer Science Division
C           Argonne National Laboratory
C           Argonne, IL  60439
C
C  Latest modification: May 31, 1989
C
C--------------------------------------------------------------------
      INTEGER J,JINT
CS    REAL
      DOUBLE PRECISION
     1    A,ARG,B,EXP40,FORTY,HALF,ONE,ONE5,P,PBAR,PP,Q,QQ,REC15,
     2    RESULT,SUMP,SUMQ,TWO25,X,XINF,XMAX,XSMALL,XX,ZERO
      DIMENSION P(15),PP(8),Q(5),QQ(6)
C--------------------------------------------------------------------
C  Mathematical constants
C--------------------------------------------------------------------
CS    DATA ONE/1.0E0/,ONE5/15.0E0/,EXP40/2.353852668370199854E17/,
CS   1     FORTY/40.0E0/,REC15/6.6666666666666666666E-2/,
CS   2     TWO25/225.0E0/,HALF/0.5E0/,ZERO/0.0E0/
      DATA ONE/1.0D0/,ONE5/15.0D0/,EXP40/2.353852668370199854D17/,
     1     FORTY/40.0D0/,REC15/6.6666666666666666666D-2/,
     2     TWO25/225.0D0/,HALF/0.5D0/,ZERO/0.0D0/
C--------------------------------------------------------------------
C  Machine-dependent constants
C--------------------------------------------------------------------
CS    DATA XSMALL/2.98E-8/,XINF/3.4E38/,XMAX/91.906E0/
      DATA XSMALL/5.55D-17/,XINF/1.79D308/,XMAX/713.987D0/
C--------------------------------------------------------------------
C  Coefficients for XSMALL .LE. ABS(ARG) .LT. 15.0
C--------------------------------------------------------------------
CS    DATA P/-1.9705291802535139930E-19,-6.5245515583151902910E-16,
CS   1       -1.1928788903603238754E-12,-1.4831904935994647675E-09,
CS   2       -1.3466829827635152875E-06,-9.1746443287817501309E-04,
CS   3       -4.7207090827310162436E-01,-1.8225946631657315931E+02,
CS   4       -5.1894091982308017540E+04,-1.0588550724769347106E+07,
CS   5       -1.4828267606612366099E+09,-1.3357437682275493024E+11,
CS   6       -6.9876779648010090070E+12,-1.7732037840791591320E+14,
CS   7       -1.4577180278143463643E+15/
      DATA P/-1.9705291802535139930D-19,-6.5245515583151902910D-16,
     1       -1.1928788903603238754D-12,-1.4831904935994647675D-09,
     2       -1.3466829827635152875D-06,-9.1746443287817501309D-04,
     3       -4.7207090827310162436D-01,-1.8225946631657315931D+02,
     4       -5.1894091982308017540D+04,-1.0588550724769347106D+07,
     5       -1.4828267606612366099D+09,-1.3357437682275493024D+11,
     6       -6.9876779648010090070D+12,-1.7732037840791591320D+14,
     7       -1.4577180278143463643D+15/
CS    DATA Q/-4.0076864679904189921E+03, 7.4810580356655069138E+06,
CS   1       -8.0059518998619764991E+09, 4.8544714258273622913E+12,
CS   2       -1.3218168307321442305E+15/
      DATA Q/-4.0076864679904189921D+03, 7.4810580356655069138D+06,
     1       -8.0059518998619764991D+09, 4.8544714258273622913D+12,
     2       -1.3218168307321442305D+15/
C--------------------------------------------------------------------
C  Coefficients for 15.0 .LE. ABS(ARG)
C--------------------------------------------------------------------
CS    DATA PP/-6.0437159056137600000E-02, 4.5748122901933459000E-01,
CS   1        -4.2843766903304806403E-01, 9.7356000150886612134E-02,
CS   2        -3.2457723974465568321E-03,-3.6395264712121795296E-04,
CS   3         1.6258661867440836395E-05,-3.6347578404608223492E-07/
      DATA PP/-6.0437159056137600000D-02, 4.5748122901933459000D-01,
     1        -4.2843766903304806403D-01, 9.7356000150886612134D-02,
     2        -3.2457723974465568321D-03,-3.6395264712121795296D-04,
     3         1.6258661867440836395D-05,-3.6347578404608223492D-07/
CS    DATA QQ/-3.8806586721556593450E+00, 3.2593714889036996297E+00,
CS   1        -8.5017476463217924408E-01, 7.4212010813186530069E-02,
CS   2        -2.2835624489492512649E-03, 3.7510433111922824643E-05/
      DATA QQ/-3.8806586721556593450D+00, 3.2593714889036996297D+00,
     1        -8.5017476463217924408D-01, 7.4212010813186530069D-02,
     2        -2.2835624489492512649D-03, 3.7510433111922824643D-05/
CS    DATA PBAR/3.98437500E-01/
      DATA PBAR/3.98437500D-01/
C--------------------------------------------------------------------
      X = ABS(ARG)
      IF (X .LT. XSMALL) THEN
C--------------------------------------------------------------------
C  Return for ABS(ARG) .LT. XSMALL
C--------------------------------------------------------------------
            RESULT = HALF * X
         ELSE IF (X .LT. ONE5) THEN
C--------------------------------------------------------------------
C  XSMALL .LE. ABS(ARG) .LT. 15.0
C--------------------------------------------------------------------
            XX = X * X
            SUMP = P(1)
            DO 50 J = 2, 15
               SUMP = SUMP * XX + P(J)
   50          CONTINUE
            XX = XX - TWO25
            SUMQ = ((((XX+Q(1))*XX+Q(2))*XX+Q(3))*XX+Q(4))
     1           *XX+Q(5)
            RESULT = (SUMP / SUMQ) * X
            IF (JINT .EQ. 2) RESULT = RESULT * EXP(-X)
         ELSE IF ((JINT .EQ. 1) .AND. (X .GT. XMAX)) THEN
                  RESULT = XINF
         ELSE
C--------------------------------------------------------------------
C  15.0 .LE. ABS(ARG)
C--------------------------------------------------------------------
            XX = ONE / X - REC15
            SUMP = ((((((PP(1)*XX+PP(2))*XX+PP(3))*XX+
     1           PP(4))*XX+PP(5))*XX+PP(6))*XX+PP(7))*XX+PP(8)
            SUMQ = (((((XX+QQ(1))*XX+QQ(2))*XX+QQ(3))*XX+
     1           QQ(4))*XX+QQ(5))*XX+QQ(6)
            RESULT = SUMP / SUMQ
            IF (JINT .NE. 1) THEN
                  RESULT = (RESULT + PBAR) / SQRT(X)
               ELSE
C--------------------------------------------------------------------
C  Calculation reformulated to avoid premature overflow
C--------------------------------------------------------------------
                  IF (X .GT. XMAX-ONE5) THEN
                        A = EXP(X-FORTY)
                        B = EXP40
                     ELSE
                        A = EXP(X)
                        B = ONE
                  END IF
                  RESULT = ((RESULT * A + PBAR * A) /
     1                  SQRT(X)) * B
C--------------------------------------------------------------------
C  Error return for ABS(ARG) .GT. XMAX
C--------------------------------------------------------------------
            END IF
      END IF
      IF (ARG .LT. ZERO) RESULT = -RESULT
      RETURN
C----------- Last line of CALCI1 -----------
      END
C
CS    REAL
      DOUBLE PRECISION
     1    FUNCTION BESI1(X)
C--------------------------------------------------------------------
C
C This long precision subprogram computes approximate values for
C   modified Bessel functions of the first kind of order one for
C   arguments ABS(ARG) .LE. XMAX  (see comments heading CALCI1).
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL
      DOUBLE PRECISION
     1    X, RESULT
C--------------------------------------------------------------------
      JINT=1
      CALL CALCI1(X,RESULT,JINT)
      BESI1=RESULT
      RETURN
C---------- Last line of BESI1 ----------
      END
C
CS    REAL
      DOUBLE PRECISION
     1    FUNCTION BESEI1(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   modified Bessel function of the first kind of order one
C   multiplied by EXP(-ABS(X)), where EXP is the
C   exponential function, ABS is the absolute value, and X
C   is any argument.
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL
      DOUBLE PRECISION
     1    X, RESULT
C--------------------------------------------------------------------
      JINT=2
      CALL CALCI1(X,RESULT,JINT)
      BESEI1=RESULT
      RETURN
C---------- Last line of BESEI1 ----------
      END
C
      SUBROUTINE CALCK0(ARG,RESULT,JINT)
C--------------------------------------------------------------------
C
C This packet computes modified Bessel functions of the second kind
C   and order zero, K0(X) and EXP(X)*K0(X), for real
C   arguments X.  It contains two function type subprograms, BESK0
C   and BESEK0, and one subroutine type subprogram, CALCK0.
C   the calling statements for the primary entries are
C
C                   Y=BESK0(X)
C   and
C                   Y=BESEK0(X)
C
C   where the entry points correspond to the functions K0(X) and
C   EXP(X)*K0(X), respectively.  The routine CALCK0 is
C   intended for internal packet use only, all computations within
C   the packet being concentrated in this routine.  The function
C   subprograms invoke CALCK0 with the statement
C          CALL CALCK0(ARG,RESULT,JINT)
C   where the parameter usage is as follows
C
C      Function                     Parameters for CALCK0
C       Call              ARG                  RESULT          JINT
C
C     BESK0(ARG)   0 .LT. ARG .LE. XMAX       K0(ARG)           1
C     BESEK0(ARG)     0 .LT. ARG           EXP(ARG)*K0(ARG)     2
C
C   The main computation evaluates slightly modified forms of near 
C   minimax rational approximations generated by Russon and Blair, 
C   Chalk River (Atomic Energy of Canada Limited) Report AECL-3461, 
C   1969.  This transportable program is patterned after the 
C   machine-dependent FUNPACK packet NATSK0, but cannot match that 
C   version for efficiency or accuracy.  This version uses rational
C   functions that theoretically approximate K-SUB-0(X) to at
C   least 18 significant decimal digits.  The accuracy achieved
C   depends on the arithmetic system, the compiler, the intrinsic
C   functions, and proper selection of the machine-dependent
C   constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta   = Radix for the floating-point system
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C   XSMALL = Argument below which BESK0 and BESEK0 may
C            each be represented by a constant and a log.
C            largest X such that  1.0 + X = 1.0  to machine
C            precision.
C   XINF   = Largest positive machine number; approximately
C            beta**maxexp
C   XMAX   = Largest argument acceptable to BESK0;  Solution to
C            equation: 
C               W(X) * (1-1/8X+9/128X**2) = beta**minexp
C            where  W(X) = EXP(-X)*SQRT(PI/2X)
C
C
C     Approximate values for some important machines are:
C
C
C                           beta       minexp       maxexp
C
C  CRAY-1        (S.P.)       2        -8193         8191
C  Cyber 180/185 
C    under NOS   (S.P.)       2         -975         1070
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)       2         -126          128
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)       2        -1022         1024
C  IBM 3033      (D.P.)      16          -65           63
C  VAX D-Format  (D.P.)       2         -128          127
C  VAX G-Format  (D.P.)       2        -1024         1023
C
C
C                          XSMALL       XINF         XMAX
C
C CRAY-1        (S.P.)    3.55E-15   5.45E+2465    5674.858
C Cyber 180/855
C   under NOS   (S.P.)    1.77E-15   1.26E+322      672.788
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)    5.95E-8    3.40E+38        85.337
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)    1.11D-16   1.79D+308      705.342
C IBM 3033      (D.P.)    1.11D-16   7.23D+75       177.852
C VAX D-Format  (D.P.)    6.95D-18   1.70D+38        86.715
C VAX G-Format  (D.P.)    5.55D-17   8.98D+307      706.728
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns the value XINF for ARG .LE. 0.0, and the
C  BESK0 entry returns the value 0.0 for ARG .GT. XMAX.
C
C
C  Intrinsic functions required are:
C
C     EXP, LOG, SQRT
C
C  Latest modification: March 19, 1990
C
C  Authors: W. J. Cody and Laura Stoltz
C           Mathematics and Computer Science Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C--------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL  
      DOUBLE PRECISION
     1    ARG,F,G,ONE,P,PP,Q,QQ,RESULT,SUMF,SUMG,SUMP,SUMQ,TEMP,
     2    X,XINF,XMAX,XSMALL,XX,ZERO
      DIMENSION P(6),Q(2),PP(10),QQ(10),F(4),G(3)
C--------------------------------------------------------------------
C  Mathematical constants
C--------------------------------------------------------------------
CS    DATA ONE/1.0E0/,ZERO/0.0E0/
      DATA ONE/1.0D0/,ZERO/0.0D0/
C--------------------------------------------------------------------
C  Machine-dependent constants
C--------------------------------------------------------------------
CS    DATA XSMALL/5.95E-8/,XINF/3.40E+38/,XMAX/ 85.337E0/
      DATA XSMALL/1.11D-16/,XINF/1.79D+308/,XMAX/705.342D0/
C--------------------------------------------------------------------
C
C     Coefficients for XSMALL .LE.  ARG  .LE. 1.0
C
C--------------------------------------------------------------------
CS    DATA   P/ 5.8599221412826100000E-04, 1.3166052564989571850E-01,  
CS   1          1.1999463724910714109E+01, 4.6850901201934832188E+02,
CS   2          5.9169059852270512312E+03, 2.4708152720399552679E+03/
CS    DATA   Q/-2.4994418972832303646E+02, 2.1312714303849120380E+04/ 
CS    DATA   F/-1.6414452837299064100E+00,-2.9601657892958843866E+02,
CS   1         -1.7733784684952985886E+04,-4.0320340761145482298E+05/
CS    DATA   G/-2.5064972445877992730E+02, 2.9865713163054025489E+04,
CS   1         -1.6128136304458193998E+06/
      DATA   P/ 5.8599221412826100000D-04, 1.3166052564989571850D-01,  
     1          1.1999463724910714109D+01, 4.6850901201934832188D+02,
     2          5.9169059852270512312D+03, 2.4708152720399552679D+03/
      DATA   Q/-2.4994418972832303646D+02, 2.1312714303849120380D+04/ 
      DATA   F/-1.6414452837299064100D+00,-2.9601657892958843866D+02,
     1         -1.7733784684952985886D+04,-4.0320340761145482298D+05/
      DATA   G/-2.5064972445877992730D+02, 2.9865713163054025489D+04,
     1         -1.6128136304458193998D+06/
C--------------------------------------------------------------------
C
C     Coefficients for  1.0 .LT. ARG
C
C--------------------------------------------------------------------
CS    DATA  PP/ 1.1394980557384778174E+02, 3.6832589957340267940E+03,
CS   1          3.1075408980684392399E+04, 1.0577068948034021957E+05,
CS   2          1.7398867902565686251E+05, 1.5097646353289914539E+05,
CS   3          7.1557062783764037541E+04, 1.8321525870183537725E+04,
CS   4          2.3444738764199315021E+03, 1.1600249425076035558E+02/
CS    DATA  QQ/ 2.0013443064949242491E+02, 4.4329628889746408858E+03,
CS   1          3.1474655750295278825E+04, 9.7418829762268075784E+04,
CS   2          1.5144644673520157801E+05, 1.2689839587977598727E+05,
CS   3          5.8824616785857027752E+04, 1.4847228371802360957E+04,
CS   4          1.8821890840982713696E+03, 9.2556599177304839811E+01/
      DATA  PP/ 1.1394980557384778174D+02, 3.6832589957340267940D+03,
     1          3.1075408980684392399D+04, 1.0577068948034021957D+05,
     2          1.7398867902565686251D+05, 1.5097646353289914539D+05,
     3          7.1557062783764037541D+04, 1.8321525870183537725D+04,
     4          2.3444738764199315021D+03, 1.1600249425076035558D+02/
      DATA  QQ/ 2.0013443064949242491D+02, 4.4329628889746408858D+03,
     1          3.1474655750295278825D+04, 9.7418829762268075784D+04,
     2          1.5144644673520157801D+05, 1.2689839587977598727D+05,
     3          5.8824616785857027752D+04, 1.4847228371802360957D+04,
     4          1.8821890840982713696D+03, 9.2556599177304839811D+01/
C--------------------------------------------------------------------
      X = ARG
      IF (X .GT. ZERO) THEN
            IF (X .LE. ONE) THEN
C--------------------------------------------------------------------
C     0.0 .LT.  ARG  .LE. 1.0
C--------------------------------------------------------------------
                  TEMP = LOG(X)
                  IF (X .LT. XSMALL) THEN
C--------------------------------------------------------------------
C     Return for small ARG
C--------------------------------------------------------------------
                        RESULT = P(6)/Q(2) - TEMP
                     ELSE
                        XX = X * X
                        SUMP = ((((P(1)*XX + P(2))*XX + P(3))*XX +
     1                         P(4))*XX + P(5))*XX + P(6)
                        SUMQ = (XX + Q(1))*XX + Q(2)
                        SUMF = ((F(1)*XX + F(2))*XX + F(3))*XX + F(4)
                        SUMG = ((XX + G(1))*XX + G(2))*XX + G(3)
                        RESULT = SUMP/SUMQ - XX*SUMF*TEMP/SUMG - TEMP
                        IF (JINT .EQ. 2) RESULT = RESULT * EXP(X)
                  END IF
               ELSE IF ((JINT .EQ. 1) .AND. (X .GT. XMAX)) THEN
C--------------------------------------------------------------------
C     Error return for ARG .GT. XMAX
C--------------------------------------------------------------------
                  RESULT = ZERO
               ELSE
C--------------------------------------------------------------------
C     1.0 .LT. ARG
C--------------------------------------------------------------------
                  XX = ONE / X
                  SUMP = PP(1)
                  DO 120 I = 2, 10
                     SUMP = SUMP*XX + PP(I)
  120             CONTINUE
                  SUMQ = XX
                  DO 140 I = 1, 9
                     SUMQ = (SUMQ + QQ(I))*XX
  140             CONTINUE
                  SUMQ = SUMQ + QQ(10)
                  RESULT = SUMP / SUMQ / SQRT(X)
                  IF (JINT .EQ. 1) RESULT = RESULT * EXP(-X)
            END IF
         ELSE
C--------------------------------------------------------------------
C     Error return for ARG .LE. 0.0
C--------------------------------------------------------------------
            RESULT = XINF
      END IF
C--------------------------------------------------------------------
C     Update error counts, etc.
C--------------------------------------------------------------------
      RETURN
C---------- Last line of CALCK0 ----------
      END
C
CS    REAL 
      DOUBLE PRECISION
     1    FUNCTION BESK0(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   modified Bessel function of the second kind of order zero
C   for arguments 0.0 .LT. ARG .LE. XMAX (see comments heading
C   CALCK0).
C
C  Authors: W. J. Cody and Laura Stoltz
C
C  Latest Modification: January 19, 1988
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL 
      DOUBLE PRECISION
     1    X, RESULT
C--------------------------------------------------------------------
      JINT = 1
      CALL CALCK0(X,RESULT,JINT)
      BESK0 = RESULT
      RETURN
C---------- Last line of BESK0 ----------
      END
C
CS    REAL
      DOUBLE PRECISION 
     1    FUNCTION BESEK0(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   modified Bessel function of the second kind of order zero
C   multiplied by the Exponential function, for arguments
C   0.0 .LT. ARG.
C
C  Authors: W. J. Cody and Laura Stoltz
C
C  Latest Modification: January 19, 1988
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL 
      DOUBLE PRECISION 
     1    X, RESULT
C--------------------------------------------------------------------
      JINT = 2
      CALL CALCK0(X,RESULT,JINT)
      BESEK0 = RESULT
      RETURN
C---------- Last line of BESEK0 ----------
      END
C
      SUBROUTINE CALCK1(ARG,RESULT,JINT)
C--------------------------------------------------------------------
C
C This packet computes modified Bessel functions of the second kind
C   and order one,  K1(X)  and  EXP(X)*K1(X), for real arguments X.
C   It contains two function type subprograms, BESK1  and  BESEK1,
C   and one subroutine type subprogram, CALCK1.  The calling
C   statements for the primary entries are
C
C                   Y=BESK1(X)
C   and
C                   Y=BESEK1(X)
C
C   where the entry points correspond to the functions K1(X) and
C   EXP(X)*K1(X), respectively.  The routine CALCK1 is intended
C   for internal packet use only, all computations within the
C   packet being concentrated in this routine.  The function
C   subprograms invoke CALCK1 with the statement
C          CALL CALCK1(ARG,RESULT,JINT)
C   where the parameter usage is as follows
C
C      Function                      Parameters for CALCK1
C        Call             ARG                  RESULT          JINT
C
C     BESK1(ARG)  XLEAST .LT. ARG .LT. XMAX    K1(ARG)          1
C     BESEK1(ARG)     XLEAST .LT. ARG       EXP(ARG)*K1(ARG)    2
C
C   The main computation evaluates slightly modified forms of near 
C   minimax rational approximations generated by Russon and Blair, 
C   Chalk River (Atomic Energy of Canada Limited) Report AECL-3461, 
C   1969.  This transportable program is patterned after the 
C   machine-dependent FUNPACK packet NATSK1, but cannot match that
C   version for efficiency or accuracy.  This version uses rational
C   functions that theoretically approximate K-SUB-1(X) to at
C   least 18 significant decimal digits.  The accuracy achieved
C   depends on the arithmetic system, the compiler, the intrinsic
C   functions, and proper selection of the machine-dependent
C   constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta   = Radix for the floating-point system
C   minexp = Smallest representable power of beta
C   maxexp = Smallest power of beta that overflows
C   XLEAST = Smallest acceptable argument, i.e., smallest machine
C            number X such that 1/X is machine representable.
C   XSMALL = Argument below which BESK1(X) and BESEK1(X) may
C            each be represented by 1/X.  A safe value is the
C            largest X such that  1.0 + X = 1.0  to machine
C            precision.
C   XINF   = Largest positive machine number; approximately
C            beta**maxexp
C   XMAX   = Largest argument acceptable to BESK1;  Solution to
C            equation:  
C               W(X) * (1+3/8X-15/128X**2) = beta**minexp
C            where  W(X) = EXP(-X)*SQRT(PI/2X)
C
C
C     Approximate values for some important machines are:
C
C                           beta       minexp       maxexp
C
C  CRAY-1        (S.P.)       2        -8193         8191
C  Cyber 180/185 
C    under NOS   (S.P.)       2         -975         1070
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)       2         -126          128
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)       2        -1022         1024
C  IBM 3033      (D.P.)      16          -65           63
C  VAX D-Format  (D.P.)       2         -128          127
C  VAX G-Format  (D.P.)       2        -1024         1023
C
C
C                         XLEAST     XSMALL      XINF       XMAX 
C
C CRAY-1                1.84E-2466  3.55E-15  5.45E+2465  5674.858 
C Cyber 180/855
C   under NOS   (S.P.)  3.14E-294   1.77E-15  1.26E+322    672.789
C IEEE (IBM/XT,
C   SUN, etc.)  (S.P.)  1.18E-38    5.95E-8   3.40E+38      85.343
C IEEE (IBM/XT,
C   SUN, etc.)  (D.P.)  2.23D-308   1.11D-16  1.79D+308    705.343
C IBM 3033      (D.P.)  1.39D-76    1.11D-16  7.23D+75     177.855
C VAX D-Format  (D.P.)  5.88D-39    6.95D-18  1.70D+38      86.721
C VAX G-Format  (D.P.)  1.12D-308   5.55D-17  8.98D+307    706.728
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns the value XINF for ARG .LE. 0.0 and the
C   BESK1 entry returns the value 0.0 for ARG .GT. XMAX.
C
C
C  Intrinsic functions required are:
C
C     LOG, SQRT, EXP
C
C
C  Authors: W. J. Cody and Laura Stoltz
C           Mathematics and Computer Science Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C  Latest modification: January 28, 1988
C
C--------------------------------------------------------------------
      INTEGER I,JINT
CS    REAL 
      DOUBLE PRECISION 
     1    ARG,F,G,ONE,P,PP,Q,QQ,RESULT,SUMF,SUMG,
     2    SUMP,SUMQ,X,XINF,XMAX,XLEAST,XSMALL,XX,ZERO
      DIMENSION P(5),Q(3),PP(11),QQ(9),F(5),G(3)
C--------------------------------------------------------------------
C  Mathematical constants
C--------------------------------------------------------------------
CS    DATA ONE/1.0E0/,ZERO/0.0E0/
      DATA ONE/1.0D0/,ZERO/0.0D0/
C--------------------------------------------------------------------
C  Machine-dependent constants
C--------------------------------------------------------------------
CS    DATA XLEAST/1.18E-38/,XSMALL/5.95E-8/,XINF/3.40E+38/,
CS   1     XMAX/85.343E+0/
      DATA XLEAST/2.23D-308/,XSMALL/1.11D-16/,XINF/1.79D+308/,
     1     XMAX/705.343D+0/
C--------------------------------------------------------------------
C  Coefficients for  XLEAST .LE.  ARG  .LE. 1.0
C--------------------------------------------------------------------
CS    DATA   P/ 4.8127070456878442310E-1, 9.9991373567429309922E+1,
CS   1          7.1885382604084798576E+3, 1.7733324035147015630E+5,
CS   2          7.1938920065420586101E+5/
CS    DATA   Q/-2.8143915754538725829E+2, 3.7264298672067697862E+4,
CS   1         -2.2149374878243304548E+6/
CS    DATA   F/-2.2795590826955002390E-1,-5.3103913335180275253E+1,
CS   1         -4.5051623763436087023E+3,-1.4758069205414222471E+5,
CS   2         -1.3531161492785421328E+6/
CS    DATA   G/-3.0507151578787595807E+2, 4.3117653211351080007E+4,
CS   2         -2.7062322985570842656E+6/
      DATA   P/ 4.8127070456878442310D-1, 9.9991373567429309922D+1,
     1          7.1885382604084798576D+3, 1.7733324035147015630D+5,
     2          7.1938920065420586101D+5/
      DATA   Q/-2.8143915754538725829D+2, 3.7264298672067697862D+4,
     1         -2.2149374878243304548D+6/
      DATA   F/-2.2795590826955002390D-1,-5.3103913335180275253D+1,
     1         -4.5051623763436087023D+3,-1.4758069205414222471D+5,
     2         -1.3531161492785421328D+6/
      DATA   G/-3.0507151578787595807D+2, 4.3117653211351080007D+4,
     2         -2.7062322985570842656D+6/
C--------------------------------------------------------------------
C  Coefficients for  1.0 .LT.  ARG
C--------------------------------------------------------------------
CS    DATA  PP/ 6.4257745859173138767E-2, 7.5584584631176030810E+0,
CS   1          1.3182609918569941308E+2, 8.1094256146537402173E+2,
CS   2          2.3123742209168871550E+3, 3.4540675585544584407E+3,
CS   3          2.8590657697910288226E+3, 1.3319486433183221990E+3,
CS   4          3.4122953486801312910E+2, 4.4137176114230414036E+1,
CS   5          2.2196792496874548962E+0/
CS    DATA  QQ/ 3.6001069306861518855E+1, 3.3031020088765390854E+2,
CS   1          1.2082692316002348638E+3, 2.1181000487171943810E+3,
CS   2          1.9448440788918006154E+3, 9.6929165726802648634E+2,
CS   3          2.5951223655579051357E+2, 3.4552228452758912848E+1,
CS   4          1.7710478032601086579E+0/
      DATA  PP/ 6.4257745859173138767D-2, 7.5584584631176030810D+0,
     1          1.3182609918569941308D+2, 8.1094256146537402173D+2,
     2          2.3123742209168871550D+3, 3.4540675585544584407D+3,
     3          2.8590657697910288226D+3, 1.3319486433183221990D+3,
     4          3.4122953486801312910D+2, 4.4137176114230414036D+1,
     5          2.2196792496874548962D+0/
      DATA  QQ/ 3.6001069306861518855D+1, 3.3031020088765390854D+2,
     1          1.2082692316002348638D+3, 2.1181000487171943810D+3,
     2          1.9448440788918006154D+3, 9.6929165726802648634D+2,
     3          2.5951223655579051357D+2, 3.4552228452758912848D+1,
     4          1.7710478032601086579D+0/
C--------------------------------------------------------------------
      X = ARG
      IF (X .LT. XLEAST) THEN
C--------------------------------------------------------------------
C  Error return for  ARG  .LT. XLEAST
C--------------------------------------------------------------------
            RESULT = XINF
         ELSE IF (X .LE. ONE) THEN
C--------------------------------------------------------------------
C  XLEAST .LE.  ARG  .LE. 1.0
C--------------------------------------------------------------------
            IF (X .LT. XSMALL) THEN
C--------------------------------------------------------------------
C  Return for small ARG
C--------------------------------------------------------------------
                  RESULT = ONE / X
               ELSE
                  XX = X * X
                  SUMP = ((((P(1)*XX + P(2))*XX + P(3))*XX + P(4))*XX
     1                   + P(5))*XX + Q(3)
                  SUMQ = ((XX + Q(1))*XX + Q(2))*XX + Q(3)
                  SUMF = (((F(1)*XX + F(2))*XX + F(3))*XX + F(4))*XX
     1                   + F(5)
                  SUMG = ((XX + G(1))*XX + G(2))*XX + G(3)
                  RESULT = (XX * LOG(X) * SUMF/SUMG + SUMP/SUMQ) / X
                  IF (JINT .EQ. 2) RESULT = RESULT * EXP(X)
            END IF
         ELSE IF ((JINT .EQ. 1) .AND. (X .GT. XMAX)) THEN
C--------------------------------------------------------------------
C  Error return for  ARG  .GT. XMAX
C--------------------------------------------------------------------
            RESULT = ZERO
         ELSE
C--------------------------------------------------------------------
C  1.0 .LT.  ARG
C--------------------------------------------------------------------
            XX = ONE / X
            SUMP = PP(1)
            DO 120 I = 2, 11
               SUMP = SUMP * XX + PP(I)
  120       CONTINUE
            SUMQ = XX
            DO 140 I = 1, 8
               SUMQ = (SUMQ + QQ(I)) * XX
  140       CONTINUE
            SUMQ = SUMQ + QQ(9)
            RESULT = SUMP / SUMQ / SQRT(X)
            IF (JINT .EQ. 1) RESULT = RESULT * EXP(-X)
      END IF
      RETURN
C---------- Last line of CALCK1 ----------
      END
C
CS    REAL 
      DOUBLE PRECISION 
     1    FUNCTION BESK1(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   modified Bessel function of the second kind of order one
C   for arguments  XLEAST .LE. ARG .LE. XMAX.
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL 
      DOUBLE PRECISION  
     1    X, RESULT
C--------------------------------------------------------------------
      JINT = 1
      CALL CALCK1(X,RESULT,JINT)
      BESK1 = RESULT
      RETURN
C---------- Last line of BESK1 ----------
      END
C
CS    REAL 
      DOUBLE PRECISION 
     1    FUNCTION BESEK1(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   modified Bessel function of the second kind of order one
C   multiplied by the exponential function, for arguments
C   XLEAST .LE. ARG .LE. XMAX.
C
C--------------------------------------------------------------------
      INTEGER JINT
CS    REAL
      DOUBLE PRECISION  
     1    X, RESULT
C--------------------------------------------------------------------
      JINT = 2
      CALL CALCK1(X,RESULT,JINT)
      BESEK1 = RESULT
      RETURN
C---------- Last line of BESEK1 ----------
      END
C
C     **************************************
C
C        BESSEL FUNCTION OF THE FIRST KIND
C
C     **************************************
C
      SUBROUTINE BESSJN(X,NMAX,BJN,DBJN)
C
      REAL*8  BJN(0:NMAX),DBJN(0:NMAX),X
C
      CALL BESJNV(X,NMAX,BJN,IERR)
      IF(IERR.NE.0) RETURN
C
      DBJN(0)=BJN(1)
      IF(ABS(X).EQ.0.D0)THEN
         DBJN(1)=0.5D0
         DO N=2,NMAX
            DBJN(N)=0.D0
         ENDDO
      ELSE
         DO N=1,NMAX
            DBJN(N)=BJN(N-1)-N*BJN(N)/X
         ENDDO
      ENDIF
      RETURN
      END
C
C     ******************************************************************
C
C                       ****** LAMBDA FUNCTION ******
C
C        MODIFIED BESSEL FUNCTION OF THE FIRST KIND, MULIPLIED BY EXP
C                             COMPLEX ARGUMENT
C
C     ******************************************************************
C
      SUBROUTINE LAMBDA(N,CX,CALAM)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      DIMENSION CALAM(0:N)
      DATA ONE/1.D0/
      DATA D55/1.D-55/
C
      XA=ABS(CX)
      NA=ABS(N)
C
      IF(NA.GE.30000.OR.N.LT.0.OR.XA.GE.173.D0) THEN
         CALAM(0)=(0.D0,0.D0)
         WRITE(6,*) 'XX LAMBDA : OUT OF RANGE : N,ABS(X) = ',N,XA
         RETURN
      ENDIF
      IF(XA.LE.1.D-8) THEN
         IF(XA.LE.1.D-77) THEN
            CALAM(0)=1.D0
            DO I=1,NA
               CALAM(I)=0.D0
            ENDDO
         ELSE
            CALAM(0)=EXP(-CX)
            CT1=0.5D0*CX
            T2=ONE
            CT3=ONE
            DO I=1,NA
               IF(ABS(CT3).LE.1.D-77*ABS(T2/CT1)) THEN
                  CALAM(I)=0.D0
               ELSE
                  CT3=CT3*CT1/T2
                  T2=T2+ONE
                  CALAM(I)=CT3*EXP(-CX)
               ENDIF
            ENDDO
         ENDIF
      ELSE
         CZ=2.D0/CX
         IF(XA.GE.10.D0) THEN
            L=40
         ELSEIF(XA.GE.0.1D0) THEN
            L=25
         ELSE
            L=10
         ENDIF
         RX=XA
         NM=MAX(NA,INT(RX))+L
         CT3=0.D0
         CT2=1.D-75
         CS=0.D0
         DO I=1,NM
            K=NM-I
            CT1=(K+1)*CT2*CZ+CT3
            IF(K.LE.NA) CALAM(K)=CT1
            CS=CS+CT1
            IF(ABS(CS).GT.1.D55) THEN
               CT1=CT1*D55
               CT2=CT2*D55
               CS=CS*D55
               DO J=K,NA
                  CALAM(J)=CALAM(J)*D55
               ENDDO
            ENDIF
            CT3=CT2
            CT2=CT1
         ENDDO
         CS=CS+CS-CT1
         DO J=0,NA
            CALAM(J)=CALAM(J)/CS
         ENDDO
      ENDIF
      RETURN
      END
