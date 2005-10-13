C     $Id$
C
C     ***** One-variable Brent Method *****
C
      double precision function FBRENT(f,ax,bx,tol)
      double precision ax,bx,f,tol
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result (.ge.0.)
c
c  output..
c
c  zeroin abscissa approximating a zero of  f  in the interval ax,bx
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  this is checked, and an error message is printed if this is not
c  satisfied.   zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
c  the  relative machine precision defined as the smallest representable
c  number such that  1.+macheps .gt. 1.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice-hall, inc. (1973).
c
      double precision  a,b,c,d,e,eps,fa,fb,fc,tol1,xm,p,q,r,s
      double precision  dabs
      parameter(eps=1.D-15,itmax=100)
C
      tol1 = eps+1.0d0
c
      a=ax
      b=bx
      fa=f(a)
      fb=f(b)
c     check that f(ax) and f(bx) have different signs
      if (fa .ne. 0.0d0 .and.
     &    fb .ne. 0.0d0) then
         if (fa * (fb/dabs(fb)) .gt. 0.0d0) then
            write(6,*) 'XX FBRENT: f(ax) and f(bx) have same sign.'
            write(6,'(A,1P4E12.4)') '  ax,bx,f(ax),f(bx)=',a,b,fa,fb
            return
         endif
      endif
      c=b
      fc=fb
C
      DO IT=1,ITMAX
         if ((fb*(fc/dabs(fc))).gt.0.0d0) THEN
            c=a
            fc=fa
            d=b-a
            e=d
         endif
C
         if (dabs(fc).lt.dabs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
         endif
         tol1=2.0d0*eps*dabs(b)+0.5d0*tol
         xm = 0.5d0*(c-b)
C
         if ((dabs(xm).le.tol1).or.(fb.eq.0.0d0)) then
            fbrent=b
            return
         endif
c
c see if a bisection is forced
c
         if ((dabs(e).ge.tol1).and.(dabs(fa).gt.dabs(fb))) then
            s=fb/fa
            if (a.eq.c) then
c
c linear interpolation
c
               p=2.0d0*xm*s
               q=1.0d0-s
            else
c
c inverse quadratic interpolation
c
               q=fa/fc
               r=fb/fc
               p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
               q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
            endif
            if (p.gt.0.0d0) then
               q=-q
            else
               p=-p
            endif
            s=e
            e=d
            if (((2.0d0*p).lt.(3.0d0*xm*q-dabs(tol1*q))).and.
     &           (p.lt.dabs(0.5d0*s*q))) then
               d=p/q
            else
               d=xm
               e=d
            endif
         else
            d=xm
            e=d
         endif
         a=b
         fa=fb
         if (dabs(d).gt.tol1) then
            b=b+d
         else
            if (xm.gt.0.0d0) then
               b=b+tol1
            else
               b=b-tol1
            endif
         endif
         fb=f(b)
      ENDDO
C
      WRITE(6,*) 'XX FBRENT EXCEEDING MAXIMUM ITERATIONS'
      fbrent=b
      return
      end
C
C     ***** Multi-variable Brent Method *****
C
      SUBROUTINE FBRENTN(FCN,N,X,FVEC,TOL,INFO,WA,LWA)
      INTEGER N,INFO,LWA
      DOUBLE PRECISION TOL
      DOUBLE PRECISION X(N),FVEC(N),WA(LWA)
      EXTERNAL FCN
C     **********
C
C     SUBROUTINE BRENT1
C
C     THE PURPOSE OF THIS SUBROUTINE IS TO FIND A ZERO OF
C     A SYSTEM OF N NONLINEAR EQUATIONS IN N VARIABLES BY A
C     METHOD DUE TO R. BRENT. THIS IS DONE BY USING THE
C     MORE GENERAL NONLINEAR EQUATION SOLVER BRENTM.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE BRENT1(FCN,N,X,FVEC,TOL,INFO,WA,LWA)
C
C     WHERE
C
C       FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
C         CALCULATES COMPONENTS OF THE FUNCTION. FCN SHOULD BE
C         DECLARED IN AN EXTERNAL STATEMENT IN THE USER CALLING
C         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N)
C         ----------
C         CALCULATE THE IFLAG-TH COMPONENT OF THE FUNCTION
C         AND RETURN THIS VALUE IN FVEC(IFLAG).
C         ----------
C         RETURN
C         END
C
C         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
C         THE USER WANTS TO TERMINATE EXECUTION OF BRENT1.
C         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF EQUATIONS AND VARIABLES.
C
C       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN
C         AN INITIAL ESTIMATE OF THE SOLUTION VECTOR. ON OUTPUT X
C         CONTAINS THE FINAL ESTIMATE OF THE SOLUTION VECTOR.
C
C       FVEC IS AN ARRAY OF LENGTH N. ON OUTPUT IT CONTAINS
C         THE FINAL RESIDUALS.
C
C       TOL IS A NONNEGATIVE INPUT VARIABLE. THE ALGORITHM CONVERGES
C         IF EITHER ALL THE RESIDUALS ARE AT MOST TOL IN MAGNITUDE,
C         OR IF THE ALGORITHM ESTIMATES THAT THE RELATIVE ERROR
C         BETWEEN X AND THE SOLUTION IS AT MOST TOL.
C
C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS. IF
C         THE USER HAS TERMINATED EXECUTION, INFO WILL BE SET TO
C         THE (NEGATIVE) VALUE OF IFLAG. SEE DESCRIPTION OF FCN.
C         OTHERWISE
C
C         INFO = 0   IMPROPER INPUT PARAMETERS.
C
C         INFO = 1   ALL RESIDUALS ARE AT MOST TOL IN MAGNITUDE.
C
C         INFO = 2   ALGORITHM ESTIMATES THAT THE RELATIVE ERROR
C                    BETWEEN X AND THE SOLUTION IS AT MOST TOL.
C
C         INFO = 3   CONDITIONS FOR INFO = 1 AND INFO = 2 BOTH HOLD.
C
C         INFO = 4   NUMBER OF FUNCTION EVALUATIONS HAS REACHED OR
C                    EXCEEDED 50*(N+3).
C
C         INFO = 5   APPROXIMATE JACOBIAN MATRIX IS SINGULAR.
C
C         INFO = 6   ITERATION IS NOT MAKING GOOD PROGRESS.
C
C         INFO = 7   ITERATION IS DIVERGING.
C
C         INFO = 8   ITERATION IS CONVERGING, BUT TOL IS TOO
C                    SMALL, OR THE CONVERGENCE IS VERY SLOW
C                    DUE TO A JACOBIAN SINGULAR NEAR THE OUTPUT
C                    X OR DUE TO BADLY SCALED VARIABLES.
C
C       WA IS A WORK ARRAY OF LENGTH LWA.
C
C       LWA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
C         N*(N+3).
C
C     SUBPROGRAMS REQUIRED
C
C       USER-SUPPLIED ...... FCN, BRENTM
C
C       FORTRAN-SUPPLIED ... DLOG
C
C     **********
      INTEGER I,MAXFEV,MOPT,NFEV
      DOUBLE PRECISION EMAX,FTOL,TEMP,XTOL,ZERO
      DOUBLE PRECISION DFLOAT
      DATA ZERO /0.D0/
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (N .LE. 0 .OR. TOL .LT. ZERO .OR. LWA .LT. N*(N+3)) THEN
         INFO=0
         RETURN
      ENDIF
C
C     DETERMINE AN OPTIMAL VALUE FOR MOPT.
C
      EMAX = ZERO
      DO I = 1, N
         TEMP = DLOG(DFLOAT(I+1))/DFLOAT(N+2*I+1)
         IF (TEMP .LT. EMAX) GO TO 20
         MOPT = I
         EMAX = TEMP
      ENDDO
   20 CONTINUE
C
C     CALL BRENTM.
C
      MAXFEV = 50*(N + 3)
      FTOL = TOL
      XTOL = TOL
      CALL FBRENTM(FCN,N,X,FVEC,FTOL,XTOL,MAXFEV,MOPT,
     1           INFO,NFEV,WA(3*N+1),N,WA(1),WA(N+1),WA(2*N+1))
   30  CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE BRENT1.
C
      END
C
      SUBROUTINE FBRENTM(FCN,N,X,FVEC,FTOL,XTOL,MAXFEV,MOPT,
     1                   INFO,NFEV,Q,LDQ,SIGMA,WA1,WA2)
      INTEGER N,MAXFEV,MOPT,INFO,NFEV,LDQ
      DOUBLE PRECISION FTOL,XTOL
      DOUBLE PRECISION X(N),FVEC(N),Q(LDQ,N),SIGMA(N),WA1(N),WA2(N)
C     **********
C
C     SUBROUTINE BRENTM
C
C     THE PURPOSE OF THIS SUBROUTINE IS TO FIND A ZERO TO
C     A SYSTEM OF N NONLINEAR EQUATIONS IN N VARIABLES BY A
C     METHOD DUE TO R. BRENT.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE BRENTM(FCN,N,X,FVEC,FTOL,XTOL,MAXFEV,MOPT,
C                         INFO,NFEV,Q,LDQ,SIGMA,WA1,WA2)
C
C     WHERE
C
C       FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
C         CALCULATES COMPONENTS OF THE FUNCTION. FCN SHOULD BE
C         DECLARED IN AN EXTERNAL STATEMENT IN THE USER CALLING
C         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N)
C         ----------
C         CALCULATE THE IFLAG-TH COMPONENT OF THE FUNCTION
C         AND RETURN THIS VALUE IN FVEC(IFLAG).
C         ----------
C         RETURN
C         END
C
C         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
C         THE USER WANTS TO TERMINATE EXECUTION OF BRENTM.
C         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF
C         EQUATIONS AND VARIABLES.
C
C       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN
C         AN ESTIMATE TO THE SOLUTION OF THE SYSTEM OF EQUATIONS.
C         ON OUTPUT X CONTAINS THE FINAL ESTIMATE TO THE SOLUTION
C         OF THE SYSTEM OF EQUATIONS.
C
C       FVEC IS AN ARRAY OF LENGTH N. ON OUTPUT IT CONTAINS
C         THE FINAL RESIDUALS.
C
C       FTOL IS A NONNEGATIVE INPUT VARIABLE. CONVERGENCE
C         OCCURS IF ALL RESIDUALS ARE AT MOST FTOL IN MAGNITUDE.
C
C       XTOL IS A NONNEGATIVE INPUT VARIABLE. CONVERGENCE
C         OCCURS IF THE RELATIVE ERROR BETWEEN TWO SUCCESSIVE
C         ITERATES IS AT MOST XTOL.
C
C       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
C         OCCURS IF THE NUMBER OF FUNCTION EVALUATIONS IS AT
C         LEAST MAXFEV BY THE END OF AN ITERATION. IN BRENTM,
C         A FUNCTION EVALUATION CORRESPONDS TO N CALLS TO FCN.
C
C       MOPT IS A POSITIVE INTEGER INPUT VARIABLE. MOPT SPECIFIES
C         THE NUMBER OF TIMES THAT THE APPROXIMATE JACOBIAN IS
C         USED DURING EACH ITERATION WHICH EMPLOYS ITERATIVE
C         REFINEMENT. IF MOPT IS 1, NO ITERATIVE REFINEMENT WILL
C         BE DONE. MAXIMUM EFFICIENCY IS USUALLY OBTAINED IF
C         MOPT MAXIMIZES LOG(K+1)/(N+2*K+1) FOR K = 1,...,N.
C
C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS. IF
C         THE USER HAS TERMINATED EXECUTION, INFO WILL BE SET TO
C         THE (NEGATIVE) VALUE OF IFLAG. SEE DESCRIPTION OF FCN.
C         OTHERWISE
C
C         INFO = 0   IMPROPER INPUT PARAMETERS.
C
C         INFO = 1   ALL RESIDUALS ARE AT MOST FTOL IN MAGNITUDE.
C
C         INFO = 2   RELATIVE ERROR BETWEEN TWO SUCCESSIVE ITERATES
C                    IS AT MOST XTOL.
C
C         INFO = 3   CONDITIONS FOR INFO = 1 AND INFO = 2 BOTH HOLD.
C
C         INFO = 4   NUMBER OF FUNCTION EVALUATIONS HAS REACHED OR
C                    EXCEEDED MAXFEV.
C
C         INFO = 5   APPROXIMATE JACOBIAN MATRIX IS SINGULAR.
C
C         INFO = 6   ITERATION IS NOT MAKING GOOD PROGRESS.
C
C         INFO = 7   ITERATION IS DIVERGING.
C
C         INFO = 8   ITERATION IS CONVERGING, BUT XTOL IS TOO
C                    SMALL, OR THE CONVERGENCE IS VERY SLOW
C                    DUE TO A JACOBIAN SINGULAR NEAR THE OUTPUT
C                    X OR DUE TO BADLY SCALED VARIABLES.
C
C       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
C         FUNCTION EVALUATIONS USED IN PRODUCING X. IN BRENTM,
C         A FUNCTION EVALUATION CORRESPONDS TO N CALLS TO FCN.
C
C       Q IS AN N BY N ARRAY. IF JAC DENOTES THE APPROXIMATE
C         JACOBIAN, THEN ON OUTPUT Q IS (A MULTIPLE OF) AN
C         ORTHOGONAL MATRIX SUCH THAT JAC*Q IS A LOWER TRIANGULAR
C         MATRIX. ONLY THE DIAGONAL ELEMENTS OF JAC*Q NEED
C         TO BE STORED, AND THESE CAN BE FOUND IN SIGMA.
C
C       LDQ IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY Q.
C
C       SIGMA IS A LINEAR ARRAY OF LENGTH N. ON OUTPUT SIGMA
C         CONTAINS THE DIAGONAL ELEMENTS OF THE MATRIX JAC*Q.
C         SEE DESCRIPTION OF Q.
C
C       WA1 AND WA2 ARE LINEAR WORK ARRAYS OF LENGTH N.
C
C     SUBPROGRAMS REQUIRED
C
C       USER-SUPPLIED ...... FCN
C
C       FORTRAN-SUPPLIED ... DABS,DMAX1,DSQRT,DSIGN
C
C     **********
      INTEGER I,IFLAG,J,K,M,NFCALL,NIER6,NIER7,NIER8,NSING
      LOGICAL CONV
      DOUBLE PRECISION DELTA,DIFIT,DIFIT1,EPS,EPSMCH,ETA,FKY,FKZ,
     1       FNORM,FNORM1,H,P05,SCALE,SKNORM,TEMP,XNORM,ZERO
      DATA ZERO,P05,SCALE /0.D0,5.D-2,1.D1/
C
C     WARNING.
C
C     THIS IS AN IBM CODE. TO RUN THIS CODE ON OTHER MACHINES IT
C     IS NECESSARY TO CHANGE THE VALUE OF THE MACHINE PRECISION
C     EPSMCH. THE MACHINE PRECISION IS THE SMALLEST FLOATING
C     POINT NUMBER EPSMCH SUCH THAT
C
C           1 + EPSMCH .GT. 1
C
C     IN WORKING PRECISION. IF IN DOUBT ABOUT THE VALUE OF
C     EPSMCH, THEN THE FOLLOWING PROGRAM SEGMENT DETERMINES
C     EPSMCH ON MOST MACHINES.
C
C     EPSMCH = 0.5D0
C   1 CONTINUE
C     IF (1.D0+EPSMCH .EQ. 1.D0) GO TO 2
C     EPSMCH = 0.5D0*EPSMCH
C     GO TO 1
C   2 CONTINUE
C     EPSMCH = 2.D0*EPSMCH
C
C     THE IBM DOUBLE PRECISION EPSMCH.
C
C      EPSMCH = 16.D0**(-13)
      EPSMCH = 1.D-15
C
      INFO = 0
      IFLAG = 0
      NFEV = 0
      NFCALL = 0
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (N .LE. 0 .OR. FTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR.
     1    MAXFEV .LE. 0 .OR. MOPT .LE. 0 .OR. LDQ .LT. N) GO TO 220
C
C     INITIALIZE SOME OF THE VARIABLES.
C
      NIER6 = -1
      NIER7 = -1
      NIER8 = 0
      FNORM = ZERO
      DIFIT = ZERO
      XNORM = ZERO
      DO I = 1, N
         XNORM = DMAX1(XNORM,DABS(X(I)))
      ENDDO
      EPS = DSQRT(EPSMCH)
      DELTA = SCALE*XNORM
      IF (XNORM .EQ. ZERO) DELTA = SCALE
C
C     ENTER THE PRINCIPAL ITERATION.
C
   20  CONTINUE
C
C     TO PRINT THE ITERATES, PLACE WRITE STATEMENTS
C     FOR THE VECTOR X HERE.
C
      NSING = N
      FNORM1 = FNORM
      DIFIT1 = DIFIT
      FNORM = ZERO
C
C     COMPUTE THE STEP H FOR THE DIVIDED DIFFERENCE WHICH
C     APPROXIMATES THE K-TH ROW OF THE JACOBIAN MATRIX.
C
      H = EPS*XNORM
      IF (H .EQ. ZERO) H = EPS
      DO J = 1, N
         DO I = 1, N
            Q(I,J) = ZERO
         ENDDO
         Q(J,J) = H
         WA1(J) = X(J)
      ENDDO
C
C     ENTER A SUBITERATION.
C
      DO K = 1, N
         IFLAG = K
         CALL FCN(N,WA1,FVEC,IFLAG)
         FKY = FVEC(K)
         NFCALL = NFCALL + 1
         NFEV = NFCALL/N
         IF (IFLAG .LT. 0) GO TO 230
         FNORM = DMAX1(FNORM,DABS(FKY))
C
C        COMPUTE THE K-TH ROW OF THE JACOBIAN MATRIX.
C
         DO J = K, N
            DO I = 1, N
               WA2(I) = WA1(I) + Q(I,J)
            ENDDO
            CALL FCN(N,WA2,FVEC,IFLAG)
            FKZ = FVEC(K)
            NFCALL = NFCALL + 1
            NFEV = NFCALL/N
            IF (IFLAG .LT. 0) GO TO 230
            SIGMA(J) = FKZ - FKY
         ENDDO
         FVEC(K) = FKY
C
C        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE K-TH ROW
C        OF THE JACOBIAN MATRIX TO A MULTIPLE OF THE K-TH UNIT VECTOR.
C
         ETA = ZERO
         DO I = K, N
            ETA = DMAX1(ETA,DABS(SIGMA(I)))
         ENDDO
         IF (ETA .EQ. ZERO) GO TO 150
         NSING = NSING - 1
         SKNORM = ZERO
         DO I = K, N
            SIGMA(I) = SIGMA(I)/ETA
            SKNORM = SKNORM + SIGMA(I)**2
         ENDDO
         SKNORM = DSQRT(SKNORM)
         IF (SIGMA(K) .LT. ZERO) SKNORM = -SKNORM
         SIGMA(K) = SIGMA(K) + SKNORM
C
C        APPLY THE TRANSFORMATION AND COMPUTE THE MATRIX Q.
C
         DO I = 1, N
            WA2(I) = ZERO
         ENDDO
         DO J = K, N
            TEMP = SIGMA(J)
            DO I = 1, N
               WA2(I) = WA2(I) + TEMP*Q(I,J)
            ENDDO
         ENDDO
         DO J = K, N
            TEMP = SIGMA(J)/(SKNORM*SIGMA(K))
            DO I = 1, N
               Q(I,J) = Q(I,J) - TEMP*WA2(I)
            ENDDO
         ENDDO
C
C        COMPUTE THE SUBITERATE.
C
         SIGMA(K) = SKNORM*ETA
         TEMP = FKY/SIGMA(K)
         IF (H*DABS(TEMP) .GT. DELTA) TEMP = DSIGN(DELTA/H,TEMP)
         DO I = 1, N
            WA1(I) = WA1(I) + TEMP*Q(I,K)
         ENDDO
  150    CONTINUE
      ENDDO
C
C     COMPUTE THE NORMS OF THE ITERATE AND CORRECTION VECTOR.
C
      XNORM = ZERO
      DIFIT = ZERO
      DO I = 1, N
         XNORM = DMAX1(XNORM,DABS(WA1(I)))
         DIFIT = DMAX1(DIFIT,DABS(X(I)-WA1(I)))
         X(I) = WA1(I)
      ENDDO
C
C     UPDATE THE BOUND ON THE CORRECTION VECTOR.
C
      DELTA = DMAX1(DELTA,SCALE*XNORM)
C
C     DETERMINE THE PROGRESS OF THE ITERATION.
C
      CONV = (FNORM .LT. FNORM1 .AND. DIFIT .LT. DIFIT1 .AND.
     1        NSING .EQ. 0)
      NIER6 = NIER6 + 1
      NIER7 = NIER7 + 1
      NIER8 = NIER8 + 1
      IF (CONV) NIER6 = 0
      IF (FNORM .LT. FNORM1 .OR. DIFIT .LT. DIFIT1) NIER7 = 0
      IF (DIFIT .GT. EPS*XNORM) NIER8 = 0
C
C     TESTS FOR CONVERGENCE.
C
      IF (FNORM .LE. FTOL) INFO = 1
      IF (DIFIT .LE. XTOL*XNORM .AND. CONV) INFO = 2
      IF (FNORM .LE. FTOL .AND. INFO .EQ. 2) INFO = 3
      IF (INFO .NE. 0) GO TO 230
C
C     TESTS FOR TERMINATION.
C
      IF (NFEV .GE. MAXFEV) INFO = 4
      IF (NSING .EQ. N) INFO = 5
      IF (NIER6 .EQ. 5) INFO = 6
      IF (NIER7 .EQ. 3) INFO = 7
      IF (NIER8 .EQ. 4) INFO = 8
      IF (INFO .NE. 0) GO TO 230
C
C     ITERATIVE REFINEMENT IS USED IF THE ITERATION IS CONVERGING.
C
      IF (.NOT. CONV .OR. DIFIT .GT. P05*XNORM .OR. MOPT .EQ. 1)
     1    GO TO 220
C
C     START ITERATIVE REFINEMENT.
C
      DO M = 2, MOPT
         FNORM1 = FNORM
         FNORM = ZERO
         DO K = 1, N
            IFLAG = K
            CALL FCN(N,WA1,FVEC,IFLAG)
            FKY = FVEC(K)
            NFCALL = NFCALL + 1
            NFEV = NFCALL/N
            IF (IFLAG .LT. 0) GO TO 230
            FNORM = DMAX1(FNORM,DABS(FKY))
C
C           ITERATIVE REFINEMENT IS TERMINATED IF IT DOES NOT
C           GIVE A REDUCTION OF THE RESIDUALS.
C
            IF (FNORM .GE. FNORM1) THEN
               FNORM = FNORM1
               GO TO 220
            ENDIF
            TEMP = FKY/SIGMA(K)
            DO I = 1, N
               WA1(I) = WA1(I) + TEMP*Q(I,K)
            ENDDO
         ENDDO
C
C        COMPUTE THE NORMS OF THE ITERATE AND CORRECTION VECTOR.
C
         XNORM = ZERO
         DIFIT = ZERO
         DO I = 1, N
            XNORM = DMAX1(XNORM,DABS(WA1(I)))
            DIFIT = DMAX1(DIFIT,DABS(X(I)-WA1(I)))
            X(I) = WA1(I)
         ENDDO
C
C        STOPPING CRITERIA FOR ITERATIVE REFINEMENT.
C
         IF (FNORM .LE. FTOL) INFO = 1
         IF (DIFIT .LE. XTOL*XNORM) INFO = 2
         IF (FNORM .LE. FTOL .AND. INFO .EQ. 2) INFO = 3
         IF (NFEV .GE. MAXFEV .AND. INFO .EQ. 0) INFO = 4
         IF (INFO .NE. 0) GO TO 230
      ENDDO
  220 CONTINUE
C
C     END OF THE ITERATIVE REFINEMENT.
C
      GO TO 20
C
C     TERMINATION, EITHER NORMAL OR USER IMPOSED.
C
  230 CONTINUE
      IF (IFLAG .LT. 0) INFO = IFLAG
      RETURN
C
C     LAST CARD OF SUBROUTINE BRENTM.
C
      END
