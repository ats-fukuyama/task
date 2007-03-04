!     ***** One-variable Brent Method *****

      real(8) function FBRENT(f,ax,bx,tol)

      implicit none
      real(8), intent(in):: ax,bx,tol
      real(8), external  :: f

!      a zero of the function  f(x)  is computed in the interval ax,bx .

!  input..

!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function subprogram which evaluates f(x) for any x in
!         the interval  ax,bx
!  tol    desired length of the interval of uncertainty of the
!         final result (.ge.0.)

!  output..

!  zeroin abscissa approximating a zero of  f  in the interval ax,bx

!      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  this is checked, and an error message is printed if this is not
!  satisfied.   zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*abs(x)+tol, where macheps  is
!  the  relative machine precision defined as the smallest representable
!  number such that  1.+macheps .gt. 1.
!      this function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice-hall, inc. (1973).

      real(8) :: a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s
      integer(4) :: it
      real(8),    parameter :: eps=1.D-15
      integer(4), parameter :: itmax=100

      tol1 = eps+1.0d0

      a=ax
      b=bx
      fa=f(a)
      fb=f(b)
!     check that f(ax) and f(bx) have different signs
      if (fa .ne. 0.0d0 .and. &
     &    fb .ne. 0.0d0) then
         if (fa * (fb/abs(fb)) .gt. 0.0d0) then
            write(6,*) 'XX FBRENT: f(ax) and f(bx) have same sign.'
            write(6,'(A,1P4E12.4)') '  ax,bx,f(ax),f(bx)=',a,b,fa,fb
            return
         endif
      endif
      c=b
      fc=fb

      DO IT=1,ITMAX
         if ((fb*(fc/abs(fc))).gt.0.0d0) THEN
            c=a
            fc=fa
            d=b-a
            e=d
         endif

         if (abs(fc).lt.abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
         endif
         tol1=2.0d0*eps*abs(b)+0.5d0*tol
         xm = 0.5d0*(c-b)

         if ((abs(xm).le.tol1).or.(fb.eq.0.0d0)) then
            fbrent=b
            return
         endif

! see if a bisection is forced

         if ((abs(e).ge.tol1).and.(abs(fa).gt.abs(fb))) then
            s=fb/fa
            if (a.eq.c) then

! linear interpolation

               p=2.0d0*xm*s
               q=1.0d0-s
            else

! inverse quadratic interpolation

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
            if (((2.0d0*p).lt.(3.0d0*xm*q-abs(tol1*q))).and. &
     &           (p.lt.abs(0.5d0*s*q))) then
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
         if (abs(d).gt.tol1) then
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

      WRITE(6,*) 'XX FBRENT EXCEEDING MAXIMUM ITERATIONS'
      fbrent=b
      return
      end function FBRENT

!     ***** Multi-variable Brent Method *****

      SUBROUTINE FBRENTN(FCN,N,X,FVEC,TOL,INFO,WA,LWA)

      IMPLICIT NONE
      INTEGER(4), INTENT(IN) ::  N,LWA
      INTEGER(4), INTENT(OUT)::  INFO
      REAL(8), INTENT(IN) :: TOL
      REAL(8), DIMENSION(N),   INTENT(IN) :: X
      REAL(8), DIMENSION(N),   INTENT(OUT):: FVEC
      REAL(8), DIMENSION(LWA), INTENT(INOUT):: WA
      EXTERNAL FCN
!     **********

!     SUBROUTINE BRENT1

!     THE PURPOSE OF THIS SUBROUTINE IS TO FIND A ZERO OF
!     A SYSTEM OF N NONLINEAR EQUATIONS IN N VARIABLES BY A
!     METHOD DUE TO R. BRENT. THIS IS DONE BY USING THE
!     MORE GENERAL NONLINEAR EQUATION SOLVER BRENTM.

!     THE SUBROUTINE STATEMENT IS

!       SUBROUTINE BRENT1(FCN,N,X,FVEC,TOL,INFO,WA,LWA)

!     WHERE

!       FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
!         CALCULATES COMPONENTS OF THE FUNCTION. FCN SHOULD BE
!         DECLARED IN AN EXTERNAL STATEMENT IN THE USER CALLING
!         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.

!         SUBROUTINE FCN(N,X,FVEC,IFLAG)
!         INTEGER N,IFLAG
!         DOUBLE PRECISION X(N),FVEC(N)
!         ----------
!         CALCULATE THE IFLAG-TH COMPONENT OF THE FUNCTION
!         AND RETURN THIS VALUE IN FVEC(IFLAG).
!         ----------
!         RETURN
!         END

!         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
!         THE USER WANTS TO TERMINATE EXECUTION OF BRENT1.
!         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.

!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF EQUATIONS AND VARIABLES.

!       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN
!         AN INITIAL ESTIMATE OF THE SOLUTION VECTOR. ON OUTPUT X
!         CONTAINS THE FINAL ESTIMATE OF THE SOLUTION VECTOR.

!       FVEC IS AN ARRAY OF LENGTH N. ON OUTPUT IT CONTAINS
!         THE FINAL RESIDUALS.

!       TOL IS A NONNEGATIVE INPUT VARIABLE. THE ALGORITHM CONVERGES
!         IF EITHER ALL THE RESIDUALS ARE AT MOST TOL IN MAGNITUDE,
!         OR IF THE ALGORITHM ESTIMATES THAT THE RELATIVE ERROR
!         BETWEEN X AND THE SOLUTION IS AT MOST TOL.

!       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS. IF
!         THE USER HAS TERMINATED EXECUTION, INFO WILL BE SET TO
!         THE (NEGATIVE) VALUE OF IFLAG. SEE DESCRIPTION OF FCN.
!         OTHERWISE

!         INFO = 0   IMPROPER INPUT PARAMETERS.

!         INFO = 1   ALL RESIDUALS ARE AT MOST TOL IN MAGNITUDE.

!         INFO = 2   ALGORITHM ESTIMATES THAT THE RELATIVE ERROR
!                    BETWEEN X AND THE SOLUTION IS AT MOST TOL.

!         INFO = 3   CONDITIONS FOR INFO = 1 AND INFO = 2 BOTH HOLD.

!         INFO = 4   NUMBER OF FUNCTION EVALUATIONS HAS REACHED OR
!                    EXCEEDED 50*(N+3).

!         INFO = 5   APPROXIMATE JACOBIAN MATRIX IS SINGULAR.

!         INFO = 6   ITERATION IS NOT MAKING GOOD PROGRESS.

!         INFO = 7   ITERATION IS DIVERGING.

!         INFO = 8   ITERATION IS CONVERGING, BUT TOL IS TOO
!                    SMALL, OR THE CONVERGENCE IS VERY SLOW
!                    DUE TO A JACOBIAN SINGULAR NEAR THE OUTPUT
!                    X OR DUE TO BADLY SCALED VARIABLES.

!       WA IS A WORK ARRAY OF LENGTH LWA.

!       LWA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
!         N*(N+3).

!     SUBPROGRAMS REQUIRED

!       USER-SUPPLIED ...... FCN, BRENTM

!       FORTRAN-SUPPLIED ... DLOG

!     **********
      INTEGER(4) :: I,MAXFEV,MOPT,NFEV
      REAL(8)    :: EMAX,FTOL,TEMP,XTOL
!      REAL(8)    :: DFLOAT
      REAL(8)    :: ZERO = 0.D0

!     CHECK THE INPUT PARAMETERS FOR ERRORS.

      IF (N .LE. 0 .OR. TOL .LT. ZERO .OR. LWA .LT. N*(N+3)) THEN
         INFO=0
         RETURN
      ENDIF

!     DETERMINE AN OPTIMAL VALUE FOR MOPT.

      EMAX = ZERO
      DO I = 1, N
         TEMP = DLOG(DFLOAT(I+1))/DFLOAT(N+2*I+1)
         IF (TEMP .LT. EMAX) GO TO 20
         MOPT = I
         EMAX = TEMP
      ENDDO
   20 CONTINUE

!     CALL BRENTM.

      MAXFEV = 50*(N + 3)
      FTOL = TOL
      XTOL = TOL
      CALL FBRENTM(FCN,N,X,FVEC,FTOL,XTOL,MAXFEV,MOPT, &
     &           INFO,NFEV,WA(3*N+1),N,WA(1),WA(N+1),WA(2*N+1))
   30  CONTINUE
      RETURN

!     LAST CARD OF SUBROUTINE BRENT1.

      END SUBROUTINE FBRENTN

      SUBROUTINE FBRENTM(FCN,N,X,FVEC,FTOL,XTOL,MAXFEV,MOPT, &
     &                   INFO,NFEV,Q,LDQ,SIGMA,WA1,WA2)

      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: N,MAXFEV,MOPT,LDQ
      INTEGER(4), INTENT(OUT):: INFO, NFEV
      REAL(8),    INTENT(IN) :: FTOL,XTOL
      REAL(8),DIMENSION(N),    INTENT(OUT):: FVEC,SIGMA
      REAL(8),DIMENSION(N),  INTENT(INOUT):: X,WA1,WA2
      REAL(8),DIMENSION(LDQ,N),INTENT(OUT):: Q
      EXTERNAL :: FCN
!     **********

!     SUBROUTINE BRENTM

!     THE PURPOSE OF THIS SUBROUTINE IS TO FIND A ZERO TO
!     A SYSTEM OF N NONLINEAR EQUATIONS IN N VARIABLES BY A
!     METHOD DUE TO R. BRENT.

!     THE SUBROUTINE STATEMENT IS

!       SUBROUTINE BRENTM(FCN,N,X,FVEC,FTOL,XTOL,MAXFEV,MOPT,
!                         INFO,NFEV,Q,LDQ,SIGMA,WA1,WA2)

!     WHERE

!       FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
!         CALCULATES COMPONENTS OF THE FUNCTION. FCN SHOULD BE
!         DECLARED IN AN EXTERNAL STATEMENT IN THE USER CALLING
!         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.

!         SUBROUTINE FCN(N,X,FVEC,IFLAG)
!         INTEGER N,IFLAG
!         DOUBLE PRECISION X(N),FVEC(N)
!         ----------
!         CALCULATE THE IFLAG-TH COMPONENT OF THE FUNCTION
!         AND RETURN THIS VALUE IN FVEC(IFLAG).
!         ----------
!         RETURN
!         END

!         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
!         THE USER WANTS TO TERMINATE EXECUTION OF BRENTM.
!         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.

!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF
!         EQUATIONS AND VARIABLES.

!       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN
!         AN ESTIMATE TO THE SOLUTION OF THE SYSTEM OF EQUATIONS.
!         ON OUTPUT X CONTAINS THE FINAL ESTIMATE TO THE SOLUTION
!         OF THE SYSTEM OF EQUATIONS.

!       FVEC IS AN ARRAY OF LENGTH N. ON OUTPUT IT CONTAINS
!         THE FINAL RESIDUALS.

!       FTOL IS A NONNEGATIVE INPUT VARIABLE. CONVERGENCE
!         OCCURS IF ALL RESIDUALS ARE AT MOST FTOL IN MAGNITUDE.

!       XTOL IS A NONNEGATIVE INPUT VARIABLE. CONVERGENCE
!         OCCURS IF THE RELATIVE ERROR BETWEEN TWO SUCCESSIVE
!         ITERATES IS AT MOST XTOL.

!       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
!         OCCURS IF THE NUMBER OF FUNCTION EVALUATIONS IS AT
!         LEAST MAXFEV BY THE END OF AN ITERATION. IN BRENTM,
!         A FUNCTION EVALUATION CORRESPONDS TO N CALLS TO FCN.

!       MOPT IS A POSITIVE INTEGER INPUT VARIABLE. MOPT SPECIFIES
!         THE NUMBER OF TIMES THAT THE APPROXIMATE JACOBIAN IS
!         USED DURING EACH ITERATION WHICH EMPLOYS ITERATIVE
!         REFINEMENT. IF MOPT IS 1, NO ITERATIVE REFINEMENT WILL
!         BE DONE. MAXIMUM EFFICIENCY IS USUALLY OBTAINED IF
!         MOPT MAXIMIZES LOG(K+1)/(N+2*K+1) FOR K = 1,...,N.

!       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS. IF
!         THE USER HAS TERMINATED EXECUTION, INFO WILL BE SET TO
!         THE (NEGATIVE) VALUE OF IFLAG. SEE DESCRIPTION OF FCN.
!         OTHERWISE

!         INFO = 0   IMPROPER INPUT PARAMETERS.

!         INFO = 1   ALL RESIDUALS ARE AT MOST FTOL IN MAGNITUDE.

!         INFO = 2   RELATIVE ERROR BETWEEN TWO SUCCESSIVE ITERATES
!                    IS AT MOST XTOL.

!         INFO = 3   CONDITIONS FOR INFO = 1 AND INFO = 2 BOTH HOLD.

!         INFO = 4   NUMBER OF FUNCTION EVALUATIONS HAS REACHED OR
!                    EXCEEDED MAXFEV.

!         INFO = 5   APPROXIMATE JACOBIAN MATRIX IS SINGULAR.

!         INFO = 6   ITERATION IS NOT MAKING GOOD PROGRESS.

!         INFO = 7   ITERATION IS DIVERGING.

!         INFO = 8   ITERATION IS CONVERGING, BUT XTOL IS TOO
!                    SMALL, OR THE CONVERGENCE IS VERY SLOW
!                    DUE TO A JACOBIAN SINGULAR NEAR THE OUTPUT
!                    X OR DUE TO BADLY SCALED VARIABLES.

!       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
!         FUNCTION EVALUATIONS USED IN PRODUCING X. IN BRENTM,
!         A FUNCTION EVALUATION CORRESPONDS TO N CALLS TO FCN.

!       Q IS AN N BY N ARRAY. IF JAC DENOTES THE APPROXIMATE
!         JACOBIAN, THEN ON OUTPUT Q IS (A MULTIPLE OF) AN
!         ORTHOGONAL MATRIX SUCH THAT JAC*Q IS A LOWER TRIANGULAR
!         MATRIX. ONLY THE DIAGONAL ELEMENTS OF JAC*Q NEED
!         TO BE STORED, AND THESE CAN BE FOUND IN SIGMA.

!       LDQ IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
!         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY Q.

!       SIGMA IS A LINEAR ARRAY OF LENGTH N. ON OUTPUT SIGMA
!         CONTAINS THE DIAGONAL ELEMENTS OF THE MATRIX JAC*Q.
!         SEE DESCRIPTION OF Q.

!       WA1 AND WA2 ARE LINEAR WORK ARRAYS OF LENGTH N.

!     SUBPROGRAMS REQUIRED

!       USER-SUPPLIED ...... FCN

!       FORTRAN-SUPPLIED ... ABS,DMAX1,DSQRT,DSIGN

!     **********
      INTEGER(4) :: I,IFLAG,J,K,M,NFCALL,NIER6,NIER7,NIER8,NSING
      LOGICAL    ::CONV
      REAL(8)    :: DELTA,DIFIT,DIFIT1,EPS,EPSMCH,ETA,FKY,FKZ, &
     &              FNORM,FNORM1,H,SKNORM,TEMP,XNORM
      REAL(8)    :: ZERO = 0.D0
      REAL(8)    :: P05  = 5.D-2
      REAL(8)    :: SCALE= 1.D1


!     WARNING.

!     THIS IS AN IBM CODE. TO RUN THIS CODE ON OTHER MACHINES IT
!     IS NECESSARY TO CHANGE THE VALUE OF THE MACHINE PRECISION
!     EPSMCH. THE MACHINE PRECISION IS THE SMALLEST FLOATING
!     POINT NUMBER EPSMCH SUCH THAT

!           1 + EPSMCH .GT. 1

!     IN WORKING PRECISION. IF IN DOUBT ABOUT THE VALUE OF
!     EPSMCH, THEN THE FOLLOWING PROGRAM SEGMENT DETERMINES
!     EPSMCH ON MOST MACHINES.

!     EPSMCH = 0.5D0
!   1 CONTINUE
!     IF (1.D0+EPSMCH .EQ. 1.D0) GO TO 2
!     EPSMCH = 0.5D0*EPSMCH
!     GO TO 1
!   2 CONTINUE
!     EPSMCH = 2.D0*EPSMCH

!     THE IBM DOUBLE PRECISION EPSMCH.

!      EPSMCH = 16.D0**(-13)
      EPSMCH = 1.D-15

      INFO = 0
      IFLAG = 0
      NFEV = 0
      NFCALL = 0

!     CHECK THE INPUT PARAMETERS FOR ERRORS.

      IF (N .LE. 0 .OR. FTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR. &
     &    MAXFEV .LE. 0 .OR. MOPT .LE. 0 .OR. LDQ .LT. N) GO TO 220

!     INITIALIZE SOME OF THE VARIABLES.

      NIER6 = -1
      NIER7 = -1
      NIER8 = 0
      FNORM = ZERO
      DIFIT = ZERO
      XNORM = ZERO
      DO I = 1, N
         XNORM = DMAX1(XNORM,ABS(X(I)))
      ENDDO
      EPS = DSQRT(EPSMCH)
      DELTA = SCALE*XNORM
      IF (XNORM .EQ. ZERO) DELTA = SCALE

!     ENTER THE PRINCIPAL ITERATION.

   20  CONTINUE

!     TO PRINT THE ITERATES, PLACE WRITE STATEMENTS
!     FOR THE VECTOR X HERE.

      NSING = N
      FNORM1 = FNORM
      DIFIT1 = DIFIT
      FNORM = ZERO

!     COMPUTE THE STEP H FOR THE DIVIDED DIFFERENCE WHICH
!     APPROXIMATES THE K-TH ROW OF THE JACOBIAN MATRIX.

      H = EPS*XNORM
      IF (H .EQ. ZERO) H = EPS
      DO J = 1, N
         DO I = 1, N
            Q(I,J) = ZERO
         ENDDO
         Q(J,J) = H
         WA1(J) = X(J)
      ENDDO

!     ENTER A SUBITERATION.

      DO K = 1, N
         IFLAG = K
         CALL FCN(N,WA1,FVEC,IFLAG)
         FKY = FVEC(K)
         NFCALL = NFCALL + 1
         NFEV = NFCALL/N
         IF (IFLAG .LT. 0) GO TO 230
         FNORM = DMAX1(FNORM,ABS(FKY))

!        COMPUTE THE K-TH ROW OF THE JACOBIAN MATRIX.

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

!        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE K-TH ROW
!        OF THE JACOBIAN MATRIX TO A MULTIPLE OF THE K-TH UNIT VECTOR.

         ETA = ZERO
         DO I = K, N
            ETA = DMAX1(ETA,ABS(SIGMA(I)))
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

!        APPLY THE TRANSFORMATION AND COMPUTE THE MATRIX Q.

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

!        COMPUTE THE SUBITERATE.

         SIGMA(K) = SKNORM*ETA
         TEMP = FKY/SIGMA(K)
         IF (H*ABS(TEMP) .GT. DELTA) TEMP = DSIGN(DELTA/H,TEMP)
         DO I = 1, N
            WA1(I) = WA1(I) + TEMP*Q(I,K)
         ENDDO
  150    CONTINUE
      ENDDO

!     COMPUTE THE NORMS OF THE ITERATE AND CORRECTION VECTOR.

      XNORM = ZERO
      DIFIT = ZERO
      DO I = 1, N
         XNORM = DMAX1(XNORM,ABS(WA1(I)))
         DIFIT = DMAX1(DIFIT,ABS(X(I)-WA1(I)))
         X(I) = WA1(I)
      ENDDO

!     UPDATE THE BOUND ON THE CORRECTION VECTOR.

      DELTA = DMAX1(DELTA,SCALE*XNORM)

!     DETERMINE THE PROGRESS OF THE ITERATION.

      CONV = (FNORM .LT. FNORM1 .AND. DIFIT .LT. DIFIT1 .AND. &
     &        NSING .EQ. 0)
      NIER6 = NIER6 + 1
      NIER7 = NIER7 + 1
      NIER8 = NIER8 + 1
      IF (CONV) NIER6 = 0
      IF (FNORM .LT. FNORM1 .OR. DIFIT .LT. DIFIT1) NIER7 = 0
      IF (DIFIT .GT. EPS*XNORM) NIER8 = 0

!     TESTS FOR CONVERGENCE.

      IF (FNORM .LE. FTOL) INFO = 1
      IF (DIFIT .LE. XTOL*XNORM .AND. CONV) INFO = 2
      IF (FNORM .LE. FTOL .AND. INFO .EQ. 2) INFO = 3
      IF (INFO .NE. 0) GO TO 230

!     TESTS FOR TERMINATION.

      IF (NFEV .GE. MAXFEV) INFO = 4
      IF (NSING .EQ. N) INFO = 5
      IF (NIER6 .EQ. 5) INFO = 6
      IF (NIER7 .EQ. 3) INFO = 7
      IF (NIER8 .EQ. 4) INFO = 8
      IF (INFO .NE. 0) GO TO 230

!     ITERATIVE REFINEMENT IS USED IF THE ITERATION IS CONVERGING.

      IF (.NOT. CONV .OR. DIFIT .GT. P05*XNORM .OR. MOPT .EQ. 1) &
     &    GO TO 220

!     START ITERATIVE REFINEMENT.

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
            FNORM = DMAX1(FNORM,ABS(FKY))

!           ITERATIVE REFINEMENT IS TERMINATED IF IT DOES NOT
!           GIVE A REDUCTION OF THE RESIDUALS.

            IF (FNORM .GE. FNORM1) THEN
               FNORM = FNORM1
               GO TO 220
            ENDIF
            TEMP = FKY/SIGMA(K)
            DO I = 1, N
               WA1(I) = WA1(I) + TEMP*Q(I,K)
            ENDDO
         ENDDO

!        COMPUTE THE NORMS OF THE ITERATE AND CORRECTION VECTOR.

         XNORM = ZERO
         DIFIT = ZERO
         DO I = 1, N
            XNORM = DMAX1(XNORM,ABS(WA1(I)))
            DIFIT = DMAX1(DIFIT,ABS(X(I)-WA1(I)))
            X(I) = WA1(I)
         ENDDO

!        STOPPING CRITERIA FOR ITERATIVE REFINEMENT.

         IF (FNORM .LE. FTOL) INFO = 1
         IF (DIFIT .LE. XTOL*XNORM) INFO = 2
         IF (FNORM .LE. FTOL .AND. INFO .EQ. 2) INFO = 3
         IF (NFEV .GE. MAXFEV .AND. INFO .EQ. 0) INFO = 4
         IF (INFO .NE. 0) GO TO 230
      ENDDO
  220 CONTINUE

!     END OF THE ITERATIVE REFINEMENT.

      GO TO 20

!     TERMINATION, EITHER NORMAL OR USER IMPOSED.

  230 CONTINUE
      IF (IFLAG .LT. 0) INFO = IFLAG
      RETURN

!     LAST CARD OF SUBROUTINE BRENTM.

      END SUBROUTINE FBRENTM
