! librkf.f90

MODULE librkf

  PRIVATE
  PUBLIC RKF

CONTAINS

  SUBROUTINE RKF(NEQ,FUNC,X0,XE,YH0,INIT,RELERR,ABSERR, &
                 YHN,ESTERR,NDE,IER,H,YL0,YLN,ERR, &
                 WORK)

!*********************************************************************
!     SUBROUTINE RKF NUMERICALLY INTEGRATES A SYSTEM OF NEQ          *
!     FIRST ORDER ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM        *
!             DY(I)/DX = F(X, Y(1),..., Y(NEQ)),                     *
!     BY THE RUNGE-KUTTA-FEHLBERG (4,5) FORMULA.                     *
!                                                                    *
!     PARAMETERS                                                     *
!  === INPUT ===                                                     *
!     (1) NEQ: NUMBER OF EQUATIONS TO BE INTEGRATED                  *
!     (2) FUNC: SUBROUTINE FUNC(X,Y,F) TO EVALUATE DERIVATIVES       *
!                F(I)=DY(I)/DX                                       *
!     (3) X0: INITIAL VALUE OF INDEPENDENT VARIABLE                  *
!     (4) XE: OUTPUT POINT AT WHICH THE SOLUTION IS DESIRED          *
!     (5) YH0(I) (I=1,..,NEQ): INITIAL VALUE AT X0                   *
!     (6) INIT: INDICATOR TO INITIALIZE THE CODE                     *
!          INIT=1..THE CODE WILL BE INITIALIZED(FIRST CALL).         *
!          INIT=2..THE CODE WILL NOT BE INITIALIZED(SUBSEQUENT CALL).*
!     (7) RELERR: RELATIVE LOCAL ERROR TOLERANCE                     *
!     (8) ABSERR: ABSOLUTE LOCAL ERROR TOLERANCE                     *
!  === OUTPUT ===                                                    *
!     (9) YHN(I) (I=1,..,NEQ): APPROXIMATE SOLUTION AT XE            *
!    (10) ESTERR(I) (I=1,..,NEQ): ESTIMATE OF THE GLOBAL ERROR IN    *
!          THE APPROXIMATE SOLUTION YHN(I)                           *
!    (11) NDE: NUMBER OF DERIVATIVE EVALUATIONS                      *
!    (12) IER: INDICATOR FOR STATUS OF INTEGRATION                   *
!          IER=0..INTEGRATION REACHED XE. INDICATES SUCCESSFUL       *
!             RETURN AND IS THE NORMAL MODE FOR CONTINUING THE       *
!             INTEGRATION.                                           *
!          IER=10000..INTEGRATION WAS NOT COMPLETED BECAUSE TOO MANY *
!             DERIVATIVES EVALUATIONS WERE NEEDED.  IF THE USER WANTS*
!             TO CONTINUE THE INTEGRATION, HE JUST CALLS RKF AGAIN.  *
!          IER=20000..INTEGRATION WAS NOT COMPLETED BECAUSE ERROR    *
!             TOLERANCE WAS INAPPROPRIATE(=IT WAS ZERO).  MUST USE   *
!             NON-ZERO ABSERR TO CONTINUE THE INTEGRATION.           *
!          IER=30000..INTEGRATION WAS NOT COMPLETED BECAUSE REQUESTED*
!             ACCURACY COULD NOT BE ACHIEVED USING SMALLEST STEPSIZE.*
!             MUST INCREASE THE ERROR TOLERANCE TO CONTINUE THE      *
!             INTEGRATION.                                           *
!  === OTHERS ===                                                    *
!    (13) H: VARIABLE TO HOLD INFORMATION INTERNAL TO RKF WHICH IS   *
!           NECESSARY FOR SUBSEQUENT CALLS.                          *
!    (14) YL0(), YLN(): ARRAY (SIZE=NEQ) TO HOLD INFORMATION INTERNAL*
!           TO RKF WHICH IS NECESSARY FOR SUBSEQUENT CALLS.          *
!    (15) ERR(): ARRAY (SIZE=NEQ) TO BE USED INSIDE RKF              *
!    (16) WORK(): TWO-DIMENTIONAL ARRAY (SIZE=(NEQ,11)) TO BE        *
!                 USED INSIDE RKF                                    *
!    COPYRIGHT: M. SUGIHARA, NOVEMBER 15, 1989    V. 1               *
!*********************************************************************
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NEQ
    REAL(dp),INTENT(INOUT):: X0,YH0(NEQ)
    REAL(dp),INTENT(IN):: XE,RELERR,ABSERR
    REAL(dp),INTENT(OUT):: YHN(NEQ),ESTERR(NEQ)
    INTEGER,INTENT(INOUT):: INIT
    INTEGER,INTENT(OUT):: NDE,IER
    REAL(dp),INTENT(INOUT):: H,YL0(NEQ),YLN(NEQ),ERR(NEQ),WORK(NEQ,11)
    EXTERNAL FUNC
!    INTEGER,PARAMETER:: ITEMAX=100000
    INTEGER,PARAMETER:: ITEMAX=1000
    REAl(dp),PARAMETER:: EPSMIN=1.0D-15
    INTEGER:: I,ITER
    REAL(dp):: TOLDY,TOL,HMIN,ERRET,ET,X

    IF (INIT .EQ. 1) THEN

!   -------------- INITIALIZATION (FIRST CALL)--------------------------

       DO I = 1,NEQ
          YL0(I) = YH0(I)
       END DO
       
!     -------- SET INITIAL STEP SIZE ---------

       CALL FUNC(X0,YL0,WORK(1:NEQ,1))
       TOLDY = 10.0D+74
       DO I = 1,NEQ
          IF (WORK(I,1) .NE. 0.0D0) THEN
             TOL = RELERR * DABS(YL0(I)) + ABSERR
             TOLDY = MIN(TOLDY, ABS(TOL / WORK(I,1)))
          END IF
       END DO
       HMIN = EPSMIN * (XE - X0)
       IF (TOLDY .EQ. 10.0D+74) THEN
          H = HMIN
       ELSE
          H = MIN((XE - X0) / 100.0D0, TOLDY ** 0.2D0)
       END IF
!    ----------------------------------------
       INIT = 2
       NDE = 1
    ELSE
       NDE = 0
    END IF
!  ------------------------------------------------------------------
    IER = 0
    X = X0

!********************* MAIN ITERATION *********************************

    DO ITER = 1,ITEMAX
       CALL FLSTEP(NEQ,FUNC,X,H,YL0,YLN,ERR,WORK)
       NDE = NDE + 6
!   ----------- COMPUTE ERROR TOLERANCES ---------------------
       ERRET = 0.0D0
       DO I = 1,NEQ
          ET = RELERR * 0.5D0 * (DABS(YL0(I)) + DABS(YLN(I))) + ABSERR
          IF (ET .EQ. 0.0D0) THEN
             GO TO 20000
          ELSE
             ERRET = MAX(ERRET, ABS(ERR(I)) / ET)
          END IF
       END DO
       IF (ERRET .GT. 1.0D0) THEN
          
!   ----------- UNSUCCESSFUL STEP ------------------------------

          IF (ERRET .GE. 59049.0D0) THEN
             H = 0.1D0 * H
          ELSE
             H = 0.9D0 * H / ERRET ** 0.2D0
          END IF
          IF (H .LE. EPSMIN) GO TO 30000
       ELSE IF ((X + H) .LT. XE) THEN
          
!   ----------- SUCCESSFUL STEP (X+H < THE END POINT) -----------

          CALL FHSTEP(NEQ,FUNC,X,H,YH0,YHN,ERR,WORK)
          NDE = NDE + 12
          X = X + H
          X0 = X
          DO I = 1,NEQ
             YL0(I) = YLN(I)
             YH0(I) = YHN(I)
          END DO
          
!       ------- CHOOSE NEXT STEP --------

          IF (ERRET .LE. 1.889568D-4) THEN
             H = 5.0D0 * H
          ELSE
             H = 0.9D0 * H / ERRET ** 0.2D0
          END IF
!       ---------------------------------
       ELSE
          
!   ------------ SUCCESSFUL STEP (X+H > =  THE END POINT) ----------

          H = XE - X
          CALL FLSTEP(NEQ,FUNC,X,H,YL0,YLN,ERR,WORK)
          CALL FHSTEP(NEQ,FUNC,X,H,YH0,YHN,ERR,WORK)
          NDE = NDE + 18
          
!       -------- ESTIMATE GLOBAL ERROR -------

          DO I = 1,NEQ
             ESTERR(I) = (YLN(I) - YHN(I)) / 31.0D0
          END DO
!       --------------------------------------
          X0 = XE
          DO I = 1,NEQ
             YL0(I) = YLN(I)
             YH0(I) = YHN(I)
          END DO
          RETURN
       ENDIF
    END DO
    
!**********************************************************************

    IER = 10000
    WRITE(6,10001) X
10001 FORMAT(' ','(SUBR.-RKF) TROUBLE(TOO MANY ITERATIONS))', &
           ' AT X = ',1PE15.7)
    RETURN
20000 CONTINUE
    IER = 20000
    WRITE(6,20001) X
20001 FORMAT(' ','(SUBR.-RKF) TROUBLE(INAPPROPRIATE ERROR TOLERANCE)', &
           ' AT X = ',1PE15.7)
    RETURN
30000 CONTINUE
    IER = 30000
    WRITE(6,30001) X
30001 FORMAT(' ','(SUBR.-RKF) TROUBLE(TOO SMALL STEP SIZE)', &
           ' AT X = ',1PE15.7)
    RETURN
  END SUBROUTINE RKF

  SUBROUTINE FLSTEP(NEQ,FUNC,X,H,Y0,YN,ERR,WORK)
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NEQ
    REAL(dp),INTENT(IN):: X,H,Y0(NEQ)
    REAL(dp),INTENT(OUT):: YN(NEQ)
    REAL(dp),INTENT(INOUT):: ERR(NEQ),WORK(NEQ,11)
    EXTERNAL FUNC
    INTEGER:: LINDEX

    LINDEX = 1
    CALL FEHL(LINDEX,NEQ,FUNC,X,H,Y0,YN,ERR,WORK(1,1),WORK(1,2), &
              WORK(1,3),WORK(1,4),WORK(1,5),WORK(1,6),WORK(1,7), &
              WORK(1,8),WORK(1,9),WORK(1,10),WORK(1,11))
    RETURN
  END SUBROUTINE FLSTEP

! --- internal one-step routine ---
    
  SUBROUTINE FHSTEP(NEQ,FUNC,X,H,Y0,YN,ERR,WORK)
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NEQ
    REAL(dp),INTENT(IN):: X,H
    REAL(dp),INTENT(OUT):: YN(:)
    REAL(dp),INTENT(INOUT):: Y0(NEQ),ERR(NEQ),WORK(NEQ,11)
    EXTERNAL FUNC
    INTEGER:: LINDEX,I
    REAL(dp):: X1,H1

    X1 = X
    H1 = 0.5D0 * H
    LINDEX = 0
    CALL FEHL(LINDEX,NEQ,FUNC,X1,H1,Y0,YN,ERR,WORK(1,1),WORK(1,2), &
              WORK(1,3),WORK(1,4),WORK(1,5),WORK(1,6),WORK(1,7), &
              WORK(1,8),WORK(1,9),WORK(1,10),WORK(1,11))
    X1 = X1 + H1
    DO I = 1,NEQ
       Y0(I) = YN(I)
    END DO
    CALL FEHL(LINDEX,NEQ,FUNC,X1,H1,Y0,YN,ERR,WORK(1,1),WORK(1,2), &
              WORK(1,3),WORK(1,4),WORK(1,5),WORK(1,6),WORK(1,7), &
              WORK(1,8),WORK(1,9),WORK(1,10),WORK(1,11))
    RETURN
  END SUBROUTINE FHSTEP

! ---       
      
  SUBROUTINE FEHL(LINDEX,NEQ,FUNC,X,H,Y0,YN,ERR,AK1,AK2,AK3,AK4, &
                  AK5,AK6,W2,W3,W4,W5,W6)
    USE task_kinds,ONLY: dp
    IMPLICIT NONE
    INTEGER,INTENT(IN):: LINDEX,NEQ
    REAL(dp),INTENT(IN):: X,H,Y0(NEQ)
    REAL(dp),INTENT(OUT):: YN(NEQ)
    REAL(dp),INTENT(OUT):: ERR(NEQ)
    REAL(dp),INTENT(OUT):: AK1(NEQ),AK2(NEQ),AK3(NEQ)
    REAL(dp),INTENT(OUT):: AK4(NEQ),AK5(NEQ),AK6(NEQ)
    REAL(dp),INTENT(OUT):: W2(NEQ),W3(NEQ),W4(NEQ),W5(NEQ),W6(NEQ)
    INTEGER:: I
!    EXTERNAL FUNC
    INTERFACE
       SUBROUTINE FUNC(X,Y,F)
         USE task_kinds,ONLY: dp
         IMPLICIT NONE
         REAL(dp),INTENT(IN):: X,Y(*)
         REAL(dp),INTENT(OUT):: F(*)
       END SUBROUTINE FUNC
    END INTERFACE

    REAL(dp),PARAMETER:: ONE = 1.D0, TWO = 2.D0, THR = 3.D0, TWL = 12.D0
    REAL(dp),PARAMETER:: &
         AL2 = ONE / 4.D0, AL3 = THR / 8.D0, AL4 = TWL / 13.D0, &
         AL5 = ONE,        AL6 = ONE / 2.D0
    REAL(dp),PARAMETER:: &
         B21 = ONE / 4.D0, B31 = THR / 32.D0, B32 = 9.0D0 / 32.D0, &
         B41 = 1932.0D0 / 2197.D0, B42 = -7200.0D0 / 2197.D0,  &
         B43 = 7296.0D0 / 2197.D0, B51 = 439.0D0 / 216.D0, B52 = -8.D0, &
         B53 = 3680.0D0 / 513.D0,  B54 = -845.0D0 / 4104.D0, &
         B61 = -8.0D0 / 27.D0, B62 = 2.D0, B63 = -3544.0D0 / 2565.D0, &
         B64 = 1859.0D0 / 4104.D0, B65 = -11.0D0 / 40.D0
    REAL(dp),PARAMETER:: &
         GA1 = 16.0D0 / 135.D0, GA3 = 6656.0D0 / 12825.D0, &
         GA4 = 28561.0D0 / 56430.D0, &
         GA5 = -9.0D0 / 50.D0, GA6 = TWO / 55.D0
    REAL(dp),PARAMETER:: &
         DA1 = ONE / 360.D0, DA3 = -128.0D0 / 4275.D0, &
         DA4 = -2197.0D0 / 75240.D0, DA5 = ONE / 50.D0, DA6 = TWO / 55.D0

    CALL FUNC(X,Y0,AK1)
    DO I = 1,NEQ
       W2(I) = Y0(I) + H * B21 * AK1(I)
    END DO
    CALL FUNC(X + AL2 * H,W2,AK2)
    DO I = 1,NEQ
       W3(I) = Y0(I) + H * (B31 * AK1(I)+B32 * AK2(I))
    END DO
    CALL FUNC(X+AL3 * H,W3,AK3)
    DO I = 1,NEQ
       W4(I) = Y0(I) + H * (B41 * AK1(I) + B42 * AK2(I) + B43 * AK3(I))
    END DO
    CALL FUNC(X + AL4 * H,W4,AK4)
    DO I = 1,NEQ
       W5(I) = Y0(I) + H * (B51 * AK1(I) + B52 * AK2(I) &
                          + B53 * AK3(I) + B54 * AK4(I))
    END DO
    CALL FUNC(X + AL5 * H,W5,AK5)
    DO I = 1,NEQ
       W6(I) = Y0(I) + H * (B61 * AK1(I) + B62 * AK2(I) + B63 * AK3(I) &
                          + B64 * AK4(I) + B65 * AK5(I))
    END DO
    CALL FUNC(X + AL6 * H,W6,AK6)
    DO I = 1,NEQ
       YN(I) = Y0(I) + H * (GA1 * AK1(I) + GA3 * AK3(I) + GA4 * AK4(I) &
                          + GA5 * AK5(I) + GA6 * AK6(I))
    END DO
    IF (LINDEX .EQ. 1) THEN
       DO I = 1,NEQ 
          ERR(I) = H * (DA1 * AK1(I) + DA3 * AK3(I) + DA4 * AK4(I) &
                      + DA5 * AK5(I) + DA6 * AK6(I))
       END DO
    ENDIF
    RETURN
  END SUBROUTINE FEHL
END MODULE librkf
