C     $Id$
C
C************************************************************************
C
      SUBROUTINE RKF(NEQ,FUNC,X0,XE,YH0,INIT,RELERR,ABSERR,
     &               YHN,ESTERR,NDE,IER,H,YL0,YLN,ERR,
     &               WORK)
C
C*********************************************************************
C     SUBROUTINE RKF NUMERICALLY INTEGRATES A SYSTEM OF NEQ          *
C     FIRST ORDER ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM        *
C             DY(I)/DX = F(X, Y(1),..., Y(NEQ)),                     *
C     BY THE RUNGE-KUTTA-FEHLBERG (4,5) FORMULA.                     *
C                                                                    *
C     PARAMETERS                                                     *
C  === INPUT ===                                                     *
C     (1) NEQ: NUMBER OF EQUATIONS TO BE INTEGRATED                  *
C     (2) FUNC: SUBROUTINE FUNC(X,Y,F) TO EVALUATE DERIVATIVES       *
C                F(I)=DY(I)/DX                                       *
C     (3) X0: INITIAL VALUE OF INDEPENDENT VARIABLE                  *
C     (4) XE: OUTPUT POINT AT WHICH THE SOLUTION IS DESIRED          *
C     (5) YH0(I) (I=1,..,NEQ): INITIAL VALUE AT X0                   *
C     (6) INIT: INDICATOR TO INITIALIZE THE CODE                     *
C          INIT=1..THE CODE WILL BE INITIALIZED(FIRST CALL).         *
C          INIT=2..THE CODE WILL NOT BE INITIALIZED(SUBSEQUENT CALL).*
C     (7) RELERR: RELATIVE LOCAL ERROR TOLERANCE                     *
C     (8) ABSERR: ABSOLUTE LOCAL ERROR TOLERANCE                     *
C  === OUTPUT ===                                                    *
C     (9) YHN(I) (I=1,..,NEQ): APPROXIMATE SOLUTION AT XE            *
C    (10) ESTERR(I) (I=1,..,NEQ): ESTIMATE OF THE GLOBAL ERROR IN    *
C          THE APPROXIMATE SOLUTION YHN(I)                           *
C    (11) NDE: NUMBER OF DERIVATIVE EVALUATIONS                      *
C    (12) IER: INDICATOR FOR STATUS OF INTEGRATION                   *
C          IER=0..INTEGRATION REACHED XE. INDICATES SUCCESSFUL       *
C             RETURN AND IS THE NORMAL MODE FOR CONTINUING THE       *
C             INTEGRATION.                                           *
C          IER=10000..INTEGRATION WAS NOT COMPLETED BECAUSE TOO MANY *
C             DERIVATIVES EVALUATIONS WERE NEEDED.  IF THE USER WANTS*
C             TO CONTINUE THE INTEGRATION, HE JUST CALLS RKF AGAIN.  *
C          IER=20000..INTEGRATION WAS NOT COMPLETED BECAUSE ERROR    *
C             TOLERANCE WAS INAPPROPRIATE(=IT WAS ZERO).  MUST USE   *
C             NON-ZERO ABSERR TO CONTINUE THE INTEGRATION.           *
C          IER=30000..INTEGRATION WAS NOT COMPLETED BECAUSE REQUESTED*
C             ACCURACY COULD NOT BE ACHIEVED USING SMALLEST STEPSIZE.*
C             MUST INCREASE THE ERROR TOLERANCE TO CONTINUE THE      *
C             INTEGRATION.                                           *
C  === OTHERS ===                                                    *
C    (13) H: VARIABLE TO HOLD INFORMATION INTERNAL TO RKF WHICH IS   *
C           NECESSARY FOR SUBSEQUENT CALLS.                          *
C    (14) YL0(), YLN(): ARRAY (SIZE=NEQ) TO HOLD INFORMATION INTERNAL*
C           TO RKF WHICH IS NECESSARY FOR SUBSEQUENT CALLS.          *
C    (15) ERR(): ARRAY (SIZE=NEQ) TO BE USED INSIDE RKF              *
C    (16) WORK(): TWO-DIMENTIONAL ARRAY (SIZE=(NEQ,11)) TO BE        *
C                 USED INSIDE RKF                                    *
C    COPYRIGHT: M. SUGIHARA, NOVEMBER 15, 1989    V. 1               *
C*********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       EXTERNAL FUNC
       DIMENSION YH0(NEQ),YHN(NEQ),YL0(NEQ),YLN(NEQ),ESTERR(NEQ),
     &           ERR(NEQ),WORK(NEQ,11)
       DATA ITEMAX / 100000 /
       DATA EPSMIN / 1.0D-15 /
C
      IF (INIT .EQ. 1) THEN
C   -------------- INITIALIZATION (FIRST CALL)--------------------------
C      WRITE(6,*) NEQ,YH0(1),YH0(2),YH0(3) 
       DO 10 I = 1,NEQ
        YL0(I) = YH0(I)
   10  CONTINUE
C     -------- SET INITIAL STEP SIZE ---------
       CALL FUNC(X0,YL0,WORK(1,1))
       TOLDY = 10.0D+74
       DO 20 I = 1,NEQ
        IF (WORK(I,1) .NE. 0.0D0) THEN
         TOL = RELERR * DABS(YL0(I)) + ABSERR
         TOLDY = MIN(TOLDY, ABS(TOL / WORK(I,1)))
        END IF
   20  CONTINUE
       HMIN = EPSMIN * (XE - X0)
       IF (TOLDY .EQ. 10.0D+74) THEN
        H = HMIN
       ELSE
        H = MIN((XE - X0) / 100.0D0, TOLDY ** 0.2D0)
       END IF
C    ----------------------------------------
       INIT = 2
       NDE = 1
      ELSE
       NDE = 0
      END IF
C  ------------------------------------------------------------------
      IER = 0
      X = X0
C
C********************* MAIN ITERATION *********************************
      DO 30 ITER = 1,ITEMAX
       CALL FLSTEP(NEQ,FUNC,X,H,YL0,YLN,ERR,WORK)
       NDE = NDE + 6
C   ----------- COMPUTE ERROR TOLERANCES ---------------------
       ERRET = 0.0D0
       DO 40 I = 1,NEQ
        ET = RELERR * 0.5D0 * (DABS(YL0(I)) + DABS(YLN(I))) + ABSERR
        IF (ET .EQ. 0.0D0) THEN
         GO TO 20000
        ELSE
         ERRET = MAX(ERRET, ABS(ERR(I)) / ET)
        END IF
   40  CONTINUE
       IF (ERRET .GT. 1.0D0) THEN
C   ----------- UNSUCCESSFUL STEP ------------------------------
        IF (ERRET .GE. 59049.0D0) THEN
         H = 0.1D0 * H
        ELSE
         H = 0.9D0 * H / ERRET ** 0.2D0
        END IF
         IF (H .LE. EPSMIN) GO TO 30000
       ELSE IF ((X + H) .LT. XE) THEN
C   ----------- SUCCESSFUL STEP (X+H < THE END POINT) -----------
        CALL FHSTEP(NEQ,FUNC,X,H,YH0,YHN,ERR,WORK)
        NDE = NDE + 12
        X = X + H
        X0 = X
        DO 50 I = 1,NEQ
         YL0(I) = YLN(I)
         YH0(I) = YHN(I)
   50   CONTINUE
C       ------- CHOOSE NEXT STEP --------
        IF (ERRET .LE. 1.889568D-4) THEN
         H = 5.0D0 * H
        ELSE
         H = 0.9D0 * H / ERRET ** 0.2D0
        END IF
C       ---------------------------------
       ELSE
C   ------------ SUCCESSFUL STEP (X+H > =  THE END POINT) ----------
        H = XE - X
        CALL FLSTEP(NEQ,FUNC,X,H,YL0,YLN,ERR,WORK)
        CALL FHSTEP(NEQ,FUNC,X,H,YH0,YHN,ERR,WORK)
        NDE = NDE + 18
C       -------- ESTIMATE GLOBAL ERROR -------
        DO 60 I = 1,NEQ
         ESTERR(I) = (YLN(I) - YHN(I)) / 31.0D0
   60   CONTINUE
C       --------------------------------------
        X0 = XE
        DO 70 I = 1,NEQ
         YL0(I) = YLN(I)
         YH0(I) = YHN(I)
   70   CONTINUE
        RETURN
       ENDIF
   30 CONTINUE
C**********************************************************************
C
      IER = 10000
      WRITE( * ,10001) X
10001 FORMAT(' ','(SUBR.-RKF) TROUBLE(TOO MANY ITERATIONS))',
     &           ' AT X = ',1PE15.7)
      RETURN
20000 CONTINUE
      IER = 20000
      WRITE( * ,20001) X
20001 FORMAT(' ','(SUBR.-RKF) TROUBLE(INAPPROPRIATE ERROR TOLERANCE)',
     &           ' AT X = ',1PE15.7)
      RETURN
30000 CONTINUE
      IER = 30000
      WRITE( * ,30001) X
30001 FORMAT(' ','(SUBR.-RKF) TROUBLE(TOO SMALL STEP SIZE)',
     &           ' AT X = ',1PE15.7)
      RETURN
      END
C
      SUBROUTINE FLSTEP(NEQ,FUNC,X,H,Y0,YN,ERR,WORK)
       IMPLICIT REAL*8(A-H,O-Z)
       EXTERNAL FUNC
       DIMENSION Y0(NEQ),YN(NEQ),ERR(NEQ),WORK(NEQ,11)
      LINDEX = 1
      CALL FEHL(LINDEX,NEQ,FUNC,X,H,Y0,YN,ERR,WORK(1,1),WORK(1,2),
     &          WORK(1,3),WORK(1,4),WORK(1,5),WORK(1,6),WORK(1,7),
     &          WORK(1,8),WORK(1,9),WORK(1,10),WORK(1,11))
      RETURN
      END
C
      SUBROUTINE FHSTEP(NEQ,FUNC,X,H,Y0,YN,ERR,WORK)
       IMPLICIT REAL*8(A-H,O-Z)
       EXTERNAL FUNC
       DIMENSION Y0(NEQ),YN(NEQ),ERR(NEQ),WORK(NEQ,11)
      X1 = X
      H1 = 0.5D0 * H
      LINDEX = 0
      CALL FEHL(LINDEX,NEQ,FUNC,X1,H1,Y0,YN,ERR,WORK(1,1),WORK(1,2),
     &          WORK(1,3),WORK(1,4),WORK(1,5),WORK(1,6),WORK(1,7),
     &          WORK(1,8),WORK(1,9),WORK(1,10),WORK(1,11))
      X1 = X1 + H1
      DO 10 I = 1,NEQ
       Y0(I) = YN(I)
   10 CONTINUE
      CALL FEHL(LINDEX,NEQ,FUNC,X1,H1,Y0,YN,ERR,WORK(1,1),WORK(1,2),
     &          WORK(1,3),WORK(1,4),WORK(1,5),WORK(1,6),WORK(1,7),
     &          WORK(1,8),WORK(1,9),WORK(1,10),WORK(1,11))
      RETURN
      END
C
      SUBROUTINE FEHL(LINDEX,NEQ,FUNC,X,H,Y0,YN,ERR,AK1,AK2,AK3,AK4,
     &                AK5,AK6,W2,W3,W4,W5,W6)
       IMPLICIT REAL*8(A-H,O-Z)
       EXTERNAL FUNC
       DIMENSION Y0(NEQ),YN(NEQ),AK1(NEQ),AK2(NEQ),AK3(NEQ),AK4(NEQ),
     &           AK5(NEQ),AK6(NEQ),ERR(NEQ),W2(NEQ),W3(NEQ),W4(NEQ),
     &           W5(NEQ),W6(NEQ)
       PARAMETER(ONE = 1, TWO = 2, THR = 3, TWL = 12)
       PARAMETER(AL2 = ONE / 4, AL3 = THR / 8, AL4 = TWL / 13,
     &           AL5 = ONE, AL6 = ONE / 2)
       PARAMETER(B21 = ONE / 4, B31 = THR / 32, B32 = 9.0D0 / 32,
     &           B41 = 1932.0D0 / 2197, B42 = -7200.0D0 / 2197,
     &           B43 = 7296.0D0 / 2197,
     &           B51 = 439.0D0 / 216, B52 = -8, B53 = 3680.0D0 / 513,
     &           B54 = -845.0D0 / 4104,
     &           B61 = -8.0D0 / 27, B62 = 2, B63 = -3544.0D0 / 2565,
     &           B64 = 1859.0D0 / 4104, B65 = -11.0D0 / 40)
       PARAMETER(GA1 = 16.0D0 / 135, GA3 = 6656.0D0 / 12825,
     &           GA4 = 28561.0D0 / 56430,
     &           GA5 = -9.0D0 / 50, GA6 = TWO / 55)
       PARAMETER(DA1 = ONE / 360, DA3 = -128.0D0 / 4275,
     &           DA4 = -2197.0D0 / 75240,
     &           DA5 = ONE / 50, DA6 = TWO / 55)
      CALL FUNC(X,Y0,AK1)
      DO 10 I = 1,NEQ
       W2(I) = Y0(I) + H * B21 * AK1(I)
   10 CONTINUE
      CALL FUNC(X + AL2 * H,W2,AK2)
      DO 20 I = 1,NEQ
       W3(I) = Y0(I) + H * (B31 * AK1(I)+B32 * AK2(I))
   20 CONTINUE
      CALL FUNC(X+AL3 * H,W3,AK3)
      DO 30 I = 1,NEQ
       W4(I) = Y0(I) + H * (B41 * AK1(I) + B42 * AK2(I) + B43 * AK3(I))
   30 CONTINUE
      CALL FUNC(X + AL4 * H,W4,AK4)
      DO 40 I = 1,NEQ
       W5(I) = Y0(I) + H * (B51 * AK1(I) + B52 * AK2(I)
     &                 + B53 * AK3(I) + B54 * AK4(I))
   40 CONTINUE
      CALL FUNC(X + AL5 * H,W5,AK5)
      DO 50 I = 1,NEQ
       W6(I) = Y0(I) + H * (B61 * AK1(I) + B62 * AK2(I) + B63 * AK3(I)
     &                 + B64 * AK4(I) + B65 * AK5(I))
   50 CONTINUE
      CALL FUNC(X + AL6 * H,W6,AK6)
      DO 60 I = 1,NEQ
       YN(I) = Y0(I) + H * (GA1 * AK1(I) + GA3 * AK3(I) + GA4 * AK4(I)
     &                 + GA5 * AK5(I) + GA6 * AK6(I))
   60 CONTINUE
      IF (LINDEX .EQ. 1) THEN
       DO 70 I = 1,NEQ
        ERR(I) = H * (DA1 * AK1(I) + DA3 * AK3(I) + DA4 * AK4(I)
     &              + DA5 * AK5(I) + DA6 * AK6(I))
   70  CONTINUE
       RETURN
      ELSE
       RETURN
      ENDIF
      END
