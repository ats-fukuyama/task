C     $Id$
      SUBROUTINE RK(NEQ,FUNC,X0,XE,N,Y0,YN,WORK)
C*********************************************************************
C     SUBROUTINE RK NUMERICALLY INTEGRATES A SYSTEM OF NEQ           *
C     FIRST ORDER ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM        *
C             DY(I)/DX = F(X, Y(1),..., Y(NEQ)),                     *
C     BY THE CLASSICAL RUNGE-KUTTA FORMULA.                          *
C                                                                    *
C     PARAMETERS                                                     *
C  === INPUT ===                                                     *
C     (1) NEQ: NUMBER OF EQUATIONS TO BE INTEGRATED                  *
C     (2) FUNC: SUBROUTINE FUNC(X,Y,F) TO EVALUATE DERIVATIVES       *
C                F(I)=DY(I)/DX                                       *
C     (3) X0: INITIAL VALUE OF INDEPENDENT VARIABLE                  *
C     (4) XE: OUTPUT POINT AT WHICH THE SOLUTION IS DESIRED          *
C     (5) N: NUMBER OF DIVISIONS                                     *
C        THE INTERVAL (X0, XE) IS DIVIDED INTO N SUBINTERVALS        *
C        WITH THE LENGTH (XE-X0)/N AND IN EACH SUBINTERVAL           *
C        THE CLASSICAL RUNGE-KUTTA FORMULA IS USED.                  *
C     (6) Y0(I) (I=1,..,NEQ): INITIAL VALUE AT X0                    *
C  === OUTPUT ===                                                    *
C     (7) YN(I) (I=1,..,NEQ): APPROXIMATE SOLUTION AT XE             *
C  === OTHER ===                                                     *
C     (8) WORK(): TWO-DIMENTIONAL ARRAY (SIZE=(NEQ,2)) TO BE         *
C                 USED INSIDE RK                                     *
C     COPYRIGHT: M. SUGIHARA, NOVEMBER 15, 1989, V. 1                *
C*********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       EXTERNAL FUNC
       DIMENSION Y0(NEQ),YN(NEQ),WORK(NEQ,2)
      H = (XE - X0) / N
      DO 10 I = 1,N
       CALL RKSTEP(NEQ,FUNC,X0,H,Y0,YN,WORK(1,1),WORK(1,2))
       X0 = X0 + H
       DO 20 J = 1,NEQ
        Y0(J) = YN(J)
   20  CONTINUE
   10 CONTINUE
      X0 = XE
      RETURN
      END
C
      SUBROUTINE RKSTEP(NEQ,FUNC,X,H,Y0,YN,AK,W)
       IMPLICIT REAL * 8(A-H,O-Z)
       PARAMETER(A2 = 0.5D+0, A3 = A2)
       PARAMETER(B2 = 0.5D+0, B3 = B2)
       PARAMETER(C1 = 1.0D+0 / 6, C2 = 1.0D+0 / 3, C3 = C2, C4 = C1)
       DIMENSION Y0(NEQ),YN(NEQ),AK(NEQ),W(NEQ)
       EXTERNAL FUNC
      CALL FUNC(X,Y0,AK)
      DO 10 I = 1,NEQ
       YN(I) = Y0(I) + H * C1 * AK(I)
   10 CONTINUE
      DO 20 I = 1,NEQ
       W(I) = Y0(I) + H * B2 * AK(I)
   20 CONTINUE
      CALL FUNC(X + A2 * H,W,AK)
      DO 30 I = 1,NEQ
       YN(I) = YN(I) + H * C2 * AK(I)
   30 CONTINUE
      DO 40 I = 1,NEQ
       W(I) = Y0(I) + H * B3 * AK(I)
   40 CONTINUE
      CALL FUNC(X + A3 * H,W,AK)
      DO 50 I = 1,NEQ
       YN(I) = YN(I) + H * C3 * AK(I)
   50 CONTINUE
      DO 60 I = 1,NEQ
       W(I) = Y0(I) + H * AK(I)
   60 CONTINUE
      CALL FUNC(X + H,W,AK)
      DO 70 I = 1,NEQ
       YN(I) = YN(I) + H * C4 * AK(I)
   70 CONTINUE
      RETURN
      END
