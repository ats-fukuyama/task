C     $Id$
C
C***********************************************************************
      SUBROUTINE SPLC(X, N, Y, DF, IOPT, C, NC, IER)
C***********************************************************************
C  COMPUTE THE COEFFICIENTS OF THE CUBIC SPLINE.                       *
C  PARAMETERS                                                          *
C    (1) X: 1-DIM. ARRAY FOR KNOWN POINTS                              *
C    (2) N: NUMBER OF KNOWN POINTS                                     *
C    (3) Y: 1-DIM. ARRAY FOR FUNCTION'S VALUES ON KNOWN POINTS         *
C    (4) DF: 1-DIM. ARRAY FOR DIFFERENTIALS AT END POINTS              *
C    (5) IOPT: 1-DIM. ARRAY SPECIFYING THE CONTENT OF DF               *
C    (6) C: 2-DIM. WORKING ARRAY                                       *
C    (7) NC: ROW SIZE OF THE ARRAY (C)                                 *
C    (8) IER: ERROR CODE                                               *
C  COPYRIGHT   T. OGUNI   JUNE 30 1989    VERSION 1.0                  *
C***********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION X(N), Y(N), DF(2), IOPT(2), C(NC,3), EC(4), D(2)
C
      IF (N.LT.2 .OR. NC.LT.N-1 .OR. IOPT(1).LT.1 .OR. IOPT(1).GT.3
     &    .OR. IOPT(2).LT.1 .OR. IOPT(2).GT.3) THEN
       IER = 2
       WRITE(*,*) '(SUBR. SPLC) INVALID ARGUMENT.',N,NC,IOPT(1),IOPT(2)
       RETURN
      ENDIF
      DO 5 I=1,N-1
       IF (X(I) .GE. X(I+1)) THEN
        IER = 1
        WRITE(*,*) '(SUBR. SPLC) X SHOULD SATISFY UPWARD ORDER.'
        RETURN
       ENDIF
    5 CONTINUE
      IER = 0
C  SET THE END CONDITIONS.
      II = 2
      KS = 1
      KE = MIN0(4,N)
      IDER = 1
      DO 70 I=1,2
       I1 = 2 * I - 1
       I2 = 2 * I
       IB = IOPT(I)
       GO TO (10, 20, 30), IB
   10  EC(I1) = 0.0D0
       EC(I2) = 2.0D0 * DF(I)
       GO TO 70
   20  D(I) = DF(I)
   25  IF (I .EQ. 2) II = N
       H = X(II) - X(II-1)
       EC(I1) = 1.0D0
       HY = Y(II) - Y(II-1)
       EC(I2) = 6.0D0 * (HY / H - D(I)) / H
       IF (I .EQ. 2) EC(I2) = - EC(I2)
       GO TO 70
   30  IF (I .NE. 1) THEN
        KS = MAX0(1,N-3)
        KE = N
        IDER = N
       ENDIF
       A2 = 0.0D0
       D(I) = 0.0D0
       DO 60 K=KS,KE
        IF (IDER .NE. K) THEN
         A1 = 1.0D0
         DO 50 J=KS,KE
          IF (J .NE. IDER .AND. J .NE. K) THEN
           X1 = X(IDER) - X(J)
           X2 = X(K) - X(J)
           A1 = A1 * X1 / X2
          ENDIF
   50    CONTINUE
         X3 = X(K) - X(IDER)
         D(I) = D(I) + A1 * Y(K) / X3
         A2 = A2 - 1.0D0 / X3
        ENDIF
   60  CONTINUE
       D(I) = D(I) + Y(IDER) * A2
       GO TO 25
   70 CONTINUE
C  SET THE ELEMENTS FOR THE SYMMETRIC TRIDIAGONAL EQUATION.
      IF (N .NE. 2) THEN
       H1 = X(2) - X(1)
       Y1 = Y(2) - Y(1)
       DO 80 I=2,N-1
        H2 = X(I+1) - X(I)
        Y2 = Y(I+1) - Y(I)
        HH = H1 + H2
        C(I,1) = H2 / HH
        C(I,2) = 1.0D0 - C(I,1)
        C(I,3) = 6.0D0 * (Y2 / H2 - Y1 / H1) / HH
        H1 = H2
   80   Y1 = Y2
      ENDIF
C  SOLVE THE EQUATION
      C(1,1) = - EC(1) * 0.5D0
      C(1,2) =   EC(2) * 0.5D0
      IF (N .NE. 2) THEN
       DO 100 K=2,N-1
        PIV = 2.0D0 + C(K,2) * C(K-1,1)
        C(K,1) = - C(K,1) / PIV
  100   C(K,2) = (C(K,3) - C(K,2) * C(K-1,2)) / PIV
      ENDIF
      DY1 = (EC(4) - EC(3) * C(N-1,2)) / (2.0D0 + EC(3) * C(N-1,1))
      DO 120 I=1,N-1
       K = N - I
       DY2 = C(K,1) * DY1 + C(K,2)
       H = X(K+1) - X(K)
       C(K,3) = (DY1 - DY2) / (6.0D0 * H)
       C(K,2) = 0.5D0 * DY2
       C(K,1) = (Y(K+1) - Y(K)) / H - (C(K,2) + C(K,3) * H) * H
  120  DY1 = DY2
C
      RETURN
      END
        SUBROUTINE SPLQ(X, N, Y, C, NC, A, B, S, IER)
C***********************************************************************
C  INTERPOLATION BY THE CUBIC SPLINE.                                  *
C  PARAMETERS                                                          *
C    (1) X: 1-DIM. ARRAY FOR KNOWN POINTS                              *
C    (2) N: NUMBER OF KNOWN POINTS                                     *
C    (3) Y: 1-DIM. ARRAY FOR FUNCTION'S VALUES ON KNOWN POINTS         *
C    (4) C: 2-DIM. WORKING ARRAY                                       *
C    (5) NC: ROW SIZE OF THE ARRAY (C)                                 *
C    (6) A: LOWER LIMIT OF THE INTERPOLATION INTERVAL                  *
C    (7) B: UPPER LIMIT OF THE INTERPOLATION INTERVAL                  *
C    (8) S: INTEGRAL OF THE FUNCTION ON THE INTERVAL                   *
C    (9) IER: ERROR CODE                                               *
C  COPYRIGHT   T. OGUNI    JUNE 30 1989    VERSION 1.0                 *
C***********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION X(N), Y(N), C(NC,3)
C
      IF (N .LT. 2 .OR. NC .LT. N-1) THEN
       IER = 2
       WRITE(*,*) '(SUBR. SPLQ) INVALID ARGUMENT. ',N, NC
       RETURN
      ENDIF
      IER = 0
      IND = 1
      IA = 1
      IFLA = 0
      IFLB = 0
      XX = A
      IF (A .GT. B) XX = B
   10 X1 = XX - X(IA)
      DO 20 I=IA,N-1
       IP = I
       X2 = XX - X(I+1)
       IF (X2 .LT. 0.0) GO TO 30
       IF (I .LT. N-1) X1 = X2
   20 CONTINUE
      IP = N - 1
      IF (X2 .GT. 0.0) IER = 1
   30 IF (X1 .LT. 0.0) THEN
       IER = 1
       X1 = - X1
       IF (IND .EQ. 1) IFLA = 1
       IF (IND .EQ. 2) IFLB = 1
      ENDIF
      IF (IND .NE. 2) THEN
       IND = 2
       IA = IP
       XA = X1
       XX = B
       IF (A .GT. B) XX = A
       GO TO 10
      ENDIF
   50 IB = IP
      XB = X1
C  COMPUTE THE INTEGRAL FROM A TO B.
      SA = Y(IA)+XA*(C(IA,1)/2.0D0+XA*(C(IA,2)/3.0D0+XA*C(IA,3)/4.0D0))
      SA = SA * XA
      IF (IFLA .EQ. 1) SA = - SA
      SAB = 0.0D0
      IF (IB-1 .GE. IA) THEN
       DO 60 I=IA,IB-1
        H = X(I+1) - X(I)
   60   SAB = SAB+H*(Y(I+1)+Y(I)-(C(I+1,2)+C(I,2))*H*H/6.0D0)/2.0D0
      ENDIF
      SB = Y(IB)+XB*(C(IB,1)/2.0D0+XB*(C(IB,2)/3.0D0+XB*C(IB,3)/4.0D0))
      SB = SB * XB
      IF (IFLB .EQ. 1) SB = - SB
      S = SB + SAB - SA
      IF (A .GT. B) S =  - S
C
      RETURN
      END
