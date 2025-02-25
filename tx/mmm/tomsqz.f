      SUBROUTINE TOMSQZ(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI,IFAIL)
C-----------------------------------------------------------------------
C TOMSQZ  written by P. Strand 27-apr-98,       elfps@elmagn.chalmers.se
C----------------------------------------------------------------------- 
C CQZHES, CQZVEC, and CQZVAL: Fortran subroutines implementing the QZ 
C algorithm for solving the generalized eigenvalue problem for complex 
C matrices. (See B.S.C Garbow, ACM TOMS 4 (1978) pp. 404-410.).
C-----------------------------------------------------------------------
C 
C ON INPUT THE GENERALIZED EIGENVALUE PROBLEM IS DEFINED THROUGH THE 
C COMPLEX MATRICES
C
C       A = cmplx (AR, AI)  AND   B = cmplx (BR, BI)
C
C WHERE LEADING DIMENSION N IS AS DEFINED IN THE CALLING ROUTINE AND 
C WHERE NA IS THE ROW  RANGE IN THE CURRENT PROBLEM. THE EIGENVALUE 
C PROBLEM IS THEN DEFINED THROUGH 
C
C       A x = w B x
C
C WHERE  THE COMPLEX EIGENVECTORS 
C
C       x = cmplx (ZVR, ZVI) 
C
C TOGETHER WITH THE COMPLEX EIGENVALUE
C
C        w = cmplx(alfr, alfi)/beta
C
C IS OUTPUT FROM THE ROUTINE
C
C IFAIL WILL BE NONZERO IF CONVERGENCE HAS NOT BEEN REACH WITHIN 50 
C ITERATIONS
C-----------------------------------------------------------------------
C DECLARATIONS FOR INPUT VARIABLES
C-----------------------------------------------------------------------

      INTEGER N, NA
      REAL AR(N,NA),AI(N,NA),BR(N,NA),BI(N,NA)

C-----------------------------------------------------------------------
C DECALRATIONS FOR OUTPUT VARIABLES
C-----------------------------------------------------------------------

      REAL ALFR(N),ALFI(N),BETA(N)
      REAL ZVR(N,NA), ZVI(N,NA)

C-----------------------------------------------------------------------
C LOCAL VARIABLES
C-----------------------------------------------------------------------

      LOGICAL WANTX
      REAL EPS1

C-----------------------------------------------------------------------
C START OF ACTUAL CODING
C-----------------------------------------------------------------------

      WANTX = .TRUE.
      EPS1  = -0.0

      CALL CQZHES(N,NA,AR,AI,BR,BI,WANTX,ZVR,ZVI)
      CALL CQZVAL(N,NA,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,WANTX,
     &            ZVR,ZVI,IFAIL)
      CALL CQZVEC(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI)
      RETURN
      END
                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE CQZHES(NM,N,AR,AI,BR,BI,MATZ,ZR,ZI)                    
C
      INTEGER I,J,K,L,N,K1,LB,L1,NM,NK1,NM1
      REAL AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ZR(NM,N),ZI(NM,N)
      REAL R,S,T,TI,U1,U2,XI,XR,YI,YR,RHO,U1I
      REAL SQRT,CABS,ABS
      LOGICAL MATZ
      COMPLEX CMPLX
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FIRST STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM WITH REAL (AND NON-
C     NEGATIVE) SUBDIAGONAL ELEMENTS AND THE OTHER TO UPPER TRIANGULAR
C     FORM USING UNITARY TRANSFORMATIONS.  IT IS USUALLY FOLLOWED BY
C     CQZVAL  AND POSSIBLY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO, AND THE
C          SUBDIAGONAL ELEMENTS HAVE BEEN MADE REAL (AND NON-NEGATIVE),
C
C        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        Z=(ZR,ZI) CONTAINS THE PRODUCT OF THE RIGHT HAND
C          TRANSFORMATIONS IF MATZ HAS BEEN SET TO .TRUE.
C          OTHERWISE, Z IS NOT REFERENCED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** INITIALIZE Z **********

      ZERO = 0.0
      ZONE = 1.0
      IF (.NOT. MATZ) GO TO 10
C
      DO 3 I = 1, N
C
         DO 2 J = 1, N
            ZR(I,J) = 0.0
            ZI(I,J) = 0.0
    2    CONTINUE
C
         ZR(I,I) = 1.0
    3 CONTINUE
C     ********** REDUCE B TO UPPER TRIANGULAR FORM WITH
C                TEMPORARILY REAL DIAGONAL ELEMENTS **********
   10 IF (N .LE. 1) GO TO 170
      NM1 = N - 1
C
      DO 100 L = 1, NM1
         L1 = L + 1
         S = 0.0
C
         DO 20 I = L, N
            S = S + ABS(BR(I,L)) + ABS(BI(I,L))
   20    CONTINUE
C
         IF (S .EQ. ZERO) GO TO 100
         RHO = 0.0
C
         DO 25 I = L, N
            BR(I,L) = BR(I,L) / S
            BI(I,L) = BI(I,L) / S
            RHO = RHO + BR(I,L)**2 + BI(I,L)**2
   25    CONTINUE
C
         R = SQRT(RHO)
         XR = CABS(CMPLX(BR(L,L),BI(L,L)))
         IF (XR .EQ. ZERO) GO TO 27
         RHO = RHO + XR * R
         U1 = -BR(L,L) / XR
         U1I = -BI(L,L) / XR
         YR = R / XR + 1.0
         BR(L,L) = YR * BR(L,L)
         BI(L,L) = YR * BI(L,L)
         GO TO 28
C
   27    BR(L,L) = R
         U1 = -1.0
         U1I = 0.0
C
   28    DO 50 J = L1, N
            T = 0.0
            TI = 0.0
C
            DO 30 I = L, N
               T = T + BR(I,L) * BR(I,J) + BI(I,L) * BI(I,J)
               TI = TI + BR(I,L) * BI(I,J) - BI(I,L) * BR(I,J)
   30       CONTINUE
C
            T = T / RHO
            TI = TI / RHO
C
            DO 40 I = L, N
               BR(I,J) = BR(I,J) - T * BR(I,L) + TI * BI(I,L)
               BI(I,J) = BI(I,J) - T * BI(I,L) - TI * BR(I,L)
   40       CONTINUE
C
            XI = U1 * BI(L,J) - U1I * BR(L,J)
            BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
            BI(L,J) = XI
   50    CONTINUE
C
         DO 80 J = 1, N
            T = 0.0
            TI = 0.0
C
            DO 60 I = L, N
               T = T + BR(I,L) * AR(I,J) + BI(I,L) * AI(I,J)
               TI = TI + BR(I,L) * AI(I,J) - BI(I,L) * AR(I,J)
   60       CONTINUE
C
            T = T / RHO
            TI = TI / RHO
C
            DO 70 I = L, N
               AR(I,J) = AR(I,J) - T * BR(I,L) + TI * BI(I,L)
               AI(I,J) = AI(I,J) - T * BI(I,L) - TI * BR(I,L)
   70       CONTINUE
C
            XI = U1 * AI(L,J) - U1I * AR(L,J)
            AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
            AI(L,J) = XI
   80    CONTINUE
C
         BR(L,L) = R * S
         BI(L,L) = 0.0
C
         DO 90 I = L1, N
            BR(I,L) = 0.0
            BI(I,L) = 0.0
   90    CONTINUE
C
  100 CONTINUE
C     ********** REDUCE A TO UPPER HESSENBERG FORM WITH REAL SUBDIAGONAL
C                ELEMENTS, WHILE KEEPING B TRIANGULAR **********
      DO 160 K = 1, NM1
         K1 = K + 1
C     ********** SET BOTTOM ELEMENT IN K-TH COLUMN OF A REAL **********
         IF (AI(N,K) .EQ. ZERO) GO TO 105
         R = CABS(CMPLX(AR(N,K),AI(N,K)))
         U1 = AR(N,K) / R
         U1I = AI(N,K) / R
         AR(N,K) = R
         AI(N,K) = 0.0
C
         DO 103 J = K1, N
            XI = U1 * AI(N,J) - U1I * AR(N,J)
            AR(N,J) = U1 * AR(N,J) + U1I * AI(N,J)
            AI(N,J) = XI
  103    CONTINUE
C
         XI = U1 * BI(N,N) - U1I * BR(N,N)
         BR(N,N) = U1 * BR(N,N) + U1I * BI(N,N)
         BI(N,N) = XI
  105    IF (K .EQ. NM1) GO TO 170
         NK1 = NM1 - K
C     ********** FOR L=N-1 STEP -1 UNTIL K+1 DO -- **********
         DO 150 LB = 1, NK1
            L = N - LB
            L1 = L + 1
C     ********** ZERO A(L+1,K) **********
            S = ABS(AR(L,K)) + ABS(AI(L,K)) + AR(L1,K)
            IF (S .EQ. ZERO) GO TO 150
            U1 = AR(L,K) / S
            U1I = AI(L,K) / S
            U2 = AR(L1,K) / S
            R = SQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1 / R
            U1I = U1I / R
            U2 = U2 / R
            AR(L,K) = R * S
            AI(L,K) = 0.0
            AR(L1,K) = 0.0
C
            DO 110 J = K1, N
               XR = AR(L,J)
               XI = AI(L,J)
               YR = AR(L1,J)
               YI = AI(L1,J)
               AR(L,J) = U1 * XR + U1I * XI + U2 * YR
               AI(L,J) = U1 * XI - U1I * XR + U2 * YI
               AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
               AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  110       CONTINUE
C
            XR = BR(L,L)
            BR(L,L) = U1 * XR
            BI(L,L) = -U1I * XR
            BR(L1,L) = -U2 * XR
C
            DO 120 J = L1, N
               XR = BR(L,J)
               XI = BI(L,J)
               YR = BR(L1,J)
               YI = BI(L1,J)
               BR(L,J) = U1 * XR + U1I * XI + U2 * YR
               BI(L,J) = U1 * XI - U1I * XR + U2 * YI
               BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
               BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  120       CONTINUE
C     ********** ZERO B(L+1,L) **********
            S = ABS(BR(L1,L1)) + ABS(BI(L1,L1)) + ABS(BR(L1,L))
            IF (S .EQ. ZERO) GO TO 150
            U1 = BR(L1,L1) / S
            U1I = BI(L1,L1) / S
            U2 = BR(L1,L) / S
            R = SQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1 / R
            U1I = U1I / R
            U2 = U2 / R
            BR(L1,L1) = R * S
            BI(L1,L1) = 0.0
            BR(L1,L) = 0.0
C
            DO 130 I = 1, L
               XR = BR(I,L1)
               XI = BI(I,L1)
               YR = BR(I,L)
               YI = BI(I,L)
               BR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               BI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               BR(I,L) = U1 * YR - U1I * YI - U2 * XR
               BI(I,L) = U1 * YI + U1I * YR - U2 * XI
  130       CONTINUE
C
            DO 140 I = 1, N
               XR = AR(I,L1)
               XI = AI(I,L1)
               YR = AR(I,L)
               YI = AI(I,L)
               AR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               AI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               AR(I,L) = U1 * YR - U1I * YI - U2 * XR
               AI(I,L) = U1 * YI + U1I * YR - U2 * XI
  140       CONTINUE
C
            IF (.NOT. MATZ) GO TO 150
C
            DO 145 I = 1, N
               XR = ZR(I,L1)
               XI = ZI(I,L1)
               YR = ZR(I,L)
               YI = ZI(I,L)
               ZR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               ZI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               ZR(I,L) = U1 * YR - U1I * YI - U2 * XR
               ZI(I,L) = U1 * YI + U1I * YR - U2 * XI
  145       CONTINUE
C
  150    CONTINUE
C
  160 CONTINUE
C
  170 RETURN
C     ********** LAST CARD OF CQZHES **********
      END
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE CQZVAL(NM,N,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,           
     X                                       MATZ,ZR,ZI,IERR)
C
      INTEGER I,J,K,L,N,EN,K1,K2,LL,L1,NA,NM,ITS,KM1,LM1,
     X        ENM2,IERR,LOR1,ENORN
      REAL AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),ALFI(N),
     X       BETA(N),ZR(NM,N),ZI(NM,N)
      REAL R,S,A1,A2,EP,SH,U1,U2,XI,XR,YI,YR,ANI,A1I,A33,A34,A43,A44,
     X       BNI,B11,B33,B44,SHI,U1I,A33I,A34I,A43I,A44I,B33I,B44I,
     X       EPSA,EPSB,EPS1,ANORM,BNORM,B3344,B3344I
      REAL SQRT,CABS,ABS
      INTEGER MAX0
      LOGICAL MATZ
      COMPLEX Z3
      COMPLEX CSQRT,CMPLX
      REAL REAL,AIMAG
C
C
C
C
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF STEPS 2 AND 3 OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,
C     AS MODIFIED IN TECHNICAL NOTE NASA TN E-7305(1973) BY WARD.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES, ONE OF THEM
C     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM,
C     THE HESSENBERG MATRIX MUST FURTHER HAVE REAL SUBDIAGONAL ELEMENTS.
C     IT REDUCES THE HESSENBERG MATRIX TO TRIANGULAR FORM USING
C     UNITARY TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM
C     OF THE OTHER MATRIX AND FURTHER MAKING ITS DIAGONAL ELEMENTS
C     REAL AND NON-NEGATIVE.  IT THEN RETURNS QUANTITIES WHOSE RATIOS
C     GIVE THE GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY
C     CQZHES  AND POSSIBLY FOLLOWED BY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER HESSENBERG MATRIX
C          WITH REAL SUBDIAGONAL ELEMENTS,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
C          EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
C          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
C          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
C          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
C          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
C          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
C          BUT LESS ACCURATE RESULTS,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE,
C
C        Z=(ZR,ZI) CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
C          BY  CQZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  IN PARTICULAR, ITS DIAGONAL HAS BEEN SET
C          REAL AND NON-NEGATIVE.  THE LOCATION BR(N,1) IS USED TO
C          STORE EPS1 TIMES THE NORM OF B FOR LATER USE BY  CQZVEC,
C
C        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE
C          DIAGONAL ELEMENTS OF THE TRIANGULARIZED A MATRIX,
C
C        BETA CONTAINS THE REAL NON-NEGATIVE DIAGONAL ELEMENTS OF THE
C          CORRESPONDING B.  THE GENERALIZED EIGENVALUES ARE THEN
C          THE RATIOS ((ALFR+I*ALFI)/BETA),
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE.,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF AR(J,J-1) HAS NOT BECOME
C                     ZERO AFTER 50 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      ZTWO = 2.0
      ZONE = 1.0
      ZERO = 0.0

      IERR = 0
C     ********** COMPUTE EPSA,EPSB **********
      ANORM = 0.0
      BNORM = 0.0
C
      DO 30 I = 1, N
         ANI = 0.0
         IF (I .NE. 1) ANI = ABS(AR(I,I-1))
         BNI = 0.0
C
         DO 20 J = I, N
            ANI = ANI + ABS(AR(I,J)) + ABS(AI(I,J))
            BNI = BNI + ABS(BR(I,J)) + ABS(BI(I,J))
   20    CONTINUE
C
         IF (ANI .GT. ANORM) ANORM = ANI
         IF (BNI .GT. BNORM) BNORM = BNI
   30 CONTINUE
C
      IF (ANORM .EQ. ZERO) ANORM = 1.0
      IF (BNORM .EQ. ZERO) BNORM = 1.0
      EP = EPS1
      IF (EP .GT. ZERO) GO TO 50
C     ********** COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO **********
      EP = 1.0
   40 EP = EP / ZTWO
      IF (ZONE + EP .GT. ZONE) GO TO 40
   50 EPSA = EP * ANORM
      EPSB = EP * BNORM
C     ********** REDUCE A TO TRIANGULAR FORM, WHILE
C                KEEPING B TRIANGULAR **********
      LOR1 = 1
      ENORN = N
      EN = N
C     ********** BEGIN QZ STEP **********
   60 IF (EN .EQ. 0) GO TO 1001
      IF (.NOT. MATZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     ********** CHECK FOR CONVERGENCE OR REDUCIBILITY.
C                FOR L=EN STEP -1 UNTIL 1 DO -- **********
   70 DO 80 LL = 1, EN
         LM1 = EN - LL
         L = LM1 + 1
         IF (L .EQ. 1) GO TO 95
         IF (ABS(AR(L,LM1)) .LE. EPSA) GO TO 90
   80 CONTINUE
C
   90 AR(L,LM1) = 0.0
C     ********** SET DIAGONAL ELEMENT AT TOP OF B REAL **********
   95 B11 = CABS(CMPLX(BR(L,L),BI(L,L)))
      IF (B11     .EQ. ZERO) GO TO 98
      U1 = BR(L,L) / B11
      U1I = BI(L,L) / B11
C
      DO 97 J = L, ENORN
         XI = U1 * AI(L,J) - U1I * AR(L,J)
         AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
         AI(L,J) = XI
         XI = U1 * BI(L,J) - U1I * BR(L,J)
         BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
         BI(L,J) = XI
   97 CONTINUE
C
      BI(L,L) = 0.0
   98 IF (L .NE. EN) GO TO 100
C     ********** 1-BY-1 BLOCK ISOLATED **********
      ALFR(EN) = AR(EN,EN)
      ALFI(EN) = AI(EN,EN)
      BETA(EN) = B11
      EN = NA
      GO TO 60
C     ********** CHECK FOR SMALL TOP OF B **********
  100 L1 = L + 1
      IF (B11 .GT. EPSB) GO TO 120
      BR(L,L) = 0.0
      S = ABS(AR(L,L)) + ABS(AI(L,L)) + ABS(AR(L1,L))
      U1 = AR(L,L) / S
      U1I = AI(L,L) / S
      U2 = AR(L1,L) / S
      R = SQRT(U1*U1+U1I*U1I+U2*U2)
      U1 = U1 / R
      U1I = U1I / R
      U2 = U2 / R
      AR(L,L) = R * S
      AI(L,L) = 0.0
C
      DO 110 J = L1, ENORN
         XR = AR(L,J)
         XI = AI(L,J)
         YR = AR(L1,J)
         YI = AI(L1,J)
         AR(L,J) = U1 * XR + U1I * XI + U2 * YR
         AI(L,J) = U1 * XI - U1I * XR + U2 * YI
         AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
         AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
         XR = BR(L,J)
         XI = BI(L,J)
         YR = BR(L1,J)
         YI = BI(L1,J)
         BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
         BR(L,J) = U1 * XR + U1I * XI + U2 * YR
         BI(L,J) = U1 * XI - U1I * XR + U2 * YI
         BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  110 CONTINUE
C
      LM1 = L
      L = L1
      GO TO 90
C     ********** ITERATION STRATEGY **********
  120 IF (ITS .EQ. 50) GO TO 1000
      IF (ITS .EQ. 10) GO TO 135
C     ********** DETERMINE SHIFT **********
      B33 = BR(NA,NA)
      B33I = BI(NA,NA)
      IF (CABS(CMPLX(B33,B33I)) .GE. EPSB) GO TO 122
      B33 = EPSB
      B33I = 0.0
  122 B44 = BR(EN,EN)
      B44I = BI(EN,EN)
      IF (CABS(CMPLX(B44,B44I)) .GE. EPSB) GO TO 124
      B44 = EPSB
      B44I = 0.0
  124 B3344 = B33 * B44 - B33I * B44I
      B3344I = B33 * B44I + B33I * B44
      A33 = AR(NA,NA) * B44 - AI(NA,NA) * B44I
      A33I = AR(NA,NA) * B44I + AI(NA,NA) * B44
      A34 = AR(NA,EN) * B33 - AI(NA,EN) * B33I
     X    - AR(NA,NA) * BR(NA,EN) + AI(NA,NA) * BI(NA,EN)
      A34I = AR(NA,EN) * B33I + AI(NA,EN) * B33
     X     - AR(NA,NA) * BI(NA,EN) - AI(NA,NA) * BR(NA,EN)
      A43 = AR(EN,NA) * B44
      A43I = AR(EN,NA) * B44I
      A44 = AR(EN,EN) * B33 - AI(EN,EN) * B33I - AR(EN,NA) * BR(NA,EN)
      A44I = AR(EN,EN) * B33I + AI(EN,EN) * B33 - AR(EN,NA) * BI(NA,EN)
      SH = A44
      SHI = A44I
      XR = A34 * A43 - A34I * A43I
      XI = A34 * A43I + A34I * A43
      IF (XR .EQ. ZERO .AND. XI .EQ. ZERO) GO TO 140
      YR = (A33 - SH) / 2.0
      YI = (A33I - SHI) / 2.0
      Z3 = CSQRT(CMPLX(YR**2-YI**2+XR,2.0*YR*YI+XI))
      U1 = REAL(Z3)
      U1I = AIMAG(Z3)
      IF (YR * U1 + YI * U1I .GE. ZERO) GO TO 125
      U1 = -U1
      U1I = -U1I
  125 Z3 = (CMPLX(SH,SHI) - CMPLX(XR,XI) / CMPLX(YR+U1,YI+U1I))
     X   / CMPLX(B3344,B3344I)
      SH = REAL(Z3)
      SHI = AIMAG(Z3)
      GO TO 140
C     ********** AD HOC SHIFT **********
  135 SH = AR(EN,NA) + AR(NA,ENM2)
      SHI = 0.0
C     ********** DETERMINE ZEROTH COLUMN OF A **********
  140 A1 = AR(L,L) / B11 - SH
      A1I = AI(L,L) / B11 - SHI
      A2 = AR(L1,L) / B11
      ITS = ITS + 1
      IF (.NOT. MATZ) LOR1 = L
C     ********** MAIN LOOP **********
      DO 260 K = L, NA
         K1 = K + 1
         K2 = K + 2
         KM1 = MAX0(K-1,L)
C     ********** ZERO A(K+1,K-1) **********
         IF (K .EQ. L) GO TO 170
         A1 = AR(K,KM1)
         A1I = AI(K,KM1)
         A2 = AR(K1,KM1)
  170    S = ABS(A1) + ABS(A1I) + ABS(A2)
         U1 = A1 / S
         U1I = A1I / S
         U2 = A2 / S
         R = SQRT(U1*U1+U1I*U1I+U2*U2)
         U1 = U1 / R
         U1I = U1I / R
         U2 = U2 / R
C
         DO 180 J = KM1, ENORN
            XR = AR(K,J)
            XI = AI(K,J)
            YR = AR(K1,J)
            YI = AI(K1,J)
            AR(K,J) = U1 * XR + U1I * XI + U2 * YR
            AI(K,J) = U1 * XI - U1I * XR + U2 * YI
            AR(K1,J) = U1 * YR - U1I * YI - U2 * XR
            AI(K1,J) = U1 * YI + U1I * YR - U2 * XI
            XR = BR(K,J)
            XI = BI(K,J)
            YR = BR(K1,J)
            YI = BI(K1,J)
            BR(K,J) = U1 * XR + U1I * XI + U2 * YR
            BI(K,J) = U1 * XI - U1I * XR + U2 * YI
            BR(K1,J) = U1 * YR - U1I * YI - U2 * XR
            BI(K1,J) = U1 * YI + U1I * YR - U2 * XI
  180    CONTINUE
C
         IF (K .EQ. L) GO TO 240
         AI(K,KM1) = 0.0
         AR(K1,KM1) = 0.0
         AI(K1,KM1) = 0.0
C     ********** ZERO B(K+1,K) **********
  240    S = ABS(BR(K1,K1)) + ABS(BI(K1,K1)) + ABS(BR(K1,K))
         U1 = BR(K1,K1) / S
         U1I = BI(K1,K1) / S
         U2 = BR(K1,K) / S
         R = SQRT(U1*U1+U1I*U1I+U2*U2)
         U1 = U1 / R
         U1I = U1I / R
         U2 = U2 / R
         IF (K .EQ. NA) GO TO 245
         XR = AR(K2,K1)
         AR(K2,K1) = U1 * XR
         AI(K2,K1) = -U1I * XR
         AR(K2,K) = -U2 * XR
C
  245    DO 250 I = LOR1, K1
            XR = AR(I,K1)
            XI = AI(I,K1)
            YR = AR(I,K)
            YI = AI(I,K)
            AR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            AI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            AR(I,K) = U1 * YR - U1I * YI - U2 * XR
            AI(I,K) = U1 * YI + U1I * YR - U2 * XI
            XR = BR(I,K1)
            XI = BI(I,K1)
            YR = BR(I,K)
            YI = BI(I,K)
            BR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            BI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            BR(I,K) = U1 * YR - U1I * YI - U2 * XR
            BI(I,K) = U1 * YI + U1I * YR - U2 * XI
  250    CONTINUE
C
         BI(K1,K1) = 0.0
         BR(K1,K) = 0.0
         BI(K1,K) = 0.0
         IF (.NOT. MATZ) GO TO 260
C
         DO 255 I = 1, N
            XR = ZR(I,K1)
            XI = ZI(I,K1)
            YR = ZR(I,K)
            YI = ZI(I,K)
            ZR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            ZI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            ZR(I,K) = U1 * YR - U1I * YI - U2 * XR
            ZI(I,K) = U1 * YI + U1I * YR - U2 * XI
  255    CONTINUE
C
  260 CONTINUE
C     ********** SET LAST A SUBDIAGONAL REAL AND END QZ STEP **********
      IF (AI(EN,NA) .EQ. ZERO) GO TO 70
      R = CABS(CMPLX(AR(EN,NA),AI(EN,NA)))
      U1 = AR(EN,NA) / R
      U1I = AI(EN,NA) / R
      AR(EN,NA) = R
      AI(EN,NA) = 0.0
C
      DO 270 J = EN, ENORN
         XI = U1 * AI(EN,J) - U1I * AR(EN,J)
         AR(EN,J) = U1 * AR(EN,J) + U1I * AI(EN,J)
         AI(EN,J) = XI
         XI = U1 * BI(EN,J) - U1I * BR(EN,J)
         BR(EN,J) = U1 * BR(EN,J) + U1I * BI(EN,J)
         BI(EN,J) = XI
  270 CONTINUE
C
      GO TO 70
C     ********** SET ERROR -- BOTTOM SUBDIAGONAL ELEMENT HAS NOT
C                BECOME NEGLIGIBLE AFTER 50 ITERATIONS **********
 1000 IERR = EN
C     ********** SAVE EPSB FOR USE BY CQZVEC **********
 1001 IF (N .GT. 1) BR(N,1) = EPSB
      RETURN
C     ********** LAST CARD OF CQZVAL **********
      END
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE CQZVEC(NM,N,AR,AI,BR,BI,ALFR,ALFI,BETA,ZR,ZI)          
C
      INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN
      REAL AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),ALFI(N),
     X       BETA(N),ZR(NM,N),ZI(NM,N)
      REAL R,T,RI,TI,XI,ALMI,ALMR,BETM,EPSB
      REAL CABS
      COMPLEX Z3
      COMPLEX CMPLX
      REAL REAL,AIMAG
C
C
C
C
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FOURTH STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES IN UPPER
C     TRIANGULAR FORM, WHERE ONE OF THEM FURTHER MUST HAVE REAL DIAGONAL
C     ELEMENTS.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM
C     AND TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
C     IT IS USUALLY PRECEDED BY  CQZHES  AND  CQZVAL.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX WITH REAL
C          DIAGONAL ELEMENTS.  IN ADDITION, LOCATION BR(N,1) CONTAINS
C          THE TOLERANCE QUANTITY (EPSB) COMPUTED AND SAVED IN  CQZVAL,
C
C        ALFR, ALFI, AND BETA ARE VECTORS WITH COMPONENTS WHOSE
C          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
C          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  CQZVAL,
C
C        Z=(ZR,ZI) CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTIONS BY  CQZHES  AND  CQZVAL, IF PERFORMED.
C          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
C          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
C
C     ON OUTPUT-
C
C        A IS UNALTERED,
C
C        B HAS BEEN DESTROYED,
C
C        ALFR, ALFI, AND BETA ARE UNALTERED,
C
C        Z CONTAINS THE EIGENVECTORS.  EACH EIGENVECTOR IS NORMALIZED
C          SO THAT THE MODULUS OF ITS LARGEST COMPONENT IS 1.0 .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      ZERO = 0.0
      IF (N .LE. 1) GO TO 1001
      EPSB = BR(N,1)
C     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********
      DO 800 NN = 2, N
         EN = N + 2 - NN
         NA = EN - 1
         ALMR = ALFR(EN)
         ALMI = ALFI(EN)
         BETM = BETA(EN)
C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
         DO 700 II = 1, NA
            I = EN - II
            R = 0.0
            RI = 0.0
            M = I + 1
C
            DO 610 J = M, EN
               T = BETM * AR(I,J) - ALMR * BR(I,J) + ALMI * BI(I,J)
               TI = BETM * AI(I,J) - ALMR * BI(I,J) - ALMI * BR(I,J)
               IF (J .EQ. EN) GO TO 605
               XI = T * BI(J,EN) + TI * BR(J,EN)
               T = T * BR(J,EN) - TI * BI(J,EN)
               TI = XI
  605          R = R + T
               RI = RI + TI
  610       CONTINUE
C
            T = ALMR * BETA(I) - BETM * ALFR(I)
            TI = ALMI * BETA(I) - BETM * ALFI(I)
            IF (T .EQ. ZERO .AND. TI .EQ. ZERO) T = EPSB
            Z3 = CMPLX(R,RI) / CMPLX(T,TI)
            BR(I,EN) = REAL(Z3)
            BI(I,EN) = AIMAG(Z3)
  700    CONTINUE
C
  800 CONTINUE
C     ********** END BACK SUBSTITUTION.
C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
C                FOR J=N STEP -1 UNTIL 2 DO -- **********
      DO 880 JJ = 2, N
         J = N + 2 - JJ
         M = J - 1
C
         DO 880 I = 1, N
C
            DO 860 K = 1, M
               ZR(I,J) = ZR(I,J) + ZR(I,K) * BR(K,J) - ZI(I,K) * BI(K,J)
               ZI(I,J) = ZI(I,J) + ZR(I,K) * BI(K,J) + ZI(I,K) * BR(K,J)
  860       CONTINUE
C
  880 CONTINUE
C     ********** NORMALIZE SO THAT MODULUS OF LARGEST
C                COMPONENT OF EACH VECTOR IS 1 **********
      DO 950 J = 1, N
         T = 0.0
C
         DO 930 I = 1, N
            R = CABS(CMPLX(ZR(I,J),ZI(I,J)))
            IF (R .GT. T) T = R
  930    CONTINUE
C
         DO 940 I = 1, N
            ZR(I,J) = ZR(I,J) / T
            ZI(I,J) = ZI(I,J) / T
  940    CONTINUE
C
  950 CONTINUE
C
 1001 RETURN
C     ********** LAST CARD OF CQZVEC **********
      END



















