C
C     ********************************************
C     **  Time splitting through interpolation  **
C     ********************************************
C
      SUBROUTINE TIMESPL(T,F,TA,FA,NTMAX,NTM,IERR)
C
C     You can choose interpolation method;
C        either spline, aitken-neville or laglange method.
C     From the aspect of trade-off between accuracy and expended time,
C        we recommend you use Aitken-Neville method.
C
C     <input>
C        T     : designated time
C        TA    : time array
C        FA    : fuction value array
C        NTMAX : maximum number
C        NTM   : maximum number of array
C     <output>
C        F     : function value for designated time
C        IERR  : error indicator
C
      IMPLICIT NONE
      INTEGER NTMAX,NTM,IERR
      REAL*8 T,F,TA(NTM),FA(NTM)
      INTEGER NT,ID,MODE,SPL,AIT,LAG
      REAL*8 DERIV4P,FXA(NTM),U(4,NTM)
C
      SPL=0
      AIT=1
      LAG=2
      MODE=AIT
C
      IF(MODE.EQ.0) THEN
C     calculate first derivative
         DO NT=1,NTMAX-3
            FXA(NT)=DERIV4P(FA(NT),FA(NT+1),FA(NT+2),FA(NT+3),
     &                      TA(NT),TA(NT+1),TA(NT+2),TA(NT+3))
         ENDDO
         DO NT=NTMAX-2,NTMAX
            FXA(NT)=DERIV4P(FA(NT),FA(NT-1),FA(NT-2),FA(NT-3),
     &                      TA(NT),TA(NT-1),TA(NT-2),TA(NT-3))
         ENDDO
         CALL SPL1D(TA,FA,FXA,U,NTMAX,ID,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TIMESPL: SPL1D ERROR',IERR
C
         CALL SPL1DF(T,F,TA,U,NTMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TIMESPL: SPL1DF ERROR',IERR
      ELSEIF(MODE.EQ.1) THEN
         CALL AITKEN(T,F,TA,FA,4,NTMAX)
      ELSEIF(MODE.EQ.2) THEN
         CALL LAGLANGE(T,F,TA,FA,NTMAX,NTM,IERR)
C         IF(IERR.NE.0) WRITE(6,*) 'XX TIMESPL: LAGLANGE ERROR',IERR
      ENDIF
C
      RETURN
      END
C
C   *******************************************
C   **    LAGRANGEAN INTERPOLATION METHOD    **
C   *******************************************
C
C     input:
C
C     SLT     : Interesting Time
C     T(NTM)  : Total Time Data
C     F1(NTM) : Functional Values
C     NTMAX   : Maximum Number of the Time Mesh from UFILE
C     NTM     : Maximum Array Width for Time
C
C     output:
C
C     FOUT    : Interpolated Value
C     IERR    : Error Indicator
C     
C     ***********************************************************
C
      SUBROUTINE LAGLANGE(SLT,FOUT,T,F1,NTMAX,NTM,IERR)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
      COMMON /COMEPS/ EPS,IERRL
      DIMENSION T(NTM),F1(NTM)
C
      M=5
C      M=1
      EPS=1.D-5
C
      FOUT=FITLAG(SLT,T,F1,M,NTMAX,NTM)
      IERR=IERRL
C
      RETURN
      END
C
      FUNCTION FITLAG(X,A,F,M,N,NTM)
****************************************************
*                                                  *
*     Compute iterated Lagrange interpolation      *
*     at X based on data given in A(I) and F(I)    *
*                                                  *
*     ==== input data for arguments  ====          *
*                                                  *
*     N...total number of data (A(I),F(I)) given   *
*         in the table                             *
*     A(I)...I-th sampling point                   *
*            A(I) must be given in ascending       *
*            order with respect to I               *
*     F(I)...data at I-th point, i.e. FUNC(A(I))   *
*     M...number of sampling points used in both   *
*         sides of X for interpolation             *
*         M must be less than or equal to 5        *
*                                                  *
*     ==== input data for common block ====        *
*                                                  *
*     EPS...absolute error tolerance               *
*                                                  *
*     ==== output data for common block ====       *
*                                                  *
*     IERR...if IERR = 0 then converged            *
*            if IERR = 1 then not converged        *
*                                                  *
*     ==== work arrays ====                        *
*                                                  *
*     B(J)...A(I) - X                              *
*     V(I,J)...triangular table for interpolation  *
*                                                  *
****************************************************
      IMPLICIT REAL*8 (A-F,H,O-Z)
      DIMENSION V(10,10),B(10)
C
      COMMON /COMEPS/ EPS,IERR
      DIMENSION A(NTM),F(NTM)
C
      IF(M.GT.5) STOP 'M MUST BE LESS THAN 6.'
C
      M2=2*M
C
C     ---- find the nearest sampling point to x ----
C                 by bisection method
C
      IS=1
      IB=N
C
 1    CONTINUE
C
      IM=(IS+IB)/2
      IF(X.LT.A(IM)) THEN
         IB=IM
      ELSEIF(X.GT.A(IM)) THEN
         IS=IM
      ELSE
         FITLAG=F(IM)
         IERR=0
         RETURN
      ENDIF
C
      IF(IS.LT.IB-1) GOTO 1
C
C     ---- set sampling points to use ----
C
      IL=IS-M+1
      IR=IB+M-1
      IF(IL.LT.1) THEN
         IL=1
      ELSEIF(IR.GT.N) THEN
         IL=N-M2+1
      ENDIF
C
      DO I=1,M2
         B(I)=A(I+IL-1)-X
         V(1,I)=F(I+IL-1)
      ENDDO
C
C     ---- sort ABS(B(I)) in ascending order ----
C
      DO I=2,M2
         K=I-1
         DO J=1,M2
            IF(ABS(B(J)).LT.ABS(B(K))) K=J
         ENDDO
         IF(K.NE.I-1) THEN
            BV=B(K)
            B(K)=B(I-1)
            B(I-1)=BV
            VV=V(1,K)
            V(1,K)=V(1,I-1)
            V(1,I-1)=VV
         ENDIF
      ENDDO
C
C     ---- compute iterated interpolation ----
C
      DO I=2,M2
C
         DO J=2,I
C            write(6,*) I,J-1,B(I),B(J-1)
            V(J,I)=(V(J-1,J-1)*B(I)
     &            - V(J-1,I  )*B(J-1))/(B(I)-B(J-1))
         ENDDO
C
         IF(ABS(V(I,I)-V(I-1,I-1)).LE.EPS) THEN
C
C     ---- converged ----
C
            FITLAG=V(I,I)
            IERR=0
            RETURN
         ENDIF
C
      ENDDO
C
C     ---- not converged ----
C
      FITLAG=V(M2,M2)
      IERR=1
      RETURN
C
      END
C
C     *****************************************************
C     **  Aitken-Neville iterated linear interpolations  **
C     *****************************************************
C
      FUNCTION AITKEN2P(X0,F1,F2,F3,X1,X2,X3)
C
      IMPLICIT NONE
      INTEGER I,J,K,M,KPM,NP1MM
      REAL*8 AITKEN2P,X0,F1,F2,F3,X1,X2,X3
      REAL*8 X(3), Y(3), VAL(3,3)
C
      M=2
      Y(1)=F3
      Y(2)=F2
      Y(3)=F1
      X(1)=X3
      X(2)=X2
      X(3)=X1
C
      DO J = 1, M
         VAL(1,J) = (Y(J) * (X0 - X(J+1)) - Y(J+1) * (X0 - X(J)))
     &             /(X(J) - X(J+1))
      ENDDO
      DO I = 2, M
         NP1MM = M + 1 - I
         DO K = 1, NP1MM
            KPM = K + I
            VAL(I,K) = (VAL(I-1,K) * (X0 - X(KPM)) - VAL(I-1,K+1)
     &                * (X0 - X(K))) / (X(K) - X(KPM))
         ENDDO
      ENDDO
      AITKEN2P = VAL(M,1)
C
      RETURN
      END
C
C     ***
C
      FUNCTION AITKEN4P(X0,F1,F2,F3,F4,F5,X1,X2,X3,X4,X5)
C
      IMPLICIT NONE
      INTEGER I,J,K,M,KPM,NP1MM
      REAL*8 AITKEN4P,X0,F1,F2,F3,F4,F5,X1,X2,X3,X4,X5
      REAL*8 X(5), Y(5), VAL(5,5)
C
      M=4
      Y(1)=F5
      Y(2)=F4
      Y(3)=F3
      Y(4)=F2
      Y(5)=F1
      X(1)=X5
      X(2)=X4
      X(3)=X3
      X(4)=X2
      X(5)=X1
C
      DO J = 1, M
         VAL(1,J) = (Y(J) * (X0 - X(J+1)) - Y(J+1) * (X0 - X(J)))
     &             /(X(J) - X(J+1))
      ENDDO
      DO I = 2, M
         NP1MM = M + 1 - I
         DO K = 1, NP1MM
            KPM = K + I
            VAL(I,K) = (VAL(I-1,K) * (X0 - X(KPM)) - VAL(I-1,K+1)
     &                * (X0 - X(K))) / (X(K) - X(KPM))
         ENDDO
      ENDDO
      AITKEN4P = VAL(M,1)
C
      RETURN
      END
C
C     ***
C
      SUBROUTINE AITKEN(X0,F0,XI,F,M,N)
C
C     <input>
C        X0    : designated point
C        XI(N) : data point array
C        F(N)  : data array
C        M     : order (up to 10)
C        N     : maximum array argument
C     <output>
C        F0    : interpolated value
C
      IMPLICIT NONE
      INTEGER I, I1, J, K, KPM, L, M, M1, M2, N, NP1MM
      REAL*8 X0, X1, F0, XI(N), F(N)
      REAL*8 X(10), Y(10), VAL(10,10), XMIN, XMAX
C
      X1 = XI(1)
      XMIN = X1
      XMAX = X1
      M2 = (M + 1)/2
      DO I = 1, N
         IF(XI(I).LE.XMIN) XMIN=XI(I)
         IF(XI(I).GE.XMAX) XMAX=XI(I)
         IF(XI(I).GE.X0) THEN
            K=I-1
            GOTO 10
         ENDIF
      ENDDO
      IF(X0.LT.XMIN) THEN
         K = 0
      ELSEIF(X0.GT.XMAX) THEN
         K = N
      ENDIF
 10   IF (K .GT. M2) THEN
         IF (N .GT. K + M2) THEN
            L = M2
         ELSE
            L = M + 1 - N + K
         ENDIF
      ELSE
         L = K
      ENDIF
      M1 = M + 1
      DO I = 1, M1
         I1 = I + K - L
         Y(I) = F(I1)
         X(I) = XI(I1)
      ENDDO
C
      DO J = 1, M
         VAL(1,J) = (Y(J) * (X0 - X(J+1)) - Y(J+1) * (X0 - X(J)))
     &             /(X(J) - X(J+1))
      ENDDO
      DO I = 2, M
         NP1MM = M + 1 - I
         DO K = 1, NP1MM
            KPM = K + I
            VAL(I,K) = (VAL(I-1,K) * (X0 - X(KPM)) - VAL(I-1,K+1)
     &                * (X0 - X(K))) / (X(K) - X(KPM))
         ENDDO
      ENDDO
      F0 = VAL(M,1)
C
      RETURN
      END
C
C     ****************************************************
C     **  Derivative calculated from adjacent 3 points  **
C     ****************************************************
C     -- This formulation has a third-order accuracy. --
C
      FUNCTION DERIV4P(F0,F1,F2,F3,X0,X1,X2,X3)
C
      IMPLICIT NONE
      REAL*8 DERIV4P,F0,F1,F2,F3,X0,X1,X2,X3
      REAL*8 DX1,DX2,DX3
C
      DX1 = X1 - X0
      DX2 = X2 - X0
      DX3 = X3 - X0
C
      DERIV4P = - F1 * DX2**2 * DX3**2 * (DX2 - DX3)
     &          - F2 * DX3**2 * DX1**2 * (DX3 - DX1)
     &          - F3 * DX1**2 * DX2**2 * (DX1 - DX2)
     &          - F0 * (DX1 - DX2) * (DX2 - DX3) * (DX3 - DX1)
     &              * (DX1 * DX2 + DX2 * DX3 + DX3 * DX1)
      DERIV4P = DERIV4P / (DX1 * DX2 * DX3
     &                  * (DX1 - DX2) * (DX2 - DX3) * (DX3 - DX1))
C
      RETURN
      END
C
C     ****************************************************
C     **  Derivative calculated from adjacent 3 points  **
C     ****************************************************
C     -- This formulation has a second-order accuracy. --
C
      FUNCTION DERIV3P(f0, f1, f2, x10, x11, x12)
C
      IMPLICIT NONE
      REAL*8 DERIV3P, f0, f1, f2, x10, x11, x12
      REAL*8 dx11, dx12
C
      dx11 = x11 - x10
      dx12 = x12 - x10
      DERIV3P = (dx12**2 * f1 - dx11**2 * f2)
     &            / (dx11 * dx12 * (dx12 - dx11))
     &            - (dx12 + dx11) / (dx11 * dx12) * f0
C
      RETURN
      END
C
C     ****************************************************
C     **  Derivative calculated from adjacent 2 points  **
C     ****************************************************
C     -- This formulation has a second-order accuracy. --
C
      FUNCTION DERIV3(NR,R,F,NRMAX,NRM)
C
      IMPLICIT NONE
      INTEGER NR,NRMAX,NRM
      REAL*8 DERIV3,R(NRM),F(NRM)
      REAL*8 DLT,DLT1,DLT2
C
      IF(NR.EQ.1) THEN
         DLT=R(NR+1)-R(NR)
         DERIV3 = (-3.D0*F(NR)+4.D0*F(NR+1)-F(NR+2))/(2.D0*DLT)
      ELSEIF(NR.EQ.NRMAX) THEN
         DLT=R(NR-1)-R(NR)
         DERIV3 = (-3.D0*F(NR)+4.D0*F(NR-1)-F(NR-2))/(2.D0*DLT)
      ELSE
         DLT1=R(NR-1)-R(NR)
         DLT2=R(NR+1)-R(NR)
         DERIV3 = (DLT2**2*F(NR-1)-DLT1**2*F(NR+1))
     &            /(DLT1*DLT2*(DLT2-DLT1))
     &            -(DLT2+DLT1)*F(NR)/(DLT1*DLT2)
      ENDIF
C
      RETURN
      END
C
C     *** Extrapolate center value ***************************
C     *  assuming the gradient is zero at the center (rho=0) *
C     ********************************************************
C
      FUNCTION FCTR(R1,R2,F1,F2)
C
      IMPLICIT NONE
      REAL*8 FCTR,R1,R2,F1,F2
C
      FCTR = (R2**2*F1-R1**2*F2)/(R2**2-R1**2)
C
      RETURN
      END
