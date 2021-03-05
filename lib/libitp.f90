!     ********************************************
!     **  Time splitting through interpolation  **
!     ********************************************

      SUBROUTINE TIMESPL(T,F,TA,FA,NTMAX,NTM,IERR)

!     You can choose interpolation method;
!        either spline, aitken-neville or laglange method.
!     From the aspect of trade-off between accuracy and expended time,
!        we recommend you use Aitken-Neville method.

!     <input>
!        T     : designated time
!        TA    : time array
!        FA    : fuction value array
!        NTMAX : maximum number
!        NTM   : maximum number of array
!     <output>
!        F     : function value for designated time
!        IERR  : error indicator

      USE task_kinds,ONLY: dp
      USE libspl1d  
      implicit none
      integer,              intent(in)  :: NTMAX, NTM
      real(dp),                 intent(in)  :: T
      real(dp), dimension(NTM), intent(in)  :: TA, FA
      integer,              intent(out) :: IERR
      real(dp),                 intent(out) :: F
      integer :: NT, ID, MODE, SPL, AIT, LAG
      real(dp)    :: DERIV4P
      real(dp), dimension(NTM)   :: FXA
      real(dp), dimension(4,NTM) :: U

      SPL=0
      AIT=1
      LAG=2
      MODE=AIT

      select case(MODE)
      case(0)
!     calculate first derivative
         DO NT = 1, NTMAX-3
            FXA(NT)=DERIV4P(FA(NT),FA(NT+1),FA(NT+2),FA(NT+3), &
     &                      TA(NT),TA(NT+1),TA(NT+2),TA(NT+3))
         ENDDO
         DO NT = NTMAX-2, NTMAX
            FXA(NT)=DERIV4P(FA(NT),FA(NT-1),FA(NT-2),FA(NT-3), &
     &                      TA(NT),TA(NT-1),TA(NT-2),TA(NT-3))
         ENDDO
         CALL SPL1D(TA,FA,FXA,U,NTMAX,ID,IERR)
         IF(IERR /= 0) WRITE(6,*) 'XX TIMESPL: SPL1D ERROR',IERR

         CALL SPL1DF(T,F,TA,U,NTMAX,IERR)
         IF(IERR /= 0) WRITE(6,*) 'XX TIMESPL: SPL1DF ERROR',IERR
      case(1)
         CALL AITKEN(T,F,TA,FA,4,NTMAX)
      case(2)
         CALL LAGLANGE(T,F,TA,FA,NTMAX,NTM,IERR)
!         IF(IERR /= 0) WRITE(6,*) 'XX TIMESPL: LAGLANGE ERROR',IERR
      end select

      RETURN
      END SUBROUTINE TIMESPL

!   *******************************************
!   **    LAGRANGEAN INTERPOLATION METHOD    **
!   *******************************************

!     input:
!
!     SLT     : Interesting Time
!     T(NTM)  : Total Time Data
!     F1(NTM) : Functional Values
!     NTMAX   : Maximum Number of the Time Mesh from UFILE
!     NTM     : Maximum Array Width for Time
!
!     output:
!
!     FOUT    : Interpolated Value
!     IERR    : Error Indicator
!     
!     ***********************************************************

      SUBROUTINE LAGLANGE(SLT,FOUT,T,F1,NTMAX,NTM,IERR)

      USE task_kinds,ONLY: dp
      implicit none
      integer,              intent(in)  :: NTMAX, NTM
      real(dp),                 intent(in)  :: SLT
      real(dp), dimension(NTM), intent(in)  :: T, F1
      integer,              intent(out) :: IERR
      real(dp),                 intent(out) :: FOUT
      integer :: M, IERRL
      real(dp) :: FITLAG, EPS
!      COMMON /COMEPS/ EPS,IERRL

      M = 5
!      M=1
      EPS = 1.D-5

      FOUT = FITLAG(SLT,T,F1,M,NTMAX,NTM,EPS,IERRL)
      IERR = IERRL

      RETURN
      END SUBROUTINE LAGLANGE

      FUNCTION FITLAG(X,A,F,M,N,NTM,EPS,IERR)
!****************************************************
!*                                                  *
!*     Compute iterated Lagrange interpolation      *
!*     at X based on data given in A(I) and F(I)    *
!*                                                  *
!*     ==== input data for arguments  ====          *
!*                                                  *
!*     N...total number of data (A(I),F(I)) given   *
!*         in the table                             *
!*     A(I)...I-th sampling point                   *
!*            A(I) must be given in ascending       *
!*            order with respect to I               *
!*     F(I)...data at I-th point, i.e. FUNC(A(I))   *
!*     M...number of sampling points used in both   *
!*         sides of X for interpolation             *
!*         M must be less than or equal to 5        *
!*                                                  *
!*     ==== input data for common block ====        *
!*                                                  *
!*     EPS...absolute error tolerance               *
!*                                                  *
!*     ==== output data for common block ====       *
!*                                                  *
!*     IERR...if IERR = 0 then converged            *
!*            if IERR = 1 then not converged        *
!*                                                  *
!*     ==== work arrays ====                        *
!*                                                  *
!*     B(J)...A(I) - X                              *
!*     V(I,J)...triangular table for interpolation  *
!*                                                  *
!****************************************************
      USE task_kinds,ONLY: dp
      implicit none
      integer,              intent(in) :: M, N, NTM
      real(dp),                 intent(in) :: X, EPS
      real(dp), dimension(NTM), intent(in) :: A, F
      integer, intent(out) :: IERR
      real(dp):: FITLAG
      integer :: I, J, K, M2, IS, IB, IM, IL, IR
      real(dp) :: BV, VV
      real(dp), dimension(10,10) :: V
      real(dp), dimension(10)    :: B

!      COMMON /COMEPS/ EPS,IERR

      IF(M > 5) STOP 'M MUST BE LESS THAN 6.'

      M2 = 2*M

!     ---- find the nearest sampling point to x ----
!                 by bisection method

      IS=1
      IB=N

 1    CONTINUE

      IM=(IS+IB)/2
      IF(X < A(IM)) THEN
         IB=IM
      ELSE IF(X > A(IM)) THEN
         IS=IM
      ELSE
         FITLAG=F(IM)
         IERR=0
         RETURN
      ENDIF

      IF(IS.LT.IB-1) GOTO 1

!     ---- set sampling points to use ----

      IL=IS-M+1
      IR=IB+M-1
      IF(IL < 1) THEN
         IL=1
      ELSE IF(IR > N) THEN
         IL=N-M2+1
      ENDIF

      DO I=1,M2
         B(I)  =A(I+IL-1)-X
         V(1,I)=F(I+IL-1)
      ENDDO

!     ---- sort ABS(B(I)) in ascending order ----

      DO I=2,M2
         K=I-1
         DO J=1,M2
            IF(ABS(B(J)).LT.ABS(B(K))) K=J
         ENDDO
         IF(K /= I-1) THEN
            BV=B(K)
            B(K)=B(I-1)
            B(I-1)=BV
            VV=V(1,K)
            V(1,K)=V(1,I-1)
            V(1,I-1)=VV
         ENDIF
      ENDDO

!     ---- compute iterated interpolation ----

      DO I=2,M2

         DO J=2,I
!            write(6,*) I,J-1,B(I),B(J-1)
            V(J,I)=(V(J-1,J-1)*B(I) - V(J-1,I  )*B(J-1))/(B(I)-B(J-1))
         ENDDO

         IF(ABS(V(I,I)-V(I-1,I-1)).LE.EPS) THEN

!     ---- converged ----

            FITLAG=V(I,I)
            IERR=0
            RETURN
         ENDIF

      ENDDO

!     ---- not converged ----

      FITLAG=V(M2,M2)
      IERR=1
      RETURN
!
      END FUNCTION FITLAG

!     *****************************************************
!     **  Aitken-Neville iterated linear interpolations  **
!     *****************************************************

      FUNCTION AITKEN2P(X0,F1,F2,F3,X1,X2,X3)

      USE task_kinds,ONLY: dp
      implicit none
      real(dp), intent(in) :: X0, F1, F2, F3, X1, X2, X3
      real(dp):: AITKEN2P
      integer :: I, J, K, M, KPM, NP1MM
      real(dp), dimension(3) :: X, Y
      real(dp), dimension(3,3) :: VAL

      M=2
      Y(1)=F3
      Y(2)=F2
      Y(3)=F1
      X(1)=X3
      X(2)=X2
      X(3)=X1

      DO J = 1, M
         VAL(1,J) = (Y(J) * (X0 - X(J+1)) - Y(J+1) * (X0 - X(J))) &
     &             /(X(J) - X(J+1))
      ENDDO
      DO I = 2, M
         NP1MM = M + 1 - I
         DO K = 1, NP1MM
            KPM = K + I
            VAL(I,K) = (VAL(I-1,K) * (X0 - X(KPM)) - VAL(I-1,K+1) &
     &                * (X0 - X(K))) / (X(K) - X(KPM))
         ENDDO
      ENDDO
      AITKEN2P = VAL(M,1)

      RETURN
      END FUNCTION AITKEN2P

!     ***

      FUNCTION AITKEN4P(X0,F1,F2,F3,F4,F5,X1,X2,X3,X4,X5)

      USE task_kinds,ONLY: dp
      implicit none
      real(dp), intent(in) :: X0, F1, F2, F3, F4, F5, X1, X2, X3, X4, X5
      real(dp):: AITKEN4P
      integer :: I, J, K, M, KPM, NP1MM
      real(dp), dimension(5) :: X, Y
      real(dp), dimension(5,5) :: VAL

      M=4
      Y(1)=F5
      Y(2)=F4
      Y(3)=F3
      Y=F2
      Y(5)=F1
      X(1)=X5
      X(2)=X4
      X(3)=X3
      X=X2
      X(5)=X1

      DO J = 1, M
         VAL(1,J) = (Y(J) * (X0 - X(J+1)) - Y(J+1) * (X0 - X(J))) &
     &             /(X(J) - X(J+1))
      ENDDO
      DO I = 2, M
         NP1MM = M + 1 - I
         DO K = 1, NP1MM
            KPM = K + I
            VAL(I,K) = (VAL(I-1,K) * (X0 - X(KPM)) - VAL(I-1,K+1) &
     &                * (X0 - X(K))) / (X(K) - X(KPM))
         ENDDO
      ENDDO
      AITKEN4P = VAL(M,1)

      RETURN
      END FUNCTION AITKEN4P

!     ***

      SUBROUTINE AITKEN(X0,F0,XI,F,M,N)

!     <input>
!        X0    : designated point
!        XI(N) : data point array
!        F(N)  : data array
!        M     : order (up to 10)
!        N     : maximum array argument
!     <output>
!        F0    : interpolated value

      USE task_kinds,ONLY: dp
      implicit none
      integer,            intent(in) :: M, N
      real(dp),               intent(in) :: X0
      real(dp), dimension(N), intent(in) :: XI, F
      real(dp), intent(out) :: F0
      integer :: I, I1, J, K, KPM, L, M1, M2, NP1MM
      real(dp) :: X1, XMIN, XMAX
      real(dp), dimension(10) :: X, Y
      real(dp), dimension(10,10) :: VAL

      X1 = XI(1)
      XMIN = X1
      XMAX = X1
      M2 = (M + 1)/2
      DO I = 1, N
         IF(XI(I) <= XMIN) XMIN=XI(I)
         IF(XI(I) >= XMAX) XMAX=XI(I)
         IF(XI(I) >= X0) THEN
            K=I-1
            GOTO 10
         ENDIF
      ENDDO
      IF(X0 < XMIN) THEN
         K = 0
      ELSEIF(X0 > XMAX) THEN
         K = N
      ENDIF
 10   IF (K > M2) THEN
         IF (N > K + M2) THEN
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

      DO J = 1, M
         VAL(1,J) = (Y(J) * (X0 - X(J+1)) - Y(J+1) * (X0 - X(J))) &
     &             /(X(J) - X(J+1))
      ENDDO
      DO I = 2, M
         NP1MM = M + 1 - I
         DO K = 1, NP1MM
            KPM = K + I
            VAL(I,K) = (VAL(I-1,K) * (X0 - X(KPM)) - VAL(I-1,K+1) &
     &                * (X0 - X(K))) / (X(K) - X(KPM))
         ENDDO
      ENDDO
      F0 = VAL(M,1)

      RETURN
      END SUBROUTINE AITKEN

!     *******************************************
!     **  Derivative calculated from 4 points  **
!     *******************************************
!     -- This formulation has a third-order accuracy. --

      FUNCTION DERIV4P(F0,F1,F2,F3,X0,X1,X2,X3)

      USE task_kinds,ONLY: dp
      implicit none
      real(dp), intent(in) :: F0, F1, F2, F3, X0, X1, X2, X3
      real(dp):: DERIV4P
      real(dp) :: DX1, DX2, DX3

      DX1 = X1 - X0
      DX2 = X2 - X0
      DX3 = X3 - X0

      DERIV4P = - F1 * DX2**2 * DX3**2 * (DX2 - DX3) &
     &          - F2 * DX3**2 * DX1**2 * (DX3 - DX1) &
     &          - F3 * DX1**2 * DX2**2 * (DX1 - DX2) & 
     &          - F0 * (DX1 - DX2) * (DX2 - DX3) * (DX3 - DX1) &
     &              * (DX1 * DX2 + DX2 * DX3 + DX3 * DX1)
      DERIV4P = DERIV4P / (DX1 * DX2 * DX3 &
     &                  * (DX1 - DX2) * (DX2 - DX3) * (DX3 - DX1))

      RETURN
      END FUNCTION DERIV4P

!     ** Array input version **

      FUNCTION DERIV4(NR,R,F,NRMAX,ID)

!     ID = 0    : NR = 0 to NRMAX  ==> NR = 1 to NRMAX+1
!          else : NR = 1 to NRMAX

      USE task_kinds,ONLY: dp
      implicit none
      integer,                    intent(in) :: NR, NRMAX,ID
!      real(dp), dimension(1:NRMAX+1), intent(in) :: R, F
      real(dp), dimension(*), intent(in) :: R, F
      real(dp):: DERIV4
      integer :: NRL, NRLMAX, NR0, NR1, NR2, NR3
      real(dp) :: DX1, DX2, DX3

      IF(ID == 0) THEN
         NRL=NR+1
         NRLMAX=NRMAX+1
      ELSE
         NRL=NR
         NRLMAX=NRMAX
      ENDIF

      IF(NRL == 1) THEN
         NR0 = NRL
         NR1 = NRL+1
         NR2 = NRL+2
         NR3 = NRL+3
      ELSEIF(NRL == 2) THEN
         NR0 = NRL
         NR1 = NRL-1
         NR2 = NRL+1
         NR3 = NRL+2
      ELSEIF(NRL == NRLMAX) THEN
         NR0 = NRL
         NR1 = NRL-1
         NR2 = NRL-2
         NR3 = NRL-3
      ELSE
         NR0 = NRL
         NR1 = NRL-2
         NR2 = NRL-1
         NR3 = NRL+1
      ENDIF

      DX1 = R(NR1) - R(NR0)
      DX2 = R(NR2) - R(NR0)
      DX3 = R(NR3) - R(NR0)

      DERIV4 = - F(NR1) * DX2**2 * DX3**2 * (DX2 - DX3) &
     &         - F(NR2) * DX3**2 * DX1**2 * (DX3 - DX1) &
     &         - F(NR3) * DX1**2 * DX2**2 * (DX1 - DX2) &
     &         - F(NR0) * (DX1 - DX2) * (DX2 - DX3) * (DX3 - DX1) &
     &                  * (DX1 * DX2 + DX2 * DX3 + DX3 * DX1)
      DERIV4 = DERIV4 / (DX1 * DX2 * DX3 &
     &                * (DX1 - DX2) * (DX2 - DX3) * (DX3 - DX1))

      RETURN
      END FUNCTION DERIV4

!     *******************************************
!     **  Derivative calculated from 3 points  **
!     *******************************************
!     -- This formulation has a second-order accuracy. --

      FUNCTION DERIV3P(f0, f1, f2, x10, x11, x12)

      USE task_kinds,ONLY: dp
      implicit none
      real(dp), intent(in) :: f0, f1, f2, x10, x11, x12
      real(dp):: DERIV3P
      real(dp) :: dx11, dx12

      dx11 = x11 - x10
      dx12 = x12 - x10
      DERIV3P = (dx12**2 * f1 - dx11**2 * f2) &
     &            / (dx11 * dx12 * (dx12 - dx11)) &
     &            - (dx12 + dx11) / (dx11 * dx12) * f0

      RETURN
      END FUNCTION DERIV3P

!     ** Array input version **

      FUNCTION DERIV3(NR,R,F,NRMAX,ID)

!     ID = 0    : NR = 0 to NRMAX
!          else : NR = 1 to NRMAX

      USE task_kinds,ONLY: dp
      implicit none
      integer,                    intent(in) :: NR, NRMAX, ID
      real(dp), dimension(1:NRMAX+1), intent(in) :: R, F
      real(dp):: DERIV3
      integer :: NR0, NR1, NR2
      real(dp) :: DLT, DLT1, DLT2

      IF(ID == 0) THEN
         IF(NR == 0) THEN
            NR0 = NR+1
            NR1 = NR+2
            NR2 = NR+3
         ELSEIF(NR == NRMAX) THEN
            NR0 = NR+1
            NR1 = NR
            NR2 = NR-1
         ELSE
            NR0 = NR+1
            NR1 = NR
            NR2 = NR+2
         ENDIF
      ELSE
         IF(NR == 1) THEN
            NR0 = NR
            NR1 = NR+1
            NR2 = NR+2
         ELSEIF(NR == NRMAX) THEN
            NR0 = NR
            NR1 = NR-1
            NR2 = NR-2
         ELSE
            NR0 = NR
            NR1 = NR-1
            NR2 = NR+1
         ENDIF
      END IF

      DLT1 = R(NR1) - R(NR0)
      DLT2 = R(NR2) - R(NR0)
      DERIV3 = (DLT2**2*F(NR1)-DLT1**2*F(NR2)) &
     &       /(DLT1*DLT2*(DLT2-DLT1)) &
     &       -(DLT2+DLT1)*F(NR0)/(DLT1*DLT2)

      RETURN
      END FUNCTION DERIV3

!!$      real(dp) FUNCTION DERIV3(NR,R,F,NRMAX,NRM,ID)
!!$
!!$!     ID = 0    : NR = 0 to NRMAX  ==> NR = 1 to NRMAX+1
!!$!          else : NR = 1 to NRMAX
!!$
!!$      USE task_kinds,ONLY: dp
!!$      implicit none
!!$      integer,              intent(in) :: NR, NRMAX, NRM, ID
!!$      real(dp), dimension(NRM), intent(in) :: R, F
!!$      integer :: NRL, NRLMAX, NR0, NR1, NR2
!!$      real(dp) :: DLT, DLT1, DLT2
!!$
!!$      IF(ID == 0) THEN
!!$         NRL=NR+1
!!$         NRLMAX=NRMAX+1
!!$      ELSE
!!$         NRL=NR
!!$         NRLMAX=NRMAX
!!$      ENDIF
!!$
!!$      IF(NRL == 1) THEN
!!$         NR0 = NRL
!!$         NR1 = NRL+1
!!$         NR2 = NRL+2
!!$      ELSEIF(NRL == NRLMAX) THEN
!!$         NR0 = NRL
!!$         NR1 = NRL-1
!!$         NR2 = NRL-2
!!$      ELSE
!!$         NR0 = NRL
!!$         NR1 = NRL-1
!!$         NR2 = NRL+1
!!$      ENDIF
!!$
!!$      DLT1 = R(NR1) - R(NR0)
!!$      DLT2 = R(NR2) - R(NR0)
!!$      DERIV3 = (DLT2**2*F(NR1)-DLT1**2*F(NR2)) &
!!$     &       /(DLT1*DLT2*(DLT2-DLT1)) &
!!$     &       -(DLT2+DLT1)*F(NR0)/(DLT1*DLT2)
!!$
!!$      RETURN
!!$      END FUNCTION DERIV3

!     *** Extrapolate center value *********************************
!     *  assuming the 1st-derivative is zero at the center (rho=0) *
!     *    based on a three-point formula.                         *
!     **************************************************************

      FUNCTION FCTR(R1,R2,F1,F2)

      USE task_kinds,ONLY: dp
      implicit none
      real(dp), intent(in) :: R1, R2, F1, F2
      real(dp):: FCTR

      FCTR = (R2**2*F1-R1**2*F2)/(R2**2-R1**2)

      RETURN
      END FUNCTION FCTR

!     *** Extrapolate center value *********************************
!     *  assuming the 1st-derivative is zero at the center (rho=0) *
!     *    based on a four-point formula.                          *
!     **************************************************************

      FUNCTION FCTR4pt(R1,R2,R3,F1,F2,F3) result(F0)

      USE task_kinds,ONLY: dp
      implicit none
      real(dp), intent(in) :: R1, R2, R3, F1, F2, F3
      real(dp) :: F0, h2, h3, RHSinv, LHS

      h2 = R2 - R1 ; h3 = R3 - R2

      RHSinv = (R1 * R2 * R3) / (R1 * R3 + R2 * R3 + R1 * R2)
      LHS =    R2 * R3 / (R1 * (R2 - R1) * (R3 - R1)) * F1 &
           & - R1 * R3 / (R2 * (R2 - R1) * (R3 - R2)) * F2 &
           & + R1 * R2 / (R3 * (R3 - R1) * (R3 - R2)) * F3
      F0 = LHS * RHSinv

      RETURN
      END FUNCTION FCTR4pt

!     *** Extrapolate center value *********************************
!     *  assuming the 2st-derivative is zero at the center (rho=0) *
!     *  (Higher order version)                                    *
!     **************************************************************

      FUNCTION FCTR2(R1,R2,R3,F1,F2,F3) result(F0)

      USE task_kinds,ONLY: dp
      implicit none
      real(dp), intent(in) :: R1, R2, R3, F1, F2, F3
      real(dp) :: F0, D1, D2, D3

      D1 = (R2 - R1) * (R3 - R1) * (R1 + R2 + R3)
      D2 = (R2 - R1) * (R3 - R2) * (R1 + R2 + R3)
      D3 = (R3 - R1) * (R3 - R2) * (R1 + R2 + R3)

      F0 =  R2 * R3 * (R2 + R3) / D1 * F1 &
     &    - R3 * R1 * (R3 + R1) / D2 * F2 &
     &    + R1 * R2 * (R1 + R2) / D3 * F3

      RETURN
      END FUNCTION FCTR2

!     *** Extrapolate center and its adjacent values ***************
!     *  assuming the 1st-derivative is zero at the center (rho=0) *
!     **************************************************************

      subroutine SCTR(R1,R2,R3,R4,F2,F3,F4,F0,F1)

      USE task_kinds,ONLY: dp
      implicit none
      real(dp), intent(in)  :: R1, R2, R3, R4, F2, F3, F4
      real(dp), intent(out) :: F0, F1
      real(dp) :: D2, D3, D4

      D2 = (R3 - R2) * (R4 - R2) * (R2*R3 + R3*R4 + R4*R2)
      D3 = (R3 - R2) * (R4 - R3) * (R2*R3 + R3*R4 + R4*R2)
      D4 = (R4 - R2) * (R4 - R3) * (R2*R3 + R3*R4 + R4*R2)

      F0 =  R3**2 * R4**2 / D2 * F2 &
     &    - R4**2 * R2**2 / D3 * F3 &
     &    + R2**2 * R3**2 / D4 * F4

      F1 =  (R3 - R1) * (R4 - R1) * (R3*R4 + R1*R4 + R1*R3) / D2 * F2 &
     &    - (R2 - R1) * (R4 - R1) * (R2*R4 + R4*R1 + R1*R2) / D3 * F3 &
     &    + (R2 - R1) * (R3 - R1) * (R2*R3 + R3*R1 + R1*R2) / D4 * F4

      RETURN
      end subroutine SCTR
