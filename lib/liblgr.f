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
C     *** Extrapolate center value assuming
C         that the gradient is zero at the center (rho=0) ***
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
C
C     *** Extrapolate edge (rho=1) or arbitrary values ***
C
      FUNCTION FEDG(R0,R1,R2,F1,F2)
C
      IMPLICIT NONE
      REAL*8 FEDG,R0,R1,R2,F1,F2
C
      FEDG = ((R2-R0)*F1-(R1-R0)*F2)/(R2-R1)
C
      RETURN
      END
C
