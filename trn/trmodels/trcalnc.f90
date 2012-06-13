MODULE trcalnc
! -----------------------------------------------------------------------
!  Subprograms included in this module are refered from models
!   which calculates neoclassical transport coefficient.
! -----------------------------------------------------------------------

  USE trcomm, ONLY: ikind, rkind

  PRIVATE
  PUBLIC ftpf  

CONTAINS

! from trcalc.f90 in TASK/TR previous version
! **********************************************************************
! -----------------------------------------------------------------------
!           TRAPPED PARTICLE FRACTION
! -----------------------------------------------------------------------
  REAL(8) FUNCTION FTPF(ID,EPS)

    IMPLICIT NONE
    
    INTEGER(4)                  :: ID, IERR, N, I
    INTEGER(4),PARAMETER        :: IMAX=20
    REAL(8)                     :: EPS, EPSC, FTLL, FTUL, OMEGA, PI, S
    REAL(8),DIMENSION(IMAX,IMAX):: TABLE
!    REAL(8)                     :: FTL, FTU
!    EXTERNAL FTL, FTU

      IF(ID.EQ.1) THEN
!  Y. R. Lin-Liu and R. L. Miller, PoP 2 1666 (1995), eqs(7)(13)(18)(19)
         PI=3.14159265358979323846D0
         EPSC=1.D-9
         I=IMAX
         FTUL=1.D0-(1.D0-1.5D0*SQRT(EPS)+0.5D0*EPS**1.5D0)/SQRT(1-EPS**2)
         CALL RMBRG(0.D0,2.D0*PI,EPSC,S,I,N,IERR,TABLE,EPS,FTL)
         FTLL=1.D0-(1.D0-EPS)**1.5D0/SQRT(1.D0+EPS)*(S/(2.D0*PI))
         OMEGA=(3.D0*SQRT(2.D0)/2.D0*0.69D0-3.D0*SQRT(2.D0)/PI)/(1.5D0-3.D0*SQRT(2.D0)/PI)
         FTPF=OMEGA*FTUL+(1.D0-OMEGA)*FTLL

      ELSEIF(ID.EQ.2) THEN
!  S. P. Hirshman et al., NF 17 611 (1977)
         FTPF=1.D0-(1.D0-EPS)**2.D0/(DSQRT(1.D0-EPS**2)*(1.D0+1.46D0*DSQRT(EPS)))

      ELSEIF(ID.EQ.3) THEN
!  Y. R. Lin-Liu and R. L. Miller, PoP 2 1666 (1995), eqs(16)(17)(18)
         PI=3.14159265358979323846D0
         FTUL=1.5D0*SQRT(EPS)
         FTLL=3.D0*SQRT(2.D0)/PI*SQRT(EPS)
         FTPF=0.75D0*FTUL+0.25D0*FTLL

      ELSEIF(ID.EQ.4) THEN
!  M. N. Rosenbluth et al., PoF 15 116 (1972)
         FTPF=1.46D0*SQRT(EPS)

      ELSE
!  Y. B. Kim et al., PoF B 3 2050 (1991) eq(C18), default
         FTPF=1.46D0*SQRT(EPS)-0.46D0*(EPS)**1.5D0
      ENDIF

      RETURN      
      END FUNCTION FTPF

! -----------------------------------------------------------------------
!           ROMBERG INTEGRATION METHOD
! -----------------------------------------------------------------------
      SUBROUTINE RMBRG(A,B,EPS,S,IMAX,N,IERR,T,ARG,F)

!     <input>
!        A     : lower bound
!        B     : upper bound
!        EPS   : stopping criterion
!        IMAX  : maximum division number
!        F     : formula of integrand
!     <output>
!        S     : integration value
!        N     : division number
!        IERR  : error indicator
!        T     : Romberg T table

      IMPLICIT NONE
      REAL(8),INTENT(IN)      ::  A, B, ARG, EPS
      REAL(8),INTENT(OUT)     ::  S
      INTEGER(4),INTENT(INOUT):: IMAX
      INTEGER(4),INTENT(OUT)  :: IERR, N
      REAL(8),DIMENSION(IMAX,IMAX),INTENT(OUT)::  T
      INTEGER(4)::  I, J, K, N2
      REAL(8)   ::  F, X, H, S1, Y1, Y2

      EXTERNAL F

      DO K=1,IMAX
         N=2**(K-1)
         N2=N/2
         H=(B-A)/N
         Y1=0
         IF(N.EQ.1) THEN
            Y2=(F(A,ARG)+F(B,ARG))/2
         ELSE
            DO I=1,N2
               X=A+(2*I-1)*H
               Y1=Y1+F(X,ARG)
            ENDDO
            Y2=Y2+Y1
            Y1=0.D0
         ENDIF
         S=H*Y2
         T(K,1)=S
         IF(K.LE.1) GOTO 10
         DO J=2,K
            T(K,J)=T(K,J-1)+(T(K,J-1)-T(K-1,J-1))/(4**(J-1)-1)
         ENDDO
         S=T(K,K)
         S1=T(K,K-1)
         IF(ABS(S-S1).LT.EPS) THEN
            IERR=0
            IMAX=K
            RETURN
         ENDIF
 10      CONTINUE
      ENDDO
      IERR=1

      RETURN
      END SUBROUTINE RMBRG

! -----------------------------------------------------------------------
!           FUNCTION FOR ROMBERG INTEGRATION
! -----------------------------------------------------------------------
      REAL(8) FUNCTION FTU(X,EPS)
      IMPLICIT NONE
      REAL(8):: EPS, X

      FTU = X/SQRT(1.D0-X*(1.D0-EPS))

      RETURN
      END FUNCTION FTU


      REAL(8) FUNCTION FTL(X,EPS)
      IMPLICIT NONE
      REAL(8):: X, EPS
      REAL(8):: H

      H = (1.D0 - EPS) / (1.D0 + EPS * COS(X))
      FTL = (1.D0 - SQRT(1.D0 - H) * (1.D0 + 0.5D0 * H)) / H**2

      RETURN
      END FUNCTION FTL

! ***********************************************************************

END MODULE trcalnc
