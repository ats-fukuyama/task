C     $Id$
C
C
C ***********************************************************
C
C                    ERROR FUNCTION
C
C ***********************************************************
C
      FUNCTION ERF0(X)
      REAL*8 ERF,X
CU    USES GAMMP
      REAL*8 GAMMP
      IF(X.LT.0.D0)THEN
        ERF0=-GAMMP(0.5D0,X**2)
      ELSE
        ERF0= GAMMP(0.5D0,X**2)
      ENDIF
      RETURN
      END
C
      FUNCTION GAMMP(A,X)
      REAL*8 A,GAMMP,X
CU    USES GCF,GSER
      REAL*8 GAMMCF,GAMSER,GLN
      IF(X.LT.0.D0.OR.A.LE.0.D0) THEN
         WRITE(6,*) 'BAD ARGUMENTS IN GAMMP'
         GAMMP=0.D0
         RETURN
      ENDIF
      IF(X.LT.A+1.D0)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.D0-GAMMCF
      ENDIF
      RETURN
      END
C
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      INTEGER ITMAX
      REAL*8 A,GAMMCF,GLN,X,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.E-7,FPMIN=1.E-30)
CU    USES GAMMLN
      INTEGER I
      REAL*8 AN,B,C,D,DEL,H,GAMMLN
C
      GLN=GAMMLN(A)
      B=X+1.D0-A
      C=1.D0/FPMIN
      D=1.D0/B
      H=D
      DO I=1,ITMAX
        AN=-I*(I-A)
        B=B+2.D0
        D=AN*D+B
        IF(ABS(D).LT.FPMIN)D=FPMIN
        C=B+AN/C
        IF(ABS(C).LT.FPMIN)C=FPMIN
        D=1.D0/D
        DEL=D*C
        H=H*DEL
        IF(ABS(DEL-1.D0).LT.EPS)GOTO 1
      ENDDO
      WRITE(6,*) 'A TOO LARGE, ITMAX TOO SMALL IN GCF'
1     GAMMCF=EXP(-X+A*LOG(X)-GLN)*H
      RETURN
      END
C
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      INTEGER ITMAX
      REAL*8 A,GAMSER,GLN,X,EPS
      PARAMETER (ITMAX=100,EPS=3.E-7)
CU    USES GAMMLN
      INTEGER N
      REAL*8 AP,DEL,SUM,GAMMLN
C
      GLN=GAMMLN(A)
      IF(X.LE.0.D0)THEN
        IF(X.LT.0.D0) WRITE(6,*) 'X < 0 IN GSER'
        GAMSER=0.D0
        RETURN
      ENDIF
      AP=A
      SUM=1.D0/A
      DEL=SUM
      DO N=1,ITMAX
        AP=AP+1.D0
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GOTO 1
      ENDDO
      WRITE(6,*) 'A TOO LARGE, ITMAX TOO SMALL IN GSER'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
C
      FUNCTION GAMMLN(XX)
      REAL*8 GAMMLN,XX
      INTEGER J
      DOUBLE PRECISION SER,STP,TMP,X,Y,COF(6)
      SAVE COF,STP
      DATA COF,STP/76.18009172947146D0,-86.50532032941677D0,
     &             24.01409824083091D0,-1.231739572450155D0,
     &            .1208650973866179D-2,-.5395239384953D-5,
     &            2.5066282746310005D0/
C
      X=XX
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*LOG(TMP)-TMP
      SER=1.000000000190015D0
      DO J=1,6
        Y=Y+1.D0
        SER=SER+COF(J)/Y
      ENDDO
      GAMMLN=TMP+LOG(STP*SER/X)
      RETURN
      END
C
C ***********************************************************
C
C                    DERIVATIVE OF ERROR FUNCTION
C
C ***********************************************************
C
      FUNCTION ERF1(U)
C
      IMPLICIT REAL*8(A-H,O-R,T-Z)
C
      PI=3.14159265358979323846D0
      U2=U**2
      IF (U2.GT.100.D0)THEN
         ERF1=0.D0
      ELSE
         ERF1=2/SQRT(PI)*EXP(-U2)
      END IF
C
      RETURN
      END
C
C ***************************************************************
C
      SUBROUTINE DPLEG( X, N, Y, IER)
C
C ***************************************************************
C
C ***************************************************************
C *   SUBROUTINE DPLEG(X,N,Y,IER)                               *
C *   DPLEG CALCULATES LEGENDRE POLYNOMIALS P(N,X)  (DOUBLE)    *
C *      INPUT PARAMETERS:                                      *
C *        X = ARGUMENT OF P (REAL*8)                           *
C *        N = MAXIMUM DEGREE OF P (INTEGER)                    *
C *      OUTPUT DATA:                                           *
C *        Y(N) = VALUE OF P(N,X) (REAL*8)                      *
C *        IER = 0  (N.GE.0) :ERROR CONDITION CODE              *
C *            = -1 (N.LT.0)                                    *
C ***************************************************************
C
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION Y(*)
C
C     TEST OF (N)
      IF(N .LT. 0) THEN
         WRITE(6,1000) N
         IER = -1
         RETURN
      ENDIF
      IER = 0
      Y(1) = 1.D0
      IF(N .EQ. 0) RETURN
      Y(2) = X
      IF(N .EQ. 1) RETURN
      DO 20 I=2,N
         W  = X * Y(I)
         WY = W - Y(I-1)
   20 Y(I+1) = WY + W - WY / DBLE(I)
      RETURN
 1000 FORMAT(1H0,
     &       '(SUBR. DPLEG) N=',I5,', N MUST BE NON-NEGATIVE.')
      END
