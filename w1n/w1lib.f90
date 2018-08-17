MODULE w1lib

CONTAINS

!     ****** ARC COSINE HYPABOLIC FUNCTION ******

  FUNCTION ACOSH(X)
    USE w1comm,ONLY: rkind
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X
    REAL(rkind):: ACOSH

    ACOSH=LOG(X+SQRT(X*X-1.D0))
    RETURN
  END FUNCTION ACOSH

!     ****** LAMBDA FUNCTION ******

  SUBROUTINE LAMBDA(N,X,ALAM)
    USE w1comm,ONLY: rkind
    IMPLICIT NONE
    INTEGER,INTENT(IN):: N
    REAL(rkind),INTENT(IN):: X
    REAL(rkind),INTENT(OUT):: ALAM(0:N)
    REAL(rkind),PARAMETER:: ONE=1.D0
    REAL(rkind),PARAMETER:: D55=1.D-55
    INTEGER:: NA,I,L,NM,K,J
    REAL(rkind):: XA,T1,T2,T3,Z,RX,S

    XA=ABS(X)
    NA=ABS(N)

    IF(NA.GE.30000.OR.XA.GE.173.D0) THEN
       WRITE(6,*) 'XX LAMBDA : OUT OF RANGE : N,X = ',N,X
       RETURN
    ENDIF
    IF(XA.LE.1.D-8) THEN
       IF(XA.LE.1.D-77) THEN
          ALAM(0)=ONE
          DO I=1,NA
             ALAM(I)=0.D0
          END DO
       ELSE
          ALAM(0)=EXP(-XA)
          T1=0.5D0*XA
          T2=ONE
          T3=ONE
          DO I=1,NA
             IF(T3.LE.1.D-77*T2/T1) THEN
                ALAM(I)=0.D0
             ELSE
                T3=T3*T1/T2
                T2=T2+ONE
                ALAM(I)=T3*EXP(-XA)
             ENDIF
          END DO
       ENDIF
    ELSE
       Z=2.D0/XA
       IF(XA.GE.10.D0) THEN
          L=40
       ELSEIF(XA.GE.0.1D0) THEN
          L=25
       ELSE
          L=10
       ENDIF
       RX=XA
       NM=MAX(NA,INT(RX))+L
       T3=0.D0
       T2=1.D-75
       S=0.D0
       DO I=1,NM
          K=NM-I
          T1=(K+1)*T2*Z+T3
          IF(K.LE.NA) ALAM(K)=T1
          S=S+T1
          IF(ABS(S).GT.1.D55) THEN
             T1=T1*D55
             T2=T2*D55
             S=S*D55
             DO J=K,NA
                ALAM(J)=ALAM(J)*D55
             END DO
          ENDIF
          T3=T2
          T2=T1
       END DO
       S=S+S-T1
       DO J=0,NA
          ALAM(J)=ALAM(J)/S
       END DO
    ENDIF
    RETURN
  END SUBROUTINE LAMBDA
END MODULE w1lib
