C     $Id$
C
      SUBROUTINE PCGPME(AL,NL,NA,LL,D,N,N2,B,X,EPS,PARM,ITER,DD,WK,IER)
C***********************************************************************
C  PRECONDITIONED MULTIPLY CONJUGATED GRADIENT POLINOMIAL MINIMIZED    *
C  METHOD FOR FINITE ELEMENT METHOD.                                   *
C                                                                      *
C  PARAMETERS:                                                         *
C   ON ENTRY:                                                          *
C     AL     NON-ZERO ELEMENTS OF EACH ROW OF THE MATRIX A EXCLUDING   *
C            DIAGONAL ELEMENTS.                                        *
C     NL     MAXIMUNM NUMBER OF NON-ZERO ELEMENTS IN EACH ROW OF A.    *
C     NA     THE LEADING DIMENSION OF THE ARRAY AL.                    *
C     LL     COLUMN INDEX OF NON-ZERO ELEMENTS OF EACH ROW OF THE      *
C            MATRIX A EXCLUDING DISGONAL ELEMENTS.                     *
C     D      DIAGONAL ELEMENTS OF THE MATRIX A.                        *
C     N      THE ORDER OF THE MATRIX A.                                *
C     N2     N+N2 IS THE LEADING DIMENSION OF X, WK, DD. USUALLY N2=62 *
C     B      THE RIGHT HAND SIDE OF THE EQUATIONS.                     *
C     X      THE INITIAL VALUES OF THE SOLUTION. ZEROS ARE PERMITTED.  *
C     EPS    THE TOLERLANCE FOR CONVERGENCE.                           *
C     PARM   PARAMETER FOR CONVERGENCE. ZERO IS PERMITTED.             *
C     ITR    THE MAXIMUM NUMBER OF ITERATIONS.                         *
C  ON RETURN:                                                          *
C     X      THE SOLUTION VECTOR.                                      *
C     EPS    THE ERROR OF THE SOLUTION ON RETURN.                      *
C     ITR    THE NUMBER OF ITERATION ON RETURN.                        *
C  OTHERS:  WORKING PARAMETERS.                                        *
C                                                                      *
C  COPYRIGHT:      USHIRO YASUNORI       NOV. 1 1991        VER. 1     *
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AL(NA,NL),D(N),B(N),X(N+N2),DD(N+N2),WK(N+N2,5),EPS(2)
      DIMENSION LL(NA,NL)
      IF(EPS(1).LE.1.D-16) EPS(1)=1.D-8
      IF(ITER.LE.0) ITER=INT(10*(DSQRT(DBLE(N))+1))
C*    CLEAR.
      IER=0
      DO 20 J=1,5
        DO 10 I=1,N+N2
   10     WK(I,J)=0.D0
   20 CONTINUE
      DO 30 I=1,N+N2
   30   DD(I)=0.D0
      DO 33 I=N+1,N+N2
   33   X(I)=0.D0
C*    P=WK(1) , R0=WK(2) , Q=WK(3) , E=WK(4) , V=WK(5) , R=B.
C*    PM=(I-AM)...(I-A1)  GENERATE.
      CALL DECOMP(N2,AL,NL,NA,LL,D,N,B,DD,PARM,IER)
      DO 35 I=1,N
   35   WK(I,1)=B(I)
      IF(IER.NE.0) GO TO 110
C*    B=PM*B.
      CALL LDUSUB(N2,AL,NL,NA,LL,N,WK,WK,DD)
C*    Q=PM*A*X.
      CALL AXSUB(N2,AL,NL,NA,LL,N,D,X,WK(1,3))
      CALL LDUSUB(N2,AL,NL,NA,LL,N,WK(1,3),WK(1,3),DD)
C*    BN=(B,B) , P=R0=R=B-Q , C1=(R,R).
      BN=0.D0
      C1=0.D0
      DO 40 I=1,N
        BN=BN+WK(I,1)**2
        B(I)=WK(I,1)-WK(I,3)
        C1=C1+B(I)**2
        WK(I,1)=B(I)
        WK(I,2)=B(I)
   40 CONTINUE
C**   PMCGPM ITERATION.
      DO 100 K=1,ITER
C*    Q=PM*A*P.
        CALL AXSUB(N2,AL,NL,NA,LL,N,D,WK(1,1),WK(1,3))
        CALL LDUSUB(N2,AL,NL,NA,LL,N,WK(1,3),WK(1,3),DD)
C*    ALP=C1/(R0,Q).
        RQ=0.D0
        DO 50 I=1,N
   50     RQ=RQ+WK(I,2)*WK(I,3)
        ALP=C1/RQ
C*    E=R-ALP*Q.
        DO 60 I=1,N
   60     WK(I,4)=B(I)-ALP*WK(I,3)
C*    V=PM*A*E.
        CALL AXSUB(N2,AL,NL,NA,LL,N,D,WK(1,4),WK(1,5))
        CALL LDUSUB(N2,AL,NL,NA,LL,N,WK(1,5),WK(1,5),DD)
C*    AMU=(E,V)/(V,V).
        AMU=0.D0
        SS=0.D0
        DO 70 I=1,N
          AMU=AMU+WK(I,4)*WK(I,5)
          SS=SS+WK(I,5)**2
   70   CONTINUE
        AMU=AMU/SS
C*    X=X+ALP*P+AMU*E , R=E-AMU*V , RN=(R,R) , C1=(R0,R).
        RN=0.D0
        C1=0.D0
        DO 80 I=1,N
          X(I)=X(I)+ALP*WK(I,1)+AMU*WK(I,4)
          B(I)=WK(I,4)-AMU*WK(I,5)
          RN=RN+B(I)**2
          C1=C1+B(I)*WK(I,2)
   80   CONTINUE
        ERR=DSQRT(RN/BN)
        IF(ERR.LE.EPS(1)) GO TO 110
C*    BETA=C1/(AMU*RQ).
        BETA=C1/(AMU*RQ)
        DO 90 I=1,N
   90     WK(I,1)=B(I)+BETA*(WK(I,1)-AMU*WK(I,3))
  100 CONTINUE
      IER=3000
  110 CONTINUE
      ITER=K
      EPS(1)=ERR
      RETURN
      END
C
      SUBROUTINE DECOMP(N2,AL,NL,NA,LL,D,N,B,W,PARM,IER)
C*    PLU DECOMPOSITION.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AL(NA,NL),D(N),W(N+N2),LL(NA,NL),B(N)
C*    AL=AL+AU.
      IER=0
      N1=N+1
      DO 20 K=1,NL
        DO 10 I=1,N
C     *******************************************
C         Modified by A. Fukuyama for arbitrary N2.
C         IF(LL(I,K).EQ.0) LL(I,K)=N1+MOD(I,64)
C     *******************************************
          IF(LL(I,K).EQ.0) LL(I,K)=N1+MOD(I,N2)
C     *******************************************
   10   CONTINUE
   20 CONTINUE
      P1=PARM
      IF(PARM.LT.1.D0) P1=4.D0
      P2=P1+1.D0
      DO 25 J=1,NL
        DO 25 I=1,N
   25     W(I)=W(I)+DABS(AL(I,J))
      DO 35 I=1,N
   35   W(I)=P2/(MAX(W(I),D(I))+D(I)*P1)
C*    D=D*W , AL=AL*W.
      DO 40 K=1,NL
        DO 30 I=1,N
   30   AL(I,K)=AL(I,K)*W(I)
   40 CONTINUE
      DO 50 I=1,N
       D(I)=D(I)*W(I)
   50   B(I)=B(I)*W(I)
C*    SORT.
      KC=4
      IF(NL.GE.6) KC=6
      IF(NL.GE.8) KC=8
      IF(NL.GE.11) KC=11
      IF(NL.GE.15) KC=15
      IF(NL.GE.20) KC=20
      IF(NL.GE.26) KC=26
      IF(NL.GE.34) KC=34
      IF(NL.GE.44) KC=44
      IF(NL.GE.56) KC=56
      IF(NL.GE.70) KC=70
      IF(NL.GE.90) KC=90
      DO 80 K=1,KC
        DO 70 KK=K+1,NL
          DO 60 I=1,N
            IF(DABS(AL(I,K)).LT.DABS(AL(I,KK))) THEN
              WK=AL(I,K)
              IW=LL(I,K)
              AL(I,K)=AL(I,KK)
              LL(I,K)=LL(I,KK)
              AL(I,KK)=WK
              LL(I,KK)=IW
            END IF
   60     CONTINUE
   70   CONTINUE
   80 CONTINUE
      RETURN
      END
C
      SUBROUTINE AXSUB(N2,AL,NL,NA,LL,N,D,X,Y)
C*    Y=A*X.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AL(NA,NL),X(N+N2),Y(N+N2),D(N),LL(NA,NL)
C******************************************************
C MODIFIED BY A. FUKUYAMA 92/12/12 OWING TO FTNOPT BUG
C******************************************************
      DO 20 I=1,N
        Y(I)=X(I)*D(I)
   20 CONTINUE
      DO 40 J=1,NL
        DO 30 I=1,N
          Y(I)=Y(I)+AL(I,J)*X(LL(I,J))
   30   CONTINUE
   40 CONTINUE
C      DO 20 I=1,N
C        Y(I)=X(I)*D(I)+AL(I,1)*X(LL(I,1))+AL(I,2)*X(LL(I,2))
C     &                +AL(I,3)*X(LL(I,3))+AL(I,4)*X(LL(I,4))
C   20 CONTINUE
C      DO 40 J=5,NL-3,4
CVOPTION PSETUP
C        DO 30 I=1,N
C          Y(I)=Y(I)+AL(I,J)*X(LL(I,J))+AL(I,J+1)*X(LL(I,J+1))
C     &             +AL(I,J+2)*X(LL(I,J+2))+AL(I,J+3)*X(LL(I,J+3))
C   30   CONTINUE
C   40 CONTINUE
C      IF(MOD(NL,4).EQ.3) THEN
C        DO 50 I=1,N
C          Y(I)=Y(I)+AL(I,NL-2)*X(LL(I,NL-2))
C     &             +AL(I,NL-1)*X(LL(I,NL-1))
C     &             +AL(I,NL)*X(LL(I,NL))
C   50   CONTINUE
C      ELSE
C       IF(MOD(NL,4).EQ.2) THEN
C        DO 60 I=1,N
C          Y(I)=Y(I)+AL(I,NL-1)*X(LL(I,NL-1))
C     &             +AL(I,NL)*X(LL(I,NL))
C   60   CONTINUE
C       ELSE
C        IF(MOD(NL,4).EQ.1) THEN
C         DO 70 I=1,N
C   70      Y(I)=Y(I)+AL(I,NL)*X(LL(I,NL))
C        END IF
C       END IF
C      END IF
      RETURN
      END
C
      SUBROUTINE LDUSUB(N2,AL,NL,NA,LL,N,P,Q,W)
C*    Q=PM*P
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AL(NA,NL),W(N+N2),P(N+N2),Q(N+N2),LL(NA,NL)
C*    J=1
      DO 10 I=1,N
   10   W(I)=P(I)-AL(I,1)*P(LL(I,1))
C*    J=2
      DO 15 I=1,N
   15   Q(I)=W(I)-AL(I,2)*W(LL(I,2))
C*    J=3
      DO 20 I=1,N
   20   W(I)=Q(I)-AL(I,3)*Q(LL(I,3))
C*    J=4
      DO 25 I=1,N
   25   Q(I)=W(I)-AL(I,4)*W(LL(I,4))
      IF (NL.EQ.4) GO TO 100
      IF (NL.EQ.5) THEN
       DO 24 I=1,N
   24   W(I)=Q(I)-AL(I,5)*Q(LL(I,5))
       GO TO 90
      ENDIF 
C*    J=5,6
      DO 26 I=1,N
   26   W(I)=Q(I)-AL(I,5)*Q(LL(I,5))-AL(I,6)*Q(LL(I,6))
      IF(NL.EQ.6) GO TO 90
      JJ=7
      JX=0
      IF(NL.GE.8) THEN
C*    J=7,8
      DO 28 I=1,N
   28   Q(I)=W(I)-AL(I,7)*W(LL(I,7))-AL(I,8)*W(LL(I,8))
      IF(NL.EQ.8) GO TO 100
      JJ=9
      JX=1
      END IF
      IF(NL.GE.11) THEN
C*    J=9,10,11
      DO 30 I=1,N
        W(I)=Q(I)-AL(I,9)*Q(LL(I,9))-AL(I,10)*Q(LL(I,10))
     &           -AL(I,11)*Q(LL(I,11))
   30 CONTINUE
      IF(NL.EQ.11) GO TO 90
      JJ=12
      JX=0
      END IF
      IF(NL.GE.15) THEN
C*    J=12,13,14,15
      DO 32 I=1,N
        Q(I)=W(I)-(AL(I,12)*W(LL(I,12))+AL(I,13)*W(LL(I,13))
     &            +AL(I,14)*W(LL(I,14))+AL(I,15)*W(LL(I,15)))
   32 CONTINUE
      IF(NL.EQ.15) GO TO 100
      JJ=16
      JX=1
      END IF
      IF(NL.GE.20) THEN
C*    J=16,17,18,19,20
      DO 34 I=1,N
        W(I)=Q(I)-(AL(I,16)*Q(LL(I,16))+AL(I,17)*Q(LL(I,17))
     &            +AL(I,18)*Q(LL(I,18))+AL(I,19)*Q(LL(I,19))
     &            +AL(I,20)*Q(LL(I,20)))
   34 CONTINUE
      IF(NL.EQ.20) GO TO 90
      JJ=21
      JX=0
      END IF
      IF(NL.GE.26) THEN
C*    J=21,22,23,24,25,26
      DO 35 I=1,N
        Q(I)=W(I)-(AL(I,21)*W(LL(I,21))+AL(I,22)*W(LL(I,22))
     &            +AL(I,23)*W(LL(I,23))+AL(I,24)*W(LL(I,24))
     &            +AL(I,25)*W(LL(I,25))+AL(I,26)*W(LL(I,26)))
   35 CONTINUE
      IF(NL.EQ.26) GO TO 100
      JJ=27
      JX=1
      END IF
      IF(NL.GE.34) THEN
C*    J=27,28,29,30,31,32,33,34
      DO 36 I=1,N
        W(I)=Q(I)-(AL(I,27)*Q(LL(I,27))+AL(I,28)*Q(LL(I,28))
     &            +AL(I,29)*Q(LL(I,29))+AL(I,30)*Q(LL(I,30))
     &            +AL(I,31)*Q(LL(I,31))+AL(I,32)*Q(LL(I,32))
     &            +AL(I,33)*Q(LL(I,33))+AL(I,34)*Q(LL(I,34)))
   36 CONTINUE
      IF(NL.EQ.34) GO TO 90
      JJ=35
      JX=0
      END IF
      IF(NL.GE.44) THEN
C*    J=35,36,37,38,39,40,41,42,43,44
      DO 37 I=1,N
        Q(I)=W(I)-(AL(I,35)*W(LL(I,35))+AL(I,36)*W(LL(I,36))
     &            +AL(I,37)*W(LL(I,37))+AL(I,38)*W(LL(I,38))
     &            +AL(I,39)*W(LL(I,39))+AL(I,40)*W(LL(I,40))
     &            +AL(I,41)*W(LL(I,41))+AL(I,42)*W(LL(I,42))
     &            +AL(I,43)*W(LL(I,43))+AL(I,44)*W(LL(I,44)))
   37 CONTINUE
      IF(NL.EQ.44) GO TO 100
      JJ=45
      JX=1
      END IF
      IF(NL.GE.56) THEN
C*    J=45,46,47,48,49,50,51,52,53,54,55,56
      DO 38 I=1,N
        W(I)=Q(I)-(AL(I,45)*Q(LL(I,45))+AL(I,46)*Q(LL(I,46))
     &            +AL(I,47)*Q(LL(I,47))+AL(I,48)*Q(LL(I,48))
     &            +AL(I,49)*Q(LL(I,49))+AL(I,50)*Q(LL(I,50))
     &            +AL(I,51)*Q(LL(I,51))+AL(I,52)*Q(LL(I,52))
     &            +AL(I,53)*Q(LL(I,53))+AL(I,54)*Q(LL(I,54))
     &            +AL(I,55)*Q(LL(I,55))+AL(I,56)*Q(LL(I,56)))
   38 CONTINUE
      IF(NL.EQ.56) GO TO 90
      JJ=57
      JX=0
      END IF
      IF(NL.GE.70) THEN
C*    J=57,58,59,60,61,62,63,64,65,66,67,68,69,70
      DO 41 I=1,N
        Q(I)=W(I)-(AL(I,57)*W(LL(I,57))+AL(I,58)*W(LL(I,58))
     &            +AL(I,59)*W(LL(I,59))+AL(I,60)*W(LL(I,60))
     &            +AL(I,61)*W(LL(I,61))+AL(I,62)*W(LL(I,62))
     &            +AL(I,63)*W(LL(I,63))+AL(I,64)*W(LL(I,64))
     &            +AL(I,65)*W(LL(I,65))+AL(I,66)*W(LL(I,66))
     &            +AL(I,67)*W(LL(I,67))+AL(I,68)*W(LL(I,68))
     &            +AL(I,69)*W(LL(I,69))+AL(I,70)*W(LL(I,70)))
   41 CONTINUE
      IF(NL.EQ.70) GO TO 100
      JJ=71
      JX=1
      END IF
      IF(NL.GE.90) THEN
C*    J=71,72,...,89,90
      DO 42 I=1,N
        W(I)=Q(I)-(AL(I,71)*Q(LL(I,71))+AL(I,72)*Q(LL(I,72))
     &            +AL(I,73)*Q(LL(I,73))+AL(I,74)*Q(LL(I,74))
     &            +AL(I,75)*Q(LL(I,75))+AL(I,76)*Q(LL(I,76))
     &            +AL(I,77)*Q(LL(I,77))+AL(I,78)*Q(LL(I,78))
     &            +AL(I,79)*Q(LL(I,79))+AL(I,80)*Q(LL(I,80))
     &            +AL(I,81)*Q(LL(I,81))+AL(I,82)*Q(LL(I,82))
     &            +AL(I,83)*Q(LL(I,83))+AL(I,84)*Q(LL(I,84))
     &            +AL(I,85)*Q(LL(I,85))+AL(I,86)*Q(LL(I,86))
     &            +AL(I,87)*Q(LL(I,87))+AL(I,88)*Q(LL(I,88))
     &            +AL(I,89)*Q(LL(I,89))+AL(I,90)*Q(LL(I,90)))
   42 CONTINUE
      IF(NL.EQ.90) GO TO 90
      JJ=91
      JX=0
      END IF
C*    JX=1
      IF(JX.EQ.1) THEN
        DO 45 I=1,N
   45   W(I)=Q(I)
      END IF
      IF(MOD(NL-JJ,2).EQ.0) THEN
      DO 50 I=1,N
   50   Q(I)=AL(I,JJ)*W(LL(I,JJ))
      J1=1
      ELSE
      DO 60 I=1,N
   60   Q(I)=AL(I,JJ)*W(LL(I,JJ))+AL(I,JJ+1)*W(LL(I,JJ+1))
      J1=2
      END IF
      DO 80 J=JJ+J1,NL,2
CVOPTION PSETUP
        DO 70 I=1,N
   70     Q(I)=Q(I)+AL(I,J)*W(LL(I,J))+AL(I,J+1)*W(LL(I,J+1))
   80 CONTINUE
      DO 85 I=1,N
   85   Q(I)=W(I)-Q(I)
      GO TO 100
   90 CONTINUE
      DO 95 I=1,N
   95   Q(I)=W(I)
  100 CONTINUE
      RETURN
C  END OF PCGPME
      END
