!     $Id$

MODULE libpcgpme

  PRIVATE
  PUBLIC pcgpme

CONTAINS

  SUBROUTINE PCGPME(AL,NL,NA,LL,D,N,N2,BIN,XIN,EPSIN,ITRIN, &
                    XOUT,EPSOUT,ITROUT,IER)

!***********************************************************************
!  PRECONDITIONED MULTIPLY CONJUGATED GRADIENT POLINOMIAL MINIMIZED    *
!  METHOD FOR FINITE ELEMENT METHOD.                                   *
!                                                                      *
!  PARAMETERS:                                                         *
!   ON ENTRY:                                                          *
!     AL     NON-ZERO ELEMENTS OF EACH ROW OF THE MATRIX A EXCLUDING   *
!            DIAGONAL ELEMENTS.                                        *
!     NL     MAXIMUNM NUMBER OF NON-ZERO ELEMENTS IN EACH ROW OF A.    *
!     NA     THE LEADING DIMENSION OF THE ARRAY AL.                    *
!     LL     COLUMN INDEX OF NON-ZERO ELEMENTS OF EACH ROW OF THE      *
!            MATRIX A EXCLUDING DISGONAL ELEMENTS.                     *
!     D      DIAGONAL ELEMENTS OF THE MATRIX A.                        *
!     N      THE ORDER OF THE MATRIX A.                                *
!     N2     N+N2 IS THE LEADING DIMENSION OF X, WK, DD. USUALLY N2=62 *
!     BIN    THE RIGHT HAND SIDE OF THE EQUATIONS.                     *
!     XIN    THE INITIAL VALUES OF THE SOLUTION. ZEROS ARE PERMITTED.  *
!     EPSIN  THE TOLERLANCE FOR CONVERGENCE.                           *
!     ITRIN  THE MAXIMUM NUMBER OF ITERATIONS.                         *
!  ON RETURN:                                                          *
!     XOUT   THE SOLUTION VECTOR.                                      *
!     EPSOUT THE ERROR OF THE SOLUTION ON RETURN.                      *
!     ITROUT THE NUMBER OF ITERATION ON RETURN.                        *
!  OTHERS:  WORKING PARAMETERS.                                        *
!                                                                      *
!     PARM   PARAMETER FOR CONVERGENCE. ZERO IS PERMITTED.             *
!  COPYRIGHT:      USHIRO YASUNORI       NOV. 1 1991        VER. 1     *
!***********************************************************************

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NA,N,N2,ITRIN
      REAL(8),INTENT(IN):: EPSIN
      INTEGER,INTENT(IN):: NL
      REAL(8),DIMENSION(NA,NL),INTENT(INOUT):: AL
      INTEGER,DIMENSION(NA,NL),INTENT(INOUT):: LL
      REAL(8),DIMENSION(N),INTENT(INOUT):: D
      REAL(8),DIMENSION(N),INTENT(IN):: BIN,XIN
      REAL(8),DIMENSION(N),INTENT(OUT):: XOUt
      REAL(8),DIMENSION(N+N2):: X
      INTEGER,INTENT(OUT):: ITROUT,IER
      REAL(8),INTENT(OUT):: EPSOUT
      REAL(8),DIMENSION(N):: B
      REAL(8),DIMENSION(N+N2):: DD,TEMP
      REAL(8),DIMENSION(N+N2,5):: WK
      INTEGER:: I,J,K,ITR
      REAL(8):: PARM,BN,C1,RQ,ALP,AMU,SS,RN,ERR,BETA,EPS

      IF(EPSIN.LE.1.D-16) THEN
         EPS=1.D-8
      ELSE
         EPS=EPSIN
      END IF
      IF(ITRIN.LE.0) THEN
         ITR=INT(10*(DSQRT(DBLE(N))+1))
      ELSE
         ITR=ITRIN
      ENDIF
      PARM=0.D0

!*    CLEAR.

      IER=0
      DO J=1,5
         DO I=1,N+N2
            WK(I,J)=0.D0
         END DO
      END DO
      DO I=1,N+N2
         DD(I)=0.D0
      END DO
      DO I=1,N
         B(I)=BIN(I)
         X(I)=XIN(I)
      END DO
      DO I=N+1,N+N2
         X(I)=0.D0
      END DO

!*    P=WK(1) , R0=WK(2) , Q=WK(3) , E=WK(4) , V=WK(5) , R=B.
!*    PM=(I-AM)...(I-A1)  GENERATE.

      CALL DECOMP(N2,AL,NL,NA,LL,D,N,B,DD,PARM,IER)

      DO I=1,N
         WK(I,1)=B(I)
      END DO
      IF(IER.NE.0) GO TO 110

!*    B=PM*B.

      DO I=1,N+N2
         TEMP(I)=WK(I,1)
      END DO
      CALL LDUSUB(N2,AL,NL,NA,LL,N,WK(1,1),TEMP,DD)
      DO I=1,N
         WK(I,1)=TEMP(I)
      END DO

!*    Q=PM*A*X.

      CALL AXSUB(N2,AL,NL,NA,LL,N,D,X,WK(1,3))
      DO I=1,N+N2
         TEMP(I)=WK(I,3)
      END DO
      CALL LDUSUB(N2,AL,NL,NA,LL,N,WK(1,3),TEMP,DD)
      DO I=1,N
         WK(I,3)=TEMP(I)
      END DO

!*    BN=(B,B) , P=R0=R=B-Q , C1=(R,R).
      BN=0.D0
      C1=0.D0
      DO I=1,N
         BN=BN+WK(I,1)**2
         B(I)=WK(I,1)-WK(I,3)
         C1=C1+B(I)**2
         WK(I,1)=B(I)
         WK(I,2)=B(I)
      END DO

!**   PMCGPM ITERATION.
      DO K=1,ITR
!*    Q=PM*A*P.
         CALL AXSUB(N2,AL,NL,NA,LL,N,D,WK(1,1),WK(1,3))
         DO I=1,N+N2
            TEMP(I)=WK(I,3)
         END DO
         CALL LDUSUB(N2,AL,NL,NA,LL,N,WK(1,3),TEMP,DD)
         DO I=1,N
            WK(I,3)=TEMP(I)
         END DO
!*    ALP=C1/(R0,Q).
         RQ=0.D0
         DO I=1,N
            RQ=RQ+WK(I,2)*WK(I,3)
         END DO
         IF(RQ.EQ.0.D0) RQ=1.D0
         ALP=C1/RQ
!*    E=R-ALP*Q.
         DO I=1,N
            WK(I,4)=B(I)-ALP*WK(I,3)
         END DO
!*    V=PM*A*E.
         CALL AXSUB(N2,AL,NL,NA,LL,N,D,WK(1,4),WK(1,5))
         DO I=1,N+N2
            TEMP(I)=WK(I,5)
         END DO
         CALL LDUSUB(N2,AL,NL,NA,LL,N,WK(1,5),TEMP,DD)
         DO I=1,N
            WK(I,5)=TEMP(I)
         END DO
!*    AMU=(E,V)/(V,V).
         AMU=0.D0
         SS=0.D0
         DO I=1,N
            AMU=AMU+WK(I,4)*WK(I,5)
            SS=SS+WK(I,5)**2
         END DO
         IF(SS.EQ.0.D0) SS=1.D0
         AMU=AMU/SS
!*    X=X+ALP*P+AMU*E , R=E-AMU*V , RN=(R,R) , C1=(R0,R).
         RN=0.D0
         C1=0.D0
         DO I=1,N
            X(I)=X(I)+ALP*WK(I,1)+AMU*WK(I,4)
            B(I)=WK(I,4)-AMU*WK(I,5)
            RN=RN+B(I)**2
            C1=C1+B(I)*WK(I,2)
         END DO
         IF(BN.EQ.0.D0) BN=1.D0
         ERR=DSQRT(RN/BN)
         IF(ERR.LE.EPS) GO TO 110
!*    BETA=C1/(AMU*RQ).
         BETA=C1/(AMU*RQ)
         DO I=1,N
            WK(I,1)=B(I)+BETA*(WK(I,1)-AMU*WK(I,3))
         END DO
      END DO
      IER=3000
  110 CONTINUE
      DO I=1,N
         XOUT(I)=X(I)
      END DO
      ITROUT=K
      EPSOUT=ERR
      RETURN
    END SUBROUTINE PCGPME
!
    SUBROUTINE DECOMP(N2,AL,NL,NA,LL,D,N,B,W,PARM,IER)

!*    PLU DECOMPOSITION.
      IMPLICIT NONE
      INTEGER,INTENT(IN):: N2,NL,NA,N
      REAL(8),INTENT(IN):: PARM
      INTEGER,DIMENSION(NA,NL),INTENT(INOUT):: LL
      REAL(8),DIMENSION(NA,NL),INTENT(INOUT):: AL
      REAL(8),DIMENSION(N),INTENT(INOUT):: D,B
      REAL(8),DIMENSION(N+N2),INTENT(INOUT):: W
      INTEGER,INTENT(OUT):: IER
      INTEGER:: N1,I,J,K,KC,KK,IW
      REAL(8):: P1,P2,WK

!*    AL=AL+AU.
      IER=0
      N1=N+1
      DO K=1,NL
         DO I=1,N
!     *******************************************
!         Modified by A. Fukuyama for arbitrary N2.
!         IF(LL(I,K).EQ.0) LL(I,K)=N1+MOD(I,64)
!     *******************************************
            IF(LL(I,K).EQ.0) LL(I,K)=N1+MOD(I,N2)
!     *******************************************
         END DO
      END DO
      P1=PARM
      IF(PARM.LT.1.D0) P1=4.D0
      P2=P1+1.D0
      DO J=1,NL
         DO I=1,N
            W(I)=W(I)+DABS(AL(I,J))
         END DO
      END DO
      DO I=1,N
         W(I)=P2/(MAX(W(I),D(I))+D(I)*P1)
      END DO
!*    D=D*W , AL=AL*W.
      DO K=1,NL
         DO I=1,N
            AL(I,K)=AL(I,K)*W(I)
         END DO
      END DO
      DO I=1,N
         D(I)=D(I)*W(I)
         B(I)=B(I)*W(I)
      END DO
!*    SORT.
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
      DO K=1,KC
         DO KK=K+1,NL
            DO I=1,N
               IF(DABS(AL(I,K)).LT.DABS(AL(I,KK))) THEN
                  WK=AL(I,K)
                  IW=LL(I,K)
                  AL(I,K)=AL(I,KK)
                  LL(I,K)=LL(I,KK)
                  AL(I,KK)=WK
                  LL(I,KK)=IW
               END IF
            END DO
         END DO
      END DO
      RETURN
    END SUBROUTINE DECOMP
!
    SUBROUTINE AXSUB(N2,AL,NL,NA,LL,N,D,X,Y)

!*    Y=A*X.
      IMPLICIT NONE
      INTEGER,INTENT(IN):: N2,NL,NA,N
      INTEGER,DIMENSION(NA,NL),INTENT(IN):: LL
      REAL(8),DIMENSION(NA,NL),INTENT(IN):: AL
      REAL(8),DIMENSION(N),INTENT(IN):: D
      REAL(8),DIMENSION(N+N2),INTENT(IN):: X
      REAL(8),DIMENSION(N+N2),INTENT(OUT):: Y
      INTEGER:: I,J

!******************************************************
! MODIFIED BY A. FUKUYAMA 92/12/12 OWING TO FTNOPT BUG
!******************************************************
      DO I=1,N
         Y(I)=X(I)*D(I)
      END DO
      DO J=1,NL
         DO I=1,N
            Y(I)=Y(I)+AL(I,J)*X(LL(I,J))
         END DO
      END DO

!      DO 20 I=1,N
!        Y(I)=X(I)*D(I)+AL(I,1)*X(LL(I,1))+AL(I,2)*X(LL(I,2))
!     &                +AL(I,3)*X(LL(I,3))+AL(I,4)*X(LL(I,4))
!   20 CONTINUE
!      DO 40 J=5,NL-3,4
!VOPTION PSETUP
!        DO 30 I=1,N
!          Y(I)=Y(I)+AL(I,J)*X(LL(I,J))+AL(I,J+1)*X(LL(I,J+1))
!     &             +AL(I,J+2)*X(LL(I,J+2))+AL(I,J+3)*X(LL(I,J+3))
!   30   CONTINUE
!   40 CONTINUE
!      IF(MOD(NL,4).EQ.3) THEN
!        DO 50 I=1,N
!          Y(I)=Y(I)+AL(I,NL-2)*X(LL(I,NL-2))
!     &             +AL(I,NL-1)*X(LL(I,NL-1))
!     &             +AL(I,NL)*X(LL(I,NL))
!   50   CONTINUE
!      ELSE
!       IF(MOD(NL,4).EQ.2) THEN
!        DO 60 I=1,N
!          Y(I)=Y(I)+AL(I,NL-1)*X(LL(I,NL-1))
!     &             +AL(I,NL)*X(LL(I,NL))
!   60   CONTINUE
!       ELSE
!        IF(MOD(NL,4).EQ.1) THEN
!         DO 70 I=1,N
!   70      Y(I)=Y(I)+AL(I,NL)*X(LL(I,NL))
!        END IF
!       END IF
!      END IF
      RETURN
    END SUBROUTINE AXSUB
!
    SUBROUTINE LDUSUB(N2,AL,NL,NA,LL,N,P,Q,W)

!*    Q=PM*P
      IMPLICIT NONE
      INTEGER,INTENT(IN):: N2,NL,NA,N
      INTEGER,DIMENSION(NA,NL),INTENT(IN):: LL
      REAL(8),DIMENSION(NA,NL),INTENT(IN):: AL
      REAL(8),DIMENSION(N+N2),INTENT(IN):: P
      REAL(8),DIMENSION(N+N2),INTENT(OUT):: W,Q
      INTEGER:: I,J,JJ,JX,J1

!*    J=1
      DO I=1,N
         W(I)=P(I)-AL(I,1)*P(LL(I,1))
      END DO
!*    J=2
      DO I=1,N
         Q(I)=W(I)-AL(I,2)*W(LL(I,2))
      END DO
      IF (NL.EQ.2) GO TO 100
!*    J=3
      DO I=1,N
         W(I)=Q(I)-AL(I,3)*Q(LL(I,3))
      END DO
      IF (NL.EQ.3) GO TO 100
!*    J=4
      DO I=1,N
         Q(I)=W(I)-AL(I,4)*W(LL(I,4))
      END DO
      IF (NL.EQ.4) GO TO 100
      IF (NL.EQ.5) THEN
         DO I=1,N
            W(I)=Q(I)-AL(I,5)*Q(LL(I,5))
         END DO
         GO TO 90
      ENDIF 
!*    J=5,6
      DO I=1,N
         W(I)=Q(I)-AL(I,5)*Q(LL(I,5))-AL(I,6)*Q(LL(I,6))
      END DO
      IF(NL.EQ.6) GO TO 90
      JJ=7
      JX=0
      IF(NL.GE.8) THEN
!*    J=7,8
      DO I=1,N
         Q(I)=W(I)-AL(I,7)*W(LL(I,7))-AL(I,8)*W(LL(I,8))
      END DO
      IF(NL.EQ.8) GO TO 100
      JJ=9
      JX=1
      END IF
      IF(NL.GE.11) THEN
!*    J=9,10,11
         DO I=1,N
            W(I)=Q(I)-AL(I,9)*Q(LL(I,9))-AL(I,10)*Q(LL(I,10)) &
                 -AL(I,11)*Q(LL(I,11))
         END DO
         IF(NL.EQ.11) GO TO 90
         JJ=12
         JX=0
      END IF
      IF(NL.GE.15) THEN
!*    J=12,13,14,15
         DO I=1,N
            Q(I)=W(I)-(AL(I,12)*W(LL(I,12))+AL(I,13)*W(LL(I,13)) &
                 +AL(I,14)*W(LL(I,14))+AL(I,15)*W(LL(I,15)))
         END DO
         IF(NL.EQ.15) GO TO 100
         JJ=16
         JX=1
      END IF
      IF(NL.GE.20) THEN
!*    J=16,17,18,19,20
         DO I=1,N
            W(I)=Q(I)-(AL(I,16)*Q(LL(I,16))+AL(I,17)*Q(LL(I,17)) &
                 +AL(I,18)*Q(LL(I,18))+AL(I,19)*Q(LL(I,19)) &
                 +AL(I,20)*Q(LL(I,20)))
         END DO
         IF(NL.EQ.20) GO TO 90
         JJ=21
         JX=0
      END IF
      IF(NL.GE.26) THEN
!*    J=21,22,23,24,25,26
         DO I=1,N
            Q(I)=W(I)-(AL(I,21)*W(LL(I,21))+AL(I,22)*W(LL(I,22)) &
                 +AL(I,23)*W(LL(I,23))+AL(I,24)*W(LL(I,24)) &
                 +AL(I,25)*W(LL(I,25))+AL(I,26)*W(LL(I,26)))
         END DO
         IF(NL.EQ.26) GO TO 100
         JJ=27
         JX=1
      END IF
      IF(NL.GE.34) THEN
!*    J=27,28,29,30,31,32,33,34
         DO I=1,N
            W(I)=Q(I)-(AL(I,27)*Q(LL(I,27))+AL(I,28)*Q(LL(I,28)) &
                 +AL(I,29)*Q(LL(I,29))+AL(I,30)*Q(LL(I,30)) &
                 +AL(I,31)*Q(LL(I,31))+AL(I,32)*Q(LL(I,32)) &
                 +AL(I,33)*Q(LL(I,33))+AL(I,34)*Q(LL(I,34)))
         END DO
         IF(NL.EQ.34) GO TO 90
         JJ=35
         JX=0
      END IF
      IF(NL.GE.44) THEN
!*    J=35,36,37,38,39,40,41,42,43,44
         DO I=1,N
            Q(I)=W(I)-(AL(I,35)*W(LL(I,35))+AL(I,36)*W(LL(I,36)) &
                 +AL(I,37)*W(LL(I,37))+AL(I,38)*W(LL(I,38)) &
                 +AL(I,39)*W(LL(I,39))+AL(I,40)*W(LL(I,40)) &
                 +AL(I,41)*W(LL(I,41))+AL(I,42)*W(LL(I,42)) &
                 +AL(I,43)*W(LL(I,43))+AL(I,44)*W(LL(I,44)))
         END DO
         IF(NL.EQ.44) GO TO 100
         JJ=45
         JX=1
      END IF
      IF(NL.GE.56) THEN
!*    J=45,46,47,48,49,50,51,52,53,54,55,56
         DO I=1,N
            W(I)=Q(I)-(AL(I,45)*Q(LL(I,45))+AL(I,46)*Q(LL(I,46)) &
                 +AL(I,47)*Q(LL(I,47))+AL(I,48)*Q(LL(I,48)) &
                 +AL(I,49)*Q(LL(I,49))+AL(I,50)*Q(LL(I,50)) &
                 +AL(I,51)*Q(LL(I,51))+AL(I,52)*Q(LL(I,52)) &
                 +AL(I,53)*Q(LL(I,53))+AL(I,54)*Q(LL(I,54)) &
                 +AL(I,55)*Q(LL(I,55))+AL(I,56)*Q(LL(I,56)))
         END DO
         IF(NL.EQ.56) GO TO 90
         JJ=57
         JX=0
      END IF
      IF(NL.GE.70) THEN
!*    J=57,58,59,60,61,62,63,64,65,66,67,68,69,70
         DO I=1,N
            Q(I)=W(I)-(AL(I,57)*W(LL(I,57))+AL(I,58)*W(LL(I,58)) &
                 +AL(I,59)*W(LL(I,59))+AL(I,60)*W(LL(I,60)) &
                 +AL(I,61)*W(LL(I,61))+AL(I,62)*W(LL(I,62)) &
                 +AL(I,63)*W(LL(I,63))+AL(I,64)*W(LL(I,64)) &
                 +AL(I,65)*W(LL(I,65))+AL(I,66)*W(LL(I,66)) &
                 +AL(I,67)*W(LL(I,67))+AL(I,68)*W(LL(I,68)) &
                 +AL(I,69)*W(LL(I,69))+AL(I,70)*W(LL(I,70)))
         END DO
         IF(NL.EQ.70) GO TO 100
         JJ=71
         JX=1
      END IF
      IF(NL.GE.90) THEN
!*    J=71,72,...,89,90
         DO I=1,N
            W(I)=Q(I)-(AL(I,71)*Q(LL(I,71))+AL(I,72)*Q(LL(I,72)) &
                 +AL(I,73)*Q(LL(I,73))+AL(I,74)*Q(LL(I,74)) &
                 +AL(I,75)*Q(LL(I,75))+AL(I,76)*Q(LL(I,76)) &
                 +AL(I,77)*Q(LL(I,77))+AL(I,78)*Q(LL(I,78)) &
                 +AL(I,79)*Q(LL(I,79))+AL(I,80)*Q(LL(I,80)) &
                 +AL(I,81)*Q(LL(I,81))+AL(I,82)*Q(LL(I,82)) &
                 +AL(I,83)*Q(LL(I,83))+AL(I,84)*Q(LL(I,84)) &
                 +AL(I,85)*Q(LL(I,85))+AL(I,86)*Q(LL(I,86)) &
                 +AL(I,87)*Q(LL(I,87))+AL(I,88)*Q(LL(I,88)) &
                 +AL(I,89)*Q(LL(I,89))+AL(I,90)*Q(LL(I,90)))
         END DO
         IF(NL.EQ.90) GO TO 90
         JJ=91
         JX=0
      END IF
!*    JX=1
      IF(JX.EQ.1) THEN
        DO  I=1,N
           W(I)=Q(I)
        END DO
      END IF
      IF(MOD(NL-JJ,2).EQ.0) THEN
         DO I=1,N
            Q(I)=AL(I,JJ)*W(LL(I,JJ))
         END DO
         J1=1
      ELSE
         DO I=1,N
            Q(I)=AL(I,JJ)*W(LL(I,JJ))+AL(I,JJ+1)*W(LL(I,JJ+1))
         END DO
         J1=2
      END IF
      DO J=JJ+J1,NL,2
!VOPTION PSETUP
         DO I=1,N
            Q(I)=Q(I)+AL(I,J)*W(LL(I,J))+AL(I,J+1)*W(LL(I,J+1))
         END DO
      END DO
      DO I=1,N
         Q(I)=W(I)-Q(I)
      END DO
      GO TO 100
   90 CONTINUE
      DO I=1,N
         Q(I)=W(I)
      END DO
  100 CONTINUE
      RETURN
    END SUBROUTINE LDUSUB
  END MODULE libpcgpme
