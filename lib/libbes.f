C     $Id$
C
C     **************************************
C
C        BESSEL FUNCTION OF THE FIRST KIND
C
C     **************************************
C
      SUBROUTINE BESSJN(X,NMAX,BJN,DBJN)
C
      INTEGER N,IACC
      REAL*8  BJN(0:NMAX),DBJN(0:NMAX),X,BIGNO,BIGNI
      PARAMETER (IACC=40,BIGNO=1.D10,BIGNI=1.D-10)
CU    USES BESSJ0,BESSJ1
      INTEGER J,JSUM,M
      REAL*8 AX,BJ,BJM,BJP,SUM,TOX,BESSJ0,BESSJ1
C
      BJN(0)=BESSJ0(X)
      BJN(1)=BESSJ1(X)
C
      AX=ABS(X)
      IF(AX.EQ.0.D0)THEN
         DO N=2,NMAX
            BJN(N)=0.D0
         ENDDO
      ELSE IF(AX.GT.DBLE(NMAX))THEN
         TOX=2.D0/AX
         BJM=BJN(0)
         BJ=BJN(1)
         DO J=1,NMAX-1
            BJP=J*TOX*BJ-BJM
            BJM=BJ
            BJ=BJP
            BJN(J+1)=BJ
         ENDDO
      ELSE
         TOX=2.D0/AX
         M=2*((NMAX+INT(SQRT(DBLE(IACC*NMAX))))/2)
         DO N=2,NMAX
            BJN(N)=0.D0
         ENDDO
         JSUM=0
         SUM=0.D0
         BJP=0.D0
         BJ=1.D0
         DO J=M,1,-1
            BJM=J*TOX*BJ-BJP
            BJP=BJ
            BJ=BJM
            IF(ABS(BJ).GT.BIGNO)THEN
               BJ=BJ*BIGNI
               BJP=BJP*BIGNI
               DO N=J+1,NMAX
                  BJN(N)=BJN(N)*BIGNI
               ENDDO
               SUM=SUM*BIGNI
            ENDIF
            IF(JSUM.NE.0) SUM=SUM+BJ
            JSUM=1-JSUM
            IF(J.GE.2.AND.J.LE.NMAX) BJN(J)=BJP
         ENDDO
         SUM=2.D0*SUM-BJ
         DO N=2,NMAX
            BJN(N)=BJN(N)/SUM
         ENDDO
      ENDIF
      DO N=3,NMAX,2
         IF(X.LT.0.D0) BJN(N)=-BJN(N)
      ENDDO
      DBJN(0)=BJN(1)
      IF(AX.EQ.0.D0)THEN
         DBJN(1)=0.5D0
         DO N=2,NMAX
            DBJN(N)=0.D0
         ENDDO
      ELSE
         DO N=1,NMAX
            DBJN(N)=BJN(N-1)-N*BJN(N)/X
         ENDDO
      ENDIF
      RETURN
      END
C
C     ************************************************
C
C        BESSEL FUNCTION OF THE FIRST KIND, ORDER 0
C
C     ************************************************
C
      FUNCTION BESSJ0(X)
C
      REAL*8 BESSJ0,X
      REAL*8 AX,XX,Z
      REAL*8 P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     &       S1,S2,S3,S4,S5,S6,Y
      SAVE   P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     &       S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,
     &-.2073370639D-5,.2093887211D-6/
      DATA Q1,Q2,Q3,Q4,Q5/-.1562499995D-1,
     &.1430488765D-3,-.6911147651D-5,.7621095161D-6,-.934945152D-7/
      DATA R1,R2,R3,R4,R5,R6/57568490574.D0,-13362590354.D0,
     &651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0/
      DATA S1,S2,S3,S4,S5,S6/57568490411.D0,1029532985.D0,9494680.718D0,
     &59272.64853D0,267.8532712D0,1.D0/
C
      IF(ABS(X).LT.8.D0)THEN
         Y=X**2
         BESSJ0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     &         /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
         AX=ABS(X)
         Z=8.D0/AX
         Y=Z**2
         XX=AX-.785398164D0
         BESSJ0=SQRT(.636619772D0/AX)
     &         *(  COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*P5))))
     &          -Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
      ENDIF
      RETURN
      END
C
C     ************************************************
C
C        BESSEL FUNCTION OF THE FIRST KIND, ORDER 1
C
C     ************************************************
C
      FUNCTION BESSJ1(X)
C
      REAL*8 BESSJ1,X
      REAL*8 AX,XX,Z
      REAL*8 P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     &       S1,S2,S3,S4,S5,S6,Y
      SAVE   P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     &       S1,S2,S3,S4,S5,S6
      DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,
     & 242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0/
      DATA S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0,
     & 18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,
     & .2457520174D-5,-.240337019D-6/
      DATA Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3,
     &     .8449199096D-5,-.88228987D-6,.105787412D-6/
C
      IF(ABS(X).LT.8.D0)THEN
        Y=X**2
        BESSJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     &          /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
         AX=ABS(X)
         Z=8.D0/AX
         Y=Z**2
         XX=AX-2.356194491D0
         BESSJ1=SQRT(.636619772D0/AX)
     &          *(  COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*P5))))
     &           -Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
     &          *SIGN(1.D0,X)
      ENDIF
      RETURN
      END
C
C     ******************************************************************
C
C                       ****** LAMBDA FUNCTION ******
C
C        MODIFIED BESSEL FUNCTION OF THE FIRST KIND, MULIPLIED BY EXP
C
C     ******************************************************************
C
      SUBROUTINE LAMBDA(N,CX,CALAM)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      DIMENSION CALAM(0:N)
      DATA ONE/1.D0/
      DATA D55/1.D-55/
C
      XA=ABS(CX)
      NA=ABS(N)
C
      IF(NA.GE.30000.OR.N.LT.0.OR.XA.GE.173.D0) THEN
         CALAM(0)=(0.D0,0.D0)
         WRITE(6,*) 'XX LAMBDA : OUT OF RANGE : N,ABS(X) = ',N,XA
         RETURN
      ENDIF
      IF(XA.LE.1.D-8) THEN
         IF(XA.LE.1.D-77) THEN
            CALAM(0)=1.D0
            DO 10 I=1,NA
               CALAM(I)=0.D0
   10       CONTINUE
         ELSE
            CALAM(0)=EXP(-CX)
            CT1=0.5D0*CX
            T2=ONE
            CT3=ONE
            DO 20 I=1,NA
               IF(ABS(CT3).LE.1.D-77*ABS(T2/CT1)) THEN
                  CALAM(I)=0.D0
               ELSE
                  CT3=CT3*CT1/T2
                  T2=T2+ONE
                  CALAM(I)=CT3*EXP(-CX)
               ENDIF
   20       CONTINUE
         ENDIF
      ELSE
         CZ=2.D0/CX
         IF(XA.GE.10.D0) THEN
            L=40
         ELSEIF(XA.GE.0.1D0) THEN
            L=25
         ELSE
            L=10
         ENDIF
         RX=XA
         NM=MAX(NA,INT(RX))+L
         CT3=0.D0
         CT2=1.D-75
         CS=0.D0
         DO 100 I=1,NM
            K=NM-I
            CT1=(K+1)*CT2*CZ+CT3
            IF(K.LE.NA) CALAM(K)=CT1
            CS=CS+CT1
            IF(ABS(CS).GT.1.D55) THEN
               CT1=CT1*D55
               CT2=CT2*D55
               CS=CS*D55
               DO 50 J=K,NA
                  CALAM(J)=CALAM(J)*D55
   50          CONTINUE
            ENDIF
            CT3=CT2
            CT2=CT1
  100    CONTINUE
         CS=CS+CS-CT1
         DO 150 J=0,NA
            CALAM(J)=CALAM(J)/CS
  150    CONTINUE
      ENDIF
      RETURN
      END
C
C     *************************************************************
C
C        MODIFIED BESSEL FUNCTION OF THE SECOND KIND
C
C     *************************************************************
C
      REAL*8 FUNCTION DKBES(N,X)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 X
      REAL*8 XABS,RECX,HALFX,W,Y,SUM,EXPX,D75,D65,
     &       AI0,AI1,AI2,AK0,AK1,AK2,NUMER,DENOM,
     &       Z
C      REAL*4 XS,YS,FK,C69,C100
      EQUIVALENCE (XABS,XS),(Y,YS)
      DATA D65,D75,C69,C100/1.0D65,1.0D75,69.D0,100.D0/
C
      XABS=X
      NABS=IABS(N)
      PI = 3.141592653589793D0
      IF (XABS) 1,1,3
    1 DKBES=D75*EXP(XABS)
      WRITE(6,2) N,X
    2 FORMAT(1H ,5X,'THE ARGUMENT OF DKBES IS INVALID.  N=',I10,
     &       '  ,  X=',D23.16)
      RETURN
C
    3 IF (NABS.LT.30000) GO TO 4
      DKBES=0.0D0
      WRITE(6,34) N,X
   34 FORMAT(1H ,5X,'THE VALUE OF DKBES IS NOT ACCURATE. N=',
     &       I10,'  ,  X=',D15.7)
      RETURN
    4 IF (XS.GT.C100) GO TO 100
      RECX=1.0D0/XABS
      Z=RECX+RECX
C                      (Z=2/X)
      IF (XS.GT.2.0) GO TO 17
      HALFX=XABS*0.5D0
      W=DLOG(HALFX)
      IF (XS.GE.0.0001) GO TO 12
C                     X LESS THAN 0.0001
      Y=HALFX*HALFX
      AI0=1.0D0+Y
      AK0=-W*AI0-5.7721566490153286D-1
     &    +4.22784335D-1*Y
      IF (NABS) 6,5,6
    5 DKBES=AK0*EXP(XABS)
      RETURN
    6 CONTINUE
      AK1=(RECX-HALFX*AK0)/AI0
      IF (NABS-1) 8,7,8
    7 AK2=AK1
      GO TO 26
    8 L=INT(-C69/W)
      IF (NABS.LE.L) GO TO 22
    9 DKBES=D75*EXP(XABS)
      WRITE(6,10) N,X
   10 FORMAT(1H ,5X,'THE VALUE OF DKBES IS OVERFLOW.  N=',
     &       I10,'  ,  X=',D23.16)
      RETURN
C              X  LESS THAN OR EQUAL TO 2.0,AND
C              GREATER THAN OR EQUAL TO 0.0001.
C              AT FIRST, I(0,X) AND I(1,X) ARE COMPUTED.
   12 CONTINUE
      AI2=0.0D0
      AI1=1.0D-75
      SUM=AI1
      K=18
      IF (XS.GE.1.0) GO TO 13
      K=14
      IF (XS.GE.0.1) GO TO 13
      K=9
   13 FK=K
C              RECURRENCE RELATION USED
C              IN DESCENDING ORDER.
      AI0=(AI1*Z)*FK+AI2
      K=K-1
      IF (K) 14,15,14
   14 SUM=SUM+AI0
      AI2=AI1
      AI1=AI0
      GO TO 13
   15 SUM=SUM+SUM+AI0
      EXPX=DEXP(XABS)
      AI0=(AI0/SUM)*EXPX
C             THE VALUE OF I(0,X) HAS BEEN OBTAIND.
C             COMPUTE K(0,X) USING POWER SERIES EXPANSION
      Y=HALFX*HALFX
C
      AK0 = ( ( ( 1.1D-17 * YS + 1.533D-15 ) * YS
     &           + 1.78593D-13 ) * YS + 1.709994D-11 ) * YS
      AK0 = ( ( AK0 + 1.31674867D-9 ) * Y + 7.9350965213D-8 )
     &     * Y
      AK0 = -W * AI0 + ( ( ( ( ( ( AK0
     & + 3.61262410320     D-6 )*Y + 1.18480393641097  D-4 )*Y
     & + 2.61478761880521  D-3 )*Y + 3.489215745643890 D-2 )*Y
     & + 2.3069608377461679D-1 )*Y + 4.2278433509846714D-1 )*Y
     & - 5.7721566490153286D-1
C
      IF (NABS) 16,5,16
   16 CONTINUE
C                  COMPUTE K(1,X) USING LOMMEL'S RELATION
      AI1 = (AI1/SUM)*EXPX
      AK1 = (RECX-AK0*AI1)/AI0
      IF (NABS.EQ.1) GO TO 7
      GO TO 22
C              X GREATER THAN 2.0
C              K(0,X) AND K(1,X) ARE COMPUTED BYE MEANS
C              OF RATIONAL EXPRESSIONS
   17 Y=RECX*0.5D0
C                      ( Y=1/2X )
      W=DEXP(-XABS) * DSQRT(3.1415926535897932D0*Y)
      IF (NABS.EQ.1) GO TO 29
      IF (XS.GE.6.0) GO TO 18
      NUMER = ( ( ( ( ( ( (
     &   9.607359468936920D-1  *Y + 2.529252282967791D+2 )*Y
     & + 7.970033517428499D+3 )*Y + 7.569208233441645D+4 )*Y
     & + 3.057001687861112D+5 )*Y + 6.316930369273204D+5 )*Y
     & + 7.490317987301367D+5 )*Y + 5.503149540522960D+5 )*Y
      NUMER = ( (   NUMER
     & + 2.642871289210102D+5 )*Y + 8.617173613052524D+4 )*Y
C
      NUMER = ( ( ( ( ( ( (         NUMER
     & + 1.959077941045564D+4 )*Y + 3.161402444515176D+3 )*Y
     & + 3.658576787591567D+2 )*Y + 3.045630347257958D+1 )*Y
     & + 1.815056853732798D 0 )*Y + 7.630840968696385D-2 )*Y
     & + 2.198080790671052D-3 )*Y
      NUMER = ( ( (   NUMER
     & + 4.108400060979278D-5 )*Y + 4.468603071316900D-7 )*Y
     & + 2.138087593931531D-9 )*Y
C
      DENOM = ( ( ( ( ( ( (                 Y
     & + 3.555555555555556D+2 )*Y + 1.513244444444444D+4 )*Y
     & + 1.956717714285714D+5 )*Y + 1.079473194053918D+6 )*Y
     & + 3.045125484052925D+6 )*Y + 4.914134724130594D+6 )*Y
     & + 4.892294125356680D+6 )*Y
      DENOM = ( ( (    DENOM
     & + 3.168987751788133D+6 )*Y + 1.388013537038700D+6 )*Y
     & + 4.227486032369927D+5 )*Y
C
      DENOM = ( ( ( ( ( ( (             DENOM
     & + 9.133105119891823D+4 )*Y + 1.418093261050334D+4 )*Y
     & + 1.593555554804435D+3 )*Y + 1.296908735314483D+2 )*Y
     & + 7.594657346992128D 0 )*Y + 3.149536504573292D-1 )*Y
     & + 8.975302543644857D-3 )*Y
      DENOM = ( (   DENOM
     & + 1.663376533185147D-4 )*Y + 1.797062622699450D-6 )*Y
     & + 8.552350375726116D-9
C
      GO TO 19
C
   18 CONTINUE
      NUMER = ( ( ( ( ( ( (
     &   9.344272845910382D-1  *Y + 7.602616403443910D+1 )*Y
     & + 7.395743295243645D+2 )*Y + 2.143299351353469D+3 )*Y
     & + 2.591022778768068D+3 )*Y + 1.552060383786633D+3 )*Y
     & + 5.070112211483554D+2 )*Y + 9.496661840062622D+1 )*Y
      NUMER = ( ( ( (    NUMER
     & + 1.035232548375690D+1 )*Y + 6.417643499456601D-1 )*Y
     & + 2.076251506438793D-2 )*Y + 2.696430527842588D-4 )*Y
C
      DENOM = ( ( ( ( ( ( (     Y + 128.0D0) * Y
     & + 1.952426666666667D+3 )*Y + 8.925379047619048D+3 )*Y
     & + 1.700072199546485D+4 )*Y + 1.598598652282462D+4 )*Y
     & + 8.186476153700395D+3 )*Y + 2.418159110016117D+3 )*Y
      DENOM = ( ( ( (     DENOM
     & + 4.239448497375430D+2 )*Y + 4.421129278670809D+1 )*Y
     & + 2.659325881907253D 0 )*Y + 8.426345399508084D-2 )*Y
     & + 1.078572211137035D-3
C
   19 AK0 = W*(1.0D0-NUMER/DENOM)
      IF (NABS.EQ.0) GO TO 5
   29 CONTINUE
C                    X IS GREATER THAN 2.0
C                    N IS GREATER THAN OR EQUAL TO 1
C                    RATIONAL APPROXIMAION IS USED TO COMPUTE
C                    K(1,X)
      IF (XS.GE.6.0) GO TO 20
C
      NUMER = ( ( ( ( ( ( (
     & + 4.612351477146149D+1  *Y + 5.662820251461119D+3 )*Y
     & + 1.267729206953107D+5 )*Y + 9.840279611657081D+5 )*Y
     & + 3.490851284631075D+6 )*Y + 6.612450857584894D+6 )*Y
     & + 7.385469550146852D+6 )*Y + 5.204090498624053D+6 )*Y
      NUMER = ( (   NUMER
     & + 2.426670289964806D+6 )*Y + 7.748520914563417D+5 )*Y
C
      NUMER = ( ( ( ( ( ( (     NUMER
     & + 1.735592435557683D+5 )*Y + 2.771333212784498D+4 )*Y
     & + 3.183321561628970D+3 )*Y + 2.636216725447433D+2 )*Y
     & + 1.565459511840883D+1 )*Y + 6.565891381439845D-1 )*Y
     & + 1.888507334035547D-2 )*Y
      NUMER = ( ( (    NUMER
     & + 3.526831364930226D-4 )*Y + 3.834684961199849D-6 )*Y
     & + 1.834777493397057D-8 )*Y
C
      DENOM = ( ( ( ( ( ( (                Y
     & + 640.0D0              )*Y + 3.242666666666667D+4 )*Y
     & + 4.565674666666667D+5 )*Y + 2.649616021768707D+6 )*Y
     & + 7.729933921057426D+6 )*Y + 1.277675028273955D+7 )*Y
     & + 1.295019033182651D+7 )*Y
      DENOM = ( ( (    DENOM      + 8.506230281115516D+6 )*Y
     & + 3.767465314819330D+6 )*Y + 1.157963565388285D+6 )*Y
C
      DENOM = ( ( ( ( ( ( (          DENOM
     & + 2.520737013090144D+5 )*Y + 3.939147947362039D+4 )*Y
     & + 4.450965515143424D+3 )*Y + 3.639711612011614D+2 )*Y
     & + 2.140312525061418D+1 )*Y + 8.908688970078743D-1 )*Y
     & + 2.547045316439757D-2 )*Y
      DENOM = ( (   DENOM
     & + 4.734225517526958D-4 )*Y + 5.128203094044773D-6 )*Y
     & + 2.446369991196075D-8
C
      GO TO 21
C
   20 CONTINUE
      NUMER = ( ( ( ( ( ( (
     &   2.727430674283096D+1  *Y + 1.138548998351567D+3 )*Y
     & + 8.476942664996772D+3 )*Y + 2.136946909366115D+4 )*Y
     & + 2.388756033308294D+4 )*Y + 1.367517138363589D+4 )*Y
     & + 4.350909161812486D+3 )*Y + 8.026823492415450D+2 )*Y
      NUMER = ( ( ( (    NUMER
     & + 8.677093809042819D+1 )*Y + 5.356622184821665D 0 )*Y
     & + 1.730209588698993D-1 )*Y + 2.247025439868823D-3 )*Y
      DENOM = ( ( ( ( ( ( (         Y + 2.304D+2    )*Y
     & + 4.183771428571428D+3 )*Y + 2.082588444444444D+4 )*Y
     & + 4.172904489795918D+4 )*Y + 4.057981194255480D+4 )*Y
     & + 2.128483799962103D+4 )*Y + 6.401009408866192D+3 )*Y
      DENOM = ( ( ( (   DENOM
     & + 1.137957228242879D+3 )*Y + 1.200020804210648D+2 )*Y
     & + 7.284240459137258D 0 )*Y + 2.325671330264231D-1 )*Y
     & + 2.996033919825096D-3
C
   21 CONTINUE
      AK1 = W*(1.0D0+NUMER/DENOM)
      IF (NABS.EQ.1) GO TO 7
C            RECURENCE RELATION IS USED IN ASCENDING
C            ORDER 10 CALCURATE K(N,X)
   22 CONTINUE
      K = 1
   23 FK=K
      IF (DABS(AK1).GE.D65) GO TO 25
      AK2=AK1*Z*FK+AK0
   24 K=K+1
      IF (K.GE.NABS) GO TO 26
      AK0=AK1
      AK1=AK2
      GO TO 23
   25 W=AK1*1.0D-10
      Y=AK0*1.0D-10
      AK2=W*Z*FK+Y
      IF (AK2.GE.D65) GO TO 9
      AK2=AK2/1.0D-10
      GO TO 24
C                 VALUE OF K(NABS,X) HAS BEEN OBTAINED
C                 IN AK2
   26 CONTINUE
      IF (X) 27,28,28
   27 K=NABS/2
      IF (K+K.NE.NABS) AK2=-AK2
   28 DKBES=AK2*EXP(XABS)
      RETURN
C
C          VALUE OF X  IS  LARGE
C
  100 AMYU=4.D0*DBLE(NABS)**2
      AMYU1=AMYU-1.D0
      AMYU2=AMYU-9.D0
      AMYU3=AMYU-25.D0
      DKBES=SQRT(PI/(2.D0*XABS))*(1.D0+AMYU1/(8.D0*XABS)
     &      +AMYU1*AMYU2/(128.D0*XABS**2)+AMYU1*AMYU2*AMYU3
     &      /(3072.D0*XABS**3))
      RETURN
      END
