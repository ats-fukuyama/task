C
C     ****** TASK/W1 (SUBROUTINE LIBRARY 3) ******
C
C              V2.10 : 91/06/08
C
      FUNCTION ERF1(U)
C
      USE libspf, ONLY: derf
      IMPLICIT NONE
      REAL(8):: U,ERF1
C
      ERF1=derf(U)
C
      RETURN
      END
C
      FUNCTION ERFC1(U)
C
      USE libspf, ONLY: derfc
      IMPLICIT NONE
      REAL(8):: U,ERFC1
C
      ERFC1=ERFC(U)
C
      RETURN
      END
C
C     ****** MACHINE DEPENDENT ROUTINE FOR CPUTIME ******
C
      SUBROUTINE FCLOCK(T)
C
      REAL*8 T
C
      CALL GUTIME(ST)
      T=ST
      RETURN
      END
C
C     ****** ARC COSINE HYPABOLIC FUNCTION ******
C
      DOUBLE PRECISION FUNCTION ACOSH(X)
C
      REAL*8 X
C
      ACOSH=LOG(X+SQRT(X*X-1.D0))
      RETURN
      END
C
C     ******* CUBIC ROOT *******
C
	FUNCTION DCBRT(X)
C
	REAL*8 DCBRT,X
C
	IF(X.GE.0.D0) THEN
	   DCBRT = X**(1.D0/3.D0)
	ELSE
	   DCBRT = -(-X)**(1.D0/3.D0)
	ENDIF
	RETURN
	END
C
C     ******* BESSEL SERIES1 *******
C
      DOUBLE PRECISION FUNCTION DJBES(N,X)
      REAL*8 X,Z,BJ,QJ,S,T1,T2,T3,ONE,D55
      DATA ONE,D55/1.D0,1.D-55/
C
      NN=IABS(N)
      XA=DABS(X)
      IF(NN.GE.30000) GOTO 900
      IF(XA.GE.30000.0) GOTO 900
      IF(XA.LT.0.00002) THEN
         T1=0.5D0*X
         IF(NN.EQ.0) THEN
            BJ=ONE-T1*T1
         ELSE
            IF(XA.LE.1.0D-77) GOTO 500
            T2=ONE
            T3=ONE
            DO 5 I=1,NN
               IF(DABS(T3).LT.1.0D-77*DABS(T2/T1)) GOTO 500
               T3=T3*T1/T2
               T2=T2+ONE
    5       CONTINUE
            BJ=T3*(ONE-T1*T1/T2)
         ENDIF
      ELSE
         IF(XA.LT.10.0) THEN
            IF(XA.GT.1.0) THEN
               L=1.4*XA+14.0
            ELSE
               L=14
            ENDIF
         ELSE
            IF(XA.GE.100.0) THEN
               L=0.073*XA+47.0
            ELSE
               L=0.27*XA+27.0
            ENDIF
         ENDIF
         Z=2.0D0/X
         NM=MAX(NN,INT(XA))+L
         T3=0.0D0
         T2=1.0D-75
         S=0.0D0
         IF(MOD(NM,2).EQ.0) NM=NM+1
         DO 100 I=1,NM,2
            K=NM-I+1
            T1=DFLOAT(K+1)*T2*Z-T3
            IF(NN.EQ.K) QJ=T1
            K=K-1
            T3=T2
            T2=T1
            T1=DFLOAT(K+1)*T2*Z-T3
            IF(NN.EQ.K) QJ=T1
            S=S+T1
            IF(DABS(S).GE.1.0D55) THEN
               T1=T1*D55
               T2=T2*D55
               S=S*D55
               IF(NN.GE.K) QJ=QJ*D55
            ENDIF
            T3=T2
            T2=T1
  100    CONTINUE
         S=S+S-T1
         BJ=QJ/S
      ENDIF
      IF(N.LT.0.AND.MOD(NN,2).NE.0) THEN
         DJBES=-BJ
      ELSE
         DJBES=BJ
      ENDIF
      RETURN
C
  900 WRITE(6,1001) N,X
  500 DJBES=0.0D0
      RETURN
C
 1001 FORMAT(1H ,5X,'THE VALUE OF DJBES IS NOT ACCURATE.  N=',I7,'
     &         , X=',D23.16)
      END
C
C     ****************************************************
C
C        FINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
C                    (-1.D0, +1.D0)
C         INTEGRAND SHOULD BE DEFINED BY FUNC(X,1-X,1+X)
C
C
      SUBROUTINE DEFT(FUNC,CS,ES,H0,EPS,ILST)
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER KL*1
      EXTERNAL FUNC
      DATA HP/1.5707 96326 79489 66192D0/
C
      EPS1=EPS**0.75
      H=H0
      X=0.D0
      CSI=HP*FUNC(X,1.D0-X,1.D0+X)
      CS=H*CSI
      CSP=0.D0
      NP=0
      NM=0
      NPMIN=1
      NMMIN=1
C
    5 IND=0
      ATP=ABS(CSI)
      ATM=ATP
      NPD=2
      IF(NP.EQ.0) NPD=1
      NMD=2
      IF(NM.EQ.0) NMD=1
C
   10 IF(IND.NE.1) THEN
         IF(NP.EQ.NPMIN+2) NPD=1
         NP=NP+NPD
         HN=DBLE(NP)*H
         HC=HP*H*COSH(HN)
         HS=HP*SINH(-HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=EXP(-HS)*CC
         XP=EXP( HS)*CC
         CT=HC*FUNC(X,XM,XP)*CC*CC
         CS=CS+CT
         AT=ATP
         ATP=ABS(CT)/H
         IF(NP.GE.NPMIN) THEN
            IF(AT+ATP.LE.EPS1*ABS(CS)) THEN
               IF(IND.EQ.-1) GO TO 100
               IND=1
            ENDIF
         ENDIF
      ENDIF
C
      IF(IND.NE.-1) THEN
         IF(NM.EQ.NMMIN+2) NMD=1
         NM=NM+NMD
         HN=DBLE(NM)*H
         HC=HP*H*COSH(HN)
         HS=HP*SINH( HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=EXP(-HS)*CC
         XP=EXP( HS)*CC
         CT=HC*FUNC(X,XM,XP)*CC*CC
         CS=CS+CT
         AT=ATM
         ATM=ABS(CT)/H
         IF(NM.GE.NMMIN) THEN
            IF(AT+ATM.LE.EPS1*ABS(CS)) THEN
               IF(IND.EQ.1) GO TO 100
               IND=-1
            ENDIF
         ENDIF
      ENDIF
      GO TO 10
C
  100 ES=ABS(CS-CSP)
      IF(ILST.EQ.2) THEN
         IF(H.GE.H0) WRITE(6,601) H,NP,NM,CS
         IF(H.LT.H0) WRITE(6,602) H,NP,NM,CS,ES
      ENDIF
      CSP=CS
      IF(ES.LE.EPS1*ABS(CS)) GO TO 200
      NMAX=MAX0(NP,NM)
      IF(NMAX.GT.1000) THEN
  110    WRITE(6,603)
         READ(5,501,END=120) KL
         IF(KL.EQ.'Q') GOTO 200
         IF(KL.NE.'C') GOTO 110
      ENDIF
  120 CONTINUE
      H=0.5D0*H
      CS=0.5D0*CS
      NPMIN=NP*2-1
      NMMIN=NM*2-1
      NP=-1
      NM=-1
      GO TO 5
C
  200 IF(ILST.EQ.1) THEN
         WRITE(6,602) H,NP,NM,CS,ES
      ENDIF
      RETURN
C
  501 FORMAT(A1)
  601 FORMAT(1H ,1PD13.5,2I8,1PD24.15)
  602 FORMAT(1H ,1PD13.5,2I8,1PD24.15,1PD14.5)
  603 FORMAT(1H ,'# SUB DEFT # C or CR : CONTINUE / Q : QUIT')
      END
C
C     ****** SOLUTION OF BAND MATRIX (GAUSSIAN ELIMINATION) ******
C
      SUBROUTINE BANDCD( A , X , N , L , LA , IERR )
C
C     COMPLEX * 16    A( LA , N ) , X( N ) , ATMP( LMAX ) , TEMP
      COMPLEX * 16    A( LA , N ) , X( N ) , TEMP
      REAL    *  8    EPS , ABS1 , ABS2
      DATA EPS/ 1.D-70 /
C
      IF( MOD(L,2) .EQ. 0 ) GO TO 9000
      LH  = (L+1)/2
      LHM = LH-1
      NM  = N -1
C
      DO 30 K = 1 , LHM
         LHMK = LH-K
         NPMK = N+1-K
         DO 30 I = 1 , LHMK
            LPMI = L+1-I
            DO 40 J = 2 , L
               A( J-1 , K ) = A( J , K )
   40       CONTINUE
            A( L    , K    ) = ( 0.D0 , 0.D0 )
            A( LPMI , NPMK ) = ( 0.D0 , 0.D0 )
   30 CONTINUE
C
C     DO 50 I = 1 , NM
C        IPIVOT = I
C        IP     = I+1
C        DO 60 K = IP , LH
C           IF( CDABS(A(1,K)) .GT. CDABS(A(1,IPIVOT)) ) IPIVOT=K
C  60    CONTINUE
C
C
      DO 50 I = 1 , NM
         IPIVOT = I
         IP     = I+1
         ABS2   = CDABS( A(1,IPIVOT) )
         DO 60 K = IP , LH
            ABS1 = CDABS( A(1,K) )
            IF( ABS1 .GT. ABS2 ) THEN
                IPIVOT = K
                ABS2 = ABS1
            ENDIF
   60    CONTINUE
C
         IF( CDABS(A(1,IPIVOT)) .LT. EPS ) GO TO 9002
         IF( IPIVOT .NE. I ) THEN
            TEMP        = X( I      )
            X( I      ) = X( IPIVOT )
            X( IPIVOT ) = TEMP
            DO 90 J = 1 , L
C              ATMP( J )          = A   ( J , I      )
               TEMP               = A   ( J , I      )
               A   ( J , I      ) = A   ( J , IPIVOT )
C              A   ( J , IPIVOT ) = ATMP( J )
               A   ( J , IPIVOT ) = TEMP
   90       CONTINUE
         END IF
C
         TEMP   = 1.D0   / A( 1 , I )
         X( I ) = X( I ) * TEMP
C
         DO 120 J = 2 , L
            A( J , I ) = A( J , I ) * TEMP
  120    CONTINUE
C
         DO 130 K = IP , LH
            TEMP   = A( 1 , K )
            X( K ) = X( K ) - X( I ) * TEMP
            DO 140 J = 2 , L
               A( J-1 , K ) = A( J , K ) - A( J , I ) * TEMP
  140       CONTINUE
C
            A( L , K ) = ( 0.D0 , 0.D0 )
  130    CONTINUE
         IF( LH .LT. N ) LH = LH + 1
   50 CONTINUE
C
      X( N ) = X( N ) / A( 1 , N )
      JJ = 2
      DO 160 I = 1 , NM
         K = N-I
         TEMP = ( 0.D0 , 0.D0 )
         DO 170 J = 2 , JJ
            TEMP = A( J , K ) * X( K-1+J ) + TEMP
  170    CONTINUE
         X( K ) = X( K ) - TEMP
         IF( JJ .LT. L ) JJ = JJ + 1
  160 CONTINUE
C
      IERR = 0
      RETURN
C
 9000 IERR = 10000
      RETURN
 9002 IERR = 30000
      RETURN
C
      END
C
C     ****** COMPLEX FFT   USING LIST VECTOR ******
C
C                                         CODING BY K.HAMAMATSU
C
      SUBROUTINE FFT2L ( A , B , T , LIST , N2 , IND , KEY , LP )
C
      IMPLICIT COMPLEX*16(A-C,E-H,O-Z)
      REAL * 8     DIV
      DIMENSION    A ( N2*2 ) , B( N2*2 ) , T( N2*LP , 2 ) ,
     &             LIST( N2*LP )
      SAVE DIV
C
      IF(LP.EQ.0) THEN
         B(1) = A(1)
         RETURN
      ENDIF
C     IF ( IND .EQ. 1 ) THEN SET TALBES AND LIST
      IF ( IND .EQ. 1 ) THEN
           CALL SETTBL ( N2 , T , B , LP )
           CALL SETLST ( LIST , N2 , LP )
           DIV = 1.D0 / (N2*2)
           IND = 0
      ENDIF
      K =  1
      L = N2
      DO 10 I = 1 , LP
           M = (I-1)*N2 + 1
           IF ( K .EQ. 1 ) THEN
                CALL FFTSUB( A , B , T(M,KEY) , LIST(M) , L , N2 )
           ELSE
                CALL FFTSUB( B , A , T(M,KEY) , LIST(M) , L , N2 )
           ENDIF
           K = K * (-1)
           L = L / 2
 10   CONTINUE
C     IF ( KEY .EQ. 2 ) THEN INVERSE TRANSFORMATION
C     ELSE    NORMAL TRANSFORMATION
      IF ( KEY .EQ. 1 ) THEN
           IF ( K .EQ. 1 ) THEN
                DO 20 I = 1 , N2*2
 20                 B( I ) = A( I ) * DIV
            ELSE
                DO 30 I = 1 , N2*2
 30                 B( I ) = B( I ) * DIV
            ENDIF
      ELSE IF ( K .EQ. 1 ) THEN
            DO 40 I = 1 , N2*2
 40             B( I ) = A ( I )
      ENDIF
      RETURN
      END
C
C     ====== TRANSFORMATION FOR FFT ======
C                                          CODING BY K.HAMAMATSU
C
      SUBROUTINE FFTSUB ( A , B , T , LIST , L , N2 )
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
      DIMENSION A( N2*2 ) , B( N2 , 2 ) , T( N2 ) , LIST( N2 )
      DO 10 I = 1 , N2
          B( I , 1 ) = A( LIST( I ) ) + A( LIST( I ) + L ) * T( I )
 10   CONTINUE
      DO 20 I = 1 , N2
          B( I , 2 ) = A( LIST( I ) ) - A( LIST( I ) + L ) * T( I )
 20   CONTINUE
      RETURN
      END
C     ====== TABLE SETTING FOR FFT ======
C                                          CODING BY K.HAMAMATSU
C
      SUBROUTINE SETTBL ( N2 , T , B , LP )
C
      COMPLEX * 16     T , B
      REAL    *  8     TR , TI , PAI
      DIMENSION       T( N2 , LP , 2 ) , B( N2 , 2 )
      DATA PAI / 3.1415926535898D0 /
C
      DO 10 I = 1 , N2
          TR = DCOS( 2.D0 * PAI * (I-1) / (N2*2) )
          TI = DSIN( 2.D0 * PAI * (I-1) / (N2*2) )
          B( I , 1 ) = DCMPLX( TR , -TI )
          B( I , 2 ) = DCMPLX( TR ,  TI )
 10   CONTINUE
      K  =  1
      NN = N2
      DO 40 L = 1 , LP
          DO 30 J = 0 , K-1
              DO 20 I = 1 , NN
                  T( I+J*NN , L , 1 ) = B( 1+NN*J , 1 )
                  T( I+J*NN , L , 2 ) = B( 1+NN*J , 2 )
 20           CONTINUE
 30       CONTINUE
          K  = K * 2
          NN = NN / 2
 40   CONTINUE
      RETURN
      END
C
C     ====== LIST SETTING FOR FFT ======
C                                          CODING BY K.HAMAMATSU
C
      SUBROUTINE SETLST ( LIST , N2 , LP )
C
      DIMENSION  LIST ( N2 , LP )
      N1 = N2
      NN =  1
      DO 30 K = 1 , LP
          M = 0
          DO 20 J = 1 , NN
              DO 10 I = 1 , N1
                  M = M + 1
                  LIST ( M , K ) = I + (J-1) * 2 * N1
 10           CONTINUE
 20       CONTINUE
          N1 = N1 / 2
          NN = NN * 2
 30   CONTINUE
      RETURN
      END
C
C     ****** INITIALIZE PLASMA DISPERSION FUNCTION ******
C
      SUBROUTINE DSPFNI
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /COMDSP/
     &           EN10, EN11, EN12, EN13, EN14, EN15, EN16, EN17, EN18,
     &     EN19, EN1A, EN1B, EN1C, EN1D, EN1E, EN21, EN22, EN23, EN24,
     &     EN25, EN26, EN27, EN28, EN29, EN2A, EN2B, EN2C, EN2D, EN2E,
     &     EN31, EN32, EN33, EN34, EN35, EN36, EN37, EN38, EN39, EN3A,
     &     EN3B, EN3C, EN3D, EN3E,  EE1,  EE2,  EE3,  EE4,  EE5,  EE6,
     &      EE7,  EE8,  EE9,  EEA,  EEB,  EEC,  EED,  EEE,
     &           FN10, FN11, FN12, FN13, FN14, FN15, FN16, FN17, FN18,
     &     FN19, FN1A, FN1B, FN1C, FN1D, FN1E, FN21, FN22, FN23, FN24,
     &     FN25, FN26, FN27, FN28, FN29, FN2A, FN2B, FN2C, FN2D, FN2E,
     &     FN31, FN32, FN33, FN34, FN35, FN36, FN37, FN38, FN39, FN3A,
     &     FN3B, FN3C, FN3D, FN3E,  FE1,  FE2,  FE3,  FE4,  FE5,  FE6,
     &      FE7,  FE8,  FE9,  FEA,  FEB,  FEC,  FED,  FEE,
     &    COEF3,COEF2,COEF1,COEF0,H2PAI, HPAI,CPAI1,CPAI2, H,H2, HW,M
      H=0.484375D0
      M=14
      CPAI1= 3.141592653589793D0
      CPAI2= 6.283185307179586D0
      COEF0= 5.641895835477563D-1
      COEF2= 1.128379167095513D0
      COEF1=-COEF2
      COEF3=-2.256758334191025D0
      H2=H*0.5D0
      HW=H+H
      H2PAI=CPAI2/H
      HPAI =CPAI1/H
      HH=H*H
      EN1E=        HH
      EN1D=  4.0D0*HH
      EN1C=  9.0D0*HH
      EN1B= 16.0D0*HH
      EN1A= 25.0D0*HH
      EN19= 36.0D0*HH
      EN18= 49.0D0*HH
      EN17= 64.0D0*HH
      EN16= 81.0D0*HH
      EN15=100.0D0*HH
      EN14=121.0D0*HH
      EN13=144.0D0*HH
      EN12=169.0D0*HH
      EN11=196.0D0*HH
      FN1E=  0.25D0*HH
      FN1D=  2.25D0*HH
      FN1C=  6.25D0*HH
      FN1B= 12.25D0*HH
      FN1A= 20.25D0*HH
      FN19= 30.25D0*HH
      FN18= 42.25D0*HH
      FN17= 56.25D0*HH
      FN16= 72.25D0*HH
      FN15= 90.25D0*HH
      FN14=110.25D0*HH
      FN13=132.25D0*HH
      FN12=156.25D0*HH
      FN11=182.25D0*HH
       EEE=DEXP(-       HH)*HW
       EED=DEXP(- 3.0D0*HH)
       EEC=DEXP(- 5.0D0*HH)
       EEB=DEXP(- 7.0D0*HH)
       EEA=DEXP(- 9.0D0*HH)
       EE9=DEXP(-11.0D0*HH)
       EE8=DEXP(-13.0D0*HH)
       EE7=DEXP(-15.0D0*HH)
       EE6=DEXP(-17.0D0*HH)
       EE5=DEXP(-19.0D0*HH)
       EE4=DEXP(-21.0D0*HH)
       EE3=DEXP(-23.0D0*HH)
       EE2=DEXP(-25.0D0*HH)
       EE1=DEXP(-27.0D0*HH)
       FEE=DEXP(-0.25D0*HH)*HW
       FED=DEXP(- 2.0D0*HH)
       FEC=DEXP(- 4.0D0*HH)
       FEB=DEXP(- 6.0D0*HH)
       FEA=DEXP(- 8.0D0*HH)
       FE9=DEXP(-10.0D0*HH)
       FE8=DEXP(-12.0D0*HH)
       FE7=DEXP(-14.0D0*HH)
       FE6=DEXP(-16.0D0*HH)
       FE5=DEXP(-18.0D0*HH)
       FE4=DEXP(-20.0D0*HH)
       FE3=DEXP(-22.0D0*HH)
       FE2=DEXP(-24.0D0*HH)
       FE1=DEXP(-26.0D0*HH)
      EN2E=EN1E*2.0D0-1.0D0
      EN2D=EN1D*2.0D0-1.0D0
      EN2C=EN1C*2.0D0-1.0D0
      EN2B=EN1B*2.0D0-1.0D0
      EN2A=EN1A*2.0D0-1.0D0
      EN29=EN19*2.0D0-1.0D0
      EN28=EN18*2.0D0-1.0D0
      EN27=EN17*2.0D0-1.0D0
      EN26=EN16*2.0D0-1.0D0
      EN25=EN15*2.0D0-1.0D0
      EN24=EN14*2.0D0-1.0D0
      EN23=EN13*2.0D0-1.0D0
      EN22=EN12*2.0D0-1.0D0
      EN21=EN11*2.0D0-1.0D0
      EN3E=(EN1E*2.0D0-3.0D0)*EN1E
      EN3D=(EN1D*2.0D0-3.0D0)*EN1D
      EN3C=(EN1C*2.0D0-3.0D0)*EN1C
      EN3B=(EN1B*2.0D0-3.0D0)*EN1B
      EN3A=(EN1A*2.0D0-3.0D0)*EN1A
      EN39=(EN19*2.0D0-3.0D0)*EN19
      EN38=(EN18*2.0D0-3.0D0)*EN18
      EN37=(EN17*2.0D0-3.0D0)*EN17
      EN36=(EN16*2.0D0-3.0D0)*EN16
      EN35=(EN15*2.0D0-3.0D0)*EN15
      EN34=(EN14*2.0D0-3.0D0)*EN14
      EN33=(EN13*2.0D0-3.0D0)*EN13
      EN32=(EN12*2.0D0-3.0D0)*EN12
      EN31=(EN11*2.0D0-3.0D0)*EN11
      FN2E=FN1E*2.0D0-1.0D0
      FN2D=FN1D*2.0D0-1.0D0
      FN2C=FN1C*2.0D0-1.0D0
      FN2B=FN1B*2.0D0-1.0D0
      FN2A=FN1A*2.0D0-1.0D0
      FN29=FN19*2.0D0-1.0D0
      FN28=FN18*2.0D0-1.0D0
      FN27=FN17*2.0D0-1.0D0
      FN26=FN16*2.0D0-1.0D0
      FN25=FN15*2.0D0-1.0D0
      FN24=FN14*2.0D0-1.0D0
      FN23=FN13*2.0D0-1.0D0
      FN22=FN12*2.0D0-1.0D0
      FN21=FN11*2.0D0-1.0D0
      FN3E=(FN1E*2.0D0-3.0D0)*FN1E
      FN3D=(FN1D*2.0D0-3.0D0)*FN1D
      FN3C=(FN1C*2.0D0-3.0D0)*FN1C
      FN3B=(FN1B*2.0D0-3.0D0)*FN1B
      FN3A=(FN1A*2.0D0-3.0D0)*FN1A
      FN39=(FN19*2.0D0-3.0D0)*FN19
      FN38=(FN18*2.0D0-3.0D0)*FN18
      FN37=(FN17*2.0D0-3.0D0)*FN17
      FN36=(FN16*2.0D0-3.0D0)*FN16
      FN35=(FN15*2.0D0-3.0D0)*FN15
      FN34=(FN14*2.0D0-3.0D0)*FN14
      FN33=(FN13*2.0D0-3.0D0)*FN13
      FN32=(FN12*2.0D0-3.0D0)*FN12
      FN31=(FN11*2.0D0-3.0D0)*FN11
C
      EN10=EN11*EE1
      EN21=EN21*EE1
      EN31=EN31*EE1
      FN10=FN11*FE1
      FN21=FN21*FE1
      FN31=FN31*FE1
      RETURN
      END
C
C     ****** PLASMA DISPERSION FUNCTION /VECTORIZED VERSION/ ******
C
      SUBROUTINE DSPFNA( Z , DZ , DDZ , NMAX )
C
C   # PLASMA DISPERSION FUNCTION   BY. T.WATANABE   /1978/6/22
C   # REVISED                      BY. T.WATANABE   /1984/4/25
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16              Z ( * ), DZ( * ), DDZ( * )
      COMMON /COMDSP/
     &           EN10, EN11, EN12, EN13, EN14, EN15, EN16, EN17, EN18,
     &     EN19, EN1A, EN1B, EN1C, EN1D, EN1E, EN21, EN22, EN23, EN24,
     &     EN25, EN26, EN27, EN28, EN29, EN2A, EN2B, EN2C, EN2D, EN2E,
     &     EN31, EN32, EN33, EN34, EN35, EN36, EN37, EN38, EN39, EN3A,
     &     EN3B, EN3C, EN3D, EN3E,  EE1,  EE2,  EE3,  EE4,  EE5,  EE6,
     &      EE7,  EE8,  EE9,  EEA,  EEB,  EEC,  EED,  EEE,
     &           FN10, FN11, FN12, FN13, FN14, FN15, FN16, FN17, FN18,
     &     FN19, FN1A, FN1B, FN1C, FN1D, FN1E, FN21, FN22, FN23, FN24,
     &     FN25, FN26, FN27, FN28, FN29, FN2A, FN2B, FN2C, FN2D, FN2E,
     &     FN31, FN32, FN33, FN34, FN35, FN36, FN37, FN38, FN39, FN3A,
     &     FN3B, FN3C, FN3D, FN3E,  FE1,  FE2,  FE3,  FE4,  FE5,  FE6,
     &      FE7,  FE8,  FE9,  FEA,  FEB,  FEC,  FED,  FEE,
     &                 COEF3, COEF2, COEF1, COEF0,
     &                 H2PAI, HPAI, CPAI1, CPAI2, H, H2, HW, M
C
      DATA IDISP/0/
C
      IF(IDISP.EQ.0) THEN
         CALL DSPFNI
         IDISP=1
      ENDIF
      LMAX = NMAX
      DO 200 LOU=1,LMAX
       XRE=DREAL ( Z (LOU ) )
       XIM=DIMAG ( Z (LOU ) )
       XXR=(XRE-XIM)*(XRE+XIM)
       XXI=XRE*XIM*2.0D0
       IF(DABS(XXI) .LT. 1.0D-36)    THEN
         SXI=0.0D0
       ELSE
         SXI=XXI*XXI
       ENDIF
       IF(DABS(XIM) .LT. 1.0D-36)    THEN
         SIM=0.0D0
       ELSE
         SIM=XIM*XIM*2.0D0
       ENDIF
       IF(DABS(XRE) .LT. 1.0D-36)    THEN
         SRE=0.0D0
       ELSE
         SRE=XRE*XRE*2.0D0
       ENDIF
       IF(XIM.LE.HPAI.AND.XXR.LE.161.18D0) THEN
        EXXPHR=DEXP(XIM*H2PAI)
        DR=DCOS(XRE*H2PAI)*2.0D0
        DI=(XIM-HPAI)*XRE*2.0D0
        DEXPV=DEXP(-XXR)
        DSINXV=DSIN(XXI)
        DSINDV=DSIN(DI)
        DCOSXV=DCOS(XXI)
        DCOSDV=DCOS(DI)
       ENDIF
       IF(XXR .LT. -50.65D0 .AND . XIM .LT. 0.0D0) THEN
        Z  ( LOU ) = DCMPLX( 1.D22 , 1.D22 )
        DZ ( LOU ) = DCMPLX( 1.D22 , 1.D22 )
        DDZ( LOU ) = DCMPLX( 1.D22 , 1.D22 )
       ELSE
        XXX=XRE/H2
CF MODIFIED ON 95/11/20
        IF(ABS(XXX).LT.256*256) THEN
           JPM=MOD(IDNINT(XRE/H2),2)
        ELSE
           JPM=0
        ENDIF
CF
        IF(JPM .EQ. 0) THEN
         D1=FN11-XXR
         D2=FN12-XXR
         D3=FN13-XXR
         D4=FN14-XXR
         D5=FN15-XXR
         D6=FN16-XXR
         D7=FN17-XXR
         D8=FN18-XXR
         D9=FN19-XXR
         DA=FN1A-XXR
         DB=FN1B-XXR
         DC=FN1C-XXR
         DD=FN1D-XXR
         DE=FN1E-XXR
         DI1=1.0D0/(D1*D1+SXI)
         DI2=1.0D0/(D2*D2+SXI)
         DI3=1.0D0/(D3*D3+SXI)
         DI4=1.0D0/(D4*D4+SXI)
         DI5=1.0D0/(D5*D5+SXI)
         DI6=1.0D0/(D6*D6+SXI)
         DI7=1.0D0/(D7*D7+SXI)
         DI8=1.0D0/(D8*D8+SXI)
         DI9=1.0D0/(D9*D9+SXI)
         DIA=1.0D0/(DA*DA+SXI)
         DIB=1.0D0/(DB*DB+SXI)
         DIC=1.0D0/(DC*DC+SXI)
         DID=1.0D0/(DD*DD+SXI)
         DIE=1.0D0/(DE*DE+SXI)
         DR1=DI1*D1
         DR2=DI2*D2
         DR3=DI3*D3
         DR4=DI4*D4
         DR5=DI5*D5
         DR6=DI6*D6
         DR7=DI7*D7
         DR8=DI8*D8
         DR9=DI9*D9
         DRA=DIA*DA
         DRB=DIB*DB
         DRC=DIC*DC
         DRD=DID*DD
         DRE=DIE*DE
         SZ0I=(((((((((((((
     &             DI1 *FE1+     DI2)*FE2+     DI3)*FE3+     DI4)*FE4+
     &             DI5)*FE5+     DI6)*FE6+     DI7)*FE7+     DI8)*FE8+
     &             DI9)*FE9+     DIA)*FEA+     DIB)*FEB+     DIC)*FEC+
     &             DID)*FED+     DIE)*FEE
         SZ0R=(((((((((((((
     &             DR1 *FE1+     DR2)*FE2+     DR3)*FE3+     DR4)*FE4+
     &             DR5)*FE5+     DR6)*FE6+     DR7)*FE7+     DR8)*FE8+
     &             DR9)*FE9+     DRA)*FEA+     DRB)*FEB+     DRC)*FEC+
     &             DRD)*FED+     DRE)*FEE
         SZ1I=(((((((((((((
     &        FN10*DI1     +FN12*DI2)*FE2+FN13*DI3)*FE3+FN14*DI4)*FE4+
     &        FN15*DI5)*FE5+FN16*DI6)*FE6+FN17*DI7)*FE7+FN18*DI8)*FE8+
     &        FN19*DI9)*FE9+FN1A*DIA)*FEA+FN1B*DIB)*FEB+FN1C*DIC)*FEC+
     &        FN1D*DID)*FED+FN1E*DIE)*FEE
         SZ1R=(((((((((((((
     &        FN10*DR1     +FN12*DR2)*FE2+FN13*DR3)*FE3+FN14*DR4)*FE4+
     &        FN15*DR5)*FE5+FN16*DR6)*FE6+FN17*DR7)*FE7+FN18*DR8)*FE8+
     &        FN19*DR9)*FE9+FN1A*DRA)*FEA+FN1B*DRB)*FEB+FN1C*DRC)*FEC+
     &        FN1D*DRD)*FED+FN1E*DRE)*FEE
         SZ2I=(((((((((((((
     &        FN21*DI1     +FN22*DI2)*FE2+FN23*DI3)*FE3+FN24*DI4)*FE4+
     &        FN25*DI5)*FE5+FN26*DI6)*FE6+FN27*DI7)*FE7+FN28*DI8)*FE8+
     &        FN29*DI9)*FE9+FN2A*DIA)*FEA+FN2B*DIB)*FEB+FN2C*DIC)*FEC+
     &        FN2D*DID)*FED+FN2E*DIE)*FEE
         SZ2R=(((((((((((((
     &        FN21*DR1     +FN22*DR2)*FE2+FN23*DR3)*FE3+FN24*DR4)*FE4+
     &        FN25*DR5)*FE5+FN26*DR6)*FE6+FN27*DR7)*FE7+FN28*DR8)*FE8+
     &        FN29*DR9)*FE9+FN2A*DRA)*FEA+FN2B*DRB)*FEB+FN2C*DRC)*FEC+
     &        FN2D*DRD)*FED+FN2E*DRE)*FEE
CF       SZ3I=(((((((((((((
CF   &        FN31*DI1     +FN32*DI2)*FE2+FN33*DI3)*FE3+FN34*DI4)*FE4+
CF   &        FN35*DI5)*FE5+FN36*DI6)*FE6+FN37*DI7)*FE7+FN38*DI8)*FE8+
CF   &        FN39*DI9)*FE9+FN3A*DIA)*FEA+FN3B*DIB)*FEB+FN3C*DIC)*FEC+
CF   &        FN3D*DID)*FED+FN3E*DIE)*FEE
CF       SZ3R=(((((((((((((
CF   &        FN31*DR1     +FN32*DR2)*FE2+FN33*DR3)*FE3+FN34*DR4)*FE4+
CF   &        FN35*DR5)*FE5+FN36*DR6)*FE6+FN37*DR7)*FE7+FN38*DR8)*FE8+
CF   &        FN39*DR9)*FE9+FN3A*DRA)*FEA+FN3B*DRB)*FEB+FN3C*DRC)*FEC+
CF   &        FN3D*DRD)*FED+FN3E*DRE)*FEE
         TZ0I=(SZ0R+SRE*SZ0I)*XIM
         TZ0R=(SZ0R-SIM*SZ0I)*XRE
         TZ2I=(SZ2R+SRE*SZ2I)*XIM
         TZ2R=(SZ2R-SIM*SZ2I)*XRE
        ELSE
         D1=EN11-XXR
         D2=EN12-XXR
         D3=EN13-XXR
         D4=EN14-XXR
         D5=EN15-XXR
         D6=EN16-XXR
         D7=EN17-XXR
         D8=EN18-XXR
         D9=EN19-XXR
         DA=EN1A-XXR
         DB=EN1B-XXR
         DC=EN1C-XXR
         DD=EN1D-XXR
         DE=EN1E-XXR
         DI1=1.0D0/(D1*D1+SXI)
         DI2=1.0D0/(D2*D2+SXI)
         DI3=1.0D0/(D3*D3+SXI)
         DI4=1.0D0/(D4*D4+SXI)
         DI5=1.0D0/(D5*D5+SXI)
         DI6=1.0D0/(D6*D6+SXI)
         DI7=1.0D0/(D7*D7+SXI)
         DI8=1.0D0/(D8*D8+SXI)
         DI9=1.0D0/(D9*D9+SXI)
         DIA=1.0D0/(DA*DA+SXI)
         DIB=1.0D0/(DB*DB+SXI)
         DIC=1.0D0/(DC*DC+SXI)
         DID=1.0D0/(DD*DD+SXI)
         DIE=1.0D0/(DE*DE+SXI)
         DR1=DI1*D1
         DR2=DI2*D2
         DR3=DI3*D3
         DR4=DI4*D4
         DR5=DI5*D5
         DR6=DI6*D6
         DR7=DI7*D7
         DR8=DI8*D8
         DR9=DI9*D9
         DRA=DIA*DA
         DRB=DIB*DB
         DRC=DIC*DC
         DRD=DID*DD
         DRE=DIE*DE
         SZ0I=(((((((((((((
     &             DI1 *EE1+     DI2)*EE2+     DI3)*EE3+     DI4)*EE4+
     &             DI5)*EE5+     DI6)*EE6+     DI7)*EE7+     DI8)*EE8+
     &             DI9)*EE9+     DIA)*EEA+     DIB)*EEB+     DIC)*EEC+
     &             DID)*EED+     DIE)*EEE
         SZ0R=(((((((((((((
     &             DR1 *EE1+     DR2)*EE2+     DR3)*EE3+     DR4)*EE4+
     &             DR5)*EE5+     DR6)*EE6+     DR7)*EE7+     DR8)*EE8+
     &             DR9)*EE9+     DRA)*EEA+     DRB)*EEB+     DRC)*EEC+
     &             DRD)*EED+     DRE)*EEE
         SZ1I=(((((((((((((
     &        EN10*DI1     +EN12*DI2)*EE2+EN13*DI3)*EE3+EN14*DI4)*EE4+
     &        EN15*DI5)*EE5+EN16*DI6)*EE6+EN17*DI7)*EE7+EN18*DI8)*EE8+
     &        EN19*DI9)*EE9+EN1A*DIA)*EEA+EN1B*DIB)*EEB+EN1C*DIC)*EEC+
     &        EN1D*DID)*EED+EN1E*DIE)*EEE
         SZ1R=(((((((((((((
     &        EN10*DR1     +EN12*DR2)*EE2+EN13*DR3)*EE3+EN14*DR4)*EE4+
     &        EN15*DR5)*EE5+EN16*DR6)*EE6+EN17*DR7)*EE7+EN18*DR8)*EE8+
     &        EN19*DR9)*EE9+EN1A*DRA)*EEA+EN1B*DRB)*EEB+EN1C*DRC)*EEC+
     &        EN1D*DRD)*EED+EN1E*DRE)*EEE
         SZ2I=(((((((((((((
     &        EN21*DI1     +EN22*DI2)*EE2+EN23*DI3)*EE3+EN24*DI4)*EE4+
     &        EN25*DI5)*EE5+EN26*DI6)*EE6+EN27*DI7)*EE7+EN28*DI8)*EE8+
     &        EN29*DI9)*EE9+EN2A*DIA)*EEA+EN2B*DIB)*EEB+EN2C*DIC)*EEC+
     &        EN2D*DID)*EED+EN2E*DIE)*EEE
         SZ2R=(((((((((((((
     &        EN21*DR1     +EN22*DR2)*EE2+EN23*DR3)*EE3+EN24*DR4)*EE4+
     &        EN25*DR5)*EE5+EN26*DR6)*EE6+EN27*DR7)*EE7+EN28*DR8)*EE8+
     &        EN29*DR9)*EE9+EN2A*DRA)*EEA+EN2B*DRB)*EEB+EN2C*DRC)*EEC+
     &        EN2D*DRD)*EED+EN2E*DRE)*EEE
CF       SZ3I=(((((((((((((
CF   &        EN31*DI1     +EN32*DI2)*EE2+EN33*DI3)*EE3+EN34*DI4)*EE4+
CF   &        EN35*DI5)*EE5+EN36*DI6)*EE6+EN37*DI7)*EE7+EN38*DI8)*EE8+
CF   &        EN39*DI9)*EE9+EN3A*DIA)*EEA+EN3B*DIB)*EEB+EN3C*DIC)*EEC+
CF   &        EN3D*DID)*EED+EN3E*DIE)*EEE
CF       SZ3R=(((((((((((((
CF   &        EN31*DR1     +EN32*DR2)*EE2+EN33*DR3)*EE3+EN34*DR4)*EE4+
CF   &        EN35*DR5)*EE5+EN36*DR6)*EE6+EN37*DR7)*EE7+EN38*DR8)*EE8+
CF   &        EN39*DR9)*EE9+EN3A*DRA)*EEA+EN3B*DRB)*EEB+EN3C*DRC)*EEC+
CF   &        EN3D*DRD)*EED+EN3E*DRE)*EEE
         DD=HW/(SRE+SIM)
         TZ0I=(SZ0R+SRE*SZ0I+DD)*XIM
         TZ0R=(SZ0R-SIM*SZ0I-DD)*XRE
         TZ2I=(SZ2R+SRE*SZ2I-DD)*XIM
         TZ2R=(SZ2R-SIM*SZ2I+DD)*XRE
        ENDIF
        SZ1I=SZ1I*XXI
CF      SZ3I=SZ3I*XXI
        IF(XIM .GT. HPAI .OR. XXR .GT. 161.18D0) THEN
         Z   ( LOU ) = DCMPLX( TZ0R , TZ0I ) * COEF0
         DZ  ( LOU ) = DCMPLX( SZ1R , SZ1I ) * COEF1
         DDZ ( LOU ) = DCMPLX( TZ2R , TZ2I ) * COEF2
CF       DDDZ( LOU ) = DCMPLX( SZ3R , SZ3I ) * COEF3
        ELSE
CSX         EXXPHR=DEXP(XIM*H2PAI)
CSX         DR=DCOS(XRE*H2PAI)*2.0D0
CSX         DI=(XIM-HPAI)*XRE*2.0D0
         IF(JPM .EQ. 0) THEN
          EXXR=DEXPV*CPAI2/((EXXPHR+DR)*EXXPHR+1.0D0)
          XN0R=(DSINXV+DSINDV*EXXPHR)*EXXR
          XN0I=(DCOSXV+DCOSDV*EXXPHR)*EXXR
         ELSE
          EXXR=DEXPV*CPAI2/((EXXPHR-DR)*EXXPHR+1.0D0)
          XN0R=(DSINXV-DSINDV*EXXPHR)*EXXR
          XN0I=(DCOSXV-DCOSDV*EXXPHR)*EXXR
         ENDIF
         XN1R=XRE*XN0R-XIM*XN0I
         XN1I=XRE*XN0I+XIM*XN0R
         X2R =XXR*2.0D0-1.0D0
         X2I =XXI*2.0D0
         XN2R=X2R*XN0R-X2I*XN0I
         XN2I=X2R*XN0I+X2I*XN0R
CF       X3R =(SRE-SIM*3.0D0-3.0D0)*XRE
CF       X3I =(SRE*3.0D0-SIM-3.0D0)*XIM
CF       XN3R=X3R*XN0R-X3I*XN0I
CF       XN3I=X3R*XN0I+X3I*XN0R
         Z   ( LOU ) = DCMPLX( TZ0R+XN0R , TZ0I+XN0I ) * COEF0
         DZ  ( LOU ) = DCMPLX( SZ1R+XN1R , SZ1I+XN1I ) * COEF1
         DDZ ( LOU ) = DCMPLX( TZ2R+XN2R , TZ2I+XN2I ) * COEF2
CF       DDDZ( LOU ) = DCMPLX( SZ3R+XN3R , SZ3I+XN3I ) * COEF3
        ENDIF
       ENDIF
 200   CONTINUE
      RETURN
      END
C
C     ****** LAMBDA FUNCTION ******
C
      SUBROUTINE LAMBDA(N,X,ALAM)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      DIMENSION ALAM(0:N)
      DATA ONE/1.D0/
      DATA D55/1.D-55/
C
      XA=ABS(X)
      NA=ABS(N)
C
      IF(NA.GE.30000.OR.XA.GE.173.D0) THEN
         WRITE(6,*) 'XX LAMBDA : OUT OF RANGE : N,X = ',N,X
         RETURN
      ENDIF
      IF(XA.LE.1.D-8) THEN
         IF(XA.LE.1.D-77) THEN
            ALAM(0)=ONE
            DO 10 I=1,NA
               ALAM(I)=0.D0
   10       CONTINUE
         ELSE
            ALAM(0)=EXP(-XA)
            T1=0.5D0*XA
            T2=ONE
            T3=ONE
            DO 20 I=1,NA
               IF(T3.LE.1.D-77*T2/T1) THEN
                  ALAM(I)=0.D0
               ELSE
                  T3=T3*T1/T2
                  T2=T2+ONE
                  ALAM(I)=T3*EXP(-XA)
               ENDIF
   20       CONTINUE
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
         DO 100 I=1,NM
            K=NM-I
            T1=(K+1)*T2*Z+T3
            IF(K.LE.NA) ALAM(K)=T1
            S=S+T1
            IF(ABS(S).GT.1.D55) THEN
               T1=T1*D55
               T2=T2*D55
               S=S*D55
               DO 50 J=K,NA
                  ALAM(J)=ALAM(J)*D55
   50          CONTINUE
            ENDIF
            T3=T2
            T2=T1
  100    CONTINUE
         S=S+S-T1
         DO 150 J=0,NA
            ALAM(J)=ALAM(J)/S
  150    CONTINUE
      ENDIF
      RETURN
      END

