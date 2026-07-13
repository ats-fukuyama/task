! libpdkf.f90

MODULE decomm
  USE wimcomm, ONLY: rkind
  REAL(rkind):: G1,G2,G3
  INTEGER:: N1
END MODULE decomm

MODULE libpdkf

  PRIVATE
  PUBLIC setup_libpdkf
  PUBLIC DSPLC
  PUBLIC DSPLF

  CONTAINS

!     *****  CALCULATE KERNEL FUNCTION  *****

  SUBROUTINE setup_libpdkf

    USE wimcomm,ONLY: rkind,TT,TMAX,NTMAX,CFK1,CFKS1,CFK2,CFKS2,PI,idebug_wim
    USE decomm,ONLY: N1
    USE libgrf
    IMPLICIT NONE
    INTEGER:: IOPT(2),NT,L,IER
    REAL(rkind):: DT
    COMPLEX(rkind):: CDF(2)
    REAL(rkind):: F(NTMAX,4)

    DT=TMAX/(NTMAX-1) 
    DO NT=1,NTMAX 
       TT(NT)=(NT-1)*DT 
    END DO

    N1=1
    CFK1(1)=(0.D0,1.D0)/SQRT(2.D0*PI)
    DO NT=2,NTMAX
       CALL EUL(TT(NT),CFK1(NT),5,L)
    END DO
    CDF(1)=(-0.5D0,0.D0) 
    IOPT(1)=2
    IOPT(2)=3
    CALL DSPLC(TT,NTMAX,CFK1,CDF,IOPT,CFKS1,NTMAX,IER)

    N1=2
    CFK2(1)=-1.D0/SQRT(2.D0*PI)
    DO NT=2,NTMAX
       CALL EUL(TT(NT),CFK2(NT),5,L)
    END DO
    CDF(1)=(0.D0,0.D0) 
    IOPT(1)=2
    IOPT(2)=3
    CALL DSPLC(TT,NTMAX,CFK2,CDF,IOPT,CFKS2,NTMAX,IER)

    IF(idebug_wim(1).NE.0) THEN
       CALL PAGES
       F(1:NTMAX,1)=REAL(CFK1(1:NTMAX))
       F(1:NTMAX,2)=IMAG(CFK1(1:NTMAX))
       F(1:NTMAX,3)=REAL(CFK2(1:NTMAX))
       F(1:NTMAX,4)=IMAG(CFK2(1:NTMAX))
       CALL GRD1D(0,TT,F,NTMAX,NTMAX,4,'CFK')
       CALL PAGEE
    END IF

    RETURN
  END SUBROUTINE setup_libpdkf

!     *****  EULER TRANSFOMATION  *****

  SUBROUTINE EUL(Z,CS,M,L)

    USE wimcomm,ONLY: rkind
    USE libde
    USE decomm,ONLY: G1,G2,G3,N1
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: Z
    COMPLEX(rkind),INTENT(OUT):: CS
    INTEGER,INTENT(IN):: M
    INTEGER,INTENT(OUT):: L
    REAL(rkind):: A(200),B(200)
    REAL(rkind),PARAMETER:: HP=1.5707963267948966192D0
    INTEGER:: ILST,K
    REAL(rkind):: H0,EPS,SR1,SI1,SR,SI,ESR,ESI,SR2,SI2,PARITY,SKR,SKI

    G2=HP
    G3=Z
    H0=0.5D0
    EPS=1.D-5
    ILST=0
     
    SR1=0.D0
    SI1=0.D0
    IF(M.NE.0) THEN 
       DO K=M-1,0,-1
          G1=DBLE(K)
          CALL DEFT(SR,ESR,H0,EPS,ILST,FUNR,'FUNR')
          CALL DEFT(SI,ESI,H0,EPS,ILST,FUNI,'FUNI')
          SR1=SR1+SR
          SI1=SI1+SI
       END DO
    ENDIF

    G1=DBLE(M)
    CALL DEFT(SR,ESR,H0,EPS,ILST,FUNR,'FUNR')
    CALL DEFT(SI,ESI,H0,EPS,ILST,FUNI,'FUNI')
    A(1)=SR
    SR2=0.5D0*SR
    B(1)=SI
    SI2=0.5D0*SI
    PARITY=-1.D0
    L=0

30  CONTINUE
    L=L+1
    IF(L.GE.200) GOTO 9000
    G1=DBLE(M+L)
    CALL DEFT(SR,ESR,H0,EPS,ILST,FUNR,'FUNR')
    CALL DEFT(SI,ESI,H0,EPS,ILST,FUNI,'FUNI')
    A(L+1)=SR*PARITY
    B(L+1)=SI*PARITY
    DO K=L,1,-1    
       A(K)=A(K+1)-A(K)
       B(K)=B(K+1)-B(K)
    END DO
    SKR=A(1)*PARITY*0.5D0**(L+1)
    SR2=SR2+SKR
    SKI=B(1)*PARITY*0.5D0**(L+1)
    SI2=SI2+SKI
    PARITY=-PARITY
    IF(DABS(SKR).GT.EPS.OR.DABS(SKI).GT.EPS) GOTO 30

    SR=(SR1+SR2)/DSQRT(4.D0*G2)
    SI=(SI1+SI2)/DSQRT(4.D0*G2)
    CS=CMPLX(SR,SI)
    RETURN

9000 WRITE(6,*) '## DIMENSION OVERFLOW IN EULER TRANSFORMATION.'
    RETURN
  END SUBROUTINE EUL

!     *****  REAL PART  *****

  FUNCTION FUNR(X,XM,XP)

    USE wimcomm,ONLY: rkind
    USE decomm,ONLY: G1,G2,G3,N1
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,XM,XP
    REAL(rkind):: FUNR,Y1,T,T2,YY,AN2

    Y1=XM
    IF(IDINT(G1).EQ.0) THEN 
       T=0.5D0*G2*XP
       T2=T*T
       YY=-0.5D0*G3*G3/T2
       IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
          IF(N1.EQ.1) THEN
             AN2=1.D0
          ELSEIF(N1.EQ.2) THEN
             AN2=T
          ELSEIF(N1.EQ.3) THEN
             AN2=T2
          ELSEIF(N1.EQ.4) THEN
             AN2=T2*T
          ELSEIF(N1.EQ.5) THEN
             AN2=G3*G3/(T2*T) 
          ENDIF
          FUNR=AN2*0.5*G2*DEXP(YY)*DCOS(T)
       ELSE
          FUNR=0.D0
       ENDIF
    ELSE
       T=G2*(X+2.D0*G1)
       T2=T*T
       YY=-0.5D0*G3*G3/T2
       IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
          IF(N1.EQ.1) THEN
             AN2=1.D0
          ELSEIF(N1.EQ.2) THEN
             AN2=T
          ELSEIF(N1.EQ.3) THEN
             AN2=T2
          ELSEIF(N1.EQ.4) THEN
             AN2=T2*T
          ELSEIF(N1.EQ.5) THEN
             AN2=G3*G3/(T2*T) 
          ENDIF
          FUNR=AN2*G2*DEXP(YY)*DCOS(T)
       ELSE
          FUNR=0.D0
       ENDIF
    ENDIF
    RETURN
  END FUNCTION FUNR

!     *****  IMAG PART  *****

  FUNCTION FUNI(X,XM,XP)

    USE wimcomm,ONLY: rkind
    USE decomm,ONLY: G1,G2,G3,N1
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X,XM,XP
    REAL(rkind):: FUNI,Y1,Y2,T,T2,YY,AN2

    Y1=X
    Y2=XM
    T=G2*(XP+2.D0*G1)
    T2=T*T
    YY=-0.5D0*G3*G3/T2
    IF(ABS(YY).LT.352.D0.AND.ABS(T).LT.1.D4) THEN
       IF(N1.EQ.1) THEN
          AN2=1.D0
       ELSEIF(N1.EQ.2) THEN 
          AN2=T
       ELSEIF(N1.EQ.3) THEN
          AN2=T2
       ELSEIF(N1.EQ.4) THEN
          AN2=T2*T
       ELSEIF(N1.EQ.5) THEN
          AN2=G3*G3/(T2*T) 
       ENDIF
       FUNI=AN2*G2*DEXP(YY)*DSIN(T)
    ELSE
       FUNI=0.D0
    ENDIF
    RETURN
  END FUNCTION FUNI


!     *****  CALCULATE CUBIC SPLINE COEFFICIENT  ***** 

      SUBROUTINE DSPLC(X,N,Y,DF,IOPT,C,NC,IER)
      IMPLICIT REAL*8(A-H,O-Z)
 
!     COMPUTE THE COEFICIENTS OF THE CUBIC SPLINE.

      COMPLEX*16 Y(N),DF(2),C(NC,3),EC(4),D(2)
      COMPLEX*16 HY,Y1,Y2,PIV,DY1,DY2
      DIMENSION X(N),IOPT(2)
      IF(N.LT.2)    GOTO 9000
      IF(NC.LT.N-1) GOTO 9000
      IF(IOPT(1).LT.1 .OR. IOPT(1).GT.3) GOTO 9000
      IF(IOPT(2).LT.1 .OR. IOPT(2).GT.3) GOTO 9000
      NM1=N-1
      DO 5 I=1,NM1
         IF(X(I).GE.X(I+1)) GOTO 9100
    5 CONTINUE
      IER=0

!     SET THE END CONDITIONS.

      II=2
      KS=1
      KE=MIN0(4,N)
      IDER=1
      DO 70 I=1,2
         I1=2*I-1
         I2=2*I
         IB=IOPT(I)
         GOTO (10,20,30),IB
   10    EC(I1)=0.D0
         EC(I2)=2.D0*DF(I)
         GOTO 70
   20    D(I)=DF(I)
   25    IF(I.EQ.2) II=N
         H=X(II)-X(II-1)
         EC(I1)=1.D0
         HY=Y(II)-Y(II-1)
         EC(I2)=6.D0*(HY/H-D(I))/H
         IF(I.EQ.2) EC(I2)=-EC(I2)
         GOTO 70
   30    IF(I.EQ.1) GOTO 40
         KS=MAX0(1,N-3)
         KE=N
         IDER=N
   40    A2=0.D0
         D(I)=0.D0
         DO 60 K=KS,KE
            IF(IDER.EQ.K) GOTO 60
            A1=1.D0
            DO 50 J=KS,KE
               IF(J.EQ.IDER .OR. J.EQ.K) GOTO 50
               X1=X(IDER)-X(J)
               X2=X(K)-X(J)
               A1=A1*X1/X2
   50       CONTINUE
            X3=X(K)-X(IDER)
            D(I)=D(I)+A1*Y(K)/X3
            A2=A2-1.D0/X3
   60    CONTINUE
         D(I)=D(I)+Y(IDER)*A2
         GOTO 25
   70 CONTINUE

!     SET THE ELEMENTS FOR THE SYMMETRIC TRIDIAGONAL EQUATION. 

      IF(N.EQ.2) GOTO 90
      H1=X(2)-X(1)
      Y1=Y(2)-Y(1)
      DO 80 I=2,NM1
         H2=X(I+1)-X(I)
         Y2=Y(I+1)-Y(I)
         HH=H1+H2
         C(I,1)=H2/HH
         C(I,2)=1.D0-C(I,1)
         C(I,3)=6.D0*(Y2/H2-Y1/H1)/HH
         H1=H2
         Y1=Y2
   80 CONTINUE

!     SOLVE THE EQUATION.

   90 C(1,1)=-EC(1)*0.5
      C(1,2)= EC(2)*0.5
      IF(N.EQ.2) GOTO 110
      DO 100 K=2,NM1
         PIV=2.D0+C(K,2)*C(K-1,1)
         C(K,1)=-C(K,1)/PIV
         C(K,2)=(C(K,3)-C(K,2)*C(K-1,2))/PIV
  100 CONTINUE
  110 DY1=(EC(4)-EC(3)*C(NM1,2))/(2.D0+EC(3)*C(NM1,1))
      DO 120 I=1,NM1
         K=N-I
         DY2=C(K,1)*DY1+C(K,2)
         H=X(K+1)-X(K)
         C(K,3)=(DY1-DY2)/(6.D0*H)
         C(K,2)=0.5D0*DY2
         C(K,1)=(Y(K+1)-Y(K))/H-(C(K,2)+C(K,3)*H)*H
         DY1=DY2
  120 CONTINUE
         GOTO 9999

!         ERROR

 9000 IER=2
      WRITE(6,1000) N,NC,IOPT(1),IOPT(2)
 1000 FORMAT(1H0,'(SUBR.DSPLC)  N =',I4,', NC =',I4, &
                 ', IOPT(1) =',I3,', IOPT(2) =',I3 &
            /' ','N, NC, IOPT(1) ANDIOPT(2) SHOULD SATISFY' &
           /' ','THE FOLLOWING INEQUALITIES.' &
           /' ' ,'    2<=N, N-1<=NC, 1<=IOPT(1)<=3, 1<=IOPT(2)<=3 ' &
           /' ' ,'RETURN WITH NO CALCULATION.' )
      GOTO 9999
 9100 IER=1
      IP1=I+1
      WRITE(6,1100) I,X(I),IP1,X(IP1)
 1100 FORMAT(1H0,'(SUBR.DSPLC)  X(',I4,')=',F10.5, &
                             ', X(',I4,')=',F10.5 &
            /' ','X SHOULD SATISFY THE FOLLOWING INEQUALITIES.' &
            /' ','    X(1)<X(2)<...<X(N)' )
 9999 RETURN
   END SUBROUTINE DSPLC

!     *****  INTERPOLATE BY CUBIC SPLINE  ***** 

      SUBROUTINE DSPLF(X,N,Y,C,NC,V,M,F,IER)
      IMPLICIT REAL*8(A-H,O-Z)

!     INTERPOLATION BY THE CUBIC SPLINE.

      COMPLEX*16 Y(N),C(NC,3),F(M)
      DIMENSION X(N),V(M)
      IF(N.LT.2)    GOTO 9000
      IF(M.LT.1)    GOTO 9000
      IF(NC.LT.N-1) GOTO 9000
      IER=0
      I=1
      DO 90 K=1,M
         V1=V(K)-X(I)
         IF(V1) 10,30,40
   10    IF(I.GT.1) GOTO 20
         IER=1
         GOTO 80
   20    I=I-1
         V1=V(K)-X(I)
         IF(V1) 10,30,80
   30    F(K)=Y(I)
         GOTO 90
   40    IF(I.LT.N) GOTO 50
         IER=1
         I=N-1
         GOTO 80
   50    V2=V(K)-X(I+1)
         IF(V2) 80,60,70
   60    I=I+1
         GOTO 30
   70    I=I+1
         V1=V2
         GOTO 40
   80    F(K)=Y(I)+V1*(C(I,1)+V1*(C(I,2)+V1*C(I,3)))
   90 CONTINUE
      GOTO 9999
 9000 IER=2
      WRITE(6,1000) N,NC,M
 1000 FORMAT(' ','(SUBR.DSPLF)  N =',I4,', NC =',I4,', M =',I4/ &
             ' ','N, NC AND M SHOULD SATISFY THE FOLLOWING INEQUALITIES.'/ &
             ' ','    2<=N, N-1<=NC, 1<=M' &
             ' ','RETURN WITH NO CALCULATION.')
 9999 RETURN
    END SUBROUTINE DSPLF

END MODULE libpdkf
