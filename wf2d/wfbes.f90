!     $Id$
!
!     **************************************
!
!        BESSEL FUNCTION OF THE FIRST KIND
!
!     **************************************

!      SUBROUTINE BESSJN(X,NMAX,BJN,DBJN)

!      INTEGER N,IACC
!      REAL*8  BJN(0:NMAX),DBJN(0:NMAX),X,BIGNO,BIGNI
!      PARAMETER (IACC=40,BIGNO=1.D10,BIGNI=1.D-10)
!U    USES BESSJ0,BESSJ1
!      INTEGER J,JSUM,M
!      REAL*8 AX,BJ,BJM,BJP,SUM,TOX,BESSJ0,BESSJ1

!      BJN(0)=BESSJ0(X)
!      BJN(1)=BESSJ1(X)

!      AX=ABS(X)
!      IF(AX.EQ.0.D0)THEN
!         DO N=2,NMAX
!            BJN(N)=0.D0
!         ENDDO
!      ELSE IF(AX.GT.DBLE(NMAX))THEN
!         TOX=2.D0/AX
!         BJM=BJN(0)
!         BJ=BJN(1)
!         DO J=1,NMAX-1
!            BJP=J*TOX*BJ-BJM
!            BJM=BJ
!            BJ=BJP
!            BJN(J+1)=BJ
!         ENDDO
!      ELSE
!         TOX=2.D0/AX
!         M=2*((NMAX+INT(SQRT(DBLE(IACC*NMAX))))/2)
!         DO N=2,NMAX
!            BJN(N)=0.D0
!         ENDDO
!         JSUM=0
!         SUM=0.D0
!         BJP=0.D0
!         BJ=1.D0
!         DO J=M,1,-1
!            BJM=J*TOX*BJ-BJP
!            BJP=BJ
!            BJ=BJM
!            IF(ABS(BJ).GT.BIGNO)THEN
!               BJ=BJ*BIGNI
!               BJP=BJP*BIGNI
!               DO N=J+1,NMAX
!                  BJN(N)=BJN(N)*BIGNI
!               ENDDO
!               SUM=SUM*BIGNI
!            ENDIF
!            IF(JSUM.NE.0) SUM=SUM+BJ
!            JSUM=1-JSUM
!            IF(J.GE.2.AND.J.LE.NMAX) BJN(J)=BJP
!         ENDDO
!         SUM=2.D0*SUM-BJ
!         DO N=2,NMAX
!            BJN(N)=BJN(N)/SUM
!         ENDDO
!      ENDIF
!      DO N=3,NMAX,2
!         IF(X.LT.0.D0) BJN(N)=-BJN(N)
!      ENDDO
!      DBJN(0)=BJN(1)
!      IF(AX.EQ.0.D0)THEN
!         DBJN(1)=0.5D0
!         DO N=2,NMAX
!            DBJN(N)=0.D0
!         ENDDO
!      ELSE
!         DO N=1,NMAX
!            DBJN(N)=BJN(N-1)-N*BJN(N)/X
!         ENDDO
!      ENDIF
!      RETURN
!      END

!     ************************************************
!
!        BESSEL FUNCTION OF THE FIRST KIND, ORDER N
!
!     ************************************************

FUNCTION BESSJ(N,X)

  implicit none
  integer :: N
  integer,parameter :: IACC=40
  real(8) :: BESSJ,X
  real(8),parameter :: BIGNO=1.D10
  real(8),parameter :: BIGNI=1.D-10
!U    USES BESSJ0,BESSJ1
  integer :: J,JSUM,M
  REAL(8) :: AX,BJ,BJM,BJP,SUM,TOX,BESSJ0,BESSJ1
  
  IF(N.LT.2) THEN
     IF(N.EQ.0) THEN
        BESSJ=BESSJ0(X)
     ELSEIF(N.EQ.1) THEN
        BESSJ=BESSJ1(X)
     ENDIF
     RETURN
  ENDIF
  
  AX=ABS(X)
  IF(AX.EQ.0.D0)THEN
     BESSJ=0.D0
  ELSE IF(AX.GT.DFLOAT(N))THEN
     TOX=2.D0/AX
     BJM=BESSJ0(AX)
     BJ=BESSJ1(AX)
     DO J=1,N-1
        BJP=J*TOX*BJ-BJM
        BJM=BJ
        BJ=BJP
     ENDDO
     BESSJ=BJ
  ELSE
     TOX=2.D0/AX
     M=2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
     BESSJ=0.D0
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
           BESSJ=BESSJ*BIGNI
           SUM=SUM*BIGNI
        ENDIF
        IF(JSUM.NE.0)SUM=SUM+BJ
        JSUM=1-JSUM
        IF(J.EQ.N)BESSJ=BJP
     ENDDO
     SUM=2.D0*SUM-BJ
     BESSJ=BESSJ/SUM
  ENDIF
  IF(X.LT.0.D0.AND.MOD(N,2).EQ.1) BESSJ=-BESSJ
  RETURN
END FUNCTION BESSJ

!     ************************************************
!
!        BESSEL FUNCTION OF THE FIRST KIND, ORDER 0
!
!     ************************************************

FUNCTION BESSJ0(X)

  implicit none
  real(8) :: BESSJ0,X
  real(8) :: AX,XX,Z
  real(8) :: P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6
  real(8) :: S1,S2,S3,S4,S5,S6,Y
  SAVE   P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,&
       &       S1,S2,S3,S4,S5,S6
  DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,&
       &-.2073370639D-5,.2093887211D-6/
  DATA Q1,Q2,Q3,Q4,Q5/-.1562499995D-1,&
       &.1430488765D-3,-.6911147651D-5,.7621095161D-6,-.934945152D-7/
  DATA R1,R2,R3,R4,R5,R6/57568490574.D0,-13362590354.D0,&
       &651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0/
  DATA S1,S2,S3,S4,S5,S6/57568490411.D0,1029532985.D0,9494680.718D0,&
       &59272.64853D0,267.8532712D0,1.D0/
  
  IF(ABS(X).LT.8.D0)THEN
     Y=X**2
     BESSJ0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))&
          &         /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
  ELSE
     AX=ABS(X)
     Z=8.D0/AX
     Y=Z**2
     XX=AX-.785398164D0
     BESSJ0=SQRT(.636619772D0/AX)&
          &         *(  COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*P5))))&
          &          -Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
  ENDIF
  RETURN
END FUNCTION BESSJ0
    
!     ************************************************
!
!        BESSEL FUNCTION OF THE FIRST KIND, ORDER 1
!
!     ************************************************

FUNCTION BESSJ1(X)

  implicit none
  real(8) :: BESSJ1,X
  real(8) :: AX,XX,Z
  real(8) :: P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6
  real(8) :: S1,S2,S3,S4,S5,S6,Y
  SAVE   P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,&
       &       S1,S2,S3,S4,S5,S6
  DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,&
       & 242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0/
  DATA S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0,&
       & 18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
  DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,&
       & .2457520174D-5,-.240337019D-6/
  DATA Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3,&
       &     .8449199096D-5,-.88228987D-6,.105787412D-6/
  
  IF(ABS(X).LT.8.D0)THEN
     Y=X**2
     BESSJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))&
          &          /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
  ELSE
     AX=ABS(X)
     Z=8.D0/AX
     Y=Z**2
     XX=AX-2.356194491D0
     BESSJ1=SQRT(.636619772D0/AX)&
          &          *(  COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*P5))))&
          &           -Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))&
          &          *SIGN(1.D0,X)
  ENDIF
  RETURN
END FUNCTION BESSJ1

!     ************************************************
!
!        BESSEL FUNCTION OF THE SECOND KIND, ORDER N
!
!     ************************************************

FUNCTION BESSY(N,X)

  implicit none
  integer :: N
  real(8) :: BESSY,X
!U    USES BESSY0,BESSY1
  integer :: J
  real(8) :: BY,BYM,BYP,TOX,BESSY0,BESSY1

  IF(N.LT.2) THEN
     IF(N.EQ.0) THEN
        BESSY=BESSY0(X)
     ELSEIF(N.EQ.1) THEN
        BESSY=BESSY1(X)
     ENDIF
     RETURN
  ENDIF
  
  TOX=2.D0/X
  BY=BESSY1(X)
  BYM=BESSY0(X)
  DO J=1,N-1
     BYP=J*TOX*BY-BYM
     BYM=BY
     BY=BYP
  ENDDO
  BESSY=BY
  RETURN
END FUNCTION BESSY

!     ************************************************
!
!        BESSEL FUNCTION OF THE FIRST KIND, ORDER 0
!
!     ************************************************

FUNCTION BESSY0(X)

  implicit none
  real(8) :: BESSY0,X
!U    USES BESSJ0
  real(8) :: XX,Z,BESSJ0
  real(8) :: P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5
  real(8) :: R1,R2,R3,R4,R5,R6,S1,S2,S3,S4,S5,S6,Y
  SAVE   P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,&
       &       R1,R2,R3,R4,R5,R6,S1,S2,S3,S4,S5,S6
  DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,&
       &                         -.2073370639D-5,.2093887211D-6/
  DATA Q1,Q2,Q3,Q4,Q5/-.1562499995D-1,.1430488765D-3,&
       &        -.6911147651D-5,.7621095161D-6,-.934945152D-7/
  DATA R1,R2,R3,R4,R5,R6/-2957821389.D0,7062834065.D0,&
       &        -512359803.6D0,10879881.29D0,-86327.92757D0,228.4622733D0/
  DATA S1,S2,S3,S4,S5,S6/40076544269.D0,745249964.8D0,&
       &         7189466.438D0,47447.26470D0, 226.1030244D0,1.D0/
  
  IF(X.LT.8.D0)THEN
     Y=X**2
     BESSY0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))&
          &         /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))&
          &         +.636619772D0*BESSJ0(X)*LOG(X)
  ELSE
     Z=8.D0/X
     Y=Z**2
     XX=X-.785398164D0
     BESSY0=SQRT(.636619772D0/X)*&
          &         (  SIN(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*P5))))&
          &         +Z*COS(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
  ENDIF
  RETURN
END FUNCTION BESSY0

!     ************************************************
!
!        BESSEL FUNCTION OF THE FIRST KIND, ORDER 1
!
!     ************************************************

FUNCTION BESSY1(X)

  implicit none
  real(8) :: BESSY1,X
!U    USES BESSJ1
  real(8) :: XX,Z,BESSJ1
  real(8) :: P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5
  real(8) :: R1,R2,R3,R4,R5,R6,S1,S2,S3,S4,S5,S6,S7,Y
  SAVE   P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,&
       &       R1,R2,R3,R4,R5,R6,S1,S2,S3,S4,S5,S6,S7
  DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,&
       &                         .2457520174D-5,-.240337019D-6/
  DATA Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3,&
       &         .8449199096D-5,-.88228987D-6,.105787412D-6/
  DATA R1,R2,R3,R4,R5,R6/&
       &     -.4900604943D13,.1275274390D13,-.5153438139D11,&
       &     .7349264551D9,-.4237922726D7,.8511937935D4/
  DATA S1,S2,S3,S4,S5,S6,S7/&
       &     .2499580570D14,.4244419664D12,.3733650367D10,&
       &     .2245904002D8,.1020426050D6,.3549632885D3,1.D0/
  
  IF(X.LT.8.D0)THEN
     Y=X**2
     BESSY1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))&
          &          /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*(S6+Y*S7))))))&
          &        +.636619772D0*(BESSJ1(X)*LOG(X)-1.D0/X)
  ELSE
     Z=8.D0/X
     Y=Z**2
     XX=X-2.356194491D0
     BESSY1=SQRT(.636619772D0/X)*&
          &          (  SIN(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*P5))))&
          &          +Z*COS(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
  ENDIF
  RETURN
END FUNCTION BESSY1
