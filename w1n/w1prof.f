C
C     ****** TASK/W1 (SUBROUTINE LIBRARY 1) ******
C
C     ****** SET PROFILES ******
C
      SUBROUTINE W1PROF
C
      USE w1comm
      IMPLICIT NONE
      INTEGER:: NX,NS
      REAL(rkind):: X,FTAUS,VALPHA,VTE,VCRIT3,VCRIT,ALNLAM,TAUS,TI
      REAL(rkind):: FH,SALPHA,T
      INTERFACE 
         FUNCTION HBEAM(X)
         USE w1comm
         REAL(rkind):: X
         REAL(rkind):: HBEAM
         END FUNCTION
      END INTERFACE
C
      IF(NSYS.EQ.0) THEN
         DO 100 NX =1,NXP
            PROFB(NX)=RR/(RR+XAM(NX))
  100    CONTINUE
      ELSE
         DO 110 NX=1,NXP
            PROFB(NX)=RR/(RR+XAM(NX))+EPSH*(XAM(NX)/RA)**2
  110    CONTINUE
      ENDIF
C
      DO 200 NS = 1 , NSMAX
      DO 200 NX = 1 , NXP
         X=XAM(NX)/RA
         PROFPN(NX,NS)=(PN  (NS)-PNS(NS))*(1.D0-X*X)**APRFPN + PNS(NS)
         PROFTR(NX,NS)=(PTPR(NS)-PTS(NS))*(1.D0-X*X)**APRFTR + PTS(NS)
C        PROFTR(NX,NS)= PTPR(NS)*EXP(-3.D0*X*X)
         PROFTP(NX,NS)=(PTPP(NS)-PTS(NS))*(1.D0-X*X)**APRFTP + PTS(NS)
C        PROFTP(NX,NS)= PTPP(NS)*EXP(-3.D0*X*X)
         PROFPU(NX,NS)=PU(NS)
 200  CONTINUE
C
      IF(NALPHA.LT.2.OR.NSM.LT.4) RETURN
C
      NSMAX=4
      PA(4)=4
      PZ(4)=2
      FTAUS  =0.75D0*PI*SQRT(PI)*(EPS0/AEE)**2*(4.D0*AMP/AEE)*(AME/AEE)
      VALPHA=SQRT(2.D0*3.5D6*AEE/(4.D0*AMP))
      DO 300 NX=1,NXP
         VTE   =SQRT(2.D0*PROFTR(NX,1)*1.D3*AEE/AME)
         VCRIT3=.75D0*SQRT(PI)*AME*(PROFPN(NX,2)/(AMP*PA(2))
     &                             +PROFPN(NX,3)/(AMP*PA(3)))
     &          /PROFPN(NX,1)
         VCRIT =VTE*VCRIT3**(1.D0/3.D0)
         ALNLAM=16.1-1.15*LOG10(PROFPN(NX,1))+2.3*LOG10(PROFTR(NX,1))
         TAUS  =FTAUS*VTE**3/(PROFPN(NX,1)*1.D20*ALNLAM)
         TI    =0.5D0*(PROFTR(NX,2)+PROFTR(NX,3))
         FH    = (TI/37.D0)+5.45D0/(3.D0+TI/(1.D0+(TI/37.5)**2.8))
         SALPHA=PROFPN(NX,2)*PROFPN(NX,3)*3.7D22*TI**(-2.D0/3.D0)
     &          *EXP(-20.D0*TI**(-1.D0/3.D0))/FH
         PROFPN(NX,4)=SALPHA*TAUS*LOG(1.D0+(VALPHA/VCRIT)**3)/3.D20
         T=3.5D3*(0.5D0*(1.D0-HBEAM(VALPHA/VCRIT)))
     &          /(LOG(1.D0+(VALPHA/VCRIT)**3)/3.D0)
         PROFTP(NX,4)=T
         PROFTR(NX,4)=T
  300 CONTINUE
      PN(4)=PROFPN(NXP/2,4)
      PNS(4)=1.D-6*PN(4)
      PTPR(4)=PROFTR(NXP/2,4)
      PTPP(4)=PROFTP(NXP/2,4)
      PTS(4) =PROFTP(NXP,4)
      WRITE(6,601) PN(4)
  601 FORMAT(1H ,'** CENTRAL ALPHA DENSITY = ',1PE12.4)
      RETURN
      END
C
C     H(X) FOR SLOWING-DOWN DISTRIBUTION
C
      FUNCTION HBEAM(X)
C
      USE w1comm
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X
      REAL(rkind):: HBEAM,SQR3
C
      SQR3=SQRT(3.D0)
      HBEAM=(LOG((1.D0+X**3)/(1.D0+X)**3)/3.D0
     &      +(ATAN((2.D0*X-1.D0)/SQR3)+PI/6.D0)/SQR3)/X**2
      RETURN
      END
