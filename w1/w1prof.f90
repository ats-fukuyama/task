MODULE w1prof

  PRIVATE
  PUBLIC w1_prof
  PUBLIC w1_profx

CONTAINS

!     ****** SET PROFILES ******

  SUBROUTINE W1_PROF
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NX,NS
    REAL(rkind):: X,FTAUS,VALPHA,VTE,VCRIT3,VCRIT,ALNLAM,TAUS,TI,FH,SALPHA,T

    IF(NSYS.EQ.0) THEN
       DO NX =1,NXPMAX
          PROFB(NX)=RR/(RR+XAM(NX))
       END DO
    ELSE
       DO NX=1,NXPMAX
          PROFB(NX)=RR/(RR+XAM(NX))+EPSH*(XAM(NX)/RA)**2
       END DO
    END IF

    DO NS = 1 , NSMAX
       DO NX = 1,NXPMAX
          X=XAM(NX)/RA
          SELECT CASE(MODELN)
          CASE(0,1)
             IF(ABS(X).GT.1.D0) THEN
                SELECT CASE(MODELN)
                CASE(0)
                   PROFPN(NX,NS)=0.D0
                CASE(1)
                   PROFPN(NX,NS)=PNS(NS)
                END SELECT
             ELSE
                PROFPN(NX,NS)=(PN  (NS)-PNS(NS))*(1.D0-X*X)**APRFPN + PNS(NS)
             END IF
             PROFTR(NX,NS)=(PTPR(NS)-PTS(NS))*(1.D0-X*X)**APRFTR + PTS(NS)
!             PROFTR(NX,NS)= PTPR(NS)*EXP(-3.D0*X*X)
             PROFTP(NX,NS)=(PTPP(NS)-PTS(NS))*(1.D0-X*X)**APRFTP + PTS(NS)
!             PROFTP(NX,NS)= PTPP(NS)*EXP(-3.D0*X*X)
             PROFPU(NX,NS)=PU(NS)
          CASE(10,11)
             IF(ABS(X).GT.1.D0) THEN
                SELECT CASE(MODELN)
                CASE(10)
                   PROFPN(NX,NS)=0.D0
                CASE(11)
                   PROFPN(NX,NS)=PNS(NS)
                END SELECT
             ELSE
                PROFPN(NX,NS)=(PN  (NS)-PNS(NS))*(0.5D0*(1.D0-X))**APRFPN &
                             +PNS(NS)
             END IF
             PROFTR(NX,NS)=(PTPR(NS)-PTS(NS))*(0.5D0*(1.D0-X))**APRFTR &
                          + PTS(NS)
             PROFTP(NX,NS)=(PTPP(NS)-PTS(NS))*(0.5D0*(1.D0-X))**APRFTP &
                          + PTS(NS)
             PROFPU(NX,NS)=PU(NS)
          END SELECT
       END DO
    END DO

    IF(NALPHA.LT.2.OR.NSMAX.LT.4) RETURN

    NSMAX=4
    PA(4)=4
    PZ(4)=2
    FTAUS  =0.75D0*PI*SQRT(PI)*(EPS0/AEE)**2*(4.D0*AMP/AEE)*(AME/AEE)
    VALPHA=SQRT(2.D0*3.5D6*AEE/(4.D0*AMP))
    DO NX=1,NXPMAX
       VTE   =SQRT(2.D0*PROFTR(NX,1)*1.D3*AEE/AME)
       VCRIT3=.75D0*SQRT(PI)*AME*(PROFPN(NX,2)/(AMP*PA(2)) &
                                 +PROFPN(NX,3)/(AMP*PA(3))) &
              /PROFPN(NX,1)
       VCRIT =VTE*VCRIT3**(1.D0/3.D0)
       ALNLAM=16.1-1.15*LOG10(PROFPN(NX,1))+2.3*LOG10(PROFTR(NX,1))
       TAUS  =FTAUS*VTE**3/(PROFPN(NX,1)*1.D20*ALNLAM)
       TI    =0.5D0*(PROFTR(NX,2)+PROFTR(NX,3))
       FH    = (TI/37.D0)+5.45D0/(3.D0+TI/(1.D0+(TI/37.5)**2.8))
       SALPHA=PROFPN(NX,2)*PROFPN(NX,3)*3.7D22*TI**(-2.D0/3.D0) &
              *EXP(-20.D0*TI**(-1.D0/3.D0))/FH
       PROFPN(NX,4)=SALPHA*TAUS*LOG(1.D0+(VALPHA/VCRIT)**3)/3.D20
       T=3.5D3*(0.5D0*(1.D0-HBEAM(VALPHA/VCRIT))) &
              /(LOG(1.D0+(VALPHA/VCRIT)**3)/3.D0)
       PROFTP(NX,4)=T
       PROFTR(NX,4)=T
    END DO
    PN(4)=PROFPN(NXPMAX/2,4)
    PNS(4)=1.D-6*PN(4)
    PTPR(4)=PROFTR(NXPMAX/2,4)
    PTPP(4)=PROFTP(NXPMAX/2,4)
    PTS(4) =PROFTP(NXPMAX,4)
    WRITE(6,601) PN(4)
601 FORMAT('** CENTRAL ALPHA DENSITY = ',1PE12.4)
    RETURN
  END SUBROUTINE W1_PROF

  SUBROUTINE W1_PROFX
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NX,NS
    REAL(rkind):: X,FTAUS,VALPHA,VTE,VCRIT3,VCRIT,ALNLAM,TAUS,TI,FH,SALPHA,T

    IF(NSYS.EQ.0) THEN
       DO NX =1,NXMAX
          PROFB(NX)=RR/(RR+XAM(NX))
       END DO
    ELSE
       DO NX=1,NXMAX
          PROFB(NX)=RR/(RR+XAM(NX))+EPSH*(XAM(NX)/RA)**2
       END DO
    END IF

    DO NS = 1,NSMAX
       DO NX = 1,NXMAX
          X=XAM(NX)/RA
          IF(ABS(X).LE.1.D0) THEN
             PROFPN(NX,NS)=(PN  (NS)-PNS(NS))*(1.D0-X*X)**APRFPN + PNS(NS)
             PROFTR(NX,NS)=(PTPR(NS)-PTS(NS))*(1.D0-X*X)**APRFTR + PTS(NS)
             PROFTP(NX,NS)=(PTPP(NS)-PTS(NS))*(1.D0-X*X)**APRFTP + PTS(NS)
!            PROFTR(NX,NS)= PTPR(NS)*EXP(-3.D0*X*X)
!            PROFTP(NX,NS)= PTPP(NS)*EXP(-3.D0*X*X)
             PROFPU(NX,NS)=PU(NS)
          ELSE
             PROFPN(NX,NS)=0.D0
             PROFTR(NX,NS)=PTS(NS)
             PROFTP(NX,NS)=PTS(NS)
             PROFPU(NX,NS)=PU(NS)
          END IF
       END DO
    END DO

    IF(NALPHA.LT.2.OR.NSMAX.LT.4) RETURN

    NSMAX=4
    PA(4)=4
    PZ(4)=2
    FTAUS  =0.75D0*PI*SQRT(PI)*(EPS0/AEE)**2*(4.D0*AMP/AEE)*(AME/AEE)
    VALPHA=SQRT(2.D0*3.5D6*AEE/(4.D0*AMP))
    DO NX=1,NXMAX
       VTE   =SQRT(2.D0*PROFTR(NX,1)*1.D3*AEE/AME)
       VCRIT3=.75D0*SQRT(PI)*AME*(PROFPN(NX,2)/(AMP*PA(2)) &
                                 +PROFPN(NX,3)/(AMP*PA(3))) &
              /PROFPN(NX,1)
       VCRIT =VTE*VCRIT3**(1.D0/3.D0)
       ALNLAM=16.1-1.15*LOG10(PROFPN(NX,1))+2.3*LOG10(PROFTR(NX,1))
       TAUS  =FTAUS*VTE**3/(PROFPN(NX,1)*1.D20*ALNLAM)
       TI    =0.5D0*(PROFTR(NX,2)+PROFTR(NX,3))
       FH    = (TI/37.D0)+5.45D0/(3.D0+TI/(1.D0+(TI/37.5)**2.8))
       SALPHA=PROFPN(NX,2)*PROFPN(NX,3)*3.7D22*TI**(-2.D0/3.D0) &
              *EXP(-20.D0*TI**(-1.D0/3.D0))/FH
       PROFPN(NX,4)=SALPHA*TAUS*LOG(1.D0+(VALPHA/VCRIT)**3)/3.D20
       T=3.5D3*(0.5D0*(1.D0-HBEAM(VALPHA/VCRIT))) &
              /(LOG(1.D0+(VALPHA/VCRIT)**3)/3.D0)
       PROFTP(NX,4)=T
       PROFTR(NX,4)=T
    END DO
    PN(4)=PROFPN(NXMAX/2,4)
    PNS(4)=1.D-6*PN(4)
    PTPR(4)=PROFTR(NXMAX/2,4)
    PTPP(4)=PROFTP(NXMAX/2,4)
    PTS(4) =PROFTP(NXMAX,4)
    WRITE(6,601) PN(4)
601 FORMAT('** CENTRAL ALPHA DENSITY = ',1PE12.4)
    RETURN
  END SUBROUTINE W1_PROFX

!     H(X) FOR SLOWING-DOWN DISTRIBUTION

  FUNCTION HBEAM(X)
    USE w1comm,ONLY: rkind,PI
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X
    REAL(rkind):: HBEAM,SQR3

    SQR3=SQRT(3.D0)
    HBEAM=(LOG((1.D0+X**3)/(1.D0+X)**3)/3.D0 &
           +(ATAN((2.D0*X-1.D0)/SQR3)+PI/6.D0)/SQR3)/X**2
    RETURN
  END FUNCTION HBEAM
END MODULE w1prof
