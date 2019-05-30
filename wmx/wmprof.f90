! wmprof.f90

MODULE wmprof

  PRIVATE
  PUBLIC wmcden,wmcmag

CONTAINS

!     ****** DENSITY PROFILE ******

  SUBROUTINE wmcden(NR,RN,RTPR,RTPP,RU)

    USE wmcomm
    USE plprof
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR
    REAL(rkind),INTENT(OUT):: RN(NSMAX),RTPR(NSMAX),RTPP(NSMAX),RU(NSMAX)
    TYPE(pl_plf_type),DIMENSION(NSMAX):: plf
    REAL(rkind):: RHON
    INTEGER:: NS,NSD,NST,NSHE,NSAF

    RHON=XRHO(NR)
    IF(RHON.LT.0.D0) RHON=0.D0

    CALL pl_prof(RHON,plf)

    DO NS=1,NSMAX
       RN(NS)  =plf(NS)%RN
       RTPR(NS)=plf(NS)%RTPR
       RTPP(NS)=plf(NS)%RTPP
       RU(NS)  =plf(NS)%RU
       PZCL(NS)=plf(NS)%RZCL
    ENDDO

!     ****** CALCULATION OF FAST ALPHA DENSITY AND TEMPERATURE ******
!     **** FAST ALPHA DENSITY IS SUBTRACTED FROM HELIUM DENSITY *****

    IF(MODELA.GE.4.AND.RHON.LT.1.D0) THEN
       NSD=2
       NST=3
       NSHE=4
       NSAF=5
       IF(NSMAX.GE.NSAF) THEN
          CALL WMCALF(RN(1),RN(NSD),RN(NST), &
                      RTPP(1),RTPP(NSD),RTPP(NST), &
                      RN(NSHE),RTPP(NSHE))
          RN(NSAF)=RN(NSAF)-RN(NSHE)
          RTPR(NSHE)=RTPP(NSHE)
!            WRITE(6,*) 'ALPHA N(NR),T(NR) =', &
!              NR,RN(NSHE)/1.D20,RT(NSHE)/(1.D3*AEE)
          IF(NR.EQ.1) THEN
             PN(NSHE)  =RN(NSHE)/1.D20
             PTPR(NSHE)=RTPR(NSHE)/(1.D3*AEE)
             PTPP(NSHE)=RTPP(NSHE)/(1.D3*AEE)
          ENDIF
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE wmcden

!     ****** MAGNETIC FIELD PROFILE ******

  SUBROUTINE wcmag(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)

    USE wmcomm
    USE plprof
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR,NTH,NHH
    REAL(rkind),INTENT(OUT):: BABS,BSUPTH,BSUPPH
    
    BSUPTH=BFLD(2,NTH,NHH,NR)
    BSUPPH=BFLD(3,NTH,NHH,NR)
    BABS=BPST(NTH,NHH,NR)

    RETURN
  END SUBROUTINE wcmag

!     ***** CALCULATE N AND T OF FAST ALPHA PARTICLES ******

  SUBROUTINE wmcalf(RNE,RND,RNT,RTE,RTD,RTT,RNA,RTA)

    USE wmcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: RNE,RND,RNT,RTE,RTD,RTT
    REAL(rkind),INTENT(OUT):: RNA,RTA
    REAL(rkind):: FTAUS,VALPHA,VTE,VCRIT3,VCRIT,ALNLAM,TAUS,TI,FH,SALPHA

    FTAUS  =0.75D0*PI*SQRT(PI)*(EPS0/AEE)**2*(4.D0*AMP/AEE)*(AME/AEE)
    VALPHA=SQRT(2.D0*3.5D6*AEE/(4.D0*AMP))

    IF(RNE.LE.0.D0) THEN
       RNA=0.D0
       RTA=3.5D6*AEE
       RETURN
    ENDIF

    VTE   =SQRT(2.D0*RTE/AME)
    VCRIT3=.75D0*SQRT(PI)*AME*(RND/(AMP*PA(2)) &
                              +RNT/(AMP*PA(3))) &
            /RNE
    VCRIT =VTE*VCRIT3**(1.D0/3.D0)
    ALNLAM=16.1D0-1.15D0*LOG10(RNE/1.D20) &
                 +2.3D0*LOG10(RTE/(AEE*1.D3))
    TAUS  =FTAUS*VTE**3/(RNE*ALNLAM)
    TI    =0.5D0*(RTD+RTT)/(AEE*1.D3)
    FH    = (TI/37.D0) &
          + 5.45D0/(3.D0+TI/(1.D0+(TI/37.5D0)**2.8D0))
    SALPHA=RND*RNT*3.7D-18*TI**(-2.D0/3.D0) &
          *EXP(-20.D0*TI**(-1.D0/3.D0))/FH
    RNA=SALPHA*TAUS*LOG(1.D0+(VALPHA/VCRIT)**3)/3.D0
    RTA=3.5D6*AEE*(0.5D0*(1.D0-HBEAM(VALPHA/VCRIT))) &
       /(LOG(1.D0+(VALPHA/VCRIT)**3)/3.D0)

    RETURN
  END SUBROUTINE wmcalf

!     ****** H(X) FOR SLOWING-DOWN DISTRIBUTION ******

  FUNCTION HBEAM(X)

    USE wmcomm,ONLY: PI,rkind
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X
    REAL(rkind):: HBEAM
    REAL(rkind):: SQR3

    SQR3=SQRT(3.D0)
    HBEAM=(LOG((1.D0+X**3)/(1.D0+X)**3)/3.D0 &
         +(ATAN((2.D0*X-1.D0)/SQR3)+PI/6.D0)/SQR3)/X**2
    RETURN
  END FUNCTION HBEAM


!     ****** XL(NR): DISTANCE FROM AXIS ******

  SUBROUTINE WMCPOS(NR,XL)

    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NR
    REAL(rkind),INTENT(OUT):: XL
    REAL(rkind):: FACT

    FACT=XRHO(NR)**2
    IF(FACT.LT.0.D0)THEN
       XL=0.D0
    ELSE
       IF(FACT.GE.1.D0)THEN
          XL=RA
       ELSE
          XL=RA*XRHO(NR)
       ENDIF
    ENDIF
    RETURN
  END SUBROUTINE WMCPOS
END MODULE wmprof
