C     $Id$
C
C     ****** DENSITY PROFILE ******
C
      SUBROUTINE WMCDEN(NR,RN,RTPR,RTPP,RU)
C
      INCLUDE 'wmcomm.inc'
c      INCLUDE 'wmxprf.inc'
C
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RU(NSM)
      DIMENSION RNPL(NSM),RTPL(NSM),RUPL(NSM)
C
      RHOL=XRHO(NR)
      IF(RHOL.LE.0.D0) RHOL=0.D0
C
      IF(MODELN.EQ.0) THEN
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               RN(NS)  =0.D0
               RTPR(NS)=PTS(NS)
               RTPP(NS)=PTS(NS)
               RU(NS)  =PUS(NS)
            ENDDO
         ELSE
            FACTN=(1.D0-RHOL**PROFN1)**PROFN2
            FACTT=(1.D0-RHOL**PROFT1)**PROFT2
            FACTU=(1.D0-RHOL**PROFU1)**PROFU2
            DO NS=1,NSMAX
               RN(NS)  =((PN(NS)  -PNS(NS))*FACTN+PNS(NS))
               RTPR(NS)=((PTPR(NS)-PTS(NS))*FACTT+PTS(NS))
               RTPP(NS)=((PTPP(NS)-PTS(NS))*FACTT+PTS(NS))
               RU(NS)  = (PU(NS)  -PUS(NS))*FACTU+PUS(NS)
               IF(RHOL.LT.RHOITB) THEN
                  FACTITB=(1.D0-(RHOL/RHOITB)**4)**2
                  RN(NS)  =RN(NS)  +PNITB(NS)*FACTITB
                  RTPR(NS)=RTPR(NS)+PTITB(NS)*FACTITB
                  RTPP(NS)=RTPP(NS)+PTITB(NS)*FACTITB
                  RU(NS)  =RU(NS)  +PUITB(NS)*FACTITB
               ENDIF
            ENDDO
         ENDIF
C
      ELSEIF(MODELN.EQ.7) THEN
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               RN(NS)  =0.D0
               RTPR(NS)=RTPRF(NRMAX+1,NS)
               RTPP(NS)=RTPRF(NRMAX+1,NS)
               RU(NS)  =PUS(NS)
            ENDDO
         ELSE
            FACTU=(1.D0-RHOL**PROFU1)**PROFU2
            DO NS=1,NSMAX
               RN(NS)  = RNPRF(NR,NS)
               RTPR(NS)= RTPRF(NR,NS)
               RTPP(NS)= RTPRF(NR,NS)
               RU(NS)  = (PU(NS)  -PUS(NS))*FACTU+PUS(NS)
            ENDDO
         ENDIF
C
      ELSEIF(MODELN.EQ.8) THEN
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               RN(NS)  =0.D0
               IF (NS.EQ.1.OR.NS.GT.1) THEN
                  RTPR(NS)=PT60(NRMAX+1,NS)*1.D-3
                  RTPP(NS)=PT60(NRMAX+1,NS)*1.D-3
               ELSE
                  RTPR(NS)=PTS(NS)
                  RTPP(NS)=PTS(NS)
               ENDIF
               RU(NS)  =PUS(NS)
            ENDDO
         ELSE
            FACTN=(1.D0-RHOL**PROFN1)**PROFN2
            FACTT=(1.D0-RHOL**PROFT1)**PROFT2
            FACTU=(1.D0-RHOL**PROFU1)**PROFU2
            DO NS=1,NSMAX
               IF (NS.EQ.1.OR.NS.GT.1) THEN
                  RN(NS)  = PN60(NR,NS)*1.D-20
                  RTPR(NS)= PT60(NR,NS)*1.D-3
                  RTPP(NS)= PT60(NR,NS)*1.D-3
               ELSE
                  RN(NS)  =((PN(NS)  -PNS(NS))*FACTN+PNS(NS))
                  RTPR(NS)=((PTPR(NS)-PTS(NS))*FACTT+PTS(NS))
                  RTPP(NS)=((PTPP(NS)-PTS(NS))*FACTT+PTS(NS))
               ENDIF
               RU(NS)  = (PU(NS)  -PUS(NS))*FACTU+PUS(NS)
            ENDDO
         ENDIF
      ELSEIF(MODELN.EQ.9) THEN
         IF(RHOL.GE.1.D0) THEN
            CALL PLDATA_GETPL(1.D0,RNPL,RTPL,RUPL)
            DO NS=1,NSMAX
               RN(NS)  =0.D0
               RTPR(NS)=RTPL(NS)
               RTPP(NS)=RTPL(NS)
               RU(NS)  =PUS(NS)
            ENDDO
         ELSE
            CALL PLDATA_GETPL(RHOL,RNPL,RTPL,RUPL)
            DO NS=1,NSMAX
               RN(NS)  =RNPL(NS)
               RTPR(NS)=RTPL(NS)
               RTPP(NS)=RTPL(NS)
               RU(NS)  =RUPL(NS)
            ENDDO
         ENDIF
      ENDIF
C      WRITE(6,'(1P5E12.4)') RHOL,RN(1),RTPR(1),RN(2),RTPR(2)
C
C     ****** CALCULATION OF FAST ALPHA DENSITY AND TEMPERATURE ******
C     **** FAST ALPHA DENSITY IS SUBTRACTED FROM HELIUM DENSITY *****
C
      IF(MODELA.GE.4.AND.RHOL.LT.1.D0) THEN
         IF(MODELN.EQ.9) THEN
            NSD=3
            NST=4
            NSHE=5
            NSAF=6
         ELSE
            NSD=2
            NST=3
            NSHE=4
            NSAF=5
         ENDIF
         IF(NSM.GE.NSAF) THEN
            CALL WMCALF(RN(1),RN(NSD),RN(NST),
     &                  RTPP(1),RTPP(NSD),RTPP(NST),
     &                  RN(NSHE),RTPP(NSHE))
            RN(NSAF)=RN(NSAF)-RN(NSHE)
            RTPR(NSHE)=RTPP(NSHE)
C            WRITE(6,*) 'ALPHA N(NR),T(NR) =',
C     &         NR,RN(NSHE)/1.D20,RT(NSHE)/(1.D3*AEE)
            IF(NR.EQ.1) THEN
               PN(NSHE)  =RN(NSHE)/1.D20
               PTPR(NSHE)=RTPR(NSHE)/(1.D3*AEE)
               PTPP(NSHE)=RTPP(NSHE)/(1.D3*AEE)
            ENDIF
         ENDIF
C
      ENDIF
C
      RETURN
      END
C
C     ****** MAGNETIC FIELD PROFILE ******
C
      SUBROUTINE WMCMAG(NR,NTH,NPH,BABS,BSUPTH,BSUPPH)
C
      INCLUDE 'wmcomm.inc'
C
      BSUPTH=BFLD(2,NTH,NPH,NR)
      BSUPPH=BFLD(3,NTH,NPH,NR)
      BABS=BPST(NTH,NPH,NR)
C
      RETURN
      END
C
C     ***** CALCULATE N AND T OF FAST ALPHA PARTICLES ******
C
      SUBROUTINE WMCALF(RNE,RND,RNT,RTE,RTD,RTT,RNA,RTA)
C
      INCLUDE 'wmcomm.inc'
C
      AME=AMP*PA(1)
      FTAUS  =0.75D0*PI*SQRT(PI)*(EPS0/AEE)**2*(4.D0*AMP/AEE)*(AME/AEE)
      VALPHA=SQRT(2.D0*3.5D6*AEE/(4.D0*AMP))
C
      IF(RNE.LE.0.D0) THEN
         RNA=0.D0
         RTA=3.5D6*AEE
         RETURN
      ENDIF
C
         VTE   =SQRT(2.D0*RTE/AME)
         VCRIT3=.75D0*SQRT(PI)*AME*(RND/(AMP*PA(2))
     &                             +RNT/(AMP*PA(3)))
     &          /RNE
         VCRIT =VTE*VCRIT3**(1.D0/3.D0)
         ALNLAM=16.1D0-1.15D0*LOG10(RNE/1.D20)
     &                +2.3D0*LOG10(RTE/(AEE*1.D3))
         TAUS  =FTAUS*VTE**3/(RNE*ALNLAM)
         TI    =0.5D0*(RTD+RTT)/(AEE*1.D3)
         FH    = (TI/37.D0)
     &         + 5.45D0/(3.D0+TI/(1.D0+(TI/37.5D0)**2.8D0))
         SALPHA=RND*RNT*3.7D-18*TI**(-2.D0/3.D0)
     &          *EXP(-20.D0*TI**(-1.D0/3.D0))/FH
         RNA=SALPHA*TAUS*LOG(1.D0+(VALPHA/VCRIT)**3)/3.D0
         RTA=3.5D6*AEE*(0.5D0*(1.D0-HBEAM(VALPHA/VCRIT)))
     &                 /(LOG(1.D0+(VALPHA/VCRIT)**3)/3.D0)
C
      RETURN
      END
C
C     ****** H(X) FOR SLOWING-DOWN DISTRIBUTION ******
C
      FUNCTION HBEAM(X)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      SQR3=SQRT(3.D0)
      PI=ACOS(-1.D0)
      HBEAM=(LOG((1.D0+X**3)/(1.D0+X)**3)/3.D0
     &      +(ATAN((2.D0*X-1.D0)/SQR3)+PI/6.D0)/SQR3)/X**2
      RETURN
      END
C
C
C     ****** XL(NR): DISTANCE FROM AXIS ******
C
      SUBROUTINE WMCPOS(NR,XL)
C
      INCLUDE 'wmcomm.inc'
C
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
      END
