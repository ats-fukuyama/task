!     ****** CALCULATE LOCAL MAGNETIC FIELD ******

      SUBROUTINE PLMAG(X,Y,Z,RHON)

!      INCLUDE '../pl/plcomm.inc'
!      INCLUDE '../pl/plcom2.inc'
      USE PLLOC1
      USE PLPRM1, ONLY : RA, RR
      USE PLPRM2, ONLY : BB
      USE PLPRM9, ONLY : MODELG

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: X, Y, Z
      REAL(8), INTENT(OUT):: RHON
      REAL(8) :: BP, BR, BT, BX, BY, BZ, PP, QL, RL, RS, &
     &           RSINT, RCOST, RSINP, RCOSP

      IF(MODELG.EQ.0) THEN
         RS   = X-RR
         RHON = RS/RA
         CALL PLQPRF(RHON,QL)
         IF(RS.LE.0.D0) QL=-QL
         BT   = BB
         BP   = RS*BT/(RR*QL)
         BX   = 0.D0
         BY   = BB
         BZ   = BP

      ELSEIF(MODELG.EQ.1) THEN
         RS =SQRT((X-RR)**2+Z**2)
         RHON=RS/RA
         IF(RS.LE.0.D0) THEN
            BX   = 0.D0
            BY   = BB
            BZ   = 0.D0
         ELSE
            CALL PLQPRF(RHON,QL)
            RSINT= Z/RS
            RCOST= (X-RR)/RS
            BT   = BB
            BP   = RS*BT/(RR*QL)
            BX   =-BP*RSINT
            BY   = BB
            BZ   = BP*RCOST
         ENDIF

      ELSEIF(MODELG.EQ.2) THEN
         RL=SQRT(X**2+Y**2)
         RS =SQRT((RL-RR)**2+Z**2)
         RHON=RS/RA
         IF(RS.LE.0.D0) THEN
            BT   = BB
            BR   = 0.D0
            BZ   = 0.D0
         ELSE
            CALL PLQPRF(RHON,QL)
            RSINP= Z/RS
            RCOSP= (RL-RR)/RS
            BT   = BB/(1.D0+RS*RCOSP/RR)
            BP   = RS*BT/(RR*QL)
            BR   =-BP*RSINP
            BZ   = BP*RCOSP
         ENDIF
         RCOST=X/RL
         RSINT=Y/RL
         BX = BR*RCOST-BT*RSINT
         BY = BR*RSINT+BT*RCOST

      ELSEIF(MODELG.EQ.3.OR.MODELG.EQ.5) THEN
         RL=SQRT(X**2+Y**2)
         PP=0.D0
         CALL GETRZ(RL,Z,PP,BR,BZ,BT,RHON)
!         WRITE(6,'(1P6E12.4)') RL,ZZ,BR,BZ,BT,RHON
         RCOST=X/RL
         RSINT=Y/RL
         BX = BR*RCOST-BT*RSINT
         BY = BR*RSINT+BT*RCOST
      ENDIF

      BABS = SQRT(BX**2+BY**2+BZ**2)

      IF(BABS.LE.0.D0) THEN
         BNX = 0.D0
         BNY = 1.D0
         BNZ = 0.D0
      ELSE
         BNX = BX/BABS
         BNY = BY/BABS
         BNZ = BZ/BABS
      ENDIF

      RETURN
      END

!     ****** CALCULATE PLASMA PROFILE ******

      SUBROUTINE PLPROF(RHON)

!      INCLUDE '../pl/plcomm.inc'
!      INCLUDE '../pl/plcom2.inc'
      USE PLCOM0
      USE PLCNS1, ONLY : AEE
!      USE PLDAT2
      USE PLLOC2
      USE PLPRM3
      USE PLPRM4, ONLY : PN, PNS ,PZCL
      USE PLPRM5
      USE PLPRM6
      USE PLPRM7, ONLY : RHOEDG
      USE PLPRM8
      USE PLPRM9, ONLY : MODELN

      IMPLICIT NONE
      REAL(8), INTENT(IN):: RHON
      REAL(8)     :: RHOL, FACTN, FACTT, FACTU, FACTITB, PTOT, PL0, PSIN, PL, &
     &               FACT, FNX, DFNX, AN, BN, FTX, DFTX, AT, BT, FUX, DFUX, AU, BU
      INTEGER(4)  :: NS
      REAL(8), DIMENSION(NSM) :: RNPL, RTPL, RUPL

      IF(RHON.LE.0.D0) THEN
         RHOL=0.D0
      ELSEIF(RHON.GE.1.D0) THEN
         RHOL=1.D0
      ELSE
         RHOL=RHON
      ENDIF

      IF(MODELN.EQ.0) THEN
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               RN(NS)  =0.D0
               RTPR(NS)=PTS(NS)
               RTPP(NS)=PTS(NS)
               RU(NS)  =PUS(NS)
               RZCL(NS)=PZCL(NS)
            ENDDO
         ELSE
            FACTN=(1.D0-RHOL**PROFN1)**PROFN2
            FACTT=(1.D0-RHOL**PROFT1)**PROFT2
            FACTU=(1.D0-RHOL**PROFU1)**PROFU2

            DO NS=1,NSMAX
               RN(NS)  =(PN(NS)  -PNS(NS))*FACTN+PNS(NS)
               RTPR(NS)=(PTPR(NS)-PTS(NS))*FACTT+PTS(NS)
               RTPP(NS)=(PTPP(NS)-PTS(NS))*FACTT+PTS(NS)
               RU(NS)  =(PU(NS)  -PUS(NS))*FACTU+PUS(NS)
               RZCL(NS)=PZCL(NS)
               IF(RHOL.LT.RHOITB) THEN
                  FACTITB=(1.D0-(RHOL/RHOITB)**4)**2
                  RN(NS)  =RN(NS)  +PNITB(NS)*FACTITB
                  RTPR(NS)=RTPR(NS)+PTITB(NS)*FACTITB
                  RTPP(NS)=RTPP(NS)+PTITB(NS)*FACTITB
                  RU(NS)  =RU(NS)  +PUITB(NS)*FACTITB
               ENDIF
            ENDDO
         ENDIF

      ELSEIF(MODELN.EQ.1) THEN
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               RN(NS)  =PNS(NS)
               RTPR(NS)=PTS(NS)
               RTPP(NS)=PTS(NS)
               RU(NS)  =PUS(NS)
               RZCL(NS)=PZCL(NS)
            ENDDO
         ELSE
            FACTN=(1.D0-RHOL**PROFN1)**PROFN2
            FACTT=(1.D0-RHOL**PROFT1)**PROFT2
            FACTU=(1.D0-RHOL**PROFU1)**PROFU2

            DO NS=1,NSMAX
               RN(NS)  =(PN(NS)  -PNS(NS))*FACTN+PNS(NS)
               RTPR(NS)=(PTPR(NS)-PTS(NS))*FACTT+PTS(NS)
               RTPP(NS)=(PTPP(NS)-PTS(NS))*FACTT+PTS(NS)
               RU(NS)  =(PU(NS)  -PUS(NS))*FACTU+PUS(NS)
               RZCL(NS)=PZCL(NS)
               IF(RHOL.LT.RHOITB) THEN
                  FACTITB=(1.D0-(RHOL/RHOITB)**4)**2
                  RN(NS)  =RN(NS)  +PNITB(NS)*FACTITB
                  RTPR(NS)=RTPR(NS)+PTITB(NS)*FACTITB
                  RTPP(NS)=RTPP(NS)+PTITB(NS)*FACTITB
                  RU(NS)  =RU(NS)  +PUITB(NS)*FACTITB
               ENDIF
            ENDDO
         ENDIF

      ELSEIF(MODELN.EQ.2) THEN
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               RN(NS)  =PNS(NS)
               RTPR(NS)=PTS(NS)
               RTPP(NS)=PTS(NS)
               RU(NS)  =PUS(NS)
               RZCL(NS)=PZCL(NS)
            ENDDO
         ELSE
            PTOT=PTOT*1.D20*1.D3/AEE
            CALL GETPP(0.D0,PL0)
            CALL GETPP(PSIN,PL)
            FACT=SQRT(PTOT*PL/PL0)
            FACTU=(1.D0-RHOL**PROFU1)**PROFU2
            DO NS=1,NSMAX
               RN(NS)  =PN(NS)  *FACT
               RTPR(NS)=PTPR(NS)*FACT
               RTPP(NS)=PTPP(NS)*FACT
               RU(NS)  =(PU(NS)  -PUS(NS))*FACTU+PUS(NS)
               RZCL(NS)=PZCL(NS)
            ENDDO
         ENDIF

      ELSEIF(MODELN.EQ.3) THEN
         IF(RHOL.GE.1.D0) THEN
            DO NS=1,NSMAX
               RN(NS)  =PNS(NS)
               RTPR(NS)=PTS(NS)
               RTPP(NS)=PTS(NS)
               RU(NS)  =PUS(NS)
               RZCL(NS)=PZCL(NS)
            ENDDO
         ELSE
            IF(RHOL.LE.RHOEDG) THEN
               FACTN=(1.D0-RHOL**PROFN1)**PROFN2
               FACTT=(1.D0-RHOL**PROFT1)**PROFT2
               FACTU=(1.D0-RHOL**PROFU1)**PROFU2
            ELSE
               FNX=(1.D0-RHOEDG**PROFN1)**PROFN2
               DFNX=-PROFN1*PROFN2*RHOEDG**(PROFN1-1.D0) &
     &             *(1.D0-RHOEDG**PROFN1)**(PROFN2-1.D0)
               AN= 3*FNX/(1.D0-RHOEDG)**2+DFNX/(1.D0-RHOEDG)
               BN=-2*FNX/(1.D0-RHOEDG)**3-DFNX/(1.D0-RHOEDG)**2
               FACTN=AN*(1-RHOL)**2+BN*(1-RHOL)**3

               FTX=(1.D0-RHOEDG**PROFT1)**PROFT2
               DFTX=-PROFT1*PROFT2*RHOEDG**(PROFT1-1.D0) &
     &             *(1.D0-RHOEDG**PROFT1)**(PROFT2-1.D0)
               AT= 3*FTX/(1.D0-RHOEDG)**2+DFTX/(1.D0-RHOEDG)
               BT=-2*FTX/(1.D0-RHOEDG)**3-DFTX/(1.D0-RHOEDG)**2
               FACTT=AT*(1-RHOL)**2+BT*(1-RHOL)**3

               FUX=(1.D0-RHOEDG**PROFU1)**PROFU2
               DFUX=-PROFU1*PROFU2*RHOEDG**(PROFU1-1.D0) &
     &             *(1.D0-RHOEDG**PROFU1)**(PROFU2-1.D0)
               AU= 3*FUX/(1.D0-RHOEDG)**2+DFUX/(1.D0-RHOEDG)
               BU=-2*FUX/(1.D0-RHOEDG)**3-DFUX/(1.D0-RHOEDG)**2
               FACTU=AU*(1-RHOL)**2+BU*(1-RHOL)**3
            ENDIF

            DO NS=1,NSMAX
               RN(NS)  =(PN(NS)  -PNS(NS))*FACTN+PNS(NS)
               RTPR(NS)=(PTPR(NS)-PTS(NS))*FACTT+PTS(NS)
               RTPP(NS)=(PTPP(NS)-PTS(NS))*FACTT+PTS(NS)
               RU(NS)  =(PU(NS)  -PUS(NS))*FACTU+PUS(NS)
               RZCL(NS)=PZCL(NS)
               IF(RHOL.LT.RHOITB) THEN
                  FACTITB=(1.D0-(RHOL/RHOITB)**4)**2
                  RN(NS)  =RN(NS)  +PNITB(NS)*FACTITB
                  RTPR(NS)=RTPR(NS)+PTITB(NS)*FACTITB
                  RTPP(NS)=RTPP(NS)+PTITB(NS)*FACTITB
                  RU(NS)  =RU(NS)  +PUITB(NS)*FACTITB
               ENDIF
            ENDDO
         ENDIF

      ELSEIF(MODELN.EQ.9) THEN
         CALL PLDATA_GETPL(RHOL,RNPL,RTPL,RUPL)
         DO NS=1,NSMAX
            RN(NS)  =RNPL(NS)
            RTPR(NS)=RTPL(NS)
            RTPP(NS)=RTPL(NS)
            RU(NS)  =RUPL(NS)
            RZCL(NS)=0.D0
         ENDDO
         IF(NSMAX.EQ.6) THEN
            RN(2)=RN(2)-RN(5)
            RN(4)=RN(4)-RN(6)
         ENDIF
      ENDIF

      RETURN
      END

!     ****** CALCULATE Q PROFILE ******

      SUBROUTINE PLQPRF(RHON,QL)

!      INCLUDE '../pl/plcomm.inc'
      USE PLCNS1, ONLY : PI, RMU0
      USE PLPRM1, ONLY : RA, RR
      USE PLPRM2
      USE PLPRM7, ONLY : RHOMIN, QMIN
      USE PLPRM9
      IMPLICIT NONE
      REAL(8), INTENT(IN)  :: RHON
      REAL(8), INTENT(OUT) :: QL
      REAL(8) :: RHOL, QSA0, QSAA

      IF(RHON.LE.0.D0) THEN
         RHOL=0.D0
      ELSE
         RHOL=RHON
      ENDIF

      IF(MODELG.LE.2) THEN
         IF(MODELQ.EQ.0) THEN
            IF(RHOL.GT.1.D0) THEN
               QL = QA*RHOL**2
            ELSEIF(RHOMIN.LE.0.D0)THEN
               QL =(Q0-QA)*(1-RHOL**2)+QA
            ELSE
               QSA0    =1/Q0-1/QMIN
               QSAA    =1/QA-1/QMIN
               IF(RHOL.LE.RHOMIN)THEN
                  QL =1/(1/Q0-QSA0*(3*RHOL**2/RHOMIN**2 &
     &                       -2*RHOL**3/RHOMIN**3))
               ELSE
                  QL =1/(1/QMIN+3*QSA0*(RHOL-RHOMIN)**2/RHOMIN**2 &
     &                  +(QSAA-3*QSA0*(1-RHOMIN)**2/RHOMIN**2) &
     &                    *(RHOL-RHOMIN)**3/(1-RHOMIN)**3)
               ENDIF
            ENDIF
         ELSEIF(MODELQ.EQ.1) THEN
            QA=2.D0*PI*RA*RA*BB/(RMU0*RIP*1.D6*RR)
            Q0=QA/(1.D0+PROFJ)
            IF(RHOL.GE.1.D0) THEN
               QL = QA*RHOL**2
            ELSEIF(RHOL.LE.1.D-30) THEN
               QL = Q0
            ELSE
               QL=QA*RHOL**2/(1.D0-(1.D0-RHOL**2)**(PROFJ+1.D0))
            ENDIF
         ENDIF
      ELSE
         CALL GETQP(RHOL,QL)
      ENDIF
      RETURN
      END

!     ****** CALCULATE BMIN ON MAG SURFACE ******

      SUBROUTINE PLBMIN(RHON,BMINL)

!      INCLUDE '../pl/plcomm.inc'
      USE PLPRM9, ONLY : MODELG
      USE PLPRM2, ONLY : BB
      USE PLPRM1, ONLY : RR
      IMPLICIT NONE
      REAL(8), INTENT(IN)  :: RHON
      REAL(8), INTENT(OUT) :: BMINL
      REAL(8) :: RS, BMINP, BMINT, QL, RRMAXL, BTL, BPL
      REAL(8) :: RSRHON

      IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
         RS=RSRHON(RHON)
         BMINT= BB
         CALL GETQP(RHON,QL)
         BMINP= RS*BMINT/(RR*QL)
         BMINL= SQRT(BMINT**2+BMINP**2)
      ELSEIF(MODELG.EQ.2) THEN
         RS=RSRHON(RHON)
         BMINT= BB/(1+RS/RR)
         CALL PLQPRF(RHON,QL)
         BMINP= RS*BMINT/((RR+RS)*QL)
         BMINL= SQRT(BMINT**2+BMINP**2)
      ELSEIF(MODELG.EQ.3) THEN
         CALL GETRMX(RHON,RRMAXL)
         BTL=BB*RR/RRMAXL
         CALL PLQPRF(RHON,QL)
         BPL=RS*BTL/(RR*QL)
         BMINL=SQRT(BTL**2+BPL**2)
      ENDIF
      RETURN
      END

!     ****** CALCULATE RRMIN AND RRMAX ON MAG SURFACE ******

      SUBROUTINE PLRRMX(RHON,RRMINL,RRMAXL)

!      INCLUDE '../pl/plcomm.inc'
      USE PLPRM9, ONLY : MODELG
      USE PLPRM1, ONLY : RR
      IMPLICIT NONE
      REAL(8), INTENT(IN):: RHON
      REAL(8), INTENT(OUT)  :: RRMINL,RRMAXL
      REAL(8)  :: RS
      REAL(8) :: RSRHON

      IF(MODELG.EQ.0.OR.MODELG.EQ.1) THEN
         RRMINL=RR
         RRMAXL=RR
      ELSEIF(MODELG.EQ.2) THEN
         RS=RSRHON(RHON)
         RRMINL=RR-RS
         RRMAXL=RR+RS
      ELSEIF(MODELG.EQ.3) THEN
         CALL GETRMN(RHON,RRMINL)
         CALL GETRMX(RHON,RRMAXL)
      ENDIF
      RETURN
      END

!     ***** AVERAGE MINOR RADIUS FOR PARABOLIC Q PROFILE *****

      REAL(8) FUNCTION RSRHON(RHON)

!      INCLUDE '../pl/plcomm.inc'
      USE PLPRM1, ONLY : RA
      USE PLPRM9, ONLY : MODELG
      IMPLICIT NONE
      REAL(8), INTENT(IN):: RHON
      REAL(8)  :: RHOL, RRMINL, RRMAXL

      IF(MODELG.LT.3) THEN
         IF(RHON.LE.0.D0) THEN
           RHOL=0.D0
         ELSE
           RHOL=RHON
         ENDIF
         RSRHON=RHOL*RA
      ELSEIF(MODELG.EQ.3) THEN
         CALL GETRMN(RHON,RRMINL)
         CALL GETRMX(RHON,RRMAXL)
         RSRHON=0.5D0*(RRMAXL-RRMINL)
      ENDIF
      RETURN
      END

!     ***** PLASMA MAGNETIC AXIS *****

      SUBROUTINE PLAXIS(RAXIS,ZAXIS)

!      INCLUDE '../pl/plcomm.inc'
      USE PLPRM9, ONLY : MODELG
      USE PLPRM1, ONLY : RR
      IMPLICIT NONE
      REAL(8), INTENT(OUT) :: RAXIS,ZAXIS


      IF(MODELG.LT.3) THEN
         RAXIS=RR
         ZAXIS=0.D0
      ELSEIF(MODELG.EQ.3) THEN
         CALL GETAXS(RAXIS,ZAXIS)
      ENDIF
      RETURN
      END

!     ***** PLASMA BOUNDARY *****

      SUBROUTINE PLRZSU(RSU,ZSU,NM,NSUMAX)

!      INCLUDE '../pl/plcomm.inc'

      USE PLCNS1, ONLY : PI
      USE PLPRM9, ONLY : MODELG
      USE PLPRM1, ONLY : RA, RKAP, RR
      IMPLICIT NONE
      REAL(8), DIMENSION(NM), INTENT(OUT) :: RSU, ZSU
      INTEGER(4), INTENT(IN) :: NM
      INTEGER(4), INTENT(OUT):: NSUMAX
      REAL(8)     :: DTH, TH
      INTEGER(4)  :: NSU


      IF(MODELG.EQ.0)THEN
         NSUMAX=5
         RSU(1)=RR-RA
         ZSU(1)=  -RA
         RSU(2)=RR+RA
         ZSU(2)=  -RA
         RSU(3)=RR+RA
         ZSU(3)=  +RA
         RSU(4)=RR-RA
         ZSU(4)=  +RA
         RSU(5)=RR-RA
         ZSU(5)=  -RA
      ELSEIF(MODELG.EQ.1.OR.MODELG.EQ.2) THEN
         NSUMAX=MIN(NM,101)
         DTH=2.D0*PI/(NSUMAX-1)
         DO NSU=1,NSUMAX
            TH=(NSU-1)*DTH
            RSU(NSU)=RR+     RA*COS(TH)
            ZSU(NSU)=   RKAP*RA*SIN(TH)
         ENDDO
      ELSEIF(MODELG.EQ.3) THEN
         CALL GETRSU(RSU,ZSU,NM,NSUMAX)
      ENDIF
      RETURN
      END
