C     $Id$
C
C     ***************************************************************
C
C        MAIN ROUTINE
C
C     ***************************************************************
C
      SUBROUTINE TXEXEC
C
      INCLUDE 'txcomm.h'
C
      CHARACTER STR1*10, STR2*10, STR3*10
C
      IF (IERR .NE. 0) THEN
         WRITE(6,*) '### ERROR(TXEXEC) : Error should be cleared.'
         RETURN
      ENDIF
C
      CALL GUDATE(NDY,NDM,NDD,NTH,NTM,NTS)
      iRTIME1 = NTH * 60 * 60 + NTM * 60 + NTS
      CALL GUTIME(gCTIME1)
C
C     ***** Core of simulation *****
C
      WRITE(6,*) 'Calculating'
      CALL TXLOOP
C
C     ******************************
C
      CALL GUDATE(NDY,NDM,NDD,NTH,NTM,NTS)
      iRTIME2 = NTH * 60 * 60 + NTM * 60 + NTS
      CALL GUTIME(gCTIME2)
      iRTIME3 = iRTIME2 - iRTIME1
      if (iRTIME3 .LT. 0) iRTIME3 = iRTIME3 + 24 * 60 * 60
      gCTIME3 = gCTIME2 - gCTIME1
C
      NSTR1 = 0
      CALL APRTOS(STR1, NSTR1, gCTIME3, 'F2')
      IF (iRTIME3 .EQ. 0) THEN
         WRITE(6,*) 'real =', iRTIME3, '(sec)   ', 
     &              'CPU = ', STR1(1:NSTR1), '(sec)'
      ELSE
         NSTR2 = 0
         CALL APRTOS(STR2, NSTR2, gCTIME3 / (iRTIME3 + 1) * 100, 'F1')
         NSTR3 = 0
         CALL APRTOS(STR3, NSTR3, gCTIME3 / (iRTIME3 + 0) * 100, 'F1')
         WRITE(6,*) 'real =', iRTIME3, '(sec)   ', 
     &              'CPU = ', STR1(1:NSTR1), '(sec)   ', 
     &              '(', STR2(1:NSTR2), '% - ', STR3(1:NSTR3), '%)'
      ENDIF
C
C     ***** Print simulation results *****
C
      CALL TXWDAT
C
      RETURN
      END
C
C     ***************************************************************
C
C        Core routine of simulation
C
C     ***************************************************************
C
      SUBROUTINE TXLOOP
C
      INCLUDE 'txcomm.h'
C
      DIMENSION XP(NQM,0:NRM)
C
      IF (MODEAV. EQ. 0) THEN
         IDIV = NTMAX + 1
      ELSE
         IDIV = NTMAX / MODEAV
      ENDIF
      TIME0 = TIME
      DIP=(rIPe-rIPs)/NTMAX
C
      DO NTDO = 1, NTMAX
         NT = NTDO
         TIME = TIME0 + DT*NT
         rIP  = rIPs  + DIP*NT
C
         DO NR = 0, NRMAX
            DO NQ = 1, NQMAX
               XN(NQ,NR) = X(NQ,NR)
            ENDDO
         ENDDO
C
         IC = 0
C
   10    CONTINUE
         IC = IC + 1
C
         DO NR = 0, NRMAX
            DO NQ = 1, NQMAX
               XP(NQ,NR) = XN(NQ,NR)
            ENDDO
         ENDDO
C
         CALL TXCALV(XP)
         CALL TXCALC
         CALL TXCALA
         CALL TXCALB
         CALL TXGLOB
C
         CALL BANDRD(BA, BX, NQMAX*(NRMAX+1), 4*NQMAX-1, 4*NQM-1, IERR)
         IF (IERR .EQ. 30000) THEN
           WRITE(6,*) '### ERROR(TXLOOP) : Matrix BA is singular at ', 
     &              NT, ' -', IC, ' step.'
            IERR = 1
            DO NR = 0, NRMAX
               DO NQ = 1, NQMAX
                  XN(NQ,NR) = XP(NQ,NR)
               ENDDO
            ENDDO
            GOTO 180
         ENDIF
C
         DO NR = 0, NRMAX
            DO NQ = 1, NQMAX
               XN(NQ,NR) = BX(NQMAX * NR + NQ)
            ENDDO
         ENDDO
C
         CALL TXCHEK(NT,IC,XN,IERR)
         IF (IERR.NE.0) THEN
            DO NR = 0, NRMAX
               DO NQ = 1, NQMAX
                  X(NQ,NR) = XN(NQ,NR)
               ENDDO
            ENDDO
            CALL TXCALV(X)
            GOTO 180
         ENDIF
C
         IF (IC .GE. ICMAX) GOTO 150
C
         DO NQ = 1, NQMAX
            SUM = 0.D0
            AVM = 0.D0
            ERR1= 0.D0
            IDISP= IDIV
            DO NR = 0, NRMAX
               SUM = SUM + XN(NQ,NR)**2
               IF(SQRT(XN(NQ,NR)**2).GT.AVM) THEN
                  AVM = SQRT(XN(NQ,NR)**2)
                  NRAVM= NR
               ENDIF
            ENDDO
            AV = SQRT(SUM)/NRMAX
            IF (AV .GT. 0.D0) THEN
               DO NR = 0, NRMAX
                  ERR1 = MAX(ERR1, ABS(XN(NQ,NR) - XP(NQ,NR))/AV)
                  IF (NT. EQ. IDISP. AND. NR. EQ. NRMAX) THEN
                     IF (NQ.EQ.1) THEN
                        WRITE(6,'((1H ,A5,A3,2H =,I3))')'#####','IC',IC
                        WRITE(6,'((1H ,A4,2H =,1PD9.2),2X,A2,2H =,I3)')
     &                 'EPS   ', EPS,'NRMAX   ', NRMAX
                     ENDIF
                     WRITE(6,'((1H ,A2,2H =,I2,2X,A2,2H =,1PD9.2,
     &                 2X,A5,2H =,1PD9.2,A1,I2,2X,A5,2H =,1PD9.2))')
     &               'NQ    ', NQ ,
     &               'AV    ', AV    ,  'AVMAX ', AVM   ,':',NRAVM,
     &               'SCMAX ', ERR1
                     IDISP = IDIV + NT
                     IF (NQ. EQ. NQMAX) THEN
                        GO TO 10
                     ENDIF
                  ELSEIF (NT. NE. IDISP. AND. 
     &                     ABS(XN(NQ,NR) - XP(NQ,NR))/AV .GT. EPS) THEN
                     GO TO 10
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
C
  150    CONTINUE
C
         DO NR = 0, NRMAX
            DO NQ = 1, NQMAX
               X(NQ,NR) = XN(NQ,NR)
            ENDDO
         ENDDO
C
         CALL TXCALV(X)
         CALL TXCALC
C
         IF ((MOD(NT, NTSTEP) .EQ. 0) .AND.
     &       (NT .NE. NTMAX)) THEN
            WRITE(6,601) NT,TIME,IC
  601       FORMAT(1H ,'NT =',I4,'   T =',1PD9.2,'   IC =',I3)
         ENDIF
C
  180    IF (MOD(NT, NGRSTP) .EQ. 0) THEN
            CALL TXSTGR
         ENDIF
C
         IF (MOD(NT, NGTSTP) .EQ. 0) THEN
            CALL TXSTGT(SNGL(TIME))
         ENDIF
C
         IF (MOD(NT, NGVSTP) .EQ. 0) THEN
            CALL TXGLOB
            CALL TXSTGV(SNGL(TIME))
         ENDIF
C
         IF (IERR .NE. 0) GOTO 190
      ENDDO
C
  190 CONTINUE
      rIPs=rIPe
C
C         DO I=0,1
      WRITE(6,602) NT,TIME,IC
  602 FORMAT(1H ,'NT =',I4,'   T =',1PD9.2,'   IC =',I3)
C         ENDDO
C
      RETURN
      END
C
C     ***************************************************************
C
C        Calculate BA, BX
C
C     ***************************************************************
C
      SUBROUTINE TXCALB
C
      INCLUDE 'txcomm.h'
C
      DO I = 1, NQMAX * 4 - 1
         DO J = 1, NQMAX * (NRMAX + 1)
            BA(I,J) = 0.D0
         ENDDO
      ENDDO
C
      DO NR = 0, NRMAX
         DO NQ = 1, NQMAX
            NC = 0
               NC1 = NLC(NC,NQ,NR)
               IC = NQMAX + (NC1 - 1) - (NQ - 1)
               IB = IC + NQMAX
               IA = IB + NQMAX
               J = NR * NQMAX + NQ
               BA(IC,J)
     &       = BA(IC,J) + CLC(NC,NQ,NR)
               BA(IB,J)
     &       = BA(IB,J) + BLC(NC,NQ,NR)
               BA(IA,J)
     &       = BA(IA,J) + ALC(NC,NQ,NR)
            ENDDO
      ENDDO
C
      DO NR = 0, NRMAX
         DO NQ = 1, NQMAX
            DO NC = 1, NLCMAX(NQ)
               NC1 = NLC(NC,NQ,NR)
               IC = NQMAX + (NC1 - 1) - (NQ - 1)
               IB = IC + NQMAX
               IA = IB + NQMAX
               J = NR * NQMAX + NQ
               BA(IC,J)
     &       = BA(IC,J) - CLC(NC,NQ,NR)
               BA(IB,J)
     &       = BA(IB,J) - BLC(NC,NQ,NR)
               BA(IA,J)
     &       = BA(IA,J) - ALC(NC,NQ,NR)
            ENDDO
         ENDDO
      ENDDO
C
      DO NQ = 1, NQMAX
         DO NR = 0, NRMAX
            BX(NQMAX * NR + NQ) = 0.D0
         ENDDO
      ENDDO
C
      NC = 0
         NR = 0
            DO NQ = 1, NQMAX
               NC1 = NLC(NC,NQ,NR)
               BX(NQMAX * NR + NQ)
     &       = BX(NQMAX * NR + NQ) + BLC(NC,NQ,NR) * X(NC1,NR  )
     &                             + ALC(NC,NQ,NR) * X(NC1,NR+1)
            ENDDO
C
         DO NR = 1, NRMAX - 1
            DO NQ = 1, NQMAX
               NC1 = NLC(NC,NQ,NR)
               BX(NQMAX * NR + NQ)
     &       = BX(NQMAX * NR + NQ) + CLC(NC,NQ,NR) * X(NC1,NR-1)
     &                             + BLC(NC,NQ,NR) * X(NC1,NR  )
     &                             + ALC(NC,NQ,NR) * X(NC1,NR+1)
            ENDDO
         ENDDO
C
         NR = NRMAX
            DO NQ = 1, NQMAX
               NC1 = NLC(NC,NQ,NR)
               BX(NQMAX * NR + NQ)
     &       = BX(NQMAX * NR + NQ) + CLC(NC,NQ,NR) * X(NC1,NR-1)
     &                             + BLC(NC,NQ,NR) * X(NC1,NR  )
            ENDDO
C
      DO NR = 0, NRMAX
         DO NQ = 1, NQMAX
            DO NC = 1, NLCMAX(NQ)
                  BX(NQMAX * NR + NQ)
     &          = BX(NQMAX * NR + NQ) + PLC(NC,NQ,NR)
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C
C     ***************************************************************
C
C        Check negative density of temperature
C
C     ***************************************************************
C
      SUBROUTINE TXCHEK(NTL,IC,XL,IER)
C
      INCLUDE 'txcomm.h'
C
      DIMENSION XL(NQM,0:NRM)
C
      IER = 0
C
      DO NR = 0, NRMAX
         IF (XL(LQe1,NR).LT.0.D0 .OR. XL(LQi1,NR).LT.0.D0) THEN
            WRITE(6,*) '### ERROR(TXLOOP) : ', 
     &           'Density become negative at ', 
     &           'NR =', NR, ', ', NTL, ' -', IC, ' step.'
            WRITE(6,*) 'ne =', SNGL(XL(LQe1,NR)), 
     &             '   nNR =', SNGL(XL(LQi1,NR))
            IER = 1
            RETURN
         ENDIF
      ENDDO
C
      DO NR = 0, NRMAX
         IF (XL(LQe5,NR).LT.0.D0 .OR. XL(LQi5,NR).LT.0.D0) THEN
            WRITE(6,*) '### ERROR(TXLOOP) : ', 
     &           'Temperature become negative at ', 
     &           'NR =', NR, ', ', NTL, ' -', IC, ' step.'
            WRITE(6,*) 'Te =', SNGL(XL(LQe5,NR)), 
     &              '   Ti =', SNGL(XL(LQi5,NR))
            IER = 1
            RETURN
         ENDIF
      ENDDO
C
      RETURN
      END
C
C     ***************************************************************
C
C        Write Data
C
C     ***************************************************************
C
      SUBROUTINE TXWDAT
C
      INCLUDE 'txcomm.h'
C
C     ***** Volume-averaged density *****
C
      rNbar = 0.D0
      DO 10 NR = 0, NRMAX
         rNbar = rNbar + 2 * PI * R(NR) * PNeI(NR) * 1.D20 * DR
   10 CONTINUE
      rNbar = rNbar / (PI * RB**2)
C
      WRITE(6,'((1H ,A,1H=,1PD9.2,3(2X,A,1H=,1PD9.2)))')
     &     'Ne(0)',    PNeI(0), 
     &     'UePhi(0)', (9*X(LQe4,0)-X(LQe4,1))/(8*PNeI(0)) / 1.D3, 
     &     'UiPhi(0)', (9*X(LQi4,0)-X(LQi4,1))/(8*PNiI(0)) / 1.D3, 
     &     'N0(RB)',   (9*X(LQn1,NRMAX-1)-X(LQn1,NRMAX-2))/8 * 1.D20, 
     &     'NB(0)',    rLINEAVE(0.D0)   / 1.D20, 
     &     'NB(0.24)', rLINEAVE(0.24D0) / 1.D20, 
     &     ' NB(0.6)', rLINEAVE(0.6D0)  / 1.D20, 
     &     '    PF',   PNeI(0) * 1.D20 / rNbar
      RETURN
      END
C
C     ***************************************************************
C
C        Write Data 2
C
C     ***************************************************************
C
      SUBROUTINE TXWDAT2
C
      INCLUDE 'txcomm.h'
C
      CALL GMNMX1(GTY(0,1),  1, NGT + 1, 1,  gPNeMIN,  gPNeMAX)
      CALL GMNMX1(GTY(0,2),  1, NGT + 1, 1,  gNB0MIN,  gNB0MAX)
      CALL GMNMX1(GTY(0,11), 1, NGT + 1, 1, gUiphMIN, gUiphMAX)
      WRITE(6,'((1H ,A,1H=,1PD9.2,2(2X,A,1H=,1PD9.2)))')
     &     'MAX(Ne(0))',     gPNeMAX / 1.E20, 
     &     'MAX(NB(0))',     gNB0MAX / 1.E20, 
     &     'MAX(UiPhi(0))', gUiphMAX / 1.E3,
     &     'MIN(Ne(0))',     gPNeMIN / 1.E20, 
     &     'MIN(NB(0))',     gNB0MIN / 1.E20, 
     &     'MIN(UiPhi(0))', gUiphMIN / 1.E3
C
      RETURN
      END
C
C     ***********************************************************
C
C        LINE AVERAGE OF rN
C
C     ***********************************************************
C
      REAL*8 FUNCTION rLINEAVE(Rho)
C
      INCLUDE 'txcomm.h'
C
      D = Rho * RA
      NY = 100
      SUM = 0.D0
      DY = SQRT(RA*RA - D*D) / NY
      DO 10 I = 0, NY
         Y = DY * I
         RL = SQRT(Y*Y + D*D)
         IR = NINT(RL * NRMAX / RA)
         SUM = SUM + PNeI(IR) * 1.D20 * DY
   10 CONTINUE
      rLINEAVE = SUM / SQRT(RA*RA - D*D)
C
      RETURN
      END
