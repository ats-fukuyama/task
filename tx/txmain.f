C     $Id$
C
C        ***********************************************************
C
C                        TOKAMAK TRANSPORT SIMULATION
C            INCLUDING RADIAL ELECTRIC FIELD AND PLASMA ROTATION
C
C             DEVELOPED BY Y. FUJI, S. TAKATSUKA AND A. FUKUYAMA
C
C                            OKAYAMA UNIVERSITY
C                            OKAYAMA 700, JAPAN
C
C                             tx2.08 : 96/12/29
C                             tx2.10 : 98/11/21
C
C        ***********************************************************
C
C        X1  = Er                                       integer
C        X2  = Etheta                                   integer
C        X3  = Ephi                                half integer
C        X4  = Btheta                                   integer
C        X5  = Bphi                                half integer
C        X6  = Ne             X11  = Ni            half integer
C        X7  = Ne * Uer       X12  = Ni * Uir           integer
C        X8  = Ne * UeTheta   X13  = Ni * UiTheta       integer
C        X9  = Ne * UePhi     X14  = Ni * UiPhi    half integer
C        X10 = Ne * Te        X15  = Ni * Ti       half integer
C        X16 = Nb                                  half integer
C        X17 = Nb * UbTheta                             integer
C        X18 = Nb * UbPhi                          half integer
C        X19 = Slow n01                            half integer
C        X20 = Fast n02                            half integer
C
C        G16 = Total n0                            half integer
C        G20 = Q                                        integer
C        G21 = Jphi                                     integer
C        G22 = JePhi     G23 = JiPhi                    integer
C        G24 = JbPhi                                    integer
C        G25 = UePerp    G26 = UePara                   integer
C        G27 = UiPerp    G28 = UiPara                   integer
C        G29 = D    G30 = G1h2    G31 = MuI        half integer
C        G32 = S    G33 = Alpha   G34 = rKappa     half integer
C        G35 = Slow N0   G36 = Fast N0             half integer
C
      INCLUDE 'txcomm.inc'
C
      CALL GSOPEN
C
C     ***** Version ID *****
C     SLID is used to identify data file.
      SLID = 'tx211.30'
      WRITE(6,*) '######## /TASK/TX V2.11.40 00/02/14 ########'
C
      CALL TXINIT
      CALL TXPARF
C
      CALL TXMAIN
C
      CALL GSCLOS
      STOP
      END
C
C     ***************************************************************
C
C        Command loop
C
C     ***************************************************************
C
      SUBROUTINE TXMAIN
C
      INCLUDE 'txcomm.inc'
C
      CHARACTER STR*40, STR2*5, KID*1, KID2*2, KID3*3
C
      ICONT = 0
C
C     *** MENU ***
C
   10 CONTINUE
      TMAX=TIME+DT*NTMAX
C
      WRITE(6,*)
     &   'R:RUN  C:CONTINUE  P,V:PARAMETER  I:INIT  Q:QUIT '//
     &   'S:SAVE  L:LOAD'
      WRITE(6,'(1H ,A,I2,3(A,1PD9.2),A,I2)')
     &   'G[T|U|V]1-', MAX(NGPRM,NGPTM,NGPVM), ':GRAPH  TIME=', TIME, 
     &   '  DT=', DT, '  TMAX =', TMAX
      WRITE(6,*) '## INPUT:'
      CALL GUFLSH
      READ(*,'(A40)',END=9000) STR
C
      IER = 0
C
   20 KID  = STR(1:1)
      KID2 = STR(1:2)
      KID3 = STR(1:3)
      CALL GUCPTL(KID)
      CALL TOUPPER(KID2)
      CALL TOUPPER(KID3)
      IF      (KID .EQ. '#') THEN
C        Comments
C           All the following characters up to a new-line character
C           is ignored.
         STR = ' '
      ELSE IF (KID .EQ. 'R') THEN
         TIME = 0.D0
         TPRE = 0.D0
         IERR = 0
         ICONT = 1
         CALL TXPROF
         CALL TXEXEC
      ELSE IF (KID .EQ. 'C') THEN
         IF (ICONT .EQ. 0) THEN
            WRITE(6,*) 'XX RUN or LOAD before CONTINUE !'
            GOTO 10
         ENDIF
         NGR=-1
         CALL TXSTGR
         CALL TXEXEC
      ELSE IF (KID .EQ. 'P') THEN
         CALL TXPARM
      ELSE IF (KID .EQ. 'V') THEN
         CALL TXVIEW
      ELSE IF (KID .EQ. 'I') THEN
         CALL TXINIT
      ELSE IF (KID .EQ. 'E'.OR.KID .EQ. 'Q') THEN
         RETURN
      ELSE IF (KID .EQ. 'S') THEN
         CALL TXSAVE
      ELSE IF (KID .EQ. 'L') THEN
         CALL TXLOAD
         IERR = 0
         NGR = -1
         ICONT = 1
         CALL TXSTGR
      ELSE IF (KID .EQ. 'G') THEN
         KID=STR(2:2)
         CALL GUCPTL(KID)
C        GS
         IF      (KID .EQ. 'S') THEN
            CALL TXGSAV
C        GL
         ELSE IF (KID .EQ. 'L') THEN
            CALL TXGLOD
            IERR = 0
C        GI
         ELSE IF (KID .EQ. 'I') THEN
            TIME = 0.D0
            TPRE = 0.D0
            NGT = -1
            NGVV = -1
            CALL TXSTGT(SNGL(TIME))
            CALL TXSTGV(SNGL(TIME))
C        GT?
         ELSE IF (KID .EQ. 'T') THEN
            KID=STR(3:3)
            CALL GUCPTL(KID)
C           GTA
            IF      (KID .EQ. 'A') THEN
               DO NGPT = 1, NGPTM
                  MODE=MODEl
                  CALL TXGRFT(NGPT,MODE)
               ENDDO
C           GTB
            ELSE IF (KID .EQ. 'B') THEN
               DO NGPT = 1, 6
                  MODE=MODEl
                  CALL TXGRFT(NGPT,MODE)
               ENDDO
               DO NGPT = 9, 10
                  MODE=MODEl
                  CALL TXGRFT(NGPT,MODE)
               ENDDO
C           GTC
            ELSE IF (KID .EQ. 'C') THEN
               DO NGPT = 1, 5
                  MODE=MODEl
                  CALL TXGRFT(NGPT,MODE)
               ENDDO
C           GT1
            ELSE
               READ(STR(3:10),'(I8)',ERR=50) NGPT
               IF (NGPT.LT.1 .OR. NGPT.GT.NGPTM) THEN
                  IER = 1
               ELSE
                  MODE=MODEl
                  CALL TXGRFT(NGPT,MODE)
               ENDIF
            ENDIF
C        GV?
         ELSE IF (KID .EQ. 'V') THEN
            KID=STR(3:3)
            CALL GUCPTL(KID)
C           GTA
            IF      (KID .EQ. 'A') THEN
               DO NGPV = 1, NGPVM
                  MODE=MODEl
                  CALL TXGRFV(NGPV,MODE)
               ENDDO
C           GVB
            ELSE IF (KID .EQ. 'B') THEN
               DO NGPV = 1, 6
                  MODE=MODEl
                  CALL TXGRFV(NGPV,MODE)
               ENDDO
               DO NGPV = 9, 10
                  MODE=MODEl
                  CALL TXGRFV(NGPV,MODE)
               ENDDO
C           GVC
            ELSE IF (KID .EQ. 'C') THEN
               DO NGPV = 1, 5
                  MODE=MODEl
                  CALL TXGRFV(NGPV,MODE)
               ENDDO
C           GV1
            ELSE
               READ(STR(3:10),'(I8)',ERR=50) NGPV
               IF (NGPV.LT.1 .OR. NGPV.GT.NGPVM) THEN
                  IER = 1
               ELSE
                  MODE=MODEl
                  CALL TXGRFV(NGPV,MODE)
               ENDIF
            ENDIF
C        GU?
         ELSE IF (KID .EQ. 'U') THEN
            KID=STR(3:3)
            CALL GUCPTL(KID)
C           GUA
            IF      (KID .EQ. 'A') THEN
               DO NQ = 1, NQMAX, 4
                  CALL PAGES
                  DO NQL=NQ,MIN(NQ+3,NQMAX)
                     CALL TXGRFQ(NQL,MOD(NQL-1,4)+1)
                  ENDDO
                  CALL PAGEE
               ENDDO
C           GU1
            ELSE
               READ(STR(3:10),'(I8)',ERR=50) NQ
               IF (NQ.LT.1 .OR. NQ.GT.NQMAX) THEN
                  IER = 1
               ELSE
                  CALL PAGES
                  DO NQL=NQ,MIN(NQ+3,NQMAX)
                     CALL TXGRFQ(NQL,MOD(NQL-1,4)+1)
                  ENDDO
                  CALL PAGEE
               ENDIF
            ENDIF
C        GA
         ELSE IF (KID .EQ. 'A') THEN
            DO NGPR = 1, NGPRM
               MODE=MODEl
               CALL TXGRFR(NGPR,MODE)
            ENDDO
C        GB
         ELSE IF (KID .EQ. 'B') THEN
            DO NGPR = 1, 6
               MODE=MODEl
               CALL TXGRFR(NGPR,MODE)
            ENDDO
            DO NGPR = 9, 10
               MODE=MODEl
               CALL TXGRFR(NGPR,MODE)
            ENDDO
C        GC
         ELSE IF (KID .EQ. 'C') THEN
            DO NGPR = 1, 5
               MODE=MODEl
               CALL TXGRFR(NGPR,MODE)
            ENDDO
C        GM
         ELSE IF (KID .EQ. 'M') THEN
  100       WRITE(6,*) '## Number of files :'
            READ(5,*,END=10,ERR=100) NGFMAX
            NGR=-1
            DO NGF=1,NGFMAX
               CALL TXLOAD
               CALL TXSTGR
            ENDDO
  110       WRITE(6,*) '## INPUT GRAPH NUMBER'
            READ(5,'(A5)',ERR=110,END=10) STR2
            CALL TOUPPER(STR2)
            IF      (STR2 .EQ. '     ') THEN
C              For space or return only 
C           Correspond to GMA
            ELSE IF (STR2(1:1) .EQ. 'A') THEN
               DO I = 1, NGPRM
                  MODE=MODEl
                  CALL TXGRFR(I,MODE)
               ENDDO
C           Correspond to GMB
            ELSE IF (STR2(1:1) .EQ. 'B') THEN
               DO I = 1, 6
                  MODE=MODEl
                  CALL TXGRFR(I,MODE)
               ENDDO
               DO I = 9, 10
                  MODE=MODEl
                  CALL TXGRFR(I,MODE)
               ENDDO
C           Correspond to GMC
            ELSE IF (STR2(1:1) .EQ. 'C') THEN
               DO I = 1, 5
                  MODE=MODEl
                  CALL TXGRFR(I,MODE)
               ENDDO
            ELSE
               READ(STR2,'(I5)',ERR=110) NGPR
               IF      (NGPR .EQ. 0) THEN
                  GOTO 10
               ELSE IF (NGPR.GE.0 .AND. NGPR.LE.NGPRM) THEN
                  MODE=MODEl
                  CALL TXGRFR(NGPR,MODE)
               ENDIF
            ENDIF
            GOTO 110
C        G1
         ELSE
            READ(STR(2:10),'(I9)',ERR=50) NGPR
            IF (NGPR.LT.0 .OR. NGPR.GT.NGPRM) THEN
               IER = 1
            ELSE
               MODE=MODEl
               CALL TXGRFR(NGPR,MODE)
            ENDIF
         ENDIF
      ELSE IF (KID3 .EQ. 'DEL') THEN
         I = NINT((DelR / DR) - 0.5D0)
         X(I,1) = X(I,1) + DelN
         X(I,2) = X(I,2) + DelN
         CALL TXCALV(X)
      ELSE IF (KID .EQ. 'D') THEN
         KID=STR(2:2)
         CALL GUCPTL(KID)
C        DT?
         IF      (KID .EQ. 'T') THEN
            KID=STR(3:3)
            CALL GUCPTL(KID)
C           DT*?
            IF      (KID .EQ. '*') THEN
               READ(STR(4:10),'(I7)',ERR=50) NDT
               IF (NDT .GT. 0) THEN
                  DT = DT * NDT
                  WRITE(6,'(A,1PD9.2)') 'DT =', DT
               ELSE
                  IER = 1
               ENDIF
C           DT/?
            ELSE IF (KID .EQ. '/') THEN
               READ(STR(4:10),'(I7)',ERR=50) NDT
               IF (NDT .GT. 0) THEN
                  DT = DT / NDT
                  WRITE(6,'(A,1PD9.2)') 'DT =', DT
               ELSE
                  IER = 1
               ENDIF
C           DT1.d-3
            ELSE
               READ(STR(3:10),'(D9.2)',ERR=50) rNEWDT
               IF (rNEWDT .GT. 1.D-30) THEN
                  DT = rNEWDT
                  WRITE(6,'(A,1PD9.2)') 'DT =', DT
               ELSE
                  IER = 1
               ENDIF
            ENDIF
C        D
         ELSE
            CALL TXWDAT
            CALL TXWDAT2
         ENDIF
      ELSE IF (KID .EQ. 'B') THEN
         KID=STR(2:2)
         CALL GUCPTL(KID)
         IF      (KID .EQ. 'P') THEN
            PNBCD =  1.D0
         ELSE IF (KID .EQ. '0') THEN
            PNBCD =  0.D0
         ELSE IF (KID .EQ. 'M') THEN
            PNBCD = -1.D0
         ELSE
            IER = 1
         ENDIF
      ELSE
         IER = 1
      ENDIF
      IF (IER .NE. 0) GOTO 50
C
      NSTR = len(STR)
      DO I = 1, NSTR - 1
         IF (STR(I:I) .EQ. ',') THEN
            DO J = I + 1, NSTR
               IF (STR(J:J) .NE. ' ') THEN
                  STR(1:NSTR) = STR(J:NSTR)
                  GOTO 20
               ENDIF
            ENDDO
         ENDIF
      ENDDO
C
      GOTO 10
C
   40 CALL TXEXEC
      GOTO 10
C
   50 WRITE(6,*) '### ERROR : Invalid Command : ', STR
      GOTO 10
C
 9000 RETURN
      END
