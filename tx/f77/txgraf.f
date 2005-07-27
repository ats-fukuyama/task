!     $Id$
!
!     ***************************************************************
!
!        GRAPHIC Command loop
!
!     ***************************************************************
!
      SUBROUTINE TXGOUT

      INCLUDE 'txcomm.inc'

      INTEGER :: MODE, NGPR, NGPT, NGPV, NQ, NQL, NGF, NGFMAX, I
      CHARACTER(5) :: STR, STR2
      CHARACTER(1) :: KID1, KID2

!     *** MENU ***

   10 WRITE(6,*) '# SELECT : Rn: Tn: Un: Vn: A,B,C: ',
     &           'S,L,M:file I:init X:exit'
      READ(5,'(A5)',ERR=50,END=9000) STR

      KID1 = STR(1:1)
      KID2 = STR(2:2)
      CALL GUCPTL(KID1)
      CALL GUCPTL(KID2)

      IF      (KID1 .EQ. 'S') THEN
         CALL TXGSAV

      ELSE IF (KID1 .EQ. 'L') THEN
         CALL TXGLOD

      ELSE IF (KID1 .EQ. 'I') THEN
         TIME = 0.D0
         TPRE = 0.D0
         NGT = -1
         NGVV = -1
         CALL TXSTGT(SNGL(TIME))
         CALL TXSTGV(SNGL(TIME))

      ELSE IF (KID1.EQ.'R') THEN
         READ(STR(2:5),*,ERR=50) NGPR
         IF (NGPR.GE.0 .AND. NGPR.LE.NGPRM) THEN
            MODE=MODEl
            CALL TXGRFR(NGPR,MODE)
         ENDIF

      ELSE IF (KID1 .EQ. 'T') THEN
         IF      (KID2 .EQ. 'A') THEN
            DO NGPT = 1, NGPTM
               MODE=MODEl
               CALL TXGRFT(NGPT,MODE)
            ENDDO
         ELSE IF (KID2 .EQ. 'B') THEN
            DO NGPT = 1, 6
               MODE=MODEl
               CALL TXGRFT(NGPT,MODE)
            ENDDO
            DO NGPT = 9, 10
               MODE=MODEl
               CALL TXGRFT(NGPT,MODE)
            ENDDO

         ELSE IF (KID2 .EQ. 'C') THEN
            DO NGPT = 1, 5
               MODE=MODEl
               CALL TXGRFT(NGPT,MODE)
            ENDDO

         ELSE
            READ(STR(2:5),*,ERR=50) NGPT
            IF (NGPT.GE.1 .AND. NGPT.LE.NGPTM) THEN
               MODE=MODEl
               CALL TXGRFT(NGPT,MODE)
            ENDIF
         ENDIF

      ELSE IF (KID1 .EQ. 'V') THEN
         IF      (KID2 .EQ. 'A') THEN
            DO NGPV = 1, NGPVM
               MODE=MODEl
               CALL TXGRFV(NGPV,MODE)
            ENDDO

         ELSE IF (KID2 .EQ. 'B') THEN
            DO NGPV = 1, 6
               MODE=MODEl
               CALL TXGRFV(NGPV,MODE)
            ENDDO
            DO NGPV = 9, 10
               MODE=MODEl
               CALL TXGRFV(NGPV,MODE)
            ENDDO

         ELSE IF (KID2 .EQ. 'C') THEN
            DO NGPV = 1, 5
               MODE=MODEl
               CALL TXGRFV(NGPV,MODE)
            ENDDO

         ELSE
            READ(STR(2:5),*,ERR=50) NGPV
            IF (NGPV.GE.1 .AND. NGPV.LE.NGPVM) THEN
               MODE=MODEl
               CALL TXGRFV(NGPV,MODE)
            ENDIF
         ENDIF

      ELSE IF (KID1 .EQ. 'U') THEN
         IF (KID2 .EQ. 'A') THEN
            DO NQ = 1, NQMAX, 4
               CALL PAGES
               DO NQL=NQ,MIN(NQ+3,NQMAX)
                  CALL TXGRFQ(NQL,MOD(NQL-1,4)+1)
               ENDDO
               CALL PAGEE
            ENDDO

         ELSE
            READ(STR(2:5),*,ERR=50) NQ
            IF (NQ.GE.1 .AND. NQ.LE.NQMAX) THEN
               CALL PAGES
               DO NQL=NQ,MIN(NQ+3,NQMAX)
                  CALL TXGRFQ(NQL,NQL-NQ+1)
               ENDDO
               CALL PAGEE
            ENDIF
         ENDIF

      ELSE IF (KID1 .EQ. 'A') THEN
         DO NGPR = 1, NGPRM
            MODE=MODEl
            CALL TXGRFR(NGPR,MODE)
         ENDDO

      ELSE IF (KID1 .EQ. 'B') THEN
         DO NGPR = 1, 6
            MODE=MODEl
            CALL TXGRFR(NGPR,MODE)
         ENDDO
         DO NGPR = 9, 10
            MODE=MODEl
            CALL TXGRFR(NGPR,MODE)
         ENDDO

      ELSE IF (KID1 .EQ. 'C') THEN
         DO NGPR = 1, 5
            MODE=MODEl
            CALL TXGRFR(NGPR,MODE)
         ENDDO

      ELSE IF (KID1 .EQ. 'M') THEN
  100    WRITE(6,*) '## Number of files :'
         READ(5,*,END=10,ERR=100) NGFMAX
         NGR=-1
         DO NGF=1,NGFMAX
            CALL TXLOAD
            CALL TXSTGR
         ENDDO
  110    WRITE(6,*) '## INPUT GRAPH NUMBER'
         READ(5,'(A5)',ERR=110,END=10) STR2
         CALL TOUPPER(STR2)
         IF (STR2 .EQ. '     ') THEN
!     For space or return only 
!     Correspond to GMA
         ELSE IF (STR2(1:1) .EQ. 'A') THEN
            DO I = 1, NGPRM
               MODE=MODEl
               CALL TXGRFR(I,MODE)
            ENDDO
!     Correspond to GMB
         ELSE IF (STR2(1:1) .EQ. 'B') THEN
            DO I = 1, 6
               MODE=MODEl
               CALL TXGRFR(I,MODE)
            ENDDO
            DO I = 9, 10
               MODE=MODEl
               CALL TXGRFR(I,MODE)
            ENDDO
!     Correspond to GMC
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

      ELSE IF (KID1.EQ.'X') THEN
         GO TO 9000

      ELSE
         WRITE(6,*) 'XX UNKNOWN GRAPHIC COMMAND'
      ENDIF
      GOTO 10

   50 WRITE(6,*) '### ERROR : Invalid Command : ', STR
      GOTO 10

 9000 RETURN
      END
!
!     ***************************************************************
!
!        Initialize graphic axes
!
!     ***************************************************************
!
      SUBROUTINE TXPRFG

      INCLUDE 'txcomm.inc'
      INTEGER :: NR

      DO NR = 0, NRMAX
         GX(NR,1) = SNGL(R(NR) / RA)
      ENDDO

      DO NR = 0, NRMAX - 1
         GX(NR,2) = SNGL(RHI(NR) / RA)
      ENDDO
      GX(NRMAX,2) = SNGL(RB/RA)
      RETURN
      END
!
!     ***************************************************************
!
!        Store GY
!
!     ***************************************************************
!
      SUBROUTINE TXSTGR

      USE physical_constants, only : AEE
      INCLUDE 'txcomm.inc'
      
      INTEGER :: NR
      REAL(8) :: Bth, Bph, BBL
      
      IF (NGR.LT.NGRM) NGR = NGR + 1

      GT(NGR) = SNGL(TIME)

! Integer

      DO NR = 0, NRMAX
         GY(NR,NGR,3)  = SNGL(UerI(NR))
         GY(NR,NGR,4)  = SNGL(UethI(NR))
         GY(NR,NGR,6)  = SNGL(UirI(NR))
         GY(NR,NGR,7)  = SNGL(UithI(NR))
         GY(NR,NGR,9)  = SNGL(ErI(NR))
         GY(NR,NGR,10) = SNGL(BthI(NR))
         GY(NR,NGR,17) = SNGL(EthI(NR))
         GY(NR,NGR,19) = SNGL(UbthI(NR))
         GY(NR,NGR,20) = SNGL(Q(NR))
      ENDDO

! Half Integer

      DO NR = 0, NRMAX - 1
         GY(NR,NGR,1)  = SNGL(PNeHI(NR) * 1.D20)
         GY(NR,NGR,2)  = SNGL((PZ * PNiHI(NR) + PZ * PNbHI(NR)
     &                         -    PNeHI(NR)) * 1.D20)
         GY(NR,NGR,5)  = SNGL(UephHI(NR))
         GY(NR,NGR,8)  = SNGL(UiphHI(NR))
         GY(NR,NGR,11) = SNGL(EphHI(NR))
         GY(NR,NGR,12) = SNGL(PNbHI(NR) * 1.D20)
         GY(NR,NGR,13) = SNGL(UbphHI(NR))
         GY(NR,NGR,14) = SNGL(PTeHI(NR))
         GY(NR,NGR,15) = SNGL(PTiHI(NR))
         GY(NR,NGR,16) = SNGL((PN01HI(NR) + PN02HI(NR)) * 1.D20)
         GY(NR,NGR,18) = SNGL(BphHI(NR))
         GY(NR,NGR,21) = SNGL(-      AEE * PNeHI(NR) * UephHI(NR)*1.D20
     &                        + PZ * AEE * PNiHI(NR) * UiphHI(NR)*1.D20
     &                        + PZ * AEE * PNbHI(NR) * UbphHI(NR)*1.D20)
         GY(NR,NGR,22) = SNGL(-      AEE * PNeHI(NR) * UephHI(NR)*1.D20)
         GY(NR,NGR,23) = SNGL(  PZ * AEE * PNiHI(NR) * UiphHI(NR)*1.D20)
         GY(NR,NGR,24) = SNGL(  PZ * AEE * PNbHI(NR) * UbphHI(NR)*1.D20)
         GY(NR,NGR,35) = SNGL(PN01HI(NR) * 1.D20)
         GY(NR,NGR,36) = SNGL(PN02HI(NR) * 1.D20)
      ENDDO

! Integer

      DO NR = 0, NRMAX
         Bth = BthI(NR)
         Bph = BphI(NR)
         BBL = SQRT(Bph**2 + Bth*2)
         GY(NR,NGR,25) = SNGL(+ Bph * UethI(NR) / BBL
     &                        - Bth * UephI(NR) / BBL)
         GY(NR,NGR,26) = SNGL(+ Bth * UethI(NR) / BBL
     &                        + Bph * UephI(NR) / BBL)
         GY(NR,NGR,27) = SNGL(+ Bph * UithI(NR) / BBL
     &                        - Bth * UiphI(NR) / BBL)
         GY(NR,NGR,28) = SNGL(+ Bth * UithI(NR) / BBL
     &                        + Bph * UiphI(NR) / BBL)
      ENDDO

! Half Integer

      DO NR = 0, NRMAX - 1
         GY(NR,NGR,29) = SNGL(Di(NR)+De(NR))
         IF (rG1h2(NR).GT.3.D0) THEN
            GY(NR,NGR,30) = 3.0
         ELSE
            GY(NR,NGR,30) = SNGL(rG1h2(NR))
         ENDIF
         GY(NR,NGR,31) = SNGL(FCDBM(NR))
         GY(NR,NGR,32) = SNGL(S(NR))
         GY(NR,NGR,33) = SNGL(Alpha(NR))
         GY(NR,NGR,34) = SNGL(rKappa(NR))
      ENDDO

      RETURN
      END
!
!     ***************************************************************
!
!        Store GTY
!
!     ***************************************************************
!
      SUBROUTINE TXSTGT(GTIME)

      USE physical_constants, only : AEE, 
     &     Phys_Constants_Initialization, PI
      INCLUDE 'txcomm.inc'

      REAL, INTENT(IN) :: GTIME
      INTEGER :: NR
      REAL(8) :: rLINEAVE, Bth, Bph, BBL, SUM1 = 0.D0, SUM2 = 0.D0, DELS

      CALL Phys_Constants_Initialization

      IF (NGT.LT.NGTM) NGT=NGT+1

      GTX(NGT) = GTIME

      GTY(NGT,1) = SNGL(PNeI(0) * 1.D20)
      GTY(NGT,2) = SNGL(rLINEAVE(0.D0))
      GTY(NGT,3) = SNGL(rLINEAVE(0.24D0))
      GTY(NGT,4) = SNGL(rLINEAVE(0.6D0))
      GTY(NGT,5) = SNGL((PZ * PNiHI(0) + PZ * PNbHI(0)
     &                                 - PNeHI(0)) * 1.D20)
      GTY(NGT,6)  = SNGL(UerI(NRMAX/2))
      GTY(NGT,7)  = SNGL(UethI(NRMAX/2))
      GTY(NGT,8)  = SNGL(UephI(0))
      GTY(NGT,9)  = SNGL(UirI(NRMAX/2))
      GTY(NGT,10) = SNGL(UithI(NRMAX/2))
      GTY(NGT,11) = SNGL(UiphI(0))
      GTY(NGT,12) = SNGL(ErI(NRMAX/2))
      GTY(NGT,13) = SNGL(BthI(NRMAX/2))
      GTY(NGT,14) = SNGL(EphHI(NRMAX-1))
      GTY(NGT,15) = SNGL(PNbHI(0) * 1.D20)
      GTY(NGT,16) = SNGL(UbphI(0))
      GTY(NGT,17) = SNGL(PTeI(0))
      GTY(NGT,18) = SNGL(PTiI(0))
      GTY(NGT,19) = SNGL((PN01HI(0) + PN02HI(0)) * 1.D20)

      GTY(NGT,20) = SNGL(Q(0))
      GTY(NGT,21) = SNGL((-      AEE * PNeI(0) * UephI(0)
     &                    + PZ * AEE * PNiI(0) * UiphI(0)
     &                    + PZ * AEE * PNbI(0) * UbphI(0)) * 1.D20)
      GTY(NGT,22) = SNGL(-      AEE * PNeI(0) * UephI(0) * 1.D20)
      GTY(NGT,23) = SNGL(  PZ * AEE * PNiI(0) * UiphI(0) * 1.D20)
      GTY(NGT,24) = SNGL(  PZ * AEE * PNbI(0) * UbphI(0) * 1.D20)

      Bth = BthHI(NRMAX/2)
      Bph = BphHI(NRMAX/2)
      BBL = SQRT(Bph**2 + Bth**2)
      GTY(NGT,25) = SNGL((+ Bph * UethI(NRMAX/2)
     &                    - Bth * UephI(NRMAX/2)) / BBL)
      GTY(NGT,26) = SNGL((+ Bth * UethI(NRMAX/2)
     &                    + Bph * UephI(NRMAX/2)) / BBL)
      GTY(NGT,27) = SNGL((+ Bph * UithI(NRMAX/2)
     &                    - Bth * UiphI(NRMAX/2)) / BBL)
      GTY(NGT,28) = SNGL((+ Bth * UithI(NRMAX/2)
     &                    + Bph * UiphI(NRMAX/2)) / BBL)
      GTY(NGT,29) = SNGL(Di(NRMAX/2)+De(NRMAX/2))
      GTY(NGT,30) = SNGL(rG1h2(NRMAX/2))
      GTY(NGT,31) = SNGL(FCDBM(NRMAX/2))
      GTY(NGT,32) = SNGL(S(NRMAX/2))
      GTY(NGT,33) = SNGL(Alpha(NRMAX/2))
      GTY(NGT,34) = SNGL(rKappa(NRMAX/2))

!      GTY(NGT,19) = SNGL((PN01HI(0) + PN02HI(0)) * 1.D20)
      GTY(NGT,38) = SNGL(PN01HI(0) * 1.D20)
      GTY(NGT,39) = SNGL(PN02HI(0) * 1.D20)

      DO NR=0,NRMAX-1
         IF(R(NR).LE.RA) THEN
            IF(R(NR+1).LE.RA) THEN
               DELS=2.D0*PI*RR*2.D0*PI*RHI(NR)*DR
            ELSE
               DELS=2.D0*PI*RR*2.D0*PI*RHI(NR)*(RA-R(NR))
            ENDIF
            SUM1=SUM1+          PNeHI(NR)*DELS
            SUM2=SUM2+rNuION(NR)*PNeHI(NR)*DELS
         ENDIF
      ENDDO
      GTY(NGT,35) = SNGL(SUM1)
      GTY(NGT,36) = SNGL(SUM2)
      IF(NGT.EQ.0.OR.ABS(SUM2).LE.0.D0) THEN
         GTY(NGT,37) = 0.0
      ELSE
         GTY(NGT,37) = SNGL(SUM1/SUM2)
      ENDIF

      RETURN
      END
!
!     ***************************************************************
!
!        Store GVY
!
!     ***************************************************************
!
      SUBROUTINE TXSTGV(GTIME)

      INCLUDE 'txcomm.inc'

      REAL, INTENT(IN) :: GTIME

      IF (NGVV.LT.NGVM) NGVV=NGVV+1

      GVX(NGVV) = GTIME

      GVY(NGVV,1) = SNGL(TS0(1))
      GVY(NGVV,2) = SNGL(TS0(2))
      GVY(NGVV,3) = SNGL(TSAV(1))
      GVY(NGVV,4) = SNGL(TSAV(2))
      GVY(NGVV,5) = SNGL(PINT)
      GVY(NGVV,6) = SNGL(POHT)
      GVY(NGVV,7)  = SNGL(PNBT)
      GVY(NGVV,8)  = SNGL(PRFT)

      GVY(NGVV,10) = SNGL(AJT)
      GVY(NGVV,11) = SNGL(AJOHT)
      GVY(NGVV,12) = SNGL(AJNBT)
      GVY(NGVV,13) = SNGL(AJOHT+AJNBT)
!      GVY(NGVV,14) = SNGL(AJBST)
      GVY(NGVV,15) = SNGL(POUT)
      GVY(NGVV,16) = SNGL(PCXT)
      GVY(NGVV,17) = SNGL(PIET)
!      GVY(NGVV,18) = SNGL(PRLT)
!      GVY(NGVV,19) = SNGL(PCONT)
      GVY(NGVV,20) = SNGL(QF)
      GVY(NGVV,21) = SNGL(TS0(1))
      GVY(NGVV,22) = SNGL(TS0(2))
      GVY(NGVV,23) = SNGL(TSAV(1))
      GVY(NGVV,24) = SNGL(TSAV(2))
      GVY(NGVV,25) = SNGL(ANS0(1))
      GVY(NGVV,26) = SNGL(ANS0(2))
      GVY(NGVV,27) = SNGL(ANSAV(1))
      GVY(NGVV,28) = SNGL(ANSAV(2))
      GVY(NGVV,29) = SNGL(WST(1))
      GVY(NGVV,30) = SNGL(WST(2))

      GVY(NGVV,31) = SNGL(WST(1))
      GVY(NGVV,32) = SNGL(WST(2))

      GVY(NGVV,33) = SNGL(TAUE1)
      GVY(NGVV,34) = SNGL(TAUE2)
      GVY(NGVV,35) = SNGL(TAUEP)
      GVY(NGVV,36) = SNGL(BETAA)
      GVY(NGVV,37) = SNGL(BETA0)
      GVY(NGVV,38) = SNGL(BETAPA)
      GVY(NGVV,39) = SNGL(BETAP0)
      GVY(NGVV,40) = SNGL(VLOOP)
      GVY(NGVV,41) = SNGL(ALI)
      GVY(NGVV,42) = SNGL(Q(0))
      GVY(NGVV,43) = SNGL(RQ1)

      RETURN
      END
!
!     ***************************************************************
!
!        Time evolution of radial profile
!
!     ***************************************************************
!
      SUBROUTINE TXGRFR(NGYR,MODE)

      INCLUDE 'txcomm.inc'

      INTEGER, INTENT(IN) :: NGYR, MODE
      INTEGER :: IND, nr
      REAL :: gDIVL
      REAL(8), DIMENSION(0:NRM,0:NGRM) :: GYL
      CHARACTER(40) :: STR

      IF (NGYR. EQ. 0) THEN
         CALL PAGES
         CALL PAGEE
         RETURN
      ENDIF

      IF (NGR .LE. -1) THEN
         WRITE(6,*) 'G', NGYR, ' has no data'
         RETURN
      ENDIF

      IF (MODEG .EQ. 2) THEN
         IND = 9
      ELSE
         IND = 0
      ENDIF

      CALL PAGES
      CALL SETCHS(0.3, 0.0)
      CALL SETLIN(0, 1, 7)

      CALL MOVE(2.0,17.7)
      CALL TEXT('[G', 2)
      CALL NUMBI(NGYR, '(I2)', 2)
      CALL TEXT(']  ', 3)
      CALL TEXT('FROM', 4)
      CALL NUMBR(GT(0), '(1PE9.2)', 9)
      CALL TEXT(' TO', 3)
      CALL NUMBR(GT(NGR), '(1PE9.2)', 9)
      CALL TEXT('  DT =', 6)
      CALL NUMBD(DT, '(1PD9.2)', 9)
      CALL TEXT('  NGRSTP = ', 11)
      CALL NUMBI(NGRSTP,'(I4)',4)

      IF (NGYR .EQ. 1) THEN
         STR = '@n$-e$=(r)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(1)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GY(0,0,1), GYL, NRM+1, NRMAX, NGR+1, gDIVL)
         CALL TXGRFRX(0, GX(0,2), GYL, NRMAX-1, NGR, STR, MODE, IND)

         STR = '@Z*n$-i$=+Z*n$-b$=-n$-e$=@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(2)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GY(0,0,2), GYL, NRM+1, NRMAX, NGR+1, gDIVL)
         CALL TXGRFRX(1, GX(0,2), GYL, NRMAX-1, NGR, STR, MODE, IND)

         STR = '@E$-r$=(r)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(9)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GY(0,0,9), GYL, NRM+1, NRMAX+1, NGR+1, gDIVL)
         CALL TXGRFRX(2, GX(0,1), GYL, NRMAX, NGR, STR, MODE, IND)

         CALL TXWPGR

      ELSEIF (NGYR .EQ. 2) THEN
         STR = '@u$-er$=(r)@'
         CALL TXGRFRX(0, GX(0,1), GY(0,0,3), NRMAX, NGR, STR, MODE, IND)

         STR = '@u$-e$#q$#$=(r)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(4)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GY(0,0,4), GYL, NRM+1, NRMAX+1, NGR+1, gDIVL)
         CALL TXGRFRX(1, GX(0,1), GYL, NRMAX, NGR, STR, MODE, IND)

         STR = '@u$-e$#f$#$=(r)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(5)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GY(0,0,5), GYL, NRM+1, NRMAX, NGR+1, gDIVL)
         CALL TXGRFRX(2, GX(0,2), GYL, NRMAX-1, NGR, STR, MODE, IND)

         CALL TXWPGR

      ELSEIF (NGYR .EQ. 3) THEN
         STR = '@u$-ir$=(r)@'
         STR = '@Uir(r)@'
         CALL TXGRFRX(0, GX(0,1), GY(0,0,6), NRMAX, NGR, STR, MODE, IND)

         STR = '@u$-i$#q$#$=(r)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(7)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GY(0,0,7), GYL, NRM+1, NRMAX+1, NGR+1, gDIVL)
         CALL TXGRFRX(1, GX(0,1), GYL, NRMAX, NGR, STR, MODE, IND)

         STR = '@u$-i$#f$#$=(r)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(8)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GY(0,0,8), GYL, NRM+1, NRMAX, NGR+1, gDIVL)
         CALL TXGRFRX(2, GX(0,2), GYL, NRMAX-1, NGR, STR, MODE, IND)

         CALL TXWPGR

      ELSEIF (NGYR .EQ. 4) THEN
         STR = '@q(r)@'
        CALL TXGRFRX(0,GX(0,1),GY(0,0,20),NRMAX  ,NGR,STR,MODE,IND)

         STR = '@E$-$#f$#$=(r)@'
        CALL TXGRFRX(1,GX(0,2),GY(0,0,11),NRMAX-1,NGR,STR,MODE,IND)

         STR = '@j$-$#f$#$=(r)@'
        CALL TXGRFRX(2,GX(0,2),GY(0,0,21),NRMAX-1,NGR,STR,MODE,IND)

         CALL TXWPGR

      ELSEIF (NGYR .EQ. 5) THEN
         STR = '@n$-b$=(r)@'
        CALL TXGRFRX(0,GX(0,2),GY(0,0,12),NRMAX-1,NGR,STR,MODE,IND)
         STR = '@u$-b$#f$#$=(r)@'
        CALL TXGRFRX(1,GX(0,1),GY(0,0,13),NRMAX  ,NGR,STR,MODE,IND)
         STR = '@B$-$#q$#$=(r)@'
        CALL TXGRFRX(2,GX(0,1),GY(0,0,10),NRMAX  ,NGR,STR,MODE,IND)
         CALL TXWPGR

      ELSEIF (NGYR .EQ. 6) THEN
         STR = '@j$-e$#f$#$=(r)@'
        CALL TXGRFRX(0,GX(0,1),GY(0,0,22),NRMAX  ,NGR,STR,MODE,IND)
         STR = '@j$-i$#f$#$=(r)@'
        CALL TXGRFRX(1,GX(0,1),GY(0,0,23),NRMAX  ,NGR,STR,MODE,IND)
         STR = '@j$-b$#f$#$=(r)@'
        CALL TXGRFRX(2,GX(0,1),GY(0,0,24),NRMAX  ,NGR,STR,MODE,IND)
         CALL TXWPGR

      ELSEIF (NGYR .EQ. 7) THEN
         STR = '@u$-er$=(r)@'
        CALL TXGRFRX(0,GX(0,1),GY(0,0,3),NRMAX  ,NGR,STR,MODE,IND)
         STR = '@u$-e$#$/136$#$=(r)@'
        CALL TXGRFRX(1,GX(0,1),GY(0,0,25),NRMAX  ,NGR,STR,MODE,IND)
         STR = '@u$-e//$=(r)@'
        CALL TXGRFRX(2,GX(0,1),GY(0,0,26),NRMAX  ,NGR,STR,MODE,IND)
         CALL TXWPGR

      ELSEIF (NGYR .EQ. 8) THEN
         STR = '@u$-ir$=(r)@'
        CALL TXGRFRX(0,GX(0,1),GY(0,0,6),NRMAX  ,NGR,STR,MODE,IND)
         STR = '@u$-i$#$/136$#$=(r)@'
        CALL TXGRFRX(1,GX(0,1),GY(0,0,27),NRMAX  ,NGR,STR,MODE,IND)
         STR = '@u$-i//$=(r)@'
        CALL TXGRFRX(2,GX(0,1),GY(0,0,28),NRMAX  ,NGR,STR,MODE,IND)
         CALL TXWPGR

      ELSEIF (NGYR .EQ. 9) THEN
         STR = '@D$-i$=(r)+D$-e$=(r)@'
        CALL TXGRFRX(0,GX(0,2),GY(0,0,29),NRMAX-1,NGR,STR,MODE,IND)
         STR = '@G$-1$=h$+2$=(r)@'
        CALL TXGRFRX(1,GX(0,2),GY(0,0,30),NRMAX-1,NGR,STR,MODE,IND)
         STR = '@F$-CDBM$=(r)@'
        CALL TXGRFRX(2,GX(0,2),GY(0,0,31),NRMAX-1,NGR,STR,MODE,IND)
         CALL TXWPGR

      ELSEIF (NGYR .EQ. 10) THEN
         STR = '@T$-e$=(r)@'
        CALL TXGRFRX(0,GX(0,2),GY(0,0,14),NRMAX-1,NGR,STR,MODE,IND)

         STR = '@T$-i$=(r)@'
        CALL TXGRFRX(1,GX(0,2),GY(0,0,15),NRMAX-1,NGR,STR,MODE,IND)

         STR = '@N$-0$=(r)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(16)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR,gDIVL)
         CALL DIVGY(GY(0,0,16),GYL,NRM+1,NRMAX,NGR+1,gDIVL)
         CALL TXGRFRX(2,GX(0,2),GYL,NRMAX-1,NGR,STR,MODE,IND)

         CALL TXWPGR

      ELSEIF (NGYR .EQ. 11) THEN
         STR = '@S(r)@'
        CALL TXGRFRX(0,GX(0,2),GY(0,0,32),NRMAX-1,NGR,STR,MODE,IND)
         STR = '@$#a$#(r)@'
        CALL TXGRFRX(1,GX(0,2),GY(0,0,33),NRMAX-1,NGR,STR,MODE,IND)
         STR = '@$#k$#(r)@'
        CALL TXGRFRX(2,GX(0,2),GY(0,0,34),NRMAX-1,NGR,STR,MODE,IND)
         CALL TXWPGR

      ELSEIF (NGYR .EQ. 12) THEN
         STR = '@E$-$#q$#$=(r)@'
        CALL TXGRFRX(0,GX(0,1),GY(0,0,17),NRMAX  ,NGR,STR,MODE,IND)
         STR = '@B$-$#f$#$=(r)@'
        CALL TXGRFRX(1,GX(0,2),GY(0,0,18),NRMAX-1,NGR,STR,MODE,IND)
         STR = '@U$-b$#q$#$=(r)@'
        CALL TXGRFRX(2,GX(0,1),GY(0,0,19),NRMAX  ,NGR,STR,MODE,IND)
         CALL TXWPGR

      ELSEIF (NGYR .EQ. 13) THEN
         STR = '@SLOW N$-0$=(r)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(16)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR,gDIVL)
         CALL DIVGY(GY(0,0,35),GYL,NRM+1,NRMAX,NGR+1,gDIVL)
         CALL TXGRFRX(0,GX(0,2),GYL,NRMAX-1,NGR,STR,MODE,IND)

         STR = '@FAST N$-0$+(r)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(16)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR,gDIVL)
         CALL DIVGY(GY(0,0,36),GYL,NRM+1,NRMAX,NGR+1,gDIVL)
         CALL TXGRFRX(1,GX(0,2),GYL,NRMAX-1,NGR,STR,MODE,IND)

         STR = '@TOTAL N$-0$=(r)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(16)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR,gDIVL)
         CALL DIVGY(GY(0,0,16),GYL,NRM+1,NRMAX,NGR+1,gDIVL)
         CALL TXGRFRX(2,GX(0,2),GYL,NRMAX-1,NGR,STR,MODE,IND)
         CALL TXWPGR

      ELSE
         WRITE(6,*) 'Unknown NGYR: NGYR = ',NGYR
      ENDIF

      CALL PAGEE

      RETURN
      END
!
!     ***************************************************************
!
!        SLAVE ROUTINE TO Write graph of GX, GY
!
!     ***************************************************************
!
      SUBROUTINE TXGRFRX(K, GXL, GYL, NRMAX1, NGMAX, STR, MODE, IND)

      INCLUDE 'txcomm.inc'

      INTEGER, INTENT(IN) :: K, NRMAX1, NGMAX, MODE, IND
      REAL, DIMENSION(0:NRM), INTENT(IN) :: GXL
      REAL, DIMENSION(0:NRM,0:NGMAX), INTENT(IN) :: GYL
      CHARACTER(*), INTENT(IN) :: STR
      REAL :: GXMAX
      REAL, DIMENSION(4) :: GPXY

      GPXY(1) =  3.0 + 12.5 * MOD(K,2)
      GPXY(2) = 12.5 + 12.5 * MOD(K,2)
      GPXY(3) = 10.5 -  8.5 * REAL(K/2)
      GPXY(4) = 17.0 -  8.5 * REAL(K/2)
      GXMAX=REAL(RB/RA)
      CALL TXGRAF(GPXY, GXL, GYL, NRM+1, NRMAX1+1, NGMAX+1,
     &            0.0, GXMAX, STR, MODE, IND)

      RETURN
      END
!
!     ***************************************************************
!
!        Write graph of GTX, GTY
!
!     ***************************************************************
!
      SUBROUTINE TXGRFT(NGYT,MODE)

      INCLUDE 'txcomm.inc'
      
      INTEGER, INTENT(IN) :: NGYT, MODE
      INTEGER :: IND
      REAL :: gDIVL
      CHARACTER(40) ::  STR
      REAL, DIMENSION(0:NGTM,NGYTM) :: GTYL

      IF (NGT .LE. 1) THEN
         WRITE(6,*) 'GT', NGYT, ' has no data'
         RETURN
      ENDIF

      IF (MODEG .EQ. 2) THEN
         IND = 9
      ELSE
         IND = 0
      ENDIF

      CALL PAGES
      CALL SETCHS(0.3, 0.0)
      CALL SETLIN(0, 1, 7)

      CALL MOVE(2.0,17.7)
      CALL TEXT('[GT', 3)
      CALL NUMBI(NGYT, '(I2)', 2)
      CALL TEXT(']  ', 3)
      CALL TEXT('FROM', 4)
      CALL NUMBR(GTX(0), '(1PE9.2)', 9)
      CALL TEXT(' TO', 3)
      CALL NUMBR(GTX(NGT), '(1PE9.2)', 9)
      CALL TEXT('  DT =', 6)
      CALL NUMBD(DT, '(1PD9.2)', 9)
      CALL TEXT('  NGTSTP = ', 11)
      CALL NUMBI(NGTSTP,'(I4)',4)

      IF (NGYT .EQ. 1) THEN
         STR = '@Ne(0),LAVE(0),(0.24),(0.6)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(1)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GTY(0,1), GTYL, NGTM+1, NGT+1, 4, gDIVL)
         CALL TXGRFTX(0, GTX, GTYL, NGTM, NGT, 4, STR,MODE, IND)

         STR = '@Z*Ni+Z*Nb-Ne(0)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(5-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GTY(0,5), GTYL, NGTM+1, NGT+1, 1, gDIVL)
         CALL TXGRFTX(1, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

         STR = '@Er(b/2)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(12-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GTY(0,12), GTYL, NGTM+1, NGT+1, 1, gDIVL)
         CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

         CALL TXWPGR

      ELSEIF (NGYT .EQ. 2) THEN
         STR = '@Uer(b/2)@'
         CALL TXGRFTX(0, GTX, GTY(0,6), NGTM, NGT, 1, STR,MODE, IND)

         STR = '@UeTheta(b/2)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(7-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GTY(0,7), GTYL, NGTM+1, NGT+1, 1, gDIVL)
         CALL TXGRFTX(1, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

         STR = '@UePhi(0)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(8-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GTY(0,8), GTYL, NGTM+1, NGT+1, 1, gDIVL)
         CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

         CALL TXWPGR

      ELSEIF (NGYT .EQ. 3) THEN
         STR = '@Uir(b/2)@'
         CALL TXGRFTX(0, GTX, GTY(0,9), NGTM, NGT, 1, STR,MODE, IND)

         STR = '@UiTheta(b/2)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(10-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GTY(0,10), GTYL, NGTM+1, NGT+1, 1, gDIVL)
         CALL TXGRFTX(1, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

         STR = '@UiPhi(0)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(11-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GTY(0,11), GTYL, NGTM+1, NGT+1, 1, gDIVL)
         CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

         CALL TXWPGR

      ELSEIF (NGYT .EQ. 4) THEN
         STR = '@Q(0)@'
         CALL TXGRFTX(0, GTX, GTY(0,20), NGTM, NGT, 1, STR,MODE, IND)

         STR = '@Ephi(b)@'
         CALL TXGRFTX(1, GTX, GTY(0,14), NGTM, NGT, 1, STR,MODE, IND)

         STR = '@Jphi(0)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(21-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GTY(0,21), GTYL, NGTM+1, NGT+1, 1, gDIVL)
         CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

         CALL TXWPGR

      ELSEIF (NGYT .EQ. 5) THEN
         STR = '@Nb(0)@'
         CALL TXGRFTX(0, GTX, GTY(0,15), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@UbPhi(0)@'
         CALL TXGRFTX(1, GTX, GTY(0,16), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@Bth(b/2)@'
         CALL TXGRFTX(2, GTX, GTY(0,13), NGTM, NGT, 1, STR,MODE, IND)
         CALL TXWPGR

      ELSEIF (NGYT .EQ. 6) THEN
         STR = '@JePhi(0)@'
         CALL TXGRFTX(0, GTX, GTY(0,22), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@JiPhi(0)@'
         CALL TXGRFTX(1, GTX, GTY(0,23), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@JbPhi(0)@'
         CALL TXGRFTX(2, GTX, GTY(0,24), NGTM, NGT, 1, STR,MODE, IND)
         CALL TXWPGR

      ELSEIF (NGYT .EQ. 7) THEN
         STR = '@Uer(b/2)@'
         CALL TXGRFTX(0, GTX, GTY(0,6), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@UePerp(b/2)@'
         CALL TXGRFTX(1, GTX, GTY(0,25), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@UePara(b/2)@'
         CALL TXGRFTX(2, GTX, GTY(0,26), NGTM, NGT, 1, STR,MODE, IND)
         CALL TXWPGR

      ELSEIF (NGYT .EQ. 8) THEN
         STR = '@Uir(b/2)@'
         CALL TXGRFTX(0, GTX, GTY(0,9), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@UiPerp(b/2)@'
         CALL TXGRFTX(1, GTX, GTY(0,27), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@UiPara(b/2)@'
         CALL TXGRFTX(2, GTX, GTY(0,28), NGTM, NGT, 1, STR,MODE, IND)
         CALL TXWPGR

      ELSEIF (NGYT .EQ. 9) THEN
         STR = '@Di(b/2)+De(b/2)@'
         CALL TXGRFTX(0, GTX, GTY(0,29), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@rG1h2(b/2)@'
         CALL TXGRFTX(1, GTX, GTY(0,30), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@FCDBM(b/2)@'
         CALL TXGRFTX(2, GTX, GTY(0,31), NGTM, NGT, 1, STR,MODE, IND)
         CALL TXWPGR

      ELSEIF (NGYT .EQ. 10) THEN
         STR = '@Te(0)@'
         CALL TXGRFTX(0, GTX, GTY(0,17), NGTM, NGT, 1, STR,MODE, IND)

         STR = '@Ti(0)@'
         CALL TXGRFTX(1, GTX, GTY(0,18), NGTM, NGT, 1, STR,MODE, IND)

         STR = '@N0(0)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(19-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GTY(0,19), GTYL, NGTM+1, NGT+1, 1, gDIVL)
         CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

         CALL TXWPGR

      ELSEIF (NGYT .EQ. 11) THEN
         STR = '@S(b/2)@'
         CALL TXGRFTX(0, GTX, GTY(0,32), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@Alpha(b/2)@'
         CALL TXGRFTX(1, GTX, GTY(0,33), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@rKappa(b/2)@'
         CALL TXGRFTX(2, GTX, GTY(0,34), NGTM, NGT, 1, STR,MODE, IND)
         CALL TXWPGR

         CALL TXWPGR

      ELSEIF (NGYT .EQ. 12) THEN
         STR = '@Ne:total@'
         CALL TXGRFTX(0, GTX, GTY(0,35), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@Gamma:total@'
         CALL TXGRFTX(1, GTX, GTY(0,36), NGTM, NGT, 1, STR,MODE, IND)
         STR = '@taup@'
         CALL TXGRFTX(2, GTX, GTY(0,37), NGTM, NGT, 1, STR,MODE, IND)
         CALL TXWPGR

      ELSEIF (NGYT .EQ. 13) THEN

         STR = '@SLOW N0(0)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(19-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GTY(0,38), GTYL, NGTM+1, NGT+1, 1, gDIVL)
         CALL TXGRFTX(0, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

         STR = '@FAST N0(0)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(19-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GTY(0,39), GTYL, NGTM+1, NGT+1, 1, gDIVL)
         CALL TXGRFTX(1, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

         STR = '@TOTAL N0(0)@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(19-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GTY(0,19), GTYL, NGTM+1, NGT+1, 1, gDIVL)
         CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

         CALL TXWPGR
      ELSE
         WRITE(6,*) 'Unknown NGYT: NGYT = ',NGYT
      ENDIF

      CALL PAGEE

      RETURN
      END
!
!     ***************************************************************
!
!        Write graph of GVX, GVY
!
!     ***************************************************************
!
      SUBROUTINE TXGRFV(NGYV,MODE)

      INCLUDE 'txcomm.inc'

      INTEGER, INTENT(IN) :: NGYV, MODE
      INTEGER :: IND
      REAL :: gDIVL
      CHARACTER(40) ::  STR
      REAL, DIMENSION(0:NGVM,NGYVM) :: GVYL

      IF (NGVV .LE. 1) THEN
         WRITE(6,*) 'GV', NGYV, ' has no data'
         RETURN
      ENDIF

      IF (MODEG .EQ. 2) THEN
         IND = 9
      ELSE
         IND = 0
      ENDIF

      CALL PAGES
      CALL SETCHS(0.3, 0.0)
      CALL SETLIN(0, 1, 7)

      CALL MOVE(2.0,17.7)
      CALL TEXT('[GV', 3)
      CALL NUMBI(NGYV, '(I2)', 2)
      CALL TEXT(']  ', 3)
      CALL TEXT('FROM', 4)
      CALL NUMBR(GVX(0), '(1PE9.2)', 9)
      CALL TEXT(' TO', 3)
      CALL NUMBR(GVX(NGVV), '(1PE9.2)', 9)
      CALL TEXT('  DT =', 6)
      CALL NUMBD(DT, '(1PD9.2)', 9)
      CALL TEXT('  NGVSTP = ', 11)
      CALL NUMBI(NGVSTP,'(I4)',4)

      IF (NGYV .EQ. 1) THEN
         STR = '@Te,Ti,<Te>,<Ti> [keV] vs t@'
         gDIVL=1.E0
         CALL DIVGY(GVY(0,1), GVYL, NGVM+1, NGVV+1, 4, gDIVL)
         CALL TXGRFTX(0, GVX, GVYL, NGVM, NGVV, 4, STR, MODE, IND)

         STR = '@PIN,POH,PNB,PRF,PNF [MW] vs t@'
            gDIVL = 1.E0
         CALL DIVGY(GVY(0,5), GVYL, NGVM+1, NGVV+1, 5, gDIVL)
         CALL TXGRFTX(1, GVX, GVYL, NGVM, NGVV, 5, STR, MODE, IND)

         STR = '@IP,IOH,INB,IOH+INB [MA] vs t@'
            gDIVL = 1.E0
         CALL DIVGY(GVY(0,10), GVYL, NGVM+1, NGVV+1, 4, gDIVL)
         CALL TXGRFTX(2, GVX, GVYL, NGVM, NGVV, 4, STR, MODE, IND)

         STR = '@POUT,PCX,PIE,PRL,PCON [MW] vs t@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(12-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL DIVGY(GVY(0,15), GVYL, NGVM+1, NGVV+1, 5, gDIVL)
         CALL TXGRFTX(3, GVX, GVYL, NGVM, NGVV, 5, STR, MODE, IND)

      ELSEIF (NGYV .EQ. 2) THEN
         STR = '@QF vs t@'
         CALL TXGRFTX(0, GVX, GVY(0,20), NGVM, NGVV, 1, STR, MODE, IND)

         STR = '@Te0,Ti0 [keV] vs t@'
            gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,21), GVYL, NGVM+1, NGVV+1, 2, gDIVL)
         CALL TXGRFTX(1, GVX, GVYL, NGVM, NGVV, 2, STR, MODE, IND)

         STR = '@<Te>,<Ti> [keV] vs t@'
             gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,23), GVYL, NGVM+1, NGVV+1, 2, gDIVL)
         CALL TXGRFTX(2, GVX, GVYL, NGVM, NGVV, 2, STR, MODE, IND)

         STR = '@Ne0,Ni0,<Ne>,<Ni> [10^20/m^3] vs t@'
            gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,25), GVYL, NGVM+1, NGVV+1, 4, gDIVL)
         CALL TXGRFTX(3, GVX, GVYL, NGVM, NGVV, 4, STR, MODE, IND)

      ELSEIF (NGYV .EQ. 3) THEN
         STR = '@WF,WB,Wi,We [MJ] vs t@'
            gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,29), GVYL, NGVM+1, NGVV+1, 4, gDIVL)
         CALL TXGRFTX(0, GVX, GVYL, NGVM, NGVV, 4, STR, MODE, IND)

         STR = '@TAUE1,TAUE2,TAUEP [s] vs t@'
            gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,33), GVYL, NGVM+1, NGVV+1, 3, gDIVL)
         CALL TXGRFTX(1, GVX, GVYL, NGVM, NGVV, 3, STR, MODE, IND)

         STR = '@BETAa,BETA0 vs t@'
            gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,36), GVYL, NGVM+1, NGVV+1, 2, gDIVL)
         CALL TXGRFTX(2, GVX, GVYL, NGVM, NGVV, 2, STR, MODE, IND)

         STR = '@BETAPa,BETAP0 vs t@'
            gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,38), GVYL, NGVM+1, NGVV+1, 2, gDIVL)
         CALL TXGRFTX(3, GVX, GVYL, NGVM, NGVV, 2, STR, MODE, IND)

      ELSEIF (NGYV .EQ. 4) THEN
         STR = '@VLOOP [V]@'
         CALL TXGRFTX(0, GVX, GVY(0,40), NGVM, NGVV, 1, STR, MODE, IND)

         STR = '@ALI@'
         CALL TXGRFTX(1, GVX, GVY(0,41), NGVM, NGVV, 1, STR, MODE, IND)

         STR = '@Q(0)@'
         CALL TXGRFTX(2, GVX, GVY(0,42), NGVM, NGVV, 1, STR, MODE, IND)

         STR = '@RQ1 [m]@'
         CALL TXGRFTX(3, GVX, GVY(0,43), NGVM, NGVV, 1, STR, MODE, IND)

      ELSEIF (NGYV .EQ. 5) THEN
         GOTO 100

      ELSEIF (NGYV .EQ. 6) THEN
         STR = '@Te,Ti,<Te>,<Ti> [keV] vs t@'
         gDIVL=1.E0
         CALL DIVGY(GVY(0,1), GVYL, NGVM+1, NGVV+1, 4, gDIVL)
         CALL TXGRFVX(0, GVX, GVYL, NGVM, NGVV, 4, STR, IND)

         STR = '@PIN,POH,PNB,PRF,PNF [MW] vs t@'
            gDIVL = 1.E0
         CALL DIVGY(GVY(0,5), GVYL, NGVM+1, NGVV+1, 5, gDIVL)
         CALL TXGRFVX(4, GVX, GVYL, NGVM, NGVV, 5, STR, IND)

         STR = '@IP,IOH,INB,IOH+INB [MA] vs t@'
            gDIVL = 1.E0
         CALL DIVGY(GVY(0,10), GVYL, NGVM+1, NGVV+1, 4, gDIVL)
         CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 4, STR, IND)

         STR = '@POUT,PCX,PIE,PRL,PCON [MW] vs t@'
         IF (MODEG .EQ. 2) THEN
            gDIVL = gDIV(12-3)
         ELSE
            gDIVL = 1.E0
         ENDIF
         CALL DIVGY(GVY(0,15), GVYL, NGVM+1, NGVV+1, 5, gDIVL)
         CALL TXGRFVX(6, GVX, GVYL, NGVM, NGVV, 5, STR, IND)

         STR = '@QF vs t@'
         CALL TXGRFVX(1, GVX, GVY(0,20), NGVM, NGVV, 1, STR, IND)

         STR = '@Te0,Ti0 [keV] vs t@'
            gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,21), GVYL, NGVM+1, NGVV+1, 2, gDIVL)
         CALL TXGRFVX(5, GVX, GVYL, NGVM, NGVV, 2, STR, IND)

         STR = '@<Te>,<Ti> [keV] vs t@'
             gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,23), GVYL, NGVM+1, NGVV+1, 2, gDIVL)
         CALL TXGRFVX(3, GVX, GVYL, NGVM, NGVV, 2, STR, IND)

         STR = '@Ne0,Ni0,<Ne>,<Ni> [10~20/m^3] vs t@'
            gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,25), GVYL, NGVM+1, NGVV+1, 4, gDIVL)
         CALL TXGRFVX(7, GVX, GVYL, NGVM, NGVV, 4, STR, IND)

      ELSEIF (NGYV .EQ. 7) THEN
         STR = '@WF,WB,Wi,We [MJ] vs t@'
            gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,29), GVYL, NGVM+1, NGVV+1, 4, gDIVL)
         CALL TXGRFVX(0, GVX, GVYL, NGVM, NGVV, 4, STR, IND)

         STR = '@TAUE1,TAUE2,TAUEP [e] vs t@'
            gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,33), GVYL, NGVM+1, NGVV+1, 3, gDIVL)
         CALL TXGRFVX(4, GVX, GVYL, NGVM, NGVV, 3, STR, IND)

         STR = '@BETAa,BETA0 vs t@'
            gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,36), GVYL, NGVM+1, NGVV+1, 2, gDIVL)
         CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 2, STR, IND)

         STR = '@BETAPa,BETAP0 vs t@'
            gDIVL = 1.E0
!         CALL APPDIV(STR, gDIVL)
         CALL DIVGY(GVY(0,38), GVYL, NGVM+1, NGVV+1, 2, gDIVL)
         CALL TXGRFVX(6, GVX, GVYL, NGVM, NGVV, 2, STR, IND)

         STR = '@VLOOP [V]@'
         CALL TXGRFVX(1, GVX, GVY(0,40), NGVM, NGVV, 1, STR, IND)

         STR = '@ALI@'
         CALL TXGRFVX(5, GVX, GVY(0,41), NGVM, NGVV, 1, STR, IND)

         STR = '@Q(0)@'
         CALL TXGRFVX(3, GVX, GVY(0,42), NGVM, NGVV, 1, STR, IND)

         STR = '@RQ1 [m]@'
         CALL TXGRFVX(7, GVX, GVY(0,43), NGVM, NGVV, 1, STR, IND)

      ELSE
         WRITE(6,*) 'Unknown NGYV: NGYV = ',NGYV
      ENDIF

  100 CALL PAGEE

      RETURN
      END
!
!     ***************************************************************
!
!        SLAVE ROUTINE TO Write graph of GTX, GTY
!
!     ***************************************************************
!
      SUBROUTINE TXGRFTX(K, GTXL, GTYL, NGTLM, NGTL, NG, STR, MODE,IND)

      INCLUDE 'txcomm.inc'
      INTEGER, INTENT(IN) :: K ,NGTLM, NGTL, NG
      REAL, DIMENSION(0:NGTLM), INTENT(IN) :: GTXL
      REAL, DIMENSION(0:NGTLM,NG), INTENT(IN) :: GTYL
      CHARACTER(*), INTENT(IN) ::  STR
      INTEGER :: MODE, IND
      REAL, DIMENSION(4) :: GPXY

      GPXY(1) =  3.0 + 12.5 * MOD(K,2)
      GPXY(2) = 12.5 + 12.5 * MOD(K,2)
      GPXY(3) = 10.5 -  8.5 * REAL(K/2)
      GPXY(4) = 17.0 -  8.5 * REAL(K/2)
      CALL TXGRAF(GPXY, GTXL, GTYL, NGTLM+1, NGTL+1, NG,
     &            GTXL(0), GTXL(NGTL), STR, MODE, IND)

      RETURN
      END
!
!     ***************************************************************
!
!        SLAVE ROUTINE TO Write graph of GVX, GVY
!
!     ***************************************************************
!
      SUBROUTINE TXGRFVX(K, GTXL, GTYL, NGTLM, NGTL, NG, STR, IND)

      INCLUDE 'txcomm.inc'
      INTEGER, INTENT(IN) :: K ,NGTLM, NGTL, NG
      REAL, DIMENSION(0:NGTLM), INTENT(IN) :: GTXL
      REAL, DIMENSION(0:NGTLM,NG), INTENT(IN) :: GTYL
      CHARACTER(*), INTENT(IN) ::  STR
      INTEGER :: IND
      REAL, DIMENSION(4) :: GPXY

      GPXY(1) =  3.0 + 12.0 * MOD(K,2)
      GPXY(2) = 12.0 + 12.0 * MOD(K,2)
      GPXY(3) = 14.0 -  4.3 * REAL(K/2)
      GPXY(4) = 17.0 -  4.3 * REAL(K/2)
      CALL TXGRAF(GPXY, GTXL, GTYL, NGTLM+1, NGTL+1, NG,
     &            GTXL(0), GTXL(NGTL), STR, 1, IND)

      RETURN
      END
!
!     ***************************************************************
!
!        Graph of each TERM
!
!     ***************************************************************
!
      SUBROUTINE TXGRFQ(NQ,ID)

      INCLUDE 'txcomm.inc'

      INTEGER, INTENT(IN) :: NQ, ID
      INTEGER :: NR, NC, NC1, NSTR, NXMAX, IND
      INTEGER, DIMENSION(NQM) :: IHIGQ
      REAL :: GXMAX
      REAL, DIMENSION(0:NRM,NCM) :: GQY
      REAL, DIMENSION(4) :: GPXY
      REAL, DIMENSION(4,4) :: GPXYA
      CHARACTER(80), DIMENSION(NQM) :: STRGQ
      CHARACTER(80) :: STR

      DATA STRGQ /
     &     '@E$-r$=@','@E$-$#q$#$=@','@E$-$#f$#$=@',
     &                '@B$-$#q$#$=@','@B$-$#f$#$=@', 
     &     '@n$-e$=@','@n$-e$= u$-er$=@',
     &     '@n$-e$=u$-e$#q$#$=@','@n$-e$#f$#$=@','@T$-e$=@',
     &     '@n$-i$=@','@n$-i$= u$-ir$=@',
     &     '@n$-i$=u$-i$#q$#$=@','@n$-i$#f$#$=@','@T$-i$=@',
     &     '@n$-b$=@',
     &     '@n$-b$=u$-b$#q$#$=@','@n$-b$#f$#$=@',
     &     '@Slow n$-0$=@', '@Fast n$-0$=@'/
!
!        1 : Equation is integer mesh.
!        2 : Equation is half integer mesh.
      DATA IHIGQ/2,1,2,1,2, 2,1,1,2,2, 2,1,1,2,2, 2,1,1, 2,2/
      DATA GPXYA/ 2.5, 9.5,10.5,17.0,
     &            2.5, 9.5, 1.5, 8.0,
     &           15.0,22.0,10.5,17.0,
     &           15.0,22.0, 1.5, 8.0/

      GPXY(1) = GPXYA(1,ID)
      GPXY(2) = GPXYA(2,ID)
      GPXY(3) = GPXYA(3,ID)
      GPXY(4) = GPXYA(4,ID)

      DO NC = 1, NLCMAX(NQ)
         NR = 0
            NC1 = NLC(NC,NQ,NR)
            GQY(NR,NC) = SNGL(+ BLC(NC,NQ,NR) * X(NC1,NR  )
     &                        + ALC(NC,NQ,NR) * X(NC1,NR+1)
     &                        + PLC(NC,NQ,NR))
         DO NR = 1, NRMAX - 1
            NC1 = NLC(NC,NQ,NR)
            GQY(NR,NC) = SNGL(+ CLC(NC,NQ,NR) * X(NC1,NR-1)
     &                        + BLC(NC,NQ,NR) * X(NC1,NR  )
     &                        + ALC(NC,NQ,NR) * X(NC1,NR+1)
     &                        + PLC(NC,NQ,NR))
         ENDDO
         NR = NRMAX
            NC1 = NLC(NC,NQ,NR)
            GQY(NR,NC) = SNGL(+ CLC(NC,NQ,NR) * X(NC1,NR-1)
     &                        + BLC(NC,NQ,NR) * X(NC1,NR  )
     &                        + PLC(NC,NQ,NR))
      ENDDO

      NSTR = 0
      CALL APITOS(STR, NSTR, NQ)
      CALL APTTOS(STR, NSTR, '@: @')
      CALL APTTOS(STR, NSTR, STRGQ(NQ))

      IF (IHIGQ(NQ) .EQ. 1) THEN
         NXMAX = NRMAX+1
      ELSE
         NXMAX = NRMAX
      ENDIF
      IF (MODEG .EQ. 2) THEN
         IND = 9
      ELSE
         IND = 0
      ENDIF
      GXMAX=REAL(RB/RA)
      CALL TXGRAF(GPXY, GX(0,IHIGQ(NQ)), GQY,
     &            NRM+1, NXMAX, NLCMAX(NQ),
     &            0.0, GXMAX, '@'//STR(1:NSTR)//'@', 3, IND)

      RETURN
      END
!
!     ***************************************************************
!
!        Write parameter to graphic screen
!
!     ***************************************************************
!
      SUBROUTINE TXWPGR

      INCLUDE 'txcomm.inc'

      INTEGER :: IFNT, NP
      REAL :: GXM, GYM, GYS
      COMMON /TXPRS1/ GXM,GYM,GYS,NP

      CALL INQFNT(IFNT)
      CALL SETFNT(32)
      CALL SETCHS(0.3, 0.0)
      CALL SETLIN(0, 1, 7)

      GXM = 13.0 + 1.0
      GYM =  8.5 + 0.2
      GYS =  0.5
      NP  = 0

      CALL TXWPSS('+'//SLID//'+')
      CALL TXWPSD('@PNBCD @', PNBCD)
      CALL TXWPSI('@NRMAX @', NRMAX)

      NP = NP + 1
      CALL TXWPSD('@BB    @', BB)
      CALL TXWPSD('@rIp   @', rIp)
      CALL TXWPSD('@FSDFIX@', FSDFIX)
      CALL TXWPSD('@FSCDBM@', FSCDBM)
      CALL TXWPSD('@FSBOHM@', FSBOHM)
      CALL TXWPSD('@FSPSCL@', FSPSCL)
      CALL TXWPSD('@PROFD @', PROFD)
      CALL TXWPSD('@FSCX  @', FSCX)
      CALL TXWPSD('@FSLC  @', FSLC)
      CALL TXWPSD('@FSNC  @', FSNC)
      CALL TXWPSD('@FSLP  @', FSLP)
      CALL TXWPSD('@FSION @', FSION)
      CALL TXWPSD('@FSD0  @', FSD0)

      GXM = GXM + 0.35 * 17
      NP = 0
      CALL TXWPSD('@PNBH  @', PNBH)
      CALL TXWPSD('@PRFH  @', PRFH)
      CALL TXWPSD('@De0   @', De0)
      CALL TXWPSD('@Di0   @', Di0)
      CALL TXWPSD('@rMue0 @', rMue0)
      CALL TXWPSD('@rMui0 @', rMui0)
      CALL TXWPSD('@Chie0 @', Chie0)
      CALL TXWPSD('@Chii0 @', Chii0)
      CALL TXWPSD('@WPM0  @', WPM0)
      CALL TXWPSD('@PTe0  @', PTe0)
      CALL TXWPSD('@PTea  @', PTea)
      CALL TXWPSD('@PTi0  @', PTi0)
      CALL TXWPSD('@PTia  @', PTia)
      CALL TXWPSD('@V0    @', V0)
      CALL TXWPSD('@rGamm0@', rGamm0)
      CALL TXWPSD('@rGASPF@', rGASPF)
      CALL TXWPSD('@Zeff  @', Zeff)

      CALL SETFNT(IFNT)
      RETURN
      END
!
!     ***************************************************************
!
!        WPgr's Sub : write Double
!
!     ***************************************************************
!
      SUBROUTINE TXWPSD(STR, VAL)

      IMPLICIT NONE
      
      REAL(8), INTENT(IN) :: VAL
      CHARACTER(*), INTENT(IN) :: STR
      INTEGER :: NP
      REAL :: GXM, GYM, GYS
      COMMON /TXPRS1/ GXM,GYM,GYS,NP

      CALL MOVE(GXM, GYM - GYS * NP)
      CALL TEXTX(STR)
      CALL NUMBD(VAL,'(1PD9.2)', 9)
      NP = NP + 1

      RETURN
      END
!
!     ***************************************************************
!
!        WPgr's Sub : write Integer
!
!     ***************************************************************
!
      SUBROUTINE TXWPSI(STR, IVAL)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: IVAL
      CHARACTER(*), INTENT(IN) :: STR
      INTEGER :: NP
      REAL :: GXM, GYM, GYS
      COMMON /TXPRS1/ GXM,GYM,GYS,NP

      CALL MOVE(GXM, GYM - GYS * NP)
      CALL TEXTX(STR)
      CALL TEXT(' = ', 3)
      CALL NUMBI(IVAL,'(I6)',6)
      NP = NP + 1

      RETURN
      END
!
!     ***************************************************************
!
!        WPgr's Sub : write Strings
!
!     ***************************************************************
!
      SUBROUTINE TXWPSS(STR)

      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: STR
      INTEGER :: NP
      REAL :: GXM, GYM, GYS
      COMMON /TXPRS1/ GXM,GYM,GYS,NP

      CALL MOVE(GXM, GYM - GYS * NP)
      CALL TEXTX(STR)
      NP = NP + 1

      RETURN
      END
!
!     ***************************************************************
!
!        Write graph
!
!     ***************************************************************
!
      SUBROUTINE TXGRAF(GPXY, GX, GY, NXM, NXMAX, NGMAX,
     &                  GXMIN, GXMAX, STR, MODE, IND)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NXM, NXMAX, NGMAX ,MODE, IND
      REAL, INTENT(IN) :: GXMIN, GXMAX
      REAL, DIMENSION(4), INTENT(IN) :: GPXY
      REAL, DIMENSION(NXMAX), INTENT(IN) :: GX
      REAL, DIMENSION(NXM,NGMAX), INTENT(IN) :: GY
      CHARACTER(*), INTENT(IN) :: STR

      INTEGER :: IFNT, NGV, NGULEN, ICL, IPAT, IMRK, ISTEP, NG
      REAL :: GX1, GX2, GY1, GY2, gSLEN, GSXMIN, GSXMAX, GXSTEP,
     &        GYMIN, GYMAX, GSYMIN, GSYMAX, GYSTEP, GYORG, 
     &        GMRK, GCHH, GXL, GYL
      INTEGER, DIMENSION(0:4) :: NLTYPE
      DATA NLTYPE/0,2,3,4,6/

      IF (MODE.LT.0 .OR. MODE.GT.3) THEN
         WRITE(6,*) '### ERROR(TXGRAF) : MODE = ', MODE
         RETURN
      ENDIF

      GX1=GPXY(1)
      GX2=GPXY(2)
      GY1=GPXY(3)
      GY2=GPXY(4)

      CALL INQFNT(IFNT)

      IF (IND .EQ. 0) THEN
         gSLEN = 1.0
      ELSE
         gSLEN = 0.2
      ENDIF

      CALL SETFNT(32)
      CALL SETCHS(0.3, 0.0)
      CALL SETLIN(0, 0, 7)
      CALL GTEXTX(GX1,GY2+0.2,STR,0)

      CALL SETFNT(44)
      CALL SETCHS(0.3, 0.0)
      CALL SETLIN(0, 0, 7)
      CALL SETLNW(0.017)

      CALL GQSCAL(GXMIN, GXMAX, GSXMIN, GSXMAX, GXSTEP)
      GSXMIN = GXMIN
      GSXMAX = GXMAX
      CALL GMNMX2(GY,NXM,1,NXMAX,1,1,NGMAX,1,GYMIN,GYMAX)
      IF (GYMAX.GT.0.0) THEN
         IF (GYMIN.GT.0.0) THEN
            GYMIN=0.0
         ENDIF
      ELSE
         GYMAX=0.0
      ENDIF
      CALL GQSCAL(GYMIN, GYMAX, GSYMIN, GSYMAX, GYSTEP)
      IF (GSYMIN .GT. GYMIN) GSYMIN = GSYMIN - GYSTEP
      IF (GSYMAX .LT. GYMAX) GSYMAX = GSYMAX + GYSTEP
      GYORG = 0.0

      CALL GDEFIN(GX1, GX2, GY1, GY2, 
     &            GSXMIN, GSXMAX, GSYMIN, GSYMAX)
      CALL GFRAME
      CALL GSCALE(GSXMIN, GXSTEP, 0.0, 0.0, gSLEN, IND)
      IF (GXSTEP.LT.0.01) THEN
         NGV=NGULEN(GXSTEP*5)
         IF (NGV.LT.0) THEN
            NGV=NGV-3200
         ELSE
            NGV=NGV+3200
         ENDIF
         CALL GVALUE(GSXMIN, GXSTEP*5, 0.0, 0.0, NGV)
      ELSE
         NGV=NGULEN(GXSTEP*2)
         IF (NGV.LT.0) THEN
            NGV=NGV-3200
         ELSE
            NGV=NGV+3200
         ENDIF
         CALL GVALUE(GSXMIN, GXSTEP*2, 0.0, 0.0, NGV)
      ENDIF
      CALL GSCALE(0.0, 0.0, GYORG, GYSTEP, gSLEN, IND)
      IF (GSYMIN.LT.0.0 .AND. GSYMAX.GT.0.0)
     &     CALL GSCALE(0.0, 0.0, 0.0, GSYMAX-GSYMIN,  2.0, 0)
      CALL GVALUE(0.0,0.0,GYORG,GYSTEP*2,NGULEN(GYSTEP*2))

! MODE = 0: Change Line Color (Last Color Fixed)

      IF (MODE .EQ. 0) THEN
         DO NG = 1, NGMAX
            ICL = 7 - MOD(NGMAX - NG, 5)
            CALL SETLIN(0, 1, ICL)
            CALL GPLOTP(GX, GY(1,NG), 1, NXMAX, 1, 0, 0, 0)
         ENDDO

! MODE = 1: Change Line Color and Style

      ELSE IF (MODE .EQ. 1) THEN
!      IF (MODE .EQ. 0.or.MODE .EQ. 1) THEN
         DO NG = 1, NGMAX
            ICL  = 7 - MOD(NG-1, 5)
            IPAT = NLTYPE(MOD(NG-1, 5))
            CALL SETLIN(0, 1, ICL)
            CALL GPLOTP(GX, GY(1,NG), 1, NXMAX, 1, 0, 0, IPAT)
         ENDDO

! MODE = 2: Change Line Color, Style and Mark
! MODE = 3: Change Line Color, Style and Mark (With Legend)

      ELSE IF (MODE .EQ. 2 .OR. MODE .EQ. 3) THEN
         IMRK = 0
         GMRK = 0.3
         CALL SETMKS(IMRK, GMRK)
         DO NG = 1, NGMAX
            ICL = 7 - MOD(NG - 1, 5)
            IMRK  = MOD(NG - 1, 5) + 1
            ISTEP = NXMAX / 10
            IPAT  = (NG - 1) / 5
            CALL SETLIN(0, 1, ICL)
            CALL GPLOTP(GX, GY(1,NG), 1, NXMAX, 1, -IMRK, ISTEP, IPAT)
         ENDDO
! Legend
         IF (MODE .EQ. 3) THEN
            CALL SETCHS(0.25,0.0)
            GCHH = 0.25
            GXL = GX2 + GCHH
            GYL = GY2 - GCHH
            DO NG = 1, NGMAX
               ICL = 7 - MOD(NG - 1, 5)
               CALL SETLIN(0, 1, ICL)
               IPAT = (NG - 1) / 5
               CALL MOVEPT(GXL + GCHH * 3, GYL + GCHH / 2, IPAT)
               CALL DRAWPT(GXL + GCHH * 7, GYL + GCHH / 2)
               CALL MOVE(GXL, GYL)
               CALL NUMBI(NG, '(I2)', 2)
               IMRK = MOD(NG - 1, 5) + 1
               IF (IMRK .NE. 0) THEN
                  CALL SETMKS(IMRK, GMRK)
                  CALL MARK(GXL + GCHH * 5, GYL + GCHH / 2)
               ENDIF
               GYL = GYL - GCHH * 2
            ENDDO
         ENDIF
! End of Legend
         IMRK = 0
         GMRK = 0.2
         CALL SETMKS(IMRK, GMRK)
      ENDIF

      CALL SETFNT(IFNT)

      RETURN
      END
!
!     ***************************************************************
!
!        Divide GY, GTY by gDIV
!
!     ***************************************************************
!
      SUBROUTINE DIVGY(GIN, GOUT, NXM, NXMAX, NYMAX, gDIVL)

      INCLUDE 'txcomm.inc'

      INTEGER, INTENT(IN) :: NXM, NXMAX, NYMAX
      REAL, INTENT(IN) :: gDIVL
      REAL, DIMENSION(NXM,NYMAX), INTENT(IN) :: GIN
      REAL, DIMENSION(NXM,NYMAX), INTENT(OUT) :: GOUT
      INTEGER :: IX, IY

      DO IY = 1, NYMAX
         DO IX = 1, NXMAX
            GOUT(IX,IY) = GIN(IX,IY) / gDIVL
         ENDDO
      ENDDO

      RETURN 
      END
!
!     ***************************************************************
!
!        Append gDIV to STRX
!
!     ***************************************************************
!
      SUBROUTINE APPDIV(STRX, gDIVL)

      IMPLICIT NONE
      REAL, INTENT(IN) :: gDIVL
      CHARACTER(*), INTENT(IN) :: STRX
      INTEGER :: NSTRX, NSTRLEN

      IF (gDIVL .EQ. 1.E0) RETURN
      NSTRX = NSTRLEN(STRX) + 1
      CALL APSTOS(STRX, NSTRX, ' [', 2)
      CALL APRTOS(STRX, NSTRX, gDIVL, 'E0')
      CALL APSTOS(STRX, NSTRX, ']@', 2)

      RETURN
      END
