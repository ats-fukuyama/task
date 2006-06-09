!     $Id$

module graphic
  implicit none
  private
  real :: GXM,GYM,GYS
  integer :: NP
  public :: TXGOUT, TX_GRAPH_SAVE, TXSTGR, TXSTGT, TXSTGV

contains
  
  !***************************************************************
  !
  !   GRAPHIC Command loop
  !
  !***************************************************************
  
  SUBROUTINE TXGOUT
    use commons
    use libraries, only : TOUPPER
    use file_io, only : TXLOAD, TXGSAV, TXGLOD

    INTEGER :: MODE, NGPR, NGPT, NGPV, NQ, NQL, NGF, NGFMAX, I, IST
    character(len=5) :: STR, STR2
    character(len=1) :: KID1, KID2

    !     *** MENU ***

    OUTER : DO
       WRITE(6,*) '# SELECT : Rn: Tn: Un: Vn: A,B,C: ', &
            &           'S,L,M:file I:init X:exit'
       READ(5,'(A5)',IOSTAT=IST) STR
       IF (IST > 0) THEN
          WRITE(6,*) '### ERROR : Invalid Command : ', STR
          CYCLE
       ELSEIF (IST < 0) THEN
          RETURN
       END IF

       KID1 = STR(1:1)
       KID2 = STR(2:2)
       CALL TOUPPER(KID1)
       CALL TOUPPER(KID2)

       SELECT CASE(KID1)
       CASE('S')
          CALL TXGSAV

       CASE('L')
          CALL TXGLOD
          CALL TXPRFG

       CASE('I')
          T_TX = 0.D0
          TPRE = 0.D0
          NGT = -1
          NGVV = -1
          CALL TXSTGT(SNGL(T_TX))
          CALL TXSTGV(SNGL(T_TX))

       CASE('R')
          IF (KID2 == 'A') THEN
             MODE=MODEL
             CALL TXGRFR(-1,MODE)
          ELSEIF (KID2 == 'B') THEN
             MODE=MODEL
             CALL TXGRFR(-3,MODE)
          ELSEIF (KID2 == 'C') THEN
             MODE=MODEL
             CALL TXGRFR(-5,MODE)
          ELSE
             READ(STR(2:5),*,IOSTAT=IST) NGPR
             IF (IST /= 0) THEN
                WRITE(6,*) '### ERROR : Invalid Command : ', STR
                CYCLE
             END IF
             IF (NGPR >= 0 .AND. NGPR <= NGPRM) THEN
                MODE=MODEL
                CALL TXGRFR(NGPR,MODE)
             END IF
          END IF

       CASE('T')
          SELECT CASE(KID2)
          CASE('A')
             DO NGPT = 1, NGPTM
                MODE=MODEL
                CALL TXGRFT(NGPT,MODE)
             END DO
          CASE('B')
             DO NGPT = 1, 6
                MODE=MODEL
                CALL TXGRFT(NGPT,MODE)
             END DO
             DO NGPT = 9, 10
                MODE=MODEL
                CALL TXGRFT(NGPT,MODE)
             END DO

          CASE('C')
             DO NGPT = 1, 5
                MODE=MODEL
                CALL TXGRFT(NGPT,MODE)
             END DO

          CASE DEFAULT
             READ(STR(2:5),*,IOSTAT=IST) NGPT
             IF (IST < 0) THEN
                WRITE(6,*) '### ERROR : Invalid Command : ', STR
                CYCLE
             END IF
             IF (NGPT >= 1 .AND. NGPT <= NGPTM) THEN
                MODE=MODEL
                CALL TXGRFT(NGPT,MODE)
             END IF
          END SELECT

       CASE('V')
          SELECT CASE(KID2)
          CASE('A')
             DO NGPV = 1, NGPVM
                MODE=MODEL
                CALL TXGRFV(NGPV,MODE)
             END DO

          CASE('B')
             DO NGPV = 1, 6
                MODE=MODEL
                CALL TXGRFV(NGPV,MODE)
             END DO
             DO NGPV = 9, 10
                MODE=MODEL
                CALL TXGRFV(NGPV,MODE)
             END DO

          CASE('C')
             DO NGPV = 1, 5
                MODE=MODEL
                CALL TXGRFV(NGPV,MODE)
             END DO

          CASE DEFAULT
             READ(STR(2:5),*,IOSTAT=IST) NGPV
             IF (IST < 0) THEN
                WRITE(6,*) '### ERROR : Invalid Command : ', STR
                CYCLE
             END IF
             IF (NGPV >= 1 .AND. NGPV <= NGPVM) THEN
                MODE=MODEL
                CALL TXGRFV(NGPV,MODE)
             END IF
          END SELECT

       CASE('U')
          IF (KID2 == 'A') THEN
             DO NQ = 1, NQMAX, 4
                CALL PAGES
                DO NQL=NQ,MIN(NQ+3,NQMAX)
                   CALL TXGRFQ(NQL,MOD(NQL-1,4)+1)
                END DO
                CALL PAGEE
             END DO

          ELSE
             READ(STR(2:5),*,IOSTAT=IST) NQ
             IF (IST < 0) THEN
                WRITE(6,*) '### ERROR : Invalid Command : ', STR
                CYCLE
             END IF
             IF (NQ >= 1 .AND. NQ <= NQMAX) THEN
                CALL PAGES
!!$              DO NQL=NQ,MIN(NQ+3,NQMAX)
!!$                 CALL TXGRFQ(NQL,NQL-NQ+1)
!!$              END DO
                CALL TXGRFQ(NQ,5)
                CALL PAGEE
             END IF
          END IF

       CASE('A')
          DO NGPR = 1, NGPRM
             MODE=MODEL
             CALL TXGRFR(NGPR,MODE)
          END DO

       CASE('B')
          DO NGPR = 1, 6
             MODE=MODEL
             CALL TXGRFR(NGPR,MODE)
          END DO
          DO NGPR = 9, 10
             MODE=MODEL
             CALL TXGRFR(NGPR,MODE)
          END DO

       CASE('C')
          DO NGPR = 1, 5
             MODE=MODEL
             CALL TXGRFR(NGPR,MODE)
          END DO

       CASE('M')
          DO
             WRITE(6,*) '## Number of files :'
             READ(5,*,IOSTAT=IST) NGFMAX
             IF (IST > 0) THEN
                CYCLE
             ELSEIF (IST < 0) THEN
                CYCLE OUTER
             ELSE
                EXIT
             END IF
          END DO
          NGR=-1
          DO NGF=1,NGFMAX
             CALL TXLOAD
             CALL TX_GRAPH_SAVE
          END DO
          DO
             WRITE(6,*) '## INPUT GRAPH NUMBER'
             READ(5,'(A5)',IOSTAT=IST) STR2
             IF (IST > 0) THEN
                CYCLE
             ELSEIF (IST < 0) THEN
                CYCLE OUTER
             ELSE
                EXIT
             END IF
          END DO
          CALL TOUPPER(STR2)
          IF (STR2 == '     ') THEN
             !     For space or return only 
             !     Correspond to GMA
          ELSE IF (STR2(1:1) == 'A') THEN
             DO I = 1, NGPRM
                MODE=MODEL
                CALL TXGRFR(I,MODE)
             END DO
             !     Correspond to GMB
          ELSE IF (STR2(1:1) == 'B') THEN
             DO I = 1, 6
                MODE=MODEL
                CALL TXGRFR(I,MODE)
             END DO
             DO I = 9, 10
                MODE=MODEL
                CALL TXGRFR(I,MODE)
             END DO
             !     Correspond to GMC
          ELSE IF (STR2(1:1) == 'C') THEN
             DO I = 1, 5
                MODE=MODEL
                CALL TXGRFR(I,MODE)
             END DO
          ELSE
             READ(STR2,'(I5)',IOSTAT=IST) NGPR
             IF (IST < 0) CYCLE
             IF      (NGPR == 0) THEN
                CYCLE
             ELSE IF (NGPR >= 0 .AND. NGPR <= NGPRM) THEN
                MODE=MODEL
                CALL TXGRFR(NGPR,MODE)
             END IF
          END IF

       CASE('X')
          RETURN

       CASE DEFAULT
          WRITE(6,*) 'XX UNKNOWN GRAPHIC COMMAND'
       END SELECT
    END DO OUTER

    RETURN
  END SUBROUTINE TXGOUT

  !***************************************************************
  !
  !   Save graphic data
  !
  !***************************************************************

  SUBROUTINE TX_GRAPH_SAVE

    use commons, only : T_TX

    !  Define radial coordinate for graph

    CALL TXPRFG

    !  Store center or edge values of variables for showing time-evolution graph

    CALL TXSTGT(SNGL(T_TX))

    !  Store global quantities for showing time-evolution graph

    CALL TXSTGV(SNGL(T_TX))

    !  Store profile data for showing graph

    CALL TXSTGR

  END SUBROUTINE TX_GRAPH_SAVE

  !***************************************************************
  !
  !   Initialize graphic axes
  !
  !***************************************************************

  SUBROUTINE TXPRFG

    use commons, only : NRMAX, GX, R, RA
    INTEGER :: NR

    !  GX(NR) : Integer

    GX(0:NRMAX) = SNGL(R(0:NRMAX) / RA)

    RETURN
  END SUBROUTINE TXPRFG

  !***************************************************************
  !
  !   Store GY
  !
  !***************************************************************

  SUBROUTINE TXSTGR

    use commons
    use physical_constants, only : AEE

    INTEGER :: NR

    IF (NGR < NGRM) NGR = NGR + 1

    GT(NGR) = SNGL(T_TX)

    GY(0:NRMAX,NGR,1)  = SNGL(PNeV(0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,2)  = SNGL((PZ*PNiV(0:NRMAX)+PZ*PNbV(0:NRMAX)-PNeV(0:NRMAX))* 1.D20)
    GY(0:NRMAX,NGR,3)  = SNGL(UerV(0:NRMAX))
    GY(0:NRMAX,NGR,4)  = SNGL(UethV(0:NRMAX))
    GY(0:NRMAX,NGR,5)  = SNGL(UephV(0:NRMAX))
    GY(0:NRMAX,NGR,6)  = SNGL(UirV(0:NRMAX))
    GY(0:NRMAX,NGR,7)  = SNGL(UithV(0:NRMAX))
    GY(0:NRMAX,NGR,8)  = SNGL(UiphV(0:NRMAX))
    GY(0:NRMAX,NGR,9)  = SNGL(ErV(0:NRMAX))
    GY(0:NRMAX,NGR,10) = SNGL(BthV(0:NRMAX))
    GY(0:NRMAX,NGR,11) = SNGL(EphV(0:NRMAX))
    GY(0:NRMAX,NGR,12) = SNGL(PNbV(0:NRMAX) * 1.D20)
    GY(0:NRMAX,NGR,13) = SNGL(UbphV(0:NRMAX))
    GY(0:NRMAX,NGR,14) = SNGL(PTeV(0:NRMAX))
    GY(0:NRMAX,NGR,15) = SNGL(PTiV(0:NRMAX))
    GY(0:NRMAX,NGR,16) = SNGL((PN01V(0:NRMAX)+PN02V(0:NRMAX))*1.D20)
    GY(0:NRMAX,NGR,17) = SNGL(EthV(0:NRMAX))
    GY(0:NRMAX,NGR,18) = SNGL(BphV(0:NRMAX))
    GY(0:NRMAX,NGR,19) = SNGL(UbthV(0:NRMAX))
    GY(0:NRMAX,NGR,20) = SNGL(Q(0:NRMAX))

    GY(0:NRMAX,NGR,21) = SNGL((-   AEE*PNeV(0:NRMAX)*UephV(0:NRMAX) &
         &                     +PZ*AEE*PNiV(0:NRMAX)*UiphV(0:NRMAX) &
         &                     +PZ*AEE*PNbV(0:NRMAX)*UbphV(0:NRMAX))*1.D20)
    GY(0:NRMAX,NGR,22) = SNGL( -   AEE*PNeV(0:NRMAX)*UephV(0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,23) = SNGL(  PZ*AEE*PNiV(0:NRMAX)*UiphV(0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,24) = SNGL(  PZ*AEE*PNbV(0:NRMAX)*UbphV(0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,25) = SNGL(( BphV(0:NRMAX)*UethV(0:NRMAX) &
         &                     -BthV(0:NRMAX)*UephV(0:NRMAX)) &
         &                    /SQRT(BphV(0:NRMAX)**2 + BthV(0:NRMAX)*2))
    GY(0:NRMAX,NGR,26) = SNGL(( BthV(0:NRMAX)*UethV(0:NRMAX) &
         &                     +BphV(0:NRMAX)*UephV(0:NRMAX)) &
         &                    /SQRT(BphV(0:NRMAX)**2 + BthV(0:NRMAX)*2))
    GY(0:NRMAX,NGR,27) = SNGL(( BphV(0:NRMAX)*UithV(0:NRMAX) &
         &                     -BthV(0:NRMAX)*UiphV(0:NRMAX)) &
         &                    /SQRT(BphV(0:NRMAX)**2 + BthV(0:NRMAX)*2))
    GY(0:NRMAX,NGR,28) = SNGL(( BthV(0:NRMAX)*UithV(0:NRMAX) &
         &                     +BphV(0:NRMAX)*UiphV(0:NRMAX)) &
         &                    /SQRT(BphV(0:NRMAX)**2 + BthV(0:NRMAX)*2))
    GY(0:NRMAX,NGR,29) = SNGL(Di(0:NRMAX)+De(0:NRMAX))
    DO NR = 0, NRMAX
       IF (rG1h2(NR) > 3.D0) THEN
          GY(NR,NGR,30) = 3.0
       ELSE
          GY(NR,NGR,30) = SNGL(rG1h2(NR))
       END IF
    END DO
    GY(0:NRMAX,NGR,31) = SNGL(FCDBM(0:NRMAX))
    GY(0:NRMAX,NGR,32) = SNGL(S(0:NRMAX))
    GY(0:NRMAX,NGR,33) = SNGL(Alpha(0:NRMAX))
    GY(0:NRMAX,NGR,34) = SNGL(rKappa(0:NRMAX))
    GY(0:NRMAX,NGR,35) = SNGL(PN01V(0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,36) = SNGL(PN02V(0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,37) = SNGL(PNiV (0:NRMAX)*1.D20)

    !  Raw data

    GY(0:NRMAX,NGR,38) = SNGL(X(LQm1,0:NRMAX))
    GY(0:NRMAX,NGR,39) = SNGL(X(LQm2,0:NRMAX))
    GY(0:NRMAX,NGR,40) = SNGL(X(LQm3,0:NRMAX))
    GY(0:NRMAX,NGR,41) = SNGL(X(LQm4,0:NRMAX))
    GY(0:NRMAX,NGR,42) = SNGL(X(LQm5,0:NRMAX))
    GY(0:NRMAX,NGR,43) = SNGL(X(LQe1,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,44) = SNGL(X(LQe2,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,45) = SNGL(X(LQe3,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,46) = SNGL(X(LQe4,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,47) = SNGL(X(LQe5,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,48) = SNGL(X(LQi1,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,49) = SNGL(X(LQi2,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,50) = SNGL(X(LQi3,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,51) = SNGL(X(LQi4,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,52) = SNGL(X(LQi5,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,53) = SNGL(X(LQb1,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,54) = SNGL(X(LQb3,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,55) = SNGL(X(LQb4,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,56) = SNGL(X(LQn1,0:NRMAX)*1.D14)
    GY(0:NRMAX,NGR,57) = SNGL(X(LQn2,0:NRMAX)*1.D14)

    !  Coefficients

    GY(0:NRMAX,NGR,58) = SNGL(rMue(0:NRMAX))
    GY(0:NRMAX,NGR,59) = SNGL(rMui(0:NRMAX))
    GY(0:NRMAX,NGR,60) = SNGL(rNuei(0:NRMAX))
    GY(0:NRMAX,NGR,61) = SNGL(rNuL(0:NRMAX))
    GY(0:NRMAX,NGR,62) = SNGL(rNube(0:NRMAX))
    GY(0:NRMAX,NGR,63) = SNGL(rNubi(0:NRMAX))
    GY(0:NRMAX,NGR,64) = SNGL(FWthe(0:NRMAX))
    GY(0:NRMAX,NGR,65) = SNGL(FWthi(0:NRMAX))
    GY(0:NRMAX,NGR,66) = SNGL(WPM(0:NRMAX))
    GY(0:NRMAX,NGR,67) = SNGL(rNuTei(0:NRMAX))
    GY(0:NRMAX,NGR,68) = SNGL(rNu0e(0:NRMAX))
    GY(0:NRMAX,NGR,69) = SNGL(rNu0i(0:NRMAX))
    GY(0:NRMAX,NGR,70) = SNGL(Chie(0:NRMAX))
    GY(0:NRMAX,NGR,71) = SNGL(Chii(0:NRMAX))
    GY(0:NRMAX,NGR,72) = SNGL(D01(0:NRMAX))
    GY(0:NRMAX,NGR,73) = SNGL(D02(0:NRMAX))
    GY(0:NRMAX,NGR,74) = SNGL(rNueNC(0:NRMAX))
    GY(0:NRMAX,NGR,75) = SNGL(rNuiNC(0:NRMAX))
    GY(0:NRMAX,NGR,76) = SNGL(rNueHL(0:NRMAX))
    GY(0:NRMAX,NGR,77) = SNGL(rNuiHL(0:NRMAX))
    GY(0:NRMAX,NGR,78) = SNGL(rNuiCX(0:NRMAX))
    GY(0:NRMAX,NGR,79) = SNGL(rNuB(0:NRMAX))
    GY(0:NRMAX,NGR,80) = SNGL(rNuION(0:NRMAX))
    GY(0:NRMAX,NGR,81) = SNGL(SiLC(0:NRMAX))
    GY(0:NRMAX,NGR,82) = SNGL(SiLCth(0:NRMAX))
    GY(0:NRMAX,NGR,83) = SNGL(SiLCph(0:NRMAX))
    GY(0:NRMAX,NGR,84) = SNGL(SNB(0:NRMAX))
    GY(0:NRMAX,NGR,85) = SNGL(PRFi(0:NRMAX))
    GY(0:NRMAX,NGR,86) = SNGL(PRFe(0:NRMAX))

    GY(0:NRMAX,NGR,87) = SNGL(PIE(0:NRMAX))
    GY(0:NRMAX,NGR,88) = SNGL(PCX(0:NRMAX))
    GY(0:NRMAX,NGR,89) = SNGL(POH(0:NRMAX))

    RETURN
  END SUBROUTINE TXSTGR

  !***************************************************************
  !
  !   Store GTY
  !
  !***************************************************************

  SUBROUTINE TXSTGT(GTIME)

    use commons
    use physical_constants, only : AEE, PI
    use output_console, only : rLINEAVE
    use libraries, only : VALINT_SUB

    REAL, INTENT(IN) :: GTIME
    REAL(8) :: BthL, BphL, BBL, PNESUM1, PNESUM2

    IF (NGT < NGTM) NGT=NGT+1

    GTX(NGT) = GTIME

    GTY(NGT,1)  = SNGL(PNeV(0) * 1.D20)
    GTY(NGT,2)  = SNGL(rLINEAVE(0.D0))
    GTY(NGT,3)  = SNGL(rLINEAVE(0.24D0))
    GTY(NGT,4)  = SNGL(rLINEAVE(0.6D0))
    GTY(NGT,5)  = SNGL((PZ * PNiV(0) + PZ * PNbV(0) - PNeV(0)) * 1.D20)
    GTY(NGT,6)  = SNGL(UerV(NRMAX/2))
    GTY(NGT,7)  = SNGL(UethV(NRMAX/2))
    GTY(NGT,8)  = SNGL(UephV(0))
    GTY(NGT,9)  = SNGL(UirV(NRMAX/2))
    GTY(NGT,10) = SNGL(UithV(NRMAX/2))
    GTY(NGT,11) = SNGL(UiphV(0))
    GTY(NGT,12) = SNGL(ErV(NRMAX/2))
    GTY(NGT,13) = SNGL(BthV(NRMAX/2))
    GTY(NGT,14) = SNGL(EphV(NRMAX))
    GTY(NGT,15) = SNGL(PNbV(0) * 1.D20)
    GTY(NGT,16) = SNGL(UbphV(0))
    GTY(NGT,17) = SNGL(PTeV(0))
    GTY(NGT,18) = SNGL(PTiV(0))
    GTY(NGT,19) = SNGL((PN01V(0) + PN02V(0)) * 1.D20)

    GTY(NGT,20) = SNGL(Q(0))
    GTY(NGT,21) = SNGL((-      AEE * PNeV(0) * UephV(0) &
         &              + PZ * AEE * PNiV(0) * UiphV(0) &
         &              + PZ * AEE * PNbV(0) * UbphV(0)) * 1.D20)
    GTY(NGT,22) = SNGL(-      AEE * PNeV(0) * UephV(0) * 1.D20)
    GTY(NGT,23) = SNGL(  PZ * AEE * PNiV(0) * UiphV(0) * 1.D20)
    GTY(NGT,24) = SNGL(  PZ * AEE * PNbV(0) * UbphV(0) * 1.D20)

    BthL = BthV(NRMAX/2)
    BphL = BphV(NRMAX/2)
    BBL = SQRT(BphL**2 + BthL**2)
    GTY(NGT,25) = SNGL((+ BphL * UethV(NRMAX/2) &
         &              - BthL * UephV(NRMAX/2)) / BBL)
    GTY(NGT,26) = SNGL((+ BthL * UethV(NRMAX/2) &
         &              + BphL * UephV(NRMAX/2)) / BBL)
    GTY(NGT,27) = SNGL((+ BphL * UithV(NRMAX/2) &
         &              - BthL * UiphV(NRMAX/2)) / BBL)
    GTY(NGT,28) = SNGL((+ BthL * UithV(NRMAX/2) &
         &              + BphL * UiphV(NRMAX/2)) / BBL)
    GTY(NGT,29) = SNGL(Di(NRMAX/2)+De(NRMAX/2))
    GTY(NGT,30) = SNGL(rG1h2(NRMAX/2))
    GTY(NGT,31) = SNGL(FCDBM(NRMAX/2))
    GTY(NGT,32) = SNGL(S(NRMAX/2))
    GTY(NGT,33) = SNGL(Alpha(NRMAX/2))
    GTY(NGT,34) = SNGL(rKappa(NRMAX/2))

    !  GTY(NGT,19) = SNGL((PN01V(0) + PN02V(0)) * 1.D20)
    GTY(NGT,38) = SNGL(PN01V(0) * 1.D20)
    GTY(NGT,39) = SNGL(PN02V(0) * 1.D20)

    CALL VALINT_SUB(PNeV,NRA,PNESUM1)
    PNESUM1 = 2.D0*PI*RR*2.D0*PI*(PNESUM1 + PNeV(NRA)*R(NRA)*(RA-R(NRA)))
    CALL VALINT_SUB(PNeV(0:NRMAX)*rNuION(0:NRMAX),NRA,PNESUM2)
    PNESUM2 = 2.D0*PI*RR*2.D0*PI*(PNESUM2 + PNeV(NRA)*rNuION(NRA)*R(NRA)*(RA-R(NRA)))

    GTY(NGT,35) = SNGL(PNESUM1)
    GTY(NGT,36) = SNGL(PNESUM2)
    IF(NGT == 0.OR.ABS(PNESUM2) <= 0.D0) THEN
       GTY(NGT,37) = 0.0
    ELSE
       GTY(NGT,37) = SNGL(PNESUM1/PNESUM2)
    END IF

    RETURN
  END SUBROUTINE TXSTGT

  !***************************************************************
  !
  !   Store GVY
  !
  !***************************************************************

  SUBROUTINE TXSTGV(GTIME)

    use commons
    REAL, INTENT(IN) :: GTIME

    IF (NGVV < NGVM) NGVV=NGVV+1

    GVX(NGVV) = GTIME

    GVY(NGVV,1)  = SNGL(TS0(1))
    GVY(NGVV,2)  = SNGL(TS0(2))
    GVY(NGVV,3)  = SNGL(TSAV(1))
    GVY(NGVV,4)  = SNGL(TSAV(2))
    GVY(NGVV,5)  = SNGL(PINT)
    GVY(NGVV,6)  = SNGL(POHT)
    GVY(NGVV,7)  = SNGL(PNBT)
    GVY(NGVV,8)  = SNGL(PRFT)

    GVY(NGVV,10) = SNGL(AJT)
    GVY(NGVV,11) = SNGL(AJOHT)
    GVY(NGVV,12) = SNGL(AJNBT)
    GVY(NGVV,13) = SNGL(AJOHT+AJNBT)
    !  GVY(NGVV,14) = SNGL(AJBST)
    GVY(NGVV,15) = SNGL(POUT)
    GVY(NGVV,16) = SNGL(PCXT)
    GVY(NGVV,17) = SNGL(PIET)
    !  GVY(NGVV,18) = SNGL(PRLT)
    !  GVY(NGVV,19) = SNGL(PCONT)
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
  END SUBROUTINE TXSTGV

  !***************************************************************
  !
  !   Time evolution of radial profile
  !
  !***************************************************************

  SUBROUTINE TXGRFR(NGYRIN,MODE)

    use commons
    INTEGER, INTENT(IN) :: MODE
    INTEGER, INTENT(IN) :: NGYRIN
    INTEGER :: IND, NG, NR, NGYR, NE
    REAL(4), DIMENSION(0:NRM,0:NGRM) :: GYL
    character(len=40) :: STR
    character(len=1) :: KSTR
    character(len=3) :: KEND
    real, dimension(4) :: GPX, GPY
    real :: GPXL, FACT

    NGYR = NGYRIN

!    IF (NGYR == 0) NGYR = -1

    IF (NGR <= -1) THEN
       WRITE(6,*) 'G', NGYR, ' has no data'
       RETURN
    END IF

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF

    DO

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

    SELECT CASE(NGYR)
    CASE(0)
       ! Draw frame
       GPX(1) =  2.0 ; GPY(1) =  1.5
       GPX(2) =  2.0 ; GPY(2) = 17.0
       GPX(3) = 24.0 ; GPY(3) = 17.0
       GPX(4) = 24.0 ; GPY(4) =  1.5
       CALL LINES2D(GPX,GPY,4)
       ! Draw lines
       FACT = (GPX(4) - GPX(1)) / R(NRMAX)
       GPXL = GPX(1)
       DO NE = 1, NEMAX-1
          GPXL = GPXL + REAL(H(NE)) * FACT
          CALL LINE2D(GPXL,1.5,GPXL,17.0)
       END DO
       WRITE(KSTR,'(I1)') 0
       WRITE(KEND,'(I3)') NRMAX
       CALL GTEXT2D(GPX(1),GPY(1)-0.5,KSTR,1,2)
       CALL GTEXT2D(GPX(4),GPY(4)-0.5,KEND,3,2)

    CASE(1) 
       STR = '@n$-e$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,1), GYL, STR, NRM, NRMAX, NGR, gDIV(1))
       CALL TXGRFRX(0, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@Z*n$-i$=+Z*n$-b$=-n$-e$=@'
       CALL APPROPGY(MODEG, GY(0,0,2), GYL, STR, NRM, NRMAX, NGR, gDIV(2))
       CALL TXGRFRX(1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@E$-r$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,9), GYL, STR, NRM, NRMAX, NGR, gDIV(9))
       CALL TXGRFRX(2, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       CALL TXWPGR

    CASE(2)
       STR = '@u$-er$=(r)@'
       CALL TXGRFRX(0, GX, GY(0,0,3), NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-e$#q$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,4), GYL, STR, NRM, NRMAX, NGR, gDIV(4))
       CALL TXGRFRX(1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-e$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,5), GYL, STR, NRM, NRMAX, NGR, gDIV(5))
       CALL TXGRFRX(2, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       CALL TXWPGR

    CASE(3)
       STR = '@u$-ir$=(r)@'
       CALL TXGRFRX(0, GX, GY(0,0,6), NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-i$#q$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,7), GYL, STR, NRM, NRMAX, NGR, gDIV(7))
       CALL TXGRFRX(1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-i$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,8), GYL, STR, NRM, NRMAX, NGR, gDIV(8))
       CALL TXGRFRX(2, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       CALL TXWPGR

    CASE(4)
       STR = '@q(r)@'
       CALL TXGRFRX(0,GX,GY(0,0,20),NRMAX,NGR,STR,MODE,IND)

       STR = '@E$-$#f$#$=(r)@'
       CALL TXGRFRX(1,GX,GY(0,0,11),NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,21), GYL, STR, NRM, NRMAX, NGR, gDIV(21))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       CALL TXWPGR

    CASE(5)
       STR = '@n$-b$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,12), GYL, STR, NRM, NRMAX, NGR, gDIV(12))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-b$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,13), GYL, STR, NRM, NRMAX, NGR, gDIV(13))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@B$-$#q$#$=(r)@'
       CALL TXGRFRX(2,GX,GY(0,0,10),NRMAX,NGR,STR,MODE,IND)

       CALL TXWPGR

    CASE(6)
       STR = '@j$-e$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,22), GYL, STR, NRM, NRMAX, NGR, gDIV(22))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-i$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,23), GYL, STR, NRM, NRMAX, NGR, gDIV(23))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-b$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,24), GYL, STR, NRM, NRMAX, NGR, gDIV(24))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       CALL TXWPGR

    CASE(7)
       STR = '@u$-er$=(r)@'
       CALL TXGRFRX(0,GX,GY(0,0,3),NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-e$#$/136$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,25), GYL, STR, NRM, NRMAX, NGR, gDIV(25))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-e//$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,26), GYL, STR, NRM, NRMAX, NGR, gDIV(26))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       CALL TXWPGR

    CASE(8)
       STR = '@u$-ir$=(r)@'
       CALL TXGRFRX(0,GX,GY(0,0,6),NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-i$#$/136$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,27), GYL, STR, NRM, NRMAX, NGR, gDIV(27))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-i//$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,28), GYL, STR, NRM, NRMAX, NGR, gDIV(28))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       CALL TXWPGR

    CASE(9)
       STR = '@D$-i$=(r)+D$-e$=(r)@'
       CALL TXGRFRX(0,GX,GY(0,0,29),NRMAX,NGR,STR,MODE,IND)
       STR = '@G$-1$=h$+2$=(r)@'
       CALL TXGRFRX(1,GX,GY(0,0,30),NRMAX,NGR,STR,MODE,IND)
       STR = '@F$-CDBM$=(r)@'
       CALL TXGRFRX(2,GX,GY(0,0,31),NRMAX,NGR,STR,MODE,IND)
       CALL TXWPGR

    CASE(10)
       STR = '@T$-e$=(r)@'
       CALL TXGRFRX(0,GX,GY(0,0,14),NRMAX,NGR,STR,MODE,IND)

       STR = '@T$-i$=(r)@'
       CALL TXGRFRX(1,GX,GY(0,0,15),NRMAX,NGR,STR,MODE,IND)

       STR = '@Chi$-e$=(r)@'
       CALL TXGRFRX(2,GX,GY(0,0,70),NRMAX,NGR,STR,MODE,IND)

       CALL TXWPGR

    CASE(11)
       STR = '@S(r)@'
       CALL TXGRFRX(0,GX,GY(0,0,32),NRMAX,NGR,STR,MODE,IND)
       STR = '@$#a$#(r)@'
       CALL TXGRFRX(1,GX,GY(0,0,33),NRMAX,NGR,STR,MODE,IND)
       STR = '@$#k$#(r)@'
       CALL TXGRFRX(2,GX,GY(0,0,34),NRMAX,NGR,STR,MODE,IND)
       CALL TXWPGR

    CASE(12)
       STR = '@E$-$#q$#$=(r)@'
       CALL TXGRFRX(0,GX,GY(0,0,17),NRMAX,NGR,STR,MODE,IND)
       STR = '@B$-$#f$#$=(r)@'
       CALL TXGRFRX(1,GX,GY(0,0,18),NRMAX,NGR,STR,MODE,IND)
       STR = '@u$-b$#q$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,19), GYL, STR, NRM, NRMAX, NGR, gDIV(19))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)
       CALL TXWPGR

    CASE(13)
       STR = '@SLOW N$-0$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,35), GYL, STR, NRM, NRMAX, NGR, gDIV(35))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@FAST N$-0$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,36), GYL, STR, NRM, NRMAX, NGR, gDIV(36))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@TOTAL N$-0$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,16), GYL, STR, NRM, NRMAX, NGR, gDIV(16))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       CALL TXWPGR

    CASE(14)
       STR = '@PIE(r)@'
       CALL TXGRFRX(0,GX,GY(0,0,87),NRMAX,NGR,STR,MODE,IND)

       STR = '@PCX(r)@'
       CALL TXGRFRX(1,GX,GY(0,0,88),NRMAX,NGR,STR,MODE,IND)

       STR = '@POH(r)@'
       CALL APPROPGY(MODEG, GY(0,0,89), GYL, STR, NRM, NRMAX, NGR, gDIV(79))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       CALL TXWPGR

    CASE(-1)
       STR = '@E$-r$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,9), GYL, STR, NRM, NRMAX, NGR, gDIV(9))
       CALL TXGRFRXS(0, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@E$-$#q$#$=(r)@'
       CALL TXGRFRXS(1,GX,GY(0,0,17),NRMAX,NGR,STR,MODE,IND)

       STR = '@E$-$#f$#$=(r)@'
       CALL TXGRFRXS(2,GX,GY(0,0,11),NRMAX,NGR,STR,MODE,IND)

       STR = '@T$-e$=(r)@'
       CALL TXGRFRXS(3,GX,GY(0,0,14),NRMAX,NGR,STR,MODE,IND)

       STR = '@n$-e$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,1), GYL, STR, NRM, NRMAX, NGR, gDIV(1))
       CALL TXGRFRXS(4, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@B$-$#q$#$=(r)@'
       CALL TXGRFRXS(5,GX,GY(0,0,10),NRMAX  ,NGR,STR,MODE,IND)

       STR = '@B$-$#f$#$=(r)@'
       CALL TXGRFRXS(6,GX,GY(0,0,18),NRMAX,NGR,STR,MODE,IND)

       STR = '@T$-i$=(r)@'
       CALL TXGRFRXS(7,GX,GY(0,0,15),NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-er$=(r)@'
       CALL TXGRFRXS(8, GX, GY(0,0,3), NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-e$#q$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,4), GYL, STR, NRM, NRMAX, NGR, gDIV(4))
       CALL TXGRFRXS(9, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-e$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,5), GYL, STR, NRM, NRMAX, NGR, gDIV(5))
       CALL TXGRFRXS(10, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@SLOW N$-0$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,35), GYL, STR, NRM, NRMAX, NGR, gDIV(35))
       CALL TXGRFRXS(11,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-ir$=(r)@'
       CALL TXGRFRXS(12, GX, GY(0,0,6), NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-i$#q$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,7), GYL, STR, NRM, NRMAX, NGR, gDIV(7))
       CALL TXGRFRXS(13, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-i$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,8), GYL, STR, NRM, NRMAX, NGR, gDIV(8))
       CALL TXGRFRXS(14, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@FAST N$-0$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,36), GYL, STR, NRM, NRMAX, NGR, gDIV(36))
       CALL TXGRFRXS(15,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    CASE(-2)
       STR = '@Z*n$-i$=+Z*n$-b$=-n$-e$=@'
       CALL APPROPGY(MODEG, GY(0,0,2), GYL, STR, NRM, NRMAX, NGR, gDIV(2))
       CALL TXGRFRXS(0, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@n$-b$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,12), GYL, STR, NRM, NRMAX, NGR, gDIV(12))
       CALL TXGRFRXS(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-b$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,13), GYL, STR, NRM, NRMAX, NGR, gDIV(13))
       CALL TXGRFRXS(2,GX,GYL,NRMAX  ,NGR,STR,MODE,IND)     

       STR = '@u$-b$#q$#$=(r)@'
       CALL TXGRFRXS(3,GX,GY(0,0,19),NRMAX,NGR,STR,MODE,IND)

       STR = '@q(r)@'
       CALL TXGRFRXS(4,GX,GY(0,0,20),NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-e$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,22), GYL, STR, NRM, NRMAX, NGR, gDIV(22))
       CALL TXGRFRXS(5,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-i$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,23), GYL, STR, NRM, NRMAX, NGR, gDIV(23))
       CALL TXGRFRXS(6,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-b$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,24), GYL, STR, NRM, NRMAX, NGR, gDIV(24))
       CALL TXGRFRXS(7,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,21), GYL, STR, NRM, NRMAX, NGR, gDIV(21))
       CALL TXGRFRXS(8,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@D$-i$=(r)+D$-e$=(r)@'
       CALL TXGRFRXS(9,GX,GY(0,0,29),NRMAX,NGR,STR,MODE,IND)

       STR = '@F$-CDBM$=(r)@'
       CALL TXGRFRXS(10,GX,GY(0,0,31),NRMAX,NGR,STR,MODE,IND)     

       STR = '@G$-1$=h$+2$=(r)@'
       CALL TXGRFRXS(11,GX,GY(0,0,30),NRMAX,NGR,STR,MODE,IND)

       STR = '@S(r)@'
       CALL TXGRFRXS(12,GX,GY(0,0,32),NRMAX,NGR,STR,MODE,IND)

       STR = '@$#a$#(r)@'
       CALL TXGRFRXS(13,GX,GY(0,0,33),NRMAX,NGR,STR,MODE,IND)

       STR = '@$#k$#(r)@'
       CALL TXGRFRXS(14,GX,GY(0,0,34),NRMAX,NGR,STR,MODE,IND)

       STR = '@n$-i$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,37), GYL, STR, NRM, NRMAX, NGR, gDIV(1))
       CALL TXGRFRXS(15,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    CASE(-3)
       STR = '@LQm1@'
       CALL APPROPGY(MODEG, GY(0,0,38), GYL, STR, NRM, NRMAX, NGR, gDIV(9))
       CALL TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQm2@'
       CALL TXGRFRXS( 1,GX,GY(0,0,39),NRMAX,NGR,STR,MODE,IND)

       STR = '@LQm3@'
       CALL TXGRFRXS( 2,GX,GY(0,0,40),NRMAX,NGR,STR,MODE,IND)

       STR = '@LQm4@'
       CALL TXGRFRXS( 3,GX,GY(0,0,41),NRMAX,NGR,STR,MODE,IND)

       STR = '@LQm5@'
       CALL TXGRFRXS( 4,GX,GY(0,0,42),NRMAX,NGR,STR,MODE,IND)

       STR = '@LQe1@'
       CALL APPROPGY(MODEG, GY(0,0,43), GYL, STR, NRM, NRMAX, NGR, gDIV(1))
       CALL TXGRFRXS( 5,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQe2@'
       CALL APPROPGY(MODEG, GY(0,0,44), GYL, STR, NRM, NRMAX, NGR, gDIV(38))
       CALL TXGRFRXS( 6,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQe3@'
       CALL APPROPGY(MODEG, GY(0,0,45), GYL, STR, NRM, NRMAX, NGR, gDIV(39))
       CALL TXGRFRXS( 7,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQe4@'
       CALL APPROPGY(MODEG, GY(0,0,46), GYL, STR, NRM, NRMAX, NGR, gDIV(40))
       CALL TXGRFRXS( 8,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQe5@'
       CALL APPROPGY(MODEG, GY(0,0,47), GYL, STR, NRM, NRMAX, NGR, gDIV(41))
       CALL TXGRFRXS( 9,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQi1@'
       CALL APPROPGY(MODEG, GY(0,0,48), GYL, STR, NRM, NRMAX, NGR, gDIV(37))
       CALL TXGRFRXS(10,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQi2@'
       CALL APPROPGY(MODEG, GY(0,0,49), GYL, STR, NRM, NRMAX, NGR, gDIV(42))
       CALL TXGRFRXS(11,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQi3@'
       CALL APPROPGY(MODEG, GY(0,0,50), GYL, STR, NRM, NRMAX, NGR, gDIV(43))
       CALL TXGRFRXS(12,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQi4@'
       CALL APPROPGY(MODEG, GY(0,0,51), GYL, STR, NRM, NRMAX, NGR, gDIV(44))
       CALL TXGRFRXS(13,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQi5@'
       CALL APPROPGY(MODEG, GY(0,0,52), GYL, STR, NRM, NRMAX, NGR, gDIV(45))
       CALL TXGRFRXS(14,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQb1@'
       CALL APPROPGY(MODEG, GY(0,0,53), GYL, STR, NRM, NRMAX, NGR, gDIV(12))
       CALL TXGRFRXS(15,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    CASE(-4)
       STR = '@LQb3@'
       CALL APPROPGY(MODEG, GY(0,0,54), GYL, STR, NRM, NRMAX, NGR, gDIV(46))
       CALL TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQb4@'
       CALL APPROPGY(MODEG, GY(0,0,55), GYL, STR, NRM, NRMAX, NGR, gDIV(47))
       CALL TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQn1@'
       CALL APPROPGY(MODEG, GY(0,0,56), GYL, STR, NRM, NRMAX, NGR, gDIV(35))
       CALL TXGRFRXS( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)     

       STR = '@LQn2@'
       CALL APPROPGY(MODEG, GY(0,0,57), GYL, STR, NRM, NRMAX, NGR, gDIV(36))
       CALL TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    CASE(-5)
       STR = '@rMue@'
       CALL TXGRFRXS( 0,GX,GY(0,0,58),NRMAX,NGR,STR,MODE,IND)

       STR = '@rMui@'
       CALL TXGRFRXS( 1,GX,GY(0,0,59),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuei@'
       CALL TXGRFRXS( 2,GX,GY(0,0,60),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuL@'
       CALL TXGRFRXS( 3,GX,GY(0,0,61),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNube@'
       CALL TXGRFRXS( 4,GX,GY(0,0,62),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNubi@'
       CALL TXGRFRXS( 5,GX,GY(0,0,63),NRMAX,NGR,STR,MODE,IND)

       STR = '@FWthe@'
       CALL TXGRFRXS( 6,GX,GY(0,0,64),NRMAX,NGR,STR,MODE,IND)

       STR = '@FWthi@'
       CALL TXGRFRXS( 7,GX,GY(0,0,65),NRMAX,NGR,STR,MODE,IND)

       STR = '@WPM@'
       CALL TXGRFRXS( 8,GX,GY(0,0,66),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuTei@'
       CALL TXGRFRXS( 9,GX,GY(0,0,67),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNu0e@'
       CALL TXGRFRXS(10,GX,GY(0,0,68),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNu0i@'
       CALL TXGRFRXS(11,GX,GY(0,0,69),NRMAX,NGR,STR,MODE,IND)

       STR = '@Chie@'
       CALL TXGRFRXS(12,GX,GY(0,0,70),NRMAX,NGR,STR,MODE,IND)

       STR = '@Chii@'
       CALL TXGRFRXS(13,GX,GY(0,0,71),NRMAX,NGR,STR,MODE,IND)

       STR = '@D01@'
       CALL TXGRFRXS(14,GX,GY(0,0,72),NRMAX,NGR,STR,MODE,IND)

       STR = '@D02@'
       CALL TXGRFRXS(15,GX,GY(0,0,73),NRMAX,NGR,STR,MODE,IND)

    CASE(-6)
       STR = '@rNueNC@'
       CALL TXGRFRXS( 0,GX,GY(0,0,74),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuiNC@'
       CALL TXGRFRXS( 1,GX,GY(0,0,75),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNueHL@'
       CALL TXGRFRXS( 2,GX,GY(0,0,76),NRMAX,NGR,STR,MODE,IND)     

       STR = '@rNuiHL@'
       CALL TXGRFRXS( 3,GX,GY(0,0,77),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuiCX@'
       CALL TXGRFRXS( 4,GX,GY(0,0,78),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuION@'
       CALL TXGRFRXS( 5,GX,GY(0,0,79),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuB@'
       CALL TXGRFRXS( 6,GX,GY(0,0,80),NRMAX,NGR,STR,MODE,IND)     

       STR = '@SiLC@'
       CALL TXGRFRXS( 7,GX,GY(0,0,81),NRMAX,NGR,STR,MODE,IND)

       STR = '@SiLCth@'
       CALL TXGRFRXS( 8,GX,GY(0,0,82),NRMAX,NGR,STR,MODE,IND)

       STR = '@SiLCph@'
       CALL TXGRFRXS( 9,GX,GY(0,0,83),NRMAX,NGR,STR,MODE,IND)     

       STR = '@PRFe@'
       CALL TXGRFRXS(10,GX,GY(0,0,84),NRMAX,NGR,STR,MODE,IND)

       STR = '@PRFi@'
       CALL TXGRFRXS(11,GX,GY(0,0,85),NRMAX,NGR,STR,MODE,IND)

       STR = '@SNB@'
       CALL TXGRFRXS(12,GX,GY(0,0,86),NRMAX,NGR,STR,MODE,IND)

    CASE DEFAULT
       WRITE(6,*) 'Unknown NGYR: NGYR = ',NGYR
    END SELECT

    CALL PAGEE

    SELECT CASE(NGYR)
    CASE(-1)
       NGYR = -2
       CYCLE
    CASE(-2)
       NGYR = 0
       EXIT
    CASE(-3)
       NGYR = -4
       CYCLE
    CASE(-4)
       NGYR = 0
       EXIT
    CASE(-5)
       NGYR = -6
       CYCLE
    CASE(-6)
       NGYR = 0
       EXIT
    CASE DEFAULT
       EXIT
    END SELECT

    END DO

    RETURN
  END SUBROUTINE TXGRFR

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GX, GY
  !
  !***************************************************************

  SUBROUTINE TXGRFRX(K, GXL, GYL, NRMAX1, NGMAX, STR, MODE, IND)

    use commons, only : RA, RB, NRM
    INTEGER, INTENT(IN) :: K, NRMAX1, NGMAX, MODE, IND
    REAL, DIMENSION(0:NRM), INTENT(IN) :: GXL
    REAL, DIMENSION(0:NRM,0:NGMAX), INTENT(IN) :: GYL
    character(len=*), INTENT(IN) :: STR
    REAL :: GXMAX
    REAL, DIMENSION(4) :: GPXY

    GPXY(1) =  3.0 + 12.5 * MOD(K,2)
    GPXY(2) = 12.5 + 12.5 * MOD(K,2)
    GPXY(3) = 10.5 -  8.5 * REAL(K/2)
    GPXY(4) = 17.0 -  8.5 * REAL(K/2)
    GXMAX=REAL(RB/RA)
    CALL TXGRAF(GPXY, GXL, GYL, NRM+1, NRMAX1+1, NGMAX+1, &
         &            0.0, GXMAX, STR, MODE, IND)

    RETURN
  END SUBROUTINE TXGRFRX

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GX, GY (small version)
  !
  !***************************************************************

  SUBROUTINE TXGRFRXS(K, GXL, GYL, NRMAX1, NGMAX, STR, MODE, IND)

    use commons, only : RA, RB, NRM
    INTEGER, INTENT(IN) :: K, NRMAX1, NGMAX, MODE, IND
    REAL, DIMENSION(0:NRM), INTENT(IN) :: GXL
    REAL, DIMENSION(0:NRM,0:NGMAX), INTENT(IN) :: GYL
    character(len=*), INTENT(IN) :: STR
    REAL :: GXMAX
    REAL, DIMENSION(4) :: GPXY

    GPXY(1) =   2.0 + 6.1  * MOD(K,4)
    GPXY(2) =   6.7 + 6.1  * MOD(K,4)
    GPXY(3) = 13.75 - 4.25 * REAL(K/4)
    GPXY(4) = 17.0  - 4.25 * REAL(K/4)
    GXMAX=REAL(RB/RA)
    CALL TXGRAF(GPXY, GXL, GYL, NRM+1, NRMAX1+1, NGMAX+1, &
         &            0.0, GXMAX, STR, MODE, IND)

    RETURN
  END SUBROUTINE TXGRFRXS

  !***************************************************************
  !
  !   Write graph of GTX, GTY
  !
  !***************************************************************

  SUBROUTINE TXGRFT(NGYT,MODE)

    use commons
    INTEGER, INTENT(IN) :: NGYT, MODE
    INTEGER :: IND
    REAL :: gDIVL
    character(len=40) ::  STR
    REAL, DIMENSION(0:NGTM,NGYTM) :: GTYL

    IF (NGT <= 1) THEN
       WRITE(6,*) 'GT', NGYT, ' has no data'
       RETURN
    END IF

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF

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

    SELECT CASE(NGYT)
    CASE(1)
       STR = '@Ne(0),LAVE(0),(0.24),(0.6)@'
       CALL APPROPGY(MODEG, GTY(0,1), GTYL, STR, NGTM, NGT, 4, gDIV(1))
       CALL TXGRFTX(0, GTX, GTYL, NGTM, NGT, 4, STR,MODE, IND)

       STR = '@Z*Ni+Z*Nb-Ne(0)@'
       CALL APPROPGY(MODEG, GTY(0,5), GTYL, STR, NGTM, NGT, 1, gDIV(5-3))
       CALL TXGRFTX(1, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       STR = '@Er(b/2)@'
       CALL APPROPGY(MODEG, GTY(0,12), GTYL, STR, NGTM, NGT, 1, gDIV(12-3))
       CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       CALL TXWPGR

    CASE(2)
       STR = '@Uer(b/2)@'
       CALL TXGRFTX(0, GTX, GTY(0,6), NGTM, NGT, 1, STR,MODE, IND)

       STR = '@UeTheta(b/2)@'
       CALL APPROPGY(MODEG, GTY(0,7), GTYL, STR, NGTM, NGT, 1, gDIV(7-3))
       CALL TXGRFTX(1, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       STR = '@UePhi(0)@'
       CALL APPROPGY(MODEG, GTY(0,8), GTYL, STR, NGTM, NGT, 1, gDIV(8-3))
       CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       CALL TXWPGR

    CASE(3)
       STR = '@Uir(b/2)@'
       CALL TXGRFTX(0, GTX, GTY(0,9), NGTM, NGT, 1, STR,MODE, IND)

       STR = '@UiTheta(b/2)@'
       CALL APPROPGY(MODEG, GTY(0,10), GTYL, STR, NGTM, NGT, 1, gDIV(10-3))
       CALL TXGRFTX(1, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       STR = '@UiPhi(0)@'
       CALL APPROPGY(MODEG, GTY(0,11), GTYL, STR, NGTM, NGT, 1, gDIV(11-3))
       CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       CALL TXWPGR

    CASE(4)
       STR = '@Q(0)@'
       CALL TXGRFTX(0, GTX, GTY(0,20), NGTM, NGT, 1, STR,MODE, IND)

       STR = '@Ephi(b)@'
       CALL TXGRFTX(1, GTX, GTY(0,14), NGTM, NGT, 1, STR,MODE, IND)

       STR = '@Jphi(0)@'
       CALL APPROPGY(MODEG, GTY(0,21), GTYL, STR, NGTM, NGT, 1, gDIV(21-3))
       CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       CALL TXWPGR

    CASE(5)
       STR = '@Nb(0)@'
       CALL APPROPGY(MODEG, GTY(0,15), GTYL, STR, NGTM, NGT, 1, gDIV(15-3))
       CALL TXGRFTX(0, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       STR = '@UbPhi(0)@'
       CALL APPROPGY(MODEG, GTY(0,16), GTYL, STR, NGTM, NGT, 1, gDIV(16-3))
       CALL TXGRFTX(1, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       STR = '@Bth(b/2)@'
       CALL TXGRFTX(2, GTX, GTY(0,13), NGTM, NGT, 1, STR,MODE, IND)

       CALL TXWPGR

    CASE(6)
       STR = '@JePhi(0)@'
       CALL APPROPGY(MODEG, GTY(0,22), GTYL, STR, NGTM, NGT, 1, gDIV(22))
       CALL TXGRFTX(0, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       STR = '@JiPhi(0)@'
       CALL APPROPGY(MODEG, GTY(0,23), GTYL, STR, NGTM, NGT, 1, gDIV(23))
       CALL TXGRFTX(1, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       STR = '@JbPhi(0)@'
       CALL APPROPGY(MODEG, GTY(0,24), GTYL, STR, NGTM, NGT, 1, gDIV(24))
       CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       CALL TXWPGR

    CASE(7)
       STR = '@Uer(b/2)@'
       CALL TXGRFTX(0, GTX, GTY(0,6), NGTM, NGT, 1, STR,MODE, IND)

       STR = '@UePerp(b/2)@'
       CALL APPROPGY(MODEG, GTY(0,25), GTYL, STR, NGTM, NGT, 1, gDIV(25))
       CALL TXGRFTX(1, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       STR = '@UePara(b/2)@'
       CALL APPROPGY(MODEG, GTY(0,26), GTYL, STR, NGTM, NGT, 1, gDIV(26))
       CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       CALL TXWPGR

    CASE(8)
       STR = '@Uir(b/2)@'
       CALL TXGRFTX(0, GTX, GTY(0,9), NGTM, NGT, 1, STR,MODE, IND)

       STR = '@UiPerp(b/2)@'
       CALL APPROPGY(MODEG, GTY(0,27), GTYL, STR, NGTM, NGT, 1, gDIV(27))
       CALL TXGRFTX(1, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       STR = '@UiPara(b/2)@'
       CALL APPROPGY(MODEG, GTY(0,28), GTYL, STR, NGTM, NGT, 1, gDIV(28))
       CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       CALL TXWPGR

    CASE(9)
       STR = '@Di(b/2)+De(b/2)@'
       CALL TXGRFTX(0, GTX, GTY(0,29), NGTM, NGT, 1, STR,MODE, IND)
       STR = '@rG1h2(b/2)@'
       CALL TXGRFTX(1, GTX, GTY(0,30), NGTM, NGT, 1, STR,MODE, IND)
       STR = '@FCDBM(b/2)@'
       CALL TXGRFTX(2, GTX, GTY(0,31), NGTM, NGT, 1, STR,MODE, IND)
       CALL TXWPGR

    CASE(10)
       STR = '@Te(0)@'
       CALL TXGRFTX(0, GTX, GTY(0,17), NGTM, NGT, 1, STR,MODE, IND)

       STR = '@Ti(0)@'
       CALL TXGRFTX(1, GTX, GTY(0,18), NGTM, NGT, 1, STR,MODE, IND)

       STR = '@N0(0)@'
       CALL APPROPGY(MODEG, GTY(0,19), GTYL, STR, NGTM, NGT, 1, gDIV(19-3))
       CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       CALL TXWPGR

    CASE(11)
       STR = '@S(b/2)@'
       CALL TXGRFTX(0, GTX, GTY(0,32), NGTM, NGT, 1, STR,MODE, IND)
       STR = '@Alpha(b/2)@'
       CALL TXGRFTX(1, GTX, GTY(0,33), NGTM, NGT, 1, STR,MODE, IND)
       STR = '@rKappa(b/2)@'
       CALL TXGRFTX(2, GTX, GTY(0,34), NGTM, NGT, 1, STR,MODE, IND)
       CALL TXWPGR

       CALL TXWPGR

    CASE(12)
       STR = '@Ne:total@'
       CALL TXGRFTX(0, GTX, GTY(0,35), NGTM, NGT, 1, STR,MODE, IND)
       STR = '@Gamma:total@'
       CALL TXGRFTX(1, GTX, GTY(0,36), NGTM, NGT, 1, STR,MODE, IND)
       STR = '@taup@'
       CALL TXGRFTX(2, GTX, GTY(0,37), NGTM, NGT, 1, STR,MODE, IND)
       CALL TXWPGR

    CASE(13)
       STR = '@SLOW N0(0)@'
       CALL APPROPGY(MODEG, GTY(0,38), GTYL, STR, NGTM, NGT, 1, gDIV(38-3))
       CALL TXGRFTX(0, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       STR = '@FAST N0(0)@'
       CALL APPROPGY(MODEG, GTY(0,39), GTYL, STR, NGTM, NGT, 1, gDIV(39-3))
       CALL TXGRFTX(1, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       STR = '@TOTAL N0(0)@'
       CALL APPROPGY(MODEG, GTY(0,19), GTYL, STR, NGTM, NGT, 1, gDIV(19-3))
       CALL TXGRFTX(2, GTX, GTYL, NGTM, NGT, 1, STR,MODE, IND)

       CALL TXWPGR
    CASE DEFAULT
       WRITE(6,*) 'Unknown NGYT: NGYT = ',NGYT
    END SELECT

    CALL PAGEE

    RETURN
  END SUBROUTINE TXGRFT

  !***************************************************************
  !
  !   Write graph of GVX, GVY
  !
  !***************************************************************

  SUBROUTINE TXGRFV(NGYV,MODE)

    use commons
    INTEGER, INTENT(IN) :: NGYV, MODE
    INTEGER :: IND
    REAL :: gDIVL
    character(len=40) ::  STR
    REAL, DIMENSION(0:NGVM,NGYVM) :: GVYL

    IF (NGVV <= 1) THEN
       WRITE(6,*) 'GV', NGYV, ' has no data'
       RETURN
    END IF

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF

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

    SELECT CASE(NGYV)
    CASE(1)
       STR = '@Te,Ti,<Te>,<Ti> [keV] vs t@'
       CALL TXGRFTX(0, GVX, GVY(0,1), NGVM, NGVV, 4, STR, MODE, IND)

       STR = '@PIN,POH,PNB,PRF,PNF [MW] vs t@'
       CALL TXGRFTX(1, GVX, GVY(0,5), NGVM, NGVV, 5, STR, MODE, IND)

       STR = '@IP,IOH,INB,IOH+INB [MA] vs t@'
       CALL TXGRFTX(2, GVX, GVY(0,10), NGVM, NGVV, 4, STR, MODE, IND)

       STR = '@POUT,PCX,PIE,PRL,PCON [MW] vs t@'
       CALL TXGRFTX(3, GVX, GVY(0,15), NGVM, NGVV, 5, STR, MODE, IND)

    CASE(2)
       STR = '@QF vs t@'
       CALL TXGRFTX(0, GVX, GVY(0,20), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@Te0,Ti0 [keV] vs t@'
       CALL TXGRFTX(1, GVX, GVY(0,21), NGVM, NGVV, 2, STR, MODE, IND)

       STR = '@<Te>,<Ti> [keV] vs t@'
       CALL TXGRFTX(2, GVX, GVY(0,23), NGVM, NGVV, 2, STR, MODE, IND)

       STR = '@Ne0,Ni0,<Ne>,<Ni> [10^20/m^3] vs t@'
       CALL TXGRFTX(3, GVX, GVY(0,25), NGVM, NGVV, 4, STR, MODE, IND)

    CASE(3)
       STR = '@WF,WB,Wi,We [MJ] vs t@'
       CALL TXGRFTX(0, GVX, GVY(0,29), NGVM, NGVV, 4, STR, MODE, IND)

       STR = '@TAUE1,TAUE2,TAUEP [s] vs t@'
       CALL TXGRFTX(1, GVX, GVY(0,33), NGVM, NGVV, 3, STR, MODE, IND)

       STR = '@BETAa,BETA0 vs t@'
       CALL TXGRFTX(2, GVX, GVY(0,36), NGVM, NGVV, 2, STR, MODE, IND)

       STR = '@BETAPa,BETAP0 vs t@'
       CALL TXGRFTX(3, GVX, GVY(0,38), NGVM, NGVV, 2, STR, MODE, IND)

    CASE(4)
       STR = '@VLOOP [V]@'
       CALL TXGRFTX(0, GVX, GVY(0,40), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@ALI@'
       CALL TXGRFTX(1, GVX, GVY(0,41), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@Q(0)@'
       CALL TXGRFTX(2, GVX, GVY(0,42), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@RQ1 [m]@'
       CALL TXGRFTX(3, GVX, GVY(0,43), NGVM, NGVV, 1, STR, MODE, IND)

    CASE(5)
       GOTO 100

    CASE(6)
       STR = '@Te,Ti,<Te>,<Ti> [keV] vs t@'
       CALL TXGRFVX(0, GVX, GVY(0,1), NGVM, NGVV, 4, STR, IND)

       STR = '@PIN,POH,PNB,PRF,PNF [MW] vs t@'
       CALL TXGRFVX(4, GVX, GVY(0,5), NGVM, NGVV, 5, STR, IND)

       STR = '@IP,IOH,INB,IOH+INB [MA] vs t@'
       CALL TXGRFVX(2, GVX, GVY(0,10), NGVM, NGVV, 4, STR, IND)

       STR = '@POUT,PCX,PIE,PRL,PCON [MW] vs t@'
       CALL TXGRFVX(6, GVX, GVY(0,15), NGVM, NGVV, 5, STR, IND)

       STR = '@QF vs t@'
       CALL TXGRFVX(1, GVX, GVY(0,20), NGVM, NGVV, 1, STR, IND)

       STR = '@Te0,Ti0 [keV] vs t@'
       CALL TXGRFVX(5, GVX, GVY(0,21), NGVM, NGVV, 2, STR, IND)

       STR = '@<Te>,<Ti> [keV] vs t@'
       CALL TXGRFVX(3, GVX, GVY(0,23), NGVM, NGVV, 2, STR, IND)

       STR = '@Ne0,Ni0,<Ne>,<Ni> [10~20/m^3] vs t@'
       CALL TXGRFVX(7, GVX, GVY(0,25), NGVM, NGVV, 4, STR, IND)

    CASE(7)
       STR = '@WF,WB,Wi,We [MJ] vs t@'
       CALL TXGRFVX(0, GVX, GVY(0,29), NGVM, NGVV, 4, STR, IND)

       STR = '@TAUE1,TAUE2,TAUEP [e] vs t@'
       CALL TXGRFVX(4, GVX, GVY(0,33), NGVM, NGVV, 3, STR, IND)

       STR = '@BETAa,BETA0 vs t@'
       CALL TXGRFVX(2, GVX, GVY(0,36), NGVM, NGVV, 2, STR, IND)

       STR = '@BETAPa,BETAP0 vs t@'
       CALL TXGRFVX(6, GVX, GVY(0,38), NGVM, NGVV, 2, STR, IND)

       STR = '@VLOOP [V]@'
       CALL TXGRFVX(1, GVX, GVY(0,40), NGVM, NGVV, 1, STR, IND)

       STR = '@ALI@'
       CALL TXGRFVX(5, GVX, GVY(0,41), NGVM, NGVV, 1, STR, IND)

       STR = '@Q(0)@'
       CALL TXGRFVX(3, GVX, GVY(0,42), NGVM, NGVV, 1, STR, IND)

       STR = '@RQ1 [m]@'
       CALL TXGRFVX(7, GVX, GVY(0,43), NGVM, NGVV, 1, STR, IND)

    CASE DEFAULT
       WRITE(6,*) 'Unknown NGYV: NGYV = ',NGYV
    END SELECT

100 CALL PAGEE

    RETURN
  END SUBROUTINE TXGRFV

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GTX, GTY
  !
  !***************************************************************

  SUBROUTINE TXGRFTX(K, GTXL, GTYL, NGTLM, NGTL, NG, STR, MODE, IND)

    INTEGER, INTENT(IN) :: K ,NGTLM, NGTL, NG
    REAL, DIMENSION(0:NGTLM), INTENT(IN) :: GTXL
    REAL, DIMENSION(0:NGTLM,NG), INTENT(IN) :: GTYL
    character(len=*), INTENT(IN) ::  STR
    INTEGER :: MODE, IND
    REAL, DIMENSION(4) :: GPXY

    GPXY(1) =  3.0 + 12.5 * MOD(K,2)
    GPXY(2) = 12.5 + 12.5 * MOD(K,2)
    GPXY(3) = 10.5 -  8.5 * REAL(K/2)
    GPXY(4) = 17.0 -  8.5 * REAL(K/2)
    CALL TXGRAF(GPXY, GTXL, GTYL, NGTLM+1, NGTL+1, NG, &
         &            GTXL(0), GTXL(NGTL), STR, MODE, IND)

    RETURN
  END SUBROUTINE TXGRFTX

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GVX, GVY
  !
  !***************************************************************

  SUBROUTINE TXGRFVX(K, GTXL, GTYL, NGTLM, NGTL, NG, STR, IND)

    INTEGER, INTENT(IN) :: K ,NGTLM, NGTL, NG
    REAL, DIMENSION(0:NGTLM), INTENT(IN) :: GTXL
    REAL, DIMENSION(0:NGTLM,NG), INTENT(IN) :: GTYL
    character(len=*), INTENT(IN) ::  STR
    INTEGER :: IND
    REAL, DIMENSION(4) :: GPXY

    GPXY(1) =  3.0 + 12.0 * MOD(K,2)
    GPXY(2) = 12.0 + 12.0 * MOD(K,2)
    GPXY(3) = 14.0 -  4.3 * REAL(K/2)
    GPXY(4) = 17.0 -  4.3 * REAL(K/2)
    CALL TXGRAF(GPXY, GTXL, GTYL, NGTLM+1, NGTL+1, NG, &
         &            GTXL(0), GTXL(NGTL), STR, 1, IND)

    RETURN
  END SUBROUTINE TXGRFVX

  !***************************************************************
  !
  !   Graph of each TERM
  !
  !***************************************************************

  SUBROUTINE TXGRFQ(NQ,ID)

    use commons
    use libraries, only : APITOS, APTTOS

    INTEGER, INTENT(IN) :: NQ, ID
    INTEGER :: NR, NC, NC1, NSTR, IND
    REAL :: GXMAX
    REAL, DIMENSION(0:NRM,NCM) :: GQY
    REAL, DIMENSION(4) :: GPXY
    REAL, DIMENSION(4,5) :: GPXYA
    character(len=80), DIMENSION(NQM) :: STRGQ
    character(len=80) :: STR

    !        Title
    DATA STRGQ /'E$-r$=','E$-$#q$#$=','E$-$#f$#$=','B$-$#q$#$=','B$-$#f$#$=',  &
         &      'n$-e$=','n$-e$=u$-er$=', &
         &      'n$-e$=u$-e$#q$#$=','n$-e$=u$-e$#f$#$=','n$-e$=T$-e$=', &
         &      'n$-i$=','n$-i$=u$-ir$=', &
         &      'n$-i$=u$-i$#q$#$=','n$-i$=u$-i$#f$#$=','n$-i$=T$-i$=', &
         &      'n$-b$=', &
         &      'n$-b$=u$-b$#q$#$=','n$-b$=u$-b$#f$#$=', &
         &      'Slow n$-0$=', 'Fast n$-0$='/

    !        Graph coordinate
    DATA GPXYA/ 2.5, 9.5, 10.5,17.0, &
         &      2.5, 9.5,  1.5, 8.0, &
         &      15.0,22.0,10.5,17.0, &
         &      15.0,22.0, 1.5, 8.0, &
         &      2.5 ,22.0, 1.5,17.5/

    GPXY(1) = GPXYA(1,ID)
    GPXY(2) = GPXYA(2,ID)
    GPXY(3) = GPXYA(3,ID)
    GPXY(4) = GPXYA(4,ID)

    DO NC = 1, NLCMAX(NQ)
       NR = 0
       NC1 = NLCR(NC,NQ,NR)
       IF(NC1 == 0) THEN
          GQY(NR,NC) = SNGL(PLC(NC,NQ,NR))
       ELSE
          GQY(NR,NC) = SNGL(  BLC(NC,NQ,NR) * X(NC1,NR  ) &
               &            + ALC(NC,NQ,NR) * X(NC1,NR+1) &
               &            + PLC(NC,NQ,NR))
       END IF
       DO NR = 1, NRMAX - 1
          NC1 = NLCR(NC,NQ,NR)
          IF(NC1 == 0) THEN
             GQY(NR,NC) = SNGL(PLC(NC,NQ,NR))
          ELSE
             GQY(NR,NC) = SNGL(  CLC(NC,NQ,NR) * X(NC1,NR-1) &
                  &            + BLC(NC,NQ,NR) * X(NC1,NR  ) &
                  &            + ALC(NC,NQ,NR) * X(NC1,NR+1) &
                  &            + PLC(NC,NQ,NR))
          END IF
       END DO
       
       NR = NRMAX
       NC1 = NLCR(NC,NQ,NR)
       IF(NC1 == 0) THEN
          GQY(NR,NC) = SNGL(PLC(NC,NQ,NR))
       ELSE
          GQY(NR,NC) = SNGL(  CLC(NC,NQ,NR) * X(NC1,NR-1) &
               &            + BLC(NC,NQ,NR) * X(NC1,NR  ) &
               &            + PLC(NC,NQ,NR))
       END IF
    END DO

    NSTR = 0
    CALL APITOS(STR, NSTR, NQ)
    CALL APTTOS(STR, NSTR, ': ')
    CALL APTTOS(STR, NSTR, STRGQ(NQ))

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF
    GXMAX=REAL(RB/RA)
    CALL TXGRAF(GPXY, GX, GQY, NRM+1, NRMAX+1, NLCMAX(NQ), &
         &      0.0, GXMAX, '@'//STR(1:NSTR)//'@', 3, IND)

    RETURN
  END SUBROUTINE TXGRFQ

  !***************************************************************
  !
  !   Write parameter to graphic screen
  !
  !***************************************************************

  SUBROUTINE TXWPGR

    use commons
    INTEGER :: IFNT

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
  END SUBROUTINE TXWPGR

  !***************************************************************
  !
  !   WPgr's Sub : write Double
  !
  !***************************************************************

  SUBROUTINE TXWPSD(STR, VAL)

    REAL(8), INTENT(IN) :: VAL
    character(len=*), INTENT(IN) :: STR

    CALL MOVE(GXM, GYM - GYS * NP)
    CALL TEXTX(STR)
    CALL NUMBD(VAL,'(1PD9.2)', 9)
    NP = NP + 1

    RETURN
  END SUBROUTINE TXWPSD

  !***************************************************************
  !
  !   WPgr's Sub : write Integer
  !
  !***************************************************************

  SUBROUTINE TXWPSI(STR, IVAL)

    INTEGER, INTENT(IN) :: IVAL
    character(len=*), INTENT(IN) :: STR

    CALL MOVE(GXM, GYM - GYS * NP)
    CALL TEXTX(STR)
    CALL TEXT(' = ', 3)
    CALL NUMBI(IVAL,'(I6)',6)
    NP = NP + 1

    RETURN
  END SUBROUTINE TXWPSI

  !***************************************************************
  !
  !   WPgr's Sub : write Strings
  !
  !***************************************************************

  SUBROUTINE TXWPSS(STR)

    character(len=*), INTENT(IN) :: STR

    CALL MOVE(GXM, GYM - GYS * NP)
    CALL TEXTX(STR)
    NP = NP + 1

    RETURN
  END SUBROUTINE TXWPSS

  !***************************************************************
  !
  !   Draw graph
  !
  !***************************************************************

  SUBROUTINE TXGRAF(GPXY, GX, GY, NXM, NXMAX, NGMAX, &
       &                  GXMIN, GXMAX, STR, MODE, IND)

    INTEGER, INTENT(IN) :: NXM, NXMAX, NGMAX ,MODE, IND
    REAL, INTENT(IN) :: GXMIN, GXMAX
    REAL, DIMENSION(4), INTENT(IN) :: GPXY
    REAL, DIMENSION(NXMAX), INTENT(IN) :: GX
    REAL, DIMENSION(NXM,NGMAX), INTENT(IN) :: GY
    character(len=*), INTENT(IN) :: STR

    INTEGER :: IFNT, NGV, NGULEN, ICL, IPAT, IMRK, ISTEP, NG
    REAL :: GX1, GX2, GY1, GY2, gSLEN, GSXMIN, GSXMAX, GXSTEP, &
         &        GYMIN, GYMAX, GSYMIN, GSYMAX, GYSTEP, GYORG,  &
         &        GMRK, GCHH, GXL, GYL
    INTEGER, DIMENSION(0:4) :: NLTYPE
    DATA NLTYPE/0,2,3,4,6/

    IF (MODE < 0 .OR. MODE > 3) THEN
       WRITE(6,*) '### ERROR(TXGRAF) : MODE = ', MODE
       RETURN
    END IF

    GX1=GPXY(1)
    GX2=GPXY(2)
    GY1=GPXY(3)
    GY2=GPXY(4)

    CALL INQFNT(IFNT)

    IF (IND == 0) THEN
       gSLEN = 1.0
    ELSE
       gSLEN = 0.2
    END IF

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
    IF (GYMAX > 0.0) THEN
       IF (GYMIN > 0.0) THEN
          GYMIN=0.0
       END IF
    ELSE
       GYMAX=0.0
    END IF
    CALL GQSCAL(GYMIN, GYMAX, GSYMIN, GSYMAX, GYSTEP)
    IF (GSYMIN > GYMIN) GSYMIN = GSYMIN - GYSTEP
    IF (GSYMAX < GYMAX) GSYMAX = GSYMAX + GYSTEP
    GYORG = 0.0

    CALL GDEFIN(GX1, GX2, GY1, GY2, GSXMIN, GSXMAX, GSYMIN, GSYMAX)
    CALL GFRAME
    CALL GSCALE(GSXMIN, GXSTEP, 0.0, 0.0, gSLEN, IND)
    IF (GXSTEP < 0.01) THEN
       NGV=NGULEN(GXSTEP*5)
       IF (NGV < 0) THEN
          NGV=NGV-3200
       ELSE
          NGV=NGV+3200
       END IF
       CALL GVALUE(GSXMIN, GXSTEP*5, 0.0, 0.0, NGV)
    ELSE
       NGV=NGULEN(GXSTEP*2)
       IF (NGV < 0) THEN
          NGV=NGV-3200
       ELSE
          NGV=NGV+3200
       END IF
       CALL GVALUE(GSXMIN, GXSTEP*2, 0.0, 0.0, NGV)
    END IF
    CALL GSCALE(0.0, 0.0, GYORG, GYSTEP, gSLEN, IND)
    IF (GSYMIN < 0.0 .AND. GSYMAX > 0.0) &
         &     CALL GSCALE(0.0, 0.0, 0.0, GSYMAX-GSYMIN,  2.0, 0)
    CALL GVALUE(0.0,0.0,GYORG,GYSTEP*2,NGULEN(GYSTEP*2))

    ! MODE = 0: Change Line Color (Last Color Fixed)

    IF (MODE == 0) THEN
       DO NG = 1, NGMAX
          ICL = 7 - MOD(NGMAX - NG, 5)
          CALL SETLIN(0, 1, ICL)
          CALL GPLOTP(GX, GY(1,NG), 1, NXMAX, 1, 0, 0, 0)
       END DO

       ! MODE = 1: Change Line Color and Style

    ELSE IF (MODE == 1) THEN
       !     IF (MODE == 0.OR.MODE == 1) THEN
       DO NG = 1, NGMAX
          ICL  = 7 - MOD(NG-1, 5)
          IPAT = NLTYPE(MOD(NG-1, 5))
          CALL SETLIN(0, 1, ICL)
          CALL GPLOTP(GX, GY(1,NG), 1, NXMAX, 1, 0, 0, IPAT)
       END DO

       ! MODE = 2: Change Line Color, Style and Mark
       ! MODE = 3: Change Line Color, Style and Mark (With Legend)

    ELSE IF (MODE == 2 .OR. MODE == 3) THEN
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
       END DO
       ! Legend
       IF (MODE == 3) THEN
          CALL SETCHS(0.25,0.0)
          GCHH = 0.25
          GXL = GX2 + GCHH
          GYL = GY2 - GCHH
          DO NG = 1, NGMAX
             ICL = 7 - MOD(NG - 1, 5)
             CALL SETLIN(0, 1, ICL)
             IPAT = (NG - 1) / 5
             CALL MOVEPT(GXL + GCHH * 3.0, GYL + GCHH / 2.0, IPAT)
             CALL DRAWPT(GXL + GCHH * 7.0, GYL + GCHH / 2.0)
             CALL MOVE(GXL, GYL)
             CALL NUMBI(NG, '(I2)', 2)
             IMRK = MOD(NG - 1, 5) + 1
             IF (IMRK /= 0) THEN
                CALL SETMKS(IMRK, GMRK)
                CALL MARK(GXL + GCHH * 5.0, GYL + GCHH / 2.0)
             END IF
             GYL = GYL - GCHH * 2.0
          END DO
       END IF
       ! End of Legend
       IMRK = 0
       GMRK = 0.2
       CALL SETMKS(IMRK, GMRK)
    END IF

    CALL SETFNT(IFNT)

    RETURN
  END SUBROUTINE TXGRAF

  !***************************************************************
  !
  !   Appropriate GY and STR for graphic routine
  !
  !***************************************************************

  SUBROUTINE APPROPGY(MODE, GIN, GOUT, STR, NXM, NXMAX, NYMAX, gDIV)
    use libraries, only : APSTOS, APRTOS

    INTEGER, INTENT(IN) :: MODE, NXM, NXMAX, NYMAX
    REAL, INTENT(IN) :: gDIV
    REAL, DIMENSION(0:NXM,0:NYMAX),INTENT(IN)  :: GIN
    REAL, DIMENSION(0:NXM,0:NYMAX),INTENT(OUT) :: GOUT
    character(len=*), INTENT(INOUT) :: STR

    INTEGER :: NSTR, POSAT
    REAL :: gDIVL

    gDIVL = 1.0
    IF(MODE == 2) gDIVL = gDIV

    ! Append gDIV to string for showing multiplication factor
    IF (gDIVL /= 1.0) THEN
       POSAT = INDEX(STR,'@',.TRUE.)
       IF (POSAT /= 0) STR = STR(1:POSAT-1)
       NSTR = LEN_TRIM(STR)+1
       CALL APSTOS(STR, NSTR, ' [', 2)
       CALL APRTOS(STR, NSTR, gDIVL, 'E0')
       CALL APSTOS(STR, NSTR, ']@', 2)
    END IF

    GOUT(0:NXMAX,0:NYMAX) = GIN(0:NXMAX,0:NYMAX) / gDIVL

    RETURN
  END SUBROUTINE APPROPGY

end module graphic
