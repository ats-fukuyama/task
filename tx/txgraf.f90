!     $Id$

module graphic
  implicit none
  private
  real :: GXM, GYM, GYS
  integer :: NP
  public :: TXGOUT, TX_GRAPH_SAVE, TXSTGR, TXSTGT, TXSTGV

  interface TXWPS
     module procedure TXWPSD
     module procedure TXWPSI
     module procedure TXWPSS
  end interface

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

    INTEGER :: MODE, NGPR, NGPT, NGPV, NQ, NQL, NGF, NGFMAX, I, IST, NGRT
    real(4), dimension(0:NRMAX,0:5,1:NGYRM) :: GYL
    character(len=5) :: STR, STR2
    character(len=1) :: KID1, KID2

    !     *** MENU ***

    MODE = MODEGL
    OUTER : DO
       WRITE(6,*) '# SELECT : Rn: Tn: Un: Vn: C: A,B: ', &
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
          SELECT CASE(KID2)
          CASE('A')
             CALL TXGRFR(-1,MODE)
          CASE('B')
             CALL TXGRFR(-3,MODE)
          CASE('C')
             CALL TXGRFR(-5,MODE)
          CASE DEFAULT
             READ(STR(2:5),*,IOSTAT=IST) NGPR
             IF (IST /= 0) THEN
                WRITE(6,*) '### ERROR : Invalid Command : ', STR
                CYCLE
             END IF
             IF (NGPR >= 0 .AND. NGPR <= NGPRM) THEN
                CALL TXGRFR(NGPR,MODE)
             END IF
          END SELECT

       CASE('V')
          SELECT CASE(KID2)
          CASE('A')
             DO NGPT = 1, NGPVM
                CALL TXGRFV(NGPT,MODE)
             END DO
          CASE('B')
             DO NGPT = 1, 6
                CALL TXGRFV(NGPT,MODE)
             END DO
             DO NGPT = 9, 10
                CALL TXGRFV(NGPT,MODE)
             END DO

          CASE DEFAULT
             READ(STR(2:5),*,IOSTAT=IST) NGPT
             IF (IST < 0) THEN
                WRITE(6,*) '### ERROR : Invalid Command : ', STR
                CYCLE
             END IF
             IF (NGPT >= 1 .AND. NGPT <= NGPVM) THEN
                CALL TXGRFV(NGPT,MODE)
             END IF
          END SELECT

       CASE('T')
          SELECT CASE(KID2)
          CASE('A')
             DO NGPV = 1, NGPTM
                CALL TXGRFT(NGPV,MODE)
             END DO

          CASE('B')
             DO NGPV = 1, 6
                CALL TXGRFT(NGPV,MODE)
             END DO
             DO NGPV = 9, 10
                CALL TXGRFT(NGPV,MODE)
             END DO

          CASE DEFAULT
             READ(STR(2:5),*,IOSTAT=IST) NGPV
             IF (IST < 0) THEN
                WRITE(6,*) '### ERROR : Invalid Command : ', STR
                CYCLE
             END IF
             IF (NGPV >= 1 .AND. NGPV <= NGPTM) THEN
                CALL TXGRFT(NGPV,MODE)
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
             CALL TXGRFR(NGPR,MODE)
          END DO

       CASE('B')
          DO NGPR = 1, 6
             CALL TXGRFR(NGPR,MODE)
          END DO
          DO NGPR = 9, 10
             CALL TXGRFR(NGPR,MODE)
          END DO

       CASE('C')
          CALL TXGRCP(MODE)

       CASE('M')
          DO
             WRITE(6,*) '## Number of files (Max = 5):'
             READ(5,*,IOSTAT=IST) NGFMAX
             IF (IST > 0) THEN
                CYCLE
             ELSEIF (IST < 0) THEN
                CYCLE OUTER
             ELSE
                EXIT
             END IF
          END DO
          IF(T_TX == 0.D0) THEN
             NGRT = 1
          ELSE
             NGRT = 0
             GYL(0:NRMAX,0,1:NGYRM) = GY(0:NRMAX,NGR,1:NGYRM)
          END IF
          NGR=0
          DO NGF=1,NGFMAX
             CALL TXGLOD
             GYL(0:NRMAX,NGF-NGRT,1:NGYRM) = GY(0:NRMAX,NGR,1:NGYRM)
             CALL TX_GRAPH_SAVE
          END DO
          NGR = NGFMAX-NGRT
          GY(0:NRMAX,0:NGR,1:NGYRM) = GYL(0:NRMAX,0:NGR,1:NGYRM)
          DO 
             DO
                WRITE(6,*) '## INPUT GRAPH NUMBER: A,B,RA,RB,RC,Rn,X:exit'
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
                   CALL TXGRFR(I,MODE)
                END DO
                !     Correspond to GMB
             ELSE IF (STR2(1:1) == 'B') THEN
                DO I = 1, 6
                   CALL TXGRFR(I,MODE)
                END DO
                DO I = 9, 10
                   CALL TXGRFR(I,MODE)
                END DO
                !     Correspond to GMC
             ELSE IF (STR2(1:2) == 'RA') THEN
                CALL TXGRFR(-1,MODE)
             ELSE IF (STR2(1:2) == 'RB') THEN
                CALL TXGRFR(-3,MODE)
             ELSE IF (STR2(1:2) == 'RC') THEN
                CALL TXGRFR(-5,MODE)
             ELSE IF (STR2(1:1) == 'X') THEN
                EXIT
             ELSE
                READ(STR2,'(I5)',IOSTAT=IST) NGPR
                IF (IST < 0) CYCLE
                IF      (NGPR == 0) THEN
                   CYCLE
                ELSE IF (NGPR >= 0 .AND. NGPR <= NGPRM) THEN
                   CALL TXGRFR(NGPR,MODE)
                END IF
             END IF
          END DO

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
         &                    /SQRT(BphV(0:NRMAX)**2 + BthV(0:NRMAX)**2))
    GY(0:NRMAX,NGR,26) = SNGL(( BthV(0:NRMAX)*UethV(0:NRMAX) &
         &                     +BphV(0:NRMAX)*UephV(0:NRMAX)) &
         &                    /SQRT(BphV(0:NRMAX)**2 + BthV(0:NRMAX)**2))
    GY(0:NRMAX,NGR,27) = SNGL(( BphV(0:NRMAX)*UithV(0:NRMAX) &
         &                     -BthV(0:NRMAX)*UiphV(0:NRMAX)) &
         &                    /SQRT(BphV(0:NRMAX)**2 + BthV(0:NRMAX)**2))
    GY(0:NRMAX,NGR,28) = SNGL(( BthV(0:NRMAX)*UithV(0:NRMAX) &
         &                     +BphV(0:NRMAX)*UiphV(0:NRMAX)) &
         &                    /SQRT(BphV(0:NRMAX)**2 + BthV(0:NRMAX)**2))
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
    GY(0:NRMAX,NGR,56) = SNGL(X(LQn1,0:NRMAX)*1.D20)
    GY(0:NRMAX,NGR,57) = SNGL(X(LQn2,0:NRMAX)*1.D20)

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
    GY(0:NRMAX,NGR,87) = SNGL(rNu0b(0:NRMAX))

    GY(0:NRMAX,NGR,88) = SNGL(PIE(0:NRMAX))
    GY(0:NRMAX,NGR,89) = SNGL(PCX(0:NRMAX))
    GY(0:NRMAX,NGR,90) = SNGL(POH(0:NRMAX))
    GY(0:NRMAX,NGR,91) = SNGL(PBr(0:NRMAX))

    GY(0:NRMAX,NGR,92) = SNGL(rNuLTe(0:NRMAX))
    GY(0:NRMAX,NGR,93) = SNGL(rNuLTi(0:NRMAX))
    GY(0:NRMAX,NGR,94) = SNGL(rNuAse(0:NRMAX))
    GY(0:NRMAX,NGR,95) = SNGL(rNuASi(0:NRMAX))

    GY(0:NRMAX,NGR,96) = SNGL(PNBe(0:NRMAX))
    GY(0:NRMAX,NGR,97) = SNGL(PNBi(0:NRMAX))
    GY(0:NRMAX,NGR,98) = SNGL(POHe(0:NRMAX))
    GY(0:NRMAX,NGR,99) = SNGL(POHi(0:NRMAX))
    GY(0:NRMAX,NGR,100) = SNGL(PEQe(0:NRMAX))
    GY(0:NRMAX,NGR,101) = SNGL(PEQi(0:NRMAX))
    GY(0:NRMAX,NGR,102) = SNGL(Vbedir(0:NRMAX))
    GY(0:NRMAX,NGR,103) = SNGL(Vbidir(0:NRMAX))

    RETURN
  END SUBROUTINE TXSTGR

  !***************************************************************
  !
  !   Store GVY
  !
  !***************************************************************

  SUBROUTINE TXSTGV(GTIME)

    use commons
    use physical_constants, only : AEE, PI
    use output_console, only : rLINEAVE
    use libraries, only : VALINT_SUB

    REAL, INTENT(IN) :: GTIME
    REAL(8) :: BthL, BphL, BBL, PNESUM1, PNESUM2

    IF (NGVV < NGVM) NGVV=NGVV+1

    GVX(NGVV) = GTIME

    GVY(NGVV,1)  = SNGL(PNeV(0) * 1.D20)
    GVY(NGVV,2)  = SNGL((PZ * PNiV(0) + PZ * PNbV(0) - PNeV(0)) * 1.D20)
    GVY(NGVV,3)  = SNGL(UerV(NRC))
    GVY(NGVV,4)  = SNGL(UethV(NRC))
    GVY(NGVV,5)  = SNGL(UephV(0))
    GVY(NGVV,6)  = SNGL(UirV(NRC))
    GVY(NGVV,7)  = SNGL(UithV(NRC))
    GVY(NGVV,8)  = SNGL(UiphV(0))
    GVY(NGVV,9)  = SNGL(ErV(NRC))
    GVY(NGVV,10) = SNGL(BthV(NRC))
    GVY(NGVV,11) = SNGL(EphV(NRMAX))
    GVY(NGVV,12) = SNGL(PNbV(0) * 1.D20)
    GVY(NGVV,13) = SNGL(UbphV(0))
    GVY(NGVV,14) = SNGL(PTeV(0))
    GVY(NGVV,15) = SNGL(PTiV(0))
    GVY(NGVV,16) = SNGL((PN01V(0) + PN02V(0)) * 1.D20)
    GVY(NGVV,17) = SNGL(EthV(NRC))
    GVY(NGVV,18) = SNGL(BphV(0))
    GVY(NGVV,19) = SNGL(UbthV(NRC))
    GVY(NGVV,20) = SNGL(Q(0))

    GVY(NGVV,21) = SNGL((-      AEE * PNeV(0) * UephV(0) &
         &               + PZ * AEE * PNiV(0) * UiphV(0) &
         &               + PZ * AEE * PNbV(0) * UbphV(0)) * 1.D20)
    GVY(NGVV,22) = SNGL( -      AEE * PNeV(0) * UephV(0)  * 1.D20)
    GVY(NGVV,23) = SNGL(   PZ * AEE * PNiV(0) * UiphV(0)  * 1.D20)
    GVY(NGVV,24) = SNGL(   PZ * AEE * PNbV(0) * UbphV(0)  * 1.D20)

    BthL = BthV(NRC)
    BphL = BphV(NRC)
    BBL = SQRT(BphL**2 + BthL**2)
    GVY(NGVV,25) = SNGL((+ BphL * UethV(NRC) - BthL * UephV(NRC)) / BBL)
    GVY(NGVV,26) = SNGL((+ BthL * UethV(NRC) + BphL * UephV(NRC)) / BBL)
    GVY(NGVV,27) = SNGL((+ BphL * UithV(NRC) - BthL * UiphV(NRC)) / BBL)
    GVY(NGVV,28) = SNGL((+ BthL * UithV(NRC) + BphL * UiphV(NRC)) / BBL)
    GVY(NGVV,29) = SNGL(Di(NRC)+De(NRC))
    GVY(NGVV,30) = SNGL(rG1h2(NRC))
    GVY(NGVV,31) = SNGL(FCDBM(NRC))
    GVY(NGVV,32) = SNGL(S(NRC))
    GVY(NGVV,33) = SNGL(Alpha(NRC))
    GVY(NGVV,34) = SNGL(rKappa(NRC))
    GVY(NGVV,35) = SNGL(PN01V(0) * 1.D20)
    GVY(NGVV,36) = SNGL(PN02V(0) * 1.D20)
    GVY(NGVV,37) = SNGL(PNiV (0) * 1.D20)

    GVY(NGVV,38)  = SNGL(rLINEAVE(0.D0))
    GVY(NGVV,39)  = SNGL(rLINEAVE(0.24D0))
    GVY(NGVV,40)  = SNGL(rLINEAVE(0.6D0))

    CALL VALINT_SUB(PNeV,NRA,PNESUM1)
    PNESUM1 = 2.D0*PI*RR*2.D0*PI*(PNESUM1 + PNeV(NRA)*R(NRA)*(RA-R(NRA)))
    CALL VALINT_SUB(PNeV(0:NRMAX)*rNuION(0:NRMAX),NRA,PNESUM2)
    PNESUM2 = 2.D0*PI*RR*2.D0*PI*(PNESUM2 + PNeV(NRA)*rNuION(NRA)*R(NRA)*(RA-R(NRA)))

    GVY(NGVV,41) = SNGL(PNESUM1)
    GVY(NGVV,42) = SNGL(PNESUM2)
    IF(NGVV == 0.OR.ABS(PNESUM2) <= 0.D0) THEN
       GVY(NGVV,43) = 0.0
    ELSE
       GVY(NGVV,43) = SNGL(PNESUM1/PNESUM2)
    END IF

    GVY(NGVV,44) = SNGL(Chie(NRC))
    GVY(NGVV,45) = SNGL(Chii(NRC))
    GVY(NGVV,46) = SNGL(PIE(NRC))
    GVY(NGVV,47) = SNGL(PCX(NRC))
    GVY(NGVV,48) = SNGL(POH(NRC))
    GVY(NGVV,49) = SNGL(PBr(NRC))

    RETURN
  END SUBROUTINE TXSTGV

  !***************************************************************
  !
  !   Store GTY
  !
  !***************************************************************

  SUBROUTINE TXSTGT(GTIME)

    use commons
    REAL, INTENT(IN) :: GTIME

    IF (NGT < NGTM) NGT=NGT+1

    GTX(NGT) = GTIME

    GTY(NGT,1)  = SNGL(TS0(1))
    GTY(NGT,2)  = SNGL(TS0(2))
    GTY(NGT,3)  = SNGL(TSAV(1))
    GTY(NGT,4)  = SNGL(TSAV(2))
    GTY(NGT,5)  = SNGL(PINT)
    GTY(NGT,6)  = SNGL(POHT)
    GTY(NGT,7)  = SNGL(PNBT)
    GTY(NGT,8)  = SNGL(PRFT)

    GTY(NGT,10) = SNGL(AJT)
    GTY(NGT,11) = SNGL(AJOHT)
    GTY(NGT,12) = SNGL(AJNBT)
    GTY(NGT,13) = SNGL(AJOHT+AJNBT)
    !  GTY(NGT,14) = SNGL(AJBST)
    GTY(NGT,15) = SNGL(POUT)
    GTY(NGT,16) = SNGL(PCXT)
    GTY(NGT,17) = SNGL(PIET)
    !  GTY(NGT,18) = SNGL(PRLT)
    !  GTY(NGT,19) = SNGL(PCONT)
    GTY(NGT,20) = SNGL(QF)
    GTY(NGT,21) = SNGL(TS0(1))
    GTY(NGT,22) = SNGL(TS0(2))
    GTY(NGT,23) = SNGL(TSAV(1))
    GTY(NGT,24) = SNGL(TSAV(2))
    GTY(NGT,25) = SNGL(ANS0(1))
    GTY(NGT,26) = SNGL(ANS0(2))
    GTY(NGT,27) = SNGL(ANSAV(1))
    GTY(NGT,28) = SNGL(ANSAV(2))
    GTY(NGT,29) = SNGL(WPT)
    GTY(NGT,30) = SNGL(WBULKT)

    GTY(NGT,31) = SNGL(WST(1))
    GTY(NGT,32) = SNGL(WST(2))

    GTY(NGT,33) = SNGL(TAUE1)
    GTY(NGT,34) = SNGL(TAUE2)
    GTY(NGT,35) = SNGL(TAUEP)
    GTY(NGT,36) = SNGL(BETAA) * 100.0
    GTY(NGT,37) = SNGL(BETA0) * 100.0
    GTY(NGT,38) = SNGL(BETAPA)
    GTY(NGT,39) = SNGL(BETAP0)
    GTY(NGT,40) = SNGL(VLOOP)
    GTY(NGT,41) = SNGL(ALI)
    GTY(NGT,42) = SNGL(Q(0))
    GTY(NGT,43) = SNGL(RQ1)

    GTY(NGT,44) = SNGL(ANF0(1))  * 100.0
    GTY(NGT,45) = SNGL(ANFAV(1)) * 100.0

    RETURN
  END SUBROUTINE TXSTGT

  !***************************************************************
  !
  !   Time evolution of radial profile
  !
  !***************************************************************

  SUBROUTINE TXGRFR(NGYRIN,MODE)

    use commons
    INTEGER, INTENT(IN) :: MODE
    INTEGER, INTENT(IN) :: NGYRIN
    INTEGER :: IND, NG, NR, NGYR, NE, IFNT
    REAL(4), DIMENSION(0:NRM,0:NGRM) :: GYL
    character(len=50) :: STR
    character(len=1) :: KSTR,KLABEL
    character(len=3) :: KEND
    real, dimension(4) :: GPX, GPY
    real :: GPXL, FACT, GYMAX

    NGYR = NGYRIN

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
    CALL INQFNT(IFNT)

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
       ! Draw frame of R
       GPX(1) =  2.5 ; GPY(1) = 10.5
       GPX(2) =  2.5 ; GPY(2) = 16.5
       GPX(3) = 24.5 ; GPY(3) = 16.5
       GPX(4) = 24.5 ; GPY(4) = 10.5
       CALL LINES2D(GPX,GPY,4)
       ! Draw lines
       FACT = (GPX(4) - GPX(1)) / R(NRMAX)
       GPXL = GPX(1)
       DO NE = 1, NEMAX-1
          GPXL = GPXL + REAL(H(NE)) * FACT
          IF(NE == NRA) CALL SETLIN(-1,-1,6)
          CALL LINE2D(GPXL,10.5,GPXL,16.5)
          IF(NE == NRA) CALL SETLIN(-1,-1,7)
       END DO
       WRITE(KSTR,'(I1)') 0
       WRITE(KEND,'(I3)') NRMAX
       KLABEL='r'
       CALL SETCHS(0.4,0.0)
       CALL GTEXT2D(GPX(1),GPY(1)-0.5,KSTR,1,2)
       CALL GTEXT2D(GPX(4),GPY(4)-0.5,KEND,3,2)
       CALL SETFNT(33)
       CALL SETCHS(0.6,0.0)
       CALL GTEXT2D(GPX(1)-1.0,0.5*(GPY(1)+GPY(2)),KLABEL,1,2)
       CALL SETFNT(IFNT)

       ! Draw frame of PSI
       GPX(1) =  2.5 ; GPY(1) =  2.0
       GPX(2) =  2.5 ; GPY(2) =  8.0
       GPX(3) = 24.5 ; GPY(3) =  8.0
       GPX(4) = 24.5 ; GPY(4) =  2.0
       CALL LINES2D(GPX,GPY,4)
       ! Draw lines
       FACT = (GPX(4) - GPX(1)) / PSI(NRMAX)
       GPXL = GPX(1)
       DO NE = 1, NEMAX-1
          GPXL = GPXL + REAL(HPSI(NE)) * FACT
          IF(NE == NRA) CALL SETLIN(-1,-1,6)
          CALL LINE2D(GPXL,2.0,GPXL,8.0)
          IF(NE == NRA) CALL SETLIN(-1,-1,7)
       END DO
       WRITE(KSTR,'(I1)') 0
       WRITE(KEND,'(I3)') NRMAX
       KLABEL='s'
       CALL SETCHS(0.4,0.0)
       CALL GTEXT2D(GPX(1),GPY(1)-0.5,KSTR,1,2)
       CALL GTEXT2D(GPX(4),GPY(4)-0.5,KEND,3,2)
       CALL SETFNT(33)
       CALL SETCHS(0.6,0.0)
       CALL GTEXT2D(GPX(1)-1.0,0.5*(GPY(1)+GPY(2)),KLABEL,1,2)
       CALL SETFNT(IFNT)
       CALL SETCHS(0.3,0.0)

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

       STR = '@$#c$#$-e$=(r)@'
       CALL TXGRFRX(2,GX,GY(0,0,70),NRMAX,NGR,STR,MODE,IND)

       CALL TXWPGR

    CASE(11)
       STR = '@s(r)@'
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
       CALL TXGRFRX(0,GX,GY(0,0,88),NRMAX,NGR,STR,MODE,IND)

       STR = '@PCX(r)@'
       CALL APPROPGY(MODEG, GY(0,0,89), GYL, STR, NRM, NRMAX, NGR, gDIV(88))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@POH(r)@'
       CALL APPROPGY(MODEG, GY(0,0,90), GYL, STR, NRM, NRMAX, NGR, gDIV(89))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@PBr(r)@'
       CALL TXGRFRX(3,GX,GY(0,0,91),NRMAX,NGR,STR,MODE,IND)
       !CALL TXWPGR

    CASE(15)
       STR = '@PNBe(r)@'
       CALL APPROPGY(MODEG, GY(0,0,96), GYL, STR, NRM, NRMAX, NGR, gDIV(96))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@POHe(r)@'
       CALL APPROPGY(MODEG, GY(0,0,98), GYL, STR, NRM, NRMAX, NGR, gDIV(98))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@PNBi(r)@'
       CALL APPROPGY(MODEG, GY(0,0,97), GYL, STR, NRM, NRMAX, NGR, gDIV(97))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@POHi(r)@'
       CALL APPROPGY(MODEG, GY(0,0,99), GYL, STR, NRM, NRMAX, NGR, gDIV(99))
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    CASE(16)
       STR = '@PEQe(r)@'
       CALL APPROPGY(MODEG, GY(0,0,100), GYL, STR, NRM, NRMAX, NGR, gDIV(100))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@PEQi(r)@'
       CALL APPROPGY(MODEG, GY(0,0,101), GYL, STR, NRM, NRMAX, NGR, gDIV(101))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@Vbedir(r)@'
       CALL APPROPGY(MODEG, GY(0,0,102), GYL, STR, NRM, NRMAX, NGR, gDIV(102))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@Vbidir(r)@'
       CALL APPROPGY(MODEG, GY(0,0,103), GYL, STR, NRM, NRMAX, NGR, gDIV(103))
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    CASE(17)
       DO NR = 0, NRMAX
          IF(GX(NR) >= 0.95) EXIT
       END DO
       NR = NR - 1

       STR = '@n$-e$=(r)@'
       CALL APPROPGY(MODEG, GY(NR:NRMAX,0,1), GYL, STR, NRM, NRMAX-NR, NGR, gDIV(1))
       GYMAX = MAXVAL(GYL(0,0:NGR))
       CALL TXGRFRX(0, GX(NR:NRMAX), GYL, NRMAX-NR, NGR, STR, MODE, IND, &
            &       GX(NR), GYMAX)

       STR = '@T$-e$=(r)@'
       CALL APPROPGY(MODEG, GY(NR:NRMAX,0,14), GYL, STR, NRM, NRMAX-NR, NGR, gDIV(14))
       GYMAX = MAXVAL(GYL(0,0:NGR))
       CALL TXGRFRX(1, GX(NR:NRMAX), GYL, NRMAX-NR, NGR, STR, MODE, IND, &
            &       GX(NR), GYMAX)

       STR = '@T$-i$=(r)@'
       CALL APPROPGY(MODEG, GY(NR:NRMAX,0,15), GYL, STR, NRM, NRMAX-NR, NGR, gDIV(15))
       GYMAX = MAXVAL(GYL(0,0:NGR))
       CALL TXGRFRX(2, GX(NR:NRMAX), GYL, NRMAX-NR, NGR, STR, MODE, IND, &
            &       GX(NR), GYMAX)

       CALL TXWPGR
    CASE(-1)
       STR = '@E$-r$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,9), GYL, STR, NRM, NRMAX, NGR, gDIV(9))
       CALL TXGRFRXS(0, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@E$-$#q$#$=(r)@'
       CALL TXGRFRXS(1,GX,GY(0,0,17),NRMAX,NGR,STR,MODE,IND)

       STR = '@E$-$#f$#$=(r)@'
       CALL TXGRFRXS(2,GX,GY(0,0,11),NRMAX,NGR,STR,MODE,IND)

       STR = '@B$-$#q$#$=(r)@'
       CALL TXGRFRXS(3,GX,GY(0,0,10),NRMAX  ,NGR,STR,MODE,IND)

       STR = '@n$-e$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,1), GYL, STR, NRM, NRMAX, NGR, gDIV(1))
       CALL TXGRFRXS(4,GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@T$-e$=(r)@'
       CALL TXGRFRXS(5,GX,GY(0,0,14),NRMAX,NGR,STR,MODE,IND)

       STR = '@T$-i$=(r)@'
       CALL TXGRFRXS(6,GX,GY(0,0,15),NRMAX,NGR,STR,MODE,IND)

       STR = '@B$-$#f$#$=(r)@'
       CALL TXGRFRXS(7,GX,GY(0,0,18),NRMAX,NGR,STR,MODE,IND)

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

       STR = '@u$-b$#q$#$=(r)@'
       CALL TXGRFRXS(2,GX,GY(0,0,19),NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-b$#f$#$=(r)@'
       CALL APPROPGY(MODEG, GY(0,0,13), GYL, STR, NRM, NRMAX, NGR, gDIV(13))
       CALL TXGRFRXS(3,GX,GYL,NRMAX  ,NGR,STR,MODE,IND)     

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

       STR = '@$#c$#$-e$=@'
       CALL TXGRFRXS( 2,GX,GY(0,0,70),NRMAX,NGR,STR,MODE,IND)

       STR = '@$#c$#$-i$=@'
       CALL TXGRFRXS( 3,GX,GY(0,0,71),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuL@'
       CALL TXGRFRXS( 4,GX,GY(0,0,61),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuLTe@'
       CALL TXGRFRXS( 5,GX,GY(0,0,92),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuLTi@'
       CALL TXGRFRXS( 6,GX,GY(0,0,93),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuei@'
       CALL APPROPGY(MODEG, GY(0,0,60), GYL, STR, NRM, NRMAX, NGR, gDIV(60))
       CALL TXGRFRXS( 7,GX,GYL      ,NRMAX,NGR,STR,MODE,IND)

       STR = '@FWthe@'
       CALL APPROPGY(MODEG, GY(0,0,64), GYL, STR, NRM, NRMAX, NGR, gDIV(64))
       CALL TXGRFRXS( 8,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@FWthi@'
       CALL APPROPGY(MODEG, GY(0,0,65), GYL, STR, NRM, NRMAX, NGR, gDIV(65))
       CALL TXGRFRXS( 9,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@WPM@'
       CALL TXGRFRXS(10,GX,GY(0,0,66),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuTei@'
       CALL TXGRFRXS(11,GX,GY(0,0,67),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNu0e@'
       CALL TXGRFRXS(12,GX,GY(0,0,68),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNu0i@'
       CALL TXGRFRXS(13,GX,GY(0,0,69),NRMAX,NGR,STR,MODE,IND)

       STR = '@D01@'
       CALL APPROPGY(MODEG, GY(0,0,72), GYL, STR, NRM, NRMAX, NGR, gDIV(72))
       CALL TXGRFRXS(14,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@D02@'
       CALL APPROPGY(MODEG, GY(0,0,73), GYL, STR, NRM, NRMAX, NGR, gDIV(73))
       CALL TXGRFRXS(15,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    CASE(-6)
       STR = '@rNueNC@'
       CALL APPROPGY(MODEG, GY(0,0,74), GYL, STR, NRM, NRMAX, NGR, gDIV(74))
       CALL TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuiNC@'
       CALL APPROPGY(MODEG, GY(0,0,75), GYL, STR, NRM, NRMAX, NGR, gDIV(75))
       CALL TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuiCX@'
       CALL TXGRFRXS( 2,GX,GY(0,0,78),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuION@'
       CALL TXGRFRXS( 3,GX,GY(0,0,80),NRMAX,NGR,STR,MODE,IND)     

       STR = '@rNuAse@'
       CALL TXGRFRXS( 4,GX,GY(0,0,94),NRMAX,NGR,STR,MODE,IND,GYMAX=8.0)

       STR = '@rNuAsi@'
       CALL TXGRFRXS( 5,GX,GY(0,0,95),NRMAX,NGR,STR,MODE,IND,GYMAX=8.0)

       STR = '@rNu0b@'
       CALL TXGRFRXS( 6,GX,GY(0,0,87),NRMAX,NGR,STR,MODE,IND)

       STR = '@SNB@'
       CALL TXGRFRXS( 7,GX,GY(0,0,84),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNube@'
       CALL TXGRFRXS( 8,GX,GY(0,0,62),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNubi@'
       CALL TXGRFRXS( 9,GX,GY(0,0,63),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuB@'
       CALL TXGRFRXS(10,GX,GY(0,0,79),NRMAX,NGR,STR,MODE,IND)

       STR = '@SiLC@'
       CALL TXGRFRXS(11,GX,GY(0,0,81),NRMAX,NGR,STR,MODE,IND)

       STR = '@SiLCth@'
       CALL TXGRFRXS(12,GX,GY(0,0,82),NRMAX,NGR,STR,MODE,IND)

       STR = '@SiLCph@'
       CALL TXGRFRXS(13,GX,GY(0,0,83),NRMAX,NGR,STR,MODE,IND)     

       STR = '@PRFe@'
       CALL TXGRFRXS(14,GX,GY(0,0,85),NRMAX,NGR,STR,MODE,IND)

       STR = '@PRFi@'
       CALL TXGRFRXS(15,GX,GY(0,0,86),NRMAX,NGR,STR,MODE,IND)

!!$       STR = '@rNueHL@'
!!$       CALL TXGRFRXS(14,GX,GY(0,0,76),NRMAX,NGR,STR,MODE,IND)     
!!$
!!$       STR = '@rNuiHL@'
!!$       CALL TXGRFRXS(15,GX,GY(0,0,77),NRMAX,NGR,STR,MODE,IND)

    CASE DEFAULT
       WRITE(6,*) 'Unknown NGYR: NGYR = ',NGYR
    END SELECT

    CALL PAGEE

    SELECT CASE(NGYR)
    CASE(-1)
       NGYR = -2  ; CYCLE
    CASE(-2)
       NGYR =  0  ; EXIT
    CASE(-3)
       NGYR = -4  ; CYCLE
    CASE(-4)
       NGYR =  0  ; EXIT
    CASE(-5)
       NGYR = -6  ; CYCLE
    CASE(-6)
       NGYR =  0  ; EXIT
    CASE DEFAULT
       EXIT
    END SELECT

    END DO

    RETURN
  END SUBROUTINE TXGRFR


  !***************************************************************
  !
  !   Comparison with radial profiles at one slice time
  !
  !***************************************************************

  SUBROUTINE TXGRCP(MODE)

    use commons
    integer, intent(in) :: MODE
    character(len=50) :: STR
    integer :: IND, IFNT, NR
    real, dimension(0:NRM,1:3) :: GYL, GYL2

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF

    CALL PAGES
    CALL SETCHS(0.3, 0.0)
    CALL SETLIN(0, 1, 7)
    CALL INQFNT(IFNT)

    CALL MOVE(2.0,17.7)
    CALL TEXT('[G', 2)
    CALL NUMBI(1, '(I2)', 2)
    CALL TEXT(']  ', 3)
    CALL TEXT('FROM', 4)
    CALL NUMBR(GT(0), '(1PE9.2)', 9)
    CALL TEXT(' TO', 3)
    CALL NUMBR(GT(NGR), '(1PE9.2)', 9)
    CALL TEXT('  DT =', 6)
    CALL NUMBD(DT, '(1PD9.2)', 9)
    CALL TEXT('  NGRSTP = ', 11)
    CALL NUMBI(NGRSTP,'(I4)',4)

    DO NR = 0, NRMAX
       GYL(NR,1) = GLOG(ETA1(NR),1.D-10,1.D0)
!       GYL(NR,2) = GLOG(ETA2(NR),1.D-10,1.D0)
       GYL(NR,2) = GLOG(ETA3(NR),1.D-10,1.D0)
!       write(6,*) R(NR)/RA,AJBS3(NR)/AJBS1(NR)
    END DO

    STR = '@LOG: ETA@'
    CALL TXGRFRS(0, GX, GYL, NRMAX, 2, STR, MODE, IND, 1)
!!$    CALL TXGRFRS(0, GX, GYL, NRA, 2, STR, MODE, IND, 1)

    GYL(0:NRMAX,1) = SNGL(AJBS1(0:NRMAX))
!    GYL(0:NRMAX,2) = SNGL(AJBS2(0:NRMAX))
    GYL(0:NRMAX,2) = SNGL(AJBS3(0:NRMAX))

    STR = '@AJBS@'
    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRM, NRMAX, 2-1, gDIV(22))
    CALL TXGRFRS(1, GX, GYL2, NRMAX, 2, STR, MODE, IND, 0)
!!$    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRM, NRA, 2-1, gDIV(22))
!!$    CALL TXGRFRS(1, GX, GYL2, NRA, 2, STR, MODE, IND, 0)

    CALL PAGEE

  END SUBROUTINE TXGRCP

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GX, GY
  !
  !***************************************************************

  SUBROUTINE TXGRFRX(K, GXL, GYL, NRMAX, NGMAX, STR, MODE, IND, GXMIN, GYMAX, ILOGIN)

    use commons, only : NRM, RA, RB
    INTEGER, INTENT(IN) :: K, NRMAX, NGMAX, MODE, IND
    REAL, DIMENSION(0:NRM), INTENT(IN) :: GXL
    REAL, DIMENSION(0:NRM,1:NGMAX+1), INTENT(IN) :: GYL
    real, intent(in), optional :: GXMIN, GYMAX
    integer, intent(in), optional :: ILOGIN
    character(len=*), INTENT(IN) :: STR
    integer :: ILOG
    REAL :: GXMAX, GXMINL
    REAL, DIMENSION(4) :: GPXY

    GPXY(1) =  3.0 + 12.5 * MOD(K,2)
    GPXY(2) = 12.5 + 12.5 * MOD(K,2)
    GPXY(3) = 10.5 -  8.5 * REAL(K/2)
    GPXY(4) = 17.0 -  8.5 * REAL(K/2)
    GXMAX=REAL(RB/RA)

    IF(PRESENT(GXMIN)) THEN
       GXMINL = GXMIN
    ELSE
       GXMINL = 0.0
    END IF
    IF(PRESENT(ILOGIN)) THEN
       ILOG = ILOGIN
    ELSE
       ILOG = 0
    END IF

    IF(PRESENT(GYMAX)) THEN
       CALL TXGRAF(GPXY, GXL, GYL, NRM+1, NRMAX+1, NGMAX+1, &
            &            GXMINL, GXMAX, STR, 0.3, MODE, IND, ILOG, GYMAX)
    ELSE
       CALL TXGRAF(GPXY, GXL, GYL, NRM+1, NRMAX+1, NGMAX+1, &
            &            GXMINL, GXMAX, STR, 0.3, MODE, IND, ILOG)
    END IF

    RETURN
  END SUBROUTINE TXGRFRX

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GX, GY (small version)
  !
  !***************************************************************

  SUBROUTINE TXGRFRXS(K, GXL, GYL, NRMAX, NGMAX, STR, MODE, IND, GYMAX, ILOGIN)

    use commons, only : NRM, RA, RB
    INTEGER, INTENT(IN) :: K, NRMAX, NGMAX, MODE, IND
    REAL, DIMENSION(0:NRM), INTENT(IN) :: GXL
    REAL, DIMENSION(0:NRM,1:NGMAX+1), INTENT(IN) :: GYL
    character(len=*), INTENT(IN) :: STR
    real, intent(in), optional :: GYMAX
    integer, intent(in), optional :: ILOGIN
    integer :: ILOG
    REAL :: GXMAX
    REAL, DIMENSION(4) :: GPXY

    GPXY(1) =   2.0 + 6.1  * MOD(K,4)
    GPXY(2) =   6.7 + 6.1  * MOD(K,4)
    GPXY(3) = 13.75 - 4.25 * REAL(K/4)
    GPXY(4) = 17.0  - 4.25 * REAL(K/4)
    GXMAX=REAL(RB/RA)

    IF(PRESENT(ILOGIN)) THEN
       ILOG = ILOGIN
    ELSE
       ILOG = 0
    END IF

    IF(PRESENT(GYMAX)) THEN
       CALL TXGRAF(GPXY, GXL, GYL, NRM+1, NRMAX+1, NGMAX+1, &
            &            0.0, GXMAX, STR, 0.26, MODE, IND, ILOG, GYMAX)
    ELSE
       CALL TXGRAF(GPXY, GXL, GYL, NRM+1, NRMAX+1, NGMAX+1, &
            &            0.0, GXMAX, STR, 0.26, MODE, IND, ILOG)
    END IF

    RETURN
  END SUBROUTINE TXGRFRXS


  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GX, GY at one slice time
  !
  !***************************************************************

  SUBROUTINE TXGRFRS(K, GXL, GYL, NRMAX, NGMAX, STR, MODE, IND, ILOGIN)

    use commons, only : NRM, RA, RB
    INTEGER, INTENT(IN) :: K, NRMAX, NGMAX, MODE, IND
    REAL, DIMENSION(0:NRM), INTENT(IN) :: GXL
    REAL, DIMENSION(0:NRM,1:NGMAX), INTENT(IN) :: GYL
    character(len=*), INTENT(IN) :: STR
    integer, intent(in), optional :: ILOGIN
    integer :: ILOG
    REAL :: GXMAX
    REAL, DIMENSION(4) :: GPXY

    GPXY(1) =  3.0 + 12.5 * MOD(K,2)
    GPXY(2) = 12.5 + 12.5 * MOD(K,2)
    GPXY(3) = 10.5 -  8.5 * REAL(K/2)
    GPXY(4) = 17.0 -  8.5 * REAL(K/2)
    GXMAX=REAL(RB/RA)

    IF(PRESENT(ILOGIN)) THEN
       ILOG = ILOGIN
    ELSE
       ILOG = 0
    END IF

    CALL TXGRAF(GPXY, GXL, GYL, NRM+1, NRMAX+1, NGMAX, &
         &            0.0, GXMAX, STR, 0.3, MODE, IND, ILOG)

    RETURN
  END SUBROUTINE TXGRFRS

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
    character(len=50) ::  STR
    REAL, DIMENSION(0:NGVM,1:4) :: GVYL

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
       STR = '@n$-e$=(0),LAVE(0),(0.24),(0.6)@'
       CALL APPROPGY(MODEG, GVY(0,38), GVYL, STR, NGVM, NGVV, 3, gDIV(1), GVY(0, 1))
       CALL TXGRFVX(0, GVX, GVYL, NGVM, NGVV, 4, STR, MODE, IND)

       STR = '@Z*n$-i$=+Z*n$-b$=-n$-e$=(0)@'
       CALL APPROPGY(MODEG, GVY(0, 2), GVYL, STR, NGVM, NGVV, 1, gDIV(2))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@E$-r$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0, 9), GVYL, STR, NGVM, NGVV, 1, gDIV(9))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(2)
       STR = '@u$-er$=(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0, 3), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-e$#q$#$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0, 4), GVYL, STR, NGVM, NGVV, 1, gDIV(4))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-e$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0, 5), GVYL, STR, NGVM, NGVV, 1, gDIV(5))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(3)
       STR = '@u$-ir$=(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0, 6), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-i$#q$#$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0, 7), GVYL, STR, NGVM, NGVV, 1, gDIV(7))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-i$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0, 8), GVYL, STR, NGVM, NGVV, 1, gDIV(8))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(4)
       STR = '@q(0)@'
       CALL TXGRFVX(0, GVX, GVY(0,20), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@E$-$#f$#$=(b)@'
       CALL TXGRFVX(1, GVX, GVY(0,11), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@j$-$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,21), GVYL, STR, NGVM, NGVV, 1, gDIV(21))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(5)
       STR = '@n$-b$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,12), GVYL, STR, NGVM, NGVV, 1, gDIV(12))
       CALL TXGRFVX(0, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-b$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,13), GVYL, STR, NGVM, NGVV, 1, gDIV(13))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@B$-$#q$#$=(a/2)@'
       CALL TXGRFVX(2, GVX, GVY(0,10), NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(6)
       STR = '@j$-e$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,22), GVYL, STR, NGVM, NGVV, 1, gDIV(22))
       CALL TXGRFVX(0, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@j$-i$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,23), GVYL, STR, NGVM, NGVV, 1, gDIV(23))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@j$-b$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,24), GVYL, STR, NGVM, NGVV, 1, gDIV(24))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(7)
       STR = '@u$-er$=(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0, 3), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-e$#$/136$#$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0,25), GVYL, STR, NGVM, NGVV, 1, gDIV(25))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-e//$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0,26), GVYL, STR, NGVM, NGVV, 1, gDIV(26))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(8)
       STR = '@u$-ir$=(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0, 6), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-i$#$/136$#$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0,27), GVYL, STR, NGVM, NGVV, 1, gDIV(27))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-i//$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0,28), GVYL, STR, NGVM, NGVV, 1, gDIV(28))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(9)
       STR = '@D$-i$=(a/2)+D$-e$=(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0,29), NGVM, NGVV, 1, STR, MODE, IND)
       STR = '@G$-1$=h$+2$=(a/2)@'
       CALL TXGRFVX(1, GVX, GVY(0,30), NGVM, NGVV, 1, STR, MODE, IND)
       STR = '@F$-CDBM$=(a/2)@'
       CALL TXGRFVX(2, GVX, GVY(0,31), NGVM, NGVV, 1, STR, MODE, IND)
       CALL TXWPGR

    CASE(10)
       STR = '@T$-e$=(0)@'
       CALL TXGRFVX(0, GVX, GVY(0,14), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@T$-i$=(0)@'
       CALL TXGRFVX(1, GVX, GVY(0,14), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@$#c$#$-e$=(0)@'
       CALL TXGRFVX(2, GVX, GVY(0,44), NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(11)
       STR = '@s(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0,32), NGVM, NGVV, 1, STR, MODE, IND)
       STR = '@$#a$#(a/2)@'
       CALL TXGRFVX(1, GVX, GVY(0,33), NGVM, NGVV, 1, STR, MODE, IND)
       STR = '@$#k$#(a/2)@'
       CALL TXGRFVX(2, GVX, GVY(0,34), NGVM, NGVV, 1, STR, MODE, IND)
       CALL TXWPGR

       CALL TXWPGR

    CASE(12)
       STR = '@E$-$#q$#$=(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0,17), NGVM, NGVV, 1, STR, MODE, IND)
       STR = '@B$-$#f$#$=(0)@'
       CALL TXGRFVX(1, GVX, GVY(0,18), NGVM, NGVV, 1, STR, MODE, IND)
       STR = '@u$-b$#q$#$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0,19), GVYL, STR, NGVM, NGVV, 1, gDIV(19))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)
       CALL TXWPGR

       CALL TXWPGR

    CASE(13)
       STR = '@SLOW N$-0$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,35), GVYL, STR, NGVM, NGVV, 1, gDIV(35))
       CALL TXGRFVX(0, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@FAST N$-0$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,36), GVYL, STR, NGVM, NGVV, 1, gDIV(36))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@TOTAL N$-0$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,16), GVYL, STR, NGVM, NGVV, 1, gDIV(16))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(14)
       STR = '@PIE(r)@'
       CALL TXGRFVX(0, GVX, GVY(0,46), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@PCX(r)@'
       CALL APPROPGY(MODEG, GVY(0,47), GVYL, STR, NGVM, NGVV, 1, gDIV(88))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@POH(r)@'
       CALL APPROPGY(MODEG, GVY(0,48), GVYL, STR, NGVM, NGVV, 1, gDIV(89))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@PBr(r)@'
       CALL TXGRFVX(3, GVX, GVY(0,49), NGVM, NGVV, 1, STR, MODE, IND)
       !CALL TXWPGR

    CASE(15)
       STR = '@n$-e$=:total@'
       CALL TXGRFVX(0, GVX, GVY(0,41), NGVM, NGVV, 1, STR, MODE, IND)
       STR = '@$#g$#:total@'
       CALL TXGRFVX(1, GVX, GVY(0,42), NGVM, NGVV, 1, STR, MODE, IND)
       STR = '@$#t$#$-p$=@'
       CALL TXGRFVX(2, GVX, GVY(0,43), NGVM, NGVV, 1, STR, MODE, IND)
       CALL TXWPGR

    CASE DEFAULT
       WRITE(6,*) 'Unknown NGYV: NGYV = ',NGYV
    END SELECT

    CALL PAGEE

    RETURN
  END SUBROUTINE TXGRFV

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
    character(len=50) ::  STR
    REAL, DIMENSION(0:NGTM,1:NGYTM) :: GTYL

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
       STR = '@Te,Ti,<Te>,<Ti> [keV] vs t@'
       CALL TXGRFVX(0, GTX, GTY(0,1), NGTM, NGT, 4, STR, MODE, IND)

       STR = '@PIN,POH,PNB,PRF,PNF [MW] vs t@'
       CALL TXGRFVX(1, GTX, GTY(0,5), NGTM, NGT, 5, STR, MODE, IND)

       STR = '@IP,IOH,INB,IOH+INB [MA] vs t@'
       CALL TXGRFVX(2, GTX, GTY(0,10), NGTM, NGT, 4, STR, MODE, IND)

       STR = '@POUT,PCX,PIE,PRL,PCON [MW] vs t@'
       CALL TXGRFVX(3, GTX, GTY(0,15), NGTM, NGT, 5, STR, MODE, IND)

    CASE(2)
       STR = '@QF vs t@'
       CALL TXGRFVX(0, GTX, GTY(0,20), NGTM, NGT, 1, STR, MODE, IND)

       STR = '@Te0,Ti0 [keV] vs t@'
       CALL TXGRFVX(1, GTX, GTY(0,21), NGTM, NGT, 2, STR, MODE, IND)

       STR = '@<Te>,<Ti> [keV] vs t@'
       CALL TXGRFVX(2, GTX, GTY(0,23), NGTM, NGT, 2, STR, MODE, IND)

       STR = '@Ne0,Ni0,<Ne>,<Ni> [10$+20$=/m$+3$=] vs t@'
       CALL TXGRFVX(3, GTX, GTY(0,25), NGTM, NGT, 4, STR, MODE, IND)

    CASE(3)
       STR = '@WF,WB,Wi,We [MJ] vs t@'
       CALL TXGRFVX(0, GTX, GTY(0,29), NGTM, NGT, 4, STR, MODE, IND)

       STR = '@TAUE1,TAUE2,TAUEP [s] vs t@'
       CALL TXGRFVX(1, GTX, GTY(0,33), NGTM, NGT, 3, STR, MODE, IND)

       STR = '@BETAa,BETA0 [%] vs t@'
       CALL TXGRFVX(2, GTX, GTY(0,36), NGTM, NGT, 2, STR, MODE, IND)

       STR = '@BETAPa,BETAP0 vs t@'
       CALL TXGRFVX(3, GTX, GTY(0,38), NGTM, NGT, 2, STR, MODE, IND)

    CASE(4)
       STR = '@VLOOP [V]@'
       CALL TXGRFVX(0, GTX, GTY(0,40), NGTM, NGT, 1, STR, MODE, IND)

       STR = '@ALI@'
       CALL TXGRFVX(1, GTX, GTY(0,41), NGTM, NGT, 1, STR, MODE, IND)

       STR = '@Q(0)@'
       CALL TXGRFVX(2, GTX, GTY(0,42), NGTM, NGT, 1, STR, MODE, IND)

       STR = '@RQ1 [m]@'
       CALL TXGRFVX(3, GTX, GTY(0,43), NGTM, NGT, 1, STR, MODE, IND)

    CASE(5)
       ! Do nothing.

    CASE(6)
       STR = '@Te,Ti,<Te>,<Ti> [keV] vs t@'
       CALL TXGRFTX(0, GTX, GTY(0, 1), NGTM, NGT, 4, STR, IND)

       STR = '@IP,IOH,INB,IOH+INB [MA] vs t@'
       CALL TXGRFTX(2, GTX, GTY(0,10), NGTM, NGT, 4, STR, IND)

       STR = '@PIN,POH,PNB,PRF,PNF [MW] vs t@'
       CALL TXGRFTX(4, GTX, GTY(0, 5), NGTM, NGT, 5, STR, IND)

       STR = '@POUT,PCX,PIE,PRL,PCON [MW] vs t@'
       CALL TXGRFTX(6, GTX, GTY(0,15), NGTM, NGT, 5, STR, IND)

       STR = '@QF vs t@'
       CALL TXGRFTX(1, GTX, GTY(0,20), NGTM, NGT, 1, STR, IND)

       STR = '@<Te>,<Ti> [keV] vs t@'
       CALL TXGRFTX(3, GTX, GTY(0,23), NGTM, NGT, 2, STR, IND)

       STR = '@Te0,Ti0 [keV] vs t@'
       CALL TXGRFTX(5, GTX, GTY(0,21), NGTM, NGT, 2, STR, IND)

       STR = '@Ne0,Ni0,<Ne>,<Ni> [10$+20$=/m$+3$=] vs t@'
       CALL TXGRFTX(7, GTX, GTY(0,25), NGTM, NGT, 4, STR, IND)

    CASE(7)
       STR = '@WF,WB,Wi,We [MJ] vs t@'
       CALL TXGRFTX(0, GTX, GTY(0,29), NGTM, NGT, 4, STR, IND)

       STR = '@BETAa,BETA0 [%] vs t@'
       CALL TXGRFTX(2, GTX, GTY(0,36), NGTM, NGT, 2, STR, IND)

       STR = '@TAUE1,TAUE2,TAUEP [s] vs t@'
       CALL TXGRFTX(4, GTX, GTY(0,33), NGTM, NGT, 3, STR, IND)

       STR = '@BETAPa,BETAP0 vs t@'
       CALL TXGRFTX(6, GTX, GTY(0,38), NGTM, NGT, 2, STR, IND)

       STR = '@VLOOP [V]@'
       CALL TXGRFTX(1, GTX, GTY(0,40), NGTM, NGT, 1, STR, IND)

       STR = '@Q(0)@'
       CALL TXGRFTX(3, GTX, GTY(0,42), NGTM, NGT, 1, STR, IND)

       STR = '@ALI@'
       CALL TXGRFTX(5, GTX, GTY(0,41), NGTM, NGT, 1, STR, IND)

       STR = '@RQ1 [m]@'
       CALL TXGRFTX(7, GTX, GTY(0,43), NGTM, NGT, 1, STR, IND)

    CASE(8)
       STR = '@Nb0,<Nb> [10$+18$=/m$+3$=] vs t@'
       CALL TXGRFTX(0, GTX, GTY(0,44), NGTM, NGT, 2, STR, IND)

    CASE DEFAULT
       WRITE(6,*) 'Unknown NGYT: NGYT = ',NGYT
    END SELECT

    CALL PAGEE

    RETURN
  END SUBROUTINE TXGRFT

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GVX, GVY
  !
  !***************************************************************

  SUBROUTINE TXGRFVX(K, GTXL, GTYL, NGTM, NGTL, NG, STR, MODE, IND)

    INTEGER, INTENT(IN) :: K, NGTM, NGTL, NG
    REAL, DIMENSION(0:NGTM), INTENT(IN) :: GTXL
    REAL, DIMENSION(0:NGTM,1:NG), INTENT(IN) :: GTYL
    character(len=*), INTENT(IN) ::  STR
    INTEGER :: MODE, IND
    REAL, DIMENSION(4) :: GPXY

    GPXY(1) =  3.0 + 12.5 * MOD(K,2)
    GPXY(2) = 12.5 + 12.5 * MOD(K,2)
    GPXY(3) = 10.5 -  8.5 * REAL(K/2)
    GPXY(4) = 17.0 -  8.5 * REAL(K/2)
    CALL TXGRAF(GPXY, GTXL, GTYL, NGTM+1, NGTL+1, NG, &
         &            GTXL(0), GTXL(NGTL), STR, 0.3, MODE, IND, 0)

    RETURN
  END SUBROUTINE TXGRFVX

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GTX, GTY
  !
  !***************************************************************

  SUBROUTINE TXGRFTX(K, GTXL, GTYL, NGTM, NGTL, NG, STR, IND)

    INTEGER, INTENT(IN) :: K, NGTM, NGTL, NG
    REAL, DIMENSION(0:NGTM), INTENT(IN) :: GTXL
    REAL, DIMENSION(0:NGTM,1:NG), INTENT(IN) :: GTYL
    character(len=*), INTENT(IN) ::  STR
    INTEGER :: IND
    REAL, DIMENSION(4) :: GPXY

    GPXY(1) =  3.0 + 12.0 * MOD(K,2)
    GPXY(2) = 12.0 + 12.0 * MOD(K,2)
    GPXY(3) = 14.0 -  4.3 * REAL(K/2)
    GPXY(4) = 17.0 -  4.3 * REAL(K/2)
    CALL TXGRAF(GPXY, GTXL, GTYL, NGTM+1, NGTL+1, NG, &
         &            GTXL(0), GTXL(NGTL), STR, 0.3, 1, IND, 0)

    RETURN
  END SUBROUTINE TXGRFTX

  !***************************************************************
  !
  !   Graph of each TERM
  !
  !***************************************************************

  SUBROUTINE TXGRFQ(NQ,ID)

    use commons
    use libraries, only : APTOS
    use physical_constants, only : EPS0, rMU0

    INTEGER, INTENT(IN) :: NQ, ID
    INTEGER :: NR, NC, NC1, NSTR, IND
    REAL :: GXMAX
    REAL, DIMENSION(0:NRM,NCM) :: GQY
    REAL, DIMENSION(4) :: GPXY
    REAL, DIMENSION(4,5) :: GPXYA
    character(len=80), DIMENSION(NQM) :: STRGQ
    character(len=80) :: STR

    !        Title
    DATA STRGQ /'$#f$#$',"A$-$#q$#$='","A$-$#f$#$='",'A$-$#f$#$=','A$-$#q$#$=',  &
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

!!$    GQY(0:10,1:NLCMAX(NQ))=100.0*GQY(0:10,1:NLCMAX(NQ))
!!$    GQY(25:NRMAX,1:NLCMAX(NQ))=100.0*GQY(25:NRMAX,1:NLCMAX(NQ))
!!$    do nc=1,nlcmax(nq)
!!$       write(6,*) nc,gqy(15,nc),gqy(nrmax,nc)
!!$    end do

    NSTR = 0
    CALL APTOS(STR,NSTR,NQ)
    CALL APTOS(STR,NSTR, ': ',2)
    CALL APTOS(STR, NSTR, STRGQ(NQ), LEN_TRIM(STRGQ(NQ)))

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF
    GXMAX=REAL(RB/RA)
    CALL TXGRAF(GPXY, GX, GQY, NRM+1, NRMAX+1, NLCMAX(NQ), &
         &      0.0, GXMAX, '@'//STR(1:NSTR)//'@', 0.3, 2, IND, 0)
!         &      0.0, GXMAX, '@'//STR(1:NSTR)//'@', 0.3, 4, IND, 0)

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
!    CALL SETFNT(32)
    CALL SETFNT(0)
    CALL SETCHS(0.3, 0.0)
    CALL SETLIN(0, 1, 7)

    GXM = 13.0 + 1.0
    GYM =  8.5 + 0.2
    GYS =  0.5
    NP  = 0

    CALL TXWPS('+'//SLID//'+')
    CALL TXWPS('@PNBCD @', PNBCD)
    CALL TXWPS('@NRMAX @', NRMAX)

    NP = NP + 1
    CALL TXWPS('@BB    @', BB)
    CALL TXWPS('@rIp   @', rIp)
    CALL TXWPS('@FSDFIX@', FSDFIX)
    CALL TXWPS('@FSCDBM@', FSCDBM)
    CALL TXWPS('@FSBOHM@', FSBOHM)
    CALL TXWPS('@FSPSCL@', FSPSCL)
    CALL TXWPS('@PROFD @', PROFD)
    CALL TXWPS('@FSCX  @', FSCX)
    CALL TXWPS('@FSLC  @', FSLC)
    CALL TXWPS('@FSNC  @', FSNC)
    CALL TXWPS('@FSLP  @', FSLP)
    CALL TXWPS('@FSION @', FSION)
    CALL TXWPS('@FSD0  @', FSD0)

    GXM = GXM + 0.35 * 17
    NP = 0
    CALL TXWPS('@PNBH  @', PNBH)
    CALL TXWPS('@PRFH  @', PRFH)
    CALL TXWPS('@De0   @', De0)
    CALL TXWPS('@Di0   @', Di0)
    CALL TXWPS('@rMue0 @', rMue0)
    CALL TXWPS('@rMui0 @', rMui0)
    CALL TXWPS('@Chie0 @', Chie0)
    CALL TXWPS('@Chii0 @', Chii0)
    CALL TXWPS('@WPM0  @', WPM0)
    CALL TXWPS('@PTe0  @', PTe0)
    CALL TXWPS('@PTea  @', PTea)
    CALL TXWPS('@PTi0  @', PTi0)
    CALL TXWPS('@PTia  @', PTia)
    CALL TXWPS('@V0    @', V0)
    CALL TXWPS('@rGamm0@', rGamm0)
    CALL TXWPS('@rGASPF@', rGASPF)
    CALL TXWPS('@Zeff  @', Zeff)

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
       &                  GXMIN, GXMAX, STR, FONT, MODE, IND, ILOG, GYMAX_IN)

    INTEGER, INTENT(IN) :: NXM, NXMAX, NGMAX, MODE, IND
    REAL, INTENT(IN) :: GXMIN, GXMAX, FONT
    REAL, DIMENSION(4), INTENT(IN) :: GPXY
    REAL, DIMENSION(1:NXMAX), INTENT(IN) :: GX
    REAL, DIMENSION(1:NXM,1:NGMAX), INTENT(IN) :: GY
    integer, intent(in) :: ILOG
    REAL, INTENT(IN), OPTIONAL :: GYMAX_IN
    character(len=*), INTENT(IN) :: STR

    INTEGER :: IFNT, NGV, NGULEN, ICL, IPAT, IMRK, ISTEP, NG
    REAL :: GX1, GX2, GY1, GY2, gSLEN, GSXMIN, GSXMAX, GXSTEP, &
         &        GYMIN, GYMAX, GSYMIN, GSYMAX, GYSTEP, GYORG,  &
         &        GMRK, GCHH, GXL, GYL
    INTEGER, DIMENSION(0:4) :: NLTYPE
    DATA NLTYPE/0,2,3,4,6/

    IF (MODE < 0 .OR. MODE > 4) THEN
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
    CALL SETCHS(FONT, 0.0)
    CALL SETLIN(0, 0, 7)
    CALL GTEXTX(GX1,GY2+0.2,STR,0)

    CALL SETFNT(44)
    CALL SETCHS(FONT, 0.0)
    CALL SETLIN(0, 0, 7)
    CALL SETLNW(0.035)

    CALL GQSCAL(GXMIN, GXMAX, GSXMIN, GSXMAX, GXSTEP)
    GSXMIN = GXMIN
    GSXMAX = GXMAX
    CALL GMNMX2(GY,NXM,1,NXMAX,1,1,NGMAX,1,GYMIN,GYMAX)
    IF(ILOG == 0) THEN
       IF (GYMAX > 0.0) THEN
          IF (GYMIN > 0.0) GYMIN=0.0
       ELSE
          GYMAX=0.0
       END IF
    END IF
    IF(PRESENT(GYMAX_IN)) GYMAX=GYMAX_IN
    CALL GQSCAL(GYMIN, GYMAX, GSYMIN, GSYMAX, GYSTEP)
    IF (GSYMIN > GYMIN) GSYMIN = GSYMIN - GYSTEP
    IF (GSYMAX < GYMAX) GSYMAX = GSYMAX + GYSTEP
    IF(GYMIN * GYMAX <= 0.0) THEN
       GYORG = 0.0
    ELSE
       GYORG = GSYMIN
    ENDIF

    CALL GDEFIN(GX1, GX2, GY1, GY2, GSXMIN, GSXMAX, GSYMIN, GSYMAX)
    CALL GFRAME
    CALL SETLNW( 0.017)
    CALL GSCALE(GSXMIN, GXSTEP, 0.0, 0.0, gSLEN, IND)
!!$    IF (GXSTEP < 0.01) THEN
!!$       NGV=NGULEN(GXSTEP*5)
!!$       IF (NGV < 0) THEN
!!$          NGV=NGV-3200
!!$       ELSE
!!$          NGV=NGV+3200
!!$       END IF
!!$       CALL GVALUE(GSXMIN, GXSTEP*5, 0.0, 0.0, NGV)
       NGV=NGULEN(GXSTEP*2)
       CALL GVALUE(GSXMIN, GXSTEP*2, 0.0, 0.0, NGV)
!!$    ELSE
!!$       NGV=NGULEN(GXSTEP*2)
!!$       IF (NGV < 0) THEN
!!$          NGV=NGV-3200
!!$       ELSE
!!$          NGV=NGV+3200
!!$       END IF
!!$       CALL GVALUE(GSXMIN, GXSTEP*2, 0.0, 0.0, NGV)
!!$    END IF

    ! Semi-Log or not

    IF(ILOG == 0) THEN
       CALL GSCALE(0.0, 0.0, GYORG, GYSTEP, gSLEN, IND)
       CALL SETLNW(-0.017)
       IF (GSYMIN < 0.0 .AND. GSYMAX > 0.0) &
            &     CALL GSCALE(0.0, 0.0, 0.0, GSYMAX-GSYMIN,  2.0, 0)
       CALL GVALUE(0.0,0.0,GYORG,GYSTEP*2,NGULEN(GYSTEP*2))
    ELSE
       CALL GSCALL(0.0, 0.0, GYORG, 4, gSLEN, IND)
       CALL SETLNW(-0.017)
       IF (GSYMIN < 0.0 .AND. GSYMAX > 0.0) &
            &     CALL GSCALL(0.0, 0.0, 0.0, 3,  2.0, 0)
       CALL GVALUL(0.0,0.0,GYORG,1,NGULEN(GYSTEP*2))
    END IF

    ! MODE = 0: Change Line Color (Last Color Fixed)

    SELECT CASE(MODE)
    CASE (0)
       DO NG = 1, NGMAX
          ICL = 7 - MOD(NGMAX - NG, 5)
          CALL SETLIN(0, 1, ICL)
          CALL GPLOTP(GX, GY(1,NG), 1, NXMAX, 1, 0, 0, 0)
       END DO

       ! MODE = 1: Change Line Color and Style
       ! MODE = 2: Change Line Color and Style (With Legend)

    CASE (1:2)
       IF (MODE == 1) THEN
          DO NG = 1, NGMAX
             ICL  = 7 - MOD(NG-1, 5)
             IPAT = NLTYPE(MOD(NG-1, 5))
             CALL SETLIN(0, 1, ICL)
             CALL GPLOTP(GX, GY(1,NG), 1, NXMAX, 1, 0, 0, IPAT)
          END DO
       ELSE
          DO NG = 1, NGMAX
             ICL = 7 - MOD(NG - 1, 5)
             ISTEP = NXMAX / 10
             IPAT  = (NG - 1) / 5
             CALL SETLIN(0, 1, ICL)
             CALL GPLOTP(GX, GY(1,NG), 1, NXMAX, 1, 0, ISTEP, IPAT)
          END DO
       END IF
       ! Legend
       IF (MODE == 2) THEN
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
             GYL = GYL - GCHH * 2.0
          END DO
       END IF
       ! End of Legend

       ! MODE = 3: Change Line Color, Style and Mark
       ! MODE = 4: Change Line Color, Style and Mark (With Legend)

    CASE (3:4)
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
       IF (MODE == 4) THEN
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
    END SELECT

    CALL SETFNT(IFNT)

    RETURN
  END SUBROUTINE TXGRAF

  !***************************************************************
  !
  !   Appropriate GY and STR for graphic routine
  !
  !***************************************************************

  SUBROUTINE APPROPGY(MODE, GIN, GOUT, STR, NXM, NXMAX, NYMAX, gDIV, GIN1)
    use libraries, only : APTOS

    INTEGER, INTENT(IN) :: MODE, NXM, NXMAX, NYMAX
    REAL, INTENT(IN) :: gDIV
    REAL, DIMENSION(0:NXM,0:NYMAX),INTENT(IN)  :: GIN
    REAL, DIMENSION(0:NXM,0:NYMAX),INTENT(OUT) :: GOUT
    REAL, DIMENSION(0:NXM),INTENT(IN), optional :: GIN1
    character(len=*), INTENT(INOUT) :: STR

    INTEGER :: NSTR, POSAT
    REAL :: gDIVL

    gDIVL = 1.0
    IF(MODE /= 0) gDIVL = gDIV

    ! Append gDIV to string for showing multiplication factor
    IF (gDIVL /= 1.0) THEN
       POSAT = INDEX(STR,'@',.TRUE.)
       IF (POSAT /= 0) STR = STR(1:POSAT-1)
       NSTR = LEN_TRIM(STR)+1
       CALL APTOS(STR, NSTR, ' [', 2)
       CALL APTOS(STR, NSTR, gDIVL, 'E0')
       CALL APTOS(STR, NSTR, ']@', 2)
    END IF

    if(present(GIN1)) then
       GOUT(0:NXMAX,0) = GIN1(0:NXMAX) / gDIVL
       GOUT(0:NXMAX,1:NYMAX) = GIN(0:NXMAX,0:NYMAX-1) / gDIVL
    else
       GOUT(0:NXMAX,0:NYMAX) = GIN(0:NXMAX,0:NYMAX) / gDIVL
    end if

    RETURN
  END SUBROUTINE APPROPGY

  !***********************************************************
  !
  !   CEILING FUNCTION FOR LOG10 PLOT
  !
  !***********************************************************

  REAL FUNCTION GLOG(X,XMIN,XMAX)

    implicit none
    real GUCLIP
    real(8), intent(in) :: X, XMIN, XMAX
    real(8) :: PLOG

    GLOG = GUCLIP(PLOG(X,XMIN,XMAX))

    RETURN
  END FUNCTION GLOG

end module graphic
