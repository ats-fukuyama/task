!***************************************************************
!
!   Command loop
!
!***************************************************************

SUBROUTINE TXMENU

  use tx_commons, only : allocate_txcomm, deallocate_txcomm, &
       &              rMU0, IERR, rMUb1, rMUb2, TMAX, DT, NTMAX, T_TX, PNBHT1, &
       &              PNBHT2, PNBHP, PNBHex, X, LQm4, NRMAX, TPRE, PNBCD, &
       &              iflag_file, &
       &              DelRho, DelN, Rho, LQe1, LQi1, achg, ICONT, IRPIN
  use tx_main, only : TXEXEC
  use tx_graphic, only : TX_GRAPH_SAVE, TXSTGR, TXGOUT, allocate_txgraf, deallocate_txgraf, &
       &                 NGR, NGRM, GT, GY, gDIV
  use tx_variables, only : TXCALV
  use tx_parameter_control, only : TXPARM_CHECK, TXPARM, TXVIEW
  use tx_interface, only : TXKLIN, TOUPPER, TXLOAD
  use tx_ripple, only : ripple_input, ripple_spl, deallocate_ripple
  use tx_ntv, only : dealloc_ntv
  use eqread_mod, only : alloc_equ
  implicit none
  INTEGER(4) :: MODE, I, IST, ier, NR, IER_RP
  character(len=80) :: LINE
  character(len=1)  :: KID, KID2

  IERR  = 0
  ICONT = 0
  IRPIN = 0 ! Read ripple file if IRPIN /= 0
  IF((PNBHP + PNBHT1 + PNBHT2 + PNBHex) /= 0) THEN
     rMUb1 = 1.D0
     rMUb2 = rMU0
     gDIV(41) = 1.e3
  END IF

  !  *** MENU ***

  DO
     ier = 0
     IF(ICONT == 0) TMAX = DT*NTMAX

     WRITE(6,'(3(A,1PD12.4))') &
          &   ' ## TIME=',T_TX,'  DT=',DT,'  NEXT TIME =',TMAX
     WRITE(6,*) '## INPUT: ', &
          &   'R:RUN  C:CONT  P,V:PARM  G:GRAPH  '// &
          &   'W:STAT  S:SAVE  L:LOAD  I:INIT '
     WRITE(6,'(11X,A)') 'F,FR:FILE  N:PTRB  M:ITG  O:OUT  Q:QUIT'
     CALL GUFLSH

     CALL TXKLIN(LINE,KID,MODE)
     IF(MODE /= 1) THEN
        CALL TXPARM_CHECK
        TMAX=T_TX+DT*NTMAX
        IF(rMUb1 == rMU0 .and. &
        & (PNBHT1 /= 0.D0 .OR. PNBHT2 /= 0.D0 .OR. PNBHP /= 0.D0 .OR. PNBHex /= 0.D0)) THEN
           rMUb1 = 1.D0
           rMUb2 = rMU0
           gDIV(41) = 1.e3
           if(allocated(X)) X(0:NRMAX,LQm4) = X(0:NRMAX,LQm4) / rMUb2
        END IF
        CYCLE
     END IF

     SELECT CASE(KID)
     CASE('R')
        IF (ICONT /= 0) THEN
           WRITE(6,*) '# Would you like to restart? [y/N]'
           READ(5,'(A1)',IOSTAT=IST) KID
           IF(IST /= 0) CYCLE
           CALL TOUPPER(KID)
           IF(KID /= 'Y') CYCLE
        END IF
        call allocate_txcomm(ier, icont)
        call allocate_txgraf(ier, icont)
        if(ier /= 0) cycle
        T_TX = 0.D0
        TPRE = 0.D0
        IERR = 0
        ICONT = 1
        CALL TXPROF
        CALL TX_GRAPH_SAVE
        if(IRPIN /= 0) CALL ripple_spl
        CALL TXEXEC
        TMAX=T_TX+DT*NTMAX
     CASE('C')
        IF (ICONT == 0) THEN
           WRITE(6,*) 'XX RUN or LOAD before CONTINUE !'
           CYCLE
        END IF
        NGR=-1
        CALL TXSTGR(NGR,GT,GY,NGRM)
        if(IRPIN /= 0) CALL ripple_spl
        CALL TXEXEC
        TMAX=T_TX+DT*NTMAX
     CASE('P')
        CALL TXPARM(KID)
        IF(KID == 'Q') EXIT
     CASE('V')
        CALL TXVIEW
     CASE('I')
        CALL TXINIT
     CASE('Q')
        EXIT
     CASE('W')
        IF (ICONT == 0) THEN
           WRITE(6,*) 'XX RUN or LOAD before CONTINUE !'
           CYCLE
        END IF
        CALL TXSTAT
     CASE('S')
        CALL TXSAVE
     CASE('L')
        CALL TXLOAD(IER)
        IF(IER /= 0) CYCLE
        IERR = 0
        NGR = -1
        ICONT = 1
        CALL TX_GRAPH_SAVE
     CASE('G')
        CALL TXGOUT
     CASE('N')
        IF (ICONT == 0) THEN
           WRITE(6,*) 'XX RUN or LOAD before CONTINUE !'
           CYCLE
        END IF
        DO NR = 0, NRMAX-1
           IF(Rho(NR) <= DelRho .AND. Rho(NR+1) >= DelRho) THEN
              I = NR
              EXIT
           END IF
        END DO
        X(LQe1,I) = X(LQe1,I) + DelN
        X(LQi1,I) = X(LQi1,I) + DelN * achg(2)
        CALL TXCALV(X)
     CASE('B')
        KID2=LINE(2:2)
        CALL TOUPPER(KID2)
        SELECT CASE(KID2)
        CASE('P')
           PNBCD =  1.D0
        CASE('0')
           PNBCD =  0.D0
        CASE('M')
           PNBCD = -1.D0
        CASE DEFAULT
           WRITE(6,*) 'XX Unknown beam command'
        END SELECT
     CASE('F')
        KID2=LINE(2:2)
        CALL TOUPPER(KID2)
        SELECT CASE(KID2)
        CASE('R')
           CALL ripple_input(IER_RP)
           if(IER_RP == 0) IRPIN = 1
        CASE DEFAULT
           CALL ascii_input

           if(iflag_file == 1) then
              if(rMUb2 == 1.D0) then
                 rMUb1 = 1.D0
                 rMUb2 = rMU0
                 gDIV(41) = 1.e3

                 if(allocated(X)) X(LQm4,0:NRMAX) = X(LQm4,0:NRMAX) / rMUb2
                 IF(T_TX /= 0.D0) CALL TXCALV(X)
              end if
           end if
        END SELECT
     CASE('O')
        IF (ICONT == 0) THEN
           WRITE(6,*) 'XX RUN or LOAD before CONTINUE !'
           CYCLE
        END IF
        CALL outfile
     CASE('M')
        IF (ICONT == 0) THEN
           WRITE(6,*) 'XX RUN or LOAD before CONTINUE !'
           CYCLE
        END IF
        CALL ITG_growthrate
     CASE('#')
        CONTINUE
     CASE DEFAULT
        WRITE(6,*) 'XX Unknown command'
     END SELECT
  END DO

  ! *** deallocate dynamic arrays before termination ***

  call deallocate_ripple
  call dealloc_ntv
  if(allocated(rho)) call deallocate_txcomm
  if(allocated(GY) ) call deallocate_txgraf
  call alloc_equ(-1)

  RETURN
END SUBROUTINE TXMENU

!***** INPUT KID or LINE *****
!           MODE=0: LINE INPUT 
!                1: KID INPUT
!                2: PARM INPUT
!                3: NEW PROMPT

SUBROUTINE TXKLIN(LINE,KID,MODE)

  use tx_parameter_control, only : TXPARL
  use tx_interface, only : TOUPPER
  implicit none

  integer(4), intent(out) :: MODE
  character(len=80), intent(out) :: LINE
  character(len=1), intent(out) :: KID
  integer(4) :: ID, I, IST

  !  ----- read line input input -----

  READ(5,'(A80)',IOSTAT=IST) LINE

  IF(IST == 0) THEN
     !  ----- parameter input -----

     ID=0
     DO I=1,80
        IF(LINE(I:I) == '=') ID=1
     END DO
     IF(ID == 1) THEN
        CALL TXPARL(LINE)
        KID=' '
        MODE=2
        RETURN
     END IF

     !  ----- command input -----

     KID=LINE(1:1)
     CALL TOUPPER(KID)
     IF(KID >= 'A'.AND.KID <= 'Z') THEN
        MODE=1
        RETURN
     END IF

     !  ----- line input -----

     KID=' '
     MODE=0

  ELSE IF(IST < 0) THEN

     !  ----- input end -----

     KID='Q'
     MODE=1

  ELSE IF(IST > 0) THEN

     !  ----- input error -----

     WRITE(6,*) 'XX INPUT ERROR !'
     MODE=3

  END IF

  RETURN
END SUBROUTINE TXKLIN
