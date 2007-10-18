!***************************************************************
!
!   Command loop
!
!***************************************************************

module menu
  use libraries, only : TOUPPER
  use parameter_control
  implicit none
  private
  public :: TXMENU

contains

  SUBROUTINE TXMENU

    use physical_constants, only : rMU0
    use commons
    use main, only : TXEXEC
    use graphic, only : TX_GRAPH_SAVE, TXSTGR, TXGOUT
    use file_io, only : TXSAVE, TXLOAD
    use variables, only : TXCALV
    use init_prof, only : TXPROF, TXINIT
    use parameter_control, only : TXPARM_CHECK
    use results, only : TXSTAT
    INTEGER :: ICONT, MODE, I, IST
    character(len=80) :: LINE
    character(len=1)  :: KID, KID2

    IERR  = 0
    ICONT = 0
    IF(PNBH /= 0) THEN
       rMUb1 = 1.D0
       rMUb2 = rMU0
    END IF

    !  *** MENU ***

    DO
       IF(ICONT == 0) TMAX = DT*NTMAX

       WRITE(6,'(3(A,1PD12.4))') &
            &   ' ## TIME=',T_TX,'  DT=',DT,'  NEXT TIME =',TMAX
       WRITE(6,*) '## INPUT: ', &
            &   'R:RUN  C:CONT  P,V:PARM  G:GR  '// &
            &   'W:STAT  S:SAVE  L:LOAD  I:INIT  Q:QUIT'
       CALL GUFLSH

       CALL TXKLIN(LINE,KID,MODE)
       IF(MODE /= 1) THEN
          CALL TXPARM_CHECK
          TMAX=T_TX+DT*NTMAX
          IF(rMUb1 == rMU0 .and. (PNBHT1 /= 0.D0 .OR. PNBHT2 /= 0.D0 .OR. PNBHP /= 0.D0)) THEN
             rMUb1 = 1.D0
             rMUb2 = rMU0
             X(LQm4,0:NRMAX) = X(LQm4,0:NRMAX) / rMUb2
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
          T_TX = 0.D0
          TPRE = 0.D0
          IERR = 0
          ICONT = 1
          CALL TXPROF
          CALL TX_GRAPH_SAVE
          CALL TXEXEC
          TMAX=T_TX+DT*NTMAX
       CASE('C')
          IF (ICONT == 0) THEN
             WRITE(6,*) 'XX RUN or LOAD before CONTINUE !'
             CYCLE
          END IF
          NGR=-1
          CALL TXSTGR
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
          CALL TXSTAT
       CASE('S')
          CALL TXSAVE
       CASE('L')
          CALL TXLOAD
          IERR = 0
          NGR = -1
          ICONT = 1
          CALL TX_GRAPH_SAVE
       CASE('G')
          CALL TXGOUT
       CASE('N')
          I = NINT((DelR / (RB / NRMAX)) - 0.5D0)
          X(I,1) = X(I,1) + DelN
          X(I,2) = X(I,2) + DelN
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
       CASE('#')
          CONTINUE
       CASE DEFAULT
          WRITE(6,*) 'XX Unknown command'
       END SELECT
    END DO

    RETURN
  END SUBROUTINE TXMENU

  !***** INPUT KID or LINE *****
  !           MODE=0: LINE INPUT 
  !                1: KID INPUT
  !                2: PARM INPUT
  !                3: NEW PROMPT

  SUBROUTINE TXKLIN(LINE,KID,MODE)

    INTEGER :: ID, I, MODE, IST
    character(len=80) :: LINE
    character(len=1)  :: KID

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
end module menu
