!***************************************************************
!
!   Command loop
!
!***************************************************************

SUBROUTINE TXMENU

  INCLUDE 'txcomm.inc'

  INTEGER :: ICONT, MODE, I
  CHARACTER(80) :: LINE, KID*1, KID2*1

  IERR=0
  ICONT = 0

!  *** MENU ***

  DO
     WRITE(6,'(3(A,1PD12.4))') &
          &   ' ## TIME=',TIME,'  DT=',DT,'  TMAX =',TMAX
     WRITE(6,*) '## INPUT: ', &
          &   'R:RUN  C:CONT  P,V:PARM  G:GRAPH  '// &
          &   'S:SAVE  L:LOAD  I,N,Bn: Q:QUIT'
     CALL GUFLSH

     CALL TXKLIN(LINE,KID,MODE)
     IF(MODE /= 1) CYCLE

     IF      (KID == 'R') THEN
        TIME = 0.D0
        TPRE = 0.D0
        IERR = 0
        ICONT = 1
        CALL TXPROF
        CALL TXEXEC
        TMAX=TIME+DT*NTMAX
     ELSE IF (KID == 'C') THEN
        IF (ICONT == 0) THEN
           WRITE(6,*) 'XX RUN or LOAD before CONTINUE !'
           CYCLE
        ENDIF
        NGR=-1
        CALL TXSTGR
        CALL TXEXEC
        TMAX=TIME+DT*NTMAX
     ELSE IF (KID == 'P') THEN
        CALL TXPARM(KID)
        IF(KID == 'Q') EXIT
     ELSE IF (KID == 'V') THEN
        CALL TXVIEW
     ELSE IF (KID == 'I') THEN
        CALL TXINIT
     ELSE IF (KID == 'Q') THEN
        EXIT
     ELSE IF (KID == 'S') THEN
        CALL TXSAVE
     ELSE IF (KID == 'L') THEN
        CALL TXLOAD
        IERR = 0
        NGR = -1
        ICONT = 1
        CALL TXSTGR
     ELSE IF (KID == 'G') THEN
        CALL TXGOUT
     ELSE IF (KID == 'N') THEN
        I = NINT((DelR / DR) - 0.5D0)
        X(I,1) = X(I,1) + DelN
        X(I,2) = X(I,2) + DelN
        CALL TXCALV(X)
     ELSE IF (KID == 'B') THEN
        KID2=LINE(2:2)
        CALL TOUPPER(KID2)
        IF      (KID2 == 'P') THEN
           PNBCD =  1.D0
        ELSE IF (KID2 == '0') THEN
           PNBCD =  0.D0
        ELSE IF (KID2 == 'M') THEN
           PNBCD = -1.D0
        ELSE
           WRITE(6,*) 'XX Unknown beam command'
        ENDIF
     ELSE IF (KID == '#') THEN
        CONTINUE
     ELSE
        WRITE(6,*) 'XX Unknown command'
     ENDIF
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
  CHARACTER LINE*80,KID*1

!  ----- read line input input -----

  READ(5,'(A80)',IOSTAT=IST) LINE

  IF(IST == 0) THEN
     !  ----- parameter input -----
     
     ID=0
     DO I=1,80
        IF(LINE(I:I) == '=') ID=1
     ENDDO
     IF(ID == 1) THEN
        CALL TXPARL(LINE)
        KID=' '
        MODE=2
        RETURN
     ENDIF

     !  ----- command input -----

     KID=LINE(1:1)
     CALL TOUPPER(KID)
     IF(KID >= 'A'.AND.KID <= 'Z') THEN
        MODE=1
        RETURN
     ENDIF

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

  ENDIF

  RETURN
END SUBROUTINE TXKLIN
