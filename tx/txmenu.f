C
C     ***************************************************************
C
C        Command loop
C
C     ***************************************************************
C
      SUBROUTINE TXMENU
C
      INCLUDE 'txcomm.inc'
C
      CHARACTER LINE*80,KID*1,KID2*1
C
      IERR=0
      ICONT = 0
C
C     *** MENU ***
C
    1 CONTINUE
C
      WRITE(6,'(3(A,1PD12.4))')
     &   ' ## TIME=',TIME,'  DT=',DT,'  TMAX =',TMAX
      WRITE(6,*) '## INPUT: ',
     &   'R:RUN  C:CONT  P,V:PARM  G:GRAPH '//
     &   'S:SAVE  L:LOAD  I,N,Bn: Q:QUIT'
      CALL GUFLSH
C
      CALL TXKLIN(LINE,KID,MODE)
      IF(MODE.NE.1) GOTO 1
C
  100 IF      (KID .EQ. 'R') THEN
         TIME = 0.D0
         TPRE = 0.D0
         IERR = 0
         ICONT = 1
         CALL TXPROF
         CALL TXEXEC
         TMAX=TIME+DT*NTMAX
      ELSE IF (KID .EQ. 'C') THEN
         IF (ICONT .EQ. 0) THEN
            WRITE(6,*) 'XX RUN or LOAD before CONTINUE !'
            GOTO 1
         ENDIF
         NGR=-1
         CALL TXSTGR
         CALL TXEXEC
         TMAX=TIME+DT*NTMAX
      ELSE IF (KID .EQ. 'P') THEN
         CALL TXPARM(KID)
         IF(KID.EQ.'Q') GOTO 9000
      ELSE IF (KID .EQ. 'V') THEN
         CALL TXVIEW
      ELSE IF (KID .EQ. 'I') THEN
         CALL TXINIT
      ELSE IF (KID .EQ. 'Q') THEN
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
         CALL TXGOUT
      ELSE IF (KID .EQ. 'N') THEN
         I = NINT((DelR / DR) - 0.5D0)
         X(I,1) = X(I,1) + DelN
         X(I,2) = X(I,2) + DelN
         CALL TXCALV(X)
      ELSE IF (KID .EQ. 'B') THEN
         KID2=LINE(2:2)
         CALL GUCPTL(KID2)
         IF      (KID2 .EQ. 'P') THEN
            PNBCD =  1.D0
         ELSE IF (KID2 .EQ. '0') THEN
            PNBCD =  0.D0
         ELSE IF (KID2 .EQ. 'M') THEN
            PNBCD = -1.D0
         ELSE
            WRITE(6,*) 'XX Unknown beam command'
         ENDIF
      ELSE IF (KID .EQ. '#') THEN
         CONTINUE
      ELSE
         WRITE(6,*) 'XX Unknown command'
      ENDIF
      GOTO 1
C
 9000 RETURN
      END
C
C     ***** INPUT KID or LINE *****
C                   MODE=0: LINE INPUT 
C                        1: KID INPUT
C                        2: PARM INPUT
C                        3: NEW PROMPT
C
      SUBROUTINE TXKLIN(LINE,KID,MODE)
C
      CHARACTER LINE*80,KID*1
C
C     ----- read line input input -----
C
      READ(5,'(A80)',ERR=2,END=3) LINE
C
C     ----- parameter input -----
C
      ID=0
      DO I=1,80
         IF(LINE(I:I).EQ.'=') ID=1
      ENDDO
      IF(ID.EQ.1) THEN
         CALL TXPARL(LINE)
         KID=' '
         MODE=2
         RETURN
      ENDIF
C
C     ----- command input -----
C
      KID=LINE(1:1)
      CALL GUCPTL(KID)
      IF(KID.GE.'A'.AND.KID.LE.'Z') THEN
         MODE=1
         RETURN
      ENDIF
C
C     ----- line input -----
C
      KID=' '
      MODE=0
      RETURN
C
C     ----- input error -----
C
    2 WRITE(6,*) 'XX INPUT ERROR !'
      MODE=3
      RETURN
C
C     ----- input end -----
C
    3 KID='Q'
      MODE=1
      RETURN
      END
