MODULE trmenu

  PUBLIC tr_menu

  PRIVATE

CONTAINS

  SUBROUTINE tr_menu

    USE trcomm, ONLY : &
         ikind,rkind
    USE trinit
    USE trsetup
    USE trloop
    USE trgout

    IMPLICIT NONE
    INTEGER(ikind)      :: IERR, MODE ! , NFL, NFLMAX, NTMOLD
    INTEGER(ikind),SAVE :: INIT=0
    CHARACTER(LEN=1) :: KID
    CHARACTER(LEN=80):: LINE

!     ------ SELECTION OF TASK TYPE ------

    IERR=0
1   IF(INIT.EQ.0) THEN
       WRITE(6,601)
601    FORMAT('## TR MENU: P,V/PARM  R/RUN  L/LOAD  ', &
              'D/DATA  H/HELP  Q/QUIT')
    ELSEIF(INIT.EQ.1) THEN
       WRITE(6,602)
602    FORMAT('## TR MENU: G/GRAPH  '/ &
              '            P,V/PARM  R/RUN  L/LOAD  ', &
              'D/DATA  H/HELP  Q/QUIT')
    ELSE
       WRITE(6,603)
603    FORMAT('## TR MENU: C/CONT  G/GRAPH  W,F/WRITE  ', &
              'S/SAVE  O/UFILEOUT  M/MDLTST'/ &
              '            P,V/PARM  R/RUN  L/LOAD  ',&
              'D/DATA  H/HELP  Q/QUIT')
    ENDIF

    CALL TASK_KLIN(LINE,KID,MODE,tr_parm)
    IF(MODE.NE.1) GOTO 1

    IF(KID.EQ.'P') THEN
       CALL tr_parm(0,'TR',ierr)
    ELSE IF(KID.EQ.'V') THEN
       CALL tr_view

    ELSE IF(KID.EQ.'L') THEN
!       CALL tr_load(ierr)
       INIT=2
    ELSE IF(KID.EQ.'S'.AND.INIT.EQ.2) THEN
!       CALL tr_save(ierr)

    ELSE IF(KID.EQ.'R') THEN
       CALL tr_idneq
       CALL tr_setup
       if(ierr.ne.0) GO TO 1
       CALL tr_loop
       INIT=2

    ELSE IF(KID.EQ.'C'.AND.INIT.EQ.2) THEN
       CALL tr_loop

    ELSE IF(KID.EQ.'G'.AND.INIT.GE.1) THEN
       CALL tr_gout

    ELSE IF(KID.EQ.'F'.AND.INIT.GE.1) THEN
!       CALL tr_fout

    ELSE IF(KID.EQ.'W'.AND.INIT.EQ.2) THEN
102    WRITE(6,*) '# SELECT ',  ': PRINT TYPE (1..9)  X/EXIT'
       READ(5,'(A1)',ERR=102,END=1) KID
       CALL GUCPTL(KID)
       IF(KID.EQ.'X') GO TO 1
       !         CALL TRPRNT(KID)
       GOTO 102

    ELSE IF(KID.EQ.'O'.AND.INIT.EQ.2) THEN
       !         CALL TRXOUT
    ELSE IF(KID.EQ.'M'.AND.INIT.EQ.2) THEN
       !         CALL TRMDLT
    ELSE IF(KID.EQ.'H') THEN
       !         CALL TRHELP('M')
    ELSE IF(KID.EQ.'Q') THEN
       GOTO 9000

    ELSE IF(KID.EQ.'X'.OR.KID.EQ.'#') THEN
       CONTINUE
    ELSE
       WRITE(6,*) 'XX TRMAIN: UNKNOWN KID'
    END IF

    GOTO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE tr_menu
END MODULE trmenu
