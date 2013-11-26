MODULE T2MENU

  PUBLIC T2_MENU

  PRIVATE

CONTAINS
  
  SUBROUTINE T2_MENU

    USE T2CNST, ONLY: i0rkind, i0ikind
    
    IMPLICIT NONE
    INTEGER(i0ikind)    :: i0err,i0ioerr, i0mode
    INTEGER(i0kind),SAVE:: i0init=0
    CHARACTER(LEN= 1)   :: a0kid
    CHARACTER(LEN=80)   :: a0line

    !C SELECTION OF TASK TYPE

    i0err = 0

    DO 
       IF(i0init.EQ.0)THEN
          WRITE(6,601)
601       FORMAT('## T2 MENU: P,V/PARM  R/RUN  L/LOAD  ', &
                 'D/DATA  H/HELP  Q/QUIT')
       ELSEIF(INIT.EQ.1) THEN
          WRITE(6,602)
602       FORMAT('## T2 MENU: G/GRAPH  '/ &
                 '            P,V/PARM  R/RUN  L/LOAD  ', &
                 'D/DATA  H/HELP  Q/QUIT')
       ELSE
          WRITE(6,603)
603       FORMAT('## T2 MENU: C/CONT  G/GRAPH  W,F/WRITE  ', &
                 'S/SAVE  O/UFILEOUT  M/MDLTST'/ &
                 '            P,V/PARM  R/RUN  L/LOAD  ',&
                 'D/DATA  H/HELP  Q/QUIT')
       ENDIF
       
       CALL TASK_KLIN(a0line,a0kid,i0mode,t2_parm)
       IF(i0mode.NE.1)  CYCLE
       
       IF(     a0kid.EQ.'P') THEN
          CALL T2_PARM(0,'T2',ierr)
          
       ELSE IF(a0kid.EQ.'V') THEN
          CALL T2_VIEW

       ELSE IF(a0kid.EQ.'L') THEN
          !CALL t2_load(i0err)
          i0init=2

       ELSE IF((a0kid.EQ.'S').AND.(i0init.EQ.2)) THEN
!          CALL t2_save(0ierr)

       ELSE IF(a0kid.EQ.'R') THEN
          CALL T2_SETUP
          IF(i0err.NE.0) CYCLE
          CALL T2_LOOP
          i0init=2
          
       ELSE IF((a0kid.EQ.'C').AND.(i0init.EQ.2)) THEN
          CALL T2_LOOP

       ELSE IF((a0kid.EQ.'G').AND.(i0init.GE.1)) THEN
          !CALL tr_gout

       ELSE IF((a0kid.EQ.'F').AND.(i0init.GE.1)) THEN
          !CALL tr_fout

       ELSE IF((a0kid.EQ.'W').AND.(i0init.EQ.2)) THEN
          DO
             WRITE(6,*) '# SELECT ',  ': PRINT TYPE (1..9)  X/EXIT'
             READ(5,'(A1)',IOSTAT=i0ioerr) a0kid
             IF(i0ioerr.LT.0) EXIT
             IF(i0ioerr.GT.0) CYCLE
             !CALL GUCPTL(a0kid)

             IF(a0kid.EQ.'X') EXIT
!             CALL TRPRNT(KID)
          END DO
          CYCLE

       ELSE IF((a0kid.EQ.'O').AND.(i0init.EQ.2)) THEN
!          CALL TRXOUT
       ELSE IF((a0kid.EQ.'M').AND.(i0init.EQ.2)) THEN
!          CALL TRMDLT
       ELSE IF(a0kid.EQ.'H') THEN
!          CALL TRHELP('M')
       ELSE IF(a0kid.EQ.'Q') THEN
          EXIT

       ELSE IF(KID.EQ.'X'.OR.KID.EQ.'#') THEN
          CONTINUE
       ELSE
          WRITE(6,*) 'XX T2MAIN: UNKNOWN KID'
       END IF

    END DO

    RETURN
  END SUBROUTINE T2_MENU
