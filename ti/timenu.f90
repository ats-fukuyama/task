!   MODULE timenu

MODULE timenu

  PUBLIC ti_menu

CONTAINS

  SUBROUTINE ti_menu

      USE bpsd
      USE tiinit
      USE tiparm
      USE tiprep
      USE tiexec
      USE tigout
      USE plload,ONLY: pl_load
      
      IMPLICIT NONE
      INTEGER(ikind)       :: IERR, MODE
      INTEGER(ikind), SAVE :: INIT=0
      CHARACTER(LEN=1) :: KID
      CHARACTER(LEN=80):: LINE

!     ------ SELECTION OF TASK TYPE ------

      IERR=0

    1 IF(INIT.EQ.0) THEN
         WRITE(6,601)
  601    FORMAT('## TI MENU: P,V/PARM  R/RUN  L/LOAD  ', &
                            'W/WRITE  H/HELP  Q/QUIT')
      ELSEIF(INIT.EQ.1) THEN
         WRITE(6,602)
  602    FORMAT('## TI MENU: G/GRAPH  '/ &
                            'P,V/PARM  R/RUN  L/LOAD  ', &
                            'W/WRITE  H/HELP  Q/QUIT')
      ELSE
         WRITE(6,603)
  603    FORMAT('## TI MENU: C/CONT  G/GRAPH  ', &
                            'P,V/PARM  R/RUN  L/LOAD  ',&
                            'W/WRITE  H/HELP  Q/QUIT')
      ENDIF

      CALL TASK_KLIN(LINE,KID,MODE,ti_parm)
      IF(MODE.NE.1) GOTO 1

      IF(KID.EQ.'P') THEN
         CALL ti_parm(0,'TI',IERR)
      ELSE IF(KID.EQ.'V') THEN
         CALL ti_view
      ELSE IF(KID.EQ.'L') THEN
         CALL pl_load(ierr)
         if(ierr.ne.0) GO TO 1
      ELSE IF(KID.EQ.'R') THEN
         CALL ti_prep(ierr)
         if(ierr.ne.0) GO TO 1
         CALL ti_exec(ierr)
         INIT=2
      ELSE IF(KID.EQ.'C'.AND.INIT.EQ.2) THEN
         CALL ti_exec(ierr)
      ELSE IF(KID.EQ.'G'.AND.INIT.GE.1) THEN
         CALL ti_gout
      ELSE IF(KID.EQ.'W'.AND.INIT.EQ.2) THEN
  102    WRITE(6,*) '# SELECT ',  ': PRINT TYPE (1..9)  N/NAMELIST X/EXIT'
         READ(5,'(A1)',ERR=102,END=1) KID
         CALL GUCPTL(KID)
         IF(KID.EQ.'H') THEN
!            CALL ti_help('W')
         ELSEIF(KID.EQ.'X') THEN
            GOTO 1
         ELSE
!            CALL ti_print(KID)
         ENDIF
         GOTO 102
      ELSE IF(KID.EQ.'H') THEN
!         CALL ti_help('M')
      ELSE IF(KID.EQ.'Q') THEN
         GOTO 9000

      ELSE IF(KID.EQ.'X'.OR.KID.EQ.'#') THEN
         CONTINUE
      ELSE
         WRITE(6,*) 'XX timenu: UNKNOWN KID'
      END IF

      GOTO 1

 9000 CONTINUE
      RETURN
    END SUBROUTINE ti_menu
  END MODULE timenu
