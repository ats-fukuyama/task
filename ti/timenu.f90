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
      USE libmpi
      USE commpi
      
      IMPLICIT NONE
      INTEGER(ikind)       :: IERR, MODE
      INTEGER(ikind), SAVE :: INIT=0
      CHARACTER(LEN=1) :: KID
      CHARACTER(LEN=80):: LINE

!     ------ SELECTION OF TASK TYPE ------

      IERR=0

    1 CONTINUE
      IF(nrank.EQ.0) THEN
         IF(INIT.EQ.0) THEN
            WRITE(6,601)
601         FORMAT('## TI MENU: P,V/PARM  R/RUN  L/LOAD  ', &
                               'W/WRITE  H/HELP  Q/QUIT')
         ELSEIF(INIT.EQ.1) THEN
            WRITE(6,602)
602         FORMAT('## TI MENU: G/GRAPH  '/ &
                               'P,V/PARM  R/RUN  L/LOAD  ', &
                               'W/WRITE  S/SAVE  Q/QUIT')
         ELSE
            WRITE(6,603)
603         FORMAT('## TI MENU: C/CONT  G/GRAPH  ', &
                               'P,V/PARM  R/RUN  L/LOAD  ',&
                               'W/WRITE  S/SAVE  Q/QUIT')
         ENDIF

         CALL TASK_KLIN(LINE,KID,MODE,ti_parm)
      END IF
      CALL mtx_broadcast1_character(KID)
      CALL mtx_broadcast1_integer(MODE)

      IF(MODE.EQ.2) CALL ti_broadcast
      IF(MODE.NE.1) GOTO 1

      IF(KID.EQ.'P') THEN           ! parameter input
         IF(nrank.eq.0) CALL ti_parm(0,'TI',IERR)
         CALL ti_broadcast
      ELSE IF(KID.EQ.'V') THEN      ! view input parameters
         IF(nrank.eq.0) CALL ti_view
      ELSE IF(KID.EQ.'L') THEN      ! load bulk component profile from TR
         CALL pl_load(ierr)
         if(ierr.ne.0) GO TO 1
      ELSE IF(KID.EQ.'R') THEN      ! run simulation
         CALL ti_prep(ierr)
         if(ierr.ne.0) GO TO 1
         CALL ti_exec(ierr)
         INIT=2
      ELSE IF(KID.EQ.'C'.AND.INIT.EQ.2) THEN  ! continue simulation
         CALL ti_exec(ierr)
      ELSE IF(KID.EQ.'G'.AND.INIT.GE.1) THEN  ! graphic output
         IF(nrank.EQ.0) CALL ti_gout
      ELSE IF(KID.EQ.'W'.AND.INIT.EQ.2) THEN  ! show simulation results in text
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
      ELSE IF(KID.EQ.'H') THEN  ! show helpful information
!         CALL ti_help('M')
      ELSE IF(KID.EQ.'S') THEN  ! save results
!         CALL ti_save(ierr)
      ELSE IF(KID.EQ.'Q') THEN  ! quit 
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
