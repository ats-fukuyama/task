! obmenu.f90

MODULE obmenu

  PRIVATE
  PUBLIC ob_menu

CONTAINS

!     ***** TASK/OB MENU *****

  SUBROUTINE OB_MENU

    USE obcomm
    USE plparm,ONLY: pl_view
    USE obparm,ONLY: ob_parm
    USE obview,ONLY: ob_view
    USE obprep,ONLY: ob_prep
    USE obcalc,ONLY: ob_calc
    USE obgout,ONLY: ob_gout
!    USE obfile,ONLY: ob_save,ob_load
    IMPLICIT NONE
    CHARACTER(LEN=1):: KID
    CHARACTER(LEN=80):: LINE
    INTEGER:: NSTAT,IERR,MODE,NID

      NSTAT=0

    1 CONTINUE
         IERR=0
         WRITE(6,601)
  601    FORMAT('## OB MENU: P,V/PARM  R/RUN  G/GRAPH  S,L/FILE', &
                '  Q/QUIT')
         CALL TASK_KLIN(LINE,KID,MODE,OB_PARM)
      IF(MODE.NE.1) GOTO 1

      IF(KID.EQ.'P') THEN
         CALL OB_PARM(0,'OB',IERR)
      ELSEIF(KID.EQ.'V') THEN
         CALL PL_VIEW
         CALL OB_VIEW
      ELSEIF(KID.EQ.'R') THEN
         CALL ob_prep(ierr)
         CALL ob_allocate
         CALL OB_CALC(ierr)
         NSTAT=1
      ELSEIF(KID.EQ.'G') THEN
         CALL OB_GOUT
      ELSEIF(KID.EQ.'S') THEN
!         CALL OB_SAVE
      ELSEIF(KID.EQ.'L') THEN
!         CALL OB_LOAD
      ELSEIF(KID.EQ.'Q') THEN
         GOTO 9000
      ELSE
         WRITE(6,*) 'XX OBMENU: UNKNOWN KID'
      ENDIF
      GOTO 1

9000  CALL ob_deallocate
      RETURN
  END SUBROUTINE OB_MENU
END MODULE obmenu
