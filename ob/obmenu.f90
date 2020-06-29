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
    CHARACTER(LEN=1):: kid
    CHARACTER(LEN=80):: line
    INTEGER:: nstat,ierr,mode,nid

      nstat=0

    1 CONTINUE
         ierr=0
         WRITE(6,601)
  601    FORMAT('## OB MENU: P,V/PARM  R/RUN  G/GRAPH  S,L/FILE', &
                '  Q/QUIT')
         CALL task_klin(line,kid,mode,ob_parm)
      IF(mode.NE.1) GOTO 1

      IF(kid.EQ.'P') THEN
         CALL ob_parm(0,'OB',ierr)
      ELSEIF(kid.EQ.'V') THEN
         CALL pl_view
         CALL ob_view
      ELSEIF(kid.EQ.'R') THEN
         CALL ob_prep(ierr)
         CALL ob_allocate
         CALL ob_calc(ierr)
         nstat=1
      ELSEIF(kid.EQ.'G') THEN
         IF(nstat.EQ.0) CALL ob_prep(ierr)
         CALL ob_gout
      ELSEIF(kid.EQ.'S') THEN
!         CALL ob_save
      ELSEIF(kid.EQ.'L') THEN
!         CALL ob_load
      ELSEIF(kid.EQ.'Q') THEN
         GOTO 9000
      ELSE
         WRITE(6,*) 'XX ob_menu: unknown kid'
      ENDIF
      GOTO 1

9000  IF(ALLOCATED(obt_in)) CALL ob_deallocate
      RETURN
  END SUBROUTINE OB_MENU
END MODULE obmenu
