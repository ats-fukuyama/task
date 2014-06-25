!
!
!
MODULE T2MENU

PRIVATE
PUBLIC T2_MENU

CONTAINS

  SUBROUTINE T2_MENU

    USE T2CNST,ONLY: ikind,rkind
    USE T2COMM,ONLY: T2COMM_DEALLOCATE,T2NGRA_DEALLOCATE,&!
         !debug
         d2xout,Xvec
    USE T2PARM,ONLY: T2_PARM,T2_VIEW
    USE T2PREP,ONLY: T2PREP_EXECUTE
    USE T2DIV, ONLY: T2_DIV
    USE T2LOOP,ONLY: T2_LOOP
    USE T2GOUT,ONLY: T2_GOUT

    IMPLICIT NONE

    INTEGER(ikind)    :: ierr,mode,ind
    CHARACTER         :: kid
    CHARACTER(LEN=80) :: line
    INTEGER(ikind)    :: init=0

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') &
         '## T2 MENU: P,V/PARM  R/RUN  C/CONT  W/FILE  G/GOUT  Q/QUIT'

    CALL TASK_KLIN(line,kid,mode,T2_PARM)
    IF(mode /= 1) GOTO 1

    IF(kid.EQ.'P') THEN
       CALL T2_PARM(0,'T2',ierr)
    ELSEIF(kid.EQ.'V') THEN
       CALL T2_VIEW
    ELSEIF(kid.EQ.'R') THEN
       CALL T2PREP_EXECUTE
       CALL T2_LOOP
       INIT=1
    ELSEIF(kid.EQ.'C') THEN
       CALL T2_LOOP
       INIT=1
    ELSEIF(kid.EQ.'W') THEN
       IF(INIT.EQ.0) THEN
          WRITE(6,*) 'XX data is not ready or destroyed'
       ELSE
!          CALL T2_WRIT
       END IF
    ELSEIF(kid.EQ.'G') THEN
       IF(INIT.EQ.0) THEN
          WRITE(6,*) 'XX data is not ready or destroyed'
       ELSE
          !print*,'AAAAA'
          !d2xout = Xvec!for_debug
          CALL T2_GOUT
       END IF
    ELSEIF(kid.EQ.'Q') THEN
       GOTO 9000
    ELSE
       WRITE(6,*) 'XX T2_MENU: UNKNOWN kid: kid = ',kid
    ENDIF
    GOTO 1

9000 CALL T2COMM_DEALLOCATE
    CALL T2NGRA_DEALLOCATE
    RETURN
  END SUBROUTINE T2_MENU
END MODULE T2MENU
