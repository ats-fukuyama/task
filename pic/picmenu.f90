!  ***** TASK/PIC MENU *****

MODULE picmenu

PRIVATE
PUBLIC pic_menu

CONTAINS

  SUBROUTINE pic_menu

    USE piccomm,ONLY: ikind,nrank,nsize,pic_deallocate
    USE picparm,ONLY: pic_parm,pic_view,pic_broadcast
    USE picprep,ONLY: pic_prep
    USE picexec,ONLY: pic_exec
    USE picgout,ONLY: pic_gout
    USE libmpi

    IMPLICIT NONE
    INTEGER(ikind)    :: ierr,mode,ind,init=0
    CHARACTER         :: kid
    CHARACTER(LEN=80) :: line

1   CONTINUE
    IF(nrank == 0) THEN
       ierr=0
       WRITE(6,'(A)') &
            '## PIC MENU: P,V/Parm  R/Run  C/Continue  G/Graf  Q/QUIT'
       CALL TASK_KLIN(line,kid,mode,pic_parm)
    END IF
    CALL mtx_barrier
    CALL mtx_broadcast1_character(kid)
    CALL mtx_broadcast1_integer(mode)
    IF(mode /= 1) GOTO 1

    CALL pic_broadcast

    IF(kid.EQ.'P') THEN
       IF(nrank == 0) CALL pic_parm(0,'PIC',ierr)
       CALL pic_broadcast
    ELSEIF(kid.EQ.'V') THEN
       IF(nrank == 0) CALL pic_view
    ELSEIF(kid.EQ.'R') THEN
       CALL pic_prep(ierr)
       if(ierr /= 0) GO TO 1
       init=1
       CALL pic_exec(ierr)
    ELSEIF(kid.EQ.'C') THEN
       IF(init /= 0) CALL pic_exec(ierr)
    ELSEIF(kid.EQ.'G') THEN
       IF(nrank == 0 .AND. init /= 0) CALL pic_gout
    ELSEIF(kid.EQ.'Q') THEN
       GOTO 9000
    ELSE
       WRITE(6,*) 'XX PICMENU: UNKNOWN kid: kid = ',kid
    ENDIF
    CALL mtx_barrier
    GOTO 1

9000 CONTINUE
    CALL pic_deallocate
    RETURN
  END SUBROUTINE pic_menu
END MODULE picmenu
