!  ***** TASK/XX MENU *****

MODULE picmenu

PRIVATE
PUBLIC pic_menu

CONTAINS

  SUBROUTINE pic_menu

    USE piccomm,ONLY: ikind,rkind,pic_allocate,pic_deallocate
    USE picparm,ONLY: pic_parm,pic_view
    USE picprep,ONLY: pic_prep
    USE picexec,ONLY: pic_exec
    USE picgout,ONLY: pic_gout

    IMPLICIT NONE
    INTEGER(ikind)    :: ierr,mode,ind
    CHARACTER         :: kid
    CHARACTER(LEN=80) :: line

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') '## PIC MENU: P,V/Parm  R/Run  C/Continue G/Graf  Q/QUIT'

    CALL TASK_KLIN(line,kid,mode,pic_parm)
    IF(mode /= 1) GOTO 1

    IF(kid.EQ.'P') THEN
       CALL pic_parm(0,'PIC',ierr)
    ELSEIF(kid.EQ.'V') THEN
       CALL pic_view
    ELSEIF(kid.EQ.'R') THEN
       CALL pic_prep(ierr)
       if(ierr == 1) GO TO 1
       CALL pic_exec(ierr)
    ELSEIF(kid.EQ.'C') THEN
       CALL pic_exec(ierr)
    ELSEIF(kid.EQ.'G') THEN
       CALL pic_gout
    ELSEIF(kid.EQ.'Q') THEN
       GOTO 9000
    ELSE
       WRITE(6,*) 'XX PICMENU: UNKNOWN kid: kid = ',kid
    ENDIF
    GOTO 1

9000 CONTINUE
    CALL pic_deallocate
    RETURN
  END SUBROUTINE pic_menu
END MODULE picmenu
