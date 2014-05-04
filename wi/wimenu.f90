!  $Id$

!  ***** TASK/WI MENU *****

MODULE wimenu

CONTAINS

  SUBROUTINE wi_menu

    USE wicomm,ONLY: ikind,rkind,wi_deallocate,alfa,beta,xmax,pn0,any, &
                     nxmax,nwmax,pnu,dx0,xwint,dxmin,xwmin,modelp
    USE wiparm,ONLY: wi_parm,wi_view
    USE wiprep,ONLY: wi_prep
    USE wiexec,ONLY: wi_exec
    USE wiscan,ONLY: wi_scan
    USE wigout,ONLY: wi_gout

    IMPLICIT NONE
    INTEGER(ikind)    :: ierr,mode,ind
    CHARACTER         :: kid
    CHARACTER(LEN=80) :: line
    INTEGER(ikind)    :: init=0
    REAL(rkind)       :: dx0_save=0.D0,dxmin_save=0.D0
    REAL(rkind)       :: xwmin_save=0.D0,xmax_save=0.D0
    REAL(rkind)       :: xwint_save=0.D0
    REAL(rkind)       :: ratea,rk0l

1   CONTINUE
    ierr=0
    WRITE(6,'(A)') '## WI MENU: P,V/PARM  R/RUN  S/SCAN  G/GRAF  Q/QUIT'

    CALL TASK_KLIN(line,kid,mode,wi_parm)
    IF(mode /= 1) GOTO 1
    IF(dx0  .NE.dx0_save   .OR. &
       xwint.NE.xwint_save .OR. &
       dxmin.NE.dxmin_save .OR. &
       xwmin.NE.xwmin_save .OR. &
       xmax .NE.xmax_save) THEN  ! x array was modified
       INIT=0
    END IF

    IF(kid.EQ.'P') THEN
       CALL wi_parm(0,'WI',ierr)
    ELSEIF(kid.EQ.'V') THEN
       CALL wi_view
    ELSEIF(kid.EQ.'R') THEN
       CALL wi_prep
       write(6,'(A,1P4E12.4)')     '## alfa,beta,pn0,any       =', &
                                       alfa,beta,pn0,any
       write(6,'(A,1PE12.4,3I12)') '## xmax,nxmax,nwmax,modelp =', &
                                       xmax,nxmax,nwmax,modelp
       SELECT CASE(modelp)
       CASE(0,1)
          write(6,'(A,1P4E12.4)')     '## dx0,pnu,dxmin,xwmin     =', &
                                          dx0,pnu,dxmin,xwmin
       CASE(2)
          write(6,'(A,1P4E12.4)')     '## dx0,xwint,dxmin,xwmin   =', &
                                          dx0,xwint,dxmin,xwmin
       END SELECT
       dx0_save  =dx0
       xwint_save=xwint
       dxmin_save=dxmin
       xwmin_save=xwmin
       xmax_save =xmax
       CALL wi_exec(1,ratea,ierr)
       INIT=1
    ELSEIF(kid.EQ.'S') THEN
       CALL wi_prep
       rk0l=beta/alfa
       write(6,'(A,1P4E12.4)')     '## alfa,beta,pn0,rk0l      =', &
                                       alfa,beta,pn0,any
       write(6,'(A,1PE12.4,3I12)') '## xmax,nxmax,nwmax,modelp =', &
                                       xmax,nxmax,nwmax,modelp
       SELECT CASE(modelp)
       CASE(0,1)
          write(6,'(A,1P4E12.4)')     '## dx0,pnu,dxmin,xwmin     =', &
                                          dx0,pnu,dxmin,xwmin
       CASE(2)
          write(6,'(A,1P4E12.4)')     '## dx0,xwint,dxmin,xwmin   =', &
                                          dx0,xwint,dxmin,xwmin
       END SELECT
       dx0_save  =dx0
       xwint_save=xwint
       dxmin_save=dxmin
       xwmin_save=xwmin
       xmax_save =xmax
       CALL wi_scan(ierr)
       INIT=1
    ELSEIF(kid.EQ.'G') THEN
       IF(INIT.EQ.0) THEN
          WRITE(6,*) 'XX data is not ready or destroyed'
       ELSE
          CALL wi_gout
       END IF
    ELSEIF(kid.EQ.'Q') THEN
       GOTO 9000
    ELSE
       WRITE(6,*) 'XX WIMENU: UNKNOWN kid: kid = ',kid
    ENDIF
    GOTO 1

9000 CALL wi_deallocate
    RETURN
  END SUBROUTINE wi_menu
END MODULE wimenu
