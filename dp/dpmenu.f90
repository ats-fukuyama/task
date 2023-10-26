! dpmenu.f90

MODULE dpmenu

PRIVATE

PUBLIC dp_menu

CONTAINS

!     ***** TASK/DP MENU *****

  SUBROUTINE dp_menu

    USE dpcomm_local
    USE plprof
    USE plprofw
    USE plparm,ONLY: pl_view
    USE dpparm,ONLY: dp_parm,dp_view
    USE dproot,ONLY: dp_root,dpgrp1,dpgrp0
    USE dpcont,ONLY: dp_cont2,dp_cont3
    USE dpcont4,ONLY: dp_cont4
    USE dptest,ONLY: dp_test
    USE libkio
    USE libchar
    IMPLICIT NONE
      
    CHARACTER(LEN=1):: KID
    CHARACTER(LEN=80):: LINE
    INtEGER:: MODE,IERR,NID

1   CONTINUE
    WRITE(6,*) '## DP MENU: P,V/PARM  ', &
               'D0,D1,D2,D3,D4,D5,D6/DISP F/ROOT T,B/TEST ?/HELP Q/QUIT'

    CALL TASK_KLIN(LINE,KID,MODE,DP_PARM)
    IF(MODE.NE.1) GOTO 1
    CALL toupper(KID)

    IF(KID.EQ.'P') THEN
       CALL DP_PARM(0,'DP',IERR)
    ELSEIF(KID.EQ.'V') THEN
       CALL PL_VIEW
       CALL DP_VIEW
    ELSEIF(KID.EQ.'D') THEN
       READ(LINE(2:),*,ERR=1,END=1) NID
       SELECT CASE(NID)
       CASE(0)
          CALL DPGRP0
       CASE(1)
          CALL DPGRP1
       CASE(2)
          CALL DP_CONT2
       CASE(3)
          CALL DP_CONT3
       CASE(4,5,6)
          CALL DP_CONT4(NID)
       CASE DEFAULT
          WRITE(6,*) 'XX DPMENU: unknown NID'
       END SELECT
    ELSEIF(KID.EQ.'F') THEN
       CALL DP_ROOT
    ELSEIF(KID.EQ.'T') THEN
       CALL dp_test
       GOTO 1
    ELSEIF(KID.EQ.'?') THEN
       CALL dp_help(0)
    ELSEIF(KID.EQ.'Q') THEN
       GOTO 9000
    ENDIF
    GOTO 1

9000 RETURN
  END SUBROUTINE dp_menu

  ! --- Help message ---

  SUBROUTINE dp_help(id)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: id

    SELECT CASE(id)
    CASE(0)
       WRITE(6,'(A)') 'Menu Help:  D=the determinant of dispersion&
            & matrix'
       WRITE(6,'(A)') 'D0: Show D for various plasma model, '
       WRITE(6,'(A)') 'D1: Draw D as a function of one parameter, omega, k'
       WRITE(6,'(A)') 'D2: Draw contour of D=0 on 2D parameter plane'
       WRITE(6,'(A)') 'D3: Draw contours of D=0 on 2D for various density'
       WRITE(6,'(A)') 'D4: Draw disp curv D=0 and plot with imag on (omega,*)'
       WRITE(6,'(A)') 'D5: Plot D=0 with imag color'
       WRITE(6,'(A)') 'D6: Plot D=0 with shifted imag with color'
       WRITE(6,'(A)') 'F : Find a root of D=0'
       WRITE(6,'(A)') 'T : show D for various plasma model for EC waves'
       WRITE(6,'(A)') 'B : show D for various plasma model for EC waves'
       WRITE(6,'(A)') ' '
       WRITE(6,'(A)') 'P : parameter input through namelist'
       WRITE(6,'(A)') 'V : show input parameters'
       WRITE(6,'(A)') '? : this helpmessage'
       WRITE(6,'(A)') 'Q : end of this run.'
       RETURN
    END SELECT
  END SUBROUTINE dp_help
       
END MODULE dpmenu
