! wmmenu.f90

!     ***** TASK/WM MENU *****

MODULE wmmenu

  PRIVATE
  PUBLIC wm_menu

CONTAINS

  SUBROUTINE wm_menu

    USE wmcomm
    USE wminit,ONLY: wm_init
    USE wmparm,ONLY: wm_parm,wm_broadcast
    USE wmview,ONLY: wm_view
    USE wmloop,ONLY: wm_loop
    USE wmdout,ONLY: wm_dout
    USE wmgout,ONLY: wm_gout
    USE wmfile,ONLY: wm_load,wm_save
    USE wmeign,ONLY: wm_am0d,wm_am1d,wm_am2d,wm_scan,wm_eign,wm_wout
    USE wmtaem,ONLY: wm_taem
    USE libkio
    USE libmpi

    IMPLICIT NONE
    INTEGER:: init=0
    CHARACTER(LEN=1):: kid
    CHARACTER(LEN=80):: line
    INTEGER:: ierr,mode,nid
    
! --- selecion of action type ----

    ierr=0

1   CONTINUE
    IF(nrank.EQ.0) THEN
       IF(init.EQ.0) THEN
          WRITE(6,601)
601       FORMAT('## WM MENU: P,V/parm R/run D0-3/amp F/root ', &
                 'T/tae L/load Q/quit')
       ELSE
          WRITE(6,602)
602       FORMAT('## WM MENU: P,V/parm R/run D0-3/amp F/root ', &
                 'G,Y/graph T/tae S,W,O,L/file ? Q/quit')
       END IF
       CALL TASK_KLIN(line,kid,mode,wm_parm)
    END IF
    CALL mtx_broadcast1_character(kid)
    CALL mtx_broadcast1_integer(mode)

    IF(mode.EQ.2) CALL wm_broadcast
    IF(mode.NE.1) GOTO 1

    SELECT CASE(kid)
    CASE('P')                                      ! paramter input
       IF(nrank.EQ.0) CALL wm_parm(0,'WM',ierr)
       CALL wm_broadcast

    CASE('V')                                      ! view paramter
       IF(nrank.EQ.0) CALL wm_view

    CASE('R')                                      ! single/multi mode calc
       CALL wm_allocate
       CALL wm_loop(ierr)
       init=1

    CASE('D')                                      ! amp calc for given source
       READ(LINE(2:),*,ERR=1,END=1) nid
       SELECT CASE(nid)
       CASE(0)                                     ! given frequency
          CALL wm_am0d(kid,line)
       CASE(1)                                     ! real/imag frequency scan
          CALL wm_am1d(kid,line)
       CASE(2)                                     ! complex frequency scan
          CALL wm_am2d(kid,line)
       CASE(3)                                     ! parameter scan
          CALL wm_scan(kid,line)
       CASE DEFAULT
          WRITE(6,*) 'XX wm_menu: unknown nid for kid=R'
       END SELECT
       init=1

    CASE('F')                                      ! find complex root
       CALL wm_eign(kid,line)
       init=1

    CASE('G')                                      ! graph output
       IF(init.EQ.1.AND.nrank.EQ.0) CALL wm_gout

    CASE('S')
       IF(nrank.EQ.0) CALL wm_save(ierr)

    CASE('L')
       IF(nrank.EQ.0) CALL wm_load(ierr)

    CASE('W')
       IF(nrank.EQ.0) CALL wm_wout

    CASE('T')
       IF(nrank.EQ.0) CALL wm_taem

    CASE('?')
       WRITE(6,*) '# KID: P:  PARAMETER INPUT by namelist'
       WRITE(6,*) '       V:  VIEW PARAMETERS'
       WRITE(6,*) '       R:  WAVE EXCITED BY EXTERNAL ANTENNA'
       WRITE(6,*) '       D0: AMPLITUDE OF INTERNALLY EXCITED WAVE'
       WRITE(6,*) '       D1: FREQUENCY SCAN OF AMPLITUDE'
       WRITE(6,*) '       D2: COMPLEX FREQUENCY SCAN OF AMPLITUDE'
       WRITE(6,*) '       D3: PARAMETER SCAN OF EIGEN VALUE'
       WRITE(6,*) '       F:  FIND EIGEN VALUE'
       WRITE(6,*) '       G:  GRAPHICS'
       WRITE(6,*) '       S:  SAVE FIELD DATA'
       WRITE(6,*) '       L:  LOAD FIELD DATA'
       WRITE(6,*) '       W:  FILE OUTPUT'
       WRITE(6,*) '       O:  FILE OUTPUT for topics'
       WRITE(6,*) '       T:  TAE frequecy plot'
       WRITE(6,*) '       ?:  HELP'
       WRITE(6,*) '       Q:  QUIT'

    CASE('Z')
!       CALL wm_debug(ierr)

    CASE('Y')
!       CALL wm_test(CEFLD3D,CPABS3D,PABST3D,PABSTT3D, &
!                    MDM,NPHM,NRM,NSM, &
!                    NTHMAX,NPHMAX,NRMAX,NSMAX, &
!                    RGMIN,RGMAX,RAXIS,NCONT,NGRAPH)

    CASE('X','#','!')
       CONTINUE

    CASE('Q')
       GOTO 9000

    CASE DEFAULT
       IF(NRANK.EQ.0) WRITE(6,*) 'XX wm_menu: unknown kid'
    END SELECT
    GO TO 1

 9000 CONTINUE
    RETURN
  END SUBROUTINE wm_menu
END MODULE wmmenu
