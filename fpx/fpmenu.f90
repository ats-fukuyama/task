! fpmenu.f90

MODULE fpmenu

  PRIVATE
  PUBLIC fp_menu

CONTAINS

!     ***** TASK/FP MENU *****

  SUBROUTINE fp_menu

    USE fpcomm
    USE fowcomm
    USE fpinit
    USE fpparm
    USE fpprep
    USE fploop
    USE fpgout
    USE fpfout
    USE plinit
    USE plparm
    USE fpfile
    USE fpcaltp
    USE fpcalte
    USE fpcaldeff
    USE fpcalchieff
    use fowprep
    use foworbit
    use fowdistribution
    use fowloop
    use fowevalNC
    USE libkio
    USE libmpi
    USE libmtx

    IMPLICIT NONE
    CHARACTER(LEN=1)::  KID
    CHARACTER(LEN=80):: LINE
    INTEGER:: ierr,NGRAPH_SAVE
    INTEGER:: mode
    REAL:: cputime1,cputime2

1   CONTINUE
    IF(nrank.EQ.0) THEN
       ierr=0
       WRITE(6,601)
601    FORMAT('## FP MENU: R:RUN C:CONT P,V:PARAM G,F:GRAPH', &
            'L,S:FILE W,Y,Z:WRITE N:NC Q:QUIT')
       CALL TASK_KLIN(LINE,KID,MODE,fp_parm)
    ENDIF
    CALL mtx_barrier
    CALL mtx_broadcast_character(KID,1)
    CALL mtx_broadcast1_integer(MODE)
    IF(MODE.NE.1) GOTO 1

    CALL fp_broadcast

    SELECT CASE(KID)
    CASE('R')
       IF(nrank.EQ.0) CALL CPU_TIME(cputime1)
       CALL OPEN_EVOLVE_DATA_OUTPUT
       CALL fp_prep(ierr)
       IF(ierr.ne.0) GO TO 1
       IF(MODEL_FOW.EQ.0) THEN
          CALL fp_loop
       ELSE
          call fow_allocate
          call fow_prep
          call fow_loop
       END IF
       IF(nrank.eq.0) THEN
          CALL CPU_TIME(cputime2)
          write(6,'(A,F12.3)') &
               '--cpu time =',cputime2-cputime1
       ENDIF
    CASE('C')
       IF(nrank.EQ.0) CALL CPU_TIME(cputime1)
       CALL fp_continue(ierr)
       IF(ierr.ne.0) GO TO 1
       IF(MODEL_FOW.EQ.0) THEN
          CALL fp_loop
       ELSE
          CALL fow_loop
       END IF
       IF(nrank.eq.0) THEN
          CALL CPU_TIME(cputime2)
          write(6,'(A,F12.3)') &
               '--cpu time =',cputime2-cputime1
       ENDIF
    CASE('P')
       if(nrank.eq.0) then
          CALL fp_parm(0,'FP',ierr)
       endif
       CALL fp_broadcast
    CASE('V')
       if(nrank.eq.0) then
          CALL pl_view
          CALL fp_view
       endif
    CASE('G')
       IF(nrank.EQ.0) CALL fp_gout
       CALL mtx_barrier
    CASE('F')
       IF(nrank.EQ.0) THEN
          NGRAPH_SAVE=NGRAPH
          NGRAPH=0
          CALL fp_gout
          NGRAPH=NGRAPH_SAVE
       ENDIF
       CALL mtx_barrier
    CASE('W')
       CALL FPWRTSNAP
    CASE('N')
       IF(MODEL_FOW.NE.0) CALL output_neoclass
    CASE('Y')
       TIMEFP=0.D0
       CALL fp_prep(ierr)
       IF(ierr.ne.0) GO TO 1
       CALL fp_coef(0)
       CALL fpsglb
       CALL fpwrtglb
       CALL fpsprf
       CALL fpwrtprf
    CASE('S')
       CALL FP_PRE_SAVE
       if(nrank.eq.0) CALL fp_save2
       CALL mtx_barrier
    CASE('L')
       CALL OPEN_EVOLVE_DATA_OUTPUT
       CALL FP_PRE_LOAD
       if(nrank.eq.0) CALL fp_load2
       CALL mtx_barrier
       CALL FP_POST_LOAD
    CASE('Z')
       CALL fp_caldeff
       CALL fp_calchieff
    CASE('Q')
       CALL CLOSE_EVOLVE_DATA_OUTPUT
       GO TO 9000
       
    CASE('X','#','!')
       CONTINUE
    CASE DEFAULT
       IF(nrank.eq.0) WRITE(6,*) 'XX fp_menu: UNKNOWN KID'
    END SELECT
    
    GO TO 1

9000 CONTINUE
    IF(id_fp_allocate.EQ.1) THEN
!       CALL fp_deallocate
!       CALL fp_deallocate_ntg1
!       CALL fp_deallocate_ntg2
       id_fp_allocate=0
    END IF
    RETURN
  END SUBROUTINE fp_menu

END MODULE fpmenu
