! fpmenu.f90

!     ***** TASK/FP MENU *****

MODULE fpmenu

CONTAINS

  SUBROUTINE fp_menu

    USE fpcomm
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
    USE fowcomm,ONLY: fow_allocate
    USE fowprep
    USE fowloop
    USE libkio
    USE libmpi
    USE libmtx

    IMPLICIT NONE
    CHARACTER(LEN=1)::  KID
    CHARACTER(LEN=80):: LINE
    INTEGER:: ierr,NGRAPH_SAVE
    INTEGER,DIMENSION(1):: mode
    REAL:: cputime1,cputime2

1   CONTINUE
    IF(nrank.EQ.0) THEN
       ierr=0
       WRITE(6,601)
601    FORMAT('## FP MENU: R:RUN C:CONT P,V:PARAM G,F:GRAPH', &
              ' L,S:FILE Y,Z:COEF W:WRITE Q:QUIT')
       CALL TASK_KLIN(LINE,KID,MODE(1),fp_parm)
    ENDIF
    CALL mtx_barrier
    CALL mtx_broadcast_character(KID,1)
    CALL mtx_broadcast_integer(MODE,1)
    IF(MODE(1).NE.1) GOTO 1

    CALL fp_broadcast

    IF(KID.EQ.'R'.OR.KID.EQ.'C') THEN
       IF(nrank.EQ.0) CALL CPU_TIME(cputime1)
       SELECT CASE(MODEL_FOW)
       CASE(0)
          IF(KID.EQ.'R') THEN
             CALL OPEN_EVOLVE_DATA_OUTPUT
             CALL fp_prep(ierr)
             IF(ierr.NE.0) GO TO 1
          ELSEIF(KID.EQ.'C')THEN
             CALL fp_continue(ierr)
             IF(ierr.NE.0) GO TO 1
          ENDIF
          CALL fp_loop
       CASE(1)
          IF(KID.EQ.'R') THEN
             CALL fp_prep(IERR)
             IF(ierr.NE.0) GO TO 1
             CALL fow_allocate
             CALL fow_prep(ierr)
             IF(ierr.NE.0) GO TO 1
          ELSEIF(KID.EQ.'C')THEN
             CALL fp_continue(ierr)
             IF(ierr.NE.0) GO TO 1
          END IF
          CALL fow_loop
       END SELECT
       IF(nrank.EQ.0) THEN
          CALL CPU_TIME(cputime2)
          WRITE(6,'(A,F12.3)') &
               '--cpu time =',cputime2-cputime1
       ENDIF
    ELSEIF (KID.EQ.'P') THEN
       IF(nrank.EQ.0) THEN
          CALL fp_parm(0,'FP',ierr)
       ENDIF
       CALL fp_broadcast
    ELSEIF (KID.EQ.'V') THEN
       IF(nrank.EQ.0) THEN
          CALL pl_view
          CALL fp_view
       ENDIF
    ELSEIF (KID.EQ.'G') THEN
       IF(nrank.EQ.0) CALL fp_gout
       CALL mtx_barrier
    ELSEIF (KID.EQ.'F') THEN
       IF(nrank.EQ.0) THEN
          NGRAPH_SAVE=NGRAPH
          NGRAPH=0
          CALL fp_gout
          NGRAPH=NGRAPH_SAVE
       ENDIF
       CALL mtx_barrier
    ELSEIF (KID.EQ.'W') THEN
       CALL FPWRTSNAP
       !         IF(NTG2.ne.0)THEN
       !            CALL fpsglb
       !            CALL fpwrtglb
       !            CALL fpsprf
       !            CALL fpwrtprf
       !         ELSE
       !            if(nrank.eq.0) WRITE(6,*) 'XX no data to write'
       !         END IF
    ELSEIF (KID.EQ.'Y') THEN
       TIMEFP=0.D0
       CALL fp_prep(ierr)
       IF(ierr.NE.0) GO TO 1
       CALL fp_coef(0)
       CALL fpsglb
       CALL fpwrtglb
       CALL fpsprf
       CALL fpwrtprf
    ELSEIF (KID.EQ.'S') THEN
       CALL FP_PRE_SAVE
       IF(nrank.EQ.0) CALL fp_save2
       CALL mtx_barrier
    ELSEIF (KID.EQ.'L') THEN
       CALL OPEN_EVOLVE_DATA_OUTPUT
       CALL FP_PRE_LOAD
       IF(nrank.EQ.0) CALL fp_load2
       CALL mtx_barrier
       CALL FP_POST_LOAD
    ELSEIF (KID.EQ.'Z') THEN
       CALL fp_caldeff
       CALL fp_calchieff
    ELSEIF (KID.EQ.'Q') THEN
       CALL CLOSE_EVOLVE_DATA_OUTPUT 
       GO TO 9000

    ELSE IF(KID.EQ.'X'.OR.KID.EQ.'#') THEN
       CONTINUE
    ELSE
       IF(nrank.EQ.0) WRITE(6,*) 'XX fp_menu: UNKNOWN KID'
    END IF
    
    GO TO 1

9000 RETURN
  END SUBROUTINE fp_menu

END MODULE fpmenu
