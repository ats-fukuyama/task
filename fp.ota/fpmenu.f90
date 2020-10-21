!     $Id: fpmenu.f90,v 1.3 2013/01/14 16:48:26 fukuyama Exp $

!     ***** TASK/FP MENU *****

      MODULE fpmenu

      contains

      SUBROUTINE fp_menu

      USE fpcomm
      USE fpinit
!      USE fpparm
      USE fpprep
      USE fploop
      USE fpgout
      USE fpfout
      USE plinit
      USE plparm
      USE fpfile
      USE fpdiffusivity
      USE fprad
      USE fowmenu
      USE libmpi
      USE libmtx

      IMPLICIT NONE
      CHARACTER(LEN=1)::  KID
      CHARACTER(LEN=80):: LINE
      integer:: ierr,NGRAPH_SAVE
      integer,DIMENSION(1):: mode
      real(8):: cputime1,cputime2

    1 CONTINUE
      IF(nrank.EQ.0) THEN
         ierr=0
         WRITE(6,601)
  601    FORMAT('## FP MENU: R:RUN C:CONT P,V:PARAM G,F:GRAPH', &
                           ' O:fow L,S:FILE Y,Z:COEF W:WRITE Q:QUIT')
         CALL TASK_KLIN(LINE,KID,MODE(1),fp_parm)
      ENDIF
      CALL mtx_barrier
      CALL mtx_broadcast_character(KID,1)
      CALL mtx_broadcast_integer(MODE,1)
      IF(MODE(1).NE.1) GOTO 1

      CALL fp_broadcast

      IF(KID.EQ.'R'.OR.KID.EQ.'C') THEN
         IF(nrank.EQ.0) CALL CPU_TIME(cputime1)
         IF(KID.EQ.'R') THEN
            CALL OPEN_EVOLVE_DATA_OUTPUT
            CALL fp_prep(ierr)
            IF(ierr.ne.0) GO TO 1
         ELSEIF(KID.eq.'C')THEN
            CALL fp_continue(ierr)
            IF(ierr.ne.0) GO TO 1
         ENDIF
         CALL fp_loop
         IF(nrank.eq.0) THEN
            CALL CPU_TIME(cputime2)
            write(6,'(A,F12.3)') &
                 '--cpu time =',cputime2-cputime1
         ENDIF
      ELSEIF (KID.EQ.'P') THEN
         if(nrank.eq.0) then
            CALL fp_parm(0,'FP',ierr)
         endif
         CALL fp_broadcast
      ELSEIF (KID.EQ.'V') THEN
         if(nrank.eq.0) then
            CALL pl_view
            CALL fp_view
         endif
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
         IF(ierr.ne.0) GO TO 1
         CALL fp_coef(0)
         CALL fpsglb
         CALL fpwrtglb
         CALL fpsprf
         CALL fpwrtprf
      ELSEIF (KID.EQ.'S') THEN
         CALL FP_PRE_SAVE
         if(nrank.eq.0) CALL fp_save2
         CALL mtx_barrier
      ELSEIF (KID.EQ.'L') THEN
         CALL OPEN_EVOLVE_DATA_OUTPUT
         CALL FP_PRE_LOAD
         if(nrank.eq.0) CALL fp_load2
         CALL mtx_barrier
         CALL FP_POST_LOAD
      ELSEIF (KID.EQ.'Z') THEN
        if(nrank.eq.0)then
         CALL fp_deff
         CALL fp_chieff

         !CALL fp_radiation
       end if
       call mtx_barrier
    ELSEIF (KID.EQ.'O') THEN
       CALL fow_menu
    ELSEIF (KID.EQ.'Q') THEN
         CALL CLOSE_EVOLVE_DATA_OUTPUT
         GO TO 9000

      ELSE IF(KID.EQ.'X'.OR.KID.EQ.'#') THEN
         CONTINUE
      ELSE
         IF(nrank.eq.0) WRITE(6,*) 'XX fp_menu: UNKNOWN KID'
      END IF

      GO TO 1

 9000 RETURN
      END SUBROUTINE fp_menu

      END MODULE fpmenu
