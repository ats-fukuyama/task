!     $Id: fpmenu.f90,v 1.3 2013/01/14 16:48:26 fukuyama Exp $

!     ***** TASK/FP MENU *****

      MODULE fpmenu

      contains

      SUBROUTINE fp_menu

      USE fpcomm
      USE fpinit
      USE fpprep
      USE fploop
      USE fpgout
      USE fpfout
      USE plinit
      USE fpfile
      USE libmpi
      USE libmtx

      IMPLICIT NONE
      CHARACTER(LEN=1)::  KID
      CHARACTER(LEN=80):: LINE
      integer:: ierr,NSA,NGRAPH_SAVE
      integer,DIMENSION(1):: mode
      REAL(4):: cputime1,cputime2

    1 CONTINUE
      IF(nrank.EQ.0) THEN
         ierr=0
         WRITE(6,601)
  601    FORMAT('## FP MENU: R:RUN C:CONT P,V:PARAM G,F:GRAPH', &
                           ' L,S:FILE Y:COEF Q:QUIT')
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
            CALL fp_prep(ierr)
            IF(MODEL_DISRUPT.ne.0.and.NRANK.eq.0)THEN
               OPEN(9,file="f1_1.dat") 
               open(10,file='time_evol.dat') 
               open(11,file='efield_e1.dat') 
               open(12,file='dndt.dat') 
               open(13,file='radial.dat') 
               open(14,file='nth-re.dat')
               open(15,file='re_pitch.dat')
               open(18,file='efield_ref.dat')
            END IF
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
         IF(NTG2.ne.0)THEN
            CALL fpsglb
            CALL fpwrtglb
            CALL fpsprf
            CALL fpwrtprf
         ELSE
            if(nrank.eq.0) WRITE(6,*) 'XX no data to write'
         END IF
      ELSEIF (KID.EQ.'Y') THEN
         TIMEFP=0.D0
         CALL fp_prep(ierr)
         IF(ierr.ne.0) GO TO 1
         DO NSA=1,NSAMAX
            CALL fp_coef(NSA)
         END DO
         CALL fpsglb
         CALL fpwrtglb
         CALL fpsprf
         CALL fpwrtprf
      ELSEIF (KID.EQ.'S') THEN
         if(nrank.eq.0) CALL fp_save2
         CALL mtx_barrier
      ELSEIF (KID.EQ.'L') THEN
         CALL FP_PRE_LOAD
         if(nrank.eq.0) CALL fp_load2
         CALL mtx_barrier
         CALL FP_POST_LOAD
      ELSEIF (KID.EQ.'Q') THEN
         IF(MODEL_DISRUPT.ne.0.and.NRANK.eq.0)THEN
            close(9)
            close(10)
            close(11)
            close(12)
            close(13)
            close(14)
            close(15)
            close(18)
         END IF
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
