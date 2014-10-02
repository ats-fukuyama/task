!     $Id$

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
      INTEGER:: IERR,NSA,NGRAPH_SAVE,NR,NTH,NP
      REAL(rkind):: P1,P2,FL1
      INTEGER,DIMENSION(1):: mode
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
         if(nrank.eq.0) CALL fp_save
         CALL mtx_barrier
      ELSEIF (KID.EQ.'L') THEN
         if(nrank.eq.0) CALL fp_load
         CALL mtx_barrier

      ELSEIF (KID.EQ.'T') THEN
         IF(nrank.EQ.0) THEN
            CALL FROPEN(49,'FNS01.DAT',1,0,'FNS',IERR)
            NSAMAX=1
            NSBMAX=1
            NSB_NSA(1)=1
            PMAX(1)=23
            NGLINE=30
            NRMAX=1
            NTHMAX=64
            NPMAX=1024
            NTMAX=0
            TIMEFP=0.D0
            NTG1=0
            NTG2=0
            CALL fp_allocate
            call fp_allocate_ntg1
            call fp_allocate_ntg2
            READ(49,*)
            READ(49,*)
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  READ(49,*) PM(NP,1),THM(NTH),P1,P2,FNS(NTH,NP,1,1),FL1
               END DO
            END DO
            CLOSE(49)
         END IF
         CALL mtx_barrier

      ELSEIF (KID.EQ.'Q') THEN
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
