C     $Id$
C
C*** /TASK/WM **********************************************************
C
C              ANALYSIS OF ICRF WAVE PROPAGATION AND ABSORPTION
C              ------------------------------------------------
C
C                                A. FUKUYAMA
C             Department of Nuclear Engineering, Kyoto University
C                      Sakyo-ku, Kyoto 606-8501, Japan
C                    Email: fukuyama@nucleng.kyoto-u.ac.jp
C                  URL: http://p-grp.nucleng.kyoto-u.ac.jp/wm/
C***********************************************************************
C
      USE plinit,ONLY:pl_init
      USE plparm,ONLY:pl_parm
      USE dpinit,ONLY:dp_init
      USE dpparm,ONLY:dp_parm
      INCLUDE 'wmcomm.inc'
C
      CALL mtx_initialize

      IF(NRANK.EQ.0) THEN
         WRITE(6,*) '##### /TASK/WM  05/06/22 #####'
         CALL GSOPEN
      ENDIF
      CALL mtx_barrier
C
      CALL PL_INIT
      CALL EQINIT
      CALL DP_INIT
      CALL WMINIT
      IF(NRANK.EQ.0) THEN
         OPEN(7,STATUS='SCRATCH')
         CALL PL_PARM(1,'plparm',IERR)
         CALL EQPARM(1,'eqparm',IERR)
         CALL DP_PARM(1,'dpparm',IERR)
         CALL WMPARM(1,'wmparm',IERR)
      ENDIF
      CALL mtx_barrier
      CALL WMPRBC
C
      CALL WMMENU
C
      IF(NRANK.EQ.0) CALL GSCLOS
      IF(NRANK.EQ.0) CLOSE(7)
      CALL mtx_barrier
C
      CALL mtx_finalize
C
      STOP
      END
