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
      INCLUDE 'wmcomm.inc'
C
      CALL MPINIT(NPROCS,MYRANK)
      IF(NPROCS.LT.NCPUMIN) THEN
         WRITE(6,*) 'XX NPROCS.LT.NCPUMIN :',NPROCS,'.LT.',NCPUMIN
         STOP
      ENDIF
      IF(NPROCS.GT.NCPUMAX) THEN
         WRITE(6,*) 'XX NPROCS.GT.NCPUMAX :',NPROCS,'.GT.',NCPUMAX
         STOP
      ENDIF
C
      IF(MYRANK.EQ.0) THEN
         WRITE(6,*) '##### /TASK/WM  04/11/08 #####'
         CALL GSOPEN
      ENDIF
      CALL MPSYNC
C
      CALL PLINIT
      CALL DPINIT
      CALL WMINIT
      IF(MYRANK.EQ.0) THEN
         OPEN(7,STATUS='SCRATCH')
         CALL PLPARM(1,'plparm',IERR)
         CALL DPPARM(1,'dpparm',IERR)
         CALL WMPARM(1,'wmparm',IERR)
      ENDIF
      CALL MPSYNC
      CALL WMPRBC
C
      CALL WMMENU
C
      IF(MYRANK.EQ.0) CALL GSCLOS
      IF(MYRANK.EQ.0) CLOSE(7)
      CALL MPSYNC
C
      CALL MPTERM
C
      STOP
      END
