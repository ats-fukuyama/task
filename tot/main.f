C     $Id$
C
C*** /TASK/TASK*********************************************************
C                                A. FUKUYAMA
C             Department of Nuclear Engineering, Kyoto University
C                      Sakyo-ku, Kyoto 606-8501, Japan
C                    Email: fukuyama@nucleng.kyoto-u.ac.jp
C                  URL: http://p-grp.nucleng.kyoto-u.ac.jp/wm/
C***********************************************************************
C
      INCLUDE '../mpi/mpicom.inc'
C
         WRITE(6,*) '##### test1 #####'
      CALL MPINIT(NPROCS,MYRANK)
         WRITE(6,*) '##### test2 #####', NPROCS, MYRANK
      IF(NPROCS.LT.NCPUMIN) THEN
         WRITE(6,*) 'XX NPROCS.LT.NCPUMIN :',NPROCS,'.LT.',NCPUMIN
         STOP
      ENDIF
      IF(NPROCS.GT.NCPUMAX) THEN
         WRITE(6,*) 'XX NPROCS.GT.NCPUMAX :',NPROCS,'.GT.',NCPUMAX
         STOP
      ENDIF
      WRITE(6,*) '##### test3 #####', NPROCS, MYRANK
C
      IF(MYRANK.EQ.0) THEN
         WRITE(6,*) '##### /TASK/TASK  04/11/08 #####'
         CALL GSOPEN
      ENDIF
      CALL MPSYNC
         WRITE(6,*) '##### test4 #####'
C
      CALL PLINIT
      CALL EQINIT
      CALL TRINIT
         WRITE(6,*) '##### test5 #####'
      CALL DPINIT
      CALL WRINIT
      CALL WMINIT
      CALL FPINIT
         WRITE(6,*) '##### test6 #####'
      IF(MYRANK.EQ.0) THEN
         OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
         CALL PLPARM(1,'plparm',IERR)
         CALL EQPARM(1,'eqparm',IERR)
         CALL TRPARM(1,'trparm',IERR)
         WRITE(6,*) '##### test7 #####'
         CALL DPPARM(1,'dpparm',IERR)
         CALL WRPARM(1,'wrparm',IERR)
         CALL WMPARM(1,'wmparm',IERR)
         CALL FPPARM(1,'fpparm',IERR)
      ENDIF
         WRITE(6,*) '##### test8 #####'
      CALL MPSYNC
      CALL WMPRBC
         WRITE(6,*) '##### test9 #####'
C
      CALL TASKMENU
C
      IF(MYRANK.EQ.0) CALL GSCLOS
      IF(MYRANK.EQ.0) CLOSE(7)
      CALL MPSYNC
C
      CALL MPTERM
C
      STOP
      END
C
C     ***** TASK MENU *****
C
      SUBROUTINE TASKMENU
C
      INCLUDE '../mpi/mpicom.inc'
C
      CHARACTER KID*2,KID1*1,KID2*1
C
    1 CONTINUE
         IF(MYRANK.EQ.0) THEN
            WRITE(6,601)
  601       FORMAT('## TASK MENU: EQ/TR/WR/WM/FP/DP/PL Q/QUIT')
            READ(5,'(A2)') KID
            KID1=KID(1:1)
            KID2=KID(2:2)
            CALL GUCPTL(KID1)
            CALL GUCPTL(KID2)
            KID=KID1//KID2
         ENDIF
         CALL MPBCKA(KID1)
         CALL MPBCKA(KID2)
         KID=KID1//KID2
C
         IF (KID.EQ.'EQ') THEN
            CALL EQMENU
         ELSE IF (KID.EQ.'TR') THEN
            CALL TRMENU
         ELSE IF (KID.EQ.'WR') THEN
            CALL WRMENU
         ELSE IF (KID.EQ.'WM') THEN
            CALL WMMENU
         ELSE IF (KID.EQ.'FP') THEN
            CALL FPMENU
         ELSE IF (KID.EQ.'DP') THEN
            CALL DPMENU
         ELSE IF (KID.EQ.'PL') THEN
            CALL PLMENU
         ELSE IF (KID1.EQ.'Q') THEN
            GOTO 9000
         ENDIF
      GO TO 1
C
 9000 RETURN
      END
