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
      USE plinit,ONLY:pl_init,pl_parm
      USE fpinit,ONLY:fp_init,fp_parm
      INCLUDE '../mpi/mpicom.inc'
C
      CALL MPINIT(nsize,nrank)
      IF(nsize.LT.NCPUMIN) THEN
         WRITE(6,*) 'XX nsize.LT.NCPUMIN :',nsize,'.LT.',NCPUMIN
         STOP
      ENDIF
      IF(nsize.GT.NCPUMAX) THEN
         WRITE(6,*) 'XX nsize.GT.NCPUMAX :',nsize,'.GT.',NCPUMAX
         STOP
      ENDIF
C
      IF(nrank.EQ.0) THEN
         WRITE(6,*) '##### /TASK/TASK  2013/01/20 #####'
         CALL GSOPEN
      ENDIF
      CALL MPSYNC
C
      CALL PL_INIT
      CALL EQINIT
      CALL TRINIT
      CALL DPINIT
      CALL WRINIT
      CALL WMINIT
      CALL FP_INIT

      IF(nrank.EQ.0) THEN
         OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
         CALL PL_PARM(1,'plparm',IERR)
         CALL EQPARM(1,'eqparm',IERR)
         CALL TRPARM(1,'trparm',IERR)
         CALL DPPARM(1,'dpparm',IERR)
         CALL WRPARM(1,'wrparm',IERR)
         CALL WMPARM(1,'wmparm',IERR)
         CALL FP_PARM(1,'fpparm',IERR)
      ENDIF
      CALL MPSYNC
      CALL WMPRBC
C
      CALL TASKMENU
C
      IF(nrank.EQ.0) CALL GSCLOS
      IF(nrank.EQ.0) CLOSE(7)
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
      USE plmenu,ONLY:pl_menu
      USE fpmenu,ONLY:fp_menu
      INCLUDE '../mpi/mpicom.inc'
C
      CHARACTER KID*2,KID1*1,KID2*1
C
    1 CONTINUE
         IF(nrank.EQ.0) THEN
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
            CALL FP_MENU
         ELSE IF (KID.EQ.'DP') THEN
            CALL DPMENU
         ELSE IF (KID.EQ.'PL') THEN
            CALL PL_MENU
         ELSE IF (KID1.EQ.'Q') THEN
            GOTO 9000
         ENDIF
      GO TO 1
C
 9000 RETURN
      END
