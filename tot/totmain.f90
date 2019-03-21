! totmain.f90

!*** /TASK/TASK*********************************************************
!                                A. FUKUYAMA
!             Department of Nuclear Engineering, Kyoto University
!                      Sakyo-ku, Kyoto 606-8501, Japan
!                    Email: fukuyama@nucleng.kyoto-u.ac.jp
!                  URL: http://p-grp.nucleng.kyoto-u.ac.jp/wm/
!***********************************************************************

PROGRAM tot

  USE plinit,ONLY:pl_init
  USE plparm,ONLY:pl_parm
  USE equnit_mod,ONLY:eq_init
  USE fpinit,ONLY:fp_init,fp_parm
  USE tiinit,ONLY:ti_init
  USE tiparm,ONLY:ti_parm,ti_broadcast
  USE commpi
  USE libmtx
  USE totmenu

  CALL mtx_initialize

  IF(nrank.EQ.0) THEN
     WRITE(6,*) '##### /TASK/TOT  2019/02/25 #####'
     CALL GSOPEN
     OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
  ENDIF

  CALL pl_init
  CALL eq_init
  CALL trinit
  CALL dpinit
  CALL wrinit
  CALL wminit
  CALL fp_init
  CALL ti_init
  

  IF(nrank.EQ.0) THEN
     CALL pl_parm(1,'plparm',IERR)
     CALL eqparm(1,'eqparm',IERR)
     CALL trparm(1,'trparm',IERR)
     CALL dpparm(1,'dpparm',IERR)
     CALL wrparm(1,'wrparm',IERR)
     CALL wmparm(1,'wmparm',IERR)
     CALL fp_parm(1,'fpparm',IERR)
     CALL ti_parm(1,'tiparm',IERR)
  ENDIF
  CALL wmprbc
  CALL ti_broadcast

  CALL tot_menu

  IF(nrank.EQ.0) THEN
     CALL GSCLOS
     CLOSE(7)
  END IF

  CALL mtx_finalize

  STOP
END PROGRAM tot
