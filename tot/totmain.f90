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
  USE equnit,ONLY:eq_init
  USE trcomm,ONLY:open_trcomm
  USE trinit,ONLY:tr_init
  USE trparm,ONLY:tr_parm
  USE ticomm,ONLY:open_ticomm_parm
  USE tiinit,ONLY:ti_init
  USE tiparm,ONLY:ti_parm,ti_broadcast
  USE fpinit,ONLY:fp_init
  USE fpparm,ONLY:fp_parm
  USE fpcomm,ONLY:fpcomm_parm
  USE dpinit,ONLY:dp_init
  USE dpparm,ONLY:dp_parm
  USE wrcomm,ONLY:open_wrcomm_parm
  USE wrinit,ONLY:wr_init
  USE wrparm,ONLY:wr_parm
  USE wminit,ONLY:wm_init
  USE wmparm,ONLY:wm_parm,wm_broadcast
  USE commpi
  USE libmtx
  USE totmenu,ONLY: tot_menu

  CALL open_trcomm
  CALL open_ticomm_parm
!  CALL open_fpcomm_parm
  CALL open_wrcomm_parm

  CALL mtx_initialize

  IF(nrank.EQ.0) THEN
     WRITE(6,*) '##### /TASK/TOT  2019/02/25 #####'
     CALL GSOPEN
     OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
  ENDIF

  CALL pl_init
  CALL eq_init
  CALL tr_init
  CALL dp_init
  CALL wr_init
  CALL wm_init
  CALL fp_init
  CALL ti_init
  

  IF(nrank.EQ.0) THEN
     CALL pl_parm(1,'plparm',IERR)
     CALL eqparm(1,'eqparm',IERR)
     CALL tr_parm(1,'trparm',IERR)
     CALL dp_parm(1,'dpparm',IERR)
     CALL wr_parm(1,'wrparm',IERR)
     CALL wm_parm(1,'wmparm',IERR)
     CALL fp_parm(1,'fpparm',IERR)
     CALL ti_parm(1,'tiparm',IERR)
  ENDIF
  CALL wm_broadcast
  CALL ti_broadcast

  CALL tot_menu

  IF(nrank.EQ.0) THEN
     CALL GSCLOS
     CLOSE(7)
  END IF

  CALL mtx_finalize

  STOP
END PROGRAM tot
