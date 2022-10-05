! wmmain.f90

!*** /TASK/WM **********************************************************
!
!              ANALYSIS OF ICRF WAVE PROPAGATION AND ABSORPTION
!              ------------------------------------------------
!
!                                A. FUKUYAMA
!             Department of Nuclear Engineering, Kyoto University
!                      Sakyo-ku, Kyoto 606-8501, Japan
!                    Email: fukuyama@nucleng.kyoto-u.ac.jp
!                  URL: http://p-grp.nucleng.kyoto-u.ac.jp/wm/
!***********************************************************************

PROGRAM wm

  USE wmcomm
  USE plinit,ONLY: pl_init
  USE equnit,ONLY: eq_init
  USE dpinit,ONLY: dp_init
  USE dpparm,ONLY: dp_parm
  USE wminit,ONLY: wm_init
  USE wmparm,ONLY: wm_parm,wm_broadcast
  USE wmmenu,ONLY: wm_menu
  USE libmtx
  IMPLICIT NONE

  INTEGER:: ierr
  EXTERNAL:: gsopen,gsclos,eqinit

! --- Initialization

  CALL mtx_initialize
      
  IF(nrank.EQ.0) THEN
     WRITE(6,*) '##### /TASK/WM  2019-02-10 #####'
     CALL GSOPEN
     OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
  END IF

  CALL pl_init
  CALL eqinit
  CALL dp_init
  CALL wm_init

  IF(nrank.EQ.0) THEN
     CALL dp_parm(1,'dpparm',IERR)
     CALL wm_parm(1,'wmparm',IERR)
  ENDIF
  CALL wm_broadcast

  CALL wm_menu

  IF(nrank.EQ.0) THEN
     CALL GSCLOS
     CLOSE(7)
  END IF
  IF(nrank.EQ.0) CLOSE(30)
  IF(nrank.EQ.1) CLOSE(31)
  IF(nrank.EQ.2) CLOSE(32)
  IF(nrank.EQ.3) CLOSE(33)
  CALL mtx_finalize

  STOP
END PROGRAM wm
