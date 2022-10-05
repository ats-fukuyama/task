! obmain.f90

!               ############# TASK/OB #############

!                          ORBIT TRACING

!                           A. Fukuyama

!                Department of Nuclear Engineering
!                       Kyoto Univerisity
!                     Kyoto 606-8501, Japan

!                      V1.00  : 2020 MAY 11
!     --------------------------------------------------------

PROGRAM ob

  USE obcomm
  USE plinit,ONLY: pl_init
  USE plparm,ONLY: pl_parm
  USE obinit,ONLY: ob_init
  USE obparm,ONLY: ob_parm
  USE obmenu,ONLY: ob_menu
  USE libmtx
  IMPLICIT NONE
  INTEGER:: ierr

      CALL mtx_initialize

      IF(NRANK.EQ.0) THEN
         WRITE(6,*) '## TASK/OB 2020/05/11'
         CALL GSOPEN
      END IF

      CALL pl_init
      CALL EQINIT
      CALL ob_init
      IF(NRANK.EQ.0) THEN
         OPEN(7,STATUS='SCRATCH')
         CALL ob_parm(1,'obparm',IERR)
      END IF

      CALL ob_menu

      IF(NRANK.EQ.0) THEN
         CLOSE(7)
         CALL GSCLOS
      END IF

      CALL mtx_finalize

      STOP
      END
