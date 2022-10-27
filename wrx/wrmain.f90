! wrmain.f

!               ############# TASK/WR #############

!                        RAY TRACING
!                     BASED ON TASK/DP

!                           A. Fukuyama

!                Department of Nuclear Engineering
!                       Kyoto Univerisity
!                     Kyoto 606-8501, Japan

!                      V1.00  : 1988 FEB 10
!                      V2.00  : 1996 FEB 04
!                      V3.00  : 1997 AUG 05
!                      V3.10  : 2000 NOV 25
!                      V4.00  : 2018 SEP 14
!     --------------------------------------------------------

PROGRAM wr

  USE wrcomm
  USE plinit,ONLY: pl_init
  USE plparm,ONLY: pl_parm
  USE dpinit,ONLY: dp_init
  USE dpparm,ONLY: dp_parm
  USE wrinit,ONLY: wr_init
  USE wrparm,ONLY: wr_parm
  USE wrmenu,ONLY: wr_menu
  USE libmtx
  IMPLICIT NONE
  INTEGER:: ierr
  EXTERNAL GSOPEN,GSCLOS,EQINIT

      CALL mtx_initialize

      IF(NRANK.EQ.0) THEN
         WRITE(6,*) '## TASK/WR 2017/12/06'
         CALL GSOPEN
      END IF

      CALL pl_init
      CALL EQINIT
      CALL dp_init
      CALL wr_init
      IF(NRANK.EQ.0) THEN
         OPEN(7,STATUS='SCRATCH')
         CALL wr_parm(1,'wrparm',IERR)
      END IF

      CALL wr_menu

      IF(NRANK.EQ.0) THEN
         CLOSE(7)
         CALL GSCLOS
      END IF

      CALL mtx_finalize

      STOP
      END
