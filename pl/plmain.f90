!     $Id$
!
!               ############# TASK/PL #############
!
!             BASIC PARAMETERS OF DEVICE AND PLASMAS
!
!                           A. Fukuyama

!                Department of Nuclear Engineering
!                       Kyoto Univerisity
!                     Kyoto 606-8501, Japan

!                     V1.00  : 1997 AUG 05
!                     V1.20  : 2009 JUL 21
!
!-----------------------------------------------------------------------

    PROGRAM plmain
      USE plcomm
      USE plinit,ONLY: pl_init
      USE plparm,ONLY: pl_parm
      USE plmenu,ONLY: pl_menu

      IMPLICIT none
      INTEGER(ikind)  :: ierr

      WRITE(6,*) '## TASK/PL 2019/09/01'
      OPEN(7,STATUS='SCRATCH')
      CALL gsopen
      CALL pl_init
      CALL pl_parm(1,'plparm',IERR)

      CALL pl_menu

      CALL gsclos
      CLOSE(7)
      STOP
    END PROGRAM plmain
