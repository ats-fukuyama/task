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
!                     V1.10  : 2000 NOV 25
!
!-----------------------------------------------------------------------

!      INCLUDE 'plcomm.inc'
      implicit none
      integer(4)  :: ier

      WRITE(6,*) '## TASK/PL 2004/11/08'
      OPEN(7,STATUS='SCRATCH')
      CALL PLINIT
      CALL PLPARM(1,'plparm',IER)

      CALL PLMENU

      CLOSE(7)
      STOP
      END
