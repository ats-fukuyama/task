!     *************************** TASK.TR ***************************
!     *                                                             *
!     *  RADIAL TRANSPORT IN A TOKAMAK                              *
!     *                                                             *
!     *  PROGRAMMED BY                                              *
!     *                                                             *
!     *        A. FUKUYAMA and M. HONDA                             *
!     *        DEPARTMENT OF NUCLEAR ENGINEERING                    *
!     *        GRADUATE SCHOOL OF ENGINEERING                       *
!     *        KYOTO UNIVERSITY                                     *
!     *        KYOTO 606-8501, JAPAN                                *
!     *                                                             *
!     *  ACKNOWLEDGEMENTS FOR                                       *
!     *                                                             *
!     *        T. KASAI                                             *
!     *        M. ISHIZU                                            *
!     *        S. TAKATSUKA                                         *
!     *                                                             *
!     *  PROGRAM VERSION                                            *
!     *                                                             *
!     *        V1.00 : 89/04/26                                     *
!     *        V1.10 : 89/05/20                                     *
!     *        V1.20 : 90/06/12                                     *
!     *        V1.30 : 90/08/18                                     *
!     *        V1.40 : 91/07/24                                     *
!     *        V2.00 : 93/06/18                                     *
!     *        V2.01 : 93/07/28                                     *
!     *        V2.02 : 93/08/05                                     *
!     *        V2.03 : 93/08/17                                     *
!     *        V2.10 : 97/08/23                                     *
!     *        V2.20 : 98/09/21                                     *
!     *        V2.21 : 00/04/19                                     *
!     *                                                             *
!     ***************************************************************

      USE TRCOMM, ONLY : GTCPU1, NFM, NGM, NRM, NSM, NTM
      use bpsd_mod
      use equnit_mod
      IMPLICIT NONE
      REAL(4)   :: GTCPU2
      INTEGER(4):: IERR

!     ------ INITIALIZATION ------

      WRITE(6,600)
  600 FORMAT(' ')
      WRITE(6,601) NRM,NSM,NFM,NGM,NTM
  601 FORMAT(' ','***** TASK/TR  04/12/10 *****', &
     &       '*** NRM, NSM, NFM, NGM, NTM  ***'/ &
     &       ' ',33X,I5,I5,I5,I5,I6)
      OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')

      CALL GSOPEN
      CALL GUTIME(GTCPU1)
!      bpsd_debug_flag = .FALSE.
!      call bpsd_init
      call eq_init
      CALL TRINIT
      CALL eq_parm(1,'eqparm',IERR)
      CALL TRPARM(1,'trparm',IERR)

      CALL TRMENU

      CALL GSCLOS
      CLOSE(7)
      CALL GUTIME(GTCPU2)
      WRITE(6,666) GTCPU2-GTCPU1
  666 FORMAT(' ','#      CPU TIME :   ',F8.3,' SEC ')
      STOP
      END PROGRAM
