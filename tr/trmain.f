C     $Id$
C
C     *************************** TASK.TR ***************************
C     *                                                             *
C     *  RADIAL TRANSPORT IN A TOKAMAK                              *
C     *                                                             *
C     *  PROGRAMMED BY                                              *
C     *                                                             *
C     *        A. FUKUYAMA                                          *
C     *        DEPARTMENT OF NUCLEAR ENGINEERING                    *
C     *        GRADUATE SCHOOL OF ENGINEERING                       *
C     *        KYOTO UNIVERSITY                                     *
C     *        KYOTO 606-8501, JAPAN                                *
C     *                                                             *
C     *  ACKNOWLEDGEMENTS FOR                                       *
C     *                                                             *
C     *        T. KASAI                                             *
C     *        M. ISHIZU                                            *
C     *        S. TAKATSUKA                                         *
C     *                                                             *
C     *  PROGRAM VERSION                                            *
C     *                                                             *
C     *        V1.00 : 89/04/26                                     *
C     *        V1.10 : 89/05/20                                     *
C     *        V1.20 : 90/06/12                                     *
C     *        V1.30 : 90/08/18                                     *
C     *        V1.40 : 91/07/24                                     *
C     *        V2.00 : 93/06/18                                     *
C     *        V2.01 : 93/07/28                                     *
C     *        V2.02 : 93/08/05                                     *
C     *        V2.03 : 93/08/17                                     *
C     *        V2.10 : 97/08/23                                     *
C     *        V2.20 : 98/09/21                                     *
C     *        V2.21 : 00/04/19                                     *
C     *                                                             *
C     ***************************************************************
C
      INCLUDE 'trcomm.inc'
C
C     ------ INITIALIZATION ------
C
      WRITE(6,600)
  600 FORMAT(' ')
      WRITE(6,601) NRM,NSM,NFM,NGM,NTM
  601 FORMAT(' ','***** TASK/TR  04/12/10 *****',
     &       '*** NRM, NSM, NFM, NGM, NTM  ***'/
     &       ' ',33X,I5,I5,I5,I5,I6)
      OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
C
      CALL GSOPEN
      CALL GUTIME(GTCPU1)
      CALL TRINIT
      CALL TRPARM(1,'trparm',IERR)
C
      CALL TRMENU
C
      CALL GSCLOS
      CLOSE(7)
      CALL GUTIME(GTCPU2)
      WRITE(6,666) GTCPU2-GTCPU1
  666 FORMAT(' ','#      CPU TIME :   ',F8.3,' SEC ')
      STOP
      END
