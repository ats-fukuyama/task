C     $Id$
C
C ************************************************************
C
C                    3D FOKKER-PLANCK CODE
C
C                  TASK/FP  V2.0  1992/12/12
C                           V2.1  1993/09/16
C                           V2.2  1997/03/18
C                           V3.0  1997/08/05
C
C                        PROGRAMMED BY
C                         A. FUKUYAMA
C                      OKAYAMA UNIVERSITY
C
C ************************************************************
C
      INCLUDE 'fpcomm.inc'
C
      CHARACTER KPNAME*80
C
      WRITE(6,*) '***** TASK/FP 2004/12/10 *****'
      CALL GSOPEN
      OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
C
      CALL PLINIT
      CALL EQINIT
      CALL FPINIT
      CALL PLPARM(1,'plparm',IERR)
      CALL EQPARM(1,'eqparm',IERR)
      CALL FPPARM(1,'fpparm',IERR)
C
      CALL FPMENU
C
      CLOSE(7)
      CALL GSCLOS
      STOP
      END
