C     $Id$
C
C   ************************************************
C   **  Grad-Shafranov Equation (Fixed Boundary)  **
C   ************************************************
C
      INCLUDE 'eqcomc.inc'
C
      WRITE(6,*) '## TASK/EQ 2004/11/08'
      CALL GSOPEN
      OPEN(7,STATUS='SCRATCH')
      CALL EQINIT
      CALL EQPARM(1,'eqparm',IERR)
C
      CALL EQMENU
C
      CLOSE(7)
      CALL GSCLOS
      STOP
      END
