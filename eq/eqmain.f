C     $Id$
C
C   ************************************************
C   **  Grad-Shafranov Equation (Fixed Boundary)  **
C   ************************************************
C
C   PSIN: 0 on axis, 1 on boundary
C
      INCLUDE 'eqcomc.inc'
C
      CHARACTER KPNAME*80
C
      WRITE(6,*) '## TASK/EQ 2004/11/08'
      CALL GSOPEN
      OPEN(7,STATUS='SCRATCH')
      CALL EQINIT
      KPNAME='eqparm'
      CALL EQPARF(KPNAME)
C
      CALL EQMENU
C
      CLOSE(7)
      CALL GSCLOS
      STOP
      END
