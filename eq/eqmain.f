C     $Id$
C
C   ************************************************
C   **  Grad-Shafranov Equation (Fixed Boundary)  **
C
C   sigma  : normalized radial coordinate [0,1]
C   theta  : periodic poloidal coodinated [0,2*Pi)
C   PSI(sigma,theta) : poloidal magnetic flux
C                      PSI0 at (0,theta), 0 at (1,theta)
C   PSITA  : toroidal magnetic flux on surface R=RA
C   PSIRZ(R,Z)       : poloidal magnetic flux
C                      PSI0 at (Raxis,Zaxis), 0 on surface
C
C   ************************************************
C
      USE plinit,ONLY: pl_init
      USE plparm,ONLY: pl_parm
      INCLUDE 'eqcomc.inc'
C
      WRITE(6,*) '## TASK/EQ 2009/09/01'
      CALL GSOPEN
      OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
C
      CALL PL_INIT
      CALL EQINIT
      CALL PL_PARM(1,'plparm',IERR)
      CALL EQPARM(1,'eqparm',IERR)
C
      CALL EQMENU
C
      CLOSE(7)
      CALL GSCLOS
      STOP
      END
