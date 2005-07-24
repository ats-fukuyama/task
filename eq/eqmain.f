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
      INCLUDE 'eqcomc.inc'
C
      WRITE(6,*) '## TASK/EQ 2005/07/24'
      CALL GSOPEN
      OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
C
      CALL EQINIT
      CALL EQPARM(1,'eqparm',IERR)
C
      CALL EQMENU
C
      CLOSE(7)
      CALL GSCLOS
      STOP
      END
