C
C   ************************************************
C   **         Boundary shape function            **
C   ************************************************
C
      FUNCTION EQFBND(X)
C      
      INCLUDE 'eqcomc.h'
C
      EQFBND=ZBRF*COS(X+RDLT*SIN(X))-RKAP*SIN(X)
      RETURN
      END
C
C   ************************************************ 
C   **          　　P(psi)                        **
C   ************************************************
C
      REAL*8 FUNCTION PPSI(PSIN)
C
      INCLUDE 'eqcomc.h'
C
      IF(PSIN.LE.0.D0) THEN
         PPSI=0.D0
         RETURN
      ENDIF
      PPSI=PP0*(1.D0-(1.D0-PSIN)**PROFR0)**PROFP0
     &    +PP1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFP1
      PSIITB=1.D0-RHOITB**2
      IF(PSIN.GT.PSIITB) THEN
         PPSI=PPSI
     &       +PP2*(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**PROFP2
      ENDIF
      PPSI=PPSI*1.D6
      RETURN
      END
C
      REAL*8 FUNCTION DPPSI(PSIN)
C
      INCLUDE 'eqcomc.h'
C
      IF(PSIN.LE.0.D0) THEN
         DPPSI=0.D0
         RETURN
      ENDIF
      DPPSI=PP0*PROFP0*(1.D0-(1.D0-PSIN)**PROFR0)**(PROFP0-1.D0)
     &                        *PROFR0*(1.D0-PSIN)**(PROFR0-1.D0)
     &     +PP1*PROFP1*(1.D0-(1.D0-PSIN)**PROFR1)**(PROFP1-1.D0)
     &                        *PROFR1*(1.D0-PSIN)**(PROFR1-1.D0)
      PSIITB=1.D0-RHOITB**2
      IF(PSIN.GT.PSIITB) THEN
         DPPSI=DPPSI
     &     +PP2*PROFP2
     &     *(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**(PROFP2-1.D0)
     &     *PROFR2*((1.D0-PSIN)/(1.D0-PSIITB))**(PROFR2-1.D0)
     &     /(1.D0-PSIITB)
      ENDIF
      DPPSI=DPPSI*1.D6/PSI0
      RETURN
      END
C
C   ************************************************ 
C   **          　　J0(psi)                        **
C   ************************************************
C
      REAL*8 FUNCTION HJPSID(PSIN)
C
      INCLUDE 'eqcomc.h'
C
      IF(PSIN.LE.0.D0) THEN
         HJPSID=0.D0
         RETURN
      ENDIF
      HJPSID=-PJ0*(1.D0-(1.D0-PSIN)**PROFR0)**(PROFJ0+1.D0)
     &           *PSI0/(PROFR0*(PROFJ0+1.D0))
     &       -PJ1*(1.D0-(1.D0-PSIN)**PROFR1)**(PROFJ1+1.D0)
     &           *PSI0/(PROFR1*(PROFJ1+1.D0))
     &       -PJ2*(1.D0-(1.D0-PSIN)**PROFR2)**(PROFJ2+1.D0)
     &           *PSI0/(PROFR2*(PROFJ2+1.D0))
      RETURN
      END
C
      REAL*8 FUNCTION HJPSI(PSIN)
C
      INCLUDE 'eqcomc.h'
C
      IF(PSIN.LE.0.D0) THEN
         HJPSI=0.D0
         RETURN
      ENDIF
      HJPSI=-PJ0*(1.D0-(1.D0-PSIN)**PROFR0)**PROFJ0
     &                *(1.D0-PSIN)**(PROFR0-1.D0)
     &      -PJ1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFJ1
     &                *(1.D0-PSIN)**(PROFR1-1.D0)
     &      -PJ2*(1.D0-(1.D0-PSIN)**PROFR2)**PROFJ2
     &               *(1.D0-PSIN)**(PROFR2-1.D0)
      RETURN
      END
