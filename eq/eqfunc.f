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
      REAL*8 FUNCTION PPSI(PSIN1)
C
      INCLUDE 'eqcomc.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
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
      REAL*8 FUNCTION DPPSI(PSIN1)
C
      INCLUDE 'eqcomc.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
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
      REAL*8 FUNCTION HJPSID(PSIN1)
C
      INCLUDE 'eqcomc.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
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
      REAL*8 FUNCTION HJPSI(PSIN1)
C
      INCLUDE 'eqcomc.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
      HJPSI=-PJ0*(1.D0-(1.D0-PSIN)**PROFR0)**PROFJ0
     &                *(1.D0-PSIN)**(PROFR0-1.D0)
     &      -PJ1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFJ1
     &                *(1.D0-PSIN)**(PROFR1-1.D0)
     &      -PJ2*(1.D0-(1.D0-PSIN)**PROFR2)**PROFJ2
     &               *(1.D0-PSIN)**(PROFR2-1.D0)
      RETURN
      END
C
C   ************************************************ 
C   **          　 　T(psi)                       **
C   ************************************************
C
      REAL*8 FUNCTION TPSI(PSIN1)
C
      INCLUDE 'eqcomc.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C      TPSI=1.D0+0.34935013D0*PSIN-7.9768279D0*PSIN**2
C     &       +65.503148D0*PSIN**3-247.45634D0*PSIN**4
C     &       +436.68766D0*PSIN**5-292.12510D0*PSIN**6
C      PSIITB=0.65D0
C      IF(TPSI.GT.1.D0) THEN
C         TPSI=1.D0
C      ENDIF
C      IF(PSIN.GT.PSIITB) THEN
C         TPSI=TS+(T0*0.5D0-TS)*(1.D0-PSIN**2)
C      ENDIF
      TPSI=PTS+(PT0-PTS)*(1.D0-(1.D0-PSIN)**PROFR0)**PROFT0
     &              +PT1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFT1
      PSIITB=1.D0-RHOITB**2
C      OLDTPS=TPSI
      IF(PSIN.GT.PSIITB) THEN
        TPSI=TPSI
     &      +PT2*(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**PROFT2
      ENDIF
C      write(6,*) OLDTPS,TPSI-OLDTPS
      TPSI=TPSI*1.D3*AEE/AMP
      RETURN
      END
C
      REAL*8 FUNCTION DTPSI(PSIN1)
C
      INCLUDE 'eqcomc.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C      DTPSI=0.34935013D0-7.9768279D0*2.D0*PSIN
C     &     +65.503148D0*3.D0*PSIN**2-247.45634D0*4.D0*PSIN**3
C     &     +436.68766D0*5.D0*PSIN**4-292.12510D0*6.D0*PSIN**5
C      PSIITB=0.65D0
C      IF(PSIN.GT.PSIITB) THEN
C         DTPSI=-2.D0*(T0*0.5D0-TS)*PSIN
C      ENDIF
      DTPSI=(PT0-PTS)*PROFT0*(1.D0-(1.D0-PSIN)**PROFR0)**(PROFT0-1.D0)
     &                        *PROFR0*(1.D0-PSIN)**(PROFR0-1.D0)
     &           +PT1*PROFT1*(1.D0-(1.D0-PSIN)**PROFR1)**(PROFT1-1.D0)
     &                        *PROFR1*(1.D0-PSIN)**(PROFR1-1.D0)
      PSIITB=1.D0-RHOITB**2
      IF(PSIN.GT.PSIITB) THEN
         DTPSI=DTPSI
     &           +PT2*PROFT2
     &     *(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**(PROFT2-1.D0)
     &     *PROFR2*((1.D0-PSIN)/(1.D0-PSIITB))**(PROFR2-1.D0)
     &     /(1.D0-PSIITB)
      ENDIF
      DTPSI=DTPSI*1.D3*AEE/AMP
      RETURN
      END
C
C   ************************************************ 
C   **          　 　OMGPS(psi)                   **
C   ************************************************
C
      REAL*8 FUNCTION OMGPS(PSIN1)
C
      INCLUDE 'eqcomc.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
      OMGPS=PV0*(1.D0-(1.D0-PSIN)**PROFR0)**PROFV0
     &     +PV1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFV1
      PSIITB=1.D0-RHOITB**2
      IF(PSIN.GT.PSIITB) THEN
         OMGPS=OMGPS
     &     +PV2*(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**PROFV2
      ENDIF
      OMGPS=OMGPS/RAXIS
      RETURN
      END
C
C     DOMGPS not normalized and not multiplyed by OMGPS0
C
      REAL*8 FUNCTION DOMGPS(PSIN1)
C
      INCLUDE 'eqcomc.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
      DOMGPS=PV0*PROFV0*(1.D0-(1.D0-PSIN)**PROFR0)**(PROFV0-1.D0)
     &          *PROFR0*(1.D0-PSIN)**(PROFR0-1.D0)
     &      +PV1*PROFV1*(1.D0-(1.D0-PSIN)**PROFR1)**(PROFV1-1.D0)
     &                         *PROFR1*(1.D0-PSIN)**(PROFR1-1.D0)
      PSIITB=1.D0-RHOITB**2
      IF(PSIN.GT.PSIITB) THEN
         DOMGPS=DOMGPS
     &     +PV2*PROFV2
     &     *(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**(PROFV2-1.D0)
     &     *PROFR2*((1.D0-PSIN)/(1.D0-PSIITB))**(PROFR2-1.D0)
     &     /(1.D0-PSIITB)
      ENDIF
      DOMGPS=DOMGPS/RAXIS
      RETURN
      END
