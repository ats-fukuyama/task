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
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MODELF.EQ.0) THEN
         PPSI=PP0*(1.D0-(1.D0-PSIN)**PROFR0)**PROFP0
     &       +PP1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFP1
         PSIITB=1.D0-RHOITB**2
         IF(PSIN.GT.PSIITB) THEN
            PPSI=PPSI
     &          +PP2*(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**PROFP2
         ENDIF
         PPSI=PPSI*1.D6
      ELSEIF(MODELF.EQ.2) THEN
         CALL SPL1DF(PSIN,PPSI,PSITR,UPPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX PPSI: SPL1DF : IERR=',IERR
         PPSI=PPSI*1.D6
      ENDIF
      RETURN
      END
C
      REAL*8 FUNCTION DPPSI(PSIN1)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MODELF.EQ.0) THEN
         DPPSI=PP0*PROFP0*(1.D0-(1.D0-PSIN)**PROFR0)**(PROFP0-1.D0)
     &                           *PROFR0*(1.D0-PSIN)**(PROFR0-1.D0)
     &        +PP1*PROFP1*(1.D0-(1.D0-PSIN)**PROFR1)**(PROFP1-1.D0)
     &                           *PROFR1*(1.D0-PSIN)**(PROFR1-1.D0)
         PSIITB=1.D0-RHOITB**2
         IF(PSIN.GT.PSIITB) THEN
            DPPSI=DPPSI
     &           +PP2*PROFP2
     &   *(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**(PROFP2-1.D0)
     &           *PROFR2*((1.D0-PSIN)/(1.D0-PSIITB))**(PROFR2-1.D0)
     &        /(1.D0-PSIITB)
         ENDIF
         DPPSI=DPPSI*1.D6/PSI0
      ELSEIF(MODELF.EQ.2) THEN
         CALL SPL1DD(PSIN,PPSI,DPPSI,PSITR,UPPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX DPPSI: SPL1DD : IERR=',IERR
         DPPSI=DPPSI*1.D6/PSI0
      ENDIF
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
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MODELF.EQ.0) THEN
         HJPSID=-PJ0*(1.D0-(1.D0-PSIN)**PROFR0)**(PROFJ0+1.D0)
     &              *PSI0/(PROFR0*(PROFJ0+1.D0))
     &          -PJ1*(1.D0-(1.D0-PSIN)**PROFR1)**(PROFJ1+1.D0)
     &              *PSI0/(PROFR1*(PROFJ1+1.D0))
     &          -PJ2*(1.D0-(1.D0-PSIN)**PROFR2)**(PROFJ2+1.D0)
     &              *PSI0/(PROFR2*(PROFJ2+1.D0))
         HJPSID=HJPSID*1.D6
      ELSEIF(MODELF.EQ.2) THEN
         CALL SPL1DI(PSIN,HJPSID,PSITR,UJPSI,UJPSI0,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX HJPSID: SPL1DF : IERR=',IERR
         HJPSID=HJPSID*1.D6
      ENDIF
      RETURN
      END
C
      REAL*8 FUNCTION HJPSI(PSIN1)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MODELF.EQ.0) THEN
         HJPSI=-PJ0*(1.D0-(1.D0-PSIN)**PROFR0)**PROFJ0
     &                   *(1.D0-PSIN)**(PROFR0-1.D0)
     &         -PJ1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFJ1
     &                   *(1.D0-PSIN)**(PROFR1-1.D0)
     &         -PJ2*(1.D0-(1.D0-PSIN)**PROFR2)**PROFJ2
     &                   *(1.D0-PSIN)**(PROFR2-1.D0)
         HJPSI=HJPSI*1.D6
      ELSEIF(MODELF.EQ.2) THEN
         CALL SPL1DF(PSIN,HJPSI,PSITR,UJPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX HJPSI: SPL1DF : IERR=',IERR
         HJPSI=HJPSI*1.D6
      ENDIF
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
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MODELF.EQ.0) THEN
         TPSI=PTS+(PT0-PTS)*(1.D0-(1.D0-PSIN)**PROFR0)**PROFT0
     &                 +PT1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFT1
         PSIITB=1.D0-RHOITB**2
         IF(PSIN.GT.PSIITB) THEN
            TPSI=TPSI
     &          +PT2*(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**PROFT2
         ENDIF
         TPSI=TPSI*1.D3*AEE
      ELSEIF(MODELF.EQ.2) THEN
         CALL SPL1DF(PSIN,TPSI,PSITR,UTPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TPSI: SPL1DF : IERR=',IERR
         TPSI=TPSI*1.D3*AEE
      ENDIF
      RETURN
      END
C
      REAL*8 FUNCTION DTPSI(PSIN1)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MODELF.EQ.0) THEN
         DTPSI=(PT0-PTS)
     &            *PROFT0*(1.D0-(1.D0-PSIN)**PROFR0)**(PROFT0-1.D0)
     &            *PROFR0*(1.D0-PSIN)**(PROFR0-1.D0)
     &        +PT1*PROFT1*(1.D0-(1.D0-PSIN)**PROFR1)**(PROFT1-1.D0)
     &            *PROFR1*(1.D0-PSIN)**(PROFR1-1.D0)
         PSIITB=1.D0-RHOITB**2
         IF(PSIN.GT.PSIITB) THEN
            DTPSI=DTPSI
     &           +PT2*PROFT2
     &        *(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**(PROFT2-1.D0)
     &        *PROFR2*((1.D0-PSIN)/(1.D0-PSIITB))**(PROFR2-1.D0)
     &         /(1.D0-PSIITB)
         ENDIF
         DTPSI=DTPSI*1.D3*AEE
      ELSEIF(MODELF.EQ.2) THEN
         CALL SPL1DD(PSIN,TPSI,DTPSI,PSITR,UTPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX DTPSI: SPL1DD : IERR=',IERR
         TPSI=TPSI*1.D3*AEE
      ENDIF
      RETURN
      END
C
C   ************************************************ 
C   **          　 　OMGPS(psi)                   **
C   ************************************************
C
      REAL*8 FUNCTION OMGPSI(PSIN1)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MODELF.EQ.0) THEN
         VPSI=PV0*(1.D0-(1.D0-PSIN)**PROFR0)**PROFV0
     &       +PV1*(1.D0-(1.D0-PSIN)**PROFR1)**PROFV1
         PSIITB=1.D0-RHOITB**2
         IF(PSIN.GT.PSIITB) THEN
            VPSI=VPSI
     &         +PV2*(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**PROFV2
         ENDIF
         OMGPSI=VPSI/RAXIS
      ELSEIF(MODELF.EQ.2) THEN
         CALL SPL1DF(PSIN,VPSI,PSITR,UVTPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX OMGPSI: SPL1DF : IERR=',IERR
         OMGPSI=VPSI/RAXIS
      ENDIF
      RETURN
      END
C
      REAL*8 FUNCTION DOMGPSI(PSIN1)
C
      INCLUDE 'eqcomc.h'
      INCLUDE 'eqcom4.h'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MODELF.EQ.0) THEN
         DVPSI=PV0*PROFV0*(1.D0-(1.D0-PSIN)**PROFR0)**(PROFV0-1.D0)
     &            *PROFR0*(1.D0-PSIN)**(PROFR0-1.D0)
     &        +PV1*PROFV1*(1.D0-(1.D0-PSIN)**PROFR1)**(PROFV1-1.D0)
     &                           *PROFR1*(1.D0-PSIN)**(PROFR1-1.D0)
         PSIITB=1.D0-RHOITB**2
         IF(PSIN.GT.PSIITB) THEN
            DVPSI=DVPSI
     &        +PV2*PROFV2
     &        *(1.D0-((1.D0-PSIN)/(1.D0-PSIITB))**PROFR2)**(PROFV2-1.D0)
     &       *PROFR2*((1.D0-PSIN)/(1.D0-PSIITB))**(PROFR2-1.D0)
     &        /(1.D0-PSIITB)
         ENDIF
         DOMGPSI=DVPSI/RAXIS
      ELSE IF(MODELF.EQ.2) THEN
         CALL SPL1DD(PSIN,VPSI,DVPSI,PSITR,UVTPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX DOMGPSI: SPL1DD : IERR=',IERR
         DOMGPSI=DVPSI/RAXIS
      ENDIF
      RETURN
      END
