C
C   ************************************************
C   **         Boundary shape function            **
C   ************************************************
C
      FUNCTION EQFBND(X)
C      
      INCLUDE 'eqcomc.inc'
C
      EQFBND=ZBRF*COS(X+RDLT*SIN(X))-RKAP*SIN(X)
      RETURN
      END
C
C   ************************************************ 
C   **                  P(psi)                    **
C   ************************************************
C
      REAL*8 FUNCTION PPSI(PSIN1)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MDLEQF.LT.5) THEN
         PPSI=PP0*(1.D0-PSIN**PROFR0)**PROFP0
     &       +PP1*(1.D0-PSIN**PROFR1)**PROFP1
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            PPSI=PPSI
     &          +PP2*(1.D0-(PSIN/PSIITB)**PROFR2)**PROFP2
         ENDIF
         PPSI=PPSI*1.D6
      ELSE
         CALL SPL1DF(PSIN,PPSI,PSITRX,UPPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX PPSI: SPL1DF : IERR=',IERR
         IF(PPSI.LT.0.D0) PPSI=0.D0
         PPSI=PPSI*1.D6
      ENDIF
C      write(6,*) PPSI
      RETURN
      END
C
      REAL*8 FUNCTION DPPSI(PSIN1)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MDLEQF.LT.5) THEN
         DPPSI=-PP0*PROFP0*(1.D0-PSIN**PROFR0)**(PROFP0-1.D0)
     &                            *PROFR0*PSIN**(PROFR0-1.D0)
     &         -PP1*PROFP1*(1.D0-PSIN**PROFR1)**(PROFP1-1.D0)
     &                            *PROFR1*PSIN**(PROFR1-1.D0)
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            DPPSI=DPPSI
     &           -PP2*PROFP2
     &           *(1.D0-(PSIN/PSIITB)**PROFR2)**(PROFP2-1.D0)
     &                   *PROFR2*(PSIN/PSIITB)**(PROFR2-1.D0)
     &           /PSIITB
         ENDIF
         DPPSI=DPPSI*1.D6
      ELSE
         CALL SPL1DD(PSIN,PPSI,DPPSI,PSITRX,UPPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX DPPSI: SPL1DD : IERR=',IERR
         DPPSI=DPPSI*1.D6
      ENDIF
      RETURN
      END
C
C   ************************************************ 
C   **                  F(psi)=B*R                **
C   ************************************************
C
      REAL*8 FUNCTION FPSI(PSIN1)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MDLEQF.LT.5) THEN
         FPSI=BB*RR
     &       +FF0*(1.D0-PSIN**PROFR0)**PROFF0
     &       +FF1*(1.D0-PSIN**PROFR1)**PROFF1
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            FPSI=FPSI
     &          +FF2*(1.D0-(PSIN/PSIITB)**PROFR2)**PROFF2
         ENDIF
      ELSE
         CALL SPL1DF(PSIN,FPSI,PSITRX,UFPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX FPSI: SPL1DF : IERR=',IERR
      ENDIF
      RETURN
      END
C
      REAL*8 FUNCTION DFPSI(PSIN1)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MDLEQF.LT.5) THEN
         DFPSI=-FF0*PROFF0*(1.D0-PSIN**PROFR0)**(PROFF0-1.D0)
     &                            *PROFR0*PSIN**(PROFR0-1.D0)
     &         -FF1*PROFF1*(1.D0-PSIN**PROFR1)**(PROFF1-1.D0)
     &                            *PROFR1*PSIN**(PROFR1-1.D0)
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            DFPSI=DFPSI
     &           -FF2*PROFF2
     &           *(1.D0-(PSIN/PSIITB)**PROFR2)**(PROFF2-1.D0)
     &                   *PROFR2*(PSIN/PSIITB)**(PROFR2-1.D0)
     &           /PSIITB
         ENDIF
      ELSE
         CALL SPL1DD(PSIN,FPSI,DFPSI,PSITRX,UFPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX DFPSI: SPL1DD : IERR=',IERR
      ENDIF
      RETURN
      END
C
C   ************************************************ 
C   **                 Q(psi)                     **
C   ************************************************
C
      REAL*8 FUNCTION QPSI(PSIN1)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MDLEQF.LT.5) THEN
         QPSI=QQS
     &       +(QQ0-QQS)*(1.D0-PSIN**PROFR0)**PROFQ0
     &       +QQ1*(1.D0-PSIN**PROFR1)**PROFQ1
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            QPSI=QPSI
     &          +QQ2*(1.D0-(PSIN/PSIITB)**PROFR2)**PROFQ2
         ENDIF
      ELSE
         CALL SPL1DF(PSIN,QPSI,QSITRX,UQPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX QPSI: SPL1DF : IERR=',IERR
         IF(QPSI.LT.0.D0) QPSI=0.D0
      ENDIF
C      write(6,*) QPSI
      RETURN
      END
C
      REAL*8 FUNCTION DQPSI(PSIN1)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MDLEQF.LT.5) THEN
         DQPSI=-(QQ0-QQS)*PROFQ0*(1.D0-PSIN**PROFR0)**(PROFQ0-1.D0)
     &                            *PROFR0*PSIN**(PROFR0-1.D0)
     &         -QQ1*PROFQ1*(1.D0-PSIN**PROFR1)**(PROFQ1-1.D0)
     &                            *PROFR1*PSIN**(PROFR1-1.D0)
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            DQPSI=DQPSI
     &           +QQ2*PROFQ2
     &           *(1.D0-(PSIN/PSIITB)**PROFR2)**(PROFQ2-1.D0)
     &                   *PROFR2*(PSIN/PSIITB)**(PROFR2-1.D0)
     &           /PSIITB
         ENDIF
         DQPSI=DQPSI
      ELSE
         CALL SPL1DD(PSIN,QPSI,DQPSI,PSITRX,UQPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX DQPSI: SPL1DD : IERR=',IERR
         DQPSI=DQPSI
      ENDIF
      RETURN
      END
C
C   ************************************************ 
C   **                  J0(psi)                   **
C   ************************************************
C
      REAL*8 FUNCTION HJPSID(PSIN1)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MDLEQF.LT.5) THEN
         HJPSID=-PJ0*(1.D0-PSIN**PROFR0)**(PROFJ0+1.D0)
     &              /(PROFR0*(PROFJ0+1.D0))
     &          -PJ1*(1.D0-PSIN**PROFR1)**(PROFJ1+1.D0)
     &              /(PROFR1*(PROFJ1+1.D0))
     &          -PJ2*(1.D0-PSIN**PROFR2)**(PROFJ2+1.D0)
     &              /(PROFR2*(PROFJ2+1.D0))
         HJPSID=HJPSID*1.D6
      ELSE
         CALL SPL1DI(PSIN,HJPSID,PSITRX,UJPSI,UJPSI0,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX HJPSID: SPL1DF : IERR=',IERR
         HJPSID=HJPSID*1.D6
      ENDIF
      RETURN
      END
C
      REAL*8 FUNCTION HJPSI(PSIN1)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MDLEQF.LT.5) THEN
         HJPSI=PJ0*(1.D0-PSIN**PROFR0)**PROFJ0
     &                   *PSIN**(PROFR0-1.D0)
     &        +PJ1*(1.D0-PSIN**PROFR1)**PROFJ1
     &                   *PSIN**(PROFR1-1.D0)
     &        +PJ2*(1.D0-PSIN**PROFR2)**PROFJ2
     &                   *PSIN**(PROFR2-1.D0)
         HJPSI=HJPSI*1.D6
      ELSE
         CALL SPL1DF(PSIN,HJPSI,PSITRX,UJPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX HJPSI: SPL1DF : IERR=',IERR
         HJPSI=HJPSI*1.D6
      ENDIF
      RETURN
      END
C
C   ************************************************ 
C   **                  T(psi)                    **
C   ************************************************
C
      REAL*8 FUNCTION TPSI(PSIN1)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MDLEQF.LT.5) THEN
         TPSI=PTS+(PT0-PTS)*(1.D0-PSIN**PROFR0)**PROFT0
     &                 +PT1*(1.D0-PSIN**PROFR1)**PROFT1
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            TPSI=TPSI
     &          +PT2*(1.D0-(PSIN/PSIITB)**PROFR2)**PROFT2
         ENDIF
         TPSI=TPSI*1.D3*AEE
      ELSE
         CALL SPL1DF(PSIN,TPSI,PSITRX,UTPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX TPSI: SPL1DF : IERR=',IERR
         TPSI=TPSI*1.D3*AEE
      ENDIF
      RETURN
      END
C
      REAL*8 FUNCTION DTPSI(PSIN1)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MDLEQF.LT.5) THEN
         DTPSI=(PT0-PTS)
     &            *PROFT0*(1.D0-PSIN**PROFR0)**(PROFT0-1.D0)
     &            *PROFR0*PSIN**(PROFR0-1.D0)
     &        +PT1*PROFT1*(1.D0-PSIN**PROFR1)**(PROFT1-1.D0)
     &            *PROFR1*PSIN**(PROFR1-1.D0)
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            DTPSI=DTPSI
     &           +PT2*PROFT2
     &        *(1.D0-(PSIN/PSIITB)**PROFR2)**(PROFT2-1.D0)
     &        *PROFR2*(PSIN/PSIITB)**(PROFR2-1.D0)
     &         /PSIITB
         ENDIF
         DTPSI=DTPSI*1.D3*AEE
      ELSEIF(MDLEQF.EQ.2) THEN
         CALL SPL1DD(PSIN,TPSI,DTPSI,PSITRX,UTPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX DTPSI: SPL1DD : IERR=',IERR
         DTPSI=DTPSI*1.D3*AEE
      ENDIF
      RETURN
      END
C
C   ************************************************ 
C   **               OMGPSI(psi)                  **
C   ************************************************
C
      REAL*8 FUNCTION OMGPSI(PSIN1)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MDLEQF.LT.5) THEN
         VPSI=PV0*(1.D0-PSIN**PROFR0)**PROFV0
     &       +PV1*(1.D0-PSIN**PROFR1)**PROFV1
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            VPSI=VPSI
     &         +PV2*(1.D0-(PSIN/PSIITB)**PROFR2)**PROFV2
         ENDIF
         OMGPSI=VPSI/RAXIS
      ELSE
         CALL SPL1DF(PSIN,VPSI,PSITRX,UVTPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX OMGPSI: SPL1DF : IERR=',IERR
         OMGPSI=VPSI/RAXIS
      ENDIF
      RETURN
      END
C
      REAL*8 FUNCTION DOMGPSI(PSIN1)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      IF(PSIN1.LE.0.D0) THEN
         PSIN=0.D0
      ELSEIF(PSIN1.GE.1.D0) THEN
         PSIN=1.D0
      ELSE
         PSIN=PSIN1
      ENDIF
C
      IF(MDLEQF.LT.5) THEN
         DVPSI=-PV0*PROFV0*(1.D0-PSIN**PROFR0)**(PROFV0-1.D0)
     &             *PROFR0*PSIN**(PROFR0-1.D0)
     &         -PV1*PROFV1*(1.D0-PSIN**PROFR1)**(PROFV1-1.D0)
     &                            *PROFR1*PSIN**(PROFR1-1.D0)
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            DVPSI=DVPSI
     &        -PV2*PROFV2
     &        *(1.D0-(PSIN/PSIITB)**PROFR2)**(PROFV2-1.D0)
     &       *PROFR2*(PSIN/PSIITB)**(PROFR2-1.D0)
     &        /PSIITB
         ENDIF
         DOMGPSI=DVPSI/RAXIS
      ELSE
         CALL SPL1DD(PSIN,VPSI,DVPSI,PSITRX,UVTPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX DOMGPSI: SPL1DD : IERR=',IERR
         DOMGPSI=DVPSI/RAXIS
      ENDIF
      RETURN
      END
