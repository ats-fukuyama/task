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
C   ** Plasma pressure                            **
C   **                  P(psin)                   **
C   ************************************************
C
      SUBROUTINE EQPPSI(PSIN1,PPSI,DPPSI)
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
         ARG0=FPOW(PSIN,PROFR0)
         ARG1=FPOW(PSIN,PROFR1)
         PPSI=  PP0       *FPOW(1.D0-ARG0,PROFP0)
     &         +PP1       *FPOW(1.D0-ARG1,PROFP1)
         DPPSI=-PP0*PROFP0*FPOW(1.D0-ARG0,PROFP0-1.D0)
     &             *PROFR0*FPOW(     PSIN,PROFR0-1.D0)
     &         -PP1*PROFP1*FPOW(1.D0-ARG1,PROFP1-1.D0)
     &             *PROFR1*FPOW(     PSIN,PROFR1-1.D0)
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            ARG2=FPOW(PSIN/PSIITB,PROFR2)
            PPSI= PPSI
     &           +PP2       *FPOW(1.D0-ARG2,  PROFP2)
            DPPSI=DPPSI
     &           -PP2*PROFP2*FPOW(1.D0-ARG2,  PROFP2-1.D0)
     &               *PROFR2*FPOW(PSIN/PSIITB,PROFR2-1.D0)
     &           /PSIITB
         ENDIF
         PPSI=PPSI*1.D6
         DPPSI=DPPSI*1.D6
      ELSE
         CALL SPL1DD(PSIN,PPSI,DPPSI,PSITRX,UPPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQPPSI: SPL1DD : IERR=',IERR
         IF(PPSI.LT.0.D0) PPSI=0.D0
      ENDIF
      PPSI=PPSI*1.D6
      DPPSI=DPPSI*1.D6
      RETURN
      END
C
C   ************************************************ 
C   ** Poloidal currenet                          **
C   **                  F(psin)=B*R               **
C   ************************************************
C
      SUBROUTINE EQFPSI(PSIN1,FPSI,DFPSI)
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
         ARG0=FPOW(PSIN,PROFR0)
         ARG1=FPOW(PSIN,PROFR1)
         FPSI=  BB*RR
     &         +FF0       *FPOW(1.D0-ARG0,PROFF0)
     &         +FF1       *FPOW(1.D0-ARG1,PROFF1)
         DFPSI=-FF0*PROFF0*FPOW(1.D0-ARG0,PROFF0-1.D0)
     &             *PROFR0*FPOW(     PSIN,PROFR0-1.D0)
     &         -FF1*PROFF1*FPOW(1.D0-ARG1,PROFF1-1.D0)
     &             *PROFR1*FPOW(     PSIN,PROFR1-1.D0)
C
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            ARG2=FPOW(PSIN/PSIITB,PROFR2)
            FPSI= FPSI
     &           +FF2       *FPOW(1.D0-ARG2,  PROFF2)
            DFPSI=DFPSI
     &           -FF2*PROFF2*FPOW(1.D0-ARG2,  PROFF2-1.D0)
     &               *PROFR2*FPOW(PSIN/PSIITB,PROFR2-1.D0)
     &               /PSIITB
         ENDIF
      ELSE
         CALL SPL1DD(PSIN,FPSI,DFPSI,PSITRX,UFPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQFPSI: SPL1DD : IERR=',IERR
      ENDIF
      RETURN
      END
C
C   ************************************************ 
C   ** Safety factor                              **
C   **                 Q(psin)                    **
C   ************************************************
C
      SUBROUTINE EQQPSI(PSIN1,QPSI,DQPSI)
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
         ARG0=FPOW(PSIN,PROFR0)
         ARG1=FPOW(PSIN,PROFR1)
         QPSI=   QQS
     &         +(QQ0-QQS)       *FPOW(1.D0-ARG0,PROFQ0)
     &         + QQ1            *FPOW(1.D0-ARG1,PROFQ1)
         DQPSI=-(QQ0-QQS)*PROFQ0*FPOW(1.D0-ARG0,PROFQ0-1.D0)
     &                   *PROFR0*FPOW(     PSIN,PROFR0-1.D0)
     &         -QQ1      *PROFQ1*FPOW(1.D0-ARG1,PROFQ1-1.D0)
     &                   *PROFR1*FPOW(     PSIN,PROFR1-1.D0)
C
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            ARG2=FPOW(PSIN/PSIITB,PROFR2)
            QPSI= QPSI
     &           +QQ2       *FPOW(1.D0-ARG2,  PROFQ2)
            DQPSI=DQPSI
     &           +QQ2*PROFQ2*FPOW(1.D0-ARG2,  PROFQ2-1.D0)
     &               *PROFR2*FPOW(PSIN/PSIITB,PROFR2-1.D0)
     &               /PSIITB
         ENDIF
      ELSE
         CALL SPL1DD(PSIN,QPSI,DQSI,QSITRX,UQPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQQPSI: SPL1DD : IERR=',IERR
      ENDIF
      RETURN
      END
C
C   ************************************************ 
C   ** Plasma current density                     **
C   **                  J(psin)                   **
C   ************************************************
C
      SUBROUTINE EQJPSI(PSIN1,HJPSID,HJPSI)
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
         ARG0=FPOW(PSIN,PROFR0)
         ARG1=FPOW(PSIN,PROFR1)
         ARG2=FPOW(PSIN,PROFR2)
         HJPSID=-PJ0*FPOW(1.D0-ARG0,PROFJ0+1.D0)
     &                    /(PROFR0*(PROFJ0+1.D0))
     &          -PJ1*FPOW(1.D0-ARG1,PROFJ1+1.D0)
     &                    /(PROFR1*(PROFJ1+1.D0))
     &          -PJ2*FPOW(1.D0-ARG2,PROFJ2+1.D0)
     &                    /(PROFR2*(PROFJ2+1.D0))
         HJPSI=  PJ0*FPOW(1.D0-ARG0,PROFJ0)
     &              *FPOW(PSIN,PROFR0-1.D0)
     &          +PJ1*FPOW(1.D0-ARG1,PROFJ1)
     &              *FPOW(PSIN,PROFR1-1.D0)
     &          +PJ2*FPOW(1.D0-ARG2,PROFJ2)
     &              *FPOW(PSIN,PROFR2-1.D0)
      ELSE
         CALL SPL1DI(PSIN,HJPSID,PSITRX,UJPSI,UJPSI0,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQJPSI: SPL1DI : IERR=',IERR
         CALL SPL1DF(PSIN,HJPSI,PSITRX,UJPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQJPSI: SPL1DF : IERR=',IERR
      ENDIF
      HJPSID=HJPSID*1.D6
      HJPSI=HJPSI*1.D6
      RETURN
      END
C
C   ************************************************ 
C   ** Plasma tempertature                        **
C   **                  T(psin)                   **
C   ************************************************
C
      SUBROUTINE EQTPSI(PSIN1,TPSI,DTPSI)
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
         ARG0=FPOW(PSIN,PROFR0)
         ARG1=FPOW(PSIN,PROFR1)
         TPSI= PTS
     &        +(PT0-PTS)       *FPOW(1.D0-ARG0,PROFT0)
     &             +PT1        *FPOW(1.D0-ARG1,PROFT1)
         DTPSI=(PT0-PTS)*PROFT0*FPOW(1.D0-ARG0,PROFT0-1.D0)
     &                  *PROFR0*FPOW(     PSIN,PROFR0-1.D0)
     &        +PT1      *PROFT1*FPOW(1.D0-ARG1,PROFT1-1.D0)
     &                  *PROFR1*FPOW(     PSIN,PROFR1-1.D0)
C
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            ARG2=FPOW(PSIN/PSIITB,PROFR2)
            TPSI= TPSI
     &           +PT2       *FPOW(1.D0-ARG2,  PROFT2)
            DTPSI=DTPSI
     &           +PT2*PROFT2*FPOW(1.D0-ARG2,  PROFT2-1.D0)
     &               *PROFR2*FPOW(PSIN/PSIITB,PROFR2-1.D0)
     &               /PSIITB
         ENDIF
      ELSE
         CALL SPL1DD(PSIN,TPSI,DTPSI,PSITRX,UTPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQTPSI: SPL1DD : IERR=',IERR
      ENDIF
      TPSI=  TPSI*1.D3*AEE
      DTPSI=DTPSI*1.D3*AEE
      RETURN
      END
C
C   ************************************************ 
C   ** Toroidal rotation angular velocisty        **
C   **               OPSI(psin)                   **
C   ************************************************
C
      SUBROUTINE EQOPSI(PSIN1,OMGPSI,DOMGPSI)
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
         ARG0=FPOW(PSIN,PROFR0)
         ARG1=FPOW(PSIN,PROFR1)
         VPSI=  PV0       *FPOW(1.D0-ARG0,PROFV0)
     &         +PV1       *FPOW(1.D0-ARG1,PROFV1)
         DVPSI=-PV0*PROFV0*FPOW(1.D0-ARG0,PROFV0-1.D0)
     &             *PROFR0*FPOW(     PSIN,PROFR0-1.D0)
     &         -PV1*PROFV1*FPOW(1.D0-ARG1,PROFV1-1.D0)
     &             *PROFR1*FPOW(     PSIN,PROFR1-1.D0)
C
         PSIITB=RHOITB**2
         IF(PSIN.LT.PSIITB) THEN
            ARG2=FPOW(PSIN/PSIITB,PROFR2)
            VPSI= VPSI
     &           +PV2       *FPOW(1.D0-ARG2,PROFV2)
            DVPSI=DVPSI
     &           -PV2*PROFV2*FPOW(1.D0-ARG2,PROFV2-1.D0)
     &               *PROFR2*FPOW(PSIN/PSIITB,PROFR2-1.D0)
     &               /PSIITB
         ENDIF
         OMGPSI=VPSI/RAXIS
      ELSE
         CALL SPL1DD(PSIN,VPSI,DVPSI,PSITRX,UVTPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQOPSI: SPL1DD : IERR=',IERR
      ENDIF
      OMGPSI= VPSI /RAXIS
      DOMGPSI=DVPSI/RAXIS
      RETURN
      END
C
C   ***************************************
C   **         Power function            **
C   ***************************************
C
      FUNCTION FPOW(X,Y)
C
      REAL*8 X,Y,FPOW
C
      IF(X.EQ.0.D0) THEN
         IF(Y.EQ.0.D0) THEN
            FPOW=1.D0
         ELSE
            FPOW=0.D0
         ENDIF
      ELSE
         FPOW=X**Y
      ENDIF
      RETURN
      END
