C     $Id$
C
C   *************************************
C   *****  Argument Psi conversion  *****
C   *************************************
C
C     ---- PSIN is converted to PSIPN or PSITN ----
C
      SUBROUTINE EQCNVA(PSIN,PSIL)
C
      INCLUDE 'eqcomc.inc'
C
      IF(PSIN.LE.0.D0) THEN
         PSIL=0.D0
      ELSEIF(PSIN.GE.1.D0) THEN
         PSIL=1.D0
      ELSE
         IF(MDLEQA.EQ.0) THEN
            PSIL=PSIN
         ELSE
            PSIL=EQPSITN(PSIN)
            IF(PSIL.LE.0.D0) PSIL=0.D0
            IF(PSIL.GE.1.D0) PSIL=1.D0
         ENDIF
      ENDIF
      RETURN
      END
C
C   ****************************************
C   *****  Calculate Profile Function  *****
C   ****************************************
C
      SUBROUTINE EQFUNC(PSIL,F,DF,F0,FS,F1,F2,PSIITB,
     &                  PROFR0,PROFR1,PROFR2,PROFF0,PROFF1,PROFF2)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      ARG0=FPOW(PSIL,PROFR0)
      ARG1=FPOW(PSIL,PROFR1)
      F = FS
     &  + (F0-FS)       *FPOW(1.D0-ARG0,PROFF0)
     &  + F1            *FPOW(1.D0-ARG1,PROFF1)
      DF=-(F0-FS)*PROFF0*FPOW(1.D0-ARG0,PROFF0-1.D0)
     &           *PROFR0*FPOW(     PSIL,PROFR0-1.D0)
     &   -F1     *PROFF1*FPOW(1.D0-ARG1,PROFF1-1.D0)
     &           *PROFR1*FPOW(     PSIL,PROFR1-1.D0)
      IF(PSIL.LT.PSIITB) THEN
         ARG2=FPOW(PSIL/PSIITB,PROFR2)
         F = F +F2       *FPOW(1.D0-ARG2,  PROFF2)
         DF= DF-F2*PROFF2*FPOW(1.D0-ARG2,  PROFF2-1.D0)
     &            *PROFR2*FPOW(PSIL/PSIITB,PROFR2-1.D0)
     &            /PSIITB
      ENDIF
      RETURN
      END
C
C   **************************************************
C   *****  Calculate factor for Psip derivative  *****
C   **************************************************
C
      SUBROUTINE EQFDPP(PSIN,FDN)
C
      INCLUDE 'eqcomc.inc'
C
      IF(MDLEQA.EQ.0) THEN
         FDN=-1.D0/PSI0
      ELSE
         QPVL=EQQPV(PSIN)
         FDN=QPVL/PSITA
      ENDIF
      RETURN
      END
C
C   ************************************************ 
C   ** Plasma pressure                            **
C   **                  P(psin)                   **
C   ************************************************
C
      SUBROUTINE EQPPSI(PSIN,PPSI,DPPSI)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      CALL EQCNVA(PSIN,PSIL)
C
      IF(MDLEQF.LT.5) THEN
         CALL EQFUNC(PSIL,F,DF,PP0,0.D0,PP1,PP2,PSIITB,
     &               PROFR0,PROFR1,PROFR2,PROFP0,PROFP1,PROFP2)
      ELSE
         CALL SPL1DD(PSIL,F,DF,PSITRX,UPPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQPPSI: SPL1DD : IERR=',IERR
         IF(F.LT.0.D0) F=0.D0
      ENDIF
C
      CALL EQFDPP(PSIN,FDN)
C
      PPSI  =      F*1.D6
      DPPSI = FDN*DF*1.D6
      RETURN
      END
C
C   ************************************************ 
C   ** Poloidal currenet                          **
C   **                  F(psin)=2*PI*B*R          **
C   ************************************************
C
      SUBROUTINE EQFPSI(PSIN,FPSI,DFPSI)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      CALL EQCNVA(PSIN,PSIL)
C
      IF(MDLEQF.LT.5) THEN
         FS=2.D0*PI*BB*RR
         CALL EQFUNC(PSIL,F,DF,FF0+FS,FS,FF1,FF2,PSIITB,
     &               PROFR0,PROFR1,PROFR2,PROFF0,PROFF1,PROFF2)
      ELSE
         CALL SPL1DD(PSIL,F,DF,PSITRX,UFPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQFPSI: SPL1DD : IERR=',IERR
      ENDIF
C
      CALL EQFDPP(PSIN,FDN)
C
      FPSI  =      F
      DFPSI = FDN*DF
      RETURN
      END
C
C   ************************************************ 
C   ** Safety factor                              **
C   **                 Q(psin)                    **
C   ************************************************
C
      SUBROUTINE EQQPSI(PSIN,QPSI,DQPSI)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      CALL EQCNVA(PSIN,PSIL)
C
      IF(MDLEQF.LT.5) THEN
         CALL EQFUNC(PSIL,F,DF,QQ0,QQS,QQ1,QQ2,PSIITB,
     &               PROFR0,PROFR1,PROFR2,PROFQ0,PROFQ1,PROFQ2)
      ELSE
         CALL SPL1DD(PSIL,F,DF,QSITRX,UQPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQQPSI: SPL1DD : IERR=',IERR
      ENDIF
C
      CALL EQFDPP(PSIN,FDN)
C
      QPSI  =      F
      DQPSI = FDN*DF
      RETURN
      END
C
C   ************************************************ 
C   ** Plasma current density                     **
C   **                  J(psin)                   **
C   ************************************************
C
      SUBROUTINE EQJPSI(PSIN,HJPSID,HJPSI)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      CALL EQCNVA(PSIN,PSIL)
C
      IF(MDLEQF.LT.5) THEN
         ARG0=FPOW(PSIL,PROFR0)
         ARG1=FPOW(PSIL,PROFR1)
         ARG2=FPOW(PSIL,PROFR2)
         FD=-PJ0*FPOW(1.D0-ARG0,PROFJ0+1.D0)
     &                /(PROFR0*(PROFJ0+1.D0))
     &      -PJ1*FPOW(1.D0-ARG1,PROFJ1+1.D0)
     &                /(PROFR1*(PROFJ1+1.D0))
     &      -PJ2*FPOW(1.D0-ARG2,PROFJ2+1.D0)
     &                /(PROFR2*(PROFJ2+1.D0))
         F = PJ0*FPOW(1.D0-ARG0,PROFJ0)
     &          *FPOW(PSIL,PROFR0-1.D0)
     &      +PJ1*FPOW(1.D0-ARG1,PROFJ1)
     &          *FPOW(PSIL,PROFR1-1.D0)
     &      +PJ2*FPOW(1.D0-ARG2,PROFJ2)
     &          *FPOW(PSIL,PROFR2-1.D0)
      ELSE
         CALL SPL1DI(PSIL,FD,PSITRX,UJPSI,UJPSI0,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQJPSI: SPL1DI : IERR=',IERR
         CALL SPL1DF(PSIL,F,PSITRX,UJPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQJPSI: SPL1DF : IERR=',IERR
      ENDIF
C
      CALL EQFDPP(PSIN,FDN)
C
      HJPSID= FD*1.D6/FDN
      HJPSI =  F*1.D6
      RETURN
      END
C
C   ************************************************ 
C   ** Plasma tempertature                        **
C   **                  T(psin)                   **
C   ************************************************
C
      SUBROUTINE EQTPSI(PSIN,TPSI,DTPSI)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      CALL EQCNVA(PSIN,PSIL)
C
      IF(MDLEQF.LT.5) THEN
         CALL EQFUNC(PSIL,F,DF,PT0,PTS,PT1,PT2,PSIITB,
     &               PROFR0,PROFR1,PROFR2,PROFT0,PROFT1,PROFT2)
      ELSE
         CALL SPL1DD(PSIN,F,DF,PSITRX,UTPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQTPSI: SPL1DD : IERR=',IERR
      ENDIF
C
      CALL EQFDPP(PSIN,FDN)
C
      TPSI=      F*1.D3*AEE
      DTPSI=FDN*DF*1.D3*AEE
      RETURN
      END
C
C   ************************************************ 
C   ** Toroidal rotation angular velocisty        **
C   **               OPSI(psin)                   **
C   ************************************************
C
      SUBROUTINE EQOPSI(PSIN,OMGPSI,DOMGPSI)
C
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      CALL EQCNVA(PSIN,PSIL)
C
      IF(MDLEQF.LT.5) THEN
         CALL EQFUNC(PSIL,F,DF,PV0,0.D0,PV1,PV2,PSIITB,
     &               PROFR0,PROFR1,PROFR2,PROFV0,PROFV1,PROFV2)
      ELSE
         CALL SPL1DD(PSIN,F,DF,PSITRX,UVTPSI,NTRMAX+2,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQOPSI: SPL1DD : IERR=',IERR
      ENDIF
C
      CALL EQFDPP(PSIN,FDN)
C
      OMGPSI=  F/RAXIS
      DOMGPSI=DF/RAXIS
      RETURN
      END
