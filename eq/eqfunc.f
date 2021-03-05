C     $Id$
C
C   *************************************
C   *****  Argument Psi conversion  *****
C   *************************************
C
C     ---- PSIN is converted to PSIPN or PSITN ----
C
      SUBROUTINE EQCNVA(PSIPNL,PSINL)
C
      INCLUDE 'eqcomc.inc'
C
      IF(PSIPNL.LE.0.D0) THEN
         PSINL=0.D0
      ELSEIF(PSIPNL.GE.1.D0) THEN
         PSINL=1.D0
      ELSE
         IF(MDLEQA.EQ.0) THEN
            PSINL=PSIPNL
         ELSE
            PSINL=EQPSITN(PSIPNL)
            IF(PSINL.LE.0.D0) PSINL=0.D0
            IF(PSINL.GE.1.D0) PSINL=1.D0
         ENDIF
      ENDIF
      RETURN
      END
C
C   ****************************************
C   *****  Calculate Profile Function  *****
C   ****************************************
C
      SUBROUTINE EQFUNC(PSINL,F,DF,F0,FS,F1,F2,PSIITB,
     &                  PROFR0,PROFR1,PROFR2,PROFF0,PROFF1,PROFF2)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      ARG0=FPOW(PSINL,PROFR0)
      ARG1=FPOW(PSINL,PROFR1)
      F = FS
     &  + (F0-FS)       *FPOW(1.D0-ARG0,PROFF0)
     &  + F1            *FPOW(1.D0-ARG1,PROFF1)
      DF=-(F0-FS)*PROFF0*FPOW(1.D0-ARG0,PROFF0-1.D0)
     &           *PROFR0*FPOW(     PSINL,PROFR0-1.D0)
     &   -F1     *PROFF1*FPOW(1.D0-ARG1,PROFF1-1.D0)
     &           *PROFR1*FPOW(     PSINL,PROFR1-1.D0)
      IF(PSINL.LT.PSIITB) THEN
         ARG2=FPOW(PSINL/PSIITB,PROFR2)
         F = F +F2       *FPOW(1.D0-ARG2,  PROFF2)
         DF= DF-F2*PROFF2*FPOW(1.D0-ARG2,  PROFF2-1.D0)
     &            *PROFR2*FPOW(PSINL/PSIITB,PROFR2-1.D0)
     &            /PSIITB
      ENDIF
      RETURN
      END
C
C   **************************************************
C   *****  Calculate factor for Psip derivative  *****
C   **************************************************
C
      SUBROUTINE EQFDPP(PSIPNL,FDN)
C
      INCLUDE 'eqcomc.inc'
C
      IF(MDLEQA.EQ.0) THEN
         FDN=1.D0/PSIPA
      ELSE
         QPVL=EQQPV(PSIPNL)
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
      SUBROUTINE EQPPSI(PSIPNL,PPSI,DPPSI)
C
      USE libspl1d
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      MDLEQFL=MOD(MDLEQF,10)
      IF(MDLEQFL.LT.5) THEN
         CALL EQCNVA(PSIPNL,PSINL)
         CALL EQFUNC(PSINL,F,DF,PP0,0.D0,PP1,PP2,PSIITB,
     &               PROFR0,PROFR1,PROFR2,PROFP0,PROFP1,PROFP2)
         CALL EQFDPP(PSIPNL,FDN)
      ELSE
         IF(PSIPNL.LT.0.D0) PSIPNL=0.D0
         IF(PSIPNL.GT.1.D0) PSIPNL=1.D0
         PSITNL=EQPSITN(PSIPNL)
         CALL SPL1DD(PSITNL,F,DF,PSITRX,UPPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQPPSI: SPL1DD : IERR=',IERR
         IF(F.LT.0.D0) F=0.D0
         QPVL=EQQPV(PSIPNL)
         FDN=QPVL/PSITA
      ENDIF
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
      SUBROUTINE EQFPSI(PSIPNL,FPSI,DFPSI)
C
      USE libspl1d
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      MDLEQFL=MOD(MDLEQF,10)
      IF(MDLEQFL.LT.5) THEN
         CALL EQCNVA(PSIPNL,PSINL)
         FS=2.D0*PI*BB*RR
         CALL EQFUNC(PSINL,F,DF,FF0+FS,FS,FF1,FF2,PSIITB,
     &               PROFR0,PROFR1,PROFR2,PROFF0,PROFF1,PROFF2)
         CALL EQFDPP(PSIPNL,FDN)
      ELSE
         PSITNL=EQPSITN(PSIPNL)
         CALL SPL1DD(PSITNL,F,DF,PSITRX,UFPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQFPSI: SPL1DD : IERR=',IERR
         QPVL=EQQPV(PSIPNL)
         FDN=QPVL/PSITA
      ENDIF
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
      SUBROUTINE EQQPSI(PSIPNL,QPSI)
C
      USE libspl1d
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      MDLEQFL=MOD(MDLEQF,10)
      IF(MDLEQFL.LT.5) THEN
         PSITNL=EQPSITN(PSIPNL)
         QPSI=Q0+(QA-Q0)*PSITNL
      ELSE
         PSITNL=EQPSITN(PSIPNL)
         CALL SPL1DD(PSITNL,F,DF,PSITRX,UQPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQQPSI: SPL1DD : IERR=',IERR
         QPSI  =      F
      ENDIF
      RETURN
      END
C
C   ************************************************ 
C   ** Plasma current density                     **
C   **                  J(psin)                   **
C   ************************************************
C
      SUBROUTINE EQJPSI(PSIPNL,HJPSID,HJPSI)
C
      USE libspl1d
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      MDLEQFL=MOD(MDLEQF,10)
      IF(MDLEQFL.LT.5) THEN
         CALL EQCNVA(PSIPNL,PSINL)
         ARG0=FPOW(PSINL,PROFR0)
         ARG1=FPOW(PSINL,PROFR1)
         ARG2=FPOW(PSINL,PROFR2)
         FD=-PJ0*FPOW(1.D0-ARG0,PROFJ0+1.D0)
     &                /(PROFR0*(PROFJ0+1.D0))
     &      -PJ1*FPOW(1.D0-ARG1,PROFJ1+1.D0)
     &                /(PROFR1*(PROFJ1+1.D0))
     &      -PJ2*FPOW(1.D0-ARG2,PROFJ2+1.D0)
     &                /(PROFR2*(PROFJ2+1.D0))
         F = PJ0*FPOW(1.D0-ARG0,PROFJ0)
     &          *FPOW(PSINL,PROFR0-1.D0)
     &      +PJ1*FPOW(1.D0-ARG1,PROFJ1)
     &          *FPOW(PSINL,PROFR1-1.D0)
     &      +PJ2*FPOW(1.D0-ARG2,PROFJ2)
     &          *FPOW(PSINL,PROFR2-1.D0)
         FD=FD*1.D6
         F =F *1.D6
         CALL EQFDPP(PSIPNL,FDN)
      ELSE
         PSITNL=EQPSITN(PSIPNL)
         CALL SPL1DI(PSITNL,FD,PSITRX,UJPSI,UJPSI0,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQJPSI: SPL1DI : IERR=',IERR
         CALL SPL1DF(PSITNL,F,PSITRX,UJPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQJPSI: SPL1DF : IERR=',IERR
         QPVL=EQQPV(PSIPNL)
         FDN=QPVL/PSITA
      ENDIF
C
      HJPSID= FD/FDN
      HJPSI =  F
C      WRITE(6,'(1P4E12.4)') PSIPNL,PSINL,HJPSI,HJPSID
      RETURN
      END
C
C   ************************************************ 
C   ** Plasma tempertature                        **
C   **                  T(psin)                   **
C   ************************************************
C
      SUBROUTINE EQTPSI(PSIPNL,TPSI,DTPSI)
C
      USE libspl1d
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      MDLEQFL=MOD(MDLEQF,10)
      IF(MDLEQFL.LT.5) THEN
         CALL EQCNVA(PSIPNL,PSINL)
         CALL EQFUNC(PSINL,F,DF,PT0,PTSEQ,PT1,PT2,PSIITB,
     &               PROFR0,PROFR1,PROFR2,PROFTP0,PROFTP1,PROFTP2)
         CALL EQFDPP(PSIPNL,FDN)
      ELSE
         PSITNL=EQPSITN(PSIPNL)
         CALL SPL1DD(PSITNL,F,DF,PSITRX,UTPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQTPSI: SPL1DD : IERR=',IERR
         QPVL=EQQPV(PSIPNL)
         FDN=QPVL/PSITA
      ENDIF
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
      SUBROUTINE EQOPSI(PSIPNL,OMGPSI,DOMGPSI)
C
      USE libspl1d
      INCLUDE 'eqcomc.inc'
      INCLUDE 'eqcom4.inc'
C
      MDLEQFL=MOD(MDLEQF,10)
      IF(MDLEQFL.LT.5) THEN
         CALL EQCNVA(PSIPNL,PSINL)
         CALL EQFUNC(PSINL,F,DF,PV0,0.D0,PV1,PV2,PSIITB,
     &               PROFR0,PROFR1,PROFR2,PROFV0,PROFV1,PROFV2)
         CALL EQFDPP(PSIPNL,FDN)
      ELSE
         PSITNL=EQPSITN(PSIPNL)
         CALL SPL1DD(PSITNL,F,DF,PSITRX,UVTPSI,NTRMAX,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX EQOPSI: SPL1DD : IERR=',IERR
         QPVL=EQQPV(PSIPNL)
         FDN=QPVL/PSITA
      ENDIF
C
      OMGPSI=      F/RAXIS
      DOMGPSI=FDN*DF/RAXIS
      RETURN
      END
