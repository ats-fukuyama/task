C     $Id$
C
C ************************************************************
C
C                      CALCULATION OF D AND F
C
C ************************************************************
C
      SUBROUTINE FPCOEF
C
      INCLUDE 'fpcomm.inc'
C
C     ----- Parallel electric field accleration term -----
C
      CALL FPCALE
C
C     ----- Quasi-linear wave-particle interaction term -----
C
      IF(MODELW.EQ.0) THEN
         CALL FPCALW
      ELSEIF(MODELW.EQ.1) THEN
         CALL FPCALWR
      ELSEIF(MODELW.EQ.2) THEN
         CALL FPCALWR
      ELSEIF(MODELW.EQ.3) THEN
         CALL FPCALWM
      ELSEIF(MODELW.EQ.4) THEN
         CALL FPCALWM
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELW =',MODELW
      ENDIF
C
C     ----- Collisional slowing down and diffusion term -----
C
      CALL FPCALC
C
C     ----- Sum up velocity diffusion terms -----
C
      DO NR=1,NRMAX
      DO NP=1,NPMAX+1
      DO NTH=1,NTHMAX
         DPP(NTH,NP,NR)=DCPP(NTH,NP,NR)+DWPP(NTH,NP,NR)
         DPT(NTH,NP,NR)=DCPT(NTH,NP,NR)+DWPT(NTH,NP,NR)
         FPP(NTH,NP,NR)=FEPP(NTH,NP,NR)+FCPP(NTH,NP,NR)
      ENDDO
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX+1
         DTP(NTH,NP,NR)=DCTP(NTH,NP,NR)+DWTP(NTH,NP,NR)
         DTT(NTH,NP,NR)=DCTT(NTH,NP,NR)+DWTT(NTH,NP,NR)
         FTH(NTH,NP,NR)=FETH(NTH,NP,NR)+FCTH(NTH,NP,NR)
      ENDDO
      ENDDO
      ENDDO
C
C     ----- Radial diffusion term (disabled) -----
C
C      CALL FPCALR
C
C     ----- Particle source term (disabled) -----
C
      DO NR=1,NRMAX
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX
         SP(NTH,NP,NR)=0.D0
      ENDDO
      ENDDO
      ENDDO
C
C     ****************************
C     Boundary condition at p=pmax
C     ****************************
C
      DO NR=1,NRMAX
      DO NTH=1,NTHMAX
C         DPP(NTH,NPMAX,NR)=0.5D0*DPP(NTH,NPMAX,NR)
C         DPT(NTH,NPMAX,NR)=0.5D0*DPT(NTH,NPMAX,NR)
         DPP(NTH,NPMAX+1,NR)=0.D0
         DPT(NTH,NPMAX+1,NR)=0.D0
         FPMAX=FPP(NTH,NPMAX+1,NR)
         FPP(NTH,NPMAX+1,NR)=MAX(FPMAX,0.D0)
      ENDDO
      ENDDO
C
C      DO NR=1,NRMAX
C      DO NTH=1,NTHMAX+1
C         DTP(NTH,NPMAX,NR)=0.D0
C         DTT(NTH,NPMAX,NR)=0.D0
C         FTH(NTH,NPMAX,NR)=0.D0
C       ENDDO
C       ENDDO
C
      DO NR=1,NRMAX
      DO NP=1,NPMAX+1
      DO NTH=1,NTHMAX
         WEIGHP(NTH,NP,NR)=FPWEGH(-DELP*FPP(NTH,NP,NR),
     &                             DPP(NTH,NP,NR))
      ENDDO
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX+1
         WEIGHT(NTH,NP,NR)=FPWEGH(-DELTH*PM(NP)*FTH(NTH,NP,NR),
     &                             DTT(NTH,NP,NR))
      ENDDO
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX+1
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX
         WEIGHR(NTH,NP,NR)=FPWEGH(-DELR*FRR(NTH,NP,NR),
     &                             DRR(NTH,NP,NR))
      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
C
C ************************************************************
C
C                 CALCULATION OF DRR AND VR
C
C ************************************************************
C
      SUBROUTINE FPCALR
C
      INCLUDE 'fpcomm.inc'
C
      DO 10 NR=1,NRMAX+1
         RHON=RG(NR)
         CALL PLPROF(RHON)
         RTFPL=RTPR(NSFP)/RTFP0
         FACTR=-2.D0*RG(NR)/RA
      DO 10 NP=1,NPMAX
          FACTP=1.D0/SQRT(1.D0+PG(NP)**2/RTFPL)
      DO 10 NTH=1,NTHMAX
         DRR(NTH,NP,NR)=DRR0*FACTP
         FRR(NTH,NP,NR)=DRR0*FACTP*FACTR
   10 CONTINUE
C
      IF (MODELA.EQ.1) THEN
         DO 100 NR=2,NRMAX
         DO 100 NTH=1,NTHMAX
            RL=0.5D0*(RLAMDA(NTH,NR-1)+RLAMDA(NTH,NR))
         DO 100 NP=1,NPMAX
            DRR(NTH,NP,NR)=RL*DRR(NTH,NP,NR)
            FRR(NTH,NP,NR)=RL*FRR(NTH,NP,NR)
  100    CONTINUE
C
         NR=1
         DO 110 NTH=1,NTHMAX
            RL=0.5D0*(1.D0+RLAMDA(NTH,NR))
         DO 110 NP=1,NPMAX
            DRR(NTH,NP,NR)=RL*DRR(NTH,NP,NR)
            FRR(NTH,NP,NR)=RL*FRR(NTH,NP,NR)
  110    CONTINUE
C
         NR=NRMAX+1
         DO 120 NTH=1,NTHMAX
            RL=RLAMDA(NTH,NR-1)
         DO 120 NP=1,NPMAX
            DRR(NTH,NP,NR)=RL*DRR(NTH,NP,NR)
            FRR(NTH,NP,NR)=RL*FRR(NTH,NP,NR)
  120    CONTINUE
      ENDIF
C
      RETURN
      END
C
C ************************************************************
C
C                    CALCULATION OF FE
C
C ************************************************************
C
      SUBROUTINE FPCALE
C
      INCLUDE 'fpcomm.inc'
C
      DO NR=1,NRMAX
      DO NP=1,NPMAX+1
      DO NTH=1,NTHMAX
         FEPP(NTH,NP,NR)= AEFP*E2(NR)/PTFP0*COSM(NTH)
      ENDDO
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX+1
         FETH(NTH,NP,NR)=-AEFP*E2(NR)/PTFP0*SING(NTH)
      ENDDO
      ENDDO
      ENDDO
C
      IF (MODELA.EQ.0) RETURN
C
      DO NR=1,NRMAX
         FACT=1.D0/SQRT(1.D0-EPSR(NR)**2)
C
         DO NP=1,NPMAX+1
         DO NTH=1,ITL(NR)-1
            FEPP(NTH,NP,NR)= FACT*FEPP(NTH,NP,NR)
         ENDDO
         ENDDO
C
         DO NP=1,NPMAX+1
         DO NTH=ITL(NR),ITU(NR)
            FEPP(NTH,NP,NR)= 0.D0
         ENDDO
         ENDDO
C
         DO NP=1,NPMAX+1
         DO NTH=ITU(NR)+1,NTHMAX
            FEPP(NTH,NP,NR)= FACT*FEPP(NTH,NP,NR)
         ENDDO
         ENDDO
C
         DO NP=1,NPMAX
         DO NTH=1,ITL(NR)
            FETH(NTH,NP,NR)=FACT*FETH(NTH,NP,NR)
         ENDDO
         ENDDO
C
         DO NP=1,NPMAX
         DO NTH=ITL(NR)+1,ITU(NR)
            FETH(NTH,NP,NR)= 0.D0
         ENDDO
         ENDDO
C
         DO NP=1,NPMAX
         DO NTH=ITU(NR)+1,NTHMAX+1
            FETH(NTH,NP,NR)=FACT*FETH(NTH,NP,NR)
         ENDDO
         ENDDO
C
      ENDDO
C
      RETURN
      END
C
C ************************************************
C     WEIGHTING FUNCTION FOR CONVECTION EFFECT
C ************************************************
C
      FUNCTION FPWEGH(X,Y)
C
      IMPLICIT REAL*8(A-F,H,O-Z)
C
      IF(ABS(Y).LT.1.D-70) THEN
         IF(X.GT.0.D0) THEN
            FPWEGH=0.D0
         ELSEIF(X.LT.0.D0) THEN
            FPWEGH=1.D0
         ELSE
            FPWEGH=0.5D0
         ENDIF
      ELSE
         Z=X/Y
         IF(ABS(Z).LT.1.D-5)THEN
            FPWEGH=0.5D0-Z/12.D0+Z**3/720.D0
         ELSE IF(Z.GE.100.D0)THEN
            FPWEGH=1.D0/Z
         ELSE IF(Z.LE.-100.D0)THEN
            FPWEGH=1.D0/Z+1.D0
         ELSE
            FPWEGH=1.D0/Z-1.D0/(EXP(Z)-1.D0)
         END IF
      ENDIF
      RETURN
      END
