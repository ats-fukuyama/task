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
      INCLUDE 'fpcomm.h'
C
      CALL FPCALE
C
      IF(MODELW.EQ.0) THEN
         CALL FPCALW
      ELSEIF(MODELW.EQ.1) THEN
         CALL FPCALV
      ELSEIF(MODELW.EQ.2) THEN
         CALL FPCALS
      ELSEIF(MODELW.EQ.3) THEN
         CALL FPCALS
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELW =',MODELW
      ENDIF
C
      IF (MODELC.EQ.0)THEN
         CALL FPCALC
      ELSE
         CALL FPCALN
      END IF
C
      DO 10 NR=1,NRMAX
      DO 10 NP=1,NPMAX+1
      DO 10 NTH=1,NTHMAX
         DPP(NTH,NP,NR)=DCPP(NTH,NP,NR)+DWPP(NTH,NP,NR)
         DPT(NTH,NP,NR)=DCPT(NTH,NP,NR)+DWPT(NTH,NP,NR)
         FPP(NTH,NP,NR)=FEPP(NTH,NP,NR)+FCPP(NTH,NP,NR)
   10 CONTINUE
C
      DO 20 NR=1,NRMAX
      DO 20 NP=1,NPMAX
      DO 20 NTH=1,NTHMAX+1
         DTP(NTH,NP,NR)=DCTP(NTH,NP,NR)+DWTP(NTH,NP,NR)
         DTT(NTH,NP,NR)=DCTT(NTH,NP,NR)+DWTT(NTH,NP,NR)
         FTH(NTH,NP,NR)=FETH(NTH,NP,NR)+FCTH(NTH,NP,NR)
   20 CONTINUE
C
C      CALL FPCALR
C
      DO 30 NR=1,NRMAX
      DO 30 NP=1,NPMAX
      DO 30 NTH=1,NTHMAX
         SP(NTH,NP,NR)=0.D0
   30 CONTINUE
C
C     ****************************
C     Boundary condition at p=pmax
C     ****************************
C
      DO 50 NR=1,NRMAX
      DO 50 NTH=1,NTHMAX
C         DPP(NTH,NPMAX,NR)=0.5D0*DPP(NTH,NPMAX,NR)
C         DPT(NTH,NPMAX,NR)=0.5D0*DPT(NTH,NPMAX,NR)
         DPP(NTH,NPMAX+1,NR)=0.D0
         DPT(NTH,NPMAX+1,NR)=0.D0
         FPMAX=FPP(NTH,NPMAX+1,NR)
         FPP(NTH,NPMAX+1,NR)=MAX(FPMAX,0.D0)
   50 CONTINUE
C      DO 60 NR=1,NRMAX
C      DO 60 NTH=1,NTHMAX+1
C         DTP(NTH,NPMAX,NR)=0.D0
C         DTT(NTH,NPMAX,NR)=0.D0
C         FTH(NTH,NPMAX,NR)=0.D0
C   60 CONTINUE
C
C
      DO 70 NR=1,NRMAX
      DO 70 NP=1,NPMAX+1
      DO 70 NTH=1,NTHMAX
         WEIGHP(NTH,NP,NR)=FPWEGH(-DELP*FPP(NTH,NP,NR),
     &                             DPP(NTH,NP,NR))
   70 CONTINUE
C
      DO 80 NR=1,NRMAX
      DO 80 NP=1,NPMAX
      DO 80 NTH=1,NTHMAX+1
         WEIGHT(NTH,NP,NR)=FPWEGH(-DELTH*PM(NP)*FTH(NTH,NP,NR),
     &                             DTT(NTH,NP,NR))
   80 CONTINUE
C
      DO 90 NR=1,NRMAX+1
      DO 90 NP=1,NPMAX
      DO 90 NTH=1,NTHMAX
         WEIGHR(NTH,NP,NR)=FPWEGH(-DELR*FRR(NTH,NP,NR),
     &                             DRR(NTH,NP,NR))
   90 CONTINUE

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
      INCLUDE 'fpcomm.h'
C
      DO 10 NR=1,NRMAX+1
         PSIN=RG(NR)**2
         CALL PLPROF(PSIN)
         RNL=RN(1)/RNE0
         TEL=RTPR(1)/TE0
         FACTR=-2.D0*RG(NR)/RA
      DO 10 NP=1,NPMAX
          FACTP=1.D0/SQRT(1.D0+PG(NP)**2/TEL)
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
      INCLUDE 'fpcomm.h'
C
      DO 10 NR=1,NRMAX
      DO 10 NP=1,NPMAX+1
      DO 10 NTH=1,NTHMAX
         FEPP(NTH,NP,NR)= AEE*E2(NR)/PTH0*COSM(NTH)
   10 CONTINUE
C
      DO 20 NR=1,NRMAX
      DO 20 NP=1,NPMAX
      DO 20 NTH=1,NTHMAX+1
         FETH(NTH,NP,NR)=-AEE*E2(NR)/PTH0*SING(NTH)
   20 CONTINUE
C
      IF (MODELA.EQ.0) RETURN
C
      DO 200 NR=1,NRMAX
         FACT=1.D0/SQRT(1.D0-EPSR(NR)**2)
C
         DO 110 NP=1,NPMAX+1
         DO 110 NTH=1,ITL(NR)-1
            FEPP(NTH,NP,NR)= FACT*FEPP(NTH,NP,NR)
  110    CONTINUE
C
         DO 120 NP=1,NPMAX+1
         DO 120 NTH=ITL(NR),ITU(NR)
            FEPP(NTH,NP,NR)= 0.D0
  120    CONTINUE
C
         DO 130 NP=1,NPMAX+1
         DO 130 NTH=ITU(NR)+1,NTHMAX
            FEPP(NTH,NP,NR)= FACT*FEPP(NTH,NP,NR)
  130    CONTINUE
C
         DO 140 NP=1,NPMAX
         DO 140 NTH=1,ITL(NR)
            FETH(NTH,NP,NR)=FACT*FETH(NTH,NP,NR)
  140    CONTINUE
C
         DO 150 NP=1,NPMAX
         DO 150 NTH=ITL(NR)+1,ITU(NR)
            FETH(NTH,NP,NR)= 0.D0
  150    CONTINUE
C
         DO 160 NP=1,NPMAX
         DO 160 NTH=ITU(NR)+1,NTHMAX+1
            FETH(NTH,NP,NR)=FACT*FETH(NTH,NP,NR)
  160    CONTINUE
C
  200 CONTINUE
C
      RETURN
      END
C
C ************************************************************
C
C       CALCULATION OF DC AND FC (NONRELATIVISTIC VERSION)
C
C ************************************************************
C
      SUBROUTINE FPCALC
C
      INCLUDE 'fpcomm.h'
      EXTERNAL FPFN1R,FPFN2R,FPFN3R,FPFN4R,FPFN5R,FPFN6R
C
      DO 100 NR=1,NRMAX
         DO 10 NP=1,NPMAX+1
         DO 10 NTH=1,NTHMAX
            DCPT(NTH,NP,NR)=0.D0
   10    CONTINUE
         DO 20 NP=1,NPMAX
         DO 20 NTH=1,NTHMAX+1
            DCTP(NTH,NP,NR)=0.D0
            FCTH(NTH,NP,NR)=0.D0
   20    CONTINUE
  100 CONTINUE
C
      DO 5000 NR=1,NRMAX
         RNUL=RNU0*RNE(NR)
         PTHL=PTH0/(SQRT(2.D0)*PTH(NR))
C
         IF(MODELR.EQ.0) THEN
C
            DO 1000 NP=1,NPMAX+1
               IF(NP.EQ.1) THEN
                  DCPPL=(2.D0/(3.D0*SQRT(PI)))*RNUL*PTHL
                  FCPPL=0.D0
               ELSE
                  U=PG(NP)*PTHL
                  DCPPL= 0.5D0*RNUL*PTHL*(ERF1(U)/U**3-ERF2(U)/U**2)
                  FCPPL=-RNUL*PTHL*PTHL*(ERF1(U)/U**2-ERF2(U)/U)
               ENDIF
            DO 1000 NTH=1,NTHMAX
               DCPP(NTH,NP,NR)=DCPPL
               FCPP(NTH,NP,NR) =FCPPL
 1000       CONTINUE
C
            DO 2000 NP=1,NPMAX
               U=PM(NP)*PTHL
               DCTTL= 0.25D0*RNUL*PTHL
     &                      *((2.D0/U-1.D0/U**3)*ERF1(U)+ERF2(U)/U**2)
     &               +0.5D0*RNUL*PTHL*ZEFF/U
            DO 2000 NTH=1,NTHMAX+1
               DCTT(NTH,NP,NR)=DCTTL
 2000       CONTINUE
C
         ELSE
            NRX=NR
            THETAL=THETA(NR)
            RNUR=(RNUL*SQRT(THETA0)**3)/(3.D0*THETAL*DKBSR(NR))
            DO 3000 NP=1,NPMAX+1
               IF(NP.EQ.1) THEN
                  DCPPL=(2.D0/(3.D0*SQRT(PI)))*RNUL*PTHL
                  FCPPL=0.D0
               ELSE
                  PX=PG(NP)
                  PV=PX/SQRT(1.D0+THETA0*PX*PX)
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R)
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R)
                  DCPPL=RNUR*(RINT1/PV**3
     &                       +RINT2)
                  CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R)
                  CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R)
                  CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R)
                  FCPPL=-RNUR*(3.D0*RINT4/PV**2
     &                        -THETA0*RINT5/PV**2
     &                        +2.D0*THETA0*RINT6*PV)
               ENDIF
               DO 3000 NTH=1,NTHMAX
                  DCPP(NTH,NP,NR)=DCPPL
                  FCPP(NTH,NP,NR)=FCPPL
 3000       CONTINUE
C
            DO 4000 NP=1,NPMAX
               PX=PM(NP)
               PV=PX/SQRT(1.D0+THETA0*PX*PX)
               CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R)
               CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R)
               CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R)
               DCTTL=RNUR*(1.5D0*RINT3/PV
     &                    -0.5D0*RINT1/PV**3
     &                    +RINT2)
     &                    +0.5D0*RNUL*ZEFF/PV
C
            DO 4000 NTH=1,NTHMAX+1
                  DCTT(NTH,NP,NR)=DCTTL
 4000       CONTINUE
         ENDIF
C
 5000 CONTINUE
C
      IF (MODELA.EQ.0) RETURN
C
      DO 8000 NR=1,NRMAX
         DO 6000 NTH=1,NTHMAX
            FACT=RLAMDA(NTH,NR)
         DO 6000 NP=1,NPMAX+1
            DCPP(NTH,NP,NR)=FACT*DCPP(NTH,NP,NR)
            FCPP(NTH,NP,NR)=FACT*FCPP(NTH,NP,NR)
 6000    CONTINUE
C
         DO 7000 NTH=1,NTHMAX+1
            FACT=RLAMDC(NTH,NR)
         DO 7000 NP=1,NPMAX
            DCTT(NTH,NP,NR)=FACT*DCTT(NTH,NP,NR)
 7000    CONTINUE
 8000 CONTINUE
C
      DO 9000 NR=1,NRMAX
      DO 9000 NP=1,NPMAX+1
         DCPP(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                    *( DCPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                      +DCPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                      +DCPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                      +DCPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
C
         FCPP(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                    *( FCPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                      +FCPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                      +FCPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                      +FCPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
         DCPP(ITU(NR),NP,NR)=DCPP(ITL(NR),NP,NR)
         FCPP(ITU(NR),NP,NR)=FCPP(ITL(NR),NP,NR)
 9000 CONTINUE
C
      RETURN
      END
C
C ***************************************************************
C
C                       SET OF INTEGRAND
C
C ***************************************************************
C
      REAL*8 FUNCTION FPFN1R(X,XM,XP)
C
      INCLUDE 'fpcomm.h'
C
      XX=XM
      XX=X
      A=0.5D0*PX
      X1=A*XP
      B=X1**4/(1.D0+X1**2*THETA0)
      FPFN1R=A*B*FPRMXW(X1)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPFN2R(X)
C
      INCLUDE 'fpcomm.h'
C
      A=PX
      X1=A*(X+1.D0)
      B=X1*SQRT(1.D0+X1**2*THETA0)
      FPFN2R=A*B*FPRMXW(X1)
C
      RETURN
      END
C
C ==============================================================
C
      REAL*8 FUNCTION FPFN3R(X,XM,XP)
C
      INCLUDE 'fpcomm.h'
C
      XX=X
      XX=XM
      A=0.5D0*PX
      X1=A*XP
      FPFN3R=A*X1**2*FPRMXW(X1)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPFN4R(X,XM,XP)
C
      INCLUDE 'fpcomm.h'
C
      XX=X
      XX=XM
      A=0.5D0*PX
      X1=A*XP
      B=X1**2/SQRT(1.D0+X1**2*THETA0)
      FPFN4R=A*B*FPRMXW(X1)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPFN5R(X,XM,XP)
C
      INCLUDE 'fpcomm.h'
C
      XX=X
      XX=XM
      A=0.5D0*PX
      X1=A*XP
      B=X1**4/(SQRT(1.D0+X1**2*THETA0))**3
      FPFN5R=A*B*FPRMXW(X1)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPFN6R(X)
C
      INCLUDE 'fpcomm.h'
C
      A=PX
      X1=A*(X+1.D0)
      FPFN6R=A*X1*FPRMXW(X1)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPRMXW(X)
C
      INCLUDE 'fpcomm.h'
C
      EX=(1.D0-SQRT(1.D0+X**2*THETA0))/THETA(NRX)
      IF (EX.LT.-100.D0)THEN
         FPRMXW=0.D0
      ELSE
         FPRMXW=EXP(EX)
      ENDIF
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
