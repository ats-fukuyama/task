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
         CALL FPCALV
      ELSEIF(MODELW.EQ.2) THEN
         CALL FPCALS
      ELSEIF(MODELW.EQ.3) THEN
         CALL FPCALS
      ELSEIF(MODELW.EQ.4) THEN
         CALL FPCALV
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
C ************************************************************
C
C       CALCULATION OF DC AND FC
C
C ************************************************************
C
      SUBROUTINE FPCALC
C
      INCLUDE 'fpcomm.inc'
C
      DO NR=1,NRMAX
         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            DCPP(NTH,NP,NR)=0.D0
            DCPT(NTH,NP,NR)=0.D0
            FCPP(NTH,NP,NR)=0.D0
         ENDDO
         ENDDO
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            DCTP(NTH,NP,NR)=0.D0
            DCTT(NTH,NP,NR)=0.D0
            FCTH(NTH,NP,NR)=0.D0
         ENDDO
         ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         DO NS=1,NSMAX
            IF(MODELC.LT.0) THEN
               IF(NS.EQ.NSFP) THEN
                  CALL FPCALC_L(NR,NS)
               ENDIF
            ELSEIF(MODELC.EQ.0) THEN
               CALL FPCALC_L(NR,NS)
            ELSEIF(MODELC.EQ.1) THEN
               IF(NS.EQ.NSFP) THEN
                  CALL FPCALC_NL(NR,NS)
               ELSE
                  CALL FPCALC_L(NR,NS)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
C
      IF (MODELA.NE.0) THEN
C
         DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               FACT=RLAMDA(NTH,NR)
               DO NP=1,NPMAX+1
                  DCPP(NTH,NP,NR)=FACT*DCPP(NTH,NP,NR)
                  FCPP(NTH,NP,NR)=FACT*FCPP(NTH,NP,NR)
               ENDDO
            ENDDO
C
            DO NTH=1,NTHMAX+1
               FACT=RLAMDC(NTH,NR)
               DO NP=1,NPMAX
                  DCTT(NTH,NP,NR)=FACT*DCTT(NTH,NP,NR)
               ENDDO
            ENDDO
         ENDDO
C
         DO NR=1,NRMAX
            DO NP=1,NPMAX+1
               DCPP(ITL(NR),NP,NR)
     &              =RLAMDA(ITL(NR),NR)/4.D0
     &                    *( DCPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                      +DCPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                      +DCPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                      +DCPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
C
               FCPP(ITL(NR),NP,NR)
     &              =RLAMDA(ITL(NR),NR)/4.D0
     &                    *( FCPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                      +FCPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                      +FCPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                      +FCPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
               DCPP(ITU(NR),NP,NR)=DCPP(ITL(NR),NP,NR)
               FCPP(ITU(NR),NP,NR)=FCPP(ITL(NR),NP,NR)
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
C
C ************************************************************
C
C       CALCULATION OF LINEAR COLLISIONAL OPERATOR
C
C ************************************************************
C
      SUBROUTINE FPCALC_L(NR,NS)
C
      INCLUDE 'fpcomm.inc'
      COMMON /FPFNV1/ PNFD,TMC2FD
      EXTERNAL FPFN0R,FPFN1R,FPFN2R,FPFN3R,FPFN4R,FPFN5R,FPFN6R
C
      AMFD=PA(NS)*AMP
      PTFPL=PTFP(NR)
      VTFPL=VTFP(NR)
      RNNL=RNFD(NR,NS)/RNFP0
      RNUFL=RNUF(NR,NS)*RNNL
      RNUDL=RNUD(NR,NS)*RNNL
      PTFDL=PTFD(NR,NS)
      VTFDL=VTFD(NR,NS)
C
C     ----- Non-Relativistic -----
C
      IF(MODELR.EQ.0) THEN
         DO NP=1,NPMAX+1
            IF(NP.EQ.1) THEN
               DCPPL=RNUDL*(2.D0/(3.D0*SQRT(PI)))
     &                    *(VTFP0/(SQRT(2.D0)*VTFD(NR,NS)))
               DCPPL=0.D0
               FCPPL=0.D0
            ELSE
               PFPL=PG(NP)*PTFP0
               VFPL=PFPL/AMFP
               V=VFPL/VTFP0
               U=VFPL/(SQRT(2.D0)*VTFD(NR,NS))
               DCPPL= 0.5D0*RNUDL/U   *(ERF0(U)/U**2-ERF1(U)/U)
               FCPPL=-      RNUFL/U**2*(ERF0(U)-U*ERF1(U))
            ENDIF
            DO NTH=1,NTHMAX
               DCPP(NTH,NP,NR)=DCPP(NTH,NP,NR)+DCPPL
               FCPP(NTH,NP,NR)=FCPP(NTH,NP,NR)+FCPPL
            ENDDO
         ENDDO
C
         DO NP=1,NPMAX
            PFPL=PM(NP)*PTFP0
            VFPL=PFPL/AMFP
            V=VFPL/VTFP0
            U=VFPL/(SQRT(2.D0)*VTFD(NR,NS))
C            DCTTL= 0.25D0*RNUDL/V 
C     &                   *((2.D0-1.D0/U**2)*ERF0(U)+ERF1(U)/U)
            DCTTL= 0.25D0*RNUDL/U
     &                   *((2.D0-1.D0/U**2)*ERF0(U)+ERF1(U)/U)
            IF(NS.EQ.NSFP.AND.MODELC.LT.0) THEN
               DCTTL=DCTTL+0.5D0*ZEFF*RNUDL/U
            ENDIF
            DO NTH=1,NTHMAX+1
               DCTT(NTH,NP,NR)=DCTT(NTH,NP,NR)+DCTTL
            ENDDO
         ENDDO
C
C     ----- Relativistic -----
C
      ELSE
         DO NP=1,NPMAX+1
            IF(NP.EQ.1) THEN
               DCPPL=RNUDL*(2.D0/(3.D0*SQRT(PI)))
     &                    *(VTFP0/(SQRT(2.D0)*VTFD(NR,NS)))
               FCPPL=0.D0
            ELSE
               PFPL=PG(NP)*PTFP0
               VFPL=PFPL/SQRT(AMFP**2+PTFPL**2/VC**2)
               VFDL=VFPL
               PFDL=AMFD*VFDL/SQRT(1.D0-VFDL**2/VC**2)
               PNFDL=PFDL/PTFDL
               PNFD=PNFDL
               TMC2FD=(PTFDL/(AMFD*VC))**2
               CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R)
               CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R)
               CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R)
               WRITE(6,'(I5,1P5E12.4)') NR,PNFD,PTFDL,AMFD,VC,TMC2FD
               DCPPL= 4.D0*PI*RNUDL*VTFP0/3.D0
     &              *(RINT1*VTFDL**2/(RINT0*VTFPL**3)
     &               +RINT2/(RINT0*VTFDL))
               WRITE(6,'(I5,1P4E12.4)') NR,RINT0,RINT1,RINT2,DCPPL
               CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R)
               CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R)
               CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R)
               FCPPL=-4.D0*PI*RNUDL*VTFP0**2/3.D0
     &              *(3.D0*RINT4*VTFDL   /(RINT0*VTFPL**2)
     &               -     RINT5*VTFDL**3/(RINT0*VTFPL**2*VC**2)
     &               +2.D0*RINT6*VTFPL   /(RINT0*VC**2))
            ENDIF
            DO NTH=1,NTHMAX
               DCPP(NTH,NP,NR)=DCPP(NTH,NP,NR)+DCPPL
               FCPP(NTH,NP,NR)=FCPP(NTH,NP,NR)+FCPPL
            ENDDO
         ENDDO
C
         DO NP=1,NPMAX
            PFPL=PM(NP)*PTFP0
            VFPL=PFPL/SQRT(AMFP**2+PTFPL**2/VC**2)
            VFDL=VFPL
            PFDL=AMFD*VFDL/SQRT(1.D0-VFDL**2/VC**2)
            PNFDL=PFDL/PTFDL
            PNFD=PNFDL
            TMC2FD=PTFDL**2/(AMFD*VC)**2
            CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R)
            CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R)
            CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R)
            CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R)
            DCTTL= 4.D0*PI*RNUDL*VTFP0/3.D0
     &           *(1.5D0*RINT3*VTFPL**2/(RINT0*VTFPL**3)
     &            -0.5D0*RINT1*VTFDL**2/(RINT0*VTFPL**3)
     &            +RINT2/(RINT0*VTFDL))
            IF(MODELC.EQ.-1) THEN
               V=VFPL/VTFP0
               DCTTL=DCTTL+0.5D0*ZEFF*RNUDL/V
            ENDIF
C
            DO NTH=1,NTHMAX+1
               DCTT(NTH,NP,NR)=DCTT(NTH,NP,NR)+DCTTL
            ENDDO
         ENDDO
      ENDIF
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
      REAL*8 FUNCTION FPFN0R(X)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFD,TMC2FD
C
      A=PNFD
      PN=A*X
      FPFN0R=A*PN**2*FPRMXW(PN)
C
      RETURN
      END
C
      REAL*8 FUNCTION FPFN1R(X,XM,XP)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFD,TMC2FD
C
      XX=XM
      XX=X
      A=0.5D0*PNFD
      PN=A*XP
      B=PN**4/(1.D0+PN**2*TMC2FD)
      FPFN1R=A*B*FPRMXW(PN)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPFN2R(X)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFD,TMC2FD
C
      A=PNFD
      PN=A*(X+1.D0)
      B=PN*SQRT(1.D0+PN**2*TMC2FD)
      FPFN2R=A*B*FPRMXW(PN)
C
      RETURN
      END
C
C ==============================================================
C
      REAL*8 FUNCTION FPFN3R(X,XM,XP)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFD,TMC2FD
C
      XX=X
      XX=XM
      A=0.5D0*PNFD
      PN=A*XP
      FPFN3R=A*PN**2*FPRMXW(PN)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPFN4R(X,XM,XP)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFD,TMC2FD
C
      XX=X
      XX=XM
      A=0.5D0*PNFD
      PN=A*XP
      B=PN**2/SQRT(1.D0+PN**2*TMC2FD)
      FPFN4R=A*B*FPRMXW(PN)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPFN5R(X,XM,XP)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFD,TMC2FD
C
      XX=X
      XX=XM
      A=0.5D0*PNFD
      PN=A*XP
      B=PN**4/(SQRT(1.D0+PN**2*TMC2FD))**3
      FPFN5R=A*B*FPRMXW(PN)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPFN6R(X)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFD,TMC2FD
C
      A=PNFD
      PN=A*(X+1.D0)
      FPFN6R=A*PN*FPRMXW(PN)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPRMXW(PN)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFD,TMC2FD
C
      EX=(1.D0-SQRT(1.D0+PN**2*TMC2FD))/TMC2FD
C      WRITE(6,'(A,1P3E12.4)') 'FPRMXW: ',PN,TMC2FD,EX
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
