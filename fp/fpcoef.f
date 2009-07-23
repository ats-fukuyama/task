C     $Id$
C
C ************************************************************
C
C                      CALCULATION OF D AND F
C
C ************************************************************
C
      SUBROUTINE FPCOEF(NSA)
C
      INCLUDE 'fpcomm.inc'

      ISAVE=0
C
C      DO NSA=1,NSAMAX
      DO NR=1,NRMAX
         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=0.D0
            DPT(NTH,NP,NR,NSA)=0.D0
            FPP(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            DTP(NTH,NP,NR,NSA)=0.D0
            DTT(NTH,NP,NR,NSA)=0.D0
            FTH(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
      ENDDO
C      ENDDO

C      DO NSA=1,NSAMAX
C
C     ----- Parallel electric field accleration term -----
C
      CALL FPCALE(NSA)

C
C     ----- Quasi-linear wave-particle interaction term -----
C
      IF(MODELW(NSA).EQ.0) THEN
         CALL FPCALW(NSA)
      ELSEIF(MODELW(NSA).EQ.1) THEN
         CALL FPCALWR(NSA)
      ELSEIF(MODELW(NSA).EQ.2) THEN
         CALL FPCALWR(NSA)
      ELSEIF(MODELW(NSA).EQ.3) THEN
         CALL FPCALWM(NSA)
      ELSEIF(MODELW(NSA).EQ.4) THEN
         CALL FPCALWM(NSA)
      ELSE
         WRITE(6,*) 'XX UNKNOWN MODELW =',MODELW(NSA)
      ENDIF
C
C     ----- Collisional slowing down and diffusion term -----
C
      CALL FPCALC(NSA)

C
C     ----- Sum up velocity diffusion terms -----
C
      DO NR=1,NRMAX
         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=0.D0
            DPT(NTH,NP,NR,NSA)=0.D0
            FPP(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            DTP(NTH,NP,NR,NSA)=0.D0
            DTT(NTH,NP,NR,NSA)=0.D0
            FTH(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO

         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=DCPP(NTH,NP,NR,NSA)+DWPP(NTH,NP,NR,NSA)
            DPT(NTH,NP,NR,NSA)=DCPT(NTH,NP,NR,NSA)+DWPT(NTH,NP,NR,NSA)
            FPP(NTH,NP,NR,NSA)=FEPP(NTH,NP,NR,NSA)+FCPP(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
C
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            DTP(NTH,NP,NR,NSA)=DCTP(NTH,NP,NR,NSA)+DWTP(NTH,NP,NR,NSA)
            DTT(NTH,NP,NR,NSA)=DCTT(NTH,NP,NR,NSA)+DWTT(NTH,NP,NR,NSA)
            FTH(NTH,NP,NR,NSA)=FETH(NTH,NP,NR,NSA)+FCTH(NTH,NP,NR,NSA)
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
      SELECT CASE(MODELS(NSA))
      CASE(0)
         DO NR=1,NRMAX
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  SP(NTH,NP,NR,NSA)=0.D0
               ENDDO
            ENDDO
         ENDDO
      CASE(1)
         PSP=SQRT(2.D0*AMFP(NSA)*SPENG(NSA)*AEE)/PTFP0(NSA)
         SUML=0.D0
         DO NP=1,NPMAX-1
            IF(PG(NP).LE.PSP.AND.PG(NP+1).GT.PSP) THEN
               DO NTH=1,NTHMAX
                  ANGSP=PI*SPANG(NSA)/180.D0
                  IF(THG(NTH).LE.ANGSP.AND.THG(NTH+1).GT.ANGSP) THEN
                     DO NR=1,NRMAX
                        SPL=EXP(-(RM(NR)-SPR0(NSA))**2/SPRW(NSA)**2)
                        SP(NTH,NP,NR,NSA)=SPL*RLAMDA(NTH,NR)
                        SUML=SUML
     &                      +SPL*VOL(NTH,NP)*RLAMDA(NTH,NR)*VOLR(NR)
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         DO NP=1,NPMAX-1
            DO NTH=1,NTHMAX
               DO NR=1,NRMAX
                  SP(NTH,NP,NR,NSA)=SPTOT(NSA)*DELT
     &                             *SP(NTH,NP,NR,NSA)/SUML
               ENDDO
            ENDDO
         ENDDO
         NP=NPMAX
         DO NTH=1,NTHMAX
            DO NR=1,NRMAX
               SP(NTH,NP,NR,NSA)=0.D0
            ENDDO
         ENDDO
      CASE(2)
         PSP=SQRT(2.D0*AMFP(NSA)*3.5D6*AEE)/PTFP0(NSA)
         SUML=0.D0
         DO NP=1,NPMAX-1
            IF(PG(NP).LE.PSP.AND.PG(NP+1).GT.PSP) THEN
               DO NR=1,NRMAX
                  SPL=EXP(-(RM(NR)-SPR0(NSA))**2/SPRW(NSA)**2)
                  DO NTH=1,NTHMAX
                     SP(NTH,NP,NR,NSA)=SPL*RLAMDA(NTH,NR)
                     SUML=SUML+SPL*VOL(NTH,NP)*RLAMDA(NTH,NR)*VOLR(NR)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
         DO NP=1,NPMAX-1
            DO NTH=1,NTHMAX
               DO NR=1,NRMAX
                  SP(NTH,NP,NR,NSA)=SPTOT(NSA)*DELT
     &                             *SP(NTH,NP,NR,NSA)/SUML
               ENDDO
            ENDDO
         ENDDO
         NP=NPMAX
         DO NTH=1,NTHMAX
            DO NR=1,NRMAX
               SP(NTH,NP,NR,NSA)=0.D0
            ENDDO
         ENDDO
      END SELECT
C
C     ****************************
C     Boundary condition at p=pmax
C     ****************************
C
      DO NR=1,NRMAX
      DO NTH=1,NTHMAX
C         DPP(NTH,NPMAX,NR,NSA)=0.5D0*DPP(NTH,NPMAX,NR,NSA)
C         DPT(NTH,NPMAX,NR,NSA)=0.5D0*DPT(NTH,NPMAX,NR,NSA)
C         DPP(NTH,NPMAX+1,NR,NSA)=0.D0
C         DPT(NTH,NPMAX+1,NR,NSA)=0.D0
         FPMAX=FPP(NTH,NPMAX+1,NR,NSA)
         FPP(NTH,NPMAX+1,NR,NSA)=MAX(FPMAX,0.D0)
      ENDDO
      ENDDO
C
C      DO NR=1,NRMAX
C      DO NTH=1,NTHMAX+1
C         DTP(NTH,NPMAX,NR,NSA)=0.D0
C         DTT(NTH,NPMAX,NR,NSA)=0.D0
C         FTH(NTH,NPMAX,NR,NSA)=0.D0
C       ENDDO
C       ENDDO
C

C      ENDDO
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
      SUBROUTINE FPCALR(NSA)
C
      INCLUDE 'fpcomm.inc'
C
      DO NR=1,NRMAX+1
         RHON=RG(NR)
         CALL PLPROF(RHON)
         RTFPL=RTPR(NS_NSA(NSA))/RTFP0(NSA)
         FACTR=-2.D0*RG(NR)/RA
      DO NP=1,NPMAX
         FACTP=1.D0/SQRT(1.D0+PG(NP)**2/RTFPL)
      DO NTH=1,NTHMAX
         DRR(NTH,NP,NR,NSA)=DRR0*FACTP
         FRR(NTH,NP,NR,NSA)=DRR0*FACTP*FACTR
      ENDDO
      ENDDO
      ENDDO
C
      IF (MODELA.EQ.1) THEN
         DO NR=2,NRMAX
         DO NTH=1,NTHMAX
            RL=0.5D0*(RLAMDA(NTH,NR-1)+RLAMDA(NTH,NR))
         DO NP=1,NPMAX
            DRR(NTH,NP,NR,NSA)=RL*DRR(NTH,NP,NR,NSA)
            FRR(NTH,NP,NR,NSA)=RL*FRR(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
         ENDDO
C
         NR=1
         DO NTH=1,NTHMAX
            RL=0.5D0*(1.D0+RLAMDA(NTH,NR))
         DO NP=1,NPMAX
            DRR(NTH,NP,NR,NSA)=RL*DRR(NTH,NP,NR,NSA)
            FRR(NTH,NP,NR,NSA)=RL*FRR(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
C
         NR=NRMAX+1
         DO NTH=1,NTHMAX
            RL=RLAMDA(NTH,NR-1)
         DO NP=1,NPMAX
            DRR(NTH,NP,NR,NSA)=RL*DRR(NTH,NP,NR,NSA)
            FRR(NTH,NP,NR,NSA)=RL*FRR(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
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
      SUBROUTINE FPCALE(NSA)
C
      INCLUDE 'fpcomm.inc'
C
      DO NR=1,NRMAX
      DO NP=1,NPMAX+1
      DO NTH=1,NTHMAX
         FEPP(NTH,NP,NR,NSA)= AEFP(NSA)*E2(NR)/PTFP0(NSA)*COSM(NTH)
      ENDDO
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX+1
         FETH(NTH,NP,NR,NSA)=-AEFP(NSA)*E2(NR)/PTFP0(NSA)*SING(NTH)
      ENDDO
      ENDDO
      ENDDO
C
      IF (MODELA.EQ.0) RETURN
C
      DO NR=1,NRMAX
         IF(MODELA.eq.1)then
C         FACT=1.D0/SQRT(1.D0-EPSR(NR)**2)
         FACT=1.D0
         DO NP=1,NPMAX+1
         DO NTH=1,ITL(NR)-1
            FEPP(NTH,NP,NR,NSA)= FACT*FEPP(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
C
         DO NP=1,NPMAX+1
         DO NTH=ITL(NR),ITU(NR)
            FEPP(NTH,NP,NR,NSA)= 0.D0
         ENDDO
         ENDDO
C
         DO NP=1,NPMAX+1
         DO NTH=ITU(NR)+1,NTHMAX
            FEPP(NTH,NP,NR,NSA)= FACT*FEPP(NTH,NP,NR,NSA)
         ENDDO
         FEPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &           *( FEPP(ITL(NR)-1,NP,NR,NSA)/RLAMDA(ITL(NR)-1,NR)
     &             +FEPP(ITL(NR)+1,NP,NR,NSA)/RLAMDA(ITL(NR)+1,NR)
     &             +FEPP(ITU(NR)-1,NP,NR,NSA)/RLAMDA(ITU(NR)-1,NR)
     &             +FEPP(ITU(NR)+1,NP,NR,NSA)/RLAMDA(ITU(NR)+1,NR))
c         FEPP(ITL(NR),NP,NR,NSA)=1.D0/4.D0
c     &           *( FEPP(ITL(NR)-1,NP,NR,NSA)
c     &             +FEPP(ITL(NR)+1,NP,NR,NSA)
c     &             +FEPP(ITU(NR)-1,NP,NR,NSA)
c     &             +FEPP(ITU(NR)+1,NP,NR,NSA))

         FEPP(ITU(NR),NP,NR,NSA)=FEPP(ITL(NR),NP,NR,NSA)
         ENDDO
C
         DO NP=1,NPMAX
         DO NTH=1,ITL(NR)
            FETH(NTH,NP,NR,NSA)=FACT*FETH(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
C
         DO NP=1,NPMAX
         DO NTH=ITL(NR)+1,ITU(NR)
            FETH(NTH,NP,NR,NSA)= 0.D0
         ENDDO
         ENDDO
C
         DO NP=1,NPMAX
         DO NTH=ITU(NR)+1,NTHMAX+1
            FETH(NTH,NP,NR,NSA)=FACT*FETH(NTH,NP,NR,NSA)
         ENDDO
         ENDDO
C
         ELSE

         DO NTH=1,NTHMAX
            DELH=2.D0*ETAM(NTH,NR)/NAVMAX
            DO NP=1,NPMAX+1
               sum11=0.D0
C
               DO NG=1,NAVMAX
                  ETAL=DELH*(NG-0.5D0)
                  X=EPSR(NR)*COS(ETAL)*RR
                  PSIB=(1.D0+EPSR(NR))/(1.D0+X/RR)
c                  IF (COSM(NTH).GE.0.D0) THEN
                     PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
c                  ELSE
c                     PCOS=-SQRT(1.D0-PSIB*SINM(NTH)**2)
c                  ENDIF
                  sum11=sum11
     &                 +FEPP(NTH,NP,NR,NSA)*ABS(COSM(NTH))/PCOS
               END DO
               FEPP(NTH,NP,NR,NSA)=SUM11*DELH/PI
            END DO
         END DO

         DO NP=1,NPMAX+1
            DO NTH=ITL(NR)+1,NTHMAX/2
               FEPP(NTH,NP,NR,NSA)
     &          =(FEPP(NTH,NP,NR,NSA)+FEPP(NTHMAX-NTH+1,NP,NR,NSA))/2.D0
               FEPP(NTHMAX-NTH+1,NP,NR,NSA)=FEPP(NTH,NP,NR,NSA)
            END DO
         FEPP(ITL(NR),NP,NR,NSA)=1.D0/4.D0
     &           *( FEPP(ITL(NR)-1,NP,NR,NSA)
     &             +FEPP(ITL(NR)+1,NP,NR,NSA)
     &             +FEPP(ITU(NR)-1,NP,NR,NSA)
     &             +FEPP(ITU(NR)+1,NP,NR,NSA))

         FEPP(ITU(NR),NP,NR,NSA)=FEPP(ITL(NR),NP,NR,NSA)
         END DO
c         DO NP=1,NPMAX+1
c         DO NTH=ITL(NR),ITU(NR)
c            FEPP(NTH,NP,NR,NSA)= 0.D0
c         ENDDO
c         ENDDO

         DO NTH=1,NTHMAX+1
            DELH=2.D0*ETAG(NTH,NR)/NAVMAX
            DO NP=1,NPMAX
               sum15=0.D0
C     
               DO NG=1,NAVMAX
                  ETAL=DELH*(NG-0.5D0)
                  X=EPSR(NR)*COS(ETAL)*RR
                  PSIB=(1.D0+EPSR(NR))/(1.D0+X/RR)
                  IF(NTH.NE.NTHMAX/2+1) THEN
                     ARG=1.D0-PSIB*SING(NTH)**2
                     IF(ARG.GT.0.D0) THEN
                        IF (COSG(NTH).GE.0.D0) THEN
                           PCOS= SQRT(ARG)
                        ELSE
                           PCOS=-SQRT(ARG)
                        ENDIF
                     ELSE
                        PCOS=0.D0
                     ENDIF
                  ENDIF
                  SUM15=SUM15
     &                 +FETH(NTH,NP,NR,NSA)/SQRT(PSIB)
               END DO
               FETH(NTH,NP,NR,NSA)=SUM15*DELH/PI
            END DO
         END DO

         DO NP=1,NPMAX
            DO NTH=ITL(NR)+1,NTHMAX/2
               FETH(NTH,NP,NR,NSA)
     &              =(FETH(NTH,NP,NR,NSA)
     &              +FETH(NTHMAX-NTH+2,NP,NR,NSA))/2.D0
               FETH(NTHMAX-NTH+2,NP,NR,NSA)
     &              =FETH(NTH,NP,NR,NSA)
            END DO
         END DO

c         DO NP=1,NPMAX
c         DO NTH=ITL(NR)+1,ITU(NR)
c            FETH(NTH,NP,NR,NSA)= 0.D0
c         ENDDO
c         ENDDO
      END IF

      ENDDO
C
      RETURN
      END
