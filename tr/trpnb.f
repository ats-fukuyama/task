C     $Id$
C
C     ***********************************************************
C
C           NEUTRAL BEAM
C
C     ***********************************************************
C
      SUBROUTINE TRPWNB
C
      INCLUDE 'trcomm.inc'
C
      IF(MDLNB.EQ.0) THEN
         DO NR=1,NRMAX
            TAUB(NR)=1.D0
         ENDDO
      ELSEIF(MDLNB.EQ.1) THEN
         CALL TRNBIA
         DO NR=1,NRMAX
            SNB(NR)=0.D0
         ENDDO
         CALL TRAJNB
      ELSEIF(MDLNB.EQ.2) THEN
         CALL TRNBIA
         CALL TRAJNB
      ELSEIF(MDLNB.EQ.3) THEN
         CALL TRNBIB
         DO NR=1,NRMAX
            SNB(NR)=0.D0
         ENDDO
         CALL TRAJNB
      ELSEIF(MDLNB.EQ.4) THEN
         CALL TRNBIB
         CALL TRAJNB
      ENDIF
      RETURN
      END
C
C     ***********************************************************
C
C           NEUTRAL BEAM : (GAUSSIAN PROFILE)
C
C     ***********************************************************
C
      SUBROUTINE TRNBIA
C
      INCLUDE 'trcomm.inc'
C
      IF(PNBTOT.LE.0.D0) RETURN
C
      SUM = 0.D0
      DO 10 NR=1,NRMAX
         SUM=SUM+DEXP(-((RA*RM(NR)-PNBR0)/PNBRW)**2)*DVRHO(NR)*DR
   10 CONTINUE
C
      PNB0=PNBTOT*1.D6/SUM
C
      DO 20 NR=1,NRMAX
         PNB(NR)=PNB0*DEXP(-((RA*RM(NR)-PNBR0)/PNBRW)**2)
         SNB(NR)=PNB(NR)/(PNBENG*RKEV)*1.D-20
   20 CONTINUE
C
      RETURN
      END
C
C     ***********************************************************
C
C           NEUTRAL BEAM : (PENCIL BEAM)
C
C     ***********************************************************
C
      SUBROUTINE TRNBIB
C
      INCLUDE 'trcomm.inc'
C
      DIMENSION AP(10),AR(10)
      DATA AP/0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,
     &        5.0D0,6.0D0,7.0D0,8.0D0,9.0D0/
      DATA AR/0.02D0,0.06D0,0.12D0,0.14D0,0.16D0,
     &        0.16D0,0.14D0,0.12D0,0.08D0,0.02D0/
C
      IF(PNBTOT.LE.0.D0) RETURN
C
      DO 100 NR=1,NRMAX
         SNB(NR) = 0.D0
  100 CONTINUE
C
      NRNBMAX=10
      DO J=1,NRNBMAX
         RWD=PNBRW*DBLE(J-1)/(RA*DBLE(NRNBMAX-1))
         RDD=AR(J)
         CALL TRNBPB(PNBRTG,RWD,RDD)
      ENDDO
C
      DO NR=1,NRMAX
         PNB(NR) = SNB(NR)*1.D20*PNBENG*RKEV
      ENDDO
      RETURN
      END
C
C     ***********************************************************
C
C           NEUTRAL BEAM CALCULATION
C
C     ***********************************************************
C
      SUBROUTINE TRNBPB(R0,RWD,RDD)
C
      INCLUDE 'trcomm.inc'
C
      IF(PNBTOT.LE.0.D0) RETURN
C
      XL = SQRT((RR+RA)**2-R0**2)
     &    *(1.D0-RWD*RWD/((RR+RA)**2-R0**2))
      ANL=RDD*PNBTOT*1.D6/(PNBENG*RKEV*1.D20)
      IB=INT(RWD/DR)+1
      I=NRMAX-IB+1
C
      SUML=0.D0
      DL=XL
C
   10 IF(I.GT.0) THEN
         IM=I+IB-1
      ELSE
         IM=1-I+IB
      ENDIF
C
C      WRITE(6,*) IB,I,IM,SUML,DL
C
      IF(IM.GE.1.AND.IM.LE.NRMAX) THEN
C
      IF(I.GT.0) THEN
         IF(IM.GT.IB) THEN
            RND=RR+RA*SQRT(RM(IM)**2-RM(IB)**2)
         ELSE
            RND=RR+RA*RM(IM)
         ENDIF
      ELSE
         IF(IM.GT.IB) THEN
            RND=RR-RA*SQRT(RM(IM)**2-RM(IB)**2)
         ELSE
            RND=RR-RA*RM(IM)
         ENDIF
      ENDIF
C
      DEX = RN(IM,1)
      TEX = RT(IM,1)
      DCX = ANC(IM)
      DFX = ANFE(IM)
      DOX = 0.D0
      ZCX = PZC(IM)
      ZFX = PZFE(IM)
      ZOX = 8.D0
C
      IF(RND-R0.GT.0.D0) THEN
         IF(I.GT.0) THEN
            IF(IM.LE.NRMAX-1) THEN
               DRR1=RM(IM+1)**2-RM(IB)**2
               DRR2=RM(IM  )**2-RM(IB)**2
               IF(DRR1.LT.0.D0) DRR1=0.D0
               IF(DRR2.LT.0.D0) DRR2=0.D0
               DRM1=(RR+RA*SQRT(DRR1))**2-R0**2
               DRM2=(RR+RA*SQRT(DRR2))**2-R0**2
               IF(DRM1.LT.0.D0) DRM1=0.D0
               IF(DRM2.LT.0.D0) DRM2=0.D0
               DL = SQRT(DRM1)-SQRT(DRM2)
            ELSE
               DL=0.D0
            ENDIF
         ELSE
            IF(IM.GE.2) THEN
               DRR1=RM(IM-1)**2-RM(IB)**2
               DRR2=RM(IM  )**2-RM(IB)**2
               IF(DRR1.LT.0.D0) DRR1=0.D0
               IF(DRR2.LT.0.D0) DRR2=0.D0
               DRM1=(RR-RA*SQRT(DRR1))**2-R0**2
               DRM2=(RR-RA*SQRT(DRR2))**2-R0**2
               IF(DRM1.LT.0.D0) DRM1=0.D0
               IF(DRM2.LT.0.D0) DRM2=0.D0
               DL = SQRT(DRM1)-SQRT(DRM2)
            ELSE
               DL=0.D0
            ENDIF
         ENDIF
      ELSE
         IF(I.GT.0) THEN
            DRR1=RM(IM+1)**2-RM(IB)**2
            IF(DRR1.LT.0.D0) DRR1=0.D0
            DRM1=(RR+RA*SQRT(DRR1))**2-R0**2
            IF(DRM1.LT.0.D0) DRM1=0.D0
            DL = 2.D0*SQRT(DRM1)
         ELSE
            DRR1=RM(IM)**2-RM(IB)**2
            IF(DRR1.LT.0.D0) DRR1=0.D0
            DRM1=(RR+RA*SQRT(DRR1))**2-R0**2
            IF(DRM1.LT.0.D0) DRM1=0.D0
            DL = 2.D0*SQRT(DRM1)
         ENDIF
      ENDIF
C
C
C
C      IF(RND-R0.GT.0.D0) THEN
C         IF(I.GT.0) THEN
C            IF(IM.NE.IB) THEN
C               DL = SQRT((RR+SQRT(RM(IM+1)**2-RM(IB)**2))**2-R0**2)
C     &             -SQRT((RR+SQRT(RM(IM  )**2-RM(IB)**2))**2-R0**2)
C            ELSE
C               DL = SQRT((RR+SQRT(RM(IM+1)**2-RM(IB)**2))**2-R0**2)
C     &             -SQRT( RR                             **2-R0**2)
C            ENDIF
C         ELSE
C            IF(IM.NE.IB.AND.IM-1.NE.IB) THEN
C               DL = SQRT((RR-SQRT(RM(IM-1)**2-RM(IB)**2))**2-R0**2)
C     &             -SQRT((RR-SQRT(RM(IM  )**2-RM(IB)**2))**2-R0**2)
C            ELSEIF(IM.NE.IB) THEN
C               DL = SQRT((RR-SQRT(RM(IM-1)**2-RM(IB)**2))**2-R0**2)
C     &             -SQRT( RR                             **2-R0**2)
C         ELSE
C               DL =-SQRT( RR                             **2-R0**2)
C     &             -SQRT((RR-SQRT(RM(IM  )**2-RM(IB)**2))**2-R0**2)
C               ENDIF
C            ENDIF
C         ELSE
C            IF(I.GT.0) THEN
C               DL = 2.D0*SQRT((RR+SQRT(RM(IM+1)-RM(IB)**2))**2-R0**2)
C            ELSE
C               DL = 2.D0*SQRT((RR+SQRT(RM(IM-1)-RM(IB)**2))**2-R0**2)
C            ENDIF
C         ENDIF
C
      CALL TRBSCS(DEX,TEX,PNBENG,DCX,DFX,DOX,
     &            ZCX,ZFX,ZOX,SGM)
C
      P1=DEX*SGM*ANL*DL
      IF(P1.LT.0.D0) P1=0.D0
      IF(P1.LE.ANL) THEN
         KL=1
      ELSE
         P1=ANL
         KL=0
      ENDIF
C
C      SNB(IM) = SNB(IM)+P1/(DVRHO(IM)*DR)
      IF(IM.GT.1) 
     &     SNB(IM-1) = SNB(IM-1)+0.25D0*P1/(DVRHO(IM-1)*DR)
      SNB(IM) = SNB(IM)+0.5D0*P1/(DVRHO(IM)*DR)
      IF(IM.LE.NRMAX) 
     &     SNB(IM+1) = SNB(IM+1)+0.25D0*P1/(DVRHO(IM+1)*DR)
C
      IF(KL.EQ.0) RETURN
      ANL=ANL-P1
C
      SUML = SUML+DL
C
      ENDIF
      IF(SUML.GE.XL) I=I+1
      IF(SUML.LT.XL) I=I-1
      IF(I.GT.0) THEN
         IM=I+IB-1
      ELSE
         IM=1-I+IB
      ENDIF
C
C      WRITE(6,*) IB,I,IM,SUML,DL
C
      IF(IM.GE.NRMAX) RETURN
      GOTO 10
C
      END
C
C     ***********************************************************
C
C           CALCULATE BEAM STOPPING CROSS-SECTIONS
C
C     ***********************************************************
C
      SUBROUTINE TRBSCS(ANEX,TEX,EB,ANCX,ANFEX,ANOX,
     &                  PZCX,PZFEX,PZOX,SGM)
C
      IMPLICIT REAL*8 (A-F,H,O-Z)
C
C     "C" DENOTES CARBON, "F" IRON, AND "O" OXIGEN.
C
      DIMENSION A(2,3,2),C(3,2,2),F(3,2,2),O(3,2,2)
C
      DATA A/ 4.40D+00, 2.30D-01, 7.46D-02, 3.16D-03,
     &       -2.55D-03, 1.32D-03,-2.49D-02,-1.15D-02,
     &        2.27D-03,-2.78D-05,-6.20D-04, 3.38D-05/
C
      DATA C/-1.49D+00, 5.18D-01,-3.36D-02,-1.19D-01,
     &        2.92D-02,-1.79D-03,-1.54D-02, 7.18D-03,
     &        3.41D-04,-1.50D-02, 3.66D-03,-2.04D-04/
C
      DATA F/-1.03D+00, 3.22D-01,-1.87D-02,-5.58D-02,
     &        1.24D-02,-7.43D-04, 1.06D-01,-3.75D-02,
     &        3.53D-03,-3.72D-03, 8.61D-04,-5.12D-05/
C
      DATA O/-1.41D+00, 4.77D-01,-3.05D-02,-1.08D-01,
     &        2.59D-02,-1.57D-03,-4.08D-04, 1.57D-03,
     &        7.35D-04,-1.38D-02, 3.33D-03,-1.86D-04/
C
      S1=0.D0
      S2=0.D0
      S3=0.D0
      S4=0.D0
      DO I=1,2
      DO J=1,3
      DO K=1,2
         S1=S1+A(I,J,K)*(LOG(EB))**(I-1)*(LOG(ANEX))**(J-1)
     &                 *(LOG(TEX))**(K-1)
      ENDDO
      ENDDO
      ENDDO
      DO I=1,3
      DO J=1,2
      DO K=1,2
         S2=S2+C(I,J,K)*(LOG(EB))**(I-1)*(LOG(ANEX))**(J-1)
     &                 *(LOG(TEX))**(K-1)
         S3=S3+F(I,J,K)*(LOG(EB))**(I-1)*(LOG(ANEX))**(J-1)
     &                 *(LOG(TEX))**(K-1)
         S4=S4+O(I,J,K)*(LOG(EB))**(I-1)*(LOG(ANEX))**(J-1)
     &                 *(LOG(TEX))**(K-1)
      ENDDO
      ENDDO
      ENDDO
C
      SGM=EXP(S1)/EB*(1.D0+(ANCX *PZCX *(PZCX -1.D0)*S2
     &                     +ANFEX*PZFEX*(PZFEX-1.D0)*S3
     &                     +ANOX *ANOX *(PZOX -1.D0)*S4)/ANEX)
C
      RETURN
      END
C
C     ***********************************************************
C
C          NEUTRAL BEAM DRIVEN CURRENT
C
C     ***********************************************************
C
      SUBROUTINE TRAJNB
C
      INCLUDE 'trcomm.inc'
C
      AMD=AMM*PA(2)
      AMT=AMM*PA(3)
      AMA=AMM*PA(4)
      AMB=AMM*PA(2)
      PZB=PZ(2)
      VB=SQRT(2.D0*PNBENG*RKEV/AMB)
C
      DO 10 NR=1,NRMAX
         ANE=RN(NR,1)
         TE =RT(NR,1)
         WB =RW(NR,1)!*1.5D0
         IF(ANE.EQ.0.D0) THEN
            P4=0.D0
            TAUS=0.D0
         ELSE
         P4 = 3.D0*SQRT(0.5D0*PI)*AME/ANE
     &       *(ABS(TE)*RKEV/AME)**1.5D0
         VCD3 = P4*RN(NR,2)*PZ(2)**2/AMD
         VCT3 = P4*RN(NR,3)*PZ(3)**2/AMT
         VCA3 = P4*RN(NR,4)*PZ(4)**2/AMA
         VC3  = VCD3+VCT3+VCA3
         VCR  = VC3**(1.D0/3.D0)
         HYB  = HY(VB/VCR)
         TAUS = 0.2D0*PA(2)*ABS(TE)**1.5D0/(PZ(2)**2*ANE*15.D0)
         TAUB(NR) = 0.5D0*TAUS*(1.D0-HYB)
         RNF(NR,1)= 2.D0*LOG(1.D0+(VB/VCR)**3)*WB
     &             /(3.D0*(1.D0-HYB)*PNBENG)
         ENDIF
C
         IF(RNF(NR,1).GT.0.D0) THEN
            RTF(NR,1)= WB/(1.5D0*RNF(NR,1))
         ELSE
            RTF(NR,1)= 0.D0
         ENDIF
         PBIN(NR)   = WB*RKEV*1.D20/TAUB(NR)
         PBCL(NR,1) =   (1.D0-HYB)*PBIN(NR)
         PBCL(NR,2) = VCD3/VC3*HYB*PBIN(NR)
         PBCL(NR,3) = VCT3/VC3*HYB*PBIN(NR)
         PBCL(NR,4) = VCA3/VC3*HYB*PBIN(NR)
   10 CONTINUE
C
      IF(PNBCD.LE.0.D0) RETURN
      IF(MDLUF.NE.0) RETURN
C
      TAUS0=6.D0*PI*SQRT(2.D0*PI)*AEPS0**2*AMB*AME
     &     /(1.D20*AEE**4*PZB**2*15.D0)
      DO 20 NR=1,NRMAX
         ANE=RN(NR,1)
         TE =RT(NR,1)
         EPS = EPSRHO(NR)
         VE  = SQRT(ABS(TE)*RKEV/AME)
       IF(ANE.EQ.0.D0) THEN
          TAUS=0.D0
          ZEFFM=0.D0
          XB=0.D0
          AJNB(NR)=0.D0
       ELSE
         TAUS=TAUS0*VE**3/ANE
         ZEFFM = (PZ(2)  *PZ(2)  *RN(NR,2)/PA(2)
     &           +PZ(3)  *PZ(3)  *RN(NR,3)/PA(3)
     &           +PZ(4)  *PZ(4)  *RN(NR,3)/PA(4)
     &           +PZC(NR) *PZC(NR) *ANC(NR)/12.D0
     &           +PZFE(NR)*PZFE(NR)*ANFE(NR)/52.D0)/ANE
         EC  = 14.8D0*TE*PA(2)*ZEFFM**(2.D0/3.D0)
         VCR = VB*SQRT(ABS(EC)/PNBENG)
         P2  = (1.55D0+0.85D0/ZEFF(NR))*SQRT(EPS)
     &        -(0.2D0+1.55D0/ZEFF(NR))*EPS
         XB  = VB/VCR
         ZN  = 0.8D0*ZEFF(NR)/PA(2)
         P3  = XB*XB/(4.D0+3.D0*ZN+XB*XB*(XB+1.39D0+0.61D0*ZN**0.7D0))
C
         AJNB(NR) = PNBCD*2.D0*AEE*PZB*TAUS/(AMB*VCR)
     &            *(1.D0-PZB*(1.D0-P2)/ZEFF(NR))*P3*PBIN(NR)
       ENDIF
C
   20 CONTINUE
C
      RETURN
      END
C
C     ***********************************************************
C
      FUNCTION HY(V)
C
      INCLUDE 'trcomm.inc'
C
      HY = 2.D0*(LOG((V**3+1.D0)/(V+1.D0)**3)/6.D0
     &      +(ATAN((2.D0*V-1.D0)/SQRT(3.D0))+PI/6.D0)
     &       /SQRT(3.D0))/V**2
      RETURN
      END
