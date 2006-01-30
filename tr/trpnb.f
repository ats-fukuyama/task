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
      DO NR=1,NRMAX
         SUM=SUM+DEXP(-((RA*RM(NR)-PNBR0)/PNBRW)**2)*DVRHO(NR)*DR
      ENDDO
C
      PNB0=PNBTOT*1.D6/SUM
C
      DO NR=1,NRMAX
         PNB(NR)=PNB0*DEXP(-((RA*RM(NR)-PNBR0)/PNBRW)**2)
         SNB(NR)=PNB(NR)/(PNBENG*RKEV)*1.D-20
      ENDDO
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
      DIMENSION AR(10)
      DATA AR/0.02D0,0.06D0,0.12D0,0.14D0,0.16D0,
     &        0.16D0,0.14D0,0.12D0,0.06D0,0.02D0/
C
      IF(PNBTOT.LE.0.D0) RETURN
C
      DO NR=1,NRMAX
         SNB(NR) = 0.D0
         DO J=1,10
            GPNB(NR        ,J)=0.0
            GPNB(NR+  NRMAX,J)=0.0
            GPNB(NR+2*NRMAX,J)=0.0
            GPNB(NR+3*NRMAX,J)=0.0
         ENDDO
      ENDDO
C
      NRNBMAX=10
      VY=PNBVY/RA
      DRTG=(PNBRTG-(RR-RA))/(NRNBMAX/2)
      DO J=1,NRNBMAX
C      DO J=1,1
         RDD=AR(J)
C         RDD=1.D0
         DRTG=2.D0*PNBRW/NRNBMAX
         RTG(J)=PNBRTG-PNBRW+0.5D0*DRTG+DRTG*(J-1)
         CALL TRNBPB(J,RTG(J),VY,RDD)
c$$$         DO NR=1,NRMAX
c$$$            PNB(NR) = SNB(NR)*1.D20*PNBENG*RKEV
c$$$         ENDDO
c$$$         CALL TRSUMD(PNB,DVRHO,NRMAX,PNBTOT)
c$$$         write(6,*) PNBTOT*DR*1.D-6
         STOP
      ENDDO
C      STOP
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
      SUBROUTINE TRNBPB(J,R0,VY,RDD)
C
C     J   : J-th NBI
C     R0  : tangential radius of NBI beam (m)
C     VY  : vertical position of NBI (set to zero at mid-plane)
C     RDD : NBI deposition rate
C
C     KL  : judging condition parameter
C           0 : beam energy has decreased at stopping condition
C           1 : beam energy has not decreased at stopping condition yet
C           2 : median center of beam line
C
      INCLUDE 'trcomm.inc'
C
      IF(PNBTOT.LE.0.D0) RETURN
C
C  XL : maximum distance from injection to wall
      XL=2.D0*SQRT((RR+RA)**2-R0**2)
C  ANL : beam intensity [num/s]
C        (the number of beam particles per unit length
C         divided by thier velocity)
      ANL=RDD*PNBTOT*1.D6/(PNBENG*RKEV*1.D20)
C  IB : radial grid point of NB deposition position
      IB=INT(VY/DR)+1
C  I : radial grid point of NB
      I=NRMAX-IB+1
C
C  SUML : total distance in direction of NB
      SUML=0.D0
C  ANL0 : stored ANL for truncation of calculation
      ANL0=ANL
C  DL : arbitrary minute length in direction of NB
      DL=XL/1000
C  COST : cosine between mid-plane and NB
      COST=SQRT((RR+RA)**2-R0**2)/(RR+RA)
C  IDL : serial number of DL
      IDL=0
C
C  IM : radial grid point turned back at magnetic axis
   10 IF(I.GT.0) THEN
         IM=I+IB-1
         write(6,*) "PART1",I,IM,IB
      ELSE
         IF(ABS(I).LE.50) THEN
            IM=-I+IB-1
         ELSEIF(ABS(I).LE.100) THEN
            IM=I-IB+2+2*NRMAX
         ELSE
            IM=-I+IB-1-2*NRMAX
         ENDIF
         write(6,*) "PART2",I,IM,IB
         if(i.lt.-150) stop
      ENDIF
C
C      WRITE(6,'(3(1X,I3),2F15.7)') IB,I,IM,SUML,DL
C
      IF(I.NE.0.AND.(IM.GE.1.AND.IM.LE.NRMAX)) THEN
C
      P1SUM = 0.D0
      IF(I.GT.0) THEN
         RADIUS1 = RR+RG(IM)*RA
C      ELSEIF(I.EQ.0) THEN
C         RADIUS1 = RR
      ELSE
         IF(ABS(I).LE.100) THEN
            RADIUS1 = RR+(DR-RG(IM))*RA 
        ELSE
            RADIUS1 = RR+RG(IM)*RA
         ENDIF
      ENDIF
      IF(I.EQ.-3*NRMAX) THEN
         NLMAX(J)=IDL
         RETURN
      ENDIF
C
      DEX = RN(IM,1)*1.D1
      TEX = RT(IM,1)
      DCX = ANC(IM)*1.D1
      DFX = ANFE(IM)*1.D1
      DOX = 0.D0
      ZCX = PZC(IM)
      ZFX = PZFE(IM)
      ZOX = 8.D0
C
C  Accumulate P1 inside one grid
 20   IDL=IDL+1
      GBL (IDL)  =GUCLIP(SUML)
      GBAN(IDL,J)=GUCLIP(ANL)
C
      CALL TRBSCS(DEX,TEX,PNBENG,DCX,DFX,DOX,ZCX,ZFX,ZOX,SGM)
C
      P1=(DEX*0.1D0)*SGM*ANL*DL
      GBP1(IDL,J)=GUCLIP(P1)
      IF(P1.LT.0.D0) P1=0.D0
      IF(ANL.GT.ANL0*1.D-3) THEN
         KL=1
      ELSE
         P1=ANL
         KL=0
      ENDIF
C
      SUML=SUML+DL
      IF(KL.EQ.1) THEN
         RADIUSG=SQRT((SUML-DL)**2+(RR+RA)*(RR+RA-2.D0*(SUML-DL)*COST))
         GBR (IDL,J)=GUCLIP(RADIUSG)
         GBRH(IDL,J)=GUCLIP(ABS((RADIUSG-RR)/RA))
         RADIUS2=SQRT(SUML**2+(RR+RA)*(RR+RA-2.D0*SUML*COST))
C      write(6,'(2I4,5F13.7)') I,IM,ABS(RADIUS1-RADIUS2),DR*RA
C     &     -ABS(RADIUS1-RADIUS2),RADIUS2-R0,RADIUS1,RADIUS2
      write(6,'(2I4,5F13.7)') I,IM,DR*RA-ABS(RADIUS1-RADIUS2),
     &     RADIUS1,RADIUS2,RADIUS2-R0,SUML
C         write(6,'(2I4,5F13.7)') I,IM,ANL,P1,P1SUM,DR*RA-ABS(RADIUS1
C     &        -RADIUS2),RADIUS2-R0
         IF(RADIUS2-R0.GT.1.D-6) THEN
            IF(DR*RA-ABS(RADIUS1-RADIUS2).GT.1.D-6) THEN
C  inside the grid
               ANL=ANL-P1
               P1SUM=P1SUM+P1
               GOTO 20
            ELSE
C  run off the grid
               P1=P1SUM
               IDL=IDL-1
               SUML=SUML-DL
            ENDIF
         ELSE
C  innermost grid
            write(6,*) "PASS"
            ANL=ANL-P1
            P1=P1SUM
            KL=2
         ENDIF
      ENDIF
C
      SNB(IM) = SNB(IM)+P1/(DVRHO(IM)*DR)
c$$$      IF(IM.GT.1) SNB(IM-1) = SNB(IM-1)+0.25D0*P1/(DVRHO(IM-1)*DR)
c$$$      IF(IM.GT.1.AND.IM.LT.NRMAX) THEN
c$$$         SNB(IM  ) = SNB(IM  )+0.5D0 *P1/(DVRHO(IM  )*DR)
c$$$      ELSEIF(IM.EQ.1.OR.IM.EQ.NRMAX) THEN
c$$$         SNB(IM  ) = SNB(IM  )+       P1/(DVRHO(IM  )*DR)
c$$$      ENDIF
c$$$      IF(IM.LT.NRMAX) SNB(IM+1) = SNB(IM+1)+0.25D0*P1/(DVRHO(IM+1)*DR)
C
C  for graphics
      IF(I.GT.0) THEN
         GPNB(IM        ,J)=GUCLIP(P1/(DVRHO(IM)*DR)*1.D20*PNBENG*RKEV)
      ELSEIF(I.LE.-1.AND.I.GE.-NRMAX) THEN
         GPNB(IM+  NRMAX,J)=GUCLIP(P1/(DVRHO(IM)*DR)*1.D20*PNBENG*RKEV)
      ELSEIF(I.LE.-NRMAX-1.AND.I.GE.-2*NRMAX) THEN
         GPNB(IM+2*NRMAX,J)=GUCLIP(P1/(DVRHO(IM)*DR)*1.D20*PNBENG*RKEV)
      ELSE
         GPNB(IM+3*NRMAX,J)=GUCLIP(P1/(DVRHO(IM)*DR)*1.D20*PNBENG*RKEV)
      ENDIF
C
C      WRITE(6,'(3(1X,I4),4F15.7)') J,I,IM,SUML,DL,ANL,P1
      IF(KL.EQ.0) THEN
         NLMAX(J)=IDL-1
         RETURN
      ENDIF
      IF(KL.EQ.2) THEN
         IF(I.GT.0) THEN
            I=-2*NRMAX-ABS(I)
         ELSE
            I=-2*NRMAX+ABS(I)-1
         ENDIF
         GOTO 10
      ENDIF
C
      ENDIF
C
      IF(SUML.LT.XL) THEN
         I=I-1
C         IF(I.EQ.0) I=I-1
      ELSE
         NLMAX(J)=IDL-1
         RETURN
      ENDIF
C
C      WRITE(6,'(3(1X,I3),4F15.7)') IB,I,IM,SUML,DL,ANL,P1
C
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
C     R.K. Janev, C.D. BOLEY and D.E. Post, Nucl. Fusion 29 (1989) 2125
C     < Input >
C        ANEX  : ELECTRON DENSITY     [10^19]
C        TEX   : ELECTRON TEMPERATURE [keV]
C        EB    : BEAM ENERGY          [keV/u]
C        ANCX  : CARBON DENSITY       [10^19]
C        ANFEX : IRON DENSITY         [10^19]
C        ANOX  : OXIGEN DENSITY       [10^19]
C        PZCX  : CARBON CHARGE NUMBER
C        PZFEX : IRON CHARGE NUMBER
C        PZOX  : OXIGEN CHARGE NUMBER
C     < Output >
C        SGM   : BEAM STOPPING CROSS-SECTION [10^-20 m^2]
C
C     NOTE: "C" DENOTES CARBON, "F" IRON AND "O" OXIGEN.
C
      IMPLICIT NONE
      REAL*8 ANEX,TEX,EB,ANCX,ANFEX,ANOX,PZCX,PZFEX,PZOX,SGM
      INTEGER I,J,K
      REAL*8 S1,S2,S3,S4
      REAL*8 A(2,3,2),C(3,2,2),F(3,2,2),O(3,2,2)
C
C     DATA FORMAT FOR A           DATA FORMAT FOR B
C         111  211  121  221          111  211  311  121
C         131  231  112  212          221  321  112  212
C         122  222  132  232          312  122  222  322
C
      DATA A/ 4.40D+00, 2.30D-01, 7.46D-02,-2.55D-03,
     &        3.16D-03, 1.32D-03,-2.49D-02,-1.15D-02,
     &        2.27D-03,-6.20D-04,-2.78D-05, 3.38D-05/
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
     &                     +ANOX *PZOX *(PZOX -1.D0)*S4)/ANEX)
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
      PAB=PA(2)
      PZB=PZ(2)
      VB=SQRT(2.D0*PNBENG*RKEV/AMB)
C
      IF(MDLUF.NE.0) THEN
      DO NR=1,NRMAX
         ANE=RN(NR,1)
         TE =RT(NR,1)
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
         TAUS = 0.2D0*PAB*ABS(TE)**1.5D0
     &         /(PZ(2)**2*ANE*COULOG(1,2,ANE,TE))
         TAUB(NR) = 0.5D0*TAUS*(1.D0-HYB)
         ENDIF
      ENDDO
      ELSE
      DO NR=1,NRMAX
         ANE=RN(NR,1)
         TE =RT(NR,1)
         WB =RW(NR,1)*1.5D0
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
         TAUS = 0.2D0*PAB*ABS(TE)**1.5D0
     &         /(PZ(2)**2*ANE*COULOG(1,2,ANE,TE))
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
      ENDDO
C
      IF(PNBCD.LE.0.D0) RETURN
C
C     D. R .Mikkelsen and C. E. Singer, J. Plasma Phys. 4 237 (1983)
C        xi_0 corresponds to PNBCD
C        H(r)*P_b/V_p corresponds to PBIN(NR)
C     
      TAUS0=6.D0*PI*SQRT(2.D0*PI)*EPS0**2*AMB*AME
     &     /(1.D20*AEE**4*PZB**2*COULOG(1,2,ANE,TE))
      DO NR=1,NRMAX
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
         EC  = 14.8D0*TE*PAB*ZEFFM**(2.D0/3.D0)
         VCR = VB*SQRT(ABS(EC)/PNBENG)
         P2  = (1.55D0+0.85D0/ZEFF(NR))*SQRT(EPS)
     &        -(0.2D0+1.55D0/ZEFF(NR))*EPS
         XB  = VB/VCR
         ZN  = 0.8D0*ZEFF(NR)/PAB
         P3  = XB*XB/(4.D0+3.D0*ZN+XB*XB*(XB+1.39D0+0.61D0*ZN**0.7D0))
C
         AJNB(NR) = PNBCD*2.D0*AEE*PZB*TAUS/(AMB*VCR)
     &            *(1.D0-PZB*(1.D0-P2)/ZEFF(NR))*P3*PBIN(NR)
      ENDIF
      ENDDO
C
      ENDIF
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
