C
C     ****** TASK/W1 (SUBROUTINE LIBRARY 2) ******
C
C     ******* LOCAL PLASMA PARAMETERS *******
C
      SUBROUTINE W1DSPA(NALPHA)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      PARAMETER (NHM=2*NHARMM+1,NHM1=NHARMM+1)
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM3/ PA(ISM),PZ(ISM),PN(ISM),PTPP(ISM),PTPR(ISM),
     &                PNS(ISM),PTS(ISM),PU(ISM),PZCL(ISM)
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1PRF1/ PROFB(NXPM),PROFPN(NXPM,ISM),PROFPU(NXPM,ISM),
     &                PROFTR(NXPM,ISM),PROFTP(NXPM,ISM)
      COMMON /W1MAT1/ CD0(4,NXPM),CD1(2,NXPM),CD2(4,NXPM)
      COMMON /W1MAT2/ CM0(4,NXPM,ISM),CM1(2,NXPM,ISM),CM2(4,NXPM,ISM)
      COMMON /W1XDAT/ XA(NXQ),XAM(NXTM)
      COMMON /W1ZETA/ CZ(NXPM,NHM),CDZ(NXPM,NHM),CDDZ(NXPM,NHM),
     &                GZ(NXPM,NHM)
      COMMON /W1QCTL/ XDMAX,DXD,NDMAX,NHARM,IHARM(ISM),MATL,NMODEL
C
C      CM0 , CD0         CM1 , CD1          CM2 , CD2
C      1   2   0         0   0   1          1   2   0
C    -(2)  3   0         0   0   2        -(2)  3   0
C      0   0   4         1 -(2)  0          0   0   4
C
      DIMENSION RKPP(NXPM),ALAM(0:NHM1)
      DATA CI/(0.D0,1.D0)/
C
      RT2= SQRT(2.D0)
      RW  = 2.D6*PI*RF
      VALPHA=SQRT(2.D0*3.5D6*AEE/(4.D0*AMM))
C
      IF(NMODEL.LE.1) THEN
         RKPR=RKZ
         RNPR=VC*RKPR/RW
         DO 5 NX=1,NXP
            RKPP(NX)=0.D0
    5    CONTINUE
      ELSE
         DO 20 NX=1,NXP
            RKPR=RKZ
            RNPR=VC*RKPR/RW
            CDT =1.D0
            CDX =0.D0
            DO 10 IS=1,ISMAX
               FWP= 1.D20*AEE*AEE*PZ(IS)*PZ(IS)/(AMM*PA(IS)*EPS0)
               FWC= AEE*PZ(IS)*BB/(AMM*PA(IS))
               WP = FWP*PROFPN(NX,IS)
               WC = FWC*PROFB(NX)
               CDT=CDT-WP/(RW*RW-WC*WC)
               CDX=CDX-WP*WC/((RW*RW-WC*WC)*RW)
   10       CONTINUE
            RNPP2=((CDT-RNPR*RNPR)**2-CDX**2)/(CDT-RNPR*RNPR)
            RKPP(NX)= RW*SQRT(MAX(RNPP2,0.D0))/VC
   20    CONTINUE
      ENDIF
C
      DO 1000 IS=1,ISMAX
         FWP= 1.D20*AEE*AEE*PZ(IS)*PZ(IS)/(AMM*PA(IS)*EPS0*RW*RW)
         FWC= AEE*PZ(IS)*BB/(AMM*PA(IS))
         FVT= AEE*1.D3/(AMM*PA(IS))
         DO 100 NX=1,NXP
            CM0(1,NX,IS) = 0.D0
            CM0(2,NX,IS) = 0.D0
            CM0(3,NX,IS) = 0.D0
            CM0(4,NX,IS) = 0.D0
            CM1(1,NX,IS) = 0.D0
            CM1(2,NX,IS) = 0.D0
            CM2(1,NX,IS) = 0.D0
            CM2(2,NX,IS) = 0.D0
            CM2(3,NX,IS) = 0.D0
            CM2(4,NX,IS) = 0.D0
  100    CONTINUE
C
         IF(MOD(NALPHA,2).EQ.1.AND.IS.EQ.4) THEN
            DO 200 NX=1,NXP
               RKPR= RKZ
               WP  = FWP*PROFPN(NX,IS)
               WC  = FWC*PROFB(NX)
               VTE   =SQRT(2.D0*PROFTR(NX,1)*1.D3*AEE/AME)
               VCRIT3=.75D0*SQRT(PI)*AME*(PROFPN(NX,2)/(AMM*PA(2))
     &                                   +PROFPN(NX,3)/(AMM*PA(3)))
     &                /PROFPN(NX,1)
               VCRIT =VTE*VCRIT3**(1.D0/3.D0)
            DO 200 NC=1,2*ABS(IHARM(4))+1
               NN  = NC-ABS(IHARM(4))-1
               VRES=(RW-NN*WC)/RKPR
               IF(VALPHA**2-VRES**2.LE.0.D0
     &           .OR.RKPP(NX).LE.0.D0) GOTO 200
               VLIM=SQRT(VALPHA**2-VRES**2)
               Y1=RKPP(NX)*VCRIT/WC
               Y2=VRES/VCRIT
               Y3=VLIM/VCRIT
               COEF=4.5D0*PI*WP*RW
     &              /(ABS(RKPR)*VCRIT*LOG(1.D0+(VALPHA/VCRIT)**3))
               CM0(1,NX,IS)=CM0(1,NX,IS)
     &                  +CI*COEF*W1FALF(Y1,Y2,Y3,NN,1)
               CM0(2,NX,IS)=CM0(2,NX,IS)
     &                  -   COEF*W1FALF(Y1,Y2,Y3,NN,2)
               CM1(1,NX,IS)=CM1(1,NX,IS)
     &                  +CI*COEF*W1FALF(Y1,Y2,Y3,NN,3)*RW/(RKPP(NX)*VC)
               CM0(3,NX,IS)=CM0(3,NX,IS)
     &                  +CI*COEF*W1FALF(Y1,Y2,Y3,NN,4)
               CM1(2,NX,IS)=CM1(2,NX,IS)
     &                  +   COEF*W1FALF(Y1,Y2,Y3,NN,5)*RW/(RKPP(NX)*VC)
               CM0(4,NX,IS)=CM0(4,NX,IS)
     &                  +CI*COEF*W1FALF(Y1,Y2,Y3,NN,6)
  200       CONTINUE
         ELSE
            IF(NMODEL.LE.3.OR.IHARM(IS).LT.0.OR.IHARM(IS).GE.3) THEN
               IF(NMODEL.LE.3.OR.IHARM(IS).LT.0) THEN
                  IHMIN=0
               ELSE
                  IHMIN=3
               ENDIF
            DO 300 NX=1,NXP
               RKPR= RKZ
               WC  = FWC*PROFB(NX)
               UD  = SQRT(FVT*PROFPU(NX,IS))
               AKPR= RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,IS))
*VDIR NOVECTOR
            DO 300 NC=1,2*ABS(IHARM(IS))+1
               NN= NC-ABS(IHARM(IS))-1
               ARG = (RW-NN*WC)/AKPR
               GZ(NX,NC)= ARG
               CZ(NX,NC)= ARG
  300       CONTINUE
C
            DO 400 NC=1,2*ABS(IHARM(IS))+1
               CALL DSPFNA(CZ(1,NC),CDZ(1,NC),CDDZ(1,NC),NXP)
  400       CONTINUE
C
            DO 500 NX=1,NXP
               RKPR= RKZ
               WP  = FWP*PROFPN(NX,IS)
               WC  = FWC*PROFB(NX)
               UD  = SQRT(FVT*PROFPU(NX,IS))
               AKPR= RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,IS))
               RT  = PROFTP(NX,IS)/PROFTR(NX,IS)
               RL  = RKPP(NX)*RKPP(NX)*FVT*PROFTP(NX,IS)/(WC*WC)
               CALL LAMBDA(ABS(IHARM(IS))+1,RL,ALAM)
            DO 500 NC=1,2*ABS(IHARM(IS))+1
               NN=NC-ABS(IHARM(IS))-1
               IF(ABS(NN).LT.IHMIN) GOTO 500
               CA1 = RW*CZ(NX,NC)/AKPR+0.5D0*(1.D0-RT)*CDZ(NX,NC)
               CA2 = 0.5D0*(RT*RW/WC-NN*(RT-1.D0))*CDZ(NX,NC)
               CA3 =-(RW-NN*WC*(1.D0-1.D0/RT))*GZ(NX,NC)*CDZ(NX,NC)
     &               /AKPR
               ALAMC=ALAM(ABS(NN))
               ALAMM=ALAM(ABS(NN-1))
               ALAMP=ALAM(ABS(NN+1))
               ALAMD=0.5D0*(ALAMM+ALAMP)-ALAMC
               ALAMN=0.5D0*(ALAMM-ALAMP)
               CM0(1,NX,IS)=CM0(1,NX,IS)
     &                     +WP*NN*ALAMN*CA1
               CM0(2,NX,IS)=CM0(2,NX,IS)
     &                     +WP*CI*NN*ALAMD*CA1
               CM1(1,NX,IS)=CM1(1,NX,IS)
     &                     -WP*ALAMN*CA2*RW/(RKPR*VC)
               CM0(3,NX,IS)=CM0(3,NX,IS)
     &                     +WP*(NN*ALAMN-2.0D0*RL*ALAMD)*CA1
               CM1(2,NX,IS)=CM1(2,NX,IS)
     &                     +WP*CI*ALAMD*CA2*RW/(RKPR*VC)
               CM0(4,NX,IS)=CM0(4,NX,IS)
     &                     +WP*ALAMC*CA3
  500       CONTINUE
         ENDIF
         IF(NMODEL.GT.3.AND.IHARM(IS).GE.0) THEN
            DO 600 NX = 1 , NXP
               RKPR = RKZ
               WC   = FWC*PROFB(NX)
               UD   = SQRT(FVT*PROFPU(NX,IS))
               AKPR = RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,IS))
C
               CZ(NX,1) = (RW-(-2.D0)*WC-RKPR*UD)/AKPR
               CZ(NX,2) = (RW-(-1.D0)*WC-RKPR*UD)/AKPR
               CZ(NX,3) = (RW           -RKPR*UD)/AKPR
               CZ(NX,4) = (RW-( 1.D0)*WC-RKPR*UD)/AKPR
               CZ(NX,5) = (RW-( 2.D0)*WC-RKPR*UD)/AKPR
  600       CONTINUE
C
            DO 700 NC=1,5
               CALL DSPFNA(CZ(1,NC),CDZ(1,NC),CDDZ(1,NC),NXP)
  700       CONTINUE
C
            DO 800 NX=1,NXP
               RKPR= RKZ
               RNPR= VC*RKPR/RW
               WP  = FWP*PROFPN(NX,IS)
               WC  = FWC*PROFB(NX)
               UD  = SQRT(FVT*PROFPU(NX,IS))
               AKPR= RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,IS))
               TPR = FVT*PROFTR(NX,IS)
               RNT = PROFTP(NX,IS)/PROFTR(NX,IS)
               TP  = WP*FVT*PROFTP(NX,IS)*(RW/(VC*WC))**2
C
               CC1 = (RW-RKPR*UD)/AKPR
               CC2 = 0.5D0*(1.D0-RNT)
               CC3 =-RKPR*UD/AKPR
               CC4 = 0.5D0*RNT*RW/WC
               CC5 = RW/AKPR
               CC6 =-WC/AKPR
     &               *(1.D0-(1.D0+UD*UD/TPR)/RNT)
               CC7 =-UD*RW/(2.D0*RKPR*TPR)
               CC8 = UD*WC*(1.D0-2.D0/RNT)/(2.D0*RKPR*TPR)
               CC9 = 0.5D0*RW/AKPR
               CC0 =-0.5D0*WC*(1.D0-1.D0/RNT)/AKPR
C
               CW1 = CC1*CZ(NX,1) + CC2*CDZ(NX,1)
               CW2 = CC1*CZ(NX,2) + CC2*CDZ(NX,2)
               CW3 = CC1*CZ(NX,3) + CC2*CDZ(NX,3)
               CW4 = CC1*CZ(NX,4) + CC2*CDZ(NX,4)
               CW5 = CC1*CZ(NX,5) + CC2*CDZ(NX,5)
C
               CF2 =-1.D0*CC3*CZ(NX,2) + (CC4-CC2)*CDZ(NX,2)
               CF3 =                     (CC4    )*CDZ(NX,3)
               CF4 = 1.D0*CC3*CZ(NX,4) + (CC4+CC2)*CDZ(NX,4)
C
               CG2 = (CC5-CC6)*CZ(NX,2) + (CC7-CC8)*CDZ(NX,2)
     &             + (CC9-CC0)*CDDZ(NX,2)
               CG3 = (CC5    )*CZ(NX,3) + (CC7    )*CDZ(NX,3)
     &             + (CC9    )*CDDZ(NX,3)
               CG4 = (CC5+CC6)*CZ(NX,4) + (CC7+CC8)*CDZ(NX,4)
     &             + (CC9+CC0)*CDDZ(NX,4)
C
               CM0(1,NX,IS)=CM0(1,NX,IS)
     &                     +    WP*0.5D0*(CW4+CW2)
               CM0(2,NX,IS)=CM0(2,NX,IS)
     &                     + CI*WP*0.5D0*(CW4-CW2)
               CM0(3,NX,IS)=CM0(3,NX,IS)
     &                     +    WP*0.5D0*(CW4+CW2)
               CM0(4,NX,IS)=CM0(4,NX,IS)
     &                     +    WP      * CG3
               CM1(1,NX,IS)=CM1(1,NX,IS)
     &                     -    WP*0.5D0*(CF4-CF2)/RNPR
               CM1(2,NX,IS)=CM1(2,NX,IS)
     &                     + CI*WP*0.5D0*(CF4-2.D0*CF3+CF2)/RNPR
               CM2(1,NX,IS)=CM2(1,NX,IS)
     &                     +    TP*0.5D0*(CW5-CW4-CW2+CW1)
               CM2(2,NX,IS)=CM2(2,NX,IS)
     &                     + CI*TP*0.5D0*(CW5-2.D0*CW4+2.D0*CW2-CW1)
               CM2(3,NX,IS)=CM2(3,NX,IS)
     &                     +    TP*0.5D0
     &                      *(CW5-3.D0*CW4+4.D0*CW3-3.D0*CW2+CW1)
               CM2(4,NX,IS)=CM2(4,NX,IS)
     &                     +    TP*0.5D0*(CG4-2.D0*CG3+CG2)
C
  800       CONTINUE
         ENDIF
         ENDIF
 1000 CONTINUE
C
      DO 1100 NX = 1 , NXP
         CD0( 1 , NX ) =  1.D0 - RNPR*RNPR
         CD0( 2 , NX ) =  0.D0
         CD0( 3 , NX ) =  1.D0 - RNPR*RNPR
         CD0( 4 , NX ) =  1.D0
         CD1( 1 , NX ) =  RNPR
         CD1( 2 , NX ) =  0.D0
         CD2( 1 , NX ) =  0.D0
         CD2( 2 , NX ) =  0.D0
         CD2( 3 , NX ) = -1.D0
         CD2( 4 , NX ) = -1.D0
 1100 CONTINUE
C
      DO 1200 IS = 1 , ISMAX
      DO 1200 NX = 1 , NXP
         CD0( 1 , NX ) = CD0( 1 , NX ) + CM0( 1 , NX , IS )
         CD0( 2 , NX ) = CD0( 2 , NX ) + CM0( 2 , NX , IS )
         CD0( 3 , NX ) = CD0( 3 , NX ) + CM0( 3 , NX , IS )
         CD0( 4 , NX ) = CD0( 4 , NX ) + CM0( 4 , NX , IS )
         CD1( 1 , NX ) = CD1( 1 , NX ) + CM1( 1 , NX , IS )
         CD1( 2 , NX ) = CD1( 2 , NX ) + CM1( 2 , NX , IS )
         CD2( 1 , NX ) = CD2( 1 , NX ) + CM2( 1 , NX , IS )
         CD2( 2 , NX ) = CD2( 2 , NX ) + CM2( 2 , NX , IS )
         CD2( 3 , NX ) = CD2( 3 , NX ) + CM2( 3 , NX , IS )
         CD2( 4 , NX ) = CD2( 4 , NX ) + CM2( 4 , NX , IS )
 1200 CONTINUE
C
      RETURN
      END
C
C     ****** FUNCTION FOR ALPHA PARTICLE ABSORPTION ******
C
      DOUBLE PRECISION FUNCTION W1FALF(Y1,Y2,Y3,NN,ID)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /W1ALF1/ FY1,FY2,FY3,NY1,NY2
C
      EXTERNAL W1FNA1,W1FNA2
C
      IF(Y3.GT.2.D0) THEN
         A=2.D0*ACOSH(0.5D0*Y3)
         FY1=Y1
         FY2=Y2
         FY3=A
         NY1=NN
         NY2=ID
         H0=0.5D0
         EPS=1.D-6
         ILST=0
         CALL DEFT(W1FNA1,CS,ES,H0,EPS,ILST)
         W1FALF=0.5D0*A*CS/SINH(0.5D0*A)
      ELSE
         FY1=Y1
         FY2=Y2
         FY3=Y3
         NY1=NN
         NY2=ID
         H0=0.5D0
         EPS=1.D-6
         ILST=0
         CALL DEFT(W1FNA2,CS,ES,H0,EPS,ILST)
         W1FALF=0.5D0*Y3*CS
      ENDIF
      RETURN
      END
C
C     ****** SLAVE FUNCTION FOR DE INTEGRATEION ******
C
      DOUBLE PRECISION FUNCTION W1FNA1(X,XP,XM)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /W1ALF1/ FY1,FY2,FY3,NY1,NY2
C
      Y1=FY1
      Y2=FY2
      A =FY3
      NN=NY1
      ID=NY2
      XX=SINH(0.5D0*A*XP)/SINH(0.5D0*A)
      DUMMY=X
      DUMMY=XM
      ARG=Y1*XX
      IF(ID.EQ.1) THEN
         BES=NN*NN*XX*DJBES(NN,ARG)**2/(Y1*Y1)
      ELSEIF(ID.EQ.2) THEN
         BES=NN*XX*XX*DJBES(NN,ARG)
     &       *(DJBES(NN-1,ARG)-DJBES(NN+1,ARG))*0.5D0/Y1
      ELSEIF(ID.EQ.3) THEN
         BES=NN*Y2*XX*DJBES(NN,ARG)**2/Y1
      ELSEIF(ID.EQ.4) THEN
         BES=0.25D0*XX*XX*XX*(DJBES(NN-1,ARG)-DJBES(NN+1,ARG))**2
      ELSEIF(ID.EQ.5) THEN
         BES=Y2*XX*XX*DJBES(NN,ARG)
     &       *(DJBES(NN-1,ARG)-DJBES(NN+1,ARG))*0.5D0
      ELSEIF(ID.EQ.6) THEN
         BES=Y2*Y2*XX*DJBES(NN,ARG)**2
      ENDIF
      Y=SQRT(XX*XX+Y2*Y2)
      W1FNA1=COSH(0.5D0*A*XP)*Y*BES/(1.D0+Y**3)**2
      RETURN
      END
C
C     ****** SLAVE FUNCTION FOR DE INTEGRATEION ******
C
      DOUBLE PRECISION FUNCTION W1FNA2(X,XP,XM)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /W1ALF1/ FY1,FY2,FY3,NY1,NY2
C
      Y1=FY1
      Y2=FY2
      Y3=FY3
      NN=NY1
      ID=NY2
      XX=0.5D0*Y3*XP
      DUMMY=X
      DUMMY=XM
      ARG=Y1*XX
      IF(ID.EQ.1) THEN
         BES=NN*NN*XX*DJBES(NN,ARG)**2/(Y1*Y1)
      ELSEIF(ID.EQ.2) THEN
         BES=NN*XX*XX*DJBES(NN,ARG)
     &       *(DJBES(NN-1,ARG)-DJBES(NN+1,ARG))*0.5D0/Y1
      ELSEIF(ID.EQ.3) THEN
         BES=NN*Y2*XX*DJBES(NN,ARG)**2/Y1
      ELSEIF(ID.EQ.4) THEN
         BES=0.25D0*XX*XX*XX*(DJBES(NN-1,ARG)-DJBES(NN+1,ARG))**2
      ELSEIF(ID.EQ.5) THEN
         BES=Y2*XX*XX*DJBES(NN,ARG)
     &       *(DJBES(NN-1,ARG)-DJBES(NN+1,ARG))*0.5D0
      ELSEIF(ID.EQ.6) THEN
         BES=Y2*Y2*XX*DJBES(NN,ARG)**2
      ENDIF
      Y=SQRT(XX*XX+Y2*Y2)
      W1FNA2=Y*BES/(1.D0+Y**3)**2
      RETURN
      END
