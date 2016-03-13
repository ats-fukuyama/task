C
C     ****** TASK/W1 (SUBROUTINE LIBRARY 1) ******
C
C     ****** POWER ABSORPTION AS A FUNCTION OF KZ ******
C
      SUBROUTINE W1CLPW(NZ,NSYM)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM4/ RZ,DZ,NZP
      COMMON /W1BND1/ CGIN(3,5),CGOT(3,5),CFJY1,CFJY2,CFJZ1,CFJZ2
      COMMON /W1EF2D/ CE2DA(NZPM,NXTM,3)
      COMMON /W1PWR2/ PAKT(NZPM,3),PANTK(NZPM),PAK(NZPM,ISM)
      COMMON /W1XDAT/ XA(NXQ),XAM(NXTM)
      COMMON /W1PWR0/ PABS(NXPM,ISM),FLUX(NXTM)
      COMMON /W1PWR1/ PABSX(NXPM,ISM),FLUXX(NXTM),PABSXZ(ISM),PABSTT,
     &                RANT1(IAM),RANT2(IAM),XANT1(IAM),XANT2(IAM),
     &                PANT1,PANT2,PANT,PWALL1,PWALL2,PWALL,
     &                PIBW1,PIBW2,PIBW,PERR
      IF(NSYM.NE.0.AND.NZ.NE.1.AND.NZ.NE.(NZP/2+1)) THEN
         FACT = 2.0*RZ
      ELSEIF(NSYM.EQ.-1.AND.NZ.EQ.(NZP/2+1)) THEN
         FACT = 0.0
      ELSE
         FACT =     RZ
      ENDIF
C
      DO 1000 IS=1,ISMAX
      DO 1000 NX=1,NXP
         PABSX(NX,IS)=PABSX(NX,IS)+PABS(NX,IS)*FACT
 1000 CONTINUE
C
      DO 2000 NX=1,NXP
         FLUXX(NX)=FLUXX(NX)+FLUX(NX)*FACT
 2000 CONTINUE
C
      DO 3000 IS=1,ISMAX
         PAK(NZ,IS)=0.D0
      DO 3000 NX=1,NXP
         PAK(NZ,IS)=PAK(NZ,IS)+PABS(NX,IS)*RZ
 3000 CONTINUE
C
      NXANT1=NXP+NXV+NXV/2+1
      NXANT2=NXP    +NXV/2
      PANTK(NZ)=-RZ*DCONJG(CE2DA(NZ,NXANT1,2))*CFJY1
     &          -RZ*DCONJG(CE2DA(NZ,NXANT2,2))*CFJY2
     &          -RZ*DCONJG(CE2DA(NZ,NXANT1,3))*CFJZ2
     &          -RZ*DCONJG(CE2DA(NZ,NXANT2,3))*CFJZ2
C
      PAKT(NZ,3)=0.
      DO 7600 IS=1,ISMAX
         PAKT(NZ,3)=PAKT(NZ,3)+PAK(NZ,IS)
 7600 CONTINUE
C
      PAKT(NZ,2)=0.
      DO 7700 I=1,3
         PAKT(NZ,2)=PAKT(NZ,2)+( FLUXX(NXP+NXV+NXV)
     &                          -FLUXX(1          )
     &                          +FLUXX(NXP        )
     &                          -FLUXX(NXP+1      ))
 7700 CONTINUE
      PAKT(NZ,1)=PAKT(NZ,2)+PAKT(NZ,3)
C
      RETURN
      END
C
C     ****** INITIALIZE POWER ABSORPTION ******
C
      SUBROUTINE W1PWRI
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PWR1/ PABSX(NXPM,ISM),FLUXX(NXTM),PABSXZ(ISM),PABSTT,
     &                RANT1(IAM),RANT2(IAM),XANT1(IAM),XANT2(IAM),
     &                PANT1,PANT2,PANT,PWALL1,PWALL2,PWALL,
     &                PIBW1,PIBW2,PIBW,PERR
      COMMON /W1CDDT/ AJCDX(NXPM),AJCDK(NZPM),AJCDT,ZEFF,WVYSIZ,ICDTYP
C
      DO 10 NX=1,NXT
         FLUXX(NX)=0.
   10 CONTINUE
C
      DO 20 IS=1,ISMAX
      DO 20 NX=1,NXP
         PABSX(NX,IS)=0.
   20 CONTINUE
C
      DO 30 NX=1,NXP
         AJCDX(NX)=0.D0
   30 CONTINUE
C
      RETURN
      END
C
C     ******* POWER ABSORPTION AND ENERGY FLUXX *******
C
      SUBROUTINE W1PWRS
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM3/ PA(ISM),PZ(ISM),PN(ISM),PTPP(ISM),PTPR(ISM),
     &                PNS(ISM),PTS(ISM),PU(ISM),PZCL(ISM)
      COMMON /W1PRM4/ RZ,DZ,NZP
      COMMON /W1ANT1/ AJYL(IAM),AJYH(IAM),AJZL(IAM),AJZH(IAM),
     &                ALYL(IAM),ALYH(IAM),APYL(IAM),APYH(IAM)
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1EF2D/ CE2DA(NZPM,NXTM,3)
      COMMON /W1PWR1/ PABSX(NXPM,ISM),FLUXX(NXTM),PABSXZ(ISM),PABSTT,
     &                RANT1(IAM),RANT2(IAM),XANT1(IAM),XANT2(IAM),
     &                PANT1,PANT2,PANT,PWALL1,PWALL2,PWALL,
     &                PIBW1,PIBW2,PIBW,PERR
      COMMON /W1XDAT/ XA(NXQ),XAM(NXTM)
      COMMON /W1ZDAT/ CJ1(NZPM),CJ2(NZPM),CJ3(NZPM),CJ4(NZPM),
     &                ZA(NZPM),AKZ(NZPM),NZANT1(IAM),NZANT2(IAM),NANT
      COMMON /W1CDDT/ AJCDX(NXPM),AJCDK(NZPM),AJCDT,ZEFF,WVYSIZ,ICDTYP
C
      PABSTT=0.D0
      AJCDT=0.D0
      DO 1410 IS=1,ISMAX
         PABSXZ(IS)=0.
 1410 CONTINUE
C
      DO 1430 IS=1,ISMAX
      DO 1430 NX=1,NXP
         PABSXZ(IS)  =PABSXZ(IS)  +PABSX(NX,IS)
 1430 CONTINUE
C
      DO 1435 IS=1,ISMAX
         PABSX(  1,IS)=2.D0*PABSX(  1,IS)/(XA(2)-XA(1))
         PABSX(NXP,IS)=2.D0*PABSX(NXP,IS)/(XA(NXP)-XA(NXP-1))
      DO 1435 NX=2,NXP-1
         PABSX( NX,IS)=2.D0*PABSX( NX,IS)/(XA(NX+1)-XA(NX-1))
 1435 CONTINUE
C
      DO 1440 IS=1,ISMAX
         PABSTT=PABSTT+PABSXZ(IS)
 1440 CONTINUE
C
      DO 1450 NX=1,NXP
         AJCDT =AJCDT +AJCDX(NX)
 1450 CONTINUE
      AJCDX(  1   )=2.D0*AJCDX(  1   )/(XA(2)-XA(1))
      AJCDX(NXP   )=2.D0*AJCDX(NXP   )/(XA(NXP)-XA(NXP-1))
      DO 1460 NX=2,NXP-1
         AJCDX( NX   )=2.D0*AJCDX( NX   )/(XA(NX+1)-XA(NX-1))
 1460 CONTINUE
C
      NXANT1=NXP+NXV+NXV/2+1
      NXANT2=NXP    +NXV/2
      PANT1 = 0.
      PANT2 = 0.
      DO 1500 NZ = 1 , NZP
         PANT1 = PANT1 + DCONJG(CE2DA(NZ,NXANT1,2))*CJ1(NZ)
     &                 + DCONJG(CE2DA(NZ,NXANT1,3))*CJ3(NZ)
         PANT2 = PANT2 + DCONJG(CE2DA(NZ,NXANT2,2))*CJ2(NZ)
     &                 + DCONJG(CE2DA(NZ,NXANT2,3))*CJ4(NZ)
 1500 CONTINUE
      PANT1 = - PANT1 * DZ
      PANT2 = - PANT2 * DZ
      PANT  = PANT1 + PANT2
C
      PWALL1 = -FLUXX( NXP+NXV+1 )
      PWALL2 =  FLUXX( NXP+NXV   )
      PWALL  =  PWALL1 + PWALL2
C
C     PIBW1 = FLUXX( NXP+NXV+NXV ) - FLUXX(     1 )
      PIBW1 = 0.D0
C     PIBW2 = FLUXX( NXP         ) - FLUXX( NXP+1 )
      PIBW2 = 0.D0
      PIBW  = PIBW1 + PIBW2
C
      PERR  = PANT - ( PABSTT + PWALL + PIBW )
C
      DO 1600 IA = 1 , IAMAX
         PHASE  = APYH(IA)*PI/180.D0
         CPHASE = DCMPLX( DCOS( PHASE ) , DSIN( PHASE ) )
         CJYH   = AJYH( IA ) * CPHASE
         PHASE  = APYL(IA)*PI/180.D0
         CPHASE = DCMPLX( DCOS( PHASE ) , DSIN( PHASE ) )
         CJYL   = AJYL( IA ) * CPHASE
         IF( ABS( CJYH ) .GT. 1.E-6 ) THEN
            RANT1(IA) = -DBLE( CE2DA(NZANT1(IA),NXANT1,2)/ CJYH)
            XANT1(IA) =  IMAG( CE2DA(NZANT1(IA),NXANT1,2)/ CJYH)
         ELSE
            RANT1(IA) = 0.
            XANT1(IA) = 0.
         ENDIF
         IF( ABS( CJYL ) .GT. 1.E-6 ) THEN
            RANT2(IA) = -REAL( CE2DA(NZANT2(IA),NXANT2,2)/ CJYL)
            XANT2(IA) =  IMAG( CE2DA(NZANT2(IA),NXANT2,2)/ CJYL)
         ELSE
            RANT2(IA) = 0.
            XANT2(IA) = 0.
         ENDIF
 1600 CONTINUE
      RETURN
      END
C
C     ****** CALCULATE DRIVEN CURRENT ******
C
      SUBROUTINE W1CLCD(NZ)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM3/ PA(ISM),PZ(ISM),PN(ISM),PTPP(ISM),PTPR(ISM),
     &                PNS(ISM),PTS(ISM),PU(ISM),PZCL(ISM)
      COMMON /W1PRM6/ IELEC(ISM)
      COMMON /W1CTRL/ DRF,DRKZ,DXFACT,DXWDTH,APRFPN,APRFTR,APRFTP,
     &                NPRINT,NFILE,NGRAPH,NLOOP,NSYM,NFLR,NALPHA,
     &                NSYS,NDISP,NXABS
      COMMON /W1XDAT/ XA(NXQ),XAM(NXTM)
      COMMON /W1PRF1/ PROFB(NXPM),PROFPN(NXPM,ISM),PROFPU(NXPM,ISM),
     &                PROFTR(NXPM,ISM),PROFTP(NXPM,ISM)
      COMMON /W1CDDT/ AJCDX(NXPM),AJCDK(NZPM),AJCDT,ZEFF,WVYSIZ,ICDTYP
      COMMON /W1PWR0/ PABS(NXPM,ISM),FLUX(NXTM)
C
      AJCDK(NZ)=0.D0
C
      DO 100 NX=1,NXP
         XR=XAM(NX)/RR
      DO 100 IS=1,ISMAX
      IF(IELEC(IS).EQ.1) THEN
         VTE=SQRT(PROFTR(NX,IS)*AEE*1.D3/AME)
         VPH=2.D0*PI*RF*1.D6/RKZ
         W=ABS(VPH/VTE)
         IF(ABS(VPH).GT.VC) THEN
            EFCD=0.D0
         ELSE
            IF(WVYSIZ.LE.0.D0) THEN
               EFCD=W1CDEF(W,ZEFF,XR,0.D0,ICDTYP)
            ELSE
               Y0=0.00D0*WVYSIZ/RR
               Y1=0.25D0*WVYSIZ/RR
               Y2=0.50D0*WVYSIZ/RR
               Y3=0.75D0*WVYSIZ/RR
               E0=W1CDEF(W,ZEFF,XR,Y0,ICDTYP)
               E1=W1CDEF(W,ZEFF,XR,Y1,ICDTYP)
               E2=W1CDEF(W,ZEFF,XR,Y2,ICDTYP)
               E3=W1CDEF(W,ZEFF,XR,Y3,ICDTYP)
C     *** WEIGHTING WITH PARABOLIC ELECTRIC FIELD PROFILE ***
               EFCD=(256*E0+450*E1+288*E2+98*E3)/1092.D0
C     *** WEIGHTING WITH PARABOLIC POWER DENSITY PROFILE ***
C              EFCD=(16*E0+30*E1+24*E2+14*E3)/84.D0
            ENDIF
         ENDIF
         RLNLMD=16.1D0 - 1.15D0*LOG10(PROFPN(NX,1))
     &                 + 2.30D0*LOG10(PROFTR(NX,1))
         AJCD=0.384D0*PROFTR(NX,IS)*EFCD/(PROFPN(NX,1)*RLNLMD)
     &        *PABS(NX,IS)
         IF(VPH.LT.0.D0) AJCD=-AJCD
         AJCDX(NX)=AJCDX(NX)+AJCD
         AJCDK(NZ)=AJCDK(NZ)+AJCD
      ENDIF
  100 CONTINUE
      RETURN
      END
C
C     ****** CURRENT DRIVE EFFICIENCY ******
C
C      WT = V / VT : PHASE VELOCITY NORMALIZED BY THERMAL VELOCITY
C      Z  = ZEFF   : EFFECTIVE Z
C      XR = X / RR : NORMALIZED X
C      YR = Y / RR : NORMALIZED Y
C      ID : 0 : LANDAU DAMPING
C           1 : TTMP
C
      DOUBLE PRECISION FUNCTION W1CDEF(WT,Z,XR,YR,ID)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      R=SQRT(XR*XR+YR*YR)
      IF(ID.EQ.0) THEN
         D=3.D0/Z
         C=3.83D0
         A=0.D0
         RM=1.38D0
         RC=0.389D0
      ELSE
         D=11.91D0/(0.678D0+Z)
         C=4.13D0
         A=12.3D0
         RM=2.48D0
         RC=0.0987D0
      ENDIF
      IF(WT.LE.1.D-20) THEN
         W=1.D-20
      ELSE
         W=WT
      ENDIF
      EFF0=D/W+C/Z**0.707D0+4.D0*W*W/(5.D0+Z)
      EFF1=1.D0-R**0.77D0*SQRT(3.5D0**2+W*W)/(3.5D0*R**0.77D0+W)
C
      Y2=(R+XR)/(1.D0+R)
      IF(Y2.LT.0.D0) Y2=0.D0
      Y1=SQRT(Y2)
      EFF2=1.D0+A*(Y1/W)**3
C
      IF(Y2.LE.1.D-20) THEN
         YT=(1.D0-Y2)*WT*WT/1.D-60
      ELSE
         YT=(1.D0-Y2)*WT*WT/Y2
      ENDIF
      IF(YT.GE.0.D0.AND.RC*YT.LT.40.D0) THEN
         ARG=(RC*YT)**RM
         IF(ARG.LE.100.D0) THEN
            EFF3=1.D0-MIN(EXP(-ARG),1.D0)
         ELSE
            EFF3=1.D0
         ENDIF
      ELSE
         EFF3=1.D0
      ENDIF
C
      W1CDEF=EFF0*EFF1*EFF2*EFF3
      RETURN
      END
C
C     ****** CALCULATION OF FORCE AND HELICITY ******
C
      SUBROUTINE W1HELD(MODE)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-H,O-Z)
C
      INCLUDE 'w1comm.f'
      PARAMETER (NZPM=2**NZLM)
      PARAMETER (NXTM=NXPM+2*NXVM,NXM=NXPM*6+10,NXQ=NXPM+1)
      COMMON /W1PRM1/ NXP,NXV,NXT,ISMAX,IAMAX
      COMMON /W1PRM2/ RF,RKZ,BB,RR,RA,RD,RB,WALLR
      COMMON /W1PRM3/ PA(ISM),PZ(ISM),PN(ISM),PTPP(ISM),PTPR(ISM),
     &                PNS(ISM),PTS(ISM),PU(ISM),PZCL(ISM)
      COMMON /W1PRF1/ PROFB(NXPM),PROFPN(NXPM,ISM),PROFPU(NXPM,ISM),
     &                PROFTR(NXPM,ISM),PROFTP(NXPM,ISM)
      COMMON /W1CNST/ AEE,AME,AMM,VC,EPS0,AMYU0,PI
      COMMON /W1BAH1/ CEF(NXTM,3),CBF(NXTM,3),CAF(NXTM,4)
      COMMON /W1BAH2/ RHL(NXTM,7),CJF(NXPM,ISM,7),FHL(NXPM,ISM,2)
      COMMON /W1BAH3/ AHLT(ISM,4),AHL(NXPM,ISM,4)
      COMMON /W1BAH4/ CED(NXTM,3),CJD(NXPM,ISM,4)
      COMMON /W1XDAT/ XA(NXQ),XAM(NXTM)
      COMMON /W1EPOL/ CSKX(NXPM*3),CSPX(NXPM*3),CSPZ(NXPM*3)
      COMMON /W1BND2/ CA(NXM)
      COMMON /W1MAT1/ CD0(4,NXPM),CD1(2,NXPM),CD2(4,NXPM)
      COMMON /W1MAT2/ CM0(4,NXPM,ISM),CM1(2,NXPM,ISM),CM2(4,NXPM,ISM)
      COMMON /W1CDDT/ AJCDX(NXPM),AJCDK(NZPM),AJCDT,ZEFF,WVYSIZ,ICDTYP
C
      DATA CI  /( 0.D0 , 1.0D0 )/
C
      RKV=2.D6*PI*RF/VC
      RNPR=RKZ/RKV
      RW=2.D6*PI*RF
C
      DO 100 NX = 1 , NXP
         DX = RKV * ( XA( NX+1 ) - XA( NX ) )
         K  = (NX-1)*MODE + 2
         KK = (NX-1)*MODE/2
C
         DO 30 I=1,3
            CEF(NX,I)=0.D0
            CBF(NX,I)=0.D0
            CAF(NX,I)=0.D0
            CED(NX,I)=0.D0
   30    CONTINUE
         CAF(NX,4)=0.D0
C
         DO 90 IS=1,ISMAX
            DO 70 I=1,7
               CJF(NX,IS,I)=0.D0
   70       CONTINUE
            DO 80 I=1,3
               CJD(NX,IS,I)=0.D0
   80       CONTINUE
   90    CONTINUE
C
         DO 100 I=1,MODE
            CAL=CA(K+I)
            II=(I-1)/2+1
            CSKXL=CSKX(KK+II)
            CSPXL=CSPX(KK+II)
            CSPZL=CSPZ(KK+II)
            IF(MOD(I,2).EQ.0) THEN
               CSKXL=-CSKXL
               CSPZL=-CSPZL
            ENDIF
C
            CPH  = CDEXP(0.5D0*CI*DX*CSKXL)
            CEX  = CAL*CPH*CSPXL
            CEY  = CAL*CPH
            CEZ  = CAL*CPH*CSPZL
            CPT  = CI*(CSKXL*CEX+RNPR*CEZ)/(CSKXL**2+RNPR**2)
C
            CEF(NX,1) = CEF(NX,1)+CEX
            CEF(NX,2) = CEF(NX,2)+CEY
            CEF(NX,3) = CEF(NX,3)+CEZ
C
            CED(NX,1) = CED(NX,1)+CI*CSKXL*CEX
            CED(NX,2) = CED(NX,2)+CI*CSKXL*CEY
            CED(NX,3) = CED(NX,3)+CI*CSKXL*CEZ
C
            CBF(NX,1) = CBF(NX,1)-RNPR*CEY
            CBF(NX,2) = CBF(NX,2)+RNPR*CEX-CSKXL*CEZ
            CBF(NX,3) = CBF(NX,3)+CSKXL*CEY
C
            CAF(NX,1) = CAF(NX,1)-CI*CEX+CSKXL*CPT
            CAF(NX,2) = CAF(NX,2)-CI*CEY
            CAF(NX,3) = CAF(NX,3)-CI*CEZ+RNPR*CPT
            CAF(NX,4) = CAF(NX,4)+CPT
C
            DO 100 IS=1,ISMAX 
               CJX=-CI*( CM0(1,NX,IS)         *CEX
     &                  +CM0(2,NX,IS)         *CEY
     &                  +CM1(1,NX,IS)   *CSKXL*CEZ)
               CDX=-CI*(+CM2(1,NX,IS)*CI*CSKXL*CEX
     &                  +CM2(2,NX,IS)*CI*CSKXL*CEY)
               CJY=-CI*(-CM0(2,NX,IS)         *CEX
     &                  +CM0(3,NX,IS)         *CEY)
               CDY=-CI*(+CM1(2,NX,IS)*CI      *CEZ
     &                  -CM2(2,NX,IS)*CI*CSKXL*CEX
     &                  +CM2(3,NX,IS)*CI*CSKXL*CEY)
               CJZ=-CI*( CM0(4,NX,IS)         *CEZ
     &                  -CM1(2,NX,IS)   *CSKXL*CEY)
               CDZ=-CI*(+CM1(1,NX,IS)*CI      *CEX
     &                  +CM2(4,NX,IS)*CI*CSKXL*CEZ)
               CVX=-CI*( CM0(1,NX,IS)            *CEX
     &                  +CM0(2,NX,IS)            *CEY
     &                  +CM1(1,NX,IS)*CSKXL      *CEZ
     &                  +CM2(1,NX,IS)*CSKXL*CSKXL*CEX
     &                  +CM2(2,NX,IS)*CSKXL*CSKXL*CEY)
               CVZ=-CI*( CM0(4,NX,IS)            *CEZ
     &                  +CM1(1,NX,IS)*CSKXL      *CEX
     &                  -CM1(2,NX,IS)*CSKXL      *CEY
     &                  +CM2(4,NX,IS)*CSKXL*CSKXL*CEZ)
C
               CJF(NX,IS,1)=CJF(NX,IS,1)+CJX
               CJF(NX,IS,2)=CJF(NX,IS,2)+CJY
               CJF(NX,IS,3)=CJF(NX,IS,3)+CJZ
               CJF(NX,IS,4)=CJF(NX,IS,4)+CVX
               CJF(NX,IS,5)=CJF(NX,IS,5)+CI*CSKXL*CVX
               CJF(NX,IS,6)=CJF(NX,IS,6)+CVZ
               CJF(NX,IS,7)=CJF(NX,IS,7)+CI*CSKXL*CVZ
               CJD(NX,IS,1)=CJD(NX,IS,1)+CDX
               CJD(NX,IS,2)=CJD(NX,IS,2)+CDY
               CJD(NX,IS,3)=CJD(NX,IS,3)+CDZ
  100 CONTINUE
C
      DO 200 NX=1,NXP
         RHL(NX,1)=DCONJG(CAF(NX,1))*CBF(NX,1)
     &            +DCONJG(CAF(NX,2))*CBF(NX,2)
     &            +DCONJG(CAF(NX,3))*CBF(NX,3)
         RHL(NX,2)=DCONJG(CAF(NX,2))*CAF(NX,3)*(-CI)
     &            -DCONJG(CAF(NX,3))*CAF(NX,2)*(-CI)
     &            +2.D0*DCONJG(CAF(NX,4))*CBF(NX,1)
         RHL(NX,3)=DCONJG(CEF(NX,1))*CBF(NX,1)
     &            +DCONJG(CEF(NX,2))*CBF(NX,2)
     &            +DCONJG(CEF(NX,3))*CBF(NX,3)
C
         RHL(NX,4)=0.D0
         RHL(NX,5)=0.D0
C
         DO 200 IS=1,ISMAX
            FWP= 1.D20*AEE*AEE*PZ(IS)*PZ(IS)/(AMM*PA(IS)*EPS0*RW*RW)
            IF(PROFPN(NX,IS).LE.0.D0) THEN
               WP2I=1.D0
            ELSE
               WP2I=1.D0/(FWP*PROFPN(NX,IS))
            ENDIF
            FHL(NX,IS,1)=RNPR*(DCONJG(CEF(NX,1))*CJF(NX,IS,1)
     &                        +DCONJG(CEF(NX,2))*CJF(NX,IS,2)
     &                        +DCONJG(CEF(NX,3))*CJF(NX,IS,3)
     &                        +DCONJG(CED(NX,1))*CJD(NX,IS,1)
     &                        +DCONJG(CED(NX,2))*CJD(NX,IS,2)
     &                        +DCONJG(CED(NX,3))*CJD(NX,IS,3))
            FHL(NX,IS,2)=CI*(DCONJG(CJF(NX,IS,4))*CED(NX,3)
     &                      +DCONJG(CJF(NX,IS,5))*CEF(NX,3))
CXX  &             +CI*(WP2I*DCONJG(CJF(NX,IS,4))*CJF(NX,IS,7)
CXX  &                 +WP2I*DCONJG(CJF(NX,IS,5))*CJF(NX,IS,6))
C
            RHL(NX,4)=RHL(NX,4)+       FHL(NX,IS,2)
            RHL(NX,5)=RHL(NX,5)-PZ(IS)*FHL(NX,IS,2)
  200 CONTINUE
C
      TAU=(4.D0*PI*(EPS0/AEE)**2*(AME/AEE)**2)/(1.D20)
      VPH=2.D0*PI*RF*1.D6/RKZ
      FACTC=3.D0*SQRT(0.5D0*PI)/(ZEFF*(0.29+0.46/(1.08+ZEFF)))
      AWEC=2.D0*PI*RF*1.D6*EPS0*AEE/(VC*AMM)
C
      DO 300 NX=1,NXP
         XR=XAM(NX)/RR
         YR=0.1D0/RR
         EP=SQRT(XR*XR+YR*YR)
         RLNLMD=16.1D0 - 1.15D0*LOG10(PROFPN(NX,1))
     &                 + 2.30D0*LOG10(PROFTR(NX,1))
         VTE=SQRT(PROFTR(NX,1)*AEE*1.D3/AME)
         TAU0=TAU*VTE**3/(PROFPN(NX,1)*RLNLMD)
         FRAC=1.D0-1.9D0*SQRT(EP)+ABS(EP)
 
            IS=1
            W=VPH/VTE
            TAU1=TAU0*W*W1CDEF(W,ZEFF,0.D0,0.D0,ICDTYP)
            TAU2=TAU0*FACTC
            TAU3=TAU0*W*W1CDEF(W,ZEFF,XR,YR,ICDTYP)
            TAU4=TAU0*FACTC*FRAC
            AHL(NX,IS,1)=-AWEC*(PZ(IS)/PA(IS))*FHL(NX,IS,1)*TAU1
            AHL(NX,IS,2)=-AWEC*(PZ(IS)/PA(IS))*FHL(NX,IS,2)*TAU2
            AHL(NX,IS,3)=-AWEC*(PZ(IS)/PA(IS))*FHL(NX,IS,1)*TAU3
            AHL(NX,IS,4)=-AWEC*(PZ(IS)/PA(IS))*FHL(NX,IS,2)*TAU4
C
         DO 300 IS=2,ISMAX
            TAUE=TAU0*FACTC*SQRT(PROFTR(NX,IS)/PROFTR(NX,1))**3
     &          *SQRT(2.D0*PA(IS)*AMM/AME)/(ZEFF*PZ(IS)**2) 
            TAUI=TAU0*FACTC*PA(IS)*AMM/(AME*PZ(IS)**2)
            TAU1=(1.D0-PZ(IS)/ZEFF)/(1.D0/TAUE+1.D0/TAUI)
            TAU2=TAU1
            TAU3=(1.D0-PZ(IS)/ZEFF)/(1.D0/TAUE+1.D0/TAUI)
            TAU4=TAU1
            AHL(NX,IS,1)=-AWEC*(PZ(IS)/PA(IS))*FHL(NX,IS,1)*TAU1
            AHL(NX,IS,2)=-AWEC*(PZ(IS)/PA(IS))*FHL(NX,IS,2)*TAU2
            AHL(NX,IS,3)=-AWEC*(PZ(IS)/PA(IS))*FHL(NX,IS,1)*TAU3
            AHL(NX,IS,4)=-AWEC*(PZ(IS)/PA(IS))*FHL(NX,IS,2)*TAU4
  300 CONTINUE
C
      DO 400 IS=1,ISMAX
         AHLT(IS,1)=0.D0
         AHLT(IS,2)=0.D0
         AHLT(IS,3)=0.D0
         AHLT(IS,4)=0.D0
      DO 400 NX=1,NXP
         AHLT(IS,1)=AHLT(IS,1)+AHL(NX,IS,1)*(XA(NX+1)-XA(NX))
         AHLT(IS,2)=AHLT(IS,2)+AHL(NX,IS,2)*(XA(NX+1)-XA(NX))
         AHLT(IS,3)=AHLT(IS,3)+AHL(NX,IS,3)*(XA(NX+1)-XA(NX))
         AHLT(IS,4)=AHLT(IS,4)+AHL(NX,IS,4)*(XA(NX+1)-XA(NX))
  400 CONTINUE
C
C     RHL(,6)=EPARA=HDOT*C/(B0*C)    V/M FOR MW
C     RHL(,7)=EPARA=FE2*NE*ME/QE     V/M FOR MW
C
      A6=1.D6/(BB*VC)
      A7=1.D6*2.D0*PI*RF*1.D6*EPS0/(VC*1.D20*AEE)
      DO 500 NX=1,NXP
         RHL(NX,6)=A6*RHL(NX,3)/PROFB(NX)
         RHL(NX,7)=A7*FHL(NX,1,2)/PROFPN(NX,1)
  500 CONTINUE
C
      RETURN
      END
