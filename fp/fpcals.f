C     $Id$
C
C ************************************************************
C
C            CALCULATION OF DW (BOUNCE AVERAGED,RAY)
C
C ************************************************************
C
      SUBROUTINE FPCALS
C
      INCLUDE 'fpcomm.h'
C
C =============  CALCULATION OF DWPP AND DWPT  ===============
C
      FACT=0.5D0
C
      DO 1000 NRDO=1,NRMAX
         NR=NRDO
         DO 101 NTH=1,NTHMAX
            IF(NTH.EQ.ITL(NR).OR.NTH.EQ.ITU(NR)) GOTO 101
         DO 100 NP=1,NPMAX+1
            IF(NP.EQ.10.AND.NTH.EQ.2) THEN
               IFLAG=1
            ELSE
               IFLAG=0
            ENDIF
            CALL FPDWAV(ETAM(NTH,NR),SINM(NTH),COSM(NTH),PG(NP),NR,
     &                  DWPPS,DWPTS,DWTPS,DWTTS)
            DWPP(NTH,NP,NR)=DWPPS
            DWPT(NTH,NP,NR)=DWPTS
  100    CONTINUE
  101    CONTINUE
C
         IF(MODELA.EQ.1) THEN
         DO 200 NP=1,NPMAX+1
C
         DO 190 NTH=ITL(NR)+1,NTHMAX/2
            DWPP(NTH,NP,NR)  =(DWPP(NTH,NP,NR)
     &                        +DWPP(NTHMAX-NTH+1,NP,NR))*FACT
            DWPT(NTH,NP,NR)  =(DWPT(NTH,NP,NR)
     &                        +DWPT(NTHMAX-NTH+1,NP,NR))*FACT
            DWPP(NTHMAX-NTH+1,NP,NR)  =DWPP(NTH,NP,NR)
            DWPT(NTHMAX-NTH+1,NP,NR)  =DWPT(NTH,NP,NR)
  190    CONTINUE
         DWPP(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                    +DWPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                    +DWPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                    +DWPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
C
         DWPT(ITL(NR),NP,NR)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWPT(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                    +DWPT(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                    +DWPT(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                    +DWPT(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
         DWPP(ITU(NR),NP,NR)  =DWPP(ITL(NR),NP,NR)
         DWPT(ITU(NR),NP,NR)  =DWPT(ITL(NR),NP,NR)
  200    CONTINUE
         ENDIF
 1000 CONTINUE
C
C =============  CALCULATION OF DWTP AND DWTT  ===============
C
      DO 2000 NRDO=1,NRMAX
         NR=NRDO
C
         DO 1101 NTH=1,NTHMAX+1
            IF(NTH.NE.NTHMAX/2+1) THEN
               DO 1100 NP=1,NPMAX
                  CALL FPDWAV(ETAG(NTH,NR),SING(NTH),COSG(NTH),PM(NP),
     &                        NR,DWPPS,DWPTS,DWTPS,DWTTS)
                  DWTP(NTH,NP,NR)=DWTPS
                  DWTT(NTH,NP,NR)=DWTTS
 1100          CONTINUE
            ELSE
               DO 1200 NP=1,NPMAX
C                  CALL FPDWAV(ETAG(NTH,NR),SING(NTH),COSG(NTH),PM(NP),
C     &                        NR,DWPPS,DWPTS,DWTPS,DWTTS)
                  DWTP(NTH,NP,NR)=0.D0
                  DWTT(NTH,NP,NR)=0.D0
 1200          CONTINUE
            ENDIF
 1101    CONTINUE
C
         IF(MODELA.EQ.1) THEN
            DO 1300 NTH=ITL(NR)+1,NTHMAX/2
            DO 1300 NP=1,NPMAX
               DWTP(NTH,NP,NR)=(DWTP(NTH,NP,NR)
     &                         +DWTP(NTHMAX-NTH+2,NP,NR))*FACT
               DWTT(NTH,NP,NR)=(DWTT(NTH,NP,NR)
     &                         +DWTT(NTHMAX-NTH+2,NP,NR))*FACT
               DWTP(NTHMAX-NTH+2,NP,NR)=DWTP(NTH,NP,NR)
               DWTT(NTHMAX-NTH+2,NP,NR)=DWTT(NTH,NP,NR)
 1300       CONTINUE
         ENDIF
 2000 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C     Calculate DW averaged over magnetic surface for (p,theta0,r)
C***********************************************************************
C
      SUBROUTINE FPDWAV(ETA,RSIN,RCOS,P,NR,DWPPS,DWPTS,DWTPS,DWTTS)
C
      INCLUDE 'fpcomm.h'
      INCLUDE '../wr/wrcom1.h'
      INCLUDE 'fpcom2.h'
C
      DELH=4.D0*ETA/NAVMAX
C
      DWPPS=0.D0
      DWPTS=0.D0
      DWTPS=0.D0
      DWTTS=0.D0
C
      DO NAV=1,NAVMAX
         ETAL=DELH*(NAV-0.5D0)-2.D0*ETA
C
         THETAL=ATAN2(RKAP*SIN(ETAL),COS(ETAL))
         RRAVE=0.5D0*(RRMAX(NR)+RRMIN(NR))
         RSAVE=0.5D0*(RRMAX(NR)-RRMIN(NR))
         X=RRAVE+RSAVE*COS(THETAL)
         Y=0.D0
         Z=      RSAVE*SIN(THETAL)*RKAP
C
         RL=SQRT(X**2+Y**2)
         ZL=Z
         DO NRAY=1,NRAYMX
C            IF(IFLAG.EQ.1) WRITE(6,*) 'NR,NCRMAX=',NR,NCRMAX(NR,NRAY)
         DO NCR=1,NCRMAX(NR,NRAY)
            RFDW=RAYIN(1,NRAY)
            RX=RCR(1,NCR,NR,NRAY)
            RY=RCR(2,NCR,NR,NRAY)
            RZ=RCR(3,NCR,NR,NRAY)
            RLCR=SQRT(RX**2+RY**2)
            ZLCR=RZ
            DELR2=(RL-RLCR)**2+(ZL-ZLCR)**2
            DELCR2=DELYEC**2
            ARG=DELR2/DELCR2
C            IF(IFLAG.EQ.1) THEN
C               WRITE(6,'(3I3)') NR,NAV,NCR
C               WRITE(6,'(1P3E12.4)') X,Y,Z
C               WRITE(6,'(1P3E12.4)') RX,RY,RZ
C               WRITE(6,'(1P3E12.4)') DELR2,DELCR2,ARG
C            ENDIF
            IF(ARG.LT.15.D0) THEN
               FACTOR=EXP(-ARG)
C               WRITE(25,*) 'NAV,NCR,DELR2,ARG=',NAV,NCR,DELR2,ARG
               CALL FPDWRP(NR,ETAL,RSIN,RCOS,PSIN,PCOS,PSI)
               CEX=CECR(1,NCR,NR,NRAY)*FACTOR
               CEY=CECR(2,NCR,NR,NRAY)*FACTOR
               CEZ=CECR(3,NCR,NR,NRAY)*FACTOR
               RKX=RKCR(1,NCR,NR,NRAY)
               RKY=RKCR(2,NCR,NR,NRAY)
               RKZ=RKCR(3,NCR,NR,NRAY)
               CALL FPDWLL(P,PSIN,PCOS,
     &                     CEX,CEY,CEZ,RKX,RKY,RKZ,RX,RY,RZ,
     &                     DWPPL,DWPTL,DWTPL,DWTTL)
               DWPPS=DWPPS+DWPPL*RCOS/PCOS
               DWPTS=DWPTS+DWPTL          /SQRT(PSI)
               DWTPS=DWTPS+DWTPL          /SQRT(PSI)
               DWTTS=DWTTS+DWTTL*PCOS/RCOS/PSI
            ENDIF
         ENDDO
         ENDDO
      ENDDO
C
      FACTOR=PWAVE*1.D6/(VC*EPS0)
     &      /(2.D0*PI*RR)
     &      /SQRT(PI*DELYEC**2/2)
     &      *DELH/(2.D0*PI)
      DWPPS=DWPPS*FACTOR
      DWPTS=DWPTS*FACTOR
      DWTPS=DWTPS*FACTOR
      DWTTS=DWTTS*FACTOR
C
      RETURN
      END
C
C***********************************************************************
C     Calculate PSIN, PCOS, PSI
C***********************************************************************
C
      SUBROUTINE FPDWRP(NR,ETAL,RSIN,RCOS,PSIN,PCOS,PSI)
C
      INCLUDE 'fpcomm.h'
C
      IF(MODELA.EQ.0) THEN
         PSI=1.D0
         PSIN=RSIN
         PCOS=RCOS
      ELSE
         PSI=(1.D0+EPSR(NR))/(1.D0+EPSR(NR)*COS(ETAL))
         PSIN=SQRT(PSI)*RSIN
         IF (RCOS.GT.0.0D0) THEN
            PCOS= SQRT(1.D0-PSI*RSIN**2)
         ELSE
            PCOS=-SQRT(1.D0-PSI*RSIN**2)
         END IF
      ENDIF
      RETURN
      END
C
C***********************************************************************
C     calculate local diffusion coefficient
C***********************************************************************
C  
      SUBROUTINE FPDWLL(P,PSIN,PCOS,
     &                  CEX,CEY,CEZ,RKX,RKY,RKZ,RX,RY,RZ,
     &                  DWPPL,DWPTL,DWTPL,DWTTL)
C
      INCLUDE 'fpcomm.h'
C
      PARAMETER(NJMAX=100)
      DIMENSION RJ(0:NJMAX),DRJ(0:NJMAX)
      DATA CI/(0.D0,1.D0)/   
C
      CALL PLMAG(RX,RY,RZ,PSIL)
C
      RW     =2.D0*PI*RFDW*1.D6
      RWC    =-AEE*BL/AME
C
      RKPARA = RKX*BX +RKY*BY +RKZ*BZ
      RKPERP = SQRT(RKX*RKX+RKY*RKY+RKZ*RKZ-RKPARA*RKPARA)
C
      U1X = (RKX-BX*RKPARA)/RKPERP
      U1Y = (RKY-BY*RKPARA)/RKPERP
      U1Z = (RKZ-BZ*RKPARA)/RKPERP
C
      U2X = (BY*RKZ-BZ*RKY)/RKPERP
      U2Y = (BZ*RKX-BX*RKZ)/RKPERP
      U2Z = (BX*RKY-BY*RKX)/RKPERP
C 
      CE1    = CEX*U1X+CEY*U1Y+CEZ*U1Z
      CE2    = CEX*U2X+CEY*U2Y+CEZ*U2Z
      CEPARA = CEX*BX +CEY*BY +CEZ*BZ
C
      CEPLUS =(CE1+CI*CE2)/SQRT(2.D0)
      CEMINUS=(CE1-CI*CE2)/SQRT(2.D0)
C
      RGAMMA =SQRT(1.D0+P*P*THETA0)
      PPARA  =PTH0*P*PCOS
      PPERP  =PTH0*P*PSIN
      VPARA  =PPARA/(AME*RGAMMA)
      VPERP  =PPERP/(AME*RGAMMA)
C 
      DWC11=0.D0
      DWC12=0.D0
      DWC21=0.D0
      DWC22=0.D0
      RKW  =RKPARA/RW
      RGZAI=RKPERP*PPERP/(RWC*AME)
C
C      WRITE(6,*) 'RKPERP,PPERP,RWC,RGZAI = ',RKPERP,PPERP,RWC,RGZAI
C      CALL BESSEL(RGZAI,RJ,NCBMAX,NJMAX+1)
C
      NHMAX=MAX(ABS(NCMIN-1),ABS(NCMAX+1),2)
      CALL BESSJN(RGZAI,NHMAX,RJ,DRJ)
C
      DO 100 NC=NCMIN,NCMAX
         
         IF (NC.LT.0) THEN
             RJN=(-1)**(-NC)*RJ(-NC)
         ELSE
             RJN=RJ(NC)
         ENDIF
         NMI=NC-1
         IF (NMI.LT.0) THEN
             RJNM=(-1)**(-NMI)*RJ(-NMI)
         ELSE
             RJNM=RJ(NMI)
         ENDIF
         NPI=NC+1
         IF (NPI.LT.0) THEN
             RJNP=(-1)**(-NPI)*RJ(-NPI)
         ELSE
             RJNP=RJ(NPI)
         ENDIF
         IF (NC.EQ.0) THEN
             CTHETA=PPERP*(CEPLUS*RJNM+CEMINUS*RJNP)/SQRT(2.D0)
     &              +CEPARA*PPARA*RJN
         ELSE
             CTHETA=(CEPLUS*RJNM+CEMINUS*RJNP)/SQRT(2.D0)
     &              +CEPARA*PPARA*(RJNM+RJNP)*RKPERP
     &              /(2*NC*RWC*AME)
         ENDIF
         RTHETA2=ABS(CTHETA)**2
         IF (NC.EQ.0) THEN
            A11=0
            A12=0
            A21=0
            A22=RTHETA2*RKW**2/(AME**2*RGAMMA**2)
	 ELSE
            A11=RTHETA2*(1.D0-RKW*VPARA)**2
            A12=RTHETA2*RKW*VPERP*(1.D0-RKW*VPARA)
            A21=RTHETA2*RKW*VPERP*(1.D0-RKW*VPARA)
            A22=RTHETA2*RKW**2*VPERP**2
	 ENDIF        
         IF(VPARA.EQ.0.D0) THEN
            DWC=0.D0  
         ELSE
            EX=-((RGAMMA-RKPARA*PPARA/(RW*AME)-NC*RWC/RW)
     &           /(PPARA*DELNPR/(AME*VC)))**2
            IF (EX.LT.-100.D0) THEN 
                DWC=0.D0
            ELSE
                DWC=0.5D0*SQRT(PI)*AEE**2*EXP(EX)/PTH0**2
     &              /(RW*ABS(PPARA)*DELNPR/(AME*VC))
            ENDIF
         ENDIF
         DWC11=DWC11+DWC*A11
         DWC12=DWC12+DWC*A12
         DWC21=DWC21+DWC*A21
         DWC22=DWC22+DWC*A22
  100 CONTINUE
C
      DWPPL=PSIN*(PSIN*DWC11+PCOS*DWC12)
     &     +PCOS*(PSIN*DWC21+PCOS*DWC22)
      DWPTL=PSIN*(PCOS*DWC11-PSIN*DWC12)
     &     +PCOS*(PCOS*DWC21-PSIN*DWC22)
      DWTPL=PCOS*(PSIN*DWC11+PCOS*DWC12)
     &     -PSIN*(PSIN*DWC21+PCOS*DWC22)
      DWTTL=PCOS*(PCOS*DWC11-PSIN*DWC12)
     &     -PSIN*(PCOS*DWC21-PSIN*DWC22)
C
      RETURN
      END
