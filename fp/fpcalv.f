C     $Id$
C
C ************************************************************
C
C            CALCULATION OF DW (BOUNCE AVERAGED)
C
C ************************************************************
C
      SUBROUTINE FPCALV
C
      INCLUDE 'fpcomm.inc'
C
C =============  CALCULATION OF DWPP AND DWPT  ===============
C
      FACT=0.5D0
C
      DO 1000 NRDO=1,NRMAX
         NR=NRDO
C
         DO 101 NTH=1,NTHMAX
            IF(NTH.EQ.ITL(NR).OR.NTH.EQ.ITU(NR)) GOTO 101
         DO 100 NP=1,NPMAX+1
            CALL FPSUMV(ETAM(NTH,NR),SINM(NTH),COSM(NTH),PG(NP),NR,
     &                 DWPPS,DWPTS,DWTPS,DWTTS)
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
                  CALL FPSUMV(ETAG(NTH,NR),SING(NTH),COSG(NTH),PM(NP),
     &                        NR,DWPPS,DWPTS,DWTPS,DWTTS)
                  DWTP(NTH,NP,NR)=DWTPS
                  DWTT(NTH,NP,NR)=DWTTS
 1100          CONTINUE
            ELSE
               PSIL=RM(NR)**2
               DO 1200 NP=1,NPMAX
                  P=PM(NP)
                  CALL FPWAVV(PSIL,ETAG(NTH,NR),SING(NTH),COSG(NTH),P,
     &                        DWPPL,DWPTL,DWTPL,DWTTL)
                  DWTP(NTH,NP,NR)=0.D0
                  DWTT(NTH,NP,NR)=DWTTL
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
C =======================================================
C
      SUBROUTINE FPSUMV(ETA,RSIN,RCOS,P,NR,DWPPS,DWPTS,DWTPS,DWTTS)
C
      INCLUDE 'fpcomm.inc'
C
      DELH=2.D0*ETA/NAVMAX
C
      DWPPS=0.D0
      DWPTS=0.D0
      DWTPS=0.D0
      DWTTS=0.D0
      PSIL=RM(NR)**2
C      
      CALL PLBMIN(PSIL,BMIN)
C
      DO 100 N=1,NAVMAX
         ETAL=DELH*(N-0.5D0)
         CALL FPSETB(PSIL,ETAL,B0TH,B0PH)
C
         IF(MODELA.EQ.0) THEN
            PSI=1.D0
            PSIN=RSIN
            PCOS=RCOS
         ELSE
            PSI=SQRT((B0TH**2+B0PH**2)/BMIN**2)
            PSIN=SQRT(PSI)*RSIN
            IF (RCOS.GT.0.0D0) THEN
               PCOS= SQRT(1.D0-PSI*RSIN**2)
            ELSE
               PCOS=-SQRT(1.D0-PSI*RSIN**2)
            END IF
         ENDIF
         CALL FPWAVV(PSIL,ETAL,PSIN,PCOS,P,DWPPL,DWPTL,DWTPL,DWTTL)
C
C         RVPARA=RCOS/PCOS
C         TAN0=RSIN/RCOS
C         TANL=PSIN/PCOS
C         DWPPS=DWPPS+DWPPL*RVPARA
C         DWPTS=DWPTS+DWPTL*RVPARA*TAN0/TANL
C         DWTPS=DWTPS+DWTPL*RVPARA*TAN0/TANL
C         DWTTS=DWTTS+DWTTL*RVPARA*TAN0**2/(TANL**2)
C
         DWPPS=DWPPS+DWPPL*RCOS/PCOS
         DWPTS=DWPTS+DWPTL          /SQRT(PSI)
         DWTPS=DWTPS+DWTPL          /SQRT(PSI)
         DWTTS=DWTTS+DWTTL*PCOS/RCOS/PSI
  100 CONTINUE
         DWPPS=DWPPS*DELH/PI
         DWPTS=DWPTS*DELH/PI
         DWTPS=DWTPS*DELH/PI
         DWTTS=DWTTS*DELH/PI
      RETURN
      END
C
C ************************************************************
C
      SUBROUTINE FPWAVV(PSIL,ETAL,PSIN,PCOS,P,DWPPL,DWPTL,DWTPL,DWTTL)
C
      INCLUDE 'fpcomm.inc'
C
      PARAMETER(NJMAX=100)
      DIMENSION RJ(0:NJMAX),DRJ(0:NJMAX)
      DATA CI/(0.D0,1.D0)/   
C
      THETA0=RTFP0*1.D3*AEE/(AMFP*VC*VC)
      CALL FPSETV(PSIL,ETAL,RKR,RKTH,RKPH,CER,CETH,CEPH)
      CALL FPSETB(PSIL,ETAL,B0TH,B0PH)
C
      RW     =2.D0*PI*RFDW*1.D6
      B2     =B0TH**2+B0PH**2
      CEPARA =(B0TH*CETH+B0PH*CEPH)/SQRT(B2)
      CEPERP =(B0PH*CETH-B0TH*CEPH)/SQRT(B2)
      CEPLUS =(CER+CI*CEPERP)/SQRT(2.D0)
      CEMINUS=(CER-CI*CEPERP)/SQRT(2.D0)
      RGAMMA =SQRT(1.D0+P*P*THETA0)
      VPARA  =PTFP0*P*PCOS/(AMFP*RGAMMA)
      VPERP  =PTFP0*P*PSIN/(AMFP*RGAMMA)
      RKPARA =(B0TH*RKTH+B0PH*RKPH)/SQRT(B2)
      RKPERP =(B0PH*RKTH-B0TH*RKPH)/SQRT(B2)
      RWC    =AEFP*SQRT(B2)/AMFP
C 
      DWC11=0.D0
      DWC12=0.D0
      DWC21=0.D0
      DWC22=0.D0
      RKW  =RKPARA/RW
      RGZAI=RKPERP*VPERP/ABS(RWC)
      NHMAX=MAX(ABS(NCMIN-1),ABS(NCMAX+1),2)
      CALL BESSJN(RGZAI,NHMAX,RJ,DRJ)
C      CALL BESSEL(RGZAI,RJ,NCBMAX,NJMAX+1)
C
      DO 100 N=NCMIN,NCMAX
         IF (N.LT.0) THEN
             RJN=(-1)**(-N)*RJ(-N)
         ELSE
             RJN=RJ(N)
         ENDIF
         NMI=N-1
         IF (NMI.LT.0) THEN
             RJNM=(-1)**(-NMI)*RJ(-NMI)
         ELSE
             RJNM=RJ(NMI)
         ENDIF
         NPI=N+1
         IF (NPI.LT.0) THEN
             RJNP=(-1)**(-NPI)*RJ(-NPI)
         ELSE
             RJNP=RJ(NPI)
         ENDIF
         IF (N.EQ.0) THEN
             CTHETA=VPERP*(CEPLUS*RJNM+CEMINUS*RJNP)/SQRT(2.D0)
     &              +CEPARA*VPARA*RJN
         ELSE
             CTHETA=(CEPLUS*RJNM+CEMINUS*RJNP)/SQRT(2.D0)
     &              +CEPARA*VPARA*(RJNM+RJNP)*RKPERP
     &              /(2*N*RWC)
         ENDIF
         RTHETA2=ABS(CTHETA)**2
         IF (N.EQ.0) THEN
            A11=0
            A12=0
            A21=0
            A22=RTHETA2*RKW**2
	 ELSE
            A11=RTHETA2*(1.D0-RKW*VPARA)**2
            A12=RTHETA2*RKW*VPERP*(1.D0-RKW*VPARA)
            A21=RTHETA2*RKW*VPERP*(1.D0-RKW*VPARA)
            A22=RTHETA2*RKW**2*VPERP**2
	 ENDIF        
         IF(VPARA.EQ.0.D0) THEN
            DWC=0.D0  
         ELSE
            EX=-((RW-RKPARA*VPARA-N*RWC/RGAMMA)
     &           /(RW*VPARA*DELNPR/VC))**2
            IF (EX.LT.-100.D0) THEN 
                DWC=0.D0
            ELSE
                DWC=0.5D0*SQRT(PI)*AEFP**2*EXP(EX)/PTFP0**2
     &              /(RW*ABS(VPARA)*DELNPR/VC)
            ENDIF
         ENDIF
         DWC11=DWC11+DWC*A11
         DWC12=DWC12+DWC*A12
         DWC21=DWC21+DWC*A21
         DWC22=DWC22+DWC*A22
  100 CONTINUE
C
      DWPPL=PSIN**2*DWC11  +PSIN*PCOS*(DWC12+DWC21)    +PCOS**2*DWC22
      DWPTL=PSIN*PCOS*DWC11-PSIN**2*DWC12+PCOS**2*DWC21-PSIN*PCOS*DWC22
      DWTPL=PSIN*PCOS*DWC11+PCOS**2*DWC12-PSIN**2*DWC21-PSIN*PCOS*DWC22
      DWTTL=PCOS**2*DWC11  -PSIN*PCOS*(DWC12+DWC21)    +PSIN**2*DWC22
C
      RETURN
      END
C
C     ****** CALCULATE LOCAL MAGNETIC FIELD ******
C
      SUBROUTINE FPSETB(PSIL,ETAL,BTL,BPL)
C
      INCLUDE 'fpcomm.inc'
C
      IF(MODELG.EQ.2) THEN
         RS=SQRT(PSIL)*RA
         X=RS*COS(ETAL)
         BTL=BB/(1.D0+X/RR)
         CALL PLQPRF(PSIL,QL)
         BPL=RS*BTL/((RR+X)*QL)
      ELSEIF(MODELG.EQ.3) THEN
         CALL PLRRMX(PSIL,RRMINL,RRMAXL)
         RS1=RRMAXL-RR
         RS2=RR-RRMINL
         RS=0.5D0*(RS1+RS2)+0.5D0*(RS1-RS2)*COS(ETAL)
         X=RS*COS(ETAL)
         CALL PLQPRF(PSIL,QL)
         BTL=BB/(1.D0+X/RR)
         BPL=RS*BTL/((RR+X)*QL)
      ENDIF
      RETURN
      END
C
C     ****** CALCULATE LOCAL WAVE ELECTRIC FIELD ******
C
      SUBROUTINE FPSETV(RHO0,ETAL,RKR,RKTH,RKPH,CER,CETH,CEPH)
C
      INCLUDE 'fpcomm.inc'
C
      Y=RHO0*RA*SIN(ETAL)
C
      ARG=(Y-REWY)**2/DREWY**2
      IF(ARG.GT.100.D0) THEN
         FACT=0.D0
      ELSE
         FACT=EXP(-ARG)
      ENDIF
C
      CER= CEWR*FACT
      CETH=CEWTH*FACT
      CEPH=CEWPH*FACT
C
      RKR= RKWR
      RKTH=RKWTH
      RKPH=RKWPH
C
      RETURN
      END

