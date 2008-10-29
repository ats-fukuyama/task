C     $Id$
C
C ************************************************************
C
C            CALCULATION OF DW (BOUNCE AVERAGED)
C
C ************************************************************
C
      SUBROUTINE FPCALWM(NSA)
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
     &                 DWPPS,DWPTS,DWTPS,DWTTS,NSA)
            DWPP(NTH,NP,NR,NSA)=DWPPS
            DWPT(NTH,NP,NR,NSA)=DWPTS
  100    CONTINUE
  101    CONTINUE
C
         IF(MODELA.EQ.1) THEN
         DO 200 NP=1,NPMAX+1
C
         DO 190 NTH=ITL(NR)+1,NTHMAX/2
            DWPP(NTH,NP,NR,NSA)  =(DWPP(NTH,NP,NR,NSA)
     &                        +DWPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
            DWPT(NTH,NP,NR,NSA)  =(DWPT(NTH,NP,NR,NSA)
     &                        +DWPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
            DWPP(NTHMAX-NTH+1,NP,NR,NSA)  =DWPP(NTH,NP,NR,NSA)
            DWPT(NTHMAX-NTH+1,NP,NR,NSA)  =DWPT(NTH,NP,NR,NSA)
  190    CONTINUE
         DWPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWPP(ITL(NR)-1,NP,NR,NSA)
     &                                   /RLAMDA(ITL(NR)-1,NR)
     &                    +DWPP(ITL(NR)+1,NP,NR,NSA)
     &                                   /RLAMDA(ITL(NR)+1,NR)
     &                    +DWPP(ITU(NR)-1,NP,NR,NSA)
     &                                   /RLAMDA(ITU(NR)-1,NR)
     &                    +DWPP(ITU(NR)+1,NP,NR,NSA)
     &                                   /RLAMDA(ITU(NR)+1,NR))
C
         DWPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWPT(ITL(NR)-1,NP,NR,NSA)
     &                                   /RLAMDA(ITL(NR)-1,NR)
     &                    +DWPT(ITL(NR)+1,NP,NR,NSA)
     &                                   /RLAMDA(ITL(NR)+1,NR)
     &                    +DWPT(ITU(NR)-1,NP,NR,NSA)
     &                                   /RLAMDA(ITU(NR)-1,NR)
     &                    +DWPT(ITU(NR)+1,NP,NR,NSA)
     &                                   /RLAMDA(ITU(NR)+1,NR))
         DWPP(ITU(NR),NP,NR,NSA)  =DWPP(ITL(NR),NP,NR,NSA)
         DWPT(ITU(NR),NP,NR,NSA)  =DWPT(ITL(NR),NP,NR,NSA)
  200    CONTINUE
         ENDIF
c         write(*,*)"test",DWPP(2,2,1,NSA),NSA
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
     &                        NR,DWPPS,DWPTS,DWTPS,DWTTS,NSA)
                  DWTP(NTH,NP,NR,NSA)=DWTPS
                  DWTT(NTH,NP,NR,NSA)=DWTTS
 1100          CONTINUE
            ELSE
               RHOL=RM(NR)
               DO 1200 NP=1,NPMAX
                  P=PM(NP)
                  CALL FPWAVV(RHOL,ETAG(NTH,NR),SING(NTH),COSG(NTH),P,
     &                        DWPPL,DWPTL,DWTPL,DWTTL,NSA)
                  DWTP(NTH,NP,NR,NSA)=0.D0
                  DWTT(NTH,NP,NR,NSA)=DWTTL
 1200          CONTINUE
            ENDIF
 1101    CONTINUE
C
         IF(MODELA.EQ.1) THEN
            DO 1300 NTH=ITL(NR)+1,NTHMAX/2
            DO 1300 NP=1,NPMAX
               DWTP(NTH,NP,NR,NSA)=(DWTP(NTH,NP,NR,NSA)
     &                         +DWTP(NTHMAX-NTH+2,NP,NR,NSA))*FACT
               DWTT(NTH,NP,NR,NSA)=(DWTT(NTH,NP,NR,NSA)
     &                         +DWTT(NTHMAX-NTH+2,NP,NR,NSA))*FACT
               DWTP(NTHMAX-NTH+2,NP,NR,NSA)=DWTP(NTH,NP,NR,NSA)
               DWTT(NTHMAX-NTH+2,NP,NR,NSA)=DWTT(NTH,NP,NR,NSA)
 1300       CONTINUE
         ENDIF
 2000 CONTINUE
C
      RETURN
      END
C
C =======================================================
C
      SUBROUTINE FPSUMV(ETA,RSIN,RCOS,P,NR,DWPPS,DWPTS,DWTPS,DWTTS,NSA)
C
      INCLUDE 'fpcomm.inc'
C
      DELH=2.D0*ETA/NAVMAX
C
      DWPPS=0.D0
      DWPTS=0.D0
      DWTPS=0.D0
      DWTTS=0.D0
      RHOL=RM(NR)
C      
      CALL PLBMIN(RHOL,BMIN)
C
      DO 100 N=1,NAVMAX
         ETAL=DELH*(N-0.5D0)
         CALL FPDWRP(NR,ETAL,RSIN,RCOS,PSIN,PCOS,PSI,NSA)
         CALL FPWAVV(RHOL,ETAL,PSIN,PCOS,P,DWPPL,DWPTL,DWTPL,DWTTL,NSA)
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
      SUBROUTINE FPWAVV(RHOL,ETAL,PSIN,PCOS,P,
     &                  DWPPL,DWPTL,DWTPL,DWTTL,NSA)
C
      INCLUDE 'fpcomm.inc'
C
      PARAMETER(NJMAX=100)
      DIMENSION RJ(0:NJMAX),DRJ(0:NJMAX)
C
c      THETA0=RTFP0*1.D3*AEE/(AMFP*VC*VC)
      THETA0(NSA)=RTFP0(NSA)*1.D3*AEE/(AMFP(NSA)*VC*VC)

      CALL FPSETV(RHOL,ETAL,RKR,RKTH,RKPH,CER,CETH,CEPH)
      CALL FPSETB(RHOL,ETAL,B0PH,B0TH)
C     BOPH=B_\phi=B_toroidal, B0TH=B_\theta=B_poloidal
C
      RW     =2.D0*PI*RFDW*1.D6
      B2     =B0TH**2+B0PH**2
C      CEPARA =(B0TH*CETH+B0PH*CEPH)/SQRT(B2)
C      CEPERP =(B0PH*CETH-B0TH*CEPH)/SQRT(B2)
      CEPARA =CEPH
      CEPERP =CETH
      RABSE = SQRT(ABS(CEPARA)**2+ABS(CEPERP)**2)
      CEPLUS =(CER+CI*CEPERP)/SQRT(2.D0)
      CEMINUS=(CER-CI*CEPERP)/SQRT(2.D0)
      RGAMMA =SQRT(1.D0+P*P*THETA0(NSA))
      VPARA  =PTFP0(NSA)*P*PCOS/(AMFP(NSA)*RGAMMA)
      VPERP  =PTFP0(NSA)*P*PSIN/(AMFP(NSA)*RGAMMA)
      RKPARA =(B0TH*RKTH+B0PH*RKPH)/SQRT(B2)
      RKPERP =(B0PH*RKTH-B0TH*RKPH)/SQRT(B2)
      RWC    =AEFP(NSA)*SQRT(B2)/AMFP(NSA)
C 
      DWC11=0.D0
      DWC12=0.D0
      DWC21=0.D0
      DWC22=0.D0
      RKW  =RKPARA/RW
      RGZAI=RKPERP*VPERP/ABS(RWC)
      NHMAX=MAX(ABS(NCMIN-1),ABS(NCMAX+1),2)

c      write(*,*)"test1",B2,RW
c      write(*,*)"test2",RABSE
      CALL BESSJN(RGZAI,NHMAX,RJ,DRJ)
C      CALL BESSEL(RGZAI,RJ,NCBMAX,NJMAX+1)
C
      NTEST=0
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
             CTHETA=VPERP*
     &           (CEPLUS*RJNM+CEMINUS*RJNP)/SQRT(2.D0)
     &              +CEPARA*VPARA*RJN
         ELSE
             CTHETA=(CEPLUS*RJNM+CEMINUS*RJNP)/SQRT(2.D0)
     &              +CEPARA*VPARA*(RJNM+RJNP)*RKPERP
     &              /(2*N*RWC)
         ENDIF
         RTHETA2=ABS(CTHETA)**2

         IF (N.EQ.0) THEN
            NTEST=1
            IF(NTEST.eq.0)THEN
            A11=0
            A12=0
            A21=0
            A22=RTHETA2*RKW**2
            ELSE
c            A11=((1.D0-RKW*VPARA)*PSIN + RKW*VPERP*PCOS )**2*RTHETA2
c            A12=((1.D0-RKW*VPARA)*PSIN + RKW*VPERP*PCOS )
c     &         *((1.D0-RKW*VPARA)*PCOS - RKW*VPERP*PSIN )*RTHETA2
c            A21=((1.D0-RKW*VPARA)*PSIN + RKW*VPERP*PCOS )
c     &         *((1.D0-RKW*VPARA)*PCOS - RKW*VPERP*PSIN )*RTHETA2
c            A22=((1.D0-RKW*VPARA)*PCOS - RKW*VPERP*PSIN )**2*RTHETA2
            A11=( RKW*PCOS )**2  *RTHETA2
            A12=-RKW**2*PCOS*PSIN*RTHETA2
            A21=-RKW**2*PCOS*PSIN*RTHETA2
            A22=( RKW*PSIN )**2  *RTHETA2
            END IF
	 ELSE
            IF(NTEST.eq.0)THEN
            A11=RTHETA2*(1.D0-RKW*VPARA)**2
            A12=RTHETA2*RKW*VPERP*(1.D0-RKW*VPARA)
            A21=RTHETA2*RKW*VPERP*(1.D0-RKW*VPARA)
            A22=RTHETA2*RKW**2*VPERP**2
            ELSE
c            A11=((1.D0-RKW*VPARA)*PSIN + RKW*VPERP*PCOS )**2*RTHETA2
c            A12=((1.D0-RKW*VPARA)*PSIN + RKW*VPERP*PCOS )
c     &         *((1.D0-RKW*VPARA)*PCOS - RKW*VPERP*PSIN )*RTHETA2
c            A21=((1.D0-RKW*VPARA)*PSIN + RKW*VPERP*PCOS )
c     &         *((1.D0-RKW*VPARA)*PCOS - RKW*VPERP*PSIN )*RTHETA2
c            A22=((1.D0-RKW*VPARA)*PCOS - RKW*VPERP*PSIN )**2*RTHETA2
            A11=(N*RWC/RW/RGAMMA*PSIN + RKW*VPERP*PCOS )**2*RTHETA2
            A12=(N*RWC/RW/RGAMMA*PSIN + RKW*VPERP*PCOS )
     &         *(N*RWC/RW/RGAMMA*PCOS - RKW*VPERP*PSIN )*RTHETA2
            A21=(N*RWC/RW/RGAMMA*PSIN + RKW*VPERP*PCOS )
     &         *(N*RWC/RW/RGAMMA*PCOS - RKW*VPERP*PSIN )*RTHETA2
            A22=(N*RWC/RW/RGAMMA*PCOS - RKW*VPERP*PSIN )**2*RTHETA2
            END IF
         ENDIF        
         IF(VPARA.EQ.0.D0) THEN
            DWC=0.D0  
         ELSE
C            EX=-((RW-RKPARA*VPARA-N*RWC/RGAMMA)
C     &           /(RW*VPARA*DELNPR/VC))**2
            EX=-( (RW-RKPARA*VPARA-N*RWC/RGAMMA)
     &           /(DELNPR*RW))**2
            IF (EX.LT.-100.D0) THEN 
                DWC=0.D0
            ELSE
C                DWC=0.5D0*SQRT(PI)*AEFP(NSA)**2*EXP(EX)/PTFP0(NSA)**2
C     &              /(RW*ABS(VPARA)*DELNPR/VC)
                DWC=0.5D0*SQRT(PI)*AEFP(NSA)**2*EXP(EX)/PTFP0(NSA)**2
     &              /ABS(RW)/DELNPR
            ENDIF
         ENDIF
c         if(DWC.ne.0.and.ratioCE.ge.1.D0)
c     &        write(6,1333)N,P,ratioCE,ratioBE,ratioCE2
c     &        ,rtest
c     &        ,RKPARA*PTFP0/(AMFP*RGAMMA)/vc*DELT*P
c     &        ,DELNPR
c     &        ,RW
c     &        ,RKPARA*VPARA
c     &        ,RWC
c/RGAMMA
c     &        ,RW/RKPARA*(1.D0-RWC/RW)/VTFP0

c     &       ,PCOS,p
c     &       ,ex

 1333    FORMAT(I2,6E14.6)
         DWC11=DWC11+DWC*A11
         DWC12=DWC12+DWC*A12
         DWC21=DWC21+DWC*A21
         DWC22=DWC22+DWC*A22

  100 CONTINUE

      ratioCE=ABS(CEPARA)/ABS(CEPERP)
      ratioCE2=ABS(CEPH)/ABS(CETH)
      ratioBE=B0PH/B0TH
c      if(ratioCE.ge.1.D0)
c     &     write(*,1333)N, B0TH, B0PH, CETH, CEPH
cCEPARA, CEPERP

C
      IF(NTEST.eq.0)THEN
      DWPPL=PSIN**2*DWC11  +PSIN*PCOS*(DWC12+DWC21)    +PCOS**2*DWC22
      DWPTL=PSIN*PCOS*DWC11-PSIN**2*DWC12+PCOS**2*DWC21-PSIN*PCOS*DWC22
      DWTPL=PSIN*PCOS*DWC11+PCOS**2*DWC12-PSIN**2*DWC21-PSIN*PCOS*DWC22
      DWTTL=PCOS**2*DWC11  -PSIN*PCOS*(DWC12+DWC21)    +PSIN**2*DWC22
      ELSE
      DWPPL = DWC11
      DWPTL = DWC12
      DWTPL = DWC21
      DWTTL = DWC22
      END IF

C
      RETURN
      END
C
C     ****** CALCULATE LOCAL MAGNETIC FIELD ******
C
      SUBROUTINE FPSETB(RHOL,ETAL,BTL,BPL)
C
C     BTL=B_toroidal, BPL=B_poloidal
C
      INCLUDE 'fpcomm.inc'
C
      IF(MODELG.EQ.2) THEN
         RS=RSRHON(RHOL)
         X=RS*COS(ETAL)
         BTL=BB/(1.D0+X/RR)
         CALL PLQPRF(RHOL,QL)
         BPL=RS*BTL/((RR+X)*QL)
      ELSEIF(MODELG.EQ.3) THEN
         CALL PLRRMX(RHOL,RRMINL,RRMAXL)
         RS1=RRMAXL-RR
         RS2=RR-RRMINL
         RS=0.5D0*(RS1+RS2)+0.5D0*(RS1-RS2)*COS(ETAL)
         X=RS*COS(ETAL)
         CALL PLQPRF(RHOL,QL)
         BTL=BB/(1.D0+X/RR)
         BPL=RS*BTL/((RR+X)*QL)
      ENDIF
      RETURN
      END
C
C     ****** CALCULATE LOCAL WAVE ELECTRIC FIELD ******
C
      SUBROUTINE FPSETV(RHOL,ETAL,RKR,RKTH,RKPH,CER,CETH,CEPH)
C
      INCLUDE 'fpcomm.inc'
C
      IF(MODELW.EQ.3) THEN
         Y=RHOL*RA*SIN(ETAL)
C
         DREWY=1.D1
         REWY=0.D0

         ARG=(Y-REWY)**2/DREWY**2
         IF(ARG.GT.100.D0) THEN
            FACT=0.D0
         ELSE
            FACT=EXP(-ARG)
         ENDIF
         CER= CEWR*FACT
         CETH=CEWTH*FACT
         CEPH=CEWPH*FACT
         RKR= RKWR
         RKTH=RKWTH
         RKPH=RKWPH
      ELSE
C
         CALL FPWMGET(RHOL,ETAL,CEWR1,CEWTH1,CEWPH1,
     &                          CKWR1,CKWTH1,CKWPH1,IERR)
         CER =SQRT(PWAVE)*CEWR1
         CETH=SQRT(PWAVE)*CEWTH1
         CEPH=SQRT(PWAVE)*CEWPH1
         RKR = DBLE(CKWR1)
         RKTH= DBLE(CKWTH1)
         RKPH= DBLE(CKWPH1)
      ENDIF
C
      RETURN
      END

