C     $Id$
C
C   ************************************************
C   **             Graphic output                 **
C   ************************************************
C
C     ****** DRAW 2D POLOIDAL GRAPH ******
C
      SUBROUTINE EQGRAP
C
      INCLUDE 'eqcomc.h'
C
      DIMENSION GX(NXM),GRX(NPSM)
      DIMENSION GYPS(NXM,1),GYJT(NXM,3),GYPP(NXM,1),GYTT(NXM,1)
      DIMENSION GRYPP(NPSM),GRYTT(NPSM)
      DIMENSION XX(NPSM),XPPPS(NPSM),XTTPS(NPSM)
C
      DIMENSION GF(NSGM,NTGM),GR(NSGM,NTGM),GZ(NSGM,NTGM)
      DIMENSION GRS(NTGMP),GZS(NTGMP)
      DIMENSION KA(4,NSGM,NTGM)
      CHARACTER KSTR*80
C
      DO 100 NTG=1,NTGMAX
      DO 100 NSG=1,NSGMAX
         GF(NSG,NTG)=GCLIP(PSI(NTG,NSG))
         GR(NSG,NTG)=GCLIP(RR+SIGM(NSG)*RHOM(NTG)*COS(THGM(NTG)))
         GZ(NSG,NTG)=GCLIP(   SIGM(NSG)*RHOM(NTG)*SIN(THGM(NTG)))
C         CALL SPL2DD(THGM(NTG),SIGM(NSG),PSIL,DPSITL,DPSISL,
C     &               THGMX,SIGMX,U,NTGPM,NTGPMAX,NSGPMAX,IERR)
C         GF(NSG,NTG)=GCLIP(DPSITL)
  100 CONTINUE
C
      DO 200 NTG=1,NTGMAX
         GRS(NTG)=GCLIP(RR+RHOM(NTG)*COS(THGM(NTG)))
         GZS(NTG)=GCLIP(   RHOM(NTG)*SIN(THGM(NTG)))
  200 CONTINUE
         GRS(NTGMAX+1)=GCLIP(RR+RHOM(1)*COS(THGM(1)))
         GZS(NTGMAX+1)=GCLIP(   RHOM(1)*SIN(THGM(1)))
C
      CALL PAGES
      CALL SETCHS(0.3,0.0)
      CALL EQGXEQ(GF,GR,GZ,GRS,GZS,NSGM,NSGMAX,NTGMAX,KA)
      CALL EQGPRM
      CALL PAGEE
C
      DO 300 NTG=1,NTGMAX
      DO 300 NSG=1,NSGMAX
C         CALL SPL2DD(THGM(NTG),SIGM(NSG),PSIL,DPSITL,DPSISL,
C     &               THGMX,SIGMX,U,NTGPM,NTGPMAX,NSGPMAX,IERR)
C        GF(NSG,NTG)=GCLIP(DPSISL)
         GF(NSG,NTG)=GCLIP(HJT(NTG,NSG))
  300 CONTINUE
C
      CALL PAGES
      CALL SETCHS(0.3,0.0)
      CALL EQGXEQ(GF,GR,GZ,GRS,GZS,NSGM,NSGMAX,NTGMAX,KA)
      CALL EQGPRM
      CALL PAGEE
C
      NTGX=NTGMAX/2+1
      DO 500 NSG=NSGMAX,1,-1
         NX=NSGMAX-NSG+1
         XAX(NX)=RMG(NSG,NTGX)
         RPSI(NX)=PSI(NTGX,NSG)
         RHJT(NX)=HJT(NTGX,NSG)
         RHJP(NX)=HJP1(NTGX,NSG)
         RPP(NX)=PP(NTGX,NSG)
         RTT(NX)=TT(NTGX,NSG)
 500  CONTINUE
      NTGX=1
      DO 510 NSG=1,NSGMAX
         NX=NSGMAX+NSG
         XAX(NX)=RMG(NSG,NTGX)
         RPSI(NX)=PSI(NTGX,NSG)
         RHJT(NX)=HJT(NTGX,NSG)
         RHJP(NX)=HJP1(NTGX,NSG)
         RPP(NX)=PP(NTGX,NSG)
         RTT(NX)=TT(NTGX,NSG)
 510  CONTINUE
C
      GX1=2.5
      GX2=12.5
      GY1=10.5
      GY2=17.5
C
      GX3=2.5
      GX4=12.5
      GY3=1.5
      GY4=8.5
C
      GX5=15.0
      GX6=25.0
      GY5=10.5
      GY6=17.5
C
      GX7=15.0
      GX8=25.0
      GY7=1.5
      GY8=8.5
C
      DO 600 NX=1,2*NSGMAX
         GX(NX)=GCLIP(XAX(NX))
         GYPS(NX,1)=GCLIP(RPSI(NX))
         GYJT(NX,1)=GCLIP(-RHJT(NX)*1.D-6)
         GYJT(NX,2)=GCLIP(-RHJP(NX)*1.D-6)
         GYJT(NX,3)=GCLIP(-(RHJT(NX)-RHJP(NX))*1.D-6)
         GYPP(NX,1)=GCLIP(RPP(NX)*1.D-6)
         GYTT(NX,1)=GCLIP(RTT(NX))
  600 CONTINUE
C
      CALL PAGES
      KSTR='/PSI(X)/'
      CALL EQGR1D(GX1,GX2,GY1,GY2,GX,GYPS(1,1),NXM,2*NSGMAX,1,KSTR,0)
      KSTR='/-HJT(X)/'
      CALL EQGR1D(GX3,GX4,GY3,GY4,GX,GYJT(1,1),NXM,2*NSGMAX,3,KSTR,0)
      KSTR='/PP(X)/'
      CALL EQGR1D(GX5,GX6,GY5,GY6,GX,GYPP(1,1),NXM,2*NSGMAX,1,KSTR,0)
      KSTR='/TT(X)/'
      CALL EQGR1D(GX7,GX8,GY7,GY8,GX,GYTT(1,1),NXM,2*NSGMAX,1,KSTR,0)
      CALL PAGEE
C
      DO 1000 NPS=1,NPSMAX
         XX(NPS)=PSIPS(NPS)
         XPPPS(NPS)=PPPS(NPS)
         XTTPS(NPS)=TTPS(NPS)
 1000 CONTINUE
C
      GPX1=2.5
      GPX2=12.5
      GPY1=10.5
      GPY2=17.5
C
      GPX3=2.5
      GPX4=12.5
      GPY3=1.5
      GPY4=8.5
C
      DO 1100 NPS=1,NPSMAX
         GRX(NPS)=GCLIP(XX(NPS))
         GRYPP(NPS)=GCLIP(XPPPS(NPS)*1.D-6)
         GRYTT(NPS)=GCLIP(XTTPS(NPS))
 1100 CONTINUE
C
      CALL PAGES
      KSTR='/PPS(PSI)/'
      CALL EQGR1D(GPX1,GPX2,GPY1,GPY2,GRX,GRYPP,NPSM,NPSMAX,1,KSTR,1)
      KSTR='/TTS(PSI)/'
      CALL EQGR1D(GPX3,GPX4,GPY3,GPY4,GRX,GRYTT,NPSM,NPSMAX,1,KSTR,1)
      CALL PAGEE
C
 9000 RETURN
      END
C
C     ***** Draw Parameters *****
C
      SUBROUTINE EQGPRM
C
      INCLUDE 'eqcomc.h'
C
      REAL*4 XPOS,YPOS,DELY
C
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(0,0,7)
C
      XPOS=16.5
      YPOS=15.0
      DELY= 0.5
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('RR    :',7)
      CALL NUMBD(RR,  '(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('RA    :',7)
      CALL NUMBD(RA, '(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('RKAP  :',7)
      CALL NUMBD(RKAP,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('RDLT  :',7)
      CALL NUMBD(RDLT,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('RIP   :',7)
      CALL NUMBD(RIP,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('BB    :',7)
      CALL NUMBD(BB,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('PP0   :',7)
      CALL NUMBD(PP0,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('PP1   :',7)
      CALL NUMBD(PP1,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('PP2   :',7)
      CALL NUMBD(PP2,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('PJ0   :',7)
      CALL NUMBD(PJ0,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('PJ1   :',7)
      CALL NUMBD(PJ1,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('PJ2   :',7)
      CALL NUMBD(PJ2,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('R0:',3)
      CALL NUMBD(PROFR0,'(F5.1)',5)
      CALL NUMBD(PROFP0,'(F5.1)',5)
      CALL NUMBD(PROFJ0,'(F5.1)',5)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('R1:',3)
      CALL NUMBD(PROFR1,'(F5.1)',5)
      CALL NUMBD(PROFP1,'(F5.1)',5)
      CALL NUMBD(PROFJ1,'(F5.1)',5)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('R2:',3)
      CALL NUMBD(PROFR2,'(F5.1)',5)
      CALL NUMBD(PROFP2,'(F5.1)',5)
      CALL NUMBD(PROFJ2,'(F5.1)',5)
      YPOS=YPOS-DELY
C
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('EPSEQ :',7)
      CALL NUMBD(EPSEQ,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('NSGMAX :',7)
      CALL NUMBI(NSGMAX,'(I8)',8)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('NTGMAX:',7)
      CALL NUMBI(NTGMAX,'(I8)',8)
      YPOS=YPOS-DELY
C
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('RAXIS :',7)
      CALL NUMBD(RAXIS,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('ZAXIS :',7)
      CALL NUMBD(ZAXIS,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('SAXIS :',7)
      CALL NUMBD(SAXIS,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('PVOL  :',7)
      CALL NUMBD(PVOL,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('RAAVE :',7)
      CALL NUMBD(RAAVE,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('BETA  :',7)
      CALL NUMBD(BETA,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('BETAP :',7)
      CALL NUMBD(BETAP,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('QAXIS :',7)
      CALL NUMBD(QAXIS,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      CALL MOVE(XPOS,YPOS)
      CALL TEXT('QSURF :',7)
      CALL NUMBD(QSURF,'(1PE11.3)',11)
      YPOS=YPOS-DELY
C
      RETURN
      END
