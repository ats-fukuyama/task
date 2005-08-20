C     $Id$
C
C   ************************************************
C   **             Graphic output                 **
C   ************************************************
C
      SUBROUTINE EQGOUT(MODE)
C
      INCLUDE 'eqcomc.inc'
C
      CHARACTER KSTR*2,K1*1,K2*1
C
    1 IF(MODE.EQ.1) THEN
        WRITE(6,*) ' ## INPUT KID : C,C1,C2 ',
     &                             'S,S1,S2,SR,ST,SD,SB M A X/EXIT'
      ELSE
        WRITE(6,*) ' ## INPUT KID : S,S1,S2,SR,ST,SD,SB X/EXIT'
      ENDIF
      READ(5,'(A2)',ERR=1,END=9000) KSTR
      K1=KSTR(1:1)
      K2=KSTR(2:2)
      CALL GUCPTL(K1)
      CALL GUCPTL(K2)
      IF(K1.EQ.'C') THEN
         IF(MODE.EQ.1) THEN
            IF(K2.EQ.' ') THEN
               CALL EQGC2D
               CALL EQGC1D
            ELSEIF(K2.EQ.'1') THEN
               CALL EQGC1D
            ELSEIF(K2.EQ.'2') THEN
               CALL EQGC2D
            ENDIF
         ELSE
            WRITE(6,*) 'XX: EQGOUT: NO DATA CREATED!'
         ENDIF
      ELSEIF(K1.EQ.'S') THEN
         IF(K2.EQ.' ') THEN
            CALL EQGS2D
            CALL EQGS1D(0)
         ELSEIF(K2.EQ.'1') THEN
            CALL EQGS1D(0)
         ELSEIF(K2.EQ.'2') THEN
            CALL EQGS2D
         ELSEIF(K2.EQ.'R') THEN
            CALL EQGS1D(1)
         ELSEIF(K2.EQ.'T') THEN
            CALL EQGS1D(2)
         ELSEIF(K2.EQ.'D') THEN
            CALL EQGSDD
         ELSEIF(K2.EQ.'B') THEN
            CALL EQGSBB
         ENDIF
      ELSEIF (K1.EQ.'M') THEN
         CALL EQGC1M
      ELSEIF (K1.EQ.'A') THEN
         CALL EQGC2D
         CALL EQGC1D
         CALL EQGS2D
         CALL EQGS1D(0)
         CALL EQGSDD
      ELSEIF (K1.EQ.'X') THEN
         GOTO 9000
      ENDIF
      GOTO 1
C
 9000 RETURN
      END
C
C     ****** DRAW CALCULATED EQ1D GRAPH ******
C
      SUBROUTINE EQGC1D
C
      INCLUDE 'eqcomc.inc'
C
      DIMENSION GX(NXM)
      DIMENSION GYPS(NXM,1),GYJT(NXM,3),GYPP(NXM,1),GYTT(NXM,1)
      DIMENSION GRX(NPSM)
      DIMENSION GRYPP(NPSM),GRYTT(NPSM)
      DIMENSION GRYTE(NPSM),GRYOM(NPSM)
      DIMENSION GVX(NRVM)
      DIMENSION GVY1(NRVM,1),GVY2(NRVM,1),GVY3(NRVM,2),GVY4(NRVM,1)
C
      CHARACTER KSTR*80
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
C     ----- Major radius dependence -----
C
      NTGX=NTGMAX/2+1
      DO NSG=NSGMAX,1,-1
         NX=NSGMAX-NSG+1
         XAX(NX)=RMG(NSG,NTGX)
         RPSI(NX)=PSI(NTGX,NSG)
         RHJT(NX)=HJT(NTGX,NSG)
         RHJP(NX)=HJP2(NTGX,NSG)
         RPP(NX)=PP(NTGX,NSG)
         RTT(NX)=TT(NTGX,NSG)
      ENDDO
      NTGX=1
      DO NSG=1,NSGMAX
         NX=NSGMAX+NSG
         XAX(NX)=RMG(NSG,NTGX)
         RPSI(NX)=PSI(NTGX,NSG)
         RHJT(NX)=HJT(NTGX,NSG)
         RHJP(NX)=HJP2(NTGX,NSG)
         RPP(NX)=PP(NTGX,NSG)
         RTT(NX)=TT(NTGX,NSG)
      ENDDO
C
      DO NX=1,2*NSGMAX
         GX(NX)=GUCLIP(XAX(NX))
         GYPS(NX,1)=GUCLIP(RPSI(NX))
         GYJT(NX,1)=GUCLIP(RHJT(NX)*1.D-6)
         GYJT(NX,2)=GUCLIP(RHJP(NX)*1.D-6)
         GYJT(NX,3)=GUCLIP((RHJT(NX)-RHJP(NX))*1.D-6)
         GYPP(NX,1)=GUCLIP(RPP(NX)*1.D-6)
         GYTT(NX,1)=GUCLIP(RTT(NX))
      ENDDO
C
      CALL PAGES
      KSTR='/PSI(X)/'
      CALL EQGR1D(GX1,GX2,GY1,GY2,GX,GYPS(1,1),NXM,2*NSGMAX,1,KSTR,0)
      KSTR='/-HJT(X)/'
      CALL EQGR1D(GX3,GX4,GY3,GY4,
     &            GX,GYJT(1,1),NXM,2*NSGMAX,3,KSTR,0)
      KSTR='/PP(X)/'
      CALL EQGR1D(GX5,GX6,GY5,GY6,GX,GYPP(1,1),NXM,2*NSGMAX,1,KSTR,0)
      KSTR='/TT(X)/'
      CALL EQGR1D(GX7,GX8,GY7,GY8,GX,GYTT(1,1),NXM,2*NSGMAX,1,KSTR,0)
      CALL PAGEE
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
      GPX5=15.0
      GPX6=25.0
      GPY5=10.5
      GPY6=17.5
C
      GPX7=15.0
      GPX8=25.0
      GPY7=1.5
      GPY8=8.5
C
C     ----- Poloidal flux dependence -----
C
      DO NPS=1,NPSMAX
         GRX(NPS)=GUCLIP(PSIPS(NPS))
         GRYPP(NPS)=GUCLIP(PPPS(NPS)*1.D-6)
         GRYTT(NPS)=GUCLIP(TTPS(NPS))
         GRYTE(NPS)=GUCLIP(TEPS(NPS))
         GRYOM(NPS)=GUCLIP(OMPS(NPS))
      ENDDO
C
      CALL PAGES
      KSTR='/PPS(PSI)/'
      CALL EQGR1D(GPX1,GPX2,GPY1,GPY2,GRX,GRYPP,NPSM,NPSMAX,1,KSTR,1)
      KSTR='/TTS(PSI)/'
      CALL EQGR1D(GPX3,GPX4,GPY3,GPY4,GRX,GRYTT,NPSM,NPSMAX,1,KSTR,1)
      KSTR='/TEMP(PSI)/'
      CALL EQGR1D(GPX5,GPX6,GPY5,GPY6,GRX,GRYTE,NPSM,NPSMAX,1,KSTR,1)
      KSTR='/OMEGA(PSI)/'
      CALL EQGR1D(GPX7,GPX8,GPY7,GPY8,GRX,GRYOM,NPSM,NPSMAX,1,KSTR,1)
      CALL PAGEE
C
C     ----- Flux average quantities -----
C
      DO NRV=1,NRVMAX
         GVX(NRV)=GUCLIP(PSIPV(NRV))
         GVY1(NRV,1)=GUCLIP(PSITV(NRV))
         GVY2(NRV,1)=GUCLIP(QPV(NRV))
         GVY3(NRV,1)=GUCLIP(AVBB2(NRV))
         GVY3(NRV,2)=GUCLIP(AVIR2(NRV))
         GVY4(NRV,1)=GUCLIP(AVRR2(NRV))
      ENDDO
C
      CALL PAGES
      KSTR='/PSIT(PSIP)/'
      CALL EQGR1D(GPX1,GPX2,GPY1,GPY2,GVX,GVY1,NRVM,NRVMAX,1,KSTR,1)
      KSTR='/QPV(PSIP)/'
      CALL EQGR1D(GPX3,GPX4,GPY3,GPY4,GVX,GVY2,NRVM,NRVMAX,1,KSTR,1)
      KSTR='/AVBB2,IR2(PSIP)/'
      CALL EQGR1D(GPX5,GPX6,GPY5,GPY6,GVX,GVY3,NRVM,NRVMAX,2,KSTR,1)
      KSTR='/AVRR2(PSIP)/'
      CALL EQGR1D(GPX7,GPX8,GPY7,GPY8,GVX,GVY4,NRVM,NRVMAX,1,KSTR,1)
      CALL PAGEE
C
      RETURN
      END
C
C     ****** DRAW CALCULATED EQ2D GRAPH ******
C
      SUBROUTINE EQGC2D
C
      INCLUDE 'eqcomc.inc'
C
      DIMENSION GF(NSGM,NTGM),GR(NSGM,NTGM),GZ(NSGM,NTGM)
      DIMENSION GRS(NTGMP),GZS(NTGMP)
      DIMENSION KA(4,NSGM,NTGM)
C
      DO NTG=1,NTGMAX
      DO NSG=1,NSGMAX
         GR(NSG,NTG)=GUCLIP(RR+SIGM(NSG)*RHOM(NTG)*COS(THGM(NTG)))
         GZ(NSG,NTG)=GUCLIP(   SIGM(NSG)*RHOM(NTG)*SIN(THGM(NTG)))
         GF(NSG,NTG)=GUCLIP(PSI(NTG,NSG))
      ENDDO
      ENDDO
C
      DO NTG=1,NTGMAX
         GRS(NTG)=GUCLIP(RR+RHOM(NTG)*COS(THGM(NTG)))
         GZS(NTG)=GUCLIP(   RHOM(NTG)*SIN(THGM(NTG)))
      ENDDO
         GRS(NTGMAX+1)=GUCLIP(RR+RHOM(1)*COS(THGM(1)))
         GZS(NTGMAX+1)=GUCLIP(   RHOM(1)*SIN(THGM(1)))
C
      CALL PAGES
      CALL SETCHS(0.3,0.0)
      CALL EQGR2D(GF,GR,GZ,GRS,GZS,NSGM,NSGMAX,NTGMAX,KA,'/PSI/')
      CALL EQGPRM
      CALL PAGEE
C
      DO NTG=1,NTGMAX
      DO NSG=1,NSGMAX
         GF(NSG,NTG)=GUCLIP(HJT(NTG,NSG))
      ENDDO
      ENDDO
C
      CALL PAGES
      CALL SETCHS(0.3,0.0)
      CALL EQGR2D(GF,GR,GZ,GRS,GZS,NSGM,NSGMAX,NTGMAX,KA,'/HJT/')
      CALL EQGPRM
      CALL PAGEE
C
      DO NTG=1,NTGMAX
      DO NSG=1,NSGMAX
         GF(NSG,NTG)=GUCLIP(RHO(NTG,NSG))
      ENDDO
      ENDDO
C
      CALL PAGES
      CALL SETCHS(0.3,0.0)
      CALL EQGR2D(GF,GR,GZ,GRS,GZS,NSGM,NSGMAX,NTGMAX,KA,'/RHOM/')
      CALL EQGPRM
      CALL PAGEE
C
      RETURN
      END
C
C     ****** DRAW SPLINED EQ1D GRAPH ******
C
      SUBROUTINE EQGS1D(MODE)
C
      INCLUDE 'eqcomq.inc'
C
      DIMENSION GX(NRM),GY(NRM,4)
C
      IF(MODE.EQ.0) THEN
         DO NR=1,NRMAX
            GX(NR)=GUCLIP(PSIP(NR))
         ENDDO
      ELSEIF(MODE.EQ.1) THEN
         DO NR=1,NRMAX
            GX(NR)=GUCLIP(RHOT(NR))
         ENDDO
      ELSEIF(MODE.EQ.2) THEN
         DO NR=1,NRMAX
            GX(NR)=GUCLIP(PSIT(NR))
         ENDDO
      ENDIF

      CALL PAGES
      CALL SETCHS(0.35,0.0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(PPS(NR)*1.D-6)
      ENDDO
      CALL EQGR1D( 3.0,13.0,10.0,16.0,GX,GY,NRM,NRMAX,1,'@PPS@',0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(TTS(NR))
      ENDDO
      CALL EQGR1D(15.0,25.0,10.0,16.0,GX,GY,NRM,NRMAX,1,'@TTS@',1)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(RLEN(NR))
      ENDDO
      CALL EQGR1D( 3.0,13.0, 2.0, 8.0,GX,GY,NRM,NRMAX,1,'@RLEN@',0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(QPS(NR))
      ENDDO
      CALL EQGR1D(15.0,25.0, 2.0, 8.0,GX,GY,NRM,NRMAX,1,'@QPS@',0)
C
      CALL PAGEE
C
      CALL PAGES
      CALL SETCHS(0.35,0.0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(VPS(NR))
      ENDDO
      CALL EQGR1D( 3.0,13.0,10.0,16.0,GX,GY,NRM,NRMAX,1,'@VPS@',0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(SPS(NR))
      ENDDO
      CALL EQGR1D(15.0,25.0,10.0,16.0,GX,GY,NRM,NRMAX,1,'@SPS@',0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(RRMIN(NR))
         GY(NR,2)=GUCLIP(RRMAX(NR))
         GY(NR,3)=GUCLIP(ZZMIN(NR))
         GY(NR,4)=GUCLIP(ZZMAX(NR))
      ENDDO
      CALL EQGR1D( 3.0,13.0, 2.0, 8.0,GX,GY,NRM,NRMAX,4,'@RRMIN/MAX@',0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(BBMIN(NR))
         GY(NR,2)=GUCLIP(BBMAX(NR))
      ENDDO
      CALL EQGR1D(15.0,25.0, 2.0, 8.0,GX,GY,NRM,NRMAX,2,'@BBMIN/MAX@',0)
C
      CALL PAGEE
C
      CALL PAGES
      CALL SETCHS(0.35,0.0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(AVERHR(NR))
      ENDDO
      CALL EQGR1D( 3.0,13.0,10.0,16.0,GX,GY,NRM,NRMAX,1,'@AVERHR@',0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(AVEIR2(NR))
      ENDDO
      CALL EQGR1D(15.0,25.0,10.0,16.0,GX,GY,NRM,NRMAX,1,'@AVEIR2@',0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(AVERH1(NR))
         GY(NR,2)=GUCLIP(AVERH2(NR))
      ENDDO
      CALL EQGR1D( 3.0,13.0, 2.0, 8.0,GX,GY,NRM,NRMAX,2,
     &            '@AVERH1,AVERH2@',0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(AVERHB(NR))
      ENDDO
      CALL EQGR1D(15.0,25.0, 2.0, 8.0,GX,GY,NRM,NRMAX,1,'@AVERHB@',0)
C
      CALL PAGEE
C
      CALL PAGES
      CALL SETCHS(0.35,0.0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(AVEBB2(NR))
      ENDDO
      CALL EQGR1D( 3.0,13.0,10.0,16.0,GX,GY,NRM,NRMAX,1,'@AVEBB2@',0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(AVEIB2(NR))
      ENDDO
      CALL EQGR1D(15.0,25.0,10.0,16.0,GX,GY,NRM,NRMAX,1,'@AVEIB2@',0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(AVEJPR(NR))*1.E-6
         GY(NR,2)=GUCLIP(AVEJTR(NR))*1.E-6
      ENDDO
      CALL EQGR1D( 3.0,13.0, 2.0, 8.0,GX,GY,NRM,NRMAX,2,
     &                                           '@AVEJPR,AVEJTR@',0)
C
      DO NR=1,NRMAX
         GY(NR,1)=GUCLIP(PSIT(NR))
      ENDDO
      CALL EQGR1D(15.0,25.0, 2.0, 8.0,GX,GY,NRM,NRMAX,1,'@PSIT@',0)
C
      CALL PAGEE
C
      RETURN
      END
C
C     ****** DRAW SPLINED EQ2D GRAPH ******
C
      SUBROUTINE EQGS2D
C
      INCLUDE 'eqcomq.inc'
C
      DIMENSION GR(NTHMP),GZ(NTHMP)
      DIMENSION GPSIRZ(NRGM,NZGM),GRG(NRGM),GZG(NZGM)
      DIMENSION KA(8,NRGM,NZGM)
      CHARACTER KTITL*80
C
      DO NR=1,NRGMAX
         GRG(NR)=GUCLIP(RG(NR))
      ENDDO
      DO NZ=1,NZGMAX
         GZG(NZ)=GUCLIP(ZG(NZ))
      ENDDO
      DO NZ=1,NZGMAX
         DO NR=1,NRGMAX
            GPSIRZ(NR,NZ)=GUCLIP(PSIRZ(NR,NZ))
         ENDDO
      ENDDO
C
      CALL GMNMX1(GRG,1,NRGMAX,1,GRMIN,GRMAX)
      CALL GMNMX1(GZG,1,NZGMAX,1,GZMIN,GZMAX)
      CALL GMNMX2(GPSIRZ,NRGM,1,NRGMAX,1,1,NZGMAX,1,GPMIN,GPMAX)
C      WRITE(6,'(A,1P2E12.4)') 'GRMIN,GRMAX=',GRMIN,GRMAX
C      WRITE(6,'(A,1P2E12.4)') 'GZMIN,GZMAX=',GZMIN,GZMAX
C      WRITE(6,'(A,1P2E12.4)') 'GPMIN,GPMAX=',GPMIN,GPMAX
      NPSTEP=20
      GPORG=GPMIN
      GPSTEP=(GPMAX-GPMIN)/NPSTEP
C
      GRLEN=GRMAX-GRMIN
      GZLEN=GZMAX-GZMIN
      IF(GRLEN.GT.GZLEN) THEN
         GPR=15.0
         GPZ=15.0*GZLEN/GRLEN
      ELSE
         GPR=15.0*GRLEN/GZLEN
         GPZ=15.0
      ENDIF
C
      CALL PAGES
      CALL MOVE(2.0,17.5)
      KTITL='/PSIRZ contour in R-Z/'
      CALL TEXTX(KTITL)
C
      CALL GQSCAL(GRMIN,GRMAX,GGRMIN,GGRMAX,GGRSTP)
      CALL GQSCAL(GZMIN,GZMAX,GGZMIN,GGZMAX,GGZSTP)
C
      CALL SETLIN(0,2,7)
      CALL GDEFIN(2.0,2.0+GPR,2.0,2.0+GPZ,
     &            GRMIN,GRMAX,GZMIN,GZMAX)
      CALL GFRAME
      CALL GSCALE(GGRMIN+GGRSTP,GGRSTP,0.0,0.0,0.1,9)
      CALL GVALUE(GGRMIN+GGRSTP,GGRSTP*2,0.0,0.0,1)
      CALL GSCALE(0.0,0.0,0.0,GGZSTP,0.1,9)
      CALL GVALUE(0.0,0.0,0.0,GGZSTP*2,1)
C
      CALL SETRGB(1.0,0.0,0.0)
      CALL SETCLP(2.0,2.0+GPR,2.0,2.0+GPZ)
      CALL CONTQ2(GPSIRZ,GRG,GZG,NRGM,NRGMAX,NZGMAX,
     &            0.5*GPSTEP,GPSTEP,NPSTEP,0,0,KA)
      CALL SETRGB(0.0,0.0,1.0)
      CALL CONTQ2(GPSIRZ,GRG,GZG,NRGM,NRGMAX,NZGMAX,
     &            -0.5*GPSTEP,-GPSTEP,NPSTEP,0,0,KA)
      CALL OFFCLP
      CALL SETRGB(0.0,0.0,0.0)
      CALL EQGPRM
      CALL PAGEE
C
      GRMIN=GUCLIP(RGMIN)
      GRMAX=GUCLIP(RGMAX)
      GZMIN=GUCLIP(ZGMIN)
      GZMAX=GUCLIP(ZGMAX)
      GRLEN=GRMAX-GRMIN
      GZLEN=GZMAX-GZMIN
      IF(GRLEN.GT.GZLEN) THEN
         GPR=15.0
         GPZ=15.0*GZLEN/GRLEN
      ELSE
         GPR=15.0*GRLEN/GZLEN
         GPZ=15.0
      ENDIF
C
      CALL PAGES
      CALL MOVE(2.0,17.5)
      KTITL='/PSI contour in R-Z/'
      CALL TEXTX(KTITL)
C
      CALL GQSCAL(GRMIN,GRMAX,GGRMIN,GGRMAX,GGRSTP)
      CALL GQSCAL(GZMIN,GZMAX,GGZMIN,GGZMAX,GGZSTP)
C
      CALL SETLIN(0,2,7)
      CALL GDEFIN(2.0,2.0+GPR,2.0,2.0+GPZ,
     &            GRMIN,GRMAX,GZMIN,GZMAX)
      CALL GFRAME
      CALL GSCALE(GGRMIN+GGRSTP,GGRSTP,0.0,0.0,0.1,9)
      CALL GVALUE(GGRMIN+GGRSTP,GGRSTP*2,0.0,0.0,1)
      CALL GSCALE(0.0,0.0,0.0,GGZSTP,0.1,9)
      CALL GVALUE(0.0,0.0,0.0,GGZSTP*2,1)
C
      CALL SETRGB(1.0,0.0,0.0)
      DO NR=1,NRPMAX
         DO NTH=1,NTHMAX
            GR(NTH)=GUCLIP(RPS(NTH,NR))
            GZ(NTH)=GUCLIP(ZPS(NTH,NR))
         ENDDO
         GR(NTHMAX+1)=GR(1)
         GZ(NTHMAX+1)=GZ(1)
         CALL GPLOTP(GR,GZ,1,NTHMAX+1,1,0,0,0)
      ENDDO
C
      CALL SETRGB(0.0,0.0,1.0)
      DO NR=NRPMAX+1,NRMAX
         DO NTH=1,NTHMAX
            GR(NTH)=GUCLIP(RPS(NTH,NR))
            GZ(NTH)=GUCLIP(ZPS(NTH,NR))
         ENDDO
         GR(NTHMAX+1)=GR(1)
         GZ(NTHMAX+1)=GZ(1)
         CALL GPLOTP(GR,GZ,1,NTHMAX+1,1,0,0,0)
      ENDDO
C
      CALL SETRGB(0.0,0.0,0.0)
      CALL EQGPRM
      CALL PAGEE
C
      RETURN
      END
C
C     ****** DRAW CALCULATED EQ2D DERIVATIVE GRAPH ******
C
      SUBROUTINE EQGSDD
C
      INCLUDE 'eqcomq.inc'
C
      DIMENSION GF(NRM,NTHM),GR(NRM,NTHM),GZ(NRM,NTHM)
      DIMENSION GRSU(NSUM),GZSU(NSUM)
      DIMENSION GRSW(NSUM),GZSW(NSUM)
      DIMENSION KA(4,NRM,NTHM)
      CHARACTER KTITL*80
C
      DO NR=1,NRMAX
      DO NTH=1,NTHMAX
         GR(NR,NTH)=GUCLIP(RPS(NTH,NR))
         GZ(NR,NTH)=GUCLIP(ZPS(NTH,NR))
      ENDDO
      ENDDO
C
      DO IND=1,8
C
         IF(IND.EQ.1) THEN
            DO NR=1,NRMAX
            DO NTH=1,NTHMAX
C               R2=(RPS(NTH,NR)-RAXIS)**2+(ZPS(NTH,NR)-ZAXIS)**2
               R2=RHOT(NR)**2
               GF(NR,NTH)=GUCLIP(DRPSI(NTH,NR)*R2)
            ENDDO
            ENDDO
            KTITL='/DRPSI*r^2/'
         ELSEIF(IND.EQ.2) THEN
            DO NR=1,NRMAX
            DO NTH=1,NTHMAX
C               R2=(RPS(NTH,NR)-RAXIS)**2+(ZPS(NTH,NR)-ZAXIS)**2
               R2=RHOT(NR)**2
               GF(NR,NTH)=GUCLIP(DZPSI(NTH,NR)*R2)
            ENDDO
            ENDDO
            KTITL='/DZPSI*r^2/'
         ELSEIF(IND.EQ.3) THEN
            DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               GF(NR,NTH)=GUCLIP(DRCHI(NTH,NR))
            ENDDO
            ENDDO
            KTITL='/DRCHI/'
         ELSEIF(IND.EQ.4) THEN
            DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               GF(NR,NTH)=GUCLIP(DZCHI(NTH,NR))
            ENDDO
            ENDDO
            KTITL='/DZCHI/'
         ELSEIF(IND.EQ.5) THEN
            DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               GF(NR,NTH)=GUCLIP(RPS(NTH,NR))
            ENDDO
            ENDDO
            KTITL='/RPS/'
         ELSEIF(IND.EQ.6) THEN
            DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               GF(NR,NTH)=GUCLIP(ZPS(NTH,NR))
            ENDDO
            ENDDO
            KTITL='/ZPS/'
         ELSEIF(IND.EQ.7) THEN
            DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               GF(NR,NTH)=GUCLIP(PSIP(NR))
            ENDDO
            ENDDO
            KTITL='/PSIP/'
         ELSEIF(IND.EQ.8) THEN
            DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               GF(NR,NTH)=GUCLIP((2.D0*PI*(NTH-1))/NTHMAX)
            ENDDO
            ENDDO
            KTITL='/CHI/'
         ENDIF
C
         GSUM=0.0
         DO NTH=1,NTHMAX
            GSUM=GSUM+GF(1,NTH)
         ENDDO
         GSUM=GSUM/NTHMAX
         DO NTH=1,NTHMAX
            GF(1,NTH)=GSUM
         ENDDO
C
         IF(IND.NE.8) THEN
            CALL GMNMX2(GF,NRM,1,NRMAX,1,1,NTHMAX,1,GFMIN,GFMAX)
            CALL GQSCAL(GFMIN,GFMAX,GGFMIN,GGFMAX,GGFSTP)
            GGFSTP=0.5*GGFSTP
            NSTEP=INT((GGFMAX-GGFMIN)/GGFSTP)+1
         ELSE
            GFMIN=0.0
            GFMAX=GUCLIP(2.D0*PI)
            GGFMIN=0.0
            NSTEP=12
            GGFSTP=GFMAX/NSTEP
         ENDIF
C
         DO NSU=1,NSUMAX
            GRSU(NSU)=GUCLIP(RSU(NSU))
            GZSU(NSU)=GUCLIP(ZSU(NSU))
            GRSW(NSU)=GUCLIP(RSW(NSU))
            GZSW(NSU)=GUCLIP(ZSW(NSU))
         ENDDO
         GRSU(NSUMAX+1)=GUCLIP(RSU(NSU))
         GZSU(NSUMAX+1)=GUCLIP(ZSU(NSU))
         GRSW(NSUMAX+1)=GUCLIP(RSW(NSU))
         GZSW(NSUMAX+1)=GUCLIP(ZSW(NSU))
C
         GRLEN=GUCLIP(RGMAX-RGMIN)
         GZLEN=GUCLIP(ZGMAX-ZGMIN)
         IF(GRLEN.GT.GZLEN) THEN
            GPR=15.0
            GPZ=15.0*GZLEN/GRLEN
         ELSE
            GPR=15.0*GRLEN/GZLEN
            GPZ=15.0
         ENDIF
C
         CALL PAGES
         CALL MOVE(2.0,17.5)
         CALL TEXTX(KTITL)
C
         CALL SETLIN(0,2,7)
         CALL GDEFIN(2.0,2.0+GPR,2.0,2.0+GPZ,
     &               REAL(RGMIN),REAL(RGMAX),
     &               REAL(ZGMIN),REAL(ZGMAX))
         CALL GFRAME
         CALL SETCLP(2.0,2.0+GPR,2.0,2.0+GPZ)
C
         IF(GFMIN*GFMAX.GE.0.) THEN
            CALL CONTP5(GF,GR,GZ,NRM,NRMAX,NTHMAX,
     &                  GGFMIN,GGFSTP,NSTEP,2,0,KA)
         ELSE
            CALL CONTP5(GF,GR,GZ,NRM,NRMAX,NTHMAX,
     &                   0.5*GGFSTP, GGFSTP,NSTEP,2,0,KA)
            CALL CONTP5(GF,GR,GZ,NRM,NRMAX,NTHMAX,
     &                  -0.5*GGFSTP,-GGFSTP,NSTEP,2,2,KA)
         ENDIF
         CALL OFFCLP
         CALL SETLIN(-1,-1,6)
         CALL GPLOTP(GRSU,GZSU,1,NSUMAX+1,1,0,0,0)
         CALL SETLIN(-1,-1,5)
         CALL GPLOTP(GRSW,GZSW,1,NSUMAX+1,1,0,0,0)
C
         CALL SETLIN(0,-1,7)
         CALL MOVE(20.0,17.0)
         CALL TEXT('MAX :',5)
         CALL NUMBR(GFMAX,'(1PE12.4)',12)
         CALL MOVE(20.0,16.5)
         CALL TEXT('MIN :',5)
         CALL NUMBR(GFMIN,'(1PE12.4)',12)
         CALL MOVE(20.0,16.0)
         CALL TEXT('STEP:',5)
         CALL NUMBR(GGFSTP,'(1PE12.4)',12)
C
         CALL PAGEE
C
      ENDDO
C
      RETURN
      END
C
C     ****** DRAW CALCULATED MAGNETIC FIELD ******
C
      SUBROUTINE EQGSBB
C
      INCLUDE 'eqcomq.inc'
C
      DIMENSION GF(NRM,NTHM),GR(NRM,NTHM),GZ(NRM,NTHM)
      DIMENSION GRSU(NSUM),GZSU(NSUM)
      DIMENSION GRSW(NSUM),GZSW(NSUM)
      DIMENSION KA(4,NRM,NTHM)
      CHARACTER KTITL*80
C
      DO NR=1,NRMAX
      DO NTH=1,NTHMAX
         GR(NR,NTH)=GUCLIP(RPS(NTH,NR))
         GZ(NR,NTH)=GUCLIP(ZPS(NTH,NR))
      ENDDO
      ENDDO
C
      DO IND=1,4
C
         IF(IND.EQ.1) THEN
            DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               GF(NR,NTH)=GUCLIP(BPR(NTH,NR))
            ENDDO
            ENDDO
            KTITL='/BPR/'
         ELSEIF(IND.EQ.2) THEN
            DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               GF(NR,NTH)=GUCLIP(BPZ(NTH,NR))
            ENDDO
            ENDDO
            KTITL='/BPZ/'
         ELSEIF(IND.EQ.3) THEN
            DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               GF(NR,NTH)=GUCLIP(BPT(NTH,NR))
            ENDDO
            ENDDO
            KTITL='/BPT/'
         ELSEIF(IND.EQ.4) THEN
            DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               GF(NR,NTH)=GUCLIP(BTP(NTH,NR))
            ENDDO
            ENDDO
            KTITL='/BTP/'
         ENDIF
C
         GSUM=0.0
         DO NTH=1,NTHMAX
            GSUM=GSUM+GF(1,NTH)
         ENDDO
         GSUM=GSUM/NTHMAX
         DO NTH=1,NTHMAX
            GF(1,NTH)=GSUM
         ENDDO
C
         CALL GMNMX2(GF,NRM,1,NRMAX,1,1,NTHMAX,1,GFMIN,GFMAX)
         CALL GQSCAL(GFMIN,GFMAX,GGFMIN,GGFMAX,GGFSTP)
         GGFSTP=0.5*GGFSTP
         NSTEP=INT((GGFMAX-GGFMIN)/GGFSTP)+1
         DO NSU=1,NSUMAX
            GRSU(NSU)=GUCLIP(RSU(NSU))
            GZSU(NSU)=GUCLIP(ZSU(NSU))
            GRSW(NSU)=GUCLIP(RSW(NSU))
            GZSW(NSU)=GUCLIP(ZSW(NSU))
         ENDDO
         GRSU(NSUMAX+1)=GUCLIP(RSU(NSU))
         GZSU(NSUMAX+1)=GUCLIP(ZSU(NSU))
         GRSW(NSUMAX+1)=GUCLIP(RSW(NSU))
         GZSW(NSUMAX+1)=GUCLIP(ZSW(NSU))
C
         GRLEN=GUCLIP(RGMAX-RGMIN)
         GZLEN=GUCLIP(ZGMAX-ZGMIN)
         IF(GRLEN.GT.GZLEN) THEN
            GPR=15.0
            GPZ=15.0*GZLEN/GRLEN
         ELSE
            GPR=15.0*GRLEN/GZLEN
            GPZ=15.0
         ENDIF
C
         CALL PAGES
         CALL MOVE(2.0,17.5)
         CALL TEXTX(KTITL)
C
         CALL SETLIN(0,2,7)
         CALL GDEFIN(2.0,2.0+GPR,2.0,2.0+GPZ,
     &               REAL(RGMIN),REAL(RGMAX),
     &               REAL(ZGMIN),REAL(ZGMAX))
         CALL GFRAME
         CALL SETCLP(2.0,2.0+GPR,2.0,2.0+GPZ)
C
         IF(GFMIN*GFMAX.GT.0.) THEN
            CALL CONTP5(GF,GR,GZ,NRM,NRMAX,NTHMAX,
     &                  GGFMIN,GGFSTP,NSTEP,2,0,KA)
         ELSE
            CALL CONTP5(GF,GR,GZ,NRM,NRMAX,NTHMAX,
     &                   0.5*GGFSTP, GGFSTP,NSTEP,2,0,KA)
            CALL CONTP5(GF,GR,GZ,NRM,NRMAX,NTHMAX,
     &                  -0.5*GGFSTP,-GGFSTP,NSTEP,2,2,KA)
         ENDIF
         CALL OFFCLP
         CALL SETLIN(-1,-1,6)
         CALL GPLOTP(GRSU,GZSU,1,NSUMAX+1,1,0,0,0)
         CALL SETLIN(-1,-1,5)
         CALL GPLOTP(GRSW,GZSW,1,NSUMAX+1,1,0,0,0)
C
         CALL SETLIN(0,-1,7)
         CALL MOVE(20.0,17.0)
         CALL TEXT('MAX :',5)
         CALL NUMBR(GFMAX,'(1PE12.4)',12)
         CALL MOVE(20.0,16.5)
         CALL TEXT('MIN :',5)
         CALL NUMBR(GFMIN,'(1PE12.4)',12)
         CALL MOVE(20.0,16.0)
         CALL TEXT('STEP:',5)
         CALL NUMBR(GGFSTP,'(1PE12.4)',12)
C
         CALL PAGEE
C
      ENDDO
C
      RETURN
      END
C
C     ***** Draw Parameters *****
C
      SUBROUTINE EQGPRM
C
      INCLUDE 'eqcomc.inc'
C
      REAL*4 XPOS,YPOS,DELY
C
      CALL SETCHS(0.3,0.0)
      CALL SETLIN(0,-1,7)
C
C      XPOS=16.5
      XPOS=20.0
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
      CALL TEXT('RIPX  :',7)
      CALL NUMBD(RIPX,'(1PE11.3)',11)
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
      CALL TEXT('PSI0  :',7)
      CALL NUMBD(PSI0,'(1PE11.3)',11)
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
      CALL TEXT('BETAT :',7)
      CALL NUMBD(BETAT,'(1PE11.3)',11)
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
C
C     ****** DRAW CALCULATED EQ1D MULTI GRAPH ******
C
      SUBROUTINE EQGC1M
C
      INCLUDE 'eqcomc.inc'
C
      DIMENSION GX(NXM),GYPS(NXM,5)
C
      CHARACTER KSTR*80
C
      GX1=2.5
      GX2=25.0
      GY1=1.5
      GY2=17.5
C
C     ----- Major radius dependence -----
C
      NTGX=NTGMAX/2+1
      DO NSG=NSGMAX,1,-1
         NX=NSGMAX-NSG+1
         XAX(NX)=RMG(NSG,NTGX)
         RPSINT(NX,1)=PSIO(NTGX,NSG,1)
         RPSINT(NX,2)=PSIO(NTGX,NSG,2)
         RPSINT(NX,3)=PSIO(NTGX,NSG,3)
         RPSINT(NX,4)=PSIO(NTGX,NSG,4)
         RPSINT(NX,5)=PSIO(NTGX,NSG,5)
      ENDDO
      NTGX=1
      DO NSG=1,NSGMAX
         NX=NSGMAX+NSG
         XAX(NX)=RMG(NSG,NTGX)
         RPSINT(NX,1)=PSIO(NTGX,NSG,1)
         RPSINT(NX,2)=PSIO(NTGX,NSG,2)
         RPSINT(NX,3)=PSIO(NTGX,NSG,3)
         RPSINT(NX,4)=PSIO(NTGX,NSG,4)
         RPSINT(NX,5)=PSIO(NTGX,NSG,5)
      ENDDO
C
      DO NX=1,2*NSGMAX
         GX(NX)=GUCLIP(XAX(NX))
         GYPS(NX,1)=GUCLIP(RPSINT(NX,1))
         GYPS(NX,2)=GUCLIP(RPSINT(NX,2))
         GYPS(NX,3)=GUCLIP(RPSINT(NX,3))
         GYPS(NX,4)=GUCLIP(RPSINT(NX,4))
         GYPS(NX,5)=GUCLIP(RPSINT(NX,5))
      ENDDO
C
      CALL PAGES
      KSTR='/PSI(X)/'
      CALL EQGR1D(GX1,GX2,GY1,GY2,GX,GYPS(1,1),NXM,2*NSGMAX,5,KSTR,0)
      CALL PAGEE
C
      RETURN
      END
C
