      IMPLICIT REAL*8 (A-F,H,O-Z)
      PARAMETER (NXM=65,NZM=65)
      DIMENSION GX(NXM),GY(NXM,2),GJY(NXM),GYC(NXM,2)
      DIMENSION RX(NXM),ZX(NZM),PSI(NXM,NZM),PSIN(NXM),HJ(NXM)
      DIMENSION PSIZAX(NXM),PSIZC(NXM)
      CHARACTER TXT*80
C
      CALL GSOPEN
C
      NXMAX = 64
      NZMAX = 64 
C
C     *********************************************
C     *            Given Values                   *
C     *********************************************
C
C     *** CONSTANTS ***
C
C        PI    : Pi
C        RMU0  : Permeability of free space
C        BLTZ  : Boltzmann constant
C        AMP   : Proton mass
C        RGAS  : Gas constant
C
      PI     = 2.D0*ASIN(1.D0)
      RMU0   = 4.D0*PI*1.D-7
      BLTZ   = 1.38066D-23
      AMP    = 1.6726231D-27
      RGAS   = BLTZ/AMP
C
C     *** CONFIGURATION PARAMETERS ***
C
C        RR    : Plasma major radius                             (m)
C        RA    : Plasma minor radius                             (m)
C        EPSA  : Cross-section parameter ; -1 < EPSA < 1
C        RKAP  : Plasma shape elongation
C        ZE    : Plasma range of z axis                          (m)
C
      RR     = 3.D0
      RA     = 1.D0
      EPSA   = 0.D0
      RKAP   = 1.6D0
      ZE     = RA*RKAP
C
C     *** NUBERICAL SOLUTIONS etc.***
C
      HM     = 2.1D1*1.D6
      PSI0   = 2.4052573068990517D-3
      PP1    = 7.38D2*1.D6
C      PP1    = 7.55D2*1.D6
      OTC    = 1.5D-1
      RAXIS  = 3.211204904635547D0
      ZAXIS  = 3.8654003453737973D-12
C
C     *** CHANGE PARAMETER ***
C
      NCHG   = 1
C
C     *********************************************
C     *            ANALYTIC PSI,J                 *
C     *********************************************
C
      DX = (2.D0*RA)/(NXMAX-1)
      DZ = (2.D0*ZE)/(NZMAX-1)
      CC=((EPSA-1.D0)*RAXIS**2/(8.D0*RR**2))+(RGAS/(2.D0*RR**2*OTC))
     &     *(EXP(OTC*RAXIS**2/(2.D0*RGAS))-1.D0)
C
      DO NX=1,NXMAX
      DO NZ=1,NZMAX
         IF (NX.LE.NXMAX/2) THEN
            RX(NX)=(RR-RA)+DX*DBLE(NX-1)
         ELSE
            RX(NX)=(RR+DBLE(DX/2))+DX*DBLE(NX-NXMAX/2-1)
         ENDIF
         IF (NX.LE.NZMAX/2) THEN
            ZX(NZ)=-DZ*(DBLE(NZMAX/2-1)+0.5D0)+DZ*DBLE(NZ-1)
         ELSE
            ZX(NZ)=DBLE(DZ/2)+DZ*DBLE(NZ-NZMAX/2-1)
         ENDIF
         PSI(NX,NZ)=CC*PP1*RX(NX)**2*RR**2
     &        -(HM*ZX(NZ)**2/(2.D0*RR**2))
     &        +PP1*((EPSA-1.D0)/4.D0)*(ZX(NZ)**2/RR**2
     &        -RX(NX)**2/(4.D0*RR**2))*RX(NX)**2*RR**2
     &        +(RGAS**2/OTC**2)*PP1*(1.D0+(OTC*RX(NX)**2/(2.D0*RGAS))
     &        -EXP(OTC*RX(NX)**2/(2.D0*RGAS)))
      ENDDO
         PSIZAX(NX)=CC*PP1*RX(NX)**2*RR**2
     &        -(HM*ZAXIS**2/(2.D0*RR**2))
     &        +PP1*((EPSA-1.D0)/4.D0)*(ZAXIS**2/RR**2
     &        -RX(NX)**2/(4.D0*RR**2))*RX(NX)**2*RR**2
     &        +(RGAS**2/OTC**2)*PP1*(1.D0+(OTC*RX(NX)**2/(2.D0*RGAS))
     &        -EXP(OTC*RX(NX)**2/(2.D0*RGAS)))
         HJ(NX)=(HM/(RR**2*RX(NX))+RX(NX)*PP1*RMU0
     &           *EXP(RX(NX)**2*OTC/(2.D0*RGAS)))*1.D-6
      ENDDO
C
C     *********************************************
C     *            NUMERICAL PSI                  *
C     *********************************************
C
      PSIN(1)  = -9.5943710774369713D-5
      PSIN(2)  = -2.7781853147859463D-4
      PSIN(3)  = -4.4987665092485276D-4
      PSIN(4)  = -6.1232317983507167D-4
      PSIN(5)  = -7.6537226430740732D-4
      PSIN(6)  = -9.0924692169933881D-4
      PSIN(7)  = -1.0441788874719422D-3
      PSIN(8)  = -1.1704084723640209D-3
      PSIN(9)  = -1.2881844293105824D-3
      PSIN(10) = -1.3977638295624946D-3
      PSIN(11) = -1.4994119475035205D-3
      PSIN(12) = -1.5934021536977157D-3
      PSIN(13) = -1.6800158157349323D-3
      PSIN(14) = -1.7595422064754263D-3
      PSIN(15) = -1.8322784193269666D-3
      PSIN(16) = -1.8985292902203341D-3
      PSIN(17) = -1.9586073259828074D-3
      PSIN(18) = -2.0128326388460667D-3
      PSIN(19) = -2.0615328868676790D-3
      PSIN(20) = -2.1050432200985472D-3
      PSIN(21) = -2.1437062324002696D-3
      PSIN(22) = -2.1778719189200739D-3
      PSIN(23) = -2.2078976393925167D-3
      PSIN(24) = -2.2341480877064495D-3
      PSIN(25) = -2.2569952686557486D-3
      PSIN(26) = -2.2768184837097474D-3
      PSIN(27) = -2.2940043295359409D-3
      PSIN(28) = -2.3089467173143237D-3
      PSIN(29) = -2.3220469319012516D-3
      PSIN(30) = -2.3337137830265524D-3
      PSIN(31) = -2.3443640269135653D-3
      PSIN(32) = -2.3544239539133178D-3
      PSIN(33) = -2.3643985387783673D-3
      PSIN(34) = -2.3742213958465248D-3
      PSIN(35) = -2.3834158728508212D-3
      PSIN(36) = -2.3915353301910693D-3
      PSIN(37) = -2.3981234959896489D-3
      PSIN(38) = -2.4027161136010049D-3
      PSIN(39) = -2.4048413862981803D-3
      PSIN(40) = -2.4040201414376567D-3
      PSIN(41) = -2.3997658993229019D-3
      PSIN(42) = -2.3915849015107735D-3
      PSIN(43) = -2.3789761188577589D-3
      PSIN(44) = -2.3614312480809536D-3
      PSIN(45) = -2.3384347010830666D-3
      PSIN(46) = -2.3094635892908437D-3
      PSIN(47) = -2.2739877042886230D-3
      PSIN(48) = -2.2314694955250074D-3
      PSIN(49) = -2.1813640455922889D-3
      PSIN(50) = -2.1231190434163461D-3
      PSIN(51) = -2.0561747555963223D-3
      PSIN(52) = -1.9799639960711936D-3
      PSIN(53) = -1.8939120942495311D-3
      PSIN(54) = -1.7974368617112080D-3
      PSIN(55) = -1.6899485575705255D-3
      PSIN(56) = -1.5708498525763877D-3
      PSIN(57) = -1.4395357920149266D-3
      PSIN(58) = -1.2953937574722421D-3
      PSIN(59) = -1.1378034275088542D-3
      PSIN(60) = -9.6613673729269269D-4
      PSIN(61) = -7.7975783723351704D-4
      PSIN(62) = -5.7802305065841413D-4
      PSIN(63) = -3.6028083056526983D-4
      PSIN(64) = -1.2587171548875036D-4
C
C     *********************************************
C     *             CORRECTING PSI                *
C     *********************************************
C
      DELP=(PSIZAX(1)*RMU0*1.D-6+PSI0)+PSIN(1)
      DO NX=1,NXMAX
         PSIZC(NX)=(((PSIZAX(NX)*RMU0*1.D-6+PSI0)-DELP)-PSI0)
     &            /(RMU0*1.D-6)
      ENDDO
C
C     *********************************************
C     *               DRAW GRAPH                  *
C     *********************************************
C
      DO 100 N=1,NXMAX
        GX(N)   = GUCLIP(RX(N))
        IF (NCHG.EQ.0) THEN
           GY(N,1) = GUCLIP(-(PSI(N,NZMAX/2)*RMU0*1.D-6+PSI0))
        ELSEIF (NCHG.EQ.1) THEN
           GY(N,1)  = GUCLIP(-(PSIZAX(N)*RMU0*1.D-6+PSI0))
           GYC(N,1) = GUCLIP(-(PSIZC(N)*RMU0*1.D-6+PSI0))
           GYC(N,2) = GUCLIP(PSIN(N))
        ENDIF
        GY(N,2) = GUCLIP(PSIN(N))
        GJY(N)  = GUCLIP(HJ(N))
  100 CONTINUE
C
      CALL PAGES
      TXT='/ANALYTIC AND NUMERICAL PSI/'
      CALL GRAPH1(3.5,23.5,2.0,17.0,GX,GY,NXM,NXMAX,2,TXT,7.5)
      CALL PAGEE
C
      IF (NCHG.EQ.1) THEN
      CALL PAGES
      TXT='/ANALYTIC AND NUMERICAL PSI(CORRECTED)/'
      CALL GRAPH1(3.5,23.5,2.0,17.0,GX,GYC,NXM,NXMAX,2,TXT,9.0)
      CALL PAGEE
      ENDIF
C
      CALL PAGES
      TXT='/ANALYTIC J/'
      CALL GRAPH1(3.5,23.5,2.0,17.0,GX,GJY,NXM,NXMAX,1,TXT,3.5)
      CALL PAGEE
      GOTO 9000
C
 9000 CALL GSCLOS
      STOP
      END
C
C     *********************************************
C     *               SUBROUTINES                 *
C     *********************************************
C
      SUBROUTINE GRAPH1(PXMIN,PXMAX,PYMIN,PYMAX,GX,GY,NXM,
     &                   NXMAX,NGMAX,TXT,POS)
C
      DIMENSION GX(NXM),GY(NXM,NGMAX)
      CHARACTER*80 TXT
C
      CALL GMNMX1(GX,1,NXMAX,1,XMIN,XMAX)
      CALL GMNMX1(GY(1,1),1,NXMAX,1,YMIN,YMAX)
      DO 100 NG=2,NGMAX
         CALL GMNMX1(GY(1,NG),1,NXMAX,1,YMIN1,YMAX1)
         YMIN=MIN(YMIN,YMIN1)
         YMAX=MAX(YMAX,YMAX1)
  100 CONTINUE
      CALL GQSCAL(XMIN,XMAX,GXMIN,GXMAX,GXSCAL)
      CALL GQSCAL(YMIN,YMAX,GYMIN,GYMAX,GYSCAL)
C
C     Origin should have a scale mark
C
      IF(GXMIN*GXMAX.GT.0.0) THEN
         GXORG=GXMIN
      ELSE
         GXORG=0.0
      ENDIF
      IF(GYMIN*GYMAX.GT.0.0) THEN
         GYORG=GYMIN
      ELSE
         GYORG=0.0
      ENDIF
C
C     symbol font has better minus sign
C
      CALL SETFNT(2)
      CALL SETCHS(0.35,0.0)
      CALL GDEFIN(PXMIN,PXMAX,PYMIN,PYMAX,GXMIN,GXMAX,GYMIN,GYMAX)
      CALL SETLIN(0,0,4)
      CALL GFRAME
C
C     NGULEN choose appropriate format of scale values
C
      CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.3,9)
      CALL GVALUE(GXORG,2*GXSCAL,0.0,0.0,NGULEN(2*GXSCAL))
      CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.3,9)
      CALL GVALUE(0.0,0.0,GYORG,2*GYSCAL,NGULEN(2*GYSCAL))
      DO 200 NG=1,NGMAX
         CALL SETLIN(-1,-1,7-MOD(NG-1,5))
         CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,0)
C         CALL GPLOTP(GX,GY(1,NG),1,NXMAX,1,0,0,MOD(NG-1,8))
  200 CONTINUE
      CALL SETLIN(-1,-1,7)
C
      CALL GTEXTX(POS,17.3,TXT,2)
C
      RETURN
      END
