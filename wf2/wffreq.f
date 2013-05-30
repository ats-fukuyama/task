C     $Id$
C
C     ***** PLOT 1D GRAPH ******************************************
C
      SUBROUTINE WFFREQ
C
      INCLUDE 'wfcomm.inc'
      PARAMETER (NNXM=201)
C
      DIMENSION WP(NSM),WC(NSM)
      DIMENSION GXL(NNXM),GYL(NNXM,8)
C
      PN1=13.D0
      PN2=20.D0
      NPNMAX=101
      RPARA=0.3D0
C
    1 WRITE(6,*)'## SELECT : (1)WUH (2)WLH (3)WCR (4)WCL (5)ALL (0)END'
      READ(5,*,ERR=1,END=9000) ISW
      IF(ISW.EQ.0) THEN
         GOTO 9000
      ELSEIF(ISW.EQ.1.OR.
     &       ISW.EQ.2.OR.
     &       ISW.EQ.3.OR.
     &       ISW.EQ.4.OR.
     &       ISW.EQ.5    ) THEN
         GOTO 2
      ELSE
         GOTO 1
      ENDIF
C
  100 CONTINUE
C
    2 WRITE(6,*) ' INPUT : PN1,PN2,NPNMAX,RPARA'
      READ(5,*,ERR=2,END=1) PN1,PN2,NPNMAX,RPARA
      IF(NPNMAX.LE.1) GOTO 1
      DPN=(PN2-PN1)/(NPNMAX-1)
      DO NX=1,NPNMAX
           PNL=PN1+DPN*(NX-1)
           PNA=10.D0**PNL
         DO IS=1,2
            WP(IS)=SQRT(PNA*PZ(IS)*PZ(IS)*AEE*AEE/(PA(IS)*AMP*EPS0))
            WC(IS)=PZ(IS)*AEE*BB/(PA(IS)*AMP)
         ENDDO
C
      IF(ISW.EQ.1) THEN
C
C       WUH : Upper hybrid resonance frequency [Hz]
C
            WUH2=WC(1)*WC(1)+WP(1)*WP(1)
            WUH =SQRT(WUH2)/(2.D0*PI)
            GXL(NX)=GCLIP(PNL)
            GYL(NX,1)=GCLIP(LOG10(WUH))
C
      ELSEIF(ISW.EQ.2) THEN
C
C       WLH : Lower hybrid resonance frequency [Hz]
C
            WLH2=ABS(WC(1))*WC(2)*(WP(1)*WP(1)+ABS(WC(1))*WC(2))
     &          /(WP(1)*WP(1)+WC(1)*WC(1))
            WLH =SQRT(WLH2)/(2.D0*PI)
            GXL(NX)=GCLIP(PNL)
            GYL(NX,2)=GCLIP(LOG10(WLH))
C
      ELSEIF(ISW.EQ.3) THEN
C
C       WCR : Cut-off frequency of Right-hand polarization [Hz]
C
            WCR=(SQRT(WC(1)*WC(1)+4.D0*WP(1)*WP(1))+ABS(WC(1))-WC(2))
     &         /(4.D0*PI)
            GXL(NX)=GCLIP(PNL)
            GYL(NX,3)=GCLIP(LOG10(WCR))
C
      ELSEIF(ISW.EQ.4) THEN
C
C       WCL : Cut-off frequency of Left-hand polarization [Hz]
C
            WCL=(SQRT(WC(1)*WC(1)+4.D0*WP(1)*WP(1))-ABS(WC(1))+WC(2))
     &         /(4.D0*PI)
            GXL(NX)=GCLIP(PNL)
            GYL(NX,4)=GCLIP(LOG10(WCL))
      ELSEIF(ISW.EQ.5) THEN
            WCE=ABS(WC(1))
            WPE=WP(1)
            WCI=WC(2)
            WPI=WP(2)
C
            WCE2=WCE*WCE
            WPE2=WPE*WPE
            WCI2=WCI*WCI
            WPI2=WPI*WPI
C
            VB = WCE2+WCI2+WPE2+WPI2
            VS = WCE2*WCI2+WPE2*WCI2+WPI2*WCE2
            VD = SQRT(VB*VB-4.D0*VS)
            WUH2=0.5D0*(VB+VD)
            WLH2=0.5D0*(VB-VD)
C
            WD1=SQRT((WCE+WCI)**2+4.D0*(WPE2+WPI2))
            WCR=0.5D0*(WD1+WCE-WCI)
            WCL=0.5D0*(WD1-WCE+WCI)
C
            RKC=VC*PI/RPARA
            RKC2=RKC**2
            VA=1.D0+(WPE2+WPI2)/RKC2
            VB=-WCE+WCI
            VS=-WCE*WCI
            VD=VB*VB-4.D0*VA*VS
            IF(VD.GT.0.D0) THEN
               WHL=(-VB+SQRT(VD))/(2.D0*VA)
            ELSE
               WHL=-VB/(2.D0*VA)
            ENDIF
            WHL=ABS(WHL)
            WRITE(6,*) PNL,RKC,WPE
C
            WUH =SQRT(WUH2)/(2.D0*PI)
            WLH =SQRT(WLH2)/(2.D0*PI)
            WCR =WCR/(2.D0*PI)
            WCL =WCL/(2.D0*PI)
            WCE =WCE/(2.D0*PI)
            WCI =WCI/(2.D0*PI)
            WPE =WPE/(2.D0*PI)
            WHL =WHL/(2.D0*PI)
C
            GXL(NX)=GCLIP(PNL)
            GYL(NX,1)=GCLIP(LOG10(WUH))
            GYL(NX,2)=GCLIP(LOG10(WLH))
            GYL(NX,3)=GCLIP(LOG10(WCR))
            GYL(NX,4)=GCLIP(LOG10(WCL))
            GYL(NX,5)=GCLIP(LOG10(WCE))
            GYL(NX,6)=GCLIP(LOG10(WCI))
            GYL(NX,7)=GCLIP(LOG10(WPE))
            GYL(NX,8)=GCLIP(LOG10(WHL))
      ELSE
         GOTO 1
      ENDIF
      ENDDO
C
      CALL PAGES
      CALL SETLIN(0,0,7)
      CALL SETCHS(0.3,0.)
      CALL GMNMX1(GXL,1,NPNMAX,1,GXMIN,GXMAX)
      CALL GMNMX2(GYL,NNXM,1,NPNMAX,1,1,4,1,GYMIN,GYMAX)
      CALL GQSCAL(GXMIN,GXMAX,GXSMN,GXSMX,GSCALX)
      CALL GQSCAL(GYMIN,GYMAX,GYSMN,GYSMX,GSCALY)
      GYSMN=5.0
      GYSMX=10.0
      CALL GDEFIN(4.,18.,3.,17.,GXSMN,GXSMX,GYSMN,GYSMX)
      CALL GFRAME
      CALL GSCALL(GXSMN,3,0.0,0,0.2,9)
      CALL GSCALL(0.0,0,GYSMN,9,0.2,9)
      CALL GVALUL(GXSMN,1,0.0,1,0)
      CALL GVALUL(0.0,0,GYSMN,1,0)
      IF(ISW.EQ.1) THEN
        CALL GPLOTP(GXL,GYL(1,1),1,NPNMAX,1,0,0,0)
      ELSEIF(ISW.EQ.2) THEN
        CALL GPLOTP(GXL,GYL(1,2),1,NPNMAX,1,0,0,0)
      ELSEIF(ISW.EQ.3) THEN
        CALL GPLOTP(GXL,GYL(1,3),1,NPNMAX,1,0,0,0)
      ELSEIF(ISW.EQ.4) THEN
        CALL GPLOTP(GXL,GYL(1,4),1,NPNMAX,1,0,0,0)
      ELSEIF(ISW.EQ.5) THEN
        CALL SETLIN(4,0,3)
        CALL GPLOTP(GXL,GYL(1,1),1,NPNMAX,1,0,0,0)
        CALL SETLIN(0,0,4)
        CALL GPLOTP(GXL,GYL(1,2),1,NPNMAX,1,0,0,0)
        CALL SETLIN(1,0,5)
        CALL GPLOTP(GXL,GYL(1,3),1,NPNMAX,1,0,0,0)
        CALL SETLIN(2,0,6)
        CALL GPLOTP(GXL,GYL(1,4),1,NPNMAX,1,0,0,0)
        CALL SETLIN(2,0,7)
        CALL GPLOTP(GXL,GYL(1,5),1,NPNMAX,1,0,0,0)
        CALL GPLOTP(GXL,GYL(1,6),1,NPNMAX,1,0,0,0)
        CALL GPLOTP(GXL,GYL(1,7),1,NPNMAX,1,0,0,0)
        CALL GPLOTP(GXL,GYL(1,8),1,NPNMAX,1,0,0,0)
      ENDIF
C
      CALL SETLIN(0,0,7)
      CALL MOVE(10.0,18.0)
      IF (ISW.EQ.1) THEN
        CALL MOVE(2.0,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' WUH [Hz] ',10) 
         CALL SETCHS(0.3,0.0)
      ELSEIF (ISW.EQ.2) THEN
        CALL MOVE(2.0,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' WLH [Hz] ',10) 
         CALL SETCHS(0.3,0.0)
      ELSEIF (ISW.EQ.3) THEN
        CALL MOVE(2.0,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' WCR [Hz] ',10) 
         CALL SETCHS(0.3,0.0)
      ELSEIF (ISW.EQ.4) THEN
        CALL MOVE(2.0,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' WCL [Hz] ',10) 
         CALL SETCHS(0.3,0.0)
      ELSEIF (ISW.EQ.5) THEN
        CALL MOVE(2.0,9.0)
         CALL SETCHS(0.3,90.0)
         CALL TEXT (' f [Hz] ',8) 
         CALL SETCHS(0.3,0.0)
      ENDIF
C
      CALL MOVE(10.0,1.5)
      CALL TEXT (' PNA [/m^3]',11)
C
      CALL SETLIN(0,0,3)
      CALL MOVE(20.0,16.5)
      CALL TEXT('WUH =',5)
      CALL NUMBD(WUH,  '(1PE11.3)',11)
      CALL SETLIN(0,0,4)
      CALL MOVE(20.0,16.0)
      CALL TEXT('WLH =',5)
      CALL NUMBD(WLH, '(1PE11.3)',11)
      CALL SETLIN(0,0,5)
      CALL MOVE(20.0,15.5)
      CALL TEXT('WCR =',5)
      CALL NUMBD(WCR, '(1PE11.3)',11)
      CALL SETLIN(0,0,6)
      CALL MOVE(20.0,15.0)
      CALL TEXT('WCL =',5)
      CALL NUMBD(WCL, '(1PE11.3)',11)
      CALL SETLIN(0,0,7)
      CALL MOVE(20.0,14.5)
      CALL TEXT('PNA =',5)
      CALL NUMBD(PNA,'(1PE11.3)',11)
      CALL MOVE(20.0,13.5)
      CALL TEXT(' BB =',5)
      CALL NUMBD(BB,'(1PE11.3)',11)
      CALL MOVE(20.0,12.5)
      CALL TEXT(' RP =',5)
      CALL NUMBD(RPARA,'(1PE11.3)',11)
      CALL PAGEE
      GOTO 100
C
 9000 RETURN
      END
