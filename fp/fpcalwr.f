C     $Id$
C
C ************************************************************
C
C            CALCULATION OF DW (BOUNCE AVERAGED,RAY)
C
C ************************************************************
C
      SUBROUTINE FPCALWR(NSA)
C
      INCLUDE 'fpcomm.inc'
      INCLUDE '../wr/wrcom1.inc'
      INCLUDE 'fpcom2.inc'
C
      DIMENSION DLA(0:NITM,NRAYMX)
C
C =============  CALCULATION OF DWPP AND DWPT  ===============
C
      FACT=0.5D0
C
      IF(MODELW.EQ.2) THEN
      DO  NRDO=1,NRMAX
         NR=NRDO
         DO  NTH=1,NTHMAX
C
         DELH=4.D0*ETAM(NTH,NR)/NAVMAX
C
         DO NAV=1,NAVMAX
            ETAL=DELH*(NAV-0.5D0)-2.D0*ETAM(NTH,NR)
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
               NITMX=NITMAX(NRAY)
               RFDW=RAYIN(1,NRAY)
C
                  DO NIT=0,NITMX
                     RXB=RXS(NIT,NRAY)
                     RYB=RYS(NIT,NRAY)
                     RZB=RZS(NIT,NRAY)
                     RRLB=SQRT(RXB**2+RYB**2)
                     RZLB=RZB
                     DLA(NIT,NRAY)=SQRT((RRLB-RL)**2+(RZLB-ZL)**2)
                  ENDDO
C
                  MINNB1=0
                  DO NIT=0,NITMX-1
                     IF(DLA(NIT+1,NRAY).LT.DLA(NIT,NRAY))THEN
                        MINNB1=NIT+1
                     ENDIF
                  ENDDO
C
                  IF(MINNB1.EQ.0) GOTO 1
                  IF(MINNB1.EQ.NITMX) GOTO 1
C
                  IF(DLA(MINNB1-1,NRAY).LT.DLA(MINNB1+1,NRAY))THEN
                     MINNB2=MINNB1-1
                  ELSE
                     MINNB2=MINNB1+1
                  ENDIF
C
                  DLAMN1=   DLA(MINNB1,NRAY)
                  RXMIN1=   RXS(MINNB1,NRAY)
                  RYMIN1=   RYS(MINNB1,NRAY)
                  RZMIN1=   RZS(MINNB1,NRAY)
                  CEXMN1=  CEXS(MINNB1,NRAY)
                  CEYMN1=  CEYS(MINNB1,NRAY)
                  CEZMN1=  CEZS(MINNB1,NRAY)
                  RKXMN1=  RKXS(MINNB1,NRAY)
                  RKYMN1=  RKYS(MINNB1,NRAY)
                  RKZMN1=  RKZS(MINNB1,NRAY)
                  RBMIN1=RAYRB1(MINNB1,NRAY)
C
                  DLAMN2=   DLA(MINNB2,NRAY)
                  RXMIN2=   RXS(MINNB2,NRAY)
                  RYMIN2=   RYS(MINNB2,NRAY)
                  RZMIN2=   RZS(MINNB2,NRAY)
                  CEXMN2=  CEXS(MINNB2,NRAY)
                  CEYMN2=  CEYS(MINNB2,NRAY)
                  CEZMN2=  CEZS(MINNB2,NRAY)
                  RKXMN2=  RKXS(MINNB2,NRAY)
                  RKYMN2=  RKYS(MINNB2,NRAY)
                  RKZMN2=  RKZS(MINNB2,NRAY)
                  RBMIN2=RAYRB1(MINNB2,NRAY)
C               
C
                  RRLMN1=SQRT(RXMIN1**2+RYMIN1**2)
                  RZLMN1=RZMIN1
                  RRLMN2=SQRT(RXMIN2**2+RYMIN2**2)
                  RZLMN2=RZMIN2
                  DEL12=SQRT((RRLMN1-RRLMN2)**2+(RZLMN1-RZLMN2)**2)
C
                  XA1=(DLAMN1**2-DLAMN2**2+DEL12**2)/(2.D0*DEL12)
                  XA2=(DLAMN2**2-DLAMN1**2+DEL12**2)/(2.D0*DEL12)             
                  XLL2=DLAMN1**2-XA1**2
                  A1=XA1/(XA1+XA2)
                  A2=XA2/(XA1+XA2)
C
                  CEX =A1*CEXMN2+A2*CEXMN1
                  CEY =A1*CEYMN2+A2*CEYMN1
                  CEZ =A1*CEZMN2+A2*CEZMN1
                  RKX =A1*RKXMN2+A2*RKXMN1
                  RKY =A1*RKYMN2+A2*RKYMN1
                  RKZ =A1*RKZMN2+A2*RKZMN1
                  RXB =A1*RXMIN2+A2*RXMIN1
                  RYB =A1*RYMIN2+A2*RYMIN1
                  RZB =A1*RZMIN2+A2*RZMIN1
                  RADB=A1*RBMIN2+A2*RBMIN1
                  IF(RADB.NE.0.D0) THEN
                     DELYEC=RADB
                  ENDIF
                  DELCR2=XLL2
                  DELRB2=DELYEC**2
                  ARG=DELCR2/DELRB2
C               WRITE(25,*) 'NAV,NCR,DELR2,ARG=',NAV,NCR,DELR2,ARG
                     ARGB (NR,NTH,NAV,NRAY)=ARG
                     CEB(1,NR,NTH,NAV,NRAY)=CEX
                     CEB(2,NR,NTH,NAV,NRAY)=CEY
                     CEB(3,NR,NTH,NAV,NRAY)=CEZ
                     RKB(1,NR,NTH,NAV,NRAY)=RKX
                     RKB(2,NR,NTH,NAV,NRAY)=RKY
                     RKB(3,NR,NTH,NAV,NRAY)=RKZ
                     RBB(1,NR,NTH,NAV,NRAY)=RXB
                     RBB(2,NR,NTH,NAV,NRAY)=RYB
                     RBB(3,NR,NTH,NAV,NRAY)=RZB
 1                CONTINUE
C
C            
C            IF(IFLAG.EQ.1) THEN
C               WRITE(6,'(3I3)') NR,NAV,NCR
C               WRITE(6,'(1P3E12.4)') X,Y,Z
C               WRITE(6,'(1P3E12.4)') RX,RY,RZ
C               WRITE(6,'(1P3E12.4)') DELR2,DELCR2,ARG
C            ENDIF
C
          ENDDO
        ENDDO
       ENDDO
      ENDDO
      ENDIF
C
      DO 1000 NRDO=1,NRMAX
         NR=NRDO
         DO 101 NTH=1,NTHMAX
            IF(NTH.EQ.ITL(NR).OR.NTH.EQ.ITU(NR)) GOTO 101
         DO 100 NP=1,NPMAX+1
C            IF(NP.EQ.10.AND.NTH.EQ.2) THEN
C               IFLAG=1 
C            ELSE
C               IFLAG=0
C            ENDIF
            CALL FPDWAV(ETAM(NTH,NR),SINM(NTH),COSM(NTH),PG(NP),NR,NTH,
     &                  DWPPS,DWPTS,DWTPS,DWTTS,NSA)
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
     &                            +DWPP(NTHMAX-NTH+1,NP,NR,NSA))*FACT
            DWPT(NTH,NP,NR,NSA)  =(DWPT(NTH,NP,NR,NSA)
     &                            +DWPT(NTHMAX-NTH+1,NP,NR,NSA))*FACT
            DWPP(NTHMAX-NTH+1,NP,NR,NSA)  =DWPP(NTH,NP,NR,NSA)
            DWPT(NTHMAX-NTH+1,NP,NR,NSA)  =DWPT(NTH,NP,NR,NSA)
  190    CONTINUE
         DWPP(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWPP(ITL(NR)-1,NP,NR,NSA)
     &                                         /RLAMDA(ITL(NR)-1,NR)
     &                    +DWPP(ITL(NR)+1,NP,NR,NSA)
     &                                         /RLAMDA(ITL(NR)+1,NR)
     &                    +DWPP(ITU(NR)-1,NP,NR,NSA)
     &                                         /RLAMDA(ITU(NR)-1,NR)
     &                    +DWPP(ITU(NR)+1,NP,NR,NSA)
     &                                         /RLAMDA(ITU(NR)+1,NR))
C
         DWPT(ITL(NR),NP,NR,NSA)=RLAMDA(ITL(NR),NR)/4.D0
     &                  *( DWPT(ITL(NR)-1,NP,NR,NSA)
     &                                         /RLAMDA(ITL(NR)-1,NR)
     &                    +DWPT(ITL(NR)+1,NP,NR,NSA)
     &                                         /RLAMDA(ITL(NR)+1,NR)
     &                    +DWPT(ITU(NR)-1,NP,NR,NSA)
     &                                         /RLAMDA(ITU(NR)-1,NR)
     &                    +DWPT(ITU(NR)+1,NP,NR,NSA)
     &                                         /RLAMDA(ITU(NR)+1,NR))
         DWPP(ITU(NR),NP,NR,NSA)  =DWPP(ITL(NR),NP,NR,NSA)
         DWPT(ITU(NR),NP,NR,NSA)  =DWPT(ITL(NR),NP,NR,NSA)
  200    CONTINUE
         ENDIF
 1000 CONTINUE
C
C =============  CALCULATION OF DWTP AND DWTT  ===============
C
      IF(MODELW.EQ.2) THEN
      DO  NRDO=1,NRMAX
         NR=NRDO
         DO  NTH=1,NTHMAX
C
         DELH=4.D0*ETAG(NTH,NR)/NAVMAX
C
         DO NAV=1,NAVMAX
            ETAL=DELH*(NAV-0.5D0)-2.D0*ETAG(NTH,NR)
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
               NITMX=NITMAX(NRAY)
               RFDW=RAYIN(1,NRAY)
C
                  DO NIT=0,NITMX
                     RXB=RXS(NIT,NRAY)
                     RYB=RYS(NIT,NRAY)
                     RZB=RZS(NIT,NRAY)
                     RRLB=SQRT(RXB**2+RYB**2)
                     RZLB=RZB
                     DLA(NIT,NRAY)=SQRT((RRLB-RL)**2+(RZLB-ZL)**2)
                  ENDDO
C
                  MINNB1=0
                  DO NIT=0,NITMX-1
                     IF(DLA(NIT+1,NRAY).LT.DLA(NIT,NRAY))THEN
                        MINNB1=NIT+1
                     ENDIF
                  ENDDO
C
                  IF(MINNB1.EQ.0) GOTO 2
                  IF(MINNB1.EQ.NITMX) GOTO 2
C
                  IF(DLA(MINNB1-1,NRAY).LT.DLA(MINNB1+1,NRAY))THEN
                     MINNB2=MINNB1-1
                  ELSE
                     MINNB2=MINNB1+1
                  ENDIF
C
                  DLAMN1=   DLA(MINNB1,NRAY)
                  RXMIN1=   RXS(MINNB1,NRAY)
                  RYMIN1=   RYS(MINNB1,NRAY)
                  RZMIN1=   RZS(MINNB1,NRAY)
                  CEXMN1=  CEXS(MINNB1,NRAY)
                  CEYMN1=  CEYS(MINNB1,NRAY)
                  CEZMN1=  CEZS(MINNB1,NRAY)
                  RKXMN1=  RKXS(MINNB1,NRAY)
                  RKYMN1=  RKYS(MINNB1,NRAY)
                  RKZMN1=  RKZS(MINNB1,NRAY)
                  RBMIN1=RAYRB1(MINNB1,NRAY)
C
                  DLAMN2=   DLA(MINNB2,NRAY)
                  RXMIN2=   RXS(MINNB2,NRAY)
                  RYMIN2=   RYS(MINNB2,NRAY)
                  RZMIN2=   RZS(MINNB2,NRAY)
                  CEXMN2=  CEXS(MINNB2,NRAY)
                  CEYMN2=  CEYS(MINNB2,NRAY)
                  CEZMN2=  CEZS(MINNB2,NRAY)
                  RKXMN2=  RKXS(MINNB2,NRAY)
                  RKYMN2=  RKYS(MINNB2,NRAY)
                  RKZMN2=  RKZS(MINNB2,NRAY)
                  RBMIN2=RAYRB1(MINNB2,NRAY)
C               
C
                  RRLMN1=SQRT(RXMIN1**2+RYMIN1**2)
                  RZLMN1=RZMIN1
                  RRLMN2=SQRT(RXMIN2**2+RYMIN2**2)
                  RZLMN2=RZMIN2
                  DEL12=SQRT((RRLMN1-RRLMN2)**2+(RZLMN1-RZLMN2)**2)
C
                  XA1=(DLAMN1**2-DLAMN2**2+DEL12**2)/(2.D0*DEL12)
                  XA2=(DLAMN2**2-DLAMN1**2+DEL12**2)/(2.D0*DEL12)             
                  XLL2=DLAMN1**2-XA1**2
                  A1=XA1/(XA1+XA2)
                  A2=XA2/(XA1+XA2)
C
                  CEX =A1*CEXMN2+A2*CEXMN1
                  CEY =A1*CEYMN2+A2*CEYMN1
                  CEZ =A1*CEZMN2+A2*CEZMN1
                  RKX =A1*RKXMN2+A2*RKXMN1
                  RKY =A1*RKYMN2+A2*RKYMN1
                  RKZ =A1*RKZMN2+A2*RKZMN1
                  RXB =A1*RXMIN2+A2*RXMIN1
                  RYB =A1*RYMIN2+A2*RYMIN1
                  RZB =A1*RZMIN2+A2*RZMIN1
                  RADB=A1*RBMIN2+A2*RBMIN1
                  IF(RADB.NE.0.D0) THEN
                     DELYEC=RADB
                  ENDIF
C
                  DELCR2=XLL2
                  DELRB2=DELYEC**2
                  ARG=DELCR2/DELRB2
C               WRITE(25,*) 'NAV,NCR,DELR2,ARG=',NAV,NCR,DELR2,ARG
                     ARGB (NR,NTH,NAV,NRAY)=ARG
                     CEB(1,NR,NTH,NAV,NRAY)=CEX
                     CEB(2,NR,NTH,NAV,NRAY)=CEY
                     CEB(3,NR,NTH,NAV,NRAY)=CEZ
                     RKB(1,NR,NTH,NAV,NRAY)=RKX
                     RKB(2,NR,NTH,NAV,NRAY)=RKY
                     RKB(3,NR,NTH,NAV,NRAY)=RKZ
                     RBB(1,NR,NTH,NAV,NRAY)=RXB
                     RBB(2,NR,NTH,NAV,NRAY)=RYB
                     RBB(3,NR,NTH,NAV,NRAY)=RZB
 2                   CONTINUE
C
C            
C            IF(IFLAG.EQ.1) THEN
C               WRITE(6,'(3I3)') NR,NAV,NCR
C               WRITE(6,'(1P3E12.4)') X,Y,Z
C               WRITE(6,'(1P3E12.4)') RX,RY,RZ
C               WRITE(6,'(1P3E12.4)') DELR2,DELCR2,ARG
C            ENDIF
C
          ENDDO
        ENDDO
       ENDDO
      ENDDO
      ENDIF
C
      DO 2000 NRDO=1,NRMAX
         NR=NRDO
C
         DO 1101 NTH=1,NTHMAX+1
            IF(NTH.NE.NTHMAX/2+1) THEN
               DO 1100 NP=1,NPMAX
                  CALL FPDWAV(ETAG(NTH,NR),SING(NTH),COSG(NTH),PM(NP),
     &                        NR,NTH,DWPPS,DWPTS,DWTPS,DWTTS,NSA)
                  DWTP(NTH,NP,NR,NSA)=DWTPS
                  DWTT(NTH,NP,NR,NSA)=DWTTS
 1100          CONTINUE
            ELSE
               DO 1200 NP=1,NPMAX
C                  CALL FPDWAV(ETAG(NTH,NR),SING(NTH),COSG(NTH),PM(NP),
C     &                        NR,NTH,DWPPS,DWPTS,DWTPS,DWTTS)
                  DWTP(NTH,NP,NR,NSA)=0.D0
                  DWTT(NTH,NP,NR,NSA)=0.D0
 1200          CONTINUE
            ENDIF
 1101    CONTINUE
C
C
C
         IF(MODELA.EQ.1) THEN
            DO 1300 NTH=ITL(NR)+1,NTHMAX/2
            DO 1300 NP=1,NPMAX
               DWTP(NTH,NP,NR,NSA)=(DWTP(NTH,NP,NR,NSA)
     &                            -DWTP(NTHMAX-NTH+2,NP,NR,NSA))*FACT
               DWTT(NTH,NP,NR,NSA)=(DWTT(NTH,NP,NR,NSA)
     &                            +DWTT(NTHMAX-NTH+2,NP,NR,NSA))*FACT
               DWTP(NTHMAX-NTH+2,NP,NR,NSA)=-DWTP(NTH,NP,NR,NSA)
               DWTT(NTHMAX-NTH+2,NP,NR,NSA)= DWTT(NTH,NP,NR,NSA)
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
      SUBROUTINE FPDWAV(ETA,RSIN,RCOS,P,NR,NTH,
     &                  DWPPS,DWPTS,DWTPS,DWTTS,NSA)
C
      INCLUDE 'fpcomm.inc'
      INCLUDE '../wr/wrcom1.inc'
      INCLUDE 'fpcom2.inc'
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
            RFDW=RAYIN(1,NRAY)
C
            IF(MODELW.EQ.1) THEN
             DO NCR=1,NCRMAX(NR,NRAY)
             RX=RCR(1,NCR,NR,NRAY)
             RY=RCR(2,NCR,NR,NRAY)
             RZ=RCR(3,NCR,NR,NRAY)
             RLCR=SQRT(RX**2+RY**2)
             ZLCR=RZ
             DELR2=(RL-RLCR)**2+(ZL-ZLCR)**2
             DELCR2=DELYEC**2
             ARG=DELR2/DELCR2
C
              IF(ARG.LT.15.D0) THEN
               FACTOR=EXP(-ARG)
C               WRITE(25,*) 'NAV,NCR,DELR2,ARG=',NAV,NCR,DELR2,ARG
C               CALL PLMAG(X,Y,Z,RHON)
C               WRITE(6,*) RHOR,RM(NR)
               CALL FPDWRP(NR,ETAL,RSIN,RCOS,PSIN,PCOS,PSI)
               CEX=CECR(1,NCR,NR,NRAY)*FACTOR
               CEY=CECR(2,NCR,NR,NRAY)*FACTOR
               CEZ=CECR(3,NCR,NR,NRAY)*FACTOR
               RKX=RKCR(1,NCR,NR,NRAY)
               RKY=RKCR(2,NCR,NR,NRAY)
               RKZ=RKCR(3,NCR,NR,NRAY)
               CALL FPDWLL(P,PSIN,PCOS,
     &                     CEX,CEY,CEZ,RKX,RKY,RKZ,RX,RY,RZ,
     &                     DWPPL,DWPTL,DWTPL,DWTTL,NSA)
               DWPPS=DWPPS+DWPPL*RCOS/PCOS
               DWPTS=DWPTS+DWPTL          /SQRT(PSI)
               DWTPS=DWTPS+DWTPL          /SQRT(PSI)
               DWTTS=DWTTS+DWTTL*PCOS/RCOS/PSI
              ENDIF
             ENDDO
            ENDIF
C
C
            IF(MODELW.EQ.2) THEN
               ARG=ARGB(NR,NTH,NAV,NRAY)
               IF(ARG.GT.0.D0.AND.ARG.LT.15.D0) THEN
                  FACTOR= EXP(-ARG)
                  CALL FPDWRP(NR,ETAL,RSIN,RCOS,PSIN,PCOS,PSI)
                  CEX=CEB(1,NR,NTH,NAV,NRAY)*FACTOR
                  CEY=CEB(2,NR,NTH,NAV,NRAY)*FACTOR
                  CEZ=CEB(3,NR,NTH,NAV,NRAY)*FACTOR
                  RKX=RKB(1,NR,NTH,NAV,NRAY)
                  RKY=RKB(2,NR,NTH,NAV,NRAY)
                  RKZ=RKB(3,NR,NTH,NAV,NRAY)
                  RX =RBB(1,NR,NTH,NAV,NRAY)
                  RY =RBB(2,NR,NTH,NAV,NRAY)
                  RZ =RBB(3,NR,NTH,NAV,NRAY)
                  CALL FPDWLL(P,PSIN,PCOS,
     &                 CEX,CEY,CEZ,RKX,RKY,RKZ,RX,RY,RZ,
     &                 DWPPL,DWPTL,DWTPL,DWTTL,NSA)
C
                  DWPPS=DWPPS+DWPPL*RCOS/PCOS
                  DWPTS=DWPTS+DWPTL          /SQRT(PSI)
                  DWTPS=DWTPS+DWTPL          /SQRT(PSI)
                  DWTTS=DWTTS+DWTTL*PCOS/RCOS/PSI
C                  WRITE(6,*) NR,NTH,DWPPS
C                  IF(IFLAG.EQ.1) THEN
C                     WRITE(6,'(3I3)') NR,NAV,NCR
C                     WRITE(6,'(1P3E12.4)') X,Y,Z
C                     WRITE(6,'(1P3E12.4)') RX,RY,RZ
C                     WRITE(6,'(1P3E12.4)') DELR2,DELCR2,ARG
C                  ENDIF
               ENDIF
            ENDIF
C
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
      INCLUDE 'fpcomm.inc'
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
     &                  DWPPL,DWPTL,DWTPL,DWTTL,NSA)
C
      INCLUDE 'fpcomm.inc'
C
      PARAMETER(NJMAX=100)
      DIMENSION RJ(0:NJMAX),DRJ(0:NJMAX)
C
      CALL PLMAG(RX,RY,RZ,RHON)
C
      RW     =2.D0*PI*RFDW*1.D6
      RWC    =AEFP(NSA)*BABS/AMFP(NSA)
C
      RKPARA = RKX*BNX +RKY*BNY +RKZ*BNZ
      RKPERP = SQRT(RKX*RKX+RKY*RKY+RKZ*RKZ-RKPARA*RKPARA)
C
         U1X = (RKX-BNX*RKPARA)/RKPERP
         U1Y = (RKY-BNY*RKPARA)/RKPERP
         U1Z = (RKZ-BNZ*RKPARA)/RKPERP
C
         U2X = (BNY*RKZ-BNZ*RKY)/RKPERP
         U2Y = (BNZ*RKX-BNX*RKZ)/RKPERP
         U2Z = (BNX*RKY-BNY*RKX)/RKPERP
C 
      CE1    = CEX*U1X+CEY*U1Y+CEZ*U1Z
      CE2    = CEX*U2X+CEY*U2Y+CEZ*U2Z
      CEPARA = CEX*BNX+CEY*BNY+CEZ*BNZ
C
      CEPLUS =(CE1+CI*CE2)/SQRT(2.D0)
      CEMINUS=(CE1-CI*CE2)/SQRT(2.D0)
C
      RGAMMA =SQRT(1.D0+P*P*THETA0(NSA))
      PPARA  =PTFP0(NSA)*P*PCOS
      PPERP  =PTFP0(NSA)*P*PSIN
      VPARA  =PPARA/(AMFP(NSA)*RGAMMA)
      VPERP  =PPERP/(AMFP(NSA)*RGAMMA)
C 
      DWC11=0.D0
      DWC12=0.D0
      DWC21=0.D0
      DWC22=0.D0
      RKW  =RKPARA/RW
      RGZAI=RKPERP*PPERP/(RWC*AMFP(NSA))
C
C      WRITE(26,*) 'RKPERP,PPERP,RWC,RGZAI = ',RKPERP,PPERP,RWC,RGZAI
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
     &              /(2*NC*RWC*AMFP(NSA))
         ENDIF
         RTHETA2=ABS(CTHETA)**2
         IF (NC.EQ.0) THEN
            A11=0
            A12=0
            A21=0
            A22=RTHETA2*RKW**2/(AMFP(NSA)**2*RGAMMA**2)
	 ELSE
            A11=RTHETA2*(1.D0-RKW*VPARA)**2
            A12=RTHETA2*RKW*VPERP*(1.D0-RKW*VPARA)
            A21=RTHETA2*RKW*VPERP*(1.D0-RKW*VPARA)
            A22=RTHETA2*RKW**2*VPERP**2
	 ENDIF        
         IF(VPARA.EQ.0.D0) THEN
            DWC=0.D0  
         ELSE
            EX=-((RGAMMA-RKPARA*PPARA/(RW*AMFP(NSA))-NC*RWC/RW)
     &           /(PPARA*DELNPR/(AMFP(NSA)*VC)))**2
            IF (EX.LT.-100.D0) THEN 
                DWC=0.D0
            ELSE
                DWC=0.5D0*SQRT(PI)*AEFP(NSA)**2*EXP(EX)/PTFP0(NSA)**2
     &              /(RW*ABS(PPARA)*DELNPR/(AMFP(NSA)*VC))
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
