C     $Id$
C     *********  /TASK/WF/GOUT  *********
C
C       GRAPHIC DATA PROCESSING PROGRAM
C             FOR FEM COMPUTATION
C
C     ***********************************
C
      SUBROUTINE WFGOUT
C
      USE libchar
      INCLUDE 'wfcomm.inc'
C
      CHARACTER KLINE*80,KWORD*(NCHM),KWD*(NCHM),KWTEMP*(NCHM)
      CHARACTER KID*1,KG1*1,KG2*1,KG3*1,KG4*1
      DIMENSION KWORD(NWDM)
      DIMENSION GZ(NNODM)
      DIMENSION GZ1(NNODM),GZ2(NNODM),GZ3(NNODM),GZ4(NNODM)
      DIMENSION GT(NGTM,3)
      DIMENSION TEMP(NNODM,2),AL(3)
C
      IF(XDMAX.EQ.XDMIN) THEN
         WRITE(6,*) 'XX NO DATA IS LOADED FOR GRAPHICS'
         GOTO 9000
      ENDIF
C
      CALL WFGINI
      CALL WFCALB
      CALL WFCALD
C
    1 WRITE(6,*) '# INPUT : E/B/A,X/Y/Z/+/-/P/,R/I/A/L/P/Xpos/Ypos'
      WRITE(6,*) '          P,P/B/N/T/G/M/S/F/C/R/U/V/D/E/H/K/I,',
     &                                            'E/I/1/2,C/Xpos/Ypos'
      WRITE(6,*) '          T,N/T/F/P,1/2  0-9  V,0-9  L,n  G,1-4  ',
     &                                            'X=EXIT'
      READ(5,'(A80)',ERR=1,END=9000) KLINE
      NWXMAX=0
C
    9 NL=0
      NWD=0
      NCH=0
   10 IF(NL.GE.80) GOTO 20
         NL=NL+1
         KID=KLINE(NL:NL)
         CALL toupper(KID)
         IF(KID.NE.','.AND.KID.NE.' ') THEN
            IF(NCH.LT.NCHM) NCH=NCH+1
            KWD(NCH:NCH)=KID
         ELSE
            IF(NCH.NE.0) THEN
               IF(NWD.LT.NWDM) NWD=NWD+1
               KWORD(NWD)=KWD(1:NCH)
               NCH=0
            ENDIF
         ENDIF
      GOTO 10
C
   20 NWMAX=NWD
   21 KWD=KWORD(1)
      KG1=KWD(1:1)
      KG2=KWD(2:2)
      IF(KG1.EQ.'0') THEN
         KLINE=KGINX(0)
         GOTO 9
      ELSEIF(KG1.EQ.'1') THEN
         KLINE=KGINX(1)
         GOTO 9
      ELSEIF(KG1.EQ.'2') THEN
         KLINE=KGINX(2)
         GOTO 9
      ELSEIF(KG1.EQ.'3') THEN
         KLINE=KGINX(3)
         GOTO 9
      ELSEIF(KG1.EQ.'4') THEN
         KLINE=KGINX(4)
         GOTO 9
      ELSEIF(KG1.EQ.'5') THEN
         KLINE=KGINX(5)
         GOTO 9
      ELSEIF(KG1.EQ.'6') THEN
         KLINE=KGINX(6)
         GOTO 9
      ELSEIF(KG1.EQ.'7') THEN
         KLINE=KGINX(7)
         GOTO 9
      ELSEIF(KG1.EQ.'8') THEN
         KLINE=KGINX(8)
         GOTO 9
      ELSEIF(KG1.EQ.'9') THEN
         KLINE=KGINX(9)
         GOTO 9
      ELSEIF(KG1.EQ.'V') THEN
         IF(    KG2.EQ.'0') THEN
            KLINE=KGINV(0)
            GOTO 9
         ELSEIF(KG2.EQ.'1') THEN
            KLINE=KGINV(1)
            GOTO 9
         ELSEIF(KG2.EQ.'2') THEN
            KLINE=KGINV(2)
            GOTO 9
         ELSEIF(KG2.EQ.'3') THEN
            KLINE=KGINV(3)
            GOTO 9
         ELSEIF(KG2.EQ.'4') THEN
            KLINE=KGINV(4)
            GOTO 9
         ELSEIF(KG2.EQ.'5') THEN
            KLINE=KGINV(5)
            GOTO 9
         ELSEIF(KG2.EQ.'6') THEN
            KLINE=KGINV(6)
            GOTO 9
         ELSEIF(KG2.EQ.'7') THEN
            KLINE=KGINV(7)
            GOTO 9
         ELSEIF(KG2.EQ.'8') THEN
            KLINE=KGINV(8)
            GOTO 9
         ELSEIF(KG2.EQ.'9') THEN
            KLINE=KGINV(9)
            GOTO 9
         ENDIF
      ELSEIF(KG1.EQ.'L') THEN
         READ(KG2,'(I1)',ERR=1,END=1) NWXMAX
         DO NW=1,NWMAX-1
            KWORD(NW)=KWORD(NW+1)
         ENDDO
         NWMAX=NWMAX-1
         IF(NWMAX.EQ.0) GOTO 1
         GOTO 21
      ELSEIF(KG1.EQ.'G') THEN
         READ(KG2,'(I1)',ERR=1,END=1) NDIM
         IF(    NDIM.EQ.1) THEN
            NGRAPH=1
         ELSEIF(NDIM.EQ.2) THEN
            NGRAPH=2
         ELSEIF(NDIM.EQ.3) THEN
            NGRAPH=3
         ELSEIF(NDIM.EQ.4) THEN
            NGRAPH=4
         ENDIF
         DO NW=1,NWMAX-1
            KWORD(NW)=KWORD(NW+1)
         ENDDO
         NWMAX=NWMAX-1
         IF(NWMAX.EQ.0) GOTO 1
         GOTO 21
      ELSEIF(KG1.EQ.'X') THEN
         GOTO 9000
      ENDIF
C
      IF(KG1.NE.'E'.AND.
     &   KG1.NE.'B'.AND.
     &   KG1.NE.'A'.AND.
     &   KG1.NE.'P'.AND.
     &   KG1.NE.'D'.AND.
     &   KG1.NE.'T') GOTO 1
C
      CALL PAGES
      DO 1000 NW=1,NWMAX
         KWD=KWORD(NW)
         KG1=KWD(1:1)
         KG2=KWD(2:2)
         KG3=KWD(3:3)
         KG4=KWD(4:4)
C
C        +++++ WAVE ELECTRIC FIELD +++++
C
         IF(    KG1.EQ.'E') THEN
            IDP=0
            IF(    KG2.EQ.'X') THEN
               CALL WFCTOG(CEF,GZ,1,KWD)
            ELSEIF(KG2.EQ.'Y') THEN
               CALL WFCTOG(CEF,GZ,2,KWD)
            ELSEIF(KG2.EQ.'Z') THEN
               CALL WFCTOG(CEF,GZ,3,KWD)
            ELSEIF(KG2.EQ.'+') THEN
               CALL WFCTOG(CEP,GZ,1,KWD)
            ELSEIF(KG2.EQ.'-') THEN
               CALL WFCTOG(CEP,GZ,2,KWD)
            ELSEIF(KG2.EQ.'P') THEN
               CALL WFCTOG(CEP,GZ,3,KWD)
            ELSE
               WRITE(6,*) 'XX UNKNOWN KG2:',KG2
               GOTO 1000
            ENDIF
C
C        +++++ WAVE MAGNETIC FIELD +++++
C
         ELSEIF(KG1.EQ.'B') THEN
            IDP=0
            IF(    KG2.EQ.'X') THEN
               CALL WFCTOG(CBF,GZ,1,KWD)
            ELSEIF(KG2.EQ.'Y') THEN
               CALL WFCTOG(CBF,GZ,2,KWD)
            ELSEIF(KG2.EQ.'Z') THEN
               CALL WFCTOG(CBF,GZ,3,KWD)
            ELSEIF(KG2.EQ.'+') THEN
               CALL WFCTOG(CBP,GZ,1,KWD)
            ELSEIF(KG2.EQ.'-') THEN
               CALL WFCTOG(CBP,GZ,2,KWD)
            ELSEIF(KG2.EQ.'P') THEN
               CALL WFCTOG(CBP,GZ,3,KWD)
            ELSE
               WRITE(6,*) 'XX UNKNOWN KG2:',KG2
               GOTO 1000
            ENDIF
C
C        +++++ WAVE POTENTIAL FIELD +++++
C
         ELSEIF(KG1.EQ.'A') THEN
            IDP=0
            IF(    KG2.EQ.'X') THEN
               CALL WFATOG(CAF,GZ,1,KWD)
            ELSEIF(KG2.EQ.'Y') THEN
               CALL WFATOG(CAF,GZ,2,KWD)
            ELSEIF(KG2.EQ.'Z') THEN
               CALL WFATOG(CAF,GZ,3,KWD)
            ELSEIF(KG2.EQ.'F') THEN
               CALL WFATOG(CAF,GZ,4,KWD)
            ELSE
               WRITE(6,*) 'XX UNKNOWN KG2:',KG2
               GOTO 1000
            ENDIF
C
C        +++++ VARIOUS PROFILES +++++
C
         ELSEIF(KG1.EQ.'P') THEN
            IDP=1
            IF(KG3.EQ.'E'.OR.KG3.EQ.'1') THEN
               IG3=1
            ELSEIF(KG3.EQ.'I'.OR.KG3.EQ.'2') THEN
               IG3=2
            ELSE
               IG3=1
            ENDIF              
C
C        +++++ ABSORBED POWER +++++
C
            IF(KG2.EQ.'P') THEN
               IDP=0
               CALL WFDTOG(PWRNS,GZ,IG3,KWD)
               KWTEMP='PP1X0.16'
               CALL WFDTOG(PWRNS,GZ1,IG3,KWTEMP)
               KWTEMP='PP1X0.15'
               CALL WFDTOG(PWRNS,GZ1,IG3,KWTEMP)
               KWTEMP='PP1X0.14'
               CALL WFDTOG(PWRNS,GZ1,IG3,KWTEMP)
               KWTEMP='PP1X0.13'
               CALL WFDTOG(PWRNS,GZ1,IG3,KWTEMP)
C
C        +++++ STATIC MAGNETIC FIELD (BABS, PSI) +++++
C
            ELSEIF(KG2.EQ.'B') THEN
               DO IN=1,NNOD
                  CALL WFBMAG(XD(IN),YD(IN),TEMP(IN,1),AL)
                  CALL WFBPSI(XD(IN),YD(IN),TEMP(IN,2))
               ENDDO
               CALL WFDTOG(TEMP,GZ,IG3,KWD)
C
C        +++++ PLASMA PRESSURE (ELECTRON, ION) +++++
C
            ELSEIF(KG2.EQ.'M') THEN
               IF(IG3.EQ.1) THEN
                  CALL WFDTOG(PPE,GZ,1,KWD)
               ELSEIF(IG3.EQ.2) THEN
                  CALL WFDTOG(PPI,GZ,1,KWD)
               ENDIF
C
C        +++++ PLASMA DENSITY (ELECTRON, ION) +++++
C
            ELSEIF(KG2.EQ.'N') THEN
               IF(IG3.EQ.1) THEN
                  CALL WFDTOG(PNE,GZ,1,KWD)
               ELSEIF(IG3.EQ.2) THEN
                  CALL WFDTOG(PNI,GZ,1,KWD)
               ENDIF
C
C        +++++ PLASMA TEMPERATURE (ELECTRON, ION) +++++
C
            ELSEIF(KG2.EQ.'T') THEN
               IF(IG3.EQ.1) THEN
                  CALL WFDTOG(PTE,GZ,1,KWD)
               ELSEIF(IG3.EQ.2) THEN
                  CALL WFDTOG(PTI,GZ,1,KWD)
               ENDIF
C
C        +++++ DRIFT VELOCITY (ELECTRON, ION) +++++
C
            ELSEIF(KG2.EQ.'G') THEN
               IF(IG3.EQ.1) THEN
                  CALL WFDTOG(VDE,GZ,1,KWD)
               ELSEIF(IG3.EQ.2) THEN
                  CALL WFDTOG(VDI,GZ,1,KWD)
               ENDIF
C
C        +++++ PLASMA POTENTIAL +++++
C
            ELSEIF(KG2.EQ.'F') THEN
               CALL WFDTOG(PHI,GZ,1,KWD)
               IDP=0
C
C        +++++ SURFACE CHARGE +++++
C
            ELSEIF(KG2.EQ.'S') THEN
               CALL WFDTOG(PSQ,GZ,1,KWD)
C
C        +++++ SPACE CHARGE +++++
C
            ELSEIF(KG2.EQ.'C') THEN
               DO IN=1,NNOD
                  TEMP(IN,1)=PZ(1)*PNE(IN)+PZ(2)*PNI(IN)
               ENDDO
               CALL WFDTOG(TEMP,GZ,1,KWD)
C
C        +++++ COLLISION FREQUENCY (ELECTRON, ION) +++++
C
            ELSEIF(KG2.EQ.'R') THEN
               CALL WFDTOG(RNUG,GZ,IG3,KWD)
C
C        +++++ PARALLEL MOBILITY (ELECTRON, ION) +++++
C
            ELSEIF(KG2.EQ.'U') THEN
               CALL WFDTOG(UPARAG,GZ,IG3,KWD)
C
C        +++++ PARPENDICULAR MOBILITY (ELECTRON, ION) +++++
C
            ELSEIF(KG2.EQ.'V') THEN
               CALL WFDTOG(UPERPG,GZ,IG3,KWD)
C
C        +++++ PARALLEL PARTICLE DIFFUSIVITY (ELECTRON, ION) +++++
C
            ELSEIF(KG2.EQ.'D') THEN
               CALL WFDTOG(DPARAG,GZ,IG3,KWD)
C
C        +++++ PARPENDICULAR PARTICLE  DIFFUSIVITY (ELECTRON, ION) +++++
C
            ELSEIF(KG2.EQ.'E') THEN
               CALL WFDTOG(DPERPG,GZ,IG3,KWD)
C
C        +++++ PARALLEL HEAT DIFFUSIVITY (ELECTRON, ION) +++++
C
            ELSEIF(KG2.EQ.'H') THEN
               CALL WFDTOG(HPARAG,GZ,IG3,KWD)
C
C        +++++ PERPENDICULAR HEAT DIFFUSIVITY (ELECTRON, ION) +++++
C
            ELSEIF(KG2.EQ.'K') THEN
               CALL WFDTOG(HPERPG,GZ,IG3,KWD)
C
C        +++++ IONIZATION FREQUECY +++++
C
            ELSEIF(KG2.EQ.'I') THEN
               CALL WFDTOG(RIONG,GZ,1,KWD)
            ELSE
               WRITE(6,*) 'XX UNKNOWN KG2:',KG2
               GOTO 1000
            ENDIF
C
C        +++++ TIME DEPENDENCE +++++
C
         ELSEIF(KG1.EQ.'T') THEN
            IF(KG3.EQ.'E'.OR.KG3.EQ.'1') THEN
               IG3=1
            ELSEIF(KG3.EQ.'I'.OR.KG3.EQ.'2') THEN
               IG3=2
            ELSE
               IG3=1
            ENDIF              
C
            IF(KG2.EQ.'N') THEN
               CALL WFTTOG(PNT,GT,IG3,KWD)
            ELSEIF(KG2.EQ.'T') THEN
               CALL WFTTOG(PTT,GT,IG3,KWD)
            ELSEIF(KG2.EQ.'F') THEN
               CALL WFTTOG(PHIT,GT,1,KWD)
            ELSEIF(KG2.EQ.'P') THEN
               CALL WFTTOG(PABST,GT,1,KWD)
            ELSE
               WRITE(6,*) 'XX UNKNOWN KG2:',KG2
               GOTO 1000
            ENDIF
C
C        +++++ DISPERSION +++++
C
         ELSEIF(KG1.EQ.'D') THEN
            IDP=3
            IF(KG3.EQ.'1') THEN
               IG3=1
            ELSEIF(KG3.EQ.'2') THEN
               IG3=2
            ELSEIF(KG3.EQ.'3') THEN
               IG3=3
            ELSEIF(KG3.EQ.'4') THEN
               IG3=4
            ELSEIF(KG3.EQ.'5') THEN
               IG3=5
            ELSEIF(KG3.EQ.'6') THEN
               IDP=2
               IG3=6
            ELSE
               WRITE(6,*) 'XX UNKNOWN KG3:',KG3
               GOTO 1000
            ENDIF              
C
            IF(KG2.EQ.'P') THEN
               CALL WFDTOG(DISPG,GZ,IG3,KWD)
            ELSEIF(KG2.EQ.'Q') THEN
               IDP=0
               CALL WFDTOG(DISPH,GZ,IG3,KWD)
            ELSE
               WRITE(6,*) 'XX UNKNOWN KG2:',KG2
               GOTO 1000
            ENDIF              
         ELSE
            WRITE(6,*) 'XX UNKNOWN KG1:',KG1
            GOTO 1000
         ENDIF
C
         IF(KG1.EQ.'T') THEN
            CALL WFGPFT(NW,NWMAX,GT,KWD)
         ELSE IF(KG1.EQ.'P') THEN
            IF(KG4.EQ.'C') THEN
               CALL WFGPFX(NW,NWMAX,GZ,KWD,IDP)
            ELSEIF(KG4.EQ.'X'.OR.
     &             KG4.EQ.'Y') THEN
               IF(NGRAPH.EQ.1) THEN
                  CALL WFGPPR(NW,NWMAX,GZ,KWD)
               ELSE
                  CALL WFGPPRX(NW,NWMAX,GZ1,GZ2,GZ3,GZ4,KWD)
               ENDIF
            ELSE
               WRITE(6,*) 'XX UNKNOWN KG4:',KG4
            ENDIF
         ELSE IF(KG1.EQ.'D') THEN
            IF(KG4.EQ.'C') THEN
               CALL WFGPFX(NW,NWMAX,GZ,KWD,IDP)
            ELSEIF(KG4.EQ.'X'.OR.
     &             KG4.EQ.'Y') THEN
               CALL WFGPPR(NW,NWMAX,GZ,KWD)
            ELSE
               WRITE(6,*) 'XX UNKNOWN KG4:',KG4
            ENDIF
         ELSE
            IF(    KG3.EQ.'R'.OR.
     &             KG3.EQ.'I'.OR.
     &             KG3.EQ.'A'.OR.
     &             KG3.EQ.'L'.OR.
     &             KG3.EQ.'P') THEN
                CALL WFGPFX(NW,NWMAX,GZ,KWD,IDP)
            ELSEIF(KG3.EQ.'X'.OR.
     &             KG3.EQ.'Y') THEN
               IF(NGRAPH.EQ.1) THEN
                  CALL WFGPPR(NW,NWMAX,GZ,KWD)
               ELSEIF(NGRAPH.EQ.2) THEN
                  CALL WFGPFR(NW,NWMAX,GZ,KWD)
               ELSEIF(NGRAPH.EQ.3) THEN
                  CALL WFGPAR(NW,NWMAX,GZ,KWD)
               ELSE
                  CALL WFGPFR(NW,NWMAX,GZ,KWD)
               ENDIF
            ELSE
               WRITE(6,*) 'XX UNKNOWN KG3:',KG3
            ENDIF
         ENDIF
 1000 CONTINUE
      CALL WFGPRM
      CALL PAGEE
      GOTO 1
C
 9000 RETURN
      END
C
C     ****** EXTRACT FROM COMPLEX DATA ******
C
      SUBROUTINE WFCTOG(CZ,GZ,ID,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CZ(3,NNODM),GZ(NNODM)
      DIMENSION X(NRM),Y(NRM)
      CHARACTER KWD*(NCHM),KID*1
C
      KID=KWD(3:3)
C
      IF(    KID.EQ.'R') THEN
         DO IN=1,NNOD
            GZ(IN)=GCLIP(DBLE(CZ(ID,IN)))
         ENDDO
      ELSEIF(KID.EQ.'I') THEN
         DO IN=1,NNOD
            GZ(IN)=GCLIP(DIMAG(CZ(ID,IN)))
         ENDDO
      ELSEIF(KID.EQ.'A') THEN
         DO IN=1,NNOD
            GZ(IN)=GCLIP(ABS(CZ(ID,IN)))
         ENDDO
      ELSEIF(KID.EQ.'L') THEN
         DO IN=1,NNOD
            GZ(IN)=GCLIP(LOG10(MAX(ABS(CZ(ID,IN)),1.D-5)))
         ENDDO
      ELSEIF(KID.EQ.'P') THEN
         DO IN=1,NNOD
            GZ(IN)=GCLIP(ATAN2(DBLE(CZ(ID,IN)),DIMAG(CZ(ID,IN))))
         ENDDO
      ELSEIF(KID.EQ.'X') THEN
         READ(KWD(4:NCHM),*,ERR=9000) YPOS
         DX=(XDMAX-XDMIN)/(NRMAX-1)
         DO N=1,NRMAX
            X(N)=XDMIN+DX*(N-1)
            CALL WFFEP(X(N),YPOS,IE)
            IF(IE.EQ.0) THEN
               GZ(        N)=GCLIP(X(N))
               GZ(  NRMAX+N)=0.0
               GZ(2*NRMAX+N)=0.0
               GZ(3*NRMAX+N)=0.0
            ELSE
               CALL FIELDC(IE,X(N),YPOS,CZ,3,ID,ZR,ZI)
               GZ(        N)=GCLIP(X(N))
               GZ(  NRMAX+N)=GCLIP(ZR)
               GZ(2*NRMAX+N)=GCLIP(ZI)
               GZ(3*NRMAX+N)=GCLIP(SQRT(ZR**2+ZI**2))
            ENDIF
         ENDDO
      ELSEIF(KID.EQ.'Y') THEN
         READ(KWD(4:NCHM),*,ERR=9000) XPOS
         DY=(YDMAX-YDMIN)/(NRMAX-1)
         DO N=1,NRMAX
            Y(N)=YDMIN+DY*(N-1)
            CALL WFFEP(XPOS,Y(N),IE)
            IF(IE.EQ.0) THEN
               GZ(        N)=GCLIP(Y(N))
               GZ(  NRMAX+N)=0.0
               GZ(2*NRMAX+N)=0.0
               GZ(3*NRMAX+N)=0.0
            ELSE
               CALL FIELDC(IE,XPOS,Y(N),CZ,3,ID,ZR,ZI)
               GZ(        N)=GCLIP(Y(N))
               GZ(  NRMAX+N)=GCLIP(ZR)
               GZ(2*NRMAX+N)=GCLIP(ZI)
               GZ(3*NRMAX+N)=GCLIP(SQRT(ZR**2+ZI**2))
            ENDIF
         ENDDO
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID2:',KID
      ENDIF
 9000 RETURN
      END
C
C     ****** EXTRACT FROM COMPLEX DATA ******
C
      SUBROUTINE WFATOG(CZ,GZ,ID,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CZ(4,NNODM),GZ(NNODM)
      DIMENSION X(NRM),Y(NRM)
      CHARACTER KWD*(NCHM),KID*1
C
      KID=KWD(3:3)
C
      IF(    KID.EQ.'R') THEN
         DO IN=1,NNOD
            GZ(IN)=GCLIP(DBLE(CZ(ID,IN)))
         ENDDO
      ELSEIF(KID.EQ.'I') THEN
         DO IN=1,NNOD
            GZ(IN)=GCLIP(DIMAG(CZ(ID,IN)))
         ENDDO
      ELSEIF(KID.EQ.'A') THEN
         DO IN=1,NNOD
            GZ(IN)=GCLIP(ABS(CZ(ID,IN)))
         ENDDO
      ELSEIF(KID.EQ.'X') THEN
         READ(KWD(4:NCHM),*,ERR=9000) YPOS
         DX=(XDMAX-XDMIN)/(NRMAX-1)
         DO N=1,NRMAX
            X(N)=XDMIN+DX*(N-1)
            CALL WFFEP(X(N),YPOS,IE)
            IF(IE.EQ.0) THEN
               GZ(        N)=GCLIP(X(N))
               GZ(  NRMAX+N)=0.0
               GZ(2*NRMAX+N)=0.0
               GZ(3*NRMAX+N)=0.0
            ELSE
               CALL FIELDC(IE,X(N),YPOS,CZ,4,ID,ZR,ZI)
               GZ(        N)=GCLIP(X(N))
               GZ(  NRMAX+N)=GCLIP(ZR)
               GZ(2*NRMAX+N)=GCLIP(ZI)
               GZ(3*NRMAX+N)=GCLIP(SQRT(ZR**2+ZI**2))
            ENDIF
         ENDDO
      ELSEIF(KID.EQ.'Y') THEN
         READ(KWD(4:NCHM),*,ERR=9000) XPOS
         DY=(YDMAX-YDMIN)/(NRMAX-1)
         DO N=1,NRMAX
            Y(N)=YDMIN+DY*(N-1)
            CALL WFFEP(XPOS,Y(N),IE)
            IF(IE.EQ.0) THEN
               GZ(        N)=GCLIP(Y(N))
               GZ(  NRMAX+N)=0.0
               GZ(2*NRMAX+N)=0.0
               GZ(3*NRMAX+N)=0.0
            ELSE
               CALL FIELDC(IE,XPOS,Y(N),CZ,4,ID,ZR,ZI)
               GZ(        N)=GCLIP(Y(N))
               GZ(  NRMAX+N)=GCLIP(ZR)
               GZ(2*NRMAX+N)=GCLIP(ZI)
               GZ(3*NRMAX+N)=GCLIP(SQRT(ZR**2+ZI**2))
            ENDIF
         ENDDO
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID2:',KID
      ENDIF
 9000 RETURN
      END
C
C     ****** EXTRACT FROM DBLE DATA ******
C
      SUBROUTINE WFDTOG(DZ,GZ,ID,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION DZ(NNODM,ID),GZ(NNODM)
      DIMENSION X(NRM),Y(NRM)
      CHARACTER KWD*(NCHM),KID*1
C
      IF(KWD(1:1).EQ.'P'.OR.KWD(1:1).EQ.'D') THEN
         IK=4
      ELSE
         IK=3
      ENDIF
      KID=KWD(IK:IK)
C
      IF(    KID.EQ.'C') THEN
         DO IN=1,NNOD
            GZ(IN)=GCLIP(DZ(IN,ID))
         ENDDO
C         DO IN=1,NNOD,5
C            WRITE(6,'(1P5E12.4)') (DZ(IN+I,ID),I=0,4)
C         ENDDO
      ELSEIF(KID.EQ.'X') THEN
         READ(KWD(IK+1:NCHM),*,ERR=9000) YPOS
         DX=(XDMAX-XDMIN)/(NRMAX-1)
         DO N=1,NRMAX
            X(N)=XDMIN+DX*(N-1)
            CALL WFFEP(X(N),YPOS,IE)
            IF(IE.EQ.0) THEN
               GZ(      N)=GCLIP(X(N))
               GZ(NRMAX+N)=0.0
            ELSE
               CALL FIELDD(IE,X(N),YPOS,DZ,ID,Z)
               GZ(      N)=GCLIP(X(N))
               GZ(NRMAX+N)=GCLIP(Z)
            ENDIF
         ENDDO
      ELSEIF(KID.EQ.'Y') THEN
         READ(KWD(IK+1:NCHM),*,ERR=9000) XPOS
         DY=(YDMAX-YDMIN)/(NRMAX-1)
         DO N=1,NRMAX
            Y(N)=YDMIN+DY*(N-1)
            CALL WFFEP(XPOS,Y(N),IE)
            IF(IE.EQ.0) THEN
               GZ(      N)=GCLIP(Y(N))
               GZ(NRMAX+N)=0.0
            ELSE
               CALL FIELDD(IE,XPOS,Y(N),DZ,ID,Z)
               GZ(      N)=GCLIP(Y(N))
               GZ(NRMAX+N)=GCLIP(Z)
            ENDIF
         ENDDO
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID3:',KID
      ENDIF
 9000 RETURN
      END
C
C     ****** EXTRACT FROM DBLE DATA ******
C
      SUBROUTINE WFTTOG(DGT,GT,ID,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION DGT(NGTM,3,ID),GT(NGTM,3)
      CHARACTER KWD*(NCHM),KID*1
C
      KID=KWD(3:3)
C
      IF(KID.EQ.'L') THEN
         DO I=1,3
            DO NGT=1,NGTMAX
               GT(NGT,I)=GCLIP(LOG10(ABS(DGT(NGT,I,ID))))
            ENDDO
         ENDDO
      ELSE
         DO I=1,3
            DO NGT=1,NGTMAX
               GT(NGT,I)=GCLIP(DGT(NGT,I,ID))
            ENDDO
         ENDDO
      ENDIF
      RETURN
      END
C
C     ****** CALCULATE MAGNETIC FIELD AND POLARIZED FIELD ******
C
      SUBROUTINE WFCALB
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION SS(3)
      REAL*8 C(3)
      DIMENSION XE(3),YE(3),A(3),B(3),AL(3)
C      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
C
      COEFB=-CI*VC/(2.D0*PI*RF*1.D6)
C
      DO 10 NN=1,NNOD
         CBF(1,NN)=0.D0
         CBF(2,NN)=0.D0
         CBF(3,NN)=0.D0
         SA(NN)   =0.D0
   10 CONTINUE
      DO 100 IDO=1,NELM
         IE=IDO
         CALL WFNPOS(IE,XE,YE)
         XG=(XE(1)+XE(2)+XE(3))/3.D0
         YG=(YE(1)+YE(2)+YE(3))/3.D0
         R=RR+XG
         CALL WFABC(IE,A,B,C,S)
         DO 20 J=1,3
            SS(J)=A(J)+B(J)*XG+C(J)*YG
   20    CONTINUE
         CEXG =(0.D0,0.D0)
         CEYG =(0.D0,0.D0)
         CEZG =(0.D0,0.D0)
         CDEXY=(0.D0,0.D0)
         CDEYX=(0.D0,0.D0)
         CDEZX=(0.D0,0.D0)
         CDEZY=(0.D0,0.D0)
         DO 30 K=1,3
            NO=IELM(K,IE)
            CEXG =CEXG +SS(K)*CEF(1,NO)
            CEYG =CEYG +SS(K)*CEF(2,NO)
            CEZG =CEZG +SS(K)*CEF(3,NO)
            CDEXY=CDEXY+ C(K)*CEF(1,NO)
            CDEYX=CDEYX+ B(K)*CEF(2,NO)
            CDEZX=CDEZX+ B(K)*CEF(3,NO)
            CDEZY=CDEZY+ C(K)*CEF(3,NO)
   30    CONTINUE
C
         IF(MODELS.EQ.1.OR.MODELS.EQ.2) THEN
            CBX =COEFB*(        CDEZY  -CI*NPHI*CEYG/R)
            CBY =COEFB*(CI*NPHI*CEXG/R -        CDEZX)
            CBZ =COEFB*(        CDEYX  -        CDEXY)
         ELSE
            CBX =COEFB*(       CDEZY-CI*RKZ*CEYG)
            CBY =COEFB*(CI*RKZ*CEXG -      CDEZX)
            CBZ =COEFB*(       CDEYX-      CDEXY)
         ENDIF
C
         S3=S/3.D0
         DO 40 K=1,3
            NO=IELM(K,IE)
            CBF(1,NO)=CBF(1,NO)+S3*CBX
            CBF(2,NO)=CBF(2,NO)+S3*CBY
            CBF(3,NO)=CBF(3,NO)+S3*CBZ
            SA(NO)=SA(NO)+S3
   40    CONTINUE
  100 CONTINUE
C
      DO 50 NN=1,NNOD
         IF(MOD(KNODW(NN),4).EQ.1) THEN
            IF(NPHI.EQ.0) THEN
               CBF(1,NN)=(0.D0,0.D0)
               CBF(2,NN)=CBF(2,NN)/SA(NN)
               CBF(3,NN)=(0.D0,0.D0)
            ELSEIF(NPHI.EQ.1.OR.NPHI.EQ.-1) THEN
               CBF(1,NN)=CBF(1,NN)/SA(NN)
               CBF(2,NN)=(0.D0,0.D0)
               CBF(3,NN)=CBF(3,NN)/SA(NN)
            ELSE
               CBF(1,NN)=(0.D0,0.D0)
               CBF(2,NN)=(0.D0,0.D0)
               CBF(3,NN)=(0.D0,0.D0)
            ENDIF
         ELSE
            CBF(1,NN)=CBF(1,NN)/SA(NN)
            CBF(2,NN)=CBF(2,NN)/SA(NN)
            CBF(3,NN)=CBF(3,NN)/SA(NN)
         ENDIF
   50 CONTINUE
C
      DO 110 N=1,NNOD
         CALL WFBMAG(XD(N),YD(N),BABS,AL)
         SUM=SQRT(AL(1)*AL(1)+AL(3)*AL(3))
         IF(SUM.EQ.0.D0) THEN
            CE1=CEF(3,N)
            CE2=CEF(1,N)
            CB1=CBF(3,N)
            CB2=CBF(1,N)
         ELSE
            CE1=(AL(3)*CEF(1,N)-AL(1)*CEF(3,N))/SUM
            CE2=SUM*CEF(2,N)-AL(2)*(AL(1)*CEF(1,N)+AL(3)*CEF(3,N))/SUM
            CB1=(AL(3)*CBF(1,N)-AL(1)*CBF(3,N))/SUM
            CB2=SUM*CBF(2,N)-AL(2)*(AL(1)*CBF(1,N)+AL(3)*CBF(3,N))/SUM
         ENDIF
         CEP(1,N)=CE1-CI*CE2
         CEP(2,N)=CE1+CI*CE2
         CEP(3,N)=AL(1)*CEF(1,N)+AL(2)*CEF(2,N)+AL(3)*CEF(3,N)
         CBP(1,N)=CB1-CI*CB2
         CBP(2,N)=CB1+CI*CB2
         CBP(3,N)=AL(1)*CBF(1,N)+AL(2)*CBF(2,N)+AL(3)*CBF(3,N)
  110 CONTINUE
C
      RETURN
      END
C
C     ****** CALCULATE DISPERSION ******
C
      SUBROUTINE WFCALD
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION WP(NSM),WC(NSM)
      DIMENSION RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
      DIMENSION AL(3)
C
      WW=2.D0*PI*RF*1.D6
C
      DO IS=1,NSMAX
         WP(IS)=PZ(IS)*PZ(IS)*AEE*AEE*1.D20/(PA(IS)*AMP*EPS0*WW*WW)
         WC(IS)=PZ(IS)*AEE/(PA(IS)*AMP*WW)
      ENDDO
C
      DO IN=1,NNOD
         CALL WFBMAG(XD(IN),YD(IN),BABS,AL)
         CALL WFSDEN(IN,RN,RTPR,RTPP,RZCL)
C
         CDT=1.D0
         CDP=1.D0
         CDX=0.D0
         DO IS=1,NSMAX
            CWP=WP(IS)*RN(IS)/DCMPLX(1.D0,RZCL(IS))
            CWC=WC(IS)*BABS/(DCMPLX(1.D0,RZCL(IS)))
            CDT=CDT+   CWP    /(1.D0-CWC**2)
            CDX=CDX+CI*CWP*CWC/(1.D0-CWC**2)
            CDP=CDP+   CWP
         ENDDO
C
         WPE= WP(1)*RN(1)
         WPI= WP(2)*RN(2)
         WCE= WC(1)*BABS
         WCI= WC(2)*BABS
C
         DISPG(IN,1)= 1.D0-WPE-WPI
         DISPG(IN,2)= (1.D0-WCE*WCE)*(1.D0-WCI*WCI)
     &               -(1.D0-WCI*WCI)*WPE
     &               -(1.D0-WCE*WCE)*WPI
         DISPG(IN,3)= (1.D0+WCE)*(1.D0+WCI)
     &               -(1.D0+WCI)*WPE
     &               -(1.D0+WCE)*WPI
         DISPG(IN,4)= (1.D0-WCE)*(1.D0-WCI)
     &               -(1.D0-WCI)*WPE
     &               -(1.D0-WCE)*WPI
         DISPG(IN,5)= DISPG(IN,1)*DISPG(IN,2)*DISPG(IN,3)*DISPG(IN,4)
     &               *(1.D0-WCE*WCE)*(1.D0-WCI*WCI)
         DISPG(IN,6)=-(1.D0-WPE-WPI)
     &               /(1.D0-WPE/(1.D0-WCE*WCE)-WPI/(1.D0-WCI*WCI))
C
         DISPH(IN,1)= 1.D0-WPE-WPI
         DISPH(IN,2)= 1.D0-WPE/(1.D0-WCE*WCE)-WPI/(1.D0-WCI*WCI)
         DISPH(IN,3)= 1.D0-WPE/(1.D0+WCE)-WPI/(1.D0+WCI)
         DISPH(IN,4)= 1.D0-WPE/(1.D0-WCE)-WPE/(1.D0-WCI)
         DISPH(IN,5)= DISPG(IN,1)*DISPG(IN,2)*DISPG(IN,3)*DISPG(IN,4)
         DISPH(IN,6)= WPE/(WCE*WCE)
C
      ENDDO
      RETURN
      END
