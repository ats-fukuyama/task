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
      CHARACTER KLINE*80,KWORD*(NCHM),KWD*(NCHM),KID*1,KTAIL*7
      CHARACTER KG1*1,KG2*1
      DIMENSION KWORD(NWDM)
C
      CALL WFVLIM
      IF(XNDMAX.EQ.XNDMIN) THEN
         WRITE(6,*) 'XX NO DATA IS LOADED FOR GRAPHICS'
         GOTO 9000
      ENDIF
C
    1 WRITE(6,*) '# INPUT : E/B/D,X/Y/Z/+/-/P/F,[R/I/A][Xx/Yy/Zz/B6]/',
     &                     'Xyz/Yxz/Zxy/A9 0-9 V,0-9'
      WRITE(6,*) '          P/F,1/2,CXx/CYy/CZz/Xyz/Yxz/Zxy  X=EXIT'
      CALL GUFLSH
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
         IF(KID.NE.' ') THEN
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
      IF(KG1.EQ.'M') THEN
         DO NCH=3,9
            IF(KWD(NCH:NCH).EQ.' ') GOTO 30
         ENDDO
         NCH=NCH-1
   30    KTAIL=KWD(3:NCH-1)//' '
C         WRITE(6,*) 'KTAIL = /',KTAIL,'/'
C
         IF(KG2.EQ.'E') THEN
            KLINE='EXR'//KTAIL//'EXI'//KTAIL//
     &            'EYR'//KTAIL//'EYI'//KTAIL//
     &            'EZR'//KTAIL//'EZI'//KTAIL//
     &            'P1C'//KTAIL//'P2C'//KTAIL
            GOTO 9
         ELSEIF(KG2.EQ.'D') THEN
            KLINE='DXR'//KTAIL//'DXI'//KTAIL//
     &            'DYR'//KTAIL//'DYI'//KTAIL//
     &            'DZR'//KTAIL//'DZI'//KTAIL//
     &            'P1C'//KTAIL//'P2C'//KTAIL
            GOTO 9
         ELSEIF(KG2.EQ.'B') THEN
            KLINE='BXR'//KTAIL//'BXI'//KTAIL//
     &            'BYR'//KTAIL//'BYI'//KTAIL//
     &            'BZR'//KTAIL//'BZI'//KTAIL//
     &            'P1C'//KTAIL//'P2C'//KTAIL
            GOTO 9
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2:',KG2
            GOTO 1
         ENDIF
      ELSEIF(KG1.EQ.'1'.OR.KG1.EQ.'2'.OR.KG1.EQ.'3'.OR.
     &       KG1.EQ.'4'.OR.KG1.EQ.'5'.OR.KG1.EQ.'6'.OR.
     &       KG1.EQ.'7'.OR.KG1.EQ.'8'.OR.KG1.EQ.'9'.OR.
     &       KG1.EQ.'0') THEN
         READ(KG1,'(I1)',ERR=1,END=1) I
         KLINE=KGINX(I)
         GOTO 9
      ELSEIF(KG1.EQ.'V') THEN
         IF(KG2.EQ.'1'.OR.KG2.EQ.'2'.OR.KG2.EQ.'3'.OR.
     &      KG2.EQ.'4'.OR.KG2.EQ.'5'.OR.KG2.EQ.'6'.OR.
     &      KG2.EQ.'7'.OR.KG2.EQ.'8'.OR.KG2.EQ.'9'.OR.
     &      KG2.EQ.'0') THEN
            READ(KG2,'(I1)',ERR=1,END=1) I
            KLINE=KGINV(I)
            GOTO 9
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID2:',KG2
            GOTO 1
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
         READ(KG2,'(I1)',ERR=1,END=1) NID
         IF(    NID.EQ.0) THEN
            NGRAPH=0
         ELSEIF(NID.EQ.1) THEN
            NGRAPH=1
         ELSEIF(NID.EQ.2) THEN
            NGRAPH=2
         ELSEIF(NID.EQ.3) THEN
            NGRAPH=3
         ELSEIF(NID.EQ.4) THEN
            NGRAPH=4
         ELSEIF(NID.EQ.5) THEN
            NGRAPH=5
         ELSEIF(NID.EQ.6) THEN
            NGRAPH=6
         ELSE
            WRITE(6,*) 'XX UNKNOWN NGRAPH:',NID
            GOTO 1
         ENDIF
         DO NW=1,NWMAX-1
            KWORD(NW)=KWORD(NW+1)
         ENDDO
         NWMAX=NWMAX-1
         IF(NWMAX.EQ.0) GOTO 1
         GOTO 21
      ELSEIF(KWD(1:1).EQ.'X') THEN
         GOTO 9000
      ENDIF
C
      IF(NWMAX.EQ.0) GOTO 1
C
      IF(NGRAPH.GE.1) CALL PAGES
      DO 1000 NW=1,NWMAX
         KWD=KWORD(NW)
         KID=KWD(1:1)
         IF(    KID.EQ.'E') THEN
             KID=KWD(2:2)
            IF(    KID.EQ.'X') THEN
               CALL WFCTOG(CEF,3,1,KWD)
            ELSEIF(KID.EQ.'Y') THEN
               CALL WFCTOG(CEF,3,2,KWD)
            ELSEIF(KID.EQ.'Z') THEN
               CALL WFCTOG(CEF,3,3,KWD)
            ELSEIF(KID.EQ.'+') THEN
               CALL WFCTOG(CEP,3,1,KWD)
            ELSEIF(KID.EQ.'-') THEN
               CALL WFCTOG(CEP,3,2,KWD)
            ELSEIF(KID.EQ.'P') THEN
               CALL WFCTOG(CEP,3,3,KWD)
            ELSE
               WRITE(6,*) 'XX UNKNOWN KID2:',KID
               GOTO 1000
            ENDIF
         ELSEIF(    KID.EQ.'D') THEN
             KID=KWD(2:2)
            IF(    KID.EQ.'X') THEN
               CALL WFCTOG(CEP,3,-1,KWD)
            ELSEIF(KID.EQ.'Y') THEN
               CALL WFCTOG(CEP,3,-2,KWD)
            ELSEIF(KID.EQ.'Z') THEN
               CALL WFCTOG(CEP,3,-3,KWD)
            ELSEIF(KID.EQ.'+') THEN
               CALL WFCTOG(CEP,3,1,KWD)
            ELSEIF(KID.EQ.'-') THEN
               CALL WFCTOG(CEP,3,2,KWD)
            ELSEIF(KID.EQ.'P') THEN
               CALL WFCTOG(CEP,3,3,KWD)
            ELSE
               WRITE(6,*) 'XX UNKNOWN KID2:',KID
               GOTO 1000
            ENDIF
         ELSEIF(KID.EQ.'B') THEN
            KID=KWD(2:2)
            IF(    KID.EQ.'X') THEN
               CALL WFCTOG(CBF,3,1,KWD)
            ELSEIF(KID.EQ.'Y') THEN
               CALL WFCTOG(CBF,3,2,KWD)
            ELSEIF(KID.EQ.'Z') THEN
               CALL WFCTOG(CBF,3,3,KWD)
            ELSEIF(KID.EQ.'+') THEN
               CALL WFCTOG(CBP,3,1,KWD)
            ELSEIF(KID.EQ.'-') THEN
               CALL WFCTOG(CBP,3,2,KWD)
            ELSEIF(KID.EQ.'P') THEN
               CALL WFCTOG(CBP,3,3,KWD)
            ELSE
               WRITE(6,*) 'XX UNKNOWN KID2:',KID
               GOTO 1000
            ENDIF
         ELSEIF(KID.EQ.'P') THEN
            KID=KWD(2:2)
            READ(KID,*,ERR=1000) ID
            IF(ID.EQ.0) THEN
               CALL WFDTOG(PABSTN,1,KWD)
            ELSEIF(ID.GT.NSMAX) THEN
               CALL WFDTOG(PABSKN,ID-NSMAX,KWD)
            ELSEIF(ID.GE.1) THEN
               CALL WFDTOG(PABSSN,ID,KWD)
            ELSE
               WRITE(6,*) 'XX UNKNOWN KID2:',KID
               GOTO 1000
            ENDIF
         ELSEIF(KID.EQ.'F') THEN
            KID=KWD(2:2)
            IF(    KID.EQ.'X') THEN
               CALL WFDTOG(PFV,1,KWD)
            ELSEIF(KID.EQ.'Y') THEN
               CALL WFDTOG(PFV,2,KWD)
            ELSEIF(KID.EQ.'Z') THEN
               CALL WFDTOG(PFV,3,KWD)
            ELSE
               WRITE(6,*) 'XX UNKNOWN KID2:',KID
               GOTO 1000
            ENDIF
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID1:',KID
            GOTO 1000
         ENDIF
C
         KID=KWD(3:3)
         IF(    KID.EQ.'R'.OR.
     &          KID.EQ.'I'.OR.
     &          KID.EQ.'A'.OR.
     &          KID.EQ.'C') THEN
            IF(NGRAPH.EQ.0) THEN
               CALL WFGWFC(KWD)
            ELSEIF(NGRAPH.EQ.1) THEN
               CALL WFGPPC(NW,NWMAX,KWD)
            ELSEIF(NGRAPH.EQ.2) THEN
               CALL WFGPFC(NW,NWMAX,KWD)
            ELSEIF(NGRAPH.GE.3) THEN
               CALL WFGPBC(NW,NWMAX,KWD)
            ENDIF
         ELSEIF(KID.EQ.'X'.OR.
     &          KID.EQ.'Y'.OR.
     &          KID.EQ.'Z'.OR.
     &          KID.EQ.'B') THEN
            IF(NGRAPH.EQ.0) THEN
               CALL WFGWFR(KWD)
            ELSE
               CALL WFGPFR(NW,NWMAX,KWD)
            ENDIF
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID3:',KID
         ENDIF
 1000 CONTINUE
      IF(NGRAPH.GE.1) THEN
         CALL WFGPRM
         CALL PAGEE
      ENDIF
      GOTO 1
C
 9000 RETURN
      END
C
C     ****** EXTRACT FROM COMPLEX DATA ******
C
      SUBROUTINE WFCTOG(CZ,IDM,ID,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION CZ(IDM,NNM),CE(IDM)
      CHARACTER KWD*(NCHM),KID*1,KIDC*1
C
      KID=KWD(3:3)
      IE=0
C
      IF(KID.EQ.'R'.OR.
     &   KID.EQ.'I'.OR.
     &   KID.EQ.'A') THEN
         KIDC=KWD(4:4)
C
         IF(KIDC.EQ.'X') THEN
            READ(KWD(5:NCHM),*,ERR=9000) XPOS
            DZ=(ZNDMAX-ZNDMIN)/(NGYMAX-1)
            DY=(YNDMAX-YNDMIN)/(NGXMAX-1)
            DO NGX=1,NGXMAX
               G2X(NGX)=GUCLIP(YNDMIN+DY*(NGX-1))
            ENDDO
            DO NGY=1,NGYMAX
               G2Y(NGY)=GUCLIP(ZNDMIN+DZ*(NGY-1))
            ENDDO
            DO NGY=1,NGYMAX
               Z=ZNDMIN+DZ*(NGY-1)
               DO NGX=1,NGXMAX
                  Y=YNDMIN+DY*(NGX-1)
                  CALL FEP(XPOS,Y,Z,IE)
                  IF(IE.EQ.0) THEN
                     GZ(NGX,NGY)=0.0
                  ELSE
                     IF(ID.LT.0) THEN
                        CALL FIELDE(IE,XPOS,Y,Z,CE)
                        I=-ID
                        ZR=DBLE(CE(I))
                        ZI=IMAG(CE(I))
                     ELSE
                        CALL FIELDC(IE,XPOS,Y,Z,CZ,IDM,ID,ZR,ZI)
                     ENDIF
                     IF(KID.EQ.'R') THEN
                        GZ(NGX,NGY)=GUCLIP(ZR)
                     ELSE IF(KID.EQ.'I') THEN
                        GZ(NGX,NGY)=GUCLIP(ZI)
                     ELSE IF(KID.EQ.'A') THEN
                        GZ(NGX,NGY)=GUCLIP(SQRT(ZR**2+ZI**2))
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(KIDC.EQ.'Y') THEN
            READ(KWD(5:NCHM),*,ERR=9000) YPOS
            DZ=(ZNDMAX-ZNDMIN)/(NGYMAX-1)
            DX=(XNDMAX-XNDMIN)/(NGXMAX-1)
            DO NGX=1,NGXMAX
               G2X(NGX)=GUCLIP(XNDMIN+DX*(NGX-1))
            ENDDO
            DO NGY=1,NGYMAX
               G2Y(NGY)=GUCLIP(ZNDMIN+DZ*(NGY-1))
            ENDDO
            DO NGY=1,NGYMAX
               Z=ZNDMIN+DZ*(NGY-1)
               DO NGX=1,NGXMAX
                  X=XNDMIN+DX*(NGX-1)
                  CALL FEP(X,YPOS,Z,IE)
                  IF(IE.EQ.0) THEN
                     GZ(NGX,NGY)=0.0
                  ELSE
                     IF(ID.LT.0) THEN
                        CALL FIELDE(IE,X,YPOS,Z,CE)
                        I=-ID
                        ZR=DBLE(CE(I))
                        ZI=IMAG(CE(I))
                     ELSE
                        CALL FIELDC(IE,X,YPOS,Z,CZ,IDM,ID,ZR,ZI)
                     ENDIF
                     IF(KID.EQ.'R') THEN
                        GZ(NGX,NGY)=GUCLIP(ZR)
                     ELSE IF(KID.EQ.'I') THEN
                        GZ(NGX,NGY)=GUCLIP(ZI)
                     ELSE IF(KID.EQ.'A') THEN
                        GZ(NGX,NGY)=GUCLIP(SQRT(ZR**2+ZI**2))
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(KIDC.EQ.'Z') THEN
            READ(KWD(5:NCHM),*,ERR=9000) ZPOS
            DY=(YNDMAX-YNDMIN)/(NGYMAX-1)
            DX=(XNDMAX-XNDMIN)/(NGXMAX-1)
            DO NGX=1,NGXMAX
               G2X(NGX)=GUCLIP(XNDMIN+DX*(NGX-1))
            ENDDO
            DO NGY=1,NGYMAX
               G2Y(NGY)=GUCLIP(YNDMIN+DY*(NGY-1))
            ENDDO
            DO NGY=1,NGYMAX
               Y=YNDMIN+DY*(NGY-1)
               DO NGX=1,NGXMAX
                  X=XNDMIN+DX*(NGX-1)
                  CALL FEP(X,Y,ZPOS,IE)
                  IF(IE.EQ.0) THEN
                     GZ(NGX,NGY)=0.0
                  ELSE
                     IF(ID.LT.0) THEN
                        CALL FIELDE(IE,X,Y,ZPOS,CE)
                        I=-ID
                        ZR=DBLE(CE(I))
                        ZI=IMAG(CE(I))
                     ELSE
                        CALL FIELDC(IE,X,Y,ZPOS,CZ,IDM,ID,ZR,ZI)
                     ENDIF
                     IF(KID.EQ.'R') THEN
                        GZ(NGX,NGY)=GUCLIP(ZR)
                     ELSE IF(KID.EQ.'I') THEN
                        GZ(NGX,NGY)=GUCLIP(ZI)
                     ELSE IF(KID.EQ.'A') THEN
                        GZ(NGX,NGY)=GUCLIP(SQRT(ZR**2+ZI**2))
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(KIDC.EQ.'A') THEN
            READ(KWD(5:NCHM),*,ERR=9000) X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3
            DX2=(X2-X1)/(NGXMAX-1)
            DY2=(Y2-Y1)/(NGXMAX-1)
            DZ2=(Z2-Z1)/(NGXMAX-1)
            DX3=(X3-X1)/(NGYMAX-1)
            DY3=(Y3-Y1)/(NGYMAX-1)
            DZ3=(Z3-Z1)/(NGYMAX-1)
C
            IF(ABS(DX2).GE.ABS(DY2)) THEN
               IF(ABS(DX2).GE.ABS(DZ2)) THEN
                  DO NGX=1,NGXMAX
                     G2X(NGX)=GUCLIP(X1+DX2*(NGX-1))
                  ENDDO
               ELSE
                  DO NGX=1,NGXMAX
                     G2X(NGX)=GUCLIP(Z1+DZ2*(NGX-1))
                  ENDDO
               ENDIF
            ELSE
               IF(ABS(DY2).GE.ABS(DZ2)) THEN
                  DO NGX=1,NGXMAX
                     G2X(NGX)=GUCLIP(Y1+DY2*(NGX-1))
                  ENDDO
               ELSE
                  DO NGX=1,NGXMAX
                     G2X(NGX)=GUCLIP(Z1+DZ2*(NGX-1))
                  ENDDO
               ENDIF
            ENDIF
            IF(ABS(DX3).GE.ABS(DY3)) THEN
               IF(ABS(DX3).GE.ABS(DZ3)) THEN
                  DO NGY=1,NGYMAX
                     G2Y(NGY)=GUCLIP(X1+DX3*(NGY-1))
                  ENDDO
               ELSE
                  DO NGY=1,NGYMAX
                     G2Y(NGY)=GUCLIP(Z1+DZ3*(NGY-1))
                  ENDDO
               ENDIF
            ELSE
               IF(ABS(DY3).GE.ABS(DZ3)) THEN
                  DO NGY=1,NGYMAX
                     G2Y(NGY)=GUCLIP(Y1+DY2*(NGY-1))
                  ENDDO
               ELSE
                  DO NGY=1,NGYMAX
                     G2Y(NGY)=GUCLIP(Z1+DZ2*(NGY-1))
                  ENDDO
               ENDIF
            ENDIF
C
            DO NGY=1,NGYMAX
               X0=X1+DX3*(NGY-1)
               Y0=Y1+DY3*(NGY-1)
               Z0=Z1+DZ3*(NGY-1)
               DO NGX=1,NGXMAX
                  X=X0+DX2*(NGX-1)
                  Y=Y0+DY2*(NGX-1)
                  Z=Z0+DZ2*(NGX-1)
                  CALL FEP(X,Y,Z,IE)
                  IF(IE.EQ.0) THEN
C                     WRITE(6,'(2I5,1P,3E12.4)') NGX,NGY,X,Y,Z
                     GZ(NGX,NGY)=0.0
                  ELSE
                     IF(ID.LT.0) THEN
                        CALL FIELDE(IE,X,Y,Z,CE)
                        I=-ID
                        ZR=DBLE(CE(I))
                        ZI=IMAG(CE(I))
                     ELSE
                        CALL FIELDC(IE,X,Y,Z,CZ,IDM,ID,ZR,ZI)
C                        WRITE(6,'(2I5,1P,5E12.4)') NGX,NGY,X,Y,Z,ZR,ZI
                     ENDIF
                     IF(KID.EQ.'R') THEN
                        GZ(NGX,NGY)=GUCLIP(ZR)
                     ELSE IF(KID.EQ.'I') THEN
                        GZ(NGX,NGY)=GUCLIP(ZI)
                     ELSE IF(KID.EQ.'A') THEN
                        GZ(NGX,NGY)=GUCLIP(SQRT(ZR**2+ZI**2))
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ELSEIF(KID.EQ.'X') THEN
         READ(KWD(4:NCHM),*,ERR=9000) YPOS,ZPOS
         DX=(XNDMAX-XNDMIN)/(NGVMAX-1)
         DO NGV=1,NGVMAX
            X=XNDMIN+DX*(NGV-1)
            CALL FEP(X,YPOS,ZPOS,IE)
            IF(IE.EQ.0) THEN
               GX(NGV)=GUCLIP(X)
               GV(NGV,1)=0.0
               GV(NGV,2)=0.0
               GV(NGV,3)=0.0
            ELSE
               IF(ID.LT.0) THEN
                  CALL FIELDE(IE,X,YPOS,ZPOS,CE)
                  I=-ID
                  ZR=DBLE(CE(I))
                  ZI=IMAG(CE(I))
               ELSE
                  CALL FIELDC(IE,X,YPOS,ZPOS,CZ,IDM,ID,ZR,ZI)
               ENDIF
               GX(NGV)=GUCLIP(X)
               GV(NGV,1)=GUCLIP(ZR)
               GV(NGV,2)=GUCLIP(ZI)
               GV(NGV,3)=GUCLIP(SQRT(ZR**2+ZI**2))
            ENDIF
         ENDDO
      ELSEIF(KID.EQ.'Y') THEN
         READ(KWD(4:NCHM),*,ERR=9000) XPOS,ZPOS
         DY=(YNDMAX-YNDMIN)/(NGVMAX-1)
         DO NGV=1,NGVMAX
            Y=YNDMIN+DY*(NGV-1)
            CALL FEP(XPOS,Y,ZPOS,IE)
            IF(IE.EQ.0) THEN
               GX(NGV)=GUCLIP(Y)
               GV(NGV,1)=0.0
               GV(NGV,2)=0.0
               GV(NGV,3)=0.0
            ELSE
               IF(ID.LT.0) THEN
                  CALL FIELDE(IE,XPOS,Y,ZPOS,CE)
                  I=-ID
                  ZR=DBLE(CE(I))
                  ZI=IMAG(CE(I))
               ELSE
                  CALL FIELDC(IE,XPOS,Y,ZPOS,CZ,IDM,ID,ZR,ZI)
               ENDIF
               GX(NGV)=GUCLIP(Y)
               GV(NGV,1)=GUCLIP(ZR)
               GV(NGV,2)=GUCLIP(ZI)
               GV(NGV,3)=GUCLIP(SQRT(ZR**2+ZI**2))
            ENDIF
         ENDDO
      ELSEIF(KID.EQ.'Z') THEN
         READ(KWD(4:NCHM),*,ERR=9000) XPOS,YPOS
         DZ=(ZNDMAX-ZNDMIN)/(NGVMAX-1)
         DO NGV=1,NGVMAX
            Z=ZNDMIN+DZ*(NGV-1)
            CALL FEP(XPOS,YPOS,Z,IE)
            IF(IE.EQ.0) THEN
               GX(NGV)=GUCLIP(Z)
               GV(NGV,1)=0.0
               GV(NGV,2)=0.0
               GV(NGV,3)=0.0
            ELSE
               IF(ID.LT.0) THEN
                  CALL FIELDE(IE,XPOS,YPOS,Z,CE)
                  I=-ID
                  ZR=DBLE(CE(I))
                  ZI=IMAG(CE(I))
               ELSE
                  CALL FIELDC(IE,XPOS,YPOS,Z,CZ,IDM,ID,ZR,ZI)
               ENDIF
               GX(NGV)=GUCLIP(Z)
               GV(NGV,1)=GUCLIP(ZR)
               GV(NGV,2)=GUCLIP(ZI)
               GV(NGV,3)=GUCLIP(SQRT(ZR**2+ZI**2))
            ENDIF
         ENDDO
      ELSEIF(KID.EQ.'B') THEN
         READ(KWD(4:NCHM),*,ERR=9000) X1,Y1,Z1,X2,Y2,Z2
         DX=(X2-X1)/(NGVMAX-1)
         DY=(Y2-Y1)/(NGVMAX-1)
         DZ=(Z2-Z1)/(NGVMAX-1)
C
         IF(ABS(DX).GE.ABS(DY)) THEN
            IF(ABS(DX).GE.ABS(DZ)) THEN
               W1=X1
               DW=DX
            ELSE
               W1=Z1
               DW=DZ
            ENDIF
         ELSE
            IF(ABS(DY).GE.ABS(DZ)) THEN
               W1=Y1
               DW=DY
            ELSE
               W1=Z1
               DW=DZ
            ENDIF
         ENDIF
C
         DO NGV=1,NGVMAX
            X=X1+DX*(NGV-1)
            Y=Y1+DY*(NGV-1)
            Z=Z1+DZ*(NGV-1)
            GX(NGV)=GUCLIP(W1+DW*(NGV-1))
            CALL FEP(X,Y,Z,IE)
            IF(IE.EQ.0) THEN
               GV(NGV,1)=0.0
               GV(NGV,2)=0.0
               GV(NGV,3)=0.0
            ELSE
               IF(ID.LT.0) THEN
                  CALL FIELDE(IE,X,Y,Z,CE)
                  I=-ID
                  ZR=DBLE(CE(I))
                  ZI=IMAG(CE(I))
               ELSE
                  CALL FIELDC(IE,X,Y,Z,CZ,IDM,ID,ZR,ZI)
               ENDIF
               GV(NGV,1)=GUCLIP(ZR)
               GV(NGV,2)=GUCLIP(ZI)
               GV(NGV,3)=GUCLIP(SQRT(ZR**2+ZI**2))
            ENDIF
         ENDDO
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID3:',KID
      ENDIF
 9000 RETURN
      END
C
C     ****** EXTRACT FROM DOUBLE DATA ******
C
      SUBROUTINE WFDTOG(DV,ID,KWD)
C
      INCLUDE 'wfcomm.inc'
C
      DIMENSION DV(NNM,ID)
      CHARACTER KWD*(NCHM),KID*1,KIDC*1
C
      KID=KWD(3:3)
      IE=0
C
      IF(KID.EQ.'C') THEN
         KIDC=KWD(4:4)
C
         IF(KIDC.EQ.'X') THEN
            READ(KWD(5:NCHM),*,ERR=9000) XPOS
            DZ=(ZNDMAX-ZNDMIN)/(NGYMAX-1)
            DY=(YNDMAX-YNDMIN)/(NGXMAX-1)
            DO NGX=1,NGXMAX
               G2X(NGX)=GUCLIP(YNDMIN+DY*(NGX-1))
            ENDDO
            DO NGY=1,NGYMAX
               G2Y(NGY)=GUCLIP(ZNDMIN+DZ*(NGY-1))
            ENDDO
            DO NGY=1,NGYMAX
               Z=ZNDMIN+DZ*(NGY-1)
               DO NGX=1,NGXMAX
                  Y=YNDMIN+DY*(NGX-1)
                  CALL FEP(XPOS,Y,Z,IE)
                  IF(IE.EQ.0) THEN
                     GZ(NGX,NGY)=0.0
C                     WRITE(6,'(A,1P3E12.4,I5)') 
C     &                    'XX DTOG: X,Y,Z,IE=',XPOS,Y,Z,IE
                  ELSE
                     CALL FIELDD(IE,XPOS,Y,Z,DV,ID,V)
                     GZ(NGX,NGY)=GUCLIP(V)
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(KIDC.EQ.'Y') THEN
            READ(KWD(5:NCHM),*,ERR=9000) YPOS
            DZ=(ZNDMAX-ZNDMIN)/(NGYMAX-1)
            DX=(XNDMAX-XNDMIN)/(NGXMAX-1)
            DO NGX=1,NGXMAX
               G2X(NGX)=GUCLIP(XNDMIN+DX*(NGX-1))
            ENDDO
            DO NGY=1,NGYMAX
               G2Y(NGY)=GUCLIP(ZNDMIN+DZ*(NGY-1))
            ENDDO
            DO NGY=1,NGYMAX
               Z=ZNDMIN+DZ*(NGY-1)
               DO NGX=1,NGXMAX
                  X=XNDMIN+DX*(NGX-1)
                  CALL FEP(X,YPOS,Z,IE)
                  IF(IE.EQ.0) THEN
C                     WRITE(6,'(A,1P3E12.4,I5)') 
C     &                    'XX DTOG: X,Y,Z,IE=',X,YPOS,Z,IE
                     GZ(NGX,NGY)=0.0
                  ELSE
                     CALL FIELDD(IE,X,YPOS,Z,DV,ID,V)
                     GZ(NGX,NGY)=GUCLIP(V)
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(KIDC.EQ.'Z') THEN
            READ(KWD(5:NCHM),*,ERR=9000) ZPOS
            DY=(YNDMAX-YNDMIN)/(NGYMAX-1)
            DX=(XNDMAX-XNDMIN)/(NGXMAX-1)
            DO NGX=1,NGXMAX
               G2X(NGX)=GUCLIP(XNDMIN+DX*(NGX-1))
            ENDDO
            DO NGY=1,NGYMAX
               G2Y(NGY)=GUCLIP(YNDMIN+DY*(NGY-1))
            ENDDO
            DO NGY=1,NGYMAX
               Y=YNDMIN+DY*(NGY-1)
               DO NGX=1,NGXMAX
                  X=XNDMIN+DX*(NGX-1)
                  CALL FEP(X,Y,ZPOS,IE)
                  IF(IE.EQ.0) THEN
                     GZ(NGX,NGY)=0.0
C                     WRITE(6,'(A,1P3E12.4,I5)') 
C     &                    'XX DTOG: X,Y,Z,IE=',X,Y,ZPOS,IE
                  ELSE
                     CALL FIELDD(IE,X,Y,ZPOS,DV,ID,V)
                     GZ(NGX,NGY)=GUCLIP(V)
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(KIDC.EQ.'A') THEN
            READ(KWD(5:NCHM),*,ERR=9000) X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3
            DX2=(X2-X1)/(NGXMAX-1)
            DY2=(Y2-Y1)/(NGXMAX-1)
            DZ2=(Z2-Z1)/(NGXMAX-1)
            DX3=(X3-X1)/(NGYMAX-1)
            DY3=(Y3-Y1)/(NGYMAX-1)
            DZ3=(Z3-Z1)/(NGYMAX-1)
C
            IF(ABS(DX2).GE.ABS(DY2)) THEN
               IF(ABS(DX2).GE.ABS(DZ2)) THEN
                  DO NGX=1,NGXMAX
                     G2X(NGX)=GUCLIP(X1+DX2*(NGX-1))
                  ENDDO
               ELSE
                  DO NGX=1,NGXMAX
                     G2X(NGX)=GUCLIP(Z1+DZ2*(NGX-1))
                  ENDDO
               ENDIF
            ELSE
               IF(ABS(DY2).GE.ABS(DZ2)) THEN
                  DO NGX=1,NGXMAX
                     G2X(NGX)=GUCLIP(Y1+DY2*(NGX-1))
                  ENDDO
               ELSE
                  DO NGX=1,NGXMAX
                     G2X(NGX)=GUCLIP(Z1+DZ2*(NGX-1))
                  ENDDO
               ENDIF
            ENDIF
            IF(ABS(DX3).GE.ABS(DY3)) THEN
               IF(ABS(DX3).GE.ABS(DZ3)) THEN
                  DO NGY=1,NGYMAX
                     G2Y(NGY)=GUCLIP(X1+DX3*(NGY-1))
                  ENDDO
               ELSE
                  DO NGY=1,NGYMAX
                     G2Y(NGY)=GUCLIP(Z1+DZ3*(NGY-1))
                  ENDDO
               ENDIF
            ELSE
               IF(ABS(DY3).GE.ABS(DZ3)) THEN
                  DO NGY=1,NGYMAX
                     G2Y(NGY)=GUCLIP(Y1+DY2*(NGY-1))
                  ENDDO
               ELSE
                  DO NGY=1,NGYMAX
                     G2Y(NGY)=GUCLIP(Z1+DZ2*(NGY-1))
                  ENDDO
               ENDIF
            ENDIF
C
            DO NGY=1,NGYMAX
               X0=X1+DX3*(NGY-1)
               Y0=Y1+DY3*(NGY-1)
               Z0=Z1+DZ3*(NGY-1)
               DO NGX=1,NGXMAX
                  X=X0+DX2*(NGX-1)
                  Y=Y0+DY2*(NGX-1)
                  Z=Z0+DZ2*(NGX-1)
                  CALL FEP(X,Y,Z,IE)
                  IF(IE.EQ.0) THEN
                     GZ(NGX,NGY)=0.0
                  ELSE
                     CALL FIELDD(IE,X,Y,Z,DV,ID,V)
                     GZ(NGX,NGY)=GUCLIP(V)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ELSEIF(KID.EQ.'X') THEN
         READ(KWD(4:NCHM),*,ERR=9000) YPOS,ZPOS
         DX=(XNDMAX-XNDMIN)/(NGVMAX-1)
         DO NGV=1,NGVMAX
            X=XNDMIN+DX*(NGV-1)
            CALL FEP(X,YPOS,ZPOS,IE)
            IF(IE.EQ.0) THEN
               GX(NGV)=GUCLIP(X)
               GV(NGV,1)=0.0
               GV(NGV,2)=0.0
               GV(NGV,3)=0.0
            ELSE
               CALL FIELDD(IE,X,YPOS,ZPOS,DV,ID,V)
               GX(NGV)=GUCLIP(X)
               GV(NGV,1)=GUCLIP(Z)
            ENDIF
         ENDDO
      ELSEIF(KID.EQ.'Y') THEN
         READ(KWD(4:NCHM),*,ERR=9000) XPOS,ZPOS
         DY=(YNDMAX-YNDMIN)/(NGVMAX-1)
         DO NGV=1,NGVMAX
            Y=YNDMIN+DY*(NGV-1)
            CALL FEP(XPOS,Y,ZPOS,IE)
            IF(IE.EQ.0) THEN
               GX(NGV)=GUCLIP(Y)
               GV(NGV,1)=0.0
               GV(NGV,2)=0.0
               GV(NGV,3)=0.0
            ELSE
               CALL FIELDD(IE,XPOS,Y,ZPOS,DV,ID,V)
               GX(NGV)=GUCLIP(Y)
               GV(NGV,1)=GUCLIP(V)
            ENDIF
         ENDDO
      ELSEIF(KID.EQ.'Z') THEN
         READ(KWD(4:NCHM),*,ERR=9000) XPOS,YPOS
         DZ=(ZNDMAX-ZNDMIN)/(NGVMAX-1)
         DO NGV=1,NGVMAX
            Z=ZNDMIN+DZ*(NGV-1)
            CALL FEP(XPOS,YPOS,Z,IE)
            IF(IE.EQ.0) THEN
               GX(NGV)=GUCLIP(Z)
               GV(NGV,1)=0.0
               GV(NGV,2)=0.0
               GV(NGV,3)=0.0
            ELSE
               CALL FIELDD(IE,XPOS,YPOS,Z,DV,ID,V)
               GX(NGV)=GUCLIP(Z)
               GV(NGV,1)=GUCLIP(V)
            ENDIF
         ENDDO
      ELSEIF(KID.EQ.'B') THEN
         READ(KWD(4:NCHM),*,ERR=9000) X1,Y1,Z1,X2,Y2,Z2
         DX=(X2-X1)/(NGVMAX-1)
         DY=(Y2-Y1)/(NGVMAX-1)
         DZ=(Z2-Z1)/(NGVMAX-1)
C
         IF(ABS(DX).GE.ABS(DY)) THEN
            IF(ABS(DX).GE.ABS(DZ)) THEN
               W1=X1
               DW=DX
            ELSE
               W1=Z1
               DW=DZ
            ENDIF
         ELSE
            IF(ABS(DY).GE.ABS(DZ)) THEN
               W1=Y1
               DW=DY
            ELSE
               W1=Z1
               DW=DZ
            ENDIF
         ENDIF
C
         DO NGV=1,NGVMAX
            X=X1+DX*(NGV-1)
            Y=Y1+DY*(NGV-1)
            Z=Z1+DZ*(NGV-1)
            GX(NGV)=GUCLIP(W1+DW*(NGV-1))
            CALL FEP(X,Y,Z,IE)
            IF(IE.EQ.0) THEN
               GV(NGV,1)=0.0
            ELSE
               CALL FIELDD(IE,X,Y,Z,DV,ID,V)
               GV(NGV,1)=GUCLIP(V)
            ENDIF
         ENDDO
      ELSE
         WRITE(6,*) 'XX UNKNOWN KID3:',KID
      ENDIF
 9000 RETURN
      END
C
C     ****** WRITE 2D PROFILE IN TEXT FILE ******
C
      SUBROUTINE WFGWFC(KWD)
C
      INCLUDE 'wfcomm.inc'
C
      CHARACTER KWD*(NCHM)
C
      NFD=22
      WRITE(NFD,'(A79)') KWD
      WRITE(NFD,'(I8)') 2
      WRITE(NFD,'(2I8)') NGXMAX,NGYMAX
      WRITE(NFD,'(1P5E15.7)') (G2X(NGX),NGX=1,NGXMAX)
      WRITE(NFD,'(1P5E15.7)') (G2Y(NGY),NGY=1,NGYMAX)
      WRITE(NFD,'(1P5E15.7)') ((GZ(NGX,NGY),NGX=1,NGXMAX),NGY=1,NGYMAX)
C
      RETURN
      END
C
C     ****** WRITE 1D PROFILE IN TEXT FILE ******
C
      SUBROUTINE WFGWFR(KWD)
C
      INCLUDE 'wfcomm.inc'
C
      CHARACTER KWD*(NCHM)
C
      IF(KWD(1:1).EQ.'E'.OR.
     &   KWD(1:1).EQ.'D'.OR.
     &   KWD(1:1).EQ.'B'.OR.
     &   KWD(1:1).EQ.'A') THEN
         NGMAX=3
      ELSE
         NGMAX=1
      ENDIF
C
      NFD=22
      WRITE(NFD,'(A79)') KWD
      WRITE(NFD,'(I8)') 1
      WRITE(NFD,'(2I8)') NGVMAX,NGMAX
      WRITE(NFD,'(1P5E15.7)') (GX(NGV),NGV=1,NGVMAX)
      WRITE(NFD,'(1P5E15.7)') ((GV(NGV,NG),NGV=1,NGVMAX),NG=1,NGMAX)
C
      RETURN
      END
