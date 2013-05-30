C     $Id$
C
C     ****** Define 2-D Node Array (POLYGON) ******
C
      SUBROUTINE DFNODP
C
      INCLUDE 'wfcomm.inc'
      PARAMETER (NFQM=11,NVQM=11,NPQM=101)
      COMMON /WFDVP1/ XPQ(NPQM),YPQ(NPQM),IPQ(NPQM)
      COMMON /WFDVP2/ NPQ(NVQM,NFQM),NVQMAX(NFQM)
      COMMON /WFDVP3/ NVPQ(NVQM,NFQM),NFPQ(NVQM,NFQM)
      COMMON /WFDVP4/ DELF1(NFQM),DELF2(NFQM)
      COMMON /WFDVP5/ DELFQ(NVQM,NFQM)
C
C     NFQM    : Maximum number of polygon
C     NVQM    : Maximum number of vertex of a polygon
C     NPQM    : Maximum number of point
C
C     XPQ,YPQ : Position of point NPQ
C     NPQ     : Point number at vertex NPFQ of polygon NFQ
C     NVMAX   : Number of vertex of poly NFQ
C     NFPQ    : Polygon number sharing the edge NVQ-(NVQ+1) of poly NFQ
C     NVPQ    : Vertex number sharing the edge NVQ-(NVQ+1) of poly NFQ
C     DELF1   : Mesh size parallel to the line 1-2 of poly NFQ
C     DELF2   : Mesh size perpendicular to the line 1-2 of poly NFQ
C     DELFQ   : Mesh size along the line NVQ-(NVQ+1) of poly NFQ
C
      CHARACTER KID*1,KLINE*80
      DATA INITPL/0/
C
    1 IF(INITPL.EQ.0) THEN
         DO NP=1,NPQM
            XPQ(NP)=0.D0
            YPQ(NP)=0.D0
            IPQ(NP)=0
         ENDDO
         DO NF=1,NFQM
            DO NV=1,NVQM
               NPQ(NV,NF)=0
               NVPQ(NV,NF)=0
               NFPQ(NV,NF)=0
            ENDDO
            NVQMAX(NF)=0
            DELF1(NF)=0.01D0
            DELF2(NF)=0.01D0
         ENDDO
         INITPL=1
      ENDIF
      IMODFY=0
C     
      WRITE(6,*) '## DIVP: P np xp yp :define point'
      WRITE(6,*) '##       F nf np1 np2 np3 ...: define polygon'
      WRITE(6,*) '##       D nf df1 df2: mesh size'
      WRITE(6,*) '##       V: view  C: clear  S:save  X:exit'
    2 READ(5,'(A80)',ERR=1,END=9000) KLINE
      KID=KLINE(1:1)
      CALL GUCPTL(KID)
C
      IF(KID.EQ.'P') THEN
         KLINE(1:1)=' '
         READ(KLINE,*,ERR=8,END=8) NP,XP,YP
         IF(NP.GE.1.AND.NP.LE.NPQM) THEN
            XPQ(NP)=XP
            YPQ(NP)=YP
            IPQ(NP)=1
            IMODFY=1
            WRITE(6,*) '## OK'
         ELSE
            WRITE(6,*) 'XX np is out of range'
         ENDIF
         GOTO 2
C
      ELSEIF(KID.EQ.'F') THEN
         NF=0
         NV=-1
         NC=1
   10    NC=NC+1
         IF(NC.GT.80) GOTO 15
         KID=KLINE(NC:NC)
         IF(KID.EQ.' '.OR.KID.EQ.',') GOTO 10
         NC1=NC
   12    NC=NC+1
         KID=KLINE(NC:NC)
         IF(KID.NE.' '.AND.KID.NE.',') GOTO 12
         READ(KLINE(NC1:NC-1),*,ERR=8,END=8) NP
         IF(NV.EQ.-1) THEN
            IF(NP.GE.1.AND.NP.LE.NFQM)THEN
               NV=0
               NF=NP
               GOTO 10
            ELSE
               WRITE(6,*) 'XX nf is out of range'
            ENDIF
         ELSEIF(NV+1.GE.1.AND.NV+1.LE.NVQM) THEN
            IF(NF.GE.1.AND.NF.LE.NFQM)THEN
               IF(NP.GE.1.AND.NP.LE.NPQM)THEN
                  NV=NV+1
                  NPQ(NV,NF)=NP
                  GOTO 10
               ELSE
                  WRITE(6,*) 'XX np is out of range'
               ENDIF
            ELSE
               WRITE(6,*) 'XX nf is out of range'
            ENDIF
         ELSE
            WRITE(6,*) 'XX nv is out of range'
         ENDIF
         GOTO 2
C
   15    IF(NF.GE.1.AND.NF.LE.NFQM) THEN
            IF(NV.GE.3.AND.NV.LE.NVQM) THEN
               NVQMAX(NF)=NV
               IMODFY=1
               WRITE(6,*) '## OK'
            ELSE
               WRITE(6,*) 'XX nv is out of range'
            ENDIF
         ELSE
            WRITE(6,*) 'XX nf is out of range'
         ENDIF
         GOTO 2
C
      ELSEIF(KID.EQ.'D') THEN
         KLINE(1:1)=' '
         READ(KLINE,*,ERR=8,END=8) NF,D1,D2
         IF(NF.GE.1.AND.NF.LE.NFQM) THEN
            IF(NVQMAX(NF).GT.0) THEN
               DELF1(NF)=D1
               DELF2(NF)=D2
               IMODFY=1
               WRITE(6,*) '## OK'
            ELSE
               WRITE(6,*) 'XX polygon(nf) is not defined'
            ENDIF
         ELSE
            WRITE(6,*) 'XX nf is out of range'
         ENDIF
         GOTO 2
C
      ELSEIF(KID.EQ.'V') THEN
         DO NP=1,NPQM
            IF(IPQ(NP).EQ.1) WRITE(6,601) NP,XPQ(NP),YPQ(NP)
  601       FORMAT(' NP:',I4,2F12.4)
         ENDDO
         DO NF=1,NFQM
            IF(NVQMAX(NF).GT.0) THEN
               WRITE(6,602) NF,(NPQ(NV,NF),NV=1,NVQMAX(NF))
  602          FORMAT(' NF:',I4,(I4))
               WRITE(6,603) DELF1(NF),DELF2(NF)
  603          FORMAT('    ',10X,2F12.4)
            ENDIF
         ENDDO
         GOTO 2
C
      ELSEIF(KID.EQ.'C') THEN
         WRITE(6,*) '## OK'
         INITPL=1
C
      ELSEIF(KID.EQ.'S') THEN
         CALL WFANAP(IERR)
         IF(IERR.EQ.0) THEN
            WRITE(6,*) '## OK'
            IMODFY=0
         ELSE
            WRITE(6,*) 'XX WFANAP: error detected.'
         ENDIF
C
      ELSEIF(KID.EQ.'X') THEN
         IF(IMODFY.NE.0) THEN
            WRITE(6,*) '## Data has been modified.'
            WRITE(6,*) '   Are you sure to exit? [Y/N]' 
            READ(5,'(A1)',ERR=8,end=8) KID
            CALL GUCPTL(KID)
            IF(KID.EQ.'Y') GOTO 9000
         ELSE
            GOTO 9000
         ENDIF
         GOTO 2
      ELSE
         WRITE(6,*) 'XX Unknown command'
      ENDIF
      GOTO 1
C
    8 WRITE(6,*) 'XX Input error detected'
      GOTO 1
C
 9000 RETURN
      END
C
C     ****** Define 2-D Node Array (POLYGON) ******
C
      SUBROUTINE WFANAP(IERR)
C
      INCLUDE 'wfcomm.inc'
      PARAMETER (NFQM=11,NVQM=11,NPQM=101)
      COMMON /WFDVP1/ XPQ(NPQM),YPQ(NPQM),IPQ(NPQM)
      COMMON /WFDVP2/ NPQ(NVQM,NFQM),NVQMAX(NFQM)
      COMMON /WFDVP3/ NVPQ(NVQM,NFQM),NFPQ(NVQM,NFQM)
      COMMON /WFDVP4/ DELF1(NFQM),DELF2(NFQM)
      COMMON /WFDVP5/ DELFQ(NVQM,NFQM)
      DATA EPSL/1.D-12/
C
C     +++++ CHECK TWO POINTS WITH SAME POSITION +++++
C
      IERR=0
      DO NP=1,NPQM
         IF(IPQ(NP).EQ.1) THEN
            DO NP1=NP+1,NPQM
               IF(IPQ(NP1).EQ.1) THEN
                  IF((XPQ(NP).EQ.XPQ(NP1)).AND.
     &               (YPQ(NP).EQ.YPQ(NP1))) THEN
                        WRITE(6,*) 'XX Points ',NP,' and ',NP1,
     &                             ' are located same position.'
                        WRITE(6,*) '   X,Y=',XPQ(NP),YPQ(NP)
                        IERR=1
                        GOTO 9000
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO
C
C     +++++ LOOK FOR SHARED LINE +++++
C
      DO NF=1,NFQM
      DO NV=1,NVQMAX(NF)
         NP=NPQ(NV,NF)
         NPP=NPQ(MOD(NV,NVQMAX(NF))+1,NF)
         DO NF1=NF+1,NFQM
         DO NV1=1,NVQMAX(NF1)
            NP1=NPQ(NV1,NF1)
            IF(NP1.EQ.NPP)THEN
               NV2=MOD(NV1,NVQMAX(NF1))+1
               NP2=NPQ(NV2,NF1)
               IF(NP.EQ.NP2)THEN
                  NFPQ(NV,NF)=NF1
                  NVPQ(NV,NF)=NV1
                  NFPQ(NV1,NF1)=NF
                  NVPQ(NV1,NF1)=NV
                  GOTO 100
               ELSE
                  WRITE(6,*) 'XX These edges are inconsistent.'
                  WRITE(6,*) '   NP, NV, NF =',NP, NV, NF
                  WRITE(6,*) '   NP1,NV1,NF1=',NP1,NV1,NF1
                  ERR=2
                  GOTO 9000
               ENDIF
            ENDIF
         ENDDO
         ENDDO
         WRITE(6,*) '** This vertex is isolated: NV,NF=',NV,NF
  100    CONTINUE
      ENDDO
      ENDDO
C
C     +++++ DEFINE DELFQ FROM DELF1,DELF2 +++++
C
      DO NF=1,NFQM
         DO NV=1,NVQMAX(NF)
            DELFQ(NV,NF)=0.D0
         ENDDO
      ENDDO
C
      DO NF=1,NFQM
      IF(NVQMAX(NF).GE.3) THEN
         DELFL1=DELF1(NF)
         DELFL2=DELF2(NF)
C
C     --- CALCULATE DIRECTIONAL COSINE OF EDGE1 
C
         NV=1
         NP1=NPQ(NV,NF)
         NP2=NPQ(NV+1,NF)
         DX=XPQ(NP2)-XPQ(NP1)
         DY=YPQ(NP2)-YPQ(NP1)
         DL=SQRT(DX*DX+DY*DY)
         DX=DX/DL
         DY=DY/DL
C
C     --- SET DELFQ
C
         NL=(DL-EPSL)/DELFL1+1
         DELL=DL/NL
         DEL=DELFQ(NV,NF)
         IF(DEL.EQ.0.D0) THEN
            DEL=DELL
         ELSE
            IF(DELL.LT.DEL) THEN
               DEL=DELL
               DELFQ(NVPQ(NV,NF),NFPQ(NV,NF))=DELL
            ENDIF
         ENDIF
         DELFQ(NV,NF)=DEL
C
         DO NV=2,NVQMAX(NF)
            NP1=NPQ(NV,NF)
            NP2=NPQ(MOD(NV,NVQMAX(NF))+1,NF)
            DX=XPQ(NP2)-XPQ(NP1)
            DY=YPQ(NP2)-YPQ(NP1)
            DL=SQRT(DX*DX+DY*DY)
            DX=DX/DL
            DY=DY/DL
            COSL=DX*DX1+DY*DY1
            DELFL=DELFL1*COSL**2+DELFL2*(1.D0-COSL**2)
C
            NL=(DL-EPSL)/DELFL+1
            DELL=DL/NL
            DEL=DELFQ(NV,NF)
            IF(DEL.EQ.0.D0) THEN
               DEL=DELL
            ELSE
               IF(DELL.LT.DEL) THEN
                  DEL=DELL
                  DELFQ(NVPQ(NV,NF),NFPQ(NV,NF))=DELL
               ENDIF
            ENDIF
            DELFQ(NV,NF)=DEL
         ENDDO
      ENDIF
      ENDDO
C
C     +++++ DEFINE NODES AND ELEMENTS +++++
C
      NN=0
      NE=0
      DO NF=1,NFQM
      IF(NVQMAX(NF).GE.3) THEN
         NVL1=1
         NVR1=2
         DELH=DELFQ(1,NF)
         DO NV=2,(NVQMAX(NF)+1)/2
           NVL2=NVQMAX(NF)+2-NV
           NVR2=NV+1
           NVL=NVL2
           NVR=NVR1
           DELL=DELFQ(NVL,NF)
           DELR=DELFQ(NVR,NF)
           NPL1=NPQ(NVL1,NF)
           NPR1=NPQ(NVR1,NF)
           NPL2=NPQ(NVL2,NF)
           NPR2=NPQ(NVR2,NF)
           XPL1=XPQ(NPL1)
           YPL1=YPQ(NPL1)
           XPR1=XPQ(NPR1)
           YPR1=YPQ(NPR1)
           XPL2=XPQ(NPL2)
           YPL2=YPQ(NPL2)
           XPR2=XPQ(NPR2)
           YPR2=YPQ(NPR2)
           DLENL=SQRT((XPL2-XPL1)**2+(YPL2-YPL1)**2)
           DLENR=SQRT((XPR2-XPR1)**2+(YPR2-YPR1)**2)
           NDELL=(DLENL+EPSL)/DELL
           NDELR=(DLENR+EPSL)/DELR
           IF(NDELL.GT.NDELR) THEN
              INDV=-1
              NDELV=NDELR
           ELSEIF(NDELL.LT.NDELR) THEN
              INDV=+1
              NDELV=NDELL
           ELSE
              INDV= 0
              NDELV=NDELL
           ENDIF
C
           DO MV=1,NDELV-1
              XL1=XPL1+(XPL2-XPL1)*(MV-1)/(NDELV-1)
              YL1=YPL1+(YPL2-YPL1)*(MV-1)/(NDELV-1)
              XL2=XPL1+(XPL2-XPL1)* MV   /(NDELV-1)
              YL2=YPL1+(YPL2-YPL1)* MV   /(NDELV-1)
              XR1=XPR1+(XPR2-XPR1)*(MV-1)/(NDELV-1)
              YR1=YPR1+(YPR2-YPR1)*(MV-1)/(NDELV-1)
              XR2=XPR1+(XPR2-XPR1)* MV   /(NDELV-1)
              YR2=YPR1+(YPR2-YPR1)* MV   /(NDELV-1)
C
              DLEN1=SQRT((XR1-XL1)**2+(YR1-YL1)**2)
              NDEL1=(DLEN1+EPSL)/DELH
              IF(MV.EQ.1) THEN
                 DO MH=1,NDEL1
                    NN=NN+1
                    XD(NN)=XL1+(XR1-XL1)*(MH-1)/NDEL1
                    YD(NN)=YL1+(YR1-YL1)*(MH-1)/NDEL1
                 ENDDO
              ENDDO
              NP1=NP-NDEL1
              NP2=NP
              DLEN2=SQRT((XR2-XL2)**2+(YR2-YL2)**2)
              NDEL2=(DLEN2+EPSL)/DELH
                 DO MH=1,NDEL2
                    NN=NN+1
                    XD(NN)=XL2+(XR2-XL2)*(MH-1)/NDEL2
                    YD(NN)=YL2+(YR2-YL2)*(MH-1)/NDEL2
                 ENDDO
C
              NP1INC=1
              NP2INC=1
  100         IF(NP1INC.EQ.NDEL1.AND.NP2INC.EQ.NDEL2) GOTO 110
                 IF(NP1INC.EQ.NDEL1) THEN
                    NE=NE+1
                    IELM(1,NE)=NP1+NP1INC
                    IELM(2,NE)=NP2+NP2INC+1
                    IELM(3,NE)=NP2+NP2INC
                    NP2INC=NP2INC+1
                 ELSEIF(NP2INC.EQ.NDEL2) THEN
                    NE=NE+1
                    IELM(1,NE)=NP1+NP1INC
                    IELM(2,NE)=NP1+NP1INC+1
                    IELM(3,NE)=NP2+NP2INC
                    NP1INC=NP1INC+1
                 ELSE
C
C
C
      IERR=0
      RETURN
C
 9000 RETURN
C
      END
