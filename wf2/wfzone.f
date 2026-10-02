C     $Id$
C
C     ********** F.E.M. ZONING UTILITY **********
C
      SUBROUTINE WFZONE
C
      USE libchar
      INCLUDE 'wfcomm.inc'
      CHARACTER KID*1
C
    1 WRITE(6,601) 
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL toupper(KID)
C
      IF(KID.EQ.'Z') THEN
         CALL WFXZON
         CALL WFSETZ
      ELSEIF(KID.EQ.'G') THEN
         CALL WFGZON
      ELSEIF(KID.EQ.'W') THEN
         CALL WFLZON
      ELSEIF(KID.EQ.'S') THEN
         CALL WFWZON
      ELSEIF(KID.EQ.'L') THEN
         CALL WFRZON
      ELSEIF(KID.EQ.'P') THEN
         CALL WFPARM
         CALL WFSETZ
      ELSEIF(KID.EQ.'V') THEN
         CALL WFVIEW
      ELSEIF(KID.EQ.'X') THEN
         GOTO 9000
      ENDIF
      GOTO 1
C
 9000 RETURN
  601 FORMAT(1H ,'## INPUT: Z/ZONE  G/DRAW  P,V/PARM  S/SAVE  L/LOAD  ',
     &           'W/LIST  X/EXIT')
      END
C
C     ****** Input Element Data ******
C
      SUBROUTINE WFXZON
C
      USE libchar
      INCLUDE 'wfcomm.inc'
      CHARACTER KID*1
C
    1 WRITE(6,*) '## SELECT: D/ELMENT  E,M/MEDIUM  B/BOUNDARY  ',
     &           'W,P/VARPHI  V/VIEW  X/EXIT'
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL toupper(KID)
C
      IF(KID.EQ.'D') THEN
   10    WRITE(6,*) '## INPUT: NDMAX'
         READ(5,*,ERR=10,END=1) NDMAX
         IF(NDMAX.LE.0) GOTO 1
   11    CONTINUE
            WRITE(6,*) '## INPUT: ND,ID,XMIN,XMAX,YMIN,YMAX ',
     &                 '(ND=0 for exit)'
            READ(5,*,ERR=11,END=10) ND,ID,X1,X2,Y1,Y2
            IF(ND.LE.0) GOTO 1
            IF(ND.LE.NDMAX) THEN
               IDZONE(ND)=ID
               XYZONE(1,ND)=X1
               XYZONE(2,ND)=X2
               XYZONE(3,ND)=Y1
               XYZONE(4,ND)=Y2
            ELSE
               WRITE(6,*) 'XX ND IS LARGER THEN NDMAX'
            ENDIF
         GOTO 11
      ELSEIF(KID.EQ.'E') THEN
   20    WRITE(6,*) '## INPUT: NMMAX (NM=1:Vacuum)'
         READ(5,*,ERR=20,END=1) NMMAX
         IF(NMMAX.LE.1) GOTO 1
   21    CONTINUE
            WRITE(6,*) '## INPUT: NM,EPSD ',
     &                 '(NM=0 for exit)'
            READ(5,*,ERR=21,END=1) NM,EPSD
            IF(NM.LE.1) GOTO 1
            IF(NM.LE.NMMAX) THEN
               EPSDM(NM)=EPSD
               RMUDM(NM)=1.D0
               SIGDM(NM)=0.D0
            ELSE
               WRITE(6,*) 'XX NM IS LARGER THEN NMMAX'
            ENDIF
         GOTO 21
      ELSEIF(KID.EQ.'M') THEN
   25    WRITE(6,*) '## INPUT: NMMAX (NM=1:Vacuum)'
         READ(5,*,ERR=25,END=1) NMMAX
         IF(NMMAX.LE.1) GOTO 1
   26    CONTINUE
            WRITE(6,*) '## INPUT: NM,EPSD,RMUD,SIGD ',
     &                 '(NM=0 for exit)'
            READ(5,*,ERR=26,END=1) NM,EPSD,RMUD,SIGD
            IF(NM.LE.1) GOTO 1
            IF(NM.LE.NMMAX) THEN
               EPSDM(NM)=EPSD
               RMUDM(NM)=RMUD
               SIGDM(NM)=SIGD
            ELSE
               WRITE(6,*) 'XX NM IS LARGER THEN NMMAX'
            ENDIF
         GOTO 26
      ELSEIF(KID.EQ.'A') THEN
   30    WRITE(6,*) '## INPUT: NBWMAX'
         READ(5,*,ERR=30,END=1) NBWMAX
         IF(NBWMAX.LE.0) GOTO 1
   31    CONTINUE
            WRITE(6,*) '## INPUT: NBW,ID,XMIN,XMAX,YMIN,YMAX ',
     &                 '(NBW=0 for exit)'
            READ(5,*,ERR=31,END=1) NBW,ID,X1,X2,Y1,Y2
            IF(NBW.LE.0) GOTO 1
            IF(NBW.LE.NBWMAX) THEN
               IDBDYW(NBW)=ID
               XYBDYW(1,NBW)=X1
               XYBDYW(2,NBW)=X2
               XYBDYW(3,NBW)=Y1
               XYBDYW(4,NBW)=Y2
            ELSE
               WRITE(6,*) 'XX NBW IS LARGER THEN NBWMAX'
            ENDIF
         GOTO 31
      ELSEIF(KID.EQ.'B') THEN
   40    WRITE(6,*) '## INPUT: NBPMAX'
         READ(5,*,ERR=40,END=1) NBPMAX
         IF(NBPMAX.LE.0) GOTO 1
   41    CONTINUE
            WRITE(6,*) '## INPUT: NBP,ID,XMIN,XMAX,YMIN,YMAX ',
     &                 '(NBP=0 for exit)'
            READ(5,*,ERR=41,END=1) NBP,ID,X1,X2,Y1,Y2
            IF(NBP.LE.0) GOTO 1
            IF(NBP.LE.NBPMAX) THEN
               IDBDYP(NBP)=ID
               XYBDYP(1,NBP)=X1
               XYBDYP(2,NBP)=X2
               XYBDYP(3,NBP)=Y1
               XYBDYP(4,NBP)=Y2
            ELSE
               WRITE(6,*) 'XX NBP IS LARGER THEN NBPMAX'
            ENDIF
         GOTO 41
      ELSEIF(KID.EQ.'W') THEN
   50    WRITE(6,*) '## INPUT: NVWMAX'
         READ(5,*,ERR=50,END=1) NVWMAX
         IF(NVWMAX.LE.0) GOTO 1
   51    CONTINUE
            WRITE(6,*) '## INPUT: NVW,XMIN,XMAX,YMIN,YMAX ',
     &                 '(NVW=0 for exit)'
            READ(5,*,ERR=51,END=1) NVW,X1,X2,Y1,Y2
            IF(NVW.LE.0) GOTO 1
            IF(NVW.LE.NVWMAX) THEN
               XYVARW(1,NVW)=X1
               XYVARW(2,NVW)=X2
               XYVARW(3,NVW)=Y1
               XYVARW(4,NVW)=Y2
            ELSE
               WRITE(6,*) 'XX NVW IS LARGER THEN NVWMAX'
            ENDIF
         GOTO 51
      ELSEIF(KID.EQ.'P') THEN
   60    WRITE(6,*) '## INPUT: NVPMAX'
         READ(5,*,ERR=60,END=1) NVPMAX
         IF(NVPMAX.LE.0) GOTO 1
   61    CONTINUE
            WRITE(6,*) '## INPUT: NVP,XMIN,XMAX,YMIN,YMAX ',
     &                 '(NVP=0 for exit)'
            READ(5,*,ERR=61,END=1) NVP,X1,X2,Y1,Y2
            IF(NVP.LE.0) GOTO 1
            IF(NVP.LE.NVPMAX) THEN
               XYVARP(1,NVP)=X1
               XYVARP(2,NVP)=X2
               XYVARP(3,NVP)=Y1
               XYVARP(4,NVP)=Y2
            ELSE
               WRITE(6,*) 'XX NVP IS LARGER THEN NVPMAX'
            ENDIF
         GOTO 61
      ELSEIF(KID.EQ.'V') THEN
         DO ND=1,NDMAX
            WRITE(6,601) ND,IDZONE(ND),(XYZONE(I,ND),I=1,4)
         ENDDO
         DO NM=1,NMMAX
            WRITE(6,602) NM,EPSDM(NM),RMUDM(NM),SIGDM(NM)
         ENDDO
         DO NBW=1,NBWMAX
            WRITE(6,603) NBW,IDBDYW(NBW),(XYBDYW(I,NBW),I=1,4)
         ENDDO
         DO NBP=1,NBPMAX
            WRITE(6,604) NBP,IDBDYP(NBP),(XYBDYP(I,NBP),I=1,4)
         ENDDO
         DO NVW=1,NVWMAX
            WRITE(6,605) NVW,(XYVARW(I,NVW),I=1,4)
         ENDDO
         DO NVP=1,NVPMAX
            WRITE(6,606) NVP,(XYVARP(I,NVP),I=1,4)
         ENDDO
  601    FORMAT(' ND: ',I3,I3,4F12.4)
  602    FORMAT(' NM: ',I3,3X,1P3E12.4)
  603    FORMAT(' NBW:',I3,I3,4F12.4)
  604    FORMAT(' NBP:',I3,I3,4F12.4)
  605    FORMAT(' NVW:',I3,3X,4F12.4)
  606    FORMAT(' NVP:',I3,3X,4F12.4)
      ELSEIF(KID.EQ.'X') THEN
         GOTO 9000
      ENDIF
      GOTO 1
C
 9000 RETURN
      END
C
C
C     ****** Draw Zone Data ******
C
      SUBROUTINE WFGZON
C
      INCLUDE 'wfcomm.inc'
C
      CALL WFGVEW(1,1,GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &                GXMIN,GXMAX,GYMIN,GYMAX)
C
      CALL PAGES
      CALL MOVE(1.0,17.0)
      CALL TEXT('KNODW',5)
      CALL WFPRME
      CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
      CALL WFGELM
      CALL WFGNOD(0)
      CALL PAGEE
C
      CALL PAGES

      CALL MOVE(1.0,17.0)
      CALL TEXT('KNODP',5)
      CALL WFPRME
      CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
      CALL WFGELM
      CALL WFGNOD(1)
      CALL PAGEE
      RETURN
      END
C
C     ****** List Element Data ******
C
      SUBROUTINE WFLZON
C
      INCLUDE 'wfcomm.inc'
C
      WRITE(6,100) NNOD,
     &             (IN,KNODW(IN),KNODP(IN),XD(IN),YD(IN),IN=1,NNOD)
  100 FORMAT(1H ,'NODE DATA',7X,'NNOD=',I6/
     &       1H ,2(3X,'IN KNODW/P',8X,'XD',11X,'YD',2X)/
     &      (1H ,2(I6,I3,I3,1P2E13.4)))
C
      WRITE(6,200) NELM,(IE,KELM(IE),(IELM(J,IE),J=1,3),
     &                    IE=1,NELM)
  200 FORMAT(1H ,'ELEMENT DATA',5X,'NELM=',I6/
     &       1H ,2(3X,'IE ','KELM',17X)/
     &      (1H ,2(I6,I3,'(',3I5,')')))
C
      WRITE(6,300) NBDYP,(IBDYP(IB),IB=1,NBDYP)
  300 FORMAT(1H ,'BOUNDARY DATA',4X,'NBDYP=',I6/
     &      (1H ,10(I6)))
C
      WRITE(6,400) NBDYW,(IBDYW(I),I=1,NBDYW)
  400 FORMAT(1H ,'BOUNDARY DATA',4X,'NBDYW=',I6/
     &      (1H ,10(I6)))
C
      RETURN
      END
C
C     ****** Set ZONE INFO ******
C
      SUBROUTINE WFSETZ
C
      INCLUDE 'wfcomm.inc'
      DIMENSION IBDYL(NBDYM)
C
      DO IE=1,NELM
         KELM(IE)=0
      ENDDO
C
      CALL SETBDY(0,IBDYL,NBDYL)
      NBDYW=NBDYL
      DO NB=1,NBDYW
         IBDYW(NB)=IBDYL(NB)
      ENDDO
C      DO NB=1,NBDYW
C      WRITE(6,'(A,2I4,1P2E12.4)') 'NB,IBDYW,XD,YD=',NB,IBDYW(NB),
C     &                             XD(IBDYW(NB)),YD(IBDYW(NB))
C      ENDDO
C
      DO ND=1,NDMAX
         XMIN=XYZONE(1,ND)
         XMAX=XYZONE(2,ND)
         YMIN=XYZONE(3,ND)
         YMAX=XYZONE(4,ND)
         DO IE=1,NELM
            IN1=IELM(1,IE)
            IN2=IELM(2,IE)
            IN3=IELM(3,IE)
            X1=XD(IN1)
            Y1=YD(IN1)
            X2=XD(IN2)
            Y2=YD(IN2)
            X3=XD(IN3)
            Y3=YD(IN3)
            IF((X1-XMIN)*(X1-XMAX).LE.0.D0.AND.
     &         (Y1-YMIN)*(Y1-YMAX).LE.0.D0.AND.
     &         (X2-XMIN)*(X2-XMAX).LE.0.D0.AND.
     &         (Y2-YMIN)*(Y2-YMAX).LE.0.D0.AND.
     &         (X3-XMIN)*(X3-XMAX).LE.0.D0.AND.
     &         (Y3-YMIN)*(Y3-YMAX).LE.0.D0) THEN
               KELM(IE)=IDZONE(ND)
            ENDIF
         ENDDO
      ENDDO
C
      DO IN=1,NNOD
         KNODW(IN)=0
         KNODP(IN)=0
      ENDDO
      DO IB=1,NBDYW
         KNODW(IBDYW(IB))=2
         KNODP(IBDYW(IB))=1
      ENDDO
C
      IF(MODELS.EQ.1) THEN
         IN=IBDYW(NBDYW)
         IF(XD(IN).LT.1.D-6) THEN
            IND=1
         ELSE
            IND=0
         ENDIF
         DO IB=1,NBDYW
            IN=IBDYW(IB)
            IF(XD(IN).LT.1.D-6) THEN
               IF(IND.EQ.0) THEN
                  KNODW(IN)=3
                  IND=1
               ELSE
                  KNODW(IN)=1
                  KNODP(IN)=0
               ENDIF
            ELSE
               IF(IND.EQ.1) THEN
                  IF(IB.GE.1) THEN
                     KNODW(IBDYW(IB-1))=3
                     KNODP(IBDYW(IB-1))=1
                  ENDIF
                  IND=0
               ENDIF
            ENDIF
         ENDDO
      ENDIF
C
      DO NBW=1,NBWMAX
         XMIN=XYBDYW(1,NBW)
         XMAX=XYBDYW(2,NBW)
         YMIN=XYBDYW(3,NBW)
         YMAX=XYBDYW(4,NBW)
         DO IN=1,NNOD
            X1=XD(IN)
            Y1=YD(IN)
            IF((X1-XMIN)*(X1-XMAX).LE.0.D0.AND.
     &         (Y1-YMIN)*(Y1-YMAX).LE.0.D0) THEN
               IF(KNODW(IN).EQ.2) THEN
                  KNODW(IN)=2*IDBDYW(NBW)
               ELSEIF(KNODW(IN).EQ.3) THEN
                  KNODW(IN)=2*IDBDYW(NBW)+1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
C
      CALL SETBDY(0,IBDYL,NBDYL)
      NBDYP=NBDYL
      DO NB=1,NBDYP
         IBDYP(NB)=IBDYL(NB)
      ENDDO
C
      DO IE=1,NELM
         IF(KELM(IE).NE.0) THEN
            IN=IELM(1,IE)
            IF(KNODP(IN).EQ.0) KNODP(IN)=10
            IN=IELM(2,IE)
            IF(KNODP(IN).EQ.0) KNODP(IN)=10
            IN=IELM(3,IE)
            IF(KNODP(IN).EQ.0) KNODP(IN)=10
         ENDIF
      ENDDO
C
      DO NBP=1,NBPMAX
         XMIN=XYBDYP(1,NBP)
         XMAX=XYBDYP(2,NBP)
         YMIN=XYBDYP(3,NBP)
         YMAX=XYBDYP(4,NBP)
         DO IN=1,NNOD
            X1=XD(IN)
            Y1=YD(IN)
            IF((X1-XMIN)*(X1-XMAX).LE.0.D0.AND.
     &         (Y1-YMIN)*(Y1-YMAX).LE.0.D0) THEN
               KNODP(IN)=IDBDYP(NBP)
            ENDIF
         ENDDO
      ENDDO
C
C      DO IB=1,NBDYP
C         IN=IBDYP(IB)
C         IF(KNODW(IN).GE.0.AND.KNODP(IN).EQ.0) THEN
C            WRITE(6,601) IN,XD(IN),YD(IN)
C  601       FORMAT('XX BOUNDARY DEF ERROR: IN,XD,YD =',I5,2F10.4)
C         ENDIF
C      ENDDO
C
C
      NELMP=0
      DO IE=1,NELM
         IF(KELM(IE).EQ.0) THEN
            NELMP=NELMP+1
            IELMP(1,NELMP)=IELM(1,IE)
            IELMP(2,NELMP)=IELM(2,IE)
            IELMP(3,NELMP)=IELM(3,IE)
            IF(KNODP(IELM(1,IE)).GE.10) 
     &           WRITE(6,*) 'XX IELMP DEF ERROR: IN,KNODP=',
     &           IELM(1,IE),KNODP(IELM(1,IE))
            IF(KNODP(IELM(2,IE)).GE.10) 
     &           WRITE(6,*) 'XX IELMP DEF ERROR: IN,KNODP=',
     &           IELM(2,IE),KNODP(IELM(2,IE))
            IF(KNODP(IELM(3,IE)).GE.10) 
     &           WRITE(6,*) 'XX IELMP DEF ERROR: IN,KNODP=',
     &           IELM(3,IE),KNODP(IELM(3,IE))
         ENDIF
      ENDDO
C
      RETURN
      END

C
C     ****** Set Bounary Node Array ******
C
      SUBROUTINE SETBDY(KEMAX,IBDYL,NBDYL)
C
      INCLUDE 'wfcomm.inc'
      DIMENSION IBDYL(NBDYM)
C
      DO IE=1,NELM
         IE0=IE
         IF(KELM(IE).LE.KEMAX) THEN
            IN1=IELM(1,IE0)
            IN2=IELM(2,IE0)
            IN3=IELM(3,IE0)
            CALL EFINDK(IE0,IN1,IN2,IE1,KEMAX)
            IF(IE1.EQ.0) THEN
               INS=IN1
               GOTO 20
            ENDIF
            CALL EFINDK(IE0,IN2,IN3,IE1,KEMAX)
            IF(IE1.EQ.0) THEN
               INS=IN2
               GOTO 20
            ENDIF
            CALL EFINDK(IE0,IN3,IN1,IE1,KEMAX)
            IF(IE1.EQ.0) THEN
               INS=IN3
               GOTO 20
            ENDIF
         ENDIF
      ENDDO
      GOTO 9100
C
   20 INPRE=0
      IN=INS
      IB=1
      IBDYL(IB)=IN
C
C      WRITE(6,*) 'IB,IN,IE =',IB,IN,IE
C
C     FIND NODE IN1 IN ELEMENT IE INCLUDING IN
C                   BUT NOT SHARED BY ANY OTHER ELEMENT
C
   30 CONTINUE
         DO I=1,MAX(IE0,NELM-IE0)
         DO J=1,2
           IF(J.EQ.1) THEN
               IF(IB.EQ.1) THEN
                  IE=IE0-I+1
               ELSE
                  IE=IE0-I
               ENDIF
            ELSE
               IE=IE0+I
            ENDIF

            CALL GUFLSH

            IF(IE.GE.1.AND.IE.LE.NELM) THEN
               IF(KELM(IE).LE.KEMAX) THEN
               IF(IELM(1,IE).EQ.IN) THEN
                  IN1=IELM(2,IE)
                  CALL EFINDK(IE,IN,IN1,IE1,KEMAX)
                  IF(IE1.EQ.0) GOTO 100
                  IN1=IELM(3,IE)
                  CALL EFINDK(IE,IN,IN1,IE1,KEMAX)
                  IF(IE1.EQ.0) GOTO 100
               ELSEIF(IELM(2,IE).EQ.IN) THEN
                  IN1=IELM(3,IE)
                  CALL EFINDK(IE,IN,IN1,IE1,KEMAX)
                  IF(IE1.EQ.0) GOTO 100
                  IN1=IELM(1,IE)
                  CALL EFINDK(IE,IN,IN1,IE1,KEMAX)
                  IF(IE1.EQ.0) GOTO 100
               ELSEIF(IELM(3,IE).EQ.IN) THEN
                  IN1=IELM(1,IE)
                  CALL EFINDK(IE,IN,IN1,IE1,KEMAX)
                  IF(IE1.EQ.0) GOTO 100
                  IN1=IELM(2,IE)
                  CALL EFINDK(IE,IN,IN1,IE1,KEMAX)
                  IF(IE1.EQ.0) GOTO 100
               ENDIF
            ENDIF
            ENDIF
         ENDDO
         ENDDO
            IE=IE0
            IF(IELM(1,IE).EQ.IN) THEN
               IN1=IELM(2,IE)
               IF(IN1.NE.INPRE) THEN
                  CALL EFINDK(IE,IN,IN1,IE1,KEMAX)
                  IF(IE1.EQ.0) GOTO 100
               ELSE
                  IN1=IELM(3,IE)
                  CALL EFINDK(IE,IN,IN1,IE1,KEMAX)
                  IF(IE1.EQ.0) GOTO 100
               ENDIF
            ELSEIF(IELM(2,IE).EQ.IN) THEN
               IN1=IELM(3,IE)
               IF(IN1.NE.INPRE) THEN
                  CALL EFINDK(IE,IN,IN1,IE1,KEMAX)
                  IF(IE1.EQ.0) GOTO 100
               ELSE
                  IN1=IELM(1,IE)
                  CALL EFINDK(IE,IN,IN1,IE1,KEMAX)
                  IF(IE1.EQ.0) GOTO 100
               ENDIF
            ELSEIF(IELM(3,IE).EQ.IN) THEN
               IN1=IELM(1,IE)
               IF(IN1.NE.INPRE) THEN
                  CALL EFINDK(IE,IN,IN1,IE1,KEMAX)
                  IF(IE1.EQ.0) GOTO 100
               ELSE
                  IN1=IELM(2,IE)
                  CALL EFINDK(IE,IN,IN1,IE1,KEMAX)
                  IF(IE1.EQ.0) GOTO 100
               ENDIF
            ENDIF
         GOTO 9200
C
  100    CONTINUE
         IF(IB.GE.3.AND.IN1.EQ.INS) GOTO 1000
         INPRE=IN
         IN=IN1
         IF(IB.GE.NBDYM) GOTO 9300
         IB=IB+1
         IBDYL(IB)=IN
         IE0=IE
C         WRITE(6,*) 'IB,IN,IE =',IB,IN,IE
      GOTO 30
C
 1000 NBDYL=IB
      RETURN
C
 9100 WRITE(6,*) 'XX ERROR IN SETBDY: NO ELEMENT ON BOUNDARY'
      NBDYL=0
      RETURN
C
 9200 WRITE(6,*) 'XX ERROR IN SETBDY: CANNOT FIND NEW IN: IN,IE=',
     &           IN,IE
      NBDYL=0
      RETURN
C
 9300 WRITE(6,*) 'XX ERROR IN SETBDY: IB EXCEEDS NBDYM: IB=',IB
      NBDYL=0
      RETURN
      END      
C
C     ******* FIND ELEMENT INCLUDING NODES NJ,NK *******
C
      SUBROUTINE EFINDK(IES,N1,N2,IE,KEMAX)
C
      INCLUDE 'wfcomm.inc'
C
      IF(IES.LT.0.OR.IES.GT.NELM+1) GOTO 9000
      DO 10 I=1,MAX(NELM-IES,IES)
         IDELT=I
         DO 20 J=1,2
            IDELT=-IDELT
            IE    =IES+IDELT
            IF(IE.GE.1.AND.IE.LE.NELM) THEN
            IF(KELM(IE).LE.KEMAX) THEN
               DO 30 K=1,3
                  NE1=IELM(K,IE)
                  IF(NE1.EQ.N1) THEN
                     DO 40 L=1,2
                        LL=MOD(K+L-1,3)+1
                        NE2=IELM(LL,IE)
                        IF(NE2.EQ.N2) RETURN
   40                CONTINUE
                  ENDIF
   30          CONTINUE
            ENDIF
            ENDIF
   20    CONTINUE
   10 CONTINUE
C
      IE=0
      RETURN
C
 9000 IE=0
      RETURN
      END
