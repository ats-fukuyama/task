C     $Id$
C
C     ######### /TASK/WF/WFANT ########
C
C      ANTENNA DATA GENERATION PROGRAM
C
C     #################################
C
      SUBROUTINE WFANT
C
      INCLUDE 'wfcomm.inc'
      CHARACTER KID*1
C
    1 WRITE(6,601)
  601 FORMAT(1H ,'## INPUT: A/ANT  G/DRAW  P,V/PARM  S/SAVE  L/LOAD  ',
     &           'W/LIST  X/EXIT')
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL GUCPTL(KID)
C
      IF(KID.EQ.'A') THEN
         CALL WFDEFA
      ELSEIF(KID.EQ.'G') THEN
         CALL WFPLTA
      ELSEIF(KID.EQ.'P') THEN
         CALL WFPARM
      ELSEIF(KID.EQ.'V') THEN
         CALL WFVIEW
      ELSEIF(KID.EQ.'S') THEN
         CALL WFWANT
      ELSEIF(KID.EQ.'L') THEN
         CALL WFRANT
      ELSEIF(KID.EQ.'W') THEN
         DO NA=1,NAMAX
         WRITE(6,610) (N,XJ(N,NA),YJ(N,NA),JAELM(N,NA),N=1,JNUM(NA))
  610    FORMAT(1H ,I5,2F12.5,I10)
         ENDDO
      ELSEIF(KID.EQ.'X') THEN
         GOTO 9000
      ENDIF
      GOTO 1
C
 9000 RETURN
      END
C
C     ****** Define Antenna Data ******
C
      SUBROUTINE WFDEFA
C
      INCLUDE 'wfcomm.inc'
      CHARACTER KID*1
C
      DEGN=PI/180.D0
      IF(MODELB.EQ.5) THEN
         RKAP=BKAP
         RDEL=BDEL
      ENDIF
C
    1 WRITE(6,*) '## ENTER NUMBER OF ANTENNAS (1-',NAM,')'
      READ(5,*,ERR=1,END=9000) NAMAX
      IF(NAMAX.EQ.0) GOTO 9000
      IF(NAMAX.LT.1.OR.NAMAX.GT.NAM) GOTO 1
C
      NA=1
  100 CONTINUE
         JNUM0(NA)=0
C
    2    WRITE(6,*) '## ANTENNA NUMBER = ',NA
         WRITE(6,*) '## TYPE:  C/CIRCLE  A/ARC  P/POINTS  ',
     &              'H/Helical  L/Helicon  X/EXIT'
         READ(5,'(A1)',ERR=2,END=1) KID
         CALL GUCPTL(KID)
C
         IF(KID.EQ.'C') THEN
    3       WRITE(6,602) RD,RKAP,NJMAX
  602       FORMAT(1H ,'## RD,RKAP,NJMAX = ',2F10.3,I5)
            READ(5,*,ERR=3,END=2) RD,RKAP,NJMAX
C
            DTHETA=2.D0*PI/(NJMAX-1)
            DO 30 NJ=1,NJMAX
               THETA=DTHETA*(NJ-1)+1.D-5
               XJ0(NJ,NA)=RD*COS(THETA)
               YJ0(NJ,NA)=RD*SIN(THETA)*RKAP
   30       CONTINUE
            JNUM0(NA)=NJMAX
            RTJ0(1,NA)=0.D0
            RTJ0(2,NA)=0.D0
            NTYPJ0(NA)=-1
C
         ELSEIF(KID.EQ.'A') THEN
    4       WRITE(6,603) THJ1,THJ2,RD,RKAP,NJMAX
  603       FORMAT(1H ,'## THJ1,THJ2 = ',2F10.3/
     &             1H ,'## RD,RKAP,NJMAX = ',2F10.3,I5)
            READ(5,*,ERR=4,END=2) THJ1,THJ2,RD,RKAP,NJMAX
C
            THETA=DEGN*THJ1
            XJ0(1,NA)=1.5D0*RD*COS(THETA)
            YJ0(1,NA)=1.5D0*RD*SIN(THETA)*RKAP
            THJ0(1,NA)=THJ1
            THJ0(1,NA)=THJ2
            DTHETA=(THJ2-THJ1)/(NJMAX-3)
            DO 40 NJ=2,NJMAX-1
               THETA=DEGN*(DTHETA*(NJ-2)+THJ1)
               XJ0(NJ,NA)=RD*COS(THETA)
               YJ0(NJ,NA)=RD*SIN(THETA)*RKAP
   40       CONTINUE
            THETA=DEGN*THJ2
            XJ0(NJMAX,NA)=1.5D0*RD*COS(THETA)
            YJ0(NJMAX,NA)=1.5D0*RD*SIN(THETA)*RKAP
            JNUM0(NA)=NJMAX
            NTYPJ0(NA)=-1
C
         ELSEIF(KID.EQ.'P') THEN
    5       WRITE(6,*) '## NUMBER OF POINTS : NJMAX=',NJMAX
            READ(5,*,ERR=5,END=2) NJMAX
            DO 50 NJ=1,NJMAX
    6          WRITE(6,*) '## NO.',NJ,': X,Y ?'
               READ(5,*,ERR=6,END=5) X,Y
               XJ0(NJ,NA)=X
               YJ0(NJ,NA)=Y
   50       CONTINUE
            JNUM0(NA)=NJMAX
            NTYPJ0(NA)=-1
C
         ELSEIF(KID.EQ.'H') THEN
    7       WRITE(6,604) ZJH1,ZJH2,RD,RTJ1,RTJ2,NTYPJH
  604       FORMAT(1H ,'## ZJH1,ZJH2 = ',2F10.3/
     &             1H ,'## RD,RTJ1,RTJ2,NTYPJH = ',3F10.3,I5)
            READ(5,*,ERR=7,END=2) ZJH1,ZJH2,RD,RTJ1,RTJ2,NTYPJH
C
            RDX=MAX(ABS(XDMIN),ABS(XDMAX))+0.01D0
C
            RTJ0(1,NA)=RTJ1
            RTJ0(2,NA)=RTJ2
            IF(NTYPJH.EQ.0) THEN
               NJMAX=2
               XJ0(1,NA)=RD
               YJ0(1,NA)=ZJH1
               XJ0(2,NA)=RD
               YJ0(2,NA)=ZJH2
            ELSEIF(NTYPJH.EQ.1) THEN
               NJMAX=3
               XJ0(1,NA)=RDX
               YJ0(1,NA)=ZJH1
               XJ0(2,NA)=RD
               YJ0(2,NA)=ZJH1
               XJ0(3,NA)=RD
               YJ0(3,NA)=ZJH2
            ELSEIF(NTYPJH.EQ.2) THEN
               NJMAX=3
               XJ0(1,NA)=RD
               YJ0(1,NA)=ZJH1
               XJ0(2,NA)=RD
               YJ0(2,NA)=ZJH2
               XJ0(3,NA)=RDX
               YJ0(3,NA)=ZJH2
            ELSEIF(NTYPJH.EQ.3) THEN
               NJMAX=4
               XJ0(1,NA)=RDX
               YJ0(1,NA)=ZJH1
               XJ0(2,NA)=RD
               YJ0(2,NA)=ZJH1
               XJ0(3,NA)=RD
               YJ0(3,NA)=ZJH2
               XJ0(4,NA)=RDX
               YJ0(4,NA)=ZJH2
            ELSEIF(NTYPJH.EQ.10) THEN
               NJMAX=2
               XJ0(1,NA)=RD
               YJ0(1,NA)=ZJH1
               XJ0(2,NA)=RD
               YJ0(2,NA)=ZJH2
            ENDIF
            JNUM0(NA)=NJMAX
            NTYPJ0(NA)=NTYPJH
C
         ELSEIF(KID.EQ.'L') THEN
    8       WRITE(6,605) ZJH1,ZJH2,ZJH3,ZJH4,RTJ1,RTJ2,RD
  605       FORMAT(1H ,'## ZJH1,ZJH2,ZJH3,ZJH4 = ',4F10.3/
     &             1H ,'## RTJ1,RTJ2,RD        = ',3F10.3)
            READ(5,*,ERR=8,END=2) ZJH1,ZJH2,ZJH3,ZJH4,
     &                            RTJ1,RTJ2,RD
C
            RDX=MAX(ABS(XDMIN),ABS(XDMAX))+0.01D0
C
            NTYPJH=11
            NAMAX=4
            DO NA=2,NAMAX
               AJ(NA)=AJ(1)
               APH(NA)=APH(1)
            END DO

            NA=1
            JNUM0(NA)=1
            XJ0(1,NA)=RD
            YJ0(1,NA)=ZJH1
            RTJ0(1,NA)=RTJ1
            RTJ0(2,NA)=RTJ2
            NTYPJ0(NA)=NTYPJH

            NA=2
            JNUM0(NA)=1
            XJ0(1,NA)=RD
            YJ0(1,NA)=ZJH2
            RTJ0(1,NA)=RTJ1
            RTJ0(2,NA)=RTJ2
            NTYPJ0(NA)=NTYPJH

            NA=3
            JNUM0(NA)=2
            XJ0(1,NA)=RD
            YJ0(1,NA)=ZJH1
            XJ0(2,NA)=RD
            YJ0(2,NA)=ZJH2
            RTJ0(1,NA)=RTJ1
            RTJ0(2,NA)=RTJ2
            NTYPJ0(NA)=NTYPJH

            NA=4
            JNUM0(NA)=4
            XJ0(1,NA)=RDX
            YJ0(1,NA)=ZJH3
            XJ0(2,NA)=RD
            YJ0(2,NA)=ZJH3
            XJ0(3,NA)=RD
            YJ0(3,NA)=ZJH4
            XJ0(4,NA)=RDX
            YJ0(4,NA)=ZJH4
            RTJ0(1,NA)=RTJ1
            RTJ0(2,NA)=RTJ2
            NTYPJ0(NA)=NTYPJH

         ELSEIF(KID.EQ.'X') THEN
            GOTO 1
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID: ',KID
            GOTO 2
         ENDIF            
C
         NA=NA+1
         IF(NA.LE.NAMAX) GOTO 100

      CALL MODANT(IERR)
      IF(IERR.NE.0) GOTO 9000
C
 9000 RETURN
      END
C
C     ****** Draw Antenna ******
C
      SUBROUTINE WFPLTA
C
      INCLUDE 'wfcomm.inc'
C
      CALL WFGVEW(1,1,GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &                GXMIN,GXMAX,GYMIN,GYMAX)
C
      CALL PAGES
      CALL WFPRMJ
      CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,
     &            GXMIN,GXMAX,GYMIN,GYMAX)
C
      CALL SETLIN(0,0,4)
      IF(NDRAWA.LE.1) THEN
         CALL WFGBDY
      ELSE
         NTEMP=NDRAWD
         NDRAWD=NDRAWA-1
         CALL WFGELM
         NDRAWD=NTEMP
      ENDIF
C
      CALL SETLIN(0,0,5)
      CALL WFGPLA
C
      CALL SETLIN(0,0,6)
      CALL WFGANT
C
      CALL PAGEE
      RETURN
      END
C
C     ****** Draw Antenna Paramters ******
C
      SUBROUTINE WFPRMJ
C
      INCLUDE 'wfcomm.inc'
C
      GXMIN=20.
      GYMAX=17.
      CALL SETCHS(0.3,0.)
      GDY=1.5*0.3
C
      GXL=GXMIN
      GYL=GYMAX
      CALL MOVE(GXL,GYL)
      CALL TEXT('NNOD =',6)
      CALL NUMBI(NNOD,'(I5)',5)
C
      GYL=GYL-GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('ELMT =',6)
      CALL NUMBI(NELM,'(I5)',5)
C
      GYL=GYL-GDY
      CALL MOVE(GXL,GYL)
      CALL TEXT('JNUM =',6)
      GXL=GXL+6*0.3
      DO 10 NA=1,NAMAX
         CALL MOVE(GXL,GYL)
         CALL NUMBI(JNUM(NA),'(I5)',5)
         GYL=GYL-GDY
   10 CONTINUE
      RETURN
      END
