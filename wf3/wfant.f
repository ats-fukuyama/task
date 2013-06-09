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
C
      WRITE(6,*) '--- SETBDY start ---'
      CALL SETBDY(IERR)
      IF(IERR.NE.0) RETURN
C
      WRITE(6,*) '--- SETSID start ---'
      CALL SETSID(IERR)
      IF(IERR.NE.0) RETURN
C
    1 WRITE(6,601)
  601 FORMAT(' ','## INPUT: A/ANT  G/DRAW  P,V/PARM  S/SAVE  L/LOAD  ',
     &           'W/LIST  X/EXIT')
      READ(5,'(A1)',ERR=1,END=9000) KID
      CALL GUCPTL(KID)
C
      IF(KID.EQ.'A') THEN
         CALL WFDEFA
      ELSEIF(KID.EQ.'G') THEN
         CALL WFPLTA
      ELSEIF(KID.EQ.'P') THEN
         CALL WFPARM(KID)
      ELSEIF(KID.EQ.'V') THEN
         CALL WFVIEW
      ELSEIF(KID.EQ.'S') THEN
         CALL WFWANT
      ELSEIF(KID.EQ.'L') THEN
         CALL WFRANT
      ELSEIF(KID.EQ.'W') THEN
         DO NA=1,NAMAX
         WRITE(6,610) (N,XJ0(N,NA),YJ0(N,NA),ZJ0(N,NA),
     &                 N=1,JNUM0(NA))
  610    FORMAT(' ',I5,3F12.5)
         ENDDO
         DO NA=1,NAMAX
         WRITE(6,611) (N,XJ(N,NA),YJ(N,NA),ZJ(N,NA),JELMT(N,NA),
     &                 N=1,JNUM(NA))
  611    FORMAT(' ',I5,3F12.5,I10)
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
C
    1 WRITE(6,*) '## ENTER NUMBER OF ANTENNAS (1-',NAM,')'
      READ(5,*,ERR=1,END=9000) NAMAX
      IF(NAMAX.EQ.0) GOTO 9000
      IF(NAMAX.LT.1.OR.NAMAX.GT.NAM) GOTO 1
C
      DO 100 NA=1,NAMAX
         JNUM0(NA)=0
C
    2    WRITE(6,*) '## ANTENNA NUMBER = ',NA
         WRITE(6,*) '## TYPE: C/CIRCLE  A/ARC  S/SPIRAL P/POINTS'//
     &              '  X/EXIT'
         READ(5,'(A1)',ERR=2,END=1) KID
         CALL GUCPTL(KID)
C
         IF(KID.EQ.'C') THEN
    3       WRITE(6,602) RD,NJMAX,ZANT
  602       FORMAT(' ','## RD,NJMAX,ZANT = ',F10.3,I5,F10.3)
            READ(5,*,ERR=3,END=2) RD,NJMAX,ZANT
C
            DTHETA=2.D0*PI/(NJMAX-1)
            DO 30 NJ=1,NJMAX
               THETA=DTHETA*(NJ-1)+1.D-5
               XJ0(NJ,NA)=RD*COS(THETA)
               YJ0(NJ,NA)=RD*SIN(THETA)
               ZJ0(NJ,NA)=ZANT
   30       CONTINUE
            JNUM0(NA)=NJMAX
C
         ELSEIF(KID.EQ.'A') THEN
    4       WRITE(6,603) THETJ1,THETJ2,RD,NJMAX,ZANT
  603       FORMAT(' ','## THETJ1,THETJ2 = ',2F10.3/
     &             ' ','## RD,NJMAX = ',F10.3,I5/
     &             ' ','## ZANT = ',F10.3)
            READ(5,*,ERR=4,END=2) THETJ1,THETJ2,RD,NJMAX
C
            THETA=DEGN*THETJ1
            XJ0(1,NA)=1.5D0*RD*COS(THETA)
            YJ0(1,NA)=1.5D0*RD*SIN(THETA)
            ZJ0(1,NA)=ZANT
            DTHETA=(THETJ2-THETJ1)/(NJMAX-3)
            DO 40 NJ=2,NJMAX-1
               THETA=DEGN*(DTHETA*(NJ-2)+THETJ1)
               XJ0(NJ,NA)=RD*COS(THETA)
               YJ0(NJ,NA)=RD*SIN(THETA)
               ZJ0(NJ,NA)=ZANT
   40       CONTINUE
            THETA=DEGN*THETJ2
            XJ0(NJMAX,NA)=1.5D0*RD*COS(THETA)
            YJ0(NJMAX,NA)=1.5D0*RD*SIN(THETA)
            ZJ0(NJMAX,NA)=ZANT
            JNUM0(NA)=NJMAX
C
         ELSEIF(KID.EQ.'S') THEN
    5       WRITE(6,604) THETS1,THETS2,RD1,RD2,NJMAX,ZANT,ZWALL
  604       FORMAT(' ','## THETS1,THETJS2 = ',2F10.3/
     &             ' ','## RD1,RD2,NJMAX = ',2F10.3,I5/
     &             ' ','## ZANT,ZWALL = ',2F10.3)
            READ(5,*,ERR=5,END=2) THETS1,THETS2,RD1,RD2,NJMAX,
     &                            ZANT,ZWALL
C
            THETA=DEGN*THETS1
            XJ0(1,NA)=RD1*COS(THETA)
            YJ0(1,NA)=RD1*SIN(THETA)
            ZJ0(1,NA)=ZWALL
            DTHETA=(THETS2-THETS1)/(NJMAX-3)
            DRD   =(RD2   -RD1   )/(NJMAX-3)
            DO 45 NJ=2,NJMAX-1
               THETA=DEGN*(DTHETA*(NJ-2)+THETS1)
               RDL  =      DRD   *(NJ-2)+RD1
               XJ0(NJ,NA)=RDL*COS(THETA)
               YJ0(NJ,NA)=RDL*SIN(THETA)
               ZJ0(NJ,NA)=ZANT
   45       CONTINUE
            THETA=DEGN*THETS2
            XJ0(NJMAX,NA)=RD2*COS(THETA)
            YJ0(NJMAX,NA)=RD2*SIN(THETA)
            ZJ0(NJMAX,NA)=ZWALL
            JNUM0(NA)=NJMAX
C
         ELSEIF(KID.EQ.'P') THEN
    6       WRITE(6,*) '## NUMBER OF POINTS : NJMAX=',NJMAX
            READ(5,*,ERR=6,END=2) NJMAX
            DO 50 NJ=1,NJMAX
    7          WRITE(6,*) '## NO.',NJ,': X,Y,Z ?'
               READ(5,*,ERR=7,END=6) X,Y,Z
               XJ0(NJ,NA)=X
               YJ0(NJ,NA)=Y
               ZJ0(NJ,NA)=Z
   50       CONTINUE
            JNUM0(NA)=NJMAX
         ELSEIF(KID.EQ.'X') THEN
            GOTO 1
         ELSE
            WRITE(6,*) 'XX UNKNOWN KID: ',KID
            GOTO 2
         ENDIF            
C
  100 CONTINUE
      CALL MODANT(IERR)
      IF(IERR.NE.0) GOTO 9000
C
 9000 RETURN
      END
