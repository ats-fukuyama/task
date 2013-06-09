C     $Id$
C
C             ########## /TASK/WF ##########
C
C             WAVE PROPAGATION AND ABSORPTION
C
C                  HOT/COLD PLASMA
C                  3-DIMENSIONAL INHOMOGENEITY
C                  FIRST ORDER FINITE ELEMENT METHOD
C                  SCALAR AND VECTOR POTENTIAL
C
C                  CODED BY A. FUKUYAMA
C
C                  DEPARTMENT OF NUCLEAR ENGINEERING
C                  KYOTO UNIVERSITY
C                  KYOTO 606-8501, JAPAN
C
C                  V1.0   : 1996 JULY (CARTESIAN COORDINATES)
C                  V1.01  : 1996 SEPT (CYLINDRICAL COORDINATES)
C                  V1.1   : 1996 SEPT (TRANSLATIONAL THREE-DIMENSIONAL)
C                  V2.0   : 2003 JAN  (FULL THREE-DIMENSIONAL)
C                  V3.0   : 2003 JUN  (SIDE ELEMENTS)
C
C     ************************************************************
C
      INCLUDE 'wfcomm.inc'
C
      CHARACTER KID*1,LINE*80
C
      WRITE(6,*) '## /TASK/WFX  V3.01  03/07/28 ###'
      CALL GSOPEN
      CALL SETAIF
      CALL WFINIT
      CALL WFPARF
C
    1 CONTINUE
         WRITE(6,601)
  601    FORMAT(' ','## INPUT: P,V:PARM  D:DIV  A:ANT  W,C:WAVE',
     &                   '  G:GRAPH  N:NAS  S,L:FILE  Q:QUIT')
         CALL GUFLSH
         CALL WFKLIN(LINE,KID,MODE)
         IF(MODE.NE.1) GOTO 1
C
         IF (KID.EQ.'P') THEN
            CALL WFPARM(KID)
            GOTO 1
         ELSEIF (KID.EQ.'V') THEN
            CALL WFVIEW
         ELSEIF (KID.EQ.'D') THEN
            CALL WFDIV
         ELSEIF (KID.EQ.'A') THEN
            IF(NNMAX.EQ.0) CALL WFRELM(ID)
            CALL WFANT
         ELSEIF (KID.EQ.'C') THEN
            CALL WFWPRE(IERR)
         ELSEIF (KID.EQ.'W') THEN
            CALL WFWAVE
         ELSEIF (KID.EQ.'G') THEN
            CALL WFGOUT
         ELSEIF (KID.EQ.'N') THEN
            CALL WFNAS
         ELSEIF (KID.EQ.'S') THEN
            CALL WFWFLD
         ELSEIF (KID.EQ.'L') THEN
            CALL WFRFLD
         ELSEIF (KID.EQ.'?') THEN
            CALL WFINFO
         ELSEIF (KID.EQ.'Q') THEN
            GO TO 9000
         END IF
         KID=' '
      GO TO 1
C
 9000 CALL GSCLOS
      IF(NFOPEN.NE.0) CLOSE(26)
      STOP
      END
C
C     ***** INPUT KID or LINE *****
C                   MODE=0: LINE INPUT 
C                        1: KID INPUT
C                        2: PARM INPUT
C                        3: NEW PROMPT
C
      SUBROUTINE WFKLIN(LINE,KID,MODE)
C
      CHARACTER LINE*80,KID*1
C
      READ(5,'(A80)',ERR=2,END=3) LINE
C
      ID=0
      DO I=1,80
         IF(LINE(I:I).EQ.'=') ID=1
      ENDDO
      IF(ID.EQ.1) THEN
         CALL WFPARL(LINE)
         MODE=2
         RETURN
      ENDIF
C
      KID=LINE(1:1)
      CALL GUCPTL(KID)
      IF(KID.EQ.'P'.OR.
     &   KID.EQ.'V'.OR.
     &   KID.EQ.'D'.OR.
     &   KID.EQ.'A'.OR.
     &   KID.EQ.'N'.OR.
     &   KID.EQ.'S'.OR.
     &   KID.EQ.'L'.OR.
     &   KID.EQ.'W'.OR.
     &   KID.EQ.'G'.OR.
     &   KID.EQ.'?'.OR.
     &   KID.EQ.'Q') THEN
         MODE=1
         RETURN
      ENDIF
C
      KID=' '
      MODE=0
      RETURN
C
    2 WRITE(6,*) 'XX INPUT ERROR !'
      MODE=3
      RETURN
C
    3 KID='Q'
      MODE=1
      RETURN
      END
C
C     ***** DEBUG INFORMATION ROUTINE *****
C
      SUBROUTINE WFINFO
C
      INCLUDE 'wfcomm.inc'
      CHARACTER KID*1
C
 8001 CONTINUE
         WRITE(6,*) '## INPUT: E:element  N:node  B:boundary  '//
     &              'A,M:ant  F:FEP  D:EFINDK  X:end'
         READ(5,'(A1)',ERR=8001,END=9000) KID
         CALL GUCPTL(KID)
C
         IF(KID.EQ.'E') THEN
 8002       WRITE(6,*) '## INPUT: Element number '
            READ(5,*,ERR=8002,END=8001) IE
            IF(IE.EQ.0) GOTO 8001
            WRITE(6,'(A,5I8)') '  NE,KN=',IE,
     &               KNELM(1,IE),KNELM(2,IE),KNELM(3,IE),KNELM(4,IE)
            WRITE(6,'(A,5I8)') '  KA,ND=',KAELM(IE),
     &               NDELM(1,IE),NDELM(2,IE),NDELM(3,IE),NDELM(4,IE)
            DO I=1,4
               IN=NDELM(I,IE)
               RND=SQRT(XND(IN)**2+YND(IN)**2)
               WRITE(6,'(A,I5,1P4E12.4,I3,I5)') 'IN,X,Y,Z,R,KA,IM =',
     &           IN,XND(IN),YND(IN),ZND(IN),RND,KANOD(IN),IMLEN(IN)
            ENDDO
            GOTO 8002
C
         ELSEIF(KID.EQ.'N') THEN
 8003       WRITE(6,*) '## INPUT: Node number '
            READ(5,*,ERR=8003,END=8001) NN
            IF(NN.EQ.0) GOTO 8001
            WRITE(6,'(A,1P3E12.4,2I5)') '   XND,YND,ZND,KA,IM =',
     &           XND(NN),YND(NN),ZND(NN),KANOD(NN),IMLEN(NN)
            DO NE=1,NEMAX
               DO IN=1,4
                  IF(NDELM(IN,NE).EQ.NN) 
     &                 WRITE(6,'(A,3I8)') '   NE,IN,KN=',
     &                                        NE,IN,KNELM(IN,NE)
               ENDDO
            ENDDO
            DO NSF=1,NSFMAX
               DO IN=1,3
                  IF(NDSRF(IN,NSF).EQ.NN) 
     &                 WRITE(6,'(A,6I8)') '   NSF,IN,KN,ND=',
     &                                        NSF,IN,KNSRF(IN,NSF),
     &                       NDSRF(1,NSF),NDSRF(2,NSF),NDSRF(3,NSF)
               ENDDO
            ENDDO
            GOTO 8003
C
         ELSEIF(KID.EQ.'B') THEN
 8004       WRITE(6,*) '## INPUT: Boundary number '
            READ(5,*,ERR=8004,END=8001) IB
            IF(IB.EQ.0) GOTO 8001
            WRITE(6,*) '   NDBDY =',NDBDY(IB)
            GOTO 8004
C
         ELSEIF(KID.EQ.'A') THEN
 8005       WRITE(6,*) '## INPUT: Antenna number, Segment number '
            READ(5,*,ERR=8005,END=8001) IA,IS
            IF(IA.EQ.0) GOTO 8001
            WRITE(6,*) '   XJ0,YJ0,ZJ0 =',
     &                     XJ0(IS,IA),YJ0(IS,IA),ZJ0(IS,IA)
            GOTO 8005
C
         ELSEIF(KID.EQ.'M') THEN
 8006       WRITE(6,*) '## INPUT: Antenna number, Segment number '
            READ(5,*,ERR=8006,END=8001) IA,IS
            IF(IA.EQ.0) GOTO 8001
            WRITE(6,*) '   JELMT,XJ,YJ,ZJ =',
     &                     JELMT(IS,IA),XJ(IS,IA),YJ(IS,IA),ZJ(IS,IA)
            GOTO 8006
C
         ELSEIF(KID.EQ.'F') THEN
 8007       WRITE(6,*) '## INPUT: IE,X,Y,Z'
            READ(5,*,ERR=8007,END=8001) IE,X,Y,Z
            IF(IE.LT.0) GOTO 8001
            IF(IE.GT.NEMAX) IE=NEMAX
            CALL FEP(X,Y,Z,IE)
            WRITE(6,*) '   IE =',IE
            IF(IE.NE.0) THEN
               WRITE(6,*) '  NDELM=',
     &                NDELM(1,IE),NDELM(2,IE),NDELM(3,IE),NDELM(4,IE)
            ENDIF
            GOTO 8007
C
         ELSEIF(KID.EQ.'D') THEN
 8008       WRITE(6,*) '## INPUT: IE,IN1,IN2,IN3'
            READ(5,*,ERR=8008,END=8001) IE,IN1,IN2,IN3
            IF(IE.EQ.0) GOTO 8001
            IDEBUGS=IDEBUG
            IDEBUG=1
            CALL EFINDK(IE,IN1,IN2,IN3,L)
            IDEBUG=IDEBUGS
            WRITE(6,*) 'IE,L = ',IE,L
            GOTO 8008
C
         ELSEIF(KID.EQ.'V') THEN
            DO NN=1,NNMAX
               IF(KANOD(NN).GT.0) THEN
                  WRITE(6,'(I5,1P3E12.4)')
     &                 NN,XND(NN),YND(NN),ZND(NN)
               ENDIF
            ENDDO
C
         ELSEIF(KID.EQ.'X') THEN
            GOTO 9000
         ENDIF
      GOTO 8001
C
 9000 RETURN
      END
