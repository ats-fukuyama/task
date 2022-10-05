C     $Id$
C
C             ########## /TASK/TF ##########
C
C             WAVE PROPAGATION AND ABSORPTION
C
C                  HOT/COLD PLASMA
C                  2-DIMENSIONAL INHOMOGENEITY
C                  FIRST ORDER FINITE ELEMENT METHOD
C                  SCALAR AND VECTOR POTENTIAL
C
C                  CODED BY A. FUKUYAMA
C
C                  Ddepartment of Nuclear Engineering
C                  Kyoto Univeristy
C                  Kyoto 606-8501, Japan
C
C                  V1.0   : 1996 JUL (CARTESIAN COORDINATES)
C                  V1.01  : 1996 SEP (CYLINDRICAL COORDINATES)
C                  V1.40  : 1997 JUL (ZONING ADDED)
C                  V1.50  : 1997 AUG (ZONING ADDED)
C
C     ************************************************************
C
      USE wfhout
      INCLUDE 'wfcomm.inc'
C
      CHARACTER KID*1,LINE*80
C
      WRITE(6,*) '## /TASK/WF  V1.73  2003/08/28 ##'
      CALL GSOPEN
      CALL SETAIF
      CALL WFINIT
      CALL WFPARF
      OPEN (7,STATUS='SCRATCH')
      NEVOL=0
      CALL WFRELM
      IF(NNOD.NE.0) THEN
         CALL WFRZON
         CALL WFRANT
      ENDIF
C
    1 WRITE(6,601)
  601 FORMAT(1H ,'## INPUT: P,V:PARM  D:DIV  Z:ZONE  A:ANT  W:WAVE',
     &                   '  T:EVOL  G,H:GRAPH  '/
     &       1H ,'          R:RUN  C:CONT  S:SAVE  L:LOAD',
     &                   '  B:FREQ  Q:QUIT')
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
      ELSEIF (KID.EQ.'Z') THEN
         IF(NNOD.EQ.0) THEN
            CALL WFRELM
         ENDIF
         CALL WFZONE
      ELSEIF (KID.EQ.'A') THEN
         IF(NNOD.EQ.0) THEN
            CALL WFRELM
            CALL WFRZON
         ENDIF
         CALL WFANT
      ELSEIF (KID.EQ.'W') THEN
         NEVOL=0
         T=0.D0
         CALL WFSTUP
         CALL WFWAVE
      ELSEIF(KID.EQ.'T') THEN 
         NEVOL=1
         T=0.D0
         CALL WFSTUP
         CALL WFEVIN
         CALL WFEVST
         DO NT=1,NTMAX
            T=T+DT
            CALL WFEVOL
            CALL WFEVST
         ENDDO
      ELSEIF (KID.EQ.'B') THEN
         CALL WFFREQ
      ELSEIF (KID.EQ.'G') THEN
         CALL WFGOUT
      ELSEIF (KID.EQ.'H') THEN
         CALL wf_hout(xd,yd,nnod,nsmax,nzmax)
      ELSEIF(KID.EQ.'R') THEN
         NEVOL=2
         T=0.D0
         CALL WFSTUP
         CALL WFEVIN
         CALL WFEVST
         DO NT=1,NTMAX
            IF(MOD(NT-1,NTSTEP).EQ.0) CALL WFWAVE
            T=T+DT
            CALL WFEVOL
            CALL WFEVST
         ENDDO
      ELSEIF(KID.EQ.'C') THEN
         IF(NEVOL.EQ.0) THEN
            WRITE(6,*) 'XX BEFORE CONTINUE, T, F, OR R'
         ELSE IF(NEVOL.EQ.3) THEN
            DO NT=1,NTMAX
               T=T+DT
C               CALL TFEVOL
C               CALL TFEVST
            ENDDO
         ELSE
            DO NT=1,NTMAX
               IF(NEVOL.GE.2.AND.
     &            MOD(NT-1,NTSTEP).EQ.0) CALL WFWAVE
               T=T+DT
               CALL WFEVOL
               CALL WFEVST
            ENDDO
         ENDIF
      ELSEIF (KID.EQ.'S') THEN
         CALL WFWFLD
      ELSEIF (KID.EQ.'L') THEN
         CALL WFRFLD
      ELSEIF (KID.EQ.'Q') THEN
         GO TO 9000
      END IF
      KID=' '
      GO TO 1
C
 9000 CALL GSCLOS
      IF(NFOPEN.NE.0) CLOSE(26)
      CLOSE(7)
      STOP
      END
C
C     ***** INPUT KID or LINE *****
C
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
C     --- when the line includes '=' ---
C
      ID=0
      DO I=1,80
         IF(LINE(I:I).EQ.'=') ID=1
      ENDDO
C
      IF(ID.EQ.1) THEN
         CALL WFPARL(LINE)
         MODE=2
         RETURN
      ENDIF
C
C     --- when the first character is a valid command char ---
C
      KID=LINE(1:1)
      CALL WFCPTL(KID)
      IF(KID.EQ.'P'.OR.
     &   KID.EQ.'D'.OR.
     &   KID.EQ.'V'.OR.
     &   KID.EQ.'Z'.OR.
     &   KID.EQ.'A'.OR.
     &   KID.EQ.'W'.OR.
     &   KID.EQ.'F'.OR.
     &   KID.EQ.'T'.OR.
     &   KID.EQ.'B'.OR.
     &   KID.EQ.'G'.OR.
     &   KID.EQ.'H'.OR.
     &   KID.EQ.'R'.OR.
     &   KID.EQ.'C'.OR.
     &   KID.EQ.'M'.OR.
     &   KID.EQ.'S'.OR.
     &   KID.EQ.'L'.OR.
     &   KID.EQ.'Q') THEN
         MODE=1
         RETURN
      ENDIF
C
      KID=' '
      MODE=0
      RETURN
C
C     --- when input error occurs ---
C
    2 WRITE(6,*) 'XX INPUT ERROR !'
      MODE=3
      RETURN
C
C     --- when EOF detected ---
C
    3 KID='Q'
      MODE=1
      RETURN
      END
C
C     ***** CAPITALIZE *****
C
      SUBROUTINE WFCPTL(KID)
C
      CHARACTER KID*1
C
      ID=ICHAR(KID)
      IF(ID.GE.97.AND.ID.LE.122) THEN
         KID=CHAR(ID-32)
      ENDIF
      RETURN
      END
