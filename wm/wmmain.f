C     $Id$
C
C*** /TASK/WM **********************************************************
C
C              ANALYSIS OF ICRF WAVE PROPAGATION AND ABSORPTION
C              ------------------------------------------------
C
C                                A. FUKUYAMA
C             Department of Nuclear Engineering, Kyoto University
C                      Sakyo-ku, Kyoto 606-8501, Japan
C                    Email: fukuyama@nucleng.kyoto-u.ac.jp
C                  URL: http://p-grp.nucleng.kyoto-u.ac.jp/wm/
C***********************************************************************
C
      INCLUDE 'wmcomm.h'
      CHARACTER KID*1,LINE*80
C
      CALL MPINIT(NPROCS,MYRANK)
      IF(NPROCS.LT.NCPUMIN) THEN
         WRITE(6,*) 'XX NPROCS.LT.NCPUMIN :',NPROCS,'.LT.',NCPUMIN
         STOP
      ENDIF
      IF(NPROCS.GT.NCPUMAX) THEN
         WRITE(6,*) 'XX NPROCS.GT.NCPUMAX :',NPROCS,'.GT.',NCPUMAX
         STOP
      ENDIF
C
      IF(MYRANK.EQ.0) THEN
         OPEN(33,STATUS='SCRATCH',FORM='FORMATTED')
         WRITE(6,*) '######## /TASK/WM V3.42 00/10/31 ########'
         CALL GSOPEN
      ENDIF
      CALL MPSYNC
C
      CALL WMINIT
      IF(MYRANK.EQ.0) CALL WMPARF
      CALL MPSYNC
      CALL WMPRBC
C
      MODE=0
    1 CONTINUE
         IF(MYRANK.EQ.0) THEN
            WRITE(6,*)
     &      '## INPUT: P,V/PARM  W/WAVE  A,F,C/AMP  E,S/SCAN  ',
     &      'G/GRAPH  T/TAE  Q/QUIT'
            CALL WMKLIN(LINE,KID,MODE)
         ENDIF
         CALL MPBCIA(MODE)
         IF(MODE.EQ.1) THEN
            CALL MPBCKA(KID)
            GOTO 2
         ENDIF
         IF(MODE.EQ.2) THEN
            CALL WMPRBC
            GOTO 1
         ENDIF
         IF(MODE.EQ.3) GOTO 1
C
    2    CONTINUE
         IF (KID.EQ.'P') THEN
            IF(MYRANK.EQ.0) CALL WMPARM(KID)
            CALL MPBCKA(KID)
            CALL WMPRBC
         ELSE IF(KID.EQ.'V') THEN
            IF(MYRANK.EQ.0) CALL WMVIEW
            CALL MPSYNC
            KID=' '
C
C        *** WAVE CALCULATION ***
C
         ELSEIF (KID.EQ.'W') THEN
            MODEEG=0
            CALL WMSETG(IERR)
               IF(IERR.NE.0) GOTO 1
            CALL WMSETJ(IERR)
               IF(IERR.NE.0) GOTO 1
C
            CALL WMSOLV
            CALL WMEFLD
            CALL WMBFLD
            CALL WMPABS
            IF(MYRANK.EQ.0) THEN
               CALL WMPFLX
               CALL WMPANT
               CALL WMPOUT
               IF(MODELW.EQ.1) CALL WMDOUT(IERR)
            ENDIF
            CALL MPSYNC
            KID=' '
         ELSEIF (KID.EQ.'A') THEN
            RF=DBLE(CRF)
            RFI=DIMAG(CRF)
            CALL DIAMIN(RF,RFI,AMPL)
            IF(MYRANK.EQ.0) 
     &      WRITE(6,'(A,1P3E12.4)') '      RF,RFI,AMPL=',RF,RFI,AMPL
            KID=' '
C
C        *** AMPLITUDE SURVEY ***
C
         ELSE IF (KID.EQ.'F') THEN
            CALL WMAM1D(KID)
         ELSE IF (KID.EQ.'C') THEN
            CALL WMAM2D(KID)
C
C        *** EIGENMODE ***
C
         ELSE IF (KID.EQ.'E') THEN
            CALL WMEIGN(KID)
         ELSE IF (KID.EQ.'S') THEN
            CALL WMSCAN(KID)
C
C        *** GRAPHICS ***
C
         ELSE IF (KID.EQ.'G') THEN
            IF(MYRANK.EQ.0) CALL WMGOUT
            CALL MPSYNC
            KID=' '
C
C        *** TAE FREQUENCY ***
C
         ELSE IF (KID.EQ.'T') THEN
            MODEEG=0
            CALL WMSETG(IERR)
               IF(IERR.NE.0) GOTO 1
            CALL WMSETJ(IERR)
               IF(IERR.NE.0) GOTO 1
            IF(MYRANK.EQ.0) CALL WMTAE
            CALL MPSYNC
            KID=' '
C
         ELSE IF(KID.EQ.'H') THEN
            WRITE(6,*) '# KID:  P: PARAMETER INPUT (VARNAME = VALUE)'
            WRITE(6,*) '        V: VIEW PARAMETERS'
            WRITE(6,*) '        W: WAVE EXCITED BY EXTERNAL ANTENNA'
            WRITE(6,*) '        A: AMPLITUDE OF INTERNALLY EXCITED WAVE'
            WRITE(6,*) '        F: REAL FREQUENCY SCAN OF AMPLITUDE'
            WRITE(6,*) '        C: COMPLEX FREQUENCY SCAN OF AMPLITUDE'
            WRITE(6,*) '        E: EIGEN VALUE SEARCH'
            WRITE(6,*) '        S: PARAMETER SCAN OF EIGEN VALUE'
            WRITE(6,*) '        G: GRAPHICS'
            WRITE(6,*) '        H: HELP'
            WRITE(6,*) '        Q: QUIT'
            KID=' '
C
         ELSE IF(KID.EQ.'Q') THEN
            GOTO 9000
C
         ELSE IF(KID.EQ.'X'.OR.KID.EQ.'#') THEN
            KID=' '
         ELSE
            IF(MYRANK.EQ.0) WRITE(6,*) 'XX WMMAIN: UNKNOWN KID'
            KID=' '
         END IF
C
      IF(KID.NE.' '.AND.KID.NE.'Q') GOTO 2
      GO TO 1
C
 9000 IF(MYRANK.EQ.0) CALL GSCLOS
      IF(MYRANK.EQ.0) CLOSE(33)
      CALL MPSYNC
C
      CALL MPTERM
C
      STOP
      END
C
C     ***** INPUT KID or LINE *****
C                   MODE=0: LINE INPUT 
C                        1: KID INPUT
C                        2: PARM INPUT
C                        3: NEW PROMPT
C
      SUBROUTINE WMKLIN(LINE,KID,MODE)
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
         CALL WMPARL(LINE)
         MODE=2
         RETURN
      ENDIF
C
      KID=LINE(1:1)
      CALL GUCPTL(KID)
      IF(KID.EQ.'P'.OR.
     &   KID.EQ.'V'.OR.
     &   KID.EQ.'W'.OR.
     &   KID.EQ.'A'.OR.
     &   KID.EQ.'F'.OR.
     &   KID.EQ.'C'.OR.
     &   KID.EQ.'E') THEN
         MODE=1
         RETURN
      ENDIF
      IF(KID.EQ.'S'.OR.
     &   KID.EQ.'G'.OR.
     &   KID.EQ.'T'.OR.
     &   KID.EQ.'H'.OR.
     &   KID.EQ.'X'.OR.
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
