C     $Id: wmmenu.f,v 1.15 2013/01/20 23:24:03 fukuyama Exp $
C
C     ***** TASK/WM MENU *****
C
      SUBROUTINE WMMENU
C
C      use plfile_prof_mod
      use wmtest
      INCLUDE 'wmcomm.inc'
C
      EXTERNAL WMPARM
      CHARACTER KID*1,LINE*80
      SAVE INIT,NFILEINI
      DATA INIT/0/
C
      IF(INIT.EQ.0) THEN
         NFILEINI=0
         INIT=1
      ENDIF
C
    1 CONTINUE
         IF(NRANK.EQ.0) THEN
            WRITE(6,601)
  601       FORMAT('## WM MENU: P,V/PARM R/RUN D0-3/AMP F/ROOT ',
     &      'G,Y/GRAPH T/TAE O/OUT S,W/SAVE L/LOAD Q/QUIT')
            CALL TASK_KLIN(LINE,KID,MODE,WMPARM)
         ENDIF
         CALL mtx_broadcast1_integer(MODE)
         IF(MODE.EQ.2) CALL WMPRBC
      IF(MODE.NE.1) GOTO 1
C
    2 CONTINUE
         CALL mtx_broadcast1_character(KID)
C
         IF (KID.EQ.'P') THEN
            IF(NRANK.EQ.0) CALL WMPARM(0,'WM',IERR)
            CALL mtx_barrier
            CALL WMPRBC
            KID=' '
         ELSE IF(KID.EQ.'V') THEN
            IF(NRANK.EQ.0) CALL WMVIEW
            CALL mtx_barrier
            KID=' '
C
C        *** WAVE CALCULATION ***
C
         ELSEIF (KID.EQ.'R') THEN
C            CALL plfile_prof_read(modeln,modelq,ierr)
            
            NPH0_SV  = NPH0
            NPHMAX_SV = NPHMAX
            NHHMAX_SV = NHHMAX

            IF(NHHMAX.GT.1) THEN
              NHHMAX=NPHMAX/NHC
              NPHMAX=NHC
            ELSE
              NHC=1
            END IF
            IF (NPHMAX .EQ. 1 .and. NHHMAX .NE. 1)THEN
               CALL WMEXEC(IERR)
            ELSE
               CALL WM_LOOP(IERR)
            ENDIF
            CALL mtx_barrier
            IF(IERR.NE.0) GOTO 1
            KID=' '
C
C        *** AMPLITUDE SURVEY ***
C
      ELSEIF(KID.EQ.'D') THEN
C         CALL plfile_prof_read(modeln,modelq,ierr)
         READ(LINE(2:),*,ERR=1,END=1) NID
         IF(NID.EQ.0) THEN
            CALL WMAM0D(KID,LINE)
         ELSEIF(NID.EQ.1) THEN
            CALL WMAM1D(KID,LINE)
         ELSEIF(NID.EQ.2) THEN
            CALL WMAM2D(KID,LINE)
         ELSEIF(NID.EQ.3) THEN
            CALL WMSCAN(KID,LINE)
         ELSE
            WRITE(6,*) 'XX WMMENU: unknown NID'
         ENDIF
C
C        *** FIND ROOT ***
C
         ELSE IF (KID.EQ.'F') THEN
C            CALL plfile_prof_read(modeln,modelq,ierr)
            CALL WMEIGN(KID,LINE)
C
C        *** GRAPHICS ***
C
         ELSE IF (KID.EQ.'G') THEN
            IF(NRANK.EQ.0) CALL WMGOUT
            CALL mtx_barrier
            KID=' '
C
C        *** FILE OUTPUT ***
C
         ELSE IF (KID.EQ.'S') THEN
            IF(NRANK.EQ.0) THEN
               CALL WMSAVE
            ENDIF
            CALL mtx_barrier
            KID=' '
         ELSE IF (KID.EQ.'L') THEN
            NPH0_SV  = NPH0
            NPHMAX_SV = NPHMAX
            NHHMAX_SV = NHHMAX

            IF(NHHMAX.GT.1) THEN
              NHHMAX=NPHMAX/NHC
              NPHMAX=NHC
            ELSE
              NHC=1
            END IF

            IF(NRANK.EQ.0) THEN
               CALL WMLOAD
            ENDIF
            CALL mtx_barrier
            KID=' '
         ELSE IF (KID.EQ.'W') THEN
            IF(NRANK.EQ.0) THEN
               IF(NFILEINI.EQ.0) THEN
                  REWIND(26)
                  NFILEINI=1
               ENDIF
               CALL WMWOUT
            ENDIF
            CALL mtx_barrier
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
            IF(NRANK.EQ.0) CALL WMTAE
            CALL mtx_barrier
            KID=' '
C
C        *** Pabs(r,s) output for TOPICS ***
C
         ELSE IF (KID.EQ.'O') THEN
            CALL WM_TOPICS(IERR)
            IF(IERR.NE.0) GOTO 1
            CALL mtx_barrier
            KID=' '
C
         ELSE IF(KID.EQ.'H') THEN
            WRITE(6,*) '# KID: P: PARAMETER INPUT (VARNAME = VALUE)'
            WRITE(6,*) '       V:  VIEW PARAMETERS'
            WRITE(6,*) '       R:  WAVE EXCITED BY EXTERNAL ANTENNA'
            WRITE(6,*) '       D0: AMPLITUDE OF INTERNALLY EXCITED WAVE'
            WRITE(6,*) '       D1: FREQUENCY SCAN OF AMPLITUDE'
            WRITE(6,*) '       D2: COMPLEX FREQUENCY SCAN OF AMPLITUDE'
            WRITE(6,*) '       D3: PARAMETER SCAN OF EIGEN VALUE'
            WRITE(6,*) '       F:  FIND EIGEN VALUE'
            WRITE(6,*) '       G: GRAPHICS'
            WRITE(6,*) '       W: FILE 26 OUTPUT: FIELD DATA'
            WRITE(6,*) '       H: HELP'
            WRITE(6,*) '       Q: QUIT'
            KID=' '
C
         ELSEIF (KID.EQ.'Z') THEN
            CALL WMDEBUG(IERR)
            CALL mtx_barrier
            IF(IERR.NE.0) GOTO 1
            KID=' '
C
         ELSEIF (KID.EQ.'Y') THEN
            CALL WM_TEST(CEFLD3D,CPABS3D,PABST3D,PABSTT3D,
     &                   MDM,NPHM,NRM,NSM,
     &                   NTHMAX,NPHMAX,NRMAX,NSMAX,
     &                   RGMIN,RGMAX,RAXIS,NCONT,NGRAPH)
            KID=' '
C
         ELSE IF(KID.EQ.'X') THEN
            KID=' '
C
         ELSE IF(KID.EQ.'Q') THEN
            GOTO 9000
         ELSE
            IF(NRANK.EQ.0) WRITE(6,*) 'XX WMMENU: UNKNOWN KID'
            KID=' '
         END IF
C
      IF(KID.NE.' '.AND.KID.NE.'Q') GOTO 2
      GO TO 1
C
 9000 RETURN
      END
