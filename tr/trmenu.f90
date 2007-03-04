!     ***** TASK/TR MENU *****

      SUBROUTINE TRMENU

      USE TRCOMM, ONLY : KUFDCG, KUFDEV, MDLUF, MDLXP, NGR, NGT, NRMAX, NSMAX, NT, NTMAX, NTMAX_SAVE, &
     &                   ALLOCATE_TRCOMM
      USE TRCOM1, ONLY : KDIRX
      IMPLICIT NONE
      INTEGER(4)       :: IERR, MODE, NFL, NFLMAX, NTMOLD
      INTEGER(4), SAVE :: INIT=0
      CHARACTER(LEN=1) :: KID
      CHARACTER(LEN=80):: LINE
      EXTERNAL TRPARM

!     ------ SELECTION OF TASK TYPE ------

      IERR=0

    1 IF(INIT.EQ.0) THEN
         WRITE(6,601)
  601    FORMAT('## TR MENU: P,V,U/PARM  R/RUN  L/LOAD  ','D/DATA  H/HELP  Q/QUIT')
      ELSEIF(INIT.EQ.1) THEN
         WRITE(6,602)
  602    FORMAT('## TR MENU: G/GRAPH  '/ &
     &          '            P,V,U/PARM  R/RUN  L/LOAD  ','D/DATA  H/HELP  Q/QUIT')
      ELSE
         WRITE(6,603)
  603    FORMAT('## TR MENU: C/CONT  G/GRAPH  W/WRITE  ','S/SAVE  O/UFILEOUT  M/MDLTST'/ &
     &          '            P,V,U/PARM  R/RUN  L/LOAD  ','D/DATA  H/HELP  Q/QUIT')
      ENDIF

      CALL TASK_KLIN(LINE,KID,MODE,TRPARM)
      IF(MODE.NE.1) GOTO 1

      IF(KID.EQ.'P') THEN
         CALL TRPARM(0,'TR',IERR)
      ELSE IF(KID.EQ.'V') THEN
         CALL TRVIEW(0)
      ELSE IF(KID.EQ.'U') THEN
         CALL TRVIEW(1)

      ELSE IF(KID.EQ.'L') THEN
         CALL TRLOAD
         CALL PLDATA_SETN(NRMAX,NSMAX)
         CALL PLDATA_CLEAR
         INIT=2
      ELSE IF(KID.EQ.'S'.AND.INIT.EQ.2) THEN
         CALL TRSAVE

      ELSE IF(KID.EQ.'R') THEN
         CALL ALLOCATE_TRCOMM(IERR)
           IF(IERR.NE.0) GOTO 1
         IF(MDLUF.NE.0.AND.MDLXP.NE.0) CALL IPDB_OPEN(KUFDEV, KUFDCG)
         CALL PLDATA_SETN(NRMAX,NSMAX)
         CALL PLDATA_CLEAR
         IF(MDLUF.NE.0) CALL UFILE_INTERFACE(KDIRX,KUFDEV,KUFDCG,0)
         CALL TR_EQS_SELECT(0)
         IF(MDLUF.EQ.1) THEN
            IF(INIT.EQ.2.AND.NT.NE.0) THEN
               NT=0
               NTMAX=NTMAX_SAVE
            ENDIF
            CALL TR_UFILE_CONTROL(1)
         ELSEIF(MDLUF.EQ.2) THEN
            CALL TR_UFILE_CONTROL(2)
         ELSEIF(MDLUF.EQ.3) THEN
            IF(INIT.EQ.2.AND.NT.NE.0) THEN
               NT=0
               NTMAX=NTMAX_SAVE
            ENDIF
            CALL TR_UFILE_CONTROL(3)
         ELSE
            CALL TR_UFILE_CONTROL(0)
         ENDIF
         CALL TRPROF
         CALL TRLOOP

         INIT=2
         NTMOLD=NTMAX
!
      ELSE IF(KID.EQ.'C'.AND.INIT.EQ.2) THEN
         IF(MDLUF.EQ.1) THEN
            NT=NTMOLD
            NTMAX=NTMAX_SAVE+NTMOLD
         ELSE
            NT=0
         ENDIF
         CALL TRLOOP
         NTMOLD=NTMAX

      ELSE IF(KID.EQ.'G'.AND.INIT.GE.1) THEN
         CALL TRGOUT

      ELSE IF(KID.EQ.'W'.AND.INIT.EQ.2) THEN
!         write(6,*)  "J0=",AJ(1)*1.D-6
  102    WRITE(6,*) '# SELECT ',  ': PRINT TYPE (1..9)  N/NAMELIST H/HELP  X/EXIT'
         READ(5,'(A1)',ERR=102,END=1) KID
         CALL GUCPTL(KID)
         IF(KID.EQ.'H') THEN
            CALL TRHELP('W')
         ELSEIF(KID.EQ.'X') THEN
            GOTO 1
         ELSE
            CALL TRPRNT(KID)
         ENDIF
         GOTO 102

      ELSE IF(KID.EQ.'D') THEN
! 200607 start delete
!         NFLMAX=0
!    4    WRITE(6,*) '## HOW MANY DATA FILES ?'
!         READ(5,*,END=1,ERR=4) NFLMAX
!         DO NFL=1,NFLMAX
!            CALL TRLOAD
!            NGR=NFL-1
!            NGT=NFL-1
!            CALL TRCALC(IERR)
!            CALL TRGLOB
!            CALL TRATOT
!            CALL TRATOG
!         ENDDO
!         INIT=1
! 200607 end delete

      ELSE IF(KID.EQ.'O'.AND.INIT.EQ.2) THEN
         CALL TRXOUT
      ELSE IF(KID.EQ.'M'.AND.INIT.EQ.2) THEN
         CALL TRMDLT
      ELSE IF(KID.EQ.'H') THEN
         CALL TRHELP('M')
      ELSE IF(KID.EQ.'Q') THEN
         GOTO 9000

      ELSE IF(KID.EQ.'X'.OR.KID.EQ.'#') THEN
         CONTINUE
      ELSE
         WRITE(6,*) 'XX TRMAIN: UNKNOWN KID'
      END IF

      GOTO 1

 9000 IF(MDLUF.NE.0.AND.MDLXP.NE.0) CALL IPDB_CLOSE
      RETURN
      END SUBROUTINE TRMENU
