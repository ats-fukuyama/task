! trmenu.f90

MODULE trmenu

  PRIVATE
  PUBLIC tr_menu

CONTAINS

!     ***** TASK/TR MENU *****

  SUBROUTINE tr_menu

    USE trcomm
    USE trparm
    USE trview
    USE trprep
    USE trloop
    USE trfile
    USE trgout
    USE trfout
    USE trpnf
      
    USE libfio
    USE libkio
    USE libchar
    
    IMPLICIT NONE
    INTEGER, SAVE:: INIT=0
    CHARACTER(LEN=1) :: KID
    CHARACTER(LEN=80):: LINE
    INTEGER:: IERR,MODE,NTMOLD
    INTEGER:: NR,NS,NF,NTYPE,id_loop

!     ------ SELECTION OF TASK TYPE ------

    IERR=0

1   IF(INIT.EQ.0) THEN
       WRITE(6,601)
601    FORMAT('## TR MENU: P,V,U/PARM  R/RUN  L/LOAD  ', &
              'D/DATA  H/HELP  Q/QUIT')
    ELSEIF(INIT.EQ.1) THEN
       WRITE(6,602)
602    FORMAT('## TR MENU: G/GRAPH  ', &
              '            P,V,U/PARM  R/RUN  L/LOAD  '/ &
              'D/DATA  H/HELP  Q/QUIT')
    ELSE
       WRITE(6,603)
603    FORMAT('## TR MENU: C/CONT  G/GRAPH  W,F/WRITE ', &
              'S/SAVE  P,V,U/PARM  R/RUN  L/LOAD  '/ &
              'D/DATA  H/HELP  Q/QUIT')
    ENDIF

    CALL TASK_KLIN(LINE,KID,MODE,tr_parm)
    IF(MODE.NE.1) GOTO 1

    IF(KID.EQ.'P') THEN
       CALL tr_parm(0,'TR',IERR)

    ELSE IF(KID.EQ.'V') THEN
       CALL tr_view(0)

    ELSE IF(KID.EQ.'U') THEN
       CALL tr_view(1)

    ELSE IF(KID.EQ.'L') THEN
       CALL tr_prep(ierr)
       if(ierr.ne.0) GO TO 1
       CALL tr_load
       INIT=2
       id_loop=0
       NT=0

    ELSE IF(KID.EQ.'S'.AND.INIT.EQ.2) THEN
       CALL tr_save

    ELSE IF(KID.EQ.'R') THEN
       id_loop=0
       CALL tr_prep(ierr)
       if(ierr.ne.0) GO TO 1

       CALL tr_loop(ierr)
       IF(ierr.NE.0) id_loop=1
         
       INIT=2
       NTMOLD=NTMAX

    ELSE IF(KID.EQ.'C'.AND.INIT.EQ.2) THEN
       NT=0
       IF(id_loop.EQ.0) THEN
          CALL tr_loop(ierr)
          IF(ierr.NE.0) id_loop=1
          NTMOLD=NTMAX
       END IF
       
    ELSE IF(KID.EQ.'G'.AND.INIT.GE.1) THEN
       CALL tr_gout

    ELSE IF(KID.EQ.'F'.AND.INIT.GE.1) THEN
       CALL tr_fout

    ELSE IF(KID.EQ.'W'.AND.INIT.EQ.2) THEN
102    CONTINUE
       WRITE(6,*) &
            '# SELECT ',  ': PRINT TYPE (1..9)  N/NAMELIST H/HELP  X/EXIT'
       READ(5,'(A1)',ERR=102,END=1) KID
       CALL toupper(KID)
       IF(KID.EQ.'H') THEN
          CALL TRHELP('W')
       ELSEIF(KID.EQ.'X') THEN
          GOTO 1
       ELSE
          CALL TRPRNT(KID)
       ENDIF
       GOTO 102

    ELSE IF(KID.EQ.'D') THEN
4      CONTINUE
       WRITE(6,*) '## profile data output: 1:type, 0:end'
       READ(5,*,ERR=4,END=1) NTYPE
       IF(NTYPE.EQ.0) GO TO 1
       SELECT CASE(NTYPE)
       CASE(1)
          CALL FWOPEN(26,'trdata1',0,1,'trdata',IERR)
       CASE(2)
          CALL FWOPEN(26,'trdata2',0,1,'trdata',IERR)
       CASE(3)
          CALL FWOPEN(26,'trdata3',0,1,'trdata',IERR)
       CASE(4)
          CALL FWOPEN(26,'trdata4',0,1,'trdata',IERR)
       CASE(5)
          CALL FWOPEN(26,'trdata5',0,1,'trdata',IERR)
       CASE(6)
          CALL FWOPEN(26,'trdata6',0,1,'trdata',IERR)
       CASE DEFAULT
          WRITE(6,*) 'XX unknown ntype'
          GO TO 4
       END SELECT
       IF(IERR.NE.0) GO TO 4
       WRITE(26) NRMAX,NSMAX,NFMAX
       WRITE(26) (RM(NR),RG(NR),NR=1,NRMAX)
       WRITE(26) ((RN(NR,NS),RT(NR,NS),NR=1,NRMAX),NS=1,NSMAX)
       WRITE(26) ((RW(NR,NF),RNF(NR,NF),RTF(NR,NF),NR=1,NRMAX),NF=1,NFMAX)
       CLOSE(26)
       WRITE(6,'(A,I1)') '## Data saved in trdata',NTYPE
       GO TO 4
            
    ELSE IF(KID.EQ.'B') THEN
7777   CONTINUE
       WRITE(6,*) '## INPUT T'
       READ(5,*,END=7779) T
       WRITE(6,*) '   sigma=',SIGMAM(T,T)
       GOTO 7777
7779   CONTINUE
       
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

9000 CONTINUE
    RETURN
  END SUBROUTINE tr_menu
END MODULE trmenu
