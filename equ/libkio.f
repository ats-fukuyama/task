C     $Id$
C
      module libkio_mod
      use libchar_mod
      public
      contains
C
C     ***** INPUT KID or LINE *****
C                   MODE=0: LINE INPUT 
C                        1: KID INPUT
C                        2: PARM INPUT
C                        3: NEW PROMPT
C
      SUBROUTINE TASK_KLIN(LINE,KID,MODE,XXPARM)
C
      USE libchar
      EXTERNAL XXPARM
      CHARACTER LINE*80,KID*1
C
      READ(5,'(A80)',ERR=2,END=3) LINE
C
      ID=0
      DO I=1,80
         IF(LINE(I:I).EQ.'=') ID=1
      ENDDO
      IF(ID.EQ.1) THEN
         CALL XXPARM(2,LINE,IERR)
         MODE=2
         RETURN
      ENDIF
C
      KID=LINE(1:1)
      CALL toupper(KID)
      IF((KID.GE.'A'.AND.KID.LE.'Z').OR.
     &    KID.EQ.'?'.OR.KID.EQ.'#') THEN
         MODE=1
         RETURN
      ENDIF
C
      KID=' '
      MODE=0
      RETURN
C
    2 WRITE(6,'(A)') 'XX INPUT ERROR !'
      MODE=3
      RETURN
C
    3 KID='Q'
      MODE=1
      RETURN
      END SUBROUTINE TASK_KLIN
C
C     ****** INPUT PARAMETERS ******
C
      SUBROUTINE TASK_PARM(MODE,KWD,KIN,XXNLIN,XXPLST,IERR)
C
C     MODE=0 : standard namelist input  KIN: not used
C     MODE=1 : namelist file input      KIN: file name
C     MODE=2 : namelist line input      KIN: data string
C
C     KWD:     namelist name (not longer than 6 chars)
C     XXNLIN:  namelist subroutine XXNLIN(NID,IST,IERR)
C                      input  NID:  device number
C                      output IST:  io status
C                             IERR: error indicater
C     XXPLST:  error prompt subroutine XXPLST
C                      indicating valid namelist variables
C
C     IERR=0 : normal end
C     IERR=1 : namelist standard input end of file
C     IERR=2 : namelist file does not exist
C     IERR=3 : namelist file open error
C     IERR=4 : namelist file read error
C     IERR=5 : namelist file abormal end of file
C     IERR=6 : namelist line input error
C     IERR=7 : unknown MODE
C
      EXTERNAL XXNLIN,XXPLST
      LOGICAL LEX
      CHARACTER KWD*(*),KIN*(*),LINE*80,KNLINE*94,KNL*6
C
      IERR=0
      KNL=KWD
      LINE=KIN
C
      IF(MODE.EQ.0) THEN
    1    CONTINUE
         CALL KTRIM(KNL,KL)
         WRITE(6,'(A,A,A)') '## INPUT ',KNL(1:KL),' :'
         CALL XXNLIN(5,IST,IERR)
         IF(IERR.EQ.8) THEN
            CALL XXPLST
            GOTO 1
         ENDIF
         IF(IERR.EQ.9) IERR=1
C
      ELSEIF(MODE.EQ.1) THEN
C
         INQUIRE(FILE=LINE,EXIST=LEX,ERR=9800)
         IF(.NOT.LEX) THEN
            IERR=2
            RETURN
         ENDIF
         OPEN(25,FILE=LINE,IOSTAT=IST,STATUS='OLD',ERR=9100)
         CALL XXNLIN(25,IST,IERR)
         IF(IERR.EQ.8) GOTO 9800
         IF(IERR.EQ.9) GOTO 9900
         CLOSE(25)
         CALL KTRIM(LINE,KL)
         WRITE(6,'(A,A,A)') 
     &        '## FILE (',LINE(1:KL),') IS ASSIGNED FOR PARM INPUT'
C
      ELSEIF(MODE.EQ.2) THEN
         CALL KTRIM(KNL,KL1)
         CALL KTRIM(LINE,KL2)
         KNLINE=' &'//KNL(1:KL1)//' '//LINE(1:KL2)//' &END'
         WRITE(7,'(A90)') KNLINE
         REWIND(7)
         CALL XXNLIN(7,IST,IERR)
         REWIND(7)
         IF(IERR.EQ.8.OR.IERR.EQ.9) THEN
            WRITE(6,'(A)') '## PARM INPUT ERROR.'
            CALL XXPLST
            IERR=6
            RETURN
         ENDIF
         WRITE(6,'(A)') '## PARM INPUT ACCEPTED.'
      ELSE
         WRITE(6,'(A,I4)') 'XX XXPARM : UNKNOWN MODE =',MODE
         IERR=7
         RETURN
      ENDIF
      RETURN
C
 9100 WRITE(6,'(A,I6)') 'XX PARM FILE OPEN ERROR : IOSTAT = ',IST
      IERR=3
      RETURN
 9800 WRITE(6,'(A,I6)') 'XX PARM FILE READ ERROR : IOSTAT = ',IST
      IERR=4
      RETURN
 9900 WRITE(6,'(A)') 'XX PARM FILE EOF ERROR'
      IERR=5
      RETURN
      END SUBROUTINE TASK_PARM
c
      end module libkio_mod
