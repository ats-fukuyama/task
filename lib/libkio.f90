!     ***** INPUT KID or LINE *****
!                   MODE=0: LINE INPUT
!                        1: KID INPUT
!                        2: PARM INPUT
!                        3: NEW PROMPT

      SUBROUTINE TASK_KLIN(LINE,KID,MODE,XXPARM)

      IMPLICIT NONE
      INTEGER(4), INTENT(OUT)  :: MODE
      CHARACTER(LEN=80),INTENT(OUT) :: LINE
      CHARACTER(LEN=1), INTENT(OUT) :: KID
      INTEGER(4)  :: I, ID, IERR, IKID
      EXTERNAL XXPARM

      READ(5,'(A80)',ERR=2,END=3) LINE

      ID=0
      DO I=1,80
         IF(LINE(I:I).EQ.'=') ID=1
      ENDDO
      IF(ID.EQ.1) THEN
         CALL XXPARM(2,LINE,IERR)
         MODE=2
         RETURN
      ENDIF

      KID=LINE(1:1)
      IKID=ICHAR(KID)
      IF(IKID.GE.97.AND.IKID.LE.122) IKID=IKID-32
      KID=CHAR(IKID)

      IF((KID.GE.'A'.AND.KID.LE.'Z').OR. &
     &    KID.EQ.'?'.OR.KID.EQ.'#') THEN
         MODE=1
         RETURN
      ENDIF

      KID=' '
      MODE=0
      RETURN

    2 WRITE(6,'(A)') 'XX INPUT ERROR !'
      MODE=3
      RETURN

    3 KID='Q'
      MODE=1
      RETURN
      END SUBROUTINE TASK_KLIN

!     ****** INPUT PARAMETERS ******

      SUBROUTINE TASK_PARM(MODE,KWD,KIN,XXNLIN,XXPLST,IERR)

!     MODE=0 : standard namelist input  KIN: not used
!     MODE=1 : namelist file input      KIN: file name
!     MODE=2 : namelist line input      KIN: data string

!     KWD:     namelist name (not longer than 6 chars)
!     XXNLIN:  namelist subroutine XXNLIN(NID,IST,IERR)
!                      input  NID:  device number
!                      output IST:  io status
!                             IERR: error indicater
!     XXPLST:  error prompt subroutine XXPLST
!                      indicating valid namelist variables

!     IERR=0 : normal end
!     IERR=1 : namelist standard input end of file
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE

      IMPLICIT NONE

      INTEGER(4), INTENT(IN)   :: MODE
      INTEGER(4), INTENT(OUT)  :: IERR
      CHARACTER(LEN=*) , INTENT(IN) :: KWD, KIN
      CHARACTER(LEN=80):: LINE
      CHARACTER(LEN=94):: KNLINE
      CHARACTER(LEN=6) :: KNL
      LOGICAL          :: LEX
      INTEGER(4)       :: IST, KL, KL1, KL2
      EXTERNAL XXNLIN,XXPLST

      IERR=0
      KNL=KWD
      LINE=KIN

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

      ELSEIF(MODE.EQ.1) THEN

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
         WRITE(6,'(A,A,A)')  &
     &        '## FILE (',LINE(1:KL),') IS ASSIGNED FOR PARM INPUT'

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
            write(6,*) KNLINE
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
