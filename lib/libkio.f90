! libkio.f90

MODULE libkio

  PRIVATE
  PUBLIC TASK_KLIN
  PUBLIC TASK_PARM

CONTAINS

!     ***** INPUT KID or LINE *****
!                   MODE=0: LINE INPUT
!                        1: KID INPUT (first one char: a-z,A-Z,#,?,!)
!                        2: PARM INPUT
!                        3: ERROR INPUT

      SUBROUTINE TASK_KLIN(LINE,KID,MODE,XXPARM)

      USE libchar
      IMPLICIT NONE
      INTEGER, INTENT(OUT)  :: MODE
      CHARACTER(LEN=80),INTENT(OUT) :: LINE
      CHARACTER(LEN=1), INTENT(OUT) :: KID
      INTEGER  :: I, ID, IERR, IKID
      EXTERNAL XXPARM

      READ(5,'(A80)',ERR=2,END=3) LINE    ! read one line

      KID=LINE(1:1)                  ! if first char is '!'
      IF(KID.EQ.'!') THEN           ! comment input
         MODE=0
         RETURN
      END IF

      ID=0                                ! if "=" is included, namelist
      DO I=1,80
         IF(LINE(I:I).EQ.'=') ID=1
      ENDDO
      IF(ID.EQ.1) THEN
         CALL XXPARM(2,LINE,IERR)
         MODE=2
         RETURN
      ENDIF

      KID=LINE(1:1)                  ! if first char is lower-case, captalized
      CALL toupper(KID)

      IF((KID.GE.'A'.AND.KID.LE.'Z').OR. &
     &    KID.EQ.'?'.OR.KID.EQ.'#'.OR.KID.EQ.'!') THEN  ! one char input
         MODE=1
         RETURN
      ENDIF

      KID=' '                                           ! line input
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

      INTEGER, INTENT(IN)   :: MODE
      INTEGER, INTENT(OUT)  :: IERR
      CHARACTER(LEN=*) , INTENT(IN) :: KWD, KIN
      CHARACTER(LEN=80):: LINE
      CHARACTER(LEN=94):: KNLINE
      CHARACTER(LEN=6) :: KNL
      LOGICAL          :: LEX
      INTEGER       :: IST, KL, KL1, KL2
      EXTERNAL XXNLIN,XXPLST

      IERR=0
      KNL=KWD
      LINE=KIN

      IF(MODE.EQ.0) THEN
    1    CONTINUE
         WRITE(6,'(A,A,A)') '## INPUT ',TRIM(KNL),' :'
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
         OPEN(25,FILE=LINE,IOSTAT=IST,STATUS='OLD',FORM='FORMATTED',ERR=9100)
         CALL XXNLIN(25,IST,IERR)
         IF(IERR.EQ.8) GOTO 9800
         IF(IERR.EQ.9) GOTO 9900
         CLOSE(25)
         WRITE(6,'(A,A,A)')  &
     &        '## FILE (',TRIM(LINE),') IS ASSIGNED FOR PARM INPUT'

      ELSEIF(MODE.EQ.2) THEN
         KNLINE=' &'//TRIM(KNL)//' '//TRIM(LINE)//' /'
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
         WRITE(6,'(A,A6,A2,A)') '## PARM INPUT ACCEPTED: ',KWD,': ',TRIM(LINE)
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

END MODULE libkio
