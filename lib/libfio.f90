!     ***** OPEN FILE FOR READ *****

      SUBROUTINE FROPEN(NFL,KNAMFL,MODEF,MODEP,KPR,IERR)

!     INPUT:
!        NFL    : FILE DEVICE NUMBER
!        KNAMFL : FILE NAME
!        MODEF  : 0 : UNFORMATTED
!                 1 : FORMATTED
!        MODEP  : 0 : WITHOUT PROMPT
!                 1 : WITH FILE NAME INPUT
!        KPR    : PROMPT

!     OUTPUT:
!        IERR   : ERROR CODE
!                 0 : NO ERROR
!                 1 : EOF IN READ FILE NAME
!                 2 : CANCEL WITH BLANK FILE NAME
!                 3 : UNDEFINED MODEP
!                 4 : FILE NAME ERROR
!                 5 : UNDEFINED MODEF
!                 6 : OLD FILE OPEN ERROR
!                 7 : OLD FILE NOT FOUND
!                 8 : EMPTY FILE NAME

      IMPLICIT NONE

      INTEGER(4),       INTENT(IN)   :: NFL, MODEF, MODEP
      INTEGER(4),       INTENT(OUT)  :: IERR
      CHARACTER(LEN=80),INTENT(INOUT):: KNAMFL
      CHARACTER(LEN=*), INTENT(IN)   :: KPR
      INTEGER(4)        :: KL, IST
      CHARACTER(LEN=80) :: KNAM
      LOGICAL           :: LEX

      KNAM=KNAMFL
      CALL KTRIM(KNAM,KL)
      IF(MODEP.EQ.0) THEN
         IF(KL.EQ.0) GOTO  9008
      ELSEIF(MODEP.EQ.1) THEN
    1    WRITE(6,*) '#',KPR,'> INPUT : LOAD FILE NAME : ',KNAM(1:KL)
         READ(5,'(A80)',ERR=1,END=9001) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMFL=KNAM
         IF(KNAMFL(1:2).EQ.'  ') GOTO 9002
      ELSE
         WRITE(6,*) 'XX FROPEN: UNKNOWN MODEP : MODEP=',MODEP
         GOTO 9003
      ENDIF

      INQUIRE(FILE=KNAMFL,EXIST=LEX,ERR=9004)
      KNAM=KNAMFL
      CALL KTRIM(KNAM,KL)
      IF(LEX) THEN
         IF(MODEF.EQ.0) THEN
            OPEN(NFL,FILE=KNAMFL,IOSTAT=IST,STATUS='OLD',ERR=20, &
     &           FORM='UNFORMATTED')
         ELSEIF(MODEF.EQ.1) THEN
            OPEN(NFL,FILE=KNAMFL,IOSTAT=IST,STATUS='OLD',ERR=20, &
     &           FORM='FORMATTED')
         ELSE
            WRITE(6,*) 'XX FROPEN: UNKNOWN MODEF : MODEF=',MODEF
            GOTO 9005
         ENDIF
         KNAM=KNAMFL
         WRITE(6,*) '# OLD FILE (',KNAM(1:KL),') IS ASSIGNED FOR INPUT.'
         GOTO 9000

   20    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 9006
      ELSE
         WRITE(6,*) 'XX FILE (',KNAM(1:KL),') NOT FOUND'
         GOTO 9007
      ENDIF

 9000 IERR=0
      RETURN

 9001 IERR=1
      RETURN

 9002 IERR=2
      RETURN

 9003 IERR=3
      RETURN

 9004 IERR=4
      RETURN

 9005 IERR=5
      RETURN

 9006 IERR=6
      RETURN

 9007 IERR=7
      RETURN

 9008 IERR=8
      RETURN
      END SUBROUTINE FROPEN

!     ***** OPEN FILE FOR WRITE *****

      SUBROUTINE FWOPEN(NFL,KNAMFL,MODEF,MODEP,KPR,IERR)

!     INPUT:
!        NFL    : FILE DEVICE NUMBER
!        KNAMFL : FILE NAME
!        MODEF  : 0 : UNFORMATTED
!                 1 : FORMATTED
!        MODEP  : 0 : WITHOUT PROMPT, ALWAYS OVERWRITE
!                 1 : WITHOUT PROMPT, CONFIRM, IF FILE EXISTS
!                 2 : WITHOUT PROMPT, ASK NEW NAME, IF FILE EXISTS
!                 3 : WITHOUT PROMPT, ERROR, IF FILE EXISTS
!                 4 : WITH FILE NAME INPUT, ALWAYS OVERWRITE
!                 5 : WITH FILE NAME INPUT, CONFIRM, IF FILE EXISTS
!                 6 : WITH FILE NAME INPUT, ASK NEW NAME, IF FILE EXISTS
!                 7 : WITH FILE NAME INPUT, ERROR, IF FILE EXISTS

!     OUTPUT:
!        IERR   : ERROR CODE
!                 0 : NO ERROR
!                 1 : EOF IN READ FILE NAME
!                 2 : CANCEL WITH BLANK FILE NAME
!                 3 : UNDEFINED MODEP
!                 4 : FILE NAME ERROR
!                 5 : UNDEFINED MODEF
!                 6 : OLD FILE OPEN ERROR
!                 7 : OLD FILE EXISTS
!                 8 : EMPTY FILE NAME

      IMPLICIT NONE

      INTEGER(4),       INTENT(IN)   :: NFL, MODEF, MODEP
      INTEGER(4),       INTENT(OUT)  :: IERR
      CHARACTER(LEN=80),INTENT(INOUT):: KNAMFL
      CHARACTER(LEN=*), INTENT(IN)   :: KPR
      INTEGER(4)        :: MODEPI, MODEPII, KL, IST
      CHARACTER(LEN=80) :: KNAM
      CHARACTER(LEN=1)  :: KID
      LOGICAL           :: LEX

      MODEPI=MODEP

      KNAM=KNAMFL
      CALL KTRIM(KNAM,KL)

 1000 IF(MODEPI.LE.3) THEN
         IF(KL.EQ.0) GOTO 9008
      ELSE
    1    WRITE(6,*) '#',KPR,'> INPUT : SAVE FILE NAME : ',KNAM(1:KL)
         READ(5,'(A80)',ERR=1,END=9001) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMFL=KNAM
         IF(KNAM(1:2).EQ.'  ') GOTO 9002
      ENDIF

      INQUIRE(FILE=KNAMFL,EXIST=LEX,ERR=9004)
      KNAM=KNAMFL
      CALL KTRIM(KNAM,KL)

      IF(LEX) THEN
         MODEPII=MOD(MODEPI,4)
         IF(MODEPII.EQ.0) THEN
            WRITE(6,*) '# OLD FILE (',KNAM(1:KL), &
     &                 ') WILL BE OVERWRITTEN'
         ELSEIF(MODEPII.EQ.1) THEN
    3       WRITE(6,*) '# OLD FILE (',KNAM(1:KL), &
     &                 ') IS GOING TO BE OVERWRITTEN'
            WRITE(6,*) '  ARE YOU SURE ? (Y/N)'
            READ(5,'(A1)',ERR=3,END=9001) KID
            CALL GUCPTL(KID)
            IF(KID.EQ.'N') GOTO 9007
         ELSEIF(MODEPII.EQ.2) THEN
            MODEPI=1
            GOTO 1000
         ELSEIF(MODEPII.EQ.3) THEN
            WRITE(6,*) 'XX FWOPEN: FILE ALREADY EXISTS.'
            GOTO 9007
         ELSE
            WRITE(6,*) 'XX FWOPEN: UNKNOWN MODEP : MODEP=',MODEP
            GOTO 9003
         ENDIF

         IF(MODEF.EQ.0) THEN
            OPEN(NFL,FILE=KNAMFL,IOSTAT=IST,STATUS='OLD',ERR=10, &
     &           FORM='UNFORMATTED')
         ELSEIF(MODEF.EQ.1) THEN
            OPEN(NFL,FILE=KNAMFL,IOSTAT=IST,STATUS='OLD',ERR=10, &
     &           FORM='FORMATTED')
         ELSE
            WRITE(6,*) 'XX FEOPEN: UNKNOWN MODEF : MODEF=',MODEF
            GOTO 9005
         ENDIF
         WRITE(6,*) '# OLD FILE (',KNAM(1:KL), &
     &                 ') IS ASSIGNED FOR OUTPUT.'
         GOTO 9000

   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 9006
      ELSE
         IF(MODEF.EQ.0) THEN
            OPEN(21,FILE=KNAMFL,IOSTAT=IST,STATUS='NEW',ERR=20, &
     &           FORM='UNFORMATTED')
         ELSEIF(MODEF.EQ.1) THEN
            OPEN(21,FILE=KNAMFL,IOSTAT=IST,STATUS='NEW',ERR=20, &
     &           FORM='FORMATTED')
         ELSE
            WRITE(6,*) 'XX FEOPEN: UNKNOWN MODEF : MODEF=',MODEF
            GOTO 9005
         ENDIF
         WRITE(6,*) '# NEW FILE (',KNAM(1:KL),') IS CREATED FOR OUTPUT.'
         GOTO 9000

   20    WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 9006
      ENDIF

 9000 IERR=0
      RETURN

 9001 IERR=1
      RETURN

 9002 IERR=2
      RETURN

 9003 IERR=3
      RETURN

 9004 IERR=4
      RETURN

 9005 IERR=5
      RETURN

 9006 IERR=6
      RETURN

 9007 IERR=7
      RETURN

 9008 IERR=8
      RETURN
      END SUBROUTINE FWOPEN
