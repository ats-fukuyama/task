C
C     ***** OPEN FILE FOR READ *****
C
      SUBROUTINE FROPEN(NFL,KNAMFL,MODEF,MODEP,KPR,IERR)
C
C     INPUT:
C        NFL    : FILE DEVICE NUMBER
C        KNAMFL : FILE NAME
C        MODEF  : 0 : UNFORMATTED
C                 1 : FORMATTED
C        MODEP  : 0 : WITHOUT PROMPT
C                 1 : WITH FILE NAME INPUT
C        KPR    : PROMPT
C
C     OUTPUT:
C        IERR   : ERROR CODE
C                 0 : NO ERROR
C                 1 : EOF IN READ FILE NAME
C                 2 : CANCEL WITH BLANK FILE NAME
C                 3 : UNDEFINED MODEP
C                 4 : FILE NAME ERROR
C                 5 : UNDEFINED MODEF
C                 6 : OLD FILE OPEN ERROR
C                 7 : OLD FILE NOT FOUND
C
      CHARACTER KNAMFL*80,KNAM*80,KPR*(*)
      LOGICAL LEX
C
      IF(MODEP.EQ.1) THEN
         KNAM=KNAMFL
         CALL KTRIM(KNAM,KL)
    1    WRITE(6,*) '#',KPR,'> INPUT : LOAD FILE NAME : ',KNAM(1:KL)
         READ(5,'(A80)',ERR=1,END=9001) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMFL=KNAM
         IF(KNAMFL(1:2).EQ.'  ') GOTO 9002
      ELSEIF(MODEP.NE.0) THEN
         WRITE(6,*) 'XX FROPEN: UNKNOWN MODEP : MODEP=',MODEP
         GOTO 9003
      ENDIF
C
      INQUIRE(FILE=KNAMFL,EXIST=LEX,ERR=9004)
      KNAM=KNAMFL
      CALL KTRIM(KNAM,KL)
      IF(LEX) THEN
         IF(MODEF.EQ.0) THEN
            OPEN(NFL,FILE=KNAMFL,IOSTAT=IST,STATUS='OLD',ERR=20,
     &           FORM='UNFORMATTED')
         ELSEIF(MODEF.EQ.1) THEN
            OPEN(NFL,FILE=KNAMFL,IOSTAT=IST,STATUS='OLD',ERR=20,
     &           FORM='FORMATTED')
         ELSE
            WRITE(6,*) 'XX FROPEN: UNKNOWN MODEF : MODEF=',MODEF
            GOTO 9005
         ENDIF
         KNAM=KNAMFL
         WRITE(6,*) '# OLD FILE (',KNAM(1:KL),') IS ASSIGNED FOR INPUT.'
         GOTO 9000
C
   20    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 9006
      ELSE
         WRITE(6,*) 'XX FILE (',KNAM(1:KL),') NOT FOUND'
         GOTO 9007
      ENDIF
C
 9000 IERR=0
      RETURN
C
 9001 IERR=1
      RETURN
C
 9002 IERR=2
      RETURN
C
 9003 IERR=3
      RETURN
C
 9004 IERR=4
      RETURN
C
 9005 IERR=5
      RETURN
C
 9006 IERR=6
      RETURN
C
 9007 IERR=7
      RETURN
      END
C
C     ***** OPEN FILE FOR WRITE *****
C
      SUBROUTINE FWOPEN(NFL,KNAMFL,MODEF,MODEP,KPR,IERR)
C
C     INPUT:
C        NFL    : FILE DEVICE NUMBER
C        KNAMFL : FILE NAME
C        MODEF  : 0 : UNFORMATTED
C                 1 : FORMATTED
C        MODEP  : 0 : WITHOUT PROMPT
C                 1 : WITH FILE NAME INPUT
C                 2 : ASK NEW NAME, IF FILE EXISTS
C                 3 : CONFIRM, IF FILE EXISTS
C                 4 : ERROR, IF FILE EXISTS
C
C     OUTPUT:
C        IERR   : ERROR CODE
C                 0 : NO ERROR
C                 1 : EOF IN READ FILE NAME
C                 2 : CANCEL WITH BLANK FILE NAME
C                 3 : UNDEFINED MODEP
C                 4 : FILE NAME ERROR
C                 5 : UNDEFINED MODEF
C                 6 : OLD FILE OPEN ERROR
C                 7 : OLD FILE EXISTS
C
      CHARACTER KNAMFL*80,KNAM*80,KPR*(*),KID*1
      LOGICAL LEX
C
      MODEPI=MODEP
C
 1000 IF(MODEPI.EQ.1) THEN
         KNAM=KNAMFL
         CALL KTRIM(KNAM,KL)
    1    WRITE(6,*) '#',KPR,'> INPUT : SAVE FILE NAME : ',KNAM(1:KL)
         READ(5,'(A80)',ERR=1,END=9001) KNAM
         IF(KNAM(1:2).NE.'/ ') KNAMFL=KNAM
         IF(KNAM(1:2).EQ.'  ') GOTO 9002
      ENDIF
C
      INQUIRE(FILE=KNAMFL,EXIST=LEX,ERR=9004)
      KNAM=KNAMFL
      CALL KTRIM(KNAM,KL)
C
      IF(LEX) THEN
         IF(MODEPI.EQ.1) THEN
    2       WRITE(6,*) '# OLD FILE (',KNAM(1:KL),
     &                 ') IS GOING TO BE OVERWRITTEN'
            WRITE(6,*) '  ARE YOU SURE ? (Y/N)'
            READ(5,'(A1)',ERR=2,END=9001) KID
            CALL GUCPTL(KID)
            IF(KID.EQ.'N') GOTO 1000
         ELSEIF(MODEPI.EQ.2) THEN
            MODEPI=1
            GOTO 1000
         ELSEIF(MODEPI.EQ.3) THEN
    3       WRITE(6,*) '# OLD FILE (',KNAM(1:KL),
     &                 ') IS GOING TO BE OVERWRITTEN'
            WRITE(6,*) '  ARE YOU SURE ? (Y/N)'
            READ(5,'(A1)',ERR=3,END=9001) KID
            CALL GUCPTL(KID)
            IF(KID.EQ.'N') GOTO 9007
         ELSEIF(MODEPI.EQ.4) THEN
            WRITE(6,*) 'XX FWOPEN: FILE ALREADY EXISTS.'
            GOTO 9007
         ELSEIF(MODEPI.NE.0) THEN
            WRITE(6,*) 'XX FWOPEN: UNKNOWN MODEP : MODEP=',MODEP
            GOTO 9003
         ENDIF
C
         IF(MODEF.EQ.0) THEN
            OPEN(NFL,FILE=KNAMFL,IOSTAT=IST,STATUS='OLD',ERR=10,
     &           FORM='UNFORMATTED')
         ELSEIF(MODEF.EQ.1) THEN
            OPEN(NFL,FILE=KNAMFL,IOSTAT=IST,STATUS='OLD',ERR=10,
     &           FORM='FORMATTED')
         ELSE
            WRITE(6,*) 'XX FEOPEN: UNKNOWN MODEF : MODEF=',MODEF
            GOTO 9005
         ENDIF
         WRITE(6,*) '# OLD FILE (',KNAM(1:KL),
     &                 ') IS ASSIGNED FOR OUTPUT.'
         GOTO 9000
C
   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 9006
      ELSE
         IF(MODEF.EQ.0) THEN
            OPEN(21,FILE=KNAMFL,IOSTAT=IST,STATUS='NEW',ERR=20,
     &           FORM='UNFORMATTED')
         ELSEIF(MODEF.EQ.1) THEN
            OPEN(21,FILE=KNAMFL,IOSTAT=IST,STATUS='NEW',ERR=20,
     &           FORM='FORMATTED')
         ELSE
            WRITE(6,*) 'XX FEOPEN: UNKNOWN MODEF : MODEF=',MODEF
            GOTO 9005
         ENDIF
         WRITE(6,*) '# NEW FILE (',KNAM(1:KL),') IS CREATED FOR OUTPUT.'
         GOTO 9000
C
   20    WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 9006
      ENDIF
C
 9000 IERR=0
      RETURN
C
 9001 IERR=1
      RETURN
C
 9002 IERR=2
      RETURN
C
 9003 IERR=3
      RETURN
C
 9004 IERR=4
      RETURN
C
 9005 IERR=5
      RETURN
C
 9006 IERR=6
      RETURN
C
 9007 IERR=7
      RETURN
      END
