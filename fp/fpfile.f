C     $Id$
C
C     ***********************************************************
C
C           SAVE VELOCITY DISTRIBUTION DATA
C
C     ***********************************************************
C
      SUBROUTINE FPSAVE
C
      INCLUDE 'fpcomm.inc'
C
      CHARACTER*72 KNAM
C      CHARACTER*1 KID
      LOGICAL LEX
C
    1 WRITE(6,*) '# INPUT : FPDATA FILE NAME : ',KNAMFP
      READ(5,'(A72)',ERR=1,END=900) KNAM
      IF(KNAM(1:2).NE.'/ ') KNAMFP=KNAM
C
      INQUIRE(FILE=KNAMFP,EXIST=LEX)
      IF(LEX) THEN
C         WRITE(6,*) '# OLD FILE IS GOING TO BE OVERWRITTEN.  ',
C     &              'ARE YOU SURE {Y/N} ?'
C         READ(5,502) KID
C  502    FORMAT(A1)
C         CALL GUCPTL(KID)
C         IF(KID.NE.'Y') GOTO 1
         OPEN(21,FILE=KNAMFP,IOSTAT=IST,STATUS='OLD',ERR=10,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (',KNAMFP,') IS ASSIGNED FOR OUTPUT.'
         GOTO 30
   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 1
      ELSE
         OPEN(21,FILE=KNAMFP,IOSTAT=IST,STATUS='NEW',ERR=20,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# NEW FILE (',KNAMFP,') IS CREATED FOR OUTPUT.'
         GOTO 30
   20    WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 1
      ENDIF
C
   30 WRITE(21) RNFP0,RTFP0,PTFP0
      WRITE(21) DELR,DELP,DELTH
      WRITE(21) NRMAX,NPMAX,NTHMAX
      WRITE(21) (((F(NTH,NP,NR),NTH=1,NTHMAX),NP=1,NPMAX),NR=1,NRMAX)
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'
C
  900 RETURN
      END
C
C     ***********************************************************
C
C           LOAD VELOCITY DISTRIBUTION DATA
C
C     ***********************************************************
C
      SUBROUTINE FPLOAD
C
      INCLUDE 'fpcomm.inc'
C
      CHARACTER*72 KNAM
      LOGICAL LEX
C
    1 WRITE(6,*) '# INPUT : FPDATA FILE NAME : ',KNAMFP
      READ(5,'(A72)',ERR=1,END=900) KNAM
      IF(KNAM(1:2).NE.'/ ') KNAMFP=KNAM
C
      INQUIRE(FILE=KNAMFP,EXIST=LEX)
      IF(LEX) THEN
         OPEN(21,FILE=KNAMFP,IOSTAT=IST,STATUS='OLD',ERR=10,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (',KNAMFP,') IS ASSIGNED FOR INPUT.'
         GOTO 30
   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 1
      ELSE
         WRITE(6,*) 'XX FILE (',KNAMFP,') NOT FOUND'
         GOTO 1
      ENDIF
C
   30 READ(21) RNFP0,RTFP0,PTFP0
      READ(21) DELR,DELP,DELTH
      READ(21) NRMAX,NPMAX,NTHMAX
      READ(21) (((F(NTH,NP,NR),NTH=1,NTHMAX),NP=1,NPMAX),NR=1,NRMAX)
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
C
  900 RETURN
      END
