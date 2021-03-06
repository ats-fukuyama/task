! libchar.f90

MODULE libchar

  PRIVATE
  PUBLIC kkindex
  PUBLIC ksplit
  PUBLIC ktrim
  PUBLIC kextr
  PUBLIC kmatch
  PUBLIC kkint
  PUBLIC toupper
  PUBLIC tolower
  PUBLIC linsep
  PUBLIC lenx
  
CONTAINS
  
!     ***** FIND CHARACTER POSITION *****
!     ===== this routine is no longer necessary =====
!     ===== use 'INDEX'                         =====

      INTEGER FUNCTION KKINDEX(KKLINEX,KID)

      IMPLICIT NONE
      INTEGER, PARAMETER :: KSLEN=80
      CHARACTER(LEN=*), INTENT(IN) :: KKLINEX
      CHARACTER(LEN=1), INTENT(IN) :: KID
      CHARACTER(LEN=80)            :: KKLINE
      INTEGER :: I

      KKLINE=KKLINEX
      DO I=1,KSLEN
         IF(KKLINE(I:I).EQ.KID) GOTO 100
      ENDDO
      KKINDEX=0
      RETURN

  100 KKINDEX=I
      RETURN
      END FUNCTION KKINDEX

!     ***** SEPARATE STRING AT CHAR *****

      SUBROUTINE KSPLIT(KKLINE,KID,KKLINE1,KKLINE2)

      IMPLICIT NONE
      CHARACTER(LEN=*),  INTENT(IN)  :: KKLINE
      CHARACTER(LEN=1),  INTENT(IN)  :: KID
      CHARACTER(LEN=80), INTENT(OUT) :: KKLINE1, KKLINE2
      INTEGER :: I

      I = INDEX(KKLINE,KID)
      IF(I == 0) THEN
         KKLINE1 = KKLINE
         KKLINE2 = ' '
      ELSEIF(I == 1) THEN
         KKLINE1 = ' '
         KKLINE2 = KKLINE(I+1:)
      ELSE
         KKLINE1 = KKLINE(1:I-1)
         KKLINE2 = KKLINE(I+1:)
      END IF

      END SUBROUTINE KSPLIT

!     ***** TRIM STRING AND RETURN LENGTH *****
!     ===== this routine is no longer necessary =====
!     ===== use 'LEN_TRIM'                      =====

      SUBROUTINE KTRIM(KKLINEX,KL)

      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(INOUT):: KKLINEX
      INTEGER, INTENT(OUT)        :: KL
      INTEGER, PARAMETER          :: KSLEN=80
      CHARACTER(LEN=80)    :: KKLINE, KKLINE1
      INTEGER           :: I, ITOP, IBOTTOM

      KKLINE=KKLINEX
      DO I=1,KSLEN
         IF(KKLINE(I:I).NE.' ') GOTO 100
      ENDDO
  100 ITOP=I

      DO I=KSLEN,1,-1
         IF(KKLINE(I:I).NE.' ') GOTO 200
      ENDDO
  200 IBOTTOM=I

      IF(ITOP.GT.KSLEN) THEN
         KKLINE=' '
         KL=0
      ELSE IF(IBOTTOM.LT.1) THEN
         KKLINE=' '
         KL=0
      ELSE
         KKLINE1=KKLINE(ITOP:IBOTTOM)
         KKLINE=KKLINE1
         KL=IBOTTOM-ITOP+1
      ENDIF
      KKLINEX=KKLINE
      RETURN
      END SUBROUTINE KTRIM

!     ***** EXTRACT STRING FROM KSTRING *****

      SUBROUTINE KEXTR(KKLINEX)

      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(INOUT) :: KKLINEX
      CHARACTER(LEN=80) :: KKLINE
      CHARACTER(LEN=1)  :: KID
      INTEGER        :: INCL

      KKLINE=KKLINEX
      KID=KKLINE(1:1)
      INCL=INDEX(KKLINE(2:),KID)
      IF(INCL <= 1) THEN
         KKLINEX=KKLINE(2:)
      ELSE
         KKLINEX=KKLINE(2:INCL)
      ENDIF
      RETURN
      END SUBROUTINE KEXTR

!     ***** RETURN TRUE WHEN STRINGS MATCH *****

      LOGICAL FUNCTION KMATCH(KKLINE1X,KKLINE2X)

      IMPLICIT NONE
      INTEGER, PARAMETER :: KSLEN=80
      CHARACTER(LEN=*),INTENT(IN) :: KKLINE1X, KKLINE2X
      CHARACTER(LEN=80)           :: KKLINE1,  KKLINE2
      INTEGER                  :: I

      KKLINE1=KKLINE1X
      KKLINE2=KKLINE2X
      DO I=1,KSLEN
         IF(KKLINE1(I:I).NE.KKLINE2(I:I)) GO TO 100
      ENDDO
      KMATCH=.TRUE.
      RETURN

  100 KMATCH=.FALSE.
      RETURN
      END FUNCTION KMATCH

      SUBROUTINE KKINT(I,K)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: I
      CHARACTER(LEN=*),INTENT(OUT):: K
      CHARACTER(LEN=8):: LINE
      INTEGER:: L

      WRITE(LINE,'(I8)') I
      L=0
1     L=L+1
      IF(LINE(L:L).EQ.' ') GOTO 1
      K=LINE(L:8)
      RETURN
      END SUBROUTINE KKINT

!***************************************************************
!
!   Convert Strings to Upper Case
!
!***************************************************************

SUBROUTINE TOUPPER(KTEXT)

  implicit none
  character(len=*), INTENT(INOUT) ::  KTEXT

  INTEGER :: NCHAR, I, ID

  NCHAR = LEN(KTEXT)
  DO I = 1, NCHAR
     ID=IACHAR(KTEXT(I:I))
     IF(ID >= 97 .AND. ID <= 122) ID = ID - 32
     KTEXT(I:I)=ACHAR(ID)
  END DO

  RETURN
END SUBROUTINE TOUPPER

!***************************************************************
!
!   Convert Strings to Lower Case
!
!***************************************************************

SUBROUTINE TOLOWER(KTEXT)

  implicit none
  character(len=*), INTENT(INOUT) ::  KTEXT

  INTEGER :: NCHAR, I, ID

  NCHAR = LEN(KTEXT)
  DO I = 1, NCHAR
     ID=IACHAR(KTEXT(I:I))
     IF(ID >= 65 .AND. ID <= 90) ID = ID + 32
     KTEXT(I:I)=ACHAR(ID)
  END DO

  RETURN
END SUBROUTINE TOLOWER

  ! *** separate line string by separator and space ***

SUBROUTINE linsep(ckey,csep,nk,mjs,mje,ndm)

    IMPLICIT NONE
!
!::arguments
    CHARACTER(LEN=*),INTENT(IN):: ckey  ! input line
    CHARACTER(LEN=*),INTENT(IN):: csep  ! separator
    INTEGER,INTENT(IN):: ndm            ! maximumm number of separated list
    INTEGER,INTENT(OUT):: nk            ! number of words
    INTEGER,INTENT(OUT):: mjs(ndm)      ! start char number of a word
    INTEGER,INTENT(OUT):: mje(ndm)      ! end char number of a word
!
!::local varaibales
    INTEGER:: nmax,ii,js,je,i,is
    CHARACTER(LEN=1):: ctab
!
    ctab = CHAR(9)
    nmax = LEN_TRIM(ckey)
!
    ii = 0
!
!::start
    is = 0
    DO i = 1, nmax
       is = i
       IF( ckey(i:i).EQ.ctab ) CYCLE
       IF( ckey(i:i).NE." " )  EXIT
    ENDDO
    IF( is.EQ.0 ) GOTO 100
!
!::separateor
    IF( INDEX(csep,ckey(i:i)).GT.0 ) is = is+1
    js = is
    DO i = is, nmax
       IF( INDEX(csep,ckey(i:i)).GT.0 ) THEN
          ii = ii + 1
          IF( ii.GT.ndm ) GOTO 910
          je = i - 1
          je = LEN_TRIM(ckey(1:je))
!
!--Note.    xs_units = "      ",  ==> " "
          IF( ckey(je:je).EQ.'"' .AND. ckey(i:i).EQ.'"' ) je = js+ 1
!--
          mjs(ii) = js
          mje(ii) = je
          IF( js.GT.je ) THEN
             ii = ii - 1
          ENDIF
          js = i+1
       ENDIF
    ENDDO
!
    IF( js.LE.nmax ) THEN
       ii = ii + 1
       IF( ii.GT.ndm ) GOTO 910
       mjs(ii) = js
       mje(ii) = nmax
    ENDIF
!
100 CONTINUE
    nk = ii
!
    RETURN
!
910 CONTINUE
    WRITE(6,'(/2x,"*** linsep ***  too many word  ",2i5)') ii,ndm
    STOP
END SUBROUTINE linsep

  ! *** almost equal to TRIM ***
  
  FUNCTION lenx(clin)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::  clin ! input character line
    INTEGER:: lenx                      ! lenggthof 
    INTEGER:: nmax,i,il

    nmax = LEN(clin)
    DO i = nmax, 1, -1
       il = i
       IF( clin(i:i).NE." " ) GOTO 120
    END DO
    il = 0
120 CONTINUE
    lenx = il

    RETURN
  END FUNCTION lenx
END MODULE libchar
