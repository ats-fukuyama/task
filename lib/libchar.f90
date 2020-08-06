!     ***** FIND CHARACTER POSITION *****
!     ===== this routine is no longer necessary =====
!     ===== use 'INDEX'                         =====

      INTEGER(4) FUNCTION KKINDEX(KKLINEX,KID)

      IMPLICIT NONE
      INTEGER(4), PARAMETER :: KSLEN=80
      CHARACTER(LEN=*), INTENT(IN) :: KKLINEX
      CHARACTER(LEN=1), INTENT(IN) :: KID
      CHARACTER(LEN=80)            :: KKLINE
      INTEGER(4) :: I

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
      INTEGER(4) :: I

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
      INTEGER(4), INTENT(OUT)        :: KL
      INTEGER(4), PARAMETER          :: KSLEN=80
      CHARACTER(LEN=80)    :: KKLINE, KKLINE1
      INTEGER(4)           :: I, ITOP, IBOTTOM

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
      INTEGER(4)        :: INCL

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
      INTEGER(4), PARAMETER :: KSLEN=80
      CHARACTER(LEN=*),INTENT(IN) :: KKLINE1X, KKLINE2X
      CHARACTER(LEN=80)           :: KKLINE1,  KKLINE2
      INTEGER(4)                  :: I

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

  INTEGER(4) :: NCHAR, I, ID

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

  INTEGER(4) :: NCHAR, I, ID

  NCHAR = LEN(KTEXT)
  DO I = 1, NCHAR
     ID=IACHAR(KTEXT(I:I))
     IF(ID >= 65 .AND. ID <= 90) ID = ID + 32
     KTEXT(I:I)=ACHAR(ID)
  END DO

  RETURN
END SUBROUTINE TOLOWER
