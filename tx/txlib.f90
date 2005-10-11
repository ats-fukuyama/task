!     $Id$

!***************************************************************
!
!   For no '*** MATH LIBRARY ERROR 14: DEXP(X) UNDERFLOW'
!
!***************************************************************

REAL(8) FUNCTION EXPV(X)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: X

  IF (X < -708) THEN
     EXPV = 0.D0
  ELSE
     EXPV = EXP(X)
  END IF

  RETURN
END FUNCTION EXPV

!***************************************************************
!
!   SUBROUTINE APpend Integer TO Strings
!     INPUT  : STR, NSTR, I
!              STR(NSTR(original)+1:NSTR(return))
!              NSTR : Number of STR. First, NSTR = 0.
!     OUTPUT : STR, NSTR
!
!***************************************************************

SUBROUTINE APITOS(STR, NSTR, I)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: I
  INTEGER, INTENT(INOUT) :: NSTR
  CHARACTER(*), INTENT(INOUT) :: STR

  INTEGER :: J, NSTRI
  CHARACTER(25) :: KVALUE

  WRITE(KVALUE,'(I25)') I
  DO J = 1, 25
     IF (KVALUE(J:J) /= ' ') EXIT
  END DO
  NSTRI = 25 - J + 1
  STR(NSTR+1:NSTR+NSTRI) = KVALUE(J:25)
  NSTR = NSTR + NSTRI

  RETURN
END SUBROUTINE APITOS

!***************************************************************
!
!  SUBROUTINE APpend Text TO Strings
!     INPUT  : STR, NSTR, TX
!              STR(NSTR(original+1):NSTR(return))
!              NSTR : Number of STR. First, NSTR = 0.
!              TX : Delimited text
!     OUTPUT : STR, NSTR
!
!***************************************************************

SUBROUTINE APTTOS(STR, NSTR, TX)

  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: NSTR
  CHARACTER(*), INTENT(IN) :: TX
  CHARACTER(*), INTENT(INOUT) :: STR

  INTEGER :: NTX

  NTX = LEN_TRIM(TX)
  STR(NSTR+1:NSTR+NTX) = TX(1:NTX)
  NSTR = NSTR + NTX

  RETURN
END SUBROUTINE APTTOS

!***************************************************************
!
!  SUBROUTINE APpend Strings TO Strings
!     INPUT  : STR, NSTR, INSTR, NINSTR
!              STR(NSTR(original+1):NSTR(return))
!              NSTR : Number of STR. First, NSTR = 0.
!              NINSTR : Number of INSTR
!     OUTPUT : STR, NSTR
!
!***************************************************************

SUBROUTINE APSTOS(STR, NSTR, INSTR, NINSTR)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NINSTR
  CHARACTER(*), INTENT(IN) :: INSTR
  INTEGER, INTENT(INOUT) :: NSTR
  CHARACTER(*), INTENT(INOUT) :: STR

  STR(NSTR+1:NSTR+NINSTR) = INSTR(1:NINSTR)
  NSTR = NSTR + NINSTR

  RETURN
END SUBROUTINE APSTOS

!***************************************************************
!
!  SUBROUTINE APpend Double precision real number TO Strings
!     INPUT  : STR, NSTR, D, FORM
!              STR(NSTR(original+1):NSTR(return))
!              NSTR : Number of STR. First, NSTR = 0.
!              FORM : '{D|E|F|G}n' or '*'
!     OUTPUT : STR, NSTR
!
!***************************************************************

SUBROUTINE APDTOS(STR, NSTR, D, FORM)

  IMPLICIT NONE
  REAL(8), INTENT(IN) :: D
  CHARACTER(*), INTENT(IN) :: FORM
  INTEGER, INTENT(INOUT) :: NSTR
  CHARACTER(*), INTENT(INOUT) :: STR

  INTEGER(1) :: IND
  INTEGER :: L, IS, IE, NSTRD, IST
  CHARACTER(10) :: KFORM, KVALUE*25

  L = LEN(FORM)
  IF      (L == 0) THEN
     WRITE(6,*) '### ERROR(APDTOS) : Invalid Format : "', FORM , '"'
     NSTRD = 0
     RETURN
  ELSE IF (L == 1) THEN
     IF (FORM(1:1) == '*') THEN
        WRITE(KVALUE,*) SNGL(D)
     ELSE
        WRITE(6,*) '### ERROR(APDTOS) : Invalid Format : "', FORM , '"'
        NSTRD = 0
        RETURN
     END IF
  ELSE
     READ(FORM(2:2),'(I1)',IOSTAT=IST) IND
     IF (IST > 0) THEN
        WRITE(6,*) '### ERROR(APDTOS) : Invalid Format : "', FORM , '"'
        NSTRD = 0
        RETURN
     END IF
     IF (FORM(1:1) == 'F') THEN
        WRITE(KFORM,'(A,I2,A)') '(F25.', IND, ')'
        WRITE(KVALUE,KFORM) D
     ELSE IF (FORM(1:1) == 'D' .OR. FORM(1:1) == 'E' &
          &            .OR. FORM(1:1) == 'G') THEN
        WRITE(KFORM,'(3A,I2,A)') '(1P', FORM(1:1), '25.', IND, ')'
        WRITE(KVALUE,KFORM) D
     ELSE
        WRITE(6,*) '### ERROR(APDTOS) : Invalid Format : "', FORM , '"'
        NSTRD = 0
        RETURN
     END IF
  END IF

  DO IS = 1, 25
     IF (KVALUE(IS:IS) /= ' ') EXIT
  END DO
     
  IE = IS
20 CONTINUE
  IE = IE + 1
  IF (KVALUE(IE:IE) /= ' ' .AND. IE < 25) GOTO 20
  IF (KVALUE(IE:IE) /= ' ' .AND. IE == 25) IE = 25 + 1

  IF (KVALUE(IS:IS) == '-') THEN
     IF (IS > 1 .AND. KVALUE(IS+1:IS+1) == '.') THEN
        KVALUE(IS-1:IS-1) = '-'
        KVALUE(IS  :IS  ) = '0'
        IS = IS - 1
     END IF
  ELSE IF (KVALUE(IS:IS) == '.') THEN
     IF (IS > 1) THEN
        KVALUE(IS-1:IS-1) = '0'
        IS = IS - 1
     END IF
  END IF

  NSTRD = IE - IS
  STR(NSTR+1:NSTR+NSTRD) = KVALUE(IS:IE-1)
  NSTR = NSTR + NSTRD

  RETURN
END SUBROUTINE APDTOS

!***************************************************************
!
!  SUBROUTINE APpend Real number TO Strings
!     INPUT  : STR, NSTR, GR, FORM
!              STR(NSTR(original+1):NSTR(return))
!              NSTR : Number of STR. First, NSTR = 0.
!              FORM : '{D|E|F|G}n' or '*'
!     OUTPUT : STR, NSTR
!
!***************************************************************

SUBROUTINE APRTOS(STR, NSTR, GR, FORM)

  IMPLICIT NONE
  REAL, INTENT(IN) :: GR
  CHARACTER(*), INTENT(IN) :: FORM
  INTEGER, INTENT(INOUT) :: NSTR
  CHARACTER(*), INTENT(INOUT) :: STR

  REAL(8) :: D

  D = DBLE(GR)
  CALL APDTOS(STR, NSTR, D, FORM)

  RETURN
END SUBROUTINE APRTOS

!***************************************************************
!
!   Convert Strings to Upper Case
!
!***************************************************************

SUBROUTINE TOUPPER(KTEXT)

  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT) ::  KTEXT

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
!   CPU used time
!
!***************************************************************

! Not really cpu time but wallclock time.
! Should work with all Fortran90 systems.
! Note that TIME has to be default real kind (F95 standard allows any kind).

SUBROUTINE CPU_TIME(TIME)
  REAL, INTENT(OUT) :: TIME

  INTEGER :: COUNT, COUNT_RATE

  CALL SYSTEM_CLOCK(COUNT,COUNT_RATE)
  IF (COUNT_RATE /= 0) THEN
     TIME = REAL(COUNT)/COUNT_RATE
  ELSE
     TIME = -1
  END IF
END SUBROUTINE CPU_TIME
