!     $Id$
!
!     ***************************************************************
!
!        Check negative value
!
!     ***************************************************************
!
      SUBROUTINE CHKNEG(XL, NRL, STR, XLNEG)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NRL
      REAL(8), INTENT(IN) :: XLNEG
      REAL(8), DIMENSION(0:NRL), INTENT(INOUT) :: XL
      CHARACTER(*), INTENT(IN) :: STR

      INTEGER :: NSTRLEN, I, NSTR, NSTR1
      CHARACTER(10) :: STR1

      DO I = 0, NRL
         IF (XL(I) .LT. 0.D0) THEN
            NSTR1 = 0
            CALL APITOS(STR1, NSTR1, I)
            NSTR = NSTRLEN(STR)
            WRITE(6,*) '### ERROR(CHKNEG) : ', STR(2:NSTR+1), &
     &                 '(', STR1(1:NSTR1), ') is', SNGL(XL(I))
            XL(I) = XLNEG
         END IF
      END DO

      RETURN
      END
!
!     ***************************************************************
!
!     For no '*** MATH LIBRARY ERROR 14: DEXP(X) UNDERFLOW'
!
!     ***************************************************************
!
      REAL(8) FUNCTION EXPV(X)

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: X

      IF (X .LT. -708) THEN
         EXPV = 0.D0
      ELSE
         EXPV = EXP(X)
      END IF

      RETURN
      END
!
!     ***************************************************************
!
!        Length of strings   STR(2:NSTRLEN(A)+1)
!
!     ***************************************************************
!
      INTEGER FUNCTION NSTRLEN(STR)

      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: STR

      INTEGER :: I
      CHARACTER(1) :: DELIM

      DELIM = STR(1:1)
      DO I = 2, LEN(STR)
         IF (STR(I:I) .EQ. DELIM) THEN
            NSTRLEN = I - 2
            GOTO 10
         END IF
      END DO
      NSTRLEN = LEN(STR)
   10 CONTINUE

      RETURN
      END
!
!     ***************************************************************
!
!     SUBROUTINE APpend Integer TO Strings
!        INPUT  : NSTR, I
!                 NSTR : Number of STR. First, NSTR = 0.
!        OUTPUT : STR, NSTR
!                 STR(NSTR(original)+1:NSTR(return))
!
!     ***************************************************************
!
      SUBROUTINE APITOS(STR, NSTR, I)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I
      INTEGER, INTENT(INOUT) :: NSTR
      CHARACTER(*), INTENT(INOUT) :: STR

      INTEGER :: J, NSTRI
      CHARACTER(25) :: KVALUE

      WRITE(KVALUE,'(I25)') I
      J = 0
   10 CONTINUE
         J = J + 1
      IF (KVALUE(J:J) .EQ. ' ') GOTO 10
      NSTRI = 25 - J + 1
      STR(NSTR+1:NSTR+NSTRI) = KVALUE(J:25)
      NSTR = NSTR + NSTRI

      RETURN
      END
!
!     ***************************************************************
!
!     SUBROUTINE APpend Textx TO Strings
!        INPUT  : NSTR, TX
!                 NSTR : Number of STR. First, NSTR = 0.
!                 TX : Delimited text
!        OUTPUT : STR, NSTR
!                 STR(NSTR(original+1):NSTR(return))
!
!     ***************************************************************
!
      SUBROUTINE APTTOS(STR, NSTR, TX)

      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: NSTR
      CHARACTER(*), INTENT(IN) :: TX
      CHARACTER(*), INTENT(INOUT) :: STR

      INTEGER :: NTX, NSTRLEN

      NTX = NSTRLEN(TX)
      STR(NSTR+1:NSTR+NTX) = TX(2:NTX+1)
      NSTR = NSTR + NTX

      RETURN
      END
!
!     ***************************************************************
!
!     SUBROUTINE APpend Strings TO Strings
!        INPUT  : NSTR, INSTR, NINSTR
!                 NSTR : Number of STR. First, NSTR = 0.
!                 NINSTR : Number of INSTR
!        OUTPUT : STR, NSTR
!                 STR(NSTR(original+1):NSTR(return))
!
!     ***************************************************************
!
      SUBROUTINE APSTOS(STR, NSTR, INSTR, NINSTR)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NINSTR
      CHARACTER(*), INTENT(IN) :: INSTR
      INTEGER, INTENT(INOUT) :: NSTR
      CHARACTER(*), INTENT(INOUT) :: STR

      STR(NSTR+1:NSTR+NINSTR) = INSTR(1:NINSTR)
      NSTR = NSTR + NINSTR

      RETURN
      END
!
!     ***************************************************************
!
!     SUBROUTINE APpend Double precision real number TO Strings
!        INPUT  : NSTR, D, FORM
!                 NSTR : Number of STR. First, NSTR = 0.
!                 FORM : '{D|E|F|G}n' or '*'
!        OUTPUT : STR, NSTR
!                 STR(NSTR(original+1):NSTR(return))
!
!     ***************************************************************
!
      SUBROUTINE APDTOS(STR, NSTR, D, FORM)

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: D
      CHARACTER(*), INTENT(IN) :: FORM
      INTEGER, INTENT(INOUT) :: NSTR
      CHARACTER(*), INTENT(INOUT) :: STR

      INTEGER(1) :: IND
      INTEGER :: L, IS, IE, NSTRD
      CHARACTER(10) :: KFORM, KVALUE*25

      L = LEN(FORM)
      IF      (L .EQ. 0) THEN
         GOTO 30
      ELSE IF (L .EQ. 1) THEN
         IF (FORM(1:1) .EQ. '*') THEN
            WRITE(KVALUE,*) SNGL(D)
         ELSE
            GOTO 30
         END IF
      ELSE
         READ(FORM(2:2),'(I1)',ERR=30) IND
         IF (FORM(1:1) .EQ. 'F') THEN
            WRITE(KFORM,'(A,I2,A)') '(F25.', IND, ')'
            WRITE(KVALUE,KFORM) D
         ELSE IF (FORM(1:1).EQ.'D' .OR. FORM(1:1).EQ.'E' &
     &            .OR. FORM(1:1).EQ.'G') THEN
            WRITE(KFORM,'(3A,I2,A)') '(1P', FORM(1:1), '25.', IND, ')'
            WRITE(KVALUE,KFORM) D
         ELSE
            GOTO 30
         END IF
      ENDIF

      IS = 0
   10 CONTINUE
         IS = IS + 1
      IF (KVALUE(IS:IS) .EQ. ' ') GOTO 10

      IE = IS
   20 CONTINUE
         IE = IE + 1
      IF (KVALUE(IE:IE).NE.' ' .AND. IE.LT.25) GOTO 20
      IF (KVALUE(IE:IE).NE.' ' .AND. IE.EQ.25) IE = 25 + 1

      IF (KVALUE(IS:IS) .EQ. '-') THEN
         IF (IS.GT.1 .AND. KVALUE(IS+1:IS+1).EQ.'.') THEN
            KVALUE(IS-1:IS-1) = '-'
            KVALUE(IS  :IS  ) = '0'
            IS = IS - 1
         END IF
      ELSE IF (KVALUE(IS:IS) .EQ. '.') THEN
         IF (IS .GT. 1) THEN
            KVALUE(IS-1:IS-1) = '0'
            IS = IS - 1
         END IF
      END IF

      NSTRD = IE - IS
      STR(NSTR+1:NSTR+NSTRD) = KVALUE(IS:IE-1)
      NSTR = NSTR + NSTRD

      RETURN

   30 CONTINUE
      WRITE(6,*) '### ERROR(APDTOS) : Invalid Format : "', FORM , '"'
      NSTRD = 0
      RETURN
      END
!
!     ***************************************************************
!
!     SUBROUTINE APpend Real number TO Strings
!        INPUT  : NSTR, GR, FORM
!                 NSTR : Number of STR. First, NSTR = 0.
!                 FORM : '{D|E|F|G}n' or '*'
!        OUTPUT : STR, NSTR
!                 STR(NSTR(original+1):NSTR(return))
!
!     ***************************************************************
!
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
      END
!
!     ***************************************************************
!
!        Convert Strings to Upper Case
!
!     ***************************************************************
!
      SUBROUTINE TOUPPER(KTEXT)

      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) ::  KTEXT
      
      INTEGER :: NCHAR, I, ID
      REAL(8), DIMENSION(256) :: IASC

      NCHAR = LEN(KTEXT)
      CALL CHRASC(KTEXT, IASC, NCHAR)
      DO I = 1, NCHAR
         ID = IASC(I)
         IF(ID .GE. 97 .AND. ID .LE. 122) IASC(I) = ID - 32
      END DO
      CALL ASCCHR(IASC, KTEXT, NCHAR)

      RETURN
      END
