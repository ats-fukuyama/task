C     $Id$
C
C     ***************************************************************
C
C        Check negative value
C
C     ***************************************************************
C
      SUBROUTINE CHKNEG(XL, NRL, STR, XLNEG)
C
      IMPLICIT REAL*8(A-F, H, O-Z)
C
      DIMENSION XL(0:NRL)
      CHARACTER STR*(*), STR1*10
C
      DO I = 0, NRL
         IF (XL(I) .LT. 0.D0) THEN
            NSTR1 = 0
            CALL APITOS(STR1, NSTR1, I)
            NSTR = NSTRLEN(STR)
            WRITE(6,*) '### ERROR(CHKNEG) : ', STR(2:NSTR+1),
     &                 '(', STR1(1:NSTR1), ') is', SNGL(XL(I))
            XL(I) = XLNEG
         ENDIF
      ENDDO
C
      RETURN
      END
C
C     ***************************************************************
C
C     For no '*** MATH LIBRARY ERROR 14: DEXP(X) UNDERFLOW'
C
C     ***************************************************************
C
      FUNCTION EXPV(X)
C
      IMPLICIT REAL*8(A-F, H, O-Z)
C
      IF (X .LT. -708) THEN
         EXPV = 0.D0
      ELSE
         EXPV = EXP(X)
      ENDIF
C
      RETURN
      END
C
C     ***************************************************************
C
C        Length of strings   STR(2:NSTRLEN(A)+1)
C
C     ***************************************************************
C
      FUNCTION NSTRLEN(STR)
C
      CHARACTER STR*(*), DELIM*1
C
      DELIM = STR(1:1)
      DO I = 2, LEN(STR)
         IF (STR(I:I) .EQ. DELIM) THEN
            NSTRLEN = I - 2
            GOTO 10
         ENDIF
      ENDDO
      NSTRLEN = LEN(STR)
   10 CONTINUE
C
      RETURN
      END
C
C     ***************************************************************
C
C     SUBROUTINE APpend Integer TO Strings
C        INPUT  : NSTR, I
C                 NSTR : Number of STR. First, NSTR = 0.
C        OUTPUT : STR, NSTR
C                 STR(NSTR(original)+1:NSTR(return))
C
C     ***************************************************************
C
      SUBROUTINE APITOS(STR, NSTR, I)
      IMPLICIT REAL*8(A-F, H, O-Z)
      CHARACTER STR*(*), KVALUE*25
C
      WRITE(KVALUE,'(I25)') I
      J = 0
   10 CONTINUE
         J = J + 1
      IF (KVALUE(J:J) .EQ. ' ') GOTO 10
      NSTRI = 25 - J + 1
      STR(NSTR+1:NSTR+NSTRI) = KVALUE(J:25)
      NSTR = NSTR + NSTRI
C
      RETURN
      END
C
C     ***************************************************************
C
C     SUBROUTINE APpend Textx TO Strings
C        INPUT  : NSTR, TX
C                 NSTR : Number of STR. First, NSTR = 0.
C                 TX : Delimited text
C        OUTPUT : STR, NSTR
C                 STR(NSTR(original+1):NSTR(return))
C
C     ***************************************************************
C
      SUBROUTINE APTTOS(STR, NSTR, TX)
      IMPLICIT REAL*8(A-F, H, O-Z)
      CHARACTER TX*(*), STR*(*)
C
      NTX = NSTRLEN(TX)
      STR(NSTR+1:NSTR+NTX) = TX(2:NTX+1)
      NSTR = NSTR + NTX
C
      RETURN
      END
C
C     ***************************************************************
C
C     SUBROUTINE APpend Strings TO Strings
C        INPUT  : NSTR, INSTR, NINSTR
C                 NSTR : Number of STR. First, NSTR = 0.
C                 NINSTR : Number of INSTR
C        OUTPUT : STR, NSTR
C                 STR(NSTR(original+1):NSTR(return))
C
C     ***************************************************************
C
      SUBROUTINE APSTOS(STR, NSTR, INSTR, NINSTR)
      IMPLICIT REAL*8(A-F, H, O-Z)
      CHARACTER STR*(*), INSTR*(*)
C
      STR(NSTR+1:NSTR+NINSTR) = INSTR(1:NINSTR)
      NSTR = NSTR + NINSTR
C
      RETURN
      END
C
C     ***************************************************************
C
C     SUBROUTINE APpend Double precision real number TO Strings
C        INPUT  : NSTR, D, FORM
C                 NSTR : Number of STR. First, NSTR = 0.
C                 FORM : '{D|E|F|G}n' or '*'
C        OUTPUT : STR, NSTR
C                 STR(NSTR(original+1):NSTR(return))
C
C     ***************************************************************
C
      SUBROUTINE APDTOS(STR, NSTR, D, FORM)
      IMPLICIT REAL*8(A-F, H, O-Z)
      CHARACTER STR*(*), FORM*(*), KFORM*10, KVALUE*25
C
      L = LEN(FORM)
      IF      (L .EQ. 0) THEN
         GOTO 30
      ELSE IF (L .EQ. 1) THEN
         IF (FORM(1:1) .EQ. '*') THEN
            WRITE(KVALUE,*) SNGL(D)
         ELSE
            GOTO 30
         ENDIF
      ELSE
         READ(FORM(2:2),'(I1)',ERR=30) IND
         IF (FORM(1:1) .EQ. 'F') THEN
            WRITE(KFORM,'(A,I2,A)') '(F25.', IND, ')'
            WRITE(KVALUE,KFORM) D
         ELSE IF (FORM(1:1).EQ.'D' .OR. FORM(1:1).EQ.'E'
     &            .OR. FORM(1:1).EQ.'G') THEN
            WRITE(KFORM,'(3A,I2,A)') '(1P', FORM(1:1), '25.', IND, ')'
            WRITE(KVALUE,KFORM) D
         ELSE
            GOTO 30
         ENDIF
      ENDIF
C
      IS = 0
   10 CONTINUE
         IS = IS + 1
      IF (KVALUE(IS:IS) .EQ. ' ') GOTO 10
C
      IE = IS
   20 CONTINUE
         IE = IE + 1
      IF (KVALUE(IE:IE).NE.' ' .AND. IE.LT.25) GOTO 20
      IF (KVALUE(IE:IE).NE.' ' .AND. IE.EQ.25) IE = 25 + 1
C
      IF (KVALUE(IS:IS) .EQ. '-') THEN
         IF (IS.GT.1 .AND. KVALUE(IS+1:IS+1).EQ.'.') THEN
            KVALUE(IS-1:IS-1) = '-'
            KVALUE(IS  :IS  ) = '0'
            IS = IS - 1
         ENDIF
      ELSE IF (KVALUE(IS:IS) .EQ. '.') THEN
         IF (IS .GT. 1) THEN
            KVALUE(IS-1:IS-1) = '0'
            IS = IS - 1
         ENDIF
      ENDIF
C
      NSTRD = IE - IS
      STR(NSTR+1:NSTR+NSTRD) = KVALUE(IS:IE-1)
      NSTR = NSTR + NSTRD
C
      RETURN
C
   30 CONTINUE
      WRITE(6,*) '### ERROR(APDTOS) : Invalid Format : "', FORM , '"'
      NSTRD = 0
      RETURN
      END
C
C     ***************************************************************
C
C     SUBROUTINE APpend Real number TO Strings
C        INPUT  : NSTR, GR, FORM
C                 NSTR : Number of STR. First, NSTR = 0.
C                 FORM : '{D|E|F|G}n' or '*'
C        OUTPUT : STR, NSTR
C                 STR(NSTR(original+1):NSTR(return))
C
C     ***************************************************************
C
      SUBROUTINE APRTOS(STR, NSTR, GR, FORM)
      IMPLICIT REAL*8(A-F, H, O-Z)
      CHARACTER STR*(*), FORM*(*)
C
      D = DBLE(GR)
      CALL APDTOS(STR, NSTR, D, FORM)
C
      RETURN
      END
C
C     ***************************************************************
C
C        Convert Strings to Upper Case
C
C     ***************************************************************
C
      SUBROUTINE TOUPPER(KTEXT)
C
      CHARACTER KTEXT*(*)
      DIMENSION IASC(256)
C
      NCHAR = LEN(KTEXT)
      CALL CHRASC(KTEXT, IASC, NCHAR)
      DO I = 1, NCHAR
         ID = IASC(I)
         IF(ID .GE. 97 .AND. ID .LE. 122) IASC(I) = ID - 32
      ENDDO
      CALL ASCCHR(IASC, KTEXT, NCHAR)
C
      RETURN
      END
