!     $Id$
module libraries
  implicit none
  private
  public :: EXPV, APITOS, APTTOS, APSTOS, APRTOS, TOUPPER

contains
!***************************************************************
!
!   For no '*** MATH LIBRARY ERROR 14: DEXP(X) UNDERFLOW'
!
!***************************************************************

  REAL(8) FUNCTION EXPV(X)

    REAL(8), INTENT(IN) :: X

    IF (X < -708.D0) THEN
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

    INTEGER, INTENT(IN) :: I
    INTEGER, INTENT(INOUT) :: NSTR
    character(len=*), INTENT(INOUT) :: STR

    INTEGER :: J, NSTRI
    character(len=25) :: KVALUE

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

    INTEGER, INTENT(INOUT) :: NSTR
    character(len=*), INTENT(IN) :: TX
    character(len=*), INTENT(INOUT) :: STR

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

    INTEGER, INTENT(IN) :: NINSTR
    character(len=*), INTENT(IN) :: INSTR
    INTEGER, INTENT(INOUT) :: NSTR
    character(len=*), INTENT(INOUT) :: STR

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

    REAL(8), INTENT(IN) :: D
    character(len=*), INTENT(IN) :: FORM
    INTEGER, INTENT(INOUT) :: NSTR
    character(len=*), INTENT(INOUT) :: STR

    INTEGER(1) :: IND
    INTEGER :: L, IS, IE, NSTRD, IST
    character(len=10) :: KFORM
    character(len=25) :: KVALUE

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
    DO
       IE = IE + 1
       IF (KVALUE(IE:IE) /= ' ' .AND. IE < 25) THEN
          CYCLE
       ELSE
          EXIT
       END IF
    END DO
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

    REAL, INTENT(IN) :: GR
    character(len=*), INTENT(IN) :: FORM
    INTEGER, INTENT(INOUT) :: NSTR
    character(len=*), INTENT(INOUT) :: STR

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

end module libraries

!*****************************************************************

module core_module
  use commons, only : nrmax, h
  implicit none
  public

contains

  function fem_integral(id,ne,a,actr) result(x)
!-------------------------------------------------------
!
!   Calculate "\int_{r_i}^{r_{i+1}} function(r) dr"
!      dr : mesh interval
!      a  : coefficient vector
!      u  : variable vector
!      w  : weighting vector
!
!   function(r) is classified as follows:
!      id = 0  : a * w
!      id = 1  : u * w
!      id = 2  : a * u * w
!      id = 3  :(a * u)'* w
!      id = 4  : u'* w
!      id = 5  : a * u'* w
!      id = 6  : a'* u * w
!      id = 7  : a'* u'* w
!      id = 8  : u * w'
!      id = 9  : a * u * w'
!      id = 10 :(a * u)'* w'
!      id = 11 : u'* w'
!      id = 12 : a * u'* w'
!      id = 13 : a'* u * w'
!      id = 14 : a'* u'* w'
!
!   where ' means the derivative of r
!
!  < input >
!     id       : mode select
!     ne       : current number of elements
!     nnode    : maximum of nodes
!     a(nnode) : coefficient vector, optional
!     actr     : element averaged value, optional
!  < output >
!     x(4)     : matrix of integrated values
!
!-------------------------------------------------------
    integer, intent(in) :: id, ne
    real(8), intent(in), dimension(0:nrmax), optional  :: a
    real(8), intent(in), optional :: actr
    integer :: node1, node2
    real(8) :: x(1:4), a1, a2
    
    node1 = ne-1  ; node2 = ne
    a1 = a(node1) ; a2 = a(node2)
    if(ne == 1.and.present(actr)) then
       a1 = actr
       a2 = actr
    end if

    select case(id)
    case(0)
       x(1) = h(ne) / 3.d0 * a1
       x(2) = h(ne) / 6.d0 * a2
       x(3) = h(ne) / 6.d0 * a1
       x(4) = h(ne) / 3.d0 * a2
    case(1)
       x(1) = h(ne) / 3.d0
       x(2) = h(ne) / 6.d0
       x(3) = h(ne) / 6.d0
       x(4) = h(ne) / 3.d0
    case(2)
       x(1) = ( 3.d0 * a1 +        a2) * h(ne) / 12.d0
       x(2) = (        a1 +        a2) * h(ne) / 12.d0
       x(3) = (        a1 +        a2) * h(ne) / 12.d0
       x(4) = (        a1 + 3.d0 * a2) * h(ne) / 12.d0
    case(3)
       x(1) = (-4.d0 * a1 +        a2) / 6.d0
       x(2) = (        a1 + 2.d0 * a2) / 6.d0
       x(3) = (-2.d0 * a1 -        a2) / 6.d0
       x(4) = (-       a1 + 4.d0 * a2) / 6.d0
    case(4)
       x(1) = -0.5d0
       x(2) =  0.5d0
       x(3) = -0.5d0
       x(4) =  0.5d0
    case(5)
       x(1) = (-2.d0 * a1 -        a2) / 6.d0
       x(2) = ( 2.d0 * a1 +        a2) / 6.d0
       x(3) = (-       a1 - 2.d0 * a2) / 6.d0
       x(4) = (        a1 + 2.d0 * a2) / 6.d0
    case(6)
       x(1) = (-a1 + a2) / 3.d0
       x(2) = (-a1 + a2) / 6.d0
       x(3) = (-a1 + a2) / 6.d0
       x(4) = (-a1 + a2) / 3.d0
    case(7)
       x(1) = ( a1 - a2) / (2.d0 * h(ne))
       x(2) = ( a1 - a2) / (2.d0 * h(ne))
       x(3) = (-a1 + a2) / (2.d0 * h(ne))
       x(4) = (-a1 + a2) / (2.d0 * h(ne))
    case(8)
       x(1) =-0.5d0
       x(2) =-0.5d0
       x(3) = 0.5d0
       x(4) = 0.5d0
    case(9)
       x(1) = (-2.d0 * a1 -        a2) / 6.d0
       x(2) = (-       a1 - 2.d0 * a2) / 6.d0
       x(3) = ( 2.d0 * a1 +        a2) / 6.d0
       x(4) = (        a1 + 2.d0 * a2) / 6.d0
    case(10)
       x(1) = a1 / h(ne)
       x(2) =-a2 / h(ne)
       x(3) =-a1 / h(ne)
       x(4) = a2 / h(ne)
    case(11)
       x(1) = 1.d0 / h(ne)**2
       x(2) =-1.d0 / h(ne)**2
       x(3) =-1.d0 / h(ne)**2
       x(4) = 1.d0 / h(ne)**2
    case(12)
       x(1) = ( a1 + a2) / (2.d0 * h(ne))
       x(2) = (-a1 - a2) / (2.d0 * h(ne))
       x(3) = (-a1 - a2) / (2.d0 * h(ne))
       x(4) = ( a1 + a2) / (2.d0 * h(ne))
    case(13)
       x(1) = ( a1 - a2) / (2.d0 * h(ne))
       x(2) = ( a1 - a2) / (2.d0 * h(ne))
       x(3) = (-a1 + a2) / (2.d0 * h(ne))
       x(4) = (-a1 + a2) / (2.d0 * h(ne))
    case(14)
       x(1) = (-a1 + a2) / h(ne)**3
       x(2) = ( a1 - a2) / h(ne)**3
       x(3) = ( a1 - a2) / h(ne)**3
       x(4) = (-a1 + a2) / h(ne)**3
    case default
       stop 'XX falut ID in fem_integral'
    end select

  end function fem_integral
end module core_module

