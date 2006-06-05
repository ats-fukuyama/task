!     $Id$
module core_module
  use commons, only : nrmax, h, r
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
!      id = 15 : r * a * w
!      id = 16 : r * a * u * w
!      id = 17 : r * a'* u * w
!      id = 18 : r * a * u'* w
!      id = 19 : r * a * u'* w'
!
!      id = 20 : r**2 * a * w
!      id = 21 : r**2 * u * w
!      id = 22 : r**2 * a * u * w
!      id = 23 : r**2 * a'* u * w
!      id = 24 : r**2 * a * u'* w
!      id = 25 : r**2 * a * u'* w'
!      id = 26 : r**3 * a * u * w
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
    real(8) :: x(1:4), a1, a2, r1, r2
    
    node1 = ne-1  ; node2 = ne
    a1 = a(node1) ; a2 = a(node2)
    r1 = r(node1) ; r2 = r(node2)
    if(ne == 1.and.present(actr)) then
       a1 = actr
       a2 = actr
    end if

! minus sign number indicates lumped mass matrix.
    select case(id)
    case(-2)
       x(1) = ( 2.d0 * a1 +        a2) * h(ne) / 6.d0
       x(2) = 0.d0
       x(3) = 0.d0
       x(4) = (        a1 + 2.d0 * a2) * h(ne) / 6.d0
    case(-1)
       x(1) = h(ne) / 2.d0
       x(2) = 0.d0
       x(3) = 0.d0
       x(4) = h(ne) / 2.d0
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
       x(1) = h(ne) / 3.d0 * a1
       x(2) = h(ne) / 6.d0 * a1
       x(3) = h(ne) / 6.d0 * a2
       x(4) = h(ne) / 3.d0 * a2
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
       x(1) = 1.d0 / h(ne)
       x(2) =-1.d0 / h(ne)
       x(3) =-1.d0 / h(ne)
       x(4) = 1.d0 / h(ne)
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
       x(1) = (-a1 + a2) / h(ne)**2
       x(2) = ( a1 - a2) / h(ne)**2
       x(3) = ( a1 - a2) / h(ne)**2
       x(4) = (-a1 + a2) / h(ne)**2
    case(15)
       x(1) = (3.d0 * r1 +        r2) * h(ne) / 12.d0 * a1
       x(2) = (       r1 +        r2) * h(ne) / 12.d0 * a2
       x(3) = (       r1 +        r2) * h(ne) / 12.d0 * a1
       x(4) = (       r1 + 3.d0 * r2) * h(ne) / 12.d0 * a2
    case(16)
       x(1) = (12.d0*r1*a1 + 3.d0*r2*a1 + 3.d0*r1*a2 + 2.d0*r2*a2) * h(ne) / 60.d0
       x(2) = ( 3.d0*r1*a1 + 2.d0*r2*a1 + 2.d0*r1*a2 + 3.d0*r2*a2) * h(ne) / 60.d0
       x(3) = ( 3.d0*r1*a1 + 2.d0*r2*a1 + 2.d0*r1*a2 + 3.d0*r2*a2) * h(ne) / 60.d0
       x(4) = ( 2.d0*r1*a1 + 3.d0*r2*a1 + 3.d0*r1*a2 +12.d0*r2*a2) * h(ne) / 60.d0
    case(17)
       x(1) = (3.d0 * r1 +        r2) * (-a1 + a2) / 12.d0
       x(2) = (       r1 +        r2) * (-a1 + a2) / 12.d0
       x(3) = (       r1 +        r2) * (-a1 + a2) / 12.d0
       x(4) = (       r1 + 3.d0 * r2) * (-a1 + a2) / 12.d0
    case(18)
       x(1) = (-3.d0*r1*a1 - r2*a1 - r1*a2 -      r2*a2) / 12.d0
       x(2) = ( 3.d0*r1*a1 + r2*a1 + r1*a2 +      r2*a2) / 12.d0
       x(3) = (-     r1*a1 - r2*a1 - r1*a2 - 3.d0*r2*a2) / 12.d0
       x(4) = (      r1*a1 + r2*a1 + r1*a2 + 3.d0*r2*a2) / 12.d0
    case(19)
       x(1) = ( 2.d0*r1*a1 + r2*a1 + r1*a2 + 2.d0*r2*a2) / (6.d0 * h(ne))
       x(2) = (-2.d0*r1*a1 - r2*a1 - r1*a2 - 2.d0*r2*a2) / (6.d0 * h(ne))
       x(3) = (-2.d0*r1*a1 - r2*a1 - r1*a2 - 2.d0*r2*a2) / (6.d0 * h(ne))
       x(4) = ( 2.d0*r1*a1 + r2*a1 + r1*a2 + 2.d0*r2*a2) / (6.d0 * h(ne))
    case(20)
       x(1) = (6.d0*r1*r1 + 3.d0*r1*r2 +      r2*r2) * h(ne) / 30.d0 * a1
       x(2) = (3.d0*r1*r1 + 4.d0*r1*r2 + 3.d0*r2*r2) * h(ne) / 30.d0 * a2
       x(3) = (3.d0*r1*r1 + 4.d0*r1*r2 + 3.d0*r2*r2) * h(ne) / 30.d0 * a1
       x(4) = (     r1*r1 + 3.d0*r1*r2 + 6.d0*r2*r2) * h(ne) / 30.d0 * a2
    case(21)
       x(1) = (6.d0*r1*r1 + 3.d0*r1*r2 +      r2*r2) * h(ne) / 30.d0
       x(2) = (3.d0*r1*r1 + 4.d0*r1*r2 + 3.d0*r2*r2) * h(ne) / 30.d0
       x(3) = (3.d0*r1*r1 + 4.d0*r1*r2 + 3.d0*r2*r2) * h(ne) / 30.d0
       x(4) = (     r1*r1 + 3.d0*r1*r2 + 6.d0*r2*r2) * h(ne) / 30.d0
    case(22)
       x(1) = ( 10.d0*r1*r1*a1 + 4.d0*r1*r2*a1 +      r2*r2*a1 &
            &  + 2.d0*r1*r1*a2 + 2.d0*r1*r2*a2 +      r2*r2*a2) *h(ne) / 60.d0
       x(2) = (  2.d0*r1*r1*a1 + 2.d0*r1*r2*a1 +      r2*r2*a1 &
            &  +      r1*r1*a2 + 2.d0*r1*r2*a2 + 2.d0*r2*r2*a2) *h(ne) / 60.d0
       x(3) = (  2.d0*r1*r1*a1 + 2.d0*r1*r2*a1 +      r2*r2*a1 &
            &  +      r1*r1*a2 + 2.d0*r1*r2*a2 + 2.d0*r2*r2*a2) *h(ne) / 60.d0
       x(4) = (       r1*r1*a1 + 2.d0*r1*r2*a1 + 2.d0*r2*r2*a1 &
            &  +      r1*r1*a2 + 4.d0*r1*r2*a2 +10.d0*r2*r2*a2) *h(ne) / 60.d0
    case(23)
       x(1) = (- 6.d0*r1*r1*a1 - 3.d0*r1*r2*a1 -      r2*r2*a1 &
            &  + 6.d0*r1*r1*a2 + 3.d0*r1*r2*a2 +      r2*r2*a2) / 30.d0
       x(2) = (- 3.d0*r1*r1*a1 - 4.d0*r1*r2*a1 - 3.d0*r2*r2*a1 &
            &  + 3.d0*r1*r1*a2 + 4.d0*r1*r2*a2 + 3.d0*r2*r2*a2) / 60.d0
       x(3) = (- 3.d0*r1*r1*a1 - 4.d0*r1*r2*a1 - 3.d0*r2*r2*a1 &
            &  + 3.d0*r1*r1*a2 + 4.d0*r1*r2*a2 + 3.d0*r2*r2*a2) / 60.d0
       x(4) = (-      r1*r1*a1 - 3.d0*r1*r2*a1 - 6.d0*r2*r2*a1 &
            &  +      r1*r1*a2 + 3.d0*r1*r2*a2 + 6.d0*r2*r2*a2) / 30.d0
    case(24)
       x(1) = (-12.d0*r1*r1*a1 - 6.d0*r1*r2*a1 - 2.d0*r2*r2*a1 &
            &  - 3.d0*r1*r1*a2 - 4.d0*r1*r2*a2 - 3.d0*r2*r2*a2) / 60.d0
       x(2) = ( 12.d0*r1*r1*a1 + 6.d0*r1*r2*a1 + 2.d0*r2*r2*a1 &
            &  + 3.d0*r1*r1*a2 + 4.d0*r1*r2*a2 + 3.d0*r2*r2*a2) / 60.d0
       x(3) = (- 3.d0*r1*r1*a1 - 4.d0*r1*r2*a1 - 3.d0*r2*r2*a1 &
            &  - 2.d0*r1*r1*a2 - 6.d0*r1*r2*a2 -12.d0*r2*r2*a2) / 60.d0
       x(4) = (  3.d0*r1*r1*a1 + 4.d0*r1*r2*a1 + 3.d0*r2*r2*a1 &
            &  + 2.d0*r1*r1*a2 + 6.d0*r1*r2*a2 +12.d0*r2*r2*a2) / 60.d0
    case(25)
       x(1) = (  3.d0*r1*r1*a1 + 2.d0*r1*r2*a1 +      r2*r2*a1 &
            &  +      r1*r1*a2 + 2.d0*r1*r2*a2 + 3.d0*r2*r2*a2) / (12.d0 * h(ne))
       x(2) = (- 3.d0*r1*r1*a1 - 2.d0*r1*r2*a1 -      r2*r2*a1 &
            &  -      r1*r1*a2 - 2.d0*r1*r2*a2 - 3.d0*r2*r2*a2) / (12.d0 * h(ne))
       x(3) = (- 3.d0*r1*r1*a1 - 2.d0*r1*r2*a1 -      r2*r2*a1 &
            &  -      r1*r1*a2 - 2.d0*r1*r2*a2 - 3.d0*r2*r2*a2) / (12.d0 * h(ne))
       x(4) = (  3.d0*r1*r1*a1 + 2.d0*r1*r2*a1 +      r2*r2*a1 &
            &  +      r1*r1*a2 + 2.d0*r1*r2*a2 + 3.d0*r2*r2*a2) / (12.d0 * h(ne))
    case(26)
       x(1) = ( 60.d0*r1**3*a1 + 30.d0*r1**2*r2*a1 + 12.d0*r1*r2**2*a1 +  3.d0*r2**3*a1 &
            &  +10.d0*r1**3*a2 + 12.d0*r1**2*r2*a2 +  9.d0*r1*r2**2*a2 +  4.d0*r2**3*a2) &
            & / 420.d0
       x(2) = ( 10.d0*r1**3*a1 + 12.d0*r1**2*r2*a1 +  9.d0*r1*r2**2*a1 +  4.d0*r2**3*a1 &
            &  + 4.d0*r1**3*a2 +  9.d0*r1**2*r2*a2 + 12.d0*r1*r2**2*a2 + 10.d0*r2**3*a2) &
            & / 420.d0
       x(3) = ( 10.d0*r1**3*a1 + 12.d0*r1**2*r2*a1 +  9.d0*r1*r2**2*a1 +  4.d0*r2**3*a1 &
            &  + 4.d0*r1**3*a2 +  9.d0*r1**2*r2*a2 + 12.d0*r1*r2**2*a2 + 10.d0*r2**3*a2) &
            & / 420.d0
       x(4) = (  4.d0*r1**3*a1 +  9.d0*r1**2*r2*a1 + 12.d0*r1*r2**2*a1 + 10.d0*r2**3*a1 &
            &  + 3.d0*r1**3*a2 + 12.d0*r1**2*r2*a2 + 30.d0*r1*r2**2*a2 + 60.d0*r2**3*a2) &
            & / 420.d0
    case default
       stop 'XX falut ID in fem_integral'
    end select

  end function fem_integral
end module core_module

!*****************************************************************

module libraries
  implicit none
  private
  public :: EXPV, APITOS, APTTOS, APSTOS, APRTOS, TOUPPER, DERIV3SB, &
            F33, INTG_F, VALINT_SUB, INTG_P, LORENTZ, BISECTION

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

!***************************************************************
!
!   Derivative
!
!***************************************************************

  subroutine DERIV3SB(X,R,dX,NRMAX)
    use commons, only : NRM
    integer, intent(in) :: NRMAX 
    real(8), intent(in), dimension(0:NRMAX) :: X, R
    real(8), intent(out), dimension(0:NRMAX) :: dX
    integer :: NR
    real(8) :: DERIV3

    DO NR = 0, NRMAX
       dX(NR) = DERIV3(NR,R,X,NRMAX,NRM+1,0)
    END DO

  end subroutine DERIV3SB

!***************************************************************
!
!   Integral Method
!
!***************************************************************

  REAL(8) FUNCTION INTG_F(X)

    ! Calculate \int_0^b (r * X) dr

    use commons, only : NRMAX, NEMAX
    use core_module, only : fem_integral
    real(8), dimension(0:NRMAX), intent(in) :: X
    integer :: NE
    real(8) :: SUML

    SUML = 0.D0
    DO NE = 1, NEMAX
       SUML = SUML + SUM(fem_integral(15,NE,X))
    END DO
    INTG_F = SUML
    
  END FUNCTION INTG_F

  REAL(8) FUNCTION INTG_P(X,NR)

    ! Calculate \int (r * X) dr at one mesh

    use commons, only : NRMAX
    use core_module, only : fem_integral
    integer, intent(in) :: NR
    real(8), dimension(0:NRMAX), intent(in) :: X
    integer :: NE

    IF(NR == 0) THEN
       INTG_P = 0.D0
    ELSE
       NE = NR
       INTG_P = SUM(fem_integral(15,NE,X))
    END IF
    
  END FUNCTION INTG_P

  SUBROUTINE VALINT_SUB(X,NRLMAX,VAL)

    ! Calculate \int_0^r (r * X) dr

    use commons, only : NRMAX
    use core_module, only : fem_integral
    integer, intent(in) :: NRLMAX
    real(8), dimension(0:NRMAX), intent(in) :: X
    real(8), intent(out) :: VAL
    integer :: NE, NEMAX
    real(8) :: SUML

    NEMAX = NRLMAX
    SUML = 0.D0
    DO NE = 1, NEMAX
       SUML = SUML + SUM(fem_integral(15,NE,X))
    END DO
    VAL = SUML
    
  END SUBROUTINE VALINT_SUB

!***************************************************************
!
!   Coefficient function of Sauter model
!
!***************************************************************

  REAL(8) FUNCTION F33(X,Z)

    real(8), intent(in) :: X, Z

    F33 = 1.D0-(1.D0+0.36D0/Z)*X+0.59D0/Z*X**2-0.23D0/Z*X**3

  END FUNCTION F33

!***************************************************************
!
!   Mesh generating function
!
!***************************************************************

! This function is the function made by integration of the Lorentz function
!   R  : radial coordinate
!   C  : amplitude factor
!   W  : width of flat region around R=RC
!   RC : canter radial point of fine mesh region

  REAL(8) FUNCTION LORENTZ(R,C,W,RC,AMP)
    real(8), intent(in) :: r, c, w, rc
    real(8), intent(in), optional :: AMP

    LORENTZ = R + C * (W * ATAN((R - RC) / W) + W * ATAN(RC / W))
    if(present(amp)) LORENTZ = LORENTZ / AMP

  END FUNCTION LORENTZ

!***************************************************************
!
!   Bisection method for solving the equation
!
!***************************************************************

  SUBROUTINE BISECTION(f,cl,w,rc,amp,s,val)
    real(8), external :: f
    real(8), intent(in) :: cl, w, rc, amp, s
    real(8), intent(out) :: val
    integer :: i, n
    real(8) :: a, b, c, eps, fa, fc

    a = 0.d0 ; b= 1.d0
    eps = 1.d-5
    n = log10((b - a) / eps) / log10(2.d0) + 0.5d0
    fa = f(a,cl,w,rc,amp) - s
    do i = 1, n
       c = 0.5d0 * (a + b)
       fc = f(c,cl,w,rc,amp) - s
       if(fa * fc < 0.d0) then
          b = c
       else
          a = c
       end if
    end do
    val = c

  END SUBROUTINE BISECTION

end module libraries

