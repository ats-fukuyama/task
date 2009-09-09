!     $Id$
module tx_core_module
  use tx_commons, only : nrmax, h, r, psi, hpsi, nemax
  implicit none
  real(8) :: c13 = 1.d0 / 3.d0, c16 = 1.d0 / 6.d0, c112 = 1.d0 / 12.d0, c160 = 1.d0 / 60.d0
  public

contains

  function fem_int(id,a,b,c) result(x)

!-------------------------------------------------------
!
!   Calculate "\int_0^{psi_b} function(psi) dpsi"
!      dpsi : mesh interval
!      a    : coefficient vector
!      b    : coefficient vector
!      u    : variable vector
!      w    : weighting vector
!
!   function(psi) is classified as follows:
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
!      id = 15 : psi * a * u * w
!      id = 16 : psi * a'* u * w
!      id = 17 : psi * a * u'* w
!      id = 18 : psi * a * u'* w'
!      id = 19 : psi * a * u * w'
!      id = 20 : psi * a * b'* u * w'
!
!      id = 21 : sqrt(psi) * a * w (r*a*w)
!      id = 22 : r * a * u * w
!      id = 23 : r * a * u'* w
!      id = 24 : r * a * u * w'
!      id = 25 : r * a * b'* u * w
!      id = 26 : r * psi * a * u'* w'
!      id = 27 : r * psi * a * b'* u * w'
!
!      id = 28 : a * b * u * w
!      id = 29 : a * b'* u * w
!      id = 30 : a * b * u'* w
!      id = 31 : a * b * u * w'
!      id = 32 :(a * b * u)'* w
!      id = 33 : a'* b'* u * w
!      id = 34 : a * b'* u * w'
!      id = 35 : a * b * u'* w'
!
!      id = 36 : psi * a'* b * u * w
!      id = 37 : psi * a * b'* u * w
!      id = 38 : psi * a * b * u'* w
!      id = 39 :(psi * a * b * u)'* w
!      id = 40 : psi * a'* b * u * w'
!      id = 41 : psi * a * b * u'* w'
!      id = 42 :(psi * a * b * u)'* w'
!      id = 43 : psi * a * b * (u / b)'* w'
!
!      id = 44 : psi * a * b * u * w
!      id = 45 :(psi * a * b'* u)'* w
!      id = 46 :(psi * a * b'* u)'* w'
!      id = 47 : psi * a'* b'* u * w
!      id = 48 : psi * a * b'* u'* w
!
!      id = 49 : psi * a * b'* c'* u * w
!
!      id = -1 : a * w
!      id = -2 : a * b * w
!      id = -3 :(a * b)'* w
!      id = -8 : a * w'
!      id = -9 : a * b * w'
!
!   where ' means the derivative of psi
!
!  < input >
!     id       : mode select
!     a(nrmax) : coefficient vector, optional
!     b(nrmax) : coefficient vector, optional
!  < output >
!     x(1:nemax,4) : matrix of integrated values
!
!-------------------------------------------------------

    integer(4), intent(in) :: id
    real(8), intent(in), dimension(0:nrmax), optional  :: a, b, c
    integer(4) :: ne
    real(8) :: x(1:nemax,1:4), csq15, csq25, a1, a2, b1, b2, c1, c2, p1, p2, hp
    
    select case(id)
    case(-1)
       do ne = 1, nemax
          x(ne,1) = hpsi(ne) * c13 * a(ne-1)
          x(ne,2) = hpsi(ne) * c16 * a(ne)
          x(ne,3) = hpsi(ne) * c16 * a(ne-1)
          x(ne,4) = hpsi(ne) * c13 * a(ne)
       end do
    case(-2)
       do ne = 1, nemax
          x(ne,1) = ( 3.d0 * a(ne-1) +        a(ne)) * hpsi(ne) * c112 * b(ne-1)
          x(ne,2) = (        a(ne-1) +        a(ne)) * hpsi(ne) * c112 * b(ne)
          x(ne,3) = (        a(ne-1) +        a(ne)) * hpsi(ne) * c112 * b(ne-1)
          x(ne,4) = (        a(ne-1) + 3.d0 * a(ne)) * hpsi(ne) * c112 * b(ne)
       end do
    case(-3)
       do ne = 1, nemax
          x(ne,1) = (-4.d0 * a(ne-1) +        a(ne)) * c16 * b(ne-1)
          x(ne,2) = (        a(ne-1) + 2.d0 * a(ne)) * c16 * b(ne)
          x(ne,3) = (-2.d0 * a(ne-1) -        a(ne)) * c16 * b(ne-1)
          x(ne,4) = (-       a(ne-1) + 4.d0 * a(ne)) * c16 * b(ne)
       end do
    case(-8)
       x(1:nemax,1) =-0.5d0 * a(ne-1)
       x(1:nemax,2) =-0.5d0 * a(ne)
       x(1:nemax,3) = 0.5d0 * a(ne-1)
       x(1:nemax,4) = 0.5d0 * a(ne)
    case(0)
       !  for SUPG
       csq15 = 1.d0 / sqrt(15.d0) ! Raymond and Garder
       csq25 = 2.d0 * csq15 ! NASA report
       x(1:nemax,1) = hpsi(1:nemax) * csq25
       x(1:nemax,2) = x(1:nemax,1)
       x(1:nemax,3) = x(1:nemax,1)
       x(1:nemax,4) = x(1:nemax,1)
    case(1)
       x(1:nemax,1) = hpsi(1:nemax) * c13
       x(1:nemax,2) = hpsi(1:nemax) * c16
       x(1:nemax,3) = x(1:nemax,2)
       x(1:nemax,4) = x(1:nemax,1)
    case(2)
       do ne = 1, nemax
          x(ne,1) = ( 3.d0 * a(ne-1) +        a(ne)) * hpsi(ne) * c112
          x(ne,2) = (        a(ne-1) +        a(ne)) * hpsi(ne) * c112
          x(ne,3) = x(ne,2)
          x(ne,4) = (        a(ne-1) + 3.d0 * a(ne)) * hpsi(ne) * c112
       end do
    case(3)
       do ne = 1, nemax
          x(ne,1) = (-4.d0 * a(ne-1) +        a(ne)) * c16
          x(ne,2) = (        a(ne-1) + 2.d0 * a(ne)) * c16
          x(ne,3) = (-2.d0 * a(ne-1) -        a(ne)) * c16
          x(ne,4) = (-       a(ne-1) + 4.d0 * a(ne)) * c16
       end do
    case(4)
       x(1:nemax,1) = -0.5d0
       x(1:nemax,2) =  0.5d0
       x(1:nemax,3) = x(1:nemax,1)
       x(1:nemax,4) = x(1:nemax,2)
    case(5)
       do ne = 1, nemax
          x(ne,1) = (-2.d0 * a(ne-1) -        a(ne)) / 6.d0
          x(ne,2) = ( 2.d0 * a(ne-1) +        a(ne)) / 6.d0
          x(ne,3) = (-       a(ne-1) - 2.d0 * a(ne)) / 6.d0
          x(ne,4) = (        a(ne-1) + 2.d0 * a(ne)) / 6.d0
       end do
    case(6)
       do ne = 1, nemax
          x(ne,1) = (- a(ne-1) + a(ne)) * c13
          x(ne,2) = (- a(ne-1) + a(ne)) * c16
          x(ne,3) = x(ne,2)
          x(ne,4) = x(ne,1)
       end do
    case(8)
       x(1:nemax,1) =-0.5d0
       x(1:nemax,2) = x(1:nemax,1)
       x(1:nemax,3) = 0.5d0
       x(1:nemax,4) = x(1:nemax,3)
    case(9)
       do ne = 1, nemax
          x(ne,1) = (-2.d0 * a(ne-1) -        a(ne)) * c16
          x(ne,2) = (-       a(ne-1) - 2.d0 * a(ne)) * c16
          x(ne,3) =-x(ne,1)
          x(ne,4) =-x(ne,2)
       end do
    case(10)
       do ne = 1, nemax
          x(ne,1) = a(ne-1) / hpsi(ne)
          x(ne,2) =-a(ne)   / hpsi(ne)
          x(ne,3) =-x(ne,1)
          x(ne,4) =-x(ne,2)
       end do
    case(15)
       do ne = 1, nemax
          x(ne,1) = ( 12.d0*psi(ne-1)*a(ne-1) + 3.d0*psi(ne)*a(ne-1) &
               &     + 3.d0*psi(ne-1)*a(ne)   + 2.d0*psi(ne)*a(ne)) * hpsi(ne) * c160
          x(ne,2) = (  3.d0*psi(ne-1)*a(ne-1) + 2.d0*psi(ne)*a(ne-1) &
               &     + 2.d0*psi(ne-1)*a(ne)   + 3.d0*psi(ne)*a(ne)) * hpsi(ne) * c160
          x(ne,3) = x(ne,2)
          x(ne,4) = (  2.d0*psi(ne-1)*a(ne-1) + 3.d0*psi(ne)*a(ne-1) &
               &     + 3.d0*psi(ne-1)*a(ne)   +12.d0*psi(ne)*a(ne)) * hpsi(ne) * c160
       end do
    case(16)
       do ne = 1, nemax
          x(ne,1) = (3.d0 * psi(ne-1) +        psi(ne)) * (-a(ne-1) + a(ne)) * c112
          x(ne,2) = (       psi(ne-1) +        psi(ne)) * (-a(ne-1) + a(ne)) * c112
          x(ne,3) = x(ne,2)
          x(ne,4) = (       psi(ne-1) + 3.d0 * psi(ne)) * (-a(ne-1) + a(ne)) * c112
       end do
    case(17)
       do ne = 1, nemax
          x(ne,1) = (-3.d0*psi(ne-1)*a(ne-1) -      psi(ne)*a(ne-1) &
               &     -     psi(ne-1)*a(ne)   -      psi(ne)*a(ne)) * c112
          x(ne,2) =-x(ne,1)
          x(ne,3) = (-     psi(ne-1)*a(ne-1) -      psi(ne)*a(ne-1) &
               &     -     psi(ne-1)*a(ne)   - 3.d0*psi(ne)*a(ne)) * c112
          x(ne,4) =-x(ne,3)
       end do
    case(18)
       do ne = 1, nemax
          x(ne,1) = ( 2.d0*psi(ne-1)*a(ne-1) +      psi(ne)*a(ne-1) &
               &     +     psi(ne-1)*a(ne)   + 2.d0*psi(ne)*a(ne))  / hpsi(ne) * c16
          x(ne,2) =-x(ne,1)
          x(ne,3) = (-2.d0*psi(ne-1)*a(ne-1) -      psi(ne)*a(ne-1) &
               &     -     psi(ne-1)*a(ne)   - 2.d0*psi(ne)*a(ne))  / hpsi(ne) * c16
          x(ne,4) =-x(ne,3)
       end do
    case(20)
       do ne = 1, nemax
          x(ne,1) = ( 3.d0*psi(ne-1)*a(ne-1) +      psi(ne)*a(ne-1) &
               &     +     psi(ne-1)*a(ne)   +      psi(ne)*a(ne)) &
               &  * (b(ne-1) - b(ne)) / (12.d0 * hpsi(ne))
          x(ne,2) = (      psi(ne-1)*a(ne-1) +      psi(ne)*a(ne-1) &
               &     +     psi(ne-1)*a(ne)   + 3.d0*psi(ne)*a(ne)) &
               &  * (b(ne-1) - b(ne)) / (12.d0 * hpsi(ne))
          x(ne,3) =-x(ne,1)
          x(ne,4) =-x(ne,2)
       end do
    case(28)
       do ne = 1, nemax
          x(ne,1) = (12.d0*a(ne-1)*b(ne-1) + 3.d0*a(ne)*b(ne-1) &
               &    + 3.d0*a(ne-1)*b(ne)   + 2.d0*a(ne)*b(ne)) * hpsi(ne) * c160
          x(ne,2) = ( 3.d0*a(ne-1)*b(ne-1) + 2.d0*a(ne)*b(ne-1) &
               &    + 2.d0*a(ne-1)*b(ne)   + 3.d0*a(ne)*b(ne)) * hpsi(ne) * c160
          x(ne,3) = x(ne,2)
          x(ne,4) = ( 2.d0*a(ne-1)*b(ne-1) + 3.d0*a(ne)*b(ne-1) &
               &    + 3.d0*a(ne-1)*b(ne)   +12.d0*a(ne)*b(ne)) * hpsi(ne) * c160
       end do
    case(29)
       do ne = 1, nemax
          x(ne,1) = (-3.d0*a(ne-1)*b(ne-1) + 3.d0*a(ne-1)*b(ne) &
               &     -     a(ne)  *b(ne-1) +      a(ne)  *b(ne)) * c112
          x(ne,2) = (-     a(ne-1)*b(ne-1) +      a(ne-1)*b(ne) &
               &     -     a(ne)  *b(ne-1) +      a(ne)  *b(ne)) * c112
          x(ne,3) = (-     a(ne-1)*b(ne-1) +      a(ne-1)*b(ne) &
               &     -     a(ne)  *b(ne-1) +      a(ne)  *b(ne)) * c112
          x(ne,4) = (-     a(ne-1)*b(ne-1) +      a(ne-1)*b(ne) &
               &     -3.d0*a(ne)  *b(ne-1) + 3.d0*a(ne)  *b(ne)) * c112
       end do
    case(32)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne)
          b1 = b(ne-1) ; b2 = b(ne)
          x(ne,1) = (-3.d0*a1*b1 + 3.d0*a1*b2 -      a2*b1 +      a2*b2) * c112 &
               &  + (-3.d0*a1*b1 -      a1*b2 + 3.d0*a2*b1 +      a2*b2) * c112 &
               &  + (-3.d0*a1*b1 -      a1*b2 -      a2*b1 -      a2*b2) * c112
          x(ne,2) = (-     a1*b1 +      a1*b2 -      a2*b1 +      a2*b2) * c112 &
               &  + (-     a1*b1 -      a1*b2 +      a2*b1 +      a2*b2) * c112 &
               &  + ( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) * c112
          x(ne,3) = (-     a1*b1 +      a1*b2 -      a2*b1 +      a2*b2) * c112 &
               &  + (-     a1*b1 -      a1*b2 +      a2*b1 +      a2*b2) * c112 &
               &  + (-     a1*b1 -      a1*b2 -      a2*b1 - 3.d0*a2*b2) * c112
          x(ne,4) = (-     a1*b1 +      a1*b2 - 3.d0*a2*b1 + 3.d0*a2*b2) * c112 &
               &  + (-     a1*b1 - 3.d0*a1*b2 +      a2*b1 + 3.d0*a2*b2) * c112 &
               &  + (      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) * c112
       end do
    case(34)
       do ne = 1, nemax
          x(ne,1) = ( 2.d0*a(ne-1)*b(ne-1) - 2.d0*a(ne-1)*b(ne) &
               &     +     a(ne)  *b(ne-1) -      a(ne)  *b(ne)) / hpsi(ne) * c16
          x(ne,2) = (      a(ne-1)*b(ne-1) -      a(ne-1)*b(ne) &
               &     +2.d0*a(ne)  *b(ne-1) - 2.d0*a(ne)  *b(ne)) / hpsi(ne) * c16
          x(ne,3) =-x(ne,1)
          x(ne,4) =-x(ne,2)
       end do
    case(36)
       do ne = 1, nemax
          x(ne,1) = ( 12.d0*psi(ne-1)*b(ne-1) + 3.d0*psi(ne)*b(ne-1) &
               &     + 3.d0*psi(ne-1)*b(ne)   + 2.d0*psi(ne)*b(ne)) * (-a(ne-1)+a(ne)) * c160
          x(ne,2) = (  3.d0*psi(ne-1)*b(ne-1) + 2.d0*psi(ne)*b(ne-1) &
               &     + 2.d0*psi(ne-1)*b(ne)   + 3.d0*psi(ne)*b(ne)) * (-a(ne-1)+a(ne)) * c160
          x(ne,3) = x(ne,2)
          x(ne,4) = (  2.d0*psi(ne-1)*b(ne-1) + 3.d0*psi(ne)*b(ne-1) &
               &     + 3.d0*psi(ne-1)*b(ne)   +12.d0*psi(ne)*b(ne)) * (-a(ne-1)+a(ne)) * c160
       end do
    case(37)
       do ne = 1, nemax
          x(ne,1) = (12.d0*psi(ne-1)*a(ne-1) + 3.d0*psi(ne)*a(ne-1) &
               &    + 3.d0*psi(ne-1)*a(ne)   + 2.d0*psi(ne)*a(ne)) * (-b(ne-1)+b(ne)) * c160
          x(ne,2) = ( 3.d0*psi(ne-1)*a(ne-1) + 2.d0*psi(ne)*a(ne-1) &
               &    + 2.d0*psi(ne-1)*a(ne)   + 3.d0*psi(ne)*a(ne)) * (-b(ne-1)+b(ne)) * c160
          x(ne,3) = ( 3.d0*psi(ne-1)*a(ne-1) + 2.d0*psi(ne)*a(ne-1) &
               &    + 2.d0*psi(ne-1)*a(ne)   + 3.d0*psi(ne)*a(ne)) * (-b(ne-1)+b(ne)) * c160
          x(ne,4) = ( 2.d0*psi(ne-1)*a(ne-1) + 3.d0*psi(ne)*a(ne-1) &
               &    + 3.d0*psi(ne-1)*a(ne)   +12.d0*psi(ne)*a(ne)) * (-b(ne-1)+b(ne)) * c160
       end do
    case(38)
       do ne = 1, nemax
          a1 = a(ne-1) ; b1 = b(ne-1) ; p1 = psi(ne-1)
          a2 = a(ne)   ; b2 = b(ne)   ; p2 = psi(ne)
          x(ne,1) =-( 12.d0*p1*a1*b1 + 3.d0*p2*a1*b1 + 3.d0*p1*a2*b1 + 2.d0*p2*a2*b1 &
               &     + 3.d0*p1*a1*b2 + 2.d0*p2*a1*b2 + 2.d0*p1*a2*b2 + 3.d0*p2*a2*b2) * c160
          x(ne,2) = ( 12.d0*p1*a1*b1 + 3.d0*p2*a1*b1 + 3.d0*p1*a2*b1 + 2.d0*p2*a2*b1 &
               &     + 3.d0*p1*a1*b2 + 2.d0*p2*a1*b2 + 2.d0*p1*a2*b2 + 3.d0*p2*a2*b2) * c160
          x(ne,3) =-(  3.d0*p1*a1*b1 + 2.d0*p2*a1*b1 + 2.d0*p1*a2*b1 + 3.d0*p2*a2*b1 &
               &     + 2.d0*p1*a1*b2 + 3.d0*p2*a1*b2 + 3.d0*p1*a2*b2 +12.d0*p2*a2*b2) * c160
          x(ne,4) = (  3.d0*p1*a1*b1 + 2.d0*p2*a1*b1 + 2.d0*p1*a2*b1 + 3.d0*p2*a2*b1 &
               &     + 2.d0*p1*a1*b2 + 3.d0*p2*a1*b2 + 3.d0*p1*a2*b2 +12.d0*p2*a2*b2) * c160
       end do
    case(39)
       do ne = 1, nemax
          a1 = a(ne-1) ; b1 = b(ne-1) ; p1 = psi(ne-1)
          a2 = a(ne)   ; b2 = b(ne)   ; p2 = psi(ne)
          x(ne,1) = (-48.d0*a1*b1*p1+3.d0*a2*b1*p1+3.d0*a1*b2*p1+ 2.d0*a2*b2*p1 &
               &     + 3.d0*a1*b1*p2+2.d0*a2*b1*p2+2.d0*a1*b2*p2+ 3.d0*a2*b2*p2) * c160
          x(ne,2) = (  3.d0*a1*b1*p1+2.d0*a2*b1*p1+2.d0*a1*b2*p1+ 3.d0*a2*b2*p1 &
               &     + 2.d0*a1*b1*p2+3.d0*a2*b1*p2+3.d0*a1*b2*p2+12.d0*a2*b2*p2) * c160
          x(ne,3) = (-12.d0*a1*b1*p1-3.d0*a2*b1*p1-3.d0*a1*b2*p1- 2.d0*a2*b2*p1 &
               &     - 3.d0*a1*b1*p2-2.d0*a2*b1*p2-2.d0*a1*b2*p2- 3.d0*a2*b2*p2) * c160
          x(ne,4) = (- 3.d0*a1*b1*p1-2.d0*a2*b1*p1-2.d0*a1*b2*p1- 3.d0*a2*b2*p1 &
               &     - 2.d0*a1*b1*p2-3.d0*a2*b1*p2-3.d0*a1*b2*p2+48.d0*a2*b2*p2) * c160
       end do
    case(41)
       do ne = 1, nemax
          x(ne,1) = (b(ne-1)*( 3.d0*psi(ne-1)*a(ne-1)+     psi(ne)*a(ne-1) &
               &              +     psi(ne-1)*a(ne)  +     psi(ne)*a(ne)) &
               &    +b(ne)  *(      psi(ne-1)*a(ne-1)+     psi(ne)*a(ne-1) &
               &              +     psi(ne-1)*a(ne)  +3.d0*psi(ne)*a(ne))) / hpsi(ne) * c112
          x(ne,2) =-x(ne,1)
          x(ne,3) = x(ne,2)
          x(ne,4) = x(ne,1)
       end do
    case(42)
       do ne = 1, nemax
          x(ne,1) = a(ne-1)*b(ne-1)*psi(ne-1) / hpsi(ne)
          x(ne,2) =-a(ne)  *b(ne)  *psi(ne)   / hpsi(ne)
          x(ne,3) =-x(ne,1)
          x(ne,4) =-x(ne,2)
       end do
    case(44)
       do ne = 1, nemax
          x(ne,1) = ( 10.d0*a(ne-1)*b(ne-1)*psi(ne-1)+ 2.d0*a(ne)*b(ne-1)*psi(ne-1) &
               &     + 2.d0*a(ne-1)*b(ne)  *psi(ne-1)+      a(ne)*b(ne)  *psi(ne-1) &
               &     + 2.d0*a(ne-1)*b(ne-1)*psi(ne)  +      a(ne)*b(ne-1)*psi(ne) &
               &     +      a(ne-1)*b(ne)  *psi(ne)  +      a(ne)*b(ne)  *psi(ne)) * hpsi(ne) * c160
          x(ne,2) = (  2.d0*a(ne-1)*b(ne-1)*psi(ne-1)+      a(ne)*b(ne-1)*psi(ne-1) &
               &     +      a(ne-1)*b(ne)  *psi(ne-1)+      a(ne)*b(ne)  *psi(ne-1) &
               &     +      a(ne-1)*b(ne-1)*psi(ne)  +      a(ne)*b(ne-1)*psi(ne) &
               &     +      a(ne-1)*b(ne)  *psi(ne)  + 2.d0*a(ne)*b(ne)  *psi(ne)) * hpsi(ne) * c160
          x(ne,3) = x(ne,2)
          x(ne,4) = (       a(ne-1)*b(ne-1)*psi(ne-1)+      a(ne)*b(ne-1)*psi(ne-1) &
               &     +      a(ne-1)*b(ne)  *psi(ne-1)+ 2.d0*a(ne)*b(ne)  *psi(ne-1) &
               &     +      a(ne-1)*b(ne-1)*psi(ne)  + 2.d0*a(ne)*b(ne-1)*psi(ne) &
               &     + 2.d0*a(ne-1)*b(ne)  *psi(ne)  +10.d0*a(ne)*b(ne)  *psi(ne)) * hpsi(ne) * c160
       end do
    case(45)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne)
          p1 = psi(ne-1) ; p2 = psi(ne) ; hp = hpsi(ne)
          x(ne,1) = (b1-b2)*( 9.d0*a1*p1-a2*p1-a1*p2-     a2*p2)/hp * c112
          x(ne,2) =-(b1-b2)*(      a1*p1+a2*p1+a1*p2+3.d0*a2*p2)/hp * c112
          x(ne,3) = (b1-b2)*( 3.d0*a1*p1+a2*p1+a1*p2+     a2*p2)/hp * c112
          x(ne,4) =-(b1-b2)*(-     a1*p1-a2*p1-a1*p2+9.d0*a2*p2)/hp * c112
       end do
    case(48)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne)
          p1 = psi(ne-1) ; p2 = psi(ne) ; hp = hpsi(ne)
          x(ne,1) =-(3.d0*a1*p1+a2*p1+a1*p2+     a2*p2)*(b2-b1)/hp * c112
          x(ne,2) = (3.d0*a1*p1+a2*p1+a1*p2+     a2*p2)*(b2-b1)/hp * c112
          x(ne,3) =-(     a1*p1+a2*p1+a1*p2+3.d0*a2*p2)*(b2-b1)/hp * c112
          x(ne,4) = (     a1*p1+a2*p1+a1*p2+3.d0*a2*p2)*(b2-b1)/hp * c112
       end do
    case(49)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne) ; c1 = c(ne-1) ; c2 = c(ne)
          p1 = psi(ne-1) ; p2 = psi(ne) ; hp = hpsi(ne)
          x(ne,1) = (3.d0*p1*(4.d0*a1+a2)+p2*(3.d0*a1+2.d0*a2))*(b2-b1)*(c2-c1)/hp*c160
          x(ne,2) = (p1*(3.d0*a1+2.d0*a2)+p2*(2.d0*a1+3.d0*a2))*(b2-b1)*(c2-c1)/hp*c160
          x(ne,3) = x(ne,2)
          x(ne,4) = (p1*(2.d0*a1+3.d0*a2)+3.d0*p2*(a1+4.d0*a2))*(b2-b1)*(c2-c1)/hp*c160
       end do
    case default
       write(6,*)  'XX falut ID in fem_int, id= ',id
       stop
    end select

  end function fem_int

  function fem_int_point(id,ne,a,b,c) result(x)
!-------------------------------------------------------
!
!   Calculate "\int_{psi_i}^{psi_{i+1}} function(psi) dpsi"
!      dpsi : mesh interval
!      a    : coefficient vector
!      b    : coefficient vector
!      u    : variable vector
!      w    : weighting vector
!
!   function(psi) is classified as follows:
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
!      id = 15 : psi * a * u * w
!      id = 16 : psi * a'* u * w
!      id = 17 : psi * a * u'* w
!      id = 18 : psi * a * u'* w'
!      id = 19 : psi * a * u * w'
!      id = 20 : psi * a * b'* u * w'
!
!      id = 21 : sqrt(psi) * a * w (r*a*w)
!      id = 22 : r * a * u * w
!      id = 23 : r * a * u'* w
!      id = 24 : r * a * u * w'
!      id = 25 : r * a * b'* u * w
!      id = 26 : r * psi * a * u'* w'
!      id = 27 : r * psi * a * b'* u * w'
!
!      id = 28 : a * b * u * w
!      id = 29 : a * b'* u * w
!      id = 30 : a * b * u'* w
!      id = 31 : a * b * u * w'
!      id = 32 :(a * b * u)'* w
!      id = 33 : a'* b'* u * w
!      id = 34 : a * b'* u * w'
!      id = 35 : a * b * u'* w'
!
!      id = 36 : psi * a'* b * u * w
!      id = 37 : psi * a * b'* u * w
!      id = 38 : psi * a * b * u'* w
!      id = 39 :(psi * a * b * u)'* w
!      id = 40 : psi * a'* b * u * w'
!      id = 41 : psi * a * b * u'* w'
!      id = 42 :(psi * a * b * u)'* w'
!      id = 43 : psi * a * b * (u / b)'* w'
!
!      id = 44 : psi * a * b * u * w
!      id = 45 :(psi * a * b'* u)'* w
!      id = 46 :(psi * a * b'* u)'* w'
!      id = 47 : psi * a'* b'* u * w
!      id = 48 : psi * a * b'* u'* w
!
!      id = 49 : psi * a * b'* c'* u * w
!
!      id = -1 : a * w
!      id = -2 : a * b * w
!      id = -8 : a * w'
!      id = -9 : a * b * w'
!
!   where ' means the derivative of psi
!
!  < input >
!     id       : mode select
!     ne       : current number of elements
!     a(nrmax) : coefficient vector, optional
!     b(nrmax) : coefficient vector, optional
!  < output >
!     x(4)     : matrix of integrated values
!
!  If ne = 0 is given, the coefficient "a" at the magnetic axis is forced to be zero.
!
!-------------------------------------------------------
    integer(4), intent(in) :: id, ne
    real(8), intent(in), dimension(0:nrmax), optional  :: a, b, c
    integer(4) :: nel, node1, node2, iflag
    real(8) :: x(1:4), a1, a2, r1, r2, p1, p2, b1, b2, c1, c2, hp, csq15

    iflag = 0
    if(ne == 0) then
       nel = ne + 1
       iflag = 1
    else
       nel = ne
    end if

    node1 = nel-1  ; node2 = nel
    if(present(a)) then
       a1 = a(node1) ; a2 = a(node2)
       if(present(b)) then
          b1 = b(node1) ; b2 = b(node2)
       end if
    end if
    if(iflag == 1) a1 = 0.d0
    r1 = r(node1) ; r2 = r(node2)
    p1 = psi(node1) ; p2 = psi(node2) ; hp = hpsi(nel)

    select case(id)
    case(-1)
       x(1) = hp / 3.d0 * a1
       x(2) = hp / 6.d0 * a2
       x(3) = hp / 6.d0 * a1
       x(4) = hp / 3.d0 * a2
    case(-2)
       x(1) = ( 3.d0 * a1 +        a2) * hp / 12.d0 * b1
       x(2) = (        a1 +        a2) * hp / 12.d0 * b2
       x(3) = (        a1 +        a2) * hp / 12.d0 * b1
       x(4) = (        a1 + 3.d0 * a2) * hp / 12.d0 * b2
    case(-8)
       x(1) =-0.5d0 * a1
       x(2) =-0.5d0 * a2
       x(3) = 0.5d0 * a1
       x(4) = 0.5d0 * a2
    case(-9)
       x(1) = (-2.d0 * a1 -        a2) / 6.d0 * b1
       x(2) = (-       a1 - 2.d0 * a2) / 6.d0 * b2
       x(3) = ( 2.d0 * a1 +        a2) / 6.d0 * b1
       x(4) = (        a1 + 2.d0 * a2) / 6.d0 * b2
!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    case(0)
       csq15 = 1.d0 / sqrt(15.d0)
       x(1) = hp * csq15
       x(2) = x(1)
       x(3) = x(1)
       x(4) = x(1)
    case(1)
       x(1) = hp / 3.d0
       x(2) = hp / 6.d0
       x(3) = hp / 6.d0
       x(4) = hp / 3.d0
    case(2)
       x(1) = ( 3.d0 * a1 +        a2) * hp / 12.d0
       x(2) = (        a1 +        a2) * hp / 12.d0
       x(3) = (        a1 +        a2) * hp / 12.d0
       x(4) = (        a1 + 3.d0 * a2) * hp / 12.d0
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
       x(1) = (- a1 + a2) / 3.d0
       x(2) = (- a1 + a2) / 6.d0
       x(3) = (- a1 + a2) / 6.d0
       x(4) = (- a1 + a2) / 3.d0
    case(7)
       x(1) = ( a1 - a2) / (2.d0 * hp)
       x(2) = ( a1 - a2) / (2.d0 * hp)
       x(3) = (-a1 + a2) / (2.d0 * hp)
       x(4) = (-a1 + a2) / (2.d0 * hp)
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
       x(1) = a1 / hp
       x(2) =-a2 / hp
       x(3) =-a1 / hp
       x(4) = a2 / hp
    case(11)
       x(1) = 1.d0 / hp
       x(2) =-1.d0 / hp
       x(3) =-1.d0 / hp
       x(4) = 1.d0 / hp
    case(12)
       x(1) = ( a1 + a2) / (2.d0 * hp)
       x(2) = (-a1 - a2) / (2.d0 * hp)
       x(3) = (-a1 - a2) / (2.d0 * hp)
       x(4) = ( a1 + a2) / (2.d0 * hp)
    case(13)
       x(1) = ( a1 - a2) / (2.d0 * hp)
       x(2) = ( a1 - a2) / (2.d0 * hp)
       x(3) = (-a1 + a2) / (2.d0 * hp)
       x(4) = (-a1 + a2) / (2.d0 * hp)
    case(14)
       x(1) = (-a1 + a2) / hp**2
       x(2) = ( a1 - a2) / hp**2
       x(3) = ( a1 - a2) / hp**2
       x(4) = (-a1 + a2) / hp**2
    case(15)
       x(1) = (12.d0*p1*a1 + 3.d0*p2*a1 + 3.d0*p1*a2 + 2.d0*p2*a2) * hp / 60.d0
       x(2) = ( 3.d0*p1*a1 + 2.d0*p2*a1 + 2.d0*p1*a2 + 3.d0*p2*a2) * hp / 60.d0
       x(3) = ( 3.d0*p1*a1 + 2.d0*p2*a1 + 2.d0*p1*a2 + 3.d0*p2*a2) * hp / 60.d0
       x(4) = ( 2.d0*p1*a1 + 3.d0*p2*a1 + 3.d0*p1*a2 +12.d0*p2*a2) * hp / 60.d0
    case(16)
       x(1) = (3.d0 * p1 +        p2) * (-a1 + a2) / 12.d0
       x(2) = (       p1 +        p2) * (-a1 + a2) / 12.d0
       x(3) = (       p1 +        p2) * (-a1 + a2) / 12.d0
       x(4) = (       p1 + 3.d0 * p2) * (-a1 + a2) / 12.d0
    case(17)
       x(1) = (-3.d0*p1*a1 - p2*a1 - p1*a2 -      p2*a2) / 12.d0
       x(2) = ( 3.d0*p1*a1 + p2*a1 + p1*a2 +      p2*a2) / 12.d0
       x(3) = (-     p1*a1 - p2*a1 - p1*a2 - 3.d0*p2*a2) / 12.d0
       x(4) = (      p1*a1 + p2*a1 + p1*a2 + 3.d0*p2*a2) / 12.d0
    case(18)
       x(1) = ( 2.d0*p1*a1 + p2*a1 + p1*a2 + 2.d0*p2*a2) / (6.d0 * hp)
       x(2) = (-2.d0*p1*a1 - p2*a1 - p1*a2 - 2.d0*p2*a2) / (6.d0 * hp)
       x(3) = (-2.d0*p1*a1 - p2*a1 - p1*a2 - 2.d0*p2*a2) / (6.d0 * hp)
       x(4) = ( 2.d0*p1*a1 + p2*a1 + p1*a2 + 2.d0*p2*a2) / (6.d0 * hp)
    case(19)
       x(1) = (-3.d0*p1*a1 - p2*a1 - p1*a2 -      p2*a2) / 12.d0
       x(2) = (-     p1*a1 - p2*a1 - p1*a2 - 3.d0*p2*a2) / 12.d0
       x(3) = ( 3.d0*p1*a1 + p2*a1 + p1*a2 +      p2*a2) / 12.d0
       x(4) = (      p1*a1 + p2*a1 + p1*a2 + 3.d0*p2*a2) / 12.d0
    case(20)
       x(1) = (3.d0*p1*a1 + p2*a1 + p1*a2 +      p2*a2) * (b1 - b2) / (12.d0 * hp)
       x(2) = (     p1*a1 + p2*a1 + p1*a2 + 3.d0*p2*a2) * (b1 - b2) / (12.d0 * hp)
       x(3) =-(3.d0*p1*a1 + p2*a1 + p1*a2 +      p2*a2) * (b1 - b2) / (12.d0 * hp)
       x(4) =-(     p1*a1 + p2*a1 + p1*a2 + 3.d0*p2*a2) * (b1 - b2) / (12.d0 * hp)
    case(21)
!       x(1) = (3.d0 * r1 +        r2) * hp / 12.d0 * a1
!       x(2) = (       r1 +        r2) * hp / 12.d0 * a2
!       x(3) = (       r1 +        r2) * hp / 12.d0 * a1
!       x(4) = (       r1 + 3.d0 * r2) * hp / 12.d0 * a2
       x(1) = 2.d0*h(nel)*( 15.d0*p1**2+45.d0*p1*r1*r2+48.d0*p1*p2+24.d0*r1*p2*r2 &
            &             + 8.d0*p2**2) / (105.d0*(r1+r2)**2) * a1
       x(2) = 4.d0*h(nel)*(  3.d0*p1**2+ 9.d0*p1*r1*r2+11.d0*p1*p2+ 9.d0*r1*p2*r2 &
            &             + 3.d0*p2**2) / (105.d0*(r1+r2)**2) * a2
       x(3) = 4.d0*h(nel)*(  3.d0*p1**2+ 9.d0*p1*r1*r2+11.d0*p1*p2+ 9.d0*r1*p2*r2 &
            &             + 3.d0*p2**2) / (105.d0*(r1+r2)**2) * a1
       x(4) = 2.d0*h(nel)*(  8.d0*p1**2+24.d0*p1*r1*r2+48.d0*p1*p2+45.d0*r1*p2*r2 &
            &             +15.d0*p2**2) / (105.d0*(r1+r2)**2) * a2
    case(22)
       x(1) = (12.d0*r1*a1 + 3.d0*r2*a1 + 3.d0*r1*a2 + 2.d0*r2*a2) * hp / 60.d0
       x(2) = ( 3.d0*r1*a1 + 2.d0*r2*a1 + 2.d0*r1*a2 + 3.d0*r2*a2) * hp / 60.d0
       x(3) = ( 3.d0*r1*a1 + 2.d0*r2*a1 + 2.d0*r1*a2 + 3.d0*r2*a2) * hp / 60.d0
       x(4) = ( 2.d0*r1*a1 + 3.d0*r2*a1 + 3.d0*r1*a2 +12.d0*r2*a2) * hp / 60.d0
    case(23)
       x(1) = (-3.d0*r1*a1 - r2*a1 - r1*a2 -      r2*a2) / 12.d0
       x(2) = ( 3.d0*r1*a1 + r2*a1 + r1*a2 +      r2*a2) / 12.d0
       x(3) = (-     r1*a1 - r2*a1 - r1*a2 - 3.d0*r2*a2) / 12.d0
       x(4) = (      r1*a1 + r2*a1 + r1*a2 + 3.d0*r2*a2) / 12.d0
    case(24)
       x(1) = (-3.d0*r1*a1 - r2*a1 - r1*a2 -      r2*a2) / 12.d0
       x(2) = (-     r1*a1 - r2*a1 - r1*a2 - 3.d0*r2*a2) / 12.d0
       x(3) = ( 3.d0*r1*a1 + r2*a1 + r1*a2 +      r2*a2) / 12.d0
       x(4) = (      r1*a1 + r2*a1 + r1*a2 + 3.d0*r2*a2) / 12.d0
    case(25)
       x(1) =  ((3.d0*r1+r2)*a1+(r1+     r2)*a2)*(b1-b2) / (12.d0 * hp)
       x(2) =  ((     r1+r2)*a1+(r1+3.d0*r2)*a2)*(b1-b2) / (12.d0 * hp)
       x(3) = -((3.d0*r1+r2)*a1+(r1+     r2)*a2)*(b1-b2) / (12.d0 * hp)
       x(4) = -((     r1+r2)*a1+(r1+3.d0*r2)*a2)*(b1-b2) / (12.d0 * hp)
    case(26)
       x(1) =  ((p1*r2+p2*r2+3.d0*p1*r1+p2*r1)*a1+(p1*r2+3.d0*p2*r2+p1*r1+p2*r1)*a2) &
            & / (12.d0 * hp)
       x(2) = -((p1*r2+p2*r2+3.d0*p1*r1+p2*r1)*a1+(p1*r2+3.d0*p2*r2+p1*r1+p2*r1)*a2) &
            & / (12.d0 * hp)
       x(3) = -((p1*r2+p2*r2+3.d0*p1*r1+p2*r1)*a1+(p1*r2+3.d0*p2*r2+p1*r1+p2*r1)*a2) &
            & / (12.d0 * hp)
       x(4) =  ((p1*r2+p2*r2+3.d0*p1*r1+p2*r1)*a1+(p1*r2+3.d0*p2*r2+p1*r1+p2*r1)*a2) &
            & / (12.d0 * hp)
    case(27)
       x(1) = -( (10.d0*p1*r1+2.d0*p2*r1+2.d0*p1*r2+      p2*r2)*a1 &
            &   +( 2.d0*p1*r1+     p2*r1+     p1*r2+      p2*r2)*a2)*(b1-b2) / 60.d0
       x(2) = -( ( 2.d0*p1*r1+     p2*r1+     p1*r2+      p2*r2)*a1 &
            &   +(      p1*r1+     p2*r1+     p1*r2+ 2.d0*p2*r2)*a2)*(b1-b2) / 60.d0
       x(3) = -( ( 2.d0*p1*r1+     p2*r1+     p1*r2+      p2*r2)*a1 &
            &   +(      p1*r1+     p2*r1+     p1*r2+ 2.d0*p2*r2)*a2)*(b1-b2) / 60.d0
       x(4) = -( (      p1*r1+     p2*r1+     p1*r2+ 2.d0*p2*r2)*a1 &
            &   +(      p1*r1+2.d0*p2*r1+2.d0*p1*r2+10.d0*p2*r2)*a2)*(b1-b2) / 60.d0
    case(28)
       x(1) = (12.d0*a1*b1 + 3.d0*a2*b1 + 3.d0*a1*b2 + 2.d0*a2*b2) * hp / 60.d0
       x(2) = ( 3.d0*a1*b1 + 2.d0*a2*b1 + 2.d0*a1*b2 + 3.d0*a2*b2) * hp / 60.d0
       x(3) = ( 3.d0*a1*b1 + 2.d0*a2*b1 + 2.d0*a1*b2 + 3.d0*a2*b2) * hp / 60.d0
       x(4) = ( 2.d0*a1*b1 + 3.d0*a2*b1 + 3.d0*a1*b2 +12.d0*a2*b2) * hp / 60.d0
    case(29)
       x(1) = (-3.d0*a1*b1 + 3.d0*a1*b2 -      a2*b1 +      a2*b2) / 12.d0
       x(2) = (-     a1*b1 +      a1*b2 -      a2*b1 +      a2*b2) / 12.d0
       x(3) = (-     a1*b1 +      a1*b2 -      a2*b1 +      a2*b2) / 12.d0
       x(4) = (-     a1*b1 +      a1*b2 - 3.d0*a2*b1 + 3.d0*a2*b2) / 12.d0
    case(30)
       x(1) = (-3.d0*a1*b1 -      a1*b2 -      a2*b1 -      a2*b2) / 12.d0
       x(2) = ( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) / 12.d0
       x(3) = (-     a1*b1 -      a1*b2 -      a2*b1 - 3.d0*a2*b2) / 12.d0
       x(4) = (      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) / 12.d0
    case(31)
       x(1) = (-3.d0*a1*b1 -      a1*b2 -      a2*b1 -      a2*b2) / 12.d0
       x(2) = (-     a1*b1 -      a1*b2 -      a2*b1 - 3.d0*a2*b2) / 12.d0
       x(3) = ( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) / 12.d0
       x(4) = (      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) / 12.d0
    case(32)
       x(1) = (-3.d0*a1*b1 + 3.d0*a1*b2 -      a2*b1 +      a2*b2) / 12.d0 &
            &+(-3.d0*a1*b1 -      a1*b2 + 3.d0*a2*b1 +      a2*b2) / 12.d0 &
            &+(-3.d0*a1*b1 -      a1*b2 -      a2*b1 -      a2*b2) / 12.d0
       x(2) = (-     a1*b1 +      a1*b2 -      a2*b1 +      a2*b2) / 12.d0 &
            &+(-     a1*b1 -      a1*b2 +      a2*b1 +      a2*b2) / 12.d0 &
            &+( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) / 12.d0
       x(3) = (-     a1*b1 +      a1*b2 -      a2*b1 +      a2*b2) / 12.d0 &
            &+(-     a1*b1 -      a1*b2 +      a2*b1 +      a2*b2) / 12.d0 &
            &+(-     a1*b1 -      a1*b2 -      a2*b1 - 3.d0*a2*b2) / 12.d0
       x(4) = (-     a1*b1 +      a1*b2 - 3.d0*a2*b1 + 3.d0*a2*b2) / 12.d0 &
            &+(-     a1*b1 - 3.d0*a1*b2 +      a2*b1 + 3.d0*a2*b2) / 12.d0 &
            &+(      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) / 12.d0
    case(33)
       x(1) = (-3.d0*a1*b1 -      a1*b2 + 3.d0*a2*b1 +      a2*b2) / 12.d0
       x(2) = (-     a1*b1 -      a1*b2 +      a2*b1 +      a2*b2) / 12.d0
       x(3) = (-     a1*b1 -      a1*b2 +      a2*b1 +      a2*b2) / 12.d0
       x(4) = (-     a1*b1 - 3.d0*a1*b2 +      a2*b1 + 3.d0*a2*b2) / 12.d0
    case(34)
       x(1) = ( 2.d0*a1*b1 - 2.d0*a1*b2 +      a2*b1 -      a2*b2) / (6.d0 * hp)
       x(2) = (      a1*b1 -      a1*b2 + 2.d0*a2*b1 - 2.d0*a2*b2) / (6.d0 * hp)
       x(3) = (-2.d0*a1*b1 + 2.d0*a1*b2 -      a2*b1 +      a2*b2) / (6.d0 * hp)
       x(4) = (-     a1*b1 +      a1*b2 - 2.d0*a2*b1 + 2.d0*a2*b2) / (6.d0 * hp)
    case(35)
       x(1) = (2.d0*a1*b1 + a2*b1 + a1*b2 + 2.d0*a2*b2) / (6.d0 * hp)
       x(2) =-(2.d0*a1*b1 + a2*b1 + a1*b2 + 2.d0*a2*b2) / (6.d0 * hp)
       x(3) =-(2.d0*a1*b1 + a2*b1 + a1*b2 + 2.d0*a2*b2) / (6.d0 * hp)
       x(4) = (2.d0*a1*b1 + a2*b1 + a1*b2 + 2.d0*a2*b2) / (6.d0 * hp)
    case(36)
       x(1) = (12.d0*p1*b1 + 3.d0*p2*b1 + 3.d0*p1*b2 + 2.d0*p2*b2) * (-a1+a2) / 60.d0
       x(2) = ( 3.d0*p1*b1 + 2.d0*p2*b1 + 2.d0*p1*b2 + 3.d0*p2*b2) * (-a1+a2) / 60.d0
       x(3) = ( 3.d0*p1*b1 + 2.d0*p2*b1 + 2.d0*p1*b2 + 3.d0*p2*b2) * (-a1+a2) / 60.d0
       x(4) = ( 2.d0*p1*b1 + 3.d0*p2*b1 + 3.d0*p1*b2 +12.d0*p2*b2) * (-a1+a2) / 60.d0
    case(37)
       x(1) = (12.d0*p1*a1 + 3.d0*p2*a1 + 3.d0*p1*a2 + 2.d0*p2*a2) * (-b1+b2) / 60.d0
       x(2) = ( 3.d0*p1*a1 + 2.d0*p2*a1 + 2.d0*p1*a2 + 3.d0*p2*a2) * (-b1+b2) / 60.d0
       x(3) = ( 3.d0*p1*a1 + 2.d0*p2*a1 + 2.d0*p1*a2 + 3.d0*p2*a2) * (-b1+b2) / 60.d0
       x(4) = ( 2.d0*p1*a1 + 3.d0*p2*a1 + 3.d0*p1*a2 +12.d0*p2*a2) * (-b1+b2) / 60.d0
    case(38)
       x(1) =-( 12.d0*p1*a1*b1 + 3.d0*p2*a1*b1 + 3.d0*p1*a2*b1 + 2.d0*p2*a2*b1 &
            &  + 3.d0*p1*a1*b2 + 2.d0*p2*a1*b2 + 2.d0*p1*a2*b2 + 3.d0*p2*a2*b2) / 60.d0
       x(2) = ( 12.d0*p1*a1*b1 + 3.d0*p2*a1*b1 + 3.d0*p1*a2*b1 + 2.d0*p2*a2*b1 &
            &  + 3.d0*p1*a1*b2 + 2.d0*p2*a1*b2 + 2.d0*p1*a2*b2 + 3.d0*p2*a2*b2) / 60.d0
       x(3) =-(  3.d0*p1*a1*b1 + 2.d0*p2*a1*b1 + 2.d0*p1*a2*b1 + 3.d0*p2*a2*b1 &
            &  + 2.d0*p1*a1*b2 + 3.d0*p2*a1*b2 + 3.d0*p1*a2*b2 +12.d0*p2*a2*b2) / 60.d0
       x(4) = (  3.d0*p1*a1*b1 + 2.d0*p2*a1*b1 + 2.d0*p1*a2*b1 + 3.d0*p2*a2*b1 &
            &  + 2.d0*p1*a1*b2 + 3.d0*p2*a1*b2 + 3.d0*p1*a2*b2 +12.d0*p2*a2*b2) / 60.d0
    case(39)
       x(1) = (-48.d0*a1*b1*p1+3.d0*a2*b1*p1+3.d0*a1*b2*p1+ 2.d0*a2*b2*p1 &
            &  + 3.d0*a1*b1*p2+2.d0*a2*b1*p2+2.d0*a1*b2*p2+ 3.d0*a2*b2*p2) / 60.d0
       x(2) = (  3.d0*a1*b1*p1+2.d0*a2*b1*p1+2.d0*a1*b2*p1+ 3.d0*a2*b2*p1 &
            &  + 2.d0*a1*b1*p2+3.d0*a2*b1*p2+3.d0*a1*b2*p2+12.d0*a2*b2*p2) / 60.d0
       x(3) = (-12.d0*a1*b1*p1-3.d0*a2*b1*p1-3.d0*a1*b2*p1- 2.d0*a2*b2*p1 &
            &  - 3.d0*a1*b1*p2-2.d0*a2*b1*p2-2.d0*a1*b2*p2- 3.d0*a2*b2*p2) / 60.d0
       x(4) = (- 3.d0*a1*b1*p1-2.d0*a2*b1*p1-2.d0*a1*b2*p1- 3.d0*a2*b2*p1 &
            &  - 2.d0*a1*b1*p2-3.d0*a2*b1*p2-3.d0*a1*b2*p2+48.d0*a2*b2*p2) / 60.d0
    case(40)
       x(1) = (3.d0*p1*b1+p2*b1+p1*b2+     p2*b2) * (a1-a2) / (12.d0 * hp)
       x(2) =-(3.d0*p1*b1+p2*b1+p1*b2+     p2*b2) * (a1-a2) / (12.d0 * hp)
       x(3) = (     p1*b1+p2*b1+p1*b2+3.d0*p2*b2) * (a1-a2) / (12.d0 * hp)
       x(4) =-(     p1*b1+p2*b1+p1*b2+3.d0*p2*b2) * (a1-a2) / (12.d0 * hp)
    case(41)
       x(1) = (b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
          & / (12.d0 * hp)
       x(2) =-(b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
          & / (12.d0 * hp)
       x(3) =-(b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
          & / (12.d0 * hp)
       x(4) = (b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
          & / (12.d0 * hp)
    case(42)
       x(1) = a1*b1*p1 / hp
       x(2) =-a2*b2*p2 / hp
       x(3) =-a1*b1*p1 / hp
       x(4) = a2*b2*p2 / hp
    case(43)
       x(1) = (b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
          & / (12.d0 * b1 * hp)
       x(2) =-(b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
          & / (12.d0 * b2 * hp)
       x(3) =-(b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
          & / (12.d0 * b1 * hp)
       x(4) = (b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
          & / (12.d0 * b2 * hp)
    case(44)
       x(1) = ( 10.d0*a1*b1*p1+2.d0*a2*b1*p1+2.d0*a1*b2*p1+      a2*b2*p1 &
            &  + 2.d0*a1*b1*p2+     a2*b1*p2+     a1*b2*p2+      a2*b2*p2) * hp / 60.d0
       x(2) = (  2.d0*a1*b1*p1+     a2*b1*p1+     a1*b2*p1+      a2*b2*p1 &
            &  +      a1*b1*p2+     a2*b1*p2+     a1*b2*p2+ 2.d0*a2*b2*p2) * hp / 60.d0
       x(3) = (  2.d0*a1*b1*p1+     a2*b1*p1+     a1*b2*p1+      a2*b2*p1 &
            &  +      a1*b1*p2+     a2*b1*p2+     a1*b2*p2+ 2.d0*a2*b2*p2) * hp / 60.d0
       x(4) = (       a1*b1*p1+     a2*b1*p1+     a1*b2*p1+ 2.d0*a2*b2*p1 &
            &  +      a1*b1*p2+2.d0*a2*b1*p2+2.d0*a1*b2*p2+10.d0*a2*b2*p2) * hp / 60.d0
    case(45)
       x(1) = (-9.d0*a1*p1+a2*p1+a1*p2+     a2*p2)*(b2-b1)/(12.d0*hp)
       x(2) = (      a1*p1+a2*p1+a1*p2+3.d0*a2*p2)*(b2-b1)/(12.d0*hp)
       x(3) = (-3.d0*a1*p1-a2*p1-a1*p2-     a2*p2)*(b2-b1)/(12.d0*hp)
       x(4) = (-     a1*p1-a2*p1-a1*p2+9.d0*a2*p2)*(b2-b1)/(12.d0*hp)
    case(46)
       x(1) =-a1*(b1-b2)*p1/hp**2
       x(2) = a2*(b1-b2)*p2/hp**2
       x(3) = a1*(b1-b2)*p1/hp**2
       x(4) =-a2*(b1-b2)*p2/hp**2
    case(47)
       x(1) = (3.d0*p1+     p2)*(a2-a1)*(b2-b1)/(12.d0*hp)
       x(2) = (     p1+     p2)*(a2-a1)*(b2-b1)/(12.d0*hp)
       x(3) = (     p1+     p2)*(a2-a1)*(b2-b1)/(12.d0*hp)
       x(4) = (     p1+3.d0*p2)*(a2-a1)*(b2-b1)/(12.d0*hp)
    case(48)
       x(1) =-(3.d0*a1*p1+a2*p1+a1*p2+     a2*p2)*(b2-b1)/(12.d0*hp)
       x(2) = (3.d0*a1*p1+a2*p1+a1*p2+     a2*p2)*(b2-b1)/(12.d0*hp)
       x(3) =-(     a1*p1+a2*p1+a1*p2+3.d0*a2*p2)*(b2-b1)/(12.d0*hp)
       x(4) = (     a1*p1+a2*p1+a1*p2+3.d0*a2*p2)*(b2-b1)/(12.d0*hp)
    case(49)
       x(1) = (3.d0*p1*(4.d0*a1+a2)+p2*(3.d0*a1+2.d0*a2))*(b2-b1)*(c2-c1)/(60.d0*hp)
       x(2) = (p1*(3.d0*a1+2.d0*a2)+p2*(2.d0*a1+3.d0*a2))*(b2-b1)*(c2-c1)/(60.d0*hp)
       x(3) = x(2)
       x(4) = (p1*(2.d0*a1+3.d0*a2)+3.d0*p2*(a1+4.d0*a2))*(b2-b1)*(c2-c1)/(60.d0*hp)
    case default
       write(6,*)  'XX falut ID in fem_int_point, id= ',id
       stop
    end select

  end function fem_int_point

  subroutine inv_int(ne_in,val,x1,x2,x)

    integer(4), intent(in)  :: ne_in
    real(8), intent(in)  :: val, x1, x2
    real(8), intent(out) :: x

    integer(4) :: ne
    real(8) :: f(1:4), a1, a2, suml, lhs

    ! left grid

    ne = ne_in
    a1 = x1

    f(1) = hpsi(ne) / 3.d0 * a1
    f(3) = hpsi(ne) / 6.d0 * a1

    suml = f(1) + f(3)

    f(2) = hpsi(ne) / 6.d0
    f(4) = hpsi(ne) / 3.d0

    lhs = f(2) + f(4)

    ! right grid

    ne = ne_in + 1
    a2 = x2

    f(2) = hpsi(ne) / 6.d0 * a2
    f(4) = hpsi(ne) / 3.d0 * a2

    suml = 0.5d0 * (suml + f(2) + f(4))

    f(1) = hpsi(ne) / 3.d0
    f(3) = hpsi(ne) / 6.d0

    lhs = 0.5d0 * (lhs + f(1) + f(3))

    x = (- val - suml) / lhs

  end subroutine inv_int

end module tx_core_module

!*****************************************************************

!***************************************************************
!
!   For no '*** MATH LIBRARY ERROR 14: DEXP(X) UNDERFLOW'
!
!***************************************************************

pure REAL(8) FUNCTION EXPV(X)

  implicit none
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
!   SUBROUTINE APpend Integer(4) TO Strings
!     INPUT  : STR, NSTR, I
!              STR(NSTR(original)+1:NSTR(return))
!              NSTR : Number of STR. First, NSTR = 0.
!     OUTPUT : STR, NSTR
!
!***************************************************************

SUBROUTINE APITOS(STR, NSTR, I)

  implicit none
  character(len=*), INTENT(INOUT) :: STR
  INTEGER(4),       INTENT(INOUT) :: NSTR
  INTEGER(4),       INTENT(IN)    :: I

  INTEGER(4) :: J, NSTRI
  character(len=25) :: KVALUE

  WRITE(KVALUE,'(I25)') I
  J = index(KVALUE,' ',.true.)
  NSTRI = 25 - J
  STR(NSTR+1:NSTR+NSTRI) = KVALUE(J+1:25)
  NSTR = NSTR + NSTRI

  RETURN
END SUBROUTINE APITOS

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

  implicit none
  character(len=*), INTENT(INOUT) :: STR
  INTEGER(4),       INTENT(INOUT) :: NSTR
  character(len=*), INTENT(IN)    :: INSTR
  INTEGER(4),       INTENT(IN)    :: NINSTR

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

  implicit none
  character(len=*), INTENT(INOUT) :: STR
  INTEGER(4),       INTENT(INOUT) :: NSTR
  REAL(8),          INTENT(IN)    :: D
  character(len=*), INTENT(IN)    :: FORM

  INTEGER(4) :: IND
  INTEGER(4) :: L, IS, IE, NSTRD, IST
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

  IS = index(KVALUE,' ',.true.) + 1
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

  implicit none
  character(len=*), INTENT(INOUT) :: STR
  INTEGER(4),       INTENT(INOUT) :: NSTR
  REAL(4),          INTENT(IN)    :: GR
  character(len=*), INTENT(IN)    :: FORM

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
!   Separate string at char (similar to KSPLIT in task/lib/libchar.f90)
!
!***************************************************************

SUBROUTINE KSPLIT_TX(KKLINE,KID,KKLINE1,KKLINE2)

  IMPLICIT NONE
  CHARACTER(LEN=*),  INTENT(IN)  :: KKLINE
  CHARACTER(LEN=1),  INTENT(IN)  :: KID
  CHARACTER(LEN=*), INTENT(OUT) :: KKLINE1, KKLINE2
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
  
END SUBROUTINE KSPLIT_TX

!***************************************************************
!
!   Derivative
!
!***************************************************************

SUBROUTINE DERIVS1D(R,F,NRMAX,G)

  implicit none
  real(8), dimension(0:NRMAX), intent(in)  :: R, F
  real(8), dimension(0:NRMAX), intent(out) :: G
  integer(4), intent(in) :: NRMAX
  integer(4) :: NR
  real(8) :: DR1, DR2

  NR = 0
  DR1 = R(NR+1) - R(NR)
  DR2 = R(NR+2) - R(NR)
  G(NR) = (DR2**2 * F(NR+1) - DR1**2 * F(NR+2)) / (DR1 * DR2 * (DR2 - DR1)) &
       &- (DR2 + DR1) * F(NR) / (DR1 * DR2)
  DO NR = 1, NRMAX - 1
     DR1 = R(NR-1) - R(NR)
     DR2 = R(NR+1) - R(NR)
     G(NR) = (DR2**2 * F(NR-1) - DR1**2 * F(NR+1)) / (DR1 * DR2 * (DR2 - DR1)) &
          &- (DR2 + DR1) * F(NR) / (DR1 * DR2)
  END DO
  NR = NRMAX
  DR1 = R(NR-1) - R(NR)
  DR2 = R(NR-2) - R(NR)
  G(NR) = (DR2**2 * F(NR-1) - DR1**2 * F(NR-2)) / (DR1 * DR2 * (DR2 - DR1)) &
       &- (DR2 + DR1) * F(NR) / (DR1 * DR2)

END SUBROUTINE DERIVS1D

SUBROUTINE DERIVS2D(R,F,LQ,NQMAX,NRMAX,G)

  implicit none
  real(8), dimension(0:NRMAX), intent(in)  :: R
  real(8), dimension(1:NQMAX,0:NRMAX), intent(in)  :: F
  real(8), dimension(0:NRMAX), intent(out) :: G
  integer(4), intent(in) :: LQ, NRMAX, NQMAX
  integer(4) :: NR
  real(8) :: DR1, DR2

  NR = 0
  DR1 = R(NR+1) - R(NR)
  DR2 = R(NR+2) - R(NR)
  G(NR) = (DR2**2 * F(LQ,NR+1) - DR1**2 * F(LQ,NR+2)) / (DR1 * DR2 * (DR2 - DR1)) &
       &- (DR2 + DR1) * F(LQ,NR) / (DR1 * DR2)
  DO NR = 1, NRMAX - 1
     DR1 = R(NR-1) - R(NR)
     DR2 = R(NR+1) - R(NR)
     G(NR) = (DR2**2 * F(LQ,NR-1) - DR1**2 * F(LQ,NR+1)) / (DR1 * DR2 * (DR2 - DR1)) &
          &- (DR2 + DR1) * F(LQ,NR) / (DR1 * DR2)
  END DO
  NR = NRMAX
  DR1 = R(NR-1) - R(NR)
  DR2 = R(NR-2) - R(NR)
  G(NR) = (DR2**2 * F(LQ,NR-1) - DR1**2 * F(LQ,NR-2)) / (DR1 * DR2 * (DR2 - DR1)) &
       &- (DR2 + DR1) * F(LQ,NR) / (DR1 * DR2)

END SUBROUTINE DERIVS2D

pure REAL(8) FUNCTION DERIVF(NR,R,F,LQ,NQMAX,NRMAX)

  implicit none
  real(8), dimension(0:NRMAX), intent(in)  :: R
  real(8), dimension(1:NQMAX,0:NRMAX), intent(in)  :: F
  integer(4), intent(in) :: NR, LQ, NRMAX, NQMAX
  real(8) :: DR1, DR2

  IF(NR == 0) THEN
     DR1 = R(NR+1) - R(NR)
     DR2 = R(NR+2) - R(NR)
     DERIVF = (DR2**2 * F(LQ,NR+1) - DR1**2 * F(LQ,NR+2)) / (DR1 * DR2 * (DR2 - DR1)) &
          &- (DR2 + DR1) * F(LQ,NR) / (DR1 * DR2)
  ELSE IF(NR == NRMAX) THEN
     DR1 = R(NR-1) - R(NR)
     DR2 = R(NR-2) - R(NR)
     DERIVF = (DR2**2 * F(LQ,NR-1) - DR1**2 * F(LQ,NR-2)) / (DR1 * DR2 * (DR2 - DR1)) &
          &- (DR2 + DR1) * F(LQ,NR) / (DR1 * DR2)
  ELSE
     DR1 = R(NR-1) - R(NR)
     DR2 = R(NR+1) - R(NR)
     DERIVF = (DR2**2 * F(LQ,NR-1) - DR1**2 * F(LQ,NR+1)) / (DR1 * DR2 * (DR2 - DR1)) &
          &- (DR2 + DR1) * F(LQ,NR) / (DR1 * DR2)
  END IF

  RETURN
END FUNCTION DERIVF

! *** Universal high-accuracy 1st-derivative routine ***
!
!   Mode=0 gives you higher accuracy of the derivative.
!   Mode=1 gives you faster calculation of the derivative.
!
!   Input: x(0:nmax)     : radial coordinate
!          f(0:nmax)     : function
!          nmax          : size of array
!          mode          : 0 : Spline with 1st derivatives
!                              with O(4) accuracy at boundaries (recommended)
!                          1 : Taylor expansion derivative
!                              with O(4) accuracy
!   Output: dfdx(0:nmax) : derivative of the function with respect to x

function dfdx(x,f,nmax,mode)

  implicit none
  integer(4), intent(in) :: nmax, mode
  real(8), dimension(0:nmax), intent(in) :: x, f
  real(8), dimension(0:nmax) :: dfdx

  integer(4) :: n, ierr
  real(8), dimension(0:nmax) :: deriv
  real(8), dimension(1:4,0:nmax) :: u
  real(8) :: deriv4

  dfdx(0)    = deriv4(   0,x,f,nmax,0)
  dfdx(nmax) = deriv4(nmax,x,f,nmax,0)
  if(mode == 0) then
     deriv(0)    = dfdx(0)
     deriv(nmax) = dfdx(nmax)
     call spl1d(x,f,deriv,u,nmax+1,3,ierr)
     if(ierr /= 0) stop 'dfdx: SPL1D error!'
     dfdx(1:nmax-1) = deriv(1:nmax-1)
  else
     do n = 1, nmax-1
        dfdx(n) = deriv4(n,x,f,nmax,0)
     end do
  end if

end function dfdx

!***************************************************************
!
!   Integral Method
!
!***************************************************************

REAL(8) FUNCTION INTG_F(X)

  ! Calculate \int (r * X) dpsi

  ! === ATTENTION !!! ======================================================!
  !   This function can be used only if size(X) is equivalent to NEMAX,     !
  !   denoting the maximum index of the integral domain of fem_int funtion. !
  ! ========================================================================!

  use tx_core_module, only : fem_int
  implicit none
  real(8), dimension(*), intent(in) :: X

  INTG_F = 0.5D0 * SUM(fem_int(-1,X))

END FUNCTION INTG_F

REAL(8) FUNCTION INTG_P(X,NR,ID)

  ! Integrate X at a certain ONE point (NOT the domain)

  ! Calculate \int r * X(r) dpsi or 0.5 * \int X(psi) dpsi (ID == 0) 
  !           \int     X(r) dpsi (ID == 1)

  use tx_core_module, only : fem_int_point
  implicit none
  integer(4), intent(in) :: NR, ID
  real(8), dimension(*), intent(in) :: X
  integer(4) :: NE

  IF(ID == 0) THEN
     IF(NR == 0) THEN
        INTG_P = 0.D0
     ELSE
        NE = NR
        INTG_P = 0.5D0 * SUM(fem_int_point(-1,NE,X))
     END IF
  ELSEIF(ID == 1) THEN
     IF(NR == 0) THEN
        INTG_P = 0.D0
     ELSE
        NE = NR
        INTG_P = 0.5D0 * SUM(fem_int_point(21,NE,X))
     END IF
  END IF

END FUNCTION INTG_P

SUBROUTINE VALINT_SUB(X,NRLMAX,VAL,NR_START)

  ! Integrate X in the arbitrary size domain from NR_START to NRLMAX

  ! Calculate \int_{r(NR_START)}^r(NRLMAX) (r * X) dr

  use tx_core_module, only : fem_int_point
  implicit none
  integer(4), intent(in) :: NRLMAX
  integer(4), intent(in), optional :: NR_START
  real(8), dimension(*), intent(in) :: X
  real(8), intent(out) :: VAL
  integer(4) :: NE, NEMAX, NE_START
  real(8) :: SUML

  NEMAX = NRLMAX

  IF(PRESENT(NR_START)) THEN
     NE_START = NR_START
  ELSE
     ! default
     NE_START = 1
  END IF

  SUML = 0.D0
  DO NE = NE_START, NEMAX
     SUML = SUML + 0.5D0 * SUM(fem_int_point(-1,NE,X))
  END DO
  VAL = SUML

END SUBROUTINE VALINT_SUB

! *** Integral Method By Inverting Derivative Method ***

!     This formula can be used only if intX(0    ) is known of FVAL (ID = 0)
!                                   or intX(NRMAX) is known of FVAL (ID = else).

SUBROUTINE INTDERIV3(X,R,intX,FVAL,NRMAX,ID)

  implicit none
  integer(4), intent(in) :: NRMAX, ID
  real(8), intent(in), dimension(0:NRMAX) :: X, R
  real(8), intent(in) :: FVAL
  real(8), intent(out), dimension(0:NRMAX) :: intX
  integer(4) :: NR
  real(8) :: D1, D2, D3

  IF(ID == 0) THEN
     intX(0) = FVAL
     D1 = R(1) - R(0)
     D2 = R(2) - R(0)
     D3 = R(2) - R(1)
     intX(1) = ( D1*D2*(D2-D1)*X(0)+D1*D3*(D1+D3)*X(1)) &
          &  / (-D1**2+D2**2+D3**2) + FVAL
     intX(2) = ( D1**2*D2*(D2-D1)*X(0)+D2*D3**2*(D1-D2)*X(0) &
          &     +D2**2*D3*(D1+D3)*X(1)) / (D1*(-D1**2+D2**2+D3**2)) + FVAL
     DO NR = 3, NRMAX
        D1 = R(NR-2) - R(NR-1)
        D2 = R(NR  ) - R(NR-1)
        intX(NR) = (D2/D1)**2*intX(NR-2)-((D2/D1)**2-1.D0)*intX(NR-1) &
             &    -X(NR-1)*D2/D1*(D2-D1)
     END DO
  ELSE
     intX(NRMAX) = FVAL
     D1 = R(NRMAX-1) - R(NRMAX  )
     D2 = R(NRMAX-2) - R(NRMAX  )
     D3 = R(NRMAX-2) - R(NRMAX-1)
     intX(NRMAX-1) = ( D1*D2*(D2-D1)*X(NRMAX)+D1*D3*(D1+D3)*X(NRMAX-1)) &
          &        / (-D1**2+D2**2+D3**2) + FVAL
     intX(NRMAX-2) = ( D1**2*D2*(D2-D1)*X(NRMAX  )+D2*D3**2*(D1-D2)*X(NRMAX) &
          &           +D2**2*D3*(D1+D3)*X(NRMAX-1)) / (D1*(-D1**2+D2**2+D3**2)) + FVAL
     DO NR = NRMAX - 3, 0, -1
        D1 = R(NR+2) - R(NR+1)
        D2 = R(NR  ) - R(NR+1)
        intX(NR) = (D2/D1)**2*intX(NR+2)-((D2/D1)**2-1.D0)*intX(NR+1) &
             &    -X(NR+1)*D2/D1*(D2-D1)
     END DO
  END IF

END SUBROUTINE INTDERIV3

!***************************************************************
!
!   Coefficient function of CDBM model
!
!***************************************************************

pure REAL(8) FUNCTION TRCOFS(S,ALFA,RKCV)

  implicit none
  real(8), intent(in) :: S, ALFA, RKCV
  real(8) :: SA, FS1, FS2

  IF(ALFA > 0.D0) THEN
     SA = S - ALFA
     IF(SA > 0.D0) THEN
        FS1 = (1.D0 + 9.0D0 * SQRT(2.D0) * SA**2.5D0) &
             &  / (SQRT(2.D0) * (1.D0 - 2.D0 * SA + 3.D0 * SA**2 + 2.0D0 * SA**3))
     ELSE
        FS1 = 1.D0 / SQRT(2.D0 * (1.D0 - 2.D0 * SA) * (1.D0 - 2.D0 * SA + 3.D0 * SA**2))
     ENDIF
     IF(RKCV > 0.D0) THEN
        FS2 = SQRT(RKCV)**3 / S**2
     ELSE
        FS2 = 0.D0
     ENDIF
  ELSE
     SA = ALFA - S
     IF(SA > 0.D0) THEN
        FS1 = (1.D0 + 9.0D0 * SQRT(2.D0) * SA**2.5D0) &
             &  / (SQRT(2.D0) * (1.D0 - 2.D0 * SA + 3.D0 * SA**2 + 2.0D0 * SA**3))
     ELSE
        FS1 = 1.D0 / SQRT(2.D0 * (1.D0 - 2.D0 * SA) * (1.D0 - 2.D0 * SA + 3.D0 * SA**2))
     ENDIF
     IF(RKCV < 0.D0) THEN
        FS2 = SQRT(-RKCV)**3 / S**2
     ELSE
        FS2 = 0.D0
     ENDIF
  ENDIF
  TRCOFS = MAX(FS1,FS2)

END FUNCTION TRCOFS

!***************************************************************
!
!   Coefficient function of CDIM model
!   09/06/17~ miki_m
!   given up, txcalv
!
!***************************************************************

!pure REAL(8) FUNCTION TRCOFSIM(S,ALFA)

!  implicit none
!  real(8), intent(in) :: S
!  real(8) :: OMEGAPR,  RAQPR

!  OMEGAPR = (R(NR) / RR)**2.D0 * (NCph / NCth) / (R / RA)**2.D0 * RAQPR

!  RAQPR = 2.D0 * R(0:NRMAX) * dfdx (PSI , (R / RA)**4 / Q , NMAX , 0 )
! cf. txcalv.f90 L501

!  TRCOFSIM = 3 * (OMEGAPR / 2.D0)**1.5D0 * (RR / RA)**1.5D0 / (Q(NR) * S**2.D0)

!END FUNCTION TRCOFSIM


!***************************************************************
!
!   Mesh generating function
!
!***************************************************************

! This function is the one obtained by integration of the Lorentz function
!   R  : radial coordinate
!   C  : amplitude factor
!   W  : width of flat region around R=RC
!   RC : center radial point of fine mesh region

pure REAL(8) FUNCTION LORENTZ(R,C1,C2,W1,W2,RC1,RC2,AMP)

  implicit none
  real(8), intent(in) :: r, c1, c2, w1, w2, rc1, rc2
  real(8), intent(in), optional :: AMP

  LORENTZ = R + C1 * ( W1 * ATAN((R - RC1) / W1) + W1 * ATAN(RC1 / W1)) &
       &      + C2 * ( W2 * ATAN((R - RC2) / W2) + W2 * ATAN(RC2 / W2))
  if(present(amp)) LORENTZ = LORENTZ / AMP

END FUNCTION LORENTZ

pure REAL(8) FUNCTION LORENTZ_PART(R,W1,W2,RC1,RC2,ID)

  implicit none
  real(8), intent(in) :: r, w1, w2, rc1, rc2
  integer(4), intent(in) :: ID

  IF(ID == 0) THEN
     LORENTZ_PART = W1 * ATAN((R - RC1) / W1) + W1 * ATAN(RC1 / W1)
  ELSE
     LORENTZ_PART = W2 * ATAN((R - RC2) / W2) + W2 * ATAN(RC2 / W2)
  END IF

END FUNCTION LORENTZ_PART

!***************************************************************
!
!   Bisection method for solving the equation
!
!***************************************************************

! Bisection method can solve the equation only if the solution is unique
! in the designated region.
!   f      (in)  : the function which is set to LHS in the equation
!   cl     (in)  : argument of f
!   w      (in)  : argument of f
!   rc     (in)  : argument of f
!   amp    (in)  : argument of f
!   s      (in)  : the value which is set to RHS in the equation
!   valmax (in)  : maximum value of codomain
!   val    (out) : solution
!   valmin (in)  : minimum value of codomain, optional
! i.e. we now handle the equation of "f = s" or "f - s = 0".
!   eps    : arithmetic precision

SUBROUTINE BISECTION(f,cl1,cl2,w1,w2,rc1,rc2,amp,s,valmax,val,valmin)

  implicit none
  real(8), external :: f
  real(8), intent(in) :: cl1, cl2, w1, w2, rc1, rc2, amp, s, valmax
  real(8), intent(in), optional :: valmin
  real(8), intent(out) :: val
  integer(4) :: i, n
  real(8) :: a, b, c, eps, fa, fc

  if(present(valmin)) then
     a = valmin
  else
     a = 0.d0
  end if
  b = valmax
  eps = 1.d-10
  n = log10((b - a) / eps) / log10(2.d0) + 0.5d0
  fa = f(a,cl1,cl2,w1,w2,rc1,rc2,amp) - s
  do i = 1, n
     c = 0.5d0 * (a + b)
     fc = f(c,cl1,cl2,w1,w2,rc1,rc2,amp) - s
     if(fa * fc < 0.d0) then
        b = c
     else
        a = c
     end if
  end do
  val = c

END SUBROUTINE BISECTION

!************************************************************************************
!
!   Interpolate and extrapolate data loaded from file
!     INPUT  : nmax_in  : the number of the data point of the file
!              r_in     : radial coordinate of the file
!              dat_in   : the data of the file
!              nmax_std : the number of the data point
!              r_std    : radial coordinate
!              iedge    : indication of the treatment at the axis and the outer boundary
!              ideriv   : (optional) r-derivative at the axis becomes zero when imode=4,5,6
!     OUTPUT : dat_out  : interpolated and extrapolated value of dat_in
!              nrbound  : (optional) outermost index corresponding to rho_in(nmax_in)
!
!      imode |  iaxis  iedge  |    axis      outer boundary |
!     -------------------------------------------------------
!        1   |    1      1    |     0              0        |
!        2   |    1      2    |     0              ex       |
!        3   |    1      3    |     0              val      |
!        4   |    2      1    |     ex             0        |
!        5   |    2      2    |     ex             ex       |
!        6   |    2      3    |     ex             val      |
!        7   |    3      1    |     val            0        |
!        8   |    3      2    |     val            ex       |
!        9   |    3      3    |     val            val      |
!
!     0 : zero, ex : extrapolation, val : three quarters of the nearest value
!************************************************************************************

subroutine inexpolate(nmax_in,rho_in,dat_in,nmax_std,rho_std,imode,dat_out,ideriv,nrbound,idx)
  implicit none
  integer(4), intent(in) :: nmax_in, nmax_std, imode
  integer(4), intent(in),  optional :: ideriv, idx
  integer(4), intent(out), optional :: nrbound
  real(8), dimension(1:nmax_in),  intent(in)  :: rho_in, dat_in
  real(8), dimension(0:nmax_std), intent(in)  :: rho_std
  real(8), dimension(0:nmax_std), intent(out) :: dat_out
  integer(4) :: i, iaxis, iedge, nmax, isep
  real(8) :: rhoa
  real(8), dimension(:), allocatable :: rho_tmp, dat_tmp
  real(8) :: aitken2p, fctr

  isep = 0

  if(imode < 1 .or. imode > 9) stop 'inexpolate: inappropriate imode'

  ! Check whether the outermost value of the input is over the separatrix or not
  if(rho_in(nmax_in) >= 1.d0) isep = 1
  iedge = mod((imode-1),3)+1

  ! Reshape the mesh and data arrays
  nmax = nmax_in
  if(rho_in(1) == 0.d0) then ! Originally having the value at the axis
     nmax = nmax - 1
     allocate(rho_tmp(0:nmax),dat_tmp(0:nmax))
     rho_tmp(0:nmax_in-1) = rho_in(1:nmax_in)
     dat_tmp(0:nmax_in-1) = dat_in(1:nmax_in)
     iaxis = 0
  else ! Not so
     allocate(rho_tmp(0:nmax),dat_tmp(0:nmax))
     rho_tmp(0) = 0.d0
     dat_tmp(0) = 0.d0
     rho_tmp(1:nmax_in) = rho_in(1:nmax_in)
     dat_tmp(1:nmax_in) = dat_in(1:nmax_in)
     iaxis = 1
  end if

  ! Extrapolate the value at the axis
  if(iaxis == 0) then
     dat_out(0) = dat_tmp(0)
  else
     iaxis = int((imode-1)/3)+1
     if(iaxis == 1) then
        dat_tmp(0) = 0.d0
     else if(iaxis == 2) then
        dat_tmp(0) = aitken2p(0.d0,dat_tmp(1),dat_tmp(2),dat_tmp(3), &
             &                     rho_tmp(1),rho_tmp(2),rho_tmp(3))
        if(present(ideriv)) then
           if(ideriv == 1) then
              dat_tmp(0) = fctr(rho_tmp(1),rho_tmp(2),dat_tmp(1),dat_tmp(2))
           else
              call sctr(           rho_tmp(1),rho_tmp(2),rho_tmp(3),rho_tmp(4), &
                   &                          dat_tmp(2),dat_tmp(3),dat_tmp(4), &
                   &    dat_tmp(0),dat_tmp(1))
              dat_out(1) = dat_tmp(1)
           end if
        end if
     else if(iaxis == 3) then
        dat_tmp(0) = 0.75d0 * dat_tmp(1)
     end if
     dat_out(0) = dat_tmp(0)
  end if

  ! Extrapolate the value at the outer boundary
  if(isep == 0) then
     rhoa = 1.d0
  else
     rhoa = rho_tmp(nmax)
  end if
  if(iedge == 1) then
     rho_tmp(nmax) = rhoa
     dat_tmp(nmax) = 0.d0
  else if(iedge == 2) then
     rho_tmp(nmax) = rhoa
     dat_tmp(nmax) = aitken2p(rhoa,dat_tmp(nmax-1),dat_tmp(nmax-2),dat_tmp(nmax-3), &
             &                     rho_tmp(nmax-1),rho_tmp(nmax-2),rho_tmp(nmax-3))
  else if(iedge == 3) then
     rho_tmp(nmax) = rhoa
     dat_tmp(nmax) = 0.75d0 * dat_tmp(nmax-1)
  end if

  ! Interpolate the value
  do i = 1, nmax_std
     if(rho_std(i) < rho_tmp(nmax)) then
        call aitken(rho_std(i),dat_out(i),rho_tmp,dat_tmp,2,size(dat_tmp))
     else
        if(present(nrbound)) nrbound = i - 1
        exit
     end if
  end do

  deallocate(rho_tmp,dat_tmp)

  ! Avoid negative values
  if(present(idx) .eqv. .false.) where(dat_out < 0.d0) dat_out = 0.d0

!!$  do i=0,nmax_std
!!$     write(6,*) rho_std(i),dat_out(i)
!!$  end do

end subroutine inexpolate

!***************************************************************
!
!   Gaussian (Maxwellian) distribution
!
!   (input)
!     x     : position
!     mu    : average
!     sigma : standard deviation
!
!***************************************************************

pure real(8) function fgaussian(x,mu,sigma,norm) result(f)
  real(8), parameter :: pi = 3.14159265358979323846d0
  real(8), intent(in) :: x, mu, sigma
  integer(4), intent(in), optional :: norm

  f = exp(- (x - mu)**2 / (2.d0 * sigma**2))
  if(present(norm)) f = f / (sqrt(2.d0 * pi) * sigma)

end function fgaussian
