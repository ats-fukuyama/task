!     $Id$
module core_module
  use commons, only : nrmax, h, r, nra, psi, hpsi, nemax
  implicit none
  public

contains

  function fem_int(id,ne,a,b) result(x)
!-------------------------------------------------------
!
!   Calculate "\int_{psi_i}^{psi_{i+1}} function(psi) dpsi"
!      dpsi : mesh interval
!      a    : coefficient vector
!      u    : variable vector
!      w    : weighting vector
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
!      id = 44 : a * b * w
!
!   where ' means the derivative of psi
!
!  < input >
!     id       : mode select
!     ne       : current number of elements
!     nnode    : maximum of nodes
!     a(nnode) : coefficient vector, optional
!     b(nnode) : coefficient vector, optional
!  < output >
!     x(4)     : matrix of integrated values
!
!-------------------------------------------------------
    integer, intent(in) :: id, ne
    real(8), intent(in), dimension(0:nrmax), optional  :: a, b
    integer :: node1, node2
    real(8) :: x(1:4), a1, a2, r1, r2, p1, p2, b1, b2, hp
    
    node1 = ne-1  ; node2 = ne
    if(present(a)) then
       a1 = a(node1) ; a2 = a(node2)
       if(present(b)) then
          b1 = b(node1) ; b2 = b(node2)
       end if
    end if
    r1 = r(node1) ; r2 = r(node2)
    p1 = psi(node1) ; p2 = psi(node2) ; hp = hpsi(ne)

    select case(id)
    case(-1)
       x(1) =-0.5d0 * a1
       x(2) =-0.5d0 * a2
       x(3) = 0.5d0 * a1
       x(4) = 0.5d0 * a2
    case(0)
       x(1) = hp / 3.d0 * a1
       x(2) = hp / 6.d0 * a2
       x(3) = hp / 6.d0 * a1
       x(4) = hp / 3.d0 * a2
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
!!$       x(1) = (3.d0 * r1 +        r2) * hp / 12.d0 * a1
!!$       x(2) = (       r1 +        r2) * hp / 12.d0 * a2
!!$       x(3) = (       r1 +        r2) * hp / 12.d0 * a1
!!$       x(4) = (       r1 + 3.d0 * r2) * hp / 12.d0 * a2
       x(1) = 2.d0*h(ne)*( 15.d0*p1**2+45.d0*p1*r1*r2+48.d0*p1*p2+24.d0*r1*p2*r2 &
            &             + 8.d0*p2**2) / (105.d0*(r1+r2)**2) * a1
       x(2) = 4.d0*h(ne)*(  3.d0*p1**2+ 9.d0*p1*r1*r2+11.d0*p1*p2+ 9.d0*r1*p2*r2 &
            &             + 3.d0*p2**2) / (105.d0*(r1+r2)**2) * a2
       x(3) = 4.d0*h(ne)*(  3.d0*p1**2+ 9.d0*p1*r1*r2+11.d0*p1*p2+ 9.d0*r1*p2*r2 &
            &             + 3.d0*p2**2) / (105.d0*(r1+r2)**2) * a1
       x(4) = 2.d0*h(ne)*(  8.d0*p1**2+24.d0*p1*r1*r2+48.d0*p1*p2+45.d0*r1*p2*r2 &
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
       x(1) = ( 12.d0*a1*b1 + 3.d0*a2*b1 + 3.d0*a1*b2 + 2.d0*a2*b2) * hp / 60.d0 &
            &+( 12.d0*p1*b1 + 3.d0*p2*b1 + 3.d0*p1*b2 + 2.d0*p2*b2) * (-a1+a2) / 60.d0 &
            &+( 12.d0*p1*a1 + 3.d0*p2*a1 + 3.d0*p1*a2 + 2.d0*p2*a2) * (-b1+b2) / 60.d0 &
            &-( 12.d0*p1*a1*b1 + 3.d0*p2*a1*b1 + 3.d0*p1*a2*b1 + 2.d0*p2*a2*b1 &
            &  + 3.d0*p1*a1*b2 + 2.d0*p2*a1*b2 + 2.d0*p1*a2*b2 + 3.d0*p2*a2*b2) / 60.d0
       x(2) = (  3.d0*a1*b1 + 2.d0*a2*b1 + 2.d0*a1*b2 + 3.d0*a2*b2) * hp / 60.d0 &
            &+(  3.d0*p1*b1 + 2.d0*p2*b1 + 2.d0*p1*b2 + 3.d0*p2*b2) * (-a1+a2) / 60.d0 &
            &+(  3.d0*p1*a1 + 2.d0*p2*a1 + 2.d0*p1*a2 + 3.d0*p2*a2) * (-b1+b2) / 60.d0 &
            &+( 12.d0*p1*a1*b1 + 3.d0*p2*a1*b1 + 3.d0*p1*a2*b1 + 2.d0*p2*a2*b1 &
            &  + 3.d0*p1*a1*b2 + 2.d0*p2*a1*b2 + 2.d0*p1*a2*b2 + 3.d0*p2*a2*b2) / 60.d0
       x(3) = (  3.d0*a1*b1 + 2.d0*a2*b1 + 2.d0*a1*b2 + 3.d0*a2*b2) * hp / 60.d0 &
            &+(  3.d0*p1*b1 + 2.d0*p2*b1 + 2.d0*p1*b2 + 3.d0*p2*b2) * (-a1+a2) / 60.d0 &
            &+(  3.d0*p1*a1 + 2.d0*p2*a1 + 2.d0*p1*a2 + 3.d0*p2*a2) * (-b1+b2) / 60.d0 &
            &-(  3.d0*p1*a1*b1 + 2.d0*p2*a1*b1 + 2.d0*p1*a2*b1 + 3.d0*p2*a2*b1 &
            &  + 2.d0*p1*a1*b2 + 3.d0*p2*a1*b2 + 3.d0*p1*a2*b2 +12.d0*p2*a2*b2) / 60.d0
       x(4) = (  2.d0*a1*b1 + 3.d0*a2*b1 + 3.d0*a1*b2 +12.d0*a2*b2) * hp / 60.d0 &
            &+(  2.d0*p1*b1 + 3.d0*p2*b1 + 3.d0*p1*b2 +12.d0*p2*b2) * (-a1+a2) / 60.d0 &
            &+(  2.d0*p1*a1 + 3.d0*p2*a1 + 3.d0*p1*a2 +12.d0*p2*a2) * (-b1+b2) / 60.d0 &
            &+(  3.d0*p1*a1*b1 + 2.d0*p2*a1*b1 + 2.d0*p1*a2*b1 + 3.d0*p2*a2*b1 &
            &  + 2.d0*p1*a1*b2 + 3.d0*p2*a1*b2 + 3.d0*p1*a2*b2 +12.d0*p2*a2*b2) / 60.d0
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
       x(1) = (-3.d0*a1*b1 -      a1*b2 -      a2*b1 -      a2*b2) / 12.d0 &
          &  +( 3.d0*p1*b1+p2*b1+p1*b2+     p2*b2) * (a1-a2) / (12.d0 * hp) &
          &  +(b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
          & / (12.d0 * hp) &
          &  +( 3.d0*p1*a1 + p2*a1 + p1*a2 +      p2*a2) * (b1 - b2) / (12.d0 * hp)
       x(2) = ( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) / 12.d0 &
          &  -( 3.d0*p1*b1+p2*b1+p1*b2+     p2*b2) * (a1-a2) / (12.d0 * hp) &
          &  -(b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
          & / (12.d0 * hp) &
          &  +(     p1*a1 + p2*a1 + p1*a2 + 3.d0*p2*a2) * (b1 - b2) / (12.d0 * hp)
       x(3) = (-     a1*b1 -      a1*b2 -      a2*b1 - 3.d0*a2*b2) / 12.d0 &
          &  +(     p1*b1+p2*b1+p1*b2+3.d0*p2*b2) * (a1-a2) / (12.d0 * hp) &
          &  -(b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
          & / (12.d0 * hp) &
          &  -( 3.d0*p1*a1 + p2*a1 + p1*a2 +      p2*a2) * (b1 - b2) / (12.d0 * hp)
       x(4) = (      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) / 12.d0 &
          &  -(     p1*b1+p2*b1+p1*b2+3.d0*p2*b2) * (a1-a2) / (12.d0 * hp) &
          &  +(b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
          & / (12.d0 * hp) &
          &  -(     p1*a1 + p2*a1 + p1*a2 + 3.d0*p2*a2) * (b1 - b2) / (12.d0 * hp)
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
       x(1) = ( 3.d0 * a1 +        a2) * hp / 12.d0 * b1
       x(2) = (        a1 +        a2) * hp / 12.d0 * b2
       x(3) = (        a1 +        a2) * hp / 12.d0 * b1
       x(4) = (        a1 + 3.d0 * a2) * hp / 12.d0 * b2
    case default
       stop 'XX falut ID in fem_int'
    end select

  end function fem_int

  function lump_int(id,ne,a) result(x)

    integer, intent(in) :: id, ne
    real(8), intent(in), dimension(0:nrmax), optional  :: a
    integer :: node1, node2
    real(8) :: x(1:4), a1, a2, r1, r2, p1, p2, hp
    
    node1 = ne-1  ; node2 = ne
    if(present(a)) then
       a1 = a(node1) ; a2 = a(node2)
    end if
    r1 = r(node1) ; r2 = r(node2)
    p1 = psi(node1) ; p2 = psi(node2) ; hp = hpsi(ne)

    select case(id)
    case(1)
       x(1) = 0.5d0 * hp
       x(2) = 0.d0
       x(3) = 0.d0
       x(4) = 0.5d0 * hp
    case(2)
       x(1) = ( 2.d0 * a1 +        a2) * hp / 6.d0
       x(2) = 0.d0
       x(3) = 0.d0
       x(4) = (        a1 + 2.d0 * a2) * hp / 6.d0
    case(8)
       x(1) =-1.d0 * hp
       x(2) = 0.d0
       x(3) = 0.d0
       x(4) = 1.d0 * hp
    end select

  end function lump_int

end module core_module

!*****************************************************************

module libraries
  implicit none
  private
  public :: EXPV, APTOS, TOUPPER, TRCOFS, DERIVS, DERIVF, INTDERIV3, &
       &    LORENTZ, LORENTZ_PART, BISECTION, VALINT_SUB, INTG_F, INTG_P

  interface APTOS
     module procedure APITOS
     module procedure APSTOS
     module procedure APRTOS
     module procedure APDTOS
  end interface

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

    character(len=*), INTENT(INOUT) :: STR
    INTEGER,          INTENT(INOUT) :: NSTR
    INTEGER,          INTENT(IN)    :: I

    INTEGER :: J, NSTRI
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

    character(len=*), INTENT(INOUT) :: STR
    INTEGER,          INTENT(INOUT) :: NSTR
    character(len=*), INTENT(IN)    :: INSTR
    INTEGER,          INTENT(IN)    :: NINSTR

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

    character(len=*), INTENT(INOUT) :: STR
    INTEGER,          INTENT(INOUT) :: NSTR
    REAL(8),          INTENT(IN)    :: D
    character(len=*), INTENT(IN)    :: FORM

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

    character(len=*), INTENT(INOUT) :: STR
    INTEGER,          INTENT(INOUT) :: NSTR
    REAL,             INTENT(IN)    :: GR
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

  SUBROUTINE DERIVS(R,F,G,NRMAX)

    real(8), dimension(0:NRMAX), intent(in)  :: R, F
    real(8), dimension(0:NRMAX), intent(out) :: G
    integer, intent(in) :: NRMAX
    integer :: NR
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

  END SUBROUTINE DERIVS

  REAL(8) FUNCTION DERIVF(NR,R,F,NRMAX,ID)
    integer, intent(in) :: NR, NRMAX, ID
    real(8), dimension(0:NRMAX), intent(in) :: F, R
    integer :: NR0, NR1, NR2
    real(8) :: DR, DLT1, DLT2, AITKEN4P, RETP, FETP

    IF(ID == 0) THEN
       IF(NR == 0) THEN
          NR0 = NR
          NR1 = NR+1
          NR2 = NR+2
       ELSEIF(NR == NRMAX) THEN
          NR0 = NR
          NR1 = NR-1
          NR2 = NR-2
       ELSE
          NR0 = NR
          NR1 = NR-1
          NR2 = NR+1
       ENDIF

       DLT1=R(NR1)-R(NR0)
       DLT2=R(NR2)-R(NR0)
       DERIVF = (DLT2**2*F(NR1)-DLT1**2*F(NR2))/(DLT1*DLT2*(DLT2-DLT1)) &
       &       -(DLT2+DLT1)*F(NR0)/(DLT1*DLT2)
    ELSE
       IF(NR == 0) THEN
          IF(ID == 1 .OR. ID == 3) THEN
             RETP = 2.D0 * R(NR) - R(NR+1)
             FETP = AITKEN4P(RETP,F(NR),F(NR+1),F(NR+2),F(NR+3),F(NR+4), &
                  &               R(NR),R(NR+1),R(NR+2),R(NR+3),R(NR+4))

             NR0 = NR
             NR1 = NR+1
             DLT1=R(NR1)-R(NR0)
             DLT2=RETP  -R(NR0)
             DERIVF = (DLT2**2*F(NR1)-DLT1**2*FETP)/(DLT1*DLT2*(DLT2-DLT1)) &
                  &  -(DLT2+DLT1)*F(NR0)/(DLT1*DLT2)
          ELSE
             NR0 = NR
             NR1 = NR+1
             NR2 = NR+2
             DLT1=R(NR1)-R(NR0)
             DLT2=R(NR2)-R(NR0)
             DERIVF = (DLT2**2*F(NR1)-DLT1**2*F(NR2))/(DLT1*DLT2*(DLT2-DLT1)) &
                  &  -(DLT2+DLT1)*F(NR0)/(DLT1*DLT2)
          END IF
       ELSEIF(NR == NRMAX) THEN
          IF(ID == 2 .OR. ID == 3) THEN
             RETP = 2.D0 * R(NR) - R(NR-1)
             FETP = AITKEN4P(RETP,F(NR),F(NR-1),F(NR-2),F(NR-3),F(NR-4), &
                  &               R(NR),R(NR-1),R(NR-2),R(NR-3),R(NR-4))

             NR0 = NR
             NR1 = NR-1
             DLT1=R(NR1)-R(NR0)
             DLT2=RETP  -R(NR0)
             DERIVF = (DLT2**2*F(NR1)-DLT1**2*FETP)/(DLT1*DLT2*(DLT2-DLT1)) &
                  &  -(DLT2+DLT1)*F(NR0)/(DLT1*DLT2)
          ELSE
             NR0 = NR
             NR1 = NR-1
             NR2 = NR-2
             DLT1=R(NR1)-R(NR0)
             DLT2=R(NR2)-R(NR0)
             DERIVF = (DLT2**2*F(NR1)-DLT1**2*F(NR2))/(DLT1*DLT2*(DLT2-DLT1)) &
                  &  -(DLT2+DLT1)*F(NR0)/(DLT1*DLT2)
          END IF
       ELSE
          NR0 = NR
          NR1 = NR-1
          NR2 = NR+1
          DLT1=R(NR1)-R(NR0)
          DLT2=R(NR2)-R(NR0)
          DERIVF = (DLT2**2*F(NR1)-DLT1**2*F(NR2))/(DLT1*DLT2*(DLT2-DLT1)) &
               &  -(DLT2+DLT1)*F(NR0)/(DLT1*DLT2)
       END IF
    END IF

    RETURN
  END FUNCTION DERIVF

!***************************************************************
!
!   Integral Method
!
!***************************************************************

  REAL(8) FUNCTION INTG_F(X)

    ! Calculate \int (r * X) dpsi

    use commons, only : NRMAX, NEMAX
    use core_module, only : fem_int
    real(8), dimension(0:NRMAX), intent(in) :: X
    integer :: NE
    real(8) :: SUML

    SUML = 0.D0
    DO NE = 1, NEMAX
       SUML = SUML + 0.5D0 * SUM(fem_int(0,NE,X))
    END DO
    INTG_F = SUML
    
  END FUNCTION INTG_F

  REAL(8) FUNCTION INTG_P(X,NR,ID)

    ! Calculate \int r * X(r) dpsi or 0.5 * \int X(psi) dpsi (ID == 0) 
    !           \int     X(r) dpsi (ID == 1) at one mesh

    use commons, only : NRMAX
    use core_module, only : fem_int
    integer, intent(in) :: NR, ID
    real(8), dimension(0:NRMAX), intent(in) :: X
    integer :: NE

    IF(ID == 0) THEN
       IF(NR == 0) THEN
          INTG_P = 0.D0
       ELSE
          NE = NR
          INTG_P = 0.5D0 * SUM(fem_int(0,NE,X))
       END IF
    ELSEIF(ID == 1) THEN
       IF(NR == 0) THEN
          INTG_P = 0.D0
       ELSE
          NE = NR
          INTG_P = 0.5D0 * SUM(fem_int(21,NE,X))
       END IF
    END IF
    
  END FUNCTION INTG_P

  SUBROUTINE VALINT_SUB(X,NRLMAX,VAL)

    ! Calculate \int_0^r (r * X) dr

    use commons, only : NRMAX
    use core_module, only : fem_int
    integer, intent(in) :: NRLMAX
    real(8), dimension(0:NRMAX), intent(in) :: X
    real(8), intent(out) :: VAL
    integer :: NE, NEMAX
    real(8) :: SUML

    NEMAX = NRLMAX
    SUML = 0.D0
    DO NE = 1, NEMAX
       SUML = SUML + 0.5D0 * SUM(fem_int(0,NE,X))
    END DO
    VAL = SUML
    
  END SUBROUTINE VALINT_SUB

! *** Integral Method By Inversing Derivative Method ***

!     This formula can be used only if intX(0    ) is known of FVAL (ID = 0)
!                                   or intX(NRMAX) is known of FVAL (ID = else).

  SUBROUTINE INTDERIV3(X,R,intX,FVAL,NRMAX,ID)
    integer, intent(in) :: NRMAX, ID
    real(8), intent(in), dimension(0:NRMAX) :: X, R
    real(8), intent(in) :: FVAL
    real(8), intent(out), dimension(0:NRMAX) :: intX
    integer :: NR
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

  REAL(8) FUNCTION TRCOFS(S,ALFA,RKCV)

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
!   Mesh generating function
!
!***************************************************************

! This function is the one obtained by integration of the Lorentz function
!   R  : radial coordinate
!   C  : amplitude factor
!   W  : width of flat region around R=RC
!   RC : center radial point of fine mesh region

  REAL(8) FUNCTION LORENTZ(R,C,W,RC,AMP)
    real(8), intent(in) :: r, c, w, rc
    real(8), intent(in), optional :: AMP

    LORENTZ = R + C * (W * ATAN((R - RC) / W) + W * ATAN(RC / W))
    if(present(amp)) LORENTZ = LORENTZ / AMP

  END FUNCTION LORENTZ

  REAL(8) FUNCTION LORENTZ_PART(R,W,RC)
    real(8), intent(in) :: r, w, rc

    LORENTZ_PART = W * ATAN((R - RC) / W) + W * ATAN(RC / W)

  END FUNCTION LORENTZ_PART

!***************************************************************
!
!   Bisection method for solving the equation
!
!***************************************************************

! Bisection method can solve the equation only if the solution is unique
! in the designated region.
!   f      : the function which is set to LHS in the equation
!   cl     : argument of f
!   w      : argument of f
!   rc     : argument of f
!   amp    : argument of f
!   s      : the value which is set to RHS in the equation
!   valmax : maximum value of codomain
!   val    : solution
! i.e. we now handle the equation of "f = s" or "f - s = 0".
!   eps    : arithmetic precision

  SUBROUTINE BISECTION(f,cl,w,rc,amp,s,valmax,val)
    real(8), external :: f
    real(8), intent(in) :: cl, w, rc, amp, s, valmax
    real(8), intent(out) :: val
    integer :: i, n
    real(8) :: a, b, c, eps, fa, fc

    a = 0.d0 ; b = valmax
    eps = 1.d-10
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

