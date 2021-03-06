!     $Id: txlib.f90,v 1.50 2011/04/13 05:50:25 honda Exp $
!!! Miscellaneous libraries NOT related to the physics or physical model
module tx_core_module
  use tx_commons, only : nrmax, vv, hv, nemax, SUPGstab
  implicit none
  real(8), parameter :: c13 = 1.d0 / 3.d0, c16 = 1.d0 / 6.d0, c112 = 1.d0 / 12.d0, c160 = 1.d0 / 60.d0
  real(8), parameter :: csq15 = 1.d0 / sqrt(15.d0) ! Raymond and Garder, Monthly Weather Review (1976)
  real(8) :: cp
  public

contains

  function fem_int(id,a,b,c) result(x)

!-------------------------------------------------------
!
!   Calculate "\int_0^{vv_b} function(vv) dvv"
!      dvv : mesh interval
!      a   : coefficient vector
!      b   : coefficient vector
!      u   : variable vector
!      w   : weighting vector
!
!   function(vv) is classified as follows:
!      (+-100)  : upwind SUPG
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
!      id = 20 : a * b * u * w
!      id = 21 : a * b'* u * w
!      id = 22 : a * b * u'* w
!      id = 23 : a * b * u * w'
!      id = 24 :(a * b * u)'* w
!      id = 25 : a * b'* u'* w
!      id = 26 : a * b'* u * w'
!      id = 27 : a * b * u'* w'
!
!      id = 30 : a * b * c * u * w
!      id = 31 : a * b * c'* u * w
!      id = 32 : a * b * c * u'* w
!      id = 33 : a * b * c * u * w'
!      id = 36 : a * b * c'* u * w'
!      id = 37 : a * b * c * u'* w'
!      id = 38 : a * b *(c * u)'*w' =  a * b * c * u'* w' + a * b * c'* u * w'
!
!      id = -1 : a * w
!      id = -2 : a * b * w
!      id = -3 :(a * b)'* w
!      id = -4 : a'* w
!      id = -8 : a * w'
!      id = -9 : a * b * w'
!
!      id =-20 : a * b * c * w
!      id =-22 : a * b * c'* w
!      id =-23 : a * b * c * w'
!
!      < obsolete >
!      id = 15 : vv * a * u * w
!      id = 16 : vv * a'* u * w
!      id = 17 : vv * a * u'* w
!      id = 18 : vv * a * u'* w'
!      id = 19 : vv * a * u * w'
!      id = 20 : vv * a * b'* u * w'
!
!      id = 36 : vv * a'* b * u * w
!      id = 37 : vv * a * b'* u * w
!      id = 38 : vv * a * b * u'* w
!      id = 39 :(vv * a * b * u)'* w
!      id = 40 : vv * a'* b * u * w'
!      id = 41 : vv * a * b * u'* w'
!      id = 42 :(vv * a * b * u)'* w'
!      id = 43 : vv * a * b * (u / b)'* w'
!
!      id = 44 : vv * a * b * u * w
!      id = 45 :(vv * a * b'* u)'* w
!      id = 46 :(vv * a * b'* u)'* w'
!      id = 47 : vv * a'* b'* u * w
!      id = 48 : vv * a * b'* u'* w
!
!      id = 49 : vv * a * b'* c'* u * w
!
!   where ' means the derivative of vv
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
    integer(4) :: ne, icheck
    real(8) :: x(1:nemax,1:4), a1, a2, b1, b2, c1, c2, hvlinv, &
         &     coef1, coef2, coef3, coef4, coef5, coef6, coef7
!    real(8) :: p1, p2

#ifdef debug
    icheck = 0
    select case(id)
    case(1,4,8,11,101)
       icheck = 0
    case(-101,-8,-4:-1,2:3,5:7,9:10,12:14,102,103,105,106)
       if(.not. present(a)) icheck = 1
    case(-102,-12,-9,20:27,120:122)
       if(.not. present(a) .or. .not. present(b)) icheck = 2
    case(-120,-23:-22,-20,30:33,36:38,130:132)
       if(.not. present(a) .or. .not. present(b) .or. .not. present(c)) icheck = 3
    case default
       icheck = 4
    end select
    if(icheck /= 0) then
       write(6,*) 'txlib: argument of fem_int mismatch! id=',id
       stop
    end if
#endif debug

    cp = SUPGstab * csq15

    select case(id)
    case(-1)
       do ne = 1, nemax
          x(ne,1) = hv(ne) * c13 * a(ne-1)
          x(ne,2) = hv(ne) * c16 * a(ne)
          x(ne,3) = hv(ne) * c16 * a(ne-1)
          x(ne,4) = hv(ne) * c13 * a(ne)
       end do
    case(-2)
       do ne = 1, nemax
          x(ne,1) = ( 3.d0 * a(ne-1) +        a(ne)) * hv(ne) * c112 * b(ne-1)
          x(ne,2) = (        a(ne-1) +        a(ne)) * hv(ne) * c112 * b(ne)
          x(ne,3) = (        a(ne-1) +        a(ne)) * hv(ne) * c112 * b(ne-1)
          x(ne,4) = (        a(ne-1) + 3.d0 * a(ne)) * hv(ne) * c112 * b(ne)
       end do
    case(-3)
       do ne = 1, nemax
          x(ne,1) = (-4.d0 * a(ne-1) +        a(ne)) * c16 * b(ne-1)
          x(ne,2) = (        a(ne-1) + 2.d0 * a(ne)) * c16 * b(ne)
          x(ne,3) = (-2.d0 * a(ne-1) -        a(ne)) * c16 * b(ne-1)
          x(ne,4) = (-       a(ne-1) + 4.d0 * a(ne)) * c16 * b(ne)
       end do
    case(-4)
       do ne = 1, nemax
          x(ne,1) = -0.5d0 * a(ne-1)
          x(ne,2) =  0.5d0 * a(ne)
          x(ne,3) = x(ne,1)
          x(ne,4) = x(ne,2)
       end do
    case(-8)
       do ne = 1, nemax
          x(ne,1) =-0.5d0 * a(ne-1)
          x(ne,2) =-0.5d0 * a(ne)
          x(ne,3) = 0.5d0 * a(ne-1)
          x(ne,4) = 0.5d0 * a(ne)
       end do
    case(-9)
       do ne = 1, nemax
          x(ne,1) = (-2.d0 * a(ne-1) -        a(ne)) * c16 * b(ne-1)
          x(ne,2) = (-       a(ne-1) - 2.d0 * a(ne)) * c16 * b(ne)
          x(ne,3) = ( 2.d0 * a(ne-1) +        a(ne)) * c16 * b(ne-1)
          x(ne,4) = (        a(ne-1) + 2.d0 * a(ne)) * c16 * b(ne)
       end do
    case(-12)
       do ne = 1, nemax
          x(ne,1) = 0.5d0 * ( a(ne-1) + a(ne) ) * b(ne-1) / hv(ne)
          x(ne,2) =-0.5d0 * ( a(ne-1) + a(ne) ) * b(ne  ) / hv(ne)
          x(ne,3) =-0.5d0 * ( a(ne-1) + a(ne) ) * b(ne-1) / hv(ne)
          x(ne,4) = 0.5d0 * ( a(ne-1) + a(ne) ) * b(ne  ) / hv(ne)
       end do
    case(-20)
       do ne = 1, nemax
          x(ne,1) = (12.d0*a(ne-1)*b(ne-1) + 3.d0*a(ne)*b(ne-1) &
               &    + 3.d0*a(ne-1)*b(ne)   + 2.d0*a(ne)*b(ne)) * hv(ne) * c160 * c(ne-1)
          x(ne,2) = ( 3.d0*a(ne-1)*b(ne-1) + 2.d0*a(ne)*b(ne-1) &
               &    + 2.d0*a(ne-1)*b(ne)   + 3.d0*a(ne)*b(ne)) * hv(ne) * c160 * c(ne)
          x(ne,3) = ( 3.d0*a(ne-1)*b(ne-1) + 2.d0*a(ne)*b(ne-1) &
               &    + 2.d0*a(ne-1)*b(ne)   + 3.d0*a(ne)*b(ne)) * hv(ne) * c160 * c(ne-1)
          x(ne,4) = ( 2.d0*a(ne-1)*b(ne-1) + 3.d0*a(ne)*b(ne-1) &
               &    + 3.d0*a(ne-1)*b(ne)   +12.d0*a(ne)*b(ne)) * hv(ne) * c160 * c(ne)
       end do
    case(-22)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne)
          b1 = b(ne-1) ; b2 = b(ne)
          x(ne,1) = (-3.d0*a1*b1 -      a1*b2 -      a2*b1 -      a2*b2) * c112 * c(ne-1)
          x(ne,2) = ( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) * c112 * c(ne)
          x(ne,3) = (-     a1*b1 -      a1*b2 -      a2*b1 - 3.d0*a2*b2) * c112 * c(ne-1)
          x(ne,4) = (      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) * c112 * c(ne)
       end do
    case(-23)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne)
          c1 = c(ne-1) ; c2 = c(ne)
          x(ne,1) =-( 3.d0*a1*b1 + a1*b2 + a2*b1 +      a2*b2) * c112 * c1
          x(ne,2) =-(      a1*b1 + a1*b2 + a2*b1 + 3.d0*a2*b2) * c112 * c2
          x(ne,3) = ( 3.d0*a1*b1 + a1*b2 + a2*b1 +      a2*b2) * c112 * c1
          x(ne,4) = (      a1*b1 + a1*b2 + a2*b1 + 3.d0*a2*b2) * c112 * c2
       end do
    case(-122)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne)
          b1 = b(ne-1) ; b2 = b(ne)
          c1 = c(ne-1) ; c2 = c(ne)
          x(ne,1) = (-3.d0*a1*b1 - a1*b2 - a2*b1 -      a2*b2) * c112 &
               &  + ( 2.d0*a1*b1 + a1*b2 + a2*b1 + 2.d0*a2*b2) * c16 * cp * c1
          x(ne,2) = ( 3.d0*a1*b1 + a1*b2 + a2*b1 +      a2*b2) * c112 &
               &  - ( 2.d0*a1*b1 + a1*b2 + a2*b1 + 2.d0*a2*b2) * c16 * cp * c2
          x(ne,3) = (-     a1*b1 - a1*b2 - a2*b1 - 3.d0*a2*b2) * c112 &
               &  - ( 2.d0*a1*b1 + a1*b2 + a2*b1 + 2.d0*a2*b2) * c16 * cp * c1
          x(ne,4) = (      a1*b1 + a1*b2 + a2*b1 + 3.d0*a2*b2) * c112 &
               &  + ( 2.d0*a1*b1 + a1*b2 + a2*b1 + 2.d0*a2*b2) * c16 * cp * c2
       end do
    case(0)
       !  for SUPG
       do ne = 1, nemax
          x(ne,1) = hv(ne) * cp
          x(ne,2) = x(ne,1)
          x(ne,3) = x(ne,1)
          x(ne,4) = x(ne,1)
       end do
    case(1)
       do ne = 1, nemax
          x(ne,1) = hv(ne) * c13
          x(ne,2) = hv(ne) * c16
          x(ne,3) = x(ne,2)
          x(ne,4) = x(ne,1)
       end do
    case(2)
       do ne = 1, nemax
          x(ne,1) = ( 3.d0 * a(ne-1) +        a(ne)) * hv(ne) * c112
          x(ne,2) = (        a(ne-1) +        a(ne)) * hv(ne) * c112
          x(ne,3) = x(ne,2)
          x(ne,4) = (        a(ne-1) + 3.d0 * a(ne)) * hv(ne) * c112
       end do
    case(3)
       do ne = 1, nemax
          x(ne,1) = (-4.d0 * a(ne-1) +        a(ne)) * c16
          x(ne,2) = (        a(ne-1) + 2.d0 * a(ne)) * c16
          x(ne,3) = (-2.d0 * a(ne-1) -        a(ne)) * c16
          x(ne,4) = (-       a(ne-1) + 4.d0 * a(ne)) * c16
       end do
    case(4)
       do ne = 1, nemax
          x(ne,1) = -0.5d0
          x(ne,2) =  0.5d0
          x(ne,3) = x(ne,1)
          x(ne,4) = x(ne,2)
       end do
    case(5)
       do ne = 1, nemax
          x(ne,1) = (-2.d0 * a(ne-1) -        a(ne)) * c16
          x(ne,2) = ( 2.d0 * a(ne-1) +        a(ne)) * c16
          x(ne,3) = (-       a(ne-1) - 2.d0 * a(ne)) * c16
          x(ne,4) = (        a(ne-1) + 2.d0 * a(ne)) * c16
       end do
    case(6)
       do ne = 1, nemax
          x(ne,1) = (- a(ne-1) + a(ne)) * c13
          x(ne,2) = (- a(ne-1) + a(ne)) * c16
          x(ne,3) = x(ne,2)
          x(ne,4) = x(ne,1)
       end do
    case(8)
       do ne = 1, nemax
          x(ne,1) =-0.5d0
          x(ne,2) = x(ne,1)
          x(ne,3) = 0.5d0
          x(ne,4) = x(ne,3)
       end do
    case(9)
       do ne = 1, nemax
          x(ne,1) = (-2.d0 * a(ne-1) -        a(ne)) * c16
          x(ne,2) = (-       a(ne-1) - 2.d0 * a(ne)) * c16
          x(ne,3) =-x(ne,1)
          x(ne,4) =-x(ne,2)
       end do
    case(10)
       do ne = 1, nemax
          x(ne,1) = a(ne-1) / hv(ne)
          x(ne,2) =-a(ne)   / hv(ne)
          x(ne,3) =-x(ne,1)
          x(ne,4) =-x(ne,2)
       end do
    case(11)
       do ne = 1, nemax
          x(ne,1) = 1.d0 / hv(ne)
          x(ne,2) =-x(ne,1)
          x(ne,3) = x(ne,2)
          x(ne,4) = x(ne,1)
       end do
    case(12)
       do ne = 1, nemax
          x(ne,1) = 0.5d0 * ( a(ne-1) + a(ne) ) / hv(ne)
          x(ne,2) =-x(ne,1)
          x(ne,3) = x(ne,2)
          x(ne,4) = x(ne,1)
       end do
    case(13)
       do ne = 1, nemax
          x(ne,1) = 0.5d0 * ( a(ne-1) - a(ne) ) / hv(ne)
          x(ne,2) = x(ne,1)
          x(ne,3) =-x(ne,1)
          x(ne,4) = x(ne,3)
       end do
    case(20)
       do ne = 1, nemax
          x(ne,1) = (12.d0*a(ne-1)*b(ne-1) + 3.d0*a(ne)*b(ne-1) &
               &    + 3.d0*a(ne-1)*b(ne)   + 2.d0*a(ne)*b(ne)) * hv(ne) * c160
          x(ne,2) = ( 3.d0*a(ne-1)*b(ne-1) + 2.d0*a(ne)*b(ne-1) &
               &    + 2.d0*a(ne-1)*b(ne)   + 3.d0*a(ne)*b(ne)) * hv(ne) * c160
          x(ne,3) = x(ne,2)
          x(ne,4) = ( 2.d0*a(ne-1)*b(ne-1) + 3.d0*a(ne)*b(ne-1) &
               &    + 3.d0*a(ne-1)*b(ne)   +12.d0*a(ne)*b(ne)) * hv(ne) * c160
       end do
    case(21)
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
    case(22)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne)
          b1 = b(ne-1) ; b2 = b(ne)
          x(ne,1) = (-3.d0*a1*b1 -      a1*b2 -      a2*b1 -      a2*b2) * c112
          x(ne,2) = ( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) * c112
          x(ne,3) = (-     a1*b1 -      a1*b2 -      a2*b1 - 3.d0*a2*b2) * c112
          x(ne,4) = (      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) * c112
       end do
    case(23)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne)
          b1 = b(ne-1) ; b2 = b(ne)
          x(ne,1) =-( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) * c112
          x(ne,2) =-(      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) * c112
          x(ne,3) =-x(ne,1)
          x(ne,4) =-x(ne,2)
       end do
    case(24)
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
    case(25)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne)
          b1 = b(ne-1) ; b2 = b(ne) ; hvlinv = 1.d0 / hv(ne)
          x(ne,1) = (2.d0*a1+     a2)*(b1-b2) * hvlinv * c16
          x(ne,2) =-(2.d0*a1+     a2)*(b1-b2) * hvlinv * c16
          x(ne,3) = (     a1+2.d0*a2)*(b1-b2) * hvlinv * c16
          x(ne,4) =-(     a1+2.d0*a2)*(b1-b2) * hvlinv * c16
       end do
    case(26)
       do ne = 1, nemax
          x(ne,1) = ( 2.d0*a(ne-1)*b(ne-1) - 2.d0*a(ne-1)*b(ne) &
               &     +     a(ne)  *b(ne-1) -      a(ne)  *b(ne)) / hv(ne) * c16
          x(ne,2) = (      a(ne-1)*b(ne-1) -      a(ne-1)*b(ne) &
               &     +2.d0*a(ne)  *b(ne-1) - 2.d0*a(ne)  *b(ne)) / hv(ne) * c16
          x(ne,3) =-x(ne,1)
          x(ne,4) =-x(ne,2)
       end do
    case(27)
       do ne = 1, nemax
          x(ne,1) = ( 2.d0*a(ne-1)*b(ne-1) +      a(ne-1)*b(ne) &
               &     +     a(ne)  *b(ne-1) + 2.d0*a(ne)  *b(ne)) / hv(ne) * c16
          x(ne,2) =-x(ne,1)
          x(ne,3) = x(ne,2)
          x(ne,4) = x(ne,1)
       end do
    case(30)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne) ; c1 = c(ne-1) ; c2 = c(ne)
          coef1 =      b1*(2.d0*c1+     c2)+     b2*(     c1+     c2)
          coef2 =      b1*(     c1+     c2)+     b2*(     c1+2.d0*c2)
          coef3 = 2.d0*b1*(5.d0*c1+     c2)+     b2*(2.d0*c1+     c2)
          coef4 =      b1*(     c1+2.d0*c2)+2.d0*b2*(     c1+5.d0*c2)
          x(ne,1) = (a1*coef3+a2*coef1)*hv(ne)*c160
          x(ne,2) = (a1*coef1+a2*coef2)*hv(ne)*c160
          x(ne,3) = x(ne,2)
          x(ne,4) = (a1*coef2+a2*coef4)*hv(ne)*c160
       end do
    case(31)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne) ; c1 = c(ne-1) ; c2 = c(ne)
          x(ne,1) = (12.d0*a1*b1+3.d0*a1*b2+3.d0*a2*b1+ 2.d0*a2*b2)*(c2-c1)*c160
          x(ne,2) = ( 3.d0*a1*b1+2.d0*a1*b2+2.d0*a2*b1+ 3.d0*a2*b2)*(c2-c1)*c160
          x(ne,3) = x(ne,2)
          x(ne,4) = ( 2.d0*a1*b1+3.d0*a1*b2+3.d0*a2*b1+12.d0*a2*b2)*(c2-c1)*c160
       end do
    case(32)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne) ; c1 = c(ne-1) ; c2 = c(ne)
          coef1 = 3.d0*b1*(4.d0*c1+c2)+b2*(3.d0*c1+2.d0*c2)
          coef2 = 3.d0*b1*c1+2.d0*b2*c1+2.d0*b1*c2+3.d0*b2*c2
          coef3 = b1*(2.d0*c1+3.d0*c2)+3.d0*b2*(c1+4.d0*c2)
          x(ne,1) =-(a1*coef1+a2*coef2)*c160
          x(ne,2) =-x(ne,1)
          x(ne,3) =-(a1*coef2+a2*coef3)*c160
          x(ne,4) =-x(ne,3)
       end do
    case(33)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne) ; c1 = c(ne-1) ; c2 = c(ne)
          coef1 =      b1*(3.d0*c1+2.d0*c2)+     b2*(2.d0*c1+3.d0*c2)
          coef2 = 3.d0*b1*(4.d0*c1+     c2)+     b2*(3.d0*c1+2.d0*c2)
          coef3 =      b1*(2.d0*c1+3.d0*c2)+3.d0*b2*(     c1+4.d0*c2)
          x(ne,1) =-(a1*coef2+a2*coef1)*c160
          x(ne,2) =-(a1*coef1+a2*coef3)*c160
          x(ne,3) = (a1*coef2+a2*coef1)*c160
          x(ne,4) = (a1*coef1+a2*coef3)*c160
       end do
    case(36)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne) ; c1 = c(ne-1) ; c2 = c(ne)
          hvlinv = 1.d0 / hv(ne)
          x(ne,1) = (3.d0*a1*b1+a1*b2+a2*b1+     a2*b2)*(c1-c2)*hvlinv*c112
          x(ne,2) = (     a1*b1+a1*b2+a2*b1+3.d0*a2*b2)*(c1-c2)*hvlinv*c112 
          x(ne,3) = (3.d0*a1*b1+a1*b2+a2*b1+     a2*b2)*(c2-c1)*hvlinv*c112
          x(ne,4) = (     a1*b1+a1*b2+a2*b1+3.d0*a2*b2)*(c2-c1)*hvlinv*c112 
       end do
    case(37)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne) ; c1 = c(ne-1) ; c2 = c(ne)
          hvlinv = 1.d0 / hv(ne)
          x(ne,1) = (a2*(b1+b2)*c1+a1*(3.d0*b1+b2)*c1+a1*(b1+b2)*c2+a2*(b1+3.d0*b2)*c2)*hvlinv*c112
          x(ne,2) =-x(ne,1)
          x(ne,3) = x(ne,2)
          x(ne,4) = x(ne,1)
       end do
    case(38)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne) ; c1 = c(ne-1) ; c2 = c(ne)
          hvlinv = 1.d0 / hv(ne)
          x(ne,1) =( (a2*(b1+b2)*c1+a1*(3.d0*b1+b2)*c1+a1*(b1+b2)*c2+a2*(b1+3.d0*b2)*c2) &
               &    +(3.d0*a1*b1+a1*b2+a2*b1+     a2*b2)*(c1-c2))*hvlinv*c112
          x(ne,2) =(-(a2*(b1+b2)*c1+a1*(3.d0*b1+b2)*c1+a1*(b1+b2)*c2+a2*(b1+3.d0*b2)*c2) &
               &    +(     a1*b1+a1*b2+a2*b1+3.d0*a2*b2)*(c1-c2))*hvlinv*c112 
          x(ne,3) =(-(a2*(b1+b2)*c1+a1*(3.d0*b1+b2)*c1+a1*(b1+b2)*c2+a2*(b1+3.d0*b2)*c2) &
               &    +(3.d0*a1*b1+a1*b2+a2*b1+     a2*b2)*(c2-c1))*hvlinv*c112
          x(ne,4) =( (a2*(b1+b2)*c1+a1*(3.d0*b1+b2)*c1+a1*(b1+b2)*c2+a2*(b1+3.d0*b2)*c2) &
               &    +(     a1*b1+a1*b2+a2*b1+3.d0*a2*b2)*(c2-c1))*hvlinv*c112 
       end do
    ! --- SUPG ---
    case(-101) ! (-1)+(-8)*(0)
       do ne = 1, nemax
          x(ne,1) = hv(ne) * c13 * a(ne-1) - 0.5d0 * a(ne-1) * hv(ne) * cp
          x(ne,2) = hv(ne) * c16 * a(ne)   - 0.5d0 * a(ne)   * hv(ne) * cp
          x(ne,3) = hv(ne) * c16 * a(ne-1) + 0.5d0 * a(ne-1) * hv(ne) * cp
          x(ne,4) = hv(ne) * c13 * a(ne)   + 0.5d0 * a(ne)   * hv(ne) * cp
       end do
    case(-102) ! (-2)+(-9)*(0)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne)
          b1 = b(ne-1) ; b2 = b(ne)
          x(ne,1) = ( 3.d0 * a1 +        a2) * c112 * b1 * hv(ne) &
               &  + (-2.d0 * a1 -        a2) * c16  * b1 * hv(ne) * cp
          x(ne,2) = (        a1 +        a2) * c112 * b2 * hv(ne) &
               &  + (-       a1 - 2.d0 * a2) * c16  * b2 * hv(ne) * cp
          x(ne,3) = (        a1 +        a2) * c112 * b1 * hv(ne) &
               &  + ( 2.d0 * a1 +        a2) * c16  * b1 * hv(ne) * cp
          x(ne,4) = (        a1 + 3.d0 * a2) * c112 * b2 * hv(ne) &
               &  + (        a1 + 2.d0 * a2) * c16  * b2 * hv(ne) * cp
       end do
    case(-120) ! (-20)+(-23)*(0)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne)
          c1 = c(ne-1) ; c2 = c(ne)
          x(ne,1) = (12.d0*a1*b1 + 3.d0*a2*b1 + 3.d0*a1*b2 + 2.d0*a2*b2) * c160 * c1 * hv(ne) &
               &  - ( 3.d0*a1*b1 +      a2*b1 +      a1*b2 +      a2*b2) * c112 * c1 * hv(ne) * cp
          x(ne,2) = ( 3.d0*a1*b1 + 2.d0*a2*b1 + 2.d0*a1*b2 + 3.d0*a2*b2) * c160 * c2 * hv(ne) &
               &  - (      a1*b1 +      a2*b1 +      a1*b2 + 3.d0*a2*b2) * c112 * c2 * hv(ne) * cp
          x(ne,3) = ( 3.d0*a1*b1 + 2.d0*a2*b1 + 2.d0*a1*b2 + 3.d0*a2*b2) * c160 * c1 * hv(ne) &
               &  + ( 3.d0*a1*b1 +      a2*b1 +      a1*b2 +      a2*b2) * c112 * c1 * hv(ne) * cp
          x(ne,4) = ( 2.d0*a1*b1 + 3.d0*a2*b1 + 3.d0*a1*b2 +12.d0*a2*b2) * c160 * c2 * hv(ne) &
               &  + (      a1*b1 +      a2*b1 +      a1*b2 + 3.d0*a2*b2) * c112 * c2 * hv(ne) * cp
       end do
    case(101) ! (1)+(8)*(0)
       do ne = 1, nemax
          x(ne,1) = hv(ne) * c13 - 0.5d0 * hv(ne) * cp
          x(ne,2) = hv(ne) * c16 - 0.5d0 * hv(ne) * cp
          x(ne,3) = hv(ne) * c16 + 0.5d0 * hv(ne) * cp
          x(ne,4) = hv(ne) * c13 + 0.5d0 * hv(ne) * cp
       end do
    case(102) ! (2)+(9)*(0)
       do ne = 1, nemax
          x(ne,1) = ( 3.d0 * a(ne-1) +        a(ne)) * c112 * hv(ne) &
               &  + (-2.d0 * a(ne-1) -        a(ne)) * c16  * hv(ne) * cp
          x(ne,2) = (        a(ne-1) +        a(ne)) * c112 * hv(ne) &
               &  + (-       a(ne-1) - 2.d0 * a(ne)) * c16  * hv(ne) * cp
          x(ne,3) = (        a(ne-1) +        a(ne)) * c112 * hv(ne) &
               &  + ( 2.d0 * a(ne-1) +        a(ne)) * c16  * hv(ne) * cp
          x(ne,4) = (        a(ne-1) + 3.d0 * a(ne)) * c112 * hv(ne) &
               &  + (        a(ne-1) + 2.d0 * a(ne)) * c16  * hv(ne) * cp
       end do
    case(103) ! (3)+(10)*(0)
       do ne = 1, nemax
          x(ne,1) = (-4.d0 * a(ne-1) +        a(ne)) * c16 + a(ne-1) * cp
          x(ne,2) = (        a(ne-1) + 2.d0 * a(ne)) * c16 - a(ne)   * cp 
          x(ne,3) = (-2.d0 * a(ne-1) -        a(ne)) * c16 - a(ne-1) * cp
          x(ne,4) = (-       a(ne-1) + 4.d0 * a(ne)) * c16 + a(ne)   * cp
       end do
    case(105) ! (5)+(12)*(0)
      do ne = 1, nemax
          x(ne,1) = (-2.d0 * a(ne-1) -        a(ne)) * c16 + 0.5d0 * ( a(ne-1) + a(ne) ) * cp
          x(ne,2) = ( 2.d0 * a(ne-1) +        a(ne)) * c16 - 0.5d0 * ( a(ne-1) + a(ne) ) * cp
          x(ne,3) = (-       a(ne-1) - 2.d0 * a(ne)) * c16 - 0.5d0 * ( a(ne-1) + a(ne) ) * cp
          x(ne,4) = (        a(ne-1) + 2.d0 * a(ne)) * c16 + 0.5d0 * ( a(ne-1) + a(ne) ) * cp
       end do
    case(106) ! (6)+(13)*(0)
       do ne = 1, nemax
          x(ne,1) = (- a(ne-1) + a(ne)) * c13 + 0.5d0 * ( a(ne-1) - a(ne) ) * cp
          x(ne,2) = (- a(ne-1) + a(ne)) * c16 + 0.5d0 * ( a(ne-1) - a(ne) ) * cp
          x(ne,3) = (- a(ne-1) + a(ne)) * c16 - 0.5d0 * ( a(ne-1) - a(ne) ) * cp
          x(ne,4) = (- a(ne-1) + a(ne)) * c13 - 0.5d0 * ( a(ne-1) - a(ne) ) * cp
       end do
    case(120) ! (20)+(23)*(0)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne)
          b1 = b(ne-1) ; b2 = b(ne)
          x(ne,1) = (12.d0*a1*b1 + 3.d0*a2*b1 + 3.d0*a1*b2 + 2.d0*a2*b2) * c160 * hv(ne) &
               &  - ( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) * c112 * hv(ne) * cp
          x(ne,2) = ( 3.d0*a1*b1 + 2.d0*a2*b1 + 2.d0*a1*b2 + 3.d0*a2*b2) * c160 * hv(ne) &
               &  - (      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) * c112 * hv(ne) * cp
          x(ne,3) = ( 3.d0*a1*b1 + 2.d0*a2*b1 + 2.d0*a1*b2 + 3.d0*a2*b2) * c160 * hv(ne) &
               &  + ( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) * c112 * hv(ne) * cp
          x(ne,4) = ( 2.d0*a1*b1 + 3.d0*a2*b1 + 3.d0*a1*b2 +12.d0*a2*b2) * c160 * hv(ne) &
               &  + (      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) * c112 * hv(ne) * cp
       end do
    case(121) ! (21)+(26)*(0)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne)
          b1 = b(ne-1) ; b2 = b(ne)
          x(ne,1) = (-3.d0*a1*b1 + 3.d0*a1*b2 -     a2  *b1 +      a2  *b2) * c112 &
               &  + ( 2.d0*a1*b1 - 2.d0*a1*b2 +     a2  *b1 -      a2  *b2) * c16 * cp
          x(ne,2) = (-     a1*b1 +      a1*b2 -     a2  *b1 +      a2  *b2) * c112 &
               &  + (      a1*b1 -      a1*b2 +2.d0*a2  *b1 - 2.d0*a2  *b2) * c16 * cp
          x(ne,3) = (-     a1*b1 +      a1*b2 -     a2  *b1 +      a2  *b2) * c112 &
               &  - ( 2.d0*a1*b1 - 2.d0*a1*b2 +     a2  *b1 -      a2  *b2) * c16 * cp
          x(ne,4) = (-     a1*b1 +      a1*b2 -3.d0*a2  *b1 + 3.d0*a2  *b2) * c112 &
               &  - (      a1*b1 -      a1*b2 +2.d0*a2  *b1 - 2.d0*a2  *b2) * c16 * cp
       end do
    case(122) ! (22)+(27)*(0)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne)
          b1 = b(ne-1) ; b2 = b(ne)
          x(ne,1) = (-3.d0*a1*b1 - a1*b2 - a2*b1 -      a2*b2) * c112 &
               &  + ( 2.d0*a1*b1 + a1*b2 + a2*b1 + 2.d0*a2*b2) * c16 * cp
          x(ne,2) = ( 3.d0*a1*b1 + a1*b2 + a2*b1 +      a2*b2) * c112 &
               &  - ( 2.d0*a1*b1 + a1*b2 + a2*b1 + 2.d0*a2*b2) * c16 * cp
          x(ne,3) = (-     a1*b1 - a1*b2 - a2*b1 - 3.d0*a2*b2) * c112 &
               &  - ( 2.d0*a1*b1 + a1*b2 + a2*b1 + 2.d0*a2*b2) * c16 * cp
          x(ne,4) = (      a1*b1 + a1*b2 + a2*b1 + 3.d0*a2*b2) * c112 &
               &  + ( 2.d0*a1*b1 + a1*b2 + a2*b1 + 2.d0*a2*b2) * c16 * cp
       end do
    case(130) ! (30)+(33)*(0)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne) ; c1 = c(ne-1) ; c2 = c(ne)
          coef1 =      b1*(2.d0*c1+     c2)+     b2*(     c1+     c2)
          coef2 =      b1*(     c1+     c2)+     b2*(     c1+2.d0*c2)
          coef3 = 2.d0*b1*(5.d0*c1+     c2)+     b2*(2.d0*c1+     c2)
          coef4 =      b1*(     c1+2.d0*c2)+2.d0*b2*(     c1+5.d0*c2)
          coef5 =      b1*(3.d0*c1+2.d0*c2)+     b2*(2.d0*c1+3.d0*c2)
          coef6 = 3.d0*b1*(4.d0*c1+     c2)+     b2*(3.d0*c1+2.d0*c2)
          coef7 =      b1*(2.d0*c1+3.d0*c2)+3.d0*b2*(     c1+4.d0*c2)
          x(ne,1) = (a1 * coef3 + a2 * coef1) * hv(ne) * c160 &
               &  - (a1 * coef6 + a2 * coef5) * hv(ne) * c160 * cp
          x(ne,2) = (a1 * coef1 + a2 * coef2) * hv(ne) * c160 &
               &  - (a1 * coef5 + a2 * coef7) * hv(ne) * c160 * cp
          x(ne,3) = (a1 * coef1 + a2 * coef2) * hv(ne) * c160 &
               &  + (a1 * coef6 + a2 * coef5) * hv(ne) * c160 * cp
          x(ne,4) = (a1 * coef2 + a2 * coef4) * hv(ne) * c160 &
               &  + (a1 * coef5 + a2 * coef7) * hv(ne) * c160 * cp
       end do
    case(131) ! (31)+(36)*(0)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne) ; c1 = c(ne-1) ; c2 = c(ne)
          x(ne,1) = (12.d0*a1*b1+3.d0*a1*b2+3.d0*a2*b1+ 2.d0*a2*b2)*(c2-c1) * c160 &
               &  - ( 3.d0*a1*b1+     a1*b2+     a2*b1+      a2*b2)*(c2-c1) * c112 * cp
          x(ne,2) = ( 3.d0*a1*b1+2.d0*a1*b2+2.d0*a2*b1+ 3.d0*a2*b2)*(c2-c1) * c160 &
               &  - (      a1*b1+     a1*b2+     a2*b1+ 3.d0*a2*b2)*(c2-c1) * c112 * cp
          x(ne,3) = ( 3.d0*a1*b1+2.d0*a1*b2+2.d0*a2*b1+ 3.d0*a2*b2)*(c2-c1) * c160 &
               &  + ( 3.d0*a1*b1+     a1*b2+     a2*b1+      a2*b2)*(c2-c1) * c112 * cp
          x(ne,4) = ( 2.d0*a1*b1+3.d0*a1*b2+3.d0*a2*b1+12.d0*a2*b2)*(c2-c1) * c160 &
               &  + (      a1*b1+     a1*b2+     a2*b1+ 3.d0*a2*b2)*(c2-c1) * c112 * cp
       end do
    case(132) ! (32)+(37)*(0)
       do ne = 1, nemax
          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne) ; c1 = c(ne-1) ; c2 = c(ne)
          coef1 = 3.d0*b1*(4.d0*c1+c2)+b2*(3.d0*c1+2.d0*c2)
          coef2 = 3.d0*b1*c1+2.d0*b2*c1+2.d0*b1*c2+3.d0*b2*c2
          coef3 = b1*(2.d0*c1+3.d0*c2)+3.d0*b2*(c1+4.d0*c2)
          x(ne,1) =-(a1*coef1+a2*coef2) * c160 &
               &  + (a2*(b1+b2)*c1+a1*(3.d0*b1+b2)*c1+a1*(b1+b2)*c2+a2*(b1+3.d0*b2)*c2) * c112 * cp
          x(ne,2) = (a1*coef1+a2*coef2) * c160 &
               &  - (a2*(b1+b2)*c1+a1*(3.d0*b1+b2)*c1+a1*(b1+b2)*c2+a2*(b1+3.d0*b2)*c2) * c112 * cp
          x(ne,3) =-(a1*coef2+a2*coef3)*c160 &
               &  - (a2*(b1+b2)*c1+a1*(3.d0*b1+b2)*c1+a1*(b1+b2)*c2+a2*(b1+3.d0*b2)*c2) * c112 * cp
          x(ne,4) = (a1*coef2+a2*coef3)*c160 &
               &  + (a2*(b1+b2)*c1+a1*(3.d0*b1+b2)*c1+a1*(b1+b2)*c2+a2*(b1+3.d0*b2)*c2) * c112 * cp
       end do
!!$    case(15)
!!$       do ne = 1, nemax
!!$          x(ne,1) = ( 12.d0*vv(ne-1)*a(ne-1) + 3.d0*vv(ne)*a(ne-1) &
!!$               &     + 3.d0*vv(ne-1)*a(ne)   + 2.d0*vv(ne)*a(ne)) * hv(ne) * c160
!!$          x(ne,2) = (  3.d0*vv(ne-1)*a(ne-1) + 2.d0*vv(ne)*a(ne-1) &
!!$               &     + 2.d0*vv(ne-1)*a(ne)   + 3.d0*vv(ne)*a(ne)) * hv(ne) * c160
!!$          x(ne,3) = x(ne,2)
!!$          x(ne,4) = (  2.d0*vv(ne-1)*a(ne-1) + 3.d0*vv(ne)*a(ne-1) &
!!$               &     + 3.d0*vv(ne-1)*a(ne)   +12.d0*vv(ne)*a(ne)) * hv(ne) * c160
!!$       end do
!!$    case(16)
!!$       do ne = 1, nemax
!!$          x(ne,1) = (3.d0 * vv(ne-1) +        vv(ne)) * (-a(ne-1) + a(ne)) * c112
!!$          x(ne,2) = (       vv(ne-1) +        vv(ne)) * (-a(ne-1) + a(ne)) * c112
!!$          x(ne,3) = x(ne,2)
!!$          x(ne,4) = (       vv(ne-1) + 3.d0 * vv(ne)) * (-a(ne-1) + a(ne)) * c112
!!$       end do
!!$    case(17)
!!$       do ne = 1, nemax
!!$          x(ne,1) = (-3.d0*vv(ne-1)*a(ne-1) -      vv(ne)*a(ne-1) &
!!$               &     -     vv(ne-1)*a(ne)   -      vv(ne)*a(ne)) * c112
!!$          x(ne,2) =-x(ne,1)
!!$          x(ne,3) = (-     vv(ne-1)*a(ne-1) -      vv(ne)*a(ne-1) &
!!$               &     -     vv(ne-1)*a(ne)   - 3.d0*vv(ne)*a(ne)) * c112
!!$          x(ne,4) =-x(ne,3)
!!$       end do
!!$    case(18)
!!$       do ne = 1, nemax
!!$          x(ne,1) = ( 2.d0*vv(ne-1)*a(ne-1) +      vv(ne)*a(ne-1) &
!!$               &     +     vv(ne-1)*a(ne)   + 2.d0*vv(ne)*a(ne))  / hv(ne) * c16
!!$          x(ne,2) =-x(ne,1)
!!$          x(ne,3) = x(ne,2)
!!$          x(ne,4) =-x(ne,3)
!!$       end do
!!$    case(20)
!!$       do ne = 1, nemax
!!$          x(ne,1) = ( 3.d0*vv(ne-1)*a(ne-1) +      vv(ne)*a(ne-1) &
!!$               &     +     vv(ne-1)*a(ne)   +      vv(ne)*a(ne)) &
!!$               &  * (b(ne-1) - b(ne)) / (12.d0 * hv(ne))
!!$          x(ne,2) = (      vv(ne-1)*a(ne-1) +      vv(ne)*a(ne-1) &
!!$               &     +     vv(ne-1)*a(ne)   + 3.d0*vv(ne)*a(ne)) &
!!$               &  * (b(ne-1) - b(ne)) / (12.d0 * hv(ne))
!!$          x(ne,3) =-x(ne,1)
!!$          x(ne,4) =-x(ne,2)
!!$       end do
!!$    case(36)
!!$       do ne = 1, nemax
!!$          x(ne,1) = ( 12.d0*vv(ne-1)*b(ne-1) + 3.d0*vv(ne)*b(ne-1) &
!!$               &     + 3.d0*vv(ne-1)*b(ne)   + 2.d0*vv(ne)*b(ne)) * (-a(ne-1)+a(ne)) * c160
!!$          x(ne,2) = (  3.d0*vv(ne-1)*b(ne-1) + 2.d0*vv(ne)*b(ne-1) &
!!$               &     + 2.d0*vv(ne-1)*b(ne)   + 3.d0*vv(ne)*b(ne)) * (-a(ne-1)+a(ne)) * c160
!!$          x(ne,3) = x(ne,2)
!!$          x(ne,4) = (  2.d0*vv(ne-1)*b(ne-1) + 3.d0*vv(ne)*b(ne-1) &
!!$               &     + 3.d0*vv(ne-1)*b(ne)   +12.d0*vv(ne)*b(ne)) * (-a(ne-1)+a(ne)) * c160
!!$       end do
!!$    case(37)
!!$       do ne = 1, nemax
!!$          x(ne,1) = (12.d0*vv(ne-1)*a(ne-1) + 3.d0*vv(ne)*a(ne-1) &
!!$               &    + 3.d0*vv(ne-1)*a(ne)   + 2.d0*vv(ne)*a(ne)) * (-b(ne-1)+b(ne)) * c160
!!$          x(ne,2) = ( 3.d0*vv(ne-1)*a(ne-1) + 2.d0*vv(ne)*a(ne-1) &
!!$               &    + 2.d0*vv(ne-1)*a(ne)   + 3.d0*vv(ne)*a(ne)) * (-b(ne-1)+b(ne)) * c160
!!$          x(ne,3) = ( 3.d0*vv(ne-1)*a(ne-1) + 2.d0*vv(ne)*a(ne-1) &
!!$               &    + 2.d0*vv(ne-1)*a(ne)   + 3.d0*vv(ne)*a(ne)) * (-b(ne-1)+b(ne)) * c160
!!$          x(ne,4) = ( 2.d0*vv(ne-1)*a(ne-1) + 3.d0*vv(ne)*a(ne-1) &
!!$               &    + 3.d0*vv(ne-1)*a(ne)   +12.d0*vv(ne)*a(ne)) * (-b(ne-1)+b(ne)) * c160
!!$       end do
!!$    case(38)
!!$       do ne = 1, nemax
!!$          a1 = a(ne-1) ; b1 = b(ne-1) ; p1 = vv(ne-1)
!!$          a2 = a(ne)   ; b2 = b(ne)   ; p2 = vv(ne)
!!$          x(ne,1) =-( 12.d0*p1*a1*b1 + 3.d0*p2*a1*b1 + 3.d0*p1*a2*b1 + 2.d0*p2*a2*b1 &
!!$               &     + 3.d0*p1*a1*b2 + 2.d0*p2*a1*b2 + 2.d0*p1*a2*b2 + 3.d0*p2*a2*b2) * c160
!!$          x(ne,2) = ( 12.d0*p1*a1*b1 + 3.d0*p2*a1*b1 + 3.d0*p1*a2*b1 + 2.d0*p2*a2*b1 &
!!$               &     + 3.d0*p1*a1*b2 + 2.d0*p2*a1*b2 + 2.d0*p1*a2*b2 + 3.d0*p2*a2*b2) * c160
!!$          x(ne,3) =-(  3.d0*p1*a1*b1 + 2.d0*p2*a1*b1 + 2.d0*p1*a2*b1 + 3.d0*p2*a2*b1 &
!!$               &     + 2.d0*p1*a1*b2 + 3.d0*p2*a1*b2 + 3.d0*p1*a2*b2 +12.d0*p2*a2*b2) * c160
!!$          x(ne,4) = (  3.d0*p1*a1*b1 + 2.d0*p2*a1*b1 + 2.d0*p1*a2*b1 + 3.d0*p2*a2*b1 &
!!$               &     + 2.d0*p1*a1*b2 + 3.d0*p2*a1*b2 + 3.d0*p1*a2*b2 +12.d0*p2*a2*b2) * c160
!!$       end do
!!$    case(39)
!!$       do ne = 1, nemax
!!$          a1 = a(ne-1) ; b1 = b(ne-1) ; p1 = vv(ne-1)
!!$          a2 = a(ne)   ; b2 = b(ne)   ; p2 = vv(ne)
!!$          x(ne,1) = (-48.d0*a1*b1*p1+3.d0*a2*b1*p1+3.d0*a1*b2*p1+ 2.d0*a2*b2*p1 &
!!$               &     + 3.d0*a1*b1*p2+2.d0*a2*b1*p2+2.d0*a1*b2*p2+ 3.d0*a2*b2*p2) * c160
!!$          x(ne,2) = (  3.d0*a1*b1*p1+2.d0*a2*b1*p1+2.d0*a1*b2*p1+ 3.d0*a2*b2*p1 &
!!$               &     + 2.d0*a1*b1*p2+3.d0*a2*b1*p2+3.d0*a1*b2*p2+12.d0*a2*b2*p2) * c160
!!$          x(ne,3) = (-12.d0*a1*b1*p1-3.d0*a2*b1*p1-3.d0*a1*b2*p1- 2.d0*a2*b2*p1 &
!!$               &     - 3.d0*a1*b1*p2-2.d0*a2*b1*p2-2.d0*a1*b2*p2- 3.d0*a2*b2*p2) * c160
!!$          x(ne,4) = (- 3.d0*a1*b1*p1-2.d0*a2*b1*p1-2.d0*a1*b2*p1- 3.d0*a2*b2*p1 &
!!$               &     - 2.d0*a1*b1*p2-3.d0*a2*b1*p2-3.d0*a1*b2*p2+48.d0*a2*b2*p2) * c160
!!$       end do
!!$    case(41)
!!$       do ne = 1, nemax
!!$          x(ne,1) = (b(ne-1)*( 3.d0*vv(ne-1)*a(ne-1)+     vv(ne)*a(ne-1) &
!!$               &              +     vv(ne-1)*a(ne)  +     vv(ne)*a(ne)) &
!!$               &    +b(ne)  *(      vv(ne-1)*a(ne-1)+     vv(ne)*a(ne-1) &
!!$               &              +     vv(ne-1)*a(ne)  +3.d0*vv(ne)*a(ne))) / hv(ne) * c112
!!$          x(ne,2) =-x(ne,1)
!!$          x(ne,3) = x(ne,2)
!!$          x(ne,4) = x(ne,1)
!!$       end do
!!$    case(42)
!!$       do ne = 1, nemax
!!$          x(ne,1) = a(ne-1)*b(ne-1)*vv(ne-1) / hv(ne)
!!$          x(ne,2) =-a(ne)  *b(ne)  *vv(ne)   / hv(ne)
!!$          x(ne,3) =-x(ne,1)
!!$          x(ne,4) =-x(ne,2)
!!$       end do
!!$    case(44)
!!$       do ne = 1, nemax
!!$          x(ne,1) = ( 10.d0*a(ne-1)*b(ne-1)*vv(ne-1)+ 2.d0*a(ne)*b(ne-1)*vv(ne-1) &
!!$               &     + 2.d0*a(ne-1)*b(ne)  *vv(ne-1)+      a(ne)*b(ne)  *vv(ne-1) &
!!$               &     + 2.d0*a(ne-1)*b(ne-1)*vv(ne)  +      a(ne)*b(ne-1)*vv(ne) &
!!$               &     +      a(ne-1)*b(ne)  *vv(ne)  +      a(ne)*b(ne)  *vv(ne)) * hv(ne) * c160
!!$          x(ne,2) = (  2.d0*a(ne-1)*b(ne-1)*vv(ne-1)+      a(ne)*b(ne-1)*vv(ne-1) &
!!$               &     +      a(ne-1)*b(ne)  *vv(ne-1)+      a(ne)*b(ne)  *vv(ne-1) &
!!$               &     +      a(ne-1)*b(ne-1)*vv(ne)  +      a(ne)*b(ne-1)*vv(ne) &
!!$               &     +      a(ne-1)*b(ne)  *vv(ne)  + 2.d0*a(ne)*b(ne)  *vv(ne)) * hv(ne) * c160
!!$          x(ne,3) = x(ne,2)
!!$          x(ne,4) = (       a(ne-1)*b(ne-1)*vv(ne-1)+      a(ne)*b(ne-1)*vv(ne-1) &
!!$               &     +      a(ne-1)*b(ne)  *vv(ne-1)+ 2.d0*a(ne)*b(ne)  *vv(ne-1) &
!!$               &     +      a(ne-1)*b(ne-1)*vv(ne)  + 2.d0*a(ne)*b(ne-1)*vv(ne) &
!!$               &     + 2.d0*a(ne-1)*b(ne)  *vv(ne)  +10.d0*a(ne)*b(ne)  *vv(ne)) * hv(ne) * c160
!!$       end do
!!$    case(45)
!!$       do ne = 1, nemax
!!$          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne)
!!$          p1 = vv(ne-1) ; p2 = vv(ne) ; hvlinv = 1.d0 / hv(ne)
!!$          x(ne,1) = (b1-b2)*( 9.d0*a1*p1-a2*p1-a1*p2-     a2*p2)*hvlinv * c112
!!$          x(ne,2) =-(b1-b2)*(      a1*p1+a2*p1+a1*p2+3.d0*a2*p2)*hvlinv * c112
!!$          x(ne,3) = (b1-b2)*( 3.d0*a1*p1+a2*p1+a1*p2+     a2*p2)*hvlinv * c112
!!$          x(ne,4) =-(b1-b2)*(-     a1*p1-a2*p1-a1*p2+9.d0*a2*p2)*hvlinv * c112
!!$       end do
!!$    case(48)
!!$       do ne = 1, nemax
!!$          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne)
!!$          p1 = vv(ne-1) ; p2 = vv(ne) ; hvlinv = 1.d0 / hv(ne)
!!$          x(ne,1) =-(3.d0*a1*p1+a2*p1+a1*p2+     a2*p2)*(b2-b1)*hvlinv * c112
!!$          x(ne,2) = (3.d0*a1*p1+a2*p1+a1*p2+     a2*p2)*(b2-b1)*hvlinv * c112
!!$          x(ne,3) =-(     a1*p1+a2*p1+a1*p2+3.d0*a2*p2)*(b2-b1)*hvlinv * c112
!!$          x(ne,4) = (     a1*p1+a2*p1+a1*p2+3.d0*a2*p2)*(b2-b1)*hvlinv * c112
!!$       end do
!!$    case(49)
!!$       do ne = 1, nemax
!!$          a1 = a(ne-1) ; a2 = a(ne) ; b1 = b(ne-1) ; b2 = b(ne) ; c1 = c(ne-1) ; c2 = c(ne)
!!$          p1 = vv(ne-1) ; p2 = vv(ne) ; hvlinv = 1.d0 / hv(ne)
!!$          x(ne,1) = (3.d0*p1*(4.d0*a1+a2)+p2*(3.d0*a1+2.d0*a2))*(b2-b1)*(c2-c1)*hvlinv*c160
!!$          x(ne,2) = (p1*(3.d0*a1+2.d0*a2)+p2*(2.d0*a1+3.d0*a2))*(b2-b1)*(c2-c1)*hvlinv*c160
!!$          x(ne,3) = x(ne,2)
!!$          x(ne,4) = (p1*(2.d0*a1+3.d0*a2)+3.d0*p2*(a1+4.d0*a2))*(b2-b1)*(c2-c1)*hvlinv*c160
!!$       end do
    case default
       write(6,*)  'XX falut ID in fem_int, id= ',id
       stop
    end select

  end function fem_int

  function fem_int_point(id,ne,a,b,c) result(x)
!-------------------------------------------------------
!
!   Calculate "\int_{vv_i}^{vv_{i+1}} function(vv) dvv"
!      dvv : mesh interval
!      a   : coefficient vector
!      b   : coefficient vector
!      u   : variable vector
!      w   : weighting vector
!
!   function(vv) is classified as follows:
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
!      id = 20 : a * b * u * w
!      id = 21 : a * b'* u * w
!      id = 22 : a * b * u'* w
!      id = 23 : a * b * u * w'
!      id = 24 :(a * b * u)'* w
!      id = 25 : a'* b'* u * w
!      id = 26 : a * b'* u * w'
!      id = 27 : a * b * u'* w'
!
!      id = -1 : a * w
!      id = -2 : a * b * w
!      id = -8 : a * w'
!      id = -9 : a * b * w'
!
!      < obsolete >
!
!      id = 15 : vv * a * u * w
!      id = 16 : vv * a'* u * w
!      id = 17 : vv * a * u'* w
!      id = 18 : vv * a * u'* w'
!      id = 19 : vv * a * u * w'
!      id = 20 : vv * a * b'* u * w'
!
!      id = 36 : vv * a'* b * u * w
!      id = 37 : vv * a * b'* u * w
!      id = 38 : vv * a * b * u'* w
!      id = 39 :(vv * a * b * u)'* w
!      id = 40 : vv * a'* b * u * w'
!      id = 41 : vv * a * b * u'* w'
!      id = 42 :(vv * a * b * u)'* w'
!      id = 43 : vv * a * b * (u / b)'* w'
!
!      id = 44 : vv * a * b * u * w
!      id = 45 :(vv * a * b'* u)'* w
!      id = 46 :(vv * a * b'* u)'* w'
!      id = 47 : vv * a'* b'* u * w
!      id = 48 : vv * a * b'* u'* w
!
!      id = 49 : vv * a * b'* c'* u * w
!
!   where ' means the derivative of vv
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
    real(8), intent(in), dimension(:), optional  :: a, b, c
    integer(4) :: nel, node1, node2, iflag
    real(8) :: x(1:4), a1, a2, p1, p2, b1, b2, c1, c2, hvl, hvlinv

    cp = SUPGstab * csq15

    iflag = 0
    if(ne == 0) then
       nel = ne + 1
       iflag = 1
    else
       nel = ne
    end if

    node1 = nel  ; node2 = nel+1
    if(present(a)) then
       a1 = a(node1) ; a2 = a(node2)
       if(present(b)) then
          b1 = b(node1) ; b2 = b(node2)
       end if
       if(present(c)) then
          c1 = c(node1) ; c2 = c(node2)
       end if
    end if
    if(iflag == 1) a1 = 0.d0
    p1 = vv(nel-1) ; p2 = vv(nel) ; hvl = hv(nel) ; hvlinv = 1.d0 / hv(nel)

    select case(id)
    case(-1)
       x(1) = hvl * c13 * a1
       x(2) = hvl * c16 * a2
       x(3) = hvl * c16 * a1
       x(4) = hvl * c13 * a2
    case(-2)
       x(1) = ( 3.d0 * a1 +        a2) * hvl * c112 * b1
       x(2) = (        a1 +        a2) * hvl * c112 * b2
       x(3) = (        a1 +        a2) * hvl * c112 * b1
       x(4) = (        a1 + 3.d0 * a2) * hvl * c112 * b2
    case(-8)
       x(1) =-0.5d0 * a1
       x(2) =-0.5d0 * a2
       x(3) = 0.5d0 * a1
       x(4) = 0.5d0 * a2
    case(-9)
       x(1) = (-2.d0 * a1 -        a2) * c16 * b1
       x(2) = (-       a1 - 2.d0 * a2) * c16 * b2
       x(3) = ( 2.d0 * a1 +        a2) * c16 * b1
       x(4) = (        a1 + 2.d0 * a2) * c16 * b2
!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    case(0)
       x(1) = hvl * cp
       x(2) = x(1)
       x(3) = x(1)
       x(4) = x(1)
    case(1)
       x(1) = hvl * c13
       x(2) = hvl * c16
       x(3) = hvl * c16
       x(4) = hvl * c13
    case(2)
       x(1) = ( 3.d0 * a1 +        a2) * hvl * c112
       x(2) = (        a1 +        a2) * hvl * c112
       x(3) = (        a1 +        a2) * hvl * c112
       x(4) = (        a1 + 3.d0 * a2) * hvl * c112
    case(3)
       x(1) = (-4.d0 * a1 +        a2) * c16
       x(2) = (        a1 + 2.d0 * a2) * c16
       x(3) = (-2.d0 * a1 -        a2) * c16
       x(4) = (-       a1 + 4.d0 * a2) * c16
    case(4)
       x(1) = -0.5d0
       x(2) =  0.5d0
       x(3) = -0.5d0
       x(4) =  0.5d0
    case(5)
       x(1) = (-2.d0 * a1 -        a2) * c16
       x(2) = ( 2.d0 * a1 +        a2) * c16
       x(3) = (-       a1 - 2.d0 * a2) * c16
       x(4) = (        a1 + 2.d0 * a2) * c16
    case(6)
       x(1) = (- a1 + a2) * c13
       x(2) = (- a1 + a2) * c16
       x(3) = (- a1 + a2) * c16
       x(4) = (- a1 + a2) * c13
    case(7)
       x(1) = 0.5d0 * ( a1 - a2) * hvlinv
       x(2) = 0.5d0 * ( a1 - a2) * hvlinv
       x(3) = 0.5d0 * (-a1 + a2) * hvlinv
       x(4) = 0.5d0 * (-a1 + a2) * hvlinv
    case(8)
       x(1) =-0.5d0
       x(2) =-0.5d0
       x(3) = 0.5d0
       x(4) = 0.5d0
    case(9)
       x(1) = (-2.d0 * a1 -        a2) * c16
       x(2) = (-       a1 - 2.d0 * a2) * c16
       x(3) = ( 2.d0 * a1 +        a2) * c16
       x(4) = (        a1 + 2.d0 * a2) * c16
    case(10)
       x(1) = a1 * hvlinv
       x(2) =-a2 * hvlinv
       x(3) =-a1 * hvlinv
       x(4) = a2 * hvlinv
    case(11)
       x(1) = hvlinv
       x(2) =-hvlinv
       x(3) =-hvlinv
       x(4) = hvlinv
    case(12)
       x(1) = 0.5d0 * ( a1 + a2) * hvlinv
       x(2) = 0.5d0 * (-a1 - a2) * hvlinv
       x(3) = 0.5d0 * (-a1 - a2) * hvlinv
       x(4) = 0.5d0 * ( a1 + a2) * hvlinv
    case(13)
       x(1) = 0.5d0 * ( a1 - a2) * hvlinv
       x(2) = 0.5d0 * ( a1 - a2) * hvlinv
       x(3) = 0.5d0 * (-a1 + a2) * hvlinv
       x(4) = 0.5d0 * (-a1 + a2) * hvlinv
    case(14)
       x(1) = (-a1 + a2) * hvlinv * hvlinv
       x(2) = ( a1 - a2) * hvlinv * hvlinv
       x(3) = ( a1 - a2) * hvlinv * hvlinv
       x(4) = (-a1 + a2) * hvlinv * hvlinv
    case(20)
       x(1) = (12.d0*a1*b1 + 3.d0*a2*b1 + 3.d0*a1*b2 + 2.d0*a2*b2) * hvl * c160
       x(2) = ( 3.d0*a1*b1 + 2.d0*a2*b1 + 2.d0*a1*b2 + 3.d0*a2*b2) * hvl * c160
       x(3) = ( 3.d0*a1*b1 + 2.d0*a2*b1 + 2.d0*a1*b2 + 3.d0*a2*b2) * hvl * c160
       x(4) = ( 2.d0*a1*b1 + 3.d0*a2*b1 + 3.d0*a1*b2 +12.d0*a2*b2) * hvl * c160
    case(21)
       x(1) = (-3.d0*a1*b1 + 3.d0*a1*b2 -      a2*b1 +      a2*b2) * c112
       x(2) = (-     a1*b1 +      a1*b2 -      a2*b1 +      a2*b2) * c112
       x(3) = (-     a1*b1 +      a1*b2 -      a2*b1 +      a2*b2) * c112
       x(4) = (-     a1*b1 +      a1*b2 - 3.d0*a2*b1 + 3.d0*a2*b2) * c112
    case(22)
       x(1) = (-3.d0*a1*b1 -      a1*b2 -      a2*b1 -      a2*b2) * c112
       x(2) = ( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) * c112
       x(3) = (-     a1*b1 -      a1*b2 -      a2*b1 - 3.d0*a2*b2) * c112
       x(4) = (      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) * c112
    case(23)
       x(1) = (-3.d0*a1*b1 -      a1*b2 -      a2*b1 -      a2*b2) * c112
       x(2) = (-     a1*b1 -      a1*b2 -      a2*b1 - 3.d0*a2*b2) * c112
       x(3) = ( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) * c112
       x(4) = (      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) * c112
    case(24)
       x(1) = (-3.d0*a1*b1 + 3.d0*a1*b2 -      a2*b1 +      a2*b2) * c112 &
            &+(-3.d0*a1*b1 -      a1*b2 + 3.d0*a2*b1 +      a2*b2) * c112 &
            &+(-3.d0*a1*b1 -      a1*b2 -      a2*b1 -      a2*b2) * c112
       x(2) = (-     a1*b1 +      a1*b2 -      a2*b1 +      a2*b2) * c112 &
            &+(-     a1*b1 -      a1*b2 +      a2*b1 +      a2*b2) * c112 &
            &+( 3.d0*a1*b1 +      a1*b2 +      a2*b1 +      a2*b2) * c112
       x(3) = (-     a1*b1 +      a1*b2 -      a2*b1 +      a2*b2) * c112 &
            &+(-     a1*b1 -      a1*b2 +      a2*b1 +      a2*b2) * c112 &
            &+(-     a1*b1 -      a1*b2 -      a2*b1 - 3.d0*a2*b2) * c112
       x(4) = (-     a1*b1 +      a1*b2 - 3.d0*a2*b1 + 3.d0*a2*b2) * c112 &
            &+(-     a1*b1 - 3.d0*a1*b2 +      a2*b1 + 3.d0*a2*b2) * c112 &
            &+(      a1*b1 +      a1*b2 +      a2*b1 + 3.d0*a2*b2) * c112
    case(25)
       x(1) = (-3.d0*a1*b1 -      a1*b2 + 3.d0*a2*b1 +      a2*b2) * c112
       x(2) = (-     a1*b1 -      a1*b2 +      a2*b1 +      a2*b2) * c112
       x(3) = (-     a1*b1 -      a1*b2 +      a2*b1 +      a2*b2) * c112
       x(4) = (-     a1*b1 - 3.d0*a1*b2 +      a2*b1 + 3.d0*a2*b2) * c112
    case(26)
       x(1) = ( 2.d0*a1*b1 - 2.d0*a1*b2 +      a2*b1 -      a2*b2) * hvlinv * c16
       x(2) = (      a1*b1 -      a1*b2 + 2.d0*a2*b1 - 2.d0*a2*b2) * hvlinv * c16
       x(3) = (-2.d0*a1*b1 + 2.d0*a1*b2 -      a2*b1 +      a2*b2) * hvlinv * c16
       x(4) = (-     a1*b1 +      a1*b2 - 2.d0*a2*b1 + 2.d0*a2*b2) * hvlinv * c16
    case(27)
       x(1) = (2.d0*a1*b1 + a2*b1 + a1*b2 + 2.d0*a2*b2) * hvlinv * c16
       x(2) =-(2.d0*a1*b1 + a2*b1 + a1*b2 + 2.d0*a2*b2) * hvlinv * c16
       x(3) =-(2.d0*a1*b1 + a2*b1 + a1*b2 + 2.d0*a2*b2) * hvlinv * c16
       x(4) = (2.d0*a1*b1 + a2*b1 + a1*b2 + 2.d0*a2*b2) * hvlinv * c16
!!$    case(15)
!!$       x(1) = (12.d0*p1*a1 + 3.d0*p2*a1 + 3.d0*p1*a2 + 2.d0*p2*a2) * hvl * c160
!!$       x(2) = ( 3.d0*p1*a1 + 2.d0*p2*a1 + 2.d0*p1*a2 + 3.d0*p2*a2) * hvl * c160
!!$       x(3) = ( 3.d0*p1*a1 + 2.d0*p2*a1 + 2.d0*p1*a2 + 3.d0*p2*a2) * hvl * c160
!!$       x(4) = ( 2.d0*p1*a1 + 3.d0*p2*a1 + 3.d0*p1*a2 +12.d0*p2*a2) * hvl * c160
!!$    case(16)
!!$       x(1) = (3.d0 * p1 +        p2) * (-a1 + a2) * c112
!!$       x(2) = (       p1 +        p2) * (-a1 + a2) * c112
!!$       x(3) = (       p1 +        p2) * (-a1 + a2) * c112
!!$       x(4) = (       p1 + 3.d0 * p2) * (-a1 + a2) * c112
!!$    case(17)
!!$       x(1) = (-3.d0*p1*a1 - p2*a1 - p1*a2 -      p2*a2) * c112
!!$       x(2) = ( 3.d0*p1*a1 + p2*a1 + p1*a2 +      p2*a2) * c112
!!$       x(3) = (-     p1*a1 - p2*a1 - p1*a2 - 3.d0*p2*a2) * c112
!!$       x(4) = (      p1*a1 + p2*a1 + p1*a2 + 3.d0*p2*a2) * c112
!!$    case(18)
!!$       x(1) = ( 2.d0*p1*a1 + p2*a1 + p1*a2 + 2.d0*p2*a2) * hvlinv * c16
!!$       x(2) = (-2.d0*p1*a1 - p2*a1 - p1*a2 - 2.d0*p2*a2) * hvlinv * c16
!!$       x(3) = (-2.d0*p1*a1 - p2*a1 - p1*a2 - 2.d0*p2*a2) * hvlinv * c16
!!$       x(4) = ( 2.d0*p1*a1 + p2*a1 + p1*a2 + 2.d0*p2*a2) * hvlinv * c16
!!$    case(19)
!!$       x(1) = (-3.d0*p1*a1 - p2*a1 - p1*a2 -      p2*a2) * c112
!!$       x(2) = (-     p1*a1 - p2*a1 - p1*a2 - 3.d0*p2*a2) * c112
!!$       x(3) = ( 3.d0*p1*a1 + p2*a1 + p1*a2 +      p2*a2) * c112
!!$       x(4) = (      p1*a1 + p2*a1 + p1*a2 + 3.d0*p2*a2) * c112
!!$    case(20)
!!$       x(1) = (3.d0*p1*a1 + p2*a1 + p1*a2 +      p2*a2) * (b1 - b2) * hvlinv * c112
!!$       x(2) = (     p1*a1 + p2*a1 + p1*a2 + 3.d0*p2*a2) * (b1 - b2) * hvlinv * c112
!!$       x(3) =-(3.d0*p1*a1 + p2*a1 + p1*a2 +      p2*a2) * (b1 - b2) * hvlinv * c112
!!$       x(4) =-(     p1*a1 + p2*a1 + p1*a2 + 3.d0*p2*a2) * (b1 - b2) * hvlinv * c112
!!$    case(36)
!!$       x(1) = (12.d0*p1*b1 + 3.d0*p2*b1 + 3.d0*p1*b2 + 2.d0*p2*b2) * (-a1+a2) * c160
!!$       x(2) = ( 3.d0*p1*b1 + 2.d0*p2*b1 + 2.d0*p1*b2 + 3.d0*p2*b2) * (-a1+a2) * c160
!!$       x(3) = ( 3.d0*p1*b1 + 2.d0*p2*b1 + 2.d0*p1*b2 + 3.d0*p2*b2) * (-a1+a2) * c160
!!$       x(4) = ( 2.d0*p1*b1 + 3.d0*p2*b1 + 3.d0*p1*b2 +12.d0*p2*b2) * (-a1+a2) * c160
!!$    case(37)
!!$       x(1) = (12.d0*p1*a1 + 3.d0*p2*a1 + 3.d0*p1*a2 + 2.d0*p2*a2) * (-b1+b2) * c160
!!$       x(2) = ( 3.d0*p1*a1 + 2.d0*p2*a1 + 2.d0*p1*a2 + 3.d0*p2*a2) * (-b1+b2) * c160
!!$       x(3) = ( 3.d0*p1*a1 + 2.d0*p2*a1 + 2.d0*p1*a2 + 3.d0*p2*a2) * (-b1+b2) * c160
!!$       x(4) = ( 2.d0*p1*a1 + 3.d0*p2*a1 + 3.d0*p1*a2 +12.d0*p2*a2) * (-b1+b2) * c160
!!$    case(38)
!!$       x(1) =-( 12.d0*p1*a1*b1 + 3.d0*p2*a1*b1 + 3.d0*p1*a2*b1 + 2.d0*p2*a2*b1 &
!!$            &  + 3.d0*p1*a1*b2 + 2.d0*p2*a1*b2 + 2.d0*p1*a2*b2 + 3.d0*p2*a2*b2) * c160
!!$       x(2) = ( 12.d0*p1*a1*b1 + 3.d0*p2*a1*b1 + 3.d0*p1*a2*b1 + 2.d0*p2*a2*b1 &
!!$            &  + 3.d0*p1*a1*b2 + 2.d0*p2*a1*b2 + 2.d0*p1*a2*b2 + 3.d0*p2*a2*b2) * c160
!!$       x(3) =-(  3.d0*p1*a1*b1 + 2.d0*p2*a1*b1 + 2.d0*p1*a2*b1 + 3.d0*p2*a2*b1 &
!!$            &  + 2.d0*p1*a1*b2 + 3.d0*p2*a1*b2 + 3.d0*p1*a2*b2 +12.d0*p2*a2*b2) * c160
!!$       x(4) = (  3.d0*p1*a1*b1 + 2.d0*p2*a1*b1 + 2.d0*p1*a2*b1 + 3.d0*p2*a2*b1 &
!!$            &  + 2.d0*p1*a1*b2 + 3.d0*p2*a1*b2 + 3.d0*p1*a2*b2 +12.d0*p2*a2*b2) * c160
!!$    case(39)
!!$       x(1) = (-48.d0*a1*b1*p1+3.d0*a2*b1*p1+3.d0*a1*b2*p1+ 2.d0*a2*b2*p1 &
!!$            &  + 3.d0*a1*b1*p2+2.d0*a2*b1*p2+2.d0*a1*b2*p2+ 3.d0*a2*b2*p2) * c160
!!$       x(2) = (  3.d0*a1*b1*p1+2.d0*a2*b1*p1+2.d0*a1*b2*p1+ 3.d0*a2*b2*p1 &
!!$            &  + 2.d0*a1*b1*p2+3.d0*a2*b1*p2+3.d0*a1*b2*p2+12.d0*a2*b2*p2) * c160
!!$       x(3) = (-12.d0*a1*b1*p1-3.d0*a2*b1*p1-3.d0*a1*b2*p1- 2.d0*a2*b2*p1 &
!!$            &  - 3.d0*a1*b1*p2-2.d0*a2*b1*p2-2.d0*a1*b2*p2- 3.d0*a2*b2*p2) * c160
!!$       x(4) = (- 3.d0*a1*b1*p1-2.d0*a2*b1*p1-2.d0*a1*b2*p1- 3.d0*a2*b2*p1 &
!!$            &  - 2.d0*a1*b1*p2-3.d0*a2*b1*p2-3.d0*a1*b2*p2+48.d0*a2*b2*p2) * c160
!!$    case(40)
!!$       x(1) = (3.d0*p1*b1+p2*b1+p1*b2+     p2*b2) * (a1-a2) * hvlinv * c112
!!$       x(2) =-(3.d0*p1*b1+p2*b1+p1*b2+     p2*b2) * (a1-a2) * hvlinv * c112
!!$       x(3) = (     p1*b1+p2*b1+p1*b2+3.d0*p2*b2) * (a1-a2) * hvlinv * c112
!!$       x(4) =-(     p1*b1+p2*b1+p1*b2+3.d0*p2*b2) * (a1-a2) * hvlinv * c112
!!$    case(41)
!!$       x(1) = (b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
!!$          & * hvlinv * c112
!!$       x(2) =-(b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
!!$          & * hvlinv * c112
!!$       x(3) =-(b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
!!$          & * hvlinv * c112
!!$       x(4) = (b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
!!$          & * hvlinv * c112
!!$    case(42)
!!$       x(1) = a1*b1*p1 * hvlinv
!!$       x(2) =-a2*b2*p2 * hvlinv
!!$       x(3) =-a1*b1*p1 * hvlinv
!!$       x(4) = a2*b2*p2 * hvlinv
!!$    case(43)
!!$       x(1) = (b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
!!$          & / (12.d0 * b1 * hvl)
!!$       x(2) =-(b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
!!$          & / (12.d0 * b2 * hvl)
!!$       x(3) =-(b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
!!$          & / (12.d0 * b1 * hvl)
!!$       x(4) = (b1*(3.d0*p1*a1+p2*a1+p1*a2+p2*a2)+b2*(p1*a1+p2*a1+p1*a2+3.d0*p2*a2)) &
!!$          & / (12.d0 * b2 * hvl)
!!$    case(44)
!!$       x(1) = ( 10.d0*a1*b1*p1+2.d0*a2*b1*p1+2.d0*a1*b2*p1+      a2*b2*p1 &
!!$            &  + 2.d0*a1*b1*p2+     a2*b1*p2+     a1*b2*p2+      a2*b2*p2) * hvl * c160
!!$       x(2) = (  2.d0*a1*b1*p1+     a2*b1*p1+     a1*b2*p1+      a2*b2*p1 &
!!$            &  +      a1*b1*p2+     a2*b1*p2+     a1*b2*p2+ 2.d0*a2*b2*p2) * hvl * c160
!!$       x(3) = (  2.d0*a1*b1*p1+     a2*b1*p1+     a1*b2*p1+      a2*b2*p1 &
!!$            &  +      a1*b1*p2+     a2*b1*p2+     a1*b2*p2+ 2.d0*a2*b2*p2) * hvl * c160
!!$       x(4) = (       a1*b1*p1+     a2*b1*p1+     a1*b2*p1+ 2.d0*a2*b2*p1 &
!!$            &  +      a1*b1*p2+2.d0*a2*b1*p2+2.d0*a1*b2*p2+10.d0*a2*b2*p2) * hvl * c160
!!$    case(45)
!!$       x(1) = (-9.d0*a1*p1+a2*p1+a1*p2+     a2*p2)*(b2-b1)/(12.d0*hvl)
!!$       x(2) = (      a1*p1+a2*p1+a1*p2+3.d0*a2*p2)*(b2-b1)/(12.d0*hvl)
!!$       x(3) = (-3.d0*a1*p1-a2*p1-a1*p2-     a2*p2)*(b2-b1)/(12.d0*hvl)
!!$       x(4) = (-     a1*p1-a2*p1-a1*p2+9.d0*a2*p2)*(b2-b1)/(12.d0*hvl)
!!$    case(46)
!!$       x(1) =-a1*(b1-b2)*p1 * hvlinv * hvlinv
!!$       x(2) = a2*(b1-b2)*p2 * hvlinv * hvlinv
!!$       x(3) = a1*(b1-b2)*p1 * hvlinv * hvlinv
!!$       x(4) =-a2*(b1-b2)*p2 * hvlinv * hvlinv
!!$    case(47)
!!$       x(1) = (3.d0*p1+     p2)*(a2-a1)*(b2-b1) * hvlinv * c112
!!$       x(2) = (     p1+     p2)*(a2-a1)*(b2-b1) * hvlinv * c112
!!$       x(3) = (     p1+     p2)*(a2-a1)*(b2-b1) * hvlinv * c112
!!$       x(4) = (     p1+3.d0*p2)*(a2-a1)*(b2-b1) * hvlinv * c112
!!$    case(48)
!!$       x(1) =-(3.d0*a1*p1+a2*p1+a1*p2+     a2*p2)*(b2-b1) * hvlinv * c112
!!$       x(2) = (3.d0*a1*p1+a2*p1+a1*p2+     a2*p2)*(b2-b1) * hvlinv * c112
!!$       x(3) =-(     a1*p1+a2*p1+a1*p2+3.d0*a2*p2)*(b2-b1) * hvlinv * c112
!!$       x(4) = (     a1*p1+a2*p1+a1*p2+3.d0*a2*p2)*(b2-b1) * hvlinv * c112
!!$    case(49)
!!$       x(1) = (3.d0*p1*(4.d0*a1+a2)+p2*(3.d0*a1+2.d0*a2))*(b2-b1)*(c2-c1) * hvlinv * c160
!!$       x(2) = (p1*(3.d0*a1+2.d0*a2)+p2*(2.d0*a1+3.d0*a2))*(b2-b1)*(c2-c1) * hvlinv * c160
!!$       x(3) = x(2)
!!$       x(4) = (p1*(2.d0*a1+3.d0*a2)+3.d0*p2*(a1+4.d0*a2))*(b2-b1)*(c2-c1) * hvlinv * c160
    case default
       write(6,*)  'XX falut ID in fem_int_point, id= ',id
       stop
    end select

  end function fem_int_point

!!$  subroutine inv_int(ne_in,val,x1,x2,x)
!!$
!!$    integer(4), intent(in)  :: ne_in
!!$    real(8), intent(in)  :: val, x1, x2
!!$    real(8), intent(out) :: x
!!$
!!$    integer(4) :: ne
!!$    real(8) :: f(1:4), a1, a2, suml, lhs
!!$
!!$    ! left grid
!!$
!!$    ne = ne_in
!!$    a1 = x1
!!$
!!$    f(1) = hv(ne) / 3.d0 * a1
!!$    f(3) = hv(ne) / 6.d0 * a1
!!$
!!$    suml = f(1) + f(3)
!!$
!!$    f(2) = hv(ne) / 6.d0
!!$    f(4) = hv(ne) / 3.d0
!!$
!!$    lhs = f(2) + f(4)
!!$
!!$    ! right grid
!!$
!!$    ne = ne_in + 1
!!$    a2 = x2
!!$
!!$    f(2) = hv(ne) / 6.d0 * a2
!!$    f(4) = hv(ne) / 3.d0 * a2
!!$
!!$    suml = 0.5d0 * (suml + f(2) + f(4))
!!$
!!$    f(1) = hv(ne) / 3.d0
!!$    f(3) = hv(ne) / 6.d0
!!$
!!$    lhs = 0.5d0 * (lhs + f(1) + f(3))
!!$
!!$    x = (- val - suml) / lhs
!!$
!!$  end subroutine inv_int

!***************************************************************
!
!   Integral Method
!
!***************************************************************

! === ATTENTION !!! ======================================================!
!   These functions may be used only if size(X) is equivalent to NEMAX,   !
!   denoting the maximum index of the integral domain of fem_int funtion. !

  real(8) function intg_vol(X)

    ! Calculate \int X dV

    real(8), dimension(:), intent(in) :: X

    intg_vol = sum(fem_int(-1,X))

  end function intg_vol

  real(8) function intg_area(X)

    ! Calculate \int X dS = \int X <1/R>/(2 Pi) dV
    !   Note: ait = <1/R> = 2 Pi dS/dV

    use tx_commons, only : Pi, ait
    real(8), dimension(:), intent(in) :: X

    intg_area = sum(fem_int(-2,ait,X)) / ( 2.d0 * Pi )

  end function intg_area
!                                                                         !
! ========================================================================!

  real(8) function intg_vol_p(X,NR)

    ! Integrate X at a certain ONE point (NOT the domain)

    ! Calculate \int X dV

    integer(4), intent(in) :: NR
    real(8), dimension(:), intent(in) :: X
    integer(4) :: NE

    if(NR == 0) THEN
       intg_vol_p = 0.d0
    else
       NE = NR
       intg_vol_p = sum(fem_int_point(-1,NE,X))
    end if

  end function intg_vol_p


  real(8) function intg_area_p(X,NR)

    ! Integrate X at a certain ONE point (NOT the domain)

    ! Calculate \int X dV

    use tx_commons, only : Pi, ait
    integer(4), intent(in) :: NR
    real(8), dimension(:), intent(in) :: X
    integer(4) :: NE

    if(NR == 0) THEN
       intg_area_p = 0.d0
    else
       NE = NR
       intg_area_p = sum(fem_int_point(-2,NE,ait,X)) / ( 2.d0 * Pi )
    end if

  end function intg_area_p

  subroutine sub_intg_vol(X,NRLMAX,VAL,NR_START)

    ! Integrate X in the arbitrary size domain from NR_START to NRLMAX

    ! Calculate \int_{V(NR_START)}^V(NRLMAX) X dV

    integer(4), intent(in) :: NRLMAX
    integer(4), intent(in), optional :: NR_START
    real(8), dimension(:), intent(in) :: X
    real(8), intent(out) :: VAL
    integer(4) :: NE, NEMAX, NE_START
    real(8) :: SUML

    NEMAX = NRLMAX

    if(present(NR_START)) then
       NE_START = NR_START
    else
       ! default
       NE_START = 1
    end if

    SUML = 0.D0
    do NE = NE_START, NEMAX
       SUML = SUML + sum(fem_int_point(-1,NE,X))
    end do
    VAL = SUML

  end subroutine sub_intg_vol

end module tx_core_module

!*****************************************************************

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
        WRITE(KVALUE,*) REAL(D)
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

  D = real(GR,8)
  CALL APDTOS(STR, NSTR, D, FORM)

  RETURN
END SUBROUTINE APRTOS

!***************************************************************
!
!   Convert Strings to Upper Case
!
!***************************************************************

SUBROUTINE TOUPPERX(KTEXT)

  implicit none
  character(len=*), INTENT(INOUT) ::  KTEXT

  integer(4) :: i

  do i = 1, len_trim(KTEXT)
     if ( KTEXT(i:i) >= 'a' .and. KTEXT(i:i) <= 'z' ) then
        KTEXT(i:i) = char( ichar(KTEXT(i:i))-32 )
     end if
  end do

  RETURN
END SUBROUTINE TOUPPERX

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

pure REAL(8) FUNCTION DERIVF(NR,R,F,NRMAX)

  implicit none
  integer(4), intent(in) :: NR, NRMAX
  real(8), dimension(0:NRMAX), intent(in)  :: R
  real(8), dimension(0:NRMAX), intent(in)  :: F
  real(8) :: DR1, DR2

  IF(NR == 0) THEN
     DR1 = R(NR+1) - R(NR)
     DR2 = R(NR+2) - R(NR)
     DERIVF = (DR2**2 * F(NR+1) - DR1**2 * F(NR+2)) / (DR1 * DR2 * (DR2 - DR1)) &
          &- (DR2 + DR1) * F(NR) / (DR1 * DR2)
  ELSE IF(NR == NRMAX) THEN
     DR1 = R(NR-1) - R(NR)
     DR2 = R(NR-2) - R(NR)
     DERIVF = (DR2**2 * F(NR-1) - DR1**2 * F(NR-2)) / (DR1 * DR2 * (DR2 - DR1)) &
          &- (DR2 + DR1) * F(NR) / (DR1 * DR2)
  ELSE
     DR1 = R(NR-1) - R(NR)
     DR2 = R(NR+1) - R(NR)
     DERIVF = (DR2**2 * F(NR-1) - DR1**2 * F(NR+1)) / (DR1 * DR2 * (DR2 - DR1)) &
          &- (DR2 + DR1) * F(NR) / (DR1 * DR2)
  END IF

  RETURN
END FUNCTION DERIVF

! *** Universal high-accuracy 1st-derivative routine ***
!
!   Mode=0 computes a derivative with higher accuracy.
!   Mode=1 computes a derivative faster than mode=0.
!
!   Input: x(0:nmax)     : radial coordinate
!          f(0:nmax)     : function
!          nmax          : size of array
!          mode          : 0 : Spline with 1st derivatives
!                              with O(4) accuracy at boundaries (recommended)
!                          1 : Taylor expansion derivative
!                              with O(4) accuracy
!          daxs          : optional, given dfdx(0)
!          dbnd          : optional, given dfdx(nmax)
!   Output: dfdx(0:nmax) : derivative of the function with respect to x

function dfdx(x,f,nmax,mode,daxs,dbnd)

  USE libspl1d
  USE libitp
  implicit none
  integer(4), intent(in) :: nmax, mode
  real(8), dimension(0:nmax), intent(in) :: x, f
  real(8), dimension(0:nmax) :: dfdx
  real(8), intent(in), optional :: daxs, dbnd

  integer(4) :: n, ierr
  real(8), dimension(0:nmax) :: deriv
  real(8), dimension(1:4,0:nmax) :: u

  if(present(daxs)) then
     dfdx(0)    = daxs
  else
     dfdx(0)    = deriv4(   0,x,f,nmax,0)
  end if
  if(present(dbnd)) then
     dfdx(nmax) = dbnd
  else
     dfdx(nmax) = deriv4(nmax,x,f,nmax,0)
  end if

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
!              idx      : (optional) accept negative values if idx is true
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
  USE libitp
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
        if(present(nrbound)) then
           nrbound = i - 1
           exit
        end if
        dat_out(i) = 0.d0
     end if
!!$        if(present(nrbound)) nrbound = i - 1
!!$        exit
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
  implicit none
  real(8), parameter :: pi = 3.14159265358979323846d0
  real(8), intent(in) :: x, mu, sigma
  integer(4), intent(in), optional :: norm

  f = exp(- (x - mu)**2 / (2.d0 * sigma**2))
  if(present(norm)) f = f / (sqrt(2.d0 * pi) * sigma)

end function fgaussian

!***************************************************************
!
!   Moving window averaging
!      [Numerical Recipes 3rd p.767]
!
!   (input)
!     i    : index of the position
!     f    : function of interest
!     imax : size of f
!     iend : artificial end of averaging region
!
!***************************************************************

pure real(8) function moving_average(i,f,imax,iend) result(g)
  implicit none
  integer(4), intent(in) :: i, imax
  integer(4), intent(in), optional :: iend
  real(8), dimension(0:imax), intent(in) :: f
  integer(4) :: imaxl
  real(8) :: cn

  if(present(iend)) then
     imaxl = iend + 2
     if(i >= iend) then
        g = f(i)
        return
     end if
  else
     imaxl = imax
  end if

  if(i == 0 .or. i == imaxl) then
     g = f(i)
  else if(i == 1 .or. i == imaxl - 1) then
     cn = 1.d0 / 3.d0
     g = sum(f(i-1:i+1)) * cn
  else
     cn = 1.d0 / 5.d0
     g = sum(f(i-2:i+2)) * cn
  end if

end function moving_average

!***************************************************************
!
!   Replace the value by the interpolated one
!
!   (input)
!     index  : index of the position at which "val" is evaluated
!              0 < index < NRMAX
!     xarray : radial coordinate
!     varray : values
!   (output)
!     val    : interpolated value at xarray(index+1)
!
!***************************************************************

subroutine replace_interpolate_value(val,index,xarray,varray)
  use mod_spln
  USE libspl1d
  implicit none
  integer(4), intent(in) :: index
  real(8), dimension(:), intent(in) :: xarray, varray ! 0:NRMAX
  real(8), intent(out) :: val
  integer(4) :: ierr, iarray
  real(8), dimension(size(xarray)-1) :: xl, vl, dl
  real(8), dimension(size(xarray)) :: vl_akima
  real(8), dimension(4,size(xarray)-1) :: ul
  character(len=8) :: choose = 'akima'

  iarray = size(xarray)
  xl(1:index) = xarray(1:index)
  xl(index+1:iarray-1) = xarray(index+2:iarray)
  vl(1:index) = varray(1:index)
  vl(index+1:iarray-1) = varray(index+2:iarray)

  if( choose == 'akima' ) then
     ! Akima interpolation method
     call spln( vl_akima, xarray, iarray, vl, xl, iarray-1, 0 )
     val = vl_akima(index+1)
  else
     ! C-spline interpolation method
     call spl1d(xl,vl,dl,ul,iarray-1,0,ierr)
     call spl1df(xarray(index+1),val,xl,ul,iarray-1,ierr)
  end if

end subroutine replace_interpolate_value
