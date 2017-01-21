module mod_spln
  implicit none

contains

! Quasi-Hermite interpolation, what is called "Akima" method
! [Suuchi Keisan Handbook, Ohm sha, 1990, pp.708-717]
!=======================================================================
  subroutine spln( y, x, n, y0, x0, n0, iop, ispln )
!=======================================================================
! input arguments
    integer, intent(in)  :: n, n0, iop
    integer, intent(in), optional :: ispln
    real(8), intent(in)  :: x(n), y0(n0), x0(n0)
! output arguments
    real(8), intent(out) :: y(n)
! y : interpolated value in spln, differential value in dakima
! iop : < option >
!      1 : positive value
!     -1 : negative value
!      2 : monotone increase
!     -2 : monotone decrease
!      3 : reset to  1 if all y0 are positive,
!          reset to -1 if all y0 are negative, rset to 0 else
!      4 : if all y0 are positive(negative), all y are positive(negative)
!          negative(positive) value is reset in 0.0
!   +  0 : inside is y_i(1) and outside is y_i(n_o)
!   +100 : extrapolate inside & outside (by spline)
!   +110 : extrapolate inside (by spline) and outside is y_i(n_o)
!   +120 : extrapolate outside (by spline) and inside is y_i(1)
!   +130 : extrapolate inside & outside (by linear interpolation)
!   +140 : extrapolate inside (by linear interpolation) and outside is y_i(n_o)
!   +150 : extrapolate outside (by linear interpolation) and inside is y_i(1)
!  +1000 : data of random row
! input  : iop = lop * 1000 + kop * 10 + jop
! local variables
    real(8)    dd0, dd1, dx0, dx1, dx2, dy0, dy1, dy2, ft &
         &   , w(4,n0), y00, y01, z, zy
    integer    ifl, j0, j1, kz, m, m0, m1, nm
    integer    jop, kop, lop, mop
! iop = lop * 1000 + kop * 10 + jop 
! kop <> 10, 11, 13, 14 : y(i) = y0(1)  if y(i) < y0(1) (i = 1~n)
! kop == 10, 11         : y(i) is set up by spline interpolation
!                         used x0(1,2,3) if y(i) < y0(1) (i = 1~n)
! kop == 13, 14         : y(i) is set up by linear interpolation
!                         used x0(1,2) if y(i) < y0(1) (i = 1~n)
! kop <> 10, 12, 13, 15 : y(i) = y0(n0) if y(i) > y0(n0) (i = 1~n)
! kop == 10, 12         : y(i) is set up by spline interpolation
!                         used x0(n0-2,n0-1,n0) if y(i) > y0(n0) (i = 1~n)
! kop == 13, 15         : y(i) is set up by linear interpolation
!                         used x0(n0-1,n0) if y(i) > y0(n0) (i = 1~n)
    integer    kc(-1:1)
    integer    n1
    real(8)    cm1, cm2, cm21, cm3, cm4, cm43, dft, dx0inv
    real(8)    fdif, fint
    integer    kspln
    integer :: iuse_spln = 0 ! = 0 : Akima, = 2 : original
!=======================================================================
    if( present(ispln) ) iuse_spln = ispln
    fint = 1.0_8
    goto 1000
    entry dakima( y, x, n, y0, x0, n0, iop )
    fint = 0.0_8

1000 continue
    fdif = 1.0_8 - fint

!:: set flag
    if( n == 0 ) return
    jop = mod( iop, 10 )
    if( jop >= 8 ) jop = jop - 10
    kop = ( iop - jop ) / 10
    lop = kop / 100
    kop = kop - lop * 100
    if( lop == 0 ) then
       lop = n0
    else
       lop = 2
    endif

! check n0
    kspln = iuse_spln
    if( iuse_spln == 0 ) then
! Akima method requires at least 5 points.
       if( n0 <= 4 ) stop 'Akima method requires at least 5 points.'
    elseif( iuse_spln == 1 ) then
       if( n0 == 3 .or. n0 == 4 ) kspln = 2
       if( n0 <= 2 ) stop 'Akima method requires at least 5 points.'
    elseif( iuse_spln == 2 ) then
       if( n0 <= 2 ) stop 'spline method requires at least 3 points.'
    endif

    if( kspln == 0 ) then
!:: Akima spline interpolation
       do m = 1, n0-1
! Interim derivative, share the array to save memory
          w(1,m) = ( y0(m+1) - y0(m) ) / ( x0(m+1) - x0(m) )
       enddo
       w(1,n0) = 2.0_8 * w(1,n0-1) - w(1,n0-2)

! first order coefficients: w(2,m)
       m = 1
       cm1  = 3.0_8 * w(1,m) - 2.0_8 * w(1,m+1)
       cm2  = 2.0_8 * w(1,m) -         w(1,m+1)
       cm43 = abs( w(1,m+1) - w(1,m) )
       cm21 = abs( cm2 - cm1 )
       if( cm43 * cm21 == 0.0_8 ) then
          w(2,m) = 0.5_8 * ( w(1,m) + cm2 )
       else
          w(2,m) = ( cm43 * cm2 + cm21 * w(1,m) ) / ( cm43 + cm21 )
       endif

       m = 2
       cm1  = cm2
       cm43 = abs( w(1,m+1) - w(1,m) )
       cm21 = abs( w(1,m-1) - cm1 )
       if( cm43 * cm21 == 0.0_8 ) then
          w(2,m) = 0.5_8 * ( w(1,m) + w(1,m-1) )
       else
          w(2,m) = ( cm43 * w(1,m-1) + cm21 * w(1,m) ) / ( cm43 + cm21 )
       endif

       do m = 3, n0-1
          cm43 = abs( w(1,m+1) - w(1,m  ) )
          cm21 = abs( w(1,m-1) - w(1,m-2) )
          if( cm43 == 0.0_8 .and. cm21 == 0.0_8 ) then
             w(2,m) = 0.5_8 * ( w(1,m) + w(1,m-1) )
          else
             w(2,m) = ( cm43 * w(1,m-1) + cm21 * w(1,m) ) &
                  & / ( cm43 + cm21 )
          endif
       enddo

       m = n0
       cm3  = w(1,n0)
       cm4  = 3.0_8 * w(1,m-1) - 2.0_8 * w(1,m-2)
       cm43 = abs( cm4 - cm3 )
       cm21 = abs( w(1,m-1) - w(1,m-2) )
       if( cm43 * cm21 == 0.0_8 ) then
          w(2,m) = 0.5_8 * ( cm3 + w(1,m-1) )
       else
          w(2,m) = ( cm43 * w(1,m-1) + cm21 * cm3 ) / ( cm43 + cm21 )
       endif
       w(1,m) = y0(m) ! overridden
       w(3,m) = 0.0_8
       w(4,m) = 0.0_8
       do m = 1, n0-1
          dx0inv = 1.0_8 / (x0(m+1) - x0(m) )
! second and third order coefficients: w(3,m), w(4,m)
          w(3,m) = ( 3.0_8*w(1,m) - 2.0_8*w(2,m)   -       w(2,m+1) ) &
               & * dx0inv
          w(4,m) = (       w(2,m) +       w(2,m+1) - 2.0_8*w(1,m)   ) &
               & * dx0inv**2
! constant term: w(1,m), overridden
          w(1,m) = y0(m)
       enddo
    else
!:: spline interpolation
       nm  = n0 - 2
       dx1 = x0(2) - x0(1)
       dx2 = x0(3) - x0(2)
       dy1 = ( y0(2) - y0(1) ) / dx1
       dy2 = ( y0(3) - y0(2) ) / dx2
       dd1 = (-dx1 * dy2 + ( dx2 + 2.0_8 * dx1 ) * dy1 ) / ( dx2 + dx1 )

       do m = 1, nm
          dx0    = dx1
          dy0    = dy1
          dd0    = dd1
          dx1    = x0(m+2) - x0(m+1)
          dy1    = ( y0(m+2) - y0(m+1) ) / dx1
          dd1    = ( dx0 * dy1 + dx1 * dy0 ) / ( dx1 + dx0 )
          w(1,m) = y0(m)
          w(2,m) = dd0
          w(3,m) = ( 3.0_8 * dy0 - dd1 - 2.0_8 * dd0 ) / dx0
          w(4,m) = (-2.0_8 * dy0 + dd1 + dd0 ) / ( dx0 * dx0 )
       enddo

       dd0 = dd1
       dd1 = ( ( 2.0_8 * dx1 + dx0 ) * dy1 - dx1 * dy0 ) / ( dx1 + dx0 )
       m      = n0 - 1
       w(1,m) = y0(m)
       w(2,m) = dd0
       w(3,m) = ( 3.0_8 * dy1 - dd1 - 2.0_8 * dd0 ) / dx1
       w(4,m) = (-2.0_8 * dy1 + dd1 + dd0 ) / ( dx1 * dx1 )

       w(1,n0) = y0(n0)
       w(2,n0) = dd1
       w(3,n0) = 0.0_8
       w(4,n0) = 0.0_8
    endif

!-----------------------------------------------------------------------
    j1 = 2
    do m = 1, n
       j1 = min( j1, lop )
       if( x(m) <= x0(1) ) then
          if( kop /= 10 .and. kop /= 11 .and. kop /= 13 .and. kop /= 14 ) then
             y(m) = fint * y0(1) ! interpolation
             if( fdif == 1.0_8 .and. x(m) == x0(1) ) y(m) = w(2,1) ! differential
          elseif( kop == 10 .or. kop == 11 ) then
             j0 = 1
             z  = x(m) - x0(j0)
             y(m) = fint * ( w(1,j0) + z * ( w(2,j0) &
                & + z * ( w(3,j0) + z * w(4,j0) ) ) ) & ! interpolation 
                & + fdif * ( w(2,j0) + z * ( 2.0_8 * w(3,j0) &
                & + 3.0_8 * z * w(4,j0) ) )   ! differential
          elseif( kop == 13 .or. kop == 14 ) then
             dx0inv = 1.0_8 / ( x0(1) - x0(2) )
             y(m) = fint * (  ( y0(1) - y0(2) ) * dx0inv &
                &  * ( x(m) - x0(2) ) + y0(2) ) & ! interpolation
                &  + fdif * ( -( y0(1) - y0(2) ) * dx0inv &
                &  * x0(2) )                    ! differential
          endif
       elseif( x0(n0) <= x(m) ) then
          if( kop /= 10 .and. kop /= 12 .and. kop /= 13 .and. kop /= 15 ) then
             y(m) = fint * y0(n0) ! interpolation
             if( fdif == 1.0_8 .and. x(m) == x0(n0) ) y(m) = w(2,n0) ! differential
          elseif( kop == 10 .or. kop == 12 ) then
             j0 = n0 - 1
             z  = x(m) - x0(j0)
             y(m) = fint * ( w(1,j0) + z * ( w(2,j0) &
                & + z * ( w(3,j0) + z * w(4,j0) ) ) ) & ! interpolation
                & + fdif * ( w(2,j0) + z * ( 2.0_8 * w(3,j0) &
                & + 3.0_8 * z * w(4,j0) ) )   ! differential
          elseif( kop == 13 .or. kop == 15 ) then
             n1 = n0 - 1
             dx0inv = 1.0_8 / ( x0(n0) - x0(n1) )
             y(m) = fint * (   ( y0(n0) - y0(n1) ) * dx0inv &
                & * ( x(m) - x0(n1) ) + y0(n1) ) &! interpolation 
                & + fdif *  ( -( y0(n0) - y0(n1) ) * dx0inv &
                & * x0(n1) )                     ! differential
          endif
       else
          do while( x0(j1) <= x(m) .and. j1 < n0 )
             j1 = j1 + 1
          enddo
          j0 = j1 - 1
          z  = x(m) - x0(j0)
          y(m) = fint * ( w(1,j0) + z * ( w(2,j0) &
             & + z * ( w(3,j0) + z * w(4,j0) ) ) ) & ! interpolation
             & + fdif * ( w(2,j0) + z * ( 2.0_8 * w(3,j0) &
             & + 3.0_8 * z * w(4,j0) ) )   ! differential
       endif
    enddo
!-----------------------------------------------------------------------

!:: correction processing
    if( jop == 3 .or. jop == 4 ) then
! -- check positive value or negative value for all y0(input)
       kc(-1:1) = 0
       do m = 1, n0
          if( y0(m) < 0.0_8 ) then
             kc(-1) = kc(-1) + 1
          elseif( y0(m) == 0.0_8 ) then
             kc(0) = kc(0) + 1
          else
             kc(1) = kc(1) + 1
          endif
       enddo
       mop = 0
       if( kc(-1)+kc(0) == n0 ) then
          mop = -1
       elseif( kc(1)+kc(0) == n0 ) then
          mop =  1
       endif
! reset y(1:n) to 0.0
       if( jop == 3 ) then
          jop = mop
       else
          if( mop == 1  ) where(  y(1:n) < 0.0_8 ) y(1:n) = 0.0_8
          if( mop == -1 ) where(  y(1:n) > 0.0_8 ) y(1:n) = 0.0_8
          jop = 0
       endif
    endif

    where( abs( y(1:n) ) < 1.0d-38 ) y(1:n) = 0.0_8

    if( jop == 0 ) return

    m1 = 1

    do j1 = 2, n0
       if( x(m1) > x0(j1) ) cycle
       m0 = m1
       do while( x(m1) < x0(j1) .and. m1 < n )
          m1 = m1 + 1
       enddo
       j0 = j1 - 1

       ifl = 0
       kz = abs( jop )
       if( kz == 1 ) then
          do m = m0, m1
             if( x(m) < x0(j0) .or. x(m) >= x0(j1) ) cycle
             zy = jop * y(m)
             if( zy <= 0.0_8 ) ifl = 1
          enddo
       elseif( kz == 2 ) then
          do m = m0, m1
             if( x(m) < x0(j0) .or. x(m) >= x0(j1) ) cycle
             if( m == 1 ) cycle
             zy = jop * ( y(m) - y(m-1) )
             if( zy <= 0.0_8 ) ifl = 1
          enddo
       endif
       if( ifl == 0 ) cycle

       y00 = y0(j0)
       y01 = y0(j1)
       if( jop == 1 ) then
          y00 = max( y00, 0.0_8 )
          y01 = max( y01, 0.0_8 )
       elseif( jop == -1 ) then
          y00 = min( y00, 0.0_8 )
          y01 = min( y01, 0.0_8 )
       endif

       do m = m0, m1
          if( x(m) < x0(j0) .or. x(m) >= x0(j1) ) cycle
          dx0inv = 1.0_8 / ( x0(j1) - x0(j0) )
          ft   = ( x(m) - x0(j0) ) * dx0inv
          dft  =        - x0(j0)   * dx0inv
          y(m) = fint * ( ( 1.0_8 - ft  ) * y00 + ft  * y01 ) & ! interpolation
             & + fdif * ( ( 1.0_8 - dft ) * y00 + dft * y01 ) ! differential
       enddo

    enddo

  end subroutine spln

end module mod_spln
