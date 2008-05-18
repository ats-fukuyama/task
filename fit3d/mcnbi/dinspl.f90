SUBROUTINE    dinspl(x,y,dy,n,c,d,e,icon)

  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: n
  INTEGER(4),INTENT(OUT):: icon
  REAL(8),INTENT(IN) :: x(n), y(n), dy(n)
  REAL(8),INTENT(OUT):: c(n), d(n), e(n)
  INTEGER(4):: i, i0, i1, n1, n3
  REAL(8):: dw, dw1, h, h1, yh, yh1, tm    !mcode(6)

 1000 CONTINUE
      IF( n .lt. 2 )     GO TO  9000
      n1 = n - 1
      DO  i = 1, n1
        IF( x(i) .ge. x(i+1) )     GO TO  9000
      END DO
 1100 CONTINUE
      icon = 0
      d(1) = dy(1)
      d(n) = dy(2)
      IF( n .eq. 2 )     GO TO  1300
      h1   = x(2) - x(1)
      yh1  = ( y(2) - y(1) ) / h1
      h    = x(3) - x(2)
      yh   = ( y(3) - y(2) ) / h
      d(2) = 0.5d0 / ( h1 + h )
      c(2) = 6.0d0 * ( yh - yh1 ) - h1 * d(1)
      IF(n.gt.3) GO TO 1200
      d(2)=(c(2)-h*d(3))*d(2)
      GO TO 1300
 1200 CONTINUE
      DO  i = 3, n1
        h1   = h
        yh1  = yh
        h    = x(i+1) - x(i)
        yh   = ( y(i+1) - y(i) ) / h
        tm   = h1 * d(i-1)
        d(i) = 1.0d0 / ( 2.0d0 * ( h1 + h ) - h1 * tm )
        e(i) = tm
        c(i) = 6.0d0 * ( yh - yh1 ) - tm * c(i-1)
      END DO
      d(n1) = ( c(n1) - h * d(n) ) * d(n1)
      n3    = n1 - 2
      DO  i = 1, n3
        i0 = n1 - i
        i1 = n  - i
        d(i0) = c(i0) * d(i0) - d(i1) * e(i1)
      END DO
 1300 CONTINUE
      DO  i = 1, n1
        dw   = d(i)
        dw1  = d(i+1)
        h    = x(i+1) - x(i)
        c(i) = ( y(i+1) - y(i) ) / h - h * ( dw1 + 2.0d0 * dw ) / 6.0d0
        e(i) = ( dw1 - dw ) / ( h * 6.0d0 )
        d(i) = dw * 0.5d0
      END DO
      c(n) = ( y(n) - y(n1) ) / h + h * ( d(n) + d(n1) ) / 3.0d0
      d(n) = d(n) * 0.5d0
      e(n) = 0.0d0
 8000 continue
      RETURN
 9000 CONTINUE
      icon = 30000
      GO TO  8000
END
