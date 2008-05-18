      SUBROUTINE ranu2(ix,a,n,icon)

  IMPLICIT NONE
  INTEGER(4),INTENT(INOUT) :: ix, n
  INTEGER(4),INTENT(OUT):: icon
  REAL(4),INTENT(OUT) :: a(n)
  INTEGER(4):: i
  REAL(8) :: x

      IF(n.lt.1.or.ix.lt.0) GO TO 8000
      x=DFLOAT(ix)/2147483648.D0
      DO i=1,n
        x=x*32771.d0+0.574890473391860723d0
        x=x-INT(x)
        a(i)=x
      END DO
      ix=x*2147483648.d0
      icon=0
      RETURN

 8000 icon=30000
      RETURN
      END
