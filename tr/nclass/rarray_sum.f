      FUNCTION RARRAY_SUM(n,x,incx)
!***********************************************************************
!  RARRAY_SUM is the sum of the elements of the array x
!  W.A.Houlberg 12/98
!  Input:
!  n-number of elements to be summed
!  x-array to be summed
!  incx-increment in sx index
!***********************************************************************
      USE trcomm,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind)    RARRAY_SUM
!Declaration of input variables
      INTEGER        incx,                    n
      REAL(rkind)    x(*)
!Declaration of local variables
      INTEGER        i,                       ix
      RARRAY_SUM=0.D0
      ix=1
      DO i=1,n
        RARRAY_SUM=RARRAY_SUM+x(ix)
        ix=ix+incx
      ENDDO   
      RETURN
      END
