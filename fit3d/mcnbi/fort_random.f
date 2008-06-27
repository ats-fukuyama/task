!     $Id$

      SUBROUTINE fort_random(n,ran,init,ix,iy,ierr)

      IMPLICIT none
      INTEGER,INTENT(IN):: ix,iy,n
      INTEGER,INTENT(INOUT):: init
      INTEGER,INTENT(OUT):: ierr
      REAL(4),DIMENSION(n),INTENT(OUT):: ran
      INTEGER:: i,ns
      INTEGER,DIMENSION(:),ALLOCATABLE:: seed

      IF(init.eq.0) THEN
         CALL random_seed(size=ns)
         ALLOCATE(seed(ns))
         seed=ix+iy*(/ (i-1,i=1,ns) /)
         CALL random_seed(put=seed)
         DEALLOCATE(seed)
         init=1
      ENDIF
      CALL random_number(ran)
      ierr=0
      RETURN
      END
