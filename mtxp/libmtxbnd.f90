!     $Id$

! -----------------------------------------------------------------------
!
!  Description: Solves a linear system with BANDRD (direct method)
!
! -----------------------------------------------------------------------

      MODULE libmtx

      IMPLICIT NONE
      PUBLIC mtx_initialize
      PUBLIC mtx_setup
      PUBLIC mtx_set_matrix
      PUBLIC mtx_set_vector
      PUBLIC mtx_solve
      PUBLIC mtx_get_vector
      PUBLIC mtx_cleanup
      PUBLIC mtx_finalize
      PUBLIC mtx_barrier
      PRIVATE

      INTEGER:: imax,jmax,joffset,ierr
      REAL(8),DIMENSION(:),POINTER:: x,b
      REAL(8),DIMENSION(:,:),POINTER:: A

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      CONTAINS

      SUBROUTINE mtx_initialize(irank,isize)
      INTEGER,INTENT(OUT):: irank,isize

      WRITE(6,*) '## libmtxbnd initialized'
      irank=0
      isize=1
      RETURN
      END SUBROUTINE mtx_initialize

      SUBROUTINE mtx_setup(i_max,j_width,i_start,i_end)

      INTEGER,INTENT(IN):: i_max           ! total matrix size
      INTEGER,INTENT(IN):: j_width         ! band matrix width
      INTEGER,INTENT(OUT):: i_start,i_end  ! allocated range of lines 

      WRITE(6,*) '## libmtxbnd 
      imax=i_max
      jmax=j_width
      ALLOCATE(A(jmax,imax))
      ALLOCATE(b(imax),x(imax))
      i_start=1
      i_end=imax
      joffset=(jmax+1)/2

      RETURN
      END SUBROUTINE mtx_setup

      SUBROUTINE mtx_set_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      REAL(8),INTENT(IN):: v    ! value to be inserted

      A(j-i+joffset,i)=v
      RETURN
      END SUBROUTINE mtx_set_matrix
      
      SUBROUTINE mtx_set_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      b(j)=v
      RETURN
      END SUBROUTINE mtx_set_vector
      
      SUBROUTINE mtx_solve(its)
      INTEGER,INTENT(OUT):: its
      INTEGER:: i

      DO i=1,imax
         x(i)=b(i)
      ENDDO
         
      CALL BANDRD(A,X,imax,jmax,jmax,ierr)
      IF(ierr.ne.0) then
         WRITE(6,'(A,I5)') 'XX BANDRD in mtx_solve: ierr=',ierr
         its=-1
      ELSE
         its=0
      ENDIF
      RETURN
      END SUBROUTINE mtx_solve

      SUBROUTINE mtx_get_vector(j,v)

      INTEGER,INTENT(IN):: j
      REAL(8),INTENT(OUT):: v

      v=X(j)
      RETURN
      END SUBROUTINE mtx_get_vector

      SUBROUTINE mtx_cleanup

      DEALLOCATE(x,b)
      DEALLOCATE(A)
      RETURN
      END SUBROUTINE mtx_cleanup

      SUBROUTINE mtx_finalize
      RETURN
      END SUBROUTINE mtx_finalize

      SUBROUTINE mtx_barrier
      RETURN
      END SUBROUTINE mtx_barrier

      END MODULE libmtx
