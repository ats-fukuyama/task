!     $Id$

! -----------------------------------------------------------------------
!
!  Description: Solves a linear system with BANDRD (direct method)
!
! -----------------------------------------------------------------------

      MODULE libmtx

      use libmpi

      PUBLIC mtx_initialize
      PUBLIC mtx_finalize
      PUBLIC mtx_barrier
      PUBLIC mtx_broadcast_character
      PUBLIC mtx_broadcast_integer
      PUBLIC mtx_broadcast_real8
      PUBLIC mtx_broadcast_complex8
      PUBLIC mtx_gather_integer
      PUBLIC mtx_gather_real8
      PUBLIC mtx_allgather_integer
      PUBLIC mtx_gatherv_real8
      PUBLIC mtx_allgatherv_real8

      PUBLIC mtx_setup
      PUBLIC mtx_set_matrix
      PUBLIC mtx_set_source
      PUBLIC mtx_set_vector
      PUBLIC mtx_solve
      PUBLIC mtx_get_vector
      PUBLIC mtx_gather_vector
      PUBLIC mtx_cleanup

      PRIVATE
      INTEGER:: imax,jmax,joffset,ierr
      REAL(8),DIMENSION(:),POINTER:: x,b
      REAL(8),DIMENSION(:,:),POINTER:: A

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      CONTAINS

      SUBROUTINE mtx_setup(imax_,istart_,iend_,jwidth_)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(OUT):: istart_,iend_  ! allocated range of lines 
      INTEGER,INTENT(IN):: jwidth_         ! band matrix width
      INTEGER:: i,j

      imax=imax_
      jmax=jwidth_
      ALLOCATE(A(jmax,imax))
      ALLOCATE(b(imax),x(imax))
      istart_=1
      iend_=imax
      joffset=(jmax+1)/2

      DO i=1,imax
         DO j=1,jmax
            A(J,i)=0.D0
         ENDDO
         b(i)=0.D0
         x(i)=0.d0
      ENDDO

      RETURN
      END SUBROUTINE mtx_setup

      SUBROUTINE mtx_set_matrix(i,j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      REAL(8),INTENT(IN):: v    ! value to be inserted

      A(j-i+joffset,i)=v
      RETURN
      END SUBROUTINE mtx_set_matrix
      
      SUBROUTINE mtx_set_source(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      b(j)=v
      RETURN
      END SUBROUTINE mtx_set_source
      
      SUBROUTINE mtx_set_vector(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      x(j)=v
      RETURN
      END SUBROUTINE mtx_set_vector
      
      SUBROUTINE mtx_add_vector(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      x(j)=x(j)+v
      RETURN
      END SUBROUTINE mtx_add_vector
      
      SUBROUTINE mtx_split_operation
      IMPLICIT NONE
      RETURN
      END SUBROUTINE mtx_split_operation
      
      SUBROUTINE mtx_solve(itype,tolerance,its, &
           methodKSP,methodPC,damping_factor,emax,emin,max_steps)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: itype     ! not used
      REAL(8),INTENT(IN):: tolerance ! not used
      INTEGER,INTENT(OUT):: its
      INTEGER,OPTIONAL:: methodKSP,methodPC,max_steps
      REAL(8),OPTIONAL:: damping_factor,emax,emin
      INTEGER:: i,j

      DO i=1,imax
         x(i)=b(i)
      ENDDO

!      do i=1,imax
!         write(21,'(i5/(1P5E12.4))') i,(A(j,i),j=1,jmax)
!      enddo
         
      CALL BANDRD(A,x,imax,jmax,jmax,ierr)
      IF(ierr.ne.0) then
         WRITE(6,'(A,I5)') 'XX BANDRD in mtx_solve: ierr=',ierr
         its=-1
      ELSE
         its=0
      ENDIF
      DO i=1,imax
         DO j=1,jmax
            A(J,i)=0.D0
         ENDDO
         b(i)=0.D0
      ENDDO

      RETURN
      END SUBROUTINE mtx_solve

      SUBROUTINE mtx_get_vector(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j
      REAL(8),INTENT(OUT):: v
      v=x(j)
      RETURN
      END SUBROUTINE mtx_get_vector

      SUBROUTINE mtx_gather_vector(v)
      IMPLICIT NONE
      REAL(8),DIMENSION(imax),INTENT(OUT):: v
      INTEGER:: i

      DO i=1,imax
         v(i)=x(i)
      ENDDO
      RETURN
      END SUBROUTINE mtx_gather_vector

      SUBROUTINE mtx_cleanup
      IMPLICIT NONE

      DEALLOCATE(x,b)
      DEALLOCATE(A)
      RETURN
      END SUBROUTINE mtx_cleanup

      END MODULE libmtx
