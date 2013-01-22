!     $Id$

! -----------------------------------------------------------------------
!
!  Description: Solves a linear system with BANDRD (direct method)
!
! -----------------------------------------------------------------------

      MODULE libmtxc

      use libmpi
      PRIVATE

      PUBLIC mtxc_setup
      PUBLIC mtxc_set_matrix
      PUBLIC mtxc_set_source
      PUBLIC mtxc_set_vector
      PUBLIC mtxc_solve
      PUBLIC mtxc_get_vector
      PUBLIC mtxc_gather_vector
      PUBLIC mtxc_cleanup

      INTEGER:: imax,jmax,joffset,ierr
      COMPLEX(8),DIMENSION(:),POINTER:: x,b
      COMPLEX(8),DIMENSION(:,:),POINTER:: A

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      CONTAINS

      SUBROUTINE mtxc_setup(imax_,istart_,iend_,jwidth,nzmax)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(OUT):: istart_,iend_  ! allocated range of lines 
      INTEGER,OPTIONAL,INTENT(IN):: jwidth ! band matrix width
      INTEGER,OPTIONAL,INTENT(IN):: nzmax  ! number of nonzero components
      INTEGER:: i,j

      imax=imax_
      IF(PRESENT(jwidth)) THEN
         jmax=jwidth
      ELSE
         jmax=2*imax_-1
      ENDIF
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
      END SUBROUTINE mtxc_setup

      SUBROUTINE mtxc_set_matrix(i,j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      COMPLEX(8),INTENT(IN):: v ! value to be inserted

      A(j-i+joffset,i)=v
      RETURN
      END SUBROUTINE mtxc_set_matrix
      
      SUBROUTINE mtxc_set_source(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j    ! vector positon j=row
      COMPLEX(8),INTENT(IN):: v ! value to be inserted

      b(j)=v
      RETURN
      END SUBROUTINE mtxc_set_source
      
      SUBROUTINE mtxc_set_vector(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j    ! vector positon j=row
      COMPLEX(8),INTENT(IN):: v ! value to be inserted

      x(j)=v
      RETURN
      END SUBROUTINE mtxc_set_vector
      
      SUBROUTINE mtxc_solve(itype,tolerance,its, &
           methodKSP,methodPC,damping_factor,emax,emin,max_steps)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: itype     ! not used
      REAL(8),INTENT(IN):: tolerance ! not used
      INTEGER,INTENT(OUT):: its      ! number of iterations
      INTEGER,OPTIONAL:: methodKSP,methodPC,max_steps
      REAL(8),OPTIONAL:: damping_factor,emax,emin
      INTEGER:: i,j

      DO i=1,imax
         x(i)=b(i)
      ENDDO

!      do i=1,imax
!         write(21,'(i5/(1P5E12.4))') i,(A(j,i),j=1,jmax)
!      enddo
         
      CALL BANDCD(A,x,imax,jmax,jmax,ierr)
      IF(ierr.ne.0) then
         WRITE(6,'(A,I5)') 'XX BANDRD in mtxc_solve: ierr=',ierr
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
      END SUBROUTINE mtxc_solve

      SUBROUTINE mtxc_get_vector(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j
      COMPLEX(8),INTENT(OUT):: v
      v=x(j)
      RETURN
      END SUBROUTINE mtxc_get_vector

      SUBROUTINE mtxc_gather_vector(v)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(imax),INTENT(OUT):: v
      INTEGER:: i

      DO i=1,imax
         v(i)=x(i)
      ENDDO
      RETURN
      END SUBROUTINE mtxc_gather_vector

      SUBROUTINE mtxc_cleanup
      IMPLICIT NONE

      DEALLOCATE(x,b)
      DEALLOCATE(A)
      RETURN
      END SUBROUTINE mtxc_cleanup

      END MODULE libmtxc
