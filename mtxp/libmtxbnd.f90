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
      PUBLIC mtx_add_matrix
      PUBLIC mtx_set_source
      PUBLIC mtx_add_source
      PUBLIC mtx_set_vector
      PUBLIC mtx_add_vector
      PUBLIC mtx_split_operation
      PUBLIC mtx_solve
      PUBLIC mtx_get_vector
      PUBLIC mtx_gather_vector
      PUBLIC mtx_cleanup
      PUBLIC mtx_finalize
      PUBLIC mtx_barrier
      PUBLIC mtx_broadcast_character
      PUBLIC mtx_broadcast_integer
      PUBLIC mtx_broadcast_real8
      PUBLIC mtx_broadcast_complex8
      PUBLIC mtx_gather_integer
      PUBLIC mtx_allgather_integer
      PUBLIC mtx_gatherv_real8
      PUBLIC mtx_allgatherv_real8
      PRIVATE

      INTEGER:: imax,jmax,joffset,ierr
      REAL(8),DIMENSION(:),POINTER:: x,b
      REAL(8),DIMENSION(:,:),POINTER:: A

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      CONTAINS

      SUBROUTINE mtx_initialize(rank_,size_)
      INTEGER,INTENT(OUT):: rank_,size_

      WRITE(6,'(A)') '# mtx_initialize: libmtxbnd'
      rank_=0
      size_=1
      RETURN
      END SUBROUTINE mtx_initialize

      SUBROUTINE mtx_setup(imax_,istart_,iend_,jwidth_)

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
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      REAL(8),INTENT(IN):: v    ! value to be inserted

      A(j-i+joffset,i)=v
      RETURN
      END SUBROUTINE mtx_set_matrix
      
      SUBROUTINE mtx_add_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      REAL(8),INTENT(IN):: v    ! value to be inserted

      A(j-i+joffset,i)=A(j-i+joffset,i)+v
      RETURN
      END SUBROUTINE mtx_add_matrix
      
      SUBROUTINE mtx_set_source(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      b(j)=v
      RETURN
      END SUBROUTINE mtx_set_source
      
      SUBROUTINE mtx_add_source(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      b(j)=b(j)+v
      RETURN
      END SUBROUTINE mtx_add_source
      
      SUBROUTINE mtx_set_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      x(j)=v
      RETURN
      END SUBROUTINE mtx_set_vector
      
      SUBROUTINE mtx_add_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      x(j)=x(j)+v
      RETURN
      END SUBROUTINE mtx_add_vector
      
      SUBROUTINE mtx_split_operation
      RETURN
      END SUBROUTINE mtx_split_operation
      
      SUBROUTINE mtx_solve(itype,tolerance,its)
      INTEGER,INTENT(IN):: itype     ! not used
      REAL(8),INTENT(IN):: tolerance ! not used
      INTEGER,INTENT(OUT):: its
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
      INTEGER,INTENT(IN):: j
      REAL(8),INTENT(OUT):: v
      v=x(j)
      RETURN
      END SUBROUTINE mtx_get_vector

      SUBROUTINE mtx_gather_vector(v)
      REAL(8),DIMENSION(imax),INTENT(OUT):: v
      INTEGER:: i
      DO i=1,imax
         v(i)=x(i)
      ENDDO
      RETURN
      END SUBROUTINE mtx_gather_vector

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

      SUBROUTINE mtx_broadcast_character(kdata,n)
      CHARACTER(LEN=n),INTENT(INOUT):: kdata
      INTEGER,INTENT(IN):: n      
      RETURN
      END SUBROUTINE mtx_broadcast_character

      SUBROUTINE mtx_broadcast_integer(idata,n)
      INTEGER,DIMENSION(n),INTENT(INOUT):: idata
      INTEGER,INTENT(IN):: n      
      RETURN
      END SUBROUTINE mtx_broadcast_integer

      SUBROUTINE mtx_broadcast_real8(vdata,n)
      REAL(8),DIMENSION(n),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n
            RETURN
      END SUBROUTINE mtx_broadcast_real8

      SUBROUTINE mtx_broadcast_complex8(vdata,n)
      COMPLEX(8),DIMENSION(n),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n
            RETURN
      END SUBROUTINE mtx_broadcast_complex8

      SUBROUTINE mtx_gather_integer(idata,itot,ntot)

      INTEGER,INTENT(IN):: idata
      INTEGER,INTENT(INOUT):: ntot
      INTEGER,DIMENSION(ntot),INTENT(OUT):: itot

      itot=idata
      RETURN
      END SUBROUTINE mtx_gather_integer

      SUBROUTINE mtx_allgather_integer(idata,itot,ntot)

      INTEGER,INTENT(IN):: idata
      INTEGER,INTENT(INOUT):: ntot
      INTEGER,DIMENSION(ntot),INTENT(OUT):: itot

      itot=idata
      RETURN
      END SUBROUTINE mtx_allgather_integer

      SUBROUTINE mtx_gatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa)

      INTEGER,INTENT(IN):: ndata
      INTEGER,INTENT(INOUT):: ntot
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,DIMENSION(1):: ilena,iposa
      INTEGER:: n

      DO n=1,ndata
         vtot(n)=vdata(n)
      ENDDO
      RETURN
      END SUBROUTINE mtx_gatherv_real8

      SUBROUTINE mtx_allgatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa)

      INTEGER,INTENT(IN):: ndata
      INTEGER,INTENT(INOUT):: ntot
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,DIMENSION(1):: ilena,iposa
      INTEGER:: n

      DO n=1,ndata
         vtot(n)=vdata(n)
      ENDDO
      RETURN
      END SUBROUTINE mtx_allgatherv_real8

      END MODULE libmtx
