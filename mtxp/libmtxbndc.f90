!     $Id$

! -----------------------------------------------------------------------
!
!  Description: Solves a linear system with BANDRD (direct method)
!
! -----------------------------------------------------------------------

      MODULE libmtxc

      use libmpi
      PRIVATE

      PUBLIC mtx_initialize
      PUBLIC mtx_finalize
      PUBLIC mtx_set_communicator
      PUBLIC mtx_reset_communicator

      PUBLIC mtx_setup
      PUBLIC mtx_set_matrix
      PUBLIC mtx_set_source
      PUBLIC mtx_set_vector
      PUBLIC mtx_solve
      PUBLIC mtx_get_vector
      PUBLIC mtx_gather_vector
      PUBLIC mtx_cleanup

      TYPE(mtx_mpi_type):: mtx_global
      INTEGER:: ncomm,nrank,nsize

      INTEGER:: imax,jmax,joffset,ierr
      COMPLEX(8),DIMENSION(:),POINTER:: x,b
      COMPLEX(8),DIMENSION(:,:),POINTER:: A

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      CONTAINS

      SUBROUTINE mtx_initialize(nrank_,nsize_)
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: nrank_,nsize_
      INTEGER:: ierr


      CALL MPI_Init(ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_initialize: MPI_Init: ierr=',ierr
      ncomm=0
      CALL mtx_set_communicator_global(ncomm,nrank,nsize)
      nrank_=nrank
      nsize_=nsize
      mtx_global%comm=ncomm
      mtx_global%rank=nrank
      mtx_global%size=nsize
      mtx_global%rankg=0
      mtx_global%sizeg=1
      mtx_global%rankl=nrank
      mtx_global%sizel=nsize
      return
      END SUBROUTINE mtx_initialize

!-----

      SUBROUTINE mtx_finalize
      IMPLICIT NONE
      INTEGER:: ierr

      CALL MPI_Finalize(ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_finalize: MPI_Finalize: ierr=',ierr
      END SUBROUTINE mtx_finalize

!-----

      SUBROUTINE mtx_set_communicator(mtx_mpi,nrank_,nsize_)
        IMPLICIT NONE
        TYPE(mtx_mpi_type),INTENT(IN):: mtx_mpi
        INTEGER,INTENT(OUT):: nrank_,nsize_

        CALL mtx_set_communicator_local(mtx_mpi)
        ncomm=mtx_mpi%comm
        nrank=mtx_mpi%rank
        nsize=mtx_mpi%size
        nrank_=nrank
        nsize_=nsize
        return
      END SUBROUTINE mtx_set_communicator

!-----

      SUBROUTINE mtx_reset_communicator(nrank_,nsize_)
        IMPLICIT NONE
        INTEGER,INTENT(OUT):: nrank_,nsize_

        CALL mtx_reset_communicator_local
        ncomm=mtx_global%comm
        nrank=mtx_global%rank
        nsize=mtx_global%size
        nrank_=nrank
        nsize_=nsize
        return
      END SUBROUTINE mtx_reset_communicator

      SUBROUTINE mtx_setup(imax_,istart_,iend_,jwidth,nzmax)
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
      END SUBROUTINE mtx_setup

      SUBROUTINE mtx_set_matrix(i,j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      COMPLEX(8),INTENT(IN):: v ! value to be inserted

      A(j-i+joffset,i)=v
      RETURN
      END SUBROUTINE mtx_set_matrix
      
      SUBROUTINE mtx_set_source(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j    ! vector positon j=row
      COMPLEX(8),INTENT(IN):: v ! value to be inserted

      b(j)=v
      RETURN
      END SUBROUTINE mtx_set_source
      
      SUBROUTINE mtx_set_vector(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j    ! vector positon j=row
      COMPLEX(8),INTENT(IN):: v ! value to be inserted

      x(j)=v
      RETURN
      END SUBROUTINE mtx_set_vector
      
      SUBROUTINE mtx_solve(itype,tolerance,its, &
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
      COMPLEX(8),INTENT(OUT):: v
      v=x(j)
      RETURN
      END SUBROUTINE mtx_get_vector

      SUBROUTINE mtx_gather_vector(v)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(imax),INTENT(OUT):: v
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


      END MODULE libmtxc
