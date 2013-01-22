!     $Id$

! -----------------------------------------------------------------------
!
!  Description: Solves a linear system with BANDRD/CD (direct method)
!
! -----------------------------------------------------------------------

      MODULE libmtx

      use libmpi
      use commpi
      PRIVATE

      PUBLIC mtx_initialize
      PUBLIC mtx_finalize

      PUBLIC mtx_setup
      PUBLIC mtx_set_matrix
      PUBLIC mtx_set_source
      
      PUBLIC mtx_set_vector
      PUBLIC mtx_solve
      PUBLIC mtx_get_vector
      PUBLIC mtx_gather_vector
      PUBLIC mtx_cleanup

      PUBLIC mtxc_setup
      PUBLIC mtxc_set_matrix
      PUBLIC mtxc_set_source
      PUBLIC mtxc_set_vector
      PUBLIC mtxc_solve
      PUBLIC mtxc_get_vector
      PUBLIC mtxc_gather_vector
      PUBLIC mtxc_cleanup

      TYPE(mtx_mpi_type):: mtx_global

      INTEGER:: imax,jmax,joffset,ierr
      REAL(8),DIMENSION(:),POINTER:: x,b
      REAL(8),DIMENSION(:,:),POINTER:: A
      COMPLEX(8),DIMENSION(:),POINTER:: xc,bc
      COMPLEX(8),DIMENSION(:,:),POINTER:: Ac

      CONTAINS

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of common section
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      SUBROUTINE mtx_initialize
      IMPLICIT NONE
      INTEGER:: ncomm_in,ierr

      ncomm_in=0
      CALL mtx_set_communicator_global(ncomm_in)
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

      RETURN
      END SUBROUTINE mtx_finalize

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of bandrd section
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

!-----

      SUBROUTINE mtx_set_matrix(i,j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      REAL(8),INTENT(IN):: v    ! value to be inserted

      A(j-i+joffset,i)=v
      RETURN
      END SUBROUTINE mtx_set_matrix
      
!-----

      SUBROUTINE mtx_set_source(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      b(j)=v
      RETURN
      END SUBROUTINE mtx_set_source
      
!-----

      SUBROUTINE mtx_set_vector(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      x(j)=v
      RETURN
      END SUBROUTINE mtx_set_vector
      
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

!-----

      SUBROUTINE mtx_get_vector(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j
      REAL(8),INTENT(OUT):: v
      v=x(j)
      RETURN
      END SUBROUTINE mtx_get_vector

!-----

      SUBROUTINE mtx_gather_vector(v)
      IMPLICIT NONE
      REAL(8),DIMENSION(imax),INTENT(OUT):: v
      INTEGER:: i

      DO i=1,imax
         v(i)=x(i)
      ENDDO
      RETURN
      END SUBROUTINE mtx_gather_vector

!-----

      SUBROUTINE mtx_cleanup
      IMPLICIT NONE

      DEALLOCATE(x,b)
      DEALLOCATE(A)
      RETURN
      END SUBROUTINE mtx_cleanup

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of bandcd section
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

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
      ALLOCATE(Ac(jmax,imax))
      ALLOCATE(bc(imax),xc(imax))
      istart_=1
      iend_=imax
      joffset=(jmax+1)/2

      DO i=1,imax
         DO j=1,jmax
            Ac(J,i)=0.D0
         ENDDO
         bc(i)=0.D0
         xc(i)=0.d0
      ENDDO

      RETURN
      END SUBROUTINE mtxc_setup

!-----

      SUBROUTINE mtxc_set_matrix(i,j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      COMPLEX(8),INTENT(IN):: v ! value to be inserted

      Ac(j-i+joffset,i)=v
      RETURN
      END SUBROUTINE mtxc_set_matrix
      
!-----

      SUBROUTINE mtxc_set_source(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j    ! vector positon j=row
      COMPLEX(8),INTENT(IN):: v ! value to be inserted

      bc(j)=v
      RETURN
      END SUBROUTINE mtxc_set_source
      
!-----

      SUBROUTINE mtxc_set_vector(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j    ! vector positon j=row
      COMPLEX(8),INTENT(IN):: v ! value to be inserted

      xc(j)=v
      RETURN
      END SUBROUTINE mtxc_set_vector
      
!-----

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
         xc(i)=bc(i)
      ENDDO

!      do i=1,imax
!         write(21,'(i5/(1P5E12.4))') i,(A(j,i),j=1,jmax)
!      enddo
         
      CALL BANDCD(Ac,xc,imax,jmax,jmax,ierr)
      IF(ierr.ne.0) then
         WRITE(6,'(A,I5)') 'XX BANDCD in mtxc_solve: ierr=',ierr
         its=-1
      ELSE
         its=0
      ENDIF
      DO i=1,imax
         DO j=1,jmax
            Ac(j,i)=0.D0
         ENDDO
         bc(i)=0.D0
      ENDDO

      RETURN
      END SUBROUTINE mtxc_solve

      SUBROUTINE mtxc_get_vector(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j
      COMPLEX(8),INTENT(OUT):: v
      v=xc(j)
      RETURN
      END SUBROUTINE mtxc_get_vector

      SUBROUTINE mtxc_gather_vector(v)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(imax),INTENT(OUT):: v
      INTEGER:: i

      DO i=1,imax
         v(i)=xc(i)
      ENDDO
      RETURN
      END SUBROUTINE mtxc_gather_vector

      SUBROUTINE mtxc_cleanup
      IMPLICIT NONE

      DEALLOCATE(xc,bc)
      DEALLOCATE(Ac)
      RETURN
      END SUBROUTINE mtxc_cleanup

      END MODULE libmtx
