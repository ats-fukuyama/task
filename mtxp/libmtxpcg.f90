!     $Id$

! -----------------------------------------------------------------------

!  Description: Solves a linear system with PCGPME method (iterative method)

! -----------------------------------------------------------------------

      MODULE libmtx

      USE task_kinds,ONLY: dp
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
      PUBLIC mtx_get_vector_j
      PUBLIC mtx_get_vector
      PUBLIC mtx_gather_vector
      PUBLIC mtx_cleanup

      PUBLIC mtxc_setup
      PUBLIC mtxc_set_matrix
      PUBLIC mtxc_set_source
      PUBLIC mtxc_set_vector
      PUBLIC mtxc_solve
      PUBLIC mtxc_get_vector_j
      PUBLIC mtxc_get_vector
      PUBLIC mtxc_gather_vector
      PUBLIC mtxc_cleanup

      TYPE(mtx_mpi_type):: mtx_global

      INTEGER:: imax,jmax,nzmax_save,nzcount,idebug_save
      INTEGER,DIMENSION(:),ALLOCATABLE:: iC,jC
      REAL(dp),DIMENSION(:),ALLOCATABLE:: x,b
      REAL(dp),DIMENSION(:),ALLOCATABLE:: vC

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
!                 Beginning of pcgpme section
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE mtx_setup(imax_,istart_,iend_,jwidth,nzmax,idebug)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(OUT):: istart_,iend_  ! allocated range of lines 
      INTEGER,OPTIONAL,INTENT(IN):: jwidth ! band matrix width
      INTEGER,OPTIONAL,INTENT(IN):: nzmax  ! number of nonzero components
      INTEGER,OPTIONAL,INTENT(IN):: idebug ! debug level
      INTEGER:: i

      IF(PRESENT(idebug)) THEN
         idebug_save=idebug
      ELSE
         idebug_save=0
      END IF

      imax=imax_
      IF(PRESENT(nzmax)) THEN
         nzmax_save=nzmax
      ELSE
         nzmax_save=imax*27
      ENDIF
      ALLOCATE(vC(nzmax_save),iC(nzmax_save),jC(nzmax_save))
      ALLOCATE(b(imax),x(imax))
      nzcount=0
      istart_=1
      iend_=imax
      DO i=1,imax
         b(i)=0.D0
         x(i)=0.D0
      END DO
      RETURN
      END SUBROUTINE mtx_setup

!-----

      SUBROUTINE mtx_set_matrix(i,j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      REAL(dp),INTENT(IN):: v    ! value to be inserted

      IF(ABS(v) > 0.D0) THEN
         nzcount=nzcount+1
         IF(nzcount > nzmax_save) THEN
            WRITE(6,'(A,I6)') 'XX mtx_set_matrix: nzcount > nzmax: ',nzcount
            STOP
         END IF
         iC(nzcount)=i
         jC(nzcount)=j
         vC(nzcount)=v
      END IF
            
      RETURN
      END SUBROUTINE mtx_set_matrix
      
!-----

      SUBROUTINE mtx_set_source(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(dp),INTENT(IN):: v ! value to be inserted

      b(j)=v
!      write(18,'(I5,1P2E12.4)') j,v,b(j)
      RETURN
      END SUBROUTINE mtx_set_source
      
!-----

      SUBROUTINE mtx_set_vector(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(dp),INTENT(IN):: v ! value to be inserted

      x(j)=v
      RETURN
      END SUBROUTINE mtx_set_vector
      
      SUBROUTINE mtx_solve(itype,tolerance,its, &
           methodKSP,methodPC,damping_factor,emax,emin,max_steps)
      USE libpcgpme
      IMPLICIT NONE
      INTEGER,INTENT(IN):: itype     ! not used
      REAL(dp),INTENT(IN):: tolerance ! tolerance
      INTEGER,INTENT(OUT):: its      ! number of iteration
      INTEGER,OPTIONAL:: max_steps   ! maximum number of iteration
      INTEGER,OPTIONAL:: methodKSP,methodPC
      REAL(dp),OPTIONAL:: damping_factor,emax,emin
      REAL(dp),DIMENSION(:),ALLOCATABLE:: D,xin
      INTEGER,DIMENSION(:),ALLOCATABLE:: NL
      REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AL
      INTEGER,DIMENSION(:,:),ALLOCATABLE:: LL
      INTEGER:: i,j,nz,NLMAX,itrin,IER
      REAL(dp):: EPSOUT

      ALLOCATE(D(imax),NL(imax),xin(imax))
      DO i=1,imax
         xin(i)=x(i)
      END DO

!     calculate maximum number of each row

      DO i=1,imax
         NL(i)=0
      END DO
      DO nz=1,nzcount
         i=iC(nz)
         j=jC(nz)
         IF(i /= j) NL(i)=NL(i)+1
      END DO
      NLMAX=NL(1)
      DO i=2,imax
         IF(NL(i) > NLMAX) NLMAX=NL(i)
      END DO

!     define max_stetps

      IF(PRESENT(max_steps)) THEN
         itrin=max_steps
      ELSE
         itrin=MAX(2*NLMAX,imax)
      END IF

!     create coefficient array

      ALLOCATE(AL(imax,NLMAX),LL(imax,NLMAX))
      DO i=1,imax
         NL(i)=0
         DO j=1,NLMAX
            LL(i,j)=0
            AL(i,j)=0.D0
         END DO
      END DO
      DO nz=1,nzcount
         i=iC(nz)
         j=jC(nz)
         IF(i == j) THEN
            D(i)=vC(nz)
         ELSE
            NL(i)=NL(i)+1
            AL(i,NL(i))=vC(nz)
            LL(i,NL(i))=j
         END IF
      END DO

!      DO i=1,imax
!         write(16,'(I5,5(I5,1PE12.4))') &
!              i,NL(i),D(i),(LL(i,j),AL(i,j),j=1,NL(i))
!      END DO
!      DO i=1,imax
!         write(17,'(I5,1P2E12.4)') i,b(i),x(i)
!      END DO

      CALL PCGPME(AL,NLMAX,imax,LL,D,imax,64,b,xin,tolerance,itrin, &
                  x,EPSOUT,its,IER)
      IF(IER /= 0) &
           WRITE(6,'(A,I5)') 'XX pcgpme: IER =',IER
      IF(idebug_save.GE.1) &
           write(6,'(A,I5,A,1PE12.4)') &
           '-- pcgpme: iteration =',its,'  convergence =',EPSOUT
      
      DEALLOCATE(AL,LL,D,NL,xin)
         
      RETURN
      END SUBROUTINE mtx_solve

!-----

      SUBROUTINE mtx_get_vector_j(j,v)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: j
        REAL(dp),INTENT(OUT):: v
        v=x(j)
        RETURN
      END SUBROUTINE mtx_get_vector_j

!-----

      SUBROUTINE mtx_get_vector(v)
        IMPLICIT NONE
        REAL(dp),DIMENSION(imax),INTENT(OUT):: v
        INTEGER:: i

        DO i=1,imax
           v(i)=x(i)
        ENDDO
        RETURN
      END SUBROUTINE mtx_get_vector

!-----

      SUBROUTINE mtx_gather_vector(v)
        IMPLICIT NONE
        REAL(dp),DIMENSION(imax),INTENT(OUT):: v
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
      DEALLOCATE(iC,jC,vC)
      RETURN
      END SUBROUTINE mtx_cleanup

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Complex inteface for compatibility
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      SUBROUTINE mtxc_setup(imax_,istart_,iend_,jwidth,nzmax,idebug)

      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(OUT):: istart_,iend_  ! allocated range of lines 
      INTEGER,OPTIONAL,INTENT(IN):: jwidth ! band matrix width
      INTEGER,OPTIONAL,INTENT(IN):: nzmax  ! number of nonzero components
      INTEGER,OPTIONAL,INTENT(IN):: idebug ! debug level

      istart_=0
      iend_=0
      RETURN
      END SUBROUTINE mtxc_setup
      
!-----

      SUBROUTINE mtxc_set_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      COMPLEX(dp),INTENT(IN):: v    ! value to be inserted

      RETURN
      END SUBROUTINE mtxc_set_matrix
      
!-----

      SUBROUTINE mtxc_set_source(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      COMPLEX(dp),INTENT(IN):: v ! value to be inserted

      RETURN
      END SUBROUTINE mtxc_set_source
      
!-----

      SUBROUTINE mtxc_set_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(dp),INTENT(IN):: v ! value to be inserted

      RETURN
      END SUBROUTINE mtxc_set_vector
      
!-----

      SUBROUTINE mtxc_solve(itype,tolerance,its, &
           methodKSP,methodPC,damping_factor,emax,emin,max_steps)
      INTEGER,INTENT(IN):: itype     ! info level
      REAL(dp),INTENT(IN):: tolerance
      INTEGER,INTENT(OUT):: its
      INTEGER,OPTIONAL:: methodKSP,methodPC,max_steps
      REAL(dp),OPTIONAL:: damping_factor,emax,emin

      its=0
      RETURN
      END SUBROUTINE mtxc_solve

!-----

      SUBROUTINE mtxc_get_vector_j(j,v)
        INTEGER,INTENT(IN):: j
        COMPLEX(dp),INTENT(OUT):: v
        v=0.D0
        RETURN
      END SUBROUTINE mtxc_get_vector_j

!-----

      SUBROUTINE mtxc_get_vector(v)
        COMPLEX(dp),DIMENSION(imax),INTENT(OUT):: v
        v(1:imax)=0.D0
        RETURN
      END SUBROUTINE mtxc_get_vector

!-----

      SUBROUTINE mtxc_gather_vector(v)
        COMPLEX(dp),DIMENSION(imax),INTENT(OUT):: v
        v(1:imax)=0.D0
        RETURN
      END SUBROUTINE mtxc_gather_vector

!-----

      SUBROUTINE mtxc_cleanup

      RETURN
      END SUBROUTINE mtxc_cleanup

      END MODULE libmtx
