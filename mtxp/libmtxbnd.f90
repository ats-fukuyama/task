!     $Id: libmtxbnd.f90,v 1.20 2013/12/24 09:29:19 fukuyama Exp $

! -----------------------------------------------------------------------
!
!  Description: Solves a linear system with BANDRD/CD (direct method)
!
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

      INTEGER:: imax,jmax,nzzmax
      INTEGER:: joffset,ierr,mode,irmax,icmin,icmax,irc,idebug_save
      REAL(dp),DIMENSION(:),POINTER:: x,b
      REAL(dp),DIMENSION(:,:),POINTER:: A
      INTEGER,DIMENSION(:),POINTER:: ir,ic
      REAL(dp),DIMENSION(:),POINTER:: drc
      COMPLEX(dp),DIMENSION(:),POINTER:: xc,bc
      COMPLEX(dp),DIMENSION(:,:),POINTER:: Ac
      COMPLEX(dp),DIMENSION(:),POINTER:: drcc

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

      SUBROUTINE mtx_setup(imax_,istart_,iend_,jwidth,nzmax,idebug)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(OUT):: istart_,iend_  ! allocated range of lines 
      INTEGER,OPTIONAL,INTENT(IN):: jwidth ! band matrix width
      INTEGER,OPTIONAL,INTENT(IN):: nzmax  ! number of nonzero components
      INTEGER,OPTIONAL,INTENT(IN):: idebug ! debug level
      INTEGER:: i,j

      IF(PRESENT(idebug)) THEN
         idebug_save=idebug
      ELSE
         idebug_save=0
      END IF

      imax=imax_
      IF(PRESENT(jwidth)) THEN
         jmax=jwidth
         mode=1
      ELSE
         IF(PRESENT(nzmax)) THEN
            nzzmax=nzmax
            jmax=nzzmax/imax
         ELSE
            jmax=2*imax_-1
            nzzmax=imax*jmax
         END IF
         mode=2
      END IF
      IF(mode.EQ.1) THEN
         ALLOCATE(A(jmax,imax))
         DO i=1,imax
            DO j=1,jmax
               A(J,i)=0.D0
            ENDDO
         ENDDO
      ELSE
         ALLOCATE(ir(nzzmax),ic(nzzmax),drc(nzzmax))
      END IF
      irmax=0
      icmin=0
      icmax=0
      irc=0

      ALLOCATE(b(imax),x(imax))
      istart_=1
      iend_=imax
      joffset=(jmax+1)/2

      DO i=1,imax
         b(i)=0.D0
         x(i)=0.d0
      ENDDO

      RETURN
      END SUBROUTINE mtx_setup

!-----

      SUBROUTINE mtx_set_matrix(i,j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      REAL(dp),INTENT(IN):: v    ! value to be inserted

      IF(mode.EQ.1) THEN
         A(j-i+joffset,i)=v
      ELSE
         irc=irc+1
         IF(irc.gt.nzzmax) THEN
            WRITE(6,'(A,I10)') &
                 'XX libmtxbnd: mtx_set_matrix: irc overflow: irc=',irc
         ELSE
            ir(irc)=i
            ic(irc)=j
            drc(irc)=v
            IF(idebug_save.EQ.2) &
                 write(6,'(3I5,1PE12.4)') irc,ir(irc),ic(irc),drc(irc)
         END IF
      END IF
      irmax=max(irmax,i)
      icmin=min(icmin,j-i)
      icmax=max(icmax,j-i)
      RETURN
      END SUBROUTINE mtx_set_matrix
      
!-----

      SUBROUTINE mtx_set_source(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(dp),INTENT(IN):: v ! value to be inserted

      b(j)=v
      IF(idebug_save.EQ.2) &
           write(6,'(3I5,1PE12.4)') j,j,0,v
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
      USE libbnd
      IMPLICIT NONE
      INTEGER,INTENT(IN):: itype     ! not used
      REAL(dp),INTENT(IN):: tolerance ! not used
      INTEGER,INTENT(OUT):: its
      INTEGER,OPTIONAL:: methodKSP,methodPC,max_steps
      REAL(dp),OPTIONAL:: damping_factor,emax,emin
      INTEGER:: i,j

      DO i=1,imax
         x(i)=b(i)
      ENDDO

!      do i=1,imax
!         write(21,'(i5/(1P5E12.4))') i,(A(j,i),j=1,jmax)
!      enddo

      IF(MODE.EQ.2) THEN
         jmax=2*max(-icmin,icmax)+1
         joffset=(jmax+1)/2
         ALLOCATE(A(jmax,imax))
         DO i=1,imax
            DO j=1,jmax
               A(J,i)=0.D0
            ENDDO
         ENDDO
         DO i=1,irc
            A(ic(i)-ir(i)+joffset,ir(i))=drc(i)
!            write(6,'(3I5,1PE12.4)') i,ir(i),ic(i),drc(i)
         END DO
!         WRITE(6,'(A,4I8)') &
!              '-- libmtxbnd: mtx_solve: irmax,icmin,icmax,irc=', &
!                                        irmax,icmin,icmax,irc
      END IF
         
      CALL BANDRD(A,x,imax,jmax,jmax,ierr)
      IF(ierr.ne.0) then
         WRITE(6,'(A,I5)') 'XX BANDRD in mtx_solve: ierr=',ierr
         its=-1
      ELSE
         its=0
      ENDIF
      IF(mode.EQ.1) THEN
         DO i=1,imax
            DO j=1,jmax
               A(J,i)=0.D0
            ENDDO
         END DO
      ELSE
         DEALLOCATE(A)
      END IF
      DO i=1,imax
         b(i)=0.D0
      ENDDO

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
      IF(mode.EQ.1) THEN
         DEALLOCATE(A)
      ELSE
         DEALLOCATE(ir,ic,drc)
      END IF
      RETURN
      END SUBROUTINE mtx_cleanup

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of bandcd section
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      SUBROUTINE mtxc_setup(imax_,istart_,iend_,jwidth,nzmax,idebug)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(OUT):: istart_,iend_  ! allocated range of lines 
      INTEGER,OPTIONAL,INTENT(IN):: jwidth ! band matrix width
      INTEGER,OPTIONAL,INTENT(IN):: nzmax  ! number of nonzero components
      INTEGER,OPTIONAL,INTENT(IN):: idebug ! debug level
      INTEGER:: i,j

      IF(PRESENT(idebug)) THEN
         idebug_save=idebug
      ELSE
         idebug_save=0
      END IF

      imax=imax_
      IF(PRESENT(jwidth)) THEN
         jmax=jwidth
         mode=1
      ELSE
         IF(PRESENT(nzmax)) THEN
            nzzmax=nzmax
            jmax=nzmax/imax
         ELSE
            jmax=2*imax_-1
            nzzmax=imax*jmax
         END IF
         mode=2
      END IF
      IF(mode.EQ.1) THEN
         ALLOCATE(Ac(jmax,imax))
         DO i=1,imax
            DO j=1,jmax
               Ac(J,i)=0.D0
            ENDDO
         ENDDO
      ELSE
         ALLOCATE(ir(nzzmax),ic(nzzmax),drcc(nzzmax))
      END IF
      irmax=0
      icmin=0
      icmax=0
      irc=0

      ALLOCATE(bc(imax),xc(imax))
      istart_=1
      iend_=imax
      joffset=(jmax+1)/2

      DO i=1,imax
         bc(i)=0.D0
         xc(i)=0.d0
      ENDDO

      RETURN
      END SUBROUTINE mtxc_setup

!-----

      SUBROUTINE mtxc_set_matrix(i,j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      COMPLEX(dp),INTENT(IN):: v ! value to be inserted

      IF(mode.EQ.1) THEN
         Ac(j-i+joffset,i)=v
      ELSE
         irc=irc+1
         IF(irc.gt.nzzmax) THEN
            WRITE(6,'(A,I10)') &
                 'XX libmtxbnd: mtx_set_matrix: irc overflow: irc=',irc
         ELSE
            ir(irc)=i
            ic(irc)=j
            drcc(irc)=v
            IF(idebug_save.EQ.2) &
                 write(6,'(3I5,1P2E12.4)') irc,ir(irc),ic(irc),drcc(irc)
         END IF
      END IF
      irmax=max(irmax,i)
      icmin=min(icmin,j-i)
      icmax=max(icmax,j-i)
      RETURN
      END SUBROUTINE mtxc_set_matrix
      
!-----

      SUBROUTINE mtxc_set_source(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j    ! vector positon j=row
      COMPLEX(dp),INTENT(IN):: v ! value to be inserted

      bc(j)=v
      IF(idebug_save.EQ.2) &
           write(6,'(3I5,1P2E12.4)') j,j,0,v
      RETURN
      END SUBROUTINE mtxc_set_source
      
!-----

      SUBROUTINE mtxc_set_vector(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j    ! vector positon j=row
      COMPLEX(dp),INTENT(IN):: v ! value to be inserted

      xc(j)=v
      RETURN
      END SUBROUTINE mtxc_set_vector
      
!-----

      SUBROUTINE mtxc_solve(itype,tolerance,its, &
           methodKSP,methodPC,damping_factor,emax,emin,max_steps)
      USE libbnd
      IMPLICIT NONE
      INTEGER,INTENT(IN):: itype     ! not used
      REAL(dp),INTENT(IN):: tolerance ! not used
      INTEGER,INTENT(OUT):: its      ! number of iterations
      INTEGER,OPTIONAL:: methodKSP,methodPC,max_steps
      REAL(dp),OPTIONAL:: damping_factor,emax,emin
      INTEGER:: i,j

      WRITE(6,*) '@@@ point 31'
      DO i=1,imax
         xc(i)=bc(i)
      ENDDO

!      do i=1,imax
!         write(21,'(i5/(1P5E12.4))') i,(A(j,i),j=1,jmax)
!      enddo
         
      WRITE(6,*) '@@@ point 32'
      IF(MODE.EQ.2) THEN
         jmax=2*max(-icmin,icmax)+1
         joffset=(jmax+1)/2
!         write(6,'(A,2I10)') 'imax,jmax=',imax,jmax
         ALLOCATE(Ac(jmax,imax))
         DO i=1,imax
            DO j=1,jmax
               Ac(J,i)=0.D0
            ENDDO
         ENDDO
         DO i=1,irc
            Ac(ic(i)-ir(i)+joffset,ir(i))=drcc(i)
!            write(6,'(3I5,1PE12.4)') i,ir(i),ic(i),drcc(i)
         END DO
!         WRITE(6,'(A,4I8)') &
!              '-- libmtxbnd: mtxc_solve: irmax,icmin,icmax,irc=', &
!                                         irmax,icmin,icmax,irc
      END IF
         
      WRITE(6,*) '@@@ point 33'
      CALL BANDCD(Ac,xc,imax,jmax,jmax,ierr)
      WRITE(6,*) '@@@ point 34'
      IF(ierr.ne.0) then
         WRITE(6,'(A,I5)') 'XX BANDCD in mtxc_solve: ierr=',ierr
         its=-1
      ELSE
         its=0
      ENDIF
      IF(mode.EQ.1) THEN
         DO i=1,imax
            DO j=1,jmax
               Ac(J,i)=0.D0
            ENDDO
         END DO
      ELSE
         DEALLOCATE(Ac)
      END IF
      DO i=1,imax
         bc(i)=0.D0
      ENDDO

      RETURN
      END SUBROUTINE mtxc_solve

      SUBROUTINE mtxc_get_vector_j(j,v)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: j
      COMPLEX(dp),INTENT(OUT):: v
      v=xc(j)
      RETURN
      END SUBROUTINE mtxc_get_vector_j

      SUBROUTINE mtxc_get_vector(v)
      IMPLICIT NONE
      COMPLEX(dp),DIMENSION(imax),INTENT(OUT):: v
      INTEGER:: i

      DO i=1,imax
         v(i)=xc(i)
      ENDDO
      RETURN
      END SUBROUTINE mtxc_get_vector

      SUBROUTINE mtxc_gather_vector(v)
      IMPLICIT NONE
      COMPLEX(dp),DIMENSION(imax),INTENT(OUT):: v
      INTEGER:: i

      DO i=1,imax
         v(i)=xc(i)
      ENDDO
      RETURN
      END SUBROUTINE mtxc_gather_vector

      SUBROUTINE mtxc_cleanup
      IMPLICIT NONE

      IF(mode.EQ.1) THEN
         DEALLOCATE(Ac)
      ELSE
         DEALLOCATE(ir,ic,drcc)
      END IF
      DEALLOCATE(xc,bc)
      RETURN
      END SUBROUTINE mtxc_cleanup

      END MODULE libmtx
