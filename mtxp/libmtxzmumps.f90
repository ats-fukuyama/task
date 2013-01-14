!     $Id$

! -----------------------------------------------------------------------
!
!  Description: Solves a linear system with ZMUMPS (direct method)
!
! -----------------------------------------------------------------------

      MODULE libmtxc

      USE mpi
      USE libmpi
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

      INCLUDE 'zmumps_struc.h'
      TYPE (ZMUMPS_STRUC) id
      INTEGER,DIMENSION(:),POINTER:: istartx,iendx,isizex,nz_tot
      COMPLEX(8),DIMENSION(:),POINTER:: b,b_loc
      INTEGER:: imax,istart,iend,jmax,nzcount,nzmax_save

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
      ncomm=MPI_COMM_WORLD
      CALL mtx_set_communicator_global(ncomm,nrank,nsize)
      nsize_=nsize
      nrank_=nrank
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

      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(OUT):: istart_,iend_  ! allocated range of lines 
      INTEGER,OPTIONAL,INTENT(IN):: jwidth ! band matrix width
      INTEGER,OPTIONAL,INTENT(IN):: nzmax  ! number of nonzero components
      INTEGER:: i,iwork1,iwork2

!     ----- define a communicator -----      
      id%COMM=ncomm
!     ----- unsymmetric matrix -----      
      id%SYM=0
!     ----- host working -----      
      id%PAR=1
!     ----- initialize MUMPS -----
      id%JOB=-1

      CALL ZMUMPS(id)

      imax=imax_
      IF(PRESENT(jwidth)) THEN
         jmax=jwidth
      ELSE
         jmax=2*imax_-1
      ENDIF

      iwork1 = imax/nsize
      iwork2 = mod(imax,nsize)
      istart =  nrank   *iwork1 + min(nrank,  iwork2) + 1
      iend   = (nrank+1)*iwork1 + min(nrank+1,iwork2)
      
      istart_=istart
      iend_=iend

      IF(PRESENT(nzmax)) THEN
         nzmax_save=nzmax
      ELSE
         nzmax_save=(iend-istart+1)*jmax
      ENDIF

!      if(nrank.eq.0) then
!         write(6,'(A,3I10)') 'nrank,nsize,nzmax=',nrank,nsize,nzmax
!         write(6,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif
!      if(nrank.eq.1) then
!         write(21,'(A,3I10)') 'nrank,nsize,nzmax=',nrank,nsize,nzmax
!         write(21,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif
!      if(nrank.eq.2) then
!         write(22,'(A,3I10)') 'nrank,nsize,nzmax=',nrank,nsize,nzmax
!         write(22,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif
!      if(nrank.eq.3) then
!         write(23,'(A,3I10)') 'nrank,nsize,nzmax=',nrank,nsize,nzmax
!         write(23,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif

      ALLOCATE(id%IRN_loc(nzmax_save))
      ALLOCATE(id%JCN_loc(nzmax_save))
      ALLOCATE(id%A_loc(nzmax_save))
      ALLOCATE(id%RHS(imax))
      nzcount=0
      DO i=1,imax
         id%RHS(i)=0.d0
      ENDDO
      ALLOCATE(istartx(0:nsize-1),iendx(0:nsize-1),isizex(0:nsize-1))
      ALLOCATE(nz_tot(0:nsize-1))
      ALLOCATE(b(imax),b_loc(iend-istart+1))
      write(6,*) 'point 3'
      RETURN
      END SUBROUTINE mtx_setup
      
      SUBROUTINE mtx_set_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      COMPLEX(8),INTENT(IN):: v    ! value to be inserted

      IF(i.GE.istart.AND.i.LE.iend) THEN
         nzcount=nzcount+1
         IF(nzcount > nzmax_save) THEN
            WRITE(6,*) "XX mtx_set_matrix: nzcount > nzmax_save:", &
                       nzcount,nzmax_save
            WRITE(6,*) " at component (i,j): ",i,j
            STOP
         END IF
         id%A_loc(nzcount)=v
         id%IRN_loc(nzcount)=i
         id%JCN_loc(nzcount)=j
      ELSE
         write(6,'(A)') &
              'XX libmtxzmumps:mtx_set_matrix: i : out of range'
         write(6,'(A,4I10)') '   nrank,istart,iend,i=',nrank,istart,iend,i
      ENDIF
      return
      END SUBROUTINE mtx_set_matrix
      
      SUBROUTINE mtx_set_source(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      COMPLEX(8),INTENT(IN):: v ! value to be inserted

      IF(j.GE.istart.AND.j.LE.iend) THEN
         b_loc(j-istart+1)=v
      ELSE
         write(6,'(A)') &
              'XX libmtxzmumps:mtx_set_source: j : out of range'
         write(6,'(A,4I10)') '   nrank,istart,iend,j=',nrank,istart,iend,j
      ENDIF
      RETURN
      END SUBROUTINE mtx_set_source
      
      SUBROUTINE mtx_set_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      return
      END SUBROUTINE mtx_set_vector
      
      SUBROUTINE mtx_solve(itype,tolerance,its, &
           methodKSP,methodPC,damping_factor,emax,emin,max_steps)
      INTEGER,INTENT(IN):: itype     ! info level
      REAL(8),INTENT(IN):: tolerance
      INTEGER,INTENT(OUT):: its
      INTEGER,OPTIONAL:: methodKSP,methodPC,max_steps
      REAL(8),OPTIONAL:: damping_factor,emax,emin
      INTEGER:: i,isum

      call mtx_allgather1_integer(istart-1,istartx)
      call mtx_allgather1_integer(iend,iendx)
      do i=0,nsize-1
         isizex(i)=iendx(i)-istartx(i)
      enddo

      call mtx_allgather1_integer(nzcount,nz_tot)
      isum=0
      do i=0,nsize-1
         isum=isum+nz_tot(i)
      enddo

!      if(nrank.eq.0) then
!         write(21,'(A,3I10)') 'imax,isum,nzcount=',imax,isum,nzcount
!         do i=0,nsize-1
!            write(21,'(A,4I10)') 'nrank,istartx,iendx,isizex=', &
!                  i,istartx(i),iendx(i),isizex(i)
!         enddo
!      endif
!      if(nrank.eq.1) then
!         write(22,'(A,3I10)') 'imax,isum,nzcount=',imax,isum,nzcount
!         do i=0,nsize-1
!            write(22,'(A,4I10)') 'nrank,istartx,iendx,isizex=', &
!                  i,istartx(i),iendx(i),isizex(i)
!         enddo
!      endif

      call mtx_allgatherv_complex8(b_loc,iend-istart+1,b,imax,isizex,istartx)
      do i=1,imax
         id%RHS(i)=b(i)
      enddo

      id%N=imax
      id%NZ=isum
      id%NZ_loc=nzcount
      id%NRHS=1
      id%LRHS=imax
!     ----- distributed matrix -----
      id%ICNTL(18)=3
!     ----- error output level contrall -----
      IF(itype.EQ.0) THEN
         id%ICNTL(1)=0
         id%ICNTL(2)=0
         id%ICNTL(3)=0
         id%ICNTL(4)=0
      ELSE IF(itype.EQ.1) THEN
         id%ICNTL(1)=6
         id%ICNTL(2)=0
         id%ICNTL(3)=0
         id%ICNTL(4)=1
      ELSE IF(itype.GE.2.AND.itype.LE.4) THEN
         id%ICNTL(1)=6
         id%ICNTL(2)=0
         id%ICNTL(3)=6
         id%ICNTL(4)=itype
      ELSE
         id%ICNTL(1)=6
         id%ICNTL(2)=0
         id%ICNTL(3)=6
         id%ICNTL(4)=2
      END IF

!      if(nrank.eq.0) then
!         do i=1,nzcount
!            write(21,'(A,2I10,1PE12.4)') 'i,j,A=', &
!                 id%IRN_loc(i),id%JCN_loc(i),id%A_loc(i)
!         enddo
!         do i=1,imax
!            write(21,'(A,I10,1PE12.4)') 'i,b=',i,id%RHS(i)
!         enddo
!         do i=istart,iend
!            write(21,'(A,I10,1PE12.4)') 'i,b_loc=',i,b_loc(i)
!         enddo
!      endif
!      if(nrank.eq.1) then
!         do i=1,nzcount
!            write(22,'(A,2I10,1PE12.4)') 'i,j,A=', &
!                 id%IRN_loc(i),id%JCN_loc(i),id%A_loc(i)
!         enddo
!         do i=1,imax
!            write(22,'(A,I10,1PE12.4)') 'i,b=',i,b(i)
!         enddo
!         do i=istart,iend
!            write(22,'(A,I10,1PE12.4)') 'i,b_loc=',i,b_loc(i)
!         enddo
!      endif

      id%JOB=6
      CALL zmumps(id)

      its=0
      RETURN
      END SUBROUTINE mtx_solve

      SUBROUTINE mtx_get_vector(j,v)

      INTEGER,INTENT(IN):: j
      COMPLEX(8),INTENT(OUT):: v

      v=id%RHS(j)
      RETURN
      END SUBROUTINE mtx_get_vector

      SUBROUTINE mtx_gather_vector(v)

      COMPLEX(8),DIMENSION(imax),INTENT(OUT):: v
      INTEGER:: j

      DO j=1,imax
         v(j)=id%RHS(j)
      ENDDO
      RETURN
      END SUBROUTINE mtx_gather_vector

      SUBROUTINE mtx_cleanup

      id%JOB = -2
      CALL ZMUMPS(id)

      DEALLOCATE(id%IRN_loc)
      DEALLOCATE(id%JCN_loc)
      DEALLOCATE(id%A_loc)
      DEALLOCATE(id%RHS)

      DEALLOCATE(istartx,iendx,isizex)
      DEALLOCATE(nz_tot)
      DEALLOCATE(b,b_loc)

      RETURN
      END SUBROUTINE mtx_cleanup

      END MODULE libmtxc
