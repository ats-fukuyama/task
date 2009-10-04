!     $Id$

      MODULE libmtx

      IMPLICIT NONE
      PUBLIC mtx_initialize
      PUBLIC mtx_setup
      PUBLIC mtx_set_matrix
      PUBLIC mtx_set_source
      PUBLIC mtx_set_vector
      PUBLIC mtx_solve
      PUBLIC mtx_get_vector
      PUBLIC mtx_gather_vector
      PUBLIC mtx_cleanup
      PUBLIC mtx_finalize
      PUBLIC mtx_barrier
      PUBLIC mtx_broadcast_character   ! broadcast from rank=0
      PUBLIC mtx_broadcast_integer
      PUBLIC mtx_broadcast_real8
      PUBLIC mtx_broadcast_complex8
      PUBLIC mtx_gather_integer
      PUBLIC mtx_allgather_integer
      PUBLIC mtx_gatherv_real8
      PUBLIC mtx_allgatherv_real8
      PRIVATE

#include "finclude/petsc.h"
      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC) id
      INTEGER,DIMENSION(:),POINTER:: istartx,iendx,isizex,nz_tot
      REAL(8),DIMENSION(:),POINTER:: b,b_loc
      INTEGER:: rank,size,ierr
      INTEGEr:: imax,istart,iend,jwidth,nzcount

      CONTAINS

      SUBROUTINE mtx_initialize(rank_,size_)
      INTEGER,INTENT(OUT):: rank_,size_

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
      call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
      rank_=rank
      size_=size
!      if(rank.eq.0) write(6,'(A,I10)') '# mtx_initialize: size=',size
      return

      END SUBROUTINE mtx_initialize

      SUBROUTINE mtx_setup(imax_,istart_,iend_,jwidth_)

      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(IN):: jwidth_         ! band matrix width
      INTEGER,INTENT(OUT):: istart_,iend_  ! allocated range of lines 
      INTEGER:: i,nzmax,iwork1,iwork2

!     ----- define a communicator -----      
      id%COMM=PETSC_COMM_WORLD
!     ----- unsymmetric matrix -----      
      id%SYM=0
!     ----- host working -----      
      id%PAR=1
!     ----- initialize MUMPS -----
      id%JOB=-1

      CALL DMUMPS(id)

      imax=imax_
      jwidth=jwidth_

      iwork1 = imax/size
      iwork2 = mod(imax,size)
      istart =  rank   *iwork1 + min(rank,  iwork2) + 1
      iend   = (rank+1)*iwork1 + min(rank+1,iwork2)
      
      istart_=istart
      iend_=iend
      nzmax=(iend-istart+1)*jwidth

!      if(rank.eq.0) then
!         write(6,'(A,3I10)') 'size,rank,nzmax=',size,rank,nzmax
!         write(6,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif
!      if(rank.eq.1) then
!         write(21,'(A,3I10)') 'size,rank,nzmax=',size,rank,nzmax
!         write(21,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif
!      if(rank.eq.2) then
!         write(22,'(A,3I10)') 'size,rank,nzmax=',size,rank,nzmax
!         write(22,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif
!      if(rank.eq.3) then
!         write(23,'(A,3I10)') 'size,rank,nzmax=',size,rank,nzmax
!         write(23,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif

      ALLOCATE(id%IRN_loc(nzmax))
      ALLOCATE(id%JCN_loc(nzmax))
      ALLOCATE(id%A_loc(nzmax))
      ALLOCATE(id%RHS(imax))
      nzcount=0
      DO i=1,imax
         id%RHS(i)=0.d0
      ENDDO
      ALLOCATE(istartx(0:size-1),iendx(0:size-1),isizex(0:size-1))
      ALLOCATE(nz_tot(0:size-1))
      ALLOCATE(b(imax),b_loc(iend-istart+1))
      RETURN
      END SUBROUTINE mtx_setup
      
      SUBROUTINE mtx_set_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      REAL(8),INTENT(IN):: v    ! value to be inserted

      IF(i.GE.istart.AND.i.LE.iend) THEN
         nzcount=nzcount+1
         id%A_loc(nzcount)=v
         id%IRN_loc(nzcount)=i
         id%JCN_loc(nzcount)=j
      ELSE
         write(6,'(A)') 
     &        'XX libmtxdmumps:mtx_set_matrix: i : out of range'
         write(6,'(A,4I10)') '   rank,istart,iend,i=',rank,istart,iend,i
      ENDIF
      return
      END SUBROUTINE mtx_set_matrix
      
      SUBROUTINE mtx_set_source(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      IF(j.GE.istart.AND.j.LE.iend) THEN
         b_loc(j-istart+1)=v
      ELSE
         write(6,'(A)') 
     &        'XX libmtxdmumps:mtx_set_source: j : out of range'
         write(6,'(A,4I10)') '   rank,istart,iend,j=',rank,istart,iend,j
      ENDIF
      RETURN
      END SUBROUTINE mtx_set_source
      
      SUBROUTINE mtx_set_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(8),INTENT(IN):: v ! value to be inserted

      return
      END SUBROUTINE mtx_set_vector
      
      SUBROUTINE mtx_solve(itype,tolerance,its)
      INTEGER,INTENT(IN):: itype     ! info level
      REAL(8),INTENT(IN):: tolerance
      INTEGER,INTENT(OUT):: its
      INTEGER:: i,isum

      call mtx_allgather_integer(istart-1,istartx,size)
      call mtx_allgather_integer(iend,iendx,size)
      do i=0,size-1
         isizex(i)=iendx(i)-istartx(i)
      enddo

      call mtx_allgather_integer(nzcount,nz_tot,size)
      isum=0
      do i=0,size-1
         isum=isum+nz_tot(i)
      enddo

!      if(rank.eq.0) then
!         write(21,'(A,3I10)') 'imax,isum,nzcount=',imax,isum,nzcount
!         do i=0,size-1
!            write(21,'(A,4I10)') 'rank,istartx,iendx,isizex=',
!     &            i,istartx(i),iendx(i),isizex(i)
!         enddo
!      endif
!      if(rank.eq.1) then
!         write(22,'(A,3I10)') 'imax,isum,nzcount=',imax,isum,nzcount
!         do i=0,size-1
!            write(22,'(A,4I10)') 'rank,istartx,iendx,isizex=',
!     &            i,istartx(i),iendx(i),isizex(i)
!         enddo
!      endif

      call mtx_allgatherv_real8(b_loc,iend-istart+1,b,imax,
     &                          isizex,istartx)
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
      IF(itype.GE.1.AND.itype.LE.4) THEN
         id%ICNTL(4)=itype
      ELSE
         id%ICNTL(4)=0
      ENDIF

!      if(rank.eq.0) then
!         do i=1,nzcount
!            write(21,'(A,2I10,1PE12.4)') 'i,j,A=', 
!     &           id%IRN_loc(i),id%JCN_loc(i),id%A_loc(i)
!         enddo
!         do i=1,imax
!            write(21,'(A,I10,1PE12.4)') 'i,b=',i,id%RHS(i)
!         enddo
!         do i=istart,iend
!            write(21,'(A,I10,1PE12.4)') 'i,b_loc=',i,b_loc(i)
!         enddo
!      endif
!      if(rank.eq.1) then
!         do i=1,nzcount
!            write(22,'(A,2I10,1PE12.4)') 'i,j,A=', 
!     &           id%IRN_loc(i),id%JCN_loc(i),id%A_loc(i)
!         enddo
!         do i=1,imax
!            write(22,'(A,I10,1PE12.4)') 'i,b=',i,b(i)
!         enddo
!         do i=istart,iend
!            write(22,'(A,I10,1PE12.4)') 'i,b_loc=',i,b_loc(i)
!         enddo
!      endif

      id%JOB=6
      CALL dmumps(id)

      its=0
      RETURN
      END SUBROUTINE mtx_solve

      SUBROUTINE mtx_get_vector(j,v)

      INTEGER,INTENT(IN):: j
      REAL(8),INTENT(OUT):: v

      v=id%RHS(j)
      RETURN
      END SUBROUTINE mtx_get_vector

      SUBROUTINE mtx_gather_vector(v)

      REAL(8),DIMENSION(imax),INTENT(OUT):: v
      INTEGER:: j

      DO j=1,imax
         v(j)=id%RHS(j)
      ENDDO
      RETURN
      END SUBROUTINE mtx_gather_vector

      SUBROUTINE mtx_cleanup

      id%JOB = -2
      CALL DMUMPS(id)

      DEALLOCATE(id%IRN_loc)
      DEALLOCATE(id%JCN_loc)
      DEALLOCATE(id%A_loc)
      DEALLOCATE(id%RHS)

      DEALLOCATE(istartx,iendx,isizex)
      DEALLOCATE(nz_tot)
      DEALLOCATE(b,b_loc)

      RETURN
      END SUBROUTINE mtx_cleanup

      SUBROUTINE mtx_finalize

      call PetscFinalize(ierr)
      END SUBROUTINE mtx_finalize

      SUBROUTINE mtx_barrier
      call MPI_Barrier(PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_barrier

      SUBROUTINE mtx_broadcast_character(kdata,n)

      CHARACTER(LEN=n),INTENT(INOUT):: kdata
      INTEGER,INTENT(IN):: n
      
      call mpi_bcast(kdata,n,MPI_CHARACTER,
     &               0,PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_broadcast_character

      SUBROUTINE mtx_broadcast_integer(idata,n)

      INTEGER,DIMENSION(n),INTENT(INOUT):: idata
      INTEGER,INTENT(IN):: n
      
      call mpi_bcast(idata,n,MPI_INTEGER,
     &               0,PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_broadcast_integer

      SUBROUTINE mtx_broadcast_real8(vdata,n)

      REAL(8),DIMENSION(n),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n
      
      call mpi_bcast(vdata,n,MPI_DOUBLE_PRECISION,
     &               0,PETSC_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE mtx_broadcast_real8

      SUBROUTINE mtx_broadcast_complex8(vdata,n)

      COMPLEX(8),DIMENSION(n),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n
      
      call mpi_bcast(vdata,n,MPI_DOUBLE_COMPLEX,
     &               0,PETSC_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE mtx_broadcast_complex8

      SUBROUTINE mtx_gather_integer(idata,itot,ntot)

      INTEGER,INTENT(IN):: idata
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(ntot),INTENT(OUT):: itot

      call MPI_GATHER(idata,1,MPI_INTEGER,
     &                itot,1,MPI_INTEGER,
     &                0,PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_gather_integer

      SUBROUTINE mtx_allgather_integer(idata,itot,ntot)

      INTEGER,INTENT(IN):: idata
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(ntot),INTENT(OUT):: itot

      call MPI_ALLGATHER(idata,1,MPI_INTEGER,
     &                   itot,1,MPI_INTEGER,
     &                   PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_allgather_integer

      SUBROUTINE mtx_gatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa)

      INTEGER,INTENT(IN):: ndata
      INTEGER,INTENT(IN):: ntot
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,DIMENSION(size),INTENT(IN):: ilena,iposa

      call MPI_GATHERV(vdata,ndata,MPI_DOUBLE_PRECISION,
     &                 vtot,ilena,iposa,MPI_DOUBLE_PRECISION,
     &                 0,PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_gatherv_real8

      SUBROUTINE mtx_allgatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa)

      INTEGER,INTENT(IN):: ndata
      INTEGER,INTENT(IN):: ntot
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,DIMENSION(size),INTENT(IN):: ilena,iposa

      call MPI_ALLGATHERV(vdata,ndata,MPI_DOUBLE_PRECISION,
     &                    vtot,ilena,iposa,MPI_DOUBLE_PRECISION,
     &                    PETSC_COMM_WORLD,ierr)
      RETURN
      END SUBROUTINE mtx_allgatherv_real8

      END MODULE libmtx
