!     $Id$

      MODULE libmtx

      use libmpi

      IMPLICIT NONE

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
      PUBLIC mtx_allgatherv_complex8

      PUBLIC mtx_setup
      PUBLIC mtx_set_matrix
      PUBLIC mtx_set_source
      PUBLIC mtx_set_vector
      PUBLIC mtx_solve
      PUBLIC mtx_get_vector
      PUBLIC mtx_gather_vector
      PUBLIC mtx_cleanup

      PRIVATE

      INCLUDE 'zmumps_struc.h'
      TYPE (ZMUMPS_STRUC) id
      INTEGER,DIMENSION(:),POINTER:: istartx,iendx,isizex,nz_tot
      COMPLEX(8),DIMENSION(:),POINTER:: b,b_loc
      INTEGEr:: imax,istart,iend,jwidth,nzcount

      CONTAINS

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

      CALL ZMUMPS(id)

      imax=imax_
      jwidth=jwidth_

      iwork1 = imax/nsize
      iwork2 = mod(imax,nsize)
      istart =  rank   *iwork1 + min(rank,  iwork2) + 1
      iend   = (rank+1)*iwork1 + min(rank+1,iwork2)
      
      istart_=istart
      iend_=iend
      nzmax=(iend-istart+1)*jwidth

!      if(rank.eq.0) then
!         write(6,'(A,3I10)') 'nsize,rank,nzmax=',nsize,rank,nzmax
!         write(6,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif
!      if(rank.eq.1) then
!         write(21,'(A,3I10)') 'nsize,rank,nzmax=',nsize,rank,nzmax
!         write(21,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif
!      if(rank.eq.2) then
!         write(22,'(A,3I10)') 'nsize,rank,nzmax=',nsize,rank,nzmax
!         write(22,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif
!      if(rank.eq.3) then
!         write(23,'(A,3I10)') 'nsize,rank,nzmax=',nsize,rank,nzmax
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
      ALLOCATE(istartx(0:nsize-1),iendx(0:nsize-1),isizex(0:nsize-1))
      ALLOCATE(nz_tot(0:nsize-1))
      ALLOCATE(b(imax),b_loc(iend-istart+1))
      RETURN
      END SUBROUTINE mtx_setup
      
      SUBROUTINE mtx_set_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      COMPLEX(8),INTENT(IN):: v    ! value to be inserted

      IF(i.GE.istart.AND.i.LE.iend) THEN
         nzcount=nzcount+1
         id%A_loc(nzcount)=v
         id%IRN_loc(nzcount)=i
         id%JCN_loc(nzcount)=j
      ELSE
         write(6,'(A)') &
              'XX libmtxdmumps:mtx_set_matrix: i : out of range'
         write(6,'(A,4I10)') '   rank,istart,iend,i=',rank,istart,iend,i
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
              'XX libmtxdmumps:mtx_set_source: j : out of range'
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

      call mtx_allgather_integer(istart-1,istartx)
      call mtx_allgather_integer(iend,iendx)
      do i=0,nsize-1
         isizex(i)=iendx(i)-istartx(i)
      enddo

      call mtx_allgather_integer(nzcount,nz_tot)
      isum=0
      do i=0,nsize-1
         isum=isum+nz_tot(i)
      enddo

!      if(rank.eq.0) then
!         write(21,'(A,3I10)') 'imax,isum,nzcount=',imax,isum,nzcount
!         do i=0,nsize-1
!            write(21,'(A,4I10)') 'rank,istartx,iendx,isizex=', &
!                  i,istartx(i),iendx(i),isizex(i)
!         enddo
!      endif
!      if(rank.eq.1) then
!         write(22,'(A,3I10)') 'imax,isum,nzcount=',imax,isum,nzcount
!         do i=0,nsize-1
!            write(22,'(A,4I10)') 'rank,istartx,iendx,isizex=', &
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

!      if(rank.eq.0) then
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
!      if(rank.eq.1) then
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

      END MODULE libmtx
