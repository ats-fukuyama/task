!     $Id$

! -----------------------------------------------------------------------
!
!  Common interface for DMUMPS and ZMUMPS
!
! -----------------------------------------------------------------------

      MODULE libmtx

      USE task_kinds,ONLY: dp
      USE mpi
      USE libmpi
      USE commpi

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

      INCLUDE 'dmumps_struc.h'
      TYPE (DMUMPS_STRUC) id
      REAL(dp),DIMENSION(:),POINTER:: b,b_loc

      INCLUDE 'zmumps_struc.h'
      TYPE (ZMUMPS_STRUC) idc
      COMPLEX(dp),DIMENSION(:),POINTER:: bc,bc_loc

      INTEGER,DIMENSION(:),POINTER:: istartx,iendx,isizex,nz_tot
      INTEGER:: imax,istart,iend,irange,nzcount,nzmax_save,idebug_save

      CONTAINS

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Common section
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      SUBROUTINE mtx_initialize
      IMPLICIT NONE
      INTEGER:: ierr

      CALL MPI_Init(ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_initialize: MPI_Init: ierr=',ierr
      ncomm=MPI_COMM_WORLD
      CALL mtx_set_communicator_global(ncomm)
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

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 DMUMPS section
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      SUBROUTINE mtx_setup(imax_,istart_,iend_,jwidth,nzmax,idebug)

      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(OUT):: istart_,iend_  ! allocated range of lines 
      INTEGER,OPTIONAL,INTENT(IN):: jwidth ! band matrix width
      INTEGER,OPTIONAL,INTENT(IN):: nzmax  ! number of nonzero components
      INTEGER,OPTIONAL,INTENT(IN):: idebug ! debug level
      INTEGER:: i,jmax,iwork1,iwork2

      IF(PRESENT(idebug)) THEN
         idebug_save=idebug
      ELSE
         idebug_save=0
      END IF

!     ----- define a communicator -----      
      id%COMM=ncomm
!     ----- unsymmetric matrix -----      
      id%SYM=0
!     ----- host working -----      
      id%PAR=1
!     ----- initialize MUMPS -----
      id%JOB=-1

      CALL DMUMPS(id)

      imax=imax_
      IF(PRESENT(jwidth)) THEN
         jmax=jwidth
      ELSE
         jmax=2*imax-1
      ENDIF

      iwork1 = imax/nsize
      iwork2 = mod(imax,nsize)
      istart =  nrank   *iwork1 + min(nrank,  iwork2) + 1
      iend   = (nrank+1)*iwork1 + min(nrank+1,iwork2)
      
      istart_=istart
      iend_=iend
      irange=iend-istart+1
      IF(PRESENT(nzmax)) THEN
         nzmax_save=nzmax
      ELSE
         nzmax_save=(iend-istart+1)*jmax
      ENDIF

!      if(nrank.eq.0) then
!         write(6,'(A,3I10)') 'nrank,nsize,nzmax=',nrank,nsize,nzmax_l
!         write(6,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif
!      if(nrank.eq.1) then
!         write(21,'(A,3I10)') 'nrank,nsize,nzmax=',nrank,nsize,nzmax_l
!         write(21,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif
!      if(nrank.eq.2) then
!         write(22,'(A,3I10)') 'nrank,nsize,nzmax=',nrank,nsize,nzmax_l
!         write(22,'(A,3I10)') 'imax,istart,iend=',imax,istart,iend
!      endif
!      if(nrank.eq.3) then
!         write(23,'(A,3I10)') 'nrank,nsize,nzmax=',nrank,nsize,nzmax_l
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
      b_loc(1:iend-istart+1)=0.D0
      RETURN
      END SUBROUTINE mtx_setup
      
!-----

      SUBROUTINE mtx_set_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      REAL(dp),INTENT(IN):: v    ! value to be inserted

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
              'XX libmtxdmumps:mtx_set_matrix: i : out of range'
         write(6,'(A,4I10)') '   nrank,istart,iend,i=',nrank,istart,iend,i
      ENDIF
      return
      END SUBROUTINE mtx_set_matrix
      
!-----

      SUBROUTINE mtx_set_source(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(dp),INTENT(IN):: v ! value to be inserted

      IF(j.GE.istart.AND.j.LE.iend) THEN
         b_loc(j-istart+1)=v
      ELSE
         write(6,'(A)') &
              'XX libmtxdmumps:mtx_set_source: j : out of range'
         write(6,'(A,4I10)') '   nrank,istart,iend,j=',nrank,istart,iend,j
      ENDIF
      RETURN
      END SUBROUTINE mtx_set_source
      
!-----

      SUBROUTINE mtx_set_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(dp),INTENT(IN):: v ! value to be inserted

      return
      END SUBROUTINE mtx_set_vector
      
!-----

      SUBROUTINE mtx_solve(itype,tolerance,its, &
           methodKSP,methodPC,damping_factor,emax,emin,max_steps)
      INTEGER,INTENT(IN):: itype     ! info level
      REAL(dp),INTENT(IN):: tolerance
      INTEGER,INTENT(OUT):: its
      INTEGER,OPTIONAL:: methodKSP,methodPC,max_steps
      REAL(dp),OPTIONAL:: damping_factor,emax,emin
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

      call mtx_allgatherv_real8(b_loc,iend-istart+1,b,imax,isizex,istartx)
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
      CALL dmumps(id)

      CALL mtx_broadcast_real8(id%RHS,imax)
      its=0
      RETURN
      END SUBROUTINE mtx_solve

!-----

      SUBROUTINE mtx_get_vector_j(j,v)

        INTEGER,INTENT(IN):: j
        REAL(dp),INTENT(OUT):: v

        v=id%RHS(j)
      RETURN
      END SUBROUTINE mtx_get_vector_j

!-----

      SUBROUTINE mtx_get_vector(v)

        REAL(dp),DIMENSION(irange),INTENT(OUT):: v
        INTEGER:: j

        DO j=1,irange
           v(j)=id%RHS(istart+j-1)
        ENDDO
        RETURN
      END SUBROUTINE mtx_get_vector

!-----

      SUBROUTINE mtx_gather_vector(v)

      REAL(dp),DIMENSION(imax),INTENT(OUT):: v
      INTEGER:: j

      DO j=1,imax
         v(j)=id%RHS(j)
      ENDDO
      RETURN
      END SUBROUTINE mtx_gather_vector

!-----

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

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 ZMUMPS section
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      SUBROUTINE mtxc_setup(imax_,istart_,iend_,jwidth,nzmax,idebug)

      INTEGER,INTENT(IN):: imax_           ! total matrix size
      INTEGER,INTENT(OUT):: istart_,iend_  ! allocated range of lines 
      INTEGER,OPTIONAL,INTENT(IN):: jwidth ! band matrix width
      INTEGER,OPTIONAL,INTENT(IN):: nzmax  ! number of nonzero components
      INTEGER,OPTIONAL,INTENT(IN):: idebug ! debug level
      INTEGER:: i,iwork1,iwork2

      IF(PRESENT(idebug)) THEN
         idebug_save=idebug
      ELSE
         idebug_save=0
      END IF

!     ----- define a communicator -----      
      idc%COMM=ncomm
!     ----- unsymmetric matrix -----      
      idc%SYM=0
!     ----- host working -----      
      idc%PAR=1
!     ----- initialize MUMPS -----
      idc%JOB=-1

      CALL ZMUMPS(idc)

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
      irange=iend-istart

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

      ALLOCATE(idc%IRN_loc(nzmax_save))
      ALLOCATE(idc%JCN_loc(nzmax_save))
      ALLOCATE(idc%A_loc(nzmax_save))
      ALLOCATE(idc%RHS(imax))
      nzcount=0
      DO i=1,imax
         idc%RHS(i)=0.d0
      ENDDO
      ALLOCATE(istartx(0:nsize-1),iendx(0:nsize-1),isizex(0:nsize-1))
      ALLOCATE(nz_tot(0:nsize-1))
      ALLOCATE(bc(imax),bc_loc(iend-istart+1))
      bc_loc(1:iend-istart+1)=0.D0
      RETURN
      END SUBROUTINE mtxc_setup
      
!-----

      SUBROUTINE mtxc_set_matrix(i,j,v)
      INTEGER,INTENT(IN):: i,j  ! matrix position i=line, j=row
      COMPLEX(dp),INTENT(IN):: v    ! value to be inserted

      IF(i.GE.istart.AND.i.LE.iend) THEN
         nzcount=nzcount+1
         IF(nzcount > nzmax_save) THEN
            WRITE(6,*) "XX mtxc_set_matrix: nzcount > nzmax_save:", &
                       nzcount,nzmax_save
            WRITE(6,*) " at component (i,j): ",i,j
            STOP
         END IF
         idc%A_loc(nzcount)=v
         idc%IRN_loc(nzcount)=i
         idc%JCN_loc(nzcount)=j
      ELSE
         write(6,'(A)') &
              'XX libmtxzmumps:mtxc_set_matrix: i : out of range'
         write(6,'(A,4I10)') '   nrank,istart,iend,i=',nrank,istart,iend,i
      ENDIF
      return
      END SUBROUTINE mtxc_set_matrix
      
!-----

      SUBROUTINE mtxc_set_source(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      COMPLEX(dp),INTENT(IN):: v ! value to be inserted

      IF(j.GE.istart.AND.j.LE.iend) THEN
         bc_loc(j-istart+1)=v
      ELSE
         write(6,'(A)') &
              'XX libmtxzmumps:mtxc_set_source: j : out of range'
         write(6,'(A,4I10)') '   nrank,istart,iend,j=',nrank,istart,iend,j
      ENDIF
      RETURN
      END SUBROUTINE mtxc_set_source
      
!-----

      SUBROUTINE mtxc_set_vector(j,v)
      INTEGER,INTENT(IN):: j ! vector positon j=row
      REAL(dp),INTENT(IN):: v ! value to be inserted

      return
      END SUBROUTINE mtxc_set_vector
      
!-----

      SUBROUTINE mtxc_solve(itype,tolerance,its, &
           methodKSP,methodPC,damping_factor,emax,emin,max_steps)
      INTEGER,INTENT(IN):: itype     ! info level
      REAL(dp),INTENT(IN):: tolerance
      INTEGER,INTENT(OUT):: its
      INTEGER,OPTIONAL:: methodKSP,methodPC,max_steps
      REAL(dp),OPTIONAL:: damping_factor,emax,emin
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

      call mtx_allgatherv_complex8(bc_loc,iend-istart+1,bc,imax,isizex,istartx)
      do i=1,imax
         idc%RHS(i)=bc(i)
      enddo

      idc%N=imax
      idc%NZ=isum
      idc%NZ_loc=nzcount
      idc%NRHS=1
      idc%LRHS=imax
!     ----- distributed matrix -----
      idc%ICNTL(18)=3
!     ----- error output level contrall -----
      IF(itype.EQ.0) THEN
         idc%ICNTL(1)=0
         idc%ICNTL(2)=0
         idc%ICNTL(3)=0
         idc%ICNTL(4)=0
      ELSE IF(itype.EQ.1) THEN
         idc%ICNTL(1)=6
         idc%ICNTL(2)=0
         idc%ICNTL(3)=0
         idc%ICNTL(4)=1
      ELSE IF(itype.GE.2.AND.itype.LE.4) THEN
         idc%ICNTL(1)=6
         idc%ICNTL(2)=0
         idc%ICNTL(3)=6
         idc%ICNTL(4)=itype
      ELSE
         idc%ICNTL(1)=6
         idc%ICNTL(2)=0
         idc%ICNTL(3)=6
         idc%ICNTL(4)=2
      END IF

!      if(nrank.eq.0) then
!         do i=1,nzcount
!            write(21,'(A,2I10,1PE12.4)') 'i,j,A=', &
!                 idc%IRN_loc(i),idc%JCN_loc(i),idc%A_loc(i)
!         enddo
!         do i=1,imax
!            write(21,'(A,I10,1PE12.4)') 'i,bc=',i,idc%RHS(i)
!         enddo
!         do i=istart,iend
!            write(21,'(A,I10,1PE12.4)') 'i,bc_loc=',i,bc_loc(i)
!         enddo
!      endif
!      if(nrank.eq.1) then
!         do i=1,nzcount
!            write(22,'(A,2I10,1PE12.4)') 'i,j,A=', &
!                 idc%IRN_loc(i),idc%JCN_loc(i),idc%A_loc(i)
!         enddo
!         do i=1,imax
!            write(22,'(A,I10,1PE12.4)') 'i,bc=',i,bc(i)
!         enddo
!         do i=istart,iend
!            write(22,'(A,I10,1PE12.4)') 'i,bc_loc=',i,bc_loc(i)
!         enddo
!      endif

      idc%JOB=6
      CALL zmumps(idc)

      CALL mtx_broadcast_complex8(idc%RHS,imax)
      its=0
      RETURN
      END SUBROUTINE mtxc_solve

!-----

      SUBROUTINE mtxc_get_vector_j(j,v)

        INTEGER,INTENT(IN):: j
        COMPLEX(dp),INTENT(OUT):: v

        v=idc%RHS(j)
        RETURN
      END SUBROUTINE mtxc_get_vector_j

!-----

      SUBROUTINE mtxc_get_vector(v)

        COMPLEX(dp),DIMENSION(irange),INTENT(OUT):: v
        INTEGER:: j

        DO j=1,irange
           v(j)=idc%RHS(j+istart-1)
        ENDDO
        RETURN
      END SUBROUTINE mtxc_get_vector

!-----

      SUBROUTINE mtxc_gather_vector(v)

        COMPLEX(dp),DIMENSION(imax),INTENT(OUT):: v
        INTEGER:: j

        DO j=1,imax
           v(j)=idc%RHS(j)
        ENDDO
        RETURN
      END SUBROUTINE mtxc_gather_vector

!-----

      SUBROUTINE mtxc_cleanup

      idc%JOB = -2
      CALL ZMUMPS(idc)

      DEALLOCATE(idc%IRN_loc)
      DEALLOCATE(idc%JCN_loc)
      DEALLOCATE(idc%A_loc)
      DEALLOCATE(idc%RHS)

      DEALLOCATE(istartx,iendx,isizex)
      DEALLOCATE(nz_tot)
      DEALLOCATE(bc,bc_loc)

      RETURN
      END SUBROUTINE mtxc_cleanup

    END MODULE libmtx
