  program mpi_sample1

    implicit none
    include 'mpif.h'

  INTEGER :: ierr, myid, nodes, i, n
  REAL(8) :: total
  REAL(8),DIMENSION(2) :: x
  REAL(8),ALLOCATABLE,DIMENSION(:) :: ista
  ALLOCATE(ista(MPI_STATUS_SIZE))
  x(1) = 1.d0
  x(2) = 1.d0
  CALL mpi_init(ierr)
  CALL mpi_comm_rank(mpi_comm_world,myid,ierr)
  CALL mpi_comm_size(mpi_comm_world,nodes,ierr)
  IF(myid.NE.0) THEN
     CALL MPI_Send(x, 1, MPI_REAL8, 0, myid, MPI_COMM_WORLD, ierr)
  ELSE
     total=x(1)
     DO i=1,nodes-1
        CALL MPI_Recv(x, 1, MPI_REAL8, i, i, MPI_COMM_WORLD, ista, ierr)
        total=total+x(1)
     ENDDO
     WRITE(*,*) total
     CALL mpi_get_count(ista,MPI_REAL8,n,ierr)
     write(*,*) n
  ENDIF
  CALL mpi_finalize(ierr)
  DEALLOCATE(ista)
  end program
