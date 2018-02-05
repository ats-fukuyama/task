!     $Id$

      MODULE commpi
!
!       This module indicates the status of mpi communicator
!       Do not change the variables in this module directly.
!       In order to change the communicator, 
!       use the routine "mtx_set_communicator" 
!
        INTEGER,PUBLIC:: ncomm,nrank,nsize
      END MODULE commpi
        

      MODULE libmpi

      USE mpi
      PRIVATE

      TYPE mtx_mpi_type
         integer:: comm   ! communicator
         integer:: rank   ! rank of processor
         integer:: size   ! number of processors
         integer:: rankg  ! rank of the group (color)
         integer:: sizeg  ! number of groups
         integer:: rankl  ! rank of the processors in the group (key) = rank
         integer:: sizel  ! number of processors in the group = size
      END TYPE mtx_mpi_type
      PUBLIC mtx_mpi_type

      PUBLIC mtx_mpi

      PUBLIC mtx_set_communicator_global
      PUBLIC mtx_set_communicator
      PUBLIC mtx_reset_communicator
      PUBLIC mtx_comm_split2D
      PUBLIC mtx_comm_split3D
      PUBLIC mtx_comm_free
      PUBLIC mtx_barrier

      PUBLIC mtx_broadcast1_character
      PUBLIC mtx_broadcast1_logical
      PUBLIC mtx_broadcast1_integer
      PUBLIC mtx_broadcast1_real4
      PUBLIC mtx_broadcast1_real8
      PUBLIC mtx_broadcast1_complex8
      PUBLIC mtx_broadcast_character
      PUBLIC mtx_broadcast_logical
      PUBLIC mtx_broadcast_integer
      PUBLIC mtx_broadcast_real4
      PUBLIC mtx_broadcast_real8
      PUBLIC mtx_broadcast_complex8
      PUBLIC mtx_broadcast2D_integer
      PUBLIC mtx_broadcast2D_real4
      PUBLIC mtx_broadcast2D_real8
      PUBLIC mtx_broadcast2D_complex8

      PUBLIC mtx_gather1_integer
      PUBLIC mtx_gather1_real4
      PUBLIC mtx_gather1_real8
      PUBLIC mtx_gather1_complex8
      PUBLIC mtx_gather_integer
      PUBLIC mtx_gather_real4
      PUBLIC mtx_gather_real8
      PUBLIC mtx_gather_complex8
      PUBLIC mtx_gatherv_integer
      PUBLIC mtx_gatherv_real4
      PUBLIC mtx_gatherv_real8
      PUBLIC mtx_gatherv_complex8
      PUBLIC mtx_allgather1_integer
      PUBLIC mtx_allgather1_real4
      PUBLIC mtx_allgather1_real8
      PUBLIC mtx_allgather1_complex8
      PUBLIC mtx_allgather_integer
      PUBLIC mtx_allgather_real4
      PUBLIC mtx_allgather_real8
      PUBLIC mtx_allgather_complex8
      PUBLIC mtx_allgatherv_integer
      PUBLIC mtx_allgatherv_real4
      PUBLIC mtx_allgatherv_real8
      PUBLIC mtx_allgatherv_complex8

      PUBLIC mtx_reduce1_integer
      PUBLIC mtx_reduce1_real4
      PUBLIC mtx_reduce1_real8
      PUBLIC mtx_reduce1_complex8
      PUBLIC mtx_reduce_integer
      PUBLIC mtx_reduce_real4
      PUBLIC mtx_reduce_real8
      PUBLIC mtx_reduce_complex8
      PUBLIC mtx_allreduce1_integer
      PUBLIC mtx_allreduce1_real4
      PUBLIC mtx_allreduce1_real8
      PUBLIC mtx_allreduce1_complex8
      PUBLIC mtx_allreduce_integer
      PUBLIC mtx_allreduce_real4
      PUBLIC mtx_allreduce_real8
      PUBLIC mtx_allreduce_complex8

      PUBLIC mtx_sendrecv_integer
      PUBLIC mtx_sendrecv_real4
      PUBLIC mtx_sendrecv_real8
      PUBLIC mtx_sendrecv_complex8

      PUBLIC mtx_scatterv_real8
      PUBLIC mtx_send_real8
      PUBLIC mtx_recv_real8

      PUBLIC mtx_abort

      TYPE(mtx_mpi_type):: mtx_global
      INTEGER:: ncomm,nrank,nsize

      CONTAINS

!-----

      SUBROUTINE mtx_mpi(ncomm_,nrank_,nsize_)
        USE commpi
        IMPLICIT NONE
        INTEGER,INTENT(IN):: ncomm_,nrank_,nsize_

        ncomm=ncomm_
        nrank=nrank_
        nsize=nsize_
        return
      END SUBROUTINE mtx_mpi

!-----

      SUBROUTINE mtx_set_communicator_global(ncomm_in)
        USE commpi,ncomm_=>ncomm,nrank_=>nrank,nsize_=>nsize
        IMPLICIT NONE
        INTEGER,INTENT(IN):: ncomm_in
        INTEGER:: ierr

        ncomm=ncomm_in
        call MPI_Comm_rank(ncomm,nrank,ierr)
        IF(ierr.NE.0) WRITE(6,*) &
             'XX mtx_set_communicator_global: MPI_Comm_rank: ierr=',ierr
        call MPI_Comm_size(ncomm,nsize,ierr)
        IF(ierr.NE.0) WRITE(6,*) &
             'XX mtx_set_communicator_global: MPI_Comm_size: ierr=',ierr
        mtx_global%comm=ncomm
        mtx_global%rank=nrank
        mtx_global%size=nsize
        mtx_global%rankg=0
        mtx_global%sizeg=1
        mtx_global%rankl=nrank
        mtx_global%sizel=nsize
        ncomm_=ncomm
        nrank_=nrank
        nsize_=nsize
        return
      END SUBROUTINE mtx_set_communicator_global

!-----

      SUBROUTINE mtx_set_communicator(mtx_mpi)
        USE commpi,ncomm_=>ncomm,nrank_=>nrank,nsize_=>nsize
        IMPLICIT NONE
        TYPE(mtx_mpi_type),INTENT(IN):: mtx_mpi

        ncomm=mtx_mpi%comm
        nrank=mtx_mpi%rank
        nsize=mtx_mpi%size
        ncomm_=ncomm
        nrank_=nrank
        nsize_=nsize
        return
      END SUBROUTINE mtx_set_communicator

!-----

      SUBROUTINE mtx_reset_communicator
        USE commpi,ncomm_=>ncomm,nrank_=>nrank,nsize_=>nsize
        IMPLICIT NONE

        ncomm=mtx_global%comm
        nrank=mtx_global%rank
        nsize=mtx_global%size
        ncomm_=ncomm
        nrank_=nrank
        nsize_=nsize
        return
      END SUBROUTINE mtx_reset_communicator

!-----

      SUBROUTINE mtx_comm_split2D(n1,n2,commx1,commx2)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: n1,n2 ! number of groups
        TYPE(mtx_mpi_type),intent(OUT):: commx1,commx2
        integer:: ierr

        IF(n1*n2 > nsize) THEN
           WRITE(6,*) 'XX mtx_comm_split2D: n1*n2 > nsize: ',n1,n2,nsize
           STOP
        ENDIF
        commx1%sizeg=n1       ! number of groups
        commx1%sizel=n2       ! number of processors in the group
        commx1%rankg=nrank/n2 ! colors
        commx1%rankl=nrank - commx1%rankg*n2 ! keys

        commx2%sizeg=n2
        commx2%sizel=n1
        commx2%rankg=MOD(nrank+1,n2) ! colorr
        commx2%rankl=nrank/n2        ! keyr

        CALL MPI_Comm_split(ncomm,commx1%rankg,commx1%rankl, &
                            commx1%comm,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split2D: MPI_Comm_split_1: ierr=", ierr
        CALL MPI_Comm_split(ncomm,commx2%rankg,commx2%rankl, &
                            commx2%comm,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split2D: MPI_Comm_split_2: ierr=", ierr

        CALL MPI_Comm_rank(commx1%comm,commx1%rank,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_rank_1: ierr=", ierr
        CALL MPI_Comm_rank(commx2%comm,commx2%rank,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_rank_2: ierr=", ierr
        CALL MPI_Comm_size(commx1%comm,commx1%size,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_size_1: ierr=", ierr
        CALL MPI_Comm_size(commx2%comm,commx2%size,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_size_2: ierr=", ierr
        RETURN
      END SUBROUTINE mtx_comm_split2D

!-----

      SUBROUTINE mtx_comm_split3D(n1,n2,n3,commx1,commx2,commx3,commx4,commx5,commx6)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: n1,n2,n3 ! number of groups
        TYPE(mtx_mpi_type),intent(OUT):: commx1,commx2,commx3,commx4,commx5,commx6
        integer:: ierr, na, nb

        IF(n1*n2*n3 > nsize) THEN
           WRITE(6,*) 'XX mtx_comm_split3D: n1*n2*n3 > nsize: ',n1,n2,n3,nsize
           STOP
        ENDIF
!       enable to communicate for each NSA 
        commx1%sizeg=n2*n3            ! number of groups
        commx1%sizel=n1               ! number of processors in a group
        commx1%rankg=MOD(NRANK,n2*n3) ! colors 
        commx1%rankl=nrank/(n2*n3)    ! keys = RANK

!       enable to communicate for each NR
        commx2%sizeg=n1*n3
        commx2%sizel=n2
        na = NRANK/(n2*n3)
        nb = na * (n2*n3)
        commx2%rankg=mod(NRANK,n3) + nb ! colorr
        commx2%rankl= (NRANK-nb)/n3     ! keyr

!       enable to communicate for each NP
        commx3%sizeg=n1*n2
        commx3%sizel=n3
        commx3%rankg=NRANK/n3   ! colorp
        commx3%rankl=MOD(nrank,n3) ! keyp

!       enable to communicate for each NP and NR
        commx4%sizeg=n1
        commx4%sizel=n2*n3
        commx4%rankg=NRANK/(n2*n3)   ! colorrp
        commx4%rankl=MOD(nrank,n2*n3) ! keyrp

!       enable to communicate for each NP and NSA
        commx5%sizeg=n2
        commx5%sizel=n1*n3
        commx5%rankg=MOD(NRANK,n2*n3)/n3   ! colorsp
        commx5%rankl=MOD(nrank,n3) + NRANK/(n2*n3)*n3 ! keysp

!       enable to communicate for each NR and NSA
        commx6%sizeg=n3
        commx6%sizel=n1*n2
        commx6%rankg=MOD(NRANK,n3)   ! colorsp
        commx6%rankl=NRANK/n3 ! keysp
!
        CALL MPI_Comm_split(ncomm,commx1%rankg,commx1%rankl, &
                            commx1%comm,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split3D: MPI_Comm_split_1: ierr=", ierr
!
        CALL MPI_Comm_split(ncomm,commx2%rankg,commx2%rankl, &
                            commx2%comm,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split3D: MPI_Comm_split_2: ierr=", ierr
!
        CALL MPI_Comm_split(ncomm,commx3%rankg,commx3%rankl, &
                            commx3%comm,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split3D: MPI_Comm_split_3: ierr=", ierr
!
        CALL MPI_Comm_split(ncomm,commx4%rankg,commx4%rankl, &
                            commx4%comm,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split3D: MPI_Comm_split_4: ierr=", ierr
!
        CALL MPI_Comm_split(ncomm,commx5%rankg,commx5%rankl, &
                            commx5%comm,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split3D: MPI_Comm_split_5: ierr=", ierr
!
        CALL MPI_Comm_split(ncomm,commx6%rankg,commx6%rankl, &
                            commx6%comm,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split3D: MPI_Comm_split_6: ierr=", ierr

!!!!
        CALL MPI_Comm_rank(commx1%comm,commx1%rank,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_rank_1: ierr=", ierr
        CALL MPI_Comm_rank(commx2%comm,commx2%rank,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_rank_2: ierr=", ierr
        CALL MPI_Comm_rank(commx3%comm,commx3%rank,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_rank_3: ierr=", ierr
        CALL MPI_Comm_rank(commx4%comm,commx4%rank,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_rank_4: ierr=", ierr
        CALL MPI_Comm_rank(commx5%comm,commx5%rank,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_rank_5: ierr=", ierr
        CALL MPI_Comm_rank(commx6%comm,commx6%rank,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_rank_6: ierr=", ierr
!
        CALL MPI_Comm_size(commx1%comm,commx1%size,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_size_1: ierr=", ierr
        CALL MPI_Comm_size(commx2%comm,commx2%size,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_size_2: ierr=", ierr
        CALL MPI_Comm_size(commx3%comm,commx3%size,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_size_3: ierr=", ierr
        CALL MPI_Comm_size(commx4%comm,commx4%size,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_size_4: ierr=", ierr
        CALL MPI_Comm_size(commx5%comm,commx5%size,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_size_5: ierr=", ierr
        CALL MPI_Comm_size(commx6%comm,commx6%size,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_split: MPI_Comm_size_6: ierr=", ierr

        RETURN
      END SUBROUTINE mtx_comm_split3D
!-----

      SUBROUTINE mtx_comm_free(commx)

        IMPLICIT NONE
        TYPE(mtx_mpi_type),intent(INOUT):: commx
        integer:: ierr

        CALL MPI_COMM_FREE(commx%comm,ierr)
        IF(ierr.ne.0) WRITE(6,*) &
             "XX mtx_comm_free: MPI_Comm_free: ierr=", ierr
        commx%comm=0
        RETURN
      END SUBROUTINE mtx_comm_free

!-----

      SUBROUTINE mtx_barrier
        IMPLICIT NONE
        INTEGER:: ierr

        call MPI_Barrier(ncomm,ierr)
        IF(ierr.NE.0) WRITE(6,*) &
             'XX mtx_barrier: MPI_Barrier: ierr=',ierr
        RETURN
      END SUBROUTINE mtx_barrier

!-----

      SUBROUTINE mtx_broadcast1_character(v)
        IMPLICIT NONE
        CHARACTER(LEN=1),INTENT(INOUT):: v

        CALL mtx_broadcast_character(v,1)

        RETURN
      END SUBROUTINE mtx_broadcast1_character

!-----

      SUBROUTINE mtx_broadcast1_logical(v)
        IMPLICIT NONE
        LOGICAL,INTENT(INOUT):: v
        LOGICAL,DIMENSION(1):: vdata

        vdata(1)=v
        CALL mtx_broadcast_logical(vdata,1)
        v=vdata(1)
        RETURN
      END SUBROUTINE mtx_broadcast1_logical

!-----

      SUBROUTINE mtx_broadcast1_integer(v)
        IMPLICIT NONE
        INTEGER,INTENT(INOUT):: v
        INTEGER,DIMENSION(1):: vdata

        vdata(1)=v
        CALL mtx_broadcast_integer(vdata,1)
        v=vdata(1)
        RETURN
      END SUBROUTINE mtx_broadcast1_integer

!-----

      SUBROUTINE mtx_broadcast1_real4(v)
        IMPLICIT NONE
        REAL(4),INTENT(INOUT):: v
        REAL(4),DIMENSION(1):: vdata

        vdata(1)=v
        CALL mtx_broadcast_real4(vdata,1)
        v=vdata(1)
        RETURN
      END SUBROUTINE mtx_broadcast1_real4

!-----

      SUBROUTINE mtx_broadcast1_real8(v)
        IMPLICIT NONE
        REAL(8),INTENT(INOUT):: v
        REAL(8),DIMENSION(1):: vdata

        vdata(1)=v
        CALL mtx_broadcast_real8(vdata,1)
        v=vdata(1)
        RETURN
      END SUBROUTINE mtx_broadcast1_real8

!-----

      SUBROUTINE mtx_broadcast1_complex8(v)
        IMPLICIT NONE
        COMPLEX(8),INTENT(INOUT):: v
        COMPLEX(8),DIMENSION(1):: vdata

        vdata(1)=v
        CALL mtx_broadcast_complex8(vdata,1)
        v=vdata(1)
        RETURN
      END SUBROUTINE mtx_broadcast1_complex8

!-----

      SUBROUTINE mtx_broadcast_character(vdata,ndata)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: ndata
        CHARACTER(LEN=ndata),INTENT(INOUT):: vdata
        INTEGER:: ierr
      
        call MPI_BCAST(vdata,ndata,MPI_CHARACTER,0,ncomm,ierr)
        IF(ierr.NE.0) WRITE(6,*) &
             'XX mtx_broadcast_character: MPI_BCAST: ierr=',ierr
        RETURN
      END SUBROUTINE mtx_broadcast_character

!-----

      SUBROUTINE mtx_broadcast_logical(vdata,ndata)
        IMPLICIT NONE
        LOGICAL,DIMENSION(ndata),INTENT(INOUT):: vdata
        INTEGER,INTENT(IN):: ndata
        INTEGER:: ierr
      
        call MPI_BCAST(vdata,ndata,MPI_LOGICAL,0,ncomm,ierr)
        IF(ierr.NE.0) WRITE(6,*) &
             'XX mtx_broadcast_logical: MPI_BCAST: ierr=',ierr
        RETURN
      END SUBROUTINE mtx_broadcast_logical

!-----

      SUBROUTINE mtx_broadcast_integer(vdata,ndata)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: ndata
      INTEGER:: ierr

      call MPI_BCAST(vdata,ndata,MPI_INTEGER,0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_broadcast_integer: MPI_BCAST: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_broadcast_integer

!-----

      SUBROUTINE mtx_broadcast_real4(vdata,ndata)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: ndata
      INTEGER:: ierr
      
      call MPI_BCAST(vdata,ndata,MPI_REAL,0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_broadcast_real4: MPI_BCAST: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_broadcast_real4

!-----

      SUBROUTINE mtx_broadcast_real8(vdata,ndata)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: ndata
      INTEGER:: ierr
      
      call MPI_BCAST(vdata,ndata,MPI_DOUBLE_PRECISION,0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_broadcast_real8: MPI_BCAST: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_broadcast_real8

!-----

      SUBROUTINE mtx_broadcast_complex8(vdata,ndata)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: ndata
      INTEGER:: ierr
      
      call MPI_BCAST(vdata,ndata,MPI_DOUBLE_COMPLEX,0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_broadcast_complex8: MPI_BCAST: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_broadcast_complex8

!-----

      SUBROUTINE mtx_broadcast2d_integer(vdata,n1,m1,m2)
      IMPLICIT NONE
      INTEGER,DIMENSION(n1,m2),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n1,m1,m2
!      INTEGER,DIMENSION(m1*m2):: tdata
      INTEGER,DIMENSION(:),allocatable:: tdata
      INTEGER:: i,i1,i2,ierr

      ALLOCATE(tdata(m1*m2))
      i=0
      DO i2=1,m2
         DO i1=1,m1
            i=i+1
            tdata(i)=vdata(i1,i2)
         END DO
      END DO
      CALL MPI_BCAST(tdata,i,MPI_INTEGER,0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_broadcast2d_integer: MPI_BCAST: ierr=',ierr
      i=0
      DO i2=1,m2
         DO i1=1,m1
            i=i+1
            vdata(i1,i2)=tdata(i)
         END DO
      END DO
      DEALLOCATE(tdata)
      RETURN
      END SUBROUTINE mtx_broadcast2D_integer

!-----

      SUBROUTINE mtx_broadcast2D_real4(vdata,n1,m1,m2)
      IMPLICIT NONE
      REAL(4),DIMENSION(n1,m2),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n1,m1,m2
!      REAL(4),DIMENSION(m1*m2):: tdata
      REAL(4),DIMENSION(:),allocatable:: tdata
      INTEGER:: i,i1,i2,ierr

      ALLOCATE(tdata(m1*m2))
      i=0
      DO i2=1,m2
         DO i1=1,m1
            i=i+1
            tdata(i)=vdata(i1,i2)
         END DO
      END DO
      call MPI_BCAST(tdata,i,MPI_REAL,0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_broadcast2D_real4: MPI_BCAST: ierr=',ierr
      i=0
      DO i2=1,m2
         DO i1=1,m1
            i=i+1
            vdata(i1,i2)=tdata(i)
         END DO
      END DO
      DEALLOCATE(tdata)
      RETURN
      END SUBROUTINE mtx_broadcast2D_real4

!-----

      SUBROUTINE mtx_broadcast2D_real8(vdata,n1,m1,m2)
      IMPLICIT NONE
      REAL(8),DIMENSION(n1,m2),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n1,m1,m2
!      REAL(8),DIMENSION(m1*m2):: tdata
      REAL(8),DIMENSION(:),allocatable:: tdata
      INTEGER:: i,i1,i2,ierr

      ALLOCATE(tdata(m1*m2))
      i=0
      DO i2=1,m2
         DO i1=1,m1
            i=i+1
            tdata(i)=vdata(i1,i2)
         END DO
      END DO
      call MPI_BCAST(tdata,i,MPI_DOUBLE_PRECISION,0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_broadcast2D_real8: MPI_BCAST: ierr=',ierr
      i=0
      DO i2=1,m2
         DO i1=1,m1
            i=i+1
            vdata(i1,i2)=tdata(i)
         END DO
      END DO
      DEALLOCATE(tdata)
      RETURN
      END SUBROUTINE mtx_broadcast2D_real8

!-----

      SUBROUTINE mtx_broadcast2D_complex8(vdata,n1,m1,m2)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(n1,m2),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n1,m1,m2
!      COMPLEX(8),DIMENSION(m1*m2):: tdata
      COMPLEX(8),DIMENSION(:),allocatable:: tdata
      INTEGER:: i,i1,i2,ierr

      ALLOCATE(tdata(m1*m2))
      i=0
      DO i2=1,m2
         DO i1=1,m1
            i=i+1
            tdata(i)=vdata(i1,i2)
         END DO
      END DO
      CALL MPI_BCAST(tdata,i,MPI_DOUBLE_COMPLEX,0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_broadcast2D_complex8: MPI_BCAST: ierr=',ierr
      i=0
      DO i2=1,m2
         DO i1=1,m1
            i=i+1
            vdata(i1,i2)=tdata(i)
         END DO
      END DO
      DEALLOCATE(tdata)
      RETURN
      END SUBROUTINE mtx_broadcast2D_complex8

!-----

      SUBROUTINE mtx_gather1_integer(vdata,vtot)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: vdata
      INTEGER,DIMENSION(nsize),INTENT(OUT):: vtot
      INTEGER,DIMENSION(1):: tdata

      tdata(1)=vdata
      CALL mtx_gather_integer(tdata,1,vtot)
      RETURN
      END SUBROUTINE mtx_gather1_integer

!-----

      SUBROUTINE mtx_gather1_real4(vdata,vtot)
      IMPLICIT NONE
      REAL(4),INTENT(IN):: vdata
      REAL(4),DIMENSION(nsize),INTENT(OUT):: vtot
      REAL(4),DIMENSION(1):: tdata

      tdata(1)=vdata
      CALL mtx_gather_real4(tdata,1,vtot)
      RETURN
      END SUBROUTINE mtx_gather1_real4

!-----

      SUBROUTINE mtx_gather1_real8(vdata,vtot)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: vdata
      REAL(8),DIMENSION(nsize),INTENT(OUT):: vtot
      REAL(8),DIMENSION(1):: tdata

      tdata(1)=vdata
      CALL mtx_gather_real8(tdata,1,vtot)
      RETURN
      END SUBROUTINE mtx_gather1_real8

!-----

      SUBROUTINE mtx_gather1_complex8(vdata,vtot)
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN):: vdata
      COMPLEX(8),DIMENSION(nsize),INTENT(OUT):: vtot
      COMPLEX(8),DIMENSION(1):: tdata

      tdata(1)=vdata
      CALL mtx_gather_complex8(tdata,1,vtot)
      RETURN
      END SUBROUTINE mtx_gather1_complex8

!-----

      SUBROUTINE mtx_gather_integer(vdata,ndata,vtot)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      INTEGER,DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: ierr

      call MPI_GATHER(vdata,ndata,MPI_INTEGER, &
                      vtot,ndata,MPI_INTEGER, &
                      0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_integer: MPI_GATHER: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_gather_integer

!-----

      SUBROUTINE mtx_gather_real4(vdata,ndata,vtot)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(4),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: ierr

      call MPI_GATHER(vdata,ndata,MPI_REAL, &
                      vtot,ndata,MPI_REAL, &
                      0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_real4: MPI_GATHER: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_gather_real4

!-----

      SUBROUTINE mtx_gather_real8(vdata,ndata,vtot)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(8),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: ierr

      call MPI_GATHER(vdata,ndata,MPI_DOUBLE_PRECISION, &
                      vtot,ndata,MPI_DOUBLE_PRECISION, &
                      0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_real8: MPI_GATHER: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_gather_real8

!-----

      SUBROUTINE mtx_gather_complex8(vdata,ndata,vtot)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      COMPLEX(8),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: ierr

      call MPI_GATHER(vdata,ndata,MPI_DOUBLE_COMPLEX, &
                      vtot,ndata,MPI_DOUBLE_COMPLEX, &
                      0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_complex8: MPI_GATHER: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_gather_complex8

!-----

      SUBROUTINE mtx_gatherv_integer(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      INTEGER,DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize):: ilena,iposa
      INTEGER:: ierr

      CALL MPI_GATHERV(vdata,ndata,MPI_INTEGER, &
                       vtot,ilena,iposa,MPI_INTEGER, &
                       0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gatherv_integer: MPI_GATHERV: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_gatherv_integer

!-----

      SUBROUTINE mtx_gatherv_real4(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(4),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize):: ilena,iposa
      INTEGER:: ierr

      CALL MPI_GATHERV(vdata,ndata,MPI_REAL, &
                       vtot,ilena,iposa,MPI_REAL, &
                       0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gatherv_real4: MPI_GATHERV: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_gatherv_real4

!-----

      SUBROUTINE mtx_gatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize):: ilena,iposa
      INTEGER:: ierr

      CALL MPI_GATHERV(vdata,ndata,MPI_DOUBLE_PRECISION, &
                       vtot,ilena,iposa,MPI_DOUBLE_PRECISION, &
                       0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gatherv_real8: MPI_GATHERV: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_gatherv_real8

!-----

      SUBROUTINE mtx_gatherv_complex8(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      COMPLEX(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize):: ilena,iposa
      INTEGER:: ierr

      CALL MPI_GATHERV(vdata,ndata,MPI_DOUBLE_COMPLEX, &
                       vtot,ilena,iposa,MPI_DOUBLE_COMPLEX, &
                       0,ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gatherv_complex8: MPI_GATHERV: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_gatherv_complex8

!-----

      SUBROUTINE mtx_allgather1_integer(vdata,vtot)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: vdata
      INTEGER,DIMENSION(nsize),INTENT(OUT):: vtot
      INTEGER,DIMENSION(1):: tdata

      tdata(1)=vdata
      CALL mtx_allgather_integer(tdata,1,vtot)
      RETURN
      END SUBROUTINE mtx_allgather1_integer

!-----

      SUBROUTINE mtx_allgather1_real4(vdata,vtot)
      IMPLICIT NONE
      REAL(4),INTENT(IN):: vdata
      REAL(4),DIMENSION(nsize),INTENT(OUT):: vtot
      REAL(4),DIMENSION(1):: tdata

      tdata(1)=vdata
      CALL mtx_allgather_real4(tdata,1,vtot)
      RETURN
      END SUBROUTINE mtx_allgather1_real4

!-----

      SUBROUTINE mtx_allgather1_real8(vdata,vtot)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: vdata
      REAL(8),DIMENSION(nsize),INTENT(OUT):: vtot
      REAL(8),DIMENSION(1):: tdata

      tdata(1)=vdata
      CALL mtx_allgather_real8(tdata,1,vtot)
      RETURN
      END SUBROUTINE mtx_allgather1_real8

!-----

      SUBROUTINE mtx_allgather1_complex8(vdata,vtot)
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN):: vdata
      COMPLEX(8),DIMENSION(nsize),INTENT(OUT):: vtot
      COMPLEX(8),DIMENSION(1):: tdata

      tdata(1)=vdata
      CALL mtx_allgather_complex8(tdata,1,vtot)
      RETURN
      END SUBROUTINE mtx_allgather1_complex8

!-----

      SUBROUTINE mtx_allgather_integer(vdata,ndata,vtot)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      INTEGER,DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: ierr

      call MPI_ALLGATHER(vdata,ndata,MPI_INTEGER, &
                         vtot,ndata,MPI_INTEGER, &
                         ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allgather_integer: MPI_ALLGATHER: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allgather_integer

!-----

      SUBROUTINE mtx_allgather_real4(vdata,ndata,vtot)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(4),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: ierr

      call MPI_ALLGATHER(vdata,ndata,MPI_REAL, &
                         vtot,ndata,MPI_REAL, &
                         ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allgather_real4: MPI_ALLGATHER: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allgather_real4

!-----

      SUBROUTINE mtx_allgather_real8(vdata,ndata,vtot)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(8),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: ierr

      call MPI_ALLGATHER(vdata,ndata,MPI_DOUBLE_PRECISION, &
                         vtot,ndata,MPI_DOUBLE_PRECISION, &
                         ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allgather_real8: MPI_ALLGATHER: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allgather_real8

!-----

      SUBROUTINE mtx_allgather_complex8(vdata,ndata,vtot)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      COMPLEX(8),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: ierr

      CALL MPI_ALLGATHER(vdata,ndata,MPI_DOUBLE_COMPLEX, &
                         vtot,ndata,MPI_DOUBLE_COMPLEX, &
                         ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allgather_complex8: MPI_ALLGATHER: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allgather_complex8

!-----

      SUBROUTINE mtx_allgatherv_integer(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      INTEGER,DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize):: ilena,iposa
      INTEGER:: ierr

      CALL MPI_ALLGATHERV(vdata,ndata,MPI_INTEGER, &
                          vtot,ilena,iposa,MPI_INTEGER, &
                          ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allgatherv_integer: MPI_ALLGATHERV: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allgatherv_integer

!-----

      SUBROUTINE mtx_allgatherv_real4(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(4),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize):: ilena,iposa
      INTEGER:: ierr

      CALL MPI_ALLGATHERV(vdata,ndata,MPI_REAL, &
                          vtot,ilena,iposa,MPI_REAL, &
                          ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allgatherv_real4: MPI_ALLGATHERV: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allgatherv_real4

!-----

      SUBROUTINE mtx_allgatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize):: ilena,iposa
      INTEGER:: ierr

      CALL MPI_ALLGATHERV(vdata,ndata,MPI_DOUBLE_PRECISION, &
                          vtot,ilena,iposa,MPI_DOUBLE_PRECISION, &
                          ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allgatherv_real8: MPI_ALLGATHERV: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allgatherv_real8

!-----

      SUBROUTINE mtx_allgatherv_complex8(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      COMPLEX(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize):: ilena,iposa
      INTEGER:: ierr

      CALL MPI_ALLGATHERV(vdata,ndata,MPI_DOUBLE_COMPLEX, &
                          vtot,ilena,iposa,MPI_DOUBLE_COMPLEX, &
                          ncomm,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allgatherv_complex8: MPI_ALLGATHERV: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allgatherv_complex8

!-----

      SUBROUTINE mtx_reduce1_integer(vdata_,nop,vreduce_,vloc_)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: vdata_
      INTEGER,INTENT(IN):: nop
      INTEGER,INTENT(OUT):: vreduce_
      INTEGER,INTENT(OUT):: vloc_
      INTEGER,DIMENSION(1):: vdata
      INTEGER,DIMENSION(1):: vreduce
      INTEGER,DIMENSION(1):: vloc

      vdata(1)=vdata_
      CALL mtx_reduce_integer(vdata,1,nop,vreduce,vloc)
      vreduce_=vreduce(1)
      vloc_=vloc(1)
      RETURN
      END SUBROUTINE mtx_reduce1_integer
      
!-----

      SUBROUTINE mtx_reduce1_real4(vdata_,nop,vreduce_,vloc_)
      IMPLICIT NONE
      REAL(4),INTENT(IN):: vdata_
      INTEGER,INTENT(IN):: nop
      REAL(4),INTENT(OUT):: vreduce_
      INTEGER,INTENT(OUT):: vloc_
      REAL(4),DIMENSION(1):: vdata
      REAL(4),DIMENSION(1):: vreduce
      INTEGER,DIMENSION(1):: vloc

      vdata(1)=vdata_
      CALL mtx_reduce_real4(vdata,1,nop,vreduce,vloc)
      vreduce_=vreduce(1)
      vloc_=vloc(1)
      RETURN
      END SUBROUTINE mtx_reduce1_real4
      
!-----

      SUBROUTINE mtx_reduce1_real8(vdata_,nop,vreduce_,vloc_)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: vdata_
      INTEGER,INTENT(IN):: nop
      REAL(8),INTENT(OUT):: vreduce_
      INTEGER,INTENT(OUT):: vloc_
      REAL(8),DIMENSION(1):: vdata
      REAL(8),DIMENSION(1):: vreduce
      INTEGER,DIMENSION(1):: vloc

      vdata(1)=vdata_
      CALL mtx_reduce_real8(vdata,1,nop,vreduce,vloc)
      vreduce_=vreduce(1)
      vloc_=vloc(1)
      RETURN
      END SUBROUTINE mtx_reduce1_real8
      
!-----

      SUBROUTINE mtx_reduce1_complex8(vdata_,nop,vreduce_,vloc_)
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN):: vdata_
      INTEGER,INTENT(IN):: nop
      COMPLEX(8),INTENT(OUT):: vreduce_
      INTEGER,INTENT(OUT):: vloc_
      COMPLEX(8),DIMENSION(1):: vdata
      COMPLEX(8),DIMENSION(1):: vreduce
      INTEGER,DIMENSION(1):: vloc

      vdata(1)=vdata_
      CALL mtx_reduce_complex8(vdata,1,nop,vreduce,vloc)
      vreduce_=vreduce(1)
      vloc_=vloc(1)
      RETURN
      END SUBROUTINE mtx_reduce1_complex8
      
!-----

      SUBROUTINE mtx_reduce_integer(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
!      INTEGER,DIMENSION(2,ndata):: d_send, d_recv
      INTEGER,DIMENSION(:,:),ALLOCATABLE:: d_send, d_recv
      INTEGER:: ierr,i

      DO i=1,ndata
         vloc(i)=0
      END DO
      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_INTEGER, &
                         MPI_MAX,0,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_INTEGER, &
                         MPI_MIN,0,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_INTEGER, &
                         MPI_SUM,0,ncomm,ierr)
      CASE(4,5)! MAX/MINOC
         ALLOCATE(d_send(2,ndata),d_recv(2,ndata))
         DO i=1,ndata
            d_send(1,i)=vdata(i)
            d_send(2,i)=nrank
         END DO
         SELECT CASE(NOP)
         CASE(4) ! MAXLOC
            CALL MPI_REDUCE(d_send,d_recv,ndata,MPI_2INTEGER, &
                            MPI_MAXLOC,0,ncomm,ierr)
         CASE(5) ! MINLOC
            CALL MPI_REDUCE(d_send,d_recv,ndata,MPI_2INTEGER, &
                            MPI_MINLOC,0,ncomm,ierr)
         END SELECT
         DO i=1,ndata
            vreduce(i)=d_recv(1,i)
            vloc(i)=d_recv(2,i)
         END DO
         DEALLOCATE(d_send,d_recv)
      END SELECT
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_reduce_integer: MPI_REDUCE: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_reduce_integer
      
!-----

      SUBROUTINE mtx_reduce_real4(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      REAL(4),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
!      REAL(4),DIMENSION(2,ndata):: d_send, d_recv
      REAL(4),DIMENSION(:,:),ALLOCATABLE:: d_send, d_recv
      INTEGER:: ierr, i

      DO i=1,ndata
         vloc(i)=0
      END DO
      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_REAL, &
                         MPI_MAX,0,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_REAL, &
                         MPI_MIN,0,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_REAL, &
                         MPI_SUM,0,ncomm,ierr)
      CASE(4,5)! MAX/MINLOC
         ALLOCATE(d_send(2,ndata),d_recv(2,ndata))
         DO i=1,ndata
            d_send(1,i)=vdata(i)
            d_send(2,i)=nrank*1.0
         END DO
         SELECT CASE(NOP)
         CASE(4) ! MAXLOC
            CALL MPI_REDUCE(d_send,d_recv,ndata,MPI_2REAL, &
                            MPI_MAXLOC,0,ncomm,ierr)
         CASE(5) ! MINLOC
            CALL MPI_REDUCE(d_send,d_recv,ndata,MPI_2REAL, &
                            MPI_MINLOC,0,ncomm,ierr)
         END SELECT
         DO i=1,ndata
            vreduce(i)=d_recv(1,i)
            vloc(i)=int(d_recv(2,i))
         END DO
         DEALLOCATE(d_send,d_recv)
      END SELECT
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_reduce_real4: MPI_REDUCE: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_reduce_real4
      
!-----

      SUBROUTINE mtx_reduce_real8(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      REAL(8),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
!      REAL(8),DIMENSION(2,ndata):: d_send, d_recv
      REAL(8),DIMENSION(:,:),ALLOCATABLE:: d_send, d_recv
      INTEGER:: ierr, i

      DO i=1,ndata
         vloc(i)=0
      END DO
      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_DOUBLE_PRECISION, &
                         MPI_MAX,0,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_DOUBLE_PRECISION, &
                         MPI_MIN,0,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,0,ncomm,ierr)
      CASE(4,5)! MAX/MINLOC
         ALLOCATE(d_send(2,ndata),d_recv(2,ndata))
         DO i=1,ndata
            d_send(1,i)=vdata(i)
            d_send(2,i)=nrank*1.D0
         END DO
         SELECT CASE(NOP)
         CASE(4) ! MAXLOC
            CALL MPI_REDUCE(d_send,d_recv,ndata,MPI_2DOUBLE_PRECISION, &
                            MPI_MAXLOC,0,ncomm,ierr)
         CASE(5) ! MINLOC
            CALL MPI_REDUCE(d_send,d_recv,ndata,MPI_2DOUBLE_PRECISION, &
                            MPI_MINLOC,0,ncomm,ierr)
         END SELECT
         DO i=1,ndata
            vreduce(i)=d_recv(1,i)
            vloc(i)=int(d_recv(2,i))
         END DO
         DEALLOCATE(d_send,d_recv)
      END SELECT
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_reduce_real8: MPI_REDUCE: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_reduce_real8
      
!-----

      SUBROUTINE mtx_reduce_complex8(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: ndata,nop
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      COMPLEX(8),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
!      COMPLEX(8),DIMENSION(2,ndata):: d_send, d_recv
      COMPLEX(8),DIMENSION(:,:),ALLOCATABLE:: d_send, d_recv
      INTEGER:: ierr,i

      DO i=1,ndata
         vloc(i)=0
      END DO
      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_DOUBLE_COMPLEX, &
                         MPI_MAX,0,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_DOUBLE_COMPLEX, &
                         MPI_MIN,0,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_DOUBLE_COMPLEX, &
                         MPI_SUM,0,ncomm,ierr)
      CASE(4,5)! MAX/MINLOC
         ALLOCATE(d_send(2,ndata),d_recv(2,ndata))
         DO i=1,ndata
            d_send(1,i)=abs(vdata(i))
            d_send(2,i)=nrank*1.D0
         END DO
         SELECT CASE(NOP)
         CASE(4) ! MAXLOC
            CALL MPI_REDUCE(d_send,d_recv,ndata,MPI_2DOUBLE_PRECISION, &
                            MPI_MAXLOC,0,ncomm,ierr)
         CASE(5) ! MINLOC
            CALL MPI_REDUCE(d_send,d_recv,ndata,MPI_2DOUBLE_PRECISION, &
                            MPI_MINLOC,0,ncomm,ierr)
         END SELECT
         DO i=1,ndata
            vreduce(i)=d_recv(1,i)
            vloc(i)=int(d_recv(2,i))
         END DO
         DEALLOCATE(d_send,d_recv)
      END SELECT
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_reduce_complex8: MPI_REDUCE: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_reduce_complex8
      
!-----

      SUBROUTINE mtx_allreduce1_integer(vdata_,nop,vreduce_,vloc_)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: vdata_
      INTEGER,INTENT(IN):: nop
      INTEGER,INTENT(OUT):: vreduce_
      INTEGER,INTENT(OUT):: vloc_
      INTEGER,DIMENSION(1):: vdata
      INTEGER,DIMENSION(1):: vreduce
      INTEGER,DIMENSION(1):: vloc

      vdata(1)=vdata_
      CALL mtx_allreduce_integer(vdata,1,nop,vreduce,vloc)
      vreduce_=vreduce(1)
      vloc_=vloc(1)
      RETURN
      END SUBROUTINE mtx_allreduce1_integer
      
!-----

      SUBROUTINE mtx_allreduce1_real4(vdata_,nop,vreduce_,vloc_)
      IMPLICIT NONE
      REAL(4),INTENT(IN):: vdata_
      INTEGER,INTENT(IN):: nop
      REAL(4),INTENT(OUT):: vreduce_
      INTEGER,INTENT(OUT):: vloc_
      REAL(4),DIMENSION(1):: vdata
      REAL(4),DIMENSION(1):: vreduce
      INTEGER,DIMENSION(1):: vloc

      vdata(1)=vdata_
      CALL mtx_allreduce_real4(vdata,1,nop,vreduce,vloc)
      vreduce_=vreduce(1)
      vloc_=vloc(1)
      RETURN
      END SUBROUTINE mtx_allreduce1_real4
      
!-----

      SUBROUTINE mtx_allreduce1_real8(vdata_,nop,vreduce_,vloc_)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: vdata_
      INTEGER,INTENT(IN):: nop
      REAL(8),INTENT(OUT):: vreduce_
      INTEGER,INTENT(OUT):: vloc_
      REAL(8),DIMENSION(1):: vdata
      REAL(8),DIMENSION(1):: vreduce
      INTEGER,DIMENSION(1):: vloc

      vdata(1)=vdata_
      CALL mtx_allreduce_real8(vdata,1,nop,vreduce,vloc)
      vreduce_=vreduce(1)
      vloc_=vloc(1)
      RETURN
      END SUBROUTINE mtx_allreduce1_real8
      
!-----

      SUBROUTINE mtx_allreduce1_complex8(vdata_,nop,vreduce_,vloc_)
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN):: vdata_
      INTEGER,INTENT(IN):: nop
      COMPLEX(8),INTENT(OUT):: vreduce_
      INTEGER,INTENT(OUT):: vloc_
      COMPLEX(8),DIMENSION(1):: vdata
      COMPLEX(8),DIMENSION(1):: vreduce
      INTEGER,DIMENSION(1):: vloc

      vdata(1)=vdata_
      CALL mtx_allreduce_complex8(vdata,1,nop,vreduce,vloc)
      vreduce_=vreduce(1)
      vloc_=vloc(1)
      RETURN
      END SUBROUTINE mtx_allreduce1_complex8
      
!-----

      SUBROUTINE mtx_allreduce_integer(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
!      INTEGER,DIMENSION(2,ndata):: d_send, d_recv
      INTEGER,DIMENSION(:,:),allocatable:: d_send, d_recv
      INTEGER:: ierr,i

      DO i=1,ndata
         vloc(i)=0
      END DO
      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_INTEGER, &
                            MPI_MAX,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_INTEGER, &
                            MPI_MIN,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_INTEGER, &
                            MPI_SUM,ncomm,ierr)
      CASE(4,5)! MAX/MINOC
         ALLOCATE(d_send(2,ndata),d_recv(2,ndata))
         DO i=1,ndata
            d_send(1,i)=vdata(i)
            d_send(2,i)=nrank
         END DO
         SELECT CASE(NOP)
         CASE(4) ! MAXLOC
            CALL MPI_ALLREDUCE(d_send,d_recv,ndata,MPI_2INTEGER, &
                               MPI_MAXLOC,ncomm,ierr)
         CASE(5) ! MINLOC
            CALL MPI_ALLREDUCE(d_send,d_recv,ndata,MPI_2INTEGER, &
                               MPI_MINLOC,ncomm,ierr)
         END SELECT
         DO i=1,ndata
            vreduce(i)=d_recv(1,i)
            vloc(i)=d_recv(2,i)
         END DO
         DEALLOCATE(d_send,d_recv)
      END SELECT
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allreduce_integer: MPI_allREDUCE: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allreduce_integer
      
!-----

      SUBROUTINE mtx_allreduce_real4(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      REAL(4),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
!      REAL(4),DIMENSION(2,ndata):: d_send, d_recv
      REAL(4),DIMENSION(:,:),allocatable:: d_send, d_recv
      INTEGER:: ierr,i

      DO i=1,ndata
         vloc(i)=0
      END DO
      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_REAL, &
                            MPI_MAX,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_REAL, &
                            MPI_MIN,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_REAL, &
                            MPI_SUM,ncomm,ierr)
      CASE(4,5)! MAX/MINLOC
         ALLOCATE(d_send(2,ndata),d_recv(2,ndata))
         DO i=1,ndata
            d_send(1,i)=vdata(i)
            d_send(2,i)=nrank*1.0
         END DO
         SELECT CASE(NOP)
         CASE(4) ! MAXLOC
            CALL MPI_ALLREDUCE(d_send,d_recv,ndata,MPI_2REAL, &
                               MPI_MAXLOC,ncomm,ierr)
         CASE(5) ! MINLOC
            CALL MPI_ALLREDUCE(d_send,d_recv,ndata,MPI_2REAL, &
                               MPI_MINLOC,ncomm,ierr)
         END SELECT
         DO i=1,ndata
            vreduce(i)=d_recv(1,i)
            vloc(i)=int(d_recv(2,i))
         END DO
         DEALLOCATE(d_send,d_recv)
      END SELECT
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allreduce_real4: MPI_ALLREDUCE: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allreduce_real4
      
!-----

      SUBROUTINE mtx_allreduce_real8(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      REAL(8),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
!      REAL(8),DIMENSION(2,ndata):: d_send, d_recv
      REAL(8),DIMENSION(:,:),allocatable:: d_send, d_recv
      INTEGER:: ierr,i

      DO i=1,ndata
         vloc(i)=0
      END DO
      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_DOUBLE_PRECISION, &
                            MPI_MAX,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_DOUBLE_PRECISION, &
                            MPI_MIN,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_DOUBLE_PRECISION, &
                            MPI_SUM,ncomm,ierr)
      CASE(4,5)! MAX/MINLOC
         allocate(d_send(2,ndata),d_recv(2,ndata))
         DO i=1,ndata
            d_send(1,i)=vdata(i)
            d_send(2,i)=nrank*1.D0
         END DO
         SELECT CASE(NOP)
         CASE(4) ! MAXLOC
            CALL MPI_ALLREDUCE(d_send,d_recv,ndata,MPI_2DOUBLE_PRECISION, &
                               MPI_MAXLOC,ncomm,ierr)
         CASE(5) ! MINLOC
            CALL MPI_ALLREDUCE(d_send,d_recv,ndata,MPI_2DOUBLE_PRECISION, &
                               MPI_MINLOC,ncomm,ierr)
         END SELECT
         DO i=1,ndata
            vreduce(i)=d_recv(1,i)
            vloc(i)=int(d_recv(2,i))
         END DO
         DEALLOCATE(d_send,d_recv)
      END SELECT
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allreduce_real8: MPI_ALLREDUCE: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allreduce_real8
      
!-----

      SUBROUTINE mtx_allreduce_complex8(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      COMPLEX(8),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
!      COMPLEX(8),DIMENSION(2,ndata):: d_send, d_recv
      COMPLEX(8),DIMENSION(:,:),allocatable:: d_send, d_recv
      INTEGER:: ierr,i

      DO i=1,ndata
         vloc(i)=0
      END DO
      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_DOUBLE_COMPLEX, &
                            MPI_MAX,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_DOUBLE_COMPLEX, &
                            MPI_MIN,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_DOUBLE_COMPLEX, &
                            MPI_SUM,ncomm,ierr)
      CASE(4,5)! MAX/MINLOC
         ALLOCATE(d_send(2,ndata),d_recv(2,ndata))
         DO i=1,ndata
            d_send(1,i)=abs(vdata(i))
            d_send(2,i)=nrank*1.D0
         END DO
         SELECT CASE(NOP)
         CASE(4) ! MAXLOC
            CALL MPI_ALLREDUCE(d_send,d_recv,ndata,MPI_DOUBLE_PRECISION, &
                               MPI_MAXLOC,ncomm,ierr)
         CASE(5) ! MINLOC
            CALL MPI_ALLREDUCE(d_send,d_recv,ndata,MPI_DOUBLE_PRECISION, &
                               MPI_MINLOC,ncomm,ierr)
         END SELECT
         DO i=1,ndata
            vreduce(i)=d_recv(1,i)
            vloc(i)=int(d_recv(2,i))
         END DO
         DEALLOCATE(d_send,d_recv)
      END SELECT
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allreduce_complex8: MPI_ALLREDUCE: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allreduce_complex8

!-----

      SUBROUTINE mtx_sendrecv_integer(sendbuf,sendcount,dest, &
                                    recvbuf,recvcount,source)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: sendcount,recvcount
        INTEGER,INTENT(IN),DIMENSION(sendcount):: sendbuf
        INTEGER,INTENT(OUT),DIMENSION(recvcount):: recvbuf
        INTEGER,INTENT(IN):: dest,source
        INTEGER:: dest_,source_,ierr,i
        integer:: istatus(MPI_STATUS_SIZE)

        IF(nsize.eq.1) then
           DO i=1,sendcount
              recvbuf(i)=sendbuf(i)
           END DO
        ELSE
           IF(dest.ge.nsize) THEN
              dest_=MPI_PROC_NULL
           ELSE
              dest_=dest
           END IF
           IF(source.lt.0) THEN
              source_=MPI_PROC_NULL
           ELSE
              source_=source
           END IF

           CALL MPI_SENDRECV( &
                sendbuf, sendcount, MPI_INTEGER, dest_,  0,&
                recvbuf, recvcount, MPI_INTEGER, source_,0,&
                ncomm, ISTATUS, IERR)
        END IF
        RETURN
      END SUBROUTINE mtx_sendrecv_integer

!-----

      SUBROUTINE mtx_sendrecv_real4(sendbuf,sendcount,dest, &
                                    recvbuf,recvcount,source)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: sendcount,recvcount
        REAL(4),INTENT(IN),DIMENSION(sendcount):: sendbuf
        REAL(4),INTENT(OUT),DIMENSION(recvcount):: recvbuf
        INTEGER,INTENT(IN):: dest,source
        INTEGER:: dest_,source_,ierr,i
        integer:: istatus(MPI_STATUS_SIZE)

        IF(nsize.eq.1) then
           DO i=1,sendcount
              recvbuf(i)=sendbuf(i)
           END DO
        ELSE
           IF(dest.ge.nsize) THEN
              dest_=MPI_PROC_NULL
           ELSEIF(dest.lt.0) THEN
              dest_=MPI_PROC_NULL
           ELSE
              dest_=dest
           END IF
           IF(source.lt.0) THEN
              source_=MPI_PROC_NULL
           ELSEIF(source.ge.nsize) THEN
              source_=MPI_PROC_NULL
           ELSE
              source_=source
           END IF

           CALL MPI_SENDRECV( &
                sendbuf, sendcount, MPI_REAL, dest_,  0,&
                recvbuf, recvcount, MPI_REAL, source_,0,&
                ncomm, ISTATUS, IERR)
        END IF
        RETURN
      END SUBROUTINE mtx_sendrecv_real4

!-----

      SUBROUTINE mtx_sendrecv_real8(sendbuf,sendcount,dest, &
                                    recvbuf,recvcount,source)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: sendcount,recvcount
        REAL(8),INTENT(IN),DIMENSION(sendcount):: sendbuf
        REAL(8),INTENT(OUT),DIMENSION(recvcount):: recvbuf
        INTEGER,INTENT(IN):: dest,source
        INTEGER:: dest_,source_,ierr,i
        integer:: istatus(MPI_STATUS_SIZE)

        IF(nsize.eq.1) then
           DO i=1,sendcount
              recvbuf(i)=sendbuf(i)
           END DO
        ELSE
           IF(dest.ge.nsize) THEN
              dest_=MPI_PROC_NULL
           ELSEIF(dest.lt.0) THEN
              dest_=MPI_PROC_NULL
           ELSE
              dest_=dest
           END IF

           IF(source.lt.0) THEN
              source_=MPI_PROC_NULL
           ELSEIF(source.ge.nsize) THEN
              source_=MPI_PROC_NULL
           ELSE
              source_=source
           END IF

           CALL MPI_SENDRECV( &
                sendbuf, sendcount, MPI_DOUBLE_PRECISION, dest_,  0,&
                recvbuf, recvcount, MPI_DOUBLE_PRECISION, source_,0,&
                ncomm, ISTATUS, IERR)
        END IF
        RETURN
      END SUBROUTINE mtx_sendrecv_real8

!-----

      SUBROUTINE mtx_sendrecv_complex8(sendbuf,sendcount,dest, &
                                    recvbuf,recvcount,source)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: sendcount,recvcount
        COMPLEX(8),INTENT(IN),DIMENSION(sendcount):: sendbuf
        COMPLEX(8),INTENT(OUT),DIMENSION(recvcount):: recvbuf
        INTEGER,INTENT(IN):: dest,source
        INTEGER:: dest_,source_,ierr,i
        integer:: istatus(MPI_STATUS_SIZE)

        IF(nsize.eq.1) then
           DO i=1,sendcount
              recvbuf(i)=sendbuf(i)
           END DO
        ELSE
           IF(dest.ge.nsize) THEN
              dest_=MPI_PROC_NULL
           ELSE
              dest_=dest
           END IF
           IF(source.lt.0) THEN
              source_=MPI_PROC_NULL
           ELSE
              source_=source
           END IF

           CALL MPI_SENDRECV( &
                sendbuf, sendcount, MPI_DOUBLE_COMPLEX, dest_,  0,&
                recvbuf, recvcount, MPI_DOUBLE_COMPLEX, source_,0,&
                ncomm, ISTATUS, IERR)
        END IF
        RETURN
      END SUBROUTINE mtx_sendrecv_complex8
!-----------------------------------------------------------
      SUBROUTINE mtx_scatterv_real8(sendbuf,sendbuf_size,recvcount,recvbuf)

      IMPLICIT NONE
      integer:: ierr
      integer,intent(in):: sendbuf_size, recvcount
      double precision,dimension(sendbuf_size),intent(in)::sendbuf
      double precision,dimension(recvcount),intent(out)::recvbuf

      call MPI_Scatter(sendbuf, recvcount, MPI_DOUBLE_PRECISION, &
                       recvbuf, recvcount, MPI_DOUBLE_PRECISION, &
                       0, ncomm, ierr)

      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_scatterv_real8: MPI_SCATTER: ierr=',ierr
      RETURN

      END SUBROUTINE mtx_scatterv_real8
!-----------------------------------------------------------
      SUBROUTINE mtx_send_real8(sendbuf,sendcount,dest,tag)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: sendcount
      REAL(8),INTENT(IN),DIMENSION(sendcount):: sendbuf
      INTEGER,INTENT(IN):: dest,tag
      INTEGER:: ierr

      CALL MPI_SEND(sendbuf, sendcount, MPI_DOUBLE_PRECISION, &
           dest, tag, ncomm, ierr)


      END SUBROUTINE mtx_send_real8
!-----------------------------------------------------------
      SUBROUTINE mtx_recv_real8(recvbuf,recvcount,source,tag)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: recvcount
      REAL(8),INTENT(OUT),DIMENSION(recvcount):: recvbuf
      INTEGER,INTENT(IN):: source,tag
      INTEGER:: ierr
      integer:: istatus(MPI_STATUS_SIZE)

      CALL MPI_RECV(recvbuf, recvcount, MPI_DOUBLE_PRECISION, &
           source, tag, ncomm, istatus, ierr)


      END SUBROUTINE MTX_RECV_REAL8
!-----------------------------------------------------------
      SUBROUTINE mtx_abort(ierr_g)

      IMPLICIT NONE
      integer:: ierr, ierr_g

      CALL MPI_ABORT(ncomm,ierr_g,ierr)

      END SUBROUTINE mtx_abort
!-----
      END MODULE libmpi

