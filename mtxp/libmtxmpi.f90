!     $Id$

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

      PUBLIC mtx_set_communicator_global
      PUBLIC mtx_set_communicator_local
      PUBLIC mtx_reset_communicator_local
      PUBLIC mtx_comm_split2D
      PUBLIC mtx_comm_free
      PUBLIC mtx_barrier
      PUBLIC mtx_broadcast1_character
      PUBLIC mtx_broadcast1_logical
      PUBLIC mtx_broadcast1_integer
      PUBLIC mtx_broadcast1_real8
      PUBLIC mtx_broadcast1_complex8
      PUBLIC mtx_broadcast_character
      PUBLIC mtx_broadcast_logical
      PUBLIC mtx_broadcast_integer
      PUBLIC mtx_broadcast_real8
      PUBLIC mtx_broadcast_complex8
      PUBLIC mtx_broadcast2D_integer
      PUBLIC mtx_broadcast2D_real8
      PUBLIC mtx_broadcast2D_complex8
      PUBLIC mtx_gather1_integer
      PUBLIC mtx_gather1_real8
      PUBLIC mtx_gather1_complex8
      PUBLIC mtx_gather_integer
      PUBLIC mtx_gather_real8
      PUBLIC mtx_gather_complex8
      PUBLIC mtx_gatherv_integer
      PUBLIC mtx_gatherv_real8
      PUBLIC mtx_gatherv_complex8
      PUBLIC mtx_allgather1_integer
      PUBLIC mtx_allgather1_real8
      PUBLIC mtx_allgather1_complex8
      PUBLIC mtx_allgather_integer
      PUBLIC mtx_allgather_real8
      PUBLIC mtx_allgather_complex8
      PUBLIC mtx_allgatherv_integer
      PUBLIC mtx_allgatherv_real8
      PUBLIC mtx_allgatherv_complex8
      PUBLIC mtx_reduce1_integer
      PUBLIC mtx_reduce1_real8
      PUBLIC mtx_reduce1_complex8
      PUBLIC mtx_reduce_integer
      PUBLIC mtx_reduce_real8
      PUBLIC mtx_reduce_complex8
      PUBLIC mtx_allreduce1_integer
      PUBLIC mtx_allreduce1_real8
      PUBLIC mtx_allreduce1_complex8
      PUBLIC mtx_allreduce_integer
      PUBLIC mtx_allreduce_real8
      PUBLIC mtx_allreduce_complex8

      TYPE(mtx_mpi_type):: mtx_global
      INTEGER:: ncomm,nrank,nsize

      CONTAINS

!-----

      SUBROUTINE mtx_set_communicator_global(ncomm_,nrank_,nsize_)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: ncomm_
        INTEGER,INTENT(OUT):: nrank_,nsize_
        INTEGER:: ierr

        ncomm=ncomm_
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
        nrank_=nrank
        nsize_=nsize
        return
      END SUBROUTINE mtx_set_communicator_global

!-----

      SUBROUTINE mtx_set_communicator_local(mtx_mpi)
        IMPLICIT NONE
        TYPE(mtx_mpi_type),INTENT(IN):: mtx_mpi

        ncomm=mtx_mpi%comm
        nrank=mtx_mpi%rank
        nsize=mtx_mpi%size
        return
      END SUBROUTINE mtx_set_communicator_local

!-----

      SUBROUTINE mtx_reset_communicator_local
        IMPLICIT NONE

        ncomm=mtx_global%comm
        nrank=mtx_global%rank
        nsize=mtx_global%size
        return
      END SUBROUTINE mtx_reset_communicator_local

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
        commx1%sizeg=n2
        commx2%sizeg=n1
        commx1%sizel=n1
        commx2%sizel=n2
        commx1%rankg=nrank/n1
        commx1%rankl=nrank - commx1%rankg*n1
        commx2%rankl=nrank/n1
        commx2%rankg=MOD(nrank,n1)

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

      SUBROUTINE mtx_broadcast1_real8(v)
        IMPLICIT NONE
        REAL(8),INTENT(INOUT):: v
        REAL(8),DIMENSION(1):: vdata
        INTEGER:: ierr

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
        INTEGER:: ierr

        vdata(1)=v
        CALL mtx_broadcast_complex8(vdata,1)
        v=vdata(1)
        RETURN
      END SUBROUTINE mtx_broadcast1_complex8

!-----

      SUBROUTINE mtx_broadcast_character(vdata,ndata)
        IMPLICIT NONE
        CHARACTER(LEN=ndata),INTENT(INOUT):: vdata
        INTEGER,INTENT(IN):: ndata
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
      INTEGER,DIMENSION(m1*m2):: tdata
      INTEGER:: i,i1,i2,ierr

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
      RETURN
      END SUBROUTINE mtx_broadcast2D_integer

!-----

      SUBROUTINE mtx_broadcast2D_real8(vdata,n1,m1,m2)
      IMPLICIT NONE
      REAL(8),DIMENSION(n1,m2),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n1,m1,m2
      REAL(8),DIMENSION(m1*m2):: tdata
      INTEGER:: i,i1,i2,ierr

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
      RETURN
      END SUBROUTINE mtx_broadcast2D_real8

!-----

      SUBROUTINE mtx_broadcast2D_complex8(vdata,n1,m1,m2)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(n1,m2),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n1,m1,m2
      COMPLEX(8),DIMENSION(m1*m2):: tdata
      INTEGER:: i,i1,i2,ierr

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
      RETURN
      END SUBROUTINE mtx_broadcast2D_complex8

!-----

      SUBROUTINE mtx_gather1_integer(vdata,vtot)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: vdata
      INTEGER,DIMENSION(nsize),INTENT(OUT):: vtot
      INTEGER,DIMENSION(1):: tdata
      INTEGER:: ierr

      tdata(1)=vdata
      CALL mtx_gather_integer(tdata,1,vtot)
      RETURN
      END SUBROUTINE mtx_gather1_integer

!-----

      SUBROUTINE mtx_gather1_real8(vdata,vtot)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: vdata
      REAL(8),DIMENSION(nsize),INTENT(OUT):: vtot
      REAL(8),DIMENSION(1):: tdata
      INTEGER:: ierr

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
      INTEGER:: ierr

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
      INTEGER:: ierr

      tdata(1)=vdata
      CALL mtx_allgather_integer(tdata,1,vtot)
      RETURN
      END SUBROUTINE mtx_allgather1_integer

!-----

      SUBROUTINE mtx_allgather1_real8(vdata,vtot)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: vdata
      REAL(8),DIMENSION(nsize),INTENT(OUT):: vtot
      REAL(8),DIMENSION(1):: tdata
      INTEGER:: ierr

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
      INTEGER:: ierr

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
      INTEGER:: ierr

      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_INTEGER, &
                         MPI_MAX,0,ncomm,ierr)
         CALL MPI_REDUCE(vdata,vloc,ndata,MPI_INTEGER, &
                         MPI_MAXLOC,0,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_INTEGER, &
                         MPI_MIN,0,ncomm,ierr)
         CALL MPI_REDUCE(vdata,vloc,ndata,MPI_INTEGER, &
                         MPI_MINLOC,0,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_INTEGER, &
                         MPI_SUM,0,ncomm,ierr)
      END SELECT
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_reduce_integer: MPI_REDUCE: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_reduce_integer
      
!-----

      SUBROUTINE mtx_reduce_real8(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      REAL(8),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
      INTEGER:: ierr

      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_DOUBLE_PRECISION, &
                         MPI_MAX,0,ncomm,ierr)
         CALL MPI_REDUCE(vdata,vloc,ndata,MPI_DOUBLE_PRECISION, &
                         MPI_MAXLOC,0,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_DOUBLE_PRECISION, &
                         MPI_MIN,0,ncomm,ierr)
         CALL MPI_REDUCE(vdata,vloc,ndata,MPI_DOUBLE_PRECISION, &
                         MPI_MINLOC,0,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,0,ncomm,ierr)
      END SELECT
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_reduce_real8: MPI_REDUCE: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_reduce_real8
      
!-----

      SUBROUTINE mtx_reduce_complex8(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      COMPLEX(8),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
      INTEGER:: ierr

      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_DOUBLE_COMPLEX, &
                         MPI_MAX,0,ncomm,ierr)
         CALL MPI_REDUCE(vdata,vloc,ndata,MPI_DOUBLE_COMPLEX, &
                         MPI_MAXLOC,0,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_DOUBLE_COMPLEX, &
                         MPI_MIN,0,ncomm,ierr)
         CALL MPI_REDUCE(vdata,vloc,ndata,MPI_DOUBLE_COMPLEX, &
                         MPI_MINLOC,0,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_REDUCE(vdata,vreduce,ndata,MPI_DOUBLE_COMPLEX, &
                         MPI_SUM,0,ncomm,ierr)
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
      INTEGER:: ierr

      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_INTEGER, &
                            MPI_MAX,ncomm,ierr)
         CALL MPI_ALLREDUCE(vdata,vloc,ndata,MPI_INTEGER, &
                            MPI_MAXLOC,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_INTEGER, &
                            MPI_MIN,ncomm,ierr)
         CALL MPI_ALLREDUCE(vdata,vloc,ndata,MPI_INTEGER, &
                            MPI_MINLOC,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_INTEGER, &
                            MPI_SUM,ncomm,ierr)
      END SELECT
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allreduce_integer: MPI_allREDUCE: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allreduce_integer
      
!-----

      SUBROUTINE mtx_allreduce_real8(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      REAL(8),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
      INTEGER:: ierr

      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_DOUBLE_PRECISION, &
                            MPI_MAX,ncomm,ierr)
         CALL MPI_ALLREDUCE(vdata,vloc,ndata,MPI_DOUBLE_PRECISION, &
                            MPI_MAXLOC,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_DOUBLE_PRECISION, &
                            MPI_MIN,ncomm,ierr)
         CALL MPI_ALLREDUCE(vdata,vloc,ndata,MPI_DOUBLE_PRECISION, &
                            MPI_MINLOC,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_DOUBLE_PRECISION, &
                            MPI_SUM,ncomm,ierr)
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
      INTEGER:: ierr

      SELECT CASE(NOP)
      CASE(1)! MAX
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_DOUBLE_COMPLEX, &
                            MPI_MAX,ncomm,ierr)
         CALL MPI_ALLREDUCE(vdata,vloc,ndata,MPI_DOUBLE_COMPLEX, &
                            MPI_MAXLOC,ncomm,ierr)
      CASE(2)! MIN
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_DOUBLE_COMPLEX, &
                            MPI_MIN,ncomm,ierr)
         CALL MPI_ALLREDUCE(vdata,vloc,ndata,MPI_DOUBLE_COMPLEX, &
                            MPI_MINLOC,ncomm,ierr)
      CASE(3)! SUM
         CALL MPI_ALLREDUCE(vdata,vreduce,ndata,MPI_DOUBLE_COMPLEX, &
                            MPI_SUM,ncomm,ierr)
      END SELECT
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allreduce_complex8: MPI_ALLREDUCE: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_allreduce_complex8
      
      END MODULE libmpi

