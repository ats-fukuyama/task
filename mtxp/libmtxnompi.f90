!     $Id$

      MODULE libmpi

      PRIVATE

      TYPE mtx_mpi_type
         integer:: comm   ! communicator
         integer:: rank   ! rank of processor
         integer:: size   ! number of processors
         integer:: rankg  ! rank of the group (color)
         integer:: sizeg  ! number of groups
         integer:: rankl  ! rank of processor in the group (key)
         integer:: sizel  ! number of processors in the groups
      END TYPE mtx_mpi_type
      PUBLIC mtx_mpi_type

      PUBLIC mtx_set_communicator_global
      PUBLIC mtx_set_communicator_local
      PUBLIC mtx_reset_communicator_local
      PUBLIC mtx_comm_split
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

        ncomm=ncomm_
        nsize=1
        nrank=0
        mtx_global%comm=ncomm
        mtx_global%rank=0
        mtx_global%size=1
        mtx_global%rankg=0
        mtx_global%sizeg=1
        mtx_global%rankl=0
        mtx_global%sizel=1
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

      SUBROUTINE mtx_comm_split(ngroup,mtx_mpi)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: ngroup ! number of groups
        TYPE(mtx_mpi_type),intent(OUT):: mtx_mpi

        IF(ngroup > nsize) THEN
           WRITE(6,*) 'XX mtx_comm_split: ngroup > nsize: ',ngroup,nsize
           STOP
        ENDIF
        mtx_mpi%sizel=ngroup
        mtx_mpi%sizeg=nsize/ngroup
        mtx_mpi%rankg=nrank/mtx_mpi%sizeg
        mtx_mpi%rankl=nrank-mtx_mpi%rankg*mtx_mpi%sizeg
        mtx_mpi%comm=ncomm
        mtx_mpi%rank=nrank
        mtx_mpi%size=nsize
        RETURN
      END SUBROUTINE mtx_comm_split

!-----

      SUBROUTINE mtx_comm_free(mtx_mpi)
        IMPLICIT NONE
        TYPE(mtx_mpi_type),intent(INOUT):: mtx_mpi

        mtx_mpi%comm=0
        RETURN
      END SUBROUTINE mtx_comm_free

!-----

      SUBROUTINE mtx_barrier
        IMPLICIT NONE
        RETURN
      END SUBROUTINE mtx_barrier

!-----

      SUBROUTINE mtx_broadcast1_character(v)
        IMPLICIT NONE
        CHARACTER(LEN=1),INTENT(INOUT):: v
        RETURN
      END SUBROUTINE mtx_broadcast1_character

!-----

      SUBROUTINE mtx_broadcast1_logical(v)
        IMPLICIT NONE
        LOGICAL,INTENT(INOUT):: v
        RETURN
      END SUBROUTINE mtx_broadcast1_logical

!-----

      SUBROUTINE mtx_broadcast1_integer(v)
        IMPLICIT NONE
        INTEGER,INTENT(INOUT):: v
        RETURN
      END SUBROUTINE mtx_broadcast1_integer

!-----

      SUBROUTINE mtx_broadcast1_real8(v)
        IMPLICIT NONE
        REAL(8),INTENT(INOUT):: v
        RETURN
      END SUBROUTINE mtx_broadcast1_real8

!-----

      SUBROUTINE mtx_broadcast1_complex8(v)
        IMPLICIT NONE
        COMPLEX(8),INTENT(INOUT):: v
        RETURN
      END SUBROUTINE mtx_broadcast1_complex8

!-----

      SUBROUTINE mtx_broadcast_character(vdata,ndata)
        IMPLICIT NONE
        CHARACTER(LEN=ndata),INTENT(INOUT):: vdata
        INTEGER,INTENT(IN):: ndata
        RETURN
      END SUBROUTINE mtx_broadcast_character

!-----

      SUBROUTINE mtx_broadcast_logical(vdata,ndata)
        IMPLICIT NONE
        LOGICAL,DIMENSION(ndata),INTENT(INOUT):: vdata
        INTEGER,INTENT(IN):: ndata
        RETURN
      END SUBROUTINE mtx_broadcast_logical

!-----

      SUBROUTINE mtx_broadcast_integer(vdata,ndata)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: ndata
      RETURN
      END SUBROUTINE mtx_broadcast_integer

!-----

      SUBROUTINE mtx_broadcast_real8(vdata,ndata)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: ndata
      RETURN
      END SUBROUTINE mtx_broadcast_real8

!-----

      SUBROUTINE mtx_broadcast_complex8(vdata,ndata)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: ndata
      RETURN
      END SUBROUTINE mtx_broadcast_complex8

!-----

      SUBROUTINE mtx_broadcast2d_integer(vdata,n1,m1,m2)
      IMPLICIT NONE
      INTEGER,DIMENSION(n1,m2),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n1,m1,m2
      RETURN
      END SUBROUTINE mtx_broadcast2D_integer

!-----

      SUBROUTINE mtx_broadcast2D_real8(vdata,n1,m1,m2)
      IMPLICIT NONE
      REAL(8),DIMENSION(n1,m2),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n1,m1,m2
      RETURN
      END SUBROUTINE mtx_broadcast2D_real8

!-----

      SUBROUTINE mtx_broadcast2D_complex8(vdata,n1,m1,m2)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(n1,m2),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n1,m1,m2
      RETURN
      END SUBROUTINE mtx_broadcast2D_complex8

!-----

      SUBROUTINE mtx_gather1_integer(vdata,vtot)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: vdata
      INTEGER,DIMENSION(nsize),INTENT(OUT):: vtot

      vtot(1)=vdata
      RETURN
      END SUBROUTINE mtx_gather1_integer

!-----

      SUBROUTINE mtx_gather1_real8(vdata,vtot)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: vdata
      REAL(8),DIMENSION(nsize),INTENT(OUT):: vtot

      vtot(1)=vdata
      RETURN
      END SUBROUTINE mtx_gather1_real8

!-----

      SUBROUTINE mtx_gather1_complex8(vdata,vtot)
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN):: vdata
      COMPLEX(8),DIMENSION(nsize),INTENT(OUT):: vtot

      vtot(1)=vdata
      RETURN
      END SUBROUTINE mtx_gather1_complex8

!-----

      SUBROUTINE mtx_gather_integer(vdata,ndata,vtot)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      INTEGER,DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      RETURN
      END SUBROUTINE mtx_gather_integer

!-----

      SUBROUTINE mtx_gather_real8(vdata,ndata,vtot)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(8),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      RETURN
      END SUBROUTINE mtx_gather_real8

!-----

      SUBROUTINE mtx_gather_complex8(vdata,ndata,vtot)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      COMPLEX(8),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      RETURN
      END SUBROUTINE mtx_gather_complex8

!-----

      SUBROUTINE mtx_gatherv_integer(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      INTEGER,DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize),INTENT(OUT):: ilena,iposa
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      ilena(1)=ndata
      iposa(1)=1
      RETURN
      END SUBROUTINE mtx_gatherv_integer

!-----

      SUBROUTINE mtx_gatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize),INTENT(OUT):: ilena,iposa
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      ilena(1)=ndata
      iposa(1)=1
      RETURN
      END SUBROUTINE mtx_gatherv_real8

!-----

      SUBROUTINE mtx_gatherv_complex8(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      COMPLEX(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize),INTENT(OUT):: ilena,iposa
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      ilena(1)=ndata
      iposa(1)=1
      RETURN
      END SUBROUTINE mtx_gatherv_complex8

!-----

      SUBROUTINE mtx_allgather1_integer(vdata,vtot)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: vdata
      INTEGER,DIMENSION(nsize),INTENT(OUT):: vtot

      vtot(1)=vdata
      RETURN
      END SUBROUTINE mtx_allgather1_integer

!-----

      SUBROUTINE mtx_allgather1_real8(vdata,vtot)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: vdata
      REAL(8),DIMENSION(nsize),INTENT(OUT):: vtot

      vtot(1)=vdata
      RETURN
      END SUBROUTINE mtx_allgather1_real8

!-----

      SUBROUTINE mtx_allgather1_complex8(vdata,vtot)
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN):: vdata
      COMPLEX(8),DIMENSION(nsize),INTENT(OUT):: vtot

      vtot(1)=vdata
      RETURN
      END SUBROUTINE mtx_allgather1_complex8

!-----

      SUBROUTINE mtx_allgather_integer(vdata,ndata,vtot)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      INTEGER,DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      RETURN
      END SUBROUTINE mtx_allgather_integer

!-----

      SUBROUTINE mtx_allgather_real8(vdata,ndata,vtot)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(8),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      RETURN
      END SUBROUTINE mtx_allgather_real8

!-----

      SUBROUTINE mtx_allgather_complex8(vdata,ndata,vtot)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      COMPLEX(8),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      RETURN
      END SUBROUTINE mtx_allgather_complex8

!-----

      SUBROUTINE mtx_allgatherv_integer(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      INTEGER,DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize),INTENT(OUT):: ilena,iposa
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      ilena(1)=ndata
      iposa(1)=1
      RETURN
      END SUBROUTINE mtx_allgatherv_integer

!-----

      SUBROUTINE mtx_allgatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize),INTENT(OUT):: ilena,iposa
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      ilena(1)=ndata
      iposa(1)=1
      RETURN
      END SUBROUTINE mtx_allgatherv_real8

!-----

      SUBROUTINE mtx_allgatherv_complex8(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      COMPLEX(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize),INTENT(OUT):: ilena,iposa
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      ilena(1)=ndata
      iposa(1)=1
      RETURN
      END SUBROUTINE mtx_allgatherv_complex8

!-----

      SUBROUTINE mtx_reduce1_integer(vdata,nop,vreduce,vloc)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: vdata
      INTEGER,INTENT(IN):: nop
      INTEGER,INTENT(OUT):: vreduce
      INTEGER,INTENT(OUT):: vloc

      vreduce=vdata
      vloc=1
      RETURN
      END SUBROUTINE mtx_reduce1_integer
      
!-----

      SUBROUTINE mtx_reduce1_real8(vdata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: nop
      REAL(8),INTENT(OUT):: vreduce
      INTEGER,INTENT(OUT):: vloc

      vreduce=vdata
      vloc=1
      RETURN
      END SUBROUTINE mtx_reduce1_real8
      
!-----

      SUBROUTINE mtx_reduce1_complex8(vdata,nop,vreduce,vloc)
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: nop
      COMPLEX(8),INTENT(OUT):: vreduce
      INTEGER,INTENT(OUT):: vloc

      vreduce=vdata
      vloc=1
      RETURN
      END SUBROUTINE mtx_reduce1_complex8
      
!-----

      SUBROUTINE mtx_reduce_integer(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
      INTEGER:: i

      DO i=1,ndata
         vreduce(i)=vdata(i)
         vloc(i)=1
      END DO
      RETURN
      END SUBROUTINE mtx_reduce_integer
      
!-----

      SUBROUTINE mtx_reduce_real8(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      REAL(8),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
      INTEGER:: i

      DO i=1,ndata
         vreduce(i)=vdata(i)
         vloc(i)=1
      END DO
      RETURN
      END SUBROUTINE mtx_reduce_real8
      
!-----

      SUBROUTINE mtx_reduce_complex8(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      COMPLEX(8),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
      INTEGER:: i

      DO i=1,ndata
         vreduce(i)=vdata(i)
         vloc(i)=1
      END DO
      RETURN
      END SUBROUTINE mtx_reduce_complex8
      
!-----

      SUBROUTINE mtx_allreduce1_integer(vdata,nop,vreduce,vloc)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: vdata
      INTEGER,INTENT(IN):: nop
      INTEGER,INTENT(OUT):: vreduce
      INTEGER,INTENT(OUT):: vloc

      vreduce=vdata
      vloc=1
      RETURN
      END SUBROUTINE mtx_allreduce1_integer
      
!-----

      SUBROUTINE mtx_allreduce1_real8(vdata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: nop
      REAL(8),INTENT(OUT):: vreduce
      INTEGER,INTENT(OUT):: vloc
      INTEGER:: i

      vreduce=vdata
      vloc=1
      RETURN
      END SUBROUTINE mtx_allreduce1_real8
      
!-----

      SUBROUTINE mtx_allreduce1_complex8(vdata,nop,vreduce,vloc)
      IMPLICIT NONE
      COMPLEX(8),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: nop
      COMPLEX(8),INTENT(OUT):: vreduce
      INTEGER,INTENT(OUT):: vloc

      vreduce=vdata
      vloc=1
      RETURN
      END SUBROUTINE mtx_allreduce1_complex8
      
!-----

      SUBROUTINE mtx_allreduce_integer(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      INTEGER,DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
      INTEGER:: i

      DO i=1,ndata
         vreduce(i)=vdata(i)
         vloc(i)=1
      END DO
      RETURN
      END SUBROUTINE mtx_allreduce_integer
      
!-----

      SUBROUTINE mtx_allreduce_real8(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      REAL(8),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
      INTEGER:: i

      DO i=1,ndata
         vreduce(i)=vdata(i)
         vloc(i)=1
      END DO
      RETURN
      END SUBROUTINE mtx_allreduce_real8
      
!-----

      SUBROUTINE mtx_allreduce_complex8(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      COMPLEX(8),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
      INTEGER:: i

      DO i=1,ndata
         vreduce(i)=vdata(i)
         vloc(i)=1
      END DO
      RETURN
      END SUBROUTINE mtx_allreduce_complex8
      
      END MODULE libmpi

