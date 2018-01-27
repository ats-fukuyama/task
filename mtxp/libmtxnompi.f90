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

      PRIVATE

      TYPE mtx_mpi_type
         integer:: comm   ! communicator
         integer:: rank   ! rank of processor
         integer:: size   ! number of processors
         integer:: rankg  ! rank of the group (color)
         integer:: sizeg  ! number of groups
         integer:: rankl  ! rank of processor in the group (key) = rank
         integer:: sizel  ! number of processors in the groups =size
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
      PUBLIC mtx_sendrecv_comple8

      PUBLIC mtx_scatterv_real8
      PUBLIC mtx_send_real8
      PUBLIC mtx_recv_real8

      PUBLIC mtx_abort

      TYPE(mtx_mpi_type):: mtx_global
      INTEGER:: ncomm,nrank,nsize

      INTEGER:: real8_buffer_count
      REAL(8),DIMENSION(:),ALLOCATABLE:: real8_buffer


      CONTAINS

!-----

      SUBROUTINE mtx_mpi(ncomm_,nrank_,nsize_)
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

        ncomm=ncomm_in
        nsize=1
        nrank=0
        mtx_global%comm=ncomm
        mtx_global%rank=0
        mtx_global%size=1
        mtx_global%rankg=0
        mtx_global%sizeg=1
        mtx_global%rankl=0
        mtx_global%sizel=1
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

        IF(n1*n2 > nsize) THEN
           WRITE(6,*) 'XX mtx_comm_split2D: n1*n2 > nsize: ',n1,n2,nsize
           STOP
        ENDIF

        commx1%sizeg=n2       ! number of groups
        commx1%sizel=n1       ! number of processors in the group
        commx1%rankg=MOD(nrank,n2) ! color1
        commx1%rankl=nrank/n2        ! key1

        commx2%sizeg=n1       ! number of groups
        commx2%sizel=n2       ! number of processors in the group
        commx2%rankg=nrank/n2 ! color2
        commx2%rankl=MOD(nrank,n2) ! key2

        commx1%comm=ncomm
        commx1%rank=nrank
        commx1%size=nsize
        commx2%comm=ncomm
        commx2%rank=nrank
        commx2%size=nsize
        RETURN
      END SUBROUTINE mtx_comm_split2D

!-----

      SUBROUTINE mtx_comm_split3D(n1,n2,n3,commx1,commx2,commx3, &
                                           commx23,commx31,commx12)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: n1,n2,n3 ! number of groups
        TYPE(mtx_mpi_type),intent(OUT):: commx1,commx2,commx3, &
                                         commx23,commx31,commx12
        integer:: na, nb

        IF(n1*n2*n3 > nsize) THEN
           WRITE(6,*) 'XX mtx_comm_split3D: n1*n2*n3 > nsize: ',n1,n2,n3,nsize
           STOP
        ENDIF
!       enable to communicate for each n1
        commx1%sizeg=n2*n3            ! number of groups
        commx1%sizel=n1               ! number of processors in a group
        commx1%rankg=MOD(NRANK,n2*n3) ! colors 
        commx1%rankl=nrank/(n2*n3)    ! keys = RANK

!       enable to communicate for each n2
        commx2%sizeg=n1*n3
        commx2%sizel=n2
        na = NRANK/(n2*n3)
        nb = na * (n2*n3)
        commx2%rankg=mod(NRANK,n3) + nb ! colorr
        commx2%rankl= (NRANK-nb)/n3     ! keyr

!       enable to communicate for each n3
        commx3%sizeg=n1*n2
        commx3%sizel=n3
        commx3%rankg=NRANK/n3   ! colorp
        commx3%rankl=MOD(nrank,n3) ! keyp

!       enable to communicate for each n2 and n3
        commx23%sizeg=n1
        commx23%sizel=n2*n3
        commx23%rankg=NRANK/(n2*n3)   ! colorrp
        commx23%rankl=MOD(nrank,n2*n3) ! keyrp

!       enable to communicate for each n3 and n1
        commx31%sizeg=n2
        commx31%sizel=n1*n3
        commx31%rankg=MOD(NRANK,n2*n3)/n3   ! colorsp
        commx31%rankl=MOD(nrank,n3) + NRANK/(n2*n3)*n3 ! keysp

!       enable to communicate for each n1 and n2
        commx12%sizeg=n3
        commx12%sizel=n1*n2
        commx12%rankg=MOD(NRANK,n3)   ! colorsp
        commx12%rankl=NRANK/n3 ! keysp

        commx1%comm=ncomm
        commx1%rank=nrank
        commx1%size=nsize
        commx2%comm=ncomm
        commx2%rank=nrank
        commx2%size=nsize
        commx3%comm=ncomm
        commx3%rank=nrank
        commx3%size=nsize

        commx23%comm=ncomm
        commx23%rank=nrank
        commx23%size=nsize
        commx31%comm=ncomm
        commx31%rank=nrank
        commx31%size=nsize
        commx12%comm=ncomm
        commx12%rank=nrank
        commx12%size=nsize
        RETURN
      END SUBROUTINE mtx_comm_split3D
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

      SUBROUTINE mtx_broadcast1_real4(v)
        IMPLICIT NONE
        REAL(4),INTENT(INOUT):: v
        RETURN
      END SUBROUTINE mtx_broadcast1_real4

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
        INTEGER,INTENT(IN):: ndata
        CHARACTER(LEN=ndata),INTENT(INOUT):: vdata
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

      SUBROUTINE mtx_broadcast_real4(vdata,ndata)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: ndata
      RETURN
      END SUBROUTINE mtx_broadcast_real4

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

      SUBROUTINE mtx_broadcast2D_real4(vdata,n1,m1,m2)
      IMPLICIT NONE
      REAL(4),DIMENSION(n1,m2),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n1,m1,m2
      RETURN
      END SUBROUTINE mtx_broadcast2D_real4

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

      SUBROUTINE mtx_gather1_real4(vdata,vtot)
      IMPLICIT NONE
      REAL(4),INTENT(IN):: vdata
      REAL(4),DIMENSION(nsize),INTENT(OUT):: vtot

      vtot(1)=vdata
      RETURN
      END SUBROUTINE mtx_gather1_real4

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

      SUBROUTINE mtx_gather_real4(vdata,ndata,vtot)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(4),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      RETURN
      END SUBROUTINE mtx_gather_real4

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

      SUBROUTINE mtx_gatherv_real4(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(4),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize),INTENT(OUT):: ilena,iposa
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      ilena(1)=ndata
      iposa(1)=1
      RETURN
      END SUBROUTINE mtx_gatherv_real4

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

      SUBROUTINE mtx_allgather1_real4(vdata,vtot)
      IMPLICIT NONE
      REAL(4),INTENT(IN):: vdata
      REAL(4),DIMENSION(nsize),INTENT(OUT):: vtot

      vtot(1)=vdata
      RETURN
      END SUBROUTINE mtx_allgather1_real4

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

      SUBROUTINE mtx_allgather_real4(vdata,ndata,vtot)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(4),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      RETURN
      END SUBROUTINE mtx_allgather_real4

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

      SUBROUTINE mtx_allgatherv_real4(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(4),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntot
      INTEGER,DIMENSION(nsize),INTENT(OUT):: ilena,iposa
      INTEGER:: i

      DO i=1,ndata
         vtot(i)=vdata(i)
      END DO
      ilena(1)=ndata
      iposa(1)=1
      RETURN
      END SUBROUTINE mtx_allgatherv_real4

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

      SUBROUTINE mtx_reduce1_real4(vdata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(4),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: nop
      REAL(4),INTENT(OUT):: vreduce
      INTEGER,INTENT(OUT):: vloc

      vreduce=vdata
      vloc=1
      RETURN
      END SUBROUTINE mtx_reduce1_real4
      
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

      SUBROUTINE mtx_reduce_real4(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      REAL(4),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
      INTEGER:: i

      DO i=1,ndata
         vreduce(i)=vdata(i)
         vloc(i)=1
      END DO
      RETURN
      END SUBROUTINE mtx_reduce_real4
      
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

      SUBROUTINE mtx_allreduce1_real4(vdata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(4),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: nop
      REAL(4),INTENT(OUT):: vreduce
      INTEGER,INTENT(OUT):: vloc
      INTEGER:: i

      vreduce=vdata
      vloc=1
      RETURN
      END SUBROUTINE mtx_allreduce1_real4
      
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

      SUBROUTINE mtx_allreduce_real4(vdata,ndata,nop,vreduce,vloc)
      IMPLICIT NONE
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata,nop
      REAL(4),DIMENSION(ndata),INTENT(OUT):: vreduce
      INTEGER,DIMENSION(ndata),INTENT(OUT):: vloc
      INTEGER:: i

      DO i=1,ndata
         vreduce(i)=vdata(i)
         vloc(i)=1
      END DO
      RETURN
      END SUBROUTINE mtx_allreduce_real4
      
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
      
!-----

      SUBROUTINE mtx_sendrecv_integer(sendbuf,sendcount,dest, &
                                      recvbuf,recvcount,source)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: sendcount,recvcount
        INTEGER,INTENT(IN),DIMENSION(sendcount):: sendbuf
        INTEGER,INTENT(OUT),DIMENSION(recvcount):: recvbuf
        INTEGER,INTENT(IN):: dest,source
        INTEGER:: i

        DO i=1,sendcount
           recvbuf(i)=sendbuf(i)
        END DO
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
        INTEGER:: i

        DO i=1,sendcount
           recvbuf(i)=sendbuf(i)
        END DO
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
        INTEGER:: i

        DO i=1,sendcount
           recvbuf(i)=sendbuf(i)
        END DO
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
        INTEGER:: i

        DO i=1,sendcount
           recvbuf(i)=sendbuf(i)
        END DO
        RETURN
      END SUBROUTINE mtx_sendrecv_complex8

!-----------------------------------------------------------
      SUBROUTINE mtx_scatterv_real8(sendbuf,sendbuf_size,recvcount,recvbuf)

      IMPLICIT NONE
      integer,intent(in):: sendbuf_size, recvcount
      double precision,dimension(sendbuf_size),intent(in)::sendbuf
      double precision,dimension(recvcount),intent(out)::recvbuf
      integer:: i

      DO i=1,recvcount
         recvbuf(i)=sendbuf(i)
      END DO
      RETURN

      END SUBROUTINE mtx_scatterv_real8
!-----------------------------------------------------------
      SUBROUTINE mtx_send_real8(sendbuf,sendcount,dest,tag)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: sendcount
      REAL(8),INTENT(IN),DIMENSION(sendcount):: sendbuf
      INTEGER,INTENT(IN):: dest,tag
      INTEGER:: i

      ALLOCATE(real8_buffer(sendcount))

      real8_buffer_count=sendcount
      DO i=1,sendcount
         real8_buffer(i)=sendbuf(i)
      END DO

      END SUBROUTINE mtx_send_real8
!-----------------------------------------------------------
      SUBROUTINE mtx_recv_real8(recvbuf,recvcount,source,tag)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: recvcount
      REAL(8),INTENT(OUT),DIMENSION(recvcount):: recvbuf
      INTEGER,INTENT(IN):: source,tag
      INTEGER:: i

      DO i=1,MIN(recvcount,real8_buffer_count)
         recvbuf(i)=real8_buffer(i)
      END DO
      DEALLOCATE(real8_buffer)

      END SUBROUTINE MTX_RECV_REAL8
!-----------------------------------------------------------
      SUBROUTINE mtx_abort(ierr_g)

      IMPLICIT NONE
      integer:: ierr_g

      ierr_g=0

      END SUBROUTINE mtx_abort
!-----
    END MODULE libmpi

