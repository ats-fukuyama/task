!     $Id$

      MODULE libmpi

      integer:: rank,nsize

      CONTAINS

!-----

      SUBROUTINE mtx_initialize(rank_,nsize_,ncom)
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: rank_,nsize_,ncom

      WRITE(6,'(A)') '# mtx_initialize: libmtxbnd'
      ncom=0
      rank=0
      nsize=1
      rank_=0
      nsize_=1
      RETURN
      END SUBROUTINE mtx_initialize

!-----

      SUBROUTINE mtx_finalize
      RETURN
      END SUBROUTINE mtx_finalize

!-----

      SUBROUTINE mtx_barrier
      RETURN
      END SUBROUTINE mtx_barrier

!-----

      SUBROUTINE mtx_broadcast_character(kdata,n)
      IMPLICIT NONE
      CHARACTER(LEN=n),INTENT(INOUT):: kdata
      INTEGER,INTENT(IN):: n      
      RETURN
      END SUBROUTINE mtx_broadcast_character

!-----

      SUBROUTINE mtx_broadcast_integer(idata,n)
      IMPLICIT NONE
      INTEGER,DIMENSION(n),INTENT(INOUT):: idata
      INTEGER,INTENT(IN):: n      
      RETURN
      END SUBROUTINE mtx_broadcast_integer

!-----

      SUBROUTINE mtx_broadcast_real8(vdata,n)
      IMPLICIT NONE
      REAL(8),DIMENSION(n),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n
            RETURN
      END SUBROUTINE mtx_broadcast_real8

!-----

      SUBROUTINE mtx_broadcast_complex8(vdata,n)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(n),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n
            RETURN
      END SUBROUTINE mtx_broadcast_complex8

!-----

      SUBROUTINE mtx_broadcast_integer_2d(idata,n1,m1,m2)
      IMPLICIT none
      INTEGER,DIMENSION(n1,m2),INTENT(INOUT):: idata
      INTEGER,INTENT(IN):: n1,m1,m2
      RETURN
      END SUBROUTINE mtx_broadcast_integer_2D

!-----

      SUBROUTINE mtx_broadcast_real8_2D(vdata,n1,m1,m2)
      IMPLICIT NONE
      REAL(8),DIMENSION(n1,m2),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n1,m1,m2
      RETURN
      END SUBROUTINE mtx_broadcast_real8_2D

!-----

      SUBROUTINE mtx_broadcast_complex8_2D(vdata,n1,m1,m2)
      IMPLICIT NONE
      COMPLEX(8),DIMENSION(n1,m2),INTENT(INOUT):: vdata
      INTEGER,INTENT(IN):: n1,m1,m2
      RETURN
      END SUBROUTINE mtx_broadcast_complex8_2D

!-----

      SUBROUTINE mtx_gather_integer(idata,itot)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: idata
      INTEGER,DIMENSION(1),INTENT(OUT):: itot

      itot(1)=idata
      RETURN
      END SUBROUTINE mtx_gather_integer

!-----

      SUBROUTINE mtx_gather_real8(ddata,dtot)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: ddata
      REAL(8),DIMENSION(1),INTENT(OUT):: dtot

      dtot(1)=ddata
      RETURN
      END SUBROUTINE mtx_gather_real8

!-----

      SUBROUTINE mtx_allgather_integer(idata,itot)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: idata
      INTEGER,DIMENSION(1),INTENT(OUT):: itot

      itot(1)=idata
      RETURN
      END SUBROUTINE mtx_allgather_integer

!-----

      SUBROUTINE mtx_gatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: ndata
      INTEGER,INTENT(INOUT):: ntot
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,DIMENSION(1):: ilena,iposa
      INTEGER:: n

      DO n=1,ndata
         vtot(n)=vdata(n)
      ENDDO
      RETURN
      END SUBROUTINE mtx_gatherv_real8

!-----

      SUBROUTINE mtx_allgatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa,ncom)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: ndata,ncom
      INTEGER,INTENT(INOUT):: ntot
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,DIMENSION(nsize),INTENT(IN):: ilena,iposa
      INTEGER:: n

      DO n=1,ndata
         vtot(n)=vdata(n)
      ENDDO
      RETURN
      END SUBROUTINE mtx_allgatherv_real8

! ----------

      SUBROUTINE mtx_allgatherv_complex8(vdata,ndata,vtot,ntot,ilena,iposa,ncom)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: ndata,ncom
      INTEGER,INTENT(INOUT):: ntot
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      COMPLEX(8),DIMENSION(ntot),INTENT(OUT):: vtot
      INTEGER,DIMENSION(nsize),INTENT(IN):: ilena,iposa
      INTEGER:: n

      DO n=1,ndata
         vtot(n)=vdata(n)
      ENDDO
      RETURN
      END SUBROUTINE mtx_allgatherv_complex8

!-----

      SUBROUTINE mtx_reduce_real8(din,NOP,dout,ncom)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: din
      INTEGER,INTENT(IN):: NOP,ncom
      REAL(8),INTENT(OUT):: dout

      dout=din
      RETURN
      END SUBROUTINE mtx_reduce_real8
      
!-----

      SUBROUTINE mtx_broadcast_integer1(idata)
      IMPLICIT NONE
      INTEGER,INTENT(INOUT):: idata
      
      RETURN
      END SUBROUTINE mtx_broadcast_integer1

!-----

      SUBROUTINE mtx_reduce_v_real8(din,ncount,NOP,dout)
      IMPLICIT NONE
      REAL(8),dimension(ncount),INTENT(IN):: din
      INTEGER,INTENT(IN):: NOP, ncount
      REAL(8),dimension(ncount),INTENT(OUT):: dout
      INTEGER:: n

      DO n=1,ncount
         dout(n)=din(n)
      END DO

      RETURN
      END SUBROUTINE mtx_reduce_v_real8

!-----

      SUBROUTINE mtx_maxloc_real8(din,ncount,dout,loc)
      IMPLICIT NONE
      REAL(8),dimension(ncount),INTENT(IN):: din
      REAL(8),dimension(ncount),INTENT(OUT):: dout
      INTEGER,INTENT(IN):: ncount
      INTEGER,dimension(ncount),INTENT(OUT):: loc
      INTEGER:: n

      DO n = 1, ncount
         dout(n) = din(n)
         loc(n) = 1
      END DO

      RETURN
      END SUBROUTINE mtx_maxloc_real8
      END MODULE libmpi
