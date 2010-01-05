!     $Id$

      MODULE libmpi

      integer:: rank,size

      CONTAINS

!-----

      SUBROUTINE mtx_initialize(rank_,size_)
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: rank_,size_

      WRITE(6,'(A)') '# mtx_initialize: libmtxbnd'
      rank=0
      size=1
      rank_=0
      size_=1
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

      SUBROUTINE mtx_allgatherv_real8(vdata,ndata,vtot,ntot,ilena,iposa)
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
      END SUBROUTINE mtx_allgatherv_real8

      END MODULE libmpi
