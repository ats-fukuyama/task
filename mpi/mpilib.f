C     $Id$
C
      SUBROUTINE MPSYNC
C
      call mtx_barrier
C
      RETURN
      END
C
      SUBROUTINE MPBCDN(vdata,ndata)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: ndata
      REAL(8),DIMENSION(ndata),INTENT(INOUT):: vdata
C
      call mtx_broadcast_real8(vdata,ndata)
C
      RETURN
      END
C
      SUBROUTINE MPBCRN(vdata,ndata)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: ndata
      REAL(4),DIMENSION(ndata),INTENT(INOUT):: vdata
C
      call mtx_broadcast_real4(vdata,ndata)
C
      RETURN
      END
C
      SUBROUTINE MPBCIN(vdata,ndata)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: ndata
      INTEGER,DIMENSION(ndata),INTENT(INOUT):: vdata
C
      call mtx_broadcast_integer(vdata,ndata)

      RETURN
      END
C
      SUBROUTINE MPBCKN(vdata,ndata)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: ndata
      CHARACTER(LEN=ndata),INTENT(INOUT):: vdata
C
      call mtx_broadcast_character(vdata,ndata)
C
      RETURN
      END
C
      SUBROUTINE MPBCCN(vdata,ndata)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: ndata
      COMPLEX(8),DIMENSION(ndata),INTENT(INOUT):: vdata
C
      call mtx_broadcast_complex8(vdata,ndata)
C
      RETURN
      END
C
      SUBROUTINE MPBCDA(v)
      IMPLICIT NONE
      REAL(8),INTENT(INOUT):: v
C
      call mtx_broadcast1_real8(v)
C
      RETURN
      END
C
      SUBROUTINE MPBCRA(v)
      IMPLICIT NONE
      REAL(4),INTENT(INOUT):: v
C
      call mtx_broadcast1_real4(v)
C
      RETURN
      END
C
      SUBROUTINE MPBCIA(v)
      IMPLICIT NONE
      INTEGER,INTENT(INOUT):: v
C
      call mtx_broadcast1_integer(v)
C
      RETURN
      END
C
      SUBROUTINE MPBCKA(v)
      IMPLICIT NONE
      CHARACTER(LEN=1),INTENT(INOUT):: v
C
      call mtx_broadcast1_character(v)
C
      RETURN
      END
C
      SUBROUTINE MPBCLA(v)
      IMPLICIT NONE
      LOGICAL,INTENT(INOUT):: v
C
      call mtx_broadcast1_logical(v)
C
      RETURN
      END
C
      SUBROUTINE MPBCCA(v)
      IMPLICIT NONE
      COMPLEX(8),INTENT(INOUT):: v
C
      call mtx_broadcast1_complex8(v)
C
      RETURN
      END
