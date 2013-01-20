C     $Id$
C
      SUBROUTINE MPINIT(nsize_,nrank_)
C
      USE mpi
      INCLUDE '../mpi/mpilib.inc'
      INTEGER,INTENT(OUT):: nsize_,nrank_
      INTEGER::ierr
C
      CALL MPI_Init(ierr)
      IF(ierr.NE.0) WRITE(6,*)
     &     'XX MPINIT: MPI_Init: ierr=',ierr
      ncomm=MPI_COMM_WORLD
      CALL mtx_set_communicator_global(ncomm,nrank,nsize)
      nsize_=nsize
      nrank_=nrank
C
      RETURN
      END
C
      SUBROUTINE MPTERM
C
      USE mpi
      INCLUDE '../mpi/mpilib.inc'
      INTEGER::ierr
C
      CALL MPI_Finalize(ierr)
      IF(ierr.NE.0) WRITE(6,*)
     &     'XX MPTERM: MPI_Finalaize: ierr=',ierr
      RETURN
      END
C
      SUBROUTINE MPSYNC
C
      INCLUDE '../mpi/mpilib.inc'
C
      call mtx_barrier
C
      RETURN
      END
C
      SUBROUTINE MPSETI(NMAX,NRANK_,ista,iend)
C
      include '../mpi/mpilib.inc'
      INTEGER,INTENT(IN):: NMAX,NRANK_
      INTEGER,INTENT(OUT):: ista,iend
      INTEGER:: iwork1,iwork2
C
      iwork1 = NMAX/nsize
      iwork2 = mod(NMAX,nsize)
      ista =  nrank   *iwork1 + min(nrank,  iwork2) + 1
      iend = (nrank+1)*iwork1 + min(nrank+1,iwork2)
      RETURN
      END
C
      SUBROUTINE MPBCDN(vdata,ndata)
C
      include '../mpi/mpilib.inc'
      INTEGER,INTENT(IN):: ndata
      REAL(8),DIMENSION(ndata),INTENT(INOUT):: vdata
C
      call mtx_broadcast_real8(vdata,ndata)
C
      RETURN
      END
C
      SUBROUTINE MPBCRN(vdata,ndata)
C
      include '../mpi/mpilib.inc'
      INTEGER,INTENT(IN):: ndata
      REAL(4),DIMENSION(ndata),INTENT(INOUT):: vdata
C
      call mtx_broadcast_real4(vdata,ndata)
C
      RETURN
      END
C
      SUBROUTINE MPBCIN(vdata,ndata)
C
      include '../mpi/mpilib.inc'
      INTEGER,INTENT(IN):: ndata
      INTEGER,DIMENSION(ndata),INTENT(INOUT):: vdata
C
      call mtx_broadcast_integer(vdata,ndata)

      RETURN
      END
C
      SUBROUTINE MPBCKN(vdata,ndata)
C
      include '../mpi/mpilib.inc'
      INTEGER,INTENT(IN):: ndata
      CHARACTER(LEN=ndata),INTENT(INOUT):: vdata
C
      call mtx_broadcast_character(vdata,ndata)
C
      RETURN
      END
C
      SUBROUTINE MPBCCN(vdata,ndata)
C
      include '../mpi/mpilib.inc'
      INTEGER,INTENT(IN):: ndata
      COMPLEX(8),DIMENSION(ndata),INTENT(INOUT):: vdata
C
      call mtx_broadcast_complex8(vdata,ndata)
C
      RETURN
      END
C
      SUBROUTINE MPBCDA(v)
C
      include '../mpi/mpilib.inc'
      REAL(8),INTENT(INOUT):: v
C
      call mtx_broadcast1_real8(v)
C
      RETURN
      END
C
      SUBROUTINE MPBCRA(v)
C
      include '../mpi/mpilib.inc'
      REAL(4),INTENT(INOUT):: v
C
      call mtx_broadcast1_real4(v)
C
      RETURN
      END
C
      SUBROUTINE MPBCIA(v)
C
      include '../mpi/mpilib.inc'
      INTEGER,INTENT(INOUT):: v
C
      call mtx_broadcast1_integer(v)
C
      RETURN
      END
C
      SUBROUTINE MPBCKA(v)
C
      include '../mpi/mpilib.inc'
      CHARACTER(LEN=1),INTENT(INOUT):: v
C
      call mtx_broadcast1_character(v)
C
      RETURN
      END
C
      SUBROUTINE MPBCLA(v)
C
      include '../mpi/mpilib.inc'
      LOGICAL,INTENT(INOUT):: v
C
      call mtx_broadcast1_logical(v)
C
      RETURN
      END
C
      SUBROUTINE MPBCCA(v)
C
      include '../mpi/mpilib.inc'
      COMPLEX(8),INTENT(INOUT):: v
C
      call mtx_broadcast1_complex8(v)
C
      RETURN
      END
C
      SUBROUTINE MPGTDN(vdata,ndata,vtot)
C
      include '../mpi/mpilib.inc'
      INTEGER,INTENT(IN):: ndata
      REAL(8),DIMENSION(ndata),INTENT(INOUT):: vdata
      REAL(8),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
C
      CALL mtx_gather_real8(vdata,ndata,vtot)
C
      RETURN
      END
C
      SUBROUTINE MPGTRN(vdata,ndata,vtot)
C
      include '../mpi/mpilib.inc'
      INTEGER,INTENT(IN):: ndata
      REAL(4),DIMENSION(ndata),INTENT(INOUT):: vdata
      REAL(4),DIMENSION(ndata*nsize),INTENT(OUT):: vtot
C
      CALL mtx_gather_real4(vdata,ndata,vtot)
C
      RETURN
      END
C
      SUBROUTINE MPGTRV(vdata,ndata,vtot,ntot,ntotm)
C
      include '../mpi/mpilib.inc'
      REAL(4),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      REAL(4),DIMENSION(ntotm),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntotm
      INTEGER,INTENT(OUT):: ntot
      INTEGER,DIMENSION(nsize):: ilena,iposa
      INTEGER:: i
C
      CALL mtx_allgather1_integer(ndata,ilena)
C
      ntot=0
      DO i=1,nsize
         ilena(i)=ntot
         ntot=ntot+ilena(I)
      END DO
C     
      CALL mtx_gatherv_real4(vdata,ndata,vtot,ntot,ilena,iposa)
C
      RETURN
      END
C
      SUBROUTINE MPGTCV(vdata,ndata,vtot,ntot,ntotm)
C
      include '../mpi/mpilib.inc'
      COMPLEX(8),DIMENSION(ndata),INTENT(IN):: vdata
      INTEGER,INTENT(IN):: ndata
      COMPLEX(8),DIMENSION(ntotm),INTENT(OUT):: vtot
      INTEGER,INTENT(IN):: ntotm
      INTEGER,INTENT(OUT):: ntot
      INTEGER,DIMENSION(nsize):: ilena,iposa
      INTEGER:: i
C
      CALL mtx_allgather1_integer(ndata,ilena)
C
      ntot=0
      DO i=1,nsize
         ilena(i)=ntot
         ntot=ntot+ilena(I)
      END DO
C     
      CALL mtx_gatherv_complex8(vdata,ndata,vtot,ntot,ilena,iposa)
C
      RETURN
      END
