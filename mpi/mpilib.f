C     $Id$
C
      SUBROUTINE MPINIT(nprocs1,myrank1)
C
      include '../mpi/mpilib.inc'
C
      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world,nprocs,ierr)
      call mpi_comm_rank(mpi_comm_world,myrank,ierr)
      nprocs1=nprocs
      myrank1=myrank
C
      RETURN
      END
C
      SUBROUTINE MPTERM
C
      include '../mpi/mpilib.inc'
C
      call mpi_finalize(ierr)
C
      RETURN
      END
C
      SUBROUTINE MPSYNC
C
      include '../mpi/mpilib.inc'
C
      call mpi_barrier(mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPSETI(NMAX,NRANK,ista,iend)
C
      include '../mpi/mpilib.inc'
C
      iwork1 = NMAX/nprocs
      iwork2 = mod(NMAX,nprocs)
      ista =  NRANK   *iwork1 + min(NRANK,  iwork2) + 1
      iend = (NRANK+1)*iwork1 + min(NRANK+1,iwork2)
      RETURN
      END
C
      SUBROUTINE MPBCDN(dtmp,NDTMP)
C
      include '../mpi/mpilib.inc'
      REAL*8 dtmp(NDTMP)
C
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_bcast(dtmp,NDTMP,mpi_double_precision,
     &               0,mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPBCRN(rtmp,NRTMP)
C
      include '../mpi/mpilib.inc'
      DIMENSION rtmp(NRTMP)
C
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_bcast(rtmp,NRTMP,mpi_real,
     &               0,mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPBCIN(itmp,NITMP)
C
      include '../mpi/mpilib.inc'
      DIMENSION itmp(NITMP)
C
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_bcast(itmp,NITMP,mpi_integer,
     &               0,mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPBCKN(ktmp,NKTMP)
C
      include '../mpi/mpilib.inc'
      CHARACTER ktmp*(*)
C
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_bcast(ktmp,NKTMP,mpi_character,
     &               0,mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPBCCN(ctmp,NDTMP)
C
      include '../mpi/mpilib.inc'
      COMPLEX*16 ctmp(NDTMP)
C
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_bcast(ctmp,NDTMP,mpi_double_complex,
     &               0,mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPBCDA(D)
C
      include '../mpi/mpilib.inc'
      REAL*8 D
C
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_bcast(D,1,mpi_double_precision,
     &               0,mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPBCRA(R)
C
      include '../mpi/mpilib.inc'
C
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_bcast(R,1,mpi_real,
     &               0,mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPBCIA(I)
C
      include '../mpi/mpilib.inc'
C
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_bcast(I,1,mpi_integer,
     &               0,mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPBCKA(K)
C
      include '../mpi/mpilib.inc'
      CHARACTER K*1
C
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_bcast(K,1,mpi_character,
     &               0,mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPBCLA(L)
C
      include '../mpi/mpilib.inc'
      LOGICAL L
C
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_bcast(L,1,mpi_logical,
     &               0,mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPBCCA(C)
C
      include '../mpi/mpilib.inc'
      COMPLEX*16 C
C
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_bcast(C,1,mpi_double_complex,
     &               0,mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPGTDN(dtmp,NDTMP)
C
      include '../mpi/mpilib.inc'
C
      REAL*8 dtmp(NDTMP)
      DIMENSION istatus(mpi_status_size)
      DIMENSION ireq(0:NCPUMAX-1)
C
      if(myrank.eq.0)then
         do irank=1,nprocs-1
            call MPSETI(NDTMP,irank,ista,iend)
            if (ista.le.iend) then
               call mpi_irecv(dtmp(ista),iend-ista+1,
     &              mpi_double_precision,
     &              irank,0,mpi_comm_world,ireq(irank),ierr)
            endif
         enddo
         do irank=1,nprocs-1
            call MPSETI(NDTMP,irank,ista,iend)
            if (ista.le.iend) then
               call mpi_wait(ireq(irank),istatus,ierr)
            endif
         enddo
      else
         call MPSETI(NDTMP,myrank,ista,iend)
         if (ista.le.iend) then
            call mpi_isend(dtmp(ista),iend-ista+1,mpi_double_precision,
     &           0,0,mpi_comm_world,ireq1,ierr)
            call mpi_wait(ireq1,istatus,ierr)
         endif
      endif
C
      call mpi_barrier(mpi_comm_world,ierr)
      RETURN
      END
C
      SUBROUTINE MPGTRN(rtmp,NRTMP)
C
      include '../mpi/mpilib.inc'
C
      DIMENSION rtmp(NRTMP)
      DIMENSION istatus(mpi_status_size)
      DIMENSION ireq(0:NCPUMAX-1)
C
      if(myrank.eq.0)then
         do irank=1,nprocs-1
            call MPSETI(NRTMP,irank,ista,iend)
            call mpi_irecv(rtmp(ista),iend-ista+1,mpi_real,irank,0,
     &                     mpi_comm_world,ireq(irank),ierr)
         enddo
         do irank=1,nprocs-1
            call mpi_wait(ireq(irank),istatus,ierr)
         enddo
      else
         call MPSETI(NRTMP,myrank,ista,iend)
         call mpi_isend(rtmp(ista),iend-ista+1,mpi_real,0,0,
     &                  mpi_comm_world,ireq1,ierr)
         call mpi_wait(ireq1,istatus,ierr)
      endif
C
      call mpi_barrier(mpi_comm_world,ierr)
      RETURN
      END
C
      SUBROUTINE MPGTRV(GX,NV,GXTOT,NVTOT,NM)
C
      include '../mpi/mpilib.inc'
C
      DIMENSION GX(NV),GXTOT(NM)
      dimension ircnt(NCPUMAX)
      dimension idisp(NCPUMAX)
C
      call mpi_gather(NV,1,mpi_integer,
     &                ircnt,1,mpi_integer,
     &                0,mpi_comm_world,ierr)
C
      if(myrank.eq.0)then
         nrecv=0
         do i=1,nprocs
            idisp(i)=nrecv
            nrecv=nrecv+ircnt(i)
         enddo
         NVTOT=nrecv
      endif
C
      CALL MPBCIN(ircnt,nprocs)
      CALL MPBCIN(idisp,nprocs)
C
      call mpi_gatherv(GX,NV,mpi_real,
     &                 GXTOT,ircnt,idisp ,mpi_real,
     &                 0,mpi_comm_world,ierr)
C
      RETURN
      END
C
      SUBROUTINE MPGTCV(CX,NV,CXTOT,NVTOT,NM)
C
      include '../mpi/mpilib.inc'
C
      complex*16 CX(NV),CXTOT(NM)
      dimension ircnt(NCPUMAX)
      dimension idisp(NCPUMAX)
C
      call mpi_gather(NV,1,mpi_integer,
     &                ircnt,1,mpi_integer,
     &                0,mpi_comm_world,ierr)
C
      if(myrank.eq.0)then
         nrecv=0
         do i=1,nprocs
            idisp(i)=nrecv
            nrecv=nrecv+ircnt(i)
         enddo
         NVTOT=nrecv
      endif
C
      CALL MPBCIN(ircnt,nprocs)
      CALL MPBCIN(idisp,nprocs)
C
      call mpi_gatherv(CX,NV,mpi_double_complex,
     &                 CXTOT,ircnt,idisp ,mpi_double_complex,
     &                 0,mpi_comm_world,ierr)
C
      RETURN
      END
