      MODULE FPMPI

      USE MPI
      USE fpcomm

      CONTAINS

!------------------------------------------------------

      SUBROUTINE mtx_comm_split_s(N_partition_s,color,key,PETSC_COMM_LOCAL_S)

      USE MPI
      IMPLICIT NONE
      integer,intent(IN):: N_partition_s
      integer,intent(OUT):: color, key, PETSC_COMM_LOCAL_S
      integer:: N_PS, ierr

      N_PS=NPROCS/N_partition_s
      
      color = NRANK / N_PS ! number of belonging sub group
      key = NRANK - color*N_PS ! rank in sub group
      CALL MPI_COMM_SPLIT(ncomw, color, key, PETSC_COMM_LOCAL_S, ierr)
      IF(ierr.ne.0) WRITE(*,*) "mtx_split_s, ierr=", ierr
      call MPI_Comm_rank(PETSC_COMM_LOCAL_S,nranks,ierr)
      call MPI_Comm_size(PETSC_COMM_LOCAL_S,nprocss,ierr) 

      END SUBROUTINE mtx_comm_split_s
      
!---------------------------------------------------

      SUBROUTINE mtx_comm_split_r(N_partition_r,color,key,PETSC_COMM_LOCAL_R)

      USE MPI
      IMPLICIT NONE
      integer,intent(IN):: N_partition_r
      integer,intent(OUT):: color, key, PETSC_COMM_LOCAL_R
      integer:: N_PR, ierr

      N_PR=NPROCS/N_partition_r

      key = NRANK/N_partition_r ! rank in sub group
      color = mod(NRANK+1,N_partition_r) ! number of belonging sub group
!      color = NRANK / N_PR 
!      key = NRANK - color*N_PR 
      CALL MPI_COMM_SPLIT(ncomw, color, key, PETSC_COMM_LOCAL_R, ierr)
      IF(ierr.ne.0) WRITE(*,*) "mtx_split_r, ierr=", ierr
      call MPI_Comm_rank(PETSC_COMM_LOCAL_R,nrankr,ierr)
      call MPI_Comm_size(PETSC_COMM_LOCAL_R,nprocsr,ierr) 

      END SUBROUTINE mtx_comm_split_r

!---------------------------------------------------

      SUBROUTINE mtx_gatherv_real8_local(vdata,ndata,vtot,ntot,ilena,iposa,ncom)

      IMPLICIT NONE

      INTEGER,INTENT(IN):: ndata ! mtxlen
      INTEGER,INTENT(IN):: ncom ! communicator
      INTEGER,INTENT(INOUT):: ntot ! nrmax
      REAL(8),DIMENSION(ndata),INTENT(IN):: vdata ! work
      REAL(8),DIMENSION(ntot),INTENT(OUT):: vtot ! workg
      INTEGER,DIMENSION(nprocs):: ilena,iposa
      INTEGER:: ierr

      call MPI_GATHERV(vdata,ndata,MPI_DOUBLE_PRECISION,vtot,ilena,iposa, &
                       MPI_DOUBLE_PRECISION,0,ncom,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gatherv_real8: MPI_GATHERV: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_gatherv_real8_local

!---------------------------------------------------

      SUBROUTINE mtx_allgather_f(ncom)

      IMPLICIT NONE
      integer,intent(IN):: ncom
      double precision,dimension(nthmax,npmax,NRSTART:NREND,NSASTART:NSAEND):: sendf
      integer:: ierr, NSW, NRW, NTH, NP, NR, NSB, sendcount

      DO NSB = NSASTART, NSAEND
      DO NR = NRSTART, NREND
      DO NP = 1, NPMAX
      DO NTH = 1, NTHMAX
         sendf(NTH,NP,NR,NSB) = FNSP(NTH,NP,NR,NSB)
      END DO
      END DO
      END DO
      END DO

      NSW = NSAEND-NSASTART+1
      NRW = NREND-NRSTART+1
      sendcount = NTHMAX*NPMAX*NRW*NSW
      CALL MPI_ALLGATHER(sendf, sendcount, MPI_DOUBLE_PRECISION, &
                         FNSB , sendcount, MPI_DOUBLE_PRECISION, &
                         ncom, ierr)


      IF(ierr.ne.0) WRITE(*,*) "mtx_allgather_f, ierr=", ierr

      END SUBROUTINE mtx_allgather_f
!---------------------------------------------------

      SUBROUTINE mtx_allgather_integer_sav(isend,n,nsaw,ncom)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: isend, nsaw, n, ncom
!      INTEGER,DIMENSION(nprocs,nsaw),INTENT(OUT):: irecv
      INTEGER:: ierr

      call MPI_ALLGATHER(isend,1,MPI_INTEGER,savpos(1,n),1,MPI_INTEGER, &
                         ncom,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_allgather_integer_sav: MPI_ALLGATHER: ierr=',ierr

      RETURN
      END SUBROUTINE mtx_allgather_integer_sav
!-----
      SUBROUTINE mtx_gatherv_real8_sav(vsend,nscnt,vreturn,n,nsa,ncom)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: nscnt
      REAL(8),DIMENSION(nrstart:nrend,NSAMAX),INTENT(IN):: vsend
      REAL(8),DIMENSION(nrmax*nsamax):: vrecv
      REAL(8),DIMENSION(nrmax,nsamax),INTENT(OUT):: vreturn
      INTEGER,DIMENSION(nprocs):: nlen, npos
      INTEGER:: ierr, ncom, n, nsa, nn, ns, nr, nsw, nse
      integer,dimension(nprocs):: idisp

      DO nn=1,nprocs
         idisp(nn)=savpos(nn,n)-1
      END DO

      call MPI_GATHERV(vsend(NRSTART,NSA), nscnt, MPI_DOUBLE_PRECISION, &
                       vrecv  ,savlen,idisp, MPI_DOUBLE_PRECISION, &
                       0,ncom,ierr)
      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gatherv_real8: MPI_GATHERV: ierr=',ierr
      IF(NRANK.eq.0)THEN
         nsw = NSAEND-NSASTART+1
         nse = NSAMAX/NSW
         DO NS=n,NSAMAX,nsw
            DO NR=1,NRMAX
               vreturn(NR,NS)=vrecv( (NS-1)*NRMAX+NR )
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE mtx_gatherv_real8_sav

!-----
      SUBROUTINE mtx_maxloc_real8_local(din,ncount,rank,ncom,dout,loc)

      IMPLICIT NONE
      REAL(8),dimension(nsastart:nsaend),INTENT(IN):: din
      INTEGER,INTENT(IN):: ncount, ncom
      REAL(8),dimension(2,ncount):: d_send, d_recv
      REAL(8),dimension(nsastart:nsaend),INTENT(OUT):: dout
      INTEGER,dimension(nsastart:nsaend),INTENT(OUT):: loc
      INTEGER:: ierr, i, rank

      DO i = 1, ncount
         d_send(1,i) = din(i+NSASTART-1)
         d_send(2,i) = rank*1.D0
      END DO

      CALL MPI_ALLREDUCE(d_send,d_recv,ncount,MPI_2DOUBLE_PRECISION &
                              ,MPI_MAXLOC,ncom,ierr)

      DO i = 1, ncount 
         dout(i+NSASTART-1) = d_recv(1,i)
         loc(i+NSASTART-1) = idint(d_recv(2,i))
      END DO

      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_maxloc_real8: MPI_ALLREDUCE_MAXLOC: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_maxloc_real8_local
!---------------------------------------------------

      SUBROUTINE mtx_gather_real8_deps(dsend,nsw,ncom,drecv)

      IMPLICIT NONE

      INTEGER,INTENT(IN):: nsw ! mtxlen
      INTEGER,INTENT(IN):: ncom ! communicator
      double precision,DIMENSION(nsastart:nsaend),INTENT(IN):: dsend ! work
      double precision,DIMENSION(nsamax),INTENT(OUT):: drecv ! workg
      INTEGER:: ierr

      call MPI_GATHER(dsend,nsw,MPI_DOUBLE_PRECISION, & 
                      drecv,nsw,MPI_DOUBLE_PRECISION,0,ncom,ierr)

      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_real8_deps: MPI_GATHER: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_gather_real8_deps

!------------------------------------------

      SUBROUTINE mtx_gather_integer_deps(isend,nsw,ncom,irecv)

      IMPLICIT NONE

      INTEGER,INTENT(IN):: nsw ! mtxlen
      INTEGER,INTENT(IN):: ncom ! communicator
      integer,DIMENSION(nsastart:nsaend),INTENT(IN):: isend ! work
      integer,DIMENSION(nsamax),INTENT(OUT):: irecv ! workg
      INTEGER:: ierr

      call MPI_GATHER(isend,nsw,MPI_INTEGER, & 
                      irecv,nsw,MPI_INTEGER,0,ncom,ierr)

      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_integer_deps: MPI_GATHER: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_gather_integer_deps

!------------------------------------------

      SUBROUTINE mtx_gather_fns_rs(nsa,ncom,local,global)

      IMPLICIT NONE
      integer,intent(in):: nsa, ncom
      REAL(8),dimension(NTHMAX,NPMAX,NRSTART-1:NREND+1,NSASTART:NSAEND),INTENT(IN):: local
      REAL(8),dimension(NTHMAX,NPMAX,NRMAX,NSAMAX),INTENT(OUT):: global
      INTEGER:: ierr, nsend, nrecv

      nsend = NTHMAX*NPMAX*(NREND-NRSTART+1) 
      nrecv = nsend
      call MPI_GATHER(local(1,1,nrstart,nsa),nsend,MPI_DOUBLE_PRECISION, &
                     global(1,1,nrstart,nsa),nrecv,MPI_DOUBLE_PRECISION, &
                     0, ncom, ierr)

      IF(ierr.NE.0) WRITE(6,*) &
           'XX mtx_gather_fns_rs: MPI_GATHER: ierr=',ierr
      RETURN
      END SUBROUTINE mtx_gather_fns_rs

!-----

      END MODULE FPMPI
