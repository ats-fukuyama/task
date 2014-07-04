      MODULE FPMPI

      USE libmpi
      USE fpcomm

      CONTAINS

      SUBROUTINE fp_gatherv_real8_sav(vsend,nscnt,vreturn,n,nsa)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: nscnt
      REAL(8),DIMENSION(nrstart:nrend,NSAMAX),INTENT(IN):: vsend
      double precision,dimension(nrstart:nrend):: vtemp
!      REAL(8),DIMENSION(nrmax*nsamax):: vrecv
      REAL(8),DIMENSION(nrmax*N_partition_s):: vrecv
      REAL(8),DIMENSION(nrmax,nsamax),INTENT(OUT):: vreturn
      INTEGER,DIMENSION(nsize):: nlen, npos
      INTEGER:: n, nsa, nn, ns, nr, nsw, nse
      integer,dimension(nsize):: idisp

      DO nn=1,nsize
         idisp(nn)=savpos(nn,n)-1
      END DO
      DO NR=NRSTART,NREND
         vtemp(NR)=vsend(nr,nsa)
      END DO
!      call mtx_gather_real8(vtemp,nscnt,vrecv)
      call mtx_gather_real8(vtemp,nscnt,vrecv)

      IF(NRANK.eq.0)THEN
         nsw = NSAEND-NSASTART+1
         DO NS=1,n_partition_s
            NSE=(NS-1)*NSW + N
            DO NR=1,NRMAX
               vreturn(NR,NSE)=vrecv( (NS-1)*NRMAX+NR )
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE fp_gatherv_real8_sav
!-----
      SUBROUTINE update_fnsb 

      IMPLICIT NONE
      integer:: nsend, nth, np, nr, nsa
      double precision,dimension(nthmax,npstart:npend,nrstart:nrend,nsastart:nsaend):: dsend

      DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX
                  dsend(nth,np,nr,nsa)=FNSP(nth,np,nr,nsa)
               END DO
            END DO
         END DO
      END DO

      nsend=NTHMAX*(NPEND-NPSTART+1)*(NREND-NRSTART+1)*(NSAEND-NSASTART+1)
      CALL mtx_allgather_real8(dsend,nsend,FNSB(1:NTHMAX,NPSTART,NRSTART,1)) 

      END SUBROUTINE update_fnsb
!-----
      SUBROUTINE update_fns

      USE libmtx
      IMPLICIT NONE
      integer:: nsend, nth, np, nr, nsa, nsw, nswi,N
!      double precision,dimension(nthmax,npstart:npend,nrstart:nrend):: dsend
!      double precision,dimension(nthmax,npmax,nrmax,n_partition_s):: drecv
      double precision,dimension(nthmax,npstart:npend,nrstart:nrend):: dsend
      double precision,dimension(nthmax,npmax,nrmax):: drecv
      double precision,dimension(nthmax,npmax,nrmax,nsastart:nsaend):: dsend2

      CALL mtx_set_communicator(comm_nrnp)
      DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            DO NP=NPSTART,NPEND
               DO NTH=1,NTHMAX
                  dsend(nth,np,nr)=fnsp(nth,np,nr,nsa)
               END DO
            END DO
         END DO
         nsend=NTHMAX*(NPEND-NPSTART+1)*(NREND-NRSTART+1)
         CALL mtx_gather_real8(dsend,nsend,drecv) 
         DO NR=1,NRMAX
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  dsend2(nth,np,nr,nsa)=drecv(nth,np,nr)
               END DO
            END DO
         END DO
      END DO
      CALL mtx_set_communicator(comm_nsa)
      nsend=NTHMAX*NPMAX*NRMAX*(NSAEND-NSASTART+1)
      CALL mtx_gather_real8(dsend2,nsend,FNS) 
      CALL mtx_reset_communicator 

      END SUBROUTINE update_fns
!-----

      SUBROUTINE source_allreduce(array)

      IMPLICIT NONE
      DOUBLE PRECISION,dimension(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX), &
           INTENT(INOUT):: array
      DOUBLE PRECISION,dimension(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX):: &
           sendbuf, recvbuf
      INTEGER,dimension(NTHMAX,NPSTART:NPEND,NRSTART:NREND,NSAMAX):: &
           lloc
      INTEGER:: ierr, ncount

         ncount = NTHMAX*(NPEND-NPSTART+1)*(NREND-NRSTART+1)*NSAMAX
         sendbuf(:,:,:,:)=0.D0
         sendbuf(:,:,:,:)=array(:,:,:,:)

         CALL mtx_allreduce_real8(sendbuf, ncount, 1, recvbuf,lloc)

         array(:,:,:,:) = recvbuf(:,:,:,:)

      RETURN

      END SUBROUTINE source_allreduce
!-----
      SUBROUTINE p_theta_integration(vdata)

      IMPLICIT NONE
      double precision,intent(inout):: vdata
      double precision:: vrecv
      integer::vloc

      CALL mtx_allreduce1_real8(vdata,3,vrecv,vloc) 
      vdata = vrecv

      END SUBROUTINE p_theta_integration
!-----
      SUBROUTINE fpl_comm(vdata,vrecv)

      IMPLICIT NONE
      double precision,dimension(NPSTART:NPEND),intent(in)::vdata
      double precision,dimension(NPMAX),intent(out)::vrecv
      integer:: ndata

      ndata = (NPEND-NPSTART+1)
      CALL mtx_allgather_real8(vdata,ndata,vrecv)

      END SUBROUTINE fpl_comm
!-----
      END MODULE FPMPI
