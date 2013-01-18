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
!      call mtx_gatherv_real8(vtemp,nscnt,vrecv,nrmax*nsamax,savlen,idisp)
      call mtx_gather_real8(vtemp,nscnt,vrecv)

      IF(NRANK.eq.0)THEN
         nsw = NSAEND-NSASTART+1
!         DO NS=n,NSAMAX,nsw
         DO NS=1,n_partition_s
            NSE=(NS-1)*NSW + N
            DO NR=1,NRMAX
!               vreturn(NR,NS)=vrecv( (NS-1)*NRMAX+NR )
               vreturn(NR,NSE)=vrecv( (NS-1)*NRMAX+NR )
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE fp_gatherv_real8_sav
!-----
      SUBROUTINE update_fnsb 

      IMPLICIT NONE
      integer:: nsend, nth, np, nr, nsa, nsw, nswi
      double precision,dimension(nthmax,npmax,nrstart:nrend,nsastart:nsaend):: dsend

      DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  dsend(nth,np,nr,nsa)=FNSP(nth,np,nr,nsa)
               END DO
            END DO
         END DO
      END DO

      nsw=NSAEND-NSASTART+1
      nsend=NTHMAX*NPMAX*(NREND-NRSTART+1)*NSW
!      CALL mtx_set_communicator(comm_nr,nrank,nsize)
      DO NSWI=1,NSW
         NSA=NSASTART-1+NSWI
         CALL mtx_allgather_real8(dsend,nsend,FNSB(1:NTHMAX,1:NPMAX,NRSTART,1)) 
      END DO
!      CALL mtx_reset_communicator(nrank,nsize)

      END SUBROUTINE update_fnsb
!-----
      SUBROUTINE update_fns

      IMPLICIT NONE
      integer:: nsend, nth, np, nr, nsa, nsw, nswi,N
      double precision,dimension(nthmax,npmax,nrstart:nrend):: dsend
      double precision,dimension(nthmax,npmax,nrmax,n_partition_s):: drecv


      nsw=NSAEND-NSASTART+1
      nsend=NTHMAX*NPMAX*(NREND-NRSTART+1)
      DO NSWI=1,NSW
         NSA=NSASTART-1+NSWI
         DO NR=NRSTART,NREND
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  dsend(nth,np,nr)=FNSP(nth,np,nr,nsa)
               END DO
            END DO
         END DO
         CALL mtx_gather_real8(dsend,nsend,drecv) 
         IF(NRANK.eq.0)THEN
            N=0
            DO NSA=NSWI,NSAMAX,NSW
               N=N+1
               DO NR=1,NRMAX
                  DO NP=1,NPMAX
                     DO NTH=1,NTHMAX
                        FNS(NTH,NP,NR,NSA)=drecv(NTH,NP,NR,N)
                     END DO
                  END DO
               END DO
            END DO
         END IF
      END DO

      END SUBROUTINE update_fns
!-----

      SUBROUTINE source_allreduce(array)

      IMPLICIT NONE
      DOUBLE PRECISION,dimension(NTHMAX,NPMAX,NRSTART:NREND+1,NSAMAX), &
           INTENT(INOUT):: array
      DOUBLE PRECISION,dimension(NTHMAX,NPMAX,NRSTART:NREND+1,NSAMAX):: &
           sendbuf, recvbuf
      INTEGER,dimension(NTHMAX,NPMAX,NRSTART:NREND+1,NSAMAX):: &
           lloc
      INTEGER:: ierr, ncount

         ncount = NTHMAX*NPMAX*(NREND-NRSTART+2)*NSAMAX
         sendbuf(:,:,:,:)=0.D0
         sendbuf(:,:,:,:)=array(:,:,:,:)

         CALL mtx_allreduce_real8(sendbuf, ncount, 1, recvbuf,lloc)

         array(:,:,:,:) = recvbuf(:,:,:,:)

      RETURN

      END SUBROUTINE source_allreduce
!-----
      END MODULE FPMPI
