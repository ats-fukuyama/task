      MODULE FPMPI

      USE libmpi
      USE fpcomm

      CONTAINS

      SUBROUTINE fp_gatherv_real8_sav(vsend,nscnt,vreturn,n,nsa)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: nscnt
      REAL(8),DIMENSION(nrstart:nrend,NSAMAX),INTENT(IN):: vsend
      REAL(8),DIMENSION(nrmax*nsamax):: vrecv
      REAL(8),DIMENSION(nrmax,nsamax),INTENT(OUT):: vreturn
      INTEGER,DIMENSION(nsize):: nlen, npos
      INTEGER:: n, nsa, nn, ns, nr, nsw, nse
      integer,dimension(nsize):: idisp

      DO nn=1,nsize
         idisp(nn)=savpos(nn,n)-1
      END DO

      call mtx_gatherv_real8(vsend(NRSTART,NSA),nscnt,vrecv,nscnt,savlen,idisp)

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
      END SUBROUTINE fp_gatherv_real8_sav

!-----

      SUBROUTINE source_allreduce(array,ncom)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: ncom
      DOUBLE PRECISION,dimension(NTHMAX,NPMAX,NRSTART:NREND+1,NSAMAX), &
           INTENT(INOUT):: array
      DOUBLE PRECISION,dimension(NTHMAX,NPMAX,NRSTART:NREND+1,NSAMAX):: &
           sendbuf, recvbuf
      INTEGER,dimension(NTHMAX,NPMAX,NRSTART:NREND+1,NSAMAX):: &
           lloc
      INTEGER:: ierr, ncount

         ncount = NTHMAX*NPMAX*(NREND-NRSTART+2)*NSAMAX

         CALL mtx_allreduce_real8(sendbuf, ncount, 1, recvbuf,lloc, NSASTART)

         array(:,:,:,:) = recvbuf(:,:,:,:)

      RETURN

      END SUBROUTINE source_allreduce
!-----
      END MODULE FPMPI
